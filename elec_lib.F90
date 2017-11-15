module elec_lib
  use props
  use ptcl_props
  implicit none

  real(8), allocatable :: ne_pl(:,:,:), ne_mi(:,:,:), &
                          nte_pl(:,:,:), nte_mi(:,:,:)

  contains

  ! *** Electron Source Term ***
  subroutine elecSrc(ne, ni, nte, nm, src)
    real(8), intent(in)  :: ni, ne, nte, nm
    real(8), intent(out) :: src
    real(8) :: Te, k_ir, k_si

    ! rates and coefficients
    Te = get_Te(nte, ne)
    k_ir = get_k_ir(Te)
    k_si = get_k_si(Te)

    t_m = min(t_m, 1.0 / (get_mue(Te) * ne + ni * mui))

    ! evaluate source terms
    src =   k_ir * ninf    * ne &
          - beta * ni * ne &
          + k_si * nm * ne &
          + k_mp * nm**2
  end subroutine

  ! *** Electron Energy Source Term ***
  subroutine elecEnrgSrc(g, i, j, ph, ne, ni, nte, nm, src)
    type(grid), intent(in) :: g
    integer, intent(in)    :: i, j
    real(8), intent(in)    :: ph(:,:), ne(:,:), ni(:,:), nte(:,:), nm(:,:)
    real(8), intent(out)   :: src
    real(8) :: Ex, Ey, Te, k_ir, k_ex, k_sc, k_si, nu, flxe_x(2), flxe_y(2)

    call elecFlx(g, i, j, ph, ne, ni, nte, flxe_x, flxe_y)

    flxe_x(1) = 0.5 * sum(flxe_x)
    flxe_y(1) = 0.5 * sum(flxe_y)

    Ex = 0
    if (g%nx > 1) then
      if (g%type_x(i-1,j-1) == -1) then
        Ex = -((ph(i+1,j) - ph(i,j)) / g%dx(i) &
              ) / 2.0
      else if(g%type_x(i-1,j-1) == 1) then
        Ex = -((ph(i,j) - ph(i-1,j)) / g%dx(i-1) &
              ) / 2.0
      else
        Ex = -((ph(i+1,j) - ph(i,j)) / g%dx(i) &
               +(ph(i,j) - ph(i-1,j)) / g%dx(i-1) &
              ) / 2.0
      end if
    end if

    Ey = 0
    if (g%ny > 1) then
      if (g%type_y(i-1,j-1) == -1) then
        Ey = -((ph(i,j+1) - ph(i,j)) / g%dy(j) &
              ) / 2.0
      else if(g%type_y(i-1,j-1) == 1) then
        Ey = -((ph(i,j) - ph(i,j-1)) / g%dy(j-1) &
              ) / 2.0
      else
        Ey = -((ph(i,j+1) - ph(i,j)) / g%dy(j) &
               +(ph(i,j) - ph(i,j-1)) / g%dy(j-1) &
              ) / 2.0
      end if
    end if

    ! rates and coefficients
    Te = get_Te(nte(i,j), ne(i,j))
    k_ir = get_k_ir(Te)
    k_sc = get_k_sc(Te)
    k_si = get_k_si(Te)
    k_ex = get_k_ex(Te)
    nu   = get_nu(Te)

    ! evaluate source term
    src = -flxe_x(1) * Ex - flxe_y(1) * Ey &
          -nte(i,j) * nu * me/mi           &
          -h_ir * k_ir * ninf    * ne(i,j) &
          -h_si * k_si * nm(i,j) * ne(i,j) &
          -h_ex * k_ex * ninf    * ne(i,j) &
          -h_sc * k_sc * nm(i,j) * ne(i,j)
  end subroutine

  ! *** Divergence of Electron Flux ***
  subroutine elecDFlx(g, i, j, ph, ne, ni, nte, dflx)
    type(grid), intent(in) :: g
    integer, intent(in) :: i, j
    real(8), intent(in) :: ph(:,:), ne(:,:), ni(:,:), nte(:,:)
    real(8), intent(out) :: dflx
    real(8) :: flx_x(2), flx_y(2), dfdx = 0, dfdy = 0

    call elecFlx(g, i, j, ph, ne, ni, nte, flx_x, flx_y)

    if (g%nx > 1) dfdx = (flx_x(2) - flx_x(1)) / g%dlx(i-1)
    if (g%ny > 1) then
      if (cyl) then
        dfdy = (flx_y(2) * (g%r(j+1) + g%r(j)) / 2.0   &
               - flx_y(1) * (g%r(j) + g%r(j-1)) / 2.0) &
               / g%dly(j-1) / g%r(j)
      else
        dfdy = (flx_y(2) - flx_y(1)) / g%dly(j-1)
      end if
    end if

    dflx = dfdx + dfdy
  end subroutine

  ! *** Divergence of Electron Energy Flux ***
  subroutine elecEnrgDFlx(g, i, j, ph, ne, ni, nte, dflx)
    type(grid), intent(in) :: g
    integer, intent(in) :: i, j
    real(8), intent(in) :: ph(:,:), ne(:,:), ni(:,:), nte(:,:)
    real(8), intent(out) :: dflx
    real(8) :: flx_x(2), flx_y(2), dfdx = 0, dfdy = 0

    call elecEnrgFlx(g, i, j, ph, ne, ni, nte, flx_x, flx_y)

    if (g%nx > 1) dfdx = (flx_x(2) - flx_x(1)) / g%dlx(i-1)
    if (g%ny > 1) then
      if (cyl) then
        dfdy = (flx_y(2) * (g%r(j+1) + g%r(j)) / 2.0   &
               - flx_y(1) * (g%r(j) + g%r(j-1)) / 2.0) &
               / g%dly(j-1) / g%r(j)
      else
        dfdy = (flx_y(2) - flx_y(1)) / g%dly(j-1)
      end if
    end if

    dflx = dfdx + dfdy
  end subroutine

  ! *** Electron Flux ***
  subroutine elecFlx(g, i, j, ph, ne, ni, nte, flx_x, flx_y)
    type(grid), intent(in) :: g
    integer, intent(in) :: i, j
    real(8), intent(in) :: ph(:,:), ne(:,:), ni(:,:), nte(:,:)
    real(8), intent(out) :: flx_x(2), flx_y(2)
    real(8) :: a, Te(3), mu(2), D(2), ve, Ex(2), Ey(2), flxi

    flx_x = 0
    flx_y = 0

    ! X-dir fields:
    Ex(1) = -(ph(i,j) - ph(i-1,j)) / g%dx(i-1)
    Ex(2) = -(ph(i+1,j) - ph(i,j)) / g%dx(i)

    ! X-dir Fluxes:
    ! - center -
    if (g%type_x(i-1,j-1) == 0) then
      ! rates and coefficients
      Te(1) = get_Te(nte(i-1,j), ne(i-1,j))
      Te(2) = get_Te(nte(i,j),   ne(i,j))
      Te(3) = get_Te(nte(i+1,j), ne(i+1,j))

      mu(1) = 0.5 * (get_mue(Te(1)) + get_mue(Te(2)))
      mu(2) = 0.5 * (get_mue(Te(2)) + get_mue(Te(3)))

      D(1) = 0.5 * (get_De(Te(1)) + get_De(Te(2)))
      D(2) = 0.5 * (get_De(Te(2)) + get_De(Te(3)))

      call getFlx(flx_x(1), Ex(1), g%dx(i-1), -1, mu(1), D(1), &
                  ne(i-1,j), ne(i,j))
      call getFlx(flx_x(2), Ex(2), g%dx(i), -1, mu(2), D(2), &
                  ne(i,j), ne(i+1,j))

    ! - left -
    else if (g%type_x(i-1,j-1) < 0) then
      ! rates and coefficients
      Te(2) = get_Te(nte(i,j),   ne(i,j))
      Te(3) = get_Te(nte(i+1,j), ne(i+1,j))

      mu(2) = 0.5 * (get_mue(Te(2)) + get_mue(Te(3)))
      D(2) = 0.5 * (get_De(Te(2)) + get_De(Te(3)))

      call getFlx(flx_x(2), Ex(2), g%dx(i), -1, mu(2), D(2), &
                  ne(i,j), ne(i+1,j))

      ! - electrode -
      if (g%type_x(i-1,j-1) == -2) then
        if (Ex(1) > 0) then
          a = 1
        else
          a = 0
        end if

        mu(1) = get_mue(Te(2))
        ve = sqrt((16.0 * e * ph0 * Te(2)) / (3.0 * pi * me)) * t0 / x0

        flxi = (1 - a) * mui * Ex(1) * ni(i,j) - 0.25 * vi * ni(i,j)

        flx_x(1) = - a * mu(1) * Ex(1) * ne(i,j) &
                   - 0.25 * ve * ne(i,j) &
                   - gam * flxi

      ! - vacuum -
      else if (g%type_x(i-1,j-1) == -1) then
        flx_x(1) = 0
      end if

    ! - right -
    else if (g%type_x(i-1,j-1) > 0) then
      ! rates and coefficients
      Te(1) = get_Te(nte(i-1,j), ne(i-1,j))
      Te(2) = get_Te(nte(i,j),   ne(i,j))

      mu(1) = 0.5 * (get_mue(Te(1)) + get_mue(Te(2)))
      D(1) = 0.5 * (get_De(Te(1)) + get_De(Te(2)))

      call getFlx(flx_x(1), Ex(1), g%dx(i-1), -1, mu(1), D(1), &
                  ne(i-1,j), ne(i,j))

      ! - electrode -
      if (g%type_x(i-1,j-1) == 2) then
        if (-Ex(2) > 0) then
          a = 1
        else
          a = 0
        end if

        mu(2) = get_mue(Te(2))
        ve = sqrt((16.0 * e * ph0 * Te(2)) / (3.0 * pi * me)) * t0 / x0

        flxi = (1 - a) * mui * Ex(2) * ni(i,j) + 0.25 * vi * ni(i,j)

        flx_x(2) = - a * mu(2) * Ex(2) * ne(i,j) &
                   + 0.25 * ve * ne(i,j) &
                   - gam * flxi

      ! - vacuum -
      else if (g%type_x(i-1,j-1) == 1) then
        flx_x(2) = 0
      end if
    end if

    ! Y-dir Fluxes
    if (g%ny > 1) then
      Ey(1) = -(ph(i,j) - ph(i,j-1)) / g%dy(j-1)
      Ey(2) = -(ph(i,j+1) - ph(i,j)) / g%dy(j)

      ! - center -
      if (g%type_y(i-1,j-1) == 0) then
        ! rates and coefficients
        Te(1) = get_Te(nte(i,j-1), ne(i,j-1))
        Te(2) = get_Te(nte(i,j),   ne(i,j))
        Te(3) = get_Te(nte(i,j+1), ne(i,j+1))

        mu(1) = 0.5 * (get_mue(Te(1)) + get_mue(Te(2)))
        mu(2) = 0.5 * (get_mue(Te(2)) + get_mue(Te(3)))

        D(1) = 0.5 * (get_De(Te(1)) + get_De(Te(2)))
        D(2) = 0.5 * (get_De(Te(2)) + get_De(Te(3)))

        call getFlx(flx_y(1), Ey(1), g%dy(j-1), -1, mu(1), D(1), &
                    ne(i,j-1), ne(i,j))
        call getFlx(flx_y(2), Ey(2), g%dy(j), -1, mu(2), D(2), &
                    ne(i,j), ne(i,j+1))

      ! - left -
      else if (g%type_y(i-1,j-1) < 0) then
        ! rates and coefficients
        Te(2) = get_Te(nte(i,j),   ne(i,j))
        Te(3) = get_Te(nte(i,j+1), ne(i,j+1))

        mu(2) = 0.5 * (get_mue(Te(2)) + get_mue(Te(3)))
        D(2) = 0.5 * (get_De(Te(2)) + get_De(Te(3)))

        call getFlx(flx_y(2), Ey(2), g%dy(j), -1, mu(2), D(2), &
                    ne(i,j), ne(i,j+1))

        flx_y(1) = 0

      ! - right -
      else if (g%type_y(i-1,j-1) > 0) then
        ! rates and coefficients
        Te(1) = get_Te(nte(i,j-1), ne(i,j-1))
        Te(2) = get_Te(nte(i,j),   ne(i,j))

        mu(1) = 0.5 * (get_mue(Te(1)) + get_mue(Te(2)))
        D(1) = 0.5 * (get_De(Te(1)) + get_De(Te(2)))

        call getFlx(flx_y(1), Ey(1), g%dy(j-1), -1, mu(1), D(1), &
                    ne(i,j-1), ne(i,j))

        if (rwall) then
          if (-Ey(2) > 0) then
            a = 1
          else
            a = 0
          end if

          mu(2) = get_mue(Te(2))
          ve = sqrt((16.0 * e * ph0 * Te(2)) / (3.0 * pi * me)) * t0 / x0

          ! Flux at j + 1/2
          flx_y(2) = - a * mu(2) * Ey(2) * ne(i,j) &
                     + 0.25 * ve * ne(i,j)

        else
          flx_y(2) = 0
        end if
      end if
    end if
  end subroutine

  ! *** Electron Energy Flux ***
  subroutine elecEnrgFlx(g, i, j, ph, ne, ni, nte, flx_x, flx_y)
    type(grid), intent(in) :: g
    integer, intent(in) :: i, j
    real(8), intent(in) :: ph(:,:), ne(:,:), ni(:,:), nte(:,:)
    real(8), intent(out) :: flx_x(2), flx_y(2)
    real(8) :: a, Te(3), mu(2), D(2), ve, Ex(2), Ey(2), flxi

    flx_x = 0
    flx_y = 0

    ! X-dir fields:
    Ex(1) = -(ph(i,j) - ph(i-1,j)) / g%dx(i-1)
    Ex(2) = -(ph(i+1,j) - ph(i,j)) / g%dx(i)

    ! X-dir Fluxes:
    ! - center -
    if (g%type_x(i-1,j-1) == 0) then
      ! rates and coefficients
      Te(1) = get_Te(nte(i-1,j), ne(i-1,j))
      Te(2) = get_Te(nte(i,j),   ne(i,j))
      Te(3) = get_Te(nte(i+1,j), ne(i+1,j))

      mu(1) = 0.5 * (get_mut(Te(1)) + get_mut(Te(2)))
      mu(2) = 0.5 * (get_mut(Te(2)) + get_mut(Te(3)))

      D(1) = 0.5 * (get_Dt(Te(1)) + get_Dt(Te(2)))
      D(2) = 0.5 * (get_Dt(Te(2)) + get_Dt(Te(3)))

      call getFlx(flx_x(1), Ex(1), g%dx(i-1), -1, mu(1), D(1), &
                  nte(i-1,j), nte(i,j))
      call getFlx(flx_x(2), Ex(2), g%dx(i), -1, mu(2), D(2), &
                  nte(i,j), nte(i+1,j))

    ! - left -
    else if (g%type_x(i-1,j-1) < 0) then
      ! rates and coefficients
      Te(2) = get_Te(nte(i,j),   ne(i,j))
      Te(3) = get_Te(nte(i+1,j), ne(i+1,j))

      mu(2) = 0.5 * (get_mut(Te(2)) + get_mut(Te(3)))
      D(2) = 0.5 * (get_Dt(Te(2)) + get_Dt(Te(3)))

      call getFlx(flx_x(2), Ex(2), g%dx(i), -1, mu(2), D(2), &
                  nte(i,j), nte(i+1,j))

      ! - electrode -
      if (g%type_x(i-1,j-1) == -2) then
        if (Ex(1) > 0) then
          a = 1
        else
          a = 0
        end if

        mu(1) = get_mut(Te(2))
        ve = sqrt((16.0 * e * ph0 * Te(2)) / (3.0 * pi * me)) * t0 / x0

        flxi = (1 - a) * mui * Ex(1) * ni(i,j) - 0.25 * vi * ni(i,j)

        flx_x(1) = - a * mu(1) * Ex(1) * nte(i,j) &
                   - 1.0/3.0 * ve * nte(i,j) &
                   - gam * Te(2) * flxi

      ! - vacuum -
      else if (g%type_x(i-1,j-1) == -1) then
        flx_x(1) = 0
      end if

    ! - right -
    else if (g%type_x(i-1,j-1) > 0) then
      ! rates and coefficients
      Te(1) = get_Te(nte(i-1,j), ne(i-1,j))
      Te(2) = get_Te(nte(i,j),   ne(i,j))

      mu(1) = 0.5 * (get_mut(Te(1)) + get_mut(Te(2)))
      D(1) = 0.5 * (get_Dt(Te(1)) + get_Dt(Te(2)))

      call getFlx(flx_x(1), Ex(1), g%dx(i-1), -1, mu(1), D(1), &
                  nte(i-1,j), nte(i,j))

      ! - electrode -
      if (g%type_x(i-1,j-1) == 2) then
        if (-Ex(2) > 0) then
          a = 1
        else
          a = 0
        end if

        mu(2) = get_mut(Te(2))
        ve = sqrt((16.0 * e * ph0 * Te(2)) / (3.0 * pi * me)) * t0 / x0

        flxi = (1 - a) * mui * Ex(2) * ni(i,j) + 0.25 * vi * ni(i,j)

        flx_x(2) = - a * mu(2) * Ex(2) * nte(i,j) &
                   + 1.0/3.0 * ve * nte(i,j) &
                   - gam * Te(2) * flxi

      ! - vacuum -
      else if (g%type_x(i-1,j-1) == 1) then
        flx_x(2) = 0
      end if
    end if

    ! Y-dir Fluxes
    if (g%ny > 1) then
      Ey(1) = -(ph(i,j) - ph(i,j-1)) / g%dy(j-1)
      Ey(2) = -(ph(i,j+1) - ph(i,j)) / g%dy(j)

      ! - center -
      if (g%type_y(i-1,j-1) == 0) then
        ! rates and coefficients
        Te(1) = get_Te(nte(i,j-1), ne(i,j-1))
        Te(2) = get_Te(nte(i,j),   ne(i,j))
        Te(3) = get_Te(nte(i,j+1), ne(i,j+1))

        mu(1) = 0.5 * (get_mut(Te(1)) + get_mut(Te(2)))
        mu(2) = 0.5 * (get_mut(Te(2)) + get_mut(Te(3)))

        D(1) = 0.5 * (get_Dt(Te(1)) + get_Dt(Te(2)))
        D(2) = 0.5 * (get_Dt(Te(2)) + get_Dt(Te(3)))

        ! Flux at j - 1/2
        call getFlx(flx_y(1), Ey(1), g%dy(j-1), -1, mu(1), D(1), &
                    nte(i,j-1), nte(i,j))
        call getFlx(flx_y(2), Ey(2), g%dy(j), -1, mu(2), D(2), &
                    nte(i,j), nte(i,j+1))

      ! - left -
      else if (g%type_y(i-1,j-1) < 0) then
        ! rates and coefficients
        Te(2) = get_Te(nte(i,j),   ne(i,j))
        Te(3) = get_Te(nte(i,j+1), ne(i,j+1))

        mu(2) = 0.5 * (get_mut(Te(2)) + get_mut(Te(3)))
        D(2) = 0.5 * (get_Dt(Te(2)) + get_Dt(Te(3)))

        call getFlx(flx_y(2), Ey(2), g%dy(j), -1, mu(2), D(2), &
                    nte(i,j), nte(i,j+1))

        flx_y(1) = 0

      ! - right -
      else if (g%type_y(i-1,j-1) > 0) then
        ! rates and coefficients
        Te(1) = get_Te(nte(i,j-1), ne(i,j-1))
        Te(2) = get_Te(nte(i,j),   ne(i,j))

        mu(1) = 0.5 * (get_mut(Te(1)) + get_mut(Te(2)))
        D(1) = 0.5 * (get_Dt(Te(1)) + get_Dt(Te(2)))

        call getFlx(flx_y(1), Ey(1), g%dy(j-1), -1, mu(1), D(1), &
                    nte(i,j-1), nte(i,j))

        if (rwall) then
          if (-Ey(2) > 0) then
            a = 1
          else
            a = 0
          end if

          mu(2) = get_mut(Te(2))
          ve = sqrt((16.0 * e * ph0 * Te(2)) / (3.0 * pi * me)) * t0 / x0

          ! Flux at j + 1/2
          flx_y(2) = - a * mu(2) * Ey(2) * nte(i,j) &
                     + 1.0/3.0 * ve * nte(i,j)
        else
          flx_y(2) = 0
        end if
      end if
    end if
  end subroutine
end module
