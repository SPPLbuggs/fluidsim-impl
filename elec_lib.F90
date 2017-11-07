module elec_lib
    use props
    use ptcl_props
    implicit none
    
    contains
    
    ! *** Electron Continuity ***
    subroutine elecEqn(g, i, j, ph, ni, ne, nte, nm, b)
        type(grid), intent(in) :: g
        integer, intent(in)    :: i, j
        real(8), intent(in)    :: ph(:,:), ni(:,:), ne(:,:), nte(:,:), nm(:,:)
        real(8), intent(out)   :: b
        real(8) :: dflx, Te, src, k_ir, k_si, flx_x(2), flx_y(2)
        
        call elecDFlx(g, i, j, ph, ne, ni, nte, dflx, flx_x, flx_y)
        
        ! rates and coefficients
        Te = get_Te(nte(i,j), ne(i,j))
        k_ir = get_k_ir(Te)
        k_si = get_k_si(Te)

        ! evaluate source terms
        src =   k_ir * ninf    * ne(i,j) &
              - beta * ni(i,j) * ne(i,j) &
              + k_si * nm(i,j) * ne(i,j) &
              + k_mp * nm(i,j)**2
        
        ! evaluate expression
        b = -dflx + src
    end subroutine
    
    ! *** Electron Energy Continuity ***
    subroutine elecEnrgEqn(g, i, j, ph, ne, ni, nte, nm, b)
        type(grid), intent(in) :: g
        integer, intent(in)    :: i, j
        real(8), intent(in)    :: ph(:,:), ne(:,:), ni(:,:), nte(:,:), nm(:,:)
        real(8), intent(out)   :: b
        real(8) :: Ex, Ey, Te, dflx, src(3), k_ir, k_ex, k_sc, k_si, nu, flxe_x(2), flxe_y(2)
        
        call elecDFlx(g, i, j, ph, ne, ni, nte, dflx, flxe_x, flxe_y)
        call elecDEnrgFlx(g, i, j, ph, ne, ni, nte, dflx)
        
        flxe_x(1) = 0.5 * sum(flxe_x)
        Ex = -((ph(i+1,j) - ph(i,j)) / g%dx(i) &
               +(ph(i,j) - ph(i-1,j)) / g%dx(i-1) &
              ) / 2.0
        
        flxe_y(1) = 0.5 * sum(flxe_y)
        Ey = 0
        
        ! rates and coefficients
        Te = get_Te(nte(i,j), ne(i,j))
        k_ir = get_k_ir(Te)
        k_sc = get_k_sc(Te)
        k_si = get_k_si(Te)
        k_ex = get_k_ex(Te)
        nu   = get_nu(Te)

        ! evaluate source terms   
        ! -e flux_e . E
        src(1) = -(flxe_x(1) * Ex + flxe_y(1) * Ey)

        ! -me/mg nu_e (Te - Tg)
        src(2) = -nte(i,j) * nu * me/mi

        ! reactions
        src(3) = - h_ir * k_ir * ninf * ne(i,j) &
                 - h_si * k_si * nm(i,j) * ne(i,j) &
                 - h_ex * k_ex * ninf * ne(i,j)    &
                 - h_sc * k_sc * nm(i,j) * ne(i,j)
        
        ! evaluate expression
        b = -dflx + sum(src)
    
    end subroutine

    ! *** Calculate Divergence of Electron Flux ***
    subroutine elecDFlx(g, i, j, ph, ne, ni, nte, dflx, flx_x, flx_y)
        type(grid), intent(in) :: g
        integer, intent(in) :: i, j
        real(8), intent(in) :: ph(:,:), ne(:,:), ni(:,:), nte(:,:)
        real(8), intent(out) :: dflx, flx_x(2), flx_y(2)
        real(8) :: a, Te(3), mu(2), D(2), ve, Ex(2), flxi

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
            
            ! Flux at i - 1/2
            call getFlx(flx_x(1), Ex(1), g%dx(i-1), -1, mu(1), D(1), &
                          ne(i-1,j), ne(i,j))

            ! Flux at i + 1/2
            call getFlx(flx_x(2), Ex(2), g%dx(i), -1, mu(2), D(2), &
                          ne(i,j), ne(i+1,j))
        
        ! - left -
        else if (g%type_x(i-1,j-1) < 0) then
            ! rates and coefficients
            Te(2) = get_Te(nte(i,j),   ne(i,j))
            Te(3) = get_Te(nte(i+1,j), ne(i+1,j))
            
            mu(2) = 0.5 * (get_mue(Te(2)) + get_mue(Te(3)))
            D(2) = 0.5 * (get_De(Te(2)) + get_De(Te(3)))

            ! Flux at i + 1/2
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
                
                ! Flux at i - 1/2
                flxi = (1 - a) * mui * Ex(1) * ni(i,j) - 0.25 * vi * ni(i,j)
                
                flx_x(1) = - a * mu(1) * Ex(1) * ne(i,j) &
                           - 0.25 * ve * ne(i,j) &
                           - gam * flxi


            ! - vacuum -
            else if (g%type_x(i-1,j-1) == -1) then
                ! Flux at i - 1/2
                flx_x(1) = 0
            end if
        
        ! - right -
        else if (g%type_x(i-1,j-1) > 0) then
            ! rates and coefficients
            Te(1) = get_Te(nte(i-1,j), ne(i-1,j))
            Te(2) = get_Te(nte(i,j),   ne(i,j))
            
            mu(1) = 0.5 * (get_mue(Te(1)) + get_mue(Te(2)))        
            D(1) = 0.5 * (get_De(Te(1)) + get_De(Te(2)))
            
            ! Flux at i - 1/2
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
                
                ! Flux at i + 1/2
                flxi = (1 - a) * mui * Ex(2) * ni(i,j) + 0.25 * vi * ni(i,j)
                
                flx_x(2) = - a * mu(2) * Ex(2) * ne(i,j) &
                           + 0.25 * ve * ne(i,j) &
                           - gam * flxi

            ! - vacuum -
            else if (g%type_x(i-1,j-1) == 1) then
                ! Flux at i + 1/2
                flx_x(2) = 0
            end if
        end if
        
        dflx = (flx_x(2) - flx_x(1)) / g%dlx(i-1)
    end subroutine
    
    ! *** Calculate Divergence of Electron Energy Flux ***
    subroutine elecDEnrgFlx(g, i, j, ph, ne, ni, nte, dflx)
        type(grid), intent(in) :: g
        integer, intent(in) :: i, j
        real(8), intent(in) :: ph(:,:), ne(:,:), ni(:,:), nte(:,:)
        real(8), intent(out) :: dflx
        real(8) :: a, flx_x(2), flx_y(2), Te(3), mu(2), D(2), ve, Ex(2), flxi

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
            
            ! Flux at i - 1/2
            call getFlx(flx_x(1), Ex(1), g%dx(i-1), -1, mu(1), D(1), &
                          nte(i-1,j), nte(i,j))

            ! Flux at i + 1/2
            call getFlx(flx_x(2), Ex(2), g%dx(i), -1, mu(2), D(2), &
                          nte(i,j), nte(i+1,j))
        
        ! - left -
        else if (g%type_x(i-1,j-1) < 0) then
            ! rates and coefficients
            Te(2) = get_Te(nte(i,j),   ne(i,j))
            Te(3) = get_Te(nte(i+1,j), ne(i+1,j))
            
            mu(2) = 0.5 * (get_mut(Te(2)) + get_mut(Te(3)))
            D(2) = 0.5 * (get_Dt(Te(2)) + get_Dt(Te(3)))

            ! Flux at i + 1/2
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
                
                ! Flux at i - 1/2
                flxi = (1 - a) * mui * Ex(1) * ni(i,j) - 0.25 * vi * ni(i,j)
                
                flx_x(1) = - a * mu(1) * Ex(1) * nte(i,j) &
                           - 1.0/3.0 * ve * nte(i,j) &
                           - gam * Te(2) * flxi


            ! - vacuum -
            else if (g%type_x(i-1,j-1) == -1) then
                ! Flux at i - 1/2
                flx_x(1) = 0
            end if
        
        ! - right -
        else if (g%type_x(i-1,j-1) > 0) then
            ! rates and coefficients
            Te(1) = get_Te(nte(i-1,j), ne(i-1,j))
            Te(2) = get_Te(nte(i,j),   ne(i,j))
            
            mu(1) = 0.5 * (get_mut(Te(1)) + get_mut(Te(2)))        
            D(1) = 0.5 * (get_Dt(Te(1)) + get_Dt(Te(2)))
            
            ! Flux at i - 1/2
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
                
                ! Flux at i + 1/2
                flxi = (1 - a) * mui * Ex(2) * ni(i,j) + 0.25 * vi * ni(i,j)
                
                flx_x(2) = - a * mu(2) * Ex(2) * nte(i,j) &
                           + 1.0/3.0 * ve * nte(i,j) &
                           - gam * Te(2) * flxi

            ! - vacuum -
            else if (g%type_x(i-1,j-1) == 1) then
                ! Flux at i + 1/2
                flx_x(2) = 0
            end if
        end if
        
        dflx = (flx_x(2) - flx_x(1)) / g%dlx(i-1)
    end subroutine
end module
    
