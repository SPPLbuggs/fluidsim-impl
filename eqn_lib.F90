module eqn_lib
  use props
  use lapl_lib
  use elec_lib
  use ion_lib
  implicit none

contains

  subroutine eqn_init(g)
    type(grid), intent(in) :: g
    allocate(ph_pl(g%bx+2, g%by+2, 1), ph_mi(g%bx+2, g%by+2, 1), &
             ne_pl(g%bx+2, g%by+2, 2), ne_mi(g%bx+2, g%by+2, 2), &
             ni_pl(g%bx+2, g%by+2, 1), ni_mi(g%bx+2, g%by+2, 1), &
             nm_pl(g%bx+2, g%by+2, 1), nm_mi(g%bx+2, g%by+2, 1))

    ph_pl  = 0
    ne_pl  = n_init
    ni_pl  = n_init
    ne_pl(:,:,2) = n_init / ph0 / 100.
    nm_pl  = n_init
    ph_mi = 0
  end subroutine

  ! *** Implicit Function Evals ***
  subroutine phEval(g, i, j, n, m, dof, ph, b_temp)
    type(grid), intent(in) :: g
    integer, intent(in) :: i, j, n, m, dof
    real(8), intent(in) :: ph(n,m,dof)
    real(8), intent(out) :: b_temp(dof)

    call laplEqn(g, i, j, ph(:,:,1), ne_pl(:,:,1), ni_pl(:,:,1), &
                 ne_pl(:,:,2), b_temp(1))
  end subroutine

  subroutine neEval(g, i, j, n, m, dof, ne, b_temp)
    type(grid), intent(in) :: g
    integer, intent(in) :: i, j, n, m, dof
    real(8), intent(in) :: ne(n,m,dof)
    real(8), intent(out) :: b_temp(dof)
    real(8) :: dflx, src

    dflx = 0
    src = 0
    call elecDFlx(g, i, j, ph_pl(:,:,1), ne(:,:,1), ni_pl(:,:,1), &
                  ne_pl(:,:,2), dflx)
    call elecSrc(ne_mi(i,j,1), ni_pl(i,j,1), ne_pl(i,j,2), &
                 nm_pl(i,j,1), src)

    b_temp(1) = ne(i,j,1) - ne_mi(i,j,1) + g%dt * (dflx - src)
  end subroutine

  subroutine niEval(g, i, j, n, m, dof, ni, b_temp)
    type(grid), intent(in) :: g
    integer, intent(in) :: i, j, n, m, dof
    real(8), intent(in) :: ni(n,m,dof)
    real(8), intent(out) :: b_temp(dof)
    real(8) :: dflx, src

    dflx = 0
    call ionDFlx(g, i, j, ph_pl(:,:,1), ni(:,:,1), dflx)
    call ionSrc(ne_mi(i,j,1), ni_mi(i,j,1), ne_pl(i,j,2), &
                nm_pl(i,j,1), src)

    b_temp(1) = ni(i,j,1) - ni_mi(i,j,1) + g%dt * (dflx - src)
  end subroutine

  subroutine nteEval(g, i, j, n, m, dof, nte, b_temp)
    type(grid), intent(in) :: g
    integer, intent(in) :: i, j, n, m, dof
    real(8), intent(in) :: nte(n,m,dof)
    real(8), intent(out) :: b_temp(dof)
    real(8) :: dflx, src

    dflx = 0
    src = 0
    call elecEnrgDFlx(g, i, j, ph_pl(:,:,1), ne_pl(:,:,1), ni_pl(:,:,1), &
                      nte(:,:,1), dflx)
    call elecEnrgSrc(g, i, j, ph_pl(:,:,1), ne_pl(:,:,1), ni_pl(:,:,1), &
                      nte(:,:,1), nm_pl(:,:,1), src)

    b_temp(1) = nte(i,j,1) - ne_mi(i,j,2) + g%dt * (dflx - src)
  end subroutine

  subroutine nmEval(g, i, j, n, m, dof, nm, b_temp)
    type(grid), intent(in) :: g
    integer, intent(in) :: i, j, n, m, dof
    real(8), intent(in) :: nm(n,m,dof)
    real(8), intent(out) :: b_temp(dof)
    real(8) :: dflx, src

    dflx = 0
    src = 0
    call metaDFlx(g, i, j, nm(:,:,1), dflx)
    call metaSrc(ne_pl(i,j,1), ne_pl(i,j,2), nm_mi(i,j,1), src)

    b_temp(1) = nm(i,j,1) - nm_mi(i,j,1) + g%dt * (dflx - src)
  end subroutine

  ! *** Explicit Function Evals ***
  subroutine neEval_ex(g, i, j, n, m, dof, ne, b_temp)
    type(grid), intent(in) :: g
    integer, intent(in) :: i, j, n, m, dof
    real(8), intent(in) :: ne(n,m,dof)
    real(8), intent(out) :: b_temp(dof)
    real(8) :: dflx(2), src(2)

    dflx = 0
    src = 0
    call elecDFlx(g, i, j, ph_pl(:,:,1), ne(:,:,1), ni_pl(:,:,1), &
                  ne(:,:,2), dflx(1))
    call elecSrc(ne(i,j,1), ni_pl(i,j,1), ne(i,j,2), &
                 nm_pl(i,j,1), src(1))

    call elecEnrgDFlx(g, i, j, ph_pl(:,:,1), ne(:,:,1), ni_pl(:,:,1), &
                      ne(:,:,2), dflx(2))
    call elecEnrgSrc(g, i, j, ph_pl(:,:,1), ne(:,:,1), ni_pl(:,:,1), &
                      ne(:,:,2), nm_pl(:,:,1), src(2))

    b_temp(1) = src(1) - dflx(1)
    b_temp(2) = src(2) - dflx(2)
  end subroutine

  subroutine niEval_ex(g, i, j, n, m, dof, ni, b_temp)
    type(grid), intent(in) :: g
    integer, intent(in) :: i, j, n, m, dof
    real(8), intent(in) :: ni(n,m,dof)
    real(8), intent(out) :: b_temp(dof)
    real(8) :: dflx, src

    dflx = 0
    call ionDFlx(g, i, j, ph_pl(:,:,1), ni(:,:,1), dflx)
    call ionSrc(ne_pl(i,j,1), ni(i,j,1), ne_pl(i,j,2), &
                nm_pl(i,j,1), src)

    b_temp(1) = src - dflx
  end subroutine

  subroutine nmEval_ex(g, i, j, n, m, dof, nm, b_temp)
    type(grid), intent(in) :: g
    integer, intent(in) :: i, j, n, m, dof
    real(8), intent(in) :: nm(n,m,dof)
    real(8), intent(out) :: b_temp(dof)
    real(8) :: dflx, src

    dflx = 0
    src = 0
    call metaDFlx(g, i, j, nm(:,:,1), dflx)
    call metaSrc(ne_pl(i,j,1), ne_pl(i,j,2), nm(i,j,1), src)

    b_temp(1) = src - dflx
  end subroutine
end module
