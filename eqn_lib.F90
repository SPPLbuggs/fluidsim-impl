module eqn_lib
    use props
    use lapl_lib
    use elec_lib
    use ion_lib
    implicit none
    
    real(8), allocatable :: ph_pl(:,:,:), &
                            ne_pl(:,:,:), ne_mi(:,:,:), &
                            ni_pl(:,:,:), ni_mi(:,:,:)
    
    contains
    
    subroutine eqn_init(g)
    type(grid), intent(in) :: g
    
        allocate(ph_pl(g%bx+2, g%by+2, 1), &
                 ne_pl(g%bx+2, g%by+2, 1), ne_mi(g%bx+2, g%by+2, 1), &
                 ni_pl(g%bx+2, g%by+2, 1), ni_mi(g%bx+2, g%by+2, 1))
        ph_pl = 0
        ne_pl = 1
        ne_mi = 1
        ni_pl = 1
        ni_mi = 1
    end subroutine
    
    subroutine laplEval(g, i, j, n, m, dof, ph, b_temp)
        type(grid), intent(in) :: g
        integer, intent(in) :: i, j, n, m, dof
        real(8), intent(in) :: ph(n,m,dof)
        real(8), intent(out) :: b_temp(dof)
        
        call laplEqn(g, i, j, ph(:,:,1), ne_pl(:,:,1), ni_pl(:,:,1), b_temp(1))
    
    end subroutine
    
    subroutine elecEval(g, i, j, n, m, dof, ne, b_temp)
        type(grid), intent(in) :: g
        integer, intent(in) :: i, j, n, m, dof
        real(8), intent(in) :: ne(n,m,dof)
        real(8), intent(out) :: b_temp(dof)
        real(8) :: dflx_pl, dflx_mi
        
        call elecEqn(g, i, j, ph_pl(:,:,1), ne(:,:,1), dflx_pl)
        call elecEqn(g, i, j, ph_pl(:,:,1), ne_mi(:,:,1), dflx_mi)
        
        b_temp(1) = ne(i,j,1) - ne_mi(i,j,1) + 0.5 * dflx_pl + 0.5 * dflx_mi
    end subroutine
    
    subroutine ionEval(g, i, j, n, m, dof, ni, b_temp)
        type(grid), intent(in) :: g
        integer, intent(in) :: i, j, n, m, dof
        real(8), intent(in) :: ni(n,m,dof)
        real(8), intent(out) :: b_temp(dof)
        real(8) :: dflx_pl, dflx_mi
        
        call ionEqn(g, i, j, ph_pl(:,:,1), ni(:,:,1), dflx_pl)
        call ionEqn(g, i, j, ph_pl(:,:,1), ni_mi(:,:,1), dflx_mi)
        
        b_temp(1) = ni(i,j,1) - ni_mi(i,j,1) + 0.5 * dflx_pl + 0.5 * dflx_mi
    end subroutine
end module
