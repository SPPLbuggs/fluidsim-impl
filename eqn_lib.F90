module eqn_lib
    use props
    use lapl_lib
    use elec_lib
    use ion_lib
    implicit none
    
    real(8), allocatable :: f_pl(:,:,:), f_mi(:,:,:)
    
    contains
    
    subroutine eqn_init(g)
    type(grid), intent(in) :: g
        allocate(f_pl(g%bx+2, g%by+2, 3), f_mi(g%bx+2, g%by+2, 3))
        f_pl = n_init
        f_pl(:,:,1) = 0
        
        if (rx == 0) then
            f_pl(1,:,2:3) = n_zero
        end if
        
        if (rx == px - 1) then
            f_pl(g%bx+2,:,2:3) = n_zero
        end if
        
        if (ry == 0) then
            f_pl(:,1,2:3) = n_zero
        end if
        
        if (ry == py - 1) then
            f_pl(:,g%by+2,2:3) = n_zero
        end if
    end subroutine
    
    subroutine fEval(g, i, j, n, m, dof, f, b_temp)
        type(grid), intent(in) :: g
        integer, intent(in) :: i, j, n, m, dof
        real(8), intent(in) :: f(n,m,dof)
        real(8), intent(out) :: b_temp(dof)
        real(8) :: src, dflxe_pl, dflxi_pl!, dflxe_mi, dflxi_mi
        
        call laplEqn(g, i, j, f(:,:,1), f(:,:,2), f(:,:,3), b_temp(1))
        call elecEqn(g, i, j, f(:,:,1), f(:,:,2), dflxe_pl)
        call ionEqn(g, i, j, f(:,:,1), f(:,:,3), dflxi_pl)
        !call elecEqn(g, i, j, f_mi(:,:,1), f_mi(:,:,2), dflxe_mi)
        !call ionEqn(g, i, j, f_mi(:,:,1), f_mi(:,:,3), dflxi_mi)
        call calc_src(g, i, j, f_mi(:,:,1), f_mi(:,:,2), src)
        
        b_temp(2) = f(i,j,2) - f_mi(i,j,2) - src + dflxe_pl
        b_temp(3) = f(i,j,3) - f_mi(i,j,3) - src + dflxi_pl
    end subroutine
end module
