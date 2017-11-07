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
        allocate(f_pl(g%bx+2, g%by+2, 5), f_mi(g%bx+2, g%by+2, 5))
        f_pl = n_init
        f_pl(:,:,1) = 0
        f_pl(:,:,4) = n_init / ph0 / 100.
    end subroutine
    
    subroutine fEval(g, i, j, n, m, dof, f, b_temp)
        type(grid), intent(in) :: g
        integer, intent(in) :: i, j, n, m, dof
        real(8), intent(in) :: f(n,m,dof)
        real(8), intent(out) :: b_temp(dof)
        real(8) :: src(3), dflxe_pl, dflxi_pl, dflxt_pl, dflxm_pl, &
                   dflxe_mi, dflxi_mi, dflxt_mi, dflxm_mi
        
        call laplEqn(g, i, j, f(:,:,1), f(:,:,2), f(:,:,3), b_temp(1))
        
        call elecEqn(g, i, j, f(:,:,1), f(:,:,2), f_mi(:,:,2), f_mi(:,:,4), f_mi(:,:,3), dflxe_pl)
        call elecEnrgEqn(g, i, j, f(:,:,1), f(:,:,4), f_mi(:,:,2), f_mi(:,:,4), f_mi(:,:,3), dflxt_pl)
        call ionEqn(g, i, j, f(:,:,1), f(:,:,3), dflxi_pl)
        call metaEqn(g, i, j, f(:,:,5), dflxm_pl)
        
        !call elecEqn(g, i, j, f_mi(:,:,1), f_mi(:,:,2), f_mi(:,:,2), f_mi(:,:,4), f_mi(:,:,3), dflxe_mi)
        !call elecEnrgEqn(g, i, j, f_mi(:,:,1), f_mi(:,:,4), f_mi(:,:,2), f_mi(:,:,4), f_mi(:,:,3), dflxt_mi)
        !call ionEqn(g, i, j, f_mi(:,:,1), f_mi(:,:,3), dflxi_mi)
        !call metaEqn(g, i, j, f_mi(:,:,5), dflxm_mi)
        
        call calc_src(g, i, j, f_mi(:,:,1), f_mi(:,:,2), f_mi(:,:,3), f_mi(:,:,4), f_mi(:,:,5), src)
        
        b_temp(2) = f(i,j,2) - f_mi(i,j,2) - src(1) + dflxe_pl !+ 0.5 * (dflxe_pl + dflxe_mi)
        b_temp(3) = f(i,j,3) - f_mi(i,j,3) - src(1) + dflxi_pl !+ 0.5 * (dflxi_pl + dflxi_mi)
        b_temp(4) = f(i,j,4) - f_mi(i,j,4) - src(2) + dflxt_pl !+ 0.5 * (dflxt_pl + dflxt_mi)
        b_temp(5) = f(i,j,5) - f_mi(i,j,5) - src(3) + dflxm_pl !+ 0.5 * (dflxm_pl + dflxm_mi)
    end subroutine
end module















