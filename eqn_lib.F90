module eqn_lib
    use props
    use lapl_lib
    use elec_lib
    use ion_lib
    implicit none
    
    real(8), allocatable :: ph_pl(:,:,:), ph_mi(:,:,:), &
                            ne_pl(:,:,:), ne_mi(:,:,:), &
                            ni_pl(:,:,:), ni_mi(:,:,:), &
                            nte_pl(:,:,:), nte_mi(:,:,:), &
                            nm_pl(:,:,:), nm_mi(:,:,:)
    
    contains
    
    subroutine eqn_init(g)
    type(grid), intent(in) :: g
        allocate(ph_pl(g%bx+2, g%by+2, 1), ph_mi(g%bx+2, g%by+2, 1), &
                 ne_pl(g%bx+2, g%by+2, 1), ne_mi(g%bx+2, g%by+2, 1), &
                 ni_pl(g%bx+2, g%by+2, 1), ni_mi(g%bx+2, g%by+2, 1), &
                 nte_pl(g%bx+2, g%by+2, 1), nte_mi(g%bx+2, g%by+2, 1), &
                 nm_pl(g%bx+2, g%by+2, 1), nm_mi(g%bx+2, g%by+2, 1))
        
        ph_pl  = 0
        ne_pl  = n_init
        ni_pl  = n_init
        nte_pl = n_init / ph0 / 100.
        nm_pl  = n_init
    end subroutine
    
    subroutine phEval(g, i, j, n, m, dof, ph, b_temp)
        type(grid), intent(in) :: g
        integer, intent(in) :: i, j, n, m, dof
        real(8), intent(in) :: ph(n,m,dof)
        real(8), intent(out) :: b_temp(dof)
        
        call laplEqn(g, i, j, ph(:,:,1), ne_pl(:,:,1), ni_pl(:,:,1), &
                     nte_pl(:,:,1), b_temp(1))
    end subroutine
    
    subroutine neEval(g, i, j, n, m, dof, ne, b_temp)
        type(grid), intent(in) :: g
        integer, intent(in) :: i, j, n, m, dof
        real(8), intent(in) :: ne(n,m,dof)
        real(8), intent(out) :: b_temp(dof)
        
        call elecEqn(g, i, j, ph_pl(:,:,1), ni_pl(:,:,1), ne(:,:,1), &
                     nte_pl(:,:,1), nm_pl(:,:,1), b_temp(1))
    end subroutine
    
    subroutine niEval(g, i, j, n, m, dof, ni, b_temp)
        type(grid), intent(in) :: g
        integer, intent(in) :: i, j, n, m, dof
        real(8), intent(in) :: ni(n,m,dof)
        real(8), intent(out) :: b_temp(dof)
        
        call ionEqn(g, i, j, ph_pl(:,:,1), ni(:,:,1), ne_pl(:,:,1), &
                    nte_pl(:,:,1), nm_pl(:,:,1), b_temp(1))
    end subroutine
    
    subroutine nteEval(g, i, j, n, m, dof, nte, b_temp)
        type(grid), intent(in) :: g
        integer, intent(in) :: i, j, n, m, dof
        real(8), intent(in) :: nte(n,m,dof)
        real(8), intent(out) :: b_temp(dof)
        
        call elecEnrgEqn(g, i, j, ph_pl(:,:,1), ne_pl(:,:,1), ni_pl(:,:,1), &
                         nte(:,:,1), nm_pl(:,:,1), b_temp(1))
    end subroutine
    
    subroutine nmEval(g, i, j, n, m, dof, nm, b_temp)
        type(grid), intent(in) :: g
        integer, intent(in) :: i, j, n, m, dof
        real(8), intent(in) :: nm(n,m,dof)
        real(8), intent(out) :: b_temp(dof)
        
        call metaEqn(g, i, j, nm(:,:,1), ne_pl(:,:,1), nte_pl(:,:,1), b_temp(1))
    end subroutine
end module















