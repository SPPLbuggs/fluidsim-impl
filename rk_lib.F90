module rk_lib
    use props
    implicit none
    
    abstract interface
        subroutine subIn(g, i, j, n, m, dof, f, b_temp)
            use props
            type(grid), intent(in) :: g
            integer, intent(in) :: i, j, n, m, dof
            real(8), intent(in) :: f(n,m,dof)
            real(8), intent(out) :: b_temp(dof)
        end subroutine
    end interface
    
contains
! *** Runge-Kutta Timestepping ***
    subroutine rk_step(g, order, f_pl, f_mi, fEval, f_min)
        type(grid), intent(in) :: g
        integer, intent(in)    :: order
        real(8), intent(in)    :: f_min
        real(8), intent(inout) :: f_mi(:,:,:), f_pl(:,:,:)
        procedure(subIn)       :: feval
        integer :: stage, i, j
        real(8), allocatable :: f_or(:,:,:), k(:,:,:)
        real(8) :: err_cur = 1.0
        
        allocate(f_or(g%bx+2, g%by+2,1), k(g%bx+2, g%by+2, order))
        
        stage = 0
        f_or = f_mi
        
        do
            stage = stage + 1
            if (stage > order) exit
            
            do j = 2, g%by+1
                do i = 2, g%bx+1
                    call fEval(g, i, j, g%bx+2, g%by+2, 1, f_mi, k(i,j,stage))
                end do
            end do
            
            call rk_int(g, f_pl(:,:,1), f_mi(:,:,1), f_or(:,:,1), k, stage, g%dt, err_cur, order, f_min)
            
        end do
    end subroutine

    subroutine rk_int(g, f_pl, f_mi, f_or, k, stage, dt, nerr_n, order, f_min)
        type(grid), intent(in) :: g
        integer, intent(in) :: stage, order
        real(8), intent(in) :: k(:,:,:), f_min, dt
        real(8), intent(inout) :: f_pl(:,:), f_mi(:,:), f_or(:,:), nerr_n
        real(8) :: err_n(g%bx+2, g%by+2), abs_tol = 1e-4, rel_tol = 1e-4
        integer :: i,j
        
        if (order == 1) then
        ! euler scheme
            do j = 2, g%by+1
                do i = 2, g%bx+1
                    f_pl(i,j) = f_or(i,j) + k(i,j,1) * dt
                end do
            end do
            
            call comm_real(g%bx, g%by, f_pl)
            
        else if (order == 2) then
        ! 2nd order
            if (stage == 1) then
                do j = 2, g%by+1
                    do i = 2, g%bx+1
                        f_pl(i,j) = f_or(i,j) + k(i,j,1) * dt / 2.0
                        f_mi(i,j) = f_or(i,j) + k(i,j,1) * dt
                    end do
                end do
                
            else
                do j = 2, g%by+1
                    do i = 2, g%bx+1
                        f_pl(i,j) = f_pl(i,j) + k(i,j,2) * dt / 2.0
                    end do
                end do
            end if
            
            f_pl = max(f_pl, f_min)
            f_mi = max(f_mi, f_min)
            
            call comm_real(g%bx, g%by, f_pl)
            call comm_real(g%bx, g%by, f_mi)
        
        else if (order == 4) then
            if (stage == 1) then
                do j = 2, g%by+1
                    do i = 2, g%bx+1
                        f_pl(i,j) = f_or(i,j) + k(i,j,1) * dt / 6.0
                    end do
                end do
                
                f_pl = max(f_pl, f_min)
                call comm_real(g%bx, g%by, f_pl)
            else if (stage == 2) then
                do j = 2, g%by+1
                    do i = 2, g%bx+1
                        f_pl(i,j) = f_pl(i,j) + k(i,j,2) * dt / 3.0
                        f_mi(i,j) = f_or(i,j) + k(i,j,2) * dt / 2.0
                    end do
                end do
                
                f_pl = max(f_pl, f_min)
                f_mi = max(f_mi, f_min)
                
                call comm_real(g%bx, g%by, f_pl)
                call comm_real(g%bx, g%by, f_mi)
            else if (stage == 3) then
                do j = 2, g%by+1
                    do i = 2, g%bx+1
                        f_pl(i,j) = f_pl(i,j) + k(i,j,3) * dt / 3.0
                        f_mi(i,j) = f_or(i,j) + k(i,j,3) * dt
                    end do
                end do
                
                f_pl = max(f_pl, f_min)
                f_mi = max(f_mi, f_min)
                
                call comm_real(g%bx, g%by, f_pl)
                call comm_real(g%bx, g%by, f_mi)
            else
                do j = 2, g%by+1
                    do i = 2, g%bx+1
                        f_pl(i,j) = f_pl(i,j) + k(i,j,4) * dt / 6.0
                    end do
                end do
                
                f_pl = max(f_pl, f_min)
                call comm_real(g%bx, g%by, f_pl)
            end if
            
        ! merson 4("5") adaptive time-stepping
        else if (order == 5) then
            if (stage == 1) then
                do j = 2, g%by+1
                    do i = 2, g%bx+1
                        f_mi(i,j) = f_or(i,j) + k(i,j,1) * dt / 3.0
                    end do
                end do
                
                f_mi(:,:) = max(f_mi, f_min)
                call comm_real(g%bx, g%by, f_mi)
                
            else if (stage == 2) then
                do j = 2, g%by+1
                    do i = 2, g%bx+1
                        f_mi(i,j) = f_or(i,j) + dt * ( k(i,j,1) + k(i,j,2) ) / 6.0
                     end do
                end do
                   
                f_mi = max(f_mi, f_min)
                
                call comm_real(g%bx, g%by, f_mi)
                
            else if (stage == 3) then
                do j = 2, g%by+1
                    do i = 2, g%bx+1
                        f_mi(i,j) = f_or(i,j) + dt * (k(i,j,1) + k(i,j,3) * 3.0) / 8.0
                    end do
                end do
                
                f_mi = max(f_mi, f_min)
                call comm_real(g%bx, g%by, f_mi)
                
            else if (stage == 4) then
                do j = 2, g%by+1
                    do i = 2, g%bx+1
                        f_mi(i,j) = f_or(i,j) + dt * (k(i,j,1) &
                                    - k(i,j,3) * 3.0 + k(i,j,4) * 4.0 ) / 2.0
                    end do
                end do
                
                f_mi = max(f_mi, f_min)
                call comm_real(g%bx, g%by, f_mi)
            else
                do j = 2, g%by+1
                    do i = 2, g%bx+1
                        f_pl(i,j) = f_or(i,j) + dt * (k(i,j,1) &
                                    + k(i,j,4) * 4.0 + k(i,j,5)) / 6.0
                    end do
                end do
                
                f_pl = max(f_pl, f_min)
                call comm_real(g%bx, g%by, f_pl)
                
                err_n = 0
                
                err_n = abs(dt * (k(:,:,1) * 2.0 / 30.0 - k(:,:,3) * 3.0 / 10.0 &
                            + k(:,:,4) * 4.0 / 15.0 - k(:,:,5) / 30.0 ))
                
                nerr_n = maxval(err_n/(abs_tol+rel_tol*abs(f_or)))
                
                call MPI_Allreduce(MPI_In_Place, nerr_n, 1, etype, &
                                   MPI_Max, comm, ierr)
            end if
        end if
    end subroutine
end module
