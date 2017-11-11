module rk_lib
    use props
    implicit none
    
    real(8) :: err_mi = 1.0
    real(8), allocatable :: f_or(:,:,:), k(:,:,:,:)
    
    
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
! *** Runge-Kutta Initialization ***
    subroutine rk_init(g)
        type(grid), intent(in) :: g
        
        allocate(f_or(g%bx+2, g%by+2, 5), &
                 k(g%bx+2, g%by+2, 5, 5))
        
    end subroutine

! *** Runge-Kutta Timestepping ***
    subroutine rk_step(g, order, dof, f_pl, f_mi, fEval, f_min)
        type(grid), intent(inout) :: g
        integer, intent(in)    :: order, dof
        real(8), intent(in)    :: f_min(:)
        real(8), intent(inout) :: f_mi(:,:,:), f_pl(:,:,:)
        procedure(subIn)       :: feval
        integer :: stage, i, j
        real(8) :: scfac, err_pl, err_n(dof)
        
        stage = 0
        f_or(:,:,1:dof) = f_mi
        
        do
            stage = stage + 1
            if (stage > order) exit
            
            do j = 2, g%by+1
                do i = 2, g%bx+1
                    call fEval(g, i, j, g%bx+2, g%by+2, dof, f_mi, k(i,j,1:dof,stage))
                end do
            end do
            
            err_pl = 0
            call rk_int(g, dof, f_pl, f_mi, f_or(:,:,1:dof), k(:,:,1:dof,:), &
                        stage, g%dt, err_n, order, f_min)
            
            if ((order == 5) .and. (stage == 5)) then
                err_pl = sqrt(sum(err_n**2))
                scfac = 0.8 * err_pl**(-0.7 / 4.) &
                        * err_mi**( 0.4 / 4.)
                scfac = min(2.5d0, max(3d-1, scfac))
                g%dt = scfac * g%dt
                
                g%dt = min(max(g%dt, 1d-12), 2d-3)
                
                call MPI_Bcast(g%dt, 1, etype, 0, comm, ierr)
               
                if (g%dt <= 1.1d-12) then
                    write(*,*) 'minimum timestep reached; finishing simulation'
                    write(*,'(es10.2)') err_pl
                    stop
                end if
                
                if (err_pl .le. 1.0) then
                    err_mi = err_pl
                else
                    stage = 0
                    
                    f_pl = f_or(:,:,1:dof)
                    f_mi = f_or(:,:,1:dof)
                end if
            end if
        end do
    end subroutine

    subroutine rk_int(g, dof, f_pl, f_mi, f_or, k, stage, dt, nerr_n, order, f_min)
        type(grid), intent(in) :: g
        integer, intent(in) :: dof, stage, order
        real(8), intent(in) :: k(:,:,:,:), f_min(:), dt
        real(8), intent(inout) :: f_pl(:,:,:), f_mi(:,:,:), f_or(:,:,:), nerr_n(:)
        real(8) :: err_n(g%bx+2, g%by+2), abs_tol = 1e-4, rel_tol = 1e-4
        integer :: i,j,d
        
        ! Explicit Euler
        if (order == 1) then
            do d = 1, dof
                do j = 2, g%by+1
                    do i = 2, g%bx+1
                        f_pl(i,j,d) = f_or(i,j,d) + k(i,j,d,1) * dt
                    end do
                end do
                
                call comm_real(g%bx, g%by, f_pl(:,:,d))
            end do
            
        ! 2nd order
        else if (order == 2) then
            if (stage == 1) then
                do d = 1, dof
                    do j = 2, g%by+1
                        do i = 2, g%bx+1
                            f_pl(i,j,d) = f_or(i,j,d) + k(i,j,d,1) * dt / 2.0
                            f_mi(i,j,d) = f_or(i,j,d) + k(i,j,d,1) * dt
                        end do
                    end do
                    
                    f_pl(:,:,d) = max(f_pl(:,:,d), f_min(d))
                    f_mi(:,:,d) = max(f_mi(:,:,d), f_min(d))
                    
                    call comm_real(g%bx, g%by, f_pl(:,:,d))
                    call comm_real(g%bx, g%by, f_mi(:,:,d))
                end do
            else
                do d = 1, dof
                    do j = 2, g%by+1
                        do i = 2, g%bx+1
                            f_pl(i,j,d) = f_pl(i,j,d) + k(i,j,d,2) * dt / 2.0
                        end do
                    end do
                    
                    f_pl(:,:,d) = max(f_pl(:,:,d), f_min(d))
                    
                    call comm_real(g%bx, g%by, f_pl(:,:,d))
                end do
            end if
        
        ! 4th order
        else if (order == 4) then
            if (stage == 1) then
                do d = 1, dof
                    do j = 2, g%by+1
                        do i = 2, g%bx+1
                            f_pl(i,j,d) = f_or(i,j,d) + k(i,j,d,1) * dt / 6.0
                        end do
                    end do
                
                    f_pl(:,:,d) = max(f_pl(:,:,d), f_min(d))
                    call comm_real(g%bx, g%by, f_pl(:,:,d))
                end do
            else if (stage == 2) then
                do d = 1, dof
                    do j = 2, g%by+1
                        do i = 2, g%bx+1
                            f_pl(i,j,d) = f_pl(i,j,d) + k(i,j,d,2) * dt / 3.0
                            f_mi(i,j,d) = f_or(i,j,d) + k(i,j,d,2) * dt / 2.0
                        end do
                    end do
                    
                    f_pl(:,:,d) = max(f_pl(:,:,d), f_min(d))
                    f_mi(:,:,d) = max(f_mi(:,:,d), f_min(d))
                    
                    call comm_real(g%bx, g%by, f_pl(:,:,d))
                    call comm_real(g%bx, g%by, f_mi(:,:,d))
                end do
            else if (stage == 3) then
                do d = 1, dof
                    do j = 2, g%by+1
                        do i = 2, g%bx+1
                            f_pl(i,j,d) = f_pl(i,j,d) + k(i,j,d,3) * dt / 3.0
                            f_mi(i,j,d) = f_or(i,j,d) + k(i,j,d,3) * dt
                        end do
                    end do
                    
                    f_pl(:,:,d) = max(f_pl(:,:,d), f_min(d))
                    f_mi(:,:,d) = max(f_mi(:,:,d), f_min(d))
                    
                    call comm_real(g%bx, g%by, f_pl(:,:,d))
                    call comm_real(g%bx, g%by, f_mi(:,:,d))
                end do
            else
                do d = 1, dof
                    do j = 2, g%by+1
                        do i = 2, g%bx+1
                            f_pl(i,j,d) = f_pl(i,j,d) + k(i,j,d,4) * dt / 6.0
                        end do
                    end do
                    
                    f_pl(:,:,d) = max(f_pl(:,:,d), f_min(d))
                    call comm_real(g%bx, g%by, f_pl(:,:,d))
                end do
            end if
            
        ! merson 4("5") adaptive time-stepping
        else if (order == 5) then
            if (stage == 1) then
                do d = 1, dof
                    do j = 2, g%by+1
                        do i = 2, g%bx+1
                            f_mi(i,j,d) = f_or(i,j,d) + k(i,j,d,1) * dt / 3.0
                        end do
                    end do
                    
                    f_mi(:,:,d) = max(f_mi(:,:,d), f_min(d))
                    call comm_real(g%bx, g%by, f_mi(:,:,d))
                end do
                
            else if (stage == 2) then
                do d = 1, dof
                    do j = 2, g%by+1
                        do i = 2, g%bx+1
                            f_mi(i,j,d) = f_or(i,j,d) + dt * ( k(i,j,d,1) + k(i,j,d,2) ) / 6.0
                         end do
                    end do
                       
                    f_mi(:,:,d) = max(f_mi(:,:,d), f_min(d))
                    
                    call comm_real(g%bx, g%by, f_mi(:,:,d))
                end do
            else if (stage == 3) then
                do d = 1, dof
                    do j = 2, g%by+1
                        do i = 2, g%bx+1
                            f_mi(i,j,d) = f_or(i,j,d) + dt * (k(i,j,d,1) + k(i,j,d,3) * 3.0) / 8.0
                        end do
                    end do
                    
                    f_mi(:,:,d) = max(f_mi(:,:,d), f_min(d))
                    call comm_real(g%bx, g%by, f_mi(:,:,d))
                end do
            else if (stage == 4) then
                do d = 1, dof
                    do j = 2, g%by+1
                        do i = 2, g%bx+1
                            f_mi(i,j,d) = f_or(i,j,d) + dt * (k(i,j,d,1) &
                                        - k(i,j,d,3) * 3.0 + k(i,j,d,4) * 4.0 ) / 2.0
                        end do
                    end do
                    
                    f_mi(:,:,d) = max(f_mi(:,:,d), f_min(d))
                    call comm_real(g%bx, g%by, f_mi(:,:,d))
                end do
            else
                do d = 1, dof
                    do j = 2, g%by+1
                        do i = 2, g%bx+1
                            f_pl(i,j,d) = f_or(i,j,d) + dt * (k(i,j,d,1) &
                                        + k(i,j,d,4) * 4.0 + k(i,j,d,5)) / 6.0
                        end do
                    end do
                    
                    f_pl(:,:,d) = max(f_pl(:,:,d), f_min(d))
                    call comm_real(g%bx, g%by, f_pl(:,:,d))
                    
                
                    err_n = 0
                    
                    err_n = abs(dt * (k(:,:,d,1) * 2.0 / 30.0 - k(:,:,d,3) * 3.0 / 10.0 &
                                + k(:,:,d,4) * 4.0 / 15.0 - k(:,:,d,5) / 30.0 ))
                    
                    nerr_n(d) = maxval(err_n/(abs_tol+rel_tol*abs(f_or(:,:,d))))
                    
                    call MPI_Allreduce(MPI_In_Place, nerr_n(d), 1, etype, &
                                       MPI_Max, comm, ierr)
               end do
            end if
        end if
    end subroutine
end module
