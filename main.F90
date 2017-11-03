    program main
    use props
    use petsc_lib
    use eqn_lib
    implicit none
    
    type(grid) :: g
    integer :: i, iter, ts = 0, nx, ny, dof, narg, conv
    real(8) :: l, w, dt, t_fin, t_pr, t_sv, sim_start, time, relErr
    character(80):: arg, path
    
    call PetscInitialize(petsc_null_character, ierr)
    comm = PETSC_COMM_WORLD
    
    call MPI_Comm_rank(comm, my_id, ierr)
    call MPI_Comm_size(comm, nproc, ierr)
    
    call cpu_time(sim_start)
    
    ! Default properties
    nx = 50
    ny = 1
    px = 1
    py = 1
    dof = 1
    l  = 1
    w  = 1
    dt = 1e-2
    t_fin = 10
    t_pr = t_fin/100.
    t_sv = t_fin/100.
    
    ! Check for -help
    narg = iargc()
    if (mod(narg,2) .ne. 0) then
        if (my_id == 0) then
            write(*,*)
            write(*,*) 'Usage:   mpiexec -n <nproc> ./main <options>'
            write(*,*) 'Options: -nx <nx>, -ny <ny>, -px <px>, -py <py>'
            write(*,*) '         -l <l>, -w <w>, -t <t>, -unif <T/F>'
            write(*,*)
        end if
        call MPI_Finalize(ierr)
        stop
    end if
    
    ! Read input arguments
    do i = 1, narg/2
        call getarg(2 * (i - 1) + 1, arg)
        select case (arg)
            case ('-nx')
                call getarg(2 * (i - 1) + 2, arg)
                read(arg,*) nx
            case ('-ny')
                call getarg(2 * (i - 1) + 2, arg)
                read(arg,*) ny
            case ('-px')
                call getarg(2 * (i - 1) + 2, arg)
                read(arg,*) px
            case ('-py')
                call getarg(2 * (i - 1) + 2, arg)
                read(arg,*) py
            case ('-l')
                call getarg(2 * (i - 1) + 2, arg)
                read(arg,*) l
            case ('-w')
                call getarg(2 * (i - 1) + 2, arg)
                read(arg,*) w
            case ('-t')
                call getarg(2 * (i - 1) + 2, arg)
                read(arg,*) t_fin
            case ('-dt')
                call getarg(2 * (i - 1) + 2, arg)
                read(arg,*) dt
            case ('-unif')
                call getarg(2 * (i - 1) + 2, arg)
                read(arg,*) unif
        end select
    end do
    
    if (px * py .ne. nproc) then
        if (my_id == 0) then
            write(*,*) 'Error: px * py must equal nproc. Stop.'
            write(*,*) '       Enter ./main -help for usage.'
            call MPI_Abort(comm,2,ierr)
        end if
    end if
    if (mod(nx,px) .ne. 0) then
        if (my_id == 0) then
            write(*,*) 'Error: px needs to divide nx. Stop.'
            write(*,*) '       Enter ./main -help for usage.'
            call MPI_Abort(comm,3,ierr)
        end if
    end if
    if (mod(ny,py) .ne. 0) then
        if (my_id == 0) then
            write(*,*) 'Error: py needs to divide ny. Stop.'
            write(*,*) '       Enter ./main -help for usage.'
            call MPI_Abort(comm,4,ierr)
        end if
    end if
    
    path = 'Output/'
    call g_init(g, nx, ny, px, py, dof, l, w, trim(path))
    
    allocate(f(g%bx+2, g%by+2, 2), f_mi(g%bx+2, g%by+2, 2))
    f(:,:,1) = 0
    f(:,:,2) = 1
    f_mi = f
    
    g%t  = 0
    g%dt = dt
    
    ! Create PETSc Objects
    call petsc_create(g, A1, b1, x1)
    call petsc_create(g, A2, b2, x2)
    
    do
        ts = ts + 1
        g%t = g%t + g%dt
        if (g%t >= t_fin) exit
        
        ! Update boundary conditions
        if (rx == 0) f(1,:,1) = 1!sin(2.0 * 3.14159 * g%t / 10.0)
        if (ry == 0) f(:,1,1) = 1!sin(2.0 * 3.14159 * g%t / 10.0)

        ! Solve first equation
        do iter = 1, 100
            ! Assemble jacobian and RHS
            call assem_Ab(g, A1, b1, f(:,:,1:1), feval1)
            if (ts == 1) then
                call MatSetOption(A1, Mat_New_Nonzero_Locations, PETSc_False, ierr)
            end if
            
            !call view(A1, b1)
            
            ! Check norm of residual
            call VecNorm(b1, norm_2, relErr, ierr)
            if (relErr < 1d-8) exit
            
            ! Solve system:
            call KSPSetOperators(ksp, A1, A1, ierr)
            call KSPSolve(ksp, b1, x1, ierr)
            
            call KSPGetConvergedReason(ksp, conv, ierr)
            if ((my_id == 0) .and. (conv .ne. 2)) write(*,3) conv
            3 format('KSP did not converge; Reason = ',i0)
            
            ! Update variables with solution:
            call upd_soln(g, x1, f(:,:,1:1))
        end do
        
        ! Solve second equation
        do iter = 1, 100
            ! Assemble jacobian and RHS
            call assem_Ab(g, A2, b2, f(:,:,2:2), feval2)
            if (ts == 1) then
                call MatSetOption(A2, Mat_New_Nonzero_Locations, PETSc_False, ierr)
            end if
            
            !call view(A2, b2)
            
            ! Check norm of residual
            call VecNorm(b2, norm_2, relErr, ierr)
            if (relErr < 1d-8) exit
            
            ! Solve system:
            call KSPSetOperators(ksp, A2, A2, ierr)
            call KSPSolve(ksp, b2, x2, ierr)
            
            call KSPGetConvergedReason(ksp, conv, ierr)
            if ((my_id == 0) .and. (conv .ne. 2)) write(*,3) conv
            
            ! Update variables with solution:
            call upd_soln(g, x2, f(:,:,2:2))
        end do
        
        f_mi = f
        
        if ((t_pr <= g%t) .and. (my_id == 0)) then
            call cpu_time(time)
            write(*,*)
            write(*,11) float(ts), g%t, g%dt, time - sim_start
            t_pr = t_pr + t_fin/100.
        end if
    
        if (t_sv <= g%t) then
            call savedat(trim(path)//'f1.dat', f(:,:,1))
            call savedat(trim(path)//'f2.dat', f(:,:,2))
            
            call MPI_File_Open(comm, trim(path)//'time.dat', &
                MPI_MODE_WRonly + MPI_Mode_Append,  info, fh, ierr)
            if (my_id == 0) call MPI_File_Write(fh, g%t, 1, etype, stat, ierr)
            call MPI_File_Close(fh, ierr)
            
            t_sv  = t_sv + t_fin/100.
        end if
    end do
    
    if (my_id == 0) then
        call cpu_time(time)
        write(*,*)
        write(*,9) int(time - sim_start) / 3600, &
                    mod(int(time - sim_start)/60,60)
        write(*,*)
    end if
    
    call petsc_destroy(A1, b1, x1)
    call PetscFinalize(ierr)

11 format('Timestep:', es9.2, '  Time:', es9.2, '  dT:', es9.2, '  Dur:', f6.2, ' sec')
9  format('Simulation finished in ', i0, ' hr ', i0, ' min')

    end program
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
