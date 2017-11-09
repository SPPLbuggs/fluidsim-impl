program main
    use props
    use petsc_lib
    use eqn_lib
    use rk_lib
    use circ_lib
    implicit none
    
    type(grid) :: g
    integer :: ts = 0, ts1, ts2, nx, ny, dof
    real(8) :: l, w, ew, vl, r0, dt, t_fin, t_pr, t_sv, t_sv0, &
               sim_start, time1, time2
    character(80):: path
    logical :: assem(5) = .True.
    
    ! Initialize PETSc and MPI
    call PetscInitialize(petsc_null_character, ierr)
    comm = PETSC_COMM_WORLD
    call MPI_Comm_rank(comm, my_id, ierr)
    call MPI_Comm_size(comm, nproc, ierr)
    
    call cpu_time(sim_start)
    call cpu_time(time1)
    ts1 = 0
    
    ! Default properties
    nx = 100
    ny = 1
    px = 1
    py = 1
    dof = 1
    l  = 1e-2 / x0
    w  = 1e-2 / x0
    ew = 2e-2 / x0
    dt = 1e-4
    t_fin = 10
    t_pr = 1e-2
    t_sv = 1e-3
    t_sv0 = 1e-3
    vl = 500 / ph0
    r0 = 1e6
    
    ! Read input arguments
    call read_in
    
    ! Initialize grid and arrays
    path = 'Output/'
    call g_init(g, nx, ny, px, py, dof, l, w, ew, trim(path))
    call eqn_init(g)
    call circ_init(vl, r0)
    
    g%t  = 0
    g%dt = dt
    
    do
        ts = ts + 1
        g%dt = min(g%dt*1.001, 2d-3)
        g%t = g%t + g%dt
        if (g%t >= t_fin) exit
        
        ! Update boundary conditions
        if (rx == 0) ph_pl(1,:,1) = Vd_pl
        if (ry == 0) ph_pl(:,1,1) = Vd_pl
        
        ! Solve ne system
        t_m = 1e9
        ne_mi = ne_pl
        call petsc_step(g, A2, b2, x2, ne_pl, neEval, (/ n_zero /), assem(2))
        
        ! Solve ni system
        ni_mi = ni_pl
        call petsc_step(g, A3, b3, x3, ni_pl, niEval, (/ n_zero /), assem(3))
        
        ! Solve ni system
        nte_mi = nte_pl
        call petsc_step(g, A4, b4, x4, nte_pl, nteEval, (/ n_zero / ph0 / 100. /), assem(4))
        
        ! Solve nm system
        nm_mi = nm_pl
        call petsc_step(g, A5, b5, x5, nm_pl, nmEval, (/ n_zero /), assem(5))
        
        ! Solve ph, ne, ni, nte system
        ph_mi = ph_pl
        call petsc_step(g, A1, b1, x1, ph_pl, phEval, (/ -1d0 /), assem(1))
        
        ! Solve external circuit system
        call circ_step(g)
        
        ! Print out some information
        if ((t_pr <= g%t) .and. (my_id == 0)) then
            call cpu_time(time2)
            ts2 = ts
            write(*,*)
            write(*,11) float(ts), g%t, (time2 - time1) / g%dt / float(ts2-ts1) / 60.
            write(*,12)  g%dt, t_m
            write(*,13)  Vd_pl * ph0, Id * e / t0
            t_pr = t_pr + 1e-2
            call cpu_time(time1)
            ts1 = ts
        end if
        
        ! Save data
        if (t_sv <= g%t) then
            call savedat(trim(path)//'f1.dat', ph_pl(:,:,1) * ph0)
            call savedat(trim(path)//'f2.dat', ne_pl(:,:,1) / x0**3)
            call savedat(trim(path)//'f3.dat', ni_pl(:,:,1) / x0**3)
            call savedat(trim(path)//'f4.dat', nte_pl(:,:,1) / x0**3 * ph0 / 1.5)
            call savedat(trim(path)//'f5.dat', nm_pl(:,:,1) / x0**3)
            
            call MPI_File_Open(comm, trim(path)//'time.dat', &
                MPI_MODE_WRonly + MPI_Mode_Append,  info, fh, ierr)
            if (my_id == 0) call MPI_File_Write(fh, g%t, 1, etype, stat, ierr)
            call MPI_File_Close(fh, ierr)
            
            call MPI_File_Open(comm, trim(path)//'vd.dat', &
                MPI_MODE_WRonly + MPI_Mode_Append,  info, fh, ierr)
            if (my_id == 0) call MPI_File_Write(fh, Vd_pl*ph0, 1, etype, stat, ierr)
            call MPI_File_Close(fh, ierr)
            
            call MPI_File_Open(comm, trim(path)//'id.dat', &
                MPI_MODE_WRonly + MPI_Mode_Append,  info, fh, ierr)
            if (my_id == 0) call MPI_File_Write(fh, Id*e/t0, 1, etype, stat, ierr)
            call MPI_File_Close(fh, ierr)
            
            t_sv  = t_sv + t_sv0
            t_sv0 = t_sv0 * 1.01
        end if
    end do
    
    if (my_id == 0) then
        call cpu_time(time1)
        write(*,*)
        write(*,9) int(time1 - sim_start) / 3600, &
                    mod(int(time1 - sim_start)/60,60)
        write(*,*)
    end if
    
    call petsc_destroy(A1, b1, x1)
    call petsc_destroy(A2, b2, x2)
    call petsc_destroy(A3, b3, x3)
    call petsc_destroy(A4, b4, x4)
    call petsc_destroy(A5, b5, x5)
    call PetscFinalize(ierr)

11 format('Timestep:', es9.2, '  Time:', es9.2, '  time/us:', f6.2, ' min')
12 format('  dT:', es9.2, '  tm:', es9.2)
13 format('  Vd:', es9.2, '  Id:', es9.2)
9  format('Simulation finished in ', i0, ' hr ', i0, ' min')
    
contains
    
    subroutine read_in
        integer :: i, narg
        character(80) :: arg
        
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
                    l = l / x0
                case ('-w')
                    call getarg(2 * (i - 1) + 2, arg)
                    read(arg,*) w
                    w = w / x0
                case ('-t')
                    call getarg(2 * (i - 1) + 2, arg)
                    read(arg,*) t_fin
                case ('-dt')
                    call getarg(2 * (i - 1) + 2, arg)
                    read(arg,*) dt
                case ('-v')
                    call getarg(2 * (i - 1) + 2, arg)
                    read(arg,*) vl
                    vl = vl / ph0
                case ('-r')
                    call getarg(2 * (i - 1) + 2, arg)
                    read(arg,*) r0
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
    end subroutine
end program
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
