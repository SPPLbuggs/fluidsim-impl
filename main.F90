program main
    use props
    use petsc_lib
    use eqn_lib
    implicit none
    
    type(grid) :: g
    integer :: ts = 0, nx, ny, dof
    real(8) :: l, w, ew, vl, dt, t_fin, t_pr, t_sv, sim_start, time
    character(80):: path
    logical :: assem(3) = .True.
    
    ! Initialize PETSc and MPI
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
    l  = 1e-2 / x0
    w  = 1e-2 / x0
    ew = 5e-3 / x0
    dt = 1e-3
    t_fin = 10
    t_pr = t_fin/100.
    t_sv = t_fin/100.
    vl = 500 / ph0
    
    ! Read input arguments
    call read_in
    
    ! Initialize grid and arrays
    path = 'Output/'
    call g_init(g, nx, ny, px, py, dof, l, w, ew, trim(path))
    call eqn_init(g)
    
    g%t  = 0
    g%dt = dt
    
    do
        ts = ts + 1
        g%t = g%t + g%dt
        if (g%t >= t_fin) exit
        
        ! Update boundary conditions
        if (rx == 0) ph_pl(1,:,1) = vl!sin(2.0 * 3.14159 * g%t / 10.0)
        if (ry == 0) ph_pl(:,1,1) = vl!sin(2.0 * 3.14159 * g%t / 10.0)
        
        ! Solve laplace equation
        call petsc_step(g, A1, b1, x1, ph_pl, laplEval, assem(1))
        
        ! Solve electron equation
        ne_mi = ne_pl
        call petsc_step(g, A2, b2, x2, ne_pl, elecEval, assem(2))
        ne_pl = max(ne_pl, n_zero)
        
        ! Solve ion equation
        ni_mi = ni_pl
        call petsc_step(g, A3, b3, x3, ni_pl, ionEval, assem(3))
        ni_pl = max(ni_pl, n_zero)
        
        ! Print out some information
        if ((t_pr <= g%t) .and. (my_id == 0)) then
            call cpu_time(time)
            write(*,*)
            write(*,11) float(ts), g%t, g%dt, time - sim_start
            t_pr = t_pr + t_fin/100.
        end if
        
        ! Save data
        if (t_sv <= g%t) then
            call savedat(trim(path)//'f1.dat', ph_pl(:,:,1) * ph0)
            call savedat(trim(path)//'f2.dat', ne_pl(:,:,1) / x0**3)
            call savedat(trim(path)//'f3.dat', ni_pl(:,:,1) / x0**3)
            
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
    call petsc_destroy(A2, b2, x2)
    call petsc_destroy(A3, b3, x3)
    call PetscFinalize(ierr)

11 format('Timestep:', es9.2, '  Time:', es9.2, '  dT:', es9.2, '  Dur:', f6.2, ' sec')
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
