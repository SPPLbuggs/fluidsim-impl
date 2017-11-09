! *** External Circuit Module ***
    module circ_lib
    use props
    use lapl_lib
    use elec_lib
    use ion_lib
    implicit none
    
    real(8) :: Vmax, Vd_pl, Id, Vsrc, Cap, Res
    
    ! variables
    public  :: Vd_pl, Id
    private :: Vmax, Vsrc, Cap, Res
    
    ! subroutines
    public  :: circ_step, circ_init
    !private :: 
    
    contains

! *** External Circuit Timestepping ***
    subroutine circ_step(g)
        type(grid), intent(in) :: g
        integer :: i, j, stage
        real(8) :: b, Rtemp, Vd_mi, Vd_or, &
                   flxi_x(2), flxe_x(2), flxi_y(2), flxe_y(2)
        
        Vd_mi = Vd_pl
        Vd_or = Vd_pl
        
        if (g%t < 1d0 / 16d0) then
            Vsrc = Vmax * sin(8d0 * pi * g%t)**3
        else
            Vsrc = Vmax
        end if
        
        do stage = 1, 4
            Id = 0
            
            if (rx == 0) then
                ph_pl(1,:,1) = Vd_mi
                
                i = 2
                do j = 2, g%by + 1
                
                    call ionFlx(g, i, j, ph_pl(:,:,1), ni_pl(:,:,1), flxi_x, flxi_y)
                    call elecFlx(g, i, j, ph_pl(:,:,1), ne_pl(:,:,1), ni_pl(:,:,1), &
                                 nte_pl(:,:,1), flxe_x, flxe_y)
                    
                    if (g%ny > 1) then
                        Id = Id + (flxi_x(1) - flxe_x(1)) * g%dy(j-1)
                    else
                        Id = Id + (flxi_x(1) - flxe_x(1)) * g%w**2 * pi
                    end if
                    
                end do
            end if
            
            call MPI_Allreduce(MPI_In_Place, Id, 1, etype, MPI_Sum, comm, ierr)
            
            if (g%t < 7d-2) then
                Rtemp = min(1e4 * e / (ph0 * t0), Res)
            else
                Rtemp = Res
            end if
            
            b = -1.0 / Cap * (Id - (Vsrc - Vd_mi)/Rtemp)
            
            if (stage == 1) then
                Vd_pl = Vd_or + b * g%dt / 6.0
            else if (stage == 2) then
                Vd_pl = Vd_pl + b * g%dt / 3.0
                Vd_mi = Vd_or + b * g%dt / 2.0
            else if (stage == 3) then
                Vd_pl = Vd_pl + b * g%dt / 3.0
                Vd_mi = Vd_or + b * g%dt
            else
                Vd_pl = Vd_pl + b * g%dt / 6.0
            end if
        end do
    end subroutine

! *** External Circuit Initialization ***
    subroutine circ_init(Vset, R0)
    real(8), intent(in) :: Vset, R0
    
    Vmax  = Vset
    Vsrc  = 0
    Vd_pl = 0
    Cap = 1e-12 * ph0 / e
    Res = R0 * e / (ph0 * t0)
    
    end subroutine
    end module
