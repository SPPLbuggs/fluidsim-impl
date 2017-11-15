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

  contains

! *** External Circuit Timestepping ***
  subroutine circ_step(g)
    type(grid), intent(in) :: g
    integer :: i, j
    real(8) :: num, den, mue, Te, ve, a, temp, dr

    i = 2
    num = 0
    den = 0
    Id = 0
    if (rx == 0) then
      do j = 2, g%by+1
        Te = get_Te(ne_pl(i,j,2), ne_pl(i,j,1))
        mue = get_mue(Te)
        ve = sqrt((16.0 * e * ph0 * Te) / (3.0 * pi * me)) * t0 / x0

        a = 0
        if ((Vd_pl - ph_pl(i,j,1)) < 0) a = 1

        temp = (a * (1+gam) * mui * ni_pl(i,j,1) &
               + (1 - a) * mue * ne_pl(i,j,1)) / g%dx(i-1)

        if (g%ny > 1) then
          if (cyl) then
            dr = g%dy(j-1) * g%r(j-1) * 2 * pi
          else
            dr = g%dy(j-1)
          end if
        else
          dr = g%w**2 * pi
        end if

        den = den + g%dt / Cap * temp * dr

        num = num + g%dt / Cap * (ph_pl(i,j,1) * temp &
              + (1 + gam) / 4.0 * vi * ni_pl(i,j,1) &
              - 1.0 / 4.0 * ve * ne_pl(i,j,1)) * dr

        Id = Id + ((Vd_pl - ph_pl(i,j,1)) * temp &
             - (1 + gam) / 4.0 * vi * ni_pl(i,j,1) &
             + 1.0 / 4.0 * ve * ne_pl(i,j,1)) * dr
      end do
    end if

    call MPI_Allreduce(MPI_In_Place, num, 1, etype, MPI_Sum, comm, ierr)
    call MPI_Allreduce(MPI_In_Place, den, 1, etype, MPI_Sum, comm, ierr)

    num = num + Vd_pl + g%dt / (Cap * Res) * Vsrc
    den = den + 1 + g%dt / (Cap * Res)

    Vd_pl = num / den

    do j = 2, g%by+1
      if (g%type_x(i-1,j-1) == -2) ph_pl(i-1,j,1) = Vd_pl
    end do
  end subroutine

  ! *** External Circuit Initialization ***
  subroutine circ_init(Vset, R0)
    real(8), intent(in) :: Vset, R0

    Vmax  = Vset
    Vsrc  = Vset
    Vd_pl = Vset
    Cap = 1e-12 * ph0 / e
    Res = R0 * e / (ph0 * t0)
  end subroutine
end module
