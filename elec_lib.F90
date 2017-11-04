module elec_lib
    use props
    use ptcl_props
    implicit none
    
    contains
    
    subroutine elecEqn(g, i, j, ph, ne, b)
        type(grid), intent(in) :: g
        integer, intent(in)  :: i, j
        real(8), intent(in)  :: ph(:,:), ne(:,:)
        real(8), intent(out) :: b
        real(8) :: dfdx = 0, dfdy = 0, vx(2) = 0, vy(2) = 0, mue, De
        integer :: q = -1
        
        mue = get_mue(1.5 / ph0)
        De = get_De(1.5 / ph0)
        
        if (g%nx > 1) then
            call calc_flux(ne(i-1:i,j), ph(i-1:i,j), q, mue, De, g%dx(i-1), vx(1))
            call calc_flux(ne(i:i+1,j), ph(i:i+1,j), q, mue, De, g%dx(i), vx(2))        
            dfdx = (vx(2) - vx(1)) / g%dlx(i-1)
        end if
        
        if (g%ny > 1) then
            call calc_flux(ne(i,j-1:j), ph(i,j-1:j), q, mue, De, g%dy(j-1), vy(1))
            call calc_flux(ne(i,j:j+1), ph(i,j:j+1), q, mue, De, g%dy(j), vy(2))
            dfdy = (vy(2) - vy(1)) / g%dly(j-1)
        end if
        
        b = g%dt * (dfdx + dfdy)
    
    end subroutine
end module
    
