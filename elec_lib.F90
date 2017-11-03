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
        real(8) :: dfdx = 0, dfdy = 0, vx(2) = 0, vy(2) = 0, mu, Te
        integer :: q = -1
        
        mu = 0.1
        Te = 0.1
        
        if (g%nx > 1) then
            if (g%type_x(i-1,j-1) < 0) then
                vx(1) = 0
                call calc_flux(ne(i:i+1,j), ph(i:i+1,j), q, mu, Te, g%dx(i), vx(2))
                
            else if (g%type_x(i-1,j-1) > 0) then
                call calc_flux(ne(i-1:i,j), ph(i-1:i,j), q, mu, Te, g%dx(i-1), vx(1))
                vx(2) = 0
                
            else
                call calc_flux(ne(i-1:i,j), ph(i-1:i,j), q, mu, Te, g%dx(i-1), vx(1))
                call calc_flux(ne(i:i+1,j), ph(i:i+1,j), q, mu, Te, g%dx(i), vx(2))
            end if
            
            dfdx = (vx(2) - vx(1)) / g%dlx(i-1)
        end if
        
        if (g%ny > 1) then
            if (g%type_y(i-1,j-1) < 0) then
                vy(1) = 0
                call calc_flux(ne(i,j:j+1), ph(i,j:j+1), q, mu, Te, g%dy(j), vy(2))
                
            else if (g%type_y(i-1,j-1) > 0) then
                call calc_flux(ne(i,j-1:j), ph(i,j-1:j), q, mu, Te, g%dy(j-1), vy(1))
                vy(2) = 0
                
            else
                call calc_flux(ne(i,j-1:j), ph(i,j-1:j), q, mu, Te, g%dy(j-1), vy(1))
                call calc_flux(ne(i,j:j+1), ph(i,j:j+1), q, mu, Te, g%dy(j), vy(2))
            end if
            
            dfdy = (vy(2) - vy(1)) / g%dly(j-1)
        end if
        
        b = g%dt * (dfdx + dfdy)
    
    end subroutine
end module
    
