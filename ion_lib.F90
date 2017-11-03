module ion_lib
    use props
    use ptcl_props
    implicit none
    
    contains
    
    subroutine ionEqn(g, i, j, ph, ni, b)
        type(grid), intent(in) :: g
        integer, intent(in)  :: i, j
        real(8), intent(in)  :: ph(:,:), ni(:,:)
        real(8), intent(out) :: b
        real(8) :: dfdx = 0, dfdy = 0, vx(2) = 0, vy(2) = 0, mu, Ti
        integer :: q = 1
        
        mu = 0.1
        Ti = 0.1
        
        if (g%nx > 1) then
            if (g%type_x(i-1,j-1) < 0) then
                vx(1) = 0
                call calc_flux(ni(i:i+1,j), ph(i:i+1,j), q, mu, Ti, g%dx(i), vx(2))
                
            else if (g%type_x(i-1,j-1) > 0) then
                call calc_flux(ni(i-1:i,j), ph(i-1:i,j), q, mu, Ti, g%dx(i-1), vx(1))
                vx(2) = 0
                
            else
                call calc_flux(ni(i-1:i,j), ph(i-1:i,j), q, mu, Ti, g%dx(i-1), vx(1))
                call calc_flux(ni(i:i+1,j), ph(i:i+1,j), q, mu, Ti, g%dx(i), vx(2))
            end if
            
            dfdx = (vx(2) - vx(1)) / g%dlx(i-1)
        end if
        
        if (g%ny > 1) then
            if (g%type_y(i-1,j-1) < 0) then
                vy(1) = 0
                call calc_flux(ni(i,j:j+1), ph(i,j:j+1), q, mu, Ti, g%dy(j), vy(2))
                
            else if (g%type_y(i-1,j-1) > 0) then
                call calc_flux(ni(i,j-1:j), ph(i,j-1:j), q, mu, Ti, g%dy(j-1), vy(1))
                vy(2) = 0
                
            else
                call calc_flux(ni(i,j-1:j), ph(i,j-1:j), q, mu, Ti, g%dy(j-1), vy(1))
                call calc_flux(ni(i,j:j+1), ph(i,j:j+1), q, mu, Ti, g%dy(j), vy(2))
            end if
            
            dfdy = (vy(2) - vy(1)) / g%dly(j-1)
        end if
        
        b = g%dt * (dfdx + dfdy)
    
    end subroutine
end module
    
