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
        real(8) :: dfdx = 0, dfdy = 0, vx(2) = 0, vy(2) = 0, D, mu, Ex, Ey
        integer :: q = 1, a
        
        mu = mui
        D = Di
        
        ! X-dir Fluxes
        if (g%nx > 1) then
            if (g%type_x(i-1,j-1) < 0) then                
                call calc_flux(ni(i:i+1,j), ph(i:i+1,j), q, mu, D, &
                               g%dx(i), vx(2))
            
                if (g%type_x(i-1,j-1) == -2) then
                    Ex = -(ph(i,j) - ph(i-1,j)) / g%dx(i-1)
                    if (Ex > 0) then
                        a = 1
                    else
                        a = 0
                    end if

                    vx(1) = (1 - a) * mui * Ex * ni(i,j) - 0.25 * vi * ni(i,j)
                else
                    vx(1) = 0
                end if
                
            else if (g%type_x(i-1,j-1) > 0) then
                call calc_flux(ni(i-1:i,j), ph(i-1:i,j), q, mu, D, &
                               g%dx(i-1), vx(1))
                
                if (g%type_x(i-1,j-1) == 2) then
                    Ex = -(ph(i+1,j) - ph(i,j)) / g%dx(i)
                    if (Ex < 0) then
                        a = 1
                    else
                        a = 0
                    end if

                    vx(2) = (1 - a) * mui * Ex * ni(i,j) + 0.25 * vi * ni(i,j)
                else
                    vx(2) = 0
                end if
                
            else
                call calc_flux(ni(i-1:i,j), ph(i-1:i,j), q, mu, D, &
                               g%dx(i-1), vx(1))
                call calc_flux(ni(i:i+1,j), ph(i:i+1,j), q, mu, D, &
                               g%dx(i), vx(2))
            end if
            
            dfdx = (vx(2) - vx(1)) / g%dlx(i-1)
        end if
        
        ! Y-Dir Fluxes
        if (g%ny > 1) then
            if (g%type_y(i-1,j-1) < 0) then                
                call calc_flux(ni(i,j:j+1), ph(i,j:j+1), q, mu, D, &
                               g%dy(j), vy(2))
            
                if (g%type_y(i-1,j-1) == -2) then
                    Ey = -(ph(i,j) - ph(i,j-1)) / g%dy(j-1)
                    if (Ey > 0) then
                        a = 1
                    else
                        a = 0
                    end if

                    vy(1) = (1 - a) * mui * Ey * ni(i,j) - 0.25 * vi * ni(i,j)
                else
                    vy(1) = 0
                end if
                
            else if (g%type_y(i-1,j-1) > 0) then
                call calc_flux(ni(i,j-1:j), ph(i,j-1:j), q, mu, D, &
                               g%dy(j-1), vy(1))
                
                if (g%type_y(i-1,j-1) == 2) then
                    Ey = -(ph(i,j+1) - ph(i,j)) / g%dy(j)
                    if (Ey < 0) then
                        a = 1
                    else
                        a = 0
                    end if

                    vx(2) = (1 - a) * mui * Ey * ni(i,j) + 0.25 * vi * ni(i,j)
                else
                    vx(2) = 0
                end if
                
            else
                call calc_flux(ni(i,j-1:j), ph(i,j-1:j), q, mu, D, &
                               g%dy(j-1), vy(1))
                call calc_flux(ni(i,j:j+1), ph(i,j:j+1), q, mu, D, &
                               g%dy(j), vy(2))
            end if
            
            dfdy = (vy(2) - vy(1)) / g%dly(j-1)
        end if
        
        b = g%dt * (dfdx + dfdy)
    end subroutine
    
    subroutine metaEqn(g, i, j, nm, b)
        type(grid), intent(in) :: g
        integer, intent(in)  :: i, j
        real(8), intent(in)  :: nm(:,:)
        real(8), intent(out) :: b
        real(8) :: dfdx = 0, dfdy = 0, vx(2) = 0, vy(2) = 0, D
        
        D = Dm
        
        ! X-dir Fluxes
        if (g%nx > 1) then
            if (g%type_x(i-1,j-1) < 0) then                
                vx(2) = -D * (nm(i+1,j) - nm(i,j)) / g%dx(i)
            
                if (g%type_x(i-1,j-1) == -2) then
                    vx(1) = -0.25 * vi * nm(i,j)
                else
                    vx(1) = 0
                end if
                
            else if (g%type_x(i-1,j-1) > 0) then
                vx(1) = -D * (nm(i,j) - nm(i-1,j)) / g%dx(i-1)
                
                if (g%type_x(i-1,j-1) == 2) then
                    vx(2) = 0.25 * vi * nm(i,j)
                else
                    vx(2) = 0
                end if
                
            else
                vx(1) = -D * (nm(i,j) - nm(i-1,j)) / g%dx(i-1)
                vx(2) = -D * (nm(i+1,j) - nm(i,j)) / g%dx(i)
            end if
            
            dfdx = (vx(2) - vx(1)) / g%dlx(i-1)
        end if
        
        ! Y-Dir Fluxes
        if (g%ny > 1) then
            if (g%type_y(i-1,j-1) < 0) then                
                vy(2) = -D * (nm(i,j+1) - nm(i,j)) / g%dy(j)
            
                if (g%type_y(i-1,j-1) == -2) then
                    vy(1) = -0.25 * vi * nm(i,j)
                else
                    vy(1) = 0
                end if
                
            else if (g%type_y(i-1,j-1) > 0) then
                vy(1) = -D * (nm(i,j) - nm(i,j-1)) / g%dy(j-1)
                
                if (g%type_y(i-1,j-1) == 2) then
                    vy(2) = 0.25 * vi * nm(i,j)
                else
                    vy(2) = 0
                end if
                
            else
                vy(1) = -D * (nm(i,j) - nm(i,j-1)) / g%dy(j-1)
                vy(2) = -D * (nm(i,j+1) - nm(i,j)) / g%dy(i)
            end if
            
            dfdy = (vy(2) - vy(1)) / g%dly(j-1)
        end if
        
        b = g%dt * (dfdx + dfdy)
    end subroutine
end module
    
