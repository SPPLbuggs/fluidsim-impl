module elec_lib
    use props
    use ptcl_props
    implicit none
    
    contains
    
    subroutine elecEqn(g, i, j, ph, ne, ne_mi, nte_mi, ni_mi, b)
        type(grid), intent(in) :: g
        integer, intent(in)  :: i, j
        real(8), intent(in)  :: ph(:,:), ne(:,:), ne_mi(:,:), nte_mi(:,:), ni_mi(:,:)
        real(8), intent(out) :: b
        real(8) :: dfdx = 0, dfdy = 0, vx(2) = 0, vy(2) = 0, Te(3), mu(2), D(2), flxi, Ex, Ey, ve
        integer :: a, q = -1
        
        if (g%nx > 1) then
            if (g%type_x(i-1,j-1) < 0) then
                Te(2) = get_Te(nte_mi(i,j),   ne_mi(i,j))
                Te(3) = get_Te(nte_mi(i+1,j), ne_mi(i+1,j))
                
                mu(2) = get_mue(0.5 * (Te(2) + Te(3)))
                D(2) = get_De(0.5 * (Te(2) + Te(3)))
                
                call calc_flux(ne(i:i+1,j), ph(i:i+1,j), q, mu(2), D(2), &
                               g%dx(i), vx(2))
            
                if (g%type_x(i-1,j-1) == -2) then
                    Ex = -(ph(i,j) - ph(i-1,j)) / g%dx(i-1)
                    if (Ex > 0) then
                        a = 1
                    else
                        a = 0
                    end if
                    
                    mu(1) = get_mue(Te(2))
                    ve = sqrt((16.0 * e * ph0 * Te(2)) / (3.0 * pi * me)) * t0 / x0
                    
                    ! Flux at i - 1/2
                    flxi = (1 - a) * mui * Ex * ni_mi(i,j) - 0.25 * vi * ni_mi(i,j)
                    
                    vx(1) = - a * mu(1) * Ex * ne(i,j) &
                            - 0.25 * ve * ne(i,j) &
                            - gam * flxi
                else
                    vx(1) = 0
                end if
                
            else if (g%type_x(i-1,j-1) > 0) then
                Te(1) = get_Te(nte_mi(i-1,j), ne_mi(i-1,j))
                Te(2) = get_Te(nte_mi(i,j),   ne_mi(i,j))
                
                mu(1) = get_mue(0.5 * (Te(1) + Te(2)))
                D(1) = get_De(0.5 * (Te(1) + Te(2)))
                
                call calc_flux(ne(i-1:i,j), ph(i-1:i,j), q, mu(1), D(1), &
                               g%dx(i-1), vx(1))
                
                if (g%type_x(i-1,j-1) == 2) then
                    Ex = -(ph(i+1,j) - ph(i,j)) / g%dx(i)
                    if (Ex < 0) then
                        a = 1
                    else
                        a = 0
                    end if
                    
                    mu(2) = get_mue(Te(2))
                    ve = sqrt((16.0 * e * ph0 * Te(2)) / (3.0 * pi * me)) * t0 / x0
                    
                    ! Flux at i - 1/2
                    flxi = (1 - a) * mui * Ex * ni_mi(i,j) + 0.25 * vi * ni_mi(i,j)
                    
                    vx(2) = - a * mu(2) * Ex * ne(i,j) &
                            + 0.25 * ve * ne(i,j) &
                            - gam * flxi
                else
                    vx(2) = 0
                end if
                
            else
                Te(1) = get_Te(nte_mi(i-1,j), ne_mi(i-1,j))
                Te(2) = get_Te(nte_mi(i,j),   ne_mi(i,j))
                Te(3) = get_Te(nte_mi(i+1,j), ne_mi(i+1,j))
                
                mu(1) = get_mue(0.5 * (Te(1) + Te(2)))
                mu(2) = get_mue(0.5 * (Te(2) + Te(3)))
                
                D(1) = get_De(0.5 * (Te(1) + Te(2)))
                D(2) = get_De(0.5 * (Te(2) + Te(3)))
                
                call calc_flux(ne(i-1:i,j), ph(i-1:i,j), q, mu(1), D(1), &
                               g%dx(i-1), vx(1))
                call calc_flux(ne(i:i+1,j), ph(i:i+1,j), q, mu(2), D(2), &
                               g%dx(i), vx(2))
            end if
            dfdx = (vx(2) - vx(1)) / g%dlx(i-1)
        end if
        
        if (g%ny > 1) then
            Te(1) = get_Te(nte_mi(i,j-1), ne_mi(i,j-1))
            Te(2) = get_Te(nte_mi(i,j),   ne_mi(i,j))
            Te(3) = get_Te(nte_mi(i,j+1), ne_mi(i,j+1))
            
            mu(1) = get_mue(0.5 * (Te(1) + Te(2)))
            mu(2) = get_mue(0.5 * (Te(2) + Te(3)))
            
            D(1) = get_De(0.5 * (Te(1) + Te(2)))
            D(2) = get_De(0.5 * (Te(2) + Te(3)))
            
            call calc_flux(ne(i,j-1:j), ph(i,j-1:j), q, mu(1), D(1), g%dy(j-1), vy(1))
            call calc_flux(ne(i,j:j+1), ph(i,j:j+1), q, mu(1), D(1), g%dy(j), vy(2))
            dfdy = (vy(2) - vy(1)) / g%dly(j-1)
        end if
        
        b = g%dt * (dfdx + dfdy)
    end subroutine
    
    subroutine elecEnrgEqn(g, i, j, ph, nte, ne_mi, nte_mi, ni_mi, b)
        type(grid), intent(in) :: g
        integer, intent(in)  :: i, j
        real(8), intent(in)  :: ph(:,:), nte(:,:), ne_mi(:,:), nte_mi(:,:), ni_mi(:,:)
        real(8), intent(out) :: b
        real(8) :: dfdx = 0, dfdy = 0, vx(2) = 0, vy(2) = 0, Te(3), mu(2), D(2), flxi, Ex, Ey, ve
        integer :: a, q = -1
        
        if (g%nx > 1) then
            if (g%type_x(i-1,j-1) < 0) then
                Te(2) = get_Te(nte_mi(i,j),   ne_mi(i,j))
                Te(3) = get_Te(nte_mi(i+1,j), ne_mi(i+1,j))
                
                mu(2) = get_mut(0.5 * (Te(2) + Te(3)))
                D(2) = get_Dt(0.5 * (Te(2) + Te(3)))
                
                call calc_flux(nte(i:i+1,j), ph(i:i+1,j), q, mu(2), D(2), &
                               g%dx(i), vx(2))
            
                if (g%type_x(i-1,j-1) == -2) then
                    Ex = -(ph(i,j) - ph(i-1,j)) / g%dx(i-1)
                    if (Ex > 0) then
                        a = 1
                    else
                        a = 0
                    end if
                    
                    mu(1) = get_mut(Te(2))
                    ve = sqrt((16.0 * e * ph0 * Te(2)) / (3.0 * pi * me)) * t0 / x0
                    
                    ! Flux at i - 1/2
                    flxi = (1 - a) * mui * Ex * ni_mi(i,j) - 0.25 * vi * ni_mi(i,j)
                    
                    vx(1) = - a * mu(1) * Ex * nte(i,j) &
                            - 1.0/3.0 * ve * nte(i,j) &
                            - gam * Te(2) * flxi
                else
                    vx(1) = 0
                end if
                
            else if (g%type_x(i-1,j-1) > 0) then
                Te(1) = get_Te(nte_mi(i-1,j), ne_mi(i-1,j))
                Te(2) = get_Te(nte_mi(i,j),   ne_mi(i,j))
                
                mu(1) = get_mut(0.5 * (Te(1) + Te(2)))
                D(1) = get_Dt(0.5 * (Te(1) + Te(2)))
                
                call calc_flux(nte(i-1:i,j), ph(i-1:i,j), q, mu(1), D(1), &
                               g%dx(i-1), vx(1))
                
                if (g%type_x(i-1,j-1) == 2) then
                    Ex = -(ph(i+1,j) - ph(i,j)) / g%dx(i)
                    if (Ex < 0) then
                        a = 1
                    else
                        a = 0
                    end if
                    
                    mu(2) = get_mut(Te(2))
                    ve = sqrt((16.0 * e * ph0 * Te(2)) / (3.0 * pi * me)) * t0 / x0
                    
                    ! Flux at i - 1/2
                    flxi = (1 - a) * mui * Ex * ni_mi(i,j) + 0.25 * vi * ni_mi(i,j)
                    
                    vx(2) = - a * mu(2) * Ex * nte(i,j) &
                            + 1.0/3.0 * ve * nte(i,j) &
                            - gam * Te(2) * flxi
                else
                    vx(2) = 0
                end if
                
            else
                Te(1) = get_Te(nte_mi(i-1,j), ne_mi(i-1,j))
                Te(2) = get_Te(nte_mi(i,j),   ne_mi(i,j))
                Te(3) = get_Te(nte_mi(i+1,j), ne_mi(i+1,j))
                
                mu(1) = get_mut(0.5 * (Te(1) + Te(2)))
                mu(2) = get_mut(0.5 * (Te(2) + Te(3)))
                
                D(1) = get_Dt(0.5 * (Te(1) + Te(2)))
                D(2) = get_Dt(0.5 * (Te(2) + Te(3)))
                
                call calc_flux(nte(i-1:i,j), ph(i-1:i,j), q, mu(1), D(1), &
                               g%dx(i-1), vx(1))
                call calc_flux(nte(i:i+1,j), ph(i:i+1,j), q, mu(2), D(2), &
                               g%dx(i), vx(2))
            end if
            dfdx = (vx(2) - vx(1)) / g%dlx(i-1)
        end if
        
        if (g%ny > 1) then
            Te(1) = get_Te(nte_mi(i,j-1), ne_mi(i,j-1))
            Te(2) = get_Te(nte_mi(i,j),   ne_mi(i,j))
            Te(3) = get_Te(nte_mi(i,j+1), ne_mi(i,j+1))
            
            mu(1) = get_mut(0.5 * (Te(1) + Te(2)))
            mu(2) = get_mut(0.5 * (Te(2) + Te(3)))
            
            D(1) = get_Dt(0.5 * (Te(1) + Te(2)))
            D(2) = get_Dt(0.5 * (Te(2) + Te(3)))
            
            call calc_flux(nte(i,j-1:j), ph(i,j-1:j), q, mu(1), D(1), g%dy(j-1), vy(1))
            call calc_flux(nte(i,j:j+1), ph(i,j:j+1), q, mu(1), D(1), g%dy(j), vy(2))
            dfdy = (vy(2) - vy(1)) / g%dly(j-1)
        end if
        
        b = g%dt * (dfdx + dfdy)
    end subroutine
end module
    
