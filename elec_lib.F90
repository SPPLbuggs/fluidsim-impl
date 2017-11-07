module elec_lib
    use props
    use ptcl_props
    implicit none
    
    contains
    
    subroutine elecEqn(g, i, j, ph, ne, ne_mi, nte_mi, b)
        type(grid), intent(in) :: g
        integer, intent(in)  :: i, j
        real(8), intent(in)  :: ph(:,:), ne(:,:), ne_mi(:,:), nte_mi(:,:)
        real(8), intent(out) :: b
        real(8) :: dfdx = 0, dfdy = 0, vx(2) = 0, vy(2) = 0, Te(3), mu(2), D(2)
        integer :: q = -1
        
        if (g%nx > 1) then
            Te(1) = get_Te(nte_mi(i-1,j), ne_mi(i-1,j))
            Te(2) = get_Te(nte_mi(i,j),   ne_mi(i,j))
            Te(3) = get_Te(nte_mi(i+1,j), ne_mi(i+1,j))
            
            mu(1) = get_mue(0.5 * (Te(1) + Te(2)))
            mu(2) = get_mue(0.5 * (Te(2) + Te(3)))
            
            D(1) = get_De(0.5 * (Te(1) + Te(2)))
            D(2) = get_De(0.5 * (Te(2) + Te(3)))
            
            call calc_flux(ne(i-1:i,j), ph(i-1:i,j), q, mu(1), D(1), g%dx(i-1), vx(1))
            call calc_flux(ne(i:i+1,j), ph(i:i+1,j), q, mu(1), D(1), g%dx(i), vx(2))        
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
    
    subroutine elecEnrgEqn(g, i, j, ph, nte, ne_mi, nte_mi, b)
        type(grid), intent(in) :: g
        integer, intent(in)  :: i, j
        real(8), intent(in)  :: ph(:,:), nte(:,:), ne_mi(:,:), nte_mi(:,:)
        real(8), intent(out) :: b
        real(8) :: dfdx = 0, dfdy = 0, vx(2) = 0, vy(2) = 0, Te(3), mu(2), D(2)
        integer :: q = -1
        
        if (g%nx > 1) then
            Te(1) = get_Te(nte_mi(i-1,j), ne_mi(i-1,j))
            Te(2) = get_Te(nte_mi(i,j),   ne_mi(i,j))
            Te(3) = get_Te(nte_mi(i+1,j), ne_mi(i+1,j))
            
            mu(1) = get_mut(0.5 * (Te(1) + Te(2)))
            mu(2) = get_mut(0.5 * (Te(2) + Te(3)))
            
            D(1) = get_Dt(0.5 * (Te(1) + Te(2)))
            D(2) = get_Dt(0.5 * (Te(2) + Te(3)))
            
            call calc_flux(nte(i-1:i,j), ph(i-1:i,j), q, mu(1), D(1), g%dx(i-1), vx(1))
            call calc_flux(nte(i:i+1,j), ph(i:i+1,j), q, mu(1), D(1), g%dx(i), vx(2))        
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
    
