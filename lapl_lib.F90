module lapl_lib
    use props
    use elec_lib
    use ion_lib
    implicit none
    
    contains

    subroutine laplEqn(g, i, j, ph, ne, ni, b)
        type(grid), intent(in) :: g
        integer, intent(in)  :: i, j
        real(8), intent(in)  :: ph(:,:), ne(:,:), ni(:,:)
        real(8), intent(out) :: b
        real(8) :: dfdx = 0, dfdy = 0, dflxe = 0, dflxi = 0
        
        if (g%nx > 1) then
            ! Left symmetry boundary
            if (g%type_x(i-1,j-1) == -1) then
                dfdx = (ph(i+1,j) - ph(i,j)) / g%dx(i) / g%dlx(i-1)
            
            ! Right symmetry boundary
            else if (g%type_x(i-1,j-1) == 1) then
                dfdx = (ph(i-1,j) - ph(i,j)) / g%dx(i-1) / g%dlx(i-1)
            
            ! Domain center
            else
                dfdx = ((ph(i+1,j) - ph(i,j)) / g%dx(i) &
                       -(ph(i,j) - ph(i-1,j)) / g%dx(i-1)) &
                       / g%dlx(i-1)
           end if
        end if
        
        if (g%ny > 1) then
            ! Left symmetry boundary
            if (g%type_y(i-1,j-1) == -1) then
                dfdy = (ph(i,j+1) - ph(i,j)) / g%dy(j) / g%dly(j-1)
            
            ! Right symmetry boundary
            else if (g%type_y(i-1,j-1) == 1) then
                dfdy = (ph(i,j-1) - ph(i,j)) / g%dy(j-1) / g%dly(j-1)
            
            ! Domain center
            else
                dfdy = ((ph(i,j+1) - ph(i,j)) / g%dy(j) &
                       -(ph(i,j) - ph(i,j-1)) / g%dy(j-1)) &
                       / g%dly(j-1)
           end if
        end if
        
        call elecEqn(g, i, j, ph, ne, dflxe)
        call ionEqn(g, i, j, ph, ni, dflxi)
        
        b = dfdx + dfdy + (ni(i,j) - ne(i,j) + dflxe - dflxi) * 100.
    end subroutine
end module