    module eqn_lib
    use props
    implicit none
    
    real(8), allocatable :: f(:,:,:), f_mi(:,:,:)
    
    contains
    
    subroutine feval1(g, i, j, n, m, dof, f1, b_temp)
    type(grid), intent(in) :: g
    integer, intent(in) :: i, j, n, m, dof
    real(8), intent(in) :: f1(n,m,dof)
    real(8), intent(out) :: b_temp(dof)
    
    call feqn1(g, i, j, f1(:,:,1), f(:,:,2), b_temp(1))
    
    end subroutine
    
    subroutine feval2(g, i, j, n, m, dof, f2, b_temp)
    type(grid), intent(in) :: g
    integer, intent(in) :: i, j, n, m, dof
    real(8), intent(in) :: f2(n,m,dof)
    real(8), intent(out) :: b_temp(dof)
    real(8) :: dflx_pl, dflx_mi
    
    call feqn2(g, i, j, f(:,:,1), f2(:,:,1), dflx_pl)
    call feqn2(g, i, j, f(:,:,1), f_mi(:,:,2), dflx_mi)
    
    b_temp = f2(i,j,1) - f_mi(i,j,2) + 0.5 * dflx_pl + 0.5 * dflx_mi
    
    end subroutine
    
    subroutine feqn1(g, i, j, f1, f2, b)
    type(grid), intent(in) :: g
    integer, intent(in)  :: i, j
    real(8), intent(in)  :: f1(:,:), f2(:,:)
    real(8), intent(out) :: b
    real(8) :: dfdx = 0, dfdy = 0, dflx = 0
    
    if (g%nx > 1) then
        dfdx = ((f1(i+1,j) - f1(i,j)) / g%dx(i) &
               -(f1(i,j) - f1(i-1,j)) / g%dx(i-1)) &
               / g%dlx(i-1)
    end if
    
    if (g%ny > 1) then
        dfdy = ((f1(i,j+1) - f1(i,j)) / g%dy(j) &
               -(f1(i,j) - f1(i,j-1)) / g%dy(j-1)) &
               / g%dly(j-1)
    end if
    
    call feqn2(g, i, j, f1, f2, dflx)
    
    b = dfdx + dfdy + (1.0 - f2(i,j) + dflx) * 100.
    
    end subroutine
    
    subroutine feqn2(g, i, j, f1, f2, b)
    type(grid), intent(in) :: g
    integer, intent(in)  :: i, j
    real(8), intent(in)  :: f1(:,:), f2(:,:)
    real(8), intent(out) :: b
    real(8) :: dfdx = 0, dfdy = 0, vx(2) = 0, vy(2) = 0, mu, Te
    
    mu = 0.1
    Te = 0.1
    
    if (g%nx > 1) then
        if (g%type_x(i-1,j-1) == -1) then
            vx(1) = 0
            call calc_flux(f2(i:i+1,j), f1(i:i+1,j), mu, Te, g%dx(i), vx(2))
            
        else if (g%type_x(i-1,j-1) == 1) then
            call calc_flux(f2(i-1:i,j), f1(i-1:i,j), mu, Te, g%dx(i-1), vx(1))
            vx(2) = 0
            
        else
            call calc_flux(f2(i-1:i,j), f1(i-1:i,j), mu, Te, g%dx(i-1), vx(1))
            call calc_flux(f2(i:i+1,j), f1(i:i+1,j), mu, Te, g%dx(i), vx(2))
        end if
        
        dfdx = (vx(2) - vx(1)) / g%dlx(i-1)
    end if
    
    if (g%ny > 1) then
        if (g%type_y(i-1,j-1) == -1) then
            vy(1) = 0
            call calc_flux(f2(i,j:j+1), f1(i,j:j+1), mu, Te, g%dy(j), vy(2))
            
        else if (g%type_y(i-1,j-1) == 1) then
            call calc_flux(f2(i,j-1:j), f1(i,j-1:j), mu, Te, g%dy(j-1), vy(1))
            vy(2) = 0
            
        else
            call calc_flux(f2(i,j-1:j), f1(i,j-1:j), mu, Te, g%dy(j-1), vy(1))
            call calc_flux(f2(i,j:j+1), f1(i,j:j+1), mu, Te, g%dy(j), vy(2))
        end if
        
        dfdy = (vy(2) - vy(1)) / g%dly(j-1)
    end if
    
    b = g%dt * (dfdx + dfdy)
    
    end subroutine
    
    subroutine calc_flux2(n, p, mu, Te, dx, flx)
    real(8), intent(in)  :: n(2), p(2), mu, Te, dx
    real(8), intent(out) :: flx
    
    flx = 0.5 * mu * (n(2) + n(1)) * (p(2) - p(1) &
                    - Te * log(n(2) / n(1))) / dx
    end subroutine
    
    subroutine calc_flux(n, p, mu, Te, dx, flx)
    real(8), intent(in)  :: n(2), p(2), mu, Te, dx
    real(8), intent(out) :: flx
    real(8) :: v, tol, D, arg
    
    tol = 1e-12
    v = mu * (p(2) - p(1)) / dx
    D = mu * Te
    arg = v * dx / D
    
    if (abs(v) < tol) then
        flx = D * (n(1) - n(2)) / dx
    else if (arg > 0) then
        flx = v * (n(1) - n(2) * exp(-arg)) / (1.0 - exp(-arg))
    else
        flx = v * (n(2) - n(1) * exp( arg)) / (1.0 - exp( arg))
    end if
    end subroutine
    end module
    
    
    
    
    
