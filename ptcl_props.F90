module ptcl_props
    use props
    implicit none
    
contains
    
    subroutine calc_flux2(n, ph, q, mu, Te, dx, flx)
        integer, intent(in)  :: q
        real(8), intent(in)  :: n(2), ph(2), mu, Te, dx
        real(8), intent(out) :: flx
        
        flx = 0.5 * mu * (n(2) + n(1)) * (-q * (ph(2) - ph(1)) &
                    - Te * log(n(2) / n(1))) / dx
    end subroutine
    
    subroutine calc_flux(n, ph, q, mu, Te, dx, flx)
        integer, intent(in)  :: q
        real(8), intent(in)  :: n(2), ph(2), mu, Te, dx
        real(8), intent(out) :: flx
        real(8) :: v, tol, D, arg
        
        tol = 1e-12
        v = -q * mu * (ph(2) - ph(1)) / dx
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
