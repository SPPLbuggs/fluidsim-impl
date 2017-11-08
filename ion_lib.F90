module ion_lib
    use props
    use ptcl_props
    implicit none
    
    contains
    
    ! *** Ion Source Term ***
    subroutine ionSrc(ne, ni, nte, nm, src)
        real(8), intent(in)  :: ni, ne, nte, nm
        real(8), intent(out) :: src
        real(8) :: Te, k_ir, k_si
        
        ! rates and coefficients
        Te = get_Te(nte, ne)
        k_ir = get_k_ir(Te)
        k_si = get_k_si(Te)

        ! evaluate source terms
        src =   k_ir * ninf    * ne &
              - beta * ni * ne &
              + k_si * nm * ne &
              + k_mp * nm**2
    end subroutine

    ! *** Metastable Source Term ***
    subroutine metaSrc(ne, nte, nm, src)
        real(8), intent(in)  :: ne, nte, nm
        real(8), intent(out) :: src
        real(8) :: Te, k_ex, k_sc, k_si, nu
        
        ! rates and coefficients
        Te = get_Te(nte, ne)
        k_sc = get_k_sc(Te)
        k_si = get_k_si(Te)
        k_ex = get_k_ex(Te)
        nu   = get_nu(Te)

        ! evaluate source term
        src =   k_ex * ninf * ne &
              - k_si * nm   * ne &
              - k_sc * nm   * ne &
              - k_r  * nm   * ne &
              - 2d0  * k_mp * nm**2 &
              - k_2q * ninf * nm  &
              - k_3q * ninf**2 * nm
    end subroutine
    
    ! *** Divergence of Ion Flux ***
    subroutine ionDFlx(g, i, j, ph, ni, dflx)
        type(grid), intent(in) :: g
        integer, intent(in) :: i, j
        real(8), intent(in) :: ph(:,:), ni(:,:)
        real(8), intent(out) :: dflx
        real(8) :: flx_x(2), flx_y(2)
        
        call ionFlx(g, i, j, ph, ni, flx_x, flx_y)
        
        dflx = (flx_x(2) - flx_x(1)) / g%dlx(i-1)
    end subroutine
    
    ! *** Divergence of Metastable Flux ***
    subroutine metaDFlx(g, i, j, nm, dflx)
        type(grid), intent(in) :: g
        integer, intent(in) :: i, j
        real(8), intent(in) :: nm(:,:)
        real(8), intent(out) :: dflx
        real(8) :: flx_x(2), flx_y(2)
        
        call metaFlx(g, i, j, nm, flx_x, flx_y)
        
        dflx = (flx_x(2) - flx_x(1)) / g%dlx(i-1)
    end subroutine
    
    ! *** Ion Flux ***
    subroutine ionFlx(g, i, j, ph, ni, flx_x, flx_y)
        type(grid), intent(in) :: g
        integer, intent(in)    :: i, j
        real(8), intent(in)    :: ph(:,:), ni(:,:)
        real(8), intent(out)   :: flx_x(2), flx_y(2)
        real(8) :: a, Ex(2)
        
        flx_x = 0
        flx_y = 0
        
        ! X-dir fields:
        Ex(1) = -(ph(i,j) - ph(i-1,j)) / g%dx(i-1)
        Ex(2) = -(ph(i+1,j) - ph(i,j)) / g%dx(i)
        
        ! X-dir Fluxes:
        ! - center -
        if (g%type_x(i-1,j-1) == 0) then
            ! Flux at i - 1/2
            call getFlx(flx_x(1), Ex(1), g%dx(i-1), 1, mui, Di, &
                          ni(i-1,j), ni(i,j))

            ! Flux at i + 1/2
            call getFlx(flx_x(2), Ex(2), g%dx(i), 1, mui, Di, &
                          ni(i,j), ni(i+1,j))
        
        ! - left -
        else if (g%type_x(i-1,j-1) < 0) then
            ! Flux at i + 1/2
            call getFlx(flx_x(2), Ex(2), g%dx(i), 1, mui, Di, &
                          ni(i,j), ni(i+1,j))
            
            ! - electrode -
            if (g%type_x(i-1,j-1) == -2) then
                if (Ex(1) < 0) then
                    a = 1
                else
                    a = 0
                end if
                
                ! Flux at i - 1/2
                flx_x(1) = a * mui * Ex(1) * ni(i,j) - 0.25 * vi * ni(i,j)

            ! - vacuum -
            else if (g%type_x(i-1,j-1) == -1) then
                ! Flux at i - 1/2
                flx_x(1) = 0
            end if
        
        ! - right -
        else if (g%type_x(i-1,j-1) > 0) then
            ! Flux at i - 1/2
            call getFlx(flx_x(1), Ex(1), g%dx(i-1), 1, mui, Di, &
                          ni(i-1,j), ni(i,j))
            
            ! - electrode -
            if (g%type_x(i-1,j-1) == 2) then
                if (Ex(2) > 0) then
                    a = 1
                else
                    a = 0
                end if
                
                ! Flux at i + 1/2
                flx_x(2) = a * mui * Ex(2) * ni(i,j) + 0.25 * vi * ni(i,j)
                
            ! - vacuum -
            else if (g%type_x(i-1,j-1) == 1) then
                ! Flux at i + 1/2
                flx_x(2) = 0
            end if
        end if
    end subroutine
    
    ! *** Metastable Flux ***
    subroutine metaFlx(g, i, j, nm, flx_x, flx_y)
        type(grid), intent(in) :: g
        integer, intent(in) :: i, j
        real(8), intent(in) :: nm(:,:)
        real(8), intent(out) :: flx_x(2), flx_y(2)
        
        flx_x = 0
        flx_y = 0
        
        ! X-dir Fluxes:
        ! - center -
        if (g%type_x(i-1,j-1) == 0) then
            ! Flux at i - 1/2
            flx_x(1) = -Dm * (nm(i,j) - nm(i-1,j)) / g%dx(i-1)
            
            ! Flux at i + 1/2
            flx_x(2) = -Dm * (nm(i+1,j) - nm(i,j)) / g%dx(i)
        
        ! - left -
        else if (g%type_x(i-1,j-1) < 0) then
            ! Flux at i + 1/2
            flx_x(2) = -Dm * (nm(i+1,j) - nm(i,j)) / g%dx(i)
            
            ! - electrode -
            if (g%type_x(i-1,j-1) == -2) then
                
                ! Flux at i - 1/2
                flx_x(1) = - 0.25 * vi * nm(i,j)

            ! - vacuum -
            else if (g%type_x(i-1,j-1) == -1) then
                ! Flux at i - 1/2
                flx_x(1) = 0
            end if
        
        ! - right -
        else if (g%type_x(i-1,j-1) > 0) then
            ! Flux at i - 1/2
            flx_x(1) = -Dm * (nm(i,j) - nm(i-1,j)) / g%dx(i-1)
            
            ! - electrode -
            if (g%type_x(i-1,j-1) == 2) then
                
                ! Flux at i + 1/2
                flx_x(2) = 0.25 * vi * nm(i,j)
                
            ! - vacuum -
            else if (g%type_x(i-1,j-1) == 1) then
                ! Flux at i + 1/2
                flx_x(2) = 0
            end if
        end if
        
        ! Y-dir Fluxes
        if (g%ny > 1) then
            ! - center -
            if (g%type_y(i-1,j-1) == 0) then
                ! Flux at j - 1/2
                flx_y(1) = -Dm * (nm(i,j) - nm(i,j-1)) / g%dy(j-1)
                
                ! Flux at j + 1/2
                flx_y(2) = -Dm * (nm(i,j+1) - nm(i,j)) / g%dy(j)
            
            ! - left -
            else if (g%type_y(i-1,j-1) < 0) then
                ! Flux at j + 1/2
                flx_y(2) = -Dm * (nm(i,j+1) - nm(i,j)) / g%dy(j)

                ! Flux at j - 1/2
                flx_y(1) = 0
            
            ! - right -
            else if (g%type_y(i-1,j-1) > 0) then
                ! Flux at j - 1/2
                flx_y(1) = -Dm * (nm(i,j) - nm(i,j-1)) / g%dy(j-1)

                ! Flux at j + 1/2
                flx_y(2) = 0
            end if
        end if
    end subroutine
end module
    
