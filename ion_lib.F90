module ion_lib
    use props
    use ptcl_props
    implicit none
    
    contains
    
    ! *** Ion Continuity ***
    subroutine ionEqn(g, i, j, ph, ni, ne, nte, nm, b)
        type(grid), intent(in) :: g
        integer, intent(in)    :: i, j
        real(8), intent(in)    :: ph(:,:), ni(:,:), ne(:,:), nte(:,:), nm(:,:)
        real(8), intent(out)   :: b
        real(8) :: dflx, Te, src, k_ir, k_si
        
        call ionDFlx(g, i, j, ph, ni, dflx)
        
        ! rates and coefficients
        Te = get_Te(nte(i,j), ne(i,j))
        k_ir = get_k_ir(Te)
        k_si = get_k_si(Te)

        ! evaluate source terms
        src =   k_ir * ninf    * ne(i,j) &
              - beta * ni(i,j) * ne(i,j) &
              + k_si * nm(i,j) * ne(i,j) &
              + k_mp * nm(i,j)**2
        
        ! evaluate expression
        b = -dflx + src
    end subroutine

    ! *** Metastable Continuity ***
    subroutine metaEqn(g, i, j, nm, ne, nte, b)
        type(grid), intent(in) :: g
        integer, intent(in)    :: i, j
        real(8), intent(in)    :: nm(:,:), ne(:,:), nte(:,:)
        real(8), intent(out)   :: b
        real(8) :: Te, dflx, src, k_ex, k_sc, k_si, nu
        
        call metaDFlx(g, i, j, nm, dflx)
        
        ! rates and coefficients
        Te = get_Te(nte(i,j), ne(i,j))
        k_sc = get_k_sc(Te)
        k_si = get_k_si(Te)
        k_ex = get_k_ex(Te)
        nu   = get_nu(Te)

        ! evaluate source term
        src =   k_ex * ninf    * ne(i,j) &
              - k_si * nm(i,j) * ne(i,j) &
              - k_sc * nm(i,j) * ne(i,j) &
              - k_r  * nm(i,j) * ne(i,j) &
              - 2d0  * k_mp    * nm(i,j)**2 &
              - k_2q * ninf    * nm(i,j) &
              - k_3q * ninf**2 * nm(i,j)
        
        ! evaluate expression
        b = -dflx + src
    end subroutine
    
    ! *** Calculate Divergence of Ion Flux***
    subroutine ionDFlx(g, i, j, ph, ni, dflx)
        type(grid), intent(in) :: g
        integer, intent(in)    :: i, j
        real(8), intent(in)    :: ph(:,:), ni(:,:)
        real(8), intent(out)   :: dflx
        real(8) :: a, Ex(2), flx_x(2), flx_y(2)
        
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
        
        dflx = (flx_x(2) - flx_x(1)) / g%dlx(i-1)
    end subroutine
    
    ! *** Calculate Divergence of Metastable Flux ***
    subroutine metaDFlx(g, i, j, nm, dflx)
        type(grid), intent(in) :: g
        integer, intent(in) :: i, j
        real(8), intent(in) :: nm(:,:)
        real(8), intent(out) :: dflx
        real(8) :: flx_x(2), flx_y(2)
        
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
        
        dflx = (flx_x(2) - flx_x(1)) / g%dlx(i-1)
    end subroutine
end module
    
