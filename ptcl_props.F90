module ptcl_props
    use props
    implicit none
    
    ! electron properties
    real(8), parameter:: me = 9.10938188e-31

    ! positive ion properties
    real(8), parameter:: mi  = 1.67262178d-27 * 39.948, &
                         mui = 1.45d3 / p * 1d-4 * ph0 * t0 / x0**2, &
                         Ti  = Tg, &
                         vi  = Sqrt( (8d0 * kb*Ti) / (pi * mi) ) * t0 / x0, &
                         Di  = mui * (kb*Ti / e) / ph0

    ! metastable argon properties
    real(8), parameter:: mm   = mi, &
                         Dm   = 2.42d18 * x0 / ninf / 1d2 * t0, &
                         k_r  = 2d-7 / 1d6 / x0**3 * t0, &
                         k_mp = 6.2d-10 / 1.0d6 / x0**3 * t0, &
                         k_2q = 3d-15 / 1d6 / x0**3 * t0, &
                         k_3q = 1.1d-31 / 1.0d12 / x0**6 * t0

    ! reactions
    real(8), parameter:: H_ir  =  15.8  / ph0, &
                         H_ex  =  11.56 / ph0, &
                         H_si  =   4.14 / ph0, &
                         H_sc  = -11.56 / ph0, &
                         beta  =  1e-13 * t0 / x0**3, &
                         gam   =   0.1
    
contains
    
    subroutine getFlx(flux, E, dh, q, mu, D, n_left, n_right)
        real(8), intent(inout) :: flux
        integer, intent(in) :: q
        real(8), intent(in) :: E, dh, mu, D, n_left, n_right
        real(8) :: v, tol, arg
        
        tol = 1e-12
        v = q * mu * E   
        arg = v * dh / D
        
        if (abs(q*E) < tol) then
            flux = D * (n_left - n_right) / dh
        
        ! Positive exponentials blow up,
        !  so rewrite always as negative exp.
        !  this is the same analytical expression
        else if (arg > 0) then
            flux = v * (n_left - n_right * exp(max(-arg,-100d0))) &
                / (1d0 - exp(max(-arg, -100d0)))
        else
            flux = v * (n_right - n_left * exp(max(arg, -100d0))) &
                / (1d0  - exp( max(arg,-100d0)))
        end if
    end subroutine
    
    function get_Te(nte,ne)
        real(8):: get_Te
        real(8), intent(in) :: nte, ne
        
        if ((nte >= 0d0) .and. (ne >= 0d0)) then
            get_Te = nte / ne
        else
            write(*,*) "Error, negative density. Stop."
            stop
            get_Te = 1d-8
        end if
        
        return
    end function
    
    function get_mue(T)
        real(8):: get_mue
        real(8), intent(in):: T
        real(8):: x, &
                  a1 =  58.1133016145,    &
                  b1 =  -0.984082217962,  &
                  c1 =  -0.164770900119,  &
                  d1 =  -0.142058042584,  &
                  f1 =  -0.0637079234081, &
                  g1 =   0.0436642742558, &
                  a2 =  70.7846754548,    &
                  b2 = -22.7558138237,    &
                  c2 =  12.8366493242,    &
                  d2 =  -3.57242763244,   &
                  f2 =   0.486276623664,  &
                  g2 =  -0.0259045697422
        
        x = log(max(2.34d-1, min(1.57d2, T * ph0)))
        
        if (x < log(5.119)) then
            get_mue = exp(a1 + b1*x + c1*x**2 + d1*x**3 + f1*x**4 &
                          + g1*x**5) * x0 / ninf * t0 * ph0
        else
            get_mue = exp(a2 + b2*x + c2*x**2 + d2*x**3 + f2*x**4 &
                          + g2*x**5) * x0 / ninf * t0 * ph0
        end if
        return
    end function get_mue
    
    function get_mut(T)
        real(8):: get_mut
        real(8), intent(in):: T
        real(8):: x, &
                  a1 =  58.1354622038,    &
                  b1 =  -1.29918391109,   &
                  c1 =  -0.106559658509,  &
                  d1 =  -0.0202265418394, &
                  f1 =  -0.0286667333915, &
                  g1 =   0.022457210094,  &
                  a2 =  68.2627663433,    &
                  b2 = -19.3524000831,    &
                  c2 =  11.3251980583,    &
                  d2 =  -3.24070386509,   &
                  f2 =   0.449848156054,  &
                  g2 =  -0.0242585762405
        
        x = log(max(2.34d-1, min(1.57d2, T * ph0)))
        
        if (x < log(5.119)) then
            get_mut = exp(a1 + b1*x + c1*x**2 + d1*x**3 + f1*x**4 &
                          + g1*x**5) * x0 / ninf * t0 * ph0
        else
            get_mut = exp(a2 + b2*x + c2*x**2 + d2*x**3 + f2*x**4 &
                          + g2*x**5) * x0 / ninf * t0 * ph0
        end if
        return
    end function get_mut
    
    function get_De(T)
        real(8):: get_De
        real(8), intent(in):: T
        real(8):: x, &
                  a = 57.7067018813,      &
                  b = -0.0699892875381,   &
                  c = -0.117645585949,    &
                  d =  0.105390295278,    &
                  f = -0.102862612604,    &
                  g = -0.0469171521686,   &
                  h =  0.0584908312121,   &
                  i =  0.000578601715687, &
                  j = -0.0122860884883,   &
                  k =  0.00426793748856,  &
                  l = -0.000590000082557, &
                  m = 3.00706533201e-05
                  
        x = log(max(2.34d-1, min(1.57d2, T * ph0)))
        
        get_De = exp(a + b*x + c*x**2 + d*x**3 + f*x**4 + g*x**5 &
                     + h*x**6 + i*x**7 + j*x**8 + k*x**9 &
                     + l*x**10 + m*x**11) * x0 / ninf * t0
        return
    end function get_De

    function get_Dt(T)
        real(8) :: get_Dt
        real(8), intent(in):: T
        real(8):: x, &
                  a = 57.7214952841,     &
                  b = -0.353731464348,   &
                  c = -0.00505156731795, &
                  d =  0.161173720584,   &
                  f = -0.196610345869,   &
                  g = -0.0719115643218,  &
                  h =  0.11874085868,    &
                  i = -0.00669712724784, &
                  j = -0.0236445308504 , &
                  k =  0.00917326671764, &
                  l = -0.00135453096483, &
                  m =  7.26379684461e-05
                  
        x = log(max(2.34d-1, min(1.57d2, T * ph0)))
        
        get_Dt = exp(a + b*x + c*x**2 + d*x**3 + f*x**4 + g*x**5 &
                     + h*x**6 + i*x**7 + j*x**8 + k*x**9 &
                     + l*x**10 + m*x**11) * x0 / ninf * t0
        return
    end function get_Dt
    
    function get_k_ex(T)
        real(8):: get_k_ex
        real(8), intent(in):: T
        real(8):: x, &
                  a1 =  -54.4513453969,  &
                  b1 =   17.1472529739,  &
                  c1 =   -1.05188982824, &
                  d1 =  -10.5053010711,  &
                  f1 =    7.51502551486, &
                  g1 =   -1.44763942525, &
                  a2 = -150.266369415,   &
                  b2 =  158.67012941,    &
                  c2 =  -85.5163470769,  &
                  d2 =   23.0655081891,  &
                  f2 =   -3.0908688844,  &
                  g2 =    0.16403226902
        
        x = log(max(9.611d-1, min(1.57d2, T * ph0)))
        
        if (x < log(5.667)) then
            get_k_ex = exp(a1 + b1*x + c1*x**2 + d1*x**3 + f1*x**4 &
                          + g1*x**5) * t0 / x0**3
        else
            get_k_ex = exp(a2 + b2*x + c2*x**2 + d2*x**3 + f2*x**4 &
                          + g2*x**5) * t0 / x0**3
        end if
        return
    end function get_k_ex
    
    function get_k_ir(T)
        real(8):: get_k_ir
        real(8), intent(in):: T
        real(8):: x, &
                  a1 = -73.27080258,   &
                  b1 = -30.1312514721, &
                  c1 = 138.743842326,  &
                  d1 = -161.378167125, &
                  f1 = 81.9325486193,  &
                  g1 = -14.9080865672, &
                  a2 = -168.925304837, &
                  b2 = 177.113452524,  &
                  c2 = -92.3047166553, &
                  d2 = 24.1814225482,  &
                  f2 = -3.15463258104, &
                  g2 = 0.16331188192
        
        x = log(max(1.37d0, min(1.57d2, T * ph0)))
        
        if (x < log(7.5)) then
            get_k_ir = exp(a1 + b1*x + c1*x**2 + d1*x**3 + f1*x**4 &
                          + g1*x**5) * t0 / x0**3
        else
            get_k_ir = exp(a2 + b2*x + c2*x**2 + d2*x**3 + f2*x**4 &
                          + g2*x**5) * t0 / x0**3
        end if
        return
    end function get_k_ir
    
    function get_nu(T)
        real(8):: get_nu
        real(8), intent(in):: T
        real(8):: x, &
                  a = -32.275912575,      &
                  b =   1.45173283977,    &
                  c =   0.00936933121094, &
                  d =   0.129397015353,   &
                  f =  -0.0414865809044,  &
                  g =  -0.0582934303409,  &
                  h =   0.0309832277826,  &
                  i =  -0.00542014733763,  &
                  j =   0.000325615321708
                  
        x = log(max(2.34d-1, min(1.57d2, T * ph0)))
        
        get_nu = exp(a + b*x + c*x**2. + d*x**3. + f*x**4. + g*x**5. &
                     + h*x**6. + i*x**7 + j*x**8.) / x0**3 * ninf * t0
        return
    end function get_nu
    
    function get_k_sc(T)
        real(8):: get_k_sc
        real(8), intent(in):: T
        real(8):: x, &
                  a = -21.4827864151,   &
                  b =   0.457356923276, &
                  c =  -0.555439231606, &
                  d =   1.27257798891,  &
                  f =  -0.67840685073,  &
                  g =   0.10591014464
                  
        x = log(min(16d0, max(5d-1, T * ph0)))
        
        get_k_sc = exp(a + b*x + c*x**2. + d*x**3. + f*x**4. &
                       + g*x**5. ) / 1.0d6 / x0**3 * t0
        return
    end function get_k_sc

    function get_k_si(T)
        real(8):: get_k_si
        real(8), intent(in):: T
        real(8):: x, &
                  a = -43.1347385848, &
                  b =  43.9905424566, &
                  c = -28.1169537586, &
                  d =   8.28853856817, &
                  f =  -0.931626144207
                  
        x = log(min(16d0, max(5d-1, T * ph0)))
        
        get_k_si = exp(a + b*x + c*x**2. + d*x**3. + f*x**4.) / 1.0d6 / x0**3 * t0
        return
    end function get_k_si
end module
