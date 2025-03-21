module ice_enthalpy 
    ! Module contains the ice temperature and basal mass balance (grounded) solution

    use yelmo_defs, only : wp, pi  
    use solver_tridiagonal, only : solve_tridiag 
    use thermodynamics, only : calc_bmb_grounded, calc_advec_vertical_column, &
                               convert_to_enthalpy, convert_from_enthalpy_column

    !use interp1D 

    implicit none
    
    private
    public :: calc_temp_column
    public :: calc_temp_bedrock_column
    public :: calc_enth_column
    public :: calc_dzeta_terms
    public :: calc_zeta_twolayers
    public :: calc_zeta_combined
    public :: get_cts_index

contains 
    
    subroutine calc_temp_column(enth,T_ice,omega,bmb_grnd,Q_ice_b,H_cts,T_pmp,cp,kt,advecxy,uz, &
                                Q_strn,Q_b,Q_rock,T_srf,T_shlf,H_ice,H_w,f_grnd,zeta_aa,zeta_ac, &
                                dzeta_a,dzeta_b,omega_max,T0,rho_ice,rho_w,L_ice,sec_year,dt)
        ! Thermodynamics solver for a given column of ice 
        ! Note zeta=height, k=1 base, k=nz surface 
        ! Note: nz = number of vertical boundaries (including zeta=0.0 and zeta=1.0), 
        ! temperature is defined for cell centers, plus a value at the surface and the base
        ! so nz_ac = nz_aa - 1 

        ! For notes on implicit form of advection terms, see eg http://farside.ph.utexas.edu/teaching/329/lectures/node90.html
        
        implicit none 

        real(wp), intent(INOUT) :: enth(:)        ! nz_aa [J kg] Ice column enthalpy
        real(wp), intent(INOUT) :: T_ice(:)       ! nz_aa [K] Ice column temperature
        real(wp), intent(INOUT) :: omega(:)       ! nz_aa [-] Ice column water content fraction
        real(wp), intent(INOUT) :: bmb_grnd       ! [m a-1] Basal mass balance (melting is negative)
        real(wp), intent(OUT)   :: Q_ice_b        ! [mW m-2] Ice basal heat flux (positive up)
        real(wp), intent(OUT)   :: H_cts          ! [m] cold-temperate transition surface (CTS) height
        real(wp), intent(IN)    :: T_pmp(:)       ! nz_aa [K] Pressure melting point temp.
        real(wp), intent(IN)    :: cp(:)          ! nz_aa [J kg-1 K-1] Specific heat capacity
        real(wp), intent(IN)    :: kt(:)          ! nz_aa [J a-1 m-1 K-1] Heat conductivity 
        real(wp), intent(IN)    :: advecxy(:)     ! nz_aa [K a-1] Horizontal heat advection 
        real(wp), intent(IN)    :: uz(:)          ! nz_ac [m a-1] Vertical velocity 
        real(wp), intent(IN)    :: Q_strn(:)      ! nz_aa [J a-1 m-3] Internal strain heat production in ice
        real(wp), intent(IN)    :: Q_b            ! [mW m-2] Basal frictional heat production
        real(wp), intent(IN)    :: Q_rock         ! [mW m-2] Bedrock heat flux (positive up)
        real(wp), intent(IN)    :: T_srf          ! [K] Surface temperature 
        real(wp), intent(IN)    :: T_shlf         ! [K] Marine-shelf interface temperature
        real(wp), intent(IN)    :: H_ice          ! [m] Ice thickness 
        real(wp), intent(IN)    :: H_w            ! [m] Basal water layer thickness 
        real(wp), intent(IN)    :: f_grnd         ! [--] Grounded fraction
        real(wp), intent(IN)    :: zeta_aa(:)     ! nz_aa [--] Vertical sigma coordinates (zeta==height), layer centered aa-nodes
        real(wp), intent(IN)    :: zeta_ac(:)     ! nz_ac [--] Vertical height axis temperature (0:1), layer edges ac-nodes
        real(wp), intent(IN)    :: dzeta_a(:)     ! nz_aa [--] Solver discretization helper variable ak
        real(wp), intent(IN)    :: dzeta_b(:)     ! nz_aa [--] Solver discretization helper variable bk
        real(wp), intent(IN)    :: omega_max      ! [-] Maximum allowed water fraction inside ice, typically omega_max=0.02 
        real(wp), intent(IN)    :: T0             ! [K or degreesCelcius] Reference melting temperature  
        real(wp), intent(IN)    :: rho_ice 
        real(wp), intent(IN)    :: rho_w
        real(wp), intent(IN)    :: L_ice
        real(wp), intent(IN)    :: sec_year 
        real(wp), intent(IN)    :: dt             ! [a] Time step  

        ! Local variables 
        integer  :: k, nz_aa, nz_ac
        real(wp) :: H_w_predicted
        real(wp) :: dz, dz1, dz2, d2Tdz2 
        real(wp) :: T00, T01, T02, zeta_now  
        real(wp) :: T_excess
        real(wp) :: Q_ice_b_now, Q_b_now, Q_rock_now 
        real(wp) :: melt_internal
        real(wp) :: val_base, val_srf 
        logical  :: is_basal_flux  
        logical  :: is_surf_flux 

        real(wp), allocatable :: kappa_aa(:)    ! aa-nodes
        real(wp), allocatable :: Q_strn_now(:)  ! aa-nodes
        
        real(wp), parameter   :: T_ref        = 273.15_wp   
        real(wp), parameter   :: T_min_lim    = 200.00_wp 
        real(wp), parameter   :: bmb_grnd_lim = 10.0_wp     ! [m/yr]
        
        nz_aa = size(zeta_aa,1)

        allocate(kappa_aa(nz_aa))
        allocate(Q_strn_now(nz_aa))

        ! Get bedrock heat flux in proper units 
        Q_rock_now = Q_rock*1e-3*sec_year   ! [mW m-2] => [J m-2 a-1]

        ! Get basal frictional heat flux in proper units 
        Q_b_now = Q_b*1e-3*sec_year         ! [mW m-2] => [J m-2 a-1]

        ! Step 0: Calculate diffusivity on cell centers (aa-nodes)

        kappa_aa = kt / (rho_ice*cp)
        
        ! Convert units of Q_strn [J a-1 m-3] => [K a-1]
        Q_strn_now = Q_strn/(rho_ice*cp)

        ! === Surface boundary condition =====================

        val_srf      =  min(T_srf,T0)    
        is_surf_flux = .FALSE. 

        ! === Basal boundary condition =====================

        if (f_grnd .lt. 1.0) then
            ! Floating or partially floating ice - set temperature equal 
            ! to basal temperature at pressure melting point, or marine freezing temp,
            ! or weighted average between the two.
            
            val_base      = (f_grnd*T_pmp(1) + (1.0-f_grnd)*T_shlf)
            is_basal_flux = .FALSE. 

        else 
            ! Grounded ice 

            ! Determine expected basal water thickness [m] for this timestep,
            ! using basal mass balance from previous time step (good guess)
            H_w_predicted = H_w - (bmb_grnd*(rho_w/rho_ice))*dt
            
            ! == Assign grounded basal boundary conditions ==

            if (H_w_predicted .gt. 0.0_wp) then 
                ! Temperate at bed 
                ! Hold basal temperature at pressure melting point

                val_base      = T_pmp(1)
                is_basal_flux = .FALSE.

            else if ( T_ice(1) .lt. T_pmp(1) .or. H_w_predicted .lt. 0.0_wp ) then
                ! Frozen at bed, or about to become frozen 

                ! backward Euler flux basal boundary condition
                val_base      = -(Q_b_now + Q_rock_now) / kt(1)
                is_basal_flux = .TRUE. 
                
            else 
                ! Temperate at bed 
                ! Hold basal temperature at pressure melting point

                val_base      = T_pmp(1)
                is_basal_flux = .FALSE. 
                
            end if   ! melting or frozen

        end if  ! floating or grounded 

        ! === Solver =============================
     
        call calc_temp_column_internal(T_ice,kappa_aa,uz,advecxy,Q_strn_now,val_base,val_srf,H_ice, &
                                                zeta_aa,zeta_ac,dzeta_a,dzeta_b,T_ref,dt, &
                                                is_basal_flux,is_surf_flux)

        ! Now calculate internal melt (only allow melting, no accretion)
    
        melt_internal = 0.0 

        do k = nz_aa-1, 2, -1 
            ! Descend from surface to base layer (center of layer)

            ! Store temperature difference above pressure melting point (excess energy)
            T_excess = max(T_ice(k)-T_pmp(k),0.0)

            ! Calculate basal mass balance as sum of all water produced in column,
            ! reset temperature to pmp  
            if (T_excess .gt. 0.0) then 
                melt_internal = melt_internal + T_excess * H_ice*(zeta_ac(k)-zeta_ac(k-1))*cp(k) / (L_ice * dt) 
                T_ice(k)      = T_pmp(k)
            end if 
            
        end do 

        ! Make sure base is below pmp too (mass/energy balance handled via bmb_grnd calculation externally)
        if (T_ice(1) .gt. T_pmp(1)) T_ice(1) = T_pmp(1)

        ! Also ensure there are no crazy low values (can happen with bad initialization)
        where (T_ice(:) .lt. T_min_lim) T_ice(:) = T_min_lim
        
        ! ajr: some kind of limit on the temp gradient near the base is
        ! probably a good idea. But the below method wouldn't work well, if 
        ! there are only eg nz=3 points in the vertical. Needs more thought...
        ! ! Finally check temperature gradient near the base to avoid potential instabilities 
        ! if (T_ice(1) .eq. T_pmp(1) .and. T_ice(3) .gt. T_ice(2)) then 
        !     ! Change in slope of T_ice profile detected near base - likely
        !     ! the temperature solver did not converge properly. Correct
        !     ! temperature just above the base with linear interpolation

        !     dz = (zeta_aa(2)-zeta_aa(1))
        !     T_ice(2) = T_ice(1) + dz*((T_ice(3)-T_ice(1))/(zeta_aa(3)-zeta_aa(1)))
        ! end if

        ! Also set omega to constant value where ice is temperate just for some consistency 
        omega = 0.0 
!             where (T_ice .ge. T_pmp) omega = omega_max 

        ! Finally, get enthalpy too 
        call convert_to_enthalpy(enth,T_ice,omega,T_pmp,cp,L_ice)

        ! Calculate heat flux at ice base as temperature gradient * conductivity [J a-1 m-2]
        if (H_ice .gt. 0.0_wp) then 

            ! 1st order, upwind gradient dTdz 
            ! Works, but can cause oscillations in H_w 
            dz = H_ice * (zeta_aa(2)-zeta_aa(1))
            Q_ice_b_now = -kt(1) * (T_ice(2) - T_ice(1)) / dz 

        else 

            Q_ice_b_now = 0.0  

        end if 
        
        ! Calculate Q_ice_b for global output 
        Q_ice_b = Q_ice_b_now*1e3/sec_year      ! [J a-1 m-2] => [mW m-2]

        
        if (f_grnd .gt. 0.0) then 
            ! Calculate basal mass balance (valid for grounded ice only)

            call calc_bmb_grounded(bmb_grnd,T_ice(1)-T_pmp(1),Q_ice_b_now,Q_b_now,Q_rock_now,rho_ice,L_ice)

            ! Include internal melting in bmb_grnd (not allowed for floating ice, but not expected either)
            bmb_grnd = bmb_grnd - melt_internal 

            ! Safety: limit bmb_grnd to reasonable values to avoid problems
            ! (grounded ice melt should be much less than this limit, eg 10 m/yr)
            if (bmb_grnd .gt. bmb_grnd_lim)  bmb_grnd =  bmb_grnd_lim
            if (bmb_grnd .lt. -bmb_grnd_lim) bmb_grnd = -bmb_grnd_lim

        else 
            ! No grounded bmb allowed for floating ice 

            bmb_grnd = 0.0_wp 

        end if

        ! Finally, calculate the CTS height 
        H_cts = calc_cts_height(enth,T_ice,omega,T_pmp,cp,H_ice,zeta_aa)

        return 

    end subroutine calc_temp_column

    subroutine calc_temp_bedrock_column(enth,temp,Q_rock,cp,kt,Q_ice_b,Q_geo,T_srf,H_rock, &
                                                zeta_aa,zeta_ac,dzeta_a,dzeta_b,rho_rock,sec_year,dt)
        ! Thermodynamics solver for a given column of ice 
        ! Note zeta=height, k=1 base, k=nz surface 
        ! Note: nz = number of vertical boundaries (including zeta=0.0 and zeta=1.0), 
        ! temperature is defined for cell centers, plus a value at the surface and the base
        ! so nz_ac = nz_aa - 1 

        ! For notes on implicit form of advection terms, see eg http://farside.ph.utexas.edu/teaching/329/lectures/node90.html
        
        implicit none 

        real(wp), intent(INOUT) :: enth(:)        ! nz_aa [J kg] Ice column enthalpy
        real(wp), intent(INOUT) :: temp(:)        ! nz_aa [K] Ice column temperature
        real(wp), intent(OUT)   :: Q_rock         ! [mW m-2] Bed surface heat flux (positive up)
        real(wp), intent(IN)    :: cp             ! nz_aa [J kg-1 K-1] Specific heat capacity
        real(wp), intent(IN)    :: kt             ! nz_aa [J a-1 m-1 K-1] Heat conductivity 
        real(wp), intent(IN)    :: Q_ice_b        ! [mW m-2] Ice basal heat flux (positive up)
        real(wp), intent(IN)    :: Q_geo          ! [mW m-2] Bedrock heat flux (positive up)
        real(wp), intent(IN)    :: T_srf          ! [K] Surface temperature 
        real(wp), intent(IN)    :: H_rock         ! [m] Bedrock thickness 
        real(wp), intent(IN)    :: zeta_aa(:)     ! nz_aa [--] Vertical sigma coordinates (zeta==height), layer centered aa-nodes
        real(wp), intent(IN)    :: zeta_ac(:)     ! nz_ac [--] Vertical height axis temperature (0:1), layer edges ac-nodes
        real(wp), intent(IN)    :: dzeta_a(:)     ! nz_aa [--] Solver discretization helper variable ak
        real(wp), intent(IN)    :: dzeta_b(:)     ! nz_aa [--] Solver discretization helper variable bk
        real(wp), intent(IN)    :: rho_rock
        real(wp), intent(IN)    :: sec_year 
        real(wp), intent(IN)    :: dt             ! [a] Time step 

        ! Local variables 
        integer  :: k, nz_aa, nz_ac
        real(wp) :: Q_rock_now 
        real(wp) :: Q_ice_b_now 
        real(wp) :: Q_geo_now 
        real(wp) :: dz 
        real(wp) :: val_base, val_srf 
        logical  :: is_basal_flux  
        logical  :: is_surf_flux  

        real(wp), allocatable :: kappa_aa(:)      ! aa-nodes
        
        real(wp), allocatable, target :: zeros(:)
        real(wp), pointer     :: advecxy(:)
        real(wp), pointer     :: Q_strn(:)
        real(wp), pointer     :: uz(:)
        
        real(wp), parameter   :: T_ref = 273.15_wp   

        nz_aa = size(zeta_aa,1)
        nz_ac = size(zeta_ac,1) 

        ! Get basal heat flux from ice in proper units 
        Q_ice_b_now = Q_ice_b*1e-3*sec_year     ! [mW m-2] => [J m-2 a-1]

        ! Get deep bedrock heat flux in proper units 
        Q_geo_now = Q_geo*1e-3*sec_year   ! [mW m-2] => [J m-2 a-1]

        ! Step 0: Calculate diffusivity on cell centers (aa-nodes)
        
        allocate(kappa_aa(nz_aa))

        kappa_aa = kt / (rho_rock*cp)
        
        ! Set unused terms to zero (via pointers to save allocations)
        allocate(zeros(nz_ac))
        zeros = 0.0_wp 
        
        advecxy => zeros(1:nz_aa) 
        Q_strn  => zeros(1:nz_aa) 
        uz      => zeros(1:nz_ac) 

        ! === Surface boundary condition =====================

        ! if ( T_srf .lt. T_pmp_srf ) then
        !         ! Frozen at bed

        !     ! backward Euler flux surface boundary condition
        !     val_srf      =  -Q_ice_b_now / kt
        !     is_surf_flux = .TRUE.   

        ! else 
            ! Temperate at bed surface
            ! Hold bed surface temperature at ice base temperature

            val_srf      = T_srf
            is_surf_flux = .FALSE. 
            
        ! end if
        
        ! === Basal boundary condition =====================

        ! backward Euler flux basal boundary condition
        val_base = -Q_geo_now / kt
        is_basal_flux = .TRUE. 
        
        ! === Solver =============================
     
        call calc_temp_column_internal(temp,kappa_aa,uz,advecxy,Q_strn,val_base,val_srf,H_rock, &
                                                zeta_aa,zeta_ac,dzeta_a,dzeta_b,T_ref,dt, &
                                                is_basal_flux,is_surf_flux)

        ! Convert to enthalpy too 
        call convert_to_enthalpy(enth,temp,omega=0.0_wp,T_pmp=0.0_wp,cp=cp,L=0.0_wp)

        ! Calculate heat flux at bedrock surface as temperature gradient * conductivity [J a-1 m-2]

        ! 1st order, upwind gradient dTdz 
        ! Works, but can cause oscillations in H_w 
        dz = H_rock * (zeta_aa(nz_aa)-zeta_aa(nz_aa-1))
        Q_rock_now = -kt * (temp(nz_aa) - temp(nz_aa-1)) / dz 

        ! Set global units to [mW m-2]
        Q_rock = Q_rock_now *1e3 / sec_year 
        
        return 

    end subroutine calc_temp_bedrock_column

    subroutine calc_temp_column_internal(temp,kappa,uz,advecxy,Q_strn,val_base,val_srf,thickness, &
                                                zeta_aa,zeta_ac,dzeta_a,dzeta_b,T_ref,dt,is_basal_flux,is_surf_flux)
        ! Thermodynamics solver 1D (eg, for column of ice or bedrock)
        ! Note zeta= column height, k=1 base, k=nz surface 
        ! Note: nz = number of vertical boundaries (including zeta=0.0 and zeta=1.0), 
        ! temperature is defined for cell centers, plus a value at the surface and the base
        ! so nz_ac = nz_aa - 1 

        ! For notes on implicit form of advection terms, see eg http://farside.ph.utexas.edu/teaching/329/lectures/node90.html
        
        implicit none 

        real(wp), intent(INOUT) :: temp(:)        ! nz_aa [K] Ice column temperature
        real(wp), intent(IN)    :: kappa(:)       ! nz_aa [] Diffusivity
        real(wp), intent(IN)    :: uz(:)          ! nz_ac [m a-1] Vertical velocity 
        real(wp), intent(IN)    :: advecxy(:)     ! nz_aa [K a-1] Horizontal heat advection 
        real(wp), intent(IN)    :: Q_strn(:)      ! nz_aa [K a-1] Internal strain heat production in ice
        real(wp), intent(IN)    :: val_base       ! [K or flux] Basal boundary condition
        real(wp), intent(IN)    :: val_srf        ! [K or flux] Surface temperature 
        real(wp), intent(IN)    :: thickness      ! [m] Total column thickness 
        real(wp), intent(IN)    :: zeta_aa(:)     ! nz_aa [--] Vertical sigma coordinates (zeta==height), layer centered aa-nodes
        real(wp), intent(IN)    :: zeta_ac(:)     ! nz_ac [--] Vertical height axis temperature (0:1), layer edges ac-nodes
        real(wp), intent(IN)    :: dzeta_a(:)     ! nz_aa [--] Solver discretization helper variable ak
        real(wp), intent(IN)    :: dzeta_b(:)     ! nz_aa [--] Solver discretization helper variable bk
        real(wp), intent(IN)    :: T_ref          ! [K] Reference temperature to scale calculation
        real(wp), intent(IN)    :: dt             ! [a] Time step 
        logical,  intent(IN)    :: is_basal_flux  ! Is basal bc flux condition (True) or Neumann (False)
        logical,  intent(IN)    :: is_surf_flux   ! Is surf  bc flux condition (True) or Neumann (False)
        
        ! Local variables 
        integer  :: k, nz_aa
        real(wp) :: fac, fac_a, fac_b, uz_aa, dz, dzeta
        real(wp) :: h1, h2, afac_a, afac_b, afac_mid
        real(wp) :: kappa_a, kappa_b, dz1, dz2 
        real(wp), allocatable :: subd(:)      ! nz_aa 
        real(wp), allocatable :: diag(:)      ! nz_aa  
        real(wp), allocatable :: supd(:)      ! nz_aa 
        real(wp), allocatable :: rhs(:)       ! nz_aa 
        real(wp), allocatable :: solution(:)  ! nz_aa
        
        nz_aa = size(zeta_aa,1)

        allocate(subd(nz_aa))
        allocate(diag(nz_aa))
        allocate(supd(nz_aa))
        allocate(rhs(nz_aa))
        allocate(solution(nz_aa))

        ! == Ice base ==

        if (is_basal_flux) then 
            ! Impose basal flux (Neumann condition)

            ! Calculate dz for the bottom layer between the basal boundary
            ! (ac-node) and the centered (aa-node) temperature point above
            ! Note: zeta_aa(1) == zeta_ac(1) == bottom boundary 
            dz = thickness * (zeta_aa(2) - zeta_aa(1))

            ! backward Euler flux basal boundary condition
            subd(1) =  0.0_wp
            diag(1) = -1.0_wp
            supd(1) =  1.0_wp
            rhs(1)  = val_base * dz
                
        else 
            ! Impose basal temperature (Dirichlet condition)  

            subd(1) = 0.0_wp
            diag(1) = 1.0_wp
            supd(1) = 0.0_wp
            rhs(1)  = (val_base - T_ref)

        end if 

        ! == Ice interior layers 2:nz_aa-1 ==

        do k = 2, nz_aa-1

            ! Get kappa for the lower and upper ac-nodes using harmonic mean from aa-nodes
            
            dz1 = zeta_ac(k)-zeta_aa(k-1)
            dz2 = zeta_aa(k)-zeta_ac(k)
            call calc_wtd_harmonic_mean(kappa_a,kappa(k-1),kappa(k),dz1,dz2)

            dz1 = zeta_ac(k+1)-zeta_aa(k)
            dz2 = zeta_aa(k+1)-zeta_ac(k+1)
            call calc_wtd_harmonic_mean(kappa_b,kappa(k),kappa(k+1),dz1,dz2)

            ! Get diffusion factors
            fac_a   = -kappa_a*dzeta_a(k)*dt/thickness**2
            fac_b   = -kappa_b*dzeta_b(k)*dt/thickness**2

            ! Get implicit vertical advection term, ac => aa nodes
            uz_aa   = 0.5*(uz(k)+uz(k+1))

if (.TRUE.) then
    ! Use simple centered-difference scheme for vertical advection

            ! Vertical distance for centered difference advection scheme
            dzeta   =  (zeta_aa(k+1)-zeta_aa(k-1))
            dz      = thickness * dzeta

            subd(k) = fac_a - uz_aa * dt/dz
            supd(k) = fac_b + uz_aa * dt/dz
            diag(k) = 1.0_wp - fac_a - fac_b
            rhs(k)  = (temp(k)-T_ref) - dt*advecxy(k) + dt*Q_strn(k)

else
    ! Use centered-difference advection scheme accounting for variable layer thickness
    ! Note: this can get unstable in some rare cases, for now using centered differences as above. 
    
            ! Get grid spacing between neighboring points
            h1 = thickness * (zeta_aa(k)-zeta_aa(k-1))
            h2 = thickness * (zeta_aa(k+1)-zeta_aa(k))

            ! Get advective factor terms for centered difference with 
            ! uneven layers 
            afac_a   = -h2/(h1*(h1+h2))
            afac_mid = -(h1-h2)/(h1*h2)
            afac_b   = +h1/(h2*(h1+h2))

            subd(k) = fac_a + uz_aa * dt*afac_a
            supd(k) = fac_b + uz_aa * dt*afac_b
            diag(k) = 1.0_wp - fac_a - fac_b + uz_aa * dt*afac_mid
            rhs(k)  = (temp(k)-T_ref) - dt*advecxy(k) + dt*Q_strn(k)

end if 

        end do 

        ! == Column surface ==

        if (is_surf_flux) then 
            ! Impose surface flux (Neumann condition)

            ! Calculate dz for the top layer between the surface boundary
            ! (ac-node) and the centered (aa-node) temperature point above
            ! Note: zeta_aa(nz_aa) == zeta_ac(nz_ac) == top boundary 
            dz = thickness * (zeta_aa(nz_aa) - zeta_aa(nz_aa-1))
            
            ! backward Euler flux surface boundary condition
            subd(nz_aa) = -1.0_wp
            diag(nz_aa) =  1.0_wp
            supd(nz_aa) =  0.0_wp
            rhs(nz_aa)  = val_srf * dz
            
        else 
            ! Impose surface temperature (Dirichlet condition) 

            subd(nz_aa) = 0.0_wp
            diag(nz_aa) = 1.0_wp
            supd(nz_aa) = 0.0_wp
            rhs(nz_aa)  = (val_srf-T_ref)

        end if 

        ! == Call solver ==

        call solve_tridiag(subd,diag,supd,rhs,solution)

        ! Copy the solution into the temperature variable

        temp = solution + T_ref 

        return 

    end subroutine calc_temp_column_internal

    subroutine calc_enth_column(enth,T_ice,omega,bmb_grnd,Q_ice_b,H_cts,T_pmp,cp,kt,advecxy,uz, &
                                Q_strn,Q_b,Q_lith,T_srf,T_shlf,H_ice,H_w,f_grnd,zeta_aa,zeta_ac, &
                                dzeta_a,dzeta_b,cr,omega_max,T0,rho_ice,rho_w,L_ice,sec_year,dt)
        ! Thermodynamics solver for a given column of ice 
        ! Note zeta=height, k=1 base, k=nz surface 
        ! Note: nz = number of vertical boundaries (including zeta=0.0 and zeta=1.0), 
        ! temperature is defined for cell centers, plus a value at the surface and the base
        ! so nz_ac = nz_aa - 1 

        ! For notes on implicit form of advection terms, see eg http://farside.ph.utexas.edu/teaching/329/lectures/node90.html
        
        implicit none 

        real(wp), intent(INOUT) :: enth(:)        ! nz_aa [J kg] Ice column enthalpy
        real(wp), intent(INOUT) :: T_ice(:)       ! nz_aa [K] Ice column temperature
        real(wp), intent(INOUT) :: omega(:)       ! nz_aa [-] Ice column water content fraction
        real(wp), intent(INOUT) :: bmb_grnd       ! [m a-1] Basal mass balance (melting is negative)
        real(wp), intent(OUT)   :: Q_ice_b        ! [J a-1 m-2] Ice basal heat flux (positive up)
        real(wp), intent(OUT)   :: H_cts          ! [m] cold-temperate transition surface (CTS) height
        real(wp), intent(IN)    :: T_pmp(:)       ! nz_aa [K] Pressure melting point temp.
        real(wp), intent(IN)    :: cp(:)          ! nz_aa [J kg-1 K-1] Specific heat capacity
        real(wp), intent(IN)    :: kt(:)          ! nz_aa [J a-1 m-1 K-1] Heat conductivity 
        real(wp), intent(IN)    :: advecxy(:)     ! nz_aa [K a-1] Horizontal heat advection 
        real(wp), intent(IN)    :: uz(:)          ! nz_ac [m a-1] Vertical velocity 
        real(wp), intent(IN)    :: Q_strn(:)      ! nz_aa [J a-1 m-3] Internal strain heat production in ice
        real(wp), intent(IN)    :: Q_b            ! [J a-1 m-2] Basal frictional heat production
        real(wp), intent(IN)    :: Q_lith         ! [J a-1 m-2] Bedrock heat flux (positive up)
        real(wp), intent(IN)    :: T_srf          ! [K] Surface temperature 
        real(wp), intent(IN)    :: T_shlf         ! [K] Marine-shelf interface temperature
        real(wp), intent(IN)    :: H_ice          ! [m] Ice thickness 
        real(wp), intent(IN)    :: H_w            ! [m] Basal water layer thickness 
        real(wp), intent(IN)    :: f_grnd         ! [--] Grounded fraction
        real(wp), intent(IN)    :: zeta_aa(:)     ! nz_aa [--] Vertical sigma coordinates (zeta==height), layer centered aa-nodes
        real(wp), intent(IN)    :: zeta_ac(:)     ! nz_ac [--] Vertical height axis temperature (0:1), layer edges ac-nodes
        real(wp), intent(IN)    :: dzeta_a(:)     ! nz_aa [--] Solver discretization helper variable ak
        real(wp), intent(IN)    :: dzeta_b(:)     ! nz_aa [--] Solver discretization helper variable bk
        real(wp), intent(IN)    :: cr             ! [--] Conductivity ratio (kappa_water / kappa_ice)
        real(wp), intent(IN)    :: omega_max      ! [-] Maximum allowed water fraction inside ice, typically omega_max=0.02 
        real(wp), intent(IN)    :: T0             ! [K or degreesCelcius] Reference melting temperature 
        real(wp), intent(IN)    :: rho_ice 
        real(wp), intent(IN)    :: rho_w
        real(wp), intent(IN)    :: L_ice
        real(wp), intent(IN)    :: sec_year  
        real(wp), intent(IN)    :: dt             ! [a] Time step 

        ! Local variables 
        integer  :: k, nz_aa, nz_ac
        integer  :: k_cts  
        real(wp) :: H_w_predicted
        real(wp) :: dz 
        real(wp) :: omega_excess
        real(wp) :: melt_internal
        real(wp) :: val_base, val_srf 
        logical  :: is_basal_flux  

        real(wp), allocatable :: kappa_aa(:)    ! aa-nodes
        real(wp), allocatable :: Q_strn_now(:)  ! aa-nodes
        real(wp), allocatable :: enth_pmp(:)    ! aa-nodes
        
        real(wp), parameter   :: enth_ref = 273.15_wp * 2009.0_wp    ! [K] * [J kg-1 K-1]  

        nz_aa = size(zeta_aa,1)

        allocate(kappa_aa(nz_aa))
        allocate(Q_strn_now(nz_aa))
        allocate(enth_pmp(nz_aa))

        ! Get geothermal heat flux in proper units 
        ! Q_geo_now = Q_geo*1e-3*sec_year   ! [mW m-2] => [J m-2 a-1]

        ! Get enthalpy of the pressure melting point 
        enth_pmp = T_pmp*cp 
        
        ! Find height of CTS - highest temperate layer 
        k_cts = get_cts_index(enth,enth_pmp)

        ! Calculate diffusivity on cell centers (aa-nodes)
!         kappa_aa = kt / (rho_ice*cp)
        call calc_enth_diffusivity(kappa_aa,enth,enth_pmp,cp,kt,cr,rho_ice)

        ! Convert units of Q_strn [J a-1 m-3] => [J kg a-1]
        Q_strn_now = Q_strn/(rho_ice)

        ! === Surface boundary condition =====================

        val_srf =  min(T_srf,T0) * cp(nz_aa)  

        ! === Basal boundary condition =====================

        if (f_grnd .lt. 1.0) then
            ! Floating or partially floating ice - set temperature equal 
            ! to basal temperature at pressure melting point, or marine freezing temp,
            ! or weighted average between the two.
            
            val_base = (f_grnd*T_pmp(1) + (1.0-f_grnd)*T_shlf) * cp(1) 
            is_basal_flux = .FALSE. 

        else 
            ! Grounded ice 

            ! Determine expected basal water thickness [m] for this timestep,
            ! using basal mass balance from previous time step (good guess)
            H_w_predicted = H_w - (bmb_grnd*(rho_w/rho_ice))*dt 
            
            ! == Assign grounded basal boundary conditions ==

            if (H_w_predicted .gt. 0.0_wp) then 

                val_base = enth_pmp(1)
                is_basal_flux = .FALSE.

            else if ( enth(1) .lt. enth_pmp(1) .or. H_w_predicted .lt. 0.0_wp ) then
                ! Frozen at bed, or about to become frozen 

                ! backward Euler flux basal boundary condition
                val_base = (Q_b + Q_lith) / kt(1) * cp(1)
                is_basal_flux = .TRUE. 
                
            else 
                ! Temperate at bed 
                ! Hold basal temperature at pressure melting point

                val_base = enth_pmp(1)
                is_basal_flux = .FALSE. 
                
            end if   ! melting or frozen

        end if  ! floating or grounded 

        ! === Solver =============================
     
        call calc_enth_column_internal(enth,kappa_aa,uz,advecxy,Q_strn_now,val_base,val_srf,H_ice, &
                                            zeta_aa,zeta_ac,dzeta_a,dzeta_b,enth_ref,dt,k_cts,is_basal_flux)


        ! Modify enthalpy at the base in the case that a temperate layer is present above the base
        ! (water content should increase towards the base)
        ! This should come out of routine, but it helps ensure stability to check it here
        if (enth(2) .ge. enth_pmp(2)) enth(1) = enth(2)
        
        ! Get temperature and water content 
        call convert_from_enthalpy_column(enth,T_ice,omega,T_pmp,cp,L_ice)

        ! Set internal melt to zero 
        melt_internal = 0.0 

        do k = nz_aa-1, 2, -1 
            ! Descend from surface to base layer (center of layer)

            ! Store excess water above maximum allowed limit
            omega_excess = max(omega(k)-omega_max,0.0)

            ! Calculate internal melt as sum of all excess water produced in the column 
            if (omega_excess .gt. 0.0) then 
                dz = H_ice*(zeta_ac(k)-zeta_ac(k-1))
                melt_internal = melt_internal + (omega_excess*dz) / dt 
                omega(k)      = omega_max 
            end if 

        end do 

        ! Also limit basal omega to omega_max (even though it doesn't have thickness)
        if (omega(1) .gt. omega_max) omega(1) = omega_max 

        ! Finally, get enthalpy again too (to be consistent with new omega) 
        call convert_to_enthalpy(enth,T_ice,omega,T_pmp,cp,L_ice)

        ! Calculate heat flux at ice base as enthalpy gradient * rho_ice * diffusivity [J a-1 m-2]
        if (H_ice .gt. 0.0_wp) then 
            dz = H_ice * (zeta_aa(2)-zeta_aa(1))
!             Q_ice_b = kappa_aa(1) * rho_ice * (enth(2) - enth(1)) / dz    ! <== Problematic
            Q_ice_b = kt(1) * ( T_ice(2) - T_ice(1) ) / dz
        else
            Q_ice_b = 0.0_wp 
        end if 

        ! Calculate basal mass balance 
        write(*,*) "calc_enth_column:: routine needs to be updated with Q_rock etc."
        stop 

        ! call calc_bmb_grounded(bmb_grnd,T_ice(1)-T_pmp(1),Q_ice_b,Q_b,Q_lith,f_grnd,rho_ice,L_ice) 
!         call calc_bmb_grounded_enth(bmb_grnd,T_ice(1)-T_pmp(1),omega(1),Q_ice_b,Q_b,Q_lith,f_grnd,rho_ice,L_ice) 
        
        ! Include internal melting in bmb_grnd 
        bmb_grnd = bmb_grnd - melt_internal 

        ! Finally, calculate the CTS height 
        H_cts = calc_cts_height(enth,T_ice,omega,T_pmp,cp,H_ice,zeta_aa)

        return 

    end subroutine calc_enth_column

    subroutine calc_enth_column_internal(enth,kappa,uz,advecxy,Q_strn,val_base,val_srf,thickness, &
                                            zeta_aa,zeta_ac,dzeta_a,dzeta_b,enth_ref,dt,k_cts,is_basal_flux)
        ! Thermodynamics solver for a given column of ice 
        ! Note zeta=height, k=1 base, k=nz surface 
        ! Note: nz = number of vertical boundaries (including zeta=0.0 and zeta=1.0), 
        ! temperature is defined for cell centers, plus a value at the surface and the base
        ! so nz_ac = nz_aa - 1 

        ! For notes on implicit form of advection terms, see eg http://farside.ph.utexas.edu/teaching/329/lectures/node90.html
        
        implicit none 

        real(wp), intent(INOUT) :: enth(:)        ! nz_aa [J kg] Ice column temperature
        real(wp), intent(IN)    :: kappa(:)       ! nz_aa [] Diffusivity
        real(wp), intent(IN)    :: uz(:)          ! nz_ac [m a-1] Vertical velocity 
        real(wp), intent(IN)    :: advecxy(:)     ! nz_aa [J kg a-1] Horizontal heat advection 
        real(wp), intent(IN)    :: Q_strn(:)      ! nz_aa [J kg a-1] Internal strain heat production in ice
        real(wp), intent(IN)    :: val_base       ! [J kg or flux] Basal boundary condition
        real(wp), intent(IN)    :: val_srf        ! [J kg] Surface temperature 
        real(wp), intent(IN)    :: thickness      ! [m] Ice thickness 
        real(wp), intent(IN)    :: zeta_aa(:)     ! nz_aa [--] Vertical sigma coordinates (zeta==height), layer centered aa-nodes
        real(wp), intent(IN)    :: zeta_ac(:)     ! nz_ac [--] Vertical height axis temperature (0:1), layer edges ac-nodes
        real(wp), intent(IN)    :: dzeta_a(:)     ! nz_aa [--] Solver discretization helper variable ak
        real(wp), intent(IN)    :: dzeta_b(:)     ! nz_aa [--] Solver discretization helper variable bk
        real(wp), intent(IN)    :: enth_ref       ! [J kg] Reference temperature to scale calculation
        real(wp), intent(IN)    :: dt             ! [a] Time step
        integer,  intent(IN)    :: k_cts          ! Index of the CTS (highest point at pressure melting point) 
        logical,  intent(IN)    :: is_basal_flux  ! Is basal condition flux condition (True) or Neumann (False)
        ! Local variables 
        integer  :: k, nz_aa
        real(wp) :: fac, fac_a, fac_b, uz_aa, dz, dzeta
        real(wp) :: h1, h2, afac_a, afac_b, afac_mid
        real(wp) :: kappa_a, kappa_b, dz1, dz2 
        real(wp), allocatable :: subd(:)      ! nz_aa 
        real(wp), allocatable :: diag(:)      ! nz_aa  
        real(wp), allocatable :: supd(:)      ! nz_aa 
        real(wp), allocatable :: rhs(:)       ! nz_aa 
        real(wp), allocatable :: solution(:)  ! nz_aa
        
        nz_aa = size(zeta_aa,1)

        allocate(subd(nz_aa))
        allocate(diag(nz_aa))
        allocate(supd(nz_aa))
        allocate(rhs(nz_aa))
        allocate(solution(nz_aa))

        ! == Ice base ==

        if (is_basal_flux) then 
            ! Impose basal flux (Neumann condition)

            ! Calculate dz for the bottom layer between the basal boundary
            ! (ac-node) and the centered (aa-node) temperature point above
            ! Note: zeta_aa(1) == zeta_ac(1) == bottom boundary 
            dz = thickness * (zeta_aa(2) - zeta_aa(1))

            ! backward Euler flux basal boundary condition
            subd(1) =  0.0_wp
            diag(1) =  1.0_wp
            supd(1) = -1.0_wp
            rhs(1)  = val_base * dz
                
        else if (k_cts .ge. 2) then 
            ! Layer above base is also temperate (with water likely present in the ice),
            ! set K0 dE/dz = 0. To do so, set basal enthalpy equal to enthalpy above

            subd(1) =  0.0_wp
            diag(1) =  1.0_wp
            supd(1) = -1.0_wp
            rhs(1)  =  0.0_wp
              
        else
            ! Impose basal enthalpy (Dirichlet condition) 

            subd(1) = 0.0_wp
            diag(1) = 1.0_wp
            supd(1) = 0.0_wp
            rhs(1)  = (val_base - enth_ref)

        end if 

        ! == Ice interior layers 2:nz_aa-1 ==

        do k = 2, nz_aa-1

            ! Get kappa for the lower and upper ac-nodes using harmonic mean from aa-nodes
            
            dz1 = zeta_ac(k-1)-zeta_aa(k-1)
            dz2 = zeta_aa(k)-zeta_ac(k-1)
            call calc_wtd_harmonic_mean(kappa_a,kappa(k-1),kappa(k),dz1,dz2)

            dz1 = zeta_ac(k)-zeta_aa(k)
            dz2 = zeta_aa(k+1)-zeta_ac(k)
            call calc_wtd_harmonic_mean(kappa_b,kappa(k),kappa(k+1),dz1,dz2)

            ! Special treatment of diffusivity at the CTS
            if (k .eq. k_cts+1) kappa_a = kappa(k-1)
            !if (k .eq. k_cts)   kappa_b = kappa(k+1) 
            
            ! Get diffusion factors
            fac_a   = -kappa_a*dzeta_a(k)*dt/thickness**2
            fac_b   = -kappa_b*dzeta_b(k)*dt/thickness**2


            ! Get implicit vertical advection term, ac => aa nodes
            uz_aa   = 0.5*(uz(k)+uz(k+1))

if (.TRUE.) then
    ! Use simple centered-difference scheme for vertical advection

            ! Vertical distance for centered difference advection scheme
            dzeta   = (zeta_aa(k+1)-zeta_aa(k-1))
            dz      = thickness * dzeta 

            subd(k) = fac_a - uz_aa * dt/dz
            supd(k) = fac_b + uz_aa * dt/dz
            diag(k) = 1.0_wp - fac_a - fac_b
            rhs(k)  = (enth(k)-enth_ref) - dt*advecxy(k) + dt*Q_strn(k)
            
else
    ! Use centered-difference advection scheme accounting for variable layer thickness

            ! Get grid spacing between neighboring points
            h1 = thickness * (zeta_aa(k)-zeta_aa(k-1))
            h2 = thickness * (zeta_aa(k+1)-zeta_aa(k))

            ! Get advective factor terms for centered difference with 
            ! uneven layers 
            afac_a   = -h2/(h1*(h1+h2))
            afac_mid = -(h1-h2)/(h1*h2)
            afac_b   = +h1/(h2*(h1+h2))

            subd(k) = fac_a + uz_aa * dt*afac_a
            supd(k) = fac_b + uz_aa * dt*afac_b
            diag(k) = 1.0_wp - fac_a - fac_b + uz_aa * dt*afac_mid
            rhs(k)  = (enth(k)-enth_ref) - dt*advecxy(k) + dt*Q_strn(k)
            
end if 

        end do 

        ! == Ice surface ==

        subd(nz_aa) = 0.0_wp
        diag(nz_aa) = 1.0_wp
        supd(nz_aa) = 0.0_wp
        rhs(nz_aa)  = (val_srf-enth_ref)

        ! == Call solver ==

        call solve_tridiag(subd,diag,supd,rhs,solution)

        ! Copy the solution into the temperature variable

        enth = solution + enth_ref 

        return 

    end subroutine calc_enth_column_internal

    ! ========== ENTHALPY ==========================================

    subroutine calc_enth_diffusivity(kappa,enth,enth_pmp,cp,kt,cr,rho_ice)
        ! Calculate the enthalpy vertical diffusivity for use with the diffusion solver:
        ! When water is present in the layer, set kappa=kappa_therm, else kappa=kappa_cold 

        implicit none 

        real(wp), intent(OUT) :: kappa(:)         ! [nz_aa]
        real(wp), intent(IN)  :: enth(:)          ! [nz_aa]
        real(wp), intent(IN)  :: enth_pmp(:)      ! [nz_aa]
        real(wp), intent(IN)  :: cp(:)
        real(wp), intent(IN)  :: kt(:)  
        real(wp), intent(IN)  :: cr 
        real(wp), intent(IN)  :: rho_ice
        
        ! Local variables
        integer   :: k, nz
        real(wp)  :: kappa_cold       ! Cold diffusivity 
        real(wp)  :: kappa_temp       ! Temperate diffusivity 
        
        nz = size(enth)

        kappa = 0.0 

        do k = 1, nz

            ! Determine kappa_cold and kappa_temp for this level 
            kappa_cold = kt(k) / (rho_ice*cp(k))
            kappa_temp = cr * kappa_cold 

            if (enth(k) .ge. enth_pmp(k)) then
                kappa(k) = kappa_temp 
            else 
                kappa(k) = kappa_cold 
            end if 

        end do

        return 

    end subroutine calc_enth_diffusivity
    
    function calc_cts_height(enth,T_ice,omega,T_pmp,cp,H_ice,zeta) result(H_cts)
        ! Calculate the height of the cold-temperate transition surface (m)
        ! within the ice sheet. 

        implicit none 

        real(wp), intent(IN) :: enth(:) 
        real(wp), intent(IN) :: T_ice(:) 
        real(wp), intent(IN) :: omega(:) 
        real(wp), intent(IN) :: T_pmp(:) 
        real(wp), intent(IN) :: cp(:)
        real(wp), intent(IN) :: H_ice  
        real(wp), intent(IN) :: zeta(:) 
        real(wp) :: H_cts 

        ! Local variables 
        integer  :: k, k_cts, nz 
        real(wp) :: f_lin, f_lin_0, dedz0, dedz1, zeta_cts 
        real(wp), allocatable :: enth_pmp(:) 

        integer :: i, n_iter, n_prime
        real(wp), allocatable :: zeta_prime(:) 
        real(wp), allocatable :: enth_prime(:) 
        
        nz = size(enth,1) 

        allocate(enth_pmp(nz))
        allocate(enth_prime(nz)) 

        ! Get enthalpy at the pressure melting point (no water content)
        enth_pmp = T_pmp * cp

        enth_prime = enth - enth_pmp 

        ! Determine height of CTS as highest temperate layer
        k_cts = get_cts_index(enth,enth_pmp)  

        if (k_cts .eq. 0) then 
            ! No temperate ice 
            H_cts = 0.0_wp 

        else if (k_cts .eq. nz) then 
            ! Whole column is temperate
            H_cts = H_ice

        else 

            ! Assume H_cts lies at center of last temperate cell (aa-node)
!             zeta_cts = zeta(k_cts)

!             ! Assume H_cts lies on ac-node between temperate and cold layers 
!             zeta_cts = 0.5_wp*(zeta(k_cts)+zeta(k_cts+1))

            ! Perform linear interpolation between enth(k_cts) and enth(k_cts+1) to find 
            ! where enth==enth_pmp.
            f_lin_0 = ( (enth(k_cts+1)-enth(k_cts)) - (enth_pmp(k_cts+1)-enth_pmp(k_cts)) )
            if (f_lin_0 .ne. 0.0) then 
                f_lin = (enth_pmp(k_cts)-enth(k_cts)) / f_lin_0
                if (f_lin .lt. 1e-2) f_lin = 0.0 
            else 
                f_lin = 1.0
            end if 

            zeta_cts = zeta(k_cts) + f_lin*(zeta(k_cts+1)-zeta(k_cts))
            
!             ec = (zc-z0)/(z1-z0)*(e1-e0) + e0 
!              0 = (zc-z0)/(z1-z0)*(e1-e0) + e0 
!            -e0 = (zc-z0)/(z1-z0)*(e1-e0)
!            -e0*(z1-z0)/(e1-e0) = zc-z0
           
!            zc = z0 - e0*(z1-z0)/(e1-e0)
            
!             if (abs(enth_prime(k_cts)-enth_prime(k_cts+1)) .lt. 1e-3) then 
!                 zeta_cts = zeta(k_cts+1)
!             else 
!                 zeta_cts = zeta(k_cts) - enth_prime(k_cts)*(zeta(k_cts+1)-zeta(k_cts))/(enth_prime(k_cts+1)-enth_prime(k_cts))
!             end if 

!             ! Further iterate to improve estimate of H_cts 
!             n_iter = 3
!             do i = 1, n_iter 

!             end do 
            
!             n_prime = 11 

!             allocate(zeta_prime(n_prime))
!             allocate(enth_prime(n_prime))
            
!             zeta_prime(1)       = zeta(k_cts-1)
!             zeta_prime(n_prime) = zeta(k_cts+1)

!             do i = 2, n_prime-1
!                 zeta_prime(i) = ((i-2)/(n_prime-3))*(zeta(k_cts+1)-zeta(k_cts)) + zeta(k_cts)
!             end do 

!             enth_prime = interp_spline(zeta,enth-enth_pmp,zeta_prime)

!             i = minloc(abs(enth_prime),1)
!             zeta_cts = zeta_prime(i) 

! !             i = maxloc(abs(enth_prime),1,mask=enth_prime .lt. 0.0_wp)
! !             f_lin = (zeta_prime(i+1)-zeta_prime(i)) / (enth_prime(i+1)-enth_prime(i))
! !             if (abs(f_lin) .lt. 1e-3) f_lin = 0.0_wp 
! !             zeta_cts = (1.0_wp-f_lin)*zeta_prime(i) 
    
            H_cts    = H_ice*zeta_cts

        end if 

        return 


    end function calc_cts_height

    subroutine calc_dzeta_terms(dzeta_a,dzeta_b,zeta_aa,zeta_ac)
        ! zeta_aa  = depth axis at layer centers (plus base and surface values)
        ! zeta_ac  = depth axis (1: base, nz: surface), at layer boundaries
        ! Calculate ak, bk terms as defined in Hoffmann et al (2018)
        implicit none 

        real(wp), allocatable, intent(INOUT) :: dzeta_a(:)    ! nz_aa
        real(wp), allocatable, intent(INOUT) :: dzeta_b(:)    ! nz_aa
        real(wp), intent(IN) :: zeta_aa(:)    ! nz_aa 
        real(wp), intent(IN) :: zeta_ac(:)    ! nz_ac == nz_aa+1 

        ! Local variables 
        integer :: k, nz_layers, nz_aa    

        nz_aa = size(zeta_aa)

        ! Allocate dzeta_a and dzeta_b to match size of input zeta_aa vector
        if (allocated(dzeta_a)) deallocate(dzeta_a)
        if (allocated(dzeta_b)) deallocate(dzeta_b)
        allocate(dzeta_a(nz_aa))
        allocate(dzeta_b(nz_aa))
        
        ! Note: zeta_aa is calculated outside in the main program 

        ! Initialize dzeta_a/dzeta_b to zero, first and last indices will not be used (end points)
        dzeta_a = 0.0 
        dzeta_b = 0.0 
        
        do k = 2, nz_aa-1 
            dzeta_a(k) = 1.0/ ( (zeta_ac(k+1) - zeta_ac(k)) * (zeta_aa(k) - zeta_aa(k-1)) )
        enddo

        do k = 2, nz_aa-1
            dzeta_b(k) = 1.0/ ( (zeta_ac(k+1) - zeta_ac(k)) * (zeta_aa(k+1) - zeta_aa(k)) )
        end do

        return 

    end subroutine calc_dzeta_terms

    subroutine calc_wtd_harmonic_mean(var_ave,var1,var2,wt1,wt2)

        implicit none 

        real(wp), intent(OUT) :: var_ave 
        real(wp), intent(IN)  :: var1 
        real(wp), intent(IN)  :: var2 
        real(wp), intent(IN)  :: wt1 
        real(wp), intent(IN)  :: wt2 
        
        ! Local variables 
        real(wp), parameter   :: tol = 1e-5 

        !var_ave = ( wt1*(var1+tol)**(-1.0) + wt2*(var2+tol)**(-1.0) )**(-1.0)
        var_ave = ( (wt1*(var1)**(-1.0) + wt2*(var2)**(-1.0)) / (wt1+wt2) )**(-1.0)

        return 

    end subroutine calc_wtd_harmonic_mean

    subroutine calc_zeta_twolayers(zeta_pt,zeta_pc,zeta_scale,zeta_exp)
        ! Calculate the vertical layer-edge axis (vertical ac-nodes)
        ! and the vertical cell-center axis (vertical aa-nodes),
        ! including an extra zero-thickness aa-node at the base and surface

        ! This is built in two-steps, first for the basal temperate layer
        ! and second for the overlying cold layer. The height of the border
        ! is the CTS height, which will be defined for each column. The temperate layer is populated with an 
        ! evenly-spaced (linear) axis up to upper boundary, while the cold layer follows the 
        ! parameter options zeta_scale and zeta_exp. 

        implicit none 

        real(wp),   intent(INOUT) :: zeta_pt(:)
        real(wp),   intent(INOUT) :: zeta_pc(:) 
        character(*), intent(IN)  :: zeta_scale 
        real(wp),   intent(IN)    :: zeta_exp 

        ! Local variables
        integer :: k, nz_pt, nz_pc 

        integer :: nz_ac 
        real(wp), allocatable :: zeta_ac(:) 

        nz_pt  = size(zeta_pt)
        nz_pc  = size(zeta_pc) 

        ! ===== Temperate layer ===================================

        nz_ac = nz_pt - 1
        allocate(zeta_ac(nz_ac))

        ! Linear scale for cell boundaries
        do k = 1, nz_ac
            zeta_ac(k) = 0.0 + 1.0*(k-1)/real(nz_ac-1)
        end do 

        ! Get zeta_aa (between zeta_ac values, as well as at the base and surface)
        zeta_pt(1) = 0.0 
        do k = 2, nz_pt-1
            zeta_pt(k) = 0.5 * (zeta_ac(k-1)+zeta_ac(k))
        end do 
        zeta_pt(nz_pt) = 1.0 

        ! ===== Cold layer ========================================

        nz_ac = nz_pc - 1
        deallocate(zeta_ac)
        allocate(zeta_ac(nz_ac))

        ! Linear scale for cell boundaries
        do k = 1, nz_ac
            zeta_ac(k) = 0.0 + 1.0*(k-1)/real(nz_ac-1)
        end do 

        ! Scale zeta to produce different resolution through column if desired
        ! zeta_scale = ["linear","exp","wave"]
        select case(trim(zeta_scale))
            
            case("exp")
                ! Increase resolution at the base 
                zeta_ac = zeta_ac**(zeta_exp) 

            case("tanh")
                ! Increase resolution at base and surface 

                zeta_ac = tanh(1.0*pi*(zeta_ac-0.5))
                zeta_ac = zeta_ac - minval(zeta_ac)
                zeta_ac = zeta_ac / maxval(zeta_ac)

            case DEFAULT
            ! Do nothing, scale should be linear as defined above
        
        end select  
        
        ! Get zeta_aa (between zeta_ac values, as well as at the base and surface)
        zeta_pc(1) = 0.0 
        do k = 2, nz_pc-1
            zeta_pc(k) = 0.5 * (zeta_ac(k-1)+zeta_ac(k))
        end do 
        zeta_pc(nz_pc) = 1.0 

        return 

    end subroutine calc_zeta_twolayers
    
    subroutine calc_zeta_combined(zeta_aa,zeta_ac,zeta_pt,zeta_pc,H_cts,H_ice)
        ! Take two-layer axis and combine into one axis based on relative CTS height
        ! f_cts = H_cts / H_ice 

        implicit none 

        real(wp), intent(INOUT) :: zeta_aa(:) 
        real(wp), intent(INOUT) :: zeta_ac(:) 
        real(wp), intent(IN)    :: zeta_pt(:) 
        real(wp), intent(IN)    :: zeta_pc(:) 
        real(wp), intent(IN)    :: H_cts 
        real(wp), intent(IN)    :: H_ice 

        ! Local variables 
        integer  :: k 
        integer  :: nzt, nztc, nzc, nz_aa, nz_ac  
        real(wp) :: f_cts

        nz_aa = size(zeta_aa,1)
        nz_ac = size(zeta_ac,1)  ! == nz_aa-1
        nzt   = size(zeta_pt,1)
        nzc   = size(zeta_pc,1) 

        if (nzt+(nzc-1)  .ne. nz_aa) then 
            write(*,*) "calc_zeta_combined:: Error: Two-layer axis length does not match combined axis length."
            write(*,*) "nzt, nzc-1, nz_aa: ", nzt, nzc-1, nz_aa 
            stop 
        end if 

        ! Get f_cts 
        if (H_ice .gt. 0.0) then 
            f_cts = max(H_cts / H_ice,0.01)
        else 
            f_cts = 0.01 
        end if 

        zeta_aa(1:nzt) = zeta_pt(1:nzt)*f_cts
        zeta_aa(nzt+1:nzt+nzc) = f_cts + zeta_pc(2:nzc)*(1.0-f_cts)

        ! Get zeta_ac again (boundaries between zeta_aa values, as well as at the base and surface)
        zeta_ac(1) = 0.0_wp 
        do k = 2, nz_ac-1
            zeta_ac(k) = 0.5_wp * (zeta_aa(k)+zeta_aa(k+1))
        end do 
        zeta_ac(nz_ac) = 1.0_wp 

        return 

    end subroutine calc_zeta_combined

    function get_cts_index(enth,enth_pmp) result(k_cts)

        implicit none 

        real(wp), intent(IN) :: enth(:) 
        real(wp), intent(IN) :: enth_pmp(:) 
        integer :: k_cts 

        ! Local variables 
        integer :: k, nz 

        nz = size(enth,1) 

        k_cts = 1 
        do k = 1, nz 
            if (enth(k) .ge. enth_pmp(k)) then 
                k_cts = k 
            else 
                exit 
            end if 
        end do 
            
        return 

    end function get_cts_index

    function interp_linear_point(x0,x1,y0,y1,xout) result(yout)
        ! Interpolates for the y value at the desired x value, 
        ! given x and y values around the desired point.
        ! Solution outside of range x0 < x < x1 bounded by y0 < y < y1 

        implicit none

        real(wp), intent(IN)  :: x0,x1,y0,y1, xout
        real(wp) :: yout
        real(wp) :: alph

        if (xout .le. x0) then 
            yout = y0 
        else if (xout .ge. x1) then 
            yout = y1 
        else 
            alph = (xout - x0) / (x1 - x0)
            yout = y0 + alph*(y1 - y0)
        end if 

        return

    end function interp_linear_point

end module ice_enthalpy


