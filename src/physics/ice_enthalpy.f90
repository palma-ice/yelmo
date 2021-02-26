module ice_enthalpy 
    ! Module contains the ice temperature and basal mass balance (grounded) solution

    use yelmo_defs, only : prec, wp, pi, g, sec_year, rho_ice, rho_sw, rho_w, rho_l, L_ice  
    use solver_tridiagonal, only : solve_tridiag 
    use thermodynamics, only : calc_bmb_grounded, calc_bmb_grounded_enth, calc_advec_vertical_column

    !use interp1D 

    implicit none
    
    private
    public :: calc_temp_column
    public :: calc_temp_column_bedrock
    public :: calc_enth_column0
    public :: calc_enth_column
    public :: convert_to_enthalpy
    public :: convert_from_enthalpy_column
    public :: calc_dzeta_terms
    public :: calc_zeta_twolayers
    public :: calc_zeta_combined
    public :: get_cts_index

contains 
    
    subroutine calc_temp_column(enth,T_ice,omega,bmb_grnd,Q_ice_b,H_cts,T_pmp,cp,kt,advecxy,uz, &
                                Q_strn,Q_b,Q_lith,T_srf,T_shlf,H_ice,H_w,f_grnd,zeta_aa,zeta_ac, &
                                dzeta_a,dzeta_b,omega_max,T0,dt)
        ! Thermodynamics solver for a given column of ice 
        ! Note zeta=height, k=1 base, k=nz surface 
        ! Note: nz = number of vertical boundaries (including zeta=0.0 and zeta=1.0), 
        ! temperature is defined for cell centers, plus a value at the surface and the base
        ! so nz_ac = nz_aa - 1 

        ! For notes on implicit form of advection terms, see eg http://farside.ph.utexas.edu/teaching/329/lectures/node90.html
        
        implicit none 

        real(prec), intent(INOUT) :: enth(:)        ! nz_aa [J kg] Ice column enthalpy
        real(prec), intent(INOUT) :: T_ice(:)       ! nz_aa [K] Ice column temperature
        real(prec), intent(INOUT) :: omega(:)       ! nz_aa [-] Ice column water content fraction
        real(prec), intent(INOUT) :: bmb_grnd       ! [m a-1] Basal mass balance (melting is negative)
        real(prec), intent(OUT)   :: Q_ice_b        ! [mW m-2] Ice basal heat flux (positive up)
        real(prec), intent(OUT)   :: H_cts          ! [m] cold-temperate transition surface (CTS) height
        real(prec), intent(IN)    :: T_pmp(:)       ! nz_aa [K] Pressure melting point temp.
        real(prec), intent(IN)    :: cp(:)          ! nz_aa [J kg-1 K-1] Specific heat capacity
        real(prec), intent(IN)    :: kt(:)          ! nz_aa [J a-1 m-1 K-1] Heat conductivity 
        real(prec), intent(IN)    :: advecxy(:)     ! nz_aa [K a-1] Horizontal heat advection 
        real(prec), intent(IN)    :: uz(:)          ! nz_ac [m a-1] Vertical velocity 
        real(prec), intent(IN)    :: Q_strn(:)      ! nz_aa [J a-1 m-3] Internal strain heat production in ice
        real(prec), intent(IN)    :: Q_b            ! [mW m-2] Basal frictional heat production
        real(prec), intent(IN)    :: Q_lith         ! [mW m-2] Bedrock heat flux (positive up)
        real(prec), intent(IN)    :: T_srf          ! [K] Surface temperature 
        real(prec), intent(IN)    :: T_shlf         ! [K] Marine-shelf interface temperature
        real(prec), intent(IN)    :: H_ice          ! [m] Ice thickness 
        real(prec), intent(IN)    :: H_w            ! [m] Basal water layer thickness 
        real(prec), intent(IN)    :: f_grnd         ! [--] Grounded fraction
        real(prec), intent(IN)    :: zeta_aa(:)     ! nz_aa [--] Vertical sigma coordinates (zeta==height), layer centered aa-nodes
        real(prec), intent(IN)    :: zeta_ac(:)     ! nz_ac [--] Vertical height axis temperature (0:1), layer edges ac-nodes
        real(prec), intent(IN)    :: dzeta_a(:)     ! nz_aa [--] Solver discretization helper variable ak
        real(prec), intent(IN)    :: dzeta_b(:)     ! nz_aa [--] Solver discretization helper variable bk
        real(prec), intent(IN)    :: omega_max      ! [-] Maximum allowed water fraction inside ice, typically omega_max=0.02 
        real(prec), intent(IN)    :: T0             ! [K or degreesCelcius] Reference melting temperature  
        real(prec), intent(IN)    :: dt             ! [a] Time step  

        ! Local variables 
        integer    :: k, nz_aa, nz_ac
        real(prec) :: H_w_predicted
        real(prec) :: dz, dz1, dz2, d2Tdz2 
        real(prec) :: T00, T01, T02, zeta_now  
        real(prec) :: T_excess
        real(prec) :: Q_ice_b_now, Q_b_now, Q_lith_now 
        real(prec) :: melt_internal
        real(prec) :: val_base, val_srf 
        logical    :: is_basal_flux  
        logical    :: is_surf_flux 

        real(prec), allocatable :: kappa_aa(:)    ! aa-nodes
        real(prec), allocatable :: Q_strn_now(:)  ! aa-nodes
        
        real(prec), parameter   :: T_ref = 273.15_prec   


        nz_aa = size(zeta_aa,1)

        allocate(kappa_aa(nz_aa))
        allocate(Q_strn_now(nz_aa))

        ! Get lithosphere heat flux in proper units 
        Q_lith_now = Q_lith*1e-3*sec_year   ! [mW m-2] => [J m-2 a-1]

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

            if (H_w_predicted .gt. 0.0_prec) then 
                ! Temperate at bed 
                ! Hold basal temperature at pressure melting point

                val_base      = T_pmp(1)
                is_basal_flux = .FALSE.

            else if ( T_ice(1) .lt. T_pmp(1) .or. H_w_predicted .lt. 0.0_prec ) then
                ! Frozen at bed, or about to become frozen 

                ! backward Euler flux basal boundary condition
                val_base      = -(Q_b_now + Q_lith_now) / kt(1)
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

        ! Also set omega to constant value where ice is temperate just for some consistency 
        omega = 0.0 
!             where (T_ice .ge. T_pmp) omega = omega_max 

        ! Finally, get enthalpy too 
        call convert_to_enthalpy(enth,T_ice,omega,T_pmp,cp,L_ice)

        ! Calculate heat flux at ice base as temperature gradient * conductivity [J a-1 m-2]
        if (H_ice .gt. 0.0_prec) then 

            ! 1st order, upwind gradient dTdz 
            ! Works, but can cause oscillations in H_w 
            dz = H_ice * (zeta_aa(2)-zeta_aa(1))
            Q_ice_b_now = -kt(1) * (T_ice(2) - T_ice(1)) / dz 

        else 

            Q_ice_b_now = 0.0  

        end if 
        
        ! Calculate Q_ice_b for global output 
        Q_ice_b = Q_ice_b_now*1e3/sec_year      ! [J a-1 m-2] => [mW m-2]

        ! Calculate basal mass balance (valid for grounded ice only)
        call calc_bmb_grounded(bmb_grnd,T_ice(1)-T_pmp(1),Q_ice_b_now,Q_b_now,Q_lith_now,f_grnd,rho_ice)
        
        ! Include internal melting in bmb_grnd 
        bmb_grnd = bmb_grnd - melt_internal 


        ! Finally, calculate the CTS height 
        H_cts = calc_cts_height(enth,T_ice,omega,T_pmp,cp,H_ice,zeta_aa)

        return 

    end subroutine calc_temp_column

    subroutine calc_temp_column_bedrock(enth,temp,Q_lith,cp,kt,Q_ice_b,Q_geo,T_srf,H_lith, &
                                                zeta_aa,zeta_ac,dzeta_a,dzeta_b,dt)
        ! Thermodynamics solver for a given column of ice 
        ! Note zeta=height, k=1 base, k=nz surface 
        ! Note: nz = number of vertical boundaries (including zeta=0.0 and zeta=1.0), 
        ! temperature is defined for cell centers, plus a value at the surface and the base
        ! so nz_ac = nz_aa - 1 

        ! For notes on implicit form of advection terms, see eg http://farside.ph.utexas.edu/teaching/329/lectures/node90.html
        
        implicit none 

        real(prec), intent(INOUT) :: enth(:)        ! nz_aa [J kg] Ice column enthalpy
        real(prec), intent(INOUT) :: temp(:)        ! nz_aa [K] Ice column temperature
        real(prec), intent(OUT)   :: Q_lith         ! [mW m-2] Bed surface heat flux (positive up)
        real(prec), intent(IN)    :: cp(:)          ! nz_aa [J kg-1 K-1] Specific heat capacity
        real(prec), intent(IN)    :: kt(:)          ! nz_aa [J a-1 m-1 K-1] Heat conductivity 
        real(prec), intent(IN)    :: Q_ice_b        ! [mW m-2] Ice basal heat flux (positive up)
        real(prec), intent(IN)    :: Q_geo          ! [mW m-2] Bedrock heat flux (positive up)
        real(prec), intent(IN)    :: T_srf          ! [K] Surface temperature 
        real(prec), intent(IN)    :: H_lith         ! [m] Lithosphere thickness 
        real(prec), intent(IN)    :: zeta_aa(:)     ! nz_aa [--] Vertical sigma coordinates (zeta==height), layer centered aa-nodes
        real(prec), intent(IN)    :: zeta_ac(:)     ! nz_ac [--] Vertical height axis temperature (0:1), layer edges ac-nodes
        real(prec), intent(IN)    :: dzeta_a(:)     ! nz_aa [--] Solver discretization helper variable ak
        real(prec), intent(IN)    :: dzeta_b(:)     ! nz_aa [--] Solver discretization helper variable bk
        real(prec), intent(IN)    :: dt             ! [a] Time step 

        ! Local variables 
        integer    :: k, nz_aa, nz_ac
        real(prec) :: Q_lith_now 
        real(prec) :: Q_ice_b_now 
        real(prec) :: Q_geo_now 
        real(prec) :: dz 
        real(prec) :: val_base, val_srf 
        logical    :: is_basal_flux  
        logical    :: is_surf_flux  

        real(prec), allocatable :: kappa_aa(:)      ! aa-nodes
        
        real(prec), allocatable, target :: zeros(:)
        real(prec), pointer     :: advecxy(:)
        real(prec), pointer     :: Q_strn(:)
        real(prec), pointer     :: uz(:)
        
        real(prec), parameter   :: T_ref = 273.15_prec   

        nz_aa = size(zeta_aa,1)
        nz_ac = size(zeta_ac,1) 

        ! Get basal heat flux from ice in proper units 
        Q_ice_b_now = Q_ice_b*1e-3*sec_year     ! [mW m-2] => [J m-2 a-1]

        ! Get deep lithosphere heat flux in proper units 
        Q_geo_now = Q_geo*1e-3*sec_year   ! [mW m-2] => [J m-2 a-1]

        ! Step 0: Calculate diffusivity on cell centers (aa-nodes)

        kappa_aa = kt / (rho_l*cp)
        
        ! Set unused terms to zero (via pointers to save allocations)
        allocate(zeros(nz_aa))
        zeros = 0.0_wp 
        
        advecxy => zeros(1:nz_aa) 
        Q_strn  => zeros(1:nz_aa) 
        uz      => zeros(1:nz_ac) 

        ! === Surface boundary condition =====================

        ! if ( T_srf .lt. T_pmp_srf ) then
        !         ! Frozen at bed

        !     ! backward Euler flux surface boundary condition
        !     val_srf      =  -Q_ice_b_now / kt(nz_aa)
        !     is_surf_flux = .TRUE.   

        ! else 
            ! Temperate at bed surface
            ! Hold bed surface temperature at ice base temperature

            val_srf      = T_srf
            is_surf_flux = .FALSE. 
            
        ! end if
        
        ! === Basal boundary condition =====================

        ! backward Euler flux basal boundary condition
        val_base = -Q_geo_now / kt(1)
        is_basal_flux = .TRUE. 
        
        ! === Solver =============================
     
        call calc_temp_column_internal(temp,kappa_aa,uz,advecxy,Q_strn,val_base,val_srf,H_lith, &
                                                zeta_aa,zeta_ac,dzeta_a,dzeta_b,T_ref,dt, &
                                                is_basal_flux,is_surf_flux)

        ! Convert to enthalpy too 
        call convert_to_enthalpy(enth,temp,omega=0.0_wp,T_pmp=0.0_wp,cp=cp,L=0.0_wp)

        ! Calculate heat flux at lithosphere surface as temperature gradient * conductivity [J a-1 m-2]

        ! 1st order, upwind gradient dTdz 
        ! Works, but can cause oscillations in H_w 
        dz = H_lith * (zeta_aa(nz_aa)-zeta_aa(nz_aa-1))
        Q_lith_now = -kt(nz_aa) * (temp(nz_aa) - temp(nz_aa-1)) / dz 

        ! Set global units to [mW m-2]
        Q_lith = Q_lith_now *1e3 / sec_year 
        
        return 

    end subroutine calc_temp_column_bedrock

    subroutine calc_temp_column_internal(temp,kappa,uz,advecxy,Q_strn,val_base,val_srf,thickness, &
                                                zeta_aa,zeta_ac,dzeta_a,dzeta_b,T_ref,dt,is_basal_flux,is_surf_flux)
        ! Thermodynamics solver 1D (eg, for column of ice or bedrock)
        ! Note zeta= column height, k=1 base, k=nz surface 
        ! Note: nz = number of vertical boundaries (including zeta=0.0 and zeta=1.0), 
        ! temperature is defined for cell centers, plus a value at the surface and the base
        ! so nz_ac = nz_aa - 1 

        ! For notes on implicit form of advection terms, see eg http://farside.ph.utexas.edu/teaching/329/lectures/node90.html
        
        implicit none 

        real(prec), intent(INOUT) :: temp(:)        ! nz_aa [K] Ice column temperature
        real(prec), intent(IN)    :: kappa(:)       ! nz_aa [] Diffusivity
        real(prec), intent(IN)    :: uz(:)          ! nz_ac [m a-1] Vertical velocity 
        real(prec), intent(IN)    :: advecxy(:)     ! nz_aa [K a-1] Horizontal heat advection 
        real(prec), intent(IN)    :: Q_strn(:)      ! nz_aa [K a-1] Internal strain heat production in ice
        real(prec), intent(IN)    :: val_base       ! [K or flux] Basal boundary condition
        real(prec), intent(IN)    :: val_srf        ! [K or flux] Surface temperature 
        real(prec), intent(IN)    :: thickness      ! [m] Total column thickness 
        real(prec), intent(IN)    :: zeta_aa(:)     ! nz_aa [--] Vertical sigma coordinates (zeta==height), layer centered aa-nodes
        real(prec), intent(IN)    :: zeta_ac(:)     ! nz_ac [--] Vertical height axis temperature (0:1), layer edges ac-nodes
        real(prec), intent(IN)    :: dzeta_a(:)     ! nz_aa [--] Solver discretization helper variable ak
        real(prec), intent(IN)    :: dzeta_b(:)     ! nz_aa [--] Solver discretization helper variable bk
        real(prec), intent(IN)    :: T_ref          ! [K] Reference temperature to scale calculation
        real(prec), intent(IN)    :: dt             ! [a] Time step 
        logical,    intent(IN)    :: is_basal_flux  ! Is basal bc flux condition (True) or Neumann (False)
        logical,    intent(IN)    :: is_surf_flux   ! Is surf  bc flux condition (True) or Neumann (False)
        
        ! Local variables 
        integer    :: k, nz_aa
        real(prec) :: fac, fac_a, fac_b, uz_aa, dz
        real(prec) :: kappa_a, kappa_b, dz1, dz2 
        real(prec), allocatable :: subd(:)      ! nz_aa 
        real(prec), allocatable :: diag(:)      ! nz_aa  
        real(prec), allocatable :: supd(:)      ! nz_aa 
        real(prec), allocatable :: rhs(:)       ! nz_aa 
        real(prec), allocatable :: solution(:)  ! nz_aa
        
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
            subd(1) =  0.0_prec
            diag(1) = -1.0_prec
            supd(1) =  1.0_prec
            rhs(1)  = val_base * dz
                
        else 
            ! Impose basal temperature (Dirichlet condition)  

            subd(1) = 0.0_prec
            diag(1) = 1.0_prec
            supd(1) = 0.0_prec
            rhs(1)  = (val_base - T_ref)

        end if 

        ! == Ice interior layers 2:nz_aa-1 ==

        do k = 2, nz_aa-1

            ! Get implicit vertical advection term, ac => aa nodes
            uz_aa   = 0.5*(uz(k-1)+uz(k))

            ! Get kappa for the lower and upper ac-nodes using harmonic mean from aa-nodes
            
            dz1 = zeta_ac(k-1)-zeta_aa(k-1)
            dz2 = zeta_aa(k)-zeta_ac(k-1)
            call calc_wtd_harmonic_mean(kappa_a,kappa(k-1),kappa(k),dz1,dz2)

            dz1 = zeta_ac(k)-zeta_aa(k)
            dz2 = zeta_aa(k+1)-zeta_ac(k)
            call calc_wtd_harmonic_mean(kappa_b,kappa(k),kappa(k+1),dz1,dz2)

            ! Vertical distance for centered difference advection scheme
            dz      =  thickness*(zeta_aa(k+1)-zeta_aa(k-1))
            
            fac_a   = -kappa_a*dzeta_a(k)*dt/thickness**2
            fac_b   = -kappa_b*dzeta_b(k)*dt/thickness**2

            subd(k) = fac_a - uz_aa * dt/dz
            supd(k) = fac_b + uz_aa * dt/dz
            diag(k) = 1.0_prec - fac_a - fac_b
            rhs(k)  = (temp(k)-T_ref) - dt*advecxy(k) + dt*Q_strn(k)
            
        end do 

        ! == Column surface ==

        if (is_surf_flux) then 
            ! Impose surface flux (Neumann condition)

            ! Calculate dz for the top layer between the surface boundary
            ! (ac-node) and the centered (aa-node) temperature point above
            ! Note: zeta_aa(nz_aa) == zeta_ac(nz_ac) == top boundary 
            dz = thickness * (zeta_aa(nz_aa) - zeta_aa(nz_aa-1))
            
            ! backward Euler flux surface boundary condition
            subd(nz_aa) = -1.0_prec
            diag(nz_aa) =  1.0_prec
            supd(nz_aa) =  0.0_prec
            rhs(nz_aa)  = val_srf * dz
            
        else 
            ! Impose surface temperature (Dirichlet condition) 

            subd(nz_aa) = 0.0_prec
            diag(nz_aa) = 1.0_prec
            supd(nz_aa) = 0.0_prec
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
                                dzeta_a,dzeta_b,cr,omega_max,T0,dt)
        ! Thermodynamics solver for a given column of ice 
        ! Note zeta=height, k=1 base, k=nz surface 
        ! Note: nz = number of vertical boundaries (including zeta=0.0 and zeta=1.0), 
        ! temperature is defined for cell centers, plus a value at the surface and the base
        ! so nz_ac = nz_aa - 1 

        ! For notes on implicit form of advection terms, see eg http://farside.ph.utexas.edu/teaching/329/lectures/node90.html
        
        implicit none 

        real(prec), intent(INOUT) :: enth(:)        ! nz_aa [J kg] Ice column enthalpy
        real(prec), intent(INOUT) :: T_ice(:)       ! nz_aa [K] Ice column temperature
        real(prec), intent(INOUT) :: omega(:)       ! nz_aa [-] Ice column water content fraction
        real(prec), intent(INOUT) :: bmb_grnd       ! [m a-1] Basal mass balance (melting is negative)
        real(prec), intent(OUT)   :: Q_ice_b        ! [J a-1 m-2] Ice basal heat flux (positive up)
        real(prec), intent(OUT)   :: H_cts          ! [m] cold-temperate transition surface (CTS) height
        real(prec), intent(IN)    :: T_pmp(:)       ! nz_aa [K] Pressure melting point temp.
        real(prec), intent(IN)    :: cp(:)          ! nz_aa [J kg-1 K-1] Specific heat capacity
        real(prec), intent(IN)    :: kt(:)          ! nz_aa [J a-1 m-1 K-1] Heat conductivity 
        real(prec), intent(IN)    :: advecxy(:)     ! nz_aa [K a-1] Horizontal heat advection 
        real(prec), intent(IN)    :: uz(:)          ! nz_ac [m a-1] Vertical velocity 
        real(prec), intent(IN)    :: Q_strn(:)      ! nz_aa [J a-1 m-3] Internal strain heat production in ice
        real(prec), intent(IN)    :: Q_b            ! [J a-1 m-2] Basal frictional heat production
        real(prec), intent(IN)    :: Q_lith         ! [J a-1 m-2] Bedrock heat flux (positive up)
        real(prec), intent(IN)    :: T_srf          ! [K] Surface temperature 
        real(prec), intent(IN)    :: T_shlf         ! [K] Marine-shelf interface temperature
        real(prec), intent(IN)    :: H_ice          ! [m] Ice thickness 
        real(prec), intent(IN)    :: H_w            ! [m] Basal water layer thickness 
        real(prec), intent(IN)    :: f_grnd         ! [--] Grounded fraction
        real(prec), intent(IN)    :: zeta_aa(:)     ! nz_aa [--] Vertical sigma coordinates (zeta==height), layer centered aa-nodes
        real(prec), intent(IN)    :: zeta_ac(:)     ! nz_ac [--] Vertical height axis temperature (0:1), layer edges ac-nodes
        real(prec), intent(IN)    :: dzeta_a(:)     ! nz_aa [--] Solver discretization helper variable ak
        real(prec), intent(IN)    :: dzeta_b(:)     ! nz_aa [--] Solver discretization helper variable bk
        real(prec), intent(IN)    :: cr             ! [--] Conductivity ratio (kappa_water / kappa_ice)
        real(prec), intent(IN)    :: omega_max      ! [-] Maximum allowed water fraction inside ice, typically omega_max=0.02 
        real(prec), intent(IN)    :: T0             ! [K or degreesCelcius] Reference melting temperature  
        real(prec), intent(IN)    :: dt             ! [a] Time step 

        ! Local variables 
        integer    :: k, nz_aa, nz_ac
        integer    :: k_cts  
        real(prec) :: H_w_predicted
        real(prec) :: dz 
        real(prec) :: omega_excess
        real(prec) :: melt_internal
        real(prec) :: val_base, val_srf 
        logical    :: is_basal_flux  

        real(prec), allocatable :: kappa_aa(:)    ! aa-nodes
        real(prec), allocatable :: Q_strn_now(:)  ! aa-nodes
        real(prec), allocatable :: enth_pmp(:)    ! aa-nodes
        
        real(prec), parameter   :: enth_ref = 273.15_prec * 2009.0_prec    ! [K] * [J kg-1 K-1]  


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

            if (H_w_predicted .gt. 0.0_prec) then 

                val_base = enth_pmp(1)
                is_basal_flux = .FALSE.

            else if ( enth(1) .lt. enth_pmp(1) .or. H_w_predicted .lt. 0.0_prec ) then
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
        if (H_ice .gt. 0.0_prec) then 
            dz = H_ice * (zeta_aa(2)-zeta_aa(1))
!             Q_ice_b = kappa_aa(1) * rho_ice * (enth(2) - enth(1)) / dz    ! <== Problematic
            Q_ice_b = kt(1) * ( T_ice(2) - T_ice(1) ) / dz
        else
            Q_ice_b = 0.0_prec 
        end if 

        ! Calculate basal mass balance 
        call calc_bmb_grounded(bmb_grnd,T_ice(1)-T_pmp(1),Q_ice_b,Q_b,Q_lith,f_grnd,rho_ice) 
!         call calc_bmb_grounded_enth(bmb_grnd,T_ice(1)-T_pmp(1),omega(1),Q_ice_b,Q_b,Q_lith,f_grnd,rho_ice) 
        
        ! Include internal melting in bmb_grnd 
        bmb_grnd = bmb_grnd - melt_internal 

        ! Finally, calculate the CTS height 
        H_cts = calc_cts_height(enth,T_ice,omega,T_pmp,cp,H_ice,zeta_aa)

        return 

    end subroutine calc_enth_column

    subroutine calc_enth_column_internal(enth,kappa,uz,advecxy,Q_strn,val_base,val_srf,H_ice, &
                                            zeta_aa,zeta_ac,dzeta_a,dzeta_b,enth_ref,dt,k_cts,is_basal_flux)
        ! Thermodynamics solver for a given column of ice 
        ! Note zeta=height, k=1 base, k=nz surface 
        ! Note: nz = number of vertical boundaries (including zeta=0.0 and zeta=1.0), 
        ! temperature is defined for cell centers, plus a value at the surface and the base
        ! so nz_ac = nz_aa - 1 

        ! For notes on implicit form of advection terms, see eg http://farside.ph.utexas.edu/teaching/329/lectures/node90.html
        
        implicit none 

        real(prec), intent(INOUT) :: enth(:)        ! nz_aa [J kg] Ice column temperature
        real(prec), intent(IN)    :: kappa(:)       ! nz_aa [] Diffusivity
        real(prec), intent(IN)    :: uz(:)          ! nz_ac [m a-1] Vertical velocity 
        real(prec), intent(IN)    :: advecxy(:)     ! nz_aa [J kg a-1] Horizontal heat advection 
        real(prec), intent(IN)    :: Q_strn(:)      ! nz_aa [J kg a-1] Internal strain heat production in ice
        real(prec), intent(IN)    :: val_base       ! [J kg or flux] Basal boundary condition
        real(prec), intent(IN)    :: val_srf        ! [J kg] Surface temperature 
        real(prec), intent(IN)    :: H_ice          ! [m] Ice thickness 
        real(prec), intent(IN)    :: zeta_aa(:)     ! nz_aa [--] Vertical sigma coordinates (zeta==height), layer centered aa-nodes
        real(prec), intent(IN)    :: zeta_ac(:)     ! nz_ac [--] Vertical height axis temperature (0:1), layer edges ac-nodes
        real(prec), intent(IN)    :: dzeta_a(:)     ! nz_aa [--] Solver discretization helper variable ak
        real(prec), intent(IN)    :: dzeta_b(:)     ! nz_aa [--] Solver discretization helper variable bk
        real(prec), intent(IN)    :: enth_ref       ! [J kg] Reference temperature to scale calculation
        real(prec), intent(IN)    :: dt             ! [a] Time step
        integer,    intent(IN)    :: k_cts          ! Index of the CTS (highest point at pressure melting point) 
        logical,    intent(IN)    :: is_basal_flux  ! Is basal condition flux condition (True) or Neumann (False)
        ! Local variables 
        integer    :: k, nz_aa
        real(prec) :: fac, fac_a, fac_b, uz_aa, dz
        real(prec) :: kappa_a, kappa_b, dz1, dz2 
        real(prec), allocatable :: subd(:)      ! nz_aa 
        real(prec), allocatable :: diag(:)      ! nz_aa  
        real(prec), allocatable :: supd(:)      ! nz_aa 
        real(prec), allocatable :: rhs(:)       ! nz_aa 
        real(prec), allocatable :: solution(:)  ! nz_aa
        
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
            dz = H_ice * (zeta_aa(2) - zeta_aa(1))

            ! backward Euler flux basal boundary condition
            subd(1) =  0.0_prec
            diag(1) =  1.0_prec
            supd(1) = -1.0_prec
            rhs(1)  = val_base * dz
                
        else if (k_cts .ge. 2) then 
            ! Layer above base is also temperate (with water likely present in the ice),
            ! set K0 dE/dz = 0. To do so, set basal enthalpy equal to enthalpy above

            subd(1) =  0.0_prec
            diag(1) =  1.0_prec
            supd(1) = -1.0_prec
            rhs(1)  =  0.0_prec
              
        else
            ! Impose basal enthalpy (Dirichlet condition) 

            subd(1) = 0.0_prec
            diag(1) = 1.0_prec
            supd(1) = 0.0_prec
            rhs(1)  = (val_base - enth_ref)

        end if 

        ! == Ice interior layers 2:nz_aa-1 ==

        do k = 2, nz_aa-1

            ! Get implicit vertical advection term, ac => aa nodes
            uz_aa   = 0.5*(uz(k-1)+uz(k))

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
            
            ! Vertical distance for centered difference advection scheme
            dz      =  H_ice*(zeta_aa(k+1)-zeta_aa(k-1))
            
            fac_a   = -kappa_a*dzeta_a(k)*dt/H_ice**2
            fac_b   = -kappa_b*dzeta_b(k)*dt/H_ice**2

            subd(k) = fac_a - uz_aa * dt/dz
            supd(k) = fac_b + uz_aa * dt/dz
            diag(k) = 1.0_prec - fac_a - fac_b
            rhs(k)  = (enth(k)-enth_ref) - dt*advecxy(k) + dt*Q_strn(k)
            
        end do 

        ! == Ice surface ==

        subd(nz_aa) = 0.0_prec
        diag(nz_aa) = 1.0_prec
        supd(nz_aa) = 0.0_prec
        rhs(nz_aa)  = (val_srf-enth_ref)

        ! == Call solver ==

        call solve_tridiag(subd,diag,supd,rhs,solution)

        ! Copy the solution into the temperature variable

        enth = solution + enth_ref 

        return 

    end subroutine calc_enth_column_internal

    subroutine calc_enth_column0(enth,T_ice,omega,bmb_grnd,Q_ice_b,H_cts,T_pmp,cp,kt,advecxy,uz,Q_strn,Q_b, &
                                Q_lith,T_srf,T_shlf,H_ice,H_w,f_grnd,zeta_aa,zeta_ac,dzeta_a,dzeta_b,cr,omega_max,T0,dt)
        ! Thermodynamics solver for a given column of ice 
        ! Note zeta=height, k=1 base, k=nz surface 
        ! Note: nz = number of vertical boundaries (including zeta=0.0 and zeta=1.0), 
        ! temperature is defined for cell centers, plus a value at the surface and the base
        ! so nz_ac = nz_aa - 1 

        ! For notes on implicit form of advection terms, see eg http://farside.ph.utexas.edu/teaching/329/lectures/node90.html
        
        implicit none 

        real(prec), intent(INOUT) :: enth(:)        ! nz_aa [J kg] Ice column enthalpy
        real(prec), intent(INOUT) :: T_ice(:)       ! nz_aa [K] Ice column temperature
        real(prec), intent(INOUT) :: omega(:)       ! nz_aa [-] Ice column water content fraction
        real(prec), intent(INOUT) :: bmb_grnd       ! [m a-1] Basal mass balance (melting is negative)
        real(prec), intent(OUT)   :: Q_ice_b        ! [J a-1 m-2] Ice basal heat flux (positive up)
        real(prec), intent(OUT)   :: H_cts          ! [m] cold-temperate transition surface (CTS) height
        real(prec), intent(IN)    :: T_pmp(:)       ! nz_aa [K] Pressure melting point temp.
        real(prec), intent(IN)    :: cp(:)          ! nz_aa [J kg-1 K-1] Specific heat capacity
        real(prec), intent(IN)    :: kt(:)          ! nz_aa [J a-1 m-1 K-1] Heat conductivity 
        real(prec), intent(IN)    :: advecxy(:)     ! nz_aa [K a-1] Horizontal heat advection 
        real(prec), intent(IN)    :: uz(:)          ! nz_ac [m a-1] Vertical velocity 
        real(prec), intent(IN)    :: Q_strn(:)      ! nz_aa [J a-1 m-3] Internal strain heat production in ice
        real(prec), intent(IN)    :: Q_b            ! [J a-1 m-2] Basal frictional heat production
        real(prec), intent(IN)    :: Q_lith         ! [J a-1 m-2] Bedrock heat flux (positive up)
        real(prec), intent(IN)    :: T_srf          ! [K] Surface temperature 
        real(prec), intent(IN)    :: T_shlf         ! [K] Marine-shelf interface temperature
        real(prec), intent(IN)    :: H_ice          ! [m] Ice thickness 
        real(prec), intent(IN)    :: H_w            ! [m] Basal water layer thickness 
        real(prec), intent(IN)    :: f_grnd         ! [--] Grounded fraction
        real(prec), intent(IN)    :: zeta_aa(:)     ! nz_aa [--] Vertical sigma coordinates (zeta==height), layer centered aa-nodes
        real(prec), intent(IN)    :: zeta_ac(:)     ! nz_ac [--] Vertical height axis temperature (0:1), layer edges ac-nodes    
        real(prec), intent(IN)    :: dzeta_a(:)     ! nz_aa [--] Solver discretization helper variable ak
        real(prec), intent(IN)    :: dzeta_b(:)     ! nz_aa [--] Solver discretization helper variable bk
        real(prec), intent(IN)    :: cr             ! [--] Conductivity ratio (kappa_water / kappa_ice)
        real(prec), intent(IN)    :: omega_max      ! [-] Maximum allowed water fraction inside ice, typically omega_max=0.02 
        real(prec), intent(IN)    :: T0             ! [K or degreesCelcius] Reference melting temperature  
        real(prec), intent(IN)    :: dt             ! [a] Time step 
        
        ! Local variables 
        integer    :: k, nz_aa, nz_ac, k_cts 
        real(prec) :: Q_strn_now
        real(prec) :: H_w_predicted
        real(prec) :: T_excess
        real(prec) :: melt_internal   
        real(prec) :: enth_b, enth_pmp_b 
        real(prec) :: dedz 
        real(prec) :: omega_excess

!         real(prec), allocatable :: dzeta_a(:)   ! nz_aa [--] Solver discretization helper variable ak
!         real(prec), allocatable :: dzeta_b(:)   ! nz_aa [--] Solver discretization helper variable bk

        real(prec), allocatable :: fac_enth(:)  ! aa-nodes 
        real(prec), allocatable :: var(:)       ! aa-nodes 
        real(prec), allocatable :: enth_pmp(:)  ! aa-nodes
        real(prec), allocatable :: kappa_aa(:)  ! aa-nodes

        real(prec), allocatable :: subd(:)      ! nz_aa 
        real(prec), allocatable :: diag(:)      ! nz_aa  
        real(prec), allocatable :: supd(:)      ! nz_aa 
        real(prec), allocatable :: rhs(:)       ! nz_aa 
        real(prec), allocatable :: solution(:)  ! nz_aa

        real(prec) :: fac, fac_a, fac_b, uz_aa, dzeta, dz, dz1, dz2 
        real(prec) :: kappa_a, kappa_b 
        logical    :: use_enth  

        logical, parameter :: test_expl_advecz = .FALSE. 
        real(prec), allocatable :: advecz(:) 

        nz_aa = size(zeta_aa,1)
        nz_ac = size(zeta_ac,1)

!         allocate(dzeta_a(nz_aa))
!         allocate(dzeta_b(nz_aa))

        allocate(kappa_aa(nz_aa))
        allocate(fac_enth(nz_aa))
        allocate(var(nz_aa))
        allocate(enth_pmp(nz_aa))

        allocate(subd(nz_aa))
        allocate(diag(nz_aa))
        allocate(supd(nz_aa))
        allocate(rhs(nz_aa))
        allocate(solution(nz_aa))

        ! Define dzeta terms for this column
        ! Note: for constant zeta axis, this can be done once outside
        ! instead of for each column. However, it is done here to allow
        ! use of adaptive vertical axis.
!         call calc_dzeta_terms(dzeta_a,dzeta_b,zeta_aa,zeta_ac)

        ! Get geothermal heat flux in proper units 
        ! Q_geo_now = Q_geo*1e-3*sec_year   ! [mW m-2] => [J m-2 a-1]

        ! Step 0: Calculate diffusivity, set prognostic variable (T_ice or enth),
        ! and corresponding scaling factor (fac_enth)

        fac_enth = cp               ! To scale to units of [J kg]
        var      = enth             ! [J kg]

        enth_pmp = T_pmp*fac_enth 

        call calc_enth_diffusivity(kappa_aa,enth,enth_pmp,cp,kt,cr,rho_ice)

        ! Step 1: Apply vertical implicit diffusion-advection
        
        ! Step 1: Apply vertical advection (for explicit testing)
        if (test_expl_advecz) then
            allocate(advecz(nz_aa))
            advecz = 0.0
            call calc_advec_vertical_column(advecz,var,uz,H_ice,zeta_aa)
            var = var - dt*advecz
        end if

        ! == Ice base ==

        if (f_grnd .lt. 1.0) then
            ! Floating or partially floating ice - set temperature equal 
            ! to basal temperature at pressure melting point, or marine freezing temp,
            ! or weighted average between the two.

            ! Impose the weighted average of the pressure melting point and the marine freezing temp.
            subd(1) = 0.0_prec
            diag(1) = 1.0_prec
            supd(1) = 0.0_prec
            rhs(1)  = (f_grnd*T_pmp(1) + (1.0-f_grnd)*T_shlf) * fac_enth(1)

        else 
            ! Grounded ice 

            ! Determine expected basal water thickness [m] for this timestep,
            ! using basal mass balance from previous time step (good guess)
            H_w_predicted = H_w - (bmb_grnd*(rho_w/rho_ice))*dt  

            ! == Assign grounded basal boundary conditions ==

            if (T_ice(1) .lt. T_pmp(1) .or. H_w_predicted .lt. 0.0_prec) then   
                ! Frozen at bed, or about to become frozen 

                ! Calculate dzeta for the bottom layer between the basal boundary
                ! (ac-node) and the centered (aa-node) temperature point above
                ! Note: zeta_aa(1) == zeta_ac(1) == bottom boundary 
                dzeta = zeta_aa(2) - zeta_aa(1)

                ! Backward Euler flux basal boundary condition
                subd(1) =  0.0_prec
                diag(1) =  1.0_prec
                supd(1) = -1.0_prec
                rhs(1)  = ((Q_b + Q_lith) * dzeta*H_ice / kt(1)) * fac_enth(1)
                
            else 
                ! Temperate at bed 
                ! Hold basal temperature at pressure melting point

                if (T_ice(2) .ge. T_pmp(2)) then 
                    ! Layer above base is also temperate (with water likely present in the ice),
                    ! set K0 dE/dz = 0. To do so, set basal enthalpy equal to enthalpy above

                    subd(1) =  0.0_prec
                    diag(1) =  1.0_prec
                    supd(1) = -1.0_prec
                    rhs(1)  =  0.0_prec
                    
                else 
                    ! Set enthalpy/temp equal to pressure melting point value 

                    subd(1) = 0.0_prec
                    diag(1) = 1.0_prec
                    supd(1) = 0.0_prec
                    rhs(1)  = T_pmp(1) * fac_enth(1)

                end if 

            end if   ! melting or frozen

        end if  ! floating or grounded 

        ! == Ice interior layers 2:nz_aa-1 ==

        ! Find height of CTS - highest temperate layer 
        k_cts = get_cts_index(enth,T_pmp*cp)

        do k = 2, nz_aa-1

            if (test_expl_advecz) then 
            
                uz_aa = 0.0_prec 

            else 
                ! Implicit vertical advection term on aa-node
                
                uz_aa   = 0.5_prec*(uz(k-1)+uz(k))   ! ac => aa nodes
            
            end if 

            ! Convert units of Q_strn [J a-1 m-3] => [K a-1]
            Q_strn_now = Q_strn(k)/(rho_ice*cp(k))

            ! Get kappa for the lower and upper ac-nodes 
            ! Note: this is important to avoid mixing of kappa at the 
            ! CTS height (kappa_lower = kappa_temperate; kappa_upper = kappa_cold)
            ! See Blatter and Greve, 2015, Eq. 25. 
            !kappa_a = 0.5_prec*(kappa_aa(k-1) + kappa_aa(k))
            !kappa_b = 0.5_prec*(kappa_aa(k)   + kappa_aa(k+1))

            dz1 = zeta_ac(k-1)-zeta_aa(k-1)
            dz2 = zeta_aa(k)-zeta_ac(k-1)
            call calc_wtd_harmonic_mean(kappa_a,kappa_aa(k-1),kappa_aa(k),dz1,dz2)

            dz1 = zeta_ac(k)-zeta_aa(k)
            dz2 = zeta_aa(k+1)-zeta_ac(k)
            call calc_wtd_harmonic_mean(kappa_b,kappa_aa(k),kappa_aa(k+1),dz1,dz2)

            if (k .eq. k_cts+1) kappa_a = kappa_aa(k-1)
            !if (k .eq. k_cts)   kappa_b = kappa_aa(k+1) 

            ! Vertical distance for centered difference advection scheme
            dz      =  H_ice*(zeta_aa(k+1)-zeta_aa(k-1))
            
            fac_a   = -kappa_a*dzeta_a(k)*dt/H_ice**2
            fac_b   = -kappa_b*dzeta_b(k)*dt/H_ice**2

            subd(k) = fac_a - uz_aa * dt/dz
            diag(k) = 1.0_prec - fac_a - fac_b
            supd(k) = fac_b + uz_aa * dt/dz
            rhs(k)  = var(k) - dt*advecxy(k) + dt*Q_strn_now*fac_enth(k)
            
        end do 

        ! == Ice surface ==

        subd(nz_aa) = 0.0_prec
        diag(nz_aa) = 1.0_prec
        supd(nz_aa) = 0.0_prec
        rhs(nz_aa)  = min(T_srf,T0) * fac_enth(nz_aa)

        ! == Call solver ==

        call solve_tridiag(subd,diag,supd,rhs,solution)


        ! Copy the solution into the enthalpy variable,
        ! recalculate enthalpy, temperature and water content 
        
        enth  = solution

        ! Modify enthalpy at the base in the case that a temperate layer is present above the base
        ! (water content should increase towards the base)
        if (enth(2) .ge. T_pmp(2)*cp(2)) then 
            ! Temperate layer exists, interpolate enthalpy at the base. 

!             dedz    = (enth(3)-enth(2))/(zeta_aa(3)-zeta_aa(2))
!             enth(1) = enth(2) + dedz*(zeta_aa(1)-zeta_aa(2))
            
            enth(1) = enth(2)
        end if 
        
        ! Calculate heat flux at ice base as enthalpy gradient * rho_ice * diffusivity [J a-1 m-2]
        if (H_ice .gt. 0.0_prec) then 
            dz = H_ice * (zeta_aa(2)-zeta_aa(1))
            Q_ice_b = kappa_aa(1) * rho_ice * (enth(2) - enth(1)) / dz
        else
            Q_ice_b = 0.0 
        end if 

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

!         ! Calculate heat flux at ice base as enthalpy gradient * rho_ice * diffusivity [J a-1 m-2]
!         if (H_ice .gt. 0.0_prec) then 
!             dz = H_ice * (zeta_aa(2)-zeta_aa(1))
!             Q_ice_b = kappa_aa(1) * rho_ice * (enth(2) - enth(1)) / dz
!         else
!             Q_ice_b = 0.0 
!         end if 

        ! Calculate basal mass balance 
        call calc_bmb_grounded_enth(bmb_grnd,T_ice(1)-T_pmp(1),omega(1),Q_ice_b,Q_b,Q_lith,f_grnd,rho_ice)
        
        ! Include internal melting in bmb_grnd 
        bmb_grnd = bmb_grnd - melt_internal 

! ======================= Corrector step for cold ice ==========================
if (.FALSE.) then 

        ! Find height of CTS - heighest temperate layer 
        k_cts = get_cts_index(enth,T_pmp*cp)

        if (k_cts .ge. 2) then
            ! Temperate ice exists above the base, recalculate cold layers 

            ! Recalculate diffusivity (only relevant for cold points)
            call calc_enth_diffusivity(kappa_aa,enth,enth_pmp,cp,kt,cr,rho_ice)

            ! Lower boundary condition for cold ice dE/dz = 0.0 

            subd(k_cts) = 0.0_prec
            diag(k_cts) = 1.0_prec
            supd(k_cts) = 0.0_prec
            rhs(k_cts)  = enth(k_cts+1)

!             subd(k_cts+1) =  1.0_prec
!             diag(k_cts+1) = -1.0_prec
!             supd(k_cts+1) =  0.0_prec
!             rhs(k_cts+1)  =  0.0_prec
    
            ! == Cold ice interior layers k_cts:nz_aa-1 ==

            do k = k_cts+1, nz_aa-1

                ! Implicit vertical advection term on aa-node
                uz_aa   = 0.5*(uz(k-1)+uz(k))   ! ac => aa nodes

                ! Convert units of Q_strn [J a-1 m-3] => [K a-1]
                Q_strn_now = Q_strn(k)/(rho_ice*cp(k))

                ! Get kappa for the lower and upper ac-nodes 
                ! Note: this is important to avoid mixing of kappa at the 
                ! CTS height (kappa_lower = kappa_temperate; kappa_upper = kappa_cold)
                ! See Blatter and Greve, 2015, Eq. 25. 
                !kappa_a = kappa_aa(k-1)
                !kappa_b = kappa_aa(k) 

                dz1 = zeta_ac(k-1)-zeta_aa(k-1)
                dz2 = zeta_aa(k)-zeta_ac(k-1)
                call calc_wtd_harmonic_mean(kappa_a,kappa_aa(k-1),kappa_aa(k),dz1,dz2)

                dz1 = zeta_ac(k)-zeta_aa(k)
                dz2 = zeta_aa(k+1)-zeta_ac(k)
                call calc_wtd_harmonic_mean(kappa_b,kappa_aa(k),kappa_aa(k+1),dz1,dz2)

                !if (k .eq. k_cts+1) kappa_a = kappa_aa(k-1)
                if (k .eq. k_cts+1) kappa_a = 0.0_prec 

                ! Vertical distance for centered difference advection scheme
                dz      =  H_ice*(zeta_aa(k+1)-zeta_aa(k-1))
                
                fac_a   = -kappa_a*dzeta_a(k)*dt/H_ice**2
                fac_b   = -kappa_b*dzeta_b(k)*dt/H_ice**2

                subd(k) = fac_a - uz_aa * dt/dz
                supd(k) = fac_b + uz_aa * dt/dz
                diag(k) = 1.0_prec - fac_a - fac_b
                rhs(k)  = var(k) - dt*advecxy(k) + dt*Q_strn_now*fac_enth(k)
                
            end do 

            ! == Ice surface ==

            subd(nz_aa) = 0.0_prec
            diag(nz_aa) = 1.0_prec
            supd(nz_aa) = 0.0_prec
            rhs(nz_aa)  = min(T_srf,T0) * fac_enth(nz_aa)

            ! == Call solver ==

            call solve_tridiag(subd(k_cts:nz_aa),diag(k_cts:nz_aa),supd(k_cts:nz_aa), &
                                        rhs(k_cts:nz_aa),solution(k_cts:nz_aa))

            enth(k_cts+1:nz_aa) = solution(k_cts+1:nz_aa) 
            
            ! Get temperature and water content 
            call convert_from_enthalpy_column(enth,T_ice,omega,T_pmp,cp,L_ice)
        
        end if  

end if 
! ==============================================================================




        ! Finally, calculate the CTS height 
        H_cts = calc_cts_height(enth,T_ice,omega,T_pmp,cp,H_ice,zeta_aa)

        return 

    end subroutine calc_enth_column0

    subroutine calc_enth_column_poly(enth,enth_pmp,cp,kt,advecxy,uz,Q_strn, &
                                                            zeta_aa,zeta_ac,H_now,rho_ice,dt)
        ! Thermodynamics solver for a given column of ice 
        ! Note zeta=height, k=1 base, k=nz surface 
        ! Note: nz = number of vertical boundaries (including zeta=0.0 and zeta=1.0), 
        ! temperature is defined for cell centers, plus a value at the surface and the base
        ! so nz_ac = nz_aa - 1 

        ! For notes on implicit form of advection terms, see eg http://farside.ph.utexas.edu/teaching/329/lectures/node90.html
        
        implicit none 

        real(prec), intent(INOUT) :: enth(:)        ! nz_aa [J kg] Ice column enthalpy
        real(prec), intent(IN)    :: enth_pmp(:)    ! nz_aa [K] Pressure melting point temp.
        real(prec), intent(IN)    :: cp(:)          ! nz_aa [J kg-1 K-1] Specific heat capacity
        real(prec), intent(IN)    :: kt(:)          ! nz_aa [J a-1 m-1 K-1] Heat conductivity 
        real(prec), intent(IN)    :: advecxy(:)     ! nz_aa [K a-1] Horizontal heat advection 
        real(prec), intent(IN)    :: uz(:)          ! nz_ac [m a-1] Vertical velocity 
        real(prec), intent(IN)    :: Q_strn(:)      ! nz_aa [J a-1 m-3] Internal strain heat production in ice
        real(prec), intent(IN)    :: zeta_aa(:)     ! nz_aa [--] Vertical sigma coordinates (zeta==height), layer centered aa-nodes
        real(prec), intent(IN)    :: zeta_ac(:)     ! nz_ac [--] Vertical height axis temperature (0:1), layer edges ac-nodes    
        real(prec), intent(IN)    :: H_now          ! [m] Ice thickness of column
        real(prec), intent(IN)    :: rho_ice        ! [kg m-3] Ice density   
        real(prec), intent(IN)    :: dt             ! [a] Time step 
        


        return 

    end subroutine calc_enth_column_poly

    subroutine calc_enth_column_cold(enth,enth_pmp,cp,kt,advecxy,uz,Q_strn, &
                                                            zeta_aa,zeta_ac,H_now,rho_ice,dt)
        ! Thermodynamics solver for a given column of ice 
        ! Note zeta=height, k=1 base, k=nz surface 
        ! Note: nz = number of vertical boundaries (including zeta=0.0 and zeta=1.0), 
        ! temperature is defined for cell centers, plus a value at the surface and the base
        ! so nz_ac = nz_aa - 1 

        ! For notes on implicit form of advection terms, see eg http://farside.ph.utexas.edu/teaching/329/lectures/node90.html
        
        implicit none 

        real(prec), intent(INOUT) :: enth(:)        ! nz_aa [J kg] Ice column enthalpy
        real(prec), intent(IN)    :: enth_pmp(:)    ! nz_aa [K] Pressure melting point temp.
        real(prec), intent(IN)    :: cp(:)          ! nz_aa [J kg-1 K-1] Specific heat capacity
        real(prec), intent(IN)    :: kt(:)          ! nz_aa [J a-1 m-1 K-1] Heat conductivity 
        real(prec), intent(IN)    :: advecxy(:)     ! nz_aa [K a-1] Horizontal heat advection 
        real(prec), intent(IN)    :: uz(:)          ! nz_ac [m a-1] Vertical velocity 
        real(prec), intent(IN)    :: Q_strn(:)      ! nz_aa [J a-1 m-3] Internal strain heat production in ice
        real(prec), intent(IN)    :: zeta_aa(:)     ! nz_aa [--] Vertical sigma coordinates (zeta==height), layer centered aa-nodes
        real(prec), intent(IN)    :: zeta_ac(:)     ! nz_ac [--] Vertical height axis temperature (0:1), layer edges ac-nodes    
        real(prec), intent(IN)    :: H_now          ! [m] Ice thickness of column
        real(prec), intent(IN)    :: rho_ice        ! [kg m-3] Ice density   
        real(prec), intent(IN)    :: dt             ! [a] Time step 
        
        ! Local variables 
        integer    :: k, nz_aa
        real(prec) :: Q_strn_now

        real(prec), allocatable :: dzeta_a(:)   ! nz_aa [--] Solver discretization helper variable ak
        real(prec), allocatable :: dzeta_b(:)   ! nz_aa [--] Solver discretization helper variable bk

        real(prec), allocatable :: kappa_aa(:)  ! aa-nodes

        real(prec), allocatable :: subd(:)      ! nz_aa 
        real(prec), allocatable :: diag(:)      ! nz_aa  
        real(prec), allocatable :: supd(:)      ! nz_aa 
        real(prec), allocatable :: rhs(:)       ! nz_aa 
        real(prec), allocatable :: solution(:)  ! nz_aa

        real(prec) :: fac, fac_a, fac_b, uz_aa, dzeta, dz, dz1, dz2 
        real(prec) :: kappa_a, kappa_b 

        nz_aa = size(zeta_aa,1)

        allocate(dzeta_a(nz_aa))
        allocate(dzeta_b(nz_aa))

        allocate(kappa_aa(nz_aa))

        allocate(subd(nz_aa))
        allocate(diag(nz_aa))
        allocate(supd(nz_aa))
        allocate(rhs(nz_aa))
        allocate(solution(nz_aa))

        ! Define dzeta terms for this column
        ! Note: for constant zeta axis, this can be done once outside
        ! instead of for each column. However, it is done here to allow
        ! use of adaptive vertical axis.
        call calc_dzeta_terms(dzeta_a,dzeta_b,zeta_aa,zeta_ac)

        call calc_enth_diffusivity(kappa_aa,enth,enth_pmp,cp,kt,cr=0.0_prec,rho_ice=rho_ice)

        ! == Base of cold region - prescribe enthalpy of the pressure melting point ==

        subd(1) = 0.0_prec
        diag(1) = 1.0_prec
        supd(1) = 0.0_prec
        rhs(1)  = enth_pmp(1)

        ! == Ice interior layers 2:nz_aa-1 ==

        do k = 2, nz_aa-1
 
            ! Implicit vertical advection term on aa-node
            uz_aa   = 0.5_prec*(uz(k-1)+uz(k))   ! ac => aa-nodes
            
            ! Convert units of Q_strn [J a-1 m-3] => [K a-1]
            Q_strn_now = Q_strn(k)/(rho_ice*cp(k))

            ! Get kappa for the lower and upper ac-nodes 
            ! Note: this is important to avoid mixing of kappa at the 
            ! CTS height (kappa_lower = kappa_temperate; kappa_upper = kappa_cold)
            ! See Blatter and Greve, 2015, Eq. 25. 
            !kappa_a = 0.5_prec*(kappa_aa(k-1) + kappa_aa(k))
            !kappa_b = 0.5_prec*(kappa_aa(k)   + kappa_aa(k+1))

            dz1 = zeta_ac(k-1)-zeta_aa(k-1)
            dz2 = zeta_aa(k)-zeta_ac(k-1)
            call calc_wtd_harmonic_mean(kappa_a,kappa_aa(k-1),kappa_aa(k),dz1,dz2)

            dz1 = zeta_ac(k)-zeta_aa(k)
            dz2 = zeta_aa(k+1)-zeta_ac(k)
            call calc_wtd_harmonic_mean(kappa_b,kappa_aa(k),kappa_aa(k+1),dz1,dz2)

            ! Vertical distance for centered difference advection scheme
            dz      =  H_now*(zeta_aa(k+1)-zeta_aa(k-1))
            
            fac_a   = -kappa_a*dzeta_a(k)*dt/H_now**2
            fac_b   = -kappa_b*dzeta_b(k)*dt/H_now**2

            subd(k) = fac_a - uz_aa * dt/dz
            diag(k) = 1.0_prec - fac_a - fac_b
            supd(k) = fac_b + uz_aa * dt/dz
            rhs(k)  = enth(k) - dt*advecxy(k) + dt*Q_strn_now*cp(k)
            
        end do 

        ! == Surface of cold region (ice surface) ==

        subd(nz_aa) = 0.0_prec
        diag(nz_aa) = 1.0_prec
        supd(nz_aa) = 0.0_prec
        rhs(nz_aa)  = min(enth(nz_aa),enth_pmp(nz_aa))

        ! == Call solver ==

        call solve_tridiag(subd,diag,supd,rhs,solution)

        ! Copy the solution into the enthalpy variable for output
        
        enth  = solution

        return 

    end subroutine calc_enth_column_cold

    subroutine calc_enth_column_temperate(enth,enth_pmp,cp,kt,advecxy,uz,Q_strn, &
                                                            zeta_aa,zeta_ac,H_now,rho_ice,dt)
        ! Thermodynamics solver for a given column of ice 
        ! Note zeta=height, k=1 base, k=nz surface 
        ! Note: nz = number of vertical boundaries (including zeta=0.0 and zeta=1.0), 
        ! temperature is defined for cell centers, plus a value at the surface and the base
        ! so nz_ac = nz_aa - 1 

        ! For notes on implicit form of advection terms, see eg http://farside.ph.utexas.edu/teaching/329/lectures/node90.html
        
        implicit none 

        real(prec), intent(INOUT) :: enth(:)        ! nz_aa [J kg] Ice column enthalpy
        real(prec), intent(IN)    :: enth_pmp(:)    ! nz_aa [K] Pressure melting point temp.
        real(prec), intent(IN)    :: cp(:)          ! nz_aa [J kg-1 K-1] Specific heat capacity
        real(prec), intent(IN)    :: kt(:)          ! nz_aa [J a-1 m-1 K-1] Heat conductivity 
        real(prec), intent(IN)    :: advecxy(:)     ! nz_aa [K a-1] Horizontal heat advection 
        real(prec), intent(IN)    :: uz(:)          ! nz_ac [m a-1] Vertical velocity 
        real(prec), intent(IN)    :: Q_strn(:)      ! nz_aa [J a-1 m-3] Internal strain heat production in ice
        real(prec), intent(IN)    :: zeta_aa(:)     ! nz_aa [--] Vertical sigma coordinates (zeta==height), layer centered aa-nodes
        real(prec), intent(IN)    :: zeta_ac(:)     ! nz_ac [--] Vertical height axis temperature (0:1), layer edges ac-nodes    
        real(prec), intent(IN)    :: H_now          ! [m] Ice thickness of column
        real(prec), intent(IN)    :: rho_ice        ! [kg m-3] Ice density   
        real(prec), intent(IN)    :: dt             ! [a] Time step 
        
        ! Local variables 
        integer    :: k, nz_aa
        real(prec) :: Q_strn_now

        real(prec), allocatable :: dzeta_a(:)   ! nz_aa [--] Solver discretization helper variable ak
        real(prec), allocatable :: dzeta_b(:)   ! nz_aa [--] Solver discretization helper variable bk

        real(prec), allocatable :: kappa_aa(:)  ! aa-nodes

        real(prec), allocatable :: subd(:)      ! nz_aa 
        real(prec), allocatable :: diag(:)      ! nz_aa  
        real(prec), allocatable :: supd(:)      ! nz_aa 
        real(prec), allocatable :: rhs(:)       ! nz_aa 
        real(prec), allocatable :: solution(:)  ! nz_aa

        real(prec) :: fac, fac_a, fac_b, uz_aa, dzeta, dz, dz1, dz2 
        real(prec) :: kappa_a, kappa_b 

        nz_aa = size(zeta_aa,1)

        allocate(dzeta_a(nz_aa))
        allocate(dzeta_b(nz_aa))

        allocate(kappa_aa(nz_aa))

        allocate(subd(nz_aa))
        allocate(diag(nz_aa))
        allocate(supd(nz_aa))
        allocate(rhs(nz_aa))
        allocate(solution(nz_aa))

        ! Define dzeta terms for this column
        ! Note: for constant zeta axis, this can be done once outside
        ! instead of for each column. However, it is done here to allow
        ! use of adaptive vertical axis.
        call calc_dzeta_terms(dzeta_a,dzeta_b,zeta_aa,zeta_ac)

        call calc_enth_diffusivity(kappa_aa,enth,enth_pmp,cp,kt,cr=0.0_prec,rho_ice=rho_ice)

        ! == Ice base ==

        ! Temperate at bed 
        ! Hold basal temperature at pressure melting point

        if (enth(2) .ge. enth_pmp(2)) then 
            ! Layer above base is also temperate (with water likely present in the ice),
            ! set K0 dE/dz = 0. To do so, set basal enthalpy equal to enthalpy above

            subd(1) =  0.0_prec
            diag(1) =  1.0_prec
            supd(1) = -1.0_prec
            rhs(1)  =  0.0_prec
            
        else 
            ! Set enthalpy equal to pressure melting point value 

            subd(1) = 0.0_prec
            diag(1) = 1.0_prec
            supd(1) = 0.0_prec
            rhs(1)  = enth_pmp(1)

        end if 

        ! == Ice interior layers 2:nz_aa-1 ==

        do k = 2, nz_aa-1
 
            ! Implicit vertical advection term on aa-node    
            uz_aa   = 0.5_prec*(uz(k-1)+uz(k))   ! ac => aa nodes
            
            ! Convert units of Q_strn [J a-1 m-3] => [K a-1]
            Q_strn_now = Q_strn(k)/(rho_ice*cp(k))

            ! Get kappa for the lower and upper ac-nodes 
            ! Note: this is important to avoid mixing of kappa at the 
            ! CTS height (kappa_lower = kappa_temperate; kappa_upper = kappa_cold)
            ! See Blatter and Greve, 2015, Eq. 25. 
            !kappa_a = 0.5_prec*(kappa_aa(k-1) + kappa_aa(k))
            !kappa_b = 0.5_prec*(kappa_aa(k)   + kappa_aa(k+1))

            dz1 = zeta_ac(k-1)-zeta_aa(k-1)
            dz2 = zeta_aa(k)-zeta_ac(k-1)
            call calc_wtd_harmonic_mean(kappa_a,kappa_aa(k-1),kappa_aa(k),dz1,dz2)

            dz1 = zeta_ac(k)-zeta_aa(k)
            dz2 = zeta_aa(k+1)-zeta_ac(k)
            call calc_wtd_harmonic_mean(kappa_b,kappa_aa(k),kappa_aa(k+1),dz1,dz2)

            ! Vertical distance for centered difference advection scheme
            dz      =  H_now*(zeta_aa(k+1)-zeta_aa(k-1))
            
            fac_a   = -kappa_a*dzeta_a(k)*dt/H_now**2
            fac_b   = -kappa_b*dzeta_b(k)*dt/H_now**2

            subd(k) = fac_a - uz_aa * dt/dz
            diag(k) = 1.0_prec - fac_a - fac_b
            supd(k) = fac_b + uz_aa * dt/dz
            rhs(k)  = enth(k) - dt*advecxy(k) + dt*Q_strn_now*cp(k)
            
        end do 

        ! == Ice surface, temperate layer ==

        subd(nz_aa) = 0.0_prec
        diag(nz_aa) = 1.0_prec
        supd(nz_aa) = 0.0_prec
        rhs(nz_aa)  = enth_pmp(nz_aa)

        ! == Call solver ==

        call solve_tridiag(subd,diag,supd,rhs,solution)


        ! Copy the solution into the enthalpy variable
        
        enth  = solution

        ! If temperate layer exists, ensure basal boundary condition 
        ! holds dE/dz = 0 == E(1) = E(2);
        ! This is only for extra security w.r.t. the solver stability
        
        if (enth(2) .ge. enth_pmp(2)) enth(1) = enth(2)
        
        return 

    end subroutine calc_enth_column_temperate

    subroutine calc_enth_column_zoom(enth,enth_pmp,cp,kt,advecxy,uz,Q_strn, &
                                                            zeta_aa,zeta_ac,cr,H_now,rho_ice,dt)
        ! Thermodynamics solver for a given column of ice 
        ! Note zeta=height, k=1 base, k=nz surface 
        ! Note: nz = number of vertical boundaries (including zeta=0.0 and zeta=1.0), 
        ! temperature is defined for cell centers, plus a value at the surface and the base
        ! so nz_ac = nz_aa - 1 

        ! For notes on implicit form of advection terms, see eg http://farside.ph.utexas.edu/teaching/329/lectures/node90.html
        
        implicit none 

        real(prec), intent(INOUT) :: enth(:)        ! nz_aa [J kg] Ice column enthalpy
        real(prec), intent(IN)    :: enth_pmp(:)    ! nz_aa [K] Pressure melting point temp.
        real(prec), intent(IN)    :: cp(:)          ! nz_aa [J kg-1 K-1] Specific heat capacity
        real(prec), intent(IN)    :: kt(:)          ! nz_aa [J a-1 m-1 K-1] Heat conductivity 
        real(prec), intent(IN)    :: advecxy(:)     ! nz_aa [K a-1] Horizontal heat advection 
        real(prec), intent(IN)    :: uz(:)          ! nz_ac [m a-1] Vertical velocity 
        real(prec), intent(IN)    :: Q_strn(:)      ! nz_aa [J a-1 m-3] Internal strain heat production in ice
        real(prec), intent(IN)    :: zeta_aa(:)     ! nz_aa [--] Vertical sigma coordinates (zeta==height), layer centered aa-nodes
        real(prec), intent(IN)    :: zeta_ac(:)     ! nz_ac [--] Vertical height axis temperature (0:1), layer edges ac-nodes    
        real(prec), intent(IN)    :: cr 
        real(prec), intent(IN)    :: H_now          ! [m] Ice thickness of column
        real(prec), intent(IN)    :: rho_ice        ! [kg m-3] Ice density   
        real(prec), intent(IN)    :: dt             ! [a] Time step 
        
        ! Local variables 
        integer    :: k, nz_aa
        real(prec) :: Q_strn_now

        real(prec), allocatable :: dzeta_a(:)   ! nz_aa [--] Solver discretization helper variable ak
        real(prec), allocatable :: dzeta_b(:)   ! nz_aa [--] Solver discretization helper variable bk

        real(prec), allocatable :: kappa_aa(:)  ! aa-nodes

        real(prec), allocatable :: subd(:)      ! nz_aa 
        real(prec), allocatable :: diag(:)      ! nz_aa  
        real(prec), allocatable :: supd(:)      ! nz_aa 
        real(prec), allocatable :: rhs(:)       ! nz_aa 
        real(prec), allocatable :: solution(:)  ! nz_aa

        real(prec) :: fac, fac_a, fac_b, uz_aa, dzeta, dz, dz1, dz2 
        real(prec) :: kappa_a, kappa_b 

        nz_aa = size(zeta_aa,1)

        allocate(dzeta_a(nz_aa))
        allocate(dzeta_b(nz_aa))

        allocate(kappa_aa(nz_aa))

        allocate(subd(nz_aa))
        allocate(diag(nz_aa))
        allocate(supd(nz_aa))
        allocate(rhs(nz_aa))
        allocate(solution(nz_aa))

        ! Define dzeta terms for this column
        ! Note: for constant zeta axis, this can be done once outside
        ! instead of for each column. However, it is done here to allow
        ! use of adaptive vertical axis.
        call calc_dzeta_terms(dzeta_a,dzeta_b,zeta_aa,zeta_ac)

        call calc_enth_diffusivity(kappa_aa,enth,enth_pmp,cp,kt,cr,rho_ice)

        ! == Base of zoom region ==

        subd(1) = 0.0_prec
        diag(1) = 1.0_prec
        supd(1) = 0.0_prec
        rhs(1)  = enth(1)

        ! == Ice interior layers 2:nz_aa-1 ==

        do k = 2, nz_aa-1
 
            ! Implicit vertical advection term on aa-node
            uz_aa   = 0.5_prec*(uz(k-1)+uz(k))   ! ac => aa-nodes
            
            ! Convert units of Q_strn [J a-1 m-3] => [K a-1]
            Q_strn_now = Q_strn(k)/(rho_ice*cp(k))

            ! Get kappa for the lower and upper ac-nodes 
            ! Note: this is important to avoid mixing of kappa at the 
            ! CTS height (kappa_lower = kappa_temperate; kappa_upper = kappa_cold)
            ! See Blatter and Greve, 2015, Eq. 25. 
            !kappa_a = 0.5_prec*(kappa_aa(k-1) + kappa_aa(k))
            !kappa_b = 0.5_prec*(kappa_aa(k)   + kappa_aa(k+1))

            dz1 = zeta_ac(k-1)-zeta_aa(k-1)
            dz2 = zeta_aa(k)-zeta_ac(k-1)
            call calc_wtd_harmonic_mean(kappa_a,kappa_aa(k-1),kappa_aa(k),dz1,dz2)

            dz1 = zeta_ac(k)-zeta_aa(k)
            dz2 = zeta_aa(k+1)-zeta_ac(k)
            call calc_wtd_harmonic_mean(kappa_b,kappa_aa(k),kappa_aa(k+1),dz1,dz2)

            ! Vertical distance for centered difference advection scheme
            dz      =  H_now*(zeta_aa(k+1)-zeta_aa(k-1))
            
            fac_a   = -kappa_a*dzeta_a(k)*dt/H_now**2
            fac_b   = -kappa_b*dzeta_b(k)*dt/H_now**2

            subd(k) = fac_a - uz_aa * dt/dz
            diag(k) = 1.0_prec - fac_a - fac_b
            supd(k) = fac_b + uz_aa * dt/dz
            rhs(k)  = enth(k) - dt*advecxy(k) + dt*Q_strn_now*cp(k)
            
        end do 

        ! == Surface of zoom region ==

        subd(nz_aa) = 0.0_prec
        diag(nz_aa) = 1.0_prec
        supd(nz_aa) = 0.0_prec
        rhs(nz_aa)  = enth(nz_aa)

        ! == Call solver ==

        call solve_tridiag(subd,diag,supd,rhs,solution)

        ! Copy the solution into the enthalpy variable for output
        
        enth  = solution

        return 

    end subroutine calc_enth_column_zoom

    ! ========== ENTHALPY ==========================================

    elemental subroutine convert_to_enthalpy(enth,temp,omega,T_pmp,cp,L)
        ! Given temperature and water content, calculate enthalpy.

        implicit none 

        real(prec), intent(OUT) :: enth             ! [J m-3] Enthalpy 
        real(prec), intent(IN)  :: temp             ! [K] Temperature 
        real(prec), intent(IN)  :: omega            ! [-] Water content (fraction)
        real(prec), intent(IN)  :: T_pmp            ! [K] Pressure melting point
        real(prec), intent(IN)  :: cp               ! [J kg-1 K-1] Heat capacity 
        real(prec), intent(IN)  :: L                ! [J kg-1] Latent heat of fusion 
        
        enth = (1.0_prec-omega)*(cp*temp) + omega*(cp*T_pmp + L)

        return 

    end subroutine convert_to_enthalpy

    subroutine convert_from_enthalpy_column(enth,temp,omega,T_pmp,cp,L)
        ! Given enthalpy, calculate temperature and water content. 

        implicit none 

        real(prec), intent(INOUT) :: enth(:)            ! [J m-3] Enthalpy, nz_aa nodes
        real(prec), intent(OUT)   :: temp(:)            ! [K] Temperature, nz_aa nodes  
        real(prec), intent(OUT)   :: omega(:)           ! [-] Water content (fraction), nz_aa nodes 
        real(prec), intent(IN)    :: T_pmp(:)           ! [K] Pressure melting point, nz_aa nodes 
        real(prec), intent(IN)    :: cp(:)              ! [J kg-1 K-1] Heat capacity,nz_aa nodes 
        real(prec), intent(IN)    :: L                  ! [J kg-1] Latent heat of fusion
        
        ! Local variables
        integer    :: k, nz_aa  
        real(prec), allocatable :: enth_pmp(:)  

        nz_aa = size(enth,1)

        allocate(enth_pmp(nz_aa))

        ! Find pressure melting point enthalpy
        enth_pmp = T_pmp * cp 

        ! Column interior and basal layer
        ! Note: although the k=1 is a boundary value with no thickness,
        ! allow it to retain omega to maintain consistency with grid points above.
        do k = 1, nz_aa-1

            if (enth(k) .gt. enth_pmp(k)) then
                ! Temperate ice 
                
                temp(k)  = T_pmp(k)
                omega(k) = (enth(k) - enth_pmp(k)) / L 
             else
                ! Cold ice 

                temp(k)  = enth(k) / cp(k) 
                omega(k) = 0.0_prec

             end if

        end do 

        ! Surface layer 
        if (enth(nz_aa) .ge. enth_pmp(nz_aa)) then 
            ! Temperate surface, reset omega to zero and enth to pmp value 
            
            enth(nz_aa)  = enth_pmp(nz_aa)
            temp(nz_aa)  = enth(nz_aa) / cp(nz_aa)
            omega(nz_aa) = 0.0_prec 
        
        else 
            ! Cold surface, calculate T, and reset omega to zero 
            
            temp(nz_aa)  = enth(nz_aa) / cp(nz_aa)
            omega(nz_aa) = 0.0_prec 
        
        end if 
        
        return 

    end subroutine convert_from_enthalpy_column

    subroutine calc_enth_diffusivity(kappa,enth,enth_pmp,cp,kt,cr,rho_ice)
        ! Calculate the enthalpy vertical diffusivity for use with the diffusion solver:
        ! When water is present in the layer, set kappa=kappa_therm, else kappa=kappa_cold 

        implicit none 

        real(prec), intent(OUT) :: kappa(:)         ! [nz_aa]
        real(prec), intent(IN)  :: enth(:)          ! [nz_aa]
        real(prec), intent(IN)  :: enth_pmp(:)      ! [nz_aa]
        real(prec), intent(IN)  :: cp(:)
        real(prec), intent(IN)  :: kt(:)  
        real(prec), intent(IN)  :: cr 
        real(prec), intent(IN)  :: rho_ice
        
        ! Local variables
        integer     :: k, nz
        real(prec)  :: kappa_cold       ! Cold diffusivity 
        real(prec)  :: kappa_temp       ! Temperate diffusivity 
        
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

        real(prec), intent(IN) :: enth(:) 
        real(prec), intent(IN) :: T_ice(:) 
        real(prec), intent(IN) :: omega(:) 
        real(prec), intent(IN) :: T_pmp(:) 
        real(prec), intent(IN) :: cp(:)
        real(prec), intent(IN) :: H_ice  
        real(prec), intent(IN) :: zeta(:) 
        real(prec) :: H_cts 

        ! Local variables 
        integer :: k, k_cts, nz 
        real(prec) :: f_lin, f_lin_0, dedz0, dedz1, zeta_cts 
        real(prec), allocatable :: enth_pmp(:) 

        integer :: i, n_iter, n_prime
        real(prec), allocatable :: zeta_prime(:) 
        real(prec), allocatable :: enth_prime(:) 
        
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
            H_cts = 0.0_prec 

        else if (k_cts .eq. nz) then 
            ! Whole column is temperate
            H_cts = H_ice

        else 

            ! Assume H_cts lies at center of last temperate cell (aa-node)
!             zeta_cts = zeta(k_cts)

!             ! Assume H_cts lies on ac-node between temperate and cold layers 
!             zeta_cts = 0.5_prec*(zeta(k_cts)+zeta(k_cts+1))

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

! !             i = maxloc(abs(enth_prime),1,mask=enth_prime .lt. 0.0_prec)
! !             f_lin = (zeta_prime(i+1)-zeta_prime(i)) / (enth_prime(i+1)-enth_prime(i))
! !             if (abs(f_lin) .lt. 1e-3) f_lin = 0.0_prec 
! !             zeta_cts = (1.0_prec-f_lin)*zeta_prime(i) 
    
            H_cts    = H_ice*zeta_cts

        end if 

        return 


    end function calc_cts_height

    subroutine calc_dzeta_terms(dzeta_a,dzeta_b,zeta_aa,zeta_ac)
        ! zeta_aa  = depth axis at layer centers (plus base and surface values)
        ! zeta_ac  = depth axis (1: base, nz: surface), at layer boundaries
        ! Calculate ak, bk terms as defined in Hoffmann et al (2018)
        implicit none 

        real(prec), intent(INOUT) :: dzeta_a(:)    ! nz_aa
        real(prec), intent(INOUT) :: dzeta_b(:)    ! nz_aa
        real(prec), intent(IN)    :: zeta_aa(:)    ! nz_aa 
        real(prec), intent(IN)    :: zeta_ac(:)    ! nz_ac == nz_aa-1 

        ! Local variables 
        integer :: k, nz_layers, nz_aa    

        nz_aa = size(zeta_aa)

        ! Note: zeta_aa is calculated outside in the main program 

        ! Initialize dzeta_a/dzeta_b to zero, first and last indices will not be used (end points)
        dzeta_a = 0.0 
        dzeta_b = 0.0 
        
        do k = 2, nz_aa-1 
            dzeta_a(k) = 1.0/ ( (zeta_ac(k) - zeta_ac(k-1)) * (zeta_aa(k) - zeta_aa(k-1)) )
        enddo

        do k = 2, nz_aa-1
            dzeta_b(k) = 1.0/ ( (zeta_ac(k) - zeta_ac(k-1)) * (zeta_aa(k+1) - zeta_aa(k)) )
        end do

        return 

    end subroutine calc_dzeta_terms

    subroutine calc_wtd_harmonic_mean(var_ave,var1,var2,wt1,wt2)

        implicit none 

        real(prec), intent(OUT) :: var_ave 
        real(prec), intent(IN)  :: var1 
        real(prec), intent(IN)  :: var2 
        real(prec), intent(IN)  :: wt1 
        real(prec), intent(IN)  :: wt2 
        
        ! Local variables 
        real(prec), parameter   :: tol = 1e-5 

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

        real(prec),   intent(INOUT) :: zeta_pt(:)
        real(prec),   intent(INOUT) :: zeta_pc(:) 
        character(*), intent(IN)    :: zeta_scale 
        real(prec),   intent(IN)    :: zeta_exp 

        ! Local variables
        integer :: k, nz_pt, nz_pc 

        integer :: nz_ac 
        real(prec), allocatable :: zeta_ac(:) 

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

        real(prec),   intent(INOUT) :: zeta_aa(:) 
        real(prec),   intent(INOUT) :: zeta_ac(:) 
        real(prec),   intent(IN)    :: zeta_pt(:) 
        real(prec),   intent(IN)    :: zeta_pc(:) 
        real(prec),   intent(IN)    :: H_cts 
        real(prec),   intent(IN)    :: H_ice 

        ! Local variables 
        integer    :: k 
        integer    :: nzt, nztc, nzc, nz_aa, nz_ac  
        real(prec) :: f_cts

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
        zeta_ac(1) = 0.0_prec 
        do k = 2, nz_ac-1
            zeta_ac(k) = 0.5_prec * (zeta_aa(k)+zeta_aa(k+1))
        end do 
        zeta_ac(nz_ac) = 1.0_prec 

        return 

    end subroutine calc_zeta_combined

    function get_cts_index(enth,enth_pmp) result(k_cts)

        implicit none 

        real(prec), intent(IN) :: enth(:) 
        real(prec), intent(IN) :: enth_pmp(:) 
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

        real(prec), intent(IN)  :: x0,x1,y0,y1, xout
        real(prec) :: yout
        real(prec) :: alph

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


