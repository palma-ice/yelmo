module ice_enthalpy 
    ! Module contains the ice temperature and basal mass balance (grounded) solution

    use yelmo_defs, only : prec, pi, g, sec_year, rho_ice, rho_sw, rho_w, L_ice  
    use solver_tridiagonal, only : solve_tridiag 
    use thermodynamics, only : calc_bmb_grounded, calc_bmb_grounded_enth, calc_advec_vertical_column

    implicit none
    
    private
    public :: calc_enth_column 
    public :: convert_to_enthalpy
    public :: convert_from_enthalpy_column
    public :: calc_dzeta_terms
    
contains 

    subroutine calc_enth_column(enth,T_ice,omega,bmb_grnd,Q_ice_b,H_cts,T_pmp,cp,kt,advecxy,uz,Q_strn,Q_b,Q_geo, &
                    T_srf,T_shlf,H_ice,H_w,f_grnd,zeta_aa,zeta_ac,dzeta_a,dzeta_b,cr,omega_max,T0,dt,solver)
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
        real(prec), intent(IN)    :: Q_geo          ! [mW m-2] Geothermal heat flux (positive up)
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
        character(len=*), intent(IN) :: solver      ! "enth" or "temp" 
        
        ! Local variables 
        integer    :: k, nz_aa, nz_ac
        real(prec) :: Q_geo_now, ghf_conv 
        real(prec) :: Q_strn_now
        real(prec) :: H_w_predicted
        real(prec) :: T_excess
        real(prec) :: melt_internal   
        real(prec) :: enth_b, enth_pmp_b 
        real(prec) :: omega_excess

        logical, parameter      :: test_expl_advecz = .FALSE. 
        real(prec), allocatable :: advecz(:)   ! nz_aa, for explicit vertical advection solving
        
        real(prec), allocatable :: fac_enth(:)  ! aa-nodes 
        real(prec), allocatable :: var(:)       ! aa-nodes 
        real(prec), allocatable :: kappa_aa(:)  ! aa-nodes

        real(prec), allocatable :: subd(:)      ! nz_aa 
        real(prec), allocatable :: diag(:)      ! nz_aa  
        real(prec), allocatable :: supd(:)      ! nz_aa 
        real(prec), allocatable :: rhs(:)       ! nz_aa 
        real(prec), allocatable :: solution(:)  ! nz_aa
        real(prec) :: fac, fac_a, fac_b, uz_aa, dzeta, dz
        real(prec) :: kappa_a, kappa_b 
        logical    :: use_enth 

        nz_aa = size(zeta_aa,1)
        nz_ac = size(zeta_ac,1)

        allocate(kappa_aa(nz_aa))
        allocate(fac_enth(nz_aa))
        allocate(var(nz_aa))

        allocate(subd(nz_aa))
        allocate(diag(nz_aa))
        allocate(supd(nz_aa))
        allocate(rhs(nz_aa))
        allocate(solution(nz_aa))

        ! Determine which solver to use: enth or temp 
        use_enth = .TRUE. 
        if (trim(solver) .eq. "temp") use_enth = .FALSE. 

        ! Get geothermal heat flux in proper units 
        Q_geo_now = Q_geo*1e-3*sec_year   ! [mW m-2] => [J m-2 a-1]

        ! Get enthalpy to have enth, omega and T_ice all defined and consistent initially
        ! Note: in principle, these quantities should all be available and consistent
        ! when entering the routine, but it ensures that enthalpy is defined if only 
        ! T_ice and omega are known initially.
        call convert_to_enthalpy(enth,T_ice,omega,T_pmp,cp,L_ice)

        ! Step 0: Calculate diffusivity, set prognostic variable (T_ice or enth),
        ! and corresponding scaling factor (fac_enth)

        if (use_enth) then 
            ! Use enthalpy as prognostic variable 

            ! Calculate diffusivity on cell centers (aa-nodes)
            call calc_enth_diffusivity(kappa_aa,T_ice,omega,enth,T_pmp,cp,kt,rho_ice,rho_w,L_ice,cr)

            fac_enth = cp               ! To scale to units of [J kg]
            var      = enth             ! [J kg]

        else 
            ! Use temperature as prognostic variable  

            ! Calculate diffusivity on cell centers (aa-nodes)
            kappa_aa = kt / (rho_ice*cp)
        
            fac_enth = 1.0              ! Keep units of [K]
            var      = T_ice            ! [K]

        end if 

        ! Step 1: Apply vertical advection (for explicit testing)
        if (test_expl_advecz) then 
            allocate(advecz(nz_aa))
            advecz = 0.0
            call calc_advec_vertical_column(advecz,var,uz,H_ice,zeta_aa)
            var = var - dt*advecz 
        end if 

        ! Step 2: Apply vertical implicit diffusion-advection (or diffusion only if test_expl_advecz=True)
        
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

                ! maintain balance of heat sources and sinks
                ! (conductive flux, geothermal flux, and basal friction)

                ! Note: (Q_b+Q_geo_now) is generally >= 0, since defined as positive up

                ! calculate dzeta for the bottom layer between the basal boundary and the temperature point above
                dzeta = zeta_aa(2) - zeta_aa(1)

                ! backward Euler flux basal boundary condition
                subd(1) =  0.0_prec
                diag(1) =  1.0_prec
                supd(1) = -1.0_prec
                rhs(1)  = ((Q_b + Q_geo_now) * dzeta*H_ice / kt(1)) * fac_enth(1)
                
            else 
                ! Temperate at bed 
                ! Hold basal temperature at pressure melting point

                if (use_enth .and. T_ice(2) .ge. T_pmp(2)) then 
                    ! Layer above base is also temperate (with water likely present in the ice),
                    ! set K0 dE/dz = 0. To do so, set basal enthalpy equal to enthalpy above
                    ! (following MALIv6 implementation)

                    subd(1) =  0.0_prec
                    diag(1) =  1.0_prec
                    supd(1) = -1.0_prec
                    rhs(1)  =  0.0_prec
    
                    ! Testing implementation of second-order upwind derivative,
                    ! using var(3) value from previous timestep (doesn't work)
!                     subd(1) =  0.0_prec
!                     diag(1) = -3.0_prec
!                     supd(1) =  4.0_prec
!                     rhs(1)  =  var(3)

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

        do k = 2, nz_aa-1

            if (test_expl_advecz) then 
                ! No implicit vertical advection (diffusion only)
                uz_aa = 0.0 

            else
                ! With implicit vertical advection (diffusion + advection)
                uz_aa   = 0.5*(uz(k-1)+uz(k))   ! ac => aa nodes

            end if 

            ! Convert units of Q_strn [J a-1 m-3] => [K a-1]
            Q_strn_now = Q_strn(k)/(rho_ice*cp(k))

            ! Get kappa for the lower and upper ac-nodes 
            ! Note: this is important to avoid mixing of kappa at the 
            ! CTS height (kappa_lower = kappa_temperate; kappa_upper = kappa_cold)
            ! See Blatter and Greve, 2015, Eq. 25. 
            kappa_a = kappa_aa(k)
            kappa_b = kappa_aa(k+1) 

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

        call solve_tridiag(subd,diag,supd,rhs,solution)


        ! == Get variables back in consistent form (enth,T_ice,omega)

        if (use_enth) then 
            ! Copy the solution into the enthalpy variable,
            ! recalculate enthalpy, temperature and water content 
            
            enth  = solution

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
                Q_ice_b = kappa_aa(1) * rho_ice * (enth(2) - enth(1)) / dz
            else
                Q_ice_b = 0.0 
            end if 

            ! Calculate basal mass balance 
            enth_b     = enth(1)
            enth_pmp_b = T_pmp(1) * fac_enth(1)
            call calc_bmb_grounded_enth(bmb_grnd,enth_b,enth_pmp_b,Q_ice_b,Q_b,Q_geo_now,f_grnd,rho_ice)
            
        else 
            ! Copy the solution into the temperature variable,
            ! recalculate enthalpy  

            T_ice = solution 

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
            where (T_ice .ge. T_pmp) omega = 0.01 

            ! Finally, get enthalpy too 
            call convert_to_enthalpy(enth,T_ice,omega,T_pmp,cp,L_ice)

            ! Calculate heat flux at ice base as temperature gradient * conductivity [J a-1 m-2]
            if (H_ice .gt. 0.0_prec) then 
                dz = H_ice * (zeta_aa(2)-zeta_aa(1))
                Q_ice_b = kt(1) * (T_ice(2) - T_ice(1)) / dz 
            else 
                Q_ice_b = 0.0  
            end if 
            
            ! Calculate basal mass balance (valid for grounded ice only)
            call calc_bmb_grounded(bmb_grnd,T_ice(1)-T_pmp(1),Q_ice_b,Q_b,Q_geo_now,f_grnd,rho_ice)
            
        end if 
        
        ! Include internal melting in bmb_grnd 
        bmb_grnd = bmb_grnd - melt_internal 


        ! Finally, calculate the CTS height 
        H_cts = calc_cts_height(enth,T_pmp,cp,H_ice,zeta_aa)

        return 

    end subroutine calc_enth_column

    ! ========== ENTHALPY ==========================================

    elemental subroutine convert_to_enthalpy(enth,T_ice,omega,T_pmp,cp,L_ice)
        ! Given temperature and water content, calculate enthalpy.

        implicit none 

        real(prec), intent(OUT) :: enth             ! [J m-3] Ice enthalpy 
        real(prec), intent(IN)  :: T_ice            ! [K] Ice temperature 
        real(prec), intent(IN)  :: omega            ! [-] Ice water content (fraction)
        real(prec), intent(IN)  :: T_pmp            ! [K] Ice pressure melting point
        real(prec), intent(IN)  :: cp               ! [J kg-1 K-1] Heat capacity 
        real(prec), intent(IN)  :: L_ice            ! [J kg-1] Latent heat of ice 
        
        enth = (1.0_prec-omega)*(cp*T_ice) + omega*(cp*T_pmp + L_ice)

        return 

    end subroutine convert_to_enthalpy

    subroutine convert_from_enthalpy_column(enth,T_ice,omega,T_pmp,cp,L_ice)
        ! Given enthalpy, calculate temperature and water content. 

        implicit none 

        real(prec), intent(INOUT) :: enth(:)            ! [J m-3] Ice enthalpy, nz_aa nodes
        real(prec), intent(OUT)   :: T_ice(:)           ! [K] Ice temperature, nz_aa nodes  
        real(prec), intent(OUT)   :: omega(:)           ! [-] Ice water content (fraction), nz_aa nodes 
        real(prec), intent(IN)    :: T_pmp(:)           ! [K] Ice pressure melting point, nz_aa nodes 
        real(prec), intent(IN)    :: cp(:)              ! [J kg-1 K-1] Heat capacity,nz_aa nodes 
        real(prec), intent(IN)    :: L_ice              ! [J kg-1] Latent heat of ice
        
        ! Local variables
        integer    :: k, nz_aa  
        real(prec), allocatable :: enth_pmp(:)  

        nz_aa = size(enth,1)

        allocate(enth_pmp(nz_aa))

        ! Find pressure melting point enthalpy
        enth_pmp = T_pmp * cp 

        ! Ice interior and basal layer
        ! Note: although the k=1 is a boundary value with no thickness,
        ! allow it to retain omega to maintain consistency with grid points above.
        do k = 1, nz_aa-1

            if (enth(k) .gt. enth_pmp(k)) then
                ! Temperate ice 
                
                T_ice(k) = T_pmp(k)
                omega(k) = (enth(k) - enth_pmp(k)) / L_ice 
             else
                ! Cold ice 

                T_ice(k) = enth(k) / cp(k) 
                omega(k) = 0.0_prec

             end if

        end do 

        ! Surface layer 
        if (enth(nz_aa) .ge. enth_pmp(nz_aa)) then 
            ! Temperate surface, reset omega to zero and enth to pmp value 
            
            enth(nz_aa)  = enth_pmp(nz_aa)
            T_ice(nz_aa) = enth(nz_aa) / cp(nz_aa)
            omega(nz_aa) = 0.0_prec 
        
        else 
            ! Cold surface, calculate T, and reset omega to zero 
            
            T_ice(nz_aa) = enth(nz_aa) / cp(nz_aa)
            omega(nz_aa) = 0.0_prec 
        
        end if 
        
        return 

    end subroutine convert_from_enthalpy_column

    subroutine calc_enth_diffusivity(kappa,T_ice,omega,enth,T_pmp,cp,kt,rho_ice,rho_w,L_ice,cr)
        ! Calculate the enthalpy vertical diffusivity for use with the diffusion solver 
        ! Note: this routine was ported from MALIv6 (Hoffman et al., 2018), 
        ! which appears to have been ported from CISMv2.1. It appears that
        ! using the Harmonic averaging or simple mean gives comparable results.
        !--------------------------------------------------------------------
        ! Comments from MPAS-landice subroutine: 
        !
        ! Compute the enthalpy diffusivity at layer interfaces.
        ! Let d(enth)/dz = the gradient of enthalpy
        ! Can write
        !    d(enth)/dz = d(enth_T)/dz + d(enth_w)/dz,
        ! where
        !    enth_T = (1-phi_w) * rhoi*ci*T
        !    enth_w =    phi_w  * rhow*(L + ci*Tpmp)
        !    phi_w  = water fraction
        !
        ! Now let f = d(enth_T)/dz / d(enth)/dz
        !   (f -> 0 if f is computed to be negative)
        ! For cold ice, f = 1 and diffusivity = diffusivityCold
        ! For temperate ice, f ~ 0 and diffusivity = diffusivityTemperate
        ! At the interface between cold and temperate ice,
        !  f ~ 0 if the temperate ice has large phi_w, but
        !  f ~ 1 if the temperate ice has close to zero phi_w.
        ! Two ways to average:
        ! (1) arithmetic average:  diffusivity = f*diffusivityCold + (1-f)*diffusivityTemperate
        ! (2) harmonic average:    diffusivity = 1 / (f/diffusivityCold + (1-f)/diffusivityTemperate).
        ! Both methods have the same asymptotic values at f = 0 or 1,
        !  but the arithmetic average gives greater diffusivity for
        !  intermediate values.
        !
        ! Still to be determined which is more accurate.
        ! The harmonic average allows large temperature gradients between the
        !  bottom layer and the next layer up; the arithmetic average gives
        !  smoother gradients.
        ! Currently (as of Oct. 2015), the arithmetic average is the default.
        !--------------------------------------------------------------------

        implicit none 

        real(prec), intent(OUT) :: kappa(:)         ! [nz_aa]
        real(prec), intent(IN)  :: T_ice(:)         ! [nz_aa]
        real(prec), intent(IN)  :: omega(:)         ! [nz_aa]
        real(prec), intent(IN)  :: enth(:)          ! [nz_aa]
        real(prec), intent(IN)  :: T_pmp(:)         ! [nz_aa]
        real(prec), intent(IN)  :: cp(:)
        real(prec), intent(IN)  :: kt(:)  
        real(prec), intent(IN)  :: rho_ice
        real(prec), intent(IN)  :: rho_w
        real(prec), intent(IN)  :: L_ice
        real(prec), intent(IN)  :: cr 

        ! Local variables
        integer     :: k, nz_aa  
        real(prec)  :: denth, denth_temp 
        real(prec)  :: f_avg 
        real(prec)  :: kappa_cold       ! Cold diffusivity 
        real(prec)  :: kappa_temp       ! Temperate diffusivity 
        real(prec), allocatable :: enth_temp(:) 
        
        logical, parameter :: use_harmonic_avg = .FALSE. 

        nz_aa = size(enth)

        allocate(enth_temp(nz_aa))

        ! First, define enthalpy associated with temperature only for the whole column 
        ! (for cold ice enth = enth_temp, while for temperate ice enth > enth_temp) 
        !enth_temp = (1.0_prec-omega)*rho_ice*cp*T_ice
        enth_temp = (1.0_prec-omega)*cp*T_ice

        ! Compute factors relating the temperature gradient to the total enthalpy gradient.
        ! Use these factors to average the diffusivity between the cold and temperate kappa values.

        kappa = 0.0 

        do k = 1, nz_aa

            ! Determine kappa_cold and kappa_temp for this level 
            kappa_cold = kt(k) / (rho_ice*cp(k))
            kappa_temp = cr * kappa_cold 

            if (k .lt. nz_aa) then 
                denth      = enth(k+1) - enth(k)
                denth_temp = enth_temp(k+1) - enth_temp(k)   ! = denth in cold ice, < denth in temperate ice
            else 
                denth      = enth(k) - enth(k-1)
                denth_temp = enth_temp(k) - enth_temp(k-1)   ! = denth in cold ice, < denth in temperate ice
            end if 

            if (abs(denth) .lt. 1e-10_prec*L_ice .and. T_ice(k) .lt. T_pmp(k)) then 
                ! Cold ice, no gradient, assign fraction for cold ice diffusion 
                f_avg = 1.0 

            else if (abs(denth) .lt. 1e-10_prec*L_ice) then
                ! Warm ice, no gradient, assign fraction for warm ice diffusion 

                f_avg = 0.0 

            else 
                ! Gradient, determine fraction of diffusion from ice and water content 
                
                f_avg = max(0.0_prec, denth_temp/denth)
                f_avg = min(1.0_prec, f_avg)
            
            end if

            if (use_harmonic_avg) then  
                ! Take a harmonic average
                ! This gives slower cooling of temperate layers and allows
                !  large temperature gradients between cold and temperate layers
                
                kappa(k) = 1.0_prec / ((f_avg/kappa_cold) + (1.0_prec - f_avg)/kappa_temp)
         
            else   
                ! take an arithmetic average
                ! This gives faster cooling of temperate layers and smaller gradients
            
                kappa(k) = f_avg*kappa_cold + (1.0_prec - f_avg)*kappa_temp
                
            end if

        end do

        return 

    end subroutine calc_enth_diffusivity

    function calc_cts_height(enth,T_pmp,cp,H_ice,zeta) result(H_cts)
        ! Calculate the height of the cold-temperate transition surface (m)
        ! within the ice sheet. 

        implicit none 

        real(prec), intent(IN) :: enth(:) 
        real(prec), intent(IN) :: T_pmp(:) 
        real(prec), intent(IN) :: cp(:)
        real(prec), intent(IN) :: H_ice  
        real(prec), intent(IN) :: zeta(:) 
        real(prec) :: H_cts 

        ! Local variables 
        integer :: k, k0, nz 
        real(prec) :: f_lin 
        real(prec), allocatable :: enth_pmp(:) 

        nz = size(enth,1) 

        allocate(enth_pmp(nz))

        ! Get enthalpy at the pressure melting point (no water content)
        enth_pmp = T_pmp * cp

        ! Determine index of the first layer that has enthalpy below enth_pmp 
        do k = 1, nz 
            if (enth(k) .lt. enth_pmp(k)) exit 
        end do 
        
        ! Perform linear interpolation 
        if (k .gt. nz) then 
            ! Whole column is temperate
            H_cts = H_ice
        else if (k .eq. 1) then 
            ! Whole column is cold 
            H_cts = 0.0 
        else 
            ! Perform interpolation
            k0 = k-1 

            ! Get linear weight for where E(f_lin) = Epmp(f_lin)
            ! E(k0) + dE*f_lin = Epmp(k0) + dEpmp*f_lin 
            ! f_lin = (Epmp(k0)-E(k0)) / (dE - dEpmp)
            f_lin = (enth_pmp(k0)-enth(k0)) / ( (enth(k)-enth(k0)) - (enth_pmp(k)-enth_pmp(k0)) )

            H_cts = H_ice * (zeta(k0) + f_lin*(zeta(k)-zeta(k0)))

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

end module ice_enthalpy