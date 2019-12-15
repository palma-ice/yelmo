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

    subroutine calc_enth_column(enth,T_ice,omega,bmb_grnd,Q_ice_b,H_cts,T_pmp,cp,kt,advecxy,uz, &
                                Q_strn,Q_b,Q_geo,T_srf,T_shlf,H_ice,H_w,f_grnd,zeta_aa,zeta_ac, &
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
        
        ! Local variables 
        integer    :: k, nz_aa, nz_ac, k_cts
        real(prec) :: Q_geo_now, ghf_conv 
        real(prec) :: Q_strn_now
        real(prec) :: H_w_predicted
        real(prec) :: T_excess
        real(prec) :: melt_internal   
        real(prec) :: enth_b, enth_pmp_b 
        real(prec) :: dedz 
        real(prec) :: omega_excess

        real(prec), allocatable :: fac_enth(:)  ! aa-nodes 
        real(prec), allocatable :: var(:)       ! aa-nodes 
        real(prec), allocatable :: kappa_aa(:)  ! aa-nodes

        real(prec), allocatable :: subd(:)      ! nz_aa 
        real(prec), allocatable :: diag(:)      ! nz_aa  
        real(prec), allocatable :: supd(:)      ! nz_aa 
        real(prec), allocatable :: rhs(:)       ! nz_aa 
        real(prec), allocatable :: solution(:)  ! nz_aa

        real(prec) :: fac, fac_a, fac_b, uz_aa, dzeta, dz, dz1, dz2 
        real(prec) :: kappa_a, kappa_b 
        logical    :: use_enth  

        logical, parameter :: test_expl_advecz = .TRUE. 
        real(prec), allocatable :: advecz(:) 

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

        ! Get geothermal heat flux in proper units 
        Q_geo_now = Q_geo*1e-3*sec_year   ! [mW m-2] => [J m-2 a-1]

        ! Step 0: Calculate diffusivity, set prognostic variable (T_ice or enth),
        ! and corresponding scaling factor (fac_enth)

        call calc_enth_diffusivity(kappa_aa,enth,T_ice,omega,T_pmp,cp,kt,rho_ice,rho_w,L_ice,cr)

        fac_enth = cp               ! To scale to units of [J kg]
        var      = enth             ! [J kg]

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
                rhs(1)  = ((Q_b + Q_geo_now) * dzeta*H_ice / kt(1)) * fac_enth(1)
                
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

        ! Find height of CTS - heighest temperate layer 
        k_cts = 0 
        do k = 1, nz_aa-1 
            if (enth(k) .ge. T_pmp(k)*cp(k)) then
                k_cts = k 
            else 
                exit 
            end if 
        end do

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

!             ! Get slope at k=2 
!             k = 2
!             dedz = (zeta_aa(k)-zeta_aa(k+1))/((zeta_aa(k+1)-zeta_aa(k+2))*(zeta_aa(k)-zeta_aa(k+2)))*enth(k+2) &
!                  + (zeta_aa(k+2)-zeta_aa(k))/((zeta_aa(k)-zeta_aa(k+1))*(zeta_aa(k+1)-zeta_aa(k+2)))*enth(k+1) &
!                  + (2.0*zeta_aa(k)-zeta_aa(k+1)-zeta_aa(k+2))/((zeta_aa(k)-zeta_aa(k+1))*(zeta_aa(k)-zeta_aa(k+2)))*enth(k)

            dedz    = (enth(3)-enth(2))/(zeta_aa(3)-zeta_aa(2))
            !write(*,*) dedz, (enth(3)-enth(2))/(zeta_aa(3)-zeta_aa(2))

!             enth(1) = enth(2) + dedz*(zeta_aa(1)-zeta_aa(2))
            
!             enth(1) = enth(2)
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
        call calc_bmb_grounded_enth(bmb_grnd,Q_ice_b,Q_b,Q_geo_now,f_grnd,rho_ice)
        
        ! Include internal melting in bmb_grnd 
        bmb_grnd = bmb_grnd - melt_internal 

! ======================= Corrector step for cold ice ==========================
if (.FALSE.) then 

        ! Find height of CTS - heighest temperate layer 
        k_cts = 0 
        do k = 1, nz_aa-1 
            if (enth(k) .ge. T_pmp(k)*cp(k)) then
                k_cts = k 
            else 
                exit 
            end if 
        end do

        if (k_cts .ge. 2) then
            ! Temperate ice exists above the base, recalculate cold layers 

            ! Recalculate diffusivity (only relevant for cold points)
            call calc_enth_diffusivity(kappa_aa,enth,T_ice,omega,T_pmp,cp,kt,rho_ice,rho_w,L_ice,cr)

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

    subroutine calc_enth_diffusivity(kappa,enth,T_ice,omega,T_pmp,cp,kt,rho_ice,rho_w,L_ice,cr)
        ! Calculate the enthalpy vertical diffusivity for use with the diffusion solver:
        ! When water is present in the layer, set kappa=kappa_therm, else kappa=kappa_cold 

        implicit none 

        real(prec), intent(OUT) :: kappa(:)         ! [nz_aa]
        real(prec), intent(IN)  :: enth(:)          ! [nz_aa]
        real(prec), intent(IN)  :: T_ice(:)         ! [nz_aa]
        real(prec), intent(IN)  :: omega(:)         ! [nz_aa]
        real(prec), intent(IN)  :: T_pmp(:)         ! [nz_aa]
        real(prec), intent(IN)  :: cp(:)
        real(prec), intent(IN)  :: kt(:)  
        real(prec), intent(IN)  :: rho_ice
        real(prec), intent(IN)  :: rho_w
        real(prec), intent(IN)  :: L_ice
        real(prec), intent(IN)  :: cr 

        ! Local variables
        integer     :: k, nz_aa   
        real(prec)  :: enth_pmp
        real(prec)  :: kappa_cold       ! Cold diffusivity 
        real(prec)  :: kappa_temp       ! Temperate diffusivity 
        
        nz_aa = size(enth)

        kappa = 0.0 

        do k = 1, nz_aa

            ! Determine kappa_cold and kappa_temp for this level 
            kappa_cold = kt(k) / (rho_ice*cp(k))
            kappa_temp = cr * kappa_cold 

            enth_pmp = T_pmp(k)*cp(k)

            if (omega(k) .gt. 0.0) then 
!             if (enth(k) .ge. enth_pmp) then
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
        real(prec) :: f_lin 
        real(prec), allocatable :: enth_pmp(:) 

        nz = size(enth,1) 

        allocate(enth_pmp(nz))

        ! Get enthalpy at the pressure melting point (no water content)
        enth_pmp = T_pmp * cp

        ! Determine height of CTS as heighest temperate layer 
        k_cts = 0 
        do k = 1, nz 
            !if (enth(k) .ge. T_pmp(k)*cp(k)) then
            if (T_ice(k) .ge. T_pmp(k)) then 
                k_cts = k 
            else 
                exit 
            end if 
        end do 

        if (k_cts .eq. 0) then 
            ! No temperate ice 
            H_cts = 0.0_prec 

        else if (k_cts .eq. nz) then 
            ! Whole column is temperate
            H_cts = H_ice

        else 
            ! Perform linear interpolation 

!             if (T_ice(k_cts)-T_pmp(k_cts) .lt. 0.0_prec) then 
!                 write(*,*) "CTS: ", enth(k_cts)/cp(k_cts), T_ice(k_cts)-T_pmp(k_cts), H_cts
!                 stop 
!             end if 

            H_cts = H_ice * zeta(k_cts) 

        end if 

!         ! Determine index of the first layer that has enthalpy below enth_pmp 
!         do k = 1, nz 
!             if (enth(k) .lt. enth_pmp(k)) exit 
!         end do 
        
!         ! Perform linear interpolation 
!         if (k .gt. nz) then 
!             ! Whole column is temperate
!             H_cts = H_ice
!         else if (k .le. 1) then 
!             ! Whole column is cold 
!             H_cts = 0.0 
!         else 
!             ! Perform interpolation
!             k_cts = k-1 

!             ! Get linear weight for where E(f_lin) = Epmp(f_lin)
!             ! E(k_cts) + dE*f_lin = Epmp(k_cts) + dEpmp*f_lin 
!             ! f_lin = (Epmp(k_cts)-E(k_cts)) / (dE - dEpmp)
!             f_lin = (enth_pmp(k_cts)-enth(k_cts)) / ( (enth(k)-enth(k_cts)) - (enth_pmp(k)-enth_pmp(k_cts)) )
!             if (f_lin .lt. 1e-2) f_lin = 0.0 

!             H_cts = H_ice * (zeta(k_cts) + f_lin*(zeta(k)-zeta(k_cts)))

!             !H_cts = H_ice * zeta(k_cts)

!         end if 

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

end module ice_enthalpy


