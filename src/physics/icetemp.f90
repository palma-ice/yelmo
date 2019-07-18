module icetemp 
    ! Module contains the ice temperature and basal mass balance (grounded) solution

    use yelmo_defs, only : prec, pi, g, sec_year, T0, rho_ice, rho_sw, rho_w, L_ice  
    use solver_tridiagonal, only : solve_tridiag 
    use thermodynamics, only : calc_bmb_grounded, calc_advec_vertical_column, calc_advec_horizontal_column, &
                                calc_T_pmp, calc_T_base_shlf_approx, calc_temp_linear_column

    implicit none
    
    private
    public :: calc_icetemp_3D
    public :: calc_dzeta_terms

contains 

    subroutine calc_icetemp_3D(T_ice,bmb_grnd,dTdz_b,T_pmp,cp,ct,ux,uy,uz,Q_strn,Q_b,Q_geo, &
                            T_srf,H_ice,H_w,H_grnd,f_grnd,zeta_aa,zeta_ac,dzeta_a,dzeta_b,dt,dx)
        ! Solver for thermodynamics of ice 
        ! Note zeta=height, k=1 base, k=nz surface 

        implicit none 

        real(prec), intent(INOUT) :: T_ice(:,:,:)   ! [K] Ice column temperature
        real(prec), intent(INOUT) :: bmb_grnd(:,:)  ! [m a-1] Basal mass balance (melting is negative)
        real(prec), intent(OUT)   :: dTdz_b(:,:)    ! [K m-1] Ice temperature gradient at the base
        real(prec), intent(IN)    :: T_pmp(:,:,:)   ! [K] Pressure melting point temp.
        real(prec), intent(IN)    :: cp(:,:,:)      ! [J kg-1 K-1] Specific heat capacity
        real(prec), intent(IN)    :: ct(:,:,:)      ! [J a-1 m-1 K-1] Heat conductivity 
        real(prec), intent(IN)    :: ux(:,:,:)      ! [m a-1] Horizontal x-velocity 
        real(prec), intent(IN)    :: uy(:,:,:)      ! [m a-1] Horizontal y-velocity 
        real(prec), intent(IN)    :: uz(:,:,:)      ! [m a-1] Vertical velocity 
        real(prec), intent(IN)    :: Q_strn(:,:,:)  ! [K a-1] Internal strain heat production in ice
        real(prec), intent(IN)    :: Q_b(:,:)       ! [J a-1 m-2] Basal frictional heat production 
        real(prec), intent(IN)    :: Q_geo(:,:)     ! [mW m-2] Geothermal heat flux 
        real(prec), intent(IN)    :: T_srf(:,:)     ! [K] Surface temperature 
        real(prec), intent(IN)    :: H_ice(:,:)     ! [m] Ice thickness 
        real(prec), intent(IN)    :: H_w(:,:)       ! [m] Basal water layer thickness 
        real(prec), intent(IN)    :: H_grnd(:,:)    ! [--] Ice thickness above flotation 
        real(prec), intent(IN)    :: f_grnd(:,:)    ! [--] Grounded fraction
        real(prec), intent(IN)    :: zeta_aa(:)     ! [--] Vertical sigma coordinates (zeta==height), aa-nodes
        real(prec), intent(IN)    :: zeta_ac(:)     ! [--] Vertical sigma coordinates (zeta==height), ac-nodes
        real(prec), intent(IN)    :: dzeta_a(:)    ! d Vertical height axis (0:1) 
        real(prec), intent(IN)    :: dzeta_b(:)    ! d Vertical height axis (0:1) 
        real(prec), intent(IN)    :: dt             ! [a] Time step 
        real(prec), intent(IN)    :: dx             ! [a] Horizontal grid step 
        
        ! Local variables
        integer :: i, j, k, nx, ny, nz_aa, nz_ac  
        real(prec), allocatable  :: advecxy(:)   ! [K a-1 m-2] Horizontal heat advection 
        real(prec) :: T_shlf, H_grnd_lim, f_scalar, T_base  
        real(prec) :: H_ice_now 

        real(prec), allocatable :: T_ice_old(:,:,:) 
        real(prec) :: filter0(3,3), filter(3,3) 

        real(prec), parameter :: H_ice_thin = 15.0   ! [m] Threshold to define 'thin' ice

        nx    = size(T_ice,1)
        ny    = size(T_ice,2)
        nz_aa = size(zeta_aa,1)
        nz_ac = size(zeta_ac,1)

        allocate(advecxy(nz_aa))
        allocate(T_ice_old(nx,ny,nz_aa))

        ! Store original ice temperature field here for input to horizontal advection
        ! calculations 
        T_ice_old = T_ice 

        ! Initialize gaussian filter kernel 
        filter0 = gauss_values(dx,dx,sigma=2.0*dx,n=size(filter,1))

        do j = 3, ny-2
        do i = 3, nx-2 
            
            ! For floating points, calculate the approximate marine-shelf temperature 
            ! ajr, later this should come from an external model, and T_shlf would
            ! be the boundary variable directly
            if (f_grnd(i,j) .lt. 1.0) then 

                ! Calculate approximate marine freezing temp, limited to pressure melting point 
                T_shlf = calc_T_base_shlf_approx(H_ice(i,j),T_pmp(i,j,1),H_grnd(i,j))

            else 
                ! Assigned for safety 

                T_shlf   = T_pmp(i,j,1)

            end if 

            if (H_ice(i,j) .le. H_ice_thin) then 
                ! Ice is too thin or zero: prescribe linear temperature profile
                ! between temperate ice at base and surface temperature 
                ! (accounting for floating/grounded nature via T_base)

                if (f_grnd(i,j) .lt. 1.0) then 
                    ! Impose T_shlf for the basal temperature
                    T_base = T_shlf 
                else 
                    ! Impose the pressure melting point of grounded ice 
                    T_base = T_pmp(i,j,1) 
                end if 

                T_ice(i,j,:) = calc_temp_linear_column(T_srf(i,j),T_base,T_pmp(i,j,nz_aa),zeta_aa)

            else 
                ! Thick ice exists, call thermodynamic solver for the column

                ! No filtering of H_ice, take actual value
!                 H_ice_now = H_ice(i,j) 
                
                ! Filter everywhere
!                 filter = filter0
!                 H_ice_now = sum(H_ice(i-1:i+1,j-1:j+1)*filter)

                ! Filter everywhere there is ice 
                ! filter = filter0 
                ! where (H_ice(i-1:i+1,j-1:j+1) .eq. 0.0) filter = 0.0 
                ! filter = filter/sum(filter)
                ! H_ice_now = sum(H_ice(i-1:i+1,j-1:j+1)*filter)
                
                ! Filter at the margin only 
                if (count(H_ice(i-1:i+1,j-1:j+1) .eq. 0.0) .ge. 2) then
                    filter = filter0 
                    where (H_ice(i-1:i+1,j-1:j+1) .eq. 0.0) filter = 0.0 
                    filter = filter/sum(filter)
                    H_ice_now = sum(H_ice(i-1:i+1,j-1:j+1)*filter)
                else 
                    H_ice_now = H_ice(i,j) 
                end if 
                
                ! Pre-calculate the contribution of horizontal advection to column solution
                ! (use unmodified T_ice_old field as input, to avoid mixing with new solution)
                !call calc_advec_horizontal_column(advecxy,T_ice_old,H_ice,ux,uy,dx,i,j)
                call calc_advec_horizontal_column(advecxy,T_ice_old,H_ice,ux,uy,dx,i,j)
                
                call calc_temp_column(T_ice(i,j,:),bmb_grnd(i,j),dTdz_b(i,j),T_pmp(i,j,:),cp(i,j,:),ct(i,j,:), &
                            uz(i,j,:),Q_strn(i,j,:),advecxy,Q_b(i,j),Q_geo(i,j),T_srf(i,j),T_shlf,H_ice_now, &
                            H_w(i,j),f_grnd(i,j),zeta_aa,zeta_ac,dzeta_a,dzeta_b,dt)
                
            end if 

        end do 
        end do 

        ! Fill in borders 
        T_ice(2,:,:)    = T_ice(3,:,:) 
        T_ice(1,:,:)    = T_ice(3,:,:) 
        T_ice(nx-1,:,:) = T_ice(nx-2,:,:) 
        T_ice(nx,:,:)   = T_ice(nx-2,:,:) 
        
        T_ice(:,2,:)    = T_ice(:,3,:) 
        T_ice(:,1,:)    = T_ice(:,3,:) 
        T_ice(:,ny-1,:) = T_ice(:,ny-2,:) 
        T_ice(:,ny,:)   = T_ice(:,ny-2,:) 

!mmr
        print*,'holacol', sum(bmb_grnd), sum (T_ice), sum(H_ice), sum(dTdz_b)
!mmr

        
        return 

    end subroutine calc_icetemp_3D

    subroutine calc_temp_column(T_ice,bmb_grnd,dTdz_b,T_pmp,cp,ct,uz,Q_strn,advecxy,Q_b,Q_geo, &
                    T_srf,T_shlf,H_ice,H_w,f_grnd,zeta_aa,zeta_ac,dzeta_a,dzeta_b,dt)
        ! Thermodynamics solver for a given column of ice 
        ! Note zeta=height, k=1 base, k=nz surface 
        ! Note: nz = number of vertical boundaries (including zeta=0.0 and zeta=1.0), 
        ! temperature is defined for cell centers, plus a value at the surface and the base
        ! so nz_ac = nz_aa - 1 

        ! For notes on implicit form of advection terms, see eg http://farside.ph.utexas.edu/teaching/329/lectures/node90.html
        
        implicit none 

        real(prec), intent(INOUT) :: T_ice(:)     ! nz_aa [K] Ice column temperature
        real(prec), intent(INOUT) :: bmb_grnd     ! [m a-1] Basal mass balance (melting is negative)
        real(prec), intent(OUT)   :: dTdz_b       ! [K m-1] Basal temperature gradient
        real(prec), intent(IN)    :: T_pmp(:)     ! nz_aa [K] Pressure melting point temp.
        real(prec), intent(IN)    :: cp(:)        ! nz_aa [J kg-1 K-1] Specific heat capacity
        real(prec), intent(IN)    :: ct(:)        ! nz_aa [J a-1 m-1 K-1] Heat conductivity 
        real(prec), intent(IN)    :: uz(:)        ! nz_ac [m a-1] Vertical velocity 
        real(prec), intent(IN)    :: Q_strn(:)    ! nz_aa [K a-1] Internal strain heat production in ice
        real(prec), intent(IN)    :: advecxy(:)   ! nz_aa [K a-1] Horizontal heat advection 
        real(prec), intent(IN)    :: Q_b          ! [J a-1 m-2] Basal frictional heat production 
        real(prec), intent(IN)    :: Q_geo        ! [mW m-2] Geothermal heat flux 
        real(prec), intent(IN)    :: T_srf        ! [K] Surface temperature 
        real(prec), intent(IN)    :: T_shlf       ! [K] Marine-shelf interface temperature
        real(prec), intent(IN)    :: H_ice        ! [m] Ice thickness 
        real(prec), intent(IN)    :: H_w          ! [m] Basal water layer thickness 
        real(prec), intent(IN)    :: f_grnd       ! [--] Grounded fraction
        real(prec), intent(IN)    :: zeta_aa(:)  ! nz_aa [--] Vertical sigma coordinates (zeta==height), layer centered aa-nodes
        real(prec), intent(IN)    :: zeta_ac(:)  ! nz_ac [--] Vertical height axis temperature (0:1), layer edges ac-nodes
        real(prec), intent(IN)    :: dzeta_a(:)  ! nz_aa [--] Solver discretization helper variable ak
        real(prec), intent(IN)    :: dzeta_b(:)  ! nz_aa [--] Solver discretization helper variable bk
        real(prec), intent(IN)    :: dt           ! [a] Time step 

        ! Local variables 
        integer :: k, nz_aa, nz_ac
        real(prec) :: Q_geo_now, ghf_conv 
        real(prec) :: Q_strn_now
        real(prec) :: H_w_predicted
        real(prec) :: T_excess
        real(prec) :: melt_internal   

        real(prec), allocatable :: advecz(:)   ! nz_aa, for explicit vertical advection solving
        logical, parameter      :: test_expl_advecz = .FALSE. 

        real(prec), allocatable :: kappa_aa(:)
        real(prec), allocatable :: dkappadz(:)

        real(prec), allocatable :: subd(:)     ! nz_aa 
        real(prec), allocatable :: diag(:)     ! nz_aa  
        real(prec), allocatable :: supd(:)     ! nz_aa 
        real(prec), allocatable :: rhs(:)      ! nz_aa 
        real(prec), allocatable :: solution(:) ! nz_aa
        real(prec) :: fac, fac_a, fac_b, uz_aa, dzeta, dz, dz1, dz2  
        real(prec) :: kappa_a, kappa_b, dza, dzb 

        nz_aa = size(zeta_aa,1)
        nz_ac = size(zeta_ac,1)

        allocate(kappa_aa(nz_aa))
        allocate(dkappadz(nz_aa))
        
        allocate(subd(nz_aa))
        allocate(diag(nz_aa))
        allocate(supd(nz_aa))
        allocate(rhs(nz_aa))
        allocate(solution(nz_aa))

        ! Get geothermal heat flux in proper units 
        Q_geo_now = Q_geo*1e-3*sec_year   ! [mW m-2] => [J m-2 a-1]

        ! Calculate diffusivity on cell centers (aa-nodes)
        kappa_aa = ct / (rho_ice*cp)

        ! Calculate gradient in kappa (centered on aa-nodes)
        dkappadz = 0.0 
        do k = 2, nz_aa-1 
            dkappadz(k) = (kappa_aa(k+1)-kappa_aa(k-1))/(zeta_aa(k+1)-zeta_aa(k-1))
        end do 

        ! Step 1: apply vertical advection (for explicit testing)
        if (test_expl_advecz) then 
            allocate(advecz(nz_aa))
            advecz = 0.0
            call calc_advec_vertical_column(advecz,T_ice,uz,H_ice,zeta_aa)
            T_ice = T_ice - dt*advecz 
        end if 

        ! Step 2: apply vertical implicit diffusion-advection 
        
        ! == Ice base ==

        if (f_grnd .lt. 1.0) then
            ! Floating or partially floating ice - set temperature equal 
            ! to basal temperature at pressure melting point, or marine freezing temp,
            ! or weighted average between the two.

            ! Impose the weighted average of the pressure melting point and the marine freezing temp.
            subd(1) = 0.0_prec
            diag(1) = 1.0_prec
            supd(1) = 0.0_prec
            rhs(1)  = f_grnd*T_pmp(1) + (1.0-f_grnd)*T_shlf

        else 
            ! Grounded ice 

            ! Determine expected basal water thickness [m] for this timestep,
            ! using basal mass balance from previous time step (rough guess)
            H_w_predicted = H_w - (bmb_grnd*(rho_w/rho_ice))*dt 

            ! == Assign grounded basal boundary conditions ==

            if (T_ice(1) .lt. T_pmp(1) .or. H_w_predicted .lt. 0.0_prec) then   
                ! Frozen at bed, or about to become frozen 

                ! maintain balance of heat sources and sinks
                ! (conductive flux, geothermal flux, and basal friction)

                ! Note: basalHeatFlux is generally >= 0, since defined as positive up

                ! calculate dzeta for the bottom layer between the basal boundary and the temperature point above
                dzeta = zeta_aa(2) - zeta_aa(1)

                ! backward Euler flux basal boundary condition
                subd(1) =  0.0_prec
                diag(1) =  1.0_prec
                supd(1) = -1.0_prec
                rhs(1)  = (Q_b + Q_geo_now) * dzeta*H_ice / ct(1)
                
                if (.FALSE.) then 
                    ! Alternative boundary condition approach - works, but needs testing and refinement 

                ! Convert units of Q_strn [J a-1 m-3] => [K a-1]
                Q_strn_now = Q_strn(1)/(rho_ice*cp(1))

                fac      = dt * ct(1) / (rho_ice*cp(1)) / H_ice**2
               
                subd(1) =  0.0_prec
                diag(1) =  1.0_prec  + fac/((zeta_ac(2)-zeta_ac(1))*(zeta_aa(2) - zeta_aa(1))  )     !mmr  1.0_prec
                supd(1) =  -1.0_prec  * fac/((zeta_ac(2)-zeta_ac(1))*(zeta_aa(2) - zeta_aa(1))  )     !mmr -1.0_prec
                rhs(1)  =  T_ice(1) + fac* ((Q_geo_now) * H_ice / ct(1)) * 1.0_prec/( (zeta_ac(2)-zeta_ac(1))) + Q_strn_now*dt   !+ uz(1) * (Q_b_now) *H*dt / (( zeta_aa(2)-zeta_aa(1) ) * ct(1) )    ! mmr (Q_b_now) * dzetaBot*H / ct(1)
                
                end if 

            else 
                ! Temperate at bed 
                ! Hold basal temperature at pressure melting point

                subd(1) = 0.0_prec
                diag(1) = 1.0_prec
                supd(1) = 0.0_prec
                rhs(1)  = T_pmp(1) 

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

            ! Add 'vertical advection' due to the gradient in kappa 
            !uz_aa = uz_aa + dkappadz(k)

            ! Convert units of Q_strn [J a-1 m-3] => [K a-1]
            Q_strn_now = Q_strn(k)/(rho_ice*cp(k))

        if (.FALSE.) then 
            ! Stagger kappa to the lower and upper ac-nodes

            ! ac-node between k-1 and k 
            if (k .eq. 2) then 
                ! Bottom layer, kappa is kappa for now (later with bedrock kappa?)
                kappa_a = kappa_aa(1)
            else 
                ! Weighted average between lower half and upper half of point k-1 to k 
                dz1 = zeta_ac(k-1)-zeta_aa(k-1)
                dz2 = zeta_aa(k)-zeta_ac(k-1)
                kappa_a = (dz1*kappa_aa(k-1) + dz2*kappa_aa(k))/(dz1+dz2)
            end if 

            ! ac-node between k and k+1 

            ! Weighted average between lower half and upper half of point k to k+1
            dz1 = zeta_ac(k+1)-zeta_aa(k)
            dz2 = zeta_aa(k+1)-zeta_ac(k+1)
            kappa_b = (dz1*kappa_aa(k) + dz2*kappa_aa(k+1))/(dz1+dz2)

        else 
            ! ajr: simply use aa-node kappas for now
            kappa_a = kappa_aa(k) 
            kappa_b = kappa_a

        end if 

            ! Vertical distance for centered difference advection scheme
            dz      =  H_ice*(zeta_aa(k+1)-zeta_aa(k-1))
            
            fac_a   = -kappa_a*dzeta_a(k)*dt/H_ice**2
            fac_b   = -kappa_b*dzeta_b(k)*dt/H_ice**2

            subd(k) = fac_a - uz_aa * dt/dz
            supd(k) = fac_b + uz_aa * dt/dz
            diag(k) = 1.0_prec - fac_a - fac_b
            rhs(k)  = T_ice(k) + dt*Q_strn_now - dt*advecxy(k) 

        end do 

        ! == Ice surface ==

        subd(nz_aa) = 0.0_prec
        diag(nz_aa) = 1.0_prec
        supd(nz_aa) = 0.0_prec
        rhs(nz_aa)  = min(T_srf,T0)


        ! == Call solver ==

        call solve_tridiag(subd,diag,supd,rhs,solution)

        ! Copy the solution into the temperature variables
        T_ice  = solution
        
        ! === Treat basal mass balance and high temperatures ===
        
        ! First calculate internal melt (only allow melting, no accretion)
        
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
        k = 1 
        if (T_ice(k) .gt. T_pmp(k)) T_ice(k) = T_pmp(k)


        ! Get temperature gradient at ice base
        if (H_ice .gt. 0.0_prec) then 
            dz = H_ice * (zeta_aa(2)-zeta_aa(1))
            dTdz_b = (T_ice(2) - T_ice(1)) / dz 
        else 
            dTdz_b = 0.0_prec 
        end if 
        
        ! Calculate basal mass balance (valid for grounded ice only)

!mmr
!mmr        print*,'holagrnd in',sum(bmb_grnd)
!mmr

        call calc_bmb_grounded(bmb_grnd,T_ice(1)-T_pmp(1),dTdz_b,ct(1),rho_ice,Q_b,Q_geo_now,f_grnd)

!mmr
!mmr	print*,'holagrnd out',sum(bmb_grnd)
!mmr
        ! Include internal melting in bmb_grnd 
        bmb_grnd = bmb_grnd - melt_internal 

        return 

    end subroutine calc_temp_column

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


! ========== ENTHALPY ==========================================

    elemental subroutine convert_to_enthalpy_column(enth,T_ice,omega,T_pmp,cp_ice,rho_ice,rho_w,L_ice)

        implicit none 

        real(prec), intent(OUT) :: enth 
        real(prec), intent(IN)  :: T_ice 
        real(prec), intent(IN)  :: omega 
        real(prec), intent(IN)  :: T_pmp
        real(prec), intent(IN)  :: cp_ice
        real(prec), intent(IN)  :: rho_ice
        real(prec), intent(IN)  :: rho_w
        real(prec), intent(IN)  :: L_ice
        
        enth = (1.0_prec-omega)*(rho_ice*cp_ice*T_ice) + omega*(rho_w*(cp_ice*T_pmp+L_ice))

        return 

    end subroutine convert_to_enthalpy_column

    subroutine convert_from_enthalpy_column(T_ice,omega,enth,T_pmp,cp_ice,rho_ice,rho_w,L_ice)

        implicit none 

        real(prec), intent(OUT) :: T_ice(:)   ! [K] Ice temperature 
        real(prec), intent(OUT) :: omega(:)   ! Water content 
        real(prec), intent(IN)  :: enth(:)    ! [nz_aa] Enthalpy, aa-nodes 
        real(prec), intent(IN)  :: T_pmp(:)
        real(prec), intent(IN)  :: cp_ice(:)
        real(prec), intent(IN)  :: rho_ice
        real(prec), intent(IN)  :: rho_w
        real(prec), intent(IN)  :: L_ice
        
        ! Local variables
        integer    :: k, nz_aa  
        real(prec), allocatable :: enth_pmp(:)  

        nz_aa = size(enth,1)

        allocate(enth_pmp(nz_aa))

        ! Find pressure melting point enthalpy
        enth_pmp = T_pmp * rho_ice*cp_ice

        ! Basal layer 
        if (enth(1) >= enth_pmp(1)) then   
            ! Temperate ice: reset enthalpy == pressure melting point, 
            ! and set basal temperature == pressure melting point
            ! and water content to zero, since the
            ! basal layer is infinitesimally small

            !enth(1)  = enth_pmp(1)

            T_ice(1) = T_pmp(1)
            omega(1) = 0.0_prec 
            
         else
            ! Cold ice

            T_ice(1) = enth(1) / (rho_ice*cp_ice(1))
            omega(1) = 0.0_prec 

         endif

        ! Ice interior
        do k = 2, nz_aa-1

            if (enth(k) >= enth_pmp(k)) then
                ! Temperate ice 
                
                T_ice(k) = T_pmp(k)
                omega(k) = (enth(k) - enth_pmp(k)) / ((rho_w-rho_ice)*cp_ice(k)*T_pmp(k)+rho_w*L_ice)

             else
                ! Cold ice 

                T_ice(k) = enth(k) / (rho_ice*cp_ice(k))
                omega(k) = 0.0_prec

             end if

        end do 

        ! Surface temperature is prescribed, set omega at surface to zero 
        omega(nz_aa) = 0.0_prec 

        return 

    end subroutine convert_from_enthalpy_column

    subroutine calc_enth_diffusivity(kappa,T_ice,omega,enth,cp_ice, &
                                        rho_ice,rho_w,L_ice,kappa_cold,kappa_temp)
        ! Calculate the enthalpy vertical diffusivity for use
        ! with the thermodynamics solver 

        implicit none 

        !--------------------------------------------------------------------
        ! Compute the enthalpy diffusivity at layer interfaces.
        ! Let d(enth)/dz = the gradient of enthalpy
        ! Can write
        !    d(enth)/dz = d(enth_T)/dz + d(enth_w)/dz,
        ! where
        !    enth_T = (1-phi_w) * rhoi*ci*T
        !    enth_w =    phi_w  * rhow*(L + ci*Tpmp)
        !    phi_w = water fraction
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
        !
        ! At each temperature point, compute the temperature part of the enthalpy.
        ! Note: enthalpyTemp = enthalpy for cold ice, enthalpyTemp < enthalpy for temperate ice

        real(prec), intent(OUT) :: kappa(:)         ! [nz_ac]
        real(prec), intent(IN)  :: T_ice(:)         ! [nz_aa]
        real(prec), intent(IN)  :: omega(:)         ! [nz_aa]
        real(prec), intent(IN)  :: enth(:)          ! [nz_aa]
        real(prec), intent(IN)  :: cp_ice(:)
        real(prec), intent(IN)  :: rho_ice
        real(prec), intent(IN)  :: rho_w
        real(prec), intent(IN)  :: L_ice
        real(prec), intent(IN)  :: kappa_cold       ! Cold diffusivity 
        real(prec), intent(IN)  :: kappa_temp       ! Temperate diffusivity 

        ! Local variables
        integer :: k, nz_aa  
        real(prec) :: denth, denth_temp 
        real(prec) :: f_avg 
        real(prec), allocatable :: enth_temp(:) 

        logical, parameter :: use_harmonic_avg = .FALSE. 

        nz_aa = size(enth)

        allocate(enth_temp(nz_aa))

        ! First, define enthalpy associated with temperature only 
        enth_temp(1)     = enth(1)
        do k = 2, nz_aa-1
            enth_temp(k) = (1.0_prec-omega(k))*rho_ice*cp_ice(k)*T_ice(k)
        end do
        enth_temp(nz_aa) = enth(nz_aa)
      
        ! Compute factors relating the temperature gradient to the total enthalpy gradient.
        ! Use these factors to average the diffusivity between adjacent temperature points.
        ! ac-nodes vertically 
        
        kappa = 0.0 

        do k = 1, nz_aa-1

            denth      = enth(k+1) - enth(k)
            denth_temp = enth_temp(k+1) - enth_temp(k)   ! = denth in cold ice, < denth in temperate ice

            if (abs(denth) > 1.e-10_prec*rho_w*L_ice) then
                f_avg = max(0.0_prec, denth_temp/denth)
                f_avg = min(1.0_prec, f_avg)
            else
                f_avg = 0.0_prec
            endif

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

    function gauss_values(dx,dy,sigma,n) result(filt)
        ! Calculate 2D Gaussian smoothing kernel
        ! https://en.wikipedia.org/wiki/Gaussian_blur

        implicit none 

        real(prec), intent(IN) :: dx 
        real(prec), intent(IN) :: dy 
        real(prec), intent(IN) :: sigma 
        integer,    intent(IN) :: n 
        real(prec) :: filt(n,n) 

        ! Local variables 
        real(prec) :: x, y  
        integer    :: n2, i, j, i1, j1  

        if (mod(n,2) .ne. 1) then 
            write(*,*) "gauss_values:: error: n can only be odd."
            write(*,*) "n = ", n 
        end if 

        n2 = (n-1)/2 

        do j = -n2, n2 
        do i = -n2, n2 
            x = i*dx 
            y = j*dy 

            i1 = i+1+n2 
            j1 = j+1+n2 
            filt(i1,j1) = 1.0/(2.0*pi*sigma**2)*exp(-(x**2+y**2)/(2*sigma**2))

        end do 
        end do 
        
        ! Normalize to ensure sum to 1
        filt = filt / sum(filt)

        return 

    end function gauss_values

end module icetemp


