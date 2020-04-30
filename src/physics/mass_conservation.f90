module mass_conservation

    use yelmo_defs, only : sp, dp, prec, tol_underflow, g, rho_ice, rho_sw  
    use yelmo_tools, only : fill_borders_2D

    use solver_advection, only : calc_advec2D  

    implicit none 

    private
    public :: calc_ice_thickness
    public :: relax_ice_thickness

contains 

    subroutine calc_ice_thickness(H_ice,dHdt_n,H_ice_n,H_ice_pred,H_margin,f_ice,mb_applied,f_grnd,H_ocn,ux,uy,mbal,calv, &
                                    z_bed_sd,dx,dt,solver,boundaries,ice_allowed,H_min,sd_min,sd_max,calv_max,beta,pc_step)
        ! Interface subroutine to update ice thickness through application
        ! of advection, vertical mass balance terms and calving 

        implicit none 

        real(prec),       intent(INOUT) :: H_ice(:,:)           ! [m]   Ice thickness 
        real(prec),       intent(INOUT) :: dHdt_n(:,:)          ! [m/a] Advective rate of ice thickness change from previous=>current timestep 
        real(prec),       intent(INOUT) :: H_ice_n(:,:)         ! [m]   Ice thickness from previous=>current timestep 
        real(prec),       intent(INOUT) :: H_ice_pred(:,:)      ! [m]   Ice thickness from predicted timestep 
        real(prec),       intent(INOUT) :: H_margin(:,:)        ! [m]   Margin ice thickness (assuming full area coverage) 
        real(prec),       intent(INOUT) :: f_ice(:,:)           ! [m]   Ice area fraction 
        real(prec),       intent(OUT)   :: mb_applied(:,:)      ! [m/a] Actual mass balance applied to real ice points
        real(prec),       intent(IN)    :: f_grnd(:,:)          ! [--]  Grounded fraction 
        real(prec),       intent(IN)    :: H_ocn(:,:)           ! [m]   Ocean thickness (ie, depth)
        real(prec),       intent(IN)    :: ux(:,:)              ! [m/a] Depth-averaged velocity, x-direction (ac-nodes)
        real(prec),       intent(IN)    :: uy(:,:)              ! [m/a] Depth-averaged velocity, y-direction (ac-nodes)
        real(prec),       intent(IN)    :: mbal(:,:)            ! [m/a] Net mass balance; mbal = smb+bmb  !-calv 
        real(prec),       intent(IN)    :: calv(:,:)            ! [m/a] Calving rate 
        real(prec),       intent(IN)    :: z_bed_sd(:,:)        ! [m]   Standard deviation of bed topography
        real(prec),       intent(IN)    :: dx                   ! [m]   Horizontal resolution
        real(prec),       intent(IN)    :: dt                   ! [a]   Timestep 
        character(len=*), intent(IN)    :: solver               ! Solver to use for the ice thickness advection equation
        character(len=*), intent(IN)    :: boundaries           ! What boundary conditions should apply?
        logical,          intent(IN)    :: ice_allowed(:,:)     ! Mask of where ice is allowed to be greater than zero 
        real(prec),       intent(IN)    :: H_min                ! [m]   Minimum allowed ice thickness parameter
        real(prec),       intent(IN)    :: sd_min               ! [m]   Minimum stdev(z_bed) parameter
        real(prec),       intent(IN)    :: sd_max               ! [m]   Maximum stdev(z_bed) parameter
        real(prec),       intent(IN)    :: calv_max             ! [m]   Maximum grounded calving rate parameter
        real(prec),       intent(IN)    :: beta(4) 
        character(len=*), intent(IN)    :: pc_step 

        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: n  
        real(prec), allocatable :: calv_grnd(:,:) 
        real(prec), allocatable :: dHdt_advec(:,:) 
        real(prec), allocatable :: ux_tmp(:,:) 
        real(prec), allocatable :: uy_tmp(:,:) 

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        allocate(calv_grnd(nx,ny))
        calv_grnd = 0.0_prec 

        allocate(ux_tmp(nx,ny))
        allocate(uy_tmp(nx,ny))
        ux_tmp = 0.0_prec 
        uy_tmp = 0.0_prec 

        allocate(dHdt_advec(nx,ny))
        dHdt_advec = 0.0_prec 

        ! Ensure that no velocity is defined for outer boundaries of margin points
        ux_tmp = ux 
        do j = 1, ny 
        do i = 1, nx-1 
            if (H_margin(i,j) .gt. 0.0_prec .and. H_ice(i+1,j) .eq. 0.0_prec)    ux_tmp(i,j) = 0.0_prec 
            if (H_ice(i,j)    .eq. 0.0_prec .and. H_margin(i+1,j) .gt. 0.0_prec) ux_tmp(i,j) = 0.0_prec 
        end do 
        end do  

        uy_tmp = uy 
        do j = 1, ny-1 
        do i = 1, nx  
            if (H_margin(i,j) .gt. 0.0_prec .and. H_ice(i,j+1) .eq. 0.0_prec)    uy_tmp(i,j) = 0.0_prec 
            if (H_ice(i,j)    .eq. 0.0_prec .and. H_margin(i,j+1) .gt. 0.0_prec) uy_tmp(i,j) = 0.0_prec 
        end do 
        end do  
    
!         ! No margin treatment 
!         ux_tmp = ux
!         uy_tmp = uy 
        
        ! ===================================================================================
        ! First, only resolve the dynamic part (ice advection) using multistep method

        select case(trim(pc_step))
        
            case("predictor") 
                
                ! Store ice thickness from time=n
                H_ice_n   = H_ice 

                ! Store advective rate of change from saved from previous timestep (now represents time=n-1)
                dHdt_advec = dHdt_n 

                ! Determine current advective rate of change (time=n)
                call calc_advec2D(dHdt_n,H_ice,ux_tmp,uy_tmp,mbal*0.0,dx,dx,dt,solver)

                ! Calculate rate of change using weighted advective rates of change 
                dHdt_advec = beta(1)*dHdt_n + beta(2)*dHdt_advec 
                
                ! Calculate predicted ice thickness (time=n+1,pred)
                H_ice = H_ice_n + dt*dHdt_advec 

            case("corrector") ! corrector 

                ! Determine advective rate of change based on predicted H,ux/y fields (time=n+1,pred)
                call calc_advec2D(dHdt_advec,H_ice_pred,ux_tmp,uy_tmp,mbal*0.0,dx,dx,dt,solver)

                ! Calculate rate of change using weighted advective rates of change 
                dHdt_advec = beta(3)*dHdt_advec + beta(4)*dHdt_n 
                
                ! Calculate corrected ice thickness (time=n+1)
                H_ice = H_ice_n + dt*dHdt_advec 

                ! Finally, update dHdt_n with correct term to use as n-1 on next iteration
                dHdt_n = dHdt_advec 

        end select  

        ! Ensure ice thickness is greater than zero for safety 
        where(H_ice .lt. 0.0_prec) H_ice = 0.0_prec 

        ! ===================================================================================
        

        ! Next, handle mass balance in order to be able to diagnose
        ! precisely how much mass was lost/gained 
        mb_applied = mbal

        ! Ensure ice cannot form in open ocean 
        where(f_grnd .eq. 0.0 .and. H_ice .eq. 0.0)  mb_applied = 0.0  

        ! Ensure melt is limited to amount of available ice to melt  
        where((H_ice+dt*mb_applied) .lt. 0.0) mb_applied = -H_ice/dt

        ! Apply modified mass balance to update the ice thickness 
        H_ice = H_ice + dt*mb_applied
        
        if (trim(boundaries) .eq. "MISMIP3D") then 
            ! Do not use H_margin treatment for MISMIP3D, it is problematic
            ! at the domain boundaries.

                f_ice = 0.0 
                where (H_ice .gt. 0.0) f_ice = 1.0
                H_margin = 0.0 

        else 
            ! Diagnose margin ice and determine f_ice

            call calc_ice_margin(H_ice,H_margin,f_ice,f_grnd)

        end if 

        ! Determine grounded calving from highly variable terrain  
        call calc_calving_rate_grounded(calv_grnd,H_ice,f_grnd,z_bed_sd,sd_min,sd_max,calv_max,dt)

        ! Limit grounded ice thickness to below maximum threshold value
        ! based on shear stress limit 
        ! ajr: not ready yet, experimetal!!
        !call limit_grounded_margin_thickness_stress(H_ice,mb_applied,f_grnd,H_ocn,dt)

        ! Limit grounded ice thickess to above minimum and below inland neighbor at the margin
!         call limit_grounded_margin_thickness(H_ice,mb_applied,f_grnd,H_min,dt) 
        call limit_grounded_margin_thickness_flux(H_ice,mb_applied,f_grnd,mbal,ux_tmp,uy_tmp,dx,dt,H_min)
            
        ! Also ensure tiny numeric ice thicknesses are removed
        where (H_ice .lt. 1e-5) H_ice = 0.0 


        ! Artificially delete ice from locations that are not allowed
        ! according to boundary mask (ie, EISMINT, BUELER-A, open ocean)
        where (.not. ice_allowed) H_ice = 0.0 

        ! Post processing of H_ice ================

        select case(trim(boundaries))

            case("zeros")

                ! Set border values to zero
                H_ice(1,:)  = 0.0
                H_ice(nx,:) = 0.0
                H_ice(:,1)  = 0.0
                H_ice(:,ny) = 0.0

            case("EISMINT")

                ! Set border values to zero
                H_ice(1,:)  = 0.0
                H_ice(nx,:) = 0.0
                H_ice(:,1)  = 0.0
                H_ice(:,ny) = 0.0
                
            case("MISMIP3D")

                ! === MISMIP3D =====
                H_ice(1,:)    = H_ice(2,:)       ! x=0, Symmetry 
                H_ice(nx,:)   = 0.0              ! x=800km, no ice
                H_ice(:,1)    = H_ice(:,2)       ! y=-50km, Free-slip condition
                H_ice(:,ny)   = H_ice(:,ny-1)    ! y= 50km, Free-slip condition

            case("infinite")
                ! ajr: we should check setting border H values equal to inner neighbors
                
                call fill_borders_2D(H_ice,nfill=1)

            case DEFAULT 

                write(*,*) "calc_ice_thickness:: error: boundary method not recognized: "//trim(boundaries)
                stop 

        end select 
        
        return 

    end subroutine calc_ice_thickness

    subroutine limit_grounded_margin_thickness(H_ice,mb_applied,f_grnd,H_min,dt)

        implicit none 

        real(prec), intent(INOUT) :: H_ice(:,:) 
        real(prec), intent(INOUT) :: mb_applied(:,:) 
        real(prec), intent(IN)    :: f_grnd(:,:) 
        real(prec), intent(IN)    :: H_min
        real(prec), intent(IN)    :: dt
        
        ! Local variables 
        integer :: i, j, nx, ny 
        logical :: is_margin 
        real(prec) :: H_min_neighbor, H_diff 

        real(prec), allocatable :: H_ice_0(:,:) 

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        allocate(H_ice_0(nx,ny))
        H_ice_0 = H_ice 

        do j = 2, ny-1 
        do i = 2, nx-1 

            ! Determine if ice-covered point has an ice-free neighbor (ie, at the ice margin)
            is_margin = (H_ice_0(i,j) .gt. 0.0 .and. f_grnd(i,j) .gt. 0.0 &
                .and. minval([H_ice_0(i-1,j),H_ice_0(i+1,j),H_ice_0(i,j-1),H_ice_0(i,j+1)]) .eq. 0.0)

            if (is_margin .and. H_ice_0(i,j) .lt. H_min .and. mb_applied(i,j) .le. 0.0) then 
                ! If grounded margin point is too thin and mass balance not positive, 
                ! impose ablation to set ice thickness to zero

                mb_applied(i,j) = mb_applied(i,j) - H_ice_0(i,j)/dt
                H_ice(i,j)      = 0.0 

            else if (is_margin) then 
                ! Check if margin is thicker than inland neighbors

                H_min_neighbor = minval([H_ice_0(i-1,j),H_ice_0(i+1,j),H_ice_0(i,j-1),H_ice_0(i,j+1)], &
                                    mask=[H_ice_0(i-1,j),H_ice_0(i+1,j),H_ice_0(i,j-1),H_ice_0(i,j+1)] .gt. 0.0)

                if(H_ice_0(i,j) .gt. H_min_neighbor) then 
                    H_diff = H_ice_0(i,j)-0.5*H_min_neighbor
                    mb_applied(i,j) = mb_applied(i,j) - H_diff/dt 
                    H_ice(i,j)      = H_ice_0(i,j)-H_diff
                end if 

            end if
            
        end do 
        end do 

        return 

    end subroutine limit_grounded_margin_thickness

    subroutine limit_grounded_margin_thickness_flux(H_ice,mb_applied,f_grnd,mbal,ux,uy,dx,dt,H_min)
        ! Remove marginal ice that is too thin (H_ice < H_min) with mass balance
        ! and cannot be replenished from upstream flux

        implicit none 

        real(prec), intent(INOUT) :: H_ice(:,:)         ! [m] Ice thickness 
        real(prec), intent(INOUT) :: mb_applied(:,:) 
        real(prec), intent(IN) :: f_grnd(:,:)           ! [-] Grounded fraction
        real(prec), intent(IN) :: mbal(:,:)             ! [m/a] Net mass balance 
        real(prec), intent(IN) :: ux(:,:)               ! [m/a] velocity, x-direction (ac-nodes)
        real(prec), intent(IN) :: uy(:,:)               ! [m/a] velocity, y-direction (ac-nodes)
        real(prec), intent(IN) :: dx, dt 
        real(prec), intent(IN) :: H_min                 ! [m] Threshold for calving

        ! Local variables 
        integer :: i, j, nx, ny
        real(prec) :: eps_xx, eps_yy  
        logical :: test_mij, test_pij, test_imj, test_ipj
        logical :: is_margin, positive_mb 
        real(prec), allocatable :: dHdt(:,:), H_diff(:,:) 
        real(prec), allocatable :: H_ice_0(:,:) 

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        allocate(dHdt(nx,ny))
        allocate(H_diff(nx,ny))
        allocate(H_ice_0(nx,ny))
        
        ! Ice thickness above threshold
        H_diff = H_ice - H_min

        ! Diagnosed lagrangian rate of change
        dHdt = 0.0 

        do j = 2, ny
        do i = 2, nx
        
                ! Calculate strain rate locally (aa-node)
                eps_xx = (ux(i,j) - ux(i-1,j))/dx
                eps_yy = (uy(i,j) - uy(i,j-1))/dx

                ! Calculate thickness change via conservation
                dHdt(i,j) = mbal(i,j) - H_ice(i,j)*(eps_xx+eps_yy)

        end do 
        end do
        
        ! Store original ice thickness field 
        H_ice_0 = H_ice 
        
        do j = 2, ny-1
        do i = 2, nx-1

            ! Determine how many ice-free points are bordering grounded ice point

            ! Determine if grounded, ice-covered point has an ice-free neighbor (ie, at the grounded ice margin)
            is_margin = (H_ice_0(i,j) .gt. 0.0 .and. f_grnd(i,j) .gt. 0.0 &
                .and. minval([H_ice_0(i-1,j),H_ice_0(i+1,j),H_ice_0(i,j-1),H_ice_0(i,j+1)]) .eq. 0.0)

            if (is_margin .and. H_diff(i,j).lt.0.0) then 
                ! Check if current point is at the margin,
                ! and has thickness less than threshold, or if
                ! ice below H_min limit, accounting for mass flux from inland

                positive_mb = (mbal(i,j).gt.0.0)

                test_mij = ( ((H_diff(i-1,j).gt.0.0).and.(ux(i-1,j).gt.0.0)  &  ! neighbor (i-1,j) total > H_calv
                    .and.  (dHdt(i-1,j).gt.(-H_diff(i-1,j)*abs(ux(i-1,j)/dx)))) & 
                    .or.(f_grnd(i-1,j).gt.0.0.and.positive_mb ))

                test_pij = ( ((H_diff(i+1,j).gt.0.0).and.(ux(i,j).lt.0.0) & ! neighbor (i+1,j) total > H_calv
                    .and.(dHdt(i+1,j).gt.(-H_diff(i+1,j)*abs(ux(i,j)/dx)))) &
                    .or.(f_grnd(i+1,j).gt.0.0.and.positive_mb ))

                test_imj = ( ((H_diff(i,j-1).gt.0.0).and.(uy(i,j-1).gt.0.0)  &  ! neighbor (i,j-1) total > H_calv
                    .and.(dHdt(i,j-1).gt.(-H_diff(i,j-1)*abs(uy(i,j-1)/dx))))&
                    .or.(f_grnd(i,j-1).gt.0.0.and.positive_mb ))

                test_ipj = ( ((H_diff(i,j+1).gt.0.0).and.(uy(i,j).lt.0.0) & ! neighbor (i,j+1) total > H_calv
                    .and.(dHdt(i,j+1).gt.(-H_diff(i,j+1)*abs(uy(i,j)/dx))))&
                    .or.(f_grnd(i,j+1).gt.0.0.and.positive_mb ))

                if ((.not.(test_mij.or.test_pij.or.test_imj.or.test_ipj))) then
                    mb_applied(i,j) = mb_applied(i,j) - H_ice_0(i,j) / dt    
                    H_ice(i,j)      = 0.0          
                end if  

            end if

        end do
        end do

        return 

    end subroutine limit_grounded_margin_thickness_flux
    
    subroutine limit_grounded_margin_thickness_stress(H_ice,mb_applied,f_grnd,H_ocn,dt)
        ! Remove marginal ice that exceeds a stress threshold following
        ! Bassis and Walker (2012), Eq. 2.12 

        implicit none 

        real(prec), intent(INOUT) :: H_ice(:,:)             ! [m] Ice thickness 
        real(prec), intent(INOUT) :: mb_applied(:,:)        ! [m/a] Applied mass balance
        real(prec), intent(IN)    :: f_grnd(:,:)            ! [-] Grounded fraction
        real(prec), intent(IN)    :: H_ocn(:,:)             ! [m] Ocean thickness (depth)
        real(prec), intent(IN)    :: dt 

        ! Local variables 
        integer    :: i, j, nx, ny
        real(prec) :: tau_c, H_max, H_ocn_now 
        logical    :: is_margin  
        real(prec) :: rho_ice_g, rho_sw_ice, rho_ice_sw  
        real(prec), allocatable :: H_ice_0(:,:) 

        real(prec), parameter :: C0    = 1e6                ! [Pa] Depth-averaged shear stress in ice 
        real(prec), parameter :: alpha = 0.0                ! [--] Friction coefficient for Bassis and Walker (2012), Eq. 2.13
        real(prec), parameter :: r     = 0.0                ! [--] Crevasse fraction 
        
        rho_ice_g  = rho_ice * g 
        rho_sw_ice = rho_sw / rho_ice 
        rho_ice_sw = rho_ice / rho_sw 

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        allocate(H_ice_0(nx,ny))
        H_ice_0 = H_ice 
        
        do j = 2, ny-1
        do i = 2, nx-1 

            ! Determine if grounded, ice-covered point has an ice-free neighbor (ie, at the grounded ice margin)
            is_margin = (H_ice_0(i,j) .gt. 0.0 .and. f_grnd(i,j) .gt. 0.0 &
                .and. minval([H_ice_0(i-1,j),H_ice_0(i+1,j),H_ice_0(i,j-1),H_ice_0(i,j+1)]) .eq. 0.0)

            if (is_margin) then 
                ! Determine if this margin point should fail 

                ! Calculate depth of seawater (limited by ice thickness and flotation criterion)
                if (H_ocn(i,j) .gt. 0.0) then 
                    H_ocn_now = min(rho_ice_sw*H_ice_0(i,j),H_ocn(i,j))
                else 
                    H_ocn_now = 0.0 
                end if 

                ! Get depth-averaged shear-stress in ice, Bassis and Walker (2012), Eq. 2.13 vertically integrated
                ! alpha = 0.65: model S1 validated for cold ice 
                ! alpha = 0.4 : model S2 for warmer ice 
                ! alpha = 0.0 : model S3 for purely plastic yielding (default)
                tau_c = C0 + 0.5*alpha*rho_ice_g*H_ice_0(i,j)

                ! Get critical ice thickness to cause stress failure
                H_max = (1.0-r)*tau_c/rho_ice_g + sqrt(((1.0-r)*tau_c/rho_ice_g)**2 + rho_sw_ice*H_ocn_now**2)

                if (H_ice_0(i,j) .gt. H_max) then 
                    ! Critical stress exceeded, calve this ice - pass to mb_applied for now...

                    mb_applied(i,j) = mb_applied(i,j) - H_ice_0(i,j) / dt    
                    H_ice(i,j)      = 0.0  

                end if 

            end if 

        end do 
        end do 

        return 

    end subroutine limit_grounded_margin_thickness_stress
    
    subroutine calc_calving_rate_grounded(calv,H_ice,f_grnd,z_bed_sd,sd_min,sd_max,calv_max,dt)
        ! Parameterize grounded ice-margin calving as a function of 
        ! standard deviation of bedrock at each grid point.
        ! Assumes that higher variability in subgrid implies cliffs
        ! that are not represented at low resolution. 

        implicit none 

        real(prec), intent(OUT) :: calv(:,:)                ! [m/a] Calculated calving rate 
        real(prec), intent(IN)  :: H_ice(:,:)               ! [m] Ice thickness 
        real(prec), intent(IN)  :: f_grnd(:,:)              ! [-] Grounded fraction
        real(prec), intent(IN)  :: z_bed_sd(:,:)            ! [m] Standard deviation of bedrock topography
        real(prec), intent(IN)  :: sd_min                   ! [m] stdev(z_bed) at/below which calv=0
        real(prec), intent(IN)  :: sd_max                   ! [m] stdev(z_bed) at/above which calv=calv_max 
        real(prec), intent(IN)  :: calv_max                 ! [m/a] Maximum allowed calving rate
        real(prec), intent(IN)  :: dt      

        ! Local variables
        integer :: i, j, nx, ny  
        real(prec) :: f_scale 
        logical    :: is_grnd_margin 

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        calv = 0.0 
        
        if (calv_max .gt. 0.0) then 
            ! Determine grounded calving rate 

            do j = 2, ny-1
            do i = 2, nx-1 

                ! Determine if grounded, ice-covered point has an ice-free neighbor (ie, at the grounded ice margin)
                is_grnd_margin = (H_ice(i,j) .gt. 0.0 .and. f_grnd(i,j) .gt. 0.0 &
                    .and. minval([H_ice(i-1,j),H_ice(i+1,j),H_ice(i,j-1),H_ice(i,j+1)]) .eq. 0.0)

                if (is_grnd_margin) then
                    ! Grounded ice-covered point

                    f_scale = (z_bed_sd(i,j) - sd_min)/(sd_max-sd_min)
                    if (f_scale .lt. 0.0) f_scale = 0.0 
                    if (f_scale .gt. 1.0) f_scale = 1.0 

                    ! Calculate calving rate from linear function, limited
                    ! to available ice thickness 
                    calv(i,j) = min(f_scale*calv_max, H_ice(i,j)/dt) 
                    
                end if 

            end do 
            end do 

        end if 

        return 

    end subroutine calc_calving_rate_grounded

    subroutine calc_ice_margin(H_ice,H_margin,f_ice,f_grnd)
        ! Determine the area fraction of a grid cell
        ! that is ice-covered. Assume that marginal points
        ! have equal thickness to inland neighbors 

        implicit none 

        real(prec), intent(INOUT) :: H_ice(:,:)             ! [m] Ice thickness on standard grid (aa-nodes)
        real(prec), intent(INOUT) :: H_margin(:,:)          ! [m] Margin ice thickness for partially filled cells, H_margin*1.0 = H_ref*f_ice
        real(prec), intent(INOUT) :: f_ice(:,:)             ! [--] Ice covered fraction (aa-nodes)
        real(prec), intent(IN)    :: f_grnd(:,:)            ! [--] Grounded fraction (aa-nodes)

        ! Local variables 
        integer :: i, j, nx, ny, i1, i2, j1, j2  
        real(prec) :: H_neighb(4)
        logical :: mask_neighb(4)
        real(prec) :: H_ref  
        real(prec), allocatable :: H_ice_0(:,:) 

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        allocate(H_ice_0(nx,ny))

        ! Initially set fraction to one everywhere there is ice 
        ! and zero everywhere there is no ice
        f_ice = 0.0_prec  
        where (H_ice .gt. 0.0) f_ice = 1.0_prec

        ! For ice-covered points with ice-free neighbors (ie, at the floating or grounded margin),
        ! determine the fraction of grid point that should be ice covered. 

        H_ice_0 = H_ice 

        ! Reset H_margin to zero, will be diagnosed below 
        H_margin = 0.0_prec

        do j = 1, ny
        do i = 1, nx 

            i1 = max(i-1,1)
            j1 = max(j-1,1)
            i2 = min(i+1,nx)
            j2 = min(j+1,ny)

            ! Store neighbor heights 
            H_neighb = [H_ice_0(i1,j),H_ice_0(i2,j),H_ice_0(i,j1),H_ice_0(i,j2)]
            
            if (H_ice(i,j) .gt. 0.0 .and. minval(H_neighb) .eq. 0.0 .and. f_grnd(i,j) .eq. 0.0) then 
                ! This point is at the floating ice margin
!             if (H_ice(i,j) .gt. 0.0 .and. minval(H_neighb) .eq. 0.0) then 
!                 ! This point is at the ice margin

                ! Store mask of neighbors with ice 
                mask_neighb = (H_neighb .gt. 0.0)

                if (count(mask_neighb) .gt. 0) then 
                    ! This point has ice-covered neighbors (generally true)

                    ! Determine height to give to partially filled cell
                    if (f_grnd(i,j) .eq. 0.0) then 
                        ! Floating point, set H_ref = minimum of neighbors

                        H_ref = minval(H_neighb,mask=mask_neighb)

                    else 
                        ! Grounded point, set H_ref < H_mean arbitrarily (0.5 works well)
                        ! Note: H_min instead of H_mean seems to work better (tested with EISMINT2 EXPA + sliding)
                        !H_ref = 0.5*sum(H_neighb,mask=mask_neighb) / real(count(mask_neighb),prec)
                        H_ref = 0.5*minval(H_neighb,mask=mask_neighb)
                    end if
                    
                    ! Determine the cell ice fraction
                    ! Note: fraction is determined as a ratio of 
                    ! thicknesses, derived from volume conservation 
                    ! vol = H_ice*dx*dy = H_ref*area_frac 
                    ! f_ice = area_frac / (dx*dy)
                    ! f_ice = H_ice/H_ref 
                    ! Note: H_ref == 0.0 probably won't happen, but keep if-statement 
                    ! for safety 

                    if (H_ref .gt. 0.0) then 
                        f_ice(i,j) = min( H_ice(i,j) / H_ref, 1.0 ) 
                    else 
                        f_ice(i,j) = 1.0 
                    end if 

                else 
                    ! Island point, assume the cell is not full to 
                    ! ensure it is assigned as an H_margin point

                    H_ref = H_ice(i,j) 
                    f_ice(i,j) = 0.1 

                end if 

                ! Now determine if ice should be in buffer (with f_ice < 1.0)
                if (f_ice(i,j) .gt. 0.0 .and. f_ice(i,j) .lt. 1.0) then 
                    ! Ice exists, but does not fill the entire cell,
                    ! define it in H_margin

                    H_margin(i,j) = H_ice(i,j)

                end if 

            end if  

        end do 
        end do 

        return 

    end subroutine calc_ice_margin
    
    subroutine relax_ice_thickness(H_ice,f_grnd,H_ref,topo_rel,tau,dt)
        ! This routines allows ice within a given mask to be
        ! relaxed to a reference state with certain timescale tau 
        ! (if tau=0), then H_ice = H_ice_ref directly 

        implicit none 

        real(prec), intent(INOUT) :: H_ice(:,:) 
        real(prec), intent(IN)    :: f_grnd(:,:)  
        real(prec), intent(IN)    :: H_ref(:,:) 
        integer,    intent(IN)    :: topo_rel 
        real(prec), intent(IN)    :: tau
        real(prec), intent(IN)    :: dt 

        ! Local variables 
        integer    :: i, j, nx, ny 
        logical    :: apply_relax 
        real(prec) :: dHdt 

        nx = size(H_ice,1)
        ny = size(H_ice,2) 

        do j = 2, ny-1 
            do i = 2, nx-1 

                ! No relaxation to start
                apply_relax = .FALSE.

                select case(topo_rel)

                    case(1) 
                        ! Relax the shelf (floating) ice and ice-free points
                    
                        if (f_grnd(i,j) .eq. 0.0 .or. H_ref(i,j) .eq. 0.0) apply_relax = .TRUE. 
                
                    case(2) 
                        ! Relax the shelf (floating) ice and ice-free points
                        ! and the grounding-line ice too
                        
                        if (f_grnd(i,j) .eq. 0.0 .or. H_ref(i,j) .eq. 0.0) apply_relax = .TRUE. 
                        
                        if (f_grnd(i,j) .gt. 0.0 .and. &
                         (f_grnd(i-1,j) .eq. 0.0 .or. f_grnd(i+1,j) .eq. 0.0 &
                            .or. f_grnd(i,j-1) .eq. 0.0 .or. f_grnd(i,j+1) .eq. 0.0)) apply_relax = .TRUE. 
                
                    case(3)
                        ! Relax all points
                        
                        apply_relax = .TRUE. 
                
                    case DEFAULT
                        ! No relaxation

                        apply_relax = .FALSE.

                end select
                

                if (apply_relax) then 

                    if (tau .eq. 0.0) then
                        ! Impose ice thickness 

                        H_ice(i,j) = H_ref(i,j) 

                    else
                        ! Apply relaxation to reference state 

                        dHdt = (H_ref(i,j) - H_ice(i,j)) / tau 

                        H_ice(i,j) = H_ice(i,j) + dHdt*dt 

                    end if 
                end if 

            end do 
        end do 


        return 

    end subroutine relax_ice_thickness 


end module mass_conservation
