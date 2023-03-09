module mass_conservation

    use yelmo_defs, only : sp, dp, wp, TOL_UNDERFLOW, MISSING_VALUE
    use yelmo_tools, only : get_neighbor_indices, fill_borders_2D, set_boundaries_2D_aa

    use solver_advection, only : calc_advec2D  
    use velocity_general, only : set_inactive_margins 
    use topography, only : calc_H_eff

    implicit none 

    private

    public :: calc_G_advec_simple
    public :: calc_G_advec
    public :: calc_G_mbal
    public :: calc_G_calv

    public :: calc_ice_thickness_dyn
    public :: calc_ice_thickness_mbal
    public :: calc_ice_thickness_calving
    public :: apply_ice_thickness_boundaries
    public :: relax_ice_thickness

contains 
    
    subroutine calc_G_advec_simple(G_advec,H_ice,f_ice,ux,uy,mask_adv, &
                                        solver,boundaries,dx,dt,F)
        ! Interface subroutine to update ice thickness through application
        ! of advection, vertical mass balance terms and calving 

        implicit none 

        real(wp),         intent(OUT)   :: G_advec(:,:)         ! [m/yr] Tendency due to advection
        real(wp),         intent(IN)    :: H_ice(:,:)           ! [m]   Ice thickness 
        real(wp),         intent(IN)    :: f_ice(:,:)           ! [--]  Ice area fraction 
        real(wp),         intent(IN)    :: ux(:,:)              ! [m/a] Depth-averaged velocity, x-direction (ac-nodes)
        real(wp),         intent(IN)    :: uy(:,:)              ! [m/a] Depth-averaged velocity, y-direction (ac-nodes)
        integer,          intent(IN)    :: mask_adv(:,:)        ! Advection mask
        character(len=*), intent(IN)    :: solver               ! Solver to use for the ice thickness advection equation
        character(len=*), intent(IN)    :: boundaries
        real(wp),         intent(IN)    :: dx                   ! [m]   Horizontal resolution
        real(wp),         intent(IN)    :: dt                   ! [a]   Timestep 
        real(wp),         intent(IN), optional :: F(:,:) 

        ! Local variables 
        integer :: i, j, nx, ny
        real(wp), allocatable :: F_now(:,:) 
        real(wp), allocatable :: ux_tmp(:,:) 
        real(wp), allocatable :: uy_tmp(:,:) 

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        allocate(F_now(nx,ny))
        allocate(ux_tmp(nx,ny))
        allocate(uy_tmp(nx,ny))

        ! Set local velocity fields with no margin treatment intially
        ux_tmp = ux
        uy_tmp = uy
        
        F_now = 0.0_wp 
        if (present(F)) F_now = F 

        ! Ensure that no velocity is defined for outer boundaries of partially-filled margin points
        call set_inactive_margins(ux_tmp,uy_tmp,f_ice,boundaries)

        ! Determine current advective rate of change (time=n)
        call calc_advec2D(G_advec,H_ice,f_ice,ux_tmp,uy_tmp,F_now,mask_adv,dx,dx,dt,solver,boundaries)

        return 

    end subroutine calc_G_advec_simple

    subroutine calc_G_advec_interp(G_advec,H_ice,f_ice,ux,uy,mask_adv, &
                                        solver,boundaries,dx,dt,mps_to_lo,mps_to_hi)
        ! Interface subroutine to update ice thickness through application
        ! of advection, vertical mass balance terms and calving 

        ! Perform advection at lower resolution. 

        use coordinates_mapping_scrip, only : map_scrip_class, map_scrip_init, map_scrip_field, &
                                            gen_map_filename

        implicit none 

        real(wp),         intent(OUT)   :: G_advec(:,:)         ! [m/yr] Tendency due to advection
        real(wp),         intent(IN)    :: H_ice(:,:)           ! [m]   Ice thickness 
        real(wp),         intent(IN)    :: f_ice(:,:)           ! [--]  Ice area fraction 
        real(wp),         intent(IN)    :: ux(:,:)              ! [m/a] Depth-averaged velocity, x-direction (ac-nodes)
        real(wp),         intent(IN)    :: uy(:,:)              ! [m/a] Depth-averaged velocity, y-direction (ac-nodes)
        integer,          intent(IN)    :: mask_adv(:,:)        ! Advection mask
        character(len=*), intent(IN)    :: solver               ! Solver to use for the ice thickness advection equation
        character(len=*), intent(IN)    :: boundaries
        real(wp),         intent(IN)    :: dx                   ! [m]   Horizontal resolution
        real(wp),         intent(IN)    :: dt                   ! [a]   Timestep 
        type(map_scrip_class), intent(IN) :: mps_to_lo          ! Grid mapping to low resolution
        type(map_scrip_class), intent(IN) :: mps_to_hi          ! Grid mapping to high resolution


        ! Local variables 
        integer :: i, j, nx, ny
        integer :: nx_lo, ny_lo 

        real(wp), allocatable :: F_now(:,:) 
        real(wp), allocatable :: ux_tmp(:,:) 
        real(wp), allocatable :: uy_tmp(:,:) 

        ! Low resolution fields 
        real(wp), allocatable :: H_ice_lo(:,:) 

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        allocate(F_now(nx,ny))
        allocate(ux_tmp(nx,ny))
        allocate(uy_tmp(nx,ny))


        ! First, interpolate input fields to lower resolution using input mapping
        ! Perform conservative interpolation 

        H_ice_lo = MISSING_VALUE 
        call map_scrip_field(mps_to_lo,"H_ice",H_ice,H_ice_lo,method="mean", &
                                    missing_value=MISSING_VALUE,fill_method="nn")

        ! Set local velocity fields with no margin treatment intially
        ux_tmp = ux
        uy_tmp = uy
        
        F_now = 0.0_wp

        ! Ensure that no velocity is defined for outer boundaries of partially-filled margin points
        call set_inactive_margins(ux_tmp,uy_tmp,f_ice,boundaries)

        ! Determine current advective rate of change (time=n)
        call calc_advec2D(G_advec,H_ice,f_ice,ux_tmp,uy_tmp,F_now,mask_adv,dx,dx,dt,solver,boundaries)

        return 

    end subroutine calc_G_advec_interp

    subroutine calc_G_advec(G_adv,dHdt_n,H_ice_n,H_ice_pred,H_ice,f_ice,ux,uy, &
                        mask_pred_new,mask_corr_new,solver,mask_adv,boundaries, &
                        dx,dt,beta,pc_step,F)
        ! Interface subroutine to update ice thickness through application
        ! of advection, vertical mass balance terms and calving 

        implicit none 

        real(wp),         intent(OUT)   :: G_adv(:,:)           ! [m/yr] Tendency due to advection
        real(wp),         intent(INOUT) :: dHdt_n(:,:)          ! [m/a] Advective rate of ice thickness change from previous=>current timestep 
        real(wp),         intent(INOUT) :: H_ice_n(:,:)         ! [m]   Ice thickness from previous=>current timestep 
        real(wp),         intent(IN)    :: H_ice_pred(:,:)      ! [m]   Ice thickness from predicted timestep 
        real(wp),         intent(IN)    :: H_ice(:,:)           ! [m]   Ice thickness 
        real(wp),         intent(IN)    :: f_ice(:,:)           ! [--]  Ice area fraction 
        real(wp),         intent(IN)    :: ux(:,:)              ! [m/a] Depth-averaged velocity, x-direction (ac-nodes)
        real(wp),         intent(IN)    :: uy(:,:)              ! [m/a] Depth-averaged velocity, y-direction (ac-nodes)
        integer,          intent(IN)    :: mask_pred_new(:,:)   
        integer,          intent(IN)    :: mask_corr_new(:,:)  
        integer,          intent(IN)    :: mask_adv(:,:)        ! Advection mask  
        character(len=*), intent(IN)    :: solver               ! Solver to use for the ice thickness advection equation
        character(len=*), intent(IN)    :: boundaries
        real(wp),         intent(IN)    :: dx                   ! [m]   Horizontal resolution
        real(wp),         intent(IN)    :: dt                   ! [a]   Timestep 
        real(wp),         intent(IN)    :: beta(4)              ! Timestep weighting parameters
        character(len=*), intent(IN)    :: pc_step              ! Current predictor-corrector step ('predictor' or 'corrector')
        real(wp),         intent(IN), optional :: F(:,:) 
        
        ! Local variables 
        integer :: i, j, nx, ny
        integer :: im1, ip1, jm1, jp1  
        real(wp), allocatable :: F_now(:,:) 
        real(wp), allocatable :: dHdt_advec(:,:) 
        real(wp), allocatable :: ux_tmp(:,:) 
        real(wp), allocatable :: uy_tmp(:,:) 

        real(wp), parameter :: dHdt_advec_lim = 10.0_wp     ! [m/a] Hard limit on advection rate

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        allocate(F_now(nx,ny))
        allocate(ux_tmp(nx,ny))
        allocate(uy_tmp(nx,ny))
        allocate(dHdt_advec(nx,ny))

        dHdt_advec = 0.0_wp 

        F_now = 0.0_wp 
        if (present(F)) F_now = F 

        ! Set local velocity fields with no margin treatment intially
        ux_tmp = ux
        uy_tmp = uy
        
        ! ===================================================================================
        ! Resolve the dynamic part (ice advection) using multistep method

        select case(trim(pc_step))
        
            case("predictor") 
                
                ! Fill velocity field for new cells 
                call fill_vel_new_cells(ux_tmp,uy_tmp,mask_corr_new,boundaries)

                ! Ensure that no velocity is defined for outer boundaries of partially-filled margin points
                call set_inactive_margins(ux_tmp,uy_tmp,f_ice,boundaries)

                ! Store ice thickness from time=n
                H_ice_n   = H_ice 

                ! Store advective rate of change from saved from previous timestep (now represents time=n-1)
                dHdt_advec = dHdt_n 

                ! Determine current advective rate of change (time=n)
                call calc_advec2D(dHdt_n,H_ice,f_ice,ux_tmp,uy_tmp,F_now,mask_adv,dx,dx,dt,solver,boundaries)

                ! Calculate rate of change using weighted advective rates of change 
                dHdt_advec = beta(1)*dHdt_n + beta(2)*dHdt_advec 
                
                ! Calculate predicted ice thickness (time=n+1,pred)
                !H_ice = H_ice_n + dt*dHdt_advec 

            case("corrector") ! corrector 

                ! Fill velocity field for new cells 
                call fill_vel_new_cells(ux_tmp,uy_tmp,mask_pred_new,boundaries)

                ! Ensure that no velocity is defined for outer boundaries of partially-filled margin points
                call set_inactive_margins(ux_tmp,uy_tmp,f_ice,boundaries)

                ! Determine advective rate of change based on predicted H,ux/y fields (time=n+1,pred)
                call calc_advec2D(dHdt_advec,H_ice_pred,f_ice,ux_tmp,uy_tmp,F_now,mask_adv,dx,dx,dt,solver,boundaries)

                ! Calculate rate of change using weighted advective rates of change 
                dHdt_advec = beta(3)*dHdt_advec + beta(4)*dHdt_n 
                
                ! Calculate corrected ice thickness (time=n+1)
                !H_ice = H_ice_n + dt*dHdt_advec 

                ! Finally, update dHdt_n with correct term to use as n-1 on next iteration
                dHdt_n = dHdt_advec 

        end select
        
        ! Store advective tendency 
        G_adv = dHdt_advec 

        return 

    end subroutine calc_G_advec

    subroutine fill_vel_new_cells(ux,uy,mask,boundaries)

        implicit none

        real(wp), intent(INOUT) :: ux(:,:) 
        real(wp), intent(INOUT) :: uy(:,:) 
        integer,  intent(IN)    :: mask(:,:) 
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: im1, ip1, jm1, jp1 

        nx = size(mask,1)
        ny = size(mask,2) 

        do j = 1, ny
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
            
            if (mask(i,j) .eq. 2) then 
                ! This site just filled with ice, so 
                ! velocity may not be defined on borders
                ! Check borders and fill in empty velocities

                ! x-direction
                if (ux(i,j) .eq. 0.0_wp .and. ux(im1,j) .ne. 0.0_wp) then 
                    ux(i,j) = ux(im1,j) 
                else if (ux(i,j) .ne. 0.0_wp .and. ux(im1,j) .eq. 0.0_wp) then
                    ux(im1,j) = ux(i,j)
                end if 

                ! y-direction
                if (uy(i,j) .eq. 0.0_wp .and. uy(i,jm1) .ne. 0.0_wp) then 
                    uy(i,j) = uy(i,jm1) 
                else if (uy(i,j) .ne. 0.0_wp .and. uy(i,jm1) .eq. 0.0_wp) then
                    uy(i,jm1) = uy(i,j)
                end if 
                
            end if 

        end do 
        end do


        return

    end subroutine fill_vel_new_cells

    subroutine calc_G_mbal(G_mb,H_ice,f_grnd,mbal,dt)
        ! Interface subroutine to update ice thickness through application
        ! of advection, vertical mass balance terms and calving 

        implicit none 

        real(wp), intent(OUT)   :: G_mb(:,:)            ! [m/yr] Actual tendency due to mass balance
        real(wp), intent(IN)    :: H_ice(:,:)           ! [m]   Ice thickness 
        real(wp), intent(IN)    :: f_grnd(:,:)          ! [--]  Grounded fraction 
        real(wp), intent(IN)    :: mbal(:,:)            ! [m/yr] Net mass balance; mbal = smb+bmb+fmb+calv
        real(wp), intent(IN)    :: dt                   ! [a]   Timestep  

        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: im1, ip1, jm1, jp1 

        nx = size(H_ice,1)
        ny = size(H_ice,2) 

        ! ==== MASS BALANCE =====

        ! Initialize G_mb object with diagnosed mass balance everywhere
        G_mb = mbal

        ! Ensure melting is only counted where ice exists 
        where(G_mb .lt. 0.0 .and. H_ice .eq. 0.0) G_mb = 0.0 

        ! Additionally ensure ice cannot form in open ocean 
        where(f_grnd .eq. 0.0 .and. H_ice .eq. 0.0)  G_mb = 0.0  

        ! Ensure melt is limited to amount of available ice to melt  
        where((H_ice+dt*G_mb) .lt. 0.0) G_mb = -H_ice/dt

        return 

    end subroutine calc_G_mbal

    subroutine calc_G_calv(G_calv,H_ice,calv_flt,calv_grnd,dt,calv_flt_method,boundaries)
        ! Interface subroutine to update ice thickness through application
        ! of advection, vertical mass balance terms and calving 

        implicit none 

        real(wp), intent(OUT)   :: G_calv(:,:)          ! [m/yr] Actual calving rate applied to real ice points
        real(wp), intent(IN)    :: H_ice(:,:)           ! [m]   Ice thickness 
        real(wp), intent(IN)    :: calv_flt(:,:)        ! [m/a] Potential calving rate (floating)
        real(wp), intent(IN)    :: calv_grnd(:,:)       ! [m/a] Potential calving rate (grounded)
        real(wp), intent(IN)    :: dt                   ! [a]   Timestep   
        character(len=*), intent(IN) :: calv_flt_method
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: im1, ip1, jm1, jp1 
        logical :: is_margin
        logical :: kill_floating
        real(wp) :: calv_flt_now 
        real(wp) :: calv_grnd_now 


        nx = size(H_ice,1)
        ny = size(H_ice,2) 

        ! Determine whether a kill method is being applied
        kill_floating = .FALSE. 
        if (trim(calv_flt_method) .eq. "kill" .or. &
            trim(calv_flt_method) .eq. "kill-pos") then 

            kill_floating = .TRUE. 

        end if 

        ! ===== CALVING ======

        ! Combine grounded and floating calving into one field for output.
        ! It has already been scaled by area of ice in cell (f_ice).

        ! Note 1: Only allow calving at the current margin 
        ! If ice has retreated before applying calving, then H_ice is 
        ! zero and so G_calv will also be zero. But if ice has advanced,
        ! then calving should also go to zero. 

        ! Note 2: for floating ice, allow calving everywhere if kill_floating is active.

        G_calv = 0.0_wp 
        
        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
            
            is_margin = H_ice(i,j) .gt. 0.0 .and. &
                count([H_ice(im1,j),H_ice(ip1,j),H_ice(i,jm1),H_ice(i,jp1)].eq.0.0) .gt. 0

            if (is_margin .or. kill_floating) then
                calv_flt_now = calv_flt(i,j)
            else
                calv_flt_now = 0.0_wp
            end if 

            if (is_margin) then
                calv_grnd_now = calv_grnd(i,j) 
            else
                calv_grnd_now = 0.0_wp
            end if

            G_calv(i,j) = calv_flt_now + calv_grnd_now
            
            ! Limit calving rate to available ice
            if (H_ice(i,j)-dt*G_calv(i,j) .lt. 0.0) G_calv(i,j) = H_ice(i,j)/dt

        end do 
        end do

        return 

    end subroutine calc_G_calv


    subroutine calc_ice_thickness_dyn(H_ice,dHdt_n,H_ice_n,H_ice_pred,f_ice,ux,uy,mask_adv, &
                                      solver,boundaries,dx,dt,beta,pc_step)
        ! Interface subroutine to update ice thickness through application
        ! of advection, vertical mass balance terms and calving 

        implicit none 

        real(wp),         intent(INOUT) :: H_ice(:,:)           ! [m]   Ice thickness 
        real(wp),         intent(INOUT) :: dHdt_n(:,:)          ! [m/a] Advective rate of ice thickness change from previous=>current timestep 
        real(wp),         intent(INOUT) :: H_ice_n(:,:)         ! [m]   Ice thickness from previous=>current timestep 
        real(wp),         intent(INOUT) :: H_ice_pred(:,:)      ! [m]   Ice thickness from predicted timestep 
        real(wp),         intent(IN)    :: f_ice(:,:)           ! [--]  Ice area fraction 
        real(wp),         intent(IN)    :: ux(:,:)              ! [m/a] Depth-averaged velocity, x-direction (ac-nodes)
        real(wp),         intent(IN)    :: uy(:,:)              ! [m/a] Depth-averaged velocity, y-direction (ac-nodes)
        integer,          intent(IN)    :: mask_adv(:,:)        ! Advection mask
        character(len=*), intent(IN)    :: solver               ! Solver to use for the ice thickness advection equation
        character(len=*), intent(IN)    :: boundaries
        real(wp),         intent(IN)    :: dx                   ! [m]   Horizontal resolution
        real(wp),         intent(IN)    :: dt                   ! [a]   Timestep 
        real(wp),         intent(IN)    :: beta(4)              ! Timestep weighting parameters
        character(len=*), intent(IN)    :: pc_step              ! Current predictor-corrector step ('predictor' or 'corrector')

        ! Local variables 
        integer :: i, j, nx, ny
        integer :: im1, ip1, jm1, jp1  
        real(wp), allocatable :: mbal_zero(:,:) 
        real(wp), allocatable :: dHdt_advec(:,:) 
        real(wp), allocatable :: ux_tmp(:,:) 
        real(wp), allocatable :: uy_tmp(:,:) 

        real(wp), parameter :: dHdt_advec_lim = 10.0_wp     ! [m/a] Hard limit on advection rate

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        allocate(mbal_zero(nx,ny))
        mbal_zero = 0.0_wp 

        allocate(ux_tmp(nx,ny))
        allocate(uy_tmp(nx,ny))
        allocate(dHdt_advec(nx,ny))

        dHdt_advec = 0.0_wp 

        ! Set local velocity fields with no margin treatment intially
        ux_tmp = ux
        uy_tmp = uy
        
        ! Ensure that no velocity is defined for outer boundaries of partially-filled margin points
        call set_inactive_margins(ux_tmp,uy_tmp,f_ice,boundaries)

        ! ===================================================================================
        ! Resolve the dynamic part (ice advection) using multistep method

        select case(trim(pc_step))
        
            case("predictor") 
                
                ! Store ice thickness from time=n
                H_ice_n   = H_ice 

                ! Store advective rate of change from saved from previous timestep (now represents time=n-1)
                dHdt_advec = dHdt_n 

                ! Determine current advective rate of change (time=n)
                call calc_advec2D(dHdt_n,H_ice,f_ice,ux_tmp,uy_tmp,mbal_zero,mask_adv,dx,dx,dt,solver,boundaries)

                ! ajr: testing stability fix for spin-up, limit advection rate!
                ! where(dHdt_n .gt.  dHdt_advec_lim) dHdt_n = dHdt_advec_lim
                ! where(dHdt_n .lt. -dHdt_advec_lim) dHdt_n = -dHdt_advec_lim
                
                ! Calculate rate of change using weighted advective rates of change 
                dHdt_advec = beta(1)*dHdt_n + beta(2)*dHdt_advec 
                
                ! Calculate predicted ice thickness (time=n+1,pred)
                H_ice = H_ice_n + dt*dHdt_advec 

            case("corrector") ! corrector 

                ! Determine advective rate of change based on predicted H,ux/y fields (time=n+1,pred)
                call calc_advec2D(dHdt_advec,H_ice_pred,f_ice,ux_tmp,uy_tmp,mbal_zero,mask_adv,dx,dx,dt,solver,boundaries)

                ! ajr: testing stability fix for spin-up, limit advection rate!
                ! where(dHdt_advec .gt.  dHdt_advec_lim) dHdt_advec = dHdt_advec_lim
                ! where(dHdt_advec .lt. -dHdt_advec_lim) dHdt_advec = -dHdt_advec_lim
                
                ! Calculate rate of change using weighted advective rates of change 
                dHdt_advec = beta(3)*dHdt_advec + beta(4)*dHdt_n 
                
                ! Calculate corrected ice thickness (time=n+1)
                H_ice = H_ice_n + dt*dHdt_advec 

                ! Finally, update dHdt_n with correct term to use as n-1 on next iteration
                dHdt_n = dHdt_advec 

        end select

        ! Post processing of H_ice ================

        ! Also ensure tiny numeric ice thicknesses are removed
        where (H_ice .lt. TOL_UNDERFLOW) H_ice = 0.0 

        return 

    end subroutine calc_ice_thickness_dyn

    subroutine calc_ice_thickness_mbal(H_ice,f_ice,mb_applied,f_grnd,H_ocn,mbal,dx,dt)
        ! Interface subroutine to update ice thickness through application
        ! of advection, vertical mass balance terms and calving 

        implicit none 

        real(wp), intent(INOUT) :: H_ice(:,:)           ! [m]   Ice thickness 
        real(wp), intent(INOUT) :: f_ice(:,:)           ! [--]  Ice cover fraction 
        real(wp), intent(OUT)   :: mb_applied(:,:)      ! [m/a] Actual mass balance applied to real ice points
        real(wp), intent(IN)    :: f_grnd(:,:)          ! [--]  Grounded fraction 
        real(wp), intent(IN)    :: H_ocn(:,:)           ! [m]   Ocean thickness (ie, depth)
        real(wp), intent(IN)    :: mbal(:,:)            ! [m/a] Net mass balance; mbal = smb+bmb (calving separate) 
        real(wp), intent(IN)    :: dx                   ! [m]   Horizontal resolution
        real(wp), intent(IN)    :: dt                   ! [a]   Timestep  

        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: im1, ip1, jm1, jp1 

        nx = size(H_ice,1)
        ny = size(H_ice,2) 

        ! ==== MASS BALANCE =====

        ! ! Combine mass balance and calving into one field, mb_applied.
        ! ! negative mbal is scaled by f_ice, to account for available surface area. 
        ! where(mbal .lt. 0.0_wp)
        !     mb_applied = f_ice*mbal
        ! elsewhere 
        !     mb_applied = mbal
        ! end where

        ! ajr: testing, ensure f_ice scaling is already applied externally 
        mb_applied = mbal

        ! Ensure melting is only counted where ice exists 
        where(mb_applied .lt. 0.0 .and. H_ice .eq. 0.0) mb_applied = 0.0 

        ! Additionally ensure ice cannot form in open ocean 
        where(f_grnd .eq. 0.0 .and. H_ice .eq. 0.0)  mb_applied = 0.0  

        ! Ensure melt is limited to amount of available ice to melt  
        where((H_ice+dt*mb_applied) .lt. 0.0) mb_applied = -H_ice/dt

        ! Apply modified mass balance to update the ice thickness 
        H_ice = H_ice + dt*mb_applied

        ! Ensure tiny numeric ice thicknesses are removed
        where (H_ice .lt. TOL_UNDERFLOW) H_ice = 0.0 

        return 

    end subroutine calc_ice_thickness_mbal

    subroutine calc_ice_thickness_calving(H_ice,f_ice,calv_applied, &
                                       calv_flt,calv_grnd,dx,dt)
        ! Interface subroutine to update ice thickness through application
        ! of advection, vertical mass balance terms and calving 

        implicit none 

        real(wp), intent(INOUT) :: H_ice(:,:)           ! [m]   Ice thickness 
        real(wp), intent(INOUT) :: f_ice(:,:)           ! [--]  Ice cover fraction 
        real(wp), intent(OUT)   :: calv_applied(:,:)    ! [m/a] Actual calving applied to real ice points
        real(wp), intent(IN)    :: calv_flt(:,:)        ! [m/a] Potential calving rate (floating)
        real(wp), intent(IN)    :: calv_grnd(:,:)       ! [m/a] Potential calving rate (grounded)
        real(wp), intent(IN)    :: dx                   ! [m]   Horizontal resolution
        real(wp), intent(IN)    :: dt                   ! [a]   Timestep  

        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: im1, ip1, jm1, jp1 

        nx = size(H_ice,1)
        ny = size(H_ice,2) 

        ! ===== CALVING ======

        ! Combine grounded and floating calving into one field for output.
        ! It has already been scaled by area of ice in cell (f_ice).
        calv_applied = (calv_flt + calv_grnd)

        ! Ensure applied calving is limited to amount of available ice to calve
        where((H_ice-dt*calv_applied) .lt. 0.0) calv_applied = H_ice/dt

        ! Apply modified mass balance to update the ice thickness 
        H_ice = H_ice - dt*calv_applied

        ! Ensure tiny numeric ice thicknesses are removed
        where (H_ice .lt. TOL_UNDERFLOW) H_ice = 0.0 

        return 

    end subroutine calc_ice_thickness_calving

    subroutine apply_ice_thickness_boundaries(mb_resid,H_ice,f_ice,f_grnd,uxy_b,ice_allowed,boundaries,H_ice_ref, &
                                                H_min_flt,H_min_grnd,dt,reset)

        implicit none

        real(wp),           intent(INOUT)   :: mb_resid(:,:)            ! [m/yr] Residual mass balance
        real(wp),           intent(INOUT)   :: H_ice(:,:)               ! [m] Ice thickness 
        real(wp),           intent(IN)      :: f_ice(:,:)               ! [--] Fraction of ice cover
        real(wp),           intent(IN)      :: f_grnd(:,:)              ! [--] Grounded ice fraction
        real(wp),           intent(IN)      :: uxy_b(:,:)               ! [m/a] Basal sliding speed, aa-nodes
        logical,            intent(IN)      :: ice_allowed(:,:)         ! Mask of where ice is allowed to be greater than zero 
        character(len=*),   intent(IN)      :: boundaries               ! Boundary condition choice
        real(wp),           intent(IN)      :: H_ice_ref(:,:)           ! [m]  Reference ice thickness to fill with for boundaries=="fixed"
        real(wp),           intent(IN)      :: H_min_flt                ! [m] Minimum allowed floating ice thickness 
        real(wp),           intent(IN)      :: H_min_grnd               ! [m] Minimum allowed grounded ice thickness 
        real(wp),           intent(IN)      :: dt                       ! [yr] Timestep
        logical, optional,  intent(IN)      :: reset                    ! Reset mb_resid to zero? 

        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: im1, ip1, jm1, jp1 
        real(wp), allocatable :: H_ice_new(:,:)
        real(wp), allocatable :: H_tmp(:,:)
        real(wp) :: H_eff
        real(wp) :: H_max 
        logical  :: is_margin 
        logical  :: is_island 
        logical  :: is_isthmus_x 
        logical  :: is_isthmus_y 

        nx = size(H_ice,1)
        ny = size(H_ice,2) 

        allocate(H_tmp(nx,ny)) 
        allocate(H_ice_new(nx,ny)) 
        H_ice_new = H_ice 

        ! Apply special case for symmetric EISMINT domain when basal sliding is active
        ! (ensure summit thickness does not grow disproportionately)
        if (trim(boundaries) .eq. "EISMINT" .and. maxval(uxy_b) .gt. 0.0) then 
            i = (nx-1)/2 
            j = (ny-1)/2
            H_ice_new(i,j) = (H_ice(i-1,j)+H_ice(i+1,j) &
                                    +H_ice(i,j-1)+H_ice(i,j+1)) / 4.0 
        end if  
        
        ! Artificially delete ice from locations that are not allowed
        ! according to boundary mask (ie, EISMINT, BUELER-A, open ocean)
        where (.not. ice_allowed) H_ice_new = 0.0 
        
        ! Also remove ice that is very small to avoid issues
        where (H_ice_new .lt. 1e-4) H_ice_new = 0.0
        
        ! Remove margin points that are too thin ====

        H_tmp = H_ice_new 

        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
            
            is_margin = H_tmp(i,j) .gt. 0.0 .and. &
                count([H_tmp(im1,j),H_tmp(ip1,j),H_tmp(i,jm1),H_tmp(i,jp1)].eq.0.0) .gt. 0

            if (is_margin) then
                ! Ice covered point at the margin

                ! Calculate current ice thickness 
                call calc_H_eff(H_eff,H_ice_new(i,j),f_ice(i,j)) 

                ! Remove ice that is too thin 
                if (f_grnd(i,j) .eq. 0.0_wp .and. H_eff .lt. H_min_flt)  H_ice_new(i,j) = 0.0_wp 
                if (f_grnd(i,j) .gt. 0.0_wp .and. H_eff .lt. H_min_grnd) H_ice_new(i,j) = 0.0_wp 
 
            end if 

        end do 
        end do 

        ! Remove ice islands =====

        H_tmp = H_ice_new 

        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
            
            ! Check for ice islands
            is_island = H_tmp(i,j) .gt. 0.0 .and. &
                count([H_tmp(im1,j),H_tmp(ip1,j),H_tmp(i,jm1),H_tmp(i,jp1)].gt.0.0) .eq. 0

            if (is_island) then 
                ! Ice-covered island
                ! Remove ice completely. 

                H_ice_new(i,j)   = 0.0_wp 

            end if 

        end do 
        end do 

if (.FALSE.) then
        ! ajr: testing:
        ! This is needed for a bizzare case seen in 
        ! Antarctica, where a grid point with inflow from two neighbors
        ! and positive mass balance, with two ice-free ocean neighbors
        ! nextdoor, grows indefinitely, until the model is killed. A 
        ! solution is needed, but limiting the ice thickness at least 
        ! ensures the simulation continues. 
        where (H_ice_new .gt. 6e3) H_ice_new = 6e3 

        ! ajr: the above situation may be related to poor treatment of
        ! grounded ice-front boundary conditions (in ssa solver). It
        ! can also happen in Greenland at the grounded ice front.

end if 


        ! Reduce ice thickness for margin points that are thicker 
        ! than inland neighbors ====

        H_tmp = H_ice_new

        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
            
            is_margin = H_tmp(i,j) .gt. 0.0 .and. &
                count([H_tmp(im1,j),H_tmp(ip1,j),H_tmp(i,jm1),H_tmp(i,jp1)].eq.0.0) .gt. 0

            if (is_margin) then
                ! Ice covered point at the margin
                
                ! Calculate current ice thickness 
                call calc_H_eff(H_eff,H_ice_new(i,j),f_ice(i,j))

                ! Calculate maximum thickness of neighbors 
                H_max = maxval([H_tmp(im1,j),H_tmp(ip1,j), &
                                H_tmp(i,jm1),H_tmp(i,jp1)])

                if ( H_eff .gt. H_max) then 
                    H_ice_new(i,j) = H_max 
                end if
                
            end if
            
        end do 
        end do

        select case(trim(boundaries))

            case("MISMIP3D","TROUGH")

                ! Do nothing - this should be handled by the ice advection routine
                ! if the default choice ytopo.solver="impl-lis" is used.
                
                !H_ice_new(1,:)    = H_ice_new(2,:)          ! x=0, Symmetry 
                !H_ice_new(nx,:)   = 0.0                     ! x=max, no ice

                !H_ice_new(:,1)    = H_ice_new(:,2)          ! y=-50km, Free-slip condition
                !H_ice_new(:,ny)   = H_ice_new(:,ny-1)       ! y= 50km, Free-slip condition

            case("periodic","periodic-xy")

                ! Do nothing - this should be handled by the ice advection routine
                ! if the default choice ytopo.solver="impl-lis" is used.

            ! case("periodic-x") 

            !     ! Periodic x 
            !     H_ice_new(1:2,:)     = H_ice_new(nx-3:nx-2,:) 
            !     H_ice_new(nx-1:nx,:) = H_ice_new(2:3,:) 
                
            !     ! Infinite (free-slip too)
            !     H_ice_new(:,1)  = H_ice_new(:,2)
            !     H_ice_new(:,ny) = H_ice_new(:,ny-1)

            case("infinite")
                ! Set border points equal to inner neighbors 

                ! Do nothing - this should be handled by the ice advection routine
                ! if the default choice ytopo.solver="impl-lis" is used.

                !call fill_borders_2D(H_ice_new,nfill=1)

            case("fixed") 
                ! Set border points equal to prescribed values from array

                ! Do nothing - this should be handled by the ice advection routine
                ! if the default choice ytopo.solver="impl-lis" is used.

                !call fill_borders_2D(H_ice_new,nfill=1,fill=H_ice_ref)

            case DEFAULT    ! e.g., None/none, zeros, EISMINT
                ! By default, impose zero ice thickness on grid borders

                ! Set border values to zero
                H_ice_new(1,:)  = 0.0
                H_ice_new(nx,:) = 0.0

                H_ice_new(:,1)  = 0.0
                H_ice_new(:,ny) = 0.0
     
        end select
        
        ! Determine mass balance related to changes applied here

        if (reset) then  
            mb_resid = 0.0_wp 
        end if 

        if (dt .ne. 0.0) then 
            mb_resid = mb_resid + (H_ice_new - H_ice) / dt 
        else 
            mb_resid = 0.0
        end if

        ! Reset actual ice thickness to new values 
        H_ice = H_ice_new 

        return

    end subroutine apply_ice_thickness_boundaries

    subroutine relax_ice_thickness(H_ice,f_grnd,mask_grz,H_ref,topo_rel,tau,dt,boundaries)
        ! This routines allows ice within a given mask to be
        ! relaxed to a reference state with certain timescale tau 
        ! (if tau=0), then H_ice = H_ice_ref directly 

        implicit none 

        real(wp), intent(INOUT) :: H_ice(:,:) 
        real(wp), intent(IN)    :: f_grnd(:,:)  
        integer,  intent(IN)    :: mask_grz(:,:) 
        real(wp), intent(IN)    :: H_ref(:,:) 
        integer,  intent(IN)    :: topo_rel 
        real(wp), intent(IN)    :: tau
        real(wp), intent(IN)    :: dt 
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer  :: i, j, nx, ny 
        integer  :: im1, ip1, jm1, jp1
        logical  :: apply_relax 
        real(wp) :: dHdt 

        nx = size(H_ice,1)
        ny = size(H_ice,2) 

        do j = 1, ny
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
            
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
                        (f_grnd(im1,j) .eq. 0.0 .or. f_grnd(ip1,j) .eq. 0.0 &
                        .or. f_grnd(i,jm1) .eq. 0.0 .or. f_grnd(i,jp1) .eq. 0.0)) apply_relax = .TRUE. 
            
                case(3)
                    ! Relax all points
                    
                    apply_relax = .TRUE. 
                
                case(4) 
                    ! Relax all grounded grounding-zone points 

                    if (mask_grz(i,j) .eq. 0 .or. mask_grz(i,j) .eq. 1) then 

                        apply_relax = .TRUE. 

                    end if 

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
