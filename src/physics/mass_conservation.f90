module mass_conservation

    use yelmo_defs, only : sp, dp, wp, TOL_UNDERFLOW, MISSING_VALUE, io_unit_err
    use yelmo_tools, only : get_neighbor_indices, fill_borders_2D, set_boundaries_2D_aa

    use solver_advection, only : calc_advec2D  
    use velocity_general, only : set_inactive_margins 
    use topography, only : calc_H_eff

    implicit none 

    private

    public :: check_mass_conservation
    public :: apply_tendency
    public :: calc_G_advec_simple
    public :: calc_G_advec
    public :: calc_G_mbal
    public :: calc_G_calv
    public :: calc_G_boundaries
    public :: calc_G_relaxation

    public :: extend_floating_slab
    public :: calc_G_remove_fractional_ice
    public :: remove_icebergs
    
contains 
    
    subroutine check_mass_conservation(H_ice,f_ice,f_grnd,dHidt,mb_applied,mb_calv,dHidt_dyn,smb,bmb,fmb,dmb, &
                                                                mb_resid,dx,sec_year,time,dt,units,label)

        implicit none

        real(wp), intent(IN) :: H_ice(:,:)
        real(wp), intent(IN) :: f_ice(:,:)
        real(wp), intent(IN) :: f_grnd(:,:)
        real(wp), intent(IN) :: dHidt(:,:)
        real(wp), intent(IN) :: mb_applied(:,:)
        real(wp), intent(IN) :: mb_calv(:,:)
        real(wp), intent(IN) :: dHidt_dyn(:,:)
        real(wp), intent(IN) :: smb(:,:)
        real(wp), intent(IN) :: bmb(:,:)
        real(wp), intent(IN) :: fmb(:,:)
        real(wp), intent(IN) :: dmb(:,:)
        real(wp), intent(IN) :: mb_resid(:,:)
        real(wp), intent(IN) :: dx 
        real(wp), intent(IN) :: sec_year
        real(wp), intent(IN) :: time
        real(wp), intent(IN) :: dt
        character(len=*), intent(IN) :: units 
        character(len=*), intent(IN) :: label 

        ! Local variables
        integer  :: npts
        real(wp) :: tot_dHidt
        real(wp) :: tot_components
        real(wp) :: tot_mb_applied
        real(wp) :: tot_mb_calv
        real(wp) :: tot_dHidt_dyn
        real(wp) :: conv

        real(wp) :: percent_error

        real(wp), parameter :: tol_mb = 1e-6

        ! Determine conversion factor to units of interest from [m^3/yr]

        select case(trim(units))

            case("m^3/yr")

                conv = 1.0 

            case("km^3/yr")

                conv = 1e-9

            case("Sv")

                conv = 1e-6 / sec_year

            case DEFAULT

                write(io_unit_err,*) "check_mass_conservation:: Error: units not recognized."
                write(io_unit_err,*) "units = ", trim(units)
                stop

        end select

        ! Calculate totals, initially [m^3/yr] => [units]
 
        tot_dHidt       = sum(dHidt)*dx*dx      * conv
        tot_mb_applied  = sum(mb_applied)*dx*dx * conv
        tot_mb_calv     = sum(mb_calv)*dx*dx    * conv
        tot_dHidt_dyn   = sum(dHidt_dyn)*dx*dx  * conv

        ! Get total of components and percent error
        tot_components = tot_mb_applied + tot_mb_calv
        percent_error  = (tot_components - tot_dHidt) / (tot_dHidt+tol_mb) * 100.0 

        write(*,"(a8,a,2f9.3,a3,2g14.4,g10.3,a3,2g13.4,a3,g13.4)") &
                    trim(label), " mbcheck ["//trim(units)//"]: ", time, dt, " | ", &
                    tot_dHidt, tot_components, percent_error, " | ", &
                    tot_mb_applied, tot_mb_calv !, " | ", tot_dHidt_dyn

        return

    end subroutine check_mass_conservation

    subroutine apply_tendency(H_ice,mb_dot,dt,label,adjust_mb)

        implicit none

        real(wp), intent(INOUT) :: H_ice(:,:)
        real(wp), intent(INOUT) :: mb_dot(:,:)
        real(wp), intent(IN)    :: dt 
        character(len=*),  intent(IN) :: label 
        logical, optional, intent(IN) :: adjust_mb

        ! Local variables
        integer :: i, j, nx, ny
        real(wp) :: H_prev
        real(wp) :: dHdt 
        logical  :: allow_adjust_mb
        
        if (dt .gt. 0.0) then 
            ! Only apply this routine if dt > 0!

            allow_adjust_mb = .FALSE.
            if (present(adjust_mb)) allow_adjust_mb = adjust_mb 

            nx = size(H_ice,1)
            ny = size(H_ice,2)

            !!$omp parallel do collapse(2) private(i,j,H_prev,dHdt)
            do j = 1, ny 
            do i = 1, nx 

                ! Store previous ice thickness
                H_prev = H_ice(i,j) 

                ! Now update ice thickness with tendency for this timestep 
                H_ice(i,j) = H_prev + dt*mb_dot(i,j)

                ! Limit ice thickness to zero 
                if (H_ice(i,j) .lt. 0.0) H_ice(i,j) = 0.0 

                ! Ensure tiny numeric ice thicknesses are removed
                if (abs(H_ice(i,j)) .lt. TOL_UNDERFLOW) H_ice(i,j) = 0.0
                
                ! Calculate actual current rate of change
                dHdt = (H_ice(i,j) - H_prev) / dt 

                ! Update mb rate to match ice rate of change perfectly
                if (allow_adjust_mb) then
                    mb_dot(i,j) = dHdt
                end if 

            end do
            end do
            !!$omp end parallel do

        end if 

        return

    end subroutine apply_tendency

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

            ! Calculate calving rate tendency (negative == mass loss)
            G_calv(i,j) = -(calv_flt_now + calv_grnd_now)
            
            ! Limit calving rate to available ice
            if (H_ice(i,j)+dt*G_calv(i,j) .lt. 0.0) G_calv(i,j) = -H_ice(i,j)/dt
            
        end do 
        end do

        return 

    end subroutine calc_G_calv

    subroutine calc_G_boundaries(mb_resid,H_ice,f_ice,f_grnd,uxy_b,ice_allowed,boundaries, &
                                                            H_ice_ref,H_min_flt,H_min_grnd,dt)

        implicit none

        real(wp),           intent(INOUT)   :: mb_resid(:,:)            ! [m/yr] Residual mass balance
        real(wp),           intent(IN)      :: H_ice(:,:)               ! [m] Ice thickness 
        real(wp),           intent(IN)      :: f_ice(:,:)               ! [--] Fraction of ice cover
        real(wp),           intent(IN)      :: f_grnd(:,:)              ! [--] Grounded ice fraction
        real(wp),           intent(IN)      :: uxy_b(:,:)               ! [m/a] Basal sliding speed, aa-nodes
        logical,            intent(IN)      :: ice_allowed(:,:)         ! Mask of where ice is allowed to be greater than zero 
        character(len=*),   intent(IN)      :: boundaries               ! Boundary condition choice
        real(wp),           intent(IN)      :: H_ice_ref(:,:)           ! [m]  Reference ice thickness to fill with for boundaries=="fixed"
        real(wp),           intent(IN)      :: H_min_flt                ! [m] Minimum allowed floating ice thickness 
        real(wp),           intent(IN)      :: H_min_grnd               ! [m] Minimum allowed grounded ice thickness 
        real(wp),           intent(IN)      :: dt                       ! [yr] Timestep

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

        real(wp), parameter :: H_min_tol = 1e-6

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
        
        ! Remove margin points that are too thin, or points that are below tolerance ====

        H_tmp = H_ice_new 

        !!$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,is_margin,H_eff)
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

            ! Also remove very thin ice points (eg 1e-6 m thick - thicker than machine tolerance, but thinner than relevant)
            ! E.g., for a very small timestep of dt=1e-3 and an accumulation rate of 0.1 m/yr, 
            ! after one timestep, H = 1e-3*0.1 = 1e-4 m. 
            if (H_tmp(i,j) .lt. H_min_tol) H_ice_new(i,j) = 0.0_wp 

        end do 
        end do
        !!$omp end parallel do

        ! Remove ice islands =====

        H_tmp = H_ice_new 

        !!$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,is_island)
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
        !!$omp end parallel do

        ! Reduce ice thickness for margin points that are thicker 
        ! than inland neighbors ====

        H_tmp = H_ice_new

        !!$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,H_eff,H_max)
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
        !!$omp end parallel do
        
        select case(trim(boundaries))

            case("MISMIP3D","TROUGH")

                ! Ensure that x-boundary ice-thickness matches "infinite" and "zero" case.
                ! This is not fully handled by ytopo.solver="impl-lis", since
                ! subsequently the mb forcing is applied.

                H_ice_new(1,:)    = H_ice_new(2,:)          ! x=0, Symmetry 
                H_ice_new(nx,:)   = 0.0                     ! x=max, no ice

            case("periodic","periodic-xy")

                ! Do nothing - this should be handled by the ice advection routine
                ! if the default choice ytopo.solver="impl-lis" is used.

            case("infinite")
                ! Set border points equal to inner neighbors 
                ! This is not fully handled by ytopo.solver="impl-lis", since
                ! subsequently the mb forcing is applied.
                
                H_ice_new(1,:)    = H_ice_new(2,:)          ! x=0, Symmetry 
                H_ice_new(nx,:)   = 0.0                     ! x=max, no ice
                H_ice_new(:,1)  = H_ice_new(:,2)
                H_ice_new(:,ny) = H_ice_new(:,ny-1)

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
        
        ! Determine rate of mass balance related to changes applied here

        if (dt .ne. 0.0) then 
            mb_resid = (H_ice_new - H_ice) / dt 
        else 
            mb_resid = 0.0
        end if

        return

    end subroutine calc_G_boundaries

    subroutine calc_G_relaxation(dHdt,H_ice,f_grnd,mask_grz,H_ref,topo_rel,tau,dt,boundaries)
        ! This routines allows ice within a given mask to be
        ! relaxed to a reference state with certain timescale tau 
        ! (if tau=0), then H_ice = H_ice_ref directly 

        implicit none 

        real(wp), intent(OUT)   :: dHdt(:,:) 
        real(wp), intent(IN)    :: H_ice(:,:) 
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

        nx = size(H_ice,1)
        ny = size(H_ice,2) 

        dHdt = 0.0 

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

                    if (dt .gt. 0.0) then
                        dHdt(i,j)  = (H_ref(i,j) - H_ice(i,j)) / dt
                    else
                        dHdt(i,j)  = (H_ref(i,j) - H_ice(i,j)) / 1.0
                    end if

                else
                    ! Apply relaxation to reference state 

                    dHdt = (H_ref(i,j) - H_ice(i,j)) / tau 

                end if 
            end if 

        end do 
        end do 


        return 

    end subroutine calc_G_relaxation

    subroutine extend_floating_slab(H_ice,f_grnd,H_slab,n_ext)
        ! Extend ice field so that there is always 
        ! floating ice next to grounded marine margins
        ! Extended ice should be very thin, will 
        ! be assigned value H_slab. Slab will be extended
        ! n_ext points away from marine margin

        implicit none 

        real(wp), intent(INOUT) :: H_ice(:,:) 
        real(wp), intent(IN)    :: f_grnd(:,:) 
        real(wp), intent(IN)    :: H_slab       ! Typically 1 or 0.1 m. 
        integer,  intent(IN)    :: n_ext        ! Number of points to extend slab
        
        ! Local variables 
        integer :: i, j, nx, ny, iter 
        integer :: im1, ip1, jm1, jp1
        logical :: is_marine 

        logical  :: ms4(4)
        real(wp) :: Hi4(4) 
        real(wp) :: fg4(4)
        
        logical,  allocatable :: mask_slab(:,:)
        real(wp), allocatable :: H_new(:,:) 

        nx = size(H_ice,1) 
        ny = size(H_ice,2) 

        allocate(mask_slab(nx,ny)) 
        allocate(H_new(nx,ny)) 

        mask_slab = .FALSE. 
        H_new     = H_ice 

        do iter = 1, n_ext

            do j = 1, ny 
            do i = 1, nx 

                ! BC: Periodic boundary conditions
                im1 = i-1
                if (im1 == 0) then
                    im1 = nx
                end if
                ip1 = i+1
                if (ip1 == nx+1) then
                    ip1 = 1
                end if

                jm1 = j-1
                if (jm1 == 0) then
                    jm1 = ny
                end if
                jp1 = j+1
                if (jp1 == ny+1) then
                    jp1 = 1
                end if

                if ( f_grnd(i,j) .eq. 0.0 .and. H_ice(i,j) .eq. 0.0 ) then 
                    ! Floating ice-free ocean point
                    
                    ! Get neighbor values in convenient arrays
                    fg4 = [f_grnd(im1,j),f_grnd(ip1,j),f_grnd(i,jm1),f_grnd(i,jp1)]
                    Hi4 = [H_ice(im1,j),H_ice(ip1,j),H_ice(i,jm1),H_ice(i,jp1)]
                    ms4 = [mask_slab(im1,j),mask_slab(ip1,j),mask_slab(i,jm1),mask_slab(i,jp1)]

                    if ( (count(fg4 .gt. 0.0 .and. Hi4 .gt. 0.0) .gt. 0) .or. &
                         (count(ms4) .gt. 0) ) then 
                        ! At least one neighbors is either a grounded point
                        ! or an extended slab point - make this point extended slab.

                        H_new(i,j)     = H_slab 
                        
                    end if

                end if 

            end do 
            end do 

            ! Update H_ice to current array 
            H_ice = H_new 

            ! Update mask_slab
            where(H_ice .eq. H_slab) 
                mask_slab = .TRUE. 
            elsewhere
                mask_slab = .FALSE.
            end where

        end do 

        return

    end subroutine extend_floating_slab

    subroutine calc_G_remove_fractional_ice(mb_diff,H_ice,f_ice,dt)
        ! Eliminate fractional ice covered points that only 
        ! have fractional ice neighbors. 

        implicit none 

        real(wp), intent(OUT) :: mb_diff(:,:) 
        real(wp), intent(IN)  :: H_ice(:,:) 
        real(wp), intent(IN)  :: f_ice(:,:) 
        real(wp), intent(IN)  :: dt 

        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: im1, ip1, jm1, jp1 
        real(wp), allocatable :: H_new(:,:) 

        nx = size(H_ice,1) 
        ny = size(H_ice,2) 

        allocate(H_new(nx,ny)) 

        H_new = H_ice 

        !!$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1)
        do j = 1, ny 
        do i = 1, nx 

            ! BC: Periodic boundary conditions
            im1 = i-1
            if (im1 == 0) then
                im1 = nx
            end if
            ip1 = i+1
            if (ip1 == nx+1) then
                ip1 = 1
            end if

            jm1 = j-1
            if (jm1 == 0) then
                jm1 = ny
            end if
            jp1 = j+1
            if (jp1 == ny+1) then
                jp1 = 1
            end if

            if (f_ice(i,j) .gt. 0.0 .and. f_ice(i,j) .lt. 1.0) then 
                ! Fractional ice-covered point 

                if ( count([f_ice(im1,j),f_ice(ip1,j), &
                        f_ice(i,jm1),f_ice(i,jp1)] .eq. 1.0) .eq. 0) then 
                    ! No fully ice-covered neighbors available.
                    ! Point should be removed. 

                    H_new(i,j) = 0.0_wp 

                end if

            end if

        end do 
        end do
        !!$omp end parallel do

        ! Determine rate of mass balance related to changes applied here

        if (dt .ne. 0.0) then 
            mb_diff = (H_new - H_ice) / dt 
        else 
            mb_diff = 0.0
        end if

        return

    end subroutine calc_G_remove_fractional_ice

    subroutine remove_icebergs(H_ice)

        implicit none 

        real(wp), intent(INOUT) :: H_ice(:,:) 

        return

    end subroutine remove_icebergs


end module mass_conservation
