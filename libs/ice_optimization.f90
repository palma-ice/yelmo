module ice_optimization

    use yelmo_defs, only : sp, dp, wp, prec, io_unit_err, pi, missing_value, mv, tol_underflow
    use nml 

    use gaussian_filter 
    
    implicit none 

    type ice_opt_params
        logical  :: opt_cf 
        real(wp) :: cf_time_init
        real(wp) :: cf_time_end
        real(wp) :: cf_init
        real(wp) :: cf_min_par
        real(wp) :: tau_c 
        real(wp) :: H0
        logical  :: scaleH
        real(wp) :: sigma_err 
        real(wp) :: sigma_vel 
        character(len=56) :: fill_method 
        logical  :: basin_fill

        real(wp) :: rel_tau 
        real(wp) :: rel_tau1 
        real(wp) :: rel_tau2
        real(wp) :: rel_time1
        real(wp) :: rel_time2
        real(wp) :: rel_m

        logical  :: opt_tf 
        real(wp) :: tf_time_init
        real(wp) :: tf_time_end
        real(wp) :: H_grnd_lim
        real(wp) :: tf_sigma 
        real(wp) :: tau_m 
        real(wp) :: m_temp
        real(wp) :: tf_min 
        real(wp) :: tf_max
        integer  :: tf_basins(100) 

        real(wp), allocatable :: cf_min(:,:) 
        real(wp), allocatable :: cf_max(:,:) 
        
    end type 

    private 

    public :: ice_opt_params
    public :: optimize_par_load 
    public :: optimize_set_transient_param

    public :: optimize_tf_corr
    public :: optimize_tf_corr_basin
    public :: optimize_cb_ref
    public :: optimize_cb_ref_vel
    public :: optimize_cb_ref_pc12
    
    public :: update_mb_corr
    public :: guess_cb_ref
    public :: fill_nearest
    public :: fill_cb_ref
    public :: wtd_mean

    ! Obsolete:
    !public :: update_cb_ref_errscaling
    !public :: update_cb_ref_thickness_ratio

contains 
    
    subroutine optimize_par_load(opt,path_par,group)

        implicit none

        type(ice_opt_params), intent(INOUT) :: opt 
        character(len=*),     intent(IN)    :: path_par 
        character(len=*),     intent(IN)    :: group 

        ! Load optimization parameters 

        call nml_read(path_par,group,"opt_cf",      opt%opt_cf)
        call nml_read(path_par,group,"cf_time_init",opt%cf_time_init)
        call nml_read(path_par,group,"cf_time_end", opt%cf_time_end)
        call nml_read(path_par,group,"cf_init",     opt%cf_init)
        call nml_read(path_par,group,"cf_min",      opt%cf_min_par)
        call nml_read(path_par,group,"tau_c",       opt%tau_c)
        call nml_read(path_par,group,"H0",          opt%H0)
        call nml_read(path_par,group,"scaleH",      opt%scaleH)    
        call nml_read(path_par,group,"sigma_err",   opt%sigma_err)   
        call nml_read(path_par,group,"sigma_vel",   opt%sigma_vel)   
        call nml_read(path_par,group,"fill_method", opt%fill_method)
        call nml_read(path_par,group,"basin_fill",  opt%basin_fill)   
        
        call nml_read(path_par,group,"rel_tau1",    opt%rel_tau1)   
        call nml_read(path_par,group,"rel_tau2",    opt%rel_tau2)  
        call nml_read(path_par,group,"rel_time1",   opt%rel_time1)    
        call nml_read(path_par,group,"rel_time2",   opt%rel_time2) 
        call nml_read(path_par,group,"rel_m",       opt%rel_m)

        call nml_read(path_par,group,"opt_tf",      opt%opt_tf)
        call nml_read(path_par,group,"tf_time_init",opt%tf_time_init)
        call nml_read(path_par,group,"tf_time_end", opt%tf_time_end)
        call nml_read(path_par,group,"H_grnd_lim",  opt%H_grnd_lim)
        call nml_read(path_par,group,"tf_sigma",    opt%tf_sigma)
        call nml_read(path_par,group,"tau_m",       opt%tau_m)
        call nml_read(path_par,group,"m_temp",      opt%m_temp)
        call nml_read(path_par,group,"tf_min",      opt%tf_min)
        call nml_read(path_par,group,"tf_max",      opt%tf_max)
        call nml_read(path_par,group,"tf_basins",   opt%tf_basins)
        
        return

    end subroutine optimize_par_load

    subroutine optimize_set_transient_param(p,time,time1,time2,p1,p2,m)
        ! Determine value of parameter as a function of time 

        implicit none 

        real(wp), intent(OUT) :: p 
        real(wp), intent(IN)  :: time 
        real(wp), intent(IN)  :: time1 
        real(wp), intent(IN)  :: time2
        real(wp), intent(IN)  :: p1
        real(wp), intent(IN)  :: p2
        real(wp), intent(IN)  :: m         ! Non-linear exponent (m=1.0 or higher)
        

        if (time .le. time1) then 
            p = p1 
        else if (time .ge. time2) then 
            p = p2 
        else 
            ! Linear interpolation with non-linear factor m if desired
            p = p1 + (p2-p1)* ((time-time1)/(time2-time1))**m 
        end if  

        return 

    end subroutine optimize_set_transient_param

    subroutine optimize_tf_corr(tf_corr,H_ice,H_grnd,dHicedt,H_obs,H_grnd_obs,H_grnd_lim, &
                                basins,basin_fill,tau_m,m_temp,tf_min,tf_max,dx,sigma,dt)

        implicit none 

        real(wp), intent(INOUT) :: tf_corr(:,:) 
        real(wp), intent(IN)    :: H_ice(:,:)
        real(wp), intent(IN)    :: H_grnd(:,:)
        real(wp), intent(IN)    :: dHicedt(:,:)
        real(wp), intent(IN)    :: H_obs(:,:)
        real(wp), intent(IN)    :: H_grnd_obs(:,:)
        real(wp), intent(IN)    :: H_grnd_lim
        real(wp), intent(IN)    :: basins(:,:) 
        logical,  intent(IN)    :: basin_fill
        real(wp), intent(IN)    :: tau_m 
        real(wp), intent(IN)    :: m_temp
        real(wp), intent(IN)    :: tf_min 
        real(wp), intent(IN)    :: tf_max 
        real(wp), intent(IN)    :: dx 
        real(wp), intent(IN)    :: sigma
        real(wp), intent(IN)    :: dt 

        ! Local variables
        integer  :: b, nb
        integer  :: i, j, nx, ny, n 
        real(wp) :: f_damp 
        real(wp) :: tau_tgt
        real(wp) :: tf_corr_dot
        real(wp) :: tf_corr_tgt
        real(wp) :: tf_corr_bar 

        real(wp), allocatable :: H_err(:,:)
       
        logical,  allocatable :: mask(:,:)
        real(wp), allocatable :: basin_list(:)

        real(wp), parameter :: tol = 1e-5_wp

        ! Internal parameters 
        f_damp = 2.0 

        tau_tgt     = 500.0     ! [yr] Target relaxation timescale 
        tf_corr_tgt = 0.0       ! [degC] Target is no thermal forcing correction 

        nx = size(tf_corr,1)
        ny = size(tf_corr,2) 

        allocate(H_err(nx,ny))
        allocate(mask(nx,ny))

        ! Calculate ice thickness error everywhere
        H_err = H_ice - H_obs

        ! Additionally, apply a Gaussian filter to H_err to ensure smooth transitions
        ! Apply a weighted average between smoothed and original H_err, where 
        ! slow regions get more smoothed, and fast regions use more local error 
        if (sigma .gt. 0.0) then  
            call filter_gaussian(H_err,sigma,dx)
        end if

        ! Limit H_err to desired region, mainly for mainly floating points
        where (H_grnd .gt. H_grnd_lim)
            H_err = 0.0 
        end where 

        do j = 1, ny 
        do i = 1, nx 

            ! Get adjustment rate given error in ice thickness  =========

            tf_corr_dot = 1.0_wp/(tau_m*m_temp) *( (H_err(i,j) / tau_m) + f_damp*dHicedt(i,j) ) &
                    - (1.0_wp/tau_tgt)*(tf_corr(i,j)-tf_corr_tgt)

            ! Apply correction to all points =========

            tf_corr(i,j) = tf_corr(i,j) + tf_corr_dot*dt 

            ! Ensure tf_corr is not below lower or upper limit 
            if (tf_corr(i,j) .lt. tf_min) tf_corr(i,j) = tf_min 
            if (tf_corr(i,j) .gt. tf_max) tf_corr(i,j) = tf_max 

        end do 
        end do

        if (basin_fill) then
                ! Obtain the mean tf_corr value by basins and extrapolate to the basins
                ! Determine unique basin numbers that are available
                nb = MAXVAL(basins)
                allocate(basin_list(nb))
                basin_list = [(i, i=1,nb)]
                
                do b = 1, nb
                        ! Get a mask of points of interest:
                        ! 1. Points within the current basin 
                        ! 2. Points with overburden thickness near flotation,
                        !    with magnitude less than H_grnd_lim
                        ! 3. Points with observed or modeled ice thickness
                        mask =  (abs(basins-basin_list(b)) .lt. tol) .and. &
                                (H_grnd .lt. H_grnd_lim) .and. &
                                (H_obs .gt. 0.0 .or. H_ice .gt. 0.0)

                        ! How many points available 
                        n = count(mask)

                        if (n .gt. 0) then
                                ! Points are available for averaging 
                                ! Calculate average observed thickness for masked region
                                tf_corr_bar = sum(tf_corr,mask=mask)   / real(n,wp)
                        else
                                tf_corr_bar = 0.0_wp
                        end if
                        
                        where ((H_grnd .gt. H_grnd_lim) .and. (abs(basins-basin_list(b)) .lt. tol))
                                tf_corr = tf_corr_bar 
                        end where
                end do 
        else
                ! ! Finally reset tf_corr to zero where no floating ice is observed
                ! ! (ie outside the observed floating ice margin)
                ! where (H_grnd .lt. 0.0 .and. H_obs .eq. 0.0)
                !     tf_corr = 0.0
                ! end where 

                ! ! Also reset to zero for fully grounded ice 
                ! where (H_grnd .gt. 0.0)
                !     tf_corr = 0.0 
                ! end where 
        end if

        return 

    end subroutine optimize_tf_corr

    subroutine optimize_tf_corr_basin(tf_corr,H_ice,H_grnd,dHicedt,H_obs,basins,H_grnd_lim, &
                                    tau_m,m_temp,tf_min,tf_max,tf_basins,dt)

        implicit none 

        real(wp), intent(INOUT) :: tf_corr(:,:) 
        real(wp), intent(IN)    :: H_ice(:,:)
        real(wp), intent(IN)    :: H_grnd(:,:)
        real(wp), intent(IN)    :: dHicedt(:,:)
        real(wp), intent(IN)    :: H_obs(:,:)
        real(wp), intent(IN)    :: basins(:,:) 
        real(wp), intent(IN)    :: H_grnd_lim 
        real(wp), intent(IN)    :: tau_m 
        real(wp), intent(IN)    :: m_temp
        real(wp), intent(IN)    :: tf_min 
        real(wp), intent(IN)    :: tf_max 
        integer,  intent(IN)    :: tf_basins(:) 
        real(wp), intent(IN)    :: dt 

        ! Local variables
        integer  :: i, j, nx, ny, n  
        integer  :: b, nb  
        real(wp) :: f_damp 
        real(wp) :: H_obs_bar 
        real(wp) :: H_bar 
        real(wp) :: H_bar_err
        real(wp) :: dHdt_bar
        real(wp) :: tf_corr_dot

        logical,  allocatable :: mask(:,:) 
        real(wp), allocatable :: basin_list_ref(:)
        real(wp), allocatable :: basin_list(:) 
        
        real(wp), parameter :: tol = 1e-5_wp

        ! Internal parameters 
        f_damp = 2.0 

        nx = size(tf_corr,1)
        ny = size(tf_corr,2) 

        allocate(mask(nx,ny)) 

        ! Determine unique basin numbers that are available
        call unique(basin_list_ref,reshape(basins,[nx*ny]))

        ! Check if we are optimizing all basins
        if (tf_basins(1) .lt. 0) then 
            ! Optimizing all basins, set 
            ! basin list to reference list

            allocate(basin_list(size(basin_list_ref)))

            basin_list = basin_list_ref 

            nb = size(basin_list,1) 
        else

            nb = count(tf_basins .gt. 0)
            allocate(basin_list(nb))

            n = 0 
            do b = 1, size(tf_basins)
                if (tf_basins(b) .gt. 0) then 
                    n = n+1
                    basin_list(n) = tf_basins(b)
                end if
            end do

        end if 
         
        ! Loop over each basin
        do b = 1, nb 
            
            ! Get a mask of points of interest:
            ! 1. Points within the current basin 
            ! 2. Points with overburden thickness near flotation,
            !    with magnitude less than H_grnd_lim
            ! 3. Points with observed or modeled ice thickness
            mask =  (abs(basins-basin_list(b)) .lt. tol) .and. &
                    (H_grnd .lt. H_grnd_lim) .and. & 
                    (H_obs .gt. 0.0 .or. H_ice .gt. 0.0)

            ! How many points available 
            n = count(mask)
                
            if (n .gt. 0) then 
                ! Points are available for averaging 

                ! Calculate average observed thickness for masked region
                H_obs_bar = sum(H_obs,mask=mask)   / real(n,wp)

                ! Calculate average thickness and rate of change for masked region
                H_bar     = sum(H_ice,mask=mask)   / real(n,wp)
                dHdt_bar  = sum(dHicedt,mask=mask) / real(n,wp)

            else 
                ! Quantities are not defined

                H_obs_bar = missing_value 

                H_bar     = 0.0_wp
                dHdt_bar  = 0.0_wp 
            
            end if 

            
            if (H_obs_bar .ne. missing_value) then 
                ! Observed ice exists in basin, proceed with calculations
                                
                ! Get mean error for this basin
                H_bar_err = H_bar - H_obs_bar 

                ! Get adjustment rate given error in ice thickness  =========

                tf_corr_dot = &
                     1.0_wp/(tau_m*m_temp) *( (H_bar_err / tau_m) + f_damp*dHdt_bar )

                ! Apply correction to all points in basin =========

                where(basins .eq. basin_list(b))

                    tf_corr = tf_corr + tf_corr_dot*dt 

                end where 

            end if 

        end do 


        ! Ensure tf_corr is not below lower or upper limit 
        where (tf_corr .lt. tf_min) tf_corr = tf_min 
        where (tf_corr .gt. tf_max) tf_corr = tf_max 

        return 

    end subroutine optimize_tf_corr_basin

    subroutine optimize_cb_ref(cb_ref,H_ice,dHdt,z_bed,z_sl,ux,uy,H_obs,uxy_obs,H_grnd_obs, &
                                        cf_min,cf_max,dx,sigma_err,sigma_vel,tau_c,H0,scaleH,dt,fill_method,fill_dist, &
                                        cb_tgt)
        ! Update method following Lipscomb et al. (2021, tc)

        implicit none 

        real(wp), intent(INOUT) :: cb_ref(:,:) 
        real(wp), intent(IN)    :: H_ice(:,:) 
        real(wp), intent(IN)    :: dHdt(:,:) 
        real(wp), intent(IN)    :: z_bed(:,:) 
        real(wp), intent(IN)    :: z_sl(:,:) 
        real(wp), intent(IN)    :: ux(:,:) 
        real(wp), intent(IN)    :: uy(:,:) 
        real(wp), intent(IN)    :: H_obs(:,:) 
        real(wp), intent(IN)    :: uxy_obs(:,:) 
        real(wp), intent(IN)    :: H_grnd_obs(:,:) 
        real(wp), intent(IN)    :: cf_min(:,:) 
        real(wp), intent(IN)    :: cf_max(:,:) 
        real(wp), intent(IN)    :: dx 
        real(wp), intent(IN)    :: sigma_err 
        real(wp), intent(IN)    :: sigma_vel
        real(wp), intent(IN)    :: tau_c                    ! [yr]
        real(wp), intent(IN)    :: H0                       ! [m]
        logical,  intent(IN)    :: scaleH                   ! Scale ice thickness in opt with observation.
        real(wp), intent(IN)    :: dt 
        character(len=*), intent(IN) :: fill_method         ! How should missing values outside obs be filled?
        real(wp), intent(IN)    :: fill_dist                ! [km] Distance over which to smooth between nearest neighbor and minimum value
        real(wp), intent(IN), optional :: cb_tgt(:,:) 

        ! Local variables 
        integer  :: i, j, nx, ny, i1, j1 
        integer  :: im1, ip1, jm1, jp1  
        real(wp) :: f_damp   
        real(wp) :: ux_aa, uy_aa, uxy_aa
        real(wp) :: H_err_now, dHdt_now, f_vel   
        real(wp) :: xwt, ywt, xywt

        real(wp) :: f_tgt
        real(wp) :: cb_tgt_fac
        real(wp) :: tau_tgt 

        real(wp) :: cb_ref_dot 

        real(wp), allocatable   :: H_err_sm(:,:)
        real(wp), allocatable   :: H_err(:,:)
        real(wp), allocatable   :: uxy(:,:)
        real(wp), allocatable   :: uxy_err(:,:)
        real(wp), allocatable   :: cb_prev(:,:) 

        nx = size(cb_ref,1)
        ny = size(cb_ref,2)  
        
        allocate(H_err_sm(nx,ny))
        allocate(H_err(nx,ny))
        allocate(uxy(nx,ny))
        allocate(uxy_err(nx,ny))
        allocate(cb_prev(nx,ny))
        
        ! Internal parameters 
        f_damp = 2.0 

        f_tgt   = 0.05 * H0    ! [--] * [m] = [m]
        tau_tgt = tau_c        ! 500 [yr] 

        ! Store initial cb_ref solution 
        cb_prev = cb_ref 

        ! Calculate velocity magnitude and velocity error 
        uxy = calc_magnitude_from_staggered_ice(ux,uy,H_ice)
         
        uxy_err = MV 
        where(uxy_obs .ne. MV .and. uxy_obs .ne. 0.0) uxy_err = (uxy - uxy_obs)
        
        ! Calculate ice thickness error 
        H_err = H_ice - H_obs 

        ! Additionally, apply a Gaussian filter to H_err to ensure smooth transitions
        ! Apply a weighted average between smoothed and original H_err, where 
        ! slow regions get more smoothed, and fast regions use more local error 
        if (sigma_err .gt. 0.0) then

            H_err_sm = H_err  
            call filter_gaussian(H_err_sm,sigma_err,dx)

            do j = 1, ny 
            do i = 1, nx 
                f_vel = min( uxy(i,j)/sigma_vel, 1.0 )
                H_err(i,j) = (1.0-f_vel)*H_err_sm(i,j) + f_vel*H_err(i,j)  
            end do 
            end do  

        end if

        ! Initially set cf to missing value for now where no correction possible
        cb_ref = MV 

        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices
            im1 = max(i-1,1) 
            ip1 = min(i+1,nx) 
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny) 
            
            ux_aa = 0.5*(ux(i,j)+ux(im1,j))
            uy_aa = 0.5*(uy(i,j)+uy(i,jm1))
            
            uxy_aa = sqrt(ux_aa**2+uy_aa**2)

            if ( uxy(i,j) .ne. 0.0 .and. uxy_err(i,j) .ne. MV &
                            .and. (H_grnd_obs(i,j) .gt. 0.0) ) then 
                ! Update coefficient where velocity exists and 
                ! observations are not floating.

                ! Determine upstream node(s) 

                if (ux_aa .ge. 0.0) then 
                    i1 = im1
                else 
                    i1 = ip1 
                end if 

                if (uy_aa .ge. 0.0) then 
                    j1 = jm1
                else 
                    j1 = jp1  
                end if 
                
                ! Get weighted error  =========

                xywt  = abs(ux_aa)+abs(uy_aa)
                if (xywt .gt. 0.0) then 
                    xwt = abs(ux_aa) / xywt 
                    ywt = abs(uy_aa) / xywt 
                else 
                    ! Dummy weights
                    xwt   = 0.5 
                    ywt   = 0.5
                end if 

                ! Define error for ice thickness 
                H_err_now = xwt*H_err(i1,j) + ywt*H_err(i,j1) 
                dHdt_now  = xwt*dHdt(i1,j)  + ywt*dHdt(i,j1) 

                ! Determine scaling correction with respect to target cb_ref value
                if (present(cb_tgt) .and. (.not. scaleH)) then
                    cb_tgt_fac = log(cb_prev(i,j) / cb_tgt(i,j))
                else 
                    cb_tgt_fac = 0.0 
                end if 

                ! Get adjustment rate given error in ice thickness  =========
                ! jablasco
                if (scaleH) then
                        cb_ref_dot = -(cb_prev(i,j)/H_obs(i,j)) * &
                                ((H_err_now / tau_c) + f_damp*dHdt_now + (f_tgt/tau_tgt)*cb_tgt_fac)
                else
                        cb_ref_dot = -(cb_prev(i,j)/H0) * &
                                ((H_err_now / tau_c) + f_damp*dHdt_now + (f_tgt/tau_tgt)*cb_tgt_fac)
                end if
                ! Apply correction to current node =========

                cb_ref(i,j) = cb_prev(i,j) + cb_ref_dot*dt 

            end if 

        end do 
        end do 

        select case(trim(fill_method))

            case("analog")
                
                write(io_unit_err,*)
                write(io_unit_err,*) "optimize_cb_ref:: Error: &
                &fill_method='analog' is not working right now!"
                write(io_unit_err,*)
                stop 

                ! ajr: Need to adapt fill_cb_ref and fill_nearest for 2D fields
                ! of cf_min and cf_max, instead of just single values. 

                ! Fill in cb_ref for floating points using bed analogy method
                !call fill_cb_ref(cb_ref,H_ice,z_bed,z_sl,is_float_obs,cf_min,cf_max)

                ! Fill in remaining missing values with nearest neighbor or cf_min when none available
                !call fill_nearest(cb_ref,missing_value=MV,fill_value=cf_min,fill_dist=fill_dist,n=5,dx=dx)

            case("nearest")

                write(io_unit_err,*)
                write(io_unit_err,*) "optimize_cb_ref:: Error: &
                &fill_method='nearest' is not working right now!"
                write(io_unit_err,*)
                stop 

                ! ajr: Need to adapt fill_cb_ref and fill_nearest for 2D fields
                ! of cf_min and cf_max, instead of just single values. 
                
                ! Fill in remaining missing values with nearest neighbor or cf_min when none available
                !call fill_nearest(cb_ref,missing_value=MV,fill_value=cf_min,fill_dist=fill_dist,n=5,dx=dx)

            case("target")
                ! Fill in field with cg_tgt values 

                ! Ensure where obs are floating, set cb_ref = cb_tgt 
                where(H_grnd_obs .le. 0.0) cb_ref = cb_tgt 

                ! Also where no ice exists, set cb_ref = cb_tgt 
                where(H_obs .eq. 0.0) cb_ref = cb_tgt 

            case("cf_min")

                ! Ensure where obs are floating, set cb_ref = cf_min 
                where(H_grnd_obs .le. 0.0) cb_ref = cf_min 

                ! Also where no ice exists, set cb_ref = cf_min 
                where(H_obs .eq. 0.0) cb_ref = cf_min 

            case DEFAULT 

                write(io_unit_err,*)
                write(io_unit_err,*) "optimize_cb_ref:: Error: fill_method not recognized."
                write(io_unit_err,*) "fill_method = ", trim(fill_method)
                stop 

        end select
        
        ! Ensure cb_ref is not below lower or upper limit 
        where (cb_ref .lt. cf_min) cb_ref = cf_min 
        where (cb_ref .gt. cf_max) cb_ref = cf_max 

        return 

    end subroutine optimize_cb_ref

    subroutine optimize_cb_ref_vel(cb_ref,H_ice,dHdt,z_bed,z_sl,ux,uy,H_obs,uxy_obs,H_grnd_obs,duxydt, &
                                   cf_min,cf_max,f_pmp,dx,sigma_err,sigma_vel,tau_c,H0,dt,fill_method,fill_dist,cb_tgt)
        ! Update method following Lipscomb et al. but for ice velocity (2021, tc)

        implicit none 

        real(wp), intent(INOUT) :: cb_ref(:,:) 
        real(wp), intent(IN)    :: H_ice(:,:) 
        real(wp), intent(IN)    :: dHdt(:,:) 
        real(wp), intent(IN)    :: z_bed(:,:) 
        real(wp), intent(IN)    :: z_sl(:,:) 
        real(wp), intent(IN)    :: ux(:,:) 
        real(wp), intent(IN)    :: uy(:,:) 
        real(wp), intent(IN)    :: H_obs(:,:) 
        real(wp), intent(IN)    :: uxy_obs(:,:) 
        real(wp), intent(IN)    :: H_grnd_obs(:,:)
        real(wp), intent(IN)    :: duxydt(:,:) 
        real(wp), intent(IN)    :: cf_min(:,:) 
        real(wp), intent(IN)    :: cf_max(:,:)
        real(wp), intent(IN)    :: f_pmp(:,:) 
        real(wp), intent(IN)    :: dx 
        real(wp), intent(IN)    :: sigma_err 
        real(wp), intent(IN)    :: sigma_vel
        real(wp), intent(IN)    :: tau_c                  ! [yr]
        real(wp), intent(IN)    :: H0                     ! [m]
        real(wp), intent(IN)    :: dt 
        character(len=*), intent(IN) :: fill_method         ! How should missing values outside obs be filled?
        real(wp), intent(IN)    :: fill_dist                ! [km] Distance over which to smooth between nearest neighbor and minimum value
        real(wp), intent(IN), optional :: cb_tgt(:,:) 

        ! Local variables 
        integer  :: i, j, nx, ny, i1, j1 
        integer  :: im1, ip1, jm1, jp1  
        real(wp) :: f_damp   
        real(wp) :: ux_aa, uy_aa, uxy_aa
        real(wp) :: H_err_now, dHdt_now, f_vel  
        real(wp) :: uxy_err_now,duxydt_now 
        real(wp) :: xwt, ywt, xywt

        real(wp) :: f_tgt
        real(wp) :: cb_tgt_fac
        real(wp) :: tau_tgt 

        real(wp) :: cb_ref_dot 

        real(wp), allocatable   :: H_err_sm(:,:)
        real(wp), allocatable   :: H_err(:,:)
        real(wp), allocatable   :: uxy(:,:)
        real(wp), allocatable   :: uxy_err(:,:)
        real(wp), allocatable   :: cb_prev(:,:) 

        nx = size(cb_ref,1)
        ny = size(cb_ref,2)  

        allocate(H_err_sm(nx,ny))
        allocate(H_err(nx,ny))
        allocate(uxy(nx,ny))
        allocate(uxy_err(nx,ny))
        allocate(cb_prev(nx,ny))

        ! Internal parameters 
        f_damp = 2.0 

        f_tgt   = 0.05 * H0    ! [--] * [m] = [m]
        tau_tgt = tau_c        ! 500 [yr] 

        ! Store initial cb_ref solution 
        cb_prev = cb_ref 

        ! Calculate velocity magnitude and velocity error 
        uxy = calc_magnitude_from_staggered_ice(ux,uy,H_ice)
        uxy_err = MV 
        where(uxy_obs .ne. MV .and. uxy_obs .ne. 0.0) uxy_err = (uxy - uxy_obs)

        ! Calculate ice thickness error 
        H_err = H_ice - H_obs 

        ! Additionally, apply a Gaussian filter to H_err to ensure smooth transitions
        ! Apply a weighted average between smoothed and original H_err, where 
        ! slow regions get more smoothed, and fast regions use more local error 
        if (sigma_err .gt. 0.0) then

        H_err_sm = H_err  
        call filter_gaussian(H_err_sm,sigma_err,dx)

        do j = 1, ny 
        do i = 1, nx 
        f_vel = min( uxy(i,j)/sigma_vel, 1.0 )
        H_err(i,j) = (1.0-f_vel)*H_err_sm(i,j) + f_vel*H_err(i,j)  
        end do 
        end do  

        end if

        ! Initially set cf to missing value for now where no correction possible
        cb_ref = MV 

        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices
            im1 = max(i-1,1) 
            ip1 = min(i+1,nx) 
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny) 

            ux_aa = 0.5*(ux(i,j)+ux(im1,j))
            uy_aa = 0.5*(uy(i,j)+uy(i,jm1))

            uxy_aa = sqrt(ux_aa**2+uy_aa**2)

            if ( uxy(i,j) .ne. 0.0 .and. uxy_err(i,j) .ne. MV .and. (H_grnd_obs(i,j) .gt. 0.0) ) then 
                ! Update coefficient where velocity exists and 
                ! observations are not floating.

                ! Determine upstream node(s) 
                if (ux_aa .ge. 0.0) then 
                    i1 = im1
                else 
                    i1 = ip1 
                end if 

                if (uy_aa .ge. 0.0) then 
                    j1 = jm1
                else 
                    j1 = jp1  
                end if 

                ! Get weighted error  =========
                xywt  = abs(ux_aa)+abs(uy_aa)
                if (xywt .gt. 0.0) then 
                xwt = abs(ux_aa) / xywt 
                ywt = abs(uy_aa) / xywt 
                else 
                ! Dummy weights
                xwt   = 0.5 
                ywt   = 0.5
                end if 

                ! Define error for ice thickness 
                uxy_err_now = uxy_err(i,j)
                duxydt_now  = duxydt(i,j)

                ! Determine scaling correction with respect to target cb_ref value
                if (present(cb_tgt)) then
                   cb_tgt_fac = log(cb_prev(i,j) / cb_tgt(i,j))
                else 
                    cb_tgt_fac = 0.0 
                end if 

                ! Get adjustment rate given error in surface velocity =========
                if (.False.) then
                        cb_ref_dot = -(cb_prev(i,j)/H0) * &
                        ((uxy_err_now / tau_c) + f_damp*duxydt_now + (f_tgt/tau_tgt)*cb_tgt_fac)
                else
                        cb_ref_dot = -(cb_prev(i,j)/uxy_obs(i,j)) * &
                        ((uxy_err_now / tau_c) + f_damp*duxydt_now)
                end if

                ! Apply correction to current node =========

                cb_ref(i,j) = cb_prev(i,j) + cb_ref_dot*dt 

                if (.True.) then
                        cb_ref(i,j) = cb_ref(i,j)*f_pmp(i,j) + cf_min(i,j)*(1-f_pmp(i,j))
                end if

            end if 

        end do 
        end do 

        select case(trim(fill_method))

            case("analog")

                write(io_unit_err,*)
                write(io_unit_err,*) "optimize_cb_ref:: Error: &
                &fill_method='analog' is not working right now!"
                write(io_unit_err,*)
                stop 

                ! ajr: Need to adapt fill_cb_ref and fill_nearest for 2D fields
                ! of cf_min and cf_max, instead of just single values. 

                ! Fill in cb_ref for floating points using bed analogy method
                !call fill_cb_ref(cb_ref,H_ice,z_bed,z_sl,is_float_obs,cf_min,cf_max)

                ! Fill in remaining missing values with nearest neighbor or cf_min when none available
                !call fill_nearest(cb_ref,missing_value=MV,fill_value=cf_min,fill_dist=fill_dist,n=5,dx=dx)

            case("nearest")

                write(io_unit_err,*)
                write(io_unit_err,*) "optimize_cb_ref:: Error: &
                &fill_method='nearest' is not working right now!"
                write(io_unit_err,*)
                stop 

                ! ajr: Need to adapt fill_cb_ref and fill_nearest for 2D fields
                ! of cf_min and cf_max, instead of just single values. 

                ! Fill in remaining missing values with nearest neighbor or cf_min when none available
                !call fill_nearest(cb_ref,missing_value=MV,fill_value=cf_min,fill_dist=fill_dist,n=5,dx=dx)

            case("target")
                ! Fill in field with cg_tgt values 

                ! Ensure where obs are floating, set cb_ref = cb_tgt 
                where(H_grnd_obs .le. 0.0) cb_ref = cb_tgt 

                ! Also where no ice exists, set cb_ref = cb_tgt 
                where(H_obs .eq. 0.0) cb_ref = cb_tgt 

            case("cf_min")

                ! Ensure where obs are floating, set cb_ref = cf_min 
                where(H_grnd_obs .le. 0.0) cb_ref = cf_min 

                ! Also where no ice exists, set cb_ref = cf_min 
                where(H_obs .eq. 0.0) cb_ref = cf_min 

            case DEFAULT 

                write(io_unit_err,*)
                write(io_unit_err,*) "optimize_cb_ref:: Error: fill_method not recognized."
                write(io_unit_err,*) "fill_method = ", trim(fill_method)
                stop 

        end select

        ! Ensure cb_ref is not below lower or upper limit 
        where (cb_ref .lt. cf_min) cb_ref = cf_min 
        where (cb_ref .gt. cf_max) cb_ref = cf_max 

        return 

    end subroutine optimize_cb_ref_vel

    subroutine optimize_cb_ref_pc12(cb_ref,H_ice,H_ice_n,dHdt,z_bed,z_sl,ux,uy,H_obs,uxy_obs,H_grnd_obs, &
            cf_min,cf_max,dx,sigma_err,sigma_vel,tau_c,H0,dt,fill_method,fill_dist, &
            cb_tgt)
        ! Update method following Pollard & deConto (2012, tc)

        implicit none 

        real(wp), intent(INOUT) :: cb_ref(:,:) 
        real(wp), intent(IN)    :: H_ice(:,:) 
        real(wp), intent(IN)    :: H_ice_n(:,:) 
        real(wp), intent(IN)    :: dHdt(:,:) 
        real(wp), intent(IN)    :: z_bed(:,:) 
        real(wp), intent(IN)    :: z_sl(:,:) 
        real(wp), intent(IN)    :: ux(:,:) 
        real(wp), intent(IN)    :: uy(:,:) 
        real(wp), intent(IN)    :: H_obs(:,:) 
        real(wp), intent(IN)    :: uxy_obs(:,:) 
        real(wp), intent(IN)    :: H_grnd_obs(:,:) 
        real(wp), intent(IN)    :: cf_min(:,:) 
        real(wp), intent(IN)    :: cf_max(:,:) 
        real(wp), intent(IN)    :: dx 
        real(wp), intent(IN)    :: sigma_err 
        real(wp), intent(IN)    :: sigma_vel
        real(wp), intent(IN)    :: tau_c                  ! [yr]
        real(wp), intent(IN)    :: H0                     ! [m]
        real(wp), intent(IN)    :: dt 
        character(len=*), intent(IN) :: fill_method         ! How should missing values outside obs be filled?
        real(wp), intent(IN)    :: fill_dist                ! [km] Distance over which to smooth between nearest neighbor and minimum value
        real(wp), intent(IN), optional :: cb_tgt(:,:) 

        ! Local variables 
        integer  :: i, j, nx, ny, i1, j1 
        integer  :: im1, ip1, jm1, jp1  
        real(wp) :: ux_aa, uy_aa, uxy_aa
        real(wp) :: H_err_now, H_err_now_n, dHdt_now, f_vel   
        real(wp) :: dz_now, dz_n
        real(wp) :: xwt, ywt, xywt

        real(wp) :: f_tgt
        real(wp) :: cb_tgt_fac
        real(wp) :: tau_tgt 

        real(wp) :: cb_ref_dot 

        real(wp), allocatable   :: H_err_sm(:,:),H_err_n_sm(:,:)
        real(wp), allocatable   :: H_err(:,:),H_err_n(:,:)
        real(wp), allocatable   :: uxy(:,:)
        real(wp), allocatable   :: uxy_err(:,:)
        real(wp), allocatable   :: cb_prev(:,:) 

        nx = size(cb_ref,1)
        ny = size(cb_ref,2)  

        allocate(H_err_sm(nx,ny))
        allocate(H_err_n_sm(nx,ny))
        allocate(H_err(nx,ny))
        allocate(H_err_n(nx,ny))
        allocate(uxy(nx,ny))
        allocate(uxy_err(nx,ny))
        allocate(cb_prev(nx,ny))

        ! Internal parameters 
        f_tgt   = 0.05 * H0    ! [--] * [m] = [m]
        tau_tgt = tau_c        ! 500 [yr] 

        ! Store initial cb_ref solution 
        cb_prev = cb_ref 

        ! Calculate velocity magnitude and velocity error 
        uxy = calc_magnitude_from_staggered_ice(ux,uy,H_ice)

        uxy_err = MV 
        where(uxy_obs .ne. MV .and. uxy_obs .ne. 0.0) uxy_err = (uxy - uxy_obs)

        ! Calculate ice thickness error 
        H_err   = H_ice   - H_obs ! current timestep
        H_err_n = H_ice_n - H_obs ! previous timestep

        ! Additionally, apply a Gaussian filter to H_err to ensure smooth transitions
        ! Apply a weighted average between smoothed and original H_err, where 
        ! slow regions get more smoothed, and fast regions use more local error 
        if (sigma_err .gt. 0.0) then

            H_err_sm = H_err
            H_err_n_sm = H_err_n  
            call filter_gaussian(H_err_sm,sigma_err,dx)
            call filter_gaussian(H_err_n_sm,sigma_err,dx)

            do j = 1, ny 
            do i = 1, nx 
                f_vel = min( uxy(i,j)/sigma_vel, 1.0 )
                H_err(i,j)   = (1.0-f_vel)*H_err_sm(i,j)   + f_vel*H_err(i,j)  
                H_err_n(i,j) = (1.0-f_vel)*H_err_n_sm(i,j) + f_vel*H_err_n(i,j)
            end do 
            end do  

        end if

        ! Initially set cf to missing value for now where no correction possible
        cb_ref = MV 

        do j = 1, ny 
        do i = 1, nx 

        ! Get neighbor indices
        im1 = max(i-1,1) 
        ip1 = min(i+1,nx) 
        jm1 = max(j-1,1) 
        jp1 = min(j+1,ny) 

        ux_aa = 0.5*(ux(i,j)+ux(im1,j))
        uy_aa = 0.5*(uy(i,j)+uy(i,jm1))

        uxy_aa = sqrt(ux_aa**2+uy_aa**2)

        if ( uxy(i,j) .ne. 0.0 .and. uxy_err(i,j) .ne. MV &
            .and. (H_grnd_obs(i,j) .gt. 0.0) ) then 
            ! Update coefficient where velocity exists and observations are not floating.

            ! Determine upstream node(s) 
            if (ux_aa .ge. 0.0) then 
                i1 = im1
            else 
                i1 = ip1 
            end if 

            if (uy_aa .ge. 0.0) then 
                j1 = jm1
            else 
                j1 = jp1  
            end if 

            ! Get weighted error  =========
            xywt  = abs(ux_aa)+abs(uy_aa)
            if (xywt .gt. 0.0) then 
                xwt = abs(ux_aa) / xywt 
                ywt = abs(uy_aa) / xywt 
            else 
                ! Dummy weights
                xwt   = 0.5 
                ywt   = 0.5
            end if 

            ! Define error for ice thickness 
            H_err_now   = xwt*H_err(i1,j)   + ywt*H_err(i,j1)
            H_err_now_n = xwt*H_err_n(i1,j) + ywt*H_err_n(i,j1) 
            dHdt_now    = xwt*dHdt(i1,j)    + ywt*dHdt(i,j1) 

            ! Determine scaling correction with respect to target cb_ref value
            if (present(cb_tgt)) then
                cb_tgt_fac = log(cb_prev(i,j) / cb_tgt(i,j))
            else 
                cb_tgt_fac = 0.0 
            end if 

            ! Get adjustment rate given error in ice thickness  =========
            if (.False.) then
                ! current error
                dz_now   = MAX(-1.5,MIN(1.5,H_err_now/H0))
                ! previous error
                dz_n     = MAX(-1.5,MIN(1.5,H_err_now_n/H0))
            else
                ! current error
                dz_now   = MAX(-1.5,MIN(1.5,H_err_now/H_obs(i,j)))
                ! previous error
                dz_n     = MAX(-1.5,MIN(1.5,H_err_now_n/H_obs(i,j)))
            end if

            ! Only updates if the current error is smaller
            if (dz_now .lt. dz_n) then
                cb_ref(i,j) = cb_prev(i,j)*(10**(-dz_now)) 
            else
                cb_ref(i,j) = cb_prev(i,j)
            end if

        end if 

        end do 
        end do 

        select case(trim(fill_method))

        case("analog")

        write(io_unit_err,*)
        write(io_unit_err,*) "optimize_cb_ref:: Error: &
        &fill_method='analog' is not working right now!"
        write(io_unit_err,*)
        stop 

        ! ajr: Need to adapt fill_cb_ref and fill_nearest for 2D fields
        ! of cf_min and cf_max, instead of just single values. 

        ! Fill in cb_ref for floating points using bed analogy method
        !call fill_cb_ref(cb_ref,H_ice,z_bed,z_sl,is_float_obs,cf_min,cf_max)

        ! Fill in remaining missing values with nearest neighbor or cf_min when none available
        !call fill_nearest(cb_ref,missing_value=MV,fill_value=cf_min,fill_dist=fill_dist,n=5,dx=dx)

        case("nearest")

        write(io_unit_err,*)
        write(io_unit_err,*) "optimize_cb_ref:: Error: &
        &fill_method='nearest' is not working right now!"
        write(io_unit_err,*)
        stop 

        ! ajr: Need to adapt fill_cb_ref and fill_nearest for 2D fields
        ! of cf_min and cf_max, instead of just single values. 

        ! Fill in remaining missing values with nearest neighbor or cf_min when none available
        !call fill_nearest(cb_ref,missing_value=MV,fill_value=cf_min,fill_dist=fill_dist,n=5,dx=dx)

        case("target")
        ! Fill in field with cg_tgt values 

        ! Ensure where obs are floating, set cb_ref = cb_tgt 
        where(H_grnd_obs .le. 0.0) cb_ref = cb_tgt 

        ! Also where no ice exists, set cb_ref = cb_tgt 
        where(H_obs .eq. 0.0) cb_ref = cb_tgt 

        case("cf_min")

        ! Ensure where obs are floating, set cb_ref = cf_min 
        where(H_grnd_obs .le. 0.0) cb_ref = cf_min 

        ! Also where no ice exists, set cb_ref = cf_min 
        where(H_obs .eq. 0.0) cb_ref = cf_min 

        case DEFAULT 

        write(io_unit_err,*)
        write(io_unit_err,*) "optimize_cb_ref:: Error: fill_method not recognized."
        write(io_unit_err,*) "fill_method = ", trim(fill_method)
        stop 

        end select

        ! Ensure cb_ref is not below lower or upper limit 
        where (cb_ref .lt. cf_min) cb_ref = cf_min 
        where (cb_ref .gt. cf_max) cb_ref = cf_max 

        return 

    end subroutine optimize_cb_ref_pc12

    subroutine fill_cb_ref(cb_ref,H_ice,z_bed,z_sl,is_float_obs,cf_min,cf_max,rho_ice,rho_sw)
        ! Fill points that cannot be optimized with 
        ! analagous values from similar bed elevations 

        implicit none 

        real(wp), intent(INOUT) :: cb_ref(:,:) 
        real(wp), intent(IN)    :: H_ice(:,:) 
        real(wp), intent(IN)    :: z_bed(:,:) 
        real(wp), intent(IN)    :: z_sl(:,:) 
        logical,  intent(IN)    :: is_float_obs(:,:) 
        real(wp), intent(IN)    :: cf_min 
        real(wp), intent(IN)    :: cf_max
        real(wp), intent(IN)    :: rho_ice 
        real(wp), intent(IN)    :: rho_sw

        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: k, nbnd, nlev, n 
        logical :: is_float 
        real(wp) :: rho_sw_ice
        real(wp), allocatable :: z_bnd(:) 
        real(wp), allocatable :: z_lev(:) 
        real(wp), allocatable :: cf_lev(:) 
        logical,  allocatable :: mask(:,:) 

        logical,  allocatable :: mask_now(:,:) 
        real(wp), allocatable :: H_grnd(:,:) 
        real(wp), allocatable :: wts_now(:,:) 

        nx = size(cb_ref,1) 
        ny = size(cb_ref,2) 

        nbnd = 13
        nlev = nbnd-1
        allocate(z_bnd(nbnd))
        allocate(z_lev(nlev))
        allocate(cf_lev(nlev))

        allocate(mask(nx,ny)) 
        allocate(mask_now(nx,ny))
        allocate(H_grnd(nx,ny)) 
        allocate(wts_now(nx,ny)) 

        z_bnd = missing_value
        z_lev = missing_value

        ! Determine z_bed bin boundaries
        z_bnd = [-2500.0,-2000.0,-1000.0,-500.0,-400.0, &
                            -300.0,-200.0,-100.0,0.0,100.0,200.0,300.0,500.0]

        ! Calculate z_bed bin centers
        do k = 1, nlev 
            z_lev(k) = 0.5_wp*(z_bnd(k) + z_bnd(k+1))
        end do 

        ! Define mask 
        mask =  (H_ice .gt. 0.0_wp)  .and. &
                (.not. is_float_obs) .and. &
                (cb_ref .ne. mv)     .and. &
                (cb_ref .ne. cf_max)

        ! Calculate H_grnd (distance to flotation)
        H_grnd = H_ice - rho_sw_ice*max(z_sl-z_bed,0.0_wp)

        ! Determine mean values of cb_ref for each bin based 
        ! on values available for grounded ice 
        cf_lev(1) = cf_min 
        do k = 2, nlev 

            ! Restrict mask to points within current bin 
            mask_now = (z_bed .ge. z_bnd(k-1) .and. z_bed .lt. z_bnd(k) .and. mask)

            ! Set weighting as a function of distance to floating (ie, H_grnd)
            wts_now = H_grnd 
            where (wts_now .lt. 0.0_wp)    wts_now = 0.0_wp 
            where (wts_now .gt. 2000.0_wp) wts_now = 2000.0_wp 
            wts_now = 2000.0_wp - wts_now 

            ! Limit to masked area too 
            where (.not. mask_now) wts_now = 0.0_wp 

            n = count(wts_now .gt. 0.0_wp)

            if (n .gt. 0) then
                ! Analog points exist for this bin

                ! Normalize weights 
                wts_now = wts_now / sum(wts_now)
                
                ! Get weighted average of cb_ref for this bin
                cf_lev(k) = sum(cb_ref*wts_now)

            else
                ! No analog points available, assign the minimum value

                cf_lev(k) = cf_min

            end if 

        end do 
        
        ! Perform linear interpolation at points of interest 

        rho_sw_ice = rho_sw/rho_ice ! Ratio of density of seawater to ice [--]
        
        do j = 1, ny 
        do i = 1, nx 

            ! Determine if current point is floating 
            is_float = H_ice(i,j) - rho_sw_ice*max(z_sl(i,j)-z_bed(i,j),0.0_wp) .le. 0.0_wp 

            if (is_float .or. is_float_obs(i,j)) then 

                cb_ref(i,j) = interp_linear(z_lev,cf_lev,xout=z_bed(i,j))

            end if 

        end do 
        end do 
        
        return 

    end subroutine fill_cb_ref

    function interp_linear(x,y,xout) result(yout)
        ! Simple linear interpolation of a point

        implicit none 
 
        real(wp), dimension(:), intent(IN) :: x, y
        real(wp), intent(IN) :: xout
        real(wp) :: yout 
        integer :: i, j, n, nout 
        real(wp) :: alph

        n    = size(x) 

        if (xout .lt. x(1)) then
            yout = y(1)
        else if (xout .gt. x(n)) then
            yout = y(n)
        else
            do j = 1, n 
                if (x(j) .ge. xout) exit 
            end do

            if (j .eq. 1) then 
                yout = y(1) 
            else if (j .eq. n+1) then 
                yout = y(n)
            else 
                alph = (xout - x(j-1)) / (x(j) - x(j-1))
                yout = y(j-1) + alph*(y(j) - y(j-1))
            end if 
        end if 

        return 

    end function interp_linear
    
    subroutine update_mb_corr(mb_corr,H_ice,H_obs,tau)

        implicit none 

        real(wp), intent(OUT) :: mb_corr(:,:)     ! [m/a] Mass balance correction term 
        real(wp), intent(IN)  :: H_ice(:,:)       ! [m] Simulated ice thickness
        real(wp), intent(IN)  :: H_obs(:,:)       ! [m] Target observed ice thickness
        real(wp), intent(IN)  :: tau              ! [a] Relaxation time constant 

        mb_corr = -(H_ice - H_obs) / tau 

        return 

    end subroutine update_mb_corr

    subroutine guess_cb_ref(cb_ref,tau_d,uxy_obs,H_obs,H_grnd,u0,cf_min,cf_max,rho_ice,g)
        ! Use suggestion by Morlighem et al. (2013) to guess friction
        ! assuming tau_b ~ tau_d, and u_b = u_obs:
        !
        ! For a linear law, tau_b = beta * u_b, so 
        ! beta = tau_b / u_b = tau_d / (u_obs+ebs), ebs=0.1 to avoid divide by zero 
        ! beta = cb_ref/u0 * N_eff, so:
        ! cb_ref = (tau_d/(u_obs+ebs)) * (u0/N_eff)

        implicit none 

        real(wp), intent(OUT) :: cb_ref(:,:) 
        real(wp), intent(IN)  :: tau_d(:,:) 
        real(wp), intent(IN)  :: uxy_obs(:,:) 
        real(wp), intent(IN)  :: H_obs(:,:)
        real(wp), intent(IN)  :: H_grnd(:,:)
        real(wp), intent(IN)  :: u0 
        real(wp), intent(IN)  :: cf_min 
        real(wp), intent(IN)  :: cf_max  
        real(wp), intent(IN)  :: rho_ice 
        real(wp), intent(IN)  :: g 
        
        ! Local variables 
        real(wp), parameter :: ebs = 0.1          ! [m/yr] To avoid divide by zero 

        where (H_obs .eq. 0.0_prec .or. H_grnd .eq. 0.0_prec) 
            ! Set floating or ice-free points to minimum 
            cb_ref = cf_min 

        elsewhere 
            ! Apply equation 

            ! Linear law: 
            cb_ref = (tau_d / (uxy_obs + ebs)) * (u0 / (rho_ice*g*H_obs + 1.0_prec))

        end where 

        where (cb_ref .gt. cf_max) cb_ref = cf_max 

        return 

    end subroutine guess_cb_ref

    

    subroutine fill_nearest(var,missing_value,fill_value,fill_dist,n,dx)

        implicit none 

        real(wp), intent(INOUT) :: var(:,:)
        real(wp), intent(IN)    :: missing_value
        real(wp), intent(IN)    :: fill_value 
        real(wp), intent(IN)    :: fill_dist          ! [km]
        integer,    intent(IN)    :: n                  ! Average of n neighbors 
        real(wp), intent(IN)    :: dx                 ! [m] 

        ! Local variables 
        integer :: i, j, nx, ny, i1, j1, q, n_now, ij(2) 
        integer :: ntot 
        real(wp) :: dist_now, f_d 

        real(wp), allocatable :: var0(:,:) 
        real(wp), allocatable :: dist(:,:) 

        nx = size(var,1)
        ny = size(var,2) 

        allocate(var0(nx,ny)) 
        allocate(dist(nx,ny)) 

        ! Store initial field 
        var0 = var 

        ntot = 0 

        ! Loop over missing values, look for nearest non-missing neighbor
        do j = 1, ny 
        do i = 1, nx 

            if (var(i,j) .eq. missing_value) then 
                ! Find a neighbor value in var0 

                ! Populate distance matrix where necessary 
                dist = MV 
                do j1 = 1, ny 
                do i1 = 1, nx 
                    if (var0(i1,j1) .ne. MV) then 
                        dist(i1,j1) = sqrt( real( (i1-i)**2 + (j1-j)**2 ) ) * dx
                    end if 
                end do 
                end do 

                n_now    = 0 
                var(i,j) = 0.0
                dist_now = 0.0 

                do q = 1, n 
                    ! Loop over nearest neighbors to get average 

                    ! Find minimum populated neighbor 
                    ij = minloc(dist,mask=dist.ne.MV .and. var0.ne.MV)

                    ! Check if no neighbors found 
                    if (ij(1) .eq. 0) exit 

                    ! Populate with neighbor value 
                    var(i,j) = var(i,j) + var0(ij(1),ij(2))
                    dist_now = dist_now + dist(ij(1),ij(2))
                    n_now = n_now + 1 

                    ! Reset distance of neighbor to zero so it cannot be used again
                    dist(ij(1),ij(2)) = MV 
                end do 

                ! If no neighbors found, use fill value 
                if (n_now .eq. 0) var(i,j) = fill_value 

                ! Take average if multiple points used 
                if (n_now .gt. 1) then 
                    
                    ! Determine mean distance to neighbors and weighting function versus distance
                    dist_now = dist_now / real(n_now,prec) 
                    f_d      = 1.0 - min( dist_now/fill_dist, 1.0 )

                    ! Apply weighted average of mean neighbor value and fill value 
                    var(i,j) = f_d * (var(i,j) / real(n_now,prec)) + (1.0-f_d)*fill_value
                    
                end if 

                ! Add this missing point to total for diagnostics 
                ntot = ntot + 1 
            end if 

        end do
        end do 

        return 

    end subroutine fill_nearest

    
    subroutine wtd_mean(var_ave,var,wts)
        ! wts == gauss_values(dx,dy,sigma,n)

        implicit none

        real(wp), intent(OUT) :: var_ave 
        real(wp), intent(IN)  :: var(:,:) 
        real(wp), intent(IN)  :: wts(:,:) 

        ! Local variables 
        real(wp) :: wts_tot 
        real(wp) :: wts_norm(size(wts,1),size(wts,2))

        wts_tot = sum(wts) 
        if (wts_tot .gt. 0.0) then 
            wts_norm = wts / wts_tot 
        else 
            wts_norm = 0.0 
        end if 

        var_ave = sum(var*wts_norm) 

        return 

    end subroutine wtd_mean

    ! Yelmo functions duplicated here to avoid dependency

    ! From yelmo_tools.f90:
    function calc_magnitude_from_staggered_ice(u,v,H,boundaries) result(umag)
        ! Calculate the centered (aa-nodes) magnitude of a vector 
        ! from the staggered (ac-nodes) components

        implicit none 
        
        real(wp), intent(IN)  :: u(:,:), v(:,:), H(:,:) 
        real(wp) :: umag(size(u,1),size(u,2)) 
        character(len=*), intent(IN), optional :: boundaries 

        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: ip1, jp1, im1, jm1
        real(wp) :: unow, vnow 
        real(wp) :: f1, f2, H1, H2 
        
        nx = size(u,1)
        ny = size(u,2) 

        umag = 0.0_prec 

        do j = 1, ny 
        do i = 1, nx 

            im1 = max(i-1,1)
            jm1 = max(j-1,1)
            ip1 = min(i+1,nx)
            jp1 = min(j+1,ny)

            ! x-direction =====

            H1 = 0.5_prec*(H(im1,j)+H(i,j))
            H2 = 0.5_prec*(H(i,j)+H(ip1,j))

            f1 = 0.5_prec 
            f2 = 0.5_prec 
            if (H1 .eq. 0.0) f1 = 0.0_prec  
            if (H2 .eq. 0.0) f2 = 0.0_prec   

            if (f1+f2 .gt. 0.0) then 
                unow = (f1*u(im1,j) + f2*u(i,j)) / (f1+f2)
                if (abs(unow) .lt. tol_underflow) unow = 0.0_prec 
            else 
                unow = 0.0 
            end if 

            ! y-direction =====

            H1 = 0.5_prec*(H(i,jm1)+H(i,j))
            H2 = 0.5_prec*(H(i,j)+H(i,jp1))

            f1 = 0.5_prec 
            f2 = 0.5_prec 
            if (H1 .eq. 0.0) f1 = 0.0_prec  
            if (H2 .eq. 0.0) f2 = 0.0_prec   

            if (f1+f2 .gt. 0.0) then 
                vnow = (f1*v(i,jm1) + f2*v(i,j)) / (f1+f2)
                if (abs(vnow) .lt. tol_underflow) vnow = 0.0_prec 
            else 
                vnow = 0.0 
            end if 

            umag(i,j) = sqrt(unow*unow+vnow*vnow)
        end do 
        end do 

        if (present(boundaries)) then 
            ! Apply conditions at boundaries of domain 

            if (trim(boundaries) .eq. "periodic") then 

                umag(1,:)  = umag(nx-1,:) 
                umag(nx,:) = umag(2,:) 
                 
                umag(:,1)  = umag(:,ny-1)
                umag(:,ny) = umag(:,2) 
                
            end if 

        end if 

        return

    end function calc_magnitude_from_staggered_ice
    
    ! From yelmo_tools.f90:
    function gauss_values(dx,dy,sigma,n) result(filt)
        ! Calculate 2D Gaussian smoothing kernel
        ! https://en.wikipedia.org/wiki/Gaussian_blur

        implicit none 

        real(wp), intent(IN) :: dx 
        real(wp), intent(IN) :: dy 
        real(wp), intent(IN) :: sigma 
        integer,    intent(IN) :: n 
        real(wp) :: filt(n,n) 

        ! Local variables 
        real(wp) :: x, y  
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

    subroutine unique(xu,x)
        ! Return only the unique values of a vector
        ! http://rosettacode.org/wiki/Remove_duplicate_elements#Fortran

        implicit none 

        real(wp), allocatable :: xu(:)      ! The output 
        real(wp) :: x(:)                    ! The input
        
        ! Local variables 
        integer :: i, j, n
        real(wp) :: res(size(x))            ! The unique values
        logical :: found 

        real(wp), parameter :: tol = 1e-5_wp
        
        n = 1
        res(1) = x(1)
        do i=2,size(x)
            found = .FALSE.
            do j=1,n
                if (abs(res(j)-x(i)) .le. tol) then 
                   ! Found a match so start looking again
                   found = .TRUE. 
                   cycle 
                end if
            end do
            ! No match found so add it to the output
            if (.not. found) then 
                n = n + 1
                res(n) = x(i)
            end if 
        end do

        ! Store output in properly sized output vector
        if(allocated(xu)) deallocate(xu)
        allocate(xu(n))
        xu = res(1:n)

        return 

    end subroutine unique

! ============
!
! OBSOLETE OPTIMIZATION ROUTINES
!
! ============

    subroutine update_cb_ref_errscaling(cb_ref,H_ice,z_bed,ux,uy,H_obs,uxy_obs,is_float_obs, &
                                        dx,cf_min,cf_max,sigma_err,sigma_vel,err_scale,fill_dist,optvar)

        implicit none 

        real(wp), intent(INOUT) :: cb_ref(:,:) 
        real(wp), intent(IN)    :: H_ice(:,:) 
        real(wp), intent(IN)    :: z_bed(:,:) 
        real(wp), intent(IN)    :: ux(:,:) 
        real(wp), intent(IN)    :: uy(:,:) 
        real(wp), intent(IN)    :: H_obs(:,:) 
        real(wp), intent(IN)    :: uxy_obs(:,:) 
        logical,  intent(IN)    :: is_float_obs(:,:) 
        real(wp), intent(IN)    :: dx 
        real(wp), intent(IN)    :: cf_min 
        real(wp), intent(IN)    :: cf_max
        real(wp), intent(IN)    :: sigma_err 
        real(wp), intent(IN)    :: sigma_vel
        real(wp), intent(IN)    :: err_scale              ! [m] or [m/a]
        real(wp), intent(IN)    :: fill_dist              ! [km] Distance over which to smooth between nearest neighbor and minimum value
        character(len=*), intent(IN) :: optvar

        ! Local variables 
        integer  :: i, j, nx, ny, i1, j1  
        real(wp) :: dx_km, f_err, f_err_lim, f_scale, f_scale_vel   
        real(wp) :: ux_aa, uy_aa, uxy_aa
        real(wp) :: err_now, err_now_vel, err_scale_vel, f_err_vel, f_vel   

        real(wp) :: xwt, ywt, xywt   

        real(wp), allocatable   :: H_err_sm(:,:)
        real(wp), allocatable   :: H_err(:,:)
        real(wp), allocatable   :: uxy(:,:)
        real(wp), allocatable   :: uxy_err(:,:)
        real(wp), allocatable   :: cf_prev(:,:) 

        nx = size(cb_ref,1)
        ny = size(cb_ref,2) 

        dx_km = dx*1e-3  
        
        allocate(H_err_sm(nx,ny))
        allocate(H_err(nx,ny))
        allocate(uxy(nx,ny))
        allocate(uxy_err(nx,ny))
        allocate(cf_prev(nx,ny))

        ! Optimization parameters 
        f_err_lim = 1.5              ! [--] 

        ! Store initial cb_ref solution 
        cf_prev = cb_ref 

        ! Calculate ice thickness error 
        H_err = H_ice - H_obs 

        ! Calculate velocity magnitude and velocity error 
        uxy = calc_magnitude_from_staggered_ice(ux,uy,H_ice)
         
        uxy_err = MV 
        where(uxy_obs .ne. MV .and. uxy_obs .ne. 0.0) uxy_err = (uxy - uxy_obs)

        ! Additionally, apply a Gaussian filter to H_err to ensure smooth transitions
        ! Apply a weighted average between smoothed and original H_err, where 
        ! slow regions get more smoothed, and fast regions use more local error 
        if (sigma_err .gt. 0.0) then
            H_err_sm = H_err  
            call filter_gaussian(var=H_err_sm,sigma=dx_km*sigma_err,dx=dx_km)

            do j = 1, ny 
            do i = 1, nx 
                f_vel = min( uxy(i,j)/sigma_vel, 1.0 )
                H_err(i,j) = (1.0-f_vel)*H_err_sm(i,j) + f_vel*H_err(i,j)  
            end do 
            end do  

        end if 

        ! Initially set cf to missing value for now where no correction possible
        cb_ref = MV 

        do j = 3, ny-2 
        do i = 3, nx-2 

            ux_aa = 0.5*(ux(i,j)+ux(i+1,j))
            uy_aa = 0.5*(uy(i,j)+uy(i,j+1))
            
            uxy_aa = sqrt(ux_aa**2+uy_aa**2)

            if ( uxy(i,j) .ne. 0.0 .and. uxy_err(i,j) .ne. MV ) then 
                ! Update coefficient where velocity exists

                ! Determine upstream node(s) 

                if (ux_aa .ge. 0.0) then 
                    i1 = i-1 
                else 
                    i1 = i+1 
                end if 

                if (uy_aa .ge. 0.0) then 
                    j1 = j-1
                else 
                    j1 = j+1  
                end if 
                
                ! Get weighted error  =========

                xwt   = 0.5 
                ywt   = 0.5
                xywt  = abs(ux(i1,j))+abs(uy(i,j1))

                if (xywt .gt. 0.0) then 
                    xwt = abs(ux(i1,j)) / xywt 
                    ywt = abs(uy(i,j1)) / xywt 
                end if 

                ! Define error for ice thickness 
                err_now = xwt*H_err(i1,j) + ywt*H_err(i,j1) 
                
                ! Define error for surface velocity 
                err_now_vel = xwt*uxy_err(i1,j) + ywt*uxy_err(i,j1)
                err_now_vel = -err_now_vel  ! Make negative to invert relationship (higher vel, higher friction)
                
                if (trim(optvar) .eq. "vel") then
                    ! If using the velocity method, then set error to velocity error 

                    err_now = err_now_vel 
                    
                end if 

                ! Get adjustment rate given error in ice thickness  =========

                f_err   = err_now / err_scale
                f_err   = max(f_err,-f_err_lim)
                f_err   = min(f_err, f_err_lim)
                f_scale = 10.0**(-f_err) 

if (.FALSE.) then 
    ! Add contribution from velocity (experimental - not working)

                err_scale_vel = 200.0 
                f_vel         = 0.20     ! 20% velocity contribution 

                f_err_vel   = err_now_vel / err_scale_vel
                f_err_vel   = max(f_err_vel,-f_err_lim)
                f_err_vel   = min(f_err_vel, f_err_lim)
                f_scale_vel = 10.0**(-f_err_vel) 

                f_scale = (1.0-f_vel)*f_scale + f_vel*f_scale_vel 
end if 

                ! Apply correction to current node =========

                cb_ref(i,j) = f_scale * cf_prev(i,j) 

            end if 

        end do 
        end do 

        ! Fill in missing values with nearest neighbor or cf_min when none available
        call fill_nearest(cb_ref,missing_value=MV,fill_value=cf_min,fill_dist=fill_dist,n=5,dx=dx)

        ! Ensure cb_ref is not below lower or upper limit 
        where (cb_ref .lt. cf_min) cb_ref = cf_min 
        where (cb_ref .gt. cf_max) cb_ref = cf_max 

        ! Additionally, apply a Gaussian filter to cb_ref to ensure smooth transitions
        !call filter_gaussian(var=cb_ref,sigma=dx_km*0.2,dx=dx_km)     !,mask=err_z_srf .ne. 0.0)
        
        ! Ensure where obs are floating, set cb_ref = cf_min 
        !where(is_float_obs) cb_ref = cf_min 

        ! Also where no ice exists, set cb_ref = cf_min 
        !where(H_obs .eq. 0.0) cb_ref = cf_min 

        return 

    end subroutine update_cb_ref_errscaling

    subroutine update_cb_ref_thickness_ratio(cb_ref,H_ice,z_bed,ux,uy,uxy_i,uxy_b,H_obs,dx,cf_min,cf_max)

        implicit none 

        real(wp), intent(INOUT) :: cb_ref(:,:) 
        real(wp), intent(IN)    :: H_ice(:,:) 
        real(wp), intent(IN)    :: z_bed(:,:) 
        real(wp), intent(IN)    :: ux(:,:)        ! Depth-averaged velocity (ux_bar)
        real(wp), intent(IN)    :: uy(:,:)        ! Depth-averaged velocity (uy_bar)
        real(wp), intent(IN)    :: uxy_i(:,:)     ! Internal shear velocity magnitude 
        real(wp), intent(IN)    :: uxy_b(:,:)     ! Basal sliding velocity magnitude 
        real(wp), intent(IN)    :: H_obs(:,:) 
        real(wp), intent(IN)    :: dx 
        real(wp), intent(IN)    :: cf_min 
        real(wp), intent(IN)    :: cf_max

        ! Local variables 
        integer :: i, j, nx, ny, i1, j1, n 
        real(wp) :: f_err, f_vel, f_corr, dx_km 
        real(wp) :: ux_aa, uy_aa 
        real(wp) :: H_ice_now, H_obs_now 

        real(wp), allocatable   :: cf_prev(:,:) 
        real(wp) :: wts0(5,5), wts(5,5) 

        real(wp),parameter :: exp1 = 2.0

        nx = size(cb_ref,1)
        ny = size(cb_ref,2) 

        dx_km = dx*1e-3  
        
        allocate(cf_prev(nx,ny))

        ! Get Gaussian weights 
        wts0 = gauss_values(dx_km,dx_km,sigma=dx_km*1.5,n=5)

!         do i = 1, 5 
!         write(*,*) wts0(i,:) 
!         end do 
!         stop 

        ! Store initial cb_ref solution 
        cf_prev = cb_ref 

        do j = 3, ny-2 
        do i = 3, nx-2 

            if ( abs(H_ice(i,j) - H_obs(i,j)) .ne. 0.0) then 
                ! Update where thickness error exists

                ! Determine downstream point to apply changes

                ux_aa = 0.5*(ux(i,j)+ux(i+1,j))
                uy_aa = 0.5*(uy(i,j)+uy(i,j+1))
                
                if ( abs(ux_aa) .gt. abs(uy_aa) ) then 
                    ! Downstream in x-direction 
                    j1 = j 
                    if (ux_aa .lt. 0.0) then 
                        i1 = i-1 
                    else
                        i1 = i+1
                    end if 

                else 
                    ! Downstream in y-direction 
                    i1 = i 
                    if (uy_aa .lt. 0.0) then 
                        j1 = j-1
                    else
                        j1 = j+1
                    end if 

                end if 

                ! Calculate thickness error ratio 
!                 f_err = H_ice(i,j) / max(H_obs(i,j),1e-1)
                
                wts = wts0 
                !where( H_ice(i-2:i+2,j-2:j+2) .eq. 0.0) wts = 0.0 
                call wtd_mean(H_ice_now,H_ice(i-2:i+2,j-2:j+2),wts) 

                wts = wts0 
                !where( H_obs(i-2:i+2,j-2:j+2) .eq. 0.0) wts = 0.0 
                call wtd_mean(H_obs_now,H_obs(i-2:i+2,j-2:j+2),wts) 
                
!                 n = count(H_ice(i-1:i+1,j-1:j+1).gt.0.0)
!                 if (n .gt. 0) then
!                     H_ice_now = sum(H_ice(i-1:i+1,j-1:j+1),mask=H_ice(i-1:i+1,j-1:j+1).gt.0.0) / real(n,prec)
!                 else 
!                     H_ice_now = 0.0 
!                 end if 

!                 n = count(H_obs(i-1:i+1,j-1:j+1).gt.0.0)
!                 if (n .gt. 0) then
!                     H_obs_now = sum(H_obs(i-1:i+1,j-1:j+1),mask=H_obs(i-1:i+1,j-1:j+1).gt.0.0) / real(n,prec)
!                 else 
!                     H_obs_now = 0.0 
!                 end if 
                
                f_err = ( H_ice_now / max(H_obs_now,1e-1) )
                
                ! Calculate ratio of deformational velocity to sliding velocity
                f_vel = uxy_i(i,j) / max(uxy_b(i,j),1e-1) 

                ! Calculate correction factor (beta_old / beta_new)
                f_corr = ( max( f_err + f_vel*(f_err-1.0_prec), 1e-1) )**exp1

                ! Apply correction to update cb_ref
                cb_ref(i1,j1) = cf_prev(i1,j1) * f_corr**(-1.0)

            end if 

        end do 
        end do 

        ! Ensure cb_ref is not below lower or upper limit 
        where (cb_ref .lt. cf_min) cb_ref = cf_min 
        where (cb_ref .gt. cf_max) cb_ref = cf_max 

        ! Additionally, apply a Gaussian filter to cb_ref to ensure smooth transitions
        call filter_gaussian(var=cb_ref,sigma=dx_km*0.25,dx=dx_km)     !,mask=err_z_srf .ne. 0.0)
        
        ! Also where no ice exists, set cb_ref = cf_min 
        where(H_ice .eq. 0.0) cb_ref = cf_min 

        return 

    end subroutine update_cb_ref_thickness_ratio

end module ice_optimization

