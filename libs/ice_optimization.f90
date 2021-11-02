module ice_optimization

    use yelmo_defs, only : sp, dp, wp, prec, pi, missing_value, mv, tol_underflow, rho_ice, rho_sw, g 

    use gaussian_filter 
    
    implicit none 

    private 

    public :: update_tf_corr_l21
    public :: update_cf_ref_errscaling_l21
    public :: update_cf_ref_errscaling
    public :: update_cf_ref_thickness_ratio
    public :: update_mb_corr
    public :: guess_cf_ref
    public :: fill_nearest
    public :: fill_cf_ref
    public :: wtd_mean
    public :: get_opt_param

contains 
    
    subroutine update_tf_corr_l21(tf_corr,H_ice,H_grnd,dHicedt,H_obs,basins,H_grnd_lim, &
                                    tau_m,m_temp,tf_min,tf_max,dt)

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
        real(wp), intent(IN)    :: dt 

        ! Local variables
        integer :: i, j, nx, ny, n  
        integer :: b, nb  
        real(wp) :: f_damp 
        real(wp) :: H_obs_now 
        real(wp) :: H_now 
        real(wp) :: H_err_now 
        real(wp) :: dHdt_now
        real(wp) :: tf_corr_dot
        real(wp), allocatable :: basin_list(:) 
        logical,  allocatable :: mask(:,:) 

        ! Internal parameters 
        f_damp = 2.0 

        nx = size(tf_corr,1)
        ny = size(tf_corr,2) 

        allocate(mask(nx,ny)) 

        ! Determine unique basin numbers 
        call unique(basin_list,reshape(basins,[nx*ny]))
        nb = size(basin_list,1) 

        ! Loop over each basin
        do b = 1, nb 
            
            ! Get a mask of points of interest:
            ! 1. Points within the current basin 
            ! 2. Points with overburden thickness near flotation
            mask = basins .eq. basin_list(b) .and. &
                    abs(H_grnd) .lt. H_grnd_lim 

            ! Calculate averages

            n = count(mask .and. H_obs .gt. 0.0)
            if (n .gt. 0) then 
                H_obs_now = sum(H_obs,mask=mask) / real(n,wp)
            else 
                H_obs_now = missing_value 
            end if 

            n = count(mask)
            if (n .gt. 0) then 
                H_now    = sum(H_ice,mask=mask) / real(n,wp)
                dHdt_now = sum(dHicedt,mask=mask) / real(n,wp)
            else 
                H_now    = 0.0_wp
                dHdt_now = 0.0_wp 
            end if 

            if (H_obs_now .ne. missing_value) then 
                ! Observed ice exists in basin, proceed with calculations
                                
                ! Get mean error for this basin
                H_err_now = H_now - H_obs_now 

                ! Get adjustment rate given error in ice thickness  =========

                tf_corr_dot = -1.0_wp/(tau_m*m_temp) *( (H_err_now / tau_m) + f_damp*dHdt_now )

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

    end subroutine update_tf_corr_l21

    subroutine update_cf_ref_errscaling_l21(cf_ref,H_ice,dHdt,z_bed,z_sl,ux,uy,H_obs,uxy_obs,is_float_obs, &
                                        dx,cf_min,cf_max,sigma_err,sigma_vel,tau_c,H0,fill_dist,dt)
        ! Update method following Lipscomb et al. (2021, tc)

        implicit none 

        real(wp), intent(INOUT) :: cf_ref(:,:) 
        real(wp), intent(IN)    :: H_ice(:,:) 
        real(wp), intent(IN)    :: dHdt(:,:) 
        real(wp), intent(IN)    :: z_bed(:,:) 
        real(wp), intent(IN)    :: z_sl(:,:) 
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
        real(wp), intent(IN)    :: tau_c                  ! [yr]
        real(wp), intent(IN)    :: H0                     ! [m]
        real(wp), intent(IN)    :: fill_dist              ! [km] Distance over which to smooth between nearest neighbor and minimum value
        real(wp), intent(IN)    :: dt 

        ! Local variables 
        integer  :: i, j, nx, ny, i1, j1 
        integer  :: im1, ip1, jm1, jp1  
        real(wp) :: dx_km, f_damp   
        real(wp) :: ux_aa, uy_aa, uxy_aa
        real(wp) :: H_err_now, dHdt_now, f_vel   
        real(wp) :: xwt, ywt, xywt   
        real(wp) :: cf_val 

        real(wp), allocatable   :: H_err_sm(:,:)
        real(wp), allocatable   :: H_err(:,:)
        real(wp), allocatable   :: uxy(:,:)
        real(wp), allocatable   :: uxy_err(:,:)
        real(wp), allocatable   :: cf_prev(:,:) 
        real(wp), allocatable   :: cf_ref_dot(:,:)

        nx = size(cf_ref,1)
        ny = size(cf_ref,2) 

        dx_km = dx*1e-3  
        
        allocate(H_err_sm(nx,ny))
        allocate(H_err(nx,ny))
        allocate(uxy(nx,ny))
        allocate(uxy_err(nx,ny))
        allocate(cf_prev(nx,ny))
        allocate(cf_ref_dot(nx,ny)) 
        
        ! Internal parameters 
        f_damp = 2.0 

        ! Store initial cf_ref solution 
        cf_prev = cf_ref 

        ! Calculate ice thickness error 
        H_err = H_ice - H_obs 

        ! Calculate velocity magnitude and velocity error 
        uxy = calc_magnitude_from_staggered_ice(ux,uy,H_ice)
         
        uxy_err = MV 
        where(uxy_obs .ne. MV .and. uxy_obs .ne. 0.0) uxy_err = (uxy - uxy_obs)

if (.TRUE.) then 
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
end if 

        ! Initially set cf to missing value for now where no correction possible
        cf_ref = MV 

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
                            .and. (.not. is_float_obs(i,j)) ) then 
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

                ! Get adjustment rate given error in ice thickness  =========

                cf_ref_dot(i,j) = -(cf_prev(i,j)/H0)*((H_err_now / tau_c) + f_damp*dHdt_now)

                ! Apply correction to current node =========

                cf_ref(i,j) = cf_prev(i,j) + cf_ref_dot(i,j)*dt 

            end if 

        end do 
        end do 

        ! Fill in cf_ref for floating points using bed analogy method
        call fill_cf_ref(cf_ref,H_ice,z_bed,z_sl,is_float_obs,cf_min)

        ! Fill in remaining missing values with nearest neighbor or cf_min when none available
        call fill_nearest(cf_ref,missing_value=MV,fill_value=cf_min,fill_dist=fill_dist,n=5,dx=dx)

        ! Ensure cf_ref is not below lower or upper limit 
        where (cf_ref .lt. cf_min) cf_ref = cf_min 
        where (cf_ref .gt. cf_max) cf_ref = cf_max 

        ! Additionally, apply a Gaussian filter to cf_ref to ensure smooth transitions
        !call filter_gaussian(var=cf_ref,sigma=dx_km*0.2,dx=dx_km)     !,mask=err_z_srf .ne. 0.0)
        
        ! Ensure where obs are floating, set cf_ref = cf_min 
        where(is_float_obs) cf_ref = cf_min 

        ! Also where no ice exists, set cf_ref = cf_min 
        where(H_obs .eq. 0.0) cf_ref = cf_min 

        return 

    end subroutine update_cf_ref_errscaling_l21

    subroutine fill_cf_ref(cf_ref,H_ice,z_bed,z_sl,is_float_obs,cf_min)
        ! Fill points that cannot be optimized with 
        ! analagous values from similar bed elevations 

        implicit none 

        real(wp), intent(INOUT) :: cf_ref(:,:) 
        real(wp), intent(IN)    :: H_ice(:,:) 
        real(wp), intent(IN)    :: z_bed(:,:) 
        real(wp), intent(IN)    :: z_sl(:,:) 
        logical,  intent(IN)    :: is_float_obs(:,:) 
        real(wp), intent(IN)    :: cf_min 

        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: k, nlev, n 
        logical :: is_float 
        real(wp) :: rho_sw_ice
        real(wp), allocatable :: z_bnd(:) 
        real(wp), allocatable :: z_lev(:) 
        real(wp), allocatable :: cf_lev(:) 
        logical,  allocatable :: mask(:,:) 

        nx = size(cf_ref,1) 
        ny = size(cf_ref,2) 

        nlev = 10
        allocate(z_bnd(nlev+1))
        allocate(z_lev(nlev))
        allocate(cf_lev(nlev))

        allocate(mask(nx,ny)) 

        ! Determine z_bed bin boundaries and bin centers
        z_bnd = [-2000.0,-1000.0,-500.0,-400.0,-300.0,-200.0,-100.0,0.0,100.0,200.0]

        do k = 1, nlev 
            z_lev(k) = 0.5_wp*(z_bnd(k) + z_bnd(k+1))
        end do 

        ! Define mask 
        mask = (H_ice .gt. 0.0_wp) .and. (.not. is_float_obs) .and. (cf_ref .ne. mv)

        ! Determine mean values of cf_ref for each bin based 
        ! on values available for grounded ice 
        cf_lev(1) = cf_min 
        do k = 2, nlev 

            n = count(z_bed .ge. z_bnd(k-1) .and. z_bed .lt. z_bnd(k) .and. mask)

            if (n .gt. 0) then 
                cf_lev(k) = sum(cf_ref, &
                            mask=z_bed .ge. z_bnd(k-1) .and. z_bed .lt. z_bnd(k) .and. mask) &
                                    / real(n,wp)
            else 
                cf_lev(k) = cf_min 
            end if 

        end do 

        ! Perform linear interpolation at points of interest 

        rho_sw_ice = rho_sw/rho_ice ! Ratio of density of seawater to ice [--]
        
        do j = 1, ny 
        do i = 1, nx 

            ! Determine if current point is floating 
            is_float = H_ice(i,j) - rho_sw_ice*max(z_sl(i,j)-z_bed(i,j),0.0_wp) .le. 0.0_wp 

            if (is_float) then 

                cf_ref(i,j) = interp_linear(z_lev,cf_lev,xout=z_bed(i,j))

            end if 

        end do 
        end do 
        
        return 

    end subroutine fill_cf_ref

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
    
    subroutine update_cf_ref_errscaling(cf_ref,H_ice,z_bed,ux,uy,H_obs,uxy_obs,is_float_obs, &
                                        dx,cf_min,cf_max,sigma_err,sigma_vel,err_scale,fill_dist,optvar)

        implicit none 

        real(wp), intent(INOUT) :: cf_ref(:,:) 
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

        nx = size(cf_ref,1)
        ny = size(cf_ref,2) 

        dx_km = dx*1e-3  
        
        allocate(H_err_sm(nx,ny))
        allocate(H_err(nx,ny))
        allocate(uxy(nx,ny))
        allocate(uxy_err(nx,ny))
        allocate(cf_prev(nx,ny))

        ! Optimization parameters 
        f_err_lim = 1.5              ! [--] 

        ! Store initial cf_ref solution 
        cf_prev = cf_ref 

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
        cf_ref = MV 

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

                cf_ref(i,j) = f_scale * cf_prev(i,j) 

            end if 

        end do 
        end do 

        ! Fill in missing values with nearest neighbor or cf_min when none available
        call fill_nearest(cf_ref,missing_value=MV,fill_value=cf_min,fill_dist=fill_dist,n=5,dx=dx)

        ! Ensure cf_ref is not below lower or upper limit 
        where (cf_ref .lt. cf_min) cf_ref = cf_min 
        where (cf_ref .gt. cf_max) cf_ref = cf_max 

        ! Additionally, apply a Gaussian filter to cf_ref to ensure smooth transitions
        !call filter_gaussian(var=cf_ref,sigma=dx_km*0.2,dx=dx_km)     !,mask=err_z_srf .ne. 0.0)
        
        ! Ensure where obs are floating, set cf_ref = cf_min 
        !where(is_float_obs) cf_ref = cf_min 

        ! Also where no ice exists, set cf_ref = cf_min 
        !where(H_obs .eq. 0.0) cf_ref = cf_min 

        return 

    end subroutine update_cf_ref_errscaling

    subroutine update_cf_ref_thickness_ratio(cf_ref,H_ice,z_bed,ux,uy,uxy_i,uxy_b,H_obs,dx,cf_min,cf_max)

        implicit none 

        real(wp), intent(INOUT) :: cf_ref(:,:) 
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

        nx = size(cf_ref,1)
        ny = size(cf_ref,2) 

        dx_km = dx*1e-3  
        
        allocate(cf_prev(nx,ny))

        ! Get Gaussian weights 
        wts0 = gauss_values(dx_km,dx_km,sigma=dx_km*1.5,n=5)

!         do i = 1, 5 
!         write(*,*) wts0(i,:) 
!         end do 
!         stop 

        ! Store initial cf_ref solution 
        cf_prev = cf_ref 

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

                ! Apply correction to update cf_ref
                cf_ref(i1,j1) = cf_prev(i1,j1) * f_corr**(-1.0)

            end if 

        end do 
        end do 

        ! Ensure cf_ref is not below lower or upper limit 
        where (cf_ref .lt. cf_min) cf_ref = cf_min 
        where (cf_ref .gt. cf_max) cf_ref = cf_max 

        ! Additionally, apply a Gaussian filter to cf_ref to ensure smooth transitions
        call filter_gaussian(var=cf_ref,sigma=dx_km*0.25,dx=dx_km)     !,mask=err_z_srf .ne. 0.0)
        
        ! Also where no ice exists, set cf_ref = cf_min 
        where(H_ice .eq. 0.0) cf_ref = cf_min 

        return 

    end subroutine update_cf_ref_thickness_ratio

    subroutine update_mb_corr(mb_corr,H_ice,H_obs,tau)

        implicit none 

        real(wp), intent(OUT) :: mb_corr(:,:)     ! [m/a] Mass balance correction term 
        real(wp), intent(IN)  :: H_ice(:,:)       ! [m] Simulated ice thickness
        real(wp), intent(IN)  :: H_obs(:,:)       ! [m] Target observed ice thickness
        real(wp), intent(IN)  :: tau              ! [a] Relaxation time constant 

        mb_corr = -(H_ice - H_obs) / tau 

        return 

    end subroutine update_mb_corr

    subroutine guess_cf_ref(cf_ref,tau_d,uxy_obs,H_obs,H_grnd,u0,cf_min,cf_max)
        ! Use suggestion by Morlighem et al. (2013) to guess friction
        ! assuming tau_b ~ tau_d, and u_b = u_obs:
        !
        ! For a linear law, tau_b = beta * u_b, so 
        ! beta = tau_b / u_b = tau_d / (u_obs+ebs), ebs=0.1 to avoid divide by zero 
        ! beta = cf_ref/u0 * N_eff, so:
        ! cf_ref = (tau_d/(u_obs+ebs)) * (u0/N_eff)

        implicit none 

        real(wp), intent(OUT) :: cf_ref(:,:) 
        real(wp), intent(IN)  :: tau_d(:,:) 
        real(wp), intent(IN)  :: uxy_obs(:,:) 
        real(wp), intent(IN)  :: H_obs(:,:)
        real(wp), intent(IN)  :: H_grnd(:,:)
        real(wp), intent(IN)  :: u0 
        real(wp), intent(IN)  :: cf_min 
        real(wp), intent(IN)  :: cf_max  

        ! Local variables 
        real(wp), parameter :: ebs = 0.1          ! [m/yr] To avoid divide by zero 

        where (H_obs .eq. 0.0_prec .or. H_grnd .eq. 0.0_prec) 
            ! Set floating or ice-free points to minimum 
            cf_ref = cf_min 

        elsewhere 
            ! Apply equation 

            ! Linear law: 
            cf_ref = (tau_d / (uxy_obs + ebs)) * (u0 / (rho_ice*g*H_obs + 1.0_prec))

        end where 

        where (cf_ref .gt. cf_max) cf_ref = cf_max 

        return 

    end subroutine guess_cf_ref

    

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
        real(wp) :: dx_km 
        real(wp) :: dist_now, f_d 

        real(wp), allocatable :: var0(:,:) 
        real(wp), allocatable :: dist(:,:) 

        nx = size(var,1)
        ny = size(var,2) 

        allocate(var0(nx,ny)) 
        allocate(dist(nx,ny)) 

        ! Define resolution in km 
        dx_km = dx*1e-3 

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
                        dist(i1,j1) = sqrt( real( (i1-i)**2 + (j1-j)**2 ) ) * dx_km 
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

    function get_opt_param(time,time1,time2,p1,p2,m) result(p) 
        ! Determine value of parameter as a function of time 

        implicit none 

        real(wp), intent(IN) :: time 
        real(wp), intent(IN) :: time1 
        real(wp), intent(IN) :: time2
        real(wp), intent(IN) :: p1
        real(wp), intent(IN) :: p2
        real(wp), intent(IN) :: m         ! Non-linear exponent (m=1.0 or higher)
        real(wp) :: p 

        if (time .le. time1) then 
            p = p1 
        else if (time .ge. time2) then 
            p = p2 
        else 
            ! Linear interpolation with non-linear factor m if desired
            p = p1 + (p2-p1)* ((time-time1)/(time2-time1))**m 
        end if  

        return 

    end function get_opt_param

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

end module ice_optimization
