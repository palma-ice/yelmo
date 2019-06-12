

program yelmo_test

    use ncio 
    use yelmo 
    
    use gaussian_filter 

    implicit none 

    type(yelmo_class)      :: yelmo1

    character(len=256) :: outfldr, file1D, file2D, file_restart, domain 
    character(len=512) :: path_par, path_const  
    real(prec) :: time_init, time_end, time_equil, time, dtt, dt1D_out, dt2D_out   
    integer    :: n
    real(4) :: cpu_start_time, cpu_end_time 

    ! Optimization variables 
    real(prec) :: time_iter 
    integer    :: q, qmax, qmax_topo_fixed 
    logical    :: topo_fixed 
    real(prec) :: phi_min, phi_max 

    real(prec), allocatable :: dCbed(:,:) 
    real(prec), allocatable :: phi(:,:) 

    ! No-ice mask (to impose additional melting)
    logical, allocatable :: mask_noice(:,:)  

    ! Concavity field 
    real(prec), allocatable :: channels(:,:) 

    ! Start timing 
    call cpu_time(cpu_start_time)

    ! Determine the parameter file from the command line 
    call yelmo_load_command_line_args(path_par)

    ! Assume program is running from the output folder
    outfldr = "./"

    ! Define input and output locations 
    path_const = trim(outfldr)//"yelmo_const_Earth.nml"
    file1D     = trim(outfldr)//"yelmo1D.nc"
    file2D     = trim(outfldr)//"yelmo2D.nc"

    ! === Initialize ice sheet model =====

    ! General initialization of yelmo constants (used globally)
    call yelmo_global_init(path_const)

    ! Initialize data objects and load initial topography
    call yelmo_init(yelmo1,filename=path_par,grid_def="file",time=time_init)

    ! Simulation parameters
    time_init       = 0.0           ! [yr] Starting time
    time_iter       = 500.0         ! [yr] Simulation time for each iteration
    qmax            = 200           ! Total number of iterations
    qmax_topo_fixed = 0             ! Number of initial iterations that should use topo_fixed=.TRUE. 

    phi_min         =  2.0 
    phi_max         = 40.0 

    ! Prescribe key parameters here that should be set for beta optimization exercise 
    yelmo1%dyn%par%C_bed_method      = -1       ! C_Bed is set external to yelmo calculations
    yelmo1%mat%par%rf_method         = 0        ! Constant rate factor (no thermodynamics)
    yelmo1%mat%par%rf_const          = 1e-17    ! [Pa^-3 a^-1]
!     yelmo1%thrm%par%method           = "fixed"  ! No thermodynamics calculations 

    yelmo1%dyn%par%beta_method       = 3        ! 0: constant beta; 1: power law; 2: effective pressure; 3: regularized Coulomb law, 4: pism power law
    yelmo1%dyn%par%m_drag            = 3.0      ! Dragging law exponent 
    yelmo1%dyn%par%cf_stream         = 0.1      ! [--] Friction scalar, unitless for reg. coulomb law  
    
    yelmo1%dyn%par%neff_method       = 2        ! -1: external N_eff, 0: overburden pressure, 1: Leguy param., 2: Till pressure
    
    ! === Set initial boundary conditions for current time and yelmo state =====
    ! ybound: z_bed, z_sl, H_sed, H_w, smb, T_srf, bmb_shlf , Q_geo

    yelmo1%bnd%z_sl     = 0.0           ! [m]
    yelmo1%bnd%H_sed    = 0.0           ! [m]
    yelmo1%bnd%H_w      = 0.0           ! [m]
    yelmo1%bnd%Q_geo    = 50.0          ! [mW/m2]
    
    yelmo1%bnd%bmb_shlf = -20.0         ! [m.i.e./a]
    yelmo1%bnd%T_shlf   = T0            ! [K]   

    ! Impose present-day surface mass balance and present-day temperature field
    yelmo1%bnd%smb      = yelmo1%dta%pd%smb        ! [m.i.e./a]
    yelmo1%bnd%T_srf    = yelmo1%dta%pd%T_srf      ! [K]
    
    call yelmo_print_bound(yelmo1%bnd)

    time = time_init 

    ! Allocate rate of change of C_bed
    allocate(dCbed(yelmo1%grd%nx,yelmo1%grd%ny))
    dCbed = 0.0 

    allocate(phi(yelmo1%grd%nx,yelmo1%grd%ny))
    phi = phi_min  

    ! Define channel field 
    allocate(channels(yelmo1%grd%nx,yelmo1%grd%ny))

    ! Set initial guess of C_bed as a function of present-day velocity 
    call guess_C_bed(yelmo1%dyn%now%C_bed,phi,yelmo1%dta%pd%uxy_s,phi_min,phi_max,yelmo1%dyn%par%cf_stream)

    ! Initialize state variables (dyn,therm,mat)
    ! (initialize temps with robin method with a cold base)
    call yelmo_init_state(yelmo1,path_par,time=time_init,thrm_method="robin-cold")

    ! Define no-ice mask from present-day data
    allocate(mask_noice(yelmo1%grd%nx,yelmo1%grd%ny))
    mask_noice = .FALSE. 
    where(yelmo1%dta%pd%H_ice .le. 0.0) mask_noice = .TRUE. 

    ! Impose additional negative mass balance to no ice points of 2 [m.i.e./a] melting
    where(mask_noice) yelmo1%bnd%smb = yelmo1%dta%pd%smb - 2.0 

    ! Run yelmo for several years with constant boundary conditions and topo
    ! to equilibrate thermodynamics and dynamics
    call yelmo_update_equil(yelmo1,time,time_tot=100.0,topo_fixed=.FALSE.,dt=1.0,ssa_vel_max=500.0)
    
    ! Define present topo as present-day dataset for comparison 
    yelmo1%dta%pd%H_ice = yelmo1%tpo%now%H_ice 
    yelmo1%dta%pd%z_srf = yelmo1%tpo%now%z_srf 

    ! Calculate initial error 
    yelmo1%dta%pd%err_H_ice   = yelmo1%tpo%now%H_ice - yelmo1%dta%pd%H_ice 
    yelmo1%dta%pd%err_z_srf   = yelmo1%tpo%now%z_srf - yelmo1%dta%pd%z_srf 
    yelmo1%dta%pd%err_uxy_s   = yelmo1%dyn%now%uxy_s - yelmo1%dta%pd%uxy_s 
    
    ! 2D file 
    call yelmo_write_init(yelmo1,file2D,time_init=0.0,units="years")  
    call write_step_2D_opt(yelmo1,file2D,time=0.0,dCbed=dCbed,phi=phi)  

    ! Initially assume we are working with topo_fixed... (only for optimizing velocity)
    topo_fixed = .TRUE. 
    
    ! Perform loops over beta:
    ! update beta, calculate topography and velocity for 100 years, get error, try again
    do q = 1, qmax 

        ! Determine whether this iteration maintains topo_fixed conditions (only for optimizing velocity)
        if (q .gt. qmax_topo_fixed) topo_fixed = .FALSE. 

        ! Reset topography to initial state 
        yelmo1%tpo%now%H_ice = yelmo1%dta%pd%H_ice 
        yelmo1%tpo%now%z_srf = yelmo1%dta%pd%z_srf 

        ! Update C_bed based on error correction
        call update_C_bed_thickness(yelmo1%dyn%now%C_bed,dCbed,phi,yelmo1%dta%pd%err_z_srf,yelmo1%tpo%now%H_ice, &
                    yelmo1%dyn%now%ux_bar,yelmo1%dyn%now%uy_bar,yelmo1%tpo%par%dx,phi_min,phi_max,yelmo1%dyn%par%cf_stream)

        ! Run model for time_iter yrs with this C_bed configuration (no change in boundaries)
        call yelmo_update_equil(yelmo1,time,time_tot=time_iter,topo_fixed=topo_fixed,dt=5.0,ssa_vel_max=5000.0)
        
        ! Calculate error 
        yelmo1%dta%pd%err_H_ice   = yelmo1%tpo%now%H_ice - yelmo1%dta%pd%H_ice 
        yelmo1%dta%pd%err_z_srf   = yelmo1%tpo%now%z_srf - yelmo1%dta%pd%z_srf 
        yelmo1%dta%pd%err_uxy_s = yelmo1%dyn%now%uxy_s   - yelmo1%dta%pd%uxy_s 
        
        ! == MODEL OUTPUT =======================================================

        time = real(q,prec)
        
        call write_step_2D_opt(yelmo1,file2D,time=time,dCbed=dCbed,phi=phi)
        
        ! Summary 
        write(*,*) "q= ", q, maxval(abs(yelmo1%dta%pd%err_z_srf))

    end do 

    ! Finalize program
    call yelmo_end(yelmo1,time=time)

    ! Stop timing 
    call cpu_time(cpu_end_time)

    write(*,"(a,f12.3,a)") "Time  = ",(cpu_end_time-cpu_start_time)/60.0 ," min"
    write(*,"(a,f12.1,a)") "Speed = ",(1e-3*(time_end-time_init))/((cpu_end_time-cpu_start_time)/3600.0), " kiloyears / hr"

contains

    subroutine write_step_2D_opt(ylmo,filename,time,dCbed,phi)

        implicit none 
        
        type(yelmo_class), intent(IN) :: ylmo
        character(len=*),  intent(IN) :: filename
        real(prec), intent(IN) :: time
        real(prec), intent(IN) :: dCbed(:,:) 
        real(prec), intent(IN) :: phi(:,:)

        ! Local variables
        integer    :: ncid, n
        real(prec) :: time_prev 

        real(prec), allocatable :: tmp(:,:) 
        real(prec) :: rmse, err, npts  

        allocate(tmp(ylmo%grd%nx,ylmo%grd%ny))

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_size(filename,"time",ncid)
        call nc_read(filename,"time",time_prev,start=[n],count=[1],ncid=ncid) 
        if (abs(time-time_prev).gt.1e-5) n = n+1 

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        call nc_write(filename,"H_ice",ylmo%tpo%now%H_ice,units="m",long_name="Ice thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf",ylmo%tpo%now%z_srf,units="m",long_name="Surface elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_bed",ylmo%tpo%now%mask_bed,units="",long_name="Bed mask", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"C_bed",ylmo%dyn%now%C_bed,units="",long_name="Bed constant", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dCbed",dCbed,units="-",long_name="Bed constant change", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"phi",phi,units="degrees",long_name="Till friction angle", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"N_eff",ylmo%dyn%now%N_eff,units="Pa",long_name="Effective pressure", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"beta",ylmo%dyn%now%beta,units="Pa a m-1",long_name="Basal friction coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"visc_eff",ylmo%dyn%now%visc_eff,units="1",long_name="Vertically averaged viscosity", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"uxy_s",ylmo%dyn%now%uxy_s,units="m/a",long_name="Surface velocity magnitude", &
                     dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_b",ylmo%dyn%now%uxy_b,units="m/a",long_name="Basal sliding velocity magnitude", &
                     dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"f_pmp",ylmo%thrm%now%f_pmp,units="",long_name="Basal temperate fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! Boundary variables (forcing)
        call nc_write(filename,"z_bed",ylmo%bnd%z_bed,units="m",long_name="Bedrock elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"smb",ylmo%bnd%smb,units="m a-1",long_name="Annual surface mass balance (ice equiv.)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"T_srf",ylmo%bnd%T_srf,units="K",long_name="Surface temperature", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"H_w",ylmo%bnd%H_w,units="m",long_name="Basal water layer", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! Target data (not time dependent)
        if (n .eq. 1) then 
            call nc_write(filename,"pd_z_srf",ylmo%dta%pd%z_srf,units="m",long_name="Observed surface elevation (present day)", &
                          dim1="xc",dim2="yc",ncid=ncid)
            call nc_write(filename,"pd_uxy_s",ylmo%dta%pd%uxy_s,units="m",long_name="Observed surface velocity (present day)", &
                          dim1="xc",dim2="yc",ncid=ncid)
        end if 

        ! Error fields compared to targets
        call nc_write(filename,"pd_err_z_srf",ylmo%dta%pd%err_z_srf,units="m",long_name="Surface elevation error (present day)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"pd_err_uxy_s",ylmo%dta%pd%err_uxy_s,units="m",long_name="Surface velocity error (present day)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        tmp = ylmo%dta%pd%err_z_srf
        call filter_gaussian(var=tmp,sigma=2.0*ylmo%tpo%par%dx,dx=ylmo%tpo%par%dx, &
                                mask=ylmo%dta%pd%err_z_srf .ne. 0.0)
        
        call nc_write(filename,"pd_err_z_srf_sm",tmp,units="m",long_name="Smooth surface elevation error (present day)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! Get cumulative error diagnostics (z_srf)
        npts = ylmo%tpo%par%nx*ylmo%tpo%par%ny 
        err  = sum(ylmo%dta%pd%err_z_srf) / real(npts,prec)
        rmse = sqrt( sum(ylmo%dta%pd%err_z_srf**2) / real(npts,prec) ) 

        call nc_write(filename,"err_z_srf",err,units="m",long_name="Mean surface elevation error (present day)", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"rmse_z_srf",rmse,units="m",long_name="Root mean square error (present day z_srf)", &
                      dim1="time",start=[n],ncid=ncid)
        
        ! Get cumulative error diagnostics (uxy_s)
        npts = ylmo%tpo%par%nx*ylmo%tpo%par%ny 
        err  = sum(ylmo%dta%pd%err_uxy_s) / real(npts,prec)
        rmse = sqrt( sum(ylmo%dta%pd%err_uxy_s**2) / real(npts,prec) ) 

        call nc_write(filename,"err_uxy_s",err,units="m a^-1",long_name="Mean surface velocity error (present day)", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"rmse_uxy_s",rmse,units="m a^-1",long_name="Root mean square error (present day uxy_s)", &
                      dim1="time",start=[n],ncid=ncid)
        
        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine write_step_2D_opt

    subroutine guess_C_bed(C_bed,phi,uxy_s,phi_min,phi_max,cf_stream)

        implicit none 

        real(prec), intent(INOUT) :: C_bed(:,:) 
        real(prec), intent(INOUT) :: phi(:,:) 
        real(prec), intent(IN)    :: uxy_s(:,:) 
        real(prec), intent(IN)    :: phi_min
        real(prec), intent(IN)    :: phi_max
        real(prec), intent(IN)    :: cf_stream 

        ! Local variables 
        integer :: i, j, nx, ny 
        real(prec) :: logvel, logvel_max, f_scale    

        nx = size(C_bed,1)
        ny = size(C_bed,2)

        logvel_max = log(maxval(uxy_s))

        do j = 1, ny 
        do i = 1, nx 

            if (uxy_s(i,j) .gt. 0.0) then 
                ! Calculate expected till angle versus velocity 

                logvel   = max(0.0,log(uxy_s(i,j)))
                f_scale  = logvel / logvel_max
                phi(i,j) = phi_max - f_scale*(phi_max-phi_min)

            else 
                ! Set phi to the minimum 

                phi(i,j) = phi_min 

            end if 

            ! Calculate C_bed following till friction approach (Bueler and van Pelt, 2015)
            C_bed(i,j) = cf_stream*tan(phi(i,j)*pi/180)

        end do 
        end do

        return 

    end subroutine guess_C_bed

    subroutine update_C_bed_thickness(C_bed,dCbed,phi,err_z_srf,H_ice,ux,uy,dx,phi_min,phi_max,cf_stream)

        implicit none 

        real(prec), intent(INOUT) :: C_bed(:,:) 
        real(prec), intent(INOUT) :: dCbed(:,:) 
        real(prec), intent(INOUT) :: phi(:,:) 
        real(prec), intent(IN)    :: err_z_srf(:,:) 
        real(prec), intent(IN)    :: H_ice(:,:) 
        real(prec), intent(IN)    :: ux(:,:) 
        real(prec), intent(IN)    :: uy(:,:) 
        real(prec), intent(IN)    :: dx 
        real(prec), intent(IN)    :: phi_min 
        real(prec), intent(IN)    :: phi_max 
        real(prec), intent(IN)    :: cf_stream 

        ! Local variables 
        integer :: i, j, nx, ny, i1, j1 
        real(prec) :: dphi, f_dz
        real(prec) :: ux_aa, uy_aa  
        real(prec), allocatable   :: C_bed_prev(:,:) 

        real(prec), parameter :: dphi_min  = -0.5       ! [degrees]
        real(prec), parameter :: dphi_max  =  1.0       ! [degrees]
        real(prec), parameter :: err_z_fac = 100.0      ! [m] 

        nx = size(C_bed,1)
        ny = size(C_bed,2) 

        allocate(C_bed_prev(nx,ny))

        ! Store initial C_bed solution 
        C_bed_prev = C_bed 

        do j = 2, ny-1 
        do i = 2, nx-1 

            if (err_z_srf(i,j) .ne. 0.0) then 
                ! Update where elevation error exists

                ! Get adjustment rate given error in z_srf 
                f_dz = -err_z_srf(i,j) / err_z_fac 
                f_dz = max(f_dz, dphi_min)
                f_dz = min(f_dz, dphi_max)
                dphi = f_dz 
                
                ! 1. Apply change at current point 
if (.FALSE.) then 
                phi(i,j)  = phi(i,j) + dphi 
                phi(i,j)  = max(phi(i,j),phi_min)
                phi(i,j)  = min(phi(i,j),phi_max)

                C_bed(i,j) = cf_stream*tan(phi(i,j)*pi/180.0)
end if 

                ! 2. Apply change downstream (this may overlap with other changes)

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

                phi(i1,j1)  = phi(i1,j1) + dphi 
                phi(i1,j1)  = max(phi(i1,j1),phi_min)
                phi(i1,j1)  = min(phi(i1,j1),phi_max)

                C_bed(i1,j1) = cf_stream*tan(phi(i1,j1)*pi/180.0)

            end if 

        end do 
        end do 

        ! Additionally, apply a Gaussian filter to C_bed to ensure smooth transitions 
!         call filter_gaussian(var=C_bed,sigma=64.0,dx=dx) !, &
                                !mask=err_z_srf .ne. 0.0)
        
        dCbed = C_bed - C_bed_prev

        return 

    end subroutine update_C_bed_thickness

    ! Extra...

    subroutine calc_ydyn_cbed_external(dyn,tpo,thrm,bnd,channels)
        ! Update C_bed based on parameter choices

        implicit none
        
        type(ydyn_class),   intent(INOUT) :: dyn
        type(ytopo_class),  intent(IN)    :: tpo 
        type(ytherm_class), intent(IN)    :: thrm
        type(ybound_class), intent(IN)    :: bnd  
        real(prec),         intent(INOUT) :: channels(:,:) 

        integer :: i, j, nx, ny 
        integer :: i1, i2, j1, j2 
        real(prec), allocatable :: f_channel(:,:) 

        real(prec) :: channel_lim = 1e-6 

        nx = size(dyn%now%C_bed,1)
        ny = size(dyn%now%C_bed,2)
        
        allocate(f_channel(nx,ny)) 

        ! Set C_bed according to temperate character of base

        ! Smooth transition between temperate and frozen C_bed
        dyn%now%C_bed = (thrm%now%f_pmp)*dyn%par%cf_stream &
                    + (1.0_prec - thrm%now%f_pmp)*dyn%par%cf_frozen 

        if (dyn%par%streaming_margin) then 
            ! Ensure that both the margin points and the grounding line
            ! are always considered streaming, independent of their
            ! thermodynamic character (as sometimes these can incorrectly become frozen)

        
            ! Ensure any marginal point is also treated as streaming 
            do j = 1, ny 
            do i = 1, nx 

                i1 = max(i-1,1)
                i2 = min(i+1,nx)
                j1 = max(j-1,1)
                j2 = min(j+1,ny)

                if (tpo%now%H_ice(i,j) .gt. 0.0 .and. &
                    (tpo%now%H_ice(i1,j) .le. 0.0 .or. &
                     tpo%now%H_ice(i2,j) .le. 0.0 .or. &
                     tpo%now%H_ice(i,j1) .le. 0.0 .or. &
                     tpo%now%H_ice(i,j2) .le. 0.0)) then 

                    dyn%now%C_bed(i,j) = dyn%par%cf_stream

                end if 

            end do 
            end do 

            ! Also ensure that grounding line is also considered streaming
            where(tpo%now%is_grline) dyn%now%C_bed = dyn%par%cf_stream

        end if 

        ! == Until here, C_bed is defined as normally with C_bed_method=1,
        !    now refine to increase only marginal velocities 

        ! Reduce C_bed further for low elevation points
        !where(tpo%now%z_srf .lt. 1500.0) dyn%now%C_bed = 0.5*dyn%now%C_bed

        ! Next diagnose channels
        call calc_channels(channels,tpo%now%z_srf,dyn%now%ux_bar,dyn%now%uy_bar,tpo%par%dx)

        ! Finally scale C_bed according to concavity of channels 
        !f_channel = exp(-channels/channel_lim)
        !where(f_channel .lt. 0.1) f_channel = 0.1 
        !where(f_channel .gt. 2.0) f_channel = 2.0  

        f_channel = 1.0 

        dyn%now%C_bed = dyn%now%C_bed * f_channel 
        
        return 

    end subroutine calc_ydyn_cbed_external


    subroutine calc_channels(channels,z_bed,ux,uy,dx)

        implicit none 

        real(prec), intent(OUT) :: channels(:,:) 
        real(prec), intent(IN)  :: z_bed(:,:) 
        real(prec), intent(IN)  :: ux(:,:) 
        real(prec), intent(IN)  :: uy(:,:) 
        real(prec), intent(IN)  :: dx 

        ! Local variables 
        integer :: i, j, nx, ny 
        real(prec) :: ux_aa, uy_aa, uxy 
        real(prec) :: dzdx1, dzdx2 
        real(prec) :: dz2dx2, dz2dy2
        real(prec) :: theta   ! [rad] Angle of direction of ice flow
        real(prec) :: alpha   ! [rad] Angle of direction perpindicular to ice flow 

        ! Finite-difference coefficients
        real(prec), parameter :: fm2 = -1.0/12.0 
        real(prec), parameter :: fm1 =  4.0/3.0 
        real(prec), parameter :: f0  = -5.0/2.0 
        real(prec), parameter :: fp1 =  4.0/3.0
        real(prec), parameter :: fp2 = -1.0/12.0 
         
        
        
        nx = size(channels,1)
        ny = size(channels,2)

        ! Set channels to zero initially 
        channels = 0.0 

        ! Find channels based on change in elevation perpendicular to flow direction,
        ! then (to do!) negative component for along flow direction  
        do j = 3, ny-2 
        do i = 3, nx-2 

            ! Get velocity of current grid point 
            ux_aa = 0.5*(ux(i-1,j) + ux(i,j))
            uy_aa = 0.5*(uy(i,j-1) + uy(i,j))
            uxy   = sqrt(ux_aa**2 + uy_aa**2)

            if (uxy .gt. 0.0) then 

                ! Get direction perpindicular ice flow 
                alpha = atan2(uy_aa,ux_aa) - pi/2.0 

                ! Only modify areas with some velocity available 

                ! Calculate second-derivative in each direction (2nd order)
                dz2dx2 = (1.0/dx**2)*sum([fm2,fm1,f0,fp1,fp2]*z_bed(i-2:i+2,j))
                dz2dy2 = (1.0/dx**2)*sum([fm2,fm1,f0,fp1,fp2]*z_bed(i,j-2:j+2))
                
                ! Scale derivative in each direction to get approximate concavity in
                ! direction of interest 
                channels(i,j) = cos(alpha)*dz2dx2 + sin(alpha)*dz2dy2

!                 if (abs(ux_aa) .gt. abs(uy_aa)) then 
!                     ! Flow predominantly in x-direction

!                     dzdx1         = (z_bed(i,j)   - z_bed(i,j-1)) / dx 
!                     dzdx2         = (z_bed(i,j+1) - z_bed(i,j))   / dx 
!                     channels(i,j) = (dzdx2-dzdx1) / dx 

!                     !channels(i,j) = (0.5*(z_bed(i,j-1)+z_bed(i,j+1)) - z_bed(i,j)) / dx 

!                 else 
!                     ! Flow predominantly in y-direction 

!                     dzdx1         = (z_bed(i,j)   - z_bed(i-1,j)) / dx 
!                     dzdx2         = (z_bed(i+1,j) - z_bed(i,j))   / dx 
!                     channels(i,j) = (dzdx2-dzdx1) / dx 
                    
!                     !channels(i,j) = (0.5*(z_bed(i-1,j)+z_bed(i+1,j)) - z_bed(i,j)) / dx 

!                 end if 


            end if 

        end do 
        end do 

        return 

    end subroutine calc_channels

end program yelmo_test 



