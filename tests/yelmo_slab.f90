program yelmo_slab
    ! For testing an idealized slab domain

    use ncio 
    use yelmo 
    use deformation 

    implicit none 

    type(yelmo_class)     :: yelmo1

    type ctrl_par 
        character(len=56)   :: domain 
        character(len=56)   :: grid_name
        real(wp)            :: H0 
        real(wp)            :: H_stdev
        real(wp)            :: lx
        real(wp)            :: ly
        real(wp)            :: dx
        real(wp)            :: alpha
        real(wp)            :: visc_const
        real(wp)            :: beta_const

        real(wp)            :: time_init
        real(wp)            :: time_end
        real(wp)            :: dtt
        real(wp)            :: dt1D_out
        real(wp)            :: dt2D_out
    end type 

    type(ctrl_par) :: ctrl 

    character(len=56)  :: domain, grid_name  
    character(len=256) :: outfldr, file2D, file1D, file_restart
    character(len=512) :: path_par, path_const 
    
    integer  :: n  
    real(wp) :: time 
    real(wp) :: xmin, xmax, ymin, ymax 
    real(wp), allocatable :: dh(:,:) 
    
    real(8)  :: cpu_start_time
    real(8)  :: cpu_end_time
    real(8)  :: cpu_dtime
    
    ! Start timing 
    call yelmo_cpu_time(cpu_start_time)
    
    ! Assume program is running from the output folder
    outfldr = "./"

    ! Determine the parameter file from the command line 
    call yelmo_load_command_line_args(path_par)
    !path_par   = trim(outfldr)//"yelmo_TROUGH-F17.nml" 
    
    ! Define input and output locations 
    path_const = trim(outfldr)//"yelmo_const_EISMINT.nml"
    file2D     = trim(outfldr)//"yelmo2D.nc"
    file1D     = trim(outfldr)//"yelmo1D.nc"
    
    ! Define the domain, grid and experiment from parameter file
    call nml_read(path_par,"ctrl","domain",       ctrl%domain)
    call nml_read(path_par,"ctrl","H0",           ctrl%H0)    
    call nml_read(path_par,"ctrl","H_stdev",      ctrl%H_stdev)
    call nml_read(path_par,"ctrl","lx",           ctrl%lx)
    call nml_read(path_par,"ctrl","ly",           ctrl%ly)
    call nml_read(path_par,"ctrl","dx",           ctrl%dx)
    call nml_read(path_par,"ctrl","alpha",        ctrl%alpha)
    call nml_read(path_par,"ctrl","visc_const",   ctrl%visc_const)
    call nml_read(path_par,"ctrl","beta_const",   ctrl%beta_const)

    ! Timing parameters 
    call nml_read(path_par,"ctrl","time_init",    ctrl%time_init)     ! [yr] Starting time
    call nml_read(path_par,"ctrl","time_end",     ctrl%time_end)      ! [yr] Ending time
    call nml_read(path_par,"ctrl","dtt",          ctrl%dtt)           ! [yr] Main loop time step 
    call nml_read(path_par,"ctrl","dt2D_out",     ctrl%dt2D_out)      ! [yr] Frequency of 2D output 
        
    ! Define default grid name for completeness 
    ctrl%grid_name = trim(ctrl%domain)
    ctrl%dt1D_out  = ctrl%dtt 

    ! General initialization of yelmo constants (used globally)
    call yelmo_global_init(path_const)

    ! Define the domain and grid
    xmin =  0.0_prec 
    xmax =  ctrl%lx 
    ymax =  ctrl%ly/2.0_prec 
    ymin = -ctrl%ly/2.0_prec 
    call yelmo_init_grid(yelmo1%grd,grid_name,units="km",x0=xmin,dx=ctrl%dx,nx=int(xmax/ctrl%dx)+1, &
                                                         y0=ymin,dy=ctrl%dx,ny=int((ymax-ymin)/ctrl%dx)+1)

    ! Set time to initial time 

    time = ctrl%time_init 

    ! === Initialize ice sheet model =====

    ! Initialize data objects
    call yelmo_init(yelmo1,filename=path_par,grid_def="none",time=time,load_topo=.FALSE., &
                        domain=domain,grid_name=grid_name)

    ! Load boundary values

    yelmo1%bnd%z_sl     = 0.0
    yelmo1%bnd%bmb_shlf = 0.0 
    yelmo1%bnd%T_shlf   = T0  
    yelmo1%bnd%H_sed    = 0.0 

    yelmo1%bnd%T_srf    = T0                ! [K] 
    yelmo1%bnd%smb      = 0.0_prec          ! [m/yr]
    yelmo1%bnd%Q_geo    = 50.0_prec         ! [mW/m2] 

    ! Check boundary values 
    call yelmo_print_bound(yelmo1%bnd)

    ! Initialize output file 
    call yelmo_write_init(yelmo1,file2D,time_init=time,units="years")
    
    ! ===== Intialize topography =========
    
    yelmo1%bnd%z_bed = 1000.0_wp - sin(ctrl%alpha*degrees_to_radians)*(yelmo1%grd%x*1e-3)

    ! Define initial ice thickness 
    allocate(dh(yelmo1%grd%nx,yelmo1%grd%ny))
    call gen_random_normal(dh,0.0_wp,1.0_wp) 
    yelmo1%tpo%now%H_ice = ctrl%H0 + dh 

    ! Define surface elevation 
    yelmo1%tpo%now%z_srf = yelmo1%bnd%z_bed + yelmo1%tpo%now%H_ice

    ! ====================================

    ! Initialize the yelmo state (dyn,therm,mat)
    call yelmo_init_state(yelmo1,path_par,time=time,thrm_method="robin-cold")

    ! Write initial state 
    call write_step_2D(yelmo1,file2D,time=time) 

    ! 1D file 
    call write_yreg_init(yelmo1,file1D,time_init=time,units="years",mask=yelmo1%bnd%ice_allowed)
    call write_yreg_step(yelmo1,file1D,time=time)  
    
    stop 

    ! Advance timesteps
    do n = 1, ceiling((ctrl%time_end-ctrl%time_init)/ctrl%dtt)

        ! Get current time 
        time = ctrl%time_init + n*ctrl%dtt

        ! == Yelmo ice sheet ===================================================
        call yelmo_update(yelmo1,time)


        ! == MODEL OUTPUT =======================================================
        if (mod(nint(time*100),nint(ctrl%dt2D_out*100))==0) then  
            call write_step_2D(yelmo1,file2D,time=time)    
        end if 

        if (mod(nint(time*100),nint(ctrl%dt1D_out*100))==0) then 
            call write_yreg_step(yelmo1,file1D,time=time) 
        end if

        if (mod(time,10.0)==0 .and. (.not. yelmo_log)) then
            write(*,"(a,f14.4)") "yelmo:: time = ", time
        end if  

    end do 

    ! Write summary 
    write(*,*) "====== "//trim(domain)//" ======="
    write(*,*) "nz, H0 = ", yelmo1%par%nz_aa, maxval(yelmo1%tpo%now%H_ice)

    ! Finalize program
    call yelmo_end(yelmo1,time=time)

    ! Stop timing 
    call yelmo_cpu_time(cpu_end_time,cpu_start_time,cpu_dtime)
    
    write(*,"(a,f12.3,a)") "Time  = ",cpu_dtime/60.0 ," min"
    write(*,"(a,f12.1,a)") "Speed = ",(1e-3*(ctrl%time_end-ctrl%time_init))/(cpu_dtime/3600.0), " kiloyears / hr"
    
contains
    
    subroutine write_step_2D(ylmo,filename,time)

        implicit none 
        
        type(yelmo_class), intent(IN) :: ylmo
        character(len=*),  intent(IN) :: filename
        real(prec), intent(IN) :: time

        ! Local variables
        integer    :: ncid, n, i, j, nx, ny  
        real(prec) :: time_prev 

        nx = ylmo%tpo%par%nx 
        ny = ylmo%tpo%par%ny 

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_size(filename,"time",ncid)
        call nc_read(filename,"time",time_prev,start=[n],count=[1],ncid=ncid) 
        if (abs(time-time_prev).gt.1e-5) n = n+1 

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        ! Write model metrics (model speed, dt, eta)
        call yelmo_write_step_model_metrics(filename,ylmo,n,ncid)

        ! == yelmo_topography ==
        call nc_write(filename,"H_ice",ylmo%tpo%now%H_ice,units="m",long_name="Ice thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf",ylmo%tpo%now%z_srf,units="m",long_name="Surface elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_bed",ylmo%tpo%now%mask_bed,units="",long_name="Bed mask", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"N_eff",ylmo%dyn%now%N_eff,units="Pa",long_name="Effective pressure", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"dzsrfdt",ylmo%tpo%now%dzsrfdt,units="m/a",long_name="Surface elevation change", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dHicedt",ylmo%tpo%now%dHicedt,units="m/a",long_name="Ice thickness change", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"H_grnd",ylmo%tpo%now%H_grnd,units="m",long_name="Ice thickness overburden", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"bmb",ylmo%tpo%now%bmb,units="m/a",long_name="Basal mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! == yelmo_thermodynamics ==
        call nc_write(filename,"T_ice",ylmo%thrm%now%T_ice,units="K",long_name="Ice temperature", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"T_prime",ylmo%thrm%now%T_ice-ylmo%thrm%now%T_pmp,units="deg C",long_name="Homologous ice temperature", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"f_pmp",ylmo%thrm%now%f_pmp,units="1",long_name="Fraction of grid point at pmp", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"T_prime_b",ylmo%thrm%now%T_prime_b,units="deg C",long_name="Homologous basal ice temperature", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"H_w",ylmo%thrm%now%H_w,units="m",long_name="Basal water layer", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! == yelmo_material ==
!         call nc_write(filename,"visc_int",ylmo%mat%now%visc_int,units="Pa a m",long_name="Vertically integrated viscosity", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"visc",ylmo%mat%now%visc,units="Pa a",long_name="Viscosity", &
!                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"ATT_bar",ylmo%mat%now%ATT_bar,units="a^-1 Pa^-3",long_name="Vertically averaged rate factor", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"ATT",ylmo%mat%now%ATT,units="a^-1 Pa^-3",long_name="Rate factor", &
!                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        
        call nc_write(filename,"Q_ice_b",ylmo%thrm%now%Q_ice_b,units="J a-1 m-2",long_name="Basal ice heat flux", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"Q_strn",ylmo%thrm%now%Q_strn/(rho_ice*ylmo%thrm%now%cp),units="K a-1",long_name="Strain heating", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        call nc_write(filename,"Q_b",ylmo%thrm%now%Q_b,units="J a-1 m-2",long_name="Basal frictional heating", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! == yelmo_dynamics ==

        call nc_write(filename,"cf_ref",ylmo%dyn%now%cf_ref,units="--",long_name="Bed friction scalar", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"c_bed",ylmo%dyn%now%c_bed,units="Pa",long_name="Bed friction coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"beta",ylmo%dyn%now%beta,units="Pa a m-1",long_name="Basal friction coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"beta_acx",ylmo%dyn%now%beta_acx,units="Pa a m-1",long_name="Basal friction coefficient (acx)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"beta_acy",ylmo%dyn%now%beta_acy,units="Pa a m-1",long_name="Basal friction coefficient (acy)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"visc_eff_int",ylmo%dyn%now%visc_eff_int,units="Pa a m",long_name="Depth-integrated effective viscosity (SSA)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"taud_acx",ylmo%dyn%now%taud_acx,units="Pa",long_name="Driving stress, x-direction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taud_acy",ylmo%dyn%now%taud_acy,units="Pa",long_name="Driving stress, y-direction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"taub_acx",ylmo%dyn%now%taub_acx,units="Pa",long_name="Basal stress, x-direction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taub_acy",ylmo%dyn%now%taub_acy,units="Pa",long_name="Basal stress, y-direction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"ux_bar",ylmo%dyn%now%ux_bar,units="m/a",long_name="Vertically averaged velocity (x)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy_bar",ylmo%dyn%now%uy_bar,units="m/a",long_name="Vertically averaged velocity (y)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_bar",ylmo%dyn%now%uxy_bar,units="m/a",long_name="Vertically averaged velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"de",ylmo%mat%now%strn%de,units="a^-1",long_name="Vertically averaged strain rate", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        ! No time dimension::

        call nc_write(filename,"z_bed",ylmo%bnd%z_bed,units="m",long_name="Bedrock elevation", &
                      dim1="xc",dim2="yc",start=[1,1],ncid=ncid)
        
        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine write_step_2D

    subroutine gen_random_normal(ynrm,mu,sigma)
        ! Calculate a random number from a normal distribution 
        ! following the Box-Mueller algorithm 
        ! https://en.wikipedia.org/wiki/Normal_distribution#Generating_values_from_normal_distribution 

        implicit none 

        real(wp), intent(OUT) :: ynrm(:,:) 
        real(wp), intent(IN)  :: mu 
        real(wp), intent(IN)  :: sigma 

        ! Local variables 
        integer :: i, j, nx, ny 
        real(wp) :: yuni(2)
        real(wp), parameter :: pi = 4.0_wp*atan(1.0_wp)

        nx = size(ynrm,1)
        ny = size(ynrm,2) 

        do j = 1, ny 
        do i = 1, nx 

            ! Get 2 numbers from Uniform distribution between 0 and 1
            call random_number(yuni)
            
            ! Convert to normal distribution using the Box-Mueller algorithm
            ynrm(i,j) = mu + sigma * sqrt(-2.0*log(yuni(1))) * cos(2*pi*yuni(2))

        end do 
        end do 

        return 

    end subroutine gen_random_normal

end program yelmo_slab



