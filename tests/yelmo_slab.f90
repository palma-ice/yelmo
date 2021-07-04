program yelmo_slab
    ! For testing an idealized slab domain

    use ncio 
    use yelmo 
    use deformation 

    implicit none 

    type(yelmo_class)     :: yelmo1

    type ctrl_par 
        character(len=56)       :: domain 
        character(len=56)       :: grid_name
        real(wp)                :: H0 
        real(wp)                :: H_stdev
        integer                 :: nx
        integer                 :: ny
        real(wp)                :: dx
        real(wp)                :: alpha

        real(wp)                :: time_init
        ! real(wp)                :: time_end
        integer                 :: nt 
        real(wp)                :: dtt
        real(wp)                :: dt1D_out
        real(wp)                :: dt2D_out

        integer                 :: n_dtt 
        real(wp), allocatable   :: dtts(:) 
        integer                 :: n_dx 
        real(wp), allocatable   :: dxs(:) 

        character(len=512)  :: path_par 
        
    end type 

    type(ctrl_par) :: ctrl 

    character(len=256) :: outfldr, file2D, file1D, file_restart
    character(len=512) :: path_par, path_const 
    
    integer  :: n, j
    real(wp) :: time 
    real(wp) :: xmin, xmax, ymin, ymax 
    real(wp), allocatable :: dh(:,:) 
    
    integer  :: q, q1, q2  
    real(wp) :: dtt_now
    real(wp) :: stdev, factor 

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
    call nml_read(path_par,"ctrl","nx",           ctrl%nx)
    call nml_read(path_par,"ctrl","ny",           ctrl%ny)
    call nml_read(path_par,"ctrl","dx",           ctrl%dx)
    call nml_read(path_par,"ctrl","alpha",        ctrl%alpha)

    ! Timing parameters 
    call nml_read(path_par,"ctrl","time_init",    ctrl%time_init)     ! [yr] Starting time
    !call nml_read(path_par,"ctrl","time_end",     ctrl%time_end)      ! [yr] Ending time
    call nml_read(path_par,"ctrl","nt",           ctrl%nt)            ! [--] Total timesteps to run
    call nml_read(path_par,"ctrl","dtt",          ctrl%dtt)           ! [yr] Main loop time step 
    call nml_read(path_par,"ctrl","dt1D_out",     ctrl%dt1D_out)      ! [yr] Frequency of 2D output 
    call nml_read(path_par,"ctrl","dt2D_out",     ctrl%dt2D_out)      ! [yr] Frequency of 2D output 
    
    ! Set ctrl parameters for later use 
    ctrl%path_par   = path_par 
    
    ! Define default grid name for completeness 
    ctrl%grid_name = trim(ctrl%domain) 

    ! Set time to initial time 
    time = ctrl%time_init 

    ! General initialization of yelmo constants (used globally)
    call yelmo_global_init(path_const)

    ! === Initialize ice sheet model =====

    ! Define the domain and grid
    xmin =  0.0_prec 
    xmax =  (ctrl%nx-1)*ctrl%dx  
    ymax =  (ctrl%ny-1)*ctrl%dx/2.0_prec 
    ymin = -ymax
    call yelmo_init_grid(yelmo1%grd,ctrl%grid_name,units="km", &
                         x0=xmin,dx=ctrl%dx,nx=ctrl%nx, &
                         y0=ymin,dy=ctrl%dx,ny=ctrl%ny)

    ! Initialize data objects
    call yelmo_init(yelmo1,filename=path_par,grid_def="none",time=time,load_topo=.FALSE., &
                        domain=ctrl%domain,grid_name=ctrl%grid_name)

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
    
    ! ===== Intialize topography and set parameters =========
    
    yelmo1%bnd%z_bed = 1000.0_wp - ctrl%alpha*(yelmo1%grd%x)

    ! Define initial ice thickness 
    allocate(dh(yelmo1%grd%nx,yelmo1%grd%ny))
    call gen_random_normal(dh,0.0_wp,ctrl%H_stdev) 
    yelmo1%tpo%now%H_ice = ctrl%H0 + dh 

    ! Make sure all values are the same in y-direction
    do j = 1, yelmo1%grd%ny 
        yelmo1%tpo%now%H_ice(:,j) = yelmo1%tpo%now%H_ice(:,1)
    end do 

    ! Define surface elevation 
    yelmo1%tpo%now%z_srf = yelmo1%bnd%z_bed + yelmo1%tpo%now%H_ice

    ! Define reference ice thickness (for prescribing boundary values, potentially)
    yelmo1%bnd%H_ice_ref = ctrl%H0 

    ! =======================================================

    ! Initialize the yelmo state (dyn,therm,mat)
    call yelmo_init_state(yelmo1,time=time,thrm_method="robin-cold")

    ! Write initial state 
    call write_step_2D(yelmo1,file2D,time=time) 

    ! ! 1D file 
    ! call yelmo_write_reg_init(yelmo1,file1D,time_init=time,units="years",mask=yelmo1%bnd%ice_allowed)
    ! call yelmo_write_reg_step(yelmo1,file1D,time=time)  


if (ctrl%dtt .ne. 0.0) then 
    ! == Perform one simulation with an outer timestep of ctrl%dtt =======

    ! Advance timesteps
    do n = 1, ctrl%nt 

        ! Get current time 
        time = ctrl%time_init + n*ctrl%dtt

        ! == Yelmo ice sheet ===================================================
        call yelmo_update(yelmo1,time)

        if (mod(time,10.0)==0 .and. (.not. yelmo_log)) then
            write(*,"(a,f14.4)") "yelmo:: time = ", time
        end if  

    end do 

    ! Write final state 
    call write_step_2D(yelmo1,file2D,time=time) 

    ! Calculate summary 
    call calc_stdev(stdev,yelmo1%tpo%now%H_ice)
    factor = stdev / max(ctrl%H_stdev,1e-5)

    ! Write summary 
    write(*,*) "====== "//trim(ctrl%domain)//" ======="
    write(*,*) "factor", ctrl%dx, ctrl%dtt, ctrl%H_stdev, stdev, factor 

    ! Finalize program
    call yelmo_end(yelmo1,time=time)

    ! Stop timing 
    call yelmo_cpu_time(cpu_end_time,cpu_start_time,cpu_dtime)
    
    write(*,"(a,f12.3,a)") "Time  = ",cpu_dtime/60.0 ," min"
    write(*,"(a,f12.1,a)") "Speed = ",(1e-3*(time-ctrl%time_init))/(cpu_dtime/3600.0), " kiloyears / hr"
    
else 
    ! === Perform ensemble of simulations with multiple values of dx and dt =====

    ! Test different timesteps 
    
    ctrl%n_dtt = 60
    allocate(ctrl%dtts(ctrl%n_dtt)) 
    ctrl%dtts  = 0.0_wp
    do q = 1, ctrl%n_dtt 
        ctrl%dtts(q) = 10.0_dp**(log10(0.001_dp)+(log10(10.0_dp)-log10(0.001_dp))*(q-1)/real(ctrl%n_dtt-1,dp))
    end do

    ctrl%n_dx     = 12
    allocate(ctrl%dxs(ctrl%n_dx)) 
    ctrl%dxs      = 0.0_wp 
    ctrl%dxs(1:8) = [0.01_wp,0.025_wp,0.05_wp,0.1_wp,0.25_wp,0.5_wp,1.0_wp,2.5_wp,5.0_wp,10.0_wp,25.0_wp,40.0_wp]

    ! ctrl%n_dtt = 20 
    ! ctrl%dtts  = 0.0_wp
    ! do q = 1, ctrl%n_dtt 
    !     ctrl%dtts(q) = exp(log(0.005)+(log(0.3)-log(0.005))*(q-1)/(ctrl%n_dtt-1))
    ! end do 

    ! ctrl%n_dx     = 3
    ! ctrl%dxs      = 0.0_wp 
    ! ! ctrl%dxs(1:3) = [0.05_wp,0.1_wp,0.2_wp]
    ! ctrl%dxs(1:3) = [0.005_wp,0.01_wp,0.02_wp]

    write(*,*) "dtts: ", ctrl%dtts(1:ctrl%n_dtt) 
    write(*,*) "dxs:  ", ctrl%dxs(1:ctrl%n_dx) 

    open(unit=15,file=trim(outfldr)//"slab_dt_factor.txt",status="UNKNOWN")
    write(15,"(a12,a12,a12)") "dx", "dt", "factor" 

    do q1 = 1, ctrl%n_dx

        ! Reset factor to a small value 
        factor = 1e-5 

    do q2 = 1, ctrl%n_dtt 

        if (factor .ge. 3.0) then 
            ! If factor has already surpassed one, meaning 
            ! model has become unstable for ctrl%dtts(q2-1),
            ! then don't run the model further, just set 
            ! factor to high value 

            factor = 100.0_wp 

        else 
            ! factor is still small, run the model for this timestep value

            ctrl%dx  = ctrl%dxs(q1) 
            ctrl%dtt = ctrl%dtts(q2) 

            call run_yelmo_test(factor,ctrl)
            
        end if 

        write(15,"(g12.3,g12.3,g12.3)") ctrl%dxs(q1), ctrl%dtts(q2), factor 
        write(*, "(a,g12.3,g12.3,g12.3)") "factor", ctrl%dxs(q1), ctrl%dtts(q2), factor
    end do 
    end do 

    close(15) 

    ! Stop timing 
    call yelmo_cpu_time(cpu_end_time,cpu_start_time,cpu_dtime)
    
    write(*,"(a,f12.3,a)") "Time  = ",cpu_dtime/60.0 ," min"
    
end if 



contains
    
    subroutine run_yelmo_test(factor,ctrl)

        implicit none 

        real(wp),       intent(OUT) :: factor 
        type(ctrl_par), intent(IN)  :: ctrl  

        ! Local variables 
        type(yelmo_class)     :: yelmo1

        integer  :: n
        real(wp) :: time 
        real(wp) :: xmin, xmax, ymin, ymax 
        real(wp), allocatable :: dh(:,:) 
        real(wp) :: stdev

        ! Define the time 
        time = ctrl%time_init 

        ! Define the domain and grid
        xmin =  0.0_prec 
        xmax =  (ctrl%nx-1)*ctrl%dx  
        ymax =  (ctrl%ny-1)*ctrl%dx/2.0_prec 
        ymin = -ymax
        call yelmo_init_grid(yelmo1%grd,ctrl%grid_name,units="km", &
                             x0=xmin,dx=ctrl%dx,nx=ctrl%nx, &
                             y0=ymin,dy=ctrl%dx,ny=ctrl%ny)

        ! Initialize data objects
        call yelmo_init(yelmo1,filename=ctrl%path_par,grid_def="none",time=time,load_topo=.FALSE., &
                            domain=ctrl%domain,grid_name=ctrl%grid_name)

        ! ===== Intialize topography and set parameters =========
        
        yelmo1%bnd%z_bed = 1000.0_wp - ctrl%alpha*(yelmo1%grd%x)

        ! Define initial ice thickness 
        allocate(dh(yelmo1%grd%nx,yelmo1%grd%ny))
        call gen_random_normal(dh,0.0_wp,ctrl%H_stdev) 
        yelmo1%tpo%now%H_ice = ctrl%H0 + dh 

        ! Define surface elevation 
        yelmo1%tpo%now%z_srf = yelmo1%bnd%z_bed + yelmo1%tpo%now%H_ice

        ! Define reference ice thickness (for prescribing boundary values, potentially)
        yelmo1%bnd%H_ice_ref = ctrl%H0 

        ! =======================================================

        ! Initialize the yelmo state (dyn,therm,mat)
        call yelmo_init_state(yelmo1,time=time,thrm_method="robin-cold")

        ! Advance timesteps
        do n = 1, ctrl%nt 

            ! Get current time 
            time = ctrl%time_init + n*ctrl%dtt

            ! == Yelmo ice sheet ====================
            call yelmo_update(yelmo1,time)

        end do 

        ! Calculate summary 
        call calc_stdev(stdev,yelmo1%tpo%now%H_ice)
        factor = stdev / max(ctrl%H_stdev,1e-5)

        call yelmo_end(yelmo1,time) 

        return 

    end subroutine run_yelmo_test

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
        ! call yelmo_write_step_model_metrics(filename,ylmo,n,ncid)

        ! == yelmo_topography ==
        call nc_write(filename,"H_ice",ylmo%tpo%now%H_ice,units="m",long_name="Ice thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf",ylmo%tpo%now%z_srf,units="m",long_name="Surface elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        ! call nc_write(filename,"mask_bed",ylmo%tpo%now%mask_bed,units="",long_name="Bed mask", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        ! call nc_write(filename,"N_eff",ylmo%dyn%now%N_eff,units="Pa",long_name="Effective pressure", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"dzsrfdt",ylmo%tpo%now%dzsrfdt,units="m/a",long_name="Surface elevation change", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dHicedt",ylmo%tpo%now%dHicedt,units="m/a",long_name="Ice thickness change", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! call nc_write(filename,"H_grnd",ylmo%tpo%now%H_grnd,units="m",long_name="Ice thickness overburden", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! call nc_write(filename,"bmb",ylmo%tpo%now%bmb,units="m/a",long_name="Basal mass balance", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! == yelmo_thermodynamics ==
        ! call nc_write(filename,"T_ice",ylmo%thrm%now%T_ice,units="K",long_name="Ice temperature", &
        !               dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        ! call nc_write(filename,"T_prime",ylmo%thrm%now%T_ice-ylmo%thrm%now%T_pmp,units="deg C",long_name="Homologous ice temperature", &
        !               dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        ! call nc_write(filename,"f_pmp",ylmo%thrm%now%f_pmp,units="1",long_name="Fraction of grid point at pmp", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! call nc_write(filename,"T_prime_b",ylmo%thrm%now%T_prime_b,units="deg C",long_name="Homologous basal ice temperature", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! call nc_write(filename,"H_w",ylmo%thrm%now%H_w,units="m",long_name="Basal water layer", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! == yelmo_material ==
!         call nc_write(filename,"visc_int",ylmo%mat%now%visc_int,units="Pa a m",long_name="Vertically integrated viscosity", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"visc",ylmo%mat%now%visc,units="Pa a",long_name="Viscosity", &
!                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        ! call nc_write(filename,"ATT_bar",ylmo%mat%now%ATT_bar,units="a^-1 Pa^-3",long_name="Vertically averaged rate factor", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"ATT",ylmo%mat%now%ATT,units="a^-1 Pa^-3",long_name="Rate factor", &
!                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        
        ! call nc_write(filename,"Q_ice_b",ylmo%thrm%now%Q_ice_b,units="J a-1 m-2",long_name="Basal ice heat flux", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        ! call nc_write(filename,"Q_strn",ylmo%thrm%now%Q_strn/(rho_ice*ylmo%thrm%now%cp),units="K a-1",long_name="Strain heating", &
        !               dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        ! call nc_write(filename,"Q_b",ylmo%thrm%now%Q_b,units="J a-1 m-2",long_name="Basal frictional heating", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! == yelmo_dynamics ==

        ! call nc_write(filename,"cf_ref",ylmo%dyn%now%cf_ref,units="--",long_name="Bed friction scalar", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        ! call nc_write(filename,"c_bed",ylmo%dyn%now%c_bed,units="Pa",long_name="Bed friction coefficient", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"beta",ylmo%dyn%now%beta,units="Pa a m-1",long_name="Basal friction coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
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
        
        ! call nc_write(filename,"de",ylmo%mat%now%strn%de,units="a^-1",long_name="Vertically averaged strain rate", &
        !               dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

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

    subroutine calc_stdev(stdev,var)

        implicit none 

        real(wp), intent(OUT) :: stdev 
        real(wp), intent(IN)  :: var(:,:) 

        ! Local variables 
        real(wp) :: mean 
        integer  :: n 

        n = size(var,1)*size(var,2)

        ! Calculate the mean 
        mean = sum(var) / real(n,wp)

        ! Calculate standard deviation 
        stdev = sqrt(sum( (var-mean)**2 ) / real(n-1,wp))

        return 

    end subroutine calc_stdev

end program yelmo_slab



