program yelmo_slab
    ! For testing an idealized slab domain

    use nml
    use ncio 
    use yelmo 
    use deformation 

    implicit none 

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

        integer                 :: n_dx 
        integer                 :: n_dt
        real(wp)                :: dx_range(2)
        real(wp)                :: dt_range(2) 
        real(wp)                :: conv_tol 

        character(len=512)  :: path_par 
        
        ! Internal variables
        real(wp), allocatable   :: dxs(:)
        real(wp), allocatable   :: dts(:)
        real(wp), allocatable   :: factors(:)
    end type 

    type(ctrl_par) :: ctrl 

    type results_class 
        real(wp), allocatable   :: dx
        real(wp), allocatable   :: dt
        real(wp), allocatable   :: factor
        real(wp), allocatable   :: H_mean
        real(wp), allocatable   :: ux_mean
        real(wp), allocatable   :: uxb_mean
        real(wp), allocatable   :: uxs_mean
    end type 

    type(results_class) :: res
        
    character(len=256) :: outfldr, file2D, file1D, file_restart
    character(len=512) :: path_par 
    
    integer  :: q, q1, q2, qmax, qmax1

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
    
    call nml_read(path_par,"ctrl","n_dx",         ctrl%n_dx)
    call nml_read(path_par,"ctrl","n_dt",         ctrl%n_dt)
    call nml_read(path_par,"ctrl","dx_range",     ctrl%dx_range)
    call nml_read(path_par,"ctrl","dt_range",     ctrl%dt_range)
    call nml_read(path_par,"ctrl","conv_tol",     ctrl%conv_tol)
    
    ! Set ctrl parameters for later use 
    ctrl%path_par   = path_par 
    
    ! Define default grid name for completeness 
    ctrl%grid_name = trim(ctrl%domain) 
    
if (ctrl%dtt .ne. 0.0) then 
    ! == Perform one simulation with an outer timestep of ctrl%dtt =======
    
    call run_yelmo_test(res,ctrl,file2D)
    
    ! Write summary 
    write(*,*) "====== "//trim(ctrl%domain)//" ======="
    write(*, "(a,8a12)")   "factor","dx", "dt", "factor", "H", "ux_bar", "ux_b", "ux_s"
    write(*, "(a,8g12.3)") "factor", res%dx, res%dt, res%factor, res%H_mean, res%ux_mean, res%uxb_mean, res%uxs_mean
    
else 
    ! === Perform ensemble of simulations with multiple values of dx and dt =====
    ! Use bisection method to refine estimate of maximum stable timestep

    ! Define dx values to be tested
    ! (evenly spaced in log-space for dx_range)

    allocate(ctrl%dxs(ctrl%n_dx)) 
    call sample_log10_space(ctrl%dxs,ctrl%dx_range(1),ctrl%dx_range(2))

    write(*,*) "dxs:  ", ctrl%dxs

    ! Define total number of dt values, but intially set to zero
    allocate(ctrl%dts(ctrl%n_dt))
    ctrl%dts = 0.0_wp 

    ! Define factors array too
    allocate(ctrl%factors(ctrl%n_dt))
    ctrl%factors = 0.0_wp 

    
    open(unit=15,file=trim(outfldr)//"slab_dt_factor.txt",status="UNKNOWN")
    write(15,"(7a12)") "dx", "dt", "factor", "H", "ux_bar", "ux_b", "ux_s"

    do q1 = 1, ctrl%n_dx

        ! Reset all dt and factor values to zero
        ctrl%dts     = 0.0_wp 
        ctrl%factors = 0.0_wp 

        ! Initially test 5 values of dt over the desired range
        call sample_log10_space(ctrl%dts(1:5),ctrl%dt_range(1),ctrl%dt_range(2))

        do q2 = 1, 5  

            ! Assign current dx and dt values and test model
            ctrl%dx  = ctrl%dxs(q1)
            ctrl%dtt = ctrl%dts(q2) 
            call run_yelmo_test(res,ctrl)
            
            ! Store results 
            ctrl%factors(q2) = res%factor 

            ! Write current results to output table
            write(15,"(8g12.3)") res%dx, res%dt, res%factor, res%H_mean, res%ux_mean, res%uxb_mean, res%uxs_mean
            write(*, "(8g12.3)") "factor", res%dx, res%dt, res%factor, res%H_mean, res%ux_mean, res%uxb_mean, res%uxs_mean
            
        end do 

        if (maxval(ctrl%factors,mask=ctrl%dts.ne.0.0_wp) .lt. 1.0_wp) then 
            ! All dt values tested were stable, do not test more 
            ! (implies maximum stable timestep is above the tested range)

            ! Do nothing 

        else if (minval(ctrl%factors,mask=ctrl%dts.ne.0.0_wp) .gt. 1.0_wp) then
            ! No dt values tested were stable, do not test more
            ! (implies maximum stable timestep is below the tested range)

            ! Do nothing 

        else 
            ! Continue testing dt values 

            do q2 = 6, ctrl%n_dt 

                ! Determine index of current estimate of maximum stable timestep
                qmax  = maxloc(ctrl%dts,mask=ctrl%factors.lt.1.0_wp,dim=1)
                qmax1 = minloc(abs(ctrl%dts-ctrl%dts(qmax)),mask=ctrl%dts.gt.ctrl%dts(qmax),dim=1)

                ! Add another timestep between this one and the next higher 
                ctrl%dts(q2) = 0.5_wp*(ctrl%dts(qmax)+ctrl%dts(qmax1))

                ! Assign current dx and dt values and test model
                ctrl%dx  = ctrl%dxs(q1)
                ctrl%dtt = ctrl%dts(q2) 
                call run_yelmo_test(res,ctrl)
                
                ! Store results 
                ctrl%factors(q2) = res%factor 

                ! Write current results to output table
                write(15,"(8g12.3)") res%dx, res%dt, res%factor, res%H_mean, res%ux_mean, res%uxb_mean, res%uxs_mean
                write(*, "(8g12.3)") "factor", res%dx, res%dt, res%factor, res%H_mean, res%ux_mean, res%uxb_mean, res%uxs_mean
                
                ! Exit loop if dt has converged enough
                ! (checking fractional difference between current and previous timestep)
                if (abs(ctrl%dts(q2)-ctrl%dts(q2-1))/ctrl%dts(q2-1) .lt. ctrl%conv_tol) then 
                    exit 
                end if 

            end do 

        end if 

    end do 

    close(15) 

end if 

! Stop timing 
call yelmo_cpu_time(cpu_end_time,cpu_start_time,cpu_dtime)

write(*,"(a,f12.3,a)") "Time  = ",cpu_dtime/60.0 ," min"

! ===== DONE ======

contains
    
    subroutine run_yelmo_test(res,ctrl,file2D)

        implicit none 

        type(results_class),    intent(OUT) :: res 
        type(ctrl_par),         intent(IN)  :: ctrl  
        character(len=*),       intent(IN), optional :: file2D 

        ! Local variables 
        type(yelmo_class) :: yelmo1

        integer  :: n, j 
        real(wp) :: time 
        real(wp) :: xmin, xmax, ymin, ymax 
        real(wp), allocatable :: dh(:,:) 
        real(wp) :: u0, ub0 
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
        
        yelmo1%bnd%z_bed = 10000.0_wp - ctrl%alpha*(yelmo1%grd%x)

        ! Define initial ice thickness 
        allocate(dh(yelmo1%grd%nx,yelmo1%grd%ny))
        call gen_random_normal(dh,0.0_wp,ctrl%H_stdev) 
        yelmo1%tpo%now%H_ice = ctrl%H0 !+ dh 

        ! Define surface elevation 
        yelmo1%tpo%now%z_srf = yelmo1%bnd%z_bed + yelmo1%tpo%now%H_ice

        ! Define reference ice thickness (for prescribing boundary values, potentially)
        yelmo1%bnd%H_ice_ref = ctrl%H0 

        ! Initialize velocity to u0 value too
        call calc_u0(u0,ctrl%H0,yelmo1%dyn%par%beta_const, &
                        yelmo1%dyn%par%visc_const,ctrl%alpha, &
                        yelmo1%bnd%c%rho_ice,yelmo1%bnd%c%g)
        
        yelmo1%dyn%now%ux_bar = u0

        if (trim(yelmo1%dyn%par%solver) .ne. "diva") then 
            call calc_ub0(ub0,ctrl%H0,yelmo1%dyn%par%beta_const, &
                        yelmo1%dyn%par%visc_const,ctrl%alpha, &
                        yelmo1%bnd%c%rho_ice,yelmo1%bnd%c%g)
            
            yelmo1%dyn%now%ux_b = ub0 
        else 
            ub0 = 0.0_wp 
        end if 

        write(*,*) "u0, ub0 = ", u0, ub0

        ! =======================================================

        ! Initialize the yelmo state (dyn,therm,mat)
        call yelmo_init_state(yelmo1,time=time,thrm_method="robin-cold")

        yelmo1%tpo%par%topo_fixed = .TRUE. 

        write(*,*) "=== Timesteps with fixed topo. ==="

        ! Advance timesteps
        do n = 1, 200 

            ! Get current time 
            time = ctrl%time_init + n*ctrl%dtt

            ! == Yelmo ice sheet ===================================================
            call yelmo_update(yelmo1,time)

            if (mod(time,10.0)==0 .and. (.not. yelmo_log)) then
                write(*,"(a,f14.4)") "yelmo:: time = ", time
            end if  

            if (yelmo1%dyn%par%ssa_iter_now .lt. 10) then 
                ! Solution is near stable
                write(*,*) "n = ", n 
                exit 
            end if 

        end do 

        call yelmo_set_time(yelmo1,ctrl%time_init)
        yelmo1%tpo%par%topo_fixed   = .FALSE. 
        yelmo1%tpo%now%H_ice        = ctrl%H0 + dh 
        
        ! Make sure all values are the same in y-direction
        do j = 1, yelmo1%grd%ny 
            yelmo1%tpo%now%H_ice(:,j) = yelmo1%tpo%now%H_ice(:,1)
        end do 

        write(*,*) "=== Timesteps with prognostic topo. ==="

        if (present(file2D)) then 
            ! Initialize output file 
            call yelmo_write_init(yelmo1,file2D,time_init=time,units="years")
        
            ! Write initial state 
            call write_step_2D(yelmo1,file2D,time=time) 
        end if 

        ! Advance timesteps
        do n = 1, ctrl%nt 

            ! Get current time 
            time = ctrl%time_init + n*ctrl%dtt

            ! == Yelmo ice sheet ====================
            call yelmo_update(yelmo1,time)

        end do 

        if (present(file2D)) then 
            ! Write final state 
            call write_step_2D(yelmo1,file2D,time=time) 
        end if 

        ! Calculate summary 
        call calc_stdev(stdev,yelmo1%tpo%now%H_ice,ctrl%H0)
        res%factor = stdev / max(ctrl%H_stdev,1e-5)

        ! Limit factor to a reasonable value 
        res%factor = min(res%factor,1e3)

        ! Get mean values along center x-profile

        j = floor(yelmo1%grd%ny/2.0_wp)
        res%H_mean   = sum(yelmo1%tpo%now%H_ice(:,j))  / real(yelmo1%grd%nx,wp)
        res%ux_mean  = sum(yelmo1%dyn%now%ux_bar(:,j)) / real(yelmo1%grd%nx,wp)
        res%uxb_mean = sum(yelmo1%dyn%now%ux_b(:,j))   / real(yelmo1%grd%nx,wp)
        res%uxs_mean = sum(yelmo1%dyn%now%ux_s(:,j))   / real(yelmo1%grd%nx,wp)
        
        ! Also save dx and dt 
        res%dx = ctrl%dx*1e3    ! Save dx in [m]
        res%dt = ctrl%dtt 

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
        
        call nc_write(filename,"dzsdt",ylmo%tpo%now%dzsdt,units="m/a",long_name="Surface elevation change", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dHidt",ylmo%tpo%now%dHidt,units="m/a",long_name="Ice thickness change", &
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

        ! call nc_write(filename,"cb_ref",ylmo%dyn%now%cb_ref,units="--",long_name="Bed friction scalar", &
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
        
        call nc_write(filename,"ux_b",ylmo%dyn%now%ux_b,units="m/a",long_name="Basal velocity (x)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"ux_i_bar",ylmo%dyn%now%ux_i_bar,units="m/a",long_name="Vertically averaged shearing velocity (x)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"ux_i",ylmo%dyn%now%ux_i,units="m/a",long_name="Horizontal shearing velocity (x)", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        
        call nc_write(filename,"ux",ylmo%dyn%now%ux,units="m/a",long_name="Horizontal velocity (x)", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        
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

    subroutine calc_stdev(stdev,var,mean_ref)

        implicit none 

        real(wp), intent(OUT) :: stdev 
        real(wp), intent(IN)  :: var(:,:) 
        real(wp), intent(IN), optional  :: mean_ref
        ! Local variables 
        real(wp) :: mean 
        integer  :: n 

        n = size(var,1)*size(var,2)

        if (present(mean_ref)) then
            ! Use reference mean value instead of calculating it
            ! (to ensure high stdev value when mean is far from mean_ref) 
            mean = mean_ref 
        else 
            ! Calculate the mean 
            mean = sum(var) / real(n,wp)
        end if 

        ! Calculate standard deviation 
        stdev = sqrt(sum( (var-mean)**2 ) / real(n-1,wp))

        return 

    end subroutine calc_stdev

    elemental subroutine calc_u0(u0,H0,beta,eta,alpha,rho_ice,g)
        ! Prescribe initial velocity as u0 

        real(wp), intent(OUT) :: u0 
        real(wp), intent(IN)  :: H0
        real(wp), intent(IN)  :: beta 
        real(wp), intent(IN)  :: eta            ! Viscosity!
        real(wp), intent(IN)  :: alpha
        real(wp), intent(IN)  :: rho_ice 
        real(wp), intent(IN)  :: g 
        
        ! Local variables 
        real(wp) :: F2 
 
        F2 = H0/(3.0*eta)
        u0 = (rho_ice*g) * H0 * alpha * (1.0+beta*F2)/beta

        return 

    end subroutine calc_u0

    elemental subroutine calc_ub0(ub0,H0,beta,eta,alpha,rho_ice,g)
        ! Prescribe initial velocity as u0 

        real(wp), intent(OUT) :: ub0 
        real(wp), intent(IN)  :: H0
        real(wp), intent(IN)  :: beta 
        real(wp), intent(IN)  :: eta            ! Viscosity!
        real(wp), intent(IN)  :: alpha
        real(wp), intent(IN)  :: rho_ice 
        real(wp), intent(IN)  :: g 

        ub0 = (rho_ice*g) * H0 * alpha / beta

        return 

    end subroutine calc_ub0

    subroutine sample_log10_space(vals,val_min,val_max)

        implicit none 

        real(wp), intent(OUT) :: vals(:)
        real(wp), intent(IN)  :: val_min
        real(wp), intent(IN)  :: val_max
        
        ! Local variables
        integer :: n, q 

        n = size(vals)

        vals  = 0.0_wp
        
        do q = 1, n
            vals(q) = 10.0_dp**(log10(real(val_min,dp))+(log10(real(val_max,dp))-log10(real(val_min,dp)))*(q-1)/real(n-1,dp))
        end do

        return

    end subroutine sample_log10_space

end program yelmo_slab



