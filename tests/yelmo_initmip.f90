

program yelmo_test

    use omp_lib

    use nml 
    use ncio 
    use yelmo 
    
    use basal_dragging 
    use yelmo_tools, only : gauss_values
    
    use ice_optimization 
    
    implicit none 

    type(yelmo_class)       :: yelmo1 
    
    character(len=256) :: outfldr, file1D, file2D, file_restart, domain 
    character(len=512) :: path_par  
    real(wp) :: time   
    integer    :: n
    
    logical,  allocatable :: mask_noice(:,:)  

    logical, parameter  :: test_restart = .FALSE. 
    real(wp), parameter :: time_r       = 50.0_wp 
    type(yelmo_class)   :: yelmo_r
    character(len=256)  :: file1D_r, file2D_r, file_restart_r

    type ctrl_params
        character(len=56) :: run_step
        real(wp) :: time_init
        real(wp) :: time_end
        real(wp) :: time_equil      ! Only for spinup
        real(wp) :: dtt
        real(wp) :: dt1D_out
        real(wp) :: dt2D_out

        logical  :: transient_clim
        logical  :: with_ice_sheet 
        logical  :: optimize
        
        character(len=512) :: clim_nm
          
        logical :: load_cb_ref 
        character(len=256) :: file_cb_ref 

        logical :: load_bmelt
        character(len=256) :: file_bmelt 

        real(wp) :: bmb_shlf_const
        real(wp) :: dT_ann
        real(wp) :: z_sl 
    end type 

    type opt_params
        real(wp) :: cf_init
        real(wp) :: cf_min_par
        real(wp) :: cf_max_par 
        real(wp) :: tau_c 
        real(wp) :: H0
        real(wp) :: sigma_err 
        real(wp) :: sigma_vel 
        character(len=56) :: fill_method
        
        real(wp) :: rel_tau 
        real(wp) :: rel_tau1 
        real(wp) :: rel_tau2
        real(wp) :: rel_time1
        real(wp) :: rel_time2
        real(wp) :: rel_m 

        real(wp), allocatable :: cf_min(:,:) 
        real(wp), allocatable :: cf_max(:,:) 
        
    end type 

    type(ctrl_params)   :: ctl
    type(opt_params)    :: opt  

    real(8) :: cpu_start_time, cpu_end_time, cpu_dtime  
    integer  :: perr 

    ! Start timing 
    call yelmo_cpu_time(cpu_start_time)
    
    ! Determine the parameter file from the command line 
    call yelmo_load_command_line_args(path_par)

    ! Timing and other parameters 
    call nml_read(path_par,"ctrl","time_init",      ctl%time_init)          ! [yr] Starting time
    call nml_read(path_par,"ctrl","time_end",       ctl%time_end)           ! [yr] Ending time
    call nml_read(path_par,"ctrl","time_equil",     ctl%time_equil)         ! [yr] Years to equilibrate first
    call nml_read(path_par,"ctrl","dtt",            ctl%dtt)                ! [yr] Main loop time step 
    call nml_read(path_par,"ctrl","dt1D_out",       ctl%dt1D_out)           ! [yr] Frequency of 1D output 
    call nml_read(path_par,"ctrl","dt2D_out",       ctl%dt2D_out)           ! [yr] Frequency of 2D output 
    call nml_read(path_par,"ctrl","transient_clim", ctl%transient_clim)     ! Calculate transient climate? 
    call nml_read(path_par,"ctrl","with_ice_sheet", ctl%with_ice_sheet)     ! Include an active ice sheet 
    call nml_read(path_par,"ctrl","optimize",       ctl%optimize)           ! Optimize basal friction cb_ref field
    call nml_read(path_par,"ctrl","clim_nm",        ctl%clim_nm)            ! Namelist group holding climate information
    
    call nml_read(path_par,"ctrl","load_cb_ref",    ctl%load_cb_ref)        ! Load cb_ref from file? Otherwise define from till_cf_ref + inline tuning
    call nml_read(path_par,"ctrl","file_cb_ref",    ctl%file_cb_ref)        ! Filename holding cb_ref to load 

    call nml_read(path_par,"ctrl","load_bmelt",     ctl%load_bmelt)         ! Load bmelt from file?
    call nml_read(path_par,"ctrl","file_bmelt",     ctl%file_bmelt)         ! Filename holding bmelt field to load 
        
    ! Load climate (eg, clim_pd or clim_lgm)
    call nml_read(path_par,ctl%clim_nm,  "bmb_shlf_const",  ctl%bmb_shlf_const)            ! [yr] Constant imposed bmb_shlf value
    call nml_read(path_par,ctl%clim_nm,  "dT_ann",          ctl%dT_ann)                    ! [K] Temperature anomaly (atm)
    call nml_read(path_par,ctl%clim_nm,  "z_sl",            ctl%z_sl)                      ! [m] Sea level relative to present-day

    if (ctl%optimize) then 
        ! Load optimization parameters 

        call nml_read(path_par,"opt","cf_init",     opt%cf_init)
        call nml_read(path_par,"opt","cf_min",      opt%cf_min_par)
        call nml_read(path_par,"opt","cf_max",      opt%cf_max_par)
        call nml_read(path_par,"opt","tau_c",       opt%tau_c)
        call nml_read(path_par,"opt","H0",          opt%H0)
        call nml_read(path_par,"opt","sigma_err",   opt%sigma_err)   
        call nml_read(path_par,"opt","sigma_vel",   opt%sigma_vel)   
        call nml_read(path_par,"opt","fill_method", opt%fill_method)   
        
        call nml_read(path_par,"opt","rel_tau1",    opt%rel_tau1)   
        call nml_read(path_par,"opt","rel_tau2",    opt%rel_tau2)  
        call nml_read(path_par,"opt","rel_time1",   opt%rel_time1)    
        call nml_read(path_par,"opt","rel_time2",   opt%rel_time2) 
        call nml_read(path_par,"opt","rel_m",       opt%rel_m)
    end if 

    ! Assume program is running from the output folder
    outfldr = "./"

    ! Define input and output locations 
    file1D       = trim(outfldr)//"yelmo1D.nc"
    file2D       = trim(outfldr)//"yelmo2D.nc"
    file_restart = trim(outfldr)//"yelmo_restart.nc"

    ! === Initialize ice sheet model =====

    ! Initialize data objects and load initial topography
    call yelmo_init(yelmo1,filename=path_par,grid_def="file",time=ctl%time_init)

    ! Ensure optimization fields are allocated
    allocate(opt%cf_min(yelmo1%grd%nx,yelmo1%grd%ny))
    allocate(opt%cf_max(yelmo1%grd%nx,yelmo1%grd%ny))
    opt%cf_min = opt%cf_min_par 
    opt%cf_max = opt%cf_max_par 

    ! === Set initial boundary conditions for current time and yelmo state =====
    ! ybound: z_bed, z_sl, H_sed, smb, T_srf, bmb_shlf , Q_geo

    yelmo1%bnd%z_sl     = ctl%z_sl          ! [m]
    yelmo1%bnd%H_sed    = 0.0               ! [m]
    yelmo1%bnd%Q_geo    = 50.0              ! [mW/m2]
    
    ! Impose present-day surface mass balance and present-day temperature field plus any anomaly
    yelmo1%bnd%smb      = yelmo1%dta%pd%smb                 ! [m.i.e./a]
    yelmo1%bnd%T_srf    = yelmo1%dta%pd%T_srf + ctl%dT_ann  ! [K]
    
    if (ctl%load_bmelt) then

        ! Parse filenames with grid information
        call yelmo_parse_path(ctl%file_bmelt,yelmo1%par%domain,yelmo1%par%grid_name)
 
        call nc_read(ctl%file_bmelt,"bm_ac_reese",yelmo1%bnd%bmb_shlf)
        yelmo1%bnd%bmb_shlf = -yelmo1%bnd%bmb_shlf      ! Negative because bmb = -bmelt 
    else 
        yelmo1%bnd%bmb_shlf = ctl%bmb_shlf_const        ! [m.i.e./a]
    end if 

    yelmo1%bnd%T_shlf   = yelmo1%bnd%c%T0                ! [K]   

    if (ctl%dT_ann .lt. 0.0) yelmo1%bnd%T_shlf   = yelmo1%bnd%c%T0 + ctl%dT_ann*0.25_wp      ! [K] Oceanic temp anomaly
    
    call yelmo_print_bound(yelmo1%bnd)

    time = ctl%time_init 

    ! Define no-ice mask from present-day data
    allocate(mask_noice(yelmo1%grd%nx,yelmo1%grd%ny))
    mask_noice = .FALSE. 
    ! Present-day
    if (ctl%dT_ann .ge. 0.0) then
        where(yelmo1%dta%pd%H_ice .le. 0.0) mask_noice = .TRUE. 
    end if 

    ! Special treatment for Greenland
    if (trim(yelmo1%par%domain) .eq. "Greenland") then 
        
        ! Present-day
        if (ctl%dT_ann .ge. 0.0) then 
            ! Impose additional negative mass balance to no ice points of 2 [m.i.e./a] melting
            where(mask_noice) yelmo1%bnd%smb = yelmo1%dta%pd%smb - 2.0 
        end if 

    end if 

    ! Special treatment for Antarctica
    if (trim(yelmo1%par%domain) .eq. "Antarctica") then 
        
        ! Present-day
        if (ctl%dT_ann .ge. 0.0) then 
            where(mask_noice) yelmo1%bnd%bmb_shlf = -2.0                ! [m/a]        
        end if 

        ! LGM
        if (ctl%dT_ann .lt. 0.0) then 
            where(yelmo1%bnd%smb .le. 0.1) yelmo1%bnd%smb = 0.1         ! [m/a]
        end if 

        ! Present-day and LGM 
        where(yelmo1%bnd%regions .eq. 2.0) yelmo1%bnd%smb = -1.0        ! [m/a]

    end if 

    ! Initialize state variables (dyn,therm,mat)
    ! (initialize temps with robin method with a cold base)
    call yelmo_init_state(yelmo1,time=ctl%time_init,thrm_method="robin-cold")

    ! ===== basal friction optimization ======
    if (ctl%optimize) then 
        
        ! Ensure that cb_ref will be optimized (till_method == set externally) 
        yelmo1%dyn%par%till_method = -1  

        ! If not using restart, prescribe cb_ref to initial guess 
        if (.not. yelmo1%par%use_restart) then
            yelmo1%dyn%now%cb_ref = opt%cf_init 
        end if 

    else 

        ! ============================================================
        ! Load or define cb_ref 

        ! Allocate cb_ref and set it to till_cf_ref by default 
        yelmo1%dyn%now%cb_ref = yelmo1%dyn%par%till_cf_ref  

        if (ctl%load_cb_ref) then 

            ! Parse filename with grid information
            call yelmo_parse_path(ctl%file_cb_ref,yelmo1%par%domain,yelmo1%par%grid_name)

            ! Load cb_ref from specified file 
            call nc_read(ctl%file_cb_ref,"cb_ref",yelmo1%dyn%now%cb_ref)

            ! Make sure that present-day shelves have minimal cb_ref values 
            where(yelmo1%tpo%now%f_grnd .eq. 0.0) yelmo1%dyn%now%cb_ref = 0.05 

        else
            ! Define cb_ref inline 

            yelmo1%dyn%now%cb_ref = yelmo1%dyn%par%till_cf_ref 

        end if 

    end if 
    ! ============================================================
    
    if (.not. yelmo1%par%use_restart) then

        ! Run full model (tpo,dyn,thrm) without advection to clean initial topo via mass balance
        call yelmo_update_equil(yelmo1,time,time_tot=10.0_prec,dt=1.0_prec,topo_fixed=.FALSE.,tpo_solver="none")

        ! Run full model with correct solver (tpo,dyn,thrm)
        call yelmo_update_equil(yelmo1,time,time_tot=1.0_prec,dt=0.2_prec,topo_fixed=.FALSE.)

        ! Next equilibrate thermodynamics further and maintain constant ice topopgraphy (for speed)
        !call yelmo_update_equil(yelmo1,time,time_tot=1e3,dt=10.0_wp,topo_fixed=.TRUE.)

    end if 


    if (yelmo1%par%use_restart .and. yelmo1%par%restart_interpolated .eq. 1) then 

        write(*,*) "Running restart topo smoothing step..."

        ! Run model with no advection
        call yelmo_update_equil(yelmo1,time,time_tot=10.0_prec,dt=1.0_wp, &
                              tpo_solver="none",topo_fixed=.FALSE.,dyn_solver="ssa")

        ! Run thermodynamics with SSA solver very briefly to smooth it out
        call yelmo_update_equil(yelmo1,time,time_tot=10.0_prec,dt=1.0_wp, &
                                                topo_fixed=.TRUE.,dyn_solver="ssa")

        ! Run full model (tpo,dyn,thrm) with SSA solver very briefly to smooth it out
        call yelmo_update_equil(yelmo1,time,time_tot=1.0_prec,dt=0.2_wp, &
                                                topo_fixed=.FALSE.,dyn_solver="ssa")

        if (yelmo1%grd%dx .le. 8e3_wp) then 
            ! Perform additional ssa smoothing step for higher resolution simulations

            call yelmo_update_equil(yelmo1,time,time_tot=100.0_prec,dt=1.0_wp, &
                                                topo_fixed=.FALSE.,dyn_solver="ssa")

        end if 

    
    end if 

    ! Initialize output files 
    call yelmo_write_init(yelmo1,file2D,time_init=time,units="years")  
    call yelmo_write_reg_init(yelmo1,file1D,time_init=time,units="years",mask=yelmo1%bnd%ice_allowed)
        
    ! if (with_anom) then 
    !     ! Warm up the ice sheet to impose some changes 
    !     yelmo1%bnd%T_srf    = yelmo1%dta%pd%T_srf + 5.0
    !     yelmo1%bnd%bmb_shlf = -10.0 
    !     where (yelmo1%dta%pd%smb .lt. 0.0) yelmo1%bnd%smb = yelmo1%dta%pd%smb - 1.0 
    ! end if 

    ! Initialize continuous restart file 
    !call yelmo_restart_write(yelmo1,"yelmo_heavy.nc",time,init=.TRUE.)

    ! Advance timesteps
    do n = 0, ceiling((ctl%time_end-ctl%time_init)/ctl%dtt)

        ! Get current time 
        time = ctl%time_init + n*ctl%dtt

        if (test_restart) then 

            if (time .eq. time_r+ctl%dtt) then 

                file1D_r       = trim(outfldr)//"yelmo1D_r.nc"
                file2D_r       = trim(outfldr)//"yelmo2D_r.nc"
                file_restart_r = trim(outfldr)//"yelmo_restart_r.nc"

                ctl%dt2D_out = ctl%dtt 
                ctl%dt1D_out = ctl%dtt 

                ! Initialize data objects and load initial topography
                call yelmo_init(yelmo_r,filename=path_par,grid_def="file",time=time_r)

                yelmo_r%par%restart     = "yelmo_restart.nc"
                yelmo_r%par%use_restart = .TRUE. 
                yelmo_r%bnd = yelmo1%bnd

                ! Initialize state variables (dyn,therm,mat)
                ! (initialize temps with robin method with a cold base)
                call yelmo_init_topo(yelmo_r,path_par,time_r)
                call yelmo_init_state(yelmo_r,time=time_r,thrm_method="robin-cold")

                ! Write restart file to compare with expected 
                call yelmo_restart_write(yelmo_r,file_restart_r,time_r)

                ! 2D file 
                call yelmo_write_init(yelmo_r,file2D_r,time_init=time_r,units="years")  
                
                ! 1D file 
                call yelmo_write_reg_init(yelmo_r,file1D_r,time_init=time_r,units="years",mask=yelmo_r%bnd%ice_allowed)
                
                call write_step_2D(yelmo_r,file2D_r,time=time_r)
                call yelmo_write_reg_step(yelmo_r,file1D_r,time=time_r)  
            end if 

            if (time .gt. time_r) then
                ! Update restarted model
                call yelmo_update(yelmo_r,time)
            end if 

        end if 

!         ! Update temperature and smb as needed in time (ISMIP6)
!         if (time .ge. -10e6 .and. time .lt. -10e3) then 
!             ! Glacial period, impose cold climate 
!             yelmo1%bnd%T_srf = yelmo1%dta%pd%T_srf - 10.0 

!         else if (time .ge. -10e3 .and. time .lt. -6e3) then
!             ! Holocene optimum 
!             yelmo1%bnd%T_srf = yelmo1%dta%pd%T_srf + 1.0 

!         else  ! time .ge. -6e3
!             ! Entering Holocene, impose present-day temperatures 
!             yelmo1%bnd%T_srf = yelmo1%dta%pd%T_srf
!         end if 
        
        ! Spin-up procedure...
        if (ctl%with_ice_sheet .and. (time-ctl%time_init) .lt. ctl%time_equil) then 
            
            if (ctl%optimize) then 
                ! ===== basal friction optimization ==================

                ! === Optimization parameters =========
                
                ! Update model relaxation time scale and error scaling (in [m])
                call optimize_set_transient_param(opt%rel_tau,time,time1=opt%rel_time1, &
                                        time2=opt%rel_time2,p1=opt%rel_tau1,p2=opt%rel_tau2,m=opt%rel_m)
                
                ! Set model tau, and set yelmo relaxation switch (2: gl-line and shelves relaxing; 0: no relaxation)
                yelmo1%tpo%par%topo_rel_tau = opt%rel_tau 
                yelmo1%tpo%par%topo_rel     = 2
                if (time .gt. opt%rel_time2) yelmo1%tpo%par%topo_rel = 0 
                
                ! === Optimization update step =========

                ! Update cb_ref based on error metric(s) 
                call optimize_cb_ref(yelmo1%dyn%now%cb_ref,yelmo1%tpo%now%H_ice, &
                                    yelmo1%tpo%now%dHidt,yelmo1%bnd%z_bed,yelmo1%bnd%z_sl,yelmo1%dyn%now%ux_s,yelmo1%dyn%now%uy_s, &
                                    yelmo1%dta%pd%H_ice,yelmo1%dta%pd%uxy_s,yelmo1%dta%pd%H_grnd, &
                                    opt%cf_min,opt%cf_max,yelmo1%tpo%par%dx,opt%sigma_err,opt%sigma_vel,opt%tau_c,opt%H0, &
                                    dt=ctl%dtt,fill_method=opt%fill_method,fill_dist=80.0_wp)
            
            else 
                ! ===== relaxation spinup ==================

                ! Turn on relaxation for now, to let thermodynamics equilibrate
                ! without changing the topography too much. Important when 
                ! effective pressure = f(thermodynamics).

                yelmo1%tpo%par%topo_rel     = 2
                yelmo1%tpo%par%topo_rel_tau = 50.0 
                write(*,*) "timelog, tau = ", yelmo1%tpo%par%topo_rel_tau

            end if 

        else if ( (time-ctl%time_init) .eq. ctl%time_equil) then
            ! Finally, ensure all relaxation is disabled and continue as normal.

                yelmo1%tpo%par%topo_rel     = 0
                write(*,*) "timelog, relaxation off..."
              
        end if 
        ! ====================================================

        ! Update ice sheet 
        if (ctl%with_ice_sheet) call yelmo_update(yelmo1,time)

        ! if (time .lt. 130.0) then 
        !     call yelmo_update(yelmo1,time)
        ! else 
        !     call yelmo_update(yelmo1,time,"yelmo_heavy.nc")
        ! end if 

        ! == MODEL OUTPUT =======================================================

        if (mod(nint(time*100),nint(ctl%dt2D_out*100))==0) then
            call write_step_2D(yelmo1,file2D,time=time)
        end if 

        if (mod(nint(time*100),nint(ctl%dt1D_out*100))==0) then 
            call yelmo_write_reg_step(yelmo1,file1D,time=time) 
        end if 

        if (test_restart) then 

            if (time .eq. time_r) then
                ! Write restart file to load from
                call yelmo_restart_write(yelmo1,file_restart,time)
            end if 

            if (time .gt. time_r) then
                ! Write restarted output files 
                if (mod(nint(time*100),nint(ctl%dt2D_out*100))==0) then
                    call write_step_2D(yelmo_r,file2D_r,time=time)
                end if 

                if (mod(nint(time*100),nint(ctl%dt1D_out*100))==0) then 
                    call yelmo_write_reg_step(yelmo_r,file1D_r,time=time)  
                end if
            end if 

        end if 

        if (mod(time,10.0)==0 .and. (.not. yelmo_log)) then
            write(*,"(a,f14.4)") "yelmo:: time = ", time
        end if 
        
    end do 
    ! == Finished time loop == 

    ! Write a final restart file 
    if (.not. test_restart) then 
        call yelmo_restart_write(yelmo1,file_restart,time)
    end if 

    ! Finalize program
    call yelmo_end(yelmo1,time=time)
    
    ! Stop timing 
    call yelmo_cpu_time(cpu_end_time,cpu_start_time,cpu_dtime)

    write(*,"(a,f12.3,a)") "Time  = ",cpu_dtime/60.0 ," min"
    write(*,"(a,f12.1,a)") "Speed = ",(1e-3*(ctl%time_end-ctl%time_init))/(cpu_dtime/3600.0), " kiloyears / hr"

contains

    subroutine write_step_2D(ylmo,filename,time)

        implicit none 
        
        type(yelmo_class), intent(IN) :: ylmo
        character(len=*),  intent(IN) :: filename
        real(prec), intent(IN) :: time

        ! Local variables
        integer    :: ncid, n
        real(prec) :: time_prev 

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

        ! Write present-day data metrics (rmse[H],etc)
        call yelmo_write_step_pd_metrics(filename,ylmo,n,ncid)
        
        ! == ISMIP6 specific variables 
        ! http://www.climate-cryosphere.org/wiki/index.php?title=InitMIP-Greenland#Appendix_2_.E2.80.93_Naming_conventions.2C_upload_and_model_output_data.
        
        call nc_write(filename,"H_ice",ylmo%tpo%now%H_ice,units="m",long_name="Ice thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf",ylmo%tpo%now%z_srf,units="m",long_name="Surface elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_bed",ylmo%bnd%z_bed,units="m",long_name="Bedrock elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid) 
        call nc_write(filename,"Q_geo",ylmo%bnd%Q_geo,units="mW/m^2",long_name="Geothermal heat flux", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"smb",ylmo%bnd%smb,units="m/a ice equiv.",long_name="Surface mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"bmb",ylmo%tpo%now%bmb,units="m/a ice equiv.",long_name="Basal mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"fmb",ylmo%tpo%now%fmb,units="m/a ice equiv.",long_name="Margin-front mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dmb",ylmo%tpo%now%dmb,units="m/a ice equiv.",long_name="Subgrid discharge rate", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dHidt",ylmo%tpo%now%dHidt,units="m/a",long_name="Ice thickness change", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"ux_s",ylmo%dyn%now%ux_s,units="m/a",long_name="Surface velocity (x)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy_s",ylmo%dyn%now%uy_s,units="m/a",long_name="Surface velocity (y)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"uz_s",ylmo%dyn%now%uz(:,:,ylmo%par%nz_ac),units="m/a",long_name="Surface velocity (z)", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"ux_b",ylmo%dyn%now%ux_b,units="m/a",long_name="Basal velocity (x)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy_b",ylmo%dyn%now%uy_b,units="m/a",long_name="Basal velocity (y)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"uz_b",ylmo%dyn%now%uz(:,:,1),units="m/a",long_name="Basal velocity (z)", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"ux_bar",ylmo%dyn%now%ux_bar,units="m/a",long_name="Vertically averaged velocity (x)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy_bar",ylmo%dyn%now%uy_bar,units="m/a",long_name="Vertically averaged velocity (y)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"T_ice_s",ylmo%thrm%now%T_ice(:,:,ylmo%par%nz_aa),units="K",long_name="Surface ice temperature", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"T_ice_b",ylmo%thrm%now%T_ice(:,:,1),units="K",long_name="Basal homologous ice temperature", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taub",ylmo%dyn%now%taub,units="Pa",long_name="Basal dragging stress (magnitude)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"cmb",ylmo%tpo%now%cmb,units="m/yr ice equiv.",long_name="Calving mass balance rate", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"cmb_flt",ylmo%tpo%now%cmb_flt,units="m/yr ice equiv.",long_name="Floating potential calving mass balance rate", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"cmb_grnd",ylmo%tpo%now%cmb_grnd,units="m/yr ice equiv.",long_name="Grounded potential calving mass balance rate", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_ice",ylmo%tpo%now%f_ice,units="1",long_name="Total ice fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"f_ice_g",ylmo%tpo%now%f_ice*ylmo%tpo%now%f_grnd,units="1",long_name="Grounded ice fraction", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"f_ice_f",ylmo%tpo%now%f_ice*(1.0-ylmo%tpo%now%f_grnd),units="1",long_name="Floating ice fraction", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! == Additional variables 

        call nc_write(filename,"mask_bed",ylmo%tpo%now%mask_bed,units="",long_name="Bed mask", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_grz",ylmo%tpo%now%mask_grz,units="",long_name="Grounding-zone mask", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_frnt",ylmo%tpo%now%mask_frnt,units="",long_name="Ice-front mask", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_grnd",ylmo%tpo%now%f_grnd,units="1",long_name="Grounded ice fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"H_grnd",ylmo%tpo%now%H_grnd,units="m",long_name="Ice thickness above flotation (distance to flotation)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"dHidt",ylmo%tpo%now%dHidt,units="m/yr",long_name="Ice thickness rate of change", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dzsdt",ylmo%tpo%now%dzsdt,units="m/yr",long_name="Surface elevation rate of change", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"dist_grline",ylmo%tpo%now%dist_grline,units="km",long_name="Distance to grounding line", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dist_margin",ylmo%tpo%now%dist_margin,units="km",long_name="Distance to ice margin", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"f_grnd_acx",ylmo%tpo%now%f_grnd_acx,units="1",long_name="Grounded ice fraction (x)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_grnd_acy",ylmo%tpo%now%f_grnd_acy,units="1",long_name="Grounded ice fraction (y)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"f_grnd_pin",ylmo%tpo%now%f_grnd_pin,units="1",long_name="Grounded ice fraction (pinning points)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"kt",ylmo%tpo%now%kt,units="m yr-1 Pa-1",long_name="Von Mises calving scaling factor", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)


        call nc_write(filename,"uxy_s",ylmo%dyn%now%uxy_s,units="m/a",long_name="Surface velocity (magnitude)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_b",ylmo%dyn%now%uxy_b,units="m/a",long_name="Basal velocity (magnitude)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_bar",ylmo%dyn%now%uxy_bar,units="m/a",long_name="Vertically averaged velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)  
        
        call nc_write(filename,"duxydt",ylmo%dyn%now%duxydt,units="m a-2",long_name="Vertically averaged acceleration magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)  
        
        call nc_write(filename,"mb_net",ylmo%tpo%now%mb_net,units="m",long_name="Applied net mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mb_resid",ylmo%tpo%now%mb_resid,units="m",long_name="Residual mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mb_err",ylmo%tpo%now%mb_err,units="m",long_name="Residual error mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"bmb_grnd",ylmo%thrm%now%bmb_grnd,units="m/a ice equiv.",long_name="Basal mass balance (grounded)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"bmb_shlf",ylmo%bnd%bmb_shlf,units="m/a ice equiv.",long_name="Basal mass balance (shelf)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"N_eff",ylmo%dyn%now%N_eff,units="Pa",long_name="Effective pressure", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"cb_ref",ylmo%dyn%now%cb_ref,units="--",long_name="Bed friction scalar", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"c_bed",ylmo%dyn%now%c_bed,units="Pa",long_name="Bed friction coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"beta",ylmo%dyn%now%beta,units="Pa a m^-1",long_name="Basal friction coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"visc_eff_int",ylmo%dyn%now%visc_eff_int,units="Pa a m",long_name="Depth-integrated effective viscosity (SSA)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taud",ylmo%dyn%now%taud,units="Pa",long_name="Driving stress", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"enh_bar",ylmo%mat%now%enh_bar,units="1",long_name="Vertically averaged enhancement factor", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"de_eff",ylmo%dyn%now%de_eff,units="1/a",long_name="Effective strain rate", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        call nc_write(filename,"ssa_mask_acx",ylmo%dyn%now%ssa_mask_acx,units="1",long_name="Grounded ice fraction (x)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"ssa_mask_acy",ylmo%dyn%now%ssa_mask_acy,units="1",long_name="Grounded ice fraction (y)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"H_ice_dyn",ylmo%tpo%now%H_ice_dyn,units="1",long_name="Ice thickness (dynamics)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_ice_dyn",ylmo%tpo%now%f_ice_dyn,units="1",long_name="Ice fraction (dynamics)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"dzbdx",ylmo%tpo%now%dzbdx,units="1",long_name="Bed slope", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dzbdy",ylmo%tpo%now%dzbdy,units="1",long_name="Bed slope", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
                      
        ! Strain-rate and stress tensors 
        if (.TRUE.) then

            call nc_write(filename,"de",ylmo%mat%now%strn%de,units="a^-1",long_name="Effective strain rate", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
            call nc_write(filename,"te",ylmo%mat%now%strs%te,units="Pa",long_name="Effective stress", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
            call nc_write(filename,"visc_int",ylmo%mat%now%visc_int,units="Pa a m",long_name="Depth-integrated effective viscosity (SSA)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
         
            call nc_write(filename,"de2D",ylmo%mat%now%strn2D%de,units="yr^-1",long_name="Effective strain rate", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"div2D",ylmo%mat%now%strn2D%div,units="yr^-1",long_name="Divergence strain rate", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"te2D",ylmo%mat%now%strs2D%te,units="Pa",long_name="Effective stress", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            
            call nc_write(filename,"eps_eig_1",ylmo%mat%now%strn2D%eps_eig_1,units="1/yr",long_name="Eigen strain 1", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"eps_eig_2",ylmo%mat%now%strn2D%eps_eig_2,units="1/yr",long_name="Eigen strain 2", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"eps_eff",ylmo%tpo%now%eps_eff,units="yr^-1",long_name="Effective calving strain", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            
            call nc_write(filename,"tau_eig_1",ylmo%mat%now%strs2D%tau_eig_1,units="Pa",long_name="Eigen stress 1", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"tau_eig_2",ylmo%mat%now%strs2D%tau_eig_2,units="Pa",long_name="Eigen stress 2", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            call nc_write(filename,"tau_eff",ylmo%tpo%now%tau_eff,units="Pa",long_name="Effective calving stress", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            
        end if 

        call nc_write(filename,"uz_star",ylmo%dyn%now%uz_star,units="m yr-1",long_name="Advection-adjusted vertical velocity", &
                      dim1="xc",dim2="yc",dim3="zeta_ac",dim4="time",start=[1,1,1,n],ncid=ncid)
        
        call nc_write(filename,"T_rock",ylmo%thrm%now%T_rock,units="K",long_name="Bedrock temperature", &
                      dim1="xc",dim2="yc",dim3="zeta_rock",dim4="time",start=[1,1,1,n],ncid=ncid)
        
        call nc_write(filename,"Q_rock",ylmo%thrm%now%Q_rock,units="mW m-2",long_name="Bedrock surface heat flux", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"Q_ice_b",ylmo%thrm%now%Q_ice_b,units="mW m-2",long_name="Basal ice heat flux", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"Q_strn",ylmo%thrm%now%Q_strn/(ylmo%bnd%c%rho_ice*ylmo%thrm%now%cp),units="K a-1",long_name="Strain heating", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        call nc_write(filename,"Q_b",ylmo%thrm%now%Q_b,units="mW m-2",long_name="Basal frictional heating", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"beta_eff",ylmo%dyn%now%beta_eff,units="Pa a m^-1",long_name="Effective basal friction coefficient (DIVA)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"dep_time",ylmo%mat%now%dep_time,units="yr",long_name="Deposition time", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"uz",ylmo%dyn%now%uz,units="m/a",long_name="Vertical velocity (z)", &
                       dim1="xc",dim2="yc",dim3="zeta_ac",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"T_prime",ylmo%thrm%now%T_prime,units="deg C",long_name="Homologous ice temperature", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"T_prime_b",ylmo%thrm%now%T_prime_b,units="K",long_name="Basal homologous ice temperature", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_pmp",ylmo%thrm%now%f_pmp,units="1",long_name="Fraction of grid point at pmp", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"H_w",ylmo%thrm%now%H_w,units="m water equiv.",long_name="Basal water layer thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! Comparison with present-day 
        call nc_write(filename,"H_ice_pd_err",ylmo%dta%pd%err_H_ice,units="m",long_name="Ice thickness error wrt present day", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf_pd_err",ylmo%dta%pd%err_z_srf,units="m",long_name="Surface elevation error wrt present day", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_s_pd_err",ylmo%dta%pd%err_uxy_s,units="m/a",long_name="Surface velocity error wrt present day", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! Diagnostics 
        call nc_write(filename,"taud_acx",ylmo%dyn%now%taud_acx,units="Pa",long_name="Driving stress (x)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taud_acy",ylmo%dyn%now%taud_acy,units="Pa",long_name="Driving stress (y)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"taul_int_acx",ylmo%dyn%now%taul_int_acx,units="Pa m",long_name="Vertically integrated lateral stress (x)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taul_int_acy",ylmo%dyn%now%taul_int_acy,units="Pa m",long_name="Vertically integrated lateral stress (y)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
             
        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine write_step_2D

    subroutine modify_cb_ref(dyn,tpo,thrm,bnd,grd,domain,cb_ref)
        ! Modify cb_ref [unitless] with location specific tuning 

        implicit none
        
        type(ydyn_class),   intent(INOUT) :: dyn
        type(ytopo_class),  intent(IN)    :: tpo 
        type(ytherm_class), intent(IN)    :: thrm
        type(ybound_class), intent(IN)    :: bnd  
        type(ygrid_class),  intent(IN)    :: grd
        character(len=*),   intent(IN)    :: domain 
        real(prec),         intent(INOUT) :: cb_ref(:,:) 

        integer :: i, j, nx, ny 
        integer :: i1, i2, j1, j2 
        real(prec) :: f_scale 
            
            ! Additionally modify cb_ref
            if (trim(domain) .eq. "Antarctica") then


                ! Increase friction - feeding the Ronne ice shelf from the South
!                 call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*2.0,x0=-800.0, y0= 100.0,sigma=400.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
!                 call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*2.0,x0=-980.0, y0=-400.0,sigma=200.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
                
                ! Reduction friction - feeding the Ross ice shelf from the East
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.05_prec,x0= 130.0_prec, y0=-550.0_prec, sigma=100.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.05_prec,x0= 280.0_prec, y0=-760.0_prec, sigma=100.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.05_prec,x0= 380.0_prec, y0=-960.0_prec, sigma=100.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.05_prec,x0= 400.0_prec, y0=-1150.0_prec,sigma=100.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)


if (.FALSE.) then 
                ! Increase - feeding the Ronne ice shelf from the North
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*4.00_prec,x0=-700.0_prec, y0=    0.0_prec,sigma=200.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                
                ! Increase - Southeast Antarctica inland
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*2.00_prec,x0=1500.0_prec, y0= -550.0_prec,sigma=200.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*2.00_prec,x0=1700.0_prec, y0=-1000.0_prec,sigma=200.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                
                ! Reduction - South pole 
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.10_prec,x0=   0.0_prec, y0=   0.0_prec, sigma=400.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.20_prec,x0=   0.0_prec, y0= 600.0_prec, sigma=300.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.10_prec,x0= 500.0_prec, y0=-500.0_prec, sigma=400.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                
                ! Reduction - Amery ice shelf
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.05_prec,x0=1500.0_prec, y0= 650.0_prec, sigma=300.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                
                ! Reduction - feeding the Ross ice shelf from the North
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.005_prec,x0=-500.0_prec, y0=-500.0_prec, sigma=300.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)

                ! Reduction - feeding the Ross ice shelf from the East
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.05_prec,x0= 130.0_prec, y0=-550.0_prec, sigma=100.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.05_prec,x0= 280.0_prec, y0=-760.0_prec, sigma=100.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.05_prec,x0= 380.0_prec, y0=-960.0_prec, sigma=100.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.05_prec,x0= 400.0_prec, y0=-1150.0_prec,sigma=100.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)

end if 

if (.FALSE.) then
                ! Reduction 
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.25_prec,x0=-2000.0_prec,y0=1000.0_prec, sigma=100.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.20_prec,x0=-750.0_prec, y0=-900.0_prec, sigma=300.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.20_prec,x0=-600.0_prec, y0=-600.0_prec, sigma=300.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.20_prec,x0=-300.0_prec, y0=   0.0_prec, sigma=100.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.20_prec,x0=-250.0_prec, y0=-500.0_prec, sigma=100.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.20_prec,x0=-100.0_prec, y0=-600.0_prec, sigma=300.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.20_prec,x0=-100.0_prec, y0=-300.0_prec, sigma=300.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.20_prec,x0=   0.0_prec, y0=   0.0_prec, sigma=300.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.20_prec,x0= 130.0_prec, y0=-550.0_prec, sigma= 50.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.20_prec,x0= 280.0_prec, y0=-760.0_prec, sigma= 50.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.20_prec,x0= 380.0_prec, y0=-960.0_prec, sigma= 50.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.50_prec,x0= 400.0_prec, y0=   0.0_prec, sigma=300.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.20_prec,x0= 400.0_prec, y0=-1150.0_prec,sigma= 50.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.25_prec,x0= 700.0_prec, y0= -500.0_prec,sigma=400.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.20_prec,x0=1500.0_prec, y0= 650.0_prec, sigma=300.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                
                ! Increase
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*4.00_prec,x0=-600.0_prec, y0=    0.0_prec,sigma=200.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*2.00_prec,x0=1200.0_prec, y0=-1200.0_prec,sigma=300.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*1.50_prec,x0=2000.0_prec, y0=    0.0_prec,sigma=300.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
end if 

            end if 

        return 

    end subroutine modify_cb_ref

    subroutine scale_cf_gaussian(cb_ref,cf_new,x0,y0,sigma,xx,yy)

        implicit none 

        real(prec), intent(INOUT) :: cb_ref(:,:)
        real(prec), intent(IN) :: cf_new
        real(prec), intent(IN) :: x0
        real(prec), intent(IN) :: y0
        real(prec), intent(IN) :: sigma
        real(prec), intent(IN) :: xx(:,:)
        real(prec), intent(IN) :: yy(:,:)

        ! Local variables 
        integer :: nx, ny 
        real(prec), allocatable :: wts(:,:)
        
        nx = size(cb_ref,1)
        ny = size(cb_ref,2)

        allocate(wts(nx,ny))

        ! Get Gaussian weights 
        wts = 1.0/(2.0*pi*sigma**2)*exp(-((xx-x0)**2+(yy-y0)**2)/(2.0*sigma**2))
        wts = wts / maxval(wts)

        ! Scale cb_ref
        cb_ref = cb_ref*(1.0-wts) + cf_new*wts

        return 

    end subroutine scale_cf_gaussian

end program yelmo_test



