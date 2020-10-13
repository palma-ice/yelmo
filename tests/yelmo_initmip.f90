

program yelmo_test

    use ncio 
    use yelmo 
    
    use basal_dragging 
    use yelmo_tools, only : gauss_values
    
    use omp_lib

    implicit none 

    type(yelmo_class)       :: yelmo1 

    character(len=256) :: outfldr, file1D, file2D, file_restart, domain 
    character(len=512) :: path_par, path_const, clim_nm  
    real(prec) :: time_init, time_end, time_equil, time, dtt, dtt_equil, dt1D_out, dt2D_out 
    real(prec) :: dtt_equil_now 
    real(prec) :: bmb_shlf_const, dT_ann, z_sl    
    integer    :: n
    logical    :: with_anom 
    
    real(prec), allocatable :: cf_ref(:,:) 

    ! No-ice mask (to impose additional melting)
    logical, allocatable :: mask_noice(:,:)  

    ! cf_ref, bmelt  
    logical :: load_cf_ref, load_bmelt
    character(len=256) :: file_cf_ref 
    character(len=256) :: file_bmelt 

    real(8) :: cpu_start_time, cpu_end_time, cpu_dtime  
    
    ! Start timing 
    call yelmo_cpu_time(cpu_start_time)
    
    ! Determine the parameter file from the command line 
    call yelmo_load_command_line_args(path_par)

    ! Timing and other parameters 
    call nml_read(path_par,"ctrl","time_init",       time_init)                 ! [yr] Starting time
    call nml_read(path_par,"ctrl","time_end",        time_end)                  ! [yr] Ending time
    call nml_read(path_par,"ctrl","time_equil",      time_equil)                ! [yr] Years to equilibrate first
    call nml_read(path_par,"ctrl","dtt",             dtt)                       ! [yr] Main loop time step 
    call nml_read(path_par,"ctrl","dtt_equil",       dtt_equil)                 ! [yr] Timestep to use for dynamic equilibration (if dt_method=0)
    call nml_read(path_par,"ctrl","dt1D_out",        dt1D_out)                  ! [yr] Frequency of 1D output 
    call nml_read(path_par,"ctrl","dt2D_out",        dt2D_out)                  ! [yr] Frequency of 2D output 
    
    call nml_read(path_par,"ctrl","clim_nm",         clim_nm)                   ! Namelist group holding climate information
    call nml_read(path_par,"ctrl","with_anom",       with_anom)                 ! Apply anomaly at the start of the simulation (after equilibration)
    
    call nml_read(path_par,"ctrl","load_cf_ref",     load_cf_ref)               ! Load cf_ref from file? Otherwise define from cf_stream + inline tuning
    call nml_read(path_par,"ctrl","file_cf_ref",     file_cf_ref)               ! Filename holding cf_ref to load 

    call nml_read(path_par,"ctrl","load_bmelt",      load_bmelt)                ! Load bmelt from file?
    call nml_read(path_par,"ctrl","file_bmelt",      file_bmelt)                ! Filename holding bmelt field to load 

    ! Load climate (eg, clim_pd or clim_lgm)
    call nml_read(path_par,clim_nm,  "bmb_shlf_const",  bmb_shlf_const)            ! [yr] Constant imposed bmb_shlf value
    call nml_read(path_par,clim_nm,  "dT_ann",          dT_ann)                    ! [K] Temperature anomaly (atm)
    call nml_read(path_par,clim_nm,  "z_sl",            z_sl)                      ! [m] Sea level relative to present-day

    ! Assume program is running from the output folder
    outfldr = "./"

    ! Define input and output locations 
    path_const   = trim(outfldr)//"yelmo_const_Earth.nml"
    file1D       = trim(outfldr)//"yelmo1D.nc"
    file2D       = trim(outfldr)//"yelmo2D.nc"
    file_restart = trim(outfldr)//"yelmo_restart.nc"

    ! === Initialize ice sheet model =====

    ! General initialization of yelmo constants (used globally)
    call yelmo_global_init(path_const)

    ! Initialize data objects and load initial topography
    call yelmo_init(yelmo1,filename=path_par,grid_def="file",time=time_init)

    ! Determine maximum timestep to use with equilibration step 
    if (yelmo1%par%dt_method .eq. 0) then 
        dtt_equil_now = dtt_equil 
    else 
        dtt_equil_now = dtt 
    end if 
    
    ! === Set initial boundary conditions for current time and yelmo state =====
    ! ybound: z_bed, z_sl, H_sed, smb, T_srf, bmb_shlf , Q_geo

    yelmo1%bnd%z_sl     = z_sl              ! [m]
    yelmo1%bnd%H_sed    = 0.0               ! [m]
    yelmo1%bnd%Q_geo    = 50.0              ! [mW/m2]
    
    ! Impose present-day surface mass balance and present-day temperature field plus any anomaly
    yelmo1%bnd%smb      = yelmo1%dta%pd%smb             ! [m.i.e./a]
    yelmo1%bnd%T_srf    = yelmo1%dta%pd%T_srf + dT_ann  ! [K]
    
    if (load_bmelt) then

        ! Parse filenames with grid information
        call yelmo_parse_path(file_bmelt,yelmo1%par%domain,yelmo1%par%grid_name)
 
        call nc_read(file_bmelt,"bm_ac_reese",yelmo1%bnd%bmb_shlf)
        yelmo1%bnd%bmb_shlf = -yelmo1%bnd%bmb_shlf                  ! Negative because bmb = -bmelt 
    else 
        yelmo1%bnd%bmb_shlf = bmb_shlf_const    ! [m.i.e./a]
    end if 

    yelmo1%bnd%T_shlf   = T0                ! [K]   

    if (dT_ann .lt. 0.0) yelmo1%bnd%T_shlf   = T0 + dT_ann*0.25_prec  ! [K] Oceanic temp anomaly
    
    call yelmo_print_bound(yelmo1%bnd)

    time = time_init 

    ! Define no-ice mask from present-day data
    allocate(mask_noice(yelmo1%grd%nx,yelmo1%grd%ny))
    mask_noice = .FALSE. 
    ! Present-day
    if (dT_ann .ge. 0.0) then
        where(yelmo1%dta%pd%H_ice .le. 0.0) mask_noice = .TRUE. 
    end if 

    ! Special treatment for Greenland
    if (trim(yelmo1%par%domain) .eq. "Greenland") then 
        
        ! Present-day
        if (dT_ann .ge. 0.0) then 
            ! Impose additional negative mass balance to no ice points of 2 [m.i.e./a] melting
            where(mask_noice) yelmo1%bnd%smb = yelmo1%dta%pd%smb - 2.0 
        end if 

    end if 

    ! Special treatment for Antarctica
    if (trim(yelmo1%par%domain) .eq. "Antarctica") then 
        
        ! Present-day
        if (dT_ann .ge. 0.0) then 
            where(mask_noice) yelmo1%bnd%bmb_shlf = -2.0                ! [m/a]        
        end if 

        ! LGM
        if (dT_ann .lt. 0.0) then 
            where(yelmo1%bnd%smb .le. 0.1) yelmo1%bnd%smb = 0.1         ! [m/a]
        end if 

        ! Present-day and LGM 
        where(yelmo1%bnd%regions .eq. 2.0) yelmo1%bnd%smb = -1.0        ! [m/a]

    end if 

    ! Initialize state variables (dyn,therm,mat)
    ! (initialize temps with robin method with a cold base)
    call yelmo_init_state(yelmo1,path_par,time=time_init,thrm_method="robin-cold")

    ! ============================================================
    ! Load or define cf_ref 

    ! Allocate cf_ref and set it to cf_stream by default 
    allocate(cf_ref(yelmo1%grd%nx,yelmo1%grd%ny))
    cf_ref = yelmo1%dyn%par%cf_stream  

    if (load_cf_ref) then 

        ! Parse filename with grid information
        call yelmo_parse_path(file_cf_ref,yelmo1%par%domain,yelmo1%par%grid_name)

        ! Load cf_ref from specified file 
        call nc_read(file_cf_ref,"cf_ref",cf_ref)

        ! Make sure that present-day shelves have minimal cf_ref values 
        where(yelmo1%tpo%now%f_grnd .eq. 0.0) cf_ref = 0.05 

    else
        ! Define cf_ref inline 

        cf_ref = yelmo1%dyn%par%cf_stream  

    end if 

    ! ============================================================


    ! Define cf_ref initially
    call calc_ydyn_cfref_external(yelmo1%dyn,yelmo1%tpo,yelmo1%thrm,yelmo1%bnd,yelmo1%grd,yelmo1%par%domain, &
                                    mask_noice,cf_ref)

    ! Impose a colder boundary temperature for equilibration step 
    ! -5 [K] for mimicking glacial times
!     yelmo1%bnd%T_srf = yelmo1%dta%pd%T_srf - 10.0  

if (.FALSE.) then 
    ! Run yelmo for several years with constant boundary conditions and topo
    ! to equilibrate thermodynamics and dynamics
    yelmo1%par%dt_method        = 0 
    yelmo1%tpo%par%topo_rel     = 2
    yelmo1%tpo%par%topo_rel_tau = 10.0 
    write(*,*) "timelog, tau = ", yelmo1%tpo%par%topo_rel_tau
    call yelmo_update_equil(yelmo1,time,time_tot=1.0_prec,topo_fixed=.FALSE.,dt=0.05,ssa_vel_max=5000.0_prec)
    yelmo1%tpo%par%topo_rel_tau = 100.0 
    write(*,*) "timelog, tau = ", yelmo1%tpo%par%topo_rel_tau
    call yelmo_update_equil(yelmo1,time,time_tot=2.0_prec,topo_fixed=.FALSE.,dt=0.05,ssa_vel_max=5000.0_prec)
    yelmo1%tpo%par%topo_rel_tau = 1000.0 
    write(*,*) "timelog, tau = ", yelmo1%tpo%par%topo_rel_tau
    call yelmo_update_equil(yelmo1,time,time_tot=5.0_prec,topo_fixed=.FALSE.,dt=0.05,ssa_vel_max=5000.0_prec)
    yelmo1%tpo%par%topo_rel     = 0
    write(*,*) "timelog, no relaxation..."
    call yelmo_update_equil(yelmo1,time,time_tot=5.0_prec,topo_fixed=.FALSE.,dt=0.05,ssa_vel_max=5000.0_prec)
    yelmo1%par%dt_method        = 2
    write(*,*) "timelog, adaptive timestepping..."
    call yelmo_update_equil(yelmo1,time,time_tot=10.0_prec,topo_fixed=.FALSE.,dt=min(1.0_prec,dtt_equil_now),ssa_vel_max=5000.0_prec)

    ! Next equilibrate thermodynamics and maintain constant ice topopgraphy (for speed)
    call yelmo_update_equil(yelmo1,time,time_tot=time_equil,topo_fixed=.TRUE.,dt=1.0_prec,ssa_vel_max=5000.0_prec)

    ! Finally allow further dynamic equilibrium (therm + dyn) to ensure 
    ! everything is in sync
    call yelmo_update_equil(yelmo1,time,time_tot=100.0_prec,topo_fixed=.FALSE.,dt=dtt_equil_now,ssa_vel_max=5000.0_prec)

else 
    ! Just testing...

    ! Calculate dynamics and thermodynamics, constant ice thickness
    call yelmo_update_equil(yelmo1,time,time_tot=10.0,topo_fixed=.TRUE.,dt=1.0_prec,ssa_vel_max=5000.0_prec)

    ! Now let it advance one timestep with no smb/bmb
    yelmo1%bnd%smb      = 0.0_prec 
    yelmo1%bnd%bmb_shlf = 0.0_prec 
    
    call yelmo_update_equil(yelmo1,time,time_tot=1.0_prec,topo_fixed=.FALSE.,dt=0.1_prec,ssa_vel_max=5000.0_prec)

    !Re impose smb/bmb
    yelmo1%bnd%smb      = yelmo1%dta%pd%smb             ! [m.i.e./a]
    
    if (load_bmelt) then

        ! Parse filenames with grid information
        call yelmo_parse_path(file_bmelt,yelmo1%par%domain,yelmo1%par%grid_name)
 
        call nc_read(file_bmelt,"bm_ac_reese",yelmo1%bnd%bmb_shlf)
        yelmo1%bnd%bmb_shlf = -yelmo1%bnd%bmb_shlf                  ! Negative because bmb = -bmelt 
    else 
        yelmo1%bnd%bmb_shlf = bmb_shlf_const    ! [m.i.e./a]
    end if 

end if 

    ! 2D file 
    call yelmo_write_init(yelmo1,file2D,time_init=time,units="years")  
    call write_step_2D(yelmo1,file2D,time=time,cf_ref=cf_ref)
    
    ! 1D file 
    call write_yreg_init(yelmo1,file1D,time_init=time_init,units="years",mask=yelmo1%bnd%ice_allowed)
    call write_yreg_step(yelmo1%reg,file1D,time=time)  

    write(*,*) "Completed test."
    stop 

    if (with_anom) then 
        ! Warm up the ice sheet to impose some changes 
        yelmo1%bnd%T_srf    = yelmo1%dta%pd%T_srf + 5.0
        yelmo1%bnd%bmb_shlf = -10.0 
        where (yelmo1%dta%pd%smb .lt. 0.0) yelmo1%bnd%smb = yelmo1%dta%pd%smb - 1.0 
    end if 

    ! Advance timesteps
    do n = 1, ceiling((time_end-time_init)/dtt)

        ! Get current time 
        time = time_init + n*dtt

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
        
        ! Update ice sheet 
        call yelmo_update(yelmo1,time)

        ! == MODEL OUTPUT =======================================================

        if (mod(nint(time*100),nint(dt2D_out*100))==0) then
            call write_step_2D(yelmo1,file2D,time=time,cf_ref=cf_ref)
        end if 

        if (mod(nint(time*100),nint(dt1D_out*100))==0) then 
            call write_yreg_step(yelmo1%reg,file1D,time=time) 
        end if 

        if (mod(time,10.0)==0 .and. (.not. yelmo_log)) then
            write(*,"(a,f14.4)") "yelmo:: time = ", time
        end if 
        
    end do 
    ! == Finished time loop == 

    ! Write a final restart file 
    call yelmo_restart_write(yelmo1,file_restart,time)

    ! Finalize program
    call yelmo_end(yelmo1,time=time)
    
    ! Stop timing 
    call yelmo_cpu_time(cpu_end_time,cpu_start_time,cpu_dtime)

    write(*,"(a,f12.3,a)") "Time  = ",cpu_dtime/60.0 ," min"
    write(*,"(a,f12.1,a)") "Speed = ",(1e-3*(time_end-time_init))/(cpu_dtime/3600.0), " kiloyears / hr"

contains

    subroutine write_step_2D(ylmo,filename,time,cf_ref)

        implicit none 
        
        type(yelmo_class), intent(IN) :: ylmo
        character(len=*),  intent(IN) :: filename
        real(prec), intent(IN) :: time
        real(prec), intent(IN) :: cf_ref(:,:) 

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
        call nc_write(filename,"dHicedt",ylmo%tpo%now%dHicedt,units="m/a",long_name="Ice thickness change", &
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
!         call nc_write(filename,"ux_bar",ylmo%dyn%now%ux_bar,units="m/a",long_name="Vertically averaged velocity (x)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"uy_bar",ylmo%dyn%now%uy_bar,units="m/a",long_name="Vertically averaged velocity (y)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"T_ice_s",ylmo%thrm%now%T_ice(:,:,ylmo%par%nz_aa),units="K",long_name="Surface ice temperature", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"T_ice_b",ylmo%thrm%now%T_ice(:,:,1),units="K",long_name="Basal homologous ice temperature", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taub",ylmo%dyn%now%taub,units="Pa",long_name="Basal dragging stress (magnitude)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"calv",ylmo%tpo%now%calv,units="m/a ice equiv.",long_name="Calving rate", &
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
        call nc_write(filename,"f_grnd",ylmo%tpo%now%f_grnd,units="1",long_name="Grounded ice fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"H_margin",ylmo%tpo%now%H_margin,units="m",long_name="Margin ice thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"uxy_s",ylmo%dyn%now%uxy_s,units="m/a",long_name="Surface velocity (magnitude)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_b",ylmo%dyn%now%uxy_b,units="m/a",long_name="Basal velocity (magnitude)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_bar",ylmo%dyn%now%uxy_bar,units="m/a",long_name="Vertically averaged velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)  
        
        call nc_write(filename,"duxydt",ylmo%dyn%now%duxydt,units="m a-2",long_name="Vertically averaged acceleration magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)  
        
        call nc_write(filename,"mb_applied",ylmo%tpo%now%mb_applied,units="m",long_name="Applied net mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"bmb_grnd",ylmo%thrm%now%bmb_grnd,units="m/a ice equiv.",long_name="Basal mass balance (grounded)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"bmb_shlf",ylmo%bnd%bmb_shlf,units="m/a ice equiv.",long_name="Basal mass balance (shelf)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"N_eff",ylmo%dyn%now%N_eff,units="Pa",long_name="Effective pressure", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"cf_ref",ylmo%dyn%now%cf_ref,units="--",long_name="Bed friction scalar", &
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
        
        call nc_write(filename,"beta_eff",ylmo%dyn%now%beta_eff,units="Pa a m^-1",long_name="Effective basal friction coefficient (DIVA)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"beta_diva",ylmo%dyn%now%beta_diva,units="Pa a m^-1",long_name="Actual basal friction coefficient (DIVA)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"dep_time",ylmo%mat%now%dep_time,units="yr",long_name="Deposition time", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"uz",ylmo%dyn%now%uz,units="m/a",long_name="Surface velocity (z)", &
                       dim1="xc",dim2="yc",dim3="zeta_ac",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"T_prime",ylmo%thrm%now%T_ice-ylmo%thrm%now%T_pmp,units="deg C",long_name="Homologous ice temperature", &
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
!         call nc_write(filename,"taud_acx",ylmo%dyn%now%taud_acx,units="Pa",long_name="Driving stress (x)", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"taud_acy",ylmo%dyn%now%taud_acy,units="Pa",long_name="Driving stress (y)", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine write_step_2D

    subroutine calc_ydyn_cfref_external(dyn,tpo,thrm,bnd,grd,domain,mask_noice,cf_ref)
        ! Update cfref [Pa] based on parameter choices

        implicit none
        
        type(ydyn_class),   intent(INOUT) :: dyn
        type(ytopo_class),  intent(IN)    :: tpo 
        type(ytherm_class), intent(IN)    :: thrm
        type(ybound_class), intent(IN)    :: bnd  
        type(ygrid_class),  intent(IN)    :: grd
        character(len=*),   intent(IN)    :: domain 
        logical,            intent(IN)    :: mask_noice(:,:) 
        real(prec),         intent(INOUT) :: cf_ref(:,:) 

        integer :: i, j, nx, ny 
        integer :: i1, i2, j1, j2 
        real(prec) :: f_scale 
        
        real(prec), allocatable :: lambda_bed(:,:)  

        nx = size(dyn%now%cf_ref,1)
        ny = size(dyn%now%cf_ref,2)
        
        allocate(lambda_bed(nx,ny))

            ! cf_ref [unitless] is obtained as input to this routine from optimization or elsewhere

            ! =============================================================================
            ! Step 2: calculate lambda functions to scale cf_ref from default value 
            
            select case(trim(dyn%par%cb_scale))

                case("lin_zb")
                    ! Linear scaling function with bedrock elevation
                    
                    lambda_bed = calc_lambda_bed_lin(bnd%z_bed,dyn%par%cb_z0,dyn%par%cb_z1)

                case("exp_zb")
                    ! Exponential scaling function with bedrock elevation
                    
                    ! Default
                    lambda_bed = calc_lambda_bed_exp(bnd%z_bed,dyn%par%cb_z0,dyn%par%cb_z1)

                    if (trim(domain) .eq. "Antarctica") then 
                        ! Domain-specific modifications to lambda function

                        ! Increased friction in Wilkes Land (South - Southeast)
                        where (bnd%basins .ge. 12 .and. bnd%basins .le. 17) 
                            lambda_bed = calc_lambda_bed_exp(bnd%z_bed,-400.0_prec,dyn%par%cb_z1)
                        end where

                        ! Increased friction in WAIS divide area feeding the Ronne
                        where (bnd%basins .eq.  1 .or. bnd%basins .eq.  2)
                            lambda_bed = calc_lambda_bed_exp(bnd%z_bed,-2000.0_prec,dyn%par%cb_z1)
                        end where 

!                         ! Increased friction in WAIS divide area feeding the Ross
!                         where (bnd%basins .eq. 21 .or. bnd%basins .eq. 22)
!                             lambda_bed = calc_lambda_bed_exp(bnd%z_bed,-2000.0,dyn%par%cb_z1)
!                         end where

                    end if 

                case("till_const")
                    ! Constant till friction angle

                    lambda_bed = calc_lambda_till_const(dyn%par%till_phi_const)

                case("till_zb")
                    ! Linear till friction angle versus elevation

                    lambda_bed = calc_lambda_till_linear(bnd%z_bed,bnd%z_sl,dyn%par%till_phi_min,dyn%par%till_phi_max, &
                                                            dyn%par%till_phi_zmin,dyn%par%till_phi_zmax)

                case DEFAULT
                    ! No scaling

                    lambda_bed = 1.0

            end select 

            ! Set lambda_bed to lower limit for regions of noice 
            where (mask_noice) lambda_bed = dyn%par%cb_min 
            
            ! Ensure lambda_bed is not below lower limit [default range 0:1] 
            where (lambda_bed .lt. dyn%par%cb_min) lambda_bed = dyn%par%cb_min

            ! =============================================================================
            ! Step 3: calculate cf_ref [--]
            
            cf_ref = (cf_ref*lambda_bed)

            ! Use location-specific tuning functions to modify cf_ref 
            call modify_cf_ref(dyn,tpo,thrm,bnd,grd,domain,cf_ref)
            
            ! Finally store in dyn object for output
            dyn%now%cf_ref = cf_ref 

        return 

    end subroutine calc_ydyn_cfref_external

    subroutine modify_cf_ref(dyn,tpo,thrm,bnd,grd,domain,cf_ref)
        ! Modify cf_ref [unitless] with location specific tuning 

        implicit none
        
        type(ydyn_class),   intent(INOUT) :: dyn
        type(ytopo_class),  intent(IN)    :: tpo 
        type(ytherm_class), intent(IN)    :: thrm
        type(ybound_class), intent(IN)    :: bnd  
        type(ygrid_class),  intent(IN)    :: grd
        character(len=*),   intent(IN)    :: domain 
        real(prec),         intent(INOUT) :: cf_ref(:,:) 

        integer :: i, j, nx, ny 
        integer :: i1, i2, j1, j2 
        real(prec) :: f_scale 
            
            ! Additionally modify cf_ref
            if (trim(domain) .eq. "Antarctica") then


                ! Increase friction - feeding the Ronne ice shelf from the South
!                 call scale_cf_gaussian(cf_ref,dyn%par%cf_stream*2.0,x0=-800.0, y0= 100.0,sigma=400.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
!                 call scale_cf_gaussian(cf_ref,dyn%par%cf_stream*2.0,x0=-980.0, y0=-400.0,sigma=200.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
                
                ! Reduction friction - feeding the Ross ice shelf from the East
                call scale_cf_gaussian(cf_ref,dyn%par%cf_stream*0.05_prec,x0= 130.0_prec, y0=-550.0_prec, sigma=100.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cf_ref,dyn%par%cf_stream*0.05_prec,x0= 280.0_prec, y0=-760.0_prec, sigma=100.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cf_ref,dyn%par%cf_stream*0.05_prec,x0= 380.0_prec, y0=-960.0_prec, sigma=100.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cf_ref,dyn%par%cf_stream*0.05_prec,x0= 400.0_prec, y0=-1150.0_prec,sigma=100.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)


if (.FALSE.) then 
                ! Increase - feeding the Ronne ice shelf from the North
                call scale_cf_gaussian(cf_ref,dyn%par%cf_stream*4.00_prec,x0=-700.0_prec, y0=    0.0_prec,sigma=200.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                
                ! Increase - Southeast Antarctica inland
                call scale_cf_gaussian(cf_ref,dyn%par%cf_stream*2.00_prec,x0=1500.0_prec, y0= -550.0_prec,sigma=200.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cf_ref,dyn%par%cf_stream*2.00_prec,x0=1700.0_prec, y0=-1000.0_prec,sigma=200.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                
                ! Reduction - South pole 
                call scale_cf_gaussian(cf_ref,dyn%par%cf_stream*0.10_prec,x0=   0.0_prec, y0=   0.0_prec, sigma=400.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cf_ref,dyn%par%cf_stream*0.20_prec,x0=   0.0_prec, y0= 600.0_prec, sigma=300.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cf_ref,dyn%par%cf_stream*0.10_prec,x0= 500.0_prec, y0=-500.0_prec, sigma=400.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                
                ! Reduction - Amery ice shelf
                call scale_cf_gaussian(cf_ref,dyn%par%cf_stream*0.05_prec,x0=1500.0_prec, y0= 650.0_prec, sigma=300.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                
                ! Reduction - feeding the Ross ice shelf from the North
                call scale_cf_gaussian(cf_ref,dyn%par%cf_stream*0.005_prec,x0=-500.0_prec, y0=-500.0_prec, sigma=300.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)

                ! Reduction - feeding the Ross ice shelf from the East
                call scale_cf_gaussian(cf_ref,dyn%par%cf_stream*0.05_prec,x0= 130.0_prec, y0=-550.0_prec, sigma=100.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cf_ref,dyn%par%cf_stream*0.05_prec,x0= 280.0_prec, y0=-760.0_prec, sigma=100.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cf_ref,dyn%par%cf_stream*0.05_prec,x0= 380.0_prec, y0=-960.0_prec, sigma=100.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cf_ref,dyn%par%cf_stream*0.05_prec,x0= 400.0_prec, y0=-1150.0_prec,sigma=100.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)

end if 

if (.FALSE.) then
                ! Reduction 
                call scale_cf_gaussian(cf_ref,dyn%par%cf_stream*0.25_prec,x0=-2000.0_prec,y0=1000.0_prec, sigma=100.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cf_ref,dyn%par%cf_stream*0.20_prec,x0=-750.0_prec, y0=-900.0_prec, sigma=300.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cf_ref,dyn%par%cf_stream*0.20_prec,x0=-600.0_prec, y0=-600.0_prec, sigma=300.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cf_ref,dyn%par%cf_stream*0.20_prec,x0=-300.0_prec, y0=   0.0_prec, sigma=100.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cf_ref,dyn%par%cf_stream*0.20_prec,x0=-250.0_prec, y0=-500.0_prec, sigma=100.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cf_ref,dyn%par%cf_stream*0.20_prec,x0=-100.0_prec, y0=-600.0_prec, sigma=300.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cf_ref,dyn%par%cf_stream*0.20_prec,x0=-100.0_prec, y0=-300.0_prec, sigma=300.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cf_ref,dyn%par%cf_stream*0.20_prec,x0=   0.0_prec, y0=   0.0_prec, sigma=300.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cf_ref,dyn%par%cf_stream*0.20_prec,x0= 130.0_prec, y0=-550.0_prec, sigma= 50.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cf_ref,dyn%par%cf_stream*0.20_prec,x0= 280.0_prec, y0=-760.0_prec, sigma= 50.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cf_ref,dyn%par%cf_stream*0.20_prec,x0= 380.0_prec, y0=-960.0_prec, sigma= 50.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cf_ref,dyn%par%cf_stream*0.50_prec,x0= 400.0_prec, y0=   0.0_prec, sigma=300.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cf_ref,dyn%par%cf_stream*0.20_prec,x0= 400.0_prec, y0=-1150.0_prec,sigma= 50.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cf_ref,dyn%par%cf_stream*0.25_prec,x0= 700.0_prec, y0= -500.0_prec,sigma=400.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cf_ref,dyn%par%cf_stream*0.20_prec,x0=1500.0_prec, y0= 650.0_prec, sigma=300.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                
                ! Increase
                call scale_cf_gaussian(cf_ref,dyn%par%cf_stream*4.00_prec,x0=-600.0_prec, y0=    0.0_prec,sigma=200.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cf_ref,dyn%par%cf_stream*2.00_prec,x0=1200.0_prec, y0=-1200.0_prec,sigma=300.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
                call scale_cf_gaussian(cf_ref,dyn%par%cf_stream*1.50_prec,x0=2000.0_prec, y0=    0.0_prec,sigma=300.0_prec,xx=grd%x*1e-3_prec,yy=grd%y*1e-3_prec)
end if 

            end if 

        return 

    end subroutine modify_cf_ref

    subroutine scale_cf_gaussian(cf_ref,cf_new,x0,y0,sigma,xx,yy)

        implicit none 

        real(prec), intent(INOUT) :: cf_ref(:,:)
        real(prec), intent(IN) :: cf_new
        real(prec), intent(IN) :: x0
        real(prec), intent(IN) :: y0
        real(prec), intent(IN) :: sigma
        real(prec), intent(IN) :: xx(:,:)
        real(prec), intent(IN) :: yy(:,:)

        ! Local variables 
        integer :: nx, ny 
        real(prec), allocatable :: wts(:,:)
        
        nx = size(cf_ref,1)
        ny = size(cf_ref,2)

        allocate(wts(nx,ny))

        ! Get Gaussian weights 
        wts = 1.0/(2.0*pi*sigma**2)*exp(-((xx-x0)**2+(yy-y0)**2)/(2.0*sigma**2))
        wts = wts / maxval(wts)

        ! Scale cf_ref
        cf_ref = cf_ref*(1.0-wts) + cf_new*wts

        return 

    end subroutine scale_cf_gaussian

end program yelmo_test 



