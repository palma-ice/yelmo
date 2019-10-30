

program yelmo_test

    use ncio 
    use yelmo 
    
    use basal_hydro_simple 
    use basal_dragging 
    use yelmo_tools, only : gauss_values
    
    implicit none 

    type(yelmo_class)      :: yelmo1
    type(hydro_class)      :: hyd1 

    character(len=256) :: outfldr, file1D, file2D, file_restart, domain 
    character(len=512) :: path_par, path_const  
    real(prec) :: time_init, time_end, time_equil, time, dtt, dt1D_out, dt2D_out 
    real(prec) :: bmb_shlf_const, dT_ann, z_sl    
    integer    :: n
    real(4) :: cpu_start_time, cpu_end_time 

    real(prec), allocatable :: cf_ref(:,:) 

    ! No-ice mask (to impose additional melting)
    logical, allocatable :: mask_noice(:,:)  

    ! cf_ref 
    logical :: load_cf_ref
    character(len=256) :: file_cf_ref 

    ! Start timing 
    call cpu_time(cpu_start_time)

    ! Determine the parameter file from the command line 
    call yelmo_load_command_line_args(path_par)

    ! Timing and other parameters 
    call nml_read(path_par,"control","time_init",       time_init)                 ! [yr] Starting time
    call nml_read(path_par,"control","time_end",        time_end)                  ! [yr] Ending time
    call nml_read(path_par,"control","time_equil",      time_equil)                ! [yr] Years to equilibrate first
    call nml_read(path_par,"control","dtt",             dtt)                       ! [yr] Main loop time step 
    call nml_read(path_par,"control","dt1D_out",        dt1D_out)                  ! [yr] Frequency of 1D output 
    call nml_read(path_par,"control","dt2D_out",        dt2D_out)                  ! [yr] Frequency of 2D output 
    call nml_read(path_par,"control","bmb_shlf_const",  bmb_shlf_const)            ! [yr] Constant imposed bmb_shlf value
    call nml_read(path_par,"control","dT_ann",          dT_ann)                    ! [K] Temperature anomaly (atm)
    call nml_read(path_par,"control","z_sl",            z_sl)                      ! [m] Sea level relative to present-day

    call nml_read(path_par,"control","load_cf_ref",     load_cf_ref)               ! Load cf_ref from file? Otherwise define from cf_stream + inline tuning
    call nml_read(path_par,"control","file_cf_ref",     file_cf_ref)               ! Filename holding cf_ref to load 

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

    ! Also intialize simple basal hydrology object
    call hydro_init(hyd1,filename=path_par,nx=yelmo1%grd%nx,ny=yelmo1%grd%ny)
    call hydro_init_state(hyd1,yelmo1%tpo%now%H_ice,yelmo1%tpo%now%f_grnd,time)

    ! === Set initial boundary conditions for current time and yelmo state =====
    ! ybound: z_bed, z_sl, H_sed, H_w, smb, T_srf, bmb_shlf , Q_geo

    yelmo1%bnd%z_sl     = z_sl              ! [m]
    yelmo1%bnd%H_sed    = 0.0               ! [m]
    yelmo1%bnd%H_w      = hyd1%now%H_w      ! [m]
    yelmo1%bnd%Q_geo    = 50.0              ! [mW/m2]
    
    ! Impose present-day surface mass balance and present-day temperature field plus any anomaly
    yelmo1%bnd%smb      = yelmo1%dta%pd%smb             ! [m.i.e./a]
    yelmo1%bnd%T_srf    = yelmo1%dta%pd%T_srf + dT_ann  ! [K]
    
    yelmo1%bnd%bmb_shlf = bmb_shlf_const    ! [m.i.e./a]
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
        
        ! Impose additional negative mass balance to no ice points of 2 [m.i.e./a] melting
        where(mask_noice) yelmo1%bnd%smb = yelmo1%dta%pd%smb - 2.0 
    
    end if 

    ! Special treatment for Antarctica
    if (trim(yelmo1%par%domain) .eq. "Antarctica") then 
        
        ! Present-day
        if (dT_ann .ge. 0.0) then 
            where(mask_noice) yelmo1%bnd%bmb_shlf = -2.0    ! [m/a]
!             where(yelmo1%bnd%basins .ge. 23.0 .and. & 
!                   yelmo1%bnd%basins .le. 26.0) yelmo1%bnd%bmb_shlf = -2.0   ! [m/a]
!             where(yelmo1%bnd%basins .ge.  9.0 .and. & 
!                   yelmo1%bnd%basins .le. 11.0) yelmo1%bnd%bmb_shlf = -1.0   ! [m/a]
        
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

        ! Use location-specific tuning functions to modify cf_ref 
        call modify_cf_ref(yelmo1%dyn,yelmo1%tpo,yelmo1%thrm,yelmo1%bnd,yelmo1%grd,yelmo1%par%domain,cf_ref)

    end if 

    ! ============================================================


    ! Define C_bed initially
    call calc_ydyn_cbed_external(yelmo1%dyn,yelmo1%tpo,yelmo1%thrm,yelmo1%bnd,yelmo1%grd,yelmo1%par%domain, &
                                    mask_noice,cf_ref)

    ! Impose a colder boundary temperature for equilibration step 
    ! -5 [K] for mimicking glacial times
!     yelmo1%bnd%T_srf = yelmo1%dta%pd%T_srf - 10.0  
    
    ! Run yelmo for several years with constant boundary conditions and topo
    ! to equilibrate thermodynamics and dynamics
    call yelmo_update_equil(yelmo1,time,time_tot=time_equil,topo_fixed=.FALSE.,dt=1.0,ssa_vel_max=500.0)
    
    ! 2D file 
    call yelmo_write_init(yelmo1,file2D,time_init=time,units="years")  
    call write_step_2D(yelmo1,file2D,time=time,cf_ref=cf_ref)
    
    ! 1D file 
    call write_yreg_init(yelmo1,file1D,time_init=time_init,units="years",mask=yelmo1%bnd%ice_allowed)
    call write_yreg_step(yelmo1%reg,file1D,time=time)  
    
    ! Advance timesteps
    do n = 1, ceiling((time_end-time_init)/dtt)

        ! Get current time 
        time = time_init + n*dtt

!         ! Update temperature and smb as needed in time (ISMIP6)
!         if (time .ge. -10e6 .and. time .lt. -10e3) then 
!             ! Glacial period, impose cold climate 
!             yelmo1%bnd%T_srf = yelmo1%dta%pd%T_srf - 10.0 

!         else if (time .ge. -10e3 .and. time .lt. -8e3) then
!             ! Holocene optimum 
!             yelmo1%bnd%T_srf = yelmo1%dta%pd%T_srf + 1.0 

!         else if (time .ge. -8e3) then 
!             ! Entering Holocene, impose present-day temperatures 
!             yelmo1%bnd%T_srf = yelmo1%dta%pd%T_srf
!         end if 

        ! Update ice sheet 
        call yelmo_update(yelmo1,time)

        ! Update basal hydrology 
        call hydro_update(hyd1,yelmo1%tpo%now%H_ice,yelmo1%tpo%now%f_grnd, &
                    -yelmo1%thrm%now%bmb_grnd*rho_ice/rho_w,time)

        ! Pass updated boundary variables to yelmo 
        yelmo1%bnd%H_w = hyd1%now%H_w 

        ! Update C_bed (due to effective pressure)
        call calc_ydyn_cbed_external(yelmo1%dyn,yelmo1%tpo,yelmo1%thrm,yelmo1%bnd,yelmo1%grd,yelmo1%par%domain, &
                                        mask_noice,cf_ref)

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


    ! Finalize program
    call yelmo_end(yelmo1,time=time)

    ! Stop timing 
    call cpu_time(cpu_end_time)


    write(*,"(a,f12.3,a)") "Time  = ",(cpu_end_time-cpu_start_time)/60.0 ," min"
    write(*,"(a,f12.1,a)") "Speed = ",(1e-3*(time_end-time_init))/((cpu_end_time-cpu_start_time)/3600.0), " kiloyears / hr"

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

        real(prec) :: uxy_rmse, H_rmse, zsrf_rmse, loguxy_rmse 
        real(prec), allocatable :: tmp(:,:) 
        real(prec), allocatable :: tmp1(:,:) 
        
        allocate(tmp(ylmo%grd%nx,ylmo%grd%ny))
        allocate(tmp1(ylmo%grd%nx,ylmo%grd%ny))

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_size(filename,"time",ncid)
        call nc_read(filename,"time",time_prev,start=[n],count=[1],ncid=ncid) 
        if (abs(time-time_prev).gt.1e-5) n = n+1 

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        ! Write model speed 
        call nc_write(filename,"speed",ylmo%par%model_speed,units="kyr/hr",long_name="Model speed (Yelmo only)", &
                      dim1="time",start=[n],count=[1],ncid=ncid)
        
        ! initmip specific error metrics 
        tmp = ylmo%tpo%now%H_ice-ylmo%dta%pd%H_ice
        if (n .gt. 1 .or. count(tmp .ne. 0.0) .gt. 0) then 
            H_rmse = sqrt(sum(tmp**2)/count(tmp .ne. 0.0))
        else 
            H_rmse = mv 
        end if 

        ! surface elevation too 
        tmp = ylmo%dta%pd%err_z_srf
        if (n .gt. 1 .or. count(tmp .ne. 0.0) .gt. 0) then 
            zsrf_rmse = sqrt(sum(tmp**2)/count(tmp .ne. 0.0))
        else 
            zsrf_rmse = mv 
        end if 

        tmp = ylmo%dta%pd%err_uxy_s
        if (n .gt. 1 .or. count(tmp .ne. 0.0) .gt. 0) then 
            uxy_rmse = sqrt(sum(tmp**2)/count(tmp .ne. 0.0))
        else
            uxy_rmse = mv
        end if 

        tmp = ylmo%dta%pd%uxy_s 
        where(ylmo%dta%pd%uxy_s .gt. 0.0) tmp = log(tmp)
        tmp1 = ylmo%dyn%now%uxy_s 
        where(ylmo%dyn%now%uxy_s .gt. 0.0) tmp1 = log(tmp1)
        
        if (n .gt. 1 .or. count(tmp1-tmp .ne. 0.0) .gt. 0) then 
            loguxy_rmse = sqrt(sum((tmp1-tmp)**2)/count(tmp1-tmp .ne. 0.0))
        else
            loguxy_rmse = mv
        end if 
        
        call nc_write(filename,"rmse_H",H_rmse,units="m",long_name="RMSE - Ice thickness", &
                      dim1="time",start=[n],count=[1],ncid=ncid)
        call nc_write(filename,"rmse_zsrf",zsrf_rmse,units="m",long_name="RMSE - Surface elevation", &
                      dim1="time",start=[n],count=[1],ncid=ncid)
        call nc_write(filename,"rmse_uxy",uxy_rmse,units="m/a",long_name="RMSE - Surface velocity", &
                      dim1="time",start=[n],count=[1],ncid=ncid)
        call nc_write(filename,"rmse_uxy_log",loguxy_rmse,units="log(m/a)",long_name="RMSE - Log surface velocity", &
                      dim1="time",start=[n],count=[1],ncid=ncid)
        
        ! == ISMIP6 specific variables 
        ! http://www.climate-cryosphere.org/wiki/index.php?title=InitMIP-Greenland#Appendix_2_.E2.80.93_Naming_conventions.2C_upload_and_model_output_data.

        call nc_write(filename,"H_ice",ylmo%tpo%now%H_ice,units="m",long_name="Ice thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"H_margin",ylmo%tpo%now%H_margin,units="m",long_name="Margin ice thickness", &
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
        call nc_write(filename,"uz_s",ylmo%dyn%now%uz(:,:,ylmo%par%nz_ac),units="m/a",long_name="Surface velocity (z)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"ux_b",ylmo%dyn%now%ux_b,units="m/a",long_name="Basal velocity (x)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy_b",ylmo%dyn%now%uy_b,units="m/a",long_name="Basal velocity (y)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uz_b",ylmo%dyn%now%uz(:,:,1),units="m/a",long_name="Basal velocity (z)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"ux_bar",ylmo%dyn%now%ux_bar,units="m/a",long_name="Vertically averaged velocity (x)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy_bar",ylmo%dyn%now%uy_bar,units="m/a",long_name="Vertically averaged velocity (y)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"T_ice_s",ylmo%thrm%now%T_ice(:,:,ylmo%par%nz_aa),units="K",long_name="Surface ice temperature", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"T_ice_b",ylmo%thrm%now%T_ice(:,:,1),units="K",long_name="Basal ice temperature", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"taub",ylmo%dyn%now%taub,units="Pa",long_name="Basal dragging stress (magnitude)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"calv",ylmo%tpo%now%calv,units="m/a ice equiv.",long_name="Calving rate", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_ice",ylmo%tpo%now%f_ice,units="1",long_name="Total ice fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_ice_g",ylmo%tpo%now%f_ice*ylmo%tpo%now%f_grnd,units="1",long_name="Grounded ice fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_ice_f",ylmo%tpo%now%f_ice*(1.0-ylmo%tpo%now%f_grnd),units="1",long_name="Floating ice fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! == Additional variables 

        call nc_write(filename,"mask_bed",ylmo%tpo%now%mask_bed,units="",long_name="Bed mask", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_grnd",ylmo%tpo%now%f_grnd,units="1",long_name="Grounded ice fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"uxy_s",ylmo%dyn%now%uxy_s,units="m/a",long_name="Surface velocity (magnitude)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_b",ylmo%dyn%now%uxy_b,units="m/a",long_name="Basal velocity (magnitude)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_bar",ylmo%dyn%now%uxy_bar,units="m/a",long_name="Vertically averaged velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)  
        
        call nc_write(filename,"mb_applied",ylmo%tpo%now%mb_applied,units="m",long_name="Applied net mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"N_eff",ylmo%dyn%now%N_eff,units="Pa",long_name="Effective pressure", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"C_bed",ylmo%dyn%now%C_bed,units="Pa",long_name="Bed friction coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"beta",ylmo%dyn%now%beta,units="Pa a m^-1",long_name="Basal friction coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"visc_eff",ylmo%dyn%now%visc_eff,units="Pa a m",long_name="Effective viscosity (SSA)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"T_prime",ylmo%thrm%now%T_ice-ylmo%thrm%now%T_pmp,units="deg C",long_name="Homologous ice temperature", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"f_pmp",ylmo%thrm%now%f_pmp,units="1",long_name="Fraction of grid point at pmp", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"H_w",ylmo%bnd%H_w,units="m water equiv.",long_name="Basal water layer thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! Ice thickness comparison with present-day 
        call nc_write(filename,"pd_err_z_srf",ylmo%dta%pd%err_z_srf,units="m",long_name="Surface elevation error (present day)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"pd_err_uxy_s",ylmo%dta%pd%err_uxy_s,units="m",long_name="Surface velocity error (present day)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! Diagnostics 
        call nc_write(filename,"taud_acx",ylmo%dyn%now%taud_acx,units="Pa",long_name="Driving stress (x)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taud_acy",ylmo%dyn%now%taud_acy,units="Pa",long_name="Driving stress (y)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        

        ! Write cf_ref if it's the first timestep
        if (n .eq. 1) then  
            call nc_write(filename,"cf_ref",cf_ref,units="",long_name="Dragging constant coefficient", &
                          dim1="xc",dim2="yc",start=[1,1],ncid=ncid)
        end if 

        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine write_step_2D

    subroutine calc_ydyn_cbed_external(dyn,tpo,thrm,bnd,grd,domain,mask_noice,cf_ref)
        ! Update C_bed [Pa] based on parameter choices

        implicit none
        
        type(ydyn_class),   intent(INOUT) :: dyn
        type(ytopo_class),  intent(IN)    :: tpo 
        type(ytherm_class), intent(IN)    :: thrm
        type(ybound_class), intent(IN)    :: bnd  
        type(ygrid_class),  intent(IN)    :: grd
        character(len=*),   intent(IN)    :: domain 
        logical,            intent(IN)    :: mask_noice(:,:) 
        real(prec),         intent(IN)    :: cf_ref(:,:) 

        integer :: i, j, nx, ny 
        integer :: i1, i2, j1, j2 
        real(prec) :: f_scale 
        
        real(prec), allocatable :: lambda_bed(:,:)  

        nx = size(dyn%now%C_bed,1)
        ny = size(dyn%now%C_bed,2)
        
        allocate(lambda_bed(nx,ny))

            ! cf_ref [unitless] is obtained as input to this routine from optimization or elsewhere

            ! =============================================================================
            ! Step 2: calculate lambda functions to scale C_bed from default value 
            
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
                        where (bnd%basins .ge. 12 .and. &
                               bnd%basins .le. 17) lambda_bed = calc_lambda_bed_exp(bnd%z_bed,-400.0,dyn%par%cb_z1)


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
            ! Step 3: calculate C_bed [Pa]
            
            dyn%now%C_bed = (cf_ref*lambda_bed)*dyn%now%N_eff 

        return 

    end subroutine calc_ydyn_cbed_external

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

                ! Increase - feeding the Ronne ice shelf from the North
                call scale_cb_gaussian(cf_ref,dyn%par%cf_stream*4.00,x0=-700.0, y0=    0.0,sigma=200.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
                
                ! Increase - Southeast Antarctica inland
                call scale_cb_gaussian(cf_ref,dyn%par%cf_stream*2.00,x0=1500.0, y0= -550.0,sigma=200.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
                call scale_cb_gaussian(cf_ref,dyn%par%cf_stream*2.00,x0=1700.0, y0=-1000.0,sigma=200.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
                
                ! Reduction - South pole 
                call scale_cb_gaussian(cf_ref,dyn%par%cf_stream*0.10,x0=   0.0, y0=   0.0, sigma=400.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
                call scale_cb_gaussian(cf_ref,dyn%par%cf_stream*0.20,x0=   0.0, y0= 600.0, sigma=300.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
                call scale_cb_gaussian(cf_ref,dyn%par%cf_stream*0.10,x0= 500.0, y0=-500.0, sigma=400.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
                
                ! Reduction - Amery ice shelf
                call scale_cb_gaussian(cf_ref,dyn%par%cf_stream*0.05,x0=1500.0, y0= 650.0, sigma=300.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
                
                ! Reduction - feeding the Ross ice shelf from the North
                call scale_cb_gaussian(cf_ref,dyn%par%cf_stream*0.005,x0=-500.0, y0=-500.0, sigma=300.0,xx=grd%x*1e-3,yy=grd%y*1e-3)

                ! Reduction - feeding the Ross ice shelf from the East
                call scale_cb_gaussian(cf_ref,dyn%par%cf_stream*0.05,x0= 130.0, y0=-550.0, sigma=100.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
                call scale_cb_gaussian(cf_ref,dyn%par%cf_stream*0.05,x0= 280.0, y0=-760.0, sigma=100.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
                call scale_cb_gaussian(cf_ref,dyn%par%cf_stream*0.05,x0= 380.0, y0=-960.0, sigma=100.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
                call scale_cb_gaussian(cf_ref,dyn%par%cf_stream*0.05,x0= 400.0, y0=-1150.0,sigma=100.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
                
if (.FALSE.) then
                ! Reduction 
                call scale_cb_gaussian(cf_ref,dyn%par%cf_stream*0.25,x0=-2000.0,y0=1000.0, sigma=100.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
                call scale_cb_gaussian(cf_ref,dyn%par%cf_stream*0.20,x0=-750.0, y0=-900.0, sigma=300.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
                call scale_cb_gaussian(cf_ref,dyn%par%cf_stream*0.20,x0=-600.0, y0=-600.0, sigma=300.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
                call scale_cb_gaussian(cf_ref,dyn%par%cf_stream*0.20,x0=-300.0, y0=   0.0, sigma=100.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
                call scale_cb_gaussian(cf_ref,dyn%par%cf_stream*0.20,x0=-250.0, y0=-500.0, sigma=100.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
                call scale_cb_gaussian(cf_ref,dyn%par%cf_stream*0.20,x0=-100.0, y0=-600.0, sigma=300.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
                call scale_cb_gaussian(cf_ref,dyn%par%cf_stream*0.20,x0=-100.0, y0=-300.0, sigma=300.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
                call scale_cb_gaussian(cf_ref,dyn%par%cf_stream*0.20,x0=   0.0, y0=   0.0, sigma=300.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
                call scale_cb_gaussian(cf_ref,dyn%par%cf_stream*0.20,x0= 130.0, y0=-550.0, sigma= 50.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
                call scale_cb_gaussian(cf_ref,dyn%par%cf_stream*0.20,x0= 280.0, y0=-760.0, sigma= 50.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
                call scale_cb_gaussian(cf_ref,dyn%par%cf_stream*0.20,x0= 380.0, y0=-960.0, sigma= 50.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
                call scale_cb_gaussian(cf_ref,dyn%par%cf_stream*0.50,x0= 400.0, y0=   0.0, sigma=300.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
                call scale_cb_gaussian(cf_ref,dyn%par%cf_stream*0.20,x0= 400.0, y0=-1150.0,sigma= 50.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
                call scale_cb_gaussian(cf_ref,dyn%par%cf_stream*0.25,x0= 700.0, y0= -500.0,sigma=400.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
                call scale_cb_gaussian(cf_ref,dyn%par%cf_stream*0.20,x0=1500.0, y0= 650.0, sigma=300.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
                
                ! Increase
                call scale_cb_gaussian(cf_ref,dyn%par%cf_stream*4.00,x0=-600.0, y0=    0.0,sigma=200.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
                call scale_cb_gaussian(cf_ref,dyn%par%cf_stream*2.00,x0=1200.0, y0=-1200.0,sigma=300.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
                call scale_cb_gaussian(cf_ref,dyn%par%cf_stream*1.50,x0=2000.0, y0=    0.0,sigma=300.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
end if 

            end if 

        return 

    end subroutine modify_cf_ref

    subroutine scale_cb_gaussian(C_bed,cb_new,x0,y0,sigma,xx,yy)

        implicit none 

        real(prec), intent(INOUT) :: C_bed(:,:)
        real(prec), intent(IN) :: cb_new
        real(prec), intent(IN) :: x0
        real(prec), intent(IN) :: y0
        real(prec), intent(IN) :: sigma
        real(prec), intent(IN) :: xx(:,:)
        real(prec), intent(IN) :: yy(:,:)

        ! Local variables 
        integer :: nx, ny 
        real(prec), allocatable :: wts(:,:)
        
        nx = size(C_bed,1)
        ny = size(C_bed,2)

        allocate(wts(nx,ny))

        ! Get Gaussian weights 
        wts = 1.0/(2.0*pi*sigma**2)*exp(-((xx-x0)**2+(yy-y0)**2)/(2.0*sigma**2))
        wts = wts / maxval(wts)

        ! Scale C_bed
        C_bed = C_bed*(1.0-wts) + cb_new*wts

        return 

    end subroutine scale_cb_gaussian

end program yelmo_test 



