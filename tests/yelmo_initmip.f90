

program yelmo_test

    use ncio 
    use yelmo 
    
    ! External libraries
    use sealevel 
    use isostasy  
    
    use snapclim
    use marine_shelf 
    use smbpal   
    use basal_hydrology 
    use sediments 
    use geothermal
    
    implicit none 

    type(yelmo_class)      :: yelmo1
    type(ygrid_class)      :: grd1 

    type(sealevel_class)   :: sealev 
    type(snapclim_class)   :: snp1 
    type(marshelf_class)   :: mshlf1 
    type(smbpal_class)     :: smbpal1 
    type(hydro_class)      :: hyd1  
    type(sediments_class)  :: sed1 
    type(geothermal_class) :: gthrm1
    type(isos_class)       :: isos1
    
    character(len=256) :: outfldr, file1D, file2D, file_restart, domain 
    character(len=512) :: path_par, path_const  
    real(prec) :: time_init, time_end, time_equil, time, dtt, dt1D_out, dt2D_out   
    integer    :: n
    logical    :: calc_transient_climate
    real(4) :: cpu_start_time, cpu_end_time 

    ! No-ice mask (to impose additional melting)
    logical, allocatable :: mask_noice(:,:)  

    ! Concavity field 
    real(prec), allocatable :: channels(:,:) 

    ! Start timing 
    call cpu_time(cpu_start_time)

    ! Determine the parameter file from the command line 
    call yelmo_load_command_line_args(path_par)

    ! Timing and other parameters 
    call nml_read(path_par,"control","time_init",    time_init)                 ! [yr] Starting time
    call nml_read(path_par,"control","time_end",     time_end)                  ! [yr] Ending time
    call nml_read(path_par,"control","time_equil",   time_equil)                ! [yr] Years to equilibrate first
    call nml_read(path_par,"control","dtt",          dtt)                       ! [yr] Main loop time step 
    call nml_read(path_par,"control","dt1D_out",     dt1D_out)                  ! [yr] Frequency of 1D output 
    call nml_read(path_par,"control","dt2D_out",     dt2D_out)                  ! [yr] Frequency of 2D output 
    call nml_read(path_par,"control","transient",    calc_transient_climate)    ! Calculate transient climate? 

    ! Assume program is running from the output folder
    outfldr = "./"

    ! Define input and output locations 
    path_const = trim(outfldr)//"yelmo_const_Earth.nml"
    file1D     = trim(outfldr)//"yelmo1D.nc"
    file2D     = trim(outfldr)//"yelmo2D.nc"

    ! === Initialize ice sheet model =====

    ! General initialization of yelmo constants (used globally)
    call yelmo_global_init(path_const)

    ! Initialize data objects
    call yelmo_init(yelmo1,filename=path_par,grid_def="file")

    ! Initialize state variables (topo only)
    call yelmo_init_state_1(yelmo1,path_par,time=time_init)
    
    ! === Initialize external models (forcing for ice sheet) ======

    ! Store yelmo grid in grd1 and domain name as a shortcut 
    grd1  = yelmo1%grd 
    domain = yelmo1%par%domain 

    ! Initialize global sea level model (bnd%z_sl)
    call sealevel_init(sealev,path_par)

    ! Initialize bedrock model 
    call isos_init(isos1,path_par,grd1%nx,grd1%ny,grd1%dx)

    ! Initialize "climate" model (climate and ocean forcing)
    call snapclim_init(snp1,path_par,domain,grd1%name,grd1%nx,grd1%ny)

    ! Initialize surface mass balance model (bnd%smb, bnd%T_srf)
    call smbpal_init(smbpal1,path_par,x=grd1%xc,y=grd1%yc,lats=grd1%lat)
    
    ! Initialize marine melt model (bnd%bmb_shlf) - impose high basal melt rates
    call marshelf_init(mshlf1,path_par,grd1%nx,grd1%ny,domain,grd1%name,yelmo1%bnd%basins)
    
    ! Initialize basal hydrology (bnd%H_w)
    call hydro_init(hyd1,path_par,grd1%nx,grd1%ny)

    ! Load other constant boundary variables (bnd%H_sed, bnd%Q_geo)
    call sediments_init(sed1,path_par,grd1%nx,grd1%ny,domain,grd1%name)
    call geothermal_init(gthrm1,path_par,grd1%nx,grd1%ny,domain,grd1%name)

    ! === Update initial boundary conditions for current time and yelmo state =====
    ! ybound: z_bed, z_sl, H_sed, H_w, smb, T_srf, bmb_shlf , Q_geo

    call isos_init_state(isos1,z_bed=yelmo1%bnd%z_bed,z_bed_ref=yelmo1%bnd%z_bed, &
                               H_ice_ref=yelmo1%tpo%now%H_ice,z_sl=yelmo1%bnd%z_sl,time=time_init)

    call sealevel_update(sealev,year_bp=time_init)
    yelmo1%bnd%z_sl  = sealev%z_sl 

    yelmo1%bnd%H_sed = sed1%now%H 
    
    yelmo1%bnd%H_w   = 0.0   ! use hydro_init_state later

!     ! Normal snapclim call 
!     call snapclim_update(snp1,z_srf=yelmo1%tpo%now%z_srf,year_bp=time_init,domain=domain) 

!     ! Equilibrate snowpack for itm
!     if (trim(smbpal1%par%abl_method) .eq. "itm") then 
!         call smbpal_update_monthly_equil(smbpal1,snp1%now%tas,snp1%now%pr, &
!                                yelmo1%tpo%now%z_srf,yelmo1%tpo%now%H_ice,time_init,time_equil=100.0)
!     end if 

!     call smbpal_update_monthly(smbpal1,snp1%now%tas,snp1%now%pr, &
!                                yelmo1%tpo%now%z_srf,yelmo1%tpo%now%H_ice,time_init) 
!     yelmo1%bnd%smb   = smbpal1%ann%smb*conv_we_ie*1e-3    ! [mm we/a] => [m ie/a]
!     yelmo1%bnd%T_srf = smbpal1%ann%tsrf 

    ! Impose present-day surface mass balance and present-day temperature field
    yelmo1%bnd%smb   = yelmo1%dta%pd%smb_ann
    yelmo1%bnd%T_srf = yelmo1%dta%pd%t2m_ann
    
    call marshelf_calc_Tshlf(mshlf1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                         yelmo1%bnd%z_sl,depth=snp1%now%depth,to_ann=snp1%now%to_ann, &
                         dto_ann=snp1%now%to_ann - snp1%clim0%to_ann)

    call marshelf_update(mshlf1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                         yelmo1%bnd%z_sl,regions=yelmo1%bnd%regions,dx=yelmo1%grd%dx)

    yelmo1%bnd%bmb_shlf = mshlf1%now%bmb_shlf  
    yelmo1%bnd%T_shlf   = mshlf1%now%T_shlf  

    yelmo1%bnd%Q_geo    = gthrm1%now%ghf 
    
    call yelmo_print_bound(yelmo1%bnd)

    time = time_init 

    ! Define channel field 
    allocate(channels(grd1%nx,grd1%ny))

    ! Update C_bed if needed
    if (yelmo1%dyn%par%C_bed_method .eq. -1) then 
        call calc_ydyn_cbed_external(yelmo1%dyn,yelmo1%tpo,yelmo1%thrm,yelmo1%bnd,channels)
    end if 

    ! Initialize state variables (dyn,therm,mat)
    ! (initialize temps with robin method with a cold base)
    call yelmo_init_state_2(yelmo1,path_par,time=time_init,thrm_method="robin-cold")

    ! Update C_bed if needed
    if (yelmo1%dyn%par%C_bed_method .eq. -1) then 
        call calc_ydyn_cbed_external(yelmo1%dyn,yelmo1%tpo,yelmo1%thrm,yelmo1%bnd,channels)
    end if 
    
    ! Define no-ice mask from present-day data
    allocate(mask_noice(grd1%nx,grd1%ny))
    mask_noice = .FALSE. 
    where(yelmo1%dta%pd%H_ice .le. 0.0) mask_noice = .TRUE. 

    ! Impose additional negative mass balance to no ice points 2 [m/a i.e.] melting
    where(mask_noice) yelmo1%bnd%smb   = yelmo1%dta%pd%smb_ann - 2.0 

    ! Impose a colder boundary temperature for equilibration step 
    ! -5 [K] for mimicking glacial times
    yelmo1%bnd%T_srf = yelmo1%dta%pd%t2m_ann - 10.0  

    ! Run yelmo for several years with constant boundary conditions and topo
    ! to equilibrate thermodynamics and dynamics
    call yelmo_update_equil(yelmo1,time,time_tot=time_equil,topo_fixed=.FALSE.,dt=1.0,ssa_vel_max=500.0)
    
    ! 2D file 
    call yelmo_write_init(yelmo1,file2D,time_init=time,units="years")
!     call write_step_2D(yelmo1,file2D,time=time)  
!     call write_step_2D_small(yelmo1,file2D,time=time)  
    call write_step_2D_combined(yelmo1,isos1,snp1,mshlf1,smbpal1,hyd1,file2D,time=time)
    
    ! 1D file 
    call write_yreg_init(yelmo1,file1D,time_init=time_init,units="years",mask=yelmo1%bnd%ice_allowed)
    call write_yreg_step(yelmo1%reg,file1D,time=time) 
    
    ! Advance timesteps
    do n = 1, ceiling((time_end-time_init)/dtt)

        ! Get current time 
        time = time_init + n*dtt

        ! == SEA LEVEL ==========================================================
        call sealevel_update(sealev,year_bp=time)
        yelmo1%bnd%z_sl  = sealev%z_sl 

        ! Update C_bed if needed
        if (yelmo1%dyn%par%C_bed_method .eq. -1) then 
            call calc_ydyn_cbed_external(yelmo1%dyn,yelmo1%tpo,yelmo1%thrm,yelmo1%bnd,channels)
        end if 
    
        ! == Yelmo ice sheet ===================================================
        call yelmo_update(yelmo1,time)

        ! == ISOSTASY ==========================================================
        call isos_update(isos1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_sl,time)
        yelmo1%bnd%z_bed = isos1%now%z_bed

if (calc_transient_climate) then 
        ! == CLIMATE (ATMOSPHERE AND OCEAN) ====================================
!         if (mod(time,10.0)==0) then
!             ! Normal snapclim call 
!             call snapclim_update(snp1,z_srf=yelmo1%tpo%now%z_srf,year_bp=time,domain=domain)
!         end if 

!         ! == SURFACE MASS BALANCE ==============================================
!         call smbpal_update_monthly(smbpal1,snp1%now%tas,snp1%now%pr, &
!                                    yelmo1%tpo%now%z_srf,yelmo1%tpo%now%H_ice,time) 
!         yelmo1%bnd%smb   = smbpal1%ann%smb*conv_we_ie*1e-3       ! [mm we/a] => [m ie/a]
!         yelmo1%bnd%T_srf = smbpal1%ann%tsrf 

!         yelmo1%bnd%smb   = yelmo1%dta%pd%smb_ann
!         yelmo1%bnd%T_srf = yelmo1%dta%pd%t2m_ann
    
        ! == MARINE AND TOTAL BASAL MASS BALANCE ===============================
        call marshelf_calc_Tshlf(mshlf1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                         yelmo1%bnd%z_sl,depth=snp1%now%depth,to_ann=snp1%now%to_ann, &
                         dto_ann=snp1%now%to_ann - snp1%clim0%to_ann)

        call marshelf_update(mshlf1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                         yelmo1%bnd%z_sl,regions=yelmo1%bnd%regions,dx=yelmo1%grd%dx)

end if 

        yelmo1%bnd%bmb_shlf = mshlf1%now%bmb_shlf  
        yelmo1%bnd%T_shlf   = mshlf1%now%T_shlf  



        ! Update temperature and smb as needed (ISMIP6)
        if (time .ge. -10e6 .and. time .lt. -10e3) then 
            ! Glacial period, impose cold climate 
            yelmo1%bnd%T_srf = yelmo1%dta%pd%t2m_ann - 10.0 

        else if (time .ge. -10e3 .and. time .lt. -8e3) then
            ! Holocene optimum 
            yelmo1%bnd%T_srf = yelmo1%dta%pd%t2m_ann + 1.0 

        else if (time .ge. -8e3) then 
            ! Entering Holocene, impose present-day temperatures 
            yelmo1%bnd%T_srf = yelmo1%dta%pd%t2m_ann
        end if 


        ! == BASAL HYDROLOGY ====================================================
        call hydro_update(hyd1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%z_srf,yelmo1%bnd%z_sl, &
                          yelmo1%tpo%now%bmb,yelmo1%tpo%now%f_grnd,yelmo1%tpo%par%dx,yelmo1%tpo%par%dx,time)
        yelmo1%bnd%H_w   = hyd1%now%H_water 

        ! == MODEL OUTPUT =======================================================

!         if (time == 30.0) then
!             file_restart = trim(outfldr)//"yelmo_restart_yr30.0.nc" 
!             call yelmo_restart_write(yelmo1,file_restart,time=time)
!         end if 
        
        if (mod(time,dt2D_out)==0) then 
!             call write_step_2D(yelmo1,file2D,time=time)
!             call write_step_2D_small(yelmo1,file2D,time=time)
            call write_step_2D_combined(yelmo1,isos1,snp1,mshlf1,smbpal1,hyd1,file2D,time=time)
        end if 

        if (mod(time,dt1D_out)==0) then 
            call write_yreg_step(yelmo1%reg,file1D,time=time) 
        end if 

        if (mod(time,10.0)==0) then
            write(*,"(a,f14.4)") "yelmo::       time = ", time
        end if 

    end do 

!     ! Let's see if we can read a restart file 
!     call yelmo_restart_read(yelmo1,file_restart,time=time)

    ! Finalize program
    call yelmo_end(yelmo1,time=time)

    ! Stop timing 
    call cpu_time(cpu_end_time)

    write(*,"(a,f12.3,a)") "Time  = ",(cpu_end_time-cpu_start_time)/60.0 ," min"
    write(*,"(a,f12.1,a)") "Speed = ",(1e-3*(time_end-time_init))/((cpu_end_time-cpu_start_time)/3600.0), "kiloyears / hr"

contains

    subroutine write_step_2D_small(ylmo,filename,time)

        implicit none 
        
        type(yelmo_class),      intent(IN) :: ylmo
!         type(snapclim_class),   intent(IN) :: snp 
!         type(marshelf_class),   intent(IN) :: mshlf 
!         type(smbpal_class),     intent(IN) :: srf 
!         type(hydro_class),      intent(IN) :: hyd  
        !type(sediments_class),  intent(IN) :: sed 
        !type(geothermal_class), intent(IN) :: gthrm
        !type(isos_class),       intent(IN) :: isos
        
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

        ! == yelmo_topography ==
        call nc_write(filename,"H_ice",ylmo%tpo%now%H_ice,units="m",long_name="Ice thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf",ylmo%tpo%now%z_srf,units="m",long_name="Surface elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_bed",ylmo%tpo%now%mask_bed,units="",long_name="Bed mask", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"beta",ylmo%dyn%now%beta,units="Pa a m-1",long_name="Dragging coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"visc_eff",ylmo%dyn%now%visc_eff,units="Pa a^-1",long_name="Effective viscosity (SSA)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"uxy_i_bar",ylmo%dyn%now%uxy_i_bar,units="m/a",long_name="Internal shear velocity magnitude", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_b",ylmo%dyn%now%uxy_b,units="m/a",long_name="Basal sliding velocity magnitude", &
                     dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_bar",ylmo%dyn%now%uxy_bar,units="m/a",long_name="Vertically integrated velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
!         call nc_write(filename,"z_sl",ylmo%bnd%z_sl,units="m",long_name="Sea level rel. to present", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"z_bed",ylmo%bnd%z_bed,units="m",long_name="Bedrock elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine write_step_2D_small

    subroutine write_step_2D_combined(ylmo,isos,snp,mshlf,srf,hyd,filename,time)

        implicit none 
        
        type(yelmo_class),      intent(IN) :: ylmo
        type(isos_class),       intent(IN) :: isos 
        type(snapclim_class),   intent(IN) :: snp 
        type(marshelf_class),   intent(IN) :: mshlf 
        type(smbpal_class),     intent(IN) :: srf 
        type(hydro_class),      intent(IN) :: hyd  
        !type(sediments_class),  intent(IN) :: sed 
        !type(geothermal_class), intent(IN) :: gthrm
        !type(isos_class),       intent(IN) :: isos
        
        character(len=*),  intent(IN) :: filename
        real(prec), intent(IN) :: time

        ! Local variables
        integer    :: ncid, n
        real(prec) :: time_prev 

        real(prec) :: uxy_rmse, H_rmse 
        real(prec), allocatable :: tmp(:,:) 

        allocate(tmp(ylmo%grd%nx,ylmo%grd%ny))

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
        H_rmse = sqrt(sum(tmp**2)/count(tmp .ne. 0.0))

        tmp = ylmo%dyn%now%uxy_s-ylmo%dta%pd%uxy_s
        uxy_rmse = sqrt(sum(tmp**2)/count(tmp .ne. 0.0))

        call nc_write(filename,"rmse_H",H_rmse,units="m",long_name="RMSE - Ice thickness", &
                      dim1="time",start=[n],count=[1],ncid=ncid)
        call nc_write(filename,"rmse_uxy",uxy_rmse,units="m/a",long_name="RMSE - Surface velocity", &
                      dim1="time",start=[n],count=[1],ncid=ncid)
        
        ! == yelmo_topography ==
        call nc_write(filename,"H_ice",ylmo%tpo%now%H_ice,units="m",long_name="Ice thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf",ylmo%tpo%now%z_srf,units="m",long_name="Surface elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_bed",ylmo%tpo%now%mask_bed,units="",long_name="Bed mask", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"mb_applied",ylmo%tpo%now%mb_applied,units="m",long_name="Applied net mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"H_grnd",ylmo%tpo%now%H_grnd,units="m",long_name="Ice thickness overburden", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"N_eff",ylmo%tpo%now%N_eff,units="bar",long_name="Effective pressure", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"calv",ylmo%tpo%now%calv,units="m/a ice equiv.",long_name="Calving rate", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"f_grnd",ylmo%tpo%now%f_grnd,units="1",long_name="Grounded fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_ice",ylmo%tpo%now%f_ice,units="1",long_name="Ice fraction in grid cell", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"dist_grline",ylmo%tpo%now%dist_grline,units="km", &
                      long_name="Distance to nearest grounding-line point", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"C_bed",ylmo%dyn%now%C_bed,units="m a^-1 Pa^-2",long_name="Dragging constant", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"beta",ylmo%dyn%now%beta,units="Pa a m-1",long_name="Dragging coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"visc_eff",ylmo%dyn%now%visc_eff,units="Pa a^-1",long_name="Vertically averaged viscosity", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"uxy_i_bar",ylmo%dyn%now%uxy_i_bar,units="m/a",long_name="Internal shear velocity magnitude", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_b",ylmo%dyn%now%uxy_b,units="m/a",long_name="Basal sliding velocity magnitude", &
                     dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_bar",ylmo%dyn%now%uxy_bar,units="m/a",long_name="Vertically integrated velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"T_ice",ylmo%thrm%now%T_ice,units="K",long_name="Ice temperature", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        
        call nc_write(filename,"T_prime",ylmo%thrm%now%T_ice-ylmo%thrm%now%T_pmp,units="deg C",long_name="Homologous ice temperature", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"f_pmp",ylmo%thrm%now%f_pmp,units="1",long_name="Fraction of grid point at pmp", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"Q_strn",ylmo%thrm%now%Q_strn/(rho_ice*ylmo%thrm%now%cp),units="K a-1",long_name="Strain heating", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"Q_b",ylmo%thrm%now%Q_b,units="J a-1 m-2",long_name="Basal frictional heating", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"bmb_grnd",ylmo%thrm%now%bmb_grnd,units="m/a ice equiv.",long_name="Basal mass balance (grounded)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"ATT",ylmo%mat%now%ATT,units="a^-1 Pa^-3",long_name="Rate factor", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        call nc_write(filename,"enh_bar",ylmo%mat%now%enh_bar,units="1",long_name="Vertically averaged enhancement factor", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"visc_bar",ylmo%mat%now%visc_bar,units="Pa a^-1",long_name="Vertically averaged viscosity", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! Boundaries
        call nc_write(filename,"z_bed",ylmo%bnd%z_bed,units="m",long_name="Bedrock elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid) 
        call nc_write(filename,"z_sl",ylmo%bnd%z_sl,units="m",long_name="Sea level rel. to present", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"H_sed",ylmo%bnd%H_sed,units="m",long_name="Sediment thickness", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"H_w",ylmo%bnd%H_w,units="m",long_name="Basal water layer thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"smb",ylmo%bnd%smb,units="m/a ice equiv.",long_name="Surface mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"smb_errpd",ylmo%bnd%smb-ylmo%dta%pd%smb_ann,units="m/a ice equiv.",long_name="Surface mass balance error wrt present day", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"T_srf",ylmo%bnd%T_srf,units="K",long_name="Surface temperature", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"bmb_shlf",ylmo%bnd%bmb_shlf,units="m/a ice equiv.",long_name="Basal mass balance (shelf)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"Q_geo",ylmo%bnd%Q_geo,units="mW/m^2",long_name="Geothermal heat flux", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_sl",ylmo%bnd%z_sl,units="m",long_name="Sea level rel. to present", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"bmb",ylmo%tpo%now%bmb,units="m/a ice equiv.",long_name="Basal mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! External data
        call nc_write(filename,"dzbdt",isos%now%dzbdt,units="m/a",long_name="Bedrock uplift rate", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
 
        call nc_write(filename,"dT_shlf",mshlf%now%dT_shlf,units="K",long_name="Shelf temperature anomaly", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"Ta_ann",snp%now%ta_ann,units="K",long_name="Near-surface air temperature (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"Ta_sum",snp%now%ta_sum,units="K",long_name="Near-surface air temperature (sum)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"Pr_ann",snp%now%pr_ann*1e-3,units="m/a water equiv.",long_name="Precipitation (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"dTa_ann",snp%now%ta_ann-snp%clim0%ta_ann,units="K",long_name="Near-surface air temperature anomaly (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dTa_sum",snp%now%ta_sum-snp%clim0%ta_sum,units="K",long_name="Near-surface air temperature anomaly (sum)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dPr_ann",(snp%now%pr_ann-snp%clim0%pr_ann)*1e-3,units="m/a water equiv.",long_name="Precipitation anomaly (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
!         call nc_write(filename,"PDDs",srf%ann%PDDs,units="degC days",long_name="Positive degree days (annual total)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! Ice thickness comparison with present-day 
        call nc_write(filename,"H_ice_errpd",ylmo%tpo%now%H_ice-ylmo%dta%pd%H_ice,units="m",long_name="Ice thickness error wrt present day", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_s_errpd",ylmo%dyn%now%uxy_s-ylmo%dta%pd%uxy_s,units="m",long_name="Ice thickness error wrt present day", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"ssa_mask_acx",ylmo%dyn%now%ssa_mask_acx,units="1",long_name="SSA mask (acx)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"ssa_mask_acy",ylmo%dyn%now%ssa_mask_acy,units="1",long_name="SSA mask (acy)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"dzsdx",ylmo%tpo%now%dzsdx,units="m/m",long_name="Surface slope", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dzsdy",ylmo%tpo%now%dzsdy,units="m/m",long_name="Surface slope", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"f_grnd_acx",ylmo%tpo%now%f_grnd_acx,units="1",long_name="Grounded fraction (acx)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_grnd_acy",ylmo%tpo%now%f_grnd_acy,units="1",long_name="Grounded fraction (acy)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"taub_acx",ylmo%dyn%now%taub_acx,units="Pa",long_name="Basal stress (x)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taub_acy",ylmo%dyn%now%taub_acy,units="Pa",long_name="Basal stress (y)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"taud_acx",ylmo%dyn%now%taud_acx,units="Pa",long_name="Driving stress (x)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taud_acy",ylmo%dyn%now%taud_acy,units="Pa",long_name="Driving stress (y)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"ux_b",ylmo%dyn%now%ux_b,units="m/a",long_name="Basal sliding velocity (x)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy_b",ylmo%dyn%now%uy_b,units="m/a",long_name="Basal sliding velocity (y)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"ux_i_bar",ylmo%dyn%now%ux_i_bar,units="m/a",long_name="Internal shear velocity (x)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy_i_bar",ylmo%dyn%now%uy_i_bar,units="m/a",long_name="Internal shear velocity (y)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"ux_bar",ylmo%dyn%now%ux_bar,units="m/a",long_name="Vertically averaged velocity (x)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy_bar",ylmo%dyn%now%uy_bar,units="m/a",long_name="Vertically averaged velocity (y)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"beta_acx",ylmo%dyn%now%beta_acx,units="Pa a m-1",long_name="Dragging coefficient (x)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"beta_acy",ylmo%dyn%now%beta_acy,units="Pa a m-1",long_name="Dragging coefficient (y)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! Channels 
        call nc_write(filename,"channels",channels,units="m m-1",long_name="Channel diagnosis", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine write_step_2D_combined

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

        ! == yelmo_topography ==
        call nc_write(filename,"H_ice",ylmo%tpo%now%H_ice,units="m",long_name="Ice thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf",ylmo%tpo%now%z_srf,units="m",long_name="Surface elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_bed",ylmo%tpo%now%mask_bed,units="",long_name="Bed mask", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"dzsdx",ylmo%tpo%now%dzsdx,units="m/m",long_name="Surface slope", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dzsdy",ylmo%tpo%now%dzsdy,units="m/m",long_name="Surface slope", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"dzsrfdt",ylmo%tpo%now%dzsrfdt,units="m/a",long_name="Surface elevation change", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dHicedt",ylmo%tpo%now%dHicedt,units="m/a",long_name="Ice thickness change", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"f_grnd",ylmo%tpo%now%f_grnd,units="1",long_name="Grounded fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_ice",ylmo%tpo%now%f_ice,units="1",long_name="Ice-covered fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"N_eff",ylmo%tpo%now%N_eff,units="bar",long_name="Effective pressure", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"calv",ylmo%tpo%now%calv,units="m/a",long_name="Calving rate", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! == yelmo_thermodynamics ==
        call nc_write(filename,"T_ice",ylmo%thrm%now%T_ice,units="K",long_name="Ice temperature", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"T_prime",ylmo%thrm%now%T_ice-ylmo%thrm%now%T_pmp,units="deg C",long_name="Homologous ice temperature", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"f_pmp",ylmo%thrm%now%f_pmp,units="1",long_name="Fraction of grid point at pmp", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"Q_strn",ylmo%thrm%now%Q_strn/(rho_ice*ylmo%thrm%now%cp),units="K a-1",long_name="Strain heating", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"Q_b",ylmo%thrm%now%Q_b,units="J a-1 m-2",long_name="Basal frictional heating", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! == yelmo_material ==
        call nc_write(filename,"visc",ylmo%mat%now%visc,units="Pa a-1",long_name="Viscosity", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"visc_bar",ylmo%mat%now%visc_bar,units="Pa a^-1",long_name="Vertically averaged viscosity (vertically integerated)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"ATT",ylmo%mat%now%ATT,units="a^-1 Pa^-3",long_name="Rate factor", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        call nc_write(filename,"enh",ylmo%mat%now%enh,units="1",long_name="Enhancement factor", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"enh_bar",ylmo%mat%now%enh_bar,units="1",long_name="Vertically averaged enhancement factor", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! == yelmo_dynamics ==

        call nc_write(filename,"ssa_mask_acx",ylmo%dyn%now%ssa_mask_acx,units="1",long_name="SSA mask (acx)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"ssa_mask_acy",ylmo%dyn%now%ssa_mask_acy,units="1",long_name="SSA mask (acy)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"C_bed",ylmo%dyn%now%C_bed,units="1e5 a m^-1",long_name="Dragging constant", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"beta",ylmo%dyn%now%beta,units="Pa a m^-1",long_name="Dragging coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"visc_eff_dyn",ylmo%dyn%now%visc_eff,units="1",long_name="Effective viscosity from dyn module", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"sigma_horiz_sq",ylmo%dyn%now%sigma_horiz_sq,units="1",long_name="Horizontal stress components squared", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"lhs_x",ylmo%dyn%now%lhs_x,units="Pa",long_name="Shear reduction (x)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"lhs_y",ylmo%dyn%now%lhs_y,units="Pa",long_name="Shear reduction (y)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"lhs_xy",ylmo%dyn%now%lhs_xy,units="Pa",long_name="Shear reduction magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"duxdz",ylmo%dyn%now%duxdz,units="1/a",long_name="Vertical shear (x)", &
                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"duydz",ylmo%dyn%now%duydz,units="1/a",long_name="Vertical shear (y)", &
                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        
        call nc_write(filename,"taud_acx",ylmo%dyn%now%taud_acx,units="Pa",long_name="Driving stress (x)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taud_acy",ylmo%dyn%now%taud_acy,units="Pa",long_name="Driving stress (y)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
!         call nc_write(filename,"ux_i_bar",ylmo%dyn%now%ux_i_bar,units="m/a",long_name="Internal shear velocity (x)", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"uy_i_bar",ylmo%dyn%now%uy_i_bar,units="m/a",long_name="Internal shear velocity (y)", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_i_bar",ylmo%dyn%now%uxy_i_bar,units="m/a",long_name="Internal shear velocity magnitude", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"dd_ab",ylmo%dyn%now%dd_ab,units="m2 a-1",long_name="Diffusivity", &
                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"dd_ab_bar",ylmo%dyn%now%dd_ab_bar,units="m2 a-1",long_name="Diffusivity", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"ux_b",ylmo%dyn%now%ux_b,units="m/a",long_name="Basal sliding velocity (x)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"uy_b",ylmo%dyn%now%uy_b,units="m/a",long_name="Basal sliding velocity (y)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_b",ylmo%dyn%now%uxy_b,units="m/a",long_name="Basal sliding velocity magnitude", &
                     dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
!         call nc_write(filename,"ux_bar",ylmo%dyn%now%ux_bar,units="m/a",long_name="Vertically integrated velocity (x)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"uy_bar",ylmo%dyn%now%uy_bar,units="m/a",long_name="Vertically integrated velocity (y)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_bar",ylmo%dyn%now%uxy_bar,units="m/a",long_name="Vertically integrated velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"ux",ylmo%dyn%now%ux,units="m/a",long_name="Horizontal velocity (x)", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"uy",ylmo%dyn%now%uy,units="m/a",long_name="Horizontal velocity (y)", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"uxy",ylmo%dyn%now%uxy,units="m/a",long_name="Horizontal velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"uz",ylmo%dyn%now%uz,units="m/a",long_name="Vertical velocity", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        call nc_write(filename,"f_vbvs",ylmo%dyn%now%f_vbvs,units="1",long_name="Basal to surface velocity fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"f_shear_bar",ylmo%mat%now%f_shear_bar,units="1",long_name="Vertically averaged shearing fraction", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"de",ylmo%mat%now%strn%de,units="",long_name="Vertically averaged strain rate", &
!                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

!         call nc_write(filename,"de_bar",ylmo%mat%now%strn2D%de,units="",long_name="Vertically averaged strain rate", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"tau_mx",ylmo%dyn%now%tau_mx,units="Pa",long_name="Basal friction (x)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"tau_my",ylmo%dyn%now%tau_my,units="Pa",long_name="Basal friction (y)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! == yelmo_bound ==
        call nc_write(filename,"z_bed",ylmo%bnd%z_bed,units="m",long_name="Bedrock elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"z_sl",ylmo%bnd%z_sl,units="m",long_name="Sea level rel. to present", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"H_sed",ylmo%bnd%H_sed,units="m",long_name="Sediment thickness", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"H_w",ylmo%bnd%H_w,units="m",long_name="Basal water pressure", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"smb",ylmo%bnd%smb,units="m/a",long_name="Surface mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"T_srf",ylmo%bnd%T_srf,units="K",long_name="Surface temperature", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"bmb_shlf",ylmo%bnd%bmb_shlf,units="m/a",long_name="Basal mass balance (shelf)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"Q_geo",ylmo%bnd%Q_geo,units="mW/m^2",long_name="Geothermal heat flux", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        

        call nc_write(filename,"bmb",ylmo%tpo%now%bmb,units="m/a",long_name="Basal mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine write_step_2D


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
        f_channel = exp(-channels/channel_lim)
        where(f_channel .lt. 0.1) f_channel = 0.1 
        where(f_channel .gt. 2.0) f_channel = 2.0  

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



