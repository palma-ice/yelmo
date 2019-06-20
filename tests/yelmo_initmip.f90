

program yelmo_test

    use ncio 
    use yelmo 
    
    use basal_hydro_simple 

    implicit none 

    type(yelmo_class)      :: yelmo1
    type(hydro_class)      :: hyd1 

    character(len=256) :: outfldr, file1D, file2D, file_restart, domain 
    character(len=512) :: path_par, path_const  
    real(prec) :: time_init, time_end, time_equil, time, dtt, dt1D_out, dt2D_out   
    integer    :: n
    real(4) :: cpu_start_time, cpu_end_time 

    ! No-ice mask (to impose additional melting)
    logical, allocatable :: mask_noice(:,:)  

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

    yelmo1%bnd%z_sl     = 0.0           ! [m]
    yelmo1%bnd%H_sed    = 0.0           ! [m]
    yelmo1%bnd%H_w      = hyd1%now%H_w  ! [m]
    yelmo1%bnd%Q_geo    = 50.0          ! [mW/m2]
    
    yelmo1%bnd%bmb_shlf = -20.0         ! [m.i.e./a]
    yelmo1%bnd%T_shlf   = T0            ! [K]   

    ! Impose present-day surface mass balance and present-day temperature field
    yelmo1%bnd%smb      = yelmo1%dta%pd%smb        ! [m.i.e./a]
    yelmo1%bnd%T_srf    = yelmo1%dta%pd%T_srf      ! [K]
    
    call yelmo_print_bound(yelmo1%bnd)

    time = time_init 

    ! Initialize state variables (dyn,therm,mat)
    ! (initialize temps with robin method with a cold base)
    call yelmo_init_state(yelmo1,path_par,time=time_init,thrm_method="robin-cold")

    ! Define no-ice mask from present-day data
    allocate(mask_noice(yelmo1%grd%nx,yelmo1%grd%ny))
    mask_noice = .FALSE. 
    where(yelmo1%dta%pd%H_ice .le. 0.0) mask_noice = .TRUE. 

    ! Impose additional negative mass balance to no ice points of 2 [m.i.e./a] melting
    where(mask_noice) yelmo1%bnd%smb = yelmo1%dta%pd%smb - 2.0 

    ! Impose a colder boundary temperature for equilibration step 
    ! -5 [K] for mimicking glacial times
!     yelmo1%bnd%T_srf = yelmo1%dta%pd%T_srf - 10.0  

    ! Run yelmo for several years with constant boundary conditions and topo
    ! to equilibrate thermodynamics and dynamics
    call yelmo_update_equil(yelmo1,time,time_tot=time_equil,topo_fixed=.FALSE.,dt=0.5,ssa_vel_max=500.0)
    
    ! 2D file 
    call yelmo_write_init(yelmo1,file2D,time_init=time,units="years")  
    call write_step_2D(yelmo1,file2D,time=time)
    
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
                    yelmo1%thrm%now%bmb_grnd*rho_ice/rho_w,time)

        ! Pass updated boundary variables to yelmo 
        yelmo1%bnd%H_w = hyd1%now%H_w 

        ! == MODEL OUTPUT =======================================================

        if (mod(nint(time*100),nint(dt2D_out*100))==0) then
            call write_step_2D(yelmo1,file2D,time=time)
        end if 

        if (mod(nint(time*100),nint(dt1D_out*100))==0) then 
            call write_yreg_step(yelmo1%reg,file1D,time=time) 
        end if 

        if (mod(time,10.0)==0) then
            write(*,"(a,f14.4)") "yelmo::       time = ", time
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

    subroutine write_step_2D(ylmo,filename,time)

        implicit none 
        
        type(yelmo_class), intent(IN) :: ylmo
        character(len=*),  intent(IN) :: filename
        real(prec), intent(IN) :: time

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
        
        call nc_write(filename,"uxy_s",ylmo%dyn%now%uxy_s,units="m/a",long_name="Surface velocity (magnitude)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_b",ylmo%dyn%now%uxy_b,units="m/a",long_name="Basal velocity (magnitude)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_bar",ylmo%dyn%now%uxy_bar,units="m/a",long_name="Vertically averaged velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)  
        
        call nc_write(filename,"mb_applied",ylmo%tpo%now%mb_applied,units="m",long_name="Applied net mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"N_eff",ylmo%dyn%now%N_eff*1e-5,units="1e5 Pa (Bar)",long_name="Effective pressure", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"C_bed",ylmo%dyn%now%C_bed,units="",long_name="Dragging constant", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"beta",ylmo%dyn%now%beta,units="Pa a m^-1",long_name="Dragging coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"visc_eff",ylmo%dyn%now%visc_eff,units="Pa a m",long_name="Vertically integrated viscosity", &
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
        

        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine write_step_2D

end program yelmo_test 



