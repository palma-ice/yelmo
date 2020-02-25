

program yelmo_benchmarks

    use nml 
    use ncio  
    use yelmo 

    use deformation 

    use ice_benchmarks 
    use mismip3D 

    implicit none 

    type(yelmo_class)     :: yelmo1

    type(bueler_test_type) :: buel 
    
    character(len=56)  :: domain    
    character(len=256) :: outfldr, file2D, file1D, file_compare
    character(len=256) :: file_restart
    character(len=512) :: path_par, path_const 
    character(len=56)  :: experiment
    logical            :: with_ssa  
    logical    :: topo_fixed, dyn_fixed 
    character(len=256) :: topo_fixed_file 
    real(prec) :: time_init, time_end, time, dtt, dt2D_out, dt1D_out
    real(prec) :: period, dt_test 
    integer    :: n  

    character(len=56) :: grid_name
    real(prec) :: dx 
    integer    :: nx  

    real(4) :: start, finish
    
    ! Start timing 
    call cpu_time(start)

    ! Assume program is running from the output folder
    outfldr = "./"

    ! Determine the parameter file from the command line 
    call yelmo_load_command_line_args(path_par)
    !path_par   = trim(outfldr)//"yelmo_EISMINT.nml" 
    path_const = trim(outfldr)//"yelmo_const_EISMINT.nml"
    
    ! Define input and output locations 
    file1D       = trim(outfldr)//"yelmo1D.nc"
    file2D       = trim(outfldr)//"yelmo2D.nc"
    file_compare = trim(outfldr)//"yelmo_compare.nc"
    file_restart = trim(outfldr)//"yelmo_restart.nc"

    
    ! Define the domain, grid and experiment from parameter file
    call nml_read(path_par,"eismint","domain",       domain)        ! EISMINT1, EISMINT2
    call nml_read(path_par,"eismint","experiment",   experiment)    ! "fixed", "moving", "mismip", "EXPA", "EXPB", "BUELER-A"
    call nml_read(path_par,"eismint","dx",           dx)            ! [km] Grid resolution 
    call nml_read(path_par,"eismint","with_ssa",     with_ssa)      ! Include ssa in experiment?
    call nml_read(path_par,"eismint","topo_fixed",   topo_fixed)    ! Calculate the topography, or use Heiko's topo file. 
    call nml_read(path_par,"eismint","dyn_fixed",    dyn_fixed)     ! Calculate the topography, or use Heiko's topo file. 
    call nml_read(path_par,"eismint","topo_fixed_file",topo_fixed_file)     ! File containing fixed topo field of H_ice
    
    ! Timing parameters 
    call nml_read(path_par,"eismint","time_init",    time_init)     ! [yr] Starting time
    call nml_read(path_par,"eismint","time_end",     time_end)      ! [yr] Ending time
    call nml_read(path_par,"eismint","dtt",          dtt)           ! [yr] Main loop time step 
    call nml_read(path_par,"eismint","dt2D_out",     dt2D_out)      ! [yr] Frequency of 2D output 
    dt1D_out = dtt  ! Set 1D output to frequency of main loop timestep 

    ! Settings for transient EISMINT1 experiments 
    call nml_read(path_par,"eismint","period",       period)        ! [yr] for transient experiments 
    call nml_read(path_par,"eismint","dT_test",      dT_test)       ! [K] for test experiments  
    
    ! Define the model domain based on the experiment we are running
    select case(trim(experiment))

        case("mismip")
            ! EISMINT1 experiment with modified EISMINT1-mismip grid proposed by Heiko

            grid_name = "EISMINT1-mismip"
            nx = 73  

        case("BUELER-A","BUELER-B")

            grid_name = "EISMINT-EXT"
            nx = (2000.0 / dx) + 1      ! Domain width is 2000 km total (-1000 to 1000 km)
            
        case("HALFAR")

            grid_name = "HALFAR"
            nx = (60.0 / dx) + 1        ! Domain width is 60 km total (-30 to 30 km)

        case("HALFAR-MED")

            grid_name = "HALFAR-MED"
            nx = (300.0 / dx) + 1        ! Domain width is 300 km total (-150 to 150 km)

        case DEFAULT 
            ! EISMINT1, EISMINT2 and Bueler test grid setup 

            grid_name = "EISMINT"
            nx = (1500.0 / dx) + 1     ! Domain width is 1500 km total (-750 to 750 km)

    end select 

    ! === Initialize ice sheet model =====

    ! General initialization of yelmo constants (used globally)
    call yelmo_global_init(path_const)

    ! Next define grid 
    call yelmo_init_grid(yelmo1%grd,grid_name,units="km",dx=dx,nx=nx,dy=dx,ny=nx)



    ! Initialize data objects (without loading topography, which will be defined inline below)
    call yelmo_init(yelmo1,filename=path_par,grid_def="none",time=time_init,load_topo=.FALSE.,domain=domain,grid_name=grid_name)
    

    ! Update parameter values with EISMINT choices 
    yelmo1%dyn%par%use_ssa    = with_ssa 
    yelmo1%tpo%par%topo_fixed = topo_fixed 



    ! Initialize Bueler test type 
    call bueler_init(buel,yelmo1%grd%nx,yelmo1%grd%ny)


    ! === Define initial topography =====

    select case(trim(experiment))

        case("mismip")
            ! Add a marine bed for testing ice shelves following MISMIP
        
            yelmo1%bnd%z_bed  = 720.0 - 778.50*(sqrt((yelmo1%grd%x*1e-3)**2+(yelmo1%grd%y*1e-3)**2))/750.0
            yelmo1%tpo%now%H_ice  = 100.0
            yelmo1%tpo%now%z_srf  = yelmo1%bnd%z_bed + yelmo1%tpo%now%H_ice

!             where(yelmo1%bnd%z_bed .lt. 0.0) yelmo1%bnd%smb = 0.0 
        
        case DEFAULT 
            ! EISMINT1, EISMINT2, BUELER 

            yelmo1%tpo%now%z_srf  = yelmo1%bnd%z_bed + yelmo1%tpo%now%H_ice 

    end select 

    ! ==== READ STEADY-STATE TOPOGRAPHY FROM HEIKO'S RUN
    if (topo_fixed) then

        ! Consistency check 
        if (trim(experiment) .ne. "moving") then 
            write(*,*) "yelmo_benchmarks:: Error: topo_fixed=True can only be used&
            & for the EISMINT1 'moving' experiment."
            write(*,*) "experiment = ", trim(experiment)
            stop 
        end if 

        yelmo1%bnd%z_bed  = 0.0 
        call nc_read(topo_fixed_file,"Hi",yelmo1%tpo%now%H_ice)
        where(yelmo1%tpo%now%H_ice.lt.1.0) yelmo1%tpo%now%H_ice = 0.0 
        yelmo1%tpo%now%z_srf  = yelmo1%bnd%z_bed + yelmo1%tpo%now%H_ice 
        
    end if 
    
    ! Load boundary values

    yelmo1%bnd%z_sl     = 0.0
    yelmo1%bnd%bmb_shlf = 0.0  
    yelmo1%bnd%T_shlf   = T0  
    yelmo1%bnd%H_sed    = 0.0 

    select case(trim(experiment))

        case("BUELER-A")

            ! Initialize BUELER-A 
            call bueler_test_AE(buel%H_ice,buel%mbal,buel%u_b,yelmo1%grd%x,yelmo1%grd%y, &
                                L=750.0_prec,mbal0=0.3_prec,A=1e-16_prec,n=3.0_prec,rho_ice=rho_ice,g=g,mu_max=0.0_prec)

            yelmo1%bnd%T_srf = 223.15 
            yelmo1%bnd%Q_geo = 42.0 
            yelmo1%bnd%smb   = buel%mbal 
            
            ! Update boundary mask to limit ice growth 
            where(buel%H_ice .eq. 0.0) yelmo1%bnd%ice_allowed = .FALSE. 

            ! Update topography to match exact solution to start 
            yelmo1%tpo%now%H_ice  = buel%H_ice
            yelmo1%tpo%now%z_srf  = yelmo1%bnd%z_bed + yelmo1%tpo%now%H_ice
            
        case("BUELER-B")

            ! Initialize BUELER-B 
            call bueler_test_BC(buel%H_ice,buel%mbal,buel%u_b,yelmo1%grd%x,yelmo1%grd%y, &
                        time=0.0_prec,R0=750.0_prec,H0=3600.0_prec,lambda=0.0_prec,n=3.0_prec,A=1e-16_prec,rho_ice=rho_ice,g=g)

            yelmo1%bnd%T_srf = 223.15 
            yelmo1%bnd%Q_geo = 42.0 
            yelmo1%bnd%smb   = buel%mbal 
            
            ! Update boundary mask to limit ice growth 
            !where(buel%H_ice .eq. 0.0) yelmo1%bnd%ice_allowed = .FALSE. 

            ! Update topography to match exact solution to start 
            yelmo1%tpo%now%H_ice  = buel%H_ice
            yelmo1%tpo%now%z_srf  = yelmo1%bnd%z_bed + yelmo1%tpo%now%H_ice

        case("HALFAR")

            ! Initialize BUELER-B (but with HALFAR conditions)
            call bueler_test_BC(buel%H_ice,buel%mbal,buel%u_b,yelmo1%grd%x,yelmo1%grd%y, &
                        time=0.0_prec,R0=21.2132_prec,H0=707.1_prec,lambda=0.0_prec,n=3.0_prec,A=1e-16_prec,rho_ice=rho_ice,g=g)

            yelmo1%bnd%T_srf = 223.15 
            yelmo1%bnd%Q_geo = 42.0 
            yelmo1%bnd%smb   = buel%mbal 
            
            ! Update boundary mask to limit ice growth 
            !where(buel%H_ice .eq. 0.0) yelmo1%bnd%ice_allowed = .FALSE. 

            ! Update topography to match exact solution to start 
            yelmo1%tpo%now%H_ice  = buel%H_ice
            yelmo1%tpo%now%z_srf  = yelmo1%bnd%z_bed + yelmo1%tpo%now%H_ice
            
        case("HALFAR-MED")

            ! Initialize BUELER-B (but with conditions between Bueler-B and HALFAR - not too big, not too small)
            call bueler_test_BC(buel%H_ice,buel%mbal,buel%u_b,yelmo1%grd%x,yelmo1%grd%y, &
                        time=0.0_prec,R0=100.0_prec,H0=800.0_prec,lambda=0.0_prec,n=3.0_prec,A=1e-16_prec,rho_ice=rho_ice,g=g)

            yelmo1%bnd%T_srf = 223.15 
            yelmo1%bnd%Q_geo = 42.0 
            yelmo1%bnd%smb   = buel%mbal 
            
            ! Update boundary mask to limit ice growth 
            !where(buel%H_ice .eq. 0.0) yelmo1%bnd%ice_allowed = .FALSE. 

            ! Update topography to match exact solution to start 
            yelmo1%tpo%now%H_ice  = buel%H_ice
            yelmo1%tpo%now%z_srf  = yelmo1%bnd%z_bed + yelmo1%tpo%now%H_ice
        
        case("mismip") 

            ! Initialize mismip boundary values 
            call mismip3D_boundaries(yelmo1%bnd%T_srf,yelmo1%bnd%smb,yelmo1%bnd%Q_geo,experiment="Stnd") 

        case DEFAULT 
            ! EISMINT 

            ! Initialize eismint boundary values 
            call eismint_boundaries(yelmo1%bnd%T_srf,yelmo1%bnd%smb,yelmo1%bnd%Q_geo, &
                                    yelmo1%grd%x,yelmo1%grd%y,yelmo1%tpo%now%H_ice, &
                                    experiment=experiment,time=0.0_prec,period=period,dT_test=dT_test)


    end select 

    ! Call bueler_compare once to initialize comparison fields (even though it is not currently used for EISMINT sims)
    call bueler_compare(buel,yelmo1%tpo%now%H_ice,dx=yelmo1%grd%dx)

    ! Actually use a small circular ice sheet to get things rolling...
    ! Set ice thickness to a circle of low ice thickness to start
    ! (for testing only)
    if (.FALSE.) then

        yelmo1%tpo%now%H_ice  = 0.0
        where(yelmo1%bnd%smb .gt. 0.0) 
            yelmo1%tpo%now%H_ice = max(0.0, 1000.0 + (3000.0-1000.0)*(750.0-sqrt((yelmo1%grd%x*1e-3)**2+(yelmo1%grd%y*1e-3)**2))/750.0)
!             yelmo1%tpo%now%H_ice = max(0.0, 10.0 + (300.0-100.0)*(750.0-sqrt((yelmo1%grd%x*1e-3)**2+(yelmo1%grd%y*1e-3)**2))/750.0)
        end where 
        yelmo1%tpo%now%z_srf  = yelmo1%bnd%z_bed + yelmo1%tpo%now%H_ice
            
    end if 

    ! Check boundary values 
    call yelmo_print_bound(yelmo1%bnd)

    ! Initialize state variables (dyn,therm,mat)
    call yelmo_init_state(yelmo1,path_par,time=time_init,thrm_method="robin")


    ! == Write initial state ==
     
    ! 2D file 
    call yelmo_write_init(yelmo1,file2D,time_init=time_init,units="years")
    call write_step_2D(yelmo1,file2D,time=time_init)  
    
    ! 1D file 
    call write_yreg_init(yelmo1,file1D,time_init=time_init,units="years",mask=yelmo1%bnd%ice_allowed)
    call write_yreg_step(yelmo1%reg,file1D,time=time_init) 

    ! Comparison file 
    call yelmo_write_init(yelmo1,file_compare,time_init=time_init,units="years")
    call write_step_2D_bueler(yelmo1,buel,file_compare,time_init)
    
    if (dyn_fixed) then 
        ! Set yelmo parameter to fix dynamics
        yelmo1%dyn%par%solver = "fixed"
    end if 



    ! Advance timesteps
    do n = 1, ceiling((time_end-time_init)/dtt)

        ! Get current time 
        time = time_init + n*dtt
        
        ! == Yelmo ice sheet ===================================================
        call yelmo_update(yelmo1,time)
        
        ! == Update boundaries 
        select case(trim(experiment))

            case("BUELER-A")

                ! Pass - no update to perform here 

                
            case("BUELER-B")

                ! Update BUELER-B profile
                call bueler_test_BC(buel%H_ice,buel%mbal,buel%u_b,yelmo1%grd%x,yelmo1%grd%y, &
                      time=time,R0=750.0_prec,H0=3600.0_prec,lambda=0.0_prec,n=3.0_prec,A=1e-16_prec,rho_ice=rho_ice,g=g)

            case("HALFAR")

                ! Update BUELER-B (but with HALFAR conditions)
                call bueler_test_BC(buel%H_ice,buel%mbal,buel%u_b,yelmo1%grd%x,yelmo1%grd%y, &
                      time=time,R0=21.2132_prec,H0=707.1_prec,lambda=0.0_prec,n=3.0_prec,A=1e-16_prec,rho_ice=rho_ice,g=g)

            case("HALFAR-MED")

                ! Initialize BUELER-B (but with conditions between Bueler-B and HALFAR - not too big, not too small)
                call bueler_test_BC(buel%H_ice,buel%mbal,buel%u_b,yelmo1%grd%x,yelmo1%grd%y, &
                            time=time,R0=100.0_prec,H0=800.0_prec,lambda=0.0_prec,n=3.0_prec,A=1e-16_prec,rho_ice=rho_ice,g=g)
            
            case("mismip") 

                ! Initialize mismip boundary values 
                call mismip3D_boundaries(yelmo1%bnd%T_srf,yelmo1%bnd%smb,yelmo1%bnd%Q_geo,experiment="Stnd") 

            case DEFAULT 

                call eismint_boundaries(yelmo1%bnd%T_srf,yelmo1%bnd%smb,yelmo1%bnd%Q_geo, &
                                yelmo1%grd%x,yelmo1%grd%y,yelmo1%tpo%now%H_ice, &
                                experiment=experiment,time=time,period=period,dT_test=dT_test)

        end select 

        ! Make comparison if relevant 
        select case(trim(experiment))
            case("BUELER-A","BUELER-B","HALFAR","HALFAR-MED")

                call bueler_compare(buel,yelmo1%tpo%now%H_ice,dx=yelmo1%grd%dx)
                
            case DEFAULT 
            ! Pass- not a bueler test... 
        end select 

        ! == MODEL OUTPUT =======================================================
        if (mod(nint(time*100),nint(dt2D_out*100))==0) then 
            call write_step_2D(yelmo1,file2D,time=time) 
            call write_step_2D_bueler(yelmo1,buel,file_compare,time)   
        end if 

        if (mod(nint(time*100),nint(dt1D_out*100))==0) then 
            call write_yreg_step(yelmo1%reg,file1D,time=time) 
        end if 

        if (mod(time,10.0)==0 .and. (.not. yelmo_log)) then
            write(*,"(a,f14.4)") "yelmo:: time = ", time
        end if 

    end do 




    ! Write summary 
    write(*,*) "====== "//trim(domain)//"-"//trim(experiment)//" ======="
    write(*,*) "nz_aa, H0 = ", yelmo1%par%nz_aa, maxval(yelmo1%tpo%now%H_ice)

    ! Make comparison if relevant 
    select case(trim(experiment))
        case("BUELER-A","BUELER-B","HALFAR","HALFAR-MED")

            write(*,"(a,3f10.2,10g12.3)") trim(experiment), time, dtt, yelmo1%grd%dx*1e-3, &
                                            buel%rmse_H_ice, buel%err_H0, buel%err_max_H_ice, buel%err_V_ice

        case DEFAULT 
        ! Pass- not a bueler test... 
    end select 


    ! Write a restart file too
    call yelmo_restart_write(yelmo1,file_restart,time=time)

    ! Finalize program
    call yelmo_end(yelmo1,time=time)

    ! Stop timing 
    call cpu_time(finish)

    print '("Time = ",f12.3," min.")', (finish-start)/60.0 

contains
    
    subroutine write_step_2D(ylmo,filename,time,buel)

        implicit none 
        
        type(yelmo_class), intent(IN) :: ylmo
        character(len=*),  intent(IN) :: filename
        real(prec), intent(IN) :: time
        type(bueler_test_type), intent(IN), optional :: buel 

        ! Local variables
        integer    :: ncid, n, i, j, nx, ny  
        real(prec) :: time_prev 
        real(prec), allocatable :: sym(:,:) 

        nx = ylmo%tpo%par%nx 
        ny = ylmo%tpo%par%ny 

        allocate(sym(nx,ny)) 

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
        call nc_write(filename,"dt_avg",ylmo%par%dt_avg,units="yr",long_name="Average timestep", &
                      dim1="time",start=[n],count=[1],ncid=ncid)
        call nc_write(filename,"eta_avg",ylmo%par%eta_avg,units="m a-1",long_name="Average eta (maximum PC truncation error)", &
                      dim1="time",start=[n],count=[1],ncid=ncid)
        
        ! Time step limits 
        call nc_write(filename,"dt_adv",ylmo%par%dt_adv,units="a",long_name="Advective timestep", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dt_diff",ylmo%par%dt_diff,units="a",long_name="Diffusive timestep", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"dt_adv3D",ylmo%par%dt_adv3D,units="a",long_name="Advective timestep", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        ! == yelmo_topography ==
        call nc_write(filename,"H_ice",ylmo%tpo%now%H_ice,units="m",long_name="Ice thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"H_margin",ylmo%tpo%now%H_margin,units="m",long_name="Margin ice thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf",ylmo%tpo%now%z_srf,units="m",long_name="Surface elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_bed",ylmo%tpo%now%mask_bed,units="",long_name="Bed mask", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"H_grnd",ylmo%tpo%now%H_grnd,units="m",long_name="Ice thickness overburden", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mb_applied",ylmo%tpo%now%mb_applied,units="m/a",long_name="Actual ice mass balance applied", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"dzsrfdt",ylmo%tpo%now%dzsrfdt,units="m/a",long_name="Surface elevation change", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dHicedt",ylmo%tpo%now%dHicedt,units="m/a",long_name="Ice thickness change", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"dHicedx",ylmo%tpo%now%dHicedx,units="m/m",long_name="Ice thickness gradient (acx)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dHicedy",ylmo%tpo%now%dHicedy,units="m/m",long_name="Ice thickness gradient (acy)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"f_grnd",ylmo%tpo%now%f_grnd,units="1",long_name="Grounded fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_grnd_acx",ylmo%tpo%now%f_grnd_acx,units="1",long_name="Grounded fraction (acx)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_grnd_acy",ylmo%tpo%now%f_grnd_acy,units="1",long_name="Grounded fraction (acy)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_ice",ylmo%tpo%now%f_ice,units="1",long_name="Ice-covered fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"N_eff",ylmo%dyn%now%N_eff,units="Pa",long_name="Effective pressure", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"calv",ylmo%tpo%now%calv,units="m/a ice equiv.",long_name="Calving rate", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! == yelmo_thermodynamics ==
        call nc_write(filename,"enth",ylmo%thrm%now%enth,units="J m-3",long_name="Ice enthalpy", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"T_ice",ylmo%thrm%now%T_ice,units="K",long_name="Ice temperature", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"omega",ylmo%thrm%now%omega,units="--",long_name="Ice water content", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"T_prime",ylmo%thrm%now%T_ice-ylmo%thrm%now%T_pmp,units="deg C",long_name="Homologous ice temperature", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
!         call nc_write(filename,"T_pmp",ylmo%thrm%now%T_pmp,units="K",long_name="Ice pressure melting point (pmp)", &
!                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        
        call nc_write(filename,"T_prime_b",ylmo%thrm%now%T_prime_b,units="deg C",long_name="Homologous basal ice temperature", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        sym = ylmo%thrm%now%T_prime_b(nx:1:-1,:) - ylmo%thrm%now%T_prime_b
        call nc_write(filename,"T_prime_b_sym",sym,units="deg C",long_name="Homologous basal ice temperature symmetry check", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"f_pmp",ylmo%thrm%now%f_pmp,units="1",long_name="Fraction of grid point at pmp", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"Q_ice_b",ylmo%thrm%now%Q_ice_b,units="J a-1 m-2",long_name="Basal ice heat flux", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"Q_strn",ylmo%thrm%now%Q_strn/(rho_ice*ylmo%thrm%now%cp),units="K a-1",long_name="Strain heating", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        call nc_write(filename,"Q_b",ylmo%thrm%now%Q_b,units="J a-1 m-2",long_name="Basal frictional heating", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"bmb_grnd",ylmo%thrm%now%bmb_grnd,units="m/a ice equiv.",long_name="Basal mass balance (grounded)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"H_w",ylmo%thrm%now%H_w,units="m",long_name="Basal water pressure", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"H_cts",ylmo%thrm%now%H_cts,units="m",long_name="Height of CTS", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! == yelmo_material ==
!         call nc_write(filename,"enh_bar",ylmo%mat%now%enh_bar,units="1",long_name="Vertically averaged enhancement factor", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"visc_int",ylmo%mat%now%visc_int,units="Pa a m",long_name="Vertically integrated viscosity", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"ATT",ylmo%mat%now%ATT,units="a^-1 Pa^-3",long_name="Rate factor", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"visc",ylmo%mat%now%visc,units="Pa a",long_name="Viscosity", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        
!         call nc_write(filename,"dep_time",ylmo%mat%now%dep_time,units="a",long_name="Ice deposition time", &
!                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"age",time-ylmo%mat%now%dep_time,units="a",long_name="Ice age", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        ! == yelmo_dynamics ==

        call nc_write(filename,"ssa_mask_acx",ylmo%dyn%now%ssa_mask_acx,units="1",long_name="SSA mask (acx)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"ssa_mask_acy",ylmo%dyn%now%ssa_mask_acy,units="1",long_name="SSA mask (acy)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"cf_ref",ylmo%dyn%now%cf_ref,units="--",long_name="Bed friction scalar", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"c_bed",ylmo%dyn%now%c_bed,units="Pa",long_name="Bed friction coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"beta",ylmo%dyn%now%beta,units="Pa a m-1",long_name="Basal friction coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"visc_eff",ylmo%dyn%now%visc_eff,units="Pa a m",long_name="Effective viscosity (SSA)", &
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
        
!         call nc_write(filename,"ux_i_bar",ylmo%dyn%now%ux_i_bar,units="m/a",long_name="Internal shear velocity (x)", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"uy_i_bar",ylmo%dyn%now%uy_i_bar,units="m/a",long_name="Internal shear velocity (y)", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"uxy_i_bar",ylmo%dyn%now%uxy_i_bar,units="m/a",long_name="Internal shear velocity magnitude", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"ux_b",ylmo%dyn%now%ux_b,units="m/a",long_name="Basal sliding velocity (x)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy_b",ylmo%dyn%now%uy_b,units="m/a",long_name="Basal sliding velocity (y)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_b",ylmo%dyn%now%uxy_b,units="m/a",long_name="Basal sliding velocity magnitude", &
                     dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"ux_bar",ylmo%dyn%now%ux_bar,units="m/a",long_name="Vertically averaged velocity (x)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy_bar",ylmo%dyn%now%uy_bar,units="m/a",long_name="Vertically averaged velocity (y)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_bar",ylmo%dyn%now%uxy_bar,units="m/a",long_name="Vertically averaged velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"ux_s",ylmo%dyn%now%ux_s,units="m/a",long_name="Surface velocity (x)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy_s",ylmo%dyn%now%uy_s,units="m/a",long_name="Surface velocity (y)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_s",ylmo%dyn%now%uxy_s,units="m/a",long_name="Surface velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        sym = ylmo%dyn%now%uxy_bar(nx:1:-1,:) - ylmo%dyn%now%uxy_bar
        call nc_write(filename,"uxy_bar_sym",sym,units="m/a",long_name="Vertically integrated velocity magnitude symmetry check", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"qq_acx",ylmo%dyn%now%qq_acx,units="m^3/a",long_name="Ice flux (acx-nodes)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"qq_acy",ylmo%dyn%now%qq_acy,units="m^3/a",long_name="Ice flux (acy-nodes)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"qq",ylmo%dyn%now%qq,units="m^3/a",long_name="Ice flux", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
!         call nc_write(filename,"dd_ab",ylmo%dyn%now%dd_ab,units="m2 a-1",long_name="Diffusivity", &
!                        dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"dd_ab_bar",ylmo%dyn%now%dd_ab_bar,units="m2 a-1",long_name="Diffusivity", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"ux",ylmo%dyn%now%ux,units="m/a",long_name="Horizontal velocity (x)", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"uy",ylmo%dyn%now%uy,units="m/a",long_name="Horizontal velocity (y)", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"uxy",ylmo%dyn%now%uxy,units="m/a",long_name="Horizontal velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"uz",ylmo%dyn%now%uz,units="m/a",long_name="Vertical velocity", &
                      dim1="xc",dim2="yc",dim3="zeta_ac",dim4="time",start=[1,1,1,n],ncid=ncid)

        call nc_write(filename,"f_vbvs",ylmo%dyn%now%f_vbvs,units="1",long_name="Basal to surface velocity fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_shear_bar",ylmo%mat%now%f_shear_bar,units="1",long_name="Vertically averaged shearing fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"de",ylmo%mat%now%strn%de,units="a^-1",long_name="Strain rate", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        ! == yelmo_bound ==

        call nc_write(filename,"z_bed",ylmo%bnd%z_bed,units="m",long_name="Bedrock elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
!         call nc_write(filename,"z_sl",ylmo%bnd%z_sl,units="m",long_name="Sea level rel. to present", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"H_sed",ylmo%bnd%H_sed,units="m",long_name="Sediment thickness", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"smb",ylmo%bnd%smb,units="m/a ice equiv.",long_name="Surface mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"T_srf",ylmo%bnd%T_srf,units="K",long_name="Surface temperature", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
     
        call nc_write(filename,"bmb",ylmo%tpo%now%bmb,units="m/a ice equiv.",long_name="Basal mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"bmb_shlf",ylmo%bnd%bmb_shlf,units="m/a ice equiv.",long_name="Basal mass balance (shelf)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"Q_geo",ylmo%bnd%Q_geo,units="mW/m^2",long_name="Geothermal heat flux", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        

        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine write_step_2D

    subroutine write_step_2D_bueler(ylmo,buel,filename,time)

        implicit none 
        
        type(yelmo_class),      intent(IN) :: ylmo
        type(bueler_test_type), intent(IN) :: buel 
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

        call nc_write(filename,"H_ice",ylmo%tpo%now%H_ice,units="m",long_name="Ice thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"H_ice_target",buel%H_ice,units="m",long_name="Target ice thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"err_H_ice",buel%err_H_ice,units="m",long_name="Ice thickness error", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"rmse_H_ice",buel%rmse_H_ice,units="m",long_name="Ice thickness RMSE", &
                      dim1="time",start=[n],count=[1],ncid=ncid)
        
        call nc_write(filename,"err_H0",buel%err_H0,units="m",long_name="Summit error", &
                      dim1="time",start=[n],count=[1],ncid=ncid)

        call nc_write(filename,"err_max_H_ice",buel%err_max_H_ice,units="m",long_name="Maximum error", &
                      dim1="time",start=[n],count=[1],ncid=ncid)

        call nc_write(filename,"V_ice",buel%V_ice_mod,units="1e6 km^3",long_name="Ice volume", &
                      dim1="time",start=[n],count=[1],ncid=ncid)
        call nc_write(filename,"V_ice_target",buel%V_ice_mod,units="1e6 km^3",long_name="Target ice volume", &
                      dim1="time",start=[n],count=[1],ncid=ncid)
        call nc_write(filename,"err_V_ice",buel%err_V_ice,units="km^3",long_name="Ice volume error", &
                      dim1="time",start=[n],count=[1],ncid=ncid)

        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine write_step_2D_bueler 

end program yelmo_benchmarks 



! Extra code for ratefactor testing

!     ! Simulations to test rate factor feedback
!     logical :: testing_ratefactor = .FALSE. 
!     real(prec) :: mod_time_1, mod_time_2, mod_time_3 
!     real(prec) :: rad_max, rad_min, rad_now
!     real(prec) :: dT_max, dT_min, dT_now
!     real(prec) :: dsmb_max, dsmb_min, dsmb_now 
!     real(prec) :: f_lin 

!     ! Settings for transient rate factor test. 
!     mod_time_1 = 100e3 
!     mod_time_2 = 150e3 
!     mod_time_3 = 200e3
!     rad_max    = 450.0   ! Consistent with original boundary condition 
!     rad_min    = 150.0   ! Minimum equilibrium line radius to reach (ice sheet will shrink)
!     dT_min     =   0.0 
!     dT_max     =  10.0 
!     dsmb_min   =   0.0 
!     dsmb_max   =   0.5 

!     rad_now    = 450.0   ! Set radius now to default value 
!     dT_now     =   0.0 
!     dsmb_now   =   0.0 

!                 if (testing_ratefactor) then 


!                     if (time .ge. mod_time_1 .and. time .le. mod_time_2) then

!                         f_lin = (time - mod_time_1) / (mod_time_2 - mod_time_1)
!                         !rad_now = rad_max + f_lin * (rad_min - rad_max)
!                         !dT_now = dT_min + f_lin * (dT_max - dT_min)
                        
!                         dsmb_now = dsmb_min + f_lin * (dsmb_max - dsmb_min)

!                     else if (time .gt. mod_time_2 .and. time .le. mod_time_3) then

!                         f_lin = (time - mod_time_2) / (mod_time_3 - mod_time_2)
!                         !rad_now = rad_min + f_lin * (rad_max - rad_min)
!                         !dT_now = dT_max + f_lin * (dT_min - dT_max)
                        
!                         dsmb_now = dsmb_max + f_lin * (dsmb_min - dsmb_max)
                        
!                     else

!                         !dT_now   = dT_min 
!                         dsmb_now = dsmb_min 

!                     end if 

! !                     if (time .ge. mod_time_1) yelmo1%mat%par%rf_method = -1 

!                     write(*,*) "testing_ratefactor: ", time, rad_now, dT_now, dsmb_now, yelmo1%mat%par%rf_method  

!                 end if 


