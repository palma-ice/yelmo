

program yelmo_benchmarks

    use nml 
    use ncio  
    use yelmo 

    use deformation 

    use ice_benchmarks 
    use mismip3D 
    use timestepping

    implicit none 

    type(tstep_class)      :: ts
    type(yelmo_class)      :: yelmo1
    type(bueler_test_type) :: buel 
    
    character(len=56)  :: domain    
    character(len=256) :: outfldr, file2D, file1D, file_compare
    character(len=256) :: file_restart
    character(len=512) :: path_par 
    character(len=56)  :: experiment
    logical    :: with_bumps, low_z_sl 
    real(wp) :: time_init, time_end, dtt, dt2D_out, dt1D_out
    real(wp) :: period, dt_test, alpha, omega, L, amp 
    real(wp) :: bumps_L, bumps_A 
    integer    :: n  

    character(len=56) :: thrm_method_default 
    character(len=56) :: rock_method_default 

    character(len=56) :: grid_name
    real(wp) :: dx, x0  
    integer    :: nx  

    real(8) :: cpu_start_time, cpu_end_time, cpu_dtime  
    
    ! Start timing 
    call yelmo_cpu_time(cpu_start_time)

    ! Assume program is running from the output folder
    outfldr = "./"

    ! Determine the parameter file from the command line 
    call yelmo_load_command_line_args(path_par)
    !path_par   = trim(outfldr)//"yelmo_EISMINT.nml" 

    ! Define input and output locations 
    file1D       = trim(outfldr)//"yelmo1D.nc"
    file2D       = trim(outfldr)//"yelmo2D.nc"
    file_compare = trim(outfldr)//"yelmo_compare.nc"
    file_restart = trim(outfldr)//"yelmo_restart.nc"

    ! Define the domain, grid and experiment from parameter file
    call nml_read(path_par,"ctrl","domain",       domain)        ! EISMINT1, EISMINT2
    call nml_read(path_par,"ctrl","experiment",   experiment)    ! "fixed", "moving", "mismip", "EXPA", "EXPB", "BUELER-A"
    call nml_read(path_par,"ctrl","dx",           dx)            ! [km] Grid resolution 

    ! Timing parameters 
    call nml_read(path_par,"ctrl","time_init",    time_init)     ! [yr] Starting time
    call nml_read(path_par,"ctrl","time_end",     time_end)      ! [yr] Ending time
    call nml_read(path_par,"ctrl","dtt",          dtt)           ! [yr] Main loop time step 
    call nml_read(path_par,"ctrl","dt2D_out",     dt2D_out)      ! [yr] Frequency of 2D output 
    dt1D_out = dtt  ! Set 1D output to frequency of main loop timestep 

    ! Settings for transient EISMINT1 experiments 
    call nml_read(path_par,"ctrl","period",       period)        ! [yr] for transient experiments 
    call nml_read(path_par,"ctrl","dT_test",      dT_test)       ! [K] for test experiments  
    
    ! call nml_read(path_par,"ctrl","with_bumps",with_bumps)       ! Bedrock with sin bumps?
    ! call nml_read(path_par,"ctrl","bumps_L",bumps_L)             ! [km] Length scale of bumps
    ! call nml_read(path_par,"ctrl","bumps_A",bumps_A)             ! [m]  Amplitude of bumps
    with_bumps = .FALSE. 
    low_z_sl   = .FALSE. 

    ! === Initialize timestepping ===
    
    call tstep_init(ts,time_init,time_end,method="const",units="year", &
                                            time_ref=1950.0_wp,const_rel=0.0_wp,const_cal=1950.0_wp)

    ! Define the model domain based on the experiment we are running
    select case(trim(experiment))

        case("mismip")
            ! EISMINT1 experiment with modified EISMINT1-mismip grid proposed by Heiko

            grid_name = "EISMINT1-mismip"
            nx = 73  

            call yelmo_init_grid(yelmo1%grd,grid_name,units="km",dx=dx,nx=nx,dy=dx,ny=nx)

        case("BUELER-A","BUELER-B")

            grid_name = "EISMINT-EXT"
            nx = (2000.0 / dx) + 1          ! Domain width is 2000 km total (-1000 to 1000 km)
            if (mod(nx,2).eq.0) nx=nx+1     ! Make sure that nx is odd, so that we get one grid point right at dome x=0,y=0

            call yelmo_init_grid(yelmo1%grd,grid_name,units="km",dx=dx,nx=nx,dy=dx,ny=nx)

        case("HALFAR")

            grid_name = "HALFAR"
            nx = (80.0 / dx) + 1            ! Domain width is 60 km total (-30 to 30 km)
            if (mod(nx,2).eq.0) nx=nx+1     ! Make sure that nx is odd, so that we get one grid point right at dome x=0,y=0

            call yelmo_init_grid(yelmo1%grd,grid_name,units="km",dx=dx,nx=nx,dy=dx,ny=nx)

        case("HALFAR-MED")

            grid_name = "HALFAR-MED"
            nx = (300.0 / dx) + 1           ! Domain width is 300 km total (-150 to 150 km)
            if (mod(nx,2).eq.0) nx=nx+1     ! Make sure that nx is odd, so that we get one grid point right at dome x=0,y=0

            call yelmo_init_grid(yelmo1%grd,grid_name,units="km",dx=dx,nx=nx,dy=dx,ny=nx)
        
        case("HEINO")

            grid_name = "HEINO"
            nx = (500.0 / dx) + 1           ! Domain width is 500 km total (0 to 500 km)
            x0 = 0.0

            call yelmo_init_grid(yelmo1%grd,grid_name,units="km",x0=x0,dx=dx,nx=nx,y0=x0,dy=dx,ny=nx)
        
        case DEFAULT 
            ! EISMINT1, EISMINT2, dome and Bueler test grid setup 

            grid_name = "EISMINT"
            nx = (1500.0 / dx) + 1          ! Domain width is 1500 km total (-750 to 750 km)
            if (mod(nx,2).eq.0) nx=nx+1     ! Make sure that nx is odd, so that we get one grid point right at dome x=0,y=0

            call yelmo_init_grid(yelmo1%grd,grid_name,units="km",dx=dx,nx=nx,dy=dx,ny=nx)
            
    end select 

    ! === Initialize ice sheet model =====

    ! Initialize data objects (without loading topography, which will be defined inline below)
    call yelmo_init(yelmo1,filename=path_par,grid_def="none",time=ts%time,load_topo=.FALSE.,domain=domain,grid_name=grid_name)
    
    ! Initialize Bueler test type 
    call bueler_init(buel,yelmo1%grd%nx,yelmo1%grd%ny)

    ! === Define initial topography =====

    select case(trim(experiment))

        case("mismip")
            ! Add a marine bed for testing ice shelves following MISMIP
        
            yelmo1%bnd%z_bed  = 720.0 - 778.50*(sqrt((yelmo1%grd%x*1e-3)**2+(yelmo1%grd%y*1e-3)**2))/750.0
            
            !yelmo1%tpo%now%H_ice  = 100.0
            call dome_init(yelmo1%tpo%now%H_ice,yelmo1%grd%x,yelmo1%grd%y,R0=0.5_wp,H0=2000.0_wp, &
                                H0_shlf=150.0_wp,rmax_shlf=0.8_wp) 

            yelmo1%tpo%now%z_srf  = yelmo1%bnd%z_bed + yelmo1%tpo%now%H_ice

!             where(yelmo1%bnd%z_bed .lt. 0.0) yelmo1%bnd%smb = 0.0 
        
        case("dome") 
            ! Define a radially symmetric ice sheet to start, based on an ellipsoidal (square root) profile
            ! following Lipscomb et al. (2019) and CISMv2.1. 

            call dome_init(yelmo1%tpo%now%H_ice,yelmo1%grd%x,yelmo1%grd%y,R0=0.5_wp,H0=2000.0_wp)

            yelmo1%bnd%z_bed     = 0.0_wp 
            yelmo1%tpo%now%z_srf = yelmo1%bnd%z_bed + yelmo1%tpo%now%H_ice

        case DEFAULT 
            ! EISMINT1, EISMINT2, HALFAR, BUELER 

if (.TRUE.) then
            ! Use flat bed as expected by EISMINT experiments
            yelmo1%bnd%z_bed     = 0.0_wp
else
            ! Use MISMIP-like sloping bedrock
            yelmo1%bnd%z_bed  = 220.0 - 778.50*(sqrt((yelmo1%grd%x*1e-3)**2+(yelmo1%grd%y*1e-3)**2))/750.0
end if

            yelmo1%tpo%now%H_ice = 0.0_wp 
            yelmo1%tpo%now%z_srf = yelmo1%bnd%z_bed + yelmo1%tpo%now%H_ice 

    end select 

    if (with_bumps) then 

        L     = bumps_L                 ! [km] Length scale 
        amp   = bumps_A                 ! [m] Amplitude 
        alpha = 0.5*pi/180.0_wp       ! [rad] 
        omega = 2.0_wp*pi / (L*1e3)   ! [rad/m]

        !yelmo1%tpo%now%z_srf = -yelmo1%grd%x * tan(alpha)
        yelmo1%bnd%z_bed     = 0.0 + amp * sin(omega*yelmo1%grd%x) * sin(omega*yelmo1%grd%y)

    end if 

    ! Load boundary values

    yelmo1%bnd%z_sl     = 0.0
    yelmo1%bnd%bmb_shlf = 0.0  
    yelmo1%bnd%T_shlf   = yelmo1%bnd%c%T0  
    yelmo1%bnd%H_sed    = 0.0 

    if (with_bumps .or. low_z_sl) then 
        yelmo1%bnd%z_sl     = -5000.0
    end if 

    select case(trim(experiment))

        case("BUELER-A")

            ! Initialize BUELER-A 
            call bueler_test_AE(buel%H_ice,buel%mbal,buel%u_b,yelmo1%grd%x,yelmo1%grd%y, &
                                L=750.0_wp,mbal0=0.3_wp,A=1e-16_wp,n=3.0_wp, &
                                rho_ice=yelmo1%bnd%c%rho_ice,g=yelmo1%bnd%c%g,mu_max=0.0_wp)

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
                        time=0.0_wp,R0=750.0_wp,H0=3600.0_wp,lambda=0.0_wp,n=3.0_wp,A=1e-16_wp, &
                        rho_ice=yelmo1%bnd%c%rho_ice,g=yelmo1%bnd%c%g)

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
                        time=0.0_wp,R0=21.2132_wp,H0=707.1_wp,lambda=0.0_wp,n=3.0_wp,A=1e-16_wp, &
                        rho_ice=yelmo1%bnd%c%rho_ice,g=yelmo1%bnd%c%g)

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
                        time=0.0_wp,R0=100.0_wp,H0=800.0_wp,lambda=0.0_wp,n=3.0_wp,A=1e-16_wp, &
                        rho_ice=yelmo1%bnd%c%rho_ice,g=yelmo1%bnd%c%g)

            yelmo1%bnd%T_srf = 223.15 
            yelmo1%bnd%Q_geo = 42.0 
            yelmo1%bnd%smb   = buel%mbal 
            
            ! Update boundary mask to limit ice growth 
            !where(buel%H_ice .eq. 0.0) yelmo1%bnd%ice_allowed = .FALSE. 

            ! Update topography to match exact solution to start 
            yelmo1%tpo%now%H_ice  = buel%H_ice
            yelmo1%tpo%now%z_srf  = yelmo1%bnd%z_bed + yelmo1%tpo%now%H_ice
        
        case("mismip")
            
            ! To do: define calv_mask!!
            
            ! Set conditions similar to EISMINT2-EXPA with smaller radius 
            call dome_boundaries(yelmo1%bnd%T_srf,yelmo1%bnd%smb,yelmo1%bnd%Q_geo, &
                            yelmo1%grd%x,yelmo1%grd%y,yelmo1%tpo%now%H_ice, &
                            experiment="dome",time=ts%time,smb_max=0.5_wp,rad_el=1200.0_wp,period=period,dT_test=dT_test)
                
        case("mismip-stnd") 

            ! To do: define calv_mask!!

            ! Initialize mismip boundary values 
            call mismip3D_boundaries(yelmo1%bnd%T_srf,yelmo1%bnd%smb,yelmo1%bnd%Q_geo,yelmo1%bnd%calv_mask,yelmo1%bnd%c%T0,experiment="Stnd") 


        case("dome") 
            ! Boundary conditions are free to be chosen here 

            ! Set conditions similar to EISMINT2-EXPA with smaller radius 
            call dome_boundaries(yelmo1%bnd%T_srf,yelmo1%bnd%smb,yelmo1%bnd%Q_geo, &
                            yelmo1%grd%x,yelmo1%grd%y,yelmo1%tpo%now%H_ice, &
                            experiment="dome",time=ts%time,smb_max=0.3_wp,rad_el=300.0_wp,period=period,dT_test=dT_test)
                
        case DEFAULT 
            ! EISMINT 

            ! Initialize eismint boundary values 
            call eismint_boundaries(yelmo1%bnd%T_srf,yelmo1%bnd%smb,yelmo1%bnd%Q_geo, &
                                    yelmo1%grd%x,yelmo1%grd%y,yelmo1%tpo%now%H_ice, &
                                    experiment=experiment,time=0.0_wp,period=period,dT_test=dT_test)


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
    if (.FALSE.) then
        ! Set conditions similar to EISMINT2-EXPA with smaller radius 
        call dome_init(yelmo1%tpo%now%H_ice,yelmo1%grd%x,yelmo1%grd%y,R0=0.5_wp,H0=1000.0_wp,H0_shlf=200.0_wp,rmax_shlf=0.6_wp)
        yelmo1%tpo%now%z_srf  = yelmo1%bnd%z_bed + yelmo1%tpo%now%H_ice
    end if 

    ! Check boundary values 
    call yelmo_print_bound(yelmo1%bnd)

    ! Initialize state variables (dyn,therm,mat)
    call yelmo_init_state(yelmo1,time=ts%time,thrm_method="robin")

    ! For dome experiment, let it equilibrate for several thousand years
    if (trim(experiment) .eq. "dome") then 

        ! Set conditions similar to EISMINT2-EXPA with smaller radius 
        call dome_boundaries(yelmo1%bnd%T_srf,yelmo1%bnd%smb,yelmo1%bnd%Q_geo, &
                        yelmo1%grd%x,yelmo1%grd%y,yelmo1%tpo%now%H_ice, &
                        experiment="dome",time=ts%time,smb_max=0.3_wp,rad_el=300.0_wp,period=period,dT_test=dT_test)

        call yelmo_update_equil(yelmo1,ts%time,time_tot=real(5e3,wp), &
                                                dt=5.0_wp,topo_fixed=.FALSE.)  
    end if 

    ! == Write initial state ==
     
    ! 2D file 
    call yelmo_write_init(yelmo1,file2D,time_init=ts%time,units="years")
    call write_step_2D(yelmo1,file2D,time=ts%time)  
    !call yelmo_write_step(yelmo1,file2D,ts%time) 
            
    ! 1D file 
    call yelmo_write_reg_init(yelmo1,file1D,time_init=ts%time,units="years",mask=yelmo1%bnd%ice_allowed)
    call yelmo_write_reg_step(yelmo1,file1D,time=ts%time) 

    ! Comparison file 
    call yelmo_write_init(yelmo1,file_compare,time_init=ts%time,units="years")
    call write_step_2D_bueler(yelmo1,buel,file_compare,ts%time)
    
    ! Store default bedrock solver method, to be activated after several years 
    thrm_method_default = trim(yelmo1%thrm%par%method)
    rock_method_default = trim(yelmo1%thrm%par%rock_method)

    ! Advance timesteps
    do while (.not. ts%is_finished)

        ! == Update timestep ===

        call tstep_update(ts,dtt,verbose=.FALSE.)

        ! Update bnd%enh_srf to test transition of enhancement layers in time 
        if (ts%time_elapsed .lt. 0.5*(ts%time_end-ts%time_init)) then 
            yelmo1%bnd%enh_srf = 3.0 
        else 
            yelmo1%bnd%enh_srf = 1.0 
        end if 

        ! == Yelmo ice sheet ===================================================
        call yelmo_update(yelmo1,ts%time)
        
        ! == Update boundaries 
        select case(trim(experiment))

            case("BUELER-A")

                ! Pass - no update to perform here 

                
            case("BUELER-B")

                ! Update BUELER-B profile
                call bueler_test_BC(buel%H_ice,buel%mbal,buel%u_b,yelmo1%grd%x,yelmo1%grd%y, &
                      time=ts%time,R0=750.0_wp,H0=3600.0_wp,lambda=0.0_wp,n=3.0_wp,A=1e-16_wp, &
                      rho_ice=yelmo1%bnd%c%rho_ice,g=yelmo1%bnd%c%g)

            case("HALFAR")

                ! Update BUELER-B (but with HALFAR conditions)
                call bueler_test_BC(buel%H_ice,buel%mbal,buel%u_b,yelmo1%grd%x,yelmo1%grd%y, &
                      time=ts%time,R0=21.2132_wp,H0=707.1_wp,lambda=0.0_wp,n=3.0_wp,A=1e-16_wp, &
                      rho_ice=yelmo1%bnd%c%rho_ice,g=yelmo1%bnd%c%g)

            case("HALFAR-MED")

                ! Initialize BUELER-B (but with conditions between Bueler-B and HALFAR - not too big, not too small)
                call bueler_test_BC(buel%H_ice,buel%mbal,buel%u_b,yelmo1%grd%x,yelmo1%grd%y, &
                            time=ts%time,R0=100.0_wp,H0=800.0_wp,lambda=0.0_wp,n=3.0_wp,A=1e-16_wp, &
                            rho_ice=yelmo1%bnd%c%rho_ice,g=yelmo1%bnd%c%g)
            
            case("mismip")

                ! Set conditions similar to EISMINT2-EXPA with smaller radius 
                call dome_boundaries(yelmo1%bnd%T_srf,yelmo1%bnd%smb,yelmo1%bnd%Q_geo, &
                                yelmo1%grd%x,yelmo1%grd%y,yelmo1%tpo%now%H_ice, &
                                experiment="dome",time=ts%time,smb_max=0.5_wp,rad_el=1200.0_wp,period=period,dT_test=dT_test)
                
            case("mismip-stnd") 

                ! Initialize mismip boundary values 
                call mismip3D_boundaries(yelmo1%bnd%T_srf,yelmo1%bnd%smb,yelmo1%bnd%Q_geo,yelmo1%bnd%calv_mask,yelmo1%bnd%c%T0,experiment="Stnd") 

            case("dome") 
            ! Boundary conditions are free to be chosen here 

                ! Set conditions similar to EISMINT2-EXPA with smaller radius 
                call dome_boundaries(yelmo1%bnd%T_srf,yelmo1%bnd%smb,yelmo1%bnd%Q_geo, &
                                yelmo1%grd%x,yelmo1%grd%y,yelmo1%tpo%now%H_ice, &
                                experiment="dome",time=ts%time,smb_max=0.3_wp,rad_el=300.0_wp,period=period,dT_test=dT_test)

            case DEFAULT 

                call eismint_boundaries(yelmo1%bnd%T_srf,yelmo1%bnd%smb,yelmo1%bnd%Q_geo, &
                                yelmo1%grd%x,yelmo1%grd%y,yelmo1%tpo%now%H_ice, &
                                experiment=experiment,time=ts%time,period=period,dT_test=dT_test)

        end select 

        ! Make comparison if relevant 
        select case(trim(experiment))
            case("BUELER-A","BUELER-B","HALFAR","HALFAR-MED")

                call bueler_compare(buel,yelmo1%tpo%now%H_ice,dx=yelmo1%grd%dx)
                
            case DEFAULT 
            ! Pass- not a bueler test... 
        end select 

        ! == MODEL OUTPUT =======================================================
        if (mod(nint(ts%time*100),nint(dt2D_out*100))==0) then 
            call write_step_2D(yelmo1,file2D,time=ts%time) 
            !call yelmo_write_step(yelmo1,file2D,time) 

            select case(trim(experiment))
                case("BUELER-A","BUELER-B")
                    call write_step_2D_bueler(yelmo1,buel,file_compare,ts%time)
            end select
            
        end if 

        if (mod(nint(ts%time*100),nint(dt1D_out*100))==0) then 
            call yelmo_write_reg_step(yelmo1,file1D,time=ts%time) 
        end if 

        if (mod(ts%time,10.0)==0 .and. (.not. yelmo_log)) then
            write(*,"(a,f14.4)") "yelmo:: time = ", ts%time
        end if 

    end do 




    ! Write summary 
    write(*,*) "====== "//trim(domain)//"-"//trim(experiment)//" ======="
    write(*,*) "nz_aa, H0 = ", yelmo1%par%nz_aa, maxval(yelmo1%tpo%now%H_ice)

    ! Make comparison if relevant 
    select case(trim(experiment))
        case("BUELER-A","BUELER-B","HALFAR","HALFAR-MED")

            write(*,"(a,3f10.2,10g12.3)") trim(experiment), ts%time, dtt, yelmo1%grd%dx*1e-3, &
                                            buel%rmse_H_ice, buel%err_H0, buel%err_max_H_ice, buel%err_V_ice

        case DEFAULT 
        ! Pass- not a bueler test... 
    end select 


    ! Write a restart file too
    call yelmo_restart_write(yelmo1,file_restart,time=ts%time)

    ! Finalize program
    call yelmo_end(yelmo1,time=ts%time)

    ! Stop timing 
    call yelmo_cpu_time(cpu_end_time,cpu_start_time,cpu_dtime)
    
    write(*,"(a,f12.3,a)") "Time  = ",cpu_dtime/60.0 ," min"
    write(*,"(a,f12.1,a)") "Speed = ",(1e-3*(ts%time_end-ts%time_init))/(cpu_dtime/3600.0), " kiloyears / hr"
    
contains
    
    subroutine write_step_2D(ylmo,filename,time,buel)

        implicit none 
        
        type(yelmo_class), intent(IN) :: ylmo
        character(len=*),  intent(IN) :: filename
        real(wp),          intent(IN) :: time
        type(bueler_test_type), intent(IN), optional :: buel 

        ! Local variables
        integer    :: ncid, n, i, j, nx, ny   
        real(wp), allocatable :: sym(:,:) 

        integer :: q, qtot
        character(len=32), allocatable :: vnms(:)

        nx = ylmo%tpo%par%nx 
        ny = ylmo%tpo%par%ny 

        allocate(sym(nx,ny)) 

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_time_index(filename,"time",time,ncid)

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        ! Write model metrics (model speed, dt, eta)
        call yelmo_write_step_model_metrics(filename,ylmo,n,ncid)

        ! Time step limits 
        call nc_write(filename,"dt_adv",ylmo%time%dt_adv,units="a",long_name="Advective timestep", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dt_diff",ylmo%time%dt_diff,units="a",long_name="Diffusive timestep", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"dt_adv3D",ylmo%time%dt_adv3D,units="a",long_name="Advective timestep", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        qtot = 19
        allocate(vnms(qtot))
        vnms(1) = "H_ice"
        vnms(2) = "z_srf"
        vnms(3) = "mask_bed"
        vnms(4) = "H_grnd"
        vnms(5) = "mb_net"
        vnms(6) = "mb_resid"
        vnms(7) = "dzsdt"
        vnms(8) = "dHidt"
        vnms(9) = "dzsdx"
        vnms(10) = "dzsdy"
        vnms(11) = "dHidx"
        vnms(12) = "dHidy"
        vnms(13) = "f_grnd"
        vnms(14) = "f_grnd_acx"
        vnms(15) = "f_grnd_acy"
        vnms(16) = "f_ice"
        vnms(17) = "f_ice_dyn"
        vnms(18) = "N_eff"
        vnms(19) = "cmb"

        do q = 1, qtot
            call yelmo_write_var(filename,vnms(q),ylmo,n,ncid)
        end do
        
        ! == yelmo_thermodynamics ==
        call nc_write(filename,"enth",ylmo%thrm%now%enth,units="J m-3",long_name="Ice enthalpy", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"T_ice",ylmo%thrm%now%T_ice,units="K",long_name="Ice temperature", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"omega",ylmo%thrm%now%omega,units="--",long_name="Ice water content", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"T_prime",ylmo%thrm%now%T_prime,units="deg C",long_name="Homologous ice temperature", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
!         call nc_write(filename,"T_pmp",ylmo%thrm%now%T_pmp,units="K",long_name="Ice pressure melting point (pmp)", &
!                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        
        call nc_write(filename,"uz_star",ylmo%dyn%now%uz_star,units="m yr-1",long_name="Advection-adjusted vertical velocity", &
                      dim1="xc",dim2="yc",dim3="zeta_ac",dim4="time",start=[1,1,1,n],ncid=ncid)
        
        call nc_write(filename,"T_rock",ylmo%thrm%now%T_rock,units="K",long_name="Bedrock temperature", &
                      dim1="xc",dim2="yc",dim3="zeta_rock",dim4="time",start=[1,1,1,n],ncid=ncid)
        
        call nc_write(filename,"Q_rock",ylmo%thrm%now%Q_rock,units="mW m-2",long_name="Bedrock surface heat flux", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"T_prime_b",ylmo%thrm%now%T_prime_b,units="deg C",long_name="Homologous basal ice temperature", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        sym = ylmo%thrm%now%T_prime_b(nx:1:-1,:) - ylmo%thrm%now%T_prime_b
        call nc_write(filename,"T_prime_b_sym",sym,units="deg C",long_name="Homologous basal ice temperature symmetry check", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"f_pmp",ylmo%thrm%now%f_pmp,units="1",long_name="Fraction of grid point at pmp", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"Q_ice_b",ylmo%thrm%now%Q_ice_b,units="mW m-2",long_name="Basal ice heat flux", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"Q_strn",ylmo%thrm%now%Q_strn/(ylmo%bnd%c%rho_ice*ylmo%thrm%now%cp),units="K yr-1",long_name="Strain heating", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        call nc_write(filename,"dQsdT",ylmo%thrm%now%dQsdT,units="yr-1",long_name="Strain heating derivative w.r.t. temperature", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        call nc_write(filename,"Q_b",ylmo%thrm%now%Q_b,units="mW m-2",long_name="Basal frictional heating", &
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
        
        call nc_write(filename,"dep_time",ylmo%mat%now%dep_time,units="a",long_name="Ice deposition time", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
!         call nc_write(filename,"age",time-ylmo%mat%now%dep_time,units="a",long_name="Ice age", &
!                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        call nc_write(filename,"enh_bnd",ylmo%mat%now%enh_bnd,units="",long_name="Enhancement factor tracer field", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"enh",ylmo%mat%now%enh,units="",long_name="Enhancement factor", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        ! Strain-rate and stress tensors 
        if (.FALSE.) then

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
            
            sym = ylmo%mat%now%strn2D%de(nx:1:-1,:) - ylmo%mat%now%strn2D%de
            call nc_write(filename,"de2D_sym",sym,units="yr^-1",long_name="Effective strain rate symmetry resid.", &
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

        ! == yelmo_dynamics ==

        call nc_write(filename,"ssa_mask_acx",ylmo%dyn%now%ssa_mask_acx,units="1",long_name="SSA mask (acx)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"ssa_mask_acy",ylmo%dyn%now%ssa_mask_acy,units="1",long_name="SSA mask (acy)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"cb_ref",ylmo%dyn%now%cb_ref,units="--",long_name="Bed friction scalar", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"c_bed",ylmo%dyn%now%c_bed,units="Pa",long_name="Bed friction coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"beta",ylmo%dyn%now%beta,units="Pa a m-1",long_name="Basal friction coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"beta_acx",ylmo%dyn%now%beta_acx,units="Pa a m-1",long_name="Basal friction coefficient (acx)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"beta_acy",ylmo%dyn%now%beta_acy,units="Pa a m-1",long_name="Basal friction coefficient (acy)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"visc_eff",ylmo%dyn%now%visc_eff,units="Pa a m",long_name="Effective viscosity (SSA)", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        call nc_write(filename,"visc_eff_int",ylmo%dyn%now%visc_eff_int,units="Pa a m",long_name="Depth-integrated effective viscosity (SSA)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"beta_eff",ylmo%dyn%now%beta_eff,units="Pa a m-1",long_name="Effective basal friction coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"duxdz",ylmo%dyn%now%duxdz,units="1/a",long_name="Vertical shear (x)", &
                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"duydz",ylmo%dyn%now%duydz,units="1/a",long_name="Vertical shear (y)", &
                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        
        ! call nc_write(filename,"dxx",ylmo%dyn%now%strn%dxx,units="1/a",long_name="Strain dxx", &
        !                dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        ! call nc_write(filename,"dxy",ylmo%dyn%now%strn%dxy,units="1/a",long_name="Strain dxy", &
        !                dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        ! call nc_write(filename,"dyy",ylmo%dyn%now%strn%dyy,units="1/a",long_name="Strain dyy", &
        !                dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        ! call nc_write(filename,"dxz",ylmo%dyn%now%strn%dxz,units="1/a",long_name="Strain dxz", &
        !                dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        ! call nc_write(filename,"dyz",ylmo%dyn%now%strn%dyz,units="1/a",long_name="Strain dyz", &
        !                dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        
!         call nc_write(filename,"ux_i_bar",ylmo%dyn%now%ux_i_bar,units="m/a",long_name="Internal shear velocity (x)", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"uy_i_bar",ylmo%dyn%now%uy_i_bar,units="m/a",long_name="Internal shear velocity (y)", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_i_bar",ylmo%dyn%now%uxy_i_bar,units="m/a",long_name="Internal shear velocity magnitude", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

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

        call nc_write(filename,"qq_acx",ylmo%dyn%now%qq_acx,units="m^3/yr",long_name="Ice flux (acx-nodes)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"qq_acy",ylmo%dyn%now%qq_acy,units="m^3/yr",long_name="Ice flux (acy-nodes)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"qq",ylmo%dyn%now%qq,units="m^3/yr",long_name="Ice flux", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"taud_acx",ylmo%dyn%now%taud_acx,units="Pa",long_name="Driving stress (x)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taud_acy",ylmo%dyn%now%taud_acy,units="Pa",long_name="Driving stress (y)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taud",ylmo%dyn%now%taud,units="Pa",long_name="Driving stress", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"taub_acx",ylmo%dyn%now%taub_acx,units="Pa",long_name="Basal stress (x)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taub_acy",ylmo%dyn%now%taub_acy,units="Pa",long_name="Basal stress (y)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taub",ylmo%dyn%now%taub,units="Pa",long_name="Basal stress", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"ux",ylmo%dyn%now%ux,units="m/yr",long_name="Horizontal velocity (x)", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"uy",ylmo%dyn%now%uy,units="m/yr",long_name="Horizontal velocity (y)", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"uxy",ylmo%dyn%now%uxy,units="m/yr",long_name="Horizontal velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"uz",ylmo%dyn%now%uz,units="m/yr",long_name="Vertical velocity", &
                      dim1="xc",dim2="yc",dim3="zeta_ac",dim4="time",start=[1,1,1,n],ncid=ncid)

        call nc_write(filename,"f_vbvs",ylmo%dyn%now%f_vbvs,units="1",long_name="Basal to surface velocity fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_shear_bar",ylmo%mat%now%f_shear_bar,units="1",long_name="Vertically averaged shearing fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"de",ylmo%mat%now%strn%de,units="yr^-1",long_name="Strain rate", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        ! == yelmo_bound ==

        call nc_write(filename,"z_bed",ylmo%bnd%z_bed,units="m",long_name="Bedrock elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
!         call nc_write(filename,"z_sl",ylmo%bnd%z_sl,units="m",long_name="Sea level rel. to present", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"H_sed",ylmo%bnd%H_sed,units="m",long_name="Sediment thickness", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"smb",ylmo%tpo%now%smb,units="m/yr ice equiv.",long_name="Net surface mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"smb_ref",ylmo%bnd%smb,units="m/yr ice equiv.",long_name="Surface mass balance (potential)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"T_srf",ylmo%bnd%T_srf,units="K",long_name="Surface temperature", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
     
        call nc_write(filename,"bmb",ylmo%tpo%now%bmb,units="m/yr ice equiv.",long_name="Basal mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"bmb_shlf",ylmo%bnd%bmb_shlf,units="m/yr ice equiv.",long_name="Basal mass balance (shelf)", &
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
        real(wp), intent(IN) :: time
        
        ! Local variables
        integer    :: ncid, n, i, j, nx, ny  

        nx = ylmo%tpo%par%nx 
        ny = ylmo%tpo%par%ny 

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_time_index(filename,"time",time,ncid)

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
