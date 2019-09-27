

program yelmo_mismip

    use nml 
    use ncio  
    use yelmo 

    use deformation 

    use ice_benchmarks 
    use mismip3D 

    implicit none 

    type(yelmo_class)     :: yelmo1

    character(len=56)  :: domain 
    character(len=256) :: outfldr, file2D, file1D, file_compare
    character(len=512) :: path_par, path_const 
    character(len=56)  :: experiment
    logical            :: with_ssa  
    real(prec) :: time_init, time_end, time, dtt, dt2D_out, dt1D_out
    real(prec) :: period, dt_test 
    integer    :: n  

    character(len=56) :: grid_name
    real(prec) :: dx 
    integer    :: nx  

    integer :: n_att, n_att_tot, q_att, q
    real(prec), allocatable :: ATT_values(:)
    real(prec) :: ATT_time, ATT_dt 
    logical    :: is_converged, exit_loop 
    real(prec) :: err  
    real(prec) :: x_gl, x_gl_std 

    real(4) :: start, finish
    
    ! Start timing 
    call cpu_time(start)

    ! Assume program is running from the output folder
    outfldr = "./"

    ! Determine the parameter file from the command line 
    call yelmo_load_command_line_args(path_par)
    path_const = trim(outfldr)//"yelmo_const_EISMINT.nml"
    
    ! Define input and output locations 
    file1D       = trim(outfldr)//"yelmo1D.nc"
    file2D       = trim(outfldr)//"yelmo2D.nc"
    file_compare = trim(outfldr)//"yelmo_compare.nc"
    
    ! Define the domain, grid and experiment from parameter file
    call nml_read(path_par,"eismint","domain",       domain)        ! EISMINT1, EISMINT2
    call nml_read(path_par,"eismint","experiment",   experiment)    ! "linear", "overdeepened", "flat"
    call nml_read(path_par,"eismint","dx",           dx)            ! [km] Grid resolution 
    call nml_read(path_par,"eismint","with_ssa",     with_ssa)      ! Include ssa in experiment?

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
    
    ! Square domain to handle mismip sloping beds radially

    grid_name = "EISMINT1-mismip"
    nx = 73

    ! === Initialize ice sheet model =====

    ! General initialization of yelmo constants (used globally)
    call yelmo_global_init(path_const)

    ! Next define grid 
    call yelmo_init_grid(yelmo1%grd,grid_name,units="km",dx=dx,nx=nx,dy=dx,ny=nx)

    ! Initialize data objects (without loading topography, which will be defined inline below)
    call yelmo_init(yelmo1,filename=path_par,grid_def="none",time=time_init,load_topo=.FALSE.,domain=domain,grid_name=grid_name)
    
    ! Update parameter values with choices 
    yelmo1%dyn%par%use_ssa    = with_ssa 

    ! === Define initial topography =====

    select case(trim(experiment))

        case("linear")
            ! Add a marine bed for testing ice shelves following MISMIP
        
            yelmo1%bnd%z_bed  = 720.0 - 778.50*(sqrt((yelmo1%grd%x*1e-3)**2+(yelmo1%grd%y*1e-3)**2))/750.0
            yelmo1%tpo%now%H_ice  = 100.0
            yelmo1%tpo%now%z_srf  = yelmo1%bnd%z_bed + yelmo1%tpo%now%H_ice

            ! Limit ice to circular domain for symmetry improvements 
            yelmo1%bnd%ice_allowed = .TRUE. 
            where(yelmo1%bnd%z_bed .lt. -1000.0) yelmo1%bnd%ice_allowed = .FALSE. 

        case("overdeepened")
            ! Add a marine bed for testing ice shelves following MISMIP
        
            yelmo1%bnd%z_bed  = 729.0 -  (2184.8/750.0**2)*(sqrt((yelmo1%grd%x*1e-3)**2+(yelmo1%grd%y*1e-3)**2))**2 &
                                      + (1031.72/750.0**4)*(sqrt((yelmo1%grd%x*1e-3)**2+(yelmo1%grd%y*1e-3)**2))**4 &
                                      -  (151.72/750.0**6)*(sqrt((yelmo1%grd%x*1e-3)**2+(yelmo1%grd%y*1e-3)**2))**6
            yelmo1%tpo%now%H_ice  = 100.0
            yelmo1%tpo%now%z_srf  = yelmo1%bnd%z_bed + yelmo1%tpo%now%H_ice

        case DEFAULT 
            ! Flat bed, useful for diagnostics

            yelmo1%bnd%z_bed      = 0.0 
            yelmo1%tpo%now%H_ice  = 0.0
            yelmo1%tpo%now%z_srf  = yelmo1%bnd%z_bed + yelmo1%tpo%now%H_ice 
            
    end select 

    ! Load boundary values

    yelmo1%bnd%z_sl     = 0.0
    yelmo1%bnd%bmb_shlf = 0.0  
    yelmo1%bnd%T_shlf   = T0  
    yelmo1%bnd%H_sed    = 0.0 
    yelmo1%bnd%H_w      = 0.0

    select case(trim(experiment))

        case("linear","overdeepened") 
            
            ! Initialize mismip boundary values 
            call mismip3D_boundaries(yelmo1%bnd%T_srf,yelmo1%bnd%smb,yelmo1%bnd%Q_geo,experiment="Stnd-0.3") 

        case DEFAULT 
            ! EISMINT - just for testing

            ! Initialize eismint boundary values 
            call eismint_boundaries(yelmo1%bnd%T_srf,yelmo1%bnd%smb,yelmo1%bnd%Q_geo, &
                                    yelmo1%grd%x,yelmo1%grd%y,yelmo1%tpo%now%H_ice, &
                                    experiment="mismip",time=0.0_prec,period=period,dT_test=dT_test)

    end select 


    ! Define rate factor experiment ============================================================
    
!     ! Longer experiment 
!     n_att     = 17
!     allocate(ATT_values(n_att))
!     ATT_values = sec_year*1e-26*[464.16,215.44,100.00,46.416,21.544,10.0,4.6416,2.1544,1.0, &
!                                  2.1544,4.6416,10.0,21.544,46.416,100.00,215.44,464.16]
    
!     ! Shorter experiment (Pattyn, 2017)
!     n_att      = 7
!     allocate(ATT_values(n_att))
!     ATT_values = [1e-16,1e-17,1e-18,1e-19,1e-18,1e-17,1e-16]
    
    ! Shorter experiment 2 (Pattyn, 2017)
    n_att      = 3
    allocate(ATT_values(n_att))
    ATT_values = [1e-16,1e-17,1e-16]
    
    ATT_time   = 15e3
    ATT_dt     =  2e3 
    time_end   = ATT_time + n_att*ATT_dt + 50e3
    dt2D_out   = 500.0 
    
    write(*,*) "time_init = ", time_init 
    write(*,*) "time_end  = ", time_end
    do q_att = 1, n_att 
        write(*,"(i3,f10.3,g15.5)") q_att, ATT_time+(q_att-1)*ATT_dt, ATT_values(q_att)
    end do
        
    ! Ensure rate factor (and counter) is set for first value if performing RF experiment 
    q_att = 1 
    yelmo1%mat%par%rf_const = ATT_values(q_att)

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

    ! Calculate grounding line position 
    call find_x_gl_2D(x_gl,x_gl_std,yelmo1%grd%x*1e-3,yelmo1%grd%y*1e-3,yelmo1%tpo%now%f_grnd)

    ! == Write initial state ==
     
    ! 2D file 
    call yelmo_write_init(yelmo1,file2D,time_init=time_init,units="years")
    call write_step_2D(yelmo1,file2D,time=time_init,x_gl=x_gl,x_gl_std=x_gl_std)  
    
    ! 1D file 
    call write_yreg_init(yelmo1,file1D,time_init=time_init,units="years",mask=yelmo1%bnd%ice_allowed)
    call write_yreg_step(yelmo1%reg,file1D,time=time_init) 


    ! Set exit to false
    exit_loop   = .FALSE. 

    ! Advance timesteps
    do n = 1, ceiling((time_end-time_init)/dtt)

        ! Get current time 
        time = time_init + n*dtt

        ! Check for convergence 
        is_converged = .FALSE. 
        err = sqrt(sum(yelmo1%tpo%now%dHicedt**2)/yelmo1%grd%npts)
        if (err .lt. 1e-2) is_converged =.TRUE. 

        ! Modify rate factor 
        if (time .gt. ATT_time+ATT_dt) then 
            ! Ensure minimum time per step has been reached before checking convergence

            !write(*,*) "err: ", time, ATT_time, err, yelmo1%mat%par%rf_const, q_att 
        
            if (is_converged .and. q_att == n_att) then 
                ! If output timestep also reached,
                ! then time to kill simulation 
                if (mod(time,dt2D_out)==0) exit_loop = .TRUE. 
            else if (is_converged) then
                ! Time to step ATT_value 
                q_att = min(q_att+1,n_att)
                yelmo1%mat%par%rf_const = ATT_values(q_att)
                ATT_time = time 
                dt2D_out = 500.0
            end if   

        end if 


        ! == Yelmo ice sheet ===================================================
        call yelmo_update(yelmo1,time)
        
        ! == Update boundaries 
        select case(trim(experiment))

            case("linear","overdeepened") 
            
                ! Initialize mismip boundary values 
                call mismip3D_boundaries(yelmo1%bnd%T_srf,yelmo1%bnd%smb,yelmo1%bnd%Q_geo,experiment="Stnd-0.3")

            case DEFAULT 

                call eismint_boundaries(yelmo1%bnd%T_srf,yelmo1%bnd%smb,yelmo1%bnd%Q_geo, &
                                yelmo1%grd%x,yelmo1%grd%y,yelmo1%tpo%now%H_ice, &
                                experiment=experiment,time=time,period=period,dT_test=dT_test)

        end select 


        ! Calculate grounding line position 
        call find_x_gl_2D(x_gl,x_gl_std,yelmo1%grd%x*1e-3,yelmo1%grd%y*1e-3,yelmo1%tpo%now%f_grnd)

        ! == MODEL OUTPUT =======================================================
        if (mod(nint(time*100),nint(dt2D_out*100))==0) then 
            call write_step_2D(yelmo1,file2D,time=time,x_gl=x_gl,x_gl_std=x_gl_std)  
        end if 

        if (mod(nint(time*100),nint(dt1D_out*100))==0) then 
            call write_yreg_step(yelmo1%reg,file1D,time=time) 
        end if 

        write(*,"(a,f14.4)") "time = ", time

        if (mod(nint(time*100),nint((5.0*dtt)*100))==0) then
            write(*,"(a,2f14.1,g14.3,f14.1,f10.1)") "mism: ",  &
                time, maxval(yelmo1%tpo%now%H_ice), yelmo1%mat%par%rf_const, x_gl, x_gl_std 
        end if 

        if (exit_loop) exit 

    end do 

    ! Write summary 
    write(*,*) "====== "//trim(domain)//"-"//trim(experiment)//" ======="
    write(*,*) "nz_aa, H0 = ", yelmo1%par%nz_aa, maxval(yelmo1%tpo%now%H_ice)

    ! Finalize program
    call yelmo_end(yelmo1,time=time)

    ! Stop timing 
    call cpu_time(finish)

    print '("Time = ",f12.3," min.")', (finish-start)/60.0 

contains

    subroutine write_step_2D(ylmo,filename,time,x_gl,x_gl_std)

        implicit none 
        
        type(yelmo_class), intent(IN) :: ylmo
        character(len=*),  intent(IN) :: filename
        real(prec),        intent(IN) :: time
        real(prec),        intent(IN) :: x_gl 
        real(prec),        intent(IN) :: x_gl_std 

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

        ! Write model speed 
        call nc_write(filename,"speed",ylmo%par%model_speed,units="kyr/hr",long_name="Model speed (Yelmo only)", &
                      dim1="time",start=[n],count=[1],ncid=ncid)
        
        ! 1D variables of interest 
        call nc_write(filename,"x_rf",ylmo%mat%par%rf_const,units="a^-1 Pa^-3", &
            long_name="Rate factor",dim1="time",start=[n],count=[1],ncid=ncid)
        call nc_write(filename,"x_gl",x_gl,units="km", &
            long_name="Grounding line position",dim1="time",start=[n],count=[1],ncid=ncid)
        call nc_write(filename,"x_gl_std",x_gl_std,units="km", &
            long_name="Grounding line position Stdev",dim1="time",start=[n],count=[1],ncid=ncid)
            
        ! Time step limits 
        call nc_write(filename,"dt_adv",ylmo%par%dt_adv,units="a",long_name="Advective timestep", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dt_diff",ylmo%par%dt_diff,units="a",long_name="Diffusive timestep", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! == yelmo_topography ==
        call nc_write(filename,"H_ice",ylmo%tpo%now%H_ice,units="m",long_name="Ice thickness", &
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
!         call nc_write(filename,"calv",ylmo%tpo%now%calv,units="m/a ice equiv.",long_name="Calving rate", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! == yelmo_thermodynamics ==
        call nc_write(filename,"T_ice",ylmo%thrm%now%T_ice,units="K",long_name="Ice temperature", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"T_prime",ylmo%thrm%now%T_ice-ylmo%thrm%now%T_pmp,units="deg C",long_name="Homologous ice temperature", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
!         call nc_write(filename,"T_pmp",ylmo%thrm%now%T_pmp,units="K",long_name="Ice pressure melting point (pmp)", &
!                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        
!         call nc_write(filename,"enth_ice",ylmo%thrm%now%enth_ice,units="J/kg",long_name="Ice enthalpy", &
!                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
!         call nc_write(filename,"omega_ice",ylmo%thrm%now%omega_ice,units="%",long_name="Ice water content", &
!                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        
        call nc_write(filename,"f_pmp",ylmo%thrm%now%f_pmp,units="1",long_name="Fraction of grid point at pmp", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"Q_strn",ylmo%thrm%now%Q_strn/(rho_ice*ylmo%thrm%now%cp),units="K a-1",long_name="Strain heating", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"Q_b",ylmo%thrm%now%Q_b,units="J a-1 m-2",long_name="Basal frictional heating", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"dTdz_b",ylmo%thrm%now%dTdz_b,units="K/m",long_name="Basal temperature gradient (ice)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"dTrdz_b",ylmo%thrm%now%dTrdz_b,units="K/m",long_name="Surface temperature gradient (rock)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"bmb_grnd",ylmo%thrm%now%bmb_grnd,units="m/a ice equiv.",long_name="Basal mass balance (grounded)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! == yelmo_material ==
!         call nc_write(filename,"enh_bar",ylmo%mat%now%enh_bar,units="1",long_name="Vertically averaged enhancement factor", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"visc_int",ylmo%mat%now%visc_int,units="Pa a m",long_name="Vertically integrated viscosity", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"ATT",ylmo%mat%now%ATT,units="a^-1 Pa^-3",long_name="Rate factor", &
!                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"visc",ylmo%mat%now%visc,units="Pa a",long_name="Viscosity", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        
!         call nc_write(filename,"dep_time",ylmo%mat%now%dep_time,units="a",long_name="Ice deposition time", &
!                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
!         call nc_write(filename,"age",time-ylmo%mat%now%dep_time,units="a",long_name="Ice age", &
!                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        ! == yelmo_dynamics ==

        call nc_write(filename,"taud_acx",ylmo%dyn%now%taud_acx,units="Pa",long_name="Driving stress, x-direction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taud_acy",ylmo%dyn%now%taud_acy,units="Pa",long_name="Driving stress, y-direction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"taub_acx",ylmo%dyn%now%taub_acx,units="Pa",long_name="Basal stress, x-direction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taub_acy",ylmo%dyn%now%taub_acy,units="Pa",long_name="Basal stress, y-direction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"ssa_mask_acx",ylmo%dyn%now%ssa_mask_acx,units="1",long_name="SSA mask (acx)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"ssa_mask_acy",ylmo%dyn%now%ssa_mask_acy,units="1",long_name="SSA mask (acy)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"C_bed",ylmo%dyn%now%C_bed,units="Pa",long_name="Bed friction coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"beta",ylmo%dyn%now%beta,units="Pa a m-1",long_name="Basal friction coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"visc_eff",ylmo%dyn%now%visc_eff,units="Pa a m",long_name="Effective viscosity (SSA)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"sigma_horiz_sq",ylmo%dyn%now%sigma_horiz_sq,units="1",long_name="Horizontal stress components squared", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"lhs_x",ylmo%dyn%now%lhs_x,units="Pa",long_name="Shear reduction (x)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"lhs_y",ylmo%dyn%now%lhs_y,units="Pa",long_name="Shear reduction (y)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"lhs_xy",ylmo%dyn%now%lhs_xy,units="Pa",long_name="Shear reduction magnitude", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
!         call nc_write(filename,"duxdz",ylmo%dyn%now%duxdz,units="1/a",long_name="Vertical shear (x)", &
!                        dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
!         call nc_write(filename,"duydz",ylmo%dyn%now%duydz,units="1/a",long_name="Vertical shear (y)", &
!                        dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        
!         call nc_write(filename,"ux_i_bar",ylmo%dyn%now%ux_i_bar,units="m/a",long_name="Internal shear velocity (x)", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"uy_i_bar",ylmo%dyn%now%uy_i_bar,units="m/a",long_name="Internal shear velocity (y)", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"uxy_i_bar",ylmo%dyn%now%uxy_i_bar,units="m/a",long_name="Internal shear velocity magnitude", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"ux_b",ylmo%dyn%now%ux_b,units="m/a",long_name="Basal sliding velocity (x)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"uy_b",ylmo%dyn%now%uy_b,units="m/a",long_name="Basal sliding velocity (y)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"uxy_b",ylmo%dyn%now%uxy_b,units="m/a",long_name="Basal sliding velocity magnitude", &
!                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"ux_bar",ylmo%dyn%now%ux_bar,units="m/a",long_name="Vertically integrated velocity (x)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy_bar",ylmo%dyn%now%uy_bar,units="m/a",long_name="Vertically integrated velocity (y)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_bar",ylmo%dyn%now%uxy_bar,units="m/a",long_name="Vertically integrated velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"qq_acx",ylmo%dyn%now%qq_acx,units="m^3/a",long_name="Ice flux (acx-nodes)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"qq_acy",ylmo%dyn%now%qq_acy,units="m^3/a",long_name="Ice flux (acy-nodes)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"qq",ylmo%dyn%now%qq,units="m^3/a",long_name="Ice flux", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
!         call nc_write(filename,"dd_ab",ylmo%dyn%now%dd_ab,units="m2 a-1",long_name="Diffusivity", &
!                        dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
!         call nc_write(filename,"dd_ab_bar",ylmo%dyn%now%dd_ab_bar,units="m2 a-1",long_name="Diffusivity", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
!         call nc_write(filename,"ux",ylmo%dyn%now%ux,units="m/a",long_name="Horizontal velocity (x)", &
!                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
!         call nc_write(filename,"uy",ylmo%dyn%now%uy,units="m/a",long_name="Horizontal velocity (y)", &
!                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
!         call nc_write(filename,"uxy",ylmo%dyn%now%uxy,units="m/a",long_name="Horizontal velocity magnitude", &
!                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
!         call nc_write(filename,"uz",ylmo%dyn%now%uz,units="m/a",long_name="Vertical velocity", &
!                       dim1="xc",dim2="yc",dim3="zeta_ac",dim4="time",start=[1,1,1,n],ncid=ncid)

!         call nc_write(filename,"f_vbvs",ylmo%dyn%now%f_vbvs,units="1",long_name="Basal to surface velocity fraction", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"f_shear_bar",ylmo%mat%now%f_shear_bar,units="1",long_name="Vertically averaged shearing fraction", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"de",ylmo%mat%now%strn%de,units="a^-1",long_name="Strain rate", &
!                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        ! == yelmo_bound ==

        call nc_write(filename,"z_bed",ylmo%bnd%z_bed,units="m",long_name="Bedrock elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
!         call nc_write(filename,"z_sl",ylmo%bnd%z_sl,units="m",long_name="Sea level rel. to present", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"H_sed",ylmo%bnd%H_sed,units="m",long_name="Sediment thickness", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"H_w",ylmo%bnd%H_w,units="m",long_name="Basal water pressure", &
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

end program yelmo_mismip 

