program yelmo_calving

    use nml 
    use ncio  
    use yelmo 
    use lsf_module

    use calving_benchmarks
    
    implicit none

    type(yelmo_class) :: yelmo1
    type(yelmo_class) :: yelmo_ref

    type control_type    
        character(len=256) :: outfldr
        character(len=256) :: file2D, file1D
        character(len=256) :: file_restart
        character(len=512) :: path_par 
        character(len=56)  :: exp

        real(wp) :: time_init, time_end, time, dtt
        real(wp) :: dt2D_out, dt1D_out

        real(wp) :: dx

        ! Internal parameters
        character(len=56)  :: domain
        character(len=56)  :: grid_name
        real(wp) :: x0, x1
        integer  :: nx
        integer  :: ny
    end type

    type(control_type) :: ctl
    real(wp) :: time 
    integer  :: n 
    
    real(8) :: cpu_start_time, cpu_end_time, cpu_dtime  
    
    ! Start timing 
    call yelmo_cpu_time(cpu_start_time)

    ! Assume program is running from the output folder
    ctl%outfldr = "./"

    ! Determine the parameter file from the command line 
    call yelmo_load_command_line_args(ctl%path_par)
    !path_par   = trim(outfldr)//"yelmo_calving.nml" 

    ! Define input and output locations 
    ctl%file1D       = trim(ctl%outfldr)//"yelmo1D.nc"
    ctl%file2D       = trim(ctl%outfldr)//"yelmo2D.nc"
    ctl%file_restart = trim(ctl%outfldr)//"yelmo_restart.nc"

    
    ! Define the domain, grid and experiment from parameter file
    call nml_read(ctl%path_par,"ctl","exp",         ctl%exp)            ! "exp1", "exp2", "exp3", "exp4", "exp5"
    call nml_read(ctl%path_par,"ctl","dx",          ctl%dx)             ! [km] Grid resolution 

    ! Timing parameters 
    call nml_read(ctl%path_par,"ctl","time_init",   ctl%time_init)      ! [yr] Starting time
    call nml_read(ctl%path_par,"ctl","time_end",    ctl%time_end)       ! [yr] Ending time
    call nml_read(ctl%path_par,"ctl","dtt",         ctl%dtt)            ! [yr] Main loop time step 
    call nml_read(ctl%path_par,"ctl","dt2D_out",    ctl%dt2D_out)       ! [yr] Frequency of 2D output 
    ctl%dt1D_out = ctl%dtt  ! Set 1D output to frequency of main loop timestep 

    ! Now set internal parameters ===

    ! Define domain and grid size based on experiment
    select case(trim(ctl%exp))
        case("exp1","exp2")
            ctl%domain = "circular"
            ctl%x0 = -900.0
            ctl%x1 =  900.0
        case("exp3","exp4","exp5")
            ctl%domain = "thule"
            ctl%x0 = -1000.0
            ctl%x1 =  1000.0
        case("advection")
            ctl%domain = "advection"
            ctl%x0     = 0.0
            ctl%x1     = 50.0
        case DEFAULT
            write(*,*) "ctl.exp = ",trim(ctl%domain), " not recognized."
            stop
    end select

    ! Get grid size
    ctl%nx = (ctl%x1-ctl%x0) / ctl%dx + 1
    ctl%ny = ctl%nx

    ! Get grid name
    write(ctl%grid_name,"(a,i2,a2)") trim(ctl%domain)//"-",int(ctl%dx),"KM"
   
    ! === Initialize ice sheet model =====

    ! First, define grid 
    call yelmo_init_grid(yelmo1%grd,ctl%grid_name,units="km",dx=ctl%dx,nx=ctl%nx,dy=ctl%dx,ny=ctl%nx)

    ! Initialize data objects (without loading topography, which will be defined inline below)
    call yelmo_init(yelmo1,filename=ctl%path_par,grid_def="none",time=ctl%time_init, &
                    load_topo=.FALSE.,domain=ctl%domain,grid_name=ctl%grid_name)

    ! === Define initial topography ===
    call calvmip_init(yelmo1%bnd%z_bed,yelmo1%grd%x,yelmo1%grd%y,yelmo1%par%domain)

    if (.not. yelmo1%par%use_restart) then
        ! If no restart, set ice thickness to zero
        yelmo1%tpo%now%H_ice = 0.0
        yelmo1%tpo%now%z_srf = yelmo1%bnd%z_bed 
        select case(trim(ctl%exp))
            case("advection")
            call LSFadvection(yelmo1%tpo%now%lsf,yelmo1%bnd%z_bed,yelmo1%tpo%par%dx)
        case DEFAULT 
            call LSFinit(yelmo1%tpo%now%lsf,yelmo1%tpo%now%H_ice,yelmo1%tpo%now%z_srf,yelmo1%tpo%par%dx)
        end select
    end if

    

    ! === Define additional boundary conditions =====

    yelmo1%bnd%z_sl     = 0.0
    yelmo1%bnd%bmb_shlf = 0.0  
    yelmo1%bnd%T_shlf   = yelmo1%bnd%c%T0  
    yelmo1%bnd%H_sed    = 0.0 

    yelmo1%bnd%T_srf    = 223.15 
    yelmo1%bnd%Q_geo    = 42.0

    select case(trim(ctl%exp))
        case("advection")
            yelmo1%bnd%smb      = 0.0
        case DEFAULT 
            yelmo1%bnd%smb      = 0.3
    end select    

    ! Check boundary values 
    call yelmo_print_bound(yelmo1%bnd)

    ! Initialize state variables (dyn,therm,mat)
    call yelmo_init_state(yelmo1,time=ctl%time_init,thrm_method="robin")

    ! == Write initial state ==
 
    ! 2D file 
    call yelmo_write_init(yelmo1,ctl%file2D,time_init=ctl%time_init,units="years")
    call yelmo_write_step(yelmo1,ctl%file2D,time=ctl%time_init)  
    
    ! 1D file 
    call yelmo_write_reg_init(yelmo1,ctl%file1D,time_init=ctl%time_init,units="years",mask=yelmo1%bnd%ice_allowed)
    call yelmo_write_reg_step(yelmo1,ctl%file1D,time=ctl%time_init) 

    ! Store default parameters
    yelmo_ref = yelmo1

    ! Set calving mask if needed
    if (trim(yelmo1%tpo%par%calv_flt_method) .eq. "kill-pos") then
        call set_calving_mask(yelmo1%bnd%calv_mask,yelmo1%grd%x,yelmo1%grd%y,r_lim=750e3_wp)
    end if

    ! Advance timesteps
    do n = 1, ceiling((ctl%time_end-ctl%time_init)/ctl%dtt)

        ! Get current time 
        time = ctl%time_init + n*ctl%dtt
        
        ! if (time .lt. 10e3) then
        !     yelmo1%dyn%par%solver = "sia"
        ! else
        !     yelmo1%dyn%par%solver = yelmo_ref%dyn%par%solver
        ! end if

        ! == Yelmo ice sheet ===================================================
        call yelmo_update(yelmo1,time)
        
        ! == MODEL OUTPUT =======================================================
        if (mod(nint(time*100),nint(ctl%dt2D_out*100))==0) then 
            call yelmo_write_step(yelmo1,ctl%file2D,time=time)   
        end if 

        if (mod(nint(time*100),nint(ctl%dt1D_out*100))==0) then 
            call yelmo_write_reg_step(yelmo1,ctl%file1D,time=time) 
        end if 

    end do

    ! Write a restart file too
    call yelmo_restart_write(yelmo1,ctl%file_restart,time=time)

    ! Finalize program
    call yelmo_end(yelmo1,time=time)

    ! Stop timing 
    call yelmo_cpu_time(cpu_end_time,cpu_start_time,cpu_dtime)
    
    write(*,"(a,f12.3,a)") "Time  = ",cpu_dtime/60.0 ," min"
    write(*,"(a,f12.1,a)") "Speed = ",(1e-3*(ctl%time_end-ctl%time_init))/(cpu_dtime/3600.0), " kiloyears / hr"
    
contains



end program yelmo_calving
