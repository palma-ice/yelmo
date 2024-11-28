program yelmo_calving

    use nml 
    use ncio  
    use yelmo 

    use calving_benchmarks
    
    implicit none

    type(yelmo_class)     :: yelmo1
    
    type control_type
        character(len=56)  :: domain    
        character(len=256) :: outfldr
        character(len=256) :: file2D, file1D
        character(len=256) :: file_restart
        character(len=512) :: path_par 
        character(len=56)  :: exp

        real(wp) :: time_init, time_end, time, dtt
        real(wp) :: dt2D_out, dt1D_out

        real(wp) :: dx
    end type

    type(control_type) :: ctl

    
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

    ! Define domain based on experiment
    select case(trim(ctl%exp))
        case("exp1","exp2")
            ctl%domain = "circular"
        case("exp3","exp4","exp5")
            ctl%domain = "thule"
        case DEFAULT
            write(*,*) "ctl.exp = ",trim(ctl%domain), " not recognized."
            stop
    end select



contains



end program yelmo_calving