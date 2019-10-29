

module yelmo_ice

    use nml 
    use ncio 
    
    use yelmo_defs
    use yelmo_grid, only : yelmo_init_grid
    use yelmo_timesteps, only : set_adaptive_timestep, set_adaptive_timestep_fe_sbe, limit_adaptive_timestep, &
            yelmo_timestep_write_init, yelmo_timestep_write
    use yelmo_io 

    use yelmo_topography
    use yelmo_dynamics
    use yelmo_material
    use yelmo_thermodynamics
    use yelmo_boundaries
    use yelmo_data 
    use yelmo_regions 

    use velocity_hybrid_pd12 

    implicit none 

    private
    public :: yelmo_init, yelmo_init_topo, yelmo_init_state
    public :: yelmo_update, yelmo_update_equil, yelmo_end 
    public :: yelmo_print_bound, yelmo_set_time

contains

    subroutine yelmo_update(dom,time)
        ! Advance yelmo by calling yelmo_step N times until new time is reached,
        ! using the predictor-corrector method (Cheng et al., 2017) 

        type(yelmo_class), intent(INOUT) :: dom
        real(prec), intent(IN) :: time

        type(ytopo_class) :: tpo1 
        real(prec), allocatable :: H_ice_prev(:,:) 
        real(prec), allocatable :: dHdt_prev(:,:) 

        ! Local variables 
        real(prec) :: dt_now 
        real(prec) :: time_now, time_start 
        integer    :: n, nstep   
        logical    :: iter_exit 
        real(4)    :: cpu_start_time 
        real(prec), parameter :: time_tol = 1e-5

        real(prec) :: H_mean, T_mean 
        real(prec) :: dt_save(100) 
        real(prec) :: dt_adv 

        ! Load last model time (from dom%tpo, should be equal to dom%thrm)
        time_now = dom%tpo%par%time

        ! Initialize cpu timing 
        call cpu_time(cpu_start_time) 
        time_start = time_now 
        
        ! Determine maximum number of time steps to be iterated through   
        nstep  = (time-time_now) / dom%par%dtmin 

        dt_save = missing_value 


        ! Store local copy of ytopo object to use for predictor step 
        tpo1 = dom%tpo 

        ! Iterate of topo dynamics updates
        do n = 1, nstep

            !write(*,*) "xx1", time, time_now, dom%par%pc_dt, dom%par%pc_eta 

            if (write_timestep_log) then 
                ! Timestep file 
                call yelmo_timestep_write(write_timestep_log_file,time_now,dom%par%pc_dt,dom%par%pc_eta)
            end if 

            ! Determine current time step 
            if (dom%tpo%par%topo_fixed) then 
                ! No topo calcs, so nstep=1
                dt_now = time-dom%tpo%par%time
            else
                ! Use last calculated adaptive timestep 
                dt_now = dom%par%pc_dt 

                ! Finally, ensure timestep is within prescribed limits
                call limit_adaptive_timestep(dt_now,time_now,time,dom%par%dtmin,dom%par%dtmax)

            end if 

            dt_save(n) = dt_now 

            ! Advance the local time variable
            time_now = time_now + dt_now
            if (time-time_now .lt. time_tol) time_now = time 
            
!             if (yelmo_write_log) then 
!                 write(*,"(a,1f14.4,3g14.4)") "timestepping: ", time_now, dt_now, minval(dom%par%dt_adv), minval(dom%par%dt_diff)
!             end if 
            
            ! Step 1: Perform predictor step with temporary topography object 
            ! Calculate topography (elevation, ice thickness, calving, etc.)
            call calc_ytopo(tpo1,dom%dyn,dom%thrm,dom%bnd,time_now,topo_fixed=dom%tpo%par%topo_fixed)


            ! Step 2: Update other variables using predicted ice thickness 
            
            ! Calculate material (ice properties, viscosity, etc.)
            call calc_ymat(dom%mat,tpo1,dom%dyn,dom%thrm,dom%bnd,time_now)

            ! Calculate thermodynamics (temperatures and enthalpy)
            call calc_ytherm(dom%thrm,tpo1,dom%dyn,dom%mat,dom%bnd,time_now)
            
            ! Calculate dynamics (velocities and stresses)
            call calc_ydyn(dom%dyn,tpo1,dom%mat,dom%thrm,dom%bnd,time_now)
            
            ! Step 3: Finally, calculate corrector step with actual topography object 
            ! Calculate topography (elevation, ice thickness, calving, etc.)
            call calc_ytopo(dom%tpo,dom%dyn,dom%thrm,dom%bnd,time_now,topo_fixed=dom%tpo%par%topo_fixed)

            ! Calculate new adaptive timestep from predicted and corrected ice thicknesses 
!             call set_adaptive_timestep_fe_sbe(dom%par%pc_dt,dom%par%pc_eta,tpo1%now%H_ice,dom%tpo%now%H_ice, &
!                                                     dom%par%pc_ebs,time_now,time,dom%par%dtmin,dom%par%dtmax)

            call set_adaptive_timestep_fe_sbe(dom%par%pc_dt,dom%par%pc_eta,tpo1%now%H_ice,dom%tpo%now%H_ice, &
                                                    dom%par%pc_ebs,time_now,time,dom%par%dtmin,dom%par%dtmax, &
                                                    dom%dyn%now%ux_bar,dom%dyn%now%uy_bar,dom%tpo%par%dx,dom%par%cfl_max)

            ! Make sure model is still running well
            call yelmo_check_kill(dom,time_now)

            ! Check if it is time to exit adaptive iterations
            ! (if current outer time step has been reached)
            iter_exit = .FALSE. 
            if (abs(time_now - time) .lt. time_tol) iter_exit = .TRUE. 

            ! If current time step has been reached, exit loop 
            if (iter_exit) exit

        end do 

        ! Update regional calculations (for now entire domain with ice)
        call calc_yregions(dom%reg,dom%tpo,dom%dyn,dom%thrm,dom%mat,dom%bnd,mask=dom%bnd%ice_allowed)

        ! Compare with data 
        call ydata_compare(dom%dta,dom%tpo,dom%dyn,dom%thrm,dom%bnd)

        ! Calculate model speed [model-kyr / hr]
        call yelmo_calc_speed(dom%par%model_speed,dom%par%model_speeds,time_start,time_now,real(cpu_start_time,prec))

        ! Write some diagnostics to make sure something useful is happening 
        if (yelmo_write_log) then

            n = count(dom%tpo%now%H_ice.gt.0.0)
            if (n .gt. 0.0) then 
                H_mean = sum(dom%tpo%now%H_ice,mask=dom%tpo%now%H_ice.gt.0.0)/real(n)
                T_mean = sum(dom%thrm%now%T_ice(:,:,dom%thrm%par%nz_aa),mask=dom%tpo%now%H_ice.gt.0.0)/real(n)
            else 
                H_mean = 0.0_prec 
                T_mean = 0.0_prec
            end if 

            n = count(dt_save .ne. missing_value)

            !write(*,*) "xx2", time, time_now, dom%par%pc_dt, dom%par%pc_eta 

            write(*,"(a,f12.2,f8.1,2f10.1,20f7.2)") "yelmo:: [time,speed,H,T,dt]:", time_now, dom%par%model_speed, &
                                H_mean, T_mean, dt_save(1:n)
            
        end if 

        return

    end subroutine yelmo_update

    subroutine yelmo_update_original(dom,time)
        ! Advance yelmo by calling yelmo_step N times until new time is reached 

        type(yelmo_class), intent(INOUT) :: dom
        real(prec), intent(IN) :: time

        ! Local variables 
        real(prec) :: dt_now 
        real(prec) :: time_now, time_start 
        integer    :: n, nstep   
        integer    :: ntt
        logical    :: iter_exit 
        real(4)    :: cpu_start_time 
        real(prec), parameter :: time_tol = 1e-5

        real(prec) :: H_mean, T_mean 
        real(prec) :: dt_save(100) 

        ! Load last model time (from dom%tpo, should be equal to dom%thrm)
        time_now = dom%tpo%par%time

        ! Initialize cpu timing 
        call cpu_time(cpu_start_time) 
        time_start = time_now 
        
        ! Determine maximum number of time steps to be iterated through   
        nstep  = (time-time_now) / dom%par%dtmin 

        ! Reset number of thermodynamics timestep skips
        ntt = 0 

        dt_save = missing_value 

        ! Iterate of topo dynamics updates
        do n = 1, nstep

            ! Determine current time step 
            if (dom%tpo%par%topo_fixed) then 
                ! No topo calcs, so nstep=1
                dt_now = time-dom%tpo%par%time
            else
                ! Use adaptive time step
                call set_adaptive_timestep(dt_now,dom%par%dt_adv,dom%par%dt_diff,dom%par%dt_adv3D,time_now,time, &
                                    dom%dyn%now%ux,dom%dyn%now%uy,dom%dyn%now%uz,dom%dyn%now%ux_bar,dom%dyn%now%uy_bar, &
                                    dom%dyn%now%dd_ab_bar,dom%tpo%now%H_ice,dom%tpo%now%dHicedt,dom%par%zeta_ac, &
                                    dom%tpo%par%dx,dom%par%dtmin,dom%par%dtmax,dom%par%cfl_max,dom%par%cfl_diff_max) 
                
            end if 

            dt_save(n) = dt_now 

            ! Advance the local time variable
            time_now = time_now + dt_now
            if (time-time_now .lt. time_tol) time_now = time 
            
!             if (yelmo_write_log) then 
!                 write(*,"(a,1f14.4,3g14.4)") "timestepping: ", time_now, dt_now, minval(dom%par%dt_adv), minval(dom%par%dt_diff)
!             end if 
            
            ! Calculate topography (elevation, ice thickness, calving, etc.)
            call calc_ytopo(dom%tpo,dom%dyn,dom%thrm,dom%bnd,time_now,topo_fixed=dom%tpo%par%topo_fixed)

            ! Calculate dynamics (velocities and stresses)
            call calc_ydyn(dom%dyn,dom%tpo,dom%mat,dom%thrm,dom%bnd,time_now)
            
            ! Check if it is time to exit adaptive iterations
            ! (if current outer time step has been reached)
            iter_exit = .FALSE. 
            if (abs(time_now - time) .lt. time_tol) iter_exit = .TRUE. 

            ! Advance thermodynamics timestep counter 
            ntt = ntt + 1

            if (ntt .ge. dom%par%ntt .or. iter_exit) then 
                ! Call thermodynamics every ntt iteration, and for the last adaptive timestep 

                ! Calculate material (ice properties, viscosity, etc.)
                call calc_ymat(dom%mat,dom%tpo,dom%dyn,dom%thrm,dom%bnd,time_now)

                ! Calculate thermodynamics (temperatures and enthalpy)
                call calc_ytherm(dom%thrm,dom%tpo,dom%dyn,dom%mat,dom%bnd,time_now)
                
                ! Reset thermodynamics timestep counter back to zero
                ntt = 0 

            end if 

            ! Make sure model is still running well
            call yelmo_check_kill(dom,time_now)

            ! If current time step has been reached, exit loop 
            if (iter_exit) exit

        end do 

        ! Update regional calculations (for now entire domain with ice)
        call calc_yregions(dom%reg,dom%tpo,dom%dyn,dom%thrm,dom%mat,dom%bnd,mask=dom%bnd%ice_allowed)

        ! Compare with data 
        call ydata_compare(dom%dta,dom%tpo,dom%dyn,dom%thrm,dom%bnd)

        ! Calculate model speed [model-kyr / hr]
        call yelmo_calc_speed(dom%par%model_speed,dom%par%model_speeds,time_start,time_now,real(cpu_start_time,prec))

        ! Write some diagnostics to make sure something useful is happening 
        if (yelmo_write_log) then

            n = count(dom%tpo%now%H_ice.gt.0.0)
            if (n .gt. 0.0) then 
                H_mean = sum(dom%tpo%now%H_ice,mask=dom%tpo%now%H_ice.gt.0.0)/real(n)
                T_mean = sum(dom%thrm%now%T_ice(:,:,dom%thrm%par%nz_aa),mask=dom%tpo%now%H_ice.gt.0.0)/real(n)
            else 
                H_mean = 0.0_prec 
                T_mean = 0.0_prec
            end if 

            n = count(dt_save .ne. missing_value)

            write(*,"(a,f12.2,f8.1,2f10.1,20f7.2)") "yelmo:: [time,speed,H,T,dt]:", time_now, dom%par%model_speed, &
                                H_mean, T_mean, dt_save(1:n)
            
        end if 

        return

    end subroutine yelmo_update_original

    subroutine yelmo_update_equil(dom,time,time_tot,dt,topo_fixed,ssa_vel_max)
        ! Iterate yelmo solutions to equilibrate without updating boundary conditions

        type(yelmo_class), intent(INOUT) :: dom
        real(prec), intent(IN) :: time            ! [yr] Current time
        real(prec), intent(IN) :: time_tot        ! [yr] Equilibration time 
        real(prec), intent(IN) :: dt              ! Local dt to be used for all modules
        logical,    intent(IN) :: topo_fixed      ! Should topography be fixed? 
        real(prec), intent(IN) :: ssa_vel_max     ! Local vel limit to be used, if == 0.0, no ssa used

        ! Local variables 
        real(prec) :: time_now  
        integer    :: n, nstep 
        logical    :: use_ssa         ! Should ssa be active?  
        logical    :: dom_topo_fixed
        logical    :: dom_use_ssa 
        real(prec) :: dom_dtmax 
        integer    :: dom_ntt 
        real(prec) :: dom_ssa_vel_max 

        ! Only run equilibration if time_tot > 0 

        if (time_tot .gt. 0.0) then 

            ! Consistency check
            use_ssa = .FALSE. 
            if (ssa_vel_max .gt. 0.0 .and. dom%dyn%par%mix_method .ne. -2) use_ssa = .TRUE. 

            ! Save original model choices 
            dom_topo_fixed  = dom%tpo%par%topo_fixed 
            dom_use_ssa     = dom%dyn%par%use_ssa 
            dom_dtmax       = dom%par%dtmax
            dom_ntt         = dom%par%ntt 
            dom_ssa_vel_max = dom%dyn%par%ssa_vel_max

            ! Set model choices equal to equilibration choices 
            dom%tpo%par%topo_fixed  = topo_fixed 
            dom%dyn%par%use_ssa     = use_ssa 
            dom%par%dtmax           = dt 
            dom%par%ntt             = 1 
            dom%dyn%par%ssa_vel_max = ssa_vel_max

            ! Set model time to input time 
            call yelmo_set_time(dom,time)

            write(*,*) 
            write(*,*) "Starting equilibration steps, time to run [yrs]: ", time_tot 

            do n = 1, ceiling(time_tot/dt)

                time_now = time + n*dt
                call yelmo_update(dom,time_now)

            end do

            ! Restore original model choices 
            dom%tpo%par%topo_fixed  = dom_topo_fixed 
            dom%dyn%par%use_ssa     = dom_use_ssa 
            dom%par%dtmax           = dom_dtmax 
            dom%par%ntt             = dom_ntt 
            dom%dyn%par%ssa_vel_max = dom_ssa_vel_max

            write(*,*) 
            write(*,*) "Equilibration complete."
            write(*,*) 

        end if 

        ! Reset model time back to input time 
        call yelmo_set_time(dom,time)

        return

    end subroutine yelmo_update_equil
    
    subroutine yelmo_init(dom,filename,grid_def,time,load_topo,domain,grid_name)
        ! Initialize a yelmo domain, including the grid itself, 
        ! and all sub-components (topo,dyn,mat,therm,bound,data)

        implicit none

        type(yelmo_class) :: dom 
        character(len=*),  intent(IN) :: filename 
        character(len=*),  intent(IN) :: grid_def 
        real(prec),        intent(IN) :: time 
        logical, optional, intent(IN) :: load_topo 
        character(len=*),  intent(IN), optional :: domain
        character(len=*),  intent(IN), optional :: grid_name 

        ! Local variables  
        integer :: i, nz1, nz2 
        integer :: k 

        ! == yelmo == 

        ! Load the default yelmo parameters, then the domain specific parameters
        call yelmo_par_load(dom%par,filename,domain,grid_name)
        
        ! Define the grid for the current domain 
        select case(grid_def)

            case("none")
                ! Do nothing - grid has already been defined externally 

            case("name")
                ! Use grid_name to load grid parameters 

                call yelmo_init_grid(dom%grd,dom%par%grid_name)

            case("file")
                ! Load grid information from netcdf file 

                call yelmo_init_grid(dom%grd,dom%par%grid_path,dom%par%grid_name)

            case DEFAULT 

                write(*,*) "yelmo_init:: Error: parameter grid_def not recognized: "//trim(grid_def)
                stop 

        end select 

        ! Check that grid has been defined properly 
        if (.not. allocated(dom%grd%x)) then 
            write(*,*) "yelmo_init:: Error: ygrid has not been properly defined yet."
            write(*,*) "(Either use yelmo_init_grid externally with desired grid parameters &
                       &or set grid_def=['name','file'])"
            stop 
        end if 

        ! Define size of zeta_ac axis
        dom%par%nz_ac = dom%par%nz_aa-1 
        
        ! Allocate z-axes
        if (allocated(dom%par%zeta_aa)) deallocate(dom%par%zeta_aa)
        if (allocated(dom%par%zeta_ac)) deallocate(dom%par%zeta_ac)
        allocate(dom%par%zeta_aa(dom%par%nz_aa)) 
        allocate(dom%par%zeta_ac(dom%par%nz_ac))

        ! Calculate zeta_aa and zeta_ac 
        call calc_zeta(dom%par%zeta_aa,dom%par%zeta_ac,dom%par%zeta_scale,dom%par%zeta_exp)

        ! Allocate timestep arrays 
        if (allocated(dom%par%dt_adv))   deallocate(dom%par%dt_adv)
        if (allocated(dom%par%dt_diff))  deallocate(dom%par%dt_diff)
        if (allocated(dom%par%dt_adv3D)) deallocate(dom%par%dt_adv3D)
        allocate(dom%par%dt_adv(dom%grd%nx,dom%grd%ny))
        allocate(dom%par%dt_diff(dom%grd%nx,dom%grd%ny))
        allocate(dom%par%dt_adv3D(dom%grd%nx,dom%grd%ny,dom%par%nz_aa))
        
        dom%par%dt_adv   = 0.0 
        dom%par%dt_diff  = 0.0 
        dom%par%dt_adv3D = 0.0 

        dom%par%pc_dt    = dom%par%dtmin  
        dom%par%pc_eta   = 1.0 

        write(*,*) "yelmo_init:: yelmo initialized."
        
        ! == topography ==

        call ytopo_par_load(dom%tpo%par,filename,dom%grd%nx,dom%grd%ny,dom%grd%dx,init=.TRUE.)

        call ytopo_alloc(dom%tpo%now,dom%tpo%par%nx,dom%tpo%par%ny)
        
        write(*,*) "yelmo_init:: topography initialized."
        
        ! == dynamics == 

        call ydyn_par_load(dom%dyn%par,filename,dom%par%zeta_aa,dom%par%zeta_ac,dom%grd%nx,dom%grd%ny,dom%grd%dx,init=.TRUE.)

        call ydyn_alloc(dom%dyn%now,dom%dyn%par%nx,dom%dyn%par%ny,dom%dyn%par%nz_aa,dom%dyn%par%nz_ac)
        
        write(*,*) "yelmo_init:: dynamics initialized."
        
        ! == material == 

        call ymat_par_load(dom%mat%par,filename,dom%par%zeta_aa,dom%par%zeta_ac,dom%grd%nx,dom%grd%ny,dom%grd%dx,init=.TRUE.)

        call ymat_alloc(dom%mat%now,dom%mat%par%nx,dom%mat%par%ny,dom%mat%par%nz_aa,dom%mat%par%nz_ac)
        
        write(*,*) "yelmo_init:: material initialized."
        
        ! == thermodynamics == 
        
        call ytherm_par_load(dom%thrm%par,filename,dom%par%zeta_aa,dom%par%zeta_ac,dom%grd%nx,dom%grd%ny,dom%grd%dx,init=.TRUE.)

        call ytherm_alloc(dom%thrm%now,dom%thrm%par%nx,dom%thrm%par%ny,dom%thrm%par%nz_aa,dom%thrm%par%nz_ac,dom%thrm%par%nzr)
        
        write(*,*) "yelmo_init:: thermodynamics initialized."
        
        ! === Ensure consistency with specific parameters ===

        ! For particular case related to eismint,
        ! ensure that bmb is not used in mass conservation or vertical velocity 
        dom%dyn%par%use_bmb = dom%tpo%par%use_bmb

        ! For the case of SIA-only inland, ensure ssa will not be calculated
        ! and the shelves will be deleted
        if (dom%dyn%par%mix_method .eq. -2) then 
            !dom%tpo%par%calv_method = "kill" 
            dom%dyn%par%ssa_vel_max = 0.0 
        end if 

        ! Modify grid boundary treatment according to the experiment parameter 
        select case(trim(dom%par%experiment))

            case("EISMINT")

                dom%tpo%par%boundaries = "EISMINT"
                dom%dyn%par%boundaries = "EISMINT"
                
            case("MISMIP3D","TROUGH-F17") 

                dom%tpo%par%boundaries = "MISMIP3D"
                dom%dyn%par%boundaries = "MISMIP3D"

                ! Consistency check 
                if (trim(dom%tpo%par%solver) .ne. "impl-upwind") then 
                    write(*,*) "yelmo_init:: Error: the mass conservation solver for MISMIP3D experiments &
                    &must be 'impl-upwind' for stability. The 'expl' solver has not yet been designed to &
                    &handle ice advected at the border point nx-1, and thus oscillations can be produced. &
                    &Please set 'solver=impl-upwind'."
                    stop 
                end if 
                
        end select 

        ! == boundary and data == 
        
        ! Allocate the yelmo data objects (arrays, etc)
        call ybound_alloc(dom%bnd,dom%grd%nx,dom%grd%ny)

        ! Load region/basin masks
        call ybound_load_masks(dom%bnd,filename,"yelmo_masks",dom%par%domain,dom%par%grid_name)
        
        ! Update the ice_allowed mask based on domain definition 
        call ybound_define_ice_allowed(dom%bnd,dom%par%domain)
        
        write(*,*) "yelmo_init:: boundary intialized (loaded masks)."
        
        call ydata_par_load(dom%dta%par,filename,dom%par%domain,dom%par%grid_name,init=.TRUE.)
        call ydata_alloc(dom%dta%pd,dom%grd%nx,dom%grd%ny,dom%par%nz_aa)

        ! Load data objects 
        call ydata_load(dom%dta,dom%bnd%ice_allowed)

        write(*,*) "yelmo_init:: data intialized (loaded data if desired)."
        
        ! == topography ==

        ! Determine how to manage initial topography (H_ice,z_bed)
        call yelmo_init_topo(dom,filename,time,load_topo)

        write(*,*) "yelmo_init:: topo intialized (loaded data if desired)."
        

        ! Set bnd%H_ice_ref to present-day ice thickness by default 
        dom%bnd%H_ice_ref = dom%dta%pd%H_ice 

        write(*,*) 
        write(*,*) "yelmo_init:: Initialization complete for domain: "// &
                   trim(dom%par%domain) 

        if (write_timestep_log) then 
            ! Timestep file 
            call yelmo_timestep_write_init(write_timestep_log_file,time,dom%par%pc_ebs)
            call yelmo_timestep_write(write_timestep_log_file,time,dom%par%pc_dt,dom%par%pc_eta)
        end if 

        return

    end subroutine yelmo_init 

    subroutine yelmo_init_topo(dom,filename,time,load_topo)
        ! This subroutine is the first step to intializing 
        ! the state variables. It initializes only the topography
        ! to facilitate calculation of boundary variables (eg, T_srf),
        ! which should be initialized externally afterwards.
        ! The state variables are either calculated directly or
        ! loaded from a restart file. 
            
        implicit none 

        type(yelmo_class), intent(INOUT) :: dom
        character(len=*),  intent(IN)    :: filename
        real(prec),        intent(IN)    :: time 
        logical, optional, intent(IN)    :: load_topo 

        ! Local variables 
        logical :: load_topo_from_par 

        ! Initialize variables
         
        if (dom%par%use_restart) then 
            ! Load variables from a restart file

            call yelmo_restart_read_topo(dom,trim(dom%par%restart),time)

        else
            ! Determine whether topography has already been defined externally or
            ! to load initial topography data from file based on parameter choices 
            load_topo_from_par = .TRUE. 
            if (present(load_topo)) then 
                load_topo_from_par = load_topo 
            end if 

            if (load_topo_from_par) then 

                ! Load bedrock elevation from file 
                call ybound_load_z_bed(dom%bnd,filename,"yelmo_init_topo",dom%par%domain,dom%par%grid_name)
                
                ! Load ice thickness from file
                call ytopo_load_H_ice(dom%tpo,filename,"yelmo_init_topo",dom%par%domain,dom%par%grid_name,dom%bnd%ice_allowed)

            end if 

            ! Calculate topographic information (masks, etc)
            call calc_ytopo(dom%tpo,dom%dyn,dom%thrm,dom%bnd,time,topo_fixed=.TRUE.)

            ! Update regional calculations (for now entire domain with ice)
            call calc_yregions(dom%reg,dom%tpo,dom%dyn,dom%thrm,dom%mat,dom%bnd,mask=dom%bnd%ice_allowed)

        end if 

        return 

    end subroutine yelmo_init_topo

    subroutine yelmo_init_state(dom,filename,time,thrm_method)
        ! This subroutine is the second step to intializing 
        ! the state variables. It initializes ice temperatures,
        ! material properties and dynamics. It is called after the topography
        ! has already been loaded and the boundary variables (eg: T_srf,smb,Q_geo)
        ! are already initialized externally.
        ! The state variables are either calculated directly or
        ! loaded from a restart file. 
            
        implicit none 

        type(yelmo_class), intent(INOUT) :: dom
        character(len=*),  intent(IN)    :: filename
        real(prec),        intent(IN)    :: time  
        character(len=*),  intent(IN)    :: thrm_method 

        ! Local variables 
        character(len=256) :: dom_thrm_method 

        ! Store original model choices locally 
        dom_thrm_method = dom%thrm%par%method 

        ! Impose initialization choices 
        dom%thrm%par%method = thrm_method 
         
        ! Initialize variables

        if (dom%par%use_restart) then 
            ! Load variables from a restart file 

            call yelmo_restart_read(dom,trim(dom%par%restart),time)

        else 

            ! Consistency check 
            if (trim(thrm_method) .ne. "linear" .and. trim(thrm_method) .ne. "robin" &
                .and. trim(thrm_method) .ne. "robin-cold") then 
                write(*,*) "yelmo_init_state:: Error: temperature initialization must be &
                           &'linear', 'robin' or 'robin-cold' in order to properly prescribe &
                           &initial temperatures."
                stop 
            end if

            ! Run topo to make sure all fields are synchronized (masks, etc)
            call calc_ytopo(dom%tpo,dom%dyn,dom%thrm,dom%bnd,time,topo_fixed=.TRUE.)
        
            ! Calculate initial thermodynamic information
            dom%thrm%par%time = time - dom%par%dtmax
            call calc_ytherm(dom%thrm,dom%tpo,dom%dyn,dom%mat,dom%bnd,time)

            ! Calculate material information (with no dynamics), and set initial ice dep_time values
            
            dom%mat%par%time     = time - dom%par%dtmax
            dom%mat%now%dep_time = dom%mat%par%time

            call calc_ymat(dom%mat,dom%tpo,dom%dyn,dom%thrm,dom%bnd,time)

            ! Calculate initial dynamic state
            ! (normally dynamics is called right after topo, but it needs thermodynamic information,
            ! thus here it is called after initializing ytherm and ymat variables)

            ! Impose [high] beta value in case it hasn't been initialized (in the case of cb_method=-1/beta_method=-1)
            ! This will be overwritten when C_bed/beta are calculated internally
            dom%dyn%now%C_bed = 1e5
            dom%dyn%now%beta  = 1e5 
            
            ! Call dynamics 
            call calc_ydyn(dom%dyn,dom%tpo,dom%mat,dom%thrm,dom%bnd,time)

            ! Calculate material information again with updated dynamics
        
            call calc_ymat(dom%mat,dom%tpo,dom%dyn,dom%thrm,dom%bnd,time)

        end if 


        ! Re-run topo again to make sure all fields are synchronized (masks, etc)
        call calc_ytopo(dom%tpo,dom%dyn,dom%thrm,dom%bnd,time,topo_fixed=.TRUE.)
        
        ! Update regional calculations (for now entire domain with ice)
        call calc_yregions(dom%reg,dom%tpo,dom%dyn,dom%thrm,dom%mat,dom%bnd,mask=dom%bnd%ice_allowed)
        
        ! Restore original model choices 
        ! Impose initialization choices 
        dom%thrm%par%method = dom_thrm_method 
        
        return 

    end subroutine yelmo_init_state 

    subroutine yelmo_par_load(par,filename,domain,grid_name)

        type(yelmo_param_class), intent(OUT) :: par
        character(len=*), intent(IN) :: filename
        character(len=*), intent(IN), optional :: domain
        character(len=*), intent(IN), optional :: grid_name 

        call nml_read(filename,"yelmo","domain",        par%domain)
        call nml_read(filename,"yelmo","grid_name",     par%grid_name)
        call nml_read(filename,"yelmo","grid_path",     par%grid_path)
        call nml_read(filename,"yelmo","experiment",    par%experiment)
        call nml_read(filename,"yelmo","restart",       par%restart)
        call nml_read(filename,"yelmo","zeta_scale",    par%zeta_scale)
        call nml_read(filename,"yelmo","zeta_exp",      par%zeta_exp)
        call nml_read(filename,"yelmo","nz_aa",         par%nz_aa)
        call nml_read(filename,"yelmo","dtmin",         par%dtmin)
        call nml_read(filename,"yelmo","dtmax",         par%dtmax)
        call nml_read(filename,"yelmo","ntt",           par%ntt)
        call nml_read(filename,"yelmo","cfl_max",       par%cfl_max)
        call nml_read(filename,"yelmo","cfl_diff_max",  par%cfl_diff_max)
        call nml_read(filename,"yelmo","pc_ebs",        par%pc_ebs)

        ! Overwrite parameter values with argument definitions if available
        if (present(domain))     par%domain    = trim(domain)
        if (present(grid_name))  par%grid_name = trim(grid_name)
        
        ! Parse filenames with grid information
        call yelmo_parse_path(par%grid_path,par%domain,par%grid_name)
        call yelmo_parse_path(par%restart,par%domain,par%grid_name)

        ! dtmin must be greater than zero 
        if (par%dtmin .eq. 0.0) then 
            write(*,*) "yelmo_par_load:: dtmin must be greater than zero."
            write(*,*) "dtmin = ", par%dtmin 
            stop 
        end if 

        if (par%dtmin .gt. par%dtmax) then 
            write(*,*) "yelmo_par_load:: dtmin must be less than or equal to dtmax."
            write(*,*) "dtmin = ", par%dtmin
            write(*,*) "dtmax = ", par%dtmax
        end if 

        ! Ensure ntt is at least 1 
        if (par%ntt .lt. 1) par%ntt = 1 

        ! Set restart flag based on 'restart' parameter 
        if (trim(par%restart) .eq. "None" .or. &
            trim(par%restart) .eq. "none" .or. &
            trim(par%restart) .eq. "no") then 
            par%use_restart = .FALSE. 
        else 
            par%use_restart = .TRUE. 
        end if 

        return

    end subroutine yelmo_par_load
    
    subroutine yelmo_end(dom,time)
        ! Deallocate yelmo objects 

        implicit none 

        type(yelmo_class), intent(INOUT) :: dom 
        real(prec), intent(IN) :: time 

        call ytopo_dealloc(dom%tpo%now)
        call ydyn_dealloc(dom%dyn%now)
        call ymat_dealloc(dom%mat%now)
        call ytherm_dealloc(dom%thrm%now)
        call ybound_dealloc(dom%bnd)
        call ydata_dealloc(dom%dta%pd)
        
        if (yelmo_write_log) then 
            write(*,*)
            write(*,*) "yelmo_end:: time = ", time  
            write(*,*) "yelmo_end:: Finished."
            write(*,*) 
        end if 
        
        return 

    end subroutine yelmo_end

    subroutine yelmo_print_bound(bnd)

        implicit none 

        type(ybound_class), intent(IN) :: bnd

        write(*,*) "ybound:: z_bed:  ",   minval(bnd%z_bed),    maxval(bnd%z_bed)
        write(*,*) "ybound:: z_sl:  ",    minval(bnd%z_sl),     maxval(bnd%z_sl)
        write(*,*) "ybound:: H_sed: ",    minval(bnd%H_sed),    maxval(bnd%H_sed)
        write(*,*) "ybound:: H_w: ",      minval(bnd%H_w),      maxval(bnd%H_w)
        write(*,*) "ybound:: smb: ",      minval(bnd%smb),      maxval(bnd%smb)
        write(*,*) "ybound:: T_srf: ",    minval(bnd%T_srf),    maxval(bnd%T_srf)
        write(*,*) "ybound:: T_shlf: ",   minval(bnd%T_shlf),   maxval(bnd%T_shlf)
        write(*,*) "ybound:: bmb_shlf: ", minval(bnd%bmb_shlf), maxval(bnd%bmb_shlf)
        write(*,*) "ybound:: Q_geo: ",    minval(bnd%Q_geo),    maxval(bnd%Q_geo)
        
        return 

    end subroutine yelmo_print_bound 

    subroutine yelmo_set_time(dom,time)

        implicit none 

        type(yelmo_class), intent(INOUT) :: dom
        real(prec), intent(IN) :: time
        
        dom%tpo%par%time      = time 
        dom%dyn%par%time      = time 
        dom%mat%par%time      = time 
        dom%thrm%par%time     = time 
        
        return 

    end subroutine yelmo_set_time 

    subroutine calc_zeta(zeta_aa,zeta_ac,zeta_scale,zeta_exp)
        ! Calculate the vertical layer-edge axis (vertical ac-nodes)
        ! and the vertical cell-center axis (vertical aa-nodes),
        ! including an extra zero-thickness aa-node at the base and surface

        implicit none 

        real(prec), intent(INOUT)  :: zeta_aa(:) 
        real(prec), intent(INOUT)  :: zeta_ac(:) 
        character(*), intent(IN)   :: zeta_scale 
        real(prec),   intent(IN)   :: zeta_exp 

        ! Local variables
        integer :: k, nz_aa, nz_ac 

        nz_aa  = size(zeta_aa)
        nz_ac  = size(zeta_ac)   ! == nz_aa - 1 

        ! Initially define a linear zeta scale 
        ! Base = 0.0, Surface = 1.0 
        do k = 1, nz_ac
            zeta_ac(k) = 0.0 + 1.0*(k-1)/real(nz_ac-1)
        end do 

        ! Scale zeta to produce different resolution through column if desired
        ! zeta_scale = ["linear","exp","wave"]
        select case(trim(zeta_scale))
            
            case("exp")
                ! Increase resolution at the base 
                zeta_ac = zeta_ac**(zeta_exp) 

            case("tanh")
                ! Increase resolution at base and surface 

                zeta_ac = tanh(1.0*pi*(zeta_ac-0.5))
                zeta_ac = zeta_ac - minval(zeta_ac)
                zeta_ac = zeta_ac / maxval(zeta_ac)

            case DEFAULT
            ! Do nothing, scale should be linear as defined above
        
        end select  
        
        ! Get zeta_aa (between zeta_ac values, as well as at the base and surface)
        zeta_aa(1) = 0.0 
        do k = 2, nz_aa-1
            zeta_aa(k) = 0.5 * (zeta_ac(k-1)+zeta_ac(k))
        end do 
        zeta_aa(nz_aa) = 1.0 

        return 

    end subroutine calc_zeta
    
    subroutine yelmo_check_kill(dom,time)
        
        implicit none 

        type(yelmo_class), intent(IN) :: dom
        real(prec),        intent(IN) :: time
        
        ! Local variables 
        logical :: kill_it  
        real(prec), parameter :: H_lim = 1e5   ! [m] 
        real(prec), parameter :: u_lim = 1e5   ! [m/a]

        kill_it = .FALSE. 

        if ( maxval(abs(dom%tpo%now%H_ice)) .ge. H_lim .or. &
             maxval(abs(dom%tpo%now%H_ice-dom%tpo%now%H_ice)) .ne. 0.0 ) kill_it = .TRUE. 

        if ( maxval(abs(dom%dyn%now%uxy_bar)) .ge. u_lim ) kill_it = .TRUE. 

        if (kill_it) then 
            ! Model has probably crashed, kill it. 

            write(*,*) 
            write(*,*) 
            write(*,"(a)") "yelmo_check_kill:: Error: model has likely crashed, very high value of H_ice or uxy_bar found."
            write(*,"(a11,f15.3)")   "timestep = ", time 
            write(*,"(a16,2g14.4)") "range(H_ice):   ", minval(dom%tpo%now%H_ice), maxval(dom%tpo%now%H_ice)
            write(*,"(a16,2g14.4)") "range(uxy_bar): ", minval(dom%dyn%now%uxy_bar), maxval(dom%dyn%now%uxy_bar)
            write(*,*) 
            write(*,*) 
            write(*,*) "Stopping model."
            write(*,*) 

            stop "yelmo_check_kill error, see log."

        end if 

        return 

    end subroutine yelmo_check_kill 

end module yelmo_ice


