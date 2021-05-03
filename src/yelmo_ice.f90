

module yelmo_ice

    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit

    use nml 
    use ncio 
    
    use yelmo_defs
    use yelmo_grid, only : yelmo_init_grid, calc_zeta
    use yelmo_timesteps, only : ytime_init, set_adaptive_timestep, set_adaptive_timestep_pc, set_pc_mask, calc_pc_eta,  &
                                calc_pc_tau_fe_sbe,calc_pc_tau_ab_sam, calc_pc_tau_heun, limit_adaptive_timestep, &
                                yelmo_timestep_write_init, yelmo_timestep_write, calc_adv3D_timestep1
    use yelmo_io 

    use yelmo_topography
    use yelmo_dynamics
    use yelmo_material
    use yelmo_thermodynamics
    use yelmo_boundaries
    use yelmo_data 
    use yelmo_regions 

    implicit none 

    private
    public :: yelmo_init, yelmo_init_topo, yelmo_init_state
    public :: yelmo_update, yelmo_update_equil, yelmo_end 
    public :: yelmo_print_bound, yelmo_set_time

contains

    subroutine yelmo_update(dom,time)
        ! Advance yelmo by calling yelmo_step n-times until new time is reached,
        ! using the predictor-corrector method (Cheng et al., 2017) 

        type(yelmo_class), intent(INOUT) :: dom
        real(prec), intent(IN) :: time

        ! Local variables 
        type(yelmo_class)  :: dom_ref 

        real(prec) :: dt_now, dt_max  
        real(prec) :: time_now 
        integer    :: n, nstep, n_now, n_dtmin 
        integer    :: n_lim 
        real(prec), parameter :: time_tol = 1e-5

        real(8)    :: cpu_time0, cpu_time1 
        real(prec) :: model_time0, model_time1 
        real(prec) :: speed  

        real(prec) :: H_mean, T_mean 
        real(prec), allocatable :: dt_save(:) 
        real(prec) :: dt_adv_min, dt_pi
        real(prec) :: eta_now, rho_now 
        integer    :: iter_redo, pc_k, iter_redo_tot 
        real(prec) :: ab_zeta 
        logical, allocatable :: pc_mask(:,:) 

        character(len=1012) :: kill_txt

        integer :: ij(2) 

        ! Determine which predictor-corrector (pc) method we are using for timestepping,
        ! assign scheme order and weights 
        select case(trim(dom%par%pc_method))

            case("FE-SBE")
                
                ! Order of the method 
                pc_k = 2 

                dom%tpo%par%dt_beta(1) = 1.0_prec 
                dom%tpo%par%dt_beta(2) = 0.0_prec 
                
                dom%tpo%par%dt_beta(3) = 1.0_prec 
                dom%tpo%par%dt_beta(4) = 0.0_prec 
            
            case("AB-SAM")
                
                ! Order of the method 
                pc_k = 2 

                ! Assume ratio zeta=dt_n/dt_nm1=1.0 to start
                dom%tpo%par%dt_zeta  = 1.0_prec
                
                dom%tpo%par%dt_beta(1) = 1.0_prec + dom%tpo%par%dt_zeta/2.0_prec 
                dom%tpo%par%dt_beta(2) = -dom%tpo%par%dt_zeta/2.0_prec 

                dom%tpo%par%dt_beta(3) = 0.5_prec 
                dom%tpo%par%dt_beta(4) = 0.5_prec 
            
            case("HEUN")
                
                ! Order of the method 
                pc_k = 2 

                dom%tpo%par%dt_beta(1) = 1.0_prec 
                dom%tpo%par%dt_beta(2) = 0.0_prec 
                
                dom%tpo%par%dt_beta(3) = 0.5_prec 
                dom%tpo%par%dt_beta(4) = 0.5_prec 
            
            case("RALSTON")
                
                write(error_unit,*) "This method does not work yet - the truncation error is incorrect."
                stop 

                ! Order of the method 
                pc_k = 2 

                dom%tpo%par%dt_beta(1) = 2.0_prec / 3.0_prec
                dom%tpo%par%dt_beta(2) = 0.0_prec 
                
                dom%tpo%par%dt_beta(3) = 0.25_prec 
                dom%tpo%par%dt_beta(4) = 0.75_prec 
            
            case DEFAULT 

                write(error_unit,*) "yelmo_udpate:: Error: pc_method does not match available options [FE-SBE, AB-SAM, HEUN]."
                write(error_unit,*) "pc_method = ", trim(dom%par%pc_method)
                stop 

        end select 
        
        ! Calculate timestepping factors for thermodynamics horizontal advection term too 
        select case(trim(dom%thrm%par%dt_method))

            case("FE")
                ! Forward Euler 

                dom%thrm%par%dt_beta(1) = 1.0_prec  
                dom%thrm%par%dt_beta(2) = 0.0_prec  

            case("AB") 
                ! Adams-Bashforth (with default weights assuming dt_n/dt_nm1=1.0)

                ! Assume ratio zeta=dt_n/dt_nm1=1.0 to start
                dom%thrm%par%dt_zeta = dom%tpo%par%dt_zeta

                dom%thrm%par%dt_beta(1) = 1.0_prec + dom%thrm%par%dt_zeta/2.0_prec 
                dom%thrm%par%dt_beta(2) = -dom%thrm%par%dt_zeta/2.0_prec 

            case("SAM") 
                ! Semi-implicit Adams–Moulton

                dom%thrm%par%dt_beta(1) = 0.5
                dom%thrm%par%dt_beta(2) = 0.5 

            case DEFAULT 

                write(error_unit,*) "yelmo_update:: Error: thermodynamics dt_method does not match available options [FE, SBE, AB, SAM]."
                write(error_unit,*) "thrm:: dt_method = ", trim(dom%thrm%par%dt_method)
                stop 

        end select 

        ! Load last model time (from dom%tpo, should be equal to dom%thrm)
        time_now = dom%tpo%par%time

        ! Determine maximum number of time steps to be iterated through   
        nstep   = ceiling( (time-time_now) / dom%par%dt_min )
        n_now   = 0  ! Number of timesteps saved 

        iter_redo_tot = 0   ! Number of times total this loop 
        allocate(dt_save(nstep))
        dt_save = missing_value 

        allocate(pc_mask(dom%grd%nx,dom%grd%ny))

        ! Iteration of yelmo component updates until external timestep is reached
        do n = 1, nstep

            ! Initialize cpu timing for this iteration
            call yelmo_cpu_time(cpu_time0) 
            model_time0 = time_now 

            ! Store initial state of yelmo object in case a reset is necessary due to instability
            dom_ref = dom 

            ! Update dt_max as a function of the total timestep 
            dt_max = max(time-time_now,0.0_prec)

            ! === Diagnose different adaptive timestep limits ===

            ! Calculate adaptive time step from CFL constraints 
            call set_adaptive_timestep(dt_adv_min,dom%time%dt_adv,dom%time%dt_diff,dom%time%dt_adv3D, &
                                dom%dyn%now%ux,dom%dyn%now%uy,dom%dyn%now%uz,dom%dyn%now%ux_bar,dom%dyn%now%uy_bar, &
                                dom%dyn%now%dd_ab_bar,dom%tpo%now%H_ice,dom%tpo%now%dHicedt,dom%par%zeta_ac, &
                                dom%tpo%par%dx,dom%par%dt_min,dt_max,dom%par%cfl_max,dom%par%cfl_diff_max) 
            
            ! Calculate adaptive timestep using proportional-integral (PI) methods
            call set_adaptive_timestep_pc(dt_pi,dom%time%pc_dt,dom%time%pc_eta,dom%par%pc_eps,dom%par%dt_min,dt_max, &
                                    dom%dyn%now%ux_bar,dom%dyn%now%uy_bar,dom%tpo%par%dx,pc_k,dom%par%pc_controller)

            ! Determine current time step to be used based on method of choice 
            select case(dom%par%dt_method) 

                case(0) 
                    ! No internal timestep, so nstep=1 
                    
                    dt_now = dt_max 

                case(1) 
                    ! Use CFL-based adaptive timestep

                    dt_now = dt_adv_min 

                case(2) 
                    ! Use PI adaptive timestep

                    dt_now = dt_pi
                    
                case DEFAULT 

                    write(error_unit,*) "yelmo_update:: Error: dt_method not recognized."
                    write(error_unit,*) "dt_method = ", dom%par%dt_method 
                    stop 

            end select 


            do iter_redo=1, dom%par%pc_n_redo 
                ! Prepare to potentially perform several iterations of the same timestep.
                ! If at the end of one iteration, the truncation error is too high, then 
                ! redo the timestep with a lower dt. 
                ! Repeat n_redo times or until error is reduced. 
                ! Note: with eg, n_redo=5, iter_redo=n_redo is rarely met, 
                ! so it is a good choice (not too high either allowing too many iterations).

                ! Advance the local time variable
                time_now   = time_now + dt_now
                if (abs(time-time_now) .lt. time_tol) time_now = time 
                
                ! Calculate dt_zeta (ratio of current to previous timestep)
                dom%tpo%par%dt_zeta  = dt_now / dom%time%pc_dt(1) 
                dom%thrm%par%dt_zeta = dom%tpo%par%dt_zeta

                if (trim(dom%par%pc_method) .eq. "AB-SAM") then 
                    ! Update the predictor weights, since they depend on the timestep 

                    dom%tpo%par%dt_beta(1) = 1.0_prec + dom%tpo%par%dt_zeta/2.0_prec 
                    dom%tpo%par%dt_beta(2) = -dom%tpo%par%dt_zeta/2.0_prec 

                end if 

                if (trim(dom%thrm%par%dt_method) .eq. "AB") then 
                    ! Update the predictor weights, since they depend on the timestep 

                    dom%thrm%par%dt_beta(1) = 1.0_prec + dom%thrm%par%dt_zeta/2.0_prec 
                    dom%thrm%par%dt_beta(2) = -dom%thrm%par%dt_zeta/2.0_prec 

                end if 

                if (dom%par%pc_filter_vel) then 

                    ! Modify ux/y_bar to use the average between the current and previous velocity solutions
                    dom%dyn%now%ux_bar = 0.5_prec*dom%dyn%now%ux_bar + 0.5_prec*dom%dyn%now%ux_bar_prev
                    dom%dyn%now%uy_bar = 0.5_prec*dom%dyn%now%uy_bar + 0.5_prec*dom%dyn%now%uy_bar_prev
                    
                end if 
                
                dom%tpo%par%pc_step = "predictor" 

                ! Step 1: Perform predictor step for topography
                ! (Update elevation, ice thickness, calving, etc.)
                call calc_ytopo(dom%tpo,dom%dyn,dom%thrm,dom%bnd,time_now,topo_fixed=dom%tpo%par%topo_fixed)
                
                ! Store predicted ice thickness for later use 
                ! Do it here to ensure all changes to H_ice are accounted for (mb, calving, etc)
                dom%tpo%now%H_ice_pred = dom%tpo%now%H_ice 

                ! Step 2: Update other variables using predicted ice thickness 
                
                ! Calculate dynamics (velocities and stresses) 
                call calc_ydyn(dom%dyn,dom%tpo,dom%mat,dom%thrm,dom%bnd,time_now)

                ! Calculate material (ice properties, viscosity, etc.)
                call calc_ymat(dom%mat,dom%tpo,dom%dyn,dom%thrm,dom%bnd,time_now)

                ! Calculate thermodynamics (temperatures and enthalpy)
                call calc_ytherm(dom%thrm,dom%tpo,dom%dyn,dom%mat,dom%bnd,time_now)            

                dom%tpo%par%pc_step = "corrector" 
                dom%tpo%par%time    = dom_ref%tpo%par%time

                if (dom%par%pc_filter_vel) then 
                    
                    ! Modify ux/y_bar to use the average between the current and previous velocity solutions
                    dom%dyn%now%ux_bar = 0.5_prec*dom%dyn%now%ux_bar + 0.5_prec*dom%dyn%now%ux_bar_prev
                    dom%dyn%now%uy_bar = 0.5_prec*dom%dyn%now%uy_bar + 0.5_prec*dom%dyn%now%uy_bar_prev
                    
                end if 
                
                ! Step 3: Finally, calculate topography corrector step
                ! (elevation, ice thickness, calving, etc.)
                call calc_ytopo(dom%tpo,dom%dyn,dom%thrm,dom%bnd,time_now,topo_fixed=dom%tpo%par%topo_fixed)    
                
                ! Store corrected ice thickness for later use 
                ! Do it here to ensure all changes to H_ice are accounted for (mb, calving, etc)
                dom%tpo%now%H_ice_corr = dom%tpo%now%H_ice 

                if (dom%par%pc_use_H_pred) then 
                    ! Experimental: continue using H_ice_pred as main H_ice variable 
                    ! (ie, only use H_ice_corr to help calculate truncation error) 
                    ! This gives great results for EISMINT, grl, etc.

                    dom%tpo%now%H_ice = dom%tpo%now%H_ice_pred 

                    ! Also recalculate z_srf and masks for consistency
                    call calc_ytopo(dom%tpo,dom%dyn,dom%thrm,dom%bnd,time_now,topo_fixed=.TRUE.)    
                    
                end if 

                ! Determine truncation error for ice thickness 
                select case(trim(dom%par%pc_method))
                    ! No default case necessary, handled earlier 

                    case("FE-SBE")
                        
                        ! FE-SBE truncation error 
                        call calc_pc_tau_fe_sbe(dom%time%pc_tau,dom%tpo%now%H_ice_corr,dom%tpo%now%H_ice_pred,dt_now)

                    case("AB-SAM")
                        
                        ! AB-SAM truncation error 
                        call calc_pc_tau_ab_sam(dom%time%pc_tau,dom%tpo%now%H_ice_corr,dom%tpo%now%H_ice_pred,dt_now, &
                                                                                                dom%tpo%par%dt_zeta)

                    case("HEUN")

                        ! HEUN truncation error (same as FE-SBE)
                        call calc_pc_tau_heun(dom%time%pc_tau,dom%tpo%now%H_ice_corr,dom%tpo%now%H_ice_pred,dt_now)

                    case("RALSTON")

                        call calc_pc_tau_fe_sbe(dom%time%pc_tau,dom%tpo%now%H_ice_corr,dom%tpo%now%H_ice_pred,dt_now)

                end select 

                ! Calculate eta for this timestep 
                call set_pc_mask(pc_mask,dom%tpo%now%H_ice_pred,dom%tpo%now%H_ice_corr,dom%tpo%now%f_grnd)
                eta_now = calc_pc_eta(dom%time%pc_tau,mask=pc_mask)

                ! Save masked pc_tau for output too 
                dom%time%pc_tau_masked = dom%time%pc_tau 
                where( .not. pc_mask) dom%time%pc_tau_masked = 0.0_prec 

                ij = maxloc(abs(dom%time%pc_tau_masked))

                !write(*,"(a,f12.5,f12.5,f12.5,2i4,2f10.2)") &
                !    "test: ", time_now, dt_now, eta_now, ij(1), ij(2), &
                !    dom%tpo%now%H_ice_pred(ij(1),ij(2)), &
                !    dom%tpo%now%H_ice_corr(ij(1),ij(2))
                
                ! Check if this timestep should be rejected:
                ! If the redo iteration is not the last allowed and the timestep is still larger  
                ! than the minimum, then if eta > tolerance or checkerboard found in tau,
                ! then redo iteration: reject this timestep and try again with a smaller timestep
                if ( (iter_redo .lt. dom%par%pc_n_redo .and. dt_now .gt. dom%par%dt_min) &
                     .and. eta_now .gt. dom%par%pc_tol ) then

                    ! Calculate timestep reduction to apply
                    !rho_now = 0.7_prec
                    !rho_now = (dom%par%pc_tol / eta_now)
                    rho_now = 0.7_prec*(1.0_prec+(eta_now-dom%par%pc_tol)/10.0_prec)**(-1.0_prec) 

                    ! Reset yelmo and time variables to beginning of timestep
                    dom      = dom_ref 
                    time_now = dom_ref%tpo%par%time
                    dt_now   = max(dt_now*rho_now,dom%par%dt_min)
                    
                else
                    ! Timestep converged properly or total redo iterations completed,
                    ! or no further timestep reduction is possible:
                    ! Exit the iteration loop for this timestep (ie, move on)

                    exit 
                
                end if 

            end do   ! End iteration loop 
            
            ! Calculate model speed for this iteration
            call yelmo_cpu_time(cpu_time1)
            call yelmo_calc_speed(speed,model_time0,time_now,cpu_time0,cpu_time1)

            ! Collect how many times the redo-iteration loop had to run 
            ! (not counting the first pass, which is not a redo)
            iter_redo_tot = iter_redo_tot + (iter_redo-1) 
            
            ! Update dt and eta vectors for last N timesteps (first index becomes latest value)
            dom%time%pc_dt = cshift(dom%time%pc_dt,shift=-1)
            dom%time%pc_dt(1) = dt_now 

            dom%time%pc_eta = cshift(dom%time%pc_eta,shift=-1)
            dom%time%pc_eta(1) = eta_now
        
            ! Save the current timestep and other data for log and for running mean 
            n_now = n_now + 1 
            dt_save(n_now) = dt_now 
            call yelmo_calc_running_stats(dom%time%dt_avg,dom%time%dts,dt_now,stat="mean")
            
            call yelmo_calc_running_stats(dom%time%model_speed,dom%time%model_speeds,speed,stat="mean")
            call yelmo_calc_running_stats(dom%time%eta_avg,dom%time%etas,dom%time%pc_eta(1),stat="mean")
            call yelmo_calc_running_stats(dom%time%ssa_iter_avg,dom%time%ssa_iters,real(dom%dyn%par%ssa_iter_now,prec),stat="mean")
            
            ! Extra diagnostic field, not necessary for normal runs
            call yelmo_calc_running_stats_2D(dom%time%pc_tau_max,dom%time%pc_taus,dom%time%pc_tau_masked,stat="max")

            if (dom%par%log_timestep) then 
                ! Write timestep file if desired

                call yelmo_timestep_write(dom%time%log_timestep_file,time_now,dt_now,dt_adv_min,dt_pi, &
                            dom%time%pc_eta(1),dom%time%pc_tau_masked,speed,dom%tpo%par%speed,dom%dyn%par%speed, &
                            dom%dyn%par%ssa_iter_now,iter_redo_tot)
            
            end if 

            ! Make sure model is still running well
            call yelmo_check_kill(dom,time_now)

            ! Additionally check if minimum timestep is reached continuously

            ! Set limit for check to be the last 50 timesteps or, if it is smaller,
            ! the total number of timesteps in this call of yelmo_update
            n_lim = min(50,nstep)

            if (n .ge. n_lim) then

                ! Check how many of the last 50 timesteps are dt = dt_min 
                n_dtmin = count(abs(dt_save((n-n_lim+1):n)-dom%par%dt_min) .lt. dom%par%dt_min*1e-3)
                
                ! If all n_lim timesteps are at minimum, kill program
                if (n_dtmin .ge. n_lim) then 
                    
                    write(kill_txt,"(a,i10,a,i10)") &
                        "Too many iterations of dt_min called continuously for this timestep.", &
                                            n_dtmin, " of ", n 

                    call yelmo_check_kill(dom,time_now,kill_request=kill_txt)

                end if 

            end if 

            ! Check if it is time to exit adaptive iterations
            ! (if current outer time step has been reached)
            if (abs(time_now - time) .lt. time_tol) exit 

        end do 

        ! Update regional calculations (for now entire domain with ice)
        call calc_yregions(dom%reg,dom%tpo,dom%dyn,dom%thrm,dom%mat,dom%bnd,mask=dom%bnd%ice_allowed)

        ! Compare with data 
        call ydata_compare(dom%dta,dom%tpo,dom%dyn,dom%mat,dom%thrm,dom%bnd,dom%par%domain)

        ! Write some diagnostics to make sure something useful is happening 
        if (yelmo_log) then

            n = count(dom%tpo%now%H_ice.gt.0.0)
            if (n .gt. 0.0) then 
                H_mean = sum(dom%tpo%now%H_ice,mask=dom%tpo%now%H_ice.gt.0.0)/real(n)
                T_mean = sum(dom%thrm%now%T_ice(:,:,dom%thrm%par%nz_aa),mask=dom%tpo%now%H_ice.gt.0.0)/real(n)
            else 
                H_mean = 0.0_prec 
                T_mean = 0.0_prec
            end if 

            n       = count(dt_save .ne. missing_value)
            n_dtmin = count( abs(dt_save(1:n)-dom%par%dt_min) .lt. dom%par%dt_min*1e-3 )

            write(*,"(a,f13.2,f10.2,f10.1,f8.1,2G10.3,1i6)") &
                        !"yelmo:: [time,speed,H,T,max(dt),min(dt),n(dt==dt_min)]:", &
                        "yelmo:: timelog:", &
                            time_now, dom%time%model_speed, H_mean, T_mean,  &
                                            maxval(dt_save(1:n)), minval(dt_save(1:n)), n_dtmin
            
        end if 

        return

    end subroutine yelmo_update

    subroutine yelmo_update_equil(dom,time,time_tot,dt,topo_fixed,dyn_solver)
        ! Iterate yelmo solutions to equilibrate without updating boundary conditions

        type(yelmo_class), intent(INOUT) :: dom
        real(prec),        intent(IN)    :: time              ! [yr] Current time
        real(prec),        intent(IN)    :: time_tot          ! [yr] Equilibration time 
        real(prec),        intent(IN)    :: dt                ! Local dt to be used for all modules
        logical,           intent(IN)    :: topo_fixed        ! Should topography be fixed? 
        character(len=*),  intent(IN), optional :: dyn_solver
        
        ! Local variables 
        type(yelmo_class) :: dom_ref 
        real(prec) :: time_now  
        integer    :: n, nstep 
        
        ! Only run equilibration if time_tot > 0 

        if (time_tot .gt. 0.0) then 

            ! Save original model configuration 
            dom_ref = dom 

            ! Set new, temporary parameter values from arguments 
            dom%tpo%par%topo_fixed = topo_fixed 

            if (present(dyn_solver)) dom%dyn%par%solver = dyn_solver 

            ! Ensure during equilibration that at least 5 ssa iterations
            ! are allowed, for solvers that depend on ssa. Not strictly
            ! necessary, but potentially helps to get things going safely. 
            dom%dyn%par%ssa_iter_max = max(dom_ref%dyn%par%ssa_iter_max,5)
            
            ! Do not log timesteps for equilibration period,
            ! since time will be inconsistent.
            dom%par%log_timestep     = .FALSE. 

            ! Allow at least n=10 timestep redo iterations. Not strictly
            ! necessary, but potentially helps to get things going safely. 
            dom%par%pc_n_redo  = max(10,dom_ref%par%pc_n_redo)

            ! Set model time to input time 
            call yelmo_set_time(dom,time)

            write(*,*) 
            write(*,*) "Starting equilibration steps, time to run [yrs]: ", time_tot 

            do n = 1, ceiling(time_tot/dt)

                time_now = time + n*dt
                call yelmo_update(dom,time_now)

            end do

            ! Restore original model choices
            dom%par      = dom_ref%par 
            dom%tpo%par  = dom_ref%tpo%par
            dom%dyn%par  = dom_ref%dyn%par 
            dom%thrm%par = dom_ref%thrm%par  
            
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

                write(error_unit,*) "yelmo_init:: Error: parameter grid_def not recognized: "//trim(grid_def)
                stop 

        end select 

        ! Check that grid has been defined properly 
        if (.not. allocated(dom%grd%x)) then 
            write(error_unit,*) "yelmo_init:: Error: ygrid has not been properly defined yet."
            write(error_unit,*) "(Either use yelmo_init_grid externally with desired grid parameters &
                       &or set grid_def=['name','file'])"
            stop 
        end if 

        ! Calculate zeta_aa and zeta_ac 
        call calc_zeta(dom%par%zeta_aa,dom%par%zeta_ac,dom%par%nz_ac,dom%par%nz_aa, &
                                                    dom%par%zeta_scale,dom%par%zeta_exp)

        ! Initialize ytime information here too 
        call ytime_init(dom%time,dom%grd%nx,dom%grd%ny,dom%par%nz_aa,dom%par%dt_min,dom%par%pc_eps)

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

        call ymat_alloc(dom%mat%now,dom%mat%par%nx,dom%mat%par%ny,dom%mat%par%nz_aa,dom%mat%par%nz_ac,dom%mat%par%n_iso)
        
        write(*,*) "yelmo_init:: material initialized."
        
        ! == thermodynamics == 
        
        call ytherm_par_load(dom%thrm%par,filename,dom%par%zeta_aa,dom%par%zeta_ac,dom%grd%nx,dom%grd%ny,dom%grd%dx,init=.TRUE.)

        call ytherm_alloc(dom%thrm%now,dom%thrm%par%nx,dom%thrm%par%ny,dom%thrm%par%nz_aa,dom%thrm%par%nzr_aa)
        
        write(*,*) "yelmo_init:: thermodynamics initialized."
        
        ! === Ensure consistency with specific parameters ===

        ! For particular case related to eismint,
        ! ensure that bmb is not used in mass conservation or vertical velocity 
        dom%dyn%par%use_bmb = dom%tpo%par%use_bmb

        ! Modify grid boundary treatment according to the experiment parameter 
        select case(trim(dom%par%experiment))

            case("EISMINT")

                dom%tpo%par%boundaries = "EISMINT"
                dom%dyn%par%boundaries = "EISMINT"
                
            case("MISMIP3D","TROUGH-F17","MISMIP+") 

                dom%tpo%par%boundaries = "MISMIP3D"
                dom%dyn%par%boundaries = "MISMIP3D"

                ! Consistency check 
                if (trim(dom%tpo%par%solver) .ne. "impl-upwind") then 
                    write(error_unit,*) "yelmo_init:: Error: the mass conservation solver for MISMIP3D experiments &
                    &must be 'impl-upwind' for stability. The 'expl' solver has not yet been designed to &
                    &handle ice advected at the border point nx-1, and thus oscillations can be produced. &
                    &Please set 'solver=impl-upwind'."
                    stop 
                end if 
            
            case("ISMIPHOM","periodic","periodic-xy") 
                ! Periodic boundary conditions in x and y, eg: X_1 = X_n-1; X_n = X_2

                dom%tpo%par%boundaries = "periodic"
                dom%dyn%par%boundaries = "periodic"

            case("slab") 

                dom%tpo%par%boundaries = "periodic" 
                dom%dyn%par%boundaries = "periodic"
                
            case("infinite") 
                ! Set border points equal to interior neighbors 
                ! (ajr: not fully implemented yet)

                dom%tpo%par%boundaries = "infinite"
                dom%dyn%par%boundaries = "infinite"
            
            case("periodic-x") 
                ! Periodic boundary conditions in x-direction,
                ! infinite in y-direction
                dom%tpo%par%boundaries = "periodic-x"
                dom%dyn%par%boundaries = "periodic-x"
                
        end select 

        ! == boundary == 
        
        ! Allocate the yelmo data objects (arrays, etc)
        call ybound_alloc(dom%bnd,dom%grd%nx,dom%grd%ny)

        ! Load region/basin masks
        call ybound_load_masks(dom%bnd,filename,"yelmo_masks",dom%par%domain,dom%par%grid_name)
        
        ! Update the ice_allowed mask based on domain definition 
        call ybound_define_ice_allowed(dom%bnd,dom%par%domain)
        
        write(*,*) "yelmo_init:: boundary intialized (loaded masks, set ref. topography)."
        
        ! == data == 
        
        call ydata_par_load(dom%dta%par,filename,dom%par%domain,dom%par%grid_name,init=.TRUE.)
        call ydata_alloc(dom%dta%pd,dom%grd%nx,dom%grd%ny,dom%par%nz_aa,dom%dta%par%pd_age_n_iso)

        ! Load data objects   
        call ydata_load(dom%dta,dom%bnd%ice_allowed)

        ! Set H_ice_ref and z_bed_ref to present-day ice thickness by default 
        dom%bnd%H_ice_ref = dom%dta%pd%H_ice 
        dom%bnd%z_bed_ref = dom%dta%pd%z_bed

        write(*,*) "yelmo_init:: data intialized (loaded data if desired)."
        
        ! == topography ==

        ! Determine how to manage initial topography (H_ice,z_bed)
        call yelmo_init_topo(dom,filename,time,load_topo)

        write(*,*) "yelmo_init:: topo intialized (loaded data if desired)."
        
        write(*,*) 
        write(*,*) "yelmo_init:: Initialization complete for domain: "// &
                   trim(dom%par%domain) 

        if (dom%par%log_timestep) then 
            ! Timestep file 
            call yelmo_timestep_write_init(dom%time%log_timestep_file,time,dom%grd%xc,dom%grd%yc,dom%par%pc_eps)
            call yelmo_timestep_write(dom%time%log_timestep_file,time,0.0_prec,0.0_prec,dom%time%pc_dt(1), &
                            dom%time%pc_eta(1),dom%time%pc_tau_masked,0.0_prec,0.0_prec,0.0_prec,dom%dyn%par%ssa_iter_now,0)
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
        logical :: init_topo_load 
        character(len=1028) :: init_topo_path
        character(len=56)   :: init_topo_names(3)
        integer             :: init_topo_state
        real(prec)          :: z_bed_f_sd 

        ! Initialize variables
         
        if (dom%par%use_restart) then 
            ! Load variables from a restart file

            call yelmo_restart_read_topo(dom,trim(dom%par%restart),time)

        else

            ! Load parameters related to topography initiaization 
            call nml_read(filename,"yelmo_init_topo","init_topo_load",  init_topo_load)
            call nml_read(filename,"yelmo_init_topo","init_topo_path",  init_topo_path)
            call nml_read(filename,"yelmo_init_topo","init_topo_names", init_topo_names)
            call nml_read(filename,"yelmo_init_topo","init_topo_state", init_topo_state)
            call nml_read(filename,"yelmo_init_topo","z_bed_f_sd",      z_bed_f_sd)

            call yelmo_parse_path(init_topo_path,dom%par%domain,dom%par%grid_name)
            
            ! Override parameter choice if command-line argument present 
            if (present(load_topo)) init_topo_load = load_topo 

            if (init_topo_load) then 
                ! =========================================
                ! Load topography data from netcdf file 

                call nc_read(init_topo_path,init_topo_names(1), dom%tpo%now%H_ice, missing_value=mv)
                call nc_read(init_topo_path,init_topo_names(2), dom%bnd%z_bed, missing_value=mv) 

                if (trim(init_topo_names(3)) .ne. ""     .and. &
                    trim(init_topo_names(3)) .ne. "none" .and. &
                    trim(init_topo_names(3)) .ne. "None") then 

                    call nc_read(init_topo_path,init_topo_names(3),dom%bnd%z_bed_sd)

                    ! Apply scaling to adjust z_bed depending on standard deviation
                    dom%bnd%z_bed = dom%bnd%z_bed + z_bed_f_sd*dom%bnd%z_bed_sd 

                else
                    dom%bnd%z_bed_sd = 0.0_prec 
                end if 

                ! Clean up ice thickness field 
                where (.not. dom%bnd%ice_allowed)  dom%tpo%now%H_ice = 0.0_prec 
                where(dom%tpo%now%H_ice  .lt. 1.0) dom%tpo%now%H_ice = 0.0_prec 

                ! Additionally modify initial topographic state 
                select case(init_topo_state)

                    case(0) 

                        ! Pass, use topography as loaded 

                    case(1) 
                        ! Remove ice, but do not adjust bedrock 

                        dom%tpo%now%H_ice = 0.0_prec 

                    case(2)
                        ! Remove ice, set bedrock to isostatically rebounded state 

                        dom%bnd%z_bed     = dom%bnd%z_bed + (rho_ice/rho_a)*dom%tpo%now%H_ice
                        dom%tpo%now%H_ice = 0.0_prec

                    case DEFAULT 

                        write(error_unit,*) "yelmo_init_topo:: Error: init_topo_state choice not recognized."
                        write(error_unit,*) "init_topo_state = ", init_topo_state 
                        stop 

                end select

            end if 

            ! Calculate topographic information (masks, etc)
            call calc_ytopo(dom%tpo,dom%dyn,dom%thrm,dom%bnd,time,topo_fixed=.TRUE.)
            
            ! Update regional calculations (for entire domain)
            call calc_yregions(dom%reg,dom%tpo,dom%dyn,dom%thrm,dom%mat,dom%bnd,mask=dom%bnd%ice_allowed)

        end if 

        write(*,*) "yelmo_init_topo:: range(z_bed):     ", minval(dom%bnd%z_bed),     maxval(dom%bnd%z_bed)
        write(*,*) "yelmo_init_topo:: range(z_bed_sd):  ", minval(dom%bnd%z_bed_sd),  maxval(dom%bnd%z_bed_sd)
        write(*,*) "yelmo_init_topo:: range(H_ice):     ", minval(dom%tpo%now%H_ice), maxval(dom%tpo%now%H_ice)
        write(*,*) "yelmo_init_topo:: scaling fac z_bed_f_sd: ", z_bed_f_sd 

        return 

    end subroutine yelmo_init_topo

    subroutine yelmo_init_state(dom,time,thrm_method)
        ! This subroutine is the second step to intializing 
        ! the state variables. It initializes ice temperatures,
        ! material properties and dynamics. It is called after the topography
        ! has already been loaded and the boundary variables (eg: T_srf,smb,Q_geo)
        ! are already initialized externally.
        ! The state variables are either calculated directly or
        ! loaded from a restart file. 
            
        implicit none 

        type(yelmo_class), intent(INOUT) :: dom
        real(prec),        intent(IN)    :: time  
        character(len=*),  intent(IN)    :: thrm_method 

        ! Local variables 
        integer :: q 
        character(len=256) :: dom_thrm_method 
        character(len=256) :: dom_thrm_rock_method 
        
        ! Store original model choices locally 
        dom_thrm_method      = dom%thrm%par%method 
        dom_thrm_rock_method = dom%thrm%par%rock_method 

        ! Impose initialization choices 
        dom%thrm%par%method      = thrm_method 
        dom%thrm%par%rock_method = "equil" 

        ! Initialize variables

        if (dom%par%use_restart) then 
            ! Load variables from a restart file 

            call yelmo_restart_read(dom,trim(dom%par%restart),time)

        else 

            ! Consistency check 
            if (trim(thrm_method) .ne. "linear" .and. trim(thrm_method) .ne. "robin" &
                .and. trim(thrm_method) .ne. "robin-cold") then 
                write(error_unit,*) "yelmo_init_state:: Error: temperature initialization must be &
                           &'linear', 'robin' or 'robin-cold' in order to properly prescribe &
                           &initial temperatures."
                stop 
            end if

            ! Run topo to make sure all fields are synchronized (masks, etc)
            call calc_ytopo(dom%tpo,dom%dyn,dom%thrm,dom%bnd,time,topo_fixed=.TRUE.)
            
            ! Calculate initial thermodynamic information
            dom%thrm%par%time = dble(time) - dom%par%dt_min
            call calc_ytherm(dom%thrm,dom%tpo,dom%dyn,dom%mat,dom%bnd,time)

            ! Calculate material information (with no dynamics), and set initial ice dep_time values
            
            dom%mat%par%time     = dble(time) - dom%par%dt_min
            dom%mat%now%dep_time = dom%mat%par%time

            call calc_ymat(dom%mat,dom%tpo,dom%dyn,dom%thrm,dom%bnd,time)

            ! Calculate initial dynamic state
            ! (normally dynamics is called right after topo, but it needs thermodynamic information,
            ! thus here it is called after initializing ytherm and ymat variables)

            ! Impose [high] beta value in case it hasn't been initialized (eg, in the case of cb_method=-1/beta_method=-1)
            ! This will be overwritten when cf_ref/beta are calculated internally
            if (maxval(dom%dyn%now%beta) .eq. 0.0_prec) then 
                dom%dyn%now%cf_ref = 1.0
                dom%dyn%now%c_bed  = dom%dyn%now%cf_ref*1e5
                dom%dyn%now%beta   = dom%dyn%now%c_bed
            end if
             
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
        dom%thrm%par%method      = dom_thrm_method 
        dom%thrm%par%rock_method = dom_thrm_rock_method

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
        call nml_read(filename,"yelmo","log_timestep",  par%log_timestep)
        call nml_read(filename,"yelmo","disable_kill",  par%disable_kill)
        call nml_read(filename,"yelmo","zeta_scale",    par%zeta_scale)
        call nml_read(filename,"yelmo","zeta_exp",      par%zeta_exp)
        call nml_read(filename,"yelmo","nz_aa",         par%nz_aa)
        call nml_read(filename,"yelmo","dt_method",     par%dt_method)
        call nml_read(filename,"yelmo","dt_min",        par%dt_min)
        call nml_read(filename,"yelmo","cfl_max",       par%cfl_max)
        call nml_read(filename,"yelmo","cfl_diff_max",  par%cfl_diff_max)
        call nml_read(filename,"yelmo","pc_method",     par%pc_method)
        call nml_read(filename,"yelmo","pc_controller", par%pc_controller)
        call nml_read(filename,"yelmo","pc_filter_vel", par%pc_filter_vel)
        call nml_read(filename,"yelmo","pc_use_H_pred", par%pc_use_H_pred)
        call nml_read(filename,"yelmo","pc_n_redo",     par%pc_n_redo)
        call nml_read(filename,"yelmo","pc_tol",        par%pc_tol)
        call nml_read(filename,"yelmo","pc_eps",        par%pc_eps)

        ! Overwrite parameter values with argument definitions if available
        if (present(domain))     par%domain    = trim(domain)
        if (present(grid_name))  par%grid_name = trim(grid_name)
        
        ! Parse filenames with grid information
        call yelmo_parse_path(par%grid_path,par%domain,par%grid_name)
        call yelmo_parse_path(par%restart,par%domain,par%grid_name)

        ! dt_min must be greater than zero 
        if (par%dt_min .eq. 0.0) then 
            write(*,*) "yelmo_par_load:: dt_min must be greater than zero."
            write(*,*) "dt_min = ", par%dt_min 
            stop 
        end if 

        ! Set restart flag based on 'restart' parameter 
        if (trim(par%restart) .eq. "None" .or. &
            trim(par%restart) .eq. "none" .or. &
            trim(par%restart) .eq. "no") then 
            par%use_restart = .FALSE. 
        else 
            par%use_restart = .TRUE. 
        end if 

        if (par%pc_eps .gt. par%pc_tol) then 
            write(error_unit,*) "yelmo_par_load:: error: pc_eps must be greater than pc_tol."
            write(error_unit,*) "pc_eps, pc_tol: ", par%pc_eps, par%pc_tol 
            stop 
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
        
        if (yelmo_log) then 
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

    subroutine yelmo_check_kill(dom,time,kill_request)
            
        use ieee_arithmetic

        implicit none 

        type(yelmo_class), intent(IN) :: dom
        real(prec),        intent(IN) :: time
        character(len=*), optional, intent(IN) :: kill_request 

        ! Local variables 
        integer :: i, j, k 
        logical :: kill_it, kill_it_H, kill_it_vel, kill_it_nan, kill_it_eta   
        character(len=512) :: kill_msg 
        real(prec) :: pc_eta_avg 
        character(len=3) :: pc_iter_str(10) 

        real(prec), parameter :: H_lim = 1e4   ! [m] 
        real(prec), parameter :: u_lim = 1e4   ! [m/a]

        kill_it_H   = .FALSE. 
        kill_it_vel = .FALSE. 
        kill_it_nan = .FALSE. 
        kill_it_eta = .FALSE. 

        pc_iter_str = "" 
        pc_iter_str(1)  = "n"
        pc_iter_str(2)  = "n-1"
        pc_iter_str(3)  = "n-2"
        pc_iter_str(4)  = "n-3"
        pc_iter_str(5)  = "n-4"
        pc_iter_str(6)  = "n-5"
        pc_iter_str(7)  = "n-6"
        pc_iter_str(8)  = "n-7"
        pc_iter_str(9)  = "n-8"
        pc_iter_str(10) = "n-9"

        if ( maxval(abs(dom%tpo%now%H_ice)) .ge. H_lim .or. &
             maxval(abs(dom%tpo%now%H_ice-dom%tpo%now%H_ice)) .ne. 0.0 ) then 
            
            kill_it_H = .TRUE. 
            kill_msg  = "Ice thickness too high or invalid."

        end if 

        if ( maxval(abs(dom%dyn%now%uxy_bar)) .ge. u_lim .or. &
             maxval(abs(dom%dyn%now%uxy_bar-dom%dyn%now%uxy_bar)) .ne. 0.0 ) then 

            kill_it_vel = .TRUE. 
            kill_msg    = "Depth-averaged velocity too fast or invalid."

        end if 

        ! Additionally check for NANs using intrinsic ieee_arithmetic module 
        do j = 1, dom%grd%ny 
        do i = 1, dom%grd%nx 
            if (ieee_is_nan(dom%dyn%now%uxy_bar(i,j)) .or. ieee_is_nan(dom%tpo%now%H_ice(i,j))) then 
                kill_it_nan = .TRUE. 
                write(kill_msg,*) "** NANs detected ** ... i, j: ", i, j 
                exit 
            end if 
        end do 
        end do 

        pc_eta_avg = sum(dom%time%pc_eta) / real(size(dom%time%pc_eta,1),prec) 

        if (pc_eta_avg .gt. 10.0*dom%par%pc_tol) then 
            kill_it_eta = .TRUE. 
            write(kill_msg,"(a,g12.4,a,10g12.4)") "mean[pc_eta] > [10*pc_tol]: pc_eta_avg = ", pc_eta_avg, &
                                                                                " | pc_eta: ", dom%time%pc_eta
        end if 

        ! Determine if model should be killed 
        kill_it = kill_it_H .or. kill_it_vel .or. kill_it_nan .or. kill_it_eta 

        ! Definitely kill the model if it was requested externally
        if (present(kill_request)) then 
            kill_it = .TRUE. 
            kill_msg = trim(kill_request)
        end if 

        if (kill_it .and. (.not. dom%par%disable_kill)) then 
            ! Model is not running properly, kill it. 

            write(error_unit,*) 
            write(error_unit,*) 
            write(error_unit,"(a)") "yelmo_check_kill:: Error: model is not running properly:"
            write(error_unit,"(a)") trim(kill_msg) 
            write(error_unit,*) 
            write(error_unit,"(a11,f15.3)")  "timestep    = ", time
            write(error_unit,*) 
            write(error_unit,"(a,2g12.4)")   "pc_eps, tol = ", dom%par%pc_eps, dom%par%pc_tol 
            write(error_unit,"(a,g12.4)")    "pc_eta_avg  = ", pc_eta_avg
            
            write(error_unit,"(a4,1x,2a12)") "iter", "pc_dt", "pc_eta"
            do k = 1, size(dom%time%pc_eta,1)
                write(error_unit,"(a4,1x,2g12.4)") trim(pc_iter_str(k)), dom%time%pc_dt(k), dom%time%pc_eta(k) 
            end do 

            write(error_unit,*) 
            write(error_unit,"(a16,2g14.4)") "range(H_ice):   ", minval(dom%tpo%now%H_ice), maxval(dom%tpo%now%H_ice)
            write(error_unit,"(a16,2g14.4)") "range(uxy_bar): ", minval(dom%dyn%now%uxy_bar), maxval(dom%dyn%now%uxy_bar)
            write(error_unit,*) 

            call yelmo_restart_write(dom,"yelmo_killed.nc",time=time) 

            write(error_unit,*) "Restart file written: "//"yelmo_killed.nc"
            write(error_unit,*) 
            write(error_unit,"(a,f15.3,a)") "Time =", time, ": stopping model (killed)." 
            write(error_unit,*) 

            stop "yelmo_check_kill error, see log."

        end if 

        return 

    end subroutine yelmo_check_kill

end module yelmo_ice


