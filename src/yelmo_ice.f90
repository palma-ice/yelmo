

module yelmo_ice
    
    use nml 
    use ncio 
    
    use yelmo_defs
    use yelmo_grid, only : yelmo_init_grid, calc_zeta
    use yelmo_timesteps, only : ytime_init, set_pc_beta_coefficients, set_adaptive_timestep, set_adaptive_timestep_pc,   &
                                set_pc_mask, calc_pc_eta, calc_pc_tau_fe_sbe,calc_pc_tau_ab_sam, calc_pc_tau_heun,  &
                                limit_adaptive_timestep, yelmo_timestep_write_init, yelmo_timestep_write, calc_adv3D_timestep1
    use yelmo_tools, only : smooth_gauss_2D, adjust_topography_gradients
    use yelmo_io 

    use yelmo_topography
    use lsf_module, only : LSFinit
    use yelmo_dynamics
    use yelmo_material
    use yelmo_thermodynamics
    use yelmo_boundaries
    use yelmo_data 
    use yelmo_regions 

    use topography, only : remove_englacial_lakes
    use mass_conservation, only : calc_G_boundaries, check_mass_conservation, apply_tendency
    use variable_io, only : load_var_io_table
    !$  use omp_lib

    implicit none 

    private
    public :: yelmo_init, yelmo_init_topo, yelmo_init_state
    public :: yelmo_update, yelmo_update_equil, yelmo_end
    public :: yelmo_print_bound, yelmo_set_time

contains
    
    subroutine yelmo_update(dom,time,file_diagnostics)
        ! Advance yelmo by calling yelmo_step n-times until new time is reached,
        ! using the predictor-corrector method (Cheng et al., 2017) 

        type(yelmo_class), intent(INOUT) :: dom
        real(wp),          intent(IN)    :: time
        character(len=*),  intent(IN), optional :: file_diagnostics

        ! Local variables 
        type(yelmo_class)  :: dom_ref 

        real(wp) :: dt_now, dt_max, dt_max_0
        real(wp) :: time_now 
        integer  :: n, nstep, n_now, n_dtmin 
        integer  :: n_lim 
        real(wp), parameter :: time_tol = 1e-5

        real(8)  :: cpu_time0, cpu_time1 
        real(wp) :: model_time0, model_time1 
        real(wp) :: speed  

        real(wp) :: H_mean, T_mean 
        real(wp), allocatable :: dt_save(:) 
        real(wp) :: dt_adv_min, dt_pi
        real(wp) :: eta_now, rho_now 
        integer  :: iter_redo, iter_redo_tot 
        real(wp) :: ab_zeta 
        logical, allocatable :: pc_mask(:,:) 

        character(len=1012) :: kill_txt

        integer  :: ij(2) 
        real(wp) :: max_dt_used 
        real(wp) :: min_dt_used 

        logical, parameter :: update_others_pc  = .FALSE. 
        logical, parameter :: very_verbose      = .FALSE. 
        logical, parameter :: check_mb          = .FALSE. 

        !$ logical, parameter :: l_write_timer=.false.
        !$ real(8) :: time1, time2

        ! Safety: check status of model object, 
        ! Has it been initialized?
        if (.not. allocated(dom%tpo%now%H_ice)) then 
            write(io_unit_err,*) 
            write(io_unit_err,*) "yelmo_update:: Error: Yelmo object does not appear to be initialized/allocated."
            write(io_unit_err,*) "is_allocated(dom%tpo%now%H_ice) = ", allocated(dom%tpo%now%H_ice)
            stop 
        end if 

        ! Assume ratio zeta=dt_n/dt_nm1=1.0 to start
        dom%tpo%par%dt_zeta  = 1.0_wp
        dom%thrm%par%dt_zeta = dom%tpo%par%dt_zeta

        ! Determine which predictor-corrector (pc) method we are using for timestepping,
        ! assign scheme order and weights 
        call set_pc_beta_coefficients(dom%tpo%par%dt_beta,dom%tpo%par%pc_k, &
                                            dom%tpo%par%dt_zeta,dom%par%pc_method)

        ! Calculate timestepping factors for thermodynamics horizontal advection term too 
        call set_pc_beta_coefficients(dom%thrm%par%dt_beta,dom%thrm%par%pc_k, &
                                            dom%thrm%par%dt_zeta,dom%thrm%par%dt_method)

        ! Load last model time (from dom%tpo, should be equal to dom%thrm)
        time_now = dom%tpo%par%time

        ! Determine maximum number of time steps to be iterated through 
        ! Note: ensure at least one step will be performed, so that even 
        ! if the time is already up-to-date, masks etc can be calculated
        ! if desired. 
        nstep   = ceiling( (time-time_now) / dom%par%dt_min )
        nstep   = max(nstep,1)  
        n_now   = 0  ! Number of timesteps saved 

        ! Get initial estimate of maximum timestep
        dt_max_0 = max(time-time_now,0.0_wp)

        iter_redo_tot = 0   ! Number of times total this loop 
        allocate(dt_save(nstep))
        dt_save = missing_value 

        dom%time%eta_avg      = missing_value
        dom%time%ssa_iter_avg = missing_value

        ! Initialize rate averages
        call calc_ytopo_rates(dom%tpo,dom%bnd,time,dt=0.0_wp,step="init",check_mb=check_mb)

        allocate(pc_mask(dom%grd%nx,dom%grd%ny))

        ! Calculate filtered bedrock elevations adjusted for sea level on top (ie, water depth)
        !$ time1 = omp_get_wtime()
        dom%tpo%now%z_bed_filt = dom%bnd%z_bed - dom%bnd%z_sl
        if (dom%tpo%par%zb_sigma .gt. 0.0) then 
          call smooth_gauss_2D(dom%tpo%now%z_bed_filt,dom%tpo%par%dx, &
                               dom%tpo%par%zb_sigma / dom%tpo%par%dx)
        end if
        !$ time2 = omp_get_wtime()
        !$ if(l_write_timer) print *,'TIME smooth_gauss_2D',time2-time1
        
        ! ajr restart check:
        !call yelmo_restart_write(dom,"./yelmo_restart_init_loop.nc",time)

        ! Iteration of yelmo component updates until external timestep is reached
        do n = 1, nstep

            ! Initialize cpu timing for this iteration
            call yelmo_cpu_time(cpu_time0) 
            model_time0 = time_now 

            ! Store initial state of yelmo object in case a reset is necessary due to instability
            dom_ref = dom 

            ! Update dt_max as a function of the total timestep 
            dt_max = max(time-time_now,0.0_wp)

            ! === Diagnose different adaptive timestep limits ===

            ! Calculate adaptive time step from CFL constraints 
            call set_adaptive_timestep(dt_adv_min,dom%time%dt_adv,dom%time%dt_diff,dom%time%dt_adv3D, &
                                dom%dyn%now%ux,dom%dyn%now%uy,dom%dyn%now%uz,dom%dyn%now%ux_bar,dom%dyn%now%uy_bar, &
                                dom%tpo%now%H_ice,dom%tpo%now%dHidt,dom%par%zeta_ac, &
                                dom%tpo%par%dx,dom%par%dt_min,dt_max,dom%par%cfl_max,dom%par%cfl_diff_max) 
            
            ! Calculate adaptive timestep using proportional-integral (PI) methods
            call set_adaptive_timestep_pc(dt_pi,dom%time%pc_dt,dom%time%pc_eta,dom%par%pc_eps,dom%par%dt_min,dt_max, &
                                    dom%dyn%now%ux_bar,dom%dyn%now%uy_bar,dom%tpo%par%dx,dom%tpo%par%pc_k,dom%par%pc_controller)

            ! ajr restart check:
            ! write(*,*) "Set timestep: ", n, time_now, dt_pi, dt_max
            ! write(*,*) "    ", dom%time%pc_dt
            ! write(*,*) "    ", dom%time%pc_eta
            ! write(*,*) "    ", dom%par%pc_eps

            ! Determine current time step to be used based on method of choice 
            select case(dom%par%dt_method) 

                case(0) 
                    ! No internal timestep, so nstep=1 
                    
                    dt_now = dt_max 

                case(1) 
                    ! Use CFL-based adaptive timestep

                    dt_now = dt_adv_min 

                case(2) 
                    ! Use minimum of PI adaptive timestep

                    dt_now = dt_pi

                case DEFAULT 

                    write(io_unit_err,*) "yelmo_update:: Error: dt_method not recognized."
                    write(io_unit_err,*) "dt_method = ", dom%par%dt_method 
                    stop 

            end select 

            if (.not. dom%time%pc_active) then 
                ! Override timestep choice and simply use a very small value for 
                ! first timestep. 

                dt_now = dom%par%dt_min

            end if 

            ! Finally, override all cases if time is already up-to-date.
            ! This is done here, to let all PC timestepping algorithms etc.
            ! update themselves while dt_now>0. 
            if (dt_max .eq. 0.0) then 

                dt_now = 0.0_wp 

            end if 

            do iter_redo=1, dom%par%pc_n_redo 
                ! Prepare to potentially perform several iterations of the same timestep.
                ! If at the end of one iteration, the truncation error is too high, then 
                ! redo the timestep with a lower dt. 
                ! Repeat n_redo times or until error is reduced. 
                ! Note: with eg, n_redo=5, iter_redo=n_redo is rarely met, 
                ! so it is a good choice (not too high either allowing too many iterations).

                ! Advance the local time variable
                time_now = time_now + dt_now
                if (abs(time-time_now) .lt. time_tol) time_now = time 
                
                ! Calculate dt_zeta (ratio of current to previous timestep)
                dom%tpo%par%dt_zeta  = dt_now / dom%time%pc_dt(1) 
                dom%thrm%par%dt_zeta = dom%tpo%par%dt_zeta


                if (trim(dom%par%pc_method) .eq. "AB-SAM") then
                    ! Update the predictor weights, since they depend on the timestep 

                    if (.not. dom%time%pc_active) then 
                        ! Only FE-SBE is available 

                        call set_pc_beta_coefficients(dom%tpo%par%dt_beta,dom%tpo%par%pc_k, &
                                                            dom%tpo%par%dt_zeta,"FE-SBE")

                    else 

                        call set_pc_beta_coefficients(dom%tpo%par%dt_beta,dom%tpo%par%pc_k, &
                                                            dom%tpo%par%dt_zeta,dom%par%pc_method)

                    end if 

                end if 

                if (trim(dom%thrm%par%dt_method) .eq. "AB") then 
                    ! Update the predictor weights, since they depend on the timestep 

                    if (.not. dom%time%pc_active) then 
                        ! Only FE is available 

                        call set_pc_beta_coefficients(dom%thrm%par%dt_beta,dom%thrm%par%pc_k, &
                                                                dom%thrm%par%dt_zeta,"FE")

                    else 
                        ! Update the weights like normal
                        
                        call set_pc_beta_coefficients(dom%thrm%par%dt_beta,dom%thrm%par%pc_k, &
                                                    dom%thrm%par%dt_zeta,dom%thrm%par%dt_method)
                        
                    end if 

                end if 

                ! Step 1: Perform predictor step for topography
                ! Get predicted new ice thickness and store it for later use
                !$ time1 = omp_get_wtime()
                ! call calc_ytopo_rk4(dom%tpo,dom%dyn,dom%mat,dom%thrm,dom%bnd,time,dom%tpo%par%topo_fixed)
                call calc_ytopo_pc(dom%tpo,dom%dyn,dom%mat,dom%thrm,dom%bnd,dom%dta,time_now,dom%tpo%par%topo_fixed,"predictor")
                !$ time2 = omp_get_wtime()
                !$ if(l_write_timer) print *,'TIME ytopo predictor',time2-time1

                ! Step 2: Calculate dynamics for predicted ice thickness 

                !$ time1 = omp_get_wtime()
                call calc_ydyn(dom%dyn,dom%tpo,dom%mat,dom%thrm,dom%bnd,time_now)
                !$ time2 = omp_get_wtime()
                !$ if(l_write_timer) print *,'TIME ydyn',time2-time1

                if (dom%par%pc_filter_vel) then 
                    
                    ! Modify ux/y_bar to use the average between the current and previous velocity solutions
                    dom%dyn%now%ux_bar = 0.5_wp*dom%dyn%now%ux_bar + 0.5_wp*dom%dyn%now%ux_bar_prev
                    dom%dyn%now%uy_bar = 0.5_wp*dom%dyn%now%uy_bar + 0.5_wp*dom%dyn%now%uy_bar_prev
                    
                end if 

                if (update_others_pc) then
                    ! Now, using old topography still, update additional fields.

                    ! Calculate material (ice properties, viscosity, etc.)
                    call calc_ymat(dom%mat,dom%tpo,dom%dyn,dom%thrm,dom%bnd,time_now)

                    ! Calculate thermodynamics (temperatures and enthalpy)
                    call calc_ytherm(dom%thrm,dom%tpo,dom%dyn,dom%mat,dom%bnd,time_now)
                end if 

                ! Step 3: Perform corrector step for topography
                ! Get corrected ice thickness and store it for later use
                
                !$ time1 = omp_get_wtime()
                ! Call corrector step for topography
                call calc_ytopo_pc(dom%tpo,dom%dyn,dom%mat,dom%thrm,dom%bnd,dom%dta,time_now,dom%tpo%par%topo_fixed,"corrector")
                !$ time2 = omp_get_wtime()
                !$ if(l_write_timer) print *,'TIME ytopo corrector',time2-time1

                ! Step 4: Determine truncation error for ice thickness

                !$ time1 = omp_get_wtime()

                if (.TRUE.) then 
                    ! not rk4...

                    select case(trim(dom%par%pc_method))
                        ! No default case necessary, handled earlier 

                        case("FE-SBE")
                        
                            ! FE-SBE truncation error 
                            call calc_pc_tau_fe_sbe(dom%time%pc_tau,dom%tpo%now%corr%H_ice,dom%tpo%now%pred%H_ice,dt_now)

                        case("AB-SAM")

                            ! AB-SAM truncation error 
                            call calc_pc_tau_ab_sam(dom%time%pc_tau,dom%tpo%now%corr%H_ice,dom%tpo%now%pred%H_ice,dt_now, &
                                                                                                dom%tpo%par%dt_zeta)

                        case("HEUN")

                            ! HEUN truncation error (same as FE-SBE)
                            call calc_pc_tau_heun(dom%time%pc_tau,dom%tpo%now%corr%H_ice,dom%tpo%now%pred%H_ice,dt_now)

                        case("RALSTON")

                            call calc_pc_tau_fe_sbe(dom%time%pc_tau,dom%tpo%now%corr%H_ice,dom%tpo%now%pred%H_ice,dt_now)

                    end select 

                else 
                    ! rk4 
                    dom%time%pc_tau = dom%tpo%rk4%tau 
                end if 

                ! Calculate eta for this timestep 
                call set_pc_mask(pc_mask,dom%time%pc_tau,dom%tpo%now%corr%H_ice,dom%tpo%now%pred%H_ice,dom%bnd%z_bed, &
                                dom%bnd%z_sl,dom%bnd%c%rho_ice,dom%bnd%c%rho_sw,dom%par%pc_eps, &
                                dom%tpo%par%boundaries,dom%tpo%par%margin_flt_subgrid)
                eta_now = calc_pc_eta(dom%time%pc_tau,mask=pc_mask)

                ! Save masked pc_tau for output too 
                dom%time%pc_tau_masked = dom%time%pc_tau 
                where( .not. pc_mask) dom%time%pc_tau_masked = 0.0_wp 

                ij = maxloc(abs(dom%time%pc_tau_masked))

                !write(*,"(a,f12.5,f12.5,f12.5,2i4,2f10.2)") &
                !    "test: ", time_now, dt_now, eta_now, ij(1), ij(2), &
                !    dom%tpo%now%pred%H_ice(ij(1),ij(2)), &
                !    dom%tpo%now%corr%H_ice(ij(1),ij(2))
                
                ! Check if this timestep should be rejected:
                ! If the redo iteration is not the last allowed and the timestep is still larger  
                ! than the minimum, then if eta > tolerance or checkerboard found in tau,
                ! then redo iteration: reject this timestep and try again with a smaller timestep
                if ( (iter_redo .lt. dom%par%pc_n_redo .and. dt_now .gt. dom%par%dt_min) &
                     .and. eta_now .gt. dom%par%pc_tol ) then

                    ! Calculate timestep reduction to apply
                    !rho_now = 0.7_wp
                    !rho_now = (dom%par%pc_tol / eta_now)
                    rho_now = 0.7_wp*(1.0_wp+(eta_now-dom%par%pc_tol)/10.0_wp)**(-1.0_wp) 

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
                !$ time2 = omp_get_wtime()
                !$ if(l_write_timer) print *,'TIME truncation error',time2-time1

            end do   ! End iteration loop 
            

            ! === Predictor-corrector completed successfully ===

            !$ time1 = omp_get_wtime()
            if (.not. update_others_pc) then
                ! Now, using old topography still, update additional fields.

                ! Calculate material (ice properties, viscosity, etc.)
                call calc_ymat(dom%mat,dom%tpo,dom%dyn,dom%thrm,dom%bnd,time_now)

                ! Calculate thermodynamics (temperatures and enthalpy)
                call calc_ytherm(dom%thrm,dom%tpo,dom%dyn,dom%mat,dom%bnd,time_now)

            end if 

            ! Update topography accounting for advective changes
            ! and mass balance changes and calving.

            call calc_ytopo_pc(dom%tpo,dom%dyn,dom%mat,dom%thrm,dom%bnd,dom%dta,time_now,dom%tpo%par%topo_fixed,"advance",use_H_pred=dom%par%pc_use_H_pred)

            ! Update time averaging of instantaneous rates
            call calc_ytopo_rates(dom%tpo,dom%bnd,time_now,dt_now,step="step",check_mb=check_mb)

            ! if (time_now .ge. 10.0) then 
            !     write(*,*) time_now
            !     write(*,*) maxval(dom%tpo%now%pred%H_ice)
            !     write(*,*) maxval(dom%tpo%now%corr%H_ice)
            !     write(*,*) maxval(dom%tpo%now%H_ice)
            !     stop 
            ! end if 
                
            ! === Done, Yelmo fields are fully consistent with time=time_now 

            ! Output some diagnostic info if desired
            if (very_verbose) then 

                write(*,"(a,f13.2,f10.2,G10.3,i5)") &
                    "yelmo:: times: ", &
                    time_now, dt_now, eta_now, dom%dyn%par%ssa_iter_now

            end if

            ! Output diagnostic file if desired 
            if (present(file_diagnostics)) then 
                ! Write step of continuous restart file 
                ! (initialized externally)
                call yelmo_restart_write(dom,file_diagnostics,time_now,init=.FALSE.)
            end if 

            ! Calculate model speed for this iteration
            call yelmo_cpu_time(cpu_time1)
            call yelmo_calc_speed(speed,model_time0,time_now,cpu_time0,cpu_time1)

            ! Collect how many times the redo-iteration loop had to run 
            ! (not counting the first pass, which is not a redo)
            iter_redo_tot = iter_redo_tot + (iter_redo-1) 
            
            ! Update dt and eta vectors for last N timesteps (first index becomes latest value)
            dom%time%pc_dt = cshift(dom%time%pc_dt,shift=-1)
            dom%time%pc_dt(1) = max(dt_now,dom%par%dt_min) 

            dom%time%pc_eta = cshift(dom%time%pc_eta,shift=-1)
            dom%time%pc_eta(1) = eta_now
            
            ! Activate pc method if not already active
            if (.not. dom%time%pc_active) dom%time%pc_active = .TRUE. 

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
                            dom%dyn%par%ssa_iter_now,iter_redo-1)
            
            end if 

            ! Make sure model is still running well
            call yelmo_check_kill(dom,time_now)

            !$ time2 = omp_get_wtime()
            !$ if(l_write_timer) print *,'TIME end',time2-time1

            ! Additionally check if minimum timestep is reached continuously
            ! when the minimum is set to a very small value.

            ! Set limit for check to be the last 50 timesteps or, if it is smaller,
            ! the total number of timesteps in this call of yelmo_update
            n_lim = min(50,nstep)

            if (n .ge. n_lim .and. dom%par%dt_min .le. 1e-2) then
                ! Currently n timesteps have been made at the limit of n_lim
                ! and the dt_min value is quite low. So the model may be stuck
                ! advancing in time. 
                
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

        if (dt_max_0 .gt. 0.0_wp) then
            ! Finalize averaging of instantaneous rates 
            ! (only if a timestep was calculated, otherwise maintain rates that were there)
            call calc_ytopo_rates(dom%tpo,dom%bnd,time,dt_max_0,step="final",overwrite=.TRUE.,check_mb=check_mb)
        end if

        ! Update regional calculations (for entire domain and subdomains)
        call yelmo_regions_update(dom)

        ! Compare with data 
        call ydata_compare(dom%dta,dom%tpo,dom%dyn,dom%mat,dom%thrm,dom%bnd,dom%par%domain)

        ! Write some diagnostics to make sure something useful is happening 
        if (yelmo_log) then

            n = count(dom%tpo%now%H_ice.gt.0.0)
            if (n .gt. 0.0) then 
                H_mean = sum(dom%tpo%now%H_ice,mask=dom%tpo%now%H_ice.gt.0.0)/real(n)
                T_mean = sum(dom%thrm%now%T_ice(:,:,dom%thrm%par%nz_aa),mask=dom%tpo%now%H_ice.gt.0.0)/real(n)
            else 
                H_mean = 0.0_wp 
                T_mean = 0.0_wp
            end if 

            n       = count(dt_save .ne. missing_value)
            n_dtmin = count( abs(dt_save(1:n)-dom%par%dt_min) .lt. dom%par%dt_min*1e-3 )

            if (n .gt. 0) then 
                max_dt_used = maxval(dt_save(1:n))
                min_dt_used = minval(dt_save(1:n))
            else 
                max_dt_used = 0.0 
                min_dt_used = 0.0 
            end if 

            write(*,"(a,f13.2,f10.2,f10.1,f8.1,2G10.3,1i6)") &
                        "yelmo:: timelog:", &
                            time_now, dom%time%model_speed, H_mean, T_mean,  &
                                            max_dt_used, min_dt_used, n_dtmin
            

            ! write(*,*) "time2: ", time, time_now, dom%tpo%par%time, dom%tpo%par%time_calv, &
            !                                     dom%thrm%par%time, dom%mat%par%time, dom%dyn%par%time


            ! Check mass conservation if desired (uncomment)
            ! call check_mass_conservation(dom%tpo%now%H_ice,dom%tpo%now%f_ice,dom%tpo%now%f_grnd,dom%tpo%now%dHidt, &
            !                 dom%tpo%now%mb_applied,dom%tpo%now%calv,dom%tpo%now%mb_dyn,dom%bnd%smb,dom%tpo%now%bmb, &
            !                 dom%tpo%now%fmb,dom%tpo%now%mb_resid,dom%grd%dx,dom%bnd%c%sec_year,time_now,dt_max_0, &
            !                 units="km^3/yr",label="final")

        end if 


        ! Finally, update z_bed relaxation rate to high resolution bedrock topography
        ! This rate should be passed to the isostasy module as needed.
        call yelmo_update_z_bed_restart_rate(dom,time)

        ! ! ajr: diagnostics 
        ! if (time .gt. 100.0 .and. dom%dyn%par%ssa_iter_now .ge. 5) then 
        !     stop 
        ! end if 

        return

    end subroutine yelmo_update

    subroutine yelmo_update_equil(dom,time,time_tot,dt,topo_fixed,tpo_solver,dyn_solver)
        ! Iterate yelmo solutions to equilibrate without updating boundary conditions

        type(yelmo_class), intent(INOUT) :: dom
        real(wp),          intent(IN)    :: time              ! [yr] Current time
        real(wp),          intent(IN)    :: time_tot          ! [yr] Equilibration time 
        real(wp),          intent(IN)    :: dt                ! Local dt to be used for all modules
        logical,           intent(IN)    :: topo_fixed        ! Should topography be fixed? 
        character(len=*),  intent(IN), optional :: tpo_solver
        character(len=*),  intent(IN), optional :: dyn_solver
        
        ! Local variables 
        type(yelmo_class) :: dom_ref 
        real(wp) :: time_now  
        integer  :: n, nstep 
        
        ! Only run equilibration if time_tot > 0 

        if (time_tot .gt. 0.0) then 

            ! Save original model configuration 
            dom_ref = dom 

            ! Set new, temporary parameter values from arguments
            dom%tpo%par%topo_fixed = topo_fixed 

            if (present(tpo_solver)) dom%tpo%par%solver = tpo_solver 
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
    
    subroutine yelmo_init(dom,filename,grid_def,time,load_topo,domain,grid_name,group)
        ! Initialize a yelmo domain, including the grid itself, 
        ! and all sub-components (topo,dyn,mat,therm,bound,data)

        !$ use omp_lib 

        implicit none

        type(yelmo_class) :: dom 
        character(len=*),  intent(IN) :: filename 
        character(len=*),  intent(IN) :: grid_def 
        real(wp),          intent(IN) :: time 
        logical, optional, intent(IN) :: load_topo 
        character(len=*),  intent(IN), optional :: domain
        character(len=*),  intent(IN), optional :: grid_name 
        character(len=*),  intent(IN), optional :: group

        ! Local variables
        integer :: n_threads 
        character(len=10) :: n_threads_str 
        character(len=32) :: nml_group
        
        ! Make sure we know the namelist group for the yelmo block
        if (present(group)) then
            nml_group = trim(group)
        else
            nml_group = "yelmo"         ! Default parameter blcok name
        end if

        ! ==== GLOBAL INIT CHECKS ============================================
        
        ! Check openmp status - set global variable to use as a switch 
        yelmo_use_omp = .FALSE. 
        !$ yelmo_use_omp = .TRUE.

        ! Output some information about openmp status 
        if (yelmo_use_omp) then 
            
            n_threads = 1
            !$ n_threads = omp_get_max_threads() 

            write(n_threads_str,"(i10)") n_threads 
            n_threads_str = adjustl(n_threads_str)

            write(*,*) "yelmo_global_init:: openmp is active, Yelmo will run on "//trim(n_threads_str)//" thread(s)."
            
        else 
            
            n_threads = 1
            write(*,*) "yelmo_global_init:: openmp is not active, Yelmo will run on 1 thread."

        end if 

        write(*,*) "yelmo_global_init:: yelmo_log = ", yelmo_log

        ! ==== END GLOBAL INIT CHECKS ============================================


        ! == yelmo == 

        ! Load the default yelmo parameters, then the domain specific parameters
        call yelmo_par_load(dom%par,filename,nml_group,domain,grid_name)
        
        ! Define physical constants
        call ybound_define_physical_constants(dom%bnd%c,dom%par%phys_const,domain,grid_name)

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

                write(io_unit_err,*) "yelmo_init:: Error: parameter grid_def not recognized: "//trim(grid_def)
                stop 

        end select 

        ! Check that grid has been defined properly 
        if (.not. allocated(dom%grd%x)) then 
            write(io_unit_err,*) "yelmo_init:: Error: ygrid has not been properly defined yet."
            write(io_unit_err,*) "(Either use yelmo_init_grid externally with desired grid parameters &
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

        call ytopo_par_load(dom%tpo%par,filename,dom%par%nml_ytopo,dom%par%nml_ycalv,dom%grd%nx,dom%grd%ny,dom%grd%dx,init=.TRUE.)

        call ytopo_alloc(dom%tpo%now,dom%tpo%par%nx,dom%tpo%par%ny)
        
        write(*,*) "yelmo_init:: topography initialized."
        
        ! == dynamics == 

        call ydyn_par_load(dom%dyn%par,filename,dom%par%nml_ydyn,dom%par%nml_ytill,dom%par%nml_yneff, &
                            dom%par%zeta_aa,dom%par%zeta_ac,dom%grd%nx,dom%grd%ny,dom%grd%dx,init=.TRUE.)

        call ydyn_alloc(dom%dyn%now,dom%dyn%par%nx,dom%dyn%par%ny,dom%dyn%par%nz_aa,dom%dyn%par%nz_ac)
        dom%dyn%par%init_state_set = .FALSE. 

        write(*,*) "yelmo_init:: dynamics initialized."
        
        ! == material == 

        call ymat_par_load(dom%mat%par,filename,dom%par%nml_ymat,dom%par%zeta_aa,dom%par%zeta_ac,dom%grd%nx,dom%grd%ny,dom%grd%dx,init=.TRUE.)

        call ymat_alloc(dom%mat%now,dom%mat%par%nx,dom%mat%par%ny,dom%mat%par%nz_aa,dom%mat%par%nz_ac,dom%mat%par%n_iso)
        
        write(*,*) "yelmo_init:: material initialized."
        
        ! == thermodynamics == 
        
        call ytherm_par_load(dom%thrm%par,filename,dom%par%nml_ytherm,dom%par%zeta_aa,dom%par%zeta_ac,dom%grd%nx,dom%grd%ny,dom%grd%dx,init=.TRUE.)

        call ytherm_alloc(dom%thrm%now,dom%thrm%par%nx,dom%thrm%par%ny,dom%thrm%par%nz_aa,dom%thrm%par%nz_ac,dom%thrm%par%nzr_aa)
        
        write(*,*) "yelmo_init:: thermodynamics initialized."
        
        ! === Yelmo IO tables ===
        
        ! Load variable io tables
        call load_var_io_table(dom%io%tpo,"input/yelmo-variables-ytopo.md")
        call load_var_io_table(dom%io%dyn,"input/yelmo-variables-ydyn.md")
        call load_var_io_table(dom%io%mat,"input/yelmo-variables-ymat.md")
        call load_var_io_table(dom%io%thrm,"input/yelmo-variables-ytherm.md")
        call load_var_io_table(dom%io%bnd,"input/yelmo-variables-ybound.md")
        call load_var_io_table(dom%io%dta,"input/yelmo-variables-ydata.md")

        write(*,*) "yelmo_init:: variable io tables loaded."
        
        ! === Ensure consistency with specific parameters ===

        ! For particular case related to eismint,
        ! ensure that bmb is not used in mass conservation or vertical velocity 
        dom%dyn%par%use_bmb = dom%tpo%par%use_bmb

        ! Modify grid boundary treatment according to the experiment parameter 
        select case(trim(dom%par%experiment))

            case("EISMINT")

                dom%tpo%par%boundaries  = "zeros"
                dom%dyn%par%boundaries  = "zeros"
                dom%thrm%par%boundaries = "zeros"

            case("MISMIP3D","MISMIP+") 

                dom%tpo%par%boundaries  = "MISMIP3D"
                dom%dyn%par%boundaries  = "MISMIP3D"
                dom%thrm%par%boundaries = "MISMIP3D"

                ! Consistency check 
                if (trim(dom%tpo%par%solver) .ne. "impl-lis") then 
                    write(io_unit_err,*) "yelmo_init:: Error: the mass conservation solver for MISMIP3D experiments &
                    &must be 'impl-lis' for stability. The 'expl' solver has not yet been designed to &
                    &handle ice advected at the border point nx-1, and thus oscillations can be produced. &
                    &Please set 'solver=impl-lis'."
                    stop 
                end if 
            
            case("TROUGH-F17")

                dom%tpo%par%boundaries  = "TROUGH"
                dom%dyn%par%boundaries  = "TROUGH"
                dom%thrm%par%boundaries = "TROUGH"

            case("SLAB")

                dom%tpo%par%boundaries  = "infinite"
                dom%dyn%par%boundaries  = "periodic"
                dom%thrm%par%boundaries = "periodic"
                
            case("ISMIPHOM","slab","periodic","periodic-xy") 
                ! Periodic boundary conditions in x and y, eg: X_1 = X_n-1; X_n = X_2

                dom%tpo%par%boundaries  = "periodic"
                dom%dyn%par%boundaries  = "periodic"
                dom%thrm%par%boundaries = "periodic"
            
            case("periodic-x") 
                ! Periodic boundary conditions in x-direction,
                ! infinite in y-direction
                dom%tpo%par%boundaries  = "periodic-x"
                dom%dyn%par%boundaries  = "periodic-x"
                dom%thrm%par%boundaries = "periodic-x"

            case("infinite") 
                ! Set border points equal to interior neighbors 

                dom%tpo%par%boundaries  = "infinite"
                dom%dyn%par%boundaries  = "infinite"
                dom%thrm%par%boundaries = "infinite"

            case DEFAULT
                ! zeros by default - safest option

                dom%tpo%par%boundaries  = "zeros"
                dom%dyn%par%boundaries  = "zeros"
                dom%thrm%par%boundaries = "zeros"

        end select 

        ! == boundary == 
        
        ! Allocate the yelmo data objects (arrays, etc)
        call ybound_alloc(dom%bnd,dom%grd%nx,dom%grd%ny)

        ! Load region/basin masks
        call ybound_load_masks(dom%bnd,filename,dom%par%nml_masks,dom%par%domain,dom%par%grid_name)
        
        ! Update the ice_allowed mask based on domain definition 
        call ybound_define_ice_allowed(dom%bnd,dom%par%domain)
        
        ! Define the advection mask (by default, for now allow advection everywhere)
        dom%tpo%now%mask_adv = 1

        write(*,*) "yelmo_init:: boundary initialized (loaded masks, set ref. topography)."
        
        ! == regions ==

        ! Initialize region: global domain
        ! Writing to file will only take place if user-program calls yelmo_regions_write(),
        ! which uses the flag specified below. Note currently this assumes outfldr="./" too.
        call yelmo_region_init(dom%reg,"global",mask=dom%bnd%ice_allowed,write_to_file=.TRUE.)

        ! Initialize regional averaging domains too (global region + zero subdomains for now)
        ! If regional subdomains are desired, this call will be made explicitly outside the program
        ! by the user. Ie,
        ! call yelmo_regions_init(dom,n=2)
        ! call yelmo_region_init(dom%reg(1),"Region1",mask1)
        ! call yelmo_region_init(dom%reg(2),"Region2",mask2)

        call yelmo_regions_init(dom,n=0)

        
        write(*,*) "yelmo_init:: regions initialized."
        
        ! == data == 
        
        call ydata_par_load(dom%dta%par,filename,dom%par%nml_data,dom%par%domain,dom%par%grid_name,init=.TRUE.)
        call ydata_alloc(dom%dta%pd,dom%grd%nx,dom%grd%ny,dom%par%nz_aa,dom%dta%par%pd_age_n_iso)

        ! Load data objects   
        call ydata_load(dom%dta,dom%bnd,filename,dom%tpo%par%grad_lim_zb,dom%grd%dx,dom%tpo%par%boundaries)

        ! Set H_ice_ref and z_bed_ref to present-day ice thickness by default 
        dom%bnd%H_ice_ref = dom%dta%pd%H_ice 
        dom%bnd%z_bed_ref = dom%dta%pd%z_bed

        write(*,*) "yelmo_init:: data intialized (loaded data if desired)."
        
        ! == topography ==

        ! Determine how to manage initial topography (H_ice,z_bed)
        call yelmo_init_topo(dom,filename,dom%par%nml_init_topo,time,load_topo)

        write(*,*) "yelmo_init:: topo intialized (loaded data if desired)."
        
        write(*,*) 
        write(*,*) "yelmo_init:: Initialization complete for domain: "// &
                   trim(dom%par%domain) 

        if (dom%par%log_timestep) then 
            ! Timestep file 
            call yelmo_timestep_write_init(dom%time%log_timestep_file,time,dom%grd%xc,dom%grd%yc,dom%par%pc_eps)
            call yelmo_timestep_write(dom%time%log_timestep_file,time,0.0_wp,0.0_wp,dom%time%pc_dt(1), &
                            dom%time%pc_eta(1),dom%time%pc_tau_masked,0.0_wp,0.0_wp,0.0_wp,dom%dyn%par%ssa_iter_now,0)
        end if 

        return

    end subroutine yelmo_init

    subroutine yelmo_init_topo(dom,filename,group,time,load_topo)
        ! This subroutine is the first step to intializing 
        ! the state variables. It initializes only the topography
        ! to facilitate calculation of boundary variables (eg, T_srf),
        ! which should be initialized externally afterwards.
        ! The state variables are either calculated directly or
        ! loaded from a restart file. 
            
        implicit none 

        type(yelmo_class), intent(INOUT) :: dom
        character(len=*),  intent(IN)    :: filename
        character(len=*),  intent(IN)    :: group       ! Usually "yelmo_init_topo"
        real(wp),          intent(IN)    :: time 
        logical, optional, intent(IN)    :: load_topo 

        ! Local variables 
        logical :: init_topo_load 
        character(len=1028) :: init_topo_path
        character(len=56)   :: init_topo_names(4)
        integer             :: init_topo_state
        real(wp)            :: z_bed_f_sd 
        real(wp)            :: smooth_H_ice
        real(wp)            :: smooth_z_bed

        real(wp), allocatable :: H_ice(:,:) 
        real(wp), allocatable :: z_bed(:,:) 
        real(wp), allocatable :: z_bed_sd(:,:) 
        real(wp), allocatable :: z_srf(:,:) 

        real(wp), allocatable :: dzb(:,:)
        real(wp), allocatable :: dzb_restart(:,:)
        
        ! Local copies of tpo and bnd and tme
        type(ytopo_class)  :: tpo_restart 
        type(ybound_class) :: bnd_restart 
        type(ytime_class)  :: tme_restart 

        ! Allocate local arrays
        allocate(H_ice(dom%grd%nx,dom%grd%ny))
        allocate(z_bed(dom%grd%nx,dom%grd%ny))
        allocate(z_bed_sd(dom%grd%nx,dom%grd%ny))
        allocate(z_srf(dom%grd%nx,dom%grd%ny))
        allocate(dzb(dom%grd%nx,dom%grd%ny))
        allocate(dzb_restart(dom%grd%nx,dom%grd%ny))
        
        ! Set to zero to start 
        H_ice    = 0.0_wp 
        z_bed    = 0.0_wp 
        z_bed_sd = 0.0_wp 
        z_srf    = 0.0_wp 

        ! Step 1: load topography variables from a file, if desired.
        ! Manipulate in local arrays, then store in the main dom object. 

        ! Load parameters related to topography initiaization 
        call nml_read(filename,group,"init_topo_load",  init_topo_load)
        call nml_read(filename,group,"init_topo_path",  init_topo_path)
        call nml_read(filename,group,"init_topo_names", init_topo_names)
        call nml_read(filename,group,"init_topo_state", init_topo_state)
        call nml_read(filename,group,"z_bed_f_sd",      z_bed_f_sd)
        call nml_read(filename,group,"smooth_H_ice",    smooth_H_ice)
        call nml_read(filename,group,"smooth_z_bed",    smooth_z_bed)

        call yelmo_parse_path(init_topo_path,dom%par%domain,dom%par%grid_name)
            
        ! Override parameter choice if command-line argument present 
        if (present(load_topo)) init_topo_load = load_topo 

        if (init_topo_load) then 
            ! =========================================
            ! Load topography data from netcdf file 

            call nc_read(init_topo_path,init_topo_names(1), H_ice, missing_value=mv)
            call nc_read(init_topo_path,init_topo_names(2), z_bed, missing_value=mv) 

            if (trim(init_topo_names(3)) .ne. ""     .and. &
                trim(init_topo_names(3)) .ne. "none" .and. &
                trim(init_topo_names(3)) .ne. "None") then 

                call nc_read(init_topo_path,init_topo_names(3),z_bed_sd)

                ! Apply scaling to adjust z_bed depending on standard deviation
                z_bed = z_bed + z_bed_f_sd*z_bed_sd 

            else
                z_bed_sd = 0.0_wp 
            end if 

            ! If desired and available, read surface elevation field
            ! too, in order to correct for englacial lakes.
            if (trim(init_topo_names(4)) .ne. ""     .and. &
                trim(init_topo_names(4)) .ne. "none" .and. &
                trim(init_topo_names(4)) .ne. "None") then 

                call nc_read(init_topo_path,init_topo_names(4),z_srf)

                ! Note: this routine uses z_sl, that is likely still set to zero
                ! here. This routine is mainly for fixing present-day datasets,
                ! so this is probably ok.
                call remove_englacial_lakes(H_ice,z_bed,z_srf,dom%bnd%z_sl,dom%bnd%c%rho_ice,dom%bnd%c%rho_sw)

                write(*,*) "yelmo_init_topo:: removed englacial lakes."

            end if

            ! Smooth ice thickness field, if desired 
            if (smooth_H_ice .ge. 1.0_wp) then 
                call smooth_gauss_2D(H_ice,dx=dom%grd%dx,f_sigma=smooth_H_ice)
            end if 

            ! Clean up ice thickness field 
            where (H_ice .lt. 0.1_wp) H_ice = 0.0_wp 
            where (H_ice .ge. 0.1_wp .and. H_ice .lt. 10.0_wp) H_ice = 10.0_wp 

            ! Smooth z_bed field, if desired 
            if (smooth_z_bed .ge. 1.0_wp) then 
                call smooth_gauss_2D(z_bed,dx=dom%grd%dx,f_sigma=smooth_z_bed)
            end if 

            ! Adjust bedrock topography and ice thickness for smoothness
            call adjust_topography_gradients(z_bed,H_ice,dom%tpo%par%grad_lim_zb,dom%grd%dx,dom%tpo%par%boundaries)

            ! Additionally modify initial topographic state 
            select case(init_topo_state)

                case(0) 

                    ! Pass, use topography as loaded 

                case(1) 
                    ! Remove ice, but do not adjust bedrock 

                    H_ice = 0.0_wp 

                case(2)
                    ! Remove ice, set bedrock to isostatically rebounded state 

                    z_bed = z_bed + (dom%bnd%c%rho_ice/dom%bnd%c%rho_a)*H_ice
                    H_ice = 0.0_wp

                case DEFAULT 

                    write(io_unit_err,*) "yelmo_init_topo:: Error: init_topo_state choice not recognized."
                    write(io_unit_err,*) "init_topo_state = ", init_topo_state 
                    stop 

            end select

            ! Store in Yelmo object 
            dom%tpo%now%H_ice = H_ice 
            dom%bnd%z_bed     = z_bed 
            dom%bnd%z_bed_sd  = z_bed_sd 
                
            ! Finally, calculate and apply all additional (generally artificial) ice thickness adjustments
            ! Set minimum ice thickness to 1m for safety to start.
            call calc_G_boundaries(dom%tpo%now%mb_resid,dom%tpo%now%H_ice,dom%tpo%now%f_ice,dom%tpo%now%f_grnd, &
                                                dom%dyn%now%uxy_b,dom%bnd%ice_allowed,dom%tpo%par%boundaries,dom%bnd%H_ice_ref, &
                                                H_min_flt=1.0_wp,H_min_grnd=1.0_wp,dt=1.0_wp)
            ! Apply rate and update ice thickness
            call apply_tendency(dom%tpo%now%H_ice,dom%tpo%now%mb_resid,dt=1.0_wp,label="init")
        end if 

        ! Step 2: load topo and bnd variables from a restart file if desired 
        ! Store fields in temporary objects and determine which fields to pass
        ! to the main dom object.

        if (dom%par%use_restart) then 
            ! Load variables from a restart file. Note: this will
            ! overwrite all information stored in yelmo object from above.

            ! Intialize tpo and bnd and tme objects locally 
            tpo_restart = dom%tpo 
            bnd_restart = dom%bnd 
            tme_restart = dom%time 

            call yelmo_restart_read_topo_bnd(tpo_restart,bnd_restart,tme_restart,dom%par%restart_interpolated, &
                                             dom%grd,dom%par%domain,dom%par%grid_name,dom%par%restart,time)
            ! Now determine which values should be used from restart.
            ! Replace fields in restart objects that should not be used. 

            if (.not. dom%par%restart_H_ice) then 
                ! Use field from default initialization
                tpo_restart%now%H_ice = dom%tpo%now%H_ice

                write(*,*) "yelmo_init_topo: H_ice taken from input file, not restart file."
                write(*,*) "yelmo.restart: ", trim(dom%par%restart)
                write(*,*) "yelmo_init_topo.init_topo_path: ", trim(init_topo_path)
            end if 

            if ( (.not. dom%par%restart_z_bed) .or. dom%par%restart_interpolated .eq. 1) then 
                ! Use fields from default initialization, but make sure the current
                ! state of isostatic rebound is well represented. 

                ! Determine isostatic offset in each case 
                dzb         = dom%bnd%z_bed - dom%bnd%z_bed_ref 
                dzb_restart = bnd_restart%z_bed - bnd_restart%z_bed_ref 
                
                ! Set the restart field equal to that loaded via parameters but offset
                ! with the correct offset from the restart file.

                ! Modify main z_bed object to account for isostatic uplift of restart field
                dom%bnd%z_bed = dom%bnd%z_bed - dzb + dzb_restart
                
                ! Determine remaining differences between (z_bed from parameters)
                ! and (z_bed from restart file) - these differences should be equivalent
                ! to the high resolution information not contained in the low-to-high
                ! resolution field from the restart file. 
                bnd_restart%z_bed_corr = dom%bnd%z_bed - bnd_restart%z_bed
                bnd_restart%restart_relax_init = time
                
                if (dom%par%restart_relax .eq. 0.0) then
                    ! Impose high resolution changes to new topography directly

                    bnd_restart%z_bed       = dom%bnd%z_bed
                    bnd_restart%dzbdt_corr  = 0.0 
                else 

                    ! Calculate the desired rate of change based on relaxation time
                    bnd_restart%dzbdt_corr  = bnd_restart%z_bed_corr / dom%par%restart_relax

                    ! Note, now, do nothing: do not modify bnd_restart%z_bed.
                    ! So, to start with, the bedrock topography will still be fully
                    ! consistent with the simulation being loaded from the restart file. 
                    ! Pass dzbdt_corr to isostasy routine to slow incorporate
                    ! high-resolution information after initializing all other fields. 
                
                end if 

                ! Finally, store the variability field loaded from parameter choices too
                bnd_restart%z_bed_sd = dom%bnd%z_bed_sd

                write(*,*) "yelmo_init_topo: z_bed taken from input file, not restart file. But it has been &
                            &corrected to reflect isostatic offset from z_bed_ref in restart file."
                write(*,*) "yelmo.restart: ", trim(dom%par%restart)
                write(*,*) "yelmo_init_topo.init_topo_path: ", trim(init_topo_path)

                
            end if

            ! Replace several boundary fields like masks from the main dom object that 
            ! were loaded via parameter file choices. These should not be loaded 
            ! from the restart file, especially if the restart file was at, e.g., 
            ! a lower resolution 
            bnd_restart%ice_allowed = dom%bnd%ice_allowed 
            bnd_restart%calv_mask   = dom%bnd%calv_mask 
            bnd_restart%H_ice_ref   = dom%bnd%H_ice_ref 
            bnd_restart%z_bed_ref   = dom%bnd%z_bed_ref 

            bnd_restart%domain_mask = dom%bnd%domain_mask 
            bnd_restart%tau_relax   = dom%bnd%tau_relax 

            bnd_restart%basins      = dom%bnd%basins 
            bnd_restart%basin_mask  = dom%bnd%basin_mask 
            bnd_restart%regions     = dom%bnd%regions 
            bnd_restart%region_mask = dom%bnd%region_mask 
            
            ! Finally populate the main dom object with the desired restart fields
            dom%tpo  = tpo_restart 
            dom%bnd  = bnd_restart 
            dom%time = tme_restart

            ! And ensure pc is already active
            dom%time%pc_active = .TRUE. 

        end if 

        ! Step 3: update remaining topogaphic info to be consistent with initial fields 
        ! Here several fields in ytopo will be overwritten and topo will be fully consistent with itself.

        ! Run topo and masks to make sure all fields are synchronized (masks, etc)
        ! Note: thermodynamic state has not been loaded yet, so mask_bed produced here
        ! will not contain regions of temperate ice. masks should be updated again
        ! after loaded remaining fields.
        !call calc_ytopo_rk4(dom%tpo,dom%dyn,dom%mat,dom%thrm,dom%bnd,time,topo_fixed=.TRUE.)
        call calc_ytopo_pc(dom%tpo,dom%dyn,dom%mat,dom%thrm,dom%bnd,dom%dta,time,topo_fixed=.TRUE.,pc_step="none",use_H_pred=dom%par%pc_use_H_pred)

        ! Update regional calculations (for entire domain and subdomains)
        call yelmo_regions_update(dom)

        ! Summary for log file: 

        write(*,*) "yelmo_init_topo:: range(z_bed):     ", minval(dom%bnd%z_bed),     maxval(dom%bnd%z_bed)
        write(*,*) "yelmo_init_topo:: range(z_bed_sd):  ", minval(dom%bnd%z_bed_sd),  maxval(dom%bnd%z_bed_sd)
        write(*,*) "yelmo_init_topo:: range(H_ice):     ", minval(dom%tpo%now%H_ice), maxval(dom%tpo%now%H_ice) 
        write(*,*) "yelmo_init_topo:: scaling fac z_bed_f_sd: ", z_bed_f_sd  
        
        ! ajr: diagnostics!!
        !call yelmo_restart_write(dom,"./yelmo_check_z_bed.nc",time=0.0_wp,init=.TRUE.)
        !stop 

        ! Finally lets initialize the LSF mask
        call LSFinit(dom%tpo%now%lsf,dom%tpo%now%H_ice,dom%bnd%z_bed,dom%tpo%par%dx)

        return 

    end subroutine yelmo_init_topo

    subroutine yelmo_update_z_bed_restart_rate(dom,time)
        
        implicit none 

        type(yelmo_class), intent(INOUT) :: dom 
        real(wp), intent(IN) :: time 

        ! Local variables
        real(wp) :: time_elapsed 

        time_elapsed = time - dom%bnd%restart_relax_init

        if (time_elapsed .gt. dom%par%restart_relax) then 
            dom%bnd%dzbdt_corr = 0.0 
        end if 

        if (maxval(abs(dom%bnd%dzbdt_corr)) .gt. 0.0) then 
            write(*,*) "restart z_bed_corr rate: ", time, time_elapsed, &
                                        minval(dom%bnd%dzbdt_corr), maxval(dom%bnd%dzbdt_corr)
        end if 

        return

    end subroutine yelmo_update_z_bed_restart_rate

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
        real(wp),          intent(IN)    :: time  
        character(len=*),  intent(IN)    :: thrm_method 

        ! Local variables 
        integer :: q 
        character(len=256) :: dom_thrm_method 
        character(len=256) :: dom_thrm_rock_method 
        
        ! Initialize variables

        if (dom%par%use_restart) then 
            ! Load variables from a restart file 

            call yelmo_restart_read(dom,trim(dom%par%restart),time)
            
            ! And ensure pc is already active
            dom%time%pc_active = .TRUE. 
            
        else 

            ! Consistency check 
            if (trim(thrm_method) .ne. "linear" .and. trim(thrm_method) .ne. "robin" &
                .and. trim(thrm_method) .ne. "robin-cold") then 
                write(io_unit_err,*) "yelmo_init_state:: Error: temperature initialization must be &
                           &'linear', 'robin' or 'robin-cold' in order to properly prescribe &
                           &initial temperatures."
                stop 
            end if
            
            ! Store original model choices locally 
            dom_thrm_method      = dom%thrm%par%method 
            dom_thrm_rock_method = dom%thrm%par%rock_method 

            ! Impose initialization choices 
            dom%thrm%par%method      = thrm_method 
            dom%thrm%par%rock_method = "equil" 

            ! Run topo and masks to make sure all fields are synchronized (masks, etc)
            !call calc_ytopo_rk4(dom%tpo,dom%dyn,dom%mat,dom%thrm,dom%bnd,time,topo_fixed=.TRUE.)
            call calc_ytopo_pc(dom%tpo,dom%dyn,dom%mat,dom%thrm,dom%bnd,dom%dta,time,topo_fixed=.TRUE.,pc_step="none",use_H_pred=dom%par%pc_use_H_pred)

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
            ! This will be overwritten when cb_ref/beta are calculated internally
            if (maxval(dom%dyn%now%beta) .eq. 0.0_wp) then 
                if (maxval(dom%dyn%now%cb_ref) .eq. 0.0_wp) then 
                    dom%dyn%now%cb_ref = 1.0
                end if 
                dom%dyn%now%c_bed  = dom%dyn%now%cb_ref*1e5
                dom%dyn%now%beta   = dom%dyn%now%c_bed
            end if
            
            call calc_ydyn(dom%dyn,dom%tpo,dom%mat,dom%thrm,dom%bnd,time)

            ! Calculate material information again with updated dynamics
        
            call calc_ymat(dom%mat,dom%tpo,dom%dyn,dom%thrm,dom%bnd,time)

            ! Restore original model choices 
            ! Impose initialization choices 
            dom%thrm%par%method      = dom_thrm_method 
            dom%thrm%par%rock_method = dom_thrm_rock_method

        end if 

        ! Re-run topo again to make sure all fields are synchronized (masks, etc)
        !call calc_ytopo_rk4(dom%tpo,dom%dyn,dom%mat,dom%thrm,dom%bnd,time,topo_fixed=.TRUE.)
        call calc_ytopo_pc(dom%tpo,dom%dyn,dom%mat,dom%thrm,dom%bnd,dom%dta,time,topo_fixed=.TRUE.,pc_step="none",use_H_pred=dom%par%pc_use_H_pred)

        ! Update regional calculations (for entire domain and subdomains)
        call yelmo_regions_update(dom)

        return 

    end subroutine yelmo_init_state

    subroutine yelmo_par_load(par,filename,group,domain,grid_name)

        type(yelmo_param_class), intent(OUT) :: par
        character(len=*),        intent(IN)  :: filename
        character(len=*),        intent(IN)  :: group               ! Usually "yelmo"
        character(len=*),        intent(IN), optional :: domain
        character(len=*),        intent(IN), optional :: grid_name 

        call nml_read(filename,group,"domain",        par%domain)
        call nml_read(filename,group,"grid_name",     par%grid_name)
        call nml_read(filename,group,"grid_path",     par%grid_path)
        call nml_read(filename,group,"phys_const",    par%phys_const)
        call nml_read(filename,group,"experiment",    par%experiment)

        call nml_read(filename,group,"nml_ytopo",     par%nml_ytopo)
        call nml_read(filename,group,"nml_ycalv",     par%nml_ycalv)
        call nml_read(filename,group,"nml_ydyn",      par%nml_ydyn)
        call nml_read(filename,group,"nml_ytill",     par%nml_ytill)
        call nml_read(filename,group,"nml_yneff",     par%nml_yneff)
        call nml_read(filename,group,"nml_ymat",      par%nml_ymat)
        call nml_read(filename,group,"nml_ytherm",    par%nml_ytherm)
        call nml_read(filename,group,"nml_masks",     par%nml_masks)
        call nml_read(filename,group,"nml_init_topo", par%nml_init_topo)
        call nml_read(filename,group,"nml_data",      par%nml_data)
        
        call nml_read(filename,group,"restart",       par%restart)
        call nml_read(filename,group,"restart_z_bed", par%restart_z_bed)
        call nml_read(filename,group,"restart_H_ice", par%restart_H_ice)
        call nml_read(filename,group,"restart_relax", par%restart_relax)
        call nml_read(filename,group,"log_timestep",  par%log_timestep)
        call nml_read(filename,group,"disable_kill",  par%disable_kill)
        call nml_read(filename,group,"zeta_scale",    par%zeta_scale)
        call nml_read(filename,group,"zeta_exp",      par%zeta_exp)
        call nml_read(filename,group,"nz_aa",         par%nz_aa)
        call nml_read(filename,group,"dt_method",     par%dt_method)
        call nml_read(filename,group,"dt_min",        par%dt_min)
        call nml_read(filename,group,"cfl_max",       par%cfl_max)
        call nml_read(filename,group,"cfl_diff_max",  par%cfl_diff_max)
        call nml_read(filename,group,"pc_method",     par%pc_method)
        call nml_read(filename,group,"pc_controller", par%pc_controller)
        call nml_read(filename,group,"pc_use_H_pred", par%pc_use_H_pred)
        call nml_read(filename,group,"pc_filter_vel", par%pc_filter_vel)
        call nml_read(filename,group,"pc_corr_vel",   par%pc_corr_vel)
        call nml_read(filename,group,"pc_n_redo",     par%pc_n_redo)
        call nml_read(filename,group,"pc_tol",        par%pc_tol)
        call nml_read(filename,group,"pc_eps",        par%pc_eps)
        
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
            write(io_unit_err,*) "yelmo_par_load:: error: pc_eps must be greater than pc_tol."
            write(io_unit_err,*) trim(filename), " : ", trim(group)
            write(io_unit_err,*) "pc_eps, pc_tol: ", par%pc_eps, par%pc_tol 
            stop 
        end if

        return

    end subroutine yelmo_par_load
    
    subroutine yelmo_end(dom,time)
        ! Deallocate yelmo objects 

        implicit none 

        type(yelmo_class), intent(INOUT) :: dom 
        real(wp), intent(IN) :: time 

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
        real(wp), intent(IN) :: time
        
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
        real(wp),          intent(IN) :: time
        character(len=*), optional, intent(IN) :: kill_request 

        ! Local variables 
        integer :: i, j, k 
        logical :: kill_it, kill_it_H, kill_it_vel, kill_it_temp, kill_it_nan, kill_it_eta   
        character(len=512) :: kill_msg 
        real(wp) :: pc_eta_avg 
        character(len=3) :: pc_iter_str(10) 

        real(wp), parameter :: H_lim = 1e4   ! [m] 
        real(wp), parameter :: u_lim = 1e4   ! [m/a]

        kill_it_H    = .FALSE.
        kill_it_vel  = .FALSE.
        kill_it_temp = .FALSE. 
        kill_it_nan  = .FALSE. 
        kill_it_eta  = .FALSE. 

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

        if (maxval(abs(dom%thrm%now%T_ice-dom%thrm%now%T_ice)) .ne. 0.0 ) then 

            kill_it_temp = .TRUE. 
            kill_msg     = "Temperature field invalid."

        end if

        ! Additionally check for NANs using intrinsic ieee_arithmetic module 
        do j = 1, dom%grd%ny 
        do i = 1, dom%grd%nx 
            
            if (ieee_is_nan(dom%dyn%now%uxy_bar(i,j))) then 
                kill_it_nan = .TRUE. 
                write(kill_msg,*) "** NANs detected - uxy_bar ** ... i, j: ", i, j 
                exit 
            end if 

            if (ieee_is_nan(dom%tpo%now%H_ice(i,j))) then 
                kill_it_nan = .TRUE. 
                write(kill_msg,*) "** NANs detected - H_ice ** ... i, j: ", i, j 
                exit 
            end if 

            do k = 1, dom%thrm%par%nz_aa

                if (ieee_is_nan(dom%thrm%now%T_ice(i,j,k))) then 
                    kill_it_nan = .TRUE. 
                    write(kill_msg,*) "** NANs detected - T_ice ** ... i, j, k: ", i, j, k
                    exit 
                end if 

            end do

        end do 
        end do 

        pc_eta_avg = sum(dom%time%pc_eta) / real(size(dom%time%pc_eta,1),prec) 

        if (pc_eta_avg .gt. 10.0*dom%par%pc_tol) then 
            kill_it_eta = .TRUE. 
            write(kill_msg,"(a,g12.4,a,10g12.4)") "mean[pc_eta] > [10*pc_tol]: pc_eta_avg = ", pc_eta_avg, &
                                                                                " | pc_eta: ", dom%time%pc_eta
        end if 

        ! Determine if model should be killed 
        kill_it = kill_it_H .or. kill_it_vel .or. kill_it_temp .or. kill_it_nan .or. kill_it_eta 

        ! Definitely kill the model if it was requested externally
        if (present(kill_request)) then 
            kill_it = .TRUE. 
            kill_msg = trim(kill_request)
        end if 

        if (kill_it .and. (.not. dom%par%disable_kill)) then 
            ! Model is not running properly, kill it. 

            write(io_unit_err,*) 
            write(io_unit_err,*) 
            write(io_unit_err,"(a)") "yelmo_check_kill:: Error: model is not running properly:"
            write(io_unit_err,"(a)") trim(kill_msg) 
            write(io_unit_err,*) 
            write(io_unit_err,"(a11,f15.3)")  "timestep    = ", time
            write(io_unit_err,*) 
            write(io_unit_err,"(a,2g12.4)")   "pc_eps, tol = ", dom%par%pc_eps, dom%par%pc_tol 
            write(io_unit_err,"(a,g12.4)")    "pc_eta_avg  = ", pc_eta_avg
            
            write(io_unit_err,"(a4,1x,2a12)") "iter", "pc_dt", "pc_eta"
            do k = 1, size(dom%time%pc_eta,1)
                write(io_unit_err,"(a4,1x,2g12.4)") trim(pc_iter_str(k)), dom%time%pc_dt(k), dom%time%pc_eta(k) 
            end do 

            write(io_unit_err,*) 
            write(io_unit_err,"(a16,2g14.4)") "range(H_ice):   ", minval(dom%tpo%now%H_ice),   maxval(dom%tpo%now%H_ice)
            write(io_unit_err,"(a16,2g14.4)") "range(uxy_bar): ", minval(dom%dyn%now%uxy_bar), maxval(dom%dyn%now%uxy_bar)
            write(io_unit_err,"(a16,2g14.4)") "range(T_ice):   ", minval(dom%thrm%now%T_ice),  maxval(dom%thrm%now%T_ice)
            write(io_unit_err,*) 

            call yelmo_restart_write(dom,"yelmo_killed.nc",time=time) 

            write(io_unit_err,*) "Restart file written: "//"yelmo_killed.nc"
            write(io_unit_err,*) 
            write(io_unit_err,"(a,f15.3,a)") "Time =", time, ": stopping model (killed)." 
            write(io_unit_err,*) 

            stop "yelmo_check_kill error, see log."

        end if 

        return 

    end subroutine yelmo_check_kill

end module yelmo_ice


