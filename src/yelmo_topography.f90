
module yelmo_topography

    use nml  
    use ncio
    
    use yelmo_defs
    use yelmo_tools 
    
    use mass_conservation
    use calving
    use topography 

    implicit none
    
    ! Key for matching bed types given by mask_bed 
    integer, parameter :: mask_bed_ocean   = 0 
    integer, parameter :: mask_bed_land    = 1
    integer, parameter :: mask_bed_frozen  = 2
    integer, parameter :: mask_bed_stream  = 3
    integer, parameter :: mask_bed_grline  = 4
    integer, parameter :: mask_bed_float   = 5
    integer, parameter :: mask_bed_island  = 6
    integer, parameter :: mask_bed_partial = 7

    private
    public :: calc_ytopo
    public :: ytopo_par_load, ytopo_alloc, ytopo_dealloc
     
    ! Integers
    public :: mask_bed_ocean  
    public :: mask_bed_land  
    public :: mask_bed_frozen
    public :: mask_bed_stream
    public :: mask_bed_grline
    public :: mask_bed_float 
    public :: mask_bed_island
    
contains

    subroutine calc_ytopo(tpo,dyn,mat,thrm,bnd,time,topo_fixed)
        ! Calculate adjustments to surface elevation, bedrock elevation
        ! and ice thickness 

        implicit none 

        type(ytopo_class),  intent(INOUT) :: tpo
        type(ydyn_class),   intent(IN)    :: dyn
        type(ymat_class),   intent(IN)    :: mat
        type(ytherm_class), intent(IN)    :: thrm  
        type(ybound_class), intent(IN)    :: bnd 
        real(prec),         intent(IN)    :: time  
        logical,            intent(IN)    :: topo_fixed 

        ! Local variables 
        real(prec) :: dx, dt, dt_calv   
        integer :: i, j, nx, ny  
        real(prec), allocatable :: mbal(:,:) 
        real(prec), allocatable :: calv_sd(:,:) 
        logical :: reset_mb_resid 
        logical :: call_calving 

        real(8)    :: cpu_time0, cpu_time1
        real(prec) :: model_time0, model_time1 
        real(prec) :: speed 

        real(wp), parameter :: dt_calv_min = 0.0_wp 

        nx = size(tpo%now%H_ice,1)
        ny = size(tpo%now%H_ice,2)

        allocate(mbal(nx,ny))
        allocate(calv_sd(nx,ny))

        ! Initialize time if necessary 
        if (tpo%par%time .gt. dble(time)) then 
            tpo%par%time = dble(time) 
        end if 
        
        if (tpo%par%time_calv .gt. dble(time)) then 
            tpo%par%time_calv = dble(time)  
        end if 

        ! Get time step
        dt      = dble(time) - tpo%par%time 
        dt_calv = dble(time) - tpo%par%time_calv 

        call_calving = .FALSE. 
        if (dt_calv .ge. dt_calv_min) call_calving = .TRUE. 

        ! Store initial cpu time and model time for metrics later
        call yelmo_cpu_time(cpu_time0)
        model_time0 = tpo%par%time 

        ! Combine basal mass balance into one field accounting for 
        ! grounded/floating fraction of grid cells 
        call calc_bmb_total(tpo%now%bmb,thrm%now%bmb_grnd,bnd%bmb_shlf,tpo%now%H_grnd, &
                            tpo%now%f_grnd,tpo%par%bmb_gl_method,tpo%par%diffuse_bmb_shlf)
        
        ! Combine frontal mass balance into one field, and 
        ! calculate as needed 
        call calc_fmb_total(tpo%now%fmb,bnd%fmb_shlf,bnd%bmb_shlf,tpo%now%H_ice, &
                        tpo%now%H_grnd,tpo%now%f_ice,tpo%par%fmb_method,tpo%par%fmb_scale,tpo%par%dx)

        ! 1. Perform topography calculations ------------------

        if ( .not. topo_fixed .and. dt .gt. 0.0 ) then 

            ! === Step 0: define temporary variable for total column mass balance (without calving)
           
            mbal = bnd%smb + tpo%now%bmb + tpo%now%fmb 
            
            if (.not. tpo%par%use_bmb) then
                ! WHEN RUNNING EISMINT1 ensure bmb and fmb are not accounted for here !!!
                mbal = bnd%smb  
            end if 

            ! === Step 1: ice thickness evolution from dynamics alone ===

            ! Calculate the ice thickness conservation from dynamics only -----
            call calc_ice_thickness_dyn(tpo%now%H_ice,tpo%now%dHdt_n,tpo%now%H_ice_n,tpo%now%H_ice_pred, &
                                        tpo%now%f_ice,tpo%now%f_grnd,dyn%now%ux_bar,dyn%now%uy_bar, &
                                        solver=tpo%par%solver,dx=tpo%par%dx,dt=dt,beta=tpo%par%dt_beta,pc_step=tpo%par%pc_step) 

            ! Update ice fraction mask 
            call calc_ice_fraction(tpo%now%f_ice,tpo%now%H_ice,bnd%z_bed,bnd%z_sl,tpo%par%margin_flt_subgrid)
            
            ! === Step 2: ice thickness evolution from vertical column mass balance ===

            ! Apply mass-conservation step (mbal)
            call calc_ice_thickness_mbal(tpo%now%H_ice,tpo%now%f_ice,tpo%now%mb_applied, &
                                         tpo%now%f_grnd,bnd%z_sl-bnd%z_bed,mbal,tpo%par%dx,dt)

            ! If subgrid ice margin is not being used, then tiny ice 
            ! thicknesses should be removed before calving step.             
            if (.not. tpo%par%margin_flt_subgrid) then 

                call apply_ice_thickness_boundaries(tpo%now%mb_resid,tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd, &
                                                dyn%now%uxy_b,bnd%ice_allowed,tpo%par%boundaries,bnd%H_ice_ref, &
                                                H_min_flt=5.0_wp,H_min_grnd=0.0_wp,dt=dt,reset=.TRUE.)

                reset_mb_resid = .FALSE.
            else
                reset_mb_resid = .TRUE. 
            end if 

            ! Update ice fraction mask 
            call calc_ice_fraction(tpo%now%f_ice,tpo%now%H_ice,bnd%z_bed,bnd%z_sl,tpo%par%margin_flt_subgrid)
            
            ! === Step 3: ice thickness evolution from calving ===
            
            ! Diagnose strains and stresses relevant to calving 

            ! eps_eff = effective strain = eigencalving e+*e- following Levermann et al. (2012)
            call calc_eps_eff(tpo%now%eps_eff,mat%now%strn2D%eps_eig_1,mat%now%strn2D%eps_eig_2,tpo%now%f_ice)

            ! tau_eff = effective stress ~ von Mises stress following Lipscomb et al. (2019)
            call calc_tau_eff(tpo%now%tau_eff,mat%now%strs2D%tau_eig_1,mat%now%strs2D%tau_eig_2,tpo%now%f_ice,tpo%par%w2)

            ! For now, mask eps_eff and tau_eff to floating points 
            ! for easier visualization. If quantities will be used
            ! for grounded ice in the future, then this masking can 
            ! be removed. 
            where (tpo%now%f_grnd .eq. 1.0_wp) tpo%now%eps_eff = 0.0_wp 
            where (tpo%now%f_grnd .eq. 1.0_wp) tpo%now%tau_eff = 0.0_wp 
            
            ! Diagnose potential floating-ice calving rate [m/yr]

if (.TRUE.) then 
            select case(trim(tpo%par%calv_flt_method))
                
                case("vm-l19","eigen")
                    ! Consistency check

                    if (.not. tpo%par%margin_flt_subgrid) then 

                        write(io_unit_err,*) ""
                        write(io_unit_err,*) ""
                        write(io_unit_err,*) "calv_flt_method=['vm-l19','eigen'] must be used with margin_flt_subgrid=True."
                        write(io_unit_err,*) "calv_flt_method    = ", trim(tpo%par%calv_flt_method)
                        write(io_unit_err,*) "margin_flt_subgrid = ", tpo%par%margin_flt_subgrid
                        stop "Program stopped."

                    end if 
                    
            end select 
end if
            
            select case(trim(tpo%par%calv_flt_method))

                case("zero","none")

                    tpo%now%calv_flt = 0.0 

                case("simple") 
                    ! Use simple threshold method

                    call calc_calving_rate_simple(tpo%now%calv_flt,tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd, &
                                                    tpo%par%calv_H_lim,tpo%par%calv_tau)
                    
                case("flux") 
                    ! Use threshold+flux method from Peyaud et al. (2007), ie, GRISLI,
                    ! but reformulated. 

                    call calc_calving_rate_flux(tpo%now%calv_flt,tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd,mbal, &
                                                dyn%now%ux_bar,dyn%now%uy_bar,tpo%par%dx,tpo%par%calv_H_lim,tpo%par%calv_tau)
                    
                case("flux-grisli")
                    ! Use threshold+flux method from Peyaud et al. (2007), ie, GRISLI

                    call calc_calving_rate_flux_grisli(tpo%now%calv_flt,tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd,mbal, &
                                                dyn%now%ux_bar,dyn%now%uy_bar,tpo%par%dx,tpo%par%calv_H_lim,tpo%par%calv_tau)
                    
                case("vm-l19")
                    ! Use von Mises calving as defined by Lipscomb et al. (2019)

                    ! Next, diagnose calving
                    call calc_calving_rate_vonmises_l19(tpo%now%calv_flt,tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd, &
                                                                                tpo%now%tau_eff,tpo%par%dx,tpo%par%kt)
                case("eigen")
                    ! Use Eigen calving as defined by Levermann et al. (2012)

                    ! Next, diagnose calving
                    call calc_calving_rate_eigen(tpo%now%calv_flt,tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd, &
                                                                                tpo%now%eps_eff,tpo%par%dx,tpo%par%k2)

                case("kill") 
                    ! Delete all floating ice (using characteristic time parameter)
                    call calc_calving_rate_kill(tpo%now%calv_flt,tpo%now%H_ice,tpo%now%f_grnd.eq.0.0_prec,tpo%par%calv_tau,dt)

                case("kill-pos")
                    ! Delete all floating ice beyond a given location (using characteristic time parameter)

                    call calc_calving_rate_kill(tpo%now%calv_flt,tpo%now%H_ice, &
                                                    ( tpo%now%f_grnd .eq. 0.0_wp .and. &
                                                      tpo%now%H_ice  .gt. 0.0_wp .and. &
                                                      bnd%calv_mask ), tau=0.0_wp, dt=dt )

                case DEFAULT 

                    write(*,*) "calc_ytopo:: Error: floating calving method not recognized."
                    write(*,*) "calv_flt_method = ", trim(tpo%par%calv_flt_method)
                    stop 

            end select
            
            select case(trim(tpo%par%calv_flt_method))

                case("vm-l19","eigen")
                    
                    ! Scale calving with 'thin' calving rate to ensure 
                    ! small ice thicknesses are removed.
                    call apply_thin_calving_rate(tpo%now%calv_flt,tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd,tpo%par%calv_thin)

                    ! Adjust calving so that any excess is distributed to upstream neighbors
                    call calc_calving_residual(tpo%now%calv_flt,tpo%now%H_ice,tpo%now%f_ice,dt)
            
            end select 

            ! Additionally ensure higher calving rate for floating tongues of
            ! one grid-point width.
            call calc_calving_tongues(tpo%now%calv_flt,tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd,tpo%par%calv_tau)
            
            ! Diagnose potential grounded-ice calving rate [m/yr]

            select case(trim(tpo%par%calv_grnd_method))

                case("zero","none")

                    tpo%now%calv_grnd = 0.0 

                case("stress-b12") 
                    ! Use simple threshold method

                    call calc_calving_ground_rate_stress_b12(tpo%now%calv_grnd,tpo%now%H_ice,tpo%now%f_ice, &
                                                    tpo%now%f_grnd,bnd%z_bed,bnd%z_sl-bnd%z_bed,tpo%par%calv_tau)

                case DEFAULT 

                    write(*,*) "calc_ytopo:: Error: grounded calving method not recognized."
                    write(*,*) "calv_grnd_method = ", trim(tpo%par%calv_grnd_method)
                    stop 

            end select
            
            ! Additionally include parameterized grounded calving 
            ! to account for grid resolution 
            call calc_calving_ground_rate_stdev(calv_sd,tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd, &
                                            bnd%z_bed_sd,tpo%par%sd_min,tpo%par%sd_max,tpo%par%calv_max,tpo%par%calv_tau)
            tpo%now%calv_grnd = tpo%now%calv_grnd + calv_sd 

            ! Apply calving step 
            if (call_calving) then 
                call calc_ice_thickness_calving(tpo%now%H_ice,tpo%now%f_ice,tpo%now%calv, &
                                                tpo%now%f_grnd,bnd%z_sl-bnd%z_bed, &
                                                tpo%now%calv_flt,tpo%now%calv_grnd,tpo%par%dx,dt)            

            else
                ! Calving was diagnosed but not applied 

                tpo%now%calv = 0.0_wp 

            end if
            
            ! Update ice fraction mask 
            call calc_ice_fraction(tpo%now%f_ice,tpo%now%H_ice,bnd%z_bed,bnd%z_sl,tpo%par%margin_flt_subgrid)
            
            ! Finally apply all additional (generally artificial) ice thickness adjustments 
            ! and store changes in residual mass balance field. 
            call apply_ice_thickness_boundaries(tpo%now%mb_resid,tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd, &
                                                dyn%now%uxy_b,bnd%ice_allowed,tpo%par%boundaries,bnd%H_ice_ref, &
                                                tpo%par%H_min_flt,tpo%par%H_min_grnd,dt,reset_mb_resid)


            ! Save the rate of change of ice thickness in output variable [m/a]
            tpo%now%dHicedt = (tpo%now%H_ice - tpo%now%H_ice_n) / dt 

            ! If desired, finally relax solution to reference state
            if (tpo%par%topo_rel .ne. 0) then 

                select case(trim(tpo%par%topo_rel_field))

                    case("H_ref")
                        ! Relax towards reference ice thickness field H_ref

                        call relax_ice_thickness(tpo%now%H_ice,tpo%now%f_grnd,bnd%H_ice_ref, &
                                                    tpo%par%topo_rel,tpo%par%topo_rel_tau,dt)
                    case("H_ice_n")
                        ! Relax towards previous iteration ice thickness 
                        ! (ie slow down changes)
                        ! ajr: needs testing, not sure if this works well or helps anything.

                        call relax_ice_thickness(tpo%now%H_ice,tpo%now%f_grnd,tpo%now%H_ice_n, &
                                                    tpo%par%topo_rel,tpo%par%topo_rel_tau,dt)
                    
                    case DEFAULT 

                        write(*,*) "calc_ytopo:: Error: topo_rel_field not recognized."
                        write(*,*) "topo_rel_field = ", trim(tpo%par%topo_rel_field)
                        stop 

                end select

                ! Again apply ice thickness boundaries to ensure relaxed fields are consistent 
                ! with desired limitations. 
                call apply_ice_thickness_boundaries(tpo%now%mb_resid,tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd, &
                                                    dyn%now%uxy_b,bnd%ice_allowed,tpo%par%boundaries,bnd%H_ice_ref, &
                                                    tpo%par%H_min_flt,tpo%par%H_min_grnd,dt,reset=.FALSE.)

            end if

        end if 

        ! Final update of ice fraction mask (or define it now for fixed topography)
        call calc_ice_fraction(tpo%now%f_ice,tpo%now%H_ice,bnd%z_bed,bnd%z_sl,tpo%par%margin_flt_subgrid)
        
        ! Calculate grounding overburden ice thickness 
        call calc_H_grnd(tpo%now%H_grnd,tpo%now%H_ice,tpo%now%f_ice,bnd%z_bed,bnd%z_sl)

        ! 2. Calculate additional topographic properties ------------------

        ! Store previous surface elevation on predictor step for calculating
        ! rate of change of surface elevation.
        if (trim(tpo%par%pc_step) .eq. "predictor") then 
            tpo%now%z_srf_n = tpo%now%z_srf 
        end if 

        ! Calculate the surface elevation
        ! Note: two functions that should give the same results
        !call calc_z_srf(tpo%now%z_srf,tpo%now%H_ice,tpo%now%f_ice,tpo%now%H_grnd,bnd%z_bed,bnd%z_sl)
        call calc_z_srf_max(tpo%now%z_srf,tpo%now%H_ice,tpo%now%f_ice,bnd%z_bed,bnd%z_sl)
         
        select case(tpo%par%surf_gl_method)
            ! Choose method to treat grounding line points when calculating surface elevation

            case(0)
                ! Binary (grounded elevation or floating elevation via archemedes)
                
                ! Do nothing - surface elevation has already been properly diagnosed.
                
            case(1)
                ! Subgrid z_srf calculations at the grounding line 

                call calc_z_srf_gl_subgrid_area(tpo%now%z_srf,tpo%now%f_grnd,tpo%now%H_ice,tpo%now%f_ice,bnd%z_bed,bnd%z_sl,tpo%par%gl_sep_nx)

        end select 

        ! Determine rate of surface elevation change 
        if (dt .gt. 0.0) then 
            tpo%now%dzsrfdt = (tpo%now%z_srf-tpo%now%z_srf_n) / dt 
        end if 
        
        ! Calculate the ice thickness gradient (on staggered acx/y nodes)
        call calc_gradient_ac_ice(tpo%now%dHicedx,tpo%now%dHicedy,tpo%now%H_ice,tpo%now%f_ice,tpo%par%dx, &
                                                tpo%par%margin2nd,tpo%par%grad_lim,tpo%par%boundaries,zero_outside=.TRUE.)
        
        ! Calculate the surface slope
        call calc_gradient_ac(tpo%now%dzsdx,tpo%now%dzsdy,tpo%now%z_srf,tpo%par%dx)

        ! ajr: experimental, doesn't seem to work properly yet! ===>
        ! Modify surface slope gradient at the grounding line if desired 
!         call calc_gradient_ac_gl(tpo%now%dzsdx,tpo%now%dzsdy,tpo%now%z_srf,tpo%now%H_ice, &
!                                       tpo%now%f_grnd_acx,tpo%now%f_grnd_acy,tpo%par%dx,method=2,grad_lim=tpo%par%grad_lim)
        
        ! 3. Calculate new masks ------------------------------

        ! Calculate the grounded fraction and grounding line mask of each grid cell
        select case(tpo%par%gl_sep)

            case(1) 
                ! Binary f_grnd, linear f_grnd_acx/acy based on H_grnd

                call calc_f_grnd_subgrid_linear(tpo%now%f_grnd,tpo%now%f_grnd_acx,tpo%now%f_grnd_acy,tpo%now%H_grnd)

            case(2)
                ! Grounded area f_grnd, average to f_grnd_acx/acy 

                call calc_f_grnd_subgrid_area(tpo%now%f_grnd,tpo%now%f_grnd_acx,tpo%now%f_grnd_acy, &
                                                                            tpo%now%H_grnd,tpo%par%gl_sep_nx)
            
            case(3) 
                ! Grounded area using analytical solutions of Leguy et al. (2021)

                call determine_grounded_fractions(tpo%now%f_grnd,tpo%now%f_grnd_acx,tpo%now%f_grnd_acy, &
                                                                            tpo%now%f_grnd_ab,tpo%now%H_grnd)

        end select
        
        ! Calculate grounded fraction due to pinning points 
        call calc_f_grnd_pinning_points(tpo%now%f_grnd_pin,tpo%now%H_ice,tpo%now%f_ice, &
                                                            bnd%z_bed,bnd%z_bed_sd,bnd%z_sl)

!         ! Filter f_grnd to avoid lakes of one grid point inside of grounded ice 
!         ! ajr: note, this should be improved to treat ac-nodes too 
!         call filter_f_grnd(tpo%now%f_grnd)

        ! Calculate the grounding line mask 
        call calc_grline(tpo%now%is_grline,tpo%now%is_grz,tpo%now%f_grnd)

        ! Calculate the bed mask
        call gen_mask_bed(tpo%now%mask_bed,tpo%now%f_ice,thrm%now%f_pmp,tpo%now%f_grnd,tpo%now%is_grline)

        ! Calculate distance to ice margin (really slow if always on)
        !tpo%now%dist_margin = distance_to_margin(tpo%now%H_ice,tpo%par%dx)

        ! Calculate distance to grounding line (really slow if always on)
        !tpo%now%dist_grline = distance_to_grline(tpo%now%is_grline,tpo%now%f_grnd,tpo%par%dx)


        ! Determine ice thickness for use exclusively with the dynamics solver

        if (dyn%par%ssa_lat_bc .eq. "slab") then 
            ! Calculate extended ice thickness fields for use with dynamics solver
            ! Note: should not be used with MISMIP, TROUGH, etc. 

            tpo%now%H_ice_dyn = tpo%now%H_ice
            where (tpo%now%f_ice .lt. 1.0) tpo%now%H_ice_dyn = 1.0_wp
            
            ! Calculate the ice fraction mask for use with the dynamics solver
            call calc_ice_fraction(tpo%now%f_ice_dyn,tpo%now%H_ice_dyn,bnd%z_bed,bnd%z_sl,flt_subgrid=.FALSE.)
            
        else
            ! Set standard ice thickness field for use with dynamics 
            tpo%now%H_ice_dyn = tpo%now%H_ice
            tpo%now%f_ice_dyn = tpo%now%f_ice 
        end if 

            


        ! Store predicted/corrected ice thickness for later use 
        ! Do it here to ensure all changes to H_ice are accounted for (mb, calving, etc)
        if (trim(tpo%par%pc_step) .eq. "predictor") then 
            tpo%now%H_ice_pred = tpo%now%H_ice 
        else
            tpo%now%H_ice_corr = tpo%now%H_ice 
        end if 

        ! ================================

        ! Calculate computational performance (model speed in kyr/hr)
        call yelmo_cpu_time(cpu_time1)
        model_time1 = tpo%par%time 
        call yelmo_calc_speed(speed,model_time0,model_time1,cpu_time0,cpu_time1)

        ! Store the speed variable in predictor or corrector speed variable
        if (trim(tpo%par%pc_step) .eq. "predictor") then 
            tpo%par%speed_pred = speed
        else 
            tpo%par%speed_corr = speed 

            ! If corrector step, then also calculate the speed of both 
            ! predictor+corrector calls: mean of the predictor and corrector speeds
            ! divided by two, since two calls were made. 
            tpo%par%speed = 0.5_wp * (0.5_wp*(tpo%par%speed_pred+tpo%par%speed_corr))
            
        end if 

        if (trim(tpo%par%pc_step) .eq. "corrector") then 
            ! Advance ytopo timestep on corrector step 

            tpo%par%time = dble(time)
            
            if (call_calving) tpo%par%time_calv = dble(time)
        end if 

!         if (yelmo_log) then 

!             if (count(tpo%now%H_ice.gt.0.0) .gt. 0) then 
!                 write(*,"(a,f14.4,f10.4,f10.2)") "calc_ytopo::  time = ", tpo%par%time, dt, &
!                     sum(tpo%now%H_ice,mask=tpo%now%H_ice.gt.0.0)/real(count(tpo%now%H_ice.gt.0.0))
!             else 
!                 write(*,"(a,f14.4,f10.4,f10.2)") "calc_ytopo::  time = ", tpo%par%time, dt, 0.0
!             end if 

!         end if 

        return 

    end subroutine calc_ytopo

    elemental subroutine gen_mask_bed(mask,f_ice,f_pmp,f_grnd,is_grline)
        ! Generate an output mask for model conditions at bed
        ! based on input masks 
        ! 0: ocean, 1: land, 2: sia, 3: streams, grline: 4, floating: 5, islands: 6
        ! 7: partially-covered ice cell.

        implicit none 

        integer,    intent(OUT) :: mask 
        real(prec), intent(IN)  :: f_ice, f_pmp, f_grnd
        logical,    intent(IN)  :: is_grline

        if (is_grline) then
            ! Grounding line

            mask = mask_bed_grline

        else if (f_ice .eq. 0.0) then 
            ! Ice-free points 

            if (f_grnd .gt. 0.0) then
                ! Ice-free land

                mask = mask_bed_land

            else
                ! Ice-free ocean

                mask = mask_bed_ocean

            end if 

        else if (f_ice .gt. 0.0 .and. f_ice .lt. 1.0) then 
            ! Partially ice-covered points 

            mask = mask_bed_partial

        else
            ! Fully ice-covered points 

            if (f_grnd .gt. 0.0) then
                ! Grounded ice-covered points 

                if (f_pmp .gt. 0.5) then 
                    ! Temperate points

                    mask = mask_bed_stream 

                else
                    ! Frozen points 

                    mask = mask_bed_frozen 

                end if 

            else
                ! Floating ice-covered points 

                mask = mask_bed_float

            end if 

        end if 

        return 

    end subroutine gen_mask_bed

    subroutine ytopo_par_load(par,filename,nx,ny,dx,init)

        type(ytopo_param_class), intent(OUT) :: par
        character(len=*),        intent(IN)  :: filename
        integer,                  intent(IN)  :: nx, ny 
        real(prec),               intent(IN)  :: dx  
        logical, optional,       intent(IN)  :: init 

        ! Local variables
        logical :: init_pars 

        init_pars = .FALSE.
        if (present(init)) init_pars = .TRUE. 
 
        ! Store parameter values in output object
        call nml_read(filename,"ytopo","solver",            par%solver,           init=init_pars)
        call nml_read(filename,"ytopo","surf_gl_method",    par%surf_gl_method,   init=init_pars)
        call nml_read(filename,"ytopo","calv_flt_method",   par%calv_flt_method,  init=init_pars)
        call nml_read(filename,"ytopo","calv_grnd_method",  par%calv_grnd_method, init=init_pars)
        call nml_read(filename,"ytopo","bmb_gl_method",     par%bmb_gl_method,    init=init_pars)
        call nml_read(filename,"ytopo","fmb_method",        par%fmb_method,       init=init_pars)
        call nml_read(filename,"ytopo","margin2nd",         par%margin2nd,        init=init_pars)
        call nml_read(filename,"ytopo","margin_flt_subgrid",par%margin_flt_subgrid,init=init_pars)
        call nml_read(filename,"ytopo","use_bmb",           par%use_bmb,          init=init_pars)
        call nml_read(filename,"ytopo","topo_fixed",        par%topo_fixed,       init=init_pars)
        call nml_read(filename,"ytopo","topo_rel",          par%topo_rel,         init=init_pars)
        call nml_read(filename,"ytopo","topo_rel_tau",      par%topo_rel_tau,     init=init_pars)
        call nml_read(filename,"ytopo","topo_rel_field",    par%topo_rel_field,   init=init_pars)
        call nml_read(filename,"ytopo","calv_H_lim",        par%calv_H_lim,       init=init_pars)
        call nml_read(filename,"ytopo","calv_tau",          par%calv_tau,         init=init_pars)
        call nml_read(filename,"ytopo","calv_thin",         par%calv_thin,        init=init_pars)
        call nml_read(filename,"ytopo","H_min_grnd",        par%H_min_grnd,       init=init_pars)
        call nml_read(filename,"ytopo","H_min_flt",         par%H_min_flt,        init=init_pars)
        call nml_read(filename,"ytopo","sd_min",            par%sd_min,           init=init_pars)
        call nml_read(filename,"ytopo","sd_max",            par%sd_max,           init=init_pars)
        call nml_read(filename,"ytopo","calv_max",          par%calv_max,         init=init_pars)
        call nml_read(filename,"ytopo","grad_lim",          par%grad_lim,         init=init_pars)
        call nml_read(filename,"ytopo","gl_sep",            par%gl_sep,           init=init_pars)
        call nml_read(filename,"ytopo","gl_sep_nx",         par%gl_sep_nx,        init=init_pars)
        call nml_read(filename,"ytopo","diffuse_bmb_shlf",  par%diffuse_bmb_shlf, init=init_pars)
        call nml_read(filename,"ytopo","fmb_scale",         par%fmb_scale,        init=init_pars)
        call nml_read(filename,"ytopo","kt",                par%kt,               init=init_pars)
        call nml_read(filename,"ytopo","w2",                par%w2,               init=init_pars)
        call nml_read(filename,"ytopo","k2",                par%k2,               init=init_pars)
        
        ! === Set internal parameters =====

        par%nx  = nx 
        par%ny  = ny 
        par%dx  = dx 
        par%dy  = dx 

        ! Define how boundaries of grid should be treated 
        ! This should only be modified by the dom%par%experiment variable
        ! in yelmo_init. By default set boundaries to zero 
        par%boundaries = "zeros" 
        
        ! Define current time as unrealistic value
        par%time      = 1000000000   ! [a] 1 billion years in the future 
        par%time_calv = par%time 

        ! Intialize timestepping parameters to Forward Euler (beta2=beta4=0: no contribution from previous timestep)
        par%dt_zeta     = 1.0 
        par%dt_beta(1)  = 1.0 
        par%dt_beta(2)  = 0.0 
        par%dt_beta(3)  = 1.0 
        par%dt_beta(4)  = 0.0 

        ! Set some additional values to start out right
        par%pc_step    = "predictor"
        par%speed_pred = 0.0_wp 
        par%speed_corr = 0.0_wp 

        return

    end subroutine ytopo_par_load
    
    subroutine ytopo_alloc(now,nx,ny)

        implicit none 

        type(ytopo_state_class), intent(INOUT) :: now 
        integer, intent(IN) :: nx, ny  

        call ytopo_dealloc(now)

        allocate(now%H_ice(nx,ny))
        allocate(now%z_srf(nx,ny))
        allocate(now%dzsrfdt(nx,ny))
        allocate(now%dHicedt(nx,ny))
        allocate(now%bmb(nx,ny))
        allocate(now%fmb(nx,ny))
        allocate(now%mb_applied(nx,ny))
        allocate(now%mb_resid(nx,ny))
        
        allocate(now%eps_eff(nx,ny))
        allocate(now%tau_eff(nx,ny))
        allocate(now%calv(nx,ny))
        allocate(now%calv_flt(nx,ny))
        allocate(now%calv_grnd(nx,ny))
        
        allocate(now%dzsdx(nx,ny))
        allocate(now%dzsdy(nx,ny))

        allocate(now%dHicedx(nx,ny))
        allocate(now%dHicedy(nx,ny))
        
        allocate(now%H_eff(nx,ny))
        allocate(now%H_grnd(nx,ny))

        ! Masks 
        allocate(now%f_grnd(nx,ny))
        allocate(now%f_grnd_acx(nx,ny))
        allocate(now%f_grnd_acy(nx,ny))
        allocate(now%f_grnd_ab(nx,ny))
        allocate(now%f_ice(nx,ny))

        allocate(now%f_grnd_pin(nx,ny))

        allocate(now%dist_margin(nx,ny))
        allocate(now%dist_grline(nx,ny))
        
        allocate(now%mask_bed(nx,ny))
        allocate(now%is_grline(nx,ny))
        allocate(now%is_grz(nx,ny))

        allocate(now%dHdt_n(nx,ny))
        allocate(now%H_ice_n(nx,ny))
        allocate(now%H_ice_pred(nx,ny))
        allocate(now%H_ice_corr(nx,ny))
        
        allocate(now%z_srf_n(nx,ny))
        
        allocate(now%H_ice_dyn(nx,ny))
        allocate(now%f_ice_dyn(nx,ny))

        now%H_ice       = 0.0 
        now%z_srf       = 0.0  
        now%dzsrfdt     = 0.0 
        now%dHicedt     = 0.0
        now%bmb         = 0.0  
        now%fmb         = 0.0
        now%mb_applied  = 0.0 
        now%mb_resid    = 0.0
        now%eps_eff     = 0.0
        now%tau_eff     = 0.0
        now%calv        = 0.0
        now%calv_flt    = 0.0
        now%calv_grnd   = 0.0
        now%dzsdx       = 0.0 
        now%dzsdy       = 0.0 
        now%dHicedx     = 0.0 
        now%dHicedy     = 0.0
        now%H_eff       = 0.0 
        now%H_grnd      = 0.0  
        now%f_grnd      = 0.0  
        now%f_grnd_acx  = 0.0  
        now%f_grnd_acy  = 0.0  
        now%f_grnd_ab   = 0.0
        now%f_grnd_pin  = 0.0
        now%f_ice       = 0.0  
        now%dist_margin = 0.0
        now%dist_grline = 0.0 
        
        now%mask_bed    = 0.0 
        now%is_grline   = .FALSE. 
        now%is_grz      = .FALSE. 
         
        now%dHdt_n      = 0.0  
        now%H_ice_n     = 0.0 
        now%H_ice_pred  = 0.0 
        now%H_ice_corr  = 0.0 
        
        now%z_srf_n     = 0.0 

        now%H_ice_dyn   = 0.0 
        now%f_ice_dyn   = 0.0 
        
        return 

    end subroutine ytopo_alloc

    subroutine ytopo_dealloc(now)

        implicit none 

        type(ytopo_state_class), intent(INOUT) :: now

        if (allocated(now%H_ice))       deallocate(now%H_ice)
        if (allocated(now%z_srf))       deallocate(now%z_srf)
        
        if (allocated(now%dzsrfdt))     deallocate(now%dzsrfdt)
        if (allocated(now%dHicedt))     deallocate(now%dHicedt)
        if (allocated(now%bmb))         deallocate(now%bmb)
        if (allocated(now%fmb))         deallocate(now%fmb)
        if (allocated(now%mb_applied))  deallocate(now%mb_applied)
        if (allocated(now%mb_resid))    deallocate(now%mb_resid)
        
        if (allocated(now%eps_eff))     deallocate(now%eps_eff)
        if (allocated(now%tau_eff))     deallocate(now%tau_eff)
        if (allocated(now%calv))        deallocate(now%calv)
        if (allocated(now%calv_flt))    deallocate(now%calv_flt)
        if (allocated(now%calv_grnd))   deallocate(now%calv_grnd)
            
        if (allocated(now%dzsdx))       deallocate(now%dzsdx)
        if (allocated(now%dzsdy))       deallocate(now%dzsdy)
        if (allocated(now%dHicedx))     deallocate(now%dHicedx)
        if (allocated(now%dHicedy))     deallocate(now%dHicedy)
        
        if (allocated(now%H_eff))       deallocate(now%H_eff)
        if (allocated(now%H_grnd))      deallocate(now%H_grnd)

        if (allocated(now%f_grnd))      deallocate(now%f_grnd)
        if (allocated(now%f_grnd_acx))  deallocate(now%f_grnd_acx)
        if (allocated(now%f_grnd_acy))  deallocate(now%f_grnd_acy)
        if (allocated(now%f_grnd_ab))   deallocate(now%f_grnd_ab)
        if (allocated(now%f_grnd_pin))  deallocate(now%f_grnd_pin)

        if (allocated(now%f_ice))       deallocate(now%f_ice)

        if (allocated(now%dist_margin)) deallocate(now%dist_margin)
        if (allocated(now%dist_grline)) deallocate(now%dist_grline)
        
        if (allocated(now%mask_bed))    deallocate(now%mask_bed)
        if (allocated(now%is_grline))   deallocate(now%is_grline)
        if (allocated(now%is_grz))      deallocate(now%is_grz)

        if (allocated(now%dHdt_n))      deallocate(now%dHdt_n)
        if (allocated(now%H_ice_n))     deallocate(now%H_ice_n)
        if (allocated(now%H_ice_pred))  deallocate(now%H_ice_pred)
        if (allocated(now%H_ice_corr))  deallocate(now%H_ice_corr)
        
        if (allocated(now%z_srf_n))     deallocate(now%z_srf_n)
        
        if (allocated(now%H_ice_dyn))   deallocate(now%H_ice_dyn)
        if (allocated(now%f_ice_dyn))   deallocate(now%f_ice_dyn)
        
        return 

    end subroutine ytopo_dealloc
    
end module yelmo_topography
