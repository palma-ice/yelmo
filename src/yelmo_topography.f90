
module yelmo_topography

    use nml  
    use ncio
    
    use yelmo_defs
    use yelmo_tools 
    
    use mass_conservation
    use calving
    use topography 
    use discharge

    use runge_kutta 
    
    implicit none
    
    private
    public :: calc_ytopo_pc 
    public :: calc_ytopo_diagnostic 
    public :: calc_ytopo_rates
    public :: ytopo_par_load
    public :: ytopo_alloc
    public :: ytopo_dealloc
    
    ! Integers
    public :: mask_bed_ocean  
    public :: mask_bed_land  
    public :: mask_bed_frozen
    public :: mask_bed_stream
    public :: mask_bed_grline
    public :: mask_bed_float 
    public :: mask_bed_island
    
contains
    
    subroutine calc_ytopo_pc(tpo,dyn,mat,thrm,bnd,time,topo_fixed,pc_step,use_H_pred)

        implicit none 

        type(ytopo_class),  intent(INOUT) :: tpo
        type(ydyn_class),   intent(IN)    :: dyn
        type(ymat_class),   intent(IN)    :: mat
        type(ytherm_class), intent(IN)    :: thrm  
        type(ybound_class), intent(IN)    :: bnd 
        real(wp),           intent(IN)    :: time
        logical,            intent(IN)    :: topo_fixed  
        character(len=*),   intent(IN)    :: pc_step 
        logical, optional,  intent(IN)    :: use_H_pred

        ! Local variables 
        integer  :: i, j, nx, ny
        real(wp) :: dt  
        real(wp), allocatable :: dHidt_now(:,:) 
        real(wp), allocatable :: H_prev(:,:)

        logical, parameter :: use_rk4 = .FALSE. 

        nx = size(tpo%now%H_ice,1)
        ny = size(tpo%now%H_ice,2)

        allocate(dHidt_now(nx,ny))
        allocate(H_prev(nx,ny))

        ! Initialize time if necessary 
        if (tpo%par%time .gt. dble(time)) then 
            tpo%par%time = dble(time) 
        end if 
        
        ! Initialize rk4 object if necessary 
        if (.not. allocated(tpo%rk4%tau)) then 
            call rk4_2D_init(tpo%rk4,nx,ny,dt_min=1e-3_wp)
        end if 

        ! Get time step
        dt = dble(time) - tpo%par%time 

        ! Step 0: Get some diagnostic quantities for mass balance calculation --------

        ! Get ice thickness entering routine
        H_prev = tpo%now%H_ice 

        ! Step 1: Go through predictor-corrector-advance steps

        if ( .not. topo_fixed .and. dt .gt. 0.0 ) then 

            ! Ice thickness evolution from dynamics alone

            select case(trim(pc_step))

                case("predictor") 
                    ! Determine predicted ice thickness 

                    ! Store dynamic rate of change from previous timestep,
                    ! along with ice thickness and surface elevation
                    ! (the latter only for calculating rate of change later)
                    tpo%now%dHidt_dyn_n = tpo%now%dHidt_dyn
                    tpo%now%H_ice_n     = tpo%now%H_ice 
                    tpo%now%z_srf_n     = tpo%now%z_srf 

                    ! Get ice-fraction mask for current ice thickness  
                    call calc_ice_fraction(tpo%now%f_ice,tpo%now%H_ice,bnd%z_bed,bnd%z_sl,bnd%c%rho_ice, &
                                            bnd%c%rho_sw,tpo%par%boundaries,tpo%par%margin_flt_subgrid)

if (use_rk4) then
                    call rk4_2D_step(tpo%rk4,tpo%now%H_ice,tpo%now%f_ice,dHidt_now,dyn%now%ux_bar,dyn%now%uy_bar, &
                                                tpo%now%mask_adv,tpo%par%dx,dt,tpo%par%solver,tpo%par%boundaries)

else
                    call calc_G_advec_simple(dHidt_now,tpo%now%H_ice,tpo%now%f_ice,dyn%now%ux_bar,dyn%now%uy_bar, &
                                                 tpo%now%mask_adv,tpo%par%solver,tpo%par%boundaries,tpo%par%dx,dt)
                 
end if 
                    
                    ! Calculate rate of change using weighted advective rates of change 
                    ! depending on timestepping method chosen 
                    tpo%now%dHidt_dyn = tpo%par%dt_beta(1)*dHidt_now + tpo%par%dt_beta(2)*tpo%now%dHidt_dyn_n 

                    ! Apply rate and update ice thickness (predicted)
                    tpo%now%H_ice = tpo%now%H_ice_n
                    call apply_tendency(tpo%now%H_ice,tpo%now%dHidt_dyn,dt,"dyn_pred",adjust_mb=.FALSE.)

                case("corrector") 

                    ! Set current thickness to predicted thickness
                    tpo%now%H_ice = tpo%now%pred%H_ice 

                    ! Get ice-fraction mask for predicted ice thickness  
                    call calc_ice_fraction(tpo%now%f_ice,tpo%now%H_ice,bnd%z_bed,bnd%z_sl,bnd%c%rho_ice, &
                                            bnd%c%rho_sw,tpo%par%boundaries,tpo%par%margin_flt_subgrid)

if (use_rk4) then
                    call rk4_2D_step(tpo%rk4,tpo%now%H_ice,tpo%now%f_ice,dHidt_now,dyn%now%ux_bar,dyn%now%uy_bar, &
                                                tpo%now%mask_adv,tpo%par%dx,dt,tpo%par%solver,tpo%par%boundaries)
else
                    call calc_G_advec_simple(dHidt_now,tpo%now%H_ice,tpo%now%f_ice,dyn%now%ux_bar,dyn%now%uy_bar, &
                                                tpo%now%mask_adv,tpo%par%solver,tpo%par%boundaries,tpo%par%dx,dt)
                 
end if

                    ! Calculate rate of change using weighted advective rates of change 
                    ! depending on timestepping method chosen 
                    tpo%now%dHidt_dyn = tpo%par%dt_beta(3)*dHidt_now + tpo%par%dt_beta(4)*tpo%now%dHidt_dyn_n 
                    
                    ! Apply rate and update ice thickness (corrected)
                    tpo%now%H_ice = tpo%now%H_ice_n
                    call apply_tendency(tpo%now%H_ice,tpo%now%dHidt_dyn,dt,"dyn_corr",adjust_mb=.FALSE.)
                    
            end select

            ! Note: at this point, mass has only been advected (moved around). In principle,
            ! this is fully conservative and the net Δmb=0. However, due to the predictor-corrector
            ! mixing of dHdt with previous timesteps/iterations, some small amounts of negative
            ! ice thickness can arise. The quantities are much smaller than other mb quantities
            ! and localized at the margin, so it should not be problematic. They should be captured
            ! and corrected in the apply_tendency routine below, where adjust_mb=.TRUE. ensures
            ! the ice thickness stays >= 0. Alternatively, adjust_mb=.TRUE. can be imposed above,
            ! but this implies that the advective dHdt fields are less precise, potentially
            ! impacting the pc-stability.

            select case(trim(pc_step))

                case("predictor","corrector")
                    ! For either predictor or corrector step, also calculate all mass balance changes

                    ! === smb =====

                    call calc_G_mbal(tpo%now%smb,tpo%now%H_ice,tpo%now%f_grnd,bnd%smb,dt)

                    ! Apply rate and update ice thickness
                    call apply_tendency(tpo%now%H_ice,tpo%now%smb,dt,"smb",adjust_mb=.TRUE.)
                    
                    ! === bmb =====

                    ! Calculate grounded fraction on aa-nodes
                    ! (only to be used with basal mass balance, later all
                    !  f_grnd arrays will be calculated according to use choices)
                    call determine_grounded_fractions(tpo%now%f_grnd_bmb,H_grnd=tpo%now%H_grnd, &
                                                                        boundaries=tpo%par%boundaries)
                    
                    ! Combine basal mass balance into one field accounting for 
                    ! grounded/floating fraction of grid cells 
                    call calc_bmb_total(tpo%now%bmb_ref,thrm%now%bmb_grnd,bnd%bmb_shlf,tpo%now%H_ice, &
                                        tpo%now%H_grnd,tpo%now%f_grnd_bmb,tpo%par%gz_Hg0,tpo%par%gz_Hg1, &
                                        tpo%par%gz_nx,tpo%par%bmb_gl_method,tpo%par%boundaries)

                    if (tpo%par%use_bmb) then
                        call calc_G_mbal(tpo%now%bmb,tpo%now%H_ice,tpo%now%f_grnd,tpo%now%bmb_ref,dt)
                    else
                        ! Mainly for when running EISMINT1
                        tpo%now%bmb = 0.0
                    end if

                    ! Apply rate and update ice thickness
                    call apply_tendency(tpo%now%H_ice,tpo%now%bmb,dt,"bmb",adjust_mb=.TRUE.)

                    ! === fmb =====

                    ! Calculate frontal mass balance
                    call calc_fmb_total(tpo%now%fmb_ref,bnd%fmb_shlf,bnd%bmb_shlf,tpo%now%H_ice, &
                                    tpo%now%H_grnd,tpo%now%f_ice,tpo%par%fmb_method,tpo%par%fmb_scale, &
                                    bnd%c%rho_ice,bnd%c%rho_sw,tpo%par%dx,tpo%par%boundaries)

                    if (tpo%par%use_bmb) then
                        call calc_G_mbal(tpo%now%fmb,tpo%now%H_ice,tpo%now%f_grnd,tpo%now%fmb_ref,dt)
                    else
                        ! Mainly for when running EISMINT1
                        tpo%now%fmb = 0.0
                    end if

                    ! Apply rate and update ice thickness
                    call apply_tendency(tpo%now%H_ice,tpo%now%fmb,dt,"fmb",adjust_mb=.TRUE.)
                    
                    ! === dmb =====

                    call calc_mb_discharge(tpo%now%dmb_ref,tpo%now%H_ice,tpo%now%z_srf,bnd%z_bed_sd,tpo%now%dist_grline, &
                                tpo%now%dist_margin,tpo%now%f_ice,tpo%par%dmb_method,tpo%par%dx,tpo%par%dmb_alpha_max, &
                                tpo%par%dmb_tau,tpo%par%dmb_sigma_ref,tpo%par%dmb_m_d,tpo%par%dmb_m_r)
                    
                    call calc_G_mbal(tpo%now%dmb,tpo%now%H_ice,tpo%now%f_grnd,tpo%now%dmb_ref,dt)

                    ! Apply rate and update ice thickness
                    call apply_tendency(tpo%now%H_ice,tpo%now%dmb,dt,"dmb",adjust_mb=.TRUE.)
                    
                    ! === mb_net =====

                    tpo%now%mb_net = tpo%now%smb + tpo%now%bmb + tpo%now%fmb + tpo%now%dmb 

                    ! Calculate and apply calving
                    call calc_ytopo_calving(tpo,dyn,mat,thrm,bnd,dt)

                    ! Get ice-fraction mask for ice thickness  
                    call calc_ice_fraction(tpo%now%f_ice,tpo%now%H_ice,bnd%z_bed,bnd%z_sl,bnd%c%rho_ice, &
                                            bnd%c%rho_sw,tpo%par%boundaries,tpo%par%margin_flt_subgrid)

                    ! If desired, finally relax solution to reference state
                    if (tpo%par%topo_rel .ne. 0) then 

                        if (tpo%par%topo_rel .eq. -1) then
                            ! Use externally defined tau_relax field
                            tpo%now%tau_relax = bnd%tau_relax
                        else
                            ! Define the relaxation timescale field, if needed
                            call set_tau_relax(tpo%now%tau_relax,tpo%now%H_ice,tpo%now%f_grnd,tpo%now%mask_grz,bnd%H_ice_ref, &
                                               tpo%par%topo_rel,tpo%par%topo_rel_tau,tpo%par%boundaries)
                        end if

                        select case(trim(tpo%par%topo_rel_field))

                            case("H_ref")
                                ! Relax towards reference ice thickness field H_ref

                                call calc_G_relaxation(tpo%now%mb_relax,tpo%now%H_ice,bnd%H_ice_ref,tpo%now%tau_relax,dt)

                            case("H_ice_n")
                                ! Relax towards previous iteration ice thickness 
                                ! (ie slow down changes)
                                ! ajr: needs testing, not sure if this works well or helps anything.

                                call calc_G_relaxation(tpo%now%mb_relax,tpo%now%H_ice,tpo%now%H_ice_n,tpo%now%tau_relax,dt)

                            case DEFAULT 

                                write(*,*) "calc_ytopo:: Error: topo_rel_field not recognized."
                                write(*,*) "topo_rel_field = ", trim(tpo%par%topo_rel_field)
                                stop 

                        end select

                        ! Apply rate and update ice thickness
                        call apply_tendency(tpo%now%H_ice,tpo%now%mb_relax,dt,"relax",adjust_mb=.TRUE.)

                        ! Add relaxation tendency to mb_net for proper accounting of mass change
                        tpo%now%mb_net = tpo%now%mb_net + tpo%now%mb_relax

                        ! Get ice-fraction mask for ice thickness  
                        call calc_ice_fraction(tpo%now%f_ice,tpo%now%H_ice,bnd%z_bed,bnd%z_sl,bnd%c%rho_ice, &
                                                bnd%c%rho_sw,tpo%par%boundaries,tpo%par%margin_flt_subgrid)

                    end if

                    ! Finally, apply all additional (generally artificial) ice thickness adjustments 
                    ! and store changes in residual mass balance field. 
                    call calc_G_boundaries(tpo%now%mb_resid,tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd, &
                                            dyn%now%uxy_b,bnd%ice_allowed,tpo%par%boundaries,bnd%H_ice_ref, &
                                            tpo%par%H_min_flt,tpo%par%H_min_grnd,dt)

                    ! Apply rate and update ice thickness
                    call apply_tendency(tpo%now%H_ice,tpo%now%mb_resid,dt,"resid",adjust_mb=.TRUE.)

                    ! Add residual tendency to mb_net for proper accounting of mass change
                    tpo%now%mb_net = tpo%now%mb_net + tpo%now%mb_resid

                    ! Get ice-fraction mask for ice thickness  
                    call calc_ice_fraction(tpo%now%f_ice,tpo%now%H_ice,bnd%z_bed,bnd%z_sl,bnd%c%rho_ice, &
                                            bnd%c%rho_sw,tpo%par%boundaries,tpo%par%margin_flt_subgrid)

            end select 

            select case(trim(pc_step))

                case("predictor") 

                    ! Save current predictor fields, 
                    ! proceed with predictor fields for calculating dynamics.
                    tpo%now%pred%H_ice      = tpo%now%H_ice 
                    tpo%now%pred%dHidt_dyn  = tpo%now%dHidt_dyn
                    tpo%now%pred%mb_net     = tpo%now%mb_net 
                    tpo%now%pred%mb_relax   = tpo%now%mb_relax 
                    tpo%now%pred%mb_resid   = tpo%now%mb_resid 
                    tpo%now%pred%smb        = tpo%now%smb
                    tpo%now%pred%bmb        = tpo%now%bmb
                    tpo%now%pred%fmb        = tpo%now%fmb
                    tpo%now%pred%dmb        = tpo%now%dmb
                    tpo%now%pred%cmb        = tpo%now%cmb 
                    tpo%now%pred%cmb_flt    = tpo%now%cmb_flt 
                    tpo%now%pred%cmb_grnd   = tpo%now%cmb_grnd 
                    
                case("corrector")
                    ! Determine corrected ice thickness 

                    ! Save current corrector fields
                    tpo%now%corr%H_ice      = tpo%now%H_ice 
                    tpo%now%corr%dHidt_dyn  = tpo%now%dHidt_dyn
                    tpo%now%corr%mb_net     = tpo%now%mb_net 
                    tpo%now%corr%mb_relax   = tpo%now%mb_relax 
                    tpo%now%corr%mb_resid   = tpo%now%mb_resid 
                    tpo%now%corr%smb        = tpo%now%smb
                    tpo%now%corr%bmb        = tpo%now%bmb
                    tpo%now%corr%fmb        = tpo%now%fmb
                    tpo%now%corr%dmb        = tpo%now%dmb
                    tpo%now%corr%cmb        = tpo%now%cmb 
                    tpo%now%corr%cmb_flt    = tpo%now%cmb_flt 
                    tpo%now%corr%cmb_grnd   = tpo%now%cmb_grnd
                    
                    ! Restore main ice thickness field to original 
                    ! value at the beginning of the timestep for 
                    ! calculation of remaining quantities (thermo, material)
                    tpo%now%H_ice = tpo%now%H_ice_n 

                case("advance")
                    ! Now let's actually advance the ice thickness field

                    if (.not. present(use_H_pred)) then 
                        write(*,*) "calc_ytopo_pc:: Error: &
                        & For step='advance', the argument use_H_pred&
                        & must be provided."
                        stop 
                    end if 

                    ! Determine which ice thickness to use going forward
                    if (use_H_pred) then 

                        ! Load predictor fields in current state variables
                        tpo%now%H_ice       = tpo%now%pred%H_ice 
                        tpo%now%dHidt_dyn   = tpo%now%pred%dHidt_dyn
                        tpo%now%mb_net      = tpo%now%pred%mb_net 
                        tpo%now%mb_relax    = tpo%now%pred%mb_relax 
                        tpo%now%mb_resid    = tpo%now%pred%mb_resid 
                        tpo%now%smb         = tpo%now%pred%smb 
                        tpo%now%bmb         = tpo%now%pred%bmb 
                        tpo%now%fmb         = tpo%now%pred%fmb 
                        tpo%now%dmb         = tpo%now%pred%dmb 
                        tpo%now%cmb         = tpo%now%pred%cmb 
                        tpo%now%cmb_flt     = tpo%now%pred%cmb_flt 
                        tpo%now%cmb_grnd    = tpo%now%pred%cmb_grnd 
                        
                    else
                        ! Load corrector fields in current state variables
                        tpo%now%H_ice       = tpo%now%corr%H_ice 
                        tpo%now%dHidt_dyn   = tpo%now%corr%dHidt_dyn
                        tpo%now%mb_net      = tpo%now%corr%mb_net 
                        tpo%now%mb_relax    = tpo%now%corr%mb_relax 
                        tpo%now%mb_resid    = tpo%now%corr%mb_resid 
                        tpo%now%smb         = tpo%now%corr%smb 
                        tpo%now%bmb         = tpo%now%corr%bmb 
                        tpo%now%fmb         = tpo%now%corr%fmb 
                        tpo%now%dmb         = tpo%now%corr%dmb 
                        tpo%now%cmb         = tpo%now%corr%cmb 
                        tpo%now%cmb_flt     = tpo%now%corr%cmb_flt 
                        tpo%now%cmb_grnd    = tpo%now%corr%cmb_grnd 
                        
                    end if
                    
            end select

            ! Determine rates of change
            tpo%now%dHidt = (tpo%now%H_ice - tpo%now%H_ice_n) / dt 
            tpo%now%dzsdt = (tpo%now%z_srf - tpo%now%z_srf_n) / dt 

            ! Determine mass balance error by comparing mass_in - mass_out to dHidt
            tpo%now%mb_err = tpo%now%dHidt - (tpo%now%mb_net + tpo%now%cmb)

        end if 

        ! Update fields and masks
        call calc_ytopo_diagnostic(tpo,dyn,mat,thrm,bnd)

        
        if (trim(pc_step) .eq. "advance") then 
            ! Advance timestep here whether topo_fixed was true or not...
            
            ! Update ytopo time to current time 
            tpo%par%time = dble(time)
            
        end if

        return

    end subroutine calc_ytopo_pc

    subroutine calc_ytopo_calving(tpo,dyn,mat,thrm,bnd,dt)

        implicit none 

        type(ytopo_class),  intent(INOUT) :: tpo
        type(ydyn_class),   intent(IN)    :: dyn
        type(ymat_class),   intent(IN)    :: mat
        type(ytherm_class), intent(IN)    :: thrm  
        type(ybound_class), intent(IN)    :: bnd 
        real(wp),           intent(IN)    :: dt

        ! Local variables 
        integer :: i, j, nx, ny 
        real(wp), allocatable :: mbal_now(:,:) 
        real(wp), allocatable :: cmb_sd(:,:) 

        nx = size(tpo%now%H_ice,1) 
        ny = size(tpo%now%H_ice,2) 

        allocate(mbal_now(nx,ny)) 
        allocate(cmb_sd(nx,ny)) 


        ! Make sure current ice mask is correct
        call calc_ice_fraction(tpo%now%f_ice,tpo%now%H_ice,bnd%z_bed,bnd%z_sl,bnd%c%rho_ice, &
                                    bnd%c%rho_sw,tpo%par%boundaries,tpo%par%margin_flt_subgrid)

        ! === CALVING ===

        ! == Diagnose strains and stresses relevant to calving ==

        ! eps_eff = effective strain = eigencalving e+*e- following Levermann et al. (2012)
        call calc_eps_eff(tpo%now%eps_eff,dyn%now%strn2D%eps_eig_1,dyn%now%strn2D%eps_eig_2,tpo%now%f_ice,tpo%par%boundaries)

        ! tau_eff = effective stress ~ von Mises stress following Lipscomb et al. (2019)
        call calc_tau_eff(tpo%now%tau_eff,mat%now%strs2D%tau_eig_1,mat%now%strs2D%tau_eig_2,tpo%now%f_ice,tpo%par%w2,tpo%par%boundaries)

        ! == Determine thickness threshold for calving spatially ==

        call define_calving_thickness_threshold(tpo%now%H_calv,tpo%now%z_bed_filt,tpo%par%Hc_ref, &
                                            tpo%par%Hc_deep,tpo%par%zb_deep_0,tpo%par%zb_deep_1)

        ! Define factor for calving stress spatially
        call define_calving_stress_factor(tpo%now%kt,tpo%now%z_bed_filt,tpo%par%kt_ref, &
                                            tpo%par%kt_deep,tpo%par%zb_deep_0,tpo%par%zb_deep_1)

        ! == Calculate potential floating calving rate ==

        select case(trim(tpo%par%calv_flt_method))

            case("zero","none")

                tpo%now%cmb_flt = 0.0 

            case("threshold") 
                ! Use threshold method

                call calc_calving_rate_threshold(tpo%now%cmb_flt,tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd, &
                                                 tpo%now%H_calv,tpo%par%calv_tau,tpo%par%boundaries)
                
            case("vm-l19")
                ! Use von Mises calving as defined by Lipscomb et al. (2019)

                ! Next, diagnose calving
                call calc_calving_rate_vonmises_l19(tpo%now%cmb_flt,tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd, &
                                                                        tpo%now%tau_eff,tpo%par%dx,tpo%now%kt,tpo%par%boundaries)

                ! Scale calving with 'thin' calving rate to ensure 
                ! small ice thicknesses are removed.
                call apply_calving_rate_thin(tpo%now%cmb_flt,tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd,tpo%par%calv_thin,tpo%par%Hc_ref_thin,tpo%par%boundaries)

            case("eigen")
                ! Use Eigen calving as defined by Levermann et al. (2012)

                ! Next, diagnose calving
                call calc_calving_rate_eigen(tpo%now%cmb_flt,tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd, &
                                                                        tpo%now%eps_eff,tpo%par%dx,tpo%par%k2,tpo%par%boundaries)

                ! Scale calving with 'thin' calving rate to ensure 
                ! small ice thicknesses are removed.
                call apply_calving_rate_thin(tpo%now%cmb_flt,tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd,tpo%par%calv_thin,tpo%par%Hc_ref_thin,tpo%par%boundaries)

            case("kill") 
                ! Delete all floating ice (using characteristic time parameter)
                ! Make sure dt is a postive number

                call calc_calving_rate_kill(tpo%now%cmb_flt,tpo%now%H_ice,tpo%now%f_grnd.eq.0.0_wp, &
                                                                                    tpo%par%calv_tau,dt)

            case("kill-pos")
                ! Delete all floating ice beyond a given location (using characteristic time parameter)

                call calc_calving_rate_kill(tpo%now%cmb_flt,tpo%now%H_ice, &
                                                ( tpo%now%f_grnd .eq. 0.0_wp .and. &
                                                  tpo%now%H_ice  .gt. 0.0_wp .and. &
                                                  bnd%calv_mask ), tau=0.0_wp, dt=dt )

            case DEFAULT 

                write(*,*) "calc_ytopo:: Error: floating calving method not recognized."
                write(*,*) "calv_flt_method = ", trim(tpo%par%calv_flt_method)
                stop 

        end select
        
        ! Additionally ensure higher calving rate for floating tongues of
        ! one grid-point width, or lower calving rate for embayed points.

        select case(trim(tpo%par%calv_flt_method))

            case("zero","none","kill","kill-pos")
                
                ! Do nothing for these methods

            case DEFAULT 

                call calc_calving_rate_tongues(tpo%now%cmb_flt,tpo%now%H_ice,tpo%now%f_ice, &
                                                tpo%now%f_grnd,tpo%par%calv_tau,tpo%par%boundaries)
        
        end select 

        
        ! Apply rate and update ice thickness
        call apply_tendency(tpo%now%H_ice,tpo%now%cmb_flt,dt,"calv_flt",adjust_mb=.TRUE.)


        ! Diagnose potential grounded-ice calving rate [m/yr]

        select case(trim(tpo%par%calv_grnd_method))

            case("zero","none")

                tpo%now%cmb_grnd = 0.0 

            case("stress-b12") 
                ! Use simple threshold method

                call calc_calving_ground_rate_stress_b12(tpo%now%cmb_grnd,tpo%now%H_ice,tpo%now%f_ice, &
                                                tpo%now%f_grnd,bnd%z_bed,bnd%z_sl-bnd%z_bed,tpo%par%calv_tau, &
                                                bnd%c%rho_ice,bnd%c%rho_sw,bnd%c%g,tpo%par%boundaries)

            case DEFAULT 

                write(*,*) "calc_ytopo:: Error: grounded calving method not recognized."
                write(*,*) "calv_grnd_method = ", trim(tpo%par%calv_grnd_method)
                stop 

        end select
        
        ! Additionally include parameterized grounded calving 
        ! to account for grid resolution 
        select case(trim(tpo%par%calv_grnd_method))

            case("zero","none")
                
                ! Do nothing for these methods

            case DEFAULT 

                call calc_calving_ground_rate_stdev(cmb_sd,tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd, &
                                bnd%z_bed_sd,tpo%par%sd_min,tpo%par%sd_max,tpo%par%calv_grnd_max,tpo%par%calv_tau,tpo%par%boundaries)
                tpo%now%cmb_grnd = tpo%now%cmb_grnd + cmb_sd 

        end select
        
        ! Apply rate and update ice thickness
        call apply_tendency(tpo%now%H_ice,tpo%now%cmb_grnd,dt,"calv_grnd",adjust_mb=.TRUE.)



        ! Finally, get the total combined calving mass balance
        tpo%now%cmb = tpo%now%cmb_flt + tpo%now%cmb_grnd 


        ! Update ice fraction mask 
        call calc_ice_fraction(tpo%now%f_ice,tpo%now%H_ice,bnd%z_bed,bnd%z_sl,bnd%c%rho_ice, &
                                bnd%c%rho_sw,tpo%par%boundaries,tpo%par%margin_flt_subgrid)


        ! Treat fractional points that are not connected to full ice-covered points
        call calc_G_remove_fractional_ice(mbal_now,tpo%now%H_ice,tpo%now%f_ice,dt)

        ! Apply rate and update ice thickness
        call apply_tendency(tpo%now%H_ice,mbal_now,dt,"frac",adjust_mb=.TRUE.)

        ! Add this rate to calving tendency
        tpo%now%cmb = tpo%now%cmb + mbal_now

        call calc_ice_fraction(tpo%now%f_ice,tpo%now%H_ice,bnd%z_bed,bnd%z_sl,bnd%c%rho_ice, &
                                bnd%c%rho_sw,tpo%par%boundaries,tpo%par%margin_flt_subgrid)

        return

    end subroutine calc_ytopo_calving

    subroutine calc_ytopo_diagnostic(tpo,dyn,mat,thrm,bnd)
        ! Calculate adjustments to surface elevation, bedrock elevation
        ! and ice thickness 

        implicit none 

        type(ytopo_class),  intent(INOUT) :: tpo
        type(ydyn_class),   intent(IN)    :: dyn
        type(ymat_class),   intent(IN)    :: mat
        type(ytherm_class), intent(IN)    :: thrm  
        type(ybound_class), intent(IN)    :: bnd 

        ! Final update of ice fraction mask (or define it now for fixed topography)
        call calc_ice_fraction(tpo%now%f_ice,tpo%now%H_ice,bnd%z_bed,bnd%z_sl,bnd%c%rho_ice, &
                                bnd%c%rho_sw,tpo%par%boundaries,tpo%par%margin_flt_subgrid)

        ! Calculate grounding overburden ice thickness 
        call calc_H_grnd(tpo%now%H_grnd,tpo%now%H_ice,tpo%now%f_ice,bnd%z_bed,bnd%z_sl,bnd%c%rho_ice,bnd%c%rho_sw)

        ! Calculate the surface elevation to be consistent with current H_ice field
        call calc_z_srf_max(tpo%now%z_srf,tpo%now%H_ice,tpo%now%f_ice,bnd%z_bed,bnd%z_sl,bnd%c%rho_ice,bnd%c%rho_sw)
        
        ! Calculate the ice base elevation too
        ! Define z_base as the elevation at the base of the ice sheet 
        ! This is used for the basal derivative instead of bedrock so
        ! that it is valid for both grounded and floating ice. Note, 
        ! for grounded ice, z_base==z_bed.  
        tpo%now%z_base = tpo%now%z_srf - tpo%now%H_ice 

        ! 2. Calculate additional topographic properties ------------------

        ! Calculate the ice thickness gradient (on staggered acx/y nodes)
        !call calc_gradient_ac(tpo%now%dHidx,tpo%now%dHidy,tpo%now%H_ice,tpo%par%dx)
        ! call calc_gradient_ac_ice(tpo%now%dHidx,tpo%now%dHidy,tpo%now%H_ice,tpo%now%f_ice,tpo%par%dx, &
        !                                         tpo%par%margin2nd,tpo%par%grad_lim,tpo%par%boundaries,zero_outside=.TRUE.)
        
        ! Calculate the surface slope
        ! call calc_gradient_ac(tpo%now%dzsdx,tpo%now%dzsdy,tpo%now%z_srf,tpo%par%dx)


        ! New routines 
        call calc_gradient_acx(tpo%now%dzsdx,tpo%now%z_srf,tpo%now%f_ice,tpo%par%dx,tpo%par%grad_lim,tpo%par%margin2nd,zero_outside=.FALSE.,boundaries=tpo%par%boundaries)
        call calc_gradient_acy(tpo%now%dzsdy,tpo%now%z_srf,tpo%now%f_ice,tpo%par%dy,tpo%par%grad_lim,tpo%par%margin2nd,zero_outside=.FALSE.,boundaries=tpo%par%boundaries)
        
        call calc_gradient_acx(tpo%now%dHidx,tpo%now%H_ice,tpo%now%f_ice,tpo%par%dx,tpo%par%grad_lim,tpo%par%margin2nd,zero_outside=.TRUE.,boundaries=tpo%par%boundaries)
        call calc_gradient_acy(tpo%now%dHidy,tpo%now%H_ice,tpo%now%f_ice,tpo%par%dy,tpo%par%grad_lim,tpo%par%margin2nd,zero_outside=.TRUE.,boundaries=tpo%par%boundaries)
        
        call calc_gradient_acx(tpo%now%dzbdx,tpo%now%z_base,tpo%now%f_ice,tpo%par%dx,tpo%par%grad_lim,tpo%par%margin2nd,zero_outside=.FALSE.,boundaries=tpo%par%boundaries)
        call calc_gradient_acy(tpo%now%dzbdy,tpo%now%z_base,tpo%now%f_ice,tpo%par%dy,tpo%par%grad_lim,tpo%par%margin2nd,zero_outside=.FALSE.,boundaries=tpo%par%boundaries)
        
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
                                                                tpo%now%H_grnd,tpo%par%gz_nx,tpo%par%boundaries)
            
            case(3) 
                ! Grounded area using analytical solutions of Leguy et al. (2021)

                call determine_grounded_fractions(tpo%now%f_grnd,tpo%now%f_grnd_acx,tpo%now%f_grnd_acy, &
                                                                tpo%now%f_grnd_ab,tpo%now%H_grnd,tpo%par%boundaries)

        end select
        
        ! Calculate grounded fraction due to pinning points 
        call calc_f_grnd_pinning_points(tpo%now%f_grnd_pin,tpo%now%H_ice,tpo%now%f_ice, &
                                                bnd%z_bed,bnd%z_bed_sd,bnd%z_sl,bnd%c%rho_ice,bnd%c%rho_sw)

        ! Calculate the grounding-line distance
        call calc_distance_to_grounding_line(tpo%now%dist_grline,tpo%now%f_grnd,tpo%par%dx,tpo%par%boundaries)

        ! Define the grounding-zone mask too 
        call calc_grounding_line_zone(tpo%now%mask_grz,tpo%now%dist_grline,tpo%par%dist_grz)

        ! Calculate distance to the ice margin
        call calc_distance_to_ice_margin(tpo%now%dist_margin,tpo%now%f_ice,tpo%par%dx,tpo%par%boundaries)

        ! Calculate the general bed mask
        call gen_mask_bed(tpo%now%mask_bed,tpo%now%f_ice,thrm%now%f_pmp, &
                                            tpo%now%f_grnd,tpo%now%mask_grz)


        ! Calculate the ice-front mask (mainly for use in dynamics)
        call calc_ice_front(tpo%now%mask_frnt,tpo%now%f_ice,tpo%now%f_grnd,bnd%z_bed,bnd%z_sl,tpo%par%boundaries)

        ! Determine ice thickness for use exclusively with the dynamics solver
        select case(trim(dyn%par%ssa_lat_bc))

            case("slab")
                ! Calculate extended ice thickness fields over whole domain
                ! for use with dynamics solver.
                ! Note: should not be used with MISMIP, TROUGH, etc. 

                tpo%now%H_ice_dyn = tpo%now%H_ice
                where (tpo%now%f_ice .lt. 1.0) tpo%now%H_ice_dyn = 1.0_wp
                
                ! Calculate the ice fraction mask for use with the dynamics solver
                call calc_ice_fraction(tpo%now%f_ice_dyn,tpo%now%H_ice_dyn,bnd%z_bed,bnd%z_sl,bnd%c%rho_ice, &
                                bnd%c%rho_sw,tpo%par%boundaries,flt_subgrid=.FALSE.)

            case("slab-ext")
                ! Calculate extended ice thickness fields n_ext points
                ! away from grounded margin. 

                tpo%now%H_ice_dyn = tpo%now%H_ice
                where(tpo%now%H_ice_dyn .gt. 0.0 .and. tpo%now%H_ice_dyn .lt. 1.0) &
                        tpo%now%H_ice_dyn = 1.0_wp 

                call extend_floating_slab(tpo%now%H_ice_dyn,tpo%now%f_grnd,H_slab=1.0_wp,n_ext=4)

                ! Calculate the ice fraction mask for use with the dynamics solver
                call calc_ice_fraction(tpo%now%f_ice_dyn,tpo%now%H_ice_dyn,bnd%z_bed,bnd%z_sl,bnd%c%rho_ice, &
                                bnd%c%rho_sw,tpo%par%boundaries,flt_subgrid=.FALSE.)

            case DEFAULT 
                ! No modification of ice thickness for dynamics solver 
                ! Set standard ice thickness field for use with dynamics 
            
                tpo%now%H_ice_dyn = tpo%now%H_ice
                tpo%now%f_ice_dyn = tpo%now%f_ice 

        end select
        
        return 

    end subroutine calc_ytopo_diagnostic

    subroutine calc_ytopo_rates(tpo,bnd,time,dt,step,check_mb,overwrite)
        ! Calculate average rates over outer timestep

        implicit none

        type(ytopo_class),  intent(INOUT) :: tpo 
        type(ybound_class), intent(IN)    :: bnd
        real(wp),           intent(IN)    :: time
        real(wp),           intent(IN)    :: dt 
        character(len=*),   intent(IN)    :: step 
        logical,            intent(IN)    :: check_mb 
        logical, optional,  intent(IN)    :: overwrite 
        
        real(wp), parameter :: tol_dt = 1e-3

        select case(trim(step))

            case("init")
                ! Initialization of averaging fields - set to zero

                tpo%now%rates%dzsdt         = 0.0
                tpo%now%rates%dHidt         = 0.0
                tpo%now%rates%dHidt_dyn     = 0.0
                tpo%now%rates%mb_net        = 0.0
                tpo%now%rates%mb_relax      = 0.0
                tpo%now%rates%mb_resid      = 0.0
                tpo%now%rates%mb_err        = 0.0
                tpo%now%rates%smb           = 0.0
                tpo%now%rates%bmb           = 0.0
                tpo%now%rates%fmb           = 0.0
                tpo%now%rates%dmb           = 0.0
                tpo%now%rates%cmb           = 0.0
                tpo%now%rates%cmb_flt       = 0.0
                tpo%now%rates%cmb_grnd      = 0.0

                tpo%now%rates%dt_tot = 0.0 

            case("step")
                ! Add current step to total

                tpo%now%rates%dzsdt         = tpo%now%rates%dzsdt       + tpo%now%dzsdt*dt
                tpo%now%rates%dHidt         = tpo%now%rates%dHidt       + tpo%now%dHidt*dt
                tpo%now%rates%dHidt_dyn     = tpo%now%rates%dHidt_dyn   + tpo%now%dHidt_dyn*dt
                tpo%now%rates%mb_net        = tpo%now%rates%mb_net      + tpo%now%mb_net*dt
                tpo%now%rates%mb_relax      = tpo%now%rates%mb_relax    + tpo%now%mb_relax*dt
                tpo%now%rates%mb_resid      = tpo%now%rates%mb_resid    + tpo%now%mb_resid*dt
                tpo%now%rates%mb_err        = tpo%now%rates%mb_err      + tpo%now%mb_err*dt
                tpo%now%rates%smb           = tpo%now%rates%smb         + tpo%now%smb*dt
                tpo%now%rates%bmb           = tpo%now%rates%bmb         + tpo%now%bmb*dt
                tpo%now%rates%fmb           = tpo%now%rates%fmb         + tpo%now%fmb*dt
                tpo%now%rates%dmb           = tpo%now%rates%dmb         + tpo%now%dmb*dt
                tpo%now%rates%cmb           = tpo%now%rates%cmb         + tpo%now%cmb*dt
                tpo%now%rates%cmb_flt       = tpo%now%rates%cmb_flt     + tpo%now%cmb_flt*dt
                tpo%now%rates%cmb_grnd      = tpo%now%rates%cmb_grnd    + tpo%now%cmb_grnd*dt

                tpo%now%rates%dt_tot = tpo%now%rates%dt_tot + dt  
                
            case("final")
                ! Divide by total time to get average rate

                if (tpo%now%rates%dt_tot .gt. 0.0) then 

                    tpo%now%rates%dzsdt         = tpo%now%rates%dzsdt / tpo%now%rates%dt_tot
                    tpo%now%rates%dHidt         = tpo%now%rates%dHidt / tpo%now%rates%dt_tot
                    tpo%now%rates%dHidt_dyn     = tpo%now%rates%dHidt_dyn / tpo%now%rates%dt_tot
                    tpo%now%rates%mb_net        = tpo%now%rates%mb_net / tpo%now%rates%dt_tot
                    tpo%now%rates%mb_relax      = tpo%now%rates%mb_relax / tpo%now%rates%dt_tot
                    tpo%now%rates%mb_resid      = tpo%now%rates%mb_resid / tpo%now%rates%dt_tot
                    tpo%now%rates%mb_err        = tpo%now%rates%mb_err / tpo%now%rates%dt_tot
                    tpo%now%rates%smb           = tpo%now%rates%smb / tpo%now%rates%dt_tot
                    tpo%now%rates%bmb           = tpo%now%rates%bmb / tpo%now%rates%dt_tot
                    tpo%now%rates%fmb           = tpo%now%rates%fmb / tpo%now%rates%dt_tot
                    tpo%now%rates%dmb           = tpo%now%rates%dmb / tpo%now%rates%dt_tot
                    tpo%now%rates%cmb           = tpo%now%rates%cmb / tpo%now%rates%dt_tot
                    tpo%now%rates%cmb_flt       = tpo%now%rates%cmb_flt / tpo%now%rates%dt_tot
                    tpo%now%rates%cmb_grnd      = tpo%now%rates%cmb_grnd / tpo%now%rates%dt_tot

                    ! Check that dt_tot matches outer dt value
                    if ( abs(dt - tpo%now%rates%dt_tot) .gt. tol_dt) then
                        write(*,*) "calc_ytopo_rates: dt, dt_tot : ", dt, tpo%now%rates%dt_tot
                    end if 

                end if

                if (present(overwrite)) then
                    if (overwrite) then
                        ! Overwrite the instantaneous rates with averaged rates for output

                        tpo%now%dzsdt       = tpo%now%rates%dzsdt
                        tpo%now%dHidt       = tpo%now%rates%dHidt
                        tpo%now%dHidt_dyn   = tpo%now%rates%dHidt_dyn
                        tpo%now%mb_net      = tpo%now%rates%mb_net
                        tpo%now%mb_relax    = tpo%now%rates%mb_relax
                        tpo%now%mb_resid    = tpo%now%rates%mb_resid
                        tpo%now%mb_err      = tpo%now%rates%mb_err
                        tpo%now%smb         = tpo%now%rates%smb
                        tpo%now%bmb         = tpo%now%rates%bmb
                        tpo%now%fmb         = tpo%now%rates%fmb
                        tpo%now%dmb         = tpo%now%rates%dmb
                        tpo%now%cmb         = tpo%now%rates%cmb
                        tpo%now%cmb_flt     = tpo%now%rates%cmb_flt
                        tpo%now%cmb_grnd    = tpo%now%rates%cmb_grnd

                    end if
                end if 

            case DEFAULT 

                write(io_unit_err,*) "calc_ytopo_rates:: Error: step name not recognized."
                write(io_unit_err,*) "step = ", trim(step)
                stop 

        end select

        if (check_mb) then 
            ! Perform mass balance check to make sure that mass is conserved

            call check_mass_conservation(tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd,tpo%now%dHidt, &
                        tpo%now%mb_net,tpo%now%cmb,tpo%now%dHidt_dyn,tpo%now%smb,tpo%now%bmb, &
                        tpo%now%fmb,tpo%now%dmb,tpo%now%mb_resid,tpo%par%dx,bnd%c%sec_year,time,dt, &
                        units="km^3/yr",label=step)
                        
        end if 

        return

    end subroutine calc_ytopo_rates
    
    subroutine ytopo_par_load(par,filename,group,nx,ny,dx,init)

        type(ytopo_param_class), intent(OUT) :: par
        character(len=*),        intent(IN)  :: filename
        character(len=*),        intent(IN)  :: group       ! Usually "ytopo"
        integer,                 intent(IN)  :: nx, ny 
        real(wp),                intent(IN)  :: dx  
        logical, optional,       intent(IN)  :: init 

        ! Local variables
        logical :: init_pars 

        init_pars = .FALSE.
        if (present(init)) init_pars = .TRUE. 
 
        ! Store parameter values in output object
        call nml_read(filename,group,"solver",            par%solver,           init=init_pars)
        call nml_read(filename,group,"surf_gl_method",    par%surf_gl_method,   init=init_pars)
        call nml_read(filename,group,"calv_flt_method",   par%calv_flt_method,  init=init_pars)
        call nml_read(filename,group,"calv_grnd_method",  par%calv_grnd_method, init=init_pars)
        call nml_read(filename,group,"bmb_gl_method",     par%bmb_gl_method,    init=init_pars)
        call nml_read(filename,group,"fmb_method",        par%fmb_method,       init=init_pars)
        call nml_read(filename,group,"dmb_method",        par%dmb_method,       init=init_pars)
        call nml_read(filename,group,"margin2nd",         par%margin2nd,        init=init_pars)
        call nml_read(filename,group,"margin_flt_subgrid",par%margin_flt_subgrid,init=init_pars)
        call nml_read(filename,group,"use_bmb",           par%use_bmb,          init=init_pars)
        call nml_read(filename,group,"topo_fixed",        par%topo_fixed,       init=init_pars)
        call nml_read(filename,group,"topo_rel",          par%topo_rel,         init=init_pars)
        call nml_read(filename,group,"topo_rel_tau",      par%topo_rel_tau,     init=init_pars)
        call nml_read(filename,group,"topo_rel_field",    par%topo_rel_field,   init=init_pars)
        call nml_read(filename,group,"calv_tau",          par%calv_tau,         init=init_pars)
        call nml_read(filename,group,"calv_thin",         par%calv_thin,        init=init_pars)
        call nml_read(filename,group,"H_min_grnd",        par%H_min_grnd,       init=init_pars)
        call nml_read(filename,group,"H_min_flt",         par%H_min_flt,        init=init_pars)
        call nml_read(filename,group,"sd_min",            par%sd_min,           init=init_pars)
        call nml_read(filename,group,"sd_max",            par%sd_max,           init=init_pars)
        call nml_read(filename,group,"calv_grnd_max",     par%calv_grnd_max,    init=init_pars)
        call nml_read(filename,group,"grad_lim",          par%grad_lim,         init=init_pars)
        call nml_read(filename,group,"grad_lim_zb",       par%grad_lim_zb,      init=init_pars)
        call nml_read(filename,group,"dist_grz",          par%dist_grz,         init=init_pars)
        call nml_read(filename,group,"gl_sep",            par%gl_sep,           init=init_pars)
        call nml_read(filename,group,"gz_nx",             par%gz_nx,            init=init_pars)
        call nml_read(filename,group,"gz_Hg0",            par%gz_Hg0,           init=init_pars)
        call nml_read(filename,group,"gz_Hg1",            par%gz_Hg1,           init=init_pars)
        call nml_read(filename,group,"fmb_scale",         par%fmb_scale,        init=init_pars)
        call nml_read(filename,group,"k2",                par%k2,               init=init_pars)
        call nml_read(filename,group,"w2",                par%w2,               init=init_pars)
        call nml_read(filename,group,"kt_ref",            par%kt_ref,           init=init_pars)
        call nml_read(filename,group,"kt_deep",           par%kt_deep,          init=init_pars)
        call nml_read(filename,group,"Hc_ref",            par%Hc_ref,           init=init_pars)
        call nml_read(filename,group,"Hc_ref_thin",       par%Hc_ref_thin,      init=init_pars)
        call nml_read(filename,group,"Hc_deep",           par%Hc_deep,          init=init_pars)
        call nml_read(filename,group,"zb_deep_0",         par%zb_deep_0,        init=init_pars)
        call nml_read(filename,group,"zb_deep_1",         par%zb_deep_1,        init=init_pars)
        call nml_read(filename,group,"zb_sigma",          par%zb_sigma,         init=init_pars)
        call nml_read(filename,group,"dmb_alpha_max",     par%dmb_alpha_max,    init=init_pars)
        call nml_read(filename,group,"dmb_tau",           par%dmb_tau,          init=init_pars)
        call nml_read(filename,group,"dmb_sigma_ref",     par%dmb_sigma_ref,    init=init_pars)
        call nml_read(filename,group,"dmb_m_d",           par%dmb_m_d,          init=init_pars)
        call nml_read(filename,group,"dmb_m_r",           par%dmb_m_r,          init=init_pars)
        
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

        call ytopo_pc_alloc(now%pred,nx,ny)
        call ytopo_pc_alloc(now%corr,nx,ny)

        ! Rates (for timestep averages)
        allocate(now%rates%dzsdt(nx,ny))
        allocate(now%rates%dHidt(nx,ny))
        allocate(now%rates%dHidt_dyn(nx,ny))
        allocate(now%rates%mb_net(nx,ny))
        allocate(now%rates%mb_relax(nx,ny))
        allocate(now%rates%mb_resid(nx,ny))
        allocate(now%rates%mb_err(nx,ny))
        allocate(now%rates%smb(nx,ny))
        allocate(now%rates%bmb(nx,ny))
        allocate(now%rates%fmb(nx,ny))
        allocate(now%rates%dmb(nx,ny))
        allocate(now%rates%cmb(nx,ny))
        allocate(now%rates%cmb_flt(nx,ny))
        allocate(now%rates%cmb_grnd(nx,ny))
        
        ! Remaining ytopo fields...

        allocate(now%H_ice(nx,ny))
        allocate(now%z_srf(nx,ny))
        allocate(now%z_base(nx,ny))

        allocate(now%dzsdt(nx,ny))
        allocate(now%dHidt(nx,ny))
        allocate(now%dHidt_dyn(nx,ny))
        allocate(now%mb_net(nx,ny))
        allocate(now%mb_relax(nx,ny))
        allocate(now%mb_resid(nx,ny))
        allocate(now%mb_err(nx,ny))
        allocate(now%smb(nx,ny))
        allocate(now%bmb(nx,ny))
        allocate(now%fmb(nx,ny))
        allocate(now%dmb(nx,ny))
        allocate(now%cmb(nx,ny))
        allocate(now%cmb_flt(nx,ny))
        allocate(now%cmb_grnd(nx,ny))
        
        allocate(now%bmb_ref(nx,ny))
        allocate(now%fmb_ref(nx,ny))
        allocate(now%dmb_ref(nx,ny))

        allocate(now%mask_adv(nx,ny))
        
        allocate(now%eps_eff(nx,ny))
        allocate(now%tau_eff(nx,ny))
        
        allocate(now%dzsdx(nx,ny))
        allocate(now%dzsdy(nx,ny))
        allocate(now%dHidx(nx,ny))
        allocate(now%dHidy(nx,ny))
        allocate(now%dzbdx(nx,ny))
        allocate(now%dzbdy(nx,ny))

        allocate(now%H_eff(nx,ny))
        allocate(now%H_grnd(nx,ny))
        allocate(now%H_calv(nx,ny))
        allocate(now%kt(nx,ny))
        allocate(now%z_bed_filt(nx,ny))

        ! Masks 
        allocate(now%f_grnd(nx,ny))
        allocate(now%f_grnd_acx(nx,ny))
        allocate(now%f_grnd_acy(nx,ny))
        allocate(now%f_grnd_ab(nx,ny))
        allocate(now%f_ice(nx,ny))

        allocate(now%f_grnd_bmb(nx,ny))
        allocate(now%f_grnd_pin(nx,ny))

        allocate(now%dist_margin(nx,ny))
        allocate(now%dist_grline(nx,ny))

        allocate(now%mask_bed(nx,ny))
        allocate(now%mask_grz(nx,ny))
        allocate(now%mask_frnt(nx,ny))
        
        allocate(now%dHidt_dyn_n(nx,ny))
        allocate(now%H_ice_n(nx,ny))
        allocate(now%z_srf_n(nx,ny))
        
        allocate(now%H_ice_dyn(nx,ny))
        allocate(now%f_ice_dyn(nx,ny))

        allocate(now%tau_relax(nx,ny))

        now%rates%dzsdt         = 0.0
        now%rates%dHidt         = 0.0
        now%rates%dHidt_dyn     = 0.0
        now%rates%mb_net        = 0.0
        now%rates%mb_relax      = 0.0
        now%rates%mb_resid      = 0.0
        now%rates%mb_err        = 0.0
        now%rates%smb           = 0.0
        now%rates%bmb           = 0.0
        now%rates%fmb           = 0.0
        now%rates%dmb           = 0.0
        now%rates%cmb           = 0.0
        now%rates%cmb_flt       = 0.0
        now%rates%cmb_grnd      = 0.0
        
        now%H_ice       = 0.0 
        now%z_srf       = 0.0
        now%z_base      = 0.0  
        now%dzsdt       = 0.0 
        now%dHidt       = 0.0
        now%dHidt_dyn   = 0.0
        now%mb_net      = 0.0 
        now%mb_relax    = 0.0
        now%mb_resid    = 0.0
        now%mb_err      = 0.0
        now%smb         = 0.0 
        now%bmb         = 0.0  
        now%fmb         = 0.0
        now%dmb         = 0.0
        now%cmb         = 0.0
        now%cmb_flt     = 0.0
        now%cmb_grnd    = 0.0
        
        now%bmb_ref     = 0.0  
        now%fmb_ref     = 0.0
        now%dmb_ref     = 0.0
        
        now%mask_adv    = 0

        now%eps_eff     = 0.0
        now%tau_eff     = 0.0
        
        now%dzsdx       = 0.0 
        now%dzsdy       = 0.0 
        now%dHidx       = 0.0 
        now%dHidy       = 0.0
        now%dzbdx       = 0.0 
        now%dzbdy       = 0.0 
        now%H_eff       = 0.0 
        now%H_grnd      = 0.0  
        now%H_calv      = 0.0  
        now%kt          = 0.0  
        now%z_bed_filt  = 0.0  

        now%f_grnd      = 0.0  
        now%f_grnd_acx  = 0.0  
        now%f_grnd_acy  = 0.0  
        now%f_grnd_ab   = 0.0
        now%f_grnd_bmb  = 0.0
        now%f_grnd_pin  = 0.0
        now%f_ice       = 0.0  
        now%dist_margin = 0.0
        now%dist_grline = 0.0 

        now%mask_bed    = 0 
        now%mask_grz    = 0 
        now%mask_frnt   = 0

        now%dHidt_dyn_n = 0.0  
        now%H_ice_n     = 0.0 
        now%z_srf_n     = 0.0 

        now%H_ice_dyn   = 0.0 
        now%f_ice_dyn   = 0.0 
        
        now%tau_relax   = 0.0

        return 

    end subroutine ytopo_alloc

    subroutine ytopo_dealloc(now)

        implicit none 

        type(ytopo_state_class), intent(INOUT) :: now

        call ytopo_pc_dealloc(now%pred)
        call ytopo_pc_dealloc(now%corr)

        if (allocated(now%rates%dzsdt))         deallocate(now%rates%dzsdt)
        if (allocated(now%rates%dHidt))         deallocate(now%rates%dHidt)
        if (allocated(now%rates%dHidt_dyn))     deallocate(now%rates%dHidt_dyn)
        if (allocated(now%rates%mb_net))        deallocate(now%rates%mb_net)
        if (allocated(now%rates%mb_relax))      deallocate(now%rates%mb_relax)
        if (allocated(now%rates%mb_resid))      deallocate(now%rates%mb_resid)
        if (allocated(now%rates%mb_err))        deallocate(now%rates%mb_err)
        if (allocated(now%rates%smb))           deallocate(now%rates%smb)
        if (allocated(now%rates%bmb))           deallocate(now%rates%bmb)
        if (allocated(now%rates%fmb))           deallocate(now%rates%fmb)
        if (allocated(now%rates%dmb))           deallocate(now%rates%dmb)
        if (allocated(now%rates%cmb))           deallocate(now%rates%cmb)
        if (allocated(now%rates%cmb_flt))       deallocate(now%rates%cmb_flt)
        if (allocated(now%rates%cmb_grnd))      deallocate(now%rates%cmb_grnd)
        
        if (allocated(now%H_ice))       deallocate(now%H_ice)
        if (allocated(now%z_srf))       deallocate(now%z_srf)
        if (allocated(now%z_base))      deallocate(now%z_base)

        if (allocated(now%dzsdt))       deallocate(now%dzsdt)
        if (allocated(now%dHidt))       deallocate(now%dHidt)
        if (allocated(now%dHidt_dyn))   deallocate(now%dHidt_dyn)
        if (allocated(now%mb_net))      deallocate(now%mb_net)
        if (allocated(now%mb_relax))    deallocate(now%mb_relax)
        if (allocated(now%mb_resid))    deallocate(now%mb_resid)
        if (allocated(now%mb_err))      deallocate(now%mb_err)
        if (allocated(now%smb))         deallocate(now%smb)
        if (allocated(now%bmb))         deallocate(now%bmb)
        if (allocated(now%fmb))         deallocate(now%fmb)
        if (allocated(now%dmb))         deallocate(now%dmb)
        if (allocated(now%cmb))         deallocate(now%cmb)
        if (allocated(now%cmb_flt))     deallocate(now%cmb_flt)
        if (allocated(now%cmb_grnd))    deallocate(now%cmb_grnd)
        
        if (allocated(now%bmb_ref))     deallocate(now%bmb_ref)
        if (allocated(now%fmb_ref))     deallocate(now%fmb_ref)
        if (allocated(now%dmb_ref))     deallocate(now%dmb_ref)
        
        if (allocated(now%mask_adv))    deallocate(now%mask_adv)
        
        if (allocated(now%eps_eff))     deallocate(now%eps_eff)
        if (allocated(now%tau_eff))     deallocate(now%tau_eff)
        
        if (allocated(now%dzsdx))       deallocate(now%dzsdx)
        if (allocated(now%dzsdy))       deallocate(now%dzsdy)
        if (allocated(now%dHidx))       deallocate(now%dHidx)
        if (allocated(now%dHidy))       deallocate(now%dHidy)
        if (allocated(now%dzbdx))       deallocate(now%dzbdx)
        if (allocated(now%dzbdy))       deallocate(now%dzbdy)
        
        if (allocated(now%H_eff))       deallocate(now%H_eff)
        if (allocated(now%H_grnd))      deallocate(now%H_grnd)
        if (allocated(now%H_calv))      deallocate(now%H_calv)
        if (allocated(now%kt))          deallocate(now%kt)
        if (allocated(now%z_bed_filt))  deallocate(now%z_bed_filt)

        if (allocated(now%f_grnd))      deallocate(now%f_grnd)
        if (allocated(now%f_grnd_acx))  deallocate(now%f_grnd_acx)
        if (allocated(now%f_grnd_acy))  deallocate(now%f_grnd_acy)
        if (allocated(now%f_grnd_ab))   deallocate(now%f_grnd_ab)
        if (allocated(now%f_grnd_bmb))  deallocate(now%f_grnd_bmb)
        if (allocated(now%f_grnd_pin))  deallocate(now%f_grnd_pin)

        if (allocated(now%f_ice))       deallocate(now%f_ice)

        if (allocated(now%dist_margin)) deallocate(now%dist_margin)
        if (allocated(now%dist_grline)) deallocate(now%dist_grline)

        if (allocated(now%mask_bed))    deallocate(now%mask_bed)
        if (allocated(now%mask_grz))    deallocate(now%mask_grz)
        if (allocated(now%mask_frnt))   deallocate(now%mask_frnt)
        
        if (allocated(now%dHidt_dyn_n)) deallocate(now%dHidt_dyn_n)
        if (allocated(now%H_ice_n))     deallocate(now%H_ice_n)
        if (allocated(now%z_srf_n))     deallocate(now%z_srf_n)
        
        if (allocated(now%H_ice_dyn))   deallocate(now%H_ice_dyn)
        if (allocated(now%f_ice_dyn))   deallocate(now%f_ice_dyn)
        
        if (allocated(now%tau_relax))   deallocate(now%tau_relax)
        
        return 

    end subroutine ytopo_dealloc
    
    subroutine ytopo_pc_alloc(pc,nx,ny)

        implicit none

        type(ytopo_pc_class), intent(INOUT) :: pc 
        integer, intent(IN) :: nx, ny  

        ! First deallocate everything for safety
        call ytopo_pc_dealloc(pc)

        ! Allocate fields 
        allocate(pc%H_ice(nx,ny))
        allocate(pc%dHidt_dyn(nx,ny))
        allocate(pc%mb_net(nx,ny))
        allocate(pc%mb_relax(nx,ny))
        allocate(pc%mb_resid(nx,ny))
        allocate(pc%smb(nx,ny))
        allocate(pc%bmb(nx,ny))
        allocate(pc%fmb(nx,ny))
        allocate(pc%dmb(nx,ny))
        allocate(pc%cmb(nx,ny))      
        allocate(pc%cmb_flt(nx,ny))
        allocate(pc%cmb_grnd(nx,ny))
        
        ! Initialize to zero
        pc%H_ice        = 0.0
        pc%dHidt_dyn    = 0.0
        pc%mb_net       = 0.0
        pc%mb_relax     = 0.0
        pc%mb_resid     = 0.0
        pc%smb          = 0.0
        pc%bmb          = 0.0
        pc%fmb          = 0.0
        pc%dmb          = 0.0
        pc%cmb          = 0.0      
        pc%cmb_flt      = 0.0
        pc%cmb_grnd     = 0.0
        
        return

    end subroutine ytopo_pc_alloc

    subroutine ytopo_pc_dealloc(pc)

        implicit none

        type(ytopo_pc_class), intent(INOUT) :: pc 
        
        if (allocated(pc%H_ice))        deallocate(pc%H_ice)
        if (allocated(pc%dHidt_dyn))    deallocate(pc%dHidt_dyn)
        if (allocated(pc%mb_net))       deallocate(pc%mb_net)
        if (allocated(pc%mb_relax))     deallocate(pc%mb_relax)
        if (allocated(pc%mb_resid))     deallocate(pc%mb_resid)
        if (allocated(pc%smb))          deallocate(pc%smb)
        if (allocated(pc%bmb))          deallocate(pc%bmb)
        if (allocated(pc%fmb))          deallocate(pc%fmb)
        if (allocated(pc%cmb))          deallocate(pc%cmb)
        if (allocated(pc%cmb_flt))      deallocate(pc%cmb_flt)
        if (allocated(pc%cmb_grnd))     deallocate(pc%cmb_grnd)
        
        return

    end subroutine ytopo_pc_dealloc

end module yelmo_topography
