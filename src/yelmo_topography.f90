
module yelmo_topography

    use nml  
    use ncio
    
    use yelmo_defs
    use yelmo_tools 
    
    use mass_conservation
    use calving_aa
    use calving_ac
    use lsf_module
    use topography 
    use discharge

    use runge_kutta 
    use derivatives
    use distances

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
    
    subroutine calc_ytopo_pc(tpo,dyn,mat,thrm,bnd,dta,time,topo_fixed,pc_step,use_H_pred)

        implicit none 

        type(ytopo_class),  intent(INOUT) :: tpo
        type(ydyn_class),   intent(IN)    :: dyn
        type(ymat_class),   intent(INOUT)    :: mat
        type(ytherm_class), intent(IN)    :: thrm  
        type(ybound_class), intent(IN)    :: bnd
        type(ydata_class),  intent(IN)    :: dta 
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
                    tpo%now%lsf_n       = tpo%now%lsf

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
                    ! Limit dynamic rate of change for stability (typically < 100 m/yr)
                    tpo%now%H_ice = tpo%now%H_ice_n
                    call apply_tendency(tpo%now%H_ice,tpo%now%dHidt_dyn,dt,"dyn_pred",adjust_mb=.TRUE.,mb_lim=tpo%par%dHdt_dyn_lim)

                case("corrector") 

                    ! Set current thickness to predicted thickness
                    tpo%now%H_ice = tpo%now%pred%H_ice 
                    tpo%now%lsf   = tpo%now%pred%lsf

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
                    ! Limit dynamic rate of change for stability (typically < 100 m/yr)
                    tpo%now%H_ice = tpo%now%H_ice_n
                    tpo%now%lsf   = tpo%now%lsf_n
                    call apply_tendency(tpo%now%H_ice,tpo%now%dHidt_dyn,dt,"dyn_corr",adjust_mb=.TRUE.,mb_lim=tpo%par%dHdt_dyn_lim)
                    
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
                                        tpo%par%gz_nx,tpo%par%bmb_gl_method,dta%pd%mask_bed,tpo%par%boundaries)

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

                    ! === calving ===
                    ! Calculate and apply calving
                    if (tpo%par%use_lsf) then
                        ! Level-set function as calving
                        call calc_ytopo_calving_lsf(tpo,dyn,mat,thrm,bnd,dt,H_prev,time)
                    else
                        ! Mass balance calving
                        call calc_ytopo_calving(tpo,dyn,mat,thrm,bnd,dt)
                    end if

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
                    tpo%now%pred%lsf        = tpo%now%lsf 
                    
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
                    tpo%now%corr%lsf        = tpo%now%lsf
                    
                    ! Restore main ice thickness field to original 
                    ! value at the beginning of the timestep for 
                    ! calculation of remaining quantities (thermo, material)
                    tpo%now%H_ice = tpo%now%H_ice_n 
                    tpo%now%lsf   = tpo%now%lsf_n 

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
                        tpo%now%lsf         = tpo%now%pred%lsf 
                        
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
                        tpo%now%lsf         = tpo%now%corr%lsf 
                        
                    end if
                    
            end select

            ! Determine rates of change
            tpo%now%dHidt  = (tpo%now%H_ice - tpo%now%H_ice_n) / dt 
            tpo%now%dzsdt  = (tpo%now%z_srf - tpo%now%z_srf_n) / dt
            tpo%now%dlsfdt = (tpo%now%lsf   - tpo%now%lsf_n) / dt

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

        call define_calving_thickness_threshold(tpo%now%H_calv,tpo%now%z_bed_filt,tpo%par%Hc_ref_flt, &
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

    subroutine calc_ytopo_calving_lsf(tpo,dyn,mat,thrm,bnd,dt,H_prev,time_now)
        ! Calving computed as a flux. LSF mask is updated with velocity - calving front velocity.
        ! Points in ocean domain (LSF > 0) will be deleted with a melt equal to the ice thickness.

        implicit none
    
        type(ytopo_class),  intent(INOUT) :: tpo
        type(ydyn_class),   intent(IN)    :: dyn
        type(ymat_class),   intent(INOUT)    :: mat
        type(ytherm_class), intent(IN)    :: thrm
        type(ybound_class), intent(IN)    :: bnd
        real(wp),           intent(IN)    :: dt
        real(wp),           intent(IN)    :: H_prev(:,:)
        real(wp), optional, intent(IN)    :: time_now
    
        ! Local variables
        integer  :: i, j, nx, ny
        integer  :: im1, ip1, jm1, jp1
        real(wp) :: dt_kill
        real(wp), allocatable :: mbal_now(:,:)
        !real(wp), allocatable :: u_acx_fill(:,:), v_acy_fill(:,:)
        integer  :: BC

        ! Make sure dt is not zero
        dt_kill = dt 
        if (dt_kill .eq. 0.0) dt_kill = 1.0_wp
    
        nx = size(tpo%now%H_ice,1)
        ny = size(tpo%now%H_ice,2)

        ! Set boundary condition code
        BC = boundary_code(tpo%par%boundaries)

        allocate(mbal_now(nx,ny))

        ! === Floating calving laws ===
        
        ! Initialize the calving rates
        tpo%now%cmb_flt_x = 0.0_wp
        tpo%now%cmb_flt_y = 0.0_wp
        tpo%now%cmb_flt   = 0.0_wp

        select case(trim(tpo%par%calv_flt_method))
    
            case("zero","none")
                ! Do nothing. No calving.
                
            case("equil")
                ! For an equilibrated ice sheet calving rates should be opposite to ice velocity
                tpo%now%cmb_flt_x = -1*dyn%now%ux_bar
                tpo%now%cmb_flt_y = -1*dyn%now%uy_bar
    
            case("threshold")
                call calc_calving_threshold_lsf(tpo%now%cmb_flt_x,tpo%now%cmb_flt_y,dyn%now%ux_bar,dyn%now%uy_bar,tpo%now%H_ice,tpo%par%Hc_ref_flt,tpo%now%f_ice,tpo%par%boundaries)
        
            case("vm-m16")
                call calc_calving_rate_vonmises_m16(tpo%now%cmb_flt_x,tpo%now%cmb_flt_y,dyn%now%ux_bar,dyn%now%uy_bar,mat%now%strs2D%tau_eig_1,tpo%par%tau_ice,tpo%now%f_ice,tpo%par%boundaries)
                
            ! TO DO: Add new laws
    
            ! === CalvMIP laws ===
            case("exp1","exp3")
                call calvmip_exp1(tpo%now%cmb_flt_x,tpo%now%cmb_flt_y,dyn%now%ux_bar,dyn%now%uy_bar,tpo%now%lsf,tpo%par%dx,tpo%par%boundaries) 
    
            case("exp2","exp4")
                call calvmip_exp2(tpo%now%cmb_flt_x,tpo%now%cmb_flt_y,dyn%now%ux_bar,dyn%now%uy_bar,time_now,tpo%par%boundaries)
            
            case("exp5")
                call calvmip_exp5_aa(tpo%now%cmb_flt_x,tpo%now%cmb_flt_y,dyn%now%ux_bar,dyn%now%uy_bar,tpo%now%H_ice,tpo%par%Hc_ref_flt,tpo%now%f_ice,tpo%par%boundaries)

            case DEFAULT
    
                write(*,*) "calc_ytopo:: Error: floating calving method not recognized."
                write(*,*) "calv_flt_method = ", trim(tpo%par%calv_flt_method)
                stop
    
        end select
    
        ! === Marine terminating calving laws ===

        ! Initialize the calving rates
        tpo%now%cmb_grnd_x = 0.0_wp
        tpo%now%cmb_grnd_y = 0.0_wp
        tpo%now%cmb_grnd   = 0.0_wp

        select case(trim(tpo%par%calv_grnd_method))

            case("zero","none")
                ! Do nothing. No calving.

            case("equil")
                ! For an equilibrated ice sheet calving rates should be opposite to ice velocity
                tpo%now%cmb_grnd_x = -1*dyn%now%ux_bar
                tpo%now%cmb_grnd_y = -1*dyn%now%uy_bar

            case("threshold")
                ! Ice thickness threshold.
                call calc_calving_threshold_lsf(tpo%now%cmb_grnd_x,tpo%now%cmb_grnd_y,dyn%now%ux_bar,dyn%now%uy_bar,tpo%now%H_ice,tpo%par%Hc_ref_grnd,tpo%now%f_ice,tpo%par%boundaries)
        
            case("vm-m16")
                call calc_calving_rate_vonmises_m16(tpo%now%cmb_grnd_x,tpo%now%cmb_grnd_y,dyn%now%ux_bar,dyn%now%uy_bar,mat%now%strs2D%tau_eig_1,tpo%par%tau_ice,tpo%now%f_ice,tpo%par%boundaries)    

            ! Add new laws
            ! MICI should be a marine terminating calving law (only for grounding-line points?)

            case DEFAULT
    
                write(*,*) "calc_ytopo:: Error: grounded calving method not recognized."
                write(*,*) "calv_grnd_method = ", trim(tpo%par%calv_grnd_method)
                stop
    
        end select
        
        ! === Land terminating calving laws ===
        ! For the moment we will assume no calving laws for land-terminating ice points.
        ! Only deformation.
    
        ! === Merge all calving law ===
        ! Merge all calving-rates into a single velocity field.
        ! Using ac-nodes for indices now.
        tpo%now%cr_acx = 0.0_wp
        tpo%now%cr_acy = 0.0_wp
        
        do j=1,ny
        do i=1,nx
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)
            ! x-ac node
            if (tpo%now%f_grnd_acx(i,j) .eq. 0.0) then
                ! Floating point
                tpo%now%cr_acx(i,j) = tpo%now%cmb_flt_x(i,j)
            else
                if(0.5*(bnd%z_bed(i,j)+bnd%z_bed(ip1,j)) .gt. 0.0_wp) then
                    ! Point above sea level. Do not allow to move here. (check)
                    tpo%now%cr_acx(i,j) = -1*dyn%now%ux_bar(i,j)
                else
                    ! Marine-terminating point.
                    tpo%now%cr_acx(i,j) = tpo%now%cmb_grnd_x(i,j)
                end if
            end if
                
            ! y-ac node
            if (tpo%now%f_grnd_acy(i,j) .eq. 0.0) then
                ! Floating point
                tpo%now%cr_acy(i,j) = tpo%now%cmb_flt_y(i,j)
            else
                if(0.5*(bnd%z_bed(i,j)+bnd%z_bed(i,jp1)) .gt. 0.0_wp) then
                    ! Point above sea level. Do nothing for now.
                    tpo%now%cr_acy(i,j) = -1*dyn%now%uy_bar(i,j)
                else
                    ! Marine-terminating point.
                    tpo%now%cr_acy(i,j) = tpo%now%cmb_grnd_y(i,j)
                end if
            end if
    
        end do
        end do

        ! === LSF advection ===
        ! Store previous lsf mask. Necessary to avoid compute it two times.
        tpo%now%lsf_n = tpo%now%lsf
        call LSFupdate(tpo%now%dlsfdt,tpo%now%lsf,tpo%now%cr_acx,tpo%now%cr_acy,dyn%now%ux_bar,dyn%now%uy_bar, &
                       tpo%now%mask_adv,tpo%par%dx,tpo%par%dy,dt,tpo%par%solver)

        ! LSF should not affect points above sea level
        where(bnd%z_bed .gt. 0.0_wp) tpo%now%lsf = -1.0_wp

        ! === Calving ===
        ! Apply calving as a melt rate equal to ice thickness where lsf is positive
        tpo%now%cmb = 0.0_wp
        do j=1,ny
        do i=1,nx
            call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)
            ! Compute the mean calving rate in every aa node
            tpo%now%cmb_flt(i,j) = ((0.5*(tpo%now%cmb_flt_x(im1,j)+tpo%now%cmb_flt_x(i,j)))**2 + &
                                    (0.5*(tpo%now%cmb_flt_y(i,jm1)+tpo%now%cmb_flt_y(i,j)))**2)**0.5

            ! Redefine LSF    
            if(tpo%now%lsf(i,j) .gt. 0.0_wp) then
                ! Calve ice outside LSF mask (cmb = H_ice)
                tpo%now%cmb(i,j) =  -(tpo%now%H_ice(i,j) / dt_kill)

                ! reset LSF border
                if ((tpo%now%lsf(im1,j) .gt. 0.0_wp) .and. (tpo%now%lsf(ip1,j) .gt. 0.0_wp) .and. &
                    (tpo%now%lsf(i,jm1) .gt. 0.0_wp) .and. (tpo%now%lsf(i,jp1) .gt. 0.0_wp)) then
                    tpo%now%lsf(i,j) = 1.0_wp
                end if
            else
                ! reset LSF border
                if ((tpo%now%lsf(im1,j) .le. 0.0_wp) .and. (tpo%now%lsf(ip1,j) .le. 0.0_wp) .and. &
                    (tpo%now%lsf(i,jm1) .le. 0.0_wp) .and. (tpo%now%lsf(i,jp1) .le. 0.0_wp)) then
                    tpo%now%lsf(i,j) = -1.0_wp
                end if
            end if
        end do
        end do

        ! Apply rate and update ice thickness
        call apply_tendency(tpo%now%H_ice,tpo%now%cmb,dt,"calving_lsf",adjust_mb=.TRUE.)

        ! Update ice fraction mask 
        call calc_ice_fraction(tpo%now%f_ice,tpo%now%H_ice,bnd%z_bed,bnd%z_sl,bnd%c%rho_ice, &
            bnd%c%rho_sw,tpo%par%boundaries,tpo%par%margin_flt_subgrid)

        ! Treat fractional points that are not connected to full ice-covered points
        call calc_G_remove_fractional_ice(mbal_now,tpo%now%H_ice,tpo%now%f_ice,dt)

        ! Apply rate and update ice thickness
        mbal_now = 0.0_wp
        call apply_tendency(tpo%now%H_ice,mbal_now,dt,"frac",adjust_mb=.TRUE.)

        ! Add this rate to calving tendency
        tpo%now%cmb = tpo%now%cmb + mbal_now

        call calc_ice_fraction(tpo%now%f_ice,tpo%now%H_ice,bnd%z_bed,bnd%z_sl,bnd%c%rho_ice, &
                        bnd%c%rho_sw,tpo%par%boundaries,tpo%par%margin_flt_subgrid)

        ! reset LSF function after dt_lsf (only if dt_lsf is positive)
        if (tpo%par%dt_lsf .gt. 0.0) then
            if (mod(nint(time_now*100),nint(tpo%par%dt_lsf*100))==0) then
                where(tpo%now%lsf .gt. 0.0) tpo%now%lsf = 1.0
                where(tpo%now%lsf .le. 0.0) tpo%now%lsf = -1.0
            end if
        end if

        ! if there is no ice (for example due to oceanic melt) ensure that point is now ocean in the lsf mask
        select case(trim(tpo%par%calv_flt_method))
            case("equil")
                    ! Do nothing here
            case DEFAULT
                    where(tpo%now%H_ice .le. 0.0 .and. tpo%now%lsf .lt. 0.0 .and. bnd%z_bed .lt. 0.0) tpo%now%lsf = 1.0_wp
        end select 

        return
    
    end subroutine calc_ytopo_calving_lsf

    subroutine calc_ytopo_diagnostic(tpo,dyn,mat,thrm,bnd)
        ! Calculate adjustments to surface elevation, bedrock elevation
        ! and ice thickness 

        implicit none 

        type(ytopo_class),  intent(INOUT) :: tpo
        type(ydyn_class),   intent(IN)    :: dyn
        type(ymat_class),   intent(IN)    :: mat
        type(ytherm_class), intent(IN)    :: thrm  
        type(ybound_class), intent(IN)    :: bnd 

        ! Local variables
        character(len=256) :: bcx, bcy 

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

if (.TRUE.) then
        ! New routines 
        call calc_gradient_acx(tpo%now%dzsdx,tpo%now%z_srf,tpo%now%f_ice,tpo%par%dx,tpo%par%grad_lim,tpo%par%margin2nd,zero_outside=.FALSE.,boundaries=tpo%par%boundaries)
        call calc_gradient_acy(tpo%now%dzsdy,tpo%now%z_srf,tpo%now%f_ice,tpo%par%dy,tpo%par%grad_lim,tpo%par%margin2nd,zero_outside=.FALSE.,boundaries=tpo%par%boundaries)
        
        call calc_gradient_acx(tpo%now%dHidx,tpo%now%H_ice,tpo%now%f_ice,tpo%par%dx,tpo%par%grad_lim,tpo%par%margin2nd,zero_outside=.TRUE.,boundaries=tpo%par%boundaries)
        call calc_gradient_acy(tpo%now%dHidy,tpo%now%H_ice,tpo%now%f_ice,tpo%par%dy,tpo%par%grad_lim,tpo%par%margin2nd,zero_outside=.TRUE.,boundaries=tpo%par%boundaries)
        
        call calc_gradient_acx(tpo%now%dzbdx,tpo%now%z_base,tpo%now%f_ice,tpo%par%dx,tpo%par%grad_lim,tpo%par%margin2nd,zero_outside=.FALSE.,boundaries=tpo%par%boundaries)
        call calc_gradient_acy(tpo%now%dzbdy,tpo%now%z_base,tpo%now%f_ice,tpo%par%dy,tpo%par%grad_lim,tpo%par%margin2nd,zero_outside=.FALSE.,boundaries=tpo%par%boundaries)
else
        bcx = trim(tpo%par%boundaries)
        if (trim(bcx) .eq. "periodic-x") bcx = "periodic"
        bcy = trim(tpo%par%boundaries)
        if (trim(bcy) .eq. "periodic-y") bcy = "periodic"
        
        call calc_dvdx_2D(tpo%now%dzsdx_aa,tpo%now%z_srf,tpo%par%dx,tpo%now%f_ice .gt. 0.0_wp,bcx,tpo%par%grad_lim)
        call calc_dvdy_2D(tpo%now%dzsdy_aa,tpo%now%z_srf,tpo%par%dy,tpo%now%f_ice .gt. 0.0_wp,bcy,tpo%par%grad_lim)
        
        call calc_dvdx_2D(tpo%now%dHidx_aa,tpo%now%H_ice,tpo%par%dx,tpo%now%f_ice .gt. 0.0_wp,bcx,tpo%par%grad_lim)
        call calc_dvdy_2D(tpo%now%dHidy_aa,tpo%now%H_ice,tpo%par%dy,tpo%now%f_ice .gt. 0.0_wp,bcy,tpo%par%grad_lim)
        
        call calc_dvdx_2D(tpo%now%dzbdx_aa,tpo%now%z_base,tpo%par%dx,tpo%now%f_ice .gt. 0.0_wp,bcx,tpo%par%grad_lim)
        call calc_dvdy_2D(tpo%now%dzbdy_aa,tpo%now%z_base,tpo%par%dy,tpo%now%f_ice .gt. 0.0_wp,bcy,tpo%par%grad_lim)
        
        ! Stagger to acx and acy nodes

        tpo%now%dzsdx = stagger_aa_acx(tpo%now%dzsdx_aa)
        tpo%now%dzsdy = stagger_aa_acy(tpo%now%dzsdy_aa)
        
        tpo%now%dHidx = stagger_aa_acx(tpo%now%dHidx_aa)
        tpo%now%dHidy = stagger_aa_acy(tpo%now%dHidy_aa)
        
        tpo%now%dzbdx = stagger_aa_acx(tpo%now%dzbdx_aa)
        tpo%now%dzbdy = stagger_aa_acy(tpo%now%dzbdy_aa)

end if

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

if (.FALSE.) then
        if (tpo%par%dmb_method .gt. 0) then
            ! Calculate the grounding-line distance
            call calc_distance_to_grounding_line(tpo%now%dist_grline,tpo%now%f_grnd,tpo%par%dx, &
                                                        tpo%par%boundaries,calc_distances=.TRUE.)

            ! Define the grounding-zone mask too 
            call calc_grounding_line_zone(tpo%now%mask_grz,tpo%now%dist_grline,tpo%par%dist_grz)

            ! Calculate distance to the ice margin
            call calc_distance_to_ice_margin(tpo%now%dist_margin,tpo%now%f_ice,tpo%par%dx, &
                                                        tpo%par%boundaries,calc_distances=.TRUE.)
        else
            ! Calculate the grounding-line distance
            call calc_distance_to_grounding_line(tpo%now%dist_grline,tpo%now%f_grnd,tpo%par%dx, &
                                                        tpo%par%boundaries,calc_distances=.FALSE.)

            ! Define the grounding-zone mask too 
            call calc_grounding_line_zone(tpo%now%mask_grz,tpo%now%dist_grline,tpo%par%dist_grz)

                ! Calculate distance to the ice margin
            call calc_distance_to_ice_margin(tpo%now%dist_margin,tpo%now%f_ice,tpo%par%dx, &
                                                        tpo%par%boundaries,calc_distances=.FALSE.)

        end if
else
        ! Calculate the grounding-line distance
        call calc_distance_to_grounding_line(tpo%now%dist_grline,tpo%now%f_grnd,tpo%par%dx, &
                                                    tpo%par%boundaries,calc_distances=.FALSE.)

        ! Define the grounding-zone mask too 
        call calc_grounding_line_zone(tpo%now%mask_grz,tpo%now%dist_grline,tpo%par%dist_grz)

        call compute_distance_to_mask(tpo%now%dist_grline, tpo%now%mask_grz)

end if

        ! Calculate the general bed mask
        call gen_mask_bed(tpo%now%mask_bed,tpo%now%f_ice,thrm%now%f_pmp, &
                                            tpo%now%f_grnd,tpo%now%mask_grz.eq.0)


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
                tpo%now%rates%dlsfdt        = 0.0

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
                tpo%now%rates%dlsfdt        = tpo%now%rates%dlsfdt      + tpo%now%dlsfdt*dt


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
                    tpo%now%rates%dlsfdt        = tpo%now%rates%dlsfdt / tpo%now%rates%dt_tot
                    
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
                        tpo%now%dlsfdt      = tpo%now%rates%dlsfdt

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
    
    subroutine ytopo_par_load(par,filename,group_ytopo,group_ycalv,nx,ny,dx,init)

        type(ytopo_param_class), intent(OUT) :: par
        character(len=*),        intent(IN)  :: filename
        character(len=*),        intent(IN)  :: group_ytopo ! Usually "ytopo"
        character(len=*),        intent(IN)  :: group_ycalv ! calving group
        integer,                 intent(IN)  :: nx, ny 
        real(wp),                intent(IN)  :: dx  
        logical, optional,       intent(IN)  :: init 

        ! Local variables
        logical :: init_pars 

        init_pars = .FALSE.
        if (present(init)) init_pars = .TRUE. 
 
        ! Store parameter values in output object
        call nml_read(filename,group_ytopo,"solver",            par%solver,           init=init_pars)
        call nml_read(filename,group_ytopo,"surf_gl_method",    par%surf_gl_method,   init=init_pars)
        call nml_read(filename,group_ytopo,"grad_lim",          par%grad_lim,         init=init_pars)
        call nml_read(filename,group_ytopo,"grad_lim_zb",       par%grad_lim_zb,      init=init_pars)
        call nml_read(filename,group_ytopo,"dHdt_dyn_lim",      par%dHdt_dyn_lim,     init=init_pars) 
        call nml_read(filename,group_ytopo,"margin2nd",         par%margin2nd,        init=init_pars)
        call nml_read(filename,group_ytopo,"margin_flt_subgrid",par%margin_flt_subgrid,init=init_pars)
        call nml_read(filename,group_ytopo,"use_bmb",           par%use_bmb,          init=init_pars)
        call nml_read(filename,group_ytopo,"topo_fixed",        par%topo_fixed,       init=init_pars)
        call nml_read(filename,group_ytopo,"topo_rel",          par%topo_rel,         init=init_pars)
        call nml_read(filename,group_ytopo,"topo_rel_tau",      par%topo_rel_tau,     init=init_pars)
        call nml_read(filename,group_ytopo,"topo_rel_field",    par%topo_rel_field,   init=init_pars)
        ! Grounding line
        call nml_read(filename,group_ytopo,"bmb_gl_method",     par%bmb_gl_method,    init=init_pars)
        call nml_read(filename,group_ytopo,"gl_sep",            par%gl_sep,           init=init_pars)
        call nml_read(filename,group_ytopo,"gz_nx",             par%gz_nx,            init=init_pars)
        ! pmpt method
        call nml_read(filename,group_ytopo,"dist_grz",          par%dist_grz,         init=init_pars)
        call nml_read(filename,group_ytopo,"gz_Hg0",            par%gz_Hg0,           init=init_pars)
        call nml_read(filename,group_ytopo,"gz_Hg1",            par%gz_Hg1,           init=init_pars)
        ! dmb
        call nml_read(filename,group_ytopo,"dmb_method",        par%dmb_method,       init=init_pars)
        call nml_read(filename,group_ytopo,"dmb_alpha_max",     par%dmb_alpha_max,    init=init_pars)
        call nml_read(filename,group_ytopo,"dmb_tau",           par%dmb_tau,          init=init_pars)
        call nml_read(filename,group_ytopo,"dmb_sigma_ref",     par%dmb_sigma_ref,    init=init_pars)
        call nml_read(filename,group_ytopo,"dmb_m_d",           par%dmb_m_d,          init=init_pars)
        call nml_read(filename,group_ytopo,"dmb_m_r",           par%dmb_m_r,          init=init_pars)
        ! fmb
        call nml_read(filename,group_ytopo,"fmb_method",        par%fmb_method,       init=init_pars)
        call nml_read(filename,group_ytopo,"fmb_scale",         par%fmb_scale,        init=init_pars)

        ! === read calving routine ===
        call nml_read(filename,group_ycalv,"use_lsf",           par%use_lsf,            init=init_pars)
        call nml_read(filename,group_ycalv,"dt_lsf",            par%dt_lsf,             init=init_pars)        
        call nml_read(filename,group_ycalv,"calv_flt_method",   par%calv_flt_method,    init=init_pars)
        call nml_read(filename,group_ycalv,"calv_grnd_method",  par%calv_grnd_method,   init=init_pars)
        ! ?
        call nml_read(filename,group_ycalv,"H_min_grnd",        par%H_min_grnd,         init=init_pars)
        call nml_read(filename,group_ycalv,"H_min_flt",         par%H_min_flt,          init=init_pars)
        call nml_read(filename,group_ycalv,"sd_min",            par%sd_min,             init=init_pars)
        call nml_read(filename,group_ycalv,"sd_max",            par%sd_max,             init=init_pars)
        call nml_read(filename,group_ycalv,"calv_grnd_max",     par%calv_grnd_max,      init=init_pars) 
        !
        call nml_read(filename,group_ycalv,"calv_tau",          par%calv_tau,           init=init_pars)
        call nml_read(filename,group_ycalv,"calv_thin",         par%calv_thin,          init=init_pars)
        call nml_read(filename,group_ycalv,"k2",                par%k2,                 init=init_pars)
        call nml_read(filename,group_ycalv,"w2",                par%w2,                 init=init_pars)
        call nml_read(filename,group_ycalv,"kt_ref",            par%kt_ref,             init=init_pars)
        call nml_read(filename,group_ycalv,"kt_deep",           par%kt_deep,            init=init_pars)
        call nml_read(filename,group_ycalv,"tau_ice",           par%tau_ice,            init=init_pars)
        ! Threshold method
        call nml_read(filename,group_ycalv,"Hc_ref_flt",        par%Hc_ref_flt,       init=init_pars)
        call nml_read(filename,group_ycalv,"Hc_ref_grnd",       par%Hc_ref_grnd,      init=init_pars)
        call nml_read(filename,group_ycalv,"Hc_ref_thin",       par%Hc_ref_thin,      init=init_pars)
        call nml_read(filename,group_ycalv,"Hc_deep",           par%Hc_deep,          init=init_pars)
        call nml_read(filename,group_ycalv,"zb_deep_0",         par%zb_deep_0,        init=init_pars)
        call nml_read(filename,group_ycalv,"zb_deep_1",         par%zb_deep_1,        init=init_pars)
        call nml_read(filename,group_ycalv,"zb_sigma",          par%zb_sigma,         init=init_pars)
        
        ! === Set internal parameters ====
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
        allocate(now%rates%dlsfdt(nx,ny))
        
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
        allocate(now%cmb_flt_x(nx,ny))
        allocate(now%cmb_flt_y(nx,ny))
        allocate(now%cmb_grnd(nx,ny))
        allocate(now%cmb_grnd_x(nx,ny))
        allocate(now%cmb_grnd_y(nx,ny))
        allocate(now%cr_acx(nx,ny))
        allocate(now%cr_acy(nx,ny))

        allocate(now%lsf(nx,ny))       
        allocate(now%dlsfdt(nx,ny))
        
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

        allocate(now%dzsdx_aa(nx,ny))
        allocate(now%dzsdy_aa(nx,ny))
        allocate(now%dHidx_aa(nx,ny))
        allocate(now%dHidy_aa(nx,ny))
        allocate(now%dzbdx_aa(nx,ny))
        allocate(now%dzbdy_aa(nx,ny))

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
        allocate(now%lsf_n(nx,ny))
        
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
        now%rates%dlsfdt        = 0.0

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
        now%cmb_flt_x   = 0.0
        now%cmb_flt_y   = 0.0
        now%cmb_grnd    = 0.0
        now%lsf         = 1.0 ! init to 0.0?       
        now%dlsfdt      = 0.0
        
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

        now%dzsdx_aa    = 0.0 
        now%dzsdy_aa    = 0.0 
        now%dHidx_aa    = 0.0 
        now%dHidy_aa    = 0.0
        now%dzbdx_aa    = 0.0 
        now%dzbdy_aa    = 0.0

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
        now%lsf_n     = 0.0 

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
        if (allocated(now%rates%dlsfdt))        deallocate(now%rates%dlsfdt)
        
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
        if (allocated(now%cmb_flt_x))   deallocate(now%cmb_flt_x)
        if (allocated(now%cmb_flt_y))   deallocate(now%cmb_flt_y)
        if (allocated(now%cmb_grnd))    deallocate(now%cmb_grnd)
        if (allocated(now%cmb_grnd_x))  deallocate(now%cmb_grnd_x)
        if (allocated(now%cmb_grnd_y))  deallocate(now%cmb_grnd_y)
        if (allocated(now%cr_acx))      deallocate(now%cr_acx)
        if (allocated(now%cr_acy))      deallocate(now%cr_acy)
        if (allocated(now%lsf))         deallocate(now%lsf)       
        if (allocated(now%dlsfdt))      deallocate(now%dlsfdt)
        
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
        
        if (allocated(now%dzsdx_aa))       deallocate(now%dzsdx_aa)
        if (allocated(now%dzsdy_aa))       deallocate(now%dzsdy_aa)
        if (allocated(now%dHidx_aa))       deallocate(now%dHidx_aa)
        if (allocated(now%dHidy_aa))       deallocate(now%dHidy_aa)
        if (allocated(now%dzbdx_aa))       deallocate(now%dzbdx_aa)
        if (allocated(now%dzbdy_aa))       deallocate(now%dzbdy_aa)
        
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
        if (allocated(now%lsf_n))       deallocate(now%lsf_n)
        
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
        allocate(pc%lsf(nx,ny))
        
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
        pc%lsf          = 0.0            
        
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
        if (allocated(pc%lsf))          deallocate(pc%lsf)
        
        return

    end subroutine ytopo_pc_dealloc

end module yelmo_topography
