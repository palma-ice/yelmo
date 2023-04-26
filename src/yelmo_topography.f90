
module yelmo_topography

    use nml  
    use ncio
    
    use yelmo_defs
    use yelmo_tools 
    
    use mass_conservation
    use calving
    use topography 

    use runge_kutta 

    use grid_calcs 

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
    public :: calc_ytopo_rk4
    public :: calc_ytopo_pc 
    public :: calc_ytopo_diagnostic 
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
    
    subroutine calc_ytopo_rk4(tpo,dyn,mat,thrm,bnd,time,topo_fixed)

        implicit none 

        type(ytopo_class),  intent(INOUT) :: tpo
        type(ydyn_class),   intent(IN)    :: dyn
        type(ymat_class),   intent(IN)    :: mat
        type(ytherm_class), intent(IN)    :: thrm  
        type(ybound_class), intent(IN)    :: bnd 
        real(wp),           intent(IN)    :: time
        logical,            intent(IN)    :: topo_fixed  

        ! Local variables 
        integer  :: i, j, nx, ny
        real(wp) :: dt  
        real(wp), allocatable :: mbal(:,:) 
        real(wp), allocatable :: dHdt_now(:,:) 
        real(wp), allocatable :: G_mb(:,:)  
        
        nx = size(tpo%now%H_ice,1)
        ny = size(tpo%now%H_ice,2)

        allocate(mbal(nx,ny))
        allocate(dHdt_now(nx,ny))
        allocate(G_mb(nx,ny))

        ! Initialize time if necessary 
        if (tpo%par%time .gt. dble(time)) then 
            tpo%par%time = dble(time) 
        end if 
        
        ! Initialize rk4 object if necessary 
        if (.not. allocated(tpo%rk4%tau)) then 
            call rk4_2D_init(tpo%rk4,nx,ny,dt_min=1e-3)
        end if 
        
        ! Get time step
        dt = dble(time) - tpo%par%time 

        ! Step 0: Get some diagnostic quantities for mass balance calculation --------

        ! Calculate grounded fraction on aa-nodes
        ! (only to be used with basal mass balance, later all
        !  f_grnd arrays will be calculated according to use choices)
        call determine_grounded_fractions(tpo%now%f_grnd_bmb,H_grnd=tpo%now%H_grnd)
        
        ! Combine basal mass balance into one field accounting for 
        ! grounded/floating fraction of grid cells 
        call calc_bmb_total(tpo%now%bmb,thrm%now%bmb_grnd,bnd%bmb_shlf,tpo%now%H_ice,tpo%now%H_grnd, &
                            tpo%now%f_grnd_bmb,tpo%par%bmb_gl_method,tpo%par%diffuse_bmb_shlf)
        
        ! Combine frontal mass balance into one field, and 
        ! calculate as needed 
        call calc_fmb_total(tpo%now%fmb,bnd%fmb_shlf,bnd%bmb_shlf,tpo%now%H_ice, &
                        tpo%now%H_grnd,tpo%now%f_ice,tpo%par%fmb_method,tpo%par%fmb_scale, &
                        bnd%c%rho_ice,bnd%c%rho_sw,tpo%par%dx)


        ! Define temporary variable for total column mass balance (without calving)
        mbal = bnd%smb + tpo%now%bmb + tpo%now%fmb
        
        ! WHEN RUNNING EISMINT1 ensure bmb and fmb are not accounted for here !!!
        if (.not. tpo%par%use_bmb) then
            mbal = bnd%smb
        end if


        ! Step 1: Go through predictor-corrector-advance steps

        if ( .not. topo_fixed .and. dt .gt. 0.0 ) then 

            ! Store ice thickness and surface elevation from previous timestep
            tpo%now%H_ice_n = tpo%now%H_ice 
            tpo%now%z_srf_n = tpo%now%z_srf 

            
            ! Get ice-fraction mask for current ice thickness  
            call calc_ice_fraction(tpo%now%f_ice,tpo%now%H_ice,bnd%z_bed,bnd%z_sl, &
                                    bnd%c%rho_ice,bnd%c%rho_sw,tpo%par%margin_flt_subgrid)

            ! Calculate ice thickness evolution from dynamics alone
            ! call rk23_2D_step(tpo%rk4,tpo%now%H_ice,tpo%now%f_ice,dHdt_now,dyn%now%ux_bar,dyn%now%uy_bar, &
            !                                     tpo%par%dx,dt,tpo%par%solver,tpo%par%boundaries)
            call rk4_2D_step(tpo%rk4,tpo%now%H_ice,tpo%now%f_ice,dHdt_now,dyn%now%ux_bar,dyn%now%uy_bar, &
                                                tpo%now%mask_adv,tpo%par%dx,dt,tpo%par%solver,tpo%par%boundaries)

            ! Calculate rate of change using weighted advective rates of change 
            ! depending on timestepping method chosen 
            tpo%now%dHdt_pred  = dHdt_now
            tpo%now%dHdt_corr  = dHdt_now

            ! Calculate predicted ice thickness
            !tpo%now%H_ice = tpo%now%H_ice_n + dt*tpo%now%dHdt_pred

            ! Diagnose actual mass balance (forcing) tendency
            call calc_G_mbal(G_mb,tpo%now%H_ice,tpo%now%f_grnd,mbal,dt)

            ! Store for output too 
            tpo%now%mb_applied = G_mb 

            ! Now update ice thickness with all tendencies for this timestep 
            tpo%now%H_ice = tpo%now%H_ice + dt*G_mb  

            ! Ensure tiny numeric ice thicknesses are removed
            where (tpo%now%H_ice .lt. TOL_UNDERFLOW) tpo%now%H_ice = 0.0 

            ! Calculate and apply calving
            call calc_ytopo_calving(tpo,dyn,mat,thrm,bnd,mbal,dt)

            ! If desired, finally relax solution to reference state
            if (tpo%par%topo_rel .ne. 0) then 

                select case(trim(tpo%par%topo_rel_field))

                    case("H_ref")
                        ! Relax towards reference ice thickness field H_ref

                        call relax_ice_thickness(tpo%now%H_ice,tpo%now%f_grnd,tpo%now%mask_grz,bnd%H_ice_ref, &
                                                    tpo%par%topo_rel,tpo%par%topo_rel_tau,dt,tpo%par%boundaries)
                    
                    case("H_ice_n")
                        ! Relax towards previous iteration ice thickness 
                        ! (ie slow down changes)
                        ! ajr: needs testing, not sure if this works well or helps anything.

                        call relax_ice_thickness(tpo%now%H_ice,tpo%now%f_grnd,tpo%now%mask_grz,tpo%now%H_ice_n, &
                                                    tpo%par%topo_rel,tpo%par%topo_rel_tau,dt,tpo%par%boundaries)
                    
                    case DEFAULT 

                        write(*,*) "calc_ytopo:: Error: topo_rel_field not recognized."
                        write(*,*) "topo_rel_field = ", trim(tpo%par%topo_rel_field)
                        stop 

                end select

            end if

            ! Get ice-fraction mask for ice thickness  
            call calc_ice_fraction(tpo%now%f_ice,tpo%now%H_ice,bnd%z_bed,bnd%z_sl, &
                                    bnd%c%rho_ice,bnd%c%rho_sw,tpo%par%margin_flt_subgrid)
            
            ! Finally apply all additional (generally artificial) ice thickness adjustments 
            ! and store changes in residual mass balance field. 
            ! call apply_ice_thickness_boundaries(tpo%now%mb_resid,tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd, &
            !                                     dyn%now%uxy_b,bnd%ice_allowed,tpo%par%boundaries,bnd%H_ice_ref, &
            !                                     tpo%par%H_min_flt,tpo%par%H_min_grnd,dt,reset=.TRUE.)
            
            tpo%now%H_ice_pred = tpo%now%H_ice 
            tpo%now%H_ice_corr = tpo%now%H_ice 
            
        end if 

        ! Update fields and masks
        call calc_ytopo_diagnostic(tpo,dyn,mat,thrm,bnd)

        ! Determine rates of change
        if ( .not. topo_fixed .and. dt .gt. 0.0 ) then 

            tpo%now%dHidt = (tpo%now%H_ice - tpo%now%H_ice_n) / dt 
            tpo%now%dzsdt = (tpo%now%z_srf - tpo%now%z_srf_n) / dt 
        
        end if 

        ! Update ytopo time to current time 
        tpo%par%time = dble(time)

        return

    end subroutine calc_ytopo_rk4

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
        real(wp), allocatable :: mbal(:,:) 
        real(wp), allocatable :: dHdt_now(:,:) 
        real(wp), allocatable :: G_mb(:,:)  
        
        nx = size(tpo%now%H_ice,1)
        ny = size(tpo%now%H_ice,2)

        allocate(mbal(nx,ny))
        allocate(dHdt_now(nx,ny))
        allocate(G_mb(nx,ny))

        ! Initialize time if necessary 
        if (tpo%par%time .gt. dble(time)) then 
            tpo%par%time = dble(time) 
        end if 
        
        ! Initialize rk4 object if necessary 
        if (.not. allocated(tpo%rk4%tau)) then 
            call rk4_2D_init(tpo%rk4,nx,ny,dt_min=1e-3)
        end if 

        ! Get time step
        dt = dble(time) - tpo%par%time 

        ! Step 0: Get some diagnostic quantities for mass balance calculation --------

        ! Calculate grounded fraction on aa-nodes
        ! (only to be used with basal mass balance, later all
        !  f_grnd arrays will be calculated according to use choices)
        call determine_grounded_fractions(tpo%now%f_grnd_bmb,H_grnd=tpo%now%H_grnd)
        
        ! Combine basal mass balance into one field accounting for 
        ! grounded/floating fraction of grid cells 
        call calc_bmb_total(tpo%now%bmb,thrm%now%bmb_grnd,bnd%bmb_shlf,tpo%now%H_ice,tpo%now%H_grnd, &
                            tpo%now%f_grnd_bmb,tpo%par%bmb_gl_method,tpo%par%diffuse_bmb_shlf)
        
        ! Combine frontal mass balance into one field, and 
        ! calculate as needed 
        call calc_fmb_total(tpo%now%fmb,bnd%fmb_shlf,bnd%bmb_shlf,tpo%now%H_ice, &
                        tpo%now%H_grnd,tpo%now%f_ice,tpo%par%fmb_method,tpo%par%fmb_scale, &
                        bnd%c%rho_ice,bnd%c%rho_sw,tpo%par%dx)


        ! Define temporary variable for total column mass balance (without calving)
        mbal = bnd%smb + tpo%now%bmb + tpo%now%fmb
        
        ! WHEN RUNNING EISMINT1 ensure bmb and fmb are not accounted for here !!!
        if (.not. tpo%par%use_bmb) then
            mbal = bnd%smb
        end if


        ! Step 1: Go through predictor-corrector-advance steps

        if ( .not. topo_fixed .and. dt .gt. 0.0 ) then 

            ! Ice thickness evolution from dynamics alone

            select case(trim(pc_step))

                case("predictor") 
                    ! Determine predicted ice thickness 

                    ! Store ice thickness from previous timestep 
                    tpo%now%H_ice_n = tpo%now%H_ice 

                    ! Store previous surface elevation too (for calculating rate of change)
                    tpo%now%z_srf_n = tpo%now%z_srf 

                    ! Get ice-fraction mask for current ice thickness  
                    call calc_ice_fraction(tpo%now%f_ice,tpo%now%H_ice,bnd%z_bed,bnd%z_sl, &
                                            bnd%c%rho_ice,bnd%c%rho_sw,tpo%par%margin_flt_subgrid)
    
if (.TRUE.) then
                    call calc_G_advec_simple(dHdt_now,tpo%now%H_ice,tpo%now%f_ice,dyn%now%ux_bar,dyn%now%uy_bar, &
                                                 tpo%now%mask_adv,tpo%par%solver,tpo%par%boundaries,tpo%par%dx,dt)

else
                    
                    call rk4_2D_step(tpo%rk4,tpo%now%H_ice,tpo%now%f_ice,dHdt_now,dyn%now%ux_bar,dyn%now%uy_bar, &
                                                tpo%now%mask_adv,tpo%par%dx,dt,tpo%par%solver,tpo%par%boundaries)

end if 
                    
                    ! Calculate rate of change using weighted advective rates of change 
                    ! depending on timestepping method chosen 
                    tpo%now%dHdt_pred = tpo%par%dt_beta(1)*dHdt_now + tpo%par%dt_beta(2)*tpo%now%dHdt_n 

                    ! Calculate predicted ice thickness
                    tpo%now%H_ice = tpo%now%H_ice_n + dt*tpo%now%dHdt_pred

                case("corrector") 

                    ! Set current thickness to predicted thickness
                    tpo%now%H_ice = tpo%now%H_ice_pred 

                    ! Get ice-fraction mask for predicted ice thickness  
                    call calc_ice_fraction(tpo%now%f_ice,tpo%now%H_ice,bnd%z_bed,bnd%z_sl, &
                                            bnd%c%rho_ice,bnd%c%rho_sw,tpo%par%margin_flt_subgrid)

if (.TRUE.) then
                    call calc_G_advec_simple(dHdt_now,tpo%now%H_ice,tpo%now%f_ice,dyn%now%ux_bar,dyn%now%uy_bar, &
                                                tpo%now%mask_adv,tpo%par%solver,tpo%par%boundaries,tpo%par%dx,dt)

else

                    call rk4_2D_step(tpo%rk4,tpo%now%H_ice,tpo%now%f_ice,dHdt_now,dyn%now%ux_bar,dyn%now%uy_bar, &
                                                tpo%now%mask_adv,tpo%par%dx,dt,tpo%par%solver,tpo%par%boundaries)

end if

                    ! Calculate rate of change using weighted advective rates of change 
                    ! depending on timestepping method chosen 
                    tpo%now%dHdt_corr = tpo%par%dt_beta(3)*dHdt_now + tpo%par%dt_beta(4)*tpo%now%dHdt_n 
                    
                    ! Calculate corrected ice thickness
                    tpo%now%H_ice = tpo%now%H_ice_n + dt*tpo%now%dHdt_corr
                    
                    ! Store dHdt_corr as dHdt_n now for use with next timestep 
                    tpo%now%dHdt_n = tpo%now%dHdt_corr 

            end select

            ! Diagnose actual mass balance (forcing) tendency
            call calc_G_mbal(G_mb,tpo%now%H_ice_n,tpo%now%f_grnd,mbal,dt)

            ! Store for output too 
            tpo%now%mb_applied = G_mb 

            ! Now update ice thickness with all tendencies for this timestep 
            tpo%now%H_ice = tpo%now%H_ice + dt*G_mb  

            ! Ensure tiny numeric ice thicknesses are removed
            where (tpo%now%H_ice .lt. TOL_UNDERFLOW) tpo%now%H_ice = 0.0 

            ! Calculate and apply calving
            call calc_ytopo_calving(tpo,dyn,mat,thrm,bnd,mbal,dt)

            ! Get ice-fraction mask for ice thickness  
            call calc_ice_fraction(tpo%now%f_ice,tpo%now%H_ice,bnd%z_bed,bnd%z_sl, &
                                    bnd%c%rho_ice,bnd%c%rho_sw,tpo%par%margin_flt_subgrid)
            
            ! Finally apply all additional (generally artificial) ice thickness adjustments 
            ! and store changes in residual mass balance field. 
            call apply_ice_thickness_boundaries(tpo%now%mb_resid,tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd, &
                                                dyn%now%uxy_b,bnd%ice_allowed,tpo%par%boundaries,bnd%H_ice_ref, &
                                                tpo%par%H_min_flt,tpo%par%H_min_grnd,dt,reset=.TRUE.)

            ! If desired, finally relax solution to reference state
            if (tpo%par%topo_rel .ne. 0) then 

                select case(trim(tpo%par%topo_rel_field))

                    case("H_ref")
                        ! Relax towards reference ice thickness field H_ref

                        call relax_ice_thickness(tpo%now%H_ice,tpo%now%f_grnd,tpo%now%mask_grz,bnd%H_ice_ref, &
                                                        tpo%par%topo_rel,tpo%par%topo_rel_tau,dt,tpo%par%boundaries)
                    
                    case("H_ice_n")
                        ! Relax towards previous iteration ice thickness 
                        ! (ie slow down changes)
                        ! ajr: needs testing, not sure if this works well or helps anything.

                        call relax_ice_thickness(tpo%now%H_ice,tpo%now%f_grnd,tpo%now%mask_grz,tpo%now%H_ice_n, &
                                                        tpo%par%topo_rel,tpo%par%topo_rel_tau,dt,tpo%par%boundaries)
                    
                    case DEFAULT 

                        write(*,*) "calc_ytopo:: Error: topo_rel_field not recognized."
                        write(*,*) "topo_rel_field = ", trim(tpo%par%topo_rel_field)
                        stop 

                end select

                ! Get ice-fraction mask for ice thickness  
                call calc_ice_fraction(tpo%now%f_ice,tpo%now%H_ice,bnd%z_bed,bnd%z_sl, &
                                        bnd%c%rho_ice,bnd%c%rho_sw,tpo%par%margin_flt_subgrid)
            
                ! Again apply ice thickness boundaries to ensure relaxed fields are consistent 
                ! with desired limitations. 
                call apply_ice_thickness_boundaries(tpo%now%mb_resid,tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd, &
                                                    dyn%now%uxy_b,bnd%ice_allowed,tpo%par%boundaries,bnd%H_ice_ref, &
                                                    tpo%par%H_min_flt,tpo%par%H_min_grnd,dt,reset=.FALSE.)

            end if

            ! Get ice-fraction mask for ice thickness  
            call calc_ice_fraction(tpo%now%f_ice,tpo%now%H_ice,bnd%z_bed,bnd%z_sl, &
                                    bnd%c%rho_ice,bnd%c%rho_sw,tpo%par%margin_flt_subgrid)
            

            select case(trim(pc_step))

                case("predictor") 

                    ! Save predictor field, proceed with predictor field
                    ! in the main H_ice variable for calculating dynamics.
                    tpo%now%H_ice_pred = tpo%now%H_ice 

                    ! Compare previous and current ice field
                    tpo%now%mask_pred_new = 0 
                    where(tpo%now%H_ice_pred .gt. 0.0) tpo%now%mask_pred_new = 1
                    where(tpo%now%H_ice_pred .gt. 0.0 .and. tpo%now%H_ice_n .eq. 0.0)
                        tpo%now%mask_pred_new = 2
                    end where
            

                case("corrector")
                    ! Determine corrected ice thickness 

                    ! Save corrector field too
                    tpo%now%H_ice_corr = tpo%now%H_ice 

                    ! Compare previous and current ice field
                    tpo%now%mask_corr_new = 0 
                    where(tpo%now%H_ice_corr .gt. 0.0) tpo%now%mask_corr_new = 1
                    where(tpo%now%H_ice_corr .gt. 0.0 .and. tpo%now%H_ice_n .eq. 0.0)
                        tpo%now%mask_corr_new = 2
                    end where
                    
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
                        tpo%now%H_ice = tpo%now%H_ice_pred
                    else
                        tpo%now%H_ice = tpo%now%H_ice_corr
                    end if 

                    ! Compare previous and current ice field
                    tpo%now%mask_new = 0 
                    where(tpo%now%H_ice .gt. 0.0) tpo%now%mask_new = 1
                    where(tpo%now%H_ice .gt. 0.0 .and. tpo%now%H_ice_n .eq. 0.0)
                        tpo%now%mask_new = 2
                    end where
                    
            end select

        end if 

        ! Update fields and masks
        call calc_ytopo_diagnostic(tpo,dyn,mat,thrm,bnd)

        ! Determine rates of change
        if ( .not. topo_fixed .and. dt .gt. 0.0 ) then 

            tpo%now%dHidt = (tpo%now%H_ice - tpo%now%H_ice_n) / dt 
            tpo%now%dzsdt = (tpo%now%z_srf - tpo%now%z_srf_n) / dt 
        
        end if 

        if (trim(pc_step) .eq. "advance") then 
            ! Advance timestep here whether topo_fixed was true or not...
            
            ! Update ytopo time to current time 
            tpo%par%time = dble(time)
            
        end if

        return

    end subroutine calc_ytopo_pc

    subroutine calc_ytopo_calving(tpo,dyn,mat,thrm,bnd,mbal,dt)

        implicit none 

        type(ytopo_class),  intent(INOUT) :: tpo
        type(ydyn_class),   intent(IN)    :: dyn
        type(ymat_class),   intent(IN)    :: mat
        type(ytherm_class), intent(IN)    :: thrm  
        type(ybound_class), intent(IN)    :: bnd 
        real(wp),           intent(IN)    :: mbal(:,:) 
        real(wp),           intent(IN)    :: dt

        ! Local variables 
        integer :: nx, ny 
        real(wp), allocatable :: calv_sd(:,:) 

        nx = size(tpo%now%H_ice,1) 
        ny = size(tpo%now%H_ice,2) 

        allocate(calv_sd(nx,ny)) 


        ! Make sure current ice mask is correct
        call calc_ice_fraction(tpo%now%f_ice,tpo%now%H_ice,bnd%z_bed,bnd%z_sl, &
                            bnd%c%rho_ice,bnd%c%rho_sw,tpo%par%margin_flt_subgrid)
        
        ! === CALVING ===

        ! Diagnose strains and stresses relevant to calving 

        ! eps_eff = effective strain = eigencalving e+*e- following Levermann et al. (2012)
        call calc_eps_eff(tpo%now%eps_eff,dyn%now%strn2D%eps_eig_1,dyn%now%strn2D%eps_eig_2,tpo%now%f_ice,tpo%par%boundaries)

        ! tau_eff = effective stress ~ von Mises stress following Lipscomb et al. (2019)
        call calc_tau_eff(tpo%now%tau_eff,mat%now%strs2D%tau_eig_1,mat%now%strs2D%tau_eig_2,tpo%now%f_ice,tpo%par%w2,tpo%par%boundaries)

        select case(trim(tpo%par%calv_flt_method))

            case("zero","none")

                tpo%now%calv_flt = 0.0 

            case("simple") 
                ! Use simple threshold method

                call calc_calving_rate_simple(tpo%now%calv_flt,tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd, &
                                                tpo%par%calv_H_lim,tpo%par%calv_tau,tpo%par%boundaries)
                
            case("flux") 
                ! Use threshold+flux method from Peyaud et al. (2007), ie, GRISLI,
                ! but reformulated. 

                call calc_calving_rate_flux(tpo%now%calv_flt,tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd,mbal, &
                                dyn%now%ux_bar,dyn%now%uy_bar,tpo%par%dx,tpo%par%calv_H_lim,tpo%par%calv_tau,tpo%par%boundaries)
                
            case("flux-grisli")
                ! Use threshold+flux method from Peyaud et al. (2007), ie, GRISLI

                call calc_calving_rate_flux_grisli(tpo%now%calv_flt,tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd,mbal, &
                                dyn%now%ux_bar,dyn%now%uy_bar,tpo%par%dx,tpo%par%calv_H_lim,tpo%par%calv_tau,tpo%par%boundaries)
                
            case("vm-l19")
                ! Use von Mises calving as defined by Lipscomb et al. (2019)

                ! Next, diagnose calving
                call calc_calving_rate_vonmises_l19(tpo%now%calv_flt,tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd, &
                                                                        tpo%now%tau_eff,tpo%par%dx,tpo%par%kt,tpo%par%boundaries)
            case("eigen")
                ! Use Eigen calving as defined by Levermann et al. (2012)

                ! Next, diagnose calving
                call calc_calving_rate_eigen(tpo%now%calv_flt,tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd, &
                                                                        tpo%now%eps_eff,tpo%par%dx,tpo%par%k2,tpo%par%boundaries)

            case("kill") 
                ! Delete all floating ice (using characteristic time parameter)
                ! Make sure dt is a postive number

                call calc_calving_rate_kill(tpo%now%calv_flt,tpo%now%H_ice,tpo%now%f_grnd.eq.0.0_wp, &
                                                                                    tpo%par%calv_tau,dt)

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
                call apply_thin_calving_rate(tpo%now%calv_flt,tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd,tpo%par%calv_thin,tpo%par%boundaries)

                ! Adjust calving so that any excess is distributed to upstream neighbors
                !call calc_calving_residual(tpo%now%calv_flt,tpo%now%H_ice,tpo%now%f_ice,dt)
        
        end select 

        ! Additionally ensure higher calving rate for floating tongues of
        ! one grid-point width.
        call calc_calving_tongues(tpo%now%calv_flt,tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd,tpo%par%calv_tau,tpo%par%boundaries)
        
        ! Diagnose potential grounded-ice calving rate [m/yr]

        select case(trim(tpo%par%calv_grnd_method))

            case("zero","none")

                tpo%now%calv_grnd = 0.0 

            case("stress-b12") 
                ! Use simple threshold method

                call calc_calving_ground_rate_stress_b12(tpo%now%calv_grnd,tpo%now%H_ice,tpo%now%f_ice, &
                                                tpo%now%f_grnd,bnd%z_bed,bnd%z_sl-bnd%z_bed,tpo%par%calv_tau, &
                                                bnd%c%rho_ice,bnd%c%rho_sw,bnd%c%g,tpo%par%boundaries)

            case DEFAULT 

                write(*,*) "calc_ytopo:: Error: grounded calving method not recognized."
                write(*,*) "calv_grnd_method = ", trim(tpo%par%calv_grnd_method)
                stop 

        end select
        
        ! Additionally include parameterized grounded calving 
        ! to account for grid resolution 
        call calc_calving_ground_rate_stdev(calv_sd,tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd, &
                                bnd%z_bed_sd,tpo%par%sd_min,tpo%par%sd_max,tpo%par%calv_max,tpo%par%calv_tau,tpo%par%boundaries)
        tpo%now%calv_grnd = tpo%now%calv_grnd + calv_sd 


        ! Diagnose actual calving (forcing) tendency
        call calc_G_calv(tpo%now%calv,tpo%now%H_ice,tpo%now%calv_flt,tpo%now%calv_grnd,dt, &
                                                                trim(tpo%par%calv_flt_method),tpo%par%boundaries)

        ! Now update ice thickness with all tendencies for this timestep 
        tpo%now%H_ice = tpo%now%H_ice - dt*tpo%now%calv  

        ! Ensure tiny numeric ice thicknesses are removed
        where (tpo%now%H_ice .lt. TOL_UNDERFLOW) tpo%now%H_ice = 0.0 

        ! Update ice fraction mask 
        call calc_ice_fraction(tpo%now%f_ice,tpo%now%H_ice,bnd%z_bed,bnd%z_sl, &
                                bnd%c%rho_ice,bnd%c%rho_sw,tpo%par%margin_flt_subgrid)
    
        ! Treat fractional points that are not connected to full ice-covered points
        call remove_fractional_ice(tpo%now%H_ice,tpo%now%f_ice)
        call calc_ice_fraction(tpo%now%f_ice,tpo%now%H_ice,bnd%z_bed,bnd%z_sl, &
                                bnd%c%rho_ice,bnd%c%rho_sw,tpo%par%margin_flt_subgrid)
        
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
        call calc_ice_fraction(tpo%now%f_ice,tpo%now%H_ice,bnd%z_bed,bnd%z_sl, &
                                bnd%c%rho_ice,bnd%c%rho_sw,tpo%par%margin_flt_subgrid)
        
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
                                                                            tpo%now%H_grnd,tpo%par%gl_sep_nx)
            
            case(3) 
                ! Grounded area using analytical solutions of Leguy et al. (2021)

                call determine_grounded_fractions(tpo%now%f_grnd,tpo%now%f_grnd_acx,tpo%now%f_grnd_acy, &
                                                                            tpo%now%f_grnd_ab,tpo%now%H_grnd)

        end select
        
        ! Calculate grounded fraction due to pinning points 
        call calc_f_grnd_pinning_points(tpo%now%f_grnd_pin,tpo%now%H_ice,tpo%now%f_ice, &
                                                bnd%z_bed,bnd%z_bed_sd,bnd%z_sl,bnd%c%rho_ice,bnd%c%rho_sw)

        ! Calculate the grounding-line distance
        call calc_distance_to_grounding_line(tpo%now%dist_grline,tpo%now%f_grnd,tpo%par%dx)

        ! Define the grounding-zone mask too 
        call calc_grounding_line_zone(tpo%now%mask_grz,tpo%now%dist_grline,tpo%par%dist_grz)

        ! ajr: do not calculate distance to margin unless it is needed (costs some computational time)
        ! Calculate distance to the ice margin
        !call calc_distance_to_ice_margin(tpo%now%dist_margin,tpo%now%f_ice,tpo%par%dx)

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
                call calc_ice_fraction(tpo%now%f_ice_dyn,tpo%now%H_ice_dyn,bnd%z_bed,bnd%z_sl, &
                                                    bnd%c%rho_ice,bnd%c%rho_sw,flt_subgrid=.FALSE.)
            
            case("slab-ext")
                ! Calculate extended ice thickness fields n_ext points
                ! away from grounded margin. 

                tpo%now%H_ice_dyn = tpo%now%H_ice
                where(tpo%now%H_ice_dyn .gt. 0.0 .and. tpo%now%H_ice_dyn .lt. 1.0) &
                        tpo%now%H_ice_dyn = 1.0_wp 

                call extend_floating_slab(tpo%now%H_ice_dyn,tpo%now%f_grnd,H_slab=1.0_wp,n_ext=4)

                ! Calculate the ice fraction mask for use with the dynamics solver
                call calc_ice_fraction(tpo%now%f_ice_dyn,tpo%now%H_ice_dyn,bnd%z_bed,bnd%z_sl, &
                                                    bnd%c%rho_ice,bnd%c%rho_sw,flt_subgrid=.FALSE.)
                
            case DEFAULT 
                ! No modification of ice thickness for dynamics solver 
                ! Set standard ice thickness field for use with dynamics 
            
                tpo%now%H_ice_dyn = tpo%now%H_ice
                tpo%now%f_ice_dyn = tpo%now%f_ice 

        end select

        return 

    end subroutine calc_ytopo_diagnostic

    elemental subroutine gen_mask_bed(mask,f_ice,f_pmp,f_grnd,mask_grz)
        ! Generate an output mask for model conditions at bed
        ! based on input masks 
        ! 0: ocean, 1: land, 2: sia, 3: streams, grline: 4, floating: 5, islands: 6
        ! 7: partially-covered ice cell.

        implicit none 

        integer,    intent(OUT) :: mask 
        real(wp), intent(IN)  :: f_ice, f_pmp, f_grnd
        integer,    intent(IN)  :: mask_grz

        if (mask_grz .eq. 0) then
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
        real(wp),               intent(IN)  :: dx  
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
        call nml_read(filename,"ytopo","grad_lim_zb",       par%grad_lim_zb,      init=init_pars)
        call nml_read(filename,"ytopo","dist_grz",          par%dist_grz,         init=init_pars)
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
        allocate(now%z_base(nx,ny))
        allocate(now%dzsdt(nx,ny))
        allocate(now%dHidt(nx,ny))
        allocate(now%bmb(nx,ny))
        allocate(now%fmb(nx,ny))
        allocate(now%mb_applied(nx,ny))
        allocate(now%mb_resid(nx,ny))
        
        allocate(now%G_advec(nx,ny))
        allocate(now%mask_adv(nx,ny))
        
        allocate(now%mask_new(nx,ny))
        allocate(now%mask_pred_new(nx,ny))
        allocate(now%mask_corr_new(nx,ny))
        
        allocate(now%eps_eff(nx,ny))
        allocate(now%tau_eff(nx,ny))
        allocate(now%calv(nx,ny))
        allocate(now%calv_flt(nx,ny))
        allocate(now%calv_grnd(nx,ny))
        
        allocate(now%dzsdx(nx,ny))
        allocate(now%dzsdy(nx,ny))
        allocate(now%dHidx(nx,ny))
        allocate(now%dHidy(nx,ny))
        allocate(now%dzbdx(nx,ny))
        allocate(now%dzbdy(nx,ny))

        allocate(now%H_eff(nx,ny))
        allocate(now%H_grnd(nx,ny))

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
        
        allocate(now%dHdt_n(nx,ny))
        allocate(now%dHdt_pred(nx,ny))
        allocate(now%dHdt_corr(nx,ny))
        allocate(now%H_ice_n(nx,ny))
        allocate(now%H_ice_pred(nx,ny))
        allocate(now%H_ice_corr(nx,ny))
        
        allocate(now%z_srf_n(nx,ny))
        
        allocate(now%H_ice_dyn(nx,ny))
        allocate(now%f_ice_dyn(nx,ny))

        now%H_ice       = 0.0 
        now%z_srf       = 0.0
        now%z_base      = 0.0  
        now%dzsdt       = 0.0 
        now%dHidt       = 0.0
        now%bmb         = 0.0  
        now%fmb         = 0.0
        now%mb_applied  = 0.0 
        now%mb_resid    = 0.0

        now%G_advec     = 0.0
        now%mask_adv    = 0

        now%mask_new      = 0 
        now%mask_pred_new = 0 
        now%mask_corr_new = 0 

        now%eps_eff     = 0.0
        now%tau_eff     = 0.0
        now%calv        = 0.0
        now%calv_flt    = 0.0
        now%calv_grnd   = 0.0
        now%dzsdx       = 0.0 
        now%dzsdy       = 0.0 
        now%dHidx       = 0.0 
        now%dHidy       = 0.0
        now%dzbdx       = 0.0 
        now%dzbdy       = 0.0 
        now%H_eff       = 0.0 
        now%H_grnd      = 0.0  

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
        now%dHdt_n      = 0.0  
        now%dHdt_pred   = 0.0
        now%dHdt_corr   = 0.0
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
        if (allocated(now%z_base))      deallocate(now%z_base)

        if (allocated(now%dzsdt))       deallocate(now%dzsdt)
        if (allocated(now%dHidt))       deallocate(now%dHidt)
        if (allocated(now%bmb))         deallocate(now%bmb)
        if (allocated(now%fmb))         deallocate(now%fmb)
        if (allocated(now%mb_applied))  deallocate(now%mb_applied)
        if (allocated(now%mb_resid))    deallocate(now%mb_resid)
        
        if (allocated(now%G_advec))     deallocate(now%G_advec)
        if (allocated(now%mask_adv))    deallocate(now%mask_adv)
        
        if (allocated(now%mask_new))      deallocate(now%mask_new)
        if (allocated(now%mask_pred_new)) deallocate(now%mask_pred_new)
        if (allocated(now%mask_corr_new)) deallocate(now%mask_corr_new)
        
        if (allocated(now%eps_eff))     deallocate(now%eps_eff)
        if (allocated(now%tau_eff))     deallocate(now%tau_eff)
        if (allocated(now%calv))        deallocate(now%calv)
        if (allocated(now%calv_flt))    deallocate(now%calv_flt)
        if (allocated(now%calv_grnd))   deallocate(now%calv_grnd)
            
        if (allocated(now%dzsdx))       deallocate(now%dzsdx)
        if (allocated(now%dzsdy))       deallocate(now%dzsdy)
        if (allocated(now%dHidx))       deallocate(now%dHidx)
        if (allocated(now%dHidy))       deallocate(now%dHidy)
        if (allocated(now%dzbdx))       deallocate(now%dzbdx)
        if (allocated(now%dzbdy))       deallocate(now%dzbdy)
        
        if (allocated(now%H_eff))       deallocate(now%H_eff)
        if (allocated(now%H_grnd))      deallocate(now%H_grnd)

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
        
        if (allocated(now%dHdt_n))      deallocate(now%dHdt_n)
        if (allocated(now%dHdt_pred))   deallocate(now%dHdt_pred)
        if (allocated(now%dHdt_corr))   deallocate(now%dHdt_corr)
        if (allocated(now%H_ice_n))     deallocate(now%H_ice_n)
        if (allocated(now%H_ice_pred))  deallocate(now%H_ice_pred)
        if (allocated(now%H_ice_corr))  deallocate(now%H_ice_corr)
        
        if (allocated(now%z_srf_n))     deallocate(now%z_srf_n)
        
        if (allocated(now%H_ice_dyn))   deallocate(now%H_ice_dyn)
        if (allocated(now%f_ice_dyn))   deallocate(now%f_ice_dyn)
        
        return 

    end subroutine ytopo_dealloc
    
end module yelmo_topography
