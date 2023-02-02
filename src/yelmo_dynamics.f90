
module yelmo_dynamics

    use nml 
    use ncio 

    use yelmo_defs
    use yelmo_tools, only : calc_magnitude_from_staggered, calc_vertical_integrated_2D

    use deformation, only : calc_jacobian_vel_3D, calc_strain_rate_tensor_jac, &
                                                calc_strain_rate_tensor_jac_quad3D

    use velocity_general

    use velocity_sia 

    use velocity_ssa
    use solver_ssa_ac 
    ! use velocity_ssa_aa
    ! use solver_ssa_aa

    use velocity_l1l2 

    use velocity_diva
    ! use velocity_diva_ab 
    ! use solver_ssa_ab

    use basal_dragging  
    use grounding_line_flux 

    ! Note: 3D arrays defined such that first index (k=1) == base, and max index (k=nk) == surface 
    
    implicit none
      
    private

    public :: ydyn_par_load, ydyn_alloc, ydyn_dealloc
    public :: calc_ydyn
    public :: calc_ydyn_neff
    
contains

    subroutine calc_ydyn(dyn,tpo,mat,thrm,bnd,time)
        ! Velocity is a steady-state solution to a given set of boundary conditions (topo, material, etc).
        ! However, the time is passed to be able to treat relaxation conditions for stability, etc. 

        implicit none
        
        type(ydyn_class),   intent(INOUT) :: dyn
        type(ytopo_class),  intent(IN)    :: tpo 
        type(ymat_class),   intent(IN)    :: mat
        type(ytherm_class), intent(IN)    :: thrm 
        type(ybound_class), intent(IN)    :: bnd   
        real(prec),         intent(IN)    :: time 

        ! Local variables
        integer :: i, j, k, nx, ny, nz_aa, nz_ac   
        real(prec) :: dt 
        real(8)    :: cpu_time0, cpu_time1
        real(prec) :: model_time0, model_time1 

        real(prec), allocatable :: uxy_prev(:,:) 

        logical, parameter :: write_ssa_diagnostics = .FALSE.

        nx    = dyn%par%nx 
        ny    = dyn%par%ny 
        nz_aa = dyn%par%nz_aa 
        nz_ac = dyn%par%nz_ac 
        
        allocate(uxy_prev(nx,ny)) 


        ! Initialize time if necessary 
        if (dyn%par%time .gt. time) then 
            dyn%par%time = time
        end if 

        ! Store initial cpu time and model time for metrics later
        call yelmo_cpu_time(cpu_time0)
        model_time0 = dyn%par%time 

        ! Get time step
        dt = time - dyn%par%time 

        ! Store previous (n-1) depth-averaged horizontal velocity components
        ! (for use with higher-order ice thickness timestepping) 
        dyn%now%ux_bar_prev = dyn%now%ux_bar 
        dyn%now%uy_bar_prev = dyn%now%uy_bar 
        
        ! Store initial uxy_bar solution 
        uxy_prev = dyn%now%uxy_bar 
        
        ! ===== Calculate general variables ========================================

        ! Calculate driving stress 
        call calc_driving_stress(dyn%now%taud_acx,dyn%now%taud_acy,tpo%now%H_ice_dyn,tpo%now%f_ice_dyn,tpo%now%dzsdx, &
                                        tpo%now%dzsdy,dyn%par%dx,dyn%par%taud_lim,dyn%par%boundaries)

        if (dyn%par%taud_gl_method .ne. 0) then 
            ! Calculate driving stress specifically at the grounding line 
            ! via method of choice. 

            call calc_driving_stress_gl(dyn%now%taud_acx,dyn%now%taud_acy, &
                        tpo%now%H_ice_dyn,tpo%now%z_srf,bnd%z_bed,bnd%z_sl,tpo%now%H_grnd, &
                                      tpo%now%f_grnd,tpo%now%f_grnd_acx,tpo%now%f_grnd_acy, &
                                      dyn%par%dx,dyn%par%taud_gl_method,beta_gl_stag=1)

        end if 

        ! Calculate lateral boundary stress 
        call calc_lateral_bc_stress_2D(dyn%now%taul_int_acx,dyn%now%taul_int_acy,tpo%now%mask_frnt, &
                        tpo%now%H_ice,tpo%now%f_ice,tpo%now%z_srf,bnd%z_sl,rho_ice,rho_sw,dyn%par%boundaries)

        ! Calculate effective pressure 
        call calc_ydyn_neff(dyn,tpo,thrm,bnd)

        ! Calculate cb_tgt (cb_ref target value) - same as cb_ref, but always calculated,
        ! even if till_method=-1
        call calc_cb_ref(dyn%now%cb_tgt,bnd%z_bed,bnd%z_bed_sd,bnd%z_sl,bnd%H_sed,dyn%par%till_f_sed,dyn%par%till_sed_min, &
                                            dyn%par%till_sed_max,dyn%par%till_cf_ref,dyn%par%till_cf_min,dyn%par%till_z0, &
                                            dyn%par%till_z1,dyn%par%till_n_sd,dyn%par%till_scale,till_method=1)

        ! Update bed roughness coefficients cb_ref and c_bed (which are independent of velocity)
        call calc_cb_ref(dyn%now%cb_ref,bnd%z_bed,bnd%z_bed_sd,bnd%z_sl,bnd%H_sed,dyn%par%till_f_sed,dyn%par%till_sed_min, &
                                            dyn%par%till_sed_max,dyn%par%till_cf_ref,dyn%par%till_cf_min,dyn%par%till_z0, &
                                            dyn%par%till_z1,dyn%par%till_n_sd,dyn%par%till_scale,dyn%par%till_method)

        ! Finally calculate c_bed, which is simply c_bed = f(N_eff,cb_ref)
        call calc_c_bed(dyn%now%c_bed,dyn%now%cb_ref,dyn%now%N_eff,dyn%par%till_is_angle)

        ! ===== Calculate the 3D velocity field and helper variables =======================
        ! The variables to be obtained from these routines are:
        !     ux(:,:,:)           ! [m/a]
        !     uy(:,:,:)           ! [m/a]
        !     uz(:,:,:)           ! [m/a]
        !     ux_bar(:,:)         ! [m/a]
        !     uy_bar(:,:)         ! [m/a]
        !     ux_b(:,:)           ! [m/a]
        !     uy_b(:,:)           ! [m/a]
        !     ux_i(:,:,:)         ! [m/a]
        !     uy_i(:,:,:)         ! [m/a]
        !     taub_acx(:,:)       ! [Pa]
        !     taub_acy(:,:)       ! [Pa]
        !     beta(:,:)           ! [Pa a/m]
        !     beta_acx(:,:)       ! [Pa a/m]
        !     beta_acy(:,:)       ! [Pa a/m]
        !     beta_eff(:,:)       ! [Pa a/m]
        !     visc_eff(:,:,:)     ! [Pa a m]
        !     visc_eff_int(:,:)   ! [Pa a m]
        !     ssa_mask_acx(:,:)   ! [-]
        !     ssa_mask_acy(:,:)   ! [-]
        !     ssa_err_acx(:,:)
        !     ssa_err_acy(:,:)
        ! If a given solver does not use/calculate the variable, it is set to zero. 
        ! For the rest of Yelmo, at least these variables should be populated:
        ! ux, uy, uz, ux_bar, uy_bar, ux_b, uy_b, taub_acx, taub_acy, beta 

        select case(dyn%par%solver)

            case("fixed") 
                ! Do nothing - dynamics is fixed 

            case("sia") 
                ! SIA only 

                call calc_ydyn_hybrid(dyn,tpo,mat,thrm,bnd,use_sia=.TRUE.,use_ssa=.FALSE.)

            case("ssa") 
                ! SSA only 

                call calc_ydyn_hybrid(dyn,tpo,mat,thrm,bnd,use_sia=.FALSE.,use_ssa=.TRUE.)

            case("hybrid") 
                ! SIA+SSA

                call calc_ydyn_hybrid(dyn,tpo,mat,thrm,bnd,use_sia=.TRUE.,use_ssa=.TRUE.)

            case("diva","diva-noslip") 
                ! Depth-integrated variational approximation (DIVA) - Goldberg (2011); Lipscomb et al. (2019)

                call calc_ydyn_diva(dyn,tpo,mat,thrm,bnd)
            
            case("l1l2","l1l2-noslip")
                ! L1L2 solver

                call calc_ydyn_l1l2(dyn,tpo,mat,thrm,bnd)
            
            case DEFAULT 

                write(*,*) "calc_ydyn:: Error: ydyn solver not recognized." 
                write(*,*) "solver should be one of: ['fixed','hybrid','diva']"
                write(*,*) "solver = ", trim(dyn%par%solver) 
                stop 

        end select 

        ! Limit velocity values to avoid potential underflow errors 
        where (abs(dyn%now%ux) .lt. TOL_UNDERFLOW) dyn%now%ux = 0.0_wp 
        where (abs(dyn%now%uy) .lt. TOL_UNDERFLOW) dyn%now%uy = 0.0_wp 
        
        where (abs(dyn%now%ux_bar) .lt. TOL_UNDERFLOW) dyn%now%ux_bar = 0.0_wp 
        where (abs(dyn%now%uy_bar) .lt. TOL_UNDERFLOW) dyn%now%uy_bar = 0.0_wp 
        
        ! ===== Velocity Jacobian and strain rate tensor ===========================

        call calc_jacobian_vel_3D(dyn%now%jvel, dyn%now%ux, dyn%now%uy, dyn%now%uz, tpo%now%H_ice, tpo%now%f_ice, &
                                    tpo%now%f_grnd, tpo%now%dzsdx, tpo%now%dzsdy,tpo%now%dzbdx, tpo%now%dzbdy,   &
                                    dyn%par%zeta_aa, dyn%par%zeta_ac, dyn%par%dx, dyn%par%dy, dyn%par%boundaries)

        ! call calc_strain_rate_tensor_jac(dyn%now%strn, dyn%now%strn2D, dyn%now%jvel, tpo%now%H_ice, tpo%now%f_ice, tpo%now%f_grnd,  &
        !                                    dyn%par%zeta_aa, dyn%par%zeta_ac, dyn%par%dx, dyn%par%dy, mat%par%de_max, dyn%par%boundaries)
        call calc_strain_rate_tensor_jac_quad3D(dyn%now%strn, dyn%now%strn2D, dyn%now%jvel, tpo%now%H_ice, tpo%now%f_ice, tpo%now%f_grnd,  &
                                           dyn%par%zeta_aa, dyn%par%zeta_ac, dyn%par%dx, dyn%par%dy, mat%par%de_max, dyn%par%boundaries)
        
        ! ===== Additional diagnostic variables ====================================
        
        ! Diagnose ice flux 
        call calc_ice_flux(dyn%now%qq_acx,dyn%now%qq_acy,dyn%now%ux_bar,dyn%now%uy_bar,tpo%now%H_ice, &
                            dyn%par%dx,dyn%par%dy)
        dyn%now%qq        = calc_magnitude_from_staggered(dyn%now%qq_acx,dyn%now%qq_acy,tpo%now%f_ice,dyn%par%boundaries)

        dyn%now%taub      = calc_magnitude_from_staggered(dyn%now%taub_acx,dyn%now%taub_acy,tpo%now%f_ice,dyn%par%boundaries)
        dyn%now%taud      = calc_magnitude_from_staggered(dyn%now%taud_acx,dyn%now%taud_acy,tpo%now%f_ice,dyn%par%boundaries)

        dyn%now%uxy_b     = calc_magnitude_from_staggered(dyn%now%ux_b,dyn%now%uy_b,tpo%now%f_ice,dyn%par%boundaries)
        dyn%now%uxy_i_bar = calc_magnitude_from_staggered(dyn%now%ux_i_bar,dyn%now%uy_i_bar,tpo%now%f_ice,dyn%par%boundaries)
        dyn%now%uxy_bar   = calc_magnitude_from_staggered(dyn%now%ux_bar,dyn%now%uy_bar,tpo%now%f_ice,dyn%par%boundaries)

        do k = 1, nz_aa
            dyn%now%uxy(:,:,k) = calc_magnitude_from_staggered(dyn%now%ux(:,:,k),dyn%now%uy(:,:,k),tpo%now%f_ice,dyn%par%boundaries)
        end do 

        ! Store surface velocities for easy access too 
        dyn%now%ux_s  = dyn%now%ux(:,:,nz_aa)
        dyn%now%uy_s  = dyn%now%uy(:,:,nz_aa)
        dyn%now%uxy_s = dyn%now%uxy(:,:,nz_aa)

        ! Determine ratio of basal to surface velocity
        dyn%now%f_vbvs = calc_vel_ratio(uxy_base=dyn%now%uxy_b,uxy_srf=dyn%now%uxy_s)

        ! Finally, determine rate of velocity change 
        if (dt .ne. 0.0_wp) then 
            dyn%now%duxydt = (dyn%now%uxy_bar - uxy_prev) / dt 
        else 
            dyn%now%duxydt = 0.0_wp 
        end if 

        ! Advance ydyn timestep 
        dyn%par%time = time

        ! Calculate computational performance (model speed in kyr/hr)
        call yelmo_cpu_time(cpu_time1)
        model_time1 = dyn%par%time 
        call yelmo_calc_speed(dyn%par%speed,model_time0,model_time1,cpu_time0,cpu_time1)

        if (write_ssa_diagnostics) then 
            ! Write diagnostic output every timestep
            ! This could also be called internal to each Picard iteration,
            ! but this lets us see a snapshot after each full dynamics solve.
            call yelmo_write_init_ssa("ssa.nc",dyn%par%nx,dyn%par%ny,time)
            call write_step_2D_ssa(tpo,dyn,"ssa.nc",time)
        end if 

        return

    end subroutine calc_ydyn
    
    subroutine calc_ydyn_hybrid(dyn,tpo,mat,thrm,bnd,use_sia,use_ssa)
        ! Velocity is a steady-state solution to a given set of boundary conditions (topo, material, etc)

        implicit none
        
        type(ydyn_class),   intent(INOUT) :: dyn
        type(ytopo_class),  intent(IN)    :: tpo 
        type(ymat_class),   intent(IN)    :: mat
        type(ytherm_class), intent(IN)    :: thrm 
        type(ybound_class), intent(IN)    :: bnd   
        logical,            intent(IN)    :: use_sia 
        logical,            intent(IN)    :: use_ssa 

        ! Local variables
        integer :: iter, n_iter
        integer :: i, j, k, nx, ny, nz_aa, nz_ac   
        
        type(ssa_param_class) :: ssa_par 

        ! For vertical velocity calculation 
        real(wp), allocatable :: bmb(:,:)

        nx    = dyn%par%nx 
        ny    = dyn%par%ny 
        nz_aa = dyn%par%nz_aa 
        nz_ac = dyn%par%nz_ac 
        
        allocate(bmb(nx,ny))
        
        ! ===== Calculate 3D horizontal velocity solution via SIA + SSA algorithm ===================

        ! 1. Calculate SIA solution =====

        if (use_sia) then 
            ! Calculate SIA as normal 

            call calc_velocity_sia(dyn%now%ux_i,dyn%now%uy_i,dyn%now%ux_i_bar,dyn%now%uy_i_bar,tpo%now%H_ice, &
                                    tpo%now%f_ice,dyn%now%taud_acx,dyn%now%taud_acy,mat%now%ATT,dyn%par%zeta_aa, &
                                    dyn%par%dx,mat%par%n_glen,rho_ice,g,dyn%par%boundaries)

        else 
            ! Set all SIA terms to zero 

            dyn%now%ux_i     = 0.0_wp  
            dyn%now%uy_i     = 0.0_wp  
            dyn%now%ux_i_bar = 0.0_wp 
            dyn%now%uy_i_bar = 0.0_wp 

        end if 

        ! 2. Calculate SSA solution =====

        ! Define grid points with ssa active (uses beta from previous timestep)
        ! call set_ssa_masks(dyn%now%ssa_mask_acx,dyn%now%ssa_mask_acy,dyn%now%beta_acx,dyn%now%beta_acy, &
        !                    tpo%now%H_ice_dyn,tpo%now%f_ice_dyn,tpo%now%f_grnd_acx,tpo%now%f_grnd_acy,dyn%par%ssa_beta_max,use_ssa=.TRUE.)
        call set_ssa_masks(dyn%now%ssa_mask_acx,dyn%now%ssa_mask_acy,tpo%now%mask_frnt,tpo%now%H_ice, &
                                tpo%now%f_ice,tpo%now%f_grnd,use_ssa=.TRUE.,lateral_bc=dyn%par%ssa_lat_bc)

        if (use_ssa .and. dyn%par%use_ssa .and. &
                maxval(dyn%now%ssa_mask_acx+dyn%now%ssa_mask_acy) .gt. 0) then 
            ! Calculate SSA as normal 

            ! Set diva parameters from Yelmo settings 
            ssa_par%ssa_lis_opt    = dyn%par%ssa_lis_opt 
            ssa_par%boundaries     = dyn%par%boundaries  
            ssa_par%ssa_lateral_bc = dyn%par%ssa_lat_bc  
            ssa_par%visc_method    = dyn%par%visc_method 
            ssa_par%visc_const     = dyn%par%visc_const 
            ssa_par%beta_method    = dyn%par%beta_method 
            ssa_par%beta_const     = dyn%par%beta_const 
            ssa_par%beta_q         = dyn%par%beta_q 
            ssa_par%beta_u0        = dyn%par%beta_u0 
            ssa_par%beta_gl_scale  = dyn%par%beta_gl_scale 
            ssa_par%beta_gl_stag   = dyn%par%beta_gl_stag 
            ssa_par%beta_gl_f      = dyn%par%beta_gl_f 
            ssa_par%H_grnd_lim     = dyn%par%H_grnd_lim 
            ssa_par%beta_min       = dyn%par%beta_min
            ssa_par%eps_0          = dyn%par%eps_0  
            ssa_par%ssa_vel_max    = dyn%par%ssa_vel_max 
            ssa_par%ssa_iter_max   = dyn%par%ssa_iter_max 
            ssa_par%ssa_iter_rel   = dyn%par%ssa_iter_rel 
            ssa_par%ssa_iter_conv  = dyn%par%ssa_iter_conv 
            ssa_par%ssa_write_log  = yelmo_log

            call calc_velocity_ssa(dyn%now%ux_b,dyn%now%uy_b,dyn%now%taub_acx,dyn%now%taub_acy, &
                                      dyn%now%visc_eff,dyn%now%visc_eff_int,dyn%now%ssa_mask_acx,dyn%now%ssa_mask_acy, &
                                      dyn%now%ssa_err_acx,dyn%now%ssa_err_acy,dyn%par%ssa_iter_now,dyn%now%beta, &
                                      dyn%now%beta_acx,dyn%now%beta_acy,dyn%now%c_bed,dyn%now%taud_acx,dyn%now%taud_acy, &
                                      dyn%now%taul_int_acx,dyn%now%taul_int_acy, &
                                      tpo%now%H_ice_dyn,tpo%now%f_ice_dyn,tpo%now%H_grnd,tpo%now%f_grnd,tpo%now%f_grnd_acx,tpo%now%f_grnd_acy, &
                                      tpo%now%mask_frnt, &
                                      mat%now%ATT,dyn%par%zeta_aa,bnd%z_sl,bnd%z_bed,tpo%now%z_srf,dyn%par%dx,dyn%par%dy,mat%par%n_glen,ssa_par)
            ! call calc_velocity_ssa_aa(dyn%now%ux_b,dyn%now%uy_b,dyn%now%taub_acx,dyn%now%taub_acy, &
            !                           dyn%now%visc_eff,dyn%now%visc_eff_int,dyn%now%ssa_mask_acx,dyn%now%ssa_mask_acy, &
            !                           dyn%now%ssa_err_acx,dyn%now%ssa_err_acy,dyn%par%ssa_iter_now,dyn%now%beta, &
            !                           dyn%now%beta_acx,dyn%now%beta_acy,dyn%now%c_bed,dyn%now%taud_acx,dyn%now%taud_acy, &
            !                           tpo%now%H_ice_dyn,tpo%now%f_ice_dyn,tpo%now%H_grnd,tpo%now%f_grnd,tpo%now%f_grnd_acx,tpo%now%f_grnd_acy, &
            !                           mat%now%ATT,dyn%par%zeta_aa,bnd%z_sl,bnd%z_bed,tpo%now%z_srf,dyn%par%dx,dyn%par%dy,mat%par%n_glen,ssa_par)

        else 
            ! Set all SSA terms to zero 

            dyn%now%ux_b     = 0.0_wp 
            dyn%now%uy_b     = 0.0_wp 
            dyn%now%taub_acx = 0.0_wp 
            dyn%now%taub_acy = 0.0_wp 

        end if 

        ! Additionally, check if using SIA only, then apply SIA sliding as desired 
        if ( (use_sia .and. .not. use_ssa) .and. dyn%par%cb_sia .gt. 0.0) then 
            ! Calculate basal velocity from Weertman sliding law (Greve 1997)
                    
            ! call calc_velocity_basal_sia_00(dyn%now%ux_b,dyn%now%uy_b,dyn%now%taub_acx,dyn%now%taub_acy, &
            !                                 tpo%now%H_ice,tpo%now%dzsdx,tpo%now%dzsdy,thrm%now%f_pmp, &
            !                                 dyn%par%zeta_aa,dyn%par%dx,dyn%par%cb_sia,rho_ice,g)
            
        end if 

        ! 3. Join SIA and SSA solutions (SIA+SSA) =====

        ! Calculate the 3D horizontal velocity field (sum of shear and basal sliding)
        do k = 1, nz_aa 
            dyn%now%ux(:,:,k) = dyn%now%ux_i(:,:,k) + dyn%now%ux_b 
            dyn%now%uy(:,:,k) = dyn%now%uy_i(:,:,k) + dyn%now%uy_b 
        end do 

        ! Calculate the depth-averaged velocity too (sum of shear and basal sliding)
        dyn%now%ux_bar = dyn%now%ux_i_bar + dyn%now%ux_b 
        dyn%now%uy_bar = dyn%now%uy_i_bar + dyn%now%uy_b 


        ! 4. Set other variables to zero that are not treated with this solver =====
        dyn%now%duxdz     = 0.0_wp 
        dyn%now%duydz     = 0.0_wp 
        dyn%now%beta_eff  = 0.0_wp 

        ! ===== Calculate the vertical velocity through continuity ============================

        if (dyn%par%use_bmb) then 
            bmb = tpo%now%bmb 
        else 
            bmb = 0.0 
        end if 

        call calc_uz_3D(dyn%now%uz,dyn%now%uz_star,dyn%now%ux,dyn%now%uy,tpo%now%H_ice_dyn,tpo%now%f_ice_dyn,tpo%now%f_grnd, &
                        bnd%smb,bmb,tpo%now%dHidt,tpo%now%dzsdt,tpo%now%dzsdx,tpo%now%dzsdy, &
                        tpo%now%dzbdx,tpo%now%dzbdy,dyn%par%zeta_aa,dyn%par%zeta_ac,dyn%par%dx,dyn%par%dy,dyn%par%boundaries)
        
        return

    end subroutine calc_ydyn_hybrid

    subroutine calc_ydyn_diva(dyn,tpo,mat,thrm,bnd)
        ! Velocity is a steady-state solution to a given set of boundary conditions (topo, material, etc)

        implicit none
        
        type(ydyn_class),   intent(INOUT) :: dyn
        type(ytopo_class),  intent(IN)    :: tpo 
        type(ymat_class),   intent(IN)    :: mat
        type(ytherm_class), intent(IN)    :: thrm 
        type(ybound_class), intent(IN)    :: bnd   

        ! Local variables
        integer :: iter, n_iter
        integer :: i, j, k, nx, ny, nz_aa, nz_ac   
        logical :: no_slip 

        type(diva_param_class) :: diva_par 

        ! For vertical velocity calculation 
        real(prec), allocatable :: bmb(:,:)

        ! Determine whether basal sliding is allowed 
        if (trim(dyn%par%solver) .eq. "diva-noslip") then 
            no_slip = .TRUE. 
        else 
            no_slip = .FALSE. 
        end if 

        nx    = dyn%par%nx 
        ny    = dyn%par%ny 
        nz_aa = dyn%par%nz_aa 
        nz_ac = dyn%par%nz_ac 
        
        allocate(bmb(nx,ny))
        
        ! ===== Calculate 3D horizontal velocity solution via DIVA algorithm ===================

        ! Define grid points with ssa active (uses beta from previous timestep)
        ! call set_ssa_masks(dyn%now%ssa_mask_acx,dyn%now%ssa_mask_acy,dyn%now%beta_acx,dyn%now%beta_acy, &
        !                    tpo%now%H_ice_dyn,tpo%now%f_ice_dyn,tpo%now%f_grnd_acx,tpo%now%f_grnd_acy,dyn%par%ssa_beta_max,use_ssa=.TRUE.)
        call set_ssa_masks(dyn%now%ssa_mask_acx,dyn%now%ssa_mask_acy,tpo%now%mask_frnt,tpo%now%H_ice, &
                                tpo%now%f_ice,tpo%now%f_grnd,use_ssa=.TRUE.,lateral_bc=dyn%par%ssa_lat_bc)
        
        ! ajr: add these two statements for testing 2D flow (no flow in y-direction)
        ! Should consider whether this should be made into a parameter option of some kind,
        ! but it needs to be generalized, here it is hard-coded for diva only...
        ! dyn%now%ssa_mask_acy = 0_wp 
        ! dyn%now%uy_bar = 0.0_wp 

        ! Set diva parameters from Yelmo settings 
        diva_par%ssa_lis_opt    = dyn%par%ssa_lis_opt 
        diva_par%boundaries     = dyn%par%boundaries 
        diva_par%ssa_lateral_bc = dyn%par%ssa_lat_bc 
        diva_par%no_slip        = no_slip 
        diva_par%visc_method    = dyn%par%visc_method 
        diva_par%visc_const     = dyn%par%visc_const 
        diva_par%beta_method    = dyn%par%beta_method 
        diva_par%beta_const     = dyn%par%beta_const 
        diva_par%beta_q         = dyn%par%beta_q 
        diva_par%beta_u0        = dyn%par%beta_u0 
        diva_par%beta_gl_scale  = dyn%par%beta_gl_scale 
        diva_par%beta_gl_stag   = dyn%par%beta_gl_stag 
        diva_par%beta_gl_f      = dyn%par%beta_gl_f 
        diva_par%H_grnd_lim     = dyn%par%H_grnd_lim 
        diva_par%beta_min       = dyn%par%beta_min 
        diva_par%eps_0          = dyn%par%eps_0 
        diva_par%ssa_vel_max    = dyn%par%ssa_vel_max 
        diva_par%ssa_iter_max   = dyn%par%ssa_iter_max 
        diva_par%ssa_iter_rel   = dyn%par%ssa_iter_rel 
        diva_par%ssa_iter_conv  = dyn%par%ssa_iter_conv 
        diva_par%ssa_write_log  = yelmo_log

        call calc_velocity_diva(dyn%now%ux,dyn%now%uy,dyn%now%ux_bar,dyn%now%uy_bar, &
                                dyn%now%ux_b,dyn%now%uy_b,dyn%now%ux_i,dyn%now%uy_i, &
                                dyn%now%taub_acx,dyn%now%taub_acy,dyn%now%beta,dyn%now%beta_acx, &
                                dyn%now%beta_acy,dyn%now%beta_eff,dyn%now%de_eff,dyn%now%visc_eff, &
                                dyn%now%visc_eff_int,    &
                                dyn%now%duxdz,dyn%now%duydz,dyn%now%ssa_mask_acx,dyn%now%ssa_mask_acy,      &
                                dyn%now%ssa_err_acx,dyn%now%ssa_err_acy,dyn%par%ssa_iter_now,dyn%now%c_bed, &
                                dyn%now%taud_acx,dyn%now%taud_acy,dyn%now%taul_int_acx,dyn%now%taul_int_acy, &
                                tpo%now%H_ice_dyn,tpo%now%f_ice_dyn,tpo%now%H_grnd,   &
                                tpo%now%f_grnd,tpo%now%f_grnd_acx,tpo%now%f_grnd_acy,tpo%now%mask_frnt,mat%now%ATT, &
                                dyn%par%zeta_aa,bnd%z_sl,bnd%z_bed,tpo%now%z_srf,dyn%par%dx,dyn%par%dy,mat%par%n_glen,diva_par)
        ! call calc_velocity_diva_ab(dyn%now%ux,dyn%now%uy,dyn%now%ux_bar,dyn%now%uy_bar, &
        !                         dyn%now%ux_b,dyn%now%uy_b,dyn%now%ux_i,dyn%now%uy_i, &
        !                         dyn%now%taub_acx,dyn%now%taub_acy,dyn%now%beta,dyn%now%beta_acx, &
        !                         dyn%now%beta_acy, &
        !                         dyn%now%ux_bar_ab,dyn%now%uy_bar_ab, &
        !                         dyn%now%beta_eff,dyn%now%de_eff,dyn%now%visc_eff, &
        !                         dyn%now%visc_eff_int,    &
        !                         dyn%now%duxdz,dyn%now%duydz,dyn%now%ssa_mask_acx,dyn%now%ssa_mask_acy,      &
        !                         dyn%now%ssa_err_acx,dyn%now%ssa_err_acy,dyn%par%ssa_iter_now,dyn%now%c_bed, &
        !                         dyn%now%taud_acx,dyn%now%taud_acy,tpo%now%H_ice_dyn,tpo%now%f_ice_dyn,tpo%now%H_grnd,   &
        !                         tpo%now%f_grnd,tpo%now%f_grnd_acx,tpo%now%f_grnd_acy,mat%now%ATT, &
        !                         dyn%par%zeta_aa,bnd%z_sl,bnd%z_bed,tpo%now%z_srf,dyn%par%dx,dyn%par%dy,mat%par%n_glen,diva_par)
         
        ! Integrate from 3D shear velocity field to get depth-averaged field
        dyn%now%ux_i_bar = calc_vertical_integrated_2D(dyn%now%ux_i,dyn%par%zeta_aa)
        dyn%now%uy_i_bar = calc_vertical_integrated_2D(dyn%now%uy_i,dyn%par%zeta_aa)
          
        ! ===== Calculate the vertical velocity through continuity ============================

        if (dyn%par%use_bmb) then 
            bmb = tpo%now%bmb 
        else 
            bmb = 0.0 
        end if 

        call calc_uz_3D(dyn%now%uz,dyn%now%uz_star,dyn%now%ux,dyn%now%uy,tpo%now%H_ice_dyn,tpo%now%f_ice_dyn,tpo%now%f_grnd, &
                        bnd%smb,bmb,tpo%now%dHidt,tpo%now%dzsdt,tpo%now%dzsdx,tpo%now%dzsdy, &
                        tpo%now%dzbdx,tpo%now%dzbdy,dyn%par%zeta_aa,dyn%par%zeta_ac,dyn%par%dx,dyn%par%dy,dyn%par%boundaries)
        
        return

    end subroutine calc_ydyn_diva

    subroutine calc_ydyn_l1l2(dyn,tpo,mat,thrm,bnd)
        ! Velocity is a steady-state solution to a given set of boundary conditions (topo, material, etc)

        implicit none
        
        type(ydyn_class),   intent(INOUT) :: dyn
        type(ytopo_class),  intent(IN)    :: tpo 
        type(ymat_class),   intent(IN)    :: mat
        type(ytherm_class), intent(IN)    :: thrm 
        type(ybound_class), intent(IN)    :: bnd   

        ! Local variables
        integer :: iter, n_iter
        integer :: i, j, k, nx, ny, nz_aa, nz_ac   
        logical :: no_slip 

        type(l1l2_param_class) :: l1l2_par 

        ! For vertical velocity calculation 
        real(prec), allocatable :: bmb(:,:)

        ! Determine whether basal sliding is allowed 
        if (trim(dyn%par%solver) .eq. "l1l2-noslip") then 
            no_slip = .TRUE. 
        else 
            no_slip = .FALSE. 
        end if 
        
        nx    = dyn%par%nx 
        ny    = dyn%par%ny 
        nz_aa = dyn%par%nz_aa 
        nz_ac = dyn%par%nz_ac 
        
        allocate(bmb(nx,ny))
        
        ! ===== Calculate 3D horizontal velocity solution via DIVA algorithm ===================

        ! Define grid points with ssa active (uses beta from previous timestep)
        ! call set_ssa_masks(dyn%now%ssa_mask_acx,dyn%now%ssa_mask_acy,dyn%now%beta_acx,dyn%now%beta_acy, &
        !                    tpo%now%H_ice_dyn,tpo%now%f_ice_dyn,tpo%now%f_grnd_acx,tpo%now%f_grnd_acy,dyn%par%ssa_beta_max,use_ssa=.TRUE.)
        call set_ssa_masks(dyn%now%ssa_mask_acx,dyn%now%ssa_mask_acy,tpo%now%mask_frnt,tpo%now%H_ice, &
                                tpo%now%f_ice,tpo%now%f_grnd,use_ssa=.TRUE.,lateral_bc=dyn%par%ssa_lat_bc)
        
        ! Set diva parameters from Yelmo settings 
        l1l2_par%ssa_lis_opt    = dyn%par%ssa_lis_opt 
        l1l2_par%ssa_lateral_bc = dyn%par%ssa_lat_bc 
        l1l2_par%boundaries     = dyn%par%boundaries 
        l1l2_par%no_slip        = no_slip 
        l1l2_par%visc_method    = dyn%par%visc_method 
        l1l2_par%visc_const     = dyn%par%visc_const 
        l1l2_par%beta_method    = dyn%par%beta_method 
        l1l2_par%beta_const     = dyn%par%beta_const 
        l1l2_par%beta_q         = dyn%par%beta_q 
        l1l2_par%beta_u0        = dyn%par%beta_u0 
        l1l2_par%beta_gl_scale  = dyn%par%beta_gl_scale 
        l1l2_par%beta_gl_stag   = dyn%par%beta_gl_stag 
        l1l2_par%beta_gl_f      = dyn%par%beta_gl_f 
        l1l2_par%H_grnd_lim     = dyn%par%H_grnd_lim 
        l1l2_par%beta_min       = dyn%par%beta_min 
        l1l2_par%eps_0          = dyn%par%eps_0 
        l1l2_par%ssa_vel_max    = dyn%par%ssa_vel_max 
        l1l2_par%ssa_iter_max   = dyn%par%ssa_iter_max 
        l1l2_par%ssa_iter_rel   = dyn%par%ssa_iter_rel 
        l1l2_par%ssa_iter_conv  = dyn%par%ssa_iter_conv 
        l1l2_par%ssa_write_log  = yelmo_log

        call calc_velocity_l1l2(dyn%now%ux,dyn%now%uy,dyn%now%ux_bar,dyn%now%uy_bar, &
                                dyn%now%ux_b,dyn%now%uy_b,dyn%now%ux_i,dyn%now%uy_i, &
                                dyn%now%taub_acx,dyn%now%taub_acy,dyn%now%beta,dyn%now%beta_acx, &
                                dyn%now%beta_acy,dyn%now%visc_eff,dyn%now%visc_eff_int, &
                                dyn%now%ssa_mask_acx,dyn%now%ssa_mask_acy,dyn%now%ssa_err_acx, &
                                dyn%now%ssa_err_acy,dyn%par%ssa_iter_now,dyn%now%c_bed, &
                                dyn%now%taud_acx,dyn%now%taud_acy,dyn%now%taul_int_acx,dyn%now%taul_int_acy, &
                                tpo%now%H_ice_dyn,tpo%now%f_ice_dyn,tpo%now%H_grnd, &
                                tpo%now%f_grnd,tpo%now%f_grnd_acx,tpo%now%f_grnd_acy,tpo%now%mask_frnt,mat%now%ATT, &
                                dyn%par%zeta_aa,dyn%par%zeta_ac,bnd%z_sl,bnd%z_bed,tpo%now%z_srf,dyn%par%dx,dyn%par%dy,mat%par%n_glen,l1l2_par)
        
        ! Integrate from 3D shear velocity field to get depth-averaged field
        dyn%now%ux_i_bar = calc_vertical_integrated_2D(dyn%now%ux_i,dyn%par%zeta_aa)
        dyn%now%uy_i_bar = calc_vertical_integrated_2D(dyn%now%uy_i,dyn%par%zeta_aa)
          
        ! ===== Calculate the vertical velocity through continuity ============================

        if (dyn%par%use_bmb) then 
            bmb = tpo%now%bmb 
        else 
            bmb = 0.0 
        end if 

        call calc_uz_3D(dyn%now%uz,dyn%now%uz_star,dyn%now%ux,dyn%now%uy,tpo%now%H_ice_dyn,tpo%now%f_ice_dyn,tpo%now%f_grnd, &
                        bnd%smb,bmb,tpo%now%dHidt,tpo%now%dzsdt,tpo%now%dzsdx,tpo%now%dzsdy, &
                        tpo%now%dzbdx,tpo%now%dzbdy,dyn%par%zeta_aa,dyn%par%zeta_ac,dyn%par%dx,dyn%par%dy,dyn%par%boundaries)
        
        return

    end subroutine calc_ydyn_l1l2

!     subroutine calc_ydyn_ssa(dyn,tpo,thrm,mat,bnd)
!         ! Calculate the ssa solution via a linearized Picard iteration
!         ! over beta, visc and velocity
!         ! ajr: this routine is not used now, but 
!         ! it contains code that was testing prescribed gl-flux
!         ! parameterizations, which did not work yet.

!         implicit none
        
!         type(ydyn_class),   intent(INOUT) :: dyn
!         type(ytopo_class),  intent(IN)    :: tpo 
!         type(ytherm_class), intent(IN)    :: thrm 
!         type(ymat_class),   intent(IN)    :: mat 
!         type(ybound_class), intent(IN)    :: bnd   

!         ! Local variables
!         integer :: iter, i, j, nx, ny
!         real(prec) :: H_mid   
!         real(prec), allocatable :: ux_b_prev(:,:) 
!         real(prec), allocatable :: uy_b_prev(:,:) 
!         integer,    allocatable :: ssa_mask_acx(:,:) 
!         integer,    allocatable :: ssa_mask_acy(:,:) 

!         real(prec), allocatable :: beta_acx_prev(:,:) 
!         real(prec), allocatable :: beta_acy_prev(:,:) 

!         real(prec) :: L2_norm 

!         logical :: is_converged
!         logical :: write_ssa_diagnostics

!         is_converged          = .FALSE. 
!         write_ssa_diagnostics = .FALSE. 

!         nx    = dyn%par%nx 
!         ny    = dyn%par%ny

!         allocate(ux_b_prev(nx,ny))
!         allocate(uy_b_prev(nx,ny))
        
!         allocate(ssa_mask_acx(nx,ny))
!         allocate(ssa_mask_acy(nx,ny))

!         allocate(beta_acx_prev(nx,ny))
!         allocate(beta_acy_prev(nx,ny))
        
!         beta_acx_prev = dyn%now%beta_acx 
!         beta_acy_prev = dyn%now%beta_acy 

! !             if (tpo%now%f_grnd(18,3) .gt. 0.0) then 
! !                 write_ssa_diagnostics = .TRUE.

! !                 call yelmo_write_init_ssa("yelmo_ssa.nc",time_init=1.0) 
! !             end if 
        
!         if (write_ssa_diagnostics) then 
!             call yelmo_write_init_ssa("yelmo_ssa.nc",nx,ny,time_init=1.0_wp)
!         end if 

!         ! Store original ssa mask 
!         ssa_mask_acx = dyn%now%ssa_mask_acx
!         ssa_mask_acy = dyn%now%ssa_mask_acy
        
!         ! Initially set error very high 
!         dyn%now%ssa_err_acx = 1.0_wp 
!         dyn%now%ssa_err_acy = 1.0_wp 

!         do iter = 1, dyn%par%ssa_iter_max

!             ! Store previous solution 
!             ux_b_prev = dyn%now%ux_b 
!             uy_b_prev = dyn%now%uy_b 

!             !   1. Calculate basal drag coefficient beta (beta, beta_acx, beta_acy) 

! !             call calc_ydyn_beta(dyn,tpo,mat,bnd)

!             call calc_beta(dyn%now%beta,dyn%now%c_bed,dyn%now%ux_b,dyn%now%uy_b,tpo%now%H_ice,tpo%now%f_ice,tpo%now%H_grnd, &
!                             tpo%now%f_grnd,bnd%z_bed,bnd%z_sl,dyn%par%beta_method, &
!                                 dyn%par%beta_const,dyn%par%beta_q,dyn%par%beta_u0,dyn%par%beta_gl_scale, &
!                                 dyn%par%beta_gl_f,dyn%par%H_grnd_lim,dyn%par%beta_min,dyn%par%boundaries)

!             ! Stagger beta
!             call stagger_beta(dyn%now%beta_acx,dyn%now%beta_acy,dyn%now%beta,tpo%now%H_ice,tpo%now%f_ice,dyn%now%ux_b,dyn%now%uy_b,tpo%now%f_grnd, &
!                             tpo%now%f_grnd_acx,tpo%now%f_grnd_acy,dyn%par%beta_gl_stag,dyn%par%beta_min,dyn%par%boundaries)

!             !   2. Calculate effective viscosity
            
!             ! Note: disable shear contribution to viscosity for this solver, for mixed terms use hybrid-pd12 option.
!             ! Note: Here visc_eff_int is calculated using ux_b and uy_b (ssa velocity), not ux_bar/uy_bar as in hybrid-pd12. 
!             dyn%now%visc_eff_int = calc_visc_eff_2D(dyn%now%ux_b,dyn%now%uy_b,dyn%now%duxdz_bar*0.0,dyn%now%duydz_bar*0.0, &
!                                                     tpo%now%H_ice,mat%now%ATT,dyn%par%zeta_aa,dyn%par%dx,dyn%par%dy,mat%par%n_glen)
            
!             !   X. Prescribe grounding-line flux 
! if (.FALSE.) then
!             ! Testing prescribed grounding-line flux/vel - experimental!!!

!             ! Calculate the analytical grounding-line flux 
!             call calc_grounding_line_flux(dyn%now%qq_gl_acx,dyn%now%qq_gl_acy,tpo%now%H_ice,mat%now%ATT_bar, &
!                         dyn%now%c_bed,dyn%now%ux_b,dyn%now%uy_b,tpo%now%f_grnd,tpo%now%f_grnd_acx,tpo%now%f_grnd_acy, &
!                         mat%par%n_glen,dyn%par%beta_q,Q0=0.61_wp,f_drag=0.6_wp,glf_method="power")

!             ! Where qq_gl is present, prescribe velocity and set mask to -1

!             ! Restore original ssa mask (without grounding line flags)
!             dyn%now%ssa_mask_acx = ssa_mask_acx
!             dyn%now%ssa_mask_acy = ssa_mask_acy
            
!             write(*,*) "glf"

!             ! acx nodes 
!             do j = 1, ny 
!             do i = 1, nx-1

!                 H_mid = 0.5*(tpo%now%H_ice(i,j)+tpo%now%H_ice(i+1,j))
                
!                 if (dyn%now%qq_gl_acx(i,j) .ne. 0.0 .and. H_mid .gt. 0.0) then 
!                     ! Prescribe velocity at this point 

!                     if (j == 3) then 
!                         write(*,*) "glf", i, dyn%now%ux_b(i,j), dyn%now%qq_gl_acx(i,j) / H_mid
!                     end if 
                    
! !                     dyn%now%ux_b(i,j) = dyn%now%qq_gl_acx(i,j) / H_mid 
! !                     dyn%now%ssa_mask_acx(i,j) = -1

!                 end if 

!             end do 
!             end do 

!             ! acy nodes 
!             do j = 1, ny-1 
!             do i = 1, nx

!                 H_mid = 0.5*(tpo%now%H_ice(i,j)+tpo%now%H_ice(i,j+1))
                
!                 if (dyn%now%qq_gl_acy(i,j) .ne. 0.0 .and. H_mid .gt. 0.0) then 
!                     ! Prescribe velocity at this point 

!                     dyn%now%uy_b(i,j) = dyn%now%qq_gl_acy(i,j) / H_mid 
!                     dyn%now%ssa_mask_acy(i,j) = -1
                    
!                 end if 

!             end do 
!             end do
! end if 

!             !   3. Calculate SSA solution

! if (.TRUE.) then 
!             if (iter .gt. 1) then
!                 ! Update ssa mask based on convergence with previous step to reduce calls 
!                 call update_ssa_mask_convergence(dyn%now%ssa_mask_acx,dyn%now%ssa_mask_acy, &
!                                                 dyn%now%ssa_err_acx,dyn%now%ssa_err_acy,err_lim=real(1e-3,prec)) 
!             end if 
! end if 

!             ! Call ssa solver to determine ux_b/uy_b, where ssa_mask_acx/y are > 0
!             call calc_vxy_ssa_matrix(dyn%now%ux_b,dyn%now%uy_b,L2_norm,dyn%now%beta_acx,dyn%now%beta_acy,dyn%now%visc_eff_int, &
!                                      dyn%now%ssa_mask_acx,dyn%now%ssa_mask_acy,tpo%now%H_ice,tpo%now%f_ice,dyn%now%taud_acx, &
!                                      dyn%now%taud_acy,tpo%now%H_grnd,bnd%z_sl,bnd%z_bed,dyn%par%dx,dyn%par%dy, &
!                                      dyn%par%ssa_vel_max,dyn%par%boundaries,dyn%par%ssa_lis_opt)

!             ! Apply relaxation to keep things stable
!             call relax_ssa(dyn%now%ux_b,dyn%now%uy_b,ux_b_prev,uy_b_prev,rel=dyn%par%ssa_iter_rel)

!             ! Check for convergence
!             is_converged = check_vel_convergence_l2rel(dyn%now%ux_b,dyn%now%uy_b,ux_b_prev,uy_b_prev, &
!                                         dyn%now%ssa_mask_acx.gt.0,dyn%now%ssa_mask_acy.gt.0, &
!                                         dyn%par%ssa_iter_conv,iter,dyn%par%ssa_iter_max,yelmo_log,use_L2_norm=.FALSE.)

!             ! Calculate an L1 error metric over matrix for diagnostics
!             call check_vel_convergence_l1rel_matrix(dyn%now%ssa_err_acx,dyn%now%ssa_err_acy,dyn%now%ux_b,dyn%now%uy_b, &
!                                                                                             ux_b_prev,uy_b_prev)

            
!             if (write_ssa_diagnostics) then  
!                 !call write_step_2D_ssa(tpo,dyn,"yelmo_ssa.nc",ux_b_prev,uy_b_prev,time=real(iter,prec))    
!             end if 

!             ! Exit iterations if ssa solution has converged
!             if (is_converged) exit 

!         end do 
!         ! == END iterations ==

! !         if (write_ssa_diagnostics) then 
! !             stop 
! !         end if 

!         return 

!     end subroutine calc_ydyn_ssa
    
    subroutine calc_ydyn_neff(dyn,tpo,thrm,bnd)
        ! Update N_eff based on parameter choices

        implicit none
        
        type(ydyn_class),   intent(INOUT) :: dyn
        type(ytopo_class),  intent(IN)    :: tpo 
        type(ytherm_class), intent(IN)    :: thrm
        type(ybound_class), intent(IN)    :: bnd  

        ! Local variables 
        real(prec), allocatable :: H_w(:,:) 

        ! Allocate local H_w variable to represent water layer thickness if not available
        allocate(H_w(dyn%par%nx,dyn%par%ny)) 

        ! Determine whether to use actual water layer thickness or parameterized layer thickness
        if (dyn%par%neff_set_water) then
            ! Set water to maximum thickness for temperate ice
            H_w = thrm%par%H_w_max * thrm%now%f_pmp  
        else 
            ! Use boundary water thickness field
            H_w = thrm%now%H_w 
        end if 
        
        ! Calculate effective pressure N_eff [Pa]
        select case(dyn%par%neff_method)

            case(-1) 
                ! Do nothing, effective pressure is calculated externally 

            case(0)
                ! Constant value (to scale friction coefficients)

                dyn%now%N_eff = dyn%par%neff_const 

            case(1)
                ! Effective pressure == overburden pressure 

                dyn%now%N_eff = calc_effective_pressure_overburden(tpo%now%H_ice_dyn,tpo%now%f_ice_dyn,tpo%now%f_grnd)

            case(2) 
                ! Effective pressure diminishes with marine character
                ! following Leguy et al. (2014) 

                dyn%now%N_eff = calc_effective_pressure_marine(tpo%now%H_ice_dyn,tpo%now%f_ice_dyn,bnd%z_bed,bnd%z_sl, &
                                                                                            H_w,p=dyn%par%neff_p)

            case(3)
                ! Effective pressure as basal till pressure
                ! following van Pelt and Bueler (2015)

                call calc_effective_pressure_till(dyn%now%N_eff,H_w,tpo%now%H_ice_dyn,tpo%now%f_ice_dyn,tpo%now%f_grnd, &
                                  thrm%par%H_w_max,dyn%par%neff_N0,dyn%par%neff_delta,dyn%par%neff_e0,dyn%par%neff_Cc) 

            case(4) 
                ! Calculate two-valued effective pressure using till parameter neff_delta 

                call calc_effective_pressure_two_value(dyn%now%N_eff,thrm%now%f_pmp,tpo%now%H_ice_dyn,tpo%now%f_ice_dyn, &
                                                                                        tpo%now%f_grnd,dyn%par%neff_delta)
                
            case DEFAULT 

                write(*,*) "ydyn_calc_Neff:: Error: neff_method not recognized, must be one of [-1,0,1,2,3]."
                write(*,*) "neff_method = ", dyn%par%neff_method 
                stop 

        end select 
        

        return 

    end subroutine calc_ydyn_neff

    subroutine ydyn_par_load(par,filename,zeta_aa,zeta_ac,nx,ny,dx,init)

        type(ydyn_param_class), intent(OUT) :: par
        character(len=*),       intent(IN)  :: filename
        real(prec),             intent(IN)  :: zeta_aa(:)
        real(prec),             intent(IN)  :: zeta_ac(:)
        integer,                intent(IN)  :: nx, ny 
        real(prec),             intent(IN)  :: dx   
        logical, optional,      intent(IN)  :: init 

        ! Local variables 
        logical :: init_pars 

        init_pars = .FALSE.
        if (present(init)) init_pars = .TRUE. 
        
        call nml_read(filename,"ydyn","solver",             par%solver,             init=init_pars)
        call nml_read(filename,"ydyn","visc_method",        par%visc_method,        init=init_pars)
        call nml_read(filename,"ydyn","visc_const",         par%visc_const,         init=init_pars)
        call nml_read(filename,"ydyn","beta_method",        par%beta_method,        init=init_pars)
        call nml_read(filename,"ydyn","beta_const",         par%beta_const,         init=init_pars)
        call nml_read(filename,"ydyn","beta_q",             par%beta_q,             init=init_pars)
        call nml_read(filename,"ydyn","beta_u0",            par%beta_u0,            init=init_pars)
        call nml_read(filename,"ydyn","beta_gl_scale",      par%beta_gl_scale,      init=init_pars)
        call nml_read(filename,"ydyn","beta_gl_stag",       par%beta_gl_stag,       init=init_pars)
        call nml_read(filename,"ydyn","beta_gl_f",          par%beta_gl_f,          init=init_pars)
        call nml_read(filename,"ydyn","taud_gl_method",     par%taud_gl_method,     init=init_pars)
        call nml_read(filename,"ydyn","H_grnd_lim",         par%H_grnd_lim,         init=init_pars)
        call nml_read(filename,"ydyn","n_sm_beta",          par%n_sm_beta,          init=init_pars)
        call nml_read(filename,"ydyn","beta_min",           par%beta_min,           init=init_pars)
        call nml_read(filename,"ydyn","eps_0",              par%eps_0,              init=init_pars)
        call nml_read(filename,"ydyn","ssa_lis_opt",        par%ssa_lis_opt,        init=init_pars)
        call nml_read(filename,"ydyn","ssa_lat_bc",         par%ssa_lat_bc,         init=init_pars)
        call nml_read(filename,"ydyn","ssa_beta_max",       par%ssa_beta_max,       init=init_pars)
        call nml_read(filename,"ydyn","ssa_vel_max",        par%ssa_vel_max,        init=init_pars)
        call nml_read(filename,"ydyn","ssa_iter_max",       par%ssa_iter_max,       init=init_pars)
        call nml_read(filename,"ydyn","ssa_iter_rel",       par%ssa_iter_rel,       init=init_pars)
        call nml_read(filename,"ydyn","ssa_iter_conv",      par%ssa_iter_conv,      init=init_pars)
        call nml_read(filename,"ydyn","taud_lim",           par%taud_lim,           init=init_pars)
        call nml_read(filename,"ydyn","cb_sia",             par%cb_sia,             init=init_pars)
        
        call nml_read(filename,"ytill","method",            par%till_method,        init=init_pars)
        call nml_read(filename,"ytill","scale",             par%till_scale,         init=init_pars)
        call nml_read(filename,"ytill","is_angle",          par%till_is_angle,      init=init_pars)
        call nml_read(filename,"ytill","n_sd",              par%till_n_sd,          init=init_pars)
        call nml_read(filename,"ytill","f_sed",             par%till_f_sed,         init=init_pars)
        call nml_read(filename,"ytill","sed_min",           par%till_sed_min,       init=init_pars)
        call nml_read(filename,"ytill","sed_max",           par%till_sed_max,       init=init_pars)
        call nml_read(filename,"ytill","z0",                par%till_z0,            init=init_pars)
        call nml_read(filename,"ytill","z1",                par%till_z1,            init=init_pars)
        call nml_read(filename,"ytill","cf_min",            par%till_cf_min,        init=init_pars)
        call nml_read(filename,"ytill","cf_ref",            par%till_cf_ref,        init=init_pars)
        
        call nml_read(filename,"yneff","method",            par%neff_method,        init=init_pars)
        call nml_read(filename,"yneff","const",             par%neff_const,         init=init_pars)
        call nml_read(filename,"yneff","p",                 par%neff_p,             init=init_pars)
        call nml_read(filename,"yneff","set_water",         par%neff_set_water,     init=init_pars)
        call nml_read(filename,"yneff","N0",                par%neff_N0,            init=init_pars)
        call nml_read(filename,"yneff","delta",             par%neff_delta,         init=init_pars)
        call nml_read(filename,"yneff","e0",                par%neff_e0,            init=init_pars)
        call nml_read(filename,"yneff","Cc",                par%neff_Cc,            init=init_pars)

        ! === Set internal parameters ======

        par%nx    = nx 
        par%ny    = ny 
        par%dx    = dx 
        par%dy    = dx 
        par%nz_aa = size(zeta_aa,1)  
        par%nz_ac = size(zeta_ac,1)

        if (allocated(par%zeta_aa)) deallocate(par%zeta_aa)
        allocate(par%zeta_aa(par%nz_aa))
        par%zeta_aa = zeta_aa 
        
        if (allocated(par%zeta_ac)) deallocate(par%zeta_ac)
        allocate(par%zeta_ac(par%nz_ac))
        par%zeta_ac = zeta_ac 
        
        ! Define how boundaries of grid should be treated 
        ! This should only be modified by the dom%par%experiment variable
        ! in yelmo_init. By default set boundaries to zero 
        par%boundaries = "zeros" 
        
        ! Set use_bmb to decide whether bmb should enter into the calculation of vertical velocity from continuity
        ! By default, this is true, but it will automatically be set to match ydyn%par%use_bmb 
        par%use_bmb = .TRUE. 

        ! By default use ssa too, unless otherwise desired (can be set externally)
        par%use_ssa  = .TRUE. 

        ! Define current time as unrealistic value
        par%time = 1000000000   ! [a] 1 billion years in the future
        
        ! Set ssa_iter_now = 1 to start 
        par%ssa_iter_now = 1 

        return

    end subroutine ydyn_par_load

    subroutine ydyn_alloc(now,nx,ny,nz_aa,nz_ac)

        implicit none 

        type(ydyn_state_class), intent(INOUT) :: now 
        integer, intent(IN) :: nx, ny, nz_aa, nz_ac   

        call ydyn_dealloc(now)

        allocate(now%ux_bar_ab(nx,ny)) 
        allocate(now%uy_bar_ab(nx,ny))

        allocate(now%ux(nx,ny,nz_aa)) 
        allocate(now%uy(nx,ny,nz_aa)) 
        allocate(now%uxy(nx,ny,nz_aa)) 
        allocate(now%uz(nx,ny,nz_ac)) 
        allocate(now%uz_star(nx,ny,nz_ac)) 

        allocate(now%ux_bar(nx,ny)) 
        allocate(now%uy_bar(nx,ny))
        allocate(now%uxy_bar(nx,ny))
        
        allocate(now%ux_bar_prev(nx,ny)) 
        allocate(now%uy_bar_prev(nx,ny))

        allocate(now%ux_b(nx,ny)) 
        allocate(now%uy_b(nx,ny))
        allocate(now%uxy_b(nx,ny))

        allocate(now%ux_s(nx,ny)) 
        allocate(now%uy_s(nx,ny))
        allocate(now%uxy_s(nx,ny))
        
        allocate(now%ux_i(nx,ny,nz_aa)) 
        allocate(now%uy_i(nx,ny,nz_aa))
        allocate(now%ux_i_bar(nx,ny)) 
        allocate(now%uy_i_bar(nx,ny))
        allocate(now%uxy_i_bar(nx,ny))

        allocate(now%duxydt(nx,ny))

        allocate(now%duxdz(nx,ny,nz_aa)) 
        allocate(now%duydz(nx,ny,nz_aa))
        allocate(now%duxdz_bar(nx,ny)) 
        allocate(now%duydz_bar(nx,ny))

        allocate(now%taud_acx(nx,ny)) 
        allocate(now%taud_acy(nx,ny)) 
        allocate(now%taud(nx,ny)) 
        
        allocate(now%taub_acx(nx,ny)) 
        allocate(now%taub_acy(nx,ny)) 
        allocate(now%taub(nx,ny)) 

        allocate(now%taul_int_acx(nx,ny)) 
        allocate(now%taul_int_acy(nx,ny)) 

        allocate(now%qq_gl_acx(nx,ny)) 
        allocate(now%qq_gl_acy(nx,ny)) 

        allocate(now%qq_acx(nx,ny)) 
        allocate(now%qq_acy(nx,ny)) 
        allocate(now%qq(nx,ny)) 

        allocate(now%de_eff(nx,ny,nz_aa))  
        allocate(now%visc_eff(nx,ny,nz_aa))  
        allocate(now%visc_eff_int(nx,ny))

        allocate(now%cb_tgt(nx,ny))
        allocate(now%cb_ref(nx,ny))
        allocate(now%c_bed(nx,ny)) 
        
        allocate(now%N_eff(nx,ny))

        allocate(now%beta_acx(nx,ny))
        allocate(now%beta_acy(nx,ny))
        allocate(now%beta(nx,ny))
        allocate(now%beta_eff(nx,ny))

        allocate(now%f_vbvs(nx,ny)) 
        
        allocate(now%ssa_mask_acx(nx,ny)) 
        allocate(now%ssa_mask_acy(nx,ny)) 
        allocate(now%ssa_err_acx(nx,ny)) 
        allocate(now%ssa_err_acy(nx,ny)) 
        
        allocate(now%jvel%dxx(nx,ny,nz_aa))
        allocate(now%jvel%dxy(nx,ny,nz_aa))
        allocate(now%jvel%dxz(nx,ny,nz_aa))
        allocate(now%jvel%dyx(nx,ny,nz_aa))
        allocate(now%jvel%dyy(nx,ny,nz_aa))
        allocate(now%jvel%dyz(nx,ny,nz_aa))
        allocate(now%jvel%dzx(nx,ny,nz_ac))
        allocate(now%jvel%dzy(nx,ny,nz_ac))
        allocate(now%jvel%dzz(nx,ny,nz_ac))

        allocate(now%strn%dxx(nx,ny,nz_aa))
        allocate(now%strn%dyy(nx,ny,nz_aa))
        allocate(now%strn%dxy(nx,ny,nz_aa))
        allocate(now%strn%dxz(nx,ny,nz_aa))
        allocate(now%strn%dyz(nx,ny,nz_aa))
        allocate(now%strn%div(nx,ny,nz_aa))
        allocate(now%strn%de(nx,ny,nz_aa))
        allocate(now%strn%f_shear(nx,ny,nz_aa))
        
        allocate(now%strn2D%dxx(nx,ny))
        allocate(now%strn2D%dyy(nx,ny))
        allocate(now%strn2D%dxy(nx,ny))
        allocate(now%strn2D%dxz(nx,ny))
        allocate(now%strn2D%dyz(nx,ny))
        allocate(now%strn2D%div(nx,ny))
        allocate(now%strn2D%de(nx,ny))
        allocate(now%strn2D%f_shear(nx,ny))
        allocate(now%strn2D%eps_eig_1(nx,ny))
        allocate(now%strn2D%eps_eig_2(nx,ny))
        
        ! Set all variables to zero intially
        now%ux_bar_ab         = 0.0 
        now%uy_bar_ab         = 0.0

        now%ux                = 0.0 
        now%uy                = 0.0 
        now%uxy               = 0.0 
        now%uz                = 0.0 
        now%uz_star           = 0.0

        now%ux_bar            = 0.0 
        now%uy_bar            = 0.0
        now%uxy_bar           = 0.0

        now%ux_bar_prev       = 0.0 
        now%uy_bar_prev       = 0.0
        
        now%ux_b              = 0.0 
        now%uy_b              = 0.0
        now%uxy_b             = 0.0

        now%ux_s              = 0.0 
        now%uy_s              = 0.0
        now%uxy_s             = 0.0
        
        now%ux_i              = 0.0 
        now%uy_i              = 0.0
        now%ux_i_bar          = 0.0 
        now%uy_i_bar          = 0.0
        now%uxy_i_bar         = 0.0
        
        now%duxydt            = 0.0

        now%duxdz             = 0.0 
        now%duydz             = 0.0
        now%duxdz_bar         = 0.0 
        now%duydz_bar         = 0.0

        now%taud_acx          = 0.0 
        now%taud_acy          = 0.0 
        now%taud              = 0.0 
        
        now%taub_acx          = 0.0 
        now%taub_acy          = 0.0 
        now%taub              = 0.0 
        
        now%taul_int_acx    = 0.0 
        now%taul_int_acy    = 0.0 

        now%qq_gl_acx         = 0.0 
        now%qq_gl_acy         = 0.0 
        
        now%qq_acx            = 0.0 
        now%qq_acy            = 0.0 
        now%qq                = 0.0 
        
        now%de_eff            = 0.0 
        now%visc_eff          = 1e3  
        now%visc_eff_int      = 1e3  
            
        now%cb_tgt            = 0.0
        now%cb_ref            = 0.0
        now%c_bed             = 0.0 
        
        now%N_eff             = 0.0 

        now%beta_acx          = 0.0 
        now%beta_acy          = 0.0 
        now%beta              = 0.0         
        now%beta_eff          = 0.0 

        now%f_vbvs            = 0.0 

        now%ssa_mask_acx      = 0.0 
        now%ssa_mask_acy      = 0.0 
        now%ssa_err_acx       = 0.0 
        now%ssa_err_acy       = 0.0 
        
        now%jvel%dxx            = 0.0
        now%jvel%dxy            = 0.0
        now%jvel%dxz            = 0.0
        now%jvel%dyx            = 0.0
        now%jvel%dyy            = 0.0
        now%jvel%dyz            = 0.0
        now%jvel%dzx            = 0.0
        now%jvel%dzy            = 0.0
        now%jvel%dzz            = 0.0
        
        now%strn%dxx     = 0.0 
        now%strn%dyy     = 0.0 
        now%strn%dxy     = 0.0 
        now%strn%dxz     = 0.0
        now%strn%dyz     = 0.0
        now%strn%div     = 0.0
        now%strn%de      = 0.0
        now%strn%f_shear = 0.0 
        
        now%strn2D%dxx   = 0.0 
        now%strn2D%dyy   = 0.0 
        now%strn2D%dxy   = 0.0
        now%strn2D%dxz   = 0.0
        now%strn2D%dyz   = 0.0
        now%strn2D%div   = 0.0
        now%strn2D%de    = 0.0 
        now%strn2D%eps_eig_1 = 0.0 
        now%strn2D%eps_eig_2 = 0.0 
        
        return 

    end subroutine ydyn_alloc

    subroutine ydyn_dealloc(now)

        implicit none 

        type(ydyn_state_class), intent(INOUT) :: now

        if (allocated(now%ux_bar_ab))       deallocate(now%ux_bar_ab) 
        if (allocated(now%uy_bar_ab))       deallocate(now%uy_bar_ab)
        
        if (allocated(now%ux))              deallocate(now%ux) 
        if (allocated(now%uy))              deallocate(now%uy) 
        if (allocated(now%uxy))             deallocate(now%uxy) 
        if (allocated(now%uz))              deallocate(now%uz) 
        if (allocated(now%uz_star))         deallocate(now%uz_star) 

        if (allocated(now%ux_bar))          deallocate(now%ux_bar) 
        if (allocated(now%uy_bar))          deallocate(now%uy_bar)
        if (allocated(now%uxy_bar))         deallocate(now%uxy_bar)
        
        if (allocated(now%ux_bar_prev))     deallocate(now%ux_bar_prev) 
        if (allocated(now%uy_bar_prev))     deallocate(now%uy_bar_prev)
        
        if (allocated(now%ux_b))            deallocate(now%ux_b) 
        if (allocated(now%uy_b))            deallocate(now%uy_b)
        if (allocated(now%uxy_b))           deallocate(now%uxy_b)
        
        if (allocated(now%ux_s))            deallocate(now%ux_s) 
        if (allocated(now%uy_s))            deallocate(now%uy_s)
        if (allocated(now%uxy_s))           deallocate(now%uxy_s)
        
        if (allocated(now%ux_i))            deallocate(now%ux_i) 
        if (allocated(now%uy_i))            deallocate(now%uy_i)

        if (allocated(now%ux_i_bar))        deallocate(now%ux_i_bar) 
        if (allocated(now%uy_i_bar))        deallocate(now%uy_i_bar)
        if (allocated(now%uxy_i_bar))       deallocate(now%uxy_i_bar)
        
        if (allocated(now%duxydt ))         deallocate(now%duxydt)
        
        if (allocated(now%duxdz))           deallocate(now%duxdz) 
        if (allocated(now%duydz))           deallocate(now%duydz)
        if (allocated(now%duxdz_bar))       deallocate(now%duxdz_bar) 
        if (allocated(now%duydz_bar))       deallocate(now%duydz_bar)

        if (allocated(now%taud_acx))        deallocate(now%taud_acx) 
        if (allocated(now%taud_acy))        deallocate(now%taud_acy) 
        if (allocated(now%taud))            deallocate(now%taud) 
        
        if (allocated(now%taub_acx))        deallocate(now%taub_acx) 
        if (allocated(now%taub_acy))        deallocate(now%taub_acy) 
        if (allocated(now%taub))            deallocate(now%taub) 
        
        if (allocated(now%taul_int_acx))  deallocate(now%taul_int_acx) 
        if (allocated(now%taul_int_acy))  deallocate(now%taul_int_acy) 
        
        if (allocated(now%qq_gl_acx))       deallocate(now%qq_gl_acx) 
        if (allocated(now%qq_gl_acy))       deallocate(now%qq_gl_acy) 
        
        if (allocated(now%qq_acx))          deallocate(now%qq_acx) 
        if (allocated(now%qq_acy))          deallocate(now%qq_acy) 
        if (allocated(now%qq))              deallocate(now%qq) 
        
        if (allocated(now%de_eff))          deallocate(now%de_eff) 
        if (allocated(now%visc_eff))        deallocate(now%visc_eff) 
        if (allocated(now%visc_eff_int))    deallocate(now%visc_eff_int) 
        
        if (allocated(now%cb_tgt))          deallocate(now%cb_tgt) 
        if (allocated(now%cb_ref))          deallocate(now%cb_ref) 
        if (allocated(now%c_bed))           deallocate(now%c_bed) 
        
        if (allocated(now%N_eff))           deallocate(now%N_eff)
        
        if (allocated(now%beta_acx))        deallocate(now%beta_acx) 
        if (allocated(now%beta_acy))        deallocate(now%beta_acy) 
        if (allocated(now%beta))            deallocate(now%beta)         
        if (allocated(now%beta_eff))        deallocate(now%beta_eff) 

        if (allocated(now%f_vbvs))          deallocate(now%f_vbvs) 

        if (allocated(now%ssa_mask_acx))    deallocate(now%ssa_mask_acx) 
        if (allocated(now%ssa_mask_acy))    deallocate(now%ssa_mask_acy) 
        if (allocated(now%ssa_err_acx))     deallocate(now%ssa_err_acx) 
        if (allocated(now%ssa_err_acy))     deallocate(now%ssa_err_acy) 

        if (allocated(now%jvel%dxx))         deallocate(now%jvel%dxx)
        if (allocated(now%jvel%dxy))         deallocate(now%jvel%dxy)
        if (allocated(now%jvel%dxz))         deallocate(now%jvel%dxz)
        if (allocated(now%jvel%dyx))         deallocate(now%jvel%dyx)
        if (allocated(now%jvel%dyy))         deallocate(now%jvel%dyy)
        if (allocated(now%jvel%dyz))         deallocate(now%jvel%dyz)
        if (allocated(now%jvel%dzx))         deallocate(now%jvel%dzx)
        if (allocated(now%jvel%dzy))         deallocate(now%jvel%dzy)
        if (allocated(now%jvel%dzz))         deallocate(now%jvel%dzz)

        if (allocated(now%strn%dxx))        deallocate(now%strn%dxx)
        if (allocated(now%strn%dyy))        deallocate(now%strn%dyy)
        if (allocated(now%strn%dxy))        deallocate(now%strn%dxy)
        if (allocated(now%strn%dxz))        deallocate(now%strn%dxz)
        if (allocated(now%strn%dyz))        deallocate(now%strn%dyz)
        if (allocated(now%strn%div))        deallocate(now%strn%div)
        if (allocated(now%strn%de))         deallocate(now%strn%de)
        if (allocated(now%strn%f_shear))    deallocate(now%strn%f_shear)
        
        if (allocated(now%strn2D%dxx))      deallocate(now%strn2D%dxx)
        if (allocated(now%strn2D%dyy))      deallocate(now%strn2D%dyy)
        if (allocated(now%strn2D%dxy))      deallocate(now%strn2D%dxy)
        if (allocated(now%strn2D%dxz))      deallocate(now%strn2D%dxz)
        if (allocated(now%strn2D%dyz))      deallocate(now%strn2D%dyz)
        if (allocated(now%strn2D%div))      deallocate(now%strn2D%div)
        if (allocated(now%strn2D%de))       deallocate(now%strn2D%de)
        if (allocated(now%strn2D%eps_eig_1)) deallocate(now%strn2D%eps_eig_1)
        if (allocated(now%strn2D%eps_eig_2)) deallocate(now%strn2D%eps_eig_2)
        
        return 

    end subroutine ydyn_dealloc
    
    subroutine ydyn_set_borders(ux,uy,boundaries)

        implicit none 

        real(prec),       intent(INOUT) :: ux(:,:) 
        real(prec),       intent(INOUT) :: uy(:,:) 
        character(len=*), intent(IN)    :: boundaries 

        ! Local variables 
        integer :: nx, ny 

        nx = size(ux,1)
        ny = size(ux,2) 

        ! Post processing of velocity field ================

        if (.TRUE.) then 
            ! ajr: do not use yet, not well tested 

        select case(trim(boundaries))

            case("zeros","EISMINT")

                ! Border values are zero by default, do nothing 

            case("periodic") 

                ux(1,:)  = ux(nx-1,:) 
                ux(nx,:) = ux(2,:) 
                ux(:,1)  = ux(:,ny-1)
                ux(:,ny) = ux(:,2) 

                uy(1,:)  = uy(nx-1,:) 
                uy(nx,:) = uy(2,:) 
                uy(:,1)  = uy(:,ny-1)
                uy(:,ny) = uy(:,2) 

            case("MISMIP3D")

                ! === MISMIP3D =====

                ! x=0, dome - zero velocity 
                ux(1,:)    = 0.0       
                uy(1,:)    = 0.0 

                ! x=800km, no ice - zero by default 
                ux(nx,:)   = 0.0 
                uy(nx,:)   = 0.0 

                ! y=-50km, free-slip condition, no tangential velocity   
                uy(:,1)    = 0.0 

                ! y=50km, free-slip condition, no tangential velocity  
                uy(:,ny)     = 0.0 

            case("infinite")
                ! ajr: we should check setting border H values equal to inner neighbors
                
                write(*,*) "calc_ice_thickness:: error: boundary method not implemented yet: "//trim(boundaries)
                write(*,*) "TO DO!"
                stop 

            case DEFAULT 

                write(*,*) "calc_ice_thickness:: error: boundary method not recognized: "//trim(boundaries)
                stop 

        end select 
        
        end if 

        return 

    end subroutine ydyn_set_borders

    subroutine yelmo_write_init_ssa(filename,nx,ny,time_init)

        implicit none 

        character(len=*),  intent(IN) :: filename 
        integer,           intent(IN) :: nx 
        integer,           intent(IN) :: ny
        real(prec),        intent(IN) :: time_init

        ! Initialize netcdf file and dimensions
        call nc_create(filename)
        call nc_write_dim(filename,"xc",     x=0.0_wp,dx=1.0_wp,nx=nx,units="gridpoints")
        call nc_write_dim(filename,"yc",     x=0.0_wp,dx=1.0_wp,nx=ny,units="gridpoints")
        call nc_write_dim(filename,"time",   x=time_init,dx=1.0_wp,nx=1,units="iter",unlimited=.TRUE.)

        return

    end subroutine yelmo_write_init_ssa

    subroutine write_step_2D_ssa(tpo,dyn,filename,time)

        implicit none 
        
        type(ytopo_class), intent(IN) :: tpo 
        type(ydyn_class),  intent(IN) :: dyn 
        character(len=*),  intent(IN) :: filename 
        real(prec), intent(IN) :: time

        ! Local variables
        integer    :: ncid, n, i, j, nx, ny  
        real(prec) :: time_prev 

        nx = tpo%par%nx 
        ny = tpo%par%ny 

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_size(filename,"time",ncid)
        call nc_read(filename,"time",time_prev,start=[n],count=[1],ncid=ncid) 
        if (abs(time-time_prev).gt.1e-5) n = n+1 

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        ! == yelmo_topography ==
        call nc_write(filename,"H_ice",tpo%now%H_ice,units="m",long_name="Ice thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf",tpo%now%z_srf,units="m",long_name="Surface elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_bed",tpo%now%mask_bed,units="",long_name="Bed mask", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"dzsdt",tpo%now%dzsdt,units="m/a",long_name="Surface elevation change", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dHidt",tpo%now%dHidt,units="m/a",long_name="Ice thickness change", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"H_grnd",tpo%now%H_grnd,units="m",long_name="Ice thickness overburden", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"f_grnd",tpo%now%f_grnd,units="1",long_name="Grounded fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_grnd_acx",tpo%now%f_grnd_acx,units="1",long_name="Grounded fraction (acx)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_grnd_acy",tpo%now%f_grnd_acy,units="1",long_name="Grounded fraction (acy)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_ice",tpo%now%f_ice,units="1",long_name="Ice-covered fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"calv",tpo%now%calv,units="m/a",long_name="Calving rate", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! == yelmo_dynamics ==

        call nc_write(filename,"ssa_mask_acx",dyn%now%ssa_mask_acx,units="1",long_name="SSA mask (acx)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"ssa_mask_acy",dyn%now%ssa_mask_acy,units="1",long_name="SSA mask (acy)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"c_bed",dyn%now%c_bed,units="Pa",long_name="Dragging coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"N_eff",dyn%now%N_eff,units="Pa",long_name="Effective pressure", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"beta",dyn%now%beta,units="Pa a m^-1",long_name="Dragging coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"beta_acx",dyn%now%beta_acx,units="Pa a m^-1",long_name="Dragging coefficient (acx)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"beta_acy",dyn%now%beta_acy,units="Pa a m^-1",long_name="Dragging coefficient (acy)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dyn_visc_eff_int",dyn%now%visc_eff_int,units="Pa a",long_name="Vertically integrated effective viscosity", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"taud_acx",dyn%now%taud_acx,units="Pa",long_name="Driving stress, x-direction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taud_acy",dyn%now%taud_acy,units="Pa",long_name="Driving stress, y-direction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
!         call nc_write(filename,"ux_i_bar",dyn%now%ux_i_bar,units="m/a",long_name="Internal shear velocity (x)", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"uy_i_bar",dyn%now%uy_i_bar,units="m/a",long_name="Internal shear velocity (y)", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"uxy_i_bar",dyn%now%uxy_i_bar,units="m/a",long_name="Internal shear velocity magnitude", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"taul_int_acx",dyn%now%taul_int_acx,units="Pa m",long_name="Vertically integrated lateral boundary stress, x-direction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taul_int_acy",dyn%now%taul_int_acy,units="Pa m",long_name="Vertically integrated lateral boundary stress, y-direction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"ux_bar",dyn%now%ux_bar,units="m/a",long_name="Depth-averaged velocity (x)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy_bar",dyn%now%uy_bar,units="m/a",long_name="Depth-averaged velocity (y)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_bar",dyn%now%uxy_bar,units="m/a",long_name="Depth-averaged velocity magnitude", &
                     dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"ux_bar_diff",dyn%now%ux_bar-dyn%now%ux_bar_prev,units="m/a",long_name="Depth-averaged velocity difference (x)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy_bar_diff",dyn%now%uy_bar-dyn%now%uy_bar_prev,units="m/a",long_name="Depth-averaged velocity difference (y)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"ssa_err_acx",dyn%now%ssa_err_acx,units="1",long_name="SSA L1 error metric (x)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"ssa_err_acy",dyn%now%ssa_err_acy,units="1",long_name="SSA L1 error metric (y)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
!         call nc_write(filename,"ux",dyn%now%ux,units="m/a",long_name="Horizontal velocity (x)", &
!                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
!         call nc_write(filename,"uy",dyn%now%uy,units="m/a",long_name="Horizontal velocity (y)", &
!                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
!         call nc_write(filename,"uxy",dyn%now%uxy,units="m/a",long_name="Horizontal velocity magnitude", &
!                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
!         call nc_write(filename,"uz",dyn%now%uz,units="m/a",long_name="Vertical velocity", &
!                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

!         call nc_write(filename,"f_vbvs",dyn%now%f_vbvs,units="1",long_name="Basal to surface velocity fraction", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"f_shear_bar",mat%now%f_shear_bar,units="1",long_name="Vertically averaged shearing fraction", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine write_step_2D_ssa
    
end module yelmo_dynamics





! if (.FALSE.) then
! ! Testing exotic mixing solutions for treating the grounding line  
!             ! Set dyn1 equal to previous solution 
!             dyn1 = dyn 

!             ! Determine ssa mask for points near grounding line
!             dyn1%now%ssa_mask_acx = -1.0 
!             dyn1%now%ssa_mask_acy = -1.0  
             
!             do j = 1, ny 
!             do i = 1, nx-1 

!                 is_grz_mid = tpo%now%is_grz(i,j) .or. tpo%now%is_grz(i+1,j)
!                 if (dyn%now%ssa_mask_acx(i,j) .gt. 0.0 .and. is_grz_mid) then 
!                     dyn1%now%ssa_mask_acx(i,j) = 1.0 
!                 end if 

!             end do 
!             end do 

!             do j = 1, ny-1 
!             do i = 1, nx 

!                 is_grz_mid = tpo%now%is_grz(i,j) .or. tpo%now%is_grz(i,j+1)
!                 if (dyn%now%ssa_mask_acy(i,j) .gt. 0.0 .and. is_grz_mid) then 
!                     dyn1%now%ssa_mask_acy(i,j) = 1.0 
!                 end if 

!             end do 
!             end do 
            
!             ! Now populate dyn2 
!             dyn2 = dyn1 

!             ! Modify dyn1 parameters concerning beta 
!             dyn1%par%taud_gl_method = 1 
!             dyn1%par%beta_gl_sep    = 0     ! No subgrid grounding line treatment 
!             dyn1%par%beta_gl_scale  = 0     ! No special scaling at gl 
!             dyn1%par%beta_gl_stag   = 1     ! Upstream scaling 

!             ! Calculate driving stress 
!             call calc_driving_stress(dyn1%now%taud_acx,dyn1%now%taud_acy,tpo%now%H_ice,tpo%now%z_srf,bnd%z_bed,bnd%z_sl, &
!                      tpo%now%H_grnd,tpo%now%f_grnd,tpo%now%f_grnd_acx,tpo%now%f_grnd_acy,dyn1%par%dx,dyn1%par%taud_lim, &
!                      method=dyn1%par%taud_gl_method,beta_gl_stag=dyn1%par%beta_gl_stag)

!             call calc_ydyn_ssa(dyn1,tpo,thrm,mat,bnd)

!             ! Set dyn2 equal to previous solution 
!             !dyn2 = dyn 

!             ! Modify dyn1 parameters concerning beta 
!             dyn2%par%taud_gl_method = 1 
!             dyn2%par%beta_gl_sep    = 0     ! No subgrid grounding line treatment 
!             dyn2%par%beta_gl_scale  = 0     ! No special scaling at gl 
!             dyn2%par%beta_gl_stag   = 2     ! Downstream scaling 
            
!             ! Calculate driving stress 
!             call calc_driving_stress(dyn2%now%taud_acx,dyn2%now%taud_acy,tpo%now%H_ice,tpo%now%z_srf,bnd%z_bed,bnd%z_sl, &
!                      tpo%now%H_grnd,tpo%now%f_grnd,tpo%now%f_grnd_acx,tpo%now%f_grnd_acy,dyn2%par%dx,dyn1%par%taud_lim, &
!                      method=dyn2%par%taud_gl_method,beta_gl_stag=dyn2%par%beta_gl_stag)

!             call calc_ydyn_ssa(dyn2,tpo,thrm,mat,bnd)
            
!             ! Get weighted-average of the two solutions 

!             !dyn%now%taud_acx = tpo%now%f_grnd_acx*dyn1%now%taud_acx + (1.0-tpo%now%f_grnd_acx)*dyn2%now%taud_acx
!             !dyn%now%taud_acy = tpo%now%f_grnd_acy*dyn1%now%taud_acy + (1.0-tpo%now%f_grnd_acy)*dyn2%now%taud_acy
            
!             do j = 1, ny 
!             do i = 1, nx 
!                 if (tpo%now%f_grnd_acx(i,j) .gt. 0.0 .and. tpo%now%f_grnd_acx(i,j) .lt. 1.0) then 
!                     dyn%now%ux_b(i,j) = tpo%now%f_grnd_acx(i,j)*dyn1%now%ux_b(i,j) &
!                                         + (1.0-tpo%now%f_grnd_acx(i,j))*dyn2%now%ux_b(i,j)
!                     dyn%now%ux_b(i,j) = dyn1%now%ux_b(i,j)
!                     dyn%now%ssa_mask_acx(i,j) = -1.0 
!                 end if 

!             end do 
!             end do 

!             do j = 1, ny 
!             do i = 1, nx 
!                 if (tpo%now%f_grnd_acy(i,j) .gt. 0.0 .and. tpo%now%f_grnd_acy(i,j) .lt. 1.0) then 
!                     dyn%now%uy_b(i,j) = tpo%now%f_grnd_acy(i,j)*dyn1%now%uy_b(i,j) &
!                                         + (1.0-tpo%now%f_grnd_acy(i,j))*dyn2%now%uy_b(i,j)
!                     dyn%now%uy_b(i,j) = dyn1%now%uy_b(i,j)
!                     dyn%now%ssa_mask_acy(i,j) = -1.0 
!                 end if 

!             end do 
!             end do 
! end if 

