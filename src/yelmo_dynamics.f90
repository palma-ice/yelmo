
module yelmo_dynamics

    use nml 
    use ncio 

    use yelmo_defs
    use yelmo_tools, only : stagger_aa_acx, stagger_aa_acy, stagger_ac_aa, &
        calc_magnitude_from_staggered_ice, calc_vertical_integrated_2D, smooth_gauss_2D, &
        regularize2D !, limit_gradient

    use velocity_diva
    use velocity_hybrid_pd12 
    use velocity_sia 
    use solver_ssa_sico5
    use basal_dragging  
    use grounding_line_flux 

    ! Note: 3D arrays defined such that first index (k=1) == base, and max index (k=nk) == surface 
    
    implicit none
      
    private

    public :: ydyn_par_load, ydyn_alloc, ydyn_dealloc
    public :: calc_ydyn
    public :: calc_ydyn_neff, calc_ydyn_cfref
    
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
        real(prec) :: dt 

        ! Initialize time if necessary 
        if (dyn%par%time .gt. time) then 
            dyn%par%time = time
        end if 

        ! Get time step
        dt = time - dyn%par%time 

        ! Store previous (n-1) depth-averaged horizontal velocity components
        ! (for use with higher-order ice thickness timestepping) 
        dyn%now%ux_bar_prev = dyn%now%ux_bar 
        dyn%now%uy_bar_prev = dyn%now%uy_bar 
        
        ! ===== Calculate the horizontal velocity components =====
        ! These calculations are done assuming that the final
        ! 3D horizontal velocity fields (ux/uy) will be comprised
        ! of an internal shear contribution (ux_i/uy_i) and a plug
        ! flow contribution represented by basal velocity (ux_b/uy_b)

        select case(dyn%par%solver)

            case("fixed") 
                ! Do nothing - dynamics is fixed 

            case("hybrid")
                ! Classic SIA+SSA, including alternative SIA/SSA combinations and purely SIA or SSA

                call calc_ydyn_adhoc(dyn,tpo,mat,thrm,bnd,dt)

            case("hybrid-pd12")
                ! Variational approach of Pollard and de Conto (2012) - in progress!

                write(*,*) "Not working - to do."
                stop 

                call calc_ydyn_pd12(dyn,tpo,mat,thrm,bnd)

            case("diva") 
                ! Depth-integrated variational approximation (DIVA) - Goldberg (2011); Lipscomb et al. (2019)

                call calc_ydyn_diva(dyn,tpo,mat,thrm,bnd,dt)

            case DEFAULT 

                write(*,*) "calc_ydyn:: Error: ydyn solver not recognized." 
                write(*,*) "solver should be one of: ['fixed','hybrid','diva']"
                write(*,*) "solver = ", trim(dyn%par%solver) 
                stop 

        end select 

        ! Advance ydyn timestep 
        dyn%par%time = time

        return

    end subroutine calc_ydyn
    
    subroutine calc_ydyn_diva(dyn,tpo,mat,thrm,bnd,dt)
        ! Velocity is a steady-state solution to a given set of boundary conditions (topo, material, etc)
        ! so no time step is passed here. 

        implicit none
        
        type(ydyn_class),   intent(INOUT) :: dyn
        type(ytopo_class),  intent(IN)    :: tpo 
        type(ymat_class),   intent(IN)    :: mat
        type(ytherm_class), intent(IN)    :: thrm 
        type(ybound_class), intent(IN)    :: bnd  
        real(prec),         intent(IN)    :: dt 

        ! Local variables
        integer :: iter, n_iter
        integer :: i, j, k, nx, ny, nz_aa, nz_ac   
        real(prec), allocatable :: uxy_prev(:,:) 

        type(diva_param_class) :: diva_par 

        ! For vertical velocity calculation 
        real(prec), allocatable :: bmb(:,:)

        logical :: calc_ssa 

        type(ydyn_class) :: dyn1 
        type(ydyn_class) :: dyn2 
        logical          :: is_grz_mid 
        real(prec)       :: H_mid 

        nx    = dyn%par%nx 
        ny    = dyn%par%ny 
        nz_aa = dyn%par%nz_aa 
        nz_ac = dyn%par%nz_ac 
        
        allocate(bmb(nx,ny))
        allocate(uxy_prev(nx,ny)) 

        ! ===== Calculate general variables ==============================

        ! Store initial uxy_bar solution 
        uxy_prev = dyn%now%uxy_bar 
        
        ! Calculate driving stress 
        call calc_driving_stress(dyn%now%taud_acx,dyn%now%taud_acy,tpo%now%H_ice,tpo%now%dzsdx,tpo%now%dzsdy, &
                                                                                            dyn%par%dx,dyn%par%taud_lim)

        ! Calculate effective pressure 
        call calc_ydyn_neff(dyn,tpo,thrm,bnd)

        ! Update bed roughness coefficients cf_ref and c_bed (which are independent of velocity)
        call calc_ydyn_cfref(dyn,tpo,thrm,bnd)
        call calc_c_bed(dyn%now%c_bed,dyn%now%cf_ref,dyn%now%N_eff)

        ! ===== Calculate 3D horizontal velocity solution via DIVA algorithm ===================

        ! Define grid points with ssa active (uses beta from previous timestep)
        call set_ssa_masks(dyn%now%ssa_mask_acx,dyn%now%ssa_mask_acy,dyn%now%beta_acx,dyn%now%beta_acy, &
                           tpo%now%H_ice,tpo%now%f_grnd_acx,tpo%now%f_grnd_acy,dyn%par%ssa_beta_max,use_ssa=.TRUE.)

        ! Set diva parameters from Yelmo settings 
        diva_par%ssa_solver_opt = dyn%par%ssa_solver_opt 
        diva_par%boundaries     = dyn%par%boundaries 
        diva_par%diva_no_slip   = dyn%par%diva_no_slip 
        diva_par%beta_method    = dyn%par%beta_method 
        diva_par%beta_const     = dyn%par%beta_const 
        diva_par%beta_q         = dyn%par%beta_q 
        diva_par%beta_u0        = dyn%par%beta_u0 
        diva_par%beta_gl_scale  = dyn%par%beta_gl_scale 
        diva_par%beta_gl_stag   = dyn%par%beta_gl_stag 
        diva_par%beta_gl_f      = dyn%par%beta_gl_f 
        diva_par%H_grnd_lim     = dyn%par%H_grnd_lim 
        diva_par%beta_min       = dyn%par%beta_min 
        diva_par%ssa_vel_max    = dyn%par%ssa_vel_max 
        diva_par%ssa_iter_max   = dyn%par%ssa_iter_max 
        diva_par%ssa_iter_rel   = dyn%par%ssa_iter_rel 
        diva_par%ssa_iter_conv  = dyn%par%ssa_iter_conv 
        diva_par%ssa_write_log  = yelmo_log

        call calc_velocity_diva(dyn%now%ux,dyn%now%uy,dyn%now%ux_i,dyn%now%uy_i,dyn%now%ux_bar,dyn%now%uy_bar, &
                                dyn%now%ux_b,dyn%now%uy_b,dyn%now%duxdz,dyn%now%duydz,dyn%now%taub_acx,dyn%now%taub_acy, &
                                dyn%now%visc_eff,dyn%now%visc_eff_int,dyn%now%ssa_mask_acx,dyn%now%ssa_mask_acy, &
                                dyn%now%ssa_err_acx,dyn%now%ssa_err_acy,dyn%now%beta,dyn%now%beta_acx,dyn%now%beta_acy, &
                                dyn%now%beta_eff,dyn%now%beta_diva,dyn%now%c_bed,dyn%now%taud_acx,dyn%now%taud_acy,tpo%now%H_ice, &
                                tpo%now%H_grnd,tpo%now%f_grnd,tpo%now%f_grnd_acx,tpo%now%f_grnd_acy,mat%now%ATT, &
                                dyn%par%zeta_aa,bnd%z_sl,bnd%z_bed,dyn%par%dx,dyn%par%dy,mat%par%n_glen,diva_par)

        ! ===== Calculate the vertical velocity through continuity ============================

        if (dyn%par%use_bmb) then 
            bmb = tpo%now%bmb 
        else 
            bmb = 0.0 
        end if 

        call calc_uz_3D(dyn%now%uz,dyn%now%ux,dyn%now%uy,tpo%now%H_ice,bnd%z_bed,tpo%now%z_srf, &
                        bnd%smb,bmb,tpo%now%dHicedt,tpo%now%dzsrfdt,dyn%par%zeta_aa,dyn%par%zeta_ac,dyn%par%dx,dyn%par%dy)
        
        ! ===== Additional diagnostic variables =======================

        ! Integrate from 3D shear velocity field to get depth-averaged field
        dyn%now%ux_i_bar = calc_vertical_integrated_2D(dyn%now%ux_i,dyn%par%zeta_aa)
        dyn%now%uy_i_bar = calc_vertical_integrated_2D(dyn%now%uy_i,dyn%par%zeta_aa)
                
        ! Diagnose ice flux 
        call calc_ice_flux(dyn%now%qq_acx,dyn%now%qq_acy,dyn%now%ux_bar,dyn%now%uy_bar,tpo%now%H_ice, &
                            dyn%par%dx,dyn%par%dy)
        dyn%now%qq        = calc_magnitude_from_staggered_ice(dyn%now%qq_acx,dyn%now%qq_acy,tpo%now%H_ice)

        dyn%now%taub      = calc_magnitude_from_staggered_ice(dyn%now%taub_acx,dyn%now%taub_acy,tpo%now%H_ice)
        dyn%now%taud      = calc_magnitude_from_staggered_ice(dyn%now%taud_acx,dyn%now%taud_acy,tpo%now%H_ice)

        dyn%now%uxy_b     = calc_magnitude_from_staggered_ice(dyn%now%ux_b,dyn%now%uy_b,tpo%now%H_ice)
        dyn%now%uxy_i_bar = calc_magnitude_from_staggered_ice(dyn%now%ux_i_bar,dyn%now%uy_i_bar,tpo%now%H_ice)
        dyn%now%uxy_bar   = calc_magnitude_from_staggered_ice(dyn%now%ux_bar,dyn%now%uy_bar,tpo%now%H_ice)

        do k = 1, nz_aa
            dyn%now%uxy(:,:,k) = calc_magnitude_from_staggered_ice(dyn%now%ux(:,:,k),dyn%now%uy(:,:,k),tpo%now%H_ice)
        end do 

        ! Store surface velocities for easy access too 
        dyn%now%ux_s  = dyn%now%ux(:,:,nz_aa)
        dyn%now%uy_s  = dyn%now%uy(:,:,nz_aa)
        dyn%now%uxy_s = dyn%now%uxy(:,:,nz_aa)

        ! Determine ratio of basal to surface velocity
        dyn%now%f_vbvs = calc_vel_ratio(uxy_base=dyn%now%uxy_b,uxy_srf=dyn%now%uxy_s)

        ! Finally, determine rate of velocity change 
        if (dt .ne. 0.0_prec) then 
            dyn%now%duxydt = (dyn%now%uxy_bar - uxy_prev) / dt 
        else 
            dyn%now%duxydt = 0.0_prec 
        end if 

        return

    end subroutine calc_ydyn_diva

    subroutine calc_ydyn_adhoc(dyn,tpo,mat,thrm,bnd,dt)
        ! Velocity is a steady-state solution to a given set of boundary conditions (topo, material, etc)
        ! so no time step is passed here. 

        implicit none
        
        type(ydyn_class),   intent(INOUT) :: dyn
        type(ytopo_class),  intent(IN)    :: tpo 
        type(ymat_class),   intent(IN)    :: mat
        type(ytherm_class), intent(IN)    :: thrm 
        type(ybound_class), intent(IN)    :: bnd  
        real(prec),         intent(IN)    :: dt 

        ! Local variables
        integer :: iter, n_iter
        integer :: i, j, k, nx, ny, nz_aa, nz_ac   
        real(prec), allocatable :: H_ice_acx(:,:) 
        real(prec), allocatable :: H_ice_acy(:,:) 
        real(prec), allocatable :: ux_b_prev(:,:) 
        real(prec), allocatable :: uy_b_prev(:,:) 

        real(prec), allocatable :: uxy_prev(:,:) 

        ! For vertical velocity calculation 
        real(prec), allocatable :: bmb(:,:)

        logical :: calc_ssa 

        type(ydyn_class) :: dyn1 
        type(ydyn_class) :: dyn2 
        logical          :: is_grz_mid 
        real(prec)       :: H_mid 

        nx    = dyn%par%nx 
        ny    = dyn%par%ny 
        nz_aa = dyn%par%nz_aa 
        nz_ac = dyn%par%nz_ac 
        
        allocate(H_ice_acx(nx,ny))
        allocate(H_ice_acy(nx,ny))
        allocate(bmb(nx,ny))

        allocate(ux_b_prev(nx,ny))
        allocate(uy_b_prev(nx,ny))
        
        allocate(uxy_prev(nx,ny)) 

        ! Store initial uxy_bar solution 
        uxy_prev = dyn%now%uxy_bar 

        ! Stagger the ice thickness, aa=>ac nodes
        H_ice_acx = stagger_aa_acx(tpo%now%H_ice)
        H_ice_acy = stagger_aa_acy(tpo%now%H_ice)

        ! ===== Calculate general variables ==============================

        ! Calculate driving stress 
        call calc_driving_stress(dyn%now%taud_acx,dyn%now%taud_acy,tpo%now%H_ice,tpo%now%dzsdx,tpo%now%dzsdy, &
                                                                                            dyn%par%dx,dyn%par%taud_lim)

        
!         if (dyn%par%tau_gl_method .ne. 0) then
!             ! Additionally treat the driving stress at the grounding line

!             call calc_driving_stress_gl(dyn%now%taud_acx,dyn%now%taud_acy,tpo%now%H_ice,tpo%now%z_srf,bnd%z_bed,bnd%z_sl, &
!                     tpo%now%H_grnd,tpo%now%f_grnd,tpo%now%f_grnd_acx,tpo%now%f_grnd_acy,dyn%par%dx, &
!                     method=dyn%par%taud_gl_method,beta_gl_stag=dyn%par%beta_gl_stag)
        
!         end if 
        
        ! Calculate effective pressure 
        call calc_ydyn_neff(dyn,tpo,thrm,bnd)

        ! Update bed roughness coefficients cf_ref and c_bed (which are independent of velocity)
        call calc_ydyn_cfref(dyn,tpo,thrm,bnd)
        call calc_c_bed(dyn%now%c_bed,dyn%now%cf_ref,dyn%now%N_eff)

        ! ===== Calculate shear (ie, SIA) velocity solution ===========

        select case(trim(dyn%par%sia_solver))

            case("vel") 
                ! Use classic-style SIA solver (solve directly for velocity)

                ! Calculate diffusivity constant on ab-nodes
!                 if (yelmo_use_omp) then
!                     ! ajr: do not use as it is, really slow: 
!                     call calc_dd_ab_3D_omp(dyn%now%dd_ab,tpo%now%H_ice,dyn%now%taud_acx,dyn%now%taud_acy, &
!                                             mat%now%ATT,dyn%par%zeta_aa,dyn%par%dx,mat%par%n_glen,rho_ice,g)
!                 else 
                    call calc_dd_ab_3D_serial(dyn%now%dd_ab,tpo%now%H_ice,dyn%now%taud_acx,dyn%now%taud_acy, &
                                            mat%now%ATT,dyn%par%zeta_aa,dyn%par%dx,mat%par%n_glen,rho_ice,g)
!                 end if 

                dyn%now%dd_ab_bar = calc_vertical_integrated_2D(dyn%now%dd_ab,dyn%par%zeta_aa)
                

                ! Calculate the 3D horizontal shear velocity fields
                call calc_uxy_sia_3D(dyn%now%ux_i,dyn%now%uy_i,dyn%now%dd_ab, &
                                        dyn%now%taud_acx,dyn%now%taud_acy)
                
                ! Calculate the depth-averaged horizontal shear velocity fields too
                call calc_uxy_sia_2D(dyn%now%ux_i_bar,dyn%now%uy_i_bar,dyn%now%dd_ab, &
                                        dyn%now%taud_acx,dyn%now%taud_acy,dyn%par%zeta_aa)

                ! Or, simply integrate from 3D velocity field to get depth-averaged field
!                 dyn%now%ux_i_bar = calc_vertical_integrated_2D(dyn%now%ux_i,dyn%par%zeta_aa)
!                 dyn%now%uy_i_bar = calc_vertical_integrated_2D(dyn%now%uy_i,dyn%par%zeta_aa)
                
                ! Set terms from shear solver to zero that are not calculated here
                dyn%now%duxdz     = 0.0 
                dyn%now%duydz     = 0.0 
                dyn%now%dd_ab     = 0.0 
                dyn%now%duxdz_bar = 0.0 
                dyn%now%duydz_bar = 0.0 

            case("shear")
                ! Use 3D shear solver for SIA calculations (solve for shear and integrate to get velocity)
                ! ie, following Pollard and de Conto (2012) 

                ! Ensure PD12 stretching terms are zero  
                dyn%now%lhs_x          = 0.0 
                dyn%now%lhs_y          = 0.0 
                dyn%now%sigma_horiz_sq = 0.0 
                
                ! Calculate the 3D vertical shear fields 
                call calc_shear_3D(dyn%now%duxdz,dyn%now%duydz,dyn%now%dd_ab, &
                                   dyn%now%taud_acx,dyn%now%taud_acy,mat%now%ATT, &
                                   dyn%now%lhs_x,dyn%now%lhs_y,dyn%now%sigma_horiz_sq, &
                                   dyn%par%zeta_aa,mat%par%n_glen,boundaries=dyn%par%boundaries)

                ! Calculate the vertically averaged shear (for potential use in ssa viscosity) 
                dyn%now%duxdz_bar = calc_vertical_integrated_2D(dyn%now%duxdz,dyn%par%zeta_aa)
                dyn%now%duydz_bar = calc_vertical_integrated_2D(dyn%now%duydz,dyn%par%zeta_aa)
                
                ! Integrate vertically to get 3D horizontal shear velocity fields
                dyn%now%ux_i = calc_vertical_integrated_3D_ice(dyn%now%duxdz,H_ice_acx,dyn%par%zeta_aa)
                dyn%now%uy_i = calc_vertical_integrated_3D_ice(dyn%now%duydz,H_ice_acy,dyn%par%zeta_aa)

                ! Calculate depth-averaged horizontal shear velocity fields
                dyn%now%ux_i_bar  = calc_vertical_integrated_2D(dyn%now%ux_i,dyn%par%zeta_aa) 
                dyn%now%uy_i_bar  = calc_vertical_integrated_2D(dyn%now%uy_i,dyn%par%zeta_aa) 

                ! Diagnostic diffusivity field not available
                dyn%now%dd_ab     = 0.0_prec
                dyn%now%dd_ab_bar = 0.0_prec 

            case("none")

                ! Set shear velocity terms to zero 
                dyn%now%ux_i      = 0.0 
                dyn%now%uy_i      = 0.0 
                dyn%now%ux_i_bar  = 0.0 
                dyn%now%uy_i_bar  = 0.0

                ! Set terms from shear solver to zero that are not calculated here
                dyn%now%duxdz     = 0.0 
                dyn%now%duydz     = 0.0 
                dyn%now%dd_ab     = 0.0 
                dyn%now%duxdz_bar = 0.0 
                dyn%now%duydz_bar = 0.0 

                ! Ensure PD12 stretching terms are zero  
                dyn%now%lhs_x          = 0.0 
                dyn%now%lhs_y          = 0.0 
                dyn%now%sigma_horiz_sq = 0.0 
            
            case("fixed") 
                ! Pass - do nothing with SIA, use whatever solution is available in the fields

            case DEFAULT 

                write(*,*) "calc_ydyn_adhoc:: Error: sia solver not recognized: ", trim(dyn%par%sia_solver)
                stop 

        end select 

        ! ===== Calculate sliding velocity solution ===================

        ! Define grid points with ssa active (uses beta from previous timestep)
        call set_ssa_masks(dyn%now%ssa_mask_acx,dyn%now%ssa_mask_acy,dyn%now%beta_acx,dyn%now%beta_acy, &
                           tpo%now%H_ice,tpo%now%f_grnd_acx,tpo%now%f_grnd_acy,dyn%par%ssa_beta_max,dyn%par%use_ssa)

        ! Determine whether SSA solver should be called 
        calc_ssa = .FALSE. 
        if (dyn%par%use_ssa .and. maxval(dyn%now%ssa_mask_acx+dyn%now%ssa_mask_acy) .gt. 0) calc_ssa = .TRUE.    
            
        if (calc_ssa) then
            ! == Iterate over strain rate, viscosity and ssa velocity solutions until convergence ==
            ! Note: ssa solution is defined for ux_b/uy_b fields here, not ux_bar/uy_bar as in PD12
              
            call calc_ydyn_ssa(dyn,tpo,thrm,mat,bnd)
            
        else 
            ! No ssa calculations performed, set basal velocity fields to zeor 

            dyn%now%ux_b = 0.0 
            dyn%now%uy_b = 0.0 

        end if 

        ! ===== Combine sliding and shear into hybrid fields ==========

        select case(dyn%par%mix_method)
            ! Determine how to mix shearing (ux_i/uy_i) and sliding (ux_b/uy_b)
            ! Generally purely floating ice is only given by the SSA solution;
            ! Grounded or partially grounded ice is given by the hybrid solution 

            case(-2)
                ! Purely sia model
                ! (correspondingly, floating ice is killed in yelmo_topography)
                
                if (dyn%par%cb_sia .gt. 0.0) then 
                    ! Calculate basal velocity from Weertman sliding law (Greve 1997)
                    
                    call calc_uxy_b_sia(dyn%now%ux_b,dyn%now%uy_b,tpo%now%H_ice,tpo%now%dzsdx,tpo%now%dzsdy, &
                                thrm%now%f_pmp,dyn%par%zeta_aa,dyn%par%dx,dyn%par%cb_sia,rho_ice,g)
                
                else 
                    ! Otherwise no basal sliding in SIA-only mode
                
                    dyn%now%ux_b   = 0.0_prec 
                    dyn%now%uy_b   = 0.0_prec
                    
                end if 

                ! SIA solution everywhere, potentially with additional parameterized basal sliding
                dyn%now%ux_bar = dyn%now%ux_i_bar + dyn%now%ux_b 
                dyn%now%uy_bar = dyn%now%uy_i_bar + dyn%now%uy_b 

            case(-1)
                ! Purely ssa model

                dyn%now%ux_i     = 0.0_prec 
                dyn%now%uy_i     = 0.0_prec 
                dyn%now%ux_i_bar = 0.0_prec 
                dyn%now%uy_i_bar = 0.0_prec 
                
                dyn%now%ux_bar = dyn%now%ux_b 
                dyn%now%uy_bar = dyn%now%uy_b 

            case(0)
                ! Binary mixing (either shearing or sliding)

                ! TO DO 

                write(*,*) "mix_method=0 not implemented yet."
                stop 

            case(1)
                ! Shear when not sliding, otherwise shear+sliding 
                ! (ie, sia or sia+ssa)

                ! Hybrid solution everywhere 
                dyn%now%ux_bar = dyn%now%ux_i_bar + dyn%now%ux_b 
                dyn%now%uy_bar = dyn%now%uy_i_bar + dyn%now%uy_b 

            case(2)
                ! Weighted mixing 

                ! TO DO 

                write(*,*) "mix_method=2 not implemented yet."
                stop 

            case DEFAULT

                write(*,*) "mix_method not recognized."
                write(*,*) "mix_method = ", dyn%par%mix_method
                stop 

        end select

        ! ===== Calculate 3D velocity fields ========================== 

        ! Fill in horizontal velocity as the sum of shear and basal velocity 
        do k = 1, nz_aa 
            dyn%now%ux(:,:,k)  = dyn%now%ux_i(:,:,k) + dyn%now%ux_b 
            dyn%now%uy(:,:,k)  = dyn%now%uy_i(:,:,k) + dyn%now%uy_b 
        end do 

        ! Calculate the vertical velocity through continuity

        if (dyn%par%use_bmb) then 
            bmb = tpo%now%bmb 
        else 
            bmb = 0.0 
        end if 

        call calc_uz_3D(dyn%now%uz,dyn%now%ux,dyn%now%uy,tpo%now%H_ice,bnd%z_bed,tpo%now%z_srf, &
                        bnd%smb,bmb,tpo%now%dHicedt,tpo%now%dzsrfdt,dyn%par%zeta_aa,dyn%par%zeta_ac,dyn%par%dx,dyn%par%dy)
        
        ! ===== Additional diagnostic variables =======================

        ! 3b. Diagnose the shear reduction terms of PD12 
        call calc_shear_reduction(dyn%now%lhs_x,dyn%now%lhs_y,dyn%now%ux_b,dyn%now%uy_b,dyn%now%visc_eff_int,dyn%par%dx)
        
        ! 3c. Diagnose the effective horizontal stress squared
        dyn%now%sigma_horiz_sq = calc_stress_eff_horizontal_squared(dyn%now%ux_bar,dyn%now%uy_bar, &
                                            mat%now%ATT_bar,dyn%par%dx,dyn%par%dy,mat%par%n_glen)

        ! Calculate basal stress (input to thermodynamics)
        call calc_basal_stress(dyn%now%taub_acx,dyn%now%taub_acy,dyn%now%beta_acx,dyn%now%beta_acy, &
                               dyn%now%ux_b,dyn%now%uy_b)

        ! Diagnose ice flux 
        call calc_ice_flux(dyn%now%qq_acx,dyn%now%qq_acy,dyn%now%ux_bar,dyn%now%uy_bar,tpo%now%H_ice, &
                            dyn%par%dx,dyn%par%dy)
        dyn%now%qq        = calc_magnitude_from_staggered_ice(dyn%now%qq_acx,dyn%now%qq_acy,tpo%now%H_ice)

        dyn%now%taub      = calc_magnitude_from_staggered_ice(dyn%now%taub_acx,dyn%now%taub_acy,tpo%now%H_ice)
        dyn%now%taud      = calc_magnitude_from_staggered_ice(dyn%now%taud_acx,dyn%now%taud_acy,tpo%now%H_ice)

        dyn%now%lhs_xy    = calc_magnitude_from_staggered_ice(dyn%now%lhs_x,dyn%now%lhs_y,tpo%now%H_ice)
        dyn%now%uxy_b     = calc_magnitude_from_staggered_ice(dyn%now%ux_b,dyn%now%uy_b,tpo%now%H_ice)
        dyn%now%uxy_i_bar = calc_magnitude_from_staggered_ice(dyn%now%ux_i_bar,dyn%now%uy_i_bar,tpo%now%H_ice)
        dyn%now%uxy_bar   = calc_magnitude_from_staggered_ice(dyn%now%ux_bar,dyn%now%uy_bar,tpo%now%H_ice)

        do k = 1, nz_aa
            dyn%now%uxy(:,:,k) = calc_magnitude_from_staggered_ice(dyn%now%ux(:,:,k),dyn%now%uy(:,:,k),tpo%now%H_ice)
        end do 

        ! Store surface velocities for easy access too 
        dyn%now%ux_s  = dyn%now%ux(:,:,nz_aa)
        dyn%now%uy_s  = dyn%now%uy(:,:,nz_aa)
        dyn%now%uxy_s = dyn%now%uxy(:,:,nz_aa)

        ! Determine ratio of basal to surface velocity
        dyn%now%f_vbvs = calc_vel_ratio(uxy_base=dyn%now%uxy_b,uxy_srf=dyn%now%uxy_s)

        ! Finally, determine rate of velocity change 
        if (dt .ne. 0.0_prec) then 
            dyn%now%duxydt = (dyn%now%uxy_bar - uxy_prev) / dt 
        else 
            dyn%now%duxydt = 0.0_prec 
        end if 

        return

    end subroutine calc_ydyn_adhoc

    subroutine calc_ydyn_pd12(dyn,tpo,mat,thrm,bnd)
        ! Velocity is a steady-state solution to a given set of boundary conditions (topo, material, etc)
        ! so no time step is passed here. 

        implicit none
        
        type(ydyn_class),   intent(INOUT) :: dyn
        type(ytopo_class),  intent(IN)    :: tpo 
        type(ymat_class),   intent(IN)    :: mat
        type(ytherm_class), intent(IN)    :: thrm 
        type(ybound_class), intent(IN)    :: bnd  

        ! Local variables
        integer :: iter, n_iter
        integer :: k, nx, ny, nz_aa, nz_ac    
        real(prec), allocatable :: H_ice_acx(:,:) 
        real(prec), allocatable :: H_ice_acy(:,:) 

        real(prec), allocatable :: ux_bar_old(:,:) 
        real(prec), allocatable :: uy_bar_old(:,:) 
        real(prec), allocatable :: ux_bar_prev(:,:) 
        real(prec), allocatable :: uy_bar_prev(:,:) 
        
        ! For vertical velocity calculation 
        real(prec), allocatable :: bmb(:,:)

        real(prec) :: L2_norm 

        logical :: calc_ssa 
        logical :: is_converged
        logical :: write_ssa_diagnostics
        
        nx    = dyn%par%nx 
        ny    = dyn%par%ny 
        nz_aa = dyn%par%nz_aa 
        nz_ac = dyn%par%nz_ac 

        allocate(H_ice_acx(nx,ny))
        allocate(H_ice_acy(nx,ny))
        allocate(bmb(nx,ny))

        ! Stagger the ice thickness, Aa=>Ab nodes
        H_ice_acx = stagger_aa_acx(tpo%now%H_ice)
        H_ice_acy = stagger_aa_acy(tpo%now%H_ice)

        ! Update bed roughness coefficients cf_ref and c_bed (which are independent of velocity)
        call calc_ydyn_cfref(dyn,tpo,thrm,bnd)
        call calc_c_bed(dyn%now%c_bed,dyn%now%cf_ref,dyn%now%N_eff)

        ! Calculate driving stress
        call calc_driving_stress(dyn%now%taud_acx,dyn%now%taud_acy,tpo%now%H_ice,tpo%now%dzsdx,tpo%now%dzsdy, &
                                                                                            dyn%par%dx,dyn%par%taud_lim)
    
!         ! Additionally calculate driving stress at the grounding line
!         call calc_driving_stress_gl(dyn%now%taud_acx,dyn%now%taud_acy,tpo%now%H_ice,tpo%now%z_srf,bnd%z_bed,bnd%z_sl, &
!                  tpo%now%H_grnd,tpo%now%f_grnd,tpo%now%f_grnd_acx,tpo%now%f_grnd_acy,dyn%par%dx, &
!                  method=dyn%par%taud_gl_method,beta_gl_stag=dyn%par%beta_gl_stag)

        ! Calculate 2D diffusivity too (for timestepping and diagnostics)
!         call calc_diffusivity_2D(dyn%now%dd_ab_bar,tpo%now%H_ice,tpo%now%dzsdx,tpo%now%dzsdy, &
!                                  mat%now%ATT,dyn%par%zeta_aa,dyn%par%dx,mat%par%n_glen,rho_ice,g)
!         dyn%now%dd_ab_bar  = calc_vertical_integrated_2D(dyn%now%dd_ab,dyn%par%zeta_aa) 
        
        ! == Iterate over strain rate, viscosity and velocity solutions until convergence ==

        ! Initially set the mixing terms to zero 
        dyn%now%lhs_x          = 0.0 
        dyn%now%lhs_y          = 0.0 
        dyn%now%sigma_horiz_sq = 0.0 
        dyn%now%duxdz_bar      = 0.0 
        dyn%now%duydz_bar      = 0.0 

        is_converged  = .FALSE. 

        ! Store old solution (from previous time step) to be able to apply relaxation 
        ! to avoid too fast propagation of waves
        ux_bar_old = dyn%now%ux_bar
        uy_bar_old = dyn%now%uy_bar
        
        write_ssa_diagnostics = .FALSE. 

!         if (tpo%now%f_grnd(18,3) .gt. 0.0) then 
!             write_ssa_diagnostics = .TRUE.

!             call yelmo_write_init_ssa("yelmo_ssa.nc",time_init=1.0) 
!         end if 

        do iter = 1, dyn%par%ssa_iter_max

            ! Store previous solution 
            ux_bar_prev = dyn%now%ux_bar 
            uy_bar_prev = dyn%now%uy_bar 
            
            ! 1. Calculate effective stress due to stretching terms, and driving stress reduction

            dyn%now%sigma_horiz_sq = calc_stress_eff_horizontal_squared(dyn%now%ux_bar,dyn%now%uy_bar, &
                                            mat%now%ATT_bar,dyn%par%dx,dyn%par%dy,mat%par%n_glen)

            call calc_shear_reduction(dyn%now%lhs_x,dyn%now%lhs_y,dyn%now%ux_b,dyn%now%uy_b,dyn%now%visc_eff_int,dyn%par%dx)
            
            ! 2. Calculate the vertical shear fields 
            ! (accounting for effective stress with stretching and reduced driving stress)
            
            call calc_shear_3D(dyn%now%duxdz,dyn%now%duydz,dyn%now%dd_ab, &
                               dyn%now%taud_acx,dyn%now%taud_acy,mat%now%ATT, &
                               dyn%now%lhs_x,dyn%now%lhs_y,dyn%now%sigma_horiz_sq, &
                               dyn%par%zeta_aa,mat%par%n_glen,boundaries=dyn%par%boundaries)

            ! Calculate the vertically averaged shear for use in ssa viscosity 
            dyn%now%duxdz_bar = calc_vertical_integrated_2D(dyn%now%duxdz,dyn%par%zeta_aa)
            dyn%now%duydz_bar = calc_vertical_integrated_2D(dyn%now%duydz,dyn%par%zeta_aa)
            
            ! 3. Calculate shear velocity values ux_i/uy_i 

            ! Calculate the 3D shear velocity field
            dyn%now%ux_i = calc_vertical_integrated_3D_ice(dyn%now%duxdz,H_ice_acx,dyn%par%zeta_aa)
            dyn%now%uy_i = calc_vertical_integrated_3D_ice(dyn%now%duydz,H_ice_acy,dyn%par%zeta_aa)
            
            ! Calculate the vertically averaged shear velocity field 
            dyn%now%ux_i_bar  = calc_vertical_integrated_2D(dyn%now%ux_i,dyn%par%zeta_aa) 
            dyn%now%uy_i_bar  = calc_vertical_integrated_2D(dyn%now%uy_i,dyn%par%zeta_aa) 
            
!             if (iter .eq. 1) then 
!                 ! Set the hybrid solution equal to the shear solution initially 

!                 dyn%now%ux_bar = dyn%now%ux_i_bar 
!                 dyn%now%uy_bar = dyn%now%uy_i_bar 

!             end if 

            ! 4. Calculate basal velocity from previous u_bar and ux_i_bar
            
            call calc_vel_basal(dyn%now%ux_b,dyn%now%uy_b,dyn%now%ux_bar,dyn%now%uy_bar, &
                                dyn%now%ux_i_bar,dyn%now%uy_i_bar)
            
            ! 5. Calculate effective viscosity, including shear terms
            
            ! ---------------------------------------------------------------------
            ! Stable viscosity solutions for SSA solver:

!             dyn%now%visc_eff_int = 1e10 

            dyn%now%visc_eff_int = calc_visc_eff_2D(dyn%now%ux_bar,dyn%now%uy_bar,dyn%now%duxdz_bar*0.0,dyn%now%duydz_bar*0.0, &
                                                    tpo%now%H_ice,mat%now%ATT,dyn%par%zeta_aa,dyn%par%dx,dyn%par%dy,mat%par%n_glen)
            
            ! Ensure viscosity is relatively smooth
!             call regularize2D(dyn%now%visc_eff_int,tpo%now%H_ice,tpo%par%dx)

            ! ---------------------------------------------------------------------
            
            ! 6. Calculate basal drag coefficient beta (beta, beta_acx, beta_acy) 

!             call calc_ydyn_beta(dyn,tpo,mat,bnd)
            
            write(*,*) "Need to update interface to basal_dragging::calc_beta routine here!"
            stop 

            ! 7. Calculate SSA solution if needed

            ! Define grid points with ssa active
            call set_ssa_masks(dyn%now%ssa_mask_acx,dyn%now%ssa_mask_acy,dyn%now%beta_acx,dyn%now%beta_acy, &
                               tpo%now%H_ice,tpo%now%f_grnd_acx,tpo%now%f_grnd_acy,dyn%par%ssa_beta_max,dyn%par%use_ssa)

            ! Determine whether SSA solver should be called 
            calc_ssa = .FALSE. 
            if (dyn%par%use_ssa .and. maxval(dyn%now%ssa_mask_acx+dyn%now%ssa_mask_acy) .gt. 0) calc_ssa = .TRUE.    

            if (calc_ssa) then
                ! Call ssa solver to determine ux_bar/uy_bar, where ssa_mask_acx/y are > 0
                
                call calc_vxy_ssa_matrix(dyn%now%ux_bar,dyn%now%uy_bar,L2_norm,dyn%now%beta_acx,dyn%now%beta_acy,dyn%now%visc_eff_int, &
                                     dyn%now%ssa_mask_acx,dyn%now%ssa_mask_acy,tpo%now%H_ice,dyn%now%taud_acx, &
                                     dyn%now%taud_acy,tpo%now%H_grnd,bnd%z_sl,bnd%z_bed,dyn%par%dx,dyn%par%dy, &
                                     dyn%par%ssa_vel_max,dyn%par%boundaries,dyn%par%ssa_solver_opt)

            end if 
             
            ! 7. Ensure that shear solution is applied where no SSA is calculated 

            where (dyn%now%ssa_mask_acx .eq. 0) dyn%now%ux_bar = dyn%now%ux_i_bar 
            where (dyn%now%ssa_mask_acy .eq. 0) dyn%now%uy_bar = dyn%now%uy_i_bar  
            
            ! Apply relaxation to keep things stable
            call relax_ssa(dyn%now%ux_bar,dyn%now%uy_bar,ux_bar_prev,uy_bar_prev,rel=dyn%par%ssa_iter_rel)

!             ! Check convergence of ssa solution 
!             is_converged = check_vel_convergence(dyn%now%ux_bar,dyn%now%uy_bar,ux_bar_prev,uy_bar_prev, &
!                                         dyn%par%ssa_iter_conv,iter,dyn%par%ssa_iter_max,yelmo_log)

            if (write_ssa_diagnostics) then  
                call write_step_2D_ssa(tpo,dyn,"yelmo_ssa.nc",ux_bar_prev,uy_bar_prev,time=real(iter,prec))    
            end if 

            ! Exit iterations if ssa solution has converged
            if (is_converged) exit 

        end do 
        ! == END iterations ==

!         if (write_ssa_diagnostics) then 
!             stop 
!         end if 

        ! ===== Calculate 3D velocity fields ====== 

        ! Re-calculate basal velocity from current solution u_bar and ux_i_bar
            
        call calc_vel_basal(dyn%now%ux_b,dyn%now%uy_b,dyn%now%ux_bar,dyn%now%uy_bar, &
                            dyn%now%ux_i_bar,dyn%now%uy_i_bar)
        
        ! Fill in horizontal velocity as the sum of shear and basal velocity 
        do k = 1, nz_aa 
            dyn%now%ux(:,:,k)  = dyn%now%ux_i(:,:,k) + dyn%now%ux_b 
            dyn%now%uy(:,:,k)  = dyn%now%uy_i(:,:,k) + dyn%now%uy_b 
        end do 

        ! Calculate the vertical velocity through continuity

        if (dyn%par%use_bmb) then 
            bmb = tpo%now%bmb 
        else 
            bmb = 0.0 
        end if 

        call calc_uz_3D(dyn%now%uz,dyn%now%ux,dyn%now%uy,tpo%now%H_ice,bnd%z_bed,tpo%now%z_srf, &
                        bnd%smb,bmb,tpo%now%dHicedt,tpo%now%dzsrfdt,dyn%par%zeta_aa,dyn%par%zeta_ac,dyn%par%dx,dyn%par%dy)
        
        ! ===== Additional diagnostic variables ==========

        ! Calculate basal stress (input to thermodynamics)
        call calc_basal_stress(dyn%now%taub_acx,dyn%now%taub_acy,dyn%now%beta_acx,dyn%now%beta_acy, &
                               dyn%now%ux_b,dyn%now%uy_b)

        ! Diagnose ice flux 
        call calc_ice_flux(dyn%now%qq_acx,dyn%now%qq_acy,dyn%now%ux_bar,dyn%now%uy_bar,tpo%now%H_ice, &
                            dyn%par%dx,dyn%par%dy)
        dyn%now%qq        = calc_magnitude_from_staggered_ice(dyn%now%qq_acx,dyn%now%qq_acy,tpo%now%H_ice)

        dyn%now%taub      = calc_magnitude_from_staggered_ice(dyn%now%taub_acx,dyn%now%taub_acy,tpo%now%H_ice)
        dyn%now%taud      = calc_magnitude_from_staggered_ice(dyn%now%taud_acx,dyn%now%taud_acy,tpo%now%H_ice)

        dyn%now%lhs_xy    = calc_magnitude_from_staggered_ice(dyn%now%lhs_x,dyn%now%lhs_y,tpo%now%H_ice)
        dyn%now%uxy_b     = calc_magnitude_from_staggered_ice(dyn%now%ux_b,dyn%now%uy_b,tpo%now%H_ice)
        dyn%now%uxy_i_bar = calc_magnitude_from_staggered_ice(dyn%now%ux_i_bar,dyn%now%uy_i_bar,tpo%now%H_ice)
        dyn%now%uxy_bar   = calc_magnitude_from_staggered_ice(dyn%now%ux_bar,dyn%now%uy_bar,tpo%now%H_ice)

        do k = 1, nz_aa
            dyn%now%uxy(:,:,k) = calc_magnitude_from_staggered_ice(dyn%now%ux(:,:,k),dyn%now%uy(:,:,k),tpo%now%H_ice)
        end do 

        ! Store surface velocities for easy access too 
        dyn%now%ux_s  = dyn%now%ux(:,:,nz_aa)
        dyn%now%uy_s  = dyn%now%uy(:,:,nz_aa)
        dyn%now%uxy_s = dyn%now%uxy(:,:,nz_aa)

        ! Determine ratio of basal to surface velocity
        dyn%now%f_vbvs = calc_vel_ratio(uxy_base=dyn%now%uxy_b,uxy_srf=dyn%now%uxy_s)

        return

    end subroutine calc_ydyn_pd12

    subroutine calc_ydyn_ssa(dyn,tpo,thrm,mat,bnd)
        ! Calculate the ssa solution via a linearized Picard iteration
        ! over beta, visc and velocity

        implicit none
        
        type(ydyn_class),   intent(INOUT) :: dyn
        type(ytopo_class),  intent(IN)    :: tpo 
        type(ytherm_class), intent(IN)    :: thrm 
        type(ymat_class),   intent(IN)    :: mat 
        type(ybound_class), intent(IN)    :: bnd   

        ! Local variables
        integer :: iter, i, j, nx, ny
        real(prec) :: H_mid   
        real(prec), allocatable :: ux_b_prev(:,:) 
        real(prec), allocatable :: uy_b_prev(:,:) 
        integer,    allocatable :: ssa_mask_acx(:,:) 
        integer,    allocatable :: ssa_mask_acy(:,:) 

        real(prec), allocatable :: beta_acx_prev(:,:) 
        real(prec), allocatable :: beta_acy_prev(:,:) 

        real(prec) :: L2_norm 

        logical :: is_converged
        logical :: write_ssa_diagnostics

        is_converged          = .FALSE. 
        write_ssa_diagnostics = .FALSE. 

        nx    = dyn%par%nx 
        ny    = dyn%par%ny

        allocate(ux_b_prev(nx,ny))
        allocate(uy_b_prev(nx,ny))
        
        allocate(ssa_mask_acx(nx,ny))
        allocate(ssa_mask_acy(nx,ny))

        allocate(beta_acx_prev(nx,ny))
        allocate(beta_acy_prev(nx,ny))
        
        beta_acx_prev = dyn%now%beta_acx 
        beta_acy_prev = dyn%now%beta_acy 

!             if (tpo%now%f_grnd(18,3) .gt. 0.0) then 
!                 write_ssa_diagnostics = .TRUE.

!                 call yelmo_write_init_ssa("yelmo_ssa.nc",time_init=1.0) 
!             end if 
        
        if (write_ssa_diagnostics) then 
            call yelmo_write_init_ssa("yelmo_ssa.nc",nx,ny,time_init=1.0_prec)
        end if 

        ! Store original ssa mask 
        ssa_mask_acx = dyn%now%ssa_mask_acx
        ssa_mask_acy = dyn%now%ssa_mask_acy
        
        ! Initially set error very high 
        dyn%now%ssa_err_acx = 1.0_prec 
        dyn%now%ssa_err_acy = 1.0_prec 

        do iter = 1, dyn%par%ssa_iter_max

            ! Store previous solution 
            ux_b_prev = dyn%now%ux_b 
            uy_b_prev = dyn%now%uy_b 

            !   1. Calculate basal drag coefficient beta (beta, beta_acx, beta_acy) 

!             call calc_ydyn_beta(dyn,tpo,mat,bnd)

            call calc_beta(dyn%now%beta,dyn%now%c_bed,dyn%now%ux_b,dyn%now%uy_b,tpo%now%H_ice,tpo%now%H_grnd, &
                            tpo%now%f_grnd,bnd%z_bed,bnd%z_sl,dyn%par%beta_method, &
                                dyn%par%beta_const,dyn%par%beta_q,dyn%par%beta_u0,dyn%par%beta_gl_scale, &
                                dyn%par%beta_gl_f,dyn%par%H_grnd_lim,dyn%par%beta_min,dyn%par%boundaries)

            ! Stagger beta
            call stagger_beta(dyn%now%beta_acx,dyn%now%beta_acy,dyn%now%beta,tpo%now%f_grnd, &
                            tpo%now%f_grnd_acx,tpo%now%f_grnd_acy,dyn%par%beta_gl_stag,dyn%par%boundaries)

            !   2. Calculate effective viscosity
            
            ! Note: disable shear contribution to viscosity for this solver, for mixed terms use hybrid-pd12 option.
            ! Note: Here visc_eff_int is calculated using ux_b and uy_b (ssa velocity), not ux_bar/uy_bar as in hybrid-pd12. 
            dyn%now%visc_eff_int = calc_visc_eff_2D(dyn%now%ux_b,dyn%now%uy_b,dyn%now%duxdz_bar*0.0,dyn%now%duydz_bar*0.0, &
                                                    tpo%now%H_ice,mat%now%ATT,dyn%par%zeta_aa,dyn%par%dx,dyn%par%dy,mat%par%n_glen)
            
            ! Ensure viscosity is relatively smooth
!             call regularize2D(dyn%now%visc_eff_int,tpo%now%H_ice,tpo%par%dx)

            !   X. Prescribe grounding-line flux 
if (.FALSE.) then
            ! Testing prescribed grounding-line flux/vel - experimental!!!

            ! Calculate the analytical grounding-line flux 
            call calc_grounding_line_flux(dyn%now%qq_gl_acx,dyn%now%qq_gl_acy,tpo%now%H_ice,mat%now%ATT_bar, &
                        dyn%now%c_bed,dyn%now%ux_b,dyn%now%uy_b,tpo%now%f_grnd,tpo%now%f_grnd_acx,tpo%now%f_grnd_acy, &
                        mat%par%n_glen,dyn%par%beta_q,Q0=0.61_prec,f_drag=0.6_prec,gl_flux_method="power")

            ! Where qq_gl is present, prescribe velocity and set mask to -1

            ! Restore original ssa mask (without grounding line flags)
            dyn%now%ssa_mask_acx = ssa_mask_acx
            dyn%now%ssa_mask_acy = ssa_mask_acy
            
            write(*,*) "glf"

            ! acx nodes 
            do j = 1, ny 
            do i = 1, nx-1

                H_mid = 0.5*(tpo%now%H_ice(i,j)+tpo%now%H_ice(i+1,j))
                
                if (dyn%now%qq_gl_acx(i,j) .ne. 0.0 .and. H_mid .gt. 0.0) then 
                    ! Prescribe velocity at this point 

                    if (j == 3) then 
                        write(*,*) "glf", i, dyn%now%ux_b(i,j), dyn%now%qq_gl_acx(i,j) / H_mid
                    end if 
                    
!                     dyn%now%ux_b(i,j) = dyn%now%qq_gl_acx(i,j) / H_mid 
!                     dyn%now%ssa_mask_acx(i,j) = -1

                end if 

            end do 
            end do 

            ! acy nodes 
            do j = 1, ny-1 
            do i = 1, nx

                H_mid = 0.5*(tpo%now%H_ice(i,j)+tpo%now%H_ice(i,j+1))
                
                if (dyn%now%qq_gl_acy(i,j) .ne. 0.0 .and. H_mid .gt. 0.0) then 
                    ! Prescribe velocity at this point 

                    dyn%now%uy_b(i,j) = dyn%now%qq_gl_acy(i,j) / H_mid 
                    dyn%now%ssa_mask_acy(i,j) = -1
                    
                end if 

            end do 
            end do
end if 

            !   3. Calculate SSA solution

if (.TRUE.) then 
            if (iter .gt. 1) then
                ! Update ssa mask based on convergence with previous step to reduce calls 
                call update_ssa_mask_convergence(dyn%now%ssa_mask_acx,dyn%now%ssa_mask_acy, &
                                                dyn%now%ssa_err_acx,dyn%now%ssa_err_acy,err_lim=real(1e-3,prec)) 
            end if 
end if 

            ! Call ssa solver to determine ux_b/uy_b, where ssa_mask_acx/y are > 0
            call calc_vxy_ssa_matrix(dyn%now%ux_b,dyn%now%uy_b,L2_norm,dyn%now%beta_acx,dyn%now%beta_acy,dyn%now%visc_eff_int, &
                                     dyn%now%ssa_mask_acx,dyn%now%ssa_mask_acy,tpo%now%H_ice,dyn%now%taud_acx, &
                                     dyn%now%taud_acy,tpo%now%H_grnd,bnd%z_sl,bnd%z_bed,dyn%par%dx,dyn%par%dy, &
                                     dyn%par%ssa_vel_max,dyn%par%boundaries,dyn%par%ssa_solver_opt)

            ! Apply relaxation to keep things stable
            call relax_ssa(dyn%now%ux_b,dyn%now%uy_b,ux_b_prev,uy_b_prev,rel=dyn%par%ssa_iter_rel)

            ! Check for convergence
            is_converged = check_vel_convergence_l2rel(dyn%now%ux_b,dyn%now%uy_b,ux_b_prev,uy_b_prev, &
                                        dyn%now%ssa_mask_acx.gt.0.0_prec,dyn%now%ssa_mask_acy.gt.0.0_prec, &
                                        dyn%par%ssa_iter_conv,iter,dyn%par%ssa_iter_max,yelmo_log,use_L2_norm=.FALSE.)

            ! Calculate an L1 error metric over matrix for diagnostics
            call check_vel_convergence_l1rel_matrix(dyn%now%ssa_err_acx,dyn%now%ssa_err_acy,dyn%now%ux_b,dyn%now%uy_b, &
                                                                                            ux_b_prev,uy_b_prev)

            
            if (write_ssa_diagnostics) then  
                call write_step_2D_ssa(tpo,dyn,"yelmo_ssa.nc",ux_b_prev,uy_b_prev,time=real(iter,prec))    
            end if 

            ! Exit iterations if ssa solution has converged
            if (is_converged) exit 

        end do 
        ! == END iterations ==

!         if (write_ssa_diagnostics) then 
!             stop 
!         end if 

        return 

    end subroutine calc_ydyn_ssa 

!     subroutine calc_ydyn_beta(dyn,tpo,mat,bnd)
!         ! Update beta based on parameter choices

!         implicit none
        
!         type(ydyn_class),   intent(INOUT) :: dyn
!         type(ytopo_class),  intent(IN)    :: tpo 
!         type(ymat_class),   intent(IN)    :: mat 
!         type(ybound_class), intent(IN)    :: bnd   

!         ! Local variables 
!         integer :: i, j, nx, ny 
!         real(prec), allocatable :: logbeta(:,:) 
        
!         nx = size(dyn%now%beta,1)
!         ny = size(dyn%now%beta,2)
!         allocate(logbeta(nx,ny))
!         logbeta = 0.0_prec 

!         ! 0. Calculate C_bed [Pa]
!         dyn%now%c_bed = dyn%now%cf_ref * dyn%now%N_eff 

!         ! 1. Apply beta method of choice 
!         select case(dyn%par%beta_method)

!             case(-1)
!                 ! beta (aa-nodes) has been defined externally - do nothing

!             case(0)
!                 ! Constant beta everywhere

!                 dyn%now%beta = dyn%par%beta_const 

!             case(1)
!                 ! Calculate beta from a linear law (simply set beta=c_bed/u0)

!                 dyn%now%beta = dyn%now%c_bed * (1.0_prec / dyn%par%beta_u0)

!             case(2)
!                 ! Calculate beta from the quasi-plastic power-law as defined by Bueler and van Pelt (2015)

!                 call calc_beta_aa_power_plastic(dyn%now%beta,dyn%now%ux_b,dyn%now%uy_b,dyn%now%c_bed,dyn%par%beta_q,dyn%par%beta_u0)
                
!             case(3)
!                 ! Calculate beta from regularized Coulomb law (Joughin et al., GRL, 2019)

!                 call calc_beta_aa_reg_coulomb(dyn%now%beta,dyn%now%ux_b,dyn%now%uy_b,dyn%now%c_bed,dyn%par%beta_q,dyn%par%beta_u0)
                
!             case DEFAULT 
!                 ! Not recognized 

!                 write(*,*) "calc_ydyn_beta:: Error: beta_method not recognized."
!                 write(*,*) "beta_method = ", dyn%par%beta_method
!                 stop 

!         end select 

!         ! 1a. Ensure beta is relatively smooth 
! !         call regularize2D(dyn%now%beta,tpo%now%H_ice,tpo%par%dx)
! !         call limit_gradient(dyn%now%beta,tpo%now%H_ice,tpo%par%dx,log=.TRUE.)

!         ! 2. Scale beta as it approaches grounding line 
!         select case(dyn%par%beta_gl_scale) 

!             case(0) 
!                 ! Apply fractional parameter at grounding line, no scaling when beta_gl_f=1.0

!                 call scale_beta_gl_fraction(dyn%now%beta,tpo%now%f_grnd,dyn%par%beta_gl_f)

!             case(1) 
!                 ! Apply H_grnd scaling, reducing beta linearly towards zero at the grounding line 

!                 call scale_beta_gl_Hgrnd(dyn%now%beta,tpo%now%H_grnd,dyn%par%H_grnd_lim)

!             case(2) 
!                 ! Apply scaling according to thickness above flotation (Zstar approach of Gladstone et al., 2017)
!                 ! norm==.TRUE., so that zstar-scaling is bounded between 0 and 1, and thus won't affect 
!                 ! choice of c_bed value that is independent of this scaling. 
                
!                 call scale_beta_gl_zstar(dyn%now%beta,tpo%now%H_ice,bnd%z_bed,bnd%z_sl,norm=.TRUE.)

!             case DEFAULT 
!                 ! No scaling

!                 write(*,*) "calc_ydyn_beta:: Error: beta_gl_scale not recognized."
!                 write(*,*) "beta_gl_scale = ", dyn%par%beta_gl_scale
!                 stop 

!         end select 

!         ! 3. Apply grounding-line sub-element parameterization (sep, ie, subgrid method)
!         select case(dyn%par%beta_gl_sep)

!             case(-1) 
!                 ! Apply ac-node subgrid treatment, therefore do nothing here 
                
!                 !  Simply set beta to zero where purely floating
!                 where (tpo%now%f_grnd .eq. 0.0) dyn%now%beta = 0.0 
                
!                 if (dyn%par%beta_gl_stag .ne. 3) then 
!                     write(*,*) "calc_ydyn_beta:: Error: beta_gl_stag must equal 3 for beta_gl_sep=-1."
!                     write(*,*) "beta_gl_sep  = ", dyn%par%beta_gl_sep
!                     write(*,*) "beta_gl_stag = ", dyn%par%beta_gl_stag
!                     stop 
!                 end if 

!                 if (tpo%par%gl_sep .ne. 1) then 
!                     write(*,*) "calc_ydyn_beta:: Error: gl_sep must equal 1 for beta_gl_sep=-1."
!                     write(*,*) "beta_gl_sep  = ", dyn%par%beta_gl_sep
!                     write(*,*) "gl_sep       = ", tpo%par%gl_sep
!                     stop 
!                 end if 
                
!             case(0)
!                 ! No subgrid weighting at the grounding line,
!                 ! simply set beta to zero where purely floating  

!                 where (tpo%now%f_grnd .eq. 0.0) dyn%now%beta = 0.0 

!             case(1)
!                 ! Apply aa-node subgrid grounded fraction 
!                 dyn%now%beta = dyn%now%beta * tpo%now%f_grnd 

!             case DEFAULT 

!                 write(*,*) "calc_ydyn_beta:: Error: beta_gl_sep not recognized."
!                 write(*,*) "beta_gl_sep = ", dyn%par%beta_gl_sep
!                 stop 

!         end select 

!         ! 4. Apply smoothing if desired (only for points with beta > 0)
!         if (dyn%par%n_sm_beta .gt. 0) then
!             logbeta = 0.0_prec  
!             where(dyn%now%beta.gt.0.0_prec) logbeta = log10(dyn%now%beta)
!             call smooth_gauss_2D(logbeta,logbeta.gt.0.0,dyn%par%dx,dyn%par%n_sm_beta,logbeta.gt.0.0)
!             where(dyn%now%beta.gt.0.0_prec) dyn%now%beta = 10.0_prec**logbeta
!         end if 

!         ! Apply additional condition for particular experiments
!         if (trim(dyn%par%boundaries) .eq. "EISMINT") then 
!             ! Redefine beta at the summit to reduce singularity
!             ! in symmetric EISMINT experiments with sliding active
!             i = (dyn%par%nx-1)/2 
!             j = (dyn%par%ny-1)/2
!             dyn%now%beta(i,j) = (dyn%now%beta(i-1,j)+dyn%now%beta(i+1,j) &
!                                     +dyn%now%beta(i,j-1)+dyn%now%beta(i,j+1)) / 4.0 
!         else if (trim(dyn%par%boundaries) .eq. "MISMIP3D") then 
!             ! Redefine beta at the summit to reduce singularity
!             ! in MISMIP symmetric experiments
!             dyn%now%beta(1,:) = dyn%now%beta(2,:) 
!         end if 

!         ! Finally ensure that beta is higher than the lower allowed limit
!         where(dyn%now%beta .gt. 0.0 .and. dyn%now%beta .lt. dyn%par%beta_min) dyn%now%beta = dyn%par%beta_min 

!         ! ================================================================
!         ! Note: At this point the beta_aa field is available with beta=0 
!         ! for floating points and beta > 0 for non-floating points
!         ! ================================================================
        
!         ! 5. Apply staggering method with particular care for the grounding line 
!         select case(dyn%par%beta_gl_stag) 

!             case(0) 
!                 ! Apply pure staggering everywhere (ac(i) = 0.5*(aa(i)+aa(i+1))
                
!                 call stagger_beta_aa_mean(dyn%now%beta_acx,dyn%now%beta_acy,dyn%now%beta)

!             case(1) 
!                 ! Apply upstream beta_aa value at ac-node with at least one neighbor H_grnd_aa > 0

!                 call stagger_beta_aa_upstream(dyn%now%beta_acx,dyn%now%beta_acy,dyn%now%beta,tpo%now%f_grnd)

!             case(2) 
!                 ! Apply downstream beta_aa value (==0.0) at ac-node with at least one neighbor H_grnd_aa > 0

!                 call stagger_beta_aa_downstream(dyn%now%beta_acx,dyn%now%beta_acy,dyn%now%beta,tpo%now%f_grnd)

!             case(3)
!                 ! Apply subgrid scaling fraction at the grounding line when staggering 

!                 ! Note: now subgrid treatment is handled on aa-nodes above (using beta_gl_sep)

!                 call stagger_beta_aa_subgrid(dyn%now%beta_acx,dyn%now%beta_acy,dyn%now%beta,tpo%now%f_grnd, &
!                                                 tpo%now%f_grnd_acx,tpo%now%f_grnd_acy)

! !                 call stagger_beta_aa_subgrid_1(dyn%now%beta_acx,dyn%now%beta_acy,dyn%now%beta,tpo%now%H_grnd, &
! !                                             tpo%now%f_grnd,tpo%now%f_grnd_acx,tpo%now%f_grnd_acy)
                
!             case DEFAULT 

!                 write(*,*) "calc_ydyn_beta:: Error: beta_gl_stag not recognized."
!                 write(*,*) "beta_gl_stag = ", dyn%par%beta_gl_stag
!                 stop 

!         end select 

!         return 

!     end subroutine calc_ydyn_beta 

    subroutine calc_ydyn_cfref(dyn,tpo,thrm,bnd)
        ! Update cf_ref [--] based on parameter choices

        implicit none
        
        type(ydyn_class),   intent(INOUT) :: dyn
        type(ytopo_class),  intent(IN)    :: tpo 
        type(ytherm_class), intent(IN)    :: thrm
        type(ybound_class), intent(IN)    :: bnd  

        integer :: i, j, nx, ny 
        integer :: i1, i2, j1, j2 
        real(prec) :: f_scale 
        real(prec), allocatable :: cf_ref(:,:) 
        real(prec), allocatable :: lambda_bed(:,:)  

        nx = size(dyn%now%cf_ref,1)
        ny = size(dyn%now%cf_ref,2)
        
        allocate(cf_ref(nx,ny))
        allocate(lambda_bed(nx,ny))
        
        if (dyn%par%cb_method .eq. -1) then 
            ! Do nothing - cf_ref defined externally

        else 
            ! Calculate cf_ref following parameter choices 

            ! =============================================================================
            ! Step 1: calculate the cf_ref field only determined by 
            ! cf_frozen, cf_stream and temperate character of the bed 

            if (dyn%par%cb_with_pmp) then 
                ! Smooth transition between temperate and frozen cf_ref

                cf_ref = (thrm%now%f_pmp)*dyn%par%cf_stream &
                           + (1.0_prec - thrm%now%f_pmp)*dyn%par%cf_frozen 

            else 
                ! Only use cf_stream everywhere

                cf_ref = dyn%par%cf_stream

            end if 

            if (dyn%par%cb_with_pmp .and. dyn%par%cb_margin_pmp) then 
                ! Ensure that both the margin points and the grounding line
                ! are always considered streaming, independent of their
                ! thermodynamic character (as sometimes these can incorrectly become frozen)

            
                ! Ensure any marginal point is also treated as streaming 
                do j = 1, ny 
                do i = 1, nx 

                    i1 = max(i-1,1)
                    i2 = min(i+1,nx)
                    j1 = max(j-1,1)
                    j2 = min(j+1,ny)

                    if (tpo%now%H_ice(i,j) .gt. 0.0 .and. &
                        (tpo%now%H_ice(i1,j) .le. 0.0 .or. &
                         tpo%now%H_ice(i2,j) .le. 0.0 .or. &
                         tpo%now%H_ice(i,j1) .le. 0.0 .or. &
                         tpo%now%H_ice(i,j2) .le. 0.0)) then 
                        
                        cf_ref(i,j) = dyn%par%cf_stream

                    end if 

                end do 
                end do 

                ! Also ensure that grounding line is also considered streaming
                ! Note: this was related to cold ocean temps at floating-grounded interface,
                ! which is likely solved. Left here for safety. ajr, 2019-07-24
                where(tpo%now%is_grline) cf_ref = dyn%par%cf_stream

            end if

            ! =============================================================================
            ! Step 2: calculate lambda functions to scale cf_ref from default value 
            
            !------------------------------------------------------------------------------
            ! lambda_bed: scaling as a function of bedrock elevation

            select case(trim(dyn%par%cb_scale))

                case("lin_zb")
                    ! Linear scaling function with bedrock elevation
                    
                    lambda_bed = calc_lambda_bed_lin(bnd%z_bed,dyn%par%cb_z0,dyn%par%cb_z1)

                case("exp_zb")
                    ! Exponential scaling function with bedrock elevation
                    
                    lambda_bed = calc_lambda_bed_exp(bnd%z_bed,dyn%par%cb_z0,dyn%par%cb_z1)

                case("till_const")
                    ! Constant till friction angle

                    lambda_bed = calc_lambda_till_const(dyn%par%till_phi_const)

                case("till_zb")
                    ! Linear till friction angle versus elevation

                    lambda_bed = calc_lambda_till_linear(bnd%z_bed,bnd%z_sl,dyn%par%till_phi_min,dyn%par%till_phi_max, &
                                                            dyn%par%till_phi_zmin,dyn%par%till_phi_zmax)

                case DEFAULT
                    ! No scaling

                    lambda_bed = 1.0_prec

            end select 
            
            ! Ensure lambda_bed is not below lower limit [default range 0:1] 
            where (lambda_bed .lt. dyn%par%cb_min) lambda_bed = dyn%par%cb_min


            ! =============================================================================
            ! Step 3: calculate cf_ref [--]
            
            dyn%now%cf_ref = (cf_ref*lambda_bed)
            
        end if 

        return 

    end subroutine calc_ydyn_cfref

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
            H_w = dyn%par%neff_w_max * thrm%now%f_pmp  
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

                dyn%now%N_eff = calc_effective_pressure_overburden(tpo%now%H_ice,tpo%now%f_grnd)

            case(2) 
                ! Effective pressure diminishes with marine character
                ! following Leguy et al. (2014) 

                dyn%now%N_eff = calc_effective_pressure_marine(tpo%now%H_ice,bnd%z_bed,bnd%z_sl,H_w,p=dyn%par%neff_p)

            case(3)
                ! Effective pressure as basal till pressure
                ! following van Pelt and Bueler (2015)

                dyn%now%N_eff = calc_effective_pressure_till(H_w,tpo%now%H_ice,tpo%now%f_grnd,dyn%par%neff_w_max, &
                                            dyn%par%neff_N0,dyn%par%neff_delta,dyn%par%neff_e0,dyn%par%neff_Cc) 

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
        call nml_read(filename,"ydyn","sia_solver",         par%sia_solver,         init=init_pars)
        call nml_read(filename,"ydyn","ssa_solver_opt",     par%ssa_solver_opt,     init=init_pars)
        call nml_read(filename,"ydyn","mix_method",         par%mix_method,         init=init_pars)
        call nml_read(filename,"ydyn","calc_diffusivity",   par%calc_diffusivity,   init=init_pars)
        call nml_read(filename,"ydyn","diva_no_slip",       par%diva_no_slip,       init=init_pars)
        call nml_read(filename,"ydyn","beta_method",        par%beta_method,        init=init_pars)
        call nml_read(filename,"ydyn","beta_const",         par%beta_const,         init=init_pars)
        call nml_read(filename,"ydyn","beta_q",             par%beta_q,             init=init_pars)
        call nml_read(filename,"ydyn","beta_u0",            par%beta_u0,            init=init_pars)
        call nml_read(filename,"ydyn","beta_gl_scale",      par%beta_gl_scale,      init=init_pars)
        call nml_read(filename,"ydyn","beta_gl_sep",        par%beta_gl_sep,        init=init_pars)
        call nml_read(filename,"ydyn","beta_gl_stag",       par%beta_gl_stag,       init=init_pars)
        call nml_read(filename,"ydyn","beta_gl_f",          par%beta_gl_f,          init=init_pars)
        call nml_read(filename,"ydyn","taud_gl_method",     par%taud_gl_method,     init=init_pars)
        call nml_read(filename,"ydyn","H_grnd_lim",         par%H_grnd_lim,         init=init_pars)
        call nml_read(filename,"ydyn","H_sed_sat",          par%H_sed_sat,          init=init_pars)
        call nml_read(filename,"ydyn","cb_method",          par%cb_method,          init=init_pars)
        call nml_read(filename,"ydyn","cb_with_pmp",        par%cb_with_pmp,        init=init_pars)
        call nml_read(filename,"ydyn","cb_margin_pmp",      par%cb_margin_pmp,      init=init_pars)
        call nml_read(filename,"ydyn","cb_scale",           par%cb_scale,           init=init_pars)
        call nml_read(filename,"ydyn","cb_z0",              par%cb_z0,              init=init_pars)
        call nml_read(filename,"ydyn","cb_z1",              par%cb_z1,              init=init_pars)
        call nml_read(filename,"ydyn","cb_min",             par%cb_min,             init=init_pars)
        call nml_read(filename,"ydyn","cf_frozen",          par%cf_frozen,          init=init_pars)
        call nml_read(filename,"ydyn","cf_stream",          par%cf_stream,          init=init_pars)
        call nml_read(filename,"ydyn","n_sm_beta",          par%n_sm_beta,          init=init_pars)
        call nml_read(filename,"ydyn","beta_min",           par%beta_min,           init=init_pars)
        call nml_read(filename,"ydyn","ssa_beta_max",       par%ssa_beta_max,       init=init_pars)
        call nml_read(filename,"ydyn","ssa_vel_max",        par%ssa_vel_max,        init=init_pars)
        call nml_read(filename,"ydyn","ssa_iter_max",       par%ssa_iter_max,       init=init_pars)
        call nml_read(filename,"ydyn","ssa_iter_rel",       par%ssa_iter_rel,       init=init_pars)
        call nml_read(filename,"ydyn","ssa_iter_conv",      par%ssa_iter_conv,      init=init_pars)
        call nml_read(filename,"ydyn","taud_lim",           par%taud_lim,           init=init_pars)
        call nml_read(filename,"ydyn","cb_sia",             par%cb_sia,             init=init_pars)
        
        call nml_read(filename,"ydyn_till","till_phi_const",par%till_phi_const,     init=init_pars)
        call nml_read(filename,"ydyn_till","till_phi_min",  par%till_phi_min,       init=init_pars)
        call nml_read(filename,"ydyn_till","till_phi_max",  par%till_phi_max,       init=init_pars)
        call nml_read(filename,"ydyn_till","till_phi_zmin", par%till_phi_zmin,      init=init_pars)
        call nml_read(filename,"ydyn_till","till_phi_zmax", par%till_phi_zmax,      init=init_pars)
        
        call nml_read(filename,"ydyn_neff","neff_method",   par%neff_method,        init=init_pars)
        call nml_read(filename,"ydyn_neff","neff_const",    par%neff_const,         init=init_pars)
        call nml_read(filename,"ydyn_neff","neff_p",        par%neff_p,             init=init_pars)
        call nml_read(filename,"ydyn_neff","neff_set_water",par%neff_set_water,     init=init_pars)
        call nml_read(filename,"ydyn_neff","neff_w_max",    par%neff_w_max,         init=init_pars)
        call nml_read(filename,"ydyn_neff","neff_N0",       par%neff_N0,            init=init_pars)
        call nml_read(filename,"ydyn_neff","neff_delta",    par%neff_delta,         init=init_pars)
        call nml_read(filename,"ydyn_neff","neff_e0",       par%neff_e0,            init=init_pars)
        call nml_read(filename,"ydyn_neff","neff_Cc",       par%neff_Cc,            init=init_pars)

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

        ! Specifically deactivate ssa for mix_method=-2 
        if (par%mix_method .eq. -2) par%use_ssa = .FALSE. 

        ! Define current time as unrealistic value
        par%time = 1000000000   ! [a] 1 billion years in the future
        
        return

    end subroutine ydyn_par_load

    subroutine ydyn_alloc(now,nx,ny,nz_aa,nz_ac)

        implicit none 

        type(ydyn_state_class), intent(INOUT) :: now 
        integer, intent(IN) :: nx, ny, nz_aa, nz_ac   

        call ydyn_dealloc(now)

        allocate(now%ux(nx,ny,nz_aa)) 
        allocate(now%uy(nx,ny,nz_aa)) 
        allocate(now%uxy(nx,ny,nz_aa)) 
        allocate(now%uz(nx,ny,nz_ac)) 

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

        allocate(now%dd_ab(nx,ny,nz_aa))
        allocate(now%dd_ab_bar(nx,ny))

        allocate(now%sigma_horiz_sq(nx,ny))
        allocate(now%lhs_x(nx,ny)) 
        allocate(now%lhs_y(nx,ny)) 
        allocate(now%lhs_xy(nx,ny)) 
        
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

        allocate(now%qq_gl_acx(nx,ny)) 
        allocate(now%qq_gl_acy(nx,ny)) 

        allocate(now%qq_acx(nx,ny)) 
        allocate(now%qq_acy(nx,ny)) 
        allocate(now%qq(nx,ny)) 

        allocate(now%visc_eff(nx,ny,nz_aa))  
        allocate(now%visc_eff_int(nx,ny))

        allocate(now%cf_ref(nx,ny))
        allocate(now%c_bed(nx,ny)) 
        
        allocate(now%N_eff(nx,ny))

        allocate(now%beta_acx(nx,ny))
        allocate(now%beta_acy(nx,ny))
        allocate(now%beta(nx,ny))
        
        allocate(now%beta_eff(nx,ny))
        allocate(now%beta_diva(nx,ny))

        allocate(now%f_vbvs(nx,ny)) 
        
        allocate(now%ssa_mask_acx(nx,ny)) 
        allocate(now%ssa_mask_acy(nx,ny)) 
        allocate(now%ssa_err_acx(nx,ny)) 
        allocate(now%ssa_err_acy(nx,ny)) 
        
        ! Set all variables to zero intially
        now%ux                = 0.0 
        now%uy                = 0.0 
        now%uxy               = 0.0 
        now%uz                = 0.0 

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

        now%dd_ab             = 0.0
        now%dd_ab_bar         = 0.0
        
        now%sigma_horiz_sq    = 0.0
        now%lhs_x             = 0.0 
        now%lhs_y             = 0.0 
        now%lhs_xy            = 0.0 
        
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
        
        now%qq_gl_acx         = 0.0 
        now%qq_gl_acy         = 0.0 
        
        now%qq_acx            = 0.0 
        now%qq_acy            = 0.0 
        now%qq                = 0.0 
        
        now%visc_eff          = 1e3  
        now%visc_eff_int      = 1e3  
        
        now%cf_ref            = 0.0
        now%c_bed             = 0.0 
        
        now%N_eff             = 0.0 

        now%beta_acx          = 0.0 
        now%beta_acy          = 0.0 
        now%beta              = 0.0 
        
        now%beta_eff          = 0.0 
        now%beta_diva         = 0.0 

        now%f_vbvs            = 0.0 

        now%ssa_mask_acx      = 0.0 
        now%ssa_mask_acy      = 0.0 
        now%ssa_err_acx       = 0.0 
        now%ssa_err_acy       = 0.0 

        return 

    end subroutine ydyn_alloc 

    subroutine ydyn_dealloc(now)

        implicit none 

        type(ydyn_state_class), intent(INOUT) :: now

        if (allocated(now%ux))              deallocate(now%ux) 
        if (allocated(now%uy))              deallocate(now%uy) 
        if (allocated(now%uxy))             deallocate(now%uxy) 
        if (allocated(now%uz))              deallocate(now%uz) 

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
        
        if (allocated(now%dd_ab))           deallocate(now%dd_ab)
        if (allocated(now%dd_ab_bar))       deallocate(now%dd_ab_bar)
        
        if (allocated(now%sigma_horiz_sq))  deallocate(now%sigma_horiz_sq)
        if (allocated(now%lhs_x))           deallocate(now%lhs_x) 
        if (allocated(now%lhs_y))           deallocate(now%lhs_y) 
        if (allocated(now%lhs_xy))          deallocate(now%lhs_xy) 
        
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
        
        if (allocated(now%qq_gl_acx))       deallocate(now%qq_gl_acx) 
        if (allocated(now%qq_gl_acy))       deallocate(now%qq_gl_acy) 
        
        if (allocated(now%qq_acx))          deallocate(now%qq_acx) 
        if (allocated(now%qq_acy))          deallocate(now%qq_acy) 
        if (allocated(now%qq))              deallocate(now%qq) 
        
        if (allocated(now%visc_eff))        deallocate(now%visc_eff) 
        if (allocated(now%visc_eff_int))    deallocate(now%visc_eff_int) 
        
        if (allocated(now%cf_ref))          deallocate(now%cf_ref) 
        if (allocated(now%c_bed))           deallocate(now%c_bed) 
        
        if (allocated(now%N_eff))           deallocate(now%N_eff)
        
        if (allocated(now%beta_acx))        deallocate(now%beta_acx) 
        if (allocated(now%beta_acy))        deallocate(now%beta_acy) 
        if (allocated(now%beta))            deallocate(now%beta) 
        
        if (allocated(now%beta_eff))        deallocate(now%beta_eff) 
        if (allocated(now%beta_diva))       deallocate(now%beta_diva) 

        if (allocated(now%f_vbvs))          deallocate(now%f_vbvs) 

        if (allocated(now%ssa_mask_acx))    deallocate(now%ssa_mask_acx) 
        if (allocated(now%ssa_mask_acy))    deallocate(now%ssa_mask_acy) 
        if (allocated(now%ssa_err_acx))     deallocate(now%ssa_err_acx) 
        if (allocated(now%ssa_err_acy))     deallocate(now%ssa_err_acy) 

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
        call nc_write_dim(filename,"xc",     x=0.0_prec,dx=1.0_prec,nx=nx,units="gridpoints")
        call nc_write_dim(filename,"yc",     x=0.0_prec,dx=1.0_prec,nx=ny,units="gridpoints")
        call nc_write_dim(filename,"time",   x=time_init,dx=1.0_prec,nx=1,units="iter",unlimited=.TRUE.)

        return

    end subroutine yelmo_write_init_ssa 

    subroutine write_step_2D_ssa(tpo,dyn,filename,ux_b_prev,uy_b_prev,time)

        implicit none 
        
        type(ytopo_class), intent(IN) :: tpo 
        type(ydyn_class),  intent(IN) :: dyn 
        character(len=*),  intent(IN) :: filename
        real(prec), intent(IN) :: ux_b_prev(:,:) 
        real(prec), intent(IN) :: uy_b_prev(:,:) 
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
        
        call nc_write(filename,"dzsrfdt",tpo%now%dzsrfdt,units="m/a",long_name="Surface elevation change", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dHicedt",tpo%now%dHicedt,units="m/a",long_name="Ice thickness change", &
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
        
        call nc_write(filename,"sigma_horiz_sq",dyn%now%sigma_horiz_sq,units="1",long_name="Horizontal stress components squared", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"lhs_x",dyn%now%lhs_x,units="Pa",long_name="Shear reduction (x)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"lhs_y",dyn%now%lhs_y,units="Pa",long_name="Shear reduction (y)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"lhs_xy",dyn%now%lhs_xy,units="Pa",long_name="Shear reduction magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"ux_i_bar",dyn%now%ux_i_bar,units="m/a",long_name="Internal shear velocity (x)", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"uy_i_bar",dyn%now%uy_i_bar,units="m/a",long_name="Internal shear velocity (y)", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"uxy_i_bar",dyn%now%uxy_i_bar,units="m/a",long_name="Internal shear velocity magnitude", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"ux_b",dyn%now%ux_b,units="m/a",long_name="Basal sliding velocity (x)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy_b",dyn%now%uy_b,units="m/a",long_name="Basal sliding velocity (y)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_b",dyn%now%uxy_b,units="m/a",long_name="Basal sliding velocity magnitude", &
                     dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"ux_b_diff",dyn%now%ux_b-ux_b_prev,units="m/a",long_name="Basal sliding velocity difference (x)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy_b_diff",dyn%now%uy_b-uy_b_prev,units="m/a",long_name="Basal sliding velocity difference (y)", &
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

