module velocity_diva_ab

    use yelmo_defs, only  : sp, dp, prec, wp, pi, rho_ice, rho_sw, rho_w, g, TOL_UNDERFLOW
    use yelmo_tools, only : stagger_aa_ab, stagger_aa_ab_ice, stagger_ab_aa_ice, & 
                    stagger_nodes_aa_ab_ice, stagger_nodes_acx_ab_ice, stagger_nodes_acy_ab_ice, &
                    staggerdiff_nodes_acx_ab_ice, staggerdiff_nodes_acy_ab_ice, &
                    staggerdiffcross_nodes_acx_ab_ice, staggerdiffcross_nodes_acy_ab_ice, &
                    integrate_trapezoid1D_1D, integrate_trapezoid1D_pt, minmax, &
                    set_boundaries_2D_aa, set_boundaries_3D_aa, set_boundaries_3D_acx, &
                    set_boundaries_3D_acy

    use basal_dragging 
    use solver_ssa_ab
    use velocity_general, only : set_inactive_margins, &
                        picard_calc_error, picard_calc_error_angle, picard_relax, &
                        picard_calc_convergence_l1rel_matrix, picard_calc_convergence_l2 

    use grid_calcs 

    implicit none 

    type diva_param_class

        character(len=256) :: ssa_lis_opt 
        character(len=256) :: ssa_lateral_bc
        character(len=256) :: boundaries 
        character(len=256) :: glf_method 
        logical  :: no_slip 
        integer  :: visc_method
        real(wp) :: visc_const
        integer  :: beta_method
        real(wp) :: beta_const
        real(wp) :: beta_q                  ! Friction law exponent
        real(wp) :: beta_u0                 ! [m/a] Friction law velocity threshold 
        integer  :: beta_gl_scale           ! Beta grounding-line scaling method (beta => 0 at gl?)
        integer  :: beta_gl_stag            ! Beta grounding-line staggering method 
        real(wp) :: beta_gl_f               ! Fraction of beta at gl 
        real(wp) :: H_grnd_lim 
        real(wp) :: beta_min                ! Minimum allowed value of beta
        real(wp) :: eps_0 
        real(wp) :: ssa_vel_max
        integer  :: ssa_iter_max 
        real(wp) :: ssa_iter_rel 
        real(wp) :: ssa_iter_conv 
        logical  :: ssa_write_log 

        real(wp) :: glf_Q0                  ! Q0=0.61_wp
        real(wp) :: glf_f_drag              ! f_drag=0.6_wp
    end type

    private
    public :: diva_param_class 
    public :: calc_velocity_diva_ab

contains
    
    ! Variables/params needed for gl-flux methods
    ! qq_gl_acx 
    ! qq_gl_acy 
    ! ATT_bar 
    ! 
    
    subroutine calc_velocity_diva_ab(ux,uy,ux_bar,uy_bar,ux_b,uy_b,ux_i,uy_i,taub_acx,taub_acy, &
                                  beta,beta_acx,beta_acy, &
                                  ux_bar_ab,uy_bar_ab, &
                                  beta_eff,de_eff,visc_eff,visc_eff_int,duxdz,duydz, &
                                  ssa_mask_acx,ssa_mask_acy,ssa_err_acx,ssa_err_acy,ssa_iter_now, &
                                  c_bed,taud_acx,taud_acy,H_ice,f_ice,H_grnd,f_grnd, &
                                  f_grnd_acx,f_grnd_acy,ATT,zeta_aa,z_sl,z_bed,z_srf,dx,dy,n_glen,par)
        ! This subroutine is used to solve the horizontal velocity system (ux,uy)
        ! following the Depth-Integrated Viscosity Approximation (DIVA),
        ! as outlined by Lipscomb et al. (2019). Method originally 
        ! proposed by Goldberg (2011), algorithm by Arthern et al (2015), 
        ! updated by Lipscomb et al. (2019).

        implicit none 

        real(wp), intent(INOUT) :: ux(:,:,:)          ! [m/a]
        real(wp), intent(INOUT) :: uy(:,:,:)          ! [m/a]
        real(wp), intent(INOUT) :: ux_bar(:,:)        ! [m/a]
        real(wp), intent(INOUT) :: uy_bar(:,:)        ! [m/a]
        real(wp), intent(INOUT) :: ux_b(:,:)          ! [m/a]
        real(wp), intent(INOUT) :: uy_b(:,:)          ! [m/a]
        real(wp), intent(INOUT) :: ux_i(:,:,:)        ! [m/a]
        real(wp), intent(INOUT) :: uy_i(:,:,:)        ! [m/a]
        real(wp), intent(INOUT) :: taub_acx(:,:)      ! [Pa]
        real(wp), intent(INOUT) :: taub_acy(:,:)      ! [Pa]
        real(wp), intent(INOUT) :: beta(:,:)          ! [Pa a/m]
        real(wp), intent(INOUT) :: beta_acx(:,:)      ! [Pa a/m]
        real(wp), intent(INOUT) :: beta_acy(:,:)      ! [Pa a/m]
        real(wp), intent(INOUT) :: ux_bar_ab(:,:)     ! [m/a]
        real(wp), intent(INOUT) :: uy_bar_ab(:,:)     ! [m/a]
        real(wp), intent(INOUT) :: beta_eff(:,:)      ! [Pa a/m]
        real(wp), intent(INOUT) :: de_eff(:,:,:)      ! [1/a]
        real(wp), intent(INOUT) :: visc_eff(:,:,:)    ! [Pa a]
        real(wp), intent(INOUT) :: visc_eff_int(:,:)  ! [Pa a m]
        real(wp), intent(INOUT) :: duxdz(:,:,:)       ! [1/a]
        real(wp), intent(INOUT) :: duydz(:,:,:)       ! [1/a]
        integer,  intent(INOUT) :: ssa_mask_acx(:,:)  ! [-]
        integer,  intent(INOUT) :: ssa_mask_acy(:,:)  ! [-]
        real(wp), intent(OUT)   :: ssa_err_acx(:,:)
        real(wp), intent(OUT)   :: ssa_err_acy(:,:)
        integer,  intent(OUT)   :: ssa_iter_now 
        real(wp), intent(IN)    :: c_bed(:,:)         ! [Pa]
        real(wp), intent(IN)    :: taud_acx(:,:)      ! [Pa]
        real(wp), intent(IN)    :: taud_acy(:,:)      ! [Pa]
        real(wp), intent(IN)    :: H_ice(:,:)         ! [m]
        real(wp), intent(IN)    :: f_ice(:,:)         ! [--]
        real(wp), intent(IN)    :: H_grnd(:,:)        ! [m]
        real(wp), intent(IN)    :: f_grnd(:,:)        ! [-]
        real(wp), intent(IN)    :: f_grnd_acx(:,:)    ! [-]
        real(wp), intent(IN)    :: f_grnd_acy(:,:)    ! [-]
        real(wp), intent(IN)    :: ATT(:,:,:)         ! [a^-1 Pa^-n_glen]
        real(wp), intent(IN)    :: zeta_aa(:)         ! [-]
        real(wp), intent(IN)    :: z_sl(:,:)          ! [m]
        real(wp), intent(IN)    :: z_bed(:,:)         ! [m]
        real(wp), intent(IN)    :: z_srf(:,:)         ! [m]
        real(wp), intent(IN)    :: dx                 ! [m]
        real(wp), intent(IN)    :: dy                 ! [m]
        real(wp), intent(IN)    :: n_glen 
        type(diva_param_class), intent(IN) :: par       ! List of parameters that should be defined

        ! Local variables 
        integer :: i, j, k, nx, ny, nz_aa, nz_ac, iter 
        logical :: is_converged

        real(wp), allocatable :: ux_bar_nm1(:,:) 
        real(wp), allocatable :: uy_bar_nm1(:,:)  
        real(wp), allocatable :: F2(:,:)              ! [Pa^-1 a^-1 m == (Pa a/m)^-1]
        integer,  allocatable :: ssa_mask_acx_ref(:,:)
        integer,  allocatable :: ssa_mask_acy_ref(:,:)

        real(wp), allocatable :: corr_nm1(:) 
        real(wp), allocatable :: corr_nm2(:) 
        
        real(wp), allocatable :: H_ice_ab(:,:) 
        real(wp), allocatable :: taudx_ab(:,:) 
        real(wp), allocatable :: taudy_ab(:,:) 
        real(wp), allocatable :: taubx_ab(:,:) 
        real(wp), allocatable :: tauby_ab(:,:) 
        real(wp), allocatable :: ux_b_ab(:,:) 
        real(wp), allocatable :: uy_b_ab(:,:) 
        
        ! For glf methods 
        logical :: with_glf 
        real(wp), allocatable :: ATT_bar(:,:) 
        real(wp), allocatable :: qq_gl_acx(:,:) 
        real(wp), allocatable :: qq_gl_acy(:,:) 
        
        real(wp) :: corr_theta
        real(wp) :: corr_rel 
        real(wp) :: L2_norm 
        real(wp) :: ssa_resid 

        integer  :: ntot, ip1, jp1  

        logical, parameter :: write_ssa_diagnostics      = .FALSE. 
        logical, parameter :: write_ssa_diagnostics_stop = .FALSE.   ! Stop simulation after completing iterations?

        nx    = size(ux,1)
        ny    = size(ux,2)
        nz_aa = size(ux,3)

        ! Prepare local variables 
        allocate(ux_bar_nm1(nx,ny))
        allocate(uy_bar_nm1(nx,ny))
        allocate(F2(nx,ny))

        allocate(ssa_mask_acx_ref(nx,ny))
        allocate(ssa_mask_acy_ref(nx,ny))

        allocate(corr_nm1(2*nx*ny))
        allocate(corr_nm2(2*nx*ny))

        allocate(H_ice_ab(nx,ny))
        allocate(taudx_ab(nx,ny))
        allocate(taudy_ab(nx,ny))
        allocate(taubx_ab(nx,ny))
        allocate(tauby_ab(nx,ny))
        allocate(ux_b_ab(nx,ny))
        allocate(uy_b_ab(nx,ny))

        with_glf = .FALSE. 
        if (trim(par%glf_method) .eq. "power") with_glf = .TRUE.

        ! Store original ssa mask before iterations
        ssa_mask_acx_ref = ssa_mask_acx
        ssa_mask_acy_ref = ssa_mask_acy
            
        ! Initially set error very high 
        ssa_err_acx = 1.0_wp 
        ssa_err_acy = 1.0_wp 
        
        corr_nm1 = 0.0_wp 
        corr_nm2 = 0.0_wp 

        ! Get ice thickness on ab-nodes 
        call map_a_to_b_2D(H_ice_ab,H_ice) 

        ! Get driving stress on ab-nodes 
        call map_cx_to_b_2D(taudx_ab,taud_acx)
        call map_cy_to_b_2D(taudy_ab,taud_acy)
        
        if (write_ssa_diagnostics) then 
            call ssa_diagnostics_write_init("yelmo_ssa.nc",nx,ny,time_init=1.0_wp)
        end if 

        ! Calculate basal stress again initially since it is only defined locally
        call calc_basal_stress_ab(taubx_ab,tauby_ab,beta_eff,ux_bar_ab,uy_bar_ab)

        ! Calculate basal velocity initially too
        call calc_vel_basal_ab(ux_b_ab,uy_b_ab,ux_bar_ab,uy_bar_ab,F2,taubx_ab,tauby_ab,par%no_slip)


        do iter = 1, par%ssa_iter_max 

            ! Store solution from previous iteration (nm1 == n minus 1) 
            ux_bar_nm1 = ux_bar_ab 
            uy_bar_nm1 = uy_bar_ab 
            
            ! =========================================================================================
            ! Step 1: Calculate fields needed by ssa solver (visc_eff_int, beta_eff)

            ! Calculate the 3D vertical shear fields using viscosity estimated from the previous iteration 
            call calc_vertical_shear_3D_ab(duxdz,duydz,taubx_ab,tauby_ab,visc_eff,zeta_aa,par%boundaries)

            ! Calculate 3D effective viscosity
            select case(par%visc_method)

                case(0)
                    ! Impose constant viscosity value 

                    visc_eff = par%visc_const 

                case(1) 
                    ! Calculate effective viscosity, using velocity solution from previous iteration
                    ! Use default staggering stencil (to ab-node quadrature points, then average to aa-nodes)

                    call calc_visc_eff_3D_ab(visc_eff,ux_bar_ab,uy_bar_ab,duxdz,duydz,ATT,zeta_aa, &
                                                            dx,dy,n_glen,par%eps_0,par%boundaries)

                case DEFAULT 

                    write(*,*) "calc_velocity_diva_ab:: Error: visc_method not recognized."
                    write(*,*) "visc_method = ", par%visc_method 
                    stop 

            end select
            
            ! Calculate depth-integrated effective viscosity
            ! Note L19 uses eta_bar*H in the ssa equation. Yelmo uses eta_int=eta_bar*H directly.
            call calc_visc_eff_int(visc_eff_int,visc_eff,H_ice_ab,zeta_aa,par%boundaries)
            
            ! Calculate beta (at the ice base)
            ! Right now only with power law for ab-node velocities
            call calc_beta_aa_power_plastic_ab(beta,ux_b_ab,uy_b_ab,c_bed,par%beta_q,par%beta_u0)

            ! Calculate F-integeral (F2) on ab-nodes 
            call calc_F_integral(F2,visc_eff,H_ice_ab,zeta_aa,n=2.0_wp)
            
            ! Calculate effective beta 
            call calc_beta_eff(beta_eff,beta,F2,zeta_aa,no_slip=par%no_slip)
            
            ! Set to zero, not needed:
            beta_acx = 0.0_wp 
            beta_acy = 0.0_wp 

            ! =========================================================================================
            ! Step 2: Call the SSA solver to obtain new estimate of ux_bar/uy_bar
            
            ! -------------------------------------------------------
            ! amc: adaptive mesh coarsening 
            ! On high res. grid, determine which points are elegible 
            ! for amc. 

            ! If enough points exist...
            ! Interpolate quantities onto low resolution grid 
            ! (ux_bar,uy_bar,beta_eff_acx,beta_eff_acy,visc_eff_int,ssa_mask_acx
            !  ssa_mask_acy,H_ice,taud_acx,taud_acy,H_grnd,z_sl,z_bed)

            ! Call solver for low resolution solution 

            ! Copy low resolution solution to high-res amc-elegible points 
            ! Interpolate high-res amc-elegible points in between those actually solved 

            ! Mask high-res amc-points to prescribe solution there 

            ! Solve high-res as normal 
            ! ------------------------------------------------------- 
            
            ! Call ssa solver
            call calc_vxy_ssa_matrix(ux_bar_ab,uy_bar_ab,L2_norm,beta_eff,visc_eff_int,  &
                                ssa_mask_acx,H_ice,f_ice,taudx_ab,taudy_ab,H_grnd,z_sl,z_bed, &
                                z_srf,dx,dy,par%ssa_vel_max,par%boundaries,par%ssa_lateral_bc,par%ssa_lis_opt)

            ! Get ux_bar/uy_bar on ac-nodes 
            call map_b_to_cx_2D(ux_bar,ux_bar_ab)
            call map_b_to_cy_2D(uy_bar,uy_bar_ab)
            
! ajr: use adapative ssa relaxation or use parameter value
! For Antarctica, the adaptive method can give some strange
! convergence issues. It has been disabled for now (2022-02-09).
if (.FALSE.) then
            ! Calculate errors 
            corr_nm2 = corr_nm1 
            call picard_calc_error(corr_nm1,ux_bar,uy_bar,ux_bar_nm1,uy_bar_nm1)

            ! Calculate error angle 
            call picard_calc_error_angle(corr_theta,corr_nm1,corr_nm2) 

            if (corr_theta .le. pi/8.0_wp) then 
                corr_rel = 2.5_wp 
            else if (corr_theta .lt. 19.0_wp*pi/20.0_wp) then 
                corr_rel = 1.0_wp 
            else 
                corr_rel = 0.5_wp 
            end if 
else
            corr_rel = par%ssa_iter_rel
end if
            !write(*,*) "pic: ", iter, corr_theta, corr_rel

            ! Apply relaxation to keep things stable
            !call relax_ssa(ux_bar,uy_bar,ux_bar_nm1,uy_bar_nm1,rel=par%ssa_iter_rel)
            call picard_relax(ux_bar_ab,uy_bar_ab,ux_bar_nm1,uy_bar_nm1,rel=corr_rel)
            
            ! Check for convergence
            ! is_converged = check_vel_convergence_l2rel(ux_bar,uy_bar,ux_bar_nm1,uy_bar_nm1,ssa_mask_acx.gt.0,     &
            !                                            ssa_mask_acy.gt.0,par%ssa_iter_conv,iter,par%ssa_iter_max, &
            !                                            par%ssa_write_log,use_L2_norm=.FALSE.,L2_norm=L2_norm)
            call picard_calc_convergence_l2(is_converged,ssa_resid,ux_bar_ab,uy_bar_ab,ux_bar_nm1,uy_bar_nm1, &
                                                ssa_mask_acx.gt.0,ssa_mask_acy.gt.0,par%ssa_iter_conv,  &
                                                iter,par%ssa_iter_max,par%ssa_write_log)

            ! Calculate an L1 error metric over matrix for diagnostics
            call picard_calc_convergence_l1rel_matrix(ssa_err_acx,ssa_err_acy,ux_bar_ab,uy_bar_ab,ux_bar_nm1,uy_bar_nm1)

            ! Store current total iterations for output
            ssa_iter_now = iter 

            if (write_ssa_diagnostics) then  
                call ssa_diagnostics_write_step("yelmo_ssa.nc",ux_bar,uy_bar,L2_norm,beta_eff,visc_eff_int, &
                                        ssa_mask_acx,ssa_err_acx,ssa_err_acy,H_ice,f_ice,taudx_ab,taudy_ab, &
                                        H_grnd,z_sl,z_bed,z_srf,ux_bar_nm1,uy_bar_nm1,time=real(iter,wp))    
            end if 

            ! =========================================================================================
            ! Update additional fields based on output of solver
             
            ! Calculate basal stress 
            call calc_basal_stress_ab(taubx_ab,tauby_ab,beta_eff,ux_bar_ab,uy_bar_ab)

            ! Calculate basal velocity from depth-averaged solution and basal stress
            call calc_vel_basal_ab(ux_b_ab,uy_b_ab,ux_bar_ab,uy_bar_ab,F2,taubx_ab,tauby_ab,par%no_slip)

            ! Exit iterations if ssa solution has converged
            if (is_converged) exit 
            
        end do 

        ! Iterations are finished, finalize calculations of 3D velocity field 

        if (write_ssa_diagnostics .and. write_ssa_diagnostics_stop) then 
            stop 
        end if 

        ! Calculate velocity on ac-nodes 
        call map_b_to_cx_2D(ux_b,ux_b_ab)
        call map_b_to_cy_2D(uy_b,uy_b_ab)
        call map_b_to_cx_2D(ux_bar,ux_bar_ab)
        call map_b_to_cy_2D(uy_bar,uy_bar_ab)
        call map_b_to_cx_2D(taub_acx,taubx_ab)
        call map_b_to_cy_2D(taub_acy,tauby_ab)
        
        ! Calculate the 3D horizontal velocity field
        call calc_vel_horizontal_3D(ux,uy,ux_b,uy_b,taub_acx,taub_acy,visc_eff,H_ice,f_ice,zeta_aa,par%boundaries)

        ! Also calculate the shearing contribution
        do k = 1, nz_aa 
            ux_i(:,:,k) = ux(:,:,k) - ux_b 
            uy_i(:,:,k) = uy(:,:,k) - uy_b 
        end do
        

        if (par%visc_method .eq. 0) then 
            ! Diagnose viscosity for visc_method=0 (not used prognostically)
            call calc_visc_eff_3D_ab(visc_eff,ux_bar_ab,uy_bar_ab,duxdz,duydz,ATT,zeta_aa, &
                                                            dx,dy,n_glen,par%eps_0,par%boundaries)
        end if 

        return 

    end subroutine calc_velocity_diva_ab

    subroutine calc_vel_horizontal_3D(ux,uy,ux_b,uy_b,taub_acx,taub_acy,visc_eff,H_ice,f_ice,zeta_aa,boundaries)
        ! Caluculate the 3D horizontal velocity field (ux,uy)
        ! following L19, Eq. 29 

        implicit none 

        real(wp), intent(OUT) :: ux(:,:,:) 
        real(wp), intent(OUT) :: uy(:,:,:) 
        real(wp), intent(IN)  :: ux_b(:,:) 
        real(wp), intent(IN)  :: uy_b(:,:) 
        real(wp), intent(IN)  :: taub_acx(:,:) 
        real(wp), intent(IN)  :: taub_acy(:,:)
        real(wp), intent(IN)  :: visc_eff(:,:,:)       
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: zeta_aa(:) 
        character(len=*), intent(IN) :: boundaries 

        ! Local variables
        integer :: i, j, k, ip1, jp1, nx, ny, nz_aa  
        real(wp), allocatable :: F1(:,:,:) 
        real(wp) :: F1_ac

        nx    = size(ux,1)
        ny    = size(ux,2) 
        nz_aa = size(ux,3) 

        allocate(F1(nx,ny,nz_aa))
        !allocate(F1_ac(nz_aa))

        F1 = 0.0_wp 

        ! First calculate F1 array on aa-nodes 
        ! (performing integral before staggering seems to improve result slightly)
        ! Note: L19 define the F1 integral as purely going from the base to the surface,
        ! whereas here F1 is calculated from the base to each point in the vertical. So, 
        ! it is not technically "F1" as defined by L19, Eq. 30, except at the surface.
        do j = 1, ny 
        do i = 1, nx 
            if (f_ice(i,j) .eq. 1.0) then 
                F1(i,j,:) = integrate_trapezoid1D_1D((H_ice(i,j)/visc_eff(i,j,:))*(1.0-zeta_aa),zeta_aa)
            end if  
        end do
        end do  

        ! Next calculate 3D horizontal velocity components 
        do k = 1, nz_aa
        do j = 1, ny 
        do i = 1, nx 

            ip1 = min(i+1,nx)
            jp1 = min(j+1,ny) 

            ! === x direction ===============================================

            ! Stagger F1 column to ac-nodes 
            if (f_ice(i,j) .eq. 1.0 .and. f_ice(ip1,j) .lt. 1.0) then 
                F1_ac = F1(i,j,k) 
            else if (f_ice(i,j) .lt. 1.0 .and. f_ice(ip1,j) .eq. 1.0) then
                F1_ac = F1(ip1,j,k)
            else 
                F1_ac = 0.5_wp*(F1(i,j,k) + F1(ip1,j,k))
            end if 

            ! Calculate velocity column 
            ux(i,j,k) = ux_b(i,j) + taub_acx(i,j)*F1_ac 

            ! === y direction ===============================================

            ! Stagger F1 column to ac-nodes 
            if (f_ice(i,j) .eq. 1.0 .and. f_ice(i,jp1) .lt. 1.0) then 
                F1_ac = F1(i,j,k) 
            else if (f_ice(i,j) .lt. 1.0 .and. f_ice(i,jp1) .eq. 1.0) then
                F1_ac = F1(i,jp1,k)
            else 
                F1_ac = 0.5_wp*(F1(i,j,k) + F1(i,jp1,k))
            end if 

            ! Calculate velocity column
            uy(i,j,k) = uy_b(i,j) + taub_acy(i,j)*F1_ac  

        end do 
        end do  
        end do

        ! Apply boundary conditions as needed 
        call set_boundaries_3D_acx(ux,boundaries)
        call set_boundaries_3D_acy(uy,boundaries)
        
        return 

    end subroutine calc_vel_horizontal_3D

    subroutine calc_vertical_shear_3D_ab(duxdz,duydz,taubx,tauby,visc_eff,zeta_aa,boundaries)
        ! Calculate vertical shear terms (L19, Eq. 36)

        implicit none 

        real(wp), intent(OUT) :: duxdz(:,:,:)           ! [1/yr],   ab-nodes horizontal, aa-nodes vertical 
        real(wp), intent(OUT) :: duydz(:,:,:)           ! [1/yr],   ab-nodes horizontal, aa-nodes vertical 
        real(wp), intent(IN)  :: taubx(:,:)             ! [Pa],     ab-nodes
        real(wp), intent(IN)  :: tauby(:,:)             ! [Pa],     ab-nodes
        real(wp), intent(IN)  :: visc_eff(:,:,:)        ! [Pa a m], ab-nodes
        real(wp), intent(IN)  :: zeta_aa(:)           ! [-]
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer  :: i, j, k, nx, ny, nz_aa 
        integer  :: ip1, jp1 

        ! real(wp), parameter :: visc_min = 1e5_wp 

        nx    = size(visc_eff,1)
        ny    = size(visc_eff,2)
        nz_aa = size(zeta_aa,1) 
        
        !$omp parallel do
        do k = 1, nz_aa 
        do j = 1, ny
        do i = 1, nx 

            ! Get staggering indices limited to grid size
            ip1 = min(i+1,nx)
            jp1 = min(j+1,ny) 
            
            ! Calculate shear strain (ab-nodes)
            if (visc_eff(i,j,k) .ne. 0.0_wp) then 
                duxdz(i,j,k) = (taubx(i,j)/visc_eff(i,j,k)) * (1.0_wp-zeta_aa(k))
                duydz(i,j,k) = (tauby(i,j)/visc_eff(i,j,k)) * (1.0_wp-zeta_aa(k))
            else 
                duxdz(i,j,k) = 0.0_wp 
                duydz(i,j,k) = 0.0_wp
            end if 

        end do 
        end do 
        end do 
        !$omp end parallel do

        return 

    end subroutine calc_vertical_shear_3D_ab

    subroutine calc_visc_eff_3D_ab(visc_eff,ux,uy,duxdz,duydz,ATT,zeta_aa,dx,dy,n_glen,eps_0,boundaries)
        ! Calculate 3D effective viscosity following L19, Eq. 2
        ! Use of eps_0 ensures non-zero positive viscosity value everywhere 
        ! Note: viscosity is first calculated on ab-nodes, then 
        ! unstaggered back to aa-nodes. This ensures more stability for 
        ! visc_eff (less likely to blow up for low strain rates). 

        implicit none 
        
        real(wp), intent(OUT) :: visc_eff(:,:,:)        ! ab-nodes
        real(wp), intent(IN)  :: ux(:,:)                ! [m/yr] Vertically averaged horizontal velocity, x-component
        real(wp), intent(IN)  :: uy(:,:)                ! [m/yr] Vertically averaged horizontal velocity, y-component
        real(wp), intent(IN)  :: duxdz(:,:,:)           ! [1/yr] Vertical shearing, x-component
        real(wp), intent(IN)  :: duydz(:,:,:)           ! [1/yr] Vertical shearing, x-component
        real(wp), intent(IN)  :: ATT(:,:,:)             ! aa-nodes
        real(wp), intent(IN)  :: zeta_aa(:)             ! Vertical axis (sigma-coordinates from 0 to 1)
        real(wp), intent(IN)  :: dx
        real(wp), intent(IN)  :: dy
        real(wp), intent(IN)  :: n_glen   
        real(wp), intent(IN)  :: eps_0                  ! [1/yr] Regularization constant (minimum strain rate, ~1e-6)
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer  :: i, j, k
        integer  :: ip1, jp1, im1, jm1 
        integer  :: nx, ny, nz  
        real(wp) :: inv_2dx, inv_2dy 
        real(wp) :: p1, p2, eps_0_sq  
        real(wp) :: eps_sq                              ! [1/yr^2]
        real(wp) :: dudx_ab
        real(wp) :: dvdy_ab
        real(wp) :: dudy_ab
        real(wp) :: dvdx_ab 
        real(wp) :: duxdz_ab
        real(wp) :: duydz_ab
        real(wp) :: eps_sq_ab
        real(wp) :: ATT_ab

        real(wp), parameter :: visc_min = 1e5_wp        ! Just for safety 

        nx = size(visc_eff,1)
        ny = size(visc_eff,2)
        nz = size(visc_eff,3)
        
        ! Calculate scaling factors
        inv_2dx = 1.0_wp / (2.0_wp*dx) 
        inv_2dy = 1.0_wp / (2.0_wp*dy) 

        ! Calculate exponents 
        p1 = (1.0_wp - n_glen)/(2.0_wp*n_glen)
        p2 = -1.0_wp/n_glen

        ! Calculate squared minimum strain rate 
        eps_0_sq = eps_0*eps_0 

        ! Calculate visc_eff on quadrature points (ab-nodes),
        ! and then average to get grid-centered value (aa-node)
        do j = 1, ny 
        do i = 1, nx 

            im1 = max(i-1,1)    
            ip1 = min(i+1,nx) 
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny) 

            ! Get strain rate terms (ab-nodes)
            dudx_ab = (ux(ip1,j)-ux(im1,j)) * inv_2dx
            dvdy_ab = (uy(i,jp1)-uy(i,jm1)) * inv_2dy 
            
            dudy_ab = (ux(i,jp1)-ux(i,jm1)) * inv_2dy 
            dvdx_ab = (uy(ip1,j)-uy(im1,j)) * inv_2dx
            
            ! Loop over column
            do k = 1, nz 

                ! Get vertical shear strain rate terms
                duxdz_ab = duxdz(i,j,k) 
                duydz_ab = duydz(i,j,k) 
                
                ! Calculate the total effective strain rate from L19, Eq. 21 
                eps_sq_ab = dudx_ab**2 + dvdy_ab**2 + dudx_ab*dvdy_ab + 0.25_wp*(dudy_ab+dvdx_ab)**2 &
                            + 0.25_wp*duxdz_ab**2 + 0.25_wp*duydz_ab**2 + eps_0_sq
 
                ! Get the rate factor on ab-node too
                ATT_ab = 0.25_wp*(ATT(i,j,k)+ATT(ip1,j,k)+ATT(i,jp1,k)+ATT(ip1,jp1,k))

                ! Calculate effective viscosity on ab-nodes
                visc_eff(i,j,k) = 0.5_wp*(eps_sq_ab)**(p1) * ATT_ab**(p2)
                   
            end do 

        end do  
        end do 

        return 

    end subroutine calc_visc_eff_3D_ab

    subroutine calc_visc_eff_int(visc_eff_int,visc_eff,H_ice,zeta_aa,boundaries)

        implicit none 

        real(wp), intent(OUT) :: visc_eff_int(:,:)
        real(wp), intent(IN)  :: visc_eff(:,:,:)
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: zeta_aa(:)
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer :: i, j, nx, ny
        integer :: im1, ip1, jm1, jp1  
        real(wp) :: visc_eff_mean 
        real(wp) :: wt 

        nx = size(visc_eff_int,1)
        ny = size(visc_eff_int,2)

        do j = 1, ny 
        do i = 1, nx

            ! Calculate the vertically averaged viscosity for this point
            visc_eff_mean = integrate_trapezoid1D_pt(visc_eff(i,j,:),zeta_aa) 
            
            if (H_ice(i,j) .gt. 0.0) then 
                visc_eff_int(i,j) = visc_eff_mean*H_ice(i,j) 
            else
                visc_eff_int(i,j) = visc_eff_mean 
            end if 

        end do 
        end do 

        return

    end subroutine calc_visc_eff_int

    function calc_staggered_margin(var0,var1,f0,f1) result(var_mid)
        ! Calculate a staggered point but taking upstream point at the margin 

        implicit none 

        real(wp), intent(IN) :: var0 
        real(wp), intent(IN) :: var1 
        real(wp), intent(IN) :: f0 
        real(wp), intent(IN) :: f1
        real(wp) :: var_mid 

        if (f0 .eq. 1.0_wp .and. f1 .lt. 1.0_wp) then 
            ! At the margin 

            var_mid = var0 

        else if (f0 .lt. 1.0_wp .and. f1 .eq. 1.0_wp) then 
            ! At the margin 

            var_mid = var1 

        else 
            ! Inland ice, or ice free, simple average 

            var_mid = 0.5_wp*(var0+var1) 

        end if 

        return 

    end function calc_staggered_margin

    subroutine calc_F_integral(F_int,visc,H_ice,zeta_aa,n)
        ! Useful integrals, following Arthern et al. (2015) Eq. 7,
        ! and Lipscomb et al. (2019), Eq. 30
        ! F_n = int_zb_zs{ 1/visc * ((s-z)/H)**n dz}

        implicit none 

        real(wp), intent(OUT) :: F_int(:,:) 
        real(wp), intent(IN)  :: visc(:,:,:)
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: zeta_aa(:)
        real(wp), intent(IN)  :: n  

        ! Local variables 
        integer :: i, j, nx, ny, nz_aa, np
        real(wp) :: H_eff  

        nx    = size(visc,1)
        ny    = size(visc,2) 
        nz_aa = size(visc,3)

        ! Vertically integrate at each point
        do j = 1, ny 
        do i = 1, nx

            if (H_ice(i,j) .gt. 0.0) then 
                ! Viscosity should be nonzero here, perform integration 

                H_eff = H_ice(i,j)
                F_int(i,j) = integrate_trapezoid1D_pt( (H_eff/visc(i,j,:) )*(1.0_wp-zeta_aa)**n,zeta_aa)

            else 

                F_int(i,j) = 0.0_wp 

            end if 

        end do 
        end do 

        return

    end subroutine calc_F_integral
    
    subroutine calc_beta_eff(beta_eff,beta,F2,zeta_aa,no_slip)
        ! Calculate the depth-averaged horizontal velocity (ux_bar,uy_bar)

        ! Note: L19 staggers the F-integral F2, then solves for beta 

        implicit none 
        
        real(wp), intent(OUT) :: beta_eff(:,:)    ! ab-nodes
        real(wp), intent(IN)  :: beta(:,:)        ! ab-nodes
        real(wp), intent(IN)  :: F2(:,:)          ! ab-nodes
        real(wp), intent(IN)  :: zeta_aa(:)       ! aa-nodes
        logical,    intent(IN)  :: no_slip 

        ! Local variables 
        integer    :: i, j, nx, ny

        nx = size(beta_eff,1)
        ny = size(beta_eff,2)

        if (no_slip) then 
            ! No basal sliding allowed, impose beta_eff derived from viscosity 
            ! following L19, Eq. 35 (or G11, Eq. 42)

            beta_eff = 1.0_wp / F2 

        else 
            ! Basal sliding allowed, calculate beta_eff 
            ! following L19, Eq. 33 (or G11, Eq. 41)

            beta_eff = beta / (1.0_wp+beta*F2)

        end if 

        return 

    end subroutine calc_beta_eff

    subroutine calc_vel_basal_ab(ux_b,uy_b,ux_bar,uy_bar,F2,taubx,tauby,no_slip)
        ! Calculate basal sliding following Goldberg (2011), Eq. 34
        ! (or it can also be obtained from L19, Eq. 32 given ub*beta=taub)

        implicit none
        
        real(wp), intent(OUT) :: ux_b(:,:) 
        real(wp), intent(OUT) :: uy_b(:,:)
        real(wp), intent(IN)  :: ux_bar(:,:) 
        real(wp), intent(IN)  :: uy_bar(:,:)
        real(wp), intent(IN)  :: F2(:,:)
        real(wp), intent(IN)  :: taubx(:,:) 
        real(wp), intent(IN)  :: tauby(:,:)
        logical,    intent(IN)  :: no_slip

        ! Local variables 
        integer    :: i, j, nx, ny 
        integer    :: ip1, jp1 
        real(wp) :: F2_ac 

        nx = size(ux_b,1)
        ny = size(ux_b,2) 

        if (no_slip) then 
            ! Set basal velocity to zero 
            ! (this comes out naturally more or less with beta_eff set as above, 
            !  but ensuring basal velocity is zero adds stability)
            
            ux_b = 0.0_wp 
            uy_b = 0.0_wp 

        else 
            ! Calculate basal velocity normally 

            do j = 1, ny 
            do i = 1, nx 

                ! Calculate basal velocity component 
                ux_b(i,j) = ux_bar(i,j) - taubx(i,j)*F2(i,j) 
                uy_b(i,j) = uy_bar(i,j) - tauby(i,j)*F2(i,j) 

            end do 
            end do  

        end if 

        return
        
    end subroutine calc_vel_basal_ab

    subroutine calc_basal_stress_ab(taubx,tauby,beta_eff,ux_bar,uy_bar)
        ! Calculate the basal stress resulting from sliding (friction times velocity)
        ! Note: calculated on ac-nodes.
        ! taub [Pa] 
        ! beta [Pa a m-1]
        ! u    [m a-1]
        ! taub = beta*u (here defined with taub in the same direction as u)

        implicit none 

        real(wp), intent(OUT) :: taubx(:,:)             ! [Pa] Basal stress (ab nodes)
        real(wp), intent(OUT) :: tauby(:,:)             ! [Pa] Basal stress (ab nodes)
        real(wp), intent(IN)  :: beta_eff(:,:)          ! [Pa a m-1] Effective basal friction (ab nodes)
        real(wp), intent(IN)  :: ux_bar(:,:)            ! [m a-1] depth-ave velocity (ab nodes)
        real(wp), intent(IN)  :: uy_bar(:,:)            ! [m a-1] depth-ave velocity (ab nodes)
        
        ! Local variables 
        integer :: i, j, nx, ny 

        nx = size(beta_eff,1)
        ny = size(beta_eff,2) 

        do j = 1, ny 
        do i = 1, nx 

            taubx(i,j) = beta_eff(i,j) * ux_bar(i,j) 
            if (abs(taubx(i,j)) .lt. TOL_UNDERFLOW) taubx(i,j) = 0.0_wp 

            tauby(i,j) = beta_eff(i,j) * uy_bar(i,j) 
            if (abs(tauby(i,j)) .lt. TOL_UNDERFLOW) tauby(i,j) = 0.0_wp 

        end do 
        end do  

        return 

    end subroutine calc_basal_stress_ab

    elemental subroutine limit_vel(u,u_lim)
        ! Apply a velocity limit (for stability)

        implicit none 

        real(wp), intent(INOUT) :: u 
        real(wp), intent(IN)    :: u_lim

        real(wp), parameter :: tol = TOL_UNDERFLOW
        
        u = min(u, u_lim)
        u = max(u,-u_lim)

        ! Also avoid underflow errors 
        if (abs(u) .lt. tol) u = 0.0 

        return 

    end subroutine limit_vel
    
    subroutine calc_beta_aa_power_plastic_ab(beta,ux,uy,c_bed,q,u_0)
        ! Calculate basal friction coefficient (beta) that
        ! enters the SSA solver as a function of basal velocity
        ! using a power-law form following Bueler and van Pelt (2015)
        
        implicit none
        
        real(wp), intent(OUT) :: beta(:,:)          ! ab-nodes
        real(wp), intent(IN)  :: ux(:,:)            ! ab-nodes
        real(wp), intent(IN)  :: uy(:,:)            ! ab-nodes
        real(wp), intent(IN)  :: c_bed(:,:)         ! aa-nodes
        real(wp), intent(IN)  :: q
        real(wp), intent(IN)  :: u_0                ! [m/yr]  

        ! Local variables
        integer  :: i, j, nx, ny
        integer  :: im1, ip1, jm1, jp1 
        real(wp) :: uxy

        real(wp), parameter :: ub_min    = 1e-3_wp          ! [m/yr] Minimum velocity is positive small value to avoid divide by zero
        real(wp), parameter :: ub_sq_min = ub_min**2

        real(wp), allocatable :: c_bed_ab(:,:) 

        nx = size(beta,1)
        ny = size(beta,2)
        
        allocate(c_bed_ab(nx,ny))

        ! Get staggered c_bed 
        call map_a_to_b_2D(c_bed_ab,c_bed) 

        ! Initially set friction to zero everywhere
        beta = 0.0_wp 
        
        do j = 1, ny
        do i = 1, nx

            ! Get neighbor indices 
            im1 = max(i-1,1)
            ip1 = min(i+1,nx)
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny)
            

            if (q .eq. 1.0_wp) then 
                ! Linear law, no f(ub) term

                beta(i,j) = c_bed_ab(i,j) / u_0

            else
                ! Non-linear law with f(ub) term 

                uxy = sqrt(ux(i,j)**2 + uy(i,j)**2 + ub_sq_min)
                
                if (q .eq. 0.0_wp) then 
                    ! Plastic law

                    beta(i,j) = c_bed_ab(i,j) * (1.0_wp / uxy)
                    
                else

                    beta(i,j) = c_bed_ab(i,j) * (uxy / u_0)**q * (1.0_wp / uxy)
                
                end if 
                
            end if 

        end do
        end do

        return
        
    end subroutine calc_beta_aa_power_plastic_ab

end module velocity_diva_ab
