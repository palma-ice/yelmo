module velocity_l1l2

    use yelmo_defs ,only  : wp, prec, rho_ice, rho_sw, rho_w, g, tol_underflow, pi
    use yelmo_tools, only : stagger_aa_ab, stagger_aa_ab_ice, stagger_ab_aa_ice, & 
                    stagger_nodes_aa_ab_ice, stagger_nodes_acx_ab_ice, stagger_nodes_acy_ab_ice, &
                    staggerdiffx_nodes_aa_ab_ice, staggerdiffy_nodes_aa_ab_ice, &
                    staggerdiff_nodes_acx_ab_ice, staggerdiff_nodes_acy_ab_ice, &
                    staggerdiffcross_nodes_acx_ab_ice, staggerdiffcross_nodes_acy_ab_ice, &
                    staggerdiffcross_aa_acx_ice, staggerdiffcross_aa_acy_ice, &
                    integrate_trapezoid1D_1D, integrate_trapezoid1D_pt, minmax, &
                    calc_vertical_integrated_2D, calc_vertical_integrated_3D

    use basal_dragging 
    use solver_ssa_ac
    use velocity_general, only : set_inactive_margins, &
                        picard_calc_error, picard_calc_error_angle, picard_relax, &
                        picard_calc_convergence_l1rel_matrix, picard_calc_convergence_l2 
    
    implicit none 

    type l1l2_param_class

        character(len=256) :: ssa_lis_opt 
        character(len=256) :: ssa_lateral_bc
        character(len=256) :: boundaries 
        logical    :: no_slip 
        integer    :: visc_method
        real(prec) :: visc_const
        integer    :: beta_method
        real(prec) :: beta_const
        real(prec) :: beta_q                ! Friction law exponent
        real(prec) :: beta_u0               ! [m/a] Friction law velocity threshold 
        integer    :: beta_gl_scale         ! Beta grounding-line scaling method (beta => 0 at gl?)
        integer    :: beta_gl_stag          ! Beta grounding-line staggering method 
        real(prec) :: beta_gl_f             ! Fraction of beta at gl 
        real(prec) :: H_grnd_lim 
        real(prec) :: beta_min              ! Minimum allowed value of beta
        real(prec) :: eps_0 
        real(prec) :: ssa_vel_max
        integer    :: ssa_iter_max 
        real(prec) :: ssa_iter_rel 
        real(prec) :: ssa_iter_conv 
        logical    :: ssa_write_log 

    end type

    private
    public :: l1l2_param_class 
    public :: calc_velocity_l1l2

contains 
    
    subroutine calc_velocity_l1l2(ux,uy,ux_bar,uy_bar,ux_b,uy_b,ux_i,uy_i,taub_acx,taub_acy, &
                                  beta,beta_acx,beta_acy,visc_eff,visc_eff_int,  &
                                  ssa_mask_acx,ssa_mask_acy,ssa_err_acx,ssa_err_acy,ssa_iter_now, &
                                  c_bed,taud_acx,taud_acy,taul_int_acx,taul_int_acy, &
                                  H_ice,f_ice,H_grnd,f_grnd, &
                                  f_grnd_acx,f_grnd_acy,mask_frnt,ATT,zeta_aa,zeta_ac,z_sl,z_bed,z_srf,dx,dy,n_glen,par)
        ! This subroutine is used to solve the horizontal velocity system (ux,uy)
        ! following the L1L2 solver formulation (ie, depth-integrated solver
        ! that reduces to SIA in the limit of zero sliding), as outlined 
        ! by Perego et al. (2012) and following the implementation
        ! in CISMv2.

        implicit none 

        real(prec), intent(INOUT) :: ux(:,:,:)          ! [m/a]
        real(prec), intent(INOUT) :: uy(:,:,:)          ! [m/a]
        real(prec), intent(INOUT) :: ux_bar(:,:)        ! [m/a]
        real(prec), intent(INOUT) :: uy_bar(:,:)        ! [m/a]
        real(prec), intent(INOUT) :: ux_b(:,:)          ! [m/a]
        real(prec), intent(INOUT) :: uy_b(:,:)          ! [m/a]
        real(prec), intent(INOUT) :: ux_i(:,:,:)        ! [m/a]
        real(prec), intent(INOUT) :: uy_i(:,:,:)        ! [m/a]
        real(prec), intent(INOUT) :: taub_acx(:,:)      ! [Pa]
        real(prec), intent(INOUT) :: taub_acy(:,:)      ! [Pa]
        real(prec), intent(INOUT) :: beta(:,:)          ! [Pa a/m]
        real(prec), intent(INOUT) :: beta_acx(:,:)      ! [Pa a/m]
        real(prec), intent(INOUT) :: beta_acy(:,:)      ! [Pa a/m]
        real(prec), intent(INOUT) :: visc_eff(:,:,:)    ! [Pa a]
        real(prec), intent(OUT)   :: visc_eff_int(:,:)  ! [Pa a m]
        integer,    intent(OUT)   :: ssa_mask_acx(:,:)  ! [-]
        integer,    intent(OUT)   :: ssa_mask_acy(:,:)  ! [-]
        real(prec), intent(OUT)   :: ssa_err_acx(:,:)
        real(prec), intent(OUT)   :: ssa_err_acy(:,:)
        integer,    intent(OUT)   :: ssa_iter_now 
        real(prec), intent(IN)    :: c_bed(:,:)         ! [Pa]
        real(prec), intent(IN)    :: taud_acx(:,:)      ! [Pa]
        real(prec), intent(IN)    :: taud_acy(:,:)      ! [Pa]
        real(wp),   intent(IN)    :: taul_int_acx(:,:)  ! [Pa m]
        real(wp),   intent(IN)    :: taul_int_acy(:,:)  ! [Pa m]
        real(prec), intent(IN)    :: H_ice(:,:)         ! [m]
        real(prec), intent(IN)    :: f_ice(:,:)         ! [--]
        real(prec), intent(IN)    :: H_grnd(:,:)        ! [m]
        real(prec), intent(IN)    :: f_grnd(:,:)        ! [-]
        real(prec), intent(IN)    :: f_grnd_acx(:,:)    ! [-]
        real(prec), intent(IN)    :: f_grnd_acy(:,:)    ! [-]
        integer,    intent(IN)    :: mask_frnt(:,:)     ! [-]
        real(prec), intent(IN)    :: ATT(:,:,:)         ! [a^-1 Pa^-n_glen]
        real(prec), intent(IN)    :: zeta_aa(:)         ! [-]
        real(prec), intent(IN)    :: zeta_ac(:)         ! [-]
        real(prec), intent(IN)    :: z_sl(:,:)          ! [m]
        real(prec), intent(IN)    :: z_bed(:,:)         ! [m]
        real(prec), intent(IN)    :: z_srf(:,:)         ! [m]
        real(prec), intent(IN)    :: dx                 ! [m]
        real(prec), intent(IN)    :: dy                 ! [m]
        real(prec), intent(IN)    :: n_glen 
        type(l1l2_param_class), intent(IN) :: par       ! List of parameters that should be defined

        ! Local variables 
        integer :: i, j, k, nx, ny, nz_aa, nz_ac, iter 
        logical :: is_converged

        real(prec), allocatable :: ux_b_nm1(:,:) 
        real(prec), allocatable :: uy_b_nm1(:,:)  
        integer,    allocatable :: ssa_mask_acx_ref(:,:)
        integer,    allocatable :: ssa_mask_acy_ref(:,:)

        real(prec), allocatable :: visc_eff_ab(:,:,:) 

        real(wp), allocatable :: corr_nm1(:) 
        real(wp), allocatable :: corr_nm2(:) 
        
        real(wp) :: corr_theta
        real(wp) :: corr_rel 
        real(wp) :: L2_norm 
        real(wp) :: ssa_resid 

        integer :: ij(2) 

        logical, parameter :: write_ssa_diagnostics      = .FALSE. 
        logical, parameter :: write_ssa_diagnostics_stop = .FALSE.   ! Stop simulation after completing iterations?

        type(linear_solver_class) :: lgs_prev 
        type(linear_solver_class) :: lgs_now
        
        nx    = size(ux,1)
        ny    = size(ux,2)
        nz_aa = size(ux,3)

        ! Prepare local variables 
        allocate(ux_b_nm1(nx,ny))
        allocate(uy_b_nm1(nx,ny))

        allocate(ssa_mask_acx_ref(nx,ny))
        allocate(ssa_mask_acy_ref(nx,ny))

        allocate(visc_eff_ab(nx,ny,nz_aa)) 

        allocate(corr_nm1(2*nx*ny))
        allocate(corr_nm2(2*nx*ny))

        ! Store original ssa mask before iterations
        ssa_mask_acx_ref = ssa_mask_acx
        ssa_mask_acy_ref = ssa_mask_acy
            
        ! Initially set error very high 
        ssa_err_acx = 1.0_prec 
        ssa_err_acy = 1.0_prec 
        
        corr_nm1 = 0.0_wp 
        corr_nm2 = 0.0_wp 

        ! Ensure dynamically inactive cells have no velocity at 
        ! outer margins before starting iterations
        call set_inactive_margins(ux_b,uy_b,f_ice,par%boundaries)

        if (write_ssa_diagnostics) then 
            call ssa_diagnostics_write_init("yelmo_ssa.nc",nx,ny,time_init=1.0_wp)
        end if 

        ! Initialize linear solver variables for current and previous iteration
        call linear_solver_init(lgs_now,nx,ny)
        lgs_prev = lgs_now 

        do iter = 1, par%ssa_iter_max 

            ! Store solution from previous iteration (nm1 == n minus 1) 
            ux_b_nm1 = ux_b 
            uy_b_nm1 = uy_b 
            
            ! =========================================================================================
            ! Step 1: Calculate fields needed by ssa solver (visc_eff_int, beta)

            ! Calculate 3D effective viscosity (and staggered onto ab-nodes)
            select case(par%visc_method)

                case(0)
                    ! Impose constant viscosity value 

                    visc_eff    = par%visc_const 
                    visc_eff_ab = par%visc_const 

                case(1) 
                    ! Calculate 3D effective viscosity, using velocity solution from previous iteration
                    
                    call calc_visc_eff_3D(visc_eff,ux_b,uy_b,taud_acx,taud_acy,ATT, &
                                                H_ice,f_ice,zeta_aa,dx,dy,n_glen,par%eps_0,par%boundaries)

                    do k = 1, nz_aa 
                        visc_eff_ab(:,:,k) = stagger_aa_ab_ice(visc_eff(:,:,k),H_ice,f_ice)
                    end do 

                    ! call calc_visc_eff_3D_1(visc_eff,visc_eff_ab,ux_b,uy_b,taud_acx,taud_acy,ATT, &
                    !                             H_ice,f_ice,zeta_aa,dx,dy,n_glen,par%eps_0,par%boundaries)

                case DEFAULT 

                    write(*,*) "calc_velocity_l1l2:: Error: visc_method not recognized."
                    write(*,*) "visc_method = ", par%visc_method 
                    stop 
                    
            end select
                  
            ! Calculate depth-integrated effective viscosity
            ! Note L19 uses eta_bar*H in the ssa equation. Yelmo uses eta_int=eta_bar*H directly.
            call calc_visc_eff_int(visc_eff_int,visc_eff,H_ice,f_ice,zeta_aa,par%boundaries)
            
            ! Smooth the viscosity at the ice margins if it is too low
            ! to avoid singularities (mainly for EISMINT/dome experiments)
            !call smooth_visc_eff_int_margin(visc_eff_int,H_ice)

            ! Calculate beta (at the ice base)
            call calc_beta(beta,c_bed,ux_b,uy_b,H_ice,f_ice,H_grnd,f_grnd,z_bed,z_sl,par%beta_method, &
                                par%beta_const,par%beta_q,par%beta_u0,par%beta_gl_scale,par%beta_gl_f, &
                                par%H_grnd_lim,par%beta_min,par%boundaries)

            ! Stagger beta and beta_eff 
            call stagger_beta(beta_acx,beta_acy,beta,H_ice,f_ice,ux_b,uy_b, &
                        f_grnd,f_grnd_acx,f_grnd_acy,par%beta_gl_stag,par%beta_min,par%boundaries)

            ! =========================================================================================
            ! Step 2: determine the basal velocity ux_b/uy_b 

            if (par%no_slip) then 
                ! Simply set ux_b/uy_b equal to zero, as no sliding is allowed 

                ux_b = 0.0_prec 
                uy_b = 0.0_prec 

            else 
                ! Call the SSA solver to obtain new estimate of ux_b/uy_b

if (.FALSE.) then 
                if (iter .gt. 1) then
                    ! Update ssa mask based on convergence with previous step to reduce area being solved 
                    call update_ssa_mask_convergence(ssa_mask_acx,ssa_mask_acy,ssa_err_acx,ssa_err_acy,err_lim=real(1e-5,prec))
                    !call update_ssa_mask_convergence(ssa_mask_acx,ssa_mask_acy,ssa_err_acx,ssa_err_acy,err_lim=par%ssa_iter_conv*1e-2)  
                end if 
end if 
                
                ! Populate ssa matrices Ax=b
                call linear_solver_matrix_ssa_ac_csr_2D(lgs_now,ux_b,uy_b,beta_acx,beta_acy,visc_eff_int,  &
                                    ssa_mask_acx,ssa_mask_acy,mask_frnt,H_ice,f_ice,taud_acx,taud_acy, &
                                    taul_int_acx,taul_int_acy,dx,dy,par%boundaries,par%ssa_lateral_bc)

                ! Solve linear equation
                call linear_solver_matrix_solve(lgs_now,par%ssa_lis_opt)
                
                ! Save L2_norm locally
                L2_norm = lgs_now%L2_rel_norm 

                ! Store velocity solution
                call linear_solver_save_velocity(ux_b,uy_b,lgs_now,par%ssa_vel_max)


! ajr: note that for MISMIP3D experiments to converge well,
! a value of par%ssa_iter_conv <= 1e-3 is needed if also
! using the adaptive relaxation scheme with corr_theta below.
! If using a constant relaxation of par%ssa_iter_rel=0.7,
! then par%ssa_iter_conv = 1e-2 is sufficient for proper performance.
! For Antarctica, the adaptive method can give some strange
! convergence issues. It has been disabled for now (2022-02-09).
if (.FALSE.) then
                ! Calculate errors 
                corr_nm2 = corr_nm1 
                call picard_calc_error(corr_nm1,ux_b,uy_b,ux_b_nm1,uy_b_nm1)

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

                ! Apply relaxation to keep things stable
                !call relax_ssa(ux_b,uy_b,ux_b_nm1,uy_b_nm1,rel=par%ssa_iter_rel)
                call picard_relax(ux_b,uy_b,ux_b_nm1,uy_b_nm1,rel=corr_rel)
                
            end if 

            ! Check for convergence
            ! is_converged = check_vel_convergence_l2rel(ux_b,uy_b,ux_b_nm1,uy_b_nm1,ssa_mask_acx.gt.0,     &
            !                                            ssa_mask_acy.gt.0,par%ssa_iter_conv,iter,par%ssa_iter_max, &
            !                                            par%ssa_write_log,use_L2_norm=.FALSE.,L2_norm=L2_norm)
            call picard_calc_convergence_l2(is_converged,ssa_resid,ux_b,uy_b,ux_b_nm1,uy_b_nm1, &
                                                ssa_mask_acx.gt.0,ssa_mask_acy.gt.0,par%ssa_iter_conv,  &
                                                iter,par%ssa_iter_max,par%ssa_write_log)

            ! Calculate an L1 error metric over matrix for diagnostics
            call picard_calc_convergence_l1rel_matrix(ssa_err_acx,ssa_err_acy,ux_b,uy_b,ux_b_nm1,uy_b_nm1)

            ! Store current total iterations for output
            ssa_iter_now = iter 

            if (write_ssa_diagnostics) then  
                call ssa_diagnostics_write_step("yelmo_ssa.nc",ux_b,uy_b,L2_norm,beta_acx,beta_acy,visc_eff_int, &
                                        ssa_mask_acx,ssa_mask_acy,ssa_err_acx,ssa_err_acy,H_ice,f_ice,taud_acx,taud_acy, &
                                        taul_int_acx,taul_int_acy,H_grnd,z_sl,z_bed,z_srf,ux_b_nm1,uy_b_nm1,time=real(iter,wp))    
            end if 

            ! =========================================================================================
            ! Update additional fields based on output of solver
            
            ! Calculate basal stress 
            call calc_basal_stress(taub_acx,taub_acy,beta_acx,beta_acy,ux_b,uy_b)

            ! Exit iterations if ssa solution has converged
            if (is_converged) exit 
            
        end do 

        ! Iterations are finished, finalize calculations of 3D velocity field 

        if (write_ssa_diagnostics .and. write_ssa_diagnostics_stop) then 
            stop 
        end if 

        ! Calculate the 3D horizontal velocity field
        call calc_vel_horizontal_3D(ux,uy,ux_b,uy_b,taud_acx,taud_acy,visc_eff,ATT,H_ice, &
                                    f_ice,zeta_aa,dx,dy,n_glen,par%eps_0,par%boundaries)
        ! call calc_vel_horizontal_3D_0(ux,uy,ux_b,uy_b,taud_acx,taud_acy,visc_eff,visc_eff_ab,ATT,H_ice, &
        !                             f_ice,zeta_aa,zeta_ac,dx,dy,n_glen,par%eps_0,par%boundaries)

        ! Limit the velocity generally =====================
        call limit_vel(ux,par%ssa_vel_max)
        call limit_vel(uy,par%ssa_vel_max)
        
        ! Calculate depth-averaged horizontal velocity 
        ux_bar = calc_vertical_integrated_2D(ux,zeta_aa)
        uy_bar = calc_vertical_integrated_2D(uy,zeta_aa)
        
        ! Also calculate the shearing contribution
        do k = 1, nz_aa 
            ux_i(:,:,k) = ux(:,:,k) - ux_b 
            uy_i(:,:,k) = uy(:,:,k) - uy_b 
        end do

        return 

    end subroutine calc_velocity_l1l2

    subroutine calc_vel_horizontal_3D(ux,uy,ux_b,uy_b,taud_acx,taud_acy, &
                        visc_eff,ATT,H_ice,f_ice,zeta_aa, &
                        dx,dy,n_glen,eps_0,boundaries)
        ! Caluculate the 3D horizontal velocity field (ux,uy)
        ! for the L1L2 solver following Perego et al. (2012)
        ! and the blueprint by Lipscomb et al. (2019) in CISM

        implicit none 

        real(wp), intent(OUT) :: ux(:,:,:) 
        real(wp), intent(OUT) :: uy(:,:,:) 
        real(wp), intent(IN)  :: ux_b(:,:) 
        real(wp), intent(IN)  :: uy_b(:,:) 
        real(wp), intent(IN)  :: taud_acx(:,:) 
        real(wp), intent(IN)  :: taud_acy(:,:)
        real(wp), intent(IN)  :: visc_eff(:,:,:)      ! on aa-nodes 
        real(wp), intent(IN)  :: ATT(:,:,:)  
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: zeta_aa(:) 
        real(wp), intent(IN)  :: dx
        real(wp), intent(IN)  :: dy
        real(wp), intent(IN)  :: n_glen   
        real(wp), intent(IN)  :: eps_0                ! [1/a] Regularization constant (minimum strain rate, ~1e-8)
        character(len=*), intent(IN) :: boundaries 

        ! Local variables
        integer  :: i, j, k, nx, ny, nz_aa, nz_ac   
        integer  :: ip1, jp1, im1, jm1 
        real(wp) :: inv_4dx, inv_4dy 
        real(wp) :: zeta_ac1, zeta_ac0 
        real(wp) :: dw1dx, dw2dy, dw3dx, dw4dy
        real(wp) :: depth_ab 
        real(wp) :: fact_ac 

        real(wp) :: dudx_ab4(4) 
        real(wp) :: dvdy_ab4(4) 
        real(wp) :: dudy_ab4(4)
        real(wp) :: dvdx_ab4(4)
        real(wp) :: eps_par_sq4(4)
        real(wp) :: eps_par4(4)
        real(wp) :: visc_eff_ab4(4)
        real(wp) :: tau_par_ab4_up(4)
        real(wp) :: tau_par_ab4_dn(4)
        real(wp) :: tau_par_ab4(4)

        real(wp) :: dw1dx_ab(4)
        real(wp) :: dw2dy_ab(4)
        real(wp) :: dw3dx_ab(4)
        real(wp) :: dw4dy_ab(4) 

        real(wp) :: H_ice_ab4(4)
        real(wp) :: ATT_ab4_up(4) 
        real(wp) :: ATT_ab4_dn(4) 
        
        real(wp) :: tau_xz_ab4_up(4)
        real(wp) :: tau_xz_ab4_dn(4)
        real(wp) :: tau_xz_ab4(4)
        real(wp) :: tau_yz_ab4_up(4)
        real(wp) :: tau_yz_ab4_dn(4)
        real(wp) :: tau_yz_ab4(4)
        real(wp) :: tau_eff_sq_ab4(4)
        real(wp) :: ATT_ab4(4)
        real(wp) :: fact_ab4(4)

        real(wp) :: wt_ab(4) 
        real(wp) :: wt 

        real(wp), allocatable :: visc_eff_int3D(:,:,:) 
        real(wp), allocatable :: dudx_aa(:,:)
        real(wp), allocatable :: dvdy_aa(:,:)
        real(wp), allocatable :: dudy_aa(:,:)
        real(wp), allocatable :: dvdx_aa(:,:)
        real(wp), allocatable :: tau_par(:,:,:)

        real(wp), allocatable :: work1_aa(:,:)
        real(wp), allocatable :: work2_aa(:,:)
        real(wp), allocatable :: work3_aa(:,:)
        real(wp), allocatable :: work4_aa(:,:)

        real(wp), allocatable :: dw1dx_aa(:,:)
        real(wp), allocatable :: dw2dy_aa(:,:)
        real(wp), allocatable :: dw3dx_aa(:,:)
        real(wp), allocatable :: dw4dy_aa(:,:)
        
        real(wp), allocatable :: tau_xz(:,:,:) 
        real(wp), allocatable :: tau_yz(:,:,:) 
        real(wp), allocatable :: fact_ab(:,:) 

        real(wp) :: p1, eps_0_sq 
        real(wp) :: dzeta 

        nx    = size(ux,1)
        ny    = size(ux,2) 
        nz_aa = size(ux,3) 
        nz_ac = nz_aa + 1

        ! Allocate local arrays 
        allocate(visc_eff_int3D(nx,ny,nz_aa)) 
        allocate(dudx_aa(nx,ny))
        allocate(dvdy_aa(nx,ny))
        allocate(dudy_aa(nx,ny))
        allocate(dvdx_aa(nx,ny))
        allocate(tau_par(nx,ny,nz_aa))

        allocate(work1_aa(nx,ny))
        allocate(work2_aa(nx,ny))
        allocate(work3_aa(nx,ny))
        allocate(work4_aa(nx,ny))

        allocate(dw1dx_aa(nx,ny))
        allocate(dw2dy_aa(nx,ny))
        allocate(dw3dx_aa(nx,ny))
        allocate(dw4dy_aa(nx,ny))

        allocate(tau_xz(nx,ny,nz_aa))
        allocate(tau_yz(nx,ny,nz_aa))
        allocate(fact_ab(nx,ny))

        ! Calculate scaling factors
        inv_4dx = 1.0_prec / (4.0_prec*dx) 
        inv_4dy = 1.0_prec / (4.0_prec*dy) 

        ! Calculate an exponent 
        p1 = (n_glen-1.0_wp)/2.0_wp

        ! Calculate squared minimum strain rate 
        eps_0_sq = eps_0*eps_0 

        wt_ab = 1.0 
        wt = sum(wt_ab)
        wt_ab = wt_ab / wt 

        ! Step 0: 
        ! Compute the integral of visc_eff from each layer to the surface (P12, Eq. 28)
        do j = 1, ny 
        do i = 1, nx 

            ! Start at the surface
            zeta_ac1 = zeta_aa(nz_aa)
            zeta_ac0 = 0.5_wp*(zeta_aa(nz_aa)+zeta_aa(nz_aa-1))
            visc_eff_int3D(i,j,nz_aa) = visc_eff(i,j,nz_aa) * (zeta_ac1-zeta_ac0)*H_ice(i,j)

            ! Integrate down to near the base 
            do k = nz_aa-1, 2, -1 
                zeta_ac1 = 0.5_wp*(zeta_aa(k+1)+zeta_aa(k))
                zeta_ac0 = 0.5_wp*(zeta_aa(k)+zeta_aa(k-1))
                visc_eff_int3D(i,j,k) = visc_eff_int3D(i,j,k+1) &
                                        + visc_eff(i,j,k) * (zeta_ac1-zeta_ac0)*H_ice(i,j)
            end do 
            
            ! Get basal value
            zeta_ac1 = 0.5_wp*(zeta_aa(2)+zeta_aa(1))
            zeta_ac0 = zeta_aa(1)
            visc_eff_int3D(i,j,1) = visc_eff_int3D(i,j,2) &
                                        + visc_eff(i,j,1) * (zeta_ac1-zeta_ac0)*H_ice(i,j)
        
        end do 
        end do  

        ! Step 1: compute basal strain rates on ab-nodes and viscosity       
        do j = 1, ny 
        do i = 1, nx 

            im1 = max(i-1,1) 
            ip1 = min(i+1,nx) 
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny) 

            ! Calculate effective strain components from horizontal stretching on ab-nodes
            call staggerdiff_nodes_acx_ab_ice(dudx_ab4,ux_b,f_ice,i,j,dx)
            call staggerdiff_nodes_acy_ab_ice(dvdy_ab4,uy_b,f_ice,i,j,dy)

            ! Calculate cross terms on ab-nodes
            call staggerdiffcross_nodes_acx_ab_ice(dudy_ab4,ux_b,f_ice,i,j,dy)
            call staggerdiffcross_nodes_acy_ab_ice(dvdx_ab4,uy_b,f_ice,i,j,dx)
            
            ! Store on aa-nodes for later use 
            dudx_aa(i,j) = sum(dudx_ab4*wt_ab)
            dvdy_aa(i,j) = sum(dvdy_ab4*wt_ab)
            dudy_aa(i,j) = sum(dudy_ab4*wt_ab)
            dvdx_aa(i,j) = sum(dvdx_ab4*wt_ab)
            

            ! Calculate the 'parallel' effective strain rate from P12, Eq. 17
            eps_par_sq4 = dudx_ab4**2 + dvdy_ab4**2 + dudx_ab4*dvdy_ab4 &
                        + 0.25_wp*(dudy_ab4+dvdx_ab4)**2 + eps_0_sq
            eps_par4    = sqrt(eps_par_sq4) 

            ! Compute the 'parallel' shear stress for each layer (tau_parallel)
            do k = 1, nz_aa 

                ! call stagger_nodes_aa_ab_ice(visc_eff_ab4,visc_eff(:,:,k),f_ice,i,j)
                ! tau_par_ab4 = 2.0_wp * visc_eff_ab4 * eps_par4
                ! tau_par(i,j,k) = sum(tau_par_ab4*wt_ab)

                tau_par(i,j,k) = 2.0_wp * visc_eff(i,j,k) * sum(eps_par4*wt_ab)

            end do 

            ! Note: above, eps_par and thus tau_par should be zero when no_slip=True 
            ! and the basal velocity components are zero (effectively zero because eps_0 ensures nonzero value). 

        end do 
        end do

        ! Loop over layers 
        do k = 1, nz_aa

            ! Calculate working arrays for this layer (terms in parentheses in
            ! Perego et al, 2012, Eq. 27)
            work1_aa = 2.0_wp*visc_eff_int3D(:,:,k) * (2.d0*dudx_aa + dvdy_aa) 
            work2_aa = 2.0_wp*visc_eff_int3D(:,:,k) * 0.5*(dudy_aa+dvdx_aa)
            work3_aa = 2.0_wp*visc_eff_int3D(:,:,k) * 0.5*(dudy_aa+dvdx_aa)
            work4_aa = 2.0_wp*visc_eff_int3D(:,:,k) * (dudx_aa + 2.d0*dvdy_aa)  
            
            ! Calculate derivatives of work arrays, first on ab-nodes
            ! then averaged to aa-nodes. This helps with solver stability
            ! (eg, to avoid singularities in vel solution for ISMIPHOM-C)
            do j = 1, ny 
            do i = 1, nx 
                 
                ! Get derivatives on ab-nodes 
                call staggerdiffx_nodes_aa_ab_ice(dw1dx_ab,work1_aa,f_ice,i,j,dx)
                call staggerdiffy_nodes_aa_ab_ice(dw2dy_ab,work2_aa,f_ice,i,j,dy)
                call staggerdiffx_nodes_aa_ab_ice(dw3dx_ab,work3_aa,f_ice,i,j,dx)
                call staggerdiffy_nodes_aa_ab_ice(dw4dy_ab,work4_aa,f_ice,i,j,dy)
                
                ! Get derivatives on aa-nodes 
                dw1dx_aa(i,j) = sum(dw1dx_ab*wt_ab)
                dw2dy_aa(i,j) = sum(dw2dy_ab*wt_ab)
                dw3dx_aa(i,j) = sum(dw3dx_ab*wt_ab)
                dw4dy_aa(i,j) = sum(dw4dy_ab*wt_ab)
                
            end do 
            end do 

            ! Calcaluate vertical shear stress
            do j = 1, ny 
            do i = 1, nx 
                
                ! Stagger to acx-nodes 
                dw1dx = 0.5_wp*(dw1dx_aa(i,j)+dw1dx_aa(ip1,j))
                dw2dy = 0.5_wp*(dw2dy_aa(i,j)+dw2dy_aa(ip1,j))
                
                ! Stagger to acy-nodes 
                dw1dx = 0.5_wp*(dw3dx_aa(i,j)+dw3dx_aa(i,jp1))
                dw2dy = 0.5_wp*(dw4dy_aa(i,j)+dw4dy_aa(i,jp1))
                
                ! Calculate shear stress on ac-nodes
                tau_xz(i,j,k) = -(1.0_wp-zeta_aa(k))*taud_acx(i,j) + dw1dx + dw2dy
                tau_yz(i,j,k) = -(1.0_wp-zeta_aa(k))*taud_acy(i,j) + dw3dx + dw4dy

            end do 
            end do 

        end do 

        ux = 0.0_wp 
        uy = 0.0_wp

        ! Assign basal velocity value 
        ux(:,:,1)    = ux_b 
        uy(:,:,1)    = uy_b 
        fact_ab(:,:) = 0.0_prec 

        ! Loop over layers starting from first layer above the base to surface 
        do k = 2, nz_aa

            dzeta = zeta_aa(k) - zeta_aa(k-1) 

            ! Calculate tau_perp, tau_eff and factor to calculate velocities,
            ! all on ab-nodes 
            do i = 1, nx 
            do j = 1, ny 

                im1 = max(i-1,1) 
                ip1 = min(i+1,nx) 
                jm1 = max(j-1,1) 
                jp1 = min(j+1,ny) 

                ! Calculate effective stress on horizontal ab-nodes and vertical ac-node
                call stagger_nodes_acx_ab_ice(tau_xz_ab4_up,tau_xz(:,:,k),  f_ice,i,j)
                call stagger_nodes_acx_ab_ice(tau_xz_ab4_dn,tau_xz(:,:,k-1),f_ice,i,j)
                tau_xz_ab4 = 0.5_wp*(tau_xz_ab4_up+tau_xz_ab4_dn)

                call stagger_nodes_acy_ab_ice(tau_yz_ab4_up,tau_yz(:,:,k),  f_ice,i,j)
                call stagger_nodes_acy_ab_ice(tau_yz_ab4_dn,tau_yz(:,:,k-1),f_ice,i,j)
                tau_yz_ab4 = 0.5_wp*(tau_yz_ab4_up+tau_yz_ab4_dn)
                
                call stagger_nodes_aa_ab_ice(tau_par_ab4_up,tau_par(:,:,k),  f_ice,i,j)
                call stagger_nodes_aa_ab_ice(tau_par_ab4_dn,tau_par(:,:,k-1),f_ice,i,j)
                tau_par_ab4 = 0.5_wp*(tau_par_ab4_up+tau_par_ab4_dn) 

                tau_eff_sq_ab4 = tau_xz_ab4**2 + tau_yz_ab4**2 + tau_par_ab4**2

                ! Calculate factor to get velocity components
                call stagger_nodes_aa_ab_ice(ATT_ab4_up,ATT(:,:,k),f_ice,i,j)
                call stagger_nodes_aa_ab_ice(ATT_ab4_dn,ATT(:,:,k-1),f_ice,i,j)
                ATT_ab4 = 0.5_wp*(ATT_ab4_up+ATT_ab4_dn)

                ! Get ice thickness
                call stagger_nodes_aa_ab_ice(H_ice_ab4,H_ice,f_ice,i,j)

                ! Calculate multiplicative factor on ab-nodes
                if (p1 .ne. 0.0_wp) then 
                    fact_ab4 = 2.0_prec * ATT_ab4 * (dzeta*H_ice_ab4) * tau_eff_sq_ab4**p1
                else
                    fact_ab4 = 2.0_prec * ATT_ab4 * (dzeta*H_ice_ab4)
                end if 

                ! Calculate 3D horizontal velocity components on acx/acy nodes

                ! stagger factor to acx-nodes to calculate velocity
                if (f_ice(i,j) .eq. 1.0 .or. f_ice(ip1,j) .eq. 1.0) then 
                    fact_ac   = 0.5_prec*(fact_ab4(1)+fact_ab4(4))
                    ux(i,j,k) = ux(i,j,k-1) &
                                + fact_ac*0.5_wp*(tau_xz(i,j,k)+tau_xz(i,j,k-1))
                end if 

                ! stagger factor to acy-nodes to calculate velocity
                if (f_ice(i,j) .eq. 1.0 .or. f_ice(i,jp1) .eq. 1.0) then
                    fact_ac   = 0.5_prec*(fact_ab4(1)+fact_ab4(2)) 
                    uy(i,j,k) = uy(i,j,k-1) &
                                + fact_ac*0.5_wp*(tau_yz(i,j,k)+tau_yz(i,j,k-1))
                end if 

            end do 
            end do 
               
        end do 

        ! Apply boundary conditions as needed 
        if (trim(boundaries) .eq. "periodic") then

            ux(1,:,:)    = ux(nx-2,:,:) 
            ux(nx-1,:,:) = ux(2,:,:) 
            ux(nx,:,:)   = ux(3,:,:) 
            ux(:,1,:)    = ux(:,ny-1,:)
            ux(:,ny,:)   = ux(:,2,:) 

            uy(1,:,:)    = uy(nx-1,:,:) 
            uy(nx,:,:)   = uy(2,:,:) 
            uy(:,1,:)    = uy(:,ny-2,:)
            uy(:,ny-1,:) = uy(:,2,:) 
            uy(:,ny,:)   = uy(:,3,:)

        else if (trim(boundaries) .eq. "periodic-x") then 
            
            ux(1,:,:)    = ux(nx-2,:,:) 
            ux(nx-1,:,:) = ux(2,:,:) 
            ux(nx,:,:)   = ux(3,:,:) 
            ux(:,1,:)    = ux(:,2,:)
            ux(:,ny,:)   = ux(:,ny-1,:) 

            uy(1,:,:)    = uy(nx-1,:,:) 
            uy(nx,:,:)   = uy(2,:,:) 
            uy(:,1,:)    = uy(:,2,:)
            uy(:,ny-1,:) = uy(:,ny-2,:) 
            uy(:,ny,:)   = uy(:,ny-1,:)

        else if (trim(boundaries) .eq. "infinite") then 
            
            ux(1,:,:)    = ux(2,:,:) 
            ux(nx-1,:,:) = ux(nx-2,:,:) 
            ux(nx,:,:)   = ux(nx-1,:,:) 
            ux(:,1,:)    = ux(:,2,:)
            ux(:,ny,:)   = ux(:,ny-1,:) 

            uy(1,:,:)    = uy(2,:,:) 
            uy(nx,:,:)   = uy(nx-1,:,:) 
            uy(:,1,:)    = uy(:,2,:)
            uy(:,ny-1,:) = uy(:,ny-2,:) 
            uy(:,ny,:)   = uy(:,ny-1,:)

        end if 

        return 

    end subroutine calc_vel_horizontal_3D

    subroutine calc_visc_eff_3D(visc_eff,ux_b,uy_b,taud_acx,taud_acy,ATT,H_ice,f_ice,zeta_aa,dx,dy,n_glen,eps_0,boundaries)
        ! Caluculate the 3D effective viscosity field
        ! for the L1L2 solver following Perego et al. (2012)
        ! and the blueprint by Lipscomb et al. (2019) in CISM

        implicit none 

        real(prec), intent(OUT) :: visc_eff(:,:,:)        
        real(prec), intent(IN)  :: ux_b(:,:) 
        real(prec), intent(IN)  :: uy_b(:,:) 
        real(prec), intent(IN)  :: taud_acx(:,:) 
        real(prec), intent(IN)  :: taud_acy(:,:)
        real(prec), intent(IN)  :: ATT(:,:,:)  
        real(prec), intent(IN)  :: H_ice(:,:)
        real(prec), intent(IN)  :: f_ice(:,:)
        real(prec), intent(IN)  :: zeta_aa(:) 
        real(prec), intent(IN)  :: dx
        real(prec), intent(IN)  :: dy
        real(prec), intent(IN)  :: n_glen   
        real(prec), intent(IN)  :: eps_0                ! [1/a] Regularization constant (minimum strain rate, ~1e-8)
        character(len=*), intent(IN) :: boundaries 

        ! Local variables
        integer :: i, j, k, nx, ny, nz_aa  
        integer    :: ip1, jp1, im1, jm1 
        real(prec) :: inv_4dx, inv_4dy  
         
        real(prec) :: eps_0_sq 
        real(prec) :: a, b, c, rootA, rootB, np  
        real(prec) :: wt 
        integer    :: q 

        real(wp) :: dudx_ab(4)
        real(wp) :: dvdy_ab(4)
        real(wp) :: dudy_ab(4)
        real(wp) :: dvdx_ab(4) 
        real(wp) :: ATT_ab(4) 
        real(wp) :: eps_par_sq(4)
        real(wp) :: eps_par_ab(4)
        real(wp) :: taudx_ab(4) 
        real(wp) :: taudy_ab(4) 
        real(wp) :: taud_ab(4) 
        real(wp) :: tau_perp_ab(4) 
        real(wp) :: tau_par_ab(4) 
        real(wp) :: wt_ab(4) 
        real(wp) :: visc_eff_ab(4) 

        integer :: n_iter 

        nx    = size(ux_b,1)
        ny    = size(ux_b,2) 
        nz_aa = size(zeta_aa,1) 

        ! Consistency check 
        if (n_glen .ne. 3.0_prec) then 
            write(*,*) "calc_visc_eff_3D:: Error: currently, the L1L2 solver &
            & with dynamic viscosity can only be used with n_glen=3."
            stop 
        end if 

        ! Calculate scaling factors
        inv_4dx = 1.0_prec / (4.0_prec*dx) 
        inv_4dy = 1.0_prec / (4.0_prec*dy) 

        ! Calculate squared minimum strain rate 
        eps_0_sq = eps_0*eps_0 

        wt_ab = 1.0 
        wt = sum(wt_ab)
        wt_ab = wt_ab / wt 
              
        do j = 1, ny 
        do i = 1, nx 

            im1 = max(i-1,1) 
            ip1 = min(i+1,nx) 
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny) 

            ! Calculate effective strain components from horizontal stretching on ab-nodes
            call staggerdiff_nodes_acx_ab_ice(dudx_ab,ux_b,f_ice,i,j,dx)
            call staggerdiff_nodes_acy_ab_ice(dvdy_ab,uy_b,f_ice,i,j,dy)

            ! Calculate cross terms on ab-nodes
            call staggerdiffcross_nodes_acx_ab_ice(dudy_ab,ux_b,f_ice,i,j,dy)
            call staggerdiffcross_nodes_acy_ab_ice(dvdx_ab,uy_b,f_ice,i,j,dx)
            
            ! Calculate the 'parallel' effective strain rate from P12, Eq. 17
            eps_par_sq = dudx_ab**2 + dvdy_ab**2 + dudx_ab*dvdy_ab &
                        + 0.25_prec*(dudy_ab+dvdx_ab)**2 + eps_0_sq
            eps_par_ab = sqrt(eps_par_sq) 

            ! Get current magnitude of driving stress on ab-nodes 
            call stagger_nodes_acx_ab_ice(taudx_ab,taud_acx,f_ice,i,j)
            call stagger_nodes_acy_ab_ice(taudy_ab,taud_acy,f_ice,i,j)
            taud_ab = sqrt(taudx_ab**2 + taudy_ab**2)

            ! Now calculate viscosity at each layer 
            ! using the root-finding method of CISM
            ! Note this method is only valid for n_glen = 3!!!
            ! effstrain = A * (tau_parallel^2 + tau_perp^2)^{(n-1)/2} * tau_parallel
            ! y = A * (x^2 + tau^2)^{(n-1)/2} * x 

            do k = 1, nz_aa 

! ajr: Although logically I would choose to use the ab-node values 
! of ATT, calculate viscosity at each ab node and then average, this 
! seems to reduce stability of the model. Rather it seems to work 
! better by only calculating the effective strain rate at each ab-node,
! and center it, then multiply with the centered ATT value to get visc. 
! So that is why the central ATT value is used below. This should be 
! investigated further in the future perhaps.
if (.TRUE.) then  
                ! Get the rate factor on ab-nodes too
                call stagger_nodes_aa_ab_ice(ATT_ab,ATT(:,:,k),f_ice,i,j)
else
                ! Just use the aa-node central value of ATT 
                ATT_ab = ATT(i,j,k)

end if
                
                tau_perp_ab = taud_ab*(1.0_prec-zeta_aa(k))

if (.TRUE.) then 
    ! CISM root equation for n_glen=3 only 

                do q = 1, 4 
                    call calc_glen3_root(tau_par_ab(q),tau_perp_ab(q),eps_par_ab(q),ATT_ab(q))
                end do 

else 
    ! Root finding code (more expensive, but works for ISMIPHOM)
    ! Crashed for a random Antarctica simulation - needs testing! 

                ! To do: implement for solving all four ab nodes...

                ! a  = tau_perp_ab**2 
                ! b  = eps_par_ab / ATT_ab 
                ! np = (n_glen-1)/2.0_wp 

                ! !write(*,*) 'newton', a, b, np 
                ! call solve_secant(tau_par_ab,n_iter,10e3,1e-3,50,funY,a,b,np,.FALSE.) 
                ! !call solve_newton(tau_par_ab,n_iter,10e3,1e-3,50,funY,funYp,a,b,np,.FALSE.)
                ! !stop

end if 
                
                visc_eff_ab = 1.0_prec / (2.0_prec*ATT_ab*(tau_par_ab**2+tau_perp_ab**2)) 
                
                visc_eff(i,j,k) = sum(visc_eff_ab*wt_ab)

            end do 

        end do 
        end do  

        ! Extrapolate viscosity to bordering ice-free or partially ice-covered cells
        do j=1, ny
        do i=1, nx

            ! Get neighbor indices
            im1 = max(i-1,1) 
            ip1 = min(i+1,nx) 
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny) 
            
            if ( f_ice(i,j) .lt. 1.0 .and. &
                count([f_ice(im1,j),f_ice(ip1,j),f_ice(i,jm1),f_ice(i,jp1)] .eq. 1.0_wp) .gt. 0 ) then 
                ! Ice-free (or partially ice-free) with ice-covered neighbors

                visc_eff(i,j,:) = 0.0 
                wt = 0.0 

                if (f_ice(im1,j).eq.1.0) then 
                    visc_eff(i,j,:) = visc_eff(i,j,:) + visc_eff(im1,j,:) 
                    wt = wt + 1.0 
                end if 
                if (f_ice(ip1,j).eq.1.0) then 
                    visc_eff(i,j,:) = visc_eff(i,j,:) + visc_eff(ip1,j,:) 
                    wt = wt + 1.0 
                end if 
                if (f_ice(i,jm1).eq.1.0) then 
                    visc_eff(i,j,:) = visc_eff(i,j,:) + visc_eff(i,jm1,:) 
                    wt = wt + 1.0 
                end if 
                if (f_ice(i,jp1).eq.1.0) then 
                    visc_eff(i,j,:) = visc_eff(i,j,:) + visc_eff(i,jp1,:) 
                    wt = wt + 1.0 
                end if 
                
                if (wt .gt. 0.0) then 
                    visc_eff(i,j,:) = visc_eff(i,j,:) / wt 

                end if 

            end if 

        end do 
        end do 

        ! Apply boundary conditions as needed 
        if (trim(boundaries) .eq. "periodic") then

            visc_eff(1,:,:)    = visc_eff(nx-1,:,:) 
            visc_eff(nx-1,:,:) = visc_eff(2,:,:) 
            visc_eff(:,1,:)    = visc_eff(:,ny-1,:)
            visc_eff(:,ny,:)   = visc_eff(:,2,:) 

        else if (trim(boundaries) .eq. "periodic-x") then 
            
            visc_eff(1,:,:)    = visc_eff(nx-1,:,:) 
            visc_eff(nx-1,:,:) = visc_eff(2,:,:) 
            visc_eff(:,1,:)    = visc_eff(:,2,:)
            visc_eff(:,ny,:)   = visc_eff(:,ny-1,:) 

        else if (trim(boundaries) .eq. "infinite") then 
            
            visc_eff(1,:,:)    = visc_eff(2,:,:) 
            visc_eff(nx,:,:)   = visc_eff(nx-1,:,:) 
            visc_eff(:,1,:)    = visc_eff(:,2,:)
            visc_eff(:,ny,:)   = visc_eff(:,ny-1,:) 

        end if 

        ! Treat the corners to avoid extremes
        visc_eff(1,1,:)   = 0.5*(visc_eff(2,1,:)+visc_eff(1,2,:))
        visc_eff(1,ny,:)  = 0.5*(visc_eff(2,ny,:)+visc_eff(1,ny-1,:))
        visc_eff(nx,1,:)  = 0.5*(visc_eff(nx,2,:)+visc_eff(nx-1,1,:))
        visc_eff(nx,ny,:) = 0.5*(visc_eff(nx-1,ny,:)+visc_eff(nx,ny-1,:))

        return 

    end subroutine calc_visc_eff_3D

    subroutine calc_glen3_root(tau_par,tau_perp,eps_par,ATT)
        ! Calculate tau_par using analytical root for quadratic equation
        ! (only valid when n_glen=3!)

        implicit none 

        real(wp), intent(OUT) :: tau_par 
        real(wp), intent(IN)  :: tau_perp 
        real(wp), intent(IN)  :: eps_par 
        real(wp), intent(IN)  :: ATT 
        
        ! Local variables 
        real(wp) :: a, b, c 
        real(wp) :: rootA, rootB 

        a = tau_perp**2 
        b = -eps_par / ATT 
        c = sqrt(b**2/4.0_prec + a**3/27.0_prec) 

        rootA = (-b/2.0_prec + c)**(1.0_prec/3.0_prec)

        if (a**3/(27.0_prec) > 1.d-6 * (b**2/4.0_prec)) then
            rootB = -(b/2.0_prec + c)**(1.0_prec/3.0_prec)
        else    ! b/2 + c is small; compute solution to first order without subtracting two large, nearly equal numbers
            rootB = -a / (3.0_prec*(abs(b))**(1.0_prec/3.0_prec))
        end if

        tau_par = rootA + rootB

        return 

    end subroutine calc_glen3_root

    subroutine calc_visc_eff_int(visc_eff_int,visc_eff,H_ice,f_ice,zeta_aa,boundaries)

        implicit none 

        real(wp), intent(OUT) :: visc_eff_int(:,:)
        real(wp), intent(IN)  :: visc_eff(:,:,:)
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: zeta_aa(:)
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer :: i, j, nx, ny
        integer :: im1, ip1, jm1, jp1  
        real(wp) :: H_now
        real(wp) :: visc_eff_mean 
        real(wp) :: wt 

        nx = size(visc_eff_int,1)
        ny = size(visc_eff_int,2)


        do j = 1, ny 
        do i = 1, nx

            ! Calculate the vertically averaged viscosity for this point
            visc_eff_mean = integrate_trapezoid1D_pt(visc_eff(i,j,:),zeta_aa) 
            
            if (f_ice(i,j) .eq. 1.0) then 
                visc_eff_int(i,j) = visc_eff_mean*H_ice(i,j) 
            else
                !visc_eff_int(i,j) = visc_eff_mean 
                visc_eff_int(i,j) = 0.0_wp
            end if 

        end do 
        end do 

        ! Now extrapolate to ice-free or partially ice-free neighbors
        do j = 1, ny 
        do i = 1, nx

            im1 = max(i-1,1) 
            ip1 = min(i+1,nx) 
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny) 

            if ( f_ice(i,j) .lt. 1.0 .and. &
                count([f_ice(im1,j),f_ice(ip1,j),f_ice(i,jm1),f_ice(i,jp1)] .eq. 1.0) .gt. 0) then 
                ! Ice-free or partially ice-free point at ice margin

                visc_eff_int(i,j) = 0.0_wp 
                wt = 0.0_wp 

                if (f_ice(im1,j).eq.1.0) then 
                    visc_eff_int(i,j) = visc_eff_int(i,j) + visc_eff_int(im1,j) 
                    wt = wt + 1.0 
                end if 
                if (f_ice(ip1,j).eq.1.0) then 
                    visc_eff_int(i,j) = visc_eff_int(i,j) + visc_eff_int(ip1,j) 
                    wt = wt + 1.0 
                end if 
                if (f_ice(i,jm1).eq.1.0) then 
                    visc_eff_int(i,j) = visc_eff_int(i,j) + visc_eff_int(i,jm1) 
                    wt = wt + 1.0 
                end if 
                if (f_ice(i,jp1).eq.1.0) then 
                    visc_eff_int(i,j) = visc_eff_int(i,j) + visc_eff_int(i,jp1) 
                    wt = wt + 1.0 
                end if 
                
                if (wt .gt. 0.0) then 
                    visc_eff_int(i,j) = visc_eff_int(i,j) / wt 
                end if 

            end if 

        end do 
        end do 

        ! Apply boundary conditions as needed 
        if (trim(boundaries) .eq. "periodic") then

            visc_eff_int(1,:)    = visc_eff_int(nx-1,:) 
            visc_eff_int(nx-1,:) = visc_eff_int(2,:) 
            visc_eff_int(:,1)    = visc_eff_int(:,ny-1)
            visc_eff_int(:,ny)   = visc_eff_int(:,2) 

        else if (trim(boundaries) .eq. "periodic-x") then 
            
            visc_eff_int(1,:)    = visc_eff_int(nx-1,:) 
            visc_eff_int(nx-1,:) = visc_eff_int(2,:) 
            visc_eff_int(:,1)    = visc_eff_int(:,2)
            visc_eff_int(:,ny)   = visc_eff_int(:,ny-1) 

        else if (trim(boundaries) .eq. "infinite") then 
            
            visc_eff_int(1,:)    = visc_eff_int(2,:) 
            visc_eff_int(nx,:)   = visc_eff_int(nx-1,:) 
            visc_eff_int(:,1)    = visc_eff_int(:,2)
            visc_eff_int(:,ny)   = visc_eff_int(:,ny-1) 

        end if 

        ! Treat the corners to avoid extremes
        visc_eff_int(1,1)   = 0.5*(visc_eff_int(2,1)+visc_eff_int(1,2))
        visc_eff_int(1,ny)  = 0.5*(visc_eff_int(2,ny)+visc_eff_int(1,ny-1))
        visc_eff_int(nx,1)  = 0.5*(visc_eff_int(nx,2)+visc_eff_int(nx-1,1))
        visc_eff_int(nx,ny) = 0.5*(visc_eff_int(nx-1,ny)+visc_eff_int(nx,ny-1))

        return

    end subroutine calc_visc_eff_int

    subroutine calc_basal_stress(taub_acx,taub_acy,beta_acx,beta_acy,ux_b,uy_b)
        ! Calculate the basal stress resulting from sliding (friction times velocity)
        ! Note: calculated on ac-nodes.
        ! taub [Pa] 
        ! beta [Pa a m-1]
        ! u    [m a-1]
        ! taub = beta*u (here defined with taub in the same direction as u)

        implicit none 

        real(prec), intent(OUT) :: taub_acx(:,:)        ! [Pa] Basal stress (acx nodes)
        real(prec), intent(OUT) :: taub_acy(:,:)        ! [Pa] Basal stress (acy nodes)
        real(prec), intent(IN)  :: beta_acx(:,:)        ! [Pa a m-1] Effective basal friction (acx nodes)
        real(prec), intent(IN)  :: beta_acy(:,:)        ! [Pa a m-1] Effective basal friction (acy nodes)
        real(prec), intent(IN)  :: ux_b(:,:)            ! [m a-1] depth-ave velocity (acx nodes)
        real(prec), intent(IN)  :: uy_b(:,:)            ! [m a-1] depth-ave velocity (acy nodes)
        
        ! Local variables 
        integer :: i, j, nx, ny 

        nx = size(taub_acx,1)
        ny = size(taub_acy,2) 

        do j = 1, ny 
        do i = 1, nx 

            taub_acx(i,j) = beta_acx(i,j) * ux_b(i,j) 
            taub_acy(i,j) = beta_acy(i,j) * uy_b(i,j) 

        end do 
        end do  

        return 

    end subroutine calc_basal_stress

    elemental subroutine limit_vel(u,u_lim)
        ! Apply a velocity limit (for stability)

        implicit none 

        real(prec), intent(INOUT) :: u 
        real(prec), intent(IN)    :: u_lim

        real(prec), parameter :: tol = 1e-10
        
        u = min(u, u_lim)
        u = max(u,-u_lim)

        ! Also avoid underflow errors 
        if (abs(u) .lt. tol) u = 0.0 

        return

    end subroutine limit_vel
    

    subroutine solve_newton(x,n_iter,x_init,tol,n_max,f,fp,a,b,np,debug)
        ! Estimate the zero of f(x) using Newton's method. 
        ! Adapted from: 
        ! https://faculty.washington.edu/rjl/classes/am583s2013/notes/fortran_newton.html

        implicit none

        real(wp), intent(OUT) :: x          ! Best guess of root
        integer,  intent(OUT) :: n_iter     ! Number of iterations to reach it
        real(wp), intent(IN)  :: x_init     ! Initial guess
        real(wp), intent(IN)  :: tol        ! Tolerance to convergence
        integer,  intent(IN)  :: n_max      ! Maximum iterations allowed
        real(wp), external    :: f          ! Function to find root of 
        real(wp), external    :: fp         ! Derivative of Function
        real(wp), intent(IN)  :: a,b,np     ! Additional function parameters
        logical,  intent(IN)  :: debug      ! Print iteration information?
        
        ! Declare any local variables:
        real(wp) :: deltax, fx, fxprime
        integer    :: k

        ! Save initial guess
        x = x_init

        ! Newton iteration to find a zero of f(x) 
        n_iter = 0 

        do k = 1, n_max

            n_iter = n_iter + 1 

            ! evaluate function and its derivative:
            fx      = f(x,a,b,np)
            fxprime = fp(x,a,b,np)

            if (debug) then
                !write(*,*) n_iter, "x, f(x) = ", x, fx
                write(*,*) n_iter, "a, b, np, x, f(x) = ", a, b, np, x, fx
            end if 

            if (abs(fx) < tol) then
                exit  ! jump out of do loop
            end if

            ! Compute Newton increment x:
            deltax = fx/fxprime

            ! update x:
            x = x - deltax

        end do


        if (n_iter .eq. n_max .and. abs(fx) > tol) then
            write(*,*) "solve_newton:: Warning: no convergence."
        end if

        return 

    end subroutine solve_newton

    subroutine solve_secant(x,n_iter,x_init,tol,n_max,f,a,b,np,debug) 
        ! Estimate the zero of f(x) using the Secant method. 
        ! Adapted from: 
        ! http://jean-pierre.moreau.pagesperso-orange.fr/Fortran/secant_f90.txt
        ! https://rosettacode.org/wiki/Roots_of_a_function#Fortran

        implicit none 

        real(wp), intent(OUT) :: x          ! Best guess of root
        integer,  intent(OUT) :: n_iter     ! Number of iterations to reach it
        real(wp), intent(IN)  :: x_init     ! Initial guess
        real(wp), intent(IN)  :: tol        ! Tolerance to convergence
        integer,  intent(IN)  :: n_max      ! Maximum iterations allowed
        real(wp), external    :: f          ! Function to find root of 
        real(wp), intent(IN)  :: a,b,np     ! Additional function parameters
        logical,  intent(IN)  :: debug      ! Print iteration information?
        
        ! Local variables 
        integer  :: n 
        real(wp) :: x1, x2  
        real(wp) :: y1, y2  
        real(wp) :: d 

        ! Set x to initial guess and a slightly different value
        x1 = x_init 
        x2 = x_init*0.5_wp
        if (x_init .eq. 0.0_wp) x2 = x_init + 0.1_wp  

        n_iter = 0

        ! Start iterations
        do n = 1, n_max 

            n_iter = n_iter + 1 

            ! Calculate new value of y at x
            y1 = f(x1,a,b,np)
            y2 = f(x2,a,b,np)

            if (debug) write(*,*) n_iter, x1, x2, y1, y2 

            if (abs(y2) < tol) exit

            
            d = (x2 - x1) / (y2 - y1) * y2
            
            x1 = x2
            x2 = x2 - d

        end do 

        x = x2 

        if (n_iter .eq. n_max .and. abs(y2) > tol) then 
            write(*,*) "solve_secant:: Warning: no convergence."
        end if 

        return 

    end subroutine solve_secant

    ! Functions for root finding methods, specific to L1L2 solver: 

    real(wp) function funY(x,a,b,n)
        real(wp) :: x,a,b,n 
        funY = (x * (x**2 + a)**n)/b - 1
    end function funY
    
    real(wp) function funYp(x,a,b,n)
        real(wp) :: x,a,b,n
        funYp = (2.0_wp * n * x**2 * (x**2 + a)**(n-1.0_wp)/b) + ((x**2 + a)**n/b)
    end function funYp


end module velocity_l1l2
