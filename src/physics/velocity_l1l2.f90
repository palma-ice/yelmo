module velocity_l1l2

    use yelmo_defs ,only  : prec, wp, rho_ice, rho_sw, rho_w, g
    use yelmo_tools, only : stagger_aa_ab, stagger_aa_ab_ice, stagger_ab_aa_ice, & 
                    stagger_nodes_aa_ab_ice, stagger_nodes_acx_ab_ice, stagger_nodes_acy_ab_ice, &
                    staggerdiffx_nodes_aa_ab_ice, staggerdiffy_nodes_aa_ab_ice, &
                    staggerdiff_nodes_acx_ab_ice, staggerdiff_nodes_acy_ab_ice, &
                    staggerdiffcross_nodes_acx_ab_ice, staggerdiffcross_nodes_acy_ab_ice, &
                    staggerdiffcross_aa_acx_ice, staggerdiffcross_aa_acy_ice, &
                    integrate_trapezoid1D_1D, integrate_trapezoid1D_pt, minmax, &
                    calc_vertical_integrated_2D, calc_vertical_integrated_3D

    use basal_dragging 
    use solver_ssa_sico5 

    implicit none 

    type l1l2_param_class

        character(len=256) :: ssa_lis_opt 
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
                                  c_bed,taud_acx,taud_acy,H_ice,H_grnd,f_grnd, &
                                  f_grnd_acx,f_grnd_acy,ATT,zeta_aa,z_sl,z_bed,dx,dy,n_glen,par)
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
        real(prec), intent(IN)    :: H_ice(:,:)         ! [m]
        real(prec), intent(IN)    :: H_grnd(:,:)        ! [m]
        real(prec), intent(IN)    :: f_grnd(:,:)        ! [-]
        real(prec), intent(IN)    :: f_grnd_acx(:,:)    ! [-]
        real(prec), intent(IN)    :: f_grnd_acy(:,:)    ! [-]
        real(prec), intent(IN)    :: ATT(:,:,:)         ! [a^-1 Pa^-n_glen]
        real(prec), intent(IN)    :: zeta_aa(:)         ! [-]
        real(prec), intent(IN)    :: z_sl(:,:)          ! [m]
        real(prec), intent(IN)    :: z_bed(:,:)         ! [m]
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
        real(prec), allocatable :: f_ice(:,:)

        real(prec) :: L2_norm 

        integer :: ij(2) 

        nx    = size(ux,1)
        ny    = size(ux,2)
        nz_aa = size(ux,3)

        ! Prepare local variables 
        allocate(ux_b_nm1(nx,ny))
        allocate(uy_b_nm1(nx,ny))

        allocate(ssa_mask_acx_ref(nx,ny))
        allocate(ssa_mask_acy_ref(nx,ny))

        allocate(visc_eff_ab(nx,ny,nz_aa)) 
        allocate(f_ice(nx,ny))

        f_ice = 0.0 
        where(H_ice .gt. 0.0) f_ice = 1.0 

        ! Store original ssa mask before iterations
        ssa_mask_acx_ref = ssa_mask_acx
        ssa_mask_acy_ref = ssa_mask_acy
            
        ! Initially set error very high 
        ssa_err_acx = 1.0_prec 
        ssa_err_acy = 1.0_prec 
        
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
                    
                    call calc_visc_eff_3D(visc_eff,visc_eff_ab,ux_b,uy_b,taud_acx,taud_acy,ATT, &
                                                H_ice,zeta_aa,dx,dy,n_glen,par%eps_0,par%boundaries)

                case DEFAULT 

                    write(*,*) "calc_velocity_l1l2:: Error: visc_method not recognized."
                    write(*,*) "visc_method = ", par%visc_method 
                    stop 
                    
            end select

            ! Calculate depth-integrated effective viscosity
            ! Note L19 uses eta_bar*H in the ssa equation. Yelmo uses eta_int=eta_bar*H directly.
            visc_eff_int = calc_vertical_integrated_2D(visc_eff,zeta_aa) 
            where(H_ice .gt. 0.0_prec) visc_eff_int = visc_eff_int*H_ice 

            ! Smooth the viscosity at the ice margins if it is too low
            ! to avoid singularities (mainly for EISMINT/dome experiments)
            !call smooth_visc_eff_int_margin(visc_eff_int,H_ice)

            ! Calculate beta (at the ice base)
            call calc_beta(beta,c_bed,ux_b,uy_b,H_ice,H_grnd,f_grnd,z_bed,z_sl,par%beta_method, &
                                par%beta_const,par%beta_q,par%beta_u0,par%beta_gl_scale,par%beta_gl_f, &
                                par%H_grnd_lim,par%beta_min,par%boundaries)

            ! Stagger beta and beta_eff 
            call stagger_beta(beta_acx,beta_acy,beta,f_grnd,f_grnd_acx,f_grnd_acy,par%beta_gl_stag,par%boundaries)

            ! =========================================================================================
            ! Step 2: determine the basal velocity ux_b/uy_b 

            if (par%no_slip) then 
                ! Simply set ux_b/uy_b equal to zero, as no sliding is allowed 

                ux_b = 0.0_prec 
                uy_b = 0.0_prec 

            else 
                ! Call the SSA solver to obtain new estimate of ux_b/uy_b

if (.TRUE.) then 
                if (iter .gt. 1) then
                    ! Update ssa mask based on convergence with previous step to reduce area being solved 
                    call update_ssa_mask_convergence(ssa_mask_acx,ssa_mask_acy,ssa_err_acx,ssa_err_acy,err_lim=real(1e-5,prec))
                    !call update_ssa_mask_convergence(ssa_mask_acx,ssa_mask_acy,ssa_err_acx,ssa_err_acy,err_lim=par%ssa_iter_conv*1e-2)  
                end if 
end if 
            
                ! Call ssa solver
                call calc_vxy_ssa_matrix(ux_b,uy_b,L2_norm,beta_acx,beta_acy,visc_eff_int,  &
                                         ssa_mask_acx,ssa_mask_acy,H_ice,taud_acx,taud_acy,H_grnd,z_sl, &
                                         z_bed,dx,dy,par%ssa_vel_max,par%boundaries,par%ssa_lis_opt)


                ! Apply relaxation to keep things stable
                call relax_ssa(ux_b,uy_b,ux_b_nm1,uy_b_nm1,rel=par%ssa_iter_rel)
            
            end if 

            ! Check for convergence
            is_converged = check_vel_convergence_l2rel(ux_b,uy_b,ux_b_nm1,uy_b_nm1,ssa_mask_acx.gt.0,     &
                                                       ssa_mask_acy.gt.0,par%ssa_iter_conv,iter,par%ssa_iter_max, &
                                                       par%ssa_write_log,use_L2_norm=.FALSE.,L2_norm=L2_norm)

            ! Calculate an L1 error metric over matrix for diagnostics
            call check_vel_convergence_l1rel_matrix(ssa_err_acx,ssa_err_acy,ux_b,uy_b,ux_b_nm1,uy_b_nm1)

            ! Store current total iterations for output
            ssa_iter_now = iter 

            ! =========================================================================================
            ! Update additional fields based on output of solver
            
            ! Calculate basal stress 
            call calc_basal_stress(taub_acx,taub_acy,beta_acx,beta_acy,ux_b,uy_b)

            ! Exit iterations if ssa solution has converged
            if (is_converged) exit 
            
        end do 

        ! Iterations are finished, finalize calculations of 3D velocity field 

        ! Calculate the 3D horizontal velocity field
        call calc_vel_horizontal_3D(ux,uy,ux_b,uy_b,taud_acx,taud_acy,visc_eff,ATT,H_ice, &
                                    f_ice,zeta_aa,dx,dy,n_glen,par%eps_0,par%boundaries)
        ! call calc_vel_horizontal_3D_0(ux,uy,ux_b,uy_b,taud_acx,taud_acy,visc_eff_ab,ATT,H_ice, &
        !                                     zeta_aa,dx,dy,n_glen,par%eps_0,par%boundaries)

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
        nz_ac = nz_aa+1

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

    subroutine calc_vel_horizontal_3D_0(ux,uy,ux_b,uy_b,taud_acx,taud_acy,visc_eff_ab,ATT,H_ice,zeta_aa,dx,dy,n_glen,eps_0,boundaries)
        ! Caluculate the 3D horizontal velocity field (ux,uy)
        ! for the L1L2 solver following Perego et al. (2012)
        ! and the blueprint by Lipscomb et al. (2019) in CISM

        implicit none 

        real(prec), intent(OUT) :: ux(:,:,:) 
        real(prec), intent(OUT) :: uy(:,:,:) 
        real(prec), intent(IN)  :: ux_b(:,:) 
        real(prec), intent(IN)  :: uy_b(:,:) 
        real(prec), intent(IN)  :: taud_acx(:,:) 
        real(prec), intent(IN)  :: taud_acy(:,:)
        real(prec), intent(IN)  :: visc_eff_ab(:,:,:)   ! on ab-nodes already 
        real(prec), intent(IN)  :: ATT(:,:,:)  
        real(prec), intent(IN)  :: H_ice(:,:)
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
        real(prec) :: zeta_ac1, zeta_ac0 
        real(prec) :: H_ice_ac 
        real(prec) :: dw1dx, dw2dx, dw3dx 
        real(prec) :: dw1dy, dw2dy, dw3dy 
        real(prec) :: tau_xz_ab, tau_yz_ab 
        real(prec) :: tau_eff_sq_ab, ATT_ab, depth_ab 
        real(prec) :: fact_ac 
        real(prec), allocatable :: dudx_ab(:,:)
        real(prec), allocatable :: dvdy_ab(:,:)
        real(prec), allocatable :: dudy_ab(:,:)
        real(prec), allocatable :: dvdx_ab(:,:)
        real(prec), allocatable :: H_ice_ab(:,:) 
        real(prec), allocatable :: visc_eff_int3D_ab(:,:,:) 
        real(prec), allocatable :: tau_par_ab(:,:,:) 
        real(prec), allocatable :: work1_ab(:,:)
        real(prec), allocatable :: work2_ab(:,:)
        real(prec), allocatable :: work3_ab(:,:)
        real(prec), allocatable :: tau_xz(:,:,:) 
        real(prec), allocatable :: tau_yz(:,:,:) 
        real(prec), allocatable :: fact_ab(:,:) 

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

        real(wp), allocatable :: f_ice(:,:) 

        real(prec) :: eps_par_sq, eps_par 
        real(prec) :: p1, eps_0_sq 
        real(prec) :: dzeta 

        nx    = size(ux,1)
        ny    = size(ux,2) 
        nz_aa = size(ux,3) 

        ! Allocate local arrays 
        allocate(dudx_ab(nx,ny)) 
        allocate(dvdy_ab(nx,ny)) 
        allocate(dudy_ab(nx,ny)) 
        allocate(dvdx_ab(nx,ny)) 
        allocate(H_ice_ab(nx,ny))
        allocate(visc_eff_int3D_ab(nx,ny,nz_aa)) 
        allocate(tau_par_ab(nx,ny,nz_aa))
        allocate(work1_ab(nx,ny)) 
        allocate(work2_ab(nx,ny)) 
        allocate(work3_ab(nx,ny)) 
        allocate(tau_xz(nx,ny,nz_aa))
        allocate(tau_yz(nx,ny,nz_aa))
        allocate(fact_ab(nx,ny))

        allocate(f_ice(nx,ny)) 

        ! Define f_ice locally 
        f_ice = 0.0 
        where(H_ice .gt. 0.0) f_ice = 1.0 

        wt_ab = 1.0 
        wt = sum(wt_ab)
        wt_ab = wt_ab / wt 


        ! Calculate scaling factors
        inv_4dx = 1.0_prec / (4.0_prec*dx) 
        inv_4dy = 1.0_prec / (4.0_prec*dy) 

        ! Calculate exponent
        p1 = (n_glen - 1.0_prec) / 2.0_prec 

        ! Calculate squared minimum strain rate 
        eps_0_sq = eps_0*eps_0 

        ! Initialize integrated viscosity field
        visc_eff_int3D_ab = 0.0_prec 

        ! Step 1: compute basal strain rates on ab-nodes and viscosity       
        do j = 1, ny 
        do i = 1, nx 

            im1 = max(i-1,1) 
            ip1 = min(i+1,nx) 
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny) 

            ! Calculate effective strain components from horizontal stretching on ab-nodes
            dudx_ab(i,j) = ( (ux_b(ip1,j) - ux_b(im1,j)) + (ux_b(ip1,jp1) - ux_b(im1,jp1)) ) *inv_4dx
            dvdy_ab(i,j) = ( (uy_b(i,jp1) - uy_b(i,jm1)) + (uy_b(ip1,jp1) - uy_b(ip1,jm1)) ) *inv_4dy 

            ! Calculate of cross terms on ab-nodes
            dudy_ab(i,j) = (ux_b(i,jp1) - ux_b(i,j)) / dx 
            dvdx_ab(i,j) = (uy_b(ip1,j) - uy_b(i,j)) / dy 

            ! Calculate the 'parallel' effective strain rate from P12, Eq. 17
            eps_par_sq = dudx_ab(i,j)**2 + dvdy_ab(i,j)**2 + dudx_ab(i,j)*dvdy_ab(i,j) &
                        + 0.25_prec*(dudy_ab(i,j)+dvdx_ab(i,j))**2 + eps_0_sq
            eps_par    = sqrt(eps_par_sq) 

            ! Compute the 'parallel' shear stress for each layer (tau_parallel)
            do k = 1, nz_aa 
                tau_par_ab(i,j,k) = 2.d0 * visc_eff_ab(i,j,k) * eps_par
            end do 

            ! Compute the integral of visc_eff from the base of each layer to the surface (P12, Eq. 28)

            H_ice_ab(i,j) = 0.25_prec*(H_ice(i,j)+H_ice(ip1,j)+H_ice(i,jp1)+H_ice(ip1,jp1))
            
            ! Start at the surface
            visc_eff_int3D_ab(i,j,nz_aa) = visc_eff_ab(i,j,nz_aa) &
                                            * (zeta_aa(nz_aa)-zeta_aa(nz_aa-1))*H_ice_ab(i,j)

            ! Integrate down to near the base 
            do k = nz_aa-1, 2, -1 
                zeta_ac1 = 0.5_prec*(zeta_aa(k+1)+zeta_aa(k))
                zeta_ac0 = 0.5_prec*(zeta_aa(k)+zeta_aa(k-1))
                visc_eff_int3D_ab(i,j,k) = visc_eff_int3D_ab(i,j,k+1) &
                                        + visc_eff_ab(i,j,k) * (zeta_ac1-zeta_ac0)*H_ice_ab(i,j)
            end do 
            
            ! Get basal value
            visc_eff_int3D_ab(i,j,1) = visc_eff_int3D_ab(i,j,2) &
                                        + visc_eff_ab(i,j,1) * (zeta_aa(2)-zeta_aa(1))*H_ice_ab(i,j)
            
        end do  
        end do 

        ! Loop over layers 
        do k = 1, nz_aa

            ! Calculate working arrays for this layer 
            work1_ab = visc_eff_int3D_ab(:,:,k) * (2.d0*dudx_ab + dvdy_ab) 
            work2_ab = visc_eff_int3D_ab(:,:,k) *      (dudy_ab + dvdx_ab)
            work3_ab = visc_eff_int3D_ab(:,:,k) * (2.d0*dvdy_ab + dudx_ab) 

            ! Loop over horizontal grid points 
            do j = 1, ny 
            do i = 1, nx 

                im1 = max(i-1,1) 
                ip1 = min(i+1,nx) 
                jm1 = max(j-1,1) 
                jp1 = min(j+1,ny) 

                ! Calculate derivatives of work arrays on ac-nodes 
                dw1dx = ( 0.25*(work1_ab(i,j)+work1_ab(ip1,j)+work1_ab(ip1,jm1)+work1_ab(i,jm1)) &
                        - 0.25*(work1_ab(i,j)+work1_ab(im1,j)+work1_ab(im1,jm1)+work1_ab(i,jm1)) ) / dx 
                dw1dy = ( 0.25*(work1_ab(i,j)+work1_ab(i,jp1)+work1_ab(im1,jp1)+work1_ab(im1,j)) &
                        - 0.25*(work1_ab(i,j)+work1_ab(i,jm1)+work1_ab(im1,jm1)+work1_ab(im1,j)) ) / dy 
                
                dw2dx = ( 0.25*(work2_ab(i,j)+work2_ab(ip1,j)+work2_ab(ip1,jm1)+work2_ab(i,jm1)) &
                        - 0.25*(work2_ab(i,j)+work2_ab(im1,j)+work2_ab(im1,jm1)+work2_ab(i,jm1)) ) / dx 
                dw2dy = ( 0.25*(work2_ab(i,j)+work2_ab(i,jp1)+work2_ab(im1,jp1)+work2_ab(im1,j)) &
                        - 0.25*(work2_ab(i,j)+work2_ab(i,jm1)+work2_ab(im1,jm1)+work2_ab(im1,j)) ) / dy 
                
                dw3dx = ( 0.25*(work3_ab(i,j)+work3_ab(ip1,j)+work3_ab(ip1,jm1)+work3_ab(i,jm1)) &
                        - 0.25*(work3_ab(i,j)+work3_ab(im1,j)+work3_ab(im1,jm1)+work3_ab(i,jm1)) ) / dx 
                dw3dy = ( 0.25*(work3_ab(i,j)+work3_ab(i,jp1)+work3_ab(im1,jp1)+work3_ab(im1,j)) &
                        - 0.25*(work3_ab(i,j)+work3_ab(i,jm1)+work3_ab(im1,jm1)+work3_ab(im1,j)) ) / dy 
                
                ! Calculate shear stress on ac-nodes
                tau_xz(i,j,k) = -(1.0_prec-zeta_aa(k))*taud_acx(i,j) + 2.0_prec*dw1dx + dw2dy
                tau_yz(i,j,k) = -(1.0_prec-zeta_aa(k))*taud_acy(i,j) + dw2dx + 2.0_prec*dw3dy

            end do 
            end do  

        end do 

        ! Note: if basal sliding is zero, then the shear stress calculated 
        ! above (tau_xz,tau_yz) is identical to SIA, and tau_par is also zero.
        ! In that case, the below calculations are identical to the normal 
        ! SIA solver.

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
                
                tau_par_ab4_up = [tau_par_ab(i,j,k),tau_par_ab(im1,j,k), &
                                        tau_par_ab(im1,jm1,k),tau_par_ab(ip1,jm1,k)]
                tau_par_ab4_dn = [tau_par_ab(i,j,k-1),tau_par_ab(im1,j,k-1), &
                                        tau_par_ab(im1,jm1,k-1),tau_par_ab(ip1,jm1,k-1)]
                tau_par_ab4 = 0.5_wp*(tau_par_ab4_up+tau_par_ab4_dn) 

                tau_eff_sq_ab4 = tau_xz_ab4**2 + tau_yz_ab4**2 + tau_par_ab4**2

                ! Get rate factor on ab-nodes
                call stagger_nodes_aa_ab_ice(ATT_ab4_up,ATT(:,:,k),f_ice,i,j)
                call stagger_nodes_aa_ab_ice(ATT_ab4_dn,ATT(:,:,k-1),f_ice,i,j)
                ATT_ab4 = 0.5_wp*(ATT_ab4_up+ATT_ab4_dn)

                ! Get ice thickness on ab-nodes
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

        return 

    end subroutine calc_vel_horizontal_3D_0

    subroutine calc_visc_eff_3D(visc_eff,visc_eff_ab,ux_b,uy_b,taud_acx,taud_acy,ATT,H_ice,zeta_aa,dx,dy,n_glen,eps_0,boundaries)
        ! Caluculate the 3D effective viscosity field
        ! for the L1L2 solver following Perego et al. (2012)
        ! and the blueprint by Lipscomb et al. (2019) in CISM

        implicit none 

        real(prec), intent(OUT) :: visc_eff(:,:,:)       
        real(prec), intent(OUT) :: visc_eff_ab(:,:,:) 
        real(prec), intent(IN)  :: ux_b(:,:) 
        real(prec), intent(IN)  :: uy_b(:,:) 
        real(prec), intent(IN)  :: taud_acx(:,:) 
        real(prec), intent(IN)  :: taud_acy(:,:)
        real(prec), intent(IN)  :: ATT(:,:,:)  
        real(prec), intent(IN)  :: H_ice(:,:)
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
        real(prec) :: tau_eff_sq_ab, ATT_ab 
        real(prec) :: dudx_ab, dvdy_ab, dudy_ab, dvdx_ab
        
        real(prec) :: eps_par_sq, eps_par_ab 
        real(prec) :: eps_0_sq 
        real(prec) :: taud_ab, tau_par_ab, tau_perp_ab 
        real(prec) :: a, b, c, rootA, rootB 
        real(prec) :: wt 

        nx    = size(ux_b,1)
        ny    = size(ux_b,2) 
        nz_aa = size(zeta_aa,1) 

        ! Consistency check 
        if (n_glen .ne. 3.0_prec) then 
            write(*,*) "calc_visc_eff_3D:: Error: currently, the L1L2 solver &
            & can only be used with n_glen=3."
            stop 
        end if 

        ! Calculate scaling factors
        inv_4dx = 1.0_prec / (4.0_prec*dx) 
        inv_4dy = 1.0_prec / (4.0_prec*dy) 

        ! Calculate squared minimum strain rate 
        eps_0_sq = eps_0*eps_0 

        ! Step 1: compute basal strain rates on ab-nodes and viscosity       
        do j = 1, ny 
        do i = 1, nx 

            im1 = max(i-1,1) 
            ip1 = min(i+1,nx) 
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny) 

            ! Calculate effective strain components from horizontal stretching on ab-nodes
            dudx_ab = ( (ux_b(ip1,j) - ux_b(im1,j)) + (ux_b(ip1,jp1) - ux_b(im1,jp1)) ) *inv_4dx
            dvdy_ab = ( (uy_b(i,jp1) - uy_b(i,jm1)) + (uy_b(ip1,jp1) - uy_b(ip1,jm1)) ) *inv_4dy 

            ! Calculate of cross terms on ab-nodes
            dudy_ab = (ux_b(i,jp1) - ux_b(i,j)) / dx 
            dvdx_ab = (uy_b(ip1,j) - uy_b(i,j)) / dy 

            ! Calculate the 'parallel' effective strain rate from P12, Eq. 17
            eps_par_sq = dudx_ab**2 + dvdy_ab**2 + dudx_ab*dvdy_ab &
                        + 0.25_prec*(dudy_ab+dvdx_ab)**2 + eps_0_sq
            eps_par_ab = sqrt(eps_par_sq) 


            ! Get current magnitude of driving stress on ab-nodes 
            taud_ab = sqrt( (0.5_prec*(taud_acx(i,j)+taud_acx(i,jp1)))**2 &
                          + (0.5_prec*(taud_acy(i,j)+taud_acy(ip1,j)))**2 )

            ! Now calculate viscosity at each layer 
            ! using the root-finding method of CISM
            ! Note this method is only valid for n_glen = 3!!!
            do k = 1, nz_aa 
                
                ATT_ab = 0.25_prec*(ATT(i,j,k)+ATT(ip1,j,k)+ATT(i,jp1,k)+ATT(ip1,jp1,k))
                
                tau_perp_ab = taud_ab*(1.0_prec-zeta_aa(k))

                a = tau_perp_ab**2 
                b = -eps_par_ab / ATT_ab 
                c = sqrt(b**2/4.0_prec + a**3/27.0_prec) 

                rootA = (-b/2.0_prec + c)**(1.0_prec/3.0_prec)

                if (a**3/(27.0_prec) > 1.d-6 * (b**2/4.0_prec)) then
                    rootB = -(b/2.0_prec + c)**(1.0_prec/3.0_prec)
                else    ! b/2 + c is small; compute solution to first order without subtracting two large, nearly equal numbers
                    rootB = -a / (3.0_prec*(abs(b))**(1.0_prec/3.0_prec))
                end if

                tau_par_ab = rootA + rootB 
                visc_eff_ab(i,j,k) = 1.0_prec / (2.0_prec*ATT_ab*(tau_par_ab**2+tau_perp_ab**2)) 

            end do 

        end do 
        end do  


        ! Unstagger from ab-nodes to aa-nodes 
        ! only using contributions from ice covered neighbors
        do j = 1, ny 
        do i = 1, nx 

            im1 = max(i-1,1) 
            ip1 = min(i+1,nx) 
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny) 

            visc_eff(i,j,:) = 0.0 
            wt              = 0.0 

            if (count([H_ice(i,j),H_ice(ip1,j),H_ice(i,jp1),H_ice(ip1,jp1)].eq.0) .eq. 0) then  
                visc_eff(i,j,:) = visc_eff(i,j,:) + visc_eff_ab(i,j,:) 
                wt = wt + 1.0 
            end if 
            
            if (count([H_ice(i,j),H_ice(im1,j),H_ice(im1,jp1),H_ice(i,jp1)].eq.0) .eq. 0) then  
                visc_eff(i,j,:) = visc_eff(i,j,:) + visc_eff_ab(im1,j,:) 
                wt = wt + 1.0 
            end if 

            if (count([H_ice(i,j),H_ice(i,jm1),H_ice(ip1,jm1),H_ice(ip1,j)].eq.0) .eq. 0) then 
                visc_eff(i,j,:) = visc_eff(i,j,:) + visc_eff_ab(i,jm1,:) 
                wt = wt + 1.0 
            end if 
            
            if (count([H_ice(i,j),H_ice(im1,j),H_ice(im1,jm1),H_ice(i,jm1)].eq.0) .eq. 0) then 
                visc_eff(i,j,:) = visc_eff(i,j,:) + visc_eff_ab(im1,jm1,:) 
                wt = wt + 1.0 
            end if 
            
            if (wt .gt. 0.0) then 
                visc_eff(i,j,:) = visc_eff(i,j,:) / wt 
            else 
                ! Just get simple average for ice free points 

                ! Loop over column
                do k = 1, nz_aa 
                    visc_eff(i,j,k) = 0.25_prec*(visc_eff_ab(i,j,k)+visc_eff_ab(im1,j,k) &
                                                +visc_eff_ab(i,jm1,k)+visc_eff_ab(im1,jm1,k))
                end do 

            end if 

        end do 
        end do 
        
        ! Treat the corners to avoid extremes (ab-nodes)
        visc_eff_ab(1,1,:)   = 0.5*(visc_eff_ab(2,1,:)+visc_eff_ab(1,2,:))
        visc_eff_ab(1,ny,:)  = 0.5*(visc_eff_ab(2,ny,:)+visc_eff_ab(1,ny-1,:))
        visc_eff_ab(nx,1,:)  = 0.5*(visc_eff_ab(nx,2,:)+visc_eff_ab(nx-1,1,:))
        visc_eff_ab(nx,ny,:) = 0.5*(visc_eff_ab(nx-1,ny,:)+visc_eff_ab(nx,ny-1,:))

        ! Treat the corners to avoid extremes
        visc_eff(1,1,:)   = 0.5*(visc_eff(2,1,:)+visc_eff(1,2,:))
        visc_eff(1,ny,:)  = 0.5*(visc_eff(2,ny,:)+visc_eff(1,ny-1,:))
        visc_eff(nx,1,:)  = 0.5*(visc_eff(nx,2,:)+visc_eff(nx-1,1,:))
        visc_eff(nx,ny,:) = 0.5*(visc_eff(nx-1,ny,:)+visc_eff(nx,ny-1,:))

        return 

    end subroutine calc_visc_eff_3D

    function calc_staggered_margin(var0,var1,H0,H1) result(var_mid)
        ! Calculate a staggered point but taking upstream point at the margin 

        implicit none 

        real(prec), intent(IN) :: var0 
        real(prec), intent(IN) :: var1 
        real(prec), intent(IN) :: H0 
        real(prec), intent(IN) :: H1
        real(prec) :: var_mid 

        if (H0 .gt. 0.0_prec .and. H1 .eq. 0.0_prec) then 
            ! At the margin 

            var_mid = var0 

        else if (H0 .eq. 0.0_prec .and. H1 .gt. 0.0_prec) then 
            ! At the margin 

            var_mid = var1 

        else 
            ! Inland ice, or ice free, simple average 

            var_mid = 0.5_prec*(var0+var1) 

        end if 

        return 

    end function calc_staggered_margin

    subroutine smooth_visc_eff_int_margin(visc_eff_int,H_ice)

        implicit none 

        real(prec), intent(INOUT) :: visc_eff_int(:,:) 
        real(prec), intent(IN)    :: H_ice(:,:) 
        
        ! Local variables 
        integer    :: i, j, nx, ny 
        integer    :: im1, ip1, jm1, jp1 
        integer    :: n_ice 
        real(prec) :: H_neighb(4) 
        real(prec) :: visc_neighb(4) 
        logical    :: mask_neighb(4) 
        real(prec) :: visc_mean 
        logical    :: is_margin 
        logical    :: is_low_visc 

        real(prec), allocatable :: visc_eff_int_ref(:,:) 

        real(prec), parameter :: visc_ratio_lim = 1e-1      ! visc / mean(visc_neighb) > visc_ratio_lim 

        nx = size(visc_eff_int,1) 
        ny = size(visc_eff_int,2) 

        ! Allocate local arrays 
        allocate(visc_eff_int_ref(nx,ny)) 

        visc_eff_int_ref = visc_eff_int 

        ! Check margin for extreme viscosity values, smooth them as necessary
        do j = 1, ny 
        do i = 1, nx 

            im1 = max(i-1,1) 
            ip1 = min(i+1,nx) 
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny) 

            if (H_ice(i,j) .gt. 0.0) then 
                ! Ice-covered point 

                H_neighb    = [H_ice(im1,j),H_ice(ip1,j),H_ice(i,jm1),H_ice(i,jp1)]
                mask_neighb = H_neighb .gt. 0.0_prec 
                n_ice       = count(mask_neighb)
                is_margin   = n_ice .lt. 4
                
                if (n_ice .gt. 0) then 
                    ! Some neighbors are ice covered too 
                    visc_neighb = [visc_eff_int_ref(im1,j),visc_eff_int_ref(ip1,j),visc_eff_int_ref(i,jm1),visc_eff_int_ref(i,jp1)]
                    visc_mean   = sum(visc_neighb,mask=mask_neighb) / real(n_ice,prec)
                else 
                    ! Just set neighbor mean == to current visc to move on to next point 
                    visc_mean   = visc_eff_int_ref(i,j) 
                end if 

                ! Check if this point has a low viscosity relative to its neighbors 
                is_low_visc = visc_eff_int_ref(i,j) / visc_mean .lt. visc_ratio_lim 

                if (is_margin .and. is_low_visc) then 
                    ! If this is an ice-margin point, and the viscosity is *much*
                    ! lower than ice-covered neighbors, then smooth viscosity here. 

                    visc_eff_int(i,j) = visc_mean 

                end if 

            end if 
                    
        end do 
        end do  

        return 

    end subroutine smooth_visc_eff_int_margin

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
    
end module velocity_l1l2
