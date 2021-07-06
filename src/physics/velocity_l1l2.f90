module velocity_l1l2

    use yelmo_defs ,only  : wp, prec, rho_ice, rho_sw, rho_w, g
    use yelmo_tools, only : stagger_aa_ab, stagger_aa_ab_ice, stagger_ab_aa_ice, &
                    calc_vertical_integrated_2D, & 
                    integrate_trapezoid1D_1D, integrate_trapezoid1D_pt, minmax

    use basal_dragging 
    use solver_ssa_sico5 
    use velocity_general, only : set_inactive_margins
    
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
                                  c_bed,taud_acx,taud_acy,H_ice,f_ice,H_grnd,f_grnd, &
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
        real(prec), intent(IN)    :: f_ice(:,:)         ! [--]
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

        ! Store original ssa mask before iterations
        ssa_mask_acx_ref = ssa_mask_acx
        ssa_mask_acy_ref = ssa_mask_acy
            
        ! Initially set error very high 
        ssa_err_acx = 1.0_prec 
        ssa_err_acy = 1.0_prec 
        
        ! Ensure dynamically inactive cells have no velocity at 
        ! outer margins before starting iterations
        call set_inactive_margins(ux_b,uy_b,f_ice)

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
                                                H_ice,f_ice,zeta_aa,dx,dy,n_glen,par%eps_0,par%boundaries)

                case DEFAULT 

                    write(*,*) "calc_velocity_l1l2:: Error: visc_method not recognized."
                    write(*,*) "visc_method = ", par%visc_method 
                    stop 
                    
            end select
                  
            ! Calculate depth-integrated effective viscosity
            ! Note L19 uses eta_bar*H in the ssa equation. Yelmo uses eta_int=eta_bar*H directly.
            call calc_visc_eff_int(visc_eff_int,visc_eff,H_ice,f_ice,zeta_aa)
            
            ! Smooth the viscosity at the ice margins if it is too low
            ! to avoid singularities (mainly for EISMINT/dome experiments)
            !call smooth_visc_eff_int_margin(visc_eff_int,H_ice)

            ! Calculate beta (at the ice base)
            call calc_beta(beta,c_bed,ux_b,uy_b,H_ice,f_ice,H_grnd,f_grnd,z_bed,z_sl,par%beta_method, &
                                par%beta_const,par%beta_q,par%beta_u0,par%beta_gl_scale,par%beta_gl_f, &
                                par%H_grnd_lim,par%beta_min,par%boundaries)

            ! Stagger beta and beta_eff 
            call stagger_beta(beta_acx,beta_acy,beta,H_ice,f_ice,ux_b,uy_b, &
                        f_grnd,f_grnd_acx,f_grnd_acy,par%beta_gl_stag,par%boundaries)

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
                                         ssa_mask_acx,ssa_mask_acy,H_ice,f_ice,taud_acx,taud_acy,H_grnd,z_sl, &
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
        call calc_vel_horizontal_3D(ux,uy,ux_b,uy_b,taud_acx,taud_acy,visc_eff_ab,ATT,H_ice, &
                                            zeta_aa,dx,dy,n_glen,par%eps_0,par%boundaries)
        
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

    subroutine calc_vel_horizontal_3D(ux,uy,ux_b,uy_b,taud_acx,taud_acy,visc_eff_ab,ATT,H_ice,zeta_aa,dx,dy,n_glen,eps_0,boundaries)
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

        ! Calculate scaling factors
        inv_4dx = 1.0_prec / (4.0_prec*dx) 
        inv_4dy = 1.0_prec / (4.0_prec*dy) 

        ! Calculate an exponent 
        p1 = (n_glen-1.0_wp)/2.0_wp

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

        ux = 0.0 
        uy = 0.0 

        ! Assign basal velocity value 
        ux(:,:,1)    = ux_b 
        uy(:,:,1)    = uy_b 
        fact_ab(:,:) = 0.0_prec 

        ! Loop over layers starting from first layer above the base to surface 
        do k = 2, nz_aa

            if (k .eq. nz_aa) then 
                zeta_ac1 = zeta_aa(nz_aa)
            else 
                zeta_ac1 = 0.5_prec*(zeta_aa(k+1)+zeta_aa(k))
            end if 
            zeta_ac0 = 0.5_prec*(zeta_aa(k)+zeta_aa(k-1))
            
            dzeta = zeta_ac1 - zeta_ac0 

            ! Calculate tau_perp, tau_eff and factor to calculate velocities,
            ! all on ab-nodes 
            do i = 1, nx 
            do j = 1, ny 

                im1 = max(i-1,1) 
                ip1 = min(i+1,nx) 
                jm1 = max(j-1,1) 
                jp1 = min(j+1,ny) 

                ! Calculate effective stress 
                tau_xz_ab = 0.5_prec*(tau_xz(i,j,k)+tau_xz(i,jp1,k))
                tau_yz_ab = 0.5_prec*(tau_yz(i,j,k)+tau_yz(ip1,j,k))
                tau_eff_sq_ab = tau_par_ab(i,j,k)**2 + tau_xz_ab**2 + tau_yz_ab**2

                ! Calculate factor to get velocity components
                ATT_ab   = 0.25_prec*(ATT(i,j,k)+ATT(ip1,j,k)+ATT(i,jp1,k)+ATT(ip1,jp1,k))
                !depth_ab = H_ice_ab(i,j)*(1.0_prec-zeta_aa(k))
                depth_ab = (1.0_prec-zeta_aa(k))

                !fact_ab(i,j) = 2.0_prec * ATT_ab * depth_ab * tau_eff_sq_ab**p1 
                ! fact_ab(i,j) = 2.0_prec * ATT_ab * tau_eff_sq_ab**p1 * (dzeta*H_ice_ab(i,j))
                
                ! fact = 2.d0 * stagflwa(i,j) * tau_eff_sq**((gn-1.d0)/2.d0) * (sigma(k+1) - sigma(k))*stagthck(i,j)


                fact_ab(i,j) = fact_ab(i,j) &
                    - 2.0_prec * ATT_ab * depth_ab * tau_eff_sq_ab * (dzeta*H_ice_ab(i,j))
                
            end do 
            end do 

            ! Calculate 3D horizontal velocity components on acx/acy nodes
            do i = 1, nx 
            do j = 1, ny 

                im1 = max(i-1,1) 
                jm1 = max(j-1,1) 
                
                ! stagger factor to acx-nodes and calculate velocity
                fact_ac   = 0.5_prec*(fact_ab(i,j)+fact_ab(i,jm1))
                ux(i,j,k) = ux(i,j,1) + fact_ac*taud_acx(i,j)

                ! stagger factor to acy-nodes and calculate velocity
                fact_ac   = 0.5_prec*(fact_ab(i,j)+fact_ab(im1,j))
                uy(i,j,k) = uy(i,j,1) + fact_ac*taud_acy(i,j)

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

    subroutine calc_visc_eff_3D(visc_eff,visc_eff_ab,ux_b,uy_b,taud_acx,taud_acy,ATT,H_ice,f_ice,zeta_aa,dx,dy,n_glen,eps_0,boundaries)
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
        real(prec) :: tau_eff_sq_ab, ATT_ab 
        real(prec) :: dudx_ab, dvdy_ab, dudy_ab, dvdx_ab
        
        real(prec) :: eps_par_sq, eps_par_ab 
        real(prec) :: eps_0_sq 
        real(prec) :: taud_ab, tau_par_ab, tau_perp_ab 
        real(prec) :: a, b, c, rootA, rootB, np  
        real(prec) :: wt 

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
            ! effstrain = A * (tau_parallel^2 + tau_perp^2)^{(n-1)/2} * tau_parallel
            ! y = A * (x^2 + tau^2)^{(n-1)/2} * x 

            do k = 1, nz_aa 
                
                ATT_ab = 0.25_prec*(ATT(i,j,k)+ATT(ip1,j,k)+ATT(i,jp1,k)+ATT(ip1,jp1,k))
                
                tau_perp_ab = taud_ab*(1.0_prec-zeta_aa(k))

if (.TRUE.) then 
    ! CISM root equation for n_glen=3 only 
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

                !write(*,*) "tau_par_ab: ", a, b, tau_par_ab  
else 
    ! Root finding code (more expensive, but works for ISMIPHOM)
    ! Crashed for a random Antarctica simulation - needs testing! 

                a  = tau_perp_ab**2 
                b  = eps_par_ab / ATT_ab 
                np = (n_glen-1)/2.0_wp 

                !write(*,*) 'newton', a, b, np 
                call solve_secant(tau_par_ab,n_iter,10e3,1e-3,50,funY,a,b,np,.FALSE.) 
                !call solve_newton(tau_par_ab,n_iter,10e3,1e-3,50,funY,funYp,a,b,np,.FALSE.)
                !stop

end if 

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

            if (f_ice(i,j) .eq. 1.0) then
                ! Ice-covered point. 
                ! Only use contributions from ice-covered neighbors 

                if (count([f_ice(i,j),f_ice(ip1,j),f_ice(i,jp1),f_ice(ip1,jp1)].lt.1.0) .eq. 0) then  
                    visc_eff(i,j,:) = visc_eff(i,j,:) + visc_eff_ab(i,j,:) 
                    wt = wt + 1.0 
                end if 
                
                if (count([f_ice(i,j),f_ice(im1,j),f_ice(im1,jp1),f_ice(i,jp1)].lt.1.0) .eq. 0) then  
                    visc_eff(i,j,:) = visc_eff(i,j,:) + visc_eff_ab(im1,j,:) 
                    wt = wt + 1.0 
                end if 

                if (count([f_ice(i,j),f_ice(i,jm1),f_ice(ip1,jm1),f_ice(ip1,j)].lt.1.0) .eq. 0) then 
                    visc_eff(i,j,:) = visc_eff(i,j,:) + visc_eff_ab(i,jm1,:) 
                    wt = wt + 1.0 
                end if 
                
                if (count([f_ice(i,j),f_ice(im1,j),f_ice(im1,jm1),f_ice(i,jm1)].lt.1.0) .eq. 0) then 
                    visc_eff(i,j,:) = visc_eff(i,j,:) + visc_eff_ab(im1,jm1,:) 
                    wt = wt + 1.0 
                end if 
            
            else
                ! Ice-free point.
                ! Only use contributions from ice-free neighbors 

                if (count([f_ice(i,j),f_ice(ip1,j),f_ice(i,jp1),f_ice(ip1,jp1)].eq.1.0) .eq. 0) then  
                    visc_eff(i,j,:) = visc_eff(i,j,:) + visc_eff_ab(i,j,:) 
                    wt = wt + 1.0 
                end if 
                
                if (count([f_ice(i,j),f_ice(im1,j),f_ice(im1,jp1),f_ice(i,jp1)].eq.1.0) .eq. 0) then  
                    visc_eff(i,j,:) = visc_eff(i,j,:) + visc_eff_ab(im1,j,:) 
                    wt = wt + 1.0 
                end if 

                if (count([f_ice(i,j),f_ice(i,jm1),f_ice(ip1,jm1),f_ice(ip1,j)].eq.1.0) .eq. 0) then 
                    visc_eff(i,j,:) = visc_eff(i,j,:) + visc_eff_ab(i,jm1,:) 
                    wt = wt + 1.0 
                end if 
                
                if (count([f_ice(i,j),f_ice(im1,j),f_ice(im1,jm1),f_ice(i,jm1)].eq.1.0) .eq. 0) then 
                    visc_eff(i,j,:) = visc_eff(i,j,:) + visc_eff_ab(im1,jm1,:) 
                    wt = wt + 1.0 
                end if 
            
            end if 

            if (wt .gt. 0.0) then 
                ! Get the weighted mean of the viscosity for this aa-node 

                visc_eff(i,j,:) = visc_eff(i,j,:) / wt 

            else 
                ! Just get simple average for safety
                ! (this case should not occur)

                ! Loop over column
                do k = 1, nz_aa 
                    visc_eff(i,j,k) = 0.25_wp*(visc_eff_ab(i,j,k)+visc_eff_ab(im1,j,k) &
                                                +visc_eff_ab(i,jm1,k)+visc_eff_ab(im1,jm1,k))
                end do 

            end if 

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

        subroutine calc_visc_eff_int(visc_eff_int,visc_eff,H_ice,f_ice,zeta_aa)

        implicit none 

        real(wp), intent(OUT) :: visc_eff_int(:,:)
        real(wp), intent(IN)  :: visc_eff(:,:,:)
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: zeta_aa(:)

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
                visc_eff_int(i,j) = visc_eff_mean 
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
