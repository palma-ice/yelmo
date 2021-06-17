module velocity_diva

    use yelmo_defs, only  : sp, dp, wp, rho_ice, rho_sw, rho_w, g, TOL_UNDERFLOW
    use yelmo_tools, only : stagger_aa_ab, stagger_aa_ab_ice, stagger_ab_aa_ice, & 
                    integrate_trapezoid1D_1D, integrate_trapezoid1D_pt, minmax

    use basal_dragging 
    use solver_ssa_sico5 
    use velocity_general, only : set_inactive_margins
    
    implicit none 

    type diva_param_class

        character(len=256) :: ssa_lis_opt 
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
    public :: calc_velocity_diva

contains
    
    ! Variables/params needed for gl-flux methods
    ! qq_gl_acx 
    ! qq_gl_acy 
    ! ATT_bar 
    ! 
    
    subroutine calc_velocity_diva(ux,uy,ux_bar,uy_bar,ux_b,uy_b,ux_i,uy_i,taub_acx,taub_acy, &
                                  beta,beta_acx,beta_acy,beta_eff,visc_eff,visc_eff_int,duxdz,duydz, &
                                  ssa_mask_acx,ssa_mask_acy,ssa_err_acx,ssa_err_acy,ssa_iter_now, &
                                  c_bed,taud_acx,taud_acy,H_ice,f_ice,H_grnd,f_grnd, &
                                  f_grnd_acx,f_grnd_acy,ATT,zeta_aa,z_sl,z_bed,dx,dy,n_glen,par)
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
        real(wp), intent(OUT)   :: beta_eff(:,:)      ! [Pa a/m]
        real(wp), intent(INOUT) :: visc_eff(:,:,:)    ! [Pa a]
        real(wp), intent(OUT)   :: visc_eff_int(:,:)  ! [Pa a m]
        real(wp), intent(OUT)   :: duxdz(:,:,:)       ! [1/a]
        real(wp), intent(OUT)   :: duydz(:,:,:)       ! [1/a]
        integer,  intent(OUT)   :: ssa_mask_acx(:,:)  ! [-]
        integer,  intent(OUT)   :: ssa_mask_acy(:,:)  ! [-]
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
        real(wp), intent(IN)    :: dx                 ! [m]
        real(wp), intent(IN)    :: dy                 ! [m]
        real(wp), intent(IN)    :: n_glen 
        type(diva_param_class), intent(IN) :: par       ! List of parameters that should be defined

        ! Local variables 
        integer :: i, j, k, nx, ny, nz_aa, nz_ac, iter 
        logical :: is_converged

        real(wp), allocatable :: ux_bar_nm1(:,:) 
        real(wp), allocatable :: uy_bar_nm1(:,:)  
        real(wp), allocatable :: beta_eff_acx(:,:)
        real(wp), allocatable :: beta_eff_acy(:,:)  
        real(wp), allocatable :: F2(:,:)              ! [Pa^-1 a^-1 m == (Pa a/m)^-1]
        integer,  allocatable :: ssa_mask_acx_ref(:,:)
        integer,  allocatable :: ssa_mask_acy_ref(:,:)

        ! For glf methods 
        logical :: with_glf 
        real(wp), allocatable :: ATT_bar(:,:) 
        real(wp), allocatable :: qq_gl_acx(:,:) 
        real(wp), allocatable :: qq_gl_acy(:,:) 
        
        real(wp) :: L2_norm 
         
        integer  :: ntot, ip1, jp1  

        nx    = size(ux,1)
        ny    = size(ux,2)
        nz_aa = size(ux,3)

        ! Prepare local variables 
        allocate(ux_bar_nm1(nx,ny))
        allocate(uy_bar_nm1(nx,ny))
        allocate(beta_eff_acx(nx,ny))
        allocate(beta_eff_acy(nx,ny))
        allocate(F2(nx,ny))

        allocate(ssa_mask_acx_ref(nx,ny))
        allocate(ssa_mask_acy_ref(nx,ny))

        with_glf = .FALSE. 
        if (trim(par%glf_method) .eq. "power") with_glf = .TRUE.

        ! Store original ssa mask before iterations
        ssa_mask_acx_ref = ssa_mask_acx
        ssa_mask_acy_ref = ssa_mask_acy
            
        ! Initially set error very high 
        ssa_err_acx = 1.0_wp 
        ssa_err_acy = 1.0_wp 
        
        ! Ensure dynamically inactive cells have no velocity at 
        ! outer margins before starting iterations
        call set_inactive_margins(ux_bar,uy_bar,f_ice)

        do iter = 1, par%ssa_iter_max 

            ! Store solution from previous iteration (nm1 == n minus 1) 
            ux_bar_nm1 = ux_bar 
            uy_bar_nm1 = uy_bar 
            
            ! =========================================================================================
            ! Step 1: Calculate fields needed by ssa solver (visc_eff_int, beta_eff)

            ! Calculate the 3D vertical shear fields using viscosity estimated from the previous iteration 
            call calc_vertical_shear_3D(duxdz,duydz,taub_acx,taub_acy,visc_eff,H_ice,f_ice,zeta_aa,par%boundaries)

            ! Calculate 3D effective viscosity
            select case(par%visc_method)

                case(0)
                    ! Impose constant viscosity value 

                    visc_eff = par%visc_const 

                case(1) 
                    ! Calculate effective viscosity, using velocity solution from previous iteration
                    
                    call calc_visc_eff_3D(visc_eff,ux_bar,uy_bar,duxdz,duydz,ATT,H_ice,f_ice, &
                                                                zeta_aa,dx,dy,n_glen,par%eps_0)

                case DEFAULT 

                    write(*,*) "calc_velocity_diva:: Error: visc_method not recognized."
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

            ! Calculate F-integeral (F2) on aa-nodes 
            call calc_F_integral(F2,visc_eff,H_ice,f_ice,zeta_aa,n=2.0_wp)
            
            ! Calculate effective beta 
            call calc_beta_eff(beta_eff,beta,F2,zeta_aa,no_slip=par%no_slip)
            
            ! Stagger beta and beta_eff 
            call stagger_beta(beta_acx,beta_acy,beta,H_ice,f_ice,ux_bar,uy_bar, &
                        f_grnd,f_grnd_acx,f_grnd_acy,par%beta_gl_stag,par%boundaries)
            call stagger_beta(beta_eff_acx,beta_eff_acy,beta_eff,H_ice,f_ice,ux_bar,uy_bar, &
                        f_grnd,f_grnd_acx,f_grnd_acy,par%beta_gl_stag,par%boundaries)
            
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



if (.TRUE.) then 
            if (iter .gt. 1) then
                ! Update ssa mask based on convergence with previous step to reduce area being solved 
                call update_ssa_mask_convergence(ssa_mask_acx,ssa_mask_acy,ssa_err_acx,ssa_err_acy,err_lim=real(1e-5,wp))
                !call update_ssa_mask_convergence(ssa_mask_acx,ssa_mask_acy,ssa_err_acx,ssa_err_acy,err_lim=par%ssa_iter_conv*1e-2)  
            end if 
end if 
            
            ! Call ssa solver
            call calc_vxy_ssa_matrix(ux_bar,uy_bar,L2_norm,beta_eff_acx,beta_eff_acy,visc_eff_int,  &
                                     ssa_mask_acx,ssa_mask_acy,H_ice,f_ice,taud_acx,taud_acy,H_grnd,z_sl, &
                                     z_bed,dx,dy,par%ssa_vel_max,par%boundaries,par%ssa_lis_opt)


            ! Apply relaxation to keep things stable
            call relax_ssa(ux_bar,uy_bar,ux_bar_nm1,uy_bar_nm1,rel=par%ssa_iter_rel)
            
            ! Check for convergence
            is_converged = check_vel_convergence_l2rel(ux_bar,uy_bar,ux_bar_nm1,uy_bar_nm1,ssa_mask_acx.gt.0,     &
                                                       ssa_mask_acy.gt.0,par%ssa_iter_conv,iter,par%ssa_iter_max, &
                                                       par%ssa_write_log,use_L2_norm=.FALSE.,L2_norm=L2_norm)

            ! Calculate an L1 error metric over matrix for diagnostics
            call check_vel_convergence_l1rel_matrix(ssa_err_acx,ssa_err_acy,ux_bar,uy_bar,ux_bar_nm1,uy_bar_nm1)

            ! Store current total iterations for output
            ssa_iter_now = iter 

            ! =========================================================================================
            ! Update additional fields based on output of solver
             
            ! Calculate basal stress 
            call calc_basal_stress(taub_acx,taub_acy,beta_eff_acx,beta_eff_acy,ux_bar,uy_bar)

            ! Calculate basal velocity from depth-averaged solution and basal stress
            call calc_vel_basal(ux_b,uy_b,ux_bar,uy_bar,F2,taub_acx,taub_acy,H_ice,f_ice,par%no_slip,par%boundaries)

            ! Exit iterations if ssa solution has converged
            if (is_converged) exit 
            
        end do 

        ! Iterations are finished, finalize calculations of 3D velocity field 

        ! Calculate the 3D horizontal velocity field
        call calc_vel_horizontal_3D(ux,uy,ux_b,uy_b,taub_acx,taub_acy,visc_eff,H_ice,f_ice,zeta_aa,par%boundaries)

        ! Also calculate the shearing contribution
        do k = 1, nz_aa 
            ux_i(:,:,k) = ux(:,:,k) - ux_b 
            uy_i(:,:,k) = uy(:,:,k) - uy_b 
        end do

        return 

    end subroutine calc_velocity_diva

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
        real(wp) :: H_eff
        real(wp), allocatable :: F1(:,:,:) 
        real(wp), allocatable :: F1_ac(:) 
        
        nx    = size(ux,1)
        ny    = size(ux,2) 
        nz_aa = size(ux,3) 

        allocate(F1(nx,ny,nz_aa))
        allocate(F1_ac(nz_aa))

        ! First calculate F1 array on aa-nodes 
        ! (performing integral before staggering seems to improve result slightly)
        ! Note: L19 define the F1 integral as purely going from the base to the surface,
        ! whereas here F1 is calculated from the base to each point in the vertical. So, 
        ! it is not technically "F1" as defined by L19, Eq. 30, except at the surface.
        do j = 1, ny 
        do i = 1, nx 
            if (f_ice(i,j) .gt. 0.0) then 
                H_eff = H_ice(i,j) / f_ice(i,j) 
            else
                H_eff = H_ice(i,j)
            end if 
            F1(i,j,:) = integrate_trapezoid1D_1D((H_eff/visc_eff(i,j,:))*(1.0-zeta_aa),zeta_aa)
        end do
        end do  

        ! Next calculate 3D horizontal velocity components 
        do j = 1, ny 
        do i = 1, nx 

            ip1 = min(i+1,nx)
            jp1 = min(j+1,ny) 

            ! === x direction ===============================================

            ! Stagger F1 column to ac-nodes 
            if (f_ice(i,j) .eq. 1.0 .and. f_ice(ip1,j) .lt. 1.0) then 
                F1_ac = F1(i,j,:) 
            else if (f_ice(i,j) .lt. 1.0 .and. f_ice(ip1,j) .eq. 1.0) then
                F1_ac = F1(ip1,j,:)
            else 
                F1_ac = 0.5_wp*(F1(i,j,:) + F1(ip1,j,:))
            end if 

            ! Calculate velocity column 
            ux(i,j,:) = ux_b(i,j) + taub_acx(i,j)*F1_ac 

            ! === y direction ===============================================

            ! Stagger F1 column to ac-nodes 
            if (f_ice(i,j) .eq. 1.0 .and. f_ice(i,jp1) .lt. 1.0) then 
                F1_ac = F1(i,j,:) 
            else if (f_ice(i,j) .lt. 1.0 .and. f_ice(i,jp1) .eq. 1.0) then
                F1_ac = F1(i,jp1,:)
            else 
                F1_ac = 0.5_wp*(F1(i,j,:) + F1(i,jp1,:))
            end if 

            ! Calculate velocity column
            uy(i,j,:) = uy_b(i,j) + taub_acy(i,j)*F1_ac  

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

    subroutine calc_vertical_shear_3D(duxdz,duydz,taub_acx,taub_acy,visc_eff,H_ice,f_ice,zeta_aa,boundaries)
        ! Calculate vertical shear terms (L19, Eq. 36)

        implicit none 

        real(wp), intent(OUT) :: duxdz(:,:,:)         ! [1/a],    ac-nodes horizontal, aa-nodes vertical 
        real(wp), intent(OUT) :: duydz(:,:,:)         ! [1/a],    ac-nodes horizontal, aa-nodes vertical 
        real(wp), intent(IN)  :: taub_acx(:,:)        ! [Pa],     ac-nodes
        real(wp), intent(IN)  :: taub_acy(:,:)        ! [Pa],     ac-nodes
        real(wp), intent(IN)  :: visc_eff(:,:,:)      ! [Pa a m], aa-nodes
        real(wp), intent(IN)  :: H_ice(:,:) 
        real(wp), intent(IN)  :: f_ice(:,:) 
        real(wp), intent(IN)  :: zeta_aa(:)           ! [-]
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer :: i, j, k, nx, ny, nz_aa 
        integer :: ip1, jp1 
        real(wp) :: visc_eff_ac

        real(wp), parameter :: visc_min = 1e4_wp 

        nx    = size(duxdz,1)
        ny    = size(duxdz,2)
        nz_aa = size(duxdz,3) 
        
        do k = 1, nz_aa 
        do j = 1, ny
        do i = 1, nx 

            ! Get staggering indices limited to grid size
            ip1 = min(i+1,nx)
            jp1 = min(j+1,ny) 

            ! Calculate shear strain, acx-nodes
            visc_eff_ac  = calc_staggered_margin(visc_eff(i,j,k),visc_eff(ip1,j,k),f_ice(i,j),f_ice(ip1,j))
            visc_eff_ac  = max(visc_eff_ac,visc_min)    ! For safety 
            duxdz(i,j,k) = (taub_acx(i,j)/visc_eff_ac) * (1.0-zeta_aa(k))
            
            ! Calculate shear strain, acy-nodes
            visc_eff_ac  = calc_staggered_margin(visc_eff(i,j,k),visc_eff(i,jp1,k),f_ice(i,j),f_ice(i,jp1))
            visc_eff_ac  = max(visc_eff_ac,visc_min)    ! For safety 
            duydz(i,j,k) = (taub_acy(i,j)/visc_eff_ac) * (1.0-zeta_aa(k))

        end do 
        end do 
        end do 

        ! Apply boundary conditions as needed 
        if (trim(boundaries) .eq. "periodic") then

            duxdz(1,:,:)    = duxdz(nx-2,:,:) 
            duxdz(nx-1,:,:) = duxdz(2,:,:) 
            duxdz(nx,:,:)   = duxdz(3,:,:) 
            duxdz(:,1,:)    = duxdz(:,ny-1,:)
            duxdz(:,ny,:)   = duxdz(:,2,:) 

            duydz(1,:,:)    = duydz(nx-1,:,:) 
            duydz(nx,:,:)   = duydz(2,:,:) 
            duydz(:,1,:)    = duydz(:,ny-2,:)
            duydz(:,ny-1,:) = duydz(:,2,:) 
            duydz(:,ny,:)   = duydz(:,3,:)

        else if (trim(boundaries) .eq. "periodic-x") then 
                
            duxdz(1,:,:)    = duxdz(nx-2,:,:) 
            duxdz(nx-1,:,:) = duxdz(2,:,:) 
            duxdz(nx,:,:)   = duxdz(3,:,:) 
            duxdz(:,1,:)    = duxdz(:,2,:)
            duxdz(:,ny,:)   = duxdz(:,ny-1,:) 

            duydz(1,:,:)    = duydz(nx-1,:,:) 
            duydz(nx,:,:)   = duydz(2,:,:) 
            duydz(:,1,:)    = duydz(:,2,:)
            duydz(:,ny-1,:) = duydz(:,ny-2,:) 
            duydz(:,ny,:)   = duydz(:,ny-1,:)

        else if (trim(boundaries) .eq. "infinite") then 
                
            duxdz(1,:,:)    = duxdz(2,:,:) 
            duxdz(nx-1,:,:) = duxdz(nx-2,:,:) 
            duxdz(nx,:,:)   = duxdz(nx-1,:,:) 
            duxdz(:,1,:)    = duxdz(:,2,:)
            duxdz(:,ny,:)   = duxdz(:,ny-1,:) 

            duydz(1,:,:)    = duydz(2,:,:) 
            duydz(nx,:,:)   = duydz(nx-1,:,:) 
            duydz(:,1,:)    = duydz(:,2,:)
            duydz(:,ny-1,:) = duydz(:,ny-2,:) 
            duydz(:,ny,:)   = duydz(:,ny-1,:)

        end if 

        return 

    end subroutine calc_vertical_shear_3D

    subroutine calc_visc_eff_3D(visc_eff,ux,uy,duxdz,duydz,ATT,H_ice,f_ice,zeta_aa,dx,dy,n_glen,eps_0)
        ! Calculate 3D effective viscosity following L19, Eq. 2
        ! Use of eps_0 ensures non-zero positive viscosity value everywhere 
        ! Note: viscosity is first calculated on ab-nodes, then 
        ! unstaggered back to aa-nodes. This ensures more stability for 
        ! visc_eff (less likely to blow up for low strain rates). 

        implicit none 
        
        real(wp), intent(OUT) :: visc_eff(:,:,:)      ! aa-nodes
        real(wp), intent(IN)  :: ux(:,:)              ! [m/a] Vertically averaged horizontal velocity, x-component
        real(wp), intent(IN)  :: uy(:,:)              ! [m/a] Vertically averaged horizontal velocity, y-component
        real(wp), intent(IN)  :: duxdz(:,:,:)         ! [1/a] Vertical shearing, x-component
        real(wp), intent(IN)  :: duydz(:,:,:)         ! [1/a] Vertical shearing, x-component
        real(wp), intent(IN)  :: ATT(:,:,:)           ! aa-nodes
        real(wp), intent(IN)  :: H_ice(:,:) 
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: zeta_aa(:)           ! Vertical axis (sigma-coordinates from 0 to 1)
        real(wp), intent(IN)  :: dx
        real(wp), intent(IN)  :: dy
        real(wp), intent(IN)  :: n_glen   
        real(wp), intent(IN)  :: eps_0                ! [1/a] Regularization constant (minimum strain rate, ~1e-8)
        
        ! Local variables 
        integer  :: i, j, k
        integer  :: ip1, jp1, im1, jm1 
        integer  :: ip2, jp2, im2, jm2
        integer  :: nx, ny, nz  
        real(wp) :: inv_4dx, inv_4dy 
        real(wp) :: p1, p2, eps_0_sq  
        real(wp) :: eps_sq                          ! [1/a^2]
        real(wp) :: dudx_ab(4)
        real(wp) :: dvdy_ab(4)
        real(wp) :: dudy_ab(4)
        real(wp) :: dvdx_ab(4) 
        real(wp) :: duxdz_ab(4)
        real(wp) :: duydz_ab(4)
        real(wp) :: eps_sq_ab(4)
        real(wp) :: visc_eff_ab(4)
        real(wp) :: ATT_ab(4)
        real(wp) :: wt_ab(4)
        real(wp) :: wt  
        real(wp) :: visc_eff_now 

        real(wp) :: ux_aa, uy_aa 
        logical  :: is_margin 

        nx = size(visc_eff,1)
        ny = size(visc_eff,2)
        nz = size(visc_eff,3)
        
        ! Calculate scaling factors
        inv_4dx = 1.0_wp / (4.0_wp*dx) 
        inv_4dy = 1.0_wp / (4.0_wp*dy) 

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

            im2 = max(i-2,1) 
            ip2 = min(i+2,nx) 
            jm2 = max(j-2,1) 
            jp2 = min(j+2,ny) 

            ! Get ab-node weighting based on whether ice is present 
            wt_ab = 0.0_wp 
            if (count([f_ice(i,j),f_ice(ip1,j),f_ice(i,jp1),f_ice(ip1,jp1)].lt.1.0_wp) .eq. 0) then 
                wt_ab(1) = 1.0_wp
            end if
            if (count([f_ice(i,j),f_ice(im1,j),f_ice(i,jp1),f_ice(im1,jp1)].lt.1.0_wp) .eq. 0) then 
                wt_ab(2) = 1.0_wp 
            end if 
            if (count([f_ice(i,j),f_ice(im1,j),f_ice(i,jm1),f_ice(im1,jm1)].lt.1.0_wp) .eq. 0) then 
                wt_ab(3) = 1.0_wp
            end if 
            if (count([f_ice(i,j),f_ice(ip1,j),f_ice(i,jm1),f_ice(ip1,jm1)].lt.1.0_wp) .eq. 0) then 
                wt_ab(4) = 1.0_wp 
            end if 
            
            wt = sum(wt_ab)
            
            if (f_ice(i,j) .eq. 1.0_wp .and. wt .gt. 0.0_wp) then 
                ! Fully ice-covered point with some fully ice-covered neighbors 

                ! Normalize weighting 
                wt_ab = wt_ab / wt 

                ! === dudx ======================

                dudx_ab(1) = ( (ux(ip1,j) - ux(im1,j)) + (ux(ip1,jp1) - ux(im1,jp1)) ) *inv_4dx
                dudx_ab(2) = ( (ux(i,j)   - ux(im2,j)) + (ux(i,jp1)   - ux(im2,jp1)) ) *inv_4dx
                dudx_ab(3) = ( (ux(i,j)   - ux(im2,j)) + (ux(i,jm1)   - ux(im2,jm1)) ) *inv_4dx
                dudx_ab(4) = ( (ux(ip1,j) - ux(im1,j)) + (ux(ip1,jm1) - ux(im1,jm1)) ) *inv_4dx
                where(abs(dudx_ab) .lt. TOL_UNDERFLOW) dudx_ab = 0.0_wp
                
                ! === dvdy ======================

                dvdy_ab(1) = ( (uy(i,jp1) - uy(i,jm1)) + (uy(ip1,jp1) - uy(ip1,jm1)) ) *inv_4dy 
                dvdy_ab(2) = ( (uy(i,jp1) - uy(i,jm1)) + (uy(im1,jp1) - uy(im1,jm1)) ) *inv_4dy 
                dvdy_ab(3) = ( (uy(i,j)   - uy(i,jm2)) + (uy(im1,j)   - uy(im1,jm2)) ) *inv_4dy 
                dvdy_ab(4) = ( (uy(i,j)   - uy(i,jm2)) + (uy(ip1,j)   - uy(ip1,jm2)) ) *inv_4dy 
                where(abs(dvdy_ab) .lt. TOL_UNDERFLOW) dvdy_ab = 0.0_wp
                

                ! === dudy ======================

                dudy_ab(1) = (ux(i,jp1)   - ux(i,j))     / dy 
                dudy_ab(2) = (ux(im1,jp1) - ux(im1,j))   / dy 
                dudy_ab(3) = (ux(im1,j)   - ux(im1,jm1)) / dy 
                dudy_ab(4) = (ux(i,j)     - ux(i,jm1))   / dy 
                where(abs(dudy_ab) .lt. TOL_UNDERFLOW) dudy_ab = 0.0_wp
                
                ! === dvdx ======================

                dvdx_ab(1) = (uy(ip1,j)   - uy(i,j))     / dx 
                dvdx_ab(2) = (uy(i,j)     - uy(im1,j))   / dx 
                dvdx_ab(3) = (uy(i,jm1)   - uy(im1,jm1)) / dx
                dvdx_ab(4) = (uy(ip1,jm1) - uy(i,jm1))   / dx
                where(abs(dvdx_ab) .lt. TOL_UNDERFLOW) dvdx_ab = 0.0_wp

                ! Loop over column
                do k = 1, nz 

                    duxdz_ab(1) = 0.5_wp*(duxdz(i,j,k)   + duxdz(i,jp1,k))
                    duxdz_ab(2) = 0.5_wp*(duxdz(im1,j,k) + duxdz(im1,jp1,k))
                    duxdz_ab(3) = 0.5_wp*(duxdz(im1,j,k) + duxdz(im1,jm1,k))
                    duxdz_ab(4) = 0.5_wp*(duxdz(i,j,k)   + duxdz(i,jm1,k))
                    where(abs(duxdz_ab) .lt. TOL_UNDERFLOW) duxdz_ab = 0.0_wp

                    duydz_ab(1) = 0.5_wp*(duydz(i,j,k)   + duydz(ip1,j,k))
                    duydz_ab(2) = 0.5_wp*(duydz(i,j,k)   + duydz(im1,j,k))
                    duydz_ab(3) = 0.5_wp*(duydz(i,jm1,k) + duydz(im1,jm1,k))
                    duydz_ab(4) = 0.5_wp*(duydz(i,jm1,k) + duydz(ip1,jm1,k))
                    where(abs(duydz_ab) .lt. TOL_UNDERFLOW) duydz_ab = 0.0_wp

                    ! Calculate the total effective strain rate from L19, Eq. 21 
                    eps_sq_ab = dudx_ab**2 + dvdy_ab**2 + dudx_ab*dvdy_ab + 0.25_wp*(dudy_ab+dvdx_ab)**2 &
                                + 0.25_wp*duxdz_ab**2 + 0.25_wp*duydz_ab**2 + eps_0_sq
                    
                    ! Get the rate factor on ab-nodes too
                    ATT_ab(1) = 0.25_wp*(ATT(i,j,k)+ATT(ip1,j,k)+ATT(i,jp1,k)+ATT(ip1,jp1,k)) 
                    ATT_ab(2) = 0.25_wp*(ATT(i,j,k)+ATT(im1,j,k)+ATT(i,jp1,k)+ATT(im1,jp1,k)) 
                    ATT_ab(3) = 0.25_wp*(ATT(i,j,k)+ATT(im1,j,k)+ATT(i,jm1,k)+ATT(im1,jm1,k)) 
                    ATT_ab(4) = 0.25_wp*(ATT(i,j,k)+ATT(ip1,j,k)+ATT(i,jm1,k)+ATT(ip1,jm1,k)) 
                    
                    ! ajr:testing
                    ATT_ab = ATT(i,j,k)

                    ! Calculate effective viscosity on ab-nodes
                    visc_eff_ab = 0.5_wp*(eps_sq_ab)**(p1) * ATT_ab**(p2)

                    ! Calcualte effective viscosity on aa-nodes
                    visc_eff(i,j,k) = sum(wt_ab*visc_eff_ab)

                end do 

            else  
                ! Get simple effective viscosity ignoring staggering

                do k = 1, nz
                    visc_eff(i,j,k) = 0.5_wp*(eps_0_sq)**(p1) * ATT(i,j,k)**(p2)
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

    subroutine smooth_visc_eff_int_margin(visc_eff_int,H_ice)
        ! ajr: not currently, used, needs to be updated with f_ice if 
        ! it will be used in the future...
        ! (most likely we should delete this routine and solve visc_eff
        ! properly at the margin) 

        implicit none 

        real(wp), intent(INOUT) :: visc_eff_int(:,:) 
        real(wp), intent(IN)    :: H_ice(:,:) 
        
        ! Local variables 
        integer    :: i, j, nx, ny 
        integer    :: im1, ip1, jm1, jp1 
        integer    :: n_ice 
        real(wp) :: H_neighb(4) 
        real(wp) :: visc_neighb(4) 
        logical    :: mask_neighb(4) 
        real(wp) :: visc_mean 
        logical    :: is_margin 
        logical    :: is_low_visc 

        real(wp), allocatable :: visc_eff_int_ref(:,:) 

        real(wp), parameter :: visc_ratio_lim = 1e-1      ! visc / mean(visc_neighb) > visc_ratio_lim 

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
                mask_neighb = H_neighb .gt. 0.0_wp 
                n_ice       = count(mask_neighb)
                is_margin   = n_ice .lt. 4
                
                if (n_ice .gt. 0) then 
                    ! Some neighbors are ice covered too 
                    visc_neighb = [visc_eff_int_ref(im1,j),visc_eff_int_ref(ip1,j),visc_eff_int_ref(i,jm1),visc_eff_int_ref(i,jp1)]
                    visc_mean   = sum(visc_neighb,mask=mask_neighb) / real(n_ice,wp)
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

    subroutine calc_F_integral(F_int,visc,H_ice,f_ice,zeta_aa,n)
        ! Useful integrals, following Arthern et al. (2015) Eq. 7,
        ! and Lipscomb et al. (2019), Eq. 30
        ! F_n = int_zb_zs{ 1/visc * ((s-z)/H)**n dz}

        implicit none 

        real(wp), intent(OUT) :: F_int(:,:) 
        real(wp), intent(IN)  :: visc(:,:,:)
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: zeta_aa(:)
        real(wp), intent(IN)  :: n  

        ! Local variables 
        integer :: i, j, nx, ny, nz_aa, np
        integer :: im1, jm1, ip1, jp1
        real(wp) :: H_eff  
        real(wp) :: F_int_min 
        real(wp), parameter :: visc_min = 1e4_wp

        nx    = size(visc,1)
        ny    = size(visc,2) 
        nz_aa = size(visc,3)

        ! Determine the minimum value of F_int, to assign when H_ice == 0,
        ! since F_int should be nonzero everywhere for numerics
        F_int_min = integrate_trapezoid1D_pt((1.0_wp/visc_min)*(1.0_wp-zeta_aa)**n,zeta_aa)

        ! Initially set F_int to minimum value everywhere 
        F_int = F_int_min

        ! Vertically integrate at each point
        do j = 1, ny 
        do i = 1, nx

            im1 = max(i-1,1)
            jm1 = max(j-1,1)
            ip1 = min(i+1,nx)
            jp1 = min(j+1,ny)

            if (f_ice(i,j) .gt. 0.0) then 
                ! Viscosity should be nonzero here, perform integration 

                H_eff = H_ice(i,j) / f_ice(i,j) 
                F_int(i,j) = integrate_trapezoid1D_pt((H_eff/visc(i,j,:) )*(1.0_wp-zeta_aa)**n,zeta_aa)

            else 

                F_int(i,j) = F_int_min

            end if 

        end do 
        end do 

        return

    end subroutine calc_F_integral
    
    subroutine calc_beta_eff(beta_eff,beta,F2,zeta_aa,no_slip)
        ! Calculate the depth-averaged horizontal velocity (ux_bar,uy_bar)

        ! Note: L19 staggers the F-integral F2, then solves for beta 

        implicit none 
        
        real(wp), intent(OUT) :: beta_eff(:,:)    ! aa-nodes
        real(wp), intent(IN)  :: beta(:,:)        ! aa-nodes
        real(wp), intent(IN)  :: F2(:,:)          ! aa-nodes
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

    subroutine calc_vel_basal(ux_b,uy_b,ux_bar,uy_bar,F2,taub_acx,taub_acy,H_ice,f_ice,no_slip,boundaries)
        ! Calculate basal sliding following Goldberg (2011), Eq. 34
        ! (or it can also be obtained from L19, Eq. 32 given ub*beta=taub)

        implicit none
        
        real(wp), intent(OUT) :: ux_b(:,:) 
        real(wp), intent(OUT) :: uy_b(:,:)
        real(wp), intent(IN)  :: ux_bar(:,:) 
        real(wp), intent(IN)  :: uy_bar(:,:)
        real(wp), intent(IN)  :: F2(:,:)
        real(wp), intent(IN)  :: taub_acx(:,:) 
        real(wp), intent(IN)  :: taub_acy(:,:)
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: f_ice(:,:)
        logical,    intent(IN)  :: no_slip
        character(len=*), intent(IN) :: boundaries 

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

                ip1 = min(i+1,nx)
                jp1 = min(j+1,ny)

                ! ==== x-direction =====

                ! Stagger the F2 integral to the ac-nodes
                if (f_ice(i,j) .eq. 1.0 .and. f_ice(ip1,j) .lt. 1.0) then 
                    F2_ac = F2(i,j) 
                else if (f_ice(i,j) .lt. 1.0 .and. f_ice(ip1,j) .eq. 1.0) then
                    F2_ac = F2(ip1,j)
                else 
                    F2_ac = 0.5_wp*(F2(i,j) + F2(ip1,j))
                end if 

                ! Calculate basal velocity component 
                ux_b(i,j) = ux_bar(i,j) - taub_acx(i,j)*F2_ac 

                ! ==== y-direction =====
                
                ! Stagger the F2 integral to the ac-nodes
                if (f_ice(i,j) .eq. 1.0 .and. f_ice(i,jp1) .lt. 1.0) then 
                    F2_ac = F2(i,j) 
                else if (f_ice(i,j) .lt. 1.0 .and. f_ice(i,jp1) .eq. 1.0) then
                    F2_ac = F2(i,jp1)
                else 
                    F2_ac = 0.5_wp*(F2(i,j) + F2(i,jp1))
                end if 
                    
                ! Calculate basal velocity component 
                uy_b(i,j) = uy_bar(i,j) - taub_acy(i,j)*F2_ac 

            end do 
            end do  

            ! Apply boundary conditions as needed 
            if (trim(boundaries) .eq. "periodic") then 

                ux_b(1,:)    = ux_b(nx-2,:) 
                ux_b(nx-1,:) = ux_b(2,:) 
                ux_b(nx,:)   = ux_b(3,:) 
                ux_b(:,1)    = ux_b(:,ny-1)
                ux_b(:,ny)   = ux_b(:,2) 
                
                uy_b(1,:)    = uy_b(nx-1,:) 
                uy_b(nx,:)   = uy_b(2,:) 
                uy_b(:,1)    = uy_b(:,ny-2)
                uy_b(:,ny-1) = uy_b(:,2) 
                uy_b(:,ny)   = uy_b(:,3)

            else if (trim(boundaries) .eq. "periodic-x") then 
                
                ux_b(1,:)    = ux_b(nx-2,:) 
                ux_b(nx-1,:) = ux_b(2,:) 
                ux_b(nx,:)   = ux_b(3,:) 
                ux_b(:,1)    = ux_b(:,2)
                ux_b(:,ny)   = ux_b(:,ny-1) 

                uy_b(1,:)    = uy_b(nx-1,:) 
                uy_b(nx,:)   = uy_b(2,:) 
                uy_b(:,1)    = uy_b(:,2)
                uy_b(:,ny-1) = uy_b(:,ny-2) 
                uy_b(:,ny)   = uy_b(:,ny-1)

            else if (trim(boundaries) .eq. "infinite") then 
                
                ux_b(1,:)    = ux_b(2,:) 
                ux_b(nx-1,:) = ux_b(nx-2,:) 
                ux_b(nx,:)   = ux_b(nx-1,:) 
                ux_b(:,1)    = ux_b(:,2)
                ux_b(:,ny)   = ux_b(:,ny-1) 

                uy_b(1,:)    = uy_b(2,:) 
                uy_b(nx,:)   = uy_b(nx-1,:) 
                uy_b(:,1)    = uy_b(:,2)
                uy_b(:,ny-1) = uy_b(:,ny-2) 
                uy_b(:,ny)   = uy_b(:,ny-1)

            end if 

        end if 

        return
        
    end subroutine calc_vel_basal

    subroutine calc_basal_stress(taub_acx,taub_acy,beta_eff_acx,beta_eff_acy,ux_bar,uy_bar)
        ! Calculate the basal stress resulting from sliding (friction times velocity)
        ! Note: calculated on ac-nodes.
        ! taub [Pa] 
        ! beta [Pa a m-1]
        ! u    [m a-1]
        ! taub = beta*u (here defined with taub in the same direction as u)

        implicit none 

        real(wp), intent(OUT) :: taub_acx(:,:)        ! [Pa] Basal stress (acx nodes)
        real(wp), intent(OUT) :: taub_acy(:,:)        ! [Pa] Basal stress (acy nodes)
        real(wp), intent(IN)  :: beta_eff_acx(:,:)    ! [Pa a m-1] Effective basal friction (acx nodes)
        real(wp), intent(IN)  :: beta_eff_acy(:,:)    ! [Pa a m-1] Effective basal friction (acy nodes)
        real(wp), intent(IN)  :: ux_bar(:,:)          ! [m a-1] depth-ave velocity (acx nodes)
        real(wp), intent(IN)  :: uy_bar(:,:)          ! [m a-1] depth-ave velocity (acy nodes)
        
        ! Local variables 
        integer :: i, j, nx, ny 

        nx = size(taub_acx,1)
        ny = size(taub_acy,2) 

        do j = 1, ny 
        do i = 1, nx 

            taub_acx(i,j) = beta_eff_acx(i,j) * ux_bar(i,j) 
            if (abs(taub_acx(i,j)) .lt. TOL_UNDERFLOW) taub_acx(i,j) = 0.0_wp 

            taub_acy(i,j) = beta_eff_acy(i,j) * uy_bar(i,j) 
            if (abs(taub_acy(i,j)) .lt. TOL_UNDERFLOW) taub_acy(i,j) = 0.0_wp 

        end do 
        end do  

        return 

    end subroutine calc_basal_stress

    elemental subroutine limit_vel(u,u_lim)
        ! Apply a velocity limit (for stability)

        implicit none 

        real(wp), intent(INOUT) :: u 
        real(wp), intent(IN)    :: u_lim

        real(wp), parameter :: tol = 1e-10
        
        u = min(u, u_lim)
        u = max(u,-u_lim)

        ! Also avoid underflow errors 
        if (abs(u) .lt. tol) u = 0.0 

        return 

    end subroutine limit_vel
    
end module velocity_diva
