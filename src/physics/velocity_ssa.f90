module velocity_ssa 
    ! This module solves for the ice sheet's horizontal velocity field (ux,uy)
    ! using the shallow-shelf approximation (SSA). 

    use yelmo_defs ,only  : sp, dp, wp, prec, tol_underflow, pi, rho_ice, rho_sw, rho_w, g
    use yelmo_tools, only : stagger_aa_ab, stagger_aa_ab_ice, stagger_ab_aa_ice, & 
                    stagger_nodes_aa_ab_ice, stagger_nodes_acx_ab_ice, stagger_nodes_acy_ab_ice, &
                    staggerdiff_nodes_acx_ab_ice, staggerdiff_nodes_acy_ab_ice, &
                    staggerdiffcross_nodes_acx_ab_ice, staggerdiffcross_nodes_acy_ab_ice, &
                    integrate_trapezoid1D_1D, integrate_trapezoid1D_pt, minmax, &
                    set_boundaries_2D_aa, set_boundaries_3D_aa

    use basal_dragging 
    use solver_ssa_sico5 
    use velocity_general, only : set_inactive_margins, &
                        picard_calc_error, picard_calc_error_angle, picard_relax, &
                        picard_calc_convergence_l2 
    implicit none 

    type ssa_param_class

        character(len=256) :: ssa_lis_opt
        character(len=256) :: ssa_lateral_bc 
        character(len=256) :: boundaries  
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
    public :: ssa_param_class 
    public :: calc_velocity_ssa
      
contains 

    subroutine calc_velocity_ssa(ux_b,uy_b,taub_acx,taub_acy,visc_eff,visc_eff_int,ssa_mask_acx,ssa_mask_acy, &
                                  ssa_err_acx,ssa_err_acy,ssa_iter_now,beta,beta_acx,beta_acy,c_bed,taud_acx,taud_acy,H_ice, &
                                  f_ice,H_grnd,f_grnd,f_grnd_acx,f_grnd_acy,ATT,zeta_aa,z_sl,z_bed,z_srf,dx,dy,n_glen,par)
        ! This subroutine is used to solve the horizontal velocity system (ux,uy)
        ! following the SSA solution

        implicit none 

        real(prec), intent(INOUT) :: ux_b(:,:)          ! [m/a]
        real(prec), intent(INOUT) :: uy_b(:,:)          ! [m/a]
        real(prec), intent(INOUT) :: taub_acx(:,:)      ! [Pa]
        real(prec), intent(INOUT) :: taub_acy(:,:)      ! [Pa]
        real(prec), intent(OUT)   :: visc_eff(:,:,:)    ! [Pa a]
        real(prec), intent(OUT)   :: visc_eff_int(:,:)  ! [Pa a m]
        integer,    intent(OUT)   :: ssa_mask_acx(:,:)  ! [-]
        integer,    intent(OUT)   :: ssa_mask_acy(:,:)  ! [-]
        real(prec), intent(OUT)   :: ssa_err_acx(:,:)
        real(prec), intent(OUT)   :: ssa_err_acy(:,:)
        integer,    intent(OUT)   :: ssa_iter_now 
        real(prec), intent(INOUT) :: beta(:,:)          ! [Pa a/m]
        real(prec), intent(INOUT) :: beta_acx(:,:)      ! [Pa a/m]
        real(prec), intent(INOUT) :: beta_acy(:,:)      ! [Pa a/m]
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
        real(prec), intent(IN)    :: z_srf(:,:)         ! [m]
        real(prec), intent(IN)    :: dx                 ! [m]
        real(prec), intent(IN)    :: dy                 ! [m]
        real(prec), intent(IN)    :: n_glen 
        type(ssa_param_class), intent(IN) :: par       ! List of parameters that should be defined

        ! Local variables 
        integer :: i, j, k, nx, ny, nz_aa, nz_ac, iter
        logical :: is_converged

        real(prec), allocatable :: ux_b_nm1(:,:) 
        real(prec), allocatable :: uy_b_nm1(:,:)    
        integer,    allocatable :: ssa_mask_acx_ref(:,:)
        integer,    allocatable :: ssa_mask_acy_ref(:,:)

        real(wp), allocatable :: corr_nm1(:) 
        real(wp), allocatable :: corr_nm2(:) 
        
        real(wp) :: corr_theta
        real(wp) :: corr_rel 
        real(wp) :: L2_norm 
        real(wp) :: ssa_resid 

        logical, parameter :: write_ssa_diagnostics      = .FALSE. 
        logical, parameter :: write_ssa_diagnostics_stop = .FALSE.   ! Stop simulation after completing iterations?

        nx    = size(ux_b,1)
        ny    = size(ux_b,2)
        nz_aa = size(zeta_aa,1)

        ! Prepare local variables 
        allocate(ux_b_nm1(nx,ny))
        allocate(uy_b_nm1(nx,ny))
        allocate(ssa_mask_acx_ref(nx,ny))
        allocate(ssa_mask_acy_ref(nx,ny))

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
        call set_inactive_margins(ux_b,uy_b,f_ice)

        if (write_ssa_diagnostics) then 
            call ssa_diagnostics_write_init("yelmo_ssa.nc",nx,ny,time_init=1.0_wp)
        end if 

        do iter = 1, par%ssa_iter_max 

            ! Store solution from previous iteration (nm1 == n minus 1) 
            ux_b_nm1 = ux_b 
            uy_b_nm1 = uy_b 
            

            ! =========================================================================================
            ! Step 1: Calculate fields needed by ssa solver (visc_eff_int, beta)

            ! Calculate 3D effective viscosity
            select case(par%visc_method)

                case(0)
                    ! Impose constant viscosity value 

                    visc_eff = par%visc_const 

                case(1) 
                    ! Calculate 3D effective viscosity, using velocity solution from previous iteration
                    ! Use default staggering stencil (to ab-node quadrature points, then average to aa-nodes)

                    call calc_visc_eff_3D(visc_eff,ux_b,uy_b,ATT,H_ice,f_ice,zeta_aa, &
                                                    dx,dy,n_glen,par%eps_0,par%boundaries)

                case(2) 
                    ! Calculate 3D effective viscosity, using velocity solution from previous iteration
                    ! Use minimal staggering stencil (directly to aa-nodes)

                    call calc_visc_eff_3D_aa(visc_eff,ux_b,uy_b,ATT,H_ice,f_ice,zeta_aa, &
                                                    dx,dy,n_glen,par%eps_0,par%boundaries)

                case DEFAULT 

                    write(*,*) "calc_velocity_ssa:: Error: visc_method not recognized."
                    write(*,*) "visc_method = ", par%visc_method 
                    stop 

            end select
                    
            ! Calculate depth-integrated effective viscosity
            ! Note L19 uses eta_bar*H in the ssa equation. Yelmo uses eta_int=eta_bar*H directly.
            call calc_visc_eff_int(visc_eff_int,visc_eff,H_ice,f_ice,zeta_aa,par%boundaries)
            
            ! Calculate beta (at the ice base)
            call calc_beta(beta,c_bed,ux_b,uy_b,H_ice,f_ice,H_grnd,f_grnd,z_bed,z_sl,par%beta_method, &
                                par%beta_const,par%beta_q,par%beta_u0,par%beta_gl_scale,par%beta_gl_f, &
                                par%H_grnd_lim,par%beta_min,par%boundaries)

            ! Stagger beta
            call stagger_beta(beta_acx,beta_acy,beta,H_ice,f_ice,ux_b,uy_b, &
                        f_grnd,f_grnd_acx,f_grnd_acy,par%beta_gl_stag,par%beta_min,par%boundaries)


            ! =========================================================================================
            ! Step 2: Call the SSA solver to obtain new estimate of ux_b/uy_b

if (.FALSE.) then 
            if (iter .gt. 1) then
                ! Update ssa mask based on convergence with previous step to reduce area being solved 
                call update_ssa_mask_convergence(ssa_mask_acx,ssa_mask_acy,ssa_err_acx,ssa_err_acy,err_lim=real(1e-5,prec))
                !call update_ssa_mask_convergence(ssa_mask_acx,ssa_mask_acy,ssa_err_acx,ssa_err_acy,err_lim=par%ssa_iter_conv*1e-2)  
            end if 
end if 
            
            ! Call ssa solver
            call calc_vxy_ssa_matrix(ux_b,uy_b,L2_norm,beta_acx,beta_acy,visc_eff_int,  &
                                ssa_mask_acx,ssa_mask_acy,H_ice,f_ice,taud_acx,taud_acy,H_grnd,z_sl,z_bed, &
                                z_srf,dx,dy,par%ssa_vel_max,par%boundaries,par%ssa_lateral_bc,par%ssa_lis_opt)


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
            
            ! Check for convergence
            ! is_converged = check_vel_convergence_l2rel(ux_b,uy_b,ux_b_nm1,uy_b_nm1,ssa_mask_acx.gt.0,     &
            !                                            ssa_mask_acy.gt.0,par%ssa_iter_conv,iter,par%ssa_iter_max, &
            !                                            par%ssa_write_log,use_L2_norm=.FALSE.,L2_norm=L2_norm)
            call picard_calc_convergence_l2(is_converged,ssa_resid,ux_b,uy_b,ux_b_nm1,uy_b_nm1, &
                                                ssa_mask_acx.gt.0,ssa_mask_acy.gt.0,par%ssa_iter_conv,  &
                                                iter,par%ssa_iter_max,par%ssa_write_log)


            ! Calculate an L1 error metric over matrix for diagnostics
            call check_vel_convergence_l1rel_matrix(ssa_err_acx,ssa_err_acy,ux_b,uy_b,ux_b_nm1,uy_b_nm1)

            
            ! Store current total iterations for output
            ssa_iter_now = iter 

            if (write_ssa_diagnostics) then  
                call ssa_diagnostics_write_step("yelmo_ssa.nc",ux_b,uy_b,L2_norm,beta_acx,beta_acy,visc_eff_int, &
                                        ssa_mask_acx,ssa_mask_acy,ssa_err_acx,ssa_err_acy,H_ice,f_ice,taud_acx,taud_acy, &
                                        H_grnd,z_sl,z_bed,z_srf,ux_b_nm1,uy_b_nm1,time=real(iter,wp))    
            end if 

            ! =========================================================================================
            ! Update additional fields based on output of solver
             

            ! Exit iterations if ssa solution has converged
            if (is_converged) exit 
            
        end do 

        ! Iterations are finished, finalize calculations

        if (write_ssa_diagnostics .and. write_ssa_diagnostics_stop) then 
            stop 
        end if 

        ! Diagnose basal stress 
        call calc_basal_stress(taub_acx,taub_acy,beta_acx,beta_acy,ux_b,uy_b)
        
        return 

    end subroutine calc_velocity_ssa

    subroutine calc_visc_eff_3D(visc_eff,ux,uy,ATT,H_ice,f_ice,zeta_aa,dx,dy,n_glen,eps_0,boundaries)
        ! Calculate 3D effective viscosity following L19, Eq. 2
        ! Use of eps_0 ensures non-zero positive viscosity value everywhere 
        ! Note: viscosity is first calculated on ab-nodes, then 
        ! unstaggered back to aa-nodes. This ensures more stability for 
        ! visc_eff (less likely to blow up for low strain rates). 

        ! Note: this routine is identical to that of velocity_diva, except
        ! the shear strain terms duxdz/duydz are set to zero. Although most
        ! equations for effective viscosity in SSA are given in 2D, the 
        ! 3D rate factor implies a calculation in 3D first, then vertical integration
        ! (the integral could be performed on B=A^(-1/n), but we make use of the available 3D visc_eff
        ! variable and maintain the analogy with the DIVA solver)

        implicit none 
        
        real(wp), intent(OUT) :: visc_eff(:,:,:)        ! aa-nodes
        real(wp), intent(IN)  :: ux(:,:)                ! [m/yr] Vertically averaged horizontal velocity, x-component
        real(wp), intent(IN)  :: uy(:,:)                ! [m/yr] Vertically averaged horizontal velocity, y-component
        real(wp), intent(IN)  :: ATT(:,:,:)             ! aa-nodes
        real(wp), intent(IN)  :: H_ice(:,:) 
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: zeta_aa(:)             ! Vertical axis (sigma-coordinates from 0 to 1)
        real(wp), intent(IN)  :: dx
        real(wp), intent(IN)  :: dy
        real(wp), intent(IN)  :: n_glen   
        real(wp), intent(IN)  :: eps_0                  ! [1/yr] Regularization constant (minimum strain rate, ~1e-6)
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer  :: i, j, k
        integer  :: ip1, jp1, im1, jm1 
        integer  :: ip2, jp2, im2, jm2
        integer  :: nx, ny, nz  
        real(wp) :: inv_4dx, inv_4dy 
        real(wp) :: p1, p2, eps_0_sq  
        real(wp) :: eps_sq                              ! [1/yr^2]
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

        real(wp), parameter :: visc_min = 1e5_wp        ! Just for safety 

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

        ! Set equal weighting for all four ab-nodes
        wt_ab = 0.25_wp 

        ! Calculate visc_eff on quadrature points (ab-nodes),
        ! and then average to get grid-centered value (aa-node)
        do j = 1, ny 
        do i = 1, nx 

            if (f_ice(i,j) .eq. 1.0_wp) then 

                im1 = max(i-1,1)    
                ip1 = min(i+1,nx) 
                jm1 = max(j-1,1) 
                jp1 = min(j+1,ny) 

                ! Get strain rate terms
                call staggerdiff_nodes_acx_ab_ice(dudx_ab,ux,f_ice,i,j,dx) 
                call staggerdiff_nodes_acy_ab_ice(dvdy_ab,uy,f_ice,i,j,dy) 
                call staggerdiffcross_nodes_acx_ab_ice(dudy_ab,ux,f_ice,i,j,dy)
                call staggerdiffcross_nodes_acy_ab_ice(dvdx_ab,uy,f_ice,i,j,dx)
                
                ! Loop over column
                do k = 1, nz 

                    ! No vertical strain for this solver 
                    duxdz_ab = 0.0_wp 
                    duydz_ab = 0.0_wp 

                    ! Calculate the total effective strain rate from L19, Eq. 21 
                    eps_sq_ab = dudx_ab**2 + dvdy_ab**2 + dudx_ab*dvdy_ab + 0.25_wp*(dudy_ab+dvdx_ab)**2 &
                                + 0.25_wp*duxdz_ab**2 + 0.25_wp*duydz_ab**2 + eps_0_sq

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

                    ! Calculate effective viscosity on ab-nodes
                    visc_eff_ab = 0.5_wp*(eps_sq_ab)**(p1) * ATT_ab**(p2)

                    ! Calcualte effective viscosity on aa-nodes
                    visc_eff(i,j,k) = sum(wt_ab*visc_eff_ab)

                end do 

            else 

                do k = 1, nz 
                    visc_eff(i,j,k) = 0.0_wp 
                end do 

            end if 

        end do  
        end do 

        ! Set boundaries 
        call set_boundaries_3D_aa(visc_eff,boundaries)

        ! Final safety check (mainly for grid boundaries)
        ! Ensure non-zero visc value everywhere there is ice. 
        do j = 1, ny 
        do i = 1, nx 

            if (f_ice(i,j) .eq. 1.0_wp) then

                do k = 1, nz 
                    if (visc_eff(i,j,k) .lt. visc_min) visc_eff(i,j,k) = visc_min
                end do 

            end if 

        end do
        end do
        
        return 

    end subroutine calc_visc_eff_3D

    subroutine calc_visc_eff_3D_aa(visc_eff,ux,uy,ATT,H_ice,f_ice,zeta_aa,dx,dy,n_glen,eps_0,boundaries)
        ! Calculate 3D effective viscosity following L19, Eq. 2
        ! Use of eps_0 ensures non-zero positive viscosity value everywhere 
        ! Note: viscosity is first calculated on ab-nodes, then 
        ! unstaggered back to aa-nodes. This ensures more stability for 
        ! visc_eff (less likely to blow up for low strain rates). 

        implicit none 
        
        real(wp), intent(OUT) :: visc_eff(:,:,:)        ! aa-nodes
        real(wp), intent(IN)  :: ux(:,:)                ! [m/yr] Vertically averaged horizontal velocity, x-component
        real(wp), intent(IN)  :: uy(:,:)                ! [m/yr] Vertically averaged horizontal velocity, y-component
        real(wp), intent(IN)  :: ATT(:,:,:)             ! aa-nodes
        real(wp), intent(IN)  :: H_ice(:,:) 
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: zeta_aa(:)             ! Vertical axis (sigma-coordinates from 0 to 1)
        real(wp), intent(IN)  :: dx
        real(wp), intent(IN)  :: dy
        real(wp), intent(IN)  :: n_glen   
        real(wp), intent(IN)  :: eps_0                  ! [1/yr] Regularization constant (minimum strain rate, ~1e-6)
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer  :: i, j, k
        integer  :: ip1, jp1, im1, jm1 
        integer  :: ip2, jp2, im2, jm2
        integer  :: nx, ny, nz  
        real(wp) :: inv_4dx, inv_4dy 
        real(wp) :: p1, p2, eps_0_sq  
        real(wp) :: eps_sq                              ! [1/yr^2]

        real(wp) :: dudx_aa, dvdy_aa 
        real(wp) :: dudy_aa_1, dudy_aa_2, dudy_aa 
        real(wp) :: dvdx_aa_1, dvdx_aa_2, dvdx_aa 
        real(wp) :: duxdz_aa, duydz_aa
        real(wp) :: eps_sq_aa, ATT_aa
        
        real(wp), parameter :: visc_min = 1e5_wp        ! Just for safety 

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

        ! Calculate visc_eff on aa-nodes
        do j = 1, ny 
        do i = 1, nx 

            if (f_ice(i,j) .eq. 1.0_wp) then 

                im1 = max(i-1,1)    
                ip1 = min(i+1,nx) 
                jm1 = max(j-1,1) 
                jp1 = min(j+1,ny) 

                ! Get strain rate terms
                dudx_aa = (ux(i,j)-ux(im1,j))/dx 
                dvdy_aa = (uy(i,j)-uy(i,jm1))/dy 
                
                dudy_aa_1 = (ux(i,jp1)-ux(i,jm1))/(2.0_wp*dy)
                dudy_aa_2 = (ux(im1,jp1)-ux(im1,jm1))/(2.0_wp*dy)
                dudy_aa   = 0.5_wp*(dudy_aa_1+dudy_aa_2)

                dvdx_aa_1 = (uy(ip1,j)-uy(im1,j))/(2.0_wp*dx)
                dvdx_aa_2 = (uy(ip1,jm1)-uy(im1,jm1))/(2.0_wp*dx)
                dvdx_aa   = 0.5_wp*(dvdx_aa_1+dvdx_aa_2)

                ! Loop over column
                do k = 1, nz 

                    ! Vertical shear strain rate terms are zero
                    duxdz_aa = 0.0_wp
                    duydz_aa = 0.0_wp

                    ! Calculate the total effective strain rate from L19, Eq. 21 
                    eps_sq_aa = dudx_aa**2 + dvdy_aa**2 + dudx_aa*dvdy_aa + 0.25_wp*(dudy_aa+dvdx_aa)**2 &
                                + 0.25_wp*duxdz_aa**2 + 0.25_wp*duydz_aa**2 + eps_0_sq

                    ! Get rate factor on central node
                    ATT_aa = ATT(i,j,k)

                    ! Calculate effective viscosity on ab-nodes
                    visc_eff(i,j,k) = 0.5_wp*(eps_sq_aa)**(p1) * ATT_aa**(p2)

                end do 

            else 

                do k = 1, nz 
                    visc_eff(i,j,k) = 0.0_wp 
                end do 

            end if 

        end do  
        end do 

        ! Set boundaries 
        call set_boundaries_3D_aa(visc_eff,boundaries)

        ! Final safety check (mainly for grid boundaries)
        ! Ensure non-zero visc value everywhere there is ice. 
        do j = 1, ny 
        do i = 1, nx 

            if (f_ice(i,j) .eq. 1.0_wp) then

                do k = 1, nz 
                    if (visc_eff(i,j,k) .lt. visc_min) visc_eff(i,j,k) = visc_min
                end do 

            end if 

        end do
        end do

        return 

    end subroutine calc_visc_eff_3D_aa
    
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

        ! Apply boundary conditions as needed
        call set_boundaries_2D_aa(visc_eff_int,boundaries) 

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

        real(prec), intent(OUT) :: taub_acx(:,:)   ! [Pa] Basal stress (acx nodes)
        real(prec), intent(OUT) :: taub_acy(:,:)   ! [Pa] Basal stress (acy nodes)
        real(prec), intent(IN)  :: beta_acx(:,:)   ! [Pa a m-1] Basal friction (acx nodes)
        real(prec), intent(IN)  :: beta_acy(:,:)   ! [Pa a m-1] Basal friction (acy nodes)
        real(prec), intent(IN)  :: ux_b(:,:)       ! [m a-1] Basal velocity (acx nodes)
        real(prec), intent(IN)  :: uy_b(:,:)       ! [m a-1] Basal velocity (acy nodes)
        
        real(prec), parameter :: tol = 1e-5_prec 

        ! Calculate basal stress 
        taub_acx = beta_acx * ux_b 
        taub_acy = beta_acy * uy_b 

        ! Avoid underflows
        where(abs(taub_acx) .lt. tol) taub_acx = 0.0_prec 
        where(abs(taub_acy) .lt. tol) taub_acy = 0.0_prec 
        
        return 

    end subroutine calc_basal_stress
    
end module velocity_ssa
