module velocity_hybrid 
    ! This module solves for the ice sheet's horizontal velocity field (ux,uy)
    ! using the 'hybrid' (SIA+SSA) approximation. 

    use yelmo_defs ,only  : sp, dp, wp, prec, tol_underflow, rho_ice, rho_sw, rho_w, g
    use yelmo_tools, only : stagger_aa_ab, stagger_aa_ab_ice, integrate_trapezoid1D_pt, calc_vertical_integrated_2D
    
    use basal_dragging 
    use solver_ssa_sico5 
    use velocity_general, only : set_inactive_margins
    
    implicit none 

    type hybrid_param_class

        character(len=256) :: ssa_lis_opt 
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
    public :: hybrid_param_class 
    public :: calc_velocity_hybrid
      
contains 

    subroutine calc_velocity_hybrid(ux_b,uy_b,taub_acx,taub_acy,visc_eff,visc_eff_int,ssa_mask_acx,ssa_mask_acy, &
                                  ssa_err_acx,ssa_err_acy,ssa_iter_now,beta,beta_acx,beta_acy,c_bed,taud_acx,taud_acy,H_ice, &
                                  f_ice,H_grnd,f_grnd,f_grnd_acx,f_grnd_acy,ATT,zeta_aa,z_sl,z_bed,dx,dy,n_glen,par)
        ! This subroutine is used to solve the horizontal velocity system (ux,uy)
        ! following the hybrid method (combination of SIA and SSA solutions)
        ! Bueler and Brown (2009), Winkelmann et al (2011), Pollard and DeConto (2012)

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
        real(prec), intent(IN)    :: dx                 ! [m]
        real(prec), intent(IN)    :: dy                 ! [m]
        real(prec), intent(IN)    :: n_glen 
        type(hybrid_param_class), intent(IN) :: par       ! List of parameters that should be defined

        ! Local variables 
        integer :: i, j, k, nx, ny, nz_aa, nz_ac, iter
        logical :: is_converged

        real(prec), allocatable :: ux_b_nm1(:,:) 
        real(prec), allocatable :: uy_b_nm1(:,:)    
        integer,    allocatable :: ssa_mask_acx_ref(:,:)
        integer,    allocatable :: ssa_mask_acy_ref(:,:)

        real(prec) :: L2_norm 

        nx    = size(ux_b,1)
        ny    = size(ux_b,2)
        nz_aa = size(zeta_aa,1)

        ! Prepare local variables 
        allocate(ux_b_nm1(nx,ny))
        allocate(uy_b_nm1(nx,ny))
        allocate(ssa_mask_acx_ref(nx,ny))
        allocate(ssa_mask_acy_ref(nx,ny))

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

            ! Calculate 3D effective viscosity
            select case(par%visc_method)

                case(0)
                    ! Impose constant viscosity value 

                    visc_eff = par%visc_const 

                case(1) 
                    ! Calculate 3D effective viscosity, using velocity solution from previous iteration
                    
                    call calc_visc_eff_3D(visc_eff,ux_b,uy_b,ATT,H_ice,f_ice,zeta_aa,dx,dy,n_glen,par%eps_0)

                case DEFAULT 

                    write(*,*) "calc_velocity_hybrid:: Error: visc_method not recognized."
                    write(*,*) "visc_method = ", par%visc_method 
                    stop 

            end select
                    
            ! Calculate depth-integrated effective viscosity
            ! Note L19 uses eta_bar*H in the ssa equation. Yelmo uses eta_int=eta_bar*H directly.
            call calc_visc_eff_int(visc_eff_int,visc_eff,H_ice,f_ice,zeta_aa)
            
            ! Calculate beta (at the ice base)
            call calc_beta(beta,c_bed,ux_b,uy_b,H_ice,f_ice,H_grnd,f_grnd,z_bed,z_sl,par%beta_method, &
                                par%beta_const,par%beta_q,par%beta_u0,par%beta_gl_scale,par%beta_gl_f, &
                                par%H_grnd_lim,par%beta_min,par%boundaries)

            ! Stagger beta
            call stagger_beta(beta_acx,beta_acy,beta,H_ice,f_ice,ux_b,uy_b, &
                        f_grnd,f_grnd_acx,f_grnd_acy,par%beta_gl_stag,par%boundaries)


            ! =========================================================================================
            ! Step 2: Call the SSA solver to obtain new estimate of ux_b/uy_b

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
             

            ! Exit iterations if ssa solution has converged
            if (is_converged) exit 
            
        end do 

        ! Iterations are finished, finalize calculations


        ! Diagnose basal stress 
        call calc_basal_stress(taub_acx,taub_acy,beta_acx,beta_acy,ux_b,uy_b)
        
        return 

    end subroutine calc_velocity_hybrid

    subroutine calc_visc_eff_3D(visc_eff,ux,uy,ATT,H_ice,f_ice,zeta_aa,dx,dy,n_glen,eps_0)
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
        
        real(wp), intent(OUT) :: visc_eff(:,:,:)      ! aa-nodes
        real(wp), intent(IN)  :: ux(:,:)              ! [m/a] Vertically averaged horizontal velocity, x-component
        real(wp), intent(IN)  :: uy(:,:)              ! [m/a] Vertically averaged horizontal velocity, y-component
        real(wp), intent(IN)  :: ATT(:,:,:)           ! aa-nodes
        real(wp), intent(IN)  :: H_ice(:,:) 
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: zeta_aa(:)           ! Vertical axis (sigma-coordinates from 0 to 1)
        real(wp), intent(IN)  :: dx
        real(wp), intent(IN)  :: dy
        real(wp), intent(IN)  :: n_glen   
        real(wp), intent(IN)  :: eps_0                ! [1/a] Regularization constant (minimum strain rate, ~1e-8)
        
        ! Local variables 
        integer    :: i, j, k, nx, ny, nz
        integer    :: ip1, jp1, im1, jm1   
        real(wp) :: inv_4dx, inv_4dy 
        real(wp) :: dudx_ab, dvdy_ab
        real(wp) :: dudy, dvdx
        real(wp) :: duxdz_ab, duydz_ab  
        real(wp) :: p1, p2, eps_0_sq  
        real(wp) :: eps_sq                            ! [1/a^2]
        real(wp) :: ATT_ab
        real(wp) :: visc_eff_ab_now(4)
        real(wp) :: wt_ab(4)
        real(wp) :: wt  

        real(wp), allocatable :: visc_eff_ab(:,:,:)

        nx = size(visc_eff,1)
        ny = size(visc_eff,2)
        nz = size(visc_eff,3)
        
        ! Allocate local arrays 
        allocate(visc_eff_ab(nx,ny,nz)) 

        ! Calculate scaling factors
        inv_4dx = 1.0_prec / (4.0_prec*dx) 
        inv_4dy = 1.0_prec / (4.0_prec*dy) 

        ! Calculate exponents 
        p1 = (1.0_prec - n_glen)/(2.0_prec*n_glen)
        p2 = -1.0_prec/n_glen

        ! Calculate squared minimum strain rate 
        eps_0_sq = eps_0*eps_0 

        do j = 1, ny 
        do i = 1, nx 

            im1 = max(i-1,1) 
            ip1 = min(i+1,nx) 
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny) 

            ! Calculate effective strain components from horizontal stretching on ab-nodes
            dudx_ab = ( (ux(ip1,j) - ux(im1,j)) + (ux(ip1,jp1) - ux(im1,jp1)) ) *inv_4dx
            dvdy_ab = ( (uy(i,jp1) - uy(i,jm1)) + (uy(ip1,jp1) - uy(ip1,jm1)) ) *inv_4dy 

            ! Calculate of cross terms on ab-nodes
            dudy = (ux(i,jp1) - ux(i,j)) / dx 
            dvdx = (uy(ip1,j) - uy(i,j)) / dy 

            ! Loop over column
            do k = 1, nz 

                ! No vertical shear terms here 
                duxdz_ab = 0.0_wp 
                duydz_ab = 0.0_wp 

                ! Calculate the total effective strain rate from L19, Eq. 21 
                eps_sq = dudx_ab**2 + dvdy_ab**2 + dudx_ab*dvdy_ab + 0.25_prec*(dudy+dvdx)**2 &
                       + 0.25_prec*duxdz_ab**2 + 0.25_prec*duydz_ab**2 + eps_0_sq
                
                ATT_ab = 0.25_prec*(ATT(i,j,k)+ATT(ip1,j,k)+ATT(i,jp1,k)+ATT(ip1,jp1,k)) 
                
                ! Calculate effective viscosity on ab-nodes
                visc_eff_ab(i,j,k) = 0.5_prec*(eps_sq)**(p1) * ATT_ab**(p2)

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

                do k = 1, nz

                    visc_eff_ab_now(1) = visc_eff_ab(i,j,k) 
                    visc_eff_ab_now(2) = visc_eff_ab(im1,j,k) 
                    visc_eff_ab_now(3) = visc_eff_ab(im1,jm1,k) 
                    visc_eff_ab_now(4) = visc_eff_ab(i,jm1,k) 

                    ! Calcualte effective viscosity on aa-nodes
                    visc_eff(i,j,k) = sum(wt_ab*visc_eff_ab_now)

                end do 
                
            else 
                ! Just get simple average for ice free points 

                ! Loop over column
                do k = 1, nz 
                    visc_eff(i,j,k) = 0.25_prec*(visc_eff_ab(i,j,k)+visc_eff_ab(im1,j,k) &
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

        ! Treat the corners to avoid extremes
        visc_eff(1,1,:)   = 0.5*(visc_eff(2,1,:)+visc_eff(1,2,:))
        visc_eff(1,ny,:)  = 0.5*(visc_eff(2,ny,:)+visc_eff(1,ny-1,:))
        visc_eff(nx,1,:)  = 0.5*(visc_eff(nx,2,:)+visc_eff(nx-1,1,:))
        visc_eff(nx,ny,:) = 0.5*(visc_eff(nx-1,ny,:)+visc_eff(nx,ny-1,:))

        return 

    end subroutine calc_visc_eff_3D

    subroutine calc_visc_eff_3D_new(visc_eff,ux,uy,ATT,H_ice,f_ice,zeta_aa,dx,dy,n_glen,eps_0)
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
        
        real(wp), intent(OUT) :: visc_eff(:,:,:)      ! aa-nodes
        real(wp), intent(IN)  :: ux(:,:)              ! [m/a] Vertically averaged horizontal velocity, x-component
        real(wp), intent(IN)  :: uy(:,:)              ! [m/a] Vertically averaged horizontal velocity, y-component
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

        !real(wp), allocatable :: ATT_bar(:,:) 

        nx = size(visc_eff,1)
        ny = size(visc_eff,2)
        nz = size(visc_eff,3)
        
        ! Get vertically averaged value 
        !allocate(ATT_bar(nx,ny))
        !ATT_bar = calc_vertical_integrated_2D(ATT,zeta_aa)
            
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

                    ! No vertical shear terms here 
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
if (.FALSE.) then  
                    ! Get the rate factor on ab-nodes too
                    ATT_ab(1) = 0.25_wp*(ATT(i,j,k)+ATT(ip1,j,k)+ATT(i,jp1,k)+ATT(ip1,jp1,k)) 
                    ATT_ab(2) = 0.25_wp*(ATT(i,j,k)+ATT(im1,j,k)+ATT(i,jp1,k)+ATT(im1,jp1,k)) 
                    ATT_ab(3) = 0.25_wp*(ATT(i,j,k)+ATT(im1,j,k)+ATT(i,jm1,k)+ATT(im1,jm1,k)) 
                    ATT_ab(4) = 0.25_wp*(ATT(i,j,k)+ATT(ip1,j,k)+ATT(i,jm1,k)+ATT(ip1,jm1,k)) 
                    ! ATT_ab(1) = 0.25_wp*(ATT_bar(i,j)+ATT_bar(ip1,j)+ATT_bar(i,jp1)+ATT_bar(ip1,jp1)) 
                    ! ATT_ab(2) = 0.25_wp*(ATT_bar(i,j)+ATT_bar(im1,j)+ATT_bar(i,jp1)+ATT_bar(im1,jp1)) 
                    ! ATT_ab(3) = 0.25_wp*(ATT_bar(i,j)+ATT_bar(im1,j)+ATT_bar(i,jm1)+ATT_bar(im1,jm1)) 
                    ! ATT_ab(4) = 0.25_wp*(ATT_bar(i,j)+ATT_bar(ip1,j)+ATT_bar(i,jm1)+ATT_bar(ip1,jm1)) 
else
                    ! Just use the aa-node central value of ATT 
                    ATT_ab = ATT(i,j,k)
                    !ATT_ab = ATT_bar(i,j)
end if

                    ! Calculate effective viscosity on ab-nodes
                    visc_eff_ab = 0.5_wp*(eps_sq_ab)**(p1) * ATT_ab**(p2)

                    ! Calcualte effective viscosity on aa-nodes
                    visc_eff(i,j,k) = sum(wt_ab*visc_eff_ab)

                end do 

            else  
                ! Get simple effective viscosity ignoring staggering

                do k = 1, nz
                    visc_eff(i,j,k) = 0.5_wp*(eps_0_sq)**(p1) * ATT(i,j,k)**(p2)
                    !visc_eff(i,j,k) = 0.5_wp*(eps_0_sq)**(p1) * ATT_bar(i,j)**(p2)
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

    end subroutine calc_visc_eff_3D_new

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
        ! taub = -beta*u 

        implicit none 

        real(prec), intent(OUT) :: taub_acx(:,:)   ! [Pa] Basal stress (acx nodes)
        real(prec), intent(OUT) :: taub_acy(:,:)   ! [Pa] Basal stress (acy nodes)
        real(prec), intent(IN)  :: beta_acx(:,:)   ! [Pa a m-1] Basal friction (acx nodes)
        real(prec), intent(IN)  :: beta_acy(:,:)   ! [Pa a m-1] Basal friction (acy nodes)
        real(prec), intent(IN)  :: ux_b(:,:)       ! [m a-1] Basal velocity (acx nodes)
        real(prec), intent(IN)  :: uy_b(:,:)       ! [m a-1] Basal velocity (acy nodes)
        
        real(prec), parameter :: tol = 1e-3_prec 

        ! Calculate basal stress 
        taub_acx = -beta_acx * ux_b 
        taub_acy = -beta_acy * uy_b 

        ! Avoid underflows
        where(abs(taub_acx) .lt. tol) taub_acx = 0.0_prec 
        where(abs(taub_acy) .lt. tol) taub_acy = 0.0_prec 
        
        return 

    end subroutine calc_basal_stress
    
end module velocity_hybrid
