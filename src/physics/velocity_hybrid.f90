module velocity_hybrid 
    ! This module solves for the ice sheet's horizontal velocity field (ux,uy)
    ! using the 'hybrid' (SIA+SSA) approximation. 

    use yelmo_defs ,only  : sp, dp, prec, tol_underflow, rho_ice, rho_sw, rho_w, g
    use yelmo_tools, only : stagger_aa_ab, stagger_aa_ab_ice, integrate_trapezoid1D_pt
    
    use basal_dragging 
    use solver_ssa_sico5 

    implicit none 

    type hybrid_param_class

        character(len=256) :: ssa_solver_opt 
        character(len=256) :: boundaries  
        integer    :: beta_method
        real(prec) :: beta_const
        real(prec) :: beta_q                ! Friction law exponent
        real(prec) :: beta_u0               ! [m/a] Friction law velocity threshold 
        integer    :: beta_gl_scale         ! Beta grounding-line scaling method (beta => 0 at gl?)
        integer    :: beta_gl_stag          ! Beta grounding-line staggering method 
        real(prec) :: beta_gl_f             ! Fraction of beta at gl 
        real(prec) :: H_grnd_lim 
        real(prec) :: beta_min              ! Minimum allowed value of beta
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

    subroutine calc_velocity_hybrid(ux_b,uy_b,taub_acx,taub_acy,visc_eff_int,ssa_mask_acx,ssa_mask_acy, &
                                  ssa_err_acx,ssa_err_acy,beta,beta_acx,beta_acy,c_bed,taud_acx,taud_acy,H_ice, &
                                  H_grnd,f_grnd,f_grnd_acx,f_grnd_acy,ATT,zeta_aa,z_sl,z_bed,dx,dy,n_glen,par)
        ! This subroutine is used to solve the horizontal velocity system (ux,uy)
        ! following the hybrid method (combination of SIA and SSA solutions)
        ! Bueler and Brown (2009), Winkelmann et al (2011), Pollard and DeConto (2012)

        implicit none 

        real(prec), intent(INOUT) :: ux_b(:,:)          ! [m/a]
        real(prec), intent(INOUT) :: uy_b(:,:)          ! [m/a]
        real(prec), intent(INOUT) :: taub_acx(:,:)      ! [Pa]
        real(prec), intent(INOUT) :: taub_acy(:,:)      ! [Pa]
        real(prec), intent(OUT)   :: visc_eff_int(:,:)  ! [Pa a m]
        integer,    intent(OUT)   :: ssa_mask_acx(:,:)  ! [-]
        integer,    intent(OUT)   :: ssa_mask_acy(:,:)  ! [-]
        real(prec), intent(OUT)   :: ssa_err_acx(:,:)
        real(prec), intent(OUT)   :: ssa_err_acy(:,:)
        real(prec), intent(OUT)   :: beta(:,:)          ! [Pa a/m]
        real(prec), intent(OUT)   :: beta_acx(:,:)      ! [Pa a/m]
        real(prec), intent(OUT)   :: beta_acy(:,:)      ! [Pa a/m]
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
        
        do iter = 1, par%ssa_iter_max 

            ! Store solution from previous iteration (nm1 == n minus 1) 
            ux_b_nm1 = ux_b 
            uy_b_nm1 = uy_b 
            

            ! =========================================================================================
            ! Step 1: Calculate fields needed by ssa solver (visc_eff_int, beta)

            ! Use classic effective viscosity calculation using only horizontal stretching terms 
            call calc_visc_eff_int(visc_eff_int,ux_b,uy_b,H_ice,ATT,zeta_aa,dx,dy,n_glen)

            ! Calculate beta (at the ice base)
            call calc_beta(beta,c_bed,ux_b,uy_b,H_ice,H_grnd,f_grnd,z_bed,z_sl,par%beta_method, &
                                par%beta_const,par%beta_q,par%beta_u0,par%beta_gl_scale,par%beta_gl_f, &
                                par%H_grnd_lim,par%beta_min,par%boundaries)

            ! Stagger beta
            call stagger_beta(beta_acx,beta_acy,beta,f_grnd,f_grnd_acx,f_grnd_acy,par%beta_gl_stag,par%boundaries)


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
                                     ssa_mask_acx,ssa_mask_acy,H_ice,taud_acx,taud_acy,H_grnd,z_sl, &
                                     z_bed,dx,dy,par%ssa_vel_max,par%boundaries,par%ssa_solver_opt)


            ! Apply relaxation to keep things stable
            call relax_ssa(ux_b,uy_b,ux_b_nm1,uy_b_nm1,rel=par%ssa_iter_rel)
            
            ! Check for convergence
            is_converged = check_vel_convergence_l2rel(ux_b,uy_b,ux_b_nm1,uy_b_nm1,ssa_mask_acx.gt.0,     &
                                                       ssa_mask_acy.gt.0,par%ssa_iter_conv,iter,par%ssa_iter_max, &
                                                       par%ssa_write_log,use_L2_norm=.FALSE.,L2_norm=L2_norm)

            ! Calculate an L1 error metric over matrix for diagnostics
            call check_vel_convergence_l1rel_matrix(ssa_err_acx,ssa_err_acy,ux_b,uy_b,ux_b_nm1,uy_b_nm1)

            
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

    subroutine calc_visc_eff_int(visc,ux,uy,H_ice,ATT,zeta_aa,dx,dy,n_glen)
        ! Calculate effective viscosity eta to be used in SSA solver
        ! Pollard and de Conto (2012), Eqs. 2a/b and Eq. 4 (`visc=mu*H_ice*A**(-1/n)`)
        ! Note: calculated on same nodes as eps_sq (aa-nodes by default)
        ! Note: this is equivalent to the vertically-integrated viscosity, 
        ! since it is multiplied with H_ice 

        implicit none 
        
        real(prec), intent(OUT) :: visc(:,:)    ! [Pa a m]
        real(prec), intent(IN)  :: ux(:,:)      ! Vertically averaged horizontal velocity, x-component
        real(prec), intent(IN)  :: uy(:,:)      ! Vertically averaged horizontal velocity, y-component
        real(prec), intent(IN)  :: H_ice(:,:)   ! Ice thickness
        real(prec), intent(IN)  :: ATT(:,:,:)   ! [a^-1 Pa^-n_glen] (nx,ny,nz_aa) Rate factor
        real(prec), intent(IN)  :: zeta_aa(:)   ! Vertical axis (sigma-coordinates from 0 to 1)
        real(prec), intent(IN)  :: dx, dy
        real(prec), intent(IN)  :: n_glen
        
        ! Local variables 
        integer :: i, j, nx, ny, k, nz_aa 
        integer :: im1, ip1, jm1, jp1  
        real(prec) :: inv_4dx, inv_4dy
        real(prec) :: dudx, dudy
        real(prec) :: dvdx, dvdy 
        real(prec) :: eps_sq, mu  

        real(prec), allocatable :: visc1D(:) 
        
        real(prec), parameter :: visc_min = 1e3_prec    ! [Pa a]
        real(prec), parameter :: eps_sq_0 = 1e-6_prec   ! [a^-1] Bueler and Brown (2009), Eq. 26
        
        nx    = size(ux,1)
        ny    = size(ux,2)
        nz_aa = size(zeta_aa,1) 

        allocate(visc1D(nz_aa))

        inv_4dx = 1.0_prec / (4.0_prec*dx) 
        inv_4dy = 1.0_prec / (4.0_prec*dy) 

        ! Initialize viscosity to minimum value 
        visc = visc_min 
        
        ! Loop over domain to calculate viscosity at each aa-node
         
        do j = 1, ny
        do i = 1, nx

            im1 = max(i-1,1) 
            ip1 = min(i+1,nx) 
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny) 
            
            ! Calculate effective strain components from horizontal stretching on aa-nodes

            dudx = (ux(i,j) - ux(im1,j))/dx
            dvdy = (uy(i,j) - uy(i,jm1))/dy

            ! Calculation of cross terms on central aa-nodes (symmetrical results)
            dudy = ((ux(i,jp1)   - ux(i,jm1))    &
                  + (ux(im1,jp1) - ux(im1,jm1))) * inv_4dx 
            dvdx = ((uy(ip1,j)   - uy(im1,j))    &
                  + (uy(ip1,jm1) - uy(im1,jm1))) * inv_4dy 

            ! Calculate the total effective strain rate due to horizontal stretching terms

            eps_sq = dudx**2 + dvdy**2 + dudx*dvdy + 0.25*(dudy+dvdx)**2 + eps_sq_0
            
            ! Avoid underflow
            if (abs(eps_sq) .lt. tol_underflow) eps_sq = 0.0_prec 
            
            ! 4. Calculate the effective visocity (`eta` in Greve and Blatter, 2009)
            ! Pollard and de Conto (2012), Eqs. 2a/b and Eq. 4 (`visc=A**(-1/n_glen)*mu*H_ice`)

            mu = 0.5_prec*(eps_sq)**((1.0_prec - n_glen)/(2.0_prec*n_glen))

            do k = 1, nz_aa 
                visc1D(k) = ATT(i,j,k)**(-1.0_prec/n_glen) * mu
            end do 

            visc(i,j) = integrate_trapezoid1D_pt(visc1D,zeta_aa) 
            if (H_ice(i,j) .gt. 1.0_prec) visc(i,j) = visc(i,j) * H_ice(i,j) 
            
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
