module velocity_diva

    use yelmo_defs ,only  : prec, rho_ice, rho_sw, rho_w, g
    use yelmo_tools, only : stagger_aa_ab, stagger_aa_ab_ice, &
                    calc_vertical_integrated_2D, & 
                    integrate_trapezoid1D_1D, integrate_trapezoid1D_pt, minmax

    use basal_dragging 
    use solver_ssa_sico5 
    
    implicit none 

    private 
    public :: calc_velocity_diva

contains 

    subroutine calc_velocity_diva(ux,uy,ux_i,uy_i,ux_bar,uy_bar,ux_b,uy_b,duxdz,duydz,taub_acx,taub_acy, &
                                  visc_eff,visc_eff_bar,beta,beta_acx,beta_acy,ssa_mask_acx,ssa_mask_acy,c_bed, &
                                  taud_acx,taud_acy,H_ice,H_grnd,f_grnd,f_grnd_acx,f_grnd_acy,ATT,zeta_aa, &
                                  z_sl,z_bed,dx,dy,n_glen,vel_max,beta_gl_stag,boundaries,ssa_solver_opt)
        ! This subroutine is used to solve the horizontal velocity system (ux,uy)
        ! following the Depth-Integrated Viscosity Approximation (DIVA),
        ! as outlined by Lipscomb et al. (2019). Method originally 
        ! proposed by Goldberg (2011), algorithm by Arthern et al (2015), 
        ! updated by Lipscomb et al. (2019).

        implicit none 

        real(prec), intent(INOUT) :: ux(:,:,:)          ! [m/a]
        real(prec), intent(INOUT) :: uy(:,:,:)          ! [m/a]
        real(prec), intent(INOUT) :: ux_i(:,:,:)        ! [m/a]
        real(prec), intent(INOUT) :: uy_i(:,:,:)        ! [m/a]
        real(prec), intent(INOUT) :: ux_bar(:,:)        ! [m/a]
        real(prec), intent(INOUT) :: uy_bar(:,:)        ! [m/a]
        real(prec), intent(INOUT) :: ux_b(:,:)          ! [m/a]
        real(prec), intent(INOUT) :: uy_b(:,:)          ! [m/a]
        real(prec), intent(INOUT) :: duxdz(:,:,:)       ! [1/a]
        real(prec), intent(INOUT) :: duydz(:,:,:)       ! [1/a]
        real(prec), intent(INOUT) :: taub_acx(:,:)      ! [Pa]
        real(prec), intent(INOUT) :: taub_acy(:,:)      ! [Pa]
        real(prec), intent(INOUT) :: visc_eff(:,:,:)    ! [Pa a m]
        real(prec), intent(OUT)   :: visc_eff_bar(:,:)  ! [Pa a m]
        real(prec), intent(OUT)   :: beta(:,:)          ! [Pa a/m]
        real(prec), intent(OUT)   :: beta_acx(:,:)      ! [Pa a/m]
        real(prec), intent(OUT)   :: beta_acy(:,:)      ! [Pa a/m]
        integer,    intent(OUT)   :: ssa_mask_acx(:,:)  ! [-]
        integer,    intent(OUT)   :: ssa_mask_acy(:,:)  ! [-]
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
        real(prec), intent(IN)    :: vel_max            ! [m/a]
        integer,    intent(IN)    :: beta_gl_stag
        character(len=*), intent(IN) :: boundaries      
        character(len=*), intent(IN) :: ssa_solver_opt 

        ! Local variables 
        integer :: nx, ny, nz_aa, nz_ac
        integer :: iter, iter_max   
        real(prec), allocatable :: beta_eff(:,:) 
        real(prec), allocatable :: beta_eff_acx(:,:)
        real(prec), allocatable :: beta_eff_acy(:,:)  
        real(prec), allocatable :: eps_sq(:,:,:) 
        real(prec), allocatable :: F1(:,:) 
        real(prec), allocatable :: F2(:,:) 

        nx    = size(ux,1)
        ny    = size(ux,2)
        nz_aa = size(ux,3)

        iter_max = 2 

        ! Prepare local variables 
        allocate(beta_eff(nx,ny))
        allocate(beta_eff_acx(nx,ny))
        allocate(beta_eff_acy(nx,ny))
        allocate(eps_sq(nx,ny,nz_aa))
        allocate(F1(nx,ny))
        allocate(F2(nx,ny))

        do iter = 1, iter_max 

            ! =========================================================================================
            ! Step 1: Calculate fields needed by ssa solver (visc_eff_bar, beta_eff)

            ! Calculate the 3D vertical shear fields using viscosity estimated from the previous iteration 
            call calc_vertical_shear_3D(duxdz,duydz,taub_acx,taub_acy,visc_eff,zeta_aa)

            ! Calculate the effective strain rate using velocity solution from previous iteration
            call calc_strain_eff_squared(eps_sq,ux_bar,uy_bar,duxdz,duydz,zeta_aa,dx,dy)

            ! Calculate 3D effective viscosity and its vertical average 
            call calc_visc_eff_3D(visc_eff,ATT,eps_sq,n_glen)

            visc_eff_bar = calc_vertical_integrated_2D(visc_eff,zeta_aa) 
            
            ! Calculate F-integerals (F1,F2) on aa-nodes 
            call calc_F_integral(F1,visc_eff,H_ice,zeta_aa,n=1.0_prec)
            call calc_F_integral(F2,visc_eff,H_ice,zeta_aa,n=2.0_prec)

            ! Calculate effective beta 
            call calc_beta_eff(beta_eff,beta,ux_b,uy_b,F2,zeta_aa)

            ! Stagger beta_eff 
            call stagger_beta(beta_eff_acx,beta_eff_acy,beta_eff,f_grnd,f_grnd_acx,f_grnd_acy,beta_gl_stag)

            ! =========================================================================================
            ! Step 2: Call the SSA solver to obtain new estimate of ux_bar/uy_bar
            call calc_vxy_ssa_matrix(ux_bar,uy_bar,beta_eff_acx,beta_eff_acy,visc_eff_bar,ssa_mask_acx,ssa_mask_acy, &
                                H_ice,taud_acx,taud_acy,H_grnd,z_sl,z_bed,dx,dy,vel_max,boundaries,ssa_solver_opt)

            ! Update additional fields based on output of solver
            ! Calculate basal velocity from depth-averaged solution 
            call calc_vel_basal(ux_b,uy_b,ux_bar,uy_bar,F2,beta_acx,beta_acy)
            
            ! Calculate basal stress 
            call calc_basal_stress(taub_acx,taub_acy,beta_acx,beta_acy,ux_b,uy_b)

        end do 

        ! Iterations are finished, finalize calculations of 3D velocity field 

        ! Calculate the 3D horizontal velocity field
        call calc_vel_horizontal_3D(ux,uy,ux_b,uy_b,beta_acx,beta_acy,visc_eff,zeta_aa)

        ! Calculate the surface velocity field 

        return 

    end subroutine calc_velocity_diva 

    subroutine stagger_beta(beta_acx,beta_acy,beta,f_grnd,f_grnd_acx,f_grnd_acy,beta_gl_stag)

        implicit none 

        real(prec), intent(OUT) :: beta_acx(:,:) 
        real(prec), intent(OUT) :: beta_acy(:,:) 
        real(prec), intent(IN)  :: beta(:,:)
        real(prec), intent(IN)  :: f_grnd(:,:) 
        real(prec), intent(IN)  :: f_grnd_acx(:,:) 
        real(prec), intent(IN)  :: f_grnd_acy(:,:) 
        integer,    intent(IN)  :: beta_gl_stag 
        
        ! 5. Apply staggering method with particular care for the grounding line 
        select case(beta_gl_stag) 

            case(0) 
                ! Apply pure staggering everywhere (ac(i) = 0.5*(aa(i)+aa(i+1))
                
                call stagger_beta_aa_mean(beta_acx,beta_acy,beta)

            case(1) 
                ! Apply upstream beta_aa value at ac-node with at least one neighbor H_grnd_aa > 0

                call stagger_beta_aa_upstream(beta_acx,beta_acy,beta,f_grnd)

            case(2) 
                ! Apply downstream beta_aa value (==0.0) at ac-node with at least one neighbor H_grnd_aa > 0

                call stagger_beta_aa_downstream(beta_acx,beta_acy,beta,f_grnd)

            case(3)
                ! Apply subgrid scaling fraction at the grounding line when staggering 

                ! Note: now subgrid treatment is handled on aa-nodes above (using beta_gl_sep)

                call stagger_beta_aa_subgrid(beta_acx,beta_acy,beta,f_grnd, &
                                                f_grnd_acx,f_grnd_acy)

!                 call stagger_beta_aa_subgrid_1(beta_acx,beta_acy,beta,H_grnd, &
!                                             f_grnd,f_grnd_acx,f_grnd_acy)
                
            case DEFAULT 

                write(*,*) "stagger_beta:: Error: beta_gl_stag not recognized."
                write(*,*) "beta_gl_stag = ", beta_gl_stag
                stop 

        end select 

        return 

    end subroutine stagger_beta

    subroutine calc_vel_horizontal_3D(ux,uy,ux_b,uy_b,beta_acx,beta_acy,visc_eff,zeta_aa)
        ! Caluculate the 3D horizontal velocity field (ux,uy)
        ! following L19, Eq. 29 

        implicit none 

        real(prec), intent(OUT) :: ux(:,:,:) 
        real(prec), intent(OUT) :: uy(:,:,:) 
        real(prec), intent(IN)  :: ux_b(:,:) 
        real(prec), intent(IN)  :: uy_b(:,:) 
        real(prec), intent(IN)  :: beta_acx(:,:) 
        real(prec), intent(IN)  :: beta_acy(:,:) 
        real(prec), intent(IN)  :: visc_eff(:,:,:)       
        real(prec), intent(IN)  :: zeta_aa(:) 

        ! Local variables
        integer :: i, j, k, ip1, jp1, nx, ny, nz_aa  
        real(prec) :: tmpval_ac 
        real(prec), allocatable :: visc_eff_ac(:) 
        real(prec), allocatable :: tmpcol_ac(:) 
        
        nx    = size(ux,1)
        ny    = size(ux,2) 
        nz_aa = size(ux,3) 

        allocate(visc_eff_ac(nz_aa))
        allocate(tmpcol_ac(nz_aa))

        do j = 1, ny 
        do i = 1, nx 

            ip1 = min(i+1,nx)
            jp1 = min(j+1,ny) 

            ! === x direction ===============================================

            ! Stagger viscosity column to ac-nodes 
            visc_eff_ac = 0.5_prec*(visc_eff(i,j,:)+visc_eff(ip1,j,:))

            ! Calculate integrated term of L19, Eq. 29 
            tmpcol_ac = integrate_trapezoid1D_1D((1.0_prec/visc_eff_ac)*(1.0-zeta_aa),zeta_aa)

            ! Calculate velocity column 
            ux(i,j,:) = ux_b(i,j) + (beta_acx(i,j)*ux_b(i,j))*tmpcol_ac 

            ! === y direction ===============================================

            ! Stagger viscosity column to ac-nodes 
            visc_eff_ac = 0.5_prec*(visc_eff(i,j,:)+visc_eff(i,jp1,:))

            ! Calculate integrated term of L19, Eq. 29 
            tmpcol_ac = integrate_trapezoid1D_1D((1.0_prec/visc_eff_ac)*(1.0-zeta_aa),zeta_aa)

            ! Calculate velocity column 
            uy(i,j,:) = uy_b(i,j) + (beta_acy(i,j)*uy_b(i,j))*tmpcol_ac 

        end do 
        end do  

        return 

    end subroutine calc_vel_horizontal_3D

    subroutine calc_vertical_shear_3D(duxdz,duydz,taub_acx,taub_acy,visc_eff,zeta_aa)
        ! Calculate vertical shear terms (L19, Eq. 36)

        implicit none 

        real(prec), intent(OUT) :: duxdz(:,:,:)         ! [1/a],    ac-nodes horizontal, aa-nodes vertical 
        real(prec), intent(OUT) :: duydz(:,:,:)         ! [1/a],    ac-nodes horizontal, aa-nodes vertical 
        real(prec), intent(IN)  :: taub_acx(:,:)        ! [Pa],     ac-nodes
        real(prec), intent(IN)  :: taub_acy(:,:)        ! [Pa],     ac-nodes
        real(prec), intent(IN)  :: visc_eff(:,:,:)      ! [Pa a m], aa-nodes
        real(prec), intent(IN)  :: zeta_aa(:)           ! [-]
        
        ! Local variables 
        integer :: i, j, k, nx, ny, nz_aa 
        integer :: ip1, jp1 
        real(prec) :: visc_eff_ac

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
            visc_eff_ac  = 0.5_prec*(visc_eff(i,j,k) + visc_eff(ip1,j,k)) 
            duxdz(i,j,k) = (taub_acx(i,j)/visc_eff_ac) * (1.0-zeta_aa(k))
            
            ! Calculate shear strain, acy-nodes
            visc_eff_ac  = 0.5_prec*(visc_eff(i,j,k) + visc_eff(i,jp1,k)) 
            duydz(i,j,k) = (taub_acy(i,j)/visc_eff_ac) * (1.0-zeta_aa(k))

        end do 
        end do 
        end do 

        return 

    end subroutine calc_vertical_shear_3D 

    subroutine calc_strain_eff_squared(eps_sq,ux,uy,duxdz,duydz,zeta_aa,dx,dy)
        ! Calculate effective strain rate for DIVA solver, using 3D vertical shear 
        ! and 2D depth-averaged horizontal shear (L19, Eq. 21)

        implicit none 
        
        real(prec), intent(OUT) :: eps_sq(:,:,:)        ! [1/a]^2
        real(prec), intent(IN) :: ux(:,:)               ! [m/a] Vertically averaged horizontal velocity, x-component
        real(prec), intent(IN) :: uy(:,:)               ! [m/a] Vertically averaged horizontal velocity, y-component
        real(prec), intent(IN) :: duxdz(:,:,:)          ! [1/a] Vertical shearing, x-component
        real(prec), intent(IN) :: duydz(:,:,:)          ! [1/a] Vertical shearing, x-component
        real(prec), intent(IN) :: zeta_aa(:)            ! Vertical axis (sigma-coordinates from 0 to 1)
        real(prec), intent(IN) :: dx, dy
        
        ! Local variables 
        integer :: i, j, nx, ny, k, nz_aa 
        integer :: im1, ip1, jm1, jp1  
        real(prec) :: inv_4dx, inv_4dy
        real(prec) :: dudx, dudy
        real(prec) :: dvdx, dvdy 
        real(prec) :: duxdz_aa, duydz_aa 

        real(prec), parameter :: epsilon_sq_0 = 1e-6_prec   ! [a^-1] Bueler and Brown (2009), Eq. 26
        
        nx    = size(ux,1)
        ny    = size(ux,2)
        nz_aa = size(zeta_aa,1) 

        inv_4dx = 1.0_prec / (4.0_prec*dx) 
        inv_4dy = 1.0_prec / (4.0_prec*dy) 

        ! Initialize strain rate to zero 
        eps_sq = 0.0_prec 

        ! Loop over domain to calculate viscosity at each aa-node
         
        do j = 1, ny
        do i = 1, nx

            im1 = max(i-1,1) 
            ip1 = min(i+1,nx) 
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny) 
            
            ! Calculate effective strain components from horizontal stretching
            dudx = (ux(i,j) - ux(im1,j))/dx
            dvdy = (uy(i,j) - uy(i,jm1))/dy

            ! Calculate of cross terms on central aa-nodes (symmetrical results)
            dudy = ((ux(i,jp1)   - ux(i,jm1))    &
                  + (ux(im1,jp1) - ux(im1,jm1))) * inv_4dx 
            dvdx = ((uy(ip1,j)   - uy(im1,j))    &
                  + (uy(ip1,jm1) - uy(im1,jm1))) * inv_4dy 

            ! Loop over vertical dimension 
            do k = 1, nz_aa 

                ! Un-stagger shear terms to central aa-nodes in horizontal
                duxdz_aa = 0.5_prec*(duxdz(i,j,k) + duxdz(im1,j,k))
                duydz_aa = 0.5_prec*(duydz(i,j,k) + duydz(i,jm1,k))
                
                ! Calculate the total effective strain rate from L19, Eq. 21 
                eps_sq(i,j,k) = dudx**2 + dvdy**2 + dudx*dvdy + 0.25_prec*(dudy+dvdx)**2 &
                               + 0.25_prec*duxdz_aa**2 + 0.25_prec*duydz_aa**2 &
                               + epsilon_sq_0
            
            end do 

        end do 
        end do 

        return
        
    end subroutine calc_strain_eff_squared

    subroutine calc_visc_eff_3D(visc_eff,ATT,eps_sq,n_glen)
        ! Calculate 3D effective viscosity 
        ! following L19, Eq. 2

        implicit none 
        
        real(prec), intent(OUT) :: visc_eff(:,:,:)  ! aa-nodes
        real(prec), intent(IN)  :: ATT(:,:,:)       ! aa-nodes
        real(prec), intent(IN)  :: eps_sq(:,:,:)    ! aa-nodes
        real(prec), intent(IN)  :: n_glen   

        ! Local variables 
        integer :: i, j, k, nx, ny, nz  
        real(prec) :: mu 

        real(prec), parameter :: visc_min     = 1e3_prec 
        
        nx = size(visc_eff,1)
        ny = size(visc_eff,2)
        nz = size(visc_eff,3)
        
        do k = 1, nz 
        do j = 1, ny 
        do i = 1, nx 

            ! Calculate intermediate term
            mu = 0.5_prec*(eps_sq(i,j,k))**((1.0_prec - n_glen)/(2.0_prec*n_glen))

            ! Calculate effective viscosity 
            visc_eff(i,j,k) = ATT(i,j,k)**(-1.0_prec/n_glen) * mu

            if (visc_eff(i,j,k) .lt. visc_min) visc_eff(i,j,k) = visc_min 

        end do 
        end do  
        end do 

        return 

    end subroutine calc_visc_eff_3D 

    subroutine calc_F_integral(Fint,visc,H_ice,zeta_aa,n)
        ! Useful integrals, following Arthern et al. (2015) Eq. 7,
        ! and Lipscomb et al. (2019), Eq. 30
        ! F_n = int_zb_zs{ 1/visc * ((s-z)/H)**n dz}
        implicit none 

        real(prec), intent(OUT) :: Fint(:,:) 
        real(prec), intent(IN)  :: visc(:,:,:)
        real(prec), intent(IN)  :: H_ice(:,:)
        real(prec), intent(IN)  :: zeta_aa(:)
        real(prec), intent(IN)  :: n  

        ! Local variables 
        integer :: i, j, nx, ny
        real(prec) :: Fint_min 
        real(prec), parameter :: visc_min     = 1e3_prec

        nx = size(visc,1)
        ny = size(visc,2) 

        ! Determine the minimum value of Fint, to assign when H_ice == 0,
        ! since Fint should be nonzero everywhere for numerics
        Fint_min = integrate_trapezoid1D_pt((1.0_prec/visc_min)*(1.0_prec-zeta_aa)**n,zeta_aa)

        ! Vertically integrate at each point
        do j = 1, ny 
        do i = 1, nx 
            if (H_ice(i,j) .gt. 0.0_prec) then 
                ! Viscosity should be nonzero here, perform integration 

                Fint(i,j) = integrate_trapezoid1D_pt((1.0_prec/visc(i,j,:))*(1.0_prec-zeta_aa)**n,zeta_aa)

            else 

                Fint(i,j) = Fint_min

            end if 

        end do 
        end do 

        return

    end subroutine calc_F_integral
    
    subroutine calc_beta_eff(beta_eff,beta,ux_b,uy_b,F2,zeta_aa)
        ! Calculate the depth-averaged horizontal velocity (ux_bar,uy_bar)

        ! Note: L19 staggers the F-integral F2, then solves for beta 

        implicit none 
        
        real(prec), intent(OUT) :: beta_eff(:,:)    ! aa-nodes
        real(prec), intent(IN)  :: beta(:,:)        ! aa-nodes
        real(prec), intent(IN)  :: ux_b(:,:)        ! ac-nodes
        real(prec), intent(IN)  :: uy_b(:,:)        ! ac-nodes
        real(prec), intent(IN)  :: F2(:,:)          ! aa-nodes
        real(prec), intent(IN)  :: zeta_aa(:)       ! aa-nodes

        ! Local variables 
        integer    :: i, j, nx, ny
        integer    :: im1, jm1  
        real(prec) :: uxy_b 

        nx = size(beta_eff,1)
        ny = size(beta_eff,2)

        do j = 1, ny 
        do i = 1, nx 

            im1 = max(i-1,1)
            jm1 = max(j-1,1)

            ! Calculatae basal velocity magnitude at grid center, aa-nodes
            uxy_b = sqrt( 0.5_prec*(ux_b(i,j)+ux_b(im1,j))**2 + 0.5_prec*(ux_b(i,j)+ux_b(i,jm1))**2 )

            if (uxy_b .gt. 0.0) then 
                ! Basal sliding exists, follow L19, Eq. 33

                beta_eff(i,j) = beta(i,j) / (1.0+beta(i,j)*F2(i,j))

            else 
                ! No basal sliding, follow L19, Eq. 35 

                beta_eff(i,j) = 1.0 / F2(i,j) 

            end if 

        end do 
        end do 


        return 

    end subroutine calc_beta_eff 

    subroutine calc_vel_basal(ux_b,uy_b,ux_bar,uy_bar,F2,beta_acx,beta_acy)
        ! Calculate basal sliding following L19, Eq. 32 

        implicit none
        
        real(prec), intent(OUT) :: ux_b(:,:) 
        real(prec), intent(OUT) :: uy_b(:,:)
        real(prec), intent(IN)  :: ux_bar(:,:) 
        real(prec), intent(IN)  :: uy_bar(:,:)
        real(prec), intent(IN)  :: F2(:,:)
        real(prec), intent(IN)  :: beta_acx(:,:) 
        real(prec), intent(IN)  :: beta_acy(:,:)
        
        ! Local variables 
        integer    :: i, j, nx, ny 
        integer    :: ip1, jp1 
        real(prec) :: F2_ac 

        do j = 1, ny 
        do i = 1, nx 

            ip1 = min(i,nx)
            jp1 = min(j,ny)

            F2_ac = 0.5_prec*(F2(i,j) + F2(ip1,j))
            ux_b(i,j) = ux_bar(i,j) / (1.0_prec + beta_acx(i,j)*F2_ac)

            F2_ac = 0.5_prec*(F2(i,j) + F2(i,jp1))
            uy_b(i,j) = uy_bar(i,j) / (1.0_prec + beta_acy(i,j)*F2_ac)

        end do 
        end do  

        return
        
    end subroutine calc_vel_basal

    function calc_vertical_integrated_3D_ice(var,H_ice,sigma) result(var_int)
        ! Vertically integrate a field 3D field (nx,ny,nz)
        ! layer by layer (in the z-direction), return a 3D array
        
        implicit none

        real(prec), intent(IN) :: var(:,:,:)
        real(prec), intent(IN) :: H_ice(:,:) 
        real(prec), intent(IN) :: sigma(:)
        real(prec) :: var_int(size(var,1),size(var,2),size(var,3))

        ! Local variables 
        integer :: i, j, nx, ny

        nx = size(var,1)
        ny = size(var,2)

        do j = 1, ny
        do i = 1, nx
            var_int(i,j,:) = integrate_trapezoid1D_1D(var(i,j,:),sigma*H_ice(i,j))
        end do
        end do

        return

    end function calc_vertical_integrated_3D_ice

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
        
        ! Calculate basal stress 
        taub_acx = -beta_acx * ux_b 
        taub_acy = -beta_acy * uy_b 

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

end module velocity_diva
