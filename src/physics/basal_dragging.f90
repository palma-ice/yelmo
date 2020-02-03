module basal_dragging
    ! This module calculate beta, the basal friction coefficient,
    ! needed as input to the SSA solver. It corresponds to
    ! the equation for basal stress:
    !
    ! tau_b_acx = beta_acx * ux
    ! tau_b_acy = beta_acy * uy 
    !
    ! Note that for the stability and correctness of the SSA solver, 
    ! particularly at the grounding line, beta should be defined
    ! directly on the ac nodes (acx,acy). 

    use yelmo_defs, only : sp, dp, prec, pi, g, rho_sw, rho_ice, rho_w  

    use yelmo_tools, only : stagger_aa_acx, stagger_aa_acy

    implicit none 


    private

    ! Effective pressure
    public :: calc_effective_pressure_overburden
    public :: calc_effective_pressure_marine
    public :: calc_effective_pressure_till

    ! c_bed functions
    public :: calc_lambda_bed_lin
    public :: calc_lambda_bed_exp
    public :: calc_lambda_till_const
    public :: calc_lambda_till_linear

    ! beta gl functions 
    public :: scale_beta_gl_fraction 
    public :: scale_beta_gl_Hgrnd
    public :: scale_beta_gl_zstar
        
    ! Beta functions (aa-nodes)
    public :: calc_beta_aa_power_plastic
    public :: calc_beta_aa_reg_coulomb

    ! Beta staggering functions (aa- to ac-nodes)
    public :: stagger_beta_aa_mean
    public :: stagger_beta_aa_upstream
    public :: stagger_beta_aa_downstream    
    public :: stagger_beta_aa_subgrid
    public :: stagger_beta_aa_subgrid_1 
    
contains 

    elemental function calc_effective_pressure_overburden(H_ice,is_float) result(N_eff)
        ! Effective pressure as overburden pressure N_eff = rho*g*H_ice 

        implicit none 

        real(prec), intent(IN) :: H_ice 
        logical,    intent(IN) :: is_float 
        real(prec) :: N_eff                 ! [Pa]

        ! Calculate effective pressure [Pa] (overburden pressure)
        N_eff = (rho_ice*g*H_ice)

        ! No remaining pressure at base for floating ice
        if (is_float) N_eff = 0.0 

        return 

    end function calc_effective_pressure_overburden

    elemental function calc_effective_pressure_marine(H_ice,z_bed,z_sl,H_w,p) result(N_eff)
        ! Effective pressure as a function of connectivity to the ocean
        ! as defined by Leguy et al. (2014), Eq. 14, and modified
        ! by Robinson and Alvarez-Solas to account for basal water pressure (to do!)

        ! Note: input is for a given point, should be on central aa-nodes
        ! or shifted to ac-nodes before entering this routine 

        implicit none 

        real(prec), intent(IN) :: H_ice 
        real(prec), intent(IN) :: z_bed 
        real(prec), intent(IN) :: z_sl 
        real(prec), intent(IN) :: H_w 
        real(prec), intent(IN) :: p       ! [0:1], 0: no ocean connectivity, 1: full ocean connectivity
        real(prec) :: N_eff               ! [Pa]

        ! Local variables  
        real(prec) :: H_float     ! Maximum ice thickness to allow floating ice
        real(prec) :: p_w         ! Pressure of water at the base of the ice sheet
        real(prec) :: x 
        real(prec) :: rho_sw_ice 

        rho_sw_ice = rho_sw/rho_ice 

        ! Determine the maximum ice thickness to allow floating ice
        H_float = max(0.0_prec, rho_sw_ice*(z_sl-z_bed))

        ! Calculate basal water pressure 
        if (H_ice .eq. 0.0) then
            ! No water pressure for ice-free points

            p_w = 0.0 

        else if (H_ice .lt. H_float) then 
            ! Floating ice: water pressure equals ice pressure 

            p_w   = (rho_ice*g*H_ice)

        else
            ! Determine water pressure based on marine connectivity (Leguy et al., 2014, Eq. 14)

            x     = min(1.0_prec, H_float/H_ice)
            p_w   = (rho_ice*g*H_ice)*(1.0_prec - (1.0_prec-x)**p)

        end if 

        ! Calculate effective pressure [Pa] (overburden pressure minus basal water pressure)
        N_eff = (rho_ice*g*H_ice) - p_w 

        return 

    end function calc_effective_pressure_marine

    elemental function calc_effective_pressure_till(H_w,H_ice,f_grnd,H_w_max,N0,delta,e0,Cc) result(N_eff)
        ! Calculate the effective pressure of the till
        ! following van Pelt and Bueler (2015), Eq. 23.
        
        implicit none 
        
        real(prec), intent(IN)    :: H_w
        real(prec), intent(IN)    :: H_ice
        real(prec), intent(IN)    :: f_grnd  
        real(prec), intent(IN)    :: H_w_max            ! [m] Maximum allowed water depth 
        real(prec), intent(IN)    :: N0                 ! [Pa] Reference effective pressure 
        real(prec), intent(IN)    :: delta              ! [--] Fraction of overburden pressure for saturated till
        real(prec), intent(IN)    :: e0                 ! [--] Reference void ratio at N0 
        real(prec), intent(IN)    :: Cc                 ! [--] Till compressibility 
        real(prec)                :: N_eff              ! [Pa] Effective pressure 
        
        ! Local variables 
        real(prec) :: P0, s 
        real(prec) :: N_eff_maxHw

        if (f_grnd .eq. 0.0_prec) then 
            ! No effective pressure at base for floating ice
        
            N_eff = 0.0_prec 

        else 

            ! Get overburden pressure 
            P0 = rho_ice*g*H_ice

            ! Get ratio of water layer thickness to maximum
            s  = min(H_w/H_w_max,1.0)  

            ! Calculate the effective pressure in the till [Pa] (van Pelt and Bueler, 2015, Eq. 23-24)
            N_eff = min( N0*(delta*P0/N0)**s * 10**((e0/Cc)*(1-s)), P0 ) 

            ! Scale by f_grnd 
            N_eff = f_grnd * N_eff 

        end if  

        return 

    end function calc_effective_pressure_till

    subroutine modify_effective_pressure(N_eff,f_pmp) 
        ! Modify effective pressure in the neighborhood 

        implicit none 

        real(prec), intent(INOUT) :: N_eff 
        real(prec), intent(IN)    :: f_pmp 
        
        return 

    end subroutine modify_effective_pressure

    elemental function calc_lambda_bed_lin(z_bed,z0,z1) result(lambda)
        ! Calculate scaling function: linear 
        
        implicit none 
        
        real(prec), intent(IN)    :: z_bed  
        real(prec), intent(IN)    :: z0
        real(prec), intent(IN)    :: z1
        real(prec)                :: lambda 

        lambda = (z_bed - z0) / (z1 - z0)
        
        if (lambda .lt. 0.0) lambda = 0.0 
        if (lambda .gt. 1.0) lambda = 1.0
        
        return 

    end function calc_lambda_bed_lin

    elemental function calc_lambda_bed_exp(z_bed,z0,z1) result(lambda)
        ! Calculate scaling function: exponential 

        implicit none 
        
        real(prec), intent(IN)    :: z_bed  
        real(prec), intent(IN)    :: z0
        real(prec), intent(IN)    :: z1
        real(prec)                :: lambda 

        lambda = exp( (z_bed - z1) / (z1 - z0) )
                
        if (lambda .gt. 1.0) lambda = 1.0
        
        return 

    end function calc_lambda_bed_exp
    
    elemental function calc_lambda_till_const(phi) result(lambda)
        ! Calculate scaling function: till friction angle 
        
        implicit none 
        
        real(prec), intent(IN)    :: phi                ! [deg] Constant till angle  
        real(prec)                :: lambda 

        ! Calculate bed friction coefficient as the tangent
        lambda = tan(phi*pi/180.0)

        return 

    end function calc_lambda_till_const

    elemental function calc_lambda_till_linear(z_bed,z_sl,phi_min,phi_max,phi_zmin,phi_zmax) result(lambda)
        ! Calculate the effective pressure of the till
        ! following van Pelt and Bueler (2015), Eq. 23.
        
        implicit none 
        
        real(prec), intent(IN)    :: z_bed              ! [m] Bedrock elevation
        real(prec), intent(IN)    :: z_sl               ! [m] Sea level 
        real(prec), intent(IN)    :: phi_min            ! [deg] Minimum till angle 
        real(prec), intent(IN)    :: phi_max            ! [deg] Maximum till angle 
        real(prec), intent(IN)    :: phi_zmin           ! [m] C[z=zmin] = tan(phi_min)
        real(prec), intent(IN)    :: phi_zmax           ! [m] C[z=zmax] = tan(phi_max) 
        real(prec)                :: lambda 

        ! Local variables 
        real(prec) :: phi, f_scale, z_rel  

        ! Get bedrock elevation relative to sea level 
!         z_rel = z_sl - z_bed 
        z_rel = z_bed               ! Relative to present-day sea level 

        ! Scale till angle phi as a function bedrock elevation 
        f_scale = (z_rel - phi_zmin) / (phi_zmax - phi_zmin)
        f_scale = min(f_scale,1.0)
        f_scale = max(f_scale,0.0)

        ! Linear ramp between phi_min and phi_max 
        phi = phi_min + f_scale*(phi_max-phi_min)

        ! Calculate bed friction coefficient as the tangent
        lambda = tan(phi*pi/180.0)

        return 

    end function calc_lambda_till_linear

    ! ================================================================================
    !
    ! Beta functions (aa-nodes) 
    !
    ! ================================================================================

    subroutine calc_beta_aa_power_plastic(beta,ux_b,uy_b,C_bed,q,u_0)
        ! Calculate basal friction coefficient (beta) that
        ! enters the SSA solver as a function of basal velocity
        ! using a power-law form following Bueler and van Pelt (2015)
        
        implicit none
        
        real(prec), intent(OUT) :: beta(:,:)        ! aa-nodes
        real(prec), intent(IN)  :: ux_b(:,:)        ! ac-nodes
        real(prec), intent(IN)  :: uy_b(:,:)        ! ac-nodes
        real(prec), intent(IN)  :: C_bed(:,:)       ! aa-nodes
        real(prec), intent(IN)  :: q
        real(prec), intent(IN)  :: u_0              ! [m/a] 

        ! Local variables
        integer    :: i, j, nx, ny
        integer    :: i1, j1 
        real(prec) :: ux_b_mid, uy_b_mid, uxy_b  

        real(prec), parameter :: ub_sq_min = (1e-1_prec)**2  ! [m/a] Minimum velocity is positive small value to avoid divide by zero

        nx = size(beta,1)
        ny = size(beta,2)
        
        ! Initially set friction to zero everywhere
        beta = 0.0_prec 
        
        if (q .eq. 1.0) then 
            ! q==1: linear law, no loops needed 

            beta = C_bed * (1.0_prec / u_0)

        else 
            ! q/=1: Nonlinear law, loops needed 

            ! x-direction 
            do j = 1, ny
            do i = 1, nx

                i1 = max(i-1,1)
                j1 = max(j-1,1) 
                
                ! Calculate magnitude of basal velocity on aa-node 
                ux_b_mid  = 0.5_prec*(ux_b(i1,j)+ux_b(i,j))
                uy_b_mid  = 0.5_prec*(uy_b(i,j1)+uy_b(i,j))
                uxy_b     = (ux_b_mid**2 + uy_b_mid**2 + ub_sq_min)**0.5_prec

                ! Nonlinear beta as a function of basal velocity
                beta(i,j) = C_bed(i,j) * (uxy_b / u_0)**q * (1.0_prec / uxy_b)

            end do
            end do

        end if 
            
        return
        
    end subroutine calc_beta_aa_power_plastic

    subroutine calc_beta_aa_reg_coulomb(beta,ux_b,uy_b,C_bed,q,u_0)
        ! Calculate basal friction coefficient (beta) that
        ! enters the SSA solver as a function of basal velocity
        ! using a regularized Coulomb friction law following
        ! Joughin et al (2019), GRL, Eqs. 2a/2b
        ! Note: Calculated on aa-nodes
        ! Note: beta should be calculated for bed everywhere, 
        ! independent of floatation, which is accounted for later
        
        implicit none
        
        real(prec), intent(OUT) :: beta(:,:)        ! aa-nodes
        real(prec), intent(IN)  :: ux_b(:,:)        ! ac-nodes
        real(prec), intent(IN)  :: uy_b(:,:)        ! ac-nodes
        real(prec), intent(IN)  :: C_bed(:,:)       ! aa-nodes
        real(prec), intent(IN)  :: q
        real(prec), intent(IN)  :: u_0              ! [m/a] 

        ! Local variables
        integer    :: i, j, nx, ny
        integer    :: i1, j1 
        real(prec) :: ux_b_mid, uy_b_mid, uxy_b  

        real(prec), parameter :: ub_sq_min = (1e-1_prec)**2  ! [m/a] Minimum velocity is positive small value to avoid divide by zero

        nx = size(beta,1)
        ny = size(beta,2)
        
        ! Initially set friction to zero everywhere
        beta = 0.0_prec 

        do j = 1, ny
        do i = 1, nx

            i1 = max(i-1,1)
            j1 = max(j-1,1) 
            
            ! Calculate magnitude of basal velocity on aa-node 
            ux_b_mid  = 0.5_prec*(ux_b(i1,j)+ux_b(i,j))
            uy_b_mid  = 0.5_prec*(uy_b(i,j1)+uy_b(i,j))
            uxy_b     = (ux_b_mid**2 + uy_b_mid**2 + ub_sq_min)**0.5_prec

            ! Nonlinear beta as a function of basal velocity
            beta(i,j) = C_bed(i,j) * (uxy_b / (uxy_b+u_0))**q * (1.0_prec / uxy_b)

        end do
        end do 
            
        return
        
    end subroutine calc_beta_aa_reg_coulomb

    ! ================================================================================
    !
    ! Scaling functions 
    !
    ! ================================================================================

    subroutine scale_beta_gl_fraction(beta,f_grnd,f_gl)
        ! Applyt scalar between 0 and 1 to modify basal friction coefficient
        ! at the grounding line.
        
        implicit none
        
        real(prec), intent(INOUT) :: beta(:,:)     ! aa-nodes
        real(prec), intent(IN)    :: f_grnd(:,:)   ! aa-nodes
        real(prec), intent(IN)    :: f_gl          ! Fraction parameter      
        
        ! Local variables
        integer    :: i, j, nx, ny
        integer    :: im1, ip1, jm1, jp1 

        nx = size(f_grnd,1)
        ny = size(f_grnd,2) 

        ! Consistency check 
        if (f_gl .lt. 0.0 .or. f_gl .gt. 1.0) then 
            write(*,*) "scale_beta_gl_fraction:: Error: f_gl must be between 0 and 1."
            write(*,*) "f_gl = ", f_gl
            stop 
        end if 
        
        do j = 1, ny 
        do i = 1, nx-1

            im1 = max(1, i-1)
            ip1 = min(nx,i+1)
            
            jm1 = max(1, j-1)
            jp1 = min(ny,j+1)

            ! Check if point is at the grounding line 
            if (f_grnd(i,j) .gt. 0.0 .and. &
                (f_grnd(im1,j) .eq. 0.0 .or. f_grnd(ip1,j) .eq. 0.0 .or. &
                 f_grnd(i,jm1) .eq. 0.0 .or. f_grnd(i,jp1) .eq. 0.0) ) then 

                ! Set desired grounding-line fraction
                beta(i,j) = beta(i,j) * f_gl 

            end if 

        end do 
        end do  

        return
        
    end subroutine scale_beta_gl_fraction
    
    subroutine scale_beta_gl_Hgrnd(beta,H_grnd,H_grnd_lim)
        ! Calculate scalar between 0 and 1 to modify basal friction coefficient
        ! as ice approaches and achieves floatation, and apply.
        
        implicit none
        
        real(prec), intent(INOUT) :: beta(:,:)     ! aa-nodes
        real(prec), intent(IN)    :: H_grnd(:,:)  ! aa-nodes
        real(prec), intent(IN)    :: H_grnd_lim       

        ! Local variables
        integer    :: i, j, nx, ny
        real(prec) :: f_scale 

        nx = size(H_grnd,1)
        ny = size(H_grnd,2) 

        ! Consistency check 
        if (H_grnd_lim .le. 0.0) then 
            write(*,*) "scale_beta_aa_Hgrnd:: Error: H_grnd_lim must be positive."
            write(*,*) "H_grnd_lim = ", H_grnd_lim
            stop 
        end if 
         
        do j = 1, ny 
        do i = 1, nx

            f_scale = max( min(H_grnd(i,j),H_grnd_lim)/H_grnd_lim, 0.0) 

            beta(i,j) = beta(i,j) * f_scale 

        end do 
        end do  

        return
        
    end subroutine scale_beta_gl_Hgrnd
    
    subroutine scale_beta_gl_zstar(beta,H_ice,z_bed,z_sl,norm)
        ! Calculate scalar between 0 and 1 to modify basal friction coefficient
        ! as ice approaches and achieves floatation, and apply.
        ! Following "Zstar" approach of Gladstone et al. (2017) 
        
        implicit none
        
        real(prec), intent(INOUT) :: beta(:,:)     ! aa-nodes
        real(prec), intent(IN)    :: H_ice(:,:)     ! aa-nodes
        real(prec), intent(IN)    :: z_bed(:,:)     ! aa-nodes        
        real(prec), intent(IN)    :: z_sl(:,:)      ! aa-nodes        
        logical,    intent(IN)    :: norm           ! Normalize by H_ice? 

        ! Local variables
        integer    :: i, j, nx, ny
        real(prec) :: rho_sw_ice 
        real(prec) :: f_scale 

        rho_sw_ice = rho_sw / rho_ice 

        nx = size(H_ice,1)
        ny = size(H_ice,2) 

        ! acx-nodes 
        do j = 1, ny 
        do i = 1, nx

            if (z_bed(i,j) > z_sl(i,j)) then 
                ! Land based above sea level 
                f_scale = H_ice(i,j) 
            else
                ! Marine based 
                f_scale = max(0.0_prec, H_ice(i,j) - (z_sl(i,j)-z_bed(i,j))*rho_sw_ice)
            end if 

            if (norm .and. H_ice(i,j) .gt. 0.0) f_scale = f_scale / H_ice(i,j) 

            beta(i,j) = beta(i,j) * f_scale 
            
        end do 
        end do  

        return
        
    end subroutine scale_beta_gl_zstar

    ! ================================================================================
    !
    ! Staggering functions 
    !
    ! ================================================================================

    subroutine stagger_beta_aa_mean(beta_acx,beta_acy,beta)
        ! Stagger beta from aa-nodes to ac-nodes
        ! using simple staggering method, independent
        ! of any information about flotation, etc. 

        implicit none
        
        real(prec), intent(INOUT) :: beta_acx(:,:)   ! ac-nodes
        real(prec), intent(INOUT) :: beta_acy(:,:)   ! ac-nodes
        real(prec), intent(IN)    :: beta(:,:)       ! aa-nodes

        ! Local variables
        integer    :: i, j, nx, ny

        nx = size(beta_acx,1)
        ny = size(beta_acx,2) 

        ! === Stagger to ac-nodes === 

        ! acx-nodes
        do j = 1, ny 
        do i = 1, nx-1
            beta_acx(i,j) = 0.5_prec*(beta(i,j)+beta(i+1,j))
        end do 
        end do 
        beta_acx(nx,:) = beta_acx(nx-1,:) 
        
        ! acy-nodes 
        do j = 1, ny-1 
        do i = 1, nx
            beta_acy(i,j) = 0.5_prec*(beta(i,j)+beta(i,j+1))
        end do 
        end do 
        beta_acy(:,ny) = beta_acy(:,ny-1) 
        
        return
        
    end subroutine stagger_beta_aa_mean
    
    subroutine stagger_beta_aa_upstream(beta_acx,beta_acy,beta,f_grnd)
        ! Modify basal friction coefficient by grounded/floating binary mask
        ! (via the grounded fraction)
        ! Analagous to method "NSEP" in Seroussi et al (2014): 
        ! Friction is upstream value if a staggered node contains a floating fraction,
        ! ie, f_grnd_acx/acy < 1.0 & > 0.0 

        implicit none
        
        real(prec), intent(INOUT) :: beta_acx(:,:)   ! ac-nodes
        real(prec), intent(INOUT) :: beta_acy(:,:)   ! ac-nodes
        real(prec), intent(IN)    :: beta(:,:)       ! aa-nodes
        real(prec), intent(IN)    :: f_grnd(:,:)     ! aa-nodes    
        
        ! Local variables
        integer    :: i, j, nx, ny 
        logical    :: is_float 

        nx = size(beta_acx,1)
        ny = size(beta_acx,2) 

        ! === Stagger to ac-nodes === 

        ! acx-nodes
        do j = 1, ny 
        do i = 1, nx-1

            if (f_grnd(i,j) .eq. 0.0 .and. f_grnd(i+1,j) .eq. 0.0) then 
                beta_acx(i,j) = 0.0 
            else if (f_grnd(i,j) .gt. 0.0 .and. f_grnd(i+1,j) .eq. 0.0) then 
                beta_acx(i,j) = beta(i,j)
            else if (f_grnd(i,j) .eq. 0.0 .and. f_grnd(i+1,j) .gt. 0.0) then
                beta_acx(i,j) = beta(i+1,j)
            else 
                beta_acx(i,j) = 0.5_prec*(beta(i,j)+beta(i+1,j))
            end if 
            
        end do 
        end do 
        beta_acx(nx,:) = beta_acx(nx-1,:) 
        
        ! acy-nodes 
        do j = 1, ny-1 
        do i = 1, nx

            if (f_grnd(i,j) .eq. 0.0 .and. f_grnd(i,j+1) .eq. 0.0) then 
                beta_acy(i,j) = 0.0 
            else if (f_grnd(i,j) .gt. 0.0 .and. f_grnd(i,j+1) .eq. 0.0) then 
                beta_acy(i,j) = beta(i,j)
            else if (f_grnd(i,j) .eq. 0.0 .and. f_grnd(i,j+1) .gt. 0.0) then
                beta_acy(i,j) = beta(i,j+1)
            else 
                beta_acy(i,j) = 0.5_prec*(beta(i,j)+beta(i,j+1))
            end if 
            
        end do 
        end do 
        beta_acy(:,ny) = beta_acy(:,ny-1) 
        
        return
        
    end subroutine stagger_beta_aa_upstream
    
    subroutine stagger_beta_aa_downstream(beta_acx,beta_acy,beta,f_grnd)
        ! Modify basal friction coefficient by grounded/floating binary mask
        ! (via the grounded fraction)
        ! Analagous to method "NSEP" in Seroussi et al (2014): 
        ! Friction is zero if a staggered node contains a floating fraction,
        ! ie, f_grnd_acx/acy < 1.0 & > 0.0 

        implicit none
        
        real(prec), intent(INOUT) :: beta_acx(:,:)   ! ac-nodes
        real(prec), intent(INOUT) :: beta_acy(:,:)   ! ac-nodes
        real(prec), intent(IN)    :: beta(:,:)       ! aa-nodes
        real(prec), intent(IN)    :: f_grnd(:,:)     ! aa-nodes    
        
        ! Local variables
        integer    :: i, j, nx, ny 
        logical    :: is_float 

        nx = size(beta_acx,1)
        ny = size(beta_acx,2) 

        ! === Stagger to ac-nodes === 

        ! acx-nodes
        do j = 1, ny 
        do i = 1, nx-1

            if (f_grnd(i,j) .eq. 0.0 .or. f_grnd(i+1,j) .eq. 0.0) then 
                beta_acx(i,j) = 0.0 
            else 
                beta_acx(i,j) = 0.5_prec*(beta(i,j)+beta(i+1,j))
            end if 
            
        end do 
        end do 
        beta_acx(nx,:) = beta_acx(nx-1,:) 
        
        ! acy-nodes 
        do j = 1, ny-1 
        do i = 1, nx

            if (f_grnd(i,j) .eq. 0.0 .or. f_grnd(i,j+1) .eq. 0.0) then 
                beta_acy(i,j) = 0.0 
            else 
                beta_acy(i,j) = 0.5_prec*(beta(i,j)+beta(i,j+1))
            end if 
            
        end do 
        end do 
        beta_acy(:,ny) = beta_acy(:,ny-1) 
        
        return
        
    end subroutine stagger_beta_aa_downstream
    
    subroutine stagger_beta_aa_subgrid(beta_acx,beta_acy,beta,f_grnd,f_grnd_acx,f_grnd_acy)
        ! Modify basal friction coefficient by grounded/floating binary mask
        ! (via the grounded fraction)
        ! Analagous to method "NSEP" in Seroussi et al (2014): 
        ! Friction is zero if a staggered node contains a floating fraction,
        ! ie, f_grnd_acx/acy > 1.0 

        implicit none
        
        real(prec), intent(INOUT) :: beta_acx(:,:)      ! ac-nodes
        real(prec), intent(INOUT) :: beta_acy(:,:)      ! ac-nodes
        real(prec), intent(IN)    :: beta(:,:)          ! aa-nodes
        real(prec), intent(IN)    :: f_grnd(:,:)        ! aa-nodes     
        real(prec), intent(IN)    :: f_grnd_acx(:,:)    ! ac-nodes     
        real(prec), intent(IN)    :: f_grnd_acy(:,:)    ! ac-nodes     
        
        ! Local variables
        integer    :: i, j, nx, ny 

        nx = size(beta_acx,1)
        ny = size(beta_acx,2) 

        ! Apply simple staggering to ac-nodes

        ! acx-nodes
        do j = 1, ny 
        do i = 1, nx-1

            if (f_grnd(i,j) .gt. 0.0 .and. f_grnd(i+1,j) .eq. 0.0) then 
                ! Floating to the right 
                beta_acx(i,j) = f_grnd_acx(i,j)*beta(i,j) + (1.0-f_grnd_acx(i,j))*beta(i+1,j)
            else if (f_grnd(i,j) .eq. 0.0 .and. f_grnd(i+1,j) .gt. 0.0) then 
                ! Floating to the left 
                beta_acx(i,j) = (1.0-f_grnd_acx(i,j))*beta(i,j) + f_grnd_acx(i,j)*beta(i+1,j)
            else if (f_grnd(i,j) .gt. 0.0 .and. f_grnd(i+1,j) .gt. 0.0) then 
                ! Fully grounded, simple staggering 
                beta_acx(i,j) = 0.5*(beta(i,j) + beta(i+1,j))
            else 
                ! Fully floating 
                beta_acx(i,j) = 0.0 
            end if 

        end do 
        end do 
        beta_acx(nx,:) = beta_acx(nx-1,:) 
        
        ! acy-nodes 
        do j = 1, ny-1 
        do i = 1, nx

            if (f_grnd(i,j) .gt. 0.0 .and. f_grnd(i,j+1) .eq. 0.0) then 
                ! Floating to the top 
                beta_acy(i,j) = f_grnd_acy(i,j)*beta(i,j) + (1.0-f_grnd_acy(i,j))*beta(i,j+1)
            else if (f_grnd(i,j) .eq. 0.0 .and. f_grnd(i,j+1) .gt. 0.0) then 
                ! Floating to the bottom 
                beta_acy(i,j) = (1.0-f_grnd_acy(i,j))*beta(i,j) + f_grnd_acy(i,j)*beta(i,j+1)
            else if (f_grnd(i,j) .gt. 0.0 .and. f_grnd(i,j+1) .gt. 0.0) then 
                ! Fully grounded, simple staggering 
                beta_acy(i,j) = 0.5*(beta(i,j) + beta(i,j+1))
            else 
                ! Fully floating 
                beta_acy(i,j) = 0.0 
            end if 

        end do 
        end do 
        beta_acy(:,ny) = beta_acy(:,ny-1) 
        
        return
        
    end subroutine stagger_beta_aa_subgrid
    
    subroutine stagger_beta_aa_subgrid_1(beta_acx,beta_acy,beta,H_grnd,f_grnd_acx,f_grnd_acy,dHdt)
        ! Modify basal friction coefficient by grounded/floating binary mask
        ! (via the grounded fraction)
        ! Analagous to method "NSEP" in Seroussi et al (2014): 
        ! Friction is zero if a staggered node contains a floating fraction,
        ! ie, f_grnd_acx/acy > 1.0 

        implicit none
        
        real(prec), intent(INOUT) :: beta_acx(:,:)   ! ac-nodes
        real(prec), intent(INOUT) :: beta_acy(:,:)   ! ac-nodes
        real(prec), intent(INOUT) :: beta(:,:)       ! aa-nodes
        real(prec), intent(IN)    :: H_grnd(:,:)     ! aa-nodes    
        real(prec), intent(IN)    :: f_grnd_acx(:,:) ! ac-nodes 
        real(prec), intent(IN)    :: f_grnd_acy(:,:) ! ac-nodes
        real(prec), intent(IN)    :: dHdt(:,:)       ! aa-nodes  
        
        ! Local variables
        integer    :: i, j, nx, ny
        real(prec) :: H_grnd_now, dHdt_now  
        logical    :: is_float 

        nx = size(beta_acx,1)
        ny = size(beta_acx,2) 

        ! === Treat beta at the grounding line on aa-nodes ===

        ! Cut-off beta according to floating criterion on aa-nodes
        !where (H_grnd .le. 0.0) beta = 0.0 
        !beta = beta*f_grnd

        ! === Stagger to ac-nodes === 

        ! acx-nodes
        do j = 1, ny 
        do i = 1, nx-1

            H_grnd_now = 0.5_prec*(H_grnd(i,j)+H_grnd(i+1,j))
            !dHdt_now   = 0.5_prec*(dHdt(i,j)+dHdt(i+1,j))
            dHdt_now   = max(dHdt(i,j),dHdt(i+1,j))

            if (dHdt_now .gt. 0.0 .and. f_grnd_acx(i,j) .eq. 0.0) then 
                ! Purely floating points only considered floating (ie, domain more grounded)

                is_float = .TRUE. 

            else if (dHdt_now .le. 0.0 .and. f_grnd_acx(i,j) .lt. 1.0) then 
                ! Purely floating and partially grounded points considered floating (ie, domain more floating)
                
                is_float = .TRUE. 

            else 
                ! Grounded 

                is_float = .FALSE. 

            end if 

            if (is_float) then
                ! Consider ac-node floating

                beta_acx(i,j) = 0.0 

            else 
                ! Consider ac-node grounded 

                if (H_grnd(i,j) .gt. 0.0 .and. H_grnd(i+1,j) .le. 0.0) then 
                    beta_acx(i,j) = beta(i,j) !*f_grnd_acx(i,j)
                else if (H_grnd(i,j) .le. 0.0 .and. H_grnd(i+1,j) .gt. 0.0) then
                    beta_acx(i,j) = beta(i+1,j) !*f_grnd_acx(i,j)
                else 
                    beta_acx(i,j) = 0.5_prec*(beta(i,j)+beta(i+1,j))
                end if 
                
            end if 
            
        end do 
        end do 
        beta_acx(nx,:) = beta_acx(nx-1,:) 
        
        ! acy-nodes 
        do j = 1, ny-1 
        do i = 1, nx

!             dHdt_now   = 0.5_prec*(dHdt(i,j)+dHdt(i,j+1))
            dHdt_now   = max(dHdt(i,j),dHdt(i,j+1))
            
            if (dHdt_now .gt. 0.0 .and. f_grnd_acy(i,j) .eq. 0.0) then 
                ! Purely floating points only considered floating (ie, domain more grounded)

                is_float = .TRUE. 

            else if (dHdt_now .le. 0.0 .and. f_grnd_acy(i,j) .lt. 1.0) then 
                ! Purely floating and partially grounded points considered floating (ie, domain more floating)
                
                is_float = .TRUE. 

            else 
                ! Grounded 

                is_float = .FALSE. 

            end if 

            if (is_float) then
                ! Consider ac-node floating

                beta_acy(i,j) = 0.0 

            else 
                ! Consider ac-node grounded 

                if (H_grnd(i,j) .gt. 0.0 .and. H_grnd(i,j+1) .le. 0.0) then 
                    beta_acy(i,j) = beta(i,j) !*f_grnd_acy(i,j)
                else if (H_grnd(i,j) .le. 0.0 .and. H_grnd(i,j+1) .gt. 0.0) then
                    beta_acy(i,j) = beta(i,j+1) !*f_grnd_acy(i,j)
                else 
                    beta_acy(i,j) = 0.5_prec*(beta(i,j)+beta(i,j+1))
                end if 
                
            end if 
            
        end do 
        end do 
        beta_acy(:,ny) = beta_acy(:,ny-1) 
        
        return
        
    end subroutine stagger_beta_aa_subgrid_1
    
    ! ================================================================================
    ! ================================================================================


    function calc_l14_scalar(N_eff,uxy_b,ATT_base,m_drag,m_max,lambda_max) result(f_np)
        ! Calculate a friction scaling coefficient as a function
        ! of velocity and effective pressure, following
        ! Leguy et al. (2014), Eq. 15

        ! Note: input is for a given point, should be on central aa-nodes
        ! or shifted to ac-nodes before entering this routine 

        ! Note: this routine is untested and so far, not used (ajr, 2019-02-01)
        
        implicit none 

        real(prec), intent(IN) :: N_eff 
        real(prec), intent(IN) :: uxy_b 
        real(prec), intent(IN) :: ATT_base 
        real(prec), intent(IN) :: m_drag     
        real(prec), intent(IN) :: m_max 
        real(prec), intent(IN) :: lambda_max 
        real(prec) :: f_np 

        ! Local variables 
        real(prec) :: kappa 

        ! Calculate the velocity scaling 
        kappa = m_max / (lambda_max*ATT_base)

        ! Calculate the scaling coeffcient 
        f_np = (N_eff**m_drag / (kappa*uxy_b + N_eff**m_drag))**1/m_drag 

        if (f_np .lt. 0.0 .or. f_np .gt. 1.0) then 
            write(*,*) "calc_l14_scalar:: f_np out of bounds: f_np = ", f_np 
            stop 
        end if 

        return 

    end function calc_l14_scalar

end module basal_dragging 
