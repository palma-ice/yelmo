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
    
    use yelmo_defs, only : sp, dp, wp, prec, pi, TOL_UNDERFLOW, io_unit_err, degrees_to_radians

    use yelmo_tools, only : get_neighbor_indices

    use topography, only : calc_H_eff 

    use gaussian_quadrature, only : gq2D_class, gq2D_init, gq2D_to_nodes_aa, &
                                    gq2D_to_nodes_acx, gq2D_to_nodes_acy


    implicit none 


    private

    ! General functions for beta
    public :: calc_cb_ref 
    public :: calc_c_bed 
    public :: calc_beta 
    public :: stagger_beta 

    ! Effective pressure
    public :: calc_effective_pressure_overburden
    public :: calc_effective_pressure_marine
    public :: calc_effective_pressure_till
    public :: calc_effective_pressure_two_value 
    
    ! c_bed functions
    public :: calc_lambda_bed_lin
    public :: calc_lambda_bed_exp

    ! beta gl functions 
    public :: scale_beta_gl_fraction 
    public :: scale_beta_gl_Hgrnd
    public :: scale_beta_gl_zstar
        
    ! Beta functions (aa-nodes)
    public :: calc_beta_aa_power_plastic
    public :: calc_beta_aa_reg_coulomb

    ! Beta staggering functions (aa- to ac-nodes)
    public :: stagger_beta_aa_mean
    public :: stagger_beta_aa_gl_upstream
    public :: stagger_beta_aa_gl_downstream    
    public :: stagger_beta_aa_gl_subgrid
    public :: stagger_beta_aa_gl_subgrid_flux 

contains 
    
    subroutine calc_cb_ref(cb_ref,z_bed,z_bed_sd,z_sl,H_sed,f_sed,H_sed_min,H_sed_max, &
                                            cf_ref,cf_min,z0,z1,n_sd,till_scale,till_method)
        ! Update cb_ref [--] based on parameter choices

        implicit none
        
        real(wp), intent(OUT) :: cb_ref(:,:) 
        real(wp), intent(IN)  :: z_bed(:,:) 
        real(wp), intent(IN)  :: z_bed_sd(:,:)
        real(wp), intent(IN)  :: z_sl(:,:) 
        real(wp), intent(IN)  :: H_sed(:,:) 
        real(wp), intent(IN)  :: f_sed 
        real(wp), intent(IN)  :: H_sed_min
        real(wp), intent(IN)  :: H_sed_max  
        real(wp), intent(IN)  :: cf_ref 
        real(wp), intent(IN)  :: cf_min
        real(wp), intent(IN)  :: z0 
        real(wp), intent(IN)  :: z1
        integer,  intent(IN) :: n_sd 
        character(len=*), intent(IN) :: till_scale
        integer,  intent(IN) :: till_method
        
        ! Local variables
        integer :: q, i, j, nx, ny 
        real(wp) :: f_sd_min, f_sd_max
        real(wp) :: lambda_bed  
        real(wp), allocatable :: cb_ref_samples(:) 
        real(wp), allocatable :: f_sd(:) 
        real(wp), allocatable :: w_sd(:) 

        nx = size(cb_ref,1)
        ny = size(cb_ref,2)
        
        if (n_sd .le. 0) then 
            write(io_unit_err,*) "calc_cb_ref:: Error: ytill.n_sd must be > 0."
            write(io_unit_err,*) "ytill.n_sd = ", n_sd 
            stop 
        end if 

        allocate(cb_ref_samples(n_sd))
        allocate(f_sd(n_sd))
        allocate(w_sd(n_sd))
     
        if (n_sd .eq. 1) then 
            ! No sampling performed
            f_sd = 0.0 
            w_sd = 1.0

        else 

            ! Sample over range, e.g., +/- 1-sigma 
            f_sd_min = -1.0 
            f_sd_max =  1.0 

            do q = 1, n_sd 
                f_sd(q) = f_sd_min + (f_sd_max-f_sd_min)*real(q-1,wp)/real(n_sd-1,wp)
                w_sd(q) = (1.0_wp/sqrt(2.0_wp*pi)) * exp(-f_sd(q)**2/2.0_wp)
            end do 

            ! Normalize weighting 
            w_sd = w_sd/sum(w_sd)

        end if 

        if (till_method .eq. -1) then 
            ! Do nothing - cb_ref defined externally

        else 
            ! Calculate cb_ref following parameter choices 
            ! lambda_bed: scaling as a function of bedrock elevation

            ! Sample bedrock according to its standard deviation to account for uncertainty
            ! and possible pinning points, etc. 

            ! Finally, scale according to sediment parameterization too.

            select case(trim(till_scale))

                case("none")
                    ! No scaling with elevation, set reference value 

                    cb_ref = cf_ref
                    
                case("lin")
                    ! Linear scaling function with bedrock elevation

                    do j = 1, ny 
                    do i = 1, nx 

                        do q = 1, n_sd

                            lambda_bed = calc_lambda_bed_lin(z_bed(i,j)+f_sd(q)*z_bed_sd(i,j), &
                                                                                    z_sl(i,j),z0,z1)

                            ! Calculate cb_ref 
                            cb_ref_samples(q) = cf_ref * lambda_bed 
                            if(cb_ref_samples(q) .lt. cf_min) cb_ref_samples(q) = cf_min 

                        end do 

                        ! Average samples 
                        cb_ref(i,j) = sum(cb_ref_samples*w_sd)

                    end do 
                    end do 

                case("exp") 

                    do j = 1, ny 
                    do i = 1, nx 

                        do q = 1, n_sd

                            lambda_bed = calc_lambda_bed_exp(z_bed(i,j)+f_sd(q)*z_bed_sd(i,j), &
                                                                                    z_sl(i,j),z0,z1)

                            ! Calculate cb_ref 
                            cb_ref_samples(q) = cf_ref * lambda_bed 
                            if(cb_ref_samples(q) .lt. cf_min) cb_ref_samples(q) = cf_min 

                        end do 

                        ! Average samples 
                        cb_ref(i,j) = sum(cb_ref_samples*w_sd)

                    end do 
                    end do 

                case DEFAULT
                    ! Scaling not recognized.

                    write(io_unit_err,*) "calc_cb_ref:: Error: scaling of cb_ref with &
                    &elevation not recognized."
                    write(io_unit_err,*) "ydyn.till_scale = ", till_scale 
                    stop 
                    
            end select 
            
            ! Sediment scaling if desired
            if (f_sed .lt. 1.0) then 

                do j = 1, ny 
                do i = 1, nx 

                    ! Get linear scaling as a function of sediment thickness
                    lambda_bed = (H_sed(i,j) - H_sed_min) / (H_sed_max - H_sed_min)
                    if (lambda_bed .lt. 0.0) lambda_bed = 0.0
                    if (lambda_bed .gt. 1.0) lambda_bed = 1.0 

                    ! Get scaling factor ranging from f_sed(H_sed=H_sed_min) to 1.0(H_sed=H_sed_max)
                    lambda_bed = 1.0 - (1.0-f_sed)*lambda_bed

                    ! Apply scaling to cb_ref 
                    cb_ref(i,j) = cb_ref(i,j) * lambda_bed 

                end do
                end do

            end if 

        end if 

        return 

    end subroutine calc_cb_ref

    subroutine calc_c_bed(c_bed,cb_ref,N_eff,is_angle)

        implicit none 

        real(wp), intent(OUT) :: c_bed(:,:)         ! [Pa]
        real(wp), intent(IN)  :: cb_ref(:,:)        ! [-] or [degrees]
        real(wp), intent(IN)  :: N_eff(:,:)         ! [Pa] 
        logical,  intent(IN) :: is_angle            ! Is cb_ref a till strength angle? 

        if (is_angle) then 
            ! Transform cb_ref by tangent to make 
            ! c_bed = tau_c == critical bed stress
            ! e.g., Bueler and van Pelt (2015); Albrecht et al (2020a);
            ! Zoet and Iverson (2020)

            c_bed = tan(cb_ref*degrees_to_radians)*N_eff 

        else
            ! Treat cb_ref as a normal scalar field

            c_bed = cb_ref*N_eff 

        end if 

        return 

    end subroutine calc_c_bed

    subroutine calc_beta(beta,c_bed,ux_b,uy_b,H_ice,f_ice,H_grnd,f_grnd,z_bed,z_sl,beta_method, &
                         beta_const,beta_q,beta_u0,beta_gl_scale,beta_gl_f,H_grnd_lim, &
                         beta_min,rho_ice,rho_sw,boundaries)

        ! Update beta based on parameter choices

        implicit none
        
        real(wp), intent(INOUT) :: beta(:,:) 
        real(wp), intent(IN)    :: c_bed(:,:)  
        real(wp), intent(IN)    :: ux_b(:,:) 
        real(wp), intent(IN)    :: uy_b(:,:)  
        real(wp), intent(IN)    :: H_ice(:,:) 
        real(wp), intent(IN)    :: f_ice(:,:) 
        real(wp), intent(IN)    :: H_grnd(:,:) 
        real(wp), intent(IN)    :: f_grnd(:,:) 
        real(wp), intent(IN)    :: z_bed(:,:) 
        real(wp), intent(IN)    :: z_sl(:,:) 
        integer,  intent(IN)    :: beta_method
        real(wp), intent(IN)    :: beta_const 
        real(wp), intent(IN)    :: beta_q 
        real(wp), intent(IN)    :: beta_u0 
        integer,  intent(IN)    :: beta_gl_scale  
        real(wp), intent(IN)    :: beta_gl_f
        real(wp), intent(IN)    :: H_grnd_lim
        real(wp), intent(IN)    :: beta_min 
        real(wp), intent(IN)    :: rho_ice 
        real(wp), intent(IN)    :: rho_sw 
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer :: i, j, nx, ny 

        nx = size(beta,1)
        ny = size(beta,2)

        ! 1. Apply beta method of choice 
        select case(beta_method)

            case(-1)
                ! beta (aa-nodes) has been defined externally - do nothing
                
            case(0)
                ! Constant beta everywhere

                beta = beta_const 

            case(1)
                ! Calculate beta from a linear law (simply set beta=c_bed/u0)
                ! (use power-plastic function to ensure proper staggering)

                call calc_beta_aa_power_plastic(beta,ux_b,uy_b,c_bed,f_ice,1.0_wp,beta_u0,boundaries,simple_stagger=.FALSE.)
                
            case(2)
                ! Calculate beta from the quasi-plastic power-law as defined by Bueler and van Pelt (2015)

                call calc_beta_aa_power_plastic(beta,ux_b,uy_b,c_bed,f_ice,beta_q,beta_u0,boundaries,simple_stagger=.FALSE.)
                
            case(3)
                ! Calculate beta from regularized Coulomb law (Joughin et al., GRL, 2019)

                call calc_beta_aa_reg_coulomb(beta,ux_b,uy_b,c_bed,f_ice,beta_q,beta_u0,boundaries,simple_stagger=.FALSE.)
            
            case(4) 
                ! Calculate beta from the quasi-plastic power-law as defined by Bueler and van Pelt (2015)
                ! Use simple-staggering to aa-nodes - useful for Schoof (2006) slab test.

                call calc_beta_aa_power_plastic(beta,ux_b,uy_b,c_bed,f_ice,beta_q,beta_u0,boundaries,simple_stagger=.TRUE.)
            
            case(5)
                ! Calculate beta from regularized Coulomb law (Joughin et al., GRL, 2019)

                call calc_beta_aa_reg_coulomb(beta,ux_b,uy_b,c_bed,f_ice,beta_q,beta_u0,boundaries,simple_stagger=.TRUE.)
            
            case DEFAULT 
                ! Not recognized 

                write(*,*) "calc_beta:: Error: beta_method not recognized."
                write(*,*) "beta_method = ", beta_method
                stop 

        end select 

        ! 1a. Ensure beta is relatively smooth 
!         call regularize2D(beta,H_ice,dx)
!         call limit_gradient(beta,H_ice,dx,log=.TRUE.)

        ! 2. Scale beta as it approaches grounding line 
        select case(beta_gl_scale) 

            case(0) 
                ! Apply fractional parameter at grounding line, no scaling when beta_gl_f=1.0

                call scale_beta_gl_fraction(beta,f_grnd,beta_gl_f)

            case(1) 
                ! Apply H_grnd scaling, reducing beta linearly towards zero at the grounding line 

                call scale_beta_gl_Hgrnd(beta,H_grnd,H_grnd_lim)

            case(2) 
                ! Apply scaling according to thickness above flotation (Zstar approach of Gladstone et al., 2017)
                ! norm==.TRUE., so that zstar-scaling is bounded between 0 and 1, and thus won't affect 
                ! choice of c_bed value that is independent of this scaling. 
                
                call scale_beta_gl_zstar(beta,H_ice,f_ice,z_bed,z_sl,rho_ice,rho_sw,norm=.TRUE.)

            case(3)
                ! Apply f_grnd scaling on aa-nodes 
                ! Note: should be used with simple staggering

                beta = f_grnd*beta 

            case DEFAULT 
                ! No scaling

                write(*,*) "calc_beta:: Error: beta_gl_scale not recognized."
                write(*,*) "beta_gl_scale = ", beta_gl_scale
                stop 

        end select 

        ! 3. Ensure beta==0 for purely floating ice 
        ! Note: assume a binary f_grnd_aa, this does not affect any subgrid gl parameterization
        ! that may be applied during the staggering step.

        !  Simply set beta to zero where purely floating
        where (f_grnd .eq. 0.0) beta = 0.0 
        

        ! Apply additional condition for particular experiments
        select case(trim(boundaries))

            ! ajr: EISMINT case seems to not be necessary anymore
            ! since beta is first calculated on quadrature points
            ! case("EISMINT")
            !     ! Redefine beta at the summit to reduce singularity
            !     ! in symmetric EISMINT experiments with sliding active

            !     i = (nx-1)/2 
            !     j = (ny-1)/2
            !     beta(i,j) = 0.25*(beta(i-1,j)+beta(i+1,j)+beta(i,j-1)+beta(i,j+1))
            
            case("MISMIP3D") 

                ! Redefine beta at the summit to reduce singularity
                ! in MISMIP symmetric experiments
                beta(1,:) = beta(2,:) 

            case("infinite") 
            
                beta(1,:)  = beta(nx-1,:)
                beta(nx,:) = beta(2,:) 
                beta(:,1)  = beta(:,2)
                beta(:,ny) = beta(:,ny-1) 

        end select

        ! Finally ensure that beta for grounded ice is higher than the lower allowed limit
        where(beta .gt. 0.0 .and. beta .lt. beta_min) beta = beta_min 


        ! ================================================================
        ! Note: At this point the beta_aa field is available with beta=0 
        ! for purely floating points and beta > 0 for non-floating points
        ! ================================================================
        

        return 

    end subroutine calc_beta

    subroutine stagger_beta(beta_acx,beta_acy,beta,H_ice,f_ice,ux,uy,f_grnd,f_grnd_acx,f_grnd_acy,beta_gl_stag,beta_min,boundaries)

        implicit none 

        real(wp), intent(INOUT) :: beta_acx(:,:) 
        real(wp), intent(INOUT) :: beta_acy(:,:) 
        real(wp), intent(IN)    :: beta(:,:)
        real(wp), intent(IN)    :: H_ice(:,:)
        real(wp), intent(IN)    :: f_ice(:,:)
        real(wp), intent(IN)    :: ux(:,:)
        real(wp), intent(IN)    :: uy(:,:)
        real(wp), intent(IN)    :: f_grnd(:,:) 
        real(wp), intent(IN)    :: f_grnd_acx(:,:) 
        real(wp), intent(IN)    :: f_grnd_acy(:,:) 
        integer,  intent(IN)    :: beta_gl_stag 
        real(wp), intent(IN)    :: beta_min 
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer :: nx, ny 

        nx = size(beta_acx,1)
        ny = size(beta_acx,2) 

        if (beta_gl_stag .eq. -1) then 

            ! Do nothing - beta_acx/acy has been defined externally

        else  
            ! Calculate staggered beta_acx/acy fields from beta. 

            ! First calculate beta with a standard staggering approach 
            call stagger_beta_aa_mean(beta_acx,beta_acy,beta,f_ice,f_grnd)

            ! Modify beta at the grounding line 
            select case(beta_gl_stag) 
     
                case(-1,0) 

                    ! Do nothing, staggering has already been computed properly 
                    ! -1: beta_acx and beta_acy have been defined externally 
                    !  0: beta_acx and beta_acy have been defined with simple staggering above 

                case(1) 
                    ! Apply upstream beta_aa value at ac-node with at least one neighbor H_grnd_aa > 0

                    call stagger_beta_aa_gl_upstream(beta_acx,beta_acy,beta,f_ice,f_grnd)

                case(2) 
                    ! Apply downstream beta_aa value (==0.0) at ac-node with at least one neighbor H_grnd_aa > 0

                    call stagger_beta_aa_gl_downstream(beta_acx,beta_acy,beta,f_ice,f_grnd)

                case(3)
                    ! Apply subgrid scaling fraction at the grounding line when staggering 

                    call stagger_beta_aa_gl_subgrid(beta_acx,beta_acy,beta,f_ice,f_grnd, &
                                                    f_grnd_acx,f_grnd_acy)

                case(4)
                    ! Apply subgrid scaling fraction at the grounding line when staggering,
                    ! with subgrid weighting calculated from linearly interpolated flux. 

                    call stagger_beta_aa_gl_subgrid_flux(beta_acx,beta_acy,beta, &
                                            H_ice,f_ice,ux,uy,f_grnd,f_grnd_acx,f_grnd_acy)

                case DEFAULT 

                    write(*,*) "stagger_beta:: Error: beta_gl_stag not recognized."
                    write(*,*) "beta_gl_stag = ", beta_gl_stag
                    stop 

            end select 

        end if 

        select case(trim(boundaries))

            case("periodic")

                beta_acx(1,:)    = beta_acx(nx-2,:) 
                beta_acx(nx-1,:) = beta_acx(2,:) 
                beta_acx(nx,:)   = beta_acx(3,:) 
                beta_acx(:,1)    = beta_acx(:,ny-1)
                beta_acx(:,ny)   = beta_acx(:,2) 

                beta_acy(1,:)    = beta_acy(nx-1,:) 
                beta_acy(nx,:)   = beta_acy(2,:) 
                beta_acy(:,1)    = beta_acy(:,ny-2)
                beta_acy(:,ny-1) = beta_acy(:,2) 
                beta_acy(:,ny)   = beta_acy(:,3)

            case("infinite","MISMIP3D") 

                beta_acx(1,:)    = beta_acx(2,:) 
                beta_acx(nx-1,:) = beta_acx(nx-2,:) 
                beta_acx(nx,:)   = beta_acx(nx-2,:) 
                beta_acx(:,1)    = beta_acx(:,2)
                beta_acx(:,ny)   = beta_acx(:,ny-1) 

                beta_acy(1,:)    = beta_acy(2,:) 
                beta_acy(nx,:)   = beta_acy(nx-1,:) 
                beta_acy(:,1)    = beta_acy(:,2)
                beta_acy(:,ny-1) = beta_acy(:,ny-2) 
                beta_acy(:,ny)   = beta_acy(:,ny-2)

        end select 

        ! Finally ensure that beta for grounded ice is higher than the lower allowed limit
        where(beta_acx .gt. 0.0 .and. beta_acx .lt. beta_min) beta_acx = beta_min 
        where(beta_acy .gt. 0.0 .and. beta_acy .lt. beta_min) beta_acy = beta_min 
        
        return 

    end subroutine stagger_beta

    elemental function calc_effective_pressure_overburden(H_ice,f_ice,f_grnd,rho_ice,g) result(N_eff)
        ! Effective pressure as overburden pressure N_eff = rho*g*H_ice 

        implicit none 

        real(wp), intent(IN) :: H_ice 
        real(wp), intent(IN) :: f_ice
        real(wp), intent(IN) :: f_grnd 
        real(wp), intent(IN) :: rho_ice 
        real(wp), intent(IN) :: g
        real(wp) :: N_eff                 ! [Pa]

        ! Local variables 
        real(wp) :: H_eff 

        ! Get effective ice thickness 
        call calc_H_eff(H_eff,H_ice,f_ice,set_frac_zero=.TRUE.)

        ! Calculate effective pressure [Pa] (overburden pressure)
        ! Set N_eff to zero for purely floating points, but do not
        ! scale grounding-line points by f_grnd. This is done on
        ! the beta-staggering step.
        if (f_grnd .gt. 0.0) then
            N_eff = (rho_ice*g*H_eff)
        else
            N_eff = 0.0 
        end if 

        return 

    end function calc_effective_pressure_overburden

    elemental function calc_effective_pressure_marine(H_ice,f_ice,z_bed,z_sl,H_w,p,rho_ice,rho_sw,g) result(N_eff)
        ! Effective pressure as a function of connectivity to the ocean
        ! as defined by Leguy et al. (2014), Eq. 14, and modified
        ! by Robinson and Alvarez-Solas to account for basal water pressure (to do!)

        ! Note: input is for a given point, should be on central aa-nodes
        ! or shifted to ac-nodes before entering this routine 

        implicit none 

        real(wp), intent(IN) :: H_ice 
        real(wp), intent(IN) :: f_ice 
        real(wp), intent(IN) :: z_bed 
        real(wp), intent(IN) :: z_sl 
        real(wp), intent(IN) :: H_w 
        real(wp), intent(IN) :: p       ! [0:1], 0: no ocean connectivity, 1: full ocean connectivity
        real(wp), intent(IN) :: rho_ice 
        real(wp), intent(IN) :: rho_sw
        real(wp), intent(IN) :: g
        
        real(wp) :: N_eff               ! [Pa]

        ! Local variables 
        real(wp) :: H_eff  
        real(wp) :: H_float     ! Maximum ice thickness to allow floating ice
        real(wp) :: p_w         ! Pressure of water at the base of the ice sheet
        real(wp) :: x 
        real(wp) :: rho_sw_ice 

        rho_sw_ice = rho_sw/rho_ice 

        ! Determine the maximum ice thickness to allow floating ice
        H_float = max(0.0_wp, rho_sw_ice*(z_sl-z_bed))

        ! Get effective ice thickness 
        call calc_H_eff(H_eff,H_ice,f_ice,set_frac_zero=.TRUE.)

        ! Calculate basal water pressure 
        if (H_eff .eq. 0.0) then
            ! No water pressure for ice-free points

            p_w = 0.0 

        else if (H_eff .lt. H_float) then 
            ! Floating ice: water pressure equals ice pressure 

            p_w   = (rho_ice*g*H_eff)

        else
            ! Determine water pressure based on marine connectivity (Leguy et al., 2014, Eq. 14)

            x     = min(1.0_wp, H_float/H_eff)
            p_w   = (rho_ice*g*H_eff)*(1.0_wp - (1.0_wp-x)**p)

        end if 

        ! Calculate effective pressure [Pa] (overburden pressure minus basal water pressure)
        ! Note: this will set N_eff to zero for purely floating points, but do not
        ! scale grounding-line points by f_grnd. This is done on
        ! the beta-staggering step.
        N_eff = (rho_ice*g*H_eff) - p_w 

        return 

    end function calc_effective_pressure_marine

    elemental subroutine calc_effective_pressure_till(N_eff,H_w,H_ice,f_ice,f_grnd,H_w_max,N0,delta,e0,Cc,rho_ice,g)
        ! Calculate the effective pressure of the till
        ! following van Pelt and Bueler (2015), Eq. 23.
        
        implicit none 
        
        real(wp), intent(OUT) :: N_eff              ! [Pa] Effective pressure 
        real(wp), intent(IN)  :: H_w
        real(wp), intent(IN)  :: H_ice
        real(wp), intent(IN)  :: f_ice 
        real(wp), intent(IN)  :: f_grnd  
        real(wp), intent(IN)  :: H_w_max            ! [m] Maximum allowed water depth 
        real(wp), intent(IN)  :: N0                 ! [Pa] Reference effective pressure 
        real(wp), intent(IN)  :: delta              ! [--] Fraction of overburden pressure for saturated till
        real(wp), intent(IN)  :: e0                 ! [--] Reference void ratio at N0 
        real(wp), intent(IN)  :: Cc                 ! [--] Till compressibility 
        real(wp), intent(IN) :: rho_ice 
        real(wp), intent(IN) :: g
        
        ! Local variables  
        real(wp) :: H_eff
        real(wp) :: P0, s 
        real(wp) :: q1 

        if (f_grnd .eq. 0.0_wp) then 
            ! No effective pressure at base for floating ice
        
            N_eff = 0.0_wp 

        else 

            ! Get effective ice thickness 
            call calc_H_eff(H_eff,H_ice,f_ice,set_frac_zero=.TRUE.)

            ! Get overburden pressure 
            P0 = rho_ice*g*H_eff

            ! Get ratio of water layer thickness to maximum
            s  = min(H_w/H_w_max,1.0)  

            ! Calculate exponent in expression 
            q1 = (e0/Cc)*(1-s)

            ! Limit exponent to reasonable values to avoid an explosion 
            ! (eg s=1, e0=0.52, Cc=0.014 => q1=37,14)
            q1 = min(q1,10.0_wp) 

            ! Calculate the effective pressure in the till [Pa] (van Pelt and Bueler, 2015, Eq. 23-24)
            ! Note: do not scale grounding-line points by f_grnd. This is done on
            ! the beta-staggering step.
            N_eff = min( N0*(delta*P0/N0)**s * 10**q1, P0 ) 

        end if  

        return 

    end subroutine calc_effective_pressure_till
    
    elemental subroutine calc_effective_pressure_two_value(N_eff,f_pmp,H_ice,f_ice,f_grnd,delta,rho_ice,g)

        implicit none 

        real(wp), intent(OUT) :: N_eff
        real(wp), intent(IN)  :: f_pmp 
        real(wp), intent(IN)  :: H_ice 
        real(wp), intent(IN)  :: f_ice 
        real(wp), intent(IN)  :: f_grnd 
        real(wp), intent(IN)  :: delta
        real(wp), intent(IN) :: rho_ice 
        real(wp), intent(IN) :: g
        
        ! Local variables 
        real(wp) :: H_eff
        real(wp) :: P0, P1

        if (f_grnd .eq. 0.0_wp) then 
            ! No effective pressure at base for purely floating ice
        
            N_eff = 0.0_wp 

        else 

            ! Get effective ice thickness 
            call calc_H_eff(H_eff,H_ice,f_ice,set_frac_zero=.TRUE.)

            ! Get overburden pressure 
            P0 = rho_ice*g*H_eff

            ! Calculate reduced pressure assuming full application of delta
            P1 = P0 * delta 

            ! Calculate effective pressure as a weighted average of 
            ! P0 and P1 via f_pmp 
            N_eff = P0*(1.0_wp-f_pmp) + P1*f_pmp
            
        end if  

        return

    end subroutine calc_effective_pressure_two_value

    elemental function calc_lambda_bed_lin(z_bed,z_sl,z0,z1) result(lambda)
        ! Calculate scaling function: linear 
        
        implicit none 
        
        real(wp), intent(IN)    :: z_bed  
        real(wp), intent(IN)    :: z_sl
        real(wp), intent(IN)    :: z0
        real(wp), intent(IN)    :: z1
        real(wp)                :: lambda 

        ! Local variables 
        real(wp) :: z_rel 
        
        ! Get bedrock elevation relative to sea level 
        !z_rel = z_sl - z_bed        ! Relative to sea-level evolving in time
        z_rel = z_bed               ! Relative to present-day sea level 

        lambda = (z_rel - z0) / (z1 - z0)
        
        if (lambda .lt. 0.0) lambda = 0.0 
        if (lambda .gt. 1.0) lambda = 1.0
        
        return 

    end function calc_lambda_bed_lin

    elemental function calc_lambda_bed_exp(z_bed,z_sl,z0,z1) result(lambda)
        ! Calculate scaling function: exponential 

        implicit none 
        
        real(wp), intent(IN)    :: z_bed  
        real(wp), intent(IN)    :: z_sl
        real(wp), intent(IN)    :: z0
        real(wp), intent(IN)    :: z1
        real(wp)                :: lambda 

        ! Local variables 
        real(wp) :: z_rel 

        ! Get bedrock elevation relative to sea level 
        !z_rel = z_sl - z_bed        ! Relative to sea-level evolving in time
        z_rel = z_bed               ! Relative to present-day sea level 

        lambda = exp( (z_rel - z1) / (z1 - z0) )
                
        if (lambda .gt. 1.0) lambda = 1.0
        
        return 

    end function calc_lambda_bed_exp
    
    ! ================================================================================
    !
    ! Beta functions (aa-nodes) 
    !
    ! ================================================================================

    subroutine calc_beta_aa_power_plastic(beta,ux_b,uy_b,c_bed,f_ice,q,u_0,boundaries,simple_stagger)
        ! Calculate basal friction coefficient (beta) that
        ! enters the SSA solver as a function of basal velocity
        ! using a power-law form following Bueler and van Pelt (2015)
        
        implicit none
        
        real(wp), intent(OUT) :: beta(:,:)        ! aa-nodes
        real(wp), intent(IN)  :: ux_b(:,:)        ! ac-nodes
        real(wp), intent(IN)  :: uy_b(:,:)        ! ac-nodes
        real(wp), intent(IN)  :: c_bed(:,:)       ! aa-nodes
        real(wp), intent(IN)  :: f_ice(:,:)       ! aa-nodes
        real(wp), intent(IN)  :: q
        real(wp), intent(IN)  :: u_0              ! [m/a] 
        character(len=*), intent(IN) :: boundaries 
        logical,  intent(IN)  :: simple_stagger 
        
        ! Local variables
        integer  :: i, j, nx, ny
        integer  :: im1, ip1, jm1, jp1 
        real(wp) :: uxy_b
        real(wp) :: uxn(4) 
        real(wp) :: uyn(4) 
        real(wp) :: uxyn(4) 
        real(wp) :: cbn(4)
        real(wp) :: betan(4)

        real(wp), parameter :: ub_min    = 1e-3_wp          ! [m/yr] Minimum velocity is positive small value to avoid divide by zero
        real(wp), parameter :: ub_sq_min = ub_min**2

        type(gq2D_class) :: gq2D
        real(wp) :: dx_tmp, dy_tmp

        ! Initialize gaussian quadrature calculations
        call gq2D_init(gq2D)
        dx_tmp = 1.0
        dy_tmp = 1.0 

        nx = size(beta,1)
        ny = size(beta,2)
        
        ! Initially set friction to zero everywhere
        beta = 0.0_wp 
        
        do j = 1, ny
        do i = 1, nx

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

            if (f_ice(i,j) .eq. 1.0_wp) then 
                ! Fully ice-covered point

                if (simple_stagger) then 
                    ! Unstagger velocity components to aa-nodes 

                    ! Use central value of c_bed
                    cbn(1:4) = c_bed(i,j) 

                    ! Get velocity components on central node
                    uxn = 0.5*(ux_b(i,j)+ux_b(im1,j))
                    uyn = 0.5*(uy_b(i,j)+uy_b(i,jm1))
                
                else
                    
                    ! Get c_bed on nodes
                    
                    call gq2D_to_nodes_aa(gq2D,cbn,c_bed,dx_tmp,dy_tmp,i,j,im1,ip1,jm1,jp1)
                    !cbn(1:4) = c_bed(i,j)

                    call gq2D_to_nodes_acx(gq2D,uxn,ux_b,dx_tmp,dy_tmp,i,j,im1,ip1,jm1,jp1)
                    call gq2D_to_nodes_acy(gq2D,uyn,uy_b,dx_tmp,dy_tmp,i,j,im1,ip1,jm1,jp1)

                end if 
                
                ! Calculate magnitude of basal velocity on nodes
                uxyn      = sqrt(uxn**2 + uyn**2 + ub_sq_min)

                ! Calculate basal friction
                betan     = c_bed(i,j) * (uxyn / u_0)**q * (1.0_wp / uxyn)
                beta(i,j) = sum(betan*gq2D%wt)/gq2D%wt_tot
            else
                ! Assign minimum velocity value, no staggering for simplicity

                uxy_b     = ub_min
                beta(i,j) = c_bed(i,j) * (uxy_b / u_0)**q * (1.0_wp / uxy_b)

            end if 

        end do
        end do

        return
        
    end subroutine calc_beta_aa_power_plastic

    subroutine calc_beta_aa_reg_coulomb(beta,ux_b,uy_b,c_bed,f_ice,q,u_0,boundaries,simple_stagger)
        ! Calculate basal friction coefficient (beta) that
        ! enters the SSA solver as a function of basal velocity
        ! using a regularized Coulomb friction law following
        ! Joughin et al (2019), GRL, Eqs. 2a/2b
        ! Note: Calculated on aa-nodes
        ! Note: beta should be calculated for bed everywhere, 
        ! independent of floatation, which is accounted for later
        
        implicit none
        
        real(wp), intent(OUT) :: beta(:,:)        ! aa-nodes
        real(wp), intent(IN)  :: ux_b(:,:)        ! ac-nodes
        real(wp), intent(IN)  :: uy_b(:,:)        ! ac-nodes
        real(wp), intent(IN)  :: c_bed(:,:)       ! aa-nodes
        real(wp), intent(IN)  :: f_ice(:,:)       ! aa-nodes
        real(wp), intent(IN)  :: q
        real(wp), intent(IN)  :: u_0              ! [m/a] 
        character(len=*), intent(IN) :: boundaries 
        logical,  intent(IN)  :: simple_stagger 

        ! Local variables
        integer  :: i, j, nx, ny
        integer  :: im1, ip1, jm1, jp1
        integer  :: im1m, ip1m, jm1m, jp1m 
        real(wp) :: uxy_b
        real(wp) :: uxn(4) 
        real(wp) :: uyn(4) 
        real(wp) :: uxyn(4) 
        real(wp) :: cbn(4)
        real(wp) :: betan(4)

        real(wp), parameter :: ub_min    = 1e-3_wp          ! [m/yr] Minimum velocity is positive small value to avoid divide by zero
        real(wp), parameter :: ub_sq_min = ub_min**2

        type(gq2D_class) :: gq2D
        real(wp) :: dx_tmp, dy_tmp

        ! Initialize gaussian quadrature calculations
        call gq2D_init(gq2D)
        dx_tmp = 1.0
        dy_tmp = 1.0 

        nx = size(beta,1)
        ny = size(beta,2)
        
        ! Initially set friction to zero everywhere
        beta = 0.0_wp 

        !!$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,im1m,ip1m,jm1m,jp1m,cbn,uxn,uyn,uxyn,betan,uxy_b)
        do j = 1, ny
        do i = 1, nx

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

            if (f_ice(i,j) .eq. 1.0_wp) then 
                ! Fully ice-covered point

                if (simple_stagger) then 
                    ! Unstagger velocity components to aa-nodes 

                    ! Use central value of c_bed
                    cbn(1:4) = c_bed(i,j) 

                    ! Get velocity components on central node
                    uxn = 0.5*(ux_b(i,j)+ux_b(im1,j))
                    uyn = 0.5*(uy_b(i,j)+uy_b(i,jm1))
                
                else
                    ! Get c_bed on nodes
                    
                    call gq2D_to_nodes_aa(gq2D,cbn,c_bed,dx_tmp,dy_tmp,i,j,im1,ip1,jm1,jp1)
                    !cbn(1:4) = c_bed(i,j) 

                    call gq2D_to_nodes_acx(gq2D,uxn,ux_b,dx_tmp,dy_tmp,i,j,im1,ip1,jm1,jp1)
                    call gq2D_to_nodes_acy(gq2D,uyn,uy_b,dx_tmp,dy_tmp,i,j,im1,ip1,jm1,jp1)

                end if 
                
                ! Calculate magnitude of basal velocity on nodes
                uxyn      = sqrt(uxn**2 + uyn**2 + ub_sq_min)

                ! Calculate basal friction
                betan     = cbn * (uxyn / (uxyn+u_0))**q * (1.0_wp / uxyn)
                beta(i,j) = sum(betan*gq2D%wt)/gq2D%wt_tot

            else
                ! Assign minimum velocity value, ignore staggering for simplicity 

                uxy_b  = ub_min 
                beta(i,j) = c_bed(i,j) * (uxy_b / (uxy_b+u_0))**q * (1.0_wp / uxy_b)

            end if 

        end do
        end do
        !!$omp end parallel do
        
        return
        
    end subroutine calc_beta_aa_reg_coulomb

    ! ================================================================================
    !
    ! Scaling functions 
    !
    ! ================================================================================

    subroutine scale_beta_gl_fraction(beta,f_grnd,f_gl)
        ! Apply scalar between 0 and 1 to modify basal friction coefficient
        ! at the grounding line.
        
        implicit none
        
        real(wp), intent(INOUT) :: beta(:,:)     ! aa-nodes
        real(wp), intent(IN)    :: f_grnd(:,:)   ! aa-nodes
        real(wp), intent(IN)    :: f_gl          ! Fraction parameter      
        
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
       
        !!$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1)
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
        !!$omp end parallel do

        return
        
    end subroutine scale_beta_gl_fraction
    
    subroutine scale_beta_gl_Hgrnd(beta,H_grnd,H_grnd_lim)
        ! Calculate scalar between 0 and 1 to modify basal friction coefficient
        ! as ice approaches and achieves floatation, and apply.
        
        implicit none
        
        real(wp), intent(INOUT) :: beta(:,:)        ! aa-nodes
        real(wp), intent(IN)    :: H_grnd(:,:)      ! aa-nodes
        real(wp), intent(IN)    :: H_grnd_lim       

        ! Local variables
        integer    :: i, j, nx, ny
        real(wp) :: f_scale 

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
    
    subroutine scale_beta_gl_zstar(beta,H_ice,f_ice,z_bed,z_sl,rho_ice,rho_sw,norm)
        ! Calculate scalar between 0 and 1 to modify basal friction coefficient
        ! as ice approaches and achieves floatation, and apply.
        ! Following "Zstar" approach of Gladstone et al. (2017) 
        
        implicit none
        
        real(wp), intent(INOUT) :: beta(:,:)        ! aa-nodes
        real(wp), intent(IN)    :: H_ice(:,:)       ! aa-nodes
        real(wp), intent(IN)    :: f_ice(:,:)       ! aa-nodes
        real(wp), intent(IN)    :: z_bed(:,:)       ! aa-nodes        
        real(wp), intent(IN)    :: z_sl(:,:)        ! aa-nodes        
        real(wp), intent(IN)    :: rho_ice 
        real(wp), intent(IN)    :: rho_sw
        logical,  intent(IN)    :: norm             ! Normalize by H_ice? 
        
        ! Local variables
        integer    :: i, j, nx, ny
        real(wp) :: rho_sw_ice 
        real(wp) :: H_eff
        real(wp) :: f_scale 

        rho_sw_ice = rho_sw / rho_ice 

        nx = size(H_ice,1)
        ny = size(H_ice,2) 

        ! acx-nodes 
        do j = 1, ny 
        do i = 1, nx

            if (f_ice(i,j) .gt. 0.0_wp) then 

                ! Get effective ice thickness 
                call calc_H_eff(H_eff,H_ice(i,j),f_ice(i,j),set_frac_zero=.TRUE.)

                if (z_bed(i,j) > z_sl(i,j)) then 
                    ! Land based above sea level 
                    f_scale = H_eff 
                else
                    ! Marine based 
                    f_scale = max(0.0_wp, H_eff - (z_sl(i,j)-z_bed(i,j))*rho_sw_ice)
                end if 
                
                if (norm .and. H_eff .gt. 0.0) f_scale = f_scale / H_eff

                beta(i,j) = beta(i,j) * f_scale 
            
            end if 

        end do 
        end do  

        return
        
    end subroutine scale_beta_gl_zstar

    ! ================================================================================
    !
    ! Staggering functions 
    !
    ! ================================================================================

    subroutine stagger_beta_aa_mean(beta_acx,beta_acy,beta,f_ice,f_grnd)
        ! Stagger beta from aa-nodes to ac-nodes
        ! using simple staggering method, independent
        ! of any information about flotation, etc. 

        implicit none
        
        real(wp), intent(INOUT) :: beta_acx(:,:)   ! ac-nodes
        real(wp), intent(INOUT) :: beta_acy(:,:)   ! ac-nodes
        real(wp), intent(IN)    :: beta(:,:)       ! aa-nodes
        real(wp), intent(IN)    :: f_ice(:,:)      ! aa-nodes
        real(wp), intent(IN)    :: f_grnd(:,:)     ! aa-nodes

        ! Local variables
        integer :: i, j, nx, ny
        integer :: im1, ip1, jm1, jp1 

        nx = size(beta_acx,1)
        ny = size(beta_acx,2) 

        ! === Stagger to ac-nodes === 
        do j = 1, ny 
        do i = 1, nx

            ! Get neighbor indices 
            im1 = max(i-1,1)
            ip1 = min(i+1,nx)
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny)
            
            if (f_grnd(i,j) .eq. 0.0 .and. f_grnd(ip1,j) .eq. 0.0) then 
                ! Fully floating node 

                beta_acx(i,j) = 0.0 

            else  
                ! Other nodes 

                ! acx-nodes
                if (f_ice(i,j) .eq. 1.0 .and. f_ice(ip1,j) .lt. 1.0) then 
                    beta_acx(i,j) = beta(i,j) 
                else if (f_ice(i,j) .lt. 1.0 .and. f_ice(ip1,j) .eq. 1.0) then
                    beta_acx(i,j) = beta(ip1,j)  
                else 
                    beta_acx(i,j) = 0.5_wp*(beta(i,j)+beta(ip1,j))
                end if 

            end if 

            if (f_grnd(i,j) .eq. 0.0 .and. f_grnd(i,jp1) .eq. 0.0) then 
                ! Fully floating node 

                beta_acy(i,j) = 0.0 

            else  
                ! Other nodes 

                ! acy-nodes
                if (f_ice(i,j) .eq. 1.0 .and. f_ice(i,jp1) .lt. 1.0) then 
                    beta_acy(i,j) = beta(i,j) 
                else if (f_ice(i,j) .lt. 1.0 .and. f_ice(i,jp1) .eq. 1.0) then
                    beta_acy(i,j) = beta(i,jp1)  
                else 
                    beta_acy(i,j) = 0.5_wp*(beta(i,j)+beta(i,jp1))
                end if 
            
            end if 

        end do 
        end do 
        
        return
        
    end subroutine stagger_beta_aa_mean
    
    subroutine stagger_beta_aa_gl_upstream(beta_acx,beta_acy,beta,f_ice,f_grnd)
        ! Modify basal friction coefficient by grounded/floating binary mask
        ! (via the grounded fraction)
        ! Analagous to method "NSEP" in Seroussi et al (2014): 
        ! Friction is upstream value if a staggered node contains a floating fraction,
        ! ie, f_grnd_acx/acy < 1.0 & > 0.0 

        implicit none
        
        real(wp), intent(INOUT) :: beta_acx(:,:)    ! ac-nodes
        real(wp), intent(INOUT) :: beta_acy(:,:)    ! ac-nodes
        real(wp), intent(IN)    :: beta(:,:)        ! aa-nodes
        real(wp), intent(IN)    :: f_ice(:,:)       ! aa-nodes
        real(wp), intent(IN)    :: f_grnd(:,:)      ! aa-nodes    
        
        ! Local variables
        integer    :: i, j, nx, ny
        integer    :: im1, ip1, jm1, jp1  
        logical    :: is_float 

        nx = size(beta_acx,1)
        ny = size(beta_acx,2) 

        ! === Stagger to ac-nodes === 

        do j = 1, ny 
        do i = 1, nx

            ! Get neighbor indices 
            im1 = max(i-1,1)
            ip1 = min(i+1,nx)
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny)
            
            ! grounding line, acx-nodes
            if (f_ice(i,j) .eq. 1.0 .and. f_ice(ip1,j) .eq. 1.0) then 
                ! Both aa-nodes are ice covered 

                if (f_grnd(i,j) .gt. 0.0 .and. f_grnd(ip1,j) .eq. 0.0) then 
                    beta_acx(i,j) = beta(i,j)
                else if (f_grnd(i,j) .eq. 0.0 .and. f_grnd(ip1,j) .gt. 0.0) then
                    beta_acx(i,j) = beta(ip1,j)
                end if 
            
            end if 

            ! grounding line, acy-nodes
            if (f_ice(i,j) .eq. 1.0 .and. f_ice(i,jp1) .eq. 1.0) then 
                ! Both aa-nodes are ice covered 

                if (f_grnd(i,j) .gt. 0.0 .and. f_grnd(i,jp1) .eq. 0.0) then 
                    beta_acy(i,j) = beta(i,j)
                else if (f_grnd(i,j) .eq. 0.0 .and. f_grnd(i,jp1) .gt. 0.0) then
                    beta_acy(i,j) = beta(i,jp1)
                end if 
            
            end if 

        end do 
        end do

        return
        
    end subroutine stagger_beta_aa_gl_upstream
    
    subroutine stagger_beta_aa_gl_downstream(beta_acx,beta_acy,beta,f_ice,f_grnd)
        ! Modify basal friction coefficient by grounded/floating binary mask
        ! (via the grounded fraction)
        ! Analagous to method "NSEP" in Seroussi et al (2014): 
        ! Friction is zero if a staggered node contains a floating fraction,
        ! ie, f_grnd_acx/acy < 1.0 & > 0.0 

        implicit none
        
        real(wp), intent(INOUT) :: beta_acx(:,:)    ! ac-nodes
        real(wp), intent(INOUT) :: beta_acy(:,:)    ! ac-nodes
        real(wp), intent(IN)    :: beta(:,:)        ! aa-nodes
        real(wp), intent(IN)    :: f_ice(:,:)       ! aa-nodes
        real(wp), intent(IN)    :: f_grnd(:,:)      ! aa-nodes    
        
        ! Local variables
        integer :: i, j, nx, ny
        integer :: im1, ip1, jm1, jp1 
        logical :: is_float 

        nx = size(beta_acx,1)
        ny = size(beta_acx,2) 

        ! === Stagger to ac-nodes === 

        ! acx-nodes
        do j = 1, ny 
        do i = 1, nx

            ! Get neighbor indices 
            im1 = max(i-1,1)
            ip1 = min(i+1,nx)
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny)
            
            ! grounding line, acx-nodes
            if (f_ice(i,j) .eq. 1.0 .and. f_ice(ip1,j) .eq. 1.0) then 
                ! Both aa-nodes are ice covered 

                if (f_grnd(i,j) .gt. 0.0 .and. f_grnd(ip1,j) .eq. 0.0) then 
                    beta_acx(i,j) = beta(ip1,j)
                else if (f_grnd(i,j) .eq. 0.0 .and. f_grnd(ip1,j) .gt. 0.0) then
                    beta_acx(i,j) = beta(i,j)
                end if 
            
            end if 

            ! grounding line, acy-nodes
            if (f_ice(i,j) .eq. 1.0 .and. f_ice(i,jp1) .eq. 1.0) then 
                ! Both aa-nodes are ice covered 

                if (f_grnd(i,j) .gt. 0.0 .and. f_grnd(i,jp1) .eq. 0.0) then 
                    beta_acy(i,j) = beta(i,jp1)
                else if (f_grnd(i,j) .eq. 0.0 .and. f_grnd(i,jp1) .gt. 0.0) then
                    beta_acy(i,j) = beta(i,j)
                end if 

            end if
            
        end do 
        end do 

        return
        
    end subroutine stagger_beta_aa_gl_downstream
    
    subroutine stagger_beta_aa_gl_subgrid(beta_acx,beta_acy,beta,f_ice,f_grnd,f_grnd_acx,f_grnd_acy)
        ! Modify basal friction coefficient by grounded/floating binary mask
        ! (via the grounded fraction)
        ! Analagous to method "NSEP" in Seroussi et al (2014): 
        ! Friction is zero if a staggered node contains a floating fraction,
        ! ie, f_grnd_acx/acy > 1.0 

        implicit none
        
        real(wp), intent(INOUT) :: beta_acx(:,:)        ! ac-nodes
        real(wp), intent(INOUT) :: beta_acy(:,:)        ! ac-nodes
        real(wp), intent(IN)    :: beta(:,:)            ! aa-nodes
        real(wp), intent(IN)    :: f_ice(:,:)           ! aa-nodes
        real(wp), intent(IN)    :: f_grnd(:,:)          ! aa-nodes     
        real(wp), intent(IN)    :: f_grnd_acx(:,:)      ! ac-nodes     
        real(wp), intent(IN)    :: f_grnd_acy(:,:)      ! ac-nodes     
        
        ! Local variables
        integer  :: i, j, nx, ny 
        integer  :: im1, ip1, jm1, jp1
        real(wp) :: wt 

        nx = size(beta_acx,1)
        ny = size(beta_acx,2) 

        ! Apply simple staggering to ac-nodes
        do j = 1, ny 
        do i = 1, nx

            ! Get neighbor indices 
            im1 = max(i-1,1)
            ip1 = min(i+1,nx)
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny)
            
            ! grounding line, acx-nodes
            if (f_ice(i,j) .eq. 1.0 .and. f_ice(ip1,j) .eq. 1.0) then 
                ! Both aa-nodes are ice covered 

                ! Get weighting term as a function of grounded fraction 
                !wt = f_grnd_acx(i,j) 
                wt = f_grnd_acx(i,j)**2 

                if (f_grnd(i,j) .gt. 0.0 .and. f_grnd(ip1,j) .eq. 0.0) then 
                    ! Floating to the right 
                    beta_acx(i,j) = wt*beta(i,j) + (1.0-wt)*beta(ip1,j)
                else if (f_grnd(i,j) .eq. 0.0 .and. f_grnd(ip1,j) .gt. 0.0) then 
                    ! Floating to the left 
                    beta_acx(i,j) = (1.0-wt)*beta(i,j) + wt*beta(ip1,j)
                end if 

            end if 

            ! grounding line, acy-nodes 
            if (f_ice(i,j) .eq. 1.0 .and. f_ice(i,jp1) .eq. 1.0) then 
                ! Both aa-nodes are ice covered 

                ! Get weighting term as a function of grounded fraction 
                !wt = f_grnd_acy(i,j) 
                wt = f_grnd_acy(i,j)**2 

                if (f_grnd(i,j) .gt. 0.0 .and. f_grnd(i,jp1) .eq. 0.0) then 
                    ! Floating to the top 
                    beta_acy(i,j) = wt*beta(i,j) + (1.0-wt)*beta(i,jp1)
                else if (f_grnd(i,j) .eq. 0.0 .and. f_grnd(i,jp1) .gt. 0.0) then 
                    ! Floating to the bottom 
                    beta_acy(i,j) = (1.0-wt)*beta(i,j) + wt*beta(i,jp1)
                end if 

            end if

        end do 
        end do 

        return
        
    end subroutine stagger_beta_aa_gl_subgrid
    
    subroutine stagger_beta_aa_gl_subgrid_flux(beta_acx,beta_acy,beta,H_ice,f_ice,ux,uy,f_grnd,f_grnd_acx,f_grnd_acy)
        ! Modify basal friction coefficient by grounded fraction 
        ! weighted by linear-interpolated flux. 

        implicit none
        
        real(wp), intent(INOUT) :: beta_acx(:,:)      ! ac-nodes
        real(wp), intent(INOUT) :: beta_acy(:,:)      ! ac-nodes
        real(wp), intent(IN)    :: beta(:,:)          ! aa-nodes
        real(wp), intent(IN)    :: H_ice(:,:)         ! aa-nodes 
        real(wp), intent(IN)    :: f_ice(:,:)         ! aa-nodes 
        real(wp), intent(IN)    :: ux(:,:)            ! ac-nodes
        real(wp), intent(IN)    :: uy(:,:)            ! ac-nodes
        real(wp), intent(IN)    :: f_grnd(:,:)        ! aa-nodes     
        real(wp), intent(IN)    :: f_grnd_acx(:,:)    ! ac-nodes     
        real(wp), intent(IN)    :: f_grnd_acy(:,:)    ! ac-nodes     
        
        ! Local variables
        integer    :: i, j, nx, ny 
        integer    :: im1, ip1, jm1, jp1 
        real(wp)   :: ux_aa_a, ux_aa_b 
        real(wp)   :: uy_aa_a, uy_aa_b

        nx = size(beta_acx,1)
        ny = size(beta_acx,2) 

        ! Apply simple staggering to ac-nodes

        do j = 1, ny 
        do i = 1, nx

            im1 = max(1, i-1)
            ip1 = min(nx,i+1)
            
            jm1 = max(1, j-1)
            jp1 = min(ny,j+1)

            ! acx-nodes
            if (f_ice(i,j) .eq. 1.0 .and. f_ice(ip1,j) .eq. 1.0) then 
                ! Both aa-nodes are ice covered 

                if (f_grnd(i,j) .gt. 0.0 .and. f_grnd(ip1,j) .eq. 0.0) then 
                    ! Floating to the right

                    ux_aa_a = 0.5_wp*(ux(im1,j)+ux(i,j))
                    ux_aa_b = 0.5_wp*(ux(ip1,j)+ux(i,j))
                    
                    call calc_beta_gl_flux_weight(beta_acx(i,j),beta(i,j),beta(ip1,j), &
                                    ux_aa_a,ux_aa_b,H_ice(i,j),H_ice(ip1,j),f_grnd_acx(i,j))

                else if (f_grnd(i,j) .eq. 0.0 .and. f_grnd(ip1,j) .gt. 0.0) then 
                    ! Floating to the left 

                    ux_aa_a = 0.5_wp*(ux(ip1,j)+ux(i,j))
                    ux_aa_b = 0.5_wp*(ux(im1,j)+ux(i,j))
                    
                    call calc_beta_gl_flux_weight(beta_acx(i,j),beta(ip1,j),beta(i,j), &
                                    ux_aa_a,ux_aa_b,H_ice(ip1,j),H_ice(i,j),f_grnd_acx(i,j))

                end if 

            end if 

            ! acy-nodes
            if (f_ice(i,j) .eq. 1.0 .and. f_ice(i,jp1) .eq. 1.0) then 
                ! Both aa-nodes are ice covered 

                if (f_grnd(i,j) .gt. 0.0 .and. f_grnd(i,jp1) .eq. 0.0) then 
                    ! Floating to the top

                    uy_aa_a = 0.5_wp*(uy(i,jm1)+uy(i,j))
                    uy_aa_b = 0.5_wp*(uy(i,jp1)+uy(i,j))
                    
                    call calc_beta_gl_flux_weight(beta_acy(i,j),beta(i,j),beta(i,jp1), &
                                    uy_aa_a,uy_aa_b,H_ice(i,j),H_ice(i,jp1),f_grnd_acy(i,j))

                else if (f_grnd(i,j) .eq. 0.0 .and. f_grnd(ip1,j) .gt. 0.0) then 
                    ! Floating to the bottom 

                    uy_aa_a = 0.5_wp*(uy(i,jp1)+uy(i,j))
                    uy_aa_b = 0.5_wp*(uy(i,jm1)+uy(i,j))
                    
                    call calc_beta_gl_flux_weight(beta_acy(i,j),beta(i,jp1),beta(i,j), &
                                    uy_aa_a,uy_aa_b,H_ice(i,jp1),H_ice(i,j),f_grnd_acx(i,j))

                end if 

            end if 

        end do 
        end do 

        return
        
    end subroutine stagger_beta_aa_gl_subgrid_flux
    
    subroutine calc_beta_gl_flux_weight(beta_ac,beta_a,beta_b,ux_a,ux_b,H_a,H_b,f_grnd_ac)
        ! aa-node 'a' is grounded, aa-node 'b' is floating 
        ! Calculate weighted beta such that the weighting fraction 
        ! is scaled by u_g / u_tot, where u_g is sum of vel. for grounded segments 
        ! in cell between aa-nodes a and b, and u_tot is the sum of 
        ! all segments. Then:
        ! beta_ac = beta_a * (1-u_g/u_tot)

        ! Following "B2" approach proposed by Gladstone et al. (2010), 
        ! Eqs. 29, 30 & 31. 

        implicit none 

        real(wp), intent(OUT)   :: beta_ac
        real(wp), intent(IN)    :: beta_a
        real(wp), intent(IN)    :: beta_b
        real(wp), intent(IN)    :: ux_a
        real(wp), intent(IN)    :: ux_b 
        real(wp), intent(IN)    :: H_a 
        real(wp), intent(IN)    :: H_b  
        real(wp), intent(IN)    :: f_grnd_ac 

        ! Local variables 
        integer  :: i  
        real(wp) :: lambda 
        real(wp) :: q_a, q_b, q_now
        real(wp) :: H_now 
        real(wp) :: u_now 
        real(wp) :: uu_grnd, uu_tot 
        real(wp) :: weight 

        integer, parameter :: nseg = 100        ! 100 segments 
        
        ! First calculate flux and boundary aa-nodes 
        q_a = ux_a*H_a 
        q_b = ux_b*H_b 

        ! Set grounded and total flux sums to zero
        uu_grnd = 0.0 
        uu_tot  = 0.0 

        do i = 1, nseg 

            ! Get fraction along cell 
            lambda = real(i-1,wp)/real(nseg-1,wp)

            ! Get thickness and flux for current segment 
            H_now = H_a*lambda + H_b*(1.0_wp-lambda) 
            q_now = q_a*lambda + q_b*(1.0_wp-lambda)

            ! Recover velocity for current segment 
            if (H_now .gt. 0.0_wp) then 
                u_now = q_now / H_now 
            else
                u_now = 0.0_wp 
            end if 

            ! Add to total 
            uu_tot = uu_tot + q_now 

            ! If in grounded region, add to grounded total 
            if (lambda .lt. f_grnd_ac) then 
                uu_grnd = uu_grnd + q_now 
            end if 

        end do 

        ! Calculate grounded weight 
        if (uu_tot .gt. 0.0_wp) then 
            weight = uu_grnd / uu_tot 
        else
            ! Assume floating (no velocities or ice thicknesses found...)
            weight = 0.0_wp 
        end if 

        ! Calculate subgrid beta value 
        beta_ac = beta_a*weight 

        return 

    end subroutine calc_beta_gl_flux_weight

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

        real(wp), intent(IN) :: N_eff 
        real(wp), intent(IN) :: uxy_b 
        real(wp), intent(IN) :: ATT_base 
        real(wp), intent(IN) :: m_drag     
        real(wp), intent(IN) :: m_max 
        real(wp), intent(IN) :: lambda_max 
        real(wp) :: f_np 

        ! Local variables 
        real(wp) :: kappa 

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
