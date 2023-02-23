module yelmo_timesteps

    use yelmo_defs, only : sp, dp, wp, io_unit_err, ytime_class, MV, TOL_UNDERFLOW   
    use ncio 

    use topography, only : calc_ice_fraction, calc_H_grnd

    implicit none 

    private

    public :: set_pc_beta_coefficients

    public :: set_pc_mask
    public :: calc_pc_eta
    public :: calc_pc_tau_fe_sbe
    public :: calc_pc_tau_ab_sam
    public :: calc_pc_tau_heun
    public :: set_adaptive_timestep_pc
    
    public :: set_adaptive_timestep 
    public :: limit_adaptive_timestep

    public :: yelmo_timestep_write_init
    public :: yelmo_timestep_write

    public :: calc_adv2D_timestep1
    public :: calc_adv3D_timestep1 

    public :: check_checkerboard 

    public :: ytime_init 

contains

    subroutine set_pc_beta_coefficients(beta,pc_k,zeta,pc_method)

        implicit none 

        real(wp),         intent(OUT) :: beta(:)        ! Vector of 2 or 4 values 
        integer,          intent(OUT) :: pc_k           ! Order of the method
        real(wp),         intent(IN)  :: zeta           ! Ratio of current to previous timestep
        character(len=*), intent(IN)  :: pc_method 

        ! Local variables

        ! Using a small value of eps, like eps=0.1, may help improve stability further
        ! of the AB method. See this page for more details:
        ! https://mitgcm.readthedocs.io/en/latest/algorithm/algorithm.html#explicit-time-stepping-adams-bashforth
        real(wp), parameter :: eps = 0.0

        ! Determine which predictor-corrector (pc) method we are using for timestepping,
        ! assign scheme order and weights 

        if (size(beta,1) .eq. 4) then 
            ! ======== TWO-STEP METHODS with 4 beta values ========

            select case(trim(pc_method))

                case("FE-SBE")
                    
                    ! Order of the method 
                    pc_k = 2 

                    beta(1) = 1.0_wp 
                    beta(2) = 0.0_wp 
                    
                    beta(3) = 1.0_wp 
                    beta(4) = 0.0_wp 
                    
                case("AB-SAM")
                    
                    ! Order of the method 
                    pc_k = 2 

                    beta(1) = 1.0_wp + zeta/2.0_wp + eps
                    beta(2) = -zeta/2.0_wp - eps

                    beta(3) = 0.5_wp 
                    beta(4) = 0.5_wp 
                
                case("HEUN")
                    
                    ! Order of the method 
                    pc_k = 2 

                    beta(1) = 1.0_wp 
                    beta(2) = 0.0_wp 
                    
                    beta(3) = 0.5_wp 
                    beta(4) = 0.5_wp 
                
                case("RALSTON")
                    
                    write(io_unit_err,*) "This method does not work yet - the truncation error is incorrect."
                    stop 

                    ! Order of the method 
                    pc_k = 2 

                    beta(1) = 2.0_wp / 3.0_wp
                    beta(2) = 0.0_wp 
                    
                    beta(3) = 0.25_wp 
                    beta(4) = 0.75_wp 
                
                case DEFAULT 

                    write(io_unit_err,*) "set_pc_beta_coefficients:: &
                        &Error: two-step pc_method does not match available options [FE-SBE, AB-SAM, HEUN]."
                    write(io_unit_err,*) "pc_method = ", trim(pc_method)
                    stop 

            end select 
            
        else 
            ! ======== ONE-STEP METHODS with 2 beta values ========
            ! Note: also set pc_k here for consistency, but order 
            ! of method is not critical here since PC adaptive timestepping
            ! is not used with one-step methods. 
            
            select case(trim(pc_method))

                case("FE")
                    ! Forward Euler 

                    pc_k = 1

                    beta(1) = 1.0_wp  
                    beta(2) = 0.0_wp  

                case("AB") 
                    ! Adams-Bashforth

                    pc_k = 2
                    
                    beta(1) = 1.0_wp + zeta/2.0_wp 
                    beta(2) = -zeta/2.0_wp 

                case("SAM") 
                    ! Semi-implicit Adams–Moulton

                    pc_k = 2
                    
                    beta(1) = 0.5
                    beta(2) = 0.5 

                case DEFAULT 

                    write(io_unit_err,*) "set_pc_beta_coefficients:: &
                        &Error: one-step pc_method does not match available options [FE, SBE, AB, SAM]."
                    write(io_unit_err,*) "thrm:: dt_method = ", trim(pc_method)
                    stop 

            end select 

        end if 

        return

    end subroutine set_pc_beta_coefficients


    subroutine set_pc_mask(mask,pc_tau,H_ice_pred,H_ice_corr,z_bed,z_sl,margin_flt_subgrid)

        implicit none 

        logical, intent(OUT) :: mask(:,:) 
        real(wp), intent(IN) :: pc_tau(:,:) 
        real(wp), intent(IN) :: H_ice_pred(:,:) 
        real(wp), intent(IN) :: H_ice_corr(:,:) 
        real(wp), intent(IN) :: z_bed(:,:) 
        real(wp), intent(IN) :: z_sl(:,:) 
        logical,    intent(IN) :: margin_flt_subgrid 

        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: im1, jm1, ip1, jp1 

        real(wp), allocatable :: f_ice_pred(:,:) 
        real(wp), allocatable :: f_ice_corr(:,:) 
        real(wp), allocatable :: H_grnd_pred(:,:) 
        real(wp), allocatable :: H_grnd_corr(:,:) 

        real(wp), parameter :: H_lim = 10.0      ! [m] 
        
        nx = size(mask,1)
        ny = size(mask,2) 

        allocate(f_ice_pred(nx,ny))
        allocate(f_ice_corr(nx,ny))
        allocate(H_grnd_pred(nx,ny))
        allocate(H_grnd_corr(nx,ny))
        
        ! Get the ice area fraction mask for each ice thickness map 
        call calc_ice_fraction(f_ice_pred,H_ice_pred,z_bed,z_sl,margin_flt_subgrid)
        call calc_ice_fraction(f_ice_corr,H_ice_corr,z_bed,z_sl,margin_flt_subgrid)

        ! Get the grounded ice distance to flotation 
        call calc_H_grnd(H_grnd_pred,H_ice_pred,f_ice_pred,z_bed,z_sl)
        call calc_H_grnd(H_grnd_corr,H_ice_corr,f_ice_corr,z_bed,z_sl)

        ! Initially set all points to True 
        mask = .TRUE. 
    
if (.TRUE.) then 
        do j = 1, ny 
        do i = 1, nx

            im1 = max(i-1,1)
            jm1 = max(j-1,1)
            ip1 = min(i+1,nx)
            jp1 = min(j+1,ny)

            ! Define places that should not be checked 

            if (H_ice_pred(i,j) .lt. H_lim .or. H_ice_corr(i,j) .lt. H_lim) then 
                ! (Near) ice-free point or may be transitioning state

                mask(i,j) = .FALSE. 

            ! ajr: excluding thin points doesn't help...
            !else if (H_ice(i,j) .lt. 50.0_wp) then 
            !    ! Thin ice-covered point 
            !
            !    mask(i,j) = .FALSE. 
            !
            else
                ! Ice-covered points, further checks below

                 if (count(f_ice_pred(im1:ip1,jm1:jp1).lt.1.0) .gt. 0 .or. &
                     count(f_ice_corr(im1:ip1,jm1:jp1).lt.1.0) .gt. 0) then 
                    ! Point is at (or near) ice-margin point 

                    mask(i,j) = .FALSE. 

                else if (H_grnd_pred(i,j) .le. 0.0_wp .or. H_grnd_corr(i,j) .le. 0.0_wp) then 
                    ! Grounding-line or floating ice point 

                    mask(i,j) = .FALSE. 

!                 else if (H_grnd(i,j) .gt. 0.0_wp .and. &
!                             count(H_grnd(im1:ip1,jm1:jp1).lt.1.0_wp) .le. 0) then 
!                     ! Neighbor at the grounding line

!                     mask(i,j) = .FALSE. 

                end if 

            end if 

        end do 
        end do  

end if 

        return 

    end subroutine set_pc_mask

    function calc_pc_eta(tau,mask) result(eta)

        implicit none 

        real(wp), intent(IN) :: tau(:,:) 
        logical,    intent(IN) :: mask(:,:) 
        real(wp) :: eta 

        real(wp), parameter :: eta_tol = 1e-8 

        integer :: npts

if (.TRUE.) then
        ! Calculate eta 
        eta = maxval(abs(tau),mask=mask)

        ! Limit to non-zero value
        ! Note: Limiting minimum to above eg 1e-8 is very 
        ! important for reducing fluctuations in dt 
        eta = max(eta,eta_tol)

else
    ! ajr: testing rmse(tau) instead of max(tau)
    ! So far, this works, but leads to large areas of Antarctica on the coast with large tau values.
        npts = count(mask)
        if (npts .eq. 0) then 
            eta = eta_tol 
        else 
            eta = sqrt(sum(tau**2,mask=mask)/real(npts,wp))
            eta = max(eta,eta_tol)
        end if
end if 

        return 

    end function calc_pc_eta
    
    elemental subroutine calc_pc_tau_fe_sbe(tau,var_corr,var_pred,dt_n)
        ! Calculate truncation error for the FE-SBE timestepping method
        ! Forward Euler (FE) predictor step and Semi-implicit
        ! Backward Euler (SBE) corrector step. 
        ! Implemented followig Cheng et al (2017, GMD)
        ! Truncation error: tau = 1/2*dt_n * (var - var_pred)

        implicit none 

        real(wp), intent(OUT) :: tau
        real(wp), intent(IN)  :: var_corr
        real(wp), intent(IN)  :: var_pred
        real(wp), intent(IN)  :: dt_n 
        
        if (dt_n .eq. 0.0_wp) then 
            tau = 0.0_wp 
        else 
            tau = (1.0_wp / (2.0_wp*dt_n)) * (var_corr - var_pred)
        end if 

        ! Keep tau positive to make crash analysis simpler
        tau = abs(tau) 

        return 

    end subroutine calc_pc_tau_fe_sbe

    elemental subroutine calc_pc_tau_ab_sam(tau,var_corr,var_pred,dt_n,zeta)
        ! Calculate truncation error for the AB-SAM timestepping method
        ! Adams-Bashforth (AB) predictor step and Semi-implicit
        ! Adams–Moulton (SAM) corrector step. 
        ! Implemented followig Cheng et al (2017, GMD)
        ! Truncation error: tau = zeta * (var - var_pred) / ((3zeta + 3)*dt)

        implicit none 

        real(wp), intent(OUT) :: tau
        real(wp), intent(IN)  :: var_corr
        real(wp), intent(IN)  :: var_pred
        real(wp), intent(IN)  :: dt_n 
        real(wp), intent(IN)  :: zeta 

        if (dt_n .eq. 0.0_wp) then 
            tau = 0.0_wp 
        else 
            tau = zeta * (var_corr - var_pred) / ((3.0_wp*zeta + 3.0_wp)*dt_n)
        end if 

        ! Keep tau positive to make crash analysis simpler
        tau = abs(tau) 
        
        return 

    end subroutine calc_pc_tau_ab_sam

    elemental subroutine calc_pc_tau_heun(tau,var_corr,var_pred,dt_n)
        ! Calculate truncation error for Heun's timestepping method
        ! as derived by Marisa Montoya 

        implicit none 

        real(wp), intent(OUT) :: tau
        real(wp), intent(IN)  :: var_corr
        real(wp), intent(IN)  :: var_pred
        real(wp), intent(IN)  :: dt_n 
        
        if (dt_n .eq. 0.0_wp) then 
            tau = 0.0_wp 
        else 
            tau = (1.0_wp / (6.0_wp*dt_n)) * (var_corr - var_pred)
        end if 

        ! Keep tau positive to make crash analysis simpler
        tau = abs(tau) 
        
        return 

    end subroutine calc_pc_tau_heun

    subroutine set_adaptive_timestep_pc(dt_new,dt,eta,eps,dtmin,dtmax,ux_bar,uy_bar,dx,pc_k,controller)
        ! Calculate the timestep following algorithm for 
        ! a general predictor-corrector (pc) method.
        ! Implemented followig Cheng et al (2017, GMD)

        implicit none 

        real(wp), intent(OUT) :: dt_new               ! [yr]   Timestep (n+1)
        real(wp), intent(IN)  :: dt(:)                ! [yr]   Timesteps (n:n-2)
        real(wp), intent(IN)  :: eta(:)               ! [X/yr] Maximum truncation error (n:n-2)
        real(wp), intent(IN)  :: eps                  ! [--]   Tolerance value (eg, eps=1e-4)
        real(wp), intent(IN)  :: dtmin                ! [yr]   Minimum allowed timestep
        real(wp), intent(IN)  :: dtmax                ! [yr]   Maximum allowed timestep
        real(wp), intent(IN)  :: ux_bar(:,:)          ! [m/yr]
        real(wp), intent(IN)  :: uy_bar(:,:)          ! [m/yr]
        real(wp), intent(IN)  :: dx                   ! [m]
        integer,    intent(IN)  :: pc_k                 ! pc_k gives the order of the timestepping scheme (pc_k=2 for FE-SBE, pc_k=3 for AB-SAM)
        character(len=*), intent(IN) :: controller      ! Adaptive controller to use [PI42, H312b, H312PID]

        ! Local variables
        real(wp) :: dt_n, dt_nm1, dt_nm2          ! [yr]   Timesteps (n:n-2)
        real(wp) :: eta_n, eta_nm1, eta_nm2       ! [X/yr] Maximum truncation error (n:n-2)
        real(wp) :: rho_n, rho_nm1, rho_nm2
        real(wp) :: rhohat_n 
        real(wp) :: dt_adv 
        real(wp) :: dtmax_now
        real(wp) :: k_i 
        real(wp) :: k_p, k_d 

        ! Smoothing parameter; Söderlind and Wang (2006) method, Eq. 10
        ! Values on the order of [0.7,2.0] are reasonable. Higher kappa slows variation in dt
        real(wp), parameter :: kappa = 2.0_wp 
        
        ! Step 1: Save information needed for adapative controller algorithms 

        ! Save dt from several timesteps (potentially more available)
        dt_n    = max(dt(1),dtmin) 
        dt_nm1  = max(dt(2),dtmin) 
        dt_nm2  = max(dt(3),dtmin)

        ! Save eta from several timesteps (potentially more available)
        eta_n   = eta(1)
        eta_nm1 = eta(2)
        eta_nm2 = eta(3)

        ! Calculate rho from several timesteps 
        rho_nm1 = (dt_n   / dt_nm1) 
        rho_nm2 = (dt_nm1 / dt_nm2) 

        ! Step 2: calculate scaling for the next timestep (dt,n+1)
        select case(trim(controller))

            case("PI42")
                ! Söderlind and Wang, 2006; Cheng et al., 2017
                ! Deeper discussion in Söderlind, 2002. Note for example, 
                ! that Söderlind (2002) recommends:
                ! k*k_i >= 0.3 
                ! k*k_p >= 0.2 
                ! k*k_i + k*k_p <= 0.7 
                ! However, the default values of Söderlind and Wang (2006) and Cheng et al (2017)
                ! are outside of these bounds. 

                ! Default parameter values 
                k_i = 2.0_wp / (pc_k*5.0_wp)
                k_p = 1.0_wp / (pc_k*5.0_wp)
                
                ! Improved parameter values (reduced oscillations)
!                 k_i = 4.0_wp / (pc_k*10.0_wp)
!                 k_p = 3.0_wp / (pc_k*10.0_wp)

                ! Experimental parameter values (minimal oscillations, does not access largest timesteps)
                ! k_i = 0.5_wp / (pc_k*10.0_wp)
                ! k_p = 6.5_wp / (pc_k*10.0_wp)

                ! Default parameter values
                rho_n = calc_pi_rho_pi42(eta_n,eta_nm1,rho_nm1,eps,k_i,k_p,alpha_2=0.0_wp)

            case("H312b") 
                ! Söderlind (2003) H312b, Eq. 31+ (unlabeled) 
                
                rho_n = calc_pi_rho_H312b(eta_n,eta_nm1,eta_nm2,rho_nm1,rho_nm2,eps,k=real(pc_k,wp),b=8.0_wp)

            case("H312PID") 
                ! Söderlind (2003) H312PD, Eq. 38
                ! Note: Suggested k_i =(2/9)*1/pc_k, but lower value gives more stable solution

                !k_i = (2.0_wp/9.0_wp)*1.0_wp/real(pc_k,wp)
                k_i = 0.08_wp/real(pc_k,wp)

                rho_n = calc_pi_rho_H312PID(eta_n,eta_nm1,eta_nm2,eps,k_i)

            case("H321PID")

                k_i = 0.1  / real(pc_k,wp)
                k_p = 0.45 / real(pc_k,wp) 

                rho_n = calc_pi_rho_H321PID(eta_n,eta_nm1,eta_nm2,dt_n,dt_nm1,eps,k_i,k_p)
                
            case("PID1")

                k_i = 0.175 
                k_p = 0.075
                k_d = 0.01 

                rho_n = calc_pi_rho_PID1(eta_n,eta_nm1,eta_nm2,eps,k_i,k_p,k_d)
                
            case DEFAULT 

                write(*,*) "set_adaptive_timestep_pc:: Error: controller not recognized."
                write(*,*) "controller = ", trim(controller) 
                stop 

        end select 

        ! Scale rho_n for smoothness 
        rhohat_n = rho_n
        !rhohat_n = min(rho_n,1.1)
        !rhohat_n = 1.0_wp + kappa * atan((rho_n-1.0_wp)/kappa) ! Söderlind and Wang, 2006, Eq. 10
        
        ! Step 3: calculate the next time timestep (dt,n+1)
        dt_new = rhohat_n * dt_n

        ! Step 4: Modify timestep to fit within prescribed limits 

        ! Calculate CFL advection limit too, and limit maximum allowed timestep
        dt_adv    = minval( calc_adv2D_timestep1(ux_bar,uy_bar,dx,dx,cfl_max=1.0_wp) ) 
        dtmax_now = min(dtmax,dt_adv) 

        ! Finally, ensure timestep is within prescribed limits
        call limit_adaptive_timestep(dt_new,dtmin,dtmax_now)
        
        return 

    end subroutine set_adaptive_timestep_pc

    function calc_pi_rho_pi42(eta_n,eta_nm1,rho_nm1,eps,k_i,k_p,alpha_2) result(rho_n)

        implicit none 

        real(wp), intent(IN) :: eta_n 
        real(wp), intent(IN) :: eta_nm1 
        real(wp), intent(IN) :: rho_nm1 
        real(wp), intent(IN) :: eps 
        real(wp), intent(IN) :: k_i 
        real(wp), intent(IN) :: k_p
        real(wp), intent(IN) :: alpha_2 
        real(wp) :: rho_n 

        ! Söderlind and Wang, 2006; Cheng et al., 2017
        ! Original formulation: Söderlind, 2002, Eq. 3.12:
        rho_n   = (eps/eta_n)**(k_i+k_p) * (eps/eta_nm1)**(-k_p) * rho_nm1**(-alpha_2)

        return 

    end function calc_pi_rho_pi42

    function calc_pi_rho_H312b(eta_n,eta_nm1,eta_nm2,rho_nm1,rho_nm2,eps,k,b) result(rho_n)

        implicit none 

        real(wp), intent(IN) :: eta_n 
        real(wp), intent(IN) :: eta_nm1 
        real(wp), intent(IN) :: eta_nm2 
        real(wp), intent(IN) :: rho_nm1
        real(wp), intent(IN) :: rho_nm2
        real(wp), intent(IN) :: eps 
        real(wp), intent(IN) :: k 
        real(wp), intent(IN) :: b 
        real(wp) :: rho_n 

        ! Local variables 
        real(wp) :: beta_1, beta_2, beta_3 
        real(wp) :: alpha_2, alpha_3 

        beta_1  =  1.0_wp / (k*b)
        beta_2  =  2.0_wp / (k*b)
        beta_3  =  1.0_wp / (k*b)
        alpha_2 = -3.0_wp / b 
        alpha_3 = -1.0_wp / b 

        ! Söderlind (2003) H312b, Eq. 31+ (unlabeled) 
        rho_n   = (eps/eta_n)**beta_1 * (eps/eta_nm1)**beta_2 * (eps/eta_nm2)**beta_3 &
                            * rho_nm1**alpha_2 * rho_nm2**alpha_3 

        return 

    end function calc_pi_rho_H312b

    function calc_pi_rho_H312PID(eta_n,eta_nm1,eta_nm2,eps,k_i) result(rho_n)

        implicit none 

        real(wp), intent(IN) :: eta_n 
        real(wp), intent(IN) :: eta_nm1 
        real(wp), intent(IN) :: eta_nm2 
        real(wp), intent(IN) :: eps 
        real(wp), intent(IN) :: k_i  
        real(wp) :: rho_n 

        ! Local variables 
        real(wp) :: k_i_1, k_i_2, k_i_3

        k_i_1   = k_i / 4.0_wp 
        k_i_2   = k_i / 2.0_wp 
        k_i_3   = k_i / 4.0_wp 

        ! Söderlind (2003) H312PID, Eq. 38
        rho_n   = (eps/eta_n)**k_i_1 * (eps/eta_nm1)**k_i_2 * (eps/eta_nm2)**k_i_3

        return 

    end function calc_pi_rho_H312PID

    function calc_pi_rho_H321PID(eta_n,eta_nm1,eta_nm2,dt_n,dt_nm1,eps,k_i,k_p) result(rho_n)

        implicit none 

        real(wp), intent(IN) :: eta_n 
        real(wp), intent(IN) :: eta_nm1 
        real(wp), intent(IN) :: eta_nm2 
        real(wp), intent(IN) :: dt_n 
        real(wp), intent(IN) :: dt_nm1 
        real(wp), intent(IN) :: eps 
        real(wp), intent(IN) :: k_i  
        real(wp), intent(IN) :: k_p
        real(wp) :: rho_n 

        ! Local variables 
        real(wp) :: k_i_1, k_i_2, k_i_3

        k_i_1   =   0.75_wp*k_i + 0.50_wp*k_p 
        k_i_2   =   0.50_wp*k_i 
        k_i_3   = -(0.25_wp*k_i + 0.50_wp*k_p)

        ! Söderlind (2003) H321PID, Eq. 42
        rho_n   = (eps/eta_n)**k_i_1 * (eps/eta_nm1)**k_i_2 * (eps/eta_nm2)**k_i_3 * (dt_n / dt_nm1)

        return 

    end function calc_pi_rho_H321PID 

    function calc_pi_rho_PID1(eta_n,eta_nm1,eta_nm2,eps,k_i,k_p,k_d) result(rho_n)

        implicit none 

        real(wp), intent(IN) :: eta_n 
        real(wp), intent(IN) :: eta_nm1 
        real(wp), intent(IN) :: eta_nm2 
        real(wp), intent(IN) :: eps 
        real(wp), intent(IN) :: k_i  
        real(wp), intent(IN) :: k_p
        real(wp), intent(IN) :: k_d
        real(wp) :: rho_n 

        ! https://www.mathematik.uni-dortmund.de/~kuzmin/cfdintro/lecture8.pdf
        ! Page 20 (theoretical basis unclear/unknown)
        rho_n   = (eps/eta_n)**k_i * (eta_nm1/eta_n)**k_p * (eta_nm1**2/(eta_n*eta_nm2))**k_d

        return 

    end function calc_pi_rho_PID1

    subroutine set_adaptive_timestep(dt,dt_adv,dt_diff,dt_adv3D, &
                        ux,uy,uz,ux_bar,uy_bar,H_ice,dHicedt,zeta_ac, &
                        dx,dtmin,dtmax,cfl_max,cfl_diff_max)
        ! Determine value of adaptive timestep to be consistent with 
        ! min/max timestep range and maximum allowed step of model 
        ! to line up with control time steps

        implicit none 

        real(wp), intent(OUT) :: dt             ! [a] Current timestep 
        real(wp), intent(OUT) :: dt_adv(:,:)     ! [a] Diagnosed maximum advective timestep (vertical ave)
        real(wp), intent(OUT) :: dt_diff(:,:)    ! [a] Diagnosed maximum diffusive timestep (vertical ave) 
        real(wp), intent(OUT) :: dt_adv3D(:,:,:) ! [a] Diagnosed maximum advective timestep (3D) 
        real(wp), intent(IN)  :: ux(:,:,:)       ! [m a-1]
        real(wp), intent(IN)  :: uy(:,:,:)       ! [m a-1]
        real(wp), intent(IN)  :: uz(:,:,:)       ! [m a-1]
        real(wp), intent(IN)  :: ux_bar(:,:)     ! [m a-1]
        real(wp), intent(IN)  :: uy_bar(:,:)     ! [m a-1]
        real(wp), intent(IN)  :: H_ice(:,:)      ! [m]
        real(wp), intent(IN)  :: dHicedt(:,:)    ! [m a-1]
        real(wp), intent(IN)  :: zeta_ac(:)      ! [--] 
        real(wp), intent(IN)  :: dx, dtmin, dtmax ! [a]
        real(wp), intent(IN)  :: cfl_max
        real(wp), intent(IN)  :: cfl_diff_max
        
        ! Local variables 
        real(wp) :: dt_adv_min, dt_diff_min 
        real(wp) :: x 
        logical    :: is_unstable
        real(wp), parameter :: dtmax_cfl   = 20.0_wp 
        real(wp), parameter :: exp_cfl     =  2.0_wp 
        real(wp), parameter :: rate_lim    = 1.0_wp   ! Reduction in timestep for instability 
        real(wp), parameter :: rate_scalar = 0.05_wp  ! Reduction in timestep for instability 

        ! Timestep limits determined from CFL conditions for general advective
        ! velocity, as well as diagnosed diffusive magnitude
        ! (adapted from Bueler et al., 2007)

        dt_adv   = calc_adv2D_timestep1(ux_bar,uy_bar,dx,dx,cfl_max)
        !dt_diff  = calc_diff2D_timestep(D2D,dx,dx,cfl_diff_max) 
        dt_diff = 1000.0     ! Prescribe something just to avoid compiler warnings 
        ! ajr: diffusivity D2D is not available right now (SIA solver does not calculate it)
        ! So, diffusivity needs to be diagnosed, to be able to estimate dt_diff properly. 

!         dt_adv3D = calc_adv3D_timestep1(ux,uy,uz,dx,dx,H_ice,zeta_ac,cfl_max)
!         dt_adv3D = calc_adv3D_timestep(ux,uy,uz,H_ice,zeta_ac,dx,dx,cfl_max)
        dt_adv3D = 1000.0    ! Prescribe something just to avoid compiler warnings 

        ! Get minimum from adv and diffusive timesteps
        dt_adv_min  = minval(dt_adv)
        dt_diff_min = minval(dt_diff)
        
        ! Note: It's not clear whether dt_diff is working well, so for now
        ! it is not applied as a limit. Furthermore, although dt_adv should
        ! be consistent with the CFL limit, it does not guarantee that a 
        ! fully coupled thermodynamic model will remain stable. 

        ! Choose minimum timestep between advective and diffusive limits 
        !dt = min(dt_adv_min,dt_diff_min)
        dt = dt_adv_min 

        ! Apply additional reduction in timestep as it gets smaller
        if (.FALSE.) then 
            dt = max(dtmin,dt)          ! dt >= dtmin
            dt = min(dtmax_cfl,dt)      ! dt <= dtmax_cfl
            x = ((dt-dtmin)/(dtmax_cfl-dtmin))**exp_cfl
            dt = dtmin + (dtmax-dtmin)*x
        end if 

        ! Check if additional timestep reduction is necessary,
        ! due to checkerboard patterning related to mass conservation.
        ! Reduce if necessary 
        call check_checkerboard(is_unstable,dHicedt,rate_lim)
        if (is_unstable) dt = rate_scalar*dt

        ! Finally, ensure timestep is within prescribed limits
        call limit_adaptive_timestep(dt,dtmin,dtmax)

        return 

    end subroutine set_adaptive_timestep
    
    subroutine limit_adaptive_timestep(dt,dtmin,dtmax)
        ! Make sure that adaptive timestep is in range of dtmin < dt < dtmax 
        ! where dtmax is evolving to arrive at final timestep, eg time + dtmax = time_max 

        implicit none

        real(wp), intent(INOUT) :: dt               ! [a] Current timestep 
        real(wp), intent(IN)    :: dtmin            ! [a] Minimum allowed timestep
        real(wp), intent(IN)    :: dtmax            ! [a] Maximum allowed timestep 

        ! Local variables  
        real(wp), parameter :: n_decimal   = 6          ! Maximum decimals to treat for timestep
        real(wp), parameter :: dt_half_lim = 0.5_wp   ! Should be 0.5 or greater to make sense

        ! Check to avoid lopsided timesteps (1 big, 1 tiny) to arrive at time_max  
        if (dtmax .gt. 0.0) then 

            ! Ensure timestep is also within parameter limits 
            dt = max(dtmin,dt)  ! dt >= dtmin
            dt = min(dtmax,dt)  ! dt <= dtmax

            if (dt/dtmax .gt. dt_half_lim .and. dt .lt. dtmax) then 
                ! Current adaptive timestep is greater than ~0.5 of the total
                ! expected timestep, and another timestep will be needed to
                ! reach time_max. Therefore, set this timestep to a smaller
                ! value, ie, dt = dt_half_lim*dtmax. 

                dt = dt_half_lim*dtmax

            else if (dt/dtmax .lt. dt_half_lim) then 
                ! Round-off extra digits for neatness

                dt = real(floor(dt*10.0_wp**n_decimal)*10.0_wp**(-n_decimal), wp)
                
            end if 

        else 
            ! dt is simply zero 

            dt = dtmax 

        end if 

        return 

    end subroutine limit_adaptive_timestep

    subroutine set_to_nearest_timestep(dt)
        ! ajr: limit timestep to specific list of values (for neatness)
        ! note: routine not used or tested yet! 

        implicit none 

        real(wp), intent(INOUT) :: dt 

        ! Local variables 
        integer  :: i  

        integer, parameter :: n = 22 
        real(wp), parameter :: dt_set(n) = &
        [100.0,50.0,20.0,10.0,5.0,2.0,1.0,0.5,0.2,0.1,0.05,0.02, &
         0.01,0.005,0.002,0.001,0.0005,0.0002,0.0001,0.00005,0.00002,0.00001]

        do i = 1, n 
            if (dt_set(i) .le. dt) then 
                dt = dt_set(i)
                exit 
            end if
        end do 

        return

    end subroutine set_to_nearest_timestep



    elemental function calc_diff2D_timestep(D,dx,dy,cfl_diff_max) result(dt)
        ! Calculate maximum diffusion time step based
        ! on Courant–Friedrichs–Lewy condition
        ! Equation obtained from Bueler et al. (2007), Eq. 25:
        ! dt/2 * (1/dx^2 + 1/dy^2)*max(D) <= cfl_diff_max = 0.12 
        ! dt = cfl_diff_max * 2.0 / ((1/dx^2+1/dy^2)*max(D))

        implicit none 
        
        real(wp), intent(IN) :: D, dx, dy
        real(wp), intent(IN) :: cfl_diff_max       ! Maximum Courant number, default cfl_diff_max=0.12
        real(wp) :: dt 

        dt = (2.0*cfl_diff_max) / ((1.0/(dx**2)+1.0/(dy**2))*max(abs(D),1e-5))
        
        return 

    end function calc_diff2D_timestep
    
    function calc_adv2D_timestep1(ux,uy,dx,dy,cfl_max) result(dt)
        ! Calculate maximum advective time step based
        ! on Courant–Friedrichs–Lewy condition
        ! https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition

        ! 1D condition: C = u*dt/dx <= cfl_max 
        ! 2D condition: C = u*dt/dx + v*dt/dy <= cfl_max 
        ! thus when C = cfl_max:
        ! dt = cfl_max * 1/(u/dx+v/dx)


        implicit none 
        
        real(wp), intent(IN) :: ux(:,:)           ! acx-nodes
        real(wp), intent(IN) :: uy(:,:)           ! acy-nodes 
        real(wp), intent(IN) :: dx, dy
        real(wp), intent(IN) :: cfl_max           ! Maximum Courant number, default cfl_max=1.0
        real(wp) :: dt(size(ux,1),size(ux,2))     ! aa-nodes 

        ! Local variables  
        integer :: i, j, nx, ny 
        real(wp) :: ux_now, uy_now 

        real(wp), parameter :: eps = 1e-1         ! [m/a] Small factor to avoid divide by zero 

        nx = size(ux,1)
        ny = size(ux,2)

        do j = 2, ny-1 
        do i = 2, nx-1 

!             ux_now = abs( 0.5*(ux(i-1,j)+ux(i,j)) )
!             uy_now = abs( 0.5*(uy(i,j-1)+uy(i,j)) )

            ux_now = max(abs(ux(i-1,j)),abs(ux(i,j)))
            uy_now = max(abs(uy(i,j-1)),abs(uy(i,j)))
            
            if (abs(ux_now) .lt. TOL_UNDERFLOW) ux_now = 0.0_wp 
            if (abs(uy_now) .lt. TOL_UNDERFLOW) uy_now = 0.0_wp 

            dt(i,j) = cfl_max * 1.0 / (ux_now/dx + uy_now/dy + eps/dx)
            
            ! Underflow issues:
            !dt(i,j) = cfl_max * 1.0 / (abs(ux(i-1,j))/dx + abs(ux(i,j))/dx &
            !                           + abs(uy(i,j-1))/dy + abs(uy(i,j))/dy + eps/(dx+dy))
        end do 
        end do 

        dt(1,:)  = dt(2,:)
        dt(nx,:) = dt(nx-1,:) 
        dt(:,1)  = dt(:,2)
        dt(:,ny) = dt(:,ny-1)

        return 

    end function calc_adv2D_timestep1

    function calc_adv3D_timestep1(ux,uy,uz,dx,dy,H_ice,zeta_ac,cfl_max) result(dt)
        ! Calculate maximum advective time step based
        ! on Courant–Friedrichs–Lewy condition
        ! https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition

        ! 1D condition: C = u*dt/dx <= cfl_max 
        ! 2D condition: C = u*dt/dx + v*dt/dy <= cfl_max 
        ! thus when C = cfl_max:
        ! dt = cfl_max * 1/(u/dx+v/dx)


        implicit none 
        
        real(wp), intent(IN) :: ux(:,:,:)        ! acx-nodes
        real(wp), intent(IN) :: uy(:,:,:)        ! acy-nodes
        real(wp), intent(IN) :: uz(:,:,:)        ! acz-nodes  
        real(wp), intent(IN) :: dx, dy
        real(wp), intent(IN) :: H_ice(:,:)       ! aa-nodes 
        real(wp), intent(IN) :: zeta_ac(:)       ! ac-nodes 
        real(wp), intent(IN) :: cfl_max          ! Maximum Courant number, default cfl_max=1.0
        real(wp) :: dt(size(ux,1),size(ux,2),size(ux,3))    ! aa-nodes 

        ! Local variables  
        integer :: i, j, k, nx, ny, nz_aa  
        real(wp) :: ux_now, uy_now 
        real(wp) :: dz 

        nx    = size(ux,1)
        ny    = size(ux,2)
        nz_aa = size(zeta_ac)+1  

        ! Set a high timestep to start 
        dt = cfl_max * 1.0 / (1e-3)

        do j = 2, ny-1 
        do i = 2, nx-1 

            if (H_ice(i,j) .gt. 0.0) then 

!             ux_now = abs( 0.5*(ux(i-1,j)+ux(i,j)) )
!             uy_now = abs( 0.5*(uy(i,j-1)+uy(i,j)) )

            !ux_now = max(abs(ux(i-1,j)),abs(ux(i,j)))
            !uy_now = max(abs(uy(i,j-1)),abs(uy(i,j)))
            
            !dt(i,j) = cfl_max * 1.0 / max(ux_now/dx + uy_now/dy,1e-3)

!             dt(i,j) = cfl_max * 1.0 / max(abs(ux(i-1,j))/dx + abs(ux(i,j))/dx &
!                                         + abs(uy(i,j-1))/dy + abs(uy(i,j))/dy,1e-3)
            
                do k = 2, nz_aa-1 

                    dz = H_ice(i,j) * (zeta_ac(k)-zeta_ac(k-1))

                    dt(i,j,k) = cfl_max * 1.0 / max(abs(ux(i-1,j,k))/(2.0*dx) + abs(ux(i,j,k))/(2.0*dx) &
                                            + abs(uy(i,j-1,k))/(2.0*dy) + abs(uy(i,j,k))/(2.0*dy), &
                                            + abs(uz(i,j,k-1))/(2.0*dz) + abs(uz(i,j,k))/(2.0*dz), 1e-3)
    !                 dt(i,j,k) = cfl_max * 1.0 / max(abs(ux(i-1,j,k))/(2.0*dx) + abs(ux(i,j,k))/(2.0*dx) &
    !                                         + abs(uy(i,j-1,k))/(2.0*dy) + abs(uy(i,j,k))/(2.0*dy), 1e-3)
                
                end do 

            end if 

        end do 
        end do 

        return 

    end function calc_adv3D_timestep1
    
    elemental function calc_adv2D_timestep(ux,uy,dx,dy,cfl_max) result(dt)
        ! Calculate maximum advective time step based
        ! on Courant–Friedrichs–Lewy condition
        ! https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition

        ! 1D condition: C = u*dt/dx <= cfl_max 
        ! 2D condition: C = u*dt/dx + v*dt/dy <= cfl_max 
        ! thus when C = cfl_max:
        ! dt = cfl_max * 1/(u/dx+v/dx)


        implicit none 
        
        real(wp), intent(IN) :: ux
        real(wp), intent(IN) :: uy
        real(wp), intent(IN) :: dx, dy
        real(wp), intent(IN) :: cfl_max             ! Maximum Courant number, default cfl_max=1.0
        real(wp) :: dt 

        dt = cfl_max * 1.0 / max(abs(ux)/dx + abs(uy)/dy,1e-5)

        return 

    end function calc_adv2D_timestep
    
    subroutine calc_adv2D_velocity(ux,uy,dx,dy,dt,cfl_max)
        ! Calculate maximum velocity given a known time step
        ! on Courant–Friedrichs–Lewy condition
        ! https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition

        ! 1D condition: C = u*dt/dx <= cfl_max 
        ! 2D condition: C = u*dt/dx + v*dt/dy <= cfl_max 
        ! thus when C = cfl_max:
        ! dt = cfl_max * 1/(u/dx+v/dx)

        ! dt = cfl_max * dx/u 

        implicit none 
        
        real(wp), intent(INOUT) :: ux(:,:)
        real(wp), intent(INOUT) :: uy(:,:)
        real(wp), intent(IN)    :: dx, dy
        real(wp), intent(IN)    :: dt 
        real(wp), intent(IN)    :: cfl_max       ! Maximum Courant number, default cfl_max=1.0
        
        ! Local variables 
        integer    :: i, j, q, nx, ny 
        real(wp) :: uxy, dt_now 
        real(wp) :: X, X_max 
        real(wp) :: f_scale 

        nx = size(ux,1)
        ny = size(ux,2)

        X_max = cfl_max / dt 

        do j = 1, ny 
        do i = 1, nx 

            X = max(abs(ux(i,j))/dx + abs(uy(i,j))/dy,1e-5)

            if (X .gt. X_max) then 
                ! Reduce velocity of this point to below limit

                dt_now = cfl_max / X

                f_scale = X_max / X     ! Should be less than 1.0! 

                ux(i,j) = ux(i,j)*f_scale 
                uy(i,j) = uy(i,j)*f_scale 
                
            end if 

        end do 
        end do  


        return 

    end subroutine calc_adv2D_velocity
    
    function calc_adv3D_timestep(ux,uy,uz,H_ice,zeta_ac,dx,dy,cfl_max) result(dt)
        ! Calculate maximum advective time step based
        ! on Courant–Friedrichs–Lewy condition
        ! https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition

        ! 1D condition: C = u*dt/dx <= cfl_max 
        ! 2D condition: C = u*dt/dx + v*dt/dy <= cfl_max 
        ! thus when C = cfl_max:
        ! dt = cfl_max * 1/(u/dx+v/dx)

        ! Note: this is used by Bueler et al. (2007), but it seems 
        ! to impose a very, very small timestep, given the vertical
        ! velocity at the surface (ie, smb) essentially ends up being the
        ! limiting condition. 

        implicit none 
        
        real(wp), intent(IN) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
        real(wp), intent(IN) :: H_ice(:,:) 
        real(wp), intent(IN) :: zeta_ac(:) 
        real(wp), intent(IN) :: dx, dy
        real(wp), intent(IN) :: cfl_max             ! Maximum Courant number, default cfl_max=1.0
        real(wp) :: dt 

        ! Local variables 
        integer    :: i, j, k, nx, ny, nz_aa 
        real(wp) :: dt_check, dt_max 
        real(wp) :: dz, ux_aa, uy_aa, uz_aa  

        real(wp), parameter :: tol = 1e-5 

        nx    = size(ux,1)
        ny    = size(ux,2)
        nz_aa = size(ux,3)

        ! Start with a really high time step 
        dt_max = cfl_max * 1.0 / tol 
        dt     = dt_max

!         write(*,*) "cfl_max = ", cfl_max 
!         write(*,*) "dt_max  = ", dt_max 

!         write(*,*) "calc_adv3D_timestep:: Error: This routine is not working yet."
        
!         stop 

        ! Loop over horizontal grid points 
        do j = 2, ny 
        do i = 2, nx 

            if (H_ice(i,j) .gt. 0.0) then 

                ! Loop over the vertical layers
                do k = 1, nz_aa-1 
                    ux_aa = 0.5*(ux(i-1,j,k)+ux(i,j,k))
                    uy_aa = 0.5*(uy(i,j-1,k)+uy(i,j,k))

                    if (k .le. 1 .or. k .ge. nz_aa-1) then 
                        ! No interpolation of vertical velocity at the base or surface
                        uz_aa = uz(i,j,k) 
                    else
                        ! Interpolation to vertical aa-nodes
                        uz_aa = 0.5*(uz(i,j,k-1)+uz(i,j,k))
                    end if 

                    if (k .le. 1) then 
                        dz = 1e-5 
                    else 
                        dz = max(H_ice(i,j)*(zeta_ac(k)-zeta_ac(k-1)),1e-5)
                    end if 
                    
                    dt_check = cfl_max * 1.0 / max(abs(ux_aa)/dx + abs(uy_aa)/dy + abs(uz_aa)/dz,tol)

!                     write(*,*) i, j, k, dt_check, dt_max, ux_aa, ux_aa, uz_aa, &
!                                     abs(ux_aa)/dx + abs(ux_aa)/dy + abs(uz_aa)/dz

                    dt_max = min(dt_check,dt_max)
                end do 

            end if 

        end do 
        end do 

!         stop 

        return 

    end function calc_adv3D_timestep
    
    subroutine calc_checkerboard(var_check,var,mask)

        implicit none 

        real(wp), intent(OUT) :: var_check(:,:) 
        real(wp), intent(IN)  :: var(:,:) 
        logical,    intent(IN)  :: mask(:,:)  

        ! Local variables 
        integer :: i, j, nx, ny 

        nx = size(var,1)
        ny = size(var,2) 

        ! First assume everything is stable 
        var_check = 0.0_wp 

        do j = 2, ny-1
        do i = 2, nx-1 
            
            if (mask(i,j)) then 
                ! For points of interest, check for checkerboard pattern in var 

                if ( (var(i,j)*var(i-1,j) .lt. 0.0 .and. & 
                      var(i,j)*var(i+1,j) .lt. 0.0) .or. & 
                     (var(i,j)*var(i,j-1) .lt. 0.0 .and. & 
                      var(i,j)*var(i,j+1) .lt. 0.0) ) then 
                    ! Point has checkerboard pattern in at least one direction

                    var_check = var 

                end if 

            end if 

        end do 
        end do  

        return 

    end subroutine calc_checkerboard

    subroutine check_checkerboard(is_unstable,var,lim)

        implicit none 

        logical,    intent(OUT) :: is_unstable
        real(wp), intent(IN)  :: var(:,:) 
        real(wp), intent(IN)  :: lim 

        ! Local variables 
        integer :: i, j, nx, ny 

        nx = size(var,1)
        ny = size(var,2) 

        ! First assume everything is stable 
        is_unstable = .FALSE. 

        do j = 2, ny-1
        do i = 2, nx-1 
 
            if (abs(var(i,j)) .ge. lim) then
                ! Check for checkerboard pattern with var > lim

                if ( (var(i,j)*var(i-1,j) .lt. 0.0 .and. & 
                      var(i,j)*var(i+1,j) .lt. 0.0) .or. & 
                     (var(i,j)*var(i,j-1) .lt. 0.0 .and. & 
                      var(i,j)*var(i,j+1) .lt. 0.0) ) then 
                    ! Point has checkerboard pattern in at least one direction

                    is_unstable = .TRUE. 
                    exit

                end if 

            end if 

        end do 
        end do  

        return 

    end subroutine check_checkerboard

    subroutine yelmo_timestep_write_init(filename,time,xc,yc,pc_eps)

        implicit none 

        character(len=*),  intent(IN) :: filename
        real(wp), intent(IN) :: time
        real(wp), intent(IN) :: xc(:) 
        real(wp), intent(IN) :: yc(:)  
        real(wp), intent(IN) :: pc_eps

        ! Local variables 
        character(len=16) :: xnm 
        character(len=16) :: ynm 
        
        xnm = "xc"
        ynm = "yc" 

        call nc_create(filename)
        
        call nc_write_dim(filename,"pt",x=1,dx=1,nx=1,units="point")

        ! Add grid axis variables to netcdf file
        call nc_write_dim(filename,xnm,x=xc*1e-3,units="kilometers")
        call nc_write_attr(filename,xnm,"_CoordinateAxisType","GeoX")

        call nc_write_dim(filename,ynm,x=yc*1e-3,units="kilometers")
        call nc_write_attr(filename,ynm,"_CoordinateAxisType","GeoY")
        
        call nc_write_dim(filename,"time",x=time,dx=1.0_wp,nx=1,units="years",unlimited=.TRUE.)

        call nc_write(filename, "pc_eps", pc_eps,dim1="pt")
        
        return 

    end subroutine yelmo_timestep_write_init

    subroutine yelmo_timestep_write(filename,time,dt_now,dt_adv,dt_pi,pc_eta,pc_tau, &
                                                speed,speed_tpo,speed_dyn,ssa_iter,iter_redo)

        implicit none 

        character(len=*),  intent(IN) :: filename
        real(wp), intent(IN) :: time 
        real(wp), intent(IN) :: dt_now 
        real(wp), intent(IN) :: dt_adv
        real(wp), intent(IN) :: dt_pi 
        real(wp), intent(IN) :: pc_eta 
        real(wp), intent(IN) :: pc_tau(:,:) 
        real(wp), intent(IN) :: speed 
        real(wp), intent(IN) :: speed_tpo
        real(wp), intent(IN) :: speed_dyn 
        integer,    intent(IN) :: ssa_iter 
        integer,    intent(IN) :: iter_redo 

        ! Local variables
        integer    :: ncid, n, nx, ny 
        real(wp) :: time_prev 

        logical, parameter :: write_pc_tau_field = .FALSE.  ! Signficantly increases filesize, careful!

        nx = size(pc_tau,1)
        ny = size(pc_tau,2) 

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_size(filename,"time",ncid)
        call nc_read(filename,"time",time_prev,start=[n],count=[1],ncid=ncid) 
        if (abs(time-time_prev).gt.1e-5) n = n+1 

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        call nc_write(filename, "speed",speed,dim1="time",start=[n],count=[1],units="kyr/hr",long_name="Yelmo model speed",ncid=ncid)
        call nc_write(filename, "speed_tpo",speed_tpo,dim1="time",start=[n],count=[1],units="kyr/hr",long_name="Yelmo topo speed",ncid=ncid)
        call nc_write(filename, "speed_dyn",speed_dyn,dim1="time",start=[n],count=[1],units="kyr/hr",long_name="Yelmo dyn speed",ncid=ncid)
        
        call nc_write(filename, "dt_now",dt_now,dim1="time",start=[n],count=[1],units="yr",long_name="Timestep",ncid=ncid)
        call nc_write(filename, "dt_adv",dt_adv,dim1="time",start=[n],count=[1],units="yr",long_name="Timestep (CFL criterion)",ncid=ncid)
        
        call nc_write(filename,  "dt_pi", dt_pi,dim1="time",start=[n],count=[1],units="yr",long_name="Timestep (PI controller)",ncid=ncid)
        call nc_write(filename, "pc_eta",pc_eta,dim1="time",start=[n],count=[1],units="m/yr",long_name="eta (maximum PC truncation error)",ncid=ncid)
        
        if (write_pc_tau_field) then 
            call nc_write(filename, "pc_tau",pc_tau,dim1="xc",dim2="yc",dim3="time",start=[1,1,n],count=[nx,ny,1],units="m a**-1", &
                            long_name="Truncation error",ncid=ncid)
        end if 

        call nc_write(filename, "ssa_iter", ssa_iter, dim1="time",start=[n],units="",long_name="Picard iterations for SSA convergence",ncid=ncid)
        call nc_write(filename, "iter_redo",iter_redo,dim1="time",start=[n],units="",long_name="Number of redo iterations needed",count=[1],ncid=ncid)
        
        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine yelmo_timestep_write

    subroutine ytime_init(ytime,nx,ny,nz,dt_min,pc_eps)

        type(ytime_class), intent(INOUT) :: ytime
        integer,    intent(IN) :: nx, ny, nz 
        real(wp), intent(IN) :: dt_min, pc_eps 

        ! Allocate ytime object
        call ytime_alloc(ytime,nx,ny,nz)

        ytime%log_timestep_file = "timesteps.nc" 
        
        ! Initialize information about previous timesteps to minimum
        ! timestep value and high error (to keep low timesteps initially)
        ytime%pc_dt(:)      = dt_min  
        ytime%pc_eta(:)     = pc_eps

        ! Initialize arrays to zero 
        ytime%dt_adv        = 0.0 
        ytime%dt_diff       = 0.0 
        ytime%dt_adv3D      = 0.0 

        ytime%pc_tau        = 0.0_wp 
        ytime%pc_tau_masked = 0.0_wp 
        
        ytime%pc_taus       = 0.0_wp 
        ytime%pc_tau_max    = 0.0_wp
        
        ! Initialize averages to zero too
        ytime%model_speeds  = MV
        ytime%etas          = MV 
        ytime%ssa_iters     = MV 
        ytime%dts           = MV 

        ytime%model_speed   = MV
        ytime%dt_avg        = MV
        ytime%eta_avg       = MV
        ytime%ssa_iter_avg  = MV

        ! Initially, since we don't know pc_dt and pc_eta for
        ! previous timesteps, we should not use this information. 
        ytime%pc_active = .FALSE. 

        return

    end subroutine ytime_init
    
    subroutine ytime_alloc(ytime,nx,ny,nz)

        implicit none 

        type(ytime_class), intent(INOUT) :: ytime
        integer :: nx, ny, nz 

        ! Ensure object is deallocated first
        call ytime_dealloc(ytime) 

        ! Allocate timestep arrays 
        allocate(ytime%dt_adv(nx,ny))
        allocate(ytime%dt_diff(nx,ny))
        allocate(ytime%dt_adv3D(nx,ny,nz))
        
        ! Allocate truncation error array 
        allocate(ytime%pc_tau(nx,ny))
        allocate(ytime%pc_tau_masked(nx,ny))
        
        ! Allocate truncation error averaging arrays 
        allocate(ytime%pc_taus(nx,ny,50))
        allocate(ytime%pc_tau_max(nx,ny))

        return 

    end subroutine ytime_alloc

    subroutine ytime_dealloc(ytime)

        implicit none 

        type(ytime_class), intent(INOUT) :: ytime
        
        if (allocated(ytime%dt_adv))        deallocate(ytime%dt_adv)
        if (allocated(ytime%dt_diff))       deallocate(ytime%dt_diff)
        if (allocated(ytime%dt_adv3D))      deallocate(ytime%dt_adv3D)
        
        if (allocated(ytime%pc_tau))        deallocate(ytime%pc_tau)
        if (allocated(ytime%pc_tau_masked)) deallocate(ytime%pc_tau_masked)
        
        if (allocated(ytime%pc_taus))       deallocate(ytime%pc_taus)
        if (allocated(ytime%pc_tau_max))    deallocate(ytime%pc_tau_max)
        
        return 

    end subroutine ytime_dealloc
    
end module yelmo_timesteps 
