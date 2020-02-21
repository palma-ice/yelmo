module yelmo_timesteps

    use yelmo_defs, only : sp, dp, prec  
    use ncio 

    implicit none 

    private

    public :: set_pc_mask
    public :: calc_pc_eta
    public :: calc_pc_tau_fe_sbe
    public :: calc_pc_tau_ab_sam
    public :: set_adaptive_timestep_pc
    
    public :: set_adaptive_timestep 
    public :: limit_adaptive_timestep

    public :: yelmo_timestep_write_init
    public :: yelmo_timestep_write

    public :: calc_adv2D_timestep1
    public :: calc_adv3D_timestep1 

contains

    subroutine set_pc_mask(mask,H_ice,f_grnd)

        implicit none 

        logical, intent(OUT) :: mask(:,:) 
        real(prec), intent(IN) :: H_ice(:,:) 
        real(prec), intent(IN) :: f_grnd(:,:) 

        ! Local variables 
        integer :: i, j, nx, ny 

        nx = size(mask,1)
        ny = size(mask,2) 

        mask = .FALSE. 

        ! Limit to ice-covered, grounded points 
        where (H_ice .gt. 0.0_prec .and. f_grnd .eq. 1.0_prec) mask = .TRUE. 

        ! Set mask to false for ice margin points as well 
        do j = 2, ny-1 
        do i = 2, nx-1 
            if (mask(i,j)) then 
                if (count(H_ice(i-1:i+1,j-1:j+1).eq.0.0_prec) .gt. 0) then 
                    mask(i,j) = .FALSE.
                end if 
            end if 
        end do 
        end do

        return 

    end subroutine set_pc_mask

    function calc_pc_eta(tau,mask) result(eta)

        implicit none 

        real(prec), intent(IN) :: tau(:,:) 
        logical,    intent(IN) :: mask(:,:) 
        real(prec) :: eta 

        real(prec), parameter :: eta_tol = 1e-8 

        ! Calculate eta 
        eta = maxval(abs(tau),mask=mask)

        ! Limit to non-zero value
        ! Note: Limiting minimum to above eg 1e-8 is very 
        ! important for reducing fluctuations in dt 
        eta = max(eta,eta_tol)

        return 

    end function calc_pc_eta 

    elemental subroutine calc_pc_tau_fe_sbe(tau,var_corr,var_pred,dt_n)
        ! Calculate truncation error for the FE-SBE timestepping method
        ! Forward Euler (FE) predictor step and Semi-implicit
        ! Backward Euler (SBE) corrector step. 
        ! Implemented followig Cheng et al (2017, GMD)
        ! Truncation error: tau = 1/2*dt_n * (var - var_pred)

        implicit none 

        real(prec), intent(OUT) :: tau
        real(prec), intent(IN)  :: var_corr
        real(prec), intent(IN)  :: var_pred
        real(prec), intent(IN)  :: dt_n 
        
        if (dt_n .eq. 0.0_prec) then 
            tau = 0.0_prec 
        else 
            tau = (1.0_prec / (2.0_prec*dt_n)) * (var_corr - var_pred)
        end if 

        return 

    end subroutine calc_pc_tau_fe_sbe

    elemental subroutine calc_pc_tau_ab_sam(tau,var_corr,var_pred,dt_n,zeta)
        ! Calculate truncation error for the AB-SAM timestepping method
        ! Adams-Bashforth (AB) predictor step and Semi-implicit
        ! Adams–Moulton (SAM) corrector step. 
        ! Implemented followig Cheng et al (2017, GMD)
        ! Truncation error: tau = zeta * (var - var_pred) / ((3zeta + 3)*dt)

        implicit none 

        real(prec), intent(OUT) :: tau
        real(prec), intent(IN)  :: var_corr
        real(prec), intent(IN)  :: var_pred
        real(prec), intent(IN)  :: dt_n 
        real(prec), intent(IN)  :: zeta 

        if (dt_n .eq. 0.0_prec) then 
            tau = 0.0_prec 
        else 
            tau = zeta * (var_corr - var_pred) / ((3.0_prec*zeta + 3.0_prec)*dt_n)
        end if 

        return 

    end subroutine calc_pc_tau_ab_sam

    subroutine set_adaptive_timestep_pc(dt_new,dt,eta,eps,dtmin,dtmax,mask,ux_bar,uy_bar,dx,pc_k)
        ! Calculate the timestep following algorithm for 
        ! a general predictor-corrector (pc) method.
        ! Implemented followig Cheng et al (2017, GMD)

        implicit none 

        real(prec), intent(OUT) :: dt_new               ! [yr]   Timestep (n+1)
        real(prec), intent(IN)  :: dt(:)                ! [yr]   Timesteps (n:n-2)
        real(prec), intent(IN)  :: eta(:)               ! [X/yr] Maximum truncation error (n:n-2)
        real(prec), intent(IN)  :: eps                  ! [--]   Tolerance value (eg, eps=1e-4)
        real(prec), intent(IN)  :: dtmin                ! [yr]   Minimum allowed timestep
        real(prec), intent(IN)  :: dtmax                ! [yr]   Maximum allowed timestep
        logical,    intent(IN)  :: mask(:,:)            ! Where to calculate tau 
        real(prec), intent(IN)  :: ux_bar(:,:)          ! [m/yr]
        real(prec), intent(IN)  :: uy_bar(:,:)          ! [m/yr]
        real(prec), intent(IN)  :: dx                   ! [m]
        integer,    intent(IN)  :: pc_k                 ! pc_k gives the order of the timestepping scheme (pc_k=2 for FE-SBE, pc_k=3 for AB-SAM)

        ! Local variables
        real(prec) :: dt_n, dt_nm1, dt_nm2          ! [yr]   Timesteps (n:n-2)
        real(prec) :: eta_n, eta_nm1, eta_nm2       ! [X/yr] Maximum truncation error (n:n-2)
        real(prec) :: rho_n, rho_nm1, rho_nm2
        real(prec) :: rhohat_n 
        real(prec) :: dt_adv 
        real(prec) :: dtmax_now

        real(prec) :: k_i 

        ! Choose adaptive controller algorithm for updating timestep 
        ! PI42, H312b, H312PID
        character(len=56), parameter :: pc_adapt_method = "H312PID"

        ! Smoothing parameter; Söderlind and Wang (2006) method, Eq. 10
        ! Values on the order of [0.7,2.0] are reasonable. Higher kappa slows variation in dt
        real(prec), parameter :: kappa = 2.0_prec 
        
        ! Step 1: Save information needed for adapative controller algorithms 

        ! Save dt from several timesteps
        dt_n    = max(dt(1),dtmin) 
        dt_nm1  = max(dt(2),dtmin) 
        dt_nm2  = max(dt(3),dtmin)

        ! Save eta from several timesteps
        eta_n   = eta(1)
        eta_nm1 = eta(2)
        eta_nm2 = eta(3)

        ! Calculate rho from several timesteps 
        rho_nm1 = (dt_n / dt_nm1) 
        rho_nm2 = (dt_nm1 / dt_nm2) 

        ! Step 2: calculate scaling for the next timestep (dt,n+1)
        select case(trim(pc_adapt_method))

            case("PI42")
                ! Söderlind and Wang, 2006; Cheng et al., 2017
                
                rho_n = calc_pi_rho_pi42(eta_n,eta_nm1,rho_nm1,eps, &
                                            beta_1  =  3.0_prec / (pc_k*5.0_prec),  &
                                            beta_2  = -1.0_prec / (pc_k*5.0_prec), &
                                            alpha_2 =  0.0_prec )


            case("H312b") 
                ! Söderlind (2003) H312b, Eq. 31+ (unlabeled) 
                
                rho_n = calc_pi_rho_H312b(eta_n,eta_nm1,eta_nm2,rho_nm1,rho_nm2,eps,k=real(pc_k,prec),b=8.0_prec)

            case("H312PID") 
                ! Söderlind (2003) H312PD, Eq. 38
                ! Note: Suggested k_i =(2/9)*1/pc_k, but lower value gives more stable solution

                !k_i = (2.0_prec/9.0_prec)*1.0_prec/real(pc_k,prec)
                k_i = 0.08_prec/real(pc_k,prec)

                rho_n = calc_pi_rho_H312PID(eta_n,eta_nm1,eta_nm2,eps,k_i)


            case DEFAULT 

                write(*,*) "set_adaptive_timestep_pc:: Error: pc_adapt_method not recognized."
                write(*,*) "pc_adapt_method = ", trim(pc_adapt_method) 
                stop 

        end select 

        ! Scale rho_n for smoothness 
        rhohat_n = rho_n
        !rhohat_n = min(rho_n,1.1)
        !rhohat_n = 1.0_prec + kappa * atan((rho_n-1.0_prec)/kappa) ! Söderlind and Wang, 2006, Eq. 10
        
        ! Step 3: calculate the next time timestep (dt,n+1)
        dt_new = rhohat_n * dt_n

        ! Step 4: Modify timestep to fit within prescribed limits 

        ! Calculate CFL advection limit too, and limit maximum allowed timestep
        dt_adv    = minval( calc_adv2D_timestep1(ux_bar,uy_bar,dx,dx,cfl_max=1.0_prec) ) 
        dtmax_now = min(dtmax,dt_adv) 

        ! Finally, ensure timestep is within prescribed limits
        call limit_adaptive_timestep(dt_new,dtmin,dtmax_now)

        return 

    end subroutine set_adaptive_timestep_pc

    function calc_pi_rho_pi42(eta_n,eta_nm1,rho_nm1,eps,beta_1,beta_2,alpha_2) result(rho_n)

        implicit none 

        real(prec), intent(IN) :: eta_n 
        real(prec), intent(IN) :: eta_nm1 
        real(prec), intent(IN) :: rho_nm1 
        real(prec), intent(IN) :: eps 
        real(prec), intent(IN) :: beta_1 
        real(prec), intent(IN) :: beta_2
        real(prec), intent(IN) :: alpha_2 
        real(prec) :: rho_n 

        ! Söderlind and Wang, 2006; Cheng et al., 2017
        rho_n   = (eps/eta_n)**beta_1 * (eps/eta_nm1)**beta_2 * rho_nm1**alpha_2

        return 

    end function calc_pi_rho_pi42 

    function calc_pi_rho_H312b(eta_n,eta_nm1,eta_nm2,rho_nm1,rho_nm2,eps,k,b) result(rho_n)

        implicit none 

        real(prec), intent(IN) :: eta_n 
        real(prec), intent(IN) :: eta_nm1 
        real(prec), intent(IN) :: eta_nm2 
        real(prec), intent(IN) :: rho_nm1
        real(prec), intent(IN) :: rho_nm2
        real(prec), intent(IN) :: eps 
        real(prec), intent(IN) :: k 
        real(prec), intent(IN) :: b 
        real(prec) :: rho_n 

        ! Local variables 
        real(prec) :: beta_1, beta_2, beta_3 
        real(prec) :: alpha_2, alpha_3 

        beta_1  =  1.0_prec / (k*b)
        beta_2  =  2.0_prec / (k*b)
        beta_3  =  1.0_prec / (k*b)
        alpha_2 = -3.0_prec / b 
        alpha_3 = -1.0_prec / b 

        ! Söderlind (2003) H312b, Eq. 31+ (unlabeled) 
        rho_n   = (eps/eta_n)**beta_1 * (eps/eta_nm1)**beta_2 * (eps/eta_nm2)**beta_3 &
                            * rho_nm1**alpha_2 * rho_nm2**alpha_3 

        return 

    end function calc_pi_rho_H312b 

    function calc_pi_rho_H312PID(eta_n,eta_nm1,eta_nm2,eps,k_i) result(rho_n)

        implicit none 

        real(prec), intent(IN) :: eta_n 
        real(prec), intent(IN) :: eta_nm1 
        real(prec), intent(IN) :: eta_nm2 
        real(prec), intent(IN) :: eps 
        real(prec), intent(IN) :: k_i  
        real(prec) :: rho_n 

        ! Local variables 
        real(prec) :: k_i_1, k_i_2, k_i_3

        k_i_1   = k_i / 4.0_prec 
        k_i_2   = k_i / 2.0_prec 
        k_i_3   = k_i / 4.0_prec 

        ! Söderlind (2003) H312PID, Eq. 38
        rho_n   = (eps/eta_n)**k_i_1 * (eps/eta_nm1)**k_i_2 * (eps/eta_nm2)**k_i_3

        return 

    end function calc_pi_rho_H312PID 

    subroutine set_adaptive_timestep(dt,dt_adv,dt_diff,dt_adv3D, &
                        ux,uy,uz,ux_bar,uy_bar,D2D,H_ice,dHicedt,zeta_ac, &
                        dx,dtmin,dtmax,cfl_max,cfl_diff_max)
        ! Determine value of adaptive timestep to be consistent with 
        ! min/max timestep range and maximum allowed step of model 
        ! to line up with control time steps

        implicit none 

        real(prec), intent(OUT) :: dt             ! [a] Current timestep 
        real(prec), intent(OUT) :: dt_adv(:,:)     ! [a] Diagnosed maximum advective timestep (vertical ave)
        real(prec), intent(OUT) :: dt_diff(:,:)    ! [a] Diagnosed maximum diffusive timestep (vertical ave) 
        real(prec), intent(OUT) :: dt_adv3D(:,:,:) ! [a] Diagnosed maximum advective timestep (3D) 
        real(prec), intent(IN)  :: ux(:,:,:)       ! [m a-1]
        real(prec), intent(IN)  :: uy(:,:,:)       ! [m a-1]
        real(prec), intent(IN)  :: uz(:,:,:)       ! [m a-1]
        real(prec), intent(IN)  :: ux_bar(:,:)     ! [m a-1]
        real(prec), intent(IN)  :: uy_bar(:,:)     ! [m a-1]
        real(prec), intent(IN)  :: D2D(:,:)        ! [m2 a-1]
        real(prec), intent(IN)  :: H_ice(:,:)      ! [m]
        real(prec), intent(IN)  :: dHicedt(:,:)    ! [m a-1]
        real(prec), intent(IN)  :: zeta_ac(:)      ! [--] 
        real(prec), intent(IN)  :: dx, dtmin, dtmax ! [a]
        real(prec), intent(IN)  :: cfl_max
        real(prec), intent(IN)  :: cfl_diff_max
        
        ! Local variables 
        real(prec) :: dt_adv_min, dt_diff_min 
        real(prec) :: x 
        logical    :: is_unstable
        real(prec), parameter :: dtmax_cfl   = 20.0_prec 
        real(prec), parameter :: exp_cfl     =  2.0_prec 
        real(prec), parameter :: rate_lim    = 1.0_prec   ! Reduction in timestep for instability 
        real(prec), parameter :: rate_scalar = 0.05_prec  ! Reduction in timestep for instability 

        ! Timestep limits determined from CFL conditions for general advective
        ! velocity, as well as diagnosed diffusive magnitude
        ! (adapted from Bueler et al., 2007)

        dt_adv   = calc_adv2D_timestep1(ux_bar,uy_bar,dx,dx,cfl_max)
        dt_diff  = calc_diff2D_timestep(D2D,dx,dx,cfl_diff_max) 

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

        real(prec), intent(INOUT) :: dt               ! [a] Current timestep 
        real(prec), intent(IN)    :: dtmin            ! [a] Minimum allowed timestep
        real(prec), intent(IN)    :: dtmax            ! [a] Maximum allowed timestep 

        ! Local variables  
        real(prec), parameter :: n_decimal   = 4          ! Maximum decimals to treat for timestep
        real(prec), parameter :: dt_half_lim = 0.5_prec   ! Should be 0.5 or greater to make sense

        ! Ensure timestep is also within parameter limits 
        dt = max(dtmin,dt)  ! dt >= dtmin
        dt = min(dtmax,dt)  ! dt <= dtmax

        ! Check to avoid lopsided timesteps (1 big, 1 tiny) to arrive at time_max  
        if (dtmax .gt. 0.0) then 

            if (dt/dtmax .gt. dt_half_lim .and. dt .lt. dtmax) then 
                ! Current adaptive timestep is greater than ~0.5 of the total
                ! expected timestep, and another timestep will be needed to
                ! reach time_max. Therefore, set this timestep to a smaller
                ! value, ie, dt = dt_half_lim*dtmax. 

                dt = dt_half_lim*dtmax

            else if (dt/dtmax .lt. dt_half_lim) then 
                ! Round-off extra digits for neatness

                dt = real(floor(dt*10.0_prec**n_decimal)*10.0_prec**(-n_decimal), prec)
                
            end if 

        end if 

        return 

    end subroutine limit_adaptive_timestep

    elemental function calc_diff2D_timestep(D,dx,dy,cfl_diff_max) result(dt)
        ! Calculate maximum diffusion time step based
        ! on Courant–Friedrichs–Lewy condition
        ! Equation obtained from Bueler et al. (2007), Eq. 25:
        ! dt/2 * (1/dx^2 + 1/dy^2)*max(D) <= cfl_diff_max = 0.12 
        ! dt = cfl_diff_max * 2.0 / ((1/dx^2+1/dy^2)*max(D))

        implicit none 
        
        real(prec), intent(IN) :: D, dx, dy
        real(prec), intent(IN) :: cfl_diff_max       ! Maximum Courant number, default cfl_diff_max=0.12
        real(prec) :: dt 

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
        
        real(prec), intent(IN) :: ux(:,:)           ! acx-nodes
        real(prec), intent(IN) :: uy(:,:)           ! acy-nodes 
        real(prec), intent(IN) :: dx, dy
        real(prec), intent(IN) :: cfl_max           ! Maximum Courant number, default cfl_max=1.0
        real(prec) :: dt(size(ux,1),size(ux,2))     ! aa-nodes 

        ! Local variables  
        integer :: i, j, nx, ny 
        real(prec) :: ux_now, uy_now 

        real(prec), parameter :: eps = 1e-1         ! [m/a] Small factor to avoid divide by zero 

        nx = size(ux,1)
        ny = size(ux,2)

        do j = 2, ny-1 
        do i = 2, nx-1 

!             ux_now = abs( 0.5*(ux(i-1,j)+ux(i,j)) )
!             uy_now = abs( 0.5*(uy(i,j-1)+uy(i,j)) )

            !ux_now = max(abs(ux(i-1,j)),abs(ux(i,j)))
            !uy_now = max(abs(uy(i,j-1)),abs(uy(i,j)))
            
!             ux_now = abs(ux(i,j) - ux(i-1,j))
!             uy_now = abs(uy(i,j) - uy(i,j-1))

!             dt(i,j) = cfl_max * 1.0 / (ux_now/dx + uy_now/dy + eps/dx)

            dt(i,j) = cfl_max * 1.0 / (abs(ux(i-1,j))/dx + abs(ux(i,j))/dx &
                                       + abs(uy(i,j-1))/dy + abs(uy(i,j))/dy + eps/(dx+dy))
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
        
        real(prec), intent(IN) :: ux(:,:,:)        ! acx-nodes
        real(prec), intent(IN) :: uy(:,:,:)        ! acy-nodes
        real(prec), intent(IN) :: uz(:,:,:)        ! acz-nodes  
        real(prec), intent(IN) :: dx, dy
        real(prec), intent(IN) :: H_ice(:,:)       ! aa-nodes 
        real(prec), intent(IN) :: zeta_ac(:)       ! ac-nodes 
        real(prec), intent(IN) :: cfl_max          ! Maximum Courant number, default cfl_max=1.0
        real(prec) :: dt(size(ux,1),size(ux,2),size(ux,3))    ! aa-nodes 

        ! Local variables  
        integer :: i, j, k, nx, ny, nz_aa  
        real(prec) :: ux_now, uy_now 
        real(prec) :: dz 

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
        
        real(prec), intent(IN) :: ux
        real(prec), intent(IN) :: uy
        real(prec), intent(IN) :: dx, dy
        real(prec), intent(IN) :: cfl_max             ! Maximum Courant number, default cfl_max=1.0
        real(prec) :: dt 

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
        
        real(prec), intent(INOUT) :: ux(:,:)
        real(prec), intent(INOUT) :: uy(:,:)
        real(prec), intent(IN)    :: dx, dy
        real(prec), intent(IN)    :: dt 
        real(prec), intent(IN)    :: cfl_max       ! Maximum Courant number, default cfl_max=1.0
        
        ! Local variables 
        integer    :: i, j, q, nx, ny 
        real(prec) :: uxy, dt_now 
        real(prec) :: X, X_max 
        real(prec) :: f_scale 

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
        
        real(prec), intent(IN) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
        real(prec), intent(IN) :: H_ice(:,:) 
        real(prec), intent(IN) :: zeta_ac(:) 
        real(prec), intent(IN) :: dx, dy
        real(prec), intent(IN) :: cfl_max             ! Maximum Courant number, default cfl_max=1.0
        real(prec) :: dt 

        ! Local variables 
        integer    :: i, j, k, nx, ny, nz_aa 
        real(prec) :: dt_check, dt_max 
        real(prec) :: dz, ux_aa, uy_aa, uz_aa  

        real(prec), parameter :: tol = 1e-5 

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

        real(prec), intent(OUT) :: var_check(:,:) 
        real(prec), intent(IN)  :: var(:,:) 
        logical,    intent(IN)  :: mask(:,:)  

        ! Local variables 
        integer :: i, j, nx, ny 

        nx = size(var,1)
        ny = size(var,2) 

        ! First assume everything is stable 
        var_check = 0.0_prec 

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
        real(prec), intent(IN)  :: var(:,:) 
        real(prec), intent(IN)  :: lim 

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
        real(prec), intent(IN) :: time
        real(prec), intent(IN) :: xc(:) 
        real(prec), intent(IN) :: yc(:)  
        real(prec), intent(IN) :: pc_eps

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
        
        call nc_write_dim(filename,"time",x=time,dx=1.0_prec,nx=1,units="years",unlimited=.TRUE.)

        call nc_write(filename, "pc_eps", pc_eps,dim1="pt")
        
        return 

    end subroutine yelmo_timestep_write_init 

    subroutine yelmo_timestep_write(filename,time,dt_now,dt_adv,dt_pi,pc_eta,pc_tau)

        implicit none 

        character(len=*),  intent(IN) :: filename
        real(prec), intent(IN) :: time 
        real(prec), intent(IN) :: dt_now 
        real(prec), intent(IN) :: dt_adv
        real(prec), intent(IN) :: dt_pi 
        real(prec), intent(IN) :: pc_eta 
        real(prec), intent(IN) :: pc_tau(:,:) 

        ! Local variables
        integer    :: ncid, n, nx, ny 
        real(prec) :: time_prev 

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

        call nc_write(filename, "dt_now",dt_now,dim1="time",start=[n],count=[1],ncid=ncid)
        call nc_write(filename, "dt_adv",dt_adv,dim1="time",start=[n],count=[1],ncid=ncid)
        
        call nc_write(filename,  "dt_pi", dt_pi,dim1="time",start=[n],count=[1],ncid=ncid)
        call nc_write(filename, "pc_eta",pc_eta,dim1="time",start=[n],count=[1],ncid=ncid)
        call nc_write(filename, "pc_tau",pc_tau,dim1="xc",dim2="yc",dim3="time",start=[1,1,n],count=[nx,ny,1],ncid=ncid)
        
        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine yelmo_timestep_write 

end module yelmo_timesteps 
