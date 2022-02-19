module calving
    ! Definitions for various calving laws 

    use yelmo_defs, only : sp, dp, wp, prec, TOL_UNDERFLOW, rho_ice, rho_sw, g  
    use topography, only : calc_H_eff 

    implicit none 


    private 

    public :: apply_thin_calving_rate

    public :: calc_calving_tongues
    public :: calc_calving_residual
    public :: calc_calving_rate_kill 
    
    ! Calving related stress/strain routines 
    public :: calc_eps_eff
    public :: calc_tau_eff

    ! Floating calving routines 
    public :: calc_calving_rate_simple
    public :: calc_calving_rate_flux 
    public :: calc_calving_rate_flux_grisli
    public :: calc_calving_rate_vonmises_l19
    public :: calc_calving_rate_eigen
    
    ! Grounded calving routines 
    public :: calc_calving_ground_rate_stress_b12
    public :: calc_calving_ground_rate_stdev

contains 
    
    subroutine apply_thin_calving_rate(calv_flt,H_ice,f_ice,f_grnd,calv_thin)
        ! Adjust calving rate based on ice thickness 
        ! to ensure that thin ice (calv_thin*1yr=Xm) is removed
        ! following Pattyn (2017), Eq. 24. Typical parameters 
        ! calv_thin = 30 m/yr 
        ! H_ref     = 200 m 

        implicit none 

        real(wp), intent(INOUT) :: calv_flt(:,:) 
        real(wp), intent(IN)    :: H_ice(:,:) 
        real(wp), intent(IN)    :: f_ice(:,:) 
        real(wp), intent(IN)    :: f_grnd(:,:) 
        real(wp), intent(IN)    :: calv_thin
        
        ! Local variables
        integer  :: i, j, nx, ny, n_mrgn, n_grnd 
        integer  :: im1, ip1, jm1, jp1
        real(wp) :: H_eff 
        real(wp) :: wt 

        real(wp), parameter :: H_ref = 200.0_wp     ! [m] Thickness below which to scale calving rate

        nx = size(calv_flt,1)
        ny = size(calv_flt,2) 

        do j = 1, ny 
        do i = 1, nx 

            ! Define neighbor indices
            im1 = max(i-1,1)
            ip1 = min(i+1,nx)
            jm1 = max(j-1,1)
            jp1 = min(j+1,ny)
            
            ! Count number of grounded ice-covered neighbors
            n_grnd = count([f_ice(im1,j),f_ice(ip1,j),f_ice(i,jm1),f_ice(i,jp1)].gt.0.0 .and. &
                           [f_grnd(im1,j),f_grnd(ip1,j),f_grnd(i,jm1),f_grnd(i,jp1)].gt.0.0)

            ! Determine if point is at the floating margin with no grounded neighbors
            if (f_grnd(i,j) .eq. 0.0 .and. f_ice(i,j) .gt. 0.0 .and. n_grnd .eq. 0) then 
                ! Floating point, diagnose number of ice-free neighbors 

                n_mrgn = count([f_ice(im1,j),f_ice(ip1,j),f_ice(i,jm1),f_ice(i,jp1)].eq.0.0 )

            else 

                n_mrgn = 0 

            end if 

            if (n_mrgn .gt. 0) then 
                ! Floating ice margin point

                ! Calculate current ice thickness (H_eff = H_ice/f_ice)
                call calc_H_eff(H_eff,H_ice(i,j),f_ice(i,j))

                ! Get weighting factor based on effective ice thickness 
                wt = min(1.0_wp,H_eff/H_ref)

                ! Calculate adjusted calving rate, weighted
                ! between minimum rate and actual value 
                calv_flt(i,j) = calv_thin*(1.0_wp-wt) + calv_flt(i,j)*wt 

            end if 

        end do 
        end do  

        return 

    end subroutine apply_thin_calving_rate

    subroutine calc_calving_tongues(calv_flt,H_ice,f_ice,f_grnd,tau)
        ! Increase calving for floating margin points with 3+ calving
        ! fronts to avoid protruding ice tongues. 

        implicit none 

        real(wp), intent(INOUT) :: calv_flt(:,:) 
        real(wp), intent(IN)    :: H_ice(:,:) 
        real(wp), intent(IN)    :: f_ice(:,:) 
        real(wp), intent(IN)    :: f_grnd(:,:)  
        real(wp), intent(IN)    :: tau 

        ! Local variables 
        integer  :: i, j, nx, ny
        integer  :: im1, ip1, jm1, jp1 
        integer  :: n_mrgn, n_grnd
        real(wp) :: H_eff 
        logical  :: embayed 

        nx = size(H_ice,1)
        ny = size(H_ice,2) 

        do j = 1, ny 
        do i = 1, nx

            ! Define neighbor indices
            im1 = max(i-1,1)
            ip1 = min(i+1,nx)
            jm1 = max(j-1,1)
            jp1 = min(j+1,ny)
            
            ! Count number of grounded ice-covered neighbors
            n_grnd = count([f_ice(im1,j),f_ice(ip1,j),f_ice(i,jm1),f_ice(i,jp1)].gt.0.0 .and. &
                           [f_grnd(im1,j),f_grnd(ip1,j),f_grnd(i,jm1),f_grnd(i,jp1)].gt.0.0)

            if (f_grnd(i,j) .eq. 0.0 .and. f_ice(i,j) .gt. 0.0 .and. n_grnd .eq. 0) then 
                ! Floating point with no grounded neighbors, diagnose number of ice-free neighbors 

                n_mrgn = count([f_ice(im1,j),f_ice(ip1,j),f_ice(i,jm1),f_ice(i,jp1)].eq.0.0 )

            else 

                n_mrgn = 0 

            end if 


            if (n_mrgn .gt. 2) then 
                ! For points with more than two ice-free neighbors, increase calving rate 
                ! (this is designed to handle rare, ice peninsulas that can protrude
                !  from the main ice body)
                
                ! Calculate effective ice thickness for current cell
                if (f_ice(i,j) .gt. 0.0_prec) then 
                    H_eff = H_ice(i,j) / f_ice(i,j) 
                else
                    H_eff = H_ice(i,j) 
                end if 

                calv_flt(i,j) = calv_flt(i,j) + max(1000.0-H_eff,0.0)/tau 

            end if 


            ! Also check for points with an ice-free direct neighbor
            ! but two ice-covered neighbors in the corners. Assume that this
            ! should reduce the calving rate. 
            if (f_grnd(i,j) .eq. 0.0 .and. f_ice(i,j) .eq. 1.0) then 
                ! Floating point with full ice coverage

                embayed = .FALSE.

                ! Embayed to the right?
                if (   (f_grnd(ip1,j)   .eq. 0.0 .and. f_ice(ip1,j)   .eq. 0.0) &
                 .and. (f_grnd(ip1,jm1) .eq. 0.0 .and. f_ice(ip1,jm1) .eq. 1.0) &
                 .and. (f_grnd(ip1,jp1) .eq. 0.0 .and. f_ice(ip1,jp1) .eq. 1.0) ) then 

                    embayed = .TRUE. 

                end if 

                ! Embayed to the left?
                if (   (f_grnd(im1,j)   .eq. 0.0 .and. f_ice(im1,j)   .eq. 0.0) &
                 .and. (f_grnd(im1,jm1) .eq. 0.0 .and. f_ice(im1,jm1) .eq. 1.0) &
                 .and. (f_grnd(im1,jp1) .eq. 0.0 .and. f_ice(im1,jp1) .eq. 1.0) ) then 

                    embayed = .TRUE. 

                end if 

                ! Embayed to the top?
                if (   (f_grnd(i,jp1)   .eq. 0.0 .and. f_ice(i,jp1)   .eq. 0.0) &
                 .and. (f_grnd(im1,jp1) .eq. 0.0 .and. f_ice(im1,jp1) .eq. 1.0) &
                 .and. (f_grnd(ip1,jp1) .eq. 0.0 .and. f_ice(ip1,jp1) .eq. 1.0) ) then 

                    embayed = .TRUE. 

                end if 

                ! Embayed to the bottom?
                if (   (f_grnd(i,jp1)   .eq. 0.0 .and. f_ice(i,jp1)   .eq. 0.0) &
                 .and. (f_grnd(im1,jm1) .eq. 0.0 .and. f_ice(im1,jm1) .eq. 1.0) &
                 .and. (f_grnd(ip1,jm1) .eq. 0.0 .and. f_ice(ip1,jm1) .eq. 1.0) ) then 

                    embayed = .TRUE. 

                end if 


                ! ajr: this code needs testing in realistic setting - not activated yet!

                if (embayed) then 

                    calv_flt(i,j) = 0.5_wp * calv_flt(i,j)
                
                end if 

            end if 

        end do 
        end do 

        return 
        
    end subroutine calc_calving_tongues

    subroutine calc_calving_residual(calv,H_ice,f_ice,dt,resid_lim)

        implicit none 

        real(wp), intent(INOUT) :: calv(:,:) 
        real(wp), intent(IN)    :: H_ice(:,:) 
        real(wp), intent(IN)    :: f_ice(:,:)
        real(wp), intent(IN)    :: dt 
        real(wp), intent(IN), optional :: resid_lim 

        ! Local variables 
        integer :: i, j, nx, ny
        integer :: im1, ip1, jm1, jp1 
        real(wp) :: H_eff
        real(wp) :: wts(4)
        real(wp), allocatable :: calv_resid(:,:) 

        real(wp), parameter :: tol = 1e-3 

        nx = size(calv,1)
        ny = size(calv,2) 

        allocate(calv_resid(nx,ny))

        
        ! Diagnose residual calving 
        where( dt*calv - H_ice .gt. 0.0_wp) 
            calv_resid = (dt*calv-H_ice)/dt
            calv       = H_ice/dt 
        elsewhere
            calv_resid = 0.0_wp
        end where

        ! Limit residual calving as desired 
        if (present(resid_lim)) then 

            ! For upstream cells, limit calving to a fraction
            ! of the available ice thickness (resid_lim*H_ice)
            ! rather than the diagnosed calv_resid value. This 
            ! is useful for threshold-based calving. Usually 
            ! a limit of just 0.01 is imposed, so that the cell 
            ! becomes less than fully ice covered and dynamically 
            ! inactive. Approach designed following CISM in 
            ! glissade_calving

            where(calv_resid .gt. 0.0_wp .and. calv_resid .gt. (H_ice*resid_lim)/dt)
                calv_resid = (H_ice*resid_lim)/dt
            end where 

        end if 

        ! Avoid underflow errors
        where (abs(calv) .lt. TOL_UNDERFLOW) calv = 0.0_wp 

        ! Determine if any residual calving should migrate to neighboring inland
        ! cells if all ice in the margin cell is deleted.
        do j = 1, ny 
        do i = 1, nx 

            ! Define neighbor indices
            im1 = max(i-1,1)
            ip1 = min(i+1,nx)
            jm1 = max(j-1,1)
            jp1 = min(j+1,ny)
            
            ! Calculate effective ice thickness for current cell
            if (f_ice(i,j) .gt. 0.0_prec) then 
                H_eff = H_ice(i,j) / f_ice(i,j) 
            else
                H_eff = H_ice(i,j) 
            end if 

            if (calv_resid(i,j) .gt. 0.0_wp) then 
                ! Calving diagnosed for this point 

                wts = 0.0_wp 

                if (f_ice(ip1,j) .eq. 1.0 .and. &
                    (abs(H_ice(ip1,j)-H_eff) .lt. tol .or. H_ice(ip1,j) .lt. H_eff)) wts(1) = 1.0_wp 
                if (f_ice(i,jp1) .eq. 1.0 .and. &
                    (abs(H_ice(i,jp1)-H_eff) .lt. tol .or. H_ice(i,jp1) .lt. H_eff)) wts(2) = 1.0_wp 
                if (f_ice(im1,j) .eq. 1.0 .and. &
                    (abs(H_ice(im1,j)-H_eff) .lt. tol .or. H_ice(im1,j) .lt. H_eff)) wts(3) = 1.0_wp 
                if (f_ice(i,jm1) .eq. 1.0 .and. &
                    (abs(H_ice(i,jm1)-H_eff) .lt. tol .or. H_ice(i,jm1) .lt. H_eff)) wts(4) = 1.0_wp 
                
                if (sum(wts) .gt. 0.0_wp) then 
                    ! Normalize available weights to sum to 1
                    wts = wts / sum(wts) 
                end if 

                calv(ip1,j) = calv(ip1,j) + calv_resid(i,j)*wts(1)
                calv(i,jp1) = calv(i,jp1) + calv_resid(i,j)*wts(2)
                calv(im1,j) = calv(im1,j) + calv_resid(i,j)*wts(3)
                calv(i,jm1) = calv(i,jm1) + calv_resid(i,j)*wts(4)

                calv_resid(i,j) = calv_resid(i,j) - sum(wts)*calv_resid(i,j)
            
            end if 
        end do 
        end do 

if (.FALSE.) then 
    ! ajr: this check is not needed, generally everything sums up well.

        if (maxval(abs(calv_resid)) .gt. 1e-3) then 
            write(*,*) "calc_calving_residual:: Error: residual calving not &
            & properly accounted for."
            write(*,*) "calv_resid: ", minval(calv_resid), maxval(calv_resid)
            stop 
        end if 

end if 

        return 

    end subroutine calc_calving_residual

    subroutine calc_calving_rate_kill(calv,H_ice,mask,tau,dt)
        ! Kill all ice in a given mask using a characteristic timescale tau

        implicit none 

        real(wp), intent(OUT) :: calv(:,:)
        real(wp), intent(IN)  :: H_ice(:,:)
        logical,  intent(IN)  :: mask(:,:) 
        real(wp), intent(IN)  :: tau 
        real(wp), intent(IN)  :: dt 

        ! Local variables
        real(wp) :: dt_kill 

        ! Make sure dt is not zero
        dt_kill = dt 
        if (dt_kill .eq. 0.0) dt_kill = 1.0_wp 

        ! Kill all ice immediately 
        ! Ensure all ice calves by imposing a higher rate than ice exists

        where (mask) calv = H_ice / dt_kill * 1.1

        ! ajr: in principle, we could make use of a timescale as below,
        ! however, for most 'kill' applications, this added complexity is
        ! not helpful (ie, shelves might not be fully killed when they are
        ! expected to be). This needs further development, so far now,
        ! the lines above are active where all ice is calved immediately.
        
        ! if (tau .eq. 0.0_wp) then 
        !     ! Kill all ice immediately 

        !     where (mask) calv = H_ice / dt

        ! else 
        !     ! Kill using characteristic timescale 
        !     where (mask) calv = H_ice / tau 

        ! end if 

        return 

    end subroutine calc_calving_rate_kill

! ===================================================================
!
! Calving - floating ice 
!
! ===================================================================

    subroutine calc_calving_rate_simple(calv,H_ice,f_ice,f_grnd,H_calv,tau)
        ! Calculate the calving rate [m/a] based on a simple threshold rule
        ! H_ice < H_calv

        implicit none 

        real(wp), intent(OUT) :: calv(:,:)
        real(wp), intent(IN)  :: H_ice(:,:)                 ! [m] Ice thickness 
        real(wp), intent(IN)  :: f_ice(:,:)                 ! [--] Ice area fraction
        real(wp), intent(IN)  :: f_grnd(:,:)                ! [--] Grounded fraction
        real(wp), intent(IN)  :: H_calv                     ! [m] Calving thickness threshold
        real(wp), intent(IN)  :: tau                        ! [a] Calving timescale, ~ 1yr

        ! Local variables 
        integer  :: i, j, nx, ny
        integer  :: im1, ip1, jm1, jp1
        real(wp) :: H_eff
        integer  :: wt 

        nx = size(H_ice,1)
        ny = size(H_ice,2)
            
            
        ! Initially set calving rate to zero 
        calv = 0.0 

        do j=1,ny
        do i=1,nx

            ! Get neighbor indices 
            im1 = max(i-1,1)
            ip1 = min(i+1,nx)
            jm1 = max(j-1,1)
            jp1 = min(j+1,ny)

            if ( (f_grnd(i,j) .eq. 0.0 .and. f_ice(i,j) .gt. 0.0) ) then 
                ! Floating ice point

                ! How many calving-fronts for this cell?
                wt = 0.0 
                if (f_grnd(im1,j) .eq. 0.0 .and. f_ice(im1,j).eq.0.0) wt = wt + 1.0_wp
                if (f_grnd(ip1,j) .eq. 0.0 .and. f_ice(ip1,j).eq.0.0) wt = wt + 1.0_wp
                if (f_grnd(i,jm1) .eq. 0.0 .and. f_ice(i,jm1).eq.0.0) wt = wt + 1.0_wp
                if (f_grnd(i,jp1) .eq. 0.0 .and. f_ice(i,jp1).eq.0.0) wt = wt + 1.0_wp

                if (wt .gt. 0.0_wp) then 
                    ! Point is at calving front 

                    ! Calculate current ice thickness (H_eff = H_ice/f_ice)
                    call calc_H_eff(H_eff,H_ice(i,j),f_ice(i,j))

                    if (H_eff .lt. H_calv) then 
                        ! If ice is too thin, diagnose calving rate.
                        ! Scale by f_ice to apply to whole cell (following Lipscomb et al., 2019)
                        
                        calv(i,j) = f_ice(i,j) * ( (H_calv-H_eff) / tau )
                        
                    end if 

                end if

            end if

        end do
        end do

        return 

    end subroutine calc_calving_rate_simple
    
    subroutine calc_calving_rate_flux(calv,H_ice,f_ice,f_grnd,mbal,ux,uy,dx,H_calv,tau)
        ! Calculate the calving rate [m/a] based on a simple threshold rule
        ! H_ice < H_calv

        implicit none 

        real(wp), intent(OUT) :: calv(:,:)                  ! [m/yr] Calving rate scaled to horizontal grid size
        real(wp), intent(IN)  :: H_ice(:,:)                 ! [m] Ice thickness 
        real(wp), intent(IN)  :: f_ice(:,:)                 ! [--] Ice area fraction
        real(wp), intent(IN)  :: f_grnd(:,:)                ! [-] Grounded fraction
        real(wp), intent(IN)  :: mbal(:,:)                  ! [m/yr] Net mass balance 
        real(wp), intent(IN)  :: ux(:,:)                    ! [m/yr] velocity, x-direction (ac-nodes)
        real(wp), intent(IN)  :: uy(:,:)                    ! [m/yr] velocity, y-direction (ac-nodes)
        real(wp), intent(IN)  :: dx                         ! [m] Grid resolution
        real(wp), intent(IN)  :: H_calv                     ! [m] Threshold for calving
        real(wp), intent(IN)  :: tau                        ! [yr] Calving timescale, ~ 1yr

        ! Local variables 
        integer  :: i, j, nx, ny
        integer  :: im1, ip1, jm1, jp1
        real(wp) :: dy  
        real(wp) :: flux_xr, flux_xl, flux_yd, flux_yu
        real(wp) :: dhdt_upstream
        logical  :: positive_mb 
        real(wp), allocatable :: dHdt(:,:)
        real(wp), allocatable :: H_diff(:,:)   
        real(wp) :: H_eff
        real(wp) :: wt 

        dy = dx 

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        allocate(dHdt(nx,ny))
        allocate(H_diff(nx,ny))

        ! Calculate lagrangian rate of change and thickness relative to threshold

        dHdt   = 0.0 
        H_diff = 0.0 

        do j = 1, ny
        do i = 1, nx
            
            ! Get neighbor indices
            im1 = max(i-1,1)
            ip1 = min(i+1,nx)
            jm1 = max(j-1,1)
            jp1 = min(j+1,ny)

            ! Get effective ice thickness
            call calc_H_eff(H_eff,H_ice(i,j),f_ice(i,j))

            ! Calculate ice thickness relative to the calving threshold 
            H_diff(i,j) = H_eff - H_calv


            ! Calculate expected ice thickness change at each ice covered point

            if (f_ice(i,j) .eq. 1.0) then 
                ! Fully ice-covered points     
                
                ! Use explicit upstream advection algorithm 
                ! to calculate the flux across each boundary [m^2 a^-1]
                if (ux(i,j) .gt. 0.0) then 
                    flux_xr = ux(i,j)*H_ice(i,j)
                else
                    flux_xr = ux(i,j)*H_ice(ip1,j) 
                end if 

                if (ux(im1,j) .gt. 0.0) then 
                    flux_xl = ux(im1,j)*H_ice(im1,j)
                else
                    flux_xl = ux(im1,j)*H_ice(i,j) 
                end if 

                if (uy(i,j) .gt. 0.0) then 
                    flux_yu = uy(i,j)*H_ice(i,j)
                else
                    flux_yu = uy(i,j)*H_ice(i,jp1) 
                end if 

                if (uy(i,jm1) .gt. 0.0) then 
                    flux_yd = uy(i,jp1)*H_ice(i,jm1)
                else
                    flux_yd = uy(i,jp1)*H_ice(i,j) 
                end if 

                ! Calculate flux divergence on aa-nodes combined with mass balance
                ! (ie,thickness change via conservation)
                dHdt(i,j) = f_ice(i,j)*mbal(i,j) &
                    + (1.0 / dx) * (flux_xl - flux_xr) + (1.0 / dy) * (flux_yd - flux_yu)

            end if 
                 
        end do 
        end do
        
        ! Initialize calving to zero 
        calv = 0.0_wp 

        do j = 1, ny
        do i = 1, nx

            ! Get neighbor indices 
            im1 = max(i-1,1)
            ip1 = min(i+1,nx)
            jm1 = max(j-1,1)
            jp1 = min(j+1,ny)
            
            ! Check if mass balance of current point is positive, 
            ! and there is at least one grounded neighbor. In this
            ! case, do not allow calving.
            positive_mb = mbal(i,j).gt.0.0 .and. &
                (f_grnd(im1,j).gt.0.0 .or. f_grnd(ip1,j).gt.0.0 .or. &
                 f_grnd(i,jm1).gt.0.0 .or. f_grnd(i,jp1).gt.0.0)
            
            if ( (f_grnd(i,j) .eq. 0.0 .and. f_ice(i,j) .gt. 0.0 .and. &
                    H_diff(i,j).lt.0.0) .and. (.not. positive_mb) .and. &
                   ( (f_grnd(im1,j) .eq. 0.0 .and. f_ice(im1,j).eq.0.0) .or. &
                     (f_grnd(ip1,j) .eq. 0.0 .and. f_ice(ip1,j).eq.0.0) .or. &
                     (f_grnd(i,jm1) .eq. 0.0 .and. f_ice(i,jm1).eq.0.0) .or. &
                     (f_grnd(i,jp1) .eq. 0.0 .and. f_ice(i,jp1).eq.0.0) ) ) then 
                ! Ice-shelf floating margin with potential for calving: 
                ! floating ice point with open ocean neighbor, 
                ! that is also below the calving threshold thickness,
                ! and not with positive mass balance next to grounded ice. 

                ! Get the flux [m^2/yr] from upstream points that are 
                ! fully ice covered, with ice thickness greater than calving threshold 
                ! and velocity flowing into the margin point of interest. 

                if (f_ice(im1,j).eq.1.0.and.H_diff(im1,j).gt.0.0.and.ux(im1,j).gt.0.0) then
                    !flux_xl = max( (dHdt(im1,j)+H_diff(im1,j))/tau,0.0_wp) 
                    flux_xl = ux(im1,j)*H_diff(im1,j)
                end if 

                if (f_ice(ip1,j).eq.1.0.and.H_diff(ip1,j).gt.0.0.and.ux(i,j).lt.0.0) then 
                    !flux_xr = max( (dHdt(ip1,j)+H_diff(ip1,j))/tau,0.0_wp)
                    flux_xr = ux(i,j)*H_diff(ip1,j)
                end if

                if (f_ice(i,jm1).eq.1.0.and.H_diff(i,jm1).gt.0.0.and.uy(i,jm1).gt.0.0) then
                    !flux_yd = max( (dHdt(i,jm1)+H_diff(i,jm1))/tau,0.0_wp) 
                    flux_yd = uy(i,jm1)*H_diff(i,jm1)
                end if

                if (f_ice(i,jp1).eq.1.0.and.H_diff(i,jp1).gt.0.0.and.uy(i,j).lt.0.0) then
                    !flux_yu = max( (dHdt(i,jp1)+H_diff(i,jp1))/tau,0.0_wp)
                    flux_yu = uy(i,j)*H_diff(i,jp1)
                end if

                ! Diagnose total incoming horizontal flux to calving cell [m/yr]
                dhdt_upstream = (flux_xl+flux_xr)/dx + (flux_yd+flux_yu)/dy 


                ! Calculate effective ice thickness of calving cell
                if (f_ice(i,j) .gt. 0.0) then 
                    H_eff = H_ice(i,j) / f_ice(i,j) 
                else 
                    H_eff = H_ice(i,j)  ! == 0.0
                end if 


                ! Diagnose calving rate not accounting for upstream flux
                ! Scale by f_ice to apply to whole cell (following Lipscomb et al., 2019)
                calv(i,j) = f_ice(i,j) * (max(H_calv - H_eff,0.0_wp) / tau)
                

                ! Adjust calving rate for upstream flux 
                calv(i,j) = max(calv(i,j) - dhdt_upstream,0.0_wp)

                if (abs(calv(i,j)) .lt. TOL_UNDERFLOW) calv(i,j) = 0.0_wp 
            end if

        end do
        end do

        return 

    end subroutine calc_calving_rate_flux
    
    subroutine calc_calving_rate_flux_grisli(calv,H_ice,f_ice,f_grnd,mbal,ux,uy,dx,H_calv,tau)
        ! Calculate the calving rate [m/a] based on a simple threshold rule
        ! H_ice < H_calv

        implicit none 

        real(wp), intent(OUT) :: calv(:,:)                  ! [m/yr] Calving rate scaled to horizontal grid size
        real(wp), intent(IN)  :: H_ice(:,:)                 ! [m] Ice thickness 
        real(wp), intent(IN)  :: f_ice(:,:)                 ! [--] Ice area fraction
        real(wp), intent(IN)  :: f_grnd(:,:)                ! [-] Grounded fraction
        real(wp), intent(IN)  :: mbal(:,:)                  ! [m/yr] Net mass balance 
        real(wp), intent(IN)  :: ux(:,:)                    ! [m/yr] velocity, x-direction (ac-nodes)
        real(wp), intent(IN)  :: uy(:,:)                    ! [m/yr] velocity, y-direction (ac-nodes)
        real(wp), intent(IN)  :: dx                         ! [m] Grid resolution
        real(wp), intent(IN)  :: H_calv                     ! [m] Threshold for calving
        real(wp), intent(IN)  :: tau                        ! [yr] Calving timescale, ~ 1yr

        ! Local variables 
        integer :: i, j, nx, ny
        integer :: im1, ip1, jm1, jp1
        real(wp) :: dxx, dyy  
        logical :: flux_mij, flux_pij, flux_imj, flux_ipj
        logical :: positive_mb 
        real(wp), allocatable :: dHdt(:,:)
        real(wp), allocatable :: H_diff(:,:)  
        real(wp), allocatable :: ddiv(:,:)  
        real(wp) :: H_eff
        real(wp) :: wt 

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        allocate(dHdt(nx,ny))
        allocate(H_diff(nx,ny))
        allocate(ddiv(nx,ny)) 

        ! Calculate horizontal divergence over the grid
        do j = 1, ny
        do i = 1, nx
            
            ! Get neighbor indices
            im1 = max(i-1,1)
            jm1 = max(j-1,1)
            
            if (f_ice(i,j) .eq. 1.0) then 

                ! Calculate horizontal divergence 
                ddiv(i,j) = (ux(i,j)-ux(im1,j))/dx &
                          + (uy(i,j)-uy(i,jm1))/dx

                ! Avoid underflow errors 
                if (abs(ddiv(i,j)) .lt. 1e-8) ddiv(i,j) = 0.0_prec
            
            else 

                ddiv(i,j) = 0.0_wp

            end if 

        end do 
        end do
        
        ! Extrapolate to partially ice-covered points, then 
        ! calculate lagrangian rate of change and thickness relative to threshold

        dHdt   = 0.0 
        H_diff = 0.0 

        do j = 1, ny
        do i = 1, nx
            
            ! Get neighbor indices
            im1 = max(i-1,1)
            ip1 = min(i+1,nx)
            jm1 = max(j-1,1)
            jp1 = min(j+1,ny)

            if ( f_ice(i,j) .lt. 1.0 .and. &
                count([f_ice(im1,j),f_ice(ip1,j),f_ice(i,jm1),f_ice(i,jp1)] .eq. 1.0_wp) .gt. 0 ) then 
                ! Ice-free (or partially ice-free) with ice-covered neighbors

                ddiv(i,j) = 0.0
                wt        = 0.0 

                if (f_ice(im1,j) .eq. 1.0) then 
                    ddiv(i,j) = ddiv(i,j) + ddiv(im1,j)
                    wt = wt + 1.0 
                end if 
                if (f_ice(ip1,j) .eq. 1.0) then 
                    ddiv(i,j) = ddiv(i,j) + ddiv(ip1,j)
                    wt = wt + 1.0 
                end if
                if (f_ice(i,jm1) .eq. 1.0) then 
                    ddiv(i,j) = ddiv(i,j) + ddiv(i,jm1)
                    wt = wt + 1.0 
                end if
                if (f_ice(i,jp1) .eq. 1.0) then 
                    ddiv(i,j) = ddiv(i,j) + ddiv(i,jp1)
                    wt = wt + 1.0 
                end if

                if (wt .gt. 0.0) then 
                    ddiv(i,j) = ddiv(i,j) / wt
                end if 

            end if 

            ! Margin points 
            call calc_H_eff(H_eff,H_ice(i,j),f_ice(i,j))

            ! Calculate thickness change via conservation
            dHdt(i,j) = mbal(i,j) - H_eff*ddiv(i,j)

            ! Also calculate ice thickness relative to the calving threshold 
            H_diff(i,j) = H_eff - H_calv 

        end do 
        end do
        
        ! Initially set calving rate to zero 
        calv = 0.0 

        do j = 1, ny
        do i = 1, nx

            ! Get neighbor indices 
            im1 = max(i-1,1)
            ip1 = min(i+1,nx)
            jm1 = max(j-1,1)
            jp1 = min(j+1,ny)
            
            if ( (f_grnd(i,j) .eq. 0.0 .and. f_ice(i,j) .gt. 0.0 .and. H_diff(i,j).lt.0.0) .and. &
                   ( (f_grnd(im1,j) .eq. 0.0 .and. f_ice(im1,j).eq.0.0) .or. &
                     (f_grnd(ip1,j) .eq. 0.0 .and. f_ice(ip1,j).eq.0.0) .or. &
                     (f_grnd(i,jm1) .eq. 0.0 .and. f_ice(i,jm1).eq.0.0) .or. &
                     (f_grnd(i,jp1) .eq. 0.0 .and. f_ice(i,jp1).eq.0.0) ) ) then 
                ! Ice-shelf floating margin: floating ice point with open ocean neighbor 
                ! that is also below the calving threshold thickness 

                ! Check if current point would still be below threshold when 
                ! accounting for mass flux from inland (ie, inland ice 
                ! thickness is larger than threshold and flow is going into
                ! the margin point, and inland point has enough thickness to 
                ! maintain itself above the threshold, even accounting for ice 
                ! advected out of it). Or, inland point is grounded with positive mass balance

                positive_mb = (mbal(i,j).gt.0.0)

                flux_mij = ( ((H_diff(im1,j).gt.0.0).and.(ux(im1,j).gt.0.0)  &  ! neighbor (im1,j) total > H_calv
                    .and.  (dHdt(im1,j).gt.(-H_diff(im1,j)*abs(ux(im1,j)/dx)))) & 
                    .or.(f_grnd(im1,j).gt.0.0.and.positive_mb ))

                flux_pij = ( ((H_diff(ip1,j).gt.0.0).and.(ux(i,j).lt.0.0) & ! neighbor (ip1,j) total > H_calv
                    .and.(dHdt(ip1,j).gt.(-H_diff(ip1,j)*abs(ux(i,j)/dx)))) &
                    .or.(f_grnd(ip1,j).gt.0.0.and.positive_mb ))

                flux_imj = ( ((H_diff(i,jm1).gt.0.0).and.(uy(i,jm1).gt.0.0)  &  ! neighbor (i,jm1) total > H_calv
                    .and.(dHdt(i,jm1).gt.(-H_diff(i,jm1)*abs(uy(i,jm1)/dx))))&
                    .or.(f_grnd(i,jm1).gt.0.0.and.positive_mb ))

                flux_ipj = ( ((H_diff(i,jp1).gt.0.0).and.(uy(i,j).lt.0.0) & ! neighbor (i,jp1) total > H_calv
                    .and.(dHdt(i,jp1).gt.(-H_diff(i,jp1)*abs(uy(i,j)/dx))))&
                    .or.(f_grnd(i,jp1).gt.0.0.and.positive_mb ))

                if ( .not. (flux_mij.or.flux_pij.or.flux_imj.or.flux_ipj) ) then
                    ! This point is not sustained from neighbors, determine calving rate 
                    ! f_ice ensures rate is adjusted to size of grid cell 

                    ! Get effective ice thickness
                    call calc_H_eff(H_eff,H_ice(i,j),f_ice(i,j))

                    ! Calculate calving rate.
                    ! Scale by f_ice to apply to whole cell (following Lipscomb et al., 2019)
                    calv(i,j) = f_ice(i,j) * (max(H_calv - H_eff,0.0) / tau)

                end if  

            end if

        end do
        end do

        return 

    end subroutine calc_calving_rate_flux_grisli
    
    subroutine calc_calving_rate_vonmises_l19(calv,H_ice,f_ice,f_grnd,tau_eff,dx,kt)
        ! Calculate the 'horizontal' calving rate [m/yr] based on the 
        ! von Mises stress approach, as outlined by Lipscomb et al. (2019)
        ! Eqs. 73-75.
        ! L19: kt = 0.0025 m yr-1 Pa-1, w2=25

        implicit none 

        real(wp), intent(OUT) :: calv(:,:)
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: f_grnd(:,:)  
        real(wp), intent(IN)  :: tau_eff(:,:)
        real(wp), intent(IN)  :: dx
        real(wp), intent(IN)  :: kt

        ! Local variables 
        integer  :: i, j
        integer  :: im1, jm1, ip1, jp1
        integer  :: nx, ny 
        real(wp) :: dy   
        real(wp) :: calv_ref
        real(wp) :: H_eff 
        real(wp) :: wt
        
        real(wp), parameter :: calv_lim = 2000.0_wp     ! To avoid really high calving values

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        ! Assume square grid cells 
        dy = dx 

        calv = 0.0_wp

        do j = 1, ny
        do i = 1, nx  
            
            ! Get neighbor indices
            im1 = max(i-1,1) 
            ip1 = min(i+1,nx) 
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny) 
            
            if ( (f_grnd(i,j) .eq. 0.0 .and. f_ice(i,j) .gt. 0.0) ) then 
                ! Floating ice point, calculate calving rate at the margin 

                ! How many calving-fronts for this cell?
                wt = 0.0 
                if (f_grnd(im1,j) .eq. 0.0 .and. f_ice(im1,j).eq.0.0) wt = wt + 1.0_wp
                if (f_grnd(ip1,j) .eq. 0.0 .and. f_ice(ip1,j).eq.0.0) wt = wt + 1.0_wp
                if (f_grnd(i,jm1) .eq. 0.0 .and. f_ice(i,jm1).eq.0.0) wt = wt + 1.0_wp
                if (f_grnd(i,jp1) .eq. 0.0 .and. f_ice(i,jp1).eq.0.0) wt = wt + 1.0_wp

                if (wt .gt. 0.0_wp) then 
                    ! Point is at calving front 

                    ! Calculate lateral calving rate 
                    calv_ref = kt*tau_eff(i,j) 

                    ! Get effective ice thickness
                    call calc_H_eff(H_eff,H_ice(i,j),f_ice(i,j))

                    ! Convert to horizontal volume change, weighted
                    ! by number of exposed faces.
                    calv(i,j) = (H_eff*calv_ref) / sqrt(dx*dy)

                    ! write(*,*) "calv", i, j, tau_eff(i,j), calv_ref,  &
                    !             H_eff, H_eff/sqrt(dx*dy), calv(i,j) 
                    
                    ! Apply calving limit
                    calv(i,j) = min(calv(i,j),calv_lim)
                    
                end if 

            end if

        end do
        end do

        return 

    end subroutine calc_calving_rate_vonmises_l19
       
    subroutine calc_calving_rate_eigen(calv,H_ice,f_ice,f_grnd,eps_eff,dx,k2)
        ! Calculate the 'horizontal' calving rate [m/yr] based on the 
        ! von Mises stress approach, as outlined by Lipscomb et al. (2019)
        ! Eqs. 73-75.
        ! L19: kt = 0.0025 m yr-1 Pa-1, w2=25

        implicit none 

        real(wp), intent(OUT) :: calv(:,:)
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: f_grnd(:,:)  
        real(wp), intent(IN)  :: eps_eff(:,:)
        real(wp), intent(IN)  :: dx
        real(wp), intent(IN)  :: k2

        ! Local variables 
        integer  :: i, j
        integer  :: im1, jm1, ip1, jp1
        integer  :: nx, ny
        integer  :: n_ocean 
        logical  :: is_margin 
        real(wp) :: dy   
        real(wp) :: calv_ref
        real(wp) :: H_eff 

        real(wp) :: dxx, dyy, dxy 
        real(wp) :: eps_eig_1_now, eps_eig_2_now
        real(wp) :: eps_eff_neighb(4)
        real(wp) :: wt
        
        real(wp), parameter :: calv_lim = 2000.0_wp     ! To avoid really high calving values

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        ! Assume square grid cells 
        dy = dx 

        calv = 0.0_wp

        do j = 1, ny
        do i = 1, nx  
            
            ! Get neighbor indices
            im1 = max(i-1,1) 
            ip1 = min(i+1,nx) 
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny) 
            
            if ( (f_grnd(i,j) .eq. 0.0 .and. f_ice(i,j) .gt. 0.0) ) then 
                ! Floating ice point, calculate calving rate at the margin 

                ! How many calving-fronts for this cell?
                wt = 0.0 
                if (f_grnd(im1,j) .eq. 0.0 .and. f_ice(im1,j).eq.0.0) wt = wt + 1.0_wp
                if (f_grnd(ip1,j) .eq. 0.0 .and. f_ice(ip1,j).eq.0.0) wt = wt + 1.0_wp
                if (f_grnd(i,jm1) .eq. 0.0 .and. f_ice(i,jm1).eq.0.0) wt = wt + 1.0_wp
                if (f_grnd(i,jp1) .eq. 0.0 .and. f_ice(i,jp1).eq.0.0) wt = wt + 1.0_wp

                if (wt .gt. 0.0_wp) then 
                    ! Point is at calving front 

                    ! Calculate lateral calving rate 
                    calv_ref = max( k2*eps_eff(i,j), 0.0_wp )

                    ! Get effective ice thickness
                    call calc_H_eff(H_eff,H_ice(i,j),f_ice(i,j))

                    ! Convert to horizontal volume change, weighted
                    ! by number of exposed faces.
                    calv(i,j) = (H_eff*calv_ref) / sqrt(dx*dy)

                    ! write(*,*) "calv", tau_eff, calv_ref,  &
                    !             H_eff, wt*H_eff/sqrt(dx*dy), calv(i,j) 
                    
                    ! Apply calving limit
                    calv(i,j) = min(calv(i,j),calv_lim)

                end if 

            end if

        end do
        end do

        return 

    end subroutine calc_calving_rate_eigen

    subroutine calc_eps_eff(eps_eff,eps_eig_1,eps_eig_2,f_ice)

        implicit none 

        real(wp), intent(OUT) :: eps_eff(:,:)
        real(wp), intent(IN)  :: eps_eig_1(:,:)
        real(wp), intent(IN)  :: eps_eig_2(:,:) 
        real(wp), intent(IN)  :: f_ice(:,:) 

        ! Local variables 
        integer  :: i, j, nx, ny, n  
        integer  :: im1, jm1, ip1, jp1 
        real(wp) :: eps_eff_neighb(4)

        nx = size(eps_eff,1)
        ny = size(eps_eff,2) 

        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices
            im1 = max(i-1,1) 
            ip1 = min(i+1,nx) 
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny) 
            
            if (f_ice(i,j) .eq. 0.0_wp) then 
                ! Ice-free point, no strain

                eps_eff(i,j) = 0.0_wp 

            else if (eps_eig_1(i,j) .eq. 0.0 .and. eps_eig_2(i,j) .eq. 0.0) then 
                ! Margin point was likely just advected, no stresses available, 
                ! use maximum value of eps_eff from upstream neighbors.

                eps_eff_neighb = 0.0_wp 

                if (f_ice(im1,j).gt.0.0) eps_eff_neighb(1) = eps_eig_1(im1,j) * eps_eig_2(im1,j)
                if (f_ice(ip1,j).gt.0.0) eps_eff_neighb(2) = eps_eig_1(ip1,j) * eps_eig_2(ip1,j)
                if (f_ice(i,jm1).gt.0.0) eps_eff_neighb(3) = eps_eig_1(i,jm1) * eps_eig_2(i,jm1)
                if (f_ice(i,jp1).gt.0.0) eps_eff_neighb(4) = eps_eig_1(i,jp1) * eps_eig_2(i,jp1)

                n = count(eps_eff_neighb.ne.0.0_wp)

                if (n .gt. 0) then 
                    eps_eff(i,j) = sum(eps_eff_neighb,mask=eps_eff_neighb.ne.0.0_wp) / real(n,wp)
                else 
                    eps_eff(i,j) = 0.0_wp
                end if 

            else 
                ! Stresses are available at this margin point. 
                ! Calculate the effective strain rate directly.

                eps_eff(i,j) = eps_eig_1(i,j) * eps_eig_2(i,j)
                
            end if 
            
        end do 
        end do 

        return 

    end subroutine calc_eps_eff
    
    subroutine calc_tau_eff(tau_eff,tau_eig_1,tau_eig_2,f_ice,w2)

        implicit none 

        real(wp), intent(OUT) :: tau_eff(:,:)
        real(wp), intent(IN)  :: tau_eig_1(:,:)
        real(wp), intent(IN)  :: tau_eig_2(:,:) 
        real(wp), intent(IN)  :: f_ice(:,:) 
        real(wp), intent(IN)  :: w2 

        ! Local variables 
        integer  :: i, j, nx, ny, n  
        integer  :: im1, jm1, ip1, jp1 
        real(wp) :: tau_eff_neighb(4) 

        nx = size(tau_eff,1)
        ny = size(tau_eff,2) 
        
        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices
            im1 = max(i-1,1) 
            ip1 = min(i+1,nx) 
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny) 
            
            if (f_ice(i,j) .eq. 0.0_wp) then 
                ! Ice-free point, no stress

                tau_eff(i,j) = 0.0_wp 

            else if (tau_eig_1(i,j) .eq. 0.0 .and. tau_eig_2(i,j) .eq. 0.0) then 
                ! Margin point was likely just advected, no stresses available, 
                ! use maximum value of tau_eff from upstream neighbors.

                tau_eff_neighb = 0.0_wp 

                if (f_ice(im1,j).gt.0.0) tau_eff_neighb(1) = calc_tau_eff_now(tau_eig_1(im1,j),tau_eig_2(im1,j),w2)
                if (f_ice(ip1,j).gt.0.0) tau_eff_neighb(2) = calc_tau_eff_now(tau_eig_1(ip1,j),tau_eig_2(ip1,j),w2)
                if (f_ice(i,jm1).gt.0.0) tau_eff_neighb(3) = calc_tau_eff_now(tau_eig_1(i,jm1),tau_eig_2(i,jm1),w2)
                if (f_ice(i,jp1).gt.0.0) tau_eff_neighb(4) = calc_tau_eff_now(tau_eig_1(i,jp1),tau_eig_2(i,jp1),w2)

                n = count(tau_eff_neighb.ne.0.0_wp)

                if (n .gt. 0) then 
                    tau_eff(i,j) = sum(tau_eff_neighb,mask=tau_eff_neighb.ne.0.0_wp) / real(n,wp)
                else 
                    tau_eff(i,j) = 0.0_wp 
                end if 

            else 
                ! Stresses are available at this margin point. 
                ! Calculate the effective strain rate directly.

                tau_eff(i,j) = calc_tau_eff_now(tau_eig_1(i,j),tau_eig_2(i,j),w2)
            
            end if 

        end do 
        end do 

        return 

    end subroutine calc_tau_eff
    
    elemental function calc_tau_eff_now(teig1,teig2,w2) result(tau_eff) 

        implicit none 

        real(wp), intent(IN) :: teig1 
        real(wp), intent(IN) :: teig2
        real(wp), intent(IN) :: w2
        real(wp) :: tau_eff

        ! Local variables 
        real(wp) :: tau1, tau2

        tau1    = max(teig1,0.0_wp)
        tau2    = max(teig2,0.0_wp)
        tau_eff = sqrt(tau1**2 + (w2 * tau2)**2)

        return 

    end function calc_tau_eff_now
     
! ===================================================================
!
! Calving - grounded ice 
!
! ===================================================================

    subroutine calc_calving_ground_rate_stress_b12(calv,H_ice,f_ice,f_grnd,z_bed,H_ocn,tau)
        ! Remove marginal ice that exceeds a stress threshold following
        ! Bassis and Walker (2012), Eq. 2.12 

        implicit none 

        real(wp), intent(OUT) :: calv(:,:)              ! [m/yr] Grounded horizontal calving rate 
        real(wp), intent(IN)  :: H_ice(:,:)             ! [m] Ice thickness 
        real(wp), intent(IN)  :: f_ice(:,:)             ! [--] Ice area fraction 
        real(wp), intent(IN)  :: f_grnd(:,:)            ! [-] Grounded fraction
        real(wp), intent(IN)  :: z_bed(:,:)             ! [m] Bedrock elevation 
        real(wp), intent(IN)  :: H_ocn(:,:)             ! [m] Ocean thickness (depth)
        real(wp), intent(IN)  :: tau                    ! [yr] Calving timescale 

        ! Local variables 
        integer  :: i, j, nx, ny
        integer  :: im1, ip1, jm1, jp1 
        real(wp) :: tau_c, H_eff, H_max, H_ocn_now 
        logical  :: is_grnd_margin  
        real(wp) :: rho_ice_g, rho_sw_ice, rho_ice_sw  

        real(wp), parameter :: C0    = 1e6                ! [Pa] Depth-averaged shear stress in ice 
        real(wp), parameter :: alpha = 0.0                ! [--] Friction coefficient for Bassis and Walker (2012), Eq. 2.13
        real(wp), parameter :: r     = 0.0                ! [--] Crevasse fraction 
        
        logical  :: mask_neighb(4) 

        rho_ice_g  = rho_ice * g 
        rho_sw_ice = rho_sw / rho_ice 
        rho_ice_sw = rho_ice / rho_sw 

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        ! Intialize calving rate to zero 
        calv = 0.0 

        do j = 1, ny
        do i = 1, nx 

            ! Get neighbor indices
            im1 = max(i-1,1) 
            ip1 = min(i+1,nx) 
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny) 

            ! Check if neighbors are ice free and not at higher bedrock elevation
            mask_neighb = ([f_ice(im1,j),f_ice(ip1,j),f_ice(i,jm1),f_ice(i,jp1)] .eq. 0.0) .and. &
                          ([z_bed(im1,j),z_bed(ip1,j),z_bed(i,jm1),z_bed(i,jp1)] .le. z_bed(i,j))

            ! Determine if grounded, ice-covered point has an ice-free neighbor (ie, at the grounded ice margin)
            is_grnd_margin = (f_ice(i,j) .gt. 0.0 .and. f_grnd(i,j) .gt. 0.0 &
                                                  .and. count(mask_neighb) .gt. 0)

            if (is_grnd_margin) then 
                ! Margin point with potential for grounded stress calving

                ! Get effective ice thickness 
                call calc_H_eff(H_eff,H_ice(i,j),f_ice(i,j))

                ! Calculate depth of seawater (limited by ice thickness and flotation criterion)
                if (H_ocn(i,j) .gt. 0.0) then 
                    H_ocn_now = min(rho_ice_sw*H_eff,H_ocn(i,j))
                else 
                    ! Bedrock is above sea level, ocean depth is zero
                    H_ocn_now = 0.0 
                end if 

                ! Get depth-averaged shear-stress in ice, Bassis and Walker (2012), Eq. 2.13 vertically averaged
                ! alpha = 0.65: model S1 validated for cold ice 
                ! alpha = 0.4 : model S2 for warmer ice 
                ! alpha = 0.0 : model S3 for purely plastic yielding (default)
                tau_c = C0 + 0.5*alpha*rho_ice_g*H_eff

                ! Get critical ice thickness to cause stress failure
                H_max = (1.0-r)*tau_c/rho_ice_g + sqrt(((1.0-r)*tau_c/rho_ice_g)**2 + rho_sw_ice*H_ocn_now**2)

                if (H_eff .gt. H_max) then 
                    ! Critical stress exceeded, determine horizontal calving rate 

                    calv(i,j) = f_ice(i,j) * max(H_eff-H_max,0.0) / tau

                end if 

            end if

        end do 
        end do 

        return 

    end subroutine calc_calving_ground_rate_stress_b12
    
    subroutine calc_calving_ground_rate_stdev(calv,H_ice,f_ice,f_grnd,z_bed_sd,sd_min,sd_max,calv_max,tau)
        ! Parameterize grounded ice-margin calving as a function of 
        ! standard deviation of bedrock at each grid point.
        ! Assumes that higher variability in subgrid implies cliffs
        ! that are not represented at low resolution. 

        implicit none 

        real(wp), intent(OUT) :: calv(:,:)                ! [m/yr] Calculated calving rate 
        real(wp), intent(IN)  :: H_ice(:,:)               ! [m] Ice thickness 
        real(wp), intent(IN)  :: f_ice(:,:)               ! [-] Ice area fraction
        real(wp), intent(IN)  :: f_grnd(:,:)              ! [-] Grounded fraction
        real(wp), intent(IN)  :: z_bed_sd(:,:)            ! [m] Standard deviation of bedrock topography
        real(wp), intent(IN)  :: sd_min                   ! [m] stdev(z_bed) at/below which calv=0
        real(wp), intent(IN)  :: sd_max                   ! [m] stdev(z_bed) at/above which calv=calv_max 
        real(wp), intent(IN)  :: calv_max                 ! [m/yr] Maximum allowed calving rate
        real(wp), intent(IN)  :: tau                      ! [yr] Calving timescale       

        ! Local variables
        integer  :: i, j, nx, ny  
        integer  :: im1, ip1, jm1, jp1 
        real(wp) :: f_scale 
        real(wp) :: H_eff
        logical  :: is_grnd_margin 

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        ! Intialize calving rate to zero 
        calv = 0.0 
        
        if (calv_max .gt. 0.0) then 
            ! Determine grounded calving rate 

            do j = 1, ny
            do i = 1, nx 

                ! Get neighbor indices
                im1 = max(i-1,1) 
                ip1 = min(i+1,nx) 
                jm1 = max(j-1,1) 
                jp1 = min(j+1,ny) 

                ! Determine if grounded, ice-covered point has an ice-free neighbor (ie, at the grounded ice margin)
                is_grnd_margin = (f_ice(i,j) .gt. 0.0 .and. f_grnd(i,j) .gt. 0.0 &
                    .and. count([f_ice(im1,j),f_ice(ip1,j),f_ice(i,jm1),f_ice(i,jp1)] .eq. 0.0) .gt. 0)

                if (is_grnd_margin) then 
                    ! Margin point

                    f_scale = (z_bed_sd(i,j) - sd_min)/(sd_max-sd_min)
                    if (f_scale .lt. 0.0) f_scale = 0.0 
                    if (f_scale .gt. 1.0) f_scale = 1.0 

                    ! Get effective ice thickness 
                    call calc_H_eff(H_eff,H_ice(i,j),f_ice(i,j))

                    ! Calculate calving rate from linear function, 
                    ! limited to available (effective) ice thickness 
                    calv(i,j) = min(f_scale*calv_max, H_eff/tau) 
                    
                end if 

            end do 
            end do 

        end if 

        return 

    end subroutine calc_calving_ground_rate_stdev

end module calving
