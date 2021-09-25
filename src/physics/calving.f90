module calving
    ! Definitions for various calving laws 

    use yelmo_defs, only : sp, dp, wp, prec, TOL_UNDERFLOW, rho_ice, rho_sw, g  
    use deformation, only : calc_stress_eigen_values

    implicit none 


    private 

    public :: apply_calving     ! ajr: currently not used...
    public :: calc_calving_residual
    public :: calc_calving_rate_kill 
    
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

    subroutine apply_calving(H_ice,calv,f_grnd,H_min_flt,dt)
        ! Given a diagnosed calving rate, make additional modifications
        ! as needed and apply the calving rate to the ice thickness for this timestep

        implicit none 

        real(wp), intent(INOUT) :: H_ice(:,:) 
        real(wp), intent(INOUT) :: calv(:,:) 
        real(wp), intent(IN)    :: f_grnd(:,:) 
        real(wp), intent(IN)    :: H_min_flt
        real(wp), intent(IN)    :: dt  
        
        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: n_mrgn 

        real(wp), parameter :: tau_mrgn = 5.0         ! [a] Time scale for calving of points with many ice-free neighbors
        nx = size(H_ice,1)
        ny = size(H_ice,2) 

        do j = 2, ny-1 
        do i = 2, nx-1

            if (f_grnd(i,j) .eq. 0.0 .and. H_ice(i,j) .gt. 0.0) then 
                ! Floating point, diagnose number of ice-free neighbors 

                n_mrgn = count([H_ice(i-1,j),H_ice(i+1,j),H_ice(i,j-1),H_ice(i,j+1)].eq.0.0 )

            else 

                n_mrgn = 0 

            end if 

            if (n_mrgn .gt. 2) then 
                ! For points with more than two ice-free neighbors, increase calving rate 
                ! (this is designed to handle rare, ice peninsulas that can protrude
                !  from the main ice body)
                
                calv(i,j) = calv(i,j) + max(1000.0-H_ice(i,j),0.0)/tau_mrgn 

            end if 

            ! Additionally modify calving to remove any margin ice less than H_min_flt 
            if (n_mrgn .gt. 0 .and. H_ice(i,j) .lt. H_min_flt) calv(i,j) = H_ice(i,j)/dt
            
            ! Ensure calving is limited to amount of available ice to calve  
            if(f_grnd(i,j) .eq. 0.0 .and. (H_ice(i,j)-dt*calv(i,j)) .lt. 0.0) calv(i,j) = H_ice(i,j)/dt

            ! Apply modified mass balance to update the ice thickness 
            H_ice(i,j) = H_ice(i,j) - dt*calv(i,j)
            
        end do 
        end do 

        ! Also ensure tiny numeric ice thicknesses are removed
        where (f_grnd .eq. 0.0 .and. H_ice .lt. 1e-5) H_ice = 0.0 
            
        return 
        
    end subroutine apply_calving

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
        logical,    intent(IN)  :: mask(:,:) 
        real(wp), intent(IN)  :: tau 
        real(wp), intent(IN)  :: dt 

        if (tau .eq. 0.0_prec) then 
            ! Kill all ice immediately 

            where (mask) calv = H_ice / dt

        else 
            ! Kill using characteristic timescale 
            where (mask) calv = H_ice / tau 

        end if 

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

                ! How many calving-front neighbors?
                wt = 0.0 
                if (f_grnd(im1,j) .eq. 0.0 .and. f_ice(im1,j).eq.0.0) wt = wt + 1.0_wp
                if (f_grnd(ip1,j) .eq. 0.0 .and. f_ice(ip1,j).eq.0.0) wt = wt + 1.0_wp
                if (f_grnd(i,jm1) .eq. 0.0 .and. f_ice(i,jm1).eq.0.0) wt = wt + 1.0_wp
                if (f_grnd(i,jp1) .eq. 0.0 .and. f_ice(i,jp1).eq.0.0) wt = wt + 1.0_wp

                if (wt .gt. 0.0_wp) then 
                    ! Point is at calving front 

                    ! Calculate current ice thickness (H_eff = H_ice/f_ice)
                    ! Check f_ice==0 for safety, but this should never happen for an ice-covered point
                    if (f_ice(i,j) .gt. 0.0_prec) then 
                        H_eff = H_ice(i,j) / f_ice(i,j) 
                    else
                        H_eff = H_ice(i,j) 
                    end if 

                    if (H_eff .lt. H_calv) then 
                        ! If ice is too thin, diagnose calving rate, with
                        ! faster calving timescale for more exposed fronts
                        ! (multiply by f_ice to convert to a horizontal calving rate)
                        
                        ! calv(i,j) = f_ice(i,j) * (H_calv-H_eff) / tau                        
                        !calv(i,j) = f_ice(i,j) * (H_calv-H_eff) * wt / tau
                        calv(i,j) = (H_calv-H_eff) * wt / tau
                        
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
            if (f_ice(i,j) .gt. 0.0) then 
                H_eff = H_ice(i,j) / f_ice(i,j) 
            else 
                H_eff = H_ice(i,j) 
            end if 
            
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
                !calv(i,j) = f_ice(i,j) * max(H_calv - H_eff,0.0_wp) / tau
                calv(i,j) = max(H_calv - H_eff,0.0_wp) / tau
                

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
            if (f_ice(i,j) .gt. 0.0) then 
                H_eff = H_ice(i,j) / f_ice(i,j) 
            else 
                H_eff = H_ice(i,j)  ! == 0.0
            end if 

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

                    if (f_ice(i,j) .gt. 0.0) then 
                        H_eff = H_ice(i,j) / f_ice(i,j) 
                    else 
                        H_eff = H_ice(i,j)  ! == 0.0
                    end if 

                    !calv(i,j) = f_ice(i,j) * max(H_calv - H_eff,0.0) / tau
                    calv(i,j) = max(H_calv - H_eff,0.0) / tau

                end if  

            end if

        end do
        end do

        return 

    end subroutine calc_calving_rate_flux_grisli
    
    subroutine calc_calving_rate_vonmises_l19(calv,H_ice,f_ice,f_grnd,teig1,teig2,ATT_bar,visc_bar,dx,dy,kt,w2,n_glen)
        ! Calculate the 'horizontal' calving rate [m/yr] based on the 
        ! von Mises stress approach, as outlined by Lipscomb et al. (2019)
        ! Eqs. 73-75.
        ! L19: kt = 0.0025 m yr-1 Pa-1, w2=25

        implicit none 

        real(wp), intent(OUT) :: calv(:,:)
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: f_grnd(:,:)  
        real(wp), intent(IN)  :: teig1(:,:)
        real(wp), intent(IN)  :: teig2(:,:)
        real(wp), intent(IN)  :: ATT_bar(:,:)
        real(wp), intent(IN)  :: visc_bar(:,:)
        real(wp), intent(IN)  :: dx, dy 
        real(wp), intent(IN)  :: kt
        real(wp), intent(IN)  :: w2
        real(wp), intent(IN)  :: n_glen 

        ! Local variables 
        integer  :: i, j
        integer  :: im1, jm1, ip1, jp1
        integer  :: nx, ny
        integer  :: n_ocean 
        logical  :: is_margin 
        real(wp) :: tau1, tau2 
        real(wp) :: tau_eff 
        real(wp) :: calv_ref
        real(wp) :: H_eff 

        real(wp) :: ddiv, dxx, dyy, dxy 
        real(wp) :: txx, tyy, txy
        real(wp) :: teig1_now, teig2_now

        real(wp) :: wt
        
        nx = size(H_ice,1)
        ny = size(H_ice,2)

        calv = 0.0_wp

        do j = 1, ny
        do i = 1, nx  
            
            ! Get neighbor indices
            im1 = max(i-1,1) 
            ip1 = min(i+1,nx) 
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny) 
            
            ! Ice-shelf floating margin: floating ice point with open ocean neighbor 
            ! If this point is an ice front, calculate calving
            is_margin =  (f_grnd(i,j) .eq. 0.0 .and. H_ice(i,j) .gt. 0.0) .and. &
                       ( (f_grnd(im1,j) .eq. 0.0 .and. H_ice(im1,j).eq.0.0) .or. &
                         (f_grnd(ip1,j) .eq. 0.0 .and. H_ice(ip1,j).eq.0.0) .or. &
                         (f_grnd(i,jm1) .eq. 0.0 .and. H_ice(i,jm1).eq.0.0) .or. &
                         (f_grnd(i,jp1) .eq. 0.0 .and. H_ice(i,jp1).eq.0.0) ) 

            if (is_margin) then 
                ! Calculate calving rate here 

                ! Calculate the effective stress
                tau_eff = calc_tau_eff(teig1(i,j),teig2(i,j),w2)
                                 
                ! Calculate lateral calving rate 
                calv_ref = kt*tau_eff 

                ! Convert to horizontal volume change 
                if (f_ice(i,j) .gt. 0.0) then 
                    H_eff = H_ice(i,j) / f_ice(i,j)
                else 
                    H_eff = H_ice(i,j)
                end if 

                calv(i,j) = (H_eff*calv_ref) / sqrt(dx*dy)

            end if

        end do
        end do

        return 

    end subroutine calc_calving_rate_vonmises_l19

    elemental function calc_tau_eff(teig1,teig2,w2) result(tau_eff) 

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

    end function calc_tau_eff
        
    subroutine calc_calving_rate_eigen(calv,H_ice,f_grnd,f_ice,ux_bar,uy_bar,dx,dy,H_calv,k_calv)
        ! Calculate the calving rate [m/a] based on the "eigencalving" law
        ! from Levermann et al. (2012)

        ! Note: this routine is untested and incorrect. strain rate must be projected onto 
        ! principal direction of flow. 

        implicit none 

        real(wp), intent(OUT) :: calv(:,:)
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: f_grnd(:,:)  
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: ux_bar(:,:)
        real(wp), intent(IN)  :: uy_bar(:,:)
        real(wp), intent(IN)  :: dx, dy 
        real(wp), intent(IN)  :: H_calv 
        real(wp), intent(IN)  :: k_calv

        ! Local variables 
        integer :: i, j, nx, ny
        real(wp), allocatable :: eps_xx(:,:), eps_yy(:,:) 
        real(wp) :: eps_xx_max, eps_yy_max 
        integer :: imax_xx, jmax_xx, imax_yy, jmax_yy                                          

        logical, allocatable :: is_front(:,:)   
        integer :: n_ocean 

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        allocate(is_front(nx,ny))

        do j = 2, ny
        do i = 2, nx
            ! Calculate strain rate locally (aa-node)
            eps_xx(i,j) = (ux_bar(i,j) - ux_bar(i-1,j))/dx
            eps_yy(i,j) = (uy_bar(i,j) - uy_bar(i,j-1))/dy            
        end do
        end do
        

        calv = 0.0

        do j = 2, ny-1
        do i = 2, nx-1  
            
            if ( (f_grnd(i,j) .eq. 0.0 .and. H_ice(i,j) .gt. 0.0) .and. &
                   ( (f_grnd(i-1,j) .eq. 0.0 .and. H_ice(i-1,j).eq.0.0) .or. &
                     (f_grnd(i+1,j) .eq. 0.0 .and. H_ice(i+1,j).eq.0.0) .or. &
                     (f_grnd(i,j-1) .eq. 0.0 .and. H_ice(i,j-1).eq.0.0) .or. &
                     (f_grnd(i,j+1) .eq. 0.0 .and. H_ice(i,j+1).eq.0.0) ) ) then 
                ! Ice-shelf floating margin: floating ice point with open ocean neighbor 
                ! If this point is an ice front, check for calving

                if ((eps_xx(i,j).gt.0.0).and.(eps_yy(i,j).gt.0.0)) then                   
                    ! Divergence in both directions, apply calving law 
                    ! Flux condition + calving rate with spreading:       

                    calv(i,j) = k_calv * eps_xx(i,j)*eps_yy(i,j)                                       
                
                end if

            end if

        end do
        end do

        return 

    end subroutine calc_calving_rate_eigen

! ===================================================================
!
! Calving - grounded ice 
!
! ===================================================================

    subroutine calc_calving_ground_rate_stress_b12(calv,H_ice,f_ice,f_grnd,H_ocn,tau)
        ! Remove marginal ice that exceeds a stress threshold following
        ! Bassis and Walker (2012), Eq. 2.12 

        implicit none 

        real(wp), intent(OUT) :: calv(:,:)              ! [m/yr] Grounded horizontal calving rate 
        real(wp), intent(IN)  :: H_ice(:,:)             ! [m] Ice thickness 
        real(wp), intent(IN)  :: f_ice(:,:)             ! [--] Ice area fraction 
        real(wp), intent(IN)  :: f_grnd(:,:)            ! [-] Grounded fraction
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

            ! Determine if grounded, ice-covered point has an ice-free neighbor (ie, at the grounded ice margin)
            is_grnd_margin = (f_ice(i,j) .gt. 0.0 .and. f_grnd(i,j) .gt. 0.0 &
                .and. count([f_ice(im1,j),f_ice(ip1,j),f_ice(i,jm1),f_ice(i,jp1)] .eq. 0.0) .gt. 0)

            if (is_grnd_margin) then 
                ! Margin point

                ! Get effective ice thickness 
                if (f_ice(i,j) .gt. 0.0) then 
                    H_eff = H_ice(i,j) / f_ice(i,j)
                else 
                    H_eff = H_ice(i,j)
                end if 

                ! Calculate depth of seawater (limited by ice thickness and flotation criterion)
                if (H_ocn(i,j) .gt. 0.0) then 
                    H_ocn_now = min(rho_ice_sw*H_eff,H_ocn(i,j))
                else 
                    H_ocn_now = 0.0 
                end if 

                ! Get depth-averaged shear-stress in ice, Bassis and Walker (2012), Eq. 2.13 vertically integrated
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

                    ! Calculate calving rate from linear function, 
                    ! limited to available ice thickness 
                    calv(i,j) = min(f_scale*calv_max, H_ice(i,j)/tau) 
                    
                end if 

            end do 
            end do 

        end if 

        return 

    end subroutine calc_calving_ground_rate_stdev

end module calving
