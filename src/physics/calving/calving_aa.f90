module calving_aa
    ! Definitions for various calving laws on aa nodes. 
    ! This will be used for mass balance calving 

    use yelmo_defs, only : sp, dp, wp, prec, TOL_UNDERFLOW
    use yelmo_tools, only : get_neighbor_indices
    use topography, only : calc_H_eff 

    implicit none 


    private 

    public :: apply_calving_rate_thin

    public :: calc_calving_rate_tongues
    public :: calc_calving_rate_kill 
    
    ! Calving related stress/strain routines 
    public :: calc_eps_eff
    public :: calc_tau_eff

    ! Floating calving routines 
    public :: define_calving_thickness_threshold
    public :: calc_calving_rate_threshold

    public :: define_calving_stress_factor
    public :: calc_calving_rate_vonmises_l19

    public :: calc_calving_rate_eigen
    
    ! Grounded calving routines 
    public :: calc_calving_ground_rate_stress_b12
    public :: calc_calving_ground_rate_stdev

contains 
    
    subroutine apply_calving_rate_thin(mb_calv,H_ice,f_ice,f_grnd,calv_thin,Hc_ref_thin,boundaries)
        ! Adjust calving rate based on ice thickness 
        ! to ensure that thin ice (calv_thin*1yr=Xm) is removed
        ! following Pattyn (2017), Eq. 24. Typical parameters 
        ! calv_thin = 30 m/yr 
        ! H_ref     = 200 m 

        implicit none 

        real(wp), intent(INOUT) :: mb_calv(:,:) 
        real(wp), intent(IN)    :: H_ice(:,:) 
        real(wp), intent(IN)    :: f_ice(:,:) 
        real(wp), intent(IN)    :: f_grnd(:,:) 
        real(wp), intent(IN)    :: calv_thin
        real(wp), intent(IN)    :: Hc_ref_thin      ! [m] Thickness below which to scale calving rate
        character(len=*), intent(IN) :: boundaries 

        ! Local variables
        integer  :: i, j, nx, ny, n_mrgn, n_grnd 
        integer  :: im1, ip1, jm1, jp1
        real(wp) :: H_eff 
        real(wp) :: wt 
        
        nx = size(mb_calv,1)
        ny = size(mb_calv,2) 

        !!$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,n_grnd,n_mrgn,H_eff,wt)
        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

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
                wt = min(1.0_wp,H_eff/Hc_ref_thin)

                ! Calculate adjusted calving rate, weighted
                ! between minimum rate and actual value 
                ! -calv_thin since mb_calv has calving as negative
                mb_calv(i,j) = -calv_thin*(1.0_wp-wt) + mb_calv(i,j)*wt 

            end if 

        end do 
        end do  
        !!$omp end parallel do

        return 

    end subroutine apply_calving_rate_thin

    subroutine calc_calving_rate_tongues(mb_calv,H_ice,f_ice,f_grnd,tau,boundaries)
        ! Increase calving for floating margin points with 3+ calving
        ! fronts to avoid protruding ice tongues. 

        implicit none 

        real(wp), intent(INOUT) :: mb_calv(:,:) 
        real(wp), intent(IN)    :: H_ice(:,:) 
        real(wp), intent(IN)    :: f_ice(:,:) 
        real(wp), intent(IN)    :: f_grnd(:,:)  
        real(wp), intent(IN)    :: tau 
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer  :: i, j, nx, ny
        integer  :: im1, ip1, jm1, jp1 
        integer  :: n_mrgn, n_grnd
        real(wp) :: H_eff 
        logical  :: embayed 

        nx = size(H_ice,1)
        ny = size(H_ice,2) 

        !!$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,n_grnd,n_mrgn,H_eff,embayed)
        do j = 1, ny 
        do i = 1, nx

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

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

                mb_calv(i,j) = mb_calv(i,j) - max(1000.0-H_eff,0.0)/tau 

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
                if (   (f_grnd(i,jm1)   .eq. 0.0 .and. f_ice(i,jm1)   .eq. 0.0) &
                 .and. (f_grnd(im1,jm1) .eq. 0.0 .and. f_ice(im1,jm1) .eq. 1.0) &
                 .and. (f_grnd(ip1,jm1) .eq. 0.0 .and. f_ice(ip1,jm1) .eq. 1.0) ) then 

                    embayed = .TRUE. 

                end if 


                ! ajr: this code needs testing in realistic setting - not activated yet!

                if (embayed) then 

                    mb_calv(i,j) = 0.5_wp * mb_calv(i,j)
                
                end if 

            end if 

        end do 
        end do
        !!$omp end parallel do

        return 
        
    end subroutine calc_calving_rate_tongues

    subroutine calc_calving_rate_kill(mb_calv,H_ice,mask,tau,dt)
        ! Kill all ice in a given mask using a characteristic timescale tau

        implicit none 

        real(wp), intent(OUT) :: mb_calv(:,:)
        real(wp), intent(IN)  :: H_ice(:,:)
        logical,  intent(IN)  :: mask(:,:) 
        real(wp), intent(IN)  :: tau 
        real(wp), intent(IN)  :: dt 

        ! Local variables
        real(wp) :: dt_kill 

        ! Make sure dt is not zero
        dt_kill = dt 
        if (dt_kill .eq. 0.0) dt_kill = 1.0_wp 
        
        ! Reset calving field to zero
        mb_calv = 0.0

        ! Kill all ice immediately 
        ! Ensure all ice calves by imposing a higher rate than ice exists

        where (mask) mb_calv = - (H_ice / dt_kill * 1.1)

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
    
    elemental subroutine define_calving_thickness_threshold(H_calv,z_bed_corr,H_calv_shallow,H_calv_deep,z_shallow,z_deep)

        implicit none

        real(wp), intent(OUT) :: H_calv
        real(wp), intent(IN)  :: z_bed_corr      ! [m] Bed elevation relative to sea level height (negative value)
        real(wp), intent(IN)  :: H_calv_shallow
        real(wp), intent(IN)  :: H_calv_deep
        real(wp), intent(IN)  :: z_shallow 
        real(wp), intent(IN)  :: z_deep

        ! Local variables
        integer :: i, j, nx, ny 

        if(z_bed_corr .le. z_deep) then
            ! Deep water
            H_calv = H_calv_deep
        else if(z_bed_corr .ge. z_shallow) then
            ! Shallow water
            H_calv = H_calv_shallow
        else 
            ! Linear interpolation between thickness thresholds
            H_calv = H_calv_deep + &
                    (H_calv_shallow - H_calv_deep) * (z_bed_corr - z_deep)/(z_shallow - z_deep)
        end if

        return
        
    end subroutine define_calving_thickness_threshold

    subroutine calc_calving_rate_threshold(mb_calv,H_ice,f_ice,f_grnd,H_calv,tau,boundaries)
        ! Calculate the calving rate [m/a] based on a simple threshold rule
        ! H_ice < H_calv

        implicit none 

        real(wp), intent(OUT) :: mb_calv(:,:)
        real(wp), intent(IN)  :: H_ice(:,:)                 ! [m] Ice thickness 
        real(wp), intent(IN)  :: f_ice(:,:)                 ! [--] Ice area fraction
        real(wp), intent(IN)  :: f_grnd(:,:)                ! [--] Grounded fraction
        real(wp), intent(IN)  :: H_calv(:,:)                ! [m] Calving thickness threshold
        real(wp), intent(IN)  :: tau                        ! [a] Calving timescale, ~ 1yr
        character(len=*), intent(IN) :: boundaries

        ! Local variables 
        integer  :: i, j, nx, ny
        integer  :: im1, ip1, jm1, jp1
        real(wp) :: H_eff
        integer  :: wt 

        nx = size(H_ice,1)
        ny = size(H_ice,2)
            
            
        ! Initially set calving rate to zero 
        mb_calv = 0.0 

        !!$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,wt,H_eff)
        do j=1,ny
        do i=1,nx

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
            
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

                    if (H_eff .lt. H_calv(i,j)) then 
                        ! If ice is too thin, diagnose calving rate.
                        ! Scale by f_ice to apply to whole cell (following Lipscomb et al., 2019)
                        
                        mb_calv(i,j) = - ( f_ice(i,j) * ( (H_calv(i,j)-H_eff) / tau ) )
                        
                    end if 

                end if

            end if

        end do
        end do
        !!$omp end parallel do

        return 

    end subroutine calc_calving_rate_threshold
    
    elemental subroutine define_calving_stress_factor(kt,z_bed_corr,kt_shallow,kt_deep,z_shallow,z_deep)

        implicit none

        real(wp), intent(OUT) :: kt
        real(wp), intent(IN)  :: z_bed_corr      ! [m] Bed elevation relative to sea level height (negative value)
        real(wp), intent(IN)  :: kt_shallow
        real(wp), intent(IN)  :: kt_deep
        real(wp), intent(IN)  :: z_shallow 
        real(wp), intent(IN)  :: z_deep

        ! Local variables
        integer :: i, j, nx, ny 

        if(z_bed_corr .le. z_deep) then
            ! Deep water
            kt = kt_deep
        else if(z_bed_corr .ge. z_shallow) then
            ! Shallow water
            kt = kt_shallow
        else 
            ! Linear interpolation between thickness thresholds
            kt = kt_deep + &
                    (kt_shallow - kt_deep) * (z_bed_corr - z_deep)/(z_shallow - z_deep)
        end if

        return
        
    end subroutine define_calving_stress_factor

    subroutine calc_calving_rate_vonmises_l19(mb_calv,H_ice,f_ice,f_grnd,tau_eff,dx,kt,boundaries)
        ! Calculate the 'horizontal' calving rate [m/yr] based on the 
        ! von Mises stress approach, as outlined by Lipscomb et al. (2019)
        ! Eqs. 73-75.
        ! L19: kt = 0.0025 m yr-1 Pa-1, w2=25

        implicit none 

        real(wp), intent(OUT) :: mb_calv(:,:)
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: f_grnd(:,:)  
        real(wp), intent(IN)  :: tau_eff(:,:)
        real(wp), intent(IN)  :: dx
        real(wp), intent(IN)  :: kt(:,:)
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer  :: i, j
        integer  :: im1, jm1, ip1, jp1
        integer  :: nx, ny 
        real(wp) :: dy   
        real(wp) :: calv_ref
        real(wp) :: calv_now
        real(wp) :: H_eff 
        real(wp) :: wt
        
        real(wp), parameter :: calv_lim = 2000.0_wp     ! To avoid really high calving values

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        ! Assume square grid cells 
        dy = dx 

        mb_calv = 0.0_wp

        !!$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,wt,calv_ref,H_eff,calv_now)
        do j = 1, ny
        do i = 1, nx  
            
            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
            
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
                    calv_ref = kt(i,j)*tau_eff(i,j) 

                    ! Get effective ice thickness
                    call calc_H_eff(H_eff,H_ice(i,j),f_ice(i,j))

                    ! Convert to vertical height change spread over the whole cell
                    ! Volume calved, V = H_eff*calv_ref*dx
                    ! Area of grid cell, A = dx*dy
                    ! Height change spread over cell, 
                    !   calv = V/A = H_eff*calv_ref*dx / (dx*dy)
                    !   calv = H_eff*calv_ref / (dy)
                    ! Using sqrt(dx*dy) in the denominator accounts for the case
                    ! where dx != dy. If dx=dy, then expresion below is equivalent to above.
                    calv_now = (H_eff*calv_ref) / sqrt(dx*dy)

                    ! write(*,*) "calv", i, j, tau_eff(i,j), calv_ref,  &
                    !             H_eff, H_eff/sqrt(dx*dy), calv_now 
                    
                    ! Apply calving limit
                    calv_now = min(calv_now,calv_lim)
                    
                    ! Get calving mass balance rate
                    mb_calv(i,j) = -calv_now

                end if 

            end if

        end do
        end do
        !!$omp end parallel do

        return 

    end subroutine calc_calving_rate_vonmises_l19
       
    subroutine calc_calving_rate_eigen(mb_calv,H_ice,f_ice,f_grnd,eps_eff,dx,k2,boundaries)
        ! Calculate the 'horizontal' calving rate [m/yr] based on the 
        ! von Mises stress approach, as outlined by Lipscomb et al. (2019)
        ! Eqs. 73-75.
        ! L19: kt = 0.0025 m yr-1 Pa-1, w2=25

        implicit none 

        real(wp), intent(OUT) :: mb_calv(:,:)
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: f_grnd(:,:)  
        real(wp), intent(IN)  :: eps_eff(:,:)
        real(wp), intent(IN)  :: dx
        real(wp), intent(IN)  :: k2
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer  :: i, j
        integer  :: im1, jm1, ip1, jp1
        integer  :: nx, ny
        integer  :: n_ocean 
        logical  :: is_margin 
        real(wp) :: dy   
        real(wp) :: calv_ref
        real(wp) :: calv_now
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

        mb_calv = 0.0_wp

        !!$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,wt,calv_ref,H_eff,calv_now)
        do j = 1, ny
        do i = 1, nx  
            
            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
            
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
                    calv_now = (H_eff*calv_ref) / sqrt(dx*dy)

                    ! write(*,*) "calv", tau_eff, calv_ref,  &
                    !             H_eff, wt*H_eff/sqrt(dx*dy), calv_now
                    
                    ! Apply calving limit
                    calv_now = min(calv_now,calv_lim)

                    ! Get calving mass balance rate
                    mb_calv(i,j) = -calv_now

                end if 

            end if

        end do
        end do
        !!$omp end parallel do

        return 

    end subroutine calc_calving_rate_eigen

    subroutine calc_eps_eff(eps_eff,eps_eig_1,eps_eig_2,f_ice,boundaries)

        implicit none 

        real(wp), intent(OUT) :: eps_eff(:,:)
        real(wp), intent(IN)  :: eps_eig_1(:,:)
        real(wp), intent(IN)  :: eps_eig_2(:,:) 
        real(wp), intent(IN)  :: f_ice(:,:) 
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer  :: i, j, nx, ny, n  
        integer  :: im1, jm1, ip1, jp1 
        real(wp) :: eps_eff_neighb(4)

        nx = size(eps_eff,1)
        ny = size(eps_eff,2) 

        !!$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,n,eps_eff_neighb)
        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
            
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
        !!$omp end parallel do

        return 

    end subroutine calc_eps_eff
    
    subroutine calc_tau_eff(tau_eff,tau_eig_1,tau_eig_2,f_ice,w2,boundaries)

        implicit none 

        real(wp), intent(OUT) :: tau_eff(:,:)
        real(wp), intent(IN)  :: tau_eig_1(:,:)
        real(wp), intent(IN)  :: tau_eig_2(:,:) 
        real(wp), intent(IN)  :: f_ice(:,:) 
        real(wp), intent(IN)  :: w2 
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer  :: i, j, nx, ny, n  
        integer  :: im1, jm1, ip1, jp1 
        real(wp) :: tau_eff_neighb(4) 

        nx = size(tau_eff,1)
        ny = size(tau_eff,2) 
        
        !!$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,n,tau_eff_neighb)
        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
            
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
        !!$omp end parallel do

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

    subroutine calc_calving_ground_rate_stress_b12(mb_calv,H_ice,f_ice,f_grnd,z_bed,H_ocn,tau,rho_ice,rho_sw,g,boundaries)
        ! Remove marginal ice that exceeds a stress threshold following
        ! Bassis and Walker (2012), Eq. 2.12 

        implicit none 

        real(wp), intent(OUT) :: mb_calv(:,:)           ! [m/yr] Grounded calving mass balance
        real(wp), intent(IN)  :: H_ice(:,:)             ! [m] Ice thickness 
        real(wp), intent(IN)  :: f_ice(:,:)             ! [--] Ice area fraction 
        real(wp), intent(IN)  :: f_grnd(:,:)            ! [-] Grounded fraction
        real(wp), intent(IN)  :: z_bed(:,:)             ! [m] Bedrock elevation 
        real(wp), intent(IN)  :: H_ocn(:,:)             ! [m] Ocean thickness (depth)
        real(wp), intent(IN)  :: tau                    ! [yr] Calving timescale 
        real(wp), intent(IN)  :: rho_ice 
        real(wp), intent(IN)  :: rho_sw
        real(wp), intent(IN)  :: g
        character(len=*), intent(IN) :: boundaries 

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
        mb_calv = 0.0 

        !!$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,mask_neighb,is_grnd_margin,H_eff,H_ocn_now,tau_c,H_max)
        do j = 1, ny
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
            
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
                    ! Critical stress exceeded, determine mass balance calving rate 

                    mb_calv(i,j) = - ( f_ice(i,j) * max(H_eff-H_max,0.0) / tau )

                end if 

            end if

        end do 
        end do
        !!$omp end parallel do

        return 

    end subroutine calc_calving_ground_rate_stress_b12
    
    subroutine calc_calving_ground_rate_stdev(mb_calv,H_ice,f_ice,f_grnd,z_bed_sd,sd_min,sd_max,calv_max,tau,boundaries)
        ! Parameterize grounded ice-margin calving as a function of 
        ! standard deviation of bedrock at each grid point.
        ! Assumes that higher variability in subgrid implies cliffs
        ! that are not represented at low resolution. 

        implicit none 

        real(wp), intent(OUT) :: mb_calv(:,:)             ! [m/yr] Calculated calving rate 
        real(wp), intent(IN)  :: H_ice(:,:)               ! [m] Ice thickness 
        real(wp), intent(IN)  :: f_ice(:,:)               ! [-] Ice area fraction
        real(wp), intent(IN)  :: f_grnd(:,:)              ! [-] Grounded fraction
        real(wp), intent(IN)  :: z_bed_sd(:,:)            ! [m] Standard deviation of bedrock topography
        real(wp), intent(IN)  :: sd_min                   ! [m] stdev(z_bed) at/below which calv=0
        real(wp), intent(IN)  :: sd_max                   ! [m] stdev(z_bed) at/above which calv=calv_max 
        real(wp), intent(IN)  :: calv_max                 ! [m/yr] Maximum allowed calving rate
        real(wp), intent(IN)  :: tau                      ! [yr] Calving timescale       
        character(len=*), intent(IN) :: boundaries 

        ! Local variables
        integer  :: i, j, nx, ny  
        integer  :: im1, ip1, jm1, jp1 
        real(wp) :: f_scale 
        real(wp) :: H_eff
        logical  :: is_grnd_margin 

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        ! Intialize calving rate to zero 
        mb_calv = 0.0 
        
        if (calv_max .gt. 0.0) then 
            ! Determine grounded calving rate 

            !!$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,is_grnd_margin,f_scale,H_eff)
            do j = 1, ny
            do i = 1, nx 

                ! Get neighbor indices
                call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
            
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

                    ! Calculate mass balance calving rate from linear function, 
                    ! limited to available (effective) ice thickness 
                    mb_calv(i,j) = -( min(f_scale*calv_max, H_eff/tau) )
                    
                end if 

            end do 
            end do
            !!$omp end parallel do

        end if 

        return 

    end subroutine calc_calving_ground_rate_stdev

end module calving_aa
