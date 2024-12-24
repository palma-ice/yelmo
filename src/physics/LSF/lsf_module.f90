module lsf_module
    ! Definitions for various calving laws 

    use yelmo_defs,        only : sp, dp, wp, prec, TOL_UNDERFLOW
    use yelmo_tools,       only : get_neighbor_indices
    use topography,        only : calc_H_eff
    use solver_advection,  only : calc_advec2D
    use mass_conservation, only : calc_G_advec_simple 

    implicit none 

    private 
    
    !=== Floating calving rates===
    public :: calc_calving_rate_vonmises_m16

    !=== Marine-terminating calving rates ===
    !public :: calc_mici_crawford21

    !=== Land-terminates calving rates ===

    !=== Total CMB_flt (aesthetics) ===
    public :: calc_cmb_flt
    public :: calc_cmb_border

    !=== CalvMIP calving rates ===
    public :: calvmip_exp1
    public :: calvmip_exp2

    !=== LSF routines ===
    public :: LSFinit
    public :: LSFupdate
    public :: LSFborder

contains 

! ===================================================================
!
! Calving related stress/strain routines
!
! ===================================================================

subroutine calc_calving_rate_vonmises_m16(crx_ac,cry_ac,ux_ac,vy_ac,f_ice,f_grnd,tau_eff,dx,tau_max)
    ! Calculate the calving rate [m/yr] based on the 
    ! von Mises stress approach, as outlined by Morlighem et al. (2016)
    ! tau_max=150-750 kPa (fracture strength of ice)
    
    implicit none 

    real(wp), intent(OUT) :: crx_ac(:,:),cry_ac(:,:)
    real(wp), intent(IN)  :: ux_ac(:,:),vy_ac(:,:)
    real(wp), intent(IN)  :: f_ice(:,:),f_grnd(:,:)  
    real(wp), intent(IN)  :: tau_eff(:,:)
    real(wp), intent(IN)  :: dx
    real(wp), intent(IN)  :: tau_max

    ! Local variables 
    integer  :: i, j, nx, ny 
    !real(wp), parameter :: clv_lim = 1e6  ! To avoid really high calving values
                                           ! jablasco: set to grid size dx? CHECK
    nx = size(f_ice,1)
    ny = size(f_ice,2)

    ! Assume square grid cells 
    do j = 1, ny
    do i = 1, nx  
        
        ! Compute over the whole domain
        ! Multiplication of eigenvectors
        crx_ac(i,j) = max( ux_ac(i,j)*tau_eff(i,j)/tau_max, 0.0_wp )
        
        ! Ensure that ice free points do not exhibit calving
        if (f_ice(i,j) .eq. 0.0) then
            crx_ac(i,j) = 0.0_wp
        end if

        ! jablasco: Apply calving limit?
        !clv_rate(i,j) = min(clv_rate(i,j),clv_lim)

    end do
    end do

    return 

end subroutine calc_calving_rate_vonmises_m16

subroutine calc_cmb_flt(cmb_flt,cmb_flt_x,cmb_flt_y,boundaries)

    implicit none

    real(wp), intent(OUT) :: cmb_flt(:,:)                  ! calving on aa-nodes
    real(wp), intent(IN)  :: cmb_flt_x(:,:),cmb_flt_y(:,:) ! calving on ac_nodes
    character(len=*), intent(IN)  :: boundaries

    ! Local variables
    integer  :: i,j,im1,ip1,jm1,jp1,nx,ny
    real(wp) :: cmb_flt_x_aa,cmb_flt_y_aa

    nx = size(cmb_flt_x,1)
    ny = size(cmb_flt_x,2) 

    ! compute to aa nodes
    do j = 1, ny
    do i = 1, nx
        call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
        cmb_flt_x_aa = 0.5_wp*(cmb_flt_x(i,j) + cmb_flt_x(im1,j))
        cmb_flt_y_aa = 0.5_wp*(cmb_flt_y(i,j) + cmb_flt_y(i,jm1))
        cmb_flt(i,j) = (cmb_flt_x_aa*cmb_flt_x_aa + cmb_flt_y_aa*cmb_flt_y_aa)**0.5
    end do
    end do

    return

end subroutine calc_cmb_flt

subroutine calc_cmb_border(crx_ac,cry_ac,lsf_aa,boundaries)

    implicit none

    real(wp), intent(INOUT) :: crx_ac(:,:),cry_ac(:,:) ! [m/yr] calving-rates on ac-nodes 
    real(wp), intent(IN)    :: lsf_aa(:,:)             ! LSF mask on aa-nodes
    character(len=*), intent(IN)  :: boundaries

    ! Local variables
    integer  :: i,j,im1,ip1,jm1,jp1,nx,ny
   
    do j=1,ny
    do i=1,nx
        call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
        ! x-direction
        if (lsf_aa(i,j) .gt. 0.0 .and. lsf_aa(ip1,j) .le. 0.0) then
            crx_ac(i,j)   = crx_ac(i,j)
        else if (lsf_aa(i,j) .gt. 0.0 .and. lsf_aa(im1,j) .le. 0.0) then
            crx_ac(im1,j) = crx_ac(im1,j)
        else
            crx_ac(i,j)   = 0.0
        end if

        ! y-direction
        if (lsf_aa(i,j) .gt. 0.0 .and. lsf_aa(i,jp1) .le. 0.0) then
            cry_ac(i,j)   = cry_ac(i,j)
        else if (lsf_aa(i,j) .gt. 0.0 .and. lsf_aa(i,jm1) .le. 0.0) then
            cry_ac(i,jm1) = cry_ac(i,jm1)
        else
            cry_ac(i,j)   = 0.0
        end if
 
    end do
    end do

    return

end subroutine calc_cmb_border

! ===================================================================
!
!                      CalvMIP experiments
!
! ===================================================================

subroutine calvmip_exp1(crx_ac,cry_ac,ux_ac,vy_ac,lsf_aa,dx,boundaries)
! Experiment 1 of CalvMIP

    implicit none

    real(wp), intent(OUT) :: crx_ac(:,:),cry_ac(:,:)   ! Calving rates on ac-nodes
    real(wp), intent(IN)  :: ux_ac(:,:),vy_ac(:,:)     ! Velocities on ac-nodes
    real(wp), intent(IN)  :: lsf_aa(:,:)               ! LSF mask on aa-nodes
    real(wp), intent(IN)  :: dx                        ! Ice resolution
    character(len=*), intent(IN)  :: boundaries        ! Boundary conditions to impose

    ! Local variables
    integer  :: i, j, im1, ip1, jm1, jp1, nx, ny
    real(wp) :: r

    ! points inside the radius 750km should have no calving rate
    nx = size(ux_ac,1)
    ny = size(ux_ac,2)

    ! test for LSF improvement
    crx_ac = 0.0_wp
    cry_ac = 0.0_wp

    r = 0.0_wp
    do j = 1, ny
    do i = 1, nx

        ! aa-nodes indices
        r = sqrt((0.5*(nx+1)-i)*(0.5*(nx+1)-i) + (0.5*(ny+1)-j)*(0.5*(ny+1)-j))*dx

        if (r .ge. 750e3) then
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
            ! x-direction
            if (lsf_aa(i,j) .gt. 0.0 .and. lsf_aa(ip1,j) .le. 0.0) then
                crx_ac(i,j)   = -ux_ac(i,j)
            else if (lsf_aa(i,j) .gt. 0.0 .and. lsf_aa(im1,j) .le. 0.0) then
                crx_ac(im1,j) = -ux_ac(im1,j)
            end if

            ! y-direction
            if (lsf_aa(i,j) .gt. 0.0 .and. lsf_aa(i,jp1) .le. 0.0) then
                cry_ac(i,j)   = -vy_ac(i,j)
            else if (lsf_aa(i,j) .gt. 0.0 .and. lsf_aa(i,jm1) .le. 0.0) then
                cry_ac(i,jm1) = -vy_ac(i,jm1)
            end if

        end if

    end do
    end do

    return

end subroutine calvmip_exp1

subroutine calvmip_exp2(crx_ac,cry_ac,ux_ac,vy_ac,time,boundaries)
    ! Experiment 2 of CalvMIP

    implicit none

    real(wp), intent(OUT) :: crx_ac(:,:), cry_ac(:,:)
    real(wp), intent(IN)  :: ux_ac(:,:),vy_ac(:,:)
    real(wp), intent(IN)  :: time
    character(len=*), intent(IN)  :: boundaries             ! Boundary conditions to impose

    ! local variables
    integer :: i, j, im1, ip1, jm1, jp1, nx, ny
    real(wp), parameter   :: pi = acos(-1.0)  ! Calculate pi intrinsically
    real(wp)              :: wv   

    nx = size(ux_ac,1)
    ny = size(ux_ac,2) 

    wv = 300.0 * sin(2.0 * pi * time / 1000.0) 

    do j = 1, ny
    do i = 1, nx
        call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

        ! define over whole domain
        ! x-direction
        if(ux_ac(i,j) .gt. 0.0_wp) then
            crx_ac(i,j) = MIN(-(ux_ac(i,j)+wv),0.0_wp)
        else if (ux_ac(i,j) .lt. 0.0_wp) then
            crx_ac(i,j) = MAX(-(ux_ac(i,j)-wv),0.0_wp)
        else
            crx_ac(i,j) = 0.0_wp
        end if

        ! y-direction
        if(vy_ac(i,j) .gt. 0.0_wp) then
            cry_ac(i,j) = MIN(-(vy_ac(i,j)+wv),0.0_wp)
        else if (vy_ac(i,j) .lt. 0.0_wp) then
            cry_ac(i,j) = MAX(-(vy_ac(i,j)-wv),0.0_wp)
        else
            cry_ac(i,j) = 0.0_wp
        end if

    end do
    end do

    return

end subroutine calvmip_exp2

! ===================================================================
!
!                        LSF functions
!
! ===================================================================

subroutine LSFinit(LSF,H_ice,z_bed,dx)

    implicit none

    real(wp), intent(OUT) :: LSF(:,:)      ! LSF mask 
    real(wp), intent(IN)  :: H_ice(:,:)    ! Ice thickness
    real(wp), intent(IN)  :: z_bed(:,:)   ! bedrock elevation     
    real(wp), intent(IN)  :: dx            ! Model resolution

    ! Initialize LSF value at ocean value
    LSF = -1.0_wp  
    
    ! Assign values
    where(H_ice .gt. 0.0_wp) LSF = 1.0_wp
    where(z_bed .gt. 0.0_wp) LSF = 1.0_wp 

    call eikonal_equation(LSF,z_bed)
    LSF = LSF * dx
    !where(LSF .le. 0.0) LSF = -1.0

    return
    
end subroutine LSFinit

subroutine LSFupdate(LSFn,LSF,crx_ac,cry_ac,ux_ac,vy_ac,H_grnd,var_dot,mask_adv,dx,dy,dt,solver,boundaries)

    implicit none

    real(wp),       intent(OUT)   :: LSFn(:,:)               ! new LSF field (aa-nodes)
    real(wp),       intent(IN)    :: LSF(:,:)                ! LSF to be advected (aa-nodes)
    real(wp),       intent(IN)    :: crx_ac(:,:),cry_ac(:,:) ! [m/yr] calving rate (vertical)
    !real(wp),       intent(IN)    :: f_ice(:,:)              ! Ice fraction to be advected
    real(wp),       intent(IN)    :: ux_ac(:,:)              ! [m/a] 2D velocity, x-direction (ac-nodes)
    real(wp),       intent(IN)    :: vy_ac(:,:)              ! [m/a] 2D velocity, y-direction (ac-nodes)
    real(wp),       intent(IN)    :: H_grnd(:,:)             ! [m] floating contribution?
    real(wp),       intent(IN)    :: var_dot(:,:)            ! [dvar/dt] Source term for variable (needed?)
    integer,        intent(IN)    :: mask_adv(:,:)           ! Advection mask
    real(wp),       intent(IN)    :: dx                      ! [m] Horizontal resolution, x-direction
    real(wp),       intent(IN)    :: dy                      ! [m] Horizontal resolution, y-direction
    real(wp),       intent(IN)    :: dt                      ! [a]   Timestep
    character(len=*), intent(IN)  :: solver                  ! Solver to use for the ice thickness advection equation
    character(len=*), intent(IN)  :: boundaries              ! Boundary conditions to impose

    ! Local variables
    integer  :: i, j, im1, ip1, jm1, jp1, nx, ny
    real(wp), allocatable :: dlsf(:,:), wx(:,:), wy(:,:), mask_lsf(:,:)
    integer, allocatable  :: mask_new_adv(:,:)

    nx = size(LSF,1)
    ny = size(LSF,2)
    allocate(dlsf(nx,ny))
    allocate(wx(nx,ny))
    allocate(wy(nx,ny))
    allocate(mask_new_adv(nx,ny))
    allocate(mask_lsf(nx,ny))

    ! Advection depends on ice velocity minus the calving rate
    dlsf     = 0.0_wp
    mask_lsf = 0.0_wp
    mask_new_adv = 1
    mask_lsf = 1.0_wp ! Allow all LSF mask to be advected
    !where(LSF .gt. 0.0_wp) mask_lsf = 1.0_wp ! Ice fraction to be advected

    ! interpolate velocities to open ocean
    !call interpolatex_missing_iterative(crx_ac,ux_ac)
    !call interpolatey_missing_iterative(cry_ac,vy_ac)

    ! calving rate has opposite sign as velocity
    wx = (ux_ac + crx_ac)/(ux_ac + 1e-6)
    wy = (vy_ac + cry_ac)/(vy_ac + 1e-6)
 
    ! Compute the advected LSF field
    call calc_G_advec_simple(dlsf,lsf,mask_lsf,wx,wy, &
                             mask_new_adv,solver,boundaries,dx,dt) ! changed mask_adv for 1
    LSFn = dlsf*dt + lsf

    ! Allow for in between values only in the zero border
    !do j = 1,ny
    !do i = 1,nx
    !
    !    ! Get neighbor indices
    !    call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
    !
    !    if ((LSFn(i,j) .gt. 0.0_wp) .and. (LSFn(im1,j) .gt. 0.0_wp) .and. (LSFn(ip1,j) .gt. 0.0_wp) &
    !        .and. (LSFn(i,jm1) .gt. 0.0_wp) .and. (LSFn(i,jp1) .gt. 0.0_wp)) then
    !        LSFn(i,j) = 1.0_wp
    ! 
    !    !else if ((LSFn(i,j) .le. 0.0_wp) .and. (LSFn(im1,j) .le. 0.0_wp) .and. (LSFn(ip1,j) .le. 0.0_wp) &
    !    !    .and. (LSFn(i,jm1) .le. 0.0_wp) .and. (LSFn(i,jp1) .le. 0.0_wp)) then
    !    !    LSFn(i,j) = -1.0_wp
    !
    !    end if
    !
    !end do
    !end do

    ! compute distance to ice front
    !call eikonal_equation(LSFn,H_grnd)
    !LSFn = LSFn*dx
    !where(LSFn .le. 0.0) LSFn = -1.0

    ! saturate values to -1 to 1 (helps with stability)
    where(LSFn .gt. 1.0)  LSFn = 1.0
    where(LSFn .lt. -1.0) LSFn = -1.0

    ! LSF should not affect points above sea level (maybe in a future it can)
    where(H_grnd .gt. 0.0_wp) LSFn = 1.0_wp

    return

end subroutine LSFupdate

subroutine LSFborder(LSF,f_ice,boundaries)
    ! We will set the ice border as value zero
    ! Test to see if it helps with stability and advance + retreat
    
    implicit none

    real(wp),       intent(OUT)   :: LSF(:,:)               ! new LSF field (aa-nodes)
    real(wp),       intent(IN)    :: f_ice(:,:)              ! Ice fraction to be advected
    character(len=*), intent(IN)  :: boundaries              ! Boundary conditions to impose

    ! Local variables
    integer  :: i, j, im1, ip1, jm1, jp1, nx, ny

    nx = size(LSF,1)
    ny = size(LSF,2)
 
    do j=1,ny
    do i=1,nx
        if ((f_ice(i,j) .gt. 0.0_wp) .and. ((f_ice(im1,j) .eq. 0.0_wp) .or. (f_ice(ip1,j) .eq. 0.0_wp) &
                .or. (f_ice(i,jm1) .eq. 0.0_wp) .or. (f_ice(i,jp1) .eq. 0.0_wp))) then
                LSF(i,j) = 0.0_wp
        end if
    end do
    end do

    return

end subroutine LSFborder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                    !
!         Internal functions         !
!                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine eikonal_equation(lsf,H_grnd)
    ! internal function to determine distance to the ice front 
    ! used for LSF mask
    ! based on PICO module

    implicit none

    real(prec), intent(INOUT) :: lsf(:,:)
    real(prec), intent(IN)    :: H_grnd(:,:)

    ! Local variables
    integer :: i, j, nx, ny
    real(prec) :: loop, current_label_ice, current_label_ocn
    real(prec), allocatable :: dist_ice_front(:,:)
    real(prec), allocatable :: mask_lsf(:,:)

    nx = size(lsf,1)
    ny = size(lsf,2)
    allocate(dist_ice_front(nx,ny))
    allocate(mask_lsf(nx,ny))
    
    current_label_ice = 1.0
    current_label_ocn = -1.0
    
    ! Define border
    mask_lsf = 0.0 ! Init mask to zero
    do i = 2, nx-1
    do j = 1, ny-1
        if (lsf(i,j) .gt. 0.0 .and. ((lsf(i-1,j) .le. 0.0) .or.  (lsf(i+1,j) .le. 0.0) .or. &
                                     (lsf(i,j-1) .le. 0.0) .or.  (lsf(i,j+1) .le. 0.0))) then
            mask_lsf(i,j) = 1.0
        else if (lsf(i,j) .le. 0.0 .and. ((lsf(i-1,j) .gt. 0.0) .or.  (lsf(i+1,j) .gt. 0.0) .or. &
                                     (lsf(i,j-1) .gt. 0.0) .or.  (lsf(i,j+1) .gt. 0.0))) then
            mask_lsf(i,j) = -1.0
        end if
    end do
    end do 

    ! Assign to dists mask
    dist_ice_front = mask_lsf !0.0
    !where(lsf .le. 0.0) dist_ice_front = 1.0
    ! LSF should not affect points above sea level
    !where(H_grnd .ge. 0.0) dist_ice_front = 0.0

    ! compute distance to mask
    loop = 1.0
    do while(loop .ne. 0.0)
        loop = 0.0

        do i = 2, nx-1
        do j = 2, ny-1

            ! lsf positive points (ice)
            if (lsf(i,j) .gt. 0.0 .and. dist_ice_front(i,j) .eq. 0.0) then
                if(((dist_ice_front(i-1,j) .eq. current_label_ice)) .or. (dist_ice_front(i+1,j) .eq. current_label_ice) .or. &
                    (dist_ice_front(i,j-1) .eq. current_label_ice) .or. (dist_ice_front(i,j+1) .eq. current_label_ice)) then
                    dist_ice_front(i,j) = current_label_ice + 1.0
                    loop = 1.0
                end if
            else if (lsf(i,j) .le. 0.0 .and. dist_ice_front(i,j) .eq. 0.0) then
                if(((dist_ice_front(i-1,j) .eq. current_label_ocn)) .or. (dist_ice_front(i+1,j) .eq. current_label_ocn) .or. &
                    (dist_ice_front(i,j-1) .eq. current_label_ocn) .or. (dist_ice_front(i,j+1) .eq. current_label_ocn)) then
                    dist_ice_front(i,j) = current_label_ocn - 1.0
                    loop = 1.0
                end if
            end if
        end do
        end do

        current_label_ice = current_label_ice+1.0
        current_label_ocn = current_label_ocn-1.0

    end do

    ! jablasco: correct distance! (distance for ice front grid is 0 not 1)
    where(dist_ice_front .gt. 0.0) dist_ice_front = dist_ice_front-1.0
    !where(lsf .lt. 0.0) dist_ice_front = -1.0
    where(0.0 .lt. lsf .and. lsf .lt. 1.0) dist_ice_front = lsf ! at the ice front the number is fractional
    !where(lsf .ge. 0.0) lsf = dist_ice_front

    ! new lsf mask is dist_ice_front
    lsf = dist_ice_front

    return

end subroutine eikonal_equation

subroutine interpolatex_missing_iterative(array,mask)
   
    implicit none

    real(prec), intent(INOUT) :: array(:,:)
    real(prec), intent(IN)    :: mask(:,:)

    ! Local variables
    integer  :: i, j, nx, ny, ni
    real(wp) :: sum, count
    logical  :: has_changed

    nx = size(array,1)
    ny = size(array,2)

    do
        has_changed = .false.

        do i = 2, nx - 1
            do j = 2, ny - 1
                if (array(i, j) .eq. 0.0 .and. mask(i,j) .eq. 0.0) then
                    sum = 0.0
                    count = 0.0

                    ! Check x-neighbors
                    do ni = -1, 1
                        if (array(i + ni, j) /= 0.0) then
                            sum = sum + array(i + ni, j)
                            count = count + 1.0
                        end if
                    end do

                    ! Interpolate if neighbors are available
                    if (count > 0.0) then
                        array(i, j) = sum / count
                        has_changed = .true.
                    end if
                end if
            end do
        end do

        ! Exit loop if no more changes
        if (.not. has_changed) exit
    end do
end subroutine interpolatex_missing_iterative

subroutine interpolatey_missing_iterative(array,mask)
   
    implicit none

    real(prec), intent(INOUT) :: array(:,:)
    real(prec), intent(IN)    :: mask(:,:)

    ! Local variables
    integer  :: i, j, nx, ny, nj
    real(wp) :: sum, count
    logical  :: has_changed

    nx = size(array,1)
    ny = size(array,2)

    do
        has_changed = .false.

        do i = 2, nx - 1
            do j = 2, ny - 1
                if (array(i, j) .eq. 0.0 .and. mask(i,j) .eq. 0.0) then
                    sum = 0.0
                    count = 0.0

                    ! Check y-neighbors
                    do nj = -1, 1
                        if (array(i, j + nj) /= 0.0) then
                            sum = sum + array(i, j + nj)
                            count = count + 1.0
                        end if
                    end do

                    ! Interpolate if neighbors are available
                    if (count > 0.0) then
                        array(i, j) = sum / count
                        has_changed = .true.
                    end if
                end if
            end do
        end do

        ! Exit loop if no more changes
        if (.not. has_changed) exit
    end do
end subroutine interpolatey_missing_iterative

end module lsf_module
