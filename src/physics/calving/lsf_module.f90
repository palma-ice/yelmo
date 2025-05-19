module lsf_module

    use yelmo_defs,        only : sp, dp, wp, prec, TOL, TOL_UNDERFLOW, MISSING_VALUE, io_unit_err
    use yelmo_tools,       only : get_neighbor_indices
    use topography,        only : calc_H_eff
    use solver_advection,  only : calc_advec2D
    use mass_conservation, only : calc_G_advec_simple 

    implicit none 

    private 
    
    !=== Floating calving rates===
    !public :: calc_calving_rate_vonmises_m16

    !=== Marine-terminating calving rates ===
    !public :: calc_mici_crawford21

    !=== Land-terminates calving rates ===

    !=== LSF routines ===
    public :: LSFinit
    public :: LSFadvection
    public :: LSFupdate
    public :: LSFborder

    !=== Total CMB_flt (aesthetics) ===
    public :: calc_cmb_flt
    public :: calc_cmb_border

    ! open ocean extrapo    lation
    public :: interpolate_ocn_acx
    public :: interpolate_ocn_acy

contains 
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
        LSF = 1.0_wp  
        
        ! Assign values
        where(H_ice .gt. 0.0_wp) LSF = -1.0_wp
        where(z_bed .gt. 0.0_wp) LSF = -1.0_wp 

        !call eikonal_equation(LSF,z_bed)
        !LSF = LSF * dx
        !where(LSF .le. 0.0) LSF = -1.0

        return
        
    end subroutine LSFinit

    subroutine LSFupdate(dlsf,lsf,cr_acx,cr_acy,u_acx,v_acy,H_grnd,var_dot,mask_adv,dx,dy,dt,solver,boundaries)

        implicit none

        real(wp),       intent(INOUT) :: dlsf(:,:)               ! advected LSF field
        real(wp),       intent(INOUT) :: lsf(:,:)                ! LSF to be advected (aa-nodes)
        real(wp),       intent(INOUT) :: cr_acx(:,:),cr_acy(:,:) ! [m/yr] calving rate (vertical)
        real(wp),       intent(IN)    :: u_acx(:,:)              ! [m/a] 2D velocity, x-direction (ac-nodes)
        real(wp),       intent(IN)    :: v_acy(:,:)              ! [m/a] 2D velocity, y-direction (ac-nodes)
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
        real(wp), allocatable :: wx(:,:), wy(:,:), mask_lsf(:,:)
        real(wp), allocatable :: ux_domain_ac(:,:), vy_domain_ac(:,:)
        integer,  allocatable :: mask_new_adv(:,:)

        nx = size(lsf,1)
        ny = size(lsf,2)
        allocate(wx(nx,ny))
        allocate(wy(nx,ny))
        allocate(mask_new_adv(nx,ny))
        allocate(mask_lsf(nx,ny))
        allocate(ux_domain_ac(nx,ny))
        allocate(vy_domain_ac(nx,ny))
        
        ! Advection depends on ice velocity minus the calving rate
        dlsf     = 0.0_wp
        mask_lsf = 0.0_wp
        mask_new_adv = 1
        !mask_lsf = 1.0_wp ! Allow all LSF mask to be advected
        where(lsf .le. 0.0_wp) mask_lsf = 1.0_wp ! Ice fraction to be advected

        ! interpolate velocities into the open ocean (test)
        !ux_domain_ac = u_acx
        !vy_domain_ac = v_acy
        !call interpolatex_missing_iterative(ux_domain_ac,u_acx)
        !call interpolatey_missing_iterative(vy_domain_ac,v_acy)

        ! calving rate has opposite sign as velocity
        !wx = ux_domain_ac + cr_acx
        !wy = vy_domain_ac + cr_acy
        wx = u_acx + cr_acx
        wy = v_acy + cr_acy

        ! Compute the advected LSF field
        call calc_G_advec_simple(dlsf,lsf,mask_lsf,wx,wy, &
                                mask_new_adv,solver,boundaries,dx,dt) ! changed mask_adv for 1
        call apply_tendency_lsf(lsf,dlsf,dt,adjust_lsf=.FALSE.)
        
        ! Set border values to ocean values
        lsf(1,:)  = 1.0
        lsf(nx,:) = 1.0
        lsf(:,1)  = 1.0
        lsf(:,ny) = 1.0
        
        ! saturate values to -1 to 1 (helps with stability)
        where(lsf .gt. 1.0)  lsf = 1.0
        where(lsf .lt. -1.0) lsf = -1.0

        ! LSF should not affect points above sea level (maybe in a future it can)
        where(H_grnd .gt. 0.0_wp) lsf = -1.0_wp

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

    ! ===================================================================
    !
    ! Total CMB_flt (aesthetics)
    !
    ! ===================================================================

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

    subroutine calc_cmb_border(cr_acx,cr_acy,lsf_aa,boundaries)

        implicit none

        real(wp), intent(INOUT) :: cr_acx(:,:),cr_acy(:,:) ! [m/yr] calving-rates on ac-nodes 
        real(wp), intent(IN)    :: lsf_aa(:,:)             ! LSF mask on aa-nodes
        character(len=*), intent(IN)  :: boundaries

        ! Local variables
        integer  :: i,j,im1,ip1,jm1,jp1,nx,ny
    
        do j=1,ny
        do i=1,nx
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
            ! x-direction
            if (lsf_aa(i,j) .gt. 0.0 .and. lsf_aa(ip1,j) .le. 0.0) then
                cr_acx(i,j)   = cr_acx(i,j)
            else if (lsf_aa(i,j) .gt. 0.0 .and. lsf_aa(im1,j) .le. 0.0) then
                cr_acx(im1,j) = cr_acx(im1,j)
            else
                cr_acx(i,j)   = 0.0
            end if

            ! y-direction
            if (lsf_aa(i,j) .gt. 0.0 .and. lsf_aa(i,jp1) .le. 0.0) then
                cr_acy(i,j)   = cr_acy(i,j)
            else if (lsf_aa(i,j) .gt. 0.0 .and. lsf_aa(i,jm1) .le. 0.0) then
                cr_acy(i,jm1) = cr_acy(i,jm1)
            else
                cr_acy(i,j)   = 0.0
            end if
    
        end do
        end do

        return

    end subroutine calc_cmb_border

    ! ===================================================================
    !
    ! Internal functions
    !
    ! ===================================================================
    
    subroutine eikonal_equation(lsf,H_grnd)
        ! internal function to determine distance to the ice front 
        ! used for LSF mask
        ! based on PICO module

        implicit none

        real(prec), intent(INOUT) :: lsf(:,:)
        real(prec), intent(IN)    :: H_grnd(:,:)

        ! Local variables
        integer :: i, j, nx, ny
        real(prec) :: loop, current_label_ice !, current_label_ocn
        real(prec), allocatable :: dist_ice_front(:,:)
        real(prec), allocatable :: mask_lsf(:,:)

        nx = size(lsf,1)
        ny = size(lsf,2)
        allocate(dist_ice_front(nx,ny))
        allocate(mask_lsf(nx,ny))
        
        current_label_ice = 1.0
        !current_label_ocn = -1.0
        
        ! Define border (do not take the border which is set to zero)
        !mask_lsf = 0.0 ! Init mask to zero
        !do i = 3, nx-2
        !do j = 3, ny-2
        !    if (lsf(i,j) .gt. 0.0 .and. ((lsf(i-1,j) .le. 0.0) .or.  (lsf(i+1,j) .le. 0.0) .or. &
        !                                 (lsf(i,j-1) .le. 0.0) .or.  (lsf(i,j+1) .le. 0.0))) then
        !        mask_lsf(i,j) = 1.0
        !    !else if (lsf(i,j) .lt. 0.0 .and. ((lsf(i-1,j) .ge. 0.0) .or.  (lsf(i+1,j) .ge. 0.0) .or. &
        !    !                             (lsf(i,j-1) .ge. 0.0) .or.  (lsf(i,j+1) .ge. 0.0))) then
        !    !    mask_lsf(i,j) = -1.0
        !    end if
        !end do
        !end do 

        ! Assign to dists mask
        dist_ice_front = 0.0
        where(lsf .le. 0.0) dist_ice_front = 1.0
        ! LSF should not affect points above sea level
        where(H_grnd .ge. 0.0) dist_ice_front = 0.0

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
                !else if (lsf(i,j) .lt. 0.0 .and. dist_ice_front(i,j) .eq. 0.0) then
                !    if(((dist_ice_front(i-1,j) .eq. current_label_ocn)) .or. (dist_ice_front(i+1,j) .eq. current_label_ocn) .or. &
                !        (dist_ice_front(i,j-1) .eq. current_label_ocn) .or. (dist_ice_front(i,j+1) .eq. current_label_ocn)) then
                !        dist_ice_front(i,j) = current_label_ocn - 1.0
                !        loop = 1.0
                !    end if
                end if
            end do
            end do

            current_label_ice = current_label_ice+1.0
            !current_label_ocn = current_label_ocn-1.0

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

        return

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
                    if (array(i, j) .eq. 0.0 .and. mask(i,j) == 0.0) then
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

        return

    end subroutine interpolatey_missing_iterative

    subroutine apply_tendency_lsf(lsf,lsf_dot,dt,adjust_lsf)
            
        implicit none
            
        real(wp), intent(INOUT) :: lsf(:,:)
        real(wp), intent(INOUT) :: lsf_dot(:,:)
        real(wp), intent(IN)    :: dt
        logical, optional, intent(IN) :: adjust_lsf
            
        ! Local variables
        integer :: i, j, nx, ny 
        real(wp) :: lsf_prev 
        real(wp) :: dlsfdt
        logical  :: allow_adjust_lsf
            
        if (dt .gt. 0.0) then
            ! Only apply this routine if dt > 0!
            
            allow_adjust_lsf = .FALSE.
            if (present(adjust_lsf)) allow_adjust_lsf = adjust_lsf
            
            nx = size(lsf,1)
            ny = size(lsf,2)
            
            !!$omp parallel do collapse(2) private(i,j,lsf_prev,dlsfdt)
            do j = 1, ny
            do i = 1, nx

                ! Store previous ice thickness
                lsf_prev = lsf(i,j)
            
                ! Now update lsf with tendency for this timestep
                lsf(i,j) = lsf_prev + dt*lsf_dot(i,j)

                ! Ensure tiny numeric LSF values tend to zero are removed
                ! check effect
                !if (abs(lsf(i,j)) .lt. TOL) lsf(i,j) = 0.0

                ! Calculate actual current rate of change
                dlsfdt = (lsf(i,j) - lsf_prev) / dt

                ! Update mb rate to match ice rate of change perfectly
                if (allow_adjust_lsf) then
                    lsf_dot(i,j) = dlsfdt
                end if

            end do
            end do
            !!$omp end parallel do

        end if

        return

    end subroutine apply_tendency_lsf

    ! Simple advection test (works)
    subroutine LSFadvection(LSF,zbed,dx)
        
        implicit none
    
        real(wp), intent(OUT) :: LSF(:,:)      ! LSF mask
        real(wp), intent(IN)  :: zbed(:,:)    
        real(wp), intent(IN)  :: dx            ! Model resolution
        
        ! Internal variables
        real(wp) :: rc, xc, yc
        integer  :: i,j,nx,ny
    
        nx = size(zbed,1)
        ny = size(zbed,2)
        rc = 12.5
        xc = 25
        yc = 25 
    
        do j=1,ny
        do i=1,nx
    
        LSF(i,j) = ((0.001*(i-1)*dx-xc)**2+(0.001*(j-1)*dx-yc)**2)**0.5 -rc 
    
        end do
        end do
    
        return
    
    end subroutine LSFadvection

    subroutine interpolate_ocn_acy(mask_fill,mask_orig)
        ! Routine to extrapolate the nearest value in the y-direction.
        ! So far we assume that value 0 in mask_orig represents ice-free points

        implicit none
    
        real(prec), intent(INOUT) :: mask_fill(:,:)
        real(prec), intent(IN)    :: mask_orig(:,:)
    
        ! Local variables
        integer  :: i, j
        integer :: count, repeat
        real(wp), allocatable :: mask_change(:,:)
        allocate(mask_change(size(mask_orig,1),size(mask_orig,2)))

        ! Initialize variables
        mask_change = 1.0_wp
        where(mask_orig .eq. 0.0_wp) mask_change = 0.0_wp
        mask_fill = mask_orig

        count  = 0
        repeat = 1
        if (SUM(mask_orig) .eq. 0.0_wp) then
            ! do nothing
        else
            do while (repeat .eq. 1)
                ! Interpolate neighbor value
                do i = 2, size(mask_orig,1)-1
                    do j = 2, size(mask_orig,2)-1
                        if      (mask_change(i,j) .eq. 0.0_wp .and. mask_change(i,j+1) .eq. 1.0) then
                            mask_fill(i,j)   = mask_fill(i,j+1) 
                            mask_change(i,j) = 1.0_wp
                            count = count+1
                        else if (mask_change(i,j) .eq. 0.0_wp .and. mask_change(i,j-1) .eq. 1.0) then
                            mask_fill(i,j)   = mask_fill(i,j-1) 
                            mask_change(i,j) = 1.0_wp
                            count = count+1
                        end if
                    end do
                end do

                ! Repeat if changes occured
                if (count .eq. 0) then 
                    repeat = 0
                    count  = 0
                else
                    repeat = 1
                    count  = 0
                end if  

            end do    
        end if

        return
    
    end subroutine interpolate_ocn_acy

    subroutine interpolate_ocn_acx(mask_fill,mask_orig)
        ! Routine to extrapolate the nearest value in the x-direction.
        ! So far we assume that value 0 in mask_orig represents ice-free points

        implicit none
        
        real(prec), intent(INOUT) :: mask_fill(:,:)
        real(prec), intent(IN)    :: mask_orig(:,:)
        
        ! Local variables
        integer  :: i, j
        integer :: count, repeat
        real(wp), allocatable :: mask_change(:,:)
        allocate(mask_change(size(mask_orig,1),size(mask_orig,2)))
    
        ! Initialize variables
        mask_change = 1.0_wp
        where(mask_orig .eq. 0.0_wp) mask_change = 0.0_wp
        mask_fill = mask_orig
    
        count  = 0
        repeat = 1
        if (SUM(mask_orig) .eq. 0.0_wp) then
            ! do nothing
        else
            do while (repeat .eq. 1)
                ! Interpolate neighbor value
                do i = 2, size(mask_orig,1)-1
                    do j = 2, size(mask_orig,2)-1
                        if      (mask_change(i,j) .eq. 0.0_wp .and. mask_change(i+1,j) .eq. 1.0) then
                            mask_fill(i,j)   = mask_fill(i+1,j) 
                            mask_change(i,j) = 1.0_wp
                            count = count+1
                        else if (mask_change(i,j) .eq. 0.0_wp .and. mask_change(i-1,j) .eq. 1.0) then
                            mask_fill(i,j)   = mask_fill(i-1,j) 
                            mask_change(i,j) = 1.0_wp
                            count = count+1
                        end if
                    end do
                end do
    
                ! Repeat if changes occured
                if (count .eq. 0) then 
                    repeat = 0
                    count  = 0
                else
                    repeat = 1
                    count  = 0
                end if 

            end do    
        end if
    
        return
        
    end subroutine interpolate_ocn_acx

end module lsf_module
