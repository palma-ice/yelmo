module lsf_module

    use yelmo_defs,        only : sp, dp, wp, prec, TOL, TOL_UNDERFLOW, MISSING_VALUE, io_unit_err
    use yelmo_tools,       only : get_neighbor_indices
    use topography,        only : calc_H_eff
    use solver_advection,  only : calc_advec2D
    use mass_conservation, only : calc_G_advec_simple 

    implicit none 

    private 
    
    ! === LSF routines ===
    public :: LSFinit
    public :: LSFadvection
    public :: LSFupdate

    ! === Total CMB_flt (aesthetics) ===
    !public :: calc_cmb_flt
    !public :: calc_cmb_border

    ! === Ocean extrapolation routines ===
    !public :: interpolate_ocn_acx
    !public :: interpolate_ocn_acy
    public :: extrapolate_ocn_laplace

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
        integer,  allocatable :: mask_new_adv(:,:)

        nx = size(lsf,1)
        ny = size(lsf,2)
        allocate(wx(nx,ny))
        allocate(wy(nx,ny))
        allocate(mask_lsf(nx,ny))
        
        ! Initialize variables
        dlsf         = 0.0_wp  ! LSF change in a time dt
        wx           = 0.0_wp  ! retreat-rate x direction (ac-node)
        wy           = 0.0_wp  ! retreat-rate y direction (ac-node)
        mask_lsf     = 1.0_wp  ! Allow all LSF mask to be advected
        !mask_lsf     = 0.0_wp
        !where(lsf .lt. 1.0_wp) mask_lsf = 1.0_wp

        ! net velocity (ice velocity minus calving)
        wx = u_acx + cr_acx
        wy = v_acy + cr_acy
            
        ! Compute the advected LSF field
        call calc_G_advec_simple(dlsf,lsf,mask_lsf,wx,wy, &
                                 mask_adv,solver,boundaries,dx,dt)
        call apply_tendency_lsf(lsf,dlsf,dt,adjust_lsf=.FALSE.)
        
        ! Set border values to ocean values
        lsf(1,:)  = 1.0
        lsf(nx,:) = 1.0
        lsf(:,1)  = 1.0
        lsf(:,ny) = 1.0
        
        ! saturate values to -1 to 1 (helps with stability)
        where(lsf .gt. 1.0)  lsf = 1.0
        where(lsf .lt. -1.0) lsf = -1.0

        ! LSF should not affect points above sea level (check)
        where(H_grnd .gt. 0.0_wp) lsf = -1.0_wp

        if (.FALSE.) then
            ! plot retreat rate instead of calving rate
            cr_acx = wx
            cr_acy = wy
        end if

        return

    end subroutine LSFupdate

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

                ! Store previous value
                lsf_prev = lsf(i,j)
            
                ! Now update lsf with tendency for this timestep
                lsf(i,j) = lsf_prev + dt*lsf_dot(i,j)

                ! Calculate actual current rate of change
                dlsfdt = (lsf(i,j) - lsf_prev) / dt

                ! Update lsf rate to match ice rate of change perfectly
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
        real(wp), allocatable :: mask_change(:,:),mask_change_n(:,:)
        allocate(mask_change(size(mask_orig,1),size(mask_orig,2)))
        allocate(mask_change_n(size(mask_orig,1),size(mask_orig,2)))

        ! Initialize variables
        mask_change   = 1.0_wp
        mask_change_n = 1.0_wp

        where(mask_orig .eq. 0.0_wp) mask_change = 0.0_wp
        where(mask_orig .eq. 0.0_wp) mask_change_n = 0.0_wp
        mask_fill = mask_orig

        count  = 0
        repeat = 1

        ! Initialize the solution array
        !u_x = 0.0
        !u_x_new = 0.0
  
        ! Define boundary conditions
        !u_x(1, :) = 1.0   ! Example boundary condition: top boundary
        !u_x(nx, :) = 0.0  ! Example boundary condition: bottom boundary
        !u_x(:, 1) = 1.0   ! Example boundary condition: left boundary
        !u_x(:, ny) = 0.0  ! Example boundary condition: right boundary
  
        ! Jacobi iteration
        !iter = 0
        !error = tol + 1.0

        !do while (error > tol)
        !    error = 0.0
        !    iter = iter + 1
        !
        !    do i = 2, nx-1
        !      do j = 2, ny-1
        !        u_x_new(i, j) = 0.25 * (u_x(i+1, j) + u_x(i-1, j) + u_x(i, j+1) + u_x(i, j-1))
        !        error = error + abs(u_x_new(i, j) - u_x(i, j))
        !      end do
        !    end do
        !
        !    u_x = u_x_new
        !    write(*,*) 'Iteration:', iter, 'Error:', error
        ! end do

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

    subroutine extrapolate_ocn_acx(mask_fill,mask_orig)
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
        
    end subroutine extrapolate_ocn_acx

    subroutine extrapolate_ocn_laplace(mask_fill, mask_orig)
        ! Routine to extrapolate values using the Laplace equation.
        ! Assumes that value 0 in mask_orig represents ice-free points
        
        implicit none
        
        real(wp), intent(INOUT) :: mask_fill(:,:)
        real(wp), intent(IN) :: mask_orig(:,:)
        
        ! Local variables
        integer :: i, j, iter
        real(wp) :: error, tol
        real(wp), allocatable :: mask_new(:,:)
        
        ! Allocate memory for the temporary array
        allocate(mask_new(size(mask_orig,1), size(mask_orig,2)))
        
        ! Initialize variables
        mask_fill = mask_orig
        mask_new = mask_orig
        tol = 1e-2_wp  ! Tolerance for convergence
        error = tol + 1.0_wp
        iter = 0
        
        ! Jacobi iteration
        do while (error > tol)
            error = 0.0_wp
            iter = iter + 1
        
            do i = 2, size(mask_orig,1)-1
                do j = 2, size(mask_orig,2)-1
                    if (mask_orig(i,j) == 0.0_wp) then
                        mask_new(i,j) = 0.25_wp * (mask_fill(i+1,j) + mask_fill(i-1,j) + mask_fill(i,j+1) + mask_fill(i,j-1))
                        error = error + abs(mask_new(i,j) - mask_fill(i,j))
                    end if
                end do
            end do
        
            mask_fill = mask_new
        end do
        
        deallocate(mask_new)
        
    end subroutine extrapolate_ocn_laplace

end module lsf_module
