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
    public :: LSFupdate

    ! === Ocean extrapolation routines ===
    public :: extrapolate_ocn_acx
    public :: extrapolate_ocn_acy
    public :: extrapolate_ocn_laplace_simple

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

    subroutine LSFupdate(dlsf,lsf,cr_acx,cr_acy,u_acx,v_acy,var_dot,mask_adv,dx,dy,dt,solver)

        implicit none

        real(wp),       intent(INOUT) :: dlsf(:,:)               ! advected LSF field
        real(wp),       intent(INOUT) :: lsf(:,:)                ! LSF to be advected (aa-nodes)
        real(wp),       intent(INOUT) :: cr_acx(:,:),cr_acy(:,:) ! [m/yr] calving rate (vertical)
        real(wp),       intent(IN)    :: u_acx(:,:)              ! [m/a] 2D velocity, x-direction (ac-nodes)
        real(wp),       intent(IN)    :: v_acy(:,:)              ! [m/a] 2D velocity, y-direction (ac-nodes)
        real(wp),       intent(IN)    :: var_dot(:,:)            ! [dvar/dt] Source term for variable (needed?)
        integer,        intent(IN)    :: mask_adv(:,:)           ! Advection mask
        real(wp),       intent(IN)    :: dx                      ! [m] Horizontal resolution, x-direction
        real(wp),       intent(IN)    :: dy                      ! [m] Horizontal resolution, y-direction
        real(wp),       intent(IN)    :: dt                      ! [a]   Timestep
        character(len=*), intent(IN)  :: solver                  ! Solver to use for the ice thickness advection equation

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

        ! lsf velocity
        wx = u_acx + cr_acx
        wy = v_acy + cr_acy   

        ! Compute the advected LSF field
        if (.TRUE.) then
            call calc_advec2D(dlsf,lsf,mask_lsf,wx,wy,var_dot, &
                                mask_adv,dx,dy,dt,solver,"periodic")
            call apply_tendency_lsf(lsf,dlsf,dt,adjust_lsf=.FALSE.)
        else
            ! Simple advecter without diagonilizing. Test.
            call LSFadvec_simple(dlsf,lsf, wx, wy, dt, dx, "periodic")
        end if
        
        ! saturate values to -1 to 1 (helps with stability)
        where(lsf .gt. 1.0)  lsf = 1.0
        where(lsf .lt. -1.0) lsf = -1.0

        if (.FALSE.) then
            ! plot retreat rate instead of calving rate
            cr_acx = wx
            cr_acy = wy
        end if

        return

    end subroutine LSFupdate

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

                ! Update lsf rate to match rate of change perfectly
                if (allow_adjust_lsf) then
                    lsf_dot(i,j) = dlsfdt
                end if

            end do
            end do
            !!$omp end parallel do

        end if

        return

    end subroutine apply_tendency_lsf
        
    subroutine LSFadvec_simple(dlsf,LSF, u, v, dt, dx, boundaries)
        ! Simple LSF advection routine. Not diagonilized.
        ! Test

        implicit none
            
        ! Define input and output variables
        real(wp), intent(INOUT) :: dlsf(:,:)    ! aa-node
        real(wp), intent(INOUT) :: LSF(:,:)     ! aa-node
        real(wp), intent(IN) :: u(:,:), v(:,:)  ! ac-node
        real(wp), intent(IN) :: dt, dx
        character(len=*), intent(IN)  :: boundaries
            
        ! Local variables
        real(wp) :: dtdx
        real(wp), dimension(size(LSF,1), size(LSF,2)) :: dLSF_acx, dLSF_acy ! ac-nodes
        real(wp), dimension(size(LSF,1), size(LSF,2)) :: qx_ac, qy_ac, qx_aa, qy_aa
        integer :: i, j, im1, ip1, jm1, jp1, nx, ny
            
        dlsf     = 0.0_wp
        dtdx     = dt / dx
        dLSF_acx = 0.0_wp
        dLSF_acy = 0.0_wp
        nx       = size(LSF,1)
        ny       = size(LSF,2)
    
        do i = 1, nx
            do j = 1, ny
                ! ac-nodes
                call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
                dLSF_acx(i,j) = (LSF(ip1,j)-LSF(i,j))
                dLSF_acy(i,j) = (LSF(i,jp1)-LSF(i,j))
                qx_ac(i,j) = u(i,j) * dLSF_acx(i,j)
                qy_ac(i,j) = v(i,j) * dLSF_acy(i,j)
            end do
        end do
    
        do i = 1, nx
            do j = 1, ny
                ! Compute to aa-nodes
                call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)    
                qx_aa(i,j) = 0.5*(qx_ac(i,j) + qx_ac(im1,j)) 
                qy_aa(i,j) = 0.5*(qy_ac(i,j) + qy_ac(i,jm1))
            end do
        end do
            
        ! Update LSF
        dlsf = - dtdx * (qx_aa + qy_aa)
        LSF  = LSF + dlsf 
            
        if (.FALSE.) then
            ! Apply bounds
            where (LSF .lt. -1.0)
                LSF = -1.0
            end where
            
            where (LSF .gt. 1.0)
                LSF = 1.0
            end where
        end if

        return
    
    end subroutine LSFadvec_simple

    ! ===================================================================
    !
    ! Oceanic extrapolation routines.
    !
    ! ===================================================================
    
    subroutine extrapolate_ocn_acy(mask_fill,mask_orig,mask_ac)
        ! Routine to extrapolate the nearest value in the y-direction.
        ! So far we assume that value 0 in mask_orig represents ice-free points
    
        implicit none
        
        real(wp), intent(INOUT) :: mask_fill(:,:)
        real(wp), intent(IN)    :: mask_orig(:,:)
        real(wp), intent(IN)    :: mask_ac(:,:)
        
        ! Local variables
        integer  :: i, j
        integer :: count, repeat
        real(wp), allocatable :: mask_change(:,:),mask_change_n(:,:)
        allocate(mask_change(size(mask_orig,1),size(mask_orig,2)))
        allocate(mask_change_n(size(mask_orig,1),size(mask_orig,2)))
    
        ! Initialize variables
        mask_change   = 1.0_wp
        mask_change_n = 1.0_wp
    
        where(mask_ac .eq. 0.0_wp) mask_change = 0.0_wp
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
    
    end subroutine extrapolate_ocn_acy
    
    subroutine extrapolate_ocn_acx(mask_fill,mask_orig,mask_ac)
        ! Routine to extrapolate the nearest value in the x-direction.
        ! So far we assume that value 0 in mask_orig represents ice-free points
    
        implicit none
            
        real(wp), intent(INOUT) :: mask_fill(:,:)
        real(wp), intent(IN)    :: mask_orig(:,:)
        real(wp), intent(IN)    :: mask_ac(:,:)
    
        ! Local variables
        integer  :: i, j
        integer :: count, repeat
        real(wp), allocatable :: mask_change(:,:)
        allocate(mask_change(size(mask_orig,1),size(mask_orig,2)))
        
        ! Initialize variables
        mask_change = 1.0_wp
        where(mask_ac .eq. 0.0_wp) mask_change = 0.0_wp
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
    
    subroutine extrapolate_ocn_laplace_simple(mask_fill, mask_orig,mask_ac)
        ! Routine to extrapolate values using the Laplace equation.
        ! Assumes that value 0 in mask_ac represents ice-free points
            
        implicit none
        
        real(wp), intent(INOUT) :: mask_fill(:,:)
        real(wp), intent(IN) :: mask_orig(:,:)
        real(wp), intent(IN) :: mask_ac(:,:)
            
        ! Local variables
        integer :: i, j, iter
        real(wp) :: error, tol
        real(wp), allocatable :: mask_new(:,:)
            
        ! Allocate memory for the temporary array
        allocate(mask_new(size(mask_orig,1), size(mask_orig,2)))
            
        ! Initialize variables
        mask_fill = mask_orig
        mask_new  = mask_orig
        tol       = 1e-2_wp      ! Tolerance for convergence
        error     = tol + 1.0_wp
        iter      = 0
            
        ! Jacobi iteration
        do while (error > tol)
            error = 0.0_wp
            iter = iter + 1
            
            do i = 2, size(mask_orig,1)-1
                do j = 2, size(mask_orig,2)-1
                    if (mask_ac(i,j) .eq. 0.0_wp) then
                        mask_new(i,j) = 0.25_wp * (mask_fill(i+1,j) + mask_fill(i-1,j) + mask_fill(i,j+1) + mask_fill(i,j-1))
                        error = error + abs(mask_new(i,j) - mask_fill(i,j))
                    end if
                end do
            end do
            mask_fill = mask_new
        end do
            
        deallocate(mask_new)
    
        return
            
    end subroutine extrapolate_ocn_laplace_simple    
    
end module lsf_module
