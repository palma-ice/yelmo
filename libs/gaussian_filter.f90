!========================================================================
!
! Initial module from: 
!   https://github.com/nicjhan/gaussian-filter
!
! ajr: 
!    - sigma can be specified in the units of the grid
!    - mask is logical instead of real
!========================================================================
module gaussian_filter

use, intrinsic :: iso_fortran_env, only: error_unit

implicit none

interface filter_gaussian
    module procedure filter_gaussian_float
    module procedure filter_gaussian_dble 
end interface 

private

public :: filter_gaussian
public :: gaussian_kernel
public :: convolve 

public :: hgrad 
public :: diffuse2D, diff2D_timestep

contains

    subroutine filter_gaussian_dble(var,sigma,dx,mask,truncate)
        ! Wrapper for input as doubles
        real(kind=8), intent(inout) :: var(:,:)
        real(kind=8), intent(in)    :: sigma
        real(kind=8), intent(in), optional :: dx 
        logical,      intent(in), optional :: mask(:,:) 
        real(kind=4), intent(in), optional :: truncate

        real(kind=4), allocatable :: kernel(:,:)
        real(kind=4) :: sigmap 
        integer :: nx, ny 

        real(kind=4), allocatable :: output4(:,:)

        ! Allocate local output array 
        allocate(output4(size(var,1),size(var,2)))

        output4 = real(var) 

        call filter_gaussian_float(output4,real(sigma),real(dx),mask,truncate)
        
        ! Return a double array
        var = dble(output4)

        return 

    end subroutine filter_gaussian_dble

    subroutine filter_gaussian_float(var,sigma,dx,mask,truncate)

        real(kind=4), intent(inout) :: var(:,:)
        real(kind=4), intent(in)    :: sigma
        real(kind=4), intent(in), optional :: dx 
        logical,      intent(in), optional :: mask(:,:) 
        real(kind=4), intent(in), optional :: truncate

        real(kind=4), allocatable :: intmp(:,:), kernel(:,:), output(:,:) 
        real(kind=4) :: sigmap 
        integer :: nx, ny, nloop, i 

        ! Get sigma in terms of points
        sigmap = sigmap 
        if (present(dx)) then 
            sigmap = sigma / dx 
        end if 

        ! Get the kernel
        call gaussian_kernel(sigmap, kernel, truncate)
        
        nloop = 1
        if (size(kernel,1) > size(var,1)+1 .or. &
            size(kernel,2) > size(var,2)+1 ) then 
            nloop  = 4
            sigmap = sigmap/2.0
            call gaussian_kernel(sigmap, kernel, truncate)
        end if 

        ! Convolve as many times as necessary to acheive 
        ! desired level of smoothing 
        allocate(intmp(size(var,1),size(var,2)))
        allocate(output(size(var,1),size(var,2)))
        intmp = var 
        do i = 1, nloop
            call convolve(intmp, kernel, output, mask)
            intmp = output 
        end do 

        ! Store final solution 
        var = output 

        return 

    end subroutine filter_gaussian_float

    ! This is a copy of code found in gaussian_filter.py. These two implementations
    ! must remain equivalent for tests to pass.
    subroutine gaussian_kernel(sigma, kernel, truncate)

        real, intent(in) :: sigma
        real, intent(out), dimension(:,:), allocatable :: kernel
        real, intent(in), optional :: truncate

        real, dimension(:,:), allocatable :: x, y
        integer :: radius, i, j
        real :: trunc, s

        trunc = 4.0
        if (present(truncate)) trunc = truncate

        radius = ceiling(trunc * sigma + 0.5)
        s = sigma**2

        ! Set up meshgrid.
        allocate(x(-radius:radius, -radius:radius))
        allocate(y(-radius:radius, -radius:radius))
        do j = -radius, radius
            do i = -radius, radius
                x(i, j) = i
                y(i, j) = j
            enddo
        enddo

        ! Make kernel.
        if (allocated(kernel)) deallocate(kernel)
        allocate(kernel(-radius:radius, -radius:radius))
        kernel = 2.0*exp(-0.5 * (x**2 + y**2) / s)
        kernel = kernel / sum(kernel)

        deallocate(x)
        deallocate(y)

!         write(*,"(a,2f10.2,3i6)") "gaussian_kernel:: truncation, sigma, radius, nx, ny: ", &
!                                     trunc, sigma, radius, size(kernel,1), size(kernel,2) 
!         write(*,"(a,2g12.3)") "    kernel weight range: ", minval(kernel), maxval(kernel)

        return 

    end subroutine gaussian_kernel

    ! Set up 3x3 tiles around the input.
    subroutine tile_and_reflect(input, output)

        real, intent(in), dimension(:,:) :: input
        real, intent(out), dimension(:,:), allocatable :: output

        integer :: rows, cols

        rows = ubound(input, 1)
        cols = ubound(input, 2)

        ! Rely on automatic deallocation to clean this up. 
        allocate(output(3*rows, 3*cols))

        ! There are 3x3 tiles, we start at the top left and set the tiles up row by
        ! row.
        ! Top left is flipped left-to-right and up-to-down.
        output(:rows, :cols) = input(rows:1:-1, cols:1:-1)
        ! Top centre is flipped up-to-down
        output(:rows, cols+1:2*cols) = input(rows:1:-1, :)
        ! Top right is flipped left-to-right and up-to-down.
        output(:rows, 2*cols+1:3*cols) = input(rows:1:-1, cols:1:-1)
        ! Middle left flipped left-to-right
        output(rows+1:2*rows, :cols) = input(:, cols:1:-1)
        ! Middle centre unchanged
        output(rows+1:2*rows, cols+1:2*cols) = input(:, :)
        ! Middle right flipped left-to-right
        output(rows+1:2*rows, 2*cols+1:3*cols) = input(:, cols:1:-1)
        ! Bottom left flipped left-to-right and up-to-down
        output(2*rows+1:3*rows, :cols) = input(rows:1:-1, cols:1:-1)
        ! Bottom cente flipped up-to-down
        output(2*rows+1:3*rows, cols+1:2*cols) = input(rows:1:-1, :)
        ! Bottom right flipped left-to-right and up-to-down
        output(2*rows+1:3*rows, 2*cols+1:3*cols) = input(rows:1:-1, cols:1:-1)

    end subroutine tile_and_reflect

    ! Convolution.
    ! ajr: changed mask to be a logical array for easier usage
    subroutine convolve(input, weights, output, mask)

        real(4), intent(in), dimension(:,:) :: input, weights
        logical, intent(in), dimension(:,:), optional :: mask  ! .FALSE. points remain unchanged.
        real(4), intent(inout), dimension(:,:) :: output

        real(4), dimension(:,:), allocatable :: mask_real   ! Local mask for use with tiling routine

        ! These are allocated within tile_and_reflect, we rely on automatic
        ! deallocation at the end of the subroutine. 
        real, dimension(:, :), allocatable, target :: tiled_input
        real, dimension(:, :), allocatable, target :: tiled_mask
        real, dimension(:, :), pointer :: overlapping, overlapping_mask
        
        integer :: rows, cols, hw_row, hw_col, i, j, tj, ti
        real :: clobber_total, correction

        character(len=512) :: errmsg 

        ! First step is to tile the input.
        rows = ubound(input, 1)
        cols = ubound(input, 2)
        ! Stands for half weights row.
        hw_row = ubound(weights, 1) / 2
        hw_col = ubound(weights, 2) / 2

        ! Only one reflection is done on each side so the weights array cannot be
        ! bigger than width/height of input +1.
        errmsg = "convolve:: error: the chosen sigma is too large for the grid. &
                 &Only one reflection is done on each side of the grid, so the &
                 &size of the weights array (determined by sigma) cannot be &
                 &bigger than the grid's nx+1 or ny+1."
        call assert(ubound(weights, 1) < rows + 1,trim(errmsg))
        call assert(ubound(weights, 2) < cols + 1,trim(errmsg))

        if (present(mask)) then
            allocate(mask_real(size(mask,1),size(mask,2)))
            mask_real = 0.0
            where(mask) mask_real = 1.0
            call assert(all(shape(mask) - shape(input) == 0), &
                        'Mask and input shapes do not match')
            call tile_and_reflect(mask_real, tiled_mask)
        endif 

        ! This ensures that in the masked case, all masked points remain unchanged.
        output(:,:) = input(:,:)
        call tile_and_reflect(input, tiled_input)

        ! Very similar Python code can be found in gaussian_filter.py. 
        do j = 1, cols
            do i = 1, rows
                ! Use i, j to offset into equivalent part of the tiled arrays.
                ti = i + rows
                tj = j + cols

                ! Skip masked points. 
                if (present(mask)) then
                    if (tiled_mask(ti, tj) == 0) then
                        cycle
                    endif
                endif

                overlapping => tiled_input(ti - hw_row:ti + hw_row, &
                                           tj - hw_col:tj + hw_col)

                if (present(mask)) then
                    ! The approach taken here is to find out which parts of the
                    ! weights matrix are masked, add up the value of all these and
                    ! destribute them evenly over the rest of the matrix. The
                    ! intention is to conserve the field. 
                    overlapping_mask => tiled_mask(ti - hw_row:ti + hw_row, &
                                                   tj - hw_col:tj + hw_col)

                    ! Total value and number of weights clobbered by the mask.
                    clobber_total = sum((1 - overlapping_mask) * weights)
                    correction = clobber_total / sum(overlapping_mask)

                    ! Add correction and calculate. 
                    output(i, j) = sum((weights(:, :) + correction) * overlapping &
                                       * overlapping_mask)
                else
                    output(i, j) = sum(weights(:,:) * overlapping)
                endif
            enddo
        enddo

    end subroutine convolve

    subroutine assert(statement, msg)

        logical, intent(in) :: statement
        character(len=*), intent(in) :: msg

        if (.not. statement) then
            write(error_unit, *) msg
            stop 'Assert triggered, see stderr.'
        endif

    end subroutine assert


    ! Get the vector magnitude from two components
    elemental function calc_magnitude(u,v) result(umag)
        implicit none 
        real(4), intent(IN)  :: u, v 
        real(4) :: umag 

        umag = sqrt(u*u+v*v)

        return
    end function calc_magnitude 

    ! Horizontal gradient, x-component
    subroutine d_dx(du,u0,dx)  

        implicit none

        integer :: i, j, nx, ny
        real(4), dimension(:,:), intent(IN) :: u0
        real(4), dimension(:,:) :: du
        real(4) :: dx, inv_2dx

        inv_2dx = 1.d0 / (2.d0 * dx)
        nx = size(u0,1)
        ny = size(u0,2)

        du = 0.d0
        do i = 2, nx-1
            du(i,:) = (u0(i+1,:) - u0(i-1,:)) * inv_2dx
        end do

        ! Assign boundary values
        du(1,:)  = du(2,:)
        du(ny,:) = du(ny-1,:)

        return

    end subroutine d_dx

    ! Horizontal gradient, y-component
    subroutine d_dy(du,u0,dx)  

        implicit none

        integer :: i, j, nx, ny
        real(4), dimension(:,:), intent(IN) :: u0
        real(4), dimension(:,:) :: du
        real(4) :: dx, inv_2dx

        inv_2dx = 1.d0 / (2.d0 * dx)
        nx = size(u0,1)
        ny = size(u0,2)

        du = 0.d0
        do j = 2, ny-1
            du(:,j) = (u0(:,j+1) - u0(:,j-1)) * inv_2dx
        end do

        ! Assign boundary values
        du(:,1)  = du(:,2)
        du(:,ny) = du(:,ny-1)

        return

    end subroutine d_dy

    function hgrad(u0,dx,norm) result(du)

        implicit none

        integer :: i, j, nx, ny
        real(4), dimension(:,:), intent(IN) :: u0
        real(4), dimension(size(u0,1),size(u0,2)) :: du
        real(4), dimension(size(u0,1),size(u0,2)) :: dudx, dudy
        real(4) :: dx
        logical, optional :: norm 
        logical :: normalize 

        call d_dx(dudx,u0,dx=dx)
        call d_dy(dudy,u0,dx=dx)

        du = calc_magnitude(dudx,dudy)

        normalize = .FALSE.
        if (present(norm)) normalize = norm 
        if (normalize) then 
            du = du / maxval(abs(du))
        end if 

        return

    end function hgrad



    ! Explicit calculation of 2D Laplace: d2u/dxdy
    ! If input units are [u], returns [u/m2]
    subroutine diffuse2D(input,output,dx,mask)

        implicit none 

        real(4), dimension(:,:) :: input, output
        logical, intent(in), dimension(:,:), optional :: mask  ! .FALSE. points remain unchanged.
        real(4) :: dx, inv_dx2
        integer :: i, j, nx, ny  

        real(4), allocatable :: tmp(:,:)
        real(4), allocatable :: grad(:,:)
        logical, allocatable :: mask_local(:,:)
        real(4) :: dt, kappa  

        nx = size(input,1)
        ny = size(input,2) 

        allocate(mask_local(nx,ny))
        allocate(tmp(nx,ny),grad(nx,ny))

        ! Determine mask points
        mask_local = .TRUE. 
        if (present(mask)) mask_local = mask 

        inv_dx2 = 1.d0 / (dx*dx)

        ! Normalize input values 
        tmp = input !/ maxval(abs(input))

        ! Calculate gradient at each point
        grad = 0.0 
        do i = 2,nx-1
            do j = 2,ny-1

                if (mask_local(i,j)) then 

                    !  five-point stencil finite difference method
                    grad(i,j) = inv_dx2*(tmp(i-1,j)-2.d0*tmp(i,j)+tmp(i+1,j)) &
                              + inv_dx2*(tmp(i,j-1)-2.d0*tmp(i,j)+tmp(i,j+1))

                end if 

            end do 
        end do

!         grad(1,:)  = grad(2,:)
!         grad(nx,:) = grad(nx-1,:)
!         grad(:,1)  = grad(:,2)
!         grad(:,ny) = grad(:,ny-1)

!         output = input + dt*kappa*grad
        output = grad 

        return

    end subroutine diffuse2D 

    function diff2D_timestep(dx,dy,kappa) result(dt)
        implicit none 
        real(4) :: dx, dy, kappa 
        real(4) :: dt 

        ! 1D condition kappa*dt/dx^2 <= 1
        ! dt = dx^2/kappa 

        dt = (1.d0/(2.d0*kappa)) *(dx*dy)**2/(dx**2+dy**2)
        return 

    end function diff2D_timestep 
    
end module gaussian_filter
