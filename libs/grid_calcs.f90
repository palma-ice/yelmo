module grid_calcs

  ! Contains all the routines involved in the Arakawa grids, for mapping and
  ! calculating derivatives between all the different grids.

  use yelmo_defs, only : wp, dp, io_unit_err

  public 

contains

! ======================
! ==== Derivatives =====
! ======================
  
  ! Aa to Aa
  
  ! 2D
  subroutine ddx_a_to_a_2D(dx_a,d_a,dx)
    ! Input:  scalar on the Aa grid
    ! Output: its x-derivative on the Aa grid
    
    implicit none
    
    ! In/output variables:
    
    real(wp), intent(OUT)   :: dx_a(:,:)
    real(wp), intent(IN)    :: d_a(:,:)
    real(wp), intent(IN)    :: dx 

    ! Local variables:
    integer :: i, j, nx, ny 
    
    nx = size(d_a,1)
    ny = size(d_a,2) 

    ! Central differencing in the interior
    do j = 1, ny
    do i = 2, nx-1
        dx_a(i,j) = (d_a(i+1,j) - d_a(i-1,j)) / (2.0*dx)
    end do
    end do
    
    
    ! One-sided differencing on the boundaries
    dx_a(1,:)  = (d_a(2,:) - d_a(1,:)) / dx
    dx_a(nx,:) = (d_a(nx,:) - d_a(nx-1,:)) / dx
    
    return 

  end subroutine ddx_a_to_a_2D

  subroutine ddy_a_to_a_2D(dy_a,d_a,dx)
    ! Input:  scalar on the Aa grid
    ! Output: its y-derivative on the Aa grid
    
    implicit none
    
    ! In/output variables:
    
    real(wp), intent(OUT)   :: dy_a(:,:)
    real(wp), intent(IN)    :: d_a(:,:)
    real(wp), intent(IN)    :: dx 

    ! Local variables:
    integer :: i, j, nx, ny 
    
    nx = size(d_a,1)
    ny = size(d_a,2) 

    ! Central differencing in the interior
    do j = 2, ny-1
    do i = 1, nx
      dy_a(i,j) = (d_a(i,j+1) - d_a(i,j-1)) / (2.0*dx)
    end do
    end do
    
    ! One-sided differencing on the boundaries
    dy_a(:,1)  = (d_a(:,2) - d_a(:,1)) / dx
    dy_a(:,ny) = (d_a(:,ny) - d_a(:,ny-1)) / dx
    
    return 

  end subroutine ddy_a_to_a_2D

  subroutine ddxx_a_to_a_2D(dxx_a, d_a, dx)
    ! Input:  scalar on the Aa grid
    ! Output: its xx-derivative on the Aa grid
    
    implicit none
    
    ! In/output variables:
    
    real(wp), intent(OUT)   :: dxx_a(:,:)
    real(wp), intent(IN)    :: d_a(:,:)
    real(wp), intent(IN)    :: dx 
    
    ! Local variables:
    integer :: i, j, nx, ny
    
    nx = size(d_a,1)
    ny = size(d_a,2) 

    ! Central differencing in the interior
    do j = 1, ny
    do i = 2, nx-1
      dxx_a(i,j) = (d_a(i+1,j) + d_a(i-1,j) - 2.0_wp*d_a(i,j)) / dx**2
    end do
    end do
    
    
    ! One-sided differencing on the boundaries
    dxx_a(1,:)  = (d_a(3,:) + d_a(1,:) - 2.0_wp*d_a(2,:)) / dx**2
    dxx_a(nx,:) = (d_a(nx,:) + d_a(nx-2,:) - 2.0_wp * d_a(nx-1,:)) / dx**2
    
    return
    
  end subroutine ddxx_a_to_a_2D

  subroutine ddyy_a_to_a_2D(dyy_a, d_a, dx)
    ! Input:  scalar on the Aa grid
    ! Output: its yy-derivative on the Aa grid
    
    implicit none
    
    ! In/output variables:
    
    real(wp), intent(OUT)   :: dyy_a(:,:)
    real(wp), intent(IN)    :: d_a(:,:)
    real(wp), intent(IN)    :: dx 

    ! Local variables:
    integer :: i, j, nx, ny
    
    nx = size(d_a,1)
    ny = size(d_a,2) 

    ! Central differencing in the interior
    do i = 1, nx
    do j = 2, ny-1
      dyy_a(i,j) = (d_a(i,j+1) + d_a(i,j-1) - 2.0_wp * d_a(i,j)) / dx**2
    end do
    end do
    
    ! One-sided differencing on the boundaries
    dyy_a(:,1)  = (d_a(:,3) + d_a(:,1) - 2.0_wp*d_a(:,2)) / dx**2
    dyy_a(:,ny) = (d_a(:,ny) + d_a(:,ny-2) - 2.0_wp*d_a(:,ny-1)) / dx**2
    
    return

  end subroutine ddyy_a_to_a_2D

!   subroutine ddxy_a_to_a_2D( d_a, dxy_a)
!     ! Input:  scalar on the Aa grid
!     ! Output: its xy-derivative on the Aa grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), intent(IN)    :: d_a
!     real(wp), intent(OUT)   :: dxy_a
    
!     ! Local variables:
!     integer :: i, j, nx, ny
    
!     ! Central differencing in the interior
!     do i = 2, nx-1
!     do j = 2, ny-1
!       dxy_a(i,j) = (d_a( j+1,i+1) + d_a( j-1,i-1) - d_a( j+1,i-1) - d_a( j-1,i+1)) / (4.0_wp * dx * dx)
!     end do
!     end do
    
    
!     ! One-sided differencing on the boundaries
!     ! NO IDEA HOW TO do THIS...
!     dxy_a( 1              ,grid%i1:grid%i2) = 0.0_wp
!     dxy_a( grid%ny        ,grid%i1:grid%i2) = 0.0_wp
!     dxy_a( grid%j1:grid%j2,1              ) = 0.0_wp
!     dxy_a( grid%j1:grid%j2,grid%nx        ) = 0.0_wp
    
    
!   end subroutine ddxy_a_to_a_2D
!   ! 3D
!   subroutine ddx_a_to_a_3D( d_a, dx_a)
!     ! Input:  scalar on the Aa grid
!     ! Output: its x-derivative on the Aa grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), intent(IN)    :: d_a
!     real(wp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), intent(OUT)   :: dx_a
    
!     ! Local variables:
!     integer :: i, j, nx, ny,k
    
!     ! Central differencing in the interior
!     do i = 2, nx-1
!     do j = 1, ny
!     do k = 1, C%nZ
!       dx_a( k,j,i) = (d_a( k,j,i+1) - d_a( k,j,i-1)) / (2 * dx)
!     end do
!     end do
!     end do
    
    
!     ! One-sided differencing on the boundaries
!     dx_a( :,grid%j1:grid%j2,1      ) = (d_a( :,grid%j1:grid%j2,2      ) - d_a( :,grid%j1:grid%j2,1        )) / dx
!     dx_a( :,grid%j1:grid%j2,grid%nx) = (d_a( :,grid%j1:grid%j2,grid%nx) - d_a( :,grid%j1:grid%j2,grid%nx-1)) / dx
    
    
!   end subroutine ddx_a_to_a_3D
!   subroutine ddy_a_to_a_3D( d_a, dy_a)
!     ! Input:  scalar on the Aa grid
!     ! Output: its y-derivative on the Aa grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), intent(IN)    :: d_a
!     real(wp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), intent(OUT)   :: dy_a
    
!     ! Local variables:
!     integer :: i, j, nx, ny,k
    
!     ! Central differencing in the interior
!     do i = 1, nx
!     do j = 2, ny-1
!     do k = 1, C%nZ
!       dy_a( k,j,i) = (d_a( k,j+1,i) - d_a( k,j-1,i)) / (2 * dx)
!     end do
!     end do
!     end do
    
    
!     ! One-sided differencing on the boundaries
!     dy_a( :,1      ,grid%i1:grid%i2) = (d_a( :,2      ,grid%i1:grid%i2) - d_a( :,1        ,grid%i1:grid%i2)) / dx
!     dy_a( :,grid%ny,grid%i1:grid%i2) = (d_a( :,grid%ny,grid%i1:grid%i2) - d_a( :,grid%ny-1,grid%i1:grid%i2)) / dx
    
    
!   end subroutine ddy_a_to_a_3D
!   subroutine ddxx_a_to_a_3D( d_a, dxx_a)
!     ! Input:  scalar on the Aa grid
!     ! Output: its xx-derivative on the Aa grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), intent(IN)    :: d_a
!     real(wp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), intent(OUT)   :: dxx_a
    
!     ! Local variables:
!     integer :: i, j, nx, ny,k
    
!     ! Central differencing in the interior
!     do i = 2, nx-1
!     do j = 1, ny
!     do k = 1, C%nZ
!       dxx_a( k,j,i) = (d_a( k,j,i+1) + d_a( k,j,i-1) - 2.0_wp * d_a( k,j,i)) / dx**2
!     end do
!     end do
!     end do
    
    
!     ! One-sided differencing on the boundaries
!     dxx_a( :,grid%j1:grid%j2,1      ) = (d_a( :,grid%j1:grid%j2,3      ) + d_a( :,grid%j1:grid%j2,1        ) - 2.0_wp * d_a( :,grid%j1:grid%j2,2        )) / dx**2
!     dxx_a( :,grid%j1:grid%j2,grid%nx) = (d_a( :,grid%j1:grid%j2,grid%nx) + d_a( :,grid%j1:grid%j2,grid%nx-2) - 2.0_wp * d_a( :,grid%j1:grid%j2,grid%nx-1)) / dx**2
    
    
!   end subroutine ddxx_a_to_a_3D
!   subroutine ddyy_a_to_a_3D( d_a, dyy_a)
!     ! Input:  scalar on the Aa grid
!     ! Output: its yy-derivative on the Aa grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), intent(IN)    :: d_a
!     real(wp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), intent(OUT)   :: dyy_a
    
!     ! Local variables:
!     integer :: i, j, nx, ny,k
    
!     ! Central differencing in the interior
!     do i = 1, nx
!     do j = 2, ny-1
!     do k = 1, C%nZ
!       dyy_a( k,j,i) = (d_a( k,j+1,i) + d_a( k,j-1,i) - 2.0_wp * d_a( k,j,i)) / dx**2
!     end do
!     end do
!     end do
    
    
!     ! One-sided differencing on the boundaries
!     dyy_a( :,1      ,grid%i1:grid%i2) = (d_a( :,3      ,grid%i1:grid%i2) + d_a( :,1        ,grid%i1:grid%i2) - 2.0_wp * d_a( :,2        ,grid%i1:grid%i2)) / dx**2
!     dyy_a( :,grid%ny,grid%i1:grid%i2) = (d_a( :,grid%ny,grid%i1:grid%i2) + d_a( :,grid%ny-2,grid%i1:grid%i2) - 2.0_wp * d_a( :,grid%ny-1,grid%i1:grid%i2)) / dx**2
    
    
!   end subroutine ddyy_a_to_a_3D
!   subroutine ddxy_a_to_a_3D( d_a, dxy_a)
!     ! Input:  scalar on the Aa grid
!     ! Output: its xy-derivative on the Aa grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), intent(IN)    :: d_a
!     real(wp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), intent(OUT)   :: dxy_a
    
!     ! Local variables:
!     integer :: i, j, nx, ny,k
    
!     ! Central differencing in the interior
!     do i = 2, nx-1
!     do j = 2, ny-1
!     do k = 1, C%nZ
!       dxy_a( k,j,i) = (d_a( k,j+1,i+1) + d_a( k,j-1,i-1) - d_a( k,j+1,i-1) - d_a( k,j-1,i+1)) / (4.0_wp * dx * dx)
!     end do
!     end do
!     end do
    
    
!     ! One-sided differencing on the boundaries
!     ! NO IDEA HOW TO do THIS...
!     dxy_a( :,1              ,grid%i1:grid%i2) = 0.0_wp
!     dxy_a( :,grid%ny        ,grid%i1:grid%i2) = 0.0_wp
!     dxy_a( :,grid%j1:grid%j2,1              ) = 0.0_wp
!     dxy_a( :,grid%j1:grid%j2,grid%nx        ) = 0.0_wp
    
    
!   end subroutine ddxy_a_to_a_3D
!   ! 3D upwind, for thermodynamics
!   subroutine ddx_a_to_a_3D_upwind( d_a, dx_a, U_3D_a)
!     ! Input:  scalar on the Aa grid
!     ! Output: its x-derivative on the Aa grid, using upwind one-sided differencing
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), intent(IN)    :: d_a
!     real(wp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), intent(OUT)   :: dx_a
!     real(wp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), intent(IN)    :: U_3D_a
    
!     ! Local variables:
!     integer :: i, j, nx, ny,k
    
!     ! Upwind one-sided differencing
!     do i = 2, nx-1
!     do j = 1, ny
!     do k = 1, C%nZ
!       IF (U_3D_a( k,j,i) > 0.0_wp) THEN
!         dx_a( k,j,i) = (d_a( k,j,i  ) - d_a( k,j,i-1)) / dx
!       ELSE
!         dx_a( k,j,i) = (d_a( k,j,i+1) - d_a( k,j,i  )) / dx
!       end IF
!     end do
!     end do
!     end do
    
    
!     dx_a( :,grid%j1:grid%j2,1      ) = 0.0_wp
!     dx_a( :,grid%j1:grid%j2,grid%nx) = 0.0_wp
    
    
!   end subroutine ddx_a_to_a_3D_upwind
!   subroutine ddy_a_to_a_3D_upwind( d_a, dy_a, V_3D_a)
!     ! Input:  scalar on the Aa grid
!     ! Output: its y-derivative on the Aa grid, using upwind one-sided differencing
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), intent(IN)    :: d_a
!     real(wp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), intent(OUT)   :: dy_a
!     real(wp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), intent(IN)    :: V_3D_a
    
!     ! Local variables:
!     integer :: i, j, nx, ny,k
    
!     ! Upwind one-sided differencing
!     do i = 1, nx
!     do j = 2, ny-1
!     do k = 1, C%nZ
!       IF (V_3D_a( k,j,i) > 0.0_wp) THEN
!         dy_a( k,j,i) = (d_a( k,j  ,i) - d_a( k,j-1,i)) / dx
!       ELSE
!         dy_a( k,j,i) = (d_a( k,j+1,i) - d_a( k,j  ,i)) / dx
!       end IF
!     end do
!     end do
!     end do
    
    
!     dy_a( :,1      ,grid%i1:grid%i2) = 0.0_wp
!     dy_a( :,grid%ny,grid%i1:grid%i2) = 0.0_wp
    
    
!   end subroutine ddy_a_to_a_3D_upwind
  
!   ! Aa to Acx/Acy
  
!   ! 2D
!   subroutine ddx_a_to_cx_2D( d_a, dx_cx)
!     ! Input:  scalar on the Aa grid
!     ! Output: its x-derivative on the Acx grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), intent(IN)    :: d_a
!     real(wp), DIMENSION(       grid%ny  , grid%nx-1), intent(OUT)   :: dx_cx
    
!     ! Local variables:
!     integer :: i, j, nx, ny
    
!     do i = grid%i1, MIN(grid%nx-1,grid%i2)
!     do j = 1, ny
!       dx_cx(i,j) = (d_a(i+1,j) - d_a(i,j)) / dx
!     end do
!     end do
    
    
!   end subroutine ddx_a_to_cx_2D
!   subroutine ddy_a_to_cy_2D( d_a, dy_cy)
!     ! Input:  scalar on the Aa grid
!     ! Output: its y-derivative on the Acy grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), intent(IN)    :: d_a
!     real(wp), DIMENSION(       grid%ny-1, grid%nx  ), intent(OUT)   :: dy_cy
    
!     ! Local variables:
!     integer :: i, j, nx, ny
    
!     do i = 1, nx
!     do j = 1, ny-1
!       dy_cy(i,j) = (d_a(i,j+1) - d_a(i,j)) / dx
!     end do
!     end do
    
    
!   end subroutine ddy_a_to_cy_2D
!   subroutine ddx_a_to_cy_2D( d_a, dx_cy)
!     ! Input:  scalar on the Aa grid
!     ! Output: its x-derivative on the Acy grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), intent(IN)    :: d_a
!     real(wp), DIMENSION(       grid%ny-1, grid%nx  ), intent(OUT)   :: dx_cy
    
!     ! Local variables:
!     integer :: i, j, nx, ny
    
!     ! Central differencing in the interior
!     do i = 2, nx-1
!     do j = 1, ny-1
!       dx_cy(i,j) = (d_a(i+1,j) + d_a( j+1,i+1) - d_a(i-1,j) - d_a( j+1,i-1)) / (4.0_wp * dx)
!     end do
!     end do
    
    
!     ! One-sided differencing on the boundary
!     do j = grid%j1, MIN(grid%ny-1,grid%j2)
!       dx_cy( j,1      ) = (d_a( j,2      ) + d_a( j+1,2      ) - d_a( j,1        ) - d_a( j+1,1        )) / (2.0_wp * dx)
!       dx_cy( j,grid%nx) = (d_a( j,grid%nx) + d_a( j+1,grid%nx) - d_a( j,grid%nx-1) - d_a( j+1,grid%nx-1)) / (2.0_wp * dx)
!     end do
    
    
!   end subroutine ddx_a_to_cy_2D
!   subroutine ddy_a_to_cx_2D( d_a, dy_cx)
!     ! Input:  scalar on the Aa grid
!     ! Output: its y-derivative on the Acx grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), intent(IN)    :: d_a
!     real(wp), DIMENSION(       grid%ny  , grid%nx-1), intent(OUT)   :: dy_cx
    
!     ! Local variables:
!     integer :: i, j, nx, ny
    
!     ! Central differencing in the interior
!     do i = grid%i1, MIN(grid%nx-1,grid%i2)
!     do j = 2, ny-1
!       dy_cx(i,j) = (d_a(i,j+1) + d_a( j+1,i+1) - d_a(i,j-1) - d_a( j-1,i+1)) / (4.0_wp * dx)
!     end do
!     end do
    
    
!     ! One-sided differencing on the boundary
!     do i = grid%i1, MIN(grid%nx-1,grid%i2)
!       dy_cx( 1      ,i) = (d_a( 2,      i) + d_a( 2,      i+1) - d_a( 1,        i) - d_a( 1,        i+1)) / (2.0_wp * dx)
!       dy_cx( grid%ny,i) = (d_a( grid%ny,i) + d_a( grid%ny,i+1) - d_a( grid%ny-1,i) - d_a( grid%ny-1,i+1)) / (2.0_wp * dx)
!     end do
    
    
!   end subroutine ddy_a_to_cx_2D
!   ! 3D
!   subroutine ddx_a_to_cx_3D( d_a, dx_cx)
!     ! Input:  scalar on the Aa grid
!     ! Output: its x-derivative on the Acx grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), intent(IN)    :: d_a
!     real(wp), DIMENSION( C%nZ, grid%ny  , grid%nx-1), intent(OUT)   :: dx_cx
    
!     ! Local variables:
!     integer :: i, j, nx, ny,k
    
!     do i = grid%i1, MIN(grid%nx-1,grid%i2)
!     do j = 1, ny
!     do k = 1, C%nZ
!       dx_cx( k,j,i) = (d_a( k,j,i+1) - d_a( k,j,i)) / dx
!     end do
!     end do
!     end do
    
    
!   end subroutine ddx_a_to_cx_3D
!   subroutine ddy_a_to_cy_3D( d_a, dy_cy)
!     ! Input:  scalar on the Aa grid
!     ! Output: its y-derivative on the Acy grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), intent(IN)    :: d_a
!     real(wp), DIMENSION( C%nZ, grid%ny-1, grid%nx  ), intent(OUT)   :: dy_cy
    
!     ! Local variables:
!     integer :: i, j, nx, ny,k
    
!     do i = 1, nx
!     do j = 1, ny-1
!     do k = 1, C%nZ
!       dy_cy( k,j,i) = (d_a( k,j+1,i) - d_a( k,j,i)) / dx
!     end do
!     end do
!     end do
    
    
!   end subroutine ddy_a_to_cy_3D
!   subroutine ddx_a_to_cy_3D( d_a, dx_cy)
!     ! Input:  scalar on the Aa grid
!     ! Output: its x-derivative on the Acy grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), intent(IN)    :: d_a
!     real(wp), DIMENSION( C%nZ, grid%ny-1, grid%nx  ), intent(OUT)   :: dx_cy
    
!     ! Local variables:
!     integer :: i, j, nx, ny,k
    
!     ! Central differencing in the interior
!     do i = 2, nx-1
!     do j = 1, ny-1
!     do k = 1, C%nZ
!       dx_cy( k,j,i) = (d_a( k,j,i+1) + d_a( k,j+1,i+1) - d_a( k,j,i-1) - d_a( k,j+1,i-1)) / (4.0_wp * dx)
!     end do
!     end do
!     end do
    
    
!     ! One-sided differencing on the boundary
!     do j = grid%j1, MIN(grid%ny-1,grid%j2)
!       dx_cy( :,j,1      ) = (d_a( :,j,2      ) + d_a( :,j+1,2      ) - d_a( :,j,1        ) - d_a( :,j+1,1        )) / (2.0_wp * dx)
!       dx_cy( :,j,grid%nx) = (d_a( :,j,grid%nx) + d_a( :,j+1,grid%nx) - d_a( :,j,grid%nx-1) - d_a( :,j+1,grid%nx-1)) / (2.0_wp * dx)
!     end do
    
    
!   end subroutine ddx_a_to_cy_3D
!   subroutine ddy_a_to_cx_3D( d_a, dy_cx)
!     ! Input:  scalar on the Aa grid
!     ! Output: its y-derivative on the Acx grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), intent(IN)    :: d_a
!     real(wp), DIMENSION( C%nZ, grid%ny  , grid%nx-1), intent(OUT)   :: dy_cx
    
!     ! Local variables:
!     integer :: i, j, nx, ny,k
    
!     ! Central differencing in the interior
!     do i = grid%i1, MIN(grid%nx-1,grid%i2)
!     do j = 2, ny-1
!     do k = 1, C%nZ
!       dy_cx( k,j,i) = (d_a( k,j+1,i) + d_a( k,j+1,i+1) - d_a( k,j-1,i) - d_a( k,j-1,i+1)) / (4.0_wp * dx)
!     end do
!     end do
!     end do
    
    
!     ! One-sided differencing on the boundary
!     do i = grid%i1, MIN(grid%nx-1,grid%i2)
!       dy_cx( :,1      ,i) = (d_a( :,2,      i) + d_a( :,2,      i+1) - d_a( :,1,        i) - d_a( :,1,        i+1)) / (2.0_wp * dx)
!       dy_cx( :,grid%nx,i) = (d_a( :,grid%nx,i) + d_a( :,grid%nx,i+1) - d_a( :,grid%nx-1,i) - d_a( :,grid%nx-1,i+1)) / (2.0_wp * dx)
!     end do
    
    
!   end subroutine ddy_a_to_cx_3D
  
!   ! Acx/Acy to Aa
  
!   ! 2D
!   subroutine ddx_cx_to_a_2D( d_cx, dx_a)
!     ! Input:  scalar on the Acx grid
!     ! Output: its x-derivative on the Aa grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION(       grid%ny  , grid%nx-1), intent(IN)    :: d_cx
!     real(wp), intent(OUT)   :: dx_a
    
!     ! Local variables:
!     integer :: i, j, nx, ny
    
!     do i = 2, nx-1
!     do j = 1, ny
!       dx_a(i,j) = (d_cx(i,j) - d_cx(i-1,j)) / dx
!     end do
!     end do
    
    
!     dx_a(1,:) = dx_a( grid%j1:grid%j2,2        )
!     dx_a(nx,:) = dx_a( grid%j1:grid%j2,grid%nx-1)
    
    
!   end subroutine ddx_cx_to_a_2D
!   subroutine ddy_cy_to_a_2D( d_cy, dy_a)
!     ! Input:  scalar on the Acy grid
!     ! Output: its y-derivative on the Aa grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION(       grid%ny-1, grid%nx  ), intent(IN)    :: d_cy
!     real(wp), intent(OUT)   :: dy_a
    
!     ! Local variables:
!     integer :: i, j, nx, ny
    
!     do i = 1, nx
!     do j = 2, ny-1
!       dy_a(i,j) = (d_cy(i,j) - d_cy(i,j-1)) / dx
!     end do
!     end do
    
    
!     dy_a(:,1) = dy_a( 2        ,grid%i1:grid%i2)
!     dy_a(:,ny) = dy_a( grid%ny-1,grid%i1:grid%i2)
    
    
!   end subroutine ddy_cy_to_a_2D
!   subroutine ddy_cx_to_a_2D( d_cx, dy_a)
!     ! Input:  scalar on the Acx grid
!     ! Output: its y-derivative on the Aa grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION(       grid%ny  , grid%nx-1), intent(IN)    :: d_cx
!     real(wp), intent(OUT)   :: dy_a
    
!     ! Local variables:
!     integer :: i, j, nx, ny
    
!     ! Interior
!     do i = 2, nx-1
!     do j = 2, ny-1
!       dy_a(i,j) = (d_cx( j+1,i-1) + d_cx(i,j+1) - d_cx( j-1,i-1) - d_cx(i,j-1)) / (4.0_wp * dx)
!     end do
!     end do
    
    
!     do i = 2, nx-1
!       ! South ex. corners
!       j = 1
!       dy_a(i,j) = (d_cx( j+1,i-1) + d_cx(i,j+1) - d_cx( j  ,i-1) - d_cx( j  ,i)) / (4.0_wp * dx)
!       ! North ex. corners
!       j = grid%ny
!       dy_a(i,j) = (d_cx( j  ,i-1) + d_cx( j  ,i) - d_cx( j-1,i-1) - d_cx(i,j-1)) / (4.0_wp * dx)
!     end do
    
    
!     do j = MAX(2,grid%j1), MIN(grid%ny-1,grid%j2)
!       ! West ex. corners
!       i = 1
!       dy_a(i,j) = (d_cx( j+1,i  ) - d_cx( j-1,i  )) / (2.0_wp * dx)
!       ! East ex. corners
!       i = grid%nx
!       dy_a(i,j) = (d_cx( j+1,i-1) - d_cx( j-1,i-1)) / (2.0_wp * dx)
!     end do
    
    
!     ! Corners
!     IF (par%master) THEN
!     dy_a( 1      ,1      ) = (d_cx( 2      ,1        ) - d_cx( 1        ,1        )) / dx
!     dy_a( 1      ,grid%nx) = (d_cx( 2      ,grid%nx-1) - d_cx( 1        ,grid%nx-1)) / dx
!     dy_a( grid%ny,1      ) = (d_cx( grid%ny,1        ) - d_cx( grid%ny-1,1        )) / dx
!     dy_a( grid%ny,grid%nx) = (d_cx( grid%ny,grid%nx-1) - d_cx( grid%ny-1,grid%nx-1)) / dx
!     end IF
    
    
!   end subroutine ddy_cx_to_a_2D
!   subroutine ddx_cy_to_a_2D( d_cy, dx_a)
!     ! Input:  scalar on the Acy grid
!     ! Output: its x-derivative on the Aa grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION(       grid%ny-1, grid%nx  ), intent(IN)    :: d_cy
!     real(wp), intent(OUT)   :: dx_a
    
!     ! Local variables:
!     integer :: i, j, nx, ny
    
!     ! Interior
!     do i = 2, nx-1
!     do j = 2, ny-1
!       dx_a(i,j) = (d_cy( j-1,i+1) + d_cy(i+1,j) - d_cy( j-1,i-1) - d_cy(i-1,j)) / (4.0_wp * dx)
!     end do
!     end do
    
    
!     do j = MAX(2,grid%j1), MIN(grid%ny-1,grid%j2)
!       ! West ex. corners
!       i = 1
!       dx_a(i,j) = (d_cy( j-1,i+1) + d_cy( j  ,i+1) - d_cy( j-1,i  ) - d_cy( j  ,i  )) / (4.0_wp * dx)
!       ! East ex. corners
!       i = grid%nx
!       dx_a(i,j) = (d_cy( j-1,i  ) + d_cy( j  ,i  ) - d_cy( j-1,i-1) - d_cy( j  ,i-1)) / (4.0_wp * dx)
!     end do
    
    
!     do i = 2, nx-1
!       ! South ex. corners
!       j = 1
!       dx_a(i,j) = (d_cy( j  ,i+1) - d_cy( j  ,i-1)) / (2.0_wp * dx)
!       ! North ex. corners
!       j = grid%ny
!       dx_a(i,j) = (d_cy( j-1,i+1) - d_cy( j-1,i-1)) / (2.0_wp * dx)
!     end do
    
    
!     ! Corners
!     IF (par%master) THEN
!     dx_a( 1      ,      1) = (d_cy( 1        ,2      ) - d_cy( 1        ,1        )) / dx
!     dx_a( 1      ,grid%nx) = (d_cy( 1        ,grid%nx) - d_cy( 1        ,grid%nx-1)) / dx
!     dx_a( grid%ny,1      ) = (d_cy( grid%ny-1,2      ) - d_cy( grid%ny-1,1        )) / dx
!     dx_a( grid%ny,grid%nx) = (d_cy( grid%ny-1,grid%nx) - d_cy( grid%ny-1,grid%nx-1)) / dx
!     end IF
    
    
!   end subroutine ddx_cy_to_a_2D
!   ! 3D
!   subroutine ddx_cx_to_a_3D( d_cx, dx_a)
!     ! Input:  scalar on the Acx grid
!     ! Output: its x-derivative on the Aa grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION( C%nZ, grid%ny  , grid%nx-1), intent(IN)    :: d_cx
!     real(wp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), intent(OUT)   :: dx_a
    
!     ! Local variables:
!     integer :: i, j, nx, ny,k
    
!     do i = 2, nx-1
!     do j = 1, ny
!     do k = 1, C%nZ
!       dx_a( k,j,i) = (d_cx( k,j,i) - d_cx( k,j,i-1)) / dx
!     end do
!     end do
!     end do
    
    
!   end subroutine ddx_cx_to_a_3D
!   subroutine ddy_cy_to_a_3D( d_cy, dy_a)
!     ! Input:  scalar on the Acy grid
!     ! Output: its y-derivative on the Aa grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION( C%nZ, grid%ny-1, grid%nx  ), intent(IN)    :: d_cy
!     real(wp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), intent(OUT)   :: dy_a
    
!     ! Local variables:
!     integer :: i, j, nx, ny,k
    
!     do i = 1, nx
!     do j = 2, ny-1
!     do k = 1, C%nZ
!       dy_a( k,j,i) = (d_cy( k,j,i) - d_cy( k,j-1,i)) / dx
!     end do
!     end do
!     end do
    
    
!   end subroutine ddy_cy_to_a_3D
  
!   ! Acx/Acy to Acx/Acy
!   subroutine ddx_cx_to_cx_2D( d_cx, dx_cx)
!     ! Input:  scalar on the Acx grid
!     ! Output: its x-derivative on the Acx grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION(       grid%ny  , grid%nx-1), intent(IN)    :: d_cx
!     real(wp), DIMENSION(       grid%ny  , grid%nx-1), intent(OUT)   :: dx_cx
    
!     ! Local variables:
!     integer :: i, j, nx, ny
    
!     ! Central differencing in the interior
!     do i = 2, grid%nx-2
!     do j = 1, ny
!       dx_cx(i,j) = (d_cx(i+1,j) - d_cx(i-1,j)) / (2 * dx)
!     end do
!     end do
    
    
!     ! One-sided differencing on the boundaries
!     dx_cx( grid%j1:grid%j2,1        ) = (d_cx( grid%j1:grid%j2,2        ) - d_cx( grid%j1:grid%j2,1        )) / dx
!     dx_cx( grid%j1:grid%j2,grid%nx-1) = (d_cx( grid%j1:grid%j2,grid%nx-1) - d_cx( grid%j1:grid%j2,grid%nx-2)) / dx
    
    
!   end subroutine ddx_cx_to_cx_2D
!   subroutine ddy_cx_to_cx_2D( d_cx, dy_cx)
!     ! Input:  scalar on the Acx grid
!     ! Output: its y-derivative on the Acx grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION(       grid%ny  , grid%nx-1), intent(IN)    :: d_cx
!     real(wp), DIMENSION(       grid%ny  , grid%nx-1), intent(OUT)   :: dy_cx
    
!     ! Local variables:
!     integer :: i, j, nx, ny
    
!     ! Central differencing in the interior
!     do i = grid%i1, MIN(grid%nx-1,grid%i2)
!     do j = 2, ny-1
!       dy_cx(i,j) = (d_cx(i,j+1) - d_cx(i,j-1)) / (2 * dx)
!     end do
!     end do
    
    
!     ! One-sided differencing on the boundaries
!     dy_cx(:,1) = (d_cx( 2      ,grid%i1:grid%i2) - d_cx( 1        ,grid%i1:grid%i2)) / dx
!     dy_cx(:,ny) = (d_cx(:,ny) - d_cx( grid%ny-1,grid%i1:grid%i2)) / dx
    
    
!   end subroutine ddy_cx_to_cx_2D
!   subroutine ddx_cy_to_cy_2D( d_cy, dx_cy)
!     ! Input:  scalar on the Acy grid
!     ! Output: its x-derivative on the Acy grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION(       grid%ny-1, grid%nx  ), intent(IN)    :: d_cy
!     real(wp), DIMENSION(       grid%ny-1, grid%nx  ), intent(OUT)   :: dx_cy
    
!     ! Local variables:
!     integer :: i, j, nx, ny
    
!     ! Central differencing in the interior
!     do i = 2, nx-1
!     do j = 1, ny-1
!       dx_cy(i,j) = (d_cy(i+1,j) - d_cy(i-1,j)) / (2 * dx)
!     end do
!     end do
    
    
!     ! One-sided differencing on the boundaries
!     dx_cy( grid%j1:MIN(grid%j2,grid%ny-1),1      ) = (d_cy( grid%j1:MIN(grid%j2,grid%ny-1),2      ) - d_cy( grid%j1:MIN(grid%j2,grid%ny-1),1        )) / dx
!     dx_cy( grid%j1:MIN(grid%j2,grid%ny-1),grid%nx) = (d_cy( grid%j1:MIN(grid%j2,grid%ny-1),grid%nx) - d_cy( grid%j1:MIN(grid%j2,grid%ny-1),grid%nx-1)) / dx
    
    
!   end subroutine ddx_cy_to_cy_2D
!   subroutine ddy_cy_to_cy_2D( d_cy, dy_cy)
!     ! Input:  scalar on the Acy grid
!     ! Output: its y-derivative on the Acy grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION(       grid%ny-1, grid%nx  ), intent(IN)    :: d_cy
!     real(wp), DIMENSION(       grid%ny-1, grid%nx  ), intent(OUT)   :: dy_cy
    
!     ! Local variables:
!     integer :: i, j, nx, ny
    
!     ! Central differencing in the interior
!     do i = 1, nx
!     do j = 2, grid%ny-2
!       dy_cy(i,j) = (d_cy(i,j+1) - d_cy(i,j-1)) / (2 * dx)
!     end do
!     end do
    
    
!     ! One-sided differencing on the boundaries
!     dy_cy( 1        ,grid%i1:grid%i2) = (d_cy( 2        ,grid%i1:grid%i2) - d_cy( 1        ,grid%i1:grid%i2)) / dx
!     dy_cy( grid%ny-1,grid%i1:grid%i2) = (d_cy( grid%ny-1,grid%i1:grid%i2) - d_cy( grid%ny-2,grid%i1:grid%i2)) / dx
    
    
!   end subroutine ddy_cy_to_cy_2D
!   subroutine ddx_cx_to_cy_2D( d_cx, dx_cy)
!     ! Input:  scalar on the Acx grid
!     ! Output: its x-derivative on the Acy grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION(       grid%ny  , grid%nx-1), intent(IN)    :: d_cx
!     real(wp), DIMENSION(       grid%ny-1, grid%nx  ), intent(OUT)   :: dx_cy
    
!     ! Local variables:
!     integer :: i, j, nx, ny
    
!     ! Interior
!     do i = 2, nx-1
!     do j = 1, ny-1
!       dx_cy(i,j) = (d_cx( j,i  ) + d_cx( j+1,i  ) - d_cx(i-1,j) - d_cx( j+1,i-1)) / (2.0_wp * dx)
!     end do
!     end do
    
    
!     ! Boundaries
!     do j = grid%j1, MIN(grid%ny-1,grid%j2)
!       ! West
!       dx_cy( j,1      ) = dx_cy( j,2        )
!       ! East
!       dx_cy( j,grid%nx) = dx_cy( j,grid%nx-1)
!     end do
    
    
!   end subroutine ddx_cx_to_cy_2D
!   subroutine ddy_cx_to_cy_2D( d_cx, dy_cy)
!     ! Input:  scalar on the Acx grid
!     ! Output: its y-derivative on the Acy grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION(       grid%ny  , grid%nx-1), intent(IN)    :: d_cx
!     real(wp), DIMENSION(       grid%ny-1, grid%nx  ), intent(OUT)   :: dy_cy
    
!     ! Local variables:
!     integer :: i, j, nx, ny
    
!     ! Interior
!     do i = 2, nx-1
!     do j = 1, ny-1
!       dy_cy(i,j) = (d_cx( j+1,i-1) + d_cx( j+1,i  ) - d_cx( j  ,i-1) - d_cx( j  ,i  )) / (2.0_wp * dx)
!     end do
!     end do
    
    
!     ! Boundaries
!     do j = grid%j1, MIN(grid%ny-1,grid%j2)
!       ! West
!       i = 1
!       dy_cy(i,j) = (d_cx( j,i  ) - d_cx( j,i  )) / dx
!       ! East
!       i = grid%nx
!       dy_cy(i,j) = (d_cx(i-1,j) - d_cx(i-1,j)) / dx
!     end do
    
    
!   end subroutine ddy_cx_to_cy_2D
!   subroutine ddx_cy_to_cx_2D( d_cy, dx_cx)
!     ! Input:  scalar on the Acy grid
!     ! Output: its x-derivative on the Acx grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION(       grid%ny-1, grid%nx  ), intent(IN)    :: d_cy
!     real(wp), DIMENSION(       grid%ny  , grid%nx-1), intent(OUT)   :: dx_cx
    
!     ! Local variables:
!     integer :: i, j, nx, ny
    
!     ! Interior
!     do i = grid%i1, MIN(grid%nx-1,grid%i2)
!     do j = 2, ny-1
!       dx_cx(i,j) = (d_cy( j-1,i+1) + d_cy( j  ,i+1) - d_cy( j-1,i  ) - d_cy( j  ,i  )) / (2.0_wp * dx)
!     end do
!     end do
    
    
!     ! Boundaries
!     do i = grid%i1, MIN(grid%nx-1,grid%i2)
!       ! South
!       j = 1
!       dx_cx(i,j) = (d_cy( j  ,i+1) - d_cy( j  ,i  )) / dx
!       ! North
!       j = grid%ny
!       dx_cx(i,j) = (d_cy( j-1,i+1) - d_cy( j-1,i  )) / dx
!     end do
    
    
!   end subroutine ddx_cy_to_cx_2D
!   subroutine ddy_cy_to_cx_2D( d_cy, dy_cx)
!     ! Input:  scalar on the Acy grid
!     ! Output: its y-derivative on the Acx grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION(       grid%ny-1, grid%nx  ), intent(IN)    :: d_cy
!     real(wp), DIMENSION(       grid%ny  , grid%nx-1), intent(OUT)   :: dy_cx
    
!     ! Local variables:
!     integer :: i, j, nx, ny
    
!     ! Interior
!     do i = grid%i1, MIN(grid%nx-1,grid%i2)
!     do j = 2, ny-1
!       dy_cx(i,j) = (d_cy( j  ,i  ) + d_cy( j  ,i+1) - d_cy( j-1,i  ) - d_cy( j-1,i+1)) / (2.0_wp * dx)
!     end do
!     end do
    
    
!     ! Boundaries
!     do i = grid%i1, MIN(grid%nx-1,grid%i2)
!       ! South
!       j = 1
!       dy_cx(i,j) = (d_cy( j+1,i  ) + d_cy(j+1,i+1) - d_cy( j  ,i  ) - d_cy( j  ,i+1)) / (2.0_wp * dx)
!       ! North
!       j = grid%ny
!       dy_cx(i,j) = (d_cy( j-1,i  ) + d_cy(j-1,i+1) - d_cy( j-2,i  ) - d_cy( j-2,i+1)) / (2.0_wp * dx)
!     end do
    
    
!   end subroutine ddy_cy_to_cx_2D
  
!   ! Acx to Ab
!   subroutine ddx_cx_to_b_2D( d_cx, dx_b)
!     ! Input:  scalar on the Acx grid
!     ! Output: its x-derivative on the Ab grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION(       grid%ny  , grid%nx-1), intent(IN)    :: d_cx
!     real(wp), DIMENSION(       grid%ny-1, grid%nx-1), intent(OUT)   :: dx_b
    
!     ! Local variables:
!     integer :: i, j, nx, ny
    
!     ! Interior
!     do i = MAX(2,grid%i1), MIN(grid%nx-2,grid%i2)
!     do j = 1, ny-1
!       dx_b(i,j) = (d_cx( j+1,i+1) + d_cx( j  ,i+1) - d_cx( j+1,i-1) - d_cx( j  ,i-1)) / (4.0_wp * dx)
!     end do
!     end do
    
    
!     ! Boundaries
!     do j = grid%j1, MIN(grid%ny-1,grid%j2)
!       i = 1
!       dx_b(i,j) = (d_cx( j+1,i+1) + d_cx( j  ,i+1) - d_cx( j+1,i  ) - d_cx( j  ,i  )) / (2.0_wp * dx)
!       i = grid%nx-1
!       dx_b(i,j) = (d_cx( j+1,i  ) + d_cx( j  ,i  ) - d_cx( j+1,i-1) - d_cx( j  ,i-1)) / (2.0_wp * dx)
!     end do
    
    
!   end subroutine ddx_cx_to_b_2D
!   subroutine ddy_cy_to_b_2D( d_cy, dy_b)
!     ! Input:  scalar on the Acy grid
!     ! Output: its y-derivative on the Ab grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION(       grid%ny-1, grid%nx  ), intent(IN)    :: d_cy
!     real(wp), DIMENSION(       grid%ny-1, grid%nx-1), intent(OUT)   :: dy_b
    
!     ! Local variables:
!     integer :: i, j, nx, ny
    
!     ! Interior
!     do i = MAX(2,grid%i1), MIN(grid%nx-2,grid%i2)
!     do j = 2, grid%ny-2
!       dy_b(i,j) = (d_cy( j+1,i+1) + d_cy( j+1,i  ) - d_cy( j-1,i+1) - d_cy( j-1,i  )) / (4.0_wp * dx)
!     end do
!     end do
    
    
!     ! Boundaries
!     do i = grid%i1, MIN(grid%nx-1,grid%i2)
!       j = 1
!       dy_b(i,j) = (d_cy( j+1,i+1) + d_cy( j+1,i  ) - d_cy( j  ,i+1) - d_cy( j  ,i  )) / (2.0_wp * dx)
!       j = grid%ny-1
!       dy_b(i,j) = (d_cy( j  ,i+1) + d_cy( j  ,i  ) - d_cy( j-1,i+1) - d_cy( j-1,i  )) / (2.0_wp * dx)
!     end do
    
    
!   end subroutine ddy_cy_to_b_2D
!   subroutine ddx_cy_to_b_2D( d_cy, dx_b)
!     ! Input:  scalar on the Acy grid
!     ! Output: its x-derivative on the Ab grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION(       grid%ny-1, grid%nx  ), intent(IN)    :: d_cy
!     real(wp), DIMENSION(       grid%ny-1, grid%nx-1), intent(OUT)   :: dx_b
    
!     ! Local variables:
!     integer :: i, j, nx, ny
    
!     ! Interior
!     do i = grid%i1, MIN(grid%nx-1,grid%i2)
!     do j = 1, ny-1
!       dx_b(i,j) = (d_cy(i+1,j) - d_cy(i,j)) / dx
!     end do
!     end do
    
    
!   end subroutine ddx_cy_to_b_2D
!   subroutine ddy_cx_to_b_2D( d_cx, dy_b)
!     ! Input:  scalar on the Acx grid
!     ! Output: its y-derivative on the Ab grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION(       grid%ny  , grid%nx-1), intent(IN)    :: d_cx
!     real(wp), DIMENSION(       grid%ny-1, grid%nx-1), intent(OUT)   :: dy_b
    
!     ! Local variables:
!     integer :: i, j, nx, ny
    
!     ! Interior
!     do i = grid%i1, MIN(grid%nx-1,grid%i2)
!     do j = 1, ny-1
!       dy_b(i,j) = (d_cx(i,j+1) - d_cx(i,j)) / dx
!     end do
!     end do
    
    
!   end subroutine ddy_cx_to_b_2D
  
! ! =============================================
! ! ===== Mapping between (staggered) grids =====
! ! =============================================

!   ! Aa to Acx/Acy
  
!   ! 2D
!   subroutine map_a_to_cx_2D( d_a, d_cx)
!     ! Input:  scalar on the Aa grid
!     ! Output: the same on the Acx grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), intent(IN)    :: d_a
!     real(wp), DIMENSION(       grid%ny  , grid%nx-1), intent(OUT)   :: d_cx
    
!     ! Local variables:
!     integer :: i, j, nx, ny
    
!     do i = grid%i1, MIN(grid%nx-1,grid%i2)
!     do j = 1, ny
!       d_cx(i,j) = (d_a(i,j) + d_a(i+1,j)) / 2.0_wp
!     end do
!     end do
    
    
!   end subroutine map_a_to_cx_2D
!   subroutine map_a_to_cy_2D( d_a, d_cy)
!     ! Input:  scalar on the Aa grid
!     ! Output: the same on the Acy grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), intent(IN)    :: d_a
!     real(wp), DIMENSION(       grid%ny-1, grid%nx  ), intent(OUT)   :: d_cy
    
!     ! Local variables:
!     integer :: i, j, nx, ny
    
!     do i = 1, nx
!     do j = 1, ny-1
!       d_cy(i,j) = (d_a(i,j) + d_a(i,j+1)) / 2.0_wp
!     end do
!     end do
    
    
!   end subroutine map_a_to_cy_2D
!   ! 3D
!   subroutine map_a_to_cx_3D( d_a, d_cx)
!     ! Input:  scalar on the Aa grid
!     ! Output: the same on the Acx grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), intent(IN)    :: d_a
!     real(wp), DIMENSION( C%nZ, grid%ny  , grid%nx-1), intent(OUT)   :: d_cx
    
!     ! Local variables:
!     integer :: i, j, nx, ny,k
    
!     do i = grid%i1, MIN(grid%nx-1,grid%i2)
!     do j = 1, ny
!     do k = 1, C%nZ
!       d_cx( k,j,i) = (d_a( k,j,i) + d_a( k,j,i+1)) / 2.0_wp
!     end do
!     end do
!     end do
    
    
!   end subroutine map_a_to_cx_3D
!   subroutine map_a_to_cy_3D( d_a, d_cy)
!     ! Input:  scalar on the Aa grid
!     ! Output: the same on the Acy grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), intent(IN)    :: d_a
!     real(wp), DIMENSION( C%nZ, grid%ny-1, grid%nx  ), intent(OUT)   :: d_cy
    
!     ! Local variables:
!     integer :: i, j, nx, ny,k
    
!     do i = 1, nx
!     do j = 1, ny-1
!     do k = 1, C%nZ
!       d_cy( k,j,i) = (d_a( k,j,i) + d_a( k,j+1,i)) / 2.0_wp
!     end do
!     end do
!     end do
    
    
!   end subroutine map_a_to_cy_3D
  
!   ! Acx/Acy to Aa
  
!   ! 2D
!   subroutine map_cx_to_a_2D( d_cx, d_a)
!     ! Input:  scalar on the Acx grid
!     ! Output: the same on the Aa grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION(       grid%ny  , grid%nx-1), intent(IN)    :: d_cx
!     real(wp), intent(OUT)   :: d_a
    
!     ! Local variables:
!     integer :: i, j, nx, ny
    
!     do i = 2, nx-1
!     do j = 1, ny
!       d_a(i,j) = (d_cx(i-1,j) + d_cx(i,j)) / 2.0_wp
!     end do
!     end do
    
    
!     d_a(1,:) = d_cx( grid%j1:grid%j2,1        )
!     d_a(nx,:) = d_cx( grid%j1:grid%j2,grid%nx-1)
    
    
!   end subroutine map_cx_to_a_2D
!   subroutine map_cy_to_a_2D( d_cy, d_a)
!     ! Input:  scalar on the Acy grid
!     ! Output: the same on the Aa grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION(       grid%ny-1, grid%nx  ), intent(IN)    :: d_cy
!     real(wp), intent(OUT)   :: d_a
    
!     ! Local variables:
!     integer :: i, j, nx, ny
    
!     do i = 1, nx
!     do j = 2, ny-1
!       d_a(i,j) = (d_cy(i,j-1) + d_cy(i,j)) / 2.0_wp
!     end do
!     end do
    
    
!     d_a(:,1) = d_cy( 1        ,grid%i1:grid%i2)
!     d_a(:,ny) = d_cy( grid%ny-1,grid%i1:grid%i2)
    
    
!   end subroutine map_cy_to_a_2D
!   ! 3D
!   subroutine map_cx_to_a_3D( d_cx, d_a)
!     ! Input:  scalar on the Acx grid
!     ! Output: the same on the Aa grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION( C%nZ, grid%ny  , grid%nx-1), intent(IN)    :: d_cx
!     real(wp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), intent(OUT)   :: d_a
    
!     ! Local variables:
!     integer :: i, j, nx, ny,k
    
!     do i = 2, nx-1
!     do j = 1, ny
!     do k = 1, C%nZ
!       d_a( k,j,i) = (d_cx( k,j,i-1) + d_cx( k,j,i)) / 2.0_wp
!     end do
!     end do
!     end do
    
    
!     d_a( :,grid%j1:grid%j2,1      ) = d_cx( :,grid%j1:grid%j2,1        )
!     d_a( :,grid%j1:grid%j2,grid%nx) = d_cx( :,grid%j1:grid%j2,grid%nx-1)
    
    
!   end subroutine map_cx_to_a_3D
!   subroutine map_cy_to_a_3D( d_cy, d_a)
!     ! Input:  scalar on the Acy grid
!     ! Output: the same on the Aa grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION( C%nZ, grid%ny-1, grid%nx  ), intent(IN)    :: d_cy
!     real(wp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), intent(OUT)   :: d_a
    
!     ! Local variables:
!     integer :: i, j, nx, ny,k
    
!     do i = 1, nx
!     do j = 2, ny-1
!     do k = 1, C%nZ
!       d_a( k,j,i) = (d_cy( k,j-1,i) + d_cy( k,j,i)) / 2.0_wp
!     end do
!     end do
!     end do
    
    
!     d_a( :,1      ,grid%i1:grid%i2) = d_cy( :,1        ,grid%i1:grid%i2)
!     d_a( :,grid%ny,grid%i1:grid%i2) = d_cy( :,grid%ny-1,grid%i1:grid%i2)
    
    
!   end subroutine map_cy_to_a_3D
  
!   ! Acx/Acy to Acy/Acx
!   subroutine map_cx_to_cy_2D( d_cx, d_cy)
!     ! Input:  scalar on the Acx grid
!     ! Output: the same on the Acy grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION(       grid%ny  , grid%nx-1), intent(IN)    :: d_cx
!     real(wp), DIMENSION(       grid%ny-1, grid%nx  ), intent(OUT)   :: d_cy
    
!     ! Local variables:
!     integer :: i, j, nx, ny
    
!     ! Interior
!     do i = 2, nx-1
!     do j = 1, ny-1
!       d_cy(i,j) = (d_cx( j  ,i-1) + d_cx( j  ,i  ) + d_cx( j+1,i-1) + d_cx( j+1,i  )) / 4.0_wp
!     end do
!     end do
    
    
!     ! Boundaries
!     do j = grid%j1, MIN(grid%ny-1,grid%j2)
!       i = 1
!       d_cy(i,j) = (d_cx( j  ,i  ) + d_cx( j+1,i  )) / 2.0_wp
!       i = grid%nx
!       d_cy(i,j) = (d_cx( j  ,i-1) + d_cx( j+1,i-1)) / 2.0_wp
!     end do
    
    
!   end subroutine map_cx_to_cy_2D
!   subroutine map_cy_to_cx_2D( d_cy, d_cx)
!     ! Input:  scalar on the Acy grid
!     ! Output: the same on the Acx grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION(       grid%ny-1, grid%nx  ), intent(IN)    :: d_cy
!     real(wp), DIMENSION(       grid%ny  , grid%nx-1), intent(OUT)   :: d_cx
    
!     ! Local variables:
!     integer :: i, j, nx, ny
    
!     ! Interior
!     do i = grid%i1, MIN(grid%nx-1,grid%i2)
!     do j = 2, ny-1
!       d_cx(i,j) = (d_cy( j-1,i  ) + d_cy( j-1,i+1) + d_cy( j  ,i  ) + d_cy( j  ,i+1)) / 4.0_wp
!     end do
!     end do
    
    
!     ! Boundaries
!     do i = grid%i1, MIN(grid%nx-1,grid%i2)
!       j = 1
!       d_cx(i,j) = (d_cy( j  ,i  ) + d_cy( j  ,i+1)) / 2.0_wp
!       j = grid%ny
!       d_cx(i,j) = (d_cy( j-1,i  ) + d_cy( j-1,i+1)) / 2.0_wp
!     end do
    
    
!   end subroutine map_cy_to_cx_2D
  
!   ! Aa to Ab
!   subroutine map_a_to_b_2D( d_a, d_b)
!     ! Input:  scalar on the Aa grid
!     ! Output: the same on the Ab grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), intent(IN)    :: d_a
!     real(wp), DIMENSION(       grid%ny-1, grid%nx-1), intent(OUT)   :: d_b
    
!     ! Local variables:
!     integer :: i, j, nx, ny
    
!     do i = grid%i1, MIN(grid%nx-1,grid%i2)
!     do j = 1, ny-1
!       d_b(i,j) = (d_a(i,j) + d_a(i+1,j) + d_a(i,j+1) + d_a( j+1,i+1)) / 4.0_wp
!     end do
!     end do
    
    
!   end subroutine map_a_to_b_2D
  
!   ! Acx/Acy to Ab
!   subroutine map_cx_to_b_2D( d_cx, d_b)
!     ! Input:  scalar on the Acx grid
!     ! Output: the same on the Ab grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION(       grid%ny  , grid%nx-1), intent(IN)    :: d_cx
!     real(wp), DIMENSION(       grid%ny-1, grid%nx-1), intent(OUT)   :: d_b
    
!     ! Local variables:
!     integer :: i, j, nx, ny
    
!     do i = grid%i1, MIN(grid%nx-1,grid%i2)
!     do j = 1, ny-1
!       d_b(i,j) = (d_cx(i,j) + d_cx(i,j+1)) / 2.0_wp
!     end do
!     end do
    
    
!   end subroutine map_cx_to_b_2D
!   subroutine map_cy_to_b_2D( d_cy, d_b)
!     ! Input:  scalar on the Acy grid
!     ! Output: the same on the Ab grid
    
!     implicit none
    
!     ! In/output variables:
    
!     real(wp), DIMENSION(       grid%ny-1, grid%nx  ), intent(IN)    :: d_cy
!     real(wp), DIMENSION(       grid%ny-1, grid%nx-1), intent(OUT)   :: d_b
    
!     ! Local variables:
!     integer :: i, j, nx, ny
    
!     do i = grid%i1, MIN(grid%nx-1,grid%i2)
!     do j = 1, ny-1
!       d_b(i,j) = (d_cy(i,j) + d_cy(i+1,j)) / 2.0_wp
!     end do
!     end do
    
    
!   end subroutine map_cy_to_b_2D
  
end module grid_calcs



