
module interp2D
    
    implicit none 

    ! === coord_constants ===================================

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the precision of the coord library (sp,dp)
    integer,  parameter :: wp = sp 


    ! Missing value and aliases
    real(dp), parameter :: MISSING_VALUE_DEFAULT = -9999.0_dp 
    real(dp), parameter :: mv = MISSING_VALUE_DEFAULT
    
    ! Error distance (very large) and error index 
    real(dp), parameter :: ERR_DIST = 1E8_dp 
    integer,  parameter :: ERR_IND  = -1 

    ! Mathematical constants
    real(dp), parameter  :: pi  = 2._dp*acos(0._dp)
    real(dp), parameter  :: degrees_to_radians = pi / 180._dp  ! Conversion factor between radians and degrees
    real(dp), parameter  :: radians_to_degrees = 180._dp / pi  ! Conversion factor between degrees and radians
     
    ! ======================================================

    interface interp_bilinear 
        module procedure interp_bilinear_dble 
        module procedure interp_bilinear_points_dble
        module procedure interp_bilinear_float 
        module procedure interp_bilinear_points_float
    end interface

    interface interp_nearest 
        module procedure interp_nearest_dble, interp_nearest_int
    end interface

    interface interp_nearest_fast
        module procedure interp_nearest_fast_dble
        module procedure interp_nearest_fast_float
    end interface

    interface fill_weighted
        module procedure fill_weighted_dble 
    end interface

    interface fill_mean
        module procedure fill_mean_dble 
    end interface

    interface fill_nearest
        module procedure fill_nearest_dble, fill_nearest_int  
        module procedure fill_nearest_dble_new
        module procedure fill_nearest_float
    end interface

    interface fill_bilinear
        module procedure fill_bilinear_dble !, fill_bilinear_float
    end interface

    interface fill_poisson
        module procedure fill_poisson_float 
        module procedure fill_poisson_dble
    end interface 
        
    interface filter_poisson 
        module procedure filter_poisson_float 
        module procedure filter_poisson_dble
    end interface 

    private
    public :: interp_bilinear, interp_nearest, interp_nearest_fast
    public :: fill_weighted, fill_nearest, fill_mean
    public :: fill_poisson 
    public :: filter_poisson
    public :: diffuse, limit_gradient 

contains

    function interp_bilinear_float(x,y,z,xout,yout,missing_value,mask,fill) result(zout)
        ! Find closest x-indices and closest y-indices on original
        ! grid (assume these and next indices will bracket our point)
        ! Perform weighted interpolation 

        implicit none 

        real(sp), dimension(:) :: x, y, xout, yout 
        real(sp), dimension(:,:) :: z
        real(dp), optional :: missing_value 
        logical, dimension(:,:), optional :: mask 
        logical, optional :: fill 
        real(sp), dimension(size(xout,1),size(yout,1)) :: zout 
        
        zout = real( interp_bilinear_dble(dble(x),dble(y),dble(z),dble(xout),dble(yout), &
                                          missing_value,mask,fill), sp)

        return 

    end function interp_bilinear_float

    function interp_bilinear_points_float(is_points,x,y,z,xout,yout,missing_value,mask) result(zout)
        ! Find closest x-indices and closest y-indices on original
        ! grid (assume these and next indices will bracket our point)
        ! Perform bilinear interpolation to new points

        implicit none 

        logical :: is_points
        real(sp), dimension(:) :: x, y, xout, yout 
        real(sp), dimension(:,:) :: z
        real(dp), optional :: missing_value 
        logical, dimension(:), optional :: mask 
        real(sp), dimension(size(xout,1)) :: zout 
        
        zout = real( interp_bilinear_points_dble(is_points,dble(x),dble(y),dble(z), &
                                                 dble(xout),dble(yout),missing_value,mask), sp)

        return 

    end function interp_bilinear_points_float

    function interp_bilinear_dble(x,y,z,xout,yout,missing_value,mask,fill) result(zout)
        ! Find closest x-indices and closest y-indices on original
        ! grid (assume these and next indices will bracket our point)
        ! Perform weighted interpolation 

        implicit none 

        real(dp), dimension(:) :: x, y, xout, yout 
        real(dp), dimension(:,:) :: z
        real(dp), optional :: missing_value 
        logical, dimension(:,:), optional :: mask 
        logical, optional :: fill 
        real(dp), dimension(size(xout,1),size(yout,1)) :: zout 
        logical,  dimension(size(xout,1),size(yout,1)) :: mask_interp 
        real(dp) :: missing_val 
        logical  :: fill_missing 
        integer, dimension(size(xout,1)) :: x_idx
        integer, dimension(size(yout,1)) :: y_idx
        
        integer :: i, j, i1, j1  
        integer :: nx, ny, nx1, ny1 
        double precision :: alpha1, alpha2, p0, p1 

        nx = size(x,1)
        ny = size(y,1)

        nx1 = size(xout,1)
        ny1 = size(yout,1)

        ! Determine which points we are interested in interpolating
        mask_interp = .TRUE. 
        if (present(mask)) mask_interp = mask 

        ! Determine missing value if present 
        missing_val = MISSING_VALUE_DEFAULT
        if (present(missing_value)) missing_val = missing_value

        fill_missing = .FALSE. 
        if (present(fill)) fill_missing = fill 

        ! Get x-indices corresponding to nearest neighbor
        ! greater-than-equal-to x-value of interest
        do i1 = 1, nx1 

            if (xout(i1) .le. x(1)) then 
                x_idx(i1) = -1

            else if (xout(i1) .ge. (x(nx))) then 
                x_idx(i1) = -2 
            else 

                do i = 1, nx 
                    if (x(i) .ge. xout(i1)) exit 
                end do 

                x_idx(i1) = i
            end if 

        end do 

        ! Get y-indices corresponding to nearest neighbor
        ! greater-than-equal-to y-value of interest
        do j1 = 1, ny1 

            if (yout(j1) .le. y(1)) then 
                y_idx(j1) = -1

            else if (yout(j1) .ge. (y(ny))) then 
                y_idx(j1) = -2 
            else 

                do j = 1, ny 
                    if (y(j) .ge. yout(j1)) exit 
                end do 

                y_idx(j1) = j 
            end if 

        end do 

        ! Now loop over output grid points and perform
        ! bilinear interpolation where desired 
        zout = missing_val 
        do i1 = 1, nx1 
        do j1 = 1, ny1 

            ! Only interpolate points of interest 
            if (mask_interp(i1,j1)) then 

                i = x_idx(i1)
                j = y_idx(j1) 

                ! Only interpolate points inside the original grid (for now)
                if (i .gt. 0 .and. i-1 .gt. 0 .and. j .gt. 0 .and. j-1 .gt. 0) then 
                    
                    ! Only interpolate points with all neighbors available
                    if (count([z(i-1,j),z(i,j),z(i,j-1),z(i-1,j-1)] .eq. missing_val) .eq. 0) then
                        alpha1 = (xout(i1) - x(i-1)) / (x(i)-x(i-1))
                        p0 = z(i-1,j-1) + alpha1*(z(i,j-1)-z(i-1,j-1))
                        p1 = z(i-1,j)   + alpha1*(z(i,j)-z(i-1,j))
                        
                        alpha2 = (yout(j1) - y(j-1)) / (y(j)-y(j-1))
                        zout(i1,j1) = p0 + alpha2*(p1-p0)
                    end if 

                end if 

            end if 

        end do 
        end do
        

!         if (fill_missing) then
!             write(*,*) "Filling in...", missing_val, count(zout .eq. missing_val), nx1*ny1
!             call fill_weighted(zout,missing_val)
!             write(*,*) "Filled in... ", missing_val, count(zout .eq. missing_val), nx1*ny1
!         end if 

        return 

    end function interp_bilinear_dble

    function interp_bilinear_points_dble(is_points,x,y,z,xout,yout,missing_value,mask) result(zout)
        ! Find closest x-indices and closest y-indices on original
        ! grid (assume these and next indices will bracket our point)
        ! Perform bilinear interpolation to new points

        implicit none 

        logical :: is_points
        real(dp), dimension(:) :: x, y, xout, yout 
        real(dp), dimension(:,:) :: z
        real(dp), optional :: missing_value 
        logical, dimension(:), optional :: mask 
        real(dp), dimension(size(xout,1)) :: zout 
        logical,  dimension(size(xout,1)) :: mask_interp 
        real(dp) :: missing_val 
        logical  :: fill_missing 
        integer, dimension(size(xout,1)) :: x_idx
        integer, dimension(size(yout,1)) :: y_idx
        
        integer :: i, j, i1, j1  
        integer :: nx, ny, nx1, ny1 
        double precision :: alpha1, alpha2, p0, p1 

        nx = size(x,1)
        ny = size(y,1)

        nx1 = size(xout,1)
        ny1 = size(yout,1)

        ! Determine which points we are interested in interpolating
        mask_interp = .TRUE. 
        if (present(mask)) mask_interp = mask 

        ! Determine missing value if present 
        missing_val = MISSING_VALUE_DEFAULT
        if (present(missing_value)) missing_val = missing_value

        ! Get x-indices corresponding to nearest neighbor
        ! greater-than-equal-to x-value of interest
        do i1 = 1, nx1 

            if (xout(i1) .le. x(1)) then 
                x_idx(i1) = -1

            else if (xout(i1) .ge. (x(nx))) then 
                x_idx(i1) = -2 
            else 

                do i = 1, nx 
                    if (x(i) .ge. xout(i1)) exit 
                end do 

                x_idx(i1) = i
            end if 

        end do 

        ! Get y-indices corresponding to nearest neighbor
        ! greater-than-equal-to y-value of interest
        do j1 = 1, ny1 

            if (yout(j1) .le. y(1)) then 
                y_idx(j1) = -1

            else if (yout(j1) .ge. (y(ny))) then 
                y_idx(j1) = -2 
            else 

                do j = 1, ny 
                    if (y(j) .ge. yout(j1)) exit 
                end do 

                y_idx(j1) = j 
            end if 

        end do 

        ! Now loop over output grid points and perform
        ! bilinear interpolation where desired 
        zout = missing_val 
        do i1 = 1, nx1 

            ! Only interpolate points of interest 
            if (mask_interp(i1)) then 

                i = x_idx(i1)
                j = y_idx(i1) 

                ! Only interpolate points inside the original grid (for now)
                if (i .gt. 0 .and. i-1 .gt. 0 .and. j .gt. 0 .and. j-1 .gt. 0) then 
                    
                    ! Only interpolate points with all neighbors available
                    if (count([z(i-1,j),z(i,j),z(i,j-1),z(i-1,j-1)] .eq. missing_val) .eq. 0) then
                        alpha1 = (xout(i1) - x(i-1)) / (x(i)-x(i-1))
                        p0 = z(i-1,j-1) + alpha1*(z(i,j-1)-z(i-1,j-1))
                        p1 = z(i-1,j)   + alpha1*(z(i,j)-z(i-1,j))
                        
                        alpha2 = (yout(i1) - y(j-1)) / (y(j)-y(j-1))
                        zout(i1) = p0 + alpha2*(p1-p0)
                    end if 

                end if 

            end if 

        end do
        
        return 

    end function interp_bilinear_points_dble

    function interp_nearest_dble0(x,y,z,xout,yout,missing_value,mask) result(zout)
        ! Find closest x-indices and closest y-indices on original
        ! grid (assume these and next indices will bracket our point)
        ! Pick closest point 
        ! Don't use - this is very slow!!! Could probably be improved with vectors
        implicit none 

        real(dp), dimension(:) :: x, y, xout, yout 
        real(dp), dimension(:,:) :: z
        real(dp), optional :: missing_value 
        logical, dimension(:,:), optional :: mask 
        real(dp), dimension(size(xout,1),size(yout,1)) :: zout 
        logical,  dimension(size(xout,1),size(yout,1)) :: mask_interp 
        real(dp) :: missing_val 

        integer, dimension(size(z,1),size(z,2)) :: dist
        real(dp) :: mindist 

        integer :: i, j, i1, j1, ij(2)
        integer :: nx, ny, nx1, ny1   
        integer :: imin, jmin 

        nx = size(x,1)
        ny = size(y,1)

        nx1 = size(xout,1)
        ny1 = size(yout,1)

        ! Determine which points we are interested in interpolating
        mask_interp = .TRUE. 
        if (present(mask)) mask_interp = mask 

        ! Determine missing value if present 
        missing_val = MISSING_VALUE_DEFAULT
        if (present(missing_value)) missing_val = missing_value

        ! Now loop over output grid points and perform
        ! nearest neighbor interpolation where desired 
        zout = missing_val 

        do i1 = 1, nx1 
        do j1 = 1, ny1 

            ! Only interpolate points of interest 
            if (mask_interp(i1,j1)) then 

                dist = ERR_DIST 
                imin = minloc(dabs(x-xout(i1)),dim=1)
                if (imin .le. 1)  imin = 2
                if (imin .ge. nx) imin = nx-1

                jmin = minloc(dabs(y-yout(j1)),dim=1)
                if (jmin .le. 1)  jmin = 2
                if (jmin .ge. ny) jmin = ny-1

!                 do i = 1, nx 
!                 do j = 1, ny
                do i = imin-1,imin+1 
                do j = jmin-1,jmin+1
                        dist(i,j) = dsqrt( (x(i)-xout(i1))**2 + (y(j)-yout(j1))**2 )
                end do 
                end do 

                where(z .eq. missing_val) dist = ERR_DIST

                ij = minloc(dist) 
                i  = ij(1)
                j  = ij(2)

                zout(i1,j1) = z(i,j) 
            end if 

        end do 
        end do
        

        return 

    end function interp_nearest_dble0

    function interp_nearest_dble(x,y,z,xout,yout,missing_value,mask) result(zout)
        ! Find closest x-indices and closest y-indices on original
        ! grid (assume these and next indices will bracket our point)
        ! Pick closest point 
        ! Faster than other version that uses a bigger matrix
        implicit none 

        real(dp), dimension(:) :: x, y, xout, yout 
        real(dp), dimension(:,:) :: z
        real(dp), optional :: missing_value 
        logical, dimension(:,:), optional :: mask 
        real(dp), dimension(size(xout,1),size(yout,1)) :: zout 
        logical,  dimension(size(xout,1),size(yout,1)) :: mask_interp 
        real(dp) :: missing_val 

        integer, parameter :: nr = 3, nd = 2*nr+1 
        integer, dimension(nd,nd) :: dist
        real(dp) :: mindist 

        integer :: i, j, i1, j1, ij(2), imin, jmin  
        integer :: nx, ny, nx1, ny1   

        nx = size(x,1)
        ny = size(y,1)

        nx1 = size(xout,1)
        ny1 = size(yout,1)

        ! Determine which points we are interested in interpolating
        mask_interp = .TRUE. 
        if (present(mask)) mask_interp = mask 

        ! Determine missing value if present 
        missing_val = MISSING_VALUE_DEFAULT
        if (present(missing_value)) missing_val = missing_value

        ! Now loop over output grid points and perform
        ! nearest neighbor interpolation where desired 
        zout = missing_val 

        do i1 = 1, nx1 
        do j1 = 1, ny1 

            ! Only interpolate points of interest 
            if (mask_interp(i1,j1)) then 

                imin = minloc(dabs(x-xout(i1)),dim=1)
                if (imin .le. nr)      imin = nr+1
                if (imin .ge. nx-nr)   imin = nx-nr-1

                jmin = minloc(dabs(y-yout(j1)),dim=1)
                if (jmin .le. nr)      jmin = nr+1
                if (jmin .ge. ny-nr)   jmin = ny-nr-1  

                do i = 1, nd
                    do j = 1, nd 
                        dist(i,j) = dsqrt( (x(imin+(i-nr))-xout(i1))**2 + (y(jmin+(j-nr))-yout(j1))**2 )
                    end do 
                end do 

                where(z(imin-nr:imin+nr,jmin-nr:jmin+nr) .eq. missing_val) dist = ERR_DIST

                ij = minloc(dist) 
                i  = imin+(ij(1)-nr)
                j  = jmin+(ij(2)-nr) 

                if (z(i,j) .ne. missing_val) zout(i1,j1) = z(i,j) 
            end if 

        end do 
        end do
        

        return 

    end function interp_nearest_dble

    function interp_nearest_fast_dble(x,y,z,xout,yout,max_dist_fac,missing_value,mask) result(zout)
        ! Find closest x-indices and closest y-indices on original
        ! grid (assume these and next indices will bracket our point)
        ! Pick closest point 
        ! Faster than other version that uses a bigger matrix
        implicit none 

        real(dp), dimension(:) :: x, y, xout, yout 
        real(dp), dimension(:,:) :: z
        real(dp) :: max_dist_fac
        real(dp), optional :: missing_value
        logical, dimension(:,:), optional :: mask 
        real(dp), dimension(size(xout,1),size(yout,1)) :: zout 
        logical,  dimension(size(xout,1),size(yout,1)) :: mask_interp 
        real(dp) :: missing_val, max_dist  

        integer :: i, j, i1, j1, ij(2), imin, jmin  
        integer :: nx, ny, nx1, ny1   

        nx = size(x,1)
        ny = size(y,1)

        nx1 = size(xout,1)
        ny1 = size(yout,1)

        ! Determine which points we are interested in interpolating
        mask_interp = .TRUE. 
        if (present(mask)) mask_interp = mask 

        ! Determine missing value if present 
        missing_val = MISSING_VALUE_DEFAULT
        if (present(missing_value)) missing_val = missing_value

        ! Maximum distance to be considered a nearest neighbor 
        max_dist = dsqrt( (x(2)-x(1))**2 + (y(2)-y(1))**2 ) * max_dist_fac

        ! Now loop over output grid points and perform
        ! nearest neighbor interpolation where desired 
        zout = missing_val 

        do i1 = 1, nx1 
        do j1 = 1, ny1 

            ! Only interpolate points of interest 
            if (mask_interp(i1,j1)) then 

                imin = minloc(dabs(x-xout(i1)),dim=1)
                jmin = minloc(dabs(y-yout(j1)),dim=1)

                if ( dsqrt((x(imin)-xout(i1))**2+(y(jmin)-yout(j1))**2) .le. max_dist ) then 
                    zout(i1,j1) = z(imin,jmin)
                end if 

            end if 

        end do 
        end do
        

        return 

    end function interp_nearest_fast_dble

    function interp_nearest_fast_float(x,y,z,xout,yout,max_dist_fac,missing_value,mask) result(zout)
        ! Find closest x-indices and closest y-indices on original
        ! grid (assume these and next indices will bracket our point)
        ! Pick closest point 
        ! Faster than other version that uses a bigger matrix
        implicit none 

        real(sp), dimension(:) :: x, y, xout, yout 
        real(sp), dimension(:,:) :: z
        real(sp) :: max_dist_fac
        real(sp), optional :: missing_value
        logical, dimension(:,:), optional :: mask 
        real(sp), dimension(size(xout,1),size(yout,1)) :: zout 
        logical,  dimension(size(xout,1),size(yout,1)) :: mask_interp 
        real(sp) :: missing_val, max_dist  

        integer :: i, j, i1, j1, ij(2), imin, jmin  
        integer :: nx, ny, nx1, ny1   

        nx = size(x,1)
        ny = size(y,1)

        nx1 = size(xout,1)
        ny1 = size(yout,1)

        ! Determine which points we are interested in interpolating
        mask_interp = .TRUE. 
        if (present(mask)) mask_interp = mask 

        ! Determine missing value if present 
        missing_val = MISSING_VALUE_DEFAULT
        if (present(missing_value)) missing_val = missing_value

        ! Maximum distance to be considered a nearest neighbor 
        max_dist = sqrt( (x(2)-x(1))**2 + (y(2)-y(1))**2 ) * max_dist_fac

        ! Now loop over output grid points and perform
        ! nearest neighbor interpolation where desired 
        zout = missing_val 

        do i1 = 1, nx1 
        do j1 = 1, ny1 

            ! Only interpolate points of interest 
            if (mask_interp(i1,j1)) then 

                imin = minloc(abs(x-xout(i1)),dim=1)
                jmin = minloc(abs(y-yout(j1)),dim=1)

                if ( sqrt((x(imin)-xout(i1))**2+(y(jmin)-yout(j1))**2) .le. max_dist ) then 
                    zout(i1,j1) = z(imin,jmin)
                end if 

            end if 

        end do 
        end do
        

        return 

    end function interp_nearest_fast_float

    subroutine fill_nearest_dble_new(var,missing_value,fill_value,fill_dist,n,dx,mask)

        implicit none 

        real(8), intent(INOUT) :: var(:,:)              ! Array with missing values to be populated completely
        real(8), intent(IN)    :: missing_value         ! Specified missing value indicating points to be filled in
        real(8), intent(IN)    :: fill_value            ! Value to impose when no neighbors are found within fill_dist
        real(8), intent(IN)    :: fill_dist             ! [km] Distance to nearest neighbor at which value=fill_value
        integer, intent(IN)    :: n                     ! Average of n neighbors 
        real(8), intent(IN)    :: dx                    ! [m] 
        logical, intent(IN), optional :: mask(:,:)      ! Mask for where to apply filling 

        ! Local variables 
        integer :: i, j, nx, ny, i1, j1, q, n_now, ij(2) 
        integer :: ntot 
        real(8) :: dx_km 
        real(8) :: dist_now, f_d 

        real(8), allocatable :: var0(:,:) 
        real(8), allocatable :: dist(:,:) 
        logical, allocatable :: mask_apply(:,:) 

        nx = size(var,1)
        ny = size(var,2) 

        allocate(var0(nx,ny)) 
        allocate(dist(nx,ny)) 

        ! Define resolution in km 
        dx_km = dx*1e-3 

        ! Store initial field 
        var0 = var 

        ntot = 0 

        ! Loop over missing values, look for nearest non-missing neighbor
        do j = 1, ny 
        do i = 1, nx 

            if (var(i,j) .eq. missing_value) then 
                ! Find a neighbor value in var0 

                ! Populate distance matrix where necessary 
                dist = MV 
                do j1 = 1, ny 
                do i1 = 1, nx 
                    if (var0(i1,j1) .ne. MV) then 
                        dist(i1,j1) = sqrt( real( (i1-i)**2 + (j1-j)**2 ) ) * dx_km 
                    end if 
                end do 
                end do 

                n_now    = 0 
                var(i,j) = 0.0
                dist_now = 0.0 

                do q = 1, n 
                    ! Loop over nearest neighbors to get average 

                    ! Find minimum populated neighbor 
                    ij = minloc(dist,mask=dist.ne.MV .and. var0.ne.MV)

                    ! Check if no neighbors found 
                    if (ij(1) .eq. 0) exit 

                    ! Populate with neighbor value 
                    var(i,j) = var(i,j) + var0(ij(1),ij(2))
                    dist_now = dist_now + dist(ij(1),ij(2))
                    n_now = n_now + 1 

                    ! Reset distance of neighbor to zero so it cannot be used again
                    dist(ij(1),ij(2)) = MV 
                end do 

                ! If no neighbors found, use fill value 
                if (n_now .eq. 0) var(i,j) = fill_value 

                ! Take average if multiple points used 
                if (n_now .gt. 1) then 
                    
                    ! Determine mean distance to neighbors and weighting function versus distance
                    dist_now = dist_now / real(n_now,wp) 
                    
                    if (fill_value .ne. missing_value) then 
                        ! Apply weighted average of mean neighbor value and fill value

                        f_d      = 1.0 - min( dist_now/fill_dist, 1.0 )
                        var(i,j) = f_d * (var(i,j) / real(n_now,wp)) + (1.0-f_d)*fill_value

                    end if 

                end if 

                ! Add this missing point to total for diagnostics 
                ntot = ntot + 1 
            end if 

        end do
        end do 

        return 

    end subroutine fill_nearest_dble_new

    subroutine fill_nearest_dble(z,missing_value,fill_value,mask,n)
        implicit none 
        double precision :: z(:,:) 
        double precision :: missing_value 
        double precision, optional :: fill_value
        logical, optional :: mask(:,:) 

        integer, optional :: n 

        ! Local variables 
        integer :: nr 

        integer :: q, nx, ny, i, j
        integer, parameter :: qmax = 100 ! Iterations 
        double precision, dimension (:,:), allocatable :: neighb, weight, weight0 
        double precision :: wtot, mval 
        double precision, dimension(:,:), allocatable :: filled
        integer :: ij(2)
        logical, allocatable :: mask_apply(:,:) 

        nx = size(z,1)
        ny = size(z,2) 
        
        allocate(mask_apply(nx,ny)) 

        mask_apply = .TRUE. 
        if (present(mask)) mask_apply = mask 

        nr = 4
        if (present(n)) nr = n 

        nx = size(z,1)
        ny = size(z,2) 

        allocate(filled(nx,ny))
        allocate(neighb(2*nr+1,2*nr+1),weight(2*nr+1,2*nr+1),weight0(2*nr+1,2*nr+1))

        if (present(fill_value)) then
            where(z .eq. missing_value .and. mask_apply) z = fill_value 
        end if 

        ! Define the default neighbor weighting 
        ! All weights == 1 incorporates some numerical diffusion to make
        ! the resulting extrapolated values smoother
        weight0 = 1
        do i = 1, 2*nr+1
        do j = 1, 2*nr+1 
            weight0(i,j) = 1/max(1d-1,dsqrt(dble((nr+1)-i)**2 + dble((nr+1)-j)**2))
        end do 
        end do  

        do q = 1, qmax 

            filled = missing_value 

            do i = 1+nr, nx-nr
                do j = 1+nr, ny-nr
                    
                    if (z(i,j) .eq. missing_value .and. mask_apply(i,j)) then 

                        neighb = z(i-nr:i+nr,j-nr:j+nr)

                        weight = weight0
                        where (neighb .eq. missing_value)  weight = 0.d0
                        wtot = sum(weight)

                        ! Should only use the nearest neighbor in this case
                        if (wtot .gt. 0) then 
                            ij = maxloc(weight) 
                            filled(i,j) = neighb(ij(1),ij(2))
                        end if 

                    end if 
                end do 
            end do 

            where(filled .ne. missing_value) z = filled 
            if ( count(z(1+nr+1:nx-(nr+1),1+nr+1:ny-(nr+1)) .eq. missing_value) .eq. 0 ) exit 

!             write(*,*) "Still missing... ", count(z(1+nr:nx-nr,1+nr:ny:nr) .eq. missing_value), &
!                         " of ", nx*ny
        end do 

        if (q .ge. qmax) then 
            write(*,*) "Too many iterations to fill array of size: ", nx, ny 
            write(*,*) "Remaining missing values: ", count(z .eq. missing_value)
!             stop 
        end if 

        ! Fill in boundaries too 
        do i = 1+nr, 1, -1
            where(z(i,:) .eq. missing_value .and. mask_apply(i,:)) z(i,:) = z(i+1,:)
        end do 
        do i = nx-nr, nx
            where(z(i,:) .eq. missing_value .and. mask_apply(i,:)) z(i,:) = z(i-1,:)
        end do 
        do j = 1+nr, 1, -1
            where(z(:,j) .eq. missing_value .and. mask_apply(:,j)) z(:,j) = z(:,j+1)
        end do 
        do j = ny-nr, ny
            where(z(:,j) .eq. missing_value .and. mask_apply(:,j)) z(:,j) = z(:,j-1)
        end do 


!         if (count(z .eq. missing_value .and. mask_apply) .gt. 0) then 
!             where(z .eq. missing_value .and. mask_apply) &
                             ! z = minval(z,mask=z .ne. missing_value .and. mask_apply)
!         end if 

!         write(*,*) "Fill iterations: ",q 

        return
    end subroutine fill_nearest_dble

    subroutine fill_bilinear_dble(x,y,z,missing_value,fill_value,cont)
        ! Fill the missing values of an array with a bilinear interpolation
        ! Extrapolate when needed 
        
        implicit none 
        double precision, intent(IN)    :: x(:), y(:)
        double precision, intent(INOUT) :: z(:,:)
        double precision :: missing_value 
        double precision, optional :: fill_value
        logical, optional :: cont 

        integer :: nr 
        integer :: q, nx, ny, i, j

        double precision :: x1(3), x2(3), x3(3), x4(3)
        double precision :: xout(2) 
        integer :: i1, i2, i3, i4 
        integer :: j1, j2, j3, j4 
        integer :: ii(2)


        nx = size(z,1)
        ny = size(z,2) 

        ! Loop over array and interpolate missing values 
        do j = 1, ny 
            do i = 1, nx

                if (z(i,j) .eq. missing_value) then 
                    ! Interpolation needed for this point

                    ! Left index  
                    do i1 = i, 1, -1 
                        if (z(i1,j) .ne. missing_value) exit 
                    end do 

                    ! Right index  
                    do i2 = i, nx 
                        if (z(i2,j) .ne. missing_value) exit 
                    end do 

                    ! Down index  
                    do j1 = j, 1, -1 
                        if (z(i,j1) .ne. missing_value) exit 
                    end do 

                    do j2 = j, ny 
                        if (z(i,j2) .ne. missing_value) exit 
                    end do 

                    x1(1) = x(i1)
                    x1(2) = y(j)
                    x1(3) = z(i1,j) 

                    x2(1) = x(i2)
                    x2(2) = y(j)
                    x2(3) = z(i2,j) 

                    x3(1) = x(i)
                    x3(2) = y(j1)
                    x3(3) = z(i,j1) 

                    x4(1) = x(i)
                    x4(2) = y(j2)
                    x4(3) = z(i,j2) 

                    if (x1(3) .eq. missing_value .and. &
                        x2(3) .ne. missing_value) then 
                        x1(3) = x2(3)
                    else if (x2(3) .eq. missing_value .and. &
                             x1(3) .ne. missing_value) then 
                        x2(3) = x1(3)
                    end if 

                    if (x3(3) .eq. missing_value .and. &
                        x4(3) .ne. missing_value) then 
                        x3(3) = x4(3)
                    else if (x4(3) .eq. missing_value .and. &
                             x3(3) .ne. missing_value) then 
                        x4(3) = x3(3)
                    end if 
                    
                    if (x1(3) .ne. missing_value .and. &
                        x3(3) .ne. missing_value) then 
                        
                        z(i,j) = calc_bilinear(x1,x2,x3,x4,xout)

                    end if 

                end if

            end do 
        end do 

        return

    end subroutine fill_bilinear_dble

!     subroutine fill_bilinear_dble(x,y,z,missing_value,fill_value,cont)
!         ! Fill the missing values of an array with a bilinear interpolation
!         ! Extrapolate when needed 

!         ! ajr: TO DO: NOT FINISHED !!!
        
!         implicit none 
!         double precision :: x(:), y(:), z(:,:)
!         double precision :: missing_value 
!         double precision, optional :: fill_value
!         logical, optional :: cont 

!         integer :: nr 
!         integer :: q, nx, ny, i, j

!         double precision :: x1(3), x2(3), x3(3), x4(3)
!         double precision :: xout(2) 
!         integer :: i1, i2, i3, i4 
!         integer :: j1, j2, j3, j4 
!         integer :: ii(2)


!         nx = size(z,1)
!         ny = size(z,2) 

!         ! Loop over array and interpolate missing values 
!         do j = 1, ny 
!             do i = 1, nx

!                 if (z(i,j) .eq. missing_value) then 
!                     ! Interpolation needed for this point

!                     ! Get horizontal bracketing indices
!                     ii = get_nn_indices_line(z(:,j),i,missing_value,cont)
!                     x1(1) = x(ii(1))
!                     x1(2) = y(j)
!                     x1(3) = z(ii(1),j)
!                     x2(1) = x(ii(2))
!                     x2(2) = y(j)
!                     x2(3) = z(ii(2),j)
                    
!                     z(i,j) = calc_bilinear(x1,x2,x3,x4,xout)

!                 end if

!             end do 
!         end do 


!         return

!     end subroutine fill_bilinear_dble

    function get_nn_indices_line(z,i,missing_value,cont) result(ii)
        ! Find the indices of the nearest bracketing neighbors
        ! for point x(i) that are not missing.
        ! cont = .TRUE. if x is continous such that x(1) and x(n) 
        ! are neighbors
        implicit none 

        double precision  :: z(:) 
        integer           :: i 
        double precision  :: missing_value 
        logical, optional :: cont 
        logical           :: is_continuous
        integer :: ii(2) 

        ! Local variables 
        double precision  :: x(size(z)), xtmp(size(z)), ztmp(size(z))
        integer :: nx, n, j, iitmp(2)

        ! Determine if x is continous
        is_continuous = .FALSE.
        if (present(cont)) is_continuous = cont 

        ! Length of vector
        nx = size(x)

        ! Populate x-values 
        do j = 1, nx 
            x(j) = dble(j) 
        end do 

        ! First assume no neighbors without missing values are available
        iitmp = -1 

        ! === Get left index ===

        ! Determine how many neighbors to the left of the point to check
        n = i-1 
        if (is_continuous) n = nx 

        ! Shift vector so that point i is the last value in the vector
        xtmp = cshift(x,shift=nx-i+1)
        ztmp = cshift(z,shift=nx-i+1)

        ! If point i is not the left-most point, then check left-neighbors
        if (n .gt. 0) then
            do j = 1, n 
                if (ztmp(nx-j) .ne. missing_value) then 
                    iitmp(1) = j 
                    exit 
                end if
            end do  
        end if 

        ! === Get right index ===

        ! Determine how many neighbors to the right of the point to check
        n = nx-i 
        if (is_continuous) n = nx 

        ! Shift vector so that point i is the first value in the vector
        xtmp = cshift(x,shift=-i+1)
        ztmp = cshift(z,shift=-i+1)

        ! If point i is not the right-most point, then check right-neighbors
        if (n .gt. 0) then
            do j = 1, n 
                if (ztmp(1+j) .ne. missing_value) then 
                    iitmp(2) = j 
                    exit 
                end if
            end do  
        end if 
        
        ! Indices correspond to the shifted vectors
        ! Now find the indices of the actual vectors
        ii = iitmp 
        if (ii(1) .gt. 0) ii(1) = minloc(abs(x-xtmp(iitmp(1))),dim=1)
        if (ii(2) .gt. 0) ii(2) = minloc(abs(x-xtmp(iitmp(2))),dim=1)
        
        return 

    end function get_nn_indices_line


    subroutine fill_weighted_dble(z,missing_value,fill_value,n,mask)
        implicit none 
        double precision, intent(INOUT) :: z(:,:) 
        double precision, intent(IN)    :: missing_value 
        double precision, intent(IN), optional :: fill_value
        integer,          intent(IN), optional :: n 
        logical,          intent(IN), optional :: mask(:,:) 

        ! Local variables 
        integer :: nr 
        integer :: q, nx, ny, i, j
        integer, parameter :: qmax = 400 ! Iterations 
!         integer, parameter :: nr   = 1  ! Neighbor radius
        
        double precision, dimension (:,:), allocatable :: neighb, weight, weight0 
        double precision :: wtot, mval 
        double precision, allocatable :: filled(:,:)
        integer, allocatable :: nquad(:,:) 
        logical, allocatable :: mask_apply(:,:) 

        integer :: quadmin 

        nr = 4
        if (present(n)) nr = n 

        nx = size(z,1)
        ny = size(z,2) 

        allocate(filled(nx,ny))
        allocate(nquad(nx,ny))
        allocate(neighb(2*nr+1,2*nr+1),weight(2*nr+1,2*nr+1),weight0(2*nr+1,2*nr+1))

        allocate(mask_apply(nx,ny)) 

        mask_apply = .TRUE. 
        if (present(mask)) mask_apply = mask 

        if (present(fill_value)) then
            where(z .eq. missing_value) z = fill_value 
        end if 

        ! Define the default neighbor weighting 
        ! All weights == 1 incorporates some numerical diffusion to make
        ! the resulting extrapolated values smoother
        weight0 = 1
        do i = 1, 2*nr+1
        do j = 1, 2*nr+1 
            weight0(i,j) = 1/max(1d-1,dsqrt(dble((nr+1)-i)**2 + dble((nr+1)-j)**2))
        end do 
        end do  
!         write(*,"(10f12.3)") weight0(nr+1,:)

        quadmin = 4 

        do q = 1, qmax 

            filled = missing_value 
            nquad  = 0

            do i = 1+nr, nx-nr
                do j = 1+nr, ny-nr
                    
                    if (z(i,j) .eq. missing_value) then 

                        neighb = z(i-nr:i+nr,j-nr:j+nr)

                        weight = weight0
                        where (neighb .eq. missing_value) weight = 0.d0
                        wtot = sum(weight)

                        nquad(i,j) = count_quadrants(weight .gt. 0)
                        if (nquad(i,j) .ge. quadmin) filled(i,j) = sum(neighb*weight)/wtot
!                         write(*,*) i,j, nquad(i,j)

                        ! Total neighbors should be 9*nr, so only fill points
                        ! with at least 3*nr valid neighbors (1/3)
!                         if (count(weight .gt. 0.d0) .ge. 1) filled(i,j) = sum(neighb*weight)/wtot

                    end if 
                end do 
            end do 

            where(filled .ne. missing_value) z = filled 
            if ( count(z(1+nr+1:nx-(nr+1),1+nr+1:ny-(nr+1)) .eq. missing_value) .eq. 0 ) exit 

            quadmin = maxval(nquad)

!             write(*,*) "Still missing... ", count(z(1+nr:nx-nr,1+nr:ny:nr) .eq. missing_value), &
!                         " of ", nx*ny
        end do 

        if (q .ge. qmax) then 
            write(*,*) "Too many iterations to fill array of size: ", nx, ny 
            stop 
        end if 

        ! Fill in boundaries too 
        do i = 1+nr, 1, -1
            where(z(i,:) .eq. missing_value) z(i,:) = z(i+1,:)
        end do 
        do i = nx-nr, nx
            where(z(i,:) .eq. missing_value) z(i,:) = z(i-1,:)
        end do 
        do j = 1+nr, 1, -1
            where(z(:,j) .eq. missing_value) z(:,j) = z(:,j+1)
        end do 
        do j = ny-nr, ny
            where(z(:,j) .eq. missing_value) z(:,j) = z(:,j-1)
        end do 


!         if (count(z .eq. missing_value) .gt. 0) then 
!             where(z .eq. missing_value) z = minval(z,mask=z .ne. missing_value)
!         end if 

!         write(*,*) "Fill iterations: ",q 

        return
    end subroutine fill_weighted_dble

    ! Fill in missing values of an array with neighbor averages
    ! or with a specified fill_value
    subroutine fill_mean_dble(var,missing_value,fill_value,mask)
        implicit none 
        double precision, intent(INOUT) :: var(:,:) 
        double precision, intent(IN)    :: missing_value 
        double precision, intent(IN), optional :: fill_value
        logical,          intent(IN), optional :: mask(:,:) 

        integer :: q, nx, ny, i, j 
        integer, parameter :: qmax = 50 ! Iterations 

        double precision, dimension (3,3) :: neighb, weight
        double precision :: wtot, mval 
        double precision, allocatable :: filled(:,:)
        logical, allocatable :: mask_apply(:,:)

        nx = size(var,1)
        ny = size(var,2) 

        allocate(filled(nx,ny))
        allocate(mask_apply(nx,ny)) 

        mask_apply = .TRUE. 
        if (present(mask)) mask_apply = mask 

        if (present(fill_value)) then
            where(var .eq. missing_value) var = fill_value 
        end if 

        do q = 1, qmax 

            filled = missing_value 

            do i = 2, nx-1 
                do j = 2, ny-1 
                    
                    if (count(mask_apply(i-1:i+1,j-1:j+1)) .gt. 1) then 

                        neighb = var(i-1:i+1,j-1:j+1)

                        weight = 0.d0 
                        where (neighb .ne. missing_value) weight = 1.d0
                        wtot = sum(weight)

                        if (wtot .gt. 0.d0) then 
                            mval = sum(neighb*weight)/wtot
                            where (neighb .eq. missing_value) neighb = mval 
                        end if 

                        filled(i-1:i+1,j-1:j+1) = neighb 

                    end if 

                end do 
            end do 

            where(filled .ne. missing_value .and. mask_apply) var = filled 

!             write(*,*) q," : Missing values: ", count(var .eq. missing_value)
            if ( count(var .eq. missing_value  .and. mask_apply) .eq. 0 ) exit 
        end do 

        return
    end subroutine fill_mean_dble

    function count_quadrants(not_missing) result(n)
        ! Assuming the origin in the center of array 'not_missing',
        ! count how many quadrants contain valid neighbors

        implicit none 

        logical, dimension(:,:) :: not_missing 
        integer :: quadrants(4)
        integer :: n 
        integer :: nr, nx

        nx = size(not_missing,1)
        nr = (nx-1) / 2

        quadrants =  [count(not_missing(nr+1:nx,nr+1:nx)), &
                      count(not_missing(1:nr,nr+1:nx)), &
                      count(not_missing(1:nr,1:nr)), &
                      count(not_missing(nr+1:nx,1:nr))]
        
        n = count(quadrants .gt. 0)

        return 

    end function count_quadrants

    
    ! Interface alternatives that use the above main routines 

    function interp_nearest_int(x,y,z,xout,yout,missing_value,mask) result(zout)
        ! Find closest x-indices and closest y-indices on original
        ! grid (assume these and next indices will bracket our point)
        ! Pick closest point 
        implicit none 

        real(dp), dimension(:) :: x, y, xout, yout 
        integer, dimension(:,:) :: z
        integer, optional :: missing_value 
        logical, dimension(:,:), optional :: mask 
        integer, dimension(size(xout,1),size(yout,1)) :: zout 
        logical,  dimension(size(xout,1),size(yout,1)) :: mask_interp 

        zout = nint(interp_nearest_dble(x,y,dble(z),xout,yout,dble(missing_value),mask))

        return 

    end function interp_nearest_int

    subroutine fill_nearest_int(z,missing_value,fill_value,mask,n)

        implicit none 
        
        integer :: z(:,:) 
        integer :: missing_value 
        double precision, optional :: fill_value
        logical, optional :: mask(:,:) 
        integer, optional :: n

        double precision, dimension(size(z,1),size(z,2)) :: z_dble 

        z_dble = dble(z) 
        call fill_nearest_dble(z_dble,dble(missing_value),fill_value,mask,n)
        z = nint(z_dble)

        return 

    end subroutine fill_nearest_int

    subroutine fill_nearest_float(z,missing_value,fill_value,mask,n)

        implicit none 
        
        real(sp) :: z(:,:) 
        real(sp) :: missing_value 
        double precision, optional :: fill_value
        logical, optional :: mask(:,:) 
        integer, optional :: n

        double precision, dimension(size(z,1),size(z,2)) :: z_dble 

        z_dble = dble(z) 
        call fill_nearest_dble(z_dble,dble(missing_value),fill_value,mask,n)
        z = real(z_dble,sp)

        return 

    end subroutine fill_nearest_float

    subroutine diffuse(z,iter,missing_value,mask)

        implicit none 

        double precision :: z(:,:)
        logical, optional :: mask(:,:)
        double precision :: missing_value 
        integer :: iter, nx, ny, q, i, j 
        double precision :: ztmp(size(z,1),size(z,2))
        logical :: mask_tmp(3,3)

        nx = size(z,1)
        ny = size(z,2)

        do q = 1, iter 

            ztmp = z 
            do i = 2, nx-1 
            do j = 2, ny-1 
                mask_tmp = z(i-1:i+1,j-1:j+1) .ne. missing_value 
                if (count(mask_tmp) .gt. 0) ztmp(i,j) = sum(z(i-1:i+1,j-1:j+1),mask_tmp) / count(mask_tmp)
            end do 
            end do 

            ! Borders 
            j = 1
            do i = 1, nx 
                mask_tmp(1,1:2) = z(i,j:j+1) .ne. missing_value 
                if (count(mask_tmp(1,1:2)) .gt. 0) ztmp(i,j) = sum(z(i,j:j+1),mask_tmp(1,1:2)) / count(mask_tmp(1,1:2))
            end do 

            j = ny
            do i = 1, nx 
                mask_tmp(1,1:2) = z(i,j-1:j) .ne. missing_value 
                if (count(mask_tmp(1,1:2)) .gt. 0) ztmp(i,j) = sum(z(i,j-1:j),mask_tmp(1,1:2)) / count(mask_tmp(1,1:2))
            end do 

            i = 1
            do j = 1, ny 
                mask_tmp(1:2,1) = z(i:i+1,j) .ne. missing_value 
                if (count(mask_tmp(1:2,1)) .gt. 0) ztmp(i,j) = sum(z(i:i+1,j),mask_tmp(1:2,1)) / count(mask_tmp(1:2,1))
            end do 

            i = nx
            do j = 1, ny 
                mask_tmp(1:2,1) = z(i-1:i,j) .ne. missing_value 
                if (count(mask_tmp(1:2,1)) .gt. 0) ztmp(i,j) = sum(z(i-1:i,j),mask_tmp(1:2,1)) / count(mask_tmp(1:2,1))
            end do 


            z = ztmp 

        end do 

        return 

    end subroutine diffuse

    subroutine limit_gradient(z,dx,dy,grad_lim,iter_max,mask)
        ! Limit the gradient to below a threshold 

        implicit none 

        real(8), intent(INOUT) :: z(:,:)
        real(8), intent(IN)    :: dx, dy 
        real(8), intent(IN)    :: grad_lim 
        integer, intent(IN)    :: iter_max 
        logical, intent(IN), optional :: mask(:,:)

        ! Local variables 
        real(8), allocatable :: z0(:,:), dz(:,:)
        logical, allocatable :: mask_local(:,:) 

        integer :: i, j, nx, ny, k  
        real(8) :: hgrad(4)
        integer :: q 

        ! Get dimensions of z
        nx = size(z,1)
        ny = size(z,2)

        ! Allocate new z and store old z initially
        allocate(z0(nx,ny))
        allocate(dz(nx,ny))
        dz = 0.d0 

        ! Allocate mask and populate it 
        allocate(mask_local(nx,ny))
        mask_local = .TRUE. 
        if (present(mask)) mask_local = mask 

        ! Iterate until now gradient limits exceeded 
        do q = 1, iter_max 

            ! Store current array in old array 
            z0 = z 

            ! Loop over z, limit gradient as desired 
            do i = 2, nx-1 
                do j = 2, ny-1 

                    if (mask_local(i,j)) then 
                        ! Perform check on points of interest only 

                        ! Find neighbor of maximum gradient
                        ! (height change rel. to current point, not actual gradient)
                        hgrad(1) = (z0(i-1,j)-z0(i,j))/dx
                        hgrad(2) = (z0(i+1,j)-z0(i,j))/dx
                        hgrad(3) = (z0(i,j-1)-z0(i,j))/dy
                        hgrad(4) = (z0(i,j+1)-z0(i,j))/dy

                        k = maxloc(abs(hgrad),dim=1)
                        dz(i,j) = hgrad(k) 

                        if (abs(hgrad(k)) .gt. grad_lim) then 
                            ! Apply gradient limit to point 
                            ! note: only correct by half of gradient limit to avoid overshooting,
                            !       since correction will likely be applied in both directions

                            select case(k)
                                case(1) 
                                    z(i,j) = z0(i-1,j)-sign(grad_lim*0.8d0,hgrad(k))*dx
                                case(2) 
                                    z(i,j) = z0(i+1,j)-sign(grad_lim*0.8d0,hgrad(k))*dx
                                case(3) 
                                    z(i,j) = z0(i,j-1)-sign(grad_lim*0.8d0,hgrad(k))*dy
                                case(4) 
                                    z(i,j) = z0(i,j+1)-sign(grad_lim*0.8d0,hgrad(k))*dy
                            end select 

                        end if

                    end if 

                end do 
            end do 

            ! Output diagnostics 
            write(*,"(a,i3,2g10.3)") "Gradient iteration, range dz: ", q, minval(dz), maxval(dz)

            ! If hgrad is below limit, exit iterative loop 
            if (maxval(abs(dz)) .le. grad_lim) exit 

            
        end do 

        return 

    end subroutine limit_gradient

    function calc_bilinear(x1,x2,x3,x4,xout) result(zout)
        ! Given the points surrounding it, calculate the value at xout
        ! x3  p2  x4
        !    xout
        ! x1  p1  x2

        implicit none 

        double precision :: x1(3), x2(3), x3(3), x4(3)   ! [x,y,z]
        double precision :: xout(2)   ! [x,y]
        double precision :: zout 

        ! Local variables
        double precision :: alpha1, alpha2, p1, p2 

        alpha1 = (xout(1) - x1(1)) / (x2(1)-x1(1))  ! x-values
        p1     = x1(3) + alpha1*(x2(3)-x1(3))
        p2     = x3(3) + alpha1*(x4(3)-x3(3))

        alpha2 = (xout(2)-x1(2))/(x3(2)-x1(2))      ! y-values
        zout   = p1 + alpha2*(p2-p1)
        
        return

    end function calc_bilinear





! ========================
    
    subroutine fill_poisson_float(var2d,missing_value,method,wraplon,verbose)

        implicit none 

        real(4),          intent(INOUT) :: var2d(:,:) 
        real(4),          intent(IN)    :: missing_value 
        integer,          intent(IN)    :: method   ! 0: start with guess=0.0, 1: guess zonal mean, 2: grid mean
        logical,          intent(IN)    :: wraplon  ! Cyclical grid (wrap longitude?)
        logical,          intent(IN), optional :: verbose 

        ! Local variables 
        integer :: nx, ny 
        double precision, allocatable :: var2d_dble(:,:) 

        nx = size(var2d,1)
        ny = size(var2d,2) 

        allocate(var2d_dble(nx,ny))

        var2d_dble = var2d 
        call fill_poisson_dble(var2d_dble,dble(missing_value),method,wraplon,verbose)
        var2d = var2d_dble

        return 

    end subroutine fill_poisson_float

    subroutine fill_poisson_dble(var2d,missing_value,method,wrapx,verbose)

        ! Poisson fill code originally obtained from:
        ! https://github.com/royalosyin/Fill-Missing-Values/blob/master/poisson_grid_fill/fill_msg_grid.f90
        ! Originally ported from NCL? 

        implicit none 

        double precision, intent(INOUT) :: var2d(:,:) 
        double precision, intent(IN)    :: missing_value 
        integer,          intent(IN)    :: method   ! 0: start with guess=0.0, 1: guess zonal mean, 2: grid mean
        logical,          intent(IN)    :: wrapx    ! Cyclical grid (wrap longitude?)
        logical,          intent(IN), optional :: verbose 

        ! Local variables 
        integer :: nx, ny, n, j  
        logical, allocatable :: mask(:,:) 
        double precision     :: varmean
        
        double precision, parameter :: tol = 1d-2 

        nx = size(var2d,1)
        ny = size(var2d,2) 

        allocate(mask(nx,ny)) 
        mask = var2d .eq. missing_value

        if (count(mask) .gt. 0) then 
            
            select case(method)

                case(0) 
                    ! Fill missing values with zeros 
                    where (var2d .eq. missing_value) var2d = 0.0d0 

                case(1) 
                    ! Fill missing values with zonal mean 
                    do j = 1, ny 
                        n = count(.not. mask(:,j))
                        if (n .gt.0) then 
                            varmean = sum(var2d(:,j),mask=(.not. mask(:,j))) / dble(n)
                        else 
                            varmean = 0.0d0 
                        end if 
                        where(mask(:,j)) var2d(:,j) = varmean 
                    end do 

                case(2) 
                    ! Fill missing values with grid mean 
                    n = count(.not.mask)
                    if (n .gt. 0) then 
                        varmean = sum(var2d,mask=(.not.mask)) / dble(n)
                    else 
                        varmean = 0.0d0 
                    end if 
                    where (mask) var2d = varmean 

                case(3)
                    ! Fill with neighborhood mean 

                    !call fill_mean_dble(var2d,missing_value=missing_value)
                    call fill_weighted_dble(var2d,missing_value=missing_value,n=6)
                    !call fill_nearest_dble(var2d,missing_value=missing_value)
                    
                case DEFAULT 

                    write(*,*) "fill_poisson:: Error: guess method not recognized."

            end select 

            call filter_poisson_dble(var2d,mask,tol,missing_value, &
                                            wrapx,verbose,filling=.TRUE.)

        end if 

        return 

    end subroutine fill_poisson_dble

    subroutine filter_poisson_float(var2d,mask,tol,missing_value,wrapx,verbose,filling)

        implicit none 

        real(4), intent(INOUT) :: var2d(:,:) 
        logical, intent(IN)    :: mask(:,:) 
        real(4), intent(IN)    :: tol 
        real(4), intent(IN)    :: missing_value 
        logical, intent(IN)    :: wrapx 
        logical, intent(IN), optional :: verbose 
        logical, intent(IN), optional :: filling 

        ! Local variables 
        real(8), allocatable :: var2d_dble(:,:) 

        allocate(var2d_dble(size(var2d,1),size(var2d,2)))

        var2d_dble = var2d 
        call filter_poisson_dble(var2d_dble,mask,dble(tol),dble(missing_value),wrapx,verbose,filling)
        var2d = var2d_dble 

        return 

    end subroutine filter_poisson_float

    subroutine filter_poisson_dble(var2d,mask,tol,missing_value,wrapx,verbose,filling)

        implicit none 

        real(8), intent(INOUT) :: var2d(:,:) 
        logical, intent(IN)    :: mask(:,:) 
        real(8), intent(IN)    :: tol               ! Fraction
        real(8), intent(IN)    :: missing_value 
        logical, intent(IN)    :: wrapx 
        logical, intent(IN), optional :: verbose 
        logical, intent(IN), optional :: filling 

        ! Local variables 
        integer :: iter, i, j, nx, ny 
        integer :: im1, ip1, jm1, jp1 
        real(8) :: resid_lim, resid_max, resid
        real(8), allocatable :: var2dnew(:,:) 
        real(8) :: var_check
        real(8) :: var4(4)
        real(8) :: wt4(4) 

        integer, parameter :: iter_max = 1000 ! iter usually 10-20
        real(8), parameter :: rel      = 0.5  ! relaxation coeff (recommended 0.45-0.6)

        nx = size(var2d,1) 
        ny = size(var2d,2) 

        allocate(var2dnew(nx,ny)) 

        var2dnew = var2d 

        ! Calculate resid_lim == ~mean grid value X tol
        resid_lim = tol * &
                    abs(0.5d0*(maxval(var2d,mask=var2d.ne.missing_value) &
                               + minval(var2d,mask=var2d.ne.missing_value)))

        ! Perform iterations 
        do iter = 1, iter_max 

            resid_max = 0.0d0 

            do j = 1, ny
            do i = 1, nx 

                !if (resid(i,j) .gt. resid_lim) then 
                if (mask(i,j)) then 
                    ! Only work on desired locations 

                    jm1 = max(j-1,1) 
                    jp1 = min(j+1,ny)
                    
                    if (j.eq.1 ) jm1 = 2
                    if (j.eq.ny) jp1 = ny-1

                    im1 = i-1
                    ip1 = i+1

                    ! cyclic in x           
                    if (i.eq.1  .and. wrapx) im1 = nx 
                    if (i.eq.nx .and. wrapx) ip1 = 1
                    
                    ! not cyclic in x           
                    if (i.eq.1  .and. (.not. wrapx)) im1 = 2
                    if (i.eq.nx .and. (.not. wrapx)) ip1 = nx-1 
                
                    var4 = [var2d(im1,j),var2d(ip1,j),var2d(i,jm1),var2d(i,jp1)]
                    wt4  = 1.0d0 
                    where(var4 .eq. missing_value) wt4 = 0.0d0 
                    var_check = sum(wt4*var4) / sum(wt4) 

                    resid         = (var_check-var2d(i,j))
                    var2dnew(i,j) = var2d(i,j) + resid*rel
                    resid_max     = max(abs(resid),resid_max)

                end if 

            end do 
            end do 

            ! Update var2d 
            var2d = var2dnew 

            ! If residual is low enough, exit iterations
            if (resid_max .le. resid_lim) exit 

        end do 

        ! If verbose, print summary
        if (present(verbose) .and. present(filling)) then 
            if (verbose .and. filling) then 
                write(*,"(a35,i5,g12.3)") "fill_poisson: iter, resid_max = ", iter, resid_max 
            end if 
        else if (present(verbose)) then 
            if (verbose) then 
                write(*,"(a35,i5,g12.3)") "filter_poisson: iter, resid_max = ", iter, resid_max 
            end if 
        end if 

        return 

    end subroutine filter_poisson_dble

    subroutine filter_poisson_dble_new(var2d,mask,dx,iter_max,tol,rel,missing_value,wrapx,verbose,filling)

        implicit none 

        real(8), intent(INOUT) :: var2d(:,:) 
        logical, intent(IN)    :: mask(:,:) 
        real(8), intent(IN)    :: dx 
        integer, intent(IN)    :: iter_max 
        real(8), intent(IN)    :: tol               ! Gradient dvar/dx
        real(8), intent(IN)    :: rel
        real(8), intent(IN)    :: missing_value 
        logical, intent(IN)    :: wrapx 
        logical, intent(IN), optional :: verbose 
        logical, intent(IN), optional :: filling 

        ! Local variables 
        integer :: iter, i, j, nx, ny 
        integer :: im1, ip1, jm1, jp1 
        real(8) :: resid_lim, resid_max 
        real(8), allocatable :: resid(:,:) 
        real(8), allocatable :: var2dnew(:,:) 
        real(8) :: var_check
        real(8) :: var4(4)
        real(8) :: wt4(4) 

        nx = size(var2d,1) 
        ny = size(var2d,2) 

        allocate(resid(nx,ny))
        allocate(var2dnew(nx,ny)) 

        var2dnew = var2d 

        ! Calculate resid_lim == ~mean grid value X tol
        resid_lim = tol * &
                    abs(0.5d0*(maxval(var2d,mask=var2d.ne.missing_value) &
                               + minval(var2d,mask=var2d.ne.missing_value)))
        ! resid_lim = tol 

        resid = resid_lim*10.0 
        where(.not. mask) resid = 0.0 

        ! Perform iterations 
        do iter = 1, iter_max 

            resid_max = 0.0d0 

            do j = 1, ny
            do i = 1, nx 

                ! if (resid(i,j) .gt. resid_lim) then 
                if (mask(i,j)) then 
                    ! Only work on desired locations 

                    jm1 = max(j-1,1) 
                    jp1 = min(j+1,ny)
                    
                    if (j.eq.1 ) jm1 = 2
                    if (j.eq.ny) jp1 = ny-1

                    im1 = i-1
                    ip1 = i+1

                    ! cyclic in x           
                    if (i.eq.1  .and. wrapx) im1 = nx 
                    if (i.eq.nx .and. wrapx) ip1 = 1
                    
                    ! not cyclic in x           
                    if (i.eq.1  .and. (.not. wrapx)) im1 = 2
                    if (i.eq.nx .and. (.not. wrapx)) ip1 = nx-1 
                
                    var4 = [var2d(im1,j),var2d(ip1,j),var2d(i,jm1),var2d(i,jp1)]
                    wt4  = 1.0d0 
                    where(var4 .eq. missing_value) wt4 = 0.0d0 
                    var_check = sum(wt4*var4) / sum(wt4) 

                    resid(i,j) = (var_check-var2d(i,j))! / dx 

                    if (resid(i,j) .gt. resid_lim) then 
                        var2dnew(i,j) = var2d(i,j) + rel*resid(i,j)
                    end if 

                end if 

            end do 
            end do 

            ! Get max resid value 
            resid_max = maxval(resid)

            ! Update var2d 
            var2d = var2dnew 

            ! If residual is low enough, exit iterations
            if (resid_max .le. resid_lim) exit 

        end do 

        ! If verbose, print summary
        if (present(verbose) .and. present(filling)) then 
            if (verbose .and. filling) then 
                write(*,"(a35,i5,g12.3)") "fill_poisson: iter, resid_max = ", iter, resid_max 
            end if 
        else if (present(verbose)) then 
            if (verbose) then 
                write(*,"(a35,i5,g12.3)") "filter_poisson: iter, resid_max = ", iter, resid_max 
            end if 
        end if 

        return 

    end subroutine filter_poisson_dble_new

! =====================

end module interp2D