program test_levelset

    use yelmo 
    use ncio 

    implicit none

    type(yelmo_class)  :: yelmo1

    character(len=256) :: infldr
    character(len=256) :: outfldr
    character(len=512) :: path_par
    character(len=512) :: path_const
    character(len=256) :: file2D
    character(len=256) :: file1D
    character(len=256) :: file_restart

    character(len=56)  :: domain
    character(len=56)  :: grid_name  

    integer  :: i, j, k 
    real(wp) :: dx, dz 
    real(wp) :: xmax, ymin, ymax, zmin, zmax 

    type levelset_class

        ! General information 
        real(wp) :: missing_value

        ! Axis information 
        real(wp) :: dx, dy, dz
        integer  :: nx, ny, nz  
        real(wp), allocatable :: x(:)
        real(wp), allocatable :: y(:)
        real(wp), allocatable :: z(:)
        
        ! Variables 
        real(wp), allocatable :: H_ice(:,:) 
        real(wp), allocatable :: z_srf(:,:)
        real(wp), allocatable :: z_base(:,:)
        real(wp), allocatable :: z_bed(:,:)
        real(wp), allocatable :: z_sl(:,:)

        real(wp), allocatable :: ATT(:,:,:) 

        real(wp), allocatable :: u(:,:,:)
        real(wp), allocatable :: v(:,:,:)
        real(wp), allocatable :: w(:,:,:)

        real(wp), allocatable :: phi(:,:,:)
        integer,  allocatable :: G(:,:,:)
        
    end type

    type(levelset_class) :: lev1 

    real(wp) :: time_init 
    real(wp) :: time_end 
    real(wp) :: dt 
    integer  :: n, nt 
    real(wp) :: time 

    real(8) :: cpu_start_time, cpu_end_time, cpu_dtime  

    ! Start timing 
    call yelmo_cpu_time(cpu_start_time)

    ! Define default grid name for completeness 
    domain    = "EISMINT"
    grid_name = "EISMINT" 

    ! Assume program is running from the output folder
    infldr  = "./"
    outfldr = "./output/levelset/"

    ! Define input and output locations 
    path_const = trim(infldr)//"par/yelmo_const_"//trim(domain)//".nml"
    path_par   = trim(infldr)//"par/yelmo_"//trim(domain)//".nml"
    file2D     = trim(outfldr)//"lev2D.nc"
    file1D     = trim(outfldr)//"lev1D.nc"
    
    
    ! General initialization of yelmo constants (used globally)
    call yelmo_global_init(path_const)


    ! Define the domain and horizontal grid [km]
    dx   =  5.0 
    xmax = 1000.0
    ymax =  0.0
    ymin =  0.0
    call yelmo_init_grid(yelmo1%grd,grid_name,units="km",x0=0.0,dx=dx,nx=int(xmax/dx)+1,y0=ymin,dy=dx,ny=int((ymax-ymin)/dx)+1)

    ! Define vertical coordinates [m]
    zmin = 0.0 
    zmax = 5000.0 
    dz   = 50.0 
    
    ! Initialize levelset object 
    call levelset_init(lev1,yelmo1%grd%nx,yelmo1%grd%ny,nz=int((zmax-zmin)/dz)+1)

    ! Store axis information
    lev1%dx = dx 
    lev1%x  = yelmo1%grd%xc 

    lev1%dy = dx 
    lev1%y  = yelmo1%grd%yc 
    
    lev1%dz = dz 
    do k = 1, lev1%nz 
        lev1%z(k) = zmin + (k-1)*lev1%dz 
    end do

    ! Summary of grid 
    write(*,*) 
    write(*,*) "====="
    write(*,*) "x: ", lev1%nx, lev1%dx, minval(lev1%x), maxval(lev1%x)
    write(*,*) "y: ", lev1%ny, lev1%dy, minval(lev1%y), maxval(lev1%y)
    write(*,*) "z: ", lev1%nz, lev1%dz, minval(lev1%z), maxval(lev1%z)
    
    time_init = 100.0_wp 
    time_end  = 100.0_wp 
    dt        = 10.0_wp 
    nt        = ceiling( (time_end-time_init) / dt ) + 1

    ! Set sea level and bed elevation
    lev1%z_bed = 0.0_wp 
    lev1%z_sl  = 0.0_wp 

    ! Set rate factor 
    lev1%ATT   = 1e-16_wp


    ! Initialize output file 
    call levelset_write_init(lev1,file2D,time_init=time_init,units="years")
    
    do n = 1, nt 

        time = time_init + dt*(n-1)

        ! Initialize Halfar ice sheet profile 
        call calc_halfar(lev1%H_ice,lev1%x,lev1%y,time=time,R0=750.0_wp,H0=3600.0_wp, &
                                lambda=0.0_wp,n=3.0_wp,A=1e-16_wp,rho_ice=rho_ice,g=g)

        if (n .eq. 1) then 
            ! Calculate the initial levelset function
            call levelset_phi_set(lev1,rho_ice,rho_sw)
        else 
            ! Update levelset (TO DO)

        end if 

        ! Calculate signed distance everywhere
        call levelset_calc_signed_distance(lev1)

        ! Calculate SIA velocity profile too 
        call calc_vel_sia_2D(lev1%u(:,1,:),lev1%w(:,1,:),lev1%x,lev1%z,lev1%H_ice(:,1), &
                        lev1%z_srf(:,1),lev1%z_bed(:,1),lev1%ATT(:,1,:),rho_ice,g,n=3.0_wp)

        ! Initialize and write output
        call levelset_write_step(lev1,file2D,time=time)

        write(*,*) time, "H: ", minval(lev1%H_ice), maxval(lev1%H_ice)
        write(*,*) time, "u: ", minval(lev1%u), maxval(lev1%u)
        write(*,*) time, "w: ", minval(lev1%w), maxval(lev1%w)

    end do 

contains
    
    subroutine levelset_init(lev,nx,ny,nz)

        implicit none 

        type(levelset_class), intent(INOUT) :: lev 
        integer, intent(IN) :: nx 
        integer, intent(IN) :: ny 
        integer, intent(IN) :: nz 
        
        ! === Axes ===

        if (allocated(lev%x)) deallocate(lev%x)
        if (allocated(lev%y)) deallocate(lev%y)
        if (allocated(lev%z)) deallocate(lev%z)

        allocate(lev%x(nx))
        allocate(lev%y(ny))
        allocate(lev%z(nz))

        lev%nx = size(lev%x)
        lev%ny = size(lev%y)
        lev%nz = size(lev%z)

        ! === Variables ===

        if (allocated(lev%H_ice))   deallocate(lev%H_ice)
        if (allocated(lev%z_srf))   deallocate(lev%z_srf)
        if (allocated(lev%z_base))  deallocate(lev%z_base)
        if (allocated(lev%z_bed))   deallocate(lev%z_bed)
        if (allocated(lev%z_sl))    deallocate(lev%z_sl)
        
        if (allocated(lev%ATT))     deallocate(lev%ATT)
        if (allocated(lev%u))       deallocate(lev%u)
        if (allocated(lev%v))       deallocate(lev%v)
        if (allocated(lev%w))       deallocate(lev%w)

        if (allocated(lev%phi))     deallocate(lev%phi)
        if (allocated(lev%G))       deallocate(lev%G)
        
        allocate(lev%H_ice(nx,ny))
        allocate(lev%z_srf(nx,ny))
        allocate(lev%z_base(nx,ny))
        allocate(lev%z_bed(nx,ny))
        allocate(lev%z_sl(nx,ny))

        allocate(lev%ATT(nx,ny,nz))
        allocate(lev%u(nx,ny,nz))
        allocate(lev%v(nx,ny,nz))
        allocate(lev%w(nx,ny,nz))
        
        allocate(lev%phi(nx,ny,nz))
        allocate(lev%G(nx,ny,nz))
        
        lev%H_ice   = 0.0_wp 
        lev%z_srf   = 0.0_wp 
        lev%z_base  = 0.0_wp 
        lev%z_bed   = 0.0_wp 
        lev%z_sl    = 0.0_wp 
        lev%ATT     = 0.0_wp 
        lev%u       = 0.0_wp 
        lev%v       = 0.0_wp 
        lev%w       = 0.0_wp 
        
        lev%phi     = 1.0_wp 
        lev%G       = 2.0_wp 
        
        lev%missing_value = 1e10_wp

        return

    end subroutine levelset_init


    subroutine levelset_phi_set(lev,rho_ice,rho_sw)

        implicit none

        type(levelset_class), intent(INOUT) :: lev 
        real(wp), intent(IN) :: rho_ice 
        real(wp), intent(IN) :: rho_sw 
        
        ! Local variables
        integer :: i, j, k
        integer :: nx, ny, nz 

        nx = size(lev%phi,1)
        ny = size(lev%phi,2)
        nz = size(lev%phi,3)

        ! Assume H_ice, z_bed and z_sl are well defined. 

        ! First calculate surface elevation 
        call calc_z_srf_max(lev%z_srf,lev%H_ice,lev%z_bed,lev%z_sl,rho_ice,rho_sw)

        ! Now get ice base elevation
        lev%z_base = lev%z_srf-lev%H_ice

        ! Now use H_ice, z_bed and z_srf to define zero-surface
        ! For now set points outside ice body to arbitrary postive/negative
        ! values, later these will be replaced by the signed distance function. 

        do i = 1, nx 
        do j = 1, ny 
                
            do k = 1, nz 

                if (lev%H_ice(i,j) .eq. 0.0_wp) then 
                    ! No ice present in this column, point is outside ice body (phi > 0)

                    lev%phi(i,j,k) = 1.0_wp 

                else if ( lev%z(k) .gt. lev%z_base(i,j) .and. &
                          lev%z(k) .lt. lev%z_srf(i,j) ) then 
                    ! Point is inside ice body (phi < 0)

                    lev%phi(i,j,k) = -1.0_wp 

                else if (lev%z(k) .eq. lev%z_srf(i,j) .or. &
                         lev%z(k) .eq. lev%z_base(i,j)) then 
                    ! Point is on the zero surface (phi = 0)

                    lev%phi(i,j,k) = 0.0_wp 

                else 
                    ! Point is outside the ice body (phi > 0)
                
                    lev%phi(i,j,k) = 1.0_wp 

                end if 

            end do 
                
        end do
        end do

        return

    end subroutine levelset_phi_set

    subroutine levelset_calc_signed_distance(lev)

        implicit none

        type(levelset_class), intent(INOUT) :: lev  
        
        ! Local variables
        integer :: q

        integer, parameter :: q_max = 1000

        ! First specify check matrix, points are either
        ! 0: accepted 
        ! 1: close 
        ! 2: far 

        ! Initially set all points as far
        lev%G = 2 

        ! Now find points that surround the boundary (accepted)
        call find_accepted_cells(lev%G,lev%z_srf,lev%z_base,lev%z)

        ! Assign initial values (phi_0 = 0)
        where(lev%G .eq. 0) lev%phi = 0.0 

        ! Assign high values for other points
        where(lev%G .ne. 0) lev%phi = lev%missing_value

        ! Perform several loops until all points labeled 'far' are eliminated
        do q = 1, q_max

            ! Now find points near the boundary (close)
            call find_close_cells(lev%G)
        
            ! Calculate signed distance of points labeled 'close'
            call levelset_calc_signed_distance_close_points(lev)

            write(*,*) q, count(lev%G .eq. 0), count(lev%G .eq. 2)
            if (maxval(lev%G) .eq. 0) exit

        end do 


        return

    end subroutine levelset_calc_signed_distance

    subroutine levelset_calc_signed_distance_close_points(lev)

        implicit none

        type(levelset_class), intent(INOUT) :: lev 

        ! Local variables 
        integer :: i, j, k, l, m, n
        integer :: nx, ny, nz 
        real(wp) :: phi_temp 

        nx = size(lev%phi,1)
        ny = size(lev%phi,2)
        nz = size(lev%phi,3) 

        do l = 1, nx 
        do m = 1, ny 
        do n = 1, nz 


            if (lev%G(l,m,n) .eq. 1) then 
                
                ! Solve quadratic distance function
                call levelset_solve_quadratic(phi_temp,lev,l,m,n)
                
                if (phi_temp .ne. lev%missing_value) then 

                    ! Assign value
                    lev%phi(l,m,n) = phi_temp 

                    ! Label point as 'accepted'
                    lev%G(l,m,n) = 0 

                end if 

            end if 

        end do
        end do
        end do

        return

    end subroutine levelset_calc_signed_distance_close_points


    subroutine levelset_solve_quadratic(phi_temp,lev,l,m,n)
        ! Following algorithm 2.3 of Yang (preprint, https://arxiv.org/abs/1811.00009v1)

        implicit none 

        real(wp), intent(OUT) :: phi_temp 
        type(levelset_class), intent(IN) :: lev 
        integer, intent(IN) :: l 
        integer, intent(IN) :: m 
        integer, intent(IN) :: n 

        ! Local variables
        integer :: nx, ny, nz
        integer :: lm1, lp1, mm1, mp1, nm1, np1, lpd, mpd, npd
        integer :: id, nd, d

        real(wp) :: phi(3) 
        real(wp) :: dx(3) 
        integer  :: idx(3) 
        real(wp) :: a, b, c 

        nx = size(lev%phi,1)
        ny = size(lev%phi,2)
        nz = size(lev%phi,3) 

        ! Initially, set all values of phi buffer to missing,
        ! dx buffer to zero and and sort order indices to zeros
        phi = 9999.0 
        dx  = 0.0 
        idx = 0

        ! Define neighbor indices 
        lm1 = max(1,l-1)
        lp1 = min(nx,l+1)
        mm1 = max(1,m-1)
        mp1 = min(ny,m+1)
        nm1 = max(1,n-1)
        np1 = min(nz,n+1)
        
        nd = 0 

        ! Check for upwind neighbor, x-direction ===
        
        d = 0 
        if (lev%G(lm1,m,n) .eq. 0 .and. lev%phi(l,m,n) .gt. lev%phi(lm1,m,n)) then 
            d = -1
        end if 
        if (lev%G(lp1,m,n) .eq. 0 .and. lev%phi(l,m,n) .gt. lev%phi(lp1,m,n)) then
            if (d .eq. 0 .or. lev%phi(lp1,m,n) .lt. lev%phi(lm1,m,n)) then 
                d = 1
            end if
        end if 

        if (d .ne. 0) then 
            nd = nd + 1 
            lpd = max(min(l+d,nx),1)
            phi(nd) = lev%phi(lpd,m,n)
            dx(nd)  = abs(lev%x(l)-lev%x(lpd))
        end if 

        ! Check for upwind neighbor, y-direction ===
        
        d = 0 
        if (lev%G(l,mm1,n) .eq. 0 .and. lev%phi(l,m,n) .gt. lev%phi(l,mm1,n)) then 
            d = -1
        end if 
        if (lev%G(l,mp1,n) .eq. 0 .and. lev%phi(l,m,n) .gt. lev%phi(l,mp1,n)) then
            if (d .eq. 0 .or. lev%phi(l,mp1,n) .lt. lev%phi(l,mm1,n)) then 
                d = 1
            end if
        end if 

        if (d .ne. 0) then 
            nd = nd + 1 
            mpd = max(min(m+d,ny),1)
            phi(nd) = lev%phi(l,mpd,n)
            dx(nd)  = abs(lev%y(m)-lev%y(mpd))
        end if 

        ! Check for upwind neighbor, z-direction ===
        
        d = 0 
        if (lev%G(l,m,nm1) .eq. 0 .and. lev%phi(l,m,n) .gt. lev%phi(l,m,nm1)) then 
            d = -1
        end if 
        if (lev%G(l,m,np1) .eq. 0 .and. lev%phi(l,m,n) .gt. lev%phi(l,m,np1)) then
            if (d .eq. 0 .or. lev%phi(l,m,np1) .lt. lev%phi(l,m,nm1)) then 
                d = 1
            end if
        end if 

        if (d .ne. 0) then 
            nd = nd + 1 
            npd = max(min(n+d,nz),1)
            phi(nd) = lev%phi(l,m,npd)
            dx(nd)  = abs(lev%z(n)-lev%z(npd))
        end if 

        if (nd .gt. 0) then 
            ! Some phi values found, proceed

            ! Get order of phi values from least to greatest 
            if (nd .eq. 1) then 
                idx(1) = 1 
            else if (nd .eq. 2) then 
                idx(nd) = maxloc(phi(1:nd),1)
                idx(1)  = minloc(phi(1:nd),1)
            else ! nd .eq. 3
                idx(nd) = maxloc(phi(1:nd),1)
                idx(1)  = minloc(phi(1:nd),1)
                if ( (idx(1) .eq. 1 .and. idx(nd) .eq. 3) .or. &
                     (idx(1) .eq. 3 .and. idx(nd) .eq. 1) ) then 
                    idx(2)  = 2
                else if ( (idx(1) .eq. 1 .and. idx(nd) .eq. 2) .or. &
                          (idx(1) .eq. 2 .and. idx(nd) .eq. 1) ) then 
                    idx(2)  = 3 
                else
                    idx(2)  = 1
                end if 
            end if 

            ! Sort phi and dx 
            phi = phi(idx)
            dx  = dx(idx) 

            do id = 1, nd

                a = real(id,wp)
                b = sum(phi(1:id))
                c = sum(phi(1:id)**2 - dx(1:id)**2)

                if (b**2 - a*c .ge. 0.0) then 
                    phi_temp = b + sqrt(b**2-a*c) / a
                    if (id .lt. nd .and. phi_temp .gt. phi(id+1)) then 
                        ! Do nothing
                    else
                        ! Exit loop, smallest phi found
                        exit
                    end if
                end if 

            end do 

        else
            ! Set output phi value to a really high number
            
            phi_temp = lev%missing_value

        end if 

        return 

    end subroutine levelset_solve_quadratic


    subroutine find_accepted_cells(check,z_srf,z_base,z)
        ! Update check array, setting boundary
        ! cells to 'accepted'

        implicit none

        integer, intent(INOUT) :: check(:,:,:) 
        real(wp), intent(IN)   :: z_srf(:,:) 
        real(wp), intent(IN)   :: z_base(:,:) 
        real(wp), intent(IN)   :: z(:) 

        ! Local variables 
        integer :: i, j, k, n  
        integer :: nx, ny, nz
        integer :: im1, ip1, jm1, jp1, km1, kp1

        nx = size(check,1)
        ny = size(check,2) 
        nz = size(check,3) 

        do i = 1, nx 
        do j = 1, ny

            if (z_base(i,j) .ne. z_srf(i,j)) then 

                do k = 1, nz

                    kp1 = min(k+1,nz)

                    if (z(k)   .le. z_base(i,j) .and. &
                        z(kp1) .gt. z_base(i,j) ) then 

                        check(i,j,k) = 0  

                    else if (z(k)   .le. z_srf(i,j) .and. &
                             z(kp1) .gt. z_srf(i,j) ) then 

                        check(i,j,k) = 0 

                    end if 

                end do 

            end if 

        end do
        end do

        return 

    end subroutine find_accepted_cells

    subroutine find_close_cells(check)
        ! Update check array, setting neighbors
        ! of 'accepted' points to 'close'

        implicit none

        integer, intent(INOUT) :: check(:,:,:) 

        ! Local variables 
        integer :: i, j, k, n  
        integer :: nx, ny, nz
        integer :: im1, ip1, jm1, jp1, km1, kp1

        nx = size(check,1)
        ny = size(check,2) 
        nz = size(check,3) 

        do i = 1, nx 
        do j = 1, ny 
        do k = 1, nz

            im1 = max(i-1,1)
            ip1 = min(i+1,nx)
            jm1 = max(j-1,1)
            jp1 = min(j+1,ny)
            km1 = max(k-1,1)
            kp1 = min(k+1,nz)

            if (check(i,j,k) .ne. 0) then 

                ! Count how many direct neighbors are 'accepted'
                n =   count(check(im1:ip1,j,k) .eq. 0) &
                    + count(check(i,jm1:jp1,k) .eq. 0) &
                    + count(check(i,j,km1:kp1) .eq. 0)

                if (n .gt. 0) then 

                    check(i,j,k) = 1

                end if 

            end if 

        end do 
        end do
        end do

        return 

    end subroutine find_close_cells

    elemental subroutine calc_dist_squared(d,x0,x1)

        implicit none

        real(wp), intent(OUT) :: d
        real(wp), intent(IN)  :: x0
        real(wp), intent(IN)  :: x1

        d = (x1-x0)**2

        return

    end subroutine calc_dist_squared

    elemental subroutine calc_z_srf_max(z_srf,H_ice,z_bed,z_sl,rho_ice,rho_sw)
        ! Calculate surface elevation
        ! Adapted from Pattyn (2017), Eq. 1
        
        implicit none 

        real(prec), intent(INOUT) :: z_srf 
        real(prec), intent(IN)    :: H_ice
        real(prec), intent(IN)    :: z_bed
        real(prec), intent(IN)    :: z_sl
        real(prec), intent(IN)    :: rho_ice
        real(prec), intent(IN)    :: rho_sw

        ! Local variables
        integer :: i, j, nx, ny 
        real(prec) :: rho_ice_sw
        real(prec) :: H_eff

        rho_ice_sw = rho_ice/rho_sw ! Ratio of density of ice to seawater [--]
        
        ! Get effective ice thickness (for now assume grid ice thickness is correct)
        H_eff = H_ice

        ! Initially calculate surface elevation everywhere 
        z_srf = max(z_bed + H_eff, z_sl + (1.0-rho_ice_sw)*H_eff)
        
        return 

    end subroutine calc_z_srf_max

    subroutine calc_halfar(H_ice,x,y,time,R0,H0,lambda,n,A,rho_ice,g)
        ! Equivalent to bueler_test_BC in ice_benchmarks.f90 

        implicit none 

        real(wp), intent(OUT) :: H_ice(:,:) 
        ! real(wp), intent(OUT) :: mbal(:,:) 
        ! real(wp), intent(OUT) :: u_b(:,:) 
        real(wp), intent(IN)  :: x(:)         ! [m] 
        real(wp), intent(IN)  :: y(:)         ! [m]
        real(wp), intent(IN)  :: time         ! [a] Time relative to t0 
        real(wp), intent(IN)  :: R0 
        real(wp), intent(IN)  :: H0 
        real(wp), intent(IN)  :: lambda  
        real(wp), intent(IN)  :: n
        real(wp), intent(IN)  :: A            ! [Pa3 a m-1] 
        real(wp), intent(IN)  :: rho_ice 
        real(wp), intent(IN)  :: g   
         
        ! Local variables 
        integer    :: i, j, k, nx, ny 
        real(wp) :: r_now  
        real(wp) :: R0_meters
        real(wp) :: alpha, beta, gamma, t0, time1  
        real(wp) :: fac 
        
        real(wp), parameter :: f = 0.0        ! isostasy fraction 

        nx = size(H_ice,1)
        ny = size(H_ice,2) 

        ! Get parameter R0 in meters [km] => [m]
        R0_meters = R0 * 1e3 

        ! Calculate alpha, beta, t0 and absolute time
        alpha = (2.0 - (n+1.0)*lambda)/(5.0*n+3.0)
        beta  = (1.0 + (2.0*n+1.0)*lambda) / (5.0*n+3.0)

        gamma = halfar_gamma(A,n,rho_ice,g)
        t0    = (beta/gamma) * ((2.0*n+1.0)/(n+1.0))**n * (R0_meters**(n+1)/H0**(2.0*n+1.0))
        time1 = time + t0 

        H_ice = 0.0_prec 

        do j = 1, ny
        do i = 1, nx 

            ! Calculate the radius value as a function of xx and yy [m]
            r_now = sqrt(x(i)**2 + y(j)**2)

            ! Consider r==x for now (no y dependence)
            !r_now = x(i)

            ! Calculate the Halfar similarity solution profile (Eq. 10-11 in Bueler et al, 2005)
            fac = max(0.0, 1.0 - (((time1/t0)**(-beta))*r_now/R0_meters)**((n+1.0)/n) )
            H_ice(i,j) = H0 * (time1/t0)**(-alpha) * fac**(n/(2.0*n+1.0))

            ! Now calculate implied mass balance
            !mbal(i,j)  = (lambda/time1)*H_ice(i,j)  

        end do 
        end do 

        ! Set the basal velocity to zero everywhere
        !u_b = 0.0 

        return 

    end subroutine calc_halfar

    elemental function halfar_gamma(A,n,rho_ice,g) result(gamma)
        ! Default gamma = 9.0177e-13 m-3 s-1 

        implicit none 

        real(wp), intent(IN) :: A 
        real(wp), intent(IN) :: n 
        real(wp), intent(IN) :: rho_ice
        real(wp), intent(IN) :: g 
        real(wp) :: gamma 

        gamma = 2.0_prec * A * (rho_ice*g)**n / (n+2.0_prec)

        return 

    end function halfar_gamma

    subroutine calc_vel_sia_2D(u,w,x,z,H_ice,z_srf,z_bed,ATT,rho_ice,g,n)

        implicit none 

        real(wp), intent(OUT) :: u(:,:)
        real(wp), intent(OUT) :: w(:,:)
        real(wp), intent(IN)  :: x(:)
        real(wp), intent(IN)  :: z(:)
        real(wp), intent(IN)  :: H_ice(:)
        real(wp), intent(IN)  :: z_srf(:)
        real(wp), intent(IN)  :: z_bed(:)
        real(wp), intent(IN)  :: ATT(:,:) 
        real(wp), intent(IN)  :: rho_ice 
        real(wp), intent(IN)  :: g 
        real(wp), intent(IN)  :: n 

        ! Local variables 
        integer :: i, j, k, nx, ny, nz 
        real(wp) :: np1, np2, nm1, fac 
        real(wp) :: r_now 
        real(wp) :: dzsdx_now 
        real(wp) :: dzsdx_mag_now
        real(wp) :: d2zsdx2_now 
        real(wp) :: dHdx_now
        real(wp) :: H_now 
        real(wp) :: zs_now 
        real(wp) :: zb_now 
        real(wp) :: ATT_now 

        real(wp), parameter :: tol = 1e-9_wp 

        nx = size(x,1)
        ny = 0 
        nz = size(z,1) 

        np1 = n + 1.0_wp 
        np2 = n + 2.0_wp 
        nm1 = n - 1.0_wp 

        fac = -2.0_wp * (rho_ice*g)**n / (n+1.0_wp)

        do i = 1, nx-1 

            ! Get staggered values x(i+1/2)
            ! r_now  = 0.5_wp*(x(i)+x(i+1))
            ! H_now  = 0.5_wp*(H_ice(i)+H_ice(i+1))
            ! zs_now = 0.5_wp*(z_srf(i)+z_srf(i+1))
            ! zb_now = 0.5_wp*(z_bed(i)+z_bed(i+1))
            
            ! dHdx_now  = (H_ice(i+1)-H_ice(i))/(x(i+1)-x(i))
            ! dzsdx_now = (z_srf(i+1)-z_srf(i))/(x(i+1)-x(i))
            ! dzsdx_mag_now = abs(dzsdx_now)

            ! d2zsdx2_now = ( (H_ice(i+1)-H_now)/(x(i+1)-r_now) &
            !                 - (H_now-H_ice(i))/(r_now-x(i)) ) &
            !             / ( 0.5_wp*(r_now+x(i+1)) - 0.5_wp*(x(i)+r_now) ) 

            ! Assume no staggering!
            r_now  = x(i)
            H_now  = H_ice(i)
            zs_now = z_srf(i)
            zb_now = z_bed(i)
            
            if (i .gt. 1) then 
                dHdx_now  = (H_ice(i+1)-H_ice(i-1))/(x(i+1)-x(i-1))
                dzsdx_now = (z_srf(i+1)-z_srf(i-1))/(x(i+1)-x(i-1))
                dzsdx_mag_now = abs(dzsdx_now)

                d2zsdx2_now = ( (H_ice(i+1)-H_ice(i))/(x(i+1)-x(i)) &
                                - (H_ice(i)-H_ice(i-1))/(x(i)-x(i-1)) ) &
                            / ( 0.5_wp*(x(i)+x(i+1)) - 0.5_wp*(x(i-1)+x(i)) ) 
            else 
                dHdx_now  = 0.0_wp 
                dzsdx_now = 0.0_wp 
                dzsdx_mag_now = abs(dzsdx_now)
                
                d2zsdx2_now = 0.0_wp 
            end if 

        do k = 2, nz 
        
            ATT_now = 0.5_wp*(ATT(i,k) + ATT(i+1,k))

            if ( zb_now + z(k) .le. zs_now .and. H_now .gt. 0.0_wp) then 
                u(i,k) = fac*ATT_now &
                        * (H_now**np1 - (zs_now-z(k))**np1) &
                        * dzsdx_mag_now**nm1 * dzsdx_now

                w(i,k) = fac*ATT_now &
                        * ( ((1.0_wp/ (r_now+tol))*dzsdx_now**n + n*dzsdx_now**nm1*d2zsdx2_now) &
                            * ((1.0_wp/np2)*(H_now**np2-(zs_now-z(k))**np2) &
                                 - H_now**np1*(z(k)-zb_now))  &
                            + dzsdx_now**np1 * (H_now**np1 - (zs_now-z(k))**np1) &
                            - np1*dHdx_now*dzsdx_now**n*H_now**n*(zs_now-zb_now) &
                          )
            else
                ! Outside of ice sheet

                u(i,k) = 0.0_wp 
                w(i,k) = 0.0_wp 

            end if 

        end do
        end do 

        ! Set first grid point velocity to zero 
        u(1,:) = 0.0_wp  

        ! Set last grid point velocity to zero 
        u(i,:) = 0.0_wp 
        w(i,:) = 0.0_wp 

        ! Set basal values to zero too (no sliding)
        u(:,1) = 0.0_wp 
        w(:,1) = 0.0_wp 

        return

    end subroutine calc_vel_sia_2D


    subroutine levelset_write_init(lev,filename,time_init,units)

        implicit none 

        type(levelset_class), intent(IN) :: lev 
        character(len=*),  intent(IN) :: filename, units 
        real(wp),          intent(IN) :: time_init
        
        ! Create the empty netcdf file
        call nc_create(filename)

        ! Add grid axis variables to netcdf file
        call nc_write_dim(filename,"x",x=lev%x*1e-3,units="km")
        call nc_write_dim(filename,"y",x=lev%y*1e-3,units="km")
        call nc_write_dim(filename,"z",x=lev%z,     units="meters")
        call nc_write_dim(filename,"time",x=time_init,dx=1.0_prec,nx=1,units=trim(units),unlimited=.TRUE.)

        return

    end subroutine levelset_write_init

    subroutine levelset_write_step(lev,filename,time)

        implicit none 
        
        type(levelset_class), intent(IN) :: lev        
        character(len=*),  intent(IN) :: filename
        real(prec),        intent(IN) :: time

        ! Local variables
        integer    :: ncid, n, nx
        real(prec) :: time_prev

        nx = size(lev%x,1) 

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_size(filename,"time",ncid)
        call nc_read(filename,"time",time_prev,start=[n],count=[1],ncid=ncid) 
        if (abs(time-time_prev).gt.1e-5) n = n+1 

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        ! Write variables
        call nc_write(filename,"H_ice",lev%H_ice(:,1),units="m",long_name="Ice thickness", &
                      dim1="x",dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"z_srf",lev%z_srf(:,1),units="m",long_name="Surface elevation", &
                      dim1="x",dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"z_base",lev%z_base(:,1),units="m",long_name="Ice base elevation", &
                      dim1="x",dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"z_bed",lev%z_bed(:,1),units="m",long_name="Bedrock elevation", &
                      dim1="x",dim2="time",start=[1,n],ncid=ncid)

        call nc_write(filename,"u",lev%u(:,1,:),units="m",long_name="Velocity, x", &
                      dim1="x",dim2="z",dim3="time",start=[1,1,n],ncid=ncid)
        ! call nc_write(filename,"v",lev%v(:,1,:),units="m",long_name="Velocity, y", &
        !               dim1="x",dim2="z",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"w",lev%w(:,1,:),units="m",long_name="Velocity, z", &
                      dim1="x",dim2="z",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"phi",lev%phi(:,1,:),units="1",long_name="Levelset function", &
                      dim1="x",dim2="z",dim3="time",start=[1,1,n],ncid=ncid,missing_value=lev%missing_value)
        call nc_write(filename,"G",lev%G(:,1,:),units="1",long_name="Levelset status (0: accepted, 1: close, 2: far)", &
                      dim1="x",dim2="z",dim3="time",start=[1,1,n],ncid=ncid)

        ! Close the netcdf file
        call nc_close(ncid)

        return

    end subroutine levelset_write_step

end program test_levelset