program test_levelset

    use yelmo 

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
    real(wp) :: dx
    real(wp) :: xmax, ymin, ymax, zmin, zmax 
    
    type levelset_class 
        real(wp) :: dx, dy, dz
        integer  :: nx, ny, nz  
        real(wp), allocatable :: x(:)
        real(wp), allocatable :: y(:)
        real(wp), allocatable :: z(:)
        

    end type

    type(levelset_class) :: lev1 

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
    file2D     = trim(outfldr)//"yelmo2D.nc"
    file1D     = trim(outfldr)//"yelmo1D.nc"
    
    
    ! General initialization of yelmo constants (used globally)
    call yelmo_global_init(path_const)


    ! Define the domain and grid [km]
    dx   =  10.0 
    xmax = 800.0
    ymax =  50.0
    ymin = -50.0
    call yelmo_init_grid(yelmo1%grd,grid_name,units="km",x0=0.0,dx=dx,nx=int(xmax/dx)+1,y0=ymin,dy=dx,ny=int((ymax-ymin)/dx)+1)

    ! Define horizontal coordinates [km]
    lev1%dx = dx 
    lev1%nx = yelmo1%grd%nx 
    allocate(lev1%x(lev1%nx))
    lev1%x = yelmo1%grd%xc 

    lev1%dy = dx 
    lev1%ny = yelmo1%grd%ny 
    allocate(lev1%y(lev1%ny))
    lev1%y = yelmo1%grd%yc 
    
    ! Define vertical coordinates [m]
    zmin = -3000.0 
    zmax =  7000.0 
    lev1%dz   = 50.0 
    lev1%nz   = int((zmax-zmin)/lev1%dz)+1
    allocate(lev1%z(lev1%nz))
    do k = 1, lev1%nz 
        lev1%z(k) = zmin + (k-1)*lev1%dz 
    end do


    ! Summary of grid 
    write(*,*) 
    write(*,*) "====="
    write(*,*) "x: ", lev1%nx, lev1%dx, minval(lev1%x), maxval(lev1%x)
    write(*,*) "y: ", lev1%ny, lev1%dy, minval(lev1%y), maxval(lev1%y)
    write(*,*) "z: ", lev1%nz, lev1%dz, minval(lev1%z), maxval(lev1%z)
    


contains





end program test_levelset