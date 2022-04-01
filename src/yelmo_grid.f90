module yelmo_grid
    ! This module holds routines to manage the gridded coordinates of a given Yelmo domain 

    use yelmo_defs, only : sp, dp, prec, wp, pi, ygrid_class
    use nml 
    use ncio 

    implicit none 

    interface yelmo_init_grid
        module procedure yelmo_init_grid_fromfile
!         module procedure yelmo_init_grid_fromnml    ! Ambiguous interface with yelmo_init_grid_fromfile
        module procedure yelmo_init_grid_fromname
        module procedure yelmo_init_grid_fromgrd
        module procedure yelmo_init_grid_fromaxes
        module procedure yelmo_init_grid_fromopt
    end interface 
    
    private 

    public :: calc_zeta 

    public :: yelmo_init_grid 
    public :: yelmo_init_grid_fromfile
    public :: yelmo_init_grid_fromnml
    public :: yelmo_init_grid_fromname
    public :: yelmo_init_grid_fromaxes
    public :: yelmo_init_grid_fromopt 
    public :: yelmo_init_grid_fromgrd
    public :: yelmo_grid_write

contains 
    
    subroutine calc_zeta(zeta_aa,zeta_ac,nz_ac,nz_aa,zeta_scale,zeta_exp)
        ! Calculate the vertical axis cell-centers first (aa-nodes),
        ! including a cell centered at the base and the surface. 
        ! Then calculate the vertical axis cell-edges (ac-nodes).
        ! The base and surface ac-nodes coincide with the cell centers.
        ! There is one more cell-edge than cell-centers (nz_ac=nz_aa+1)

        implicit none 

        real(prec), allocatable, intent(INOUT) :: zeta_aa(:) 
        real(prec), allocatable, intent(INOUT) :: zeta_ac(:) 
        integer,                 intent(OUT)   :: nz_ac 
        integer,      intent(IN)   :: nz_aa 
        character(*), intent(IN)   :: zeta_scale 
        real(prec),   intent(IN)   :: zeta_exp 

        ! Local variables
        integer :: k 
        real(prec), allocatable :: tmp(:) 

        ! Define size of zeta ac-node vector
        nz_ac = nz_aa + 1 

        ! First allocate arrays 
        if (allocated(zeta_aa)) deallocate(zeta_aa)
        if (allocated(zeta_ac)) deallocate(zeta_ac)
        allocate(zeta_aa(nz_aa))
        allocate(zeta_ac(nz_ac))

        ! Initially define a linear zeta scale 
        ! Base = 0.0, Surface = 1.0 
        do k = 1, nz_aa
            zeta_aa(k) = 0.0 + 1.0*(k-1)/real(nz_aa-1)
        end do 

        ! Scale zeta to produce different resolution through column if desired
        ! zeta_scale = ["linear","exp","wave"]
        select case(trim(zeta_scale))
            
            case("exp")
                ! Increase resolution at the base 
                
                zeta_aa = zeta_aa**(zeta_exp) 

            case("exp-inv")
                ! Increase resolution at the surface 
                
                zeta_aa = 1.0_wp - zeta_aa**(zeta_exp)

                ! Reverse order 
                allocate(tmp(nz_aa))
                tmp = zeta_aa 
                do k = 1, nz_aa
                    zeta_aa(k) = tmp(nz_aa-k+1)
                end do 

            case("tanh")
                ! Increase resolution at base and surface 

                zeta_aa = tanh(1.0*pi*(zeta_aa-0.5))
                zeta_aa = zeta_aa - minval(zeta_aa)
                zeta_aa = zeta_aa / maxval(zeta_aa)

            case DEFAULT
            ! Do nothing, scale should be linear as defined above
        
        end select  
        
        ! Get zeta_ac (between zeta_aa values, as well as at the base and surface)
        zeta_ac(1) = 0.0 
        do k = 2, nz_ac-1
            zeta_ac(k) = 0.5 * (zeta_aa(k-1)+zeta_aa(k))
        end do 
        zeta_ac(nz_ac) = 1.0 

        ! =================
        ! write(*,*) "Vertical axis:"
        ! do k = nz_ac, 1, -1
        !     if (k .eq. nz_ac) then 
        !         write(*,*) k, -9999.0, zeta_ac(k)
        !     else 
        !         write(*,*) k, zeta_aa(k), zeta_ac(k) 
        !     end if  
        ! end do 
        ! stop 
        ! =================

        return 

    end subroutine calc_zeta
    
    subroutine yelmo_init_grid_fromfile(grd,filename,grid_name,xnm,ynm,lonnm,latnm)

        implicit none 

        type(ygrid_class), intent(INOUT) :: grd 
        character(len=*),  intent(IN)    :: filename  
        character(len=*),  intent(IN)    :: grid_name
        character(len=*),  intent(IN), optional :: xnm 
        character(len=*),  intent(IN), optional :: ynm 
        character(len=*),  intent(IN), optional :: lonnm 
        character(len=*),  intent(IN), optional :: latnm 
        
        ! Local variables
        integer :: i, j 
        character(len=56) :: units 
        character(len=56) :: x_name
        character(len=56) :: y_name
        character(len=56) :: lon_name
        character(len=56) :: lat_name
        character(len=56) :: crs_name 

        ! Get arguments 
        x_name = "xc"
        if (present(xnm)) x_name = trim(xnm)

        y_name = "yc"
        if (present(ynm)) y_name = trim(ynm)
        
        lon_name = "lon2D"
        if (present(lonnm)) lon_name = trim(lonnm)
        
        lat_name = "lat2D"
        if (present(latnm)) lat_name = trim(latnm)
        

        ! Deallocate grid object to start 
        call ygrid_dealloc(grd)

        ! Define the ygrid name
        grd%name = trim(grid_name)

        ! Determine grid axis sizes
        grd%nx = nc_size(filename,x_name)
        grd%ny = nc_size(filename,y_name)
        
        ! Total number of points 
        grd%npts = grd%nx * grd%ny 
        
        ! Allocate axes 
        allocate(grd%xc(grd%nx))
        allocate(grd%yc(grd%ny))
        
        ! Allocate grid arrays
        allocate(grd%x(grd%nx,grd%ny))
        allocate(grd%y(grd%nx,grd%ny))
        allocate(grd%lon(grd%nx,grd%ny))
        allocate(grd%lat(grd%nx,grd%ny))
        allocate(grd%area(grd%nx,grd%ny))
        
        ! Determine axis units 
        call nc_read_attr(filename,x_name,"units",units)

        ! Load axes from file
        call nc_read(filename,x_name,grd%xc)
        call nc_read(filename,y_name,grd%yc)

        ! Modify axis values as needed to get units of [m]
        select case(trim(units))

            case("kilometers","km") 

                grd%xc   = grd%xc *1e3
                grd%yc   = grd%yc *1e3

            case("meters","m")

                ! Pass - do nothing

            case DEFAULT 
                write(*,*) "yelmo_init_grid_fromfile:: &
                &Error grid axis units not recognized: "//trim(units)
                stop 
        
        end select 

        ! Determine grid resolution [m]
        grd%dx = grd%xc(2) - grd%xc(1) 
        grd%dy = grd%yc(2) - grd%yc(1)  
        
        ! Populate x and y 2D arrays from axis values 
        do j = 1, grd%ny 
            grd%x(:,j) = grd%xc 
        end do 
        do i = 1, grd%nx 
            grd%y(i,:) = grd%yc 
        end do 

        ! Determine cell area
        if (nc_exists_var(filename,"area")) then 
            ! Load cell area from file
            call nc_read(filename,"area", grd%area)

            select case(trim(units))
            
                case("kilometers","km") 

                    grd%area = grd%area *1e3*1e3

            end select

        else
            ! Calculate cell area directly [m^2]
            grd%area = grd%dx * grd%dy 
        end if 


        ! Load latitude and longitude values
        call nc_read(filename,lon_name,grd%lon)
        call nc_read(filename,lat_name,grd%lat)

        
        ! Determine whether this is a projected grid 
        if (minval(grd%lon) .ne. maxval(grd%lon)) then 
            grd%is_projection = .TRUE. 
        else 
            grd%is_projection = .FALSE. 
        end if 

        ! Assign default values to projection parameters, then 
        ! determine if actual projection parameters can be loaded
        grd%mtype               = "cartesian"
        grd%lambda              = 0.0 
        grd%phi                 = 0.0 
        grd%alpha               = 0.0 
        grd%scale               = 1.0 
        grd%x_e                 = 0.0 
        grd%y_n                 = 0.0
        grd%semi_major_axis     = 0.0 
        grd%inverse_flattening  = 0.0
        grd%is_sphere           = .TRUE. 
        
        if (grd%is_projection) then 
            ! Load additional projection information if available

            if ( nc_exists_attr(filename,lon_name,"grid_mapping") ) then

                ! Read grid map name (projection name, eg "polar_stereographic")
                call nc_read_attr(filename,lon_name,"grid_mapping",crs_name)

                if ( nc_exists_var(filename,crs_name) .and. trim(crs_name) .eq. "crs" ) then
                    ! If projection variable exists, load attributes 
                    ! ajr: checking crs_name=="crs" ensures that projection info
                    ! is only loaded from newer files that have the right projection
                    ! attributes defined. Older files have crs_name=="polar_stereographic"
                    ! for example, and do not contain parameters following cf-conventions. 

                    ! Load projection name (polar_stereographic,stereographic,etc)
                    call nc_read_attr(filename,crs_name,"grid_mapping_name",grd%mtype)
                    
                    select case(trim(grd%mtype))
                    
                        case("polar_stereographic")
                            call nc_read_attr(filename,crs_name,"straight_vertical_longitude_from_pole",grd%lambda)
                            !call nc_read_attr(filename,crs_name,"latitude_of_projection_origin",grd%phi_proj_orig)
                            !   ajr: latitude_of_projection_origin not needed as it is -90 or 90. 
                            call nc_read_attr(filename,crs_name,"standard_parallel",grd%phi)
                            call nc_read_attr(filename,crs_name,"false_easting",grd%x_e)
                            call nc_read_attr(filename,crs_name,"false_northing",grd%y_n) 
                            call nc_read_attr(filename,crs_name,"semi_major_axis",grd%semi_major_axis)
                            call nc_read_attr(filename,crs_name,"inverse_flattening",grd%inverse_flattening)

                            ! Define alpha too - just in case
                            grd%alpha = 90.0_wp - abs(grd%phi)

                        case("stereographic")
                            call nc_read_attr(filename,crs_name,"longitude_of_projection_origin",grd%lambda)
                            call nc_read_attr(filename,crs_name,"latitude_of_projection_origin",grd%phi)
                            call nc_read_attr(filename,crs_name,"angle_of_oblique_tangent",grd%alpha)
                            call nc_read_attr(filename,crs_name,"scale_factor_at_projection_origin",grd%scale)
                            call nc_read_attr(filename,crs_name,"false_easting",grd%x_e)
                            call nc_read_attr(filename,crs_name,"false_northing",grd%y_n) 
                            call nc_read_attr(filename,crs_name,"semi_major_axis",grd%semi_major_axis)
                            call nc_read_attr(filename,crs_name,"inverse_flattening",grd%inverse_flattening)

                        case DEFAULT 

                            write(*,*) "yelmo_init_grid_fromfile:: Warning: unsupported grid projection. &
                                    &Projection parameters not loaded."
                            write(*,*) "map type: ", trim(grd%mtype)
                    
                    end select

                    ! Check whether working on a sphere 
                    if (grd%inverse_flattening .ne. 0.0) then 
                        grd%is_sphere = .FALSE. 
                    end if 

                end if                    

            end if 

        end if  

        ! Write grid summary 
        write(*,*) "== Yelmo grid summary =="
        write(*,*) "grid_mapping:                           ",  trim(grd%mtype)
        write(*,*) "straight_vertical_longitude_from_pole:  ",  grd%lambda
        write(*,*) "standard_parallel:                      ",  grd%phi
        write(*,*) "angle_of_oblique_tangent:               ",  grd%alpha
        write(*,*) "scale_factor_at_projection_origin:      ",  grd%scale
        write(*,*) "false_easting:                          ",  grd%x_e
        write(*,*) "false_northing:                         ",  grd%y_n
        write(*,*) "semi_major_axis:                        ",  grd%semi_major_axis
        write(*,*) "inverse_flattening:                     ",  grd%inverse_flattening
        write(*,*) 
        write(*,"(a16,i8)")        "Total points = ",   grd%npts
        write(*,"(a16,2g15.5)")    "range(xc) = ",      minval(grd%xc),maxval(grd%xc)
        write(*,"(a16,2g15.5)")    "range(yc) = ",      minval(grd%yc),maxval(grd%yc)
        write(*,"(a16,2g15.5)")    "range(lon) = ",     minval(grd%lon),maxval(grd%lon)
        write(*,"(a16,2g15.5)")    "range(lat) = ",     minval(grd%lat),maxval(grd%lat)
        write(*,"(a16,2g15.5)")    "range(area) = ",    minval(grd%area),maxval(grd%area)

        return 

    end subroutine yelmo_init_grid_fromfile

    subroutine yelmo_init_grid_fromnml(grd,filename,grid_name)

        implicit none 

        type(ygrid_class), intent(INOUT) :: grd 
        character(len=*),  intent(IN)    :: filename 
        character(len=*),  intent(IN)    :: grid_name 
        
        ! Local variables 
        character(len=56) :: units 
        integer :: nx, ny 
        real(prec) :: x0, y0, dx, dy 

        ! Use the grid_name to determine grid parameters 

        call nml_read(filename,trim(grid_name),"units",units)
        call nml_read(filename,trim(grid_name),"x0",x0)
        call nml_read(filename,trim(grid_name),"y0",y0)
        call nml_read(filename,trim(grid_name),"dx",dx)
        call nml_read(filename,trim(grid_name),"dy",dy)
        call nml_read(filename,trim(grid_name),"nx",nx)
        call nml_read(filename,trim(grid_name),"ny",ny)

        ! Now call grid define routine with the desired grid parameters 
        if (x0 .ne. -9999.0) then 
            ! Call with specific x0/y0 locations
            call yelmo_init_grid_fromopt(grd,grid_name,units,x0=x0,dx=dx,nx=nx,y0=y0,dy=dy,ny=ny)
        else
            ! Call to generate grid centered on 0/0
            call yelmo_init_grid_fromopt(grd,grid_name,units,dx=dx,nx=nx,dy=dy,ny=ny)
        end if 

        ! Assign default values to projection parameters
        ! (no projection information available with this method so far)
        grd%mtype               = "cartesian"
        grd%lambda              = 0.0 
        grd%phi                 = 0.0 
        grd%alpha               = 0.0 
        grd%scale               = 1.0 
        grd%x_e                 = 0.0 
        grd%y_n                 = 0.0
        grd%semi_major_axis     = 0.0 
        grd%inverse_flattening  = 0.0
        grd%is_sphere           = .TRUE. 
 

        return 

    end subroutine yelmo_init_grid_fromnml

    subroutine yelmo_init_grid_fromname(grd,grid_name)

        implicit none 

        type(ygrid_class), intent(INOUT) :: grd  
        character(len=*),  intent(IN)    :: grid_name  
        
        ! Local variables 
        character(len=56) :: units 
        integer :: nx, ny 
        real(prec) :: x0, y0, dx, dy 

        ! Set units of grid parameters defined below
        units = "km" 

        select case(trim(grid_name))

            ! Note - all North projections now use the ESPG-3413
            ! polar stereographic projection with (lambda=-45.d0,phi=70.d0)
            ! Smaller Northern domains like Eurasia and Greenland use
            ! the same projection for consistency. 
            ! ESPG-3413 (lambda=-45.d0,phi=70.d0) is used for Greenland in 
            ! model intercomparison exercises, eg ISMIP6. 

            ! NORTH DOMAINS ======================= 

            case("NH-40KM")
                x0 = -4900.0
                y0 =  5400.0 
                dx =    40.0 
                dy =    40.0 
                nx =   221 
                ny =   221

            case("NH-20KM")
                x0 = -4900.0
                y0 =  5400.0 
                dx =    20.0 
                dy =    20.0 
                nx =   441 
                ny =   441
                
            case("NH-10KM")
                x0 = -4900.0
                y0 =  5400.0 
                dx =    10.0 
                dy =    10.0 
                nx =   881 
                ny =   881
            
            case("NH-5KM")
                x0 = -4900.0
                y0 =  5400.0 
                dx =     5.0 
                dy =     5.0 
                nx =  1761 
                ny =  1761
            
            ! EURASIA DOMAINS ======================= 

            case("EIS-40KM")
                x0 =   380.0
                y0 = -5000.0 
                dx =    40.0 
                dy =    40.0 
                nx =    89 
                ny =   161
            
            case("EIS-20KM")
                x0 =   380.0
                y0 = -5000.0 
                dx =    20.0 
                dy =    20.0 
                nx =   177 
                ny =   321
            
            case("EIS-10KM")
                x0 =   380.0
                y0 = -5000.0 
                dx =    10.0 
                dy =    10.0 
                nx =   353 
                ny =   641
            
            case("EIS-5KM")
                x0 =   380.0
                y0 = -5000.0 
                dx =     5.0 
                dy =     5.0 
                nx =   705 
                ny =  1281
            
            ! GREENLAND DOMAINS =======================

            case("GRL-40KM")
                x0 =  -720.0
                y0 = -3450.0 
                dx =    40.0 
                dy =    40.0 
                nx =    43 
                ny =    73
            
            case("GRL-20KM")
                x0 =  -720.0
                y0 = -3450.0 
                dx =    20.0 
                dy =    20.0 
                nx =    85 
                ny =   145
            
            case("GRL-10KM")
                x0 =  -720.0
                y0 = -3450.0 
                dx =    10.0 
                dy =    10.0 
                nx =   169 
                ny =   289
            
            case("GRL-5KM")
                x0 =  -720.0
                y0 = -3450.0 
                dx =     5.0 
                dy =     5.0 
                nx =   337 
                ny =   577
            
            case("GRL-2KM")
                x0 =  -720.0
                y0 = -3450.0 
                dx =     2.0 
                dy =     2.0 
                nx =   841 
                ny =  1441
            
            case("GRL-1KM")
                x0 =  -720.0
                y0 = -3450.0 
                dx =     1.0 
                dy =     1.0 
                nx =  1681 
                ny =  2881
            

            case("Bamber01-20KM")
                x0 =  -800.0
                y0 = -3400.0  
                dx =    20.0 
                dy =    20.0 
                nx =    76 
                ny =   141
            
            case("Bamber01-10KM")
                x0 =  -800.0
                y0 = -3400.0 
                dx =    10.0 
                dy =    10.0 
                nx =   151 
                ny =   281

            ! ANTARCTICA DOMAINS ======================= 

            case("ANT-80KM")
                x0 = -9999.0
                y0 = -9999.0 
                dx =    80.0 
                dy =    80.0 
                nx =    79 
                ny =    74
            
            case("ANT-40KM")
                x0 = -9999.0
                y0 = -9999.0 
                dx =    40.0 
                dy =    40.0 
                nx =   157 
                ny =   147
            
            case("ANT-20KM")
                x0 = -9999.0
                y0 = -9999.0 
                dx =    20.0 
                dy =    20.0 
                nx =   313 
                ny =   293
            
            case("ANT-10KM")
                x0 = -9999.0
                y0 = -9999.0 
                dx =    10.0 
                dy =    10.0 
                nx =   625 
                ny =   585
                   
            case("ANT-5KM")
                x0 = -9999.0
                y0 = -9999.0 
                dx =     5.0 
                dy =     5.0 
                nx =  1249 
                ny =  1169
            
            case("ANT-1KM")
                x0 = -9999.0
                y0 = -9999.0 
                dx =     1.0 
                dy =     1.0 
                nx =  6241 
                ny =  5841
                
            case DEFAULT
                write(*,*) "yelmo_init_grid_fromname:: error: grid name not recognized: "//trim(grid_name)
                stop 

        end select 

        ! Now call grid define routine with the desired grid parameters 
        if (x0 .ne. -9999.0) then 
            ! Call with specific x0/y0 locations
            call yelmo_init_grid_fromopt(grd,grid_name,units,x0=x0,dx=dx,nx=nx,y0=y0,dy=dy,ny=ny)
        else
            ! Call teo generate grid centered on 0/0
            call yelmo_init_grid_fromopt(grd,grid_name,units,dx=dx,nx=nx,dy=dy,ny=ny)
        end if 

        ! Assign default values to projection parameters
        ! (no projection information available with this method so far)
        grd%mtype               = "cartesian"
        grd%lambda              = 0.0 
        grd%phi                 = 0.0 
        grd%alpha               = 0.0 
        grd%scale               = 1.0 
        grd%x_e                 = 0.0 
        grd%y_n                 = 0.0
        grd%semi_major_axis     = 0.0 
        grd%inverse_flattening  = 0.0
        grd%is_sphere           = .TRUE. 
 

        return 

    end subroutine yelmo_init_grid_fromname

    subroutine yelmo_init_grid_fromgrd(grd,grid_name,grd_src,dx,dy)
        ! Generate a new grid (grd) with the same bounds
        ! as source grid (grd_src), but with new resolution 

        implicit none 

        type(ygrid_class), intent(INOUT) :: grd 
        character(len=*),  intent(IN)    :: grid_name 
        type(ygrid_class), intent(IN)    :: grd_src
        real(wp), intent(IN)             :: dx 
        real(wp), intent(IN)             :: dy 
        
        ! Local variables 
        integer  :: nx, ny 
        real(wp) :: x0, y0

        ! First ensure all variables are deallocated 
        call ygrid_dealloc(grd)

        ! Set grid name 
        grd%name = trim(grid_name)

        ! Determine lower-left corner (starting point for new grid)
        x0 = grd_src%xc(1)
        y0 = grd_src%yc(1) 

        ! Assign nx/ny from source grid
        nx = grd_src%nx 
        ny = grd_src%ny 

        write(*,*) "yelmo_init_grid_fromgrd:: Error: not yet implemented - see yelmo_grid.f90."
        stop 

        ! axis_init(x,x0,dx)



!         ! Or, determine nx, ny based on desired dx/dy
!         if (present(dx)) 
        
        return 

    end subroutine yelmo_init_grid_fromgrd

    subroutine yelmo_init_grid_fromaxes(grd,grid_name,xc,yc,lon,lat,area)

        implicit none 

        type(ygrid_class), intent(INOUT) :: grd 
        character(len=*),  intent(IN)    :: grid_name 
        real(prec),        intent(IN)    :: xc(:) 
        real(prec),        intent(IN)    :: yc(:) 
        real(prec),        intent(IN), optional :: lon(:,:) 
        real(prec),        intent(IN), optional :: lat(:,:) 
        real(prec),        intent(IN), optional :: area(:,:) 
        
        ! Local variables
        integer    :: i, j

        ! First ensure all variables are deallocated 
        call ygrid_dealloc(grd)

        ! Set grid name 
        grd%name = trim(grid_name)

        ! Set up the grid axis info based on axes
        
        ! x-axis
        grd%nx = size(xc) 
        allocate(grd%xc(grd%nx))
        grd%xc  = xc 

        ! y-axis
        grd%ny = size(yc) 
        allocate(grd%yc(grd%ny))
        grd%yc  = yc

        ! Total number of points 
        grd%npts = grd%nx * grd%ny 
        
        ! Determine grid resolution [m]
        grd%dx = grd%xc(2) - grd%xc(1) 
        grd%dy = grd%yc(2) - grd%yc(1)  
        
        ! Grid variables 

        allocate(grd%x(grd%nx,grd%ny))
        allocate(grd%y(grd%nx,grd%ny))
        allocate(grd%lon(grd%nx,grd%ny))
        allocate(grd%lat(grd%nx,grd%ny))
        allocate(grd%area(grd%nx,grd%ny))
        
        grd%x    = 0.0 
        grd%y    = 0.0 
        grd%lon  = 0.0 
        grd%lat  = 0.0 
        grd%area = 0.0 
        
        ! x array
        do j = 1, grd%ny 
            grd%x(:,j) = grd%xc 
        end do 

        ! y array
        do i = 1, grd%nx 
            grd%y(i,:) = grd%yc 
        end do 

        ! Grid cell area [m^2]
        if (present(area)) then 
            grd%area = area
        else 
            grd%area = grd%dx*grd%dy 
        end if 
        
        ! Optional helper variables 
        if (present(lon))  grd%lon  = lon 
        if (present(lat))  grd%lat  = lat 
        
        ! Assign default values to projection parameters
        ! (no projection information available with this method so far)
        grd%mtype               = "cartesian"
        grd%lambda              = 0.0 
        grd%phi                 = 0.0 
        grd%alpha               = 0.0 
        grd%scale               = 1.0 
        grd%x_e                 = 0.0 
        grd%y_n                 = 0.0
        grd%semi_major_axis     = 0.0 
        grd%inverse_flattening  = 0.0
        grd%is_sphere           = .TRUE. 
 

        return 

    end subroutine yelmo_init_grid_fromaxes
    
    subroutine yelmo_init_grid_fromopt(grd,grid_name,units,x0,dx,nx,y0,dy,ny,lon,lat,area)

        implicit none 

        type(ygrid_class), intent(INOUT) :: grd 
        character(len=*),  intent(IN)    :: grid_name 
        character(len=*),  intent(IN)    :: units

        real(wp),          intent(IN), optional :: x0
        real(wp),          intent(IN) :: dx
        integer,           intent(IN) :: nx
        real(wp),          intent(IN), optional :: y0
        real(wp),          intent(IN) :: dy
        integer,           intent(IN) :: ny 
        real(wp),          intent(IN), optional :: lon(:,:) 
        real(wp),          intent(IN), optional :: lat(:,:) 
        real(wp),          intent(IN), optional :: area(:,:) 
        
        ! Local variables
        integer    :: i, j

        ! First ensure all variables are deallocated 
        call ygrid_dealloc(grd)

        ! Set grid name 
        grd%name = trim(grid_name)

        ! Set up the grid axis info based on options present 
        
        ! x-axis 
        grd%nx = nx 
        allocate(grd%xc(grd%nx))
        call axis_init(grd%xc,x0=x0,dx=dx)

        ! y-axis
        grd%ny = ny
        allocate(grd%yc(grd%ny))
        call axis_init(grd%yc,x0=y0,dx=dy)

        ! Total number of points 
        grd%npts = grd%nx * grd%ny 

        ! Assign dx, dy
        grd%dx = dx  
        grd%dy = dy 
        
        ! Modify axis values depending on input units 
        select case(trim(units))

            case("kilometers","km")

                grd%xc = grd%xc*1e3 
                grd%yc = grd%yc*1e3 
                grd%dx = grd%dx*1e3
                grd%dy = grd%dy*1e3 

            case("meters","m")

                ! Pass - do nothing 

            case DEFAULT 

                write(*,*) "yelmo_init_grid:: Error: units of input grid parameters &
                            &must be one of: 'kilometers', 'km', 'meters', 'm'."
                write(*,*) "units: ", trim(units)
                stop 

        end select 

        ! Grid variables 

        allocate(grd%x(grd%nx,grd%ny))
        allocate(grd%y(grd%nx,grd%ny))
        allocate(grd%lon(grd%nx,grd%ny))
        allocate(grd%lat(grd%nx,grd%ny))
        allocate(grd%area(grd%nx,grd%ny))
        
        grd%x    = 0.0 
        grd%y    = 0.0 
        grd%lon  = 0.0 
        grd%lat  = 0.0 
        grd%area = 0.0 
        
        ! x array
        do j = 1, grd%ny 
            grd%x(:,j) = grd%xc 
        end do 

        ! y array
        do i = 1, grd%nx 
            grd%y(i,:) = grd%yc 
        end do 

        ! Grid cell area [m^2]
        if (present(area)) then 
            grd%area = area
        else 
            grd%area = grd%dx*grd%dy 
        end if 
        
        ! Optional helper variables 
        if (present(lon))  grd%lon  = lon 
        if (present(lat))  grd%lat  = lat 
        
        ! Assign default values to projection parameters
        ! (no projection information available with this method so far)
        grd%mtype               = "cartesian"
        grd%lambda              = 0.0 
        grd%phi                 = 0.0 
        grd%alpha               = 0.0 
        grd%scale               = 1.0 
        grd%x_e                 = 0.0 
        grd%y_n                 = 0.0
        grd%semi_major_axis     = 0.0 
        grd%inverse_flattening  = 0.0
        grd%is_sphere           = .TRUE. 
 

        return 

    end subroutine yelmo_init_grid_fromopt

    subroutine axis_init(x,x0,dx)

        implicit none 

        real(prec) :: x(:)
        real(prec), optional :: x0, dx
        real(prec) :: dx_tmp 
        integer :: i, nx  

        nx = size(x) 

        do i = 1, nx 
            x(i) = real(i-1,prec)
        end do 

        dx_tmp = 1.d0 
        if (present(dx)) dx_tmp = dx 
        
        x = x*dx_tmp  

        if (present(x0)) then 
            x = x + x0 
        else
            x = x + (-(nx-1.0)/2.0*dx_tmp)
        end if 

        return 
    end subroutine axis_init

    subroutine ygrid_dealloc(grd)

        implicit none 

        type(ygrid_class), intent(INOUT) :: grd 
        
        ! Axis vectors
        if (allocated(grd%xc))   deallocate(grd%xc)
        if (allocated(grd%yc))   deallocate(grd%yc)

        ! Grid arrays
        if (allocated(grd%x))    deallocate(grd%x)
        if (allocated(grd%y))    deallocate(grd%y)
        if (allocated(grd%lon))  deallocate(grd%lon)
        if (allocated(grd%lat))  deallocate(grd%lat)
        if (allocated(grd%area)) deallocate(grd%area)
        
        return 

    end subroutine ygrid_dealloc

    subroutine yelmo_grid_write(grid,fnm,domain,grid_name,create)
        ! Write grid info to netcdf file respecting coordinate conventions
        ! for easier interpretation/plotting etc. 

        implicit none 
        type(ygrid_class), intent(IN) :: grid 
        character(len=*),  intent(IN) :: fnm
        character(len=*),  intent(IN) :: domain
        character(len=*),  intent(IN) :: grid_name
        logical,           intent(IN) :: create  

        ! Local variables 
        character(len=16) :: xnm 
        character(len=16) :: ynm 
        character(len=16) :: grid_mapping_name

        xnm = "xc"
        ynm = "yc" 
        grid_mapping_name = "crs" 

        ! Create the netcdf file if desired
        if (create) then 
            
            ! Create the empty netcdf file
            call nc_create(fnm)

            ! Write some general attributes that can be useful (e.g., for interpolation etc)
            call nc_write_attr(fnm, "domain",    trim(domain))
            call nc_write_attr(fnm, "grid_name", trim(grid_name))


            ! Add grid axis variables to netcdf file
            call nc_write_dim(fnm,xnm,x=grid%xc*1e-3,units="km")

            call nc_write_dim(fnm,ynm,x=grid%yc*1e-3,units="km")
            
            if (grid%is_projection) then 
                call nc_write_attr(fnm,xnm,"standard_name","projection_x_coordinate")
                call nc_write_attr(fnm,ynm,"standard_name","projection_y_coordinate")
            end if 

        end if 

        ! Add projection information if needed
        if (grid%is_projection) then
            call nc_write_map(fnm,grid%mtype,dble(grid%lambda),phi=dble(grid%phi), &
                              alpha=dble(grid%alpha),x_e=dble(grid%x_e),y_n=dble(grid%y_n), &
                              is_sphere=grid%is_sphere,semi_major_axis=dble(grid%semi_major_axis),& 
                              inverse_flattening=dble(grid%inverse_flattening))
        end if 

        call nc_write(fnm,"x2D",grid%x*1e-3,dim1=xnm,dim2=ynm,units="km",grid_mapping=grid_mapping_name)
        call nc_write(fnm,"y2D",grid%y*1e-3,dim1=xnm,dim2=ynm,units="km",grid_mapping=grid_mapping_name)

        if (grid%is_projection) then 
            call nc_write(fnm,"lon2D",grid%lon,dim1=xnm,dim2=ynm,grid_mapping=grid_mapping_name)
            call nc_write_attr(fnm,"lon2D","units","degrees_east")
            call nc_write(fnm,"lat2D",grid%lat,dim1=xnm,dim2=ynm,grid_mapping=grid_mapping_name)
            call nc_write_attr(fnm,"lat2D","units","degrees_north")
        end if 

        call nc_write(fnm,"area",  grid%area*1e-6,  dim1=xnm,dim2=ynm,grid_mapping=grid_mapping_name,units="km^2")
        if (grid%is_projection) call nc_write_attr(fnm,"area","coordinates","lat2D lon2D")
        !call nc_write(fnm,"border",grid%border,dim1=xnm,dim2=ynm,grid_mapping=grid_mapping_name)
        !if (grid%is_projection) call nc_write_attr(fnm,"border","coordinates","lat2D lon2D")

        return

    end subroutine yelmo_grid_write

end module yelmo_grid

