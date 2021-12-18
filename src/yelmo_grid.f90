module yelmo_grid
    ! This module holds routines to manage the gridded coordinates of a given Yelmo domain 

    use yelmo_defs, only : sp, dp, prec, pi, ygrid_class
    use nml 
    use ncio 

    implicit none 

    interface yelmo_init_grid
        module procedure yelmo_init_grid_fromfile
!         module procedure yelmo_init_grid_fromnml    ! Ambiguous interface with yelmo_init_grid_fromfile
        module procedure yelmo_init_grid_fromname
        module procedure yelmo_init_grid_fromgrd
        module procedure yelmo_init_grid_fromopt
    end interface 
    
    private 
    public :: yelmo_init_grid 
    public :: yelmo_init_grid_fromfile
    public :: yelmo_init_grid_fromnml
    public :: yelmo_init_grid_fromname
    public :: yelmo_init_grid_fromopt 
    public :: yelmo_init_grid_fromgrd
    public :: yelmo_grid_write

contains 

    subroutine yelmo_init_grid_fromfile(grd,filename,grid_name)

        implicit none 

        type(ygrid_class), intent(INOUT) :: grd 
        character(len=*),  intent(IN)    :: filename  
        character(len=*),  intent(IN)    :: grid_name

        ! Local variables
        character(len=56) :: units 

        ! Deallocate grid object to start 
        call ygrid_dealloc(grd)

        ! Define the ygrid name
        grd%name = trim(grid_name)
        
        ! Determine grid axis sizes
        grd%nx = nc_size(filename,"xc")
        grd%ny = nc_size(filename,"yc")
        
        ! Total number of points 
        grd%npts = grd%nx * grd%ny 
        
        ! Allocate axes 
        allocate(grd%xc(grd%nx))
        allocate(grd%yc(grd%ny))
        
        ! Load axes from file
        call nc_read(filename,"xc",grd%xc)
        call nc_read(filename,"yc",grd%yc)

        ! Allocate grid arrays
        allocate(grd%x(grd%nx,grd%ny))
        allocate(grd%y(grd%nx,grd%ny))
        allocate(grd%lon(grd%nx,grd%ny))
        allocate(grd%lat(grd%nx,grd%ny))
        allocate(grd%area(grd%nx,grd%ny))
        
        ! Load variables 
        call nc_read(filename,"x2D",  grd%x)
        call nc_read(filename,"y2D",  grd%y)
        call nc_read(filename,"lon2D",grd%lon)
        call nc_read(filename,"lat2D",grd%lat)
        call nc_read(filename,"area", grd%area)
        
        ! Determine axis units 
        call nc_read_attr(filename,"xc","units",units)

        ! Modify axis values as needed to get units of [m]
        if (trim(units) .eq. "kilometers") then 

            grd%xc   = grd%xc *1e3
            grd%yc   = grd%yc *1e3
            grd%x    = grd%x  *1e3
            grd%y    = grd%y  *1e3

            grd%area = grd%area *1e3*1e3 

        else if (trim(units) .eq. "meters") then 
            ! Pass - do nothing

        else 
            write(*,*) "yelmo_init_grid_fromfile:: &
            &Error grid axis units not recognized: "//trim(units)
            stop 
        end if

        ! Determine grid resolution [m]
        grd%dx = grd%xc(2) - grd%xc(1) 
        grd%dy = grd%yc(2) - grd%yc(1)  
        
        ! Determine whether this is a projected grid 
        if (minval(grd%lon) .ne. maxval(grd%lon)) then 
            grd%is_projection = .TRUE. 
        else 
            grd%is_projection = .FALSE. 
        end if 

        ! Assign default values to projection parameters, then 
        ! determine if actual projection parameters can be loaded
        grd%mtype = "cartesian"
        grd%lambda = 0.0 
        grd%phi    = 0.0 
        grd%alpha  = 0.0 
        grd%scale  = 1.0 
        grd%x_e    = 0.0 
        grd%y_n    = 0.0

        if (grd%is_projection) then 
            ! Load additional projection information if available

            if (.FALSE. .and. nc_exists_attr(filename,"lon2D","grid_mapping") ) then

                ! Read grid map name (projection name, eg "polar_stereographic")
                call nc_read_attr(filename,"lon2D","grid_mapping",grd%mtype)

                if ( nc_exists_var(filename,trim(grd%mtype)) ) then
                    ! If projection variable exists, load attributes 
                    call nc_read_attr(filename,trim(grd%mtype),"straight_vertical_longitude_from_pole",grd%lambda)
                    call nc_read_attr(filename,trim(grd%mtype),"latitude_of_projection_origin",        grd%phi)
                    call nc_read_attr(filename,trim(grd%mtype),"angle_of_oblique_tangent",             grd%alpha)
                    call nc_read_attr(filename,trim(grd%mtype),"scale_factor_at_projection_origin",    grd%scale)
                    call nc_read_attr(filename,trim(grd%mtype),"false_easting",                        grd%x_e)
                    call nc_read_attr(filename,trim(grd%mtype),"false_northing",                       grd%y_n) 
                end if                    

            end if 

        end if 
        
        ! Write grid summary 
        write(*,*) "== Yelmo grid summary =="
        write(*,*) "grid_mapping: ",                            trim(grd%mtype)
        write(*,*) "straight_vertical_longitude_from_pole: ",   grd%lambda
        write(*,*) "latitude_of_projection_origin: ",           grd%phi
        write(*,*) "angle_of_oblique_tangent: ",                grd%alpha
        write(*,*) "scale_factor_at_projection_origin: ",       grd%scale
        write(*,*) "false_easting: ",                           grd%x_e
        write(*,*) "false_northing: ",                          grd%y_n

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
        grd%mtype = "cartesian"
        grd%lambda = 0.0 
        grd%phi    = 0.0 
        grd%alpha  = 0.0 
        grd%scale  = 1.0 
        grd%x_e    = 0.0 
        grd%y_n    = 0.0

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
        units = "kilometers" 

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
        grd%mtype = "cartesian"
        grd%lambda = 0.0 
        grd%phi    = 0.0 
        grd%alpha  = 0.0 
        grd%scale  = 1.0 
        grd%x_e    = 0.0 
        grd%y_n    = 0.0

        return 

    end subroutine yelmo_init_grid_fromname

    subroutine yelmo_init_grid_fromgrd(grd,grd_src,dx,dy)

        implicit none 

        type(ygrid_class), intent(INOUT) :: grd 
        type(ygrid_class), intent(IN)    :: grd_src
        real(prec), intent(IN), optional :: dx 
        real(prec), intent(IN), optional :: dy 
        
        ! Local variables 
        integer :: nx, ny 
        real(prec) :: x0, y0

        ! Assign nx/ny from source grid
        nx = grd_src%nx 
        ny = grd_src%ny 

        write(*,*) "yelmo_init_grid_fromgrd:: Error: not yet implemented - see yelmo_grid.f90."
        stop 

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
        grd%mtype = "cartesian"
        grd%lambda = 0.0 
        grd%phi    = 0.0 
        grd%alpha  = 0.0 
        grd%scale  = 1.0 
        grd%x_e    = 0.0 
        grd%y_n    = 0.0

        return 

    end subroutine yelmo_init_grid_fromaxes 
    
    subroutine yelmo_init_grid_fromopt(grd,grid_name,units,x0,dx,nx,y0,dy,ny,lon,lat,area)

        implicit none 

        type(ygrid_class), intent(INOUT) :: grd 
        character(len=*),  intent(IN)    :: grid_name 
        character(len=*),  intent(IN)    :: units

        real(prec),        intent(IN), optional :: x0
        real(prec),        intent(IN) :: dx
        integer,           intent(IN) :: nx
        real(prec),        intent(IN), optional :: y0
        real(prec),        intent(IN) :: dy
        integer,           intent(IN) :: ny 
        real(prec),        intent(IN), optional :: lon(:,:) 
        real(prec),        intent(IN), optional :: lat(:,:) 
        real(prec),        intent(IN), optional :: area(:,:) 
        
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

                write(*,*) "yelmo_init_grid:: Error: units of input grid parameters must be 'kilometers' or 'meters'."
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
        grd%mtype = "cartesian"
        grd%lambda = 0.0 
        grd%phi    = 0.0 
        grd%alpha  = 0.0 
        grd%scale  = 1.0 
        grd%x_e    = 0.0 
        grd%y_n    = 0.0

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

    subroutine yelmo_grid_write(grid,fnm,create)
        ! Write grid info to netcdf file respecting coordinate conventions
        ! for easier interpretation/plotting etc. 

        implicit none 
        type(ygrid_class), intent(IN) :: grid 
        character(len=*),  intent(IN) :: fnm
        logical,           intent(IN) :: create  

        ! Local variables 
        character(len=16) :: xnm 
        character(len=16) :: ynm 
        
        xnm = "xc"
        ynm = "yc" 

        ! Create the netcdf file if desired
        if (create) then 
            call nc_create(fnm)
        
            ! Add grid axis variables to netcdf file
            call nc_write_dim(fnm,xnm,x=grid%xc*1e-3,units="kilometers")
            call nc_write_attr(fnm,xnm,"_CoordinateAxisType","GeoX")

            call nc_write_dim(fnm,ynm,x=grid%yc*1e-3,units="kilometers")
            call nc_write_attr(fnm,ynm,"_CoordinateAxisType","GeoY")
            
        end if 

        ! Add projection information if needed
        if (grid%is_projection) then
            call nc_write_map(fnm,grid%mtype,dble(grid%lambda),phi=dble(grid%phi), &
                              alpha=dble(grid%alpha),x_e=dble(grid%x_e),y_n=dble(grid%y_n))
        end if 

        call nc_write(fnm,"x2D",grid%x*1e-3,dim1=xnm,dim2=ynm,units="kilometers",grid_mapping=grid%mtype)
        call nc_write(fnm,"y2D",grid%y*1e-3,dim1=xnm,dim2=ynm,units="kilometers",grid_mapping=grid%mtype)

        if (grid%is_projection) then 
            call nc_write(fnm,"lon2D",grid%lon,dim1=xnm,dim2=ynm,grid_mapping=grid%mtype)
            call nc_write_attr(fnm,"lon2D","_CoordinateAxisType","Lon")
            call nc_write(fnm,"lat2D",grid%lat,dim1=xnm,dim2=ynm,grid_mapping=grid%mtype)
            call nc_write_attr(fnm,"lat2D","_CoordinateAxisType","Lat")
        end if 

        call nc_write(fnm,"area",  grid%area*1e-6,  dim1=xnm,dim2=ynm,grid_mapping=grid%mtype,units="km^2")
        if (grid%is_projection) call nc_write_attr(fnm,"area","coordinates","lat2D lon2D")
        !call nc_write(fnm,"border",grid%border,dim1=xnm,dim2=ynm,grid_mapping=grid%mtype)
        !if (grid%is_projection) call nc_write_attr(fnm,"border","coordinates","lat2D lon2D")

        return

    end subroutine yelmo_grid_write


! BELOW: PREVIOUS ROUTINE USED FOR INITIALIZING YELMO GRID WHEN
! COORDINATES LIBRARY WAS USED - NOW REPLACED WITH INTERNAL FUNCTIONS

!     subroutine yelmo_init_grid(grid,grid_name,grid_in)

!         implicit none 

!         type(grid_class), intent(OUT)   :: grid  
!         character(len=*), intent(INOUT) :: grid_name      ! Overwritable if grid_in is present
!         type(grid_class), intent(IN), optional :: grid_in 

!         if (present(grid_in)) then 

!             grid = grid_in 

!             ! Ensure parameter grid_name is consistent with defined grid 
!             grid_name = grid%name 
        
!         else 
!             ! Define yelmo grid from predefined options 

!             select case(trim(grid_name))

!                 ! Note - all North projections now use the ESPG-3413
!                 ! polar stereographic projection with (lambda=-45.d0,phi=70.d0)
!                 ! Smaller Northern domains like Eurasia and Greenland use
!                 ! the same projection for consistency. 
!                 ! ESPG-3413 (lambda=-45.d0,phi=70.d0) is used for Greenland in 
!                 ! model intercomparison exercises, eg ISMIP6. 

!                 ! NORTH DOMAINS ======================= 

!                 case("NH-40KM")
!                     call grid_init(grid,name="NH-40KM",mtype="polar_stereographic",units="kilometers", &
!                             lon180=.TRUE.,x0=-4900.d0,dx=40.0d0,nx=221,y0=-5400.d0,dy=40.0d0,ny=221, &
!                             lambda=-45.d0,phi=70.d0)
                
!                 case("NH-20KM")
!                     call grid_init(grid,name="NH-20KM",mtype="polar_stereographic",units="kilometers", &
!                             lon180=.TRUE.,x0=-4900.d0,dx=20.0d0,nx=441,y0=-5400.d0,dy=20.0d0,ny=441, &
!                             lambda=-45.d0,phi=70.d0)
                
!                 case("NH-10KM")
!                     call grid_init(grid,name="NH-10KM",mtype="polar_stereographic",units="kilometers", &
!                             lon180=.TRUE.,x0=-4900.d0,dx=10.0d0,nx=881,y0=-5400.d0,dy=10.0d0,ny=881, &
!                             lambda=-45.d0,phi=70.d0)
                
!                 case("NH-5KM")
!                     call grid_init(grid,name="NH-5KM",mtype="polar_stereographic",units="kilometers", &
!                             lon180=.TRUE.,x0=-4900.d0,dx=5.0d0,nx=1761,y0=-5400.d0,dy=5.0d0,ny=1761, &
!                             lambda=-45.d0,phi=70.d0)
            
!                 ! EURASIA DOMAINS ======================= 

!                 case("EIS-40KM")
!                     call grid_init(grid,name="EIS-40KM",mtype="polar_stereographic",units="kilometers", &
!                             lon180=.TRUE.,x0=380.d0,dx=40.0d0,nx=89,y0=-5000.d0,dy=40.0d0,ny=161, &
!                             lambda=-45.d0,phi=70.d0)
                
!                 case("EIS-20KM")
!                     call grid_init(grid,name="EIS-20KM",mtype="polar_stereographic",units="kilometers", &
!                             lon180=.TRUE.,x0=380.d0,dx=20.0d0,nx=177,y0=-5000.d0,dy=20.0d0,ny=321, &
!                             lambda=-45.d0,phi=70.d0)
                
!                 case("EIS-10KM")
!                     call grid_init(grid,name="EIS-10KM",mtype="polar_stereographic",units="kilometers", &
!                             lon180=.TRUE.,x0=380.d0,dx=10.0d0,nx=353,y0=-5000.d0,dy=10.0d0,ny=641, &
!                             lambda=-45.d0,phi=70.d0)
                
!                 case("EIS-5KM")
!                     call grid_init(grid,name="EIS-5KM",mtype="polar_stereographic",units="kilometers", &
!                             lon180=.TRUE.,x0=380.d0,dx=5.0d0,nx=705,y0=-5000.d0,dy=5.0d0,ny=1281, &
!                             lambda=-45.d0,phi=70.d0)
                    
!                 ! GREENLAND DOMAINS =======================

!                 case("GRL-40KM")
!                     call grid_init(grid,name="GRL-40KM",mtype="polar_stereographic",units="kilometers", &
!                             lon180=.TRUE.,x0=-720.d0,dx=40.0d0,nx=43,y0=-3450.d0,dy=40.0d0,ny=73, &
!                             lambda=-45.d0,phi=70.d0)
                
!                 case("GRL-20KM")
!                     call grid_init(grid,name="GRL-20KM",mtype="polar_stereographic",units="kilometers", &
!                             lon180=.TRUE.,x0=-720.d0,dx=20.0d0,nx=85,y0=-3450.d0,dy=20.0d0,ny=145, &
!                             lambda=-45.d0,phi=70.d0)
                
!                 case("GRL-10KM")
!                     call grid_init(grid,name="GRL-10KM",mtype="polar_stereographic",units="kilometers", &
!                             lon180=.TRUE.,x0=-720.d0,dx=10.0d0,nx=169,y0=-3450.d0,dy=10.0d0,ny=289, &
!                             lambda=-45.d0,phi=70.d0)
                
!                 case("GRL-5KM")
!                     call grid_init(grid,name="GRL-5KM",mtype="polar_stereographic",units="kilometers", &
!                             lon180=.TRUE.,x0=-720.d0,dx=5.0d0,nx=337,y0=-3450.d0,dy=5.0d0,ny=577, &
!                             lambda=-45.d0,phi=70.d0)
                
!                 case("GRL-2KM")
!                     call grid_init(grid,name="GRL-2KM",mtype="polar_stereographic",units="kilometers", &
!                             lon180=.TRUE.,x0=-720.d0,dx=2.0d0,nx=841,y0=-3450.d0,dy=2.0d0,ny=1441, &
!                             lambda=-45.d0,phi=70.d0)
                
!                 case("GRL-1KM")
!                     call grid_init(grid,name="GRL-1KM",mtype="polar_stereographic",units="kilometers", &
!                             lon180=.TRUE.,x0=-720.d0,dx=1.0d0,nx=1681,y0=-3450.d0,dy=1.0d0,ny=2881, &
!                             lambda=-45.d0,phi=70.d0)

!                 case("Bamber01-20KM")
!                     call grid_init(grid,name="Bamber01-20KM",mtype="polar_stereographic",units="kilometers", &
!                             lon180=.TRUE.,x0=-800.d0,dx=20.d0,nx=76,y0=-3400.d0,dy=20.d0,ny=141, &
!                             lambda=-39.d0,phi=90.d0)

!                 case("Bamber01-10KM")
!                     call grid_init(grid,name="Bamber01-10KM",mtype="polar_stereographic",units="kilometers", &
!                             lon180=.TRUE.,x0=-800.d0,dx=10.d0,nx=151,y0=-3400.d0,dy=10.d0,ny=281, &
!                             lambda=-39.d0,phi=90.d0)

!                 ! ANTARCTICA DOMAINS ======================= 

!                 case("ANT-80KM")
!                     call grid_init(grid,name="ANT-80KM",mtype="polar_stereographic",units="kilometers", &
!                            lon180=.TRUE.,dx=80.d0,nx=79,dy=80.d0,ny=74,lambda=0.d0,phi=-71.d0)

!                 case("ANT-40KM")
!                     call grid_init(grid,name="ANT-40KM",mtype="polar_stereographic",units="kilometers", &
!                            lon180=.TRUE.,dx=40.d0,nx=157,dy=40.d0,ny=147,lambda=0.d0,phi=-71.d0)

!                 case("ANT-20KM")
!                     call grid_init(grid,name="ANT-20KM",mtype="polar_stereographic",units="kilometers", &
!                            lon180=.TRUE.,dx=20.d0,nx=313,dy=20.d0,ny=293,lambda=0.d0,phi=-71.d0)

!                 case("ANT-10KM")
!                     call grid_init(grid,name="ANT-10KM",mtype="polar_stereographic",units="kilometers", &
!                            lon180=.TRUE.,dx=10.d0,nx=625,dy=10.d0,ny=585,lambda=0.d0,phi=-71.d0)

!                 case("ANT-5KM")
!                     call grid_init(grid,name="ANT-5KM",mtype="polar_stereographic",units="kilometers", &
!                            lon180=.TRUE.,dx=5.d0,nx=1249,dy=5.d0,ny=1169,lambda=0.d0,phi=-71.d0)

!                 case("ANT-1KM")
!                     call grid_init(grid,name="ANT-1KM",mtype="polar_stereographic",units="kilometers", &
!                            lon180=.TRUE.,dx=1.d0,nx=6241,dy=1.d0,ny=5841,lambda=0.d0,phi=-71.d0)

!                 case DEFAULT
!                     write(*,*) "yelmo_init_grid:: error: grid name not recognized: "//trim(grid_name)
!                     stop 

!             end select

!         end if 

!         return 

!     end subroutine yelmo_init_grid
    
end module yelmo_grid 

