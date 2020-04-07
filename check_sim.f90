program check_sim 
    ! for D in tmp/yelmo991/ant-pd-ens1/* ; do ./check_sim.x $D ; done

    use ncio 

    implicit none 

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the precision of the library (sp,dp)
    integer,  parameter :: wp = sp 

    ! Missing value and aliases
    real(wp), parameter :: MISSING_VALUE_DEFAULT = real(-9999.0,wp)
    real(wp), parameter :: MISSING_VALUE = MISSING_VALUE_DEFAULT
    real(wp), parameter :: MV = MISSING_VALUE_DEFAULT
    
    integer :: narg, nt, nx, ny, n_iso
    character(len=256) :: fldr_path, file_path  
    real(wp) :: time, rmse_H, rmse_uxy, rmse_uxy_log, rmse_H2000
    real(wp), allocatable :: rmse_iso(:) 
    character(len=256) :: fldr_path_root, fldr_sim 

    real(wp), allocatable :: H_ice(:,:) 
    real(wp), allocatable :: H_ice_pd_err(:,:) 
    logical,  allocatable :: mask(:,:) 
    
    ! Get command line arguments 

    narg = command_argument_count()

    if (narg .lt. 1) then 
        write(*,*) "check_sim:: Error: folder must be provided as an argument."
        stop 
    end if 

    call get_command_argument(1,fldr_path)

    ! Determine just the simulation folder for output purposes 
    call split_string(fldr_path,"/",fldr_path_root,fldr_sim,back=.TRUE.)

    file_path = trim(fldr_path)//"/yelmo2D.nc" 

    nt = nc_size(file_path,"time")
    nx = nc_size(file_path,"xc")
    ny = nc_size(file_path,"yc")
    
    call nc_read(file_path,"time",        time,        start=[nt],count=[1])
    call nc_read(file_path,"rmse_H",      rmse_H,      start=[nt],count=[1])
    call nc_read(file_path,"rmse_uxy",    rmse_uxy,    start=[nt],count=[1])
    call nc_read(file_path,"rmse_uxy_log",rmse_uxy_log,start=[nt],count=[1])
    
    ! Allocate variables 
    allocate(H_ice(nx,ny))
    allocate(H_ice_pd_err(nx,ny))
    allocate(mask(nx,ny)) 

    ! Read in H_ice and error at nt 
    call nc_read(file_path,"H_ice",H_ice,start=[1,1,nt],count=[nx,ny,1])
    call nc_read(file_path,"H_ice_pd_err",H_ice_pd_err,start=[1,1,nt],count=[nx,ny,1])

    ! Define mask as points with ice thickness above 2000m 
    mask = H_ice .ge. 2000.0 
    call calc_rmse(rmse_H2000,H_ice_pd_err,mask,missing_value)

    n_iso = nc_size(file_path,"age_iso") 
    allocate(rmse_iso(n_iso))
    call nc_read(file_path,"rmse_iso",rmse_iso,start=[1,nt],count=[n_iso,1])

    write(*,"(a,f10.2,3f8.1,f8.2,50f8.1)") trim(fldr_sim), time, rmse_H, rmse_H2000, rmse_uxy, rmse_uxy_log, rmse_iso 

contains 

    subroutine split_string(instring,delim,string1,string2,back)
        ! split a string into 2 either side of a delimiter token

        character(*), intent(IN)  :: instring
        character(*), intent(IN)  :: delim 
        character(*), intent(OUT) :: string1
        character(*), intent(OUT) :: string2
        logical,      intent(IN), optional :: back 

        ! Local variables
        integer :: index
        character(len=len(instring)) :: string0

        string0 = trim(adjustl(instring))

        index   = scan(string0,delim,back=back)
        if (string0(index+1:index+1) .eq. " ") string0(index:index) = " "  ! Remove trailing character
        index   = scan(string0,delim,back=back)
        string1 = string0(1:index-1)
        string2 = string0(index+1:)

        return 

    end subroutine split_string

    subroutine calc_rmse(rmse,err,mask,missing_value)

        implicit none 

        real(wp), intent(OUT) :: rmse 
        real(wp), intent(IN)  :: err(:,:) 
        logical,  intent(IN)  :: mask(:,:) 
        real(wp), intent(IN)  :: missing_value

        ! Local variables 
        integer :: npts 

        npts = count(mask)

        if (npts .eq. 0) then 
            ! No points found, set missing value 

            rmse = missing_value 

        else 
            ! Calculate masked rmse 

            rmse = sqrt( sum(err**2,mask=mask) / real(npts,wp) )

        end if 

        return 

    end subroutine calc_rmse 

end program check_sim 
