program check_sim 
    
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
    
    integer :: narg, n
    character(len=256) :: fldr_path, file_path  
    real(wp) :: time, rmse_H, rmse_uxy, rmse_uxy_log 

    ! Get command line arguments 

    narg = command_argument_count()

    if (narg .lt. 1) then 
        write(*,*) "check_sim:: Error: folder must be provided as an argument."
        stop 
    end if 

    call get_command_argument(1,fldr_path)

    file_path = trim(fldr_path)//"/yelmo2D.nc" 

    n = nc_size(file_path,"time")

    call nc_read(file_path,"time",        time,        start=[n],count=[1])
    call nc_read(file_path,"rmse_H",      rmse_H,      start=[n],count=[1])
    call nc_read(file_path,"rmse_uxy",    rmse_uxy,    start=[n],count=[1])
    call nc_read(file_path,"rmse_uxy_log",rmse_uxy_log,start=[n],count=[1])
    
    write(*,*) time, rmse_H, rmse_uxy, rmse_uxy_log 
    
contains 


end program check_sim 