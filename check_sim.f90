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
    
    integer :: narg, n
    character(len=256) :: fldr_path, file_path  
    real(wp) :: time, rmse_H, rmse_uxy, rmse_uxy_log 
    character(len=256) :: fldr_path_root, fldr_sim 

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

    n = nc_size(file_path,"time")

    call nc_read(file_path,"time",        time,        start=[n],count=[1])
    call nc_read(file_path,"rmse_H",      rmse_H,      start=[n],count=[1])
    call nc_read(file_path,"rmse_uxy",    rmse_uxy,    start=[n],count=[1])
    call nc_read(file_path,"rmse_uxy_log",rmse_uxy_log,start=[n],count=[1])
    
    write(*,*) trim(fldr_sim), time, rmse_H, rmse_uxy, rmse_uxy_log 

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

end program check_sim 