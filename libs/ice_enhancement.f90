module ice_enhancement 
    ! Module to treat transient changes to the enhancement factor.
    ! The 2D enhancement factor field calculated here is applied
    ! as a boundary condition to the ice sheet model, and is then
    ! propagated as a tracer inside of the ice sheet. 

    use nml 

    implicit none 

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the precision of the library (sp,dp)
    integer,  parameter :: wp = sp 

    type series_type
        character(len=512) :: filename 
        real(wp), allocatable :: time(:), var(:)
    end type 

    type ice_enh_class 
        ! Parameter value 
        real(wp) :: enh_glac 
        real(wp) :: enh_int 

        type(series_type) :: series 
        real(wp) :: time, enh  
    end type 

    private
    public :: ice_enh_class 
    public :: ice_enh_init 
    public :: ice_enh_update 

contains 

    subroutine ice_enh_init(ienh,filename)

        implicit none 

        type(ice_enh_class), intent(INOUT) :: ienh 
        character(len=*),      intent(IN)  :: filename 

        ! Local variables 
        character(len=56) :: varname  
        integer           :: k0, k1, n 
        real(wp)          :: time0, alpha0, time1, alpha1 

        ! Determine filename from which to load sea level time series 
        call nml_read(filename,"ice_enh","enh_path",ienh%series%filename,init=.TRUE.)
        call nml_read(filename,"ice_enh","enh_glac",ienh%enh_glac,init=.TRUE.)
        call nml_read(filename,"ice_enh","enh_int", ienh%enh_int,init=.TRUE.)

        ! Read the time series from ascii file
        call read_series(ienh%series,ienh%series%filename)

        ! Default value of enh should be scaled to one for Holocene (warm periods) and 
        ! maximum value at predefined time during glacial periods (via time series definition)

        ! Define time0 for present day and time1 for LGM 
        time0 = 0.0_wp 
        time1 = -21e3 

        ! First scale time series so that values are zero at time0 
        ! and one at time1 
        k0 = minloc(abs(ienh%series%time - time0),dim=1)
        k1 = minloc(abs(ienh%series%time - time1),dim=1)
        alpha0 = ienh%series%var(k0)
        alpha1 = ienh%series%var(k1) 

        ienh%series%var = (ienh%series%var - alpha0) / (alpha1-alpha0)

        ! Limit to range of 0:1 
        where(ienh%series%var .lt. 0.0) ienh%series%var = 0.0 
        where(ienh%series%var .gt. 1.0) ienh%series%var = 1.0 

        ! Calculate enhancement factor through time as a weighted-mean between
        ! interglacial (enh_int) and glacial (enh_glac) values 
        ienh%series%var = (1.0-ienh%series%var)*ienh%enh_int + ienh%series%var*ienh%enh_glac 

        write(*,*) "ice_enh_init:: "
        write(*,*) " time0: ", time0, ienh%series%var(k0)
        write(*,*) " time1: ", time1, ienh%series%var(k1)
        write(*,*) "range(enh): ", minval(ienh%series%var), maxval(ienh%series%var)
        
        return 

    end subroutine ice_enh_init 

    subroutine ice_enh_update(ienh,time)

        implicit none 

        type(ice_enh_class), intent(INOUT) :: ienh 
        real(wp),            intent(IN)    :: time  ! [yr] 

        ! Local variables 
        real(wp) :: enh_ref 

        ienh%time = time 
        ienh%enh  = series_interp(ienh%series,time)
        
        return 

    end subroutine ice_enh_update 

    subroutine read_series(series,filename)
        ! This subroutine will read a time series of
        ! two columns [time,var] from an ascii file.
        ! Header should be commented by "#" or "!"
        implicit none 

        type(series_type) :: series 
        character(len=*)  :: filename 

        integer, parameter :: f = 190
        integer, parameter :: nmax = 10000

        integer :: i, stat, n 
        character(len=256) :: str, str1 
        real(wp) :: x(nmax), y(nmax) 

        ! Open file for reading 
        open(f,file=filename,status="old")

        ! Read the header in the first line: 
        read(f,*,IOSTAT=stat) str

        do i = 1, nmax 
            read(f,'(a100)',IOSTAT=stat) str 

            ! Exit loop if the end-of-file is reached 
            if(IS_IOSTAT_END(stat)) exit 

            str1 = adjustl(trim(str))
!            str1=str
            if ( len(trim(str1)) .gt. 0 ) then 
                if ( .not. (str1(1:1) == "!" .or. &
                            str1(1:1) == "#") ) then 
                    read(str1,*) x(i), y(i) 
                end if
            end if  
        end do 


        ! Close the file
        close(f) 

        if (i .eq. nmax) then 
            write(*,*) "read_series:: warning: "// &
                       "Maximum length of time series reached, ", nmax
            write(*,*) "Time series in the file may be longer: ", trim(filename)
        end if 

        ! Allocate the time series object and store output data
        n = i-1 
        call series_allocate(series,n)

        series%time = x(1:n) 
        series%var  = y(1:n) 

        write(*,*) "read_series:: Time series read from file: "//trim(filename)
        write(*,*) "    range time: ",minval(series%time), maxval(series%time)
        write(*,*) "    range var : ",minval(series%var),  maxval(series%var)
        
        return 

    end subroutine read_series 

    function series_interp(series,time) result(var)
        ! Wrapper for simple `interp_linear` function for series_types. 
        implicit none 

        type(series_type) :: series 
        real(wp) :: time 
        real(wp) :: var 
        integer :: nt, i 

        ! Interpolate series object
        var = interp_linear(series%time,series%var,xout=time)

        return 

    end function series_interp 

    function interp_linear(x,y,xout) result(yout)
        ! Simple linear interpolation of a point

        implicit none 
 
        real(wp), dimension(:), intent(IN) :: x, y
        real(wp), intent(IN) :: xout
        real(wp) :: yout 
        integer :: i, j, n, nout 
        real(wp) :: alph

        n    = size(x) 

        if (xout .lt. x(1)) then
            yout = y(1)
        else if (xout .gt. x(n)) then
            yout = y(n)
        else
            do j = 1, n 
                if (x(j) .ge. xout) exit 
            end do

            if (j .eq. 1) then 
                yout = y(1) 
            else if (j .eq. n+1) then 
                yout = y(n)
            else 
                alph = (xout - x(j-1)) / (x(j) - x(j-1))
                yout = y(j-1) + alph*(y(j) - y(j-1))
            end if 
        end if 

        return 

    end function interp_linear
    
    subroutine series_allocate(series,nt)

        implicit none 

        type(series_type) :: series 
        integer :: nt 

        if (allocated(series%time))  deallocate(series%time)
        if (allocated(series%var))   deallocate(series%var)

        allocate(series%time(nt))
        allocate(series%var(nt))

        ! Initialize variables to zero
        series%time  = 0.0_wp
        series%var   = 0.0_wp

        return 

    end subroutine series_allocate

end module ice_enhancement  
