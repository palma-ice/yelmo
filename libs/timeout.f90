module timeout

    use nml 
    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit

    implicit none

    integer,  parameter :: sp  = kind(1.0)
    integer,  parameter :: dp  = kind(1.0d0)
    integer,  parameter :: wp  = sp 

    real(wp), parameter :: MV  = -9999.0_wp 

    real(wp), parameter :: time_tol = 1e-6 
    logical,  parameter :: verbose  = .TRUE. 

    type timeout_class
        character(len=56)     :: method
        character(len=56)     :: label 
        real(wp), allocatable :: times(:)
        character(len=56), allocatable :: vnms(:)
    end type

    private
    public :: timeout_check
    public :: timeout_class
    public :: timeout_init 

contains

    function timeout_check(tm,time) result(out_now)
        ! Check if the current time exists within the timeout time vector.

        implicit none

        type(timeout_class), intent(IN) :: tm
        real(wp), intent(IN) :: time 
        logical :: out_now 

        ! Local variables 
        integer :: k, n 
        
        ! Get total times to check 
        n = size(tm%times)

        ! Assume this timestep is not in timeout
        out_now = .FALSE. 

        do k = 1, n 
            if (abs(time - tm%times(k)) .lt. time_tol) then 
                out_now = .TRUE. 
                exit 
            end if
        end do 

        if (verbose .and. out_now) then 
            write(*,"(a16,a,2x,g12.3)") "timeout_check ", trim(tm%label), time 
        end if

        return

    end function timeout_check

    subroutine timeout_init(tm,filename,group,label,time_init,time_end)

        implicit none 

        type(timeout_class), intent(INOUT) :: tm
        character(len=*),    intent(IN)    :: filename 
        character(len=*),    intent(IN)    :: group 
        character(len=*),    intent(IN)    :: label 
        real(wp),            intent(IN)    :: time_init 
        real(wp),            intent(IN)    :: time_end
        
        ! Local variables
        integer  :: k, n, k0, k1 
        real(wp) :: dt_const
        character(len=512) :: timeout_file 

        integer, parameter :: nmax = 100000
        real(wp) :: times(nmax)
        real(wp), allocatable :: tmp(:) 

        ! Store label for logging 
        tm%label = trim(label) 

        ! Load parameters
        call nml_read(filename,group,"method",tm%method)
        
        times = MV 

        select case(trim(tm%method))

            case("const")
                
                call nml_read(filename,group,"dt",dt_const)

                if (dt_const .gt. 0.0) then
                    
                    n = ceiling( (time_end-time_init)/dt_const ) + 1
                    
                    if (n .gt. nmax) then 
                        write(error_unit,*) "timeout_init:: Error: too many output timesteps desired. &
                        &Maximum value limited to nmax = ", nmax 
                        write(error_unit,*) "time_init = ", time_init 
                        write(error_unit,*) "time_end  = ", time_end 
                        write(error_unit,*) "dt        = ", dt_const 
                        write(error_unit,*) "n         = ", n 
                        stop
                    end if

                    do k = 1, n 
                        times(k) = time_init + (k-1)*dt_const 
                    end do

                    k0 = 1 
                    k1 = k0+n-1

                end if

            case("file","times")

                if (trim(tm%method) .eq. "file") then 
                    ! Load time information from an input file 

                    call nml_read(filename,group,"file",timeout_file)

                    ! Get a vector of times 
                    call timeout_parse_file(times,timeout_file)

                else
                    ! Load time information from a parameter
                
                    call nml_read(filename,group,"times",times,init=.FALSE.)

                end if 

                k0 = -1 
                k1 = -1 

                do k = 1, nmax 
                    if (k0 .lt. 0 .and. times(k) .ge. time_init .and. times(k) .ne. MV) then 
                        k0 = k
                    end if 
                    if (k0 .gt. 0 .and. times(k) .le. time_end  .and. times(k) .ne. MV) then
                        k1 = k 
                    end if
                end do

                if (k0 .lt. 0 .or. k1 .lt. 0) then 
                    ! No proper range of indices was found, 
                    ! prescribe first and list times as output times.
                    k0 = 1
                    k1 = 2
                    times = MV 
                    times(1) = time_init
                    times(2) = time_end 
                end if

                n = k1-k0+1

            case DEFAULT

                write(error_unit,*) "timeout_init:: Error: timeout method not recognized."
                write(error_unit,*) "timeout.method = ", trim(tm%method)
                stop 
            
        end select 

        ! Make sure last timestep is also written
        if (times(k1) .lt. time_end) then 
            n  = n+1 
            k1 = k0+n-1
            times(k1) = time_end 
        end if
                    
        ! Store final times(k0:k1) vector of length n
        ! in timeout object for later use. 
        if (allocated(tm%times)) deallocate(tm%times)
        allocate(tm%times(n))
        tm%times(1:n) = times(k0:k1)

        ! Finally remove any duplicate times
        call remove_duplicate_times(tm%times)
        n = size(tm%times)
        
        if (verbose) then 
            ! Print summary
            write(*,*) "timeout: ", trim(tm%label)
            write(*,*) "  method    = ", trim(tm%method)
            write(*,*) "  time_init = ", time_init 
            write(*,*) "  time_end  = ", time_end 
            write(*,*) "  n         = ", n 

            do k = 1, n, 5
                k0 = k 
                k1 = min(k+5-1,n)
                write(*,"(5g12.3)") tm%times(k0:k1)
            end do 

        end if 

        return

    end subroutine timeout_init

    subroutine remove_duplicate_times(times)

        implicit none 

        real(wp), intent(INOUT), allocatable :: times(:) 

        ! Local variables 
        integer :: k, n, n1 
        real(wp) :: time
        real(wp), allocatable :: times_new(:) 

        
        n = size(times)
        allocate(times_new(n))

        times_new = MV 

        n1 = 0 
        do k = 1, n
            time = times(k)
            if (minval(abs(time-times_new)) .lt. time_tol) then 
                ! Do nothing, this value is already contained in times
            else 
                ! This time does not exist yet, add it. 
                n1 = n1+1 
                times_new(n1) = time 
            end if
        end do 

        deallocate(times)
        allocate(times(n1))
        times = times_new(1:n1)

        return

    end subroutine remove_duplicate_times

    subroutine timeout_parse_file(times,filename)

        implicit none 

        real(wp),         intent(INOUT) :: times(:) 
        character(len=*), intent(IN)    :: filename 

        ! Local variables 
        integer :: q, k, ntmp, nmax, io, stat  
        character(len=256) :: line_str 
        real(wp), allocatable :: tmp(:) 

        integer, parameter :: nmax_file = 10000

        open(newunit=io,file=filename,status="old",action="read")

        k = 0 

        do q = 1, nmax_file 

            ! Read line as a string first 
            read(io,"(a)",iostat=stat) line_str 
            line_str = trim(adjustl(line_str))

            if (stat .ne. 0) exit
            if (trim(line_str) .eq. "" .or. line_str(1:1) .eq. "#") then 
                ! Skip this line, do nothing, it is either a comment or an empty line
            else
                ! This line should contain information we want

                if (index(line_str,":") .gt. 0) then 
                    ! Parse this line as a set of times time1:dtime:time2

                    call parse_time_vector(tmp,trim(line_str))

                    ntmp = size(tmp)
                    k = k+1
                    times(k:k+ntmp-1) = tmp 
                    k = k+ntmp-1

                else
                    ! Assume only one value is available

                    k = k+1
                    times(k) = string_to_wp(line_str)
                    !write(*,*) k, " ", times(k)

                end if

            end if

        end do 

        close(io) 

        return 

    end subroutine timeout_parse_file

    subroutine parse_time_vector(times,timestr)
        ! Given a string of format time0:dt:time1,
        ! generate a vector of times: [time0,time0+dt,...,time1]
        implicit none

        real(wp),         intent(INOUT), allocatable :: times(:)
        character(len=*), intent(IN)    :: timestr 

        ! Local variables 
        integer  :: k, k0, k1, n 
        real(wp) :: t0, t1, dt 

        ! Step 1: parse the string into values of t0, dt, and t1

        k0 = index(timestr,":")
        k1 = index(timestr,":",back=.TRUE.)

        if (k0 .eq. 0 .or. k1 .eq. 0) then 
            write(*,*) "parse_time_vector:: Error: timestr does not have the &
            &right format: t0:dt:t1."
            write(*,*) "timestr = ", trim(timestr)
            stop 
        end if 

        read(timestr(1:k0-1),*) t0 
        read(timestr(k0+1:k1-1),*) dt 
        read(timestr(k1+1:len_trim(timestr)),*) t1 

        n = ceiling( (t1-t0)/dt ) + 1

        if (allocated(times)) deallocate(times)
        allocate(times(n))

        do k = 1, n 
            times(k) = t0 + (k-1)*dt 
        end do

        return

    end subroutine parse_time_vector

    ! === Helper functions (borrowed from nml.f90) ===
    
    function string_to_wp(string) result(value)

        implicit none 

        character(len=*), intent(IN) :: string 
        real(wp) :: value 

        character(len=256) :: tmpstr 
        integer  :: stat, n
        real(wp) :: x 

        tmpstr = trim(adjustl(string))
        n      = len_trim(tmpstr)

        read(tmpstr(1:n),*,IOSTAT=stat) x

        value = 0
        if (stat .eq. 0) then 
            value = x 
        else
            n = len_trim(tmpstr)-1
            READ(tmpstr(1:n),*,IOSTAT=stat) x
            if (stat .ne. 0) then 
                write(*,*) "nml:: ","Error converting string to number!"
                write(*,*) "|",trim(tmpstr),"|",n,stat,x
            else
                value = x 
            end if 
        end if 

        return 

    end function string_to_wp

    subroutine string_to_vector(string,value)

        implicit none 

        character(len=*), intent(IN) :: string 
        character(len=*) :: value(:)
        character(len=256) :: tmpvec(size(value))
        character(len=256) :: tmpstr, fmt 
        integer :: stat, n, q, q1, q2, j 

        tmpstr = trim(adjustl(string))
        n      = len_trim(tmpstr)+2

        tmpvec(:) = "" 

        q1 = 1 
        do q = 1, size(tmpvec)
            q2 = index(tmpstr(q1:n)," ") + q1
            if (q2 .gt. q1 .and. q2 .le. n) then 
                tmpvec(q) = tmpstr(q1:q2-1)
                q1 = q2

                ! Make sure gaps of more than one space are properly handled
                do j = 1, 1000
                    if (tmpstr(q1:q1) == " ") q1 = q1+1
                    if (q1 .ge. n) exit 
                end do 

!                 ! Eliminate quotes
!                 q2 = len_trim(tmpvec(q))
!                 if (tmpvec(q)(1:1) == '"') tmpvec(q) = trim(adjustl(tmpvec(q)(2:q2)))
!                 q2 = len_trim(tmpvec(q))
!                 if (tmpvec(q)(q2:q2) == '"') tmpvec(q) = trim(tmpvec(q)(1:q2-1))
                ! Remove quotes around string if they exist 
                call remove_quotes_comma(tmpvec(q))
            
            end if 
        end do 
        
        value = tmpvec 

        return 

    end subroutine string_to_vector

    subroutine remove_quotes_comma(string)

        implicit none 
        character(len=*), intent(INOUT) :: string 
        integer :: i, n 

!         ! Eliminate quotes
!         n = len_trim(string)
!         if (n == 1 .and. trim(string) == '"') then 
!             string = ""
!         else if (n > 0) then 
!             if (string(1:1) == '"') string = trim(adjustl(string(2:n)))
!             n = len_trim(string)
!             if (n > 1  .and. string(n:n) == '"') string = trim(string(1:n-1))
!             if (n == 1 .and. string(n:n) == '"') string = ""
            
!         end if 

        ! Eliminate quotes
        n = len_trim(string)
        do i = 1,n 
            if (string(i:i) == '"' .or. string(i:i) == "'") string(i:i) = " "
        end do 
        string = trim(adjustl(string))

        ! Remove final comma too
        n = len_trim(string)
        if (n > 0) then 
            if (string(n:n) == ",") string(n:n) = " "
            string = trim(adjustl(string))
        end if 
        
        return 

    end subroutine remove_quotes_comma

end module timeout