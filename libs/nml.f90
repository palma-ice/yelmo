
module nml 

    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit

    implicit none 

    logical :: VERBOSE = .FALSE.            !!  should freshly read namelist be printed to screen?
    logical :: ERROR_NO_PARAM = .TRUE.      !! Should error be thrown if parameter isn't found? 

    integer, parameter :: io_unit_err = error_unit 

    interface nml_read 
        module procedure nml_read_string, nml_read_double, nml_read_float 
        module procedure nml_read_integer, nml_read_logical
        module procedure nml_read_string_vector, nml_read_double_vector
        module procedure nml_read_float_vector, nml_read_integer_vector
        module procedure nml_read_logical_vector
    end interface 

    interface nml_write
        module procedure nml_write_string 
    end interface 

    interface nml_print 
        module procedure nml_print_string, nml_print_double, nml_print_float  
        module procedure nml_print_integer, nml_print_logical
        module procedure nml_print_string_vector, nml_print_double_vector
        module procedure nml_print_float_vector, nml_print_integer_vector
        module procedure nml_print_logical_vector
    end interface 

    private 
    public :: nml_read, nml_write 
    public :: nml_print 
    public :: nml_set_verbose
    public :: nml_replace

contains

    subroutine nml_set_verbose(verb)

        implicit none 

        logical :: verb 

        VERBOSE = verb 

        write(*,*) "nml_set_verbose:: ", VERBOSE 

        return 

    end subroutine nml_set_verbose

    subroutine nml_replace(s,text,rep,outs)
        ! Adapted from FUNCTION Replace_Text:
        ! http://fortranwiki.org/fortran/show/String_Functions
        CHARACTER(len=*)           :: s,text,rep
        CHARACTER(len=*), optional :: outs
        INTEGER :: i, nt, nr

        character(len=LEN(s)+100) :: tmps  ! Temp string to hold output 

        tmps = s ; nt = LEN_TRIM(text) ; nr = LEN_TRIM(rep)
        DO
           i = INDEX(tmps,text(:nt)) ; IF (i == 0) EXIT
           tmps = tmps(:i-1) // rep(:nr) // tmps(i+nt:)
        END DO

        if (present(outs)) then
            outs = trim(tmps)
        else
            s = trim(tmps)
        end if

        return

    end subroutine nml_replace

    ! =============================================================
    !
    ! nml parameter reading functions
    !
    ! =============================================================

    ! This is the basic nml reading subroutine
    ! All interfaces use this to read the parameter, then it
    ! is converted to the correct type
    subroutine nml_read_internal(filename,group,name,value,comment)

        implicit none 

        character(len=*), intent(INOUT) :: value 
        character(len=*), intent(IN)    :: filename, group, name 
        character(len=*), intent(INOUT), optional :: comment   

        ! Local variables
        integer :: io, file_unit 
        integer :: iostat, l, ltype 
        character(len=1000) :: line, name1, value1, comment1 
        logical :: ingroup
        logical :: found_param 

        ! Open the nml filename to be read, or get file units if already open
        inquire(file=trim(filename),NUMBER=file_unit)
        if (file_unit .eq. -1) then
            ! Open a connection to the file 
            open(unit=newunit(io),file=trim(filename),status="old",iostat=iostat)
            if (iostat /= 0) then 
                write(*,*) "nml:: namelist file could not be opened: "//trim(filename)
                stop 
            end if
        else 
            ! File is already open, rewind it
            io = file_unit 
            rewind(io)
        end if 
 
        
        ingroup     = .FALSE. 
        found_param = .FALSE. 

        do l = 1, 5000
            read(io,"(a1000)",iostat=iostat) line 
            if (iostat /= 0) exit 
            call parse_line(line, ltype, name1, value1, comment1 )

            ! Check if the parameter has been found
            if (ingroup .and. ltype == 3) then 
                if (trim(name1) == trim(name)) then 
                    value       = trim(value1)
                    comment     = trim(comment1)
                    found_param = .TRUE. 
                    exit 
                end if 
            end if 

            ! Open and close group as necessary
            if (ltype == 1 .and. trim(name1) == trim(group)) ingroup = .TRUE. 
            if (ltype == 2)                                  ingroup = .FALSE. 

            if (l .eq. 5000) then 
                write(*,*) "nml:: Warning: maximum nml length of 5000 lines reached."
            end if 
        end do 

        if (file_unit .eq. -1) then
            ! Close the file that was opened here
            close(io)
        endif

        if (ERROR_NO_PARAM .and. (.not. found_param)) then 
            write(io_unit_err,*) ""
            write(io_unit_err,*) "nml:: Error: parameter not found."
            write(io_unit_err,*) "Filename:  ", trim(filename)
            write(io_unit_err,*) "Parameter: ", trim(name)
            stop "Program stopped."
        end if 

        return 

    end subroutine nml_read_internal

    subroutine nml_read_string(filename,group,name,value,comment,init)

        implicit none 

        character(len=*), intent(INOUT) :: value 
        character(len=*), intent(IN)    :: filename, group, name 
        character(len=*), intent(INOUT), optional :: comment 
        logical, optional :: init  
        logical :: init_var 
        character(len=256) :: value_str 

        init_var = .FALSE. 
        if (present(init)) init_var = init 
        if (init_var) value = "" 

        ! First find parameter value as a string 
        value_str = ""
        call nml_read_internal(filename,group,name,value_str,comment)

        if (value_str /= "") value = trim(value_str)

        if (VERBOSE) call nml_print(name,value,comment)  ! Check

        return 

    end subroutine nml_read_string

    subroutine nml_read_double(filename,group,name,value,comment,init)

        implicit none 

        double precision, intent(INOUT) :: value 
        character(len=*), intent(IN)    :: filename, group, name 
        character(len=*), intent(INOUT), optional :: comment   
        logical, optional :: init  
        logical :: init_var 
        character(len=256) :: value_str 

        init_var = .FALSE. 
        if (present(init)) init_var = init 
        if (init_var) value = 0.d0 

        ! First find parameter value as a string 
        value_str = ""
        call nml_read_internal(filename,group,name,value_str,comment)
        if (value_str /= "") value = string_to_double(value_str)

        if (VERBOSE) call nml_print(name,value,comment)  ! Check

        return 

    end subroutine nml_read_double 

    subroutine nml_read_float(filename,group,name,value,comment,init)

        implicit none 

        real(4), intent(INOUT) :: value 
        character(len=*), intent(IN)    :: filename, group, name 
        character(len=*), intent(INOUT), optional :: comment  
        logical, optional :: init  
        logical :: init_var  
        character(len=256) :: value_str 

        init_var = .FALSE. 
        if (present(init)) init_var = init 
        if (init_var) value = 0.0 

        ! First find parameter value as a string 
        value_str = ""
        call nml_read_internal(filename,group,name,value_str,comment)
        if (value_str /= "") value = real(string_to_double(value_str))

        if (VERBOSE) call nml_print(name,value,comment)  ! Check

        return 

    end subroutine nml_read_float 

    subroutine nml_read_integer(filename,group,name,value,comment,init)

        implicit none 

        integer, intent(INOUT) :: value 
        character(len=*), intent(IN)    :: filename, group, name 
        character(len=*), intent(INOUT), optional :: comment 
        logical, optional :: init  
        logical :: init_var   
        character(len=256) :: value_str 

        init_var = .FALSE. 
        if (present(init)) init_var = init 
        if (init_var) value = 0 

        ! First find parameter value as a string 
        value_str = ""
        call nml_read_internal(filename,group,name,value_str,comment)

        if (value_str /= "") then 
            value = nint(string_to_double(value_str))
        end if 

        if (VERBOSE) call nml_print(name,value,comment)  ! Check

        return 

    end subroutine nml_read_integer

    subroutine nml_read_logical(filename,group,name,value,comment,init)

        implicit none 

        logical, intent(INOUT) :: value 
        character(len=*), intent(IN)    :: filename, group, name 
        character(len=*), intent(INOUT), optional :: comment 
        logical, optional :: init  
        logical :: init_var   
        character(len=256) :: value_str 

        init_var = .FALSE. 
        if (present(init)) init_var = init 
        if (init_var) value = .FALSE. 

        ! First find parameter value as a string 
        value_str = ""
        call nml_read_internal(filename,group,name,value_str,comment)
        if (value_str /= "") value = string_to_logical(value_str)

        if (VERBOSE) call nml_print(name,value,comment)  ! Check

        return 

    end subroutine nml_read_logical 

    !! Vectors 

    subroutine nml_read_string_vector(filename,group,name,value,comment,init)

        implicit none 

        character(len=*), intent(INOUT) :: value(:) 
        character(len=*), intent(IN)    :: filename, group, name 
        character(len=*), intent(INOUT), optional :: comment 
        logical, optional :: init  
        logical :: init_var  
        character(len=256) :: value_str 

        init_var = .FALSE. 
        if (present(init)) init_var = init 
        if (init_var) value(:) = "" 

        ! First find parameter value as a string 
        value_str = ""
        call nml_read_internal(filename,group,name,value_str,comment)

        if (value_str /= "") call string_to_vector(value_str,value)

        if (VERBOSE) call nml_print(name,value,comment)  ! Check

        return 

    end subroutine nml_read_string_vector 

    subroutine nml_read_double_vector(filename,group,name,value,comment,init)

        implicit none 

        double precision, intent(INOUT) :: value(:) 
        character(len=*), intent(IN)    :: filename, group, name 
        character(len=*), intent(INOUT), optional :: comment   
        character(len=256) :: value_str, value_str_vec(size(value))
        logical, optional :: init  
        logical :: init_var  
        integer :: q 

        init_var = .FALSE. 
        if (present(init)) init_var = init 
        if (init_var) value(:) = 0.d0

        ! First find parameter value as a string 
        value_str = ""
        call nml_read_internal(filename,group,name,value_str,comment)

        if (value_str /= "") then
            call string_to_vector(value_str,value_str_vec)
            do q = 1, size(value)
                if (trim(value_str_vec(q)) /= "") then 
                    value(q) = string_to_double(trim(adjustl(value_str_vec(q))))
                end if 

            end do 
        end if 

        if (VERBOSE) call nml_print(name,value,comment)  ! Check

        return 

    end subroutine nml_read_double_vector 

    subroutine nml_read_float_vector(filename,group,name,value,comment,init)

        implicit none 

        real(4), intent(INOUT) :: value(:) 
        character(len=*), intent(IN)    :: filename, group, name 
        character(len=*), intent(INOUT), optional :: comment   
        character(len=256) :: value_str, value_str_vec(size(value)) 
        logical, optional :: init  
        logical :: init_var
        integer :: q 

        init_var = .FALSE. 
        if (present(init)) init_var = init 
        if (init_var) value(:) = 0.0

        ! First find parameter value as a string 
        value_str = ""
        call nml_read_internal(filename,group,name,value_str,comment)

        if (value_str /= "") then
            call string_to_vector(value_str,value_str_vec)
            do q = 1, size(value)
                if (trim(value_str_vec(q)) /= "") then 
                    value(q) = real(string_to_double(trim(adjustl(value_str_vec(q)))))
                end if 

            end do 
        end if 

        if (VERBOSE) call nml_print(name,value,comment)  ! Check

        return 

    end subroutine nml_read_float_vector 

    subroutine nml_read_integer_vector(filename,group,name,value,comment,init)

        implicit none 

        integer, intent(INOUT) :: value(:) 
        character(len=*), intent(IN)    :: filename, group, name 
        character(len=*), intent(INOUT), optional :: comment   
        character(len=256) :: value_str, value_str_vec(size(value))
        logical, optional :: init  
        logical :: init_var 
        integer :: q 

        init_var = .FALSE. 
        if (present(init)) init_var = init 
        if (init_var) value(:) = 0

        ! First find parameter value as a string 
        value_str = ""
        call nml_read_internal(filename,group,name,value_str,comment)

        if (value_str /= "") then
            call string_to_vector(value_str,value_str_vec)
            do q = 1, size(value)
                if (trim(value_str_vec(q)) /= "") then 
                    value(q) = nint(string_to_double(trim(adjustl(value_str_vec(q)))))
                end if 

            end do 
        end if 

        if (VERBOSE) call nml_print(name,value,comment)  ! Check

        return 

    end subroutine nml_read_integer_vector

    subroutine nml_read_logical_vector(filename,group,name,value,comment,init)

        implicit none 

        logical, intent(INOUT) :: value(:) 
        character(len=*), intent(IN)    :: filename, group, name 
        character(len=*), intent(INOUT), optional :: comment   
        character(len=256) :: value_str, value_str_vec(size(value))
        logical, optional :: init  
        logical :: init_var 
        integer :: q 

        init_var = .FALSE. 
        if (present(init)) init_var = init 
        if (init_var) value(:) = .FALSE.

        ! First find parameter value as a string 
        value_str = ""
        call nml_read_internal(filename,group,name,value_str,comment)

        if (value_str /= "") then
            call string_to_vector(value_str,value_str_vec)
            do q = 1, size(value)
                if (trim(value_str_vec(q)) /= "") then 
                    value(q) = string_to_logical(trim(adjustl(value_str_vec(q))))
                end if 

            end do 
        end if 

        if (VERBOSE) call nml_print(name,value,comment)  ! Check

        return 

    end subroutine nml_read_logical_vector

    ! =============================================================
    !
    ! nml parameter writing functions
    !
    ! =============================================================

    ! This is the basic nml writing subroutine
    ! All interfaces use this to write the parameter to file, after
    ! converting the parameter to a string
    subroutine nml_write_string(filename,group,name,value,comment)

        implicit none 

        character(len=*), intent(IN) :: filename, group, name, value 
        character(len=*), intent(IN), optional :: comment  


        return 

    end subroutine nml_write_string 

    ! =============================================================
    !
    ! nml line printing functions
    !
    ! =============================================================

    ! This is the basic routine for printing a parameter to a formatted line
    ! All other interfaces use this routine after converting to a string.
    subroutine nml_print_string(name,value,comment,io,no_quotes)

        implicit none 
        character(len=*), intent(in) :: name, value 
        character(len=*), intent(in), optional :: comment 
        integer, intent(in), optional :: io 
        integer :: io_val 
        logical, intent(in), optional :: no_quotes
        logical :: no_quotes_loc
        character(len=1000) :: line
        character(len=500)  :: comment1 
        character(len=len(value))  :: val_repr

        io_val = 6 
        if (present(io)) io_val = io 
        val_repr = value

        no_quotes_loc = .false.
        if (present(no_quotes)) no_quotes_loc = no_quotes
        if (.not.no_quotes_loc) val_repr = '"'//trim(value)//'"'

        comment1 = "" 
        if (present(comment)) comment1 = "   "//trim(comment)

        write(line,"(a)") "    "//trim(name)//" = "//trim(val_repr)//trim(comment1)
        write(io_val,*) trim(line)

        return 

    end subroutine nml_print_string 

    subroutine nml_print_double(name,value,comment,io)

        implicit none 
        double precision :: value
        character(len=*) :: name 
        character(len=*), optional :: comment
        integer, optional :: io 
        character(len=500) :: value_str  

        write(value_str,*) value 
        if (VERBOSE) call nml_print_string(name,value_str,comment,io,no_quotes=.true.)

        return 

    end subroutine nml_print_double
    
    subroutine nml_print_float(name,value,comment,io)

        implicit none 
        real(4) :: value
        character(len=*) :: name 
        character(len=*), optional :: comment
        integer, optional :: io 
        character(len=500) :: value_str  

        write(value_str,*) value 
        if (VERBOSE) call nml_print_string(name,value_str,comment,io)

        return 

    end subroutine nml_print_float
    
    subroutine nml_print_integer(name,value,comment,io)

        implicit none 
        integer :: value
        character(len=*) :: name 
        character(len=*), optional :: comment
        integer, optional :: io 
        character(len=500) :: value_str  

        write(value_str,*) value 
        if (VERBOSE) call nml_print_string(name,value_str,comment,io,no_quotes=.true.)

        return 

    end subroutine nml_print_integer

    subroutine nml_print_logical(name,value,comment,io)

        implicit none 
        logical :: value
        character(len=*) :: name 
        character(len=*), optional :: comment
        integer, optional :: io 
        character(len=500) :: value_str  

        value_str = "F"
        if (value) value_str = "T" 
        if (VERBOSE) call nml_print_string(name,value_str,comment,io,no_quotes=.true.)

        return 

    end subroutine nml_print_logical
    
    !! Vectors

    subroutine nml_print_string_vector(name,value,comment,io)

        implicit none 
        character(len=*) :: value(:)
        character(len=*) :: name 
        character(len=*), optional :: comment
        integer, optional :: io 
        character(len=500) :: value_str  
        integer :: q 

        value_str = '"'//trim(value(1))//'"'
        do q = 2, size(value)
            write(value_str,*) trim(value_str)//" "//'"'//trim(value(q))//'"'
        end do 

        if (VERBOSE) call nml_print_string(name,value_str,comment,io,no_quotes=.true.)

        return 

    end subroutine nml_print_string_vector
    
    subroutine nml_print_double_vector(name,value,comment,io)

        implicit none 
        double precision :: value(:)
        character(len=*) :: name 
        character(len=*), optional :: comment
        integer, optional :: io 
        character(len=500) :: value_str  
        integer :: q 

        value_str = ""
        do q = 1, size(value)
            write(value_str,"(a,g12.3)") trim(value_str)//" ",value(q)
        end do 

        if (VERBOSE) call nml_print_string(name,value_str,comment,io,no_quotes=.true.)

        return 

    end subroutine nml_print_double_vector
    
    subroutine nml_print_float_vector(name,value,comment,io)

        implicit none 
        real(4) :: value(:)
        character(len=*) :: name 
        character(len=*), optional :: comment
        integer, optional :: io 
        character(len=500) :: value_str  
        integer :: q 

        value_str = ""
        do q = 1, size(value)
            write(value_str,"(a,g12.3)") trim(value_str)//" ",value(q)
        end do 

        if (VERBOSE) call nml_print_string(name,value_str,comment,io,no_quotes=.true.)

        return 

    end subroutine nml_print_float_vector
    
    subroutine nml_print_integer_vector(name,value,comment,io)

        implicit none 
        integer :: value(:)
        character(len=*) :: name 
        character(len=*), optional :: comment
        integer, optional :: io 
        character(len=500) :: value_str  
        integer :: q 

        value_str = ""
        do q = 1, size(value)
            write(value_str,"(a,i12)") trim(value_str)//" ",value(q)
        end do 

        if (VERBOSE) call nml_print_string(name,value_str,comment,io,no_quotes=.true.)

        return 

    end subroutine nml_print_integer_vector
    
    subroutine nml_print_logical_vector(name,value,comment,io)

        implicit none 
        logical :: value(:)
        character(len=*) :: name 
        character(len=*), optional :: comment
        integer, optional :: io 
        character(len=500) :: value_str  
        integer :: q 

        value_str = ""
        do q = 1, size(value)
            if (value(q)) then 
                write(value_str,"(a,a1)") trim(value_str)//" ","T"
            else
                write(value_str,"(a,a1)") trim(value_str)//" ","F"
            end if 
        end do 

        if (VERBOSE) call nml_print_string(name,value_str,comment,io,no_quotes=.true.)

        return 

    end subroutine nml_print_logical_vector
    

    ! =============================================================
    !
    ! Type conversion functions
    !
    ! =============================================================

    function string_to_double(string) result(value)

        implicit none 

        character(len=*), intent(IN) :: string 
        double precision :: value 

        character(len=256) :: tmpstr 
        integer :: stat, n
        double precision :: x 

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

    end function string_to_double

    function string_to_logical(string) result(value)

        implicit none 

        character(len=*), intent(IN) :: string 
        logical :: value 

        character(len=256) :: tmpstr 
        integer :: stat, n
        double precision :: x 

        tmpstr = trim(adjustl(string))
        
        select case(trim(tmpstr))
            case("T","True","TRUE","true",".TRUE.", ".true.")
                value = .TRUE. 
            case("F","False","FALSE","false",".FALSE.", ".false.")
                value = .FALSE. 
            case DEFAULT
                write(*,*) "nml:: Error reading logical parameter."
                stop 
        end select  

        return 

    end function string_to_logical

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

    ! =============================================================
    !
    ! Helper functions
    !
    ! =============================================================

    ! This function parses a namelist line into pieces, determining
    ! whether it is a blank line (-2), comment (-1), group name (1)
    ! end-of-group (2), or a parameter line (3)
    subroutine parse_line(line,linetype,name,value,comment)

        implicit none 
        character(len=*), intent(IN)    :: line
        character(len=*), intent(INOUT) :: name, value, comment 
        integer :: linetype

        character(len=1001) :: line1 
        integer :: q, q1, q2

        name     = ""
        value    = ""
        comment  = "" 

        line1 = trim(adjustl(line))

        if (trim(line1) == "") then         ! Blank line
            linetype = -2 
            
        else if (line1(1:1) == "!") then    ! Comment line 
            linetype = -1 
            comment = trim(line1)

        else if (line1(1:1) == "&") then    ! Group name 
            linetype = 1 
            q = len_trim(line1)
            name     = line1(2:q)

        else if (line1(1:1) == "/") then    ! End of group 
            linetype = 2 

        else   ! Line must contain parameter to read
            linetype = 3

            q = index(line1,"=")
            if (q == 0) then 
                write(*,*) "nml:: Error reading namelist file."
                write(*,*) "No '=' found on parameter line."
                stop 
            end if 

            name = trim(adjustl(line1(1:q-1)))

            q1 = index(line1,"!")
            q2 = len_trim(line1)

            if (q1 > 0) then 
                comment = trim(adjustl(line1(q1:q2)))
                value   = trim(adjustl(line1(q+1:q1-1)))
            else
                value   = trim(adjustl(line1(q+1:q2)))
            end if 

            ! Remove quotes around string, and final line comma, if they exist
            call remove_quotes_comma(value)

        end if 

        return 

    end subroutine parse_line 

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



    integer function newunit(unit)
        ! This is a simple function to search for an available unit.
        ! LUN_MIN and LUN_MAX define the range of possible LUNs to check.
        ! The UNIT value is returned by the function, and also by the optional
        ! argument. This allows the function to be used directly in an OPEN
        ! statement, and optionally save the result in a local variable.
        ! If no units are available, -1 is returned.
        ! Modified from: http://fortranwiki.org/fortran/show/newunit

        implicit none 

        integer, intent(out), optional :: unit

        ! Local variables
        integer, parameter :: LUN_MIN=10, LUN_MAX=1000
        logical :: opened
        integer :: lun

        newunit=-1

        ! Search for new unit
        do lun=LUN_MIN,LUN_MAX
            inquire(unit=lun,opened=opened)
            if (.not. opened) then
                newunit=lun
                exit
            end if
        end do

        ! Assign new unit if desired
        if (present(unit)) unit=newunit

        return 

    end function newunit

end module nml 

