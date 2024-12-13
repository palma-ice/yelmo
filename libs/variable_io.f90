module variable_io

    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit

    implicit none

    type var_io_type
        integer :: id
        character(len=32)  :: varname
        character(len=32)  :: dimnames
        character(len=32)  :: units
        character(len=256) :: long_name

        integer :: ndims
        character(len=32), allocatable :: dims(:)
    end type

    private
    public :: var_io_type
    public :: find_var_io_in_file
    public :: find_var_io_in_table
    public :: load_var_io_table
    public :: var_io_print
    
contains

    subroutine find_var_io_in_file(var,varname,filename)

        implicit none

        type(var_io_type), intent(OUT) :: var
        character(len=*),  intent(IN)  :: varname
        character(len=*),  intent(IN)  :: filename

        ! Local variables 
        type(var_io_type), allocatable :: var_table(:)
        
        ! Load the variable table
        call load_var_io_table(var_table,filename)

        ! Find the current variable in the table
        call find_var_io_in_table(var,varname,var_table)

        return

    end subroutine find_var_io_in_file

    subroutine find_var_io_in_table(var,varname,var_table,with_error)

        implicit none
        
        type(var_io_type), intent(OUT) :: var
        character(len=*),  intent(IN)  :: varname 
        type(var_io_type), intent(IN)  :: var_table(:)
        logical,           intent(IN), optional :: with_error

        ! Local variables
        integer :: n 
        logical :: found

        found = .FALSE.
        var%varname = "none"

        do n = 1, size(var_table)
            if (trim(var_table(n)%varname) .eq. trim(varname)) then
                var = var_table(n)
                found = .TRUE.
                exit
            end if
        end do

        if (.not. found .and. present(with_error)) then
            if (with_error) then
                write(error_unit,*) "find_var_io_in_table:: Error: variable not found in variable table."
                write(error_unit,*) "varname = ", trim(varname)
                stop
            end if
        end if

        return

    end subroutine find_var_io_in_table

    subroutine load_var_io_table(var_table,filename)

        implicit none

        type(var_io_type), allocatable, intent(OUT) :: var_table(:)
        character(len=*), intent(IN) :: filename

        ! Local variables
        integer, parameter :: nmax = 1000
        type(var_io_type) :: vt(nmax)       ! Max 1000 variables
        integer :: n, ntot 
        character(len=1024) :: line
        integer :: io, check

        ! Open the file for reading
        open(newunit=io, file=trim(filename), status='old', action='read')

        ! Skip the header (header finishes with first line starting with "---" )
        do n = 1, 50
            read(unit=io, fmt='(a)') line
            if (line(1:3) .eq. "|--") exit
        end do

        ntot = 0 
        do n = 1, nmax
            read(unit=io, fmt='(a)',iostat=check) line

            ! Exit at end of file
            if (check .eq. -1) exit 

            if (line(1:2) .eq. "| ") then
                ! Parse the line
                call parse_line_to_variable(vt(n)%id,vt(n)%varname,vt(n)%dimnames,vt(n)%units,vt(n)%long_name,line)
                
                ! Parse the dimensions
                call parse_dims(vt(n)%dims,vt(n)%dimnames)
                vt(n)%ndims = size(vt(n)%dims)
                ntot = ntot+1
                call var_io_print(vt(n))
            end if

            if (n .eq. nmax) then
                write(error_unit,*) "load_var_table:: Error: maximum number of variables read in, but &
                &end of the table was not reached. Increase value of nmax in this routine to be safe."
                write(error_unit,*) "nmax = ", nmax 
                write(error_unit,*) "filename = ", trim(filename)
                write(error_unit,*) "Last line read in: "
                write(error_unit,*) trim(line)
                stop
            end if

        end do

        ! Close the file
        close(io)

        if (allocated(var_table)) deallocate(var_table)
        allocate(var_table(ntot))
        var_table = vt(1:ntot)

        return

    end subroutine load_var_io_table

    subroutine parse_line_to_variable(id,varname,dimnames,units,long_name,line)

        implicit none

        integer,          intent(OUT) :: id
        character(len=*), intent(OUT) :: varname
        character(len=*), intent(OUT) :: dimnames
        character(len=*), intent(OUT) :: units
        character(len=*), intent(OUT) :: long_name
        character(len=*), intent(IN)  :: line
        ! Local variables
        character(len=200) :: temp_line
        character(len=20)  :: id_str
        integer :: pos1, pos2, length
        integer :: ios 

        ! Initialize
        id_str      = ""
        varname     = ""
        dimnames    = ""
        units       = ""
        long_name   = ""

        ! Remove leading/trailing spaces
        temp_line = adjustl(trim(line))

        ! Extract id
        pos1 = index(temp_line, '|') + 1
        temp_line = temp_line(pos1:)
        pos2 = index(temp_line, '|') - 1
        id_str = adjustl(trim(temp_line(:pos2)))

        if (.not. trim(id_str) .eq. "") then
            ! Convert to integer
            read(id_str,"(i20)",iostat=ios) id
            if (ios .ne. 0) then
                write(error_unit,*) "variable_io:: parse_line_to_variable:: Error: id must be an integer."
                write(error_unit,*) "id = ", trim(id_str)
                stop
            end if
        else
            id = 0
        end if

        ! Extract var_name
        temp_line = temp_line(pos2+2:)
        pos2 = index(temp_line, '|') - 1
        varname = adjustl(trim(temp_line(:pos2)))

        ! Extract dimnames
        temp_line = temp_line(pos2+2:)
        pos2 = index(temp_line, '|') - 1
        dimnames = adjustl(trim(temp_line(:pos2)))

        ! Extract units
        temp_line = temp_line(pos2+2:)
        pos2 = index(temp_line, '|') - 1
        units = adjustl(trim(temp_line(:pos2)))

        ! Extract long_name
        temp_line = temp_line(pos2+2:)
        pos2 = index(temp_line, '|') - 1
        long_name = adjustl(trim(temp_line(:pos2)))

        return

    end subroutine parse_line_to_variable

    subroutine parse_dims(dims,dimnames)

        implicit none

        character(len=*), allocatable, intent(OUT) :: dims(:)
        character(len=*), intent(IN)  :: dimnames

        ! Local variables
        character(len=200) :: temp_line
        integer :: pos1, pos2, length
        character(len=32) :: d(10)
        integer :: n, ndims 

        ! Initialize
        d(:) = ""

        ! Remove leading/trailing spaces
        temp_line = adjustl(trim(dimnames))

            
        ndims = 0

        do n = 1, size(d)

            pos1 = index(temp_line, ',')

            if (pos1 .eq. 0) then
                ! Only one dimension currently available
                ndims = n
                d(ndims) = adjustl(trim(temp_line))
                exit

            else
                ndims = n
                d(ndims) = adjustl(trim(temp_line(1:pos1-1)))
                temp_line = temp_line(pos1+1:)
            end if

        end do

        if (allocated(dims)) deallocate(dims)
        allocate(dims(ndims))
        dims = d(1:ndims)

        return

    end subroutine parse_dims

    subroutine var_io_print(var)

        implicit none

        type(var_io_type), intent(IN) :: var

        write(*,"(i4,3a20,2x,a50)") var%id, trim(var%varname), trim(var%dimnames), trim(var%units), trim(var%long_name)

        return

    end subroutine var_io_print

end module variable_io
