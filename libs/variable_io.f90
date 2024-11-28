module variable_io

    implicit none

    type var_io_type
        integer :: id
        character(len=32)  :: varname
        character(len=32)  :: dimnames
        character(len=32)  :: units
        character(len=128) :: long_name

        integer :: ndims
        character(len=32), allocatable :: dims(:)
    end type

    private
    public :: var_io_type
    public :: load_variable_table
    public :: find_variable_in_table
    public :: var_io_print
    
contains

    subroutine find_variable_in_table(var,var_table,varname)

        implicit none
        
        type(var_io_type), intent(OUT) :: var
        type(var_io_type), intent(IN)  :: var_table(:)
        character(len=*),  intent(IN)  :: varname 

        ! Local variables
        integer :: n 
        logical :: found

        found = .FALSE.
        do n = 1, size(var_table)
            if (trim(var_table(n)%varname) .eq. trim(varname)) then
                var = var_table(n)
                found = .TRUE.
                exit
            end if
        end do

        if (.not. found) then
            write(*,*) "find_variable_in_table:: Error: variable not found in variable table."
            write(*,*) "varname = ", trim(varname)
            stop 
        end if

        return

    end subroutine find_variable_in_table

    subroutine load_variable_table(var_table,filename)

        implicit none

        type(var_io_type), allocatable, intent(OUT) :: var_table(:)
        character(len=*), intent(IN) :: filename

        ! Local variables
        integer, parameter :: nmax = 1000
        type(var_io_type) :: vt(nmax)       ! Max 1000 variables
        integer :: n, ntot 
        character(len=1024) :: line 
        integer, parameter :: io = 55

        ntot = 0

        ! Open the file for reading
        open(unit=io, file=trim(filename), status='old', action='read')

        ! Skip the first three header lines
        read(unit=io, fmt='(a)') line
        read(unit=io, fmt='(a)') line
        read(unit=io, fmt='(a)') line

        do n = 1, nmax
            read(unit=io, fmt='(a)') line

            if (line(1:1) .eq. "+") then
                ! End of table reached, exit 
                exit
            else
                ! Parse the line
                call parse_line_to_variable(vt(n)%varname,vt(n)%dimnames,vt(n)%units,vt(n)%long_name,line)
                
                ! Store the id (variable number in table)
                vt(n)%id = n 
                
                ! Parse the dimensions
                call parse_dims(vt(n)%dims,vt(n)%dimnames)
                vt(n)%ndims = size(vt(n)%dims)
                ntot = ntot+1
                call var_io_print(vt(n))
            end if

            if (n .eq. nmax) then
                write(*,*) "load_var_table:: Error: maximum number of variables read in, but &
                &end of the table was not reached. Increase value of nmax in this routine to be safe."
                write(*,*) "nmax = ", nmax 
                write(*,*) "filename = ", trim(filename)
                write(*,*) "Last line read in: "
                write(*,*) trim(line)
                stop
            end if

        end do

        if (allocated(var_table)) deallocate(var_table)
        allocate(var_table(ntot))
        var_table = vt(1:ntot)

        return

    end subroutine load_variable_table

    subroutine parse_line_to_variable(varname,dimnames,units,long_name,line)

        implicit none

        character(len=*), intent(OUT) :: varname
        character(len=*), intent(OUT) :: dimnames
        character(len=*), intent(OUT) :: units
        character(len=*), intent(OUT) :: long_name
        character(len=*), intent(IN)  :: line
        ! Local variables
        character(len=200) :: temp_line
        integer :: pos1, pos2, length

        ! Initialize
        varname     = ""
        dimnames    = ""
        units       = ""
        long_name   = ""

        ! Remove leading/trailing spaces
        temp_line = adjustl(trim(line))

        ! Extract var_name
        pos1 = index(temp_line, '|') + 1
        temp_line = temp_line(pos1:)
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
        long_name = adjustl(trim(temp_line))

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
