program test_variables

    implicit none

    type var_config_type
        character(len=32)  :: varname
        character(len=32)  :: dimnames
        character(len=32)  :: units
        character(len=128) :: long_name

        integer :: ndims
        character(len=32), allocatable :: dims(:)
    end type

    type(var_config_type), dimension(:), allocatable :: config_table

   
    call load_config_table(config_table,"input/yelmo_variables_ytopo.txt")

contains

    subroutine load_config_table(config_table,filename)

        implicit none

        type(var_config_type), allocatable, intent(OUT) :: config_table(:)
        character(len=*), intent(IN) :: filename

        ! Local variables
        integer, parameter :: nmax = 1000
        type(var_config_type) :: ct(nmax)       ! Max 1000 variables
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
                call parse_line_to_variable(ct(n)%varname,ct(n)%dimnames,ct(n)%units,ct(n)%long_name,line)
                
                ! Parse the dimensions
                call parse_dims(ct(n)%dims,ct(n)%dimnames)

                ntot = ntot+1
                write(*,"(i4,3a20,2x,a50)") ntot, trim(ct(n)%varname), trim(ct(n)%dimnames), trim(ct(n)%units), trim(ct(n)%long_name)
            end if

            if (n .eq. nmax) then
                write(*,*) "load_config_table:: Error: maximum number of variables read in, but &
                &end of the table was not reached. Increase value of nmax in this routine to be safe."
                write(*,*) "nmax = ", nmax 
                write(*,*) "filename = ", trim(filename)
                write(*,*) "Last line read in: "
                write(*,*) trim(line)
                stop
            end if

        end do

        if (allocated(config_table)) deallocate(config_table)
        allocate(config_table(ntot))
        config_table = ct(1:ntot)

        return

    end subroutine load_config_table

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

end program test_variables