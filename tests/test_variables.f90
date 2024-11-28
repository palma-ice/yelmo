program test_variables

    implicit none

    type :: var_config_type
        character(len=32)  :: varname
        character(len=32)  :: dims
        character(len=32)  :: units
        character(len=128) :: long_name
        integer :: ndims
    end type

    type(var_config_type), dimension(:), allocatable :: config_table

   
    call load_config_table(config_table,"input/yelmo_variables_ytopo.txt")

contains

    subroutine load_config_table(config_table,filename)

        implicit none

        type(var_config_type), allocatable, intent(OUT) :: config_table(:)
        character(len=*), intent(IN) :: filename

        ! Local variables
        type(var_config_type) :: ct(1000)       ! Max 1000 variables
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

        do n = 1, 1000
            read(unit=io, fmt='(a)') line

            if (line(1:1) .eq. "+") then
                ! End of table reached, exit 
                exit
            else
                ! Parse the line
                call parse_line_to_variable(ct(n)%varname,ct(n)%dims,ct(n)%units,ct(n)%long_name,line)
                ntot = ntot+1
                write(*,"(i4,3a20,2x,a50)") ntot, trim(ct(n)%varname), trim(ct(n)%dims), trim(ct(n)%units), trim(ct(n)%long_name)
            end if
        end do

        if (allocated(config_table)) deallocate(config_table)
        allocate(config_table(ntot))
        config_table = ct(1:ntot)

        return

    end subroutine load_config_table

    subroutine parse_line_to_variable(varname,dims,units,long_name,line)

        implicit none

        character(len=*), intent(OUT) :: varname
        character(len=*), intent(OUT) :: dims
        character(len=*), intent(OUT) :: units
        character(len=*), intent(OUT) :: long_name
        character(len=*), intent(IN)  :: line
        ! Local variables
        character(len=200) :: temp_line
        integer :: pos1, pos2, length

        ! Initialize
        varname     = ""
        dims        = ""
        units       = ""
        long_name   = ""

        ! Remove leading/trailing spaces
        temp_line = adjustl(trim(line))

        ! Extract var_name
        pos1 = index(temp_line, '|') + 1
        temp_line = temp_line(pos1:)
        pos2 = index(temp_line, '|') - 1
        varname = adjustl(trim(temp_line(:pos2)))

        ! Extract dims
        temp_line = temp_line(pos2+2:)
        pos2 = index(temp_line, '|') - 1
        dims = adjustl(trim(temp_line(:pos2)))

        ! Extract units
        temp_line = temp_line(pos2+2:)
        pos2 = index(temp_line, '|') - 1
        units = adjustl(trim(temp_line(:pos2)))

        ! Extract long_name
        temp_line = temp_line(pos2+2:)
        long_name = adjustl(trim(temp_line))

        return

    end subroutine parse_line_to_variable

end program test_variables