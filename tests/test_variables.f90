program test_variables

    use variable_io 

    implicit none

    type(var_io_type), allocatable :: var_table(:)
    type(var_io_type) :: var

    !call load_variable_table(var_table,"input/yelmo-variables-ytopo.txt")
    !call find_variable_in_table(var,"H_ice",var_table)

    call find_variable_in_file(var,"H_ice","input/yelmo-variables-ytopo.txt")
    
    write(*,*)
    call var_io_print(var)

end program test_variables