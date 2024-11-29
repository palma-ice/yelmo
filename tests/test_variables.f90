program test_variables

    use variable_io 

    implicit none

    type(var_io_type), allocatable :: var_table(:)
    type(var_io_type) :: var

    !call load_var_io_table(var_table,"input/yelmo-variables-ytopo.md")
    !call find_var_io_in_table(var,"H_ice",var_table)

    write(*,*)
    call find_var_io_in_file(var,"H_ice","input/yelmo-variables-ytopo.md")
    call var_io_print(var)
    write(*,*)
    call find_var_io_in_file(var,"uxy_s","input/yelmo-variables-ydyn.md")
    call var_io_print(var)
    write(*,*)
    call find_var_io_in_file(var,"visc","input/yelmo-variables-ymat.md")
    call var_io_print(var)
    write(*,*)
    call find_var_io_in_file(var,"T_ice","input/yelmo-variables-ytherm.md")
    call var_io_print(var)

end program test_variables