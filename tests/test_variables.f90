program test_variables

    use variable_io 

    implicit none

    
    type(var_io_type), dimension(:), allocatable :: var_table

   
    call load_variable_table(var_table,"input/yelmo-variables-ytopo.txt")

end program test_variables