module yelmo_c_api
  
  use yelmo
  
  implicit none

  type(yelmo_class) :: ylmo   ! persistent state

contains

  subroutine yelmo_init_wrapper(filename, grid_def, time) bind(C, name="yelmo_init")
    use iso_c_binding
    character(c_char), intent(in) :: filename(*)
    character(c_char), intent(in) :: grid_def(*)
    real(c_double), value         :: time

    ! Local variables
    character(len=1028)           :: f_filename
    character(len=56)             :: f_grid_def
    integer                       :: i

    ! Convert C string to Fortran string (null-terminated scan)
    f_filename    = trim(c_to_f_string(filename))
    f_grid_def    = trim(c_to_f_string(grid_def))

    call yelmo_init(ylmo, filename=trim(f_filename), grid_def=trim(f_grid_def), time=real(time, wp))
  end subroutine


  subroutine yelmo_init_state_wrapper(time, thrm_method) bind(C, name="yelmo_init_state")
    use iso_c_binding
    real(c_double), value         :: time
    character(c_char), intent(in) :: thrm_method(*) ! e.g., "robin-cold"

    ! Local variables
    character(len=56) :: f_thrm_method
    integer           :: i

    ! Convert C string to Fortran string (null-terminated scan)
    f_thrm_method    = trim(c_to_f_string(thrm_method))

    call yelmo_init_state(ylmo, time=real(time, wp), thrm_method=trim(f_thrm_method))
  end subroutine

  subroutine yelmo_step_wrapper(time) bind(C, name="yelmo_step")
    use iso_c_binding
    real(c_double), value :: time
    call yelmo_update(ylmo, real(time, wp))
  end subroutine

  ! Example: pull ice thickness out to Julia
  subroutine yelmo_get_H_ice(H, nx, ny) bind(C, name="yelmo_get_H_ice")
    use iso_c_binding
    integer(c_int), value   :: nx, ny
    real(c_double), intent(out) :: H(nx, ny)
    H = real(ylmo%tpo%now%H_ice, c_double)
  end subroutine

  ! Example: push a field in from Julia (e.g. SMB, basal melt)
  subroutine yelmo_set_bmb(bmb, nx, ny) bind(C, name="yelmo_set_bmb")
    use iso_c_binding
    integer(c_int), value  :: nx, ny
    real(c_double), intent(in) :: bmb(nx, ny)
    ylmo%bnd%bmb_shlf = real(bmb, wp)
  end subroutine


  function c_to_f_string(c_str) result(f_str)
    use iso_c_binding
    character(c_char), intent(in) :: c_str(*)
    character(len=1028)           :: f_str
    integer                       :: i

    f_str = " "
    do i = 1, 1028
      if (c_str(i) == c_null_char) exit
      f_str(i:i) = c_str(i)
    end do
  end function

end module yelmo_c_api