program test_roots

    use root_finder

    implicit none 

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the working precision of the library (sp,dp)
    integer,  parameter :: wp = sp 


    integer  :: n_iter, n_max
    real(wp) :: x, x_init, tol

    n_max  = 100 
    x_init = -2.0_wp 
    tol    = real(1e-5,wp)

    call solve_secant(x,n_iter,x_init,tol,n_max,Y,.FALSE.) 

    write(*,*) 
    write(*,*) "================================"
    write(*,*)
    write(*,*) "Secant method:"
    write(*,*)
    write(*,"(a,f9.6)") ' The calculated zero is X    = ', x
    write(*,"(a,f9.6)") ' The associated Y value is Y = ', Y(x)
    write(*,"(a,i4)")   ' The number of steps was:      ', n_iter 

    call solve_newton(x,n_iter,x_init,tol,n_max,Y,Yp,.FALSE.)

    write(*,*) 
    write(*,*) "================================"
    write(*,*)
    write(*,*) "Newton method:"
    write(*,*)
    write(*,"(a,f9.6)") ' The calculated zero is X    = ', x
    write(*,"(a,f9.6)") ' The associated Y value is Y = ', Y(x)
    write(*,"(a,i4)")   ' The number of steps was:      ', n_iter 

contains
                      
! real(wp) function Y(x)
!   real(wp) :: x 
!   Y = 1.0_wp+5.0_wp*x+10.0_wp*x*x+10.0_wp*x*x*x+5.0_wp*x*x*x*x+x*x*x*x*x
! end function Y

! real(wp) function Yp(x)
!   real(wp) :: x 
!   Yp = 5.0_wp+20.0_wp*x+30.0_wp*x*x+45.0_wp*x*x*x*x+x*x*x*x
! end function Yp


real(wp) function Y(x)
  real(wp) :: x 
  Y = (x+1.0_wp)**5
end function Y

real(wp) function Yp(x)
  real(wp) :: x 
  Yp = 5.0_wp*(x+1.0_wp)**4
end function Yp

end program test_roots
