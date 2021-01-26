module root_finder

    implicit none 

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the working precision of the library (sp,dp)
    integer,  parameter :: wp = sp 

    private
    public :: solve_newton 
    public :: solve_secant 

contains
    
    subroutine solve_newton(x,n_iter,x_init,tol,n_max,f,fp,debug)
        ! Estimate the zero of f(x) using Newton's method. 
        ! Adapted from: 
        ! https://faculty.washington.edu/rjl/classes/am583s2013/notes/fortran_newton.html

        implicit none

        real(wp), intent(OUT) :: x          ! Best guess of root
        integer,  intent(OUT) :: n_iter     ! Number of iterations to reach it
        real(wp), intent(IN)  :: x_init     ! Initial guess
        real(wp), intent(IN)  :: tol        ! Tolerance to convergence
        integer,  intent(IN)  :: n_max      ! Maximum iterations allowed
        real(wp), external    :: f          ! Function to find root of 
        real(wp), external    :: fp         ! Derivative of Function
        logical,  intent(IN)  :: debug      ! Print iteration information?
        
        ! Declare any local variables:
        real(wp) :: deltax, fx, fxprime
        integer    :: k

        ! Save initial guess
        x = x_init

        ! Newton iteration to find a zero of f(x) 
        n_iter = 0 

        do k = 1, n_max

            n_iter = n_iter + 1 

            ! evaluate function and its derivative:
            fx      = f(x)
            fxprime = fp(x)

            if (debug) then
                write(*,*) n_iter, "x, f(x) = ", x, fx
            end if 

            if (abs(fx) < tol) then
                exit  ! jump out of do loop
            end if

            ! Compute Newton increment x:
            deltax = fx/fxprime

            ! update x:
            x = x - deltax

        end do


        if (n_iter .eq. n_max .and. abs(fx) > tol) then
            write(*,*) "solve_newton:: Warning: no convergence."
        end if

        return 

    end subroutine solve_newton

    subroutine solve_secant(x,n_iter,x_init,tol,n_max,f,debug) 
        ! Estimate the zero of f(x) using the Secand method. 
        ! Adapted from: 
        ! http://jean-pierre.moreau.pagesperso-orange.fr/Fortran/secant_f90.txt
        ! https://rosettacode.org/wiki/Roots_of_a_function#Fortran

        implicit none 

        real(wp), intent(OUT) :: x          ! Best guess of root
        integer,  intent(OUT) :: n_iter     ! Number of iterations to reach it
        real(wp), intent(IN)  :: x_init     ! Initial guess
        real(wp), intent(IN)  :: tol        ! Tolerance to convergence
        integer,  intent(IN)  :: n_max      ! Maximum iterations allowed
        real(wp), external    :: f          ! Function to find root of 
        logical,  intent(IN)  :: debug      ! Print iteration information?
        
        ! Local variables 
        integer  :: n 
        real(wp) :: x1, x2  
        real(wp) :: y1, y2  
        real(wp) :: d 

        ! Set x to initial guess and a slightly different value
        x1 = x_init 
        x2 = x_init*0.5_wp
        if (x_init .eq. 0.0_wp) x2 = x_init + 0.1_wp  

        n_iter = 0

        ! Start iterations
        do n = 1, n_max 

            n_iter = n_iter + 1 

            ! Calculate new value of y at x
            y1 = f(x1)
            y2 = f(x2)

            if (debug) write(*,*) n_iter, x1, x2, y1, y2 

            if (abs(y2) < tol) exit

            
            d = (x2 - x1) / (y2 - y1) * y2
            
            x1 = x2
            x2 = x2 - d

        end do 

        x = x2 

        if (n_iter .eq. n_max .and. abs(y2) > tol) then 
            write(*,*) "solve_secant:: Warning: no convergence."
        end if 

        return 

    end subroutine solve_secant

end module root_finder
