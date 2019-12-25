
module interp1D

    use yelmo_defs 

    implicit none 

    interface interp_linear
        module procedure interp_linear_pt, interp_linear_vec 
    end interface 

    private
    public :: interp_linear
    public :: interp_spline
    public :: interp_align 

contains

    function interp_align(x,y,xout,missing_value,tol) result(yout)
        ! Fill in a vector of times with the nearest input 
        ! times to within a certain tolerance
        implicit none 
 
        real(prec), dimension(:), intent(IN) :: x, y
        real(prec), dimension(:), intent(IN) :: xout
        real(prec), intent(IN), optional :: missing_value, tol 
        real(prec), dimension(size(xout)) :: yout
        real(prec) :: tolerance 
        integer :: i, n, nout 

        ! Length of output vector
        nout = size(xout)

        ! Define tolerance for filling points 
!         tolerance = 0.1d0
        tolerance = minval(xout(2:nout)-xout(1:nout-1),dim=1)*0.49d0 
        if (present(tol)) tolerance = tol 

        ! Define missing values and fill vector intially
        yout = MISSING_VALUE_DEFAULT
        if (present(missing_value)) yout = missing_value

        do i = 1, nout 
            n = minloc(abs(x-xout(i)),dim=1)
            if (abs(x(n)-xout(i)) <= tolerance) yout(i) = y(n) 
        end do 
        
        return

    end function interp_align

    function interp_linear_pt(x,y,xout) result(yout)
        ! Interpolate y from ordered x to ordered xout positions

        implicit none 
 
        real(prec), dimension(:), intent(IN) :: x, y
        real(prec), intent(IN) :: xout
        real(prec) :: yout 
        integer :: i, j, n, nout 

        n    = size(x) 

        if (xout .lt. x(1)) then
            yout = y(1)
        else if (xout .gt. x(n)) then
            yout = y(n)
        else
            do j = 1, n 
                if (x(j) .ge. xout) exit 
            end do

            if (j .eq. 1) then 
                yout = y(1) 
            else if (j .eq. n+1) then 
                yout = y(n)
            else 
                yout = interp_linear_internal(x(j-1:j),y(j-1:j),xout)
            end if 
        end if 

        return 

      end function interp_linear_pt

    function interp_linear_vec(x,y,xout) result(yout)
        ! Interpolate y from ordered x to ordered xout positions

        implicit none 
 
        real(prec), dimension(:), intent(IN) :: x, y
        real(prec), dimension(:), intent(IN) :: xout
        real(prec), dimension(size(xout)) :: yout 
        integer :: i, j, n, nout 

        n    = size(x) 
        nout = size(xout)

!         write(*,*) minval(x), maxval(x), n, nout

        do i = 1, nout 
            if (xout(i) .lt. x(1)) then
                yout(i) = y(1)
!                 write(*,*) 1, xout(i)
            else if (xout(i) .gt. x(n)) then
                yout(i) = y(n)
!                 write(*,*) 2, xout(i)
            else
                do j = 1, n 
                    if (x(j) .ge. xout(i)) exit 
                end do

                if (j .eq. 1) then 
                    yout(i) = y(1) 
!                     write(*,*) 3, xout(i)
                else if (j .eq. n+1) then 
                    yout(i) = y(n)
!                     write(*,*) 4, xout(i)
                else 
                    yout(i) = interp_linear_internal(x(j-1:j),y(j-1:j),xout(i))
!                     write(*,*) 5, xout(i)
                end if 
            end if 
        end do

        return 

      end function interp_linear_vec

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   Subroutine :  interp_linear_internal
    !   Author     :  Alex Robinson
    !   Purpose    :  Interpolates for the y value at the desired x value, 
    !                 given x and y values around the desired point.
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    function interp_linear_internal(x,y,xout) result(yout)

        implicit none

        real(prec), intent(IN)  :: x(2), y(2), xout
        real(prec) :: yout
        real(prec) :: alph

        if ( xout .lt. x(1) .or. xout .gt. x(2) ) then
            write(*,*) "interp1: xout < x0 or xout > x1 !"
            write(*,*) "xout = ",xout
            write(*,*) "x0   = ",x(1)
            write(*,*) "x1   = ",x(2)
            stop
        end if

        alph = (xout - x(1)) / (x(2) - x(1))
        yout = y(1) + alph*(y(2) - y(1))

        return

    end function interp_linear_internal 

    function interp_spline(x,y,xout) result(yout)

        implicit none 
 
        real(prec), dimension(:), intent(IN) :: x, y
        real(prec), dimension(:), intent(IN) :: xout
        real(prec), dimension(size(xout)) :: yout 
        real(prec), dimension(:), allocatable :: b, c, d 
        real(prec) :: uh, dx, yh  
        integer :: i, n, nout 

        n    = size(x) 
        nout = size(xout)

        ! Get spline coefficients b, c, d
        allocate(b(n),c(n),d(n))
        call spline (x, y, b, c, d, n)

        do i = 1, nout 
            if (xout(i) .lt. x(1)) then
                dx = x(1)-xout(i)
                uh = x(1)+dx
                yh = ispline(uh,x,y,b,c,d,n)
                yout(i) = y(1) + (y(1)-yh)
                !write(*,*) x(1), xout(i), dx, uh, y(1), yh, yout(i)
            else if (xout(i) .gt. x(n)) then
                dx = xout(i)-x(n)
                uh = x(n)-dx
                yh = ispline(uh,x,y,b,c,d,n)
                yout(i) = y(n) + (y(n)-yh)
                !write(*,*) x(n), xout(i), dx, uh, y(n), yh, yout(i)
            else
                yout(i) = ispline(xout(i), x, y, b, c, d, n)
            end if 
        end do 

        return

    end function interp_spline 

subroutine spline (x, y, b, c, d, n)
!======================================================================
!  SOURCE: http://ww2.odu.edu/~agodunov/computing/programs/book2/Ch01/spline.f90
!  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
!  for cubic spline interpolation
!  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!  for  x(i) <= x <= x(i+1)
!  Alex G: January 2010
!----------------------------------------------------------------------
!  input..
!  x = the arrays of data abscissas (in strictly increasing order)
!  y = the arrays of data ordinates
!  n = size of the arrays xi() and yi() (n>=2)
!  output..
!  b, c, d  = arrays of spline coefficients
!  comments ...
!  spline.f90 program is based on fortran version of program spline.f
!  the accompanying function fspline can be used for interpolation
!======================================================================
implicit none
integer n
real(prec), dimension(:) :: x, y, b, c, d 
!real(prec) x(n), y(n), b(n), c(n), d(n)
integer i, j, gap
real(prec) h

gap = n-1
! check input
if ( n < 2 ) return
if ( n < 3 ) then
  b(1) = (y(2)-y(1))/(x(2)-x(1))   ! linear interpolation
  c(1) = 0.
  d(1) = 0.
  b(2) = b(1)
  c(2) = 0.
  d(2) = 0.
  return
end if
!
! step 1: preparation
!
d(1) = x(2) - x(1)
c(2) = (y(2) - y(1))/d(1)
do i = 2, gap
  d(i) = x(i+1) - x(i)
  b(i) = 2.0*(d(i-1) + d(i))
  c(i+1) = (y(i+1) - y(i))/d(i)
  c(i) = c(i+1) - c(i)
end do
!
! step 2: end conditions 
!
b(1) = -d(1)
b(n) = -d(n-1)
c(1) = 0.0
c(n) = 0.0
if(n /= 3) then
  c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
  c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
  c(1) = c(1)*d(1)**2/(x(4)-x(1))
  c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
end if
!
! step 3: forward elimination 
!
do i = 2, n
  h = d(i-1)/b(i-1)
  b(i) = b(i) - h*d(i-1)
  c(i) = c(i) - h*c(i-1)
end do
!
! step 4: back substitution
!
c(n) = c(n)/b(n)
do j = 1, gap
  i = n-j
  c(i) = (c(i) - d(i)*c(i+1))/b(i)
end do
!
! step 5: compute spline coefficients
!
b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
do i = 1, gap
  b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
  d(i) = (c(i+1) - c(i))/d(i)
  c(i) = 3.*c(i)
end do
c(n) = 3.0*c(n)
d(n) = d(n-1)
end subroutine spline

function ispline(u, x, y, b, c, d, n)
!======================================================================
! SOURCE: http://ww2.odu.edu/~agodunov/computing/programs/book2/Ch01/spline.f90
! function ispline evaluates the cubic spline interpolation at point z
! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
! where  x(i) <= u <= x(i+1)
!----------------------------------------------------------------------
! input..
! u       = the abscissa at which the spline is to be evaluated
! x, y    = the arrays of given data points
! b, c, d = arrays of spline coefficients computed by spline
! n       = the number of data points
! output:
! ispline = interpolated value at point u
!=======================================================================
implicit none
real(prec) ispline
integer n
! real(prec)  u, x(n), y(n), b(n), c(n), d(n)
real(prec) :: u 
real(prec), dimension(:) :: x, y, b, c, d 
integer i, j, k
real(prec) dx

! if u is ouside the x() interval take a boundary value (left or right)
if(u <= x(1)) then
  ispline = y(1)
  return
end if
if(u >= x(n)) then
  ispline = y(n)
  return
end if

!*
!  binary search for for i, such that x(i) <= u <= x(i+1)
!*
i = 1
j = n+1
do while (j > i+1)
  k = (i+j)/2
  if(u < x(k)) then
    j=k
    else
    i=k
   end if
end do
!*
!  evaluate spline interpolation
!*
dx = u - x(i)
ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
end function ispline


end module interp1D


! program test

!     use interp1D 

!     implicit none 

!     real(prec), parameter :: pi = 2.d0*acos(0.d0)

!     integer, parameter :: n    = 12
!     integer, parameter :: nout = 360
    
!     real(prec), dimension(n)    :: x, y
!     real(prec), dimension(nout) :: xout, yout

!     integer :: k, q 

!     do k = 1, n 
!         x(k) = k*30-15 
!         y(k) = sin(x(k)*2*pi/dble(nout)) + (k-6)*0.3d0
!     end do 

!     do k = 1, nout 
!         xout(k) = k 
!     end do 

!     yout = interp_spline(x,y,xout)

!     q = 1 
!     do k = 1, nout 
        
!         if (xout(k) .eq. x(q)) then 
!             write(*,"(2f10.3,2f10.3)") x(q), y(q), xout(k), yout(k) 
!             q = q+1
!         else
!             write(*,"(2a10,2f10.3)") "NA","NA", xout(k), yout(k) 
!         end if 

!     end do 

! end program test 