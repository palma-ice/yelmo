module runge_kutta

    use yelmo_defs 
    use mass_conservation, only : calc_G_advec_simple

    implicit none



    private
    public :: rk4_2D

contains


    subroutine rk4_2D(var,fvar,dvdt,ux,uy,dx,dt,solver,boundaries)

        implicit none

        real(wp), intent(INOUT) :: var(:,:) 
        real(wp), intent(INOUT) :: fvar(:,:)
        real(wp), intent(INOUT) :: dvdt(:,:) 
        real(wp), intent(IN)    :: ux(:,:) 
        real(wp), intent(IN)    :: uy(:,:) 
        real(wp), intent(IN)    :: dx
        real(wp), intent(IN)    :: dt
        character(len=*), intent(IN)    :: solver               ! Solver to use for the ice thickness advection equation
        character(len=*), intent(IN)    :: boundaries           ! Boundary conditions to impose

        ! Local variables 
        integer  :: i, j, nx, ny
        real(wp) :: dt_now 
        real(wp), allocatable :: k1(:,:) 
        real(wp), allocatable :: k2(:,:) 
        real(wp), allocatable :: k3(:,:) 
        real(wp), allocatable :: k4(:,:) 
        real(wp), allocatable :: y_now(:,:)

        real(wp), allocatable :: F(:,:) 

        nx = size(var,1) 
        ny = size(var,2) 

        allocate(k1(nx,ny))
        allocate(k2(nx,ny))
        allocate(k3(nx,ny))
        allocate(k4(nx,ny))
        allocate(y_now(nx,ny))
        
        allocate(F(nx,ny)) 

        ! Set external forcing to zero, only treating advection here
        F = 0.0_wp 

        ! ===== k1 =====

        dt_now = dt 
        y_now  = var
        
        call calc_G_advec_simple(k1,y_now,fvar,ux,uy,solver,boundaries,dx,dt_now)

        ! ===== k2 =====

        dt_now = dt / 2.0_wp 
        y_now  = var + k1*dt_now
        
        call calc_G_advec_simple(k2,y_now,fvar,ux,uy,solver,boundaries,dx,dt_now)

        ! ===== k3 =====

        dt_now = dt / 2.0_wp 
        y_now  = var + k2*dt_now
        
        call calc_G_advec_simple(k3,y_now,fvar,ux,uy,solver,boundaries,dx,dt_now)

        ! ===== k4 =====

        dt_now = dt
        y_now  = var + k3*dt_now
        
        call calc_G_advec_simple(k4,y_now,fvar,ux,uy,solver,boundaries,dx,dt_now)


        ! ===== advance to t+dt =====

        y_now = var 
        dvdt = 1.0_wp/6.0_wp * (k1 + 2.0_wp*k2 + 2.0_wp*k3 + k4)

        var  = y_now + dvdt*dt 

        return

    end subroutine rk4_2D


end module runge_kutta