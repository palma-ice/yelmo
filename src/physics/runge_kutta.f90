module runge_kutta

    use yelmo_defs 
    use mass_conservation, only : calc_G_advec_simple

    implicit none

    ! ajr: moved to yelmo_defs...
    ! type rk4_class
    !     real(wp), allocatable  :: tau(:,:) 
    !     real(wp), allocatable  :: y_np1(:,:) 
    !     real(wp), allocatable  :: y_n(:,:) 
    !     real(wp), allocatable  :: y_nm1(:,:) 
    !     real(wp), allocatable  :: y_nm2(:,:) 
    !     real(wp), allocatable  :: f_np1(:,:) 
    !     real(wp), allocatable  :: f_n(:,:) 
    !     real(wp), allocatable  :: f_nm1(:,:) 
    !     real(wp), allocatable  :: f_nm2(:,:) 
    !     real(wp) :: dt 
    !     real(wp) :: dt_nm1
    !     real(wp) :: dt_nm2
    ! end type

    private
    !public :: rk4_2D_class
    public :: rk23_2D_step
    public :: rk4_2D_step
    public :: rk4_2D_calc_truncation_error
    public :: rk4_2D_init 

contains
    
    subroutine rk23_2D_step(rk4,var,fvar,dvdt,ux,uy,dx,dt,solver,boundaries,F)
        ! Adapted from https://fncbook.github.io/fnc/ivp/adaptive.html#function-rk23 
        ! Implementation of the Bogacki–Shampine method, described in detail here:
        ! https://en.wikipedia.org/wiki/Bogacki%E2%80%93Shampine_method

        implicit none

        type(rk4_class), intent(INOUT) :: rk4 
        real(wp), intent(INOUT) :: var(:,:) 
        real(wp), intent(INOUT) :: fvar(:,:)
        real(wp), intent(INOUT) :: dvdt(:,:) 
        real(wp), intent(IN)    :: ux(:,:) 
        real(wp), intent(IN)    :: uy(:,:) 
        real(wp), intent(IN)    :: dx
        real(wp), intent(IN)    :: dt
        character(len=*), intent(IN)    :: solver               ! Solver to use for the ice thickness advection equation
        character(len=*), intent(IN)    :: boundaries           ! Boundary conditions to impose
        real(wp), intent(IN), optional :: F(:,:) 

        ! Local variables 
        integer  :: i, j, nx, ny
        real(wp) :: dt_now 
        real(wp), allocatable :: k1(:,:) 
        real(wp), allocatable :: k2(:,:) 
        real(wp), allocatable :: k3(:,:) 
        real(wp), allocatable :: k4(:,:) 
        real(wp), allocatable :: y_now(:,:)
        real(wp), allocatable :: y_new_2(:,:)
        real(wp), allocatable :: y_new_3(:,:)

        real(wp), allocatable :: F_now(:,:) 

        nx = size(var,1) 
        ny = size(var,2) 

        allocate(k1(nx,ny))
        allocate(k2(nx,ny))
        allocate(k3(nx,ny))
        allocate(k4(nx,ny))
        allocate(y_now(nx,ny))
        allocate(y_new_2(nx,ny))
        allocate(y_new_3(nx,ny))

        allocate(F_now(nx,ny)) 

        ! Set external forcing to zero, only treating advection here
        F_now = 0.0_wp 

        ! If external forcing provided, use it.
        if (present(F)) F_now = F 

        ! ===== k1 =====

        dt_now = dt 
        y_now  = var
        
        call calc_G_advec_simple(k1,y_now,fvar,ux,uy,solver,boundaries,dx,dt_now)

        ! ===== k2 =====

        dt_now = dt / 2.0_wp 
        y_now  = var + k1*dt_now
        
        call calc_G_advec_simple(k2,y_now,fvar,ux,uy,solver,boundaries,dx,dt_now)

        ! ===== k3 =====

        dt_now = dt * (3.0_wp / 4.0_wp)
        y_now  = var + k2*dt_now
        
        call calc_G_advec_simple(k3,y_now,fvar,ux,uy,solver,boundaries,dx,dt_now)

        ! 2nd order solution
        y_new_2 = var + dt*(2.0_wp*k1 + 3.0_wp*k2 + 4.0_wp*k3)/9.0_wp

        ! ===== k4 =====

        dt_now = dt
        y_now  = y_new_2
        
        call calc_G_advec_simple(k4,y_now,fvar,ux,uy,solver,boundaries,dx,dt_now)

        ! 3rd order solution
        y_new_3 = var + dt*( (7.0_wp/24.0_wp)*k1 + (1.0_wp/4.0_wp)*k2 + (1.0_wp/3.0_wp)*k3 + (1.0_wp/8.0_wp)*k4 )

        ! 2nd/3rd difference
        rk4%tau = dt*(-5.0_wp*k1/72.0_wp + k2/12.0_wp + k3/9.0_wp - k4/8.0_wp)
        !rk4%tau = (-5.0_wp*k1/72.0_wp + k2/12.0_wp + k3/9.0_wp - k4/8.0_wp)

        ! Store 2nd or 3rd order solution as current solution
        !var = y_new_2
        var = y_new_3

        ! Notes on error
        ! E = norm(err,Inf)                         # error estimate
        ! maxerr = tol*(1 + norm(u[i],Inf))     # relative/absolute blend

        ! # Accept the proposed step?
        ! if E < maxerr     # yes
        !     push!(t,t[i]+h)
        !     push!(u,unew2)
        !     i += 1
        !     k1 = k4       # use FSAL property
        ! end

        ! # Adjust step size.
        ! q = 0.8*(maxerr/E)^(1/3)   # conservative optimal step factor
        ! q = min(q,4)               # limit stepsize growth

        return

    end subroutine rk23_2D_step

    subroutine rk4_2D_step(rk4,var,fvar,dvdt,ux,uy,dx,dt,solver,boundaries)

        implicit none

        type(rk4_class), intent(INOUT) :: rk4 
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


        ! ==========

        ! Store new information in rk4 object 
        rk4%y_nm2 = rk4%y_nm1 
        rk4%y_nm1 = rk4%y_n 
        rk4%y_n   = rk4%y_np1 
        rk4%y_np1 = var 

        rk4%f_nm2 = rk4%f_nm1 
        rk4%f_nm1 = rk4%f_n 
        rk4%f_n   = rk4%f_np1 
        rk4%f_np1 = dvdt 

        rk4%dt_nm2 = rk4%dt_nm1 
        rk4%dt_nm1 = rk4%dt 
        rk4%dt     = dt 

        ! Update truncation error 
        call rk4_2D_calc_truncation_error(rk4%tau,rk4%y_np1,rk4%y_n,rk4%y_nm1,rk4%y_nm2, &
                       rk4%f_np1,rk4%f_n,rk4%f_nm1,rk4%f_nm2,rk4%dt,rk4%dt_nm1,rk4%dt_nm2)
        !rk4%tau = 1e-1_wp 

        ! ==========

        return

    end subroutine rk4_2D_step

    subroutine rk4_2D_calc_truncation_error(tau,y_np1,y_n,y_nm1,y_nm2, &
                                f_np1,f_n,f_nm1,f_nm2,dt,dt_nm1,dt_nm2)

        implicit none

        real(wp), intent(OUT) :: tau(:,:) 
        real(wp), intent(IN)  :: y_np1(:,:) 
        real(wp), intent(IN)  :: y_n(:,:) 
        real(wp), intent(IN)  :: y_nm1(:,:) 
        real(wp), intent(IN)  :: y_nm2(:,:) 
        real(wp), intent(IN)  :: f_np1(:,:) 
        real(wp), intent(IN)  :: f_n(:,:) 
        real(wp), intent(IN)  :: f_nm1(:,:) 
        real(wp), intent(IN)  :: f_nm2(:,:) 
        real(wp), intent(IN)  :: dt 
        real(wp), intent(IN)  :: dt_nm1
        real(wp), intent(IN)  :: dt_nm2
         
        
        ! Local variables 
        real(wp) :: h 
        real(wp) :: d1 
        real(wp) :: d2 
        real(wp), parameter :: d1_nudge = 1.1e-5_wp
        real(wp), parameter :: d2_nudge = 1.0e-5_wp

        h  = dt 
        d1 = dt_nm1 / h
        d2 = dt_nm2 / h

        if (d1 .eq. 0.0_wp) d1 = d1_nudge
        if (d2 .eq. 0.0_wp) d2 = d2_nudge
        if (d1 .eq. d2) d2 = d2+d2_nudge

        ! Equation 4.4, simplified for numerical testing 

        tau =    &
        - ( (1.0_wp/(1.0_wp+d2) + 1.0_wp/(1.0_wp+d1) + 1.0_wp) * y_np1 &
           + ( ((1.0_wp+d2)/d2)**2 * ((1.0_wp+d1)/d1)**2 * (1.0_wp/d2 + 1.0_wp/d1 - 1.0_wp) ) * y_n &
           + ( ((1.0_wp+d2)/(d2-d1))**2 * (1.0_wp/d1)**2 * &
                     (1.0_wp/(d2-d1) - 1.0_wp/d1 - 1.0_wp/(1.0_wp+d1)) ) * y_nm1   &
           + ( ((1.0_wp+d1)/(d1-d2))**2 * (1.0_wp/d2)**2 * &
                     (1.0_wp/(d1-d2) - 1.0_wp/d2 - 1.0_wp/(1.0_wp+d2)) ) * y_nm2 ) &
        + (h/2.0_wp) * &
          ( f_np1 &
          + (((1.0_wp+d2)/d2)**2 * ((1.0_wp+d1)/d1)**2) * f_n   &
          + (((1.0_wp+d2)/(d2-d1))**2 * (1.0_wp/d1)**2) * f_nm1 & 
          + (((1.0_wp+d1)/(d1-d2))**2 * (1.0_wp/d2)**2) * f_nm2 )

        return

    end subroutine rk4_2D_calc_truncation_error


    subroutine rk4_2D_init(rk4,nx,ny,dt_min)

        implicit none

        type(rk4_class), intent(INOUT) :: rk4
        integer,  intent(IN) :: nx 
        integer,  intent(IN) :: ny
        real(wp), intent(IN) :: dt_min 

        ! Allocate rk4 arrays 
        call rk4_2D_alloc(rk4,nx,ny)

        ! Populate with initial values 
        rk4%tau     = 1e-8_wp
        rk4%y_np1   = 0.0_wp 
        rk4%y_n     = rk4%y_np1
        rk4%y_nm1   = rk4%y_np1
        rk4%y_nm2   = rk4%y_np1
        rk4%f_np1   = 0.0_wp 
        rk4%f_n     = rk4%f_np1
        rk4%f_nm1   = rk4%f_np1
        rk4%f_nm2   = rk4%f_np1

        rk4%dt      = dt_min 
        rk4%dt_nm1  = dt_min 
        rk4%dt_nm2  = dt_min 
        
        return

    end subroutine rk4_2D_init

    subroutine rk4_2D_alloc(rk4,nx,ny)

        implicit none

        type(rk4_class), intent(INOUT) :: rk4
        integer, intent(IN) :: nx 
        integer, intent(IN) :: ny

        call rk4_2D_dealloc(rk4)
        
        allocate(rk4%tau(nx,ny))
        allocate(rk4%y_np1(nx,ny))
        allocate(rk4%y_n(nx,ny))
        allocate(rk4%y_nm1(nx,ny))
        allocate(rk4%y_nm2(nx,ny))
        allocate(rk4%f_np1(nx,ny))
        allocate(rk4%f_n(nx,ny))
        allocate(rk4%f_nm1(nx,ny))
        allocate(rk4%f_nm2(nx,ny))
        
        return

    end subroutine rk4_2D_alloc

    subroutine rk4_2D_dealloc(rk4)

        implicit none

        type(rk4_class), intent(INOUT) :: rk4

        if (allocated(rk4%tau))     deallocate(rk4%tau)
        if (allocated(rk4%y_np1))   deallocate(rk4%y_np1)
        if (allocated(rk4%y_n))     deallocate(rk4%y_n)
        if (allocated(rk4%y_nm1))   deallocate(rk4%y_nm1)
        if (allocated(rk4%y_nm2))   deallocate(rk4%y_nm2)
        if (allocated(rk4%f_np1))   deallocate(rk4%f_np1)
        if (allocated(rk4%f_n))     deallocate(rk4%f_n)
        if (allocated(rk4%f_nm1))   deallocate(rk4%f_nm1)
        if (allocated(rk4%f_nm2))   deallocate(rk4%f_nm2)
        
        return

    end subroutine rk4_2D_dealloc
    
end module runge_kutta