module ice_benchmarks 
    
    ! This module implements the verification 
    ! tests A-E outlined in Bueler et al. (2005), 
    ! and tests F-G described in Bueler et al. (2007)

    ! It also implements the EISMINT1 and EISMINT2 boundary conditions 
    
    use yelmo_defs, only : prec, pi  

    implicit none


    type bueler_test_type 

        real(prec), allocatable :: H_ice(:,:) 
        real(prec), allocatable :: mbal(:,:) 
        real(prec), allocatable :: u_b(:,:) 
        real(prec), allocatable :: mbal_c(:,:)   ! [m a-1] Compensatory mass balance
        real(prec), allocatable :: Q_c(:,:)      ! [K a-1] Compensatory heating 
        
        ! Comparison values
        real(prec), allocatable :: err_H_ice(:,:)
        real(prec) :: err_H0
        real(prec) :: err_max_H_ice
        real(prec) :: rmse_H_ice 

        real(prec) :: V_ice_mod 
        real(prec) :: V_ice_target
        real(prec) :: err_V_ice 
        
    end type 

    private 
    public :: bueler_test_type
    public :: bueler_init
    public :: bueler_compare
    public :: bueler_test_AE 
    public :: bueler_test_BC

    public :: eismint_boundaries 

contains

    subroutine bueler_init(buel,nx,ny)

        implicit none 

        type(bueler_test_type), intent(INOUT) :: buel 
        integer, intent(IN) :: nx
        integer, intent(IN) :: ny 

        allocate(buel%H_ice(nx,ny))
        allocate(buel%mbal(nx,ny))
        allocate(buel%u_b(nx,ny))
        allocate(buel%mbal_c(nx,ny))
        allocate(buel%Q_c(nx,ny))

        allocate(buel%err_H_ice(nx,ny))

        buel%H_ice      = 0.0 
        buel%mbal       = 0.0 
        buel%u_b        = 0.0 
        buel%mbal_c     = 0.0 
        buel%Q_c        = 0.0 
        buel%err_H_ice  = 0.0 

        buel%err_H0         = 0.0
        buel%err_max_H_ice  = 0.0  
        buel%rmse_H_ice     = 0.0 
        buel%err_V_ice      = 0.0 

        return 

    end subroutine bueler_init 

    subroutine bueler_compare(buel,H_ice,dx)

        implicit none 

        type(bueler_test_type), intent(INOUT) :: buel
        real(prec), intent(IN) :: H_ice(:,:)            ! [m]
        real(prec), intent(IN) :: dx                    ! [m]

        ! Local variables 
        logical, allocatable :: msk(:,:) 

        allocate(msk(size(H_ice,1),size(H_ice,2)))

        msk = .FALSE. 
        where(H_ice .gt. 0.0 .or. buel%H_ice .gt. 0.0) msk = .TRUE. 

        buel%err_H_ice     = H_ice - buel%H_ice 
        buel%err_H0        = maxval(H_ice) - maxval(buel%H_ice)
        buel%err_max_H_ice = maxval(abs(buel%err_H_ice))

        if (count(msk) .gt. 0) then 
            buel%rmse_H_ice    = sqrt(sum(buel%err_H_ice**2,mask=msk) / real(count(msk)))
        else 
            buel%rmse_H_ice    = 0.0 
        end if 
        
        buel%V_ice_mod     = sum(dx*dx*H_ice)                           *1e-9*1e-6   ! [m^3] => [1e6 km^3] 
        buel%V_ice_target  = sum(dx*dx*buel%H_ice)                      *1e-9*1e-6   ! [m^3] => [1e6 km^3]
        buel%err_V_ice     = (buel%V_ice_mod - buel%V_ice_target)       *1e3         ! [1e6 km^3] => [1e3 km^3] 

        return 

    end subroutine bueler_compare 
    
    subroutine bueler_test_AE(H_ice,mbal,u_b,xx,yy,L,mbal0,A,n,rho_ice,g,mu_max)

        implicit none 

        real(prec), intent(OUT) :: H_ice(:,:) 
        real(prec), intent(OUT) :: mbal(:,:) 
        real(prec), intent(OUT) :: u_b(:,:) 
        real(prec), intent(IN)  :: xx(:,:)      ! [m]
        real(prec), intent(IN)  :: yy(:,:)      ! [m]
        real(prec), intent(IN)  :: A 
        real(prec), intent(IN)  :: L 
        real(prec), intent(IN)  :: mbal0 
        real(prec), intent(IN)  :: n
        real(prec), intent(IN)  :: rho_ice 
        real(prec), intent(IN)  :: g 
        real(prec), intent(IN)  :: mu_max  
         
        ! Local variables 
        integer :: nx, ny 
        real(prec), allocatable :: r(:,:) 
        real(prec), allocatable :: gamma(:,:) 
        real(prec)  :: L_meters

        nx = size(H_ice,1)
        ny = size(H_ice,2) 

        allocate(r(nx,ny))
        allocate(gamma(nx,ny))

        ! Calculate the radius value as a function of xx and yy [m]
        r = sqrt(xx**2 + yy**2)

        ! Get parameter L in meters [km] => [m]
        L_meters = L * 1e3 

        ! Calculate gamma 
        gamma = bueler_gamma(A,n,rho_ice,g)

        ! First calculate the Bodvarsson-Vialov profile (Eq. 17 in Bueler et al, 2005)
        H_ice = 0.0 
        where (r .le. L_meters) &
            H_ice = (2.0**(n-1)*mbal0/gamma)**(1.0/(2.0*n+2.0))*(L_meters**(1.0+1.0/n)-r**(1.0+1.0/n))**(n/(2.0*n+2.0))

        ! Now calculate implied mass balance (constant)
        mbal  = mbal0 

        ! Set the basal velocity to zero 
        u_b = 0.0 

        return 

    end subroutine bueler_test_AE 

    subroutine bueler_test_BC(H_ice,mbal,u_b,xx,yy,time,R0,H0,lambda,n,A,rho_ice,g)

        implicit none 

        real(prec), intent(OUT) :: H_ice(:,:) 
        real(prec), intent(OUT) :: mbal(:,:) 
        real(prec), intent(OUT) :: u_b(:,:) 
        real(prec), intent(IN)  :: xx(:,:)      ! [m] 
        real(prec), intent(IN)  :: yy(:,:)      ! [m] 
        real(prec), intent(IN)  :: time         ! [a] Time relative to t0 
        real(prec), intent(IN)  :: R0 
        real(prec), intent(IN)  :: H0 
        real(prec), intent(IN)  :: lambda  
        real(prec), intent(IN)  :: n
        real(prec), intent(IN)  :: A            ! [Pa3 a m-1] 
        real(prec), intent(IN)  :: rho_ice 
        real(prec), intent(IN)  :: g   
         
        ! Local variables 
        integer    :: i, j, nx, ny 
        real(prec) :: r_now  
        real(prec) :: R0_meters
        real(prec) :: alpha, beta, gamma, t0, time1  
        
        real(prec) :: xx1(11), yy1(11), H1(11,11)

        real(prec), parameter :: f = 0.0        ! isostasy fraction 

        nx = size(H_ice,1)
        ny = size(H_ice,2) 

        ! Get parameter R0 in meters [km] => [m]
        R0_meters = R0 * 1e3 

        ! Calculate alpha, beta, t0 and absolute time
        alpha = (2.0 - (n+1.0)*lambda)/(5.0*n+3.0)
        beta  = (1.0 + (2.0*n+1.0)*lambda) / (5.0*n+3.0)
        gamma = bueler_gamma(A,n,rho_ice,g)
        t0    = (beta/gamma) * ((2.0*n+1.0)/(n+1.0))**n * (R0_meters**(n+1)/H0**(2.0*n+1.0))
        time1 = time + t0 

        H_ice = 0.0_prec 

        do j = 2, ny-1
        do i = 2, nx-1 

            ! Calculate the radius value as a function of xx and yy [m]
            r_now = sqrt(xx(i,j)**2 + yy(i,j)**2)

            ! Calculate the Halfar similarity solution profile (Eq. 10-11 in Bueler et al, 2005)
            if ( (1.0 - (((time1/t0)**(-beta))*r_now/R0_meters)**((n+1.0)/n)) .gt. 0.0_prec) &
            H_ice(i,j) = H0 * (time1/t0)**(-alpha) * (1.0 - (((time1/t0)**(-beta))*r_now/R0_meters)**((n+1.0)/n))**(n/(2.0*n+1.0))

            ! Now calculate implied mass balance
            mbal(i,j)  = (lambda/time1)*H_ice(i,j)  

        end do 
        end do 

        ! Set the basal velocity to zero everywhere
        u_b = 0.0 

        return 

    end subroutine bueler_test_BC 

    elemental function bueler_gamma(A,n,rho_ice,g) result(gamma)
        ! Default gamma = 9.0177e-13 m-3 s-1 

        implicit none 

        real(prec), intent(IN) :: A 
        real(prec), intent(IN) :: n 
        real(prec), intent(IN) :: rho_ice
        real(prec), intent(IN) :: g 
        real(prec) :: gamma 

        gamma = 2.0_prec * A * (rho_ice*g)**n / (n+2.0_prec)

        return 

    end function bueler_gamma 


    subroutine eismint_boundaries(T_srf,smb,ghf,xx,yy,H_ice,experiment,time,period,rad_el,dT_test,dsmb_test)

        implicit none 

        real(prec), intent(OUT) :: T_srf(:,:) 
        real(prec), intent(OUT) :: smb(:,:) 
        real(prec), intent(OUT) :: ghf(:,:) 
        real(prec), intent(IN)  :: xx(:,:)      ! [m] 
        real(prec), intent(IN)  :: yy(:,:)      ! [m] 
        real(prec), intent(IN)  :: H_ice(:,:) 
        character(len=*), intent(IN) :: experiment 
        real(prec), intent(IN) :: time 
        real(prec), intent(IN) :: period 
        real(prec), intent(IN), optional :: rad_el 
        real(prec), intent(IN), optional :: dT_test 
        real(prec), intent(IN), optional :: dsmb_test  

        real(prec), parameter :: x_summit = 0.0 
        real(prec), parameter :: y_summit = 0.0 
        
        ! Local variables 
        integer    :: i, j, nx, ny 
        real(prec) :: dist
        real(prec) :: R_el, s  
        real(prec) :: dT, dR_el, dsmb 
        real(prec) :: Tmin 

        nx = size(T_srf,1)
        ny = size(T_srf,2)
        
        if (period .gt. 0.0) then
            ! Transient forcing  
            dT    =  10.0*sin(2.0*pi*time/period)
            dR_el = 100.0*sin(2.0*pi*time/period)
            dsmb  =   0.2*sin(2.0*pi*time/period)
        
        else 
            ! Constant forcing 
            dT    = 0.0 
            dR_el = 0.0
            dsmb  = 0.0 

        end if 

        ! Allow modification of surface temperature via argument 
        if (present(dT_test)) dT = dT + dT_test

        ! Allow modification of surface temperature via argument 
        if (present(dsmb_test)) dsmb = dsmb + dsmb_test


        select case(trim(experiment))

            case("fixed")

                ! Surface temperature 

                do i = 1, nx
                do j = 1, ny 

                    dist = max(abs(xx(i,j)-x_summit),abs(yy(i,j)-y_summit)) *1e-3  ! [km]
                    T_srf(i,j) = 239.0 + 8e-8*dist**3 + dT 

                end do 
                end do 

                ! Surface mass balance 

                smb = 0.3 + dsmb   ! [m/a]

                ! Geothermal heat flux 

                ghf = 42.0   ! [mW/m2]

            case("moving")

                ! Surface temperature 

                T_srf = 270.0 - 0.01*H_ice + dT     ! [K]

                ! Surface mass balance 

                ! Default EISMINT1 parameter values 
                R_el = 450.0 + dR_el ! [km]
                s    = 0.01  ! [m/a / km]

                ! Allow modification of equilibrium line radius via argument 
                if (present(rad_el)) R_el = rad_el + dR_el

                do j = 1, ny 
                do i = 1, nx
                
                    dist = sqrt((xx(i,j)-x_summit)**2 + (yy(i,j)-y_summit)**2) *1e-3  ! [km]
                    smb(i,j) = min(0.5 + dsmb,s*(R_el-dist))

                end do 
                end do 

                ! Geothermal heat flux 

                ghf = 42.0   ! [mW/m2]
            
            case ("mismip")

                ! Surface temperature 

                T_srf = 270.0 - 0.01*H_ice + dT   ! [K]

!                 Surface mass balance 

!                 ! To allow floating ice to grow
!                 R_el = 1000.0 + dR_el ! [km]
!                 s    = 0.01           ! [m/a / km]
  
!                 do j = 1, ny 
!                 do i = 1, nx
                
!                     dist = sqrt((xx(i,j)-x_summit)**2 + (yy(i,j)-y_summit)**2) *1e-3  ! [km]
!                     smb(i,j) = min(0.5,s*(R_el-dist))

!                 end do 
!                 end do 
                
                smb = 0.3   ! [m/a]
                
                ! Geothermal heat flux 

                ghf = 42.0   ! [mW/m2]
            
            case("EXPA","EXPF")

                ! Surface temperature 

                Tmin  = 238.15                          ! [K] 
                if (trim(experiment) .eq. "EXPF") Tmin = 223.15
                
                s     = 1.67e-2                         ! [K km-1]
                
                do j = 1, ny 
                do i = 1, nx
                
                    dist = sqrt((xx(i,j)-x_summit)**2 + (yy(i,j)-y_summit)**2) *1e-3  ! [km]
                    T_srf(i,j) = Tmin + s*dist + dT     ! [K]

                end do 
                end do 

                ! Surface mass balance 

                ! Default EISMINT1 parameter values 
                R_el = 450.0 + dR_el ! [km]
                s    = 0.01  ! [m/a / km]

                ! Allow modification of equilibrium line radius via argument 
                if (present(rad_el)) R_el = rad_el + dR_el

                do j = 1, ny 
                do i = 1, nx
                
                    dist = sqrt((xx(i,j)-x_summit)**2 + (yy(i,j)-y_summit)**2) *1e-3  ! [km]
                    smb(i,j) = min(0.5,s*(R_el-dist))

                end do 
                end do 

                ! Geothermal heat flux 

                ghf = 42.0   ! [mW/m2]
            
            case DEFAULT 

                write(*,*) "Experiment not recognized: "//trim(experiment)
                stop 

        end select 

        return 

    end subroutine eismint_boundaries  
    
end module ice_benchmarks  
