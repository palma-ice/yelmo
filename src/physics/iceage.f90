module iceage
    ! Module to treat online age calculations 

    use yelmo_defs, only : prec, pi, g, sec_year, T0, rho_ice, rho_sw, rho_w, L_ice  
    use solver_tridiagonal, only : solve_tridiag 
    
    implicit none
    

    private 
    public :: calc_tracer_3D

contains


    subroutine calc_tracer_3D(X_ice,ux,uy,uz,H_ice,bmb,zeta_aa,zeta_ac,dzeta_a,dzeta_b,solver,impl_kappa,dt,dx,time)
        ! Solver for the age of the ice 
        ! Note zeta=height, k=1 base, k=nz surface 

        implicit none 

        real(prec), intent(INOUT) :: X_ice(:,:,:)   ! [units] Ice tracer variable 3D
        real(prec), intent(IN)    :: ux(:,:,:)      ! [m a-1] Horizontal x-velocity 
        real(prec), intent(IN)    :: uy(:,:,:)      ! [m a-1] Horizontal y-velocity 
        real(prec), intent(IN)    :: uz(:,:,:)      ! [m a-1] Vertical velocity 
        real(prec), intent(IN)    :: H_ice(:,:)     ! [m] Ice thickness 
        real(prec), intent(IN)    :: bmb(:,:)       ! [m a-1] Basal mass balance (melting is negative)
        real(prec), intent(IN)    :: zeta_aa(:)     ! [--] Vertical sigma coordinates (zeta==height), aa-nodes
        real(prec), intent(IN)    :: zeta_ac(:)     ! [--] Vertical sigma coordinates (zeta==height), ac-nodes
        real(prec), intent(IN)    :: dzeta_a(:)     ! d Vertical height axis (0:1) 
        real(prec), intent(IN)    :: dzeta_b(:)     ! d Vertical height axis (0:1)
        character(len=*), intent(IN) :: solver      ! Solver choice ("impl","expl")
        real(prec), intent(IN)    :: impl_kappa     ! [m2 a-1] Artificial diffusivity
        real(prec), intent(IN)    :: dt             ! [a] Time step 
        real(prec), intent(IN)    :: dx             ! [a] Horizontal grid step 
        real(prec), intent(IN)    :: time           ! [a] Current time 

        ! Local variables
        integer :: i, j, k, nx, ny, nz_aa, nz_ac  
        real(prec) :: X_srf, X_base 
        real(prec), allocatable  :: advecxy(:)   ! [X a-1 m-2] Horizontal advection  
        real(prec) :: dz  

        ! Testing Rybak and Huybrechts (2003) analytical case at dome
        logical,    parameter    :: use_rh2003_profile = .FALSE. 
        real(prec), parameter    :: smb_test = 0.1     ! [m a-1] 
        real(prec), allocatable  :: uz_test(:)

        nx    = size(X_ice,1)
        ny    = size(X_ice,2)
        nz_aa = size(zeta_aa,1)
        nz_ac = size(zeta_ac,1)

        allocate(advecxy(nz_aa))
        advecxy = 0.0 

        allocate(uz_test(nz_aa))
        uz_test = zeta_aa*(-smb_test) 

        do j = 3, ny-2
        do i = 3, nx-2 
            
            if (H_ice(i,j) .gt. 10.0) then 
                ! Thick ice exists, call thermodynamic solver for the column

                ! Set surface value to current time 
                X_srf = time 

                ! Apply basal boundary condition
                if (bmb(i,j) .lt. 0.0) then 
                    ! Modify basal value explicitly for basal melting
                    ! using boundary condition of Rybak and Huybrechts (2003), Eq. 3
                    ! No horizontal advection of the basal value since it has no thickness 

                    dz = H_ice(i,j)*(zeta_aa(2) - zeta_aa(1)) 
                    X_base = X_ice(i,j,1) - dt*bmb(i,j)*(X_ice(i,j,2)-X_ice(i,j,1))/dz 

                else
                    ! Leave basal value unchanged 
                    X_base = X_ice(i,j,1)

                end if 

                ! Pre-calculate the contribution of horizontal advection to column solution
                call calc_advec_horizontal_column(advecxy,X_ice,ux,uy,dx,i,j)
                
                select case(trim(solver))

                    case("expl")
                        ! Explicit, upwind solver 
                        
                        if (use_rh2003_profile) then 

                            ! Only calculate for summit
                            if (i .eq. 15 .and. j .eq. 15) then 
                                call calc_tracer_column_expl(X_ice(i,j,:),uz_test,advecxy*0.0,X_srf,X_base,H_ice(i,j),zeta_aa,dt)
                            end if 

                        else 
                            call calc_tracer_column_expl(X_ice(i,j,:),uz(i,j,:),advecxy,X_srf,X_base,H_ice(i,j),zeta_aa,dt)

                        end if 

                    case("impl")
                        ! Implicit solver vertically, upwind horizontally 
                        
                        if (use_rh2003_profile) then 

                            ! Only calculate for summit
                            if (i .eq. 15 .and. j .eq. 15) then 
                                call calc_tracer_column(X_ice(i,j,:),uz_test,advecxy*0.0,X_srf,X_base,H_ice(i,j),zeta_aa,zeta_ac, &
                                                                                                    dzeta_a,dzeta_b,impl_kappa,dt)
                            end if 

                        else 
                            call calc_tracer_column(X_ice(i,j,:),uz(i,j,:),advecxy,X_srf,X_base,H_ice(i,j),zeta_aa,zeta_ac, &
                                                                                                dzeta_a,dzeta_b,impl_kappa,dt)
                    
                        end if 

                    case DEFAULT 

                        write(*,*) "calc_tracer_3D:: Error: solver choice must be [expl,impl]."
                        write(*,*) "solver = ", trim(solver)
                        stop 

                end select 

            else ! H_ice(i,j) .le. 10.0
                ! Ice is too thin or zero, no tracing

                X_ice(i,j,:) = time 

            end if 

        end do 
        end do 

        ! Fill in borders 
        X_ice(2,:,:)    = X_ice(3,:,:) 
        X_ice(1,:,:)    = X_ice(3,:,:) 
        X_ice(nx-1,:,:) = X_ice(nx-2,:,:) 
        X_ice(nx,:,:)   = X_ice(nx-2,:,:) 
        
        X_ice(:,2,:)    = X_ice(:,3,:) 
        X_ice(:,1,:)    = X_ice(:,3,:) 
        X_ice(:,ny-1,:) = X_ice(:,ny-2,:) 
        X_ice(:,ny,:)   = X_ice(:,ny-2,:) 
        
        return 

    end subroutine calc_tracer_3D

    subroutine calc_tracer_column(X_ice,uz,advecxy,X_srf,X_base,H_ice,zeta_aa,zeta_ac,dzeta_a,dzeta_b,kappa,dt)
        ! Tracer solver for a given column of ice 
        ! Note zeta=height, k=1 base, k=nz surface 
        ! Note: nz = number of vertical boundaries (including zeta=0.0 and zeta=1.0), 
        ! temperature is defined for cell centers, plus a value at the surface and the base
        ! so nz_ac = nz_aa - 1 

        ! For notes on implicit form of advection terms, see eg http://farside.ph.utexas.edu/teaching/329/lectures/node90.html
        
        implicit none 

        real(prec), intent(INOUT) :: X_ice(:)     ! nz_aa [units] Ice column tracer variable
        real(prec), intent(IN)    :: uz(:)        ! nz_ac [m a-1] Vertical velocity 
        real(prec), intent(IN)    :: advecxy(:)   ! nz_aa [K a-1] Horizontal heat advection 
        real(prec), intent(IN)    :: X_srf        ! [units] Surface value
        real(prec), intent(IN)    :: X_base       ! [units] Basal value
        real(prec), intent(IN)    :: H_ice        ! [m] Ice thickness 
        real(prec), intent(IN)    :: zeta_aa(:)   ! nz_aa [--] Vertical sigma coordinates (zeta==height), layer centered aa-nodes
        real(prec), intent(IN)    :: zeta_ac(:)   ! nz_ac [--] Vertical height axis temperature (0:1), layer edges ac-nodes
        real(prec), intent(IN)    :: dzeta_a(:)   ! nz_aa [--] Solver discretization helper variable ak
        real(prec), intent(IN)    :: dzeta_b(:)   ! nz_aa [--] Solver discretization helper variable bk
        real(prec), intent(IN)    :: kappa        ! [m2 a-1] Diffusivity 
        real(prec), intent(IN)    :: dt           ! [a] Time step 

        ! Local variables 
        integer :: k, nz_aa, nz_ac

        real(prec), allocatable :: advecz(:)   ! nz_aa, for explicit vertical advection solving
        logical, parameter      :: test_expl_advecz = .FALSE. 

        real(prec), allocatable :: subd(:)     ! nz_aa 
        real(prec), allocatable :: diag(:)     ! nz_aa  
        real(prec), allocatable :: supd(:)     ! nz_aa 
        real(prec), allocatable :: rhs(:)      ! nz_aa 
        real(prec), allocatable :: solution(:) ! nz_aa
        real(prec) :: T_base, fac, fac_a, fac_b, uz_aa, dz  

        nz_aa = size(zeta_aa,1)
        nz_ac = size(zeta_ac,1)

        allocate(subd(nz_aa))
        allocate(diag(nz_aa))
        allocate(supd(nz_aa))
        allocate(rhs(nz_aa))
        allocate(solution(nz_aa))

        ! Step 1: apply vertical advection (for explicit testing)
        if (test_expl_advecz) then 
            allocate(advecz(nz_aa))
            advecz = 0.0
            call calc_advec_vertical_column_upwind1(advecz,X_ice,uz,H_ice,zeta_aa)
            X_ice = X_ice - dt*advecz 
        end if 

        ! Step 2: apply vertical implicit advection 
        
        ! Ice base
        subd(1) = 0.0_prec
        diag(1) = 1.0_prec
        supd(1) = 0.0_prec
        rhs(1)  = X_base 
        
        ! Ice interior layers 2:nz_aa-1
        do k = 2, nz_aa-1

            if (test_expl_advecz) then 
                ! No implicit vertical advection (diffusion only)
                uz_aa = 0.0 

            else
                ! With implicit vertical advection (diffusion + advection)
                uz_aa   = 0.5*(uz(k-1)+uz(k))   ! ac => aa nodes
            end if 

            dz      =  H_ice*(zeta_ac(k)-zeta_ac(k-1))
            
            fac     = dt * kappa / H_ice**2
            fac_a   = -fac*dzeta_a(k)
            fac_b   = -fac*dzeta_b(k)
            subd(k) = fac_a - uz_aa*dt / (2.0*dz)
            supd(k) = fac_b + uz_aa*dt / (2.0*dz)
            diag(k) = 1.0_prec - fac_a - fac_b
            rhs(k)  = X_ice(k) - dt*advecxy(k) 

        end do 

        ! Ice surface 
        subd(nz_aa) = 0.0_prec
        diag(nz_aa) = 1.0_prec
        supd(nz_aa) = 0.0_prec
        rhs(nz_aa)  = X_srf

        ! Call solver 
        call solve_tridiag(subd,diag,supd,rhs,solution)

        ! Copy the solution into the tracer variable
        X_ice  = solution

        return 

    end subroutine calc_tracer_column

    subroutine calc_tracer_column_expl(X_ice,uz,advecxy,X_srf,X_base,H_ice,zeta_aa,dt)
        ! Tracer solver for a given column of ice 
        ! Note zeta=height, k=1 base, k=nz surface 
        ! Note: nz = number of vertical boundaries (including zeta=0.0 and zeta=1.0), 
        ! temperature is defined for cell centers, plus a value at the surface and the base
        ! so nz_ac = nz_aa - 1 

        ! For notes on implicit form of advection terms, see eg http://farside.ph.utexas.edu/teaching/329/lectures/node90.html
        
        implicit none 

        real(prec), intent(INOUT) :: X_ice(:)     ! nz_aa [units] Ice column tracer variable
        real(prec), intent(IN)    :: uz(:)        ! nz_ac [m a-1] Vertical velocity 
        real(prec), intent(IN)    :: advecxy(:)   ! nz_aa [K a-1] Horizontal heat advection 
        real(prec), intent(IN)    :: X_srf        ! [units] Surface value
        real(prec), intent(IN)    :: X_base       ! [units] Basal value
        real(prec), intent(IN)    :: H_ice        ! [m] Ice thickness 
        real(prec), intent(IN)    :: zeta_aa(:)   ! nz_aa [--] Vertical sigma coordinates (zeta==height), layer centered aa-nodes
        real(prec), intent(IN)    :: dt           ! [a] Time step 

        ! Local variables 
        integer :: k, nz_aa
        real(prec), allocatable :: advecz(:)   ! nz_aa, for explicit vertical advection solving

        nz_aa = size(zeta_aa,1)
        
        allocate(advecz(nz_aa))

        ! Update base and surface values
        X_ice(1)     = X_base 
        X_ice(nz_aa) = X_srf 

        ! Calculate vertical advection 
!         call calc_advec_vertical_column_upwind1(advecz,X_ice,uz,H_ice,zeta_aa)
        call calc_advec_vertical_column_upwind2(advecz,X_ice,uz,H_ice,zeta_aa)

        ! Use advection terms to advance column tracer value 
        X_ice = X_ice - dt*advecz - dt*advecxy 

        return 

    end subroutine calc_tracer_column_expl

    subroutine calc_advec_vertical_column_upwind1(advecz,Q,uz,H_ice,zeta_aa)
        ! Calculate vertical advection term advecz, which enters
        ! advection equation as
        ! Q_new = Q - dt*advecz = Q - dt*u*dQ/dx
        ! 1st order upwind scheme 

        implicit none 

        real(prec), intent(OUT)   :: advecz(:)      ! nz_aa: bottom, cell centers, top 
        real(prec), intent(INOUT) :: Q(:)           ! nz_aa: bottom, cell centers, top 
        real(prec), intent(IN)    :: uz(:)          ! nz_ac: cell boundaries
        real(prec), intent(IN)    :: H_ice          ! Ice thickness 
        real(prec), intent(IN)    :: zeta_aa(:)     ! nz_aa, cell centers
        
        ! Local variables
        integer :: k, nz_aa   
        real(prec) :: u_aa, dz  

        nz_aa = size(zeta_aa,1)

        advecz = 0.0 

        ! Loop over internal cell centers and perform upwind advection 
        do k = 2, nz_aa-1 
            
            u_aa = 0.5_prec*(uz(k-1)+uz(k))
            
            if (u_aa < 0.0) then 
                ! Upwind negative
                dz        = H_ice*(zeta_aa(k+1)-zeta_aa(k))
                advecz(k) = uz(k)*(Q(k+1)-Q(k))/dz  

            else
                ! Upwind positive
                dz        = H_ice*(zeta_aa(k)-zeta_aa(k-1))
                advecz(k) = uz(k-1)*(Q(k)-Q(k-1))/dz

            end if 
        end do 

        return 

    end subroutine calc_advec_vertical_column_upwind1

    subroutine calc_advec_vertical_column_upwind2(advecz,Q,uz,H_ice,zeta_aa)
        ! Calculate vertical advection term advecz, which enters
        ! advection equation as
        ! Q_new = Q - dt*advecz = Q - dt*u*dQ/dx
        ! 2nd order upwind scheme 

        implicit none 

        real(prec), intent(OUT)   :: advecz(:)      ! nz_aa: bottom, cell centers, top 
        real(prec), intent(INOUT) :: Q(:)           ! nz_aa: bottom, cell centers, top 
        real(prec), intent(IN)    :: uz(:)          ! nz_ac: cell boundaries
        real(prec), intent(IN)    :: H_ice          ! Ice thickness 
        real(prec), intent(IN)    :: zeta_aa(:)     ! nz_aa, cell centers
!         real(prec), intent(IN)    :: zeta_ac(:)     ! nz_ac, cell edges
        
        ! Local variables
        integer :: k, nz_aa   
        real(prec) :: u_aa, dz  
        real(prec) :: Q0, Q1, Q2 
        real(prec) :: z1, z2, zout 

        nz_aa = size(zeta_aa,1)

        advecz = 0.0 

        ! Loop over internal cell centers and perform upwind advection 
        do k = 2, nz_aa-1 
            
            u_aa = 0.5_prec*(uz(k-1)+uz(k))
            
            if (u_aa < 0.0) then 
                ! Upwind negative
                dz = H_ice*(zeta_aa(k+1)-zeta_aa(k))

                if (k .le. nz_aa-2) then 
                    ! 2nd order 
                    Q0   = Q(k)
                    Q1   = Q(k+1)
                    z1   = H_ice*zeta_aa(k+1)
                    z2   = H_ice*zeta_aa(k+2)
                    zout = H_ice*zeta_aa(k+1)+dz 
                    Q2   = interp_linear_pt([z1,z2],[Q(k+1),Q(k+2)],zout)

                    advecz(k) = uz(k)*((4.0*Q1-Q2-3.0*Q0))/(2.0*dz)
                else
                    ! 1st order 
                    advecz(k) = uz(k)*(Q(k+1)-Q(k))/dz 

                end if 

            else
                ! Upwind positive
                dz = H_ice*(zeta_aa(k)-zeta_aa(k-1))

                if (k .ge. 3) then 
                    ! 2nd order 
                    Q0   = Q(k)
                    Q1   = Q(k-1)
                    z1   = H_ice*zeta_aa(k-1)
                    z2   = H_ice*zeta_aa(k-2)
                    zout = H_ice*zeta_aa(k-1)-dz
                    if (zout .lt. z2) then 
                        write(*,*) "calc_advec_vertical_column_upwind2:: Error: upwind2 interp step doesnt work for this spacing."
                        stop 
                    end if 
                    Q2   = interp_linear_pt([z2,z1],[Q(k-2),Q(k-1)],zout)

                    advecz(k) = uz(k-1)*(-(4.0*Q1-Q2-3.0*Q0))/(2.0*dz)
                else
                    ! 1st order 
                    advecz(k) = uz(k-1)*(Q(k)-Q(k-1))/dz
                
                end if 
                
            end if 
        end do 

        return 

    end subroutine calc_advec_vertical_column_upwind2

    subroutine calc_advec_horizontal_column(advecxy,var_ice,ux,uy,dx,i,j)
        ! Newly implemented advection algorithms (ajr)
        ! Output: [K a-1]

        ! [m-1] * [m a-1] * [K] = [K a-1]

        implicit none

        real(prec), intent(OUT) :: advecxy(:)       ! nz_aa 
        real(prec), intent(IN)  :: var_ice(:,:,:)   ! nx,ny,nz_aa  Enth, T, age, etc...
        real(prec), intent(IN)  :: ux(:,:,:)        ! nx,ny,nz
        real(prec), intent(IN)  :: uy(:,:,:)        ! nx,ny,nz
        real(prec), intent(IN)  :: dx  
        integer,    intent(IN)  :: i, j 

        ! Local variables 
        integer :: k, nx, ny, nz_aa 
        real(prec) :: ux_aa, uy_aa 
        real(prec) :: dx_inv, dx_inv2
        real(prec) :: advecx, advecy 

        ! Define some constants 
        dx_inv  = 1.0_prec / dx 
        dx_inv2 = 1.0_prec / (2.0_prec*dx)

        nx  = size(var_ice,1)
        ny  = size(var_ice,2)
        nz_aa = size(var_ice,3) 

        advecx = 0.0 
        advecy = 0.0 

        ! Loop over each point in the column
        do k = 2, nz_aa-1 

            ! Estimate direction of current flow into cell (x and y), centered in vertical layer and grid point
            ux_aa = 0.25_prec*(ux(i,j,k)+ux(i-1,j,k)+ux(i,j,k-1)+ux(i-1,j,k-1))
            uy_aa = 0.25_prec*(uy(i,j,k)+uy(i,j-1,k)+uy(i,j,k-1)+uy(i,j-1,k-1))

            ! Explicit form (to test different order approximations)
            if (ux_aa .gt. 0.0 .and. i .ge. 3) then  
                ! Flow to the right 

                ! 1st order
!                 advecx = dx_inv * 0.5*(ux(i-1,j,k)+ux(i-1,j,k-1))*(var_ice(i,j,k)-var_ice(i-1,j,k))
                ! 2nd order
                !advecx = dx_inv2 * 0.5*(ux(i-1,j,k)+ux(i-1,j,k-1))*(-(4.0*var_ice(i-1,j,k)-var_ice(i-2,j,k)-3.0*var_ice(i,j,k)))
                ! ux/uy on zeta_aa nodes:
                advecx = dx_inv2 * ux(i-1,j,k)*(-(4.0*var_ice(i-1,j,k)-var_ice(i-2,j,k)-3.0*var_ice(i,j,k)))
                
            else if (ux_aa .lt. 0.0 .and. i .le. nx-2) then 
                ! Flow to the left

                ! 1st order 
!                 advecx = dx_inv * 0.5*(ux(i,j,k)+ux(i,j,k-1))*(var_ice(i+1,j,k)-var_ice(i,j,k))
                ! 2nd order
!                 advecx = dx_inv2 * 0.5*(ux(i,j,k)+ux(i,j,k-1))*((4.0*var_ice(i+1,j,k)-var_ice(i+2,j,k)-3.0*var_ice(i,j,k)))
                ! ux/uy on zeta_aa nodes:
                advecx = dx_inv2 * ux(i,j,k)*((4.0*var_ice(i+1,j,k)-var_ice(i+2,j,k)-3.0*var_ice(i,j,k)))
                
            else 
                ! No flow 
                advecx = 0.0

            end if 

            if (uy_aa .gt. 0.0 .and. j .ge. 3) then   
                ! Flow to the right 

                ! 1st order
!                 advecy = dx_inv * 0.5*(uy(i,j-1,k)+uy(i,j-1,k-1))*(var_ice(i,j,k)-var_ice(i,j-1,k))
                ! 2nd order
!                 advecy = dx_inv2 * 0.5*(uy(i,j-1,k)+uy(i,j-1,k-1))*(-(4.0*var_ice(i,j-1,k)-var_ice(i,j-2,k)-3.0*var_ice(i,j,k)))
                ! ux/uy on zeta_aa nodes:
                advecy = dx_inv2 * uy(i,j-1,k)*(-(4.0*var_ice(i,j-1,k)-var_ice(i,j-2,k)-3.0*var_ice(i,j,k)))
                
            else if (uy_aa .lt. 0.0 .and. j .le. ny-2) then 
                ! Flow to the left

                ! 1st order 
!                 advecy = dx_inv * 0.5*(uy(i,j,k)+uy(i,j,k-1))*(var_ice(i,j+1,k)-var_ice(i,j,k))
                ! 2nd order
!                 advecy = dx_inv2 * 0.5*(uy(i,j,k)+uy(i,j,k-1))*((4.0*var_ice(i,j+1,k)-var_ice(i,j+2,k)-3.0*var_ice(i,j,k)))
                ! ux/uy on zeta_aa nodes:
                advecy = dx_inv2 * uy(i,j,k)*((4.0*var_ice(i,j+1,k)-var_ice(i,j+2,k)-3.0*var_ice(i,j,k)))
                
            else
                ! No flow 
                advecy = 0.0 

            end if 
                    
            ! Combine advection terms for total contribution 
            advecxy(k) = (advecx+advecy)

        end do 

        return 

    end subroutine calc_advec_horizontal_column
    
    function interp_linear_pt(x,y,xout) result(yout)
        ! Interpolates for the y value at the desired x value, 
        ! given x and y values around the desired point.
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

    end function interp_linear_pt 

end module iceage 

