module ice_tracer
    ! Module to treat online age calculations 

    use yelmo_defs, only : prec, mv 
    use solver_tridiagonal, only : solve_tridiag 
    
    implicit none
    

    private 
    public :: calc_tracer_3D
    public :: calc_tracer_column
    public :: calc_tracer_column_expl 
    public :: calc_isochrones

contains


    subroutine calc_tracer_3D(X_ice,X_srf,ux,uy,uz,H_ice,bmb,zeta_aa,zeta_ac,solver,impl_kappa,dt,dx,time)
        ! Solver for the age of the ice 
        ! Note zeta=height, k=1 base, k=nz surface 

        implicit none 

        real(prec), intent(INOUT) :: X_ice(:,:,:)   ! [units] Ice tracer variable 3D
        real(prec), intent(IN)    :: X_srf(:,:)     ! [units] Surface boundary variable field
        real(prec), intent(IN)    :: ux(:,:,:)      ! [m a-1] Horizontal x-velocity 
        real(prec), intent(IN)    :: uy(:,:,:)      ! [m a-1] Horizontal y-velocity 
        real(prec), intent(IN)    :: uz(:,:,:)      ! [m a-1] Vertical velocity 
        real(prec), intent(IN)    :: H_ice(:,:)     ! [m] Ice thickness 
        real(prec), intent(IN)    :: bmb(:,:)       ! [m a-1] Basal mass balance (melting is negative)
        real(prec), intent(IN)    :: zeta_aa(:)     ! [--] Vertical sigma coordinates (zeta==height), aa-nodes
        real(prec), intent(IN)    :: zeta_ac(:)     ! [--] Vertical sigma coordinates (zeta==height), ac-nodes
        character(len=*), intent(IN) :: solver      ! Solver choice ("impl","expl")
        real(prec), intent(IN)    :: impl_kappa     ! [m2 a-1] Artificial diffusivity
        real(prec), intent(IN)    :: dt             ! [a] Time step 
        real(prec), intent(IN)    :: dx             ! [a] Horizontal grid step 
        real(prec), intent(IN)    :: time           ! [a] Current time 

        ! Local variables
        integer    :: i, j, k, nx, ny, nz_aa, nz_ac  
        
        logical    :: is_margin 
        real(prec), allocatable :: X_base(:,:)  
        real(prec), allocatable :: advecxy(:)   ! [X a-1 m-2] Horizontal advection  
        real(prec), allocatable :: uz_now(:)    ! [m a-1] Vertical velocity 
        real(prec), allocatable :: X_prev(:,:,:) 
        real(prec) :: dz  

        ! These limits help with stability and generally would only affect areas
        ! where the tracers are not very interesting (ice shelves, ice margin, thin ice)
        real(prec), parameter :: H_ice_min = 100.0      ! [m] Minimum ice thickness to calculate ages
        real(prec), parameter :: uz_max    =   5.0      ! [m a-1] Maximum considered vertical velocity magnitude at surface 

        real(prec), parameter :: bmb_thinning = -1e-3   ! [m a-1] To keep solution bounded 

        nx    = size(X_ice,1)
        ny    = size(X_ice,2)
        nz_aa = size(zeta_aa,1)
        nz_ac = size(zeta_ac,1)

        allocate(X_base(nx,ny))
        allocate(advecxy(nz_aa))
        allocate(uz_now(nz_ac)) 
        allocate(X_prev(nx,ny,nz_aa))

        advecxy = 0.0_prec 
        uz_now  = 0.0_prec 

        ! Store X_ice from previous timestep 
        X_prev = X_ice 
        X_base = X_ice(:,:,1) 

        do j = 3, ny-2
        do i = 3, nx-2 
            
            if (H_ice(i,j) .gt. 0.0_prec .and. &
                    count([H_ice(i-1,j),H_ice(i+1,j),H_ice(i,j-1),H_ice(i,j+1)].eq.0.0_prec).gt.0) then 
                is_margin = .TRUE. 
            else
                is_margin = .FALSE. 
            end if 

            if ( H_ice(i,j) .gt. H_ice_min .and. (.not. is_margin) ) then 
                ! Thick ice exists, so call tracer solver for the column

                ! Get current column of uz 
                uz_now = uz(i,j,:)

                ! If uz presents very high values, renormalize it to keep things reasonable. 
                if (uz_now(nz_ac) .gt. uz_max)  uz_now =  uz_now/uz_now(nz_ac)*uz_max 
                if (uz_now(nz_ac) .lt. -uz_max) uz_now = -uz_now/uz_now(nz_ac)*uz_max 
                
                ! Pre-calculate the contribution of horizontal advection to column solution
                call calc_advec_horizontal_column(advecxy,X_ice,ux,uy,dx,i,j,ulim=5000.0_prec)
                
                ! Calculate the updated basal value of X from 
                ! basal mass balance including additional thinning term
                call calc_X_base(X_base(i,j),X_ice(i,j,:),H_ice(i,j),bmb(i,j),bmb_thinning,zeta_aa,dt)

                select case(trim(solver))

                    case("expl")
                        ! Explicit, upwind solver 
                        
                        call calc_tracer_column_expl(X_ice(i,j,:),uz_now,advecxy, &
                                    X_srf(i,j),X_base(i,j),H_ice(i,j),zeta_aa,zeta_ac,dt)

                    case("impl")
                        ! Implicit solver vertically, upwind horizontally 
                        
                        call calc_tracer_column(X_ice(i,j,:),uz_now,advecxy, &
                                    X_srf(i,j),X_base(i,j),H_ice(i,j),zeta_aa,zeta_ac,impl_kappa,dt) 

                    case DEFAULT 

                        write(*,*) "calc_tracer_3D:: Error: solver choice must be [expl,impl]."
                        write(*,*) "solver = ", trim(solver)
                        stop 

                end select 

            else ! H_ice(i,j) .le. H_ice_min
                ! Ice is too thin or zero, no tracing

                X_ice(i,j,:) = X_srf(i,j) 

            end if 

        end do 
        end do 

        ! Check for tracer inconsistencies
        call fix_tracer_violation(X_ice,X_prev,X_base,X_srf)

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

    subroutine fix_tracer_violation(X_ice,X_prev,X_base,X_srf)

        implicit none 

        real(prec), intent(INOUT) :: X_ice(:,:,:)       ! [units] Predicted values 
        real(prec), intent(IN)    :: X_prev(:,:,:)      ! [units] Previous values
        real(prec), intent(IN)    :: X_base(:,:)        ! [units] Basal boundary value
        real(prec), intent(IN)    :: X_srf(:,:)         ! [units] Surface boudnary value 
        
        ! Local variables
        integer :: i, j, nx, ny, nz   
        real(prec) :: X_min, X_max               
        logical    :: is_active 

        nx = size(X_ice,1)
        ny = size(X_ice,2) 
        nz = size(X_ice,3) 

        do j = 2, ny-1
        do i = 1, nx-1 

            if (minval(X_ice(i,j,:)) .eq. X_srf(i,j) .and. maxval(X_ice(i,j,:)) .eq. X_srf(i,j)) then 
                ! This column had X_srf imposed everywhere, so ignore it.
                is_active = .FALSE. 
            else 
                is_active = .TRUE. 
            end if 

            if (is_active) then 
                ! Check this column for errors:
                ! Column values cannot be outside of the range of
                ! previous values or the boundary values. 

                ! Determine minimum/maxium allowed values
                X_min = minval([X_base(i,j),X_srf(i,j),X_prev(i,j,2:nz-1), &
                    X_prev(i-1,j,2:nz-1),X_prev(i+1,j,2:nz-1),X_prev(i,j-1,2:nz-1),X_prev(i,j+1,2:nz-1)])
                X_max = maxval([X_base(i,j),X_srf(i,j),X_prev(i,j,2:nz-1), &
                    X_prev(i-1,j,2:nz-1),X_prev(i+1,j,2:nz-1),X_prev(i,j-1,2:nz-1),X_prev(i,j+1,2:nz-1)])


                if (minval(X_ice(i,j,:)) .lt. X_min .or. maxval(X_ice(i,j,:)) .gt. X_max) then 
                    ! Simply maintain solution from previous timestep, but update 
                    ! surface boundary condition
                    !X_ice(i,j,:)  = X_prev(i,j,:)
                    !X_ice(i,j,nz) = X_srf(i,j) 

                    X_ice(i,j,:)  = X_srf(i,j) 
                end if 

                ! Finally, also reset tracer to previous value if value is too large, or appears problematic 
                ! (likely redundant with previous check, but kept for safety)
                if (maxval(abs(X_ice(i,j,:))) .gt. 1e12) then 
                    X_ice(i,j,:) = X_srf(i,j) 
                    !X_ice(i,j,:) = X_prev(i,j,:)
                end if 

            end if 

        end do 
        end do 

        return 

    end subroutine fix_tracer_violation

    subroutine calc_tracer_column(X_ice,uz,advecxy,X_srf,X_base,H_ice,zeta_aa,zeta_ac,kappa,dt)
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
        real(prec), intent(IN)    :: kappa        ! [m2 a-1] Diffusivity 
        real(prec), intent(IN)    :: dt           ! [a] Time step 

        ! Local variables 
        integer :: k, nz_aa, nz_ac

        real(prec), allocatable :: dzeta_a(:)   ! nz_aa [--] Solver discretization helper variable ak
        real(prec), allocatable :: dzeta_b(:)   ! nz_aa [--] Solver discretization helper variable bk

        real(prec), allocatable :: advecz(:)   ! nz_aa, for explicit vertical advection solving
        logical, parameter      :: test_expl_advecz = .FALSE. 

        real(prec), allocatable :: subd(:)     ! nz_aa 
        real(prec), allocatable :: diag(:)     ! nz_aa  
        real(prec), allocatable :: supd(:)     ! nz_aa 
        real(prec), allocatable :: rhs(:)      ! nz_aa 
        real(prec), allocatable :: solution(:) ! nz_aa
        real(prec) :: T_base, fac, fac_a, fac_b, uz_aa, dz
        real(prec) :: kappa_a, kappa_b, dz1, dz2   
        real(prec) :: bmb_tot 

        real(prec), allocatable :: kappa_aa(:) 

        nz_aa = size(zeta_aa,1)
        nz_ac = size(zeta_ac,1)

        allocate(dzeta_a(nz_aa))
        allocate(dzeta_b(nz_aa))

        allocate(kappa_aa(nz_aa))

        ! Define dzeta terms for this column
        ! Note: for constant zeta axis, this can be done once outside
        ! instead of for each column. However, it is done here to allow
        ! use of adaptive vertical axis.
        call calc_dzeta_terms(dzeta_a,dzeta_b,zeta_aa,zeta_ac)
        
!         kappa_aa(1)     = 1e-4
!         kappa_aa(nz_aa) = 1.5e0 

!         do k = 2, nz_aa-1 
!             kappa_aa(k) = kappa_aa(1) + zeta_aa(k)*(kappa_aa(nz_aa)-kappa_aa(1))
!         end do 
        
        kappa_aa = kappa 

        allocate(subd(nz_aa))
        allocate(diag(nz_aa))
        allocate(supd(nz_aa))
        allocate(rhs(nz_aa))
        allocate(solution(nz_aa))

!         ! Calculate basal mass balance including additional thinning term
!         call calc_X_base(X_base,X_ice,H_ice,bmb,bmb_thinning,zeta_aa,dt)

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

            ! Stagger kappa to the lower and upper ac-nodes

            ! ac-node between k-1 and k 
            if (k .eq. 2) then 
                ! Bottom layer, kappa is kappa for now (later with bedrock kappa?)
                kappa_a = kappa_aa(1)
            else 
                ! Weighted average between lower half and upper half of point k-1 to k 
                dz1 = zeta_ac(k-1)-zeta_aa(k-1)
                dz2 = zeta_aa(k)-zeta_ac(k-1)
                kappa_a = (dz1*kappa_aa(k-1) + dz2*kappa_aa(k))/(dz1+dz2)
            end if 

            ! ac-node between k and k+1 

            ! Weighted average between lower half and upper half of point k to k+1
            dz1 = zeta_ac(k+1)-zeta_aa(k)
            dz2 = zeta_aa(k+1)-zeta_ac(k+1)
            kappa_b = (dz1*kappa_aa(k) + dz2*kappa_aa(k+1))/(dz1+dz2)

            ! Vertical distance for centered difference scheme
            dz      =  H_ice*(zeta_aa(k+1)-zeta_aa(k-1))

            fac_a   = -kappa_a*dzeta_a(k)*dt/H_ice**2
            fac_b   = -kappa_b*dzeta_b(k)*dt/H_ice**2

            subd(k) = fac_a - uz_aa * dt/dz
            supd(k) = fac_b + uz_aa * dt/dz
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

    subroutine calc_tracer_column_expl(X_ice,uz,advecxy,X_srf,X_base,H_ice,zeta_aa,zeta_ac,dt)
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
        real(prec), intent(IN)    :: zeta_ac(:)   ! nz_ac [--] Vertical sigma coordinates (zeta==height), layer boundary ac-nodes
        real(prec), intent(IN)    :: dt           ! [a] Time step 

        ! Local variables 
        integer :: k, nz_aa 
        real(prec), allocatable :: advecz(:)   ! nz_aa, for explicit vertical advection solving

        nz_aa = size(zeta_aa,1)
        
        allocate(advecz(nz_aa))

!         ! Calculate basal mass balance including additional thinning term
!         call calc_X_base(X_base,X_ice,H_ice,bmb,bmb_thinning,zeta_aa,dt)

        ! Update base and surface values
        X_ice(1)     = X_base 
        X_ice(nz_aa) = X_srf 

        ! Calculate vertical advection 
        call calc_advec_vertical_column_upwind1(advecz,X_ice,uz,H_ice,zeta_aa)
!         call calc_advec_vertical_column_upwind2(advecz,X_ice,uz,H_ice,zeta_aa)
!         call calc_advec_vertical_column_new2(advecz,X_ice,uz,H_ice,zeta_aa,zeta_ac,dt)

        ! Use advection terms to advance column tracer value 
        X_ice = X_ice - dt*advecxy - dt*advecz

        ! Ensure tiny values are removed to avoid underflow errors 
        where (abs(X_ice) .lt. 1e-5) X_ice = 0.0_prec 

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

                    !write(*,*) k, z1, z2, zout, dz 

                    if (zout .le. z2) then 
                        Q2   = interp_linear_pt([z1,z2],[Q(k+1),Q(k+2)],zout)

                        advecz(k) = uz(k)*((4.0*Q1-Q2-3.0*Q0))/(2.0*dz)

                    else 
                        ! 1st order 
                        advecz(k) = uz(k)*(Q(k+1)-Q(k))/dz 
                    end if 

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

    subroutine calc_advec_vertical_column_new2(advecz,Q,uz,H_ice,zeta_aa,zeta_ac,dt)
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
        real(prec), intent(IN)    :: zeta_ac(:)     ! nz_ac, cell edges
        real(prec), intent(IN)    :: dt             ! [a] Timestep 

        ! Local variables
        integer :: k, nz_aa   
        real(prec) :: u_aa, dz  
        real(prec) :: Q0, Q1, Q2 
        real(prec) :: z1, z2, zout 
        real(prec) :: Q_ac_up, Q_ac_dwn, Q_aa 
        real(prec) :: uz_ac_up, uz_ac_dwn 
        real(prec) :: Q_up, Q_dwn, Gdc, Gcu, Gc 
        real(prec) :: x_ac_dwn, x_ac_up

        real(prec) :: dx0, dx1, dQdz  
        real(prec) :: dim2, dim1, dip1, alpha1, beta1, gamma1 

        real(prec), parameter :: eps = 1e-6 

        nz_aa = size(zeta_aa,1)

        advecz = 0.0_prec 

        ! Loop over internal cell centers and perform upwind advection 
        do k = 1, nz_aa-1 
            
            if (k .ge. 2) then 
                
                dz = H_ice*(zeta_ac(k)-zeta_ac(k-1))

                uz_ac_up  = uz(k)
                uz_ac_dwn = uz(k-1) 

                Q_up     = Q(k+1)
                Q_ac_up  = 0.5_prec*(Q(k)+Q(k+1))
                Q_aa     = Q(k) 
                Q_ac_dwn = 0.5_prec*(Q(k)+Q(k-1))
                Q_dwn    = Q(k-1) 

                ! UNO2+ (Li, 2008)
                ! Gradient at midpoint 
                Gdc = (Q_dwn - Q_aa) / (H_ice*(zeta_aa(k-1)-zeta_aa(k)))
                Gcu = (Q_aa  - Q_up) / (H_ice*(zeta_aa(k)-zeta_aa(k+1)))
                Gc = sign(1.0_prec,Gdc)*2.0*abs(Gdc*Gcu)/(abs(Gdc)+abs(Gcu) + eps)

                Q_ac_dwn = Q_aa + sign(1.0_prec,Q_dwn-Q_aa)*(dz - abs(uz_ac_dwn)*dt)*abs(Gdc*Gcu)/(abs(Gdc)+abs(Gcu) + eps)
                Q_ac_up  = Q_aa + sign(1.0_prec,Q_aa -Q_up)*(abs(uz_ac_up) *dt - dz)*abs(Gdc*Gcu)/(abs(Gdc)+abs(Gcu) + eps)
                
                ! Get staggered upstream and downstream values 
                x_ac_dwn = zeta_aa(k) + sign(1.0_prec,uz_ac_dwn)*(dz-abs(uz_ac_dwn)*dt)/2.0_prec 
                x_ac_up  = zeta_aa(k) + sign(1.0_prec,uz_ac_up) *(dz-abs(uz_ac_up) *dt)/2.0_prec 
                
                advecz(k) = (uz_ac_dwn*Q_ac_dwn - uz_ac_up*Q_ac_up)/dz
!                 advecz(k) = (uz_ac_up*Q_ac_up - uz_ac_dwn*Q_ac_dwn)/dz

            else
                ! 1st order upstream
                dz        = H_ice*(zeta_aa(k+1)-zeta_aa(k))
                advecz(k) = uz(k)*(Q(k+1)-Q(k))/dz 

            end if 

!             if (k .ge. 3) then 

!                 u_aa = 0.5_prec*(uz(k-1)+uz(k))
                
!                 dim2 = H_ice*(zeta_aa(k) - zeta_aa(k-2))
!                 dim1 = H_ice*(zeta_aa(k) - zeta_aa(k-1))
!                 dip1 = H_ice*(zeta_aa(k+1) - zeta_aa(k))

!                 alpha1 =  (dim1*dip1)/(dim2*(dim2+dip1)*(dim2-dim1))
!                 beta1  = -(dim2*dip1)/(dim1*(dim2-dim1)*(dim1+dip1))
!                 gamma1 =  (dim2*dim1)/(dip1*(dim1+dip1)*(dim2+dip1))

!                 dQdz = alpha1*Q(k-2) + beta1*Q(k-1) - (alpha1+beta1+gamma1)*Q(k) + gamma1*Q(k+1)

!                 advecz(k) = u_aa*dQdz

!             else
!                 ! 1st order upstream
!                 dz        = H_ice*(zeta_aa(k+1)-zeta_aa(k))
!                 advecz(k) = uz(k)*(Q(k+1)-Q(k))/dz 

!             end if 

!             if (k .gt. 1) then 
!                 ! 2nd order 
                
!                 dx0 = H_ice*(zeta_aa(k)-zeta_aa(k-1))
!                 dx1 = H_ice*(zeta_aa(k+1)-zeta_aa(k)) 

!                 dQdz =        -dx1/(dx0*(dx1+dx0))*Q(k-1) &
!                        + (dx1-dx0)/(dx1*dx0)      *Q(k)   & 
!                              + dx0/(dx1*(dx1+dx0))*Q(k+1)

!                 advecz(k) = u_aa*dQdz 

!             else
!                 ! 1st order upstream
!                 dz        = H_ice*(zeta_aa(k+1)-zeta_aa(k))
!                 advecz(k) = uz(k)*(Q(k+1)-Q(k))/dz 

!             end if 

        end do 

        return 

    end subroutine calc_advec_vertical_column_new2

    subroutine calc_X_base(X_base,X_ice,H_ice,bmb,bmb_thinning,zeta_aa,dt)
        ! Calculate the basal boundary value of tracer 

        implicit none 

        real(prec), intent(OUT)   :: X_base
        real(prec), intent(IN)    :: X_ice(:)  
        real(prec), intent(IN)    :: H_ice
        real(prec), intent(IN)    :: bmb 
        real(prec), intent(IN)    :: bmb_thinning
        real(prec), intent(IN)    :: zeta_aa(:) 
        real(prec), intent(IN)    :: dt 

        ! Local variables 
        real(prec) :: bmb_tot 
        real(prec) :: dz 
        real(prec) :: f_wt 

        ! Calculate basal mass balance including additional thinning term
        bmb_tot = bmb + bmb_thinning 

        ! Step 0: Apply basal boundary condition
        if (bmb_tot .le. 0.0) then 
            ! Modify basal value explicitly for basal melting
            ! using boundary condition of Rybak and Huybrechts (2003), Eq. 3
            ! No horizontal advection of the basal value since it has no thickness 

            dz     = H_ice*(zeta_aa(2) - zeta_aa(1)) 
            !X_base = X_ice(1) - dt*bmb_tot*(X_ice(2)-X_ice(1))/dz 

            ! Calculate equation first as f_wt, to be able to limit it
            ! to a maximum value of 1.0 (ie, avoid extrapolation)
            f_wt   = min(1.0_prec, -dt*bmb_tot/dz)
            X_base = X_ice(1) + f_wt*(X_ice(2)-X_ice(1))

        else
            ! Leave basal value unchanged 
            X_base = X_ice(1)

        end if 

        return 

    end subroutine calc_X_base 

    subroutine calc_advec_horizontal_column(advecxy,var_ice,ux,uy,dx,i,j,ulim)
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
        real(prec), intent(IN)  :: ulim             ! [m/a] Maximum allowed velocity to apply

        ! Local variables 
        integer :: k, nx, ny, nz_aa 
        real(prec) :: ux_aa, uy_aa, u_now  
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

                ! Velocity on horizontal ac-node, vertical aa-node, limited to ulim 
                u_now = ux(i-1,j,k)
                u_now = sign(min(abs(u_now),ulim),u_now)

                ! 1st order
!                 advecx = dx_inv * u_now *(var_ice(i,j,k)-var_ice(i-1,j,k))
                ! 2nd order
                advecx = dx_inv2 * u_now *(-(4.0*var_ice(i-1,j,k)-var_ice(i-2,j,k)-3.0*var_ice(i,j,k)))

            else if (ux_aa .lt. 0.0 .and. i .le. nx-2) then 
                ! Flow to the left

                ! Velocity on horizontal ac-node, vertical aa-node, limited to ulim
                u_now = ux(i,j,k)
                u_now = sign(min(abs(u_now),ulim),u_now)

                ! 1st order 
!                 advecx = dx_inv * u_now *(var_ice(i+1,j,k)-var_ice(i,j,k))
                ! 2nd order
                advecx = dx_inv2 * u_now *((4.0*var_ice(i+1,j,k)-var_ice(i+2,j,k)-3.0*var_ice(i,j,k)))

            else 
                ! No flow 
                advecx = 0.0

            end if 

            if (uy_aa .gt. 0.0 .and. j .ge. 3) then   
                ! Flow to the right 

                ! Velocity on horizontal ac-node, vertical aa-node, limited to ulim
                u_now = uy(i,j-1,k)
                u_now = sign(min(abs(u_now),ulim),u_now)

                ! 1st order
!                 advecy = dx_inv * u_now*(var_ice(i,j,k)-var_ice(i,j-1,k))
                ! 2nd order
                advecy = dx_inv2 * u_now *(-(4.0*var_ice(i,j-1,k)-var_ice(i,j-2,k)-3.0*var_ice(i,j,k)))

            else if (uy_aa .lt. 0.0 .and. j .le. ny-2) then 
                ! Flow to the left

                ! Velocity on horizontal ac-node, vertical aa-node, limited to ulim
                u_now = uy(i,j,k)
                u_now = sign(min(abs(u_now),ulim),u_now)

                ! 1st order 
!                 advecy = dx_inv * u_now *(var_ice(i,j+1,k)-var_ice(i,j,k))
                ! 2nd order
                advecy = dx_inv2 * u_now *((4.0*var_ice(i,j+1,k)-var_ice(i,j+2,k)-3.0*var_ice(i,j,k)))

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

    subroutine calc_dzeta_terms(dzeta_a,dzeta_b,zeta_aa,zeta_ac)
        ! zeta_aa  = depth axis at layer centers (plus base and surface values)
        ! zeta_ac  = depth axis (1: base, nz: surface), at layer boundaries
        ! Calculate ak, bk terms as defined in Hoffmann et al (2018)
        implicit none 

        real(prec), intent(INOUT) :: dzeta_a(:)    ! nz_aa
        real(prec), intent(INOUT) :: dzeta_b(:)    ! nz_aa
        real(prec), intent(IN)    :: zeta_aa(:)    ! nz_aa 
        real(prec), intent(IN)    :: zeta_ac(:)    ! nz_ac == nz_aa-1 

        ! Local variables 
        integer :: k, nz_layers, nz_aa    

        nz_aa = size(zeta_aa)

        ! Note: zeta_aa is calculated outside in the main program 

        ! Initialize dzeta_a/dzeta_b to zero, first and last indices will not be used (end points)
        dzeta_a = 0.0 
        dzeta_b = 0.0 
        
        do k = 2, nz_aa-1 
            dzeta_a(k) = 1.0/ ( (zeta_ac(k) - zeta_ac(k-1)) * (zeta_aa(k) - zeta_aa(k-1)) )
        enddo

        do k = 2, nz_aa-1
            dzeta_b(k) = 1.0/ ( (zeta_ac(k) - zeta_ac(k-1)) * (zeta_aa(k+1) - zeta_aa(k)) )
        end do

        return 

    end subroutine calc_dzeta_terms


    subroutine calc_isochrones(depth_iso,dep_time,H_ice,age_iso,zeta,time)
        ! Calculate specific isochronal layers (depth_iso) at specified ages (age_iso)
        ! from info about the deposition time of each ice layer (dep_time)

        ! Note: to maintain generality, calculate everything based on dep_time, rather
        ! than age (ie, convert the age of each layer to a deposition time) 

        implicit none 

        real(prec), intent(OUT) :: depth_iso(:,:,:)     ! [m]   Depth of isochronal layer
        real(prec), intent(IN)  :: dep_time(:,:,:)      ! [a]   Deposition time of ice layer
        real(prec), intent(IN)  :: H_ice(:,:)           ! [m]   Ice thickness 
        real(prec), intent(IN)  :: age_iso(:)           ! [ka]  Isochronal layer ages
        real(prec), intent(IN)  :: zeta(:)              ! [-]   Vertical sigma coordinates
        real(prec), intent(IN)  :: time                 ! [a]   Current time 

        ! Local variables 
        integer :: i, j, k, k0, k1, q, nx, ny, nz, n_iso 
        real(prec) :: dep_time_iso, zeta_iso  

        nx    = size(H_ice,1)
        ny    = size(H_ice,2) 
        nz    = size(zeta) 
        n_iso = size(age_iso) 

        ! Initially set depth of isochrones to zero 
        depth_iso = 0.0 

        ! Calculate each isochrone at each location separately 
                
        do q = 1, n_iso 

            ! Get isochronal age in terms of deposition time [ka] => [a]
            dep_time_iso = 0.0 - (age_iso(q)*1e3) 

            if (dep_time_iso .lt. time) then 
                ! This layer should exist, perform interpolation at each point 

                do j = 1, ny 
                do i = 1, nx 

                    if (H_ice(i,j) .gt. 0.0) then 
                        ! Ice exists here, get isochronal depth 

                        k1 = mv 

                        do k = 1, nz
                            if (dep_time(i,j,k) .gt. dep_time_iso) then 
                                k1 = k 
                                exit 
                            end if 
                        end do 

                        if (k1 .eq. 1) then 
                            ! All ice is younger than current isochronal deposition time
                            depth_iso(i,j,q) = (1.0-zeta(1))*H_ice(i,j) 

                        else if (k1 .eq. mv) then 
                            ! All ice is older than current isochronal deposition time
                            depth_iso(i,j,q) = 0.0 

                        else 
                            ! Isochronal deposition time is inside column 
                            k0 = k1 - 1 
                            zeta_iso = interp_linear_pt(dep_time(i,j,k0:k1),zeta(k0:k1),xout=dep_time_iso)
                            depth_iso(i,j,q) = (1.0-zeta_iso)*H_ice(i,j) 

                        end if 

                    end if 

                end do 
                end do 

            end if 

        end do 
              
        return 

    end subroutine calc_isochrones

end module ice_tracer 

