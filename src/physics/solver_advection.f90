module solver_advection
    
    use yelmo_defs, only : sp, dp, wp, tol_underflow, io_unit_err
    use yelmo_tools, only : get_neighbor_indices, stagger_aa_ab, stagger_nodes_aa_ab_ice
    
    use solver_linear
    use solver_advection_sico, only : calc_adv2D_expl_sico, calc_adv2D_impl_sico

    implicit none 
    
    private 
    public :: calc_advec2D
    public :: calc_adv2D_expl
    public :: calc_adv2D_impl_upwind

contains 

    subroutine calc_advec2D(dvdt,var,f_ice,ux,uy,var_dot,mask_adv,dx,dy,dt,solver,boundaries)
        ! General routine to apply 2D advection equation to variable `var` 
        ! with source term `var_dot`. Various solvers are possible

        real(wp),       intent(OUT)   :: dvdt(:,:)              ! [dvdt] Variable rate of change
        real(wp),       intent(IN)    :: var(:,:)               ! [var]  Variable to be advected
        real(wp),       intent(IN)    :: f_ice(:,:)             ! [var]  Variable to be advected
        real(wp),       intent(IN)    :: ux(:,:)                ! [m/a] 2D velocity, x-direction (ac-nodes)
        real(wp),       intent(IN)    :: uy(:,:)                ! [m/a] 2D velocity, y-direction (ac-nodes)
        real(wp),       intent(IN)    :: var_dot(:,:)           ! [dvar/dt] Source term for variable
        integer,        intent(IN)    :: mask_adv(:,:)          ! Advection mask
        real(wp),       intent(IN)    :: dx                     ! [m]   Horizontal resolution, x-direction
        real(wp),       intent(IN)    :: dy                     ! [m]   Horizontal resolution, y-direction
        real(wp),       intent(IN)    :: dt                     ! [a]   Timestep 
        character(len=*), intent(IN)    :: solver               ! Solver to use for the ice thickness advection equation
        character(len=*), intent(IN)    :: boundaries           ! Boundary conditions to impose

        ! Local variables
        integer :: nx, ny  
        type(linear_solver_class) :: lgs
        character(len=256)        :: adv_lis_opt 

        real(wp), allocatable :: var_now(:,:) 

        nx = size(var,1)
        ny = size(var,2)

        allocate(var_now(nx,ny))

        ! Assign local variable to be modified 
        var_now = var 

        select case(trim(solver))
            ! Choose solver to use 

            case("none") 

                ! Do nothing: no advection considered. 

            case("expl")

                call calc_adv2D_expl(var_now,f_ice,ux,uy,var_dot,dx,dy,dt,boundaries)

            case("expl-new")

                call calc_adv2D_expl_new(var_now,f_ice,ux,uy,var_dot,dx,dy,dt,boundaries)

            case("expl-upwind")

                call calc_adv2D_expl_upwind(var_now,ux,uy,var_dot,dx,dy,dt,boundaries)

            case("impl-upwind") 

                call calc_adv2D_impl_upwind(var_now,ux,uy,var_dot,dx,dy,dt,boundaries,f_upwind=1.0_wp)
            
            ! Other solvers below...
            case("expl-sico")

                call calc_adv2D_expl_sico(var_now,ux,uy,var_dot,dx,dy,dt)

            case("impl-sico")

                call calc_adv2D_impl_sico(var_now,ux,uy,var_dot,dx,dy,dt,use_lis=.FALSE.)

            case("impl-sico-lis")

                call calc_adv2D_impl_sico(var_now,ux,uy,var_dot,dx,dy,dt,use_lis=.TRUE.)

            case("impl-lis")

                ! Initialize linear solver variables
                call linear_solver_init(lgs,nx,ny,nvar=1,n_terms=5)

                ! Populate advection matrices Ax=b
                call linear_solver_matrix_advection_csr_2D(lgs,var_now,ux,uy,var_dot,mask_adv,dx,dy,dt,boundaries)
                
                ! Solve linear equation
                adv_lis_opt = "-i bicg -p ilu -maxiter 1000 -tol 1.0e-12 -initx_zeros false"
                !adv_lis_opt = "-i minres -p jacobi -maxiter 1000 -tol 1.0e-12 -initx_zeros false"
                call linear_solver_matrix_solve(lgs,adv_lis_opt)
                
                !call linear_solver_print_summary(lgs,io_unit_err)

                ! Store advection solution
                call linear_solver_save_advection(var_now,lgs)

            case DEFAULT 

                write(*,*) "calc_advec2D:: Error: solver not recognized."
                write(*,*) "solver = ", trim(solver)
                stop 

        end select 
        
        ! Determine rate of change 
        dvdt = (var_now - var) / dt 

        return 

    end subroutine calc_advec2D

    subroutine linear_solver_save_advection(H,lgs)
        ! Extract ice thickness solution from lgs object. 

        implicit none 

        real(wp), intent(OUT) :: H(:,:)                 ! [X]
        type(linear_solver_class), intent(IN) :: lgs

        ! Local variables 
        integer :: i, j, nr 

        do nr = 1, lgs%nmax

            i = lgs%n2i(nr)
            j = lgs%n2j(nr)

            H(i,j) = lgs%x_value(nr)

        end do

        return

    end subroutine linear_solver_save_advection

    subroutine linear_solver_matrix_advection_csr_2D(lgs,H,ux,uy,F,mask,dx,dy,dt,boundaries)
        ! Define sparse matrices A*x=b in format 'compressed sparse row' (csr)
        ! for 2D advection equations with velocity components
        ! ux and uy defined on ac-nodes (right and top borders of i,j grid cell)
        ! and variable to be advected H defined on aa nodes.
        ! Store sparse matrices in linear_solver_class object 'lgs' for later use.

        implicit none 

        type(linear_solver_class), intent(INOUT) :: lgs
        real(wp), intent(INOUT)   :: H(:,:)         ! [X] Variable of interest (aa nodes)
        real(wp), intent(IN)      :: ux(:,:)        ! [m a-1] Horizontal velocity x-direction (ac nodes)
        real(wp), intent(IN)      :: uy(:,:)        ! [m a-1] Horizontal velocity y-direction (ac nodes)
        real(wp), intent(IN)      :: F(:,:)         ! [m a-1] Net source/sink terms (aa nodes)
        integer,  intent(IN)      :: mask(:,:)      ! Advection mask
        real(wp), intent(IN)      :: dx             ! [m] Horizontal step x-direction
        real(wp), intent(IN)      :: dy             ! [m] Horizontal step y-direction 
        real(wp), intent(IN)      :: dt             ! [a] Time step 
        character(len=*), intent(IN) :: boundaries 

        ! Local variables  
        integer  :: i, j, k, nx, ny
        integer  :: im1, ip1, jm1, jp1 
        integer  :: n, nr, nc
        real(wp) :: dt_darea
        character(len=56) :: bcs(4)

        real(wp), allocatable  :: ux_1(:,:), ux_2(:,:)
        real(wp), allocatable  :: uy_1(:,:), uy_2(:,:)
        real(wp), allocatable  :: Hx_1(:,:), Hx_2(:,:)
        real(wp), allocatable  :: Hy_1(:,:), Hy_2(:,:)

        real(wp), parameter :: WOVI = 1.0     ! Weighting parameter for the over-implicit scheme 

        nx = size(H,1)
        ny = size(H,2) 

        dt_darea = dt/(dx*dy)

        ! Boundary conditions (bcs) counterclockwise unit circle 
        ! 1: x, right-border
        ! 2: y, upper-border 
        ! 3: x, left--border 
        ! 4: y, lower-border 
        
        ! Define border conditions (only choices are: no-slip, free-slip, periodic)
        select case(trim(boundaries)) 

            case("MISMIP3D","TROUGH")

                bcs(1) = "zero"
                bcs(2) = "periodic"
                bcs(3) = "infinite"
                bcs(4) = "periodic" 

            case("infinite")

                bcs(1:4) = "infinite" 
            
            case("zeros")

                bcs(1:4) = "zero" 

            case("periodic","periodic-xy")

                bcs(1:4) = "periodic" 

            case DEFAULT 

                bcs(1:4) = "zero"

        end select 


        ! Safety check for initialization
        if (.not. allocated(lgs%x_value)) then 
            ! Object 'lgs' has not been initialized yet, do so now.

            call linear_solver_init(lgs,nx,ny,nvar=1,n_terms=5)

        end if

        ! Allocate local variables
        allocate(ux_1(nx,ny))
        allocate(ux_2(nx,ny))
        allocate(uy_1(nx,ny))
        allocate(uy_2(nx,ny))
        
        allocate(Hx_1(nx,ny))
        allocate(Hx_2(nx,ny))
        allocate(Hy_1(nx,ny))
        allocate(Hy_2(nx,ny))
        
        ! ====================================================================================
        ! Step 1: populate variables representing velocity components (ux,uy) at each
        ! cell face (left,right,bottom,top) and upstream variable (Hx,Hy) at each cell face. 

        ux_1 = 0.0
        ux_2 = 0.0
        uy_1 = 0.0
        uy_2 = 0.0
        
        Hx_1 = 0.0
        Hx_2 = 0.0
        Hy_1 = 0.0
        Hy_2 = 0.0

        do j = 1, ny
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

            ! Get velocity components on each cell face 
            ux_1(i,j) = ux(im1,j)
            ux_2(i,j) = ux(i,j)
            uy_1(i,j) = uy(i,jm1)
            uy_2(i,j) = uy(i,j)

            if (ux_1(i,j) >= 0.0) then
                Hx_1(i,j) = H(im1,j)
            else
                Hx_1(i,j) = H(i,j)
            end if

            if (ux_2(i,j) >= 0.0) then
                Hx_2(i,j) = H(i,j)
            else
                Hx_2(i,j) = H(ip1,j)
            end if

            if (uy_1(i,j) >= 0.0) then
                Hy_1(i,j) = H(i,jm1)
            else
                Hy_1(i,j) = H(i,j)
            end if

            if (uy_2(i,j) >= 0.0) then
                Hy_2(i,j) = H(i,j)
            else
                Hy_2(i,j) = H(i,jp1)
            end if

        end do
        end do

        !-------- Assembly of the system of linear equations
        !             (matrix storage: compressed sparse row CSR) --------

        ! Initialize values to zero 
        lgs%a_index = 0.0
        lgs%a_value = 0.0 
        lgs%b_value = 0.0 
        lgs%x_value = 0.0 

        lgs%a_ptr(1) = 1

        k = 0

        do n = 1, lgs%nmax

            ! Get i,j indices of current point
            i = lgs%n2i(n)
            j = lgs%n2j(n)

            nr = n   ! row counter

            ! Get neighbor indices assuming periodic domain
            ! (all other boundary conditions are treated with special cases below)
            im1 = i-1
            if (im1 .eq. 0)    im1 = nx 
            ip1 = i+1
            if (ip1 .eq. nx+1) ip1 = 1 

            jm1 = j-1
            if (jm1 .eq. 0)    jm1 = ny
            jp1 = j+1
            if (jp1 .eq. ny+1) jp1 = 1 
            

            ! Handle special cases first, otherwise populate with normal inner discretization

            if (mask(i,j) .eq. 0) then 
                ! Zero thickness imposed 

                k = k+1
                lgs%a_index(k) = nr
                lgs%a_value(k) = 1.0_wp   ! diagonal element only
                
                lgs%b_value(nr) = 0.0_wp
                lgs%x_value(nr) = 0.0_wp

            else if (mask(i,j) .eq. -1) then 
                ! Prescribed ice thickness imposed 

                k = k+1
                lgs%a_index(k) = nr
                lgs%a_value(k) = 1.0_wp   ! diagonal element only
                
                lgs%b_value(nr) = 0.0_wp
                lgs%x_value(nr) = H(i,j)
            
            else if ( (.not. trim(bcs(1)) .eq. "periodic") .and. i .eq. nx) then
                ! Right border

                if (bcs(1) .eq. "infinite") then

                    k = k+1
                    lgs%a_index(k) = lgs%ij2n(i,j)                  ! column counter for H(i,j)
                    lgs%a_value(k) =  1.0_wp
                    
                    k = k+1
                    lgs%a_index(k) = lgs%ij2n(im1,j)                ! column counter for H(im1,j)
                    lgs%a_value(k) = -1.0_wp

                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = H(i,j)

                else
                    ! Assume zero for now

                    k = k+1
                    lgs%a_index(k) = nr
                    lgs%a_value(k) = 1.0_wp   ! diagonal element only
                    
                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = 0.0_wp

                end if

            else if ( (.not. trim(bcs(2)) .eq. "periodic") .and. j .eq. ny) then
                ! Top border

                if (bcs(2) .eq. "infinite") then

                    k = k+1
                    lgs%a_index(k) = lgs%ij2n(i,j)                  ! column counter for H(i,j)
                    lgs%a_value(k) =  1.0_wp
                    
                    k = k+1
                    lgs%a_index(k) = lgs%ij2n(i,jm1)                ! column counter for H(i,jm1)
                    lgs%a_value(k) = -1.0_wp

                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = H(i,j)

                else
                    ! Assume zero for now

                    k = k+1
                    lgs%a_index(k) = nr
                    lgs%a_value(k) = 1.0_wp   ! diagonal element only
                    
                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = 0.0_wp

                end if

            else if ( (.not. trim(bcs(3)) .eq. "periodic") .and. i .eq. 1) then
                ! Left border

                if (bcs(3) .eq. "infinite") then

                    k = k+1
                    lgs%a_index(k) = lgs%ij2n(i,j)                  ! column counter for H(i,j)
                    lgs%a_value(k) =  1.0_wp
                    
                    k = k+1
                    lgs%a_index(k) = lgs%ij2n(ip1,j)                ! column counter for H(ip1,j)
                    lgs%a_value(k) = -1.0_wp

                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = H(i,j)

                    ! Initial guess == previous H

                    lgs%x_value(nr) = H(i,j)                            
                

                else
                    ! Assume zero for now

                    k = k+1
                    lgs%a_index(k) = nr
                    lgs%a_value(k) = 1.0_wp   ! diagonal element only
                    
                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = 0.0_wp

                end if

            else if ( (.not. trim(bcs(4)) .eq. "periodic") .and. j .eq. 1) then
                ! Bottom border

                if (bcs(4) .eq. "infinite") then

                    k = k+1
                    lgs%a_index(k) = lgs%ij2n(i,j)                  ! column counter for H(i,j)
                    lgs%a_value(k) =  1.0_wp
                    
                    k = k+1
                    lgs%a_index(k) = lgs%ij2n(i,jp1)                ! column counter for H(i,jp1)
                    lgs%a_value(k) = -1.0_wp

                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = H(i,j)

                else
                    ! Assume zero for now

                    k = k+1
                    lgs%a_index(k) = nr
                    lgs%a_value(k) = 1.0_wp   ! diagonal element only
                    
                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = 0.0_wp

                end if

            else
                ! Inner point

                k = k+1
                lgs%a_index(k) = lgs%ij2n(i,jm1)                    ! for H(i,jm1)
                if (uy_1(i,j) > 0.0) &
                    lgs%a_value(k) = -dt_darea*uy_1(i,j)*dx*WOVI

                k = k+1
                lgs%a_index(k) = lgs%ij2n(im1,j)                    ! for H(im1,j)
                if (ux_1(i,j) > 0.0) &
                    lgs%a_value(k) = -dt_darea*ux_1(i,j)*dy*WOVI

                k = k+1
                lgs%a_index(k) = nr                                 ! for H(i,j)
                lgs%a_value(k) = 1.0                                ! (diagonal element)
                if (uy_1(i,j) < 0.0) &
                    lgs%a_value(k) = lgs%a_value(k) &
                                    - dt_darea*uy_1(i,j)*dx*WOVI
                if (ux_1(i,j) < 0.0) &
                    lgs%a_value(k) = lgs%a_value(k) &
                                    - dt_darea*ux_1(i,j)*dy*WOVI
                if (ux_2(i,j) > 0.0) &
                    lgs%a_value(k) = lgs%a_value(k) &
                                    + dt_darea*ux_2(i,j)*dy*WOVI
                if (uy_2(i,j) > 0.0) &
                    lgs%a_value(k) = lgs%a_value(k) &
                                    + dt_darea*uy_2(i,j)*dx*WOVI

                k = k+1
                lgs%a_index(k) = lgs%ij2n(ip1,j)                    ! for H(ip1,j)
                if (ux_2(i,j) < 0.0) &
                    lgs%a_value(k) = dt_darea*ux_2(i,j)*dy*WOVI

                k = k+1
                lgs%a_index(k) = lgs%ij2n(i,jp1)                    ! for H(i,jp1)
                if (uy_2(i,j) < 0.0) &
                    lgs%a_value(k) = dt_darea*uy_2(i,j)*dx*WOVI


                ! Right-hand side 

                lgs%b_value(nr) = H(i,j) + dt*F(i,j) &
                                -(1.0-WOVI) * dt_darea &
                                     * (  ( ux_2(i,j)*Hx_2(i,j)*dy      &
                                           -ux_1(i,j)*Hx_1(i,j)*dy )    &
                                        + ( uy_2(i,j)*Hy_2(i,j)*dx      &
                                           -uy_1(i,j)*Hy_1(i,j)*dx  ) )

                ! Initial guess == previous H

                lgs%x_value(nr) = H(i,j)                            
                
            end if

            lgs%a_ptr(nr+1) = k+1   ! row is completed, store index to next row
            
        end do

        ! Done: A, x and b matrices in Ax=b have been populated 
        ! and stored in lgs object. 

        return

    end subroutine linear_solver_matrix_advection_csr_2D



! ==== ALTERNATIVES (some useful for benchmark tests) ======

    subroutine calc_adv2D_expl_new(H_ice, f_ice, ux, uy, mdot, dx, dy, dt, boundaries)
        ! Solve 2D advection equation for ice sheet thickness via explicit flux divergence:
        ! d[H]/dt = -grad[H*(ux,uy)] + mdot 
        
        implicit none 

        real(wp), intent(INOUT) :: H_ice(:,:)             ! [m] aa-nodes, Ice thickness 
        real(wp), intent(IN)    :: f_ice(:,:)
        real(wp), intent(IN)    :: ux(:,:)                ! [m a^-1] ac-nodes, Horizontal velocity, x-direction
        real(wp), intent(IN)    :: uy(:,:)                ! [m a^-1] ac-nodes, Horizontal velocity, y-direction
        real(wp), intent(IN)    :: mdot(:,:)              ! [m a^-1] aa-nodes, Source term, net rate of change at top and bottom of cell (no vertical advection) 
        real(wp), intent(IN)    :: dx                     ! [m] Horizontal grid spacing, x-direction
        real(wp), intent(IN)    :: dy                     ! [m] Horizontal grid spacing, y-direction
        real(wp), intent(IN)    :: dt                     ! [a] Timestep 
        character(len=*), intent(IN) :: boundaries 

        ! Local variables:
        integer    :: i, j, nx, ny 
        integer    :: im1, ip1, jm1, jp1
        integer    :: im2, ip2, jm2, jp2
        real(wp)   :: H_ac(2)
        real(wp)   :: wt_ab(4)
        real(wp)   :: H_ab(4)
        real(wp)   :: ux_ab(4)
        real(wp)   :: uy_ab(4)
        real(wp)   :: fx_ab(4)
        real(wp)   :: fy_ab(4) 
        real(wp)   :: wt
        real(wp)   :: ux_now, uy_now 
        real(wp)   :: dfxdx, dfydy 
        real(wp), allocatable :: dHdt(:,:)                ! [m] aa-nodes, Total change this timestep due to fluxes divergence and mdot 
        real(wp), allocatable :: fx(:,:) 
        real(wp), allocatable :: fy(:,:) 

        real(wp), parameter :: dHdt_lim = 1e3             ! [m a-1] Maximum rate of change allowed (high value for extreme changes)

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        ! Initialize dHdt 
        allocate(dHdt(nx,ny))
        dHdt = 0.0_wp 

        allocate(fx(nx,ny))
        allocate(fy(nx,ny))

        fx = 0.0_wp 
        fy = 0.0_wp 

        ! First calculate flux on aa-nodes over entire domain 
        do i = 1, nx
        do j = 1, ny

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

            ! Get ab-node weighting based on whether ice is present 
            wt_ab = 0.0_wp 
            if (count([f_ice(i,j),f_ice(ip1,j),f_ice(i,jp1),f_ice(ip1,jp1)].lt.1.0_wp) .eq. 0) then 
                wt_ab(1) = 1.0_wp
            end if
            if (count([f_ice(i,j),f_ice(im1,j),f_ice(i,jp1),f_ice(im1,jp1)].lt.1.0_wp) .eq. 0) then 
                wt_ab(2) = 1.0_wp 
            end if 
            if (count([f_ice(i,j),f_ice(im1,j),f_ice(i,jm1),f_ice(im1,jm1)].lt.1.0_wp) .eq. 0) then 
                wt_ab(3) = 1.0_wp
            end if 
            if (count([f_ice(i,j),f_ice(ip1,j),f_ice(i,jm1),f_ice(ip1,jm1)].lt.1.0_wp) .eq. 0) then 
                wt_ab(4) = 1.0_wp 
            end if 
            
            wt = sum(wt_ab)
            
            H_ab(1) = 0.25_wp*(H_ice(i,j)+H_ice(ip1,j)+H_ice(i,jp1)+H_ice(ip1,jp1)) 
            H_ab(2) = 0.25_wp*(H_ice(i,j)+H_ice(im1,j)+H_ice(i,jp1)+H_ice(im1,jp1)) 
            H_ab(3) = 0.25_wp*(H_ice(i,j)+H_ice(im1,j)+H_ice(i,jm1)+H_ice(im1,jm1)) 
            H_ab(4) = 0.25_wp*(H_ice(i,j)+H_ice(ip1,j)+H_ice(i,jm1)+H_ice(ip1,jm1)) 

            ! === ux =========

            ux_ab(1) = 0.5_wp*(ux(i,j)+ux(i,jp1))
            ux_ab(2) = 0.5_wp*(ux(im1,j)+ux(im1,jp1))
            ux_ab(3) = 0.5_wp*(ux(im1,jm1)+ux(im1,j))
            ux_ab(4) = 0.5_wp*(ux(i,jm1)+ux(i,j))
            where (abs(ux_ab) .lt. TOL_UNDERFLOW) ux_ab = 0.0_wp 

            ! === uy =========

            uy_ab(1) = 0.5_wp*(uy(i,j)+uy(ip1,j))
            uy_ab(2) = 0.5_wp*(uy(i,j)+uy(im1,j))
            uy_ab(3) = 0.5_wp*(uy(i,jm1)+uy(im1,jm1))
            uy_ab(4) = 0.5_wp*(uy(i,jm1)+uy(ip1,jm1))
            where (abs(uy_ab) .lt. TOL_UNDERFLOW) uy_ab = 0.0_wp 

            fx_ab = H_ab*ux_ab
            fy_ab = H_ab*uy_ab 

            dfxdx = 0.5_wp*( (fx_ab(1)-fx_ab(2))/dx + (fx_ab(4)-fx_ab(3))/dx )
            dfydy = 0.5_wp*( (fy_ab(1)-fy_ab(4))/dy + (fy_ab(2)-fy_ab(3))/dy )
            
            dHdt(i,j) = -dfxdx - dfydy 

            ! Limit dHdt for stability 
            if (dHdt(i,j) .lt. -dHdt_lim) dHdt(i,j) = -dHdt_lim 
            if (dHdt(i,j) .gt.  dHdt_lim) dHdt(i,j) =  dHdt_lim 
            


            ! ! Get centered fluxes (aa-nodes)
            ! fx(i,j) = 0.25_wp*(sum(H_ab*ux_ab))
            ! fy(i,j) = 0.25_wp*(sum(H_ab*uy_ab))
            
            ! ! Get centered fluxes (aa-nodes)
            ! H_ac(1) = 0.5_wp*(H_ice(i,j)+H_ice(im1,j))
            ! H_ac(2) = 0.5_wp*(H_ice(i,j)+H_ice(ip1,j))
            ! fx(i,j) = 0.5_wp*(H_ac(1)*ux(im1,j) + H_ac(2)*ux(i,j))

            ! H_ac(1) = 0.5_wp*(H_ice(i,j)+H_ice(i,jm1))
            ! H_ac(2) = 0.5_wp*(H_ice(i,j)+H_ice(i,jp1))
            ! fy(i,j) = 0.5_wp*(H_ac(1)*uy(i,jm1) + H_ac(2)*uy(i,j))
            
        end do
        end do 

if (.FALSE.) then
        ! Now calculate spatial derivatives dfx/dx / dfy/dy
        ! to obtain dh/dt
        do i = 1, nx
        do j = 1, ny

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

            !im2 = max(i-2,1)
            !ip2 = min(i+2,nx)
            !jm2 = max(j-2,1)
            !jp2 = min(j+2,ny)
            
            ! ajr: works for simple (EISMINT tests), but not for non-linear
            ! realistic ice sheets
            ! Second-order spatial derivative following
            ! https://web.media.mit.edu/~crtaylor/calculator.html
            dfxdx = (1.0_wp/(2.0_wp*dx))*(-1.0_wp*fx(im1,j)+1.0_wp*fx(ip1,j))
            dfydy = (1.0_wp/(2.0_wp*dy))*(-1.0_wp*fy(i,jm1)+1.0_wp*fy(i,jp1))

            ! ajr: works for simple (EISMINT tests), but not for non-linear
            ! realistic ice sheets
            ! Fourth order spatial derivative following
            ! Man and Tsai (2008).
            ! dfxdx = (1.0_wp/(12.0_wp*dx))*(fx(im2,j)-8.0_wp*fx(im1,j)+8.0_wp*fx(ip1,j)-fx(ip2,j))
            ! dfydy = (1.0_wp/(12.0_wp*dy))*(fy(i,jm2)-8.0_wp*fy(i,jm1)+8.0_wp*fy(i,jp1)-fy(i,jp1))


            dHdt(i,j) = -dfxdx - dfydy 

            ! Limit dHdt for stability 
            if (dHdt(i,j) .lt. -dHdt_lim) dHdt(i,j) = -dHdt_lim 
            if (dHdt(i,j) .gt.  dHdt_lim) dHdt(i,j) =  dHdt_lim 
            
        end do
        end do

end if 

        ! Update H_ice:
        H_ice = H_ice + dt*dHdt + dt*mdot 

        return 

    end subroutine calc_adv2D_expl_new

    subroutine calc_adv2D_expl(H_ice, f_ice, ux, uy, mdot, dx, dy, dt, boundaries)
        ! Solve 2D advection equation for ice sheet thickness via explicit flux divergence:
        ! d[H]/dt = -grad[H*(ux,uy)] + mdot 
        !
        ! ajr: adapted from IMAU-ICE code from Heiko Goelzer (h.goelzer@uu.nl) 2016
        ! Note: It also works using the ac-nodes directly, but is less stable, particularly
        ! for EISMINT2-expf. 

        implicit none 

        real(wp), intent(INOUT) :: H_ice(:,:)               ! [m] aa-nodes, Ice thickness 
        real(wp), intent(IN)    :: f_ice(:,:)               ! Ice fraction
        real(wp), intent(IN)    :: ux(:,:)                  ! [m a^-1] ac-nodes, Horizontal velocity, x-direction
        real(wp), intent(IN)    :: uy(:,:)                  ! [m a^-1] ac-nodes, Horizontal velocity, y-direction
        real(wp), intent(IN)    :: mdot(:,:)                ! [m a^-1] aa-nodes, Source term, net rate of change at top and bottom of cell (no vertical advection) 
        real(wp), intent(IN)    :: dx                       ! [m] Horizontal grid spacing, x-direction
        real(wp), intent(IN)    :: dy                       ! [m] Horizontal grid spacing, y-direction
        real(wp), intent(IN)    :: dt                       ! [a] Timestep 
        character(len=*), intent(IN) :: boundaries 
        
        ! Local variables:
        integer  :: i, j, nx, ny 
        integer  :: im1, ip1, jm1, jp1
        real(wp) :: H_ab(4)
        real(wp) :: flux_xr                                 ! [m^2 a^-1] ac-nodes, Flux in the x-direction to the right
        real(wp) :: flux_xl                                 ! [m^2 a^-1] ac-nodes, Flux in the x-direction to the left
        real(wp) :: flux_yu                                 ! [m^2 a^-1] ac-nodes, Flux in the y-direction upwards
        real(wp) :: flux_yd                                 ! [m^2 a^-1] ac-nodes, Flux in the y-direction downwards
        real(wp), allocatable :: dHdt(:,:)                  ! [m] aa-nodes, Total change this timestep due to fluxes divergence and mdot 
        real(wp), allocatable :: f_ice_now(:,:)
        real(wp), parameter :: dHdt_lim = 1e3               ! [m a-1] Maximum rate of change allowed (high value for extreme changes)

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        ! Initialize dHdt 
        allocate(dHdt(nx,ny))
        dHdt = 0.0_wp 

        allocate(f_ice_now(nx,ny))
        f_ice_now = 1.0 

        do j = 1, ny
        do i = 1, nx
        
            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

            ! Get the ice thickness on ab-nodes
            call stagger_nodes_aa_ab_ice(H_ab,H_ice,f_ice_now,i,j)

            ! Calculate the flux across each boundary [m^2 a^-1]
            flux_xr = ux(i,j)   * 0.5*(H_ab(1)+H_ab(4))
            flux_xl = ux(im1,j) * 0.5*(H_ab(2)+H_ab(3))
            flux_yu = uy(i,j)   * 0.5*(H_ab(1)+H_ab(2))
            flux_yd = uy(i,jm1) * 0.5*(H_ab(3)+H_ab(4))

            ! Calculate flux divergence on aa-nodes 
            dHdt(i,j) = (1.0 / dx) * (flux_xl - flux_xr) + (1.0 / dy) * (flux_yd - flux_yu)

            ! Limit dHdt for stability 
            if (dHdt(i,j) .lt. -dHdt_lim) dHdt(i,j) = -dHdt_lim 
            if (dHdt(i,j) .gt.  dHdt_lim) dHdt(i,j) =  dHdt_lim 
            
        end do
        end do

        ! Update H_ice:
        H_ice = H_ice + dt*dHdt + dt*mdot 

        return 

    end subroutine calc_adv2D_expl

    subroutine calc_adv2D_expl_upwind(H_ice, ux, uy, mdot, dx, dy, dt, boundaries)
        ! Solve 2D advection equation for ice sheet thickness via explicit flux divergence:
        ! d[H]/dt = -grad[H*(ux,uy)] + mdot 

        ! Second-order upwind routine, see eg, Winkelmann et al (2011), Eq. 18 

        implicit none 

        real(wp), intent(INOUT) :: H_ice(:,:)             ! [m] aa-nodes, Ice thickness 
        real(wp), intent(IN)    :: ux(:,:)                ! [m a^-1] ac-nodes, Horizontal velocity, x-direction
        real(wp), intent(IN)    :: uy(:,:)                ! [m a^-1] ac-nodes, Horizontal velocity, y-direction
        real(wp), intent(IN)    :: mdot(:,:)              ! [m a^-1] aa-nodes, Source term, net rate of change at top and bottom of cell (no vertical advection) 
        real(wp), intent(IN)    :: dx                     ! [m] Horizontal grid spacing, x-direction
        real(wp), intent(IN)    :: dy                     ! [m] Horizontal grid spacing, y-direction
        real(wp), intent(IN)    :: dt                     ! [a] Timestep 
        character(len=*), intent(IN) :: boundaries 
        
        ! Local variables:
        integer    :: i, j, nx, ny
        integer    :: im1, ip1, jm1, jp1 
        real(wp) :: H_now 
        real(wp) :: flux_xr                               ! [m^2 a^-1] ac-nodes, Flux in the x-direction to the right
        real(wp) :: flux_xl                               ! [m^2 a^-1] ac-nodes, Flux in the x-direction to the left
        real(wp) :: flux_yu                               ! [m^2 a^-1] ac-nodes, Flux in the y-direction upwards
        real(wp) :: flux_yd                               ! [m^2 a^-1] ac-nodes, Flux in the y-direction downwards
        real(wp), allocatable :: dHdt(:,:)                ! [m] aa-nodes, Total change this timestep due to fluxes divergence and mdot 

        real(wp), parameter :: dHdt_lim = 1e3             ! [m a-1] Maximum rate of change allowed (high value for extreme changes)

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        ! Initialize dHdt 
        allocate(dHdt(nx,ny))
        dHdt = 0.0_wp 

        do j = 1, ny
        do i = 1, nx
        
            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

            ! Calculate the flux across each boundary [m^2 a^-1]
            if (ux(i,j) .gt. 0.0) then 
                flux_xr = ux(i,j)*H_ice(i,j)
            else
                flux_xr = ux(i,j)*H_ice(ip1,j) 
            end if 

            if (ux(im1,j) .gt. 0.0) then 
                flux_xl = ux(im1,j)*H_ice(im1,j)
            else
                flux_xl = ux(im1,j)*H_ice(i,j) 
            end if 

            if (uy(i,j) .gt. 0.0) then 
                flux_yu = uy(i,j)*H_ice(i,j)
            else
                flux_yu = uy(i,j)*H_ice(i,jp1) 
            end if 

            if (uy(i,jm1) .gt. 0.0) then 
                flux_yd = uy(i,jm1)*H_ice(i,jm1)
            else
                flux_yd = uy(i,jm1)*H_ice(i,j) 
            end if 

            ! Calculate flux divergence on aa-nodes 
            dHdt(i,j) = (1.0 / dx) * (flux_xl - flux_xr) + (1.0 / dy) * (flux_yd - flux_yu)

            ! Limit dHdt for stability 
            if (dHdt(i,j) .lt. -dHdt_lim) dHdt(i,j) = -dHdt_lim 
            if (dHdt(i,j) .gt.  dHdt_lim) dHdt(i,j) =  dHdt_lim 
            
        end do
        end do

        ! Update H_ice:
        H_ice = H_ice + dt*dHdt + dt*mdot 

        return 

    end subroutine calc_adv2D_expl_upwind

    subroutine  calc_adv2D_impl_upwind(H,ux,uy,mdot,dx,dy,dt,boundaries,f_upwind)
        ! To solve the 2D adevection equation:
        ! dH/dt =
        ! M H = Frelax
        ! ajr: adapted from GRISLI (Ritz et al., 1997)

        implicit none

        real(wp), intent(INOUT) :: H(:,:)         ! Ice thickness (aa-node)
        real(wp), intent(IN)    :: ux(:,:)        ! Depth-averaged velocity - x direction (ac-node)
        real(wp), intent(IN)    :: uy(:,:)        ! Depth-averaged velocity - y direction (ac-node)
        real(wp), intent(IN)    :: mdot(:,:)      ! Total column mass balance (aa-node)
        real(wp), intent(IN)    :: dx             ! [m] x-resolution
        real(wp), intent(IN)    :: dy             ! [m] y-resolution
        real(wp), intent(IN)    :: dt             ! [a] Timestep (assumes dx=dy)
        character(len=*), intent(IN) :: boundaries 
        real(wp), intent(IN)    :: f_upwind       ! [-] Fraction of "upwind-ness" to apply (ajr: experimental!) - between 0.5 and 1.0, default f_upwind=1.0
        
        ! Local variables
        integer    :: i, j, nx, ny
        integer    :: im1, ip1, jm1, jp1 
        integer    :: iter, ierr 
        real(wp) :: dtdx, dtdx2
        real(wp) :: reste, delh, tmp 
        real(wp) :: ux_i, ux_im1, uy_j, uy_jm1
        real(wp), allocatable :: crelax(:,:)      ! diagnonale de M
        real(wp), allocatable :: arelax(:,:)      ! sous diagonale selon x
        real(wp), allocatable :: brelax(:,:)      ! sur  diagonale selon x
        real(wp), allocatable :: drelax(:,:)      ! sous diagonale selon y
        real(wp), allocatable :: erelax(:,:)      ! sur  diagonale selon y
        real(wp), allocatable :: frelax(:,:)      ! vecteur
        real(wp), allocatable :: c_west(:,:)      ! sur demi mailles Ux
        real(wp), allocatable :: c_east(:,:)      ! sur demi mailles Ux
        real(wp), allocatable :: c_north(:,:)     ! sur demi mailles Uy
        real(wp), allocatable :: c_south(:,:)     ! sur demi mailles Uy
        real(wp), allocatable :: deltaH(:,:)      ! Change in H

        logical,    parameter :: use_upwind = .TRUE.  ! Apply upwind advection scheme or central scheme?   
                                                      ! (now this is redundant with f_upwind parameter) 
        
        ! Note: f_upwind=0.6 gives good agreement with EISMINT1 summit elevation
        ! for the moving and fixed margin experiments, when using the calc_shear_3D approach.
        ! (f_upwind=0.5, ie central method, works well when using the velocity_sia approach)

        ! Note: it may be that f_upwind>0.5 is more stable for real dynamic simulations

        ! Determine array size
        nx = size(H,1)
        ny = size(H,2)

        ! Allocate local arrays
        allocate(crelax(nx,ny))
        allocate(arelax(nx,ny))
        allocate(brelax(nx,ny))
        allocate(drelax(nx,ny))
        allocate(erelax(nx,ny))
        allocate(frelax(nx,ny))
        allocate(c_west(nx,ny))
        allocate(c_east(nx,ny))
        allocate(c_north(nx,ny))
        allocate(c_south(nx,ny))

        allocate(deltaH(nx,ny))

        ! Define some helpful values
        dtdx2 = dt/(dx**2)
        dtdx  = dt/dx

        ! Initialize relaxation arrays
        arelax = 0.0
        brelax = 0.0
        drelax = 0.0
        erelax = 0.0
        crelax = 1.0
        frelax = 0.0
        deltaH = 0.0

        ! Modify coefficients depending on method (upwind, central)
        if (use_upwind) then
            ! Upwind method

            if (f_upwind .eq. 1.0) then
                ! Apply normal upwind scheme

                where (ux.ge.0.0)
                    c_west = 1.0
                    c_east = 0.0
                elsewhere
                    c_west = 0.0
                    c_east = 1.0
                end where

                where (uy.ge.0.0)
                    c_south = 1.0
                    c_north = 0.0
                elsewhere
                    c_south = 0.0
                    c_north = 1.0
                end where

            else 
                ! Apply fractional upwind scheme to reduce numerical diffusion,
                ! but benefit from upwind stability (ajr: experimental!)
                ! f_upwind = 0.5 => central difference, f_upwind = 1.0 => full upwind 

                where (ux.ge.0.0)
                    c_west = f_upwind
                    c_east = 1.0 - f_upwind
                elsewhere
                    c_west = 1.0 - f_upwind
                    c_east = f_upwind
                end where

                where (uy.ge.0.0)
                    c_south = f_upwind
                    c_north = 1.0 - f_upwind
                elsewhere
                    c_south = 1.0 - f_upwind
                    c_north = f_upwind
                end where

            end if 

        else
            ! Central method

            c_west  = 0.5
            c_east  = 0.5
            c_south = 0.5
            c_north = 0.5

        end if

        ! Populate diagonals
        do j = 1, ny
        do i = 1, nx

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

            !  sous diagonale en x
            ux_im1 = ux(im1,j)
            if (abs(ux_im1) .lt. TOL_UNDERFLOW) ux_im1 = 0.0
            arelax(i,j) = -dtdx*c_west(im1,j)*ux_im1    ! partie advective en x

            !  sur diagonale en x
            ux_i = ux(i,j)
            if (abs(ux_i) .lt. TOL_UNDERFLOW) ux_i = 0.0
            brelax(i,j) = +dtdx*c_east(i,j)*ux_i        ! partie advective

            !  sous diagonale en y
            uy_jm1 = uy(i,jm1)
            if (abs(uy_jm1) .lt. TOL_UNDERFLOW) uy_jm1 = 0.0
            drelax(i,j) = -dtdx*c_south(i,jm1)*uy_jm1   ! partie advective en y

            !  sur diagonale en y
            uy_j = uy(i,j)
            if (abs(uy_j) .lt. TOL_UNDERFLOW) uy_j = 0.0
            erelax(i,j) = +dtdx*c_north(i,j)*uy_j       ! partie advective


            ! diagonale
            crelax(i,j) = 1.0 + dtdx* &
                       ((c_west(i,j)*ux_i - c_east(im1,j)*ux_im1) &
                      +(c_south(i,j)*uy_j - c_north(i,jm1)*uy_jm1))

            ! Combine all terms
            frelax(i,j) = H(i,j) + dt*mdot(i,j)

        end do
        end do

        ! Avoid underflows 
        where (abs(arelax) .lt. tol_underflow) arelax = 0.0_wp 
        where (abs(brelax) .lt. tol_underflow) brelax = 0.0_wp 
        where (abs(drelax) .lt. tol_underflow) drelax = 0.0_wp 
        where (abs(erelax) .lt. tol_underflow) erelax = 0.0_wp 
        
        ! Initialize new H solution to zero (to get zeros at boundaries)
        H  = 0.0

        ! Initially assume convergence criterion is not satisfied 
        ierr = -1   ! convergence criterion not fulfilled
        
        do iter = 1, 1000 
            ! Relaxation loop 

            ! Calculate change in H
            do j = 1, ny
            do i = 1, nx

                ! Get neighbor indices
                call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

                reste = (((arelax(i,j)*H(im1,j) + drelax(i,j)*H(i,jm1)) &
                        + (brelax(i,j)*H(ip1,j) + erelax(i,j)*H(i,jp1))) &
                        +  crelax(i,j)*H(i,j))  - frelax(i,j)

                deltaH(i,j) = reste/crelax(i,j)

            end do
            end do

            ! Avoid underflows
            where(abs(deltaH) .lt. tol_underflow) deltaH = 0.0_wp
            
            ! Adjust H to new value
            H = H - deltaH

            ! Avoid underflows
            where(abs(H) .lt. tol_underflow) H = 0.0_wp
            
            ! Check stopping criterion (something like rmse of remaining change in H)
            delh = sqrt(sum(deltaH**2)) / (nx*ny)
            
            ! Use simple stopping criterion: maximum remaining change in H
            ! Note: this is less likely to converge given the same stopping
            ! criterion.
!             delh = maxval(abs(deltaH))
            
            if ( delh .lt. 1e-6) then
                ! Solution has converged, exit  
                ierr = 0 
                exit 
            end if 

        end do ! End of relaxation loop

        !write(6,'(10x,a,i0,5x,i2)') 'calc_adv2D_impl_upwind: iter = ', iter, ierr
        
        return

    end subroutine calc_adv2D_impl_upwind

end module solver_advection
