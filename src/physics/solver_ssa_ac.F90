module solver_ssa_ac


! Uncomment next line if petsc should be used
!#define USEPETSC 



    use yelmo_defs, only : sp, dp, wp, io_unit_err, TOL_UNDERFLOW, rho_ice, rho_sw, g 
    use yelmo_tools, only : get_neighbor_indices

    use grid_calcs  ! For staggering routines 
    use ncio        ! For diagnostic outputting only 

    implicit none 

    type linear_solver_class

        integer :: nmax
        integer :: n_terms
        integer :: n_sprs
        
        integer,  allocatable :: n2i(:)
        integer,  allocatable :: n2j(:)
        integer,  allocatable :: ij2n(:,:)
        
        integer,  allocatable :: a_ptr(:)
        integer,  allocatable :: a_index(:)
        real(wp), allocatable :: a_value(:)
        real(wp), allocatable :: b_value(:)
        real(wp), allocatable :: x_value(:)

        real(wp), allocatable :: resid(:)
        real(wp) :: L1_norm, L2_norm, L2_rel_norm

    end type

    private 
    public :: set_ssa_masks
    public :: update_ssa_mask_convergence
    public :: ssa_diagnostics_write_init
    public :: ssa_diagnostics_write_step

    public :: linear_solver_class
    public :: linear_solver_init
    public :: linear_solver_matrix_ssa_ac_csr_2D
    public :: linear_solver_matrix_solve
    public :: linear_solver_save_velocity

contains 
    
    subroutine linear_solver_init(lgs,nx,ny)
        ! Initialize the lgs object that will hold
        ! arrays A, x and b needed to solve linear
        ! equation Ax=b. Arrays are populated in 
        ! another subroutine below - see
        ! linear_solver_matrix_ssa_ac_csr_2D. 

        implicit none

        type(linear_solver_class), intent(INOUT) :: lgs
        integer, intent(IN) :: nx, ny 

        ! Local variables
        integer :: i, j, n 

        ! Define array sizes
        lgs%nmax    = 2*nx*ny 
        lgs%n_terms = 9
        lgs%n_sprs  = 2*lgs%n_terms*nx*ny 

        ! Ensure all object arrays are deallocated first
        if (allocated(lgs%n2i))     deallocate(lgs%n2i)
        if (allocated(lgs%n2j))     deallocate(lgs%n2j)
        if (allocated(lgs%ij2n))    deallocate(lgs%ij2n)
        if (allocated(lgs%a_ptr))   deallocate(lgs%a_ptr)
        if (allocated(lgs%a_index)) deallocate(lgs%a_index)
        if (allocated(lgs%a_value)) deallocate(lgs%a_value)
        if (allocated(lgs%b_value)) deallocate(lgs%b_value)
        if (allocated(lgs%x_value)) deallocate(lgs%x_value)
        
        ! Allocate arrays to proper size
        allocate(lgs%n2i(nx*ny))
        allocate(lgs%n2j(nx*ny))
        allocate(lgs%ij2n(nx,ny))
        allocate(lgs%a_value(lgs%n_sprs))
        allocate(lgs%a_index(lgs%n_sprs))
        allocate(lgs%a_ptr(lgs%nmax+1))
        allocate(lgs%b_value(lgs%nmax))
        allocate(lgs%x_value(lgs%nmax))
        
        ! Initialize array values to zero
        lgs%a_value = 0.0
        lgs%a_index = 0
        lgs%a_ptr   = 0

        lgs%b_value = 0.0
        lgs%x_value = 0.0
        
        ! Define indices for reshaping of a 2-d array (with indices i, j)
        ! to a vector (with index n)

        n = 1

        do i = 1, nx
        do j = 1, ny
            lgs%n2i(n)    = i
            lgs%n2j(n)    = j
            lgs%ij2n(i,j) = n
            n = n+1
        end do
        end do

        return

    end subroutine linear_solver_init

    subroutine linear_solver_save_velocity(ux,uy,lgs,ulim)
        ! Extract velocity solution from lgs object. 

        implicit none 

        real(wp), intent(OUT) :: ux(:,:)                ! [m yr^-1] Horizontal velocity x
        real(wp), intent(OUT) :: uy(:,:)                ! [m yr^-1] Horizontal velocity y
        type(linear_solver_class), intent(IN) :: lgs 
        real(wp), intent(IN)    :: ulim 

        ! Local variables 
        integer :: i, j, n, nr 

        do n = 1, lgs%nmax-1, 2

            i = lgs%n2i((n+1)/2)
            j = lgs%n2j((n+1)/2)

            nr = n
            ux(i,j) = lgs%x_value(nr)

            nr = n+1
            uy(i,j) = lgs%x_value(nr)

        end do

        ! Limit the velocity generally =====================
        call limit_vel(ux,ulim)
        call limit_vel(uy,ulim)

        return

    end subroutine linear_solver_save_velocity

    subroutine linear_solver_matrix_ssa_ac_csr_2D(lgs,ux,uy,beta_acx,beta_acy, &
                            N_aa,ssa_mask_acx,ssa_mask_acy,mask_frnt,H_ice,f_ice,taud_acx, &
                            taud_acy,taul_int_acx,taul_int_acy,dx,dy,boundaries,lateral_bc)
        ! Define sparse matrices A*x=b in format 'compressed sparse row' (csr)
        ! for the SSA momentum balance equations with velocity components
        ! ux and uy defined on ac-nodes (right and top borders of i,j grid cell)
        ! Store sparse matrices in linear_solver_class object 'lgs' for later use.

        implicit none 

        type(linear_solver_class), intent(INOUT) :: lgs
        real(wp), intent(IN) :: ux(:,:)                 ! [m yr^-1] Horizontal velocity x (acx-nodes)
        real(wp), intent(IN) :: uy(:,:)                 ! [m yr^-1] Horizontal velocity y (acy-nodes)
        real(wp), intent(IN) :: beta_acx(:,:)           ! [Pa yr m^-1] Basal friction (acx-nodes)
        real(wp), intent(IN) :: beta_acy(:,:)           ! [Pa yr m^-1] Basal friction (acy-nodes)
        real(wp), intent(IN) :: N_aa(:,:)               ! [Pa yr m] Vertically integrated viscosity (aa-nodes)
        integer,  intent(IN) :: ssa_mask_acx(:,:)       ! [--] Mask to determine ssa solver actions (acx-nodes)
        integer,  intent(IN) :: ssa_mask_acy(:,:)       ! [--] Mask to determine ssa solver actions (acy-nodes)
        integer,  intent(IN) :: mask_frnt(:,:)          ! [--] Ice-front mask 
        real(wp), intent(IN) :: H_ice(:,:)              ! [m]  Ice thickness (aa-nodes)
        real(wp), intent(IN) :: f_ice(:,:)
        real(wp), intent(IN) :: taud_acx(:,:)           ! [Pa] Driving stress (acx nodes)
        real(wp), intent(IN) :: taud_acy(:,:)           ! [Pa] Driving stress (acy nodes)
        real(wp), intent(IN) :: taul_int_acx(:,:)       ! [Pa m] Vertically integrated lateral stress (acx nodes)
        real(wp), intent(IN) :: taul_int_acy(:,:)       ! [Pa m] Vertically integrated lateral stress (acy nodes) 
        real(wp), intent(IN) :: dx, dy
        character(len=*), intent(IN) :: boundaries 
        character(len=*), intent(IN) :: lateral_bc

        ! Local variables
        integer  :: nx, ny
        integer  :: i, j, k, n, m 
        integer  :: nc, nr
        real(wp) :: inv_dx, inv_dxdx 
        real(wp) :: inv_dy, inv_dydy 
        real(wp) :: inv_dxdy, inv_2dxdy, inv_4dxdy 
        
        real(wp), allocatable :: N_ab(:,:)

        ! Boundary conditions (bcs) counterclockwise unit circle 
        ! 1: x, right-border
        ! 2: y, upper-border 
        ! 3: x, left--border 
        ! 4: y, lower-border 
        character(len=56) :: bcs(4)

        integer :: im1, ip1, jm1, jp1 
        real(wp) :: N_aa_now

        nx = size(H_ice,1)
        ny = size(H_ice,2) 

        ! Safety check for initialization
        if (.not. allocated(lgs%x_value)) then 
            ! Object 'lgs' has not been initialized yet, do so now.

            call linear_solver_init(lgs,nx,ny)

        end if 

        ! Consistency check: ensure beta is defined well 
        if ( ( count(ssa_mask_acx .eq. 1 .and. beta_acx .gt. 0.0) .eq. 0 ) .or. & 
             ( count(ssa_mask_acy .eq. 1 .and. beta_acy .gt. 0.0) .eq. 0 ) ) then  
            ! No points found with a non-zero beta for grounded ice,
            ! something was not well-defined/well-initialized

            write(*,*) 
            write(*,"(a)") "linear_solver_matrix_ssa_ac_csr_2D:: Error: beta appears to be zero everywhere for grounded ice."
            write(*,*) "range(beta_acx): ", minval(beta_acx), maxval(beta_acx)
            write(*,*) "range(beta_acy): ", minval(beta_acy), maxval(beta_acy)
            write(*,*) "range(ssa_mask_acx): ", minval(ssa_mask_acx), maxval(ssa_mask_acx)
            write(*,*) "range(ssa_mask_acy): ", minval(ssa_mask_acy), maxval(ssa_mask_acy)
            write(*,*) "Stopping."
            write(*,*) 
            stop 
            
        end if 

        ! Define border conditions (no-slip, free-slip, periodic)
        select case(trim(boundaries)) 

            case("MISMIP3D")

                bcs(1) = "free-slip"
                bcs(2) = "free-slip"
                bcs(3) = "no-slip"
                bcs(4) = "free-slip" 
                ! bcs(1) = "free-slip"
                ! bcs(2) = "periodic"
                ! bcs(3) = "no-slip"
                ! bcs(4) = "periodic" 

            case("periodic")

                bcs(1:4) = "periodic" 

            case("periodic-x")

                bcs(1) = "periodic"
                bcs(2) = "free-slip"
                bcs(3) = "periodic"
                bcs(4) = "free-slip"

            case("periodic-y")

                bcs(1) = "free-slip"
                bcs(2) = "periodic"
                bcs(3) = "free-slip"
                bcs(4) = "periodic"
            
            case("infinite")

                bcs(1:4) = "free-slip" 
                
            case DEFAULT 

                bcs(1:4) = "periodic"
                
        end select 

        nx = size(H_ice,1)
        ny = size(H_ice,2)
        
        allocate(N_ab(nx,ny))

        ! Define some factors

        inv_dx      = 1.0_wp / dx 
        inv_dxdx    = 1.0_wp / (dx*dx)
        inv_dy      = 1.0_wp / dy 
        inv_dydy    = 1.0_wp / (dy*dy)
        inv_dxdy    = 1.0_wp / (dx*dy)
        inv_2dxdy   = 1.0_wp / (2.0_wp*dx*dy)
        inv_4dxdy   = 1.0_wp / (4.0_wp*dx*dy)
        

        ! Calculate the staggered depth-integrated viscosity 
        ! at the grid-cell corners (ab-nodes). 
        call stagger_visc_aa_ab(N_ab,N_aa,H_ice,f_ice,boundaries)
        

        !-------- Assembly of the system of linear equations
        !             (matrix storage: compressed sparse row CSR) --------

        lgs%a_ptr(1) = 1

        k = 0

        do n=1, lgs%nmax-1, 2

            i = lgs%n2i((n+1)/2)
            j = lgs%n2j((n+1)/2)

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

            ! ------ Equations for ux ---------------------------

            nr = n   ! row counter

            ! == Treat special cases first ==

            if (ssa_mask_acx(i,j) .eq. 0) then
                ! SSA set to zero velocity for this point

                k = k+1
                lgs%a_value(k) = 1.0_wp   ! diagonal element only
                lgs%a_index(k) = nr

                lgs%b_value(nr) = 0.0_wp
                lgs%x_value(nr) = 0.0_wp
            
            else if (ssa_mask_acx(i,j) .eq. -1) then 
                ! Assign prescribed boundary velocity to this point
                ! (eg for prescribed velocity corresponding to 
                ! analytical grounding line flux, or for a regional domain)

                k = k+1
                lgs%a_value(k)  = 1.0   ! diagonal element only
                lgs%a_index(k)  = nr

                lgs%b_value(nr) = ux(i,j)
                lgs%x_value(nr) = ux(i,j)
            
            else if (i .eq. 1 .and. trim(bcs(3)) .ne. "periodic") then 
                ! Left boundary 

                select case(trim(bcs(3)))

                    case("no-slip")

                        k = k+1
                        lgs%a_value(k)  = 1.0   ! diagonal element only
                        lgs%a_index(k)  = nr

                        lgs%b_value(nr) = 0.0
                        lgs%x_value(nr) = 0.0

                    case("free-slip")

                        nc = 2*lgs%ij2n(i,j)-1          ! column counter for ux(i,j)
                        k = k+1
                        lgs%a_value(k) =  1.0_wp
                        lgs%a_index(k) = nc

                        nc = 2*lgs%ij2n(ip1,j)-1        ! column counter for ux(ip1,j)
                        k = k+1
                        lgs%a_value(k) = -1.0_wp
                        lgs%a_index(k) = nc

                        lgs%b_value(nr) = 0.0_wp
                        lgs%x_value(nr) = ux(i,j)

                end select 
                
            else if (i .eq. nx .and. trim(bcs(1)) .ne. "periodic") then 
                ! Right boundary 
                
                select case(trim(bcs(1)))

                    case("no-slip")

                        k = k+1
                        lgs%a_value(k)  = 1.0   ! diagonal element only
                        lgs%a_index(k)  = nr

                        lgs%b_value(nr) = 0.0
                        lgs%x_value(nr) = 0.0

                    case("free-slip")

                        nc = 2*lgs%ij2n(i,j)-1          ! column counter for ux(i,j)
                        k = k+1
                        lgs%a_value(k) =  1.0_wp
                        lgs%a_index(k) = nc

                        nc = 2*lgs%ij2n(nx-1,j)-1       ! column counter for ux(nx-1,j)
                        k = k+1
                        lgs%a_value(k) = -1.0_wp
                        lgs%a_index(k) = nc

                        lgs%b_value(nr) = 0.0_wp
                        lgs%x_value(nr) = ux(i,j)

                end select 

            else if (j .eq. 1 .and. trim(bcs(4)) .ne. "periodic") then 
                ! Lower boundary 

                select case(trim(bcs(4)))

                    case("no-slip")

                        k = k+1
                        lgs%a_value(k)  = 1.0_wp   ! diagonal element only
                        lgs%a_index(k)  = nr

                        lgs%b_value(nr) = 0.0_wp
                        lgs%x_value(nr) = 0.0_wp

                    case("free-slip")

                        nc = 2*lgs%ij2n(i,j)-1          ! column counter for ux(i,j)
                        k = k+1
                        lgs%a_value(k) =  1.0_wp
                        lgs%a_index(k) = nc

                        nc = 2*lgs%ij2n(i,jp1)-1       ! column counter for ux(i,jp1)
                        k = k+1
                        lgs%a_value(k) = -1.0_wp
                        lgs%a_index(k) = nc

                        lgs%b_value(nr) = 0.0_wp
                        lgs%x_value(nr) = ux(i,j)

                end select 

            else if (j .eq. ny .and. trim(bcs(2)) .ne. "periodic") then 
                ! Upper boundary 

                select case(trim(bcs(2)))

                    case("no-slip")

                        k = k+1
                        lgs%a_value(k)  = 1.0_wp   ! diagonal element only
                        lgs%a_index(k)  = nr

                        lgs%b_value(nr) = 0.0_wp
                        lgs%x_value(nr) = 0.0_wp
                        
                    case("free-slip")

                        nc = 2*lgs%ij2n(i,j)-1          ! column counter for ux(i,j)
                        k = k+1
                        lgs%a_value(k) =  1.0_wp
                        lgs%a_index(k) = nc

                        nc = 2*lgs%ij2n(i,ny-1)-1       ! column counter for ux(i,ny-1)
                        k = k+1
                        lgs%a_value(k) = -1.0_wp
                        lgs%a_index(k) = nc

                        lgs%b_value(nr) = 0.0_wp
                        lgs%x_value(nr) = ux(i,j)

                end select 

            else if (ssa_mask_acx(i,j) .eq. 3) then 
                ! Lateral boundary condition should be applied here 

                if (f_ice(i,j) .eq. 1.0 .and. f_ice(ip1,j) .lt. 1.0) then 
                    ! === Case 1: ice-free to the right ===

                    N_aa_now = N_aa(i,j)
                    
                    nc = 2*lgs%ij2n(im1,j)-1
                        ! smallest nc (column counter), for ux(im1,j)
                    k = k+1
                    lgs%a_value(k) = -4.0_wp*inv_dx*N_aa_now
                    lgs%a_index(k) = nc 

                    nc = 2*lgs%ij2n(i,jm1)
                        ! next nc (column counter), for uy(i,jm1)
                    k = k+1
                    lgs%a_value(k) = -2.0_wp*inv_dy*N_aa_now
                    lgs%a_index(k) = nc

                    nc = 2*lgs%ij2n(i,j)-1
                        ! next nc (column counter), for ux(i,j)
                    k = k+1
                    lgs%a_value(k) = 4.0_wp*inv_dx*N_aa_now
                    lgs%a_index(k) = nc

                    nc = 2*lgs%ij2n(i,j)
                        ! next nc (column counter), for uy(i,j)
                    k = k+1
                    lgs%a_value(k) = 2.0_wp*inv_dy*N_aa_now
                    lgs%a_index(k) = nc

                    ! Assign matrix values
                    lgs%b_value(nr) = taul_int_acx(i,j) 
                    lgs%x_value(nr) = ux(i,j)
                    
                else 
                    ! Case 2: ice-free to the left
                    
                    N_aa_now = N_aa(ip1,j)
                    
                    nc = 2*lgs%ij2n(i,j)-1
                        ! next nc (column counter), for ux(i,j)
                    k = k+1
                    lgs%a_value(k) = -4.0_wp*inv_dx*N_aa_now
                    lgs%a_index(k) = nc

                    nc = 2*lgs%ij2n(ip1,jm1)
                        ! next nc (column counter), for uy(ip1,jm1)
                    k  = k+1
                    lgs%a_value(k) = -2.0_wp*inv_dy*N_aa_now
                    lgs%a_index(k) = nc

                    nc = 2*lgs%ij2n(ip1,j)-1
                        ! next nc (column counter), for ux(ip1,j)
                    k = k+1
                    lgs%a_value(k) = 4.0_wp*inv_dx*N_aa_now
                    lgs%a_index(k) = nc

                    nc = 2*lgs%ij2n(ip1,j)
                        ! largest nc (column counter), for uy(ip1,j)
                    k  = k+1
                    lgs%a_value(k) = 2.0_wp*inv_dy*N_aa_now
                    lgs%a_index(k) = nc

                    ! Assign matrix values
                    lgs%b_value(nr) = taul_int_acx(i,j) 
                    lgs%x_value(nr) = ux(i,j)
                
                end if 

            else
                ! === Inner SSA solution === 

                ! -- vx terms -- 

                nc = 2*lgs%ij2n(i,j)-1          ! column counter for ux(i,j)
                k = k+1
                lgs%a_value(k) = -4.0_wp*inv_dxdx*(N_aa(ip1,j)+N_aa(i,j)) &
                                 -1.0_wp*inv_dydy*(N_ab(i,j)+N_ab(i,jm1)) &
                                 -beta_acx(i,j)
                lgs%a_index(k) = nc

                nc = 2*lgs%ij2n(ip1,j)-1        ! column counter for ux(ip1,j)
                k = k+1
                lgs%a_value(k) =  4.0_wp*inv_dxdx*N_aa(ip1,j)
                lgs%a_index(k) = nc

                nc = 2*lgs%ij2n(im1,j)-1        ! column counter for ux(im1,j)
                k = k+1
                lgs%a_value(k) =  4.0_wp*inv_dxdx*N_aa(i,j)
                lgs%a_index(k) = nc

                nc = 2*lgs%ij2n(i,jp1)-1        ! column counter for ux(i,jp1)
                k = k+1
                lgs%a_value(k) =  1.0_wp*inv_dydy*N_ab(i,j)
                lgs%a_index(k) = nc

                nc = 2*lgs%ij2n(i,jm1)-1        ! column counter for ux(i,jm1)
                k = k+1
                lgs%a_value(k) =  1.0_wp*inv_dydy*N_ab(i,jm1)
                lgs%a_index(k) = nc

                ! -- vy terms -- 
                
                nc = 2*lgs%ij2n(i,j)            ! column counter for uy(i,j)
                k = k+1
                lgs%a_value(k) = -2.0_wp*inv_dxdy*N_aa(i,j)     &
                                 -1.0_wp*inv_dxdy*N_ab(i,j)
                lgs%a_index(k) = nc

                nc = 2*lgs%ij2n(ip1,j)          ! column counter for uy(ip1,j)
                k = k+1
                lgs%a_value(k) =  2.0_wp*inv_dxdy*N_aa(ip1,j)   &
                                 +1.0_wp*inv_dxdy*N_ab(i,j)
                lgs%a_index(k) = nc
                
                nc = 2*lgs%ij2n(ip1,jm1)        ! column counter for uy(ip1,jm1)
                k = k+1
                lgs%a_value(k) = -2.0_wp*inv_dxdy*N_aa(ip1,j)   &
                                 -1.0_wp*inv_dxdy*N_ab(i,jm1)
                lgs%a_index(k) = nc
                
                nc = 2*lgs%ij2n(i,jm1)          ! column counter for uy(i,jm1)
                k = k+1
                lgs%a_value(k) =  2.0_wp*inv_dxdy*N_aa(i,j)   &
                                 +1.0_wp*inv_dxdy*N_ab(i,jm1)
                lgs%a_index(k) = nc
                

                lgs%b_value(nr) = taud_acx(i,j)
                lgs%x_value(nr) = ux(i,j)

            end if

            lgs%a_ptr(nr+1) = k+1   ! row is completed, store index to next row

            ! ------ Equations for uy ---------------------------

            nr = n+1   ! row counter

            ! == Treat special cases first ==
            
            if (ssa_mask_acy(i,j) .eq. 0) then
                ! SSA not active here, velocity set to zero

                k = k+1
                lgs%a_value(k)  = 1.0_wp   ! diagonal element only
                lgs%a_index(k)  = nr

                lgs%b_value(nr) = 0.0_wp
                lgs%x_value(nr) = 0.0_wp

            else if (ssa_mask_acy(i,j) .eq. -1) then 
                ! Assign prescribed boundary velocity to this point
                ! (eg for prescribed velocity corresponding to analytical grounding line flux)

                k = k+1
                lgs%a_value(k)  = 1.0   ! diagonal element only
                lgs%a_index(k)  = nr

                lgs%b_value(nr) = uy(i,j)
                lgs%x_value(nr) = uy(i,j)
            
            else if (j .eq. 1 .and. trim(bcs(4)) .ne. "periodic") then 
                ! lower boundary 

                select case(trim(bcs(4)))

                    case("no-slip")

                        k = k+1
                        lgs%a_value(k)  = 1.0_wp        ! diagonal element only
                        lgs%a_index(k)  = nr

                        lgs%b_value(nr) = 0.0_wp
                        lgs%x_value(nr) = 0.0_wp

                    case("free-slip")

                        nc = 2*lgs%ij2n(i,j)            ! column counter for uy(i,j)
                        k = k+1
                        lgs%a_value(k) =  1.0_wp
                        lgs%a_index(k) = nc

                        nc = 2*lgs%ij2n(i,jp1)            ! column counter for uy(i,jp1)
                        k = k+1
                        lgs%a_value(k) = -1.0_wp
                        lgs%a_index(k) = nc

                        lgs%b_value(nr) = 0.0_wp
                        lgs%x_value(nr) = uy(i,j)

                end select 

            else if (j .eq. ny .and. trim(bcs(2)) .ne. "periodic") then 
                ! Upper boundary 

                select case(trim(bcs(2)))

                    case("no-slip")

                        k = k+1
                        lgs%a_value(k)  = 1.0_wp        ! diagonal element only
                        lgs%a_index(k)  = nr

                        lgs%b_value(nr) = 0.0_wp
                        lgs%x_value(nr) = 0.0_wp

                    case("free-slip")

                        nc = 2*lgs%ij2n(i,j)            ! column counter for uy(i,j)
                        k = k+1
                        lgs%a_value(k) =  1.0_wp
                        lgs%a_index(k) = nc

                        nc = 2*lgs%ij2n(i,ny-1)         ! column counter for uy(i,ny-1)
                        k = k+1
                        lgs%a_value(k) = -1.0_wp
                        lgs%a_index(k) = nc

                        lgs%b_value(nr) = 0.0_wp
                        lgs%x_value(nr) = uy(i,j)

                end select 

            else if (i .eq. 1 .and. trim(bcs(3)) .ne. "periodic") then 
                ! Left boundary 

                select case(trim(bcs(3)))

                    case("no-slip")

                        k = k+1
                        lgs%a_value(k)  = 1.0_wp   ! diagonal element only
                        lgs%a_index(k)  = nr

                        lgs%b_value(nr) = 0.0_wp
                        lgs%x_value(nr) = 0.0_wp

                    case("free-slip")

                        nc = 2*lgs%ij2n(i,j)            ! column counter for uy(i,j)
                        k = k+1
                        lgs%a_value(k) =  1.0_wp
                        lgs%a_index(k) = nc

                        nc = 2*lgs%ij2n(ip1,j)          ! column counter for uy(ip1,j)
                        k = k+1
                        lgs%a_value(k) = -1.0_wp
                        lgs%a_index(k) = nc

                        lgs%b_value(nr) = 0.0_wp
                        lgs%x_value(nr) = uy(i,j)

                end select 
                
            else if (i .eq. nx .and. trim(bcs(1)) .ne. "periodic") then 
                ! Right boundary 

                select case(trim(bcs(1)))

                    case("no-slip")

                        k = k+1
                        lgs%a_value(k)  = 1.0_wp   ! diagonal element only
                        lgs%a_index(k)  = nr

                        lgs%b_value(nr) = 0.0_wp
                        lgs%x_value(nr) = 0.0_wp

                    case("free-slip")

                        nc = 2*lgs%ij2n(i,j)            ! column counter for uy(i,j)
                        k = k+1
                        lgs%a_value(k) =  1.0_wp
                        lgs%a_index(k) = nc

                        nc = 2*lgs%ij2n(nx-1,j)         ! column counter for uy(nx-1,j)
                        k = k+1
                        lgs%a_value(k) = -1.0_wp
                        lgs%a_index(k) = nc

                        lgs%b_value(nr) = 0.0_wp
                        lgs%x_value(nr) = uy(i,j)

                end select 

            else if (ssa_mask_acy(i,j) .eq. 3) then 
                ! Lateral boundary condition should be applied here 

                if (f_ice(i,j) .eq. 1.0 .and. f_ice(i,jp1) .lt. 1.0) then 
                    ! === Case 1: ice-free to the top ===

                    N_aa_now = N_aa(i,j)

                    nc = 2*lgs%ij2n(im1,j)-1
                        ! smallest nc (column counter), for ux(im1,j)
                    k = k+1
                    lgs%a_value(k) = -2.0_wp*inv_dx*N_aa_now
                    lgs%a_index(k) = nc

                    nc = 2*lgs%ij2n(i,jm1)
                        ! next nc (column counter), for uy(i,jm1)
                    k = k+1
                    lgs%a_value(k) = -4.0_wp*inv_dy*N_aa_now
                    lgs%a_index(k) = nc

                    nc = 2*lgs%ij2n(i,j)-1
                        ! next nc (column counter), for ux(i,j)
                    k = k+1
                    lgs%a_value(k) = 2.0_wp*inv_dx*N_aa_now
                    lgs%a_index(k) = nc

                    nc = 2*lgs%ij2n(i,j)
                        ! next nc (column counter), for uy(i,j)
                    k = k+1
                    lgs%a_value(k) = 4.0_wp*inv_dy*N_aa_now
                    lgs%a_index(k) = nc

                    ! Assign matrix values
                    lgs%b_value(nr) = taul_int_acy(i,j)
                    lgs%x_value(nr) = uy(i,j)
                    
                else
                    ! === Case 2: ice-free to the bottom ===
                    
                    N_aa_now = N_aa(i,jp1)

                    nc = 2*lgs%ij2n(im1,jp1)-1
                        ! next nc (column counter), for ux(im1,jp1)
                    k = k+1
                    lgs%a_value(k) = -2.0_wp*inv_dx*N_aa_now
                    lgs%a_index(k) = nc

                    nc = 2*lgs%ij2n(i,j)
                        ! next nc (column counter), for uy(i,j)
                    k = k+1
                    lgs%a_value(k) = -4.0_wp*inv_dy*N_aa_now
                    lgs%a_index(k) = nc

                    nc = 2*lgs%ij2n(i,jp1)-1
                        ! next nc (column counter), for ux(i,jp1)
                    k = k+1
                    lgs%a_value(k) = 2.0_wp*inv_dx*N_aa_now
                    lgs%a_index(k) = nc

                    nc = 2*lgs%ij2n(i,jp1)
                        ! next nc (column counter), for uy(i,jp1)
                    k = k+1
                    lgs%a_value(k) = 4.0_wp*inv_dy*N_aa_now
                    lgs%a_index(k) = nc

                    ! Assign matrix values
                    lgs%b_value(nr) = taul_int_acy(i,j)
                    lgs%x_value(nr) = uy(i,j)
         
                end if 

            else
                ! === Inner SSA solution === 

                ! -- vy terms -- 

                nc = 2*lgs%ij2n(i,j)        ! column counter for uy(i,j)
                k = k+1
                lgs%a_value(k) = -4.0_wp*inv_dydy*(N_aa(i,jp1)+N_aa(i,j))   &
                                 -1.0_wp*inv_dxdx*(N_ab(i,j)+N_ab(im1,j))   &
                                 -beta_acy(i,j)
                lgs%a_index(k) = nc

                nc = 2*lgs%ij2n(i,jp1)      ! column counter for uy(i,jp1)
                k = k+1
                lgs%a_value(k) =  4.0_wp*inv_dydy*N_aa(i,jp1)
                lgs%a_index(k) = nc

                nc = 2*lgs%ij2n(i,jm1)      ! column counter for uy(i,jm1)
                k = k+1
                lgs%a_value(k) =  4.0_wp*inv_dydy*N_aa(i,j)
                lgs%a_index(k) = nc
                
                nc = 2*lgs%ij2n(ip1,j)      ! column counter for uy(ip1,j)
                k = k+1
                lgs%a_value(k) =  1.0_wp*inv_dxdx*N_ab(i,j)
                lgs%a_index(k) = nc
                
                nc = 2*lgs%ij2n(im1,j)      ! column counter for uy(im1,j)
                k = k+1
                lgs%a_value(k) =  1.0_wp*inv_dxdx*N_ab(im1,j)
                lgs%a_index(k) = nc
                
                ! -- vx terms -- 

                nc = 2*lgs%ij2n(i,j)-1      ! column counter for ux(i,j)
                k = k+1
                lgs%a_value(k) = -2.0_wp*inv_dxdy*N_aa(i,j)     &
                                 -1.0_wp*inv_dxdy*N_ab(i,j)
                lgs%a_index(k) = nc

                nc = 2*lgs%ij2n(i,jp1)-1    ! column counter for ux(i,jp1)
                k = k+1
                lgs%a_value(k) =  2.0_wp*inv_dxdy*N_aa(i,jp1)     &
                                 +1.0_wp*inv_dxdy*N_ab(i,j)
                lgs%a_index(k) = nc

                nc = 2*lgs%ij2n(im1,jp1)-1  ! column counter for ux(im1,jp1)
                k = k+1
                lgs%a_value(k) = -2.0_wp*inv_dxdy*N_aa(i,jp1)     &
                                 -1.0_wp*inv_dxdy*N_ab(im1,j)
                lgs%a_index(k) = nc

                nc = 2*lgs%ij2n(im1,j)-1  ! column counter for ux(im1,j)
                k = k+1
                lgs%a_value(k) =  2.0_wp*inv_dxdy*N_aa(i,j)     &
                                 +1.0_wp*inv_dxdy*N_ab(im1,j)
                lgs%a_index(k) = nc

                lgs%b_value(nr) = taud_acy(i,j)
                lgs%x_value(nr) = uy(i,j)

            end if

            lgs%a_ptr(nr+1) = k+1   ! row is completed, store index to next row

        end do

        ! Done: A, x and b matrices in Ax=b have been populated 
        ! and stored in lgs object. 

        return

    end subroutine linear_solver_matrix_ssa_ac_csr_2D

    subroutine set_ssa_masks(ssa_mask_acx,ssa_mask_acy,mask_frnt,H_ice,f_ice,f_grnd,use_ssa,lateral_bc)
        ! Define where ssa calculations should be performed
        ! Note: could be binary, but perhaps also distinguish 
        ! grounding line/zone to use this mask for later gl flux corrections
        ! mask = -1: no ssa calculated, velocity imposed/unchanged
        ! mask = 0: no ssa calculated, velocity set to zero
        ! mask = 1: shelfy-stream ssa calculated 
        ! mask = 2: shelf ssa calculated 
        ! mask = 3: ssa lateral boundary condition applied
        ! mask = 4: ssa lateral boundary, but treated as inner ssa

        implicit none 
        
        integer,  intent(OUT) :: ssa_mask_acx(:,:) 
        integer,  intent(OUT) :: ssa_mask_acy(:,:)
        integer,  intent(IN)  :: mask_frnt(:,:)
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: f_grnd(:,:)
        logical,  intent(IN)  :: use_ssa       ! SSA is actually active now? 
        character(len=*), intent(IN) :: lateral_bc 

        ! Local variables
        integer    :: i, j, nx, ny
        integer    :: im1, ip1, jm1, jp1
        real(wp) :: H_acx, H_acy
        
        real(wp), allocatable :: mask_frnt_dyn(:,:)

        integer, parameter :: val_ice_free  = -1 
        integer, parameter :: val_flt       = 1
        integer, parameter :: val_marine    = 2
        integer, parameter :: val_grnd      = 3
        integer, parameter :: val_disabled  = 5 

        nx = size(H_ice,1)
        ny = size(H_ice,2)
        
        allocate(mask_frnt_dyn(nx,ny))

        ! Initially no active ssa points, all velocities set to zero
        ssa_mask_acx = 0
        ssa_mask_acy = 0
        
        if (use_ssa) then 

            ! Step 1: define mask_frnt for dynamics, that disables fronts as needed 

            mask_frnt_dyn = mask_frnt

                    
            ! So far all margins have been diagnosed (marine and grounded on land)
            ! Disable some regions depending on choice above. 
            select case(trim(lateral_bc))
                ! Apply the lateral boundary condition to what? 

                case("none")
                    ! Do not apply lateral bc anywhere. Ie, disable front detection.
                    ! Treat all ice points in the domain as 'inner ssa' points.

                    do j = 1, ny
                    do i = 1, nx
                    
                        if (mask_frnt(i,j) .gt. 0) mask_frnt_dyn(i,j) = val_disabled

                    end do
                    end do

                case("floating","float","slab","slab-ext")
                    ! Only apply lateral bc to floating ice fronts.
                    ! Ie, disable detection of all grounded fronts for now.
                    ! Model is generally more stable this way.
                    ! This method is also used for the 'infinite slab' approach,
                    ! where a thin ice shelf is extended everywhere over the domain. 
                    
                    do j = 1, ny
                    do i = 1, nx
                    
                        if ( mask_frnt(i,j) .eq. val_grnd .or. &
                             mask_frnt(i,j) .eq. val_marine ) mask_frnt_dyn(i,j) = val_disabled

                    end do
                    end do

                case("marine")
                    ! Only apply lateral bc to floating ice fronts and
                    ! and grounded marine fronts. Disable detection 
                    ! of ice fronts grounded above sea level.
                    
                    do j = 1, ny
                    do i = 1, nx
                    
                        if ( mask_frnt(i,j) .eq. val_grnd ) mask_frnt_dyn(i,j) = val_disabled

                    end do
                    end do

                case("all")
                    ! Apply lateral bc to all ice-sheet fronts. 

                    ! Do nothing - all fronts have been accurately diagnosed. 

                case DEFAULT
                    
                    write(io_unit_err,*) "set_ssa_masks:: error: ssa_lat_bc parameter value not recognized."
                    write(io_unit_err,*) "ydyn.ssa_lat_bc = ", lateral_bc
                    stop 

            end select
               
                    
            ! Step 2: define ssa solver masks

            do j = 1, ny
            do i = 1, nx

                ! Get neighbor indices
                im1 = max(i-1,1) 
                ip1 = min(i+1,nx) 
                jm1 = max(j-1,1) 
                jp1 = min(j+1,ny)


                ! == x-direction ===

                if (f_ice(i,j) .eq. 1.0 .or. f_ice(ip1,j) .eq. 1.0) then
                
                    ! Current ac-node is border of an ice covered cell in x-direction
                    
                    if (f_grnd(i,j) .gt. 0.0 .or. f_grnd(ip1,j) .gt. 0.0) then 
                        ! Grounded ice or grounding line (ie, shelfy-stream)
                        ssa_mask_acx(i,j) = 1
                    else 
                        ! Shelf ice 
                        ssa_mask_acx(i,j) = 2
                    end if 

                    ! Check for special case of floating ice next to ice-free land,
                    ! then set ssa mask to zero (ie, set velocity to zero)
                    if (ssa_mask_acx(i,j) .eq. 2) then 

                        if ( f_grnd(i,j) .eq. 0.0 .and. &
                               (f_grnd(ip1,j) .gt. 0.0 .and. H_ice(ip1,j) .eq. 0.0) ) then 

                            ssa_mask_acx(i,j) = 0

                        else if ( (f_grnd(i,j) .gt. 0.0 .and. H_ice(i,j) .eq. 0.0) .and. &
                                    f_grnd(ip1,j) .eq. 0.0 ) then

                            ssa_mask_acx(i,j) = 0 

                        end if 

                    end if 

                end if

                ! Overwrite above if this point should be treated via lateral boundary conditions
                if ( (mask_frnt_dyn(i,j) .gt. 0 .and. mask_frnt_dyn(ip1,j) .lt. 0) .or. &
                     (mask_frnt_dyn(i,j) .lt. 0 .and. mask_frnt_dyn(ip1,j) .gt. 0) ) then 
                    ! Lateral boundary point 

                    ssa_mask_acx(i,j) = 3 

                end if 

                ! Overwrite again if this front should be deactivated 
                if ( (mask_frnt_dyn(i,j) .eq. 5 .and. mask_frnt_dyn(ip1,j) .lt. 0) .or. &
                     (mask_frnt_dyn(i,j) .lt. 0 .and. mask_frnt_dyn(ip1,j) .eq. 5) ) then 
                    ! Deactivated lateral boundary point 

                    ssa_mask_acx(i,j) = 4 

                end if 

                ! == y-direction ===

                if (f_ice(i,j) .eq. 1.0 .or. f_ice(i,jp1) .eq. 1.0) then
                
                    ! Current ac-node is border of an ice covered cell in x-direction
                    
                    if (f_grnd(i,j) .gt. 0.0 .or. f_grnd(i,jp1) .gt. 0.0) then 
                        ! Grounded ice or grounding line (ie, shelfy-stream)
                        ssa_mask_acy(i,j) = 1
                    else 
                        ! Shelf ice 
                        ssa_mask_acy(i,j) = 2
                    end if 

                    ! Check for special case of floating ice next to ice-free land,
                    ! then set ssa mask to zero (ie, set velocity to zero)
                    if (ssa_mask_acy(i,j) .eq. 2) then 

                        if ( f_grnd(i,j) .eq. 0.0 .and. &
                               (f_grnd(i,jp1) .gt. 0.0 .and. H_ice(i,jp1) .eq. 0.0) ) then 

                            ssa_mask_acy(i,j) = 0

                        else if ( (f_grnd(i,j) .gt. 0.0 .and. H_ice(i,j) .eq. 0.0) .and. &
                                    f_grnd(i,jp1) .eq. 0.0 ) then

                            ssa_mask_acy(i,j) = 0 

                        end if 

                    end if 
                    
                end if

                ! Overwrite above if this point should be treated via lateral boundary conditions
                if ( (mask_frnt_dyn(i,j) .gt. 0 .and. mask_frnt_dyn(i,jp1) .lt. 0) .or. &
                     (mask_frnt_dyn(i,j) .lt. 0 .and. mask_frnt_dyn(i,jp1) .gt. 0) ) then 
                    ! Lateral boundary point 

                    ssa_mask_acy(i,j) = 3 

                end if 

                ! Overwrite again if this front should be deactivated 
                if ( (mask_frnt_dyn(i,j) .eq. 5 .and. mask_frnt_dyn(i,jp1) .lt. 0) .or. &
                     (mask_frnt_dyn(i,j) .lt. 0 .and. mask_frnt_dyn(i,jp1) .eq. 5) ) then 
                    ! Deactivated lateral boundary point 

                    ssa_mask_acy(i,j) = 4 

                end if 

            end do 
            end do

        end if 

        return
        
    end subroutine set_ssa_masks
    
    subroutine update_ssa_mask_convergence(ssa_mask_acx,ssa_mask_acy,err_x,err_y,err_lim)
        ! Update grounded ice ssa_masks, by prescribing vel at points that have
        ! already converged well.

        implicit none 

        integer, intent(INOUT) :: ssa_mask_acx(:,:) 
        integer, intent(INOUT) :: ssa_mask_acy(:,:) 
        real(wp), intent(IN) :: err_x(:,:) 
        real(wp), intent(IN) :: err_y(:,:) 
        real(wp), intent(IN) :: err_lim 

        ! Local variables 
        integer :: i, j, nx, ny 

        nx = size(ssa_mask_acx,1)
        ny = size(ssa_mask_acx,2) 

        ! Initially set candidate 'converged' points to -2
        where (ssa_mask_acx .eq. 1 .and. err_x .lt. err_lim)
            ssa_mask_acx = -2 
        end where 

        where (ssa_mask_acy .eq. 1 .and. err_y .lt. err_lim)
            ssa_mask_acy = -2 
        end where 
        
        ! Fill in neighbors of points that are still ssa (mask=1) to keep things clean 
        do j = 2, ny-1 
        do i = 2, nx-1 
            
            ! acx
            if (ssa_mask_acx(i,j) .eq. 1) then
                where (ssa_mask_acx(i-1:i+1,j-1:j+1) .eq. -2) ssa_mask_acx(i-1:i+1,j-1:j+1) = 1
            end if 

            ! acy 
            if (ssa_mask_acy(i,j) .eq. 1) then
                where (ssa_mask_acy(i-1:i+1,j-1:j+1) .eq. -2) ssa_mask_acy(i-1:i+1,j-1:j+1) = 1
            end if 
            
        end do 
        end do 

        ! Finally, replace temporary -2 values with -1 to prescribe ssa vel here 
        where (ssa_mask_acx .eq. -2) ssa_mask_acx = -1 
        where (ssa_mask_acy .eq. -2) ssa_mask_acy = -1 
        
        return 

    end subroutine update_ssa_mask_convergence

! === INTERNAL ROUTINES ==== 

    subroutine stagger_visc_aa_ab(visc_ab,visc,H_ice,f_ice,boundaries)

        implicit none 

        real(wp), intent(OUT) :: visc_ab(:,:) 
        real(wp), intent(IN)  :: visc(:,:) 
        real(wp), intent(IN)  :: H_ice(:,:) 
        real(wp), intent(IN)  :: f_ice(:,:) 
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer :: i, j, k
        integer :: im1, ip1, jm1, jp1 
        integer :: nx, ny 

        nx = size(visc,1)
        ny = size(visc,2)

        ! Initialisation
        visc_ab = 0.0_wp 

        ! Stagger viscosity only using contributions from neighbors that have ice  
        do i = 1, nx 
        do j = 1, ny 

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

            visc_ab(i,j) = 0.0_wp
            k=0

            if (f_ice(i,j) .eq. 1.0) then
                k = k+1                              ! floating or grounded ice
                visc_ab(i,j) = visc_ab(i,j) + visc(i,j)
            end if

            if (f_ice(ip1,j) .eq. 1.0) then
                k = k+1                                  ! floating or grounded ice
                visc_ab(i,j) = visc_ab(i,j) + visc(ip1,j)
            end if

            if (f_ice(i,jp1) .eq. 1.0) then
                k = k+1                                  ! floating or grounded ice
                visc_ab(i,j) = visc_ab(i,j) + visc(i,jp1)
            end if

            if (f_ice(ip1,jp1) .eq. 1.0) then
                k = k+1                                      ! floating or grounded ice
                visc_ab(i,j) = visc_ab(i,j) + visc(ip1,jp1)
            end if

            if (k .gt. 0) visc_ab(i,j) = visc_ab(i,j)/real(k,wp)

        end do
        end do

        return 

    end subroutine stagger_visc_aa_ab

    elemental subroutine limit_vel(u,u_lim)
        ! Apply a velocity limit (for stability)

        implicit none 

        real(wp), intent(INOUT) :: u  
        real(wp), intent(IN)    :: u_lim

        real(wp), parameter :: tol = TOL_UNDERFLOW

        u = min(u, u_lim)
        u = max(u,-u_lim)

        ! Also avoid underflow errors 
        if (abs(u) .lt. tol) u = 0.0 

        return 

    end subroutine limit_vel

    ! === DIAGNOSTIC OUTPUT ROUTINES ===

    subroutine ssa_diagnostics_write_init(filename,nx,ny,time_init)

        implicit none 

        character(len=*),  intent(IN) :: filename 
        integer,           intent(IN) :: nx 
        integer,           intent(IN) :: ny
        real(wp),          intent(IN) :: time_init

        ! Initialize netcdf file and dimensions
        call nc_create(filename)
        call nc_write_dim(filename,"xc",     x=0.0_wp,dx=1.0_wp,nx=nx,units="gridpoints")
        call nc_write_dim(filename,"yc",     x=0.0_wp,dx=1.0_wp,nx=ny,units="gridpoints")
        call nc_write_dim(filename,"time",   x=time_init,dx=1.0_wp,nx=1,units="iter",unlimited=.TRUE.)

        return

    end subroutine ssa_diagnostics_write_init

    subroutine ssa_diagnostics_write_step(filename,ux,uy,L2_norm,beta_acx,beta_acy,visc_int, &
                                        ssa_mask_acx,ssa_mask_acy,ssa_err_acx,ssa_err_acy,H_ice,f_ice,taud_acx,taud_acy, &
                                     taul_int_acx,taul_int_acy,H_grnd,z_sl,z_bed,z_srf,ux_prev,uy_prev,time)

        implicit none 
        
        character(len=*),  intent(IN) :: filename
        real(wp), intent(IN) :: ux(:,:) 
        real(wp), intent(IN) :: uy(:,:) 
        real(wp), intent(IN) :: L2_norm
        real(wp), intent(IN) :: beta_acx(:,:) 
        real(wp), intent(IN) :: beta_acy(:,:) 
        real(wp), intent(IN) :: visc_int(:,:) 
        integer,  intent(IN) :: ssa_mask_acx(:,:) 
        integer,  intent(IN) :: ssa_mask_acy(:,:) 
        real(wp), intent(IN) :: ssa_err_acx(:,:) 
        real(wp), intent(IN) :: ssa_err_acy(:,:) 
        real(wp), intent(IN) :: H_ice(:,:) 
        real(wp), intent(IN) :: f_ice(:,:) 
        real(wp), intent(IN) :: taud_acx(:,:) 
        real(wp), intent(IN) :: taud_acy(:,:) 
        real(wp), intent(IN) :: taul_int_acx(:,:) 
        real(wp), intent(IN) :: taul_int_acy(:,:) 
        real(wp), intent(IN) :: H_grnd(:,:) 
        real(wp), intent(IN) :: z_sl(:,:) 
        real(wp), intent(IN) :: z_bed(:,:) 
        real(wp), intent(IN) :: z_srf(:,:)
        real(wp), intent(IN) :: ux_prev(:,:) 
        real(wp), intent(IN) :: uy_prev(:,:) 
        real(wp), intent(IN) :: time

        ! Local variables
        integer  :: ncid, n, i, j, nx, ny  
        real(wp) :: time_prev 

        nx = size(ux,1)
        ny = size(ux,2) 

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_size(filename,"time",ncid)
        call nc_read(filename,"time",time_prev,start=[n],count=[1],ncid=ncid) 
        if (abs(time-time_prev).gt.1e-5) n = n+1 

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        ! Write the variables 

        call nc_write(filename,"ux",ux,units="m/yr",long_name="ssa velocity (acx)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy",uy,units="m/yr",long_name="ssa velocity (acy)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"ux_diff",ux-ux_prev,units="m/yr",long_name="ssa velocity difference (acx)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy_diff",uy-uy_prev,units="m/yr",long_name="ssa velocity difference (acy)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"L2_norm",L2_norm,dim1="time",start=[n],count=[1],ncid=ncid)

        call nc_write(filename,"beta_acx",beta_acx,units="Pa yr m^-1",long_name="Dragging coefficient (acx)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"beta_acy",beta_acy,units="Pa yr m^-1",long_name="Dragging coefficient (acy)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"visc_int",visc_int,units="Pa yr",long_name="Vertically integrated effective viscosity", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"ssa_mask_acx",ssa_mask_acx,units="1",long_name="SSA mask (acx)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"ssa_mask_acy",ssa_mask_acy,units="1",long_name="SSA mask (acy)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"ssa_err_acx",ssa_err_acx,units="1",long_name="SSA err (acx)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"ssa_err_acy",ssa_err_acy,units="1",long_name="SSA err (acy)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"H_ice",H_ice,units="m",long_name="Ice thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_ice",f_ice,units="1",long_name="Ice-covered fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"taud_acx",taud_acx,units="Pa",long_name="Driving stress (acx)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taud_acy",taud_acy,units="Pa",long_name="Driving stress (acy)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"taul_int_acx",taul_int_acx,units="Pa",long_name="Vertically integrated lateral stress (acx)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taul_int_acy",taul_int_acy,units="Pa",long_name="Vertically integrated lateral stress (acy)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"H_grnd",H_grnd,units="m",long_name="Ice thickness overburden", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_sl",z_sl,units="m",long_name="Sea level", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_bed",z_bed,units="m",long_name="Bedrock elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf",z_srf,units="m",long_name="Surface elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine ssa_diagnostics_write_step

! ==== GENERAL INTERFACE TO MATRIX SOLVE ROUTINES =====
    
    subroutine linear_solver_matrix_solve(lgs,lis_settings)
        ! Use LIS to solve matrix equation Ax=b. Take
        ! predefined A, x and b arrays from the lgs object,
        ! and pass them to lis-specific variables, solve
        ! for new x and return to lgs object.

        implicit none 

        type(linear_solver_class), intent(INOUT) :: lgs 
        character(len=*), intent(IN) :: lis_settings        ! LIS solver settings

! ==== PETSC SPECIFIC CODE =====
#ifdef USEPETSC
    
        ! Solve matrix equation with PETSc
        call linear_solver_matrix_solve_petsc(lgs,rtol=1e-6_wp,atol=1e-6_wp,maxits=1000)

#else
        
        ! Solve matrix equation with LIS
        call linear_solver_matrix_solve_lis(lgs,lis_settings)
#endif

        return

    end subroutine linear_solver_matrix_solve

! ==== LIS SPECIFIC CODE =====

    subroutine linear_solver_matrix_solve_lis(lgs,lis_settings)
        ! Use LIS to solve matrix equation Ax=b. Take
        ! predefined A, x and b arrays from the lgs object,
        ! and pass them to lis-specific variables, solve
        ! for new x and return to lgs object.

        implicit none 

        type(linear_solver_class), intent(INOUT) :: lgs 
        character(len=*), intent(IN) :: lis_settings        ! LIS solver settings

        ! =========================================================
        ! LIS-specific variables 

        ! Include header for lis solver fortran interface
#include "lisf.h"
        
        LIS_INTEGER :: ierr
        LIS_INTEGER :: nr
        LIS_INTEGER :: nc
        LIS_INTEGER :: nmax
        LIS_INTEGER :: lin_iter
        LIS_REAL    :: residual 
        LIS_REAL    :: solver_time  
        LIS_MATRIX  :: lgs_a
        LIS_VECTOR  :: lgs_b, lgs_x
        LIS_SOLVER  :: solver

        LIS_INTEGER :: lgs_a_index_now
        LIS_SCALAR  :: lgs_a_value_now
        LIS_SCALAR  :: lgs_b_value_now
        LIS_SCALAR  :: lgs_x_value_now
        LIS_SCALAR, allocatable :: lgs_x_value_out(:)

        ! =========================================================
        
        ! Store nmax in local LIS-specific variable
        nmax = lgs%nmax

        !-------- Settings for Lis --------
               
        call lis_initialize(ierr)           ! Important for parallel computing environments   
        call CHKERR(ierr)

        call lis_matrix_create(LIS_COMM_WORLD, lgs_a, ierr)
        call CHKERR(ierr)
        call lis_vector_create(LIS_COMM_WORLD, lgs_b, ierr)
        call lis_vector_create(LIS_COMM_WORLD, lgs_x, ierr)

        call lis_matrix_set_size(lgs_a, 0, nmax, ierr)
        call CHKERR(ierr)
        call lis_vector_set_size(lgs_b, 0, nmax, ierr)
        call lis_vector_set_size(lgs_x, 0, nmax, ierr)

        ! === Storage order: compressed sparse row (CSR) ===

        do nr=1, nmax

            do nc=lgs%a_ptr(nr), lgs%a_ptr(nr+1)-1

                ! Use temporary values with LIS data types for use with lis routines
                lgs_a_index_now = lgs%a_index(nc)
                lgs_a_value_now = lgs%a_value(nc)
                
                call lis_matrix_set_value(LIS_INS_VALUE, nr, lgs_a_index_now, &
                                                        lgs_a_value_now, lgs_a, ierr)
            end do

            ! Use temporary values with LIS data types for use with lis routines
            lgs_b_value_now = lgs%b_value(nr)
            lgs_x_value_now = lgs%x_value(nr) 

            call lis_vector_set_value(LIS_INS_VALUE, nr, lgs_b_value_now, lgs_b, ierr)
            call lis_vector_set_value(LIS_INS_VALUE, nr, lgs_x_value_now, lgs_x, ierr)

        end do 


        call lis_matrix_set_type(lgs_a, LIS_MATRIX_CSR, ierr)
        call lis_matrix_assemble(lgs_a, ierr)

        !-------- Solution of the system of linear equations with Lis --------

        call lis_solver_create(solver, ierr)

        ! ch_solver_set_option = '-i bicgsafe -p jacobi '// &
        !                         '-maxiter 100 -tol 1.0e-4 -initx_zeros false'
        call lis_solver_set_option(trim(lis_settings), solver, ierr)
        call CHKERR(ierr)

        call lis_solve(lgs_a, lgs_b, lgs_x, solver, ierr)
        call CHKERR(ierr)

        ! Get solver solution information
        call lis_solver_get_iter(solver, lin_iter, ierr)
        call lis_solver_get_time(solver,solver_time,ierr)
        
        ! Obtain the relative L2_norm == ||b-Ax|| / ||b||
        call lis_solver_get_residualnorm(solver,residual,ierr)

        ! Store in lgs object too
        lgs%L2_rel_norm = residual

        ! Print a summary
        !write(*,*) "solve_lis: [time (s), iter, L2_rel_norm] = ", solver_time, lin_iter, residual

        ! Gather x values in local array of lis-type
        allocate(lgs_x_value_out(nmax))
        lgs_x_value_out = 0.0_wp
        call lis_vector_gather(lgs_x, lgs_x_value_out, ierr)
        call CHKERR(ierr)
        
        ! Save to lgs object
        lgs%x_value = lgs_x_value_out
        
        ! Destroy all lis variables
        call lis_matrix_destroy(lgs_a, ierr)
        call CHKERR(ierr)

        call lis_vector_destroy(lgs_b, ierr)
        call lis_vector_destroy(lgs_x, ierr)
        call lis_solver_destroy(solver, ierr)
        call CHKERR(ierr)

        ! Finalize lis.
        call lis_finalize(ierr)           ! Important for parallel computing environments
        call CHKERR(ierr)

        return

    end subroutine linear_solver_matrix_solve_lis


! ==== PETSC SPECIFIC CODE =====
#ifdef USEPETSC

    subroutine linear_solver_matrix_solve_petsc(lgs,rtol,atol,maxits)
        ! Use PETSC Krylov solver to solve matrix equation Ax=b. Take
        ! predefined A, x and b arrays from the lgs object,
        ! and pass them to lis-specific variables, solve
        ! for new x and return to lgs object.

#include <petsc/finclude/petscksp.h>

        use petscksp
        use mpi

        implicit none 

        type(linear_solver_class), intent(INOUT) :: lgs 
        real(wp) :: rtol                    ! Input tolerance parameter
        real(wp) :: atol                    ! Input tolerance parameter
        integer  :: maxits                  ! Input iteration parameter 

        ! Local variables
        
        type(PetscErrorCode) :: perr
        type(PetscInt)       :: lin_iter
        type(tMat)           :: A
        type(tVec)           :: b
        type(tVec)           :: x
        type(tKSP)           :: KSP_solver
        type(PetscReal)      :: PETSc_rtol
        type(PetscReal)      :: PETSc_atol
        type(PetscInt)       :: PETSc_maxits
        
        ! Initialise PETSc MPI stuff
        call initialise_petsc()
        
        ! Convert the CSR matrix from the lgs storage structure to a PETSc matrix
        call convert_lgs_to_petsc( A, b, x, lgs)
        
        ! Set up the KSP solver
        call KSPcreate( PETSC_COMM_WORLD, KSP_solver, perr)
        call handle_error( 'KSPcreate', perr)
      
        ! Set operators. Here the matrix that defines the linear system
        ! also serves as the preconditioning matrix.
        call KSPSetOperators( KSP_solver, A, A, perr)
        call handle_error( 'KSPSetOperators', perr)

        ! Iterative solver tolerances
        PETSc_rtol   = rtol         ! The relative convergence tolerance; relative decrease in the (possibly preconditioned) residual norm
        PETSc_atol   = atol         ! The absolute convergence tolerance; absolute size of the (possibly preconditioned) residual norm 
        PETSc_maxits = maxits       ! Maximum allowed iterations
        call KSPSetTolerances( KSP_solver, PETSc_rtol, PETSc_atol, PETSC_DEFAULT_REAL, PETSc_maxits, perr)
        call handle_error( 'KSPSetTolerances', perr)
        
        ! To start from nonzero guess
        !call KSPSetInitialGuessNonzero( KSP_solver,PetscBool flg);

        ! Set runtime options, e.g.,
        !     -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
        ! These options will override those specified above as long as
        ! KSPSetFromOptions() is called _after_ any other customization routines.
        call KSPSetFromOptions( KSP_solver, perr)
        call handle_error( 'KSPSetFromOptions', perr)
        
        ! Solve the linear system
        call KSPSolve( KSP_solver, b, x, perr)
        call handle_error( 'KSPSolve', perr)
        
        ! Check convergence information
        !call KSPSetConvergenceTest( KSP_solver, PetscErrorCode (*test)(KSP ksp,PetscInt it,PetscReal rnorm, KSPConvergedReason *reason,void *ctx),void *ctx,PetscErrorCode (*destroy)(void *ctx));

        ! Find out how many iterations it took
        call KSPGetIterationNumber( KSP_solver, lin_iter, perr)
        call handle_error( 'KSPGetIterationNumber', perr)
        write(*,*) "  PETSc solved Ax=b in ", lin_iter, " iterations"
        
        ! Get the solution back to the lgs storage structure
        call convert_petsc_solution_to_lgs( lgs, x)
        call sync
        
        ! Clean up after yourself
        call KSPDestroy( KSP_solver, perr)
        call handle_error( 'KSPDestroy', perr)
        call VecDestroy( x, perr)
        call handle_error( 'VecDestroy', perr)
        call VecDestroy( b, perr)
        call handle_error( 'VecDestroy', perr)
        call MatDestroy( A, perr)
        call handle_error( 'MatDestroy', perr)
        
        ! Finalise PETSc MPI stuff
        call finalise_petsc()

        return

        contains
            ! Petsc specific routines 
            
            ! subroutine parse_solver_settings_petsc(rtol,atol,dtol,maxits,settings)
            !     -ksp_rtol <rtol> -ksp_atol <atol> -ksp_divtol <dtol> -ksp_max_it <its>

            !     implicit none

            !     real(wp), intent(OUT) :: rtol 
            !     real(wp), intent(OUT) :: atol 
            !     real(wp), intent(OUT) :: dtol 
            !     integer,  intent(OUT) :: matits 
            !     character(len=*), intent(IN) :: settings 

            !     return

            ! end subroutine parse_solver_settings_petsc

            subroutine convert_lgs_to_petsc(A,b,x,lgs)
                ! Convert lgs storage object to PETSc matrices A, b and x.
                ! Adapted from IMAU-ICE petsc_module.F90, which was
                ! mostly copied from: https://petsc.org/release/documentation/manual/getting_started/#parallel-programming
                  
                implicit none
                
                ! In- and output variables:
                type(tMat), intent(INOUT) :: A
                type(tVec), intent(INOUT) :: b
                type(tVec), intent(INOUT) :: x
                type(linear_solver_class), intent(IN) :: lgs
                
                ! Local variables
                type(PetscErrorCode) :: perr
                integer              :: k, i, j, istart, iend
                type(PetscReal)      :: v
            
                ! == Matrix A ==
                ! ==============
            
                ! Initialise the matrix object
                call MatCreate( PETSC_COMM_WORLD, A, perr)
                call handle_error( 'MatCreate', perr)
                
                ! Set the matrix type to parallel (MPI) Aij
                call MatSetType( A, 'mpiaij', perr)
                call handle_error( 'MatSetType', perr)
                
                ! Set the size, let PETSc automatically determine parallelisation domains
                !call MatSetSizes( A, PETSC_DECIDE, PETSC_DECIDE, CSR%m, CSR%n, perr)
                call MatSetSizes( A, PETSC_DECIDE, PETSC_DECIDE, lgs%nmax, lgs%nmax, perr)
                call handle_error( 'MatSetSizes', perr)
                
                ! Not entirely sure what this one does, but apparently it's really important
                call MatSetFromOptions( A, perr)
                call handle_error( 'MatSetFromOptions', perr)
                
                ! Tell PETSc how much memory needs to be allocated
                !call MatMPIAIJSetPreallocation( A, CSR%nnz_per_row_max+1, PETSC_NULL_INTEGER, CSR%nnz_per_row_max+1, PETSC_NULL_INTEGER, perr)
                call MatMPIAIJSetPreallocation( A, lgs%n_terms+1, PETSC_NULL_INTEGER, lgs%n_terms+1, PETSC_NULL_INTEGER, perr)
                call handle_error( 'MatMPIAIJSetPreallocation', perr)
                
                ! Get parallelisation domains ("ownership ranges")
                call MatGetOwnershipRange( A, istart, iend, perr)
                call handle_error( 'MatGetOwnershipRange', perr)
            
                ! Fill in matrix values
                do i = istart+1, iend ! +1 because PETSc indexes from 0
                    do k = lgs%a_ptr(i), lgs%a_ptr(i+1)-1

                        j = lgs%a_index(k)
                        v = lgs%a_value(k)

                        call MatSetValues( A, 1, i-1, 1, j-1, v, INSERT_VALUES, perr)
                        call handle_error( 'MatSetValues', perr)

                    end do
                end do
                call sync
            
                ! == Vectors b,x ==
                ! =================

                ! Create parallel vectors.
                call VecCreate( PETSC_COMM_WORLD, x, perr)
                call handle_error( 'VecCreate', perr)
                !call VecSetSizes( x, PETSC_DECIDE, CSR%n, perr)
                call VecSetSizes( x, PETSC_DECIDE, lgs%nmax, perr)
                call handle_error( 'VecSetSizes', perr)
                call VecSetFromOptions( x, perr)
                call handle_error( 'VecSetFromOptions', perr)
                
                call VecCreate( PETSC_COMM_WORLD, b, perr)
                call handle_error( 'VecCreate', perr)
                !call VecSetSizes( b, PETSC_DECIDE, CSR%m, perr)
                call VecSetSizes( b, PETSC_DECIDE, lgs%nmax, perr)
                call handle_error( 'VecSetSizes', perr)
                call VecSetFromOptions( b, perr)
                call handle_error( 'VecSetFromOptions', perr)
            
                ! Fill in vector values
                do i = istart+1, iend ! +1 because PETSc indexes from 0

                    v = lgs%b_value( i)
                    call VecSetValues( b, 1, i-1, v, INSERT_VALUES, perr)
                    call handle_error( 'VecSetValues', perr)

                    v = lgs%x_value( i)
                    call VecSetValues( x, 1, i-1, v, INSERT_VALUES, perr)
                    call handle_error( 'VecSetValues', perr)

                end do
                call sync
            
                ! Assemble matrix and vectors, using the 2-step process:
                !   MatAssemblyBegin(), MatAssemblyEnd()
                ! Computations can be done while messages are in transition
                ! by placing code between these two statements.
                
                call MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY, perr)
                call handle_error( 'MatAssemblyBegin', perr)
                call VecAssemblyBegin( b, perr)
                call handle_error( 'VecAssemblyBegin', perr)
                call VecAssemblyBegin( x, perr)
                call handle_error( 'VecAssemblyBegin', perr)
                
                call MatAssemblyEnd(   A, MAT_FINAL_ASSEMBLY, perr)
                call handle_error( 'MatAssemblyEnd', perr)
                call VecAssemblyEnd(   b, perr)
                call handle_error( 'VecAssemblyEnd', perr)
                call VecAssemblyEnd(   x, perr)
                call handle_error( 'VecAssemblyEnd', perr)
                
                return

            end subroutine convert_lgs_to_petsc
        
            subroutine convert_petsc_solution_to_lgs(lgs,x)
      
                implicit none
                
                ! In- and output variables:
                type(linear_solver_class), intent(INOUT) :: lgs
                type(tVec), intent(IN) :: x
                
                ! Local variables
                type(PetscErrorCode) :: perr
                type(PetscInt)       :: istart, iend, i
                type(PetscInt)       :: ix(1)
                type(PetscScalar)    :: v(1)
                
                ! Get parallelisation domains ("ownership ranges")
                call VecGetOwnershipRange( x, istart, iend, perr)
                call handle_error( 'VecGetOwnershipRange', perr)
                
                ! Get values
                do i = istart+1,iend
                    ix(1) = i-1
                    call VecGetValues( x, 1, ix, v, perr)
                    call handle_error( 'VecGetValues', perr)
                    lgs%x_value(i) = v(1)
                end do
                call sync
                
                return

            end subroutine convert_petsc_solution_to_lgs
  
            ! == General PETSc initialisation and finalisation
            subroutine initialise_petsc()
                ! Initialise PETSc

                implicit none

                type(PetscErrorCode) :: perr

                ! Initialise PETSc MPI stuff
                call PetscInitialize( PETSC_NULL_CHARACTER, perr)
                call handle_error( 'PetscInitialize', perr)
                
                return

            end subroutine initialise_petsc

            subroutine finalise_petsc()
                ! Finalise PETSc

                implicit none

                type(PetscErrorCode) :: perr

                ! Finalise PETSc MPI stuff
                call PetscFinalize( perr)
                call handle_error( 'PetscFinalize', perr)
                
                return

            end subroutine finalise_petsc
  
            ! == Error handler
            subroutine handle_error( routine_name, perr)

                implicit none 

                ! In/output variables:
                character(LEN=*),     intent(IN) :: routine_Name
                type(PetscErrorCode), intent(IN) :: perr

                ! Local variables
                integer :: cerr, ierr 

                if (perr /= 0) THEN
                  write(0,*) '    PETSc routine "', routine_name, '" returned error flag ', perr
                  call MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
                  stop
                end if

                return 

            end subroutine handle_error

            ! Synchronise the different processes
            subroutine sync
                ! Use MPI_BARRIER to synchronise all the processes         

                implicit none

                integer :: ierr

                call MPI_BARRIER( MPI_COMM_WORLD, ierr)

                return

            end subroutine sync

    end subroutine linear_solver_matrix_solve_petsc

#endif

end module solver_ssa_ac
