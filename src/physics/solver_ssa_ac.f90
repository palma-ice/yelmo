module solver_ssa_ac

    use yelmo_defs, only : sp, dp, wp, io_unit_err, TOL_UNDERFLOW
    use yelmo_tools, only : get_neighbor_indices

    use solver_linear
    use ncio        ! For diagnostic outputting only 

    implicit none 

    private 
    public :: set_ssa_masks
    public :: update_ssa_mask_convergence
    public :: ssa_diagnostics_write_init
    public :: ssa_diagnostics_write_step

    ! Routines that make use of the linear_solver_class object defined in the module solver_linear.F90:
    public :: linear_solver_save_velocity
    public :: linear_solver_matrix_ssa_ac_csr_2D
    
contains 
    
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

            call linear_solver_init(lgs,nx,ny,nvar=2,n_terms=9)

        end if 

        if (count(ssa_mask_acx .eq. 1) + count(ssa_mask_acy .eq. 1) .gt. 0) then 
            ! Points exist for ssa solver to treat

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

        end if 

        ! Define border conditions (only choices are: no-slip, free-slip, periodic)
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

            case("TROUGH")

                bcs(1) = "free-slip"
                bcs(2) = "periodic"
                bcs(3) = "no-slip"
                bcs(4) = "periodic" 

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
            
            case("zeros")

                bcs(1:4) = "no-slip" 
                
            case DEFAULT 

                bcs(1:4) = "no-slip" 
                
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

                if (trim(bcs(3)) .eq. "free-slip") then 
                
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

                else ! bcs(3) == "no-slip"

                    k = k+1
                    lgs%a_value(k)  = 1.0   ! diagonal element only
                    lgs%a_index(k)  = nr

                    lgs%b_value(nr) = 0.0
                    lgs%x_value(nr) = 0.0

                end if 

            else if (i .eq. nx .and. trim(bcs(1)) .ne. "periodic") then 
                ! Right boundary 
                
                if (trim(bcs(1)) .eq. "free-slip") then 
                
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

                else ! bcs(3) == "no-slip"

                    k = k+1
                    lgs%a_value(k)  = 1.0   ! diagonal element only
                    lgs%a_index(k)  = nr

                    lgs%b_value(nr) = 0.0
                    lgs%x_value(nr) = 0.0

                end if 

            else if (j .eq. 1 .and. trim(bcs(4)) .ne. "periodic") then 
                ! Lower boundary 

                if (trim(bcs(4)) .eq. "free-slip") then 

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


                else ! bcs(4) == "no-slip"

                    k = k+1
                    lgs%a_value(k)  = 1.0_wp   ! diagonal element only
                    lgs%a_index(k)  = nr

                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = 0.0_wp

                end if 

            else if (j .eq. ny .and. trim(bcs(2)) .ne. "periodic") then 
                ! Upper boundary 

                if (trim(bcs(2)) .eq. "free-slip") then 
                    
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

                else ! bcs(2) == "no-slip"

                    k = k+1
                    lgs%a_value(k)  = 1.0_wp   ! diagonal element only
                    lgs%a_index(k)  = nr

                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = 0.0_wp

                end if 

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

                if (trim(bcs(4)) .eq. "free-slip") then 

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

                else ! bcs(4) == "no-slip"

                    k = k+1
                    lgs%a_value(k)  = 1.0_wp        ! diagonal element only
                    lgs%a_index(k)  = nr

                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = 0.0_wp

                end if 

            else if (j .eq. ny .and. trim(bcs(2)) .ne. "periodic") then 
                ! Upper boundary 

                if (trim(bcs(2)) .eq. "free-slip") then 

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

                else ! bcs(2) == "no-slip"

                    k = k+1
                    lgs%a_value(k)  = 1.0_wp        ! diagonal element only
                    lgs%a_index(k)  = nr

                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = 0.0_wp

                end if 

            else if (i .eq. 1 .and. trim(bcs(3)) .ne. "periodic") then 
                ! Left boundary 

                if (trim(bcs(3)) .eq. "free-slip") then 

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

                else ! bcs(2) == "no-slip"

                    k = k+1
                    lgs%a_value(k)  = 1.0_wp        ! diagonal element only
                    lgs%a_index(k)  = nr

                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = 0.0_wp

                end if 

            else if (i .eq. nx .and. trim(bcs(1)) .ne. "periodic") then 
                ! Right boundary 

                if (trim(bcs(1)) .eq. "free-slip") then 

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

                else ! bcs(2) == "no-slip"

                    k = k+1
                    lgs%a_value(k)  = 1.0_wp        ! diagonal element only
                    lgs%a_index(k)  = nr

                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = 0.0_wp

                end if 

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

    subroutine check_base_slope(is_steep,zb0,zb1,dx,lim)

        logical,  intent(OUT) :: is_steep
        real(wp), intent(IN) :: zb0         ! [m]
        real(wp), intent(IN) :: zb1         ! [m]
        real(wp), intent(IN) :: dx          ! [m]
        real(wp), intent(IN) :: lim         ! [dx/dx] = [unitless]

        if ( abs(zb1-zb0) / dx .gt. lim ) then 
            is_steep = .TRUE. 
        else 
            is_steep = .FALSE. 
        end if 

        return

    end subroutine check_base_slope


    subroutine set_ssa_masks(ssa_mask_acx,ssa_mask_acy,mask_frnt,H_ice,f_ice, &
                                        f_grnd,z_base,z_sl,dx,use_ssa,lateral_bc)
        ! Define where ssa calculations should be performed
        ! Note: could be binary, but perhaps also distinguish 
        ! grounding line/zone to use this mask for later gl flux corrections
        ! mask = -1: no ssa calculated, velocity imposed/unchanged
        ! mask = 0: no ssa calculated, velocity set to zero
        ! mask = 1: shelfy-stream ssa calculated 
        ! mask = 2: shelf ssa calculated 
        ! mask = 3: ssa lateral boundary condition applied
        ! mask = 4: ssa lateral boundary, but treated as inner ssa

        ! Note: the parameter gradbase_max is used to check slope of ice base. 
        ! If at a given point, it is greater than this limit, the ssa solver
        ! will be disabled in this direction. gradbase_max=0.1 is a relatively
        ! high value, but is reached for points next to deep troughs in Antarctica,
        ! and next to some fjords in Greenland. Steeper slopes are present
        ! in higher-resolution topographies typically.
        
        implicit none 
        
        integer,  intent(OUT) :: ssa_mask_acx(:,:) 
        integer,  intent(OUT) :: ssa_mask_acy(:,:)
        integer,  intent(IN)  :: mask_frnt(:,:)
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: f_grnd(:,:)
        real(wp), intent(IN)  :: z_base(:,:)
        real(wp), intent(IN)  :: z_sl(:,:)
        real(wp), intent(IN)  :: dx 
        logical,  intent(IN)  :: use_ssa       ! SSA is actually active now? 
        character(len=*), intent(IN) :: lateral_bc 

        ! Local variables
        integer  :: i, j, nx, ny
        integer  :: im1, ip1, jm1, jp1
        real(wp) :: H_acx, H_acy
        logical  :: is_steep 
        logical  :: is_convergent 
        
        real(wp), allocatable :: mask_frnt_dyn(:,:)

        ! Integer values for the mask_frnt should be consistent
        ! with those defined in topography.f90:calc_ice_front().
        ! val_disabled is an internal value only used in this routine.
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

                    ! SPECIAL CASE: floating ice next to ice-free land,
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

                    ! SPECIAL CASE: floating ice next to ice-free land,
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

end module solver_ssa_ac
