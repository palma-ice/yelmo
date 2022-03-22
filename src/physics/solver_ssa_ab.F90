module solver_ssa_ab
    ! This ssa solver code was adapted from SICOPOLIS (v5-dev, svn revision 1421) module calc_vxy_m.F90. 
    
    use yelmo_defs, only : sp, dp, wp, io_unit_err, TOL_UNDERFLOW, rho_ice, rho_sw, g 

    use grid_calcs  ! For staggering routines 
    use ncio        ! For diagnostic outputting only 

    implicit none 

    private 
    public :: calc_vxy_ssa_matrix 
    public :: set_ssa_masks
    public :: ssa_diagnostics_write_init
    public :: ssa_diagnostics_write_step


contains 

    subroutine calc_vxy_ssa_matrix(vx_m,vy_m,L2_norm,beta_ab,visc_int_ab, &
                    ssa_mask_ab,H_ice,f_ice,taudx_ab,taudy_ab,H_grnd,z_sl,z_bed, &
                    z_srf,dx,dy,ulim,boundaries,lateral_bc,lis_settings)
        ! Solution of the system of linear equations for the horizontal velocities
        ! vx_m, vy_m in the shallow shelf approximation.
        ! Adapted from sicopolis version 5-dev (svn revision 1421)
        ! Uses the LIS library 
    
        implicit none

        real(wp), intent(INOUT) :: vx_m(:,:)                ! [m yr^-1] Horizontal velocity x (ab-nodes)
        real(wp), intent(INOUT) :: vy_m(:,:)                ! [m yr^-1] Horizontal velocity y (ab-nodes)
        real(wp), intent(OUT)   :: L2_norm                  ! L2 norm convergence check from solver
        real(wp), intent(IN)    :: beta_ab(:,:)             ! [Pa yr m^-1] Basal friction (ab-nodes)
        real(wp), intent(IN)    :: visc_int_ab(:,:)         ! [Pa a m] Vertically integrated viscosity (aa-nodes)
        integer,  intent(IN)    :: ssa_mask_ab(:,:)         ! [--] Mask to determine ssa solver actions (ab-nodes)
        real(wp), intent(IN)    :: H_ice(:,:)               ! [m]  Ice thickness (aa-nodes)
        real(wp), intent(IN)    :: f_ice(:,:)
        real(wp), intent(IN)    :: taudx_ab(:,:)            ! [Pa] Driving stress x (ab nodes)
        real(wp), intent(IN)    :: taudy_ab(:,:)            ! [Pa] Driving stress y (ab nodes)
        real(wp), intent(IN)    :: H_grnd(:,:)  
        real(wp), intent(IN)    :: z_sl(:,:) 
        real(wp), intent(IN)    :: z_bed(:,:) 
        real(wp), intent(IN)    :: z_srf(:,:)
        real(wp), intent(IN)    :: dx, dy
        real(wp), intent(IN)    :: ulim 
        character(len=*), intent(IN) :: boundaries 
        character(len=*), intent(IN) :: lateral_bc
        character(len=*), intent(IN) :: lis_settings

        ! Local variables
        integer  :: nx, ny
        integer  :: i, j, k, n, m 
        integer  :: i1, j1
        real(wp) :: inv_dx, inv_dxdx 
        real(wp) :: inv_dy, inv_dydy 
        real(wp) :: inv_dxdy, inv_2dxdy, inv_4dxdy 
        
        real(wp) :: rho_sw_ice, H_ice_now, beta_now, taud_now, H_ocn_now
        integer  :: IMAX, JMAX 

        integer,  allocatable :: n2i(:), n2j(:)
        integer,  allocatable :: ij2n(:,:)
        integer,  allocatable :: maske(:,:)
        logical,  allocatable :: is_grline_1(:,:) 
        logical,  allocatable :: is_grline_2(:,:) 
        logical,  allocatable :: is_front_1(:,:)
        logical,  allocatable :: is_front_2(:,:)  
        real(wp), allocatable :: visc_int_acx(:,:) 
        real(wp), allocatable :: visc_int_acy(:,:)  
        integer :: n_check 

        ! Grid boundary conditions counterclockwise unit circle 
        ! 1: x, right-border
        ! 2: y, upper-border 
        ! 3: x, left--border 
        ! 4: y, lower-border 
        character(len=56) :: boundaries_vx(4)
        character(len=56) :: boundaries_vy(4)

        integer  :: im1, ip1, jm1, jp1 

        real(wp) :: f_submerged, f_visc
        real(wp) :: vis_int_g_now
        real(wp) :: tau_bc_int 
        real(wp) :: tau_bc_sign

        real(wp), parameter :: f_submerged_min = 0.0_wp 

! Include header for lis solver fortran interface
#include "lisf.h"
        
        LIS_INTEGER :: ierr
        LIS_INTEGER :: nc, nr
        LIS_INTEGER :: lin_iter
        LIS_REAL    :: residual 
        LIS_REAL    :: solver_time 
        LIS_MATRIX  :: lgs_a
        LIS_VECTOR  :: lgs_b, lgs_x
        LIS_SOLVER  :: solver

        LIS_INTEGER :: nmax, n_sprs 
        LIS_INTEGER, allocatable, dimension(:) :: lgs_a_ptr, lgs_a_index
        LIS_SCALAR,  allocatable, dimension(:) :: lgs_a_value, lgs_b_value, lgs_x_value

        ! Define border conditions (zeros, infinite, periodic)
        select case(trim(boundaries)) 

            case("MISMIP3D")

                boundaries_vx(1) = "zeros"
                boundaries_vx(2) = "infinite"
                boundaries_vx(3) = "zeros"
                boundaries_vx(4) = "infinite"

                boundaries_vy(1:4) = "zeros" 

            case("periodic")

                boundaries_vx(1:4) = "periodic" 

                boundaries_vy(1:4) = "periodic" 

            case("periodic-x")

                boundaries_vx(1) = "periodic"
                boundaries_vx(2) = "infinite"
                boundaries_vx(3) = "periodic"
                boundaries_vx(4) = "infinite"

                boundaries_vy(1) = "periodic"
                boundaries_vy(2) = "infinite"
                boundaries_vy(3) = "periodic"
                boundaries_vy(4) = "infinite"

            case("infinite")

                boundaries_vx(1:4) = "infinite" 

                boundaries_vy(1:4) = "infinite" 
                
            case DEFAULT 

                boundaries_vx(1:4) = "zeros" 

                boundaries_vy(1:4) = "zeros" 
                
        end select 

        nx = size(H_ice,1)
        ny = size(H_ice,2)
        
        nmax   =  2*nx*ny 
        n_sprs = 26*nx*ny 

        allocate(n2i(nx*ny),n2j(nx*ny))
        allocate(ij2n(nx,ny))

        allocate(maske(nx,ny))
        allocate(is_grline_1(nx,ny))
        allocate(is_grline_2(nx,ny))
        allocate(is_front_1(nx,ny))
        allocate(is_front_2(nx,ny))
        
        allocate(visc_int_acx(nx,ny))
        allocate(visc_int_acy(nx,ny))

        ! Define some factors

        inv_dx    = 1.0_wp / dx 
        inv_dxdx  = 1.0_wp / (dx*dx)
        inv_dy    = 1.0_wp / dy 
        inv_dydy  = 1.0_wp / (dy*dy)
        inv_dxdy  = 1.0_wp / (dx*dy)
        inv_2dxdy = 1.0_wp / (2.0_wp*dx*dy)
        inv_4dxdy = 1.0_wp / (4.0_wp*dx*dy)
        
        rho_sw_ice   = rho_sw/rho_ice ! Ratio of density of seawater to ice [--]

        ! ===== Consistency checks ==========================

        ! Ensure beta is defined well 
        
        ! TO DO 

        ! Set maske and grounding line / calving front flags

        call set_masks(maske,is_front_1,is_front_2,is_grline_1,is_grline_2, &
                                H_ice, f_ice, H_grnd, z_srf, z_bed, z_sl, lateral_bc)
        
        
        ! Calculate depth-integrated viscosity on acx- and acy-nodes 

        call map_b_to_cx_2D(visc_int_acx,visc_int_ab)
        call map_b_to_cy_2D(visc_int_acy,visc_int_ab)

        ! =======================================================================
        !-------- Reshaping of a 2-d array (with indices i, j)
        !                                  to a vector (with index n) --------
        ! ajr: note, this can be done once outside this routine, but for now
        ! do it here to keep solver portable.

        n=1

        do i=1, nx
        do j=1, ny
            n2i(n)    = i
            n2j(n)    = j
            ij2n(i,j) = n
            n=n+1
        end do
        end do

        ! =======================================================================

        !-------- Assembly of the system of linear equations
        !                         (matrix storage: compressed sparse row CSR) --------

        allocate(lgs_a_value(n_sprs), lgs_a_index(n_sprs), lgs_a_ptr(nmax+1))
        allocate(lgs_b_value(nmax), lgs_x_value(nmax))

        lgs_a_value = 0.0
        lgs_a_index = 0
        lgs_a_ptr   = 0

        lgs_b_value = 0.0
        lgs_x_value = 0.0

        lgs_a_ptr(1) = 1

        k = 0

        do n=1, nmax-1, 2

            i = n2i((n+1)/2)
            j = n2j((n+1)/2)

            !  ------ Equations for vx_m (at (i+1/2,j))

            nr = n   ! row counter

            ! == Treat special cases first ==

            if (ssa_mask_ab(i,j) .eq. -1) then 
                ! Assign prescribed boundary velocity to this point
                ! (eg for prescribed velocity corresponding to 
                ! analytical grounding line flux, or for a regional domain)

                k = k+1
                lgs_a_value(k)  = 1.0   ! diagonal element only
                lgs_a_index(k)  = nr

                lgs_b_value(nr) = vx_m(i,j)
                lgs_x_value(nr) = vx_m(i,j)
           
            else if (i .eq. 1) then 
                ! Left boundary 

                select case(trim(boundaries_vx(3)))

                    case("zeros")
                        ! Assume border velocity is zero 

                        k = k+1
                        lgs_a_value(k)  = 1.0   ! diagonal element only
                        lgs_a_index(k)  = nr

                        lgs_b_value(nr) = 0.0
                        lgs_x_value(nr) = 0.0

                    case("infinite")
                        ! Infinite boundary condition, take 
                        ! value from one point inward

                        nc = 2*ij2n(i,j)-1          ! column counter for vx_m(i,j)
                        k = k+1
                        lgs_a_value(k) =  1.0_wp
                        lgs_a_index(k) = nc

                        nc = 2*ij2n(i+1,j)-1        ! column counter for vx_m(i+1,j)
                        k = k+1
                        lgs_a_value(k) = -1.0_wp
                        lgs_a_index(k) = nc

                        lgs_b_value(nr) = 0.0_wp
                        lgs_x_value(nr) = vx_m(i,j)

                    case("periodic")
                        ! Periodic boundary, take velocity from the right boundary
                        
                        nc = 2*ij2n(i,j)-1          ! column counter for vx_m(i,j)
                        k = k+1
                        lgs_a_value(k) =  1.0_wp
                        lgs_a_index(k) = nc

                        nc = 2*ij2n(nx-2,j)-1        ! column counter for vx_m(nx-2,j)
                        k = k+1
                        lgs_a_value(k) = -1.0_wp
                        lgs_a_index(k) = nc

                        lgs_b_value(nr) = 0.0_wp
                        lgs_x_value(nr) = vx_m(i,j)

                    case DEFAULT 

                        write(*,*) "calc_vxy_ssa_matrix:: Error: left-border condition not &
                        &recognized: "//trim(boundaries_vx(3))
                        write(*,*) "boundaries parameter set to: "//trim(boundaries)
                        stop

                end select 
                
            else if (i .eq. nx) then 
                ! Right boundary 
                
                select case(trim(boundaries_vx(1)))

                    case("zeros")
                        ! Assume border velocity is zero 

                        k = k+1
                        lgs_a_value(k)  = 1.0   ! diagonal element only
                        lgs_a_index(k)  = nr

                        lgs_b_value(nr) = 0.0
                        lgs_x_value(nr) = 0.0

                    case("infinite")
                        ! Infinite boundary condition, take 
                        ! value from two points inward (one point inward will also be prescribed)

                        nc = 2*ij2n(i,j)-1          ! column counter for vx_m(i,j)
                        k = k+1
                        lgs_a_value(k) =  1.0_wp
                        lgs_a_index(k) = nc

                        nc = 2*ij2n(nx-2,j)-1       ! column counter for vx_m(nx-2,j)
                        k = k+1
                        lgs_a_value(k) = -1.0_wp
                        lgs_a_index(k) = nc

                        lgs_b_value(nr) = 0.0_wp
                        lgs_x_value(nr) = vx_m(i,j)

                    case("periodic")
                        ! Periodic boundary condition, take velocity from one point
                        ! interior to the left-border, as nx-1 will be set to value
                        ! at the left-border 

                        nc = 2*ij2n(i,j)-1          ! column counter for vx_m(i,j)
                        k = k+1
                        lgs_a_value(k) =  1.0_wp
                        lgs_a_index(k) = nc

                        nc = 2*ij2n(3,j)-1          ! column counter for vx_m(3,j)
                        k = k+1
                        lgs_a_value(k) = -1.0_wp
                        lgs_a_index(k) = nc

                        lgs_b_value(nr) = 0.0_wp
                        lgs_x_value(nr) = vx_m(i,j)

                    case DEFAULT 

                        write(*,*) "calc_vxy_ssa_matrix:: Error: right-border condition not &
                        &recognized: "//trim(boundaries_vx(1))
                        write(*,*) "boundaries parameter set to: "//trim(boundaries)
                        stop

                end select 

            else if (i .eq. nx-1 .and. trim(boundaries_vx(1)) .eq. "infinite") then 
                ! Right boundary, staggered one point further inward 
                ! (only needed for periodic conditions, otherwise
                ! this point should be treated as normal)
                
                ! Periodic boundary condition, take velocity from one point
                ! interior to the left-border, as nx-1 will be set to value
                ! at the left-border 

                nc = 2*ij2n(i,j)-1          ! column counter for vx_m(i,j)
                k = k+1
                lgs_a_value(k) =  1.0_wp
                lgs_a_index(k) = nc

                nc = 2*ij2n(nx-2,j)-1       ! column counter for vx_m(nx-2,j)
                k = k+1
                lgs_a_value(k) = -1.0_wp
                lgs_a_index(k) = nc

                lgs_b_value(nr) = 0.0_wp
                lgs_x_value(nr) = vx_m(i,j)

            else if (i .eq. nx-1 .and. trim(boundaries_vx(1)) .eq. "periodic") then 
                ! Right boundary, staggered one point further inward 
                ! (only needed for periodic conditions, otherwise
                ! this point should be treated as normal)
                
                ! Periodic boundary condition, take velocity from one point
                ! interior to the left-border, as nx-1 will be set to value
                ! at the left-border 

                nc = 2*ij2n(i,j)-1          ! column counter for vx_m(i,j)
                k = k+1
                lgs_a_value(k) =  1.0_wp
                lgs_a_index(k) = nc

                nc = 2*ij2n(2,j)-1          ! column counter for vx_m(2,j)
                k = k+1
                lgs_a_value(k) = -1.0_wp
                lgs_a_index(k) = nc

                lgs_b_value(nr) = 0.0_wp
                lgs_x_value(nr) = vx_m(i,j)

            else if (j .eq. 1) then 
                ! Lower boundary 

                select case(trim(boundaries_vx(4)))

                    case("zeros")
                        ! Assume border velocity is zero 

                        k = k+1
                        lgs_a_value(k)  = 1.0   ! diagonal element only
                        lgs_a_index(k)  = nr

                        lgs_b_value(nr) = 0.0
                        lgs_x_value(nr) = 0.0

                    case("infinite")
                        ! Infinite boundary condition, take 
                        ! value from one point inward

                        nc = 2*ij2n(i,j)-1          ! column counter for vx_m(i,j)
                        k = k+1
                        lgs_a_value(k) =  1.0_wp
                        lgs_a_index(k) = nc

                        nc = 2*ij2n(i,j+1)-1        ! column counter for vx_m(i,j+1)
                        k = k+1
                        lgs_a_value(k) = -1.0_wp
                        lgs_a_index(k) = nc

                        lgs_b_value(nr) = 0.0_wp
                        lgs_x_value(nr) = vx_m(i,j)

                    case("periodic")
                        ! Periodic boundary, take velocity from the upper boundary
                        
                        nc = 2*ij2n(i,j)-1          ! column counter for vx_m(i,j)
                        k = k+1
                        lgs_a_value(k) =  1.0_wp
                        lgs_a_index(k) = nc

                        nc = 2*ij2n(i,ny-1)-1        ! column counter for vx_m(i,ny-1)
                        k = k+1
                        lgs_a_value(k) = -1.0_wp
                        lgs_a_index(k) = nc

                        lgs_b_value(nr) = 0.0_wp
                        lgs_x_value(nr) = vx_m(i,j)

                    case DEFAULT 

                        write(*,*) "calc_vxy_ssa_matrix:: Error: upper-border condition not &
                        &recognized: "//trim(boundaries_vx(4))
                        write(*,*) "boundaries parameter set to: "//trim(boundaries)
                        stop

                end select 

            else if (j .eq. ny) then 
                ! Upper boundary 

                select case(trim(boundaries_vx(2)))

                    case("zeros")
                        ! Assume border velocity is zero 

                        k = k+1
                        lgs_a_value(k)  = 1.0   ! diagonal element only
                        lgs_a_index(k)  = nr

                        lgs_b_value(nr) = 0.0
                        lgs_x_value(nr) = 0.0

                    case("infinite")
                        ! Infinite boundary condition, take 
                        ! value from one point inward

                        nc = 2*ij2n(i,j)-1          ! column counter for vx_m(i,j)
                        k = k+1
                        lgs_a_value(k) =  1.0_wp
                        lgs_a_index(k) = nc

                        nc = 2*ij2n(i,j-1)-1        ! column counter for vx_m(i,j-1)
                        k = k+1
                        lgs_a_value(k) = -1.0_wp
                        lgs_a_index(k) = nc

                        lgs_b_value(nr) = 0.0_wp
                        lgs_x_value(nr) = vx_m(i,j)

                    case("periodic")
                        ! Periodic boundary, take velocity from the lower boundary
                        
                        nc = 2*ij2n(i,j)-1          ! column counter for vx_m(i,j)
                        k = k+1
                        lgs_a_value(k) =  1.0_wp
                        lgs_a_index(k) = nc

                        nc = 2*ij2n(i,2)-1          ! column counter for vx_m(i,2)
                        k = k+1
                        lgs_a_value(k) = -1.0_wp
                        lgs_a_index(k) = nc

                        lgs_b_value(nr) = 0.0_wp
                        lgs_x_value(nr) = vx_m(i,j)

                    case DEFAULT 

                        write(*,*) "calc_vxy_ssa_matrix:: Error: upper-border condition not &
                        &recognized: "//trim(boundaries_vx(2))
                        write(*,*) "boundaries parameter set to: "//trim(boundaries)
                        stop

                end select 

            else if (ssa_mask_ab(i,j) .eq. 0) then    ! neither neighbour is floating or grounded ice,
                    ! velocity assumed to be zero

                k = k+1
                lgs_a_value(k) = 1.0_wp   ! diagonal element only
                lgs_a_index(k) = nr

                lgs_b_value(nr) = 0.0_wp
                lgs_x_value(nr) = 0.0_wp

            else
                ! === Inner SSA solution === 

                ! -- vx terms -- 

                nc = 2*ij2n(i+1,j)-1        ! column counter for vx_m(i+1,j)
                k = k+1
                lgs_a_value(k) =  4.0_wp*inv_dxdx*visc_int_acy(i+1,j)
                lgs_a_index(k) = nc

                nc = 2*ij2n(i,j)-1          ! column counter for vx_m(i,j)
                k = k+1
                lgs_a_value(k) =  4.0_wp*inv_dxdx*visc_int_acy(i+1,j) &
                                 -4.0_wp*inv_dxdx*visc_int_acy(i,j) &
                                 -1.0_wp*inv_dydy*visc_int_acx(i,j+1) &
                                 -1.0_wp*inv_dydy*visc_int_acx(i,j) &
                                 -beta_ab(i,j)
                lgs_a_index(k) = nc
                
                nc = 2*ij2n(i-1,j)-1        ! column counter for vx_m(i-1,j)
                k = k+1
                lgs_a_value(k) =  4.0_wp*inv_dxdx*visc_int_acy(i,j)
                lgs_a_index(k) = nc

                nc = 2*ij2n(i,j+1)-1        ! column counter for vx_m(i,j+1)
                k = k+1
                lgs_a_value(k) =  1.0_wp*inv_dydy*visc_int_acx(i,j+1)
                lgs_a_index(k) = nc

                nc = 2*ij2n(i,j-1)-1        ! column counter for vx_m(i,j-1)
                k = k+1
                lgs_a_value(k) =  1.0_wp*inv_dydy*visc_int_acx(i,j)
                lgs_a_index(k) = nc

                ! -- vy terms -- 

                nc = 2*ij2n(i,j+1)          ! column counter for vy_m(i,j+1)
                k = k+1
                lgs_a_value(k) = -1.0_wp*inv_2dxdy*visc_int_acy(i+1,j)     &
                                 -1.0_wp*inv_4dxdy*visc_int_acy(i,j)
                lgs_a_index(k) = nc

                nc = 2*ij2n(i,j-1)          ! column counter for vy_m(i,j-1)
                k = k+1
                lgs_a_value(k) = -1.0_wp*inv_2dxdy*visc_int_acy(i+1,j)     &
                                 +1.0_wp*inv_4dxdy*visc_int_acy(i,j)
                lgs_a_index(k) = nc

                nc = 2*ij2n(i+1,j+1)        ! column counter for vy_m(i+1,j+1)
                k = k+1
                lgs_a_value(k) =  2.0_wp*inv_2dxdy*visc_int_acy(i+1,j)     &
                                 +1.0_wp*inv_4dxdy*visc_int_acx(i,j+1)
                lgs_a_index(k) = nc

                nc = 2*ij2n(i+1,j-1)        ! column counter for vy_m(i+1,j-1)
                k = k+1
                lgs_a_value(k) = -1.0_wp*inv_2dxdy*visc_int_acy(i+1,j)     &
                                 -1.0_wp*inv_4dxdy*visc_int_acx(i,j)
                lgs_a_index(k) = nc

                nc = 2*ij2n(i-1,j+1)        ! column counter for vy_m(i-1,j+1)
                k = k+1
                lgs_a_value(k) = -1.0_wp*inv_2dxdy*visc_int_acy(i,j)     &
                                 -1.0_wp*inv_4dxdy*visc_int_acx(i,j+1)
                lgs_a_index(k) = nc

                nc = 2*ij2n(i-1,j-1)        ! column counter for vy_m(i-1,j-1)
                k = k+1
                lgs_a_value(k) =  1.0_wp*inv_2dxdy*visc_int_acy(i,j)     &
                                 +1.0_wp*inv_4dxdy*visc_int_acx(i,j)
                lgs_a_index(k) = nc

                nc = 2*ij2n(i+1,j)          ! column counter for vy_m(i+1,j)
                k = k+1
                lgs_a_value(k) =  1.0_wp*inv_4dxdy*visc_int_acx(i,j+1)     &
                                 -1.0_wp*inv_4dxdy*visc_int_acx(i,j)
                lgs_a_index(k) = nc

                nc = 2*ij2n(i-1,j)          ! column counter for vy_m(i-1,j)
                k = k+1
                lgs_a_value(k) = -1.0_wp*inv_4dxdy*visc_int_acx(i,j+1)     &
                                 +1.0_wp*inv_4dxdy*visc_int_acx(i,j)
                lgs_a_index(k) = nc

                lgs_b_value(nr) = taudx_ab(i,j)
                lgs_x_value(nr) = vx_m(i,j)

            end if

            lgs_a_ptr(nr+1) = k+1   ! row is completed, store index to next row

            !  ------ Equations for vy_m (at (i,j+1/2))

            nr = n+1   ! row counter

            ! == Treat special cases first ==

            if (ssa_mask_ab(i,j) .eq. -1) then 
                ! Assign prescribed boundary velocity to this point
                ! (eg for prescribed velocity corresponding to analytical grounding line flux)

                k = k+1
                lgs_a_value(k)  = 1.0   ! diagonal element only
                lgs_a_index(k)  = nr

                lgs_b_value(nr) = vy_m(i,j)
                lgs_x_value(nr) = vy_m(i,j)
           
            
            else if (j .eq. 1) then 
                ! lower boundary 

                select case(trim(boundaries_vy(4)))

                    case("zeros")
                        ! Assume border vy velocity is zero 

                        k = k+1
                        lgs_a_value(k)  = 1.0_wp   ! diagonal element only
                        lgs_a_index(k)  = nr

                        lgs_b_value(nr) = 0.0_wp
                        lgs_x_value(nr) = 0.0_wp

                    case("infinite")
                        ! Infinite boundary, take velocity from one point inward
                        
                        nc = 2*ij2n(i,j)            ! column counter for vy_m(i,j)
                        k = k+1
                        lgs_a_value(k) =  1.0_wp
                        lgs_a_index(k) = nc

                        nc = 2*ij2n(i,2)            ! column counter for vy_m(i,2)
                        k = k+1
                        lgs_a_value(k) = -1.0_wp
                        lgs_a_index(k) = nc

                        lgs_b_value(nr) = 0.0_wp
                        lgs_x_value(nr) = vy_m(i,j)

                    case("periodic")
                        ! Periodic boundary, take velocity from the opposite boundary
                        
                        nc = 2*ij2n(i,j)            ! column counter for vy_m(i,j)
                        k = k+1
                        lgs_a_value(k) =  1.0_wp
                        lgs_a_index(k) = nc

                        nc = 2*ij2n(i,ny-2)         ! column counter for vy_m(i,ny-2)
                        k = k+1
                        lgs_a_value(k) = -1.0_wp
                        lgs_a_index(k) = nc

                        lgs_b_value(nr) = 0.0_wp
                        lgs_x_value(nr) = vy_m(i,j)

                    case DEFAULT 

                        write(*,*) "calc_vxy_ssa_matrix:: Error: lower-border condition not &
                        &recognized: "//trim(boundaries_vy(4))
                        write(*,*) "boundaries parameter set to: "//trim(boundaries)
                        stop

                end select 

            else if (j .eq. ny) then 
                ! Upper boundary 

                select case(trim(boundaries_vy(2)))

                    case("zeros")
                        ! Assume border velocity is zero 

                        k = k+1
                        lgs_a_value(k)  = 1.0_wp   ! diagonal element only
                        lgs_a_index(k)  = nr

                        lgs_b_value(nr) = 0.0_wp
                        lgs_x_value(nr) = 0.0_wp

                    case("infinite")
                        ! Infinite boundary, take velocity from two points inward
                        ! (to account for staggering)

                        nc = 2*ij2n(i,j)            ! column counter for vy_m(i,j)
                        k = k+1
                        lgs_a_value(k) =  1.0_wp
                        lgs_a_index(k) = nc

                        nc = 2*ij2n(i,ny-2)          ! column counter for vy_m(i,ny-2)
                        k = k+1
                        lgs_a_value(k) = -1.0_wp
                        lgs_a_index(k) = nc

                        lgs_b_value(nr) = 0.0_wp
                        lgs_x_value(nr) = vy_m(i,j)

                    case("periodic")
                        ! Periodic boundary, take velocity from the right boundary
                        
                        nc = 2*ij2n(i,j)            ! column counter for vy_m(i,j)
                        k = k+1
                        lgs_a_value(k) =  1.0_wp
                        lgs_a_index(k) = nc

                        nc = 2*ij2n(i,3)            ! column counter for vy_m(i,3)
                        k = k+1
                        lgs_a_value(k) = -1.0_wp
                        lgs_a_index(k) = nc

                        lgs_b_value(nr) = 0.0_wp
                        lgs_x_value(nr) = vy_m(i,j)

                    case DEFAULT 

                        write(*,*) "calc_vxy_ssa_matrix:: Error: upper-border condition not &
                        &recognized: "//trim(boundaries_vy(2))
                        write(*,*) "boundaries parameter set to: "//trim(boundaries)
                        stop

                end select 
            
            else if (j .eq. ny-1 .and. trim(boundaries_vy(2)) .eq. "infinite") then
                ! Upper boundary, inward by one point
                ! (only needed for periodic conditions, otherwise
                ! this point should be treated as normal)
                
                ! Periodic boundary, take velocity from the lower boundary

                nc = 2*ij2n(i,j)            ! column counter for vy_m(i,j)
                k = k+1
                lgs_a_value(k) =  1.0_wp
                lgs_a_index(k) = nc

                nc = 2*ij2n(i,ny-2)         ! column counter for vy_m(i,ny-2)
                k = k+1
                lgs_a_value(k) = -1.0_wp
                lgs_a_index(k) = nc

                lgs_b_value(nr) = 0.0_wp
                lgs_x_value(nr) = vy_m(i,j)

            else if (j .eq. ny-1 .and. trim(boundaries_vy(2)) .eq. "periodic") then
                ! Upper boundary, inward by one point
                ! (only needed for periodic conditions, otherwise
                ! this point should be treated as normal)
                
                ! Periodic boundary, take velocity from the lower boundary
                
                nc = 2*ij2n(i,j)            ! column counter for vy_m(i,j)
                k = k+1
                lgs_a_value(k) =  1.0_wp
                lgs_a_index(k) = nc

                nc = 2*ij2n(i,2)            ! column counter for vy_m(i,2)
                k = k+1
                lgs_a_value(k) = -1.0_wp
                lgs_a_index(k) = nc

                lgs_b_value(nr) = 0.0_wp
                lgs_x_value(nr) = vy_m(i,j)
  
            else if (i .eq. 1) then 
                ! Left boundary 

                select case(trim(boundaries_vy(3)))

                    case("zeros")
                        ! Assume border velocity is zero 

                        k = k+1
                        lgs_a_value(k)  = 1.0_wp   ! diagonal element only
                        lgs_a_index(k)  = nr

                        lgs_b_value(nr) = 0.0_wp
                        lgs_x_value(nr) = 0.0_wp

                    case("infinite")
                        ! Infinite boundary, take velocity from one point inward

                        nc = 2*ij2n(i,j)            ! column counter for vy_m(i,j)
                        k = k+1
                        lgs_a_value(k) =  1.0_wp
                        lgs_a_index(k) = nc

                        nc = 2*ij2n(2,j)            ! column counter for vy_m(2,j)
                        k = k+1
                        lgs_a_value(k) = -1.0_wp
                        lgs_a_index(k) = nc

                        lgs_b_value(nr) = 0.0_wp
                        lgs_x_value(nr) = vy_m(i,j)

                    case("periodic")
                        ! Periodic boundary, take velocity from the right boundary
                        
                        nc = 2*ij2n(i,j)            ! column counter for vy_m(i,j)
                        k = k+1
                        lgs_a_value(k) =  1.0_wp
                        lgs_a_index(k) = nc

                        nc = 2*ij2n(nx-1,j)         ! column counter for vy_m(nx-1,j)
                        k = k+1
                        lgs_a_value(k) = -1.0_wp
                        lgs_a_index(k) = nc

                        lgs_b_value(nr) = 0.0_wp
                        lgs_x_value(nr) = vy_m(i,j)

                    case DEFAULT 

                        write(*,*) "calc_vxy_ssa_matrix:: Error: left-border condition not &
                        &recognized: "//trim(boundaries_vy(3))
                        write(*,*) "boundaries parameter set to: "//trim(boundaries)
                        stop

                end select 
                
            else if (i .eq. nx) then 
                ! Right boundary 

                select case(trim(boundaries_vy(1)))

                    case("zeros")
                        ! Assume border velocity is zero 

                        k = k+1
                        lgs_a_value(k)  = 1.0_wp   ! diagonal element only
                        lgs_a_index(k)  = nr

                        lgs_b_value(nr) = 0.0_wp
                        lgs_x_value(nr) = 0.0_wp

                    case("infinite")
                        ! Infinite boundary, take velocity from one point inward

                        nc = 2*ij2n(i,j)            ! column counter for vy_m(i,j)
                        k = k+1
                        lgs_a_value(k) =  1.0_wp
                        lgs_a_index(k) = nc

                        nc = 2*ij2n(nx-1,j)         ! column counter for vy_m(nx-1,j)
                        k = k+1
                        lgs_a_value(k) = -1.0_wp
                        lgs_a_index(k) = nc

                        lgs_b_value(nr) = 0.0_wp
                        lgs_x_value(nr) = vy_m(i,j)

                    case("periodic")
                        ! Periodic boundary, take velocity from the right boundary
                        
                        nc = 2*ij2n(i,j)            ! column counter for vy_m(i,j)
                        k = k+1
                        lgs_a_value(k) =  1.0_wp
                        lgs_a_index(k) = nc

                        nc = 2*ij2n(2,j)            ! column counter for vy_m(2,j)
                        k = k+1
                        lgs_a_value(k) = -1.0_wp
                        lgs_a_index(k) = nc

                        lgs_b_value(nr) = 0.0_wp
                        lgs_x_value(nr) = vy_m(i,j)

                    case DEFAULT 

                        write(*,*) "calc_vxy_ssa_matrix:: Error: right-border condition not &
                        &recognized: "//trim(boundaries_vy(1))
                        write(*,*) "boundaries parameter set to: "//trim(boundaries)
                        stop

                end select 

            else if (ssa_mask_ab(i,j) .eq. 0) then    ! neither neighbour is floating or grounded ice,
                    ! velocity assumed to be zero

                k = k+1
                lgs_a_value(k)  = 1.0_wp   ! diagonal element only
                lgs_a_index(k)  = nr

                lgs_b_value(nr) = 0.0_wp
                lgs_x_value(nr) = 0.0_wp

            else
                ! === Inner SSA solution ===

                ! -- vy terms -- 

                nc = 2*ij2n(i,j+1)        ! column counter for vy_m(i,j+1)
                k = k+1
                lgs_a_value(k) =  4.0_wp*inv_dxdx*visc_int_acx(i,j+1)
                lgs_a_index(k) = nc

                nc = 2*ij2n(i,j)          ! column counter for vy_m(i,j)
                k = k+1
                lgs_a_value(k) =  4.0_wp*inv_dydy*visc_int_acx(i,j+1) &
                                 -4.0_wp*inv_dydy*visc_int_acx(i,j) &
                                 -1.0_wp*inv_dxdx*visc_int_acy(i+1,j) &
                                 -1.0_wp*inv_dxdx*visc_int_acy(i,j) &
                                 -beta_ab(i,j)
                lgs_a_index(k) = nc
                
                nc = 2*ij2n(i,j-1)        ! column counter for vy_m(i,j-1)
                k = k+1
                lgs_a_value(k) =  4.0_wp*inv_dydy*visc_int_acx(i,j)
                lgs_a_index(k) = nc

                nc = 2*ij2n(i+1,j)        ! column counter for vy_m(i+1,j)
                k = k+1
                lgs_a_value(k) =  1.0_wp*inv_dxdx*visc_int_acy(i+1,j)
                lgs_a_index(k) = nc

                nc = 2*ij2n(i-1,j)        ! column counter for vy_m(i-1,j)
                k = k+1
                lgs_a_value(k) =  1.0_wp*inv_dxdx*visc_int_acy(i,j)
                lgs_a_index(k) = nc

                ! -- vx terms -- 

                nc = 2*ij2n(i+1,j)-1          ! column counter for vx_m(i+1,j)
                k = k+1
                lgs_a_value(k) = -1.0_wp*inv_2dxdy*visc_int_acx(i,j+1)     &
                                 -1.0_wp*inv_4dxdy*visc_int_acx(i,j)
                lgs_a_index(k) = nc

                nc = 2*ij2n(i-1,j)-1          ! column counter for vx_m(i-1,j)
                k = k+1
                lgs_a_value(k) = -1.0_wp*inv_2dxdy*visc_int_acx(i,j+1)     &
                                 +1.0_wp*inv_4dxdy*visc_int_acx(i,j)
                lgs_a_index(k) = nc

                nc = 2*ij2n(i+1,j+1)-1        ! column counter for vx_m(i+1,j+1)
                k = k+1
                lgs_a_value(k) =  2.0_wp*inv_2dxdy*visc_int_acx(i,j+1)     &
                                 +1.0_wp*inv_4dxdy*visc_int_acy(i+1,j)
                lgs_a_index(k) = nc

                nc = 2*ij2n(i-1,j+1)-1        ! column counter for vx_m(i-1,j+1)
                k = k+1
                lgs_a_value(k) = -1.0_wp*inv_2dxdy*visc_int_acx(i,j+1)     &
                                 -1.0_wp*inv_4dxdy*visc_int_acy(i,j)
                lgs_a_index(k) = nc

                nc = 2*ij2n(i+1,j-1)-1        ! column counter for vx_m(i+1,j-1)
                k = k+1
                lgs_a_value(k) = -1.0_wp*inv_2dxdy*visc_int_acx(i,j)     &
                                 -1.0_wp*inv_4dxdy*visc_int_acy(i+1,j)
                lgs_a_index(k) = nc

                nc = 2*ij2n(i-1,j-1)-1        ! column counter for vx_m(i-1,j-1)
                k = k+1
                lgs_a_value(k) =  1.0_wp*inv_2dxdy*visc_int_acx(i,j)     &
                                 +1.0_wp*inv_4dxdy*visc_int_acy(i,j)
                lgs_a_index(k) = nc

                nc = 2*ij2n(i,j+1)-1          ! column counter for vx_m(i,j+1)
                k = k+1
                lgs_a_value(k) =  1.0_wp*inv_4dxdy*visc_int_acy(i+1,j)     &
                                 -1.0_wp*inv_4dxdy*visc_int_acy(i,j)
                lgs_a_index(k) = nc

                nc = 2*ij2n(i,j-1)-1          ! column counter for vx_m(i,j-1)
                k = k+1
                lgs_a_value(k) = -1.0_wp*inv_4dxdy*visc_int_acy(i+1,j)     &
                                 +1.0_wp*inv_4dxdy*visc_int_acy(i,j)
                lgs_a_index(k) = nc

                lgs_b_value(nr) = taudy_ab(i,j)
                lgs_x_value(nr) = vy_m(i,j)

            end if

            lgs_a_ptr(nr+1) = k+1   ! row is completed, store index to next row

        end do

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

            do nc=lgs_a_ptr(nr), lgs_a_ptr(nr+1)-1
                call lis_matrix_set_value(LIS_INS_VALUE, nr, lgs_a_index(nc), &
                                                        lgs_a_value(nc), lgs_a, ierr)
            end do

            call lis_vector_set_value(LIS_INS_VALUE, nr, lgs_b_value(nr), lgs_b, ierr)
            call lis_vector_set_value(LIS_INS_VALUE, nr, lgs_x_value(nr), lgs_x, ierr)

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
        
        ! Obtain the relative L2_norm == ||b-Ax||_2 / ||b||_2
        call lis_solver_get_residualnorm(solver,residual,ierr)
        L2_norm = real(residual,wp) 

        ! Print a summary
        write(*,*) "calc_vxy_ssa_matrix_ab: [time (s), iter, L2] = ", solver_time, lin_iter, residual


        lgs_x_value = 0.0_wp
        call lis_vector_gather(lgs_x, lgs_x_value, ierr)
        call CHKERR(ierr)

        call lis_matrix_destroy(lgs_a, ierr)
        call CHKERR(ierr)

        call lis_vector_destroy(lgs_b, ierr)
        call lis_vector_destroy(lgs_x, ierr)
        call lis_solver_destroy(solver, ierr)
        call CHKERR(ierr)

        call lis_finalize(ierr)           ! Important for parallel computing environments
        call CHKERR(ierr)

        do n=1, nmax-1, 2

            i = n2i((n+1)/2)
            j = n2j((n+1)/2)

            nr = n
            vx_m(i,j) = lgs_x_value(nr)

            nr = n+1
            vy_m(i,j) = lgs_x_value(nr)

        end do

        deallocate(lgs_a_value, lgs_a_index, lgs_a_ptr)
        deallocate(lgs_b_value, lgs_x_value)


        ! Limit the velocity generally =====================
        call limit_vel(vx_m,ulim)
        call limit_vel(vy_m,ulim)

        return 

    end subroutine calc_vxy_ssa_matrix

    subroutine set_ssa_masks(ssa_mask_ab,H_ice,f_ice,f_grnd_ab,use_ssa)
        ! Define where ssa calculations should be performed
        ! Note: could be binary, but perhaps also distinguish 
        ! grounding line/zone to use this mask for later gl flux corrections
        ! mask = 0: no ssa calculated
        ! mask = 1: shelfy-stream ssa calculated 
        ! mask = 2: shelf ssa calculated 

        implicit none 
        
        integer,  intent(OUT) :: ssa_mask_ab(:,:)
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: f_grnd_ab(:,:)
        logical,  intent(IN)  :: use_ssa        ! SSA is actually active now? 

        ! Local variables
        integer    :: i, j, nx, ny
        integer    :: im1, ip1, jm1, jp1
        real(wp) :: H_ab
        
        nx = size(H_ice,1)
        ny = size(H_ice,2)
        
        ! Initially no active ssa points
        ssa_mask_ab = 0
        
        if (use_ssa) then 

            do j = 1, ny
            do i = 1, nx

                ! Get neighbor indices
                im1 = max(i-1,1) 
                ip1 = min(i+1,nx) 
                jm1 = max(j-1,1) 
                jp1 = min(j+1,ny)


                ! x-direction
                if (f_ice(i,j) .eq. 1.0 .or. f_ice(ip1,j) .eq. 1.0 .or. &
                    f_ice(i,jp1) .eq. 1.0 .or. f_ice(ip1,jp1) .eq. 1.0) then
                
                    ! Ice is present on ab-node
                    
                    if (f_grnd_ab(i,j) .gt. 0.0) then 
                        ! Grounded ice or grounding line (ie, shelfy-stream)
                        ssa_mask_ab(i,j) = 1
                    else 
                        ! Shelf ice 
                        ssa_mask_ab(i,j) = 2
                    end if 

                end if

            end do 
            end do

        end if 

        return
        
    end subroutine set_ssa_masks

! === INTERNAL ROUTINES ==== 
    
    subroutine set_masks(maske,front1,front2,gl1,gl2,H_ice,f_ice,H_grnd,z_srf,z_bed,z_sl,apply_lateral_bc)
        ! Define where ssa calculations should be performed
        ! Note: could be binary, but perhaps also distinguish 
        ! grounding line/zone to use this mask for later gl flux corrections
        ! mask = 0: Grounded ice 
        ! mask = 1: Ice-free land 
        ! mask = 2: Open ocean  
        ! mask = 3: Ice shelf 

        ! Note: this mask is defined on central aa-nodes 
        
        implicit none 
        
        integer,    intent(OUT) :: maske(:,:) 
        logical,    intent(OUT) :: front1(:,:)
        logical,    intent(OUT) :: front2(:,:)  
        logical,    intent(OUT) :: gl1(:,:)
        logical,    intent(OUT) :: gl2(:,:) 
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: H_grnd(:,:)
        real(wp), intent(IN)  :: z_srf(:,:)
        real(wp), intent(IN)  :: z_bed(:,:)
        real(wp), intent(IN)  :: z_sl(:,:)
        character(len=*), intent(IN) :: apply_lateral_bc 

        ! Local variables
        integer  :: i, j, nx, ny
        integer  :: im1, ip1, jm1, jp1 
        logical  :: is_float 
        real(wp) :: f_submerged
        real(wp) :: H_ocn_now 

        ! Use 'apply_lateral_bc' to determine where to apply generalized
        ! lateral bc equation.
        ! "floating" : only apply at floating ice margins 
        ! "marine"   : only apply at floating ice margins and 
        !              grounded marine ice margins (grounded ice next to open ocean)
        ! "all"      : apply at all ice margins 
        ! "none"     : do not apply boundary condition (for testing mainly)

        nx = size(maske,1)
        ny = size(maske,2)
        
        gl1 = .FALSE. 
        gl2 = .FALSE. 

        ! First determine general ice coverage mask 
        ! (land==1,ocean==2,floating_ice==3,grounded_ice==0)

        do j = 1, ny
        do i = 1, nx
            
            ! Check if this point would be floating
            is_float = H_grnd(i,j) .le. 0.0 

            if (f_ice(i,j) .eq. 1.0) then
                ! Ice-covered point

                if (is_float) then 
                    ! Ice shelf 
                    maske(i,j) = 3
                else
                    ! Grounded ice 
                    maske(i,j) = 0 
                end if 
                
            else 
                ! Ice-free point, or only partially-covered (consider ice free)

                if (is_float) then 
                    ! Open ocean 
                    maske(i,j) = 2
                else 
                    ! Ice-free land 
                    maske(i,j) = 1 
                end if 

            end if 

        end do 
        end do 
        
        !-------- Detection of the ice front/margins --------

        front1  = .false. 
        front2  = .false. 

        do j = 1, ny
        do i = 1, nx

            ! Get neighbor indices 
            im1 = max(i-1,1)
            ip1 = min(i+1,nx)
            jm1 = max(j-1,1)
            jp1 = min(j+1,ny)
            
          if (  f_ice(i,j) .eq. 1.0 .and. &
                ( f_ice(im1,j) .lt. 1.0 .or. &
                  f_ice(ip1,j) .lt. 1.0 .or. &
                  f_ice(i,jm1) .lt. 1.0 .or. &
                  f_ice(i,jp1) .lt. 1.0 ) ) then
            front1(i,j) = .TRUE. 

          end if 
          
          if (  f_ice(i,j) .lt. 1.0 .and. &
                ( f_ice(im1,j) .eq. 1.0 .or. &
                  f_ice(ip1,j) .eq. 1.0 .or. &
                  f_ice(i,jm1) .eq. 1.0 .or. &
                  f_ice(i,jp1) .eq. 1.0 ) ) then

            front2(i,j) = .TRUE. 

          end if 

          ! So far all margins have been diagnosed (marine and grounded on land)
          ! Disable some regions depending on choice above. 
          select case(trim(apply_lateral_bc))
            ! Apply the lateral boundary condition to what? 

            case("none")
                ! Do not apply lateral bc anywhere. Ie, disable front detection.
                ! This is mainly for testing, as the 'inner' ssa section
                ! is used and may take information from neighbors outside
                ! the ice sheet margin.

                if (front1(i,j)) front1(i,j) = .FALSE. 

            case("floating","float","slab","slab-ext")
                ! Only apply lateral bc to floating ice fronts.
                ! Ie, disable detection of all grounded fronts for now.
                ! Model is generally more stable this way.
                ! This method is also used for the 'infinite slab' approach,
                ! where a thin ice shelf is extended everywhere over the domain. 
                
                if ( front1(i,j) .and. maske(i,j) .eq. 0 ) front1(i,j) = .FALSE. 

            case("marine")
                ! Only apply lateral bc to floating ice fronts and
                ! and grounded marine fronts. Disable detection 
                ! of ice fronts grounded above sea level.
                
                ! Get current ocean thickness bordering ice sheet
                ! (for bedrock above sea level, this will give zero)
                f_submerged = 1.d0 - min((z_srf(i,j)-z_sl(i,j))/H_ice(i,j),1.d0)
                H_ocn_now   = H_ice(i,j)*f_submerged

                ! if ( front1(i,j) .and. maske(i,j) .eq. 0 .and. &
                !                     H_ocn_now .eq. 0.0 ) front1(i,j) = .FALSE. 
                
                if ( front1(i,j) .and. maske(i,j) .eq. 0 .and. &
                                    f_submerged .eq. 0.0_wp ) front1(i,j) = .FALSE. 
            
            case("all")
                ! Apply lateral bc to all ice-sheet fronts. 

                ! Do nothing - all fronts have been accurately diagnosed. 

            case DEFAULT
                
                write(io_unit_err,*) "set_masks:: error: ssa_lat_bc parameter value not recognized."
                write(io_unit_err,*) "ydyn.ssa_lat_bc = ", apply_lateral_bc
                stop 

           end select

        end do
        end do
        
        return
        
    end subroutine set_masks
    
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

    subroutine ssa_diagnostics_write_step(filename,ux,uy,L2_norm,beta_ab,visc_int_ab, &
                                        ssa_mask_ab,ssa_err_x,ssa_err_y,H_ice,f_ice,taudx_ab,taudy_ab, &
                                                            H_grnd,z_sl,z_bed,z_srf,ux_prev,uy_prev,time)

        implicit none 
        
        character(len=*),  intent(IN) :: filename
        real(wp), intent(IN) :: ux(:,:) 
        real(wp), intent(IN) :: uy(:,:) 
        real(wp), intent(IN) :: L2_norm
        real(wp), intent(IN) :: beta_ab(:,:)
        real(wp), intent(IN) :: visc_int_ab(:,:) 
        integer,  intent(IN) :: ssa_mask_ab(:,:)
        real(wp), intent(IN) :: ssa_err_x(:,:) 
        real(wp), intent(IN) :: ssa_err_y(:,:) 
        real(wp), intent(IN) :: H_ice(:,:) 
        real(wp), intent(IN) :: f_ice(:,:) 
        real(wp), intent(IN) :: taudx_ab(:,:) 
        real(wp), intent(IN) :: taudy_ab(:,:) 
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

        call nc_write(filename,"beta_ab",beta_ab,units="Pa yr m^-1",long_name="Dragging coefficient (ab)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"visc_int_ab",visc_int_ab,units="Pa yr",long_name="Vertically integrated effective viscosity (ab)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"ssa_mask_ab",ssa_mask_ab,units="1",long_name="SSA mask (ab)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"ssa_err_x",ssa_err_x,units="1",long_name="SSA err (x, ab)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"ssa_err_y",ssa_err_y,units="1",long_name="SSA err (y, ab)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"H_ice",H_ice,units="m",long_name="Ice thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_ice",f_ice,units="1",long_name="Ice-covered fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"taudx_ab",taudx_ab,units="Pa",long_name="Driving stress (x, ab)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taudy_ab",taudy_ab,units="Pa",long_name="Driving stress (y, ab)", &
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
    
end module solver_ssa_ab
