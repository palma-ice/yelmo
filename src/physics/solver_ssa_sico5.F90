module solver_ssa_sico5
    ! This ssa solver code was adapted from SICOPOLIS (v5-dev, svn revision 1421) module calc_vxy_m.F90. 
    
    use yelmo_defs, only : sp, dp, wp, prec, rho_ice, rho_sw, g 

    implicit none 

    private 
    public :: calc_vxy_ssa_matrix 
    public :: set_ssa_masks
    public :: update_ssa_mask_convergence
    public :: check_vel_convergence_l1rel_matrix
    public :: check_vel_convergence_l2rel
    public :: relax_ssa
    
contains 

    subroutine calc_vxy_ssa_matrix(vx_m,vy_m,L2_norm,beta_acx,beta_acy,visc_eff, &
                    ssa_mask_acx,ssa_mask_acy,H_ice,f_ice,taud_acx,taud_acy,H_grnd,z_sl,z_bed, &
                    dx,dy,ulim,boundaries,lis_settings)
        ! Solution of the system of linear equations for the horizontal velocities
        ! vx_m, vy_m in the shallow shelf approximation.
        ! Adapted from sicopolis version 5-dev (svn revision 1421)
        ! Uses the LIS library 
    
        implicit none

        real(prec), intent(INOUT) :: vx_m(:,:)            ! [m a^-1] Horizontal velocity x (acx-nodes)
        real(prec), intent(INOUT) :: vy_m(:,:)            ! [m a^-1] Horizontal velocity y (acy-nodes)
        real(prec), intent(OUT)   :: L2_norm              ! L2 norm convergence check from solver
        real(prec), intent(IN)    :: beta_acx(:,:)        ! [Pa a m^-1] Basal friction (acx-nodes)
        real(prec), intent(IN)    :: beta_acy(:,:)        ! [Pa a m^-1] Basal friction (acy-nodes)
        real(prec), intent(IN)    :: visc_eff(:,:)        ! [Pa a m] Vertically integrated viscosity (aa-nodes)
        integer,    intent(IN)    :: ssa_mask_acx(:,:)    ! [--] Mask to determine ssa solver actions (acx-nodes)
        integer,    intent(IN)    :: ssa_mask_acy(:,:)    ! [--] Mask to determine ssa solver actions (acy-nodes)
        real(prec), intent(IN)    :: H_ice(:,:)           ! [m]  Ice thickness (aa-nodes)
        real(prec), intent(IN)    :: f_ice(:,:)
        real(prec), intent(IN)    :: taud_acx(:,:)        ! [Pa] Driving stress (acx nodes)
        real(prec), intent(IN)    :: taud_acy(:,:)        ! [Pa] Driving stress (acy nodes)
        real(prec), intent(IN)    :: H_grnd(:,:)  
        real(prec), intent(IN)    :: z_sl(:,:) 
        real(prec), intent(IN)    :: z_bed(:,:) 
        real(prec), intent(IN)    :: dx, dy
        real(prec), intent(IN)    :: ulim 
        character(len=*), intent(IN) :: boundaries 
        character(len=*), intent(IN) :: lis_settings

        ! Local variables
        integer    :: nx, ny
        real(prec) :: dxi, deta
        integer    :: i, j, k, n, m 
        integer    :: i1, j1, i00, j00
        real(prec) :: inv_dxi, inv_deta, inv_dxi_deta, inv_dxi2, inv_deta2
        real(prec) :: factor_rhs_2, factor_rhs_3a, factor_rhs_3b
        real(prec) :: rho_sw_ice, H_ice_now, beta_now, taud_now, H_ocn_now    
        integer    :: IMAX, JMAX 

        integer, allocatable    :: n2i(:), n2j(:)
        integer, allocatable    :: ij2n(:,:)
        integer, allocatable    :: maske(:,:)
        logical, allocatable    :: is_grline_1(:,:) 
        logical, allocatable    :: is_grline_2(:,:) 
        logical, allocatable    :: is_front_1(:,:)
        logical, allocatable    :: is_front_2(:,:)  
        real(prec), allocatable :: vis_int_g(:,:) 
        real(prec), allocatable :: vis_int_sgxy(:,:)  
        integer :: n_check 

        ! Boundary conditions counterclockwise unit circle 
        ! 1: x, right-border
        ! 2: y, upper-border 
        ! 3: x, left-border 
        ! 4: y, lower-border 
        character(len=56) :: boundaries_vx(4)
        character(len=56) :: boundaries_vy(4)

! Include header for lis solver fortran interface
#include "lisf.h"
        
        LIS_INTEGER :: ierr
        LIS_INTEGER :: nc, nr
        LIS_REAL    :: residual 
        LIS_MATRIX  :: lgs_a
        LIS_VECTOR  :: lgs_b, lgs_x
        LIS_SOLVER  :: solver

        LIS_INTEGER :: nmax, n_sprs 
        LIS_INTEGER, allocatable, dimension(:) :: lgs_a_ptr, lgs_a_index
        LIS_SCALAR,  allocatable, dimension(:) :: lgs_a_value, lgs_b_value, lgs_x_value

        ! Define border conditions (zeros, infinite, periodic)
        select case(trim(boundaries)) 

            case("MISMIP3D")

                boundaries_vx(1) = "infinite"
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
        n_sprs = 20*nx*ny 

        allocate(n2i(nx*ny),n2j(nx*ny))
        allocate(ij2n(nx,ny))

        allocate(maske(nx,ny))
        allocate(is_grline_1(nx,ny))
        allocate(is_grline_2(nx,ny))
        allocate(is_front_1(nx,ny))
        allocate(is_front_2(nx,ny))
        
        allocate(vis_int_g(nx,ny))
        allocate(vis_int_sgxy(nx,ny))

        !--- External yelmo arguments => local sicopolis variable names ---
        dxi          = dx 
        deta         = dy 

        vis_int_g    = visc_eff 

        rho_sw_ice   = rho_sw/rho_ice ! Ratio of density of seawater to ice [--]

        ! ===== Consistency checks ==========================

        ! Ensure beta is defined well 
        if ( count(H_grnd .gt. 100.0) .gt. 0 .and. count(beta_acx .gt. 0.0 .and. H_grnd .gt. 100.0) .eq. 0 ) then  
            ! No points found with a non-zero beta for grounded ice,
            ! something was not well-defined/well-initialized

            write(*,*) 
            write(*,"(a)") "calc_vxy_ssa_matrix:: Error: beta appears to be zero everywhere for grounded ice."
            write(*,*) "range(beta_acx): ", minval(beta_acx), maxval(beta_acx)
            write(*,*) "range(beta_acy): ", minval(beta_acy), maxval(beta_acy)
            write(*,*) "range(H_grnd):   ", minval(H_grnd), maxval(H_grnd)
            write(*,*) "Stopping."
            write(*,*) 
            stop 
            
        end if 

        !-------- Abbreviations --------

        inv_dxi       = 1.0_prec/dxi
        inv_deta      = 1.0_prec/deta
        inv_dxi_deta  = 1.0_prec/(dxi*deta)
        inv_dxi2      = 1.0_prec/(dxi*dxi)
        inv_deta2     = 1.0_prec/(deta*deta)

        factor_rhs_2  = 0.5_prec*rho_ice*g*(rho_sw-rho_ice)/rho_sw
        factor_rhs_3a = 0.5_prec*rho_ice*g
        factor_rhs_3b = 0.5_prec*rho_sw*g

        ! Set maske and grounding line / calving front flags

        call set_sico_masks(maske,is_front_1,is_front_2,is_grline_1,is_grline_2, H_ice, f_ice, H_grnd)
        
        !-------- Depth-integrated viscosity on the staggered grid
        !                                       [at (i+1/2,j+1/2)] --------

        call stagger_visc_aa_ab(vis_int_sgxy,vis_int_g,H_ice,f_ice)

        !-------- Basal drag parameter (for shelfy stream) --------

        !  (now provided as input matrix)

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

            ! Set current boundary variables for later access
            beta_now = beta_acx(i,j) 
            taud_now = taud_acx(i,j) 

            ! == Treat special cases first ==

            if (ssa_mask_acx(i,j) .eq. -1) then 
                ! Assign prescribed boundary velocity to this point
                ! (eg for prescribed velocity corresponding to analytical grounding line flux)

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

                        k = k+1
                        lgs_a_value(k)  = 1.0   ! diagonal element only
                        lgs_a_index(k)  = nr

                        lgs_b_value(nr) = vx_m(2,j)
                        lgs_x_value(nr) = vx_m(2,j)

                    case("periodic")
                        ! Periodic boundary, take velocity from the right boundary
                        
                        k = k+1
                        lgs_a_value(k)  = 1.0   ! diagonal element only
                        lgs_a_index(k)  = nr

                        lgs_b_value(nr) = vx_m(nx-2,j)
                        lgs_x_value(nr) = vx_m(nx-2,j)
                
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

                        k = k+1
                        lgs_a_value(k)  = 1.0   ! diagonal element only
                        lgs_a_index(k)  = nr

                        lgs_b_value(nr) = vx_m(nx-2,j)
                        lgs_x_value(nr) = vx_m(nx-2,j)
                    
                    case("periodic")
                        ! Periodic boundary condition, take velocity from one point
                        ! interior to the left-border, as nx-1 will be set to value
                        ! at the left-border 

                        k = k+1
                        lgs_a_value(k)  = 1.0   ! diagonal element only
                        lgs_a_index(k)  = nr

                        lgs_b_value(nr) = vx_m(3,j)
                        lgs_x_value(nr) = vx_m(3,j)
                    
                    case DEFAULT 

                        write(*,*) "calc_vxy_ssa_matrix:: Error: right-border condition not &
                        &recognized: "//trim(boundaries_vx(1))
                        write(*,*) "boundaries parameter set to: "//trim(boundaries)
                        stop

                end select 

            else if (i .eq. nx-1 .and. trim(boundaries_vx(1)) .eq. "periodic") then 
                ! Right boundary, staggered one point further inward 
                ! (only needed for periodic conditions, otherwise
                ! this point should be treated as normal)
                
                ! Periodic boundary condition, take velocity from one point
                ! interior to the left-border, as nx-1 will be set to value
                ! at the left-border 

                k = k+1
                lgs_a_value(k)  = 1.0   ! diagonal element only
                lgs_a_index(k)  = nr

                lgs_b_value(nr) = vx_m(2,j)
                lgs_x_value(nr) = vx_m(2,j)

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

                        k = k+1
                        lgs_a_value(k)  = 1.0   ! diagonal element only
                        lgs_a_index(k)  = nr

                        lgs_b_value(nr) = vx_m(i,j+1)
                        lgs_x_value(nr) = vx_m(i,j+1)
                
                    case("periodic")
                        ! Periodic boundary, take velocity from the upper boundary
                    
                        k = k+1
                        lgs_a_value(k)  = 1.0   ! diagonal element only
                        lgs_a_index(k)  = nr

                        lgs_b_value(nr) = vx_m(i,ny-1)
                        lgs_x_value(nr) = vx_m(i,ny-1)

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

                        k = k+1
                        lgs_a_value(k)  = 1.0   ! diagonal element only
                        lgs_a_index(k)  = nr

                        lgs_b_value(nr) = vx_m(i,j-1)
                        lgs_x_value(nr) = vx_m(i,j-1)

                    case("periodic")
                        ! Periodic boundary, take velocity from the lower boundary
                    
                        k = k+1
                        lgs_a_value(k)  = 1.0   ! diagonal element only
                        lgs_a_index(k)  = nr

                        lgs_b_value(nr) = vx_m(i,2)
                        lgs_x_value(nr) = vx_m(i,2)

                    case DEFAULT 

                        write(*,*) "calc_vxy_ssa_matrix:: Error: upper-border condition not &
                        &recognized: "//trim(boundaries_vx(2))
                        write(*,*) "boundaries parameter set to: "//trim(boundaries)
                        stop

                end select 

            else if (  ( is_front_1(i,j).and.is_front_2(i+1,j) ) &
                      .or. &
                      ( is_front_2(i,j).and.is_front_1(i+1,j) ) &
                    ) then
                    ! one neighbour is ice-covered and the other is ice-free
                    ! (calving front, grounded ice front)

                    if (is_front_1(i,j)) then
                        i1 = i     ! ice-front marker
                    else   ! is_front_1(i+1,j)==.true.
                        i1 = i+1   ! ice-front marker 
                    end if
    
                    if (.not.( is_front_2(i1-1,j) .and. is_front_2(i1+1,j) ) ) then
                        ! discretization of the x-component of the BC

                        nc = 2*ij2n(i1-1,j)-1
                               ! smallest nc (column counter), for vx_m(i1-1,j)
                        k  = k+1
                        lgs_a_value(k) = -4.0_prec*inv_dxi*vis_int_g(i1,j)
                        lgs_a_index(k) = nc

                        nc = 2*ij2n(i1,j-1)
                               ! next nc (column counter), for vy_m(i1,j-1)
                        k  = k+1
                        lgs_a_value(k) = -2.0_prec*inv_deta*vis_int_g(i1,j)
                        lgs_a_index(k) = nc

                        nc = 2*ij2n(i1,j)-1
                               ! next nc (column counter), for vx_m(i1,j)
                        k  = k+1
                        lgs_a_value(k) = 4.0_prec*inv_dxi*vis_int_g(i1,j)
                        lgs_a_index(k) = nc

                        nc = 2*ij2n(i1,j)
                               ! largest nc (column counter), for vy_m(i1,j)
                        k  = k+1
                        lgs_a_value(k) = 2.0_prec*inv_deta*vis_int_g(i1,j)
                        lgs_a_index(k) = nc

                        !lgs_b_value(nr) = factor_rhs_2*H_ice(i1,j)*H_ice(i1,j)

                        ! =========================================================
                        ! Generalized solution for all ice fronts (floating and grounded)

                        H_ice_now = H_ice(i1,j)
                        !H_ice_now = H_ice(i1,j)/f_ice(i1,j)

                        if (z_sl(i1,j)-z_bed(i1,j) .gt. 0.0) then 
                        ! Bed below sea level 
                          H_ocn_now = min(rho_ice/rho_sw*H_ice_now, &    ! Flotation depth 
                                          z_sl(i1,j)-z_bed(i1,j))            ! Grounded depth 

                        else 
                        ! Bed above sea level 
                        H_ocn_now = 0.0 

                        end if 

                        lgs_b_value(nr) = factor_rhs_3a*H_ice_now*H_ice_now &
                                      - factor_rhs_3b*H_ocn_now*H_ocn_now

!                         if (i .eq. 80 .and. j .eq. 70) then 
!                             ! Margin point
!                             write(*,*) "front ", vx_m(i1,j), vx_m(i1-1,j), vy_m(i1,j), vy_m(i1,j-1)
!                             write(*,*) "front ",vx_m(i1,j)-vx_m(i1-1,j), vy_m(i1,j)-vy_m(i1,j-1)
!                             write(*,*) "front ",4.0_prec*inv_dxi*vis_int_g(i1,j)*(vx_m(i1,j)-vx_m(i1-1,j))
!                             write(*,*) "front ",2.0_prec*inv_deta*vis_int_g(i1,j)*(vy_m(i1,j)-vy_m(i1,j-1))
!                             write(*,*) "front ",lgs_b_value(nr), &
!                             4.0_prec*inv_dxi*vis_int_g(i1,j)*(vx_m(i1,j)-vx_m(i1-1,j)) &
!                             + 2.0_prec*inv_deta*vis_int_g(i1,j)*(vy_m(i1,j)-vy_m(i1,j-1)) &
!                             - lgs_b_value(nr)
!                             write(*,*) "front ", H_ice(i1-1,j+1), H_ice(i1,j+1), H_ice(i1+1,j+1)
!                             write(*,*) "front ", H_ice(i1-1,j),   H_ice(i1,j),   H_ice(i1+1,j)
!                             write(*,*) "front ", H_ice(i1-1,j-1), H_ice(i1,j-1), H_ice(i1+1,j-1)
!                             write(*,*) "front "  
!                         end if 

!                         if (abs(vx_m(i,j)) .gt. 4e3) then 
!                             write(*,*) "ssaxcf:", H_ice(i1,j), H_ocn_now/H_ice(i1,j), vx_m(i,j), vis_int_g(i1,j)
!                         else if (abs(vx_m(i,j)) .lt. 0.5e3) then 
!                             write(*,*) "ssaxcs:", H_ice(i1,j), H_ocn_now/H_ice(i1,j), vx_m(i,j), vis_int_g(i1,j)
!                         end if 

                        ! =========================================================
              
                        lgs_x_value(nr) = vx_m(i,j)
              
                    else    ! (is_front_2(i1-1,j)==.true.).and.(is_front_2(i1+1,j)==.true.);
                            ! velocity assumed to be zero

                        k = k+1
                        lgs_a_value(k)  = 1.0_prec   ! diagonal element only
                        lgs_a_index(k)  = nr

                        lgs_b_value(nr) = 0.0_prec
                        lgs_x_value(nr) = 0.0_prec

                    end if 

            else if ( ( (maske(i,j)==3).and.(maske(i+1,j)==1) ) &
                      .or. &
                      ( (maske(i,j)==1).and.(maske(i+1,j)==3) ) &
                    ) then
                    ! one neighbour is floating ice and the other is ice-free land;
                    ! velocity assumed to be zero

                    k = k+1
                    lgs_a_value(k)  = 1.0_prec   ! diagonal element only
                    lgs_a_index(k)  = nr

                    lgs_b_value(nr) = 0.0_prec
                    lgs_x_value(nr) = 0.0_prec

            else if (ssa_mask_acx(i,j) .eq. 0) then    ! neither neighbour is floating or grounded ice,
                    ! velocity assumed to be zero

                k = k+1
                lgs_a_value(k) = 1.0_prec   ! diagonal element only
                lgs_a_index(k) = nr

                lgs_b_value(nr) = 0.0_prec
                lgs_x_value(nr) = 0.0_prec

            else
                ! === Proceed with normal ssa solution =================
                ! inner point on the staggered grid in x-direction
                ! ie, if ( (i /= nx).and.(j /= 1).and.(j /= ny) ) then
             
                    ! inner shelfy stream or floating ice 

                    if (i .eq. 1) then  ! ajr: filler to avoid LIS errors 
                        k  = k+1
                        lgs_a_value(k) = 1.0 
                        lgs_a_index(k) = nr 
                    else 
                        nc = 2*ij2n(i-1,j)-1
                            ! smallest nc (column counter), for vx_m(i-1,j)
                        k = k+1
                        lgs_a_value(k) = 4.0_prec*inv_dxi2*vis_int_g(i,j)
                        lgs_a_index(k) = nc
                    end if 

                    nc = 2*ij2n(i,j-1)-1
                        ! next nc (column counter), for vx_m(i,j-1)
                    k = k+1
                    lgs_a_value(k) = inv_deta2*vis_int_sgxy(i,j-1)
                    lgs_a_index(k) = nc

                    nc = 2*ij2n(i,j-1)
                        ! next nc (column counter), for vy_m(i,j-1)
                    k = k+1
                    lgs_a_value(k) = inv_dxi_deta &
                                    *(2.0_prec*vis_int_g(i,j)+vis_int_sgxy(i,j-1))
                    lgs_a_index(k) = nc

                    nc = 2*ij2n(i,j)-1
                        ! next nc (column counter), for vx_m(i,j)
!                     if (nc /= nr) then   ! (diagonal element)
!                         errormsg = ' >>> calc_vxy_ssa_matrix: ' &
!                                      //'Check for diagonal element failed!'
!                         call error(errormsg)
!                     end if
                    k = k+1
                    lgs_a_value(k) = -4.0_prec*inv_dxi2 &
                                            *(vis_int_g(i+1,j)+vis_int_g(i,j)) &
                                     -inv_deta2 &
                                            *(vis_int_sgxy(i,j)+vis_int_sgxy(i,j-1)) &
                                     -beta_now
                    lgs_a_index(k) = nc

                    nc = 2*ij2n(i,j)
                        ! next nc (column counter), for vy_m(i,j)
                    k = k+1
                    lgs_a_value(k) = -inv_dxi_deta &
                                    *(2.0_prec*vis_int_g(i,j)+vis_int_sgxy(i,j))
                    lgs_a_index(k) = nc

                    nc = 2*ij2n(i,j+1)-1
                        ! next nc (column counter), for vx_m(i,j+1)
                    k = k+1
                    lgs_a_value(k) = inv_deta2*vis_int_sgxy(i,j)
                    lgs_a_index(k) = nc

                    nc = 2*ij2n(i+1,j-1)
                        ! next nc (column counter), for vy_m(i+1,j-1)
                    k  = k+1
                    lgs_a_value(k) = -inv_dxi_deta &
                                  *(2.0_prec*vis_int_g(i+1,j)+vis_int_sgxy(i,j-1))
                    lgs_a_index(k) = nc

                    nc = 2*ij2n(i+1,j)-1
                        ! next nc (column counter), for vx_m(i+1,j)
                    k = k+1
                    lgs_a_value(k) = 4.0_prec*inv_dxi2*vis_int_g(i+1,j)
                    lgs_a_index(k) = nc

                    nc = 2*ij2n(i+1,j)
                        ! largest nc (column counter), for vy_m(i+1,j)
                    k  = k+1
                    lgs_a_value(k) = inv_dxi_deta &
                                    *(2.0_prec*vis_int_g(i+1,j)+vis_int_sgxy(i,j))
                    lgs_a_index(k) = nc

                    lgs_b_value(nr) = taud_now
                    lgs_x_value(nr) = vx_m(i,j)

            end if

            lgs_a_ptr(nr+1) = k+1   ! row is completed, store index to next row

            !  ------ Equations for vy_m (at (i,j+1/2))

            nr = n+1   ! row counter

            
            ! Set current boundary variables for later access
            beta_now  = beta_acy(i,j) 
            taud_now  = taud_acy(i,j) 

            ! == Treat special cases first ==

            if (ssa_mask_acy(i,j) .eq. -1) then 
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
                        lgs_a_value(k)  = 1.0_prec   ! diagonal element only
                        lgs_a_index(k)  = nr

                        lgs_b_value(nr) = 0.0_prec
                        lgs_x_value(nr) = 0.0_prec

                    case("infinite")
                        ! Infinite boundary, take velocity from one point inward
                        
                        k = k+1
                        lgs_a_value(k)  = 1.0   ! diagonal element only
                        lgs_a_index(k)  = nr

                        lgs_b_value(nr) = vy_m(i,2)
                        lgs_x_value(nr) = vy_m(i,2)
                        
                    case("periodic")
                        ! Periodic boundary, take velocity from the opposite boundary
                        
                        k = k+1
                        lgs_a_value(k)  = 1.0   ! diagonal element only
                        lgs_a_index(k)  = nr

                        lgs_b_value(nr) = vy_m(i,ny-2)
                        lgs_x_value(nr) = vy_m(i,ny-2)
                        
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
                        lgs_a_value(k)  = 1.0_prec   ! diagonal element only
                        lgs_a_index(k)  = nr

                        lgs_b_value(nr) = 0.0_prec
                        lgs_x_value(nr) = 0.0_prec

                    case("infinite")
                        ! Infinite boundary, take velocity from two points inward
                        ! (to account for staggering)

                        k = k+1
                        lgs_a_value(k)  = 1.0   ! diagonal element only
                        lgs_a_index(k)  = nr

                        lgs_b_value(nr) = vy_m(i,ny-2)
                        lgs_x_value(nr) = vy_m(i,ny-2)
                        
                    case("periodic")
                        ! Periodic boundary, take velocity from the right boundary
                        
                        k = k+1
                        lgs_a_value(k)  = 1.0   ! diagonal element only
                        lgs_a_index(k)  = nr

                        lgs_b_value(nr) = vy_m(i,3)
                        lgs_x_value(nr) = vy_m(i,3)
                        
                    case DEFAULT 

                        write(*,*) "calc_vxy_ssa_matrix:: Error: upper-border condition not &
                        &recognized: "//trim(boundaries_vy(2))
                        write(*,*) "boundaries parameter set to: "//trim(boundaries)
                        stop

                end select 
            
            else if (j .eq. ny-1 .and. trim(boundaries_vy(2)) .eq. "periodic") then
                ! Upper boundary, inward by one point
                ! (only needed for periodic conditions, otherwise
                ! this point should be treated as normal)
                
                ! Periodic boundary, take velocity from the lower boundary
                
                k = k+1
                lgs_a_value(k)  = 1.0   ! diagonal element only
                lgs_a_index(k)  = nr

                lgs_b_value(nr) = vy_m(i,2)
                lgs_x_value(nr) = vy_m(i,2)
                    
            else if (i .eq. 1) then 
                ! Left boundary 

                select case(trim(boundaries_vy(3)))

                    case("zeros")
                        ! Assume border velocity is zero 

                        k = k+1
                        lgs_a_value(k)  = 1.0_prec   ! diagonal element only
                        lgs_a_index(k)  = nr

                        lgs_b_value(nr) = 0.0_prec
                        lgs_x_value(nr) = 0.0_prec

                    case("infinite")
                        ! Infinite boundary, take velocity from one point inward

                        k = k+1
                        lgs_a_value(k)  = 1.0   ! diagonal element only
                        lgs_a_index(k)  = nr

                        lgs_b_value(nr) = vy_m(i,2)
                        lgs_x_value(nr) = vy_m(i,2)
                        
                    case("periodic")
                        ! Periodic boundary, take velocity from the right boundary
                        
                        k = k+1
                        lgs_a_value(k)  = 1.0   ! diagonal element only
                        lgs_a_index(k)  = nr

                        lgs_b_value(nr) = vy_m(nx-1,j)
                        lgs_x_value(nr) = vy_m(nx-1,j)
                        
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
                        lgs_a_value(k)  = 1.0_prec   ! diagonal element only
                        lgs_a_index(k)  = nr

                        lgs_b_value(nr) = 0.0_prec
                        lgs_x_value(nr) = 0.0_prec

                    case("infinite")
                        ! Infinite boundary, take velocity from one point inward

                        k = k+1
                        lgs_a_value(k)  = 1.0   ! diagonal element only
                        lgs_a_index(k)  = nr

                        lgs_b_value(nr) = vy_m(nx-1,j)
                        lgs_x_value(nr) = vy_m(nx-1,j)
                        
                    case("periodic")
                        ! Periodic boundary, take velocity from the right boundary
                        
                        k = k+1
                        lgs_a_value(k)  = 1.0   ! diagonal element only
                        lgs_a_index(k)  = nr

                        lgs_b_value(nr) = vy_m(2,j)
                        lgs_x_value(nr) = vy_m(2,j)
                        
                    case DEFAULT 

                        write(*,*) "calc_vxy_ssa_matrix:: Error: right-border condition not &
                        &recognized: "//trim(boundaries_vy(1))
                        write(*,*) "boundaries parameter set to: "//trim(boundaries)
                        stop

                end select 
            
            else if (  ( is_front_1(i,j).and.is_front_2(i,j+1) ) &
                      .or. &
                      ( is_front_2(i,j).and.is_front_1(i,j+1) ) &
                    ) then
                    ! one neighbour is ice-covered and the other is ice-free
                    ! (calving front, grounded ice front)

                    if (is_front_1(i,j)) then
                        j1 = j     ! ice-front marker
                    else   ! is_front_1(i,j+1)==.true.
                        j1 = j+1   ! ice-front marker
                    end if

                    if (.not.( is_front_2(i,j1-1) .and. is_front_2(i,j1+1) ) ) then
                        ! discretization of the y-component of the BC

                        nc = 2*ij2n(i-1,j1)-1
                            ! smallest nc (column counter), for vx_m(i-1,j1)
                        k = k+1
                        lgs_a_value(k) = -2.0_prec*inv_dxi*vis_int_g(i,j1)
                        lgs_a_index(k) = nc

                        nc = 2*ij2n(i,j1-1)
                            ! next nc (column counter), for vy_m(i,j1-1)
                        k = k+1
                        lgs_a_value(k) = -4.0_prec*inv_deta*vis_int_g(i,j1)
                        lgs_a_index(k) = nc

                        nc = 2*ij2n(i,j1)-1
                            ! next nc (column counter), for vx_m(i,j1)
                        k = k+1
                        lgs_a_value(k) = 2.0_prec*inv_dxi*vis_int_g(i,j1)
                        lgs_a_index(k) = nc

                        nc = 2*ij2n(i,j1)
                            ! largest nc (column counter), for vy_m(i,j1)
                        k = k+1
                        lgs_a_value(k) = 4.0_prec*inv_deta*vis_int_g(i,j1)
                        lgs_a_index(k) = nc

!                         lgs_b_value(nr) = factor_rhs_2*H_ice(i,j1)*H_ice(i,j1)

                        ! =========================================================
                        ! Generalized solution for all ice fronts (floating and grounded)
            
                        H_ice_now = H_ice(i,j1)
                        !H_ice_now = H_ice(i,j1)/f_ice(i,j1)

                        if (z_sl(i,j1)-z_bed(i,j1) .gt. 0.0) then 
                            ! Bed below sea level 
                            H_ocn_now = min(rho_ice/rho_sw*H_ice_now, &    ! Flotation depth 
                                          z_sl(i,j1)-z_bed(i,j1))            ! Grounded depth 

                        else 
                            ! Bed above sea level 
                            H_ocn_now = 0.0 

                        end if 

                        lgs_b_value(nr) = factor_rhs_3a*H_ice_now*H_ice_now &
                                        - factor_rhs_3b*H_ocn_now*H_ocn_now  

                        ! =========================================================
              
                        lgs_x_value(nr) = vy_m(i,j)
             
                    else    ! (is_front_2(i,j1-1)==.true.).and.(is_front_2(i,j1+1)==.true.);
                            ! velocity assumed to be zero

                        k = k+1
                        lgs_a_value(k)  = 1.0_prec   ! diagonal element only
                        lgs_a_index(k)  = nr

                        lgs_b_value(nr) = 0.0_prec
                        lgs_x_value(nr) = 0.0_prec

                    end if

            else if ( ( (maske(i,j)==3).and.(maske(i,j+1)==1) ) &
                        .or. &
                        ( (maske(i,j)==1).and.(maske(i,j+1)==3) ) &
                      ) then
                    ! one neighbour is floating ice and the other is ice-free land;
                    ! velocity assumed to be zero

                    k = k+1
                    lgs_a_value(k)  = 1.0_prec   ! diagonal element only
                    lgs_a_index(k)  = nr

                    lgs_b_value(nr) = 0.0_prec
                    lgs_x_value(nr) = 0.0_prec

            else if (ssa_mask_acy(i,j) .eq. 0) then    ! neither neighbour is floating or grounded ice,
                    ! velocity assumed to be zero

                k = k+1
                lgs_a_value(k)  = 1.0_prec   ! diagonal element only
                lgs_a_index(k)  = nr

                lgs_b_value(nr) = 0.0_prec
                lgs_x_value(nr) = 0.0_prec

            else
                ! === Proceed with normal ssa solution =================
                ! inner point on the staggered grid in y-direction
                ! ie, if ( (j /= ny).and.(i /= 1).and.(i /= nx) ) then
             
                    ! inner shelfy stream or floating ice 

                    nc = 2*ij2n(i-1,j)-1
                        ! smallest nc (column counter), for vx_m(i-1,j)
                    k = k+1
                    lgs_a_value(k) = inv_dxi_deta &
                                        *(2.0_prec*vis_int_g(i,j)+vis_int_sgxy(i-1,j))
                    lgs_a_index(k) = nc

                    nc = 2*ij2n(i-1,j)
                        ! next nc (column counter), for vy_m(i-1,j)
                    k = k+1
                    lgs_a_value(k) = inv_dxi2*vis_int_sgxy(i-1,j)
                    lgs_a_index(k) = nc

                    nc = 2*ij2n(i-1,j+1)-1
                        ! next nc (column counter), for vx_m(i-1,j+1)
                    k = k+1
                    lgs_a_value(k) = -inv_dxi_deta &
                                          *(2.0_prec*vis_int_g(i,j+1)+vis_int_sgxy(i-1,j))
                    lgs_a_index(k) = nc

                    if (j .eq. 1) then  ! ajr: filler to avoid LIS errors 
                        k  = k+1
                        lgs_a_value(k) = 1.0   ! diagonal element only
                        lgs_a_index(k) = nr
                    else
                        nc = 2*ij2n(i,j-1)
                            ! next nc (column counter), for vy_m(i,j-1)
                        k = k+1
                        lgs_a_value(k) = 4.0_prec*inv_deta2*vis_int_g(i,j)
                        lgs_a_index(k) = nc
                    end if 

                    nc = 2*ij2n(i,j)-1
                        ! next nc (column counter), for vx_m(i,j)
                    k = k+1
                    lgs_a_value(k) = -inv_dxi_deta &
                                            *(2.0_prec*vis_int_g(i,j)+vis_int_sgxy(i,j))
                    lgs_a_index(k) = nc

                    nc = 2*ij2n(i,j)
                        ! next nc (column counter), for vy_m(i,j)
!                     if (nc /= nr) then   ! (diagonal element)
!                         errormsg = ' >>> calc_vxy_ssa_matrix: ' &
!                                     //'Check for diagonal element failed!'
!                         call error(errormsg)
!                     end if
                    k = k+1
                    lgs_a_value(k) = -4.0_prec*inv_deta2 &
                                        *(vis_int_g(i,j+1)+vis_int_g(i,j)) &
                                     -inv_dxi2 &
                                        *(vis_int_sgxy(i,j)+vis_int_sgxy(i-1,j)) &
                                     -beta_now
                    lgs_a_index(k) = nc

                    nc = 2*ij2n(i,j+1)-1
                        ! next nc (column counter), for vx_m(i,j+1)
                    k = k+1
                    lgs_a_value(k) = inv_dxi_deta &
                                    *(2.0_prec*vis_int_g(i,j+1)+vis_int_sgxy(i,j))
                    lgs_a_index(k) = nc

                    nc = 2*ij2n(i,j+1)
                        ! next nc (column counter), for vy_m(i,j+1)
                    k = k+1
                    lgs_a_value(k) = 4.0_prec*inv_deta2*vis_int_g(i,j+1)
                    lgs_a_index(k) = nc

                    nc = 2*ij2n(i+1,j)
                        ! largest nc (column counter), for vy_m(i+1,j)
                    k = k+1
                    lgs_a_value(k)  = inv_dxi2*vis_int_sgxy(i,j)
                    lgs_a_index(k)  = nc

                    lgs_b_value(nr) = taud_now 
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

        !call lis_solver_get_iter(solver, lin_iter, ierr)
        !write(6,'(a,i0,a)', advance='no') 'lin_iter = ', lin_iter, ', '

        !!! call lis_solver_get_time(solver,solver_time,ierr)
        !!! print *, 'calc_vxy_ssa_matrix: time (s) = ', solver_time

        ! Obtain the relative L2_norm == ||b-Ax||_2 / ||b||_2
        call lis_solver_get_residualnorm(solver,residual,ierr)
        L2_norm = real(residual,prec) 

        lgs_x_value = 0.0_prec
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

    subroutine set_ssa_masks(ssa_mask_acx,ssa_mask_acy,beta_acx,beta_acy,H_ice,f_ice,f_grnd_acx,f_grnd_acy,beta_max,use_ssa)
        ! Define where ssa calculations should be performed
        ! Note: could be binary, but perhaps also distinguish 
        ! grounding line/zone to use this mask for later gl flux corrections
        ! mask = 0: no ssa calculated
        ! mask = 1: shelfy-stream ssa calculated 
        ! mask = 2: shelf ssa calculated 

        implicit none 
        
        integer,    intent(OUT) :: ssa_mask_acx(:,:) 
        integer,    intent(OUT) :: ssa_mask_acy(:,:)
        real(prec), intent(IN)  :: beta_acx(:,:)
        real(prec), intent(IN)  :: beta_acy(:,:)
        real(prec), intent(IN)  :: H_ice(:,:)
        real(prec), intent(IN)  :: f_ice(:,:)
        real(prec), intent(IN)  :: f_grnd_acx(:,:)
        real(prec), intent(IN)  :: f_grnd_acy(:,:)
        real(prec), intent(IN)  :: beta_max
        logical,    intent(IN)  :: use_ssa       ! SSA is actually active now? 

        ! Local variables
        integer    :: i, j, nx, ny
        integer    :: im1, ip1, jm1, jp1
        real(prec) :: H_acx, H_acy
        
        nx = size(H_ice,1)
        ny = size(H_ice,2)
        
        ! Initially no active ssa points
        ssa_mask_acx = 0
        ssa_mask_acy = 0
        
        if (use_ssa) then 

            do j = 1, ny
            do i = 1, nx

                ! Get neighbor indices
                im1 = max(i-1,1) 
                ip1 = min(i+1,nx) 
                jm1 = max(j-1,1) 
                jp1 = min(j+1,ny)

                ! x-direction
                if (f_ice(i,j) .eq. 1.0 .or. f_ice(ip1,j) .eq. 1.0) then

                    ! Ice is present on ac-node
                    
                    if (f_grnd_acx(i,j) .gt. 0.0) then 
                        ! Grounded ice or grounding line (ie, shelfy-stream)
                        ssa_mask_acx(i,j) = 1
                    else 
                        ! Shelf ice 
                        ssa_mask_acx(i,j) = 2
                    end if 

                    ! Deactivate if dragging is to high and away from grounding line
                    if ( beta_acx(i,j) .ge. beta_max .and. f_grnd_acx(i,j) .eq. 1.0 ) ssa_mask_acx(i,j) = 0 
                    
                end if

                ! y-direction
                if (f_ice(i,j) .eq. 1.0 .or. f_ice(i,jp1) .eq. 1.0) then

                    ! Ice is present on ac-node
                    
                    if (f_grnd_acy(i,j) .gt. 0.0) then 
                        ! Grounded ice or grounding line (ie, shelfy-stream)
                        ssa_mask_acy(i,j) = 1
                    else 
                        ! Shelf ice 
                        ssa_mask_acy(i,j) = 2
                    end if 

                    ! Deactivate if dragging is to high and away from grounding line
                    if ( beta_acy(i,j) .ge. beta_max .and. f_grnd_acy(i,j) .eq. 1.0 ) ssa_mask_acy(i,j) = 0 
                    
                end if
                 
            end do 
            end do

            ! Final check on both masks to avoid isolated non-ssa points
            do j = 2, ny-1
            do i = 2, nx-1

                ! acx-nodes 
                if (  ssa_mask_acx(i,j) .eq. 0 .and. &
                    ssa_mask_acx(i+1,j) .gt. 0 .and. ssa_mask_acx(i-1,j) .gt. 0 .and.  &
                    ssa_mask_acx(i,j+1) .gt. 0 .and. ssa_mask_acx(i,j-1) .gt. 0 ) then 

                    if (f_grnd_acx(i,j) .gt. 0.0) then 
                        ! Grounded ice or grounding line (ie, shelfy-stream)
                        ssa_mask_acx(i,j) = 1
                    else 
                        ! Shelf ice 
                        ssa_mask_acx(i,j) = 2
                    end if 

                end if 

                ! acy-nodes 
                if (  ssa_mask_acy(i,j) .eq. 0 .and. &
                    ssa_mask_acy(i+1,j) .gt. 0 .and. ssa_mask_acy(i-1,j) .gt. 0 .and.  &
                    ssa_mask_acy(i,j+1) .gt. 0 .and. ssa_mask_acy(i,j-1) .gt. 0 ) then 

                    if (f_grnd_acy(i,j) .gt. 0.0) then 
                        ! Shelfy stream 
                        ssa_mask_acy(i,j) = 1
                    else 
                        ! Shelf 
                        ssa_mask_acy(i,j) = 2
                    end if 

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
        real(prec), intent(IN) :: err_x(:,:) 
        real(prec), intent(IN) :: err_y(:,:) 
        real(prec), intent(IN) :: err_lim 

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

    subroutine check_vel_convergence_l1rel_matrix(err_x,err_y,ux,uy,ux_prev,uy_prev)

        implicit none 

        real(prec), intent(OUT) :: err_x(:,:)
        real(prec), intent(OUT) :: err_y(:,:)
        real(prec), intent(IN)  :: ux(:,:) 
        real(prec), intent(IN)  :: uy(:,:) 
        real(prec), intent(IN)  :: ux_prev(:,:) 
        real(prec), intent(IN)  :: uy_prev(:,:)  

        ! Local variables

        real(prec), parameter :: ssa_vel_tolerance = 1e-2   ! [m/a] only consider points with velocity above this tolerance limit
        real(prec), parameter :: tol = 1e-5 

        ! Error in x-direction
        where (abs(ux) .gt. ssa_vel_tolerance) 
            err_x = 2.0_prec * abs(ux - ux_prev) / abs(ux + ux_prev + tol)
        elsewhere 
            err_x = 0.0_prec
        end where 

        ! Error in y-direction 
        where (abs(uy) .gt. ssa_vel_tolerance) 
            err_y = 2.0_prec * abs(uy - uy_prev) / abs(uy + uy_prev + tol)
        elsewhere 
            err_y = 0.0_prec
        end where 

        return 

    end subroutine check_vel_convergence_l1rel_matrix

    function check_vel_convergence_l2rel(ux,uy,ux_prev,uy_prev,mask_acx,mask_acy, &
                                        ssa_resid_tol,iter,iter_max,log,use_L2_norm,L2_norm) result(is_converged)

        implicit none 

        real(prec), intent(IN) :: ux(:,:) 
        real(prec), intent(IN) :: uy(:,:) 
        real(prec), intent(IN) :: ux_prev(:,:) 
        real(prec), intent(IN) :: uy_prev(:,:)  
        logical,    intent(IN) :: mask_acx(:,:) 
        logical,    intent(IN) :: mask_acy(:,:) 
        real(prec), intent(IN) :: ssa_resid_tol 
        integer,    intent(IN) :: iter 
        integer,    intent(IN) :: iter_max 
        logical,    intent(IN) :: log 
        logical,    intent(IN) :: use_L2_norm           ! Use externally provided value
        real(prec), optional, intent(IN) :: L2_norm     ! Optional, externally calculated relative L2_norm 
        logical :: is_converged

        ! Local variables 
        real(prec) :: ux_resid_max 
        real(prec) :: uy_resid_max 
        real(prec) :: res1, res2, resid
        integer    :: nx_check, ny_check  
        character(len=1) :: converged_txt 

        real(prec), parameter :: ssa_vel_tolerance = 1e-2   ! [m/a] only consider points with velocity above this tolerance limit

        ! Calculate residual acoording to the L2 relative error norm
        ! (as Eq. 65 in Gagliardini et al., GMD, 2013)

        if (use_L2_norm) then 
            ! Try to use external value if available 

            if (.not. present(L2_norm)) then 
                write(*,*) "check_vel_convergence_l2rel:: Error: external L2_norm value must be provided as an &
                            &argument to set use_L2_norm=.TRUE."
                stop 
            end if 

            resid = L2_norm 

        else 
            ! Calculate our own L2 norm based on velocity solution between previous
            ! and current iteration (following SICOPOLIS implementation)

            ! Count how many points should be checked for convergence
            nx_check = count(abs(ux).gt.ssa_vel_tolerance .and. mask_acx)
            ny_check = count(abs(uy).gt.ssa_vel_tolerance .and. mask_acy)

            if ( (nx_check+ny_check) .gt. 0 ) then

                res1 = sqrt( sum((ux-ux_prev)*(ux-ux_prev),mask=abs(ux).gt.ssa_vel_tolerance .and. mask_acx) &
                           + sum((uy-uy_prev)*(uy-uy_prev),mask=abs(uy).gt.ssa_vel_tolerance .and. mask_acy) )

                res2 = sqrt( sum((ux+ux_prev)*(ux+ux_prev),mask=abs(ux).gt.ssa_vel_tolerance .and. mask_acx) &
                           + sum((uy+uy_prev)*(uy+uy_prev),mask=abs(uy).gt.ssa_vel_tolerance .and. mask_acy) )
                res2 = max(res2,1e-8)

                resid = 2.0_prec*res1/res2 

            else 
                ! No points available for comparison, set residual equal to zero 

                resid = 0.0_prec 

            end if 

        end if 

        ! Check for convergence
        if (resid .le. ssa_resid_tol) then 
            is_converged = .TRUE. 
            converged_txt = "C"
        else if (iter .eq. iter_max) then 
            is_converged = .TRUE. 
            converged_txt = "X" 
        else 
            is_converged = .FALSE. 
            converged_txt = ""
        end if 

        if (log .and. is_converged) then
            ! Write summary to log if desired and iterations have completed

            ! Also calculate maximum error magnitude for perspective
            if (nx_check .gt. 0) then 
                ux_resid_max = maxval(abs(ux-ux_prev),mask=abs(ux).gt.ssa_vel_tolerance .and. mask_acx)
            else 
                ux_resid_max = 0.0 
            end if 

            if (ny_check .gt. 0) then 
                uy_resid_max = maxval(abs(uy-uy_prev),mask=abs(uy).gt.ssa_vel_tolerance .and. mask_acy)
            else 
                uy_resid_max = 0.0 
            end if 

            ! Write summary to log
            write(*,"(a,a2,i4,g12.4,a3,2i8,2g12.4)") &
                "ssa: ", trim(converged_txt), iter, resid, " | ", nx_check, ny_check, ux_resid_max, uy_resid_max 

        end if 
        
        return 

    end function check_vel_convergence_l2rel

    elemental subroutine relax_ssa(ux,uy,ux_prev,uy_prev,rel)
        ! Relax velocity solution with previous iteration 

        implicit none 

        real(prec), intent(INOUT) :: ux
        real(prec), intent(INOUT) :: uy
        real(prec), intent(IN)    :: ux_prev
        real(prec), intent(IN)    :: uy_prev
        real(prec), intent(IN)    :: rel

        ! Apply relaxation 
        ux = rel*ux + (1.0-rel)*ux_prev 
        uy = rel*uy + (1.0-rel)*uy_prev

        return 

    end subroutine relax_ssa



! === INTERNAL ROUTINES ==== 

    subroutine stagger_visc_aa_ab(visc_ab,visc,H_ice,f_ice)

        implicit none 

        real(prec), intent(OUT) :: visc_ab(:,:) 
        real(prec), intent(IN)  :: visc(:,:) 
        real(prec), intent(IN)  :: H_ice(:,:) 
        real(prec), intent(IN)  :: f_ice(:,:) 

        ! Local variables 
        integer :: i, j, k
        integer :: im1, ip1, jm1, jp1 
        integer :: nx, ny 

        nx = size(visc,1)
        ny = size(visc,2)

        ! Initialisation
        visc_ab = 0.0_prec 

        ! Stagger viscosity only using contributions from neighbors that have ice  
        do i = 1, nx 
        do j = 1, ny 

            ! Get neighbor indices
            im1 = max(i-1,1) 
            ip1 = min(i+1,nx) 
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny) 
            
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

            if (k .gt. 0) visc_ab(i,j) = visc_ab(i,j)/real(k,prec)

        end do
        end do

        return 

    end subroutine stagger_visc_aa_ab

    subroutine set_sico_masks(maske,front1,front2,gl1,gl2,H_ice,f_ice,H_grnd)
        ! Define where ssa calculations should be performed
        ! Note: could be binary, but perhaps also distinguish 
        ! grounding line/zone to use this mask for later gl flux corrections
        ! mask = 0: Grounded ice 
        ! mask = 1: Ice-free land 
        ! mask = 2: Open ocean  
        ! mask = 3: Ice shelf 

        ! Note: this mask is defined on central Aa nodes 
        
        implicit none 
        
        integer,    intent(OUT) :: maske(:,:) 
        logical,    intent(OUT) :: front1(:,:)
        logical,    intent(OUT) :: front2(:,:)  
        logical,    intent(OUT) :: gl1(:,:)
        logical,    intent(OUT) :: gl2(:,:) 
        real(prec), intent(IN)  :: H_ice(:,:)
        real(prec), intent(IN)  :: f_ice(:,:)
        real(prec), intent(IN)  :: H_grnd(:,:)
        
        ! Local variables
        integer    :: i, j, nx, ny
        integer    :: im1, ip1, jm1, jp1 
        logical    :: is_float 
        
        logical, parameter :: disable_grounded_fronts = .TRUE. 

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

            if (f_ice(i,j) .lt. 1.0) then 
                ! Ice-free point, or only partially-covered (consider ice free)

                if (is_float) then 
                    ! Open ocean 
                    maske(i,j) = 2
                else 
                    ! Ice-free land 
                    maske(i,j) = 1 
                end if 

            else 
                ! Ice-covered point

                if (is_float) then 
                    ! Ice shelf 
                    maske(i,j) = 3
                else
                    ! Grounded ice 
                    maske(i,j) = 0 
                end if 
                
            end if 

        end do 
        end do 
        
        !-------- Detection of the calving front --------

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

          if (disable_grounded_fronts) then 
            ! Disable detection of grounded fronts for now,
            ! because it is more stable this way...

            if ( front1(i,j) .and. maske(i,j) .eq. 0 ) front1(i,j) = .FALSE. 

          end if 

        end do
        end do
        
        return
        
    end subroutine set_sico_masks
    
    elemental subroutine limit_vel(u,u_lim)
        ! Apply a velocity limit (for stability)

        implicit none 

        real(prec), intent(INOUT) :: u  
        real(prec), intent(IN)    :: u_lim

        real(prec), parameter :: tol = 1e-5

        u = min(u, u_lim)
        u = max(u,-u_lim)

        ! Also avoid underflow errors 
        if (abs(u) .lt. tol) u = 0.0 

        return 

    end subroutine limit_vel

end module solver_ssa_sico5
