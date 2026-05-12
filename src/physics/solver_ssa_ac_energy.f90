module solver_ssa_ac_energy
    ! SSA momentum solver, energy formulation.
    !
    ! Assembles the symmetric positive-(semi)definite Hessian K of the
    ! discretised SSA energy density
    !
    !     W = N_aa·(2 ux^2 + 2 vy^2 + 2 ux vy + 1/2 (uy + vx)^2)   ! membrane + shear
    !       + 1/2 beta (u^2 + v^2)                                 ! basal drag
    !       + rho_i g H (u s_x + v s_y)                            ! driving (gravitational PE)
    !
    ! integrated over the C-grid. With (eta, beta, H) frozen per Picard
    ! iteration the energy is quadratic in (u, v); the Hessian K is symmetric
    ! and positive (semi-)definite, so CG + AMG can be used in place of
    ! BiCGStab on the non-symmetric residual matrix.
    !
    ! Algebraic identity: for inner ac-nodes (no boundary or mask special-cases),
    !
    !     K_inner   = - A_residual_inner · dx · dy
    !     b_inner_E = - taud_inner       · dx · dy
    !
    ! where (A_residual, b_residual) is the matrix/RHS produced by
    ! linear_solver_matrix_ssa_ac_csr_2D in solver_ssa_ac.f90. The two
    ! formulations therefore yield the same physical solution at inner cells.
    ! Boundary handling differs — see notes below.
    !
    ! Staggering convention: yelmo-fortran C-grid. ux(i,j) lives at the RIGHT
    ! face of aa-cell (i,j); uy(i,j) lives at the TOP face. (Yelmo.jl uses
    ! the opposite convention with ux on the LEFT face — stencils here are
    ! re-derived from the energy density to match the Fortran convention.)

    use yelmo_defs, only : sp, dp, wp, io_unit_err, TOL, TOL_UNDERFLOW, is_equal
    use yelmo_tools, only : boundary_code, get_neighbor_indices_bc_codes
    use solver_linear
    use solver_ssa_ac, only : stagger_visc_aa_ab

    implicit none

    private
    public :: linear_solver_matrix_ssa_ac_csr_2D_energy

contains

    subroutine linear_solver_matrix_ssa_ac_csr_2D_energy(lgs,ux,uy,beta_acx,beta_acy, &
                            N_aa,ssa_mask_acx,ssa_mask_acy,mask_frnt,H_ice,f_ice,taud_acx, &
                            taud_acy,taul_int_acx,taul_int_acy,dx,dy,beta_min,boundaries,lateral_bc)
        ! Energy-formulation analogue of linear_solver_matrix_ssa_ac_csr_2D.
        ! Same argument list so the two assemblers are interchangeable from the
        ! Picard loop. Assembles K (SPD) and b such that K · [u; v] = b.

        implicit none

        type(linear_solver_class), intent(INOUT) :: lgs
        real(wp), intent(IN) :: ux(:,:)                 ! [m yr^-1] horizontal velocity x (acx-nodes)
        real(wp), intent(IN) :: uy(:,:)                 ! [m yr^-1] horizontal velocity y (acy-nodes)
        real(wp), intent(IN) :: beta_acx(:,:)           ! [Pa yr m^-1] basal friction (acx-nodes)
        real(wp), intent(IN) :: beta_acy(:,:)           ! [Pa yr m^-1] basal friction (acy-nodes)
        real(wp), intent(IN) :: N_aa(:,:)               ! [Pa yr m] vertically integrated viscosity (aa-nodes)
        integer,  intent(IN) :: ssa_mask_acx(:,:)       ! [--] ssa solver action mask (acx-nodes)
        integer,  intent(IN) :: ssa_mask_acy(:,:)       ! [--] ssa solver action mask (acy-nodes)
        integer,  intent(IN) :: mask_frnt(:,:)          ! [--] ice-front mask
        real(wp), intent(IN) :: H_ice(:,:)              ! [m]  ice thickness (aa-nodes)
        real(wp), intent(IN) :: f_ice(:,:)
        real(wp), intent(IN) :: taud_acx(:,:)           ! [Pa] driving stress (acx-nodes)
        real(wp), intent(IN) :: taud_acy(:,:)           ! [Pa] driving stress (acy-nodes)
        real(wp), intent(IN) :: taul_int_acx(:,:)       ! [Pa m] vertically integrated lateral stress (acx-nodes)
        real(wp), intent(IN) :: taul_int_acy(:,:)       ! [Pa m] vertically integrated lateral stress (acy-nodes)
        real(wp), intent(IN) :: dx, dy
        real(wp), intent(IN) :: beta_min                ! [Pa yr m^-1] minimum allowed basal friction
        character(len=*), intent(IN) :: boundaries
        character(len=*), intent(IN) :: lateral_bc

        ! Local variables
        integer  :: nx, ny
        integer  :: i, j, k, n
        integer  :: nc, nr
        integer  :: im1, ip1, jm1, jp1
        real(wp) :: dxdy
        real(wp) :: dyodx, dxody          ! dy/dx and dx/dy (recurring stencil prefactors)
        real(wp) :: beta_now
        real(wp), allocatable :: N_ab(:,:)

        ! Boundary conditions counterclockwise unit circle:
        ! 1: x, right border; 2: y, upper; 3: x, left; 4: y, lower
        character(len=56) :: bcs(4)

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        ! Initialise the lgs object if needed (n_terms=9 matches the residual assembler)
        if (.not. allocated(lgs%x_value)) then
            call linear_solver_init(lgs,nx,ny,nvar=2,n_terms=9)
        end if

        ! Define border conditions (only choices: no-slip, free-slip, periodic)
        select case(trim(boundaries))
            case("MISMIP3D")
                bcs(1) = "free-slip"; bcs(2) = "periodic"
                bcs(3) = "no-slip";   bcs(4) = "periodic"
            case("TROUGH")
                bcs(1) = "free-slip"; bcs(2) = "periodic"
                bcs(3) = "no-slip";   bcs(4) = "periodic"
            case("periodic")
                bcs(1:4) = "periodic"
            case("periodic-x")
                bcs(1) = "periodic";  bcs(2) = "free-slip"
                bcs(3) = "periodic";  bcs(4) = "free-slip"
            case("periodic-y")
                bcs(1) = "free-slip"; bcs(2) = "periodic"
                bcs(3) = "free-slip"; bcs(4) = "periodic"
            case("infinite")
                bcs(1:4) = "free-slip"
            case("zeros")
                bcs(1:4) = "no-slip"
            case DEFAULT
                bcs(1:4) = "no-slip"
        end select

        allocate(N_ab(nx,ny))

        ! Stencil prefactors (cell area times inverse spacings)
        dxdy  = dx * dy
        dyodx = dy / dx
        dxody = dx / dy

        ! Stagger depth-integrated viscosity to ab-nodes
        call stagger_visc_aa_ab(N_ab,N_aa,H_ice,f_ice,boundaries)

        !-------- Assembly of K · [u; v] = b in compressed-sparse-row format --------

        lgs%a_ptr(1) = 1
        k = 0

        do n = 1, lgs%nmax-1, 2

            i = lgs%n2i((n+1)/2)
            j = lgs%n2j((n+1)/2)

            ! Periodic neighbour indices (other BC types handled below)
            im1 = i-1; if (im1 .eq. 0)    im1 = nx
            ip1 = i+1; if (ip1 .eq. nx+1) ip1 = 1
            jm1 = j-1; if (jm1 .eq. 0)    jm1 = ny
            jp1 = j+1; if (jp1 .eq. ny+1) jp1 = 1

            ! ============================================================
            ! Equation for ux at acx-node (i,j)
            ! ============================================================

            nr = n   ! row counter

            if (ssa_mask_acx(i,j) .eq. 0) then
                ! Dirichlet: u = 0. Penalty form (κ on diagonal => K·u = 0 at this row).
                k = k+1
                lgs%a_value(k) = 1.0_wp
                lgs%a_index(k) = nr
                lgs%b_value(nr) = 0.0_wp
                lgs%x_value(nr) = 0.0_wp

            else if (ssa_mask_acx(i,j) .eq. -1) then
                ! Dirichlet to prescribed value.
                k = k+1
                lgs%a_value(k) = 1.0_wp
                lgs%a_index(k) = nr
                lgs%b_value(nr) = ux(i,j)
                lgs%x_value(nr) = ux(i,j)

            else if (i .eq. 1 .and. trim(bcs(3)) .ne. "periodic") then
                ! Left domain boundary
                if (trim(bcs(3)) .eq. "free-slip") then
                    nc = 2*lgs%ij2n(i,j)-1                 ! ux(i,j)
                    k = k+1
                    lgs%a_value(k) =  1.0_wp; lgs%a_index(k) = nc
                    nc = 2*lgs%ij2n(ip1,j)-1               ! ux(ip1,j)
                    k = k+1
                    lgs%a_value(k) = -1.0_wp; lgs%a_index(k) = nc
                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = ux(i,j)
                else                                       ! no-slip
                    k = k+1
                    lgs%a_value(k) = 1.0_wp; lgs%a_index(k) = nr
                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = 0.0_wp
                end if

            else if (i .eq. nx .and. trim(bcs(1)) .ne. "periodic") then
                ! Right domain boundary
                if (trim(bcs(1)) .eq. "free-slip") then
                    nc = 2*lgs%ij2n(i,j)-1                 ! ux(i,j)
                    k = k+1
                    lgs%a_value(k) =  1.0_wp; lgs%a_index(k) = nc
                    nc = 2*lgs%ij2n(nx-1,j)-1              ! ux(nx-1,j)
                    k = k+1
                    lgs%a_value(k) = -1.0_wp; lgs%a_index(k) = nc
                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = ux(i,j)
                else
                    k = k+1
                    lgs%a_value(k) = 1.0_wp; lgs%a_index(k) = nr
                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = 0.0_wp
                end if

            else if (j .eq. 1 .and. trim(bcs(4)) .ne. "periodic") then
                ! Lower domain boundary
                if (trim(bcs(4)) .eq. "free-slip") then
                    nc = 2*lgs%ij2n(i,j)-1
                    k = k+1
                    lgs%a_value(k) =  1.0_wp; lgs%a_index(k) = nc
                    nc = 2*lgs%ij2n(i,jp1)-1
                    k = k+1
                    lgs%a_value(k) = -1.0_wp; lgs%a_index(k) = nc
                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = ux(i,j)
                else
                    k = k+1
                    lgs%a_value(k) = 1.0_wp; lgs%a_index(k) = nr
                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = 0.0_wp
                end if

            else if (j .eq. ny .and. trim(bcs(2)) .ne. "periodic") then
                ! Upper domain boundary
                if (trim(bcs(2)) .eq. "free-slip") then
                    nc = 2*lgs%ij2n(i,j)-1
                    k = k+1
                    lgs%a_value(k) =  1.0_wp; lgs%a_index(k) = nc
                    nc = 2*lgs%ij2n(i,ny-1)-1
                    k = k+1
                    lgs%a_value(k) = -1.0_wp; lgs%a_index(k) = nc
                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = ux(i,j)
                else
                    k = k+1
                    lgs%a_value(k) = 1.0_wp; lgs%a_index(k) = nr
                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = 0.0_wp
                end if

            else
                ! Inner ac-node OR lateral boundary (mask=3): use the energy interior stencil.
                !
                ! The membrane contributions from the ice-free side vanish naturally because
                ! N_aa(neighbour)=0 there. The calving-front stress enters as a linear
                ! boundary-work term in b only (added below if mask=3).

                beta_now = beta_acx(i,j)
                if (ssa_mask_acx(i,j) .eq. 1 .and. beta_acx(i,j) .eq. 0.0_wp) beta_now = beta_min

                ! -- ux self terms --

                nc = 2*lgs%ij2n(i,j)-1                                 ! ux(i,j)  [diagonal]
                k = k+1
                lgs%a_value(k) =  4.0_wp*dyodx*(N_aa(i,j)+N_aa(ip1,j)) &
                                + dxody*(N_ab(i,j)+N_ab(i,jm1))         &
                                + beta_now*dxdy
                lgs%a_index(k) = nc

                nc = 2*lgs%ij2n(ip1,j)-1                                ! ux(ip1,j)
                k = k+1
                lgs%a_value(k) = -4.0_wp*dyodx*N_aa(ip1,j)
                lgs%a_index(k) = nc

                nc = 2*lgs%ij2n(im1,j)-1                                ! ux(im1,j)
                k = k+1
                lgs%a_value(k) = -4.0_wp*dyodx*N_aa(i,j)
                lgs%a_index(k) = nc

                nc = 2*lgs%ij2n(i,jp1)-1                                ! ux(i,jp1)
                k = k+1
                lgs%a_value(k) = -dxody*N_ab(i,j)
                lgs%a_index(k) = nc

                nc = 2*lgs%ij2n(i,jm1)-1                                ! ux(i,jm1)
                k = k+1
                lgs%a_value(k) = -dxody*N_ab(i,jm1)
                lgs%a_index(k) = nc

                ! -- uy cross terms --

                nc = 2*lgs%ij2n(i,j)                                    ! uy(i,j)
                k = k+1
                lgs%a_value(k) =  2.0_wp*N_aa(i,j) + N_ab(i,j)
                lgs%a_index(k) = nc

                nc = 2*lgs%ij2n(ip1,j)                                  ! uy(ip1,j)
                k = k+1
                lgs%a_value(k) = -2.0_wp*N_aa(ip1,j) - N_ab(i,j)
                lgs%a_index(k) = nc

                nc = 2*lgs%ij2n(ip1,jm1)                                ! uy(ip1,jm1)
                k = k+1
                lgs%a_value(k) =  2.0_wp*N_aa(ip1,j) + N_ab(i,jm1)
                lgs%a_index(k) = nc

                nc = 2*lgs%ij2n(i,jm1)                                  ! uy(i,jm1)
                k = k+1
                lgs%a_value(k) = -2.0_wp*N_aa(i,j) - N_ab(i,jm1)
                lgs%a_index(k) = nc

                ! Right-hand side: driving stress + (optional) calving-front boundary work
                lgs%b_value(nr) = -taud_acx(i,j)*dxdy
                if (ssa_mask_acx(i,j) .eq. 3) then
                    lgs%b_value(nr) = lgs%b_value(nr) + taul_int_acx(i,j)*dy
                end if
                lgs%x_value(nr) = ux(i,j)

            end if

            lgs%a_ptr(nr+1) = k+1

            ! ============================================================
            ! Equation for uy at acy-node (i,j)
            ! ============================================================

            nr = n+1   ! row counter

            if (ssa_mask_acy(i,j) .eq. 0) then
                k = k+1
                lgs%a_value(k) = 1.0_wp; lgs%a_index(k) = nr
                lgs%b_value(nr) = 0.0_wp
                lgs%x_value(nr) = 0.0_wp

            else if (ssa_mask_acy(i,j) .eq. -1) then
                k = k+1
                lgs%a_value(k) = 1.0_wp; lgs%a_index(k) = nr
                lgs%b_value(nr) = uy(i,j)
                lgs%x_value(nr) = uy(i,j)

            else if (j .eq. 1 .and. trim(bcs(4)) .ne. "periodic") then
                if (trim(bcs(4)) .eq. "free-slip") then
                    nc = 2*lgs%ij2n(i,j)
                    k = k+1
                    lgs%a_value(k) =  1.0_wp; lgs%a_index(k) = nc
                    nc = 2*lgs%ij2n(i,jp1)
                    k = k+1
                    lgs%a_value(k) = -1.0_wp; lgs%a_index(k) = nc
                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = uy(i,j)
                else
                    k = k+1
                    lgs%a_value(k) = 1.0_wp; lgs%a_index(k) = nr
                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = 0.0_wp
                end if

            else if (j .eq. ny .and. trim(bcs(2)) .ne. "periodic") then
                if (trim(bcs(2)) .eq. "free-slip") then
                    nc = 2*lgs%ij2n(i,j)
                    k = k+1
                    lgs%a_value(k) =  1.0_wp; lgs%a_index(k) = nc
                    nc = 2*lgs%ij2n(i,ny-1)
                    k = k+1
                    lgs%a_value(k) = -1.0_wp; lgs%a_index(k) = nc
                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = uy(i,j)
                else
                    k = k+1
                    lgs%a_value(k) = 1.0_wp; lgs%a_index(k) = nr
                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = 0.0_wp
                end if

            else if (i .eq. 1 .and. trim(bcs(3)) .ne. "periodic") then
                if (trim(bcs(3)) .eq. "free-slip") then
                    nc = 2*lgs%ij2n(i,j)
                    k = k+1
                    lgs%a_value(k) =  1.0_wp; lgs%a_index(k) = nc
                    nc = 2*lgs%ij2n(ip1,j)
                    k = k+1
                    lgs%a_value(k) = -1.0_wp; lgs%a_index(k) = nc
                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = uy(i,j)
                else
                    k = k+1
                    lgs%a_value(k) = 1.0_wp; lgs%a_index(k) = nr
                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = 0.0_wp
                end if

            else if (i .eq. nx .and. trim(bcs(1)) .ne. "periodic") then
                if (trim(bcs(1)) .eq. "free-slip") then
                    nc = 2*lgs%ij2n(i,j)
                    k = k+1
                    lgs%a_value(k) =  1.0_wp; lgs%a_index(k) = nc
                    nc = 2*lgs%ij2n(nx-1,j)
                    k = k+1
                    lgs%a_value(k) = -1.0_wp; lgs%a_index(k) = nc
                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = uy(i,j)
                else
                    k = k+1
                    lgs%a_value(k) = 1.0_wp; lgs%a_index(k) = nr
                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = 0.0_wp
                end if

            else
                ! Inner ac-node OR lateral boundary (mask=3): energy interior stencil.

                beta_now = beta_acy(i,j)
                if (ssa_mask_acy(i,j) .eq. 1 .and. beta_acy(i,j) .eq. 0.0_wp) beta_now = beta_min

                ! -- uy self terms --

                nc = 2*lgs%ij2n(i,j)                                    ! uy(i,j)  [diagonal]
                k = k+1
                lgs%a_value(k) =  4.0_wp*dxody*(N_aa(i,j)+N_aa(i,jp1)) &
                                + dyodx*(N_ab(i,j)+N_ab(im1,j))         &
                                + beta_now*dxdy
                lgs%a_index(k) = nc

                nc = 2*lgs%ij2n(i,jp1)                                  ! uy(i,jp1)
                k = k+1
                lgs%a_value(k) = -4.0_wp*dxody*N_aa(i,jp1)
                lgs%a_index(k) = nc

                nc = 2*lgs%ij2n(i,jm1)                                  ! uy(i,jm1)
                k = k+1
                lgs%a_value(k) = -4.0_wp*dxody*N_aa(i,j)
                lgs%a_index(k) = nc

                nc = 2*lgs%ij2n(ip1,j)                                  ! uy(ip1,j)
                k = k+1
                lgs%a_value(k) = -dyodx*N_ab(i,j)
                lgs%a_index(k) = nc

                nc = 2*lgs%ij2n(im1,j)                                  ! uy(im1,j)
                k = k+1
                lgs%a_value(k) = -dyodx*N_ab(im1,j)
                lgs%a_index(k) = nc

                ! -- ux cross terms --

                nc = 2*lgs%ij2n(i,j)-1                                  ! ux(i,j)
                k = k+1
                lgs%a_value(k) =  2.0_wp*N_aa(i,j) + N_ab(i,j)
                lgs%a_index(k) = nc

                nc = 2*lgs%ij2n(i,jp1)-1                                ! ux(i,jp1)
                k = k+1
                lgs%a_value(k) = -2.0_wp*N_aa(i,jp1) - N_ab(i,j)
                lgs%a_index(k) = nc

                nc = 2*lgs%ij2n(im1,jp1)-1                              ! ux(im1,jp1)
                k = k+1
                lgs%a_value(k) =  2.0_wp*N_aa(i,jp1) + N_ab(im1,j)
                lgs%a_index(k) = nc

                nc = 2*lgs%ij2n(im1,j)-1                                ! ux(im1,j)
                k = k+1
                lgs%a_value(k) = -2.0_wp*N_aa(i,j) - N_ab(im1,j)
                lgs%a_index(k) = nc

                lgs%b_value(nr) = -taud_acy(i,j)*dxdy
                if (ssa_mask_acy(i,j) .eq. 3) then
                    lgs%b_value(nr) = lgs%b_value(nr) + taul_int_acy(i,j)*dx
                end if
                lgs%x_value(nr) = uy(i,j)

            end if

            lgs%a_ptr(nr+1) = k+1

        end do

        return

    end subroutine linear_solver_matrix_ssa_ac_csr_2D_energy

end module solver_ssa_ac_energy
