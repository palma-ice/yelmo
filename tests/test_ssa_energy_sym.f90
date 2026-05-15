program test_ssa_energy_sym
    ! Verify K = K^T for the energy SSA assembler when given symmetric inputs.
    !
    ! Background: the energy formulation builds a Hessian K of the discrete
    ! viscous-energy density. K is claimed to be symmetric positive-definite
    ! and is solved by CG. The existing test_ssa_energy.f90 verifies the
    ! algebraic identity K_inner = -A_residual_inner * dx*dy entry-by-entry,
    ! but does NOT verify K = K^T directly.
    !
    ! When K is asymmetric, CG still converges (since the system has a unique
    ! solution if any), but the converged solution does not minimise the
    ! intended energy and inherits a directional bias from the asymmetric
    ! coupling. CalvingMIP exp1 with ssa_solver="energy" exhibits a strong
    ! SW-heavy ice-thickness asymmetry that is absent with ssa_solver="residual".
    !
    ! This program walks every nonzero entry K(I,J) and checks K(J,I) is
    ! present with the same value. It runs three cases:
    !
    !   1. periodic boundaries, all interior cells (mask=1)
    !   2. zero (no-slip) boundaries, all interior cells
    !   3. periodic boundaries with a mask=3 calving-front strip
    !
    ! Inputs are uniform (constant N_aa, beta, taud) so the only place
    ! asymmetry can come from is the assembler itself.

    use yelmo_defs, only : wp
    use solver_linear, only : linear_solver_class
    use solver_ssa_ac_energy, only : linear_solver_matrix_ssa_ac_csr_2D_energy

    implicit none

    integer, parameter :: nx = 9, ny = 9
    real(wp), parameter :: dx = 5000.0_wp, dy = 5000.0_wp
    real(wp), parameter :: beta_min = 0.0_wp
    real(wp), parameter :: tol = 1.0e-3_wp     ! absolute tolerance on |K(I,J) - K(J,I)|

    real(wp) :: ux(nx,ny), uy(nx,ny)
    real(wp) :: beta_acx(nx,ny), beta_acy(nx,ny)
    real(wp) :: N_aa(nx,ny)
    integer  :: ssa_mask_acx(nx,ny), ssa_mask_acy(nx,ny)
    integer  :: mask_frnt(nx,ny)
    real(wp) :: H_ice(nx,ny), f_ice(nx,ny)
    real(wp) :: taud_acx(nx,ny), taud_acy(nx,ny)
    real(wp) :: taul_int_acx(nx,ny), taul_int_acy(nx,ny)

    type(linear_solver_class) :: lgs

    integer :: total_fail
    integer :: i, j

    ! ---- Uniform symmetric inputs ----
    H_ice         = 1000.0_wp
    f_ice         = 1.0_wp
    N_aa          = 1.0e10_wp
    beta_acx      = 1.0e3_wp
    beta_acy      = 1.0e3_wp
    taud_acx      = -5.0e4_wp
    taud_acy      = -5.0e4_wp
    ux            = 10.0_wp
    uy            = 10.0_wp
    taul_int_acx  = 0.0_wp
    taul_int_acy  = 0.0_wp
    mask_frnt     = 0
    ssa_mask_acx  = 1
    ssa_mask_acy  = 1

    total_fail = 0

    ! ---------- Case 1: periodic, all interior ----------
    write(*,*) "============================================================"
    write(*,*) "Case 1: periodic boundaries, all interior cells (mask=1)"
    write(*,*) "============================================================"
    call linear_solver_matrix_ssa_ac_csr_2D_energy(lgs,ux,uy,beta_acx,beta_acy,N_aa, &
                ssa_mask_acx,ssa_mask_acy,mask_frnt,H_ice,f_ice,taud_acx,taud_acy, &
                taul_int_acx,taul_int_acy,dx,dy,beta_min,"periodic","none")
    call check_K_symmetry(lgs, total_fail)

    ! ---------- Case 2: zero (no-slip) boundaries, all interior ----------
    write(*,*)
    write(*,*) "============================================================"
    write(*,*) "Case 2: zeros (no-slip) boundaries, all interior cells"
    write(*,*) "============================================================"
    ! reset mask in case the zero boundary needs explicit edge masks; for now
    ! we leave ssa_mask=1 everywhere and let the assembler treat boundaries
    ! per the "zeros" string.
    call linear_solver_matrix_ssa_ac_csr_2D_energy(lgs,ux,uy,beta_acx,beta_acy,N_aa, &
                ssa_mask_acx,ssa_mask_acy,mask_frnt,H_ice,f_ice,taud_acx,taud_acy, &
                taul_int_acx,taul_int_acy,dx,dy,beta_min,"zeros","none")
    call check_K_symmetry(lgs, total_fail)

    ! ---------- Case 3: periodic with a mask=3 calving-front strip ----------
    write(*,*)
    write(*,*) "============================================================"
    write(*,*) "Case 3: periodic boundaries, ssa_mask=3 strip at j=ny/2"
    write(*,*) "============================================================"
    ssa_mask_acx = 1
    ssa_mask_acy = 1
    j = ny/2
    do i = 1, nx
        ssa_mask_acx(i,j) = 3
        ssa_mask_acy(i,j) = 3
    end do
    call linear_solver_matrix_ssa_ac_csr_2D_energy(lgs,ux,uy,beta_acx,beta_acy,N_aa, &
                ssa_mask_acx,ssa_mask_acy,mask_frnt,H_ice,f_ice,taud_acx,taud_acy, &
                taul_int_acx,taul_int_acy,dx,dy,beta_min,"periodic","none")
    call check_K_symmetry(lgs, total_fail)

    write(*,*)
    write(*,*) "============================================================"
    if (total_fail == 0) then
        write(*,*) " ALL CASES PASS: K is symmetric in every tested config"
    else
        write(*,'(a,i0,a)') "  FAIL across ", total_fail, " case(s)"
        stop 1
    end if

contains

    subroutine check_K_symmetry(lgs, total_fail)
        ! Numerical symmetry check: for every nonzero entry K(nr, nc) with
        ! nr /= nc, look up K(nc, nr) and compare. Missing slots are
        ! treated as 0 — what matters for CG is K's *values*, not its CSR
        ! sparsity pattern, so an asymmetric pattern with both halves
        ! numerically zero is still a symmetric matrix.
        type(linear_solver_class), intent(in) :: lgs
        integer, intent(inout) :: total_fail

        integer  :: nr, nc, idx, idx2, n_asym, n_checked
        real(wp) :: K_ij, K_ji, max_diff

        n_asym    = 0
        n_checked = 0
        max_diff  = 0.0_wp

        do nr = 1, lgs%nmax
            do idx = lgs%a_ptr(nr), lgs%a_ptr(nr+1)-1
                nc = lgs%a_index(idx)
                if (nc == nr) cycle             ! diagonal is trivially symmetric
                n_checked = n_checked + 1
                K_ij = real(lgs%a_value(idx), wp)

                ! Look up K(nc, nr); treat as 0 if no slot
                K_ji = 0.0_wp
                do idx2 = lgs%a_ptr(nc), lgs%a_ptr(nc+1)-1
                    if (lgs%a_index(idx2) == nr) then
                        K_ji = real(lgs%a_value(idx2), wp)
                        exit
                    end if
                end do

                if (abs(K_ij - K_ji) > tol) then
                    n_asym = n_asym + 1
                    if (n_asym <= 5) then
                        write(*,'(a,i0,a,i0,a,3es14.6)') "  K(", nr, ",", nc, &
                                ") vs K(nc,nr): ", K_ij, K_ji, K_ij - K_ji
                    end if
                end if

                if (abs(K_ij - K_ji) > max_diff) max_diff = abs(K_ij - K_ji)
            end do
        end do

        write(*,'(a,1x,i0)')   "  off-diag entries     :", n_checked
        write(*,'(a,1x,i0)')   "  numerically asym     :", n_asym
        write(*,'(a,es12.4)')  "  max |K(I,J)-K(J,I)|  :", max_diff
        if (n_asym /= 0) then
            write(*,*) "  -> FAIL"
            total_fail = total_fail + 1
        else
            write(*,*) "  -> pass"
        end if

    end subroutine check_K_symmetry

end program test_ssa_energy_sym
