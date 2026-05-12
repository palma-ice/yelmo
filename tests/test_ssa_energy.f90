program test_ssa_energy
    ! Standalone validation: confirm that the energy SSA assembler satisfies
    ! the algebraic identity
    !
    !     K_inner = - A_residual_inner * dx * dy
    !     b_inner = - taud_inner       * dx * dy
    !
    ! against the residual assembler in solver_ssa_ac.f90, on a small synthetic
    ! grid with periodic boundaries (so every cell is "inner").
    !
    ! Pass criterion: max relative error in K and absolute error in b smaller
    ! than tol. No solver call is made — this is a pure assembler-equivalence
    ! check, independent of LIS.

    use yelmo_defs, only : wp
    use solver_linear, only : linear_solver_class
    use solver_ssa_ac, only : linear_solver_matrix_ssa_ac_csr_2D
    use solver_ssa_ac_energy, only : linear_solver_matrix_ssa_ac_csr_2D_energy

    implicit none

    integer, parameter :: nx = 8, ny = 8
    real(wp), parameter :: dx = 5000.0_wp, dy = 5000.0_wp     ! [m]
    real(wp), parameter :: beta_min = 0.0_wp
    ! Yelmo's working precision is single (wp = sp), so the two assemblers'
    ! operation-order differences produce ~1e-7 relative roundoff between the
    ! algebraically-identical entries. Tolerate that — anything larger is a
    ! real stencil/sign bug.
    real(wp), parameter :: rtol   = 1.0e-5_wp
    real(wp), parameter :: atol_K = 1.0_wp                    ! absolute fallback for tiny K entries
    real(wp), parameter :: atol_b = 1.0e-3_wp

    real(wp) :: ux(nx,ny), uy(nx,ny)
    real(wp) :: beta_acx(nx,ny), beta_acy(nx,ny)
    real(wp) :: N_aa(nx,ny)
    integer  :: ssa_mask_acx(nx,ny), ssa_mask_acy(nx,ny)
    integer  :: mask_frnt(nx,ny)
    real(wp) :: H_ice(nx,ny), f_ice(nx,ny)
    real(wp) :: taud_acx(nx,ny), taud_acy(nx,ny)
    real(wp) :: taul_int_acx(nx,ny), taul_int_acy(nx,ny)

    type(linear_solver_class) :: lgs_res, lgs_eng

    real(wp) :: dxdy
    real(wp) :: K_val, A_val, expected, diff, rel_err
    real(wp) :: max_K_rel_err, max_K_abs_err, max_b_abs_err
    integer  :: nr, nc, idx
    integer  :: i, j
    integer  :: K_fail_count, b_fail_count
    logical  :: passed

    dxdy = dx * dy

    ! ---- Build a synthetic (but plausible) input state ----
    !
    ! Use spatially varying viscosity, drag, and driving stress so that the
    ! identity check exercises every coefficient pattern.
    do j = 1, ny
    do i = 1, nx
        H_ice(i,j)        = 1000.0_wp + 50.0_wp*sin(real(i,wp)) + 30.0_wp*cos(real(j,wp))
        f_ice(i,j)        = 1.0_wp
        N_aa(i,j)         = 1.0e10_wp * (1.0_wp + 0.2_wp*sin(0.5_wp*real(i,wp))*cos(0.5_wp*real(j,wp)))
        beta_acx(i,j)     = 1.0e3_wp * (1.0_wp + 0.1_wp*real(i,wp))
        beta_acy(i,j)     = 1.0e3_wp * (1.0_wp + 0.1_wp*real(j,wp))
        taud_acx(i,j)     = -5.0e4_wp + 1.0e3_wp*real(i+j,wp)
        taud_acy(i,j)     = -3.0e4_wp + 5.0e2_wp*real(i-j,wp)
        ux(i,j)           = 10.0_wp + 0.5_wp*real(i,wp)
        uy(i,j)           = 5.0_wp  + 0.3_wp*real(j,wp)
        taul_int_acx(i,j) = 0.0_wp                ! no calving front in periodic test
        taul_int_acy(i,j) = 0.0_wp
        mask_frnt(i,j)    = 0
        ssa_mask_acx(i,j) = 1                     ! interior SSA solve everywhere
        ssa_mask_acy(i,j) = 1
    end do
    end do

    ! ---- Assemble both matrices ----
    call linear_solver_matrix_ssa_ac_csr_2D(lgs_res,ux,uy,beta_acx,beta_acy,N_aa, &
                ssa_mask_acx,ssa_mask_acy,mask_frnt,H_ice,f_ice,taud_acx,taud_acy, &
                taul_int_acx,taul_int_acy,dx,dy,beta_min,"periodic","none")

    call linear_solver_matrix_ssa_ac_csr_2D_energy(lgs_eng,ux,uy,beta_acx,beta_acy,N_aa, &
                ssa_mask_acx,ssa_mask_acy,mask_frnt,H_ice,f_ice,taud_acx,taud_acy, &
                taul_int_acx,taul_int_acy,dx,dy,beta_min,"periodic","none")

    ! ---- Compare matrices entry-by-entry ----
    !
    ! Both assemblers use the same CSR ordering (same loop, same sparsity
    ! pattern), so we can step through a_value side-by-side.
    max_K_rel_err = 0.0_wp
    max_K_abs_err = 0.0_wp
    K_fail_count  = 0

    if (lgs_res%a_ptr(lgs_res%nmax+1) /= lgs_eng%a_ptr(lgs_eng%nmax+1)) then
        write(*,*) "FAIL: differing nnz between residual and energy assembler"
        write(*,*) "  residual nnz =", lgs_res%a_ptr(lgs_res%nmax+1)-1
        write(*,*) "  energy   nnz =", lgs_eng%a_ptr(lgs_eng%nmax+1)-1
        stop 1
    end if

    do nr = 1, lgs_res%nmax
        do idx = lgs_res%a_ptr(nr), lgs_res%a_ptr(nr+1)-1
            nc = lgs_res%a_index(idx)

            ! Sanity: same column index in same slot (same sparsity pattern)
            if (lgs_eng%a_index(idx) /= nc) then
                write(*,*) "FAIL: column index mismatch at row ", nr, " slot ", idx
                stop 1
            end if

            A_val    = real(lgs_res%a_value(idx), wp)
            K_val    = real(lgs_eng%a_value(idx), wp)
            expected = -A_val * dxdy
            diff     = K_val - expected

            if (abs(expected) > 0.0_wp) then
                rel_err = abs(diff) / abs(expected)
            else
                rel_err = abs(diff)
            end if

            if (rel_err > max_K_rel_err) max_K_rel_err = rel_err
            if (abs(diff) > max_K_abs_err) max_K_abs_err = abs(diff)

            if (rel_err > rtol .and. abs(diff) > atol_K) then
                K_fail_count = K_fail_count + 1
                if (K_fail_count <= 5) then
                    write(*,'(a,i6,a,i6,a,3es14.6)') "  K mismatch nr=", nr, &
                            " nc=", nc, &
                            " (K, -A*dxdy, diff) = ", K_val, expected, diff
                end if
            end if
        end do
    end do

    ! ---- Compare RHS ----
    max_b_abs_err = 0.0_wp
    b_fail_count  = 0
    do nr = 1, lgs_res%nmax
        expected = -real(lgs_res%b_value(nr), wp) * dxdy
        diff     = real(lgs_eng%b_value(nr), wp) - expected
        if (abs(diff) > max_b_abs_err) max_b_abs_err = abs(diff)
        if (abs(diff) > atol_b) then
            b_fail_count = b_fail_count + 1
            if (b_fail_count <= 5) then
                write(*,'(a,i6,a,3es14.6)') "  b mismatch nr=", nr, &
                        " (b_eng, -b_res*dxdy, diff) = ", &
                        real(lgs_eng%b_value(nr), wp), expected, diff
            end if
        end if
    end do

    ! ---- Report ----
    write(*,*) "==============================================="
    write(*,*) " SSA energy-vs-residual algebraic-identity test"
    write(*,*) "==============================================="
    write(*,'(a,2(1x,i0))')   "  grid (nx,ny)     :", nx, ny
    write(*,'(a,2(1x,es10.3))')"  (dx,dy)          :", dx, dy
    write(*,'(a,1x,i0)')      "  total rows       :", lgs_res%nmax
    write(*,'(a,1x,i0)')      "  total nnz        :", lgs_res%a_ptr(lgs_res%nmax+1)-1
    write(*,*)
    write(*,'(a,es12.4)')     "  max K rel err    :", max_K_rel_err
    write(*,'(a,es12.4)')     "  max K abs err    :", max_K_abs_err
    write(*,'(a,es12.4)')     "  max b abs err    :", max_b_abs_err
    write(*,'(a,i0)')         "  K mismatches     :", K_fail_count
    write(*,'(a,i0)')         "  b mismatches     :", b_fail_count
    write(*,*)

    passed = (K_fail_count == 0) .and. (b_fail_count == 0)
    if (passed) then
        write(*,*) " PASS"
    else
        write(*,*) " FAIL"
        stop 1
    end if

end program test_ssa_energy
