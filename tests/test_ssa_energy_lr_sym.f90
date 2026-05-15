program test_ssa_energy_lr_sym
    ! Verify the energy SSA assembler is equivariant under L<->R (and T<->B)
    ! reflection of the C-grid. This complements test_ssa_energy_sym which
    ! checks K = K^T algebraically; here we check that K and b respect the
    ! reflection symmetry of a mirror-symmetric calving-front configuration.
    !
    ! Background: the bug fixed in commit "solver_ssa_ac_energy: flip
    ! taul_int sign at left fronts" added the boundary-work term taul_int*dy
    ! to b with a fixed +sign at every ssa_mask=3 cell. Physically the
    ! contribution is (sigma . n) . u dl and the outward normal n flips
    ! between left and right fronts, so the sign must flip too. The
    ! pre-fix assembler hands the same +taul_int*dy to both sides, breaking
    ! reflection equivariance of b. K stays mirror-symmetric on the first
    ! Picard iter but drifts as asymmetric u feeds back through visc_eff.
    !
    ! The test is direct: build a strip of ice through the middle of the
    ! grid (vertical strip for L<->R, horizontal for T<->B), turn on
    ! taul_int at the two fronts, call the assembler once, then walk every
    ! row n and check K(n, c) == sign_n * sign_c * K(n_mirror, c_mirror)
    ! and b(n) == sign_n * b(n_mirror) for every (n, c).
    !
    ! C-grid reflection P (about the vertical center line):
    !   ux(i, j) at right face of cell i  ->  ux(nx-i,   j) with sign -1
    !   uy(i, j) at top face of cell i    ->  uy(nx+1-i, j) with sign +1
    ! Mirrors for T<->B follow by exchanging i<->j and ux<->uy.

    use yelmo_defs, only : wp
    use solver_linear, only : linear_solver_class
    use solver_ssa_ac_energy, only : linear_solver_matrix_ssa_ac_csr_2D_energy

    implicit none

    integer, parameter :: nx = 11, ny = 11
    real(wp), parameter :: dx = 5000.0_wp, dy = 5000.0_wp
    real(wp), parameter :: beta_min = 0.0_wp
    real(wp), parameter :: tol_K = 1.0e-3_wp     ! abs tol on |K(n,c) - mirror|
    real(wp), parameter :: tol_b = 1.0e-3_wp     ! abs tol on |b(n) - mirror|

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

    ! ---- Uniform symmetric base inputs ----
    H_ice         = 1000.0_wp
    f_ice         = 0.0_wp                  ! Ice-free by default; cases below carve out an ice patch
    N_aa          = 1.0e10_wp
    beta_acx      = 0.0_wp                  ! Floating shelf — no basal drag
    beta_acy      = 0.0_wp
    taud_acx      = 0.0_wp                  ! Body force off so b at mask=3 isolates taul_int
    taud_acy      = 0.0_wp
    ux            = 0.0_wp
    uy            = 0.0_wp
    taul_int_acx  = 0.0_wp
    taul_int_acy  = 0.0_wp
    mask_frnt     = 0
    ssa_mask_acx  = 0
    ssa_mask_acy  = 0

    total_fail = 0

    ! ---------- Case A: vertical ice strip, L<->R reflection ----------
    ! Ice in cells i = 4..8 (mirror-symmetric about i = (nx+1)/2 = 6).
    ! LEFT front at ac-edge ux(3, j); RIGHT front at ac-edge ux(8, j).
    write(*,*)
    write(*,*) "============================================================"
    write(*,*) "Case A: vertical ice strip, check L<->R reflection of K, b"
    write(*,*) "============================================================"
    call reset_state()
    do j = 1, ny
        do i = 4, 8
            f_ice(i, j) = 1.0_wp                ! Ice cells
        end do
        ! ssa_mask_acx semantics:
        !   0 = no SSA, u=0 ;  1 = interior shelfy-stream ;  3 = lateral BC.
        ! Interior of strip: i = 4..7 are ac-edges fully inside ice.
        do i = 4, 7
            ssa_mask_acx(i, j) = 1
        end do
        ! Two fronts:
        ssa_mask_acx(3, j) = 3                  ! ice-free LEFT, ice RIGHT
        ssa_mask_acx(8, j) = 3                  ! ice LEFT, ice-free RIGHT
        ! Activate uy in the interior of the strip too so the cross-coupling
        ! between ux and uy rows is exercised. Strip spans j = 1..ny so
        ! every uy ac-node within i = 4..8 is interior.
        do i = 4, 8
            ssa_mask_acy(i, j) = 1
        end do
        ! Lateral stress: identical magnitude on both fronts. The bug shows
        ! up as b(left) == +taul ; b(right) == +taul instead of opposite.
        taul_int_acx(3, j) = 1.0e6_wp
        taul_int_acx(8, j) = 1.0e6_wp
    end do

    ! "periodic" = all four sides periodic. boundary_code only knows
    ! "periodic" and "periodic-x"; the assembler also accepts "periodic-y"
    ! and "infinite" but stagger_visc_aa_ab via boundary_code does not.
    ! Case A only needs y to be periodic (so the strip wraps top<->bottom);
    ! using fully-periodic is fine and keeps the boundary handling out of
    ! the symmetry check.
    call linear_solver_matrix_ssa_ac_csr_2D_energy(lgs, ux, uy, beta_acx, beta_acy, N_aa, &
                ssa_mask_acx, ssa_mask_acy, mask_frnt, H_ice, f_ice, taud_acx, taud_acy, &
                taul_int_acx, taul_int_acy, dx, dy, beta_min, "periodic", "none")
    call check_lr_reflection(lgs, total_fail)

    ! ---------- Case B: horizontal ice strip, T<->B reflection ----------
    ! Symmetric to Case A with i and j roles swapped; tests the uy mask=3 branch.
    write(*,*)
    write(*,*) "============================================================"
    write(*,*) "Case B: horizontal ice strip, check T<->B reflection of K, b"
    write(*,*) "============================================================"
    call reset_state()
    do i = 1, nx
        do j = 4, 8
            f_ice(i, j) = 1.0_wp
        end do
        do j = 4, 7
            ssa_mask_acy(i, j) = 1
        end do
        ssa_mask_acy(i, 3) = 3                  ! ice-free BELOW, ice ABOVE
        ssa_mask_acy(i, 8) = 3                  ! ice BELOW, ice-free ABOVE
        do j = 4, 8
            ssa_mask_acx(i, j) = 1
        end do
        taul_int_acy(i, 3) = 1.0e6_wp
        taul_int_acy(i, 8) = 1.0e6_wp
    end do

    call linear_solver_matrix_ssa_ac_csr_2D_energy(lgs, ux, uy, beta_acx, beta_acy, N_aa, &
                ssa_mask_acx, ssa_mask_acy, mask_frnt, H_ice, f_ice, taud_acx, taud_acy, &
                taul_int_acx, taul_int_acy, dx, dy, beta_min, "periodic", "none")
    call check_tb_reflection(lgs, total_fail)

    ! ---------- Case C: square ice patch, both reflections at once ----------
    ! Ice in cells (i, j) with i in [4,8] AND j in [4,8]. All four fronts
    ! (left/right/top/bottom) carry taul_int. Both reflections must hold.
    write(*,*)
    write(*,*) "============================================================"
    write(*,*) "Case C: square ice patch, both L<->R and T<->B"
    write(*,*) "============================================================"
    call reset_state()
    do j = 4, 8
        do i = 4, 8
            f_ice(i, j) = 1.0_wp
        end do
    end do
    do j = 4, 8
        ! ux fronts at i = 3 and i = 8, ux interior at i = 4..7
        ssa_mask_acx(3, j) = 3
        ssa_mask_acx(8, j) = 3
        do i = 4, 7
            ssa_mask_acx(i, j) = 1
        end do
        taul_int_acx(3, j) = 1.0e6_wp
        taul_int_acx(8, j) = 1.0e6_wp
    end do
    do i = 4, 8
        ssa_mask_acy(i, 3) = 3
        ssa_mask_acy(i, 8) = 3
        do j = 4, 7
            ssa_mask_acy(i, j) = 1
        end do
        taul_int_acy(i, 3) = 1.0e6_wp
        taul_int_acy(i, 8) = 1.0e6_wp
    end do

    call linear_solver_matrix_ssa_ac_csr_2D_energy(lgs, ux, uy, beta_acx, beta_acy, N_aa, &
                ssa_mask_acx, ssa_mask_acy, mask_frnt, H_ice, f_ice, taud_acx, taud_acy, &
                taul_int_acx, taul_int_acy, dx, dy, beta_min, "infinite", "none")
    call check_lr_reflection(lgs, total_fail)
    call check_tb_reflection(lgs, total_fail)

    write(*,*)
    write(*,*) "============================================================"
    if (total_fail == 0) then
        write(*,*) " ALL CASES PASS: K and b are reflection-symmetric"
    else
        write(*,'(a,i0,a)') "  FAIL across ", total_fail, " case(s)"
        stop 1
    end if

contains

    subroutine reset_state()
        ! Re-initialise the per-case scratch arrays so cases don't bleed into each other.
        f_ice         = 0.0_wp
        ssa_mask_acx  = 0
        ssa_mask_acy  = 0
        taul_int_acx  = 0.0_wp
        taul_int_acy  = 0.0_wp
    end subroutine reset_state

    subroutine check_lr_reflection(lgs, total_fail)
        ! For every row n (== ux at (i,j) when n odd, uy at (i,j) when n even),
        ! find its L<->R mirror row n_m and verify K and b are reflection-equivariant.
        type(linear_solver_class), intent(in) :: lgs
        integer, intent(inout) :: total_fail
        integer :: sign_lookup_var(2)            ! sign under L<->R for var=1 (ux), var=2 (uy)

        sign_lookup_var(1) = -1                  ! ux flips sign under x-reflection
        sign_lookup_var(2) = +1                  ! uy unchanged
        call check_reflection_generic(lgs, total_fail, axis="x", sign_var=sign_lookup_var)
    end subroutine check_lr_reflection

    subroutine check_tb_reflection(lgs, total_fail)
        type(linear_solver_class), intent(in) :: lgs
        integer, intent(inout) :: total_fail
        integer :: sign_lookup_var(2)

        sign_lookup_var(1) = +1                  ! ux unchanged under y-reflection
        sign_lookup_var(2) = -1                  ! uy flips sign
        call check_reflection_generic(lgs, total_fail, axis="y", sign_var=sign_lookup_var)
    end subroutine check_tb_reflection

    subroutine check_reflection_generic(lgs, total_fail, axis, sign_var)
        ! Walk every row, compute its mirror under reflection about the given axis,
        ! then compare K(n, c) <-> sign_n*sign_c*K(n_m, c_m) for every entry and
        ! b(n) <-> sign_n*b(n_m). Reports worst violations and counts failures.
        type(linear_solver_class), intent(in) :: lgs
        integer, intent(inout) :: total_fail
        character(len=*), intent(in) :: axis     ! "x" for L<->R, "y" for T<->B
        integer, intent(in) :: sign_var(2)

        integer  :: n, n_m, c, c_m, idx, idx_m
        integer  :: i_n, j_n, i_m, j_m, i_c, j_c, i_cm, j_cm
        integer  :: var_n, var_c, sign_n, sign_c
        real(wp) :: bmax, kmax
        real(wp) :: db, dK
        real(wp) :: worst_b, worst_K
        integer  :: worst_b_n, worst_K_n
        integer  :: n_b_fail, n_K_fail
        real(wp) :: val_m, val_n
        logical  :: found

        bmax = max(maxval(abs(real(lgs%b_value, wp))), 1.0_wp)
        kmax = max(maxval(abs(real(lgs%a_value, wp))), 1.0_wp)

        worst_b = 0.0_wp;  worst_b_n = 0
        worst_K = 0.0_wp;  worst_K_n = 0
        n_b_fail = 0
        n_K_fail = 0

        do n = 1, lgs%nmax
            call decode_row(lgs, n, i_n, j_n, var_n)
            call mirror_index(i_n, j_n, var_n, axis, i_m, j_m)
            if (i_m < 1 .or. i_m > nx .or. j_m < 1 .or. j_m > ny) cycle
            if (var_n == 1) then
                n_m = 2*lgs%ij2n(i_m, j_m) - 1
            else
                n_m = 2*lgs%ij2n(i_m, j_m)
            end if
            if (n_m == n) cycle                  ! self-mirror row, nothing to check
            sign_n = sign_var(var_n)

            ! --- b: b(n) = sign_n * b(n_m) ---
            db = abs(real(lgs%b_value(n), wp) - real(sign_n, wp)*real(lgs%b_value(n_m), wp)) / bmax
            if (db > tol_b) then
                n_b_fail = n_b_fail + 1
                if (n_b_fail <= 5) then
                    write(*,'(a,a1,a,i0,a,2es14.6,a,es14.6)') "  b ", axis, &
                        " row=", n, " (n, n_m): ", real(lgs%b_value(n), wp), real(lgs%b_value(n_m), wp), &
                        "  rel-err=", db
                end if
            end if
            if (db > worst_b) then
                worst_b = db
                worst_b_n = n
            end if

            ! --- K: K(n, c) = sign_n * sign_c * K(n_m, c_m) ---
            do idx = lgs%a_ptr(n), lgs%a_ptr(n+1)-1
                c = lgs%a_index(idx)
                call decode_row(lgs, c, i_c, j_c, var_c)
                call mirror_index(i_c, j_c, var_c, axis, i_cm, j_cm)
                if (i_cm < 1 .or. i_cm > nx .or. j_cm < 1 .or. j_cm > ny) cycle
                if (var_c == 1) then
                    c_m = 2*lgs%ij2n(i_cm, j_cm) - 1
                else
                    c_m = 2*lgs%ij2n(i_cm, j_cm)
                end if
                sign_c = sign_var(var_c)

                val_n = real(lgs%a_value(idx), wp)
                val_m = 0.0_wp
                found = .FALSE.
                do idx_m = lgs%a_ptr(n_m), lgs%a_ptr(n_m+1)-1
                    if (lgs%a_index(idx_m) == c_m) then
                        val_m = real(lgs%a_value(idx_m), wp)
                        found = .TRUE.
                        exit
                    end if
                end do
                ! Missing slot is treated as 0 — what matters is the value.
                dK = abs(val_n - real(sign_n*sign_c, wp)*val_m) / kmax
                if (dK > tol_K) then
                    n_K_fail = n_K_fail + 1
                    if (n_K_fail <= 5) then
                        write(*,'(a,a1,a,i0,a,i0,a,2es14.6,a,es14.6)') "  K ", axis, &
                            " (n,c)=(", n, ",", c, ") (val, mirror): ", val_n, val_m, &
                            "  rel-err=", dK
                    end if
                end if
                if (dK > worst_K) then
                    worst_K = dK
                    worst_K_n = n
                end if
            end do
        end do

        write(*,'(a,a1)')      "  reflection axis      : ", axis
        write(*,'(a,1x,i0)')   "  b violations         :", n_b_fail
        write(*,'(a,1x,i0)')   "  K violations         :", n_K_fail
        write(*,'(a,es12.4,a,i0)') "  max b rel-err        :", worst_b, "   @row=", worst_b_n
        write(*,'(a,es12.4,a,i0)') "  max K rel-err        :", worst_K, "   @row=", worst_K_n
        if (n_b_fail /= 0 .or. n_K_fail /= 0) then
            write(*,*) "  -> FAIL"
            total_fail = total_fail + 1
        else
            write(*,*) "  -> pass"
        end if
    end subroutine check_reflection_generic

    subroutine decode_row(lgs, n, i, j, var)
        ! n -> (i, j, var) where var=1 means ux (n odd), var=2 means uy (n even).
        type(linear_solver_class), intent(in) :: lgs
        integer, intent(in)  :: n
        integer, intent(out) :: i, j, var
        i   = lgs%n2i((n+1)/2)
        j   = lgs%n2j((n+1)/2)
        var = 2 - mod(n, 2)
    end subroutine decode_row

    subroutine mirror_index(i, j, var, axis, i_m, j_m)
        ! C-grid reflection. ux lives at the right face of cell i, uy at top face of cell j.
        ! See header comment for the index map.
        integer, intent(in)  :: i, j, var
        character(len=*), intent(in) :: axis
        integer, intent(out) :: i_m, j_m

        select case (trim(axis))
            case ("x")
                if (var == 1) then
                    i_m = nx - i               ! ux right-face mirror
                else
                    i_m = nx + 1 - i           ! uy aa-column mirror
                end if
                j_m = j
            case ("y")
                if (var == 2) then
                    j_m = ny - j               ! uy top-face mirror
                else
                    j_m = ny + 1 - j           ! ux aa-row mirror
                end if
                i_m = i
            case default
                ! Unknown axis: return invalid index so caller's bounds check trips.
                i_m = -1
                j_m = -1
        end select
    end subroutine mirror_index

end program test_ssa_energy_lr_sym
