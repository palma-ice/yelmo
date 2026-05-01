module fast_hydrology_k24
    ! K24 effective-pressure / water-flux model.
    ! Diagnostic only: reads ice geometry, bed, melt, sliding speed and
    ! Glen rate factor; returns water flux components (q_x, q_y), effective
    ! pressure (N) and water pressure (p_w = Po - N). Does not evolve H_w.

    use nml

    implicit none

    integer, parameter :: dp = kind(1.d0)

    ! ---------- Substrate-type enum (par%k24%substrate_type) ----------
    integer, parameter, public :: K24_SUBSTRATE_HARD  = 0
    integer, parameter, public :: K24_SUBSTRATE_SOFT  = 1
    integer, parameter, public :: K24_SUBSTRATE_MIXED = 2

    ! ---------- Flux-solver enum (par%k24%flux_solver) ----------
    integer, parameter, public :: K24_FLUX_RECURSIVE = 0
    integer, parameter, public :: K24_FLUX_TOPOSORT  = 1

    ! ---------- Loaded "constants" (runtime, namelist-overridable) ----------
    type k24_param_class
        integer  :: substrate_type
        integer  :: flux_solver
        real(dp) :: water_density
        real(dp) :: ice_density
        real(dp) :: gravity
        real(dp) :: manning_exponent
        real(dp) :: latent_heat_water
        real(dp) :: bed_thickness
        real(dp) :: manning_coefficient_exponent
        real(dp) :: bed_friction_exponent
        real(dp) :: friction_factor
        real(dp) :: till_factor
        real(dp) :: critical_discharge
        real(dp) :: initial_cavity_height
        real(dp) :: coupling_length
        real(dp) :: long_coupling_water
        real(dp) :: seconds_per_year
        real(dp) :: min_pressure_fraction
    end type

    private
    public :: k24_param_class
    public :: k24_par_load
    public :: initialize_kappa
    public :: calc_k24
    public :: update_psi_out                 ! dispatcher
    public :: update_psi_out_recursive       ! exposed for testing/comparison
    public :: update_psi_out_toposort        ! exposed for testing/comparison

contains

    ! ============================================================
    ! Namelist load
    ! ============================================================
    subroutine k24_par_load(par, filename, init)

        implicit none

        type(k24_param_class), intent(INOUT) :: par
        character(len=*),      intent(IN)    :: filename
        logical, optional,     intent(IN)    :: init

        logical :: init_pars

        init_pars = .FALSE.
        if (present(init)) init_pars = init

        ! Defaults match the values originally hardcoded in
        ! hydrology_constants (see git history of fast_hydrology.f90).
        par%substrate_type                = K24_SUBSTRATE_HARD
        par%flux_solver                   = K24_FLUX_RECURSIVE
        par%water_density                 = 1000.0_dp
        par%ice_density                   =  917.0_dp
        par%gravity                       =    9.81_dp
        par%manning_exponent              =    3.0_dp
        par%latent_heat_water             =    3.35e5_dp
        par%bed_thickness                 =    0.1_dp
        par%manning_coefficient_exponent  =    1.25_dp
        par%bed_friction_exponent         =    1.5_dp
        par%friction_factor               =    0.1_dp
        par%till_factor                   =    1.1_dp
        par%critical_discharge            =    1.0_dp
        par%initial_cavity_height         =    0.1_dp
        par%coupling_length               =    1.0e4_dp
        par%long_coupling_water           =    5.0_dp
        par%seconds_per_year              =    3.154e7_dp
        par%min_pressure_fraction         =    0.02_dp

        call nml_read(filename,"fast_hydrology_k24","substrate_type",                par%substrate_type,                init=init_pars)
        call nml_read(filename,"fast_hydrology_k24","flux_solver",                   par%flux_solver,                   init=init_pars)
        call nml_read(filename,"fast_hydrology_k24","water_density",                 par%water_density,                 init=init_pars)
        call nml_read(filename,"fast_hydrology_k24","ice_density",                   par%ice_density,                   init=init_pars)
        call nml_read(filename,"fast_hydrology_k24","gravity",                       par%gravity,                       init=init_pars)
        call nml_read(filename,"fast_hydrology_k24","manning_exponent",              par%manning_exponent,              init=init_pars)
        call nml_read(filename,"fast_hydrology_k24","latent_heat_water",             par%latent_heat_water,             init=init_pars)
        call nml_read(filename,"fast_hydrology_k24","bed_thickness",                 par%bed_thickness,                 init=init_pars)
        call nml_read(filename,"fast_hydrology_k24","manning_coefficient_exponent",  par%manning_coefficient_exponent,  init=init_pars)
        call nml_read(filename,"fast_hydrology_k24","bed_friction_exponent",         par%bed_friction_exponent,         init=init_pars)
        call nml_read(filename,"fast_hydrology_k24","friction_factor",               par%friction_factor,               init=init_pars)
        call nml_read(filename,"fast_hydrology_k24","till_factor",                   par%till_factor,                   init=init_pars)
        call nml_read(filename,"fast_hydrology_k24","critical_discharge",            par%critical_discharge,            init=init_pars)
        call nml_read(filename,"fast_hydrology_k24","initial_cavity_height",         par%initial_cavity_height,         init=init_pars)
        call nml_read(filename,"fast_hydrology_k24","coupling_length",               par%coupling_length,               init=init_pars)
        call nml_read(filename,"fast_hydrology_k24","long_coupling_water",           par%long_coupling_water,           init=init_pars)
        call nml_read(filename,"fast_hydrology_k24","seconds_per_year",              par%seconds_per_year,              init=init_pars)
        call nml_read(filename,"fast_hydrology_k24","min_pressure_fraction",         par%min_pressure_fraction,         init=init_pars)

        return

    end subroutine k24_par_load

    ! ============================================================
    ! Substrate indicator (one-time, called from hydro_init_state)
    ! ============================================================
    subroutine initialize_kappa(kappa, b, substrate_type)

        implicit none

        real(dp),   intent(OUT) :: kappa(:,:)
        real(dp),   intent(IN)  :: b(:,:)
        integer,    intent(IN)  :: substrate_type

        integer :: i, j, nx, ny

        nx = size(kappa,1)
        ny = size(kappa,2)

        select case (substrate_type)
            case (K24_SUBSTRATE_HARD)
                kappa = 0.0_dp
            case (K24_SUBSTRATE_SOFT)
                kappa = 1.0_dp
            case (K24_SUBSTRATE_MIXED)
                do j = 1, ny
                    do i = 1, nx
                        if (b(i,j) < -1000.0_dp) then
                            kappa(i,j) = 1.0_dp
                        else
                            kappa(i,j) = 0.0_dp
                        end if
                    end do
                end do
            case default
                write(*,*) "initialize_kappa:: error: substrate_type must be one of [0,1,2]."
                write(*,*) "substrate_type = ", substrate_type
                stop
        end select

        return

    end subroutine initialize_kappa

    ! ============================================================
    ! Top-level K24 driver
    ! ============================================================
    subroutine calc_k24(q_x, q_y, N, p_w, &
                       H_ice, z_bed, mask, bmb_w, uxy_b, A_glen, kappa, &
                       dx, par)

        implicit none

        real(dp), intent(OUT) :: q_x(:,:), q_y(:,:)     ! water flux components
        real(dp), intent(OUT) :: N(:,:)                 ! effective pressure
        real(dp), intent(OUT) :: p_w(:,:)               ! water pressure (Po - N)
        real(dp), intent(IN)  :: H_ice(:,:)
        real(dp), intent(IN)  :: z_bed(:,:)
        real(dp), intent(IN)  :: mask(:,:)
        real(dp), intent(IN)  :: bmb_w(:,:)             ! basal melt (water-equiv. m/a)
        real(dp), intent(IN)  :: uxy_b(:,:)             ! basal sliding speed magnitude
        real(dp), intent(IN)  :: A_glen(:,:)
        real(dp), intent(IN)  :: kappa(:,:)
        real(dp), intent(IN)  :: dx
        type(k24_param_class), intent(IN) :: par

        ! Locals (K24 internals — all scoped to this call)
        integer :: nx, ny, i, j
        real(dp) :: K, dx_sq
        real(dp), allocatable :: phi_0(:,:)
        real(dp), allocatable :: minus_grad_phi_x(:,:), minus_grad_phi_y(:,:)
        real(dp), allocatable :: abs_grad_phi(:,:)
        real(dp), allocatable :: psi_out(:,:), q_flux(:,:), Q_disc(:,:)
        real(dp), allocatable :: S_inf(:,:), H_hard(:,:), H_soft(:,:), H_conduit(:,:)
        real(dp), allocatable :: N_inf(:,:), Po(:,:), corfac(:,:)
        real(dp), allocatable :: m_dot(:,:), h_for_pot(:,:)

        nx = size(H_ice,1)
        ny = size(H_ice,2)
        dx_sq = dx*dx

        allocate(phi_0(nx,ny))
        allocate(minus_grad_phi_x(nx,ny), minus_grad_phi_y(nx,ny), abs_grad_phi(nx,ny))
        allocate(psi_out(nx,ny), q_flux(nx,ny), Q_disc(nx,ny))
        allocate(S_inf(nx,ny), H_hard(nx,ny), H_soft(nx,ny), H_conduit(nx,ny))
        allocate(N_inf(nx,ny), Po(nx,ny), corfac(nx,ny))
        allocate(m_dot(nx,ny), h_for_pot(nx,ny))

        ! Manning-Strickler coefficient
        K = (2.0_dp / 3.14159265358979323846_dp)**0.25_dp * &
            sqrt((3.14159265358979323846_dp + 2.0_dp) / (par%water_density * par%friction_factor))

        ! Caller passes melt as m/a water-equivalent; K24 uses m_dot/rho_w
        ! conventionally in the same unit space — match original program.
        m_dot = bmb_w

        ! ========== Distributed water flux ==========
        h_for_pot = H_ice
        call update_potential(phi_0, h_for_pot, z_bed, par)
        call potential_filling(phi_0, 10)

        ! Recover h consistent with filled potential (matches original program)
        do j = 1, ny
            do i = 1, nx
                h_for_pot(i,j) = (phi_0(i,j) - par%water_density * par%gravity * z_bed(i,j)) &
                                / (par%ice_density * par%gravity)
            end do
        end do

        call update_potential_gradients_smoothed(minus_grad_phi_x, minus_grad_phi_y, abs_grad_phi, &
                                                 phi_0, h_for_pot, dx, mask, par)

        call update_psi_out(psi_out, m_dot, dx_sq, &
                            minus_grad_phi_x, minus_grad_phi_y, abs_grad_phi, mask, &
                            par%flux_solver)

        ! Coordinate-system correction
        do j = 1, ny
            do i = 1, nx
                corfac(i,j) = sqrt(minus_grad_phi_x(i,j)**2 + minus_grad_phi_y(i,j)**2)
                if (corfac(i,j) > 1.0e-12_dp) then
                    corfac(i,j) = abs_grad_phi(i,j) / corfac(i,j)
                else
                    corfac(i,j) = 1.0_dp
                end if
            end do
        end do

        do j = 1, ny
            do i = 1, nx
                q_flux(i,j) = psi_out(i,j) / (corfac(i,j) * dx)
                q_flux(i,j) = min(max(q_flux(i,j), 0.0_dp), 1.0e5_dp)
            end do
        end do

        ! ========== Effective pressure ==========
        ! Recompute potential and unsmoothed gradients (matches original program)
        do j = 1, ny
            do i = 1, nx
                phi_0(i,j) = par%ice_density * par%gravity * h_for_pot(i,j) &
                           + par%water_density * par%gravity * z_bed(i,j)
            end do
        end do

        do j = 1, ny
            do i = 2, nx-1
                minus_grad_phi_x(i,j) = -(phi_0(i+1,j) - phi_0(i-1,j)) / (2.0_dp * dx)
            end do
        end do
        minus_grad_phi_x(1,:)  = minus_grad_phi_x(2,:)
        minus_grad_phi_x(nx,:) = minus_grad_phi_x(nx-1,:)

        do j = 2, ny-1
            do i = 1, nx
                minus_grad_phi_y(i,j) = -(phi_0(i,j+1) - phi_0(i,j-1)) / (2.0_dp * dx)
            end do
        end do
        minus_grad_phi_y(:,1)  = minus_grad_phi_y(:,2)
        minus_grad_phi_y(:,ny) = minus_grad_phi_y(:,ny-1)

        do j = 1, ny
            do i = 1, nx
                abs_grad_phi(i,j) = sqrt(minus_grad_phi_x(i,j)**2 + minus_grad_phi_y(i,j)**2)
            end do
        end do

        do j = 1, ny
            do i = 1, nx
                Q_disc(i,j) = q_flux(i,j) * par%coupling_length
            end do
        end do

        call update_S_inf(S_inf, K, par%manning_coefficient_exponent, par%bed_friction_exponent, &
                          abs_grad_phi, Q_disc)
        call update_H_conduit(H_conduit, H_hard, H_soft, S_inf, &
                              par%initial_cavity_height, par%till_factor, &
                              Q_disc, par%critical_discharge, kappa)
        call update_Po(Po, par%ice_density, par%gravity, h_for_pot)
        call update_N_inf(N_inf, H_conduit, S_inf, par%ice_density, par%latent_heat_water, &
                          uxy_b, par%bed_thickness, Q_disc, abs_grad_phi, par%manning_exponent, &
                          A_glen, Po, par%min_pressure_fraction)
        call update_N(N, N_inf, phi_0)

        ! Output components: signed flux from gradient direction.
        ! q_flux is the magnitude; allocate to (q_x, q_y) by the unit gradient
        ! components. Where abs_grad_phi vanishes, distribute zero.
        do j = 1, ny
            do i = 1, nx
                if (abs_grad_phi(i,j) > 1.0e-12_dp) then
                    q_x(i,j) = q_flux(i,j) * minus_grad_phi_x(i,j) / abs_grad_phi(i,j)
                    q_y(i,j) = q_flux(i,j) * minus_grad_phi_y(i,j) / abs_grad_phi(i,j)
                else
                    q_x(i,j) = 0.0_dp
                    q_y(i,j) = 0.0_dp
                end if
                p_w(i,j) = Po(i,j) - N(i,j)
            end do
        end do

        deallocate(phi_0, minus_grad_phi_x, minus_grad_phi_y, abs_grad_phi)
        deallocate(psi_out, q_flux, Q_disc)
        deallocate(S_inf, H_hard, H_soft, H_conduit)
        deallocate(N_inf, Po, corfac, m_dot, h_for_pot)

        return

    end subroutine calc_k24

    ! ============================================================
    ! Hydraulic potential
    ! ============================================================
    subroutine update_potential(phi_0, h, b, par)
        implicit none
        real(dp), intent(OUT) :: phi_0(:,:)
        real(dp), intent(IN)  :: h(:,:), b(:,:)
        type(k24_param_class), intent(IN) :: par
        integer :: i, j, nx, ny

        nx = size(phi_0,1); ny = size(phi_0,2)
        do j = 1, ny
            do i = 1, nx
                phi_0(i,j) = par%ice_density * par%gravity * h(i,j) &
                           + par%water_density * par%gravity * b(i,j)
            end do
        end do
    end subroutine update_potential

    ! ============================================================
    ! Iterative hollow-filling for spurious sinks
    ! ============================================================
    subroutine potential_filling(phi_0, iterations)
        implicit none
        real(dp), intent(INOUT) :: phi_0(:,:)
        integer,  intent(IN)    :: iterations

        real(dp), allocatable :: pot_next(:,:)
        integer :: iter, i, j, nx, ny
        real(dp) :: p, p_mean
        logical :: all_greater

        nx = size(phi_0,1); ny = size(phi_0,2)
        allocate(pot_next(nx,ny))

        do iter = 1, iterations
            pot_next = phi_0
            do j = 2, ny-1
                do i = 2, nx-1
                    p = phi_0(i,j)
                    all_greater = phi_0(i+1,j) > p .and. phi_0(i-1,j) > p .and. &
                                  phi_0(i,j+1) > p .and. phi_0(i,j-1) > p
                    if (all_greater) then
                        p_mean = (phi_0(i+1,j) + phi_0(i-1,j) + &
                                  phi_0(i,j+1) + phi_0(i,j-1)) / 4.0_dp
                        pot_next(i,j) = p_mean
                    end if
                end do
            end do
            phi_0 = pot_next
        end do

        deallocate(pot_next)
    end subroutine potential_filling

    ! ============================================================
    ! FFT-smoothed gradient field
    ! ============================================================
    subroutine update_potential_gradients_smoothed(minus_grad_phi_x, minus_grad_phi_y, &
                                                   abs_grad_phi, phi_0, h, dx, mask, par)
        implicit none
        real(dp), intent(OUT) :: minus_grad_phi_x(:,:), minus_grad_phi_y(:,:), abs_grad_phi(:,:)
        real(dp), intent(IN)  :: phi_0(:,:), h(:,:), mask(:,:), dx
        type(k24_param_class), intent(IN) :: par

        real(dp), allocatable :: kernel(:,:), tx(:,:), ty(:,:)
        real(dp) :: h_avg, scale, width, dist, kernel_sum
        integer  :: kernel_size, frb, i, j, ni, nj, nx, ny

        nx = size(phi_0,1); ny = size(phi_0,2)

        ! Unsmoothed gradients
        do j = 1, ny
            do i = 2, nx-1
                minus_grad_phi_x(i,j) = -(phi_0(i+1,j) - phi_0(i-1,j)) / (2.0_dp * dx)
            end do
        end do
        minus_grad_phi_x(1,:)  = minus_grad_phi_x(2,:)
        minus_grad_phi_x(nx,:) = minus_grad_phi_x(nx-1,:)

        do j = 2, ny-1
            do i = 1, nx
                minus_grad_phi_y(i,j) = -(phi_0(i,j+1) - phi_0(i,j-1)) / (2.0_dp * dx)
            end do
        end do
        minus_grad_phi_y(:,1)  = minus_grad_phi_y(:,2)
        minus_grad_phi_y(:,ny) = minus_grad_phi_y(:,ny-1)

        ! Kernel scale from mean active ice thickness
        h_avg = 0.0_dp
        do j = 1, ny
            do i = 1, nx
                if (mask(i,j) == 1.0_dp) h_avg = h_avg + h(i,j)
            end do
        end do
        h_avg = h_avg / max(1, count(mask == 1.0_dp))
        h_avg = max(h_avg, 10.0_dp)

        scale = h_avg * par%long_coupling_water * 2.0_dp
        width = 2.0_dp * scale
        if (width <= dx) then
            scale = dx / 2.0_dp + 1.0_dp
            width = 2.0_dp * scale
        end if

        kernel_size = 2 * nint(width / dx - 0.5_dp) + 1
        frb = (kernel_size - 1) / 2

        ! Cone kernel
        allocate(kernel(kernel_size, kernel_size))
        kernel = 0.0_dp
        do nj = 1, kernel_size
            do ni = 1, kernel_size
                dist = sqrt( (dx * real(ni - frb - 1, dp))**2 + &
                             (dx * real(nj - frb - 1, dp))**2 ) / scale
                kernel(ni,nj) = max(0.0_dp, 1.0_dp - dist / 2.0_dp)
            end do
        end do
        kernel_sum = sum(kernel)
        if (kernel_sum > 0.0_dp) kernel = kernel / kernel_sum

        ! Save unsmoothed, then smooth
        allocate(tx(nx,ny), ty(nx,ny))
        tx = minus_grad_phi_x
        ty = minus_grad_phi_y

        call imfilter_reflect_fftw(tx, nx, ny, kernel, kernel_size, frb, minus_grad_phi_x)
        call imfilter_reflect_fftw(ty, nx, ny, kernel, kernel_size, frb, minus_grad_phi_y)

        do j = 1, ny
            do i = 1, nx
                abs_grad_phi(i,j) = abs(minus_grad_phi_x(i,j)) + abs(minus_grad_phi_y(i,j))
            end do
        end do

        deallocate(kernel, tx, ty)
    end subroutine update_potential_gradients_smoothed

    ! ============================================================
    ! 2-D convolution via FFTW3 with reflect boundary conditions.
    ! Equivalent to ImageFiltering.imfilter(input, centered(kernel), "reflect").
    ! Link with -lfftw3.
    ! ============================================================
    subroutine imfilter_reflect_fftw(input, nx, ny, kernel, kernel_size, frb, output)
        use, intrinsic :: iso_c_binding
        implicit none
        include 'fftw3.f03'

        integer,  intent(IN)  :: nx, ny, kernel_size, frb
        real(dp), intent(IN)  :: input(nx, ny)
        real(dp), intent(IN)  :: kernel(kernel_size, kernel_size)
        real(dp), intent(OUT) :: output(nx, ny)

        integer :: Npx, Npy, Mfft, Nfft, i, j, ii, jj, ni, nj
        real(C_DOUBLE),            allocatable :: work_a(:,:), work_b(:,:), work_out(:,:)
        complex(C_DOUBLE_COMPLEX), allocatable :: A_hat(:,:), B_hat(:,:)
        type(C_PTR) :: plan_fwd_a, plan_fwd_b, plan_bwd

        Npx = nx + 2*frb
        Npy = ny + 2*frb

        Mfft = 1
        do while (Mfft < Npx + kernel_size - 1); Mfft = Mfft * 2; end do
        Nfft = 1
        do while (Nfft < Npy + kernel_size - 1); Nfft = Nfft * 2; end do

        allocate(work_a(Mfft, Nfft), work_b(Mfft, Nfft), work_out(Mfft, Nfft))
        allocate(A_hat(Mfft/2+1, Nfft), B_hat(Mfft/2+1, Nfft))
        work_a = 0.0_dp
        work_b = 0.0_dp

        do j = 1, Npy
            do i = 1, Npx
                ii = i - frb;  jj = j - frb
                if (ii < 1)  ii = 2 - ii
                if (ii > nx) ii = 2*nx - ii
                if (jj < 1)  jj = 2 - jj
                if (jj > ny) jj = 2*ny - jj
                work_a(i, j) = input(ii, jj)
            end do
        end do

        do nj = 1, kernel_size
            do ni = 1, kernel_size
                ii = mod(ni - frb - 1 + Mfft, Mfft) + 1
                jj = mod(nj - frb - 1 + Nfft, Nfft) + 1
                work_b(ii, jj) = kernel(ni, nj)
            end do
        end do

        plan_fwd_a = fftw_plan_dft_r2c_2d(Nfft, Mfft, work_a, A_hat, FFTW_ESTIMATE)
        call fftw_execute_dft_r2c(plan_fwd_a, work_a, A_hat)
        call fftw_destroy_plan(plan_fwd_a)

        plan_fwd_b = fftw_plan_dft_r2c_2d(Nfft, Mfft, work_b, B_hat, FFTW_ESTIMATE)
        call fftw_execute_dft_r2c(plan_fwd_b, work_b, B_hat)
        call fftw_destroy_plan(plan_fwd_b)

        A_hat = A_hat * B_hat

        plan_bwd = fftw_plan_dft_c2r_2d(Nfft, Mfft, A_hat, work_out, FFTW_ESTIMATE)
        call fftw_execute_dft_c2r(plan_bwd, A_hat, work_out)
        call fftw_destroy_plan(plan_bwd)

        work_out = work_out / real(Mfft * Nfft, dp)

        do j = 1, ny
            do i = 1, nx
                output(i, j) = work_out(i + 2*frb, j + 2*frb)
            end do
        end do

        deallocate(work_a, work_b, work_out, A_hat, B_hat)
    end subroutine imfilter_reflect_fftw

    ! ============================================================
    ! Water-flux accumulation: dispatcher + two implementations.
    !
    ! Two solvers compute the same conceptual quantity (psi_out — accumulated
    ! upstream water-potential outflow per cell) by different algorithms:
    !
    !   * _recursive — depth-first traversal with memoization. The original
    !     formulation. Simple and direct, but Fortran recursion uses the call
    !     stack and can overflow on large active regions with long flow paths.
    !
    !   * _toposort  — topological ordering of cells by inflow direction,
    !     followed by a single linear pass in topological order. Equivalent in
    !     principle. CISM uses this approach, and FastHydrology.jl observed
    !     small numerical differences from the recursive form (presumably from
    !     accumulation order in floating-point sums). Stack-safe at any grid
    !     size.
    !
    ! Both expose the same signature and write to the same psi_out shape, so
    ! callers select via par%flux_solver. Recursive is the default to preserve
    ! current numerics; the toposort variant is provided for A/B comparison
    ! and for promotion to default once validated.
    ! ============================================================
    subroutine update_psi_out(psi_out, m_dot, dx_sq, &
                              minus_grad_phi_x, minus_grad_phi_y, abs_grad_phi, mask, &
                              flux_solver)
        implicit none
        real(dp), intent(INOUT) :: psi_out(:,:)
        real(dp), intent(IN)    :: m_dot(:,:), dx_sq
        real(dp), intent(IN)    :: minus_grad_phi_x(:,:), minus_grad_phi_y(:,:)
        real(dp), intent(IN)    :: abs_grad_phi(:,:), mask(:,:)
        integer,  intent(IN)    :: flux_solver

        select case (flux_solver)
            case (K24_FLUX_RECURSIVE)
                call update_psi_out_recursive(psi_out, m_dot, dx_sq, &
                                              minus_grad_phi_x, minus_grad_phi_y, &
                                              abs_grad_phi, mask)
            case (K24_FLUX_TOPOSORT)
                call update_psi_out_toposort(psi_out, m_dot, dx_sq, &
                                             minus_grad_phi_x, minus_grad_phi_y, &
                                             abs_grad_phi, mask)
            case default
                write(*,*) "update_psi_out:: error: flux_solver must be one of [0,1]."
                write(*,*) "flux_solver = ", flux_solver
                stop
        end select
    end subroutine update_psi_out

    ! --- Recursive (DFS + memoization) --------------------------
    subroutine update_psi_out_recursive(psi_out, m_dot, dx_sq, &
                                        minus_grad_phi_x, minus_grad_phi_y, &
                                        abs_grad_phi, mask)
        implicit none
        real(dp), intent(INOUT) :: psi_out(:,:)
        real(dp), intent(IN)    :: m_dot(:,:), dx_sq
        real(dp), intent(IN)    :: minus_grad_phi_x(:,:), minus_grad_phi_y(:,:)
        real(dp), intent(IN)    :: abs_grad_phi(:,:), mask(:,:)
        integer  :: i, j, nx, ny
        real(dp) :: dummy

        nx = size(psi_out,1); ny = size(psi_out,2)

        where (mask == 1.0_dp)
            psi_out = -1.0_dp
        end where

        do j = 1, ny
            do i = 1, nx
                if (mask(i,j) == 1.0_dp) then
                    dummy = accumulate_psi_out_recursive(psi_out, i, j, m_dot, dx_sq, &
                                                         minus_grad_phi_x, minus_grad_phi_y, &
                                                         abs_grad_phi, nx, ny)
                end if
            end do
        end do
    end subroutine update_psi_out_recursive

    recursive function accumulate_psi_out_recursive(psi_out, i, j, m_dot, dx_sq, &
                                                    minus_grad_phi_x, minus_grad_phi_y, &
                                                    abs_grad_phi, nx, ny) result(psi_value)
        implicit none
        integer,  intent(IN)    :: i, j, nx, ny
        real(dp), intent(INOUT) :: psi_out(:,:)
        real(dp), intent(IN)    :: m_dot(:,:), dx_sq
        real(dp), intent(IN)    :: minus_grad_phi_x(:,:), minus_grad_phi_y(:,:)
        real(dp), intent(IN)    :: abs_grad_phi(:,:)
        real(dp) :: psi_value

        integer  :: ni, nj, dir_idx
        real(dp) :: w, abs_grad
        integer, parameter :: directions(2,4) = reshape([ -1, 0, 1, 0, 0, -1, 0, 1 ], [2,4])

        if (psi_out(i,j) >= 0.0_dp) then
            psi_value = psi_out(i,j)
            return
        end if

        psi_out(i,j) = max(0.0_dp, m_dot(i,j) * dx_sq)

        do dir_idx = 1, 4
            ni = i + directions(1, dir_idx)
            nj = j + directions(2, dir_idx)
            if (ni < 1 .or. ni > nx .or. nj < 1 .or. nj > ny) cycle

            abs_grad = abs_grad_phi(ni, nj)
            if (abs_grad <= 1.0e-12_dp) cycle

            w = -(minus_grad_phi_x(ni,nj) * real(directions(1,dir_idx), dp) + &
                  minus_grad_phi_y(ni,nj) * real(directions(2,dir_idx), dp)) / abs_grad

            if (w > 0.0_dp) then
                psi_out(i,j) = psi_out(i,j) + &
                    accumulate_psi_out_recursive(psi_out, ni, nj, m_dot, dx_sq, &
                                                 minus_grad_phi_x, minus_grad_phi_y, &
                                                 abs_grad_phi, nx, ny) * w
            end if
        end do

        psi_value = psi_out(i,j)
    end function accumulate_psi_out_recursive

    ! --- Topological-sort (iterative, stack-safe) ----------------
    ! Equivalent computation to the recursive form, but ordered:
    !   1. For each active cell, count how many active *outgoing* neighbours
    !      contribute to it (= number of in-edges to this cell from
    !      neighbours where flow goes from neighbour into us).
    !   2. Seed a queue with cells that have zero in-degree (they receive no
    !      upstream contribution).
    !   3. Process the queue in order, accumulating psi_out from already-
    !      processed upstream cells, then decrementing the in-degree of each
    !      downstream neighbour and queuing any that hit zero.
    !
    ! Edge convention (unchanged from recursive form): cell n with negative
    ! direction vector d delivers flow into cell c = n + d when the weight
    ! w_{n→c} = -(grad_x(n)*d_x + grad_y(n)*d_y) / |grad(n)| is positive.
    subroutine update_psi_out_toposort(psi_out, m_dot, dx_sq, &
                                       minus_grad_phi_x, minus_grad_phi_y, &
                                       abs_grad_phi, mask)
        implicit none
        real(dp), intent(INOUT) :: psi_out(:,:)
        real(dp), intent(IN)    :: m_dot(:,:), dx_sq
        real(dp), intent(IN)    :: minus_grad_phi_x(:,:), minus_grad_phi_y(:,:)
        real(dp), intent(IN)    :: abs_grad_phi(:,:), mask(:,:)

        integer :: i, j, nx, ny, ni, nj, dir_idx, head, tail, qsize
        real(dp) :: w, abs_grad
        integer, allocatable :: in_degree(:,:)
        integer, allocatable :: queue_i(:), queue_j(:)
        integer, parameter :: directions(2,4) = reshape([ -1, 0, 1, 0, 0, -1, 0, 1 ], [2,4])

        nx = size(psi_out,1); ny = size(psi_out,2)

        allocate(in_degree(nx,ny))
        in_degree = 0

        ! Seed psi_out with the local melt source (same as recursive base case).
        ! Inactive cells stay untouched (matches recursive behaviour, which
        ! only writes to mask==1 cells).
        do j = 1, ny
            do i = 1, nx
                if (mask(i,j) == 1.0_dp) then
                    psi_out(i,j) = max(0.0_dp, m_dot(i,j) * dx_sq)
                end if
            end do
        end do

        ! Count in-edges: for cell c, count how many neighbours n send flow
        ! into c. Mirrors the recursive iteration over directions (the inner
        ! loop computes contributions from neighbours of the current cell).
        do j = 1, ny
            do i = 1, nx
                if (mask(i,j) /= 1.0_dp) cycle
                do dir_idx = 1, 4
                    ni = i + directions(1, dir_idx)
                    nj = j + directions(2, dir_idx)
                    if (ni < 1 .or. ni > nx .or. nj < 1 .or. nj > ny) cycle

                    abs_grad = abs_grad_phi(ni, nj)
                    if (abs_grad <= 1.0e-12_dp) cycle

                    w = -(minus_grad_phi_x(ni,nj) * real(directions(1,dir_idx), dp) + &
                          minus_grad_phi_y(ni,nj) * real(directions(2,dir_idx), dp)) / abs_grad

                    if (w > 0.0_dp) in_degree(i,j) = in_degree(i,j) + 1
                end do
            end do
        end do

        ! Queue is at most one slot per active cell.
        qsize = count(mask == 1.0_dp)
        allocate(queue_i(qsize), queue_j(qsize))
        head = 1
        tail = 0

        do j = 1, ny
            do i = 1, nx
                if (mask(i,j) == 1.0_dp .and. in_degree(i,j) == 0) then
                    tail = tail + 1
                    queue_i(tail) = i
                    queue_j(tail) = j
                end if
            end do
        end do

        ! Process in topological order. For each popped cell c, scan its
        ! neighbours: any neighbour d that c sends flow into has its
        ! in_degree decremented; once zero, d is queued. Crucially, before
        ! enqueuing d we *also* add c's contribution to d's psi_out — this
        ! is the flux-accumulation step that mirrors the recursive sum.
        do while (head <= tail)
            i = queue_i(head)
            j = queue_j(head)
            head = head + 1

            do dir_idx = 1, 4
                ni = i + directions(1, dir_idx)
                nj = j + directions(2, dir_idx)
                if (ni < 1 .or. ni > nx .or. nj < 1 .or. nj > ny) cycle
                if (mask(ni,nj) /= 1.0_dp) cycle

                ! Does cell (i,j) send flow into (ni,nj)? Same edge test as
                ! above but evaluated at the *source* (i,j) with direction d
                ! pointing toward (ni,nj). Equivalently, evaluate at (i,j)
                ! with direction (-d) flipped to match the convention.
                abs_grad = abs_grad_phi(i,j)
                if (abs_grad <= 1.0e-12_dp) cycle

                ! Weight from (ni,nj)'s perspective: c is at (ni,nj)+(-d), so
                ! the recursive form would query neighbour at offset (-d)
                ! from (ni,nj). Re-use the same formula with direction (-d).
                w = -(minus_grad_phi_x(i,j) * real(-directions(1,dir_idx), dp) + &
                      minus_grad_phi_y(i,j) * real(-directions(2,dir_idx), dp)) / abs_grad

                if (w > 0.0_dp) then
                    psi_out(ni,nj) = psi_out(ni,nj) + psi_out(i,j) * w
                    in_degree(ni,nj) = in_degree(ni,nj) - 1
                    if (in_degree(ni,nj) == 0) then
                        tail = tail + 1
                        queue_i(tail) = ni
                        queue_j(tail) = nj
                    end if
                end if
            end do
        end do

        ! Cells still with in_degree > 0 are part of a cycle (rare; the
        ! recursive form silently re-uses the memoized partial result on
        ! re-entry). Their psi_out keeps the seeded melt source — slight
        ! numerical divergence from the recursive form on cyclic flow paths.

        deallocate(in_degree, queue_i, queue_j)
    end subroutine update_psi_out_toposort

    ! ============================================================
    ! Effective-pressure helpers
    ! ============================================================
    subroutine update_Po(Po, rho_i, g, h)
        implicit none
        real(dp), intent(OUT) :: Po(:,:)
        real(dp), intent(IN)  :: rho_i, g, h(:,:)
        integer :: i, j, nx, ny
        nx = size(Po,1); ny = size(Po,2)
        do j = 1, ny
            do i = 1, nx
                Po(i,j) = max(rho_i * g * h(i,j), 1.0e5_dp)
            end do
        end do
    end subroutine update_Po

    subroutine update_S_inf(S_inf, K, alpha, beta, abs_grad_phi, Q)
        implicit none
        real(dp), intent(OUT) :: S_inf(:,:)
        real(dp), intent(IN)  :: K, alpha, beta
        real(dp), intent(IN)  :: abs_grad_phi(:,:), Q(:,:)
        integer :: i, j, nx, ny
        nx = size(S_inf,1); ny = size(S_inf,2)
        do j = 1, ny
            do i = 1, nx
                S_inf(i,j) = (K**(-1.0_dp / alpha)) * &
                             (abs_grad_phi(i,j)**((1.0_dp - beta) / alpha)) * &
                             (Q(i,j)**(1.0_dp / alpha))
                S_inf(i,j) = max(S_inf(i,j), 1.0e-12_dp)
            end do
        end do
    end subroutine update_S_inf

    subroutine update_H_conduit(H_conduit, H_hard, H_soft, S_inf, H_0, F_till, Q, Q_c, kappa)
        implicit none
        real(dp), intent(OUT) :: H_conduit(:,:), H_hard(:,:), H_soft(:,:)
        real(dp), intent(IN)  :: S_inf(:,:), H_0, F_till
        real(dp), intent(IN)  :: Q(:,:), Q_c, kappa(:,:)
        integer :: i, j, nx, ny
        nx = size(H_conduit,1); ny = size(H_conduit,2)
        do j = 1, ny
            do i = 1, nx
                H_hard(i,j) = sqrt(S_inf(i,j))
                H_soft(i,j) = H_0 + (sqrt(S_inf(i,j)) / F_till - H_0) * &
                              exp(-Q(i,j) / Q_c)
                H_conduit(i,j) = (1.0_dp - kappa(i,j)) * H_hard(i,j) + &
                                 kappa(i,j) * H_soft(i,j)
                H_conduit(i,j) = max(H_conduit(i,j), 1.0e-12_dp)
            end do
        end do
    end subroutine update_H_conduit

    subroutine update_N_inf(N_inf, H_conduit, S_inf, rho_i, L_w, abs_v_b, h_b, Q, &
                            abs_grad_phi, n, A, Po, min_pressure_fraction)
        implicit none
        real(dp), intent(OUT) :: N_inf(:,:)
        real(dp), intent(IN)  :: H_conduit(:,:), S_inf(:,:)
        real(dp), intent(IN)  :: rho_i, L_w, abs_v_b(:,:), h_b, Q(:,:)
        real(dp), intent(IN)  :: abs_grad_phi(:,:), n, A(:,:), Po(:,:)
        real(dp), intent(IN)  :: min_pressure_fraction
        real(dp) :: numerator, denominator, cavity_ratio
        integer :: i, j, nx, ny
        nx = size(N_inf,1); ny = size(N_inf,2)
        do j = 1, ny
            do i = 1, nx
                numerator   = rho_i * L_w * abs_v_b(i,j) * h_b + Q(i,j) * abs_grad_phi(i,j)
                denominator = 2.0_dp * (n**(-n)) * rho_i * L_w * A(i,j)
                cavity_ratio = H_conduit(i,j) / max(S_inf(i,j), 1.0e-12_dp)
                N_inf(i,j) = ((cavity_ratio**2 * numerator / denominator)**(1.0_dp / n))
                N_inf(i,j) = max(min(N_inf(i,j), Po(i,j)), &
                                 min_pressure_fraction * Po(i,j))
            end do
        end do
    end subroutine update_N_inf

    subroutine update_N(N, N_inf, phi_0)
        implicit none
        real(dp), intent(OUT) :: N(:,:)
        real(dp), intent(IN)  :: N_inf(:,:), phi_0(:,:)
        integer :: i, j, nx, ny
        nx = size(N,1); ny = size(N,2)
        do j = 1, ny
            do i = 1, nx
                N(i,j) = max(0.0_dp, erf(sqrt(3.14159265358979323846_dp) * phi_0(i,j) / &
                                         (2.0_dp * N_inf(i,j))) * N_inf(i,j))
            end do
        end do
    end subroutine update_N

end module fast_hydrology_k24
