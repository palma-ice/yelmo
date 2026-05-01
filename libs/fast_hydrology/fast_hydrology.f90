module fast_hydrology
    ! Top-level wrapper for distributed-water-flux / effective-pressure
    ! hydrology models. Mirrors the API style of basal_hydro_simple but
    ! dispatches to model-specific implementations (currently K24, with
    ! room for additional models to be added as new method values).
    !
    ! Contract for H_w:
    !   The persistent water-layer thickness (now%H_w) may be updated
    !   externally between calls (e.g. by the host ice sheet model).
    !   The models invoked here perform horizontal transport / diagnose
    !   flux and effective pressure on the supplied state; they do not
    !   take ownership of the H_w mass balance.

    use nml
    use fast_hydrology_k24

    implicit none

    integer, parameter :: dp = kind(1.d0)
    integer, parameter :: sp = kind(1.0)

    ! Working precision matches yelmo's wp convention
    integer, parameter :: wp = sp

    ! ---------- Method enum (par%method) ----------
    integer, parameter, public :: HYDRO_METHOD_NONE     = 0   ! no-op (state frozen)
    integer, parameter, public :: HYDRO_METHOD_RESERVED = 1   ! reserved for future port
    integer, parameter, public :: HYDRO_METHOD_K24      = 2

    ! ---------- Init-method enum (par%init_method) ----------
    integer, parameter, public :: HYDRO_INIT_EXTERNAL = 0     ! H_w filled externally

    type hydro_param_class
        integer  :: init_method
        integer  :: method
        real(wp) :: dx
        real(wp) :: dy
        type(k24_param_class) :: k24
    end type

    type hydro_state_class
        real(wp) :: time
        real(wp) :: dt
        real(wp), allocatable :: H_w(:,:)        ! water layer thickness; externally mutable
        real(wp), allocatable :: p_w(:,:)        ! water pressure (diagnostic, Po - N)
        real(wp), allocatable :: q_x(:,:)        ! water flux x-component (diagnostic)
        real(wp), allocatable :: q_y(:,:)        ! water flux y-component (diagnostic)
        real(wp), allocatable :: N(:,:)          ! effective pressure (diagnostic)
        real(wp), allocatable :: kappa(:,:)      ! substrate indicator (set in init_state)
    end type

    type hydro_class
        type(hydro_param_class) :: par
        type(hydro_state_class) :: now
    end type

    private
    public :: hydro_class
    public :: hydro_init
    public :: hydro_init_state
    public :: hydro_update
    public :: wp

contains

    subroutine hydro_init(hyd, filename, nx, ny)

        implicit none

        type(hydro_class), intent(INOUT) :: hyd
        character(len=*),  intent(IN)    :: filename
        integer,           intent(IN)    :: nx, ny

        call hydro_par_load(hyd%par, filename)

        if (hyd%par%dx /= hyd%par%dy) then
            write(*,*) "hydro_init:: error: K24 requires dx == dy."
            write(*,*) "dx = ", hyd%par%dx, "  dy = ", hyd%par%dy
            stop
        end if

        call hydro_allocate(hyd%now, nx, ny)

        hyd%now%time = 1.0e10_wp
        hyd%now%dt   = 0.0_wp

        return

    end subroutine hydro_init

    subroutine hydro_init_state(hyd, z_bed, time)
        ! Populate state items that are determined once per simulation
        ! (substrate indicator) and reset the time bookkeeping. Does NOT
        ! touch H_w — the caller (ice sheet model) is responsible for
        ! initializing it externally per the API contract.

        implicit none

        type(hydro_class), intent(INOUT) :: hyd
        real(wp),          intent(IN)    :: z_bed(:,:)
        real(wp),          intent(IN)    :: time

        real(dp), allocatable :: z_bed_dp(:,:), kappa_dp(:,:)
        integer :: nx, ny

        nx = size(z_bed,1)
        ny = size(z_bed,2)

        allocate(z_bed_dp(nx,ny), kappa_dp(nx,ny))
        z_bed_dp = real(z_bed, dp)

        call initialize_kappa(kappa_dp, z_bed_dp, hyd%par%k24%substrate_type)

        hyd%now%kappa = real(kappa_dp, wp)

        hyd%now%H_w  = 0.0_wp
        hyd%now%p_w  = 0.0_wp
        hyd%now%q_x  = 0.0_wp
        hyd%now%q_y  = 0.0_wp
        hyd%now%N    = 0.0_wp

        hyd%now%time = time
        hyd%now%dt   = 0.0_wp

        deallocate(z_bed_dp, kappa_dp)

        return

    end subroutine hydro_init_state

    subroutine hydro_update(hyd, H_ice, z_bed, mask, bmb_w, uxy_b, A_glen, time)

        implicit none

        type(hydro_class), intent(INOUT) :: hyd
        real(wp),          intent(IN)    :: H_ice(:,:)
        real(wp),          intent(IN)    :: z_bed(:,:)
        real(wp),          intent(IN)    :: mask(:,:)      ! active cells == 1.0
        real(wp),          intent(IN)    :: bmb_w(:,:)     ! basal melt (water-equiv. m/a)
        real(wp),          intent(IN)    :: uxy_b(:,:)     ! basal sliding speed magnitude
        real(wp),          intent(IN)    :: A_glen(:,:)    ! Glen's-law rate factor
        real(wp),          intent(IN)    :: time

        ! dp scratch arrays for the K24 boundary
        real(dp), allocatable :: H_ice_dp(:,:), z_bed_dp(:,:), mask_dp(:,:)
        real(dp), allocatable :: bmb_w_dp(:,:), uxy_b_dp(:,:), A_glen_dp(:,:)
        real(dp), allocatable :: kappa_dp(:,:)
        real(dp), allocatable :: q_x_dp(:,:), q_y_dp(:,:), N_dp(:,:), p_w_dp(:,:)

        integer  :: nx, ny
        real(wp) :: dt_step

        select case (hyd%par%method)

            case (HYDRO_METHOD_NONE)
                ! Frozen state: do not touch any field, time, or dt.
                return

            case (HYDRO_METHOD_K24)

                dt_step = max(time - hyd%now%time, 0.0_wp)
                hyd%now%time = time
                hyd%now%dt   = dt_step

                if (dt_step == 0.0_wp) return

                nx = size(H_ice,1)
                ny = size(H_ice,2)

                allocate(H_ice_dp(nx,ny), z_bed_dp(nx,ny), mask_dp(nx,ny))
                allocate(bmb_w_dp(nx,ny), uxy_b_dp(nx,ny), A_glen_dp(nx,ny))
                allocate(kappa_dp(nx,ny))
                allocate(q_x_dp(nx,ny), q_y_dp(nx,ny), N_dp(nx,ny), p_w_dp(nx,ny))

                H_ice_dp  = real(H_ice,           dp)
                z_bed_dp  = real(z_bed,           dp)
                mask_dp   = real(mask,            dp)
                bmb_w_dp  = real(bmb_w,           dp)
                uxy_b_dp  = real(uxy_b,           dp)
                A_glen_dp = real(A_glen,          dp)
                kappa_dp  = real(hyd%now%kappa,   dp)

                call calc_k24(q_x_dp, q_y_dp, N_dp, p_w_dp, &
                              H_ice_dp, z_bed_dp, mask_dp, bmb_w_dp, uxy_b_dp, A_glen_dp, &
                              kappa_dp, &
                              real(hyd%par%dx, dp), hyd%par%k24)

                hyd%now%q_x = real(q_x_dp, wp)
                hyd%now%q_y = real(q_y_dp, wp)
                hyd%now%N   = real(N_dp,   wp)
                hyd%now%p_w = real(p_w_dp, wp)

                deallocate(H_ice_dp, z_bed_dp, mask_dp)
                deallocate(bmb_w_dp, uxy_b_dp, A_glen_dp)
                deallocate(kappa_dp)
                deallocate(q_x_dp, q_y_dp, N_dp, p_w_dp)

            case default

                write(*,*) "hydro_update:: error: method must be one of [0,2]."
                write(*,*) "method = ", hyd%par%method
                stop

        end select

        return

    end subroutine hydro_update

    subroutine hydro_par_load(par, filename, init)

        implicit none

        type(hydro_param_class), intent(INOUT) :: par
        character(len=*),        intent(IN)    :: filename
        logical, optional,       intent(IN)    :: init

        logical :: init_pars

        init_pars = .FALSE.
        if (present(init)) init_pars = init

        par%init_method = HYDRO_INIT_EXTERNAL
        par%method      = HYDRO_METHOD_NONE
        par%dx          = 0.0_wp
        par%dy          = 0.0_wp

        call nml_read(filename,"fast_hydrology","init_method", par%init_method, init=init_pars)
        call nml_read(filename,"fast_hydrology","method",      par%method,      init=init_pars)
        call nml_read(filename,"fast_hydrology","dx",          par%dx,          init=init_pars)
        call nml_read(filename,"fast_hydrology","dy",          par%dy,          init=init_pars)

        call k24_par_load(par%k24, filename, init=init_pars)

        return

    end subroutine hydro_par_load

    subroutine hydro_allocate(now, nx, ny)

        implicit none

        type(hydro_state_class), intent(INOUT) :: now
        integer,                 intent(IN)    :: nx, ny

        call hydro_deallocate(now)

        allocate(now%H_w(nx,ny))
        allocate(now%p_w(nx,ny))
        allocate(now%q_x(nx,ny))
        allocate(now%q_y(nx,ny))
        allocate(now%N(nx,ny))
        allocate(now%kappa(nx,ny))

        now%H_w   = 0.0_wp
        now%p_w   = 0.0_wp
        now%q_x   = 0.0_wp
        now%q_y   = 0.0_wp
        now%N     = 0.0_wp
        now%kappa = 0.0_wp

        return

    end subroutine hydro_allocate

    subroutine hydro_deallocate(now)

        implicit none

        type(hydro_state_class), intent(INOUT) :: now

        if (allocated(now%H_w))   deallocate(now%H_w)
        if (allocated(now%p_w))   deallocate(now%p_w)
        if (allocated(now%q_x))   deallocate(now%q_x)
        if (allocated(now%q_y))   deallocate(now%q_y)
        if (allocated(now%N))     deallocate(now%N)
        if (allocated(now%kappa)) deallocate(now%kappa)

        return

    end subroutine hydro_deallocate

end module fast_hydrology
