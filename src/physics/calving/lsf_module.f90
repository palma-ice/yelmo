module lsf_module
    ! ----------------------------------------------------------------------
    ! Level-set function (LSF) for flux-form calving.
    !
    ! Convention: lsf < 0 => ice domain, lsf > 0 => ocean. The zero
    ! level set is the calving front.
    !
    ! Public surface:
    !   - LSFinit         : initialise phi to +-1 from H_ice / z_bed / z_sl
    !   - LSFupdate       : advect phi at w = u_bar + cr (extrapolated into
    !                       the ocean), then saturate to [-1, 1]
    !   - LSFredistance   : Sussman/Osher Hamilton-Jacobi redistancing to
    !                       restore |grad phi| ~= 1 without moving the
    !                       zero level set. Replaces the older ad-hoc
    !                       neighbour-snap and periodic ±1 re-flag.
    !
    ! Mirrors the design in Yelmo.jl/src/topo/lsf.jl.
    ! ----------------------------------------------------------------------

    use yelmo_defs,        only : sp, dp, wp, prec, TOL, TOL_UNDERFLOW, MISSING_VALUE, io_unit_err
    use yelmo_tools,       only : boundary_code, get_neighbor_indices_bc_codes
    use topography,        only : calc_H_eff
    use solver_advection,  only : calc_advec2D

    implicit none

    private

    ! === LSF routines ===
    public :: LSFinit
    public :: LSFupdate
    public :: LSFredistance

    ! === Ocean extrapolation routines ===
    public :: extrapolate_ocn_acx
    public :: extrapolate_ocn_acy

contains
    ! ===================================================================
    !
    !                        LSF functions
    !
    ! ===================================================================

    subroutine LSFinit(LSF,H_ice,z_bed,z_sl,dx)

        implicit none

        real(wp), intent(OUT) :: LSF(:,:)       ! LSF mask
        real(wp), intent(IN)  :: H_ice(:,:)     ! Ice thickness
        real(wp), intent(IN)  :: z_bed(:,:)     ! Bedrock elevation
        real(wp), intent(IN)  :: z_sl(:,:)      ! Sea level
        real(wp), intent(IN)  :: dx             ! Model resolution

        ! Initialize LSF value at ocean value
        LSF = 1.0_wp

        ! Assign values
        where(H_ice .gt. 0.0_wp) LSF = -1.0_wp
        where(z_bed .gt. z_sl)   LSF = -1.0_wp

        return

    end subroutine LSFinit

    subroutine LSFupdate(dlsf,lsf,cr_acx,cr_acy,u_acx,v_acy,mask_adv,dx,dy,dt,solver,boundaries)

        implicit none

        real(wp),       intent(INOUT) :: dlsf(:,:)               ! advected LSF field
        real(wp),       intent(INOUT) :: lsf(:,:)                ! LSF to be advected (aa-nodes)
        real(wp),       intent(INOUT) :: cr_acx(:,:),cr_acy(:,:) ! [m/yr] calving rate (vertical)
        real(wp),       intent(IN)    :: u_acx(:,:)              ! [m/a] 2D velocity, x-direction (ac-nodes)
        real(wp),       intent(IN)    :: v_acy(:,:)              ! [m/a] 2D velocity, y-direction (ac-nodes)
        integer,        intent(IN)    :: mask_adv(:,:)           ! Advection mask
        real(wp),       intent(IN)    :: dx                      ! [m] Horizontal resolution, x-direction
        real(wp),       intent(IN)    :: dy                      ! [m] Horizontal resolution, y-direction
        real(wp),       intent(IN)    :: dt                      ! [a]   Timestep
        character(len=*), intent(IN)  :: solver                  ! Solver to use for the LSF advection equation
        character(len=*), intent(IN)  :: boundaries              ! Boundary condition string (passed to advection)

        ! Local variables
        integer  :: nx, ny
        real(wp), allocatable :: wx(:,:), wy(:,:), mask_lsf(:,:)
        real(wp), allocatable :: var_dot(:,:)                    ! [dvar/dt] Source term for variable. Not used in LSF.

        nx = size(lsf,1)
        ny = size(lsf,2)
        allocate(wx(nx,ny))
        allocate(wy(nx,ny))
        allocate(mask_lsf(nx,ny))
        allocate(var_dot(nx,ny))

        ! Initialize variables
        dlsf     = 0.0_wp  ! LSF change in a time dt
        mask_lsf = 1.0_wp  ! Allow all LSF mask to be advected
        var_dot  = 0.0_wp

        ! Net LSF velocity: dynamic velocity + calving retreat rate
        wx = u_acx + cr_acx
        wy = v_acy + cr_acy

        ! Extrapolate LSF velocities outside of the ice domain so that
        ! upwind advection near the front sees a non-zero front velocity.
        call extrapolate_ocn_acx(wx,wx,u_acx)
        call extrapolate_ocn_acy(wy,wy,v_acy)

        ! Compute the advected LSF field
        call calc_advec2D(dlsf,lsf,mask_lsf,wx,wy,var_dot, &
                            mask_adv,dx,dy,dt,solver,boundaries)
        call apply_tendency_lsf(lsf,dlsf,dt,adjust_lsf=.FALSE.)

        ! Saturate to [-1, 1] as a guardrail against upwind diffusion.
        ! LSFredistance is what restores |grad phi| ~= 1; this just keeps
        ! the field bounded.
        where(lsf .gt.  1.0_wp) lsf =  1.0_wp
        where(lsf .lt. -1.0_wp) lsf = -1.0_wp

        return

    end subroutine LSFupdate

    subroutine LSFredistance(lsf,dx,dy,n_iter,boundaries)
        ! Sussman/Osher Hamilton-Jacobi redistancing.
        !
        ! Solves   d phi / d tau + sgn(phi0) (|grad phi| - 1) = 0
        !
        ! discretised with the Godunov upwind scheme for |grad phi| and a
        ! smoothed sign function sgn(phi0) = phi0 / sqrt(phi0^2 + eps^2),
        ! eps = max(dx, dy). phi0 is the LSF at the start of redistancing
        ! and is held fixed for the duration of the iteration; this is
        ! what keeps the zero level set from drifting.
        !
        ! Pseudo-timestep dtau = 0.5 * min(dx, dy) satisfies the CFL limit
        ! for the explicit Godunov scheme. n_iter ~ 5 is enough for
        ! near-saturated input fields.
        !
        ! Port of lsf_redistance! in Yelmo.jl/src/topo/lsf.jl.

        implicit none

        real(wp),         intent(INOUT) :: lsf(:,:)
        real(wp),         intent(IN)    :: dx
        real(wp),         intent(IN)    :: dy
        integer,          intent(IN)    :: n_iter
        character(len=*), intent(IN)    :: boundaries

        ! Local variables
        integer  :: i, j, n, nx, ny
        integer  :: im1, ip1, jm1, jp1
        integer  :: BC
        real(wp) :: eps, dtau
        real(wp) :: a, b, c, d
        real(wp) :: s0, sgn
        real(wp) :: gx2, gy2, grad
        real(wp), allocatable :: phi0(:,:), new_phi(:,:)

        if (n_iter .le. 0) return

        nx = size(lsf,1)
        ny = size(lsf,2)
        allocate(phi0(nx,ny))
        allocate(new_phi(nx,ny))

        BC   = boundary_code(boundaries)
        eps  = max(dx,dy)
        dtau = 0.5_wp * min(dx,dy)

        ! Freeze the sign field at the start of redistancing. Using the
        ! evolving lsf here would let the zero level set drift.
        phi0 = lsf

        do n = 1, n_iter

            !$omp parallel do collapse(2) default(shared) &
            !$omp private(i,j,im1,ip1,jm1,jp1,a,b,c,d,s0,sgn,gx2,gy2,grad)
            do j = 1, ny
            do i = 1, nx

                call get_neighbor_indices_bc_codes(im1,ip1,jm1,jp1,i,j,nx,ny,BC)

                ! One-sided differences
                a = (lsf(i,j)   - lsf(im1,j)) / dx
                b = (lsf(ip1,j) - lsf(i,j)  ) / dx
                c = (lsf(i,j)   - lsf(i,jm1)) / dy
                d = (lsf(i,jp1) - lsf(i,j)  ) / dy

                ! Smoothed sign of phi0
                s0  = phi0(i,j)
                sgn = s0 / sqrt(s0*s0 + eps*eps)

                ! Godunov upwind |grad phi|^2
                if (sgn .gt. 0.0_wp) then
                    gx2 = max(max(a, 0.0_wp)**2, min(b, 0.0_wp)**2)
                    gy2 = max(max(c, 0.0_wp)**2, min(d, 0.0_wp)**2)
                else
                    gx2 = max(min(a, 0.0_wp)**2, max(b, 0.0_wp)**2)
                    gy2 = max(min(c, 0.0_wp)**2, max(d, 0.0_wp)**2)
                end if
                grad = sqrt(gx2 + gy2)

                new_phi(i,j) = lsf(i,j) - dtau * sgn * (grad - 1.0_wp)

            end do
            end do
            !$omp end parallel do

            lsf = new_phi

        end do

        deallocate(phi0)
        deallocate(new_phi)

        return

    end subroutine LSFredistance

    ! ===================================================================
    !
    ! Internal functions
    !
    ! ===================================================================

    subroutine apply_tendency_lsf(lsf,lsf_dot,dt,adjust_lsf)

        implicit none

        real(wp), intent(INOUT) :: lsf(:,:)
        real(wp), intent(INOUT) :: lsf_dot(:,:)
        real(wp), intent(IN)    :: dt
        logical, optional, intent(IN) :: adjust_lsf

        ! Local variables
        integer :: i, j, nx, ny
        real(wp) :: lsf_prev
        real(wp) :: dlsfdt
        logical  :: allow_adjust_lsf

        if (dt .gt. 0.0) then
            ! Only apply this routine if dt > 0!

            allow_adjust_lsf = .FALSE.
            if (present(adjust_lsf)) allow_adjust_lsf = adjust_lsf

            nx = size(lsf,1)
            ny = size(lsf,2)

            !$omp parallel do collapse(2) private(i,j,lsf_prev,dlsfdt)
            do j = 1, ny
            do i = 1, nx

                ! Store previous value
                lsf_prev = lsf(i,j)

                ! Now update lsf with tendency for this timestep
                lsf(i,j) = lsf_prev + dt*lsf_dot(i,j)

                ! Calculate actual current rate of change
                dlsfdt = (lsf(i,j) - lsf_prev) / dt

                ! Update lsf rate to match rate of change perfectly
                if (allow_adjust_lsf) then
                    lsf_dot(i,j) = dlsfdt
                end if

            end do
            end do
            !$omp end parallel do

        end if

        return

    end subroutine apply_tendency_lsf

    ! ===================================================================
    !
    ! Oceanic extrapolation routines.
    !
    ! ===================================================================

    subroutine extrapolate_ocn_acx(mask_fill,mask_orig,mask_ac)
        ! Fill ocean cells along the x-axis by nearest-filled-neighbour
        ! sweep. A cell is treated as "ocean" if mask_ac == 0 there.
        ! Single forward+backward pass per row: O(nx*ny) total work.

        implicit none

        real(wp), intent(INOUT) :: mask_fill(:,:)
        real(wp), intent(IN)    :: mask_orig(:,:)
        real(wp), intent(IN)    :: mask_ac(:,:)

        ! Local variables
        integer :: i, j, nx, ny
        logical, allocatable :: filled(:,:)

        nx = size(mask_orig,1)
        ny = size(mask_orig,2)
        allocate(filled(nx,ny))

        filled    = mask_ac .ne. 0.0_wp
        mask_fill = mask_orig

        if (sum(mask_orig) .eq. 0.0_wp) then
            deallocate(filled)
            return
        end if

        do j = 1, ny
            ! Forward sweep: rightward fill.
            do i = 2, nx
                if (.not. filled(i,j) .and. filled(i-1,j)) then
                    mask_fill(i,j) = mask_fill(i-1,j)
                    filled(i,j)    = .true.
                end if
            end do
            ! Backward sweep: leftward fill (catches any unfilled left tails).
            do i = nx-1, 1, -1
                if (.not. filled(i,j) .and. filled(i+1,j)) then
                    mask_fill(i,j) = mask_fill(i+1,j)
                    filled(i,j)    = .true.
                end if
            end do
        end do

        deallocate(filled)

        return

    end subroutine extrapolate_ocn_acx

    subroutine extrapolate_ocn_acy(mask_fill,mask_orig,mask_ac)
        ! Same as extrapolate_ocn_acx but along the y-axis.

        implicit none

        real(wp), intent(INOUT) :: mask_fill(:,:)
        real(wp), intent(IN)    :: mask_orig(:,:)
        real(wp), intent(IN)    :: mask_ac(:,:)

        ! Local variables
        integer :: i, j, nx, ny
        logical, allocatable :: filled(:,:)

        nx = size(mask_orig,1)
        ny = size(mask_orig,2)
        allocate(filled(nx,ny))

        filled    = mask_ac .ne. 0.0_wp
        mask_fill = mask_orig

        if (sum(mask_orig) .eq. 0.0_wp) then
            deallocate(filled)
            return
        end if

        do i = 1, nx
            ! Forward sweep: upward fill.
            do j = 2, ny
                if (.not. filled(i,j) .and. filled(i,j-1)) then
                    mask_fill(i,j) = mask_fill(i,j-1)
                    filled(i,j)    = .true.
                end if
            end do
            ! Backward sweep: downward fill.
            do j = ny-1, 1, -1
                if (.not. filled(i,j) .and. filled(i,j+1)) then
                    mask_fill(i,j) = mask_fill(i,j+1)
                    filled(i,j)    = .true.
                end if
            end do
        end do

        deallocate(filled)

        return

    end subroutine extrapolate_ocn_acy

end module lsf_module
