module lsf_module
    ! Definitions for various calving laws 

    use yelmo_defs,       only : sp, dp, wp, prec, TOL_UNDERFLOW
    use yelmo_tools,      only : get_neighbor_indices
    use topography,       only : calc_H_eff
    use solver_advection, only : calc_adv2D_impl_upwind 

    implicit none 

    private 

    ! LSF routines
    public :: LSFinit
    public :: LSFupdate

contains 

! ===================================================================
!
! Calving related stress/strain routines
!
! ===================================================================



! ===================================================================
!
!                        LSF function
!
! ===================================================================

subroutine LSFinit(LSF,f_ice)

    implicit none

    real(wp), intent(OUT) :: LSF(:,:)   ! [m/yr] Calculated calving rate 
    real(wp), intent(IN)  :: f_ice(:,:) ! [-] Ice area fraction
     
    ! Initialize LSF value at zero
    LSF = -1.0  
    ! Assign values
    where(f_ice .gt. 0.0) LSF = 1.0
    
end subroutine

subroutine LSFupdate(LSFn,LSF,CR,f_ice,ux,uy,var_dot,dx,dy,dt,boundaries)

    implicit none

    real(wp),       intent(OUT)   :: LSFn(:,:)              ! new LSF field
    real(wp),       intent(IN)    :: LSF(:,:)               ! LSF to be advected
    real(wp),       intent(IN)    :: CR(:,:)                ! [m/yr] calving rate (vertical)
    real(wp),       intent(IN)    :: f_ice(:,:)             ! Ice fraction to be advected
    real(wp),       intent(IN)    :: ux(:,:)                ! [m/a] 2D velocity, x-direction (ac-nodes)
    real(wp),       intent(IN)    :: uy(:,:)                ! [m/a] 2D velocity, y-direction (ac-nodes)
    real(wp),       intent(IN)    :: var_dot(:,:)           ! [dvar/dt] Source term for variable (needed?)
    real(wp),       intent(IN)    :: dx                     ! [m] Horizontal resolution, x-direction
    real(wp),       intent(IN)    :: dy                     ! [m] Horizontal resolution, y-direction
    real(wp),       intent(IN)    :: dt                     ! [a]   Timestep 
    character(len=*), intent(IN)  :: boundaries             ! Boundary conditions to impose

    ! Local variables
    integer  :: nx, ny
    real(wp), allocatable :: lsf_now(:,:), wx(:,:), wy(:,:)

    nx = size(LSF,1)
    ny = size(LSF,2)
    allocate(lsf_now(nx,ny))
    allocate(wx(nx,ny))
    allocate(wy(nx,ny))

    ! Advection depends on ice velocity minus the calving rate
    lsf_now = LSF
    wx = ux + CR*ux/SQRT(ux**2.0 + uy**2.0)
    wy = uy + CR*uy/SQRT(uy**2.0 + uy**2.0)

    ! calc_advec2D(dvdt,var,f_ice,ux,uy,var_dot,mask_adv,dx,dy,dt,solver,boundaries)
    ! Compute the advected LSF field
    call calc_adv2D_impl_upwind(lsf_now,wx,wy,var_dot,dx,dy,dt,boundaries,f_upwind=1.0_wp)

    LSFn = lsf_now
    where(LSFn .gt. 1.0) LSFn = 1.0
    where(LSFn .lt. -1.0) LSFn = -1.0

end subroutine LSFupdate

subroutine CalvingAdvectionOcean

    implicit none

end subroutine CalvingAdvectionOcean

end module lsf_module