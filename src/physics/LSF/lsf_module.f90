module lsf_module
    ! Definitions for various calving laws 

    use yelmo_defs,       only : sp, dp, wp, prec, TOL_UNDERFLOW
    use yelmo_tools,      only : get_neighbor_indices
    use topography,       only : calc_H_eff
    use solver_advection, only : calc_advec2D
    use mass_conservation, only : calc_G_advec_simple 

    implicit none 

    private 
    
    ! calving rates
    public :: calc_calving_rate_vonmises_m16

    ! LSF routines
    public :: LSFinit
    public :: LSFupdate

contains 

! ===================================================================
!
! Calving related stress/strain routines
!
! ===================================================================

subroutine calc_calving_rate_vonmises_m16(calv_rate,ubar,f_ice,f_grnd,tau_eff,dx,tau_max)
    ! Calculate the calving rate [m/yr] based on the 
    ! von Mises stress approach, as outlined by Morlighem et al. (2016)
    ! tau_max=150-750 kPa w2=1?0?
    
    implicit none 

    real(wp), intent(OUT) :: calv_rate(:,:)
    real(wp), intent(IN)  :: ubar(:,:)  
    real(wp), intent(IN)  :: f_ice(:,:),f_grnd(:,:)  
    real(wp), intent(IN)  :: tau_eff(:,:)
    real(wp), intent(IN)  :: dx
    real(wp), intent(IN)  :: tau_max

    ! Local variables 
    integer  :: i, j, nx, ny 
    !real(wp), parameter :: calv_lim = 1e6       ! To avoid really high calving values
                                               ! jablasco: set to grid size dx? CHECK
    nx = size(f_ice,1)
    ny = size(f_ice,2)

    ! Assume square grid cells 
    do j = 1, ny
    do i = 1, nx  
        
        ! Compute over the whole domain
        ! Multiplication of eigenvectors
        calv_rate(i,j) = max( ubar(i,j)*tau_eff(i,j)/tau_max, 0.0_wp )
        
        ! Ensure that ice free points do not exhibit calving
        !if ((f_ice(i,j) .eq. 0.0) .or. (f_grnd(i,j) .eq. 1.0)) then
        !    calv_rate(i,j) = 0.0_wp
        !end if
        if (f_ice(i,j) .eq. 0.0) then
            calv_rate(i,j) = 0.0_wp
        end if

        ! jablasco: Apply calving limit?
        !calv_rate(i,j) = min(calv_rate(i,j),calv_lim)

    end do
    end do

    return 

end subroutine calc_calving_rate_vonmises_m16

! ===================================================================
!
!                        LSF function
!
! ===================================================================

subroutine LSFinit(LSF,mask_ice)

    implicit none

    real(wp), intent(OUT) :: LSF(:,:)   ! [m/yr] Calculated calving rate 
    real(wp),  intent(IN) :: mask_ice(:,:) ! [m] Ice tickness
     
    ! Initialize LSF value at 0.0
    LSF = -1.0_wp  
    ! Assign values
    where(mask_ice .gt. 0.0_wp) LSF = 1.0_wp
    !where(f_ice .le. 1e-8) LSF = -1.0
    return
    
end subroutine

subroutine LSFupdate(LSFn,LSF,CR,f_ice,ux,uy,var_dot,mask_adv,dx,dy,dt,solver,boundaries)

    implicit none

    real(wp),       intent(OUT)   :: LSFn(:,:)              ! new LSF field
    real(wp),       intent(IN)    :: LSF(:,:)               ! LSF to be advected
    real(wp),       intent(IN)    :: CR(:,:)                ! [m/yr] calving rate (vertical)
    real(wp),       intent(IN)    :: f_ice(:,:)             ! Ice fraction to be advected
    real(wp),       intent(IN)    :: ux(:,:)                ! [m/a] 2D velocity, x-direction (ac-nodes)
    real(wp),       intent(IN)    :: uy(:,:)                ! [m/a] 2D velocity, y-direction (ac-nodes)
    real(wp),       intent(IN)    :: var_dot(:,:)           ! [dvar/dt] Source term for variable (needed?)
    integer,        intent(IN)    :: mask_adv(:,:)          ! Advection mask
    real(wp),       intent(IN)    :: dx                     ! [m] Horizontal resolution, x-direction
    real(wp),       intent(IN)    :: dy                     ! [m] Horizontal resolution, y-direction
    real(wp),       intent(IN)    :: dt                     ! [a]   Timestep 
    character(len=*), intent(IN)  :: solver                 ! Solver to use for the ice thickness advection equation
    character(len=*), intent(IN)  :: boundaries             ! Boundary conditions to impose

    ! Local variables
    integer  :: nx, ny
    real(wp), allocatable :: dlsf(:,:), wx(:,:), wy(:,:)

    nx = size(LSF,1)
    ny = size(LSF,2)
    allocate(dlsf(nx,ny))
    allocate(wx(nx,ny))
    allocate(wy(nx,ny))

    ! Advection depends on ice velocity minus the calving rate
    dlsf = 0.0
    wx = (ux - CR*ux/SQRT(ux**2.0 + uy**2.0)) ! x-component calving rate
    wy = (uy - CR*uy/SQRT(uy**2.0 + uy**2.0)) ! y-component calving rate
    !wx = (ux - 100*ux)!*ux/SQRT(ux**2.0 + uy**2.0)) ! x-component calving rate
    !wy = (uy - 100*uy)!*uy/SQRT(ux**2.0 + uy**2.0)) ! y-component calving rate
        
    ! Compute the advected LSF field
    call calc_G_advec_simple(dlsf,lsf,f_ice,wx,wy, &
                             mask_adv,solver,boundaries,dx,dt)
    !call calc_advec2D(dlsf,lsf,f_ice,wx,wy,var_dot,mask_adv,dx,dy,dt,solver,boundaries)
    
    LSFn = dlsf*dt + lsf
    where(LSFn .gt. 1.0)  LSFn = 1.0
    where(LSFn .lt. -1.0) LSFn = -1.0
    !where(LSFn .gt. 0.0)  LSFn = 1.0
    !where(LSFn .lt. 0.0) LSFn = -1.0
        

    return

end subroutine LSFupdate

!subroutine CalvingAdvectionOcean
!
!    implicit none
!
!    return
!
!end subroutine CalvingAdvectionOcean

end module lsf_module