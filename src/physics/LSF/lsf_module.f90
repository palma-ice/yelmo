module lsf_module
    ! Definitions for various calving laws 

    use yelmo_defs,        only : sp, dp, wp, prec, TOL_UNDERFLOW
    use yelmo_tools,       only : get_neighbor_indices
    use topography,        only : calc_H_eff
    use solver_advection,  only : calc_advec2D
    use mass_conservation, only : calc_G_advec_simple 

    implicit none 

    private 
    
    ! Floating calving rates
    public :: calc_calving_rate_vonmises_m16

    ! Grounded calving rates

    ! CalvMIP calving rates
    public :: calvmip_exp1_ac
    public :: calvmip_exp1_aa
    public :: calvmip_exp2_ac_norestart
    public :: calvmip_exp2_ac

    ! LSF routines
    public :: LSFinit
    public :: LSFupdate_aa
    public :: LSFupdate_ac

contains 

! ===================================================================
!
! Calving related stress/strain routines
!
! ===================================================================

subroutine calc_calving_rate_vonmises_m16(clv_rate,uxy_bar,f_ice,f_grnd,tau_eff,dx,tau_max)
    ! Calculate the calving rate [m/yr] based on the 
    ! von Mises stress approach, as outlined by Morlighem et al. (2016)
    ! tau_max=150-750 kPa (fracture strength of ice)
    
    implicit none 

    real(wp), intent(OUT) :: clv_rate(:,:)
    real(wp), intent(IN)  :: uxy_bar(:,:)  
    real(wp), intent(IN)  :: f_ice(:,:),f_grnd(:,:)  
    real(wp), intent(IN)  :: tau_eff(:,:)
    real(wp), intent(IN)  :: dx
    real(wp), intent(IN)  :: tau_max

    ! Local variables 
    integer  :: i, j, nx, ny 
    !real(wp), parameter :: clv_lim = 1e6  ! To avoid really high calving values
                                           ! jablasco: set to grid size dx? CHECK
    nx = size(f_ice,1)
    ny = size(f_ice,2)

    ! Assume square grid cells 
    do j = 1, ny
    do i = 1, nx  
        
        ! Compute over the whole domain
        ! Multiplication of eigenvectors
        clv_rate(i,j) = max( uxy_bar(i,j)*tau_eff(i,j)/tau_max, 0.0_wp )
        
        ! Ensure that ice free points do not exhibit calving
        if (f_ice(i,j) .eq. 0.0) then
            clv_rate(i,j) = 0.0_wp
        end if

        ! jablasco: Apply calving limit?
        !clv_rate(i,j) = min(clv_rate(i,j),clv_lim)

    end do
    end do

    return 

end subroutine calc_calving_rate_vonmises_m16

! ===================================================================
!
!                      CalvMIP experiments
!
! ===================================================================

subroutine calvmip_exp1_ac(crx_ac,cry_ac,ux_ac,vy_ac,f_ice,dx,boundaries)
! Experiment 1 of CalvMIP

    implicit none

    real(wp), intent(OUT) :: crx_ac(:,:),cry_ac(:,:)      ! Calving rates on ac-nodes
    real(wp), intent(IN)  :: ux_ac(:,:),vy_ac(:,:)          ! Velocities on ac-nodes
    real(wp), intent(IN)  :: f_ice(:,:)                     ! Ice fraction
    real(wp), intent(IN)  :: dx
    character(len=*), intent(IN)  :: boundaries             ! Boundary conditions to impose

    ! Local variables
    integer  :: i, j, im1, ip1, jm1, jp1, nx, ny
    real(wp) :: r, ux_aa, vy_aa

    ! points inside the radius 750km should have no calving rate
    nx = size(ux_ac,1)
    ny = size(ux_ac,2)

    r = 0.0_wp
    do j = 1, ny
    do i = 1, nx

        ! aa-nodes indices
        r = sqrt((0.5*(nx+1)-i)*(0.5*(nx+1)-i) + (0.5*(ny+1)-j)*(0.5*(ny+1)-j))*dx

        ! should be 750e3, but put 749e3 to avoid numerical errors
        if (r .le. 749e3) then
            crx_ac(i,j) = 0.0_wp
            cry_ac(i,j) = 0.0_wp

        else
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
            ! x-direction
            if (f_ice(i,j) .gt. 0.0 .and. f_ice(ip1,j) .eq. 0.0) then
                crx_ac(i,j)   = -ux_ac(i,j)
            else if (f_ice(i,j) .gt. 0.0 .and. f_ice(im1,j) .eq. 0.0) then
                crx_ac(im1,j) = -ux_ac(im1,j)
            end if

            ! y-direction
            if (f_ice(i,j) .gt. 0.0 .and. f_ice(i,jp1) .eq. 0.0) then
                cry_ac(i,j)   = -vy_ac(i,j)
            else if (f_ice(i,j) .gt. 0.0 .and. f_ice(i,jm1) .eq. 0.0) then
                cry_ac(i,jm1) = -vy_ac(i,jm1)
            end if

        end if

    end do
    end do

    return

end subroutine calvmip_exp1_ac

subroutine calvmip_exp2_ac_norestart(crx_ac,cry_ac,ux_ac,vy_ac,f_ice,time,dx,boundaries)
! Experiment 2 of CalvMIP with no restart

    implicit none

    real(wp), intent(OUT) :: crx_ac(:,:),cry_ac(:,:)      ! Calving rates on ac-nodes
    real(wp), intent(IN)  :: ux_ac(:,:),vy_ac(:,:)          ! Velocities on ac-nodes
    real(wp), intent(IN)  :: f_ice(:,:)                     ! Ice fraction
    real(wp), intent(IN)  :: dx,time
    character(len=*), intent(IN)  :: boundaries             ! Boundary conditions to impose
    ! Local variables
    integer  :: i, j, im1, ip1, jm1, jp1, nx, ny
    real(wp) :: r, ux_aa, vy_aa
    real(wp), parameter   :: pi = acos(-1.0)  ! Calculate pi intrinsically
    real(wp) :: wv

    ! points inside the radius 750km should have no calving rate
    nx = size(ux_ac,1)
    ny = size(ux_ac,2)

    r = 0.0_wp
    do j = 1, ny
    do i = 1, nx

        if (time .lt. 50e3) then

            ! aa-nodes indices
            r = sqrt((0.5*(nx+1)-i)*(0.5*(nx+1)-i) + (0.5*(ny+1)-j)*(0.5*(ny+1)-j))*dx

            ! should be 750e3, but put 749e3 to avoid numerical errors
            if (r .le. 749e3) then
                crx_ac(i,j) = 0.0_wp
                cry_ac(i,j) = 0.0_wp

            else
                call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
                ! x-direction
                if (f_ice(i,j) .gt. 0.0 .and. f_ice(ip1,j) .eq. 0.0) then
                    crx_ac(i,j)   = -ux_ac(i,j)
                else if (f_ice(i,j) .gt. 0.0 .and. f_ice(im1,j) .eq. 0.0) then
                    crx_ac(im1,j) = -ux_ac(im1,j)
                end if

                ! y-direction
                if (f_ice(i,j) .gt. 0.0 .and. f_ice(i,jp1) .eq. 0.0) then
                    cry_ac(i,j)   = -vy_ac(i,j)
                else if (f_ice(i,j) .gt. 0.0 .and. f_ice(i,jm1) .eq. 0.0) then
                    cry_ac(i,jm1) = -vy_ac(i,jm1)
                end if

            end if

        else

            wv = 300.0 * sin(2.0 * pi * time / 1000.0)
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
           
            ! define over whole domain
            ! x-direction
            if(ux_ac(i,j) .gt. 0.0_wp) then
                crx_ac(i,j) = -(ux_ac(i,j)+wv)
            else
                crx_ac(i,j) = -(ux_ac(i,j)-wv)
                !crx_ac(i,j) = -(ux_ac(im1,j)-wv)
            end if

            ! y-direction
            if(vy_ac(i,j) .gt. 0.0_wp) then
                cry_ac(i,j) = -(vy_ac(i,j)+wv)
            else
                cry_ac(i,j) = -(vy_ac(i,j)-wv)
                !cry_ac(i,jm1) = -(vy_ac(i,jm1)-wv)
            end if

            ! alex
            ! Define at the borders (neccesary?)
            ! x-direction
            !if (f_ice(i,j) .gt. 0.0 .and. f_ice(ip1,j) .eq. 0.0) then
            !    crx_ac(i,j)   = -(ux_ac(i,j)+wv)
            !else if (f_ice(i,j) .gt. 0.0 .and. f_ice(im1,j) .eq. 0.0) then
            !    crx_ac(im1,j) = -(ux_ac(im1,j)-wv)
            !end if

            ! y-direction
            !if (f_ice(i,j) .gt. 0.0 .and. f_ice(i,jp1) .eq. 0.0) then
            !    cry_ac(i,j)   = -(vy_ac(i,j)+wv)
            !else if (f_ice(i,j) .gt. 0.0 .and. f_ice(i,jm1) .eq. 0.0) then
            !    cry_ac(i,jm1) = -(vy_ac(i,jm1)-wv)
            !end if

            ! yo
            ! x-direction
            !if (ux_ac(i,j) .gt. ux_ac(ip1,j)) then
            !    crx_ac(i,j)   = -(ux_ac(i,j)+wv)
            !else if (ux_ac(i,j) .lt. ux_ac(im1,j)) then
            !    crx_ac(im1,j) = -(ux_ac(im1,j)-wv)
            !end if

            ! y-direction
            !if (vy_ac(i,j) .gt. vy_ac(i,jp1)) then
            !    cry_ac(i,j)   = -(vy_ac(i,j)+wv)
            !else if (vy_ac(i,j) .lt. vy_ac(i,jm1)) then
            !    cry_ac(i,jm1) = -(vy_ac(i,jm1)-wv)
            !end if

        end if

    end do
    end do

    return

end subroutine calvmip_exp2_ac_norestart

subroutine calvmip_exp1_aa(clv_rate,ux_ac,vy_ac,dx,boundaries)
! Experiment 1 of CalvMIP

    implicit none

    real(wp), intent(OUT) :: clv_rate(:,:)                  ! calving rates on aa-nodes
    real(wp), intent(IN)  :: ux_ac(:,:),vy_ac(:,:)          ! velocities on ac-nodes
    real(wp), intent(IN)  :: dx
    character(len=*), intent(IN)  :: boundaries             ! Boundary conditions to impose

    ! Local variables
    integer  :: i, j, im1, ip1, jm1, jp1, nx, ny
    real(wp) :: r, ux_aa, vy_aa

    ! points inside the radius 750km should have no calving rate
    nx = size(ux_ac,1)
    ny = size(ux_ac,2)

    r = 0.0_wp
    do j = 1, ny
    do i = 1, nx

        ! wrong indices?
        r = sqrt((0.5*(nx+1)-i)*(0.5*(nx+1)-i) + (0.5*(ny+1)-j)*(0.5*(ny+1)-j))*dx

        ! should be 750e3, but put 749e3 to avoid numerical errors
        if (r .le. 749e3) then
            clv_rate(i,j) = 0.0_wp

        else
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
            ux_aa = 0.5_wp*(ux_ac(i,j)+ux_ac(im1,j))
            vy_aa = 0.5_wp*(vy_ac(i,j)+vy_ac(i,jm1))
            clv_rate(i,j) = SQRT(ux_aa*ux_aa + vy_aa*vy_aa)

        end if

    end do
    end do

    return

end subroutine calvmip_exp1_aa

subroutine calvmip_exp1_exp2_aa(clv_rate,ux_ac,vy_ac,dx,t,boundaries)
! Experiment 1 of CalvMIP

    implicit none

    real(wp), intent(OUT) :: clv_rate(:,:)
    real(wp), intent(IN)  :: ux_ac(:,:),vy_ac(:,:)
    real(wp), intent(IN)  :: dx,t
    character(len=*), intent(IN)  :: boundaries             ! Boundary conditions to impose

    ! Local variables
    integer  :: i, j, im1, ip1, jm1, jp1, nx, ny
    real(wp) :: r, ux_aa, vy_aa
    real(wp), parameter   :: pi = acos(-1.0)  ! Calculate pi intrinsically
    real(wp) :: wv    

    ! points inside the radius 750km should have no calving rate
    nx = size(ux_ac,1)
    ny = size(ux_ac,2)

    r = 0.0_wp
    wv = -300.0 * sin(2.0 * pi * t / 1000.0)
    do j = 1, ny
    do i = 1, nx

        r = sqrt((0.5*nx-i)*(0.5*nx-i) + (0.5*ny-j)*(0.5*ny-j))*dx

        if (t .lt. 50e3) then

            if (r .le. 750e3) then
                clv_rate(i,j) = 0.0_wp

            else
                call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
                ux_aa = 0.5_wp*(ux_ac(i,j)+ux_ac(im1,j))
                vy_aa = 0.5_wp*(vy_ac(i,j)+vy_ac(i,jm1))
                clv_rate(i,j) = SQRT(ux_aa*ux_aa + vy_aa*vy_aa)

            end if

        else
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
            ux_aa = 0.5_wp*(ux_ac(i,j)+ux_ac(im1,j))
            vy_aa = 0.5_wp*(vy_ac(i,j)+vy_ac(i,jm1))
            clv_rate(i,j) = SQRT(ux_aa*ux_aa + vy_aa*vy_aa)-wv 
        end if
    end do
    end do

    return

end subroutine calvmip_exp1_exp2_aa

subroutine calvmip_exp2_ac(crx_ac,cry_ac,ux_ac,vy_ac,time)
    ! Experiment 2 of CalvMIP

    implicit none

    real(wp), intent(OUT) :: crx_ac(:,:), cry_ac(:,:)
    real(wp), intent(IN)  :: ux_ac(:,:),vy_ac(:,:)
    real(wp), intent(IN)  :: time

    ! local variables
    integer :: nx,ny
    real(wp), parameter   :: pi = acos(-1.0)  ! Calculate pi intrinsically
    real(wp), allocatable :: wv(:,:)   

    nx = size(ux_ac,1)
    ny = size(ux_ac,2)
    allocate(wv(nx,ny))

    wv = -300.0 * sin(2.0 * pi * time / 1000.0) 

    crx_ac = -(ux_ac - wv)
    cry_ac = -(vy_ac - wv)

    return

end subroutine calvmip_exp2_ac

! ===================================================================
!
!                        LSF functions
!
! ===================================================================

subroutine LSFinit(LSF,mask_ice,z_bed)

    implicit none

    real(wp), intent(OUT) :: LSF(:,:)      ! LSF mask 
    real(wp), intent(IN)  :: mask_ice(:,:) ! Ice thickness
    real(wp), intent(IN)  :: z_bed(:,:)    ! bedrock elevation     

    ! Initialize LSF value at ocean value
    LSF = -1.0_wp  
    
    ! Assign values
    where(mask_ice .gt. 0.0_wp) LSF = 1.0_wp
    where(z_bed .gt. 0.0_wp) LSF = 1.0_wp    

    return
    
end subroutine LSFinit

subroutine LSFupdate_aa(LSFn,LSF,CR,f_ice,ux_ac,vy_ac,var_dot,mask_adv,dx,dy,dt,solver,boundaries)

    implicit none

    real(wp),       intent(OUT)   :: LSFn(:,:)              ! new LSF field (aa-nodes)
    real(wp),       intent(IN)    :: LSF(:,:)               ! LSF to be advected (aa-nodes)
    real(wp),       intent(IN)    :: CR(:,:)                ! [m/yr] calving rate (vertical)
    real(wp),       intent(IN)    :: f_ice(:,:)             ! Ice fraction to be advected
    real(wp),       intent(IN)    :: ux_ac(:,:)             ! [m/a] 2D velocity, x-direction (ac-nodes)
    real(wp),       intent(IN)    :: vy_ac(:,:)             ! [m/a] 2D velocity, y-direction (ac-nodes)
    real(wp),       intent(IN)    :: var_dot(:,:)           ! [dvar/dt] Source term for variable (needed?)
    integer,        intent(IN)    :: mask_adv(:,:)          ! Advection mask
    real(wp),       intent(IN)    :: dx                     ! [m] Horizontal resolution, x-direction
    real(wp),       intent(IN)    :: dy                     ! [m] Horizontal resolution, y-direction
    real(wp),       intent(IN)    :: dt                     ! [a]   Timestep 
    character(len=*), intent(IN)  :: solver                 ! Solver to use for the ice thickness advection equation
    character(len=*), intent(IN)  :: boundaries             ! Boundary conditions to impose

    ! Local variables
    integer  :: i, j, im1, ip1, jm1, jp1, nx, ny
    real(wp), allocatable :: dlsf(:,:), wx(:,:), wy(:,:)
    real(wp) :: eps, uxy_aa, ux_aa, vy_aa

    nx = size(LSF,1)
    ny = size(LSF,2)
    allocate(dlsf(nx,ny))
    allocate(wx(nx,ny))
    allocate(wy(nx,ny))

    ! Advection depends on ice velocity minus the calving rate
    dlsf = 0.0
    eps  = 1e-6 ! to avoid singularities in the velocity field

    do j = 1,ny
    do i = 1,nx

        if (f_ice(i,j) .eq. 0.0_wp) then
            wx(i,j) = 0.0_wp
            wy(i,j) = 0.0_wp

        else
            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
            ux_aa  = 0.5_wp*(ux_ac(i,j)+ux_ac(im1,j))
            vy_aa  = 0.5_wp*(vy_ac(i,j)+vy_ac(i,jm1))
            uxy_aa = SQRT(ux_aa*ux_aa + vy_aa*vy_aa)

            wx(i,j) = (ux_ac(i,j)  - CR(i,j)*ux_ac(i,j)/(uxy_aa)) ! add eps?
            wy(i,j) = (vy_ac(i,j)  - CR(i,j)*vy_ac(i,j)/(uxy_aa)) ! 

        end if

    end do
    end do
        
    ! Compute the advected LSF field
    call calc_G_advec_simple(dlsf,lsf,f_ice,wx,wy, &
                             mask_adv,solver,boundaries,dx,dt)
    
    ! jablasco: check if the sum here is done correctly
    LSFn = dlsf*dt + lsf

    ! Allow for in between values only in the zero border
    do j = 1,ny
    do i = 1,nx

        ! Get neighbor indices
        call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

        if ((LSFn(i,j) .gt. 0.0_wp) .and. (LSFn(im1,j) .gt. 0.0_wp) .and. (LSFn(ip1,j) .gt. 0.0_wp) &
            .and. (LSFn(i,jm1) .gt. 0.0_wp) .and. (LSFn(i,jp1) .gt. 0.0_wp)) then
            LSFn(i,j) = 1.0_wp

        else if ((LSFn(i,j) .le. 0.0_wp) .and. (LSFn(im1,j) .le. 0.0_wp) .and. (LSFn(ip1,j) .le. 0.0_wp) &
            .and. (LSFn(i,jm1) .le. 0.0_wp) .and. (LSFn(i,jp1) .le. 0.0_wp)) then
            LSFn(i,j) = -1.0_wp        

        end if

    end do
    end do

    ! saturate values to -1 to 1 (helps with stability)
    where(LSFn .gt. 1.0)  LSFn = 1.0
    where(LSFn .lt. -1.0) LSFn = -1.0 

    return

end subroutine LSFupdate_aa

subroutine LSFupdate_ac(LSFn,LSF,crx_ac,cry_ac,f_ice,ux_ac,vy_ac,var_dot,mask_adv,dx,dy,dt,solver,boundaries)

    implicit none

    real(wp),       intent(OUT)   :: LSFn(:,:)               ! new LSF field (aa-nodes)
    real(wp),       intent(IN)    :: LSF(:,:)                ! LSF to be advected (aa-nodes)
    real(wp),       intent(IN)    :: crx_ac(:,:),cry_ac(:,:) ! [m/yr] calving rate (vertical)
    real(wp),       intent(IN)    :: f_ice(:,:)              ! Ice fraction to be advected
    real(wp),       intent(IN)    :: ux_ac(:,:)              ! [m/a] 2D velocity, x-direction (ac-nodes)
    real(wp),       intent(IN)    :: vy_ac(:,:)              ! [m/a] 2D velocity, y-direction (ac-nodes)
    real(wp),       intent(IN)    :: var_dot(:,:)            ! [dvar/dt] Source term for variable (needed?)
    integer,        intent(IN)    :: mask_adv(:,:)           ! Advection mask
    real(wp),       intent(IN)    :: dx                      ! [m] Horizontal resolution, x-direction
    real(wp),       intent(IN)    :: dy                      ! [m] Horizontal resolution, y-direction
    real(wp),       intent(IN)    :: dt                      ! [a]   Timestep
    character(len=*), intent(IN)  :: solver                  ! Solver to use for the ice thickness advection equation
    character(len=*), intent(IN)  :: boundaries              ! Boundary conditions to impose

    ! Local variables
    integer  :: i, j, im1, ip1, jm1, jp1, nx, ny
    real(wp), allocatable :: dlsf(:,:), wx(:,:), wy(:,:)

    nx = size(LSF,1)
    ny = size(LSF,2)
    allocate(dlsf(nx,ny))
    allocate(wx(nx,ny))
    allocate(wy(nx,ny))

    ! Advection depends on ice velocity minus the calving rate
    dlsf = 0.0
    ! calving rate has opposite sign as velocity
    wx = ux_ac + crx_ac
    wy = vy_ac + cry_ac
 
    ! Compute the advected LSF field
    call calc_G_advec_simple(dlsf,lsf,f_ice,wx,wy, &
                             mask_adv,solver,boundaries,dx,dt)

    ! jablasco: check if the sum here is done correctly
    LSFn = dlsf*dt + lsf

    ! Allow for in between values only in the zero border
    !do j = 1,ny
    !do i = 1,nx
!
    !    ! Get neighbor indices
    !    call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
!
    !    if ((LSFn(i,j) .gt. 0.0_wp) .and. (LSFn(im1,j) .gt. 0.0_wp) .and. (LSFn(ip1,j) .gt. 0.0_wp) &
!            .and. (LSFn(i,jm1) .gt. 0.0_wp) .and. (LSFn(i,jp1) .gt. 0.0_wp)) then
!            LSFn(i,j) = 1.0_wp
!
!        else if ((LSFn(i,j) .le. 0.0_wp) .and. (LSFn(im1,j) .le. 0.0_wp) .and. (LSFn(ip1,j) .le. 0.0_wp) &
 !           .and. (LSFn(i,jm1) .le. 0.0_wp) .and. (LSFn(i,jp1) .le. 0.0_wp)) then
  !          LSFn(i,j) = -1.0_wp
!
 !       end if
!
 !   end do
  !  end do

    ! saturate values to -1 to 1 (helps with stability)
    where(LSFn .gt. 1.0)  LSFn = 1.0
    where(LSFn .lt. -1.0) LSFn = -1.0

    return

end subroutine LSFupdate_ac

end module lsf_module
