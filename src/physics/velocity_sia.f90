module velocity_sia 
    
    use yelmo_defs, only : sp, dp, prec, yelmo_use_omp 

    implicit none 

!     ! Internal constants
!     integer,  parameter :: dp  = kind(1.d0)
!     integer,  parameter :: sp  = kind(1.0)

!     ! Choose the precision of the library (sp,dp)
!     integer,  parameter :: prec = sp 

    private
    public :: calc_velocity_sia
    public :: calc_velocity_basal_sia_00
    public :: calc_velocity_basal_sia       ! Not ready yet, needs checking that it is consistent with Weertman sliding law 
    
    ! ajr: these routines will eventually be private once hybrid-old interface is retired
    public :: calc_dd_ab_3D_serial
    public :: calc_dd_ab_3D_omp 
    public :: calc_uxy_sia_2D
    public :: calc_uxy_sia_3D 
!     public :: calc_uxy_b_sia 

contains 

    subroutine calc_velocity_sia(ux_i,uy_i,ux_i_bar,uy_i_bar,H_ice,taud_acx,taud_acy,ATT,zeta_aa,dx,n_glen,rho_ice,g)

        implicit none 

        real(prec), intent(OUT) :: ux_i(:,:,:) 
        real(prec), intent(OUT) :: uy_i(:,:,:) 
        real(prec), intent(OUT) :: ux_i_bar(:,:) 
        real(prec), intent(OUT) :: uy_i_bar(:,:) 
        real(prec), intent(IN)  :: H_ice(:,:) 
        real(prec), intent(IN)  :: taud_acx(:,:) 
        real(prec), intent(IN)  :: taud_acy(:,:) 
        real(prec), intent(IN)  :: ATT(:,:,:)
        real(prec), intent(IN)  :: zeta_aa(:) 
        real(prec), intent(IN)  :: dx 
        real(prec), intent(IN)  :: n_glen 
        real(prec), intent(IN)  :: rho_ice 
        real(prec), intent(IN)  :: g 

        ! Local variables 
        integer :: nx, ny, nz_aa 
        real(prec), allocatable :: dd_ab(:,:,:) 

        nx    = size(ux_i,1)
        ny    = size(ux_i,2)
        nz_aa = size(ux_i,3)
        
        allocate(dd_ab(nx,ny,nz_aa))

        ! Calculate diffusivity constant on ab-nodes
        call calc_dd_ab_3D_serial(dd_ab,H_ice,taud_acx,taud_acy,ATT,zeta_aa,dx,n_glen,rho_ice,g)
        
        ! ajr: Alternatively, we could use one strategy for serial computation and one for openmp
        ! computation. So far, the routine calc_dd_ab_3D_omp is really slow, so it is not recommended
!         if (yelmo_use_omp) then
!             ! ajr: do not use as it is, really slow: 
!             call calc_dd_ab_3D_omp(dd_ab,H_ice,taud_acx,taud_acy,ATT,zeta_aa,dx,n_glen,rho_ice,g)
!         else 
!             call calc_dd_ab_3D_serial(dd_ab,H_ice,taud_acx,taud_acy,ATT,zeta_aa,dx,n_glen,rho_ice,g)
!         end if 
        
        ! Calculate the 3D horizontal shear velocity fields
        call calc_uxy_sia_3D(ux_i,uy_i,dd_ab,taud_acx,taud_acy)
        
        ! Calculate the depth-averaged horizontal shear velocity fields too
        call calc_uxy_sia_2D(ux_i_bar,uy_i_bar,dd_ab,taud_acx,taud_acy,zeta_aa)

!         Or, simply integrate from 3D velocity field to get depth-averaged field (slower)
!         ux_i_bar = calc_vertical_integrated_2D(ux_i,zeta_aa)
!         uy_i_bar = calc_vertical_integrated_2D(uy_i,zeta_aa)
        
        return 

    end subroutine calc_velocity_sia

    subroutine calc_velocity_basal_sia(ux_b,uy_b,taub_acx,taub_acy,dd_ab_3D,H_ice,taud_acx,taud_acy,f_pmp,zeta_aa,dx,cf_sia,rho_ice,g)
        ! Calculate the parameterized basal velocity for use with SIA
        ! (following a Weertman-type sliding law)
        
        ! H_ice, ATT are on aa-nodes
        ! ux, uy, dzsdx, dzsdy are on ac-nodes 
        ! Intermediate values: Diffusivity calculated on B nodes
        ! Outputs are staggered (defined at boundaries of cell, ARAWAKA-C grid)

        ! Note: This routine works, but is outdated. It should use taud_acx/acy 
        ! instead of dzsdx/dy, as with calc_uxy_sia_2D/3D.
        
        implicit none

        real(prec), intent(OUT) :: ux_b(:,:)            ! [m/a] SIA basal velocity x-direction, acx-nodes
        real(prec), intent(OUT) :: uy_b(:,:)            ! [m/a] SIA basal velocity y-direction, acy-nodes
        real(prec), intent(OUT) :: taub_acx(:,:)        ! [Pa]  Basal stress, x-direction
        real(prec), intent(OUT) :: taub_acy(:,:)        ! [Pa]  Basal stress, y-direction
        real(prec), intent(IN)  :: dd_ab_3D(:,:,:)      ! [m^2/a] SIA diffusivity, ab-nodes
        real(prec), intent(IN)  :: H_ice(:,:)           ! [m]   Ice thickness
        real(prec), intent(IN)  :: taud_acx(:,:)        ! [Pa]  Driving stress, x-direction 
        real(prec), intent(IN)  :: taud_acy(:,:)        ! [Pa]  Driving stress, y-direction
        real(prec), intent(IN)  :: f_pmp(:,:)           ! [--]  Fraction of grid point at pressure melting point
        real(prec), intent(IN)  :: zeta_aa(:)           ! [--]  Height vector 0:1 
        real(prec), intent(IN)  :: dx                   ! [m]
        real(prec), intent(IN)  :: cf_sia               ! [m/a Pa-1]
        real(prec), intent(IN)  :: rho_ice              ! [kg m-3] Ice density 
        real(prec), intent(IN)  :: g                    ! [m s-2]  Gravity acceleration  

        ! Local variables
        integer :: i, j, k, nx, ny, nz_aa 
        real(prec) :: dd_acx, dd_acy 
        real(prec), allocatable :: dd_ab(:,:)           ! [m^2/a] SIA diffusivity, ab-nodes
        real(prec), allocatable :: f_pmp_ab(:,:)
        real(prec) :: H_ac 
        real(prec), parameter :: cf_sia_frozen = 0.0    ! No sliding for purely frozen points 

        nx    = size(H_ice,1)
        ny    = size(H_ice,2)
        nz_aa = size(zeta_aa,1) 

        allocate(dd_ab(nx,ny))
        allocate(f_pmp_ab(nx,ny)) 

        ! Get pressure melting point fraction on ab-nodes 
        f_pmp_ab = stagger_aa_ab(f_pmp)

        f_pmp_ab(nx,:) = f_pmp_ab(nx-1,:) 
        f_pmp_ab(:,ny) = f_pmp_ab(:,ny-1)
        
        ! Calculate vertically integrated diffusivity constant
        dd_ab = calc_vertical_integrated_2D(dd_ab_3D,zeta_aa)

        ! Scale by fraction of point at pressure melting point 
        dd_ab = (f_pmp_ab*cf_sia + (1.0-f_pmp_ab)*cf_sia_frozen) * dd_ab

        ! Stagger diffusivity constant back from ab- to ac-nodes
        ! and calculate velocity components on ac-nodes 
        ux_b = 0.0 
        do j=2,ny
        do i=1,nx
            dd_acx    = 0.5*(dd_ab(i,j-1)+dd_ab(i,j))
            ux_b(i,j) = -dd_acx*taud_acx(i,j)
        end do
        end do
        ux_b(:,1) = ux_b(:,2)

        uy_b = 0.0 
        do j=1,ny
        do i=2,nx
            dd_acy    = 0.5*(dd_ab(i-1,j)+dd_ab(i,j))
            uy_b(i,j) = -dd_acy*taud_acy(i,j)
        end do
        end do
        uy_b(1,:) = uy_b(2,:) 

        ! ajr: to do!!
        ! Diagnose basal stress 
        !call calc_basal_stress(taub_acx,taub_acy,beta_acx,beta_acy,ux_b,uy_b)
        taub_acx = 0.0_prec 
        taub_acy = 0.0_prec 

        return
        
    end subroutine calc_velocity_basal_sia
    
    subroutine calc_velocity_basal_sia_00(ux_b,uy_b,taub_acx,taub_acy,H_ice,dzsdx,dzsdy,f_pmp,zeta_aa,dx,cf_sia,rho_ice,g)
        ! Calculate the parameterized basal velocity for use with SIA
        ! (following a Weertman-type sliding law)
        
        ! H_ice, ATT are on aa-nodes
        ! ux, uy, dzsdx, dzsdy are on ac-nodes 
        ! Intermediate values: Diffusivity calculated on B nodes
        ! Outputs are staggered (defined at boundaries of cell, ARAWAKA-C grid)

        ! Note: This routine works, but is outdated. It should use taud_acx/acy 
        ! instead of dzsdx/dy, as with calc_uxy_sia_2D/3D.
        
        implicit none

        real(prec), intent(OUT) :: ux_b(:,:)            ! [m/a] SIA basal velocity x-direction, acx-nodes
        real(prec), intent(OUT) :: uy_b(:,:)            ! [m/a] SIA basal velocity y-direction, acy-nodes
        real(prec), intent(OUT) :: taub_acx(:,:)        ! [Pa]  Basal stress, x-direction
        real(prec), intent(OUT) :: taub_acy(:,:)        ! [Pa]  Basal stress, y-direction
        real(prec), intent(IN)  :: H_ice(:,:)           ! [m]   Ice thickness
        real(prec), intent(IN)  :: dzsdx(:,:)           ! [m/m] Surface gradient x-direction 
        real(prec), intent(IN)  :: dzsdy(:,:)           ! [m/m] Surface gradient y-direction
        real(prec), intent(IN)  :: f_pmp(:,:)           ! [--]  Fraction of grid point at pressure melting point
        real(prec), intent(IN)  :: zeta_aa(:)           ! [--]  Height vector 0:1 
        real(prec), intent(IN)  :: dx                   ! [m]
        real(prec), intent(IN)  :: cf_sia               ! [m/a Pa-1]
        real(prec), intent(IN)  :: rho_ice              ! [kg m-3] Ice density 
        real(prec), intent(IN)  :: g                    ! [m s-2]  Gravity acceleration  

        ! Local variables
        integer :: i, j, k, nx, ny, nz_aa 
        real(prec) :: dd_acx, dd_acy 
        real(prec), allocatable :: H_ice_ab(:,:) 
        real(prec), allocatable :: slope_ab(:,:) 
        real(prec), allocatable :: dd_ab(:,:)           ! [m^2/a] SIA diffusivity, ab-nodes
        real(prec), allocatable :: f_pmp_ab(:,:)
        real(prec) :: H_ac 
        real(prec), parameter :: cf_sia_frozen = 0.0    ! No sliding for purely frozen points 

        nx    = size(H_ice,1)
        ny    = size(H_ice,2)
        nz_aa = size(zeta_aa,1) 

        allocate(H_ice_ab(nx,ny))
        allocate(slope_ab(nx,ny))
        allocate(dd_ab(nx,ny))
        allocate(f_pmp_ab(nx,ny)) 

        ! Calculate the ice thickness onto the ab-nodes 
        H_ice_ab = stagger_aa_ab(H_ice)

        H_ice_ab(nx,:) = H_ice_ab(nx-1,:) 
        H_ice_ab(:,ny) = H_ice_ab(:,ny-1)
        
        ! Get magnitude of slope on ab-nodes
        slope_ab = 0.0 
        do j=1,ny-1
        do i=1,nx-1
            slope_ab(i,j) = sqrt( (0.5*(dzsdx(i,j)+dzsdx(i,j+1)))**2 &
                                + (0.5*(dzsdy(i,j)+dzsdy(i+1,j)))**2 )
        end do 
        end do 
        slope_ab(nx,:) = slope_ab(nx-1,:) 
        slope_ab(:,ny) = slope_ab(:,ny-1)
        
        ! Get pressure melting point fraction on ab-nodes 
        f_pmp_ab = stagger_aa_ab(f_pmp)

        f_pmp_ab(nx,:) = f_pmp_ab(nx-1,:) 
        f_pmp_ab(:,ny) = f_pmp_ab(:,ny-1)
        
        ! Calculate diffusivity coefficient 
        dd_ab = (f_pmp_ab*cf_sia + (1.0-f_pmp_ab)*cf_sia_frozen) * (rho_ice*g*H_ice_ab) * slope_ab**2

        ! Stagger diffusivity coefficient back from ab- to ac-nodes
        ! and calculate velocity components on ac-nodes 
        ux_b = 0.0 
        do j=2,ny
        do i=1,nx
            dd_acx  = 0.5*(dd_ab(i,j-1)+dd_ab(i,j))
            ux_b(i,j) = -dd_acx*dzsdx(i,j)
        end do
        end do
        ux_b(:,1) = ux_b(:,2)

        uy_b = 0.0 
        do j=1,ny
        do i=2,nx
            dd_acy  = 0.5*(dd_ab(i-1,j)+dd_ab(i,j))
            uy_b(i,j) = -dd_acy*dzsdy(i,j)
        end do
        end do
        uy_b(1,:) = uy_b(2,:) 

        ! ajr: to do!!
        ! Diagnose basal stress 
        !call calc_basal_stress(taub_acx,taub_acy,beta_acx,beta_acy,ux_b,uy_b)
        taub_acx = 0.0_prec 
        taub_acy = 0.0_prec 

        return
        
    end subroutine calc_velocity_basal_sia_00
    
    subroutine calc_dd_ab_3D_omp(dd_ab_3D,H_ice,taud_acx,taud_acy,ATT,zeta_aa,dx,n_glen,rho_ice,g)
        ! Calculate the 3D diffusivity helper field on ab-nodes, as an input 
        ! to the SIA calculation. 

        !$ use omp_lib

        implicit none
        
        real(prec), intent(OUT) :: dd_ab_3D(:,:,:)  ! nx,ny,nz_aa [m/a] Diffusivity helper, ab-nodes
        real(prec), intent(IN)  :: H_ice(:,:)       ! [m]   Ice thickness 
        real(prec), intent(IN)  :: taud_acx(:,:)    ! [Pa] Driving stress x-direction 
        real(prec), intent(IN)  :: taud_acy(:,:)    ! [Pa] Driving stress y-direction 
        real(prec), intent(IN)  :: ATT(:,:,:)       ! nx,ny,nz_aa [a-1 Pa-3] Rate factor
        real(prec), intent(IN)  :: zeta_aa(:)       ! [--]  Height axis 0:1, layer centers (aa-nodes)
        real(prec), intent(IN)  :: dx               ! [m]   Horizontal resolution 
        real(prec), intent(IN)  :: n_glen
        real(prec), intent(IN)  :: rho_ice          ! [kg m-3] Ice density 
        real(prec), intent(IN)  :: g                ! [m s-2]  Gravitational acceleration

        ! Local variables
        integer :: i, j, k, nx, ny, nz_aa
        real(prec), allocatable :: ATT_ab(:)
        real(prec), allocatable :: ATT_int_ab(:) 

        real(prec) :: H_ice_ab, sigma_tot_ab 

        nx    = size(H_ice,1)
        ny    = size(H_ice,2)
        nz_aa = size(zeta_aa,1)

        allocate(ATT_ab(nz_aa))
        allocate(ATT_int_ab(nz_aa))

        ! Reset dd_ab_3D
        dd_ab_3D = 0.0 

        !$omp parallel do 
        do j = 2, ny-1 
        do i = 2, nx-1 

            ! Calculate staggered magnitude of driving stress
            sigma_tot_ab = sqrt( (0.5*(taud_acx(i,j)+taud_acx(i,j+1)))**2 &
                               + (0.5*(taud_acy(i,j)+taud_acy(i+1,j)))**2 )
    
            ! Calculate staggered column of ATT and integrated ATT
            ATT_ab     = 0.25_prec*(ATT(i+1,j+1,:)+ATT(i+1,j,:)+ATT(i,j+1,:)+ATT(i,j,:))
            ATT_int_ab = integrate_trapezoid1D_1D(ATT_ab(:)*(1.0-zeta_aa)**n_glen,zeta_aa)

            ! Calculate staggered ice thickness 
            H_ice_ab   = 0.25_prec*(H_ice(i+1,j+1)+H_ice(i+1,j)+H_ice(i,j+1)+H_ice(i,j))

            ! Calculate quasi-diffusivity for this layer
            dd_ab_3D(i,j,:) = 2.0 * H_ice_ab * ATT_int_ab(:) * sigma_tot_ab**(n_glen-1.0) 
        
        end do 
        end do 
        !$omp end parallel do 

        ! Fill in the borders 
        dd_ab_3D(1,:,:)  = dd_ab_3D(2,:,:) 
        dd_ab_3D(nx,:,:) = dd_ab_3D(nx-1,:,:) 
        dd_ab_3D(:,1,:)  = dd_ab_3D(:,2,:)
        dd_ab_3D(:,ny,:) = dd_ab_3D(:,ny-1,:)

        return
        
    end subroutine calc_dd_ab_3D_omp

    subroutine calc_dd_ab_3D_serial(dd_ab_3D,H_ice,taud_acx,taud_acy,ATT,zeta_aa,dx,n_glen,rho_ice,g)
        ! Calculate the 3D diffusivity helper field on ab-nodes, as an input 
        ! to the SIA calculation. 

        implicit none
        
        real(prec), intent(OUT) :: dd_ab_3D(:,:,:)  ! nx,ny,nz_aa [m/a] Diffusivity helper, ab-nodes
        real(prec), intent(IN)  :: H_ice(:,:)       ! [m]   Ice thickness 
        real(prec), intent(IN)  :: taud_acx(:,:)    ! [Pa] Driving stress x-direction 
        real(prec), intent(IN)  :: taud_acy(:,:)    ! [Pa] Driving stress y-direction 
        real(prec), intent(IN)  :: ATT(:,:,:)       ! nx,ny,nz_aa [a-1 Pa-3] Rate factor
        real(prec), intent(IN)  :: zeta_aa(:)       ! [--]  Height axis 0:1, layer centers (aa-nodes)
        real(prec), intent(IN)  :: dx               ! [m]   Horizontal resolution 
        real(prec), intent(IN)  :: n_glen
        real(prec), intent(IN)  :: rho_ice          ! [kg m-3] Ice density 
        real(prec), intent(IN)  :: g                ! [m s-2]  Gravitational acceleration

        ! Local variables
        integer :: i, j, k, nx, ny, nz_aa
        real(prec), allocatable :: H_ice_ab(:,:) 
        real(prec), allocatable :: ATT_ab(:,:,:)
        real(prec), allocatable :: ATT_int_ab(:,:,:) 
        real(prec), allocatable :: slope_ab(:,:) 

        real(prec), allocatable :: sigma_tot_ab(:,:) 

        nx    = size(H_ice,1)
        ny    = size(H_ice,2)
        nz_aa = size(zeta_aa,1)

        allocate(H_ice_ab(nx,ny))
        allocate(ATT_ab(nx,ny,nz_aa))
        allocate(ATT_int_ab(nx,ny,nz_aa))
        allocate(slope_ab(nx,ny))

        allocate(sigma_tot_ab(nx,ny))

        ! Calculate the ice thickness onto the ab-nodes 
        H_ice_ab   = stagger_aa_ab(H_ice)
        
        H_ice_ab(nx,:) = H_ice_ab(nx-1,:) 
        H_ice_ab(:,ny) = H_ice_ab(:,ny-1)
        
        ! Stagger rate factor onto ab-nodes
        do k = 1, nz_aa 
            ATT_ab(:,:,k) = stagger_aa_ab(ATT(:,:,k))
        end do 

        ATT_ab(1,:,:)  = ATT_ab(2,:,:) 
        ATT_ab(nx,:,:) = ATT_ab(nx-1,:,:) 
        ATT_ab(:,1,:)  = ATT_ab(:,2,:)
        ATT_ab(:,ny,:) = ATT_ab(:,ny-1,:)

        ! Integrate up to each layer 
        ATT_int_ab = calc_rate_factor_integrated(ATT_ab,zeta_aa,n_glen)

        ! Get magnitude of driving stress 
        sigma_tot_ab = 0.0 
        do j=1,ny-1
        do i=1,nx-1
            sigma_tot_ab(i,j) = sqrt( (0.5*(taud_acx(i,j)+taud_acx(i,j+1)))**2 &
                                    + (0.5*(taud_acy(i,j)+taud_acy(i+1,j)))**2 )
        end do 
        end do 
        sigma_tot_ab(1,:)  = sigma_tot_ab(2,:) 
        sigma_tot_ab(nx,:) = sigma_tot_ab(nx-1,:) 
        sigma_tot_ab(:,1)  = sigma_tot_ab(:,2)
        sigma_tot_ab(:,ny) = sigma_tot_ab(:,ny-1)
        
        ! Reset dd_ab_3D
        dd_ab_3D = 0.0 

        ! Loop over each vertical layer 
        do k = 1, nz_aa 

            ! Calculate quasi-diffusivity for this layer
            dd_ab_3D(:,:,k) = 2.0 * H_ice_ab * ATT_int_ab(:,:,k) * sigma_tot_ab**(n_glen-1.0) 

        end do 

        return
        
    end subroutine calc_dd_ab_3D_serial

    subroutine calc_uxy_sia_2D(ux,uy,dd_ab_3D,taud_acx,taud_acy,zeta_aa)
        ! Calculate the 2D horizontal velocity field using SIA

        implicit none

        real(prec), intent(OUT) :: ux(:,:)              ! [m/a] SIA velocity x-direction, acx-nodes
        real(prec), intent(OUT) :: uy(:,:)              ! [m/a] SIA velocity y-direction, acy-nodes
        real(prec), intent(IN)  :: dd_ab_3D(:,:,:)      ! Diffusivity constant 
        real(prec), intent(IN)  :: taud_acx(:,:)        ! [Pa] Driving stress x-direction 
        real(prec), intent(IN)  :: taud_acy(:,:)        ! [Pa] Driving stress y-direction
        real(prec), intent(IN)  :: zeta_aa(:)           ! [--]  Height vector 0:1 

        ! Local variables
        integer :: i, j, k, nx, ny
        real(prec) :: dd_acx, dd_acy 
        real(prec), allocatable :: dd_ab(:,:)           ! [m^2/a] SIA diffusivity, ab-nodes

        nx    = size(ux,1)
        ny    = size(ux,2)

        allocate(dd_ab(nx,ny))

        dd_ab = calc_vertical_integrated_2D(dd_ab_3D,zeta_aa)

        ! Stagger diffusivity constant back from ab- to ac-nodes
        ! and calculate velocity components on ac-nodes 
        ux = 0.0 
        do j=2,ny
        do i=1,nx
            dd_acx  = 0.5*(dd_ab(i,j-1)+dd_ab(i,j))
            ux(i,j) = -dd_acx*taud_acx(i,j)
        end do
        end do
        ux(:,1) = ux(:,2)

        uy = 0.0 
        do j=1,ny
        do i=2,nx
            dd_acy  = 0.5*(dd_ab(i-1,j)+dd_ab(i,j))
            uy(i,j) = -dd_acy*taud_acy(i,j)
        end do
        end do
        uy(1,:) = uy(2,:) 

        return
        
    end subroutine calc_uxy_sia_2D

    subroutine calc_uxy_sia_3D(ux,uy,dd_ab_3D,taud_acx,taud_acy)
        ! Calculate the 3D horizontal velocity field using SIA

        implicit none
        
        real(prec), intent(OUT) :: ux(:,:,:)        ! nx,ny,nz_aa [m/a] SIA velocity x-direction, acx-nodes
        real(prec), intent(OUT) :: uy(:,:,:)        ! nx,ny,nz_aa [m/a] SIA velocity y-direction, acy-nodes
        real(prec), intent(IN)  :: dd_ab_3D(:,:,:)  ! Diffusivity constant
        real(prec), intent(IN)  :: taud_acx(:,:)    ! [Pa] Driving stress x-direction 
        real(prec), intent(IN)  :: taud_acy(:,:)    ! [Pa] Driving stress y-direction
        
        ! Local variables
        integer :: i, j, k, nx, ny, nz_aa
        real(prec) :: dd_acx, dd_acy  
        real(prec), allocatable :: dd_ab(:,:) 

        nx    = size(ux,1)
        ny    = size(ux,2)
        nz_aa = size(ux,3) 

        allocate(dd_ab(nx,ny)) 

        ! Reset velocity solution to zero everywhere 
        ux = 0.0 
        uy = 0.0 

        ! Loop over each vertical layer 
        do k = 1, nz_aa 

            dd_ab = dd_ab_3D(:,:,k) 

            ! Stagger diffusivity back from Ab to Ac nodes
            ! and calculate velocity components on ac-nodes 
            do j=2,ny
            do i=1,nx
                dd_acx    = 0.5*(dd_ab(i,j-1)+dd_ab(i,j))
                ux(i,j,k) = -dd_acx*taud_acx(i,j)
            end do
            end do
            ux(:,1,k) = ux(:,2,k)

            do j=1,ny
            do i=2,nx
                dd_acy    = 0.5*(dd_ab(i-1,j)+dd_ab(i,j))
                uy(i,j,k) = -dd_acy*taud_acy(i,j)
            end do
            end do
            uy(1,:,k) = uy(2,:,k)

        end do 

        return
        
    end subroutine calc_uxy_sia_3D

    function calc_rate_factor_integrated(ATT,zeta,n_glen) result(ATT_int)
        ! Greve and Blatter (2009), Chpt 5, page 82

        !$ use omp_lib

        implicit none 

        real(prec), intent(IN) :: ATT(:,:,:)
        real(prec), intent(IN) :: zeta(:) 
        real(prec), intent(IN) :: n_glen 
        real(prec) :: ATT_int(size(ATT,1),size(ATT,2),size(ATT,3))

        ! Local variables 
        integer :: i, j, nx, ny

        nx = size(ATT,1)
        ny = size(ATT,2) 

        ! Vertically integrated values of ATT to each vertical level

        !!!$omp parallel do shared(nx,ny,ATT,zeta,n_glen) private(i,j,ATT_int)
        !!!$omp parallel do 
        do j = 1, ny 
        do i = 1, nx 
            ATT_int(i,j,:) = integrate_trapezoid1D_1D(ATT(i,j,:)*(1.0-zeta)**n_glen,zeta)
        end do 
        end do 
        !!!$omp end parallel do

        return

    end function calc_rate_factor_integrated
    
    ! Functions from yelmo_tools are duplicated as private functions here
    ! to ensure velocity_sia is a stand-alone module. 
    
    function calc_vertical_integrated_2D(var,zeta) result(var_int)
        ! Vertically integrate a field 3D field (nx,ny,nz) 
        ! to the surface, return a 2D array (nx,ny)
        
        implicit none

        real(prec), intent(IN) :: var(:,:,:)
        real(prec), intent(IN) :: zeta(:)
        real(prec) :: var_int(size(var,1),size(var,2))

        ! Local variables 
        integer :: i, j, nx, ny

        nx = size(var,1)
        ny = size(var,2)

        do j = 1, ny
        do i = 1, nx
            var_int(i,j) = integrate_trapezoid1D_pt(var(i,j,:),zeta)
        end do
        end do

        return

    end function calc_vertical_integrated_2D
    
    function integrate_trapezoid1D_1D(var,sig) result(var_int)
        ! Integrate a variable from the base to each layer zeta of the ice column.
        ! Note this is designed assuming indices 1 = base, nk = surface 
        ! The value of the integral using the trapezium rule can be found using
        ! integral = (b - a)*((f(a) +f(b))/2 + Σ_1_n-1(f(k)) )/n 
        ! Returns a 1D array with integrated value at each level 

        implicit none

        real(prec), intent(IN) :: var(:)
        real(prec), intent(IN) :: sig(:)
        real(prec) :: var_int(size(var,1))

        ! Local variables 
        integer :: k, nk

        nk = size(var,1)

        ! Initial value is zero
        var_int(1:nk) = 0.0_prec 

        ! Intermediate values include sum of all previous values 
        ! Take current value as average between points
        do k = 2, nk
             var_int(k:nk) = var_int(k:nk) + 0.5_prec*(var(k)+var(k-1))*(sig(k) - sig(k-1))
        end do
        
        return

    end function integrate_trapezoid1D_1D

    function integrate_trapezoid1D_pt(var,zeta) result(var_int)
        ! Integrate a variable from the base to height zeta(nk) in the ice column.
        ! The value of the integral using the trapezium rule can be found using
        ! integral = (b - a)*((f(a) +f(b))/2 + Σ_1_n-1(f(k)) )/n 
        ! Returns a point of integrated value of var at level zeta(nk).

        implicit none

        real(prec), intent(IN) :: var(:)
        real(prec), intent(IN) :: zeta(:)
        real(prec) :: var_int

        ! Local variables 
        integer :: k, nk

        nk = size(var,1)

        ! Initial value is zero
        var_int = 0.0_prec 

        ! Intermediate values include sum of all previous values 
        ! Take current value as average between points
        do k = 2, nk
             var_int = var_int + 0.5_prec*(var(k)+var(k-1))*(zeta(k) - zeta(k-1))
        end do

        return

    end function integrate_trapezoid1D_pt

    function stagger_aa_ab(u) result(ustag)
        ! Stagger from Aa => Ab
        ! Four point average from corner aa-nodes to central ab-node 

        implicit none 

        real(prec), intent(IN)  :: u(:,:) 
        real(prec) :: ustag(size(u,1),size(u,2)) 

        ! Local variables 
        integer :: i, j, nx, ny  

        nx = size(u,1)
        ny = size(u,2) 

        ustag = 0.0_prec 

        do j = 1, ny-1 
        do i = 1, nx-1
            ustag(i,j) = 0.25_prec*(u(i+1,j+1)+u(i+1,j)+u(i,j+1)+u(i,j))
        end do 
        end do 

        return

    end function stagger_aa_ab 
    
    function stagger_aa_ab_ice(u,H_ice) result(ustag)
        ! Stagger from Aa => Ab
        ! Four point average from corner Aa nodes to central Ab node 

        implicit none 

        real(prec), intent(IN)  :: u(:,:) 
        real(prec), intent(IN)  :: H_ice(:,:) 
        real(prec) :: ustag(size(u,1),size(u,2)) 

        ! Local variables 
        integer :: i, j, nx, ny, k   

        nx = size(u,1)
        ny = size(u,2) 

        ustag = 0.0_prec 

        do j = 1, ny-1 
        do i = 1, nx-1
            k = 0 
            ustag(i,j) = 0.0 
            if (H_ice(i,j) .gt. 0.0) then 
                ustag(i,j) = ustag(i,j) + u(i,j) 
                k = k+1
            end if 

            if (H_ice(i+1,j) .gt. 0.0) then 
                ustag(i,j) = ustag(i,j) + u(i+1,j) 
                k = k+1 
            end if 
            
            if (H_ice(i,j+1) .gt. 0.0) then 
                ustag(i,j) = ustag(i,j) + u(i,j+1) 
                k = k+1 
            end if 
            
            if (H_ice(i+1,j+1) .gt. 0.0) then 
                ustag(i,j) = ustag(i,j) + u(i+1,j+1) 
                k = k+1 
            end if 
            
            if (k .gt. 0) then 
                ustag(i,j) = ustag(i,j) / real(k,prec)
            end if 

            !ustag(i,j) = 0.25_prec*(u(i+1,j+1)+u(i+1,j)+u(i,j+1)+u(i,j))
        end do 
        end do 

        return

    end function stagger_aa_ab_ice 
    
    subroutine calc_basal_stress(taub_acx,taub_acy,beta_acx,beta_acy,ux_b,uy_b)
        ! Calculate the basal stress resulting from sliding (friction times velocity)
        ! Note: calculated on ac-nodes.
        ! taub [Pa] 
        ! beta [Pa a m-1]
        ! u    [m a-1]
        ! taub = -beta*u 

        implicit none 

        real(prec), intent(OUT) :: taub_acx(:,:)   ! [Pa] Basal stress (acx nodes)
        real(prec), intent(OUT) :: taub_acy(:,:)   ! [Pa] Basal stress (acy nodes)
        real(prec), intent(IN)  :: beta_acx(:,:)   ! [Pa a m-1] Basal friction (acx nodes)
        real(prec), intent(IN)  :: beta_acy(:,:)   ! [Pa a m-1] Basal friction (acy nodes)
        real(prec), intent(IN)  :: ux_b(:,:)       ! [m a-1] Basal velocity (acx nodes)
        real(prec), intent(IN)  :: uy_b(:,:)       ! [m a-1] Basal velocity (acy nodes)
        
        real(prec), parameter :: tol = 1e-3_prec 

        ! Calculate basal stress 
        taub_acx = -beta_acx * ux_b 
        taub_acy = -beta_acy * uy_b 

        ! Avoid underflows
        where(abs(taub_acx) .lt. tol) taub_acx = 0.0_prec 
        where(abs(taub_acy) .lt. tol) taub_acy = 0.0_prec 
        
        return 

    end subroutine calc_basal_stress

end module velocity_sia 
