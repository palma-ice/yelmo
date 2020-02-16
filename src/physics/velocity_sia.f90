module velocity_sia 
    
    use yelmo_defs, only : sp, dp, prec 

    implicit none 

!     ! Internal constants
!     integer,  parameter :: dp  = kind(1.d0)
!     integer,  parameter :: sp  = kind(1.0)

!     ! Choose the precision of the library (sp,dp)
!     integer,  parameter :: prec = sp 

    private
    public :: calc_diffusivity_2D
    public :: calc_uxy_sia_2D
    public :: calc_uxy_sia_3D 
    public :: calc_uxy_b_sia 

contains 

    subroutine calc_diffusivity_2D(dd_ab,H_ice,dzsdx,dzsdy,ATT,zeta_aa,dx,n_glen,rho_ice,g)
        ! Calculate the depth-averaged horizontal diffusivity (ab-nodes)
        ! using the SIA approximation
        
        ! H_ice, ATT are on aa-nodes
        ! ux, uy, dzsdx, dzsdy are on ac-nodes 
        ! Intermediate values: Diffusivity calculated on B nodes
        ! Outputs are staggered (defined at boundaries of cell, ARAWAKA-C grid)

        ! Note: These routines would be faster if ATT_int were
        ! passed as an argument (and only calculated when ATT is updated rather
        ! than each dynamic timestep). However, for completeness, here the 
        ! subroutine takes ATT as an argument, and ATT_int is calculated internally
        ! below. 
        
        ! ajr: Note: this produces the correct value of diffusivity, needs checking though. 

        implicit none

        real(prec), intent(OUT) :: dd_ab(:,:)           ! [m^2/a] SIA diffusivity, ab-nodes
        real(prec), intent(IN)  :: H_ice(:,:)           ! [m]   Ice thickness
        real(prec), intent(IN)  :: dzsdx(:,:)           ! [m/m] Surface gradient x-direction 
        real(prec), intent(IN)  :: dzsdy(:,:)           ! [m/m] Surface gradient y-direction
        real(prec), intent(IN)  :: ATT(:,:,:)           ! nx,ny,nz_aa [a-1 Pa-3] Rate factor
        real(prec), intent(IN)  :: zeta_aa(:)           ! [--]  Height vector 0:1 
        real(prec), intent(IN)  :: dx                   ! [m]
        real(prec), intent(IN)  :: n_glen
        real(prec), intent(IN)  :: rho_ice              ! [kg m-3] Ice density 
        real(prec), intent(IN)  :: g                    ! [m s-2]  Gravity acceleration  

        ! Local variables
        integer :: i, j, k, nx, ny, nz_aa 
        real(prec), allocatable :: H_ice_ab(:,:) 
        real(prec), allocatable :: ATT_ab(:,:,:)
        real(prec), allocatable :: ATT_int_ab(:,:)
        real(prec), allocatable :: slope_ab(:,:) 
        
        nx    = size(H_ice,1)
        ny    = size(H_ice,2)
        nz_aa = size(zeta_aa,1) 

        allocate(H_ice_ab(nx,ny))
        allocate(ATT_ab(nx,ny,nz_aa))
        allocate(ATT_int_ab(nx,ny))
        allocate(slope_ab(nx,ny))

        ! Calculate the ice thickness onto the ab-nodes 
        H_ice_ab = stagger_aa_ab(H_ice)

        H_ice_ab(nx,:) = H_ice_ab(nx-1,:) 
        H_ice_ab(:,ny) = H_ice_ab(:,ny-1)
        
        ! Stagger the rate factor onto the ab-nodes 
        do k = 1, nz_aa 
            ATT_ab(:,:,k) = stagger_aa_ab(ATT(:,:,k))
        end do 

        ATT_ab(nx,:,:) = ATT_ab(nx-1,:,:) 
        ATT_ab(:,ny,:) = ATT_ab(:,ny-1,:)
        
        ! Calculate vertically integrated rate factor
        do j = 1, ny
        do i = 1, nx
            ATT_int_ab(i,j) = integrate_trapezoid1D_pt(ATT_ab(i,j,:)*(H_ice(i,j)*(1.0-zeta_aa))**(n_glen+1),zeta_aa)
        end do
        end do
        
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
        
        ! Calculate diffusivity constant
        dd_ab = 2.0 * (rho_ice*g)**n_glen * slope_ab**(n_glen-1.0) * ATT_int_ab

        return 

    end subroutine calc_diffusivity_2D 

    subroutine calc_diffusivity_2D_helper(dd_ab,H_ice,dzsdx,dzsdy,ATT,zeta_aa,dx,n_glen,rho_ice,g)
        ! Calculate the depth-averaged horizontal diffusivity (Ab nodes)
        ! using the SIA approximation
        
        ! H_ice, ATT are on aa-nodes
        ! ux, uy, dzsdx, dzsdy are on ac-nodes 
        ! Intermediate values: Diffusivity calculated on B nodes
        ! Outputs are staggered (defined at boundaries of cell, ARAWAKA-C grid)

        ! Note: These routines would be faster if ATT_int were
        ! passed as an argument (and only calculated when ATT is updated rather
        ! than each dynamic timestep). However, for completeness, here the 
        ! subroutine takes ATT as an argument, and ATT_int is calculated internally
        ! below. 
        
        ! ajr: note this produces a useful value for the SIA solver, but it is not 2D diffusivity.

        implicit none

        real(prec), intent(OUT) :: dd_ab(:,:)           ! [m^2/a] SIA diffusivity, ab-nodes
        real(prec), intent(IN)  :: H_ice(:,:)           ! [m]   Ice thickness
        real(prec), intent(IN)  :: dzsdx(:,:)           ! [m/m] Surface gradient x-direction 
        real(prec), intent(IN)  :: dzsdy(:,:)           ! [m/m] Surface gradient y-direction
        real(prec), intent(IN)  :: ATT(:,:,:)           ! nx,ny,nz_aa [a-1 Pa-3] Rate factor
        real(prec), intent(IN)  :: zeta_aa(:)           ! [--]  Height vector 0:1 
        real(prec), intent(IN)  :: dx                   ! [m]
        real(prec), intent(IN)  :: n_glen
        real(prec), intent(IN)  :: rho_ice              ! [kg m-3] Ice density 
        real(prec), intent(IN)  :: g                    ! [m s-2]  Gravity acceleration  

        ! Local variables
        integer :: i, j, k, nx, ny, nz_aa 
        real(prec), allocatable :: H_ice_ab(:,:) 
        real(prec), allocatable :: ATT_ab(:,:,:)
        real(prec), allocatable :: ATT_int_ab(:,:,:)
        real(prec), allocatable :: ATT_int_ab_mean(:,:)
        real(prec), allocatable :: slope_ab(:,:) 
        
        nx    = size(H_ice,1)
        ny    = size(H_ice,2)
        nz_aa = size(zeta_aa,1) 

        allocate(H_ice_ab(nx,ny))
        allocate(ATT_ab(nx,ny,nz_aa))
        allocate(ATT_int_ab(nx,ny,nz_aa))
        allocate(ATT_int_ab_mean(nx,ny))
        allocate(slope_ab(nx,ny))

        ! Calculate the ice thickness onto the ab-nodes 
        H_ice_ab = stagger_aa_ab(H_ice)

        H_ice_ab(nx,:) = H_ice_ab(nx-1,:) 
        H_ice_ab(:,ny) = H_ice_ab(:,ny-1)
        
        ! Stagger the rate factor onto the ab-nodes 
        do k = 1, nz_aa 
            ATT_ab(:,:,k) = stagger_aa_ab(ATT(:,:,k))
        end do 

        ATT_ab(nx,:,:) = ATT_ab(nx-1,:,:) 
        ATT_ab(:,ny,:) = ATT_ab(:,ny-1,:)
        
        ! Calculate vertically integrated rate factor
        ATT_int_ab      = calc_rate_factor_integrated(ATT_ab,zeta_aa)
        ATT_int_ab_mean = calc_vertical_integrated_2D(ATT_int_ab,zeta_aa)
        
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
        
        ! Calculate diffusivity constant
        dd_ab = 2.0 * (rho_ice*g)**n_glen * H_ice_ab**(n_glen+1.0) &
                    * slope_ab**(n_glen-1.0) * ATT_int_ab_mean

        return 

    end subroutine calc_diffusivity_2D_helper 

    subroutine calc_uxy_sia_2D(ux,uy,H_ice,dzsdx,dzsdy,ATT,zeta_aa,dx,n_glen,rho_ice,g)
        ! Calculate the depth-averaged horizontal velocity field
        ! and intermediate variables using the SIA approximation
        
        ! H_ice, ATT are on aa-nodes
        ! ux, uy, dzsdx, dzsdy are on ac-nodes 
        ! Intermediate values: Diffusivity calculated on B nodes
        ! Outputs are staggered (defined at boundaries of cell, ARAWAKA-C grid)

        ! Note: These routines would be faster if ATT_int were
        ! passed as an argument (and only calculated when ATT is updated rather
        ! than each dynamic timestep). However, for completeness, here the 
        ! subroutine takes ATT as an argument, and ATT_int is calculated internally
        ! below. 
        
        implicit none

        real(prec), intent(OUT) :: ux(:,:)              ! [m/a] SIA velocity x-direction, acx-nodes
        real(prec), intent(OUT) :: uy(:,:)              ! [m/a] SIA velocity y-direction, acy-nodes
        real(prec), intent(IN)  :: H_ice(:,:)           ! [m]   Ice thickness
        real(prec), intent(IN)  :: dzsdx(:,:)           ! [m/m] Surface gradient x-direction 
        real(prec), intent(IN)  :: dzsdy(:,:)           ! [m/m] Surface gradient y-direction
        real(prec), intent(IN)  :: ATT(:,:,:)           ! nx,ny,nz_aa [a-1 Pa-3] Rate factor
        real(prec), intent(IN)  :: zeta_aa(:)          ! [--]  Height vector 0:1 
        real(prec), intent(IN)  :: dx                   ! [m]
        real(prec), intent(IN)  :: n_glen
        real(prec), intent(IN)  :: rho_ice              ! [kg m-3] Ice density 
        real(prec), intent(IN)  :: g                    ! [m s-2]  Gravity acceleration  

        ! Local variables
        integer :: i, j, k, nx, ny, nz_aa 
        real(prec) :: dd_acx, dd_acy 
        real(prec), allocatable :: H_ice_ab(:,:) 
        real(prec), allocatable :: ATT_ab(:,:,:)
        real(prec), allocatable :: ATT_int_ab(:,:,:)
        real(prec), allocatable :: ATT_int_ab_mean(:,:) 
        real(prec), allocatable :: slope_ab(:,:) 
        real(prec), allocatable :: dd_ab(:,:)           ! [m^2/a] SIA diffusivity, ab-nodes
        real(prec) :: H_ac 

        nx    = size(H_ice,1)
        ny    = size(H_ice,2)
        nz_aa = size(zeta_aa,1) 

        allocate(H_ice_ab(nx,ny))
        allocate(ATT_ab(nx,ny,nz_aa))
        allocate(ATT_int_ab(nx,ny,nz_aa))
        allocate(ATT_int_ab_mean(nx,ny))
        allocate(slope_ab(nx,ny))
        allocate(dd_ab(nx,ny))

        ! Calculate the ice thickness onto the ab-nodes 
        H_ice_ab = stagger_aa_ab(H_ice)

        H_ice_ab(nx,:) = H_ice_ab(nx-1,:) 
        H_ice_ab(:,ny) = H_ice_ab(:,ny-1)
        
        ! Stagger the rate factor onto the ab-nodes 
        do k = 1, nz_aa 
            ATT_ab(:,:,k) = stagger_aa_ab(ATT(:,:,k))
        end do 

        ATT_ab(nx,:,:) = ATT_ab(nx-1,:,:) 
        ATT_ab(:,ny,:) = ATT_ab(:,ny-1,:)
        
        ! Calculate vertically integrated rate factor
        ATT_int_ab      = calc_rate_factor_integrated(ATT_ab,zeta_aa)
        ATT_int_ab_mean = calc_vertical_integrated_2D(ATT_int_ab,zeta_aa)
        
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
        
        ! Calculate diffusivity constant
        dd_ab = 2.0 * (rho_ice*g)**n_glen * H_ice_ab**(n_glen+1.0) &
                    * slope_ab**(n_glen-1.0) * ATT_int_ab_mean

        ! Stagger diffusivity constant back from ab- to ac-nodes
        ! and calculate velocity components on ac-nodes 
        ux = 0.0 
        do j=2,ny
        do i=1,nx
            dd_acx  = 0.5*(dd_ab(i,j-1)+dd_ab(i,j))
            ux(i,j) = -dd_acx*dzsdx(i,j)
        end do
        end do
        ux(:,1) = ux(:,2)

        uy = 0.0 
        do j=1,ny
        do i=2,nx
            dd_acy  = 0.5*(dd_ab(i-1,j)+dd_ab(i,j))
            uy(i,j) = -dd_acy*dzsdy(i,j)
        end do
        end do
        uy(1,:) = uy(2,:) 

        return
        
    end subroutine calc_uxy_sia_2D

    subroutine calc_uxy_sia_3D(ux,uy,H_ice,taud_acx,taud_acy,ATT,zeta_aa,dx,n_glen,rho_ice,g)
        ! Calculate the 3D horizontal velocity field
        ! using sia, ssa or hybrid method

        ! Note: These routines would be faster if ATT_int were
        ! passed as an argument (and only calculated when ATT is updated rather
        ! than each dynamic timestep). However, for completeness, here the 
        ! subroutine takes ATT as an argument, and ATT_int is calculated internally
        ! below. 

        implicit none
        
        real(prec), intent(OUT) :: ux(:,:,:)        ! nx,ny,nz_aa [m/a] SIA velocity x-direction, acx-nodes
        real(prec), intent(OUT) :: uy(:,:,:)        ! nx,ny,nz_aa [m/a] SIA velocity y-direction, acy-nodes
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
        real(prec) :: dd_acx, dd_acy  
        real(prec), allocatable :: H_ice_ab(:,:) 
        real(prec), allocatable :: ATT_ab(:,:,:)
        real(prec), allocatable :: ATT_int_ab(:,:,:) 
        real(prec), allocatable :: slope_ab(:,:) 
        real(prec), allocatable :: dd_ab(:,:) 

        real(prec), allocatable :: sigma_tot_ab(:,:) 

        nx    = size(ux,1)
        ny    = size(ux,2)
        nz_aa = size(zeta_aa,1)

        allocate(H_ice_ab(nx,ny))
        allocate(ATT_ab(nx,ny,nz_aa))
        allocate(ATT_int_ab(nx,ny,nz_aa))
        allocate(slope_ab(nx,ny))
        allocate(dd_ab(nx,ny))

        allocate(sigma_tot_ab(nx,ny))

        ! Calculate the ice thickness onto the ab-nodes 
        H_ice_ab   = stagger_aa_ab(H_ice)
        
        H_ice_ab(nx,:) = H_ice_ab(nx-1,:) 
        H_ice_ab(:,ny) = H_ice_ab(:,ny-1)
        
        ! Stagger rate factor onto ab-nodes
        do k = 1, nz_aa 
            ATT_ab(:,:,k) = stagger_aa_ab(ATT(:,:,k))
        end do 

        ATT_ab(nx,:,:) = ATT_ab(nx-1,:,:) 
        ATT_ab(:,ny,:) = ATT_ab(:,ny-1,:)
        
        ! Integrate up to each layer 
        ATT_int_ab = calc_rate_factor_integrated(ATT_ab,zeta_aa)

        ! Get magnitude of driving stress 
        sigma_tot_ab = 0.0 
        do j=1,ny-1
        do i=1,nx-1
            sigma_tot_ab(i,j) = sqrt( (0.5*(taud_acx(i,j)+taud_acx(i,j+1)))**2 &
                                    + (0.5*(taud_acy(i,j)+taud_acy(i+1,j)))**2 )
        end do 
        end do 
        sigma_tot_ab(nx,:) = sigma_tot_ab(nx-1,:) 
        sigma_tot_ab(:,ny) = sigma_tot_ab(:,ny-1)

        ! Reset velocity solution to zero everywhere 
        ux = 0.0 
        uy = 0.0 

        ! Loop over each vertical layer 
        do k = 1, nz_aa 

            ! Calculate quasi-diffusivity for this layer
            dd_ab = 2.0 * H_ice_ab * sigma_tot_ab**(n_glen-1.0) * ATT_int_ab(:,:,k) 
            
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

    subroutine calc_uxy_sia_3D_00(ux,uy,H_ice,dzsdx,dzsdy,ATT,zeta_aa,dx,n_glen,rho_ice,g)
        ! Calculate the 3D horizontal velocity field
        ! using sia, ssa or hybrid method

        ! Note: These routines would be faster if ATT_int were
        ! passed as an argument (and only calculated when ATT is updated rather
        ! than each dynamic timestep). However, for completeness, here the 
        ! subroutine takes ATT as an argument, and ATT_int is calculated internally
        ! below. 

        implicit none
        
        real(prec), intent(OUT) :: ux(:,:,:)        ! nx,ny,nz_aa [m/a] SIA velocity x-direction, acx-nodes
        real(prec), intent(OUT) :: uy(:,:,:)        ! nx,ny,nz_aa [m/a] SIA velocity y-direction, acy-nodes
        real(prec), intent(IN)  :: H_ice(:,:)       ! [m]   Ice thickness 
        real(prec), intent(IN)  :: dzsdx(:,:)       ! [m/m] Surface gradient x-direction 
        real(prec), intent(IN)  :: dzsdy(:,:)       ! [m/m] Surface gradient y-direction 
        real(prec), intent(IN)  :: ATT(:,:,:)       ! nx,ny,nz_aa [a-1 Pa-3] Rate factor
        real(prec), intent(IN)  :: zeta_aa(:)      ! [--]  Height axis 0:1, layer centers (aa-nodes)
        real(prec), intent(IN)  :: dx               ! [m]   Horizontal resolution 
        real(prec), intent(IN)  :: n_glen
        real(prec), intent(IN)  :: rho_ice          ! [kg m-3] Ice density 
        real(prec), intent(IN)  :: g                ! [m s-2]  Gravitational acceleration

        ! Local variables
        integer :: i, j, k, nx, ny, nz_aa
        real(prec) :: dd_acx, dd_acy  
        real(prec), allocatable :: H_ice_ab(:,:) 
        real(prec), allocatable :: ATT_ab(:,:,:)
        real(prec), allocatable :: ATT_int_ab(:,:,:) 
        real(prec), allocatable :: slope_ab(:,:) 
        real(prec), allocatable :: dd_ab(:,:) 

        nx    = size(ux,1)
        ny    = size(ux,2)
        nz_aa = size(zeta_aa,1)

        allocate(H_ice_ab(nx,ny))
        allocate(ATT_ab(nx,ny,nz_aa))
        allocate(ATT_int_ab(nx,ny,nz_aa))
        allocate(slope_ab(nx,ny))
        allocate(dd_ab(nx,ny))

        ! Calculate the ice thickness onto the ab-nodes 
        H_ice_ab   = stagger_aa_ab(H_ice)
        
        H_ice_ab(nx,:) = H_ice_ab(nx-1,:) 
        H_ice_ab(:,ny) = H_ice_ab(:,ny-1)
        
        ! Stagger rate factor onto ab-nodes
        do k = 1, nz_aa 
            ATT_ab(:,:,k) = stagger_aa_ab(ATT(:,:,k))
        end do 

        ATT_ab(nx,:,:) = ATT_ab(nx-1,:,:) 
        ATT_ab(:,ny,:) = ATT_ab(:,ny-1,:)
        
        ! Integrate up to each layer 
        ATT_int_ab = calc_rate_factor_integrated(ATT_ab,zeta_aa)

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

        ! Reset velocity solution to zero everywhere 
        ux = 0.0 
        uy = 0.0 

        ! Loop over each vertical layer 
        do k = 1, nz_aa 

            ! Calculate diffusivity for this layer
            dd_ab = 2.0 * (rho_ice*g)**n_glen * H_ice_ab**(n_glen+1.0)  &
                        * slope_ab**(n_glen-1.0) * ATT_int_ab(:,:,k)

            ! Stagger diffusivity back from Ab to Ac nodes
            ! and calculate velocity components on ac-nodes 
            do j=2,ny
            do i=1,nx
                dd_acx    = 0.5*(dd_ab(i,j-1)+dd_ab(i,j))
                ux(i,j,k) = -dd_acx*dzsdx(i,j)
            end do
            end do
            ux(:,1,k) = ux(:,2,k)

            do j=1,ny
            do i=2,nx
                dd_acy    = 0.5*(dd_ab(i-1,j)+dd_ab(i,j))
                uy(i,j,k) = -dd_acy*dzsdy(i,j)
            end do
            end do
            uy(1,:,k) = uy(2,:,k)

        end do 

        return
        
    end subroutine calc_uxy_sia_3D_00

    subroutine calc_uxy_b_sia(ux_b,uy_b,H_ice,dzsdx,dzsdy,f_pmp,zeta_aa,dx,cf_sia,rho_ice,g)
        ! Calculate the parameterized basal velocity for use with SIA
        ! (following a Weertman-type sliding law)
        
        ! H_ice, ATT are on aa-nodes
        ! ux, uy, dzsdx, dzsdy are on ac-nodes 
        ! Intermediate values: Diffusivity calculated on B nodes
        ! Outputs are staggered (defined at boundaries of cell, ARAWAKA-C grid)

        ! Note: These routines would be faster if ATT_int were
        ! passed as an argument (and only calculated when ATT is updated rather
        ! than each dynamic timestep). However, for completeness, here the 
        ! subroutine takes ATT as an argument, and ATT_int is calculated internally
        ! below. 
        
        implicit none

        real(prec), intent(OUT) :: ux_b(:,:)            ! [m/a] SIA basal velocity x-direction, acx-nodes
        real(prec), intent(OUT) :: uy_b(:,:)            ! [m/a] SIA basal velocity y-direction, acy-nodes
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

        return
        
    end subroutine calc_uxy_b_sia
    
    function calc_rate_factor_integrated(ATT,zeta) result(ATT_int)
        ! Greve and Blatter (2009), Chpt 5, page 82

        implicit none 

        real(prec), intent(IN) :: ATT(:,:,:)
        real(prec), intent(IN) :: zeta(:) 
        real(prec) :: ATT_int(size(ATT,1),size(ATT,2),size(ATT,3))

        ! Local variables 
        integer :: i, j, nx, ny
        real(prec), parameter :: n_glen = 3.0 

        nx = size(ATT,1)
        ny = size(ATT,2) 

        ! Vertically integrated values of ATT to each vertical level
        do j = 1, ny 
        do i = 1, nx 
            ATT_int(i,j,:) = integrate_trapezoid1D_1D(ATT(i,j,:)*(1.0-zeta)**n_glen,zeta)
        end do 
        end do 

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
    
end module velocity_sia 
