module velocity_sia 
    
    use yelmo_defs, only : sp, dp, prec, wp, yelmo_use_omp 

    use yelmo_tools, only : stagger_aa_ab, stagger_aa_ab_ice, stagger_ab_aa_ice, & 
                    stagger_node_aa_ab_ice, stagger_node_acx_ab_ice, stagger_node_acy_ab_ice, &
                    stagger_nodes_aa_ab_ice, stagger_nodes_acx_ab_ice, stagger_nodes_acy_ab_ice, &
                    integrate_trapezoid1D_1D, integrate_trapezoid1D_pt, minmax, &
                    calc_vertical_integrated_2D, calc_vertical_integrated_3D

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

    subroutine calc_velocity_sia(ux_i,uy_i,ux_i_bar,uy_i_bar,H_ice,taud_acx,taud_acy,ATT,zeta_aa,dx,n_glen,rho_ice,g,boundaries)

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
        character(len=*), intent(IN) :: boundaries

        ! Local variables 
        integer :: nx, ny, nz_aa 
        real(prec), allocatable :: dd_ab(:,:,:) 
        real(prec), allocatable :: f_ice(:,:) 
        real(prec), allocatable :: tau_xz(:,:,:) 
        real(prec), allocatable :: tau_yz(:,:,:) 

        nx    = size(ux_i,1)
        ny    = size(ux_i,2)
        nz_aa = size(ux_i,3)
        
        allocate(dd_ab(nx,ny,nz_aa))
        allocate(f_ice(nx,ny))

        allocate(tau_xz(nx,ny,nz_aa))
        allocate(tau_yz(nx,ny,nz_aa))

        ! Define local f_ice 
        f_ice = 0.0 
        where(H_ice .gt. 0.0) f_ice = 1.0 
        
        ! Calculate diffusivity constant on ab-nodes
        ! call calc_dd_ab_3D_serial(dd_ab,H_ice,taud_acx,taud_acy,ATT,zeta_aa,dx,n_glen,rho_ice,g)
        
        ! ajr: Alternatively, we could use one strategy for serial computation and one for openmp
        ! computation. So far, the routine calc_dd_ab_3D_omp is really slow, so it is not recommended
!         if (yelmo_use_omp) then
!             ! ajr: do not use as it is, really slow: 
!             call calc_dd_ab_3D_omp(dd_ab,H_ice,taud_acx,taud_acy,ATT,zeta_aa,dx,n_glen,rho_ice,g)
!         else 
!             call calc_dd_ab_3D_serial(dd_ab,H_ice,taud_acx,taud_acy,ATT,zeta_aa,dx,n_glen,rho_ice,g)
!         end if 
        
        ! Calculate the 3D horizontal shear velocity fields
        !call calc_uxy_sia_3D(ux_i,uy_i,dd_ab,taud_acx,taud_acy)
        
        call calc_shear_stress_3D(tau_xz,tau_yz,taud_acx,taud_acy,f_ice,zeta_aa,boundaries)

        call calc_uxy_sia_3D(ux_i,uy_i,tau_xz,tau_yz,taud_acx,taud_acy,H_ice,f_ice,ATT,n_glen,zeta_aa,boundaries)

        ! Calculate the depth-averaged horizontal shear velocity fields too
        ! call calc_uxy_sia_2D(ux_i_bar,uy_i_bar,dd_ab,taud_acx,taud_acy,zeta_aa)

        ! Or, simply integrate from 3D velocity field to get depth-averaged field (slower)
        ux_i_bar = calc_vertical_integrated_2D(ux_i,zeta_aa)
        uy_i_bar = calc_vertical_integrated_2D(uy_i,zeta_aa)
        
        return 

    end subroutine calc_velocity_sia

    subroutine calc_shear_stress_3D(tau_xz,tau_yz,taud_acx,taud_acy,f_ice,zeta_aa,boundaries)
        ! Calculate the shear stress (x/y) components at each vertical level, as an input 
        ! to the SIA calculation. 

        !$ use omp_lib

        implicit none
        
        real(wp),           intent(OUT) :: tau_xz(:,:,:)        ! nx,ny,nz_aa [m/a] Shear stress 
        real(wp),           intent(OUT) :: tau_yz(:,:,:)        ! nx,ny,nz_aa [m/a] Shear stress 
        real(wp),           intent(IN)  :: taud_acx(:,:)      ! [Pa] Driving stress x-direction 
        real(wp),           intent(IN)  :: taud_acy(:,:)      ! [Pa] Driving stress y-direction 
        real(wp),           intent(IN)  :: f_ice(:,:)         ! [--] Ice area fraction 
        real(wp),           intent(IN)  :: zeta_aa(:)         ! [--] Height axis 0:1, layer centers (aa-nodes)
        character(len=*),   intent(IN)  :: boundaries 

        ! Local variables
        integer  :: i, j, k, nx, ny, nz_aa

        nx    = size(f_ice,1)
        ny    = size(f_ice,2)
        nz_aa = size(zeta_aa,1)

        tau_xz = 0.0_wp 
        tau_yz = 0.0_wp 

        ! Calculate shear stress on ac-nodes at each level (vertical aa-nodes) 
        do k = 1, nz_aa

            !$omp parallel do 
            do j = 1, ny 
            do i = 1, nx 

                    tau_xz(i,j,k) = -(1.0_prec-zeta_aa(k))*taud_acx(i,j)
                    tau_yz(i,j,k) = -(1.0_prec-zeta_aa(k))*taud_acy(i,j)

            end do 
            end do 
            !$omp end parallel do 

        end do 

        ! Apply boundary conditions as needed 
        select case(trim(boundaries))
            
            case("periodic")

                tau_xz(1,:,:)    = tau_xz(nx-2,:,:) 
                tau_xz(nx-1,:,:) = tau_xz(2,:,:) 
                tau_xz(nx,:,:)   = tau_xz(3,:,:) 
                tau_xz(:,1,:)    = tau_xz(:,ny-1,:)
                tau_xz(:,ny,:)   = tau_xz(:,2,:) 

                tau_yz(1,:,:)    = tau_yz(nx-1,:,:) 
                tau_yz(nx,:,:)   = tau_yz(2,:,:) 
                tau_yz(:,1,:)    = tau_yz(:,ny-2,:)
                tau_yz(:,ny-1,:) = tau_yz(:,2,:) 
                tau_yz(:,ny,:)   = tau_yz(:,3,:)

            case("periodic-x") 
            
                tau_xz(1,:,:)    = tau_xz(nx-2,:,:) 
                tau_xz(nx-1,:,:) = tau_xz(2,:,:) 
                tau_xz(nx,:,:)   = tau_xz(3,:,:) 
                tau_xz(:,1,:)    = tau_xz(:,2,:)
                tau_xz(:,ny,:)   = tau_xz(:,ny-1,:) 

                tau_yz(1,:,:)    = tau_yz(nx-1,:,:) 
                tau_yz(nx,:,:)   = tau_yz(2,:,:) 
                tau_yz(:,1,:)    = tau_yz(:,2,:)
                tau_yz(:,ny-1,:) = tau_yz(:,ny-2,:) 
                tau_yz(:,ny,:)   = tau_yz(:,ny-1,:)

            case("infinite") 
            
                tau_xz(1,:,:)    = tau_xz(2,:,:) 
                tau_xz(nx-1,:,:) = tau_xz(nx-2,:,:) 
                tau_xz(nx,:,:)   = tau_xz(nx-1,:,:) 
                tau_xz(:,1,:)    = tau_xz(:,2,:)
                tau_xz(:,ny,:)   = tau_xz(:,ny-1,:) 

                tau_yz(1,:,:)    = tau_yz(2,:,:) 
                tau_yz(nx,:,:)   = tau_yz(nx-1,:,:) 
                tau_yz(:,1,:)    = tau_yz(:,2,:)
                tau_yz(:,ny-1,:) = tau_yz(:,ny-2,:) 
                tau_yz(:,ny,:)   = tau_yz(:,ny-1,:)

        end select

        return
        
    end subroutine calc_shear_stress_3D

    subroutine calc_uxy_sia_3D(ux,uy,tau_xz,tau_yz,taud_acx,taud_acy,H_ice,f_ice,ATT,n_glen,zeta_aa,boundaries)
        ! Calculate the 3D horizontal velocity field using SIA

        implicit none
        
        real(wp),           intent(OUT) :: ux(:,:,:)        ! nx,ny,nz_aa [m/a] SIA velocity x-direction, acx-nodes
        real(wp),           intent(OUT) :: uy(:,:,:)        ! nx,ny,nz_aa [m/a] SIA velocity y-direction, acy-nodes
        real(wp),           intent(IN)  :: tau_xz(:,:,:)    ! [Pa] Shear stress, x-direction
        real(wp),           intent(IN)  :: tau_yz(:,:,:)    ! [Pa] Shear stress, y-direction
        real(wp),           intent(IN)  :: taud_acx(:,:)    ! [Pa] Driving stress x-direction 
        real(wp),           intent(IN)  :: taud_acy(:,:)    ! [Pa] Driving stress y-direction
        real(wp),           intent(IN)  :: H_ice(:,:)       ! [m] Ice thickness 
        real(wp),           intent(IN)  :: f_ice(:,:)       ! [--] Ice area fraction 
        real(wp),           intent(IN)  :: ATT(:,:,:)       ! nx,ny,nz_aa [a-1 Pa-3] Rate factor
        real(wp),           intent(IN)  :: n_glen 
        real(wp),           intent(IN)  :: zeta_aa(:) 
        character(len=*),   intent(IN)  :: boundaries 

        ! Local variables
        integer  :: i, j, k, nx, ny, nz_aa
        integer  :: im1, ip1, jm1, jp1
        real(wp) :: p1    
        real(wp) :: H_ice_ac
        real(wp) :: fact_ac
        real(wp) :: depth 
        real(wp) :: dzeta 

        real(wp) :: tau_xz_ab_up
        real(wp) :: tau_xz_ab_dn
        real(wp) :: tau_xz_ab
        real(wp) :: tau_yz_ab_up
        real(wp) :: tau_yz_ab_dn
        real(wp) :: tau_yz_ab
        real(wp) :: tau_eff_sq_ab
        real(wp) :: ATT_ab_up
        real(wp) :: ATT_ab_dn
        real(wp) :: ATT_ab
        real(wp) :: H_ice_ab

        real(wp), allocatable :: fact_ab(:,:)

        nx    = size(ux,1)
        ny    = size(ux,2)
        nz_aa = size(ux,3) 

        allocate(fact_ab(nx,ny))
        fact_ab = 0.0_wp 

        ! Define exponent 
        p1 = (n_glen-1.0_wp)/2.0_wp

        ! Intialize velocity field to zero 
        ux = 0.0_wp 
        uy = 0.0_wp
        
        ! Loop over layers starting from first layer above the base to surface 
        do k = 2, nz_aa
    
            ! Get non-dimensional thickness of vertical layer between aa-nodes
            dzeta = zeta_aa(k) - zeta_aa(k-1) 

            ! Calculate tau_perp, tau_eff and factor to calculate velocities,
            ! all on ab-nodes 
            !$omp parallel do
            do j = 1, ny 
            do i = 1, nx 
                
                ! Calculate shear stress on horizontal ab-nodes 
                ! and vertical ac-node (center of dzeta layer)
                call stagger_node_acx_ab_ice(tau_xz_ab_up,tau_xz(:,:,k),  f_ice,i,j)
                call stagger_node_acx_ab_ice(tau_xz_ab_dn,tau_xz(:,:,k-1),f_ice,i,j)
                tau_xz_ab = 0.5_wp*(tau_xz_ab_up+tau_xz_ab_dn)

                call stagger_node_acy_ab_ice(tau_yz_ab_up,tau_yz(:,:,k),  f_ice,i,j)
                call stagger_node_acy_ab_ice(tau_yz_ab_dn,tau_yz(:,:,k-1),f_ice,i,j)
                tau_yz_ab = 0.5_wp*(tau_yz_ab_up+tau_yz_ab_dn)

                ! Calculate effective stress
                tau_eff_sq_ab = tau_xz_ab**2 + tau_yz_ab**2

                ! Calculate factor on ab-nodes
                call stagger_node_aa_ab_ice(ATT_ab_up,ATT(:,:,k),f_ice,i,j)
                call stagger_node_aa_ab_ice(ATT_ab_dn,ATT(:,:,k-1),f_ice,i,j)
                ATT_ab = 0.5_wp*(ATT_ab_up+ATT_ab_dn)

                ! Get ice thickness on ab-nodes
                call stagger_node_aa_ab_ice(H_ice_ab,H_ice,f_ice,i,j)

                ! Calculate multiplicative factor on ab-nodes
                if (p1 .ne. 0.0_wp) then 
                    fact_ab(i,j) = 2.0_wp * ATT_ab * (dzeta*H_ice_ab) * tau_eff_sq_ab**p1
                else
                    fact_ab(i,j) = 2.0_wp * ATT_ab * (dzeta*H_ice_ab)
                end if 

            end do 
            end do 
            !$omp end parallel do 

            ! Calculate 3D horizontal velocity components on acx/acy nodes
            !$omp parallel do
            do j = 1, ny 
            do i = 1, nx 
            
                im1 = max(i-1,1) 
                ip1 = min(i+1,nx) 
                jm1 = max(j-1,1) 
                jp1 = min(j+1,ny) 
                
                ! stagger factor to acx-nodes to calculate velocity
                if (f_ice(i,j) .eq. 1.0 .or. f_ice(ip1,j) .eq. 1.0) then 
                    fact_ac   = 0.5_wp*(fact_ab(i,j)+fact_ab(i,jm1))
                    ux(i,j,k) = ux(i,j,k-1) &
                                + fact_ac*0.5_wp*(tau_xz(i,j,k)+tau_xz(i,j,k-1))
                end if 

                ! stagger factor to acy-nodes to calculate velocity
                if (f_ice(i,j) .eq. 1.0 .or. f_ice(i,jp1) .eq. 1.0) then
                    fact_ac   = 0.5_wp*(fact_ab(i,j)+fact_ab(im1,j))
                    uy(i,j,k) = uy(i,j,k-1) &
                                + fact_ac*0.5_wp*(tau_yz(i,j,k)+tau_yz(i,j,k-1))
                end if  

            end do 
            end do 
            !$omp end parallel do 

        end do
        
        ! Apply boundary conditions as needed 
        if (trim(boundaries) .eq. "periodic") then

            ux(1,:,:)    = ux(nx-2,:,:) 
            ux(nx-1,:,:) = ux(2,:,:) 
            ux(nx,:,:)   = ux(3,:,:) 
            ux(:,1,:)    = ux(:,ny-1,:)
            ux(:,ny,:)   = ux(:,2,:) 

            uy(1,:,:)    = uy(nx-1,:,:) 
            uy(nx,:,:)   = uy(2,:,:) 
            uy(:,1,:)    = uy(:,ny-2,:)
            uy(:,ny-1,:) = uy(:,2,:) 
            uy(:,ny,:)   = uy(:,3,:)

        else if (trim(boundaries) .eq. "periodic-x") then 
            
            ux(1,:,:)    = ux(nx-2,:,:) 
            ux(nx-1,:,:) = ux(2,:,:) 
            ux(nx,:,:)   = ux(3,:,:) 
            ux(:,1,:)    = ux(:,2,:)
            ux(:,ny,:)   = ux(:,ny-1,:) 

            uy(1,:,:)    = uy(nx-1,:,:) 
            uy(nx,:,:)   = uy(2,:,:) 
            uy(:,1,:)    = uy(:,2,:)
            uy(:,ny-1,:) = uy(:,ny-2,:) 
            uy(:,ny,:)   = uy(:,ny-1,:)

        else if (trim(boundaries) .eq. "infinite") then 
            
            ux(1,:,:)    = ux(2,:,:) 
            ux(nx-1,:,:) = ux(nx-2,:,:) 
            ux(nx,:,:)   = ux(nx-1,:,:) 
            ux(:,1,:)    = ux(:,2,:)
            ux(:,ny,:)   = ux(:,ny-1,:) 

            uy(1,:,:)    = uy(2,:,:) 
            uy(nx,:,:)   = uy(nx-1,:,:) 
            uy(:,1,:)    = uy(:,2,:)
            uy(:,ny-1,:) = uy(:,ny-2,:) 
            uy(:,ny,:)   = uy(:,ny-1,:)

        end if 

        return
        
    end subroutine calc_uxy_sia_3D

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

    subroutine calc_uxy_sia_3D_0(ux,uy,dd_ab_3D,taud_acx,taud_acy)
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
        
    end subroutine calc_uxy_sia_3D_0

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
