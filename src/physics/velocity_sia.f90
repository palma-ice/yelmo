module velocity_sia 
    
    use yelmo_defs, only : sp, dp, wp, yelmo_use_omp, TOL_UNDERFLOW 
    
    use yelmo_tools, only : get_neighbor_indices, acx_to_nodes, acy_to_nodes, aa_to_nodes, &
                    acx_to_nodes_3D, acy_to_nodes_3D, aa_to_nodes_3D, &
                    integrate_trapezoid1D_1D, integrate_trapezoid1D_pt, minmax, &
                    calc_vertical_integrated_2D, stagger_aa_ab

    implicit none 

    private
    public :: calc_velocity_sia
    public :: calc_velocity_basal_sia       ! Not ready yet, needs checking that it is consistent with Weertman sliding law 
    
contains 
    
    subroutine calc_velocity_sia(ux_i,uy_i,ux_i_bar,uy_i_bar,H_ice,f_ice,taud_acx,taud_acy,ATT,zeta_aa,dx,n_glen,rho_ice,g,boundaries)

        implicit none 

        real(wp),         intent(OUT) :: ux_i(:,:,:) 
        real(wp),         intent(OUT) :: uy_i(:,:,:) 
        real(wp),         intent(OUT) :: ux_i_bar(:,:) 
        real(wp),         intent(OUT) :: uy_i_bar(:,:) 
        real(wp),         intent(IN)  :: H_ice(:,:) 
        real(wp),         intent(IN)  :: f_ice(:,:)
        real(wp),         intent(IN)  :: taud_acx(:,:) 
        real(wp),         intent(IN)  :: taud_acy(:,:) 
        real(wp),         intent(IN)  :: ATT(:,:,:)
        real(wp),         intent(IN)  :: zeta_aa(:) 
        real(wp),         intent(IN)  :: dx 
        real(wp),         intent(IN)  :: n_glen 
        real(wp),         intent(IN)  :: rho_ice 
        real(wp),         intent(IN)  :: g 
        character(len=*),   intent(IN)  :: boundaries 

        ! Local variables 
        integer :: nx, ny, nz_aa 
        real(wp), allocatable :: tau_xz(:,:,:) 
        real(wp), allocatable :: tau_yz(:,:,:) 

        nx    = size(ux_i,1)
        ny    = size(ux_i,2)
        nz_aa = size(ux_i,3)
        
        allocate(tau_xz(nx,ny,nz_aa))
        allocate(tau_yz(nx,ny,nz_aa))

        ! Calculate the 3D shear stress field on ac-nodes
        call calc_shear_stress_3D(tau_xz,tau_yz,taud_acx,taud_acy,f_ice,zeta_aa,boundaries)

        ! Calculate the 3D horizontal shear velocity fields
        call calc_uxy_sia_3D(ux_i,uy_i,tau_xz,tau_yz,taud_acx,taud_acy,H_ice,f_ice,ATT,n_glen,zeta_aa,boundaries)

        ! Integrate from 3D velocity field to get depth-averaged field
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
        
        !!$omp parallel do collapse(3) private(i,j,k)
        do j = 1, ny 
        do i = 1, nx 

                do k = 1, nz_aa

                    tau_xz(i,j,k) = -(1.0_wp-zeta_aa(k))*taud_acx(i,j)
                    tau_yz(i,j,k) = -(1.0_wp-zeta_aa(k))*taud_acy(i,j)

                end do 

        end do 
        end do 
        !!$omp end parallel do 

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
        integer  :: im1m, ip1m, jm1m, jp1m
        real(wp) :: p1    
        real(wp) :: H_ice_ac
        real(wp) :: fact_ac
        real(wp) :: depth 
        real(wp) :: dzeta 

        real(wp) :: tau_xz_n_up
        real(wp) :: tau_xz_n_dn
        real(wp) :: tau_xz_n
        real(wp) :: tau_yz_n_up
        real(wp) :: tau_yz_n_dn
        real(wp) :: tau_yz_n
        real(wp) :: tau_eff_sq_n
        real(wp) :: ATT_n_up
        real(wp) :: ATT_n_dn
        real(wp) :: ATT_n
        real(wp) :: H_ice_n
        real(wp) :: fact_n

        real(wp), allocatable :: fact_ab(:,:)

        nx    = size(ux,1)
        ny    = size(ux,2)
        nz_aa = size(ux,3) 

        allocate(fact_ab(nx,ny))

        ! Define exponent 
        p1 = (n_glen-1.0_wp)/2.0_wp

        ! Intialize velocity field to zero 
        ux = 0.0_wp 
        uy = 0.0_wp

        ! Loop over layers starting from first layer above the base to surface 
        do k = 2, nz_aa

            ! Get non-dimensional thickness of vertical layer between aa-nodes
            dzeta = zeta_aa(k) - zeta_aa(k-1) 

            ! First calculate sia factor on aa nodes for this layer
            !!$omp parallel do collapse(2) private(i,j,k,im1,ip1,jm1,jp1,tau_xz_n_up,tau_xz_n_dn,tau_xz_n,tau_yz_n_up,tau_yz_n_dn,tau_yz_n,tau_eff_sq_n) &
            !!$omp& private(ATT_n_up,ATT_n_dn,ATT_n,H_ice_n,fact_n)
            do j = 1, ny 
            do i = 1, nx 

                ! Get neighbor indices
                call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
                
                tau_xz_n_up = 0.5_wp*(tau_xz(i,j,k)+tau_xz(i,jp1,k))
                tau_xz_n_dn = 0.5_wp*(tau_xz(i,j,k)+tau_xz(i,jp1,k-1))
                tau_xz_n    = 0.5_wp*(tau_xz_n_up+tau_xz_n_dn)

                tau_yz_n_up = 0.5_wp*(tau_yz(i,j,k)+tau_yz(ip1,j,k))
                tau_yz_n_dn = 0.5_wp*(tau_yz(i,j,k)+tau_yz(ip1,j,k-1))
                tau_yz_n    = 0.5_wp*(tau_yz_n_up+tau_yz_n_dn)

                ! Calculate effective stress
                tau_eff_sq_n = tau_xz_n**2 + tau_yz_n**2

                ATT_n_up = 0.25_wp*(ATT(i,j,k)+ATT(ip1,j,k)+ATT(i,jp1,k)+ATT(ip1,jp1,k))
                ATT_n_dn = 0.25_wp*(ATT(i,j,k-1)+ATT(ip1,j,k-1)+ATT(i,jp1,k-1)+ATT(ip1,jp1,k-1))
                ATT_n = 0.5_wp*(ATT_n_up+ATT_n_dn)

                ! Get ice thickness on node
                H_ice_n = 0.25_wp*(H_ice(i,j)+H_ice(ip1,j)+H_ice(i,jp1)+H_ice(ip1,jp1))
              
                ! Calculate multiplicative factor on node
                if (p1 .ne. 0.0_wp) then 
                    fact_n = 2.0_wp * ATT_n * (dzeta*H_ice_n) * tau_eff_sq_n**p1
                else
                    fact_n = 2.0_wp * ATT_n * (dzeta*H_ice_n)
                end if 

                ! Get fact value on ab node
                fact_ab(i,j) = fact_n

            end do 
            end do 
            !!$omp end parallel do 

            ! Next calculate 3D horizontal velocity components on acx/acy nodes
            !!$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,fact_ac)
            do j = 1, ny 
            do i = 1, nx 
            
                im1 = max(i-1,1) 
                ip1 = min(i+1,nx) 
                jm1 = max(j-1,1) 
                jp1 = min(j+1,ny) 

                fact_ac = 0.5_wp*(fact_ab(i,j)+fact_ab(i,jm1))
                ux(i,j,k) = ux(i,j,k-1) &
                            + fact_ac*0.5_wp*(tau_xz(i,j,k)+tau_xz(i,j,k-1))

                fact_ac = 0.5_wp*(fact_ab(i,j)+fact_ab(im1,j))
                uy(i,j,k) = uy(i,j,k-1) &
                                + fact_ac*0.5_wp*(tau_yz(i,j,k)+tau_yz(i,j,k-1)) 

            end do
            end do
            !!$omp end parallel do 

        end do
        
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

        real(wp), intent(OUT) :: ux_b(:,:)            ! [m/a] SIA basal velocity x-direction, acx-nodes
        real(wp), intent(OUT) :: uy_b(:,:)            ! [m/a] SIA basal velocity y-direction, acy-nodes
        real(wp), intent(OUT) :: taub_acx(:,:)        ! [Pa]  Basal stress, x-direction
        real(wp), intent(OUT) :: taub_acy(:,:)        ! [Pa]  Basal stress, y-direction
        real(wp), intent(IN)  :: dd_ab_3D(:,:,:)      ! [m^2/a] SIA diffusivity, ab-nodes
        real(wp), intent(IN)  :: H_ice(:,:)           ! [m]   Ice thickness
        real(wp), intent(IN)  :: taud_acx(:,:)        ! [Pa]  Driving stress, x-direction 
        real(wp), intent(IN)  :: taud_acy(:,:)        ! [Pa]  Driving stress, y-direction
        real(wp), intent(IN)  :: f_pmp(:,:)           ! [--]  Fraction of grid point at pressure melting point
        real(wp), intent(IN)  :: zeta_aa(:)           ! [--]  Height vector 0:1 
        real(wp), intent(IN)  :: dx                   ! [m]
        real(wp), intent(IN)  :: cf_sia               ! [m/a Pa-1]
        real(wp), intent(IN)  :: rho_ice              ! [kg m-3] Ice density 
        real(wp), intent(IN)  :: g                    ! [m s-2]  Gravity acceleration  

        ! Local variables
        integer :: i, j, k, nx, ny, nz_aa 
        real(wp) :: dd_acx, dd_acy 
        real(wp), allocatable :: dd_ab(:,:)           ! [m^2/a] SIA diffusivity, ab-nodes
        real(wp), allocatable :: f_pmp_ab(:,:)
        real(wp) :: H_ac 
        real(wp), parameter :: cf_sia_frozen = 0.0    ! No sliding for purely frozen points 

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
        taub_acx = 0.0_wp 
        taub_acy = 0.0_wp 

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

        real(wp), intent(OUT) :: ux_b(:,:)            ! [m/a] SIA basal velocity x-direction, acx-nodes
        real(wp), intent(OUT) :: uy_b(:,:)            ! [m/a] SIA basal velocity y-direction, acy-nodes
        real(wp), intent(OUT) :: taub_acx(:,:)        ! [Pa]  Basal stress, x-direction
        real(wp), intent(OUT) :: taub_acy(:,:)        ! [Pa]  Basal stress, y-direction
        real(wp), intent(IN)  :: H_ice(:,:)           ! [m]   Ice thickness
        real(wp), intent(IN)  :: dzsdx(:,:)           ! [m/m] Surface gradient x-direction 
        real(wp), intent(IN)  :: dzsdy(:,:)           ! [m/m] Surface gradient y-direction
        real(wp), intent(IN)  :: f_pmp(:,:)           ! [--]  Fraction of grid point at pressure melting point
        real(wp), intent(IN)  :: zeta_aa(:)           ! [--]  Height vector 0:1 
        real(wp), intent(IN)  :: dx                   ! [m]
        real(wp), intent(IN)  :: cf_sia               ! [m/a Pa-1]
        real(wp), intent(IN)  :: rho_ice              ! [kg m-3] Ice density 
        real(wp), intent(IN)  :: g                    ! [m s-2]  Gravity acceleration  

        ! Local variables
        integer :: i, j, k, nx, ny, nz_aa 
        real(wp) :: dd_acx, dd_acy 
        real(wp), allocatable :: H_ice_ab(:,:) 
        real(wp), allocatable :: slope_ab(:,:) 
        real(wp), allocatable :: dd_ab(:,:)           ! [m^2/a] SIA diffusivity, ab-nodes
        real(wp), allocatable :: f_pmp_ab(:,:)
        real(wp) :: H_ac 
        real(wp), parameter :: cf_sia_frozen = 0.0    ! No sliding for purely frozen points 

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
        taub_acx = 0.0_wp 
        taub_acy = 0.0_wp 

        return
        
    end subroutine calc_velocity_basal_sia_00

    subroutine calc_basal_stress(taub_acx,taub_acy,beta_acx,beta_acy,ux_b,uy_b)
        ! Calculate the basal stress resulting from sliding (friction times velocity)
        ! Note: calculated on ac-nodes.
        ! taub [Pa] 
        ! beta [Pa a m-1]
        ! u    [m a-1]
        ! taub = beta*u (here defined with taub in the same direction as u)

        implicit none 

        real(wp), intent(OUT) :: taub_acx(:,:)   ! [Pa] Basal stress (acx nodes)
        real(wp), intent(OUT) :: taub_acy(:,:)   ! [Pa] Basal stress (acy nodes)
        real(wp), intent(IN)  :: beta_acx(:,:)   ! [Pa a m-1] Basal friction (acx nodes)
        real(wp), intent(IN)  :: beta_acy(:,:)   ! [Pa a m-1] Basal friction (acy nodes)
        real(wp), intent(IN)  :: ux_b(:,:)       ! [m a-1] Basal velocity (acx nodes)
        real(wp), intent(IN)  :: uy_b(:,:)       ! [m a-1] Basal velocity (acy nodes)
        
        real(wp), parameter :: tol = 1e-3_wp 

        ! Calculate basal stress 
        taub_acx = beta_acx * ux_b 
        taub_acy = beta_acy * uy_b 

        ! Avoid underflows
        where(abs(taub_acx) .lt. tol) taub_acx = 0.0_wp 
        where(abs(taub_acy) .lt. tol) taub_acy = 0.0_wp 
        
        return 

    end subroutine calc_basal_stress

end module velocity_sia
