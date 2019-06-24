module thermodynamics 
    ! This module contains some general thermodynamics subroutines
    ! that could be used by many icetemp solver approaches.
    ! Note: once icetemp is working well, this module could be 
    ! remerged into icetemp as one module. 

    use yelmo_defs, only : prec, dp, sec_year, pi, T0, g, rho_ice, rho_sw, rho_w

    implicit none 

    real(prec), parameter :: L_ice = 3.35e5       ! Specific latent heat of fusion of ice [J Kg-1]

    private  

    public :: calc_bmb_grounded
    public :: calc_advec_vertical_column
    public :: calc_advec_horizontal_column
    public :: calc_strain_heating
    public :: calc_strain_heating_sia
    public :: calc_basal_heating
    public :: calc_specific_heat_capacity
    public :: calc_thermal_conductivity
    public :: calc_T_pmp
    public :: calc_f_pmp
    public :: calc_T_base_shlf_approx
    public :: define_temp_linear_3D
    public :: calc_temp_linear_column
    public :: define_temp_robin_3D
    
contains

    elemental subroutine calc_bmb_grounded(bmb_grnd,T_prime_b,dTdz_b,kt_b,rho_ice, &
                                                 Q_b,Q_geo_now,f_grnd)
        ! Calculate everywhere there is at least some grounded ice 
        ! (centered aa node calculation)

        ! Note: calculated bmb_grounded here as if the ice point is fully grounded, 
        ! bmb_grnd and bmb_shlf will then be weighted average using f_grnd externally
        ! (to allow ice topography to evolve with different time steps)

        implicit none 
        
        real(prec), intent(OUT) :: bmb_grnd          ! [m/a ice equiv.] Basal mass balance, grounded
        real(prec), intent(IN)  :: T_prime_b         ! [K] Basal ice temp relative to pressure melting point (ie T_prime_b=0 K == temperate)
        real(prec), intent(IN)  :: dTdz_b            ! [K/m] Gradient of temperature in ice, basal layer
        real(prec), intent(IN)  :: kt_b              ! [J a-1 m-1 K-1] Heat conductivity in ice, basal layer 
        real(prec), intent(IN)  :: rho_ice           ! [kg m-3] Ice density 
        real(prec), intent(IN)  :: Q_b               ! [J a-1 m-2] Basal heat production from friction and strain heating
        real(prec), intent(IN)  :: Q_geo_now         ! [J a-1 m-2] Geothermal heat flux 
        real(prec), intent(IN)  :: f_grnd            ! [--] Grounded fraction (centered aa node)                 
        !real(prec), intent(IN)  :: kt_m              ! [J a-1 m-1 K-1] Heat conductivity in mantle (lithosphere) 

        ! Local variables
        real(prec) :: coeff 
        real(prec), parameter :: tol = 1e-10  

        ! Calculate the grounded basal mass balance following 
        ! Cuffey and Patterson (2010), Eq. 9.38 (Page 420)
        ! with an addition  of the term Q_strn_b 
        ! Note: formula only applies when base is temperate, ie
        ! when f_pmp > 0.0 
        if ( f_grnd .gt. 0.0 .and. T_prime_b .eq. 0.0_prec) then 
            ! Bed is grounded and temperate, calculate basal mass balance  

!                 if (cond_bed) then 
!                     ! Following grisli formulation: 
                    
!                     bmb_grnd = -1.0_prec/(rho_ice*L_ice)* ( Q_b + kt_b*dTdz_b - kt_m*dTrdz_b ) 

!                 else
!                     ! Classic Cuffey and Patterson (2010) formula 

!                     bmb_grnd = -1.0_prec/(rho_ice*L_ice)* ( Q_b + kt_b*dTdz_b + (Q_geo_now) ) 

!                 end if 
            
!             bmb_grnd = -1.0_prec/(rho_ice*L_ice)* ( Q_b + kt_b*dTdz_b - kt_m*dTrdz_b )
            
            ! Classic Cuffey and Patterson (2010) formula
            bmb_grnd = -1.0_prec/(rho_ice*L_ice)* ( Q_b + kt_b*dTdz_b + Q_geo_now )

        else 
            ! No basal mass change possible if bed is not temperate 

            bmb_grnd = 0.0_prec 

        end if 

        ! Limit small values to avoid underflow errors 
        if (abs(bmb_grnd) .lt. tol) bmb_grnd = 0.0_prec 

        return 

    end subroutine calc_bmb_grounded 

    subroutine calc_advec_vertical_column(advecz,Q,uz,H_ice,zeta_aa)
        ! Calculate vertical advection term advecz, which enters
        ! advection equation as
        ! Q_new = Q - dt*advecz = Q - dt*u*dQ/dx
        ! 1st order upwind scheme 

        implicit none 

        real(prec), intent(OUT)   :: advecz(:)      ! nz_aa: bottom, cell centers, top 
        real(prec), intent(INOUT) :: Q(:)           ! nz_aa: bottom, cell centers, top 
        real(prec), intent(IN)    :: uz(:)          ! nz_ac: cell boundaries
        real(prec), intent(IN)    :: H_ice          ! Ice thickness 
        real(prec), intent(IN)    :: zeta_aa(:)    ! nz_aa, cell centers
        
        ! Local variables
        integer :: k, nz_aa   
        real(prec) :: u_aa, dx  

        nz_aa = size(zeta_aa,1)

        advecz = 0.0 

        ! Loop over internal cell centers and perform upwind advection 
        do k = 2, nz_aa-1 
            
            u_aa = 0.5_prec*(uz(k-1)+uz(k))
            
            if (u_aa < 0.0) then 
                ! Upwind negative
                dx = H_ice*(zeta_aa(k+1)-zeta_aa(k))
                advecz(k) = uz(k)*(Q(k+1)-Q(k))/dx   
            else
                ! Upwind positive
                dx = H_ice*(zeta_aa(k)-zeta_aa(k-1))
                advecz(k) = uz(k-1)*(Q(k)-Q(k-1))/dx
            end if 
        end do 

        return 

    end subroutine calc_advec_vertical_column

    subroutine calc_advec_horizontal_column(advecxy,var_ice,H_ice,ux,uy,dx,i,j)
        ! Newly implemented advection algorithms (ajr)
        ! Output: [K a-1]

        ! [m-1] * [m a-1] * [K] = [K a-1]

        implicit none

        real(prec), intent(OUT) :: advecxy(:)       ! nz_aa 
        real(prec), intent(IN)  :: var_ice(:,:,:)   ! nx,ny,nz_aa  Enth, T, age, etc...
        real(prec), intent(IN)  :: H_ice(:,:)       ! nx,ny 
        real(prec), intent(IN)  :: ux(:,:,:)        ! nx,ny,nz
        real(prec), intent(IN)  :: uy(:,:,:)        ! nx,ny,nz
        real(prec), intent(IN)  :: dx  
        integer,    intent(IN)  :: i, j 

        ! Local variables 
        integer :: k, nx, ny, nz_aa 
        real(prec) :: ux_aa, uy_aa 
        real(prec) :: dx_inv, dx_inv2
        real(prec) :: advecx, advecy 

        ! Define some constants 
        dx_inv  = 1.0_prec / dx 
        dx_inv2 = 1.0_prec / (2.0_prec*dx)

        nx  = size(var_ice,1)
        ny  = size(var_ice,2)
        nz_aa = size(var_ice,3) 

        advecx  = 0.0 
        advecy  = 0.0 
        advecxy = 0.0 

        ! Loop over each point in the column
        do k = 1, nz_aa 

            ! Estimate direction of current flow into cell (x and y), centered in vertical layer and grid point
            ux_aa = 0.5_prec*(ux(i,j,k)+ux(i-1,j,k))
            uy_aa = 0.5_prec*(uy(i,j,k)+uy(i,j-1,k))
            
            ! Explicit form (to test different order approximations)
            if (ux_aa .gt. 0.0 .and. i .ge. 3) then  
                ! Flow to the right 

                ! 1st order
                !advecx = dx_inv * ux(i-1,j,k)*(-(var_ice(i-1,j,k)-var_ice(i,j,k)))
                ! 2nd order
                advecx = dx_inv2 * ux(i-1,j,k)*(-(4.0*var_ice(i-1,j,k)-var_ice(i-2,j,k)-3.0*var_ice(i,j,k)))

            else if (ux_aa .lt. 0.0 .and. i .le. nx-2) then 
                ! Flow to the left

                ! 1st order 
                !advecx = dx_inv * ux(i,j,k)*((var_ice(i+1,j,k)-var_ice(i,j,k)))
                ! 2nd order
                advecx = dx_inv2 * ux(i,j,k)*((4.0*var_ice(i+1,j,k)-var_ice(i+2,j,k)-3.0*var_ice(i,j,k)))

            else 
                ! No flow 
                advecx = 0.0

            end if 

            if (uy_aa .gt. 0.0 .and. j .ge. 3) then   
                ! Flow to the right 

                ! 1st order
                !advecy = dx_inv * uy(i,j-1,k)*(-(var_ice(i,j-1,k)-var_ice(i,j,k)))
                ! 2nd order
                advecy = dx_inv2 * uy(i,j-1,k)*(-(4.0*var_ice(i,j-1,k)-var_ice(i,j-2,k)-3.0*var_ice(i,j,k)))

            else if (uy_aa .lt. 0.0 .and. j .le. ny-2) then 
                ! Flow to the left

                ! 1st order
                !advecy = dx_inv * uy(i,j,k)*((var_ice(i,j+1,k)-var_ice(i,j,k)))
                ! 2nd order
                advecy = dx_inv2 * uy(i,j,k)*((4.0*var_ice(i,j+1,k)-var_ice(i,j+2,k)-3.0*var_ice(i,j,k)))
                
            else
                ! No flow 
                advecy = 0.0 

            end if 
                    
            ! Combine advection terms for total contribution 
            advecxy(k) = (advecx+advecy)

        end do 

        return 

    end subroutine calc_advec_horizontal_column
    
    subroutine calc_strain_heating(Q_strn,de,visc,cp,rho_ice)
        ! Calculate the general 3D internal strain heating
        ! as sum(D_ij*tau_ij)  (strain*stress)
        ! where stress has been calculated as stress_ij = 2*visc*strain
        ! Units: Q_Strn [J a-1 m-3]

        ! Note: we use the simpler approach because in the shallow
        ! model, the stress rate is simply the strain rate squared

        ! Directly from Q_strn = tr(stress*strain)
        ! (Cuffey and Patterson (2010) pag. 417, eq. 9.30)

!         Q_strn = ( strss%txx*strn%dxx &
!                  + strss%tyy*strn%dyy &
!                  + strss%tzz*strn%dzz &    ! this term is not available yet in the code, would need to be calculated
!              + 2.0*strss%txy*strn%dxy &
!              + 2.0*strss%txz*strn%dxz &
!              + 2.0*strss%tyz*strn%dyz )

        ! Simpler approach:
        ! Calculate Q_strn from effective strain rate and viscosity
        ! (Greve and Blatter (2009) eqs. 4.7 and 5.65): 
        !     Q_strn = tr(stress*strn) = tr(2*visc*strn*strn) = 2*visc*tr(strn*strn) = 4*visc*de^2
        !     with tr(strn*strn) = 2*de^2

        implicit none

        real(prec),            intent(OUT) :: Q_strn(:,:,:)      ! nx,ny,nz_aa [K a-1] Heat production
        real(prec),            intent(IN)  :: de(:,:,:)          ! nx,ny,nz_aa [a-1] Effective strain rate 
        real(prec),            intent(IN)  :: visc(:,:,:)        ! nx,ny,nz_aa [Pa a-1] Viscosity
        real(prec),            intent(IN)  :: cp(:,:,:)          ! nx,ny,nz_aa [J kg-1 K-1] Specific heat capacity
        real(prec),            intent(IN)  :: rho_ice            ! [kg m-3] Ice density 

        ! Local variables
        integer :: nz_aa 
        !real(prec), parameter :: Q_strn_max = 0.1          ! Q_strn > 0.1 is already high, it's only a safety valve.

        nz_aa = size(Q_strn,3)

        ! Calculate strain heating 
        ! following Greve and Blatter (2009), Eqs. 4.7 and 5.65
        Q_strn = 4.0*visc * de**2
        
        ! Limit strain heating to reasonable values 
        !where (Q_strn .gt. Q_strn_max) Q_strn = Q_strn_max

        return 

    end subroutine calc_strain_heating
    
    subroutine calc_strain_heating_sia(Q_strn,ux,uy,dzsdx,dzsdy,cp,H_ice,rho_ice,zeta_aa,zeta_ac)

        ! Calculate the general 3D internal strain heating
        ! as sum(D_ij*tau_ij)  (strain*stress)
        ! where stress has been calculated as stress_ij = 2*visc*strain
        ! Units: Q_strn = Q [J a-1 m-3]

        ! SIA approximation:
        ! Q_strn = rho*g*H*(duxdz*dzsdx + duydz*dzsdy)
        ! Units: [J a-1 m-3]

        implicit none

        real(prec),            intent(OUT) :: Q_strn(:,:,:)      ! nx,ny,nz_aa  [Pa m a-1 ??] Heat production
        real(prec),            intent(IN)  :: ux(:,:,:)          ! nx,ny,nz_aa  [m a-1] Velocity x-direction
        real(prec),            intent(IN)  :: uy(:,:,:)          ! nx,ny,nz_aa  [m a-1] Velocity y-direction
        real(prec),            intent(IN)  :: dzsdx(:,:)         ! nx,ny        [m m-1] Surface slope x-direction
        real(prec),            intent(IN)  :: dzsdy(:,:)         ! nx,ny        [m m-1] Surface slope y-direction
        real(prec),            intent(IN)  :: cp(:,:,:)          ! nx,ny,nz_aa  [J/kg/K] Specific heat capacity
        real(prec),            intent(IN)  :: H_ice(:,:)         ! nx,ny        [m] Ice thickness
        real(prec),            intent(IN)  :: rho_ice            ! [kg m-3] Ice density 
        real(prec),            intent(IN)  :: zeta_aa(:)        ! [-] Height axis, centered aa-nodes 
        real(prec),            intent(IN)  :: zeta_ac(:)        ! [-] Height axis, boundaries ac-nodes

        ! Local variables
        integer :: i, j, k, nx, ny, nz_aa 
        real(prec) :: ux_aa_up, ux_aa_dwn
        real(prec) :: uy_aa_up, uy_aa_dwn
        real(prec) :: duxdz, duydz
        real(prec) :: dzsdx_aa, dzsdy_aa  
        real(prec) :: dz, depth 

        real(prec), parameter :: Q_strn_max = 0.1          ! Q_strn > 0.1 is already high, it's only a safety valve. 
        
        nx    = size(Q_strn,1)
        ny    = size(Q_strn,2)
        nz_aa = size(Q_strn,3)

        ! Base and surface Q_strn == 0.0
        Q_strn = 0.0 

        do j = 2, ny 
        do i = 2, nx 

            if (H_ice(i,j) .gt. 0.0) then 
                ! Only calculate strain rate where ice exists on horizontal aa-nodes 

                do k = 2, nz_aa-1 
                    ! Loop over the vertical column of ice 

                    ux_aa_up  = 0.25*(ux(i-1,j,k)+ux(i,j,k)+ux(i-1,j,k+1)+ux(i,j,k+1))
                    ux_aa_dwn = 0.25*(ux(i-1,j,k)+ux(i,j,k)+ux(i-1,j,k-1)+ux(i,j,k-1))
                    
                    uy_aa_up  = 0.25*(uy(i,j-1,k)+uy(i,j,k)+uy(i,j-1,k+1)+uy(i,j,k+1))
                    uy_aa_dwn = 0.25*(uy(i,j-1,k)+uy(i,j,k)+uy(i,j-1,k-1)+uy(i,j,k-1))
                    
                    dz = H_ice(i,j)*(zeta_ac(k) - zeta_ac(k-1)) 

                    duxdz = (ux_aa_up-ux_aa_dwn)/dz 
                    duydz = (uy_aa_up-uy_aa_dwn)/dz 

                    dzsdx_aa = 0.5*(dzsdx(i-1,j)+dzsdx(i,j))
                    dzsdy_aa = 0.5*(dzsdy(i,j-1)+dzsdy(i,j))
                    
                    depth = H_ice(i,j)*(1.0-zeta_aa(k))
                    
                    Q_strn(i,j,k) = (-rho_ice*g*depth) * (duxdz*dzsdx_aa + duydz*dzsdy_aa)
                    
                end do 

            end if 

        end do 
        end do 

        ! [m s-1] [m] 
        ! Limit strain heating to reasonable values 
        where (Q_strn .gt. Q_strn_max) Q_strn = Q_strn_max

        return 

    end subroutine calc_strain_heating_sia

    subroutine calc_basal_heating(Q_b,ux_b,uy_b,taub_acx,taub_acy)
         ! Qb [J a-1 m-2] == [m a-1] * [J m-3]
         ! Note: grounded ice fraction f_grnd_acx/y not used here, because taub_acx/y already accounts
         ! for the grounded fraction via beta_acx/y: Q_b = tau_b*u = -beta*u*u.

        real(prec), intent(OUT) :: Q_b(:,:)               ! [J a-1 K-1] Basal heat production (friction)
        real(prec), intent(IN)  :: ux_b(:,:)              ! Basal velocity, x-component (staggered x)
        real(prec), intent(IN)  :: uy_b(:,:)              ! Basal velocity, y-compenent (staggered y)
        real(prec), intent(IN)  :: taub_acx(:,:)          ! Basal friction (staggered x)
        real(prec), intent(IN)  :: taub_acy(:,:)          ! Basal friction (staggered y)
        
        ! Local variables
        integer    :: i, j, nx, ny 
        real(prec), allocatable :: Qb_acx(:,:)
        real(prec), allocatable :: Qb_acy(:,:)

        nx = size(Q_b,1)
        ny = size(Q_b,2)

        allocate(Qb_acx(nx,ny))
        allocate(Qb_acy(nx,ny))

        ! Determine basal frictional heating values (staggered acx/acy nodes)
        Qb_acx = abs(ux_b*taub_acx)   ! [Pa m a-1] == [J a-1 m-2]
        Qb_acy = abs(uy_b*taub_acy)   ! [Pa m a-1] == [J a-1 m-2]

        Q_b = 0.0  
 
        ! Get basal frictional heating on centered nodes (aa-grid)          
        do j = 2, ny
        do i = 2, nx

             ! Average from ac-nodes to aa-node
             Q_b(i,j) = 0.25*(Qb_acx(i,j)+Qb_acx(i-1,j)+Qb_acy(i,j)+Qb_acy(i,j-1))
 
        end do 
        end do 
 
        return 
 
    end subroutine calc_basal_heating

    subroutine calc_basal_heating1(Q_b,ux_b,uy_b,taub_acx,taub_acy)
        ! Q_b [J a-1 m-2] == [m a-1] * [J m-3]
        ! Note: grounded ice fraction f_grnd_acx/y not used here, because taub_acx/y already accounts
        ! for the grounded fraction via beta_acx/y: Q_b = tau_b*u = -beta*u*u.
        ! Note: it is assumed that strain heating is zero at the very base (no thickness) layer, so
        ! it is not included here. 

        implicit none 

        real(prec), intent(OUT) :: Q_b(:,:)               ! [J a-1 K-1] Basal heat production (friction)
        real(prec), intent(IN)  :: ux_b(:,:)              ! Basal velocity, x-component (staggered x)
        real(prec), intent(IN)  :: uy_b(:,:)              ! Basal velocity, y-compenent (staggered y)
        real(prec), intent(IN)  :: taub_acx(:,:)          ! Basal friction (staggered x)
        real(prec), intent(IN)  :: taub_acy(:,:)          ! Basal friction (staggered y)
        
        ! Local variables
        integer    :: i, j, nx, ny 
        integer    :: im1, jm1 
        real(prec) :: Qb_acx_1, Qb_acx_2, Qb_acy_1, Qb_acy_2   

        nx = size(Q_b,1)
        ny = size(Q_b,2)

        ! Initially set basal heating to zero everywhere 
        Q_b = 0.0  

        ! Get basal frictional heating on centered nodes (aa-grid)
        do j = 1, ny
        do i = 1, nx

            im1 = max(i-1,1)
            jm1 = max(j-1,1)

            ! Determine basal frictional heating values (staggered acx/acy nodes)
            ! [Pa m a-1] == [J a-1 m-2]
            Qb_acx_1 = abs(ux_b(im1,j)*taub_acx(im1,j))
            Qb_acx_2 = abs(ux_b(i,j)  *taub_acx(i,j))
            Qb_acy_1 = abs(uy_b(i,jm1)*taub_acy(i,jm1))
            Qb_acy_2 = abs(uy_b(i,j)  *taub_acy(i,j))
            
            ! Average from ac-nodes to aa-node
            Q_b(i,j) = 0.25*(Qb_acx_1+Qb_acx_2+Qb_acy_1+Qb_acy_2)

        end do
        end do

        return 

    end subroutine calc_basal_heating1

    elemental function calc_specific_heat_capacity(T_ice) result(cp)

        implicit none 

        real(prec), intent(IN) :: T_ice  
        real(prec) :: cp 

        ! Specific heat capacity (Greve and Blatter, 2009, Eq. 4.39; Ritz, 1987)
        cp = (146.3 +7.253*T_ice)    ! [J kg-1 K-1]

        return 

    end function calc_specific_heat_capacity

    elemental function calc_thermal_conductivity(T_ice) result(ct)

        implicit none 

        real(prec), intent(IN) :: T_ice  
        real(prec) :: ct 

        ! Heat conductivity (Greve and Blatter, 2009, Eq. 4.37; Ritz, 1987)
        ct = 9.828*exp(-0.0057*T_ice)*sec_year  ! [W m-1 K-1 * sec_year] => [J m-1 K-1 a-1]

        return 

    end function calc_thermal_conductivity
    
    elemental function calc_T_pmp(H_ice,zeta,T0) result(T_pmp)
        ! Greve and Blatter (Chpt 4, pg 54), Eq. 4.13
        ! This gives the pressure-corrected melting point of ice
        ! where H_ice*(1-zeta) is the thickness of ice overlying the current point 
        
        implicit none 

        real(prec), intent(IN) :: H_ice  ! [m] Total ice thickness of this point
        real(prec), intent(IN) :: zeta   ! [-] Fractional height of this point within the ice
        real(prec), intent(IN) :: T0     ! [K] Reference freezing point of water (273.15 K)
        real(prec) :: T_pmp              ! [K] Pressure corrected melting point

        ! Local variables
        real(prec) :: depth   
        real(prec), parameter :: beta1 = 8.74e-4    ! [K m^-1]   beta1 = (beta*rho*g), beta=9.8e-8 [K Pa^-1]
        !real(prec), parameter :: beta1 = 8.66e-4    ! [K m^-1]   EISMINT2 value
        
        ! Get thickness of ice above current point
        depth = H_ice*(1.0-zeta)

        ! Calculate the pressure-corrected melting point
        T_pmp = T0 - beta1*depth

        ! ajr: note: should we account here for whether ice is floating or not, changing the pressure? 
        
        return 

    end function calc_T_pmp 

    elemental function calc_f_pmp(T_ice,T_pmp,gamma,f_grnd) result(f_pmp)
        ! Calculate the fraction of gridpoint at the pressure melting point (pmp),
        ! ie, when T_ice >= T_pmp. Facilitates a smooth transition between
        ! frozen and temperate ice. (Greve, 2005; Hindmarsh and Le Meur, 2001)

        implicit none 

        real(prec), intent(IN) :: T_ice
        real(prec), intent(IN) :: T_pmp
        real(prec), intent(IN) :: gamma
        real(prec), intent(IN) :: f_grnd  
        real(prec) :: f_pmp 

        if (f_grnd .eq. 0.0) then
            ! Floating points are temperate by default
            f_pmp = 1.0 

        else 
            ! Calculate the fraction at the pressure melting point 

            if (gamma .eq. 0.0) then
                ! No decay function, binary pmp fraction

                if (T_ice .ge. T_pmp) then 
                    f_pmp = 1.0
                else 
                    f_pmp = 0.0 
                end if 

            else

                ! Apply decay function 
                f_pmp = min(1.0, exp((T_ice-T_pmp)/gamma) )

                ! Ensure pure values of 0.0 and 1.0 beyond a threshold 
                if (f_pmp .lt. 1e-2)        f_pmp = 0.0 
                if (f_pmp .gt. (1.0-1e-2))  f_pmp = 1.0 

            end if 

        end if 

        return 

    end function calc_f_pmp 
    
    elemental function calc_T_base_shlf_approx(H_ice,T_pmp,H_grnd) result(T_base_shlf)
        ! Calculate the basal shelf temperature for floating ice
        ! as the estimated freezing temperature of seawater
        ! following Jenkins (1991)
        ! ajr: modified to ensure that temp approaches T_pmp as ice becomes grounded 

        implicit none 

        real(prec), intent(IN) :: H_ice 
        real(prec), intent(IN) :: T_pmp 
        real(prec), intent(IN) :: H_grnd 
        real(prec) :: T_base_shlf

        ! Local variables 
        real(prec), parameter :: a1 = - 0.0575      ! [degC / PSU]
        real(prec), parameter :: b1 =   0.0901      ! [degC]
        real(prec), parameter :: c1 =   7.61E-4     ! [degC / m]
        real(prec), parameter :: S0 =   34.75       ! [g / kg == PSU]
        real(prec) :: f_scalar, H_grnd_lim 

        T_base_shlf = a1*S0 + b1 + c1*(rho_ice/rho_sw)*H_ice + T0 

        ! Additionally ensure that the shelf temperature is approaching the pressure melting point
        ! as the grounding line is reached 
        H_grnd_lim = -100.0
        f_scalar = (H_grnd - H_grnd_lim) / H_grnd_lim
        f_scalar = min(f_scalar,1.0)
        f_scalar = max(f_scalar,0.0)
        T_base_shlf   = f_scalar*T_base_shlf + (1.0-f_scalar)*T_pmp

        ! Limit everything to the ice pressure melting point 
        T_base_shlf = min(T_base_shlf,T_pmp)

        return 

    end function calc_T_base_shlf_approx
    
    subroutine define_temp_linear_3D(T_ice,zeta_aa,H_ice,T_srf)
        ! Define a linear vertical temperature profile
       
        implicit none

        real(prec), intent(INOUT) :: T_ice(:,:,:)       ! nx,ny,nz_aa
        real(prec), intent(IN)    :: zeta_aa(:)
        real(prec), intent(IN)    :: H_ice(:,:)
        real(prec), intent(IN)    :: T_srf(:,:) 
           
        ! Local variables
        integer :: i, j, k, nx, ny, nz_aa   
        real(prec) :: T_base, T_pmp  

        nx    = size(T_ice,1)
        ny    = size(T_ice,2)
        nz_aa = size(T_ice,3)

        do j = 1, ny 
        do i = 1, nx

            if (H_ice(i,j) .gt. 0.0) then
                ! Ice is present, define linear temperature profile with frozen bed (-10 degC)
                 
                T_base       = calc_T_pmp(H_ice(i,j),zeta_aa(1),T0) - 10.0 
                T_ice(i,j,:) = calc_temp_linear_column(T_srf(i,j),T_base,T0,zeta_aa)

            else 
                ! No ice present, set equal to T_pmp (ie, T0)
                T_ice(i,j,:) = T0

            end if 

            ! Calculate enthalpy as well for all layers 
            ! TO DO !
            !do k = 1, nz_aa
            !    T_pmp = calc_T_pmp(H_ice(i,j),zeta_aa(k),T0)
            !    enth_ice(i,j,k) = enth_fct_temp_omega(T_ice(i,j,k)-T_pmp, 0.0_prec)   ! Assume zero water content
            !end do 

        end do 
        end do  
        
        return 

    end subroutine define_temp_linear_3D 

    function calc_temp_linear_column(T_srf,T_base,T0,zeta_aa) result(T_ice)

        implicit none 

        real(prec), intent(IN) :: T_srf
        real(prec), intent(IN) :: T_base
        real(prec), intent(IN) :: T0 
        real(prec), intent(IN) :: zeta_aa(:) 
        real(prec) :: T_ice(size(zeta_aa,1))

        ! Local variables 
        integer :: k 
        real(prec) :: T_srf_now 

        ! Limit surface temperature to below melting 
        T_srf_now = min(T_srf,T0) 

        ! Perform linear interpolation 
        do k = 1, size(zeta_aa,1)
            T_ice(k) = T_base  + zeta_aa(k)*(T_srf_now-T_base)
        end do 

        return 

    end function calc_temp_linear_column

    subroutine define_temp_robin_3D (T_ice,T_pmp,cp,ct,Q_geo,T_srf,H_ice,H_w,smb,bmb,f_grnd,zeta_aa,cold)
        ! Robin solution for thermodynamics for a given column of ice 
        ! Note zeta=height, k=1 base, k=nz surface 

        implicit none 

        real(prec), intent(OUT) :: T_ice(:,:,:)     ! [K] Ice column temperature
        real(prec), intent(IN)  :: T_pmp(:,:,:)     ! [K] Pressure melting point temp.
        real(prec), intent(IN)  :: cp(:,:,:)        ! [J kg-1 K-1] Specific heat capacity
        real(prec), intent(IN)  :: ct(:,:,:)        ! [J a-1 m-1 K-1] Heat conductivity 
        real(prec), intent(IN)  :: Q_geo(:,:)       ! [mW m-2] Geothermal heat flux 
        real(prec), intent(IN)  :: T_srf(:,:)       ! [K] Surface temperature 
        real(prec), intent(IN)  :: H_ice(:,:)       ! [m] Ice thickness 
        real(prec), intent(IN)  :: H_w(:,:)         ! [m] Basal water layer thickness 
        real(prec), intent(IN)  :: smb(:,:)         ! [m a-1] Surface mass balance (melting is negative)
        real(prec), intent(IN)  :: bmb(:,:)         ! [m a-1] Basal mass balance (melting is negative)
        real(prec), intent(IN)  :: f_grnd(:,:)      ! [--] Floating point or grounded?
        real(prec), intent(IN)  :: zeta_aa(:)       ! [--] Vertical zeta coordinates (zeta==height), aa-nodes
        logical,    intent(IN)  :: cold             ! False: robin as normal, True: ensure cold base 

        ! Local variable
        integer :: i, j, k, nx, ny, nz_aa  
        logical :: is_float 
        
        real(prec), allocatable :: T1(:) 

        nx    = size(T_ice,1)
        ny    = size(T_ice,2)
        nz_aa = size(zeta_aa)

        allocate(T1(nz_aa))

        do j = 1, ny
        do i = 1, nx 

            is_float = (f_grnd(i,j) .eq. 0.0)

            T_ice(i,j,:) = calc_temp_robin_column(zeta_aa,T_pmp(i,j,:),ct(i,j,:),cp(i,j,:),rho_ice,H_ice(i,j), &
                                                  T_srf(i,j),smb(i,j)+bmb(i,j),Q_geo(i,j),is_float)

            if (cold) then 
                T1(nz_aa) = T_srf(i,j)
                T1(1)     = T_pmp(i,j,1) - 10.0

                ! Intermediate layers are linearly interpolated 
                do k = 2, nz_aa-1 
                    T1(k) = T1(1)+zeta_aa(k)*(T1(nz_aa)-T1(1))
                end do 
                
                ! Average Robin solution with cold linear solution 
                T_ice(i,j,:) = 0.5_prec*(T_ice(i,j,:) + T1)
            end if 

        end do 
        end do 

        return 

    end subroutine define_temp_robin_3D

    function calc_temp_robin_column(zeta_aa,T_pmp,kt,cp,rho_ice,H_ice,T_srf,mb_net,Q_geo,is_float) result(T_ice)
        ! This function will impose a temperature solution in a given ice column.
        ! For:
        !  Grounded ice with positive net mass balance: Robin solution where possible
        !  Grounded ice with negative net mass balance: Linear profile 
        !  Floating ice: Linear profile 
        !  No or thin ice: Surface temperature 
        
        implicit none 

        real(prec), intent(IN) :: zeta_aa(:) 
        real(prec), intent(IN) :: T_pmp(:)
        real(prec), intent(IN) :: kt(:) 
        real(prec), intent(IN) :: cp(:) 
        real(prec), intent(IN) :: rho_ice 
        real(prec), intent(IN) :: H_ice 
        real(prec), intent(IN) :: T_srf 
        real(prec), intent(IN) :: mb_net
        real(prec), intent(IN) :: Q_geo 
        logical,    intent(IN) :: is_float 

        real(prec) :: T_ice(size(zeta_aa,1))

        ! Local variables 
        integer    :: k, nz_aa 
        real(prec) :: dTdz_b, z, kappa, ll   

        real(prec), parameter :: sqrt_pi   = sqrt(pi) 
        real(prec), parameter :: T_ocn     = 271.15   ! [K]
        real(prec), parameter :: H_ice_min = 0.1      ! [m] Minimum ice thickness to calculate Robin solution 
        real(prec), parameter :: mb_net_min = 1e-2    ! [m a-1] Minimum allowed net mass balance for stability
        real(prec) :: Q_geo_now, mb_now  

        nz_aa = size(T_ice,1) 
        
        Q_geo_now = Q_geo *1e-3*sec_year    ! [mW m-2] => [J a-1 m-2]

        ! Calculate temperature gradient at base 
        dTdz_b = -Q_geo_now/kt(1) 

        if (.not. is_float .and. H_ice .gt. H_ice_min .and. mb_net .gt. 0.0) then 
            ! Impose Robin solution 
            
            !mb_now = max(mb_net,mb_net_min)
            mb_now = mb_net 

            do k = 1, nz_aa 
                z     = zeta_aa(k)*H_ice              ! Ice thickness up to this layer from base 
                kappa = kt(k)/(cp(k)*rho_ice)       ! Thermal diffusivity 
                ll    = sqrt(2*kappa*H_ice/mb_now)  ! Thermal_length_scale

                ! Calculate ice temperature for this layer 
                T_ice(k) = (sqrt_pi/2.0)*ll*dTdz_b*(error_function(z/ll)-error_function(H_ice/ll)) + T_srf 
            end do 

        else if (.not. is_float .and. H_ice .gt. H_ice_min) then 
            ! Impose linear profile with temperate base

            T_ice(nz_aa) = T_srf
            T_ice(1)  = T_pmp(1)

            ! Intermediate layers are linearly interpolated 
            do k = 2, nz_aa-1 
                T_ice(k) = T_ice(1)+zeta_aa(k)*(T_ice(nz_aa)-T_ice(1))
            end do 
            
        else if (is_float) then 
            ! Floating - impose linear profile between sea and surface temps 

            T_ice(nz_aa) = T_srf
            T_ice(1)     = T_ocn 

            ! Intermediate layers are linearly interpolated 
            do k = 2, nz_aa-1 
                T_ice(k) = T_ice(1)+zeta_aa(k)*(T_ice(nz_aa)-T_ice(1))
            end do 
            
        else 
            ! Ice thickness too small 
            T_ice = T_pmp 

        end if 

        ! Ensure ice temperature everywhere is not higher than T_pmp 
        where (T_ice .gt. T_pmp) T_ice = T_pmp 

        return 

    end function calc_temp_robin_column

    function error_function(X) result(ERR)
        ! Purpose: Compute error function erf(x)
        ! Input:   x   --- Argument of erf(x)
        ! Output:  ERR --- erf(x)
        
        implicit none 

        real(prec), intent(IN)  :: X
        real(prec) :: ERR
        
        ! Local variables:
        real(prec)              :: EPS
        real(prec)              :: X2
        real(prec)              :: ER
        real(prec)              :: R
        real(prec)              :: C0
        integer                 :: k
        
        EPS = 1.0e-15
        X2  = X * X
        if (abs(X) < 3.5) then
            ER = 1.0
            R  = 1.0
            do k = 1, 50
                R  = R * X2 / (real(k, prec) + 0.5)
                ER = ER+R
                if(abs(R) < abs(ER) * EPS) then
                    C0  = 2.0 / sqrt(pi) * X * exp(-X2)
                    ERR = C0 * ER
                    EXIT
                end if
            end do
        else
            ER = 1.0
            R  = 1.0
            do k = 1, 12
                R  = -R * (real(k, prec) - 0.5) / X2
                ER = ER + R
                C0  = EXP(-X2) / (abs(X) * sqrt(pi))
                ERR = 1.0 - C0 * ER
                if(X < 0.0) ERR = -ERR
            end do
        end if

        return

    end function error_function

end module thermodynamics 
