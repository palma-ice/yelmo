

module deformation 
    ! This module contains functions related to deformation calculations:
    ! Arrenius function
    ! Flow law
    ! viscosity
    ! strain rate??

    ! Note: 3D arrays defined such that first index (k=1) == base, and max index (k=nk) == surface 

    use yelmo_defs,  only : sp, dp, wp, prec, TOL_UNDERFLOW, T0, rho_ice, g, &
                        strain_2D_class, strain_3D_class, stress_2D_class, stress_3D_class
    use yelmo_tools, only : calc_vertical_integrated_2D, integrate_trapezoid1D_1D, integrate_trapezoid1D_pt, &
                    stagger_nodes_aa_ab_ice, stagger_nodes_acx_ab_ice, stagger_nodes_acy_ab_ice, &
                    staggerdiff_nodes_acx_ab_ice, staggerdiff_nodes_acy_ab_ice, &
                    staggerdiffcross_nodes_acx_ab_ice, staggerdiffcross_nodes_acy_ab_ice, &
                    staggerdiff_nodes_acz_dx_ab_ice, staggerdiff_nodes_acz_dy_ab_ice
                    

    implicit none 
    
    private
    public :: modify_enhancement_factor_bnd 
    public :: define_enhancement_factor_3D
    public :: define_enhancement_factor_2D
    public :: calc_viscosity_glen
    public :: calc_viscosity_glen_2D
    public :: calc_visc_int
    public :: calc_rate_factor
    public :: calc_rate_factor_eismint
    public :: scale_rate_factor_water
    public :: calc_rate_factor_integrated
    public :: calc_strain_rate_tensor
    public :: calc_strain_rate_tensor_aa
    public :: calc_strain_rate_tensor_2D
    public :: calc_stress_tensor 
    public :: calc_stress_tensor_2D
    public :: calc_stress_eigen_values
    public :: strain_2D_alloc 
    public :: stress_2D_alloc

contains 

    subroutine modify_enhancement_factor_bnd(enh,f_grnd,uxy_bar,enh_stream,enh_shlf,umin,umax)
        ! enh field is initially obtained from tracer evolution,
        ! here it is updated to account for streaming and floating regimes 

        implicit none 

        real(prec), intent(INOUT) :: enh(:,:,:)         ! [--] Enhancement factor field
        real(prec), intent(IN)    :: f_grnd(:,:)        ! [--] Fraction of cell grounded
        real(prec), intent(IN)    :: uxy_bar(:,:)       ! [m/a] Depth-averaged velocity magnitude 
        real(prec), intent(IN)    :: enh_stream         ! [--] Enhancement factor for stream regions (SSA grounded)
        real(prec), intent(IN)    :: enh_shlf           ! [--] Enhancement factor for ice shelf regions (SSA floating)
        real(prec), intent(IN)    :: umin               ! [m/a] Minimum transition velocity 
        real(prec), intent(IN)    :: umax               ! [m/a] Maximum transition velocity 

        ! Local variables 
        integer    :: i, j, k, nx, ny, nz 
        real(prec) :: f_mix 

        nx = size(enh,1)
        ny = size(enh,2)
        nz = size(enh,3) 

        if (umax-umin .eq. 0.0_prec) then 
            write(*,*) "modify_enhancement_factor_bnd:: Error: umax cannot equal umin:"
            write(*,*) "umin = ", umin 
            write(*,*) "umax = ", umax 
            stop 
        end if 
            
        do j = 1, ny 
        do i = 1, nx 

            if (f_grnd(i,j) .eq. 0.0_prec) then 
                ! Floating ice, prescribe enh_shlf in column 

                enh(i,j,:) = enh_shlf 

            else 
                ! Grounded ice, determine mixing between enh_bnd for slow
                ! (ie, purely shearing) ice and fast-flowing streaming ice 

                ! Determine mixing ratio (f_mix==1 => streaming ice, f_mix==0 => paleo shearing ice)
                f_mix = (uxy_bar(i,j)-umin) / (umax-umin)
                f_mix = min(f_mix,1.0)
                f_mix = max(f_mix,0.0)

                enh(i,j,:) = f_mix*enh_stream + (1.0-f_mix)*enh(i,j,:) 

            end if 

        end do 
        end do 

        return 

    end subroutine modify_enhancement_factor_bnd

    function define_enhancement_factor_3D(f_shear,f_grnd,uxy_srf,enh_shear,enh_stream,enh_shlf) result(enh)
        ! Greve and Blatter (2009): Chapter 4, page 54 

        implicit none 

        real(prec), intent(IN) :: f_shear(:,:,:)      ! [--] Fraction of cell in shear (as opposed to longitudinal stress)
        real(prec), intent(IN) :: f_grnd(:,:)         ! [--] Fraction of cell grounded
        real(prec), intent(IN) :: uxy_srf(:,:)        ! [m/a] Surface velocity magnitude 
        real(prec), intent(IN) :: enh_shear           ! [--] Enhancement factor for shearing regions (SIA grounded)
        real(prec), intent(IN) :: enh_stream          ! [--] Enhancement factor for stream regions (SSA grounded)
        real(prec), intent(IN) :: enh_shlf            ! [--] Enhancement factor for ice shelf regions (SSA floating)
        real(prec) :: enh(size(f_shear,1),size(f_shear,2),size(f_shear,3))  ! [--] 
        
        ! Local variables
        integer :: k, nx, ny, nz_aa 
        real(prec), allocatable :: enh_ssa_tmp(:,:)
        real(prec) :: f_tmp 
        real(prec), parameter :: uxy_srf_lim = 10.0 

        nx    = size(f_shear,1)
        ny    = size(f_shear,2)
        nz_aa = size(f_shear,3)

        allocate(enh_ssa_tmp(nx,ny))

        ! First calculate the actual ssa enh factor based on f_grnd
        enh_ssa_tmp = f_grnd*enh_stream + (1.0-f_grnd)*enh_shlf
        
        ! Then mix ssa and sia (shear) inland
        ! Note that f_shear should be zero for shelves, so there enh=enh_shlf 
        
        do k = 1, nz_aa 
            enh(:,:,k) = f_shear(:,:,k)*enh_shear   + (1.0-f_shear(:,:,k))*enh_ssa_tmp
        end do 
        
        return 

    end function define_enhancement_factor_3D
    
    elemental function define_enhancement_factor_2D(f_grnd,f_shear,uxy_srf,enh_shear,enh_stream,enh_shlf) result(enh)
        ! Greve and Blatter (2009): Chapter 4, page 54 

        implicit none 

        real(prec), intent(IN) :: f_grnd        ! [--] Fraction of cell grounded
        real(prec), intent(IN) :: f_shear       ! [--] Fraction of cell in shear (as opposed to longitudinal stress)
        real(prec), intent(IN) :: uxy_srf       ! [m/a] Surface velocity magnitude 
        real(prec), intent(IN) :: enh_shear     ! [--] Enhancement factor for shearing regions (SIA grounded)
        real(prec), intent(IN) :: enh_stream    ! [--] Enhancement factor for stream regions (SSA grounded)
        real(prec), intent(IN) :: enh_shlf      ! [--] Enhancement factor for ice shelf regions (SSA floating)
        real(prec) :: enh                       ! [--] 
        
        ! Local variables
        real(prec) :: enh_ssa_tmp, f_tmp 
        real(prec), parameter :: uxy_srf_lim = 10.0 

        ! First calculate an ssa enh factor based on f_grnd, then use this factor 
        ! to further mix with sia (shear) inland.
        ! Note that f_shear should be zero for shelves, so there enh=enh_shlf 

        enh_ssa_tmp = f_grnd*enh_stream + (1.0-f_grnd)*enh_shlf
        enh         = f_shear*enh_shear + (1.0-f_shear)*enh_ssa_tmp
        
        return 

    end function define_enhancement_factor_2D

    function calc_viscosity_glen(de,ATT,H_ice,f_ice,n_glen,visc_min,eps_0) result(visc)
        ! Calculate viscosity based on Glen's flow law 
        ! ATT [a^-1 Pa^-n] is the "depth dependent ice stiffness parameter based on
        !     vertical variations in temperature, chemistry and crystal fabric" (MacAyeal, 1989, JGR)
        ! de [a^-1] is the second-invariant of the strain rate tensor 
        ! visc [Pa a] is the 3D, temperature dependent viscosity field 

        ! Equation: visc = 0.5 * ATT^(-1/n_glen) * (de)^((1-n_glen)/n_glen)
        ! ATT  => from Greve and Blatter (2009), Eq. 4.15 (written as `A(T_prime)`)
        ! de   => from Greve and Blatter (2009), Eq. 6.53
        ! visc => from Greve and Blatter (2009), Eq. 4.22 

        implicit none
        
        real(wp), intent(IN)  :: de(:,:,:)          ! [a^-1] second-invariant of the strain rate tensor
        real(wp), intent(IN)  :: ATT(:,:,:)         ! [a^-1 Pa^-3] Rate factor 
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: n_glen             ! Glen's flow law exponent
        real(wp), intent(IN)  :: visc_min           ! [Pa a] Minimum allowed viscosity (for stability, ~1e3)
        real(wp), intent(IN)  :: eps_0              ! [1/yr] Regularization constant (minimum strain rate, ~1e-6)

        real(wp) :: visc(size(ATT,1),size(ATT,2),size(ATT,3)) ! [Pa a] 3D viscosity field

        ! Local variables
        integer :: i, j, k, nx, ny, nz
        integer :: im1, ip1, jm1, jp1  
        real(wp) :: exp1, exp2
        real(wp) :: eps_0_sq, de_now 
        real(wp) :: wt

        nx = size(visc,1)
        ny = size(visc,2)
        nz = size(visc,3)

        ! Determine exponent values 
        exp1 = -1.0/n_glen
        exp2 = (1.0 - n_glen)/n_glen 

        eps_0_sq = eps_0*eps_0

        !$omp parallel do
        do k = 1, nz 
        do j = 1, ny 
        do i = 1, nx 

            if (f_ice(i,j) .eq. 1.0_wp) then 

                ! Calculate regularized strain rate 
                de_now = sqrt(de(i,j,k)**2 + eps_0_sq)

                ! Calculate viscosity at each aa-node
                visc(i,j,k) = 0.5_wp * ATT(i,j,k)**exp1 * (de_now)**exp2 

                ! Limit viscosity to above minimum value 
                if (visc(i,j,k) .lt. visc_min) visc(i,j,k) = visc_min 

            else 

                visc(i,j,k) = 0.0_wp 

            end if 

        end do 
        end do 
        end do
        !$omp end parallel do

        return
        
    end function calc_viscosity_glen

        function calc_viscosity_glen_2D(de,ATT,H_ice,f_ice,n_glen,visc_min,eps_0) result(visc)
        ! Calculate viscosity based on Glen's flow law 
        ! ATT [a^-1 Pa^-n] is the "depth dependent ice stiffness parameter based on
        !     vertical variations in temperature, chemistry and crystal fabric" (MacAyeal, 1989, JGR)
        ! de [a^-1] is the second-invariant of the strain rate tensor 
        ! visc [Pa a] is the 3D, temperature dependent viscosity field 

        ! Equation: visc = 0.5 * ATT^(-1/n_glen) * (de)^((1-n_glen)/n_glen)
        ! ATT  => from Greve and Blatter (2009), Eq. 4.15 (written as `A(T_prime)`)
        ! de   => from Greve and Blatter (2009), Eq. 6.53
        ! visc => from Greve and Blatter (2009), Eq. 4.22 

        implicit none
        
        real(wp), intent(IN)  :: de(:,:)            ! [a^-1] second-invariant of the strain rate tensor
        real(wp), intent(IN)  :: ATT(:,:,:)         ! [a^-1 Pa^-3] Rate factor 
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: n_glen             ! Glen's flow law exponent
        real(wp), intent(IN)  :: visc_min           ! [Pa a] Minimum allowed viscosity (for stability, ~1e3)
        real(wp), intent(IN)  :: eps_0              ! [1/yr] Regularization constant (minimum strain rate, ~1e-6)

        real(wp) :: visc(size(ATT,1),size(ATT,2),size(ATT,3)) ! [Pa a] 3D viscosity field

        ! Local variables
        integer :: i, j, k, nx, ny, nz
        integer :: im1, ip1, jm1, jp1  
        real(wp) :: exp1, exp2
        real(wp) :: eps_0_sq, de_now 
        real(wp) :: wt

        nx = size(visc,1)
        ny = size(visc,2)
        nz = size(visc,3)

        ! Determine exponent values 
        exp1 = -1.0/n_glen
        exp2 = (1.0 - n_glen)/n_glen 

        eps_0_sq = eps_0*eps_0

        !$omp parallel do
        do k = 1, nz 
        do j = 1, ny 
        do i = 1, nx 

            if (f_ice(i,j) .eq. 1.0_wp) then 

                ! Calculate regularized strain rate 
                de_now = sqrt(de(i,j)**2 + eps_0_sq)

                ! Calculate viscosity at each aa-node
                visc(i,j,k) = 0.5_wp * ATT(i,j,k)**exp1 * (de_now)**exp2 

                ! Limit viscosity to above minimum value 
                if (visc(i,j,k) .lt. visc_min) visc(i,j,k) = visc_min 

            else 

                visc(i,j,k) = 0.0_wp 

            end if 

        end do 
        end do 
        end do
        !$omp end parallel do

        return
        
    end function calc_viscosity_glen_2D


    subroutine calc_visc_int(visc_eff_int,visc_eff,H_ice,f_ice,zeta_aa,boundaries)

        implicit none 

        real(wp), intent(OUT) :: visc_eff_int(:,:)
        real(wp), intent(IN)  :: visc_eff(:,:,:)
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: zeta_aa(:)
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer :: i, j, nx, ny
        integer :: im1, ip1, jm1, jp1  
        real(wp) :: H_now
        real(wp) :: visc_eff_mean 
        real(wp) :: wt 

        nx = size(visc_eff_int,1)
        ny = size(visc_eff_int,2)


        do j = 1, ny 
        do i = 1, nx

            ! Calculate the vertically averaged viscosity for this point
            visc_eff_mean = integrate_trapezoid1D_pt(visc_eff(i,j,:),zeta_aa) 
            
            if (f_ice(i,j) .eq. 1.0) then 
                visc_eff_int(i,j) = visc_eff_mean*H_ice(i,j) 
            else
                !visc_eff_int(i,j) = visc_eff_mean 
                visc_eff_int(i,j) = 0.0_wp
            end if 

        end do 
        end do 

        ! Apply boundary conditions as needed 
        if (trim(boundaries) .eq. "periodic") then

            visc_eff_int(1,:)    = visc_eff_int(nx-1,:) 
            visc_eff_int(nx-1,:) = visc_eff_int(2,:) 
            visc_eff_int(:,1)    = visc_eff_int(:,ny-1)
            visc_eff_int(:,ny)   = visc_eff_int(:,2) 

        else if (trim(boundaries) .eq. "periodic-x") then 
            
            visc_eff_int(1,:)    = visc_eff_int(nx-1,:) 
            visc_eff_int(nx-1,:) = visc_eff_int(2,:) 
            visc_eff_int(:,1)    = visc_eff_int(:,2)
            visc_eff_int(:,ny)   = visc_eff_int(:,ny-1) 

        else if (trim(boundaries) .eq. "infinite") then 
            
            visc_eff_int(1,:)    = visc_eff_int(2,:) 
            visc_eff_int(nx,:)   = visc_eff_int(nx-1,:) 
            visc_eff_int(:,1)    = visc_eff_int(:,2)
            visc_eff_int(:,ny)   = visc_eff_int(:,ny-1) 

        end if 

        return

    end subroutine calc_visc_int

    elemental function calc_rate_factor(T_ice,T_pmp,enh) result(ATT)
        ! Greve and Blatter (2009): Chapter 4, page 54 
        ! Note: only valid for a Glen's flow law exponent n=3

        implicit none 

        real(prec), intent(IN) :: T_ice     ! [K]  Ice temperature
        real(prec), intent(IN) :: T_pmp     ! [K]  Pressure-corrected melting point
        real(prec), intent(IN) :: enh       !  [--] Enhancement factor 
        real(prec) :: ATT                   ! [a^-1 Pa^-3]

        ! Local variables
        real(prec) :: T_prime                         ! [K]  Temperature relative to the pressure melting point
        real(prec), parameter :: T_prime_lim = 263.15 ! [K] 
        real(prec), parameter :: A0_1 = 1.25671e-05   ! [a^-1 Pa^-3]
        real(prec), parameter :: A0_2 = 6.0422976e10  ! [a^-1 Pa^-3]
        real(prec), parameter :: Q_1  =  60.e3        ! [J mol^-1]
        real(prec), parameter :: Q_2  = 139.e3        ! [J mol^-1]
        real(prec), parameter :: R    = 8.314         ! [J mol^-1 K^-1]

        ! Calculate T_prime following Greve and Blatter (2009), Eq. 4.14 
        T_prime = T_ice - T_pmp + T0

        ! Limit T_prime to avoid under/overflows 
        T_prime = max(T_prime,220.0_wp)
        T_prime = min(T_prime,T_pmp)
        
        if (T_prime <= T_prime_lim) then 
            ATT = enh * A0_1 * exp(-Q_1/(R*T_prime))
        else 
            ATT = enh * A0_2 * exp(-Q_2/(R*T_prime))
        end if

        return 

    end function calc_rate_factor
    
    elemental function calc_rate_factor_eismint(T_ice,T_pmp,enh) result(ATT)
        ! Greve and Blatter (2009): Chapter 4, page 54 
        ! Note: only valid for a Glen's flow law exponent n=3

        implicit none 

        real(prec), intent(IN) :: T_ice     ! [K]  Ice temperature
        real(prec), intent(IN) :: T_pmp     ! [K]  Pressure-corrected melting point
        real(prec), intent(IN) :: enh       ! [--] Enhancement factor 
        real(prec) :: ATT                   ! [a^-1 Pa^-3]

        ! Local variables
        real(prec) :: T_prime                           ! [K]  Temperature relative to the pressure melting point
        real(prec), parameter :: T_prime_lim = 263.15   ! [K] 
        real(prec), parameter :: A0_1 = 1.139205e-05    ! [a^-1 Pa^-3]
        real(prec), parameter :: A0_2 = 5.459348198e10  ! [a^-1 Pa^-3]
        real(prec), parameter :: Q_1  =  60.e3          ! [J mol^-1]
        real(prec), parameter :: Q_2  = 139.e3          ! [J mol^-1]
        real(prec), parameter :: R    = 8.314           ! [J mol^-1 K^-1]

        ! Calculate T_prime following Greve and Blatter (2009), Eq. 4.14 
        T_prime = T_ice - T_pmp + T0

        if (T_prime <= T_prime_lim) then 
            ATT = enh * A0_1 * exp(-Q_1/(R*T_prime))
        else 
            ATT = enh * A0_2 * exp(-Q_2/(R*T_prime))
        end if

        return 

    end function calc_rate_factor_eismint
    
    elemental subroutine scale_rate_factor_water(ATT,omega)
        ! This routine scales the rate factor when 
        ! water content (omega) is present, ie, when
        ! the rate factor is for T_pmp 
        ! Parameterization following Greve and Blatter (2016) Eq. 14,
        ! following Lliboutry and Duval (1985):
        ! A = A(melting_temp)*(1.18125*omega*100)

        implicit none 

        real(prec), intent(INOUT) :: ATT        ! [a^-1 Pa^-3] Rate factor
        real(prec), intent(IN)    :: omega      ! [--] Water content fraction

        if (omega .gt. 0.0) ATT = ATT * (1.0+181.25*omega)

        return 

    end subroutine scale_rate_factor_water

    function calc_rate_factor_integrated(ATT,zeta_aa,n_glen) result(ATT_int)
        ! Greve and Blatter (2009), Chpt 5, page 82 

        implicit none 

        real(prec), intent(IN) :: ATT(:,:,:)
        real(prec), intent(IN) :: zeta_aa(:)
        real(prec), intent(IN) :: n_glen  
        real(prec) :: ATT_int(size(ATT,1),size(ATT,2),size(ATT,3))

        ! Local variables 
        integer :: i, j, nx, ny

        nx = size(ATT,1)
        ny = size(ATT,2) 

        ! Vertically integrated values of ATT to each vertical level
        do j = 1, ny 
        do i = 1, nx 
            ATT_int(i,j,:) = integrate_trapezoid1D_1D(ATT(i,j,:)*(1.0-zeta_aa)**n_glen,zeta_aa)
        end do 
        end do 

        return

    end function calc_rate_factor_integrated
    
    subroutine calc_strain_rate_tensor(strn, strn2D, vx, vy, vz, H_ice, f_ice, f_grnd,  &
                    zeta_aa, zeta_ac, dx, de_max, n_glen)
        ! -------------------------------------------------------------------------------
        !  Computation of all components of the strain-rate tensor, the full
        !  effective strain rate and the shear fraction.
        !  Alexander Robinson: Adapted from sicopolis5-dev::calc_dxyz 
        ! ------------------------------------------------------------------------------

        ! Note: vx, vy are staggered on ac-nodes in the horizontal, but are on the zeta_aa nodes (ie layer-centered)
        ! in the vertical. vz is centered on aa-nodes in the horizontal, but staggered on zeta_ac nodes
        ! in the vertical. 

        ! Note: first calculate each tensor component on ab-nodes, then interpolate to aa-nodes)
        ! This is a quadrature approach and is generally more stable. 
        ! The temperorary variable dd_ab(1:4) is used to hold the values 
        ! calculated at each corner (ab-node), starting from dd_ab(1)==upper-right, and
        ! moving counter-clockwise. The average of dd_ab(1:4) gives the cell-centered
        ! (aa-node) value.

        implicit none
        
        type(strain_3D_class), intent(INOUT) :: strn            ! [yr^-1] on aa-nodes (3D)
        type(strain_2D_class), intent(INOUT) :: strn2D          ! [yr^-1] on aa-nodes (2D)
        real(wp), intent(IN) :: vx(:,:,:)                       ! nx,ny,nz_aa
        real(wp), intent(IN) :: vy(:,:,:)                       ! nx,ny,nz_aa
        real(wp), intent(IN) :: vz(:,:,:)                       ! nx,ny,nz_ac
        real(wp), intent(IN) :: H_ice(:,:)
        real(wp), intent(IN) :: f_ice(:,:)
        real(wp), intent(IN) :: f_grnd(:,:)
        real(wp), intent(IN) :: zeta_aa(:) 
        real(wp), intent(IN) :: zeta_ac(:) 
        real(wp), intent(IN) :: dx
        real(wp), intent(IN) :: de_max                          ! [yr^-1] Maximum allowed effective strain rate
        real(wp), intent(IN) :: n_glen

        ! Local variables 
        integer  :: i, j, k
        integer  :: im1, ip1, jm1, jp1 
        integer  :: nx, ny, nz_aa, nz_ac  
        real(wp) :: dxi, deta, dzeta
        real(wp) :: dy  
        real(wp) :: dx_inv, dy_inv
        real(wp) :: dx_2_inv, dy_2_inv
        real(wp) :: H_ice_inv
        real(wp) :: lxy, lyx, lxz, lzx, lyz, lzy
        real(wp) :: shear_squared 
        real(wp) :: ux_aa, uy_aa 
        real(wp), allocatable :: fact_x(:,:), fact_y(:,:)
        real(wp), allocatable :: fact_z(:)

        real(wp) :: wt 
        real(wp) :: dd_ab(4) 
        real(wp) :: dd_ab_up(4)
        real(wp) :: dd_ab_dn(4)
        real(wp) :: wt_ab(4) 

        ! Define dy 
        dy = dx 

        ! Determine sizes and allocate local variables 
        nx    = size(vx,1)
        ny    = size(vx,2)
        nz_aa = size(zeta_aa,1)
        nz_ac = size(zeta_ac,1)
        
        allocate(fact_x(nx,ny))
        allocate(fact_y(nx,ny))
        allocate(fact_z(nz_aa))

        ! Change arguments to local (sicopolis) variable names 
        dxi  = dx 
        deta = dy 

        !-------- Term abbreviations --------

        dx_inv = 1.0_wp/dx
        dy_inv = 1.0_wp/dy

        dx_2_inv = 1.0_wp/dx
        dy_2_inv = 1.0_wp/dy

        fact_x   = dx_inv
        fact_y   = dy_inv

        fact_z(1) = 1.0_wp/(zeta_aa(2)-zeta_aa(1))
        do k = 2, nz_aa-1 
            fact_z(k) = 1.0_wp/(zeta_aa(k+1)-zeta_aa(k-1))
        end do
        fact_z(nz_aa) = 1.0_wp/(zeta_aa(nz_aa)-zeta_aa(nz_aa-1))

        !-------- Initialisation --------

        strn%dxx          = 0.0_wp
        strn%dyy          = 0.0_wp
        strn%dxy          = 0.0_wp
        strn%dxz          = 0.0_wp
        strn%dyz          = 0.0_wp
        strn%de           = 0.0_wp
        strn%div          = 0.0_wp 
        strn%f_shear      = 0.0_wp

        !-------- Computation --------

        !$omp parallel do
        do j=1, ny
        do i=1, nx

            ! Get neighbor indices
            im1 = max(i-1,1) 
            ip1 = min(i+1,nx) 
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny) 
            
            if (f_ice(i,j) .eq. 1.0_wp) then 

                H_ice_inv = 1.0_wp/H_ice(i,j)

                ! Get equal weighting for each ab-node
                wt_ab = 0.25_wp 

                ! ====== Loop over each column ====== 

                do k = 1, nz_aa 

                    ! === dxx =================================

                    call staggerdiff_nodes_acx_ab_ice(dd_ab,vx(:,:,k),f_ice,i,j,dx)

                    strn%dxx(i,j,k) = sum(wt_ab*dd_ab)
                    if (abs(strn%dxx(i,j,k)) .lt. TOL_UNDERFLOW) strn%dxx(i,j,k) = 0.0 
                    
                    ! === dyy =================================

                    call staggerdiff_nodes_acy_ab_ice(dd_ab,vy(:,:,k),f_ice,i,j,dy)

                    strn%dyy(i,j,k) = sum(wt_ab*dd_ab)
                    if (abs(strn%dyy(i,j,k)) .lt. TOL_UNDERFLOW) strn%dyy(i,j,k) = 0.0 

                    ! === lxy =================================

                    call staggerdiffcross_nodes_acx_ab_ice(dd_ab,vx(:,:,k),f_ice,i,j,dy)

                    lxy = sum(wt_ab*dd_ab)

                    ! === lyx =================================

                    call staggerdiffcross_nodes_acy_ab_ice(dd_ab,vy(:,:,k),f_ice,i,j,dx)

                    lyx = sum(wt_ab*dd_ab)

                    ! === dxy ================================= 

                    strn%dxy(i,j,k) = 0.5_wp*(lxy+lyx)
                    if (abs(strn%dxy(i,j,k)) .lt. TOL_UNDERFLOW) strn%dxy(i,j,k) = 0.0 
                    

                    ! ====== Vertical cross terms (lzx,lzy) ====== 

                    ! === lzx ================================

                    if (k .eq. 1) then
                        ! Basal layer 

                        call staggerdiff_nodes_acz_dx_ab_ice(dd_ab,vz(:,:,k),f_ice,i,j,dx)
                    
                    else if (k .eq. nz_aa) then
                        ! Surface layer 

                        call staggerdiff_nodes_acz_dx_ab_ice(dd_ab,vz(:,:,k+1),f_ice,i,j,dx)
                    
                    else 
                        ! Intermediate layers

                        call staggerdiff_nodes_acz_dx_ab_ice(dd_ab_up,vz(:,:,k+1),f_ice,i,j,dx)
                        call staggerdiff_nodes_acz_dx_ab_ice(dd_ab_dn,vz(:,:,k),  f_ice,i,j,dx)
                        
                        dd_ab = 0.5_wp*(dd_ab_dn+dd_ab_up)

                    end if 

                    lzx = sum(wt_ab*dd_ab)

                    ! === lzy ================================

                    if (k .eq. 1) then
                        ! Basal layer

                        call staggerdiff_nodes_acz_dy_ab_ice(dd_ab,vz(:,:,k),f_ice,i,j,dy)

                    else if (k .eq. nz_aa) then 
                        ! Surface layer

                        call staggerdiff_nodes_acz_dy_ab_ice(dd_ab,vz(:,:,k+1),f_ice,i,j,dy)

                    else 
                        ! Intermediate layers

                        call staggerdiff_nodes_acz_dy_ab_ice(dd_ab_up,vz(:,:,k+1),f_ice,i,j,dy)
                        call staggerdiff_nodes_acz_dy_ab_ice(dd_ab_dn,vz(:,:,k),  f_ice,i,j,dy)
                        
                        dd_ab = 0.5_wp*(dd_ab_dn+dd_ab_up)
                    end if 

                    lzy = sum(wt_ab*dd_ab)

                    ! ====== Shear terms (lxz,lyz) ================= 

                    ! === lxz ================================

                    if (k .eq. 1) then 
                        ! Basal layer
                        ! Gradient from first aa-node above base to base 

                        call stagger_nodes_acx_ab_ice(dd_ab_up,vx(:,:,k+1),f_ice,i,j)
                        call stagger_nodes_acx_ab_ice(dd_ab_dn,vx(:,:,k),  f_ice,i,j)

                        dd_ab = (dd_ab_up - dd_ab_dn)*fact_z(k)*H_ice_inv

                    else if (k .eq. nz_aa) then 
                        ! Surface layer
                        ! Gradient from surface to first aa-node below surface 

                        call stagger_nodes_acx_ab_ice(dd_ab_up,vx(:,:,k),  f_ice,i,j)
                        call stagger_nodes_acx_ab_ice(dd_ab_dn,vx(:,:,k-1),f_ice,i,j)

                        dd_ab = (dd_ab_up - dd_ab_dn)*fact_z(k)*H_ice_inv
                        
                    else 
                        ! Intermediate layers
                        ! Gradient from aa-node above to aa-node below

                        call stagger_nodes_acx_ab_ice(dd_ab_up,vx(:,:,k+1),f_ice,i,j)
                        call stagger_nodes_acx_ab_ice(dd_ab_dn,vx(:,:,k-1),f_ice,i,j)

                        dd_ab = (dd_ab_up - dd_ab_dn)*fact_z(k)*H_ice_inv
                        
                    end if 

                    lxz = sum(wt_ab*dd_ab)

                    ! === lyz ================================

                    if (k .eq. 1) then 
                        ! Basal layer
                        ! Gradient from first aa-node above base to base 

                        call stagger_nodes_acy_ab_ice(dd_ab_up,vy(:,:,k+1),f_ice,i,j)
                        call stagger_nodes_acy_ab_ice(dd_ab_dn,vy(:,:,k),  f_ice,i,j)

                        dd_ab = (dd_ab_up - dd_ab_dn)*fact_z(k)*H_ice_inv
                        
                    else if (k .eq. nz_aa) then 
                        ! Surface layer
                        ! Gradient from surface to first aa-node below surface 

                        call stagger_nodes_acy_ab_ice(dd_ab_up,vy(:,:,k),  f_ice,i,j)
                        call stagger_nodes_acy_ab_ice(dd_ab_dn,vy(:,:,k-1),f_ice,i,j)

                        dd_ab = (dd_ab_up - dd_ab_dn)*fact_z(k)*H_ice_inv
                        
                    else 
                        ! Intermediate layers
                        ! Gradient from aa-node above to aa-node below

                        call stagger_nodes_acy_ab_ice(dd_ab_up,vy(:,:,k+1),f_ice,i,j)
                        call stagger_nodes_acy_ab_ice(dd_ab_dn,vy(:,:,k-1),f_ice,i,j)

                        dd_ab = (dd_ab_up - dd_ab_dn)*fact_z(k)*H_ice_inv
                        
                    end if 

                    lyz = sum(wt_ab*dd_ab)

                    ! ====== Calculate cross terms from intermediate values (dxy,dxz,dyz) ====== 

                    strn%dxz(i,j,k) = 0.5_wp*(lxz+lzx)
                    strn%dyz(i,j,k) = 0.5_wp*(lyz+lzy)

                    ! Avoid extreme values
                    if (strn%dxz(i,j,k) .lt. -de_max) strn%dxz(i,j,k) = -de_max 
                    if (strn%dxz(i,j,k) .gt.  de_max) strn%dxz(i,j,k) =  de_max 

                    if (strn%dyz(i,j,k) .lt. -de_max) strn%dyz(i,j,k) = -de_max 
                    if (strn%dyz(i,j,k) .gt.  de_max) strn%dyz(i,j,k) =  de_max 
                    
                    ! Avoid underflows 
                    if (abs(strn%dxz(i,j,k)) .lt. TOL_UNDERFLOW) strn%dxz(i,j,k) = 0.0 
                    if (abs(strn%dyz(i,j,k)) .lt. TOL_UNDERFLOW) strn%dyz(i,j,k) = 0.0 
    
                    ! ====== Finished calculating individual strain rate terms ====== 
                    
                    strn%de(i,j,k) =  sqrt(  strn%dxx(i,j,k)*strn%dxx(i,j,k) &
                                           + strn%dyy(i,j,k)*strn%dyy(i,j,k) &
                                           + strn%dxx(i,j,k)*strn%dyy(i,j,k) &
                                           + strn%dxy(i,j,k)*strn%dxy(i,j,k) &
                                           + strn%dxz(i,j,k)*strn%dxz(i,j,k) &
                                           + strn%dyz(i,j,k)*strn%dyz(i,j,k) )
                    
                    if (strn%de(i,j,k) .gt. de_max) strn%de(i,j,k) = de_max 

                    ! Calculate the horizontal divergence too 
                    strn%div(i,j,k) = strn%dxx(i,j,k) + strn%dyy(i,j,k) 

                    ! Note: Using only the below should be equivalent to applying
                    ! the SIA approximation to calculate `de`
                    !strn%de(i,j,k)    =  sqrt( shear_squared(k) )

                    if (strn%de(i,j,k) .gt. 0.0) then 
                        ! Calculate the shear-based strain, stretching and the shear-fraction
                        shear_squared  =   strn%dxz(i,j,k)*strn%dxz(i,j,k) &
                                         + strn%dyz(i,j,k)*strn%dyz(i,j,k)
                        strn%f_shear(i,j,k) = sqrt(shear_squared)/strn%de(i,j,k)
                    else 
                        strn%f_shear(i,j,k) = 1.0   ! Shearing by default for low strain rates
                    end if 

                    !  ------ Modification of the shear fraction for floating ice (ice shelves)

                    if (f_grnd(i,j) .eq. 0.0) then 
                        strn%f_shear(i,j,k) = 0.0    ! Assume ice shelf is only stretching, no shear 
                    end if 

                    !  ------ Constrain the shear fraction to reasonable [0,1] interval

                    strn%f_shear(i,j,k) = min(max(strn%f_shear(i,j,k), 0.0), 1.0)

                end do 

            end if ! ice-free or ice-covered 

        end do
        end do
        !$omp end parallel do

        ! === Also calculate vertically averaged strain rate tensor ===
        
        ! Get the 2D average of strain rate in case it is needed 
        strn2D%dxx     = calc_vertical_integrated_2D(strn%dxx, zeta_aa)
        strn2D%dyy     = calc_vertical_integrated_2D(strn%dyy, zeta_aa)
        strn2D%dxy     = calc_vertical_integrated_2D(strn%dxy, zeta_aa)
        strn2D%dxz     = calc_vertical_integrated_2D(strn%dxz, zeta_aa)
        strn2D%dyz     = calc_vertical_integrated_2D(strn%dyz, zeta_aa)
        strn2D%div     = calc_vertical_integrated_2D(strn%div, zeta_aa)
        strn2D%de      = calc_vertical_integrated_2D(strn%de,  zeta_aa)
        strn2D%f_shear = calc_vertical_integrated_2D(strn%f_shear,zeta_aa) 
        
        return 

    end subroutine calc_strain_rate_tensor

    subroutine calc_strain_rate_tensor_aa(strn, strn2D, vx, vy, vz, H_ice, f_ice, f_grnd, &
                    zeta_aa, zeta_ac, dx, de_max, n_glen)
        ! -------------------------------------------------------------------------------
        !  Computation of all components of the strain-rate tensor, the full
        !  effective strain rate and the shear fraction.
        !  Alexander Robinson: Adapted from sicopolis5-dev::calc_dxyz 
        ! ------------------------------------------------------------------------------

        ! Note: vx, vy are staggered on ac-nodes in the horizontal, but are on the zeta_aa nodes (ie layer-centered)
        ! in the vertical. vz is centered on aa-nodes in the horizontal, but staggered on zeta_ac nodes
        ! in the vertical. 

        ! Note: this routine has been deprecated (2021-06-14) and will be removed, 
        ! given that the quadrature approach (calculating values on ab-nodes) appears to be more stable. 

        implicit none
        
        type(strain_3D_class), intent(INOUT) :: strn            ! [yr^-1] on aa-nodes (3D)
        type(strain_2D_class), intent(INOUT) :: strn2D          ! [yr^-1] on aa-nodes (2D)
        real(wp), intent(IN) :: vx(:,:,:)                       ! nx,ny,nz_aa
        real(wp), intent(IN) :: vy(:,:,:)                       ! nx,ny,nz_aa
        real(wp), intent(IN) :: vz(:,:,:)                       ! nx,ny,nz_ac
        real(wp), intent(IN) :: H_ice(:,:)
        real(wp), intent(IN) :: f_ice(:,:)
        real(wp), intent(IN) :: f_grnd(:,:)
        real(wp), intent(IN) :: zeta_aa(:) 
        real(wp), intent(IN) :: zeta_ac(:) 
        real(wp), intent(IN) :: dx
        real(wp), intent(IN) :: de_max                          ! [yr^-1] Maximum allowed effective strain rate
        real(wp), intent(IN) :: n_glen

        ! Local variables 
        integer  :: i, j, k
        integer  :: im1, ip1, jm1, jp1 
        integer  :: nx, ny, nz_aa, nz_ac  
        real(wp) :: dxi, deta, dzeta
        real(wp) :: dy  
        real(wp) :: dx_inv, dy_inv
        real(wp) :: H_ice_inv
        real(wp) :: lxy, lyx, lxz, lzx, lyz, lzy
        real(wp) :: shear_squared 
        real(wp) :: ux_aa, uy_aa 
        real(wp), allocatable :: fact_x(:,:), fact_y(:,:)
        real(wp), allocatable :: fact_z(:)

        real(wp) :: wt 

        ! Define dy 
        dy = dx 

        ! Determine sizes and allocate local variables 
        nx    = size(vx,1)
        ny    = size(vx,2)
        nz_aa = size(zeta_aa,1)
        nz_ac = size(zeta_ac,1)
        
        allocate(fact_x(nx,ny))
        allocate(fact_y(nx,ny))
        allocate(fact_z(nz_aa))

        ! Change arguments to local (sicopolis) variable names 
        dxi  = dx 
        deta = dy 

        !-------- Term abbreviations --------

        dx_inv = 1.0_wp/dx
        dy_inv = 1.0_wp/dy

        fact_x   = dx_inv
        fact_y   = dy_inv

        fact_z(1) = 1.0_wp/(zeta_aa(2)-zeta_aa(1))
        do k = 2, nz_aa-1 
            fact_z(k) = 1.0_wp/(zeta_aa(k+1)-zeta_aa(k-1))
        end do
        fact_z(nz_aa) = 1.0_wp/(zeta_aa(nz_aa)-zeta_aa(nz_aa-1))

        !-------- Initialisation --------

        strn%dxx          = 0.0_wp
        strn%dyy          = 0.0_wp
        strn%dxy          = 0.0_wp
        strn%dxz          = 0.0_wp
        strn%dyz          = 0.0_wp
        strn%de           = 0.0_wp
        strn%f_shear      = 0.0_wp

        !-------- Computation --------

        !$omp parallel do
        do j=1, ny
        do i=1, nx

            ! Get neighbor indices
            im1 = max(i-1,1) 
            ip1 = min(i+1,nx) 
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny) 
            
            if (f_ice(i,j) .eq. 1.0_wp) then 

                H_ice_inv = 1.0_wp/H_ice(i,j)

                
                ! ====== Loop over each column ====== 

                do k = 1, nz_aa 

                    ! ====== Horizontal gradients (dxx,dyy,lxy,lyx) ======

                    ! Full column from base to surface  
                    strn%dxx(i,j,k) = (vx(i,j,k)-vx(im1,j,k))*fact_x(i,j)
                    strn%dyy(i,j,k) = (vy(i,j,k)-vy(i,jm1,k))*fact_y(i,j)

                    ! Avoid underflows 
                    if (abs(strn%dxx(i,j,k)) .lt. TOL_UNDERFLOW) strn%dxx(i,j,k) = 0.0 
                    if (abs(strn%dyy(i,j,k)) .lt. TOL_UNDERFLOW) strn%dyy(i,j,k) = 0.0 
                    
                    lxy = (  (vx(i,jp1,k)+vx(im1,jp1,k))   &
                           - (vx(i,jm1,k)+vx(im1,jm1,k)) ) &
                            *0.25*fact_y(i,j)

                    lyx = (  (vy(ip1,j,k)+vy(ip1,jm1,k))   &
                           - (vy(im1,j,k)+vy(im1,jm1,k)) ) &
                            *0.25*fact_x(i,j)

                    ! ====== Vertical cross terms (lzx,lzy) ====== 

                    if (k .eq. 1) then
                        ! Basal layer

                        lzx = (vz(ip1,j,k)-vz(im1,j,k))*0.5*fact_x(i,j)
                        lzy = (vz(i,jp1,k)-vz(i,jm1,k))*0.5*fact_y(i,j)

                    else if (k .eq. nz_aa) then 
                        ! Surface layer

                        lzx = (vz(ip1,j,k+1)-vz(im1,j,k+1))*0.5*fact_x(i,j)
                        lzy = (vz(i,jp1,k+1)-vz(i,jm1,k+1))*0.5*fact_y(i,j)

                    else 
                        ! Intermediate layers 

                        lzx = (  (vz(ip1,j,k+1)+vz(ip1,j,k))   &
                               - (vz(im1,j,k+1)+vz(im1,j,k)) ) &
                                *0.25*fact_x(i,j)

                        lzy = (  (vz(i,jp1,k+1)+vz(i,jp1,k))   &
                               - (vz(i,jm1,k+1)+vz(i,jm1,k)) ) &
                                *0.25*fact_y(i,j)

                    end if 
                    
                    ! ====== Shear terms (lxz,lyz) ====== 

                    if (k .eq. 1) then 
                        ! Basal layer

                        ! Gradient from first aa-node above base to base 
                        lxz = (  (vx(i,j,k+1)+vx(im1,j,k+1))   &
                               - (vx(i,j,  k)+vx(im1,j,  k)) ) &
                                *0.5*fact_z(k)*H_ice_inv

                        lyz = (  (vy(i,j,k+1)+vy(i,jm1,k+1))   &
                               - (vy(i,j,  k)+vy(i,jm1,  k)) ) &
                                *0.5*fact_z(k)*H_ice_inv

                    else if (k .eq. nz_aa) then 
                        ! Surface layer

                        ! Gradient from surface to first aa-node below surface 
                        lxz = (  (vx(i,j,  k)+vx(im1,j,  k))   &
                               - (vx(i,j,k-1)+vx(im1,j,k-1)) ) &
                                *0.5*fact_z(k)*H_ice_inv

                        lyz = (  (vy(i,j,  k)+vy(i,jm1,  k))   &
                               - (vy(i,j,k-1)+vy(i,jm1,k-1)) ) &
                                *0.5*fact_z(k)*H_ice_inv

                    else 
                        ! Intermediate layers

                        ! Gradient from aa-node above to aa-node below
                        lxz = (  (vx(i,j,k+1)+vx(im1,j,k+1))   &
                               - (vx(i,j,k-1)+vx(im1,j,k-1)) ) &
                                *0.5*fact_z(k)*H_ice_inv

                        lyz = (  (vy(i,j,k+1)+vy(i,jm1,k+1))   &
                               - (vy(i,j,k-1)+vy(i,jm1,k-1)) ) &
                                *0.5*fact_z(k)*H_ice_inv 

                    end if 

                    ! ====== Calculate cross terms from intermediate values (dxy,dxz,dyz) ====== 

                    strn%dxy(i,j,k) = 0.5*(lxy+lyx)
                    strn%dxz(i,j,k) = 0.5*(lxz+lzx)
                    strn%dyz(i,j,k) = 0.5*(lyz+lzy)

                    ! Avoid underflows 
                    if (abs(strn%dxy(i,j,k)) .lt. TOL_UNDERFLOW) strn%dxy(i,j,k) = 0.0 
                    if (abs(strn%dxz(i,j,k)) .lt. TOL_UNDERFLOW) strn%dxz(i,j,k) = 0.0 
                    if (abs(strn%dyz(i,j,k)) .lt. TOL_UNDERFLOW) strn%dyz(i,j,k) = 0.0 
    
                    ! ====== Finished calculating individual strain rate terms ====== 
                    
                    strn%de(i,j,k) =  sqrt(  strn%dxx(i,j,k)*strn%dxx(i,j,k) &
                                           + strn%dyy(i,j,k)*strn%dyy(i,j,k) &
                                           + strn%dxx(i,j,k)*strn%dyy(i,j,k) &
                                           + strn%dxy(i,j,k)*strn%dxy(i,j,k) &
                                           + strn%dxz(i,j,k)*strn%dxz(i,j,k) &
                                           + strn%dyz(i,j,k)*strn%dyz(i,j,k) )
                    
                    if (strn%de(i,j,k) .gt. de_max) strn%de(i,j,k) = de_max 

                    ! Calculate the horizontal divergence too 
                    strn%div(i,j,k) = strn%dxx(i,j,k) + strn%dyy(i,j,k) 

                    ! Note: Using only the below should be equivalent to applying
                    ! the SIA approximation to calculate `de`
                    !strn%de(i,j,k)    =  sqrt( shear_squared(k) )

                    if (strn%de(i,j,k) .gt. 0.0) then 
                        ! Calculate the shear-based strain, stretching and the shear-fraction
                        shear_squared  =   strn%dxz(i,j,k)*strn%dxz(i,j,k) &
                                         + strn%dyz(i,j,k)*strn%dyz(i,j,k)
                        strn%f_shear(i,j,k) = sqrt(shear_squared)/strn%de(i,j,k)
                    else 
                        strn%f_shear(i,j,k) = 1.0   ! Shearing by default for low strain rates
                    end if 

                    !  ------ Modification of the shear fraction for floating ice (ice shelves)

                    if (f_grnd(i,j) .eq. 0.0) then 
                        strn%f_shear(i,j,k) = 0.0    ! Assume ice shelf is only stretching, no shear 
                    end if 

                    !  ------ Constrain the shear fraction to reasonable [0,1] interval

                    strn%f_shear(i,j,k) = min(max(strn%f_shear(i,j,k), 0.0), 1.0)

                end do

            end if ! ice-free or ice-covered 

        end do
        end do
        !$omp end parallel do

        ! === Also calculate vertically averaged strain rate tensor ===
        
        ! Get the 2D average of strain rate in case it is needed 
        strn2D%dxx     = calc_vertical_integrated_2D(strn%dxx, zeta_aa)
        strn2D%dyy     = calc_vertical_integrated_2D(strn%dyy, zeta_aa)
        strn2D%dxy     = calc_vertical_integrated_2D(strn%dxy, zeta_aa)
        strn2D%dxz     = calc_vertical_integrated_2D(strn%dxz, zeta_aa)
        strn2D%dyz     = calc_vertical_integrated_2D(strn%dyz, zeta_aa)
        strn2D%div     = calc_vertical_integrated_2D(strn%div, zeta_aa)
        strn2D%de      = calc_vertical_integrated_2D(strn%de,  zeta_aa)
        strn2D%f_shear = calc_vertical_integrated_2D(strn%f_shear,zeta_aa) 
        
        return 

    end subroutine calc_strain_rate_tensor_aa

    subroutine calc_strain_rate_tensor_2D(strn2D,vx,vy,H_ice,f_ice,f_grnd,dx,de_max,n_glen)
        ! Calculate the 2D (vertically averaged) strain rate tensor,
        ! assuming a constant vertical velocity profile. 

        ! ajr: to do: perform calculation on ab-nodes (quadrature)
        ! to be consistent with new formulation above for 
        ! calc_strain_rate_tensor. 

        implicit none

        type(strain_2D_class), intent(INOUT) :: strn2D          ! [yr^-1] Strain rate tensor
        real(wp), intent(IN) :: vx(:,:)                         ! [m/yr] Vertically averaged vel., x
        real(wp), intent(IN) :: vy(:,:)                         ! [m/yr] Vertically averaged vel., y
        real(wp), intent(IN) :: H_ice(:,:)                      ! [m] Ice thickness 
        real(wp), intent(IN) :: f_ice(:,:)
        real(wp), intent(IN) :: f_grnd(:,:)
        real(wp), intent(IN) :: dx                              ! [m] Resolution
        real(wp), intent(IN) :: de_max                          ! [yr^-1] Maximum allowed effective strain rate
        real(wp), intent(IN) :: n_glen 

        ! Local variables
        integer  :: i, j, k
        integer  :: im1, ip1, jm1, jp1
        integer  :: nx, ny
        real(wp) :: dy 

        real(wp) :: lxy, lyx 
        real(wp) :: wt_ab(4) 
        real(wp) :: dd_ab(4) 

        nx = size(vx,1)
        ny = size(vy,2)

        ! Assume square grid 
        dy = dx 

        strn2D%dxx      = 0.0 
        strn2D%dyy      = 0.0 
        strn2D%dxy      = 0.0 
        strn2D%dxz      = 0.0       ! Always zero in this case
        strn2D%dyz      = 0.0       ! Always zero in this case
        strn2D%de       = 0.0 
        strn2D%div      = 0.0 
        strn2D%f_shear  = 0.0       ! Always zero in this case

        !$omp parallel do
        do j=1, ny
        do i=1, nx

            ! Get neighbor indices
            im1 = max(i-1,1) 
            ip1 = min(i+1,nx) 
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny) 
            
            if (f_ice(i,j) .eq. 1.0_wp) then 

                ! Get equal weighting for each ab-node
                wt_ab = 0.25_wp 

                ! === dxx =================================

                call staggerdiff_nodes_acx_ab_ice(dd_ab,vx,f_ice,i,j,dx)

                strn2D%dxx(i,j) = sum(wt_ab*dd_ab)
                if (abs(strn2D%dxx(i,j)) .lt. TOL_UNDERFLOW) strn2D%dxx(i,j) = 0.0 
                
                ! === dyy =================================

                call staggerdiff_nodes_acy_ab_ice(dd_ab,vy,f_ice,i,j,dy)

                strn2D%dyy(i,j) = sum(wt_ab*dd_ab)
                if (abs(strn2D%dyy(i,j)) .lt. TOL_UNDERFLOW) strn2D%dyy(i,j) = 0.0 

                ! === lxy =================================

                call staggerdiffcross_nodes_acx_ab_ice(dd_ab,vx,f_ice,i,j,dy)

                lxy = sum(wt_ab*dd_ab)

                ! === lyx =================================

                call staggerdiffcross_nodes_acy_ab_ice(dd_ab,vy,f_ice,i,j,dx)

                lyx = sum(wt_ab*dd_ab)

                ! === dxy ================================= 

                strn2D%dxy(i,j) = 0.5_wp*(lxy+lyx)
                if (abs(strn2D%dxy(i,j)) .lt. TOL_UNDERFLOW) strn2D%dxy(i,j) = 0.0 
                
                ! ====== Finished calculating individual strain rate terms ====== 
                
                strn2D%de(i,j) =  sqrt(  strn2D%dxx(i,j)*strn2D%dxx(i,j) &
                                       + strn2D%dyy(i,j)*strn2D%dyy(i,j) &
                                       + strn2D%dxx(i,j)*strn2D%dyy(i,j) &
                                       + strn2D%dxy(i,j)*strn2D%dxy(i,j) )
                
                if (strn2D%de(i,j) .gt. de_max) strn2D%de(i,j) = de_max 

                ! Calculate the horizontal divergence too 
                strn2D%div(i,j) = strn2D%dxx(i,j) + strn2D%dyy(i,j) 

            end if ! ice-free or ice-covered 

        end do
        end do
        !$omp end parallel do

        return
        
    end subroutine calc_strain_rate_tensor_2D
    
    subroutine calc_stress_tensor(strs,strs2D,visc,strn,zeta_aa)
        ! Calculate the deviatoric stress tensor components [Pa]
        ! following from, eg, Thoma et al. (2014), Eq. 7.
        
        implicit none 

        type(stress_3D_class), intent(INOUT) :: strs 
        type(stress_2D_class), intent(INOUT) :: strs2D 
        real(wp),              intent(IN)    :: visc(:,:,:)
        type(strain_3D_class), intent(IN)    :: strn 
        real(wp),              intent(IN)    :: zeta_aa(:) 

        strs%txx = 2.0*visc*strn%dxx
        strs%tyy = 2.0*visc*strn%dyy
        strs%txy = 2.0*visc*strn%dxy
        strs%txz = 2.0*visc*strn%dxz
        strs%tyz = 2.0*visc*strn%dyz
        
        ! Next calculate the effective stress 
        ! analogous to the effective strain rate
        ! (or, eg, Lipscomb et al., 2019, Eq. 44)

        strs%te =  sqrt(  strs%txx*strs%txx &
                        + strs%tyy*strs%tyy &
                        + strs%txx*strs%tyy &
                        + strs%txy*strs%txy &
                        + strs%txz*strs%txz &
                        + strs%tyz*strs%tyz )

        ! === Also calculate vertically averaged stress tensor ===
        strs2D%txx = calc_vertical_integrated_2D(strs%txx,zeta_aa)
        strs2D%tyy = calc_vertical_integrated_2D(strs%tyy,zeta_aa)
        strs2D%txy = calc_vertical_integrated_2D(strs%txy,zeta_aa)
        strs2D%txz = calc_vertical_integrated_2D(strs%txz,zeta_aa)
        strs2D%tyz = calc_vertical_integrated_2D(strs%tyz,zeta_aa)
        strs2D%te  = calc_vertical_integrated_2D(strs%te, zeta_aa)
        
        ! Finally, calculate the first two eigenvectors for 2D stress tensor 
        call calc_stress_eigen_values(strs2D%tau_eig_1,strs2D%tau_eig_2, &
                                        strs2D%txx,strs2D%tyy,strs2D%txy)

        return 

    end subroutine calc_stress_tensor
    
    subroutine calc_stress_tensor_2D(strs2D,visc_bar,strn2D)
        ! Calculate the deviatoric stress tensor components [Pa]
        ! following from, eg, Thoma et al. (2014), Eq. 7.
        
        implicit none 

        type(stress_2D_class), intent(INOUT) :: strs2D 
        real(prec),            intent(IN)    :: visc_bar(:,:)
        type(strain_2D_class), intent(IN)    :: strn2D 

        strs2D%txx = 2.0*visc_bar*strn2D%dxx
        strs2D%tyy = 2.0*visc_bar*strn2D%dyy
        strs2D%txy = 2.0*visc_bar*strn2D%dxy
        
        ! Also calculate the effective stress 
        ! analogous to the effective strain rate
        ! (or, eg, Lipscomb et al., 2019, Eq. 44)

        strs2D%te =  sqrt(strs2D%txx*strs2D%txx &
                        + strs2D%tyy*strs2D%tyy &
                        + strs2D%txx*strs2D%tyy &
                        + strs2D%txy*strs2D%txy &
                        + strs2D%txz*strs2D%txz &
                        + strs2D%tyz*strs2D%tyz )

        ! Finally, calculate the first two eigenvectors for 2D stress tensor 
        call calc_stress_eigen_values(strs2D%tau_eig_1,strs2D%tau_eig_2, &
                                    strs2D%txx,strs2D%tyy,strs2D%txy)

        return 

    end subroutine calc_stress_tensor_2D
    
    elemental subroutine calc_stress_eigen_values(tau_eig_1,tau_eig_2,txx,tyy,txy)
        ! Calculate the first two eigenvectors of 2D deviatoric stress tensor 

        implicit none

        real(wp), intent(OUT) :: tau_eig_1              ! [Pa] Eigenvalue 1
        real(wp), intent(OUT) :: tau_eig_2              ! [Pa] Eigenvalue 2
        real(wp), intent(IN)  :: txx
        real(wp), intent(IN)  :: tyy
        real(wp), intent(IN)  :: txy

        ! Local variables
        real(wp) :: a, b, c, root
        real(wp) :: lambda1, lambda2 

        ! compute the eigenvalues of the vertically averaged stress tensor
        a = 1.0_wp
        b = -(txx + tyy)
        c = txx*tyy - txy*txy
        if (b*b - 4.0_wp*a*c > 0.0_wp) then   ! two real eigenvalues
            root = sqrt(b*b - 4.0_wp*a*c)
            lambda1 = (-b + root) / (2.0_wp*a)
            lambda2 = (-b - root) / (2.0_wp*a)
            if (lambda1 > lambda2) then
                tau_eig_1 = lambda1
                tau_eig_2 = lambda2
            else
                tau_eig_1 = lambda2
                tau_eig_2 = lambda1
            end if
        else 
            ! No eigenvalues, set to zero 
            tau_eig_1 = 0.0_wp 
            tau_eig_2 = 0.0_wp 
        end if  ! b^2 - 4ac > 0

        return

    end subroutine calc_stress_eigen_values


    subroutine strain_2D_alloc(strn2D,nx,ny)

        implicit none 

        type(strain_2D_class), intent(INOUT) :: strn2D 
        integer :: nx, ny 

        ! First make sure fields are deallocated
        ! TO DO 

        ! Allocate fields to desired dimensions 
        allocate(strn2D%dxx(nx,ny))
        allocate(strn2D%dyy(nx,ny))
        allocate(strn2D%dxy(nx,ny))
        allocate(strn2D%dxz(nx,ny))
        allocate(strn2D%dyz(nx,ny))
        allocate(strn2D%div(nx,ny))
        allocate(strn2D%de(nx,ny))
        allocate(strn2D%f_shear(nx,ny))

        strn2D%dxx     = 0.0 
        strn2D%dyy     = 0.0 
        strn2D%dxy     = 0.0
        strn2D%dxz     = 0.0
        strn2D%dyz     = 0.0
        strn2D%div     = 0.0
        strn2D%de      = 0.0 
        strn2D%f_shear = 0.0 

        return 

    end subroutine strain_2D_alloc

    subroutine stress_2D_alloc(strs2D,nx,ny)

        implicit none 

        type(stress_2D_class), intent(INOUT) :: strs2D 
        integer :: nx, ny 

        ! First make sure fields are deallocated
        ! TO DO 

        ! Allocate fields to desired dimensions 
        allocate(strs2D%txx(nx,ny))
        allocate(strs2D%tyy(nx,ny))
        allocate(strs2D%txy(nx,ny))
        allocate(strs2D%txz(nx,ny))
        allocate(strs2D%tyz(nx,ny))
        allocate(strs2D%te(nx,ny))
        allocate(strs2D%tau_eig_1(nx,ny))
        allocate(strs2D%tau_eig_2(nx,ny))
        
        strs2D%txx   = 0.0 
        strs2D%tyy   = 0.0 
        strs2D%txy   = 0.0
        strs2D%txz   = 0.0
        strs2D%tyz   = 0.0
        strs2D%te    = 0.0 
        strs2D%tau_eig_1 = 0.0 
        strs2D%tau_eig_2 = 0.0 
        
        return 

    end subroutine stress_2D_alloc
    
end module deformation



