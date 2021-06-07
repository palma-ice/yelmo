

module deformation 
    ! This module contains functions related to deformation calculations:
    ! Arrenius function
    ! Flow law
    ! viscosity
    ! strain rate??

    ! Note: 3D arrays defined such that first index (k=1) == base, and max index (k=nk) == surface 

    use yelmo_defs,  only : sp, dp, wp, prec, tol_underflow, T0, rho_ice, g, &
                        strain_2D_class, strain_3D_class, stress_2D_class, stress_3D_class
    use yelmo_tools, only : calc_vertical_integrated_2D, integrate_trapezoid1D_1D   

    implicit none 
    
    private
    public :: modify_enhancement_factor_bnd 
    public :: define_enhancement_factor_3D
    public :: define_enhancement_factor_2D
    public :: calc_viscosity_glen
    public :: calc_viscosity_glen_2D
    public :: calc_rate_factor
    public :: calc_rate_factor_eismint
    public :: scale_rate_factor_water
    public :: calc_rate_factor_inverted
    public :: calc_rate_factor_integrated
    public :: calc_strain_rate_tensor
    public :: calc_strain_rate_tensor_2D
    public :: calc_stress_tensor 
    public :: calc_stress_tensor_2D
    public :: calc_stress_eigen_values
    
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

    function calc_viscosity_glen(de,ATT,n_glen,visc_min) result(visc)
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
        
        real(prec), intent(IN)  :: de(:,:,:)        ! [a^-1] second-invariant of the strain rate tensor
        real(prec), intent(IN)  :: ATT(:,:,:)       ! [a^-1 Pa^-3] Rate factor 
        real(prec), intent(IN)  :: n_glen           ! Glen's flow law exponent
        real(prec), intent(IN)  :: visc_min         ! [Pa a] Minimum allowed viscosity (for stability, ~1e3)
        
        real(prec) :: visc(size(ATT,1),size(ATT,2),size(ATT,3)) ! [Pa a] 3D viscosity field

        ! Local variables
        integer :: k, nz  
        real(prec) :: exp1, exp2  
        real(prec), parameter :: de_0 = 1.0e-6      ! [a^-1] Bueler and Brown (2009), Eq. 26

        nz = size(visc,3)

        ! Determine exponent values 
        exp1 = -1.0/n_glen
        exp2 = (1.0 - n_glen)/n_glen 

        ! Calculate viscosity at each layer 
        visc = 0.5 * ATT**exp1 * (de + de_0)**exp2 

        ! Limit viscosity to above minimum value 
        where(visc .lt. visc_min) visc = visc_min 

        return
        
    end function calc_viscosity_glen

    function calc_viscosity_glen_2D(de,ATT,n_glen,visc_min) result(visc)
        ! Calculate viscosity based on Glen's flow law 
        ! ATT [a^-1 Pa^-3] is the "depth dependent ice stiffness parameter based on
        !     vertical variations in temperature, chemistry and crystal fabric" (MacAyeal, 1989, JGR)
        ! de [a^-1] is the second-invariant of the strain rate tensor 
        ! visc [Pa a] is the 2D, temperature dependent viscosity field 

        ! Equation: visc = 0.5 * ATT^(-1/n_glen) * (de)^((1-n_glen)/n_glen)
        ! ATT  => from Greve and Blatter (2009), Eq. 4.15 (written as `A(T_prime)`)
        ! de   => from Greve and Blatter (2009), Eq. 6.53
        ! visc => from Greve and Blatter (2009), Eq. 4.22 

        implicit none
           
        real(prec), intent(IN)  :: de(:,:)          ! [a^-1] second-invariant of the strain rate tensor
        real(prec), intent(IN)  :: ATT(:,:,:)       ! [a^-1 Pa^-3] Rate factor 
        real(prec), intent(IN)  :: n_glen           ! Glen's flow law exponent
        real(prec), intent(IN)  :: visc_min         ! [Pa a] Minimum allowed viscosity (for stability, ~1e3)
        
        real(prec) :: visc(size(ATT,1),size(ATT,2),size(ATT,3)) ! [Pa a] 3D viscosity field

        ! Local variables
        integer :: k, nz_aa  
        real(prec) :: exp1, exp2  
        real(prec), parameter :: de_0 = 1.0e-6    ! [a^-1] Bueler and Brown (2009), Eq. 26
        
        nz_aa = size(visc,3)

        ! Determine exponent values 
        exp1 = -1.0/n_glen
        exp2 = (1.0 - n_glen)/n_glen 

        ! Calculate viscosity at each layer 
        do k = 1, nz_aa 
            visc(:,:,k) = 0.5 * ATT(:,:,k)**exp1 * (de+de_0)**exp2
        end do 

        ! Limit viscosity to above minimum value 
        where(visc .lt. visc_min) visc = visc_min 

        return
        
    end function calc_viscosity_glen_2D

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

    elemental function calc_rate_factor_inverted(ATT,n_glen) result(BTT)
        ! BTT  => from Greve and Blatter (2009), Eq. 4.35 (written as `B(T_prime)`)
        ! ajr, currently not used, instead ATT**(1/n_glen) is used directly in subroutines 

        implicit none 

        real(prec), intent(IN) :: ATT      ! Rate factor
        real(prec), intent(IN) :: n_glen   ! [--] Glen law exponent (n=3)
        real(prec) :: BTT 

        BTT = ATT**(-1.0/n_glen) 

        return 

    end function calc_rate_factor_inverted
    
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
    
    subroutine calc_strain_rate_tensor(strn, strn2D, vx, vy, vz, H_ice, f_grnd, f_ice, &
                    zeta_aa, zeta_ac, dx, de_max, ATT_bar, n_glen)
        ! -------------------------------------------------------------------------------
        !  Computation of all components of the strain-rate tensor, the full
        !  effective strain rate and the shear fraction.
        !  Alexander Robinson: Adapted from sicopolis5-dev::calc_dxyz 
        ! ------------------------------------------------------------------------------

        ! Note: vx, vy are staggered on ac-nodes in the horizontal, but are on the zeta_aa nodes (ie layer-centered)
        ! in the vertical. vz is centered on aa-nodes in the horizontal, but staggered on zeta_ac nodes
        ! in the vertical. 

        implicit none
        
        type(strain_3D_class), intent(INOUT) :: strn            ! [yr^-1] on aa-nodes (3D)
        type(strain_2D_class), intent(INOUT) :: strn2D          ! [yr^-1] on aa-nodes (2D)
        real(wp), intent(IN) :: vx(:,:,:)                       ! nx,ny,nz_aa
        real(wp), intent(IN) :: vy(:,:,:)                       ! nx,ny,nz_aa
        real(wp), intent(IN) :: vz(:,:,:)                       ! nx,ny,nz_ac
        real(wp), intent(IN) :: H_ice(:,:)
        real(wp), intent(IN) :: f_grnd(:,:)
        real(wp), intent(IN) :: f_ice(:,:)
        real(wp), intent(IN) :: zeta_aa(:) 
        real(wp), intent(IN) :: zeta_ac(:) 
        real(wp), intent(IN) :: dx
        real(wp), intent(IN) :: de_max                          ! [yr^-1] Maximum allowed effective strain rate
        real(wp), intent(IN) :: ATT_bar(:,:)
        real(wp), intent(IN) :: n_glen

        ! Local variables 
        integer  :: i, j, k
        integer  :: im1, ip1, jm1, jp1 
        integer  :: nx, ny, nz_aa, nz_ac  
        real(wp) :: dxi, deta, dzeta
        real(wp) :: dy  
        real(wp) :: dx_inv, dy_inv
        real(wp) :: H_ice_inv
        real(wp) :: abs_v_ssa_inv, nx1, ny1
        real(wp) :: shear_x_help, shear_y_help
        real(wp) :: f_shear_help
        real(wp) :: lxy, lyx, lxz, lzx, lyz, lzy
        real(wp) :: shear_squared 
        real(wp) :: ux_aa, uy_aa 
        real(wp), allocatable :: fact_x(:,:), fact_y(:,:)
        real(wp), allocatable :: fact_z(:)

        logical :: is_margin 
        real(wp) :: ddiv_free, dxx_free, dyy_free 
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

            if (H_ice(i,j) .gt. 0.0_wp) then 
                ! Grounded or floating ice, calculate strain rate here

                H_ice_inv = 1.0_wp/H_ice(i,j)

                ! Get neighbor indices
                im1 = max(i-1,1) 
                ip1 = min(i+1,nx) 
                jm1 = max(j-1,1) 
                jp1 = min(j+1,ny) 
                
                ! ====== Loop over each column ====== 

                do k = 1, nz_aa 

                    ! ====== Horizontal gradients (dxx,dyy,lxy,lyx) ======

                    ! Full column from base to surface  
                    strn%dxx(i,j,k) = (vx(i,j,k)-vx(im1,j,k))*fact_x(i,j)
                    strn%dyy(i,j,k) = (vy(i,j,k)-vy(i,jm1,k))*fact_y(i,j)

                    ! Avoid underflows 
                    if (abs(strn%dxx(i,j,k)) .lt. tol_underflow) strn%dxx(i,j,k) = 0.0 
                    if (abs(strn%dyy(i,j,k)) .lt. tol_underflow) strn%dyy(i,j,k) = 0.0 
                    
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
                    if (abs(strn%dxy(i,j,k)) .lt. tol_underflow) strn%dxy(i,j,k) = 0.0 
                    if (abs(strn%dxz(i,j,k)) .lt. tol_underflow) strn%dxz(i,j,k) = 0.0 
                    if (abs(strn%dyz(i,j,k)) .lt. tol_underflow) strn%dyz(i,j,k) = 0.0 
    
                    ! ====== Finished calculating individual strain rate terms ====== 
                    
                    strn%de(i,j,k) =  sqrt(  strn%dxx(i,j,k)*strn%dxx(i,j,k) &
                                           + strn%dyy(i,j,k)*strn%dyy(i,j,k) &
                                           + strn%dxx(i,j,k)*strn%dyy(i,j,k) &
                                           + strn%dxy(i,j,k)*strn%dxy(i,j,k) &
                                           + strn%dxz(i,j,k)*strn%dxz(i,j,k) &
                                           + strn%dyz(i,j,k)*strn%dyz(i,j,k) )
                    
                    if (strn%de(i,j,k) .gt. de_max) strn%de(i,j,k) = de_max 

                    ! Calculate the horizontal divergence too 
                    strn%ddiv(i,j,k) = strn%dxx(i,j,k) + strn%dyy(i,j,k) 

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



if (.FALSE.) then
    ! ajr: Extrapolating to ice-free and partially ice-free neighbors, 
    ! as done further below, is more stable and convincing than 
    ! imposing the free-spreading strain rate. So this section is disabled. 

                ! Also estimate the free-spreading strain rate (Pollard et al., 2015, EPSL, Eq. B2.b)
                ! ddiv = A*(rho*g*h/4)^n = dxx + dyy
                ! assume equal spreading in both directions:
                ! dxx = dyy; de = 2*dxx
                ! dxx = de/2
                ddiv_free = ATT_bar(i,j) * (0.25*rho_ice*g*H_ice(i,j))**n_glen
                ! dxx_free  = ddiv_free / 2.0
                ! dyy_free  = dxx_free 
                if ( abs(0.5*vx(i,j,nz_aa)+vx(im1,j,nz_aa)) &
                      .gt. abs(0.5*vy(i,j,nz_aa)+vy(i,jm1,nz_aa)) ) then 
                    dxx_free  = ddiv_free
                    dyy_free  = 0.0 
                else 
                    dxx_free  = 0.0
                    dyy_free  = ddiv_free
                end if 
                

                is_margin = ( H_ice(i,j) .gt. 0.0_wp .and. f_ice(i,j) .lt. 1.0_wp .and. &
                    count([H_ice(im1,j),H_ice(ip1,j),H_ice(i,jm1),H_ice(i,jp1)] .eq. 0.0_wp) .gt. 0 )
                
                ! For partially covered grid cells at margin, or floating points
                ! ensure effective strain rate is not larger than free-spreading strain rate
                if (is_margin .or. &
                     (f_grnd(i,j) .eq. 0.0 .and. strn2D%ddiv(i,j) .gt. ddiv_free) ) then 
                    ! Overwrite above value and impose free-spreading strain 

                    strn%ddiv(i,j,:)    = ddiv_free
                    strn%dxx(i,j,:)     = dxx_free
                    strn%dyy(i,j,:)     = dyy_free
                    strn%dxy(i,j,:)     = 0.0
                    strn%dxz(i,j,:)     = 0.0 
                    strn%dyz(i,j,:)     = 0.0 

                    strn%de(i,j,k) =  sqrt(  strn%dxx(i,j,k)*strn%dxx(i,j,k) &
                                           + strn%dyy(i,j,k)*strn%dyy(i,j,k) &
                                           + strn%dxx(i,j,k)*strn%dyy(i,j,k) &
                                           + strn%dxy(i,j,k)*strn%dxy(i,j,k) &
                                           + strn%dxz(i,j,k)*strn%dxz(i,j,k) &
                                           + strn%dyz(i,j,k)*strn%dyz(i,j,k) )
                    
                    if (strn%de(i,j,k) .gt. de_max) strn%de(i,j,k) = de_max 

                    ! Calculate the horizontal divergence too 
                    strn%ddiv(i,j,k) = strn%dxx(i,j,k) + strn%dyy(i,j,k) 

                end if 
end if 



            end if ! ice-free or ice-covered 

        end do
        end do
        !$omp end parallel do


        ! === Extrapolate to ice-free neighbors === 
        ! (in case ice gets advected there)
        do j=1, ny
        do i=1, nx

            ! Get neighbor indices
            im1 = max(i-1,1) 
            ip1 = min(i+1,nx) 
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny) 
                
            if ( (H_ice(i,j) .eq. 0.0 .or. f_ice(i,j) .lt. 1.0) .and. &
                count([H_ice(im1,j),H_ice(ip1,j),H_ice(i,jm1),H_ice(i,jp1)] .gt. 0.0_wp) .gt. 0 ) then 
                ! Ice-free (or partially ice-free) with ice-covered neighbors


                strn%dxx(i,j,:) = 0.0 
                strn%dyy(i,j,:) = 0.0
                strn%dxy(i,j,:) = 0.0
                strn%dxz(i,j,:) = 0.0 
                strn%dyz(i,j,:) = 0.0

                wt = 0.0 

                if (H_ice(im1,j) .gt. 0.0 .and. f_ice(im1,j) .eq. 1.0) then 
                    strn%dxx(i,j,:) = strn%dxx(i,j,:) + strn%dxx(im1,j,:)
                    strn%dyy(i,j,:) = strn%dyy(i,j,:) + strn%dyy(im1,j,:)
                    strn%dxy(i,j,:) = strn%dxy(i,j,:) + strn%dxy(im1,j,:)
                    strn%dxz(i,j,:) = strn%dxz(i,j,:) + strn%dxz(im1,j,:)
                    strn%dyz(i,j,:) = strn%dyz(i,j,:) + strn%dyz(im1,j,:)
                    wt = wt + 1.0 
                end if 
                if (H_ice(ip1,j) .gt. 0.0 .and. f_ice(ip1,j) .eq. 1.0) then 
                    strn%dxx(i,j,:) = strn%dxx(i,j,:) + strn%dxx(ip1,j,:)
                    strn%dyy(i,j,:) = strn%dyy(i,j,:) + strn%dyy(ip1,j,:)
                    strn%dxy(i,j,:) = strn%dxy(i,j,:) + strn%dxy(ip1,j,:)
                    strn%dxz(i,j,:) = strn%dxz(i,j,:) + strn%dxz(ip1,j,:)
                    strn%dyz(i,j,:) = strn%dyz(i,j,:) + strn%dyz(ip1,j,:)
                    wt = wt + 1.0 
                end if
                if (H_ice(i,jm1) .gt. 0.0 .and. f_ice(i,jm1) .eq. 1.0) then 
                    strn%dxx(i,j,:) = strn%dxx(i,j,:) + strn%dxx(i,jm1,:)
                    strn%dyy(i,j,:) = strn%dyy(i,j,:) + strn%dyy(i,jm1,:)
                    strn%dxy(i,j,:) = strn%dxy(i,j,:) + strn%dxy(i,jm1,:)
                    strn%dxz(i,j,:) = strn%dxz(i,j,:) + strn%dxz(i,jm1,:)
                    strn%dyz(i,j,:) = strn%dyz(i,j,:) + strn%dyz(i,jm1,:)
                    wt = wt + 1.0 
                end if
                if (H_ice(i,jp1) .gt. 0.0 .and. f_ice(i,jp1) .eq. 1.0) then 
                    strn%dxx(i,j,:) = strn%dxx(i,j,:) + strn%dxx(i,jp1,:)
                    strn%dyy(i,j,:) = strn%dyy(i,j,:) + strn%dyy(i,jp1,:)
                    strn%dxy(i,j,:) = strn%dxy(i,j,:) + strn%dxy(i,jp1,:)
                    strn%dxz(i,j,:) = strn%dxz(i,j,:) + strn%dxz(i,jp1,:)
                    strn%dyz(i,j,:) = strn%dyz(i,j,:) + strn%dyz(i,jp1,:)
                    wt = wt + 1.0 
                end if

                if (wt .gt. 0.0) then 
                    strn%dxx(i,j,:) = strn%dxx(i,j,:) / wt
                    strn%dyy(i,j,:) = strn%dyy(i,j,:) / wt
                    strn%dxy(i,j,:) = strn%dxy(i,j,:) / wt
                    strn%dxz(i,j,:) = strn%dxz(i,j,:) / wt
                    strn%dyz(i,j,:) = strn%dyz(i,j,:) / wt
                end if 


                ! Obtain effective strain rate and divergence

                do k = 1, nz_aa 
                    ! ====== Finished calculating individual strain rate terms ====== 
                    
                    strn%de(i,j,k) =  sqrt(  strn%dxx(i,j,k)*strn%dxx(i,j,k) &
                                           + strn%dyy(i,j,k)*strn%dyy(i,j,k) &
                                           + strn%dxx(i,j,k)*strn%dyy(i,j,k) &
                                           + strn%dxy(i,j,k)*strn%dxy(i,j,k) &
                                           + strn%dxz(i,j,k)*strn%dxz(i,j,k) &
                                           + strn%dyz(i,j,k)*strn%dyz(i,j,k) )
                    
                    if (strn%de(i,j,k) .gt. de_max) strn%de(i,j,k) = de_max 

                    ! Calculate the horizontal divergence too 
                    strn%ddiv(i,j,k) = strn%dxx(i,j,k) + strn%dyy(i,j,k) 

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


            end if 

        end do
        end do


        ! === Also calculate vertically averaged strain rate tensor ===
        
        ! Get the 2D average of strain rate in case it is needed 
        strn2D%dxx     = calc_vertical_integrated_2D(strn%dxx, zeta_aa)
        strn2D%dyy     = calc_vertical_integrated_2D(strn%dyy, zeta_aa)
        strn2D%dxy     = calc_vertical_integrated_2D(strn%dxy, zeta_aa)
        strn2D%dxz     = calc_vertical_integrated_2D(strn%dxz, zeta_aa)
        strn2D%dyz     = calc_vertical_integrated_2D(strn%dyz, zeta_aa)
        strn2D%ddiv    = calc_vertical_integrated_2D(strn%ddiv,zeta_aa)
        strn2D%de      = calc_vertical_integrated_2D(strn%de,  zeta_aa)
        strn2D%f_shear = calc_vertical_integrated_2D(strn%f_shear,zeta_aa) 
        
        return 

    end subroutine calc_strain_rate_tensor

    subroutine calc_strain_rate_tensor_2D(strn2D,ux_bar,uy_bar,H_ice,f_ice,f_grnd,dx,dy,de_max,ATT_bar,n_glen)
        ! Calculate the 2D (vertically averaged) strain rate tensor,
        ! assuming a constant vertical velocity profile. 

        implicit none

        type(strain_2D_class), intent(OUT) :: strn2D            ! [yr^-1] Strain rate tensor
        real(wp), intent(IN) :: ux_bar(:,:)                     ! [m/yr] Vertically averaged vel., x
        real(wp), intent(IN) :: uy_bar(:,:)                     ! [m/yr] Vertically averaged vel., y
        real(wp), intent(IN) :: H_ice(:,:)                      ! [m] Ice thickness 
        real(wp), intent(IN) :: f_ice(:,:)
        real(wp), intent(IN) :: f_grnd(:,:)
        real(wp), intent(IN) :: dx, dy                          ! [m] Resolution
        real(wp), intent(IN) :: de_max                          ! [yr^-1] Maximum allowed effective strain rate
        real(wp), intent(IN) :: ATT_bar(:,:) 
        real(wp), intent(IN) :: n_glen 

        ! Local variables
        integer  :: i, j, k
        integer  :: im1, ip1, jm1, jp1
        integer  :: nx, ny
        real(wp) :: dvdx, dudy
        real(wp) :: ux_aa, uy_aa 
        logical  :: is_margin 

        real(wp) :: wt 

        real(wp) :: ddiv_free, dxx_free, dyy_free  

        nx = size(ux_bar,1)
        ny = size(ux_bar,2)

        strn2D%dxx      = 0.0 
        strn2D%dyy      = 0.0 
        strn2D%dxy      = 0.0 
        strn2D%dxz      = 0.0       ! Always zero in this case
        strn2D%dyz      = 0.0       ! Always zero in this case
        strn2D%de       = 0.0 
        strn2D%f_shear  = 0.0       ! Always zero in this case

        do j = 1, ny
        do i = 1, nx
            
            ! Get neighbor indices
            im1 = max(i-1,1) 
            ip1 = min(i+1,nx) 
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny) 
            
            if (H_ice(i,j) .gt. 0.0) then 
                ! Grounded or floating ice, calculate strain rate here (aa-nodes)

                ! aa-nodes
                strn2D%dxx(i,j) = (ux_bar(i,j) - ux_bar(im1,j))/dx
                strn2D%dyy(i,j) = (uy_bar(i,j) - uy_bar(i,jm1))/dy

                ! Calculation of cross terms on central aa-nodes (symmetrical results)
                dudy = ((ux_bar(i,jp1)   - ux_bar(i,jm1))    &
                      + (ux_bar(im1,jp1) - ux_bar(im1,jm1))) /(4.0*dy)
                dvdx = ((uy_bar(ip1,j)   - uy_bar(im1,j))    &
                      + (uy_bar(ip1,jm1) - uy_bar(im1,jm1))) /(4.0*dx)

                strn2D%dxy(i,j) = 0.5*(dudy+dvdx)

                ! Calculate the effective strain rate from the trace 
                strn2D%de(i,j) = sqrt( strn2D%dxx(i,j)**2 + strn2D%dyy(i,j)**2 &
                            + strn2D%dxx(i,j)*strn2D%dyy(i,j) + strn2D%dxy(i,j)**2 )
                
if (.FALSE.) then
    ! ajr: Extrapolating to ice-free and partially ice-free neighbors, 
    ! as done further below, is more stable and convincing than 
    ! imposing the free-spreading strain rate. So this section is disabled. 

                ! Also estimate the free-spreading strain rate (Pollard et al., 2015, EPSL, Eq. B2.b)
                ! ddiv = A*(rho*g*h/4)^n = dxx + dyy
                ! assume equal spreading in both directions:
                ! dxx = dyy; ddiv = 2*dxx
                ! dxx = ddiv/2
                ddiv_free = ATT_bar(i,j) * (0.25*rho_ice*g*H_ice(i,j))**n_glen
                ! dxx_free  = ddiv_free / 2.0
                ! dyy_free  = dxx_free 
                if ( abs(0.5*ux_bar(i,j)+ux_bar(im1,j)) &
                      .gt. abs(0.5*uy_bar(i,j)+uy_bar(i,jm1)) ) then 
                    dxx_free  = ddiv_free
                    dyy_free  = 0.0 
                else 
                    dxx_free  = 0.0
                    dyy_free  = ddiv_free
                end if 

                is_margin = ( H_ice(i,j) .gt. 0.0_wp .and. f_ice(i,j) .lt. 1.0_wp .and. &
                    count([H_ice(im1,j),H_ice(ip1,j),H_ice(i,jm1),H_ice(i,jp1)] .eq. 0.0_wp) .gt. 0 )
                
                ! For partially covered grid cells at margin, or floating points
                ! ensure effective strain rate is not larger than free-spreading strain rate
                if (is_margin .or. &
                     (f_grnd(i,j) .eq. 0.0 .and. strn2D%de(i,j) .gt. ddiv_free) ) then 
                    ! Overwrite above value and impose free-spreading strain 
                    strn2D%ddiv(i,j) = ddiv_free 
                    strn2D%dxx(i,j)  = dxx_free 
                    strn2D%dyy(i,j)  = dyy_free 
                    strn2D%dxy(i,j)  = 0.0

                    strn2D%de(i,j) = sqrt( strn2D%dxx(i,j)**2 + strn2D%dyy(i,j)**2 &
                            + strn2D%dxx(i,j)*strn2D%dyy(i,j) + strn2D%dxy(i,j)**2 )
                 
                end if 

end if 



            end if 

        end do
        end do

        ! === Extrapolate to ice-free neighbors === 
        ! (in case ice gets advected there)
        do j=1, ny
        do i=1, nx

            ! Get neighbor indices
            im1 = max(i-1,1) 
            ip1 = min(i+1,nx) 
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny) 
                
            if ( (H_ice(i,j) .eq. 0.0 .or. f_ice(i,j) .lt. 1.0) .and. &
                count([H_ice(im1,j),H_ice(ip1,j),H_ice(i,jm1),H_ice(i,jp1)] .gt. 0.0_wp) .gt. 0 ) then 
                ! Ice-free (or partially ice-free) with ice-covered neighbors


                strn2D%dxx(i,j) = 0.0 
                strn2D%dyy(i,j) = 0.0
                strn2D%dxy(i,j) = 0.0
                strn2D%dxz(i,j) = 0.0 
                strn2D%dyz(i,j) = 0.0

                wt = 0.0 

                if (H_ice(im1,j) .gt. 0.0 .and. f_ice(im1,j) .eq. 1.0) then 
                    strn2D%dxx(i,j) = strn2D%dxx(i,j) + strn2D%dxx(im1,j)
                    strn2D%dyy(i,j) = strn2D%dyy(i,j) + strn2D%dyy(im1,j)
                    strn2D%dxy(i,j) = strn2D%dxy(i,j) + strn2D%dxy(im1,j)
                    wt = wt + 1.0 
                end if 
                if (H_ice(ip1,j) .gt. 0.0 .and. f_ice(ip1,j) .eq. 1.0) then 
                    strn2D%dxx(i,j) = strn2D%dxx(i,j) + strn2D%dxx(ip1,j)
                    strn2D%dyy(i,j) = strn2D%dyy(i,j) + strn2D%dyy(ip1,j)
                    strn2D%dxy(i,j) = strn2D%dxy(i,j) + strn2D%dxy(ip1,j)
                    wt = wt + 1.0 
                end if
                if (H_ice(i,jm1) .gt. 0.0 .and. f_ice(i,jm1) .eq. 1.0) then 
                    strn2D%dxx(i,j) = strn2D%dxx(i,j) + strn2D%dxx(i,jm1)
                    strn2D%dyy(i,j) = strn2D%dyy(i,j) + strn2D%dyy(i,jm1)
                    strn2D%dxy(i,j) = strn2D%dxy(i,j) + strn2D%dxy(i,jm1)
                    wt = wt + 1.0 
                end if
                if (H_ice(i,jp1) .gt. 0.0 .and. f_ice(i,jp1) .eq. 1.0) then 
                    strn2D%dxx(i,j) = strn2D%dxx(i,j) + strn2D%dxx(i,jp1)
                    strn2D%dyy(i,j) = strn2D%dyy(i,j) + strn2D%dyy(i,jp1)
                    strn2D%dxy(i,j) = strn2D%dxy(i,j) + strn2D%dxy(i,jp1)
                    wt = wt + 1.0 
                end if

                if (wt .gt. 0.0) then 
                    strn2D%dxx(i,j) = strn2D%dxx(i,j) / wt
                    strn2D%dyy(i,j) = strn2D%dyy(i,j) / wt
                    strn2D%dxy(i,j) = strn2D%dxy(i,j) / wt
                end if 

                ! ====== Finished calculating individual strain rate terms ====== 
                
                strn2D%de(i,j) =  sqrt(  strn2D%dxx(i,j)*strn2D%dxx(i,j) &
                                       + strn2D%dyy(i,j)*strn2D%dyy(i,j) &
                                       + strn2D%dxx(i,j)*strn2D%dyy(i,j) &
                                       + strn2D%dxy(i,j)*strn2D%dxy(i,j) &
                                       + strn2D%dxz(i,j)*strn2D%dxz(i,j) &
                                       + strn2D%dyz(i,j)*strn2D%dyz(i,j) )
                
                if (strn2D%de(i,j) .gt. de_max) strn2D%de(i,j) = de_max 

                ! Calculate the horizontal divergence too 
                strn2D%ddiv(i,j) = strn2D%dxx(i,j) + strn2D%dyy(i,j) 

                ! No shearing estimated from 2D components
                strn2D%f_shear(i,j) = 0.0

            end if 

        end do
        end do

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
        call calc_stress_eigen_values(strs2D%teig1,strs2D%teig2, &
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
        call calc_stress_eigen_values(strs2D%teig1,strs2D%teig2, &
                                    strs2D%txx,strs2D%tyy,strs2D%txy)

        return 

    end subroutine calc_stress_tensor_2D
    
    elemental subroutine calc_stress_eigen_values(teig1,teig2,txx,tyy,txy)
        ! Calculate the first two eigenvectors of 2D deviatoric stress tensor 

        implicit none

        real(wp), intent(OUT) :: teig1              ! [Pa] Eigenvalue 1
        real(wp), intent(OUT) :: teig2              ! [Pa] Eigenvalue 2
        real(wp), intent(IN)  :: txx
        real(wp), intent(IN)  :: tyy
        real(wp), intent(IN)  :: txy

        ! Local variables
        real(wp) :: a, b, c, root
        real(wp) :: lambda1, lambda2 

        ! compute the eigenvalues of the vertically integrated stress tensor
        a = 1.0_wp
        b = -(txx + tyy)
        c = txx*tyy - txy*txy
        if (b*b - 4.0_wp*a*c > 0.0_wp) then   ! two real eigenvalues
            root = sqrt(b*b - 4.0_wp*a*c)
            lambda1 = (-b + root) / (2.0_wp*a)
            lambda2 = (-b - root) / (2.0_wp*a)
            if (lambda1 > lambda2) then
                teig1 = lambda1
                teig2 = lambda2
            else
                teig1 = lambda2
                teig2 = lambda1
            end if
        else 
            ! No eigenvalues, set to zero 
            teig1 = 0.0_wp 
            teig2 = 0.0_wp 
        end if  ! b^2 - 4ac > 0

        return

    end subroutine calc_stress_eigen_values


end module deformation



