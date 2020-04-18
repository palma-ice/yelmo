

module deformation 
    ! This module contains functions related to deformation calculations:
    ! Arrenius function
    ! Flow law
    ! viscosity
    ! strain rate??

    ! Note: 3D arrays defined such that first index (k=1) == base, and max index (k=nk) == surface 

    use yelmo_defs,  only : sp, dp, prec, tol_underflow, T0, strain_2D_class, strain_3D_class
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
    public :: calc_strain_rate_2D
    public :: calc_strain_rate_3D

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
        enh         = f_shear*enh_shear   + (1.0-f_shear)*enh_ssa_tmp
        
        return 

    end function define_enhancement_factor_2D

    function calc_viscosity_glen(de,ATT,n_glen,visc_min) result(visc)
        ! Calculate viscosity based on Glen's flow law 
        ! ATT [a^-1 Pa^-3] is the "depth dependent ice stiffness parameter based on
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
        real(prec), parameter :: de_0 = 1.0e-6    ! [a^-1] Bueler and Brown (2009), Eq. 26

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
    
    subroutine calc_strain_rate_2D(strn2D,ux_bar,uy_bar,dx,dy)
    
        implicit none

        real(prec), intent(IN) :: ux_bar(:,:)    ! [m/a] Vertically integrated vel., x
        real(prec), intent(IN) :: uy_bar(:,:)    ! [m/a] Vertically integrated vel., y
        real(prec), intent(IN) :: dx, dy         ! [m] Resolution
        type(strain_2D_class)  :: strn2D         ! [a^-1] Strain rate tensor

        ! Local variables
        integer :: i, j, k, nx, ny, nz
        real(prec) :: dvdx, dudy
        
        nx = size(ux_bar,1)
        ny = size(ux_bar,2)

        strn2D%dxx = 0.0 
        strn2D%dyy = 0.0 
        strn2D%dxy = 0.0 
        strn2D%de  = 0.0 

        do j = 2, ny-1
        do i = 2, nx-1
            
            ! Strain rates should be returned on central Aa nodes
            
            ! Note: in Grisli-version8-chris, the cross terms dudy/dvdx were
            ! calculated on the 'noeud mineur': Ab nodes. However, this results
            ! in asymmetrical de and visc fields in EISMINT tests. Using the older
            ! version that caculates cross terms on the Aa nodes results
            ! in symmetrical de and visc fields. 
            ! Grisli-version8-chris also did not employ the cross terms in the calulation
            ! of de, which is perhaps why no asymmetery was noticed in the results? 

            ! Aa node
            strn2D%dxx(i,j) = (ux_bar(i,j) - ux_bar(i-1,j))/dx
            strn2D%dyy(i,j) = (uy_bar(i,j) - uy_bar(i,j-1))/dy

            ! Calculation of cross terms on central Aa nodes (symmetrical results)
            dudy = ((ux_bar(i,j+1) - ux_bar(i,j-1))    &
                  + (ux_bar(i-1,j+1)   - ux_bar(i-1,j-1))) /(4.0*dy)
            dvdx = ((uy_bar(i+1,j) - uy_bar(i-1,j))    &
                  + (uy_bar(i+1,j-1)   - uy_bar(i-1,j-1)))/(4.0*dx)

            strn2D%dxy(i,j) = 0.5*(dudy+dvdx)

            ! Calculate the effective strain rate from the trace 
            strn2D%de(i,j) = strn2D%dxx(i,j)**2+strn2D%dyy(i,j)**2+strn2D%dxx(i,j)*strn2D%dyy(i,j)+strn2D%dxy(i,j)**2
            strn2D%de(i,j) = strn2D%de(i,j)**0.5
            
        end do
        end do
        
        return
        
    end subroutine calc_strain_rate_2D
    
    subroutine calc_strain_rate_3D(strn,vx,vy,vz,H_ice,f_grnd,zeta_aa,zeta_ac,dx)

        implicit none
         
        type(strain_3D_class), intent(INOUT) :: strn  ! [a^-1] 
        real(prec), intent(IN) :: vx(:,:,:) 
        real(prec), intent(IN) :: vy(:,:,:) 
        real(prec), intent(IN) :: vz(:,:,:) 
        real(prec), intent(IN) :: H_ice(:,:) 
        real(prec), intent(IN) :: f_grnd(:,:) 
        real(prec), intent(IN) :: zeta_aa(:) 
        real(prec), intent(IN) :: zeta_ac(:) 
        real(prec), intent(IN) :: dx 

        ! Local variables 
        real(prec) :: dy  

        dy = dx 

        ! Calculate 3D strain rate
        call calc_dxyz(strn,vx,vy,vz,H_ice,f_grnd,zeta_aa,zeta_ac,dx,dy)
        
        return 

    end subroutine calc_strain_rate_3D

    subroutine calc_dxyz(strn, vx, vy, vz, H_ice, f_grnd, zeta_aa, zeta_ac, dx, dy)
        ! -------------------------------------------------------------------------------
        !  Computation of all components of the strain-rate tensor, the full
        !  effective strain rate and the shear fraction.
        !  Alexander Robinson: Adapted from sicopolis5-dev
        ! ------------------------------------------------------------------------------

        ! Note: vx, vy are staggered on ac-nodes in the horizontal, but are on the zeta_aa nodes (ie layer-centered)
        ! in the vertical. vz is centered on aa-nodes in the horizontal, but staggered on zeta_ac nodes
        ! in the vertical. 

        implicit none
        
        type(strain_3D_class), intent(INOUT) :: strn                ! [a^-1] on aa-nodes (3D)
        real(prec), intent(IN) :: vx(:,:,:)                         ! nx,ny,nz_aa
        real(prec), intent(IN) :: vy(:,:,:)                         ! nx,ny,nz_aa
        real(prec), intent(IN) :: vz(:,:,:)                         ! nx,ny,nz_ac
        real(prec), intent(IN) :: H_ice(:,:)
        real(prec), intent(IN) :: f_grnd(:,:)
        real(prec), intent(IN) :: zeta_aa(:) 
        real(prec), intent(IN) :: zeta_ac(:) 
        real(prec), intent(IN) :: dx, dy

        ! Local variables 
        integer    :: i, j, k
        integer    :: im1, ip1, jm1, jp1 
        integer    :: nx, ny, nz_aa, nz_ac  
        real(prec) :: dxi, deta, dzeta
        real(prec) :: dx_inv, dy_inv
        real(prec) :: H_ice_inv
        real(prec) :: abs_v_ssa_inv, nx1, ny1
        real(prec) :: shear_x_help, shear_y_help
        real(prec) :: f_shear_help
        real(prec), allocatable :: fact_x(:,:), fact_y(:,:)
        real(prec), allocatable :: fact_z(:)
        real(prec), allocatable :: lxy(:), lyx(:), lxz(:), lzx(:), lyz(:), lzy(:)
        real(prec), allocatable :: shear_squared(:)
        
        real(prec), parameter :: epsi   = 1e-5
        real(prec), parameter :: de_max = 0.5    ! [a^-1] 

        ! Determine sizes and allocate local variables 
        nx    = size(vx,1)
        ny    = size(vx,2)
        nz_aa = size(zeta_aa,1)
        nz_ac = size(zeta_ac,1)
        
        allocate(fact_x(nx,ny))
        allocate(fact_y(nx,ny))
        allocate(fact_z(nz_aa))
        allocate(lxy(nz_aa))
        allocate(lyx(nz_aa))
        allocate(lxz(nz_aa))
        allocate(lzx(nz_aa))
        allocate(lyz(nz_aa))
        allocate(lzy(nz_aa))
        allocate(shear_squared(nz_aa))

        ! Change arguments to local (sicopolis) variable names 
        dxi  = dx 
        deta = dy 

        !-------- Term abbreviations --------

        dx_inv = 1.0/dx
        dy_inv = 1.0/dy

        fact_x   = dx_inv
        fact_y   = dy_inv

        fact_z(1) = 1.0/(zeta_aa(2)-zeta_aa(1))
        do k = 2, nz_aa-1 
            fact_z(k) = 1.0/(zeta_aa(k+1)-zeta_aa(k-1))
        end do
        fact_z(nz_aa) = 1.0/(zeta_aa(nz_aa)-zeta_aa(nz_aa-1))

        !-------- Initialisation --------

        strn%dxx          = 0.0
        strn%dyy          = 0.0
        strn%dxy          = 0.0
        strn%dxz          = 0.0
        strn%dyz          = 0.0
        strn%dzz          = 0.0
        strn%de           = 0.0
        strn%f_shear      = 0.0

        !-------- Computation --------

        !!!$omp parallel do
        do j=1, ny
        do i=1, nx

            im1 = max(i-1,1) 
            ip1 = min(i+1,nx) 
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny) 
            
            if (H_ice(i,j) .eq. 0.0) then 
                ! Ice-free, set all strain rate terms to zero 

                strn%dxx(i,j,:)     = 0.0
                strn%dyy(i,j,:)     = 0.0
                strn%dxy(i,j,:)     = 0.0
                strn%dxz(i,j,:)     = 0.0
                strn%dyz(i,j,:)     = 0.0
                strn%dzz(i,j,:)     = 0.0
                strn%de(i,j,:)      = 0.0
                strn%f_shear(i,j,:) = 0.0

            else 
                ! Grounded or floating ice present

                H_ice_inv = 1.0/(abs(H_ice(i,j))+epsi)

                ! ====== Horizontal gradients (dxx,dyy,lxy,lyx) ======

                ! Full column from base to surface 
                do k = 1, nz_aa 
                    strn%dxx(i,j,k) = (vx(i,j,k)-vx(im1,j,k))*fact_x(i,j)
                    strn%dyy(i,j,k) = (vy(i,j,k)-vy(i,jm1,k))*fact_y(i,j)

                    ! Avoid underflows 
                    if (abs(strn%dxx(i,j,k)) .lt. tol_underflow) strn%dxx(i,j,k) = 0.0 
                    if (abs(strn%dyy(i,j,k)) .lt. tol_underflow) strn%dyy(i,j,k) = 0.0 
                    
                    lxy(k) = (  (vx(i,jp1,k)+vx(im1,jp1,k))   &
                              - (vx(i,jm1,k)+vx(im1,jm1,k)) ) &
                             *0.25*fact_y(i,j)

                    lyx(k) = (  (vy(ip1,j,k)+vy(ip1,jm1,k))   &
                              - (vy(im1,j,k)+vy(im1,jm1,k)) ) &
                             *0.25*fact_x(i,j)
                end do 

                ! ====== Vertical cross terms (lzx,lzy) ====== 

                ! Basal layer
                k = 1
                lzx(k) = (vz(ip1,j,k)-vz(im1,j,k))*0.5*fact_x(i,j)
                lzy(k) = (vz(i,jp1,k)-vz(i,jm1,k))*0.5*fact_y(i,j)

                ! Intermediate layers 
                do k = 2, nz_aa-1
                    lzx(k) = (  (vz(ip1,j,k)+vz(ip1,j,k-1))   &
                              - (vz(im1,j,k)+vz(im1,j,k-1)) ) &
                             *0.25*fact_x(i,j)

                    lzy(k) = (  (vz(i,jp1,k)+vz(i,jp1,k-1))   &
                              - (vz(i,jm1,k)+vz(i,jm1,k-1)) ) &
                             *0.25*fact_y(i,j)
                end do 

                ! Surface layer 
                k = nz_aa
                lzx(k) = (vz(ip1,j,k-1)-vz(im1,j,k-1))*0.5*fact_x(i,j)
                lzy(k) = (vz(i,jp1,k-1)-vz(i,jm1,k-1))*0.5*fact_y(i,j)


                ! ====== Shear terms (lxz,lyz) ====== 

                ! Basal layer
                k = 1

                ! Gradient from first aa-node above base to base 
                lxz(k) = (  (vx(i,j,k+1)+vx(im1,j,k+1))   &
                          - (vx(i,j,  k)+vx(im1,j,  k)) ) &
                         *0.5*fact_z(k)*H_ice_inv

                lyz(k) = (  (vy(i,j,k+1)+vy(i,jm1,k+1))   &
                          - (vy(i,j,  k)+vy(i,jm1,  k)) ) &
                         *0.5*fact_z(k)*H_ice_inv

                ! Intermediate layers
                do k = 2, nz_aa-1
                     
                    ! Gradient from aa-node above to aa-node below
                    lxz(k) = (  (vx(i,j,k+1)+vx(im1,j,k+1))   &
                              - (vx(i,j,k-1)+vx(im1,j,k-1)) ) &
                             *0.5*fact_z(k)*H_ice_inv

                    lyz(k) = (  (vy(i,j,k+1)+vy(i,jm1,k+1))   &
                              - (vy(i,j,k-1)+vy(i,jm1,k-1)) ) &
                             *0.5*fact_z(k)*H_ice_inv 

                end do
                
                ! Surface layer
                ! Gradient from surface to first aa-node below surface 
                k = nz_aa
                lxz(k) = (  (vx(i,j,  k)+vx(im1,j,  k))   &
                          - (vx(i,j,k-1)+vx(im1,j,k-1)) ) &
                         *0.5*fact_z(k)*H_ice_inv

                lyz(k) = (  (vy(i,j,  k)+vy(i,jm1,  k))   &
                          - (vy(i,j,k-1)+vy(i,jm1,k-1)) ) &
                         *0.5*fact_z(k)*H_ice_inv

                ! ====== Calculate column of cross terms from intermediate values (dxy,dxz,dyz) ====== 

                strn%dxy(i,j,:) = 0.5*(lxy+lyx)
                strn%dxz(i,j,:) = 0.5*(lxz+lzx)
                strn%dyz(i,j,:) = 0.5*(lyz+lzy)

                ! Avoid underflows 
                where (abs(strn%dxy(i,j,:)) .lt. tol_underflow) strn%dxy(i,j,:) = 0.0 
                where (abs(strn%dxz(i,j,:)) .lt. tol_underflow) strn%dxz(i,j,:) = 0.0 
                where (abs(strn%dyz(i,j,:)) .lt. tol_underflow) strn%dyz(i,j,:) = 0.0 
    
                ! ====== Finished calculating individual strain rate terms ====== 
                
                ! Calculate the shear-based strain, stretching and the shear-fraction
                do k = 1, nz_aa

                    shear_squared(k) =   strn%dxz(i,j,k)*strn%dxz(i,j,k) &
                                       + strn%dyz(i,j,k)*strn%dyz(i,j,k)

                    strn%de(i,j,k)    =  sqrt(   strn%dxx(i,j,k)*strn%dxx(i,j,k) &
                                               + strn%dyy(i,j,k)*strn%dyy(i,j,k) &
                                               + strn%dxx(i,j,k)*strn%dyy(i,j,k) &
                                               + strn%dxy(i,j,k)*strn%dxy(i,j,k) &
                                               + shear_squared(k) )
                        
                    if (strn%de(i,j,k) .gt. de_max) strn%de(i,j,k) = de_max 

                    ! Note: Using only the below should be equivalent to applying
                    ! the SIA approximation to calculate `de`
                    !strn%de(i,j,k)    =  sqrt( shear_squared(k) )

                    if (strn%de(i,j,k) .gt. 0.0) then 
                        strn%f_shear(i,j,k) = sqrt(shear_squared(k))/strn%de(i,j,k)
                    else 
                        strn%f_shear(i,j,k) = 1.0   ! Shearing by default for low strain rates
                    end if 

                end do

                !  ------ Modification of the shear fraction for floating ice (ice shelves)

                if (f_grnd(i,j) .eq. 0.0) then 
                    strn%f_shear(i,j,:) = 0.0    ! Assume ice shelf is only stretching, no shear 
                end if 

                !  ------ Constrain the shear fraction to reasonable [0,1] interval

                strn%f_shear(i,j,:) = min(max(strn%f_shear(i,j,:), 0.0), 1.0)


            end if ! ice-free or ice-covered 

        end do
        end do
        !!!$omp end parallel do

        return 

    end subroutine calc_dxyz

!     subroutine calc_stress_3D(strss,visc,strn)
!         ! Note: this is not used explicitly in Yelmo,
!         ! it's just here for completeness. For SIA/SSA
!         ! it is enough to use strain everywhere. 
        
!         implicit none 

!         type(stress_3D_class), intent(INOUT) :: strss 
!         real(prec),            intent(IN)    :: visc(:,:,:)
!         type(strain_3D_class), intent(IN)    :: strn 

!         strss%txx = 2.0*visc*strn%dxx
!         strss%tyy = 2.0*visc*strn%dyy
!         strss%txy = 2.0*visc*strn%dxy
!         strss%txz = 2.0*visc*strn%dxz
!         strss%tyz = 2.0*visc*strn%dyz
!         strss%tzz = 2.0*visc*strn%dzz
        
!         return 

!     end subroutine calc_stress_3D
    
end module deformation 



