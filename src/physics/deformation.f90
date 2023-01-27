

module deformation 
    ! This module contains functions related to deformation calculations:
    ! Arrenius function
    ! Flow law
    ! viscosity
    ! strain rate??

    ! Note: 3D arrays defined such that first index (k=1) == base, and max index (k=nk) == surface 

    use yelmo_defs,  only : sp, dp, wp, prec, TOL_UNDERFLOW, T0, rho_ice, g, &
                        jacobian_3D_class, strain_2D_class, strain_3D_class, stress_2D_class, stress_3D_class
    use yelmo_tools, only : get_neighbor_indices, &
                    calc_vertical_integrated_2D, integrate_trapezoid1D_1D, integrate_trapezoid1D_pt, &
                    acx_to_nodes, acy_to_nodes, &
                    stagger_nodes_aa_ab_ice, stagger_nodes_acx_ab_ice, stagger_nodes_acy_ab_ice, &
                    staggerdiff_nodes_acx_ab_ice, staggerdiff_nodes_acy_ab_ice, &
                    staggerdiffcross_nodes_acx_ab_ice, staggerdiffcross_nodes_acy_ab_ice, &
                    staggerdiff_nodes_acz_dx_ab_ice, staggerdiff_nodes_acz_dy_ab_ice
                    
    use grid_calcs 

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

    public :: calc_jacobian_vel_3D
    public :: calc_strain_rate_tensor_jac

    public :: calc_strain_rate_tensor
    public :: calc_strain_rate_tensor_2D
    public :: calc_strain_rate_horizontal_2D
    public :: calc_stress_tensor 
    public :: calc_stress_tensor_2D
    public :: calc_2D_eigen_values
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
    
    subroutine calc_strain_rate_tensor_2D_ab(strn2D,vx,vy,H_ice,f_ice,f_grnd,dx,de_max,n_glen)
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

        ! Finally, calculate the first two eigenvectors for 2D strain rate tensor 
        call calc_2D_eigen_values(strn2D%eps_eig_1,strn2D%eps_eig_2, &
                                    strn2D%dxx,strn2D%dyy,strn2D%dxy)

        return
        
    end subroutine calc_strain_rate_tensor_2D_ab
    




        subroutine calc_jacobian_vel_3D(jvel, ux, uy, uz, H_ice, f_ice, f_grnd, dzsdx, dzsdy,  &
                                                dzbdx, dzbdy, zeta_aa, zeta_ac, dx, dy, boundaries)

        ! -------------------------------------------------------------------------------
        !  Computation of all components of the Jacobian of the velocity vector:
        !  (ux,uy,uz) == (u,v,w)
        !
        ! ------------------------------------------------------------------------------

        ! Note: vx, vy are staggered on ac-nodes in the horizontal, but are on the zeta_aa nodes (ie layer-centered)
        ! in the vertical. vz is centered on aa-nodes in the horizontal, but staggered on zeta_ac nodes
        ! in the vertical. 

        ! All tensor components are calculated in the same location as the velocity components. 
        
        implicit none
        
        type(jacobian_3D_class), intent(INOUT) :: jvel          ! [yr^-1] on ac-nodes (3D)
        real(wp), intent(IN) :: ux(:,:,:)                       ! nx,ny,nz_aa
        real(wp), intent(IN) :: uy(:,:,:)                       ! nx,ny,nz_aa
        real(wp), intent(IN) :: uz(:,:,:)                       ! nx,ny,nz_ac
        real(wp), intent(IN) :: H_ice(:,:)
        real(wp), intent(IN) :: f_ice(:,:)
        real(wp), intent(IN) :: f_grnd(:,:)
        real(wp), intent(IN) :: dzsdx(:,:) 
        real(wp), intent(IN) :: dzsdy(:,:) 
        real(wp), intent(IN) :: dzbdx(:,:) 
        real(wp), intent(IN) :: dzbdy(:,:) 
        real(wp), intent(IN) :: zeta_aa(:) 
        real(wp), intent(IN) :: zeta_ac(:) 
        real(wp), intent(IN) :: dx
        real(wp), intent(IN) :: dy
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer  :: i, j, k
        integer  :: im1, ip1, jm1, jp1 
        integer  :: nx, ny, nz_aa, nz_ac
        real(wp) :: H_now, H_now_acx, H_now_acy 
        real(wp) :: dzbdx_acy, dzbdy_acx, dzsdx_acy, dzsdy_acx
        real(wp) :: dzbdx_aa, dzbdy_aa, dzsdx_aa, dzsdy_aa
        real(wp) :: c_x, c_y, c_x_acy, c_y_acx
        real(wp) :: h1, h2 


        ! Determine sizes and allocate local variables 
        nx    = size(ux,1)
        ny    = size(ux,2)
        nz_aa = size(zeta_aa,1)
        nz_ac = size(zeta_ac,1)
        
        !-------- Initialisation --------

        jvel%dxx          = 0.0_wp
        jvel%dxy          = 0.0_wp
        jvel%dxz          = 0.0_wp
        jvel%dyx          = 0.0_wp
        jvel%dyy          = 0.0_wp
        jvel%dyz          = 0.0_wp
        jvel%dzx          = 0.0_wp
        jvel%dzy          = 0.0_wp
        jvel%dzz          = 0.0_wp

        !-------- Computation --------

        ! Step 1: Calculate all vertical derivatives, some of which are used 
        ! as correction terms for the horizontal derivatives w.r.t. sigma-coordinate transformation. 

        !$omp parallel do
        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)


            if (f_ice(i,j) .gt. 0.0) then 
                ! Ice is present at this grid point
                ! Include fractional points, since these derivatives are calculated
                ! for ac-nodes. For fractional points, ac-node may border a fully
                ! ice-covered grid point. 

                H_now     = H_ice(i,j)

                ! Get staggered ice thickness too 
                if (f_ice(i,j) .eq. 1.0 .and. f_ice(ip1,j) .lt. 1.0) then 
                    H_now_acx = H_ice(i,j)
                else if (f_ice(i,j) .lt. 1.0 .and. f_ice(ip1,j) .eq. 1.0) then 
                    H_now_acx = H_ice(ip1,j)
                else
                    H_now_acx = 0.5*(H_ice(i,j)+H_ice(ip1,j))
                end if

                if (f_ice(i,j) .eq. 1.0 .and. f_ice(i,jp1) .lt. 1.0) then 
                    H_now_acy = H_ice(i,j)
                else if (f_ice(i,j) .lt. 1.0 .and. f_ice(i,jp1) .eq. 1.0) then 
                    H_now_acy = H_ice(i,jp1)
                else
                    H_now_acy = 0.5*(H_ice(i,j)+H_ice(i,jp1))
                end if


                ! Calculate dxz, dyz on aa-nodes vertically, ac-nodes horizontally

                ! Bottom layer - upwind derivative
                k = 1
                jvel%dxz(i,j,k) = (ux(i,j,k+1)-ux(i,j,k)) / (H_now_acx*(zeta_aa(k+1)-zeta_aa(k)))
                jvel%dyz(i,j,k) = (uy(i,j,k+1)-uy(i,j,k)) / (H_now_acy*(zeta_aa(k+1)-zeta_aa(k)))

                ! Intermediate layers
                do k = 2, nz_aa-1 

                    ! Use centered, 1st order derivative for uneven layers 
                    ! see "Finite Difference Formulae for Unequal Sub-Intervals Using Lagrange’s Interpolation Formula"
                    ! by Singh and Bhadauria
                    ! http://www.m-hikari.com/ijma/ijma-password-2009/ijma-password17-20-2009/bhadauriaIJMA17-20-2009.pdf

                    !h1 = H_now_acx*(zeta_aa(k)-zeta_aa(k-1))
                    !h2 = H_now_acx*(zeta_aa(k+1)-zeta_aa(k))
                    !jvel%dxz(i,j,k) = -h2/(h1*(h1+h2))*ux(i,j,k-1) - (h1-h2)/(h1*h2)*ux(i,j,k) + h1/(h2*h1+h2)*ux(i,j,k+1)

                    !h1 = H_now_acy*(zeta_aa(k)-zeta_aa(k-1))
                    !h2 = H_now_acy*(zeta_aa(k+1)-zeta_aa(k))
                    !jvel%dyz(i,j,k) = -h2/(h1*(h1+h2))*uy(i,j,k-1) - (h1-h2)/(h1*h2)*uy(i,j,k) + h1/(h2*h1+h2)*uy(i,j,k+1)
                    
                    ! Simple upwind derivatives
                    jvel%dxz(i,j,k) = (ux(i,j,k+1)-ux(i,j,k)) / (H_now_acx*(zeta_aa(k+1)-zeta_aa(k)))
                    jvel%dyz(i,j,k) = (uy(i,j,k+1)-uy(i,j,k)) / (H_now_acy*(zeta_aa(k+1)-zeta_aa(k)))

                end do 

                ! Top layer - downwind derivative
                k = nz_aa 
                jvel%dxz(i,j,k) = (ux(i,j,k)-ux(i,j,k-1)) / (H_now_acx*(zeta_aa(k)-zeta_aa(k-1)))
                jvel%dyz(i,j,k) = (uy(i,j,k)-uy(i,j,k-1)) / (H_now_acy*(zeta_aa(k)-zeta_aa(k-1)))

            end if

            if (f_ice(i,j) .eq. 1.0) then
                ! Ice present at this point
                ! Derivatives are calculated horizontally for aa-nodes, thus, we
                ! are only concerned with fully ice-covered points now. 

                ! Calculate dzz on ac-nodes vertically, aa-nodes horizontally

                ! Bottom layer - upwind derivative
                k = 1
                jvel%dzz(i,j,k) = (uz(i,j,k+1)-uz(i,j,k)) / (H_now*(zeta_ac(k+1)-zeta_ac(k)))

                ! Intermediate layers
                do k = 2, nz_ac-1 

                    ! Use centered, 1st order derivative for uneven layers 
                    ! see "Finite Difference Formulae for Unequal Sub-Intervals Using Lagrange’s Interpolation Formula"
                    ! by Singh and Bhadauria
                    ! http://www.m-hikari.com/ijma/ijma-password-2009/ijma-password17-20-2009/bhadauriaIJMA17-20-2009.pdf

                    !h1 = H_now*(zeta_ac(k)-zeta_ac(k-1))
                    !h2 = H_now*(zeta_ac(k+1)-zeta_ac(k))
                    !jvel%dzz(i,j,k) = -h2/(h1*(h1+h2))*uz(i,j,k-1) - (h1-h2)/(h1*h2)*uz(i,j,k) + h1/(h2*h1+h2)*uz(i,j,k+1)

                    ! Simple upwind derivative
                    jvel%dzz(i,j,k) = (uz(i,j,k+1)-uz(i,j,k)) / (H_now*(zeta_ac(k+1)-zeta_ac(k)))

                end do 

                ! Top layer - downwind derivative
                k = nz_ac
                jvel%dzz(i,j,k) = (uz(i,j,k)-uz(i,j,k-1)) / (H_now*(zeta_ac(k)-zeta_ac(k-1)))

            end if

        end do
        end do 
        !$omp end parallel do

        ! Set 2: Calculate all horizontal derivatives accounting for correction terms

        !$omp parallel do
        do j = 1, ny 
        do i = 1, nx 
            
            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)


            if (f_ice(i,j) .gt. 0.0) then 
                ! Ice is present at this grid point*
                ! Include fractional points, since these derivatives are calculated
                ! for ac-nodes. For fractional points, ac-node may border a fully
                ! ice-covered grid point. 

                do k = 1, nz_aa 

                    ! === Calculate derivatives , no sigma-correction terms yet ===

                    ! Second-order, centered derivatives
                    jvel%dxx(i,j,k) = (ux(ip1,j,k)-ux(im1,j,k))/(2.0*dx)
                    jvel%dxy(i,j,k) = (ux(i,jp1,k)-ux(i,jm1,k))/(2.0*dy)
                    
                    jvel%dyx(i,j,k) = (uy(ip1,j,k)-uy(im1,j,k))/(2.0*dx)
                    jvel%dyy(i,j,k) = (uy(i,jp1,k)-uy(i,jm1,k))/(2.0*dy)

                    ! Treat special cases of ice-margin points (take upstream/downstream derivatives instead)

                    ! jvel%dxx
                    if (f_ice(i,j) .eq. 1.0 .and. f_ice(ip1,j) .lt. 1.0) then 
                        jvel%dxx(i,j,k) = (ux(i,j,k)-ux(im1,j,k))/dx
                    else if (f_ice(i,j) .lt. 1.0 .and. f_ice(ip1,j) .eq. 1.0) then 
                        jvel%dxx(i,j,k) = (ux(ip1,j,k)-ux(i,j,k))/dx
                    end if 

                    ! jvel%dxy
                    if (f_ice(i,j) .eq. 1.0 .and. f_ice(i,jp1) .lt. 1.0 .and. f_ice(i,jm1) .lt. 1.0) then 
                        jvel%dxy(i,j,k) = 0.0
                    else if (f_ice(i,j) .eq. 1.0 .and. f_ice(i,jp1) .lt. 1.0 .and. f_ice(i,jm1) .eq. 1.0) then 
                        jvel%dxy(i,j,k) = (ux(i,j,k)-ux(i,jm1,k))/dy
                    else if (f_ice(i,j) .eq. 1.0 .and. f_ice(i,jp1) .eq. 1.0 .and. f_ice(i,jm1) .lt. 1.0) then 
                        jvel%dxy(i,j,k) = (ux(i,jp1,k)-ux(i,j,k))/dy
                    end if 

                    ! jvel%dyx
                    if (f_ice(i,j) .eq. 1.0 .and. f_ice(ip1,j) .lt. 1.0 .and. f_ice(im1,j) .lt. 1.0) then 
                        jvel%dyx(i,j,k) = 0.0
                    else if (f_ice(i,j) .eq. 1.0 .and. f_ice(ip1,j) .lt. 1.0 .and. f_ice(im1,j) .eq. 1.0) then 
                        jvel%dyx(i,j,k) = (uy(i,j,k)-uy(im1,j,k))/dx
                    else if (f_ice(i,j) .eq. 1.0 .and. f_ice(ip1,j) .eq. 1.0 .and. f_ice(im1,j) .lt. 1.0) then 
                        jvel%dyx(i,j,k) = (uy(ip1,j,k)-uy(i,j,k))/dx
                    end if 

                    ! jvel%dyy
                    if (f_ice(i,j) .eq. 1.0 .and. f_ice(i,jp1) .lt. 1.0) then 
                        jvel%dyy(i,j,k) = (uy(i,j,k)-uy(i,jm1,k))/dy
                    else if (f_ice(i,j) .lt. 1.0 .and. f_ice(i,jp1) .eq. 1.0) then 
                        jvel%dyy(i,j,k) = (uy(i,jp1,k)-uy(i,j,k))/dy
                    end if 

                    ! === Calculate and apply the sigma-transformation correction terms ===

                    ! First, calculate sigma-coordinate derivative correction factors for this layer
                    ! (Greve and Blatter, 2009, Eqs. 5.131 and 5.132, also shown in 1D with Eq. 5.145)
                    ! Note: these factors can be calculated with different combinations of terms using
                    ! the gradients (dzbdx,dzsdx), (dzbdx,dHidx) or (dHidx,dzsdx), with the same result. 
                    ! (dzbdx,dzsdx) is used for convenience, since dzsdx and dzbdx are needed elsewhere.

                    ! Note that this factor should normally include H_ice in the denominator. It
                    ! is not included here, because the vertical derivative is in z-coordinates
                    ! not sigma-coordinates (zeta), so the H has already been accounted for, 
                    ! e.g.: d/dz = d/(dzeta*H)

                    c_x = - ( (1.0-zeta_aa(k))*dzbdx(i,j) + zeta_aa(k)*dzsdx(i,j))
                    c_y = - ( (1.0-zeta_aa(k))*dzbdy(i,j) + zeta_aa(k)*dzsdy(i,j))

                    ! Get cross correction terms too 
                    dzbdx_acy = 0.25*(dzbdx(i,j)+dzbdx(i,jp1)+dzbdx(im1,j)+dzbdx(im1,jp1))
                    dzsdx_acy = 0.25*(dzsdx(i,j)+dzsdx(i,jp1)+dzsdx(im1,j)+dzsdx(im1,jp1))
                    c_x_acy = - ( (1.0-zeta_aa(k))*dzbdx_acy + zeta_aa(k)*dzsdx_acy)

                    dzbdy_acx = 0.25*(dzbdy(i,j)+dzbdy(ip1,j)+dzbdy(i,jm1)+dzbdy(ip1,jm1))
                    dzsdy_acx = 0.25*(dzsdy(i,j)+dzsdy(ip1,j)+dzsdy(i,jm1)+dzsdy(ip1,jm1))
                    c_y_acx = - ( (1.0-zeta_aa(k))*dzbdy_acx + zeta_aa(k)*dzsdy_acx)

                    ! Apply the correction 

                    jvel%dxx(i,j,k) = jvel%dxx(i,j,k) + c_x*jvel%dxz(i,j,k)
                    jvel%dxy(i,j,k) = jvel%dxy(i,j,k) + c_y_acx*jvel%dxz(i,j,k)
                    
                    jvel%dyx(i,j,k) = jvel%dyx(i,j,k) + c_x_acy*jvel%dyz(i,j,k)
                    jvel%dyy(i,j,k) = jvel%dyy(i,j,k) + c_y*jvel%dyz(i,j,k)
                    
                end do

            end if 

            if (f_ice(i,j) .eq. 1.0) then 
                ! Ice present at this point
                ! Derivatives are calculated horizontally for aa-nodes, thus, we
                ! are only concerned with fully ice-covered points now. 

                ! === Now get horizontal derivatives of uz (ac-nodes vertically, aa-nodes horizontally) ===

                do k = 1, nz_ac

                    ! === Calculate derivatives , no sigma-correction terms yet ===

                    ! Second-order, centered derivatives
                    jvel%dzx(i,j,k) = (uz(ip1,j,k)-uz(im1,j,k))/(2.0*dx)
                    jvel%dzy(i,j,k) = (uz(i,jp1,k)-uz(i,jm1,k))/(2.0*dy)
                    
                    ! Treat special cases of ice-margin points (take upstream/downstream derivatives instead)

                    ! jvel%dzx
                    if (f_ice(ip1,j) .lt. 1.0 .and. f_ice(im1,j) .lt. 1.0) then 
                        jvel%dzx(i,j,k) = 0.0
                    else if (f_ice(ip1,j) .lt. 1.0 .and. f_ice(im1,j) .eq. 1.0) then 
                        jvel%dzx(i,j,k) = (uz(i,j,k)-uz(im1,j,k))/dx
                    else if (f_ice(ip1,j) .eq. 1.0 .and. f_ice(im1,j) .lt. 1.0) then 
                        jvel%dzx(i,j,k) = (uz(ip1,j,k)-uz(i,j,k))/dx
                    end if 

                    ! jvel%dzy
                    if (f_ice(i,jp1) .lt. 1.0 .and. f_ice(i,jm1) .lt. 1.0) then 
                        jvel%dzy(i,j,k) = 0.0
                    else if (f_ice(i,jp1) .lt. 1.0 .and. f_ice(i,jm1) .eq. 1.0) then 
                        jvel%dzy(i,j,k) = (uz(i,j,k)-uz(i,jm1,k))/dy
                    else if (f_ice(i,jp1) .eq. 1.0 .and. f_ice(i,jm1) .lt. 1.0) then 
                        jvel%dzy(i,j,k) = (uz(i,jp1,k)-uz(i,j,k))/dy
                    end if 
                    
                    ! === Calculate and apply the sigma-transformation correction terms ===

                    ! Recalculate correction factors on aa-nodes horizontally, ac-nodes vertically

                    dzbdx_aa = 0.5*(dzbdx(i,j)+dzbdx(im1,j))
                    dzbdy_aa = 0.5*(dzbdy(i,j)+dzbdy(i,jm1))
                    dzsdx_aa = 0.5*(dzsdx(i,j)+dzsdx(im1,j))
                    dzsdy_aa = 0.5*(dzsdy(i,j)+dzsdy(i,jm1))
                    
                    c_x = - ( (1.0-zeta_ac(k))*dzbdx_aa + zeta_ac(k)*dzsdx_aa)
                    c_y = - ( (1.0-zeta_ac(k))*dzbdy_aa + zeta_ac(k)*dzsdy_aa)
                    
                    ! Apply the correction 

                    jvel%dzx(i,j,k) = jvel%dzx(i,j,k) + c_x*jvel%dzz(i,j,k)
                    jvel%dzy(i,j,k) = jvel%dzy(i,j,k) + c_y*jvel%dzz(i,j,k)
                    
                end do 

            end if ! ice present 

        end do 
        end do 
        !$omp end parallel do


        ! Step X: fill in partially filled margin points with neighbor strain-rate values
        
        ! To do....
        
        return 

    end subroutine calc_jacobian_vel_3D

    subroutine calc_strain_rate_tensor_jac(strn, strn2D, jvel, H_ice, f_ice, f_grnd,  &
                                                    zeta_aa, zeta_ac, dx, dy, de_max, boundaries)
        ! -------------------------------------------------------------------------------
        !  Computation of all components of the strain-rate tensor, the full
        !  effective strain rate and the shear fraction.
        ! ------------------------------------------------------------------------------

        ! Note: vx, vy are staggered on ac-nodes in the horizontal, but are on the zeta_aa nodes (ie layer-centered)
        ! in the vertical. vz is centered on aa-nodes in the horizontal, but staggered on zeta_ac nodes
        ! in the vertical. 

        ! Note: first calculate each tensor component on quadrature nodes, then interpolate to aa-nodes)
        ! This is a quadrature approach and is generally more stable. 
        ! The temperorary variable ddn(1:4) is used to hold the values 
        ! calculated at each quadrature point, starting from ddn(1)==upper-right, and
        ! moving counter-clockwise. The average of ddn(1:4) gives the cell-centered
        ! (aa-node) value.

        implicit none
        
        type(strain_3D_class), intent(INOUT) :: strn            ! [yr^-1] on aa-nodes (3D)
        type(strain_2D_class), intent(INOUT) :: strn2D          ! [yr^-1] on aa-nodes (2D)
        type(jacobian_3D_class), intent(IN)  :: jvel            ! 3D velocity Jacobian: Grad([ux,uy,uz])
        real(wp), intent(IN) :: H_ice(:,:)
        real(wp), intent(IN) :: f_ice(:,:)
        real(wp), intent(IN) :: f_grnd(:,:)
        real(wp), intent(IN) :: zeta_aa(:) 
        real(wp), intent(IN) :: zeta_ac(:) 
        real(wp), intent(IN) :: dx
        real(wp), intent(IN) :: dy
        real(wp), intent(IN) :: de_max                          ! [yr^-1] Maximum allowed effective strain rate
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer  :: i, j, k
        integer  :: im1, ip1, jm1, jp1 
        integer  :: nx, ny, nz_aa, nz_ac   
        real(wp) :: lxz, lzx, lyz, lzy
        real(wp) :: shear_squared  
        real(wp), allocatable :: fact_z(:)
        
        real(wp) :: wt0
        real(wp) :: xn(4) 
        real(wp) :: yn(4) 
        real(wp) :: wtn(4)
        real(wp) :: wt1 
        
        real(wp) :: ddn(4) 
        real(wp) :: ddan(4) 
        real(wp) :: ddbn(4) 

        ! Get nodes and weighting 
        wt0 = 1.0/sqrt(3.0)
        xn  = [wt0,-wt0,-wt0, wt0]
        yn  = [wt0, wt0,-wt0,-wt0]
        wtn = [1.0,1.0,1.0,1.0]
        wt1 = sum(wtn)


        ! Determine sizes and allocate local variables 
        nx    = size(H_ice,1)
        ny    = size(H_ice,2)
        nz_aa = size(zeta_aa,1)
        nz_ac = size(zeta_ac,1)
        
        ! Calculate all strain rate tensor components on aa-nodes (horizontally and vertically)
        ! dxx = dxx
        ! dxy = 0.5*(dxy+dyx)
        ! dyy = dyy 
        ! dxz = 0.5*(dxz+dzx)
        ! dyz = 0.5*(dyz+dzy)
        ! dzz = dzz  <= Not calculated, as it is not needed 

        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

            if (f_ice(i,j) .eq. 1.0) then 
                ! Ice is present here, calculate the strain-rate tensor

                ! Loop over all aa-nodes vertically
                do k = 1, nz_aa 

if (.FALSE.) then 
    ! Use quadrature points
                    ! Get dxx on aa-nodes 
                    call acx_to_nodes(ddn,jvel%dxx(:,:,k),i,j,xn,yn,im1,ip1,jm1,jp1)
                    strn%dxx(i,j,k) = sum(ddn*wtn)/wt1

                    ! Get dxy and dyx on aa-nodes 
                    call acx_to_nodes(ddan,jvel%dxy(:,:,k),i,j,xn,yn,im1,ip1,jm1,jp1)
                    call acy_to_nodes(ddbn,jvel%dyx(:,:,k),i,j,xn,yn,im1,ip1,jm1,jp1)
                    ddn = 0.5*(ddan+ddbn)
                    strn%dxy(i,j,k) = sum(ddn*wtn)/wt1

                    ! Get dxz and dzx on aa-nodes 
                    ! (but also get dzx on aa-nodes vertically)
                    call acx_to_nodes(ddan,jvel%dxz(:,:,k),i,j,xn,yn,im1,ip1,jm1,jp1)
                    ddbn = 0.5*(jvel%dzx(i,j,k)+jvel%dzx(i,j,k+1))  ! nz_ac has one more index than nz_aa, so this is ok!
                    ddn  = 0.5*(ddan+ddbn)
                    strn%dxz(i,j,k) = sum(ddn*wtn)/wt1

                    ! Get dyz and dzy on aa-nodes 
                    ! (but also get dzy on aa-nodes vertically)
                    call acy_to_nodes(ddan,jvel%dyz(:,:,k),i,j,xn,yn,im1,ip1,jm1,jp1)
                    ddbn = 0.5*(jvel%dzy(i,j,k)+jvel%dzy(i,j,k+1))  ! nz_ac has one more index than nz_aa, so this is ok!
                    ddn  = 0.5*(ddan+ddbn)
                    strn%dyz(i,j,k) = sum(ddn*wtn)/wt1

                    ! Get dyy on aa-nodes 
                    call acy_to_nodes(ddn,jvel%dyy(:,:,k),i,j,xn,yn,im1,ip1,jm1,jp1)
                    strn%dyy(i,j,k) = sum(ddn*wtn)/wt1
else
    ! Unstagger directly to aa-nodes

                    ! Get dxx on aa-nodes 
                    strn%dxx(i,j,k) = 0.5*(jvel%dxx(i,j,k)+jvel%dxx(im1,j,k))

                    ! Get dxy and dyx on aa-nodes 
                    ddan = 0.5*(jvel%dxy(i,j,k)+jvel%dxy(im1,j,k))
                    ddbn = 0.5*(jvel%dyx(i,j,k)+jvel%dyx(i,jm1,k))
                    ddn = 0.5*(ddan+ddbn)
                    strn%dxy(i,j,k) = sum(ddn*wtn)/wt1

                    ! Get dxz and dzx on aa-nodes 
                    ! (but also get dzx on aa-nodes vertically)
                    ddan = 0.5*(jvel%dxz(i,j,k)+jvel%dxz(im1,j,k))
                    ddbn = 0.5*(jvel%dzx(i,j,k)+jvel%dzx(i,j,k+1))  ! nz_ac has one more index than nz_aa, so this is ok!
                    ddn  = 0.5*(ddan+ddbn)
                    strn%dxz(i,j,k) = sum(ddn*wtn)/wt1

                    ! Get dyz and dzy on aa-nodes 
                    ! (but also get dzy on aa-nodes vertically)
                    ddan = 0.5*(jvel%dyz(i,j,k)+jvel%dyz(i,jm1,k))
                    ddbn = 0.5*(jvel%dzy(i,j,k)+jvel%dzy(i,j,k+1))  ! nz_ac has one more index than nz_aa, so this is ok!
                    ddn  = 0.5*(ddan+ddbn)
                    strn%dyz(i,j,k) = sum(ddn*wtn)/wt1

                    ! Get dyy on aa-nodes 
                    strn%dyy(i,j,k) = 0.5*(jvel%dyy(i,j,k)+jvel%dyy(i,jm1,k))

end if

                    ! TEST - set shear strain terms to zero
                    !strn%dxz(i,j,k) =  0.0 
                    !strn%dyz(i,j,k) = 0.0 

                    ! ====== Finished calculating individual strain rate terms ====== 
                        
                    strn%de(i,j,k) =  sqrt(   strn%dxx(i,j,k)*strn%dxx(i,j,k) &
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
                        shear_squared =   strn%dxz(i,j,k)*strn%dxz(i,j,k) &
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
        strn2D%div     = calc_vertical_integrated_2D(strn%div, zeta_aa)
        strn2D%de      = calc_vertical_integrated_2D(strn%de,  zeta_aa)
        strn2D%f_shear = calc_vertical_integrated_2D(strn%f_shear,zeta_aa) 
        
        ! Finally, calculate the first two eigenvectors for 2D strain rate tensor 
        call calc_2D_eigen_values(strn2D%eps_eig_1,strn2D%eps_eig_2, &
                                    strn2D%dxx,strn2D%dyy,strn2D%dxy)

        return 

    end subroutine calc_strain_rate_tensor_jac 

    subroutine calc_strain_rate_tensor(strn, strn2D, ux, uy, uz, H_ice, f_ice, f_grnd,  &
                                                    zeta_aa, zeta_ac, dx, dy, de_max, boundaries)
        ! -------------------------------------------------------------------------------
        !  Computation of all components of the strain-rate tensor, the full
        !  effective strain rate and the shear fraction.
        ! ------------------------------------------------------------------------------

        ! Note: vx, vy are staggered on ac-nodes in the horizontal, but are on the zeta_aa nodes (ie layer-centered)
        ! in the vertical. vz is centered on aa-nodes in the horizontal, but staggered on zeta_ac nodes
        ! in the vertical. 

        ! Note: first calculate each tensor component on quadrature nodes, then interpolate to aa-nodes)
        ! This is a quadrature approach and is generally more stable. 
        ! The temperorary variable ddn(1:4) is used to hold the values 
        ! calculated at each quadrature point, starting from ddn(1)==upper-right, and
        ! moving counter-clockwise. The average of ddn(1:4) gives the cell-centered
        ! (aa-node) value.

        implicit none
        
        type(strain_3D_class), intent(INOUT) :: strn            ! [yr^-1] on aa-nodes (3D)
        type(strain_2D_class), intent(INOUT) :: strn2D          ! [yr^-1] on aa-nodes (2D)
        real(wp), intent(IN) :: ux(:,:,:)                       ! nx,ny,nz_aa
        real(wp), intent(IN) :: uy(:,:,:)                       ! nx,ny,nz_aa
        real(wp), intent(IN) :: uz(:,:,:)                       ! nx,ny,nz_ac
        real(wp), intent(IN) :: H_ice(:,:)
        real(wp), intent(IN) :: f_ice(:,:)
        real(wp), intent(IN) :: f_grnd(:,:)
        real(wp), intent(IN) :: zeta_aa(:) 
        real(wp), intent(IN) :: zeta_ac(:) 
        real(wp), intent(IN) :: dx
        real(wp), intent(IN) :: dy
        real(wp), intent(IN) :: de_max                          ! [yr^-1] Maximum allowed effective strain rate
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer  :: i, j, k
        integer  :: im1, ip1, jm1, jp1 
        integer  :: nx, ny, nz_aa, nz_ac   
        real(wp) :: H_ice_inv
        real(wp) :: lxz, lzx, lyz, lzy
        real(wp) :: shear_squared  
        real(wp), allocatable :: fact_z(:)
        
        real(wp) :: wt 
        real(wp) :: dd_ab(4) 
        real(wp) :: dd_ab_up(4)
        real(wp) :: dd_ab_dn(4)
        real(wp) :: wt_ab(4) 

        ! Determine sizes and allocate local variables 
        nx    = size(ux,1)
        ny    = size(ux,2)
        nz_aa = size(zeta_aa,1)
        nz_ac = size(zeta_ac,1)
        
        allocate(fact_z(nz_aa))

        !-------- Term abbreviations --------

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

        ! Step 1: Calculate horizontal strain rate components using separate routine
        ! that manages staggering and averaging. 

        do k = 1, nz_aa 
            call calc_strain_rate_tensor_2D(strn2D,ux(:,:,k),uy(:,:,k),H_ice,f_ice,f_grnd,dx,dy,de_max,boundaries)
            strn%dxx(:,:,k) = strn2D%dxx
            strn%dyy(:,:,k) = strn2D%dyy
            strn%dxy(:,:,k) = strn2D%dxy
        end do 

        ! Step 2: Calculate vertical strain rate components 

        !$omp parallel do
        do j=1, ny
        do i=1, nx

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

            if (f_ice(i,j) .eq. 1.0_wp) then 

                H_ice_inv = 1.0_wp/H_ice(i,j)

                ! ====== Loop over each column ====== 

                do k = 1, nz_aa

                    ! ====== Vertical cross terms (lzx,lzy) ====== 

                    ! === lzx ================================

                    if (k .eq. 1) then
                        ! Basal layer 

                        call staggerdiff_nodes_acz_dx_ab_ice(dd_ab,uz(:,:,k),f_ice,i,j,dx)
                    
                    else if (k .eq. nz_aa) then
                        ! Surface layer 

                        call staggerdiff_nodes_acz_dx_ab_ice(dd_ab,uz(:,:,k+1),f_ice,i,j,dx)
                    
                    else 
                        ! Intermediate layers

                        call staggerdiff_nodes_acz_dx_ab_ice(dd_ab_up,uz(:,:,k+1),f_ice,i,j,dx)
                        call staggerdiff_nodes_acz_dx_ab_ice(dd_ab_dn,uz(:,:,k),  f_ice,i,j,dx)
                        
                        dd_ab = 0.5_wp*(dd_ab_dn+dd_ab_up)

                    end if 

                    lzx = sum(wt_ab*dd_ab)

                    ! === lzy ================================

                    if (k .eq. 1) then
                        ! Basal layer

                        call staggerdiff_nodes_acz_dy_ab_ice(dd_ab,uz(:,:,k),f_ice,i,j,dy)

                    else if (k .eq. nz_aa) then 
                        ! Surface layer

                        call staggerdiff_nodes_acz_dy_ab_ice(dd_ab,uz(:,:,k+1),f_ice,i,j,dy)

                    else 
                        ! Intermediate layers

                        call staggerdiff_nodes_acz_dy_ab_ice(dd_ab_up,uz(:,:,k+1),f_ice,i,j,dy)
                        call staggerdiff_nodes_acz_dy_ab_ice(dd_ab_dn,uz(:,:,k),  f_ice,i,j,dy)
                        
                        dd_ab = 0.5_wp*(dd_ab_dn+dd_ab_up)
                    end if 

                    lzy = sum(wt_ab*dd_ab)

                    ! ====== Shear terms (lxz,lyz) ================= 

                    ! === lxz ================================

                    if (k .eq. 1) then 
                        ! Basal layer
                        ! Gradient from first aa-node above base to base 

                        call stagger_nodes_acx_ab_ice(dd_ab_up,ux(:,:,k+1),f_ice,i,j)
                        call stagger_nodes_acx_ab_ice(dd_ab_dn,ux(:,:,k),  f_ice,i,j)

                        dd_ab = (dd_ab_up - dd_ab_dn)*fact_z(k)*H_ice_inv

                    else if (k .eq. nz_aa) then 
                        ! Surface layer
                        ! Gradient from surface to first aa-node below surface 

                        call stagger_nodes_acx_ab_ice(dd_ab_up,ux(:,:,k),  f_ice,i,j)
                        call stagger_nodes_acx_ab_ice(dd_ab_dn,ux(:,:,k-1),f_ice,i,j)

                        dd_ab = (dd_ab_up - dd_ab_dn)*fact_z(k)*H_ice_inv
                        
                    else 
                        ! Intermediate layers
                        ! Gradient from aa-node above to aa-node below

                        call stagger_nodes_acx_ab_ice(dd_ab_up,ux(:,:,k+1),f_ice,i,j)
                        call stagger_nodes_acx_ab_ice(dd_ab_dn,ux(:,:,k-1),f_ice,i,j)

                        dd_ab = (dd_ab_up - dd_ab_dn)*fact_z(k)*H_ice_inv
                        
                    end if 

                    lxz = sum(wt_ab*dd_ab)

                    ! === lyz ================================

                    if (k .eq. 1) then 
                        ! Basal layer
                        ! Gradient from first aa-node above base to base 

                        call stagger_nodes_acy_ab_ice(dd_ab_up,uy(:,:,k+1),f_ice,i,j)
                        call stagger_nodes_acy_ab_ice(dd_ab_dn,uy(:,:,k),  f_ice,i,j)

                        dd_ab = (dd_ab_up - dd_ab_dn)*fact_z(k)*H_ice_inv
                        
                    else if (k .eq. nz_aa) then 
                        ! Surface layer
                        ! Gradient from surface to first aa-node below surface 

                        call stagger_nodes_acy_ab_ice(dd_ab_up,uy(:,:,k),  f_ice,i,j)
                        call stagger_nodes_acy_ab_ice(dd_ab_dn,uy(:,:,k-1),f_ice,i,j)

                        dd_ab = (dd_ab_up - dd_ab_dn)*fact_z(k)*H_ice_inv
                        
                    else 
                        ! Intermediate layers
                        ! Gradient from aa-node above to aa-node below

                        call stagger_nodes_acy_ab_ice(dd_ab_up,uy(:,:,k+1),f_ice,i,j)
                        call stagger_nodes_acy_ab_ice(dd_ab_dn,uy(:,:,k-1),f_ice,i,j)

                        dd_ab = (dd_ab_up - dd_ab_dn)*fact_z(k)*H_ice_inv
                        
                    end if 

                    lyz = sum(wt_ab*dd_ab)

                    ! ====== Calculate cross terms from intermediate values (dxy,dxz,dyz) ====== 

                    strn%dxz(i,j,k) = 0.5_wp*(lxz+lzx)
                    strn%dyz(i,j,k) = 0.5_wp*(lyz+lzy)

                    ! ajr: testing 
                    ! Note: when dxz=dyz=0, the model is more stable. This is clearly
                    ! not correct, but it indicates that these terms are important
                    ! for stability, as expected. Keep in mind for future work. 
                    !strn%dxz(i,j,k) = 0.0_wp
                    !strn%dyz(i,j,k) = 0.0_wp

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
        
        ! Finally, calculate the first two eigenvectors for 2D strain rate tensor 
        call calc_2D_eigen_values(strn2D%eps_eig_1,strn2D%eps_eig_2, &
                                    strn2D%dxx,strn2D%dyy,strn2D%dxy)

        return 

    end subroutine calc_strain_rate_tensor


    subroutine calc_strain_rate_tensor_2D(strn2D,ux,uy,H_ice,f_ice,f_grnd,dx,dy,de_max,boundaries)
        ! Calculate the 2D (vertically averaged) strain rate tensor,
        ! assuming a constant vertical velocity profile. 

        ! ajr: to do: perform calculation on ab-nodes (quadrature)
        ! to be consistent with new formulation above for 
        ! calc_strain_rate_tensor. 

        implicit none

        type(strain_2D_class), intent(INOUT) :: strn2D          ! [yr^-1] Strain rate tensor
        real(wp), intent(IN) :: ux(:,:)                         ! [m/yr] x-velocity, ac-nodes
        real(wp), intent(IN) :: uy(:,:)                         ! [m/yr] y-velocity, ac-nodes
        real(wp), intent(IN) :: H_ice(:,:)                      ! [m] Ice thickness 
        real(wp), intent(IN) :: f_ice(:,:)
        real(wp), intent(IN) :: f_grnd(:,:)
        real(wp), intent(IN) :: dx                              ! [m] Resolution
        real(wp), intent(IN) :: dy                              ! [m] Resolution
        real(wp), intent(IN) :: de_max                          ! [yr^-1] Maximum allowed effective strain rate 
        character(len=*), intent(IN) :: boundaries 

        ! Local variables
        integer  :: i, j, k
        integer  :: im1, ip1, jm1, jp1
        integer  :: nx, ny
        real(wp) :: wt0 
        real(wp) :: xn(4)
        real(wp) :: yn(4)
        real(wp) :: wtn(4) 
        real(wp) :: dudxn(4)
        real(wp) :: dudyn(4)
        real(wp) :: dvdxn(4)
        real(wp) :: dvdyn(4)

        real(wp), allocatable :: dudx(:,:) 
        real(wp), allocatable :: dudy(:,:) 
        real(wp), allocatable :: dvdx(:,:) 
        real(wp), allocatable :: dvdy(:,:) 
        
        nx = size(ux,1)
        ny = size(ux,2)

        allocate(dudx(nx,ny))
        allocate(dudy(nx,ny))
        allocate(dvdx(nx,ny))
        allocate(dvdy(nx,ny))


        ! === First calculate the horizontal strain rate ===


        call calc_strain_rate_horizontal_2D(dudx,dudy,dvdx,dvdy,ux,uy,f_ice,dx,dy,boundaries)


        ! === Next perform interpolations to get strain rate tensor components on aa-nodes ===


        strn2D%dxx      = 0.0 
        strn2D%dyy      = 0.0 
        strn2D%dxy      = 0.0 
        strn2D%dxz      = 0.0       ! Always zero in this case
        strn2D%dyz      = 0.0       ! Always zero in this case
        strn2D%de       = 0.0 
        strn2D%div      = 0.0 
        strn2D%f_shear  = 0.0       ! Always zero in this case

        wt0 = 1.0/sqrt(3.0)
        xn  = [wt0,-wt0,-wt0, wt0]
        yn  = [wt0, wt0,-wt0,-wt0]
        wtn = [1.0,1.0,1.0,1.0]

        !$omp parallel do
        do j=1, ny
        do i=1, nx

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

            if (f_ice(i,j) .eq. 1.0_wp) then 
                
                ! Get strain rate terms on node locations
                call acx_to_nodes(dudxn,dudx,i,j,xn,yn,im1,ip1,jm1,jp1)
                call acx_to_nodes(dudyn,dudy,i,j,xn,yn,im1,ip1,jm1,jp1)
                
                call acy_to_nodes(dvdxn,dvdx,i,j,xn,yn,im1,ip1,jm1,jp1)
                call acy_to_nodes(dvdyn,dvdy,i,j,xn,yn,im1,ip1,jm1,jp1)
                    
                ! Calculate strain rate tensor terms 
                strn2D%dxx(i,j) = sum(dudxn*wtn)/sum(wtn)
                strn2D%dyy(i,j) = sum(dvdyn*wtn)/sum(wtn)
                strn2D%dxy(i,j) = sum(0.5_wp*(dudyn+dvdxn)*wtn)/sum(wtn)

                ! Check tolerance limits
                if (abs(strn2D%dxx(i,j)) .lt. TOL_UNDERFLOW) strn2D%dxx(i,j) = 0.0 
                if (abs(strn2D%dyy(i,j)) .lt. TOL_UNDERFLOW) strn2D%dyy(i,j) = 0.0 
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

        ! Finally, calculate the first two eigenvectors for 2D strain rate tensor 
        call calc_2D_eigen_values(strn2D%eps_eig_1,strn2D%eps_eig_2, &
                                    strn2D%dxx,strn2D%dyy,strn2D%dxy)

        return
        
    end subroutine calc_strain_rate_tensor_2D
    
    subroutine calc_strain_rate_horizontal_2D(dudx,dudy,dvdx,dvdy,ux,uy,f_ice,dx,dy,boundaries)

        implicit none

        real(wp), intent(OUT) :: dudx(:,:) 
        real(wp), intent(OUT) :: dudy(:,:) 
        real(wp), intent(OUT) :: dvdx(:,:) 
        real(wp), intent(OUT) :: dvdy(:,:) 
        real(wp), intent(IN)  :: ux(:,:) 
        real(wp), intent(IN)  :: uy(:,:) 
        real(wp), intent(IN)  :: f_ice(:,:) 
        real(wp), intent(IN)  :: dx 
        real(wp), intent(IN)  :: dy
        character(len=*), intent(IN) :: boundaries         
        
        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: im1, ip1, jm1, jp1 


        nx = size(dudx,1)
        ny = size(dudx,2) 

        ! Populate strain rates over the whole domain on acx- and acy-nodes

        dudx = 0.0
        dvdy = 0.0
        dudy = 0.0
        dvdx = 0.0
        
        do j = 1, ny  
        do i = 1, nx
            
            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

            ! Calculate derivatives (second-order, centered)
            dudx(i,j) = (ux(ip1,j)-ux(im1,j))/(2.0*dx)
            dudy(i,j) = (ux(i,jp1)-ux(i,jm1))/(2.0*dy)
            dvdx(i,j) = (uy(ip1,j)-uy(im1,j))/(2.0*dx)
            dvdy(i,j) = (uy(i,jp1)-uy(i,jm1))/(2.0*dy)

            ! Treat special cases of ice-margin points (take upstream/downstream derivatives instead)

            ! dudx
            if (f_ice(i,j) .eq. 1.0 .and. f_ice(ip1,j) .lt. 1.0) then 
                dudx(i,j) = (ux(i,j)-ux(im1,j))/dx
            else if (f_ice(i,j) .lt. 1.0 .and. f_ice(ip1,j) .eq. 1.0) then 
                dudx(i,j) = (ux(ip1,j)-ux(i,j))/dx
            end if 

            ! dudy
            if (f_ice(i,j) .eq. 1.0 .and. f_ice(i,jp1) .lt. 1.0 .and. f_ice(i,jm1) .lt. 1.0) then 
                dudy(i,j) = 0.0
            else if (f_ice(i,j) .eq. 1.0 .and. f_ice(i,jp1) .lt. 1.0 .and. f_ice(i,jm1) .eq. 1.0) then 
                dudy(i,j) = (ux(i,j)-ux(i,jm1))/dy
            else if (f_ice(i,j) .eq. 1.0 .and. f_ice(i,jp1) .eq. 1.0 .and. f_ice(i,jm1) .lt. 1.0) then 
                dudy(i,j) = (ux(i,jp1)-ux(i,j))/dy
            end if 

            ! dvdy
            if (f_ice(i,j) .eq. 1.0 .and. f_ice(i,jp1) .lt. 1.0) then 
                dvdy(i,j) = (uy(i,j)-uy(i,jm1))/dy
            else if (f_ice(i,j) .lt. 1.0 .and. f_ice(i,jp1) .eq. 1.0) then 
                dvdy(i,j) = (uy(i,jp1)-uy(i,j))/dy
            end if 

            ! dvdx
            if (f_ice(i,j) .eq. 1.0 .and. f_ice(ip1,j) .lt. 1.0 .and. f_ice(im1,j) .lt. 1.0) then 
                dvdx(i,j) = 0.0
            else if (f_ice(i,j) .eq. 1.0 .and. f_ice(ip1,j) .lt. 1.0 .and. f_ice(im1,j) .eq. 1.0) then 
                dvdx(i,j) = (uy(i,j)-uy(im1,j))/dx
            else if (f_ice(i,j) .eq. 1.0 .and. f_ice(ip1,j) .eq. 1.0 .and. f_ice(im1,j) .lt. 1.0) then 
                dvdx(i,j) = (uy(ip1,j)-uy(i,j))/dx
            end if 

        end do
        end do

        ! To do... ? 
        
        ! Further special cases when using 'infinite' boundary conditions 
        ! if (trim(boundaries) .eq. "infinite") then 

        !     dudx(1,:)  = dudx(2,:) 
        !     dudx(nx,:) = dudx(nx-1,:) 
        !     dudx(:,1)  = dudx(:,2) 
        !     dudx(:,ny) = dudx(:,ny-1)
            
        ! end if
                
        return

    end subroutine calc_strain_rate_horizontal_2D

    subroutine calc_strain_rate_uz(dwdx,dwdy,dwdz,uz,H_ice,f_ice,dzbdx,dzbdy, &
                                            dzsdx,dzsdy,zeta_ac,dx,dy,boundaries)

        implicit none
        
        real(wp), intent(OUT) :: dwdx(:,:,:) 
        real(wp), intent(OUT) :: dwdy(:,:,:) 
        real(wp), intent(OUT) :: dwdz(:,:,:) 
        real(wp), intent(IN)  :: uz(:,:,:) 
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: f_ice(:,:) 
        real(wp), intent(IN)  :: dzbdx(:,:) 
        real(wp), intent(IN)  :: dzbdy(:,:)
        real(wp), intent(IN)  :: dzsdx(:,:) 
        real(wp), intent(IN)  :: dzsdy(:,:) 
        real(wp), intent(IN)  :: zeta_ac(:)
        real(wp), intent(IN)  :: dx 
        real(wp), intent(IN)  :: dy
        character(len=*), intent(IN) :: boundaries         
        
        ! Local variables 
        integer  :: i, j, k, nx, ny, nz 
        integer  :: im1, ip1, jm1, jp1 
        real(wp) :: H_inv
        real(wp) :: dzbdx_aa, dzbdy_aa 
        real(wp) :: dzsdx_aa, dzsdy_aa
        real(wp) :: c_x, c_y 

        nx = size(dwdx,1)
        ny = size(dwdx,2) 
        nz = size(dwdx,3)

        ! First calculate dwdz - it is also used as a correction term 
        ! for horizontal derivatives, due to vertical sigma coordinate

        dwdz = 0.0 

        do j = 1, ny  
        do i = 1, nx
            
            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)


            do k = 1, nz 
                ! To do
            end do


        end do 
        end do



        dwdx = 0.0
        dwdy = 0.0
        
        do j = 1, ny  
        do i = 1, nx
            
            if (f_ice(i,j) .eq. 1.0) then 
                ! Ice is present

                ! Get neighbor indices
                call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

                
                ! Calculate standard derivatives (second-order, centered)
                do k = 1, nz
                    dwdx(i,j,k) = (uz(ip1,j,k)-uz(im1,j,k))/(2.0*dx)
                    dwdy(i,j,k) = (uz(i,jp1,k)-uz(i,jm1,k))/(2.0*dy)
                end do          
                 
                ! Treat special cases of ice-margin points (take upstream/downstream derivatives instead)

                ! dwdx

                if (f_ice(ip1,j) .lt. 1.0 .and. f_ice(im1,j) .lt. 1.0) then 
                    do k = 1, nz
                        dwdx(i,j,k) = 0.0
                    end do
                else if (f_ice(ip1,j) .lt. 1.0 .and. f_ice(im1,j) .eq. 1.0) then 
                    do k = 1, nz
                        dwdx(i,j,k) = (uz(i,j,k)-uz(im1,j,k))/dx
                    end do
                else if (f_ice(ip1,j) .eq. 1.0 .and. f_ice(im1,j) .lt. 1.0) then 
                    do k = 1, nz
                        dwdx(i,j,k) = (uz(ip1,j,k)-uz(i,j,k))/dx
                    end do
                end if 

                ! dwdy

                if (f_ice(i,jp1) .lt. 1.0 .and. f_ice(i,jm1) .lt. 1.0) then 
                    do k = 1, nz
                        dwdy(i,j,k) = 0.0
                    end do
                else if (f_ice(i,jp1) .lt. 1.0 .and. f_ice(i,jm1) .eq. 1.0) then 
                    do k = 1, nz
                        dwdy(i,j,k) = (uz(i,j,k)-uz(i,jm1,k))/dy
                    end do
                else if (f_ice(i,jp1) .eq. 1.0 .and. f_ice(i,jm1) .lt. 1.0) then 
                    do k = 1, nz
                        dwdy(i,j,k) = (uz(i,jp1,k)-uz(i,j,k))/dy
                    end do
                end if 

                ! Apply sigma correction factors

                H_inv = 1.0/H_ice(i,j)

                ! Get horizontal grid-centered gradients 
                dzbdx_aa = 1
                dzbdy_aa = 1 
                dzsdx_aa = 1 
                dzsdy_aa = 1 

                do k = 1, nz

                    ! Get correction factors for this grid point on aa-nodes horizontally, ac-nodes vertically
                    c_x = -H_inv * ( (1.0-zeta_ac(k))*dzbdx_aa + zeta_ac(k)*dzsdx_aa )
                    c_y = -H_inv * ( (1.0-zeta_ac(k))*dzbdy_aa + zeta_ac(k)*dzsdy_aa )
                    
                    dwdx(i,j,k) = dwdx(i,j,k) + c_x*dwdz(i,j,k)
                    dwdy(i,j,k) = dwdy(i,j,k) + c_y*dwdz(i,j,k)

                end do

            end if 

        end do 
        end do

        return

    end subroutine calc_strain_rate_uz


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
        call calc_2D_eigen_values(strs2D%tau_eig_1,strs2D%tau_eig_2, &
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
        call calc_2D_eigen_values(strs2D%tau_eig_1,strs2D%tau_eig_2, &
                                    strs2D%txx,strs2D%tyy,strs2D%txy)

        return 

    end subroutine calc_stress_tensor_2D

    elemental subroutine calc_2D_eigen_values(eigen_1,eigen_2,txx,tyy,txy)
        ! Calculate the first two eigenvectors of 2D tensor 
        ! Given A = [txx txy 
        !            tyx tyy], and txy = tyx 
        ! calculate det[A - lambda I]
        ! == det( [txx-lambda txy 
        !          tyx tyy-lambda] )
        !  det[A-lambdaI] = (txx-lambda)* (tyy-lambda) - txy*tyx = 0 
        !    = lambda^2 -lambda txx -lambda tyy + txx * tyy - txy * tyx = 0
        !    = lambda^2 -(txx+tyy)lambda + (txx * tyy - txy * tyx) = 0 
        !
        !    a = 1 
        !    b = -(txx+tyy)
        !    c = (txx * tyy - txy * tyx)
        !
        ! Solve for roots using quadratic formula 
        !

        implicit none

        real(wp), intent(OUT) :: eigen_1              ! [Pa] Eigenvalue 1
        real(wp), intent(OUT) :: eigen_2              ! [Pa] Eigenvalue 2
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
                eigen_1 = lambda1
                eigen_2 = lambda2
            else
                eigen_1 = lambda2
                eigen_2 = lambda1
            end if
        else 
            ! No eigenvalues, set to zero 
            eigen_1 = 0.0_wp 
            eigen_2 = 0.0_wp 
        end if  ! b^2 - 4ac > 0

        return

    end subroutine calc_2D_eigen_values


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

        allocate(strn2D%eps_eig_1(nx,ny))
        allocate(strn2D%eps_eig_2(nx,ny))
        
        strn2D%dxx     = 0.0 
        strn2D%dyy     = 0.0 
        strn2D%dxy     = 0.0
        strn2D%dxz     = 0.0
        strn2D%dyz     = 0.0
        strn2D%div     = 0.0
        strn2D%de      = 0.0 
        strn2D%f_shear = 0.0 

        strn2D%eps_eig_1 = 0.0 
        strn2D%eps_eig_2 = 0.0 
        
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



