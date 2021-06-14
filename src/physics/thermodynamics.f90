module thermodynamics 
    ! This module contains some general thermodynamics subroutines
    ! that could be used by many icetemp solver approaches.
    ! Note: once icetemp is working well, this module could be 
    ! remerged into icetemp as one module. 

    use yelmo_defs, only : prec, wp, sec_year, pi, T0, g, rho_ice, &
                           rho_sw, rho_w, rho_rock, L_ice, T_pmp_beta 

    implicit none 

    private  

    public :: calc_bmb_grounded
    public :: calc_bmb_grounded_enth
    public :: calc_advec_vertical_column
    public :: calc_advec_horizontal_3D
    public :: calc_advec_horizontal_column
    public :: calc_advec_horizontal_column_quick
    public :: calc_strain_heating
    public :: calc_strain_heating_sia
    public :: calc_specific_heat_capacity
    public :: calc_thermal_conductivity
    public :: calc_T_pmp
    public :: calc_f_pmp
    public :: calc_T_base_shlf_approx
    public :: define_temp_linear_3D
    public :: define_temp_robin_3D
    public :: define_temp_linear_column
    public :: define_temp_robin_column
    
    public :: define_temp_bedrock_3D
    public :: define_temp_bedrock_column
    public :: calc_Q_bedrock
    public :: calc_Q_bedrock_column

    public :: calc_basal_water_local

    public :: calc_basal_heating_fromaa 
    public :: calc_basal_heating_fromab 
    public :: calc_basal_heating_fromac 
    public :: calc_basal_heating_interp

    public :: convert_to_enthalpy
    public :: convert_from_enthalpy_column

contains

    subroutine calc_bmb_grounded(bmb_grnd,T_prime_b,Q_ice_b_now,Q_b_now,Q_rock_now,f_grnd,rho_ice)
        ! Calculate everywhere there is at least some grounded ice 
        ! (centered aa node calculation)

        ! Note: calculated bmb_grounded here as if the ice point is fully grounded, 
        ! bmb_grnd and bmb_shlf will then be weighted average using f_grnd externally
        ! (to allow ice topography to evolve with different time steps)

        implicit none 
        
        real(prec), intent(OUT) :: bmb_grnd          ! [m/a ice equiv.] Basal mass balance, grounded
        real(prec), intent(IN)  :: T_prime_b         ! [K] Basal ice temp relative to pressure melting point (ie T_prime_b=0 K == temperate)
        real(prec), intent(IN)  :: Q_ice_b_now       ! [J a-1 m-2] Ice basal heat flux (positive up)
        real(prec), intent(IN)  :: Q_b_now           ! [J a-1 m-2] Basal heat production from friction and strain heating
        real(prec), intent(IN)  :: Q_rock_now        ! [J a-1 m-2] Geothermal heat flux 
        real(prec), intent(IN)  :: f_grnd            ! [--] Grounded fraction (centered aa node)                 
        real(prec), intent(IN)  :: rho_ice           ! [kg m-3] Ice density 
        
        ! Local variables
        real(prec) :: Q_net  
        real(prec), parameter :: tol = 1e-5  

        ! Calculate the grounded basal mass balance following 
        ! Cuffey and Patterson (2010), Eq. 9.38 (Page 420)
        if ( f_grnd .gt. 0.0) then 
            ! Bed is grounded and temperate, calculate basal mass balance  
            ! Classic Cuffey and Patterson (2010) formula
            
            ! Calculate net energy flux at the base [J a-1 m-2]
            Q_net = Q_rock_now + Q_b_now - Q_ice_b_now
            
            bmb_grnd = -Q_net /(rho_ice*L_ice)

            if (T_prime_b .lt. -1.0_prec .and. bmb_grnd .lt. 0.0_prec) then 
                ! Only allow melting for a near-temperate base 
                ! This is a safety-check to prevent strange things from
                ! happening, mainly during initialization when temperatures
                ! are not well defined.
                bmb_grnd = 0.0_prec 
            end if 

        else 
            ! No basal mass change possible if bed is not temperate 

            bmb_grnd = 0.0_prec 

        end if 

        ! Limit small values to avoid underflow errors 
        if (abs(bmb_grnd) .lt. tol) bmb_grnd = 0.0_prec 

        return 

    end subroutine calc_bmb_grounded

    elemental subroutine calc_bmb_grounded_enth(bmb_grnd,T_prime_b,omega,Q_ice_b_now,Q_b_now,Q_rock_now,f_grnd,rho_ice)
        ! Calculate everywhere there is at least some grounded ice 
        ! (centered aa node calculation)

        ! Note: calculated bmb_grounded here as if the ice point is fully grounded, 
        ! bmb_grnd and bmb_shlf will then be weighted average using f_grnd externally
        ! (to allow ice topography to evolve with different time steps)

        implicit none 
        
        real(prec), intent(OUT) :: bmb_grnd          ! [m/a ice equiv.] Basal mass balance, grounded
        real(prec), intent(IN)  :: T_prime_b         ! [K] Basal ice temp relative to pressure melting point (ie T_prime_b=0 K == temperate)
        real(prec), intent(IN)  :: omega 
        real(prec), intent(IN)  :: Q_ice_b_now       ! [J a-1 m-2] Conductive heat flux to the base (positive down)
        real(prec), intent(IN)  :: Q_b_now           ! [J a-1 m-2] Basal heat production from friction and strain heating (postive up)
        real(prec), intent(IN)  :: Q_rock_now        ! [J a-1 m-2] Geothermal heat flux (positive up)
        real(prec), intent(IN)  :: f_grnd            ! [--] Grounded fraction (centered aa node)                 
        real(prec), intent(IN)  :: rho_ice           ! [kg m-3] Ice density 
        
        ! Local variables
        real(prec) :: Q_net  
        real(prec), parameter :: tol = 1e-5  
        
        if (f_grnd .gt. 0.0) then 
            ! Grounded point 

            ! Calculate net energy flux at the base [J a-1 m-2]
            Q_net = Q_rock_now + Q_b_now - Q_ice_b_now
            
            bmb_grnd = -Q_net /(rho_ice*L_ice)
            
        else 
            ! Floating point, no grounded bmb 

            bmb_grnd = 0.0_prec 

        end if 

        ! Limit small values to avoid underflow errors 
        if (abs(bmb_grnd) .lt. tol) bmb_grnd = 0.0_prec 

        return 

    end subroutine calc_bmb_grounded_enth

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

    subroutine calc_advec_horizontal_column_quick(advecxy,var_ice,H_ice,ux,uy,dx,i,j)
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

        real(prec) :: var_w, var_e, var_s, var_n 

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
            if (ux(i,j,k) .gt. 0.0 .and. ux(i-1,j,k) .gt. 0.0 .and. i .ge. 3 .and. i .le. nx-1) then  
                ! Flow to the right 

                ! QUICK scheme 
                var_w =       (6.0/8.0)*var_ice(i-1,j,k)  &
                            + (3.0/8.0)*var_ice(i,j,k)    &
                            - (1.0/8.0)*var_ice(i-2,j,k)

                var_e =       (6.0/8.0)*var_ice(i,j,k)  &
                            + (3.0/8.0)*var_ice(i+1,j,k)    &
                            - (1.0/8.0)*var_ice(i-1,j,k)

                advecx = dx_inv * (var_e*ux(i,j,k) - var_w*ux(i-1,j,k)) 

            else if (ux(i,j,k) .lt. 0.0 .and. ux(i-1,j,k) .lt. 0.0 .and. i .le. nx-2 .and. i .ge. 2) then 
                ! Flow to the left

                ! QUICK scheme 
                var_w =       (6.0/8.0)*var_ice(i,j,k)  &
                            + (3.0/8.0)*var_ice(i-1,j,k)    &
                            - (1.0/8.0)*var_ice(i+1,j,k)

                var_e =       (6.0/8.0)*var_ice(i+1,j,k)  &
                            + (3.0/8.0)*var_ice(i,j,k)    &
                            - (1.0/8.0)*var_ice(i+2,j,k)

                advecx = dx_inv * (var_e*ux(i,j,k) - var_w*ux(i-1,j,k)) 

            else 
                ! No flow 
                advecx = 0.0

            end if 

            if (uy(i,j,k) .gt. 0.0 .and. uy(i,j-1,k) .gt. 0.0 .and. j .ge. 3 .and. j .le. ny-1) then   
                ! Flow to the right 

                ! QUICK scheme 
                var_s =       (6.0/8.0)*var_ice(i,j-1,k)  &
                            + (3.0/8.0)*var_ice(i,j,k)    &
                            - (1.0/8.0)*var_ice(i,j-2,k)

                var_n =       (6.0/8.0)*var_ice(i,j,k)  &
                            + (3.0/8.0)*var_ice(i,j+1,k)    &
                            - (1.0/8.0)*var_ice(i,j-1,k)

                advecy = dx_inv * (var_n*uy(i,j,k) - var_s*uy(i,j-1,k)) 

            else if (uy(i,j,k) .lt. 0.0 .and. uy(i,j-1,k) .lt. 0.0 .and. j .le. ny-2 .and. j .ge. 2) then 
                ! Flow to the left

                ! QUICK scheme 
                var_s =       (6.0/8.0)*var_ice(i,j,k)  &
                            + (3.0/8.0)*var_ice(i,j-1,k)    &
                            - (1.0/8.0)*var_ice(i,j+1,k)

                var_n =       (6.0/8.0)*var_ice(i,j+1,k)  &
                            + (3.0/8.0)*var_ice(i,j,k)    &
                            - (1.0/8.0)*var_ice(i,j+2,k)

                advecy = dx_inv * (var_n*uy(i,j,k) - var_s*uy(i,j-1,k)) 

            else
                ! No flow 
                advecy = 0.0 

            end if 
                    
            ! Combine advection terms for total contribution 
            advecxy(k) = (advecx+advecy)

        end do 

        return 

    end subroutine calc_advec_horizontal_column_quick
    
    subroutine calc_advec_horizontal_column(advecxy,var_ice,H_ice,z_srf,ux,uy,zeta,dx,i,j)
        ! Newly implemented advection algorithms (ajr)
        ! Output: [K a-1]

        ! [m-1] * [m a-1] * [K] = [K a-1]

        implicit none

        real(prec), intent(OUT) :: advecxy(:)       ! nz_aa 
        real(prec), intent(IN)  :: var_ice(:,:,:)   ! nx,ny,nz_aa  Enth, T, age, etc...
        real(prec), intent(IN)  :: H_ice(:,:)       ! nx,ny 
        real(prec), intent(IN)  :: z_srf(:,:)       ! nx,ny 
        real(prec), intent(IN)  :: ux(:,:,:)        ! nx,ny,nz_aa
        real(prec), intent(IN)  :: uy(:,:,:)        ! nx,ny,nz_aa
        real(prec), intent(IN)  :: zeta(:)          ! nz_aa 
        real(prec), intent(IN)  :: dx  
        integer,    intent(IN)  :: i, j 

        ! Local variables 
        integer :: k, nx, ny, nz_aa 
        real(prec) :: ux_aa, uy_aa 
        real(prec) :: dx_inv, dx_inv2
        real(prec) :: advecx, advecy, advec_rev 

        real(prec) :: c_x, c_y, dvardz 

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
            if (ux(i-1,j,k) .gt. 0.0_prec .and. ux(i,j,k) .lt. 0.0_prec .and. i .ge. 3 .and. i .le. nx-2) then 
                ! Convergent flow - take the mean 

                advecx    = dx_inv2 * ux(i-1,j,k)*(-(4.0*var_ice(i-1,j,k)-var_ice(i-2,j,k)-3.0*var_ice(i,j,k)))
                advec_rev = dx_inv2 * ux(i,j,k)*((4.0*var_ice(i+1,j,k)-var_ice(i+2,j,k)-3.0*var_ice(i,j,k)))

                advecx    = 0.5_prec * (advecx + advec_rev) 

            else if (ux_aa .gt. 0.0 .and. i .ge. 3) then  
                ! Flow to the right - inner points

                ! 2nd order
                advecx = dx_inv2 * ux(i-1,j,k)*(-(4.0*var_ice(i-1,j,k)-var_ice(i-2,j,k)-3.0*var_ice(i,j,k)))

            else if (ux_aa .gt. 0.0 .and. i .eq. 2) then  
                ! Flow to the right - border points

                ! 1st order
                advecx = dx_inv * ux(i-1,j,k)*(-(var_ice(i-1,j,k)-var_ice(i,j,k)))
                
            else if (ux_aa .lt. 0.0 .and. i .le. nx-2) then 
                ! Flow to the left

                ! 2nd order
                advecx = dx_inv2 * ux(i,j,k)*((4.0*var_ice(i+1,j,k)-var_ice(i+2,j,k)-3.0*var_ice(i,j,k)))

            else if (ux_aa .lt. 0.0 .and. i .eq. nx-1) then 
                ! Flow to the left

                ! 1st order 
                advecx = dx_inv * ux(i,j,k)*((var_ice(i+1,j,k)-var_ice(i,j,k)))
                
            else 
                ! No flow or divergent 

                advecx = 0.0

            end if 

            if (uy(i,j-1,k) .gt. 0.0_prec .and. uy(i,j,k) .lt. 0.0_prec .and. j .ge. 3 .and. j .le. ny-2) then 
                ! Convergent flow - take the mean 

                advecy    = dx_inv2 * uy(i,j-1,k)*(-(4.0*var_ice(i,j-1,k)-var_ice(i,j-2,k)-3.0*var_ice(i,j,k)))
                advec_rev = dx_inv2 * uy(i,j,k)*((4.0*var_ice(i,j+1,k)-var_ice(i,j+2,k)-3.0*var_ice(i,j,k)))
                
                advecy    = 0.5_prec * (advecy + advec_rev) 

            else if (uy_aa .gt. 0.0 .and. j .ge. 3) then   
                ! Flow to the right  - inner points

                ! 2nd order
                advecy = dx_inv2 * uy(i,j-1,k)*(-(4.0*var_ice(i,j-1,k)-var_ice(i,j-2,k)-3.0*var_ice(i,j,k)))

            else if (uy_aa .gt. 0.0 .and. j .eq. 2) then   
                ! Flow to the right - border points

                ! 1st order
                advecy = dx_inv * uy(i,j-1,k)*(-(var_ice(i,j-1,k)-var_ice(i,j,k)))
                
            else if (uy_aa .lt. 0.0 .and. j .le. ny-2) then 
                ! Flow to the left

                ! 2nd order
                advecy = dx_inv2 * uy(i,j,k)*((4.0*var_ice(i,j+1,k)-var_ice(i,j+2,k)-3.0*var_ice(i,j,k)))
            
            else if (uy_aa .lt. 0.0 .and. j .eq. ny-1) then 
                ! Flow to the left

                ! 1st order
                advecy = dx_inv * uy(i,j,k)*((var_ice(i,j+1,k)-var_ice(i,j,k)))
                  
            else
                ! No flow 
                advecy = 0.0 

            end if 
            
            ! Combine advection terms for total contribution 
            advecxy(k) = (advecx+advecy)

        end do 

        return 

    end subroutine calc_advec_horizontal_column
    
    subroutine calc_advec_horizontal_3D(advecxy,var,H_ice,z_srf,ux,uy,zeta_aa,dx,beta1,beta2)

        implicit none 

        real(prec), intent(INOUT) :: advecxy(:,:,:)     ! nz_aa 
        real(prec), intent(IN)    :: var(:,:,:)         ! nx,ny,nz_aa  Enth, T, age, etc...
        real(prec), intent(IN)    :: H_ice(:,:)         ! nx,ny 
        real(prec), intent(IN)    :: z_srf(:,:)         ! nx,ny 
        real(prec), intent(IN)    :: ux(:,:,:)          ! nx,ny,nz_aa
        real(prec), intent(IN)    :: uy(:,:,:)          ! nx,ny,nz_aa
        real(prec), intent(IN)    :: zeta_aa(:)         ! nz_aa 
        real(prec), intent(IN)    :: dx  
        real(prec), intent(IN)    :: beta1              ! Weighting term for multistep advection scheme
        real(prec), intent(IN)    :: beta2              ! Weighting term for multistep advection scheme

        ! Local variables 
        integer :: i, j 
        integer :: nx, ny, nz  

        real(prec), allocatable   :: advecxy_nm1(:,:,:) ! nz_aa, advective term from previous timestep

        nx = size(advecxy,1)
        ny = size(advecxy,2) 
        nz = size(advecxy,3) 

        allocate(advecxy_nm1(nx,ny,nz))

        ! Store previous solution for later use 
        advecxy_nm1 = advecxy 

        ! Reset current solution to zero 
        advecxy = 0.0_prec 
         
        do j = 2, ny-1
        do i = 2, nx-1 
            call calc_advec_horizontal_column(advecxy(i,j,:),var,H_ice,z_srf,ux,uy,zeta_aa,dx,i,j)
        end do 
        end do 

        ! Fill in boundaries 
        j = 1 
        do i = 2, nx-1 
            if(H_ice(i,j) .gt. 0.0_prec) advecxy(i,j,:) = advecxy(i,j+1,:) 
        end do 

        j = ny 
        do i = 2, nx-1 
            if(H_ice(i,j) .gt. 0.0_prec) advecxy(i,j,:) = advecxy(i,j-1,:) 
        end do 
        
        i = 1 
        do j = 2, ny-1 
            if(H_ice(i,j) .gt. 0.0_prec) advecxy(i,j,:) = advecxy(i+1,j,:) 
        end do 

        i = nx 
        do j = 2, ny-1 
            if(H_ice(i,j) .gt. 0.0_prec) advecxy(i,j,:) = advecxy(i-1,j,:) 
        end do 

        i = 1
        j = 1 
        if(H_ice(i,j) .gt. 0.0_prec) advecxy(i,j,:) = 0.5_prec*(advecxy(i+1,j,:)+advecxy(i,j+1,:))

        i = nx
        j = 1 
        if(H_ice(i,j) .gt. 0.0_prec) advecxy(i,j,:) = 0.5_prec*(advecxy(i-1,j,:)+advecxy(i,j+1,:))

        i = nx
        j = ny
        if(H_ice(i,j) .gt. 0.0_prec) advecxy(i,j,:) = 0.5_prec*(advecxy(i-1,j,:)+advecxy(i,j-1,:))

        i = 1
        j = ny
        if(H_ice(i,j) .gt. 0.0_prec) advecxy(i,j,:) = 0.5_prec*(advecxy(i+1,j,:)+advecxy(i,j-1,:))
        

        ! Calculate weighted average between current and previous solution following 
        ! timestepping method desired 
        advecxy = beta1*advecxy + beta2*advecxy_nm1 
         
        return 

    end subroutine calc_advec_horizontal_3D

    subroutine calc_strain_heating(Q_strn,de,visc,cp,rho_ice,beta1,beta2)
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

        real(prec), intent(INOUT) :: Q_strn(:,:,:)      ! nx,ny,nz_aa [K a-1] Heat production
        real(prec), intent(IN)    :: de(:,:,:)          ! nx,ny,nz_aa [a-1] Effective strain rate 
        real(prec), intent(IN)    :: visc(:,:,:)        ! nx,ny,nz_aa [Pa a-1] Viscosity
        real(prec), intent(IN)    :: cp(:,:,:)          ! nx,ny,nz_aa [J kg-1 K-1] Specific heat capacity
        real(prec), intent(IN)    :: rho_ice            ! [kg m-3] Ice density 
        real(prec), intent(IN)    :: beta1              ! Timestepping weighting parameter
        real(prec), intent(IN)    :: beta2              ! Timestepping weighting parameter
        
        ! Local variables
        integer :: nz_aa 

        nz_aa = size(Q_strn,3)

        ! Calculate strain heating 
        ! following Greve and Blatter (2009), Eqs. 4.7 and 5.65
        Q_strn = beta1*(4.0*visc * de**2) + beta2*Q_strn 
        
        ! Ensure Q_strn is strictly positive 
        where (Q_strn .lt. 0.0_prec) Q_strn = 0.0_prec 

        return 

    end subroutine calc_strain_heating
    
    subroutine calc_strain_heating_sia(Q_strn,ux,uy,dzsdx,dzsdy,cp,H_ice,rho_ice,zeta_aa,zeta_ac,beta1,beta2)

        ! Calculate the general 3D internal strain heating
        ! as sum(D_ij*tau_ij)  (strain*stress)
        ! where stress has been calculated as stress_ij = 2*visc*strain
        ! Units: Q_strn = Q [J a-1 m-3]

        ! SIA approximation:
        ! Q_strn = rho*g*H*(duxdz*dzsdx + duydz*dzsdy)
        ! Units: [J a-1 m-3]

        implicit none

        real(prec), intent(INOUT) :: Q_strn(:,:,:)    ! nx,ny,nz_aa  [Pa m a-1 ??] Heat production
        real(prec), intent(IN)    :: ux(:,:,:)        ! nx,ny,nz_aa  [m a-1] Velocity x-direction
        real(prec), intent(IN)    :: uy(:,:,:)        ! nx,ny,nz_aa  [m a-1] Velocity y-direction
        real(prec), intent(IN)    :: dzsdx(:,:)       ! nx,ny        [m m-1] Surface slope x-direction
        real(prec), intent(IN)    :: dzsdy(:,:)       ! nx,ny        [m m-1] Surface slope y-direction
        real(prec), intent(IN)    :: cp(:,:,:)        ! nx,ny,nz_aa  [J/kg/K] Specific heat capacity
        real(prec), intent(IN)    :: H_ice(:,:)       ! nx,ny        [m] Ice thickness
        real(prec), intent(IN)    :: rho_ice          ! [kg m-3] Ice density 
        real(prec), intent(IN)    :: zeta_aa(:)       ! [-] Height axis, centered aa-nodes 
        real(prec), intent(IN)    :: zeta_ac(:)       ! [-] Height axis, boundaries ac-nodes
        real(prec), intent(IN)    :: beta1            ! Timestepping weighting parameter
        real(prec), intent(IN)    :: beta2            ! Timestepping weighting parameter
        
        ! Local variables
        integer :: i, j, k, nx, ny, nz_aa 
        real(prec) :: ux_aa_up, ux_aa_dwn
        real(prec) :: uy_aa_up, uy_aa_dwn
        real(prec) :: duxdz, duydz
        real(prec) :: dzsdx_aa, dzsdy_aa  
        real(prec) :: dz, depth 
        real(prec) :: Q_strn_now 

        nx    = size(Q_strn,1)
        ny    = size(Q_strn,2)
        nz_aa = size(Q_strn,3)

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

                    dz = H_ice(i,j)*(zeta_ac(k+1) - zeta_ac(k)) 
                    
                    duxdz = (ux_aa_up-ux_aa_dwn)/dz 
                    duydz = (uy_aa_up-uy_aa_dwn)/dz 

                    dzsdx_aa = 0.5*(dzsdx(i-1,j)+dzsdx(i,j))
                    dzsdy_aa = 0.5*(dzsdy(i,j-1)+dzsdy(i,j))
                    
                    depth = H_ice(i,j)*(1.0-zeta_aa(k))
                    
                    Q_strn_now = (-rho_ice*g*depth) * (duxdz*dzsdx_aa + duydz*dzsdy_aa)
                    Q_strn(i,j,k) = beta1*Q_strn_now + beta2*Q_strn(i,j,k) 
                    
                    ! Ensure Q_strn is strictly positive 
                    if (Q_strn(i,j,k) .lt. 0.0_prec) Q_strn(i,j,k) = 0.0_prec 
                
                end do 

            else 
                ! No Q_strn outside of ice sheet 
                
                Q_strn(i,j,k) = 0.0_prec 

            end if 

        end do 
        end do 

        return 

    end subroutine calc_strain_heating_sia

    subroutine calc_basal_heating_fromaa(Q_b,ux_b,uy_b,taub_acx,taub_acy,beta1,beta2)
         ! Qb [J a-1 m-2] == [m a-1] * [J m-3]
         ! Note: grounded ice fraction f_grnd_acx/y not used here, because taub_acx/y already accounts
         ! for the grounded fraction via beta_acx/y: Q_b = tau_b*u = -beta*u*u.

        real(prec), intent(INOUT) :: Q_b(:,:)           ! [mW m-2] Basal heat production (friction), aa-nodes
        real(prec), intent(IN)    :: ux_b(:,:)          ! Basal velocity, x-component (acx)
        real(prec), intent(IN)    :: uy_b(:,:)          ! Basal velocity, y-compenent (acy)
        real(prec), intent(IN)    :: taub_acx(:,:)      ! Basal friction (acx)
        real(prec), intent(IN)    :: taub_acy(:,:)      ! Basal friction (acy) 
        real(prec), intent(IN)    :: beta1              ! Timestepping weighting parameter
        real(prec), intent(IN)    :: beta2              ! Timestepping weighting parameter
        
        ! Local variables
        integer    :: i, j, nx, ny, n 
        integer    :: im1, ip1, jm1, jp1 
        real(prec) :: uxy_aa, taub_aa 
        real(prec) :: Q_b_now

        nx = size(Q_b,1)
        ny = size(Q_b,2)

        ! First calculate basal frictional heating on ab-nodes 
        do j = 1, ny
        do i = 1, nx
            
            im1 = max(1,i-1)
            ip1 = min(nx,i+1)
            jm1 = max(1,j-1)
            jp1 = min(ny,j+1)
            
            uxy_aa  = sqrt( (0.5_prec*(ux_b(i,j)+ux_b(im1,j)))**2 &
                          + (0.5_prec*(uy_b(i,j)+uy_b(i,jm1)))**2 )

            taub_aa = sqrt( (0.5_prec*(taub_acx(i,j)+taub_acx(im1,j)))**2 &
                          + (0.5_prec*(taub_acy(i,j)+taub_acy(i,jm1)))**2 )
            
            Q_b_now = abs(uxy_aa*taub_aa)      ! [Pa m a-1] == [J a-1 m-2]

            ! Get weighted average of Q_b with timestepping factors
            Q_b(i,j) = beta1*Q_b_now + beta2*Q_b(i,j) 

            ! Ensure Q_b is strictly positive 
            if (Q_b(i,j) .lt. 0.0_prec) Q_b(i,j) = 0.0_prec 
            
        end do 
        end do 

        ! Finally convert to [mW m-2]
        Q_b = Q_b * 1e3 / sec_year      ! [J a-1 m-2] => [mW m-2]

        return 
 
    end subroutine calc_basal_heating_fromaa

    subroutine calc_basal_heating_fromab(Q_b,ux_b,uy_b,taub_acx,taub_acy,f_ice,beta1,beta2)
         ! Qb [J a-1 m-2] == [m a-1] * [J m-3]
         ! Note: grounded ice fraction f_grnd_acx/y not used here, because taub_acx/y already accounts
         ! for the grounded fraction via beta_acx/y: Q_b = tau_b*u = -beta*u*u.

        real(prec), intent(INOUT) :: Q_b(:,:)           ! [mW m-2] Basal heat production (friction), aa-nodes
        real(prec), intent(IN)    :: ux_b(:,:)          ! Basal velocity, x-component (acx)
        real(prec), intent(IN)    :: uy_b(:,:)          ! Basal velocity, y-compenent (acy)
        real(prec), intent(IN)    :: taub_acx(:,:)      ! Basal friction (acx)
        real(prec), intent(IN)    :: taub_acy(:,:)      ! Basal friction (acy) 
        real(prec), intent(IN)    :: f_ice(:,:)         ! [--] Ice area fraction
        real(prec), intent(IN)    :: beta1              ! Timestepping weighting parameter
        real(prec), intent(IN)    :: beta2              ! Timestepping weighting parameter
        
        ! Local variables
        integer    :: i, j, nx, ny, n 
        integer    :: im1, ip1, jm1, jp1 
        real(wp) :: uxb_ab(4) 
        real(wp) :: uyb_ab(4) 
        real(wp) :: taubx_ab(4) 
        real(wp) :: tauby_ab(4)
        real(wp) :: Qb_ab(4) 
        real(wp) :: wt_ab(4) 
        real(wp) :: Q_b_now
        real(wp) :: wt

        nx = size(Q_b,1)
        ny = size(Q_b,2)

        ! Calculate current basal frictional heating at ab-nodes
        ! following quadrature method, and then average to get 
        ! the grid center (using contributions from ice-covered nodes)
        do j = 1, ny
        do i = 1, nx
            
            im1 = max(1,i-1)
            ip1 = min(nx,i+1)
            jm1 = max(1,j-1)
            jp1 = min(ny,j+1)
            
            ! Get ab-node weighting based on whether ice is present 
            wt_ab = 0.0_wp 
            if (count([f_ice(i,j),f_ice(ip1,j),f_ice(i,jp1),f_ice(ip1,jp1)].lt.1.0_wp) .eq. 0) then 
                wt_ab(1) = 1.0_wp
            end if
            if (count([f_ice(i,j),f_ice(im1,j),f_ice(im1,jp1),f_ice(i,jp1)].lt.1.0_wp) .eq. 0) then 
                wt_ab(2) = 1.0_wp 
            end if 
            if (count([f_ice(i,j),f_ice(i,jm1),f_ice(ip1,jm1),f_ice(ip1,j)].lt.1.0_wp) .eq. 0) then 
                wt_ab(3) = 1.0_wp 
            end if 
            
            if (count([f_ice(i,j),f_ice(im1,j),f_ice(im1,jm1),f_ice(i,jm1)].lt.1.0_wp) .eq. 0) then 
                wt_ab(4) = 1.0_wp
            end if 
            
            wt = sum(wt_ab)
            
            if (f_ice(i,j) .eq. 1.0_wp .and. wt .gt. 0.0_wp) then 
                ! Fully ice-covered point with some fully ice-covered neighbors 

                ! Normalize weighting 
                wt_ab = wt_ab / wt 

                ! === ux =========

                uxb_ab(1) = 0.5_wp*(ux_b(i,j)+ux_b(i,jp1))
                uxb_ab(2) = 0.5_wp*(ux_b(im1,j)+ux_b(im1,jp1))
                uxb_ab(3) = 0.5_wp*(ux_b(im1,jm1)+ux_b(im1,j))
                uxb_ab(4) = 0.5_wp*(ux_b(i,jm1)+ux_b(i,j))

                ! === uy =========

                uyb_ab(1) = 0.5_wp*(uy_b(i,j)+uy_b(ip1,j))
                uyb_ab(2) = 0.5_wp*(uy_b(i,j)+uy_b(im1,j))
                uyb_ab(3) = 0.5_wp*(uy_b(i,jm1)+uy_b(im1,jm1))
                uyb_ab(4) = 0.5_wp*(uy_b(i,jm1)+uy_b(ip1,jm1))

                ! === taubx =========

                taubx_ab(1) = 0.5_wp*(taub_acx(i,j)+taub_acx(i,jp1))
                taubx_ab(2) = 0.5_wp*(taub_acx(im1,j)+taub_acx(im1,jp1))
                taubx_ab(3) = 0.5_wp*(taub_acx(im1,jm1)+taub_acx(im1,j))
                taubx_ab(4) = 0.5_wp*(taub_acx(i,jm1)+taub_acx(i,j))

                ! === tauby =========

                tauby_ab(1) = 0.5_wp*(taub_acy(i,j)+taub_acy(ip1,j))
                tauby_ab(2) = 0.5_wp*(taub_acy(i,j)+taub_acy(im1,j))
                tauby_ab(3) = 0.5_wp*(taub_acy(i,jm1)+taub_acy(im1,jm1))
                tauby_ab(4) = 0.5_wp*(taub_acy(i,jm1)+taub_acy(ip1,jm1))

                ! Calculate Qb at quadrature points [Pa m a-1] == [J a-1 m-2]
                Qb_ab = abs( sqrt(uxb_ab**2+uyb_ab**2) &
                              * sqrt(taubx_ab**2+tauby_ab**2) )

                ! Get weighted average using only ice-covered nodes
                Q_b_now = sum(wt_ab*Qb_ab)     

            else 

                Q_b_now = 0.0_wp 

            end if 

            ! Get weighted average of Q_b with timestepping factors
            Q_b(i,j) = beta1*Q_b_now + beta2*Q_b(i,j) 

            ! Ensure Q_b is strictly positive 
            if (Q_b(i,j) .lt. 0.0_prec) Q_b(i,j) = 0.0_prec 
            
        end do 
        end do 

        ! Finally convert to [mW m-2]
        Q_b = Q_b * 1e3 / sec_year      ! [J a-1 m-2] => [mW m-2]

        return 
 
    end subroutine calc_basal_heating_fromab

    subroutine calc_basal_heating_fromac(Q_b,ux_b,uy_b,taub_acx,taub_acy,H_ice,f_pmp,beta1,beta2)
         ! Qb [J a-1 m-2] == [m a-1] * [J m-3]
         ! Note: grounded ice fraction f_grnd_acx/y not used here, because taub_acx/y already accounts
         ! for the grounded fraction via beta_acx/y: Q_b = tau_b*u = -beta*u*u.

        real(prec), intent(INOUT) :: Q_b(:,:)           ! [mW m-2] Basal heat production (friction), aa-nodes
        real(prec), intent(IN)    :: ux_b(:,:)          ! Basal velocity, x-component (acx)
        real(prec), intent(IN)    :: uy_b(:,:)          ! Basal velocity, y-compenent (acy)
        real(prec), intent(IN)    :: taub_acx(:,:)      ! Basal friction (acx)
        real(prec), intent(IN)    :: taub_acy(:,:)      ! Basal friction (acy) 
        real(prec), intent(IN)    :: f_pmp(:,:)         ! [--] Grid-fraction at pmp (aa-nodes)
        real(prec), intent(IN)    :: H_ice(:,:)         ! [m] Ice thickness 
        real(prec), intent(IN)    :: beta1              ! Timestepping weighting parameter
        real(prec), intent(IN)    :: beta2              ! Timestepping weighting parameter
        
        ! Local variables
        integer    :: i, j, nx, ny, n 
        integer    :: im1, jm1, ip1, jp1 
        real(prec), allocatable :: Qb_acx(:,:)
        real(prec), allocatable :: Qb_acy(:,:)
        real(prec) :: Q_b_now

        nx = size(Q_b,1)
        ny = size(Q_b,2)

        allocate(Qb_acx(nx,ny))
        allocate(Qb_acy(nx,ny))

        ! Determine basal frictional heating values (staggered acx/acy nodes)
        Qb_acx = abs(ux_b*taub_acx)   ! [Pa m a-1] == [J a-1 m-2]
        Qb_acy = abs(uy_b*taub_acy)   ! [Pa m a-1] == [J a-1 m-2]


        ! Get basal frictional heating on centered nodes (aa-nodes)          
        do j = 1, ny
        do i = 1, nx

            im1 = max(1,i-1)
            ip1 = min(nx,i+1)
            jm1 = max(1,j-1)
            jp1 = min(ny,j+1)
            
            ! Average from ac-nodes to aa-node
            Q_b_now = 0.25*(Qb_acx(i,j)+Qb_acx(im1,j)+Qb_acy(i,j)+Qb_acy(i,jm1))

if (.FALSE.) then 
            ! Reduction of Q_b with f_pmp (ajr: this is an open question)
            Q_b_now = Q_b_now*f_pmp(i,j)  
end if 
            
            ! Get weighted average of Q_b with timestepping factors
            Q_b(i,j) = beta1*Q_b_now + beta2*Q_b(i,j) 

            ! Ensure Q_b is strictly positive 
            if (Q_b(i,j) .lt. 0.0_prec) Q_b(i,j) = 0.0_prec 
            
        end do 
        end do 
        
        ! Finally convert to [mW m-2]
        Q_b = Q_b * 1e3 / sec_year      ! [J a-1 m-2] => [mW m-2]

        return 
 
    end subroutine calc_basal_heating_fromac

    subroutine calc_basal_heating_interp(Q_b,ux_b,uy_b,taub_acx,taub_acy,H_ice,beta1,beta2)
         ! Qb [J a-1 m-2] == [m a-1] * [J m-3]
         ! Note: grounded ice fraction f_grnd_acx/y not used here, because beta_acx/y already accounts
         ! for the grounded fraction via beta_acx/y: Q_b = tau_b*u = -beta*u*u.

        real(prec), intent(INOUT) :: Q_b(:,:)           ! [mW m-2] Basal heat production (friction), aa-nodes
        real(prec), intent(IN)    :: ux_b(:,:)          ! Basal velocity, x-component (acx)
        real(prec), intent(IN)    :: uy_b(:,:)          ! Basal velocity, y-compenent (acy)
        real(prec), intent(IN)    :: taub_acx(:,:)      ! Basal stress (acx)
        real(prec), intent(IN)    :: taub_acy(:,:)      ! Basal stress (acy) 
        real(prec), intent(IN)    :: H_ice(:,:)         ! [m] Ice thickness
        real(prec), intent(IN)    :: beta1              ! Timestepping weighting parameter
        real(prec), intent(IN)    :: beta2              ! Timestepping weighting parameter
        
        ! Local variables
        integer    :: i, j, nx, ny, n 
        integer    :: im1, ip1, jm1, jp1 
        real(prec) :: ux_1, ux_2, ux_3, ux_4 
        real(prec) :: uy_1, uy_2, uy_3, uy_4 
        real(prec) :: taubx_1, taubx_2, taubx_3, taubx_4 
        real(prec) :: tauby_1, tauby_2, tauby_3, tauby_4
        
        real(prec) :: uxy_ab, taub_ab 
        real(prec), allocatable :: Qb_ab(:,:)
        real(prec) :: Q_b_now, wt

        integer, parameter :: nhi = 10 
        real(prec), allocatable :: uxhi(:,:)
        real(prec), allocatable :: uyhi(:,:)
        real(prec), allocatable :: taubxhi(:,:)
        real(prec), allocatable :: taubyhi(:,:)
        real(prec), allocatable :: Qbhi(:,:)
        real(prec) :: Qbhi_tot 

        integer, parameter :: ix = 36
        integer, parameter :: jx =  9
        
        nx = size(Q_b,1)
        ny = size(Q_b,2)

        allocate(Qb_ab(nx,ny))

        allocate(uxhi(nhi,nhi))
        allocate(uyhi(nhi,nhi))
        allocate(taubxhi(nhi,nhi))
        allocate(taubyhi(nhi,nhi))
        allocate(Qbhi(nhi,nhi))
        
        ! First calculate basal frictional heating on ab-nodes 
        do j = 1, ny
        do i = 1, nx
            
            if (H_ice(i,j) .gt. 0.0_prec) then 
                ! only treat ice-covered points 

                im1 = max(1,i-1)
                ip1 = min(nx,i+1)
                jm1 = max(1,j-1)
                jp1 = min(ny,j+1)
                
                ! Get variables on ab-nodes (corners of cell)
                ! and calculate hi resolution subgrid values
                ux_1 = 0.5_prec*(ux_b(i,j)  +ux_b(i,jp1))
                ux_2 = 0.5_prec*(ux_b(im1,j)+ux_b(im1,jp1))
                ux_3 = 0.5_prec*(ux_b(im1,j)+ux_b(im1,jm1))
                ux_4 = 0.5_prec*(ux_b(i,j)  +ux_b(i,jm1))

                call calc_hires_cell(uxhi,ux_1,ux_2,ux_3,ux_4)

                uy_1 = 0.5_prec*(uy_b(i,j)  +uy_b(ip1,j))
                uy_2 = 0.5_prec*(uy_b(i,j)  +uy_b(im1,j))
                uy_3 = 0.5_prec*(uy_b(i,jm1)+uy_b(im1,jm1))
                uy_4 = 0.5_prec*(uy_b(i,jm1)+uy_b(ip1,jm1))

                call calc_hires_cell(uyhi,uy_1,uy_2,uy_3,uy_4)

                taubx_1 = 0.5_prec*(taub_acx(i,j)  +taub_acx(i,jp1))
                taubx_2 = 0.5_prec*(taub_acx(im1,j)+taub_acx(im1,jp1))
                taubx_3 = 0.5_prec*(taub_acx(im1,j)+taub_acx(im1,jm1))
                taubx_4 = 0.5_prec*(taub_acx(i,j)  +taub_acx(i,jm1))

                call calc_hires_cell(taubxhi,taubx_1,taubx_2,taubx_3,taubx_4)

                tauby_1 = 0.5_prec*(taub_acy(i,j)  +taub_acy(ip1,j))
                tauby_2 = 0.5_prec*(taub_acy(i,j)  +taub_acy(im1,j))
                tauby_3 = 0.5_prec*(taub_acy(i,jm1)+taub_acy(im1,jm1))
                tauby_4 = 0.5_prec*(taub_acy(i,jm1)+taub_acy(ip1,jm1))

                call calc_hires_cell(taubyhi,tauby_1,tauby_2,tauby_3,tauby_4) 

                Qbhi = abs( sqrt(uxhi**2 + uyhi**2) * sqrt(taubxhi**2 + taubyhi**2) )  ! [Pa m a-1] == [J a-1 m-2]

                Qbhi_tot = sum(Qbhi)

                if (Qbhi_tot .gt. 0.0_prec) then 
                    
                    Q_b_now = sum(Qbhi,mask=Qbhi .gt. 0.0_prec) &
                                            / count(Qbhi .gt. 0.0_prec)

                else 

                    Q_b_now = 0.0_prec 

                end if 

            else 

                Q_b_now = 0.0_prec 

            end if 

            ! Get weighted average of Q_b with timestepping factors
            Q_b(i,j) = beta1*Q_b_now + beta2*Q_b(i,j) 

            ! Ensure Q_b is strictly positive 
            if (Q_b(i,j) .lt. 0.0_prec) Q_b(i,j) = 0.0_prec 
            
        end do 
        end do 

        ! Finally convert to [mW m-2]
        Q_b = Q_b * 1e3 / sec_year      ! [J a-1 m-2] => [mW m-2]

        return 
 
    end subroutine calc_basal_heating_interp

    subroutine calc_hires_cell(var_hi,var_1,var_2,var_3,var_4)
        ! Given the four corners of a cell in quadrants 1,2,3,4,
        ! calculate a high resolution grid of the variable
        ! through bilinear interpolation that covers the cell.

        implicit none 

        real(prec), intent(OUT) :: var_hi(:,:) 
        real(prec), intent(IN)  :: var_1,var_2,var_3,var_4
        
        ! Local variables 
        integer :: i, j, nx  
        real(prec) :: dx 
        real(prec) :: x(size(var_hi,1)), y(size(var_hi,2)) 

        nx = size(var_hi,1)

        dx = 1/real(nx,prec)

        ! Populate x,y axes for interpolation points (between 0 and 1)
        ! Note: x and y points are offset from values of 0 and 1 to be
        ! sure that each mini grid box has the same area contribution
        ! to the total grid cell.
        do i = 1, nx 
            !x(i) = 0.0_prec + real(i-1)/real(nx-1)
            x(i) = 0.5_prec*dx + real(i-1)*dx 
        end do 
        y = x 


        ! Calculate linear interpolation value 
        var_hi = 0.0_prec 

        do i = 1, nx 
        do j = 1, nx 

            var_hi(i,j) = interp_bilin_pt(var_1,var_2,var_3,var_4,x(i),y(j))

        end do 
        end do 

        return 

    end subroutine calc_hires_cell

    function interp_bilin_pt(z1,z2,z3,z4,xout,yout) result(zout)
        ! Interpolate a point given four neighbors at corners of square (0:1,0:1)
        ! z2    z1
        !    x,y
        ! z3    z4 
        ! 

        implicit none 

        real(prec), intent(IN) :: z1, z2, z3, z4 
        real(prec), intent(IN) :: xout, yout 
        real(prec) :: zout 

        ! Local variables 
        real(prec) :: x0, x1, y0, y1 
        real(prec) :: alpha1, alpha2, p0, p1 

        x0 = 0.0_prec 
        x1 = 1.0_prec
        y0 = 0.0_prec
        y1 = 1.0_prec

        alpha1  = (xout - x0) / (x1-x0)
        p0      = z3 + alpha1*(z4-z3)
        p1      = z2 + alpha1*(z1-z2)
            
        alpha2  = (yout - y0) / (y1-y0)
        zout    = p0 + alpha2*(p1-p0)

        return 

    end function interp_bilin_pt

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
    
    elemental function calc_T_pmp(H_ice,zeta,T0,beta) result(T_pmp)
        ! Greve and Blatter (Chpt 4, pg 54), Eq. 4.13
        ! This gives the pressure-corrected melting point of ice
        ! where H_ice*(1-zeta) is the thickness of ice overlying the current point 
        
        implicit none 

        real(prec), intent(IN) :: H_ice  ! [m] Total ice thickness of this point
        real(prec), intent(IN) :: zeta   ! [-] Fractional height of this point within the ice
        real(prec), intent(IN) :: T0     ! [K] Reference freezing point of water (e.g., 273.15 K or 0 C)
        real(prec), intent(IN) :: beta   ! [K Pa^-1] Melting point gradient with pressure
        real(prec) :: T_pmp              ! [K] Pressure corrected melting point

        ! Local variables
        real(prec) :: depth

!         real(prec), parameter :: beta = 9.8e-8 [K Pa^-1]      ! Greve and Blatter (2009) 
!         real(prec), parameter :: beta = 9.7e-8 [K Pa^-1]      ! EISMINT2 value (beta1 = 8.66e-4 [K m^-1])
!         real(prec), parameter :: beta = 7.9e-8 [K Pa^-1]      ! Kleiner et al. (2015)

        ! Get thickness of ice above current point
        depth = H_ice*(1.0-zeta)

        ! Calculate the pressure-corrected melting point
        T_pmp = T0 - (beta*rho_ice*g)*depth
        
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

        if (f_grnd .eq. 0.0_prec) then
            ! Floating points are temperate by default
            f_pmp = 1.0_prec 

        else 
            ! Calculate the fraction at the pressure melting point 

            if (gamma .eq. 0.0_prec) then
                ! No decay function, binary pmp fraction

                if (T_ice .ge. T_pmp) then 
                    f_pmp = 1.0_prec
                else 
                    f_pmp = 0.0_prec 
                end if 

            else

                ! Apply decay function 
                f_pmp = min(1.0_prec, exp((T_ice-T_pmp)/gamma) )

                ! Ensure pure values of 0.0 and 1.0 beyond a threshold 
                if (f_pmp .lt. 1e-2)        f_pmp = 0.0_prec 
                if (f_pmp .gt. (1.0-1e-2))  f_pmp = 1.0_prec 

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
    
    subroutine define_temp_linear_3D(enth,T_ice,omega,cp,H_ice,T_srf,zeta_aa)
        ! Define a linear vertical temperature profile
       
        implicit none

        real(prec), intent(OUT) :: enth(:,:,:)          ! [J m-3] Enthalpy 
        real(prec), intent(OUT) :: T_ice(:,:,:)         ! [K] Temperature
        real(prec), intent(OUT) :: omega(:,:,:)         ! [--] Water content
        real(prec), intent(IN)  :: cp(:,:,:)            ! Heat capacity
        real(prec), intent(IN)  :: H_ice(:,:)
        real(prec), intent(IN)  :: T_srf(:,:) 
        real(prec), intent(IN)  :: zeta_aa(:)
           
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
                 
                T_base       = calc_T_pmp(H_ice(i,j),zeta_aa(1),T0,T_pmp_beta) - 10.0 
                T_ice(i,j,:) = define_temp_linear_column(T_srf(i,j),T_base,T0,zeta_aa)

            else 
                ! No ice present, set equal to T_pmp (ie, T0)
                T_ice(i,j,:) = T0

            end if 

            ! Assume zero water content 
            omega = 0.0_wp 

            ! Calculate enthalpy too
            call convert_to_enthalpy(enth,T_ice,omega,T_pmp,cp,L_ice)
        
        end do 
        end do  
        
        return 

    end subroutine define_temp_linear_3D

    function define_temp_linear_column(T_srf,T_base,T0,zeta_aa) result(T_ice)

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

    end function define_temp_linear_column

    subroutine define_temp_robin_3D (enth,T_ice,omega,T_pmp,cp,ct,Q_rock,T_srf,H_ice,H_w,smb,bmb,f_grnd,zeta_aa,cold)
        ! Robin solution for thermodynamics for a given column of ice 
        ! Note zeta=height, k=1 base, k=nz surface 

        implicit none 

        real(prec), intent(OUT) :: enth(:,:,:)      ! [J m-3] Enthalpy 
        real(prec), intent(OUT) :: T_ice(:,:,:)     ! [K] Temperature
        real(prec), intent(OUT) :: omega(:,:,:)     ! [--] Water content
        real(prec), intent(IN)  :: T_pmp(:,:,:)     ! [K] Pressure melting point temp.
        real(prec), intent(IN)  :: cp(:,:,:)        ! [J kg-1 K-1] Specific heat capacity
        real(prec), intent(IN)  :: ct(:,:,:)        ! [J a-1 m-1 K-1] Heat conductivity 
        real(prec), intent(IN)  :: Q_rock(:,:)      ! [mW m-2] Bedrock surface heat flux 
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

            T_ice(i,j,:) = define_temp_robin_column(zeta_aa,T_pmp(i,j,:),ct(i,j,:),cp(i,j,:),rho_ice,H_ice(i,j), &
                                                  T_srf(i,j),smb(i,j)+bmb(i,j),Q_rock(i,j),is_float)

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

        ! Assume zero water content 
        omega = 0.0_wp 

        ! Calculate enthalpy too
        call convert_to_enthalpy(enth,T_ice,omega,T_pmp,cp,L_ice)

        return 

    end subroutine define_temp_robin_3D

    function define_temp_robin_column(zeta_aa,T_pmp,kt,cp,rho_ice,H_ice,T_srf,mb_net,Q_rock,is_float) result(T_ice)
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
        real(prec), intent(IN) :: Q_rock 
        logical,    intent(IN) :: is_float 

        real(prec) :: T_ice(size(zeta_aa,1))

        ! Local variables 
        integer    :: k, nz_aa 
        real(prec) :: dTdz_b, z, kappa, ll   

        real(prec), parameter :: sqrt_pi   = sqrt(pi) 
        real(prec), parameter :: T_ocn     = 271.15   ! [K]
        real(prec), parameter :: H_ice_min = 0.1      ! [m] Minimum ice thickness to calculate Robin solution 
        real(prec), parameter :: mb_net_min = 1e-2    ! [m a-1] Minimum allowed net mass balance for stability
        real(prec) :: Q_rock_now, mb_now  

        nz_aa = size(T_ice,1) 
        
        Q_rock_now = Q_rock *1e-3*sec_year    ! [mW m-2] => [J a-1 m-2]

        ! Calculate temperature gradient at base 
        dTdz_b = -Q_rock_now/kt(1) 

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

    end function define_temp_robin_column

    subroutine define_temp_bedrock_3D(enth_rock,T_rock,Q_rock,cp_rock,kt_rock,Q_geo,T_bed,H_rock,zeta_aa)
        ! Define 3D bedrock enth/temp field 

        implicit none 

        real(prec), intent(OUT) :: enth_rock(:,:,:)     ! [J m-3] 3D enthalpy field
        real(prec), intent(OUT) :: T_rock(:,:,:)        ! [K] 3D temperature field
        real(prec), intent(OUT) :: Q_rock(:,:)          ! [mW m-2] Bed surface heat flux 
        real(prec), intent(IN)  :: cp_rock              ! [J kg-1 K-1] Specific heat capacity
        real(prec), intent(IN)  :: kt_rock              ! [J a-1 m-1 K-1] Heat conductivity 
        real(prec), intent(IN)  :: Q_geo(:,:)           ! [mW m-2] Geothermal heat flux 
        real(prec), intent(IN)  :: T_bed(:,:)           ! [K] Surface temperature 
        real(prec), intent(IN)  :: H_rock               ! [m] Column thickness 
        real(prec), intent(IN)  :: zeta_aa(:)           ! [--] Vertical zeta coordinates (zeta==height), aa-nodes

        ! Local variable
        integer :: i, j, k, nx, ny, nz_aa  
        
        nx    = size(T_rock,1)
        ny    = size(T_rock,2)
        nz_aa = size(zeta_aa)

        do j = 1, ny
        do i = 1, nx 

            ! Calculate temperature profile 
            call define_temp_bedrock_column(T_rock(i,j,:),kt_rock,rho_rock,H_rock, &
                                                    T_bed(i,j),Q_geo(i,j),zeta_aa)

            ! Calculate heat flux through bed surface from lithosphere [mW m-2]
            call calc_Q_bedrock_column(Q_rock(i,j),T_rock(i,j,:),kt_rock,H_rock,zeta_aa)

        end do 
        end do 

        ! Get enthalpy too 
        call convert_to_enthalpy(enth_rock,T_rock,0.0_wp,0.0_wp,cp_rock,0.0_wp)
          
        return 

    end subroutine define_temp_bedrock_3D

    subroutine define_temp_bedrock_column(T_rock,kt_rock,rho_rock,H_rock,T_bed,Q_geo,zeta_aa)
        ! This function will impose a temperature profile in a column 
        ! of bedrock assuming equilibrium with the bed surface temperature (T_bed)
        ! and the geothermal heat flux deep in the bedrock (Q_geo) 

        implicit none 

        real(prec), intent(OUT) :: T_rock(:) 
        real(prec), intent(IN)  :: kt_rock 
        real(prec), intent(IN)  :: rho_rock 
        real(prec), intent(IN)  :: H_rock
        real(prec), intent(IN)  :: T_bed 
        real(prec), intent(IN)  :: Q_geo 
        real(prec), intent(IN)  :: zeta_aa(:) 
        
        ! Local variables
        integer :: k, nz_aa
        real(prec) :: Q_geo_now 
        real(prec) :: dTdz
        
        nz_aa = size(zeta_aa,1)

        ! Get geothermal heat flux in proper units 
        Q_geo_now = Q_geo*1e-3*sec_year   ! [mW m-2] => [J m-2 a-1]

        ! Temperature at the bed surface is equal to the boundary temperature
        T_rock(nz_aa) = T_bed 

        ! Descend into bedrock and linear gradient
        do k = nz_aa-1, 1, -1 

            ! Determine the slope based on equilibrium with geothermal heat flux 
            dTdz = -Q_geo_now / kt_rock

            ! Calculate the temperature of the lithosphere at each layer 
            T_rock(k) = T_rock(k+1) - dTdz*(H_rock*(zeta_aa(k+1)-zeta_aa(k)))

        end do 

        return 

    end subroutine define_temp_bedrock_column

    subroutine calc_Q_bedrock(Q_rock,T_rock,kt_rock,H_rock,zeta_aa)

        implicit none 

        real(wp), intent(OUT) :: Q_rock(:,:) 
        real(wp), intent(IN)  :: T_rock(:,:,:) 
        real(wp), intent(IN)  :: kt_rock
        real(wp), intent(IN)  :: H_rock 
        real(wp), intent(IN)  :: zeta_aa(:) 

        ! Local variables 
        integer  :: i, j, nx, ny  

        nx    = size(Q_rock,1)
        ny    = size(Q_rock,2) 

        do j = 1, ny 
        do i = 1, nx 

            call calc_Q_bedrock_column(Q_rock(i,j),T_rock(i,j,:),kt_rock,H_rock,zeta_aa)

        end do 
        end do  

        ! [J a-1 m-2] => [mW m-2] (same units as Q_geo by default)
        Q_rock = Q_rock *1e3 / sec_year 

        return 

    end subroutine calc_Q_bedrock

    subroutine calc_Q_bedrock_column(Q_rock,T_rock,kt_rock,H_rock,zeta_aa)

        implicit none 

        real(wp), intent(OUT) :: Q_rock 
        real(wp), intent(IN)  :: T_rock(:) 
        real(wp), intent(IN)  :: kt_rock
        real(wp), intent(IN)  :: H_rock
        real(wp), intent(IN)  :: zeta_aa(:) 

        ! Local variables 
        integer  :: nz_aa 
        real(wp) :: dz 

        nz_aa = size(zeta_aa,1) 

        ! Determine layer thickness of lithosphere at bed surface
        dz = H_rock*(zeta_aa(nz_aa)-zeta_aa(nz_aa-1))

        ! Calculate lithospheric heat flux [J a-1 m-2]
        Q_rock = -kt_rock * (T_rock(nz_aa)-T_rock(nz_aa-1)) / dz 

        ! [J a-1 m-2] => [mW m-2] (same units as Q_geo by default)
        Q_rock = Q_rock *1e3 / sec_year 

        return 

    end subroutine calc_Q_bedrock_column
    
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

    ! ==== BASAL HYDROLOGY PHYSICS ===============================

    subroutine calc_basal_water_local(H_w,dHwdt,H_ice,bmb_w,f_grnd,dt,till_rate,H_w_max)
        ! Calculate the basal water layer thickness based on a simple local 
        ! water balance: dHw/dt = bmb_w - till_rate
        implicit none 
         
        real(prec), intent(INOUT) :: H_w(:,:)         ! [m] Water layer thickness
        real(prec), intent(INOUT) :: dHwdt(:,:)       ! [m/a] Water layer thickness change
        real(prec), intent(IN)    :: H_ice(:,:)       ! [m] Ice thickness 
        real(prec), intent(IN)    :: bmb_w(:,:)       ! [m/a] Basal water mass balance
        real(prec), intent(IN)    :: f_grnd(:,:)      ! [-] Grounded fraction
        real(prec), intent(IN)    :: dt               ! [a] Timestep 
        real(prec), intent(IN)    :: till_rate        ! [m/a] Till drainage rate 
        real(prec), intent(IN)    :: H_w_max          ! [m] Maximum allowed water depth 

        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: im1, ip1, jm1, jp1 
        
        nx = size(H_ice,1)
        ny = size(H_ice,2)

        ! Store initial H_w field  
        dHwdt   = H_w 

        where (f_grnd .gt. 0.0 .and. H_ice .gt. 0.0)
            ! Grounded ice point

            ! Update mass balance of H_w
            H_w = H_w + dt*(bmb_w-till_rate)

            ! Restrict H_w to values within limits
            H_w = max(H_w,0.0)
            H_w = min(H_w,H_w_max)

        else where (f_grnd .gt. 0.0) 
            ! Ice-free land above sea level 

            H_w = 0.0 

        elsewhere
            ! Set water layer thickness to maximum layer thickness

            H_w = H_w_max 

        end where 

        ! Additionally set points at the grounding line to
        ! the maximum water thickness 
        do j = 1, ny 
        do i = 1, nx

            im1 = max(1, i-1)
            ip1 = min(nx,i+1)
            
            jm1 = max(1, j-1)
            jp1 = min(ny,j+1)

            ! Grounded point or partially floating point with floating neighbors
            if (H_ice(i,j) .gt. 0.0 .and. f_grnd(i,j) .gt. 0.0 .and. &
                (f_grnd(im1,j) .eq. 0.0 .or. f_grnd(ip1,j) .eq. 0.0 .or. &
                 f_grnd(i,jm1) .eq. 0.0 .or. f_grnd(i,jp1) .eq. 0.0) ) then 
                
                H_w(i,j) = H_w_max

            end if 

        end do 
        end do  

        ! Determine rate of change 
        if (dt .ne. 0.0_prec) then 
            dHwdt   = (dHwdt - H_w) / dt
        else 
            dHwdt   = 0.0_prec
        end if 

        return 

    end subroutine calc_basal_water_local

    ! ========== ENTHALPY ==========================================

    elemental subroutine convert_to_enthalpy(enth,temp,omega,T_pmp,cp,L)
        ! Given temperature and water content, calculate enthalpy.

        implicit none 

        real(prec), intent(OUT) :: enth             ! [J m-3] Enthalpy 
        real(prec), intent(IN)  :: temp             ! [K] Temperature 
        real(prec), intent(IN)  :: omega            ! [-] Water content (fraction)
        real(prec), intent(IN)  :: T_pmp            ! [K] Pressure melting point
        real(prec), intent(IN)  :: cp               ! [J kg-1 K-1] Heat capacity 
        real(prec), intent(IN)  :: L                ! [J kg-1] Latent heat of fusion 
        
        enth = (1.0_prec-omega)*(cp*temp) + omega*(cp*T_pmp + L)

        return 

    end subroutine convert_to_enthalpy

    subroutine convert_from_enthalpy_column(enth,temp,omega,T_pmp,cp,L)
        ! Given enthalpy, calculate temperature and water content. 

        implicit none 

        real(prec), intent(INOUT) :: enth(:)            ! [J m-3] Enthalpy, nz_aa nodes
        real(prec), intent(OUT)   :: temp(:)            ! [K] Temperature, nz_aa nodes  
        real(prec), intent(OUT)   :: omega(:)           ! [-] Water content (fraction), nz_aa nodes 
        real(prec), intent(IN)    :: T_pmp(:)           ! [K] Pressure melting point, nz_aa nodes 
        real(prec), intent(IN)    :: cp(:)              ! [J kg-1 K-1] Heat capacity,nz_aa nodes 
        real(prec), intent(IN)    :: L                  ! [J kg-1] Latent heat of fusion
        
        ! Local variables
        integer    :: k, nz_aa  
        real(prec), allocatable :: enth_pmp(:)  

        nz_aa = size(enth,1)

        allocate(enth_pmp(nz_aa))

        ! Find pressure melting point enthalpy
        enth_pmp = T_pmp * cp 

        ! Column interior and basal layer
        ! Note: although the k=1 is a boundary value with no thickness,
        ! allow it to retain omega to maintain consistency with grid points above.
        do k = 1, nz_aa-1

            if (enth(k) .gt. enth_pmp(k)) then
                ! Temperate ice 
                
                temp(k)  = T_pmp(k)
                omega(k) = (enth(k) - enth_pmp(k)) / L 
             else
                ! Cold ice 

                temp(k)  = enth(k) / cp(k) 
                omega(k) = 0.0_prec

             end if

        end do 

        ! Surface layer 
        if (enth(nz_aa) .ge. enth_pmp(nz_aa)) then 
            ! Temperate surface, reset omega to zero and enth to pmp value 
            
            enth(nz_aa)  = enth_pmp(nz_aa)
            temp(nz_aa)  = enth(nz_aa) / cp(nz_aa)
            omega(nz_aa) = 0.0_prec 
        
        else 
            ! Cold surface, calculate T, and reset omega to zero 
            
            temp(nz_aa)  = enth(nz_aa) / cp(nz_aa)
            omega(nz_aa) = 0.0_prec 
        
        end if 
        
        return 

    end subroutine convert_from_enthalpy_column

end module thermodynamics
