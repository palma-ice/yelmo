
module yelmo_thermodynamics

    use nml 
    use yelmo_defs 
    use yelmo_tools, only : smooth_gauss_2D, smooth_gauss_3D, gauss_values, fill_borders_3D, &
            stagger_aa_ab, regularize2D
    
    use thermodynamics 
    use ice_enthalpy
    use solver_advection, only : calc_advec2D  

    implicit none
    
    private
    public :: calc_ytherm 
    public :: ytherm_par_load, ytherm_alloc, ytherm_dealloc
    
contains

    subroutine calc_ytherm(thrm,tpo,dyn,mat,bnd,time)

        implicit none
        
        type(ytherm_class), intent(INOUT) :: thrm
        type(ytopo_class),  intent(IN)    :: tpo 
        type(ydyn_class),   intent(IN)    :: dyn
        type(ymat_class),   intent(IN)    :: mat
        type(ybound_class), intent(IN)    :: bnd  
        real(prec),         intent(IN)    :: time  

        ! Local variables 
        integer :: i, j, k, nx, ny  
        real(prec) :: dt 

        nx = thrm%par%nx
        ny = thrm%par%ny

        ! Initialize time if necessary 
        if (thrm%par%time .gt. time) then 
            thrm%par%time = time
        end if 

        ! Get time step and advance current time 
        dt            = time - thrm%par%time 
        thrm%par%time = time 
        
        ! === Determine some thermal properties === 

        ! Calculate the specific heat capacity of the ice
        if (thrm%par%use_const_cp) then 
            thrm%now%cp = thrm%par%const_cp 
        else  
            thrm%now%cp = calc_specific_heat_capacity(thrm%now%T_ice)
        end if 
        
        ! Calculate the heat conductivity of the ice
        if (thrm%par%use_const_kt) then 
            thrm%now%kt = thrm%par%const_kt 
        else  
            thrm%now%kt = calc_thermal_conductivity(thrm%now%T_ice)
        end if 

        ! Calculate the pressure-corrected melting point (in Kelvin)
        do k = 1, thrm%par%nz_aa
            thrm%now%T_pmp(:,:,k) = calc_T_pmp(tpo%now%H_ice,thrm%par%zeta_aa(k),T0,T_pmp_beta)
        end do 

        ! Calculate internal strain heating
        if (thrm%par%use_strain_sia) then 
            ! Calculate strain heating from SIA approximation

            call calc_strain_heating_sia(thrm%now%Q_strn,dyn%now%ux,dyn%now%uy,tpo%now%dzsdx,tpo%now%dzsdy, &
                                      thrm%now%cp,tpo%now%H_ice,rho_ice,thrm%par%zeta_aa,thrm%par%zeta_ac)
        
        else
            ! Calculate strain heating from strain rate tensor and viscosity (general approach)
            
            call calc_strain_heating(thrm%now%Q_strn,mat%now%strn%de,mat%now%visc,thrm%now%cp,rho_ice)

        end if 
        
        ! Smooth strain heating 
        if (thrm%par%n_sm_qstrn .gt. 0) then 
            call smooth_gauss_3D(thrm%now%Q_strn,tpo%now%H_ice.gt.0.0,thrm%par%dx,thrm%par%n_sm_qstrn, &
                                    tpo%now%H_ice.gt.0.0)
        end if 
        
        ! Calculate the basal frictional heating 
        call calc_basal_heating(thrm%now%Q_b,dyn%now%ux_b,dyn%now%uy_b,dyn%now%taub_acx,dyn%now%taub_acy, &
                                tpo%now%H_ice,thrm%now%T_prime_b,gamma=2.0_prec)

        ! Ensure basal frictional heating is relatively smooth
        call regularize2D(thrm%now%Q_b,tpo%now%H_ice,tpo%par%dx)
            
        ! Smooth basal frictional heating 
        if (thrm%par%n_sm_qb .gt. 0) then 
            call smooth_gauss_2D(thrm%now%Q_b,tpo%now%H_ice.gt.0.0,thrm%par%dx,thrm%par%n_sm_qb, &
                                    tpo%now%H_ice.gt.0.0)
        end if 

        if ( dt .gt. 0.0 ) then     
            ! Ice thermodynamics should evolve, perform calculations 

            select case(trim(thrm%par%method))

                case("enth","temp") 
                    ! Perform enthalpy/temperature solving via advection-diffusion equation
                    
                    call calc_ytherm_enthalpy_3D(thrm%now%enth,thrm%now%T_ice,thrm%now%omega,thrm%now%bmb_grnd,thrm%now%Q_ice_b, &
                                thrm%now%H_cts,thrm%now%T_pmp,thrm%now%cp,thrm%now%kt,dyn%now%ux,dyn%now%uy,dyn%now%uz,thrm%now%Q_strn, &
                                thrm%now%Q_b,bnd%Q_geo,bnd%T_srf,tpo%now%H_ice,thrm%now%H_w,thrm%now%dHwdt,tpo%now%H_grnd,tpo%now%f_grnd, &
                                thrm%par%zeta_aa,thrm%par%zeta_ac,thrm%par%dzeta_a,thrm%par%dzeta_b,thrm%par%enth_cr,thrm%par%omega_max, &
                                dt,thrm%par%dx,thrm%par%method,thrm%par%solver_advec)
                    
                case("robin")
                    ! Use Robin solution for ice temperature 

                    call define_temp_robin_3D(thrm%now%T_ice,thrm%now%T_pmp,thrm%now%cp,thrm%now%kt, &
                                       bnd%Q_geo,bnd%T_srf,tpo%now%H_ice,thrm%now%H_w,bnd%smb, &
                                       thrm%now%bmb_grnd,tpo%now%f_grnd,thrm%par%zeta_aa,cold=.FALSE.)

                    ! Also populate enthalpy 
                    call convert_to_enthalpy(thrm%now%enth,thrm%now%T_ice,thrm%now%omega,thrm%now%T_pmp, &
                                            thrm%now%cp,L_ice)

                case("robin-cold")
                    ! Use Robin solution for ice temperature averaged with cold linear profile
                    ! to ensure cold ice at the base

                    call define_temp_robin_3D(thrm%now%T_ice,thrm%now%T_pmp,thrm%now%cp,thrm%now%kt, &
                                       bnd%Q_geo,bnd%T_srf,tpo%now%H_ice,thrm%now%H_w,bnd%smb, &
                                       thrm%now%bmb_grnd,tpo%now%f_grnd,thrm%par%zeta_aa,cold=.TRUE.)

                    ! Also populate enthalpy 
                    call convert_to_enthalpy(thrm%now%enth,thrm%now%T_ice,thrm%now%omega,thrm%now%T_pmp, &
                                            thrm%now%cp,L_ice)

                case("linear")
                    ! Use linear solution for ice temperature

                    ! Calculate the ice temperature (eventually water content and enthalpy too)
                    call define_temp_linear_3D(thrm%now%T_ice,thrm%par%zeta_aa,tpo%now%H_ice,bnd%T_srf)

                    ! Also populate enthalpy 
                    call convert_to_enthalpy(thrm%now%enth,thrm%now%T_ice,thrm%now%omega,thrm%now%T_pmp, &
                                            thrm%now%cp,L_ice)

                case("fixed") 
                    ! Pass - do nothing, use the temperature field as it is defined

                case DEFAULT 

                    write(*,*) "ytherm:: Error: thermodynamics option not recognized: method = ", trim(thrm%par%method)
                    stop 

            end select 

!             if (trim(thrm%par%basal_water_method) .eq. "local") then 
!                 ! If updating basal water calculate it here...
                
                ! Update basal water layer thickness
                call calc_basal_water_local(thrm%now%H_w,thrm%now%dHwdt,tpo%now%H_ice,-thrm%now%bmb_grnd*(rho_ice/rho_w), &
                                        tpo%now%f_grnd,dt,thrm%par%till_rate,thrm%par%H_w_max)
            
!             end if 

        end if 

        ! Calculate homologous temperature at the base 
        thrm%now%T_prime_b = thrm%now%T_ice(:,:,1) - thrm%now%T_pmp(:,:,1) 
        
        ! Calculate gridpoint fraction at the pressure melting point
        thrm%now%f_pmp = calc_f_pmp(thrm%now%T_ice(:,:,1),thrm%now%T_pmp(:,:,1),thrm%par%gamma,tpo%now%f_grnd)

            
!         if (yelmo_log) then 
!             if (count(tpo%now%H_ice.gt.0.0) .gt. 0) then 
!                 write(*,"(a,f14.4,f10.4,f10.2)") "calc_ytherm:: time = ", thrm%par%time, dt, &
!                     sum(thrm%now%T_ice(:,:,thrm%par%nz_aa),mask=tpo%now%H_ice.gt.0.0)/real(count(tpo%now%H_ice.gt.0.0))
!             else 
!                 write(*,"(a,f14.4,f10.4,f10.2)") "calc_ytherm:: time = ", thrm%par%time, dt, 0.0 
!             end if 
!         end if 

        return

    end subroutine calc_ytherm

    subroutine calc_ytherm_enthalpy_3D(enth,T_ice,omega,bmb_grnd,Q_ice_b,H_cts,T_pmp,cp,kt,ux,uy,uz,Q_strn,Q_b,Q_geo, &
                            T_srf,H_ice,H_w,dHwdt,H_grnd,f_grnd,zeta_aa,zeta_ac,dzeta_a,dzeta_b,cr,omega_max,dt,dx,solver,solver_advec)
        ! This wrapper subroutine breaks the thermodynamics problem into individual columns,
        ! which are solved independently by calling calc_enth_column

        ! Note zeta=height, k=1 base, k=nz surface 
        
        implicit none 

        real(prec), intent(INOUT) :: enth(:,:,:)    ! [J m-3] Ice enthalpy
        real(prec), intent(INOUT) :: T_ice(:,:,:)   ! [K] Ice column temperature
        real(prec), intent(INOUT) :: omega(:,:,:)   ! [--] Ice water content
        real(prec), intent(INOUT) :: bmb_grnd(:,:)  ! [m a-1] Basal mass balance (melting is negative)
        real(prec), intent(OUT)   :: Q_ice_b(:,:)   ! [J a-1 m-2] Basal ice heat flux 
        real(prec), intent(OUT)   :: H_cts(:,:)     ! [m] Height of the cold-temperate transition surface (CTS)
        real(prec), intent(IN)    :: T_pmp(:,:,:)   ! [K] Pressure melting point temp.
        real(prec), intent(IN)    :: cp(:,:,:)      ! [J kg-1 K-1] Specific heat capacity
        real(prec), intent(IN)    :: kt(:,:,:)      ! [J a-1 m-1 K-1] Heat conductivity 
        real(prec), intent(IN)    :: ux(:,:,:)      ! [m a-1] Horizontal x-velocity 
        real(prec), intent(IN)    :: uy(:,:,:)      ! [m a-1] Horizontal y-velocity 
        real(prec), intent(IN)    :: uz(:,:,:)      ! [m a-1] Vertical velocity 
        real(prec), intent(IN)    :: Q_strn(:,:,:)  ! [K a-1] Internal strain heat production in ice
        real(prec), intent(IN)    :: Q_b(:,:)       ! [J a-1 m-2] Basal frictional heat production 
        real(prec), intent(IN)    :: Q_geo(:,:)     ! [mW m-2] Geothermal heat flux 
        real(prec), intent(IN)    :: T_srf(:,:)     ! [K] Surface temperature 
        real(prec), intent(IN)    :: H_ice(:,:)     ! [m] Ice thickness 
        real(prec), intent(IN)    :: H_w(:,:)       ! [m] Basal water layer thickness 
        real(prec), intent(IN)    :: dHwdt(:,:)     ! [m/a] Basal water layer thickness change
        real(prec), intent(IN)    :: H_grnd(:,:)    ! [--] Ice thickness above flotation 
        real(prec), intent(IN)    :: f_grnd(:,:)    ! [--] Grounded fraction
        real(prec), intent(IN)    :: zeta_aa(:)     ! [--] Vertical sigma coordinates (zeta==height), aa-nodes
        real(prec), intent(IN)    :: zeta_ac(:)     ! [--] Vertical sigma coordinates (zeta==height), ac-nodes
        real(prec), intent(IN)    :: dzeta_a(:)     ! d Vertical height axis (0:1) 
        real(prec), intent(IN)    :: dzeta_b(:)     ! d Vertical height axis (0:1) 
        real(prec), intent(IN)    :: cr             ! [--] Conductivity ratio for temperate ice (kappa_temp = enth_cr*kappa_cold)
        real(prec), intent(IN)    :: omega_max      ! [--] Maximum allowed water content fraction 
        real(prec), intent(IN)    :: dt             ! [a] Time step 
        real(prec), intent(IN)    :: dx             ! [a] Horizontal grid step 
        character(len=*), intent(IN) :: solver      ! "enth" or "temp" 
        character(len=*), intent(IN) :: solver_advec    ! "expl" or "impl-upwind"

        ! Local variables
        integer :: i, j, k, nx, ny, nz_aa, nz_ac  
        real(prec), allocatable  :: advecxy(:)   ! [K a-1 m-2] Horizontal heat advection 
        real(prec) :: T_shlf, H_grnd_lim, f_scalar, T_base  
        real(prec) :: H_ice_now 

        real(prec), allocatable :: T_ice_old(:,:,:) 
        real(prec), allocatable :: enth_old(:,:,:) 
        real(prec), allocatable :: advecxy3D(:,:,:)
        real(prec) :: filter0(3,3), filter(3,3) 

        real(prec), parameter :: H_ice_thin = 15.0   ! [m] Threshold to define 'thin' ice

        nx    = size(T_ice,1)
        ny    = size(T_ice,2)
        nz_aa = size(zeta_aa,1)
        nz_ac = size(zeta_ac,1)

        allocate(advecxy(nz_aa))
        allocate(T_ice_old(nx,ny,nz_aa))
        allocate(enth_old(nx,ny,nz_aa))
        allocate(advecxy3D(nx,ny,nz_aa))
        
        ! First perform horizontal advection (this doesn't work properly, 
        ! use column-based upwind horizontal advection below)
        !call calc_enth_horizontal_advection_3D(T_ice,ux,uy,H_ice,dx,dt,solver_advec)
        
        ! Diagnose horizontal advection 
        advecxy3D = 0.0 
!         do k = 2, nz_aa-1
!             call calc_adv2D_impl_upwind_rate(advecxy3D(:,:,k),T_ice(:,:,k),ux(:,:,k),uy(:,:,k),H_ice*0.0,dx,dx,dt,f_upwind=1.0)
!         end do 

        ! Store original ice temperature field here for input to horizontal advection
        ! calculations 
        T_ice_old = T_ice 
        enth_old  = enth 

        ! Initialize gaussian filter kernel 
        filter0 = gauss_values(dx,dx,sigma=2.0*dx,n=size(filter,1))

        do j = 3, ny-2
        do i = 3, nx-2 
            
            ! For floating points, calculate the approximate marine-shelf temperature 
            ! ajr, later this should come from an external model, and T_shlf would
            ! be the boundary variable directly
            if (f_grnd(i,j) .lt. 1.0) then 

                ! Calculate approximate marine freezing temp, limited to pressure melting point 
                T_shlf = calc_T_base_shlf_approx(H_ice(i,j),T_pmp(i,j,1),H_grnd(i,j))

            else 
                ! Assigned for safety 

                T_shlf   = T_pmp(i,j,1)

            end if 

            if (H_ice(i,j) .le. H_ice_thin) then 
                ! Ice is too thin or zero: prescribe linear temperature profile
                ! between temperate ice at base and surface temperature 
                ! (accounting for floating/grounded nature via T_base)

                if (f_grnd(i,j) .lt. 1.0) then 
                    ! Impose T_shlf for the basal temperature
                    T_base = T_shlf 
                else 
                    ! Impose the pressure melting point of grounded ice 
                    T_base = T_pmp(i,j,1) 
                end if 

                T_ice(i,j,:) = calc_temp_linear_column(T_srf(i,j),T_base,T_pmp(i,j,nz_aa),zeta_aa)

            else 
                ! Thick ice exists, call thermodynamic solver for the column

                ! No filtering of H_ice, take actual value
!                 H_ice_now = H_ice(i,j) 
                
                ! Filter everywhere
!                 filter = filter0
!                 H_ice_now = sum(H_ice(i-1:i+1,j-1:j+1)*filter)

                ! Filter everywhere there is ice 
                ! filter = filter0 
                ! where (H_ice(i-1:i+1,j-1:j+1) .eq. 0.0) filter = 0.0 
                ! filter = filter/sum(filter)
                ! H_ice_now = sum(H_ice(i-1:i+1,j-1:j+1)*filter)
                
                ! Filter at the margin only 
                if (count(H_ice(i-1:i+1,j-1:j+1) .eq. 0.0) .ge. 2) then
                    filter = filter0 
                    where (H_ice(i-1:i+1,j-1:j+1) .eq. 0.0) filter = 0.0 
                    filter = filter/sum(filter)
                    H_ice_now = sum(H_ice(i-1:i+1,j-1:j+1)*filter)
                else 
                    H_ice_now = H_ice(i,j) 
                end if 
                
                ! Pre-calculate the contribution of horizontal advection to column solution
                ! (use unmodified T_ice_old field as input, to avoid mixing with new solution)
!                 call calc_advec_horizontal_column(advecxy,T_ice_old,H_ice,ux,uy,dx,i,j)
!                 call calc_advec_horizontal_column_quick(advecxy,T_ice_old,H_ice,ux,uy,dx,i,j)
!                 do k = 1, nz_aa
!                     call calc_adv2D_expl_rate(advecxy(k),T_ice_old(:,:,k),ux(:,:,k),uy(:,:,k),dx,dx,i,j)
!                 end do 
                !advecxy = advecxy3D(i,j,:)
                !advecxy = 0.0_prec 
!                 write(*,*) "advecxy: ", i,j, maxval(abs(advecxy3D(i,j,:)-advecxy))

                !call calc_advec_horizontal_column(advecxy,enth_old,H_ice,ux,uy,dx,i,j)

                do k = 1, nz_aa
                    call calc_adv2D_expl_rate(advecxy(k),enth_old(:,:,k),ux(:,:,k),uy(:,:,k),dx,dx,i,j)
                end do 

                call calc_enth_column(enth(i,j,:),T_ice(i,j,:),omega(i,j,:),bmb_grnd(i,j),Q_ice_b(i,j),H_cts(i,j), &
                        T_pmp(i,j,:),cp(i,j,:),kt(i,j,:),advecxy,uz(i,j,:),Q_strn(i,j,:),Q_b(i,j),Q_geo(i,j),T_srf(i,j), &
                        T_shlf,H_ice_now,H_w(i,j),f_grnd(i,j),zeta_aa,zeta_ac,dzeta_a,dzeta_b,cr,omega_max,T0,dt)
                
            end if 

        end do 
        end do 

        ! Fill in borders 
        call fill_borders_3D(enth,nfill=2)
        call fill_borders_3D(T_ice,nfill=2)
        call fill_borders_3D(omega,nfill=2)

        return 

    end subroutine calc_ytherm_enthalpy_3D

    subroutine calc_enth_horizontal_advection_3D(T_ice,ux,uy,H_ice,dx,dt,solver)

        implicit none 

        real(prec),       intent(INOUT) :: T_ice(:,:,:)         ! [K]   Ice temperature/enthalpy, aa-nodes  
        real(prec),       intent(IN)    :: ux(:,:,:)            ! [m/a] 2D velocity, x-direction (ac-nodes)
        real(prec),       intent(IN)    :: uy(:,:,:)            ! [m/a] 2D velocity, y-direction (ac-nodes)
        real(prec),       intent(IN)    :: H_ice(:,:)           ! [m]   Ice thickness 
        real(prec),       intent(IN)    :: dx                   ! [m]   Horizontal resolution
        real(prec),       intent(IN)    :: dt                   ! [a]   Timestep 
        character(len=*), intent(IN)    :: solver               ! Solver to use for the ice thickness advection equation

        ! Local variables 
        integer :: i, j, k, n, nx, ny, nz 
        real(prec), allocatable :: T_dot(:,:) 

        nx = size(T_ice,1)
        ny = size(T_ice,2)
        nz = size(T_ice,3) 

        allocate(T_dot(nx,ny)) 
        T_dot = 0.0_prec 

        ! First populate boundary values so that T_ice next to ice sheet is equal to 
        ! ice sheet. 
        do j = 2, ny-1 
        do i = 2, nx-1 

            if (H_ice(i,j) .eq. 0.0 .and. &
                count([H_ice(i-1,j),H_ice(i+1,j),H_ice(i,j-1),H_ice(i,j+1)] .gt. 0.0) .gt. 0) then 
                ! Apply to ice-free points with ice-neighbors only 

                T_ice(i,j,:) = 0.0_prec 
                n = 0 

                if (H_ice(i-1,j) .gt. 0.0) then 
                    T_ice(i,j,:) = T_ice(i,j,:) + T_ice(i-1,j,:)
                    n = n+1 
                end if 
                
                if (H_ice(i+1,j) .gt. 0.0) then 
                    T_ice(i,j,:) = T_ice(i,j,:) + T_ice(i+1,j,:)
                    n = n+1 
                end if 
                
                if (H_ice(i,j-1) .gt. 0.0) then 
                    T_ice(i,j,:) = T_ice(i,j,:) + T_ice(i,j-1,:)
                    n = n+1 
                end if 
                
                if (H_ice(i,j+1) .gt. 0.0) then 
                    T_ice(i,j,:) = T_ice(i,j,:) + T_ice(i,j+1,:)
                    n = n+1 
                end if 
                
                if (n .gt. 0) then 
                    ! Get average 
                    T_ice(i,j,:) = T_ice(i,j,:) / real(n,prec) 
                else 
                    write(*,*) "calc_enth_horizontal_advection_3D:: error: something went wrong!"
                    stop 
                end if 

            end if 

        end do 
        end do 

        ! Resolve horizontal advection layer by layer 

        do k = 2, nz-1    
            call calc_advec2D(T_ice(:,:,k),ux(:,:,k),uy(:,:,k),T_dot,dx,dx,dt,solver)
        end do 

        call fill_borders_3D(T_ice,nfill=2)
        
        return 

    end subroutine calc_enth_horizontal_advection_3D

    subroutine calc_adv2D_expl_rate(dvardt,var,ux,uy,dx,dy,i,j)
        ! Solve 2D advection equation for ice sheet thickness via explicit flux divergence:
        ! d[H]/dt = -grad[H*(ux,uy)] 
        !
        ! ajr: adapted from IMAU-ICE code from Heiko Goelzer (h.goelzer@uu.nl) 2016
        ! Note: original algorithm called for interpolation of var to ab-nodes explicitly. 
        ! It also works using the ac-nodes directly, but is less stable. This can be chosen via the parameter
        ! use_ab_expl. 

        implicit none 

        real(prec), intent(OUT)   :: dvardt                 ! [m/yr] aa-nodes, Ice thickness 
        real(prec), intent(IN)    :: var(:,:)               ! [m] aa-nodes, Ice thickness 
        real(prec), intent(IN)    :: ux(:,:)                ! [m a^-1] ac-nodes, Horizontal velocity, x-direction
        real(prec), intent(IN)    :: uy(:,:)                ! [m a^-1] ac-nodes, Horizontal velocity, y-direction
        real(prec), intent(IN)    :: dx                     ! [m] Horizontal grid spacing, x-direction
        real(prec), intent(IN)    :: dy                     ! [m] Horizontal grid spacing, y-direction
        integer,    intent(IN)    :: i, j 

        ! Local variables:
        integer                 :: nx, ny 

        real(prec) :: flux_xr                  ! [m^2 a^-1] ac-nodes, Flux in the x-direction to the right
        real(prec) :: flux_xl                  ! [m^2 a^-1] ac-nodes, Flux in the x-direction to the left
        real(prec) :: flux_yu                  ! [m^2 a^-1] ac-nodes, Flux in the y-direction upwards
        real(prec) :: flux_yd                  ! [m^2 a^-1] ac-nodes, Flux in the y-direction downwards

        ! Calculate the flux across each boundary [m^2 a^-1]
        flux_xr = ux(i  ,j  ) * 0.5 * (var(i  ,j  ) + var(i+1,j  ))
        flux_xl = ux(i-1,j  ) * 0.5 * (var(i-1,j  ) + var(i  ,j  ))
        flux_yu = uy(i  ,j  ) * 0.5 * (var(i  ,j  ) + var(i  ,j+1))
        flux_yd = uy(i  ,j-1) * 0.5 * (var(i  ,j-1) + var(i  ,j  ))

        ! Calculate flux divergence on aa-node 
        dvardt = (1.0 / dx) * (flux_xl - flux_xr) + (1.0 / dy) * (flux_yd - flux_yu) 

        return 

    end subroutine calc_adv2D_expl_rate

    subroutine  calc_adv2D_impl_upwind_rate(dHdt,H_prev,ux,uy,mdot,dx,dy,dt,f_upwind)
        ! To solve the 2D adevection equation:
        ! dH/dt =
        ! M H = Frelax
        ! ajr: adapted from GRISLI (Ritz et al., 1997)

        implicit none

        real(prec), intent(OUT)   :: dHdt(:,:)      ! Variable rate of change (aa-node)
        real(prec), intent(IN)    :: H_prev(:,:)    ! Variable (aa-node)
        real(prec), intent(IN)    :: ux(:,:)        ! Depth-averaged velocity - x direction (ac-node)
        real(prec), intent(IN)    :: uy(:,:)        ! Depth-averaged velocity - y direction (ac-node)
        real(prec), intent(IN)    :: mdot(:,:)      ! Total column mass balance (aa-node)
        real(prec), intent(IN)    :: dx             ! [m] x-resolution
        real(prec), intent(IN)    :: dy             ! [m] y-resolution
        real(prec), intent(IN)    :: dt             ! [a] Timestep (assumes dx=dy)
        real(prec), intent(IN)    :: f_upwind       ! [-] Fraction of "upwind-ness" to apply (ajr: experimental!) - between 0.5 and 1.0, default f_upwind=1.0
        ! Local variables
        integer    :: i, j, nx, ny
        integer    :: iter, ierr 
        real(prec) :: dtdx, dtdx2
        real(prec) :: reste, delh, tmp 
        real(prec), allocatable :: crelax(:,:)      ! diagnonale de M
        real(prec), allocatable :: arelax(:,:)      ! sous diagonale selon x
        real(prec), allocatable :: brelax(:,:)      ! sur  diagonale selon x
        real(prec), allocatable :: drelax(:,:)      ! sous diagonale selon y
        real(prec), allocatable :: erelax(:,:)      ! sur  diagonale selon y
        real(prec), allocatable :: frelax(:,:)      ! vecteur
        real(prec), allocatable :: c_west(:,:)      ! sur demi mailles Ux
        real(prec), allocatable :: c_east(:,:)      ! sur demi mailles Ux
        real(prec), allocatable :: c_north(:,:)     ! sur demi mailles Uy
        real(prec), allocatable :: c_south(:,:)     ! sur demi mailles Uy
        real(prec), allocatable :: deltaH(:,:)      ! Change in H

        real(prec), allocatable :: H(:,:)
        
        logical,    parameter :: use_upwind = .TRUE.  ! Apply upwind advection scheme or central scheme?   
                                                      ! (now this is redundant with f_upwind parameter) 
        
        ! Note: f_upwind=0.6 gives good agreement with EISMINT1 summit elevation
        ! for the moving and fixed margin experiments, when using the calc_shear_3D approach.
        ! (f_upwind=0.5, ie central method, works well when using the velocity_sia approach)

        ! Note: it may be that f_upwind>0.5 is more stable for real dynamic simulations

        ! Determine array size
        nx = size(H,1)
        ny = size(H,2)

        ! Allocate local arrays
        allocate(crelax(nx,ny))
        allocate(arelax(nx,ny))
        allocate(brelax(nx,ny))
        allocate(drelax(nx,ny))
        allocate(erelax(nx,ny))
        allocate(frelax(nx,ny))
        allocate(c_west(nx,ny))
        allocate(c_east(nx,ny))
        allocate(c_north(nx,ny))
        allocate(c_south(nx,ny))

        allocate(deltaH(nx,ny))

        allocate(H(nx,ny))
        
        ! Define some helpful values
        dtdx2 = dt/(dx**2)
        dtdx  = dt/dx

        ! Initialize relaxation arrays
        arelax = 0.0
        brelax = 0.0
        drelax = 0.0
        erelax = 0.0
        crelax = 1.0
        frelax = 0.0
        deltaH = 0.0

        H = H_prev  

        ! Modify coefficients depending on method (upwind, central)
        if (use_upwind) then
            ! Upwind method

            if (f_upwind .eq. 1.0) then
                ! Apply normal upwind scheme

                where (ux.ge.0.0)
                    c_west = 1.0
                    c_east = 0.0
                elsewhere
                    c_west = 0.0
                    c_east = 1.0
                end where

                where (uy.ge.0.0)
                    c_south = 1.0
                    c_north = 0.0
                elsewhere
                    c_south = 0.0
                    c_north = 1.0
                end where

            else 
                ! Apply fractional upwind scheme to reduce numerical diffusion,
                ! but benefit from upwind stability (ajr: experimental!)
                ! f_upwind = 0.5 => central difference, f_upwind = 1.0 => full upwind 

                where (ux.ge.0.0)
                    c_west = f_upwind
                    c_east = 1.0 - f_upwind
                elsewhere
                    c_west = 1.0 - f_upwind
                    c_east = f_upwind
                end where

                where (uy.ge.0.0)
                    c_south = f_upwind
                    c_north = 1.0 - f_upwind
                elsewhere
                    c_south = 1.0 - f_upwind
                    c_north = f_upwind
                end where

            end if 

        else
            ! Central method

            c_west  = 0.5
            c_east  = 0.5
            c_south = 0.5
            c_north = 0.5

        end if

        ! Populate diagonals
        do j=2,ny-1
        do i=2,nx-1

            !  sous diagonale en x
            arelax(i,j) = -dtdx*c_west(i-1,j)*ux(i-1,j)    ! partie advective en x

            !  sur diagonale en x
            brelax(i,j) = +dtdx*c_east(i,j)*ux(i,j)        ! partie advective

            !  sous diagonale en y
            drelax(i,j) = -dtdx*c_south(i,j-1)*uy(i,j-1)   ! partie advective en y

            !  sur diagonale en y
            erelax(i,j) = +dtdx*c_north(i,j)*uy(i,j)       ! partie advective


            ! diagonale
            crelax(i,j) = 1.0 + dtdx* &
                       ((c_west(i,j)*ux(i,j) - c_east(i-1,j)*ux(i-1,j)) &
                      +(c_south(i,j)*uy(i,j) - c_north(i,j-1)*uy(i,j-1)))

            ! Combine all terms
            frelax(i,j) = H(i,j) + dt*mdot(i,j)

        end do
        end do

        ! Avoid underflows 
        where (abs(arelax) .lt. tol_underflow) arelax = 0.0_prec 
        where (abs(brelax) .lt. tol_underflow) brelax = 0.0_prec 
        where (abs(drelax) .lt. tol_underflow) drelax = 0.0_prec 
        where (abs(erelax) .lt. tol_underflow) erelax = 0.0_prec 
        
        ! Initialize new H solution to zero (to get zeros at boundaries)
        H  = 0.0

        ! Initially assume convergence criterion is not satisfied 
        ierr = -1   ! convergence criterion not fulfilled
        
        do iter = 1, 1000 
            ! Relaxation loop 

            ! Calculate change in H
            do j=2,ny-1
            do i=2,nx-1

                reste = (((arelax(i,j)*H(i-1,j) + drelax(i,j)*H(i,j-1)) &
                        + (brelax(i,j)*H(i+1,j) + erelax(i,j)*H(i,j+1))) &
                        +  crelax(i,j)*H(i,j))  - frelax(i,j)

                deltaH(i,j) = reste/crelax(i,j)

            end do
            end do

            ! Adjust H to new value
            H = H - deltaH

            ! Check stopping criterion (something like rmse of remaining change in H)
            where(abs(deltaH) .lt. tol_underflow) deltaH = 0.0_prec      ! Avoid underflows
            delh = sqrt(sum(deltaH**2)) / ((nx-2)*(ny-2))
            
            ! Use simple stopping criterion: maximum remaining change in H
            ! Note: this is less likely to converge given the same stopping
            ! criterion.
!             delh = maxval(abs(deltaH))
            
            if ( delh .lt. 1e-4) then
                ! Solution has converged, exit  
                ierr = 0 
                exit 
            end if 

        end do ! End of relaxation loop

        if (ierr .eq. -1) then 
            write(*,*) "calc_adv2D_impl_upwind_rate:: Error: advection did not converge."
            write(*,"(a,i10,g12.3)") "iter, delh = ", iter, delh 
            stop 
        end if 

        ! Calculate rate of change 
        dHdt = (H - H_prev) / dt 

        return

    end subroutine calc_adv2D_impl_upwind_rate

    subroutine ytherm_par_load(par,filename,zeta_aa,zeta_ac,nx,ny,dx,init)

        type(ytherm_param_class), intent(OUT) :: par
        character(len=*),         intent(IN)  :: filename
        real(prec),               intent(IN)  :: zeta_aa(:)  
        real(prec),               intent(IN)  :: zeta_ac(:)  
        integer,                  intent(IN)  :: nx, ny 
        real(prec),               intent(IN)  :: dx 
        logical, optional,        intent(IN)  :: init

        ! Local variables 
        logical :: init_pars 
        integer :: k 

        init_pars = .FALSE.
        if (present(init)) init_pars = .TRUE. 
 
        ! Store local parameter values in output object
        call nml_read(filename,"ytherm","method",         par%method,           init=init_pars)
        call nml_read(filename,"ytherm","solver_advec",   par%solver_advec,     init=init_pars)
        call nml_read(filename,"ytherm","gamma",          par%gamma,            init=init_pars)
        call nml_read(filename,"ytherm","use_strain_sia", par%use_strain_sia,   init=init_pars)
        call nml_read(filename,"ytherm","n_sm_qstrn",     par%n_sm_qstrn,       init=init_pars)
        call nml_read(filename,"ytherm","n_sm_qb",        par%n_sm_qb,          init=init_pars)
        call nml_read(filename,"ytherm","use_const_cp",   par%use_const_cp,     init=init_pars)
        call nml_read(filename,"ytherm","const_cp",       par%const_cp,         init=init_pars)
        call nml_read(filename,"ytherm","use_const_kt",   par%use_const_kt,     init=init_pars)
        call nml_read(filename,"ytherm","const_kt",       par%const_kt,         init=init_pars)
        call nml_read(filename,"ytherm","enth_cr",        par%enth_cr,          init=init_pars)
        call nml_read(filename,"ytherm","omega_max",      par%omega_max,        init=init_pars)
        call nml_read(filename,"ytherm","till_rate",      par%till_rate,        init=init_pars)
        call nml_read(filename,"ytherm","H_w_max",        par%H_w_max,          init=init_pars)
        
        ! Set internal parameters
        par%nx  = nx
        par%ny  = ny 
        par%dx  = dx
        par%dy  = dx  
        par%nz_aa = size(zeta_aa,1)     ! bottom, layer centers, top 
        par%nz_ac = size(zeta_ac,1)     ! layer borders

        if (allocated(par%zeta_aa)) deallocate(par%zeta_aa)
        allocate(par%zeta_aa(par%nz_aa))
        par%zeta_aa = zeta_aa 
        
        if (allocated(par%zeta_ac)) deallocate(par%zeta_ac)
        allocate(par%zeta_ac(par%nz_ac))
        par%zeta_ac = zeta_ac 
        
        if (allocated(par%dzeta_a)) deallocate(par%dzeta_a)
        allocate(par%dzeta_a(par%nz_aa))
        if (allocated(par%dzeta_b)) deallocate(par%dzeta_b)
        allocate(par%dzeta_b(par%nz_aa))
        
        ! Define thermodynamic zeta helper variables dzeta_a/dzeta_b
        call calc_dzeta_terms(par%dzeta_a,par%dzeta_b,par%zeta_aa,par%zeta_ac)

        ! Define current time as unrealistic value
        par%time = 1000000000   ! [a] 1 billion years in the future

        return

    end subroutine ytherm_par_load

    subroutine ytherm_alloc(now,nx,ny,nz_aa,nz_ac,nzr)

        implicit none 

        type(ytherm_state_class), intent(INOUT) :: now 
        integer, intent(IN) :: nx, ny, nz_aa, nz_ac, nzr 

        call ytherm_dealloc(now)

        allocate(now%enth(nx,ny,nz_aa))
        allocate(now%T_ice(nx,ny,nz_aa))
        allocate(now%omega(nx,ny,nz_aa))
        allocate(now%T_pmp(nx,ny,nz_aa))
        allocate(now%bmb_grnd(nx,ny))
        allocate(now%f_pmp(nx,ny))
        allocate(now%Q_strn(nx,ny,nz_aa))
        allocate(now%Q_b(nx,ny))
        allocate(now%Q_ice_b(nx,ny))
        allocate(now%cp(nx,ny,nz_aa))
        allocate(now%kt(nx,ny,nz_aa))
        allocate(now%H_cts(nx,ny))
        allocate(now%T_prime_b(nx,ny))
        allocate(now%H_w(nx,ny))
        allocate(now%dHwdt(nx,ny))

        now%enth      = 0.0
        now%T_ice     = 0.0
        now%omega     = 0.0  
        now%T_pmp     = 0.0
        now%bmb_grnd  = 0.0 
        now%f_pmp     = 0.0 
        now%Q_strn    = 0.0 
        now%Q_b       = 0.0 
        now%Q_ice_b   = 0.0 
        now%cp        = 0.0 
        now%kt        = 0.0 
        now%H_cts     = 0.0 
        now%T_prime_b = 0.0 
        now%H_w       = 0.0 
        now%dHwdt     = 0.0 

        return

    end subroutine ytherm_alloc 

    subroutine ytherm_dealloc(now)

        implicit none 

        type(ytherm_state_class), intent(INOUT) :: now

        if (allocated(now%enth))      deallocate(now%enth)
        if (allocated(now%T_ice))     deallocate(now%T_ice)
        if (allocated(now%omega))     deallocate(now%omega)
        if (allocated(now%T_pmp))     deallocate(now%T_pmp)
        if (allocated(now%bmb_grnd))  deallocate(now%bmb_grnd)
        if (allocated(now%f_pmp))     deallocate(now%f_pmp)
        if (allocated(now%Q_strn))    deallocate(now%Q_strn)
        if (allocated(now%Q_b))       deallocate(now%Q_b)
        if (allocated(now%Q_ice_b))   deallocate(now%Q_ice_b)
        if (allocated(now%cp))        deallocate(now%cp)
        if (allocated(now%kt))        deallocate(now%kt)
        if (allocated(now%H_cts))     deallocate(now%H_cts)
        if (allocated(now%T_prime_b)) deallocate(now%T_prime_b)
        if (allocated(now%H_w))       deallocate(now%H_w)
        if (allocated(now%dHwdt))     deallocate(now%dHwdt)
        
        return 

    end subroutine ytherm_dealloc 

end module yelmo_thermodynamics
