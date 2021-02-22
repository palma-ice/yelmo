
module yelmo_thermodynamics

    use nml 
    use yelmo_defs 
    use yelmo_grid, only : calc_zeta
    use yelmo_tools, only : smooth_gauss_2D, smooth_gauss_3D, gauss_values, fill_borders_2D, fill_borders_3D, &
            stagger_aa_ab
    
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
        integer :: k0, nz_aa    
        real(prec) :: dt 
        real(prec), allocatable :: H_w_now(:,:)
        real(prec), allocatable :: dTdz_b_now(:,:)
        
        nx    = thrm%par%nx
        ny    = thrm%par%ny
        nz_aa = thrm%par%ztot%nz_aa 

        ! Assign shortcut to ice basal layer index 
        k0 = thrm%par%k0 

        allocate(H_w_now(nx,ny)) 
        H_w_now = 0.0_prec 

        ! Initialize time if necessary 
        if (thrm%par%time .gt. dble(time)) then 
            thrm%par%time = dble(time)
        end if 

        ! Get time step and advance current time 
        dt            = dble(time) - thrm%par%time 
        thrm%par%time = dble(time) 

        ! === Determine some thermal properties === 

        ! == ice column == 

        ! Calculate the specific heat capacity of the ice
        if (thrm%par%use_const_cp) then 
            thrm%now%cp(:,:,k0:nz_aa)  = thrm%par%const_cp
        else  
            thrm%now%cp(:,:,k0:nz_aa)  = calc_specific_heat_capacity(thrm%now%T_tot(:,:,k0:nz_aa))
        end if 
        
        ! Calculate the heat conductivity of the ice
        if (thrm%par%use_const_kt) then 
            thrm%now%kt(:,:,k0:nz_aa)  = thrm%par%const_kt
        else  
            thrm%now%kt(:,:,k0:nz_aa)  = calc_thermal_conductivity(thrm%now%T_tot(:,:,k0:nz_aa))
        end if 

        ! Calculate the pressure-corrected melting point (in Kelvin)
        do k = k0, nz_aa  
            thrm%now%T_pmp_tot(:,:,k) = calc_T_pmp(tpo%now%H_ice,thrm%par%ztot%zeta_aa(k),T0,T_pmp_beta)
        end do 

        ! == lithosphere column == 

        ! Prescibe constant values of heat capacity, thermal conductivity 
        ! and pressure melting point withint the lithospheric column 

        thrm%now%cp(:,:,1:k0-1) = thrm%par%lith_cp 
        thrm%now%kt(:,:,1:k0-1) = thrm%par%lith_kt 
        thrm%now%T_pmp_tot(:,:,1:k0-1) = 1e8      ! [K] Set to unreachably high value
        
        ! === Calculate heat source terms (Yelmo vertical grid) === 

select case("ab")

    case("aa")
        ! Calculate the basal frictional heating (from ab-nodes)
        call calc_basal_heating_fromaa(thrm%now%Q_b,dyn%now%ux_b,dyn%now%uy_b,dyn%now%taub_acx,dyn%now%taub_acy, &
                                            beta1=thrm%par%dt_beta(1),beta2=thrm%par%dt_beta(2))

    case("ab")
        ! Calculate the basal frictional heating (from ab-nodes)
        call calc_basal_heating_fromab(thrm%now%Q_b,dyn%now%ux_b,dyn%now%uy_b,dyn%now%taub_acx,dyn%now%taub_acy, &
                                            tpo%now%H_ice,beta1=thrm%par%dt_beta(1),beta2=thrm%par%dt_beta(2))

    case("ac")
        ! ajr: old interface with scaling optional via f_pmp
        call calc_basal_heating_fromac(thrm%now%Q_b,dyn%now%ux_b,dyn%now%uy_b,dyn%now%taub_acx,dyn%now%taub_acy, &
                        tpo%now%H_ice,thrm%now%f_pmp,beta1=thrm%par%dt_beta(1),beta2=thrm%par%dt_beta(2))

    case("interp")

        call calc_basal_heating_interp(thrm%now%Q_b,dyn%now%ux_b,dyn%now%uy_b,dyn%now%taub_acx,dyn%now%taub_acy, &
                                            tpo%now%H_ice,beta1=thrm%par%dt_beta(1),beta2=thrm%par%dt_beta(2))

end select
        
        ! Smooth basal frictional heating 
        if (thrm%par%n_sm_qb .gt. 0) then 
            call smooth_gauss_2D(thrm%now%Q_b,tpo%now%H_ice.gt.0.0,thrm%par%dx,thrm%par%n_sm_qb, &
                                    tpo%now%H_ice.gt.0.0)
        end if 
        

        return

    end subroutine calc_ytherm

    subroutine calc_ytherm_enthalpy_3D(enth,T_ice,omega,bmb_grnd,Q_ice_b,H_cts,T_pmp,cp,kt,advecxy,ux,uy,uz,Q_strn,Q_b,Q_lith, &
                                        T_srf,H_ice,z_srf,H_w,dHwdt,H_grnd,f_grnd,dHdt,dzsdt,zeta_aa,zeta_ac,dzeta_a,dzeta_b, &
                                        cr,omega_max,dt,dx,solver,solver_advec)
        ! This wrapper subroutine breaks the thermodynamics problem into individual columns,
        ! which are solved independently by calling calc_enth_column

        ! Note zeta=height, k=1 base, k=nz surface 
        
        !$ use omp_lib

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
        real(prec), intent(IN)    :: advecxy(:,:,:) ! [m a-1] Horizontal x-velocity 
        real(prec), intent(IN)    :: ux(:,:,:)      ! [m a-1] Horizontal x-velocity 
        real(prec), intent(IN)    :: uy(:,:,:)      ! [m a-1] Horizontal y-velocity 
        real(prec), intent(IN)    :: uz(:,:,:)      ! [m a-1] Vertical velocity 
        real(prec), intent(IN)    :: Q_strn(:,:,:)  ! [K a-1] Internal strain heat production in ice
        real(prec), intent(IN)    :: Q_b(:,:)       ! [J a-1 m-2] Basal frictional heat production 
        real(prec), intent(IN)    :: Q_lith(:,:)    ! [mW m-2] Heat flux at bed surface from lithosphere (like Q_geo)
        real(prec), intent(IN)    :: T_srf(:,:)     ! [K] Surface temperature 
        real(prec), intent(IN)    :: H_ice(:,:)     ! [m] Ice thickness 
        real(prec), intent(IN)    :: z_srf(:,:)     ! [m] Surface elevation 
        real(prec), intent(IN)    :: H_w(:,:)       ! [m] Basal water layer thickness 
        real(prec), intent(IN)    :: dHwdt(:,:)     ! [m/a] Basal water layer thickness change
        real(prec), intent(IN)    :: H_grnd(:,:)    ! [--] Ice thickness above flotation 
        real(prec), intent(IN)    :: f_grnd(:,:)    ! [--] Grounded fraction
        real(prec), intent(IN)    :: dHdt(:,:)      ! [m/a] Ice thickness change
        real(prec), intent(IN)    :: dzsdt(:,:)     ! [m/a] Ice surface change
        real(prec), intent(IN)    :: zeta_aa(:)     ! [--] Vertical sigma coordinates (zeta==height), aa-nodes
        real(prec), intent(IN)    :: zeta_ac(:)     ! [--] Vertical sigma coordinates (zeta==height), ac-nodes
        real(prec), intent(IN)    :: dzeta_a(:)     ! nz_aa [--] Solver discretization helper variable ak
        real(prec), intent(IN)    :: dzeta_b(:)     ! nz_aa [--] Solver discretization helper variable bk
        real(prec), intent(IN)    :: cr             ! [--] Conductivity ratio for temperate ice (kappa_temp = enth_cr*kappa_cold)
        real(prec), intent(IN)    :: omega_max      ! [--] Maximum allowed water content fraction 
        real(prec), intent(IN)    :: dt             ! [a] Time step 
        real(prec), intent(IN)    :: dx             ! [a] Horizontal grid step 
        character(len=*), intent(IN) :: solver      ! "enth" or "temp" 
        character(len=*), intent(IN) :: solver_advec    ! "expl" or "impl-upwind"

        ! Local variables
        integer :: i, j, k, nx, ny, nz_aa, nz_ac  
        real(prec) :: T_shlf, H_grnd_lim, f_scalar, T_base  
        real(prec) :: H_ice_now 

        real(prec) :: filter0(3,3), filter(3,3) 

        real(prec), parameter :: H_ice_thin = 10.0   ! [m] Threshold to define 'thin' ice

        nx    = size(T_ice,1)
        ny    = size(T_ice,2)
        nz_aa = size(zeta_aa,1)
        nz_ac = size(zeta_ac,1)

        ! Initialize gaussian filter kernel for smoothing ice thickness at the margin
        filter0 = gauss_values(dx,dx,sigma=2.0*dx,n=size(filter,1))

        ! ===================================================

        ! ajr: openmp problematic here - leads to NaNs
        !!!$omp parallel do
        do j = 2, ny-1
        do i = 2, nx-1 
            
            H_ice_now = H_ice(i,j) 
 
            ! Filter the ice thickness at the margin only to avoid a large
            ! gradient in ice thickness there to rather thin ice points - 
            ! helps with stability, particularly for EISMINT2 experiments.
            if (H_ice(i,j) .gt. 0.0 .and. count(H_ice(i-1:i+1,j-1:j+1) .eq. 0.0) .ge. 2) then
                filter = filter0 
                where (H_ice(i-1:i+1,j-1:j+1) .eq. 0.0) filter = 0.0 
                filter = filter/sum(filter)
                H_ice_now = sum(H_ice(i-1:i+1,j-1:j+1)*filter)
            end if

            ! For floating points, calculate the approximate marine-shelf temperature 
            ! ajr, later this should come from an external model, and T_shlf would
            ! be the boundary variable directly
            if (f_grnd(i,j) .lt. 1.0) then 

                ! Calculate approximate marine freezing temp, limited to pressure melting point 
                T_shlf = calc_T_base_shlf_approx(H_ice_now,T_pmp(i,j,1),H_grnd(i,j))

            else 
                ! Assigned for safety 

                T_shlf   = T_pmp(i,j,1)

            end if 

            if (H_ice_now .le. H_ice_thin) then 
                ! Ice is too thin or zero: prescribe linear temperature profile
                ! between temperate ice at base and surface temperature 
                ! (accounting for floating/grounded nature via T_base)

                if (f_grnd(i,j) .lt. 1.0) then 
                    ! Impose T_shlf for the basal temperature
                    T_base = T_shlf 
                else
                    ! Impose temperature at the pressure melting point of grounded ice 
                    T_base = T_pmp(i,j,1) 
                end if 

                T_ice(i,j,:)  = calc_temp_linear_column(T_srf(i,j),T_base,T_pmp(i,j,nz_aa),zeta_aa)
                omega(i,j,:)  = 0.0_prec 
                call convert_to_enthalpy(enth(i,j,:),T_ice(i,j,:),omega(i,j,:),T_pmp(i,j,:),cp(i,j,:),L_ice)
                bmb_grnd(i,j) = 0.0_prec
                Q_ice_b(i,j)  = 0.0_prec 
                H_cts(i,j)    = 0.0_prec

            else 
                ! Thick ice exists, call thermodynamic solver for the column

                if (trim(solver) .eq. "enth") then 

                    call calc_enth_column(enth(i,j,:),T_ice(i,j,:),omega(i,j,:),bmb_grnd(i,j),Q_ice_b(i,j),H_cts(i,j), &
                            T_pmp(i,j,:),cp(i,j,:),kt(i,j,:),advecxy(i,j,:),uz(i,j,:),Q_strn(i,j,:),Q_b(i,j),Q_lith(i,j),T_srf(i,j), &
                            T_shlf,H_ice_now,H_w(i,j),f_grnd(i,j),zeta_aa,zeta_ac,dzeta_a,dzeta_b,cr,omega_max,T0,dt)
                
                else 

                    call calc_temp_column(enth(i,j,:),T_ice(i,j,:),omega(i,j,:),bmb_grnd(i,j),Q_ice_b(i,j),H_cts(i,j), &
                    T_pmp(i,j,:),cp(i,j,:),kt(i,j,:),advecxy(i,j,:),uz(i,j,:),Q_strn(i,j,:),Q_b(i,j),Q_lith(i,j),T_srf(i,j), &
                    T_shlf,H_ice_now,H_w(i,j),f_grnd(i,j),zeta_aa,zeta_ac,dzeta_a,dzeta_b,omega_max,T0,dt)                
                
                end if 

            end if 

        end do 
        end do 
        !!!$omp end parallel do

        ! Fill in borders 
        call fill_borders_3D(enth,nfill=1)
        call fill_borders_3D(T_ice,nfill=1)
        call fill_borders_3D(omega,nfill=1)
        call fill_borders_2D(bmb_grnd,nfill=1)
        
        return 

    end subroutine calc_ytherm_enthalpy_3D

    subroutine calc_ytherm_enthalpy_bedrock_3D(enth,T_lith,T_ice_b,T_pmp_b,cp,kt,H_lith, &
                                                H_ice,H_grnd,Q_geo,zeta_aa,zeta_ac,dzeta_a,dzeta_b,dt)
        ! This wrapper subroutine breaks the thermodynamics problem into individual columns,
        ! which are solved independently by calling calc_enth_column

        ! Note zeta=height, k=1 base, k=nz surface 
        
        !$ use omp_lib

        implicit none 

        real(prec), intent(INOUT) :: enth(:,:,:)    ! [J m-3] Lithosphere enthalpy
        real(prec), intent(INOUT) :: T_lith(:,:,:)  ! [K] Lithosphere temperature
        real(prec), intent(IN)    :: T_ice_b(:,:)   ! [K] Ice temperature at ice base
        real(prec), intent(IN)    :: T_pmp_b(:,:)   ! [K] Pressure melting point temp at ice base.
        real(prec), intent(IN)    :: cp(:,:,:)      ! [J kg-1 K-1] Specific heat capacity lithosphere
        real(prec), intent(IN)    :: kt(:,:,:)      ! [J a-1 m-1 K-1] Heat conductivity lithosphere
        real(prec), intent(IN)    :: H_lith(:,:)    ! [m] Lithosphere thickness 
        real(prec), intent(IN)    :: H_ice(:,:)     ! [m] Ice thickness 
        real(prec), intent(IN)    :: H_grnd(:,:)    ! [--] Ice thickness above flotation 
        real(prec), intent(IN)    :: Q_geo(:,:)     ! [mW m-2] Geothermal heat flux deep in bedrock
        real(prec), intent(IN)    :: zeta_aa(:)     ! [--] Vertical sigma coordinates (zeta==height), aa-nodes
        real(prec), intent(IN)    :: zeta_ac(:)     ! [--] Vertical sigma coordinates (zeta==height), ac-nodes
        real(prec), intent(IN)    :: dzeta_a(:)     ! nz_aa [--] Solver discretization helper variable ak
        real(prec), intent(IN)    :: dzeta_b(:)     ! nz_aa [--] Solver discretization helper variable bk
        real(prec), intent(IN)    :: dt             ! [a] Time step 

        ! Local variables
        integer :: i, j, k, nx, ny, nz_aa, nz_ac  
        real(prec) :: T_base

        nx    = size(T_lith,1)
        ny    = size(T_lith,2)
        nz_aa = size(zeta_aa,1)
        nz_ac = size(zeta_ac,1)

        ! ===================================================

        ! ajr: openmp problematic here - leads to NaNs
        !!!$omp parallel do
        do j = 1, ny
        do i = 1, nx 

            ! For floating points, calculate the approximate marine-shelf temperature 
            ! although really the temperature at the bottom of the ocean is needed
            if (H_grnd(i,j) .lt. 0.0) then 

                ! Calculate approximate marine freezing temp, limited to pressure melting point 
                T_base = calc_T_base_shlf_approx(H_ice(i,j),T_pmp_b(i,j),H_grnd(i,j))

            else 
                ! Assign ice basal temperature
                T_base   = T_ice_b(i,j)

            end if 

            if (H_ice(i,j) .gt. 0.0) then 
                ! Call thermodynamic solver for the column

                call calc_temp_column_bedrock(enth(i,j,:),T_lith(i,j,:),cp(i,j,:),kt(i,j,:),  &
                        Q_geo(i,j),T_base,H_lith(i,j),zeta_aa,zeta_ac,dzeta_a,dzeta_b,dt)
            
            else 
                ! Assume equilibrium conditions: impose linear temperature 
                ! profile following Q_geo and T_base

                T_lith(i,j,:) = calc_temp_lith_column(zeta_aa,kt(i,j,:),cp(i,j,:), &
                                                    rho_l,H_lith(i,j),T_base,Q_geo(i,j))

                ! Get enthalpy too 
                call convert_to_enthalpy(enth(i,j,:),T_lith(i,j,:),0.0_wp,0.0_wp,cp(i,j,:),0.0_wp)

            end if 

        end do 
        end do 
        !!!$omp end parallel do

        return 

    end subroutine calc_ytherm_enthalpy_bedrock_3D
    
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
        real(prec), allocatable :: dTdt(:,:) 

        nx = size(T_ice,1)
        ny = size(T_ice,2)
        nz = size(T_ice,3) 

        allocate(dTdt(nx,ny)) 
        dTdt = 0.0_prec 

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
            call calc_advec2D(dTdt,T_ice(:,:,k),ux(:,:,k),uy(:,:,k),(T_ice(:,:,k)*0.0_prec),dx,dx,dt,solver)
            T_ice(:,:,k) = T_ice(:,:,k) + dt*dTdt 
        end do 

        call fill_borders_3D(T_ice,nfill=2)
        
        return 

    end subroutine calc_enth_horizontal_advection_3D

    subroutine ytherm_par_load(par,filename,nx,ny,dx,init)

        type(ytherm_param_class), intent(OUT) :: par
        character(len=*),         intent(IN)  :: filename
        integer,                  intent(IN)  :: nx, ny 
        real(prec),               intent(IN)  :: dx 
        logical, optional,        intent(IN)  :: init

        ! Local variables 
        logical :: init_pars 
        integer :: k 

        init_pars = .FALSE.
        if (present(init)) init_pars = .TRUE. 
 
        ! Store local parameter values in output object
        call nml_read(filename,"ytherm","nz_aa",          par%zice%nz_aa,       init=init_pars)
        call nml_read(filename,"ytherm","zeta_scale",     par%zice%zeta_scale,  init=init_pars)
        call nml_read(filename,"ytherm","zeta_exp",       par%zice%zeta_exp,    init=init_pars)
        
        call nml_read(filename,"ytherm","method",         par%method,           init=init_pars)
        call nml_read(filename,"ytherm","method_lith",    par%lith_method,      init=init_pars)
        call nml_read(filename,"ytherm","dt_method",      par%dt_method,        init=init_pars)
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
        
        ! Lithosphere specific settings
        call nml_read(filename,"ytherm_lith","nz_aa",     par%zlith%nz_aa,      init=init_pars)
        call nml_read(filename,"ytherm_lith","zeta_scale",par%zlith%zeta_scale, init=init_pars)
        call nml_read(filename,"ytherm_lith","zeta_exp",  par%zlith%zeta_exp,   init=init_pars)
        
        call nml_read(filename,"ytherm_lith","H_lith",    par%H_lith,           init=init_pars)
        call nml_read(filename,"ytherm_lith","cp",        par%lith_cp,          init=init_pars)
        call nml_read(filename,"ytherm_lith","kt",        par%lith_kt,          init=init_pars)
        

        ! In case of method=="temp", prescribe some parameters
        if (trim(par%method) .eq. "temp") then  
            par%enth_cr   = 1.0_prec 
            par%omega_max = 0.0_prec 
        end if 

        ! Set internal parameters
        par%nx  = nx
        par%ny  = ny 
        par%dx  = dx
        par%dy  = dx  
        
        ! Define different vertical axes =====

        ! Set ac-node lengths
        par%zice%nz_ac  = par%zice%nz_aa - 1 
        par%zlith%nz_ac = par%zlith%nz_aa - 1 
        
        ! == Ice column == 
        call calc_zeta(par%zice%zeta_aa,par%zice%zeta_ac,par%zice%nz_aa,&
                            par%zice%zeta_scale,par%zice%zeta_exp)

        ! == Lithosphere column == (should go from -1 to 0)
        call calc_zeta(par%zlith%zeta_aa,par%zlith%zeta_ac,par%zlith%nz_aa,&
                            par%zlith%zeta_scale,par%zlith%zeta_exp)
        par%zlith%zeta_aa = par%zlith%zeta_aa - 1.0_wp 
        par%zlith%zeta_ac = par%zlith%zeta_ac - 1.0_wp 
        
        ! == Total column == 
        ! Lithosphere from -1 to 0, ice column from 0 to 1, 
        ! zlith%zeta_aa(par%zlith%nz_aa) == zice%zeta_aa(1) = 0 

        ! Get index of ice base 
        par%k0 = par%zlith%nz_aa 

        par%ztot%nz_aa = par%zice%nz_aa + par%zlith%nz_aa - 1 
        par%ztot%nz_ac = par%ztot%nz_aa - 1 

        if (allocated(par%ztot%zeta_aa)) deallocate(par%ztot%zeta_aa)
        if (allocated(par%ztot%zeta_ac)) deallocate(par%ztot%zeta_ac)
        allocate(par%ztot%zeta_aa(par%ztot%nz_aa))
        allocate(par%ztot%zeta_ac(par%ztot%nz_ac))
        
        par%ztot%zeta_aa(1:par%k0)              = par%zlith%zeta_aa 
        par%ztot%zeta_aa(par%k0:par%ztot%nz_aa) = par%zice%zeta_aa 

        par%ztot%zeta_ac(1:par%k0-1)            = par%zlith%zeta_ac 
        par%ztot%zeta_ac(par%k0:par%ztot%nz_ac) = par%zice%zeta_ac 

        par%ztot%zeta_ac(par%k0-1:par%k0) = [-1e-5,1e-5]

        if (allocated(par%ztot%dzeta_a)) deallocate(par%ztot%dzeta_a)
        if (allocated(par%ztot%dzeta_b)) deallocate(par%ztot%dzeta_b)
        allocate(par%ztot%dzeta_a(par%ztot%nz_aa))
        allocate(par%ztot%dzeta_b(par%ztot%nz_aa))
        
        par%ztot%dzeta_a = 0.0_wp 
        par%ztot%dzeta_b = 0.0_wp 
        
        call calc_dzeta_terms(par%ztot%dzeta_a,par%ztot%dzeta_b,par%ztot%zeta_aa,par%ztot%zeta_ac)

        ! write(*,*) "========="
        ! do k = par%ztot%nz_aa, 1, -1
        !     if (k .eq. par%ztot%nz_aa) then 
        !         write(*,*) k, par%ztot%zeta_aa(k), -9999.0, par%ztot%dzeta_a(k), par%ztot%dzeta_b(k)
        !     else 
        !         write(*,*) k, par%ztot%zeta_aa(k), par%ztot%zeta_ac(k), par%ztot%dzeta_a(k), par%ztot%dzeta_b(k) 
        !     end if  

        ! end do 
        ! write(*,*) "========="
        ! do k = size(par%lith_zeta_aa), 1, -1
        !     write(*,*) par%lith_zeta_aa(k), par%lith_zeta_ac(k), par%lith_dzeta_a(k), par%lith_dzeta_b(k) 
        ! end do 
        ! write(*,*) "========="
        stop 

        ! =================

        ! Define current time as unrealistic value
        par%time = 1000000000   ! [a] 1 billion years in the future

        ! Intialize timestepping parameters to Forward Euler (beta2=0: no contribution from previous timestep)
        par%dt_zeta     = 1.0 
        par%dt_beta(1)  = 1.0 
        par%dt_beta(2)  = 0.0 

        return

    end subroutine ytherm_par_load

    subroutine ytherm_alloc(now,nx,ny,nz_aa,nz_aa_ice,nz_aa_lith)

        implicit none 

        type(ytherm_state_class), intent(INOUT) :: now 
        integer, intent(IN) :: nx, ny
        integer, intent(IN) :: nz_aa
        integer, intent(IN) :: nz_aa_ice
        integer, intent(IN) :: nz_aa_lith

        call ytherm_dealloc(now)

        allocate(now%E_tot(nx,ny,nz_aa))
        allocate(now%T_tot(nx,ny,nz_aa))
        allocate(now%w_tot(nx,ny,nz_aa))
        allocate(now%T_pmp_tot(nx,ny,nz_aa))
        allocate(now%cp_tot(nx,ny,nz_aa))
        allocate(now%kt_tot(nx,ny,nz_aa))
        allocate(now%Q_strn_tot(nx,ny,nz_aa))
        allocate(now%advecxy_tot(nx,ny,nz_aa))
        
        allocate(now%E_ice(nx,ny,nz_aa_ice))
        allocate(now%T_ice(nx,ny,nz_aa_ice))
        allocate(now%w_ice(nx,ny,nz_aa_ice))
        allocate(now%T_pmp(nx,ny,nz_aa_ice))
        allocate(now%cp(nx,ny,nz_aa_ice))
        allocate(now%kt(nx,ny,nz_aa_ice))
        allocate(now%Q_strn(nx,ny,nz_aa_ice))
        allocate(now%advecxy(nx,ny,nz_aa_ice))
        
        allocate(now%E_lith(nx,ny,nz_aa_lith))
        allocate(now%T_lith(nx,ny,nz_aa_lith))

        allocate(now%Q_b(nx,ny))
        allocate(now%Q_ice_b(nx,ny))
        allocate(now%Q_lith(nx,ny))
        
        allocate(now%bmb_grnd(nx,ny))
        allocate(now%f_pmp(nx,ny))
        allocate(now%T_prime_b(nx,ny))
        
        allocate(now%H_cts(nx,ny))
        
        allocate(now%H_w(nx,ny))
        allocate(now%dHwdt(nx,ny))
        
        now%E_tot       = 0.0_wp
        now%T_tot       = 0.0_wp
        now%w_tot       = 0.0_wp
        now%T_pmp_tot   = 0.0_wp
        now%cp_tot      = 0.0_wp
        now%kt_tot      = 0.0_wp
        now%Q_strn_tot  = 0.0_wp
        now%advecxy_tot = 0.0_wp
        
        now%E_ice       = 0.0_wp
        now%T_ice       = 0.0_wp
        now%w_ice       = 0.0_wp
        now%T_pmp       = 0.0_wp
        now%cp          = 0.0_wp
        now%kt          = 0.0_wp
        now%Q_strn      = 0.0_wp
        now%advecxy     = 0.0_wp
        
        now%E_lith      = 0.0_wp
        now%T_lith      = 0.0_wp
        
        now%Q_b         = 0.0_wp
        now%Q_ice_b     = 0.0_wp
        now%Q_lith      = 0.0_wp
        
        now%bmb_grnd    = 0.0_wp
        now%f_pmp       = 0.0_wp
        now%T_prime_b   = 0.0_wp
        
        now%H_cts       = 0.0_wp
        
        now%H_w         = 0.0_wp 
        now%dHwdt       = 0.0_wp
        
        return

    end subroutine ytherm_alloc

    subroutine ytherm_dealloc(now)

        implicit none 

        type(ytherm_state_class), intent(INOUT) :: now

        if (allocated(now%E_tot))           deallocate(now%E_tot)
        if (allocated(now%T_tot))           deallocate(now%T_tot)
        if (allocated(now%w_tot))           deallocate(now%w_tot)
        if (allocated(now%T_pmp_tot))       deallocate(now%T_pmp_tot)
        if (allocated(now%cp_tot))          deallocate(now%cp_tot)
        if (allocated(now%kt_tot))          deallocate(now%kt_tot)
        if (allocated(now%Q_strn_tot))      deallocate(now%Q_strn_tot)
        if (allocated(now%advecxy_tot))     deallocate(now%advecxy_tot)
        
        if (allocated(now%E_ice))           deallocate(now%E_ice)
        if (allocated(now%T_ice))           deallocate(now%T_ice)
        if (allocated(now%w_ice))           deallocate(now%w_ice)
        if (allocated(now%T_pmp))           deallocate(now%T_pmp)
        if (allocated(now%cp))              deallocate(now%cp)
        if (allocated(now%kt))              deallocate(now%kt)
        if (allocated(now%Q_strn))          deallocate(now%Q_strn)
        if (allocated(now%advecxy))         deallocate(now%advecxy)
        
        if (allocated(now%E_lith))          deallocate(now%E_lith)
        if (allocated(now%T_lith))          deallocate(now%T_lith)

        if (allocated(now%Q_b))             deallocate(now%Q_b)
        if (allocated(now%Q_ice_b))         deallocate(now%Q_ice_b)
        if (allocated(now%Q_lith))          deallocate(now%Q_lith)
        
        if (allocated(now%bmb_grnd))        deallocate(now%bmb_grnd)
        if (allocated(now%f_pmp))           deallocate(now%f_pmp)
        if (allocated(now%T_prime_b))       deallocate(now%T_prime_b)
        
        if (allocated(now%H_cts))           deallocate(now%H_cts)
        
        if (allocated(now%H_w))             deallocate(now%H_w)
        if (allocated(now%dHwdt))           deallocate(now%dHwdt)
        
        return 

    end subroutine ytherm_dealloc

end module yelmo_thermodynamics
