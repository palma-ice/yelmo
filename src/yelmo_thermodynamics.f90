
module yelmo_thermodynamics

    use nml 
    use yelmo_defs 
    use yelmo_tools, only : smooth_gauss_2D, smooth_gauss_3D, gauss_values, fill_borders_3D
    
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
!             call calc_strain_heating_smooth(thrm%now%Q_strn,mat%now%strn%de,mat%now%visc,thrm%now%cp,rho_ice)

        end if 
        
        ! Smooth strain heating 
        if (thrm%par%n_sm_qstrn .gt. 0) then 
            call smooth_gauss_3D(thrm%now%Q_strn,tpo%now%H_ice.gt.0.0,thrm%par%dx,thrm%par%n_sm_qstrn, &
                                    tpo%now%H_ice.gt.0.0)
        end if 
        
        ! Calculate the basal frictional heating 
        call calc_basal_heating(thrm%now%Q_b,dyn%now%ux_b,dyn%now%uy_b,dyn%now%taub_acx,dyn%now%taub_acy)

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
                                thrm%now%Q_b,bnd%Q_geo,bnd%T_srf,tpo%now%H_ice,bnd%H_w,tpo%now%H_grnd,tpo%now%f_grnd,thrm%par%zeta_aa, &
                                thrm%par%zeta_ac,thrm%par%dzeta_a,thrm%par%dzeta_b,thrm%par%enth_cr,thrm%par%omega_max, &
                                dt,thrm%par%dx,thrm%par%method,thrm%par%solver_advec)
                    
                case("robin")
                    ! Use Robin solution for ice temperature 

                    call define_temp_robin_3D(thrm%now%T_ice,thrm%now%T_pmp,thrm%now%cp,thrm%now%kt, &
                                       bnd%Q_geo,bnd%T_srf,tpo%now%H_ice,bnd%H_w,bnd%smb, &
                                       thrm%now%bmb_grnd,tpo%now%f_grnd,thrm%par%zeta_aa,cold=.FALSE.)

                    ! Also populate enthalpy 
                    call convert_to_enthalpy(thrm%now%enth,thrm%now%T_ice,thrm%now%omega,thrm%now%T_pmp, &
                                            thrm%now%cp,L_ice)

                case("robin-cold")
                    ! Use Robin solution for ice temperature averaged with cold linear profile
                    ! to ensure cold ice at the base

                    call define_temp_robin_3D(thrm%now%T_ice,thrm%now%T_pmp,thrm%now%cp,thrm%now%kt, &
                                       bnd%Q_geo,bnd%T_srf,tpo%now%H_ice,bnd%H_w,bnd%smb, &
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

        end if 

        ! Calculate homologous temperature at the base 
        thrm%now%T_prime_b = thrm%now%T_ice(:,:,1) - thrm%now%T_pmp(:,:,1) 
        
        ! Calculate gridpoint fraction at the pressure melting point
        thrm%now%f_pmp = calc_f_pmp(thrm%now%T_ice(:,:,1),thrm%now%T_pmp(:,:,1),thrm%par%gamma,tpo%now%f_grnd)

!         if (yelmo_write_log) then 
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
                            T_srf,H_ice,H_w,H_grnd,f_grnd,zeta_aa,zeta_ac,dzeta_a,dzeta_b,cr,omega_max,dt,dx,solver,solver_advec)
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
        real(prec) :: filter0(3,3), filter(3,3) 

        real(prec), parameter :: H_ice_thin = 15.0   ! [m] Threshold to define 'thin' ice

        nx    = size(T_ice,1)
        ny    = size(T_ice,2)
        nz_aa = size(zeta_aa,1)
        nz_ac = size(zeta_ac,1)

        allocate(advecxy(nz_aa))
        allocate(T_ice_old(nx,ny,nz_aa))

        ! First perform horizontal advection (this doesn't work properly, 
        ! use column-based upwind horizontal advection below)
        !call calc_enth_horizontal_advection_3D(T_ice,ux,uy,dx,dt,solver_advec)

        ! Store original ice temperature field here for input to horizontal advection
        ! calculations 
        T_ice_old = T_ice 

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
                call calc_advec_horizontal_column(advecxy,T_ice_old,H_ice,ux,uy,dx,i,j)
                !advecxy = 0.0_prec 

                call calc_enth_column(enth(i,j,:),T_ice(i,j,:),omega(i,j,:),bmb_grnd(i,j),Q_ice_b(i,j),H_cts(i,j), &
                        T_pmp(i,j,:),cp(i,j,:),kt(i,j,:),advecxy,uz(i,j,:),Q_strn(i,j,:),Q_b(i,j),Q_geo(i,j),T_srf(i,j), &
                        T_shlf,H_ice_now,H_w(i,j),f_grnd(i,j),zeta_aa,zeta_ac,dzeta_a,dzeta_b,cr,omega_max,T0,dt,trim(solver))
                
            end if 

        end do 
        end do 

        ! Fill in borders 
        call fill_borders_3D(enth,nfill=2)
        call fill_borders_3D(T_ice,nfill=2)
        call fill_borders_3D(omega,nfill=2)

        return 

    end subroutine calc_ytherm_enthalpy_3D

    subroutine calc_enth_horizontal_advection_3D(T_ice,ux,uy,dx,dt,solver)

        implicit none 

        real(prec),       intent(INOUT) :: T_ice(:,:,:)         ! [K]   Ice temperature/enthalpy, aa-nodes  
        real(prec),       intent(IN)    :: ux(:,:,:)            ! [m/a] 2D velocity, x-direction (ac-nodes)
        real(prec),       intent(IN)    :: uy(:,:,:)            ! [m/a] 2D velocity, y-direction (ac-nodes)
        real(prec),       intent(IN)    :: dx                   ! [m]   Horizontal resolution
        real(prec),       intent(IN)    :: dt                   ! [a]   Timestep 
        character(len=*), intent(IN)    :: solver               ! Solver to use for the ice thickness advection equation

        ! Local variables 
        integer :: i, j, k, nx, ny, nz 
        real(prec), allocatable :: T_dot(:,:) 

        nx = size(T_ice,1)
        ny = size(T_ice,2)
        nz = size(T_ice,3) 

        allocate(T_dot(nx,ny)) 
        T_dot = 0.0_prec 

        ! Resolve horizontal advection layer by layer 

        do k = 2, nz-1    
            call calc_advec2D(T_ice(:,:,k),ux(:,:,k),uy(:,:,k),T_dot,dx,dx,dt,solver)
        end do 

        call fill_borders_3D(T_ice,nfill=2)
        
        return 

    end subroutine calc_enth_horizontal_advection_3D

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
        
        return 

    end subroutine ytherm_dealloc 

end module yelmo_thermodynamics
