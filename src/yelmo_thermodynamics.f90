
module yelmo_thermodynamics

    use nml 
    use yelmo_defs 

    use thermodynamics 
    use icetemp  
    use enthalpy 
    
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

!mmr
        print*,'hola thrm time', thrm%par%time, time
!mmr

        ! Initialize time if necessary 
        if (thrm%par%time .gt. time) then 
            thrm%par%time = time
        end if 

        ! Get time step and advance current time 
!mmr        dt            = time - thrm%par%time 
        dt            = 1.0 ! mmr trick!!!! time - thrm%par%time 
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
            thrm%now%T_pmp(:,:,k) = calc_T_pmp(tpo%now%H_ice,thrm%par%zeta_aa(k),T0)
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
            call smooth_strain_heating(thrm%now%Q_strn,tpo%now%H_ice,thrm%par%dx,thrm%par%n_sm_qstrn)
        end if 
        
        ! Calculate the basal frictional heating 
        call calc_basal_heating(thrm%now%Q_b,dyn%now%ux_b,dyn%now%uy_b,dyn%now%taub_acx,dyn%now%taub_acy)

        if ( dt .gt. 0.0 ) then    
!mmr not orig        if ( dt .ge. 0.0 ) then     
            ! Ice thermodynamics should evolve, perform calculations 

            select case(trim(thrm%par%method))

                case("active") 
                    ! Perform temperature solving via advection-diffusion equation
                    
                    call calc_icetemp_3D(thrm%now%T_ice,thrm%now%bmb_grnd,thrm%now%dTdz_b,thrm%now%T_pmp, &
                                        thrm%now%cp,thrm%now%kt,dyn%now%ux,dyn%now%uy,dyn%now%uz,thrm%now%Q_strn, &
                                        thrm%now%Q_b,bnd%Q_geo,bnd%T_srf,tpo%now%H_ice,bnd%H_w,tpo%now%H_grnd,tpo%now%f_grnd, &
                                        thrm%par%zeta_aa,thrm%par%zeta_ac,thrm%par%dzeta_a,thrm%par%dzeta_b,dt,thrm%par%dx)
                
                case("robin")
                    ! Use Robin solution for ice temperature 

                    call define_temp_robin_3D(thrm%now%T_ice,thrm%now%T_pmp,thrm%now%cp,thrm%now%kt, &
                                       bnd%Q_geo,bnd%T_srf,tpo%now%H_ice,bnd%H_w,bnd%smb, &
                                       thrm%now%bmb_grnd,tpo%now%f_grnd,thrm%par%zeta_aa,cold=.FALSE.)

                case("robin-cold")
                    ! Use Robin solution for ice temperature averaged with cold linear profile
                    ! to ensure cold ice at the base

                    call define_temp_robin_3D(thrm%now%T_ice,thrm%now%T_pmp,thrm%now%cp,thrm%now%kt, &
                                       bnd%Q_geo,bnd%T_srf,tpo%now%H_ice,bnd%H_w,bnd%smb, &
                                       thrm%now%bmb_grnd,tpo%now%f_grnd,thrm%par%zeta_aa,cold=.TRUE.)

                case("linear")
                    ! Use linear solution for ice temperature

                    ! Calculate the ice temperature (eventually water content and enthalpy too)
                    call define_temp_linear_3D(thrm%now%T_ice,thrm%par%zeta_aa,tpo%now%H_ice,bnd%T_srf)

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

!mmr
        print*,'hola yelmo_write_log'
        yelmo_write_log = .TRUE.
!mmr
        if (yelmo_write_log) then 
            if (count(tpo%now%H_ice.gt.0.0) .gt. 0) then 
                write(*,"(a,f14.4,f10.4,f10.2)") "calc_ytherm:: time = ", thrm%par%time, dt, &
!mmr                    sum(thrm%now%T_ice(:,:,thrm%par%nz_aa),mask=tpo%now%H_ice.gt.0.0)/real(count(tpo%now%H_ice.gt.0.0))
                    sum(thrm%now%T_ice(:,:,1),mask=tpo%now%H_ice.gt.0.0)/real(count(tpo%now%H_ice.gt.0.0))
            else 
                write(*,"(a,f14.4,f10.4,f10.2)") "calc_ytherm:: time = ", thrm%par%time, dt, 0.0 
            end if 
        end if 

        return

    end subroutine calc_ytherm

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
        call nml_read(filename,"ytherm","method",         par%method,     init=init_pars)
        call nml_read(filename,"ytherm","cond_bed",       par%cond_bed,   init=init_pars)
        call nml_read(filename,"ytherm","gamma",          par%gamma,      init=init_pars)
        call nml_read(filename,"ytherm","nzr",            par%nzr,        init=init_pars)
        call nml_read(filename,"ytherm","H_rock",         par%H_rock,     init=init_pars)
        
        call nml_read(filename,"ytherm","n_sm_qstrn", par%n_sm_qstrn, init=init_pars)
        call nml_read(filename,"ytherm","use_strain_sia", par%use_strain_sia, init=init_pars)
        call nml_read(filename,"ytherm","use_const_cp",   par%use_const_cp,   init=init_pars)
        call nml_read(filename,"ytherm","const_cp",       par%const_cp,       init=init_pars)
        call nml_read(filename,"ytherm","use_const_kt",   par%use_const_kt,   init=init_pars)
        call nml_read(filename,"ytherm","const_kt",       par%const_kt,       init=init_pars)
        
        call nml_read(filename,"ytherm","kt_m",           par%kt_m,       init=init_pars)
        call nml_read(filename,"ytherm","cp_m",           par%cp_m,       init=init_pars)
        call nml_read(filename,"ytherm","rho_m",          par%rho_m,      init=init_pars)
        
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
        
        if (allocated(par%zetar)) deallocate(par%zetar)
        allocate(par%zetar(par%nzr))
        
        ! Bottom = 0.0, Bedrock surface = 1.0
        ! par%zetar(par%nzr) = 1.0 is the same as bottom of ice layer, 
        ! ie par%zetar(par%nzr)==par%zeta_aa(1) represent the same location
        ! Test, have zetar end at less than 1.0, so that the surface of hte rock is separate
        ! from the ice 
        do k = 1, par%nzr 
            par%zetar(k) = 0.0 + 1.0*(k-1)/(par%nzr-1)
!             par%zetar(k) = 0.0 + 1.0*(k-1)/(par%nzr)
        end do 

        if (allocated(par%dzeta_a)) deallocate(par%dzeta_a)
        allocate(par%dzeta_a(par%nz_aa))
        if (allocated(par%dzeta_b)) deallocate(par%dzeta_b)
        allocate(par%dzeta_b(par%nz_aa))
        
        ! Define thermodynamic zeta helper variables dzeta_a/dzeta_b
        call calc_dzeta_terms(par%dzeta_a,par%dzeta_b,par%zeta_aa,par%zeta_ac)

        ! Define current time as unrealistic value
        par%time = 1000000000   ! [a] 1 billion years in the future

        ! Initialize enthalpy tables
        call enthalpy_init_tables() 

        return

    end subroutine ytherm_par_load

    subroutine ytherm_alloc(now,nx,ny,nz_aa,nz_ac,nzr)

        implicit none 

        type(ytherm_state_class), intent(INOUT) :: now 
        integer, intent(IN) :: nx, ny, nz_aa, nz_ac, nzr 

        call ytherm_dealloc(now)

        allocate(now%T_ice(nx,ny,nz_aa))
        allocate(now%T_pmp(nx,ny,nz_aa))
        
        allocate(now%cp(nx,ny,nz_aa))
        allocate(now%kt(nx,ny,nz_aa))
        allocate(now%Q_strn(nx,ny,nz_aa))
        
        allocate(now%enth_ice(nx,ny,nz_aa))
        allocate(now%omega_ice(nx,ny,nz_aa))

        allocate(now%phid(nx,ny))
        allocate(now%bmb_grnd(nx,ny))
        allocate(now%f_pmp(nx,ny))
        allocate(now%Q_b(nx,ny))
        
        allocate(now%T_rock(nx,ny,nzr))
        
        allocate(now%dTdz_b(nx,ny))
        allocate(now%dTrdz_b(nx,ny))
        allocate(now%T_prime_b(nx,ny))
        allocate(now%cts(nx,ny))
        allocate(now%T_all(nx,ny,nzr+nz_aa))

        now%T_ice     = 0.0
        now%enth_ice  = 0.0
        now%omega_ice = 0.0  
        now%T_pmp     = 0.0
        now%phid      = 0.0   
        now%bmb_grnd  = 0.0 
        now%f_pmp     = 0.0 
        now%Q_strn    = 0.0 
        now%Q_b       = 0.0 
        now%cp        = 0.0 
        now%kt        = 0.0 
        now%T_rock    = 0.0
        
        now%dTdz_b    = 0.0 
        now%dTrdz_b   = 0.0 
        now%T_prime_b = 0.0 
        now%cts       = 0.0 
        now%T_all     = 0.0 

        return 
    end subroutine ytherm_alloc 

    subroutine ytherm_dealloc(now)

        implicit none 

        type(ytherm_state_class), intent(INOUT) :: now

        if (allocated(now%T_ice))     deallocate(now%T_ice)
        if (allocated(now%T_pmp))     deallocate(now%T_pmp)
        if (allocated(now%bmb_grnd))  deallocate(now%bmb_grnd)
        if (allocated(now%f_pmp))     deallocate(now%f_pmp)
        if (allocated(now%Q_strn))    deallocate(now%Q_strn)
        if (allocated(now%Q_b))       deallocate(now%Q_b)
        if (allocated(now%cp))        deallocate(now%cp)
        if (allocated(now%kt))        deallocate(now%kt)
        
        if (allocated(now%T_rock))    deallocate(now%T_rock)
        
        if (allocated(now%dTdz_b))    deallocate(now%dTdz_b)
        if (allocated(now%dTrdz_b))   deallocate(now%dTrdz_b)
        if (allocated(now%T_prime_b)) deallocate(now%T_prime_b)
        
        if (allocated(now%T_all))     deallocate(now%T_all)
        
        return 

    end subroutine ytherm_dealloc 

end module yelmo_thermodynamics
