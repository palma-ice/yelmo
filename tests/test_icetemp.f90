program test_icetemp 

    use ncio  
    
    use yelmo_defs
    use thermodynamics 
    use ice_enthalpy  
    use ice_tracer

    use interp1D 

    implicit none 

    type zeta_class 
        real(prec), allocatable :: zeta(:)    ! [-] Sigma coordinates from 0:1 (either height or depth)
        real(prec), allocatable :: zeta_ac(:) ! nz-1 [-] sigma coordinates for internal ice layer edges
        real(prec), allocatable :: dzeta_a(:) ! [-] Sigma coordinates derivative term
        real(prec), allocatable :: dzeta_b(:) ! [-] Sigma coordinates derivative term
    end type 

    type icesheet_vectors
        real(prec), allocatable :: enth(:)    ! [J kg-1] Ice enthalpy
        real(prec), allocatable :: T_ice(:)   ! [K] Ice temperature 
        real(prec), allocatable :: omega(:)   ! [-] Ice water content (fraction)
        real(prec), allocatable :: T_pmp(:)   ! [K] Ice pressure melting point 
        
        real(prec), allocatable :: cp(:)      ! [] Ice heat capacity 
        real(prec), allocatable :: kt(:)      ! [] Ice conductivity  
        
        real(prec), allocatable :: uz(:)      ! [m a-1] Vertical velocity 
        real(prec), allocatable :: advecxy(:) ! [] Horizontal heat advection magnitude
        real(prec), allocatable :: Q_strn(:)  ! [K a-1] Strain heating 
        real(prec), allocatable :: t_dep(:)   ! [a] Deposition time 
        
        real(prec), allocatable :: enth_rock(:) ! [J kg-1] Bedrock enthalpy
        real(prec), allocatable :: T_rock(:)    ! [K] Bedrock temperature 

    end type 

    type icesheet 

        type(zeta_class) :: z           ! Ice column vertical axis
        type(zeta_class) :: zr          ! Bedrock column vertical axis
        
        type(icesheet_vectors) :: vec   ! For height coordinate systems with k=1 base and k=nz surface
        
        real(prec) :: H_w               ! [m] Water present at the ice base 
        real(prec) :: Q_ice_b           ! [mW m-2] Ice basal heat flux (positive up)
        real(prec) :: Q_rock            ! [mW m-2] Bedrock surface heat flux (positive up)
        real(prec) :: Q_geo             ! [mW m-2] Geothermal heat flux 
        real(prec) :: Q_b               ! [mW m-2] Basal heat production
        real(prec) :: H_cts             ! [m] cold-temperate transition surface (CTS) height
        real(prec) :: H_ice             ! [m] Ice thickness 
        real(prec) :: cp_rock           ! [] Bedrock heat capacity
        real(prec) :: kt_rock           ! [] Bedrock conductivity 
        real(prec) :: H_rock            ! [m] Bedrock thickness 
        
        real(prec) :: T_srf             ! [K] Ice surface temperature 
        real(prec) :: T_shlf            ! [K] Ice shelf base temperature 
        real(prec) :: smb               ! [m a**-1] Surface mass balance
        real(prec) :: bmb               ! [m a**-1] Basal mass balance
        real(prec) :: f_grnd            ! [-] Grounded fraction 
        
        character(len=56) :: age_method ! Method to use for age calculation 
        real(prec) :: age_impl_kappa    ! [m2 a-1] Artificial diffusion term for implicit age solving 
        end type 

    ! Define different icesheet objects for use in proram
    type(icesheet) :: ice1
    type(icesheet) :: robin 
    type(icesheet) :: diff 
    
    character(len=56)  :: rock_method   

    ! Local variables
    real(prec)         :: t_start, t_end, dt, time  
    integer            :: n, ntot  
    character(len=512) :: file1D 
    real(prec)         :: dt_out 
    character(len=56)  :: zeta_scale 
    integer            :: nz 
    character(len=56)  :: zeta_scale_rock 
    integer            :: nzr 
    logical            :: is_celcius 
    character(len=56)  :: age_method 
    real(prec)         :: age_impl_kappa
    character(len=12)  :: enth_solver  
    real(prec)         :: enth_cr
    real(prec)         :: omega_max 
    character(len=56)  :: experiment 

    real(prec)         :: T0_ref 

    integer            :: narg 
    character(len=12)  :: arg_nz, arg_cr
    character(len=32)  :: nz_str, cr_str, prec_str 

    integer :: iter, n_iter 

    ! General initialization of yelmo constants (used globally)
    call yelmo_global_init("par/yelmo_const_EISMINT.nml")
    
    ! ===============================================================
    ! User options 

    experiment     = "eismint"      ! "eismint", "k15expa", "k15expb", "bg15a"
    
    ! General options
    zeta_scale      = "linear"      ! "linear", "exp", "tanh"
    nz              = 22            ! [--] Number of ice sheet points (aa-nodes + base + surface)
    is_celcius      = .FALSE. 

    zeta_scale_rock = "exp-inv" 
    nzr             = 5 
    rock_method     = "equil"       ! "equil" or "active" bedrock 

    age_method      = "expl"        ! "expl" or "impl"
    age_impl_kappa  = 1.5           ! [m2 a-1] Artificial diffusion for age tracing

    enth_solver     = "temp"        ! "enth" or "temp" 
    omega_max       = 0.03          ! Maximum allowed water content (fraction)
    enth_cr         = 1e-3          ! Enthalpy solver: conductivity ratio kappa_water / kappa_ice 

    
    ! Overwrite options for nz and enth_cr if available from arguments
    narg = command_argument_count() 
    if (narg .gt. 0) then 
        call get_command_argument(1,arg_nz)
        read(arg_nz,*)  nz
    end if 

    if (narg .gt. 1) then 
        call get_command_argument(2,arg_cr)
        read(arg_cr,*)  enth_cr
    end if 

    ! ## Save simulation parameters as strings ##
    write(nz_str,*) nz 
    if (trim(zeta_scale) .ne. "linear") write(nz_str,*) trim(nz_str)//trim(zeta_scale)
    nz_str =  adjustl(nz_str)

    write(cr_str,"(e8.2)") enth_cr
    cr_str =  adjustl(cr_str)   

    prec_str = "sp" 
    if (prec .eq. dp) prec_str = "dp" 

    ! Use a more precise filename to specify cr value and dz
    !file1D = "output/test_"//trim(experiment)//".nc"
    write(file1D,*) "output/test_"//trim(experiment)//"_nz",trim(nz_str),   &
                                        "_cr",trim(cr_str),"_",trim(prec_str),".nc"
    
    ! ===============================================================

    T0_ref = T0 
    if (is_celcius) T0_ref = 0.0 

    ! Initialize icesheet object 
    call icesheet_allocate(ice1,nz,nzr,zeta_scale,zeta_scale_rock) 

    select case(trim(experiment))

        case("k15expa")

            t_start = 0.0       ! [yr]
            t_end   = 300e3     ! [yr]
            dt      = 5.0       ! [yr]
            dt_out  = 500.0     ! [yr] 

            T_pmp_beta = 7.9e-8         ! [K Pa^-1] Kleiner et al. (2015), expa

            call init_k15expa(ice1)

        case("k15expb")

            t_start = 0.0       ! [yr]
            t_end   = 1e3       ! [yr]
            dt      = 0.2       ! [yr]
            dt_out  = 20.0      ! [yr] 

            T_pmp_beta = 0.0            ! [K Pa^-1] Kleiner et al. (2015), expb
            
            call init_k15expb(ice1,smb=0.2_prec,T_srf=-3.0_prec)

        case("bg15a")

            t_start = -0.5e3    ! [yr]
            t_end   =  1e3      ! [yr]
            dt      = 0.5       ! [yr]
            dt_out  =  5.0      ! [yr] 

            T_pmp_beta = 0.0            ! [K Pa^-1] Blatter and Greve (2015), expa
            
            call init_bg15a(ice1,smb=0.2_prec,T_srf=-4.0_prec)


        case DEFAULT 
            ! EISMINT 

            t_start = 0.0       ! [yr]
            t_end   = 300e3     ! [yr]
            dt      = 0.5_prec  ! [yr]
            dt_out  = 1000.0    ! [yr] 

            !T_pmp_beta = 9.8e-8         ! [K Pa^-1] Greve and Blatter (2009) 
            T_pmp_beta = 9.7e-8         ! [K Pa^-1] EISMINT2 value (beta1 = 8.66e-4 [K m^-1])

            call init_eismint_summit(ice1,smb=0.5_prec)

    end select 

    ! Initialize time and calculate number of time steps to iterate and 
    time = t_start 
    ntot = (t_end-t_start)/dt 
    
    ! Set age method and kappa 
    ice1%age_method     = age_method 
    ice1%age_impl_kappa = age_impl_kappa

    ! Calculate the robin solution for comparison 
    robin = ice1  
    robin%vec%T_ice = calc_temp_robin_column(robin%z%zeta,robin%vec%T_pmp,robin%vec%kt,robin%vec%cp,rho_ice, &
                                       robin%H_ice,robin%T_srf,robin%smb,robin%Q_geo,is_float=robin%f_grnd.eq.0.0)

    ! Calculate initial enthalpy (robin)
    call convert_to_enthalpy(robin%vec%enth,robin%vec%T_ice,robin%vec%omega,robin%vec%T_pmp,robin%vec%cp,L_ice)

    ! Ensure zero basal water thickness to start 
    ice1%H_w = 0.0 

    ! Assume H_cts is also zero to start
    ice1%H_cts = 0.0 

    ! Initialize output file for model and write intial conditions 
    call write_init(ice1,filename=file1D,zeta=ice1%z%zeta,zeta_ac=ice1%z%zeta_ac, &
                                                zeta_r=ice1%zr%zeta,time_init=time)
    call write_step(ice1,ice1%vec,filename=file1D,time=time)

    ! Loop over time steps and perform thermodynamic calculations
    do n = 1, ntot 

        ! Get current time 
        time = t_start + n*dt 

        if (trim(experiment) .eq. "k15expa") then 
            if (time .le. 100e3) then 
                ice1%T_srf = T0_ref - 30.0
            else if (time .gt. 100e3 .and. time .le. 150e3) then 
                ice1%T_srf = T0_ref - 5.0 
            else   ! time .gt. 150e3
                ice1%T_srf = T0_ref - 30.0
            end if  
        end if 

        if (trim(experiment) .eq. "bg15a") then 

            ! === MELTING === 
            if (time .le. 50.0) then 
                ice1%T_srf = T0_ref - 4.0_prec 
            else 
                ice1%T_srf = T0_ref - 2.0_prec 
            end if 
            
            ! === FREEZING === 
!             if (time .le. 50.0) then 
!                 ice1%T_srf = T0_ref - 10.0_prec 
!             else 
!                 ice1%T_srf = T0_ref - 6.0_prec 
!             end if
            
            ! === MELTING SIN WAVE ===
            if (time .le. 0.0) then 
                ice1%T_srf = T0_ref - 2.0_prec 
            else 
                ice1%T_srf = T0_ref - 2.0_prec + 1.0_prec*sin(2.0*pi*time/100.0_prec)
            end if 
        end if 

        select case(trim(enth_solver))

            case("temp") 

                call calc_temp_column(ice1%vec%enth,ice1%vec%T_ice,ice1%vec%omega,ice1%bmb,ice1%Q_ice_b,ice1%H_cts,ice1%vec%T_pmp, &
                        ice1%vec%cp,ice1%vec%kt,ice1%vec%advecxy,ice1%vec%uz,ice1%vec%Q_strn,ice1%Q_b,ice1%Q_rock,ice1%T_srf,ice1%T_shlf, &
                        ice1%H_ice,ice1%H_w,ice1%f_grnd,ice1%z%zeta,ice1%z%zeta_ac,ice1%z%dzeta_a,ice1%z%dzeta_b,omega_max,T0_ref,dt)

            case("enth")

                call calc_enth_column(ice1%vec%enth,ice1%vec%T_ice,ice1%vec%omega,ice1%bmb,ice1%Q_ice_b,ice1%H_cts,ice1%vec%T_pmp, &
                        ice1%vec%cp,ice1%vec%kt,ice1%vec%advecxy,ice1%vec%uz,ice1%vec%Q_strn,ice1%Q_b,ice1%Q_geo,ice1%T_srf,ice1%T_shlf, &
                        ice1%H_ice,ice1%H_w,ice1%f_grnd,ice1%z%zeta,ice1%z%zeta_ac,ice1%z%dzeta_a,ice1%z%dzeta_b, &
                        enth_cr,omega_max,T0_ref,dt)

            case DEFAULT 

                write(*,*) "enth_solver not recognized: ", trim(enth_solver)

        end select 

        select case(trim(rock_method))

            case("equil")

                ! Define temperature profile in bedrock too 
                call calc_temp_bedrock_column(ice1%vec%T_rock,ice1%kt_rock,rho_rock, &
                                        ice1%H_rock,ice1%vec%T_ice(1),ice1%Q_geo,ice1%zr%zeta)

                call convert_to_enthalpy(ice1%vec%enth_rock,ice1%vec%T_rock,0.0_wp,1e8_wp,ice1%cp_rock,0.0_wp)


            case("active")

                call calc_temp_column_bedrock(ice1%vec%enth_rock,ice1%vec%T_rock,ice1%Q_rock, &
                            ice1%cp_rock,ice1%kt_rock,ice1%Q_ice_b,ice1%Q_geo,ice1%T_srf,ice1%H_rock, &
                            ice1%zr%zeta,ice1%zr%zeta_ac,ice1%zr%dzeta_a,ice1%zr%dzeta_b,dt)
                
            case DEFAULT 

                write(*,*) "rock method not recognized: ", trim(rock_method)

        end select 

        ! Update basal water thickness [m/a i.e.] => [m/a w.e.]
        ice1%H_w = max(ice1%H_w - (ice1%bmb*rho_ice/rho_w)*dt, 0.0_prec)

        if (trim(age_method) .eq. "impl") then 
            call calc_tracer_column(ice1%vec%t_dep,ice1%vec%uz,ice1%vec%advecxy*0.0,time,ice1%bmb, &
                                    ice1%H_ice,ice1%z%zeta,ice1%z%zeta_ac, &
                                    ice1%age_impl_kappa,dt)
        else 
            call calc_tracer_column_expl(ice1%vec%t_dep,ice1%vec%uz,ice1%vec%advecxy*0.0,time,ice1%bmb,ice1%H_ice,ice1%z%zeta,ice1%z%zeta_ac,dt)
        end if 

        if (mod(time,dt_out)==0) then 
            call write_step(ice1,ice1%vec,filename=file1D,time=time,T_robin=robin%vec%T_ice)
        end if 

        if (mod(time,50.0)==0) then
            write(*,"(a,f14.4)") "time = ", time
        end if 

    end do

    write(*,*)
    write(*,*) "========================="
    write(*,*) 
    write(*,*) "Program finished."
    write(*,*)
    write(*,*) "Output filename: ", trim(file1D)
    write(*,*)
    write(*,*) "========================="
    write(*,*)

contains 
    
    subroutine init_eismint_summit(ice,smb)

        implicit none 

        type(icesheet), intent(INOUT) :: ice
        real(prec),     intent(IN)    :: smb 

        ! Local variables 
        integer :: k, nz   

        nz    = size(ice%z%zeta)
    
        ! Assign point values
        ice%T_srf    = 239.0        ! [K]
        ice%T_shlf   = T0           ! [K] T_shlf not used in this idealized setup, set to T0  
        ice%smb      = smb          ! [m/a]
        ice%bmb      = 0.0          ! [m/a]
        ice%Q_geo    = 42.0         ! [mW/m2]
        ice%H_ice    = 2997.0       ! [m] Summit thickness
        ice%H_w      = 0.0          ! [m] No basal water
        ice%Q_b      = 0.0          ! [] No basal frictional heating 
        ice%f_grnd   = 1.0          ! Grounded point 

        ice%cp_rock  = 1000.0       ! [J kg-1 K-1]
        ice%kt_rock  = 9.46e7       ! [J a-1 m-1 K-1]
        ice%H_rock   = 2000.0       ! [m] 
        ice%Q_rock   = ice%Q_geo 
        
        ! EISMINT1
        ice%vec%cp      = 2009.0    ! [J kg-1 K-1]
        ice%vec%kt      = 6.67e7    ! [J a-1 m-1 K-1]
        
        ice%vec%Q_strn  = 0.0       ! [J a-1 m-3] No internal strain heating 
        ice%vec%advecxy = 0.0       ! [] No horizontal advection 

        ! Calculate pressure melting point 
        ice%vec%T_pmp = calc_T_pmp(ice%H_ice,ice%z%zeta,T0,T_pmp_beta) 

        if (is_celcius) then 
            ice%T_srf     = ice%T_srf     - T0
            ice%T_shlf    = ice%T_shlf    - T0
            ice%vec%T_pmp = ice%vec%T_pmp - T0 
        end if 

        ! Define initial temperature profile, linear 
        ! from T_srf at the surface to 10deg below freezing point at the base

        ice%vec%T_ice(nz) = ice%T_srf 
        ice%vec%T_ice(1)  = ice%vec%T_pmp(1) - 10.0 

        ! Intermediate layers are linearly interpolated 
        do k = 2, nz-1 
            ice%vec%T_ice(k) = ice%vec%T_ice(1)+ice%z%zeta(k)*(ice%vec%T_ice(nz)-ice%vec%T_ice(1))
        end do 

        ! Define linear vertical velocity profile
        ice%vec%uz = -ice%smb*ice%z%zeta_ac 
        
        ! Calculate initial enthalpy (ice1)
        call convert_to_enthalpy(ice%vec%enth,ice%vec%T_ice,ice%vec%omega,ice%vec%T_pmp,ice%vec%cp,L_ice)
        
        ! Define temperature profile in bedrock too 
        call calc_temp_bedrock_column(ice1%vec%T_rock,ice1%kt_rock,rho_rock, &
                                    ice1%H_rock,ice1%vec%T_ice(1),ice1%Q_geo,ice1%zr%zeta)

        call convert_to_enthalpy(ice%vec%enth_rock,ice%vec%T_rock,0.0_wp,1e8_wp,ice%cp_rock,0.0_wp)

        return 

    end subroutine init_eismint_summit

    subroutine init_k15expa(ice)

        implicit none 

        type(icesheet), intent(INOUT) :: ice

        ! Local variables 
        integer :: k, nz 

        nz    = size(ice%z%zeta)

        ! Assign point values
        ice%T_srf    = T0 - 30.0        ! [K]
        ice%T_shlf   = T0               ! [K] T_shlf not used in this idealized setup, set to T0  
        ice%smb      = 0.0              ! [m/a]
        ice%bmb      = 0.0              ! [m/a]
        ice%Q_geo    = 42.0             ! [mW/m2]
        ice%H_ice    = 1000.0           ! [m] Ice thickness
        ice%H_w      = 0.0              ! [m] No basal water
        ice%Q_b      = 0.0              ! [] No basal frictional heating 
        ice%f_grnd   = 1.0              ! Grounded point 

        ice%cp_rock  = 1000.0    ! [J kg-1 K-1]
        ice%kt_rock  = 9.46e7    ! [J a-1 m-1 K-1]
        ice%H_rock   = 2000.0       ! [m] 
        ice%Q_rock   = ice%Q_geo 
        
        ! EISMINT1
        ice%vec%cp      = 2009.0        ! [J kg-1 K-1]
        ice%vec%kt      = 6.6269e7      ! [J a-1 m-1 K-1]   == 2.1*sec_year  [J s-1 m-1 K-1] => J a-1 m-1 K-1]
        
        ice%vec%Q_strn  = 0.0           ! [] No internal strain heating 
        ice%vec%advecxy = 0.0           ! [] No horizontal advection 

        ! Calculate pressure melting point 
        ice%vec%T_pmp = calc_T_pmp(ice%H_ice,ice%z%zeta,T0,T_pmp_beta) 

        if (is_celcius) then 
            ice%T_srf     = ice%T_srf     - T0
            ice%T_shlf    = ice%T_shlf    - T0
            ice%vec%T_pmp = ice%vec%T_pmp - T0 
        end if 

        ! Define initial temperature profile
        ! (constant equal to surface temp)
        ice%vec%T_ice(nz) = ice%T_srf 
        ice%vec%T_ice(1)  = ice%T_srf

        ! Intermediate layers are linearly interpolated 
        do k = 2, nz-1 
            ice%vec%T_ice(k) = ice%vec%T_ice(1)+ice%z%zeta(k)*(ice%vec%T_ice(nz)-ice%vec%T_ice(1))
        end do 

        ! Define linear vertical velocity profile
        ice%vec%uz = -ice%smb*ice%z%zeta_ac 

        ! Calculate initial enthalpy (ice1)
        call convert_to_enthalpy(ice%vec%enth,ice%vec%T_ice,ice%vec%omega,ice%vec%T_pmp,ice%vec%cp,L_ice)
        
        ! Define temperature profile in bedrock too 
        call calc_temp_bedrock_column(ice1%vec%T_rock,ice1%kt_rock,rho_rock, &
                                    ice1%H_rock,ice1%vec%T_ice(1),ice1%Q_geo,ice1%zr%zeta)

        call convert_to_enthalpy(ice%vec%enth_rock,ice%vec%T_rock,0.0_wp,1e8_wp,ice%cp_rock,0.0_wp)

        return 

    end subroutine init_k15expa
    
    subroutine init_k15expb(ice,smb,T_srf)

        implicit none 

        type(icesheet), intent(INOUT) :: ice
        real(prec),     intent(IN)    :: smb 
        real(prec),     intent(IN)    :: T_srf ! [degrees Celcius]
            
        ! Local variables 
        integer :: k, nz 
        real(prec) :: ATT, gamma, T_init    
        real(prec), allocatable :: ux(:)
        real(prec), allocatable :: uy(:)  
        real(prec), allocatable :: mu(:)
        real(prec), allocatable :: eps(:) 

        nz    = size(ice%z%zeta)

        ! Assign point values
        ice%T_srf    = T_srf + T0       ! [K]
        ice%T_shlf   = T0               ! [K] T_shlf not used in this idealized setup, set to T0  
        ice%smb      = smb              ! [m/a]
        ice%bmb      = 0.0              ! [m/a]
        ice%Q_geo    = 0.0              ! [mW/m2]
        ice%H_ice    = 200.0            ! [m] Ice thickness
        ice%H_w      = 0.0              ! [m] No basal water
        ice%Q_b      = 0.0              ! [] No basal frictional heating 
        ice%f_grnd   = 1.0              ! Grounded point 

        ice%cp_rock  = 1000.0    ! [J kg-1 K-1]
        ice%kt_rock  = 9.46e7    ! [J a-1 m-1 K-1]
        ice%H_rock   = 2000.0       ! [m] 
        ice%H_rock   = 2000.0           ! [m] 
        ice%Q_rock   = ice%Q_geo 
        
        ! EISMINT1
        ice%vec%cp      = 2009.0        ! [J kg-1 K-1]
        ice%vec%kt      = 6.6269e7      ! [J a-1 m-1 K-1]   == 2.1*sec_year  [J s-1 m-1 K-1] => J a-1 m-1 K-1]

        ATT             = 5.3e-24*sec_year      ! Rate factor
        gamma           = 4.0                   ! [degrees] Bed slope 
        T_init          = T0 - 1.5 

        allocate(ux(nz))
        allocate(uy(nz))
        allocate(eps(nz))
        allocate(mu(nz))

        ux = 0.5*ATT*(rho_ice*g*sin(gamma*pi/180.0))**3 * (ice%H_ice**4 - (ice%H_ice - ice%H_ice*ice%z%zeta)**4)
        uy = 0.0 

        eps = ATT*(rho_ice*g*sin(gamma*pi/180.0))**3 *(ice%H_ice - ice%H_ice*ice%z%zeta)**3
        mu  = 0.5*ATT**(-1.0/3.0)*eps**(-2.0/3.0)

        ice%vec%Q_strn(nz)  = 0.0  
        ice%vec%Q_strn(1:nz-1)  = 4.0*mu(1:nz-1)*eps(1:nz-1)**2     ! [J a-1 m-3] Prescribed strain heating 
        ice%vec%advecxy = 0.0                                       ! [] No horizontal advection (assume constant)

        ! Write strain heating to compare basal value of ~2.6e-3 W/m-3
        !do k = nz, 1, -1 
        !    write(*,*) ice%z%zeta(k), ice%vec%Q_strn(k)/sec_year 
        !end do 

        ! Calculate pressure melting point 
        ice%vec%T_pmp = calc_T_pmp(ice%H_ice,ice%z%zeta,T0,T_pmp_beta) 

        if (is_celcius) then 
            ice%T_srf     = ice%T_srf     - T0
            ice%T_shlf    = ice%T_shlf    - T0
            ice%vec%T_pmp = ice%vec%T_pmp - T0 
        end if 

        ! Define initial temperature profile
        ! (constant equal to surface temp)
        ice%vec%T_ice(nz) = T_init
        ice%vec%T_ice(1)  = T_init 

        ! Intermediate layers are linearly interpolated 
        do k = 2, nz-1 
            ice%vec%T_ice(k) = ice%vec%T_ice(1)+ice%z%zeta(k)*(ice%vec%T_ice(nz)-ice%vec%T_ice(1))
        end do 

        ! Define constant vertical velocity profile
        ice%vec%uz = -ice%smb

        ! Calculate initial enthalpy (ice1)
        call convert_to_enthalpy(ice%vec%enth,ice%vec%T_ice,ice%vec%omega,ice%vec%T_pmp,ice%vec%cp,L_ice)
        
        ! Define temperature profile in bedrock too 
        call calc_temp_bedrock_column(ice1%vec%T_rock,ice1%kt_rock,rho_rock, &
                                    ice1%H_rock,ice1%vec%T_ice(1),ice1%Q_geo,ice1%zr%zeta)

        call convert_to_enthalpy(ice%vec%enth_rock,ice%vec%T_rock,0.0_wp,1e8_wp,ice%cp_rock,0.0_wp)

        return 

    end subroutine init_k15expb
    
    subroutine init_bg15a(ice,smb,T_srf)

        implicit none 

        type(icesheet), intent(INOUT) :: ice
        real(prec),     intent(IN)    :: smb 
        real(prec),     intent(IN)    :: T_srf ! [degrees Celcius]
            
        ! Local variables 
        integer :: k, nz 
        real(prec) :: ATT, gamma, T_init    
        real(prec), allocatable :: ux(:)
        real(prec), allocatable :: uy(:)   

        nz    = size(ice%z%zeta)

        ! Assign point values
        ice%T_srf    = T_srf + T0       ! [K]
        ice%T_shlf   = T0               ! [K] T_shlf not used in this idealized setup, set to T0  
        ice%smb      = smb              ! [m/a]
        ice%bmb      = 0.0              ! [m/a]
        ice%Q_geo    = 0.0              ! [mW/m2]
        ice%H_ice    = 200.0            ! [m] Ice thickness
        ice%H_w      = 0.0              ! [m] No basal water
        ice%Q_b      = 0.0              ! [] No basal frictional heating 
        ice%f_grnd   = 1.0              ! Grounded point 

        ice%cp_rock  = 1000.0    ! [J kg-1 K-1]
        ice%kt_rock  = 9.46e7    ! [J a-1 m-1 K-1]
        ice%H_rock   = 2000.0       ! [m] 
        ice%H_rock   = 2000.0           ! [m] 
        ice%Q_rock   = ice%Q_geo 
        
        ! EISMINT1
        ice%vec%cp      = 2009.0        ! [J kg-1 K-1]
        ice%vec%kt      = 6.6269e7      ! [J a-1 m-1 K-1]   == 2.1*sec_year  [J s-1 m-1 K-1] => J a-1 m-1 K-1]
        
        ATT             = 5.3e-24*sec_year      ! Rate factor
        gamma           = 4.0                   ! [degrees] Bed slope 
        T_init          = T0 - 1.5 

        allocate(ux(nz))
        allocate(uy(nz))

        ux = 0.0
        uy = 0.0 

        ! [J a-1 m-3] Prescribed strain heating 
        ice%vec%Q_strn  = (2.0*ATT)*(rho_ice*g*sin(gamma*pi/180.0))**4.0 &
                                * (ice%H_ice*(1.0-ice%z%zeta))**4.0
        
        ice%vec%advecxy = 0.0                                       ! [] No horizontal advection (assume constant)

        ! Write strain heating to compare basal value of ~2.6e-3 W/m-3
!         do k = nz, 1, -1 
!             write(*,*) ice%z%zeta(k), ice%vec%Q_strn(k)/sec_year 
!         end do 

        ! Calculate pressure melting point 
        ice%vec%T_pmp = calc_T_pmp(ice%H_ice,ice%z%zeta,T0,T_pmp_beta) 

        if (is_celcius) then 
            ice%T_srf     = ice%T_srf     - T0
            ice%T_shlf    = ice%T_shlf    - T0
            ice%vec%T_pmp = ice%vec%T_pmp - T0 
        end if 

        ! Define initial temperature profile
        ! (constant equal to surface temp)
        ice%vec%T_ice(nz) = T_init
        ice%vec%T_ice(1)  = T_init 

        ! Intermediate layers are linearly interpolated 
        do k = 2, nz-1 
            ice%vec%T_ice(k) = ice%vec%T_ice(1)+ice%z%zeta(k)*(ice%vec%T_ice(nz)-ice%vec%T_ice(1))
        end do 

        ! Define constant vertical velocity profile
        ice%vec%uz = -ice%smb

        ! Calculate initial enthalpy (ice1)
        call convert_to_enthalpy(ice%vec%enth,ice%vec%T_ice,ice%vec%omega,ice%vec%T_pmp,ice%vec%cp,L_ice)
        
        ! Define temperature profile in bedrock too 
        call calc_temp_bedrock_column(ice1%vec%T_rock,ice1%kt_rock,rho_rock, &
                                    ice1%H_rock,ice1%vec%T_ice(1),ice1%Q_geo,ice1%zr%zeta)

        call convert_to_enthalpy(ice%vec%enth_rock,ice%vec%T_rock,0.0_wp,1e8_wp,ice%cp_rock,0.0_wp)

        return 

    end subroutine init_bg15a
    
    subroutine icesheet_allocate(ice,nz,nzr,zeta_scale,zeta_scale_rock)
        ! Allocate the ice sheet object 

        implicit none 

        type(icesheet), intent(INOUT) :: ice 
        integer,        intent(IN)    :: nz     ! Number of ice points (aa-nodes)
        integer,        intent(IN)    :: nzr    ! Number of bedrock points (aa-nodes)
        character(*),   intent(IN)    :: zeta_scale
        character(*),   intent(IN)    :: zeta_scale_rock 

        ! Local variables 
        integer :: k, nz_ac, nzr_ac  

        ! Initialize zetas: z and zr 
        call calc_zeta(ice%z%zeta,ice%z%zeta_ac,nz_ac,nz,zeta_scale,zeta_exp=2.0_prec) 
        call calc_zeta(ice%zr%zeta,ice%zr%zeta_ac,nzr_ac,nzr,zeta_scale_rock,zeta_exp=2.0_prec) 

        ! Make sure all vectors are deallocated
        if (allocated(ice%vec%T_ice))   deallocate(ice%vec%T_ice)
        if (allocated(ice%vec%T_pmp))   deallocate(ice%vec%T_pmp)
        if (allocated(ice%vec%cp))      deallocate(ice%vec%cp)
        if (allocated(ice%vec%kt))      deallocate(ice%vec%kt)
        if (allocated(ice%vec%uz))      deallocate(ice%vec%uz)
        if (allocated(ice%vec%advecxy)) deallocate(ice%vec%advecxy)
        if (allocated(ice%vec%Q_strn))  deallocate(ice%vec%Q_strn)
        if (allocated(ice%vec%t_dep))   deallocate(ice%vec%t_dep)

        if (allocated(ice%vec%enth))    deallocate(ice%vec%enth)
        if (allocated(ice%vec%omega))   deallocate(ice%vec%omega)
        
        ! Allocate vectors with desired lengths
        allocate(ice%z%dzeta_a(nz))
        allocate(ice%z%dzeta_b(nz))
        allocate(ice%zr%dzeta_a(nzr))
        allocate(ice%zr%dzeta_b(nzr))

        allocate(ice%vec%enth(nz))
        allocate(ice%vec%T_ice(nz))
        allocate(ice%vec%omega(nz))
        allocate(ice%vec%T_pmp(nz))
        allocate(ice%vec%cp(nz))
        allocate(ice%vec%kt(nz))
        allocate(ice%vec%uz(nz_ac))
        allocate(ice%vec%advecxy(nz))
        allocate(ice%vec%Q_strn(nz))
        allocate(ice%vec%t_dep(nz))
        
        allocate(ice%vec%enth_rock(nzr))
        allocate(ice%vec%T_rock(nzr))

        ! Initialize remaining vectors to zero 
        ice%vec%enth      = 0.0 
        ice%vec%T_ice     = 0.0 
        ice%vec%omega     = 0.0 
        ice%vec%T_pmp     = 0.0 
        ice%vec%cp        = 0.0
        ice%vec%kt        = 0.0
        ice%vec%uz        = 0.0
        ice%vec%advecxy   = 0.0  
        ice%vec%Q_strn    = 0.0
        ice%vec%t_dep     = 0.0 
        
        ice%vec%enth_rock = 0.0 
        ice%vec%T_rock    = 0.0 

        ! Calculate derivative terms 
        call calc_dzeta_terms(ice%z%dzeta_a,ice%z%dzeta_b,ice%z%zeta,ice%z%zeta_ac) 
        call calc_dzeta_terms(ice%zr%dzeta_a,ice%zr%dzeta_b,ice%zr%zeta,ice%zr%zeta_ac) 

        write(*,*) "Allocated icesheet variables."

        return 

    end subroutine icesheet_allocate

    subroutine write_init(ice,filename,zeta,zeta_ac,zeta_r,time_init)

        implicit none 

        type(icesheet),   intent(IN) :: ice 
        character(len=*), intent(IN) :: filename 
        real(prec),       intent(IN) :: zeta(:)  
        real(prec),       intent(IN) :: zeta_ac(:) 
        real(prec),       intent(IN) :: zeta_r(:) 
        real(prec),       intent(IN) :: time_init

        ! Initialize netcdf file and dimensions
        call nc_create(filename)
        call nc_write_dim(filename,"zeta",    x=zeta,    units="1")
        call nc_write_dim(filename,"zeta_ac", x=zeta_ac, units="1")
        call nc_write_dim(filename,"zeta_r",  x=zeta_r,  units="1")
        call nc_write_dim(filename,"time",  x=time_init,dx=1.0_prec,nx=1,units="years",unlimited=.TRUE.)
        call nc_write_dim(filename,"pt",    x=1.0,    units="1")

        ! Write some constants 
        call nc_write(filename,"enth_cr",enth_cr,dim1="pt")

        return

    end subroutine write_init
    
    subroutine write_step(ice,vecs,filename,time,T_robin)

        implicit none 
        
        type(icesheet),         intent(IN) :: ice
        type(icesheet_vectors), intent(IN) :: vecs
        character(len=*),       intent(IN) :: filename
        real(prec),             intent(IN) :: time
        real(prec), optional,   intent(IN) :: T_robin(:) 

        ! Local variables
        integer    :: ncid, n
        real(prec) :: time_prev 
        character(len=12), parameter :: vert_dim = "zeta"

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_size(filename,"time",ncid)
        call nc_read(filename,"time",time_prev,start=[n],count=[1],ncid=ncid) 
        if (abs(time-time_prev).gt.1e-5) n = n+1 

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        ! Update variables (vectors) 
        call nc_write(filename,"enth",   vecs%enth,   units="J kg-1",   long_name="Ice enthalpy",                   dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"T_ice",  vecs%T_ice,  units="K",        long_name="Ice temperature",                dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"omega",  vecs%omega,  units="",         long_name="Ice water content (fraction)",   dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"T_pmp",  vecs%T_pmp,  units="",         long_name="Ice pressure melting point",     dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"T_prime",vecs%T_ice-vecs%T_pmp,units="K",long_name="Ice temperature",               dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        
        call nc_write(filename,"T_prime_b",vecs%T_ice(1)-vecs%T_pmp(1),units="K",long_name="Ice temperature",       dim1="time",start=[n],ncid=ncid)
        
        call nc_write(filename,"cp",     vecs%cp,     units="J kg-1 K-1",   long_name="Ice heat capacity",          dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"kt",     vecs%kt,     units="J a-1 m-1 K-1",long_name="Ice thermal conductivity",   dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"uz",     vecs%uz,     units="m a**-1",  long_name="Ice vertical velocity",          dim1="zeta_ac",dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"advecxy",vecs%advecxy,units="",         long_name="Ice horizontal advection",       dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"Q_strn", vecs%Q_strn, units="J a-1 m-3",long_name="Ice strain heating",             dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        
        call nc_write(filename,"t_dep", vecs%t_dep, units="a",long_name="Deposition time",                          dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        
        call nc_write(filename,"H_cts",  ice%H_cts, units="m",long_name="CTS height",dim1="time",start=[n],ncid=ncid)
        
        call nc_write(filename,"enth_rock", vecs%enth_rock, units="J kg-1",   long_name="Bedrock enthalpy",         dim1="zeta_r",dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"T_rock",    vecs%T_rock,    units="K",        long_name="Bedrock temperature",      dim1="zeta_r",dim2="time",start=[1,n],ncid=ncid)
        
        ! Update variables (points) 
        call nc_write(filename,"smb",     ice%smb,units="m a**-1",long_name="Surface mass balance",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"bmb",     ice%bmb,units="m a**-1",long_name="Basal mass balance",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"Q_ice_b", ice%Q_ice_b,units="mW m**-2",long_name="Ice base heat flux",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"Q_rock",  ice%Q_rock, units="mW m**-2",long_name="Bedrock surface heat flux",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"Q_b",     ice%Q_b,    units="mW m**-2",long_name="Basal frictional heating",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"Q_geo",   ice%Q_geo,  units="mW m**-2",long_name="Geothermal heat flux",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"T_srf",   ice%T_srf,units="K",long_name="Surface temperature",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"H_ice",   ice%H_ice,units="m",long_name="Ice thickness",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"H_w",     ice%H_w,units="m",long_name="Basal water thickness",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"f_grnd",  ice%f_grnd,units="1",long_name="Grounded fraction",dim1="time",start=[n],ncid=ncid)
        
        ! If available, compare with Robin analytical solution 
        if (present(T_robin)) then
            call nc_write(filename,"T_robin",  T_robin,  units="K", &
                            long_name="Ice temperature from Robin solution", &
                            dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
         
            call nc_write(filename,"T_diff",  vecs%T_ice-T_robin,  units="K", &
                            long_name="Ice temperature difference with Robin solution", &
                            dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        end if 

        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine write_step

    subroutine calc_zeta(zeta_aa,zeta_ac,nz_ac,nz_aa,zeta_scale,zeta_exp)
        ! Calculate the vertical axis cell-centers first (aa-nodes),
        ! including a cell centered at the base on the surface. 
        ! Then calculate the vertical axis cell-edges (ac-nodes).
        ! The base and surface ac-nodes coincide with the cell centers.
        ! There is one more cell-edge than cell-centers (nz_ac=nz_aa+1)

        implicit none 

        real(prec), allocatable, intent(INOUT) :: zeta_aa(:) 
        real(prec), allocatable, intent(INOUT) :: zeta_ac(:) 
        integer,                 intent(OUT)   :: nz_ac 
        integer,      intent(IN)   :: nz_aa 
        character(*), intent(IN)   :: zeta_scale 
        real(prec),   intent(IN)   :: zeta_exp 

        ! Local variables
        integer :: k 
        real(prec), allocatable :: tmp(:) 

        ! Define size of zeta ac-node vector
        nz_ac = nz_aa + 1 

        ! First allocate arrays 
        if (allocated(zeta_aa)) deallocate(zeta_aa)
        if (allocated(zeta_ac)) deallocate(zeta_ac)
        allocate(zeta_aa(nz_aa))
        allocate(zeta_ac(nz_ac))

        ! Initially define a linear zeta scale 
        ! Base = 0.0, Surface = 1.0 
        do k = 1, nz_aa
            zeta_aa(k) = 0.0 + 1.0*(k-1)/real(nz_aa-1)
        end do 

        ! Scale zeta to produce different resolution through column if desired
        ! zeta_scale = ["linear","exp","wave"]
        select case(trim(zeta_scale))
            
            case("exp")
                ! Increase resolution at the base 
                
                zeta_aa = zeta_aa**(zeta_exp) 

            case("exp-inv")
                ! Increase resolution at the surface 
                
                zeta_aa = 1.0_wp - zeta_aa**(zeta_exp)

                ! Reverse order 
                allocate(tmp(nz_aa))
                tmp = zeta_aa 
                do k = 1, nz_aa
                    zeta_aa(k) = tmp(nz_aa-k+1)
                end do 

            case("tanh")
                ! Increase resolution at base and surface 

                zeta_aa = tanh(1.0*pi*(zeta_aa-0.5))
                zeta_aa = zeta_aa - minval(zeta_aa)
                zeta_aa = zeta_aa / maxval(zeta_aa)

            case DEFAULT
            ! Do nothing, scale should be linear as defined above
        
        end select  
        
        ! Get zeta_ac (between zeta_aa values, as well as at the base and surface)
        zeta_ac(1) = 0.0 
        do k = 2, nz_ac-1
            zeta_ac(k) = 0.5 * (zeta_aa(k-1)+zeta_aa(k))
        end do 
        zeta_ac(nz_ac) = 1.0 

        ! =================
        ! write(*,*) "Vertical axis:"
        ! do k = nz_ac, 1, -1
        !     if (k .eq. nz_ac) then 
        !         write(*,*) k, -9999.0, zeta_ac(k)
        !     else 
        !         write(*,*) k, zeta_aa(k), zeta_ac(k) 
        !     end if  
        ! end do 
        ! stop 
        ! =================

        return 

    end subroutine calc_zeta
    
end program test_icetemp
