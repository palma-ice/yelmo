program test_icetemp 

    use ncio  
    
    use yelmo_defs
    use thermodynamics 
    use ice_enthalpy  
    use ice_age 

    use interp1D 

    implicit none 

    type icesheet_vectors
        real(prec), allocatable :: zeta(:)    ! [-] Sigma coordinates from 0:1 (either height or depth)
        real(prec), allocatable :: zeta_ac(:) ! nz-1 [-] sigma coordinates for internal ice layer edges
        real(prec), allocatable :: T_ice(:)   ! [K] Ice temperature 
        real(prec), allocatable :: T_pmp(:)   ! [K] Ice pressure melting point 
        real(prec), allocatable :: cp(:)      ! [] Ice heat capacity 
        real(prec), allocatable :: kt(:)      ! [] Ice conductivity  
        real(prec), allocatable :: uz(:)      ! [m a-1] Vertical velocity 
        real(prec), allocatable :: advecxy(:) ! [] Horizontal heat advection magnitude
        real(prec), allocatable :: Q_strn(:)  ! [K a-1] Strain heating 
        real(prec), allocatable :: t_dep(:)   ! [a] Deposition time 
        real(prec), allocatable :: enth(:)    ! [J kg-1] Ice enthalpy
        real(prec), allocatable :: omega(:)   ! [-] Ice water content (fraction)
    end type 

    type poly_state_class 

        integer :: nz_pt, nz_pc
        integer :: nz_aa, nz_ac 
        real(prec) :: H_cts

        real(prec), allocatable :: zeta_pt(:)       ! zeta_aa for polythermal temperate (pt) zone only 
        real(prec), allocatable :: zeta_pc(:)       ! zeta_aa for polythermal cold (pc) zone only 
        
        real(prec), allocatable :: zeta_aa(:)   ! Layer centers (aa-nodes), plus base and surface: nz_aa points 
        real(prec), allocatable :: zeta_ac(:)   ! Layer borders (ac-nodes), plus base and surface: nz_ac == nz_aa-1 points

        real(prec), allocatable :: enth(:)      ! [J m-3] Ice enthalpy 
        real(prec), allocatable :: T_ice(:)     ! [K]     Ice temp. 
        real(prec), allocatable :: omega(:)     ! [--]    Ice water content
        real(prec), allocatable :: T_pmp(:)     ! Pressure-corrected melting point
        
        real(prec), allocatable :: cp(:)        ! Specific heat capacity  
        real(prec), allocatable :: kt(:)        ! Heat conductivity  

        real(prec), allocatable :: advecxy(:) 
        real(prec), allocatable :: uz(:) 
        real(prec), allocatable :: Q_strn(:)    ! Internal heat production 
        

    end type 

    type icesheet 
        real(prec) :: H_ice             ! [m] Ice thickness 
        real(prec) :: H_w               ! [m] Water present at the ice base 
        real(prec) :: T_srf             ! [K] Ice surface temperature 
        real(prec) :: T_shlf            ! [K] Ice shelf base temperature 
        real(prec) :: smb               ! [m a**-1] Surface mass balance
        real(prec) :: bmb               ! [m a**-1] Basal mass balance
        real(prec) :: Q_ice_b           ! [J a-1 m-2] Ice basal heat flux (positive up)
        real(prec) :: H_cts             ! [m] cold-temperate transition surface (CTS) height
        real(prec) :: Q_geo             ! [mW m-2] Geothermal heat flux 
        real(prec) :: Q_b               ! [J a-1 m-2] Basal heat production
        real(prec) :: f_grnd            ! [-] Grounded fraction 
        character(len=56) :: age_method ! Method to use for age calculation 
        real(prec) :: age_impl_kappa    ! [m2 a-1] Artificial diffusion term for implicit age solving 
        type(icesheet_vectors) :: vec   ! For height coordinate systems with k=1 base and k=nz surface
        type(poly_state_class) :: poly  ! For two-layered calculations
    end type 

    ! Define different icesheet objects for use in proram
    type(icesheet) :: ice1
    type(icesheet) :: robin 
    type(icesheet) :: diff 
    
    ! Local variables
    real(prec)         :: t_start, t_end, dt, time  
    integer            :: n, ntot  
    character(len=512) :: file1D 
    real(prec)         :: dt_out 
    character(len=56)  :: zeta_scale 
    integer            :: nz 
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

    logical, parameter :: testing_poly = .TRUE.
    real(prec)         :: H_cts_prev, H_cts_ref, H_cts_now, E0, E1, dEdz 
    integer            :: k_cts, k 

    type(poly_state_class) :: poly  ! For two-layered calculations
    integer :: iter, n_iter 

    ! General initialization of yelmo constants (used globally)
    call yelmo_global_init("par/yelmo_const_EISMINT.nml")
    
    ! ===============================================================
    ! User options 

    experiment     = "bg15a"        ! "eismint", "k15expa", "k15expb", "bg15a"
    
    ! General options
    zeta_scale      = "linear"      ! "linear", "exp", "tanh"
    nz              = 22            ! [--] Number of ice sheet points (aa-nodes + base + surface)
    is_celcius      = .FALSE. 

    age_method      = "expl"        ! "expl" or "impl"
    age_impl_kappa  = 1.5           ! [m2 a-1] Artificial diffusion for age tracing

    enth_solver     = "enth"        ! "enth" or "temp" 
    omega_max       = 0.03          ! Maximum allowed water content (fraction)
    enth_cr         = 1e-3          ! Enthalpy solver: conductivity ratio kappa_water / kappa_ice 

    !file1D          = "output/test_"//trim(experiment)//"_nz30_sp.nc" 
    file1D          = "output/test_"//trim(experiment)//".nc"

    ! Overwrite options for nz and enth_cr if available from arguments
    narg = command_argument_count() 
    if (narg .gt. 0) then 
        call get_command_argument(1,arg_nz)
        call get_command_argument(2,arg_cr)

        read(arg_nz,*)  nz
        read(arg_cr,*)  enth_cr
    end if 

    if (trim(experiment) .eq. "k15expa" .or. trim(experiment) .eq. "k15expb") then 
        ! Use a more precise filename to specify cr value and dz
        write(file1D,"(a,e8.2,a,e8.2,a)") "output/test_"//trim(experiment)//"_dz", (200.0/(nz-2)), "_cr", enth_cr, ".nc"
    end if 

    ! ===============================================================

    T0_ref = T0 
    if (is_celcius) T0_ref = 0.0 

    ! Initialize icesheet object 
    call icesheet_allocate(ice1,nz=nz,zeta_scale=zeta_scale) 

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

            call init_eismint_summit(ice1,smb=0.1_prec)

    end select 

    ! Initialize polythermal data structure too 
    !call poly_init(ice1%poly,nz_pt=11,nz_pc=392,zeta_scale=zeta_scale,zeta_exp=2.0_prec)
    call poly_init(ice1%poly,nz_pt=11,nz_pc=32,zeta_scale=zeta_scale,zeta_exp=2.0_prec)

if (testing_poly) then

    ! Calculate the poly vertical axis at each grid point
    call calc_zeta_combined(ice1%poly%zeta_aa,ice1%poly%zeta_ac,ice1%poly%zeta_pt,ice1%poly%zeta_pc, &
                                                                        ice1%H_cts,ice1%H_ice)
    
    ice1%poly%cp    = interp_linear(ice1%vec%zeta,ice1%vec%cp,ice1%poly%zeta_aa)
    ice1%poly%kt    = interp_linear(ice1%vec%zeta,ice1%vec%kt,ice1%poly%zeta_aa)
    
    ice1%poly%enth  = interp_linear(ice1%vec%zeta,ice1%vec%enth,ice1%poly%zeta_aa)
    ice1%poly%T_ice = interp_linear(ice1%vec%zeta,ice1%vec%T_ice,ice1%poly%zeta_aa)
    ice1%poly%omega = interp_linear(ice1%vec%zeta,ice1%vec%omega,ice1%poly%zeta_aa)
    ice1%poly%T_pmp = interp_linear(ice1%vec%zeta,ice1%vec%T_pmp,ice1%poly%zeta_aa)

    call update_poly(ice1%poly,ice1%vec%advecxy,ice1%vec%Q_strn,ice1%vec%uz,ice1%vec%zeta, &
                                                    ice1%vec%zeta_ac,ice1%H_cts,ice1%H_ice)

else

    ! Simply set them equal for now
    ice1%poly%zeta_aa = ice1%vec%zeta
    ice1%poly%zeta_ac = ice1%vec%zeta_ac
    
    ice1%poly%cp      = ice1%vec%cp 
    ice1%poly%kt      = ice1%vec%kt 
    
    ice1%poly%enth    = ice1%vec%enth  
    ice1%poly%T_ice   = ice1%vec%T_ice 
    ice1%poly%omega   = ice1%vec%omega 
    ice1%poly%T_pmp   = ice1%vec%T_pmp 
    
    ice1%poly%advecxy = ice1%vec%advecxy 
    ice1%poly%Q_strn  = ice1%vec%Q_strn 
    ice1%poly%uz      = ice1%vec%uz 

end if 

    ! Initialize time and calculate number of time steps to iterate and 
    time = t_start 
    ntot = (t_end-t_start)/dt 
    
    ! Set age method and kappa 
    ice1%age_method     = age_method 
    ice1%age_impl_kappa = age_impl_kappa

    ! Calculate the robin solution for comparison 
    robin = ice1  
    robin%vec%T_ice = calc_temp_robin_column(robin%vec%zeta,robin%vec%T_pmp,robin%vec%kt,robin%vec%cp,rho_ice, &
                                       robin%H_ice,robin%T_srf,robin%smb,robin%Q_geo,is_float=robin%f_grnd.eq.0.0)

    ! Calculate initial enthalpy (ice1)
    call convert_to_enthalpy(ice1%vec%enth,ice1%vec%T_ice,ice1%vec%omega,ice1%vec%T_pmp,ice1%vec%cp,L_ice)
    call convert_to_enthalpy(ice1%poly%enth,ice1%poly%T_ice,ice1%poly%omega,ice1%poly%T_pmp,ice1%poly%cp,L_ice)

    ! Calculate initial enthalpy (robin)
    call convert_to_enthalpy(robin%vec%enth,robin%vec%T_ice,robin%vec%omega,robin%vec%T_pmp,robin%vec%cp,L_ice)

!     ! Write Robin solution 
!     file1D = "robin.nc"
!     call write_init(robin,filename=file1D,zeta=robin%vec%zeta,zeta_ac=robin%vec%zeta_ac,time_init=time)
!     call write_step(robin,robin%vec,filename=file1D,time=time)

    ! Initialize output file for model and write intial conditions 
    call write_init(ice1,filename=file1D,zeta=ice1%vec%zeta,zeta_ac=ice1%vec%zeta_ac,zeta_pt=ice1%poly%zeta_pt, &
                        zeta_pc=ice1%poly%zeta_pc,time_init=time)
    call write_step(ice1,ice1%vec,ice1%poly,filename=file1D,time=time)

    ! Ensure zero basal water thickness to start 
    ice1%H_w = 0.0 

    ! Assume H_cts is also zero to start
    ice1%H_cts = 0.0 
    H_cts_prev = ice1%H_cts

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

if (testing_poly) then

    ! Store H_cts from previous timestep 
    H_cts_prev = ice1%H_cts 

    ! Find height of CTS index - highest temperate layer 
    k_cts     = get_cts_index(ice1%poly%enth,ice1%poly%T_pmp*ice1%poly%cp)
    H_cts_ref = ice1%poly%zeta_aa(k_cts)*ice1%H_ice 

if (.TRUE.) then 
    ! Perform iterations to improve guess of H_cts 

!         do iter = 1, 10
            
!             ! Reset to original configuration
!             poly = ice1%poly 

!             !H_cts_now = ice1%H_cts
!             H_cts_now = H_cts_ref + (iter)*0.5 
             
!             call update_poly(poly,ice1%vec%advecxy,ice1%vec%Q_strn,ice1%vec%uz,ice1%vec%zeta, &
!                                                                 ice1%vec%zeta_ac,ice1%H_cts,ice1%H_ice)
            
!             call calc_enth_column(poly%enth,poly%T_ice,poly%omega,ice1%bmb,ice1%Q_ice_b,ice1%H_cts,poly%T_pmp, &
!                     poly%cp,poly%kt,poly%advecxy,poly%uz,poly%Q_strn,ice1%Q_b,ice1%Q_geo,ice1%T_srf,ice1%T_shlf, &
!                     ice1%H_ice,ice1%H_w,ice1%f_grnd,poly%zeta_aa,poly%zeta_ac, &
!                     enth_cr,omega_max,T0_ref,dt)

!             ! Find height of CTS index - highest temperate layer 
!             k_cts = get_cts_index(poly%enth,poly%T_pmp*poly%cp)
            
!             dEdz = (poly%enth(k_cts+2)-poly%enth(k_cts+1)) / (poly%zeta_aa(k_cts+2)-poly%zeta_aa(k_cts+1))
!             E0   = poly%enth(k_cts)-poly%T_pmp(k_cts)*poly%cp(k_cts)
!             E1   = poly%enth(k_cts+1)-poly%T_pmp(k_cts+1)*poly%cp(k_cts+1)

!             write(*,"(i10,3f10.4,3g12.5)") iter, H_cts_prev, H_cts_now, ice1%H_cts, dEdz, E0, E1 
            
!         end do  

        ice1%H_cts = H_cts_prev
        poly = ice1%poly 
            
        do iter = 1, 5
            
            H_cts_now = ice1%H_cts
            !H_cts_now = H_cts_prev + (iter)*0.5   
            call update_poly(poly,ice1%vec%advecxy,ice1%vec%Q_strn,ice1%vec%uz,ice1%vec%zeta, &
                                                                ice1%vec%zeta_ac,H_cts_now,ice1%H_ice)
            
            call calc_enth_column(poly%enth,poly%T_ice,poly%omega,ice1%bmb,ice1%Q_ice_b,ice1%H_cts,poly%T_pmp, &
                    poly%cp,poly%kt,poly%advecxy,poly%uz,poly%Q_strn,ice1%Q_b,ice1%Q_geo,ice1%T_srf,ice1%T_shlf, &
                    ice1%H_ice,ice1%H_w,ice1%f_grnd,poly%zeta_aa,poly%zeta_ac, &
                    enth_cr,omega_max,T0_ref,dt)

        end do 

else 
        H_cts_now = H_cts_prev 
        !H_cts_now = ice1%H_cts

end if 
    
        
        call update_poly(ice1%poly,ice1%vec%advecxy,ice1%vec%Q_strn,ice1%vec%uz,ice1%vec%zeta, &
                                                            ice1%vec%zeta_ac,H_cts_now,ice1%H_ice)

        H_cts_prev = ice1%H_cts 

end if 

        call calc_enth_column(ice1%poly%enth,ice1%poly%T_ice,ice1%poly%omega,ice1%bmb,ice1%Q_ice_b,ice1%H_cts,ice1%poly%T_pmp, &
                ice1%poly%cp,ice1%poly%kt,ice1%poly%advecxy,ice1%poly%uz,ice1%poly%Q_strn,ice1%Q_b,ice1%Q_geo,ice1%T_srf,ice1%T_shlf, &
                ice1%H_ice,ice1%H_w,ice1%f_grnd,ice1%poly%zeta_aa,ice1%poly%zeta_ac, &
                enth_cr,omega_max,T0_ref,dt)

!         call calc_enth_column(ice1%vec%enth,ice1%vec%T_ice,ice1%vec%omega,ice1%bmb,ice1%Q_ice_b,ice1%H_cts,ice1%vec%T_pmp, &
!                 ice1%vec%cp,ice1%vec%kt,ice1%vec%advecxy,ice1%vec%uz,ice1%vec%Q_strn,ice1%Q_b,ice1%Q_geo,ice1%T_srf,ice1%T_shlf, &
!                 ice1%H_ice,ice1%H_w,ice1%f_grnd,ice1%vec%zeta,ice1%vec%zeta_ac, &
!                 enth_cr,omega_max,T0_ref,dt)

if (testing_poly) then 

        call update_enth_1layer(ice1%vec%enth,ice1%vec%T_ice,ice1%vec%omega,ice1%vec%T_pmp,ice1%vec%cp, &
                                                ice1%vec%zeta,ice1%vec%zeta_ac,ice1%poly,L_ice)
else

        ice1%vec%enth  = ice1%poly%enth 
        ice1%vec%T_ice = ice1%poly%T_ice 
        ice1%vec%omega = ice1%poly%omega 

end if 

        ! Update basal water thickness [m/a i.e.] => [m/a w.e.]
        ice1%H_w = ice1%H_w - (ice1%bmb*rho_ice/rho_w)*dt 

        if (trim(age_method) .eq. "impl") then 
            call calc_tracer_column(ice1%vec%t_dep,ice1%vec%uz,ice1%vec%advecxy*0.0,time,ice1%bmb, &
                                    ice1%H_ice,ice1%vec%zeta,ice1%vec%zeta_ac, &
                                    ice1%age_impl_kappa,dt)
        else 
            call calc_tracer_column_expl(ice1%vec%t_dep,ice1%vec%uz,ice1%vec%advecxy*0.0,time,ice1%bmb,ice1%H_ice,ice1%vec%zeta,ice1%vec%zeta_ac,dt)
        end if 

        if (mod(time,dt_out)==0) then 
            call write_step(ice1,ice1%vec,ice1%poly,filename=file1D,time=time,T_robin=robin%vec%T_ice)
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
    write(*,*) "========================="
    write(*,*)

contains 

    subroutine update_poly(poly,advecxy,Q_strn,uz,zeta_aa,zeta_ac,H_cts,H_ice)

        implicit none

        type(poly_state_class), intent(INOUT) :: poly 
        real(prec), intent(IN) :: advecxy(:) 
        real(prec), intent(IN) :: Q_strn(:) 
        real(prec), intent(IN) :: uz(:) 
        real(prec), intent(IN) :: zeta_aa(:) 
        real(prec), intent(IN) :: zeta_ac(:) 
        real(prec), intent(IN) :: H_cts 
        real(prec), intent(IN) :: H_ice

        ! Local variables 
        integer    :: k, k_cts, nz  
        real(prec) :: H_cts_now 

        real(prec), allocatable :: p_zeta0(:) 
        real(prec), allocatable :: p_enth0(:) 
        real(prec), allocatable :: p_enth_pmp0(:) 
        
        allocate(p_zeta0(size(poly%enth,1)))
        allocate(p_enth0(size(poly%enth,1)))
        allocate(p_enth_pmp0(size(poly%enth,1)))
        
        H_cts_now = H_cts
        !H_cts_now = floor(1e1*H_cts)*1e-1
        !H_cts_now = 32.0
        
        H_cts_now = max(H_cts_now,1.0_prec)
        
        ! Update poly zeta axis 
        call calc_zeta_combined(poly%zeta_aa,poly%zeta_ac,poly%zeta_pt,poly%zeta_pc,H_cts_now,H_ice)

        ! Store original enth value and axis
        p_zeta0     = poly%zeta_aa  
        p_enth0     = poly%enth 
        p_enth_pmp0 = poly%T_pmp * poly%cp 
        
        ! Find height of CTS index - highest temperate layer 
        k_cts = get_cts_index(poly%enth,poly%T_pmp*poly%cp)
        
        ! Update enth 
        !poly%enth = interp_linear(p_zeta0,p_enth0,poly%zeta_aa)
        !call interp_enth_column(poly%enth,poly%zeta_aa,p_enth0,p_zeta0,p_enth_pmp0,H_cts_now,H_ice)
        call interp1D_bins(poly%enth,poly%zeta_aa,p_enth0,p_zeta0)

        ! Update external variables 
        poly%advecxy = interp_linear(zeta_aa,advecxy,poly%zeta_aa)
        poly%Q_strn  = interp_linear(zeta_aa,Q_strn,poly%zeta_aa)
        poly%uz      = interp_linear(zeta_ac,uz,poly%zeta_ac)

!         write(*,*) "zeta_ac:   ", minval(zeta_ac), maxval(zeta_ac)
!         write(*,*) "p_zeta_ac: ", minval(poly%zeta_ac), maxval(poly%zeta_ac)
!         write(*,*) "uz:        ", minval(uz), maxval(uz)
!         write(*,*) "p_uz:      ", minval(poly%uz), maxval(poly%uz)
        
!         stop 

        return 

    end subroutine update_poly

        subroutine update_enth_1layer(enth,T_ice,omega,T_pmp,cp,zeta_aa,zeta_ac,poly,L_ice)

        implicit none

        real(prec), intent(INOUT) :: enth(:) 
        real(prec), intent(INOUT) :: T_ice(:) 
        real(prec), intent(INOUT) :: omega(:) 
        real(prec), intent(INOUT) :: T_pmp(:)
        real(prec), intent(INOUT) :: cp(:)
        real(prec), intent(IN) :: zeta_aa(:) 
        real(prec), intent(IN) :: zeta_ac(:) 
        type(poly_state_class), intent(IN) :: poly 
        real(prec), intent(INOUT) :: L_ice
        
        !enth  = interp_linear(poly%zeta_aa,poly%enth,zeta_aa)
        call interp1D_bins(enth,zeta_aa,poly%enth,poly%zeta_aa)

        ! Get temperature and water content too
        call convert_from_enthalpy_column(enth,T_ice,omega,T_pmp,cp,L_ice)
         
        return 

    end subroutine update_enth_1layer


    subroutine init_eismint_summit(ice,smb)

        implicit none 

        type(icesheet), intent(INOUT) :: ice
        real(prec),     intent(IN)    :: smb 

        ! Local variables 
        integer :: k, nz, nz_ac   

        nz    = size(ice%vec%zeta)
        nz_ac = nz - 1 

        ! Assign point values
        ice%T_srf    = 239.0       ! [K]
        ice%T_shlf   = T0          ! [K] T_shlf not used in this idealized setup, set to T0  
        ice%smb      = smb         ! [m/a]
        ice%bmb      = 0.0         ! [m/a]
        ice%Q_geo    = 42.0        ! [mW/m2]
        ice%H_ice    = 2997.0      ! [m] Summit thickness
        ice%H_w      = 0.0         ! [m] No basal water
        ice%Q_b      = 0.0         ! [] No basal frictional heating 
        ice%f_grnd   = 1.0         ! Grounded point 

        ! EISMINT1
        ice%vec%cp      = 2009.0    ! [J kg-1 K-1]
        ice%vec%kt      = 6.67e7    ! [J a-1 m-1 K-1]
        
        ice%vec%Q_strn  = 0.0       ! [J a-1 m-3] No internal strain heating 
        ice%vec%advecxy = 0.0       ! [] No horizontal advection 

        ! Calculate pressure melting point 
        ice%vec%T_pmp = calc_T_pmp(ice%H_ice,ice%vec%zeta,T0,T_pmp_beta) 

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
            ice%vec%T_ice(k) = ice%vec%T_ice(1)+ice%vec%zeta(k)*(ice%vec%T_ice(nz)-ice%vec%T_ice(1))
        end do 

        ! Define linear vertical velocity profile
        ice%vec%uz = -ice%smb*ice%vec%zeta_ac 
        
        return 

    end subroutine init_eismint_summit 

    subroutine init_k15expa(ice)

        implicit none 

        type(icesheet), intent(INOUT) :: ice

        ! Local variables 
        integer :: k, nz 

        nz    = size(ice%vec%zeta)

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

        ! EISMINT1
        ice%vec%cp      = 2009.0        ! [J kg-1 K-1]
        ice%vec%kt      = 6.6269e7      ! [J a-1 m-1 K-1]   == 2.1*sec_year  [J s-1 m-1 K-1] => J a-1 m-1 K-1]
        
        ice%vec%Q_strn  = 0.0           ! [] No internal strain heating 
        ice%vec%advecxy = 0.0           ! [] No horizontal advection 

        ! Calculate pressure melting point 
        ice%vec%T_pmp = calc_T_pmp(ice%H_ice,ice%vec%zeta,T0,T_pmp_beta) 

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
            ice%vec%T_ice(k) = ice%vec%T_ice(1)+ice%vec%zeta(k)*(ice%vec%T_ice(nz)-ice%vec%T_ice(1))
        end do 

        ! Define linear vertical velocity profile
        ice%vec%uz = -ice%smb*ice%vec%zeta_ac 

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

        nz    = size(ice%vec%zeta)

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

        ux = 0.5*ATT*(rho_ice*g*sin(gamma*pi/180.0))**3 * (ice%H_ice**4 - (ice%H_ice - ice%H_ice*ice%vec%zeta)**4)
        uy = 0.0 

        eps = ATT*(rho_ice*g*sin(gamma*pi/180.0))**3 *(ice%H_ice - ice%H_ice*ice%vec%zeta)**3
        mu  = 0.5*ATT**(-1.0/3.0)*eps**(-2.0/3.0)

        ice%vec%Q_strn(nz)  = 0.0  
        ice%vec%Q_strn(1:nz-1)  = 4.0*mu(1:nz-1)*eps(1:nz-1)**2     ! [J a-1 m-3] Prescribed strain heating 
        ice%vec%advecxy = 0.0                                       ! [] No horizontal advection (assume constant)

        ! Write strain heating to compare basal value of ~2.6e-3 W/m-3
        do k = nz, 1, -1 
            write(*,*) ice%vec%zeta(k), ice%vec%Q_strn(k)/sec_year 
        end do 

        ! Calculate pressure melting point 
        ice%vec%T_pmp = calc_T_pmp(ice%H_ice,ice%vec%zeta,T0,T_pmp_beta) 

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
            ice%vec%T_ice(k) = ice%vec%T_ice(1)+ice%vec%zeta(k)*(ice%vec%T_ice(nz)-ice%vec%T_ice(1))
        end do 

        ! Define constant vertical velocity profile
        ice%vec%uz = -ice%smb

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

        nz    = size(ice%vec%zeta)

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
                                * (ice%H_ice*(1.0-ice%vec%zeta))**4.0
        
        ice%vec%advecxy = 0.0                                       ! [] No horizontal advection (assume constant)

        ! Write strain heating to compare basal value of ~2.6e-3 W/m-3
!         do k = nz, 1, -1 
!             write(*,*) ice%vec%zeta(k), ice%vec%Q_strn(k)/sec_year 
!         end do 

        ! Calculate pressure melting point 
        ice%vec%T_pmp = calc_T_pmp(ice%H_ice,ice%vec%zeta,T0,T_pmp_beta) 

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
            ice%vec%T_ice(k) = ice%vec%T_ice(1)+ice%vec%zeta(k)*(ice%vec%T_ice(nz)-ice%vec%T_ice(1))
        end do 

        ! Define constant vertical velocity profile
        ice%vec%uz = -ice%smb

        return 

    end subroutine init_bg15a 
    
    subroutine icesheet_allocate(ice,nz,zeta_scale)
        ! Allocate the ice sheet object 

        implicit none 

        type(icesheet), intent(INOUT) :: ice 
        integer,        intent(IN)    :: nz     ! Number of ice points (aa-nodes)
        character(*),   intent(IN)    :: zeta_scale

        ! Local variables 
        integer :: k, nz_ac 

        nz_ac = nz-1 

        ! First allocate 'up' variables (with vertical coordinate as height)

        ! Make sure all vectors are deallocated
        if (allocated(ice%vec%zeta))     deallocate(ice%vec%zeta)
        if (allocated(ice%vec%zeta_ac))  deallocate(ice%vec%zeta_ac)
        
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
        allocate(ice%vec%zeta(nz))
        allocate(ice%vec%zeta_ac(nz_ac))
        
        allocate(ice%vec%T_ice(nz))
        allocate(ice%vec%T_pmp(nz))
        allocate(ice%vec%cp(nz))
        allocate(ice%vec%kt(nz))
        allocate(ice%vec%uz(nz_ac))
        allocate(ice%vec%advecxy(nz))
        allocate(ice%vec%Q_strn(nz))
        allocate(ice%vec%t_dep(nz))
        allocate(ice%vec%enth(nz))
        allocate(ice%vec%omega(nz))

        ! Initialize zeta 
        call calc_zeta(ice%vec%zeta,ice%vec%zeta_ac,zeta_scale,zeta_exp=2.0_prec) 

        ! Initialize remaining vectors to zero 
        ice%vec%T_ice   = 0.0 
        ice%vec%T_pmp   = 0.0 
        ice%vec%cp      = 0.0
        ice%vec%kt      = 0.0
        ice%vec%uz      = 0.0
        ice%vec%advecxy = 0.0  
        ice%vec%Q_strn  = 0.0
        ice%vec%t_dep   = 0.0 
        ice%vec%enth    = 0.0 
        ice%vec%omega   = 0.0 

        write(*,*) "Allocated icesheet variables."

        return 

    end subroutine icesheet_allocate 

    subroutine poly_init(poly,nz_pt,nz_pc,zeta_scale,zeta_exp)

        implicit none 

        type(poly_state_class), intent(INOUT) :: poly 
        integer,      intent(IN) :: nz_pt
        integer,      intent(IN) :: nz_pc  
        character(*), intent(IN) :: zeta_scale 
        real(prec),   intent(IN) :: zeta_exp 

        ! Local variables 
        integer    :: k  

        poly%nz_pt  = nz_pt
        poly%nz_pc  = nz_pc
        poly%nz_aa  = poly%nz_pt + (poly%nz_pc-1) 
        poly%nz_ac  = poly%nz_aa - 1  

        ! 1D axis vectors (separate temperate and cold axes)
        allocate(poly%zeta_pt(poly%nz_pt)) 
        allocate(poly%zeta_pc(poly%nz_pc)) 

        ! 3D axis arrays (combined polythermal axis, different for each column)
        allocate(poly%zeta_aa(poly%nz_aa)) 
        allocate(poly%zeta_ac(poly%nz_ac)) 
        
        ! Variables 
        allocate(poly%enth(poly%nz_aa))
        allocate(poly%T_ice(poly%nz_aa))
        allocate(poly%omega(poly%nz_aa))
        allocate(poly%T_pmp(poly%nz_aa))
        allocate(poly%cp(poly%nz_aa))
        allocate(poly%kt(poly%nz_aa))
        
        allocate(poly%advecxy(poly%nz_aa))
        allocate(poly%Q_strn(poly%nz_aa))
        allocate(poly%uz(poly%nz_ac))
        
        ! Calculate the temperate and cold vertical axes 
        call calc_zeta_twolayers(poly%zeta_pt,poly%zeta_pc,zeta_scale,zeta_exp)


        ! Test routine to make combined axis::
!         call calc_zeta_combined(poly%zeta_aa,poly%zeta_ac,poly%zeta_pt,poly%zeta_pc,H_cts=20.0,H_ice=200.0)

        return 

    end subroutine poly_init 
    
    subroutine write_init(ice,filename,zeta,zeta_ac,zeta_pt,zeta_pc,time_init)

        implicit none 

        type(icesheet),   intent(IN) :: ice 
        character(len=*), intent(IN) :: filename 
        real(prec),       intent(IN) :: zeta(:)  
        real(prec),       intent(IN) :: zeta_ac(:) 
        real(prec),       intent(IN) :: zeta_pt(:)
        real(prec),       intent(IN) :: zeta_pc(:)  
        real(prec),       intent(IN) :: time_init

        ! Local variables
        integer :: npt_poly 

        ! Initialize netcdf file and dimensions
        call nc_create(filename)
        call nc_write_dim(filename,"zeta",    x=zeta,    units="1")
        call nc_write_dim(filename,"zeta_ac", x=zeta_ac, units="1")
        call nc_write_dim(filename,"zeta_pt", x=zeta_pt, units="1")
        call nc_write_dim(filename,"zeta_pc", x=zeta_pc, units="1")
        call nc_write_dim(filename,"time",  x=time_init,dx=1.0_prec,nx=1,units="years",unlimited=.TRUE.)
        call nc_write_dim(filename,"pt",    x=1.0,    units="1")

        ! Write the number of poly points 
        npt_poly = size(zeta_pt,1) + (size(zeta_pc,1)-1) 
        call nc_write_dim(filename,"zeta_px_aa",x=1,nx=npt_poly,dx=1,units="1")
        call nc_write_dim(filename,"zeta_px_ac",x=1,nx=npt_poly-1,dx=1,units="1")

        ! Write some constants 
        call nc_write(filename,"enth_cr",enth_cr,dim1="pt")

        return

    end subroutine write_init 
    
    subroutine write_step(ice,vecs,poly,filename,time,T_robin)

        implicit none 
        
        type(icesheet),         intent(IN) :: ice
        type(icesheet_vectors), intent(IN) :: vecs
        type(poly_state_class), intent(IN) :: poly
        character(len=*),       intent(IN) :: filename
        real(prec),             intent(IN) :: time
        real(prec), optional,   intent(IN) :: T_robin(:) 

        ! Local variables
        integer    :: ncid, n
        real(prec) :: time_prev 
        character(len=12), parameter :: vert_dim = "zeta"
        character(len=12), parameter :: vert_dim_poly = "zeta_px_aa"

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_size(filename,"time",ncid)
        call nc_read(filename,"time",time_prev,start=[n],count=[1],ncid=ncid) 
        if (abs(time-time_prev).gt.1e-5) n = n+1 

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        ! Update poly variables 
        call nc_write(filename,"pp_zeta",   poly%zeta_aa,units="J kg-1",    long_name="Vertical axis",                dim1=vert_dim_poly,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"pp_enth",   poly%enth,   units="J kg-1",    long_name="Ice enthalpy",                 dim1=vert_dim_poly,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"pp_T_ice",  poly%T_ice,  units="K",         long_name="Ice temperature",              dim1=vert_dim_poly,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"pp_omega",  poly%omega,  units="",          long_name="Ice water content (fraction)", dim1=vert_dim_poly,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"pp_T_pmp",  poly%T_pmp,  units="",          long_name="Ice pressure melting point",   dim1=vert_dim_poly,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"pp_T_prime",poly%T_ice-poly%T_pmp,units="K",long_name="Ice temperature",              dim1=vert_dim_poly,dim2="time",start=[1,n],ncid=ncid)
        
        call nc_write(filename,"pp_cp",     poly%cp,     units="J kg-1 K-1",   long_name="Ice heat capacity",       dim1=vert_dim_poly,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"pp_kt",     poly%kt,     units="J a-1 m-1 K-1",long_name="Ice thermal conductivity",dim1=vert_dim_poly,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"pp_uz",     poly%uz,     units="m a**-1",  long_name="Ice vertical velocity",   dim1="zeta_px_ac",dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"pp_advecxy",poly%advecxy,units="",         long_name="Ice horizontal advection",dim1=vert_dim_poly,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"pp_Q_strn", poly%Q_strn, units="J a-1 m-3",long_name="Ice strain heating",      dim1=vert_dim_poly,dim2="time",start=[1,n],ncid=ncid)
        
        ! Update variables (vectors) 
        call nc_write(filename,"enth",   vecs%enth,  units="J kg-1",    long_name="Ice enthalpy",               dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"T_ice",  vecs%T_ice,  units="K",        long_name="Ice temperature",            dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"omega",  vecs%omega, units="",          long_name="Ice water content (fraction)", dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"T_pmp",  vecs%T_pmp,  units="",         long_name="Ice pressure melting point", dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"T_prime",vecs%T_ice-vecs%T_pmp,units="K",long_name="Ice temperature",           dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        
        call nc_write(filename,"cp",     vecs%cp,     units="J kg-1 K-1",   long_name="Ice heat capacity",       dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"kt",     vecs%kt,     units="J a-1 m-1 K-1",long_name="Ice thermal conductivity",dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"uz",     vecs%uz,     units="m a**-1",  long_name="Ice vertical velocity",   dim1="zeta_ac",dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"advecxy",vecs%advecxy,units="",         long_name="Ice horizontal advection",dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        call nc_write(filename,"Q_strn", vecs%Q_strn, units="J a-1 m-3",long_name="Ice strain heating",      dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        
        call nc_write(filename,"t_dep", vecs%t_dep, units="a",long_name="Deposition time",      dim1=vert_dim,dim2="time",start=[1,n],ncid=ncid)
        
        call nc_write(filename,"H_cts",  ice%H_cts, units="m",long_name="CTS height",dim1="time",start=[n],ncid=ncid)
        
        ! Update variables (points) 
        call nc_write(filename,"smb",     ice%smb,units="m a**-1",long_name="Surface mass balance",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"bmb",     ice%bmb,units="m a**-1",long_name="Basal mass balance",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"Q_b",     ice%Q_b,units="",long_name="Basal heating",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"Q_geo",   ice%Q_geo,units="mW m**-2",long_name="Geothermal heat flux",dim1="time",start=[n],ncid=ncid)
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

    subroutine calc_zeta(zeta_aa,zeta_ac,zeta_scale,zeta_exp)
        ! Calculate the vertical layer-edge axis (vertical ac-nodes)
        ! and the vertical cell-center axis (vertical aa-nodes),
        ! including an extra zero-thickness aa-node at the base and surface

        implicit none 

        real(prec), intent(INOUT)  :: zeta_aa(:) 
        real(prec), intent(INOUT)  :: zeta_ac(:) 
        character(*), intent(IN)   :: zeta_scale 
        real(prec),   intent(IN)   :: zeta_exp 

        ! Local variables
        integer :: k, nz_aa, nz_ac 

        nz_aa  = size(zeta_aa)
        nz_ac  = size(zeta_ac)   ! == nz_aa - 1 

        ! Initially define a linear zeta scale 
        ! Base = 0.0, Surface = 1.0 
        do k = 1, nz_ac
            zeta_ac(k) = 0.0 + 1.0*(k-1)/real(nz_ac-1)
        end do 

        ! Scale zeta to produce different resolution through column if desired
        ! zeta_scale = ["linear","exp","wave"]
        select case(trim(zeta_scale))
            
            case("exp")
                ! Increase resolution at the base 
                zeta_ac = zeta_ac**(zeta_exp) 

            case("tanh")
                ! Increase resolution at base and surface 

                zeta_ac = tanh(1.0*pi*(zeta_ac-0.5))
                zeta_ac = zeta_ac - minval(zeta_ac)
                zeta_ac = zeta_ac / maxval(zeta_ac)

            case DEFAULT
            ! Do nothing, scale should be linear as defined above
        
        end select  
        
        ! Get zeta_aa (between zeta_ac values, as well as at the base and surface)
        zeta_aa(1) = 0.0 
        do k = 2, nz_aa-1
            zeta_aa(k) = 0.5 * (zeta_ac(k-1)+zeta_ac(k))
        end do 
        zeta_aa(nz_aa) = 1.0 

        return 

    end subroutine calc_zeta
    
    subroutine interp_enth_column(enth,zeta,enth0,zeta0,enth_pmp0,H_cts,H_ice)

        implicit none 

        real(prec), intent(OUT) :: enth(:)          ! New vertical axis (aa-nodes)
        real(prec), intent(IN)  :: zeta(:)          ! New vertical axis (aa-nodes)
        real(prec), intent(IN)  :: enth0(:)         ! Old vertical axis (aa-nodes)
        real(prec), intent(IN)  :: zeta0(:)         ! Old vertical axis (aa-nodes)
        real(prec), intent(IN)  :: enth_pmp0(:)     ! Old vertical axis (aa-nodes)
        real(prec), intent(IN)  :: H_cts 
        real(prec), intent(IN)  :: H_ice

        ! Local variables
        integer :: k, k0, k1, nz0, nz 
        real(prec) :: f_cts, enth_pmp_cts, f_lin 

        real(prec), parameter :: missing_value = -9999.0_prec 

        nz  = size(enth,1)
        nz0 = size(enth0,1) 

        f_cts = H_cts / max(H_ice,1e-5)
        if (H_ice .eq. 0.0) f_cts = 0.0 

        ! First set output enth to missing values 
        enth = missing_value 

        enth(1)  = enth0(1)
        enth(nz) = enth0(nz0)
        
        ! Loop over new vertical axis and interpolate points 
        do k = 2, nz-1 

            ! Find first index above current point, and index at or below current point
            do k1 = 2, nz0
                if (zeta0(k1) .gt. zeta(k)) exit 
            end do
            k0 = k1-1 

            if (enth0(k0) .ge. enth_pmp0(k0) .and. enth0(k1) .ge. enth_pmp0(k1)) then 
                ! Purely temperate, linear interpolation 

                enth(k) = interp_linear_internal(zeta0(k0:k1),enth0(k0:k1),zeta(k))

            else if (enth0(k0) .lt. enth_pmp0(k0) .and. enth0(k1) .lt. enth_pmp0(k1)) then
                ! Purely cold, linear interpolation 

                enth(k) = interp_linear_internal(zeta0(k0:k1),enth0(k0:k1),zeta(k))
                
            else 

                enth(k) = interp_linear_internal(zeta0(k0:k1),enth0(k0:k1),zeta(k))
                
                !enth_pmp_cts = interp_linear_internal(zeta0(k0:k1),enth_pmp0(k0:k1),f_cts)

!                 f_lin = (enth_pmp0(k0)-enth0(k0)) / ( (enth0(k1)-enth0(k0)) - (enth_pmp0(k1)-enth_pmp0(k0)) )
!                 if (f_lin .lt. 1e-2) f_lin = 0.0 
!                 f_cts = zeta0(k0) + f_lin*(zeta0(k)-zeta0(k0))
                
!                 write(*,*) k, zeta0(k0:k1),enth0(k0:k1),enth_pmp_cts,f_cts 
                
!                 if (zeta(k) .lt. f_cts) then 
!                     enth(k) = interp_linear_internal([zeta0(k0),f_cts],[enth0(k0),enth_pmp_cts],zeta(k))
!                 else 
!                     enth(k) = interp_linear_internal([f_cts,zeta0(k1)],[enth_pmp_cts,enth0(k1)],zeta(k))
!                 end if 

            end if 

        end do
        
        return 

    end subroutine interp_enth_column

    function interp_linear_internal(x,y,xout) result(yout)

        implicit none

        real(prec), intent(IN)  :: x(2), y(2), xout
        real(prec) :: yout
        real(prec) :: alph

        if ( xout .lt. x(1) .or. xout .gt. x(2) ) then
            write(*,*) "interp1: xout < x0 or xout > x1 !"
            write(*,*) "xout = ",xout
            write(*,*) "x0   = ",x(1)
            write(*,*) "x1   = ",x(2)
            stop
        end if

        alph = (xout - x(1)) / (x(2) - x(1))
        yout = y(1) + alph*(y(2) - y(1))

        return

    end function interp_linear_internal 

end program test_icetemp 
