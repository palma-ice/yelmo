

module basal_hydro_simple
    ! This module manages the variables associated with the calculation
    ! of basal hydrology (H_w, p_w) in the model 

    use nml 
    !use yelmo_defs, only : sp, dp, prec

    implicit none 

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the precision of the library (sp,dp)
    integer,  parameter :: prec = sp 

    real(prec), parameter :: sec_year  = 365.0*24.0*60.0*60.0   ! [s/a]
    real(prec), parameter :: G         = 9.81                   ! [m/s^2]
    real(prec), parameter :: pi        = 3.14159265359

    real(prec), parameter :: rho_ice   = 910.0                  ! [kg/m^3]
    real(prec), parameter :: rho_sw    = 1028.0                 ! [kg/m^3] 
    real(prec), parameter :: rho_w     = 1000.0                 ! [kg/m^3] 

    real(prec), parameter :: rho_ice_G = rho_ice*G              ! [kg/m^3 (m/s^2)]
    real(prec), parameter :: rho_sw_G  = rho_sw*G               ! [kg/m^3 (m/s^2)]
    real(prec), parameter :: rho_w_G   = rho_w*G                ! [kg/m^3 (m/s^2)]

    
    type hydro_param_class
        integer    :: init_method       ! Method to initialize hydro fields 
        integer    :: method            ! Method to treat basal water
        real(prec) :: H_w_init          ! [m] Initial water layer thickness, for init_method=1
        real(prec) :: H_w_max           ! [m] Maximum allowed water thickness       
        real(prec) :: till_rate         ! [m/a] Till infiltration (drainage rate)
    end type 

    type hydro_state_class
        real(prec) :: time, dt  
        real(prec), allocatable :: H_w(:,:)            ! [m]  Water layer thickness
        real(prec), allocatable :: p_w(:,:)            ! [Pa] Water pressure
        real(prec), allocatable :: dHwdt(:,:)          ! [m/a]  Water layer thickness change
        
    end type 

    type hydro_class
        type(hydro_param_class) :: par 
        type(hydro_state_class) :: now 
    end type

    private
    public :: hydro_class
    public :: hydro_init 
    public :: hydro_init_state 
    public :: hydro_update
!     public :: hydro_end

    public :: calc_basal_water_local

contains 

    subroutine hydro_update(hyd,H_ice,f_grnd,bmb_w,time)
        ! Update the state of the hydro variables 

        implicit none 

        type(hydro_class), intent(INOUT) :: hyd 
        real(prec), intent(IN) :: H_ice(:,:)      ! [m] Ice thickness 
        real(prec), intent(IN) :: f_grnd(:,:)     ! [-] Grounded fraction of grid cell
        real(prec), intent(IN) :: bmb_w(:,:)      ! [m/a water equiv.] Basal mass balance of water (bmb_w = -bmb_ice, since ice sheet is defined negative for melting ice)
        real(prec), intent(IN) :: time            ! [a] Current external time to advance to

        ! Local variables 
        integer :: nx, ny 
        real(prec), allocatable :: H_w_prev(:,:) 

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        allocate(H_w_prev(nx,ny))

        ! Determine current time step and time 
        hyd%now%dt   = max(time - hyd%now%time,0.0)
        hyd%now%time = time 


        if (hyd%now%dt .gt. 0.0) then 
                            
            ! Store initial H_w field 
            H_w_prev = hyd%now%H_w 

            ! Update the basal water layer thickness depending on the method desired
            select case(hyd%par%method)

                case(0)
                    ! Constant H_water at initial value 
                    
                    ! Pass, nothing happens

                case(1)
                    ! Local mass balance of H_w 

                    call calc_basal_water_local(hyd%now%H_w,hyd%now%dHwdt,H_ice,bmb_w,f_grnd, &
                                                    hyd%now%dt,hyd%par%till_rate,hyd%par%H_w_max)

                case DEFAULT 

                    write(*,*) "hydro_update:: error: method must be in one of [0,1]."
                    write(*,*) "method = ", hyd%par%method 
                    stop 

            end select 

            ! Update the water pressure 
            hyd%now%p_w = 0.0 

            ! Get rate of change of H_w field 
            hyd%now%dHwdt = (hyd%now%H_w - H_w_prev) / hyd%now%dt 

        end if 

        return 

    end subroutine hydro_update 

    subroutine hydro_init(hyd,filename,nx,ny)

        implicit none 

        type(hydro_class), intent(INOUT) :: hyd 
        character(len=*),  intent(IN)    :: filename 
        integer,           intent(IN)    :: nx, ny 

        ! Load parameter options 
        call hydro_par_load(hyd%par,filename)

        ! Initialize state object 
        call hydro_allocate(hyd%now,nx,ny)

        ! Set fields to zero for now
        ! (use hydro_init_state later to intialize properly)
        hyd%now%H_w     = 0.0 
        hyd%now%p_w     = 0.0 

        hyd%now%time    = 1e10     ! Set time way into the future, for now 
        hyd%now%dt      = 0.0 

        return 

    end subroutine hydro_init 

    subroutine hydro_init_state(hyd,H_ice,f_grnd,time)
        ! Initialize the state of the hydro fields
        ! (only after calling hydro_init)

        implicit none 

        type(hydro_class), intent(INOUT) :: hyd 
        real(prec), intent(IN) :: H_ice(:,:) 
        real(prec), intent(IN) :: f_grnd(:,:)  
        real(prec), intent(IN) :: time 
        
        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: im1, ip1, jm1, jp1 

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        ! Initialize using the method of choice (par%init_method)

        select case(hyd%par%init_method)

            case(-1) 
                ! Do nothing - initialized externally 

            case(0)
                ! Intially zero everywhere

                hyd%now%H_w = 0.0 

            case(1)
                ! Spatially constant intial value

                hyd%now%H_w = hyd%par%H_w_init 

            case DEFAULT 

                write(*,*) "hydro_init_state:: error: initialization method must be one of [0,1]."
                write(*,*) "init_method = ", hyd%par%init_method 
                stop

        end select


        ! Set floating points to the maximum water thickness
        where (f_grnd .eq. 0.0) hyd%now%H_w = hyd%par%H_w_max  


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
                
                hyd%now%H_w(i,j) = hyd%par%H_w_max 

            end if 

        end do 
        end do 
        
        ! Set the current time and time step 
        hyd%now%time = time 
        hyd%now%dt   = 0.0 

        return 

    end subroutine hydro_init_state 


    subroutine hydro_par_load(par,filename,init)

        implicit none 

        type(hydro_param_class)  :: par
        character(len=*)         :: filename 
        logical, optional        :: init

        ! Local variables 
        logical :: init_pars 

        init_pars = .FALSE.
        if (present(init)) init_pars = .TRUE. 
        
        call nml_read(filename,"basal_hydro_simple","init_method",   par%init_method,  init=init_pars)
        call nml_read(filename,"basal_hydro_simple","method",        par%method,       init=init_pars)
        call nml_read(filename,"basal_hydro_simple","H_w_init",      par%H_w_init,     init=init_pars)
        call nml_read(filename,"basal_hydro_simple","H_w_max",       par%H_w_max,      init=init_pars)
        call nml_read(filename,"basal_hydro_simple","till_rate",     par%till_rate,    init=init_pars)

        return 

    end subroutine hydro_par_load

    subroutine hydro_allocate(now,nx,ny)

        implicit none 

        type(hydro_state_class), intent(INOUT) :: now 
        integer :: nx, ny 

        ! First make sure fields are deallocated
        call hydro_deallocate(now)

        ! Allocate fields to desired dimensions
        allocate(now%H_w(nx,ny))
        allocate(now%p_w(nx,ny))
        allocate(now%dHwdt(nx,ny))

        ! Initialize all fields to zero
        now%H_w     = 0.0
        now%p_w     = 0.0
        now%dHwdt   = 0.0 

        return 

    end subroutine hydro_allocate 

    subroutine hydro_deallocate(now)

        implicit none 

        type(hydro_state_class), intent(INOUT) :: now 

        if (allocated(now%H_w))     deallocate(now%H_w)
        if (allocated(now%p_w))     deallocate(now%p_w)
        if (allocated(now%dHwdt))   deallocate(now%dHwdt)
        
        return 

    end subroutine hydro_deallocate 


    ! ==== BASAL HYDROLOGY PHYSICS ===============================

    subroutine calc_basal_water_local(H_w,dHwdt,H_ice,bmb_w,f_grnd,dt,till_rate,H_w_max)
        ! Calculate the basal water layer thickness based on a simple local 
        ! water balance: dHw/dt = bmb_w - 
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
        dHwdt = H_w 

        where (f_grnd .gt. 0.0 .and. H_ice .gt. 0.0)
            ! Grounded ice point

            ! Update mass balance of H_w
            H_w = H_w + dt*bmb_w - dt*till_rate

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
            dHwdt = (dHwdt - H_w) / dt 
        else 
            dHwdt = 0.0_prec 
        end if 

        return 

    end subroutine calc_basal_water_local

end module basal_hydro_simple 



