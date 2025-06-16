program yelmo_calving

    use nml 
    use ncio  
    use yelmo 
    use lsf_module

    use calving_benchmarks
    
    implicit none

    type(yelmo_class) :: yelmo1
    type(yelmo_class) :: yelmo_ref

    type control_type    
        character(len=256) :: outfldr
        character(len=256) :: file2D, file1D
        character(len=256) :: file_restart
        character(len=512) :: path_par 
        character(len=56)  :: exp

        real(wp) :: time_init, time_end, time, dtt
        real(wp) :: dt2D_out, dt1D_out

        real(wp) :: dx

        ! Internal parameters
        character(len=56)  :: domain
        character(len=56)  :: grid_name
        real(wp) :: x0, x1
        integer  :: nx
        integer  :: ny
    end type

    type(control_type) :: ctl
    real(wp) :: time 
    integer  :: n 
    
    real(8) :: cpu_start_time, cpu_end_time, cpu_dtime  
    
    ! Start timing 
    call yelmo_cpu_time(cpu_start_time)

    ! Assume program is running from the output folder
    ctl%outfldr = "./"

    ! Determine the parameter file from the command line 
    call yelmo_load_command_line_args(ctl%path_par)
    !path_par   = trim(outfldr)//"yelmo_calving.nml" 

    ! Define input and output locations 
    ctl%file1D       = trim(ctl%outfldr)//"yelmo1D.nc"
    ctl%file2D       = trim(ctl%outfldr)//"yelmo2D.nc"
    ctl%file_restart = trim(ctl%outfldr)//"yelmo_restart.nc"

    
    ! Define the domain, grid and experiment from parameter file
    call nml_read(ctl%path_par,"ctl","exp",         ctl%exp)            ! "exp1", "exp2", "exp3", "exp4", "exp5"
    call nml_read(ctl%path_par,"ctl","dx",          ctl%dx)             ! [km] Grid resolution 

    ! Timing parameters 
    call nml_read(ctl%path_par,"ctl","time_init",   ctl%time_init)      ! [yr] Starting time
    call nml_read(ctl%path_par,"ctl","time_end",    ctl%time_end)       ! [yr] Ending time
    call nml_read(ctl%path_par,"ctl","dtt",         ctl%dtt)            ! [yr] Main loop time step 
    call nml_read(ctl%path_par,"ctl","dt2D_out",    ctl%dt2D_out)       ! [yr] Frequency of 2D output 
    ctl%dt1D_out = ctl%dtt  ! Set 1D output to frequency of main loop timestep 

    ! Now set internal parameters ===

    ! Define domain and grid size based on experiment
    select case(trim(ctl%exp))
        case("exp1","exp2")
            ctl%domain = "circular"
            ctl%x0 = -800.0
            ctl%x1 =  800.0
        case("exp3","exp4","exp5")
            ctl%domain = "thule"
            ctl%x0 = -800.0
            ctl%x1 =  800.0
        case("advection")
            ctl%domain = "advection"
            ctl%x0     = -800.0
            ctl%x1     = 800.0
        case DEFAULT
            write(*,*) "ctl.exp = ",trim(ctl%domain), " not recognized."
            stop
    end select

    ! Get grid size
    ctl%nx = (ctl%x1-ctl%x0) / ctl%dx + 1
    ctl%ny = ctl%nx

    ! Get grid name
    write(ctl%grid_name,"(a,i2,a2)") trim(ctl%domain)//"-",int(ctl%dx),"KM"
   
    ! === Initialize ice sheet model =====

    ! First, define grid 
    call yelmo_init_grid(yelmo1%grd,ctl%grid_name,units="km",dx=ctl%dx,nx=ctl%nx,dy=ctl%dx,ny=ctl%nx)

    ! Initialize data objects (without loading topography, which will be defined inline below)
    call yelmo_init(yelmo1,filename=ctl%path_par,grid_def="none",time=ctl%time_init, &
                    load_topo=.FALSE.,domain=ctl%domain,grid_name=ctl%grid_name)

    ! === Define initial topography ===
    call calvmip_init(yelmo1%bnd%z_bed,yelmo1%grd%x,yelmo1%grd%y,yelmo1%par%domain)

    ! advection test
    if (.not. yelmo1%par%use_restart) then
        ! If no restart, set ice thickness to zero
        yelmo1%tpo%now%H_ice = 0.0
        yelmo1%tpo%now%z_srf = yelmo1%bnd%z_bed 
        select case(trim(ctl%exp))
            case("advection")
            call CircularDomain(yelmo1%tpo%now%lsf,yelmo1%bnd%z_bed,yelmo1%tpo%par%dx)
        case DEFAULT 
            call LSFinit(yelmo1%tpo%now%lsf,yelmo1%tpo%now%H_ice,yelmo1%tpo%now%z_srf,yelmo1%tpo%par%dx)
        end select
    end if

    ! === Define additional boundary conditions =====

    yelmo1%bnd%z_sl     = 0.0
    yelmo1%bnd%bmb_shlf = 0.0  
    yelmo1%bnd%T_shlf   = yelmo1%bnd%c%T0  
    yelmo1%bnd%H_sed    = 0.0 

    yelmo1%bnd%T_srf    = 223.15 
    yelmo1%bnd%Q_geo    = 42.0

    select case(trim(ctl%exp))
        case("advection")
            yelmo1%bnd%smb      = 0.0
        case DEFAULT 
            yelmo1%bnd%smb      = 0.3
    end select    

    ! Check boundary values 
    call yelmo_print_bound(yelmo1%bnd)

    ! Initialize state variables (dyn,therm,mat)
    call yelmo_init_state(yelmo1,time=ctl%time_init,thrm_method="robin")

    ! == Write initial state ==
 
    ! 2D file 
    call yelmo_write_init(yelmo1,ctl%file2D,time_init=ctl%time_init,units="years")
    call yelmo_write_step(yelmo1,ctl%file2D,time=ctl%time_init)  
    
    ! 1D file 
    call yelmo_write_reg_init(yelmo1,ctl%file1D,time_init=ctl%time_init,units="years",mask=yelmo1%bnd%ice_allowed)
    call yelmo_write_reg_step(yelmo1,ctl%file1D,time=ctl%time_init) 

    ! Store default parameters
    yelmo_ref = yelmo1

    ! Set calving mask if needed
    if (trim(yelmo1%tpo%par%calv_flt_method) .eq. "kill-pos") then
        call set_calving_mask(yelmo1%bnd%calv_mask,yelmo1%grd%x,yelmo1%grd%y,r_lim=750e3_wp)
    end if

    ! Advance timesteps
    do n = 1, ceiling((ctl%time_end-ctl%time_init)/ctl%dtt)

        ! Get current time 
        time = ctl%time_init + n*ctl%dtt
        
        ! if (time .lt. 10e3) then
        !     yelmo1%dyn%par%solver = "sia"
        ! else
        !     yelmo1%dyn%par%solver = yelmo_ref%dyn%par%solver
        ! end if

        ! == Yelmo ice sheet ===================================================
        call yelmo_update(yelmo1,time)
        
        ! == MODEL OUTPUT =======================================================
        if (mod(nint(time*100),nint(ctl%dt2D_out*100))==0) then 
            if (.FALSE.) then
                call yelmo_write_step(yelmo1,ctl%file2D,time=time)
            else
                ! CalvingMIP variables
                call write_2D_calvingmip(yelmo1,ctl%file2D,time=time)   
            end if
        end if 

        if (mod(nint(time*100),nint(ctl%dt1D_out*100))==0) then 
            if (.TRUE.) then
                call yelmo_write_reg_step(yelmo1,ctl%file1D,time=time) 
            else
                ! CalvingMIP variables
                call write_1D_calvingmip(yelmo1,ctl%file1D,time=time)
            end if
        end if 

    end do

    ! Write a restart file too
    call yelmo_restart_write(yelmo1,ctl%file_restart,time=time)

    ! Finalize program
    call yelmo_end(yelmo1,time=time)

    ! Stop timing 
    call yelmo_cpu_time(cpu_end_time,cpu_start_time,cpu_dtime)
    
    write(*,"(a,f12.3,a)") "Time  = ",cpu_dtime/60.0 ," min"
    write(*,"(a,f12.1,a)") "Speed = ",(1e-3*(ctl%time_end-ctl%time_init))/(cpu_dtime/3600.0), " kiloyears / hr"
    
contains

    subroutine CircularDomain(LSF,zbed,dx)
        
        implicit none
    
        real(wp), intent(OUT) :: LSF(:,:)      ! LSF mask
        real(wp), intent(IN)  :: zbed(:,:)    
        real(wp), intent(IN)  :: dx            ! Model resolution [m]
        
        ! Internal variables
        real(wp) :: rc
        integer  :: i,j,nx,ny
    
        nx = size(zbed,1)
        ny = size(zbed,2)
        rc = 10.0_wp ! grid points below zero
    
        do j=1,ny
        do i=1,nx
    
        LSF(i,j) = (sqrt((0.5*(nx+1)-i)**2 + (0.5*(ny+1)-j)**2) - rc)*dx*1e-3 
    
        end do
        end do
    
        return
    
    end subroutine CircularDomain

    subroutine write_1D_calvingmip(dom,filename,time)

        implicit none
    
        type(yelmo_class), intent(IN) :: dom
        character(len=*),  intent(IN) :: filename
        real(wp),          intent(IN) :: time
    
        ! Local variables
        type(yregions_class) :: reg
            
        integer  :: ncid, n
        real(wp) :: rho_ice
        real(wp) :: density_corr
        real(wp) :: m3yr_to_kgs
        real(wp) :: ismip6_correction
        real(wp) :: yr_to_sec
            
        real(wp) :: m3_km3
        real(wp) :: m2_km2 
    
        integer  :: npts_tot
        integer  :: npts_grnd
        integer  :: npts_flt
        integer  :: npts_grl
        integer  :: npts_frnt 
    
        real(wp) :: dx
        real(wp) :: dy 
        real(wp) :: smb_tot 
        real(wp) :: bmb_tot 
        real(wp) :: bmb_shlf_t 
    
        real(wp) :: A_ice_grl 
        real(wp) :: flux_grl 
        real(wp) :: A_ice_frnt 
        real(wp) :: calv_flt 
        real(wp) :: flux_frnt 
        
        logical, allocatable :: mask_tot(:,:) 
        logical, allocatable :: mask_grnd(:,:)
        logical, allocatable :: mask_flt(:,:) 
        logical, allocatable :: mask_grl(:,:) 
        logical, allocatable :: mask_frnt(:,:) 
            
        dx = dom%grd%dx 
        dy = dom%grd%dy 
    
        ! Allocate variables
    
        allocate(mask_tot(dom%grd%nx,dom%grd%ny))
        allocate(mask_grnd(dom%grd%nx,dom%grd%ny))
        allocate(mask_flt(dom%grd%nx,dom%grd%ny))
        allocate(mask_grl(dom%grd%nx,dom%grd%ny))
        allocate(mask_frnt(dom%grd%nx,dom%grd%ny))
    
        ! === Data conversion factors ========================================
    
        rho_ice             = 917.0             ! ice density kg/m3
        m3yr_to_kgs         = 3.2e-5            ! m3/yr of pure water to kg/s
        density_corr        = rho_ice/1000.0    ! ice density correction with pure water
        ismip6_correction   = m3yr_to_kgs*density_corr
        yr_to_sec           = 31556952.0
            
        m3_km3              = 1e-9 
        m2_km2              = 1e-6 
            
        ! 1. Determine regional values of variables 
    
        ! Take the global regional data object that 
        ! is calculated over the whole domain at each timestep 
        reg = dom%reg
    
        ! Assign masks of interest
    
        mask_tot  = (dom%tpo%now%H_ice .gt. 0.0) 
        mask_grnd = (dom%tpo%now%H_ice .gt. 0.0 .and. dom%tpo%now%f_grnd .gt. 0.0)
        mask_flt  = (dom%tpo%now%H_ice .gt. 0.0 .and. dom%tpo%now%f_grnd .eq. 0.0)
        mask_grl  = (dom%tpo%now%f_grnd_bmb .gt. 0.0 .and. dom%tpo%now%f_grnd_bmb .lt. 1.0)         
        mask_frnt = (dom%tpo%now%H_ice .gt. 0.0 .and. dom%tpo%now%mask_bed .eq. 4)
    
        ! Determine number of points at grl and frnt
        npts_tot  = count(mask_tot)
        npts_grnd = count(mask_grnd)
        npts_flt  = count(mask_flt)
        npts_grl  = count(mask_grl)      
        npts_frnt = count(mask_frnt) 
    
        ! Calculate additional variables of interest for ISMIP6
    
        ! To do: these variables should be defined locally!
        ! ISMIP6 boundary: jablasco
        smb_tot     = sum(dom%bnd%smb,    mask=mask_tot)*(dx*dy)    ! m^3/yr: flux
        bmb_tot     = sum(dom%tpo%now%bmb,mask=mask_tot)*(dx*dy)    ! m^3/yr: flux
        bmb_shlf_t  = sum(dom%tpo%now%bmb,mask=mask_flt)*(dx*dy)    ! m^3/yr: flux
    
        if (npts_grl .gt. 0) then
    
            A_ice_grl = count(dom%tpo%now%H_ice .gt. 0.0 .and. mask_grl)*dx*dy*m2_km2 ! [km^2]
            !flux_grl  = sum(dom%tpo%now%H_ice,mask=mask_grl)*(dx*dy)                 ! m^3/yr: flux
            flux_grl  = sum(dom%dyn%now%uxy_bar*dom%tpo%now%H_ice,mask=mask_grl)*dx   ! m^3/yr: flux
    
            ! ajr, why only *dx above? 
    
        else
    
            A_ice_grl = 0.0_wp
            flux_grl  = 0.0_wp
    
        end if
     
        ! ===== Frontal ice-shelves variables =====
    
        if (npts_frnt .gt. 0) then
    
            A_ice_frnt = count(dom%tpo%now%H_ice .gt. 0.0 .and. mask_frnt)*dx*dy*m2_km2         ! [km^2]
            calv_flt   = sum(dom%tpo%now%cmb_flt*dom%tpo%now%H_ice,mask=mask_frnt)*dx          ! m^3/yr: flux [m-1 yr-1]
            flux_frnt  = calv_flt+sum(dom%tpo%now%fmb*dom%tpo%now%H_ice,mask=mask_frnt)*dx      ! m^3/yr: flux [m-1 yr-1]
    
            ! ajr, why only *dx above? 
    
        else
    
            A_ice_frnt = 0.0_wp
            calv_flt   = 0.0_wp 
            flux_frnt  = 0.0_wp
    
        end if
    
    
        ! 2. Begin writing step 
    
        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)
    
        ! Determine current writing time step 
        n = nc_time_index(filename,"time",time,ncid)
    
        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)
    
        call nc_write(filename,"V_sl",reg%V_sl*1e-6,units="1e6 km^3",long_name="Ice volume above flotation", &
                        dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"V_sle",reg%V_sle,units="m sle",long_name="Sea-level equivalent volume", &
                        dim1="time",start=[n],ncid=ncid)
    
        ! ===== Mass ===== 
        call nc_write(filename,"lim",reg%V_ice*rho_ice*1e9,units="kg",long_name="Total ice mass", &
                    standard_name="land_ice_mass",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"limnsw",reg%V_sl*rho_ice*1e9,units="kg",long_name="Mass above flotation", &
                    standard_name="land_ice_mass_not_displacing_sea_water",dim1="time",start=[n],ncid=ncid)
    
        ! ===== Area ====
        call nc_write(filename,"iareagr",reg%A_ice_g*1e6,units="m^2",long_name="Grounded ice area", &
                    standard_name="grounded_ice_sheet_area",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"iareafl",reg%A_ice_f*1e6,units="m^2",long_name="Floating ice area", &
                    standard_name="floating_ice_shelf_area",dim1="time",start=[n],ncid=ncid)
    
        ! ==== Fluxes ====
        call nc_write(filename,"tendacabf",smb_tot*ismip6_correction,units="kg s-1",long_name="Total SMB flux", &
                    standard_name="tendency_of_land_ice_mass_due_to_surface_mass_balance",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"tendlibmassbf ",bmb_tot*ismip6_correction,units="kg s-1",long_name="Total BMB flux", &
                    standard_name="tendency_of_land_ice_mass_due_to_basal_mass_balance",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"tendlibmassbffl",bmb_shlf_t*ismip6_correction,units="kg s-1",long_name="Total BMB flux beneath floating ice", &
                    standard_name="tendency_of_land_ice_mass_due_to_basal_mass_balance",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"tendlicalvf",calv_flt*ismip6_correction,units="kg s-1",long_name="Total calving flux", &
                    standard_name="tendency_of_land_ice_mass_due_to_calving",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"tendlifmassbf",flux_frnt*ismip6_correction,units="kg s-1",long_name="Total calving and ice front melting flux", &
                    standard_name="tendency_of_land_ice_mass_due_to_calving_and_ice_front_melting",dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"tendligroundf",flux_grl*ismip6_correction,units="kg s-1",long_name="Total grounding line flux", &
                    standard_name="tendency_of_grounded_ice_mass",dim1="time",start=[n],ncid=ncid)
    
        ! Close the netcdf file
        call nc_close(ncid)
    
        return
    
    end subroutine write_1D_calvingmip

    subroutine write_2D_calvingmip(ylmo,filename,time)

        implicit none
    
        type(yelmo_class), intent(IN) :: ylmo
        character(len=*),  intent(IN) :: filename
        real(wp),          intent(IN) :: time
    
        ! Local variables
        integer  :: ncid, n, i, j
            
        ! CalvingMIP variables
        integer,  allocatable :: mask_clvmip(:,:)
        real(wp), allocatable :: ux_bar_aa(:,:)
        real(wp), allocatable :: uy_bar_aa(:,:)

        ! Allocate and initialize local arrays
        allocate(mask_clvmip(ylmo%grd%nx,ylmo%grd%ny))
        allocate(ux_bar_aa(ylmo%grd%nx,ylmo%grd%ny))
        allocate(uy_bar_aa(ylmo%grd%nx,ylmo%grd%ny)) 
    
        mask_clvmip = 0
        ux_bar_aa   = 0.0_wp
        uy_bar_aa   = 0.0_wp
    
        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)
    
        ! Determine current writing time step 
        n = nc_time_index(filename,"Time",time,ncid)
    
        ! Update the time step
        call nc_write(filename,"Time",time,dim1="Time",start=[n],count=[1],ncid=ncid)

        ! Compute mask
        where(ylmo%tpo%now%H_ice .eq. 0.0_wp) mask_clvmip = 3
        where(ylmo%tpo%now%H_ice .gt. 0.0_wp .and. ylmo%tpo%now%f_grnd .eq. 0.0_wp) mask_clvmip = 2
        where(ylmo%tpo%now%H_ice .gt. 0.0_wp .and. ylmo%tpo%now%f_grnd .gt. 0.0_wp) mask_clvmip = 1

        ! convert velocities into aa-nodes
        do i=2, ylmo%grd%nx-1
        do j=2, ylmo%grd%ny-1
            ux_bar_aa(i,j) = 0.5*(ylmo%dyn%now%ux_bar(i,j)+ylmo%dyn%now%ux_bar(i-1,j))
            uy_bar_aa(i,j) = 0.5*(ylmo%dyn%now%uy_bar(i,j)+ylmo%dyn%now%uy_bar(i,j-1))
        end do
        end do

        where(ylmo%tpo%now%H_ice .eq. 0.0_wp) ux_bar_aa = 0.0_wp
        where(ylmo%tpo%now%H_ice .eq. 0.0_wp) uy_bar_aa = 0.0_wp

        ! Write CalvingMIP variables variables


        call nc_write(filename,"xvelmean",ux_bar_aa,start=[1,1,n],units="m a-1",long_name="X velocity", &
                        standard_name="land_ice_vertical_mean_x_velocity", dim1="X",dim2="Y",dim3="Time",ncid=ncid)
        call nc_write(filename,"yvelmean",uy_bar_aa,start=[1,1,n],units="m a-1",long_name="Y velocity", &
                standard_name="land_ice_vertical_mean_y_velocity", dim1="X",dim2="Y",dim3="Time",ncid=ncid)
        call nc_write(filename,"lithk",ylmo%tpo%now%H_ice,start=[1,1,n],units="m",long_name="Ice thickness", &
                    standard_name="land_ice_thickness", dim1="X",dim2="Y",dim3="Time",ncid=ncid)
        call nc_write(filename,"mask",mask_clvmip,start=[1,1,n],units="",long_name="Ice mask", &
                        dim1="X",dim2="Y",dim3="=Time",ncid=ncid)
        call nc_write(filename,"topg",ylmo%bnd%z_bed,start=[1,1,n],units="m",long_name="Bedrock height", &
                    standard_name="bedrock_altimetry", dim1="X",dim2="Y",dim3="Time",ncid=ncid)
    
        ! Close the netcdf file
        call nc_close(ncid)
    
        return
    
    end subroutine write_2D_calvingmip

end program yelmo_calving
