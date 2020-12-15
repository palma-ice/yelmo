program yelmo_trough
    ! For mimicking Feldmann and Levermann (2017, TC) 
    ! and for running mismip+ etc. 

    use ncio 
    use yelmo 
    use deformation 

    implicit none 

    type(yelmo_class)     :: yelmo1

    character(len=56)  :: domain, grid_name  
    character(len=256) :: outfldr, file2D, file1D, file_restart
    character(len=512) :: path_par, path_const 
    character(len=56)  :: experiment, res  
    real(prec) :: time_init, time_end, time, dtt, dt1D_out, dt2D_out   
    integer    :: n  

    ! Control parameters 
    real(prec) :: dx 
    real(prec) :: lx, ly, fc, dc, wc
    real(prec) :: x_cf 
    real(prec) :: Tsrf_const, smb_const, Qgeo_const  

    real(prec) :: xmax, ymin, ymax 
    integer    :: i, j, nx, ny 
    
    real(8) :: cpu_start_time, cpu_end_time, cpu_dtime  
    
    ! Start timing 
    call yelmo_cpu_time(cpu_start_time)
    
    ! Assume program is running from the output folder
    outfldr = "./"

    ! Determine the parameter file from the command line 
    call yelmo_load_command_line_args(path_par)
    !path_par   = trim(outfldr)//"yelmo_TROUGH-F17.nml" 
    
    ! Define input and output locations 
    path_const = trim(outfldr)//"yelmo_const_TROUGH.nml"
    file2D     = trim(outfldr)//"yelmo2D.nc"
    file1D     = trim(outfldr)//"yelmo1D.nc"
    
    ! Define the domain, grid and experiment from parameter file
    call nml_read(path_par,"ctrl","domain",       domain)        ! TROUGH-F17, MISMIP+

    ! Timing parameters 
    call nml_read(path_par,"ctrl","time_init",    time_init)     ! [yr] Starting time
    call nml_read(path_par,"ctrl","time_end",     time_end)      ! [yr] Ending time
    call nml_read(path_par,"ctrl","dtt",          dtt)           ! [yr] Main loop time step 
    call nml_read(path_par,"ctrl","dt1D_out",     dt1D_out)      ! [yr] Frequency of 1D output 
    call nml_read(path_par,"ctrl","dt2D_out",     dt2D_out)      ! [yr] Frequency of 2D output 

    ! Domain parameters
    call nml_read(path_par,"ctrl","dx",           dx)            ! [km] Grid resolution ! must be multiple of xmax and ymax!!
    call nml_read(path_par,"ctrl","lx",           lx)            ! [km] Trough parameter
    call nml_read(path_par,"ctrl","ly",           ly)            ! [km] Trough parameter
    call nml_read(path_par,"ctrl","fc",           fc)            ! [km] Trough parameter
    call nml_read(path_par,"ctrl","dc",           dc)            ! [km] Trough parameter
    call nml_read(path_par,"ctrl","wc",           wc)            ! [km] Trough parameter
    call nml_read(path_par,"ctrl","x_cf",         x_cf)          ! [km] Trough parameter
    
    ! Simulation parameters 
    call nml_read(path_par,"ctrl","Tsrf_const",   Tsrf_const)    ! [degC]  Surface temperature
    call nml_read(path_par,"ctrl","smb_const",    smb_const)     ! [m/yr]  Surface mass balance
    call nml_read(path_par,"ctrl","Qgeo_const",   Qgeo_const)    ! [mW/m2] Geothermal heat flux
    
    
    ! Define default grid name for completeness 
    grid_name = trim(domain)

    ! General initialization of yelmo constants (used globally)
    call yelmo_global_init(path_const)

    ! Define the domain and grid
    xmax =  lx 
    ymax =  ly/2.0_prec 
    ymin = -ly/2.0_prec 
    call yelmo_init_grid(yelmo1%grd,grid_name,units="km",x0=0.0,dx=dx,nx=int(xmax/dx)+1,y0=ymin,dy=dx,ny=int((ymax-ymin)/dx)+1)

    ! === Initialize ice sheet model =====

    ! Initialize data objects
    call yelmo_init(yelmo1,filename=path_par,grid_def="none",time=time_init,load_topo=.FALSE., &
                        domain=domain,grid_name=grid_name)

    ! Load boundary values

    yelmo1%bnd%z_sl     = 0.0
    yelmo1%bnd%bmb_shlf = 0.0 
    yelmo1%bnd%T_shlf   = T0  
    yelmo1%bnd%H_sed    = 0.0 

    yelmo1%bnd%T_srf    = T0 + Tsrf_const   ! [K] 
    yelmo1%bnd%smb      = smb_const         ! [m/yr]
    yelmo1%bnd%Q_geo    = Qgeo_const        ! [mW/m2] 

    ! Check boundary values 
    call yelmo_print_bound(yelmo1%bnd)

    ! Initialize output file 
    call yelmo_write_init(yelmo1,file2D,time_init=time_init,units="years")
    
    ! Intialize topography 
    select case(trim(domain)) 

        case("TROUGH-F17")
            ! Feldmann and Levermann (2017) domain 

            call trough_f17_topo_init(yelmo1%bnd%z_bed,yelmo1%tpo%now%H_ice,yelmo1%tpo%now%z_srf, &
                                    yelmo1%grd%xc*1e-3,yelmo1%grd%yc*1e-3,fc,dc,wc,x_cf)
        
        case("MISMIP+") 


            call trough_mismipp_topo_init(yelmo1%bnd%z_bed,yelmo1%tpo%now%H_ice,yelmo1%tpo%now%z_srf, &
                                    yelmo1%grd%xc*1e-3,yelmo1%grd%yc*1e-3,fc,dc,wc,x_cf)
        
        case DEFAULT 

            write(*,*) "yelmo_trough:: Error: domain not recognized: "//trim(domain)
            stop 

    end select 


    ! Define calving front 
    yelmo1%bnd%calv_mask = .FALSE. 
    where (yelmo1%grd%x*1e-3 .ge. x_cf) yelmo1%bnd%calv_mask = .TRUE. 

    ! Initialize the yelmo state (dyn,therm,mat)
    call yelmo_init_state(yelmo1,path_par,time=time_init,thrm_method="robin-cold")

    ! Write initial state 
    call write_step_2D(yelmo1,file2D,time=time_init) 

    ! 1D file 
    call write_yreg_init(yelmo1,file1D,time_init=time_init,units="years",mask=yelmo1%bnd%ice_allowed)
    call write_yreg_step(yelmo1,file1D,time=time)  
    
    ! Advance timesteps
    do n = 1, ceiling((time_end-time_init)/dtt)

        ! Get current time 
        time = time_init + n*dtt

        ! == Yelmo ice sheet ===================================================
        call yelmo_update(yelmo1,time)


        ! == MODEL OUTPUT =======================================================
        if (mod(nint(time*100),nint(dt2D_out*100))==0) then  
            call write_step_2D(yelmo1,file2D,time=time)    
        end if 

        if (mod(nint(time*100),nint(dt1D_out*100))==0) then 
            call write_yreg_step(yelmo1,file1D,time=time) 
        end if

        if (mod(time,10.0)==0 .and. (.not. yelmo_log)) then
            write(*,"(a,f14.4)") "yelmo:: time = ", time
        end if  

    end do 

    ! Write summary 
    write(*,*) "====== "//trim(domain)//" ======="
    write(*,*) "nz, H0 = ", yelmo1%par%nz_aa, maxval(yelmo1%tpo%now%H_ice)

    ! Finalize program
    call yelmo_end(yelmo1,time=time)

    ! Stop timing 
    call yelmo_cpu_time(cpu_end_time,cpu_start_time,cpu_dtime)
    
    write(*,"(a,f12.3,a)") "Time  = ",cpu_dtime/60.0 ," min"
    write(*,"(a,f12.1,a)") "Speed = ",(1e-3*(time_end-time_init))/(cpu_dtime/3600.0), " kiloyears / hr"
    
contains

    subroutine trough_f17_topo_init(z_bed,H_ice,z_srf,xc,yc,fc,dc,wc,x_cf)

        implicit none 

        real(prec), intent(OUT) :: z_bed(:,:) 
        real(prec), intent(OUT) :: H_ice(:,:) 
        real(prec), intent(OUT) :: z_srf(:,:) 
        real(prec), intent(IN)  :: xc(:) 
        real(prec), intent(IN)  :: yc(:)  
        real(prec), intent(IN)  :: fc 
        real(prec), intent(IN)  :: dc 
        real(prec), intent(IN)  :: wc 
        real(prec), intent(IN)  :: x_cf 

        ! Local variables 
        integer :: i, j, nx, ny 
        real(prec) :: zb_x, zb_y 
        real(prec) :: e1, e2 

        real(prec), parameter :: zb_deep = -720.0_prec 

        nx = size(z_bed,1) 
        ny = size(z_bed,2) 

        write(*,*) "params: ", ly,fc,dc,wc,x_cf

        ! == Bedrock elevation == 
        do j = 1, ny
        do i = 1, nx 
            
            ! x-direction 
            zb_x = -150.0_prec - 0.84*abs(xc(i))

            ! y-direction 
            e1 = -2.0*(yc(j)-wc)/fc 
            e2 =  2.0*(yc(j)+wc)/fc 
            zb_y = ( dc / (1.0+exp(e1)) ) + ( dc / (1.0+exp(e2)) ) 

            ! Convolution 
            z_bed(i,j) = max(zb_x + zb_y, zb_deep)

        end do
        end do  

        ! == Ice thickness == 
        H_ice = 50.0_prec 
        do j = 1, ny 
            where(xc .gt. x_cf) H_ice(:,j) = 0.0 
        end do 

        ! == Surface elevation == 
        z_srf = z_bed + H_ice

        where(z_srf .lt. 0.0) z_srf = 0.0 

        return 

    end subroutine trough_f17_topo_init

    subroutine trough_mismipp_topo_init(z_bed,H_ice,z_srf,xc,yc,fc,dc,wc,x_cf)

        implicit none 

        real(prec), intent(OUT) :: z_bed(:,:) 
        real(prec), intent(OUT) :: H_ice(:,:) 
        real(prec), intent(OUT) :: z_srf(:,:) 
        real(prec), intent(IN)  :: xc(:) 
        real(prec), intent(IN)  :: yc(:)  
        real(prec), intent(IN)  :: fc 
        real(prec), intent(IN)  :: dc 
        real(prec), intent(IN)  :: wc 
        real(prec), intent(IN)  :: x_cf 

        ! Local variables 
        integer :: i, j, nx, ny 
        real(prec) :: zb_x, zb_y 
        real(prec) :: e1, e2 

        real(prec) :: x1 

        real(prec), parameter :: zb_deep = -720.0_prec 
        real(prec), parameter :: xbar    =   300.0_prec         ! [km] Characteristic along-flow length scale of the bedrock
        real(prec), parameter :: b0      = -150.00_prec 
        real(prec), parameter :: b2      = -728.80_prec 
        real(prec), parameter :: b4      =  343.91_prec 
        real(prec), parameter :: b6      =  -50.57_prec 


        nx = size(z_bed,1) 
        ny = size(z_bed,2) 

        write(*,*) "params: ", ly,fc,dc,wc,x_cf

        ! == Bedrock elevation == 
        do j = 1, ny
        do i = 1, nx 
            
            ! x-direction 
            x1 = xc(i) / xbar 
            zb_x = b0 + b2*x1**2 + b4*x1**4 + b6*x1**6 

            ! y-direction 
            e1 = -2.0*(yc(j)-wc)/fc 
            e2 =  2.0*(yc(j)+wc)/fc 
            zb_y = ( dc / (1.0+exp(e1)) ) + ( dc / (1.0+exp(e2)) ) 

            ! Convolution 
            z_bed(i,j) = max(zb_x + zb_y, zb_deep)

        end do
        end do  

        ! == Ice thickness == 
        H_ice = 50.0_prec 
        do j = 1, ny 
            where(xc .gt. x_cf) H_ice(:,j) = 0.0 
        end do 

        ! == Surface elevation == 
        z_srf = z_bed + H_ice

        where(z_srf .lt. 0.0) z_srf = 0.0 

        return 

    end subroutine trough_mismipp_topo_init

    subroutine write_step_2D(ylmo,filename,time)

        implicit none 
        
        type(yelmo_class), intent(IN) :: ylmo
        character(len=*),  intent(IN) :: filename
        real(prec), intent(IN) :: time

        ! Local variables
        integer    :: ncid, n, i, j, nx, ny  
        real(prec) :: time_prev 

        nx = ylmo%tpo%par%nx 
        ny = ylmo%tpo%par%ny 

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_size(filename,"time",ncid)
        call nc_read(filename,"time",time_prev,start=[n],count=[1],ncid=ncid) 
        if (abs(time-time_prev).gt.1e-5) n = n+1 

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        ! Write model metrics (model speed, dt, eta)
        call yelmo_write_step_model_metrics(filename,ylmo,n,ncid)

        ! == yelmo_topography ==
        call nc_write(filename,"H_ice",ylmo%tpo%now%H_ice,units="m",long_name="Ice thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf",ylmo%tpo%now%z_srf,units="m",long_name="Surface elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_bed",ylmo%tpo%now%mask_bed,units="",long_name="Bed mask", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"N_eff",ylmo%dyn%now%N_eff,units="Pa",long_name="Effective pressure", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"dzsrfdt",ylmo%tpo%now%dzsrfdt,units="m/a",long_name="Surface elevation change", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dHicedt",ylmo%tpo%now%dHicedt,units="m/a",long_name="Ice thickness change", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"H_grnd",ylmo%tpo%now%H_grnd,units="m",long_name="Ice thickness overburden", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"bmb",ylmo%tpo%now%bmb,units="m/a",long_name="Basal mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"dist_grline",ylmo%tpo%now%dist_grline,units="km", &
!                       long_name="Distance to nearest grounding-line point", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"is_grline",ylmo%tpo%now%is_grline,units="--",long_name="Grounding-line point", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"is_grz",ylmo%tpo%now%is_grz,units="--",long_name="Grounding-line zone", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"f_grnd",ylmo%tpo%now%f_grnd,units="1",long_name="Grounded fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_grnd_acx",ylmo%tpo%now%f_grnd_acx,units="1",long_name="Grounded fraction (acx)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_grnd_acy",ylmo%tpo%now%f_grnd_acy,units="1",long_name="Grounded fraction (acy)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_ice",ylmo%tpo%now%f_ice,units="1",long_name="Ice-covered fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"calv",ylmo%tpo%now%calv,units="m/a",long_name="Calving rate", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! == yelmo_thermodynamics ==
        call nc_write(filename,"T_ice",ylmo%thrm%now%T_ice,units="K",long_name="Ice temperature", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"T_prime",ylmo%thrm%now%T_ice-ylmo%thrm%now%T_pmp,units="deg C",long_name="Homologous ice temperature", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"f_pmp",ylmo%thrm%now%f_pmp,units="1",long_name="Fraction of grid point at pmp", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"T_prime_b",ylmo%thrm%now%T_prime_b,units="deg C",long_name="Homologous basal ice temperature", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"H_w",ylmo%thrm%now%H_w,units="m",long_name="Basal water layer", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! == yelmo_material ==
!         call nc_write(filename,"visc_int",ylmo%mat%now%visc_int,units="Pa a m",long_name="Vertically integrated viscosity", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"visc",ylmo%mat%now%visc,units="Pa a",long_name="Viscosity", &
!                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"ATT_bar",ylmo%mat%now%ATT_bar,units="a^-1 Pa^-3",long_name="Vertically averaged rate factor", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"ATT",ylmo%mat%now%ATT,units="a^-1 Pa^-3",long_name="Rate factor", &
!                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        
        call nc_write(filename,"Q_ice_b",ylmo%thrm%now%Q_ice_b,units="J a-1 m-2",long_name="Basal ice heat flux", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"Q_strn",ylmo%thrm%now%Q_strn/(rho_ice*ylmo%thrm%now%cp),units="K a-1",long_name="Strain heating", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        call nc_write(filename,"Q_b",ylmo%thrm%now%Q_b,units="J a-1 m-2",long_name="Basal frictional heating", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! == yelmo_dynamics ==

!         call nc_write(filename,"ssa_mask_acx",ylmo%dyn%now%ssa_mask_acx,units="1",long_name="SSA mask (acx)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"ssa_mask_acy",ylmo%dyn%now%ssa_mask_acy,units="1",long_name="SSA mask (acy)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"cf_ref",ylmo%dyn%now%cf_ref,units="--",long_name="Bed friction scalar", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"c_bed",ylmo%dyn%now%c_bed,units="Pa",long_name="Bed friction coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"beta",ylmo%dyn%now%beta,units="Pa a m-1",long_name="Basal friction coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"beta_acx",ylmo%dyn%now%beta_acx,units="Pa a m-1",long_name="Basal friction coefficient (acx)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"beta_acy",ylmo%dyn%now%beta_acy,units="Pa a m-1",long_name="Basal friction coefficient (acy)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"visc_eff_int",ylmo%dyn%now%visc_eff_int,units="Pa a m",long_name="Depth-integrated effective viscosity (SSA)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"taud_acx",ylmo%dyn%now%taud_acx,units="Pa",long_name="Driving stress, x-direction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taud_acy",ylmo%dyn%now%taud_acy,units="Pa",long_name="Driving stress, y-direction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"taub_acx",ylmo%dyn%now%taub_acx,units="Pa",long_name="Basal stress, x-direction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taub_acy",ylmo%dyn%now%taub_acy,units="Pa",long_name="Basal stress, y-direction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
!         call nc_write(filename,"sigma_horiz_sq",ylmo%dyn%now%sigma_horiz_sq,units="1",long_name="Horizontal stress components squared", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"lhs_x",ylmo%dyn%now%lhs_x,units="Pa",long_name="Shear reduction (x)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"lhs_y",ylmo%dyn%now%lhs_y,units="Pa",long_name="Shear reduction (y)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"lhs_xy",ylmo%dyn%now%lhs_xy,units="Pa",long_name="Shear reduction magnitude", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"duxdz",ylmo%dyn%now%duxdz,units="1/a",long_name="Vertical shear (x)", &
!                        dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
!         call nc_write(filename,"duydz",ylmo%dyn%now%duydz,units="1/a",long_name="Vertical shear (y)", &
!                        dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        
!         call nc_write(filename,"ux_i_bar",ylmo%dyn%now%ux_i_bar,units="m/a",long_name="Internal shear velocity (x)", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"uy_i_bar",ylmo%dyn%now%uy_i_bar,units="m/a",long_name="Internal shear velocity (y)", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"uxy_i_bar",ylmo%dyn%now%uxy_i_bar,units="m/a",long_name="Internal shear velocity magnitude", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"ux_b",ylmo%dyn%now%ux_b,units="m/a",long_name="Basal sliding velocity (x)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"uy_b",ylmo%dyn%now%uy_b,units="m/a",long_name="Basal sliding velocity (y)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"uxy_b",ylmo%dyn%now%uxy_b,units="m/a",long_name="Basal sliding velocity magnitude", &
!                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"ux_bar",ylmo%dyn%now%ux_bar,units="m/a",long_name="Vertically averaged velocity (x)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy_bar",ylmo%dyn%now%uy_bar,units="m/a",long_name="Vertically averaged velocity (y)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_bar",ylmo%dyn%now%uxy_bar,units="m/a",long_name="Vertically averaged velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"qq_gl_acx",ylmo%dyn%now%qq_gl_acx,units="m2/a",long_name="Grounding line flux (x)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"qq_gl_acy",ylmo%dyn%now%qq_gl_acy,units="m2/a",long_name="Grounding line flux (y)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"ux",ylmo%dyn%now%ux,units="m/a",long_name="Horizontal velocity (x)", &
!                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
!         call nc_write(filename,"uy",ylmo%dyn%now%uy,units="m/a",long_name="Horizontal velocity (y)", &
!                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
!         call nc_write(filename,"uxy",ylmo%dyn%now%uxy,units="m/a",long_name="Horizontal velocity magnitude", &
!                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
!         call nc_write(filename,"uz",ylmo%dyn%now%uz,units="m/a",long_name="Vertical velocity", &
!                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

!         call nc_write(filename,"f_vbvs",ylmo%dyn%now%f_vbvs,units="1",long_name="Basal to surface velocity fraction", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"f_shear_bar",ylmo%mat%now%f_shear_bar,units="1",long_name="Vertically averaged shearing fraction", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"de",ylmo%mat%now%strn%de,units="a^-1",long_name="Vertically averaged strain rate", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        ! == yelmo_bound ==

!         call nc_write(filename,"z_sl",ylmo%bnd%z_sl,units="m",long_name="Sea level rel. to present", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
!         call nc_write(filename,"H_sed",ylmo%bnd%H_sed,units="m",long_name="Sediment thickness", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"smb",ylmo%bnd%smb,units="m/a",long_name="Surface mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"T_srf",ylmo%bnd%T_srf,units="K",long_name="Surface temperature", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
!         call nc_write(filename,"bmb_shlf",ylmo%bnd%bmb_shlf,units="m/a",long_name="Basal mass balance (shelf)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! No time dimension::

        call nc_write(filename,"z_bed",ylmo%bnd%z_bed,units="m",long_name="Bedrock elevation", &
                      dim1="xc",dim2="yc",start=[1,1],ncid=ncid)
        
        call nc_write(filename,"Q_geo",ylmo%bnd%Q_geo,units="mW/m^2",long_name="Geothermal heat flux", &
                      dim1="xc",dim2="yc",start=[1,1],ncid=ncid)
        
        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine write_step_2D
    
end program yelmo_trough 



