program yelmo_trough
    ! For mimicking Feldmann and Levermann (2017, TC) 
    ! and for running mismip+ etc. 

    use mpi 

    use ncio 
    use yelmo 
    use deformation 

    implicit none 

    type(yelmo_class)     :: yelmo1

    character(len=56)  :: domain, grid_name  
    character(len=256) :: outfldr, file2D, file1D, file_restart
    character(len=512) :: path_par, path_const 
    character(len=56)  :: experiment, res  
    real(wp) :: time_init, time_end, time, dtt, dt1D_out, dt2D_out   
    integer    :: n

    ! Control parameters 
    real(wp) :: dx 
    real(wp) :: lx, ly, fc, dc, wc
    real(wp) :: x_cf 
    real(wp) :: Tsrf_const, smb_const, Qgeo_const  

    real(wp) :: s06_alpha, s06_H0, s06_W, s06_m
    real(wp) :: B, L  
    real(wp), allocatable :: ux_ref(:,:) 
    real(wp), allocatable :: tau_c_ref(:,:)

    real(wp) :: xmax, ymin, ymax 
    integer  :: i, j, nx, ny 
    
    real(8)  :: cpu_start_time, cpu_end_time, cpu_dtime  
    integer  :: perr 

    ! Initialize MPI
    call MPI_INIT(perr)


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
    
    ! Schoof domain parameters
    call nml_read(path_par,"ctrl_schoof","alpha",s06_alpha)      ! [m/m] Constant slope
    call nml_read(path_par,"ctrl_schoof","H0",   s06_H0)         ! [m]   Constant ice thickness
    call nml_read(path_par,"ctrl_schoof","W",    s06_W)          ! [m]   Half-width weak till
    call nml_read(path_par,"ctrl_schoof","m",    s06_m)          ! []    Exponent
    
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
    ymax =  ly/2.0_wp
    ymin = -ly/2.0_wp
    call yelmo_init_grid(yelmo1%grd,grid_name,units="km", &
                            x0=0.0_wp,dx=dx,nx=int(xmax/dx)+1, &
                            y0=ymin,dy=dx,ny=int((ymax-ymin)/dx)+1)

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

        case("RAYMOND")
            ! Raymond (2000) domain - constant slope slab

            ! ===== Intialize topography and set parameters =========
        
            yelmo1%bnd%z_bed = 10000.0_wp - s06_alpha*(yelmo1%grd%x - minval(yelmo1%grd%x))

            yelmo1%tpo%now%H_ice = s06_H0

            ! Define surface elevation 
            yelmo1%tpo%now%z_srf = yelmo1%bnd%z_bed + yelmo1%tpo%now%H_ice

            ! Define reference ice thickness (for prescribing boundary values, potentially)
            yelmo1%bnd%H_ice_ref = s06_H0 

            ! Define basal friction 
            yelmo1%dyn%now%cb_ref = 5.2e3 
            where(abs(yelmo1%grd%y) .gt. s06_W) yelmo1%dyn%now%cb_ref = 70e3 

            write(*,*) "RAYMOND: W      = ", s06_W 
            write(*,*) "RAYMOND: ATT    = ", yelmo1%mat%now%ATT(1,1,1)
            write(*,*) "RAYMOND: cb_ref = ", yelmo1%dyn%now%cb_ref(1,1)

            ! Initialiaze ux values to be safe too 
            yelmo1%dyn%now%ux_b   = 500.0 
            yelmo1%dyn%now%ux_bar = 500.0 
            yelmo1%dyn%now%ux_s   = 500.0 
            
        case("SLAB-S06")
            ! Schoof (2006) domain - constant slope slab

            ! ===== Intialize topography and set parameters =========
        
            yelmo1%bnd%z_bed = 10000.0_wp - s06_alpha*(yelmo1%grd%x - minval(yelmo1%grd%x))

            yelmo1%tpo%now%H_ice = s06_H0
            yelmo1%bnd%H_ice_ref = s06_H0 

            ! Define surface elevation 
            yelmo1%tpo%now%z_srf = yelmo1%bnd%z_bed + yelmo1%tpo%now%H_ice

            ! Calculate analytical stream function to get tau_c and ux

            allocate(ux_ref(yelmo1%grd%nx,yelmo1%grd%ny))
            allocate(tau_c_ref(yelmo1%grd%nx,yelmo1%grd%ny))
            
            call SSA_Schoof2006_analytical_solution_yelmo(ux_ref, tau_c_ref, yelmo1%grd%y, &
                                    s06_alpha,s06_H0,yelmo1%mat%par%rf_const,s06_W,s06_m, &
                                    yelmo1%mat%par%n_glen,rho_ice, g)
            
            ! Assign analytical values (tau_c as a boundary condition, ux as initial condition)
            yelmo1%dyn%now%cb_ref = tau_c_ref
            
            ! Assign initial velocity values to help achieve quicker convergence...
            yelmo1%dyn%now%ux_b   = ux_ref 
            yelmo1%dyn%now%ux_bar = ux_ref 
            yelmo1%dyn%now%ux_s   = ux_ref 
            yelmo1%dyn%now%uy_b   = 0.0_wp 
            yelmo1%dyn%now%uy_bar = 0.0_wp 
            yelmo1%dyn%now%uy_s   = 0.0_wp 

            ! Determine constant L too, for diagnostic output
            L = s06_W / ((1.0_wp+s06_m)**(1.0_wp/s06_m))

            write(*,*) "SLAB-S06: H0          = ", s06_H0 
            write(*,*) "SLAB-S06: alpha       = ", s06_alpha 
            write(*,*) "SLAB-S06: W           = ", s06_W 
            write(*,*) "SLAB-S06: L           = ", L 
            write(*,*) "SLAB-S06: m           = ", s06_m 
            write(*,*) "SLAB-S06: rho g       = ", rho_ice, g
            write(*,*) "SLAB-S06: f           = ", (rho_ice*g*s06_H0)*s06_alpha
            write(*,*) "SLAB-S06: ATT         = ", yelmo1%mat%par%rf_const
            write(*,*) "SLAB-S06: cb_ref      = ", yelmo1%dyn%now%cb_ref(1,1)
            write(*,*) "SLAB-S06: tau_c_ref   = ", tau_c_ref(1,1)
            write(*,*) "SLAB-S06: max(ux_ref) = ", maxval(ux_ref)

        case("TROUGH-F17")
            ! Feldmann and Levermann (2017) domain 

            call trough_f17_topo_init(yelmo1%bnd%z_bed,yelmo1%tpo%now%H_ice,yelmo1%tpo%now%z_srf, &
                                    yelmo1%grd%xc*1e-3,yelmo1%grd%yc*1e-3,fc,dc,wc,x_cf)
        
        case("MISMIP+") 
            ! MISMIP+ domain 

            call trough_mismipp_topo_init(yelmo1%bnd%z_bed,yelmo1%tpo%now%H_ice,yelmo1%tpo%now%z_srf, &
                                    yelmo1%grd%xc*1e-3,yelmo1%grd%yc*1e-3,fc,dc,wc,x_cf)
        
        case("SLAB-SHELF")
            ! Constant slab slope with an ice shelf

            ! call trough_f17_topo_init(yelmo1%bnd%z_bed,yelmo1%tpo%now%H_ice,yelmo1%tpo%now%z_srf, &
            !                         yelmo1%grd%xc*1e-3,yelmo1%grd%yc*1e-3,fc,dc,wc,x_cf)
            
            call slab_topo_init(yelmo1%bnd%z_bed,yelmo1%tpo%now%H_ice,yelmo1%tpo%now%z_srf, &
                                    yelmo1%grd%xc*1e-3,yelmo1%grd%yc*1e-3)


        case DEFAULT 

            write(*,*) "yelmo_trough:: Error: domain not recognized: "//trim(domain)
            stop 

    end select 


    ! Define calving front 
    call define_calving_front(yelmo1%bnd%calv_mask,yelmo1%grd%x*1e-3,x_cf)

    ! Initialize the yelmo state (dyn,therm,mat)
    call yelmo_init_state(yelmo1,time=time_init,thrm_method="robin-cold")

    ! Write initial state 
    call write_step_2D(yelmo1,file2D,time=time_init) 

    ! 1D file 
    call yelmo_write_reg_init(yelmo1,file1D,time_init=time_init,units="years",mask=yelmo1%bnd%ice_allowed)
    call yelmo_write_reg_step(yelmo1,file1D,time=time)  
    
    ! Advance timesteps
    do n = 1, ceiling((time_end-time_init)/dtt)

        ! Get current time 
        time = time_init + n*dtt

if (.FALSE.) then
        if (trim(domain) .eq. "SLAB-SHELF" .and. time .ge. 3e3) then 

            ! ! Define calving front 
            ! x_cf = 540.0_wp 
            ! call define_calving_front(yelmo1%bnd%calv_mask,yelmo1%grd%x*1e-3,x_cf)

            ! Kill all floating ice now
            yelmo1%tpo%par%calv_flt_method = "kill"

        end if 
end if 

        ! == Yelmo ice sheet ===================================================
        call yelmo_update(yelmo1,time)

        ! == MODEL OUTPUT =======================================================
        if (mod(nint(time*100),nint(dt2D_out*100))==0) then  
            call write_step_2D(yelmo1,file2D,time=time)    
        end if 

        if (mod(nint(time*100),nint(dt1D_out*100))==0) then 
            call yelmo_write_reg_step(yelmo1,file1D,time=time) 
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
    
    ! Finalize MPI
    call MPI_FINALIZE(perr)

contains
    
    subroutine define_calving_front(calv_mask,xx,x_cf)
        ! Define a calving mask in the x direction where 
        ! beyond the position x_cf ice will be calved. 

        implicit none 

        logical, intent(OUT) :: calv_mask(:,:) 
        real(wp), intent(IN) :: xx(:,:) 
        real(wp), intent(IN) :: x_cf 

        calv_mask = .FALSE. 
        where (xx .ge. x_cf) calv_mask = .TRUE. 
    
        return 

    end subroutine define_calving_front

    subroutine slab_topo_init(z_bed,H_ice,z_srf,xc,yc)

        implicit none 

        real(prec), intent(OUT) :: z_bed(:,:) 
        real(prec), intent(OUT) :: H_ice(:,:) 
        real(prec), intent(OUT) :: z_srf(:,:) 
        real(prec), intent(IN)  :: xc(:)
        real(prec), intent(IN)  :: yc(:)

        ! Local variables 
        integer :: i, j, nx, ny 

        nx = size(z_bed,1)
        ny = size(z_bed,2)
        
        ! Define bedrock as a slope
        do j = 1, ny
            z_bed(:,j) = -100.0 - xc
        end do 
        
        ! Set ice thickness to 500 m everywhere initially 
        H_ice = 500.0

        ! Remove ice from deep bed
        where(z_bed .lt. -500.0) H_ice = 0.0 

        ! Adjust for floating ice later, for now assume fully grounded
        z_srf = z_bed + H_ice 

        return 

    end subroutine slab_topo_init

    subroutine trough_f17_topo_init(z_bed,H_ice,z_srf,xc,yc,fc,dc,wc,x_cf)

        implicit none 

        real(wp), intent(OUT) :: z_bed(:,:) 
        real(wp), intent(OUT) :: H_ice(:,:) 
        real(wp), intent(OUT) :: z_srf(:,:) 
        real(wp), intent(IN)  :: xc(:) 
        real(wp), intent(IN)  :: yc(:)  
        real(wp), intent(IN)  :: fc 
        real(wp), intent(IN)  :: dc 
        real(wp), intent(IN)  :: wc 
        real(wp), intent(IN)  :: x_cf 

        ! Local variables 
        integer :: i, j, nx, ny 
        real(wp) :: zb_x, zb_y 
        real(wp) :: e1, e2 

        real(wp), parameter :: zb_deep = -720.0_wp 

        nx = size(z_bed,1) 
        ny = size(z_bed,2) 

        write(*,*) "params: ", ly,fc,dc,wc,x_cf

        ! == Bedrock elevation == 
        do j = 1, ny
        do i = 1, nx 
            
            ! x-direction 
            zb_x = -150.0_wp - 0.84*abs(xc(i))

            ! y-direction 
            e1 = -2.0*(yc(j)-wc)/fc 
            e2 =  2.0*(yc(j)+wc)/fc 
            zb_y = ( dc / (1.0+exp(e1)) ) + ( dc / (1.0+exp(e2)) ) 

            ! Convolution 
            z_bed(i,j) = max(zb_x + zb_y, zb_deep)

        end do
        end do  

        ! == Ice thickness == 
        H_ice = 50.0_wp 
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

        real(wp), intent(OUT) :: z_bed(:,:) 
        real(wp), intent(OUT) :: H_ice(:,:) 
        real(wp), intent(OUT) :: z_srf(:,:) 
        real(wp), intent(IN)  :: xc(:) 
        real(wp), intent(IN)  :: yc(:)  
        real(wp), intent(IN)  :: fc 
        real(wp), intent(IN)  :: dc 
        real(wp), intent(IN)  :: wc 
        real(wp), intent(IN)  :: x_cf 

        ! Local variables 
        integer  :: i, j, nx, ny 
        real(wp) :: zb_x, zb_y 
        real(wp) :: e1, e2 

        real(wp) :: x1 

        real(wp), parameter :: zb_deep = -720.0_wp 
        real(wp), parameter :: xbar    =  300.0_wp          ! [km] Characteristic along-flow length scale of the bedrock
        real(wp), parameter :: b0      = -150.00_wp 
        real(wp), parameter :: b2      = -728.80_wp 
        real(wp), parameter :: b4      =  343.91_wp 
        real(wp), parameter :: b6      =  -50.57_wp 


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
        H_ice = 50.0_wp 
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
        real(wp), intent(IN) :: time

        ! Local variables
        integer  :: ncid, n, i, j, nx, ny  
        real(wp) :: time_prev 

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
        
        call nc_write(filename,"mask_frnt",ylmo%tpo%now%mask_frnt,units="",long_name="Ice-front mask", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taul_int_acx",ylmo%dyn%now%taul_int_acx,units="Pa m",long_name="Vertically integrated lateral stress (x)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taul_int_acy",ylmo%dyn%now%taul_int_acy,units="Pa m",long_name="Vertically integrated lateral stress (y)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"H_ice_pred",ylmo%tpo%now%H_ice_pred,units="m",long_name="Ice thickness (predicted)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"H_ice_corr",ylmo%tpo%now%H_ice_corr,units="m",long_name="Ice thickness (corrected)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"dzsdt",ylmo%tpo%now%dzsdt,units="m/a",long_name="Surface elevation change", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dHidt",ylmo%tpo%now%dHidt,units="m/a",long_name="Ice thickness change", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"H_grnd",ylmo%tpo%now%H_grnd,units="m",long_name="Ice thickness overburden", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"bmb",ylmo%tpo%now%bmb,units="m/a",long_name="Basal mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

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

        call nc_write(filename,"ssa_mask_acx",ylmo%dyn%now%ssa_mask_acx,units="1",long_name="SSA mask (acx)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"ssa_mask_acy",ylmo%dyn%now%ssa_mask_acy,units="1",long_name="SSA mask (acy)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"cb_ref",ylmo%dyn%now%cb_ref,units="--",long_name="Bed friction scalar", &
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
        call nc_write(filename,"visc_eff",ylmo%dyn%now%visc_eff,units="Pa a m",long_name="Effective viscosity (SSA)", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        call nc_write(filename,"dzsdx",ylmo%tpo%now%dzsdx,units="m/m",long_name="Surface gradient, x-direction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dzsdy",ylmo%tpo%now%dzsdy,units="m/m",long_name="Surface gradient, y-direction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"taud_acx",ylmo%dyn%now%taud_acx,units="Pa",long_name="Driving stress, x-direction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taud_acy",ylmo%dyn%now%taud_acy,units="Pa",long_name="Driving stress, y-direction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"taub_acx",ylmo%dyn%now%taub_acx,units="Pa",long_name="Basal stress, x-direction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taub_acy",ylmo%dyn%now%taub_acy,units="Pa",long_name="Basal stress, y-direction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
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

        call nc_write(filename,"ux_b",ylmo%dyn%now%ux_b,units="m/a",long_name="Basal sliding velocity (x)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy_b",ylmo%dyn%now%uy_b,units="m/a",long_name="Basal sliding velocity (y)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_b",ylmo%dyn%now%uxy_b,units="m/a",long_name="Basal sliding velocity magnitude", &
                     dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"ux_bar",ylmo%dyn%now%ux_bar,units="m/a",long_name="Vertically averaged velocity (x)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy_bar",ylmo%dyn%now%uy_bar,units="m/a",long_name="Vertically averaged velocity (y)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_bar",ylmo%dyn%now%uxy_bar,units="m/a",long_name="Vertically averaged velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"ux_s",ylmo%dyn%now%ux_s,units="m/a",long_name="Surface velocity (x)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy_s",ylmo%dyn%now%uy_s,units="m/a",long_name="Surface velocity (y)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_s",ylmo%dyn%now%uxy_s,units="m/a",long_name="Surface velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"qq_gl_acx",ylmo%dyn%now%qq_gl_acx,units="m2/a",long_name="Grounding line flux (x)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"qq_gl_acy",ylmo%dyn%now%qq_gl_acy,units="m2/a",long_name="Grounding line flux (y)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"ux",ylmo%dyn%now%ux,units="m/a",long_name="Horizontal velocity (x)", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"uy",ylmo%dyn%now%uy,units="m/a",long_name="Horizontal velocity (y)", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"uxy",ylmo%dyn%now%uxy,units="m/a",long_name="Horizontal velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"uz",ylmo%dyn%now%uz,units="m/a",long_name="Vertical velocity", &
                      dim1="xc",dim2="yc",dim3="zeta_ac",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"uz_star",ylmo%thrm%now%uz_star,units="m/a",long_name="Advective vertical velocity", &
                      dim1="xc",dim2="yc",dim3="zeta_ac",dim4="time",start=[1,1,1,n],ncid=ncid)

        call nc_write(filename,"advecxy",ylmo%thrm%now%advecxy,units="m/a",long_name="Horizontal advection", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

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
        
        call nc_write(filename,"bmb_shlf",ylmo%bnd%bmb_shlf,units="m/a",long_name="Basal mass balance (shelf)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! No time dimension::

        call nc_write(filename,"z_bed",ylmo%bnd%z_bed,units="m",long_name="Bedrock elevation", &
                      dim1="xc",dim2="yc",start=[1,1],ncid=ncid)
        
        call nc_write(filename,"Q_geo",ylmo%bnd%Q_geo,units="mW/m^2",long_name="Geothermal heat flux", &
                      dim1="xc",dim2="yc",start=[1,1],ncid=ncid)
        
        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine write_step_2D

    ! == Analytical solution by Schoof 2006 for the "SSA_icestream" benchmark experiment
  elemental subroutine SSA_Schoof2006_analytical_solution_yelmo(u, tauc, y, tantheta, h0, A_flow, W, m,  &
                            n_glen,ice_density, grav)
      
    implicit none
    
    ! In/output variables:
    real(wp),                            intent(OUT)   :: u             ! Ice velocity in the x-direction
    real(wp),                            intent(OUT)   :: tauc          ! Till yield stress
    real(wp),                            intent(IN)    :: y             ! y-coordinate
    real(wp),                            intent(IN)    :: tantheta      ! Surface slope in the x-direction
    real(wp),                            intent(IN)    :: h0            ! Ice thickness
    real(wp),                            intent(IN)    :: A_flow        ! Ice flow factor
    real(wp),                            intent(IN)    :: W             ! Ice-stream half width (m)
    real(wp),                            intent(IN)    :: m             ! Ice stream exponent
    real(wp),                            intent(IN)    :: n_glen
    real(wp),                            intent(IN)    :: ice_density
    real(wp),                            intent(IN)    :: grav
    
    ! Local variables:
    real(dp) :: B, f, L, ua, ub, uc, ud, ue   
    real(dp) :: ux, taud, H, yy 

    ! Calculate the gravitational driving stress f
    f = ice_density * grav * h0 * tantheta
    
    ! Calculate the ice hardness factor B
    B = A_flow**(-1._dp/n_glen)
    
    ! Determine constant L (ice-stream width)
    L = W / ((1.0_dp+m)**(1.0_dp/m))

    ! Calculate the till yield stress across the stream
    tauc = f * ABS(y/L)**m
    
    taud = f 
    H    = h0
    yy   = y 
    ux = -2.0 * taud**3 * L**4 / (B**3 * H**3) * ( ((yy/L)**4 - (m+1.0)**(4.0/m))/4.0 - 3.0*( abs(yy/L)**(m+4.0) &
    - (m+1.0)**(1.0+4.0/m) )/((m+1.0)*(m+4.0)) + 3.0*( abs(yy/L)**(2.0*m+4.0) - (m+1.0)**(2.0+4.0/m) )/((m+1.0)**2*(2.0*m+4.0)) &
    - ( abs(yy/L)**(3.0*m+4.0) - (m+1.0)**(3.0+4.0/m) )/ ( (m+1.0)**3*(3.0*m+4.0)) )

if (.TRUE.) then

    u = ux 

else 

    ! Calculate the analytical solution for u
    ua = -2._dp * f**3 * L**4 / (B**3 * h0**3)
    ub = ( 1._dp / 4._dp                           ) * (   (y/L)**     4._dp  - (m+1._dp)**(       4._dp/m) )
    uc = (-3._dp / ((m+1._dp)    * (      m+4._dp))) * (ABS(y/L)**(  m+4._dp) - (m+1._dp)**(1._dp+(4._dp/m)))
    ud = ( 3._dp / ((m+1._dp)**2 * (2._dp*m+4._dp))) * (ABS(y/L)**(2*m+4._dp) - (m+1._dp)**(2._dp+(4._dp/m)))
    ue = (-1._dp / ((m+1._dp)**3 * (3._dp*m+4._dp))) * (ABS(y/L)**(3*m+4._dp) - (m+1._dp)**(3._dp+(4._dp/m)))
    u = ua * (ub + uc + ud + ue)
    
end if

    ! Outside the ice-stream, velocity is zero
    IF (ABS(y) > W) u = 0._dp
    
    end subroutine SSA_Schoof2006_analytical_solution_yelmo

  ! == Analytical solution by Schoof 2006 for the "SSA_icestream" benchmark experiment
  ELEMENTAL SUBROUTINE SSA_Schoof2006_analytical_solution(U, tauc, y, tantheta, h0, A_flow, L, m,  &
                            n_glen,ice_density, grav)
      
    IMPLICIT NONE
    
    ! In/output variables:
    REAL(wp),                            INTENT(OUT)   :: U             ! Ice velocity in the x-direction
    REAL(wp),                            INTENT(OUT)   :: tauc          ! Till yield stress
    REAL(wp),                            INTENT(IN)    :: y             ! y-coordinate
    REAL(wp),                            INTENT(IN)    :: tantheta      ! Surface slope in the x-direction
    REAL(wp),                            INTENT(IN)    :: h0            ! Ice thickness
    REAL(wp),                            INTENT(IN)    :: A_flow        ! Ice flow factor
    REAL(wp),                            INTENT(IN)    :: L             ! Ice-stream width (m), default 40e3
    REAL(wp),                            INTENT(IN)    :: m             ! Ice stream exponent
    REAL(wp),                            INTENT(IN)    :: n_glen
    REAL(wp),                            INTENT(IN)    :: ice_density
    REAL(wp),                            INTENT(IN)    :: grav
    
    ! Local variables:
    REAL(wp)                                           :: B, f, W, ua, ub, uc, ud, ue   
    
    ! Calculate the gravitational driving stress f
    f = ice_density * grav * h0 * tantheta
    
    ! Calculate the ice hardness factor B
    B = A_flow**(-1._dp/n_glen)
    
    ! Calculate the "ice stream half-width" W
    W = L * (m+1._dp)**(1._dp/m)
    
    ! Calculate the till yield stress across the stream
    tauc = f * ABS(y/L)**m
    
    ! Calculate the analytical solution for u
    ua = -2._dp * f**3 * L**4 / (B**3 * h0**3)
    ub = ( 1._dp / 4._dp                           ) * (   (y/L)**     4._dp  - (m+1._dp)**(       4._dp/m) )
    uc = (-3._dp / ((m+1._dp)    * (      m+4._dp))) * (ABS(y/L)**(  m+4._dp) - (m+1._dp)**(1._dp+(4._dp/m)))
    ud = ( 3._dp / ((m+1._dp)**2 * (2._dp*m+4._dp))) * (ABS(y/L)**(2*m+4._dp) - (m+1._dp)**(2._dp+(4._dp/m)))
    ue = (-1._dp / ((m+1._dp)**3 * (3._dp*m+4._dp))) * (ABS(y/L)**(3*m+4._dp) - (m+1._dp)**(3._dp+(4._dp/m)))
    u = ua * (ub + uc + ud + ue)
    
    ! Outside the ice-stream, velocity is zero
    IF (ABS(y) > w) U = 0._dp
    
  END SUBROUTINE SSA_Schoof2006_analytical_solution

end program yelmo_trough



