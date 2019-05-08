

program yelmo_test

    use ncio 
    use yelmo 
    
    implicit none 

    type(yelmo_class)      :: yelmo1

    character(len=256) :: outfldr, file1D, file2D, file_restart, domain 
    character(len=512) :: path_par, path_const  
    real(prec) :: time_init, time_end, time_equil, time, dtt, dt1D_out, dt2D_out   
    integer    :: n
    real(4) :: cpu_start_time, cpu_end_time 

    ! No-ice mask (to impose additional melting)
    logical, allocatable :: mask_noice(:,:)  

    ! Concavity field 
    real(prec), allocatable :: channels(:,:) 

    ! Start timing 
    call cpu_time(cpu_start_time)

    ! Determine the parameter file from the command line 
    call yelmo_load_command_line_args(path_par)

    ! Timing and other parameters 
    call nml_read(path_par,"control","time_init",    time_init)                 ! [yr] Starting time
    call nml_read(path_par,"control","time_end",     time_end)                  ! [yr] Ending time
    call nml_read(path_par,"control","time_equil",   time_equil)                ! [yr] Years to equilibrate first
    call nml_read(path_par,"control","dtt",          dtt)                       ! [yr] Main loop time step 
    call nml_read(path_par,"control","dt1D_out",     dt1D_out)                  ! [yr] Frequency of 1D output 
    call nml_read(path_par,"control","dt2D_out",     dt2D_out)                  ! [yr] Frequency of 2D output 

    ! Assume program is running from the output folder
    outfldr = "./"

    ! Define input and output locations 
    path_const = trim(outfldr)//"yelmo_const_Earth.nml"
    file1D     = trim(outfldr)//"yelmo1D.nc"
    file2D     = trim(outfldr)//"yelmo2D.nc"

    ! === Initialize ice sheet model =====

    ! General initialization of yelmo constants (used globally)
    call yelmo_global_init(path_const)

    ! Initialize data objects and load initial topography
    call yelmo_init(yelmo1,filename=path_par,grid_def="file",time=time_init)

    ! === Set initial boundary conditions for current time and yelmo state =====
    ! ybound: z_bed, z_sl, H_sed, H_w, smb, T_srf, bmb_shlf , Q_geo

    yelmo1%bnd%z_sl     = 0.0           ! [m]
    yelmo1%bnd%H_sed    = 0.0           ! [m]
    yelmo1%bnd%H_w      = 0.0           ! [m]
    yelmo1%bnd%Q_geo    = 50.0          ! [mW/m2]
    
    yelmo1%bnd%bmb_shlf = -10.0         ! [m.i.e./a]
    yelmo1%bnd%T_shlf   = T0            ! [K]   

    ! Impose present-day surface mass balance and present-day temperature field
    yelmo1%bnd%smb      = yelmo1%dta%pd%smb        ! [m.i.e./a]
    yelmo1%bnd%T_srf    = yelmo1%dta%pd%T_srf      ! [K]
    
    call yelmo_print_bound(yelmo1%bnd)

    time = time_init 

    ! Define channel field 
    allocate(channels(yelmo1%grd%nx,yelmo1%grd%ny))

    ! Update C_bed, if needed
    if (yelmo1%dyn%par%C_bed_method .eq. -1) then 
        call calc_ydyn_cbed_external(yelmo1%dyn,yelmo1%tpo,yelmo1%thrm,yelmo1%bnd,channels)
    end if 

    ! Initialize state variables (dyn,therm,mat)
    ! (initialize temps with robin method with a cold base)
    call yelmo_init_state(yelmo1,path_par,time=time_init,thrm_method="robin-cold")

    ! Update C_bed again, if needed
    if (yelmo1%dyn%par%C_bed_method .eq. -1) then 
        call calc_ydyn_cbed_external(yelmo1%dyn,yelmo1%tpo,yelmo1%thrm,yelmo1%bnd,channels)
    end if 
    
    ! Define no-ice mask from present-day data
    allocate(mask_noice(yelmo1%grd%nx,yelmo1%grd%ny))
    mask_noice = .FALSE. 
    where(yelmo1%dta%pd%H_ice .le. 0.0) mask_noice = .TRUE. 

    ! Impose additional negative mass balance to no ice points 2 [m.i.e./a] melting
    where(mask_noice) yelmo1%bnd%smb = yelmo1%dta%pd%smb - 2.0 

    ! Impose a colder boundary temperature for equilibration step 
    ! -5 [K] for mimicking glacial times
    yelmo1%bnd%T_srf = yelmo1%dta%pd%T_srf - 10.0  

    ! Run yelmo for several years with constant boundary conditions and topo
    ! to equilibrate thermodynamics and dynamics
    call yelmo_update_equil(yelmo1,time,time_tot=time_equil,topo_fixed=.FALSE.,dt=1.0,ssa_vel_max=500.0)
    
    ! 2D file 
    call yelmo_write_init(yelmo1,file2D,time_init=time,units="years")  
    call write_step_2D(yelmo1,file2D,time=time)
    
    ! 1D file 
    call write_yreg_init(yelmo1,file1D,time_init=time_init,units="years",mask=yelmo1%bnd%ice_allowed)
    call write_yreg_step(yelmo1%reg,file1D,time=time) 
    

    ! Advance timesteps
    do n = 1, ceiling((time_end-time_init)/dtt)

        ! Get current time 
        time = time_init + n*dtt

        ! Update C_bed if needed
        if (yelmo1%dyn%par%C_bed_method .eq. -1) then 
            call calc_ydyn_cbed_external(yelmo1%dyn,yelmo1%tpo,yelmo1%thrm,yelmo1%bnd,channels)
        end if 
        
        ! Update temperature and smb as needed in time (ISMIP6)
        if (time .ge. -10e6 .and. time .lt. -10e3) then 
            ! Glacial period, impose cold climate 
            yelmo1%bnd%T_srf = yelmo1%dta%pd%T_srf - 10.0 

        else if (time .ge. -10e3 .and. time .lt. -8e3) then
            ! Holocene optimum 
            yelmo1%bnd%T_srf = yelmo1%dta%pd%T_srf + 1.0 

        else if (time .ge. -8e3) then 
            ! Entering Holocene, impose present-day temperatures 
            yelmo1%bnd%T_srf = yelmo1%dta%pd%T_srf
        end if 

        ! Update ice sheet 
        call yelmo_update(yelmo1,time)

        ! == MODEL OUTPUT =======================================================

        if (mod(time,dt2D_out)==0) then 
            call write_step_2D(yelmo1,file2D,time=time)
        end if 

        if (mod(time,dt1D_out)==0) then 
            call write_yreg_step(yelmo1%reg,file1D,time=time) 
        end if 

        if (mod(time,10.0)==0) then
            write(*,"(a,f14.4)") "yelmo::       time = ", time
        end if 

    end do 
    ! == Finished time loop == 


    ! Finalize program
    call yelmo_end(yelmo1,time=time)

    ! Stop timing 
    call cpu_time(cpu_end_time)


    write(*,"(a,f12.3,a)") "Time  = ",(cpu_end_time-cpu_start_time)/60.0 ," min"
    write(*,"(a,f12.1,a)") "Speed = ",(1e-3*(time_end-time_init))/((cpu_end_time-cpu_start_time)/3600.0), " kiloyears / hr"

contains

    subroutine write_step_2D(ylmo,filename,time)

        implicit none 
        
        type(yelmo_class), intent(IN) :: ylmo
        character(len=*),  intent(IN) :: filename
        real(prec), intent(IN) :: time

        ! Local variables
        integer    :: ncid, n
        real(prec) :: time_prev 

        real(prec) :: uxy_rmse, H_rmse 
        real(prec), allocatable :: tmp(:,:) 

        allocate(tmp(ylmo%grd%nx,ylmo%grd%ny))

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_size(filename,"time",ncid)
        call nc_read(filename,"time",time_prev,start=[n],count=[1],ncid=ncid) 
        if (abs(time-time_prev).gt.1e-5) n = n+1 

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        ! Write model speed 
        call nc_write(filename,"speed",ylmo%par%model_speed,units="kyr/hr",long_name="Model speed (Yelmo only)", &
                      dim1="time",start=[n],count=[1],ncid=ncid)
        
        ! initmip specific error metrics 
        tmp = ylmo%tpo%now%H_ice-ylmo%dta%pd%H_ice
        if (count(tmp .ne. 0.0) .gt. 0) then 
            H_rmse = sqrt(sum(tmp**2)/count(tmp .ne. 0.0))
        else 
            H_rmse = mv 
        end if 

        tmp = ylmo%dyn%now%uxy_s-ylmo%dta%pd%uxy_s
        if (count(tmp .ne. 0.0) .gt. 0) then 
            uxy_rmse = sqrt(sum(tmp**2)/count(tmp .ne. 0.0))
        else
            uxy_rmse = mv
        end if 

        call nc_write(filename,"rmse_H",H_rmse,units="m",long_name="RMSE - Ice thickness", &
                      dim1="time",start=[n],count=[1],ncid=ncid)
        call nc_write(filename,"rmse_uxy",uxy_rmse,units="m/a",long_name="RMSE - Surface velocity", &
                      dim1="time",start=[n],count=[1],ncid=ncid)
        
        ! == ISMIP6 specific variables 
        ! http://www.climate-cryosphere.org/wiki/index.php?title=InitMIP-Greenland#Appendix_2_.E2.80.93_Naming_conventions.2C_upload_and_model_output_data.

        call nc_write(filename,"H_ice",ylmo%tpo%now%H_ice,units="m",long_name="Ice thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf",ylmo%tpo%now%z_srf,units="m",long_name="Surface elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_bed",ylmo%bnd%z_bed,units="m",long_name="Bedrock elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid) 
        call nc_write(filename,"Q_geo",ylmo%bnd%Q_geo,units="mW/m^2",long_name="Geothermal heat flux", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"smb",ylmo%bnd%smb,units="m/a ice equiv.",long_name="Surface mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"bmb",ylmo%tpo%now%bmb,units="m/a ice equiv.",long_name="Basal mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dHicedt",ylmo%tpo%now%dHicedt,units="m/a",long_name="Ice thickness change", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"ux_s",ylmo%dyn%now%ux_s,units="m/a",long_name="Surface velocity (x)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy_s",ylmo%dyn%now%uy_s,units="m/a",long_name="Surface velocity (y)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uz_s",ylmo%dyn%now%uz(:,:,ylmo%par%nz_ac),units="m/a",long_name="Surface velocity (z)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"ux_b",ylmo%dyn%now%ux_b,units="m/a",long_name="Basal velocity (x)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy_b",ylmo%dyn%now%uy_b,units="m/a",long_name="Basal velocity (y)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uz_b",ylmo%dyn%now%uz(:,:,1),units="m/a",long_name="Basal velocity (z)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"ux_bar",ylmo%dyn%now%ux_bar,units="m/a",long_name="Vertically averaged velocity (x)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy_bar",ylmo%dyn%now%uy_bar,units="m/a",long_name="Vertically averaged velocity (y)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"T_ice_s",ylmo%thrm%now%T_ice(:,:,ylmo%par%nz_aa),units="K",long_name="Surface ice temperature", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"T_ice_b",ylmo%thrm%now%T_ice(:,:,1),units="K",long_name="Basal ice temperature", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"taub",ylmo%dyn%now%taub,units="Pa",long_name="Basal dragging stress (magnitude)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"calv",ylmo%tpo%now%calv,units="m/a ice equiv.",long_name="Calving rate", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_ice",ylmo%tpo%now%f_ice,units="1",long_name="Total ice fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_ice_g",ylmo%tpo%now%f_ice*ylmo%tpo%now%f_grnd,units="1",long_name="Grounded ice fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_ice_f",ylmo%tpo%now%f_ice*(1.0-ylmo%tpo%now%f_grnd),units="1",long_name="Floating ice fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! == Additional variables 

        call nc_write(filename,"mask_bed",ylmo%tpo%now%mask_bed,units="",long_name="Bed mask", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"uxy_s",ylmo%dyn%now%uxy_s,units="m/a",long_name="Surface velocity (magnitude)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_b",ylmo%dyn%now%uxy_b,units="m/a",long_name="Basal velocity (magnitude)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_bar",ylmo%dyn%now%uxy_bar,units="m/a",long_name="Vertically averaged velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)  
        
        call nc_write(filename,"mb_applied",ylmo%tpo%now%mb_applied,units="m",long_name="Applied net mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"N_eff",ylmo%tpo%now%N_eff,units="bar",long_name="Effective pressure", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"C_bed",ylmo%dyn%now%C_bed,units="m a^-1 Pa^-2",long_name="Dragging constant", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"beta",ylmo%dyn%now%beta,units="Pa a m-1",long_name="Dragging coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"visc_eff",ylmo%dyn%now%visc_eff,units="Pa a^-1",long_name="Vertically averaged viscosity", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"T_prime",ylmo%thrm%now%T_ice-ylmo%thrm%now%T_pmp,units="deg C",long_name="Homologous ice temperature", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"f_pmp",ylmo%thrm%now%f_pmp,units="1",long_name="Fraction of grid point at pmp", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! Ice thickness comparison with present-day 
        call nc_write(filename,"H_ice_errpd",ylmo%tpo%now%H_ice-ylmo%dta%pd%H_ice,units="m",long_name="Ice thickness error wrt present day", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_s_errpd",ylmo%dyn%now%uxy_s-ylmo%dta%pd%uxy_s,units="m",long_name="Ice thickness error wrt present day", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! Channels 
        call nc_write(filename,"channels",channels,units="m m-1",long_name="Channel diagnosis", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine write_step_2D

    subroutine calc_ydyn_cbed_external(dyn,tpo,thrm,bnd,channels)
        ! Update C_bed based on parameter choices

        implicit none
        
        type(ydyn_class),   intent(INOUT) :: dyn
        type(ytopo_class),  intent(IN)    :: tpo 
        type(ytherm_class), intent(IN)    :: thrm
        type(ybound_class), intent(IN)    :: bnd  
        real(prec),         intent(INOUT) :: channels(:,:) 

        integer :: i, j, nx, ny 
        integer :: i1, i2, j1, j2 
        real(prec), allocatable :: f_channel(:,:) 

        real(prec) :: channel_lim = 1e-6 

        nx = size(dyn%now%C_bed,1)
        ny = size(dyn%now%C_bed,2)
        
        allocate(f_channel(nx,ny)) 

        ! Set C_bed according to temperate character of base

        ! Smooth transition between temperate and frozen C_bed
        dyn%now%C_bed = (thrm%now%f_pmp)*dyn%par%cf_stream &
                    + (1.0_prec - thrm%now%f_pmp)*dyn%par%cf_frozen 

        if (dyn%par%streaming_margin) then 
            ! Ensure that both the margin points and the grounding line
            ! are always considered streaming, independent of their
            ! thermodynamic character (as sometimes these can incorrectly become frozen)

        
            ! Ensure any marginal point is also treated as streaming 
            do j = 1, ny 
            do i = 1, nx 

                i1 = max(i-1,1)
                i2 = min(i+1,nx)
                j1 = max(j-1,1)
                j2 = min(j+1,ny)

                if (tpo%now%H_ice(i,j) .gt. 0.0 .and. &
                    (tpo%now%H_ice(i1,j) .le. 0.0 .or. &
                     tpo%now%H_ice(i2,j) .le. 0.0 .or. &
                     tpo%now%H_ice(i,j1) .le. 0.0 .or. &
                     tpo%now%H_ice(i,j2) .le. 0.0)) then 

                    dyn%now%C_bed(i,j) = dyn%par%cf_stream

                end if 

            end do 
            end do 

            ! Also ensure that grounding line is also considered streaming
            where(tpo%now%is_grline) dyn%now%C_bed = dyn%par%cf_stream

        end if 

        ! == Until here, C_bed is defined as normally with C_bed_method=1,
        !    now refine to increase only marginal velocities 

        ! Reduce C_bed further for low elevation points
        !where(tpo%now%z_srf .lt. 1500.0) dyn%now%C_bed = 0.5*dyn%now%C_bed

        ! Next diagnose channels
        call calc_channels(channels,tpo%now%z_srf,dyn%now%ux_bar,dyn%now%uy_bar,tpo%par%dx)

        ! Finally scale C_bed according to concavity of channels 
        !f_channel = exp(-channels/channel_lim)
        !where(f_channel .lt. 0.1) f_channel = 0.1 
        !where(f_channel .gt. 2.0) f_channel = 2.0  

        f_channel = 1.0 

        dyn%now%C_bed = dyn%now%C_bed * f_channel 
        
        return 

    end subroutine calc_ydyn_cbed_external


    subroutine calc_channels(channels,z_bed,ux,uy,dx)

        implicit none 

        real(prec), intent(OUT) :: channels(:,:) 
        real(prec), intent(IN)  :: z_bed(:,:) 
        real(prec), intent(IN)  :: ux(:,:) 
        real(prec), intent(IN)  :: uy(:,:) 
        real(prec), intent(IN)  :: dx 

        ! Local variables 
        integer :: i, j, nx, ny 
        real(prec) :: ux_aa, uy_aa, uxy 
        real(prec) :: dzdx1, dzdx2 
        real(prec) :: dz2dx2, dz2dy2
        real(prec) :: theta   ! [rad] Angle of direction of ice flow
        real(prec) :: alpha   ! [rad] Angle of direction perpindicular to ice flow 

        ! Finite-difference coefficients
        real(prec), parameter :: fm2 = -1.0/12.0 
        real(prec), parameter :: fm1 =  4.0/3.0 
        real(prec), parameter :: f0  = -5.0/2.0 
        real(prec), parameter :: fp1 =  4.0/3.0
        real(prec), parameter :: fp2 = -1.0/12.0 
         
        
        
        nx = size(channels,1)
        ny = size(channels,2)

        ! Set channels to zero initially 
        channels = 0.0 

        ! Find channels based on change in elevation perpendicular to flow direction,
        ! then (to do!) negative component for along flow direction  
        do j = 3, ny-2 
        do i = 3, nx-2 

            ! Get velocity of current grid point 
            ux_aa = 0.5*(ux(i-1,j) + ux(i,j))
            uy_aa = 0.5*(uy(i,j-1) + uy(i,j))
            uxy   = sqrt(ux_aa**2 + uy_aa**2)

            if (uxy .gt. 0.0) then 

                ! Get direction perpindicular ice flow 
                alpha = atan2(uy_aa,ux_aa) - pi/2.0 

                ! Only modify areas with some velocity available 

                ! Calculate second-derivative in each direction (2nd order)
                dz2dx2 = (1.0/dx**2)*sum([fm2,fm1,f0,fp1,fp2]*z_bed(i-2:i+2,j))
                dz2dy2 = (1.0/dx**2)*sum([fm2,fm1,f0,fp1,fp2]*z_bed(i,j-2:j+2))
                
                ! Scale derivative in each direction to get approximate concavity in
                ! direction of interest 
                channels(i,j) = cos(alpha)*dz2dx2 + sin(alpha)*dz2dy2

!                 if (abs(ux_aa) .gt. abs(uy_aa)) then 
!                     ! Flow predominantly in x-direction

!                     dzdx1         = (z_bed(i,j)   - z_bed(i,j-1)) / dx 
!                     dzdx2         = (z_bed(i,j+1) - z_bed(i,j))   / dx 
!                     channels(i,j) = (dzdx2-dzdx1) / dx 

!                     !channels(i,j) = (0.5*(z_bed(i,j-1)+z_bed(i,j+1)) - z_bed(i,j)) / dx 

!                 else 
!                     ! Flow predominantly in y-direction 

!                     dzdx1         = (z_bed(i,j)   - z_bed(i-1,j)) / dx 
!                     dzdx2         = (z_bed(i+1,j) - z_bed(i,j))   / dx 
!                     channels(i,j) = (dzdx2-dzdx1) / dx 
                    
!                     !channels(i,j) = (0.5*(z_bed(i-1,j)+z_bed(i+1,j)) - z_bed(i,j)) / dx 

!                 end if 


            end if 

        end do 
        end do 

        return 

    end subroutine calc_channels

end program yelmo_test 



