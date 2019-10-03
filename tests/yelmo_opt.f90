

program yelmo_test

    use ncio 
    use yelmo 
    use yelmo_tools, only : gauss_values

    use gaussian_filter 
    use basal_hydro_simple 
    use basal_dragging

    implicit none 

    type(yelmo_class)      :: yelmo1
    type(yelmo_class)      :: yelmo_ref 
    type(hydro_class)      :: hyd1 
    type(hydro_class)      :: hyd_ref 

    character(len=256) :: outfldr, file1D, file2D, file_restart, domain 
    character(len=512) :: path_par, path_const  
    real(prec) :: time_init, time_end, time_equil, time, dtt, dt1D_out, dt2D_out
    integer    :: n
    real(4)    :: cpu_start_time, cpu_end_time 

    ! Parameters 
    real(prec) :: bmb_shlf_const, dT_ann, z_sl  

    ! Optimization variables 
    real(prec) :: time_iter, time_iter_0, time_iter_1, time_iter_2
    integer    :: q, qmax, qmax_topo_fixed, qmax_iter_length_1, qmax_iter_length_2
    real(prec) :: time_tune, time_tune_0, time_tune_1, time_tune_2  
    integer    :: qmax_tune_length_1, qmax_tune_length_2
    logical    :: topo_fixed  
    real(prec) :: cf_min, cf_max, cf_init  
    integer    :: opt_method 

    integer    :: topo_rel_n, n_now  
    integer    :: topo_rel_iter(2)
    real(prec) :: topo_rel_taus(2)

    real(prec), allocatable :: cf_ref(:,:) 
    real(prec), allocatable :: cf_ref_dot(:,:) 

    ! No-ice mask (to impose additional melting)
    logical, allocatable :: mask_noice(:,:)  

    ! Start timing 
    call cpu_time(cpu_start_time)

    ! Determine the parameter file from the command line 
    call yelmo_load_command_line_args(path_par)

    call nml_read(path_par,"control","bmb_shlf_const",  bmb_shlf_const)            ! [yr] Constant imposed bmb_shlf value
    call nml_read(path_par,"control","dT_ann",          dT_ann)                    ! [K] Temperature anomaly (atm)
    call nml_read(path_par,"control","z_sl",            z_sl)                      ! [m] Sea level relative to present-day
    
    ! Choose optimization method (1: error method, 2: ratio method) 
    opt_method = 1 

    ! Simulation parameters
    time_init           = 0.0       ! [yr] Starting time
    dtt                 = 2.0       ! [yr] Time step for time loop 
    dt2D_out            = 10.0      ! [yr] 2D output writing 

    topo_rel_n          = 2
    topo_rel_iter       = [5,10]
    topo_rel_taus       = [10.0,50.0]

    cf_init    = 0.2                    ! [--]
    cf_min     = 0.001                  ! [--] 
    cf_max     = 2.0                    ! [--]

    qmax                = 200       ! Total number of iterations
    time_iter           = 10.0      ! [yr] 

!     if (opt_method .eq. 1) then 
!         ! Error method 
!         qmax                = 200       ! Total number of iterations
!         time_iter           = 500.0     ! [yr] 

!         qmax_iter_length_1  = 10        ! 1st number of iterations at which iteration length should increase
!         time_iter_1         = 1000.0    ! [yr] 
        
!         qmax_iter_length_2  = 20        ! 1st number of iterations at which iteration length should increase
!         time_iter_2         = 2000.0    ! [yr] 
        
!     else 
!         ! Ratio method 
!         qmax                = 100       ! Total number of iterations
!         time_tune           = 20.0      ! [yr]
!         time_iter           = 200.0     ! [yr] 
        
!     end if 
    
    ! Not used right now:
!     qmax_topo_fixed     = 0         ! Number of initial iterations that should use topo_fixed=.TRUE. 
!     time_iter_0         =  50.0     ! [yr] 
!     time_iter_1         = 200.0     ! [yr] 
!     time_iter_2         = 500.0     ! [yr] 
!     qmax_iter_length_1  = 10        ! 1st number of iterations at which iteration length should increase
!     qmax_iter_length_2  = 50        ! 2nd number of iterations at which iteration length should increase
    
    ! Assume program is running from the output folder
    outfldr = "./"

    ! Define input and output locations 
    path_const   = trim(outfldr)//"yelmo_const_Earth.nml"
    file1D       = trim(outfldr)//"yelmo1D.nc"
    file2D       = trim(outfldr)//"yelmo2D.nc"
    file_restart = trim(outfldr)//"yelmo_restart.nc"

    ! === Initialize ice sheet model =====

    ! General initialization of yelmo constants (used globally)
    call yelmo_global_init(path_const)

    ! Initialize data objects and load initial topography
    call yelmo_init(yelmo1,filename=path_par,grid_def="file",time=time_init)

    ! Also intialize simple basal hydrology object
    call hydro_init(hyd1,filename=path_par,nx=yelmo1%grd%nx,ny=yelmo1%grd%ny)
    call hydro_init_state(hyd1,yelmo1%tpo%now%H_ice,yelmo1%tpo%now%f_grnd,time_init)

    ! === Set initial boundary conditions for current time and yelmo state =====
    ! ybound: z_bed, z_sl, H_sed, H_w, smb, T_srf, bmb_shlf , Q_geo

    yelmo1%bnd%z_sl     = 0.0               ! [m]
    yelmo1%bnd%H_sed    = 0.0               ! [m]
    yelmo1%bnd%H_w      = hyd1%now%H_w      ! [m]
    yelmo1%bnd%Q_geo    = 50.0              ! [mW/m2]
    
    yelmo1%bnd%bmb_shlf = bmb_shlf_const    ! [m.i.e./a]
    yelmo1%bnd%T_shlf   = T0                ! [K]   

    if (dT_ann .lt. 0.0) yelmo1%bnd%T_shlf   = T0 + dT_ann*0.25_prec  ! [K] Oceanic temp anomaly
    
    ! Impose present-day surface mass balance and present-day temperature field
    yelmo1%bnd%smb      = yelmo1%dta%pd%smb             ! [m.i.e./a]
    yelmo1%bnd%T_srf    = yelmo1%dta%pd%T_srf + dT_ann  ! [K]
    
    call yelmo_print_bound(yelmo1%bnd)

    ! Define no-ice mask from present-day data
    allocate(mask_noice(yelmo1%grd%nx,yelmo1%grd%ny))
    mask_noice = .FALSE. 
    where(yelmo1%dta%pd%H_ice .le. 0.0) mask_noice = .TRUE. 

    ! Impose additional negative mass balance to no ice points of 2 [m.i.e./a] melting
    if (trim(yelmo1%par%domain) .eq. "Greenland") then 
        where(mask_noice) yelmo1%bnd%smb = yelmo1%dta%pd%smb - 2.0 
    else ! Antarctica
        where(mask_noice) yelmo1%bnd%bmb_shlf = yelmo1%bnd%bmb_shlf - 2.0 
    end if 
    
    ! Initialize cf_ref and calculate initial guess of C_bed 
    ! (requires ad-hoc initialization of N_eff too)
    allocate(cf_ref_dot(yelmo1%grd%nx,yelmo1%grd%ny))
    allocate(cf_ref(yelmo1%grd%nx,yelmo1%grd%ny))
    cf_ref_dot = 0.0 
    cf_ref     = cf_init 
    
    ! Initialize state variables (dyn,therm,mat)
    ! (initialize temps with robin method with a cold base),
    ! or from restart file, if specified 
    call yelmo_init_state(yelmo1,path_par,time=time_init,thrm_method="robin-cold")

    if (.not. yelmo1%par%use_restart) then 
        ! Run initialization steps 

        ! Note: From now on, using yelmo_update_equil_external allows for running with 
        ! interactive hydrology via hyd1 object and for passing cf_ref to be able
        ! to update C_bed as a function of N_eff interactively.

        ! ============================================================================================
        ! Step 1: Relaxtion step: run SIA model for 100 years to smooth out the input
        ! topography that will be used as a target. 

        call yelmo_update_equil_external(yelmo1,hyd1,cf_ref,time_init,time_tot=50.0,topo_fixed=.FALSE.,dt=1.0,ssa_vel_max=0.0)

        ! Define present topo as present-day dataset for comparison 
        yelmo1%dta%pd%H_ice = yelmo1%tpo%now%H_ice 
        yelmo1%dta%pd%z_srf = yelmo1%tpo%now%z_srf 

        ! ============================================================================================
        ! Step 2: Run the model for several ka in hybrid mode with topo_fixed to
        ! spin up the thermodynamics and have a reference state to reset.
        ! Store the reference state for future use.
        
        call yelmo_update_equil_external(yelmo1,hyd1,cf_ref,time_init,time_tot=20e3,topo_fixed=.TRUE.,dt=5.0,ssa_vel_max=0.0)
        call yelmo_update_equil_external(yelmo1,hyd1,cf_ref,time_init,time_tot=10e3, topo_fixed=.TRUE.,dt=1.0,ssa_vel_max=5000.0)

        ! Store the reference state
        yelmo_ref = yelmo1 
        hyd_ref   = hyd1 

        ! Initialize the 2D output file and write the initial model state 
        call yelmo_write_init(yelmo1,file2D,time_init,units="years")  
        call write_step_2D_opt(yelmo1,file2D,time_init,cf_ref,cf_ref_dot,mask_noice)  
        
        ! Initialize time variable 
        time = time_init 

        ! Write a restart file 
        call yelmo_restart_write(yelmo1,file_restart,time)
        stop "**** Done ****"

    end if 

    write(*,*) "Starting optimization..."

if (opt_method .eq. 1) then 
    ! Error method (Pollard and De Conto, 2012)

    n_now = 1 

    do q = 1, qmax 

        if (q .gt. topo_rel_iter(n_now)) then 
            ! Update relaxation parameters 
            n_now = n_now + 1 
            if (n_now .gt. topo_rel_n) then 
                ! Disable relaxation 
                yelmo1%tpo%par%topo_rel = 0 
            else
                yelmo1%tpo%par%topo_rel_tau = topo_rel_taus(n_now)
            end if 

            write(*,*) "relaxation: ", q, n_now, yelmo1%tpo%par%topo_rel, yelmo1%tpo%par%topo_rel_tau
        end if 

        ! Perform iteration loop to diagnose error for modifying C_bed 
        do n = 1, int(time_iter/dtt)
        
            time = time + dtt 

            ! Update ice sheet 
            call yelmo_update(yelmo1,time)

            ! Update C_bed 
            call calc_ydyn_cbed_external(yelmo1%dyn,yelmo1%tpo,yelmo1%thrm,yelmo1%bnd,yelmo1%grd, &
                                                                        domain,mask_noice,cf_ref)

            cf_ref_dot = 0.0_prec 

!             if (mod(nint(time*100),nint(dt2D_out*100))==0) then
!                 call write_step_2D_opt(yelmo1,file2D,time,cf_ref,cf_ref_dot,mask_noice)
!             end if 

        end do 

        ! Update cf_ref based on error metric(s) 
        call update_cf_ref_thickness_simple(cf_ref,cf_ref_dot,yelmo1%tpo%now%H_ice, &
                        yelmo1%bnd%z_bed,yelmo1%dyn%now%ux_bar,yelmo1%dyn%now%uy_bar, &
                        yelmo1%dta%pd%H_ice,yelmo1%tpo%par%dx,cf_min,cf_max)

        ! Update C_bed 
        call calc_ydyn_cbed_external(yelmo1%dyn,yelmo1%tpo,yelmo1%thrm,yelmo1%bnd,yelmo1%grd, &
                                                                        domain,mask_noice,cf_ref)

!         if (q .le. qmax_iter_length_2) then 
!             ! Reset model to the initial state (including H_w) and time, with updated C_bed field 
!             yelmo_ref%dyn%now%C_bed = yelmo1%dyn%now%C_bed 
!             yelmo1 = yelmo_ref 
!             hyd1   = hyd_ref 
!             time   = 0.0 
!             call yelmo_set_time(yelmo1,time) 
!         end if 
        
        
        ! Write the current solution 
        call write_step_2D_opt(yelmo1,file2D,time,cf_ref,cf_ref_dot,mask_noice)
        
    end do 

else 
    ! Ratio method (Le clec’h et al, 2019)

    do q = 1, qmax 

        ! Reset model to the initial state (including H_w) and time, with updated C_bed field 
        yelmo_ref%dyn%now%C_bed = yelmo1%dyn%now%C_bed 
        yelmo1 = yelmo_ref 
        hyd1   = hyd_ref 
        time   = 0.0 
        call yelmo_set_time(yelmo1,time) 
        
        ! Perform C_bed tuning step 
        do n = 1, int(time_tune)
        
            time = time + 1.0

            ! Update ice sheet 
            call yelmo_update(yelmo1,time)

            ! Update C_bed based on error metric(s) 
            call update_cf_ref_thickness_ratio(cf_ref,cf_ref_dot,yelmo1%tpo%now%H_ice, &
                            yelmo1%bnd%z_bed,yelmo1%dyn%now%ux_bar,yelmo1%dyn%now%uy_bar, &
                            yelmo1%dyn%now%uxy_i_bar,yelmo1%dyn%now%uxy_b,yelmo1%dta%pd%H_ice, &
                            yelmo1%tpo%par%dx,cf_min,cf_max=cf_max)

            ! Update C_bed 
            call calc_ydyn_cbed_external(yelmo1%dyn,yelmo1%tpo,yelmo1%thrm,yelmo1%bnd,yelmo1%grd, &
                                                                        domain,mask_noice,cf_ref)

        end do 

        ! Perform iteration loop to diagnose error for modifying C_bed 
        do n = 1, int(time_iter)
        
            time = time + 1.0

            ! Update ice sheet 
            call yelmo_update(yelmo1,time)

        end do 

        ! Write the current solution 
        call write_step_2D_opt(yelmo1,file2D,time,cf_ref,cf_ref_dot,mask_noice)
        
    end do 

end if 

!         call yelmo_update_equil_external(yelmo1,hyd1,cf_ref,time,time_tot=time_iter,topo_fixed=topo_fixed,dt=0.5,ssa_vel_max=5000.0)


    ! Finalize program
    call yelmo_end(yelmo1,time=time)

    ! Stop timing 
    call cpu_time(cpu_end_time)

    write(*,"(a,f12.3,a)") "Time  = ",(cpu_end_time-cpu_start_time)/60.0 ," min"
    write(*,"(a,f12.1,a)") "Speed = ",(1e-3*(time_end-time_init))/((cpu_end_time-cpu_start_time)/3600.0), " kiloyears / hr"

contains

    subroutine write_step_2D_opt(ylmo,filename,time,cf_ref,cf_ref_dot,mask_noice)

        implicit none 
        
        type(yelmo_class), intent(IN) :: ylmo
        character(len=*),  intent(IN) :: filename
        real(prec), intent(IN) :: time
        real(prec), intent(IN) :: cf_ref(:,:) 
        real(prec), intent(IN) :: cf_ref_dot(:,:)
        logical,    intent(IN) :: mask_noice(:,:) 

        ! Local variables
        integer    :: ncid, n
        real(prec) :: time_prev 

        real(prec) :: uxy_rmse, H_rmse, zsrf_rmse, loguxy_rmse 
        real(prec) :: rmse, err  
        real(prec), allocatable :: tmp(:,:) 
        real(prec), allocatable :: tmp1(:,:) 
        
        allocate(tmp(ylmo%grd%nx,ylmo%grd%ny))
        allocate(tmp1(ylmo%grd%nx,ylmo%grd%ny))

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
        
        ! 1D variables 
        call nc_write(filename,"V_ice",ylmo%reg%V_ice,units="km3",long_name="Ice volume", &
                              dim1="time",start=[n],count=[1],ncid=ncid)
        call nc_write(filename,"A_ice",ylmo%reg%A_ice,units="km2",long_name="Ice area", &
                              dim1="time",start=[n],count=[1],ncid=ncid)
        call nc_write(filename,"dVicedt",ylmo%reg%dVicedt,units="km3 yr-1",long_name="Rate of volume change", &
                      dim1="time",start=[n],count=[1],ncid=ncid)
        
        ! Ice limiting mask
        call nc_write(filename,"mask_noice",mask_noice,units="",long_name="No-ice mask", &
                      dim1="xc",dim2="yc",ncid=ncid)
        
        ! Variables
        call nc_write(filename,"H_ice",ylmo%tpo%now%H_ice,units="m",long_name="Ice thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf",ylmo%tpo%now%z_srf,units="m",long_name="Surface elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_bed",ylmo%tpo%now%mask_bed,units="",long_name="Bed mask", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"dHicedt",ylmo%tpo%now%dHicedt,units="m/a",long_name="Ice thickness change", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"C_bed",ylmo%dyn%now%C_bed,units="Pa",long_name="Bed friction coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"cf_ref",cf_ref,units="",long_name="Bed friction scalar", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"cf_ref_dot",cf_ref_dot,units="1/a",long_name="Bed friction scalar rate of change", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"N_eff",ylmo%dyn%now%N_eff,units="Pa",long_name="Effective pressure", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"beta",ylmo%dyn%now%beta,units="Pa a m-1",long_name="Basal friction coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"visc_eff",ylmo%dyn%now%visc_eff,units="Pa a m",long_name="Effective viscosity (SSA)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"uxy_s",ylmo%dyn%now%uxy_s,units="m/a",long_name="Surface velocity magnitude", &
                     dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_b",ylmo%dyn%now%uxy_b,units="m/a",long_name="Basal sliding velocity magnitude", &
                     dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"f_pmp",ylmo%thrm%now%f_pmp,units="",long_name="Basal temperate fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"T_prime",ylmo%thrm%now%T_ice-ylmo%thrm%now%T_pmp,units="deg C",long_name="Homologous ice temperature", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        
        call nc_write(filename,"ATT",ylmo%mat%now%ATT,units="a^-1 Pa^-3",long_name="Rate factor", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        
        ! Boundary variables (forcing)
        call nc_write(filename,"z_bed",ylmo%bnd%z_bed,units="m",long_name="Bedrock elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"smb",ylmo%bnd%smb,units="m a-1",long_name="Annual surface mass balance (ice equiv.)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"T_srf",ylmo%bnd%T_srf,units="K",long_name="Surface temperature", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"H_w",ylmo%bnd%H_w,units="m",long_name="Basal water layer", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! Target data (not time dependent)
        if (n .eq. 1) then 
            call nc_write(filename,"pd_z_srf",ylmo%dta%pd%z_srf,units="m",long_name="Observed surface elevation (present day)", &
                          dim1="xc",dim2="yc",ncid=ncid)
            call nc_write(filename,"pd_H_ice",ylmo%dta%pd%H_ice,units="m",long_name="Observed ice thickness (present day)", &
                          dim1="xc",dim2="yc",ncid=ncid)
            call nc_write(filename,"pd_uxy_s",ylmo%dta%pd%uxy_s,units="m",long_name="Observed surface velocity (present day)", &
                          dim1="xc",dim2="yc",ncid=ncid)
        end if 

        ! Error fields compared to targets
        call nc_write(filename,"pd_err_z_srf",ylmo%dta%pd%err_z_srf,units="m",long_name="Surface elevation error (present day)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"pd_err_uxy_s",ylmo%dta%pd%err_uxy_s,units="m",long_name="Surface velocity error (present day)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        tmp = ylmo%dta%pd%err_z_srf
        call filter_gaussian(var=tmp,sigma=2.0*ylmo%tpo%par%dx,dx=ylmo%tpo%par%dx, &
                                mask=ylmo%dta%pd%err_z_srf .ne. 0.0)
        
        call nc_write(filename,"pd_err_z_srf_sm",tmp,units="m",long_name="Smooth surface elevation error (present day)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! initmip specific error metrics 
        tmp = ylmo%dta%pd%err_H_ice
        if (n .gt. 1 .or. count(tmp .ne. 0.0) .gt. 0) then 
            H_rmse = sqrt(sum(tmp**2)/count(tmp .ne. 0.0))
        else 
            H_rmse = mv 
        end if 

        ! surface elevation too 
        tmp = ylmo%dta%pd%err_z_srf
        if (n .gt. 1 .or. count(tmp .ne. 0.0) .gt. 0) then 
            zsrf_rmse = sqrt(sum(tmp**2)/count(tmp .ne. 0.0))
        else 
            zsrf_rmse = mv 
        end if 

        tmp = ylmo%dta%pd%err_uxy_s
        if (n .gt. 1 .or. count(tmp .ne. 0.0) .gt. 0) then 
            uxy_rmse = sqrt(sum(tmp**2)/count(tmp .ne. 0.0))
        else
            uxy_rmse = mv
        end if 

        tmp = ylmo%dta%pd%uxy_s 
        where(ylmo%dta%pd%uxy_s .gt. 0.0) tmp = log(tmp)
        tmp1 = ylmo%dyn%now%uxy_s 
        where(ylmo%dyn%now%uxy_s .gt. 0.0) tmp1 = log(tmp1)
        
        if (n .gt. 1 .or. count(tmp1-tmp .ne. 0.0) .gt. 0) then 
            loguxy_rmse = sqrt(sum((tmp1-tmp)**2)/count(tmp1-tmp .ne. 0.0))
        else
            loguxy_rmse = mv
        end if 
        
        call nc_write(filename,"rmse_H",H_rmse,units="m",long_name="RMSE - Ice thickness", &
                      dim1="time",start=[n],count=[1],ncid=ncid,missing_value=mv)
        call nc_write(filename,"rmse_zsrf",zsrf_rmse,units="m",long_name="RMSE - Surface elevation", &
                      dim1="time",start=[n],count=[1],ncid=ncid,missing_value=mv)
        call nc_write(filename,"rmse_uxy",uxy_rmse,units="m/a",long_name="RMSE - Surface velocity", &
                      dim1="time",start=[n],count=[1],ncid=ncid,missing_value=mv)
        call nc_write(filename,"rmse_uxy_log",loguxy_rmse,units="log(m/a)",long_name="RMSE - Log surface velocity", &
                      dim1="time",start=[n],count=[1],ncid=ncid,missing_value=mv)
        
        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine write_step_2D_opt

    subroutine update_cf_ref_thickness_simple(cf_ref,cf_ref_dot,H_ice,z_bed,ux,uy,H_obs,dx,cf_min,cf_max)

        implicit none 

        real(prec), intent(INOUT) :: cf_ref(:,:) 
        real(prec), intent(INOUT) :: cf_ref_dot(:,:) 
        real(prec), intent(IN)    :: H_ice(:,:) 
        real(prec), intent(IN)    :: z_bed(:,:) 
        real(prec), intent(IN)    :: ux(:,:) 
        real(prec), intent(IN)    :: uy(:,:) 
        real(prec), intent(IN)    :: H_obs(:,:) 
        real(prec), intent(IN)    :: dx 
        real(prec), intent(IN)    :: cf_min 
        real(prec), intent(IN)    :: cf_max

        ! Local variables 
        integer :: i, j, nx, ny, i1, j1  
        real(prec) :: dx_km, f_dz, f_dz_lim, H_scale, f_scale   
        real(prec) :: ux_aa, uy_aa
        real(prec) :: H_ice_now, H_obs_now 

        real(prec), allocatable   :: cf_prev(:,:) 
        real(prec) :: wts0(5,5), wts(5,5) 

        nx = size(cf_ref,1)
        ny = size(cf_ref,2) 

        dx_km = dx*1e-3  
        
        allocate(cf_prev(nx,ny))

        ! Optimization parameters 
        H_scale  = 1000.0           ! [m] 
        f_dz_lim = 1.5              ! [--] 

        ! Get Gaussian weights 
        wts0 = gauss_values(dx_km,dx_km,sigma=dx_km*1.5,n=5)

        ! Store initial cf_ref solution 
        cf_prev = cf_ref 

        do j = 3, ny-2 
        do i = 3, nx-2 

            if ( abs(H_ice(i,j) - H_obs(i,j)) .ne. 0.0) then 
                ! Update where thickness error exists

                ! Determine downstream node

                ux_aa = 0.5*(ux(i,j)+ux(i+1,j))
                uy_aa = 0.5*(uy(i,j)+uy(i,j+1))
                
                if ( abs(ux_aa) .gt. abs(uy_aa) ) then 
                    ! Downstream in x-direction 
                    j1 = j 
                    if (ux_aa .lt. 0.0) then 
                        i1 = i-1 
                    else
                        i1 = i+1
                    end if 

                else 
                    ! Downstream in y-direction 
                    i1 = i 
                    if (uy_aa .lt. 0.0) then 
                        j1 = j-1
                    else
                        j1 = j+1
                    end if 

                end if 

                ! Get current ice thickness and obs ice thickness (smoothed)

                wts = wts0 
                where( H_ice(i-2:i+2,j-2:j+2) .eq. 0.0) wts = 0.0 
                call wtd_mean(H_ice_now,H_ice(i-2:i+2,j-2:j+2),wts) 

                wts = wts0 
                where( H_obs(i-2:i+2,j-2:j+2) .eq. 0.0) wts = 0.0 
                call wtd_mean(H_obs_now,H_obs(i-2:i+2,j-2:j+2),wts) 
                
                ! Get adjustment rate given error in z_srf
                f_dz = (H_ice_now - H_obs_now) / H_scale
                f_dz = max(f_dz,-f_dz_lim)
                f_dz = min(f_dz,f_dz_lim)
                
                f_scale = 10.0**(-f_dz) 

                cf_ref(i1,j1) = cf_prev(i1,j1)*f_scale

            end if 

        end do 
        end do 

        ! Ensure cf_ref is not below lower or upper limit 
        where (cf_ref .lt. cf_min) cf_ref = cf_min 
        where (cf_ref .gt. cf_max) cf_ref = cf_max 

        ! Additionally, apply a Gaussian filter to cf_ref to ensure smooth transitions
        call filter_gaussian(var=cf_ref,sigma=dx_km*0.2,dx=dx_km)     !,mask=err_z_srf .ne. 0.0)
        
        ! Also where no ice exists, set cf_ref = cf_min 
        where(H_obs .eq. 0.0) cf_ref = cf_min 

        ! Diagnose current rate of change of C_bed 
        cf_ref_dot = cf_ref - cf_prev

        return 

    end subroutine update_cf_ref_thickness_simple

    subroutine update_cf_ref_thickness_ratio(cf_ref,cf_ref_dot,H_ice,z_bed,ux,uy,uxy_i,uxy_b,H_obs,dx,cf_min,cf_max)

        implicit none 

        real(prec), intent(INOUT) :: cf_ref(:,:) 
        real(prec), intent(INOUT) :: cf_ref_dot(:,:) 
        real(prec), intent(IN)    :: H_ice(:,:) 
        real(prec), intent(IN)    :: z_bed(:,:) 
        real(prec), intent(IN)    :: ux(:,:)        ! Depth-averaged velocity (ux_bar)
        real(prec), intent(IN)    :: uy(:,:)        ! Depth-averaged velocity (uy_bar)
        real(prec), intent(IN)    :: uxy_i(:,:)     ! Internal shear velocity magnitude 
        real(prec), intent(IN)    :: uxy_b(:,:)     ! Basal sliding velocity magnitude 
        real(prec), intent(IN)    :: H_obs(:,:) 
        real(prec), intent(IN)    :: dx 
        real(prec), intent(IN)    :: cf_min 
        real(prec), intent(IN)    :: cf_max

        ! Local variables 
        integer :: i, j, nx, ny, i1, j1, n 
        real(prec) :: f_err, f_vel, f_corr, dx_km 
        real(prec) :: ux_aa, uy_aa 
        real(prec) :: H_ice_now, H_obs_now 

        real(prec), allocatable   :: cf_prev(:,:) 
        real(prec) :: wts0(5,5), wts(5,5) 

        real(prec),parameter :: exp1 = 2.0

        nx = size(cf_ref,1)
        ny = size(cf_ref,2) 

        dx_km = dx*1e-3  
        
        allocate(cf_prev(nx,ny))

        ! Get Gaussian weights 
        wts0 = gauss_values(dx_km,dx_km,sigma=dx_km*1.5,n=5)

!         do i = 1, 5 
!         write(*,*) wts0(i,:) 
!         end do 
!         stop 

        ! Store initial cf_ref solution 
        cf_prev = cf_ref 

        do j = 3, ny-2 
        do i = 3, nx-2 

            if ( abs(H_ice(i,j) - H_obs(i,j)) .ne. 0.0) then 
                ! Update where thickness error exists

                ! Determine downstream point to apply changes

                ux_aa = 0.5*(ux(i,j)+ux(i+1,j))
                uy_aa = 0.5*(uy(i,j)+uy(i,j+1))
                
                if ( abs(ux_aa) .gt. abs(uy_aa) ) then 
                    ! Downstream in x-direction 
                    j1 = j 
                    if (ux_aa .lt. 0.0) then 
                        i1 = i-1 
                    else
                        i1 = i+1
                    end if 

                else 
                    ! Downstream in y-direction 
                    i1 = i 
                    if (uy_aa .lt. 0.0) then 
                        j1 = j-1
                    else
                        j1 = j+1
                    end if 

                end if 

                ! Calculate thickness error ratio 
!                 f_err = H_ice(i,j) / max(H_obs(i,j),1e-1)
                
                wts = wts0 
                !where( H_ice(i-2:i+2,j-2:j+2) .eq. 0.0) wts = 0.0 
                call wtd_mean(H_ice_now,H_ice(i-2:i+2,j-2:j+2),wts) 

                wts = wts0 
                !where( H_obs(i-2:i+2,j-2:j+2) .eq. 0.0) wts = 0.0 
                call wtd_mean(H_obs_now,H_obs(i-2:i+2,j-2:j+2),wts) 
                
!                 n = count(H_ice(i-1:i+1,j-1:j+1).gt.0.0)
!                 if (n .gt. 0) then
!                     H_ice_now = sum(H_ice(i-1:i+1,j-1:j+1),mask=H_ice(i-1:i+1,j-1:j+1).gt.0.0) / real(n,prec)
!                 else 
!                     H_ice_now = 0.0 
!                 end if 

!                 n = count(H_obs(i-1:i+1,j-1:j+1).gt.0.0)
!                 if (n .gt. 0) then
!                     H_obs_now = sum(H_obs(i-1:i+1,j-1:j+1),mask=H_obs(i-1:i+1,j-1:j+1).gt.0.0) / real(n,prec)
!                 else 
!                     H_obs_now = 0.0 
!                 end if 
                
                f_err = ( H_ice_now / max(H_obs_now,1e-1) )
                
                ! Calculate ratio of deformational velocity to sliding velocity
                f_vel = uxy_i(i,j) / max(uxy_b(i,j),1e-1) 

                ! Calculate correction factor (beta_old / beta_new)
                f_corr = ( max( f_err + f_vel*(f_err-1.0_prec), 1e-1) )**exp1

                ! Apply correction to update cf_ref
                cf_ref(i1,j1) = cf_prev(i1,j1) * f_corr**(-1.0)

            end if 

        end do 
        end do 

        ! Ensure cf_ref is not below lower or upper limit 
        where (cf_ref .lt. cf_min) cf_ref = cf_min 
        where (cf_ref .gt. cf_max) cf_ref = cf_max 

        ! Additionally, apply a Gaussian filter to cf_ref to ensure smooth transitions
        call filter_gaussian(var=cf_ref,sigma=dx_km*0.25,dx=dx_km)     !,mask=err_z_srf .ne. 0.0)
        
        ! Also where no ice exists, set cf_ref = cf_min 
        where(H_ice .eq. 0.0) cf_ref = cf_min 

        ! Diagnose current rate of change of cf_ref 
        cf_ref_dot = cf_ref - cf_prev

        return 

    end subroutine update_cf_ref_thickness_ratio

    subroutine calc_ydyn_cbed_external(dyn,tpo,thrm,bnd,grd,domain,mask_noice,cf_ref)
        ! Update C_bed [Pa] based on parameter choices

        implicit none
        
        type(ydyn_class),   intent(INOUT) :: dyn
        type(ytopo_class),  intent(IN)    :: tpo 
        type(ytherm_class), intent(IN)    :: thrm
        type(ybound_class), intent(IN)    :: bnd  
        type(ygrid_class),  intent(IN)    :: grd
        character(len=*),   intent(IN)    :: domain 
        logical,            intent(IN)    :: mask_noice(:,:) 
        real(prec),         intent(IN)    :: cf_ref(:,:) 

        integer :: i, j, nx, ny 
        integer :: i1, i2, j1, j2 
        real(prec) :: f_scale 
        
        real(prec), allocatable :: lambda_bed(:,:)  

        nx = size(dyn%now%C_bed,1)
        ny = size(dyn%now%C_bed,2)
        
        allocate(lambda_bed(nx,ny))

            ! cf_ref [unitless] is obtained as input to this routine from optimization 

            ! =============================================================================
            ! Step 2: calculate lambda functions to scale C_bed from default value 
            
            select case(trim(dyn%par%cb_scale))

                case("lin_zb")
                    ! Linear scaling function with bedrock elevation
                    
                    lambda_bed = calc_lambda_bed_lin(bnd%z_bed,dyn%par%cb_z0,dyn%par%cb_z1)

                case("exp_zb")
                    ! Exponential scaling function with bedrock elevation
                    
                    ! Default
                    lambda_bed = calc_lambda_bed_exp(bnd%z_bed,dyn%par%cb_z0,dyn%par%cb_z1)

                    ! Modifications - increased friction in Wilkes Land (South)
                    where (bnd%basins .ge. 12 .and. &
                           bnd%basins .le. 17) lambda_bed = calc_lambda_bed_exp(bnd%z_bed,-400.0,dyn%par%cb_z1)

                case("till_const")
                    ! Constant till friction angle

                    lambda_bed = calc_lambda_till_const(dyn%par%till_phi_const)

                case("till_zb")
                    ! Linear till friction angle versus elevation

                    lambda_bed = calc_lambda_till_linear(bnd%z_bed,bnd%z_sl,dyn%par%till_phi_min,dyn%par%till_phi_max, &
                                                            dyn%par%till_phi_zmin,dyn%par%till_phi_zmax)

                case DEFAULT
                    ! No scaling

                    lambda_bed = 1.0

            end select 

            ! Set lambda_bed to lower limit for regions of noice 
            where (mask_noice) lambda_bed = dyn%par%cb_min 
            
            ! Ensure lambda_bed is not below lower limit [default range 0:1] 
            where (lambda_bed .lt. dyn%par%cb_min) lambda_bed = dyn%par%cb_min

            ! =============================================================================
            ! Step 3: calculate C_bed [Pa]
            
            dyn%now%C_bed = (cf_ref*lambda_bed)*dyn%now%N_eff 

        return 

    end subroutine calc_ydyn_cbed_external

    subroutine yelmo_update_equil_external(dom,hyd,cf_ref,time,time_tot,dt,topo_fixed,ssa_vel_max)
        ! Iterate yelmo solutions to equilibrate without updating boundary conditions

        type(yelmo_class), intent(INOUT) :: dom
        type(hydro_class), intent(INOUT) :: hyd 
        real(prec), intent(IN) :: cf_ref(:,:) 
        real(prec), intent(IN) :: time            ! [yr] Current time
        real(prec), intent(IN) :: time_tot        ! [yr] Equilibration time 
        real(prec), intent(IN) :: dt              ! Local dt to be used for all modules
        logical,    intent(IN) :: topo_fixed      ! Should topography be fixed? 
        real(prec), intent(IN) :: ssa_vel_max     ! Local vel limit to be used, if == 0.0, no ssa used

        ! Local variables 
        real(prec) :: time_now  
        integer    :: n, nstep 
        logical    :: use_ssa         ! Should ssa be active?  
        logical    :: dom_topo_fixed
        logical    :: dom_use_ssa 
        real(prec) :: dom_dtmax 
        integer    :: dom_ntt 
        real(prec) :: dom_ssa_vel_max 

        ! Only run equilibration if time_tot > 0 

        if (time_tot .gt. 0.0) then 

            ! Consistency check
            use_ssa = .FALSE. 
            if (ssa_vel_max .gt. 0.0) use_ssa = .TRUE. 

            ! Save original model choices 
            dom_topo_fixed  = dom%tpo%par%topo_fixed 
            dom_use_ssa     = dom%dyn%par%use_ssa 
            dom_dtmax       = dom%par%dtmax
            dom_ntt         = dom%par%ntt 
            dom_ssa_vel_max = dom%dyn%par%ssa_vel_max

            ! Set model choices equal to equilibration choices 
            dom%tpo%par%topo_fixed  = topo_fixed 
            dom%dyn%par%use_ssa     = use_ssa 
            dom%par%dtmax           = dt 
            dom%par%ntt             = 1 
            dom%dyn%par%ssa_vel_max = ssa_vel_max

            write(*,*) 
            write(*,*) "Starting equilibration steps, time to run [yrs]: ", time_tot 

            do n = 1, ceiling(time_tot/dt)

                time_now = time + n*dt
                call yelmo_update(dom,time_now)

                ! Update basal hydrology 
                call hydro_update(hyd,dom%tpo%now%H_ice,dom%tpo%now%f_grnd, &
                            -dom%thrm%now%bmb_grnd*rho_ice/rho_w,time_now)

                ! Pass updated hydrology variable to Yelmo boundary field
                dom%bnd%H_w = hyd%now%H_w 

                ! Update C_bed
                call calc_ydyn_cbed_external(yelmo1%dyn,yelmo1%tpo,yelmo1%thrm,yelmo1%bnd,yelmo1%grd, &
                                                                                domain,mask_noice,cf_ref)

            end do

            ! Restore original model choices 
            dom%tpo%par%topo_fixed  = dom_topo_fixed 
            dom%dyn%par%use_ssa     = dom_use_ssa 
            dom%par%dtmax           = dom_dtmax 
            dom%par%ntt             = dom_ntt 
            dom%dyn%par%ssa_vel_max = dom_ssa_vel_max

            write(*,*) 
            write(*,*) "Equilibration complete."
            write(*,*) 

        end if 

        ! Reset model time back to input time 
        dom%tpo%par%time      = time 
        dom%thrm%par%time     = time 

        hyd%now%time          = time 

        return

    end subroutine yelmo_update_equil_external
    

    subroutine wtd_mean(var_ave,var,wts)
        ! wts == gauss_values(dx,dy,sigma,n)

        implicit none

        real(prec), intent(OUT) :: var_ave 
        real(prec), intent(IN)  :: var(:,:) 
        real(prec), intent(IN)  :: wts(:,:) 

        ! Local variables 
        real(prec) :: wts_tot 
        real(prec) :: wts_norm(size(wts,1),size(wts,2))

        wts_tot = sum(wts) 
        if (wts_tot .gt. 0.0) then 
            wts_norm = wts / wts_tot 
        else 
            wts_norm = 0.0 
        end if 

        var_ave = sum(var*wts_norm) 

        return 

    end subroutine wtd_mean

    ! Extra...

    subroutine update_C_bed_thickness(C_bed,dCbed,phi,err_z_srf,H_ice,z_bed,ux,uy,dx,phi_min,phi_max, &
                        cf_stream,cb_z0,cb_z1,cb_min)

        implicit none 

        real(prec), intent(INOUT) :: C_bed(:,:) 
        real(prec), intent(INOUT) :: dCbed(:,:) 
        real(prec), intent(INOUT) :: phi(:,:) 
        real(prec), intent(IN)    :: err_z_srf(:,:) 
        real(prec), intent(IN)    :: H_ice(:,:) 
        real(prec), intent(IN)    :: z_bed(:,:) 
        real(prec), intent(IN)    :: ux(:,:) 
        real(prec), intent(IN)    :: uy(:,:) 
        real(prec), intent(IN)    :: dx 
        real(prec), intent(IN)    :: phi_min 
        real(prec), intent(IN)    :: phi_max 
        real(prec), intent(IN)    :: cf_stream 
        real(prec), intent(IN)    :: cb_z0
        real(prec), intent(IN)    :: cb_z1 
        real(prec), intent(IN)    :: cb_min 

        ! Local variables 
        integer :: i, j, nx, ny, i1, j1  
        real(prec) :: dphi, dx_km, f_dz 
        real(prec) :: ux_aa, uy_aa 
        real(prec) :: zsrf_rmse 

        real(prec), allocatable   :: C_bed_prev(:,:) 

        real(prec) :: dphi_min  
        real(prec) :: dphi_max 
        real(prec) :: err_z_fac 

        nx = size(C_bed,1)
        ny = size(C_bed,2) 

        allocate(C_bed_prev(nx,ny))

        ! Optimization parameters 
        dphi_min  = -0.5       ! [degrees] maximum rate of change (negative)
        dphi_max  =  1.0       ! [degrees] maximum rate of change (positive)
        
        ! Store initial C_bed solution 
        C_bed_prev = C_bed 

        ! Calculate global rmse error 
        if (count(err_z_srf .ne. 0.0) .gt. 0) then 
            zsrf_rmse = sqrt(sum(err_z_srf**2)/count(err_z_srf .ne. 0.0))
        else 
            zsrf_rmse = 0.0 
        end if 

        if (zsrf_rmse .gt. 90.0) then
            ! Maintain a faster scale
            err_z_fac = 100.0      ! [m]       Elevation-error scale for maximum
        else
            ! Slow down the optimization to reduce waves, if we are near the solution
            err_z_fac = 200.0 
        end if 

        do j = 3, ny-2 
        do i = 3, nx-2 

            if (err_z_srf(i,j) .ne. 0.0) then 
                ! Update where elevation error exists

                ! Get adjustment rate given error in z_srf 
                f_dz = -err_z_srf(i,j) / err_z_fac 
                f_dz = max(f_dz, dphi_min)
                f_dz = min(f_dz, dphi_max)
                dphi = f_dz 
                
                ! 1. Apply change at current point 
if (.FALSE.) then 
                phi(i,j)  = phi(i,j) + dphi 
                phi(i,j)  = max(phi(i,j),phi_min)
                phi(i,j)  = min(phi(i,j),phi_max)

                C_bed(i,j) = cf_stream*tan(phi(i,j)*pi/180.0)
end if 

                ! 2. Apply change downstream (this may overlap with other changes)

                ux_aa = 0.5*(ux(i,j)+ux(i+1,j))
                uy_aa = 0.5*(uy(i,j)+uy(i,j+1))
                
                if ( abs(ux_aa) .gt. abs(uy_aa) ) then 
                    ! Downstream in x-direction 
                    j1 = j 
                    if (ux_aa .lt. 0.0) then 
                        i1 = i-1 
                    else
                        i1 = i+1
                    end if 

                else 
                    ! Downstream in y-direction 
                    i1 = i 
                    if (uy_aa .lt. 0.0) then 
                        j1 = j-1
                    else
                        j1 = j+1
                    end if 

                end if 

                phi(i1,j1)  = phi(i1,j1) + dphi 
                phi(i1,j1)  = max(phi(i1,j1),phi_min)
                phi(i1,j1)  = min(phi(i1,j1),phi_max)

                C_bed(i1,j1) = cf_stream*tan(phi(i1,j1)*pi/180.0)

            end if 

        end do 
        end do 

        ! Ensure C_bed is not below lower limit 
        where (C_bed .lt. cb_min) C_bed = cb_min 

        ! Additionally, apply a Gaussian filter to C_bed to ensure smooth transitions
!         dx_km = dx*1e-3  
!         call filter_gaussian(var=C_bed,sigma=64.0,dx=dx_km)     !,mask=err_z_srf .ne. 0.0)
        
        ! Diagnose current rate of change of C_bed 
        dCbed = C_bed - C_bed_prev

        return 

    end subroutine update_C_bed_thickness

    subroutine guess_C_bed(C_bed,phi,uxy_s,phi_min,phi_max,cf_stream)

        implicit none 

        real(prec), intent(INOUT) :: C_bed(:,:) 
        real(prec), intent(INOUT) :: phi(:,:) 
        real(prec), intent(IN)    :: uxy_s(:,:) 
        real(prec), intent(IN)    :: phi_min
        real(prec), intent(IN)    :: phi_max
        real(prec), intent(IN)    :: cf_stream 

        ! Local variables 
        integer :: i, j, nx, ny 
        real(prec) :: logvel, logvel_max, f_scale    

        nx = size(C_bed,1)
        ny = size(C_bed,2)

        !logvel_max = log(maxval(uxy_s))
        logvel_max = log(100.0) 

        do j = 1, ny 
        do i = 1, nx 

            if (uxy_s(i,j) .gt. 0.0) then 
                ! Calculate expected till angle versus velocity 

                logvel   = max(0.0,log(uxy_s(i,j)))
                f_scale  = logvel / logvel_max
                if (f_scale .gt. 1.0) f_scale = 1.0 
                phi(i,j) = 0.5*phi_max - f_scale*(0.5*phi_max-2.0*phi_min)

            else 
                ! Set phi to the minimum 

                phi(i,j) = phi_min 

            end if 

            ! Calculate C_bed following till friction approach (Bueler and van Pelt, 2015)
            C_bed(i,j) = cf_stream*tan(phi(i,j)*pi/180)

        end do 
        end do

        return 

    end subroutine guess_C_bed

    subroutine calc_ydyn_cbed_external_channels(dyn,tpo,thrm,bnd,channels)
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

        if (dyn%par%cb_margin_pmp) then 
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

        ! == Until here, C_bed is defined as normally with cb_method=1,
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

    end subroutine calc_ydyn_cbed_external_channels


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



