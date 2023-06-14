

program yelmo_test

    use nml
    use ncio 
    use yelmo 

    use ice_optimization 

    implicit none 

    type(yelmo_class)      :: yelmo1
    type(yelmo_class)      :: yelmo_ref 

    character(len=256) :: outfldr, file1D, file2D, file_restart_init, file_restart, domain 
    character(len=512) :: path_par  
    real(prec) :: time_init, time_end, time_equil, time, dtt, dt1D_out, dt2D_out
    integer    :: n, q, n_now

    ! Parameters 
    real(prec) :: bmb_shlf_const, dT_ann, z_sl  

    ! Optimization variables  
    character(len=12) :: opt_method 
    real(prec) :: cf_min
    real(prec) :: cf_max
    real(prec) :: cf_init  
    
    integer    :: n_iter 
    real(prec) :: time_iter
    real(prec) :: time_iter_therm 
    real(prec) :: time_steady_end
    real(prec) :: time_tune 
    real(prec) :: time_iter_tot 

    real(prec) :: rel_time1, rel_time2, rel_tau1, rel_tau2, rel_m  
    real(prec) :: scale_ftime, scale_time1, scale_time2, scale_err1, scale_err2 
    real(prec) :: sigma_err, sigma_vel 

    real(prec) :: tau, err_scale 
    real(prec) :: tau_c 
    real(prec) :: H0 

    character(len=12) :: optvar 
    logical           :: reset_model
    character(len=56) :: cb_ref_init_method

    real(prec), allocatable :: mb_corr(:,:) 

    ! No-ice mask (to impose additional melting)
    logical, allocatable :: mask_noice(:,:)  

    real(8) :: cpu_start_time, cpu_end_time, cpu_dtime  
    
    ! Start timing 
    call yelmo_cpu_time(cpu_start_time)
    
    ! Determine the parameter file from the command line 
    call yelmo_load_command_line_args(path_par)

    ! ====================================================
    ! Load general parameters 

    call nml_read(path_par,"ctrl","opt_method", opt_method)                 ! Total number of iterations
    call nml_read(path_par,"ctrl","cb_ref_init_method", cb_ref_init_method)  ! How should cb_ref be initialized (guess, restart, none)
    call nml_read(path_par,"ctrl","sigma_err",       sigma_err)              ! [--] Smoothing radius for error to calculate correction in cb_ref (in multiples of dx)
    call nml_read(path_par,"ctrl","sigma_vel",       sigma_vel)              ! [m/a] Speed at which smoothing diminishes to zero
    call nml_read(path_par,"ctrl","cf_min",          cf_min)                 ! [--] Minimum allowed cf value 
    call nml_read(path_par,"ctrl","cf_max",          cf_max)                 ! [--] Maximum allowed cf value 
    
    call nml_read(path_par,"ctrl","bmb_shlf_const",  bmb_shlf_const)         ! [yr] Constant imposed bmb_shlf value
    call nml_read(path_par,"ctrl","dT_ann",          dT_ann)                 ! [K] Temperature anomaly (atm)
    call nml_read(path_par,"ctrl","z_sl",            z_sl)                   ! [m] Sea level relative to present-day

    ! ====================================================
    ! Load optimization-method specific parameters 

    select case(trim(opt_method))

        case("P12") 
            ! Pollard and DeConto (2012) 

            call nml_read(path_par,"opt_P12","time_init",       time_init)                 ! Total number of iterations
            call nml_read(path_par,"opt_P12","n_iter",          n_iter)                 ! Total number of iterations
            call nml_read(path_par,"opt_P12","time_iter",       time_iter)              ! [yr] Time for each iteration
            call nml_read(path_par,"opt_P12","time_iter_therm", time_iter_therm)        ! [yr] Time to run thermodynamics for each iteration
            call nml_read(path_par,"opt_P12","time_steady_end", time_steady_end)        ! [yr] Time for each iteration
            call nml_read(path_par,"opt_P12","reset_model",     reset_model)            ! Reset model to reference state between iterations?
            
            optvar = "ice"   ! ice/vel
            call nml_read(path_par,"opt_P12","rel_tau1",    rel_tau1)             ! [yr] Initial relaxation tau, fixed until rel_time1 
            call nml_read(path_par,"opt_P12","rel_tau2",    rel_tau2)             ! [yr] Final tau, reached at rel_time2, when relaxation disabled 
            call nml_read(path_par,"opt_P12","scale_ftime", scale_ftime)          ! [-]  Fraction of time_iter_tot at which to start transition from scale_err1 to scale_err2
            call nml_read(path_par,"opt_P12","scale_err1",  scale_err1)           ! [m]  Initial value for err_scale parameter in cb_ref optimization 
            call nml_read(path_par,"opt_P12","scale_err2",  scale_err2)           ! [m]  Final value for err_scale parameter reached at scale_time2   

            ! Optimization parameters 
            rel_time1           = time_iter_tot*0.4 ! [yr] Time to begin reducing tau from tau1 to tau2 
            rel_time2           = time_iter_tot*0.8 ! [yr] Time to reach tau2, and to disable relaxation 
            rel_m               = 2.0               ! [--] Non-linear exponent to scale interpolation between time1 and time2 

            scale_time1         = time_iter_tot*scale_ftime ! [yr] Time to begin increasing err_scale from scale_err1 to scale_err2 
            scale_time2         = time_iter_tot*1.0         ! [yr] Time to reach scale_H2 

            ! Determine total iteration time 
            time_iter_tot = time_iter*(n_iter-1)
            dtt                 = 5.0               ! [yr] Time step for time loop 
            dt2D_out            = time_iter         ! [yr] 2D output writing 
            cf_init             = 0.2               ! [--] Initial cf value everywhere (not too important)

            ! === Consistency checks =============================
            
            ! Ensure error scaling only increases with time
            scale_err2 = max(scale_err1,scale_err2)
            
            ! Ensure relaxation constant only increases with time 
            rel_tau2 = max(rel_tau1,rel_tau2)
            
        case("L21")
            ! Lipscomb et al. (2021) 

            call nml_read(path_par,"opt_L21","time_init",       time_init)   
            call nml_read(path_par,"opt_L21","time_end",        time_end)    
            
            call nml_read(path_par,"opt_L21","rel_tau1",        rel_tau1)             
            call nml_read(path_par,"opt_L21","rel_tau2",        rel_tau2)             
            call nml_read(path_par,"opt_L21","rel_time1",       rel_time1)           
            call nml_read(path_par,"opt_L21","rel_time2",       rel_time2)           
            
            call nml_read(path_par,"opt_L21","tau_c",  tau_c)    ! [yr] L21: Relaxation time scale for cb_ref adjustment 
            call nml_read(path_par,"opt_L21","H0",     H0)       ! [m]  L21: Error scaling

            dtt                 = 5.0               ! [yr] Time step for time loop 
            dt2D_out            = 100.0             ! [yr] 2D output writing 
            cf_init             = 0.2               ! [--] Initial cf value everywhere (not too important)

        case("L19")
            ! Le clec’h et al. (2019) - needs revising 

            write(*,*) "This method is outdated and needs checking, mainly the parameter loading."
            stop 

            ! Not used:
            ! ! Ratio method 
            ! n_iter                = 100       ! Total number of iterations
            ! time_tune           = 20.0      ! [yr]
            ! time_iter           = 200.0     ! [yr] 
    
    end select

    ! ====================================================

    ! Assume program is running from the output folder
    outfldr = "./"

    ! Define input and output locations 
    file1D       = trim(outfldr)//"yelmo1D.nc"
    file2D       = trim(outfldr)//"yelmo2D.nc"
    
    file_restart_init = trim(outfldr)//"yelmo_restart_init.nc"
    file_restart      = trim(outfldr)//"yelmo_restart.nc"

    ! === Initialize ice sheet model =====
    
    ! Initialize data objects and load initial topography
    call yelmo_init(yelmo1,filename=path_par,grid_def="file",time=time_init)

    ! Initialize mass balance correction matrix 
    allocate(mb_corr(yelmo1%grd%nx,yelmo1%grd%ny))
    mb_corr = 0.0_prec 

    ! Define no-ice mask from present-day data
    allocate(mask_noice(yelmo1%grd%nx,yelmo1%grd%ny))
    mask_noice = .FALSE. 
    !where(yelmo1%dta%pd%H_ice .le. 0.0) mask_noice = .TRUE. 

    ! === Set initial boundary conditions for current time and yelmo state =====
    ! ybound: z_bed, z_sl, H_sed, H_w, smb, T_srf, bmb_shlf , Q_geo

    yelmo1%bnd%z_sl     = 0.0                   ! [m]
    yelmo1%bnd%H_sed    = 0.0                   ! [m]
    yelmo1%bnd%Q_geo    = 50.0                  ! [mW/m2]
    
    yelmo1%bnd%bmb_shlf = bmb_shlf_const        ! [m.i.e./a]
    yelmo1%bnd%T_shlf   = T0 + dT_ann*0.25_prec ! [K]   

    ! Impose present-day surface temperature and surface mass balance fields
    yelmo1%bnd%T_srf    = yelmo1%dta%pd%T_srf + dT_ann  ! [K]
    yelmo1%bnd%smb      = yelmo1%dta%pd%smb             ! [m.i.e./a]
    
    ! Impose additional negative mass balance to no ice points of 2 [m.i.e./a] melting
    if (trim(yelmo1%par%domain) .eq. "Greenland") then 
        where(mask_noice) yelmo1%bnd%smb = yelmo1%dta%pd%smb - 2.0 
    else ! Antarctica
        where(mask_noice) yelmo1%bnd%bmb_shlf = yelmo1%bnd%bmb_shlf - 2.0 
    end if 
    
    call yelmo_print_bound(yelmo1%bnd)

    
    ! Initialize state variables (dyn,therm,mat)
    ! (initialize temps with robin method with a cold base),
    ! or from restart file, if specified 
    call yelmo_init_state(yelmo1,time=time_init,thrm_method="robin-cold")


    select case(trim(cb_ref_init_method))

        case("guess")

            ! Calculate new initial guess of cb_ref using info from dyn
            call guess_cb_ref(yelmo1%dyn%now%cb_ref,yelmo1%dyn%now%taud,yelmo1%dta%pd%uxy_s, &
                                yelmo1%dta%pd%H_ice,yelmo1%dta%pd%H_grnd,yelmo1%dyn%par%beta_u0,cf_min,cf_max)

            ! Update ice sheet to get everything in sync
            call yelmo_update_equil(yelmo1,time_init,time_tot=1.0,dt=1.0_prec,topo_fixed=.TRUE.)

        case("restart")

            ! Pass, cb_ref obtained from restart file 
            if (.not. yelmo1%par%use_restart) then 
                write(*,*) "yelmo_opt:: Error: cb_ref_init_method='restart' can only be used &
                &in conjunction with a restart file being loaded."
                stop 
            end if 

        case DEFAULT  ! "none"

            yelmo1%dyn%now%cb_ref = cf_init 
    
    end select  

    if (.not. yelmo1%par%use_restart) then 
        ! Run initialization steps 

        ! ============================================================================================
        ! Step 1: Relaxtion step: run DIVA model for a short time to smooth out the input
        ! topography that will be used as a target. 

        call yelmo_update_equil(yelmo1,time_init,time_tot=20.0_prec,dt=1.0_prec,topo_fixed=.FALSE.,dyn_solver="diva")

        ! Define present topo as present-day dataset for comparison 
        yelmo1%dta%pd%H_ice = yelmo1%tpo%now%H_ice 
        yelmo1%dta%pd%z_srf = yelmo1%tpo%now%z_srf 

        ! ============================================================================================
        ! Step 2: Run the model for several ka in standard mode with topo_fixed to
        ! spin up the thermodynamics and have a reference state to reset.
        ! Store the reference state for future use.
        
        call yelmo_update_equil(yelmo1,time_init,time_tot=20e3_prec,dt=5.0_prec,topo_fixed=.TRUE.)

        ! Write a restart file 
        call yelmo_restart_write(yelmo1,file_restart_init,time_init)
!         stop "**** Done ****"
        
    else 
        ! If restart file was used, redefine boundary conditions here since they 
        ! may have been overwritten with other information.

        ! === Set initial boundary conditions for current time and yelmo state =====
        ! ybound: z_bed, z_sl, H_sed, H_w, smb, T_srf, bmb_shlf , Q_geo

        yelmo1%bnd%z_sl     = 0.0                   ! [m]
        yelmo1%bnd%H_sed    = 0.0                   ! [m]
        yelmo1%bnd%Q_geo    = 50.0                  ! [mW/m2]
        
        yelmo1%bnd%bmb_shlf = bmb_shlf_const        ! [m.i.e./a]
        yelmo1%bnd%T_shlf   = T0 + dT_ann*0.25_prec ! [K]   

        ! Impose present-day surface temperature and surface mass balance fields
        yelmo1%bnd%T_srf    = yelmo1%dta%pd%T_srf + dT_ann  ! [K]
        yelmo1%bnd%smb      = yelmo1%dta%pd%smb             ! [m.i.e./a]
        
        ! Impose additional negative mass balance to no ice points of 2 [m.i.e./a] melting
        if (trim(yelmo1%par%domain) .eq. "Greenland") then 
            where(mask_noice) yelmo1%bnd%smb = yelmo1%dta%pd%smb - 2.0 
        else ! Antarctica
            where(mask_noice) yelmo1%bnd%bmb_shlf = yelmo1%bnd%bmb_shlf - 2.0 
        end if 
        
        call yelmo_print_bound(yelmo1%bnd)

    end if 

    
    ! Store the reference state
    yelmo_ref    = yelmo1
    
    ! Store initial optimization parameter choices 
    tau       = rel_tau1 
    err_scale = scale_err1 

    ! Initialize the 2D output file and write the initial model state 
    call yelmo_write_init(yelmo1,file2D,time_init,units="years")  
    call write_step_2D_opt(yelmo1,file2D,time_init,mb_corr,mask_noice,tau,err_scale)  
    
    write(*,*) "Starting optimization..."


    select case(opt_method)

        case("P12")
            ! Pollard and DeConto (2012) 

            ! Initialize timing variables 
            time = time_init 
            n_now = 1 

            do q = 1, n_iter 

                ! === Optimization parameters =========
                
                tau       = get_opt_param(time,time1=rel_time1,time2=rel_time2,p1=rel_tau1,p2=rel_tau2,m=rel_m)
                err_scale = get_opt_param(time,time1=scale_time1,time2=scale_time2,p1=scale_err1,p2=scale_err2,m=1.0)
                
                ! Set model tau, and set yelmo relaxation switch (2: gl-line and shelves relaxing; 0: no relaxation)
                yelmo1%tpo%par%topo_rel_tau = tau 
                yelmo1%tpo%par%topo_rel     = 2
                if (time .gt. rel_time2) yelmo1%tpo%par%topo_rel = 0 

                ! Update relaxation parameters for reference model too 
                yelmo_ref%tpo%par%topo_rel_tau = yelmo1%tpo%par%topo_rel_tau
                yelmo_ref%tpo%par%topo_rel     = yelmo1%tpo%par%topo_rel

                ! Diagnose mass balance correction term 
                call update_mb_corr(mb_corr,yelmo1%tpo%now%H_ice,yelmo1%dta%pd%H_ice,tau)
                
                ! === Update cb_ref and reset model ===================

                if (q .gt. 1) then
                    ! Perform optimization after first iteration

                    ! Update cb_ref based on error metric(s) 
                    call update_cb_ref_errscaling(yelmo1%dyn%now%cb_ref,yelmo1%tpo%now%H_ice, &
                                        yelmo1%bnd%z_bed,yelmo1%dyn%now%ux_s,yelmo1%dyn%now%uy_s, &
                                        yelmo1%dta%pd%H_ice,yelmo1%dta%pd%uxy_s,yelmo1%dta%pd%H_grnd.le.0.0_prec, &
                                        yelmo1%tpo%par%dx,cf_min,cf_max,sigma_err,sigma_vel,err_scale, &
                                        fill_dist=80.0_prec,optvar=optvar)

                end if 
                
                if (reset_model) then
                    ! Reset model to reference state with updated cb_ref 

                    yelmo_ref%dyn%now%cb_ref = yelmo1%dyn%now%cb_ref
                    call yelmo_set_time(yelmo_ref,time)

                    if (time_iter_therm .gt. 0.0) then
                        ! Run thermodynamics without updating topography 
                        call yelmo_update_equil(yelmo_ref,time,time_tot=time_iter_therm,dt=5.0_prec,topo_fixed=.TRUE.)
                    end if 

                    yelmo1 = yelmo_ref 

                else 
                    ! Continue with the model as it is

                    if (time_iter_therm .gt. 0.0) then
                        ! Run thermodynamics without updating topography  
                        call yelmo_update_equil(yelmo1,time,time_tot=time_iter_therm,dt=5.0_prec,topo_fixed=.TRUE.)
                    end if 
                     
                end if 

                ! === Update time_iter ==================
                time_end = time_iter
                if (q .eq. n_iter) time_end = time_steady_end 

                write(*,"(a,i4,f10.1,i4,f10.1,f12.1,f10.1)") "iter_par: ", q, time, &
                                    yelmo1%tpo%par%topo_rel, tau, err_scale, time_end

                ! Perform iteration loop to diagnose error for modifying c_bed 
                do n = 1, int(time_end/dtt)
                
                    time = time + dtt 

                    ! Update ice sheet 
                    call yelmo_update(yelmo1,time)

                    if (mod(nint(time*100),nint(dt2D_out*100))==0) then
                        call write_step_2D_opt(yelmo1,file2D,time,mb_corr,mask_noice,tau,err_scale)
                    end if 

                end do 

            end do 

        case("L21") 
            ! Lipscomb et al. (2021)
            
            ! Perform transient simulation with cb_ref nudging  
            do n = 1, ceiling((time_end-time_init)/dtt)

                ! Get current time 
                time = time_init + n*dtt
                
                ! === Optimization parameters =========
                
                ! Update model relaxation time scale and error scaling (in [m])
                tau = get_opt_param(time,time1=rel_time1,time2=rel_time2,p1=rel_tau1,p2=rel_tau2,m=rel_m)
                
                ! Set model tau, and set yelmo relaxation switch (2: gl-line and shelves relaxing; 0: no relaxation)
                yelmo1%tpo%par%topo_rel_tau = tau 
                yelmo1%tpo%par%topo_rel     = 2
                if (time .gt. rel_time2) yelmo1%tpo%par%topo_rel = 0 

                ! Diagnose mass balance correction term 
                call update_mb_corr(mb_corr,yelmo1%tpo%now%H_ice,yelmo1%dta%pd%H_ice,tau)
                

                ! === Optimization update step =========

                ! Update cb_ref based on error metric(s) 
                call update_cb_ref_errscaling_l21(yelmo1%dyn%now%cb_ref,yelmo1%tpo%now%H_ice, &
                                    yelmo1%tpo%now%dHicedt,yelmo1%bnd%z_bed,yelmo1%bnd%z_sl,yelmo1%dyn%now%ux_s,yelmo1%dyn%now%uy_s, &
                                    yelmo1%dta%pd%H_ice,yelmo1%dta%pd%uxy_s,yelmo1%dta%pd%H_grnd.le.0.0_prec, &
                                    yelmo1%tpo%par%dx,cf_min,cf_max,sigma_err,sigma_vel,tau_c,H0, &
                                    fill_dist=80.0_prec,dt=dtt)

                ! === Advance ice sheet =========

                ! Update ice sheet 
                call yelmo_update(yelmo1,time)

                if (mod(nint(time*100),nint(dt2D_out*100))==0) then
                    call write_step_2D_opt(yelmo1,file2D,time,mb_corr,mask_noice,tau,err_scale)
                end if 

            end do 

        case("L19")
            ! Ratio method (Le clec’h et al, 2019) - needs revising

            write(*,*) "This method may not work now - check carefully."
            stop 

            ! Initialize timing variables 
            time = time_init 
            
            do q = 1, n_iter 

                ! Reset model to the initial state (including H_w) and time, with updated c_bed field 
                yelmo_ref%dyn%now%c_bed = yelmo1%dyn%now%c_bed 
                yelmo1 = yelmo_ref 
                time   = 0.0 
                call yelmo_set_time(yelmo1,time) 
                
                ! Perform c_bed tuning step 
                do n = 1, int(time_tune)
                
                    time = time + 1.0

                    ! Update ice sheet 
                    call yelmo_update(yelmo1,time)

                    ! Update c_bed based on error metric(s) 
                    call update_cb_ref_thickness_ratio(yelmo1%dyn%now%cb_ref,yelmo1%tpo%now%H_ice, &
                                    yelmo1%bnd%z_bed,yelmo1%dyn%now%ux_bar,yelmo1%dyn%now%uy_bar, &
                                    yelmo1%dyn%now%uxy_i_bar,yelmo1%dyn%now%uxy_b,yelmo1%dta%pd%H_ice, &
                                    yelmo1%tpo%par%dx,cf_min,cf_max=cf_max)

                end do 

                ! Perform iteration loop to diagnose error for modifying c_bed 
                do n = 1, int(time_iter)
                
                    time = time + 1.0

                    ! Update ice sheet 
                    call yelmo_update(yelmo1,time)

                end do 

                ! Write the current solution 
                call write_step_2D_opt(yelmo1,file2D,time,mb_corr,mask_noice,tau,err_scale)
                
            end do 

    end select 

    
    ! Write a final restart file 
    call yelmo_restart_write(yelmo1,file_restart,time)

    ! Finalize program
    call yelmo_end(yelmo1,time=time)

    ! Stop timing 
    call yelmo_cpu_time(cpu_end_time,cpu_start_time,cpu_dtime)
    
    write(*,"(a,f12.3,a)") "Time  = ",cpu_dtime/60.0 ," min"
    write(*,"(a,f12.1,a)") "Speed = ",(1e-3*(time_end-time_init))/(cpu_dtime/3600.0), " kiloyears / hr"
    
contains

    subroutine write_step_2D_opt(ylmo,filename,time,mb_corr,mask_noice,tau,err_scale)

        implicit none 
        
        type(yelmo_class), intent(IN) :: ylmo
        character(len=*),  intent(IN) :: filename
        real(prec), intent(IN) :: time
        real(prec), intent(IN) :: mb_corr(:,:)
        logical,    intent(IN) :: mask_noice(:,:) 
        real(prec), intent(IN) :: tau 
        real(prec), intent(IN) :: err_scale 

        ! Local variables
        integer    :: ncid, n
        real(prec) :: time_prev 

        real(prec) :: rmse, err  
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

        ! Write model metrics (model speed, dt, eta)
        call yelmo_write_step_model_metrics(filename,ylmo,n,ncid)

        ! Write present-day data metrics (rmse[H],etc)
        call yelmo_write_step_pd_metrics(filename,ylmo,n,ncid)
        
        ! 1D variables 
        call nc_write(filename,"V_ice",ylmo%reg%V_ice,units="km3",long_name="Ice volume", &
                              dim1="time",start=[n],count=[1],ncid=ncid)
        call nc_write(filename,"A_ice",ylmo%reg%A_ice,units="km2",long_name="Ice area", &
                              dim1="time",start=[n],count=[1],ncid=ncid)
        call nc_write(filename,"dVicedt",ylmo%reg%dVicedt,units="km3 yr-1",long_name="Rate of volume change", &
                      dim1="time",start=[n],count=[1],ncid=ncid)
        
        call nc_write(filename,"opt_tau",tau,units="yr",long_name="Relaxation time scale (ice shelves)", &
                      dim1="time",start=[n],count=[1],ncid=ncid)
        call nc_write(filename,"opt_err_scale",err_scale,units="m",long_name="Error scaling constant", &
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
        
        call nc_write(filename,"c_bed",ylmo%dyn%now%c_bed,units="Pa",long_name="Bed friction coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"cb_ref",yelmo1%dyn%now%cb_ref,units="",long_name="Bed friction scalar", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"N_eff",ylmo%dyn%now%N_eff,units="Pa",long_name="Effective pressure", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"beta",ylmo%dyn%now%beta,units="Pa a m-1",long_name="Basal friction coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"visc_eff_int",ylmo%dyn%now%visc_eff_int,units="Pa a m",long_name="Depth-integrated effective viscosity (SSA)", &
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
        
        call nc_write(filename,"H_w",ylmo%thrm%now%H_w,units="m",long_name="Basal water layer", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"T_prime_b",ylmo%thrm%now%T_prime_b,units="m",long_name="Basal homologous ice temperature", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! Boundary variables (forcing)
        call nc_write(filename,"z_bed",ylmo%bnd%z_bed,units="m",long_name="Bedrock elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"smb",ylmo%bnd%smb,units="m a-1",long_name="Annual surface mass balance (ice equiv.)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"T_srf",ylmo%bnd%T_srf,units="K",long_name="Surface temperature", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"mb_corr",mb_corr,units="m/a",long_name="SMB correction term", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! Target data (not time dependent)
        if (n .eq. 1) then 
            call nc_write(filename,"z_srf_pd",ylmo%dta%pd%z_srf,units="m",long_name="Observed surface elevation (present day)", &
                          dim1="xc",dim2="yc",ncid=ncid)
            call nc_write(filename,"H_ice_pd",ylmo%dta%pd%H_ice,units="m",long_name="Observed ice thickness (present day)", &
                          dim1="xc",dim2="yc",ncid=ncid)
            call nc_write(filename,"uxy_s_pd",ylmo%dta%pd%uxy_s,units="m",long_name="Observed surface velocity (present day)", &
                          dim1="xc",dim2="yc",ncid=ncid)
        end if 

        ! Comparison with present-day 
        call nc_write(filename,"H_ice_pd_err",ylmo%dta%pd%err_H_ice,units="m",long_name="Ice thickness error wrt present day", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf_pd_err",ylmo%dta%pd%err_z_srf,units="m",long_name="Surface elevation error wrt present day", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_s_pd_err",ylmo%dta%pd%err_uxy_s,units="m/a",long_name="Surface velocity error wrt present day", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine write_step_2D_opt


    ! Extra...

    subroutine calc_ydyn_cbed_external_channels(dyn,tpo,thrm,bnd,channels)
        ! Update c_bed based on parameter choices

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

        nx = size(dyn%now%c_bed,1)
        ny = size(dyn%now%c_bed,2)
        
        allocate(f_channel(nx,ny)) 

        ! Set c_bed according to temperate character of base

        ! Smooth transition between temperate and frozen c_bed
        dyn%now%c_bed = (thrm%now%f_pmp)*dyn%par%cf_stream &
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

                    dyn%now%c_bed(i,j) = dyn%par%cf_stream

                end if 

            end do 
            end do 

            ! Also ensure that grounding line is also considered streaming
            where(tpo%now%is_grline) dyn%now%c_bed = dyn%par%cf_stream

        end if 

        ! == Until here, c_bed is defined as normally with cb_method=1,
        !    now refine to increase only marginal velocities 

        ! Reduce c_bed further for low elevation points
        !where(tpo%now%z_srf .lt. 1500.0) dyn%now%c_bed = 0.5*dyn%now%c_bed

        ! Next diagnose channels
        call calc_channels(channels,tpo%now%z_srf,dyn%now%ux_bar,dyn%now%uy_bar,tpo%par%dx)

        ! Finally scale c_bed according to concavity of channels 
        !f_channel = exp(-channels/channel_lim)
        !where(f_channel .lt. 0.1) f_channel = 0.1 
        !where(f_channel .gt. 2.0) f_channel = 2.0  

        f_channel = 1.0 

        dyn%now%c_bed = dyn%now%c_bed * f_channel 
        
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
                ! Only modify areas with some velocity available 

                ! Get direction perpindicular ice flow 
                alpha = atan2(uy_aa,ux_aa) - pi/2.0 

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



