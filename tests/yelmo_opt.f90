

program yelmo_test

    use ncio 
    use yelmo 
    use yelmo_tools, only : gauss_values

    use gaussian_filter 
    use basal_dragging

    implicit none 

    type(yelmo_class)      :: yelmo1
    type(yelmo_class)      :: yelmo_ref 

    character(len=256) :: outfldr, file1D, file2D, file_restart_init, file_restart, domain 
    character(len=512) :: path_par, path_const  
    real(prec) :: time_init, time_end, time_equil, time, dtt, dt1D_out, dt2D_out
    integer    :: n, q, n_now
    real(4)    :: cpu_start_time, cpu_end_time 

    ! Parameters 
    real(prec) :: bmb_shlf_const, dT_ann, z_sl  

    ! Optimization variables  
    integer    :: opt_method 
    real(prec) :: cf_min, cf_max, cf_init  
    
    integer    :: qmax 
    real(prec) :: time_iter
    real(prec) :: time_steady
    real(prec) :: time_tune 

    integer    :: iter_steps(5) 
    integer    :: topo_rels(5) 
    real(prec) :: topo_rel_taus(5)
    real(prec) :: H_scales(5) 

    real(prec) :: rel_time1, rel_time2, rel_tau1, rel_tau2, rel_q  
    real(prec) :: scale_time1, scale_time2, scale_H1, scale_H2 

    real(prec) :: tau, H_scale 

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
    dtt                 = 5.0       ! [yr] Time step for time loop 
    dt2D_out            = 500.0     ! [yr] 2D output writing 

    qmax                = 51                ! Total number of iterations
    time_iter           = 500.0             ! [yr] Time for each iteration 
    time_steady         = 50e3              ! [yr] Time to run to steady state at the end without further optimization

    ! Optimization parameters 
    rel_time1           = 10e3      ! [yr] Time to begin reducing tau from tau1 to tau2 
    rel_time2           = 20e3      ! [yr] Time to reach tau2, and to disable relaxation 
    rel_tau1            = 10.0      ! [yr] Initial relaxation tau, fixed until rel_time1 
    rel_tau2            = 500.0     ! [yr] Final tau, reached at rel_time2, when relaxation disabled 
    rel_q               = 2.0       ! [--] Non-linear exponent to scale interpolation between time1 and time2 

    scale_time1         = 15e3      ! [yr] Time to begin increasing H_scale from scale_H1 to scale_H2 
    scale_time2         = 25e3      ! [yr] Time to reach scale_H2 
    scale_H1            = 1000.0    ! [m]  Initial value for H_scale parameter in cf_ref optimization 
    scale_H2            = 2000.0    ! [m]  Final value for H_scale parameter reached at scale_time2 


!     iter_steps          = [12,16,20,25,35]
!     topo_rels           = [1,1,0,0,0]
!     topo_rel_taus       = [10.0,100.0,1000.0,0.0,0.0]
!     H_scales            = [1000.0,1000.0,1000.0,2000.0,2000.0] 

    cf_init    = 0.2                        ! [--]
    cf_min     = 1e-5                       ! [--] 
    cf_max     = 1.0                        ! [--]


! Not used:
!         ! Ratio method 
!         qmax                = 100       ! Total number of iterations
!         time_tune           = 20.0      ! [yr]
!         time_iter           = 200.0     ! [yr] 
    
    ! Consistency checks 
    if (rel_time2 .ge. qmax*time_iter) then 
        write(*,*) "Error: rel_time2 >= total time. rel_time2 must be less than the total simulation &
                    &years, so that relaxation is disabled at the end of the simulation."
        stop 
    end if 

    ! Assume program is running from the output folder
    outfldr = "./"

    ! Define input and output locations 
    path_const   = trim(outfldr)//"yelmo_const_Earth.nml"
    file1D       = trim(outfldr)//"yelmo1D.nc"
    file2D       = trim(outfldr)//"yelmo2D.nc"
    
    file_restart_init = trim(outfldr)//"yelmo_restart_init.nc"
    file_restart      = trim(outfldr)//"yelmo_restart.nc"

    ! === Initialize ice sheet model =====

    ! General initialization of yelmo constants (used globally)
    call yelmo_global_init(path_const)

    ! Initialize data objects and load initial topography
    call yelmo_init(yelmo1,filename=path_par,grid_def="file",time=time_init)

    ! === Set initial boundary conditions for current time and yelmo state =====
    ! ybound: z_bed, z_sl, H_sed, H_w, smb, T_srf, bmb_shlf , Q_geo

    yelmo1%bnd%z_sl     = 0.0               ! [m]
    yelmo1%bnd%H_sed    = 0.0               ! [m]
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
    cf_ref_dot = 0.0 
    
    ! Initialize state variables (dyn,therm,mat)
    ! (initialize temps with robin method with a cold base),
    ! or from restart file, if specified 
    call yelmo_init_state(yelmo1,path_par,time=time_init,thrm_method="robin-cold")

if (.FALSE.) then 
    ! Calculate new initial guess of cf_ref using info from dyn
    call guess_cf_ref(yelmo1%dyn%now%cf_ref,yelmo1%dyn%now%taud,yelmo1%dta%pd%uxy_s, &
                        yelmo1%dta%pd%H_ice,yelmo1%dta%pd%H_grnd,yelmo1%dyn%par%beta_u0,cf_min,cf_max)

    ! Update ice sheet to get everything in sync
    call yelmo_update_equil(yelmo1,time_init,time_tot=1.0,topo_fixed=.TRUE.,dt=1.0,ssa_vel_max=5e3)

else 

    yelmo1%dyn%now%cf_ref = cf_init 
    
end if 

    if (.not. yelmo1%par%use_restart) then 
        ! Run initialization steps 

        ! ============================================================================================
        ! Step 1: Relaxtion step: run SIA model for a short time to smooth out the input
        ! topography that will be used as a target. 

        call yelmo_update_equil(yelmo1,time_init,time_tot=20.0,topo_fixed=.FALSE.,dt=1.0,ssa_vel_max=0.0)

        ! Define present topo as present-day dataset for comparison 
        yelmo1%dta%pd%H_ice = yelmo1%tpo%now%H_ice 
        yelmo1%dta%pd%z_srf = yelmo1%tpo%now%z_srf 

        ! ============================================================================================
        ! Step 2: Run the model for several ka in hybrid mode with topo_fixed to
        ! spin up the thermodynamics and have a reference state to reset.
        ! Store the reference state for future use.
        
        call yelmo_update_equil(yelmo1,time_init,time_tot=20e3,topo_fixed=.TRUE.,dt=5.0,ssa_vel_max=5e3)

        ! Write a restart file 
        call yelmo_restart_write(yelmo1,file_restart_init,time_init)
!         stop "**** Done ****"

    end if 

    ! Store the reference state
    yelmo_ref    = yelmo1
    
    ! Store initial optimization parameter choices 
    tau     = rel_tau1 
    H_scale = scale_H1 

    ! Initialize the 2D output file and write the initial model state 
    call yelmo_write_init(yelmo1,file2D,time_init,units="years")  
    call write_step_2D_opt(yelmo1,file2D,time_init,cf_ref_dot,mask_noice,tau,H_scale)  
    
    write(*,*) "Starting optimization..."

if (opt_method .eq. 1) then 
    ! Error method (Pollard and De Conto, 2012)

    ! Initialize timing variables 
    time = time_init 
    n_now = 1 

    do q = 1, qmax 

        ! === Optimization parameters =========
        
        tau     = get_opt_param(time,time1=rel_time1,time2=rel_time2,p1=rel_tau1,p2=rel_tau2,q=rel_q)
        H_scale = get_opt_param(time,time1=scale_time1,time2=scale_time2,p1=scale_H1,p2=scale_H2,q=1.0)
        
        ! Set model tau, and set yelmo relaxation switch (1: shelves relaxing; 0: no relaxation)
        yelmo1%tpo%par%topo_rel_tau = tau 
        yelmo1%tpo%par%topo_rel = 1 
        if (time .gt. rel_time2) yelmo1%tpo%par%topo_rel = 0 

        ! === Update cf_ref and reset model ===================

        if (q .gt. 1) then
            ! Perform optimization after first iteration

            ! Update cf_ref based on error metric(s) 
            call update_cf_ref_thickness_simple(yelmo1%dyn%now%cf_ref,cf_ref_dot,yelmo1%tpo%now%H_ice, &
                            yelmo1%bnd%z_bed,yelmo1%dyn%now%ux_bar,yelmo1%dyn%now%uy_bar, &
                            yelmo1%dta%pd%H_ice,yelmo1%dta%pd%H_grnd.le.0.0_prec,yelmo1%tpo%par%dx, &
                            cf_min,cf_max,H_scale)

        end if 
        
        ! Reset model to the initial state, but with updated cf_ref field and model time
        yelmo_ref%dyn%now%cf_ref = yelmo1%dyn%now%cf_ref 
        yelmo1 = yelmo_ref  
        call yelmo_set_time(yelmo1,time) 

        ! === Update time_iter ==================
        time_end = time_iter
        if (q .eq. qmax) time_end = time_steady

        write(*,"(a,i4,f10.1,i4,f10.1,f12.1,f10.1)") "iter_par: ", q, time, &
                            yelmo1%tpo%par%topo_rel, tau, H_scale, time_end

        ! Perform iteration loop to diagnose error for modifying C_bed 
        do n = 1, int(time_end/dtt)
        
            time = time + dtt 

            ! Update ice sheet 
            call yelmo_update(yelmo1,time)

            if (mod(nint(time*100),nint(dt2D_out*100))==0) then
                call write_step_2D_opt(yelmo1,file2D,time,cf_ref_dot,mask_noice,tau,H_scale)
            end if 

        end do 

    end do 

else 
    ! Ratio method (Le clec’h et al, 2019) - needs revising

    ! Initialize timing variables 
    time = time_init 
    
    do q = 1, qmax 

        ! Reset model to the initial state (including H_w) and time, with updated C_bed field 
        yelmo_ref%dyn%now%C_bed = yelmo1%dyn%now%C_bed 
        yelmo1 = yelmo_ref 
        time   = 0.0 
        call yelmo_set_time(yelmo1,time) 
        
        ! Perform C_bed tuning step 
        do n = 1, int(time_tune)
        
            time = time + 1.0

            ! Update ice sheet 
            call yelmo_update(yelmo1,time)

            ! Update C_bed based on error metric(s) 
            call update_cf_ref_thickness_ratio(yelmo1%dyn%now%cf_ref,cf_ref_dot,yelmo1%tpo%now%H_ice, &
                            yelmo1%bnd%z_bed,yelmo1%dyn%now%ux_bar,yelmo1%dyn%now%uy_bar, &
                            yelmo1%dyn%now%uxy_i_bar,yelmo1%dyn%now%uxy_b,yelmo1%dta%pd%H_ice, &
                            yelmo1%tpo%par%dx,cf_min,cf_max=cf_max)

        end do 

        ! Perform iteration loop to diagnose error for modifying C_bed 
        do n = 1, int(time_iter)
        
            time = time + 1.0

            ! Update ice sheet 
            call yelmo_update(yelmo1,time)

        end do 

        ! Write the current solution 
        call write_step_2D_opt(yelmo1,file2D,time,cf_ref_dot,mask_noice,tau,H_scale)
        
    end do 

end if 
    
    ! Write a final restart file 
    call yelmo_restart_write(yelmo1,file_restart,time)

    ! Finalize program
    call yelmo_end(yelmo1,time=time)

    ! Stop timing 
    call cpu_time(cpu_end_time)

    write(*,"(a,f12.3,a)") "Time  = ",(cpu_end_time-cpu_start_time)/60.0 ," min"
    write(*,"(a,f12.1,a)") "Speed = ",(1e-3*(time-time_init))/((cpu_end_time-cpu_start_time)/3600.0), " kiloyears / hr"

contains

    subroutine write_step_2D_opt(ylmo,filename,time,cf_ref_dot,mask_noice,tau,H_scale)

        implicit none 
        
        type(yelmo_class), intent(IN) :: ylmo
        character(len=*),  intent(IN) :: filename
        real(prec), intent(IN) :: time
        real(prec), intent(IN) :: cf_ref_dot(:,:)
        logical,    intent(IN) :: mask_noice(:,:) 
        real(prec), intent(IN) :: tau 
        real(prec), intent(IN) :: H_scale 

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
        call nc_write(filename,"opt_H_scale",H_scale,units="m",long_name="Error scaling constant", &
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
        
        call nc_write(filename,"cf_ref",yelmo1%dyn%now%cf_ref,units="",long_name="Bed friction scalar", &
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
        
        tmp = ylmo%dta%pd%err_z_srf
        call filter_gaussian(var=tmp,sigma=2.0*ylmo%tpo%par%dx,dx=ylmo%tpo%par%dx, &
                                mask=ylmo%dta%pd%err_z_srf .ne. 0.0)
        
        call nc_write(filename,"z_srf_sm_pd_err",tmp,units="m",long_name="Smooth surface elevation error wrt present day", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine write_step_2D_opt

    subroutine guess_cf_ref(cf_ref,tau_d,uxy_obs,H_obs,H_grnd,u0,cf_min,cf_max)
        ! Use suggestion by Morlighem et al. (2013) to guess friction
        ! assuming tau_b ~ tau_d, and u_b = u_obs:
        !
        ! For a linear law, tau_b = beta * u_b, so 
        ! beta = tau_b / u_b = tau_d / (u_obs+ebs), ebs=0.1 to avoid divide by zero 
        ! beta = cf_ref/u0 * N_eff, so:
        ! cf_ref = (tau_d/(u_obs+ebs)) * (u0/N_eff)

        implicit none 

        real(prec), intent(OUT) :: cf_ref(:,:) 
        real(prec), intent(IN)  :: tau_d(:,:) 
        real(prec), intent(IN)  :: uxy_obs(:,:) 
        real(prec), intent(IN)  :: H_obs(:,:)
        real(prec), intent(IN)  :: H_grnd(:,:)
        real(prec), intent(IN)  :: u0 
        real(prec), intent(IN)  :: cf_min 
        real(prec), intent(IN)  :: cf_max  

        ! Local variables 
        real(prec), parameter :: ebs = 0.1          ! [m/yr] To avoid divide by zero 

        where (H_obs .eq. 0.0_prec .or. H_grnd .eq. 0.0_prec) 
            ! Set floating or ice-free points to minimum 
            cf_ref = cf_min 

        elsewhere 
            ! Apply equation 

            ! Linear law: 
            cf_ref = (tau_d / (uxy_obs + ebs)) * (u0 / (rho_ice*g*H_obs + 1.0_prec))

        end where 

        where (cf_ref .gt. cf_max) cf_ref = cf_max 

        return 

    end subroutine guess_cf_ref 

    subroutine update_cf_ref_thickness_simple(cf_ref,cf_ref_dot,H_ice,z_bed,ux,uy,H_obs,is_float_obs,dx,cf_min,cf_max,H_scale)

        implicit none 

        real(prec), intent(INOUT) :: cf_ref(:,:) 
        real(prec), intent(INOUT) :: cf_ref_dot(:,:) 
        real(prec), intent(IN)    :: H_ice(:,:) 
        real(prec), intent(IN)    :: z_bed(:,:) 
        real(prec), intent(IN)    :: ux(:,:) 
        real(prec), intent(IN)    :: uy(:,:) 
        real(prec), intent(IN)    :: H_obs(:,:) 
        logical,    intent(IN)    :: is_float_obs(:,:) 
        real(prec), intent(IN)    :: dx 
        real(prec), intent(IN)    :: cf_min 
        real(prec), intent(IN)    :: cf_max
        real(prec), intent(IN)    :: H_scale            ! [m] H_scale = 1000.0 m by default

        ! Local variables 
        integer :: i, j, nx, ny, i1, j1  
        real(prec) :: dx_km, f_dz, f_dz_lim, f_scale   
        real(prec) :: ux_aa, uy_aa
        real(prec) :: H_ice_now, H_obs_now 

        real(prec) :: angle 

        real(prec), allocatable   :: cf_prev(:,:) 
        real(prec) :: wts0(5,5), wts(5,5) 

        real(prec), parameter :: ulim_divide = 5.0      ! [m/a] Limit to consider we are near ice divide 

        nx = size(cf_ref,1)
        ny = size(cf_ref,2) 

        dx_km = dx*1e-3  
        
        allocate(cf_prev(nx,ny))

        ! Optimization parameters 
        !H_scale  = 1000.0           ! [m]   **Now an input parameter to change with time 
        f_dz_lim = 1.5              ! [--] 

        ! Get Gaussian weights 
        wts0 = gauss_values(dx_km,dx_km,sigma=dx_km*1.5,n=5)

        ! Store initial cf_ref solution 
        cf_prev = cf_ref 

        do j = 3, ny-2 
        do i = 3, nx-2 

            if ( H_ice(i,j) .ne. 0.0 .or. H_obs(i,j) .ne. 0.0) then 
                ! Update coefficient where ice or ice_obs exists

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

                ! Apply to current and downstream node(s) =========

!                 ux_aa = 0.5*(ux(i,j)+ux(i+1,j))
!                 uy_aa = 0.5*(uy(i,j)+uy(i,j+1))
                
!                 angle = atan2(uy_aa,ux_aa)


!                 ! First apply locally. This will be overwritten
!                 ! if point is downstream of another node.  
!                 if (i1 .ne. i .or. j1 .ne. j) then 
!                     cf_ref(i,j) = cf_prev(i,j)* ( 0.2*(f_scale-1.0)+1.0 )
!                 end if 

!                 if ( abs(ux_aa) .gt. abs(uy_aa) ) then 
!                     ! Downstream in x-direction 
!                     j1 = j 
!                     if (abs(ux_aa) .lt. ulim_divide) then 
!                         ! Near ice-divide 
!                         i1 = i 
!                     else if (ux_aa .lt. 0.0) then 
!                         i1 = i-1 
!                     else
!                         i1 = i+1
!                     end if 

!                 else 
!                     ! Downstream in y-direction 
!                     i1 = i
!                     if (abs(uy_aa) .lt. ulim_divide) then 
!                         ! Near ice-divide 
!                         j1 = j  
!                     else if (uy_aa .lt. 0.0) then 
!                         j1 = j-1
!                     else
!                         j1 = j+1
!                     end if 

!                 end if 

!                 ! Apply coefficent scaling at correct node
!                 cf_ref(i1,j1) = cf_prev(i1,j1)*f_scale

                
                ! Apply to current node 
                cf_ref(i,j) = cf_prev(i,j)*f_scale 

            end if 

        end do 
        end do 

        ! Ensure cf_ref is not below lower or upper limit 
        where (cf_ref .lt. cf_min) cf_ref = cf_min 
        where (cf_ref .gt. cf_max) cf_ref = cf_max 

        ! Additionally, apply a Gaussian filter to cf_ref to ensure smooth transitions
        !call filter_gaussian(var=cf_ref,sigma=dx_km*0.2,dx=dx_km)     !,mask=err_z_srf .ne. 0.0)
        
        ! Ensure where obs are floating, set cf_ref = cf_min 
        where(is_float_obs) cf_ref =cf_min 

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

    function get_opt_param(time,time1,time2,p1,p2,q) result(p) 
        ! Determine value of parameter as a function of time 

        implicit none 

        real(prec), intent(IN) :: time 
        real(prec), intent(IN) :: time1 
        real(prec), intent(IN) :: time2
        real(prec), intent(IN) :: p1
        real(prec), intent(IN) :: p2
        real(prec), intent(IN) :: q         ! Non-linear exponent (q=1.0 or higher)
        real(prec) :: p 

        if (time .le. time1) then 
            p = p1 
        else if (time .ge. time2) then 
            p = p2 
        else 
            ! Linear interpolation 
            p = p1 + (p2-p1)* ((time-time1)/(time2-time1))**q 
        end if  

        return 

    end function get_opt_param 

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



