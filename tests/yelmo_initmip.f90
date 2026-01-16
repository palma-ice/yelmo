

program yelmo_test

    use omp_lib

    use nml 
    use ncio 
    use timestepping
    use timeout

    use yelmo
    use basal_dragging 
    
    use ice_optimization 
    
    implicit none 

    type(tstep_class)   :: ts
    type(timeout_class) :: t1D, t2D, t2Dsm

    type(yelmo_class)   :: yelmo1 
    
    character(len=512) :: path_par
    character(len=256) :: file_restart 
    
    logical,  allocatable :: mask_noice(:,:)  

    type ctrl_params
        character(len=56) :: run_step
        character(len=56) :: tstep_method
        real(wp) :: tstep_const
        real(wp) :: time_init
        real(wp) :: time_end
        real(wp) :: time_equil      ! Only for spinup
        real(wp) :: dtt

        logical  :: with_ice_sheet 
        character(len=56) :: equil_method
        
        character(len=512) :: clim_nm
          
        logical :: load_cb_ref 
        character(len=256) :: file_cb_ref 

        logical :: load_bmelt
        character(len=256) :: file_bmelt 

        real(wp) :: bmb_shlf_const
        real(wp) :: dT_ann
        real(wp) :: z_sl 
    end type 

    type(ctrl_params)    :: ctl
    type(ice_opt_params) :: opt

    real(wp) :: dtt_now

    real(8) :: cpu_start_time, cpu_end_time, cpu_dtime

    ! Code for testing restarts
    logical, parameter  :: test_restart = .FALSE. 
    real(wp), parameter :: time_r       = 50.0_wp 
    type(yelmo_class)   :: yelmo_r
    character(len=256)  :: file1D_r, file2D_r, file_restart_r
    
    ! Start timing 
    call yelmo_cpu_time(cpu_start_time)
    
    ! Determine the parameter file from the command line 
    call yelmo_load_command_line_args(path_par)

    ! Timing and other parameters 
    call nml_read(path_par,"ctrl","tstep_method",   ctl%tstep_method)       ! Calendar choice ("const" or "rel")
    call nml_read(path_par,"ctrl","tstep_const",    ctl%tstep_const)        ! Assumed time bp for const method
    call nml_read(path_par,"ctrl","time_init",      ctl%time_init)          ! [yr] Starting time
    call nml_read(path_par,"ctrl","time_end",       ctl%time_end)           ! [yr] Ending time
    call nml_read(path_par,"ctrl","time_equil",     ctl%time_equil)         ! [yr] Years to equilibrate first
    call nml_read(path_par,"ctrl","dtt",            ctl%dtt)                ! [yr] Main loop time step 
    call nml_read(path_par,"ctrl","with_ice_sheet", ctl%with_ice_sheet)     ! Include an active ice sheet 
    call nml_read(path_par,"ctrl","equil_method",   ctl%equil_method)       ! What method should be used for spin-up?
    call nml_read(path_par,"ctrl","clim_nm",        ctl%clim_nm)            ! Namelist group holding climate information
    
    call nml_read(path_par,"ctrl","load_cb_ref",    ctl%load_cb_ref)        ! Load cb_ref from file? Otherwise define from till_cf_ref + inline tuning
    call nml_read(path_par,"ctrl","file_cb_ref",    ctl%file_cb_ref)        ! Filename holding cb_ref to load 

    call nml_read(path_par,"ctrl","load_bmelt",     ctl%load_bmelt)         ! Load bmelt from file?
    call nml_read(path_par,"ctrl","file_bmelt",     ctl%file_bmelt)         ! Filename holding bmelt field to load 
        
    ! Load climate (eg, clim_pd or clim_lgm)
    call nml_read(path_par,ctl%clim_nm,  "bmb_shlf_const",  ctl%bmb_shlf_const)            ! [yr] Constant imposed bmb_shlf value
    call nml_read(path_par,ctl%clim_nm,  "dT_ann",          ctl%dT_ann)                    ! [K] Temperature anomaly (atm)
    call nml_read(path_par,ctl%clim_nm,  "z_sl",            ctl%z_sl)                      ! [m] Sea level relative to present-day

    ! Get output times
    call timeout_init(t1D,  path_par,"t1D",  "small",  ctl%time_init,ctl%time_end)
    call timeout_init(t2Dsm,path_par,"t2Dsm","medium", ctl%time_init,ctl%time_end)
    call timeout_init(t2D,  path_par,"t2D",  "heavy",  ctl%time_init,ctl%time_end)
    
    if (trim(ctl%equil_method) .eq. "opt") then 
        ! Load optimization parameters 

        ! Initially set to zero 
        opt%tf_basins = 0 

        call optimize_par_load(opt,path_par,"opt")

    end if 

    ! Define input and output locations
    t2Dsm%filename = "yelmo2Dsm.nc"
    t2D%filename   = "yelmo2D.nc"
    file_restart   = "yelmo_restart.nc"

    ! === Initialize timestepping ===
    
    call tstep_init(ts,ctl%time_init,ctl%time_end,method=ctl%tstep_method,units="year", &
                                            time_ref=1950.0_wp,const_rel=ctl%tstep_const)

    write(*,*)
    write(*,*) "timestepping:   ",  trim(ts%method)
    if (trim(ts%method) .eq. "const") then 
        write(*,*) "time_equil: ",    ctl%time_equil
    end if 

    write(*,*) "time    = ", ts%time 
    write(*,*) "time_bp = ", ts%time_rel 
    write(*,*) 
    
    ! === Initialize ice sheet model =====

    ! Initialize data objects and load initial topography
    call yelmo_init(yelmo1,filename=path_par,grid_def="file",time=ts%time)

    ! Ensure optimization fields are allocated
    allocate(opt%cf_min(yelmo1%grd%nx,yelmo1%grd%ny))
    allocate(opt%cf_max(yelmo1%grd%nx,yelmo1%grd%ny))
    opt%cf_min = yelmo1%dyn%par%till_cf_min
    opt%cf_max = yelmo1%dyn%par%till_cf_ref

    ! === Set initial boundary conditions for current time and yelmo state =====
    ! ybound: z_bed, z_sl, H_sed, smb, T_srf, bmb_shlf , Q_geo

    yelmo1%bnd%z_sl     = ctl%z_sl          ! [m]
    yelmo1%bnd%H_sed    = 0.0               ! [m]
    yelmo1%bnd%Q_geo    = 50.0              ! [mW/m2]
    
    ! Impose present-day surface mass balance and present-day temperature field plus any anomaly
    yelmo1%bnd%smb      = yelmo1%dta%pd%smb                 ! [m.i.e./a]
    yelmo1%bnd%T_srf    = yelmo1%dta%pd%T_srf + ctl%dT_ann  ! [K]
    
    if (ctl%load_bmelt) then

        ! Parse filenames with grid information
        call yelmo_parse_path(ctl%file_bmelt,yelmo1%par%domain,yelmo1%par%grid_name)
 
        call nc_read(ctl%file_bmelt,"bm_ac_reese",yelmo1%bnd%bmb_shlf)
        yelmo1%bnd%bmb_shlf = -yelmo1%bnd%bmb_shlf      ! Negative because bmb = -bmelt 
    else 
        yelmo1%bnd%bmb_shlf = ctl%bmb_shlf_const        ! [m.i.e./a]
    end if 

    yelmo1%bnd%T_shlf   = yelmo1%bnd%c%T0                ! [K]   

    if (ctl%dT_ann .lt. 0.0) yelmo1%bnd%T_shlf   = yelmo1%bnd%c%T0 + ctl%dT_ann*0.25_wp      ! [K] Oceanic temp anomaly
    
    call yelmo_print_bound(yelmo1%bnd)

    ! Define no-ice mask from present-day data
    allocate(mask_noice(yelmo1%grd%nx,yelmo1%grd%ny))
    mask_noice = .FALSE. 
    ! Present-day
    if (ctl%dT_ann .ge. 0.0) then
        where(yelmo1%dta%pd%H_ice .le. 0.0) mask_noice = .TRUE. 
    end if 

    ! Special treatment for Greenland
    if (trim(yelmo1%par%domain) .eq. "Greenland") then 
        
        ! Present-day
        if (ctl%dT_ann .ge. 0.0) then 
            ! Impose additional negative mass balance to no ice points of 2 [m.i.e./a] melting
            where(mask_noice) yelmo1%bnd%smb = yelmo1%dta%pd%smb - 2.0 
        end if 

    end if 

    ! Special treatment for Antarctica
    if (trim(yelmo1%par%domain) .eq. "Antarctica") then 
        
        ! Present-day
        if (ctl%dT_ann .ge. 0.0) then 
            where(mask_noice) yelmo1%bnd%bmb_shlf = -2.0                ! [m/a]        
        end if 

        ! LGM
        if (ctl%dT_ann .lt. 0.0) then 
            where(yelmo1%bnd%smb .le. 0.1) yelmo1%bnd%smb = 0.1         ! [m/a]
        end if 

        ! Present-day and LGM 
        where(yelmo1%bnd%regions .eq. 2.0) yelmo1%bnd%smb = -1.0        ! [m/a]

    end if 

    ! Initialize state variables (dyn,therm,mat)
    ! (initialize temps with robin method with a cold base)
    call yelmo_init_state(yelmo1,time=ts%time,thrm_method="robin-cold")

    ! ===== basal friction optimization ======
    if (trim(ctl%equil_method) .eq. "opt") then 
        
        ! Ensure that cb_ref will be optimized (till_method == set externally) 
        yelmo1%dyn%par%till_method = -1  

        ! If not using restart, prescribe cb_ref to initial guess 
        if (.not. yelmo1%par%use_restart) then
            yelmo1%dyn%now%cb_ref = opt%cf_init 
        end if 

    else 
        ! ============================================================
        ! Load or define cb_ref 

        if (ctl%load_cb_ref) then 

            ! Parse filename with grid information
            call yelmo_parse_path(ctl%file_cb_ref,yelmo1%par%domain,yelmo1%par%grid_name)

            ! Load cb_ref from specified file 
            call nc_read(ctl%file_cb_ref,"cb_ref",yelmo1%dyn%now%cb_ref)

            ! Make sure that present-day shelves have minimal cb_ref values 
            where(yelmo1%tpo%now%f_grnd .eq. 0.0) yelmo1%dyn%now%cb_ref = 0.05 

        else
            ! Define cb_ref inline 

            yelmo1%dyn%now%cb_ref = yelmo1%dyn%par%till_cf_ref 

        end if 

    end if 
    ! ============================================================
    
    if (.not. yelmo1%par%use_restart) then

        ! Run full model (tpo,dyn,thrm) without advection to clean initial topo via mass balance
        call yelmo_update_equil(yelmo1,ts%time,time_tot=10.0_wp,dt=1.0_wp,topo_fixed=.FALSE.,tpo_solver="none")
        
        ! Run full model with correct solver (tpo,dyn,thrm)
        call yelmo_update_equil(yelmo1,ts%time,time_tot=1.0_wp,dt=0.2_wp,topo_fixed=.FALSE.)

        ! Next equilibrate thermodynamics further and maintain constant ice topopgraphy (for speed)
        !call yelmo_update_equil(yelmo1,ts%time,time_tot=1e3,dt=10.0_wp,topo_fixed=.TRUE.)

    end if 


    if (yelmo1%par%use_restart .and. yelmo1%par%restart_interpolated .eq. 1) then 

        write(*,*) "Running restart topo smoothing step..."

        ! Run model with no advection
        call yelmo_update_equil(yelmo1,ts%time,time_tot=10.0_wp,dt=1.0_wp, &
                              tpo_solver="none",topo_fixed=.FALSE.,dyn_solver="ssa")

        ! Run thermodynamics with SSA solver very briefly to smooth it out
        call yelmo_update_equil(yelmo1,ts%time,time_tot=10.0_wp,dt=1.0_wp, &
                                                topo_fixed=.TRUE.,dyn_solver="ssa")

        ! Run full model (tpo,dyn,thrm) with SSA solver very briefly to smooth it out
        call yelmo_update_equil(yelmo1,ts%time,time_tot=1.0_wp,dt=0.2_wp, &
                                                topo_fixed=.FALSE.,dyn_solver="ssa")

        if (yelmo1%grd%dx .le. 8e3_wp) then 
            ! Perform additional ssa smoothing step for higher resolution simulations

            call yelmo_update_equil(yelmo1,ts%time,time_tot=100.0_wp,dt=1.0_wp, &
                                                topo_fixed=.FALSE.,dyn_solver="ssa")

        end if 

    
    end if 

    ! Initialize output files

    call yelmo_regions_write(yelmo1,ts%time,init=.TRUE.,units="years")

    if (t2Dsm%active) then
        call yelmo_write_init(yelmo1,t2Dsm%filename,time_init=ts%time,units="years")
    end if
    if (t2D%active) then
        call yelmo_write_init(yelmo1,t2D%filename,time_init=ts%time,units="years")
    end if
    
    ! if (with_anom) then 
    !     ! Warm up the ice sheet to impose some changes 
    !     yelmo1%bnd%T_srf    = yelmo1%dta%pd%T_srf + 5.0
    !     yelmo1%bnd%bmb_shlf = -10.0 
    !     where (yelmo1%dta%pd%smb .lt. 0.0) yelmo1%bnd%smb = yelmo1%dta%pd%smb - 1.0 
    ! end if 

    ! ==== Begin main time loop =====

    dtt_now = ctl%dtt
    call tstep_print_header(ts)

    do while (.not. ts%is_finished)

        ! == Update timestep ===

        call tstep_update(ts,dtt_now)
        call tstep_print(ts)
        
!         ! Update temperature and smb as needed in time (ISMIP6)
!         if (ts%time .ge. -10e6 .and. ts%time .lt. -10e3) then 
!             ! Glacial period, impose cold climate 
!             yelmo1%bnd%T_srf = yelmo1%dta%pd%T_srf - 10.0 

!         else if (ts%time .ge. -10e3 .and. ts%time .lt. -6e3) then
!             ! Holocene optimum 
!             yelmo1%bnd%T_srf = yelmo1%dta%pd%T_srf + 1.0 

!         else  ! time .ge. -6e3
!             ! Entering Holocene, impose present-day temperatures 
!             yelmo1%bnd%T_srf = yelmo1%dta%pd%T_srf
!         end if 
        
        ! Spin-up procedure - only relevant for time_elapsed <= time_equil
        select case(trim(ctl%equil_method))
            
            case("opt")
                ! ===== basal friction optimization ==================

                if (ts%time_elapsed .le. opt%rel_time2) then 
                    ! Apply relaxation to the model 

                    ! Update model relaxation time scale and error scaling (in [m])
                    call optimize_set_transient_param(opt%rel_tau,ts%time_elapsed,time1=opt%rel_time1,time2=opt%rel_time2, &
                                                    p1=opt%rel_tau1,p2=opt%rel_tau2,m=opt%rel_m)
                    
                    ! Set model tau, and set yelmo relaxation switch (4: gl line and grounding zone relaxing; 0: no relaxation)
                    yelmo1%tpo%par%topo_rel_tau = opt%rel_tau 
                    yelmo1%tpo%par%topo_rel     = 4
                
                else 
                    ! Turn-off relaxation now

                    yelmo1%tpo%par%topo_rel = 0 

                end if 

                ! === Optimization update step =========

                if (opt%opt_cf .and. &
                        (ts%time_elapsed .ge. opt%cf_time_init .and. ts%time_elapsed .le. opt%cf_time_end) ) then
                    ! Perform cf_ref optimization
                
                    ! Update cb_ref based on error metric(s) 
                    call optimize_cb_ref(yelmo1%dyn%now%cb_ref,yelmo1%tpo%now%H_ice, &
                                                    yelmo1%tpo%now%dHidt,yelmo1%bnd%z_bed,yelmo1%bnd%z_sl,yelmo1%dyn%now%ux_s,yelmo1%dyn%now%uy_s, &
                                                    yelmo1%dta%pd%H_ice,yelmo1%dta%pd%uxy_s,yelmo1%dta%pd%H_grnd, &
                                                    opt%cf_min,opt%cf_max,yelmo1%tpo%par%dx,opt%sigma_err,opt%sigma_vel,opt%tau_c,opt%H0, &
                                                    dt=ctl%dtt,fill_method=opt%fill_method,fill_dist=opt%sigma_err,cb_tgt=yelmo1%dyn%now%cb_tgt)
                    
                end if

                if (opt%opt_tf .and. &
                        (ts%time_elapsed .ge. opt%tf_time_init .and. ts%time_elapsed .le. opt%tf_time_end) ) then
                    ! Perform tf_corr optimization

                    write(io_unit_err,*) "yelmo_initmip:: Error: thermal forcing optimization not yet defined."
                    write(io_unit_err,*) "Best solution for now: set opt_tf=False."
                    stop

                end if 

            case("relax")
                ! ===== relaxation spinup ==================

                ! Turn on relaxation for now, to let thermodynamics equilibrate
                ! without changing the topography too much. Important when 
                ! effective pressure = f(thermodynamics).

                yelmo1%tpo%par%topo_rel     = 2
                yelmo1%tpo%par%topo_rel_tau = 50.0 
                write(*,*) "timelog, tau = ", yelmo1%tpo%par%topo_rel_tau
                
                if ( ts%time_elapsed .ge. ctl%time_equil) then
                    ! Finally, ensure all relaxation is disabled and continue as normal.

                        yelmo1%tpo%par%topo_rel     = 0
                        write(*,*) "timelog, relaxation off..."
                    
                end if 
        end select 

        ! == UPDATE YELMO =======================================================


        if (ctl%with_ice_sheet) call yelmo_update(yelmo1,ts%time)


        ! == MODEL OUTPUT =======================================================

        ! if (mod(nint(ts%time*100),nint(ctl%dt1D_out*100))==0) then 
        !     call yelmo_write_reg_step(yelmo1,t1D%filename,time=ts%time) 
        ! end if 

        if (timeout_check(t1D,ts%time)) then 
            call yelmo_regions_write(yelmo1,ts%time)
        end if 

        if (timeout_check(t2Dsm,ts%time)) then 
            call yelmo_write_step(yelmo1,t2Dsm%filename,ts%time,compare_pd=.TRUE.)
        end if

        if (timeout_check(t2D,ts%time)) then
            call write_step_2D(yelmo1,t2D%filename,time=ts%time)
        end if
        
        if (mod(ts%time,10.0)==0 .and. (.not. yelmo_log)) then
            write(*,"(a,f14.4)") "yelmo:: time = ", ts%time
        end if 
        
    end do 
    ! == Finished time loop == 

    ! Write a final restart file 
    if (.not. test_restart) then 
        call yelmo_restart_write(yelmo1,file_restart,ts%time)
    end if 

    ! Finalize program
    call yelmo_end(yelmo1,time=ts%time)
    
    ! Stop timing 
    call yelmo_cpu_time(cpu_end_time,cpu_start_time,cpu_dtime)

    write(*,"(a,f12.3,a)") "Time  = ",cpu_dtime/60.0 ," min"
    write(*,"(a,f12.1,a)") "Speed = ",(1e-3*(ctl%time_end-ctl%time_init))/(cpu_dtime/3600.0), " kiloyears / hr"

contains

    subroutine test_restart_step()

        implicit none

        ! ajr: This was in the time loop. The idea is to write a restart file, then
        ! load it and run yelmo, and see if the results match. Storing it here, to have
        ! code out of the way. But ideally we can write this inside of routine like this
        ! to avoid polluting the time loop. 

        ! if (test_restart) then 

        !     if (ts%time .eq. time_r+dtt_now) then 

        !         file1D_r       = "yelmo1D_r.nc"
        !         file2D_r       = "yelmo2D_r.nc"
        !         file_restart_r = "yelmo_restart_r.nc"

        !         ctl%dt2D_out = dtt_now
        !         ctl%dt1D_out = dtt_now

        !         ! Initialize data objects and load initial topography
        !         call yelmo_init(yelmo_r,filename=path_par,grid_def="file",time=time_r)

        !         yelmo_r%par%restart     = "yelmo_restart.nc"
        !         yelmo_r%par%use_restart = .TRUE. 
        !         yelmo_r%bnd = yelmo1%bnd

        !         ! Initialize state variables (dyn,therm,mat)
        !         ! (initialize temps with robin method with a cold base)
        !         call yelmo_init_topo(yelmo_r,"ytopo",path_par,time_r)
        !         call yelmo_init_state(yelmo_r,time=time_r,thrm_method="robin-cold")

        !         ! Write restart file to compare with expected 
        !         call yelmo_restart_write(yelmo_r,file_restart_r,time_r)

        !         ! 2D file 
        !         call yelmo_write_init(yelmo_r,file2D_r,time_init=time_r,units="years")  
                
        !         ! 1D file 
        !         call yelmo_write_reg_init(yelmo_r,file1D_r,time_init=time_r,units="years",mask=yelmo_r%bnd%ice_allowed)
                
        !         call write_step_2D(yelmo_r,file2D_r,time=time_r)
        !         call yelmo_write_reg_step(yelmo_r,file1D_r,time=time_r)  
        !     end if 

        !     if (ts%time .gt. time_r) then
        !         ! Update restarted model
        !         call yelmo_update(yelmo_r,ts%time)
        !     end if 

        ! end if 

        ! -------
        ! Note: code below comes after updating yelmo...

        ! if (test_restart) then 

        !     if (ts%time .eq. time_r) then
        !         ! Write restart file to load from
        !         call yelmo_restart_write(yelmo1,file_restart,ts%time)
        !     end if 

        !     if (ts%time .gt. time_r) then
        !         ! Write restarted output files 
        !         if (mod(nint(ts%time*100),nint(ctl%dt2D_out*100))==0) then
        !             call write_step_2D(yelmo_r,file2D_r,time=ts%time)
        !         end if 

        !         if (mod(nint(ts%time*100),nint(ctl%dt1D_out*100))==0) then 
        !             call yelmo_write_reg_step(yelmo_r,file1D_r,time=ts%time)  
        !         end if
        !     end if 

        ! end if 

        return

    end subroutine test_restart_step

    subroutine write_step_2D(ylmo,filename,time)

        implicit none 
        
        type(yelmo_class), intent(IN) :: ylmo
        character(len=*),  intent(IN) :: filename
        real(wp), intent(IN) :: time

        ! Local variables
        integer  :: ncid, n
        real(wp) :: time_prev 

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_time_index(filename,"time",time,ncid)

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        ! Write model metrics (model speed, dt, eta)
        call yelmo_write_step_model_metrics(filename,ylmo,n,ncid)

        ! Write present-day data metrics (rmse[H],etc)
        call yelmo_write_step_pd_metrics(filename,ylmo,n,ncid)
        
        ! == ISMIP6 specific variables 
        ! http://www.climate-cryosphere.org/wiki/index.php?title=InitMIP-Greenland#Appendix_2_.E2.80.93_Naming_conventions.2C_upload_and_model_output_data.
        
        ! == yelmo_topography ==
        call yelmo_write_var(filename,"H_ice",ylmo,n,ncid)
        call yelmo_write_var(filename,"z_srf",ylmo,n,ncid)
        call yelmo_write_var(filename,"mask_bed",ylmo,n,ncid)
        call yelmo_write_var(filename,"mb_net",ylmo,n,ncid)
        call yelmo_write_var(filename,"mb_resid",ylmo,n,ncid)
        call yelmo_write_var(filename,"mb_err",ylmo,n,ncid)
        call yelmo_write_var(filename,"smb",ylmo,n,ncid)
        call yelmo_write_var(filename,"bmb",ylmo,n,ncid)
        call yelmo_write_var(filename,"cmb",ylmo,n,ncid)
        call yelmo_write_var(filename,"fmb",ylmo,n,ncid)
        call yelmo_write_var(filename,"dmb",ylmo,n,ncid)
        call yelmo_write_var(filename,"H_grnd",ylmo,n,ncid)
        call yelmo_write_var(filename,"N_eff",ylmo,n,ncid)
        call yelmo_write_var(filename,"f_grnd",ylmo,n,ncid)
        call yelmo_write_var(filename,"f_ice",ylmo,n,ncid)
        call yelmo_write_var(filename,"dHidt",ylmo,n,ncid)
        call yelmo_write_var(filename,"dzbdx",ylmo,n,ncid)
        call yelmo_write_var(filename,"dzbdy",ylmo,n,ncid)

        ! == yelmo_dynamics ==
        call yelmo_write_var(filename,"cb_ref",ylmo,n,ncid)
        call yelmo_write_var(filename,"N_eff",ylmo,n,ncid)
        call yelmo_write_var(filename,"c_bed",ylmo,n,ncid)
        call yelmo_write_var(filename,"beta",ylmo,n,ncid)
        call yelmo_write_var(filename,"visc_eff_int",ylmo,n,ncid)
        call yelmo_write_var(filename,"taud",ylmo,n,ncid)
        call yelmo_write_var(filename,"taub",ylmo,n,ncid)
        call yelmo_write_var(filename,"uxy_b",ylmo,n,ncid)
        call yelmo_write_var(filename,"uxy_s",ylmo,n,ncid)
        call yelmo_write_var(filename,"uz",ylmo,n,ncid)
        call yelmo_write_var(filename,"uz_star",ylmo,n,ncid)
        
        ! == yelmo_material ==
        call yelmo_write_var(filename,"enh_bar",ylmo,n,ncid)
        !call yelmo_write_var(filename,"ATT",ylmo,n,ncid)
        call yelmo_write_var(filename,"visc_int",ylmo,n,ncid)
        
        ! == yelmo_thermodynamics ==
        call yelmo_write_var(filename,"T_prime",ylmo,n,ncid)
        call yelmo_write_var(filename,"f_pmp",ylmo,n,ncid)
        call yelmo_write_var(filename,"Q_b",ylmo,n,ncid)
        call yelmo_write_var(filename,"Q_ice_b",ylmo,n,ncid)
        call yelmo_write_var(filename,"bmb_grnd",ylmo,n,ncid)
        call yelmo_write_var(filename,"H_w",ylmo,n,ncid)
        call yelmo_write_var(filename,"T_rock",ylmo,n,ncid)
        
        !call yelmo_write_var(filename,"Q_strn",ylmo,n,ncid)
        call nc_write(filename,"Q_strn",ylmo%thrm%now%Q_strn/(ylmo%bnd%c%rho_ice*ylmo%thrm%now%cp),units="K a-1",long_name="Strain heating", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)


        ! == yelmo_boundaries ==
        call yelmo_write_var(filename,"z_bed",ylmo,n,ncid)
        call yelmo_write_var(filename,"z_sl",ylmo,n,ncid)
        call yelmo_write_var(filename,"smb_ref",ylmo,n,ncid)
        call yelmo_write_var(filename,"T_srf",ylmo,n,ncid)
        call yelmo_write_var(filename,"bmb_shlf",ylmo,n,ncid)
        call yelmo_write_var(filename,"Q_geo",ylmo,n,ncid)
        
        ! == yelmo_data (comparison with present-day) ==
        call yelmo_write_var(filename,"pd_err_H_ice",ylmo,n,ncid)
        call yelmo_write_var(filename,"pd_err_z_srf",ylmo,n,ncid)
        call yelmo_write_var(filename,"pd_err_uxy_s",ylmo,n,ncid)
        
        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine write_step_2D

    subroutine modify_cb_ref(dyn,tpo,thrm,bnd,grd,domain,cb_ref)
        ! Modify cb_ref [unitless] with location specific tuning 

        implicit none
        
        type(ydyn_class),   intent(INOUT) :: dyn
        type(ytopo_class),  intent(IN)    :: tpo 
        type(ytherm_class), intent(IN)    :: thrm
        type(ybound_class), intent(IN)    :: bnd  
        type(ygrid_class),  intent(IN)    :: grd
        character(len=*),   intent(IN)    :: domain 
        real(wp),           intent(INOUT) :: cb_ref(:,:) 

        integer  :: i, j, nx, ny 
        integer  :: i1, i2, j1, j2 
        real(wp) :: f_scale 
            
            ! Additionally modify cb_ref
            if (trim(domain) .eq. "Antarctica") then


                ! Increase friction - feeding the Ronne ice shelf from the South
!                 call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*2.0,x0=-800.0, y0= 100.0,sigma=400.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
!                 call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*2.0,x0=-980.0, y0=-400.0,sigma=200.0,xx=grd%x*1e-3,yy=grd%y*1e-3)
                
                ! Reduction friction - feeding the Ross ice shelf from the East
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.05_wp,x0= 130.0_wp, y0=-550.0_wp, sigma=100.0_wp,xx=grd%x*1e-3_wp,yy=grd%y*1e-3_wp)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.05_wp,x0= 280.0_wp, y0=-760.0_wp, sigma=100.0_wp,xx=grd%x*1e-3_wp,yy=grd%y*1e-3_wp)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.05_wp,x0= 380.0_wp, y0=-960.0_wp, sigma=100.0_wp,xx=grd%x*1e-3_wp,yy=grd%y*1e-3_wp)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.05_wp,x0= 400.0_wp, y0=-1150.0_wp,sigma=100.0_wp,xx=grd%x*1e-3_wp,yy=grd%y*1e-3_wp)


if (.FALSE.) then 
                ! Increase - feeding the Ronne ice shelf from the North
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*4.00_wp,x0=-700.0_wp, y0=    0.0_wp,sigma=200.0_wp,xx=grd%x*1e-3_wp,yy=grd%y*1e-3_wp)
                
                ! Increase - Southeast Antarctica inland
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*2.00_wp,x0=1500.0_wp, y0= -550.0_wp,sigma=200.0_wp,xx=grd%x*1e-3_wp,yy=grd%y*1e-3_wp)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*2.00_wp,x0=1700.0_wp, y0=-1000.0_wp,sigma=200.0_wp,xx=grd%x*1e-3_wp,yy=grd%y*1e-3_wp)
                
                ! Reduction - South pole 
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.10_wp,x0=   0.0_wp, y0=   0.0_wp, sigma=400.0_wp,xx=grd%x*1e-3_wp,yy=grd%y*1e-3_wp)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.20_wp,x0=   0.0_wp, y0= 600.0_wp, sigma=300.0_wp,xx=grd%x*1e-3_wp,yy=grd%y*1e-3_wp)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.10_wp,x0= 500.0_wp, y0=-500.0_wp, sigma=400.0_wp,xx=grd%x*1e-3_wp,yy=grd%y*1e-3_wp)
                
                ! Reduction - Amery ice shelf
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.05_wp,x0=1500.0_wp, y0= 650.0_wp, sigma=300.0_wp,xx=grd%x*1e-3_wp,yy=grd%y*1e-3_wp)
                
                ! Reduction - feeding the Ross ice shelf from the North
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.005_wp,x0=-500.0_wp, y0=-500.0_wp, sigma=300.0_wp,xx=grd%x*1e-3_wp,yy=grd%y*1e-3_wp)

                ! Reduction - feeding the Ross ice shelf from the East
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.05_wp,x0= 130.0_wp, y0=-550.0_wp, sigma=100.0_wp,xx=grd%x*1e-3_wp,yy=grd%y*1e-3_wp)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.05_wp,x0= 280.0_wp, y0=-760.0_wp, sigma=100.0_wp,xx=grd%x*1e-3_wp,yy=grd%y*1e-3_wp)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.05_wp,x0= 380.0_wp, y0=-960.0_wp, sigma=100.0_wp,xx=grd%x*1e-3_wp,yy=grd%y*1e-3_wp)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.05_wp,x0= 400.0_wp, y0=-1150.0_wp,sigma=100.0_wp,xx=grd%x*1e-3_wp,yy=grd%y*1e-3_wp)

end if 

if (.FALSE.) then
                ! Reduction 
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.25_wp,x0=-2000.0_wp,y0=1000.0_wp, sigma=100.0_wp,xx=grd%x*1e-3_wp,yy=grd%y*1e-3_wp)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.20_wp,x0=-750.0_wp, y0=-900.0_wp, sigma=300.0_wp,xx=grd%x*1e-3_wp,yy=grd%y*1e-3_wp)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.20_wp,x0=-600.0_wp, y0=-600.0_wp, sigma=300.0_wp,xx=grd%x*1e-3_wp,yy=grd%y*1e-3_wp)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.20_wp,x0=-300.0_wp, y0=   0.0_wp, sigma=100.0_wp,xx=grd%x*1e-3_wp,yy=grd%y*1e-3_wp)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.20_wp,x0=-250.0_wp, y0=-500.0_wp, sigma=100.0_wp,xx=grd%x*1e-3_wp,yy=grd%y*1e-3_wp)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.20_wp,x0=-100.0_wp, y0=-600.0_wp, sigma=300.0_wp,xx=grd%x*1e-3_wp,yy=grd%y*1e-3_wp)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.20_wp,x0=-100.0_wp, y0=-300.0_wp, sigma=300.0_wp,xx=grd%x*1e-3_wp,yy=grd%y*1e-3_wp)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.20_wp,x0=   0.0_wp, y0=   0.0_wp, sigma=300.0_wp,xx=grd%x*1e-3_wp,yy=grd%y*1e-3_wp)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.20_wp,x0= 130.0_wp, y0=-550.0_wp, sigma= 50.0_wp,xx=grd%x*1e-3_wp,yy=grd%y*1e-3_wp)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.20_wp,x0= 280.0_wp, y0=-760.0_wp, sigma= 50.0_wp,xx=grd%x*1e-3_wp,yy=grd%y*1e-3_wp)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.20_wp,x0= 380.0_wp, y0=-960.0_wp, sigma= 50.0_wp,xx=grd%x*1e-3_wp,yy=grd%y*1e-3_wp)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.50_wp,x0= 400.0_wp, y0=   0.0_wp, sigma=300.0_wp,xx=grd%x*1e-3_wp,yy=grd%y*1e-3_wp)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.20_wp,x0= 400.0_wp, y0=-1150.0_wp,sigma= 50.0_wp,xx=grd%x*1e-3_wp,yy=grd%y*1e-3_wp)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.25_wp,x0= 700.0_wp, y0= -500.0_wp,sigma=400.0_wp,xx=grd%x*1e-3_wp,yy=grd%y*1e-3_wp)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*0.20_wp,x0=1500.0_wp, y0= 650.0_wp, sigma=300.0_wp,xx=grd%x*1e-3_wp,yy=grd%y*1e-3_wp)
                
                ! Increase
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*4.00_wp,x0=-600.0_wp, y0=    0.0_wp,sigma=200.0_wp,xx=grd%x*1e-3_wp,yy=grd%y*1e-3_wp)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*2.00_wp,x0=1200.0_wp, y0=-1200.0_wp,sigma=300.0_wp,xx=grd%x*1e-3_wp,yy=grd%y*1e-3_wp)
                call scale_cf_gaussian(cb_ref,dyn%par%till_cf_ref*1.50_wp,x0=2000.0_wp, y0=    0.0_wp,sigma=300.0_wp,xx=grd%x*1e-3_wp,yy=grd%y*1e-3_wp)
end if 

            end if 

        return 

    end subroutine modify_cb_ref

    subroutine scale_cf_gaussian(cb_ref,cf_new,x0,y0,sigma,xx,yy)

        implicit none 

        real(wp), intent(INOUT) :: cb_ref(:,:)
        real(wp), intent(IN) :: cf_new
        real(wp), intent(IN) :: x0
        real(wp), intent(IN) :: y0
        real(wp), intent(IN) :: sigma
        real(wp), intent(IN) :: xx(:,:)
        real(wp), intent(IN) :: yy(:,:)

        ! Local variables 
        integer :: nx, ny 
        real(wp), allocatable :: wts(:,:)
        
        nx = size(cb_ref,1)
        ny = size(cb_ref,2)

        allocate(wts(nx,ny))

        ! Get Gaussian weights 
        wts = 1.0/(2.0*pi*sigma**2)*exp(-((xx-x0)**2+(yy-y0)**2)/(2.0*sigma**2))
        wts = wts / maxval(wts)

        ! Scale cb_ref
        cb_ref = cb_ref*(1.0-wts) + cf_new*wts

        return 

    end subroutine scale_cf_gaussian

end program yelmo_test



