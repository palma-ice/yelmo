

program yelmo_mismip

    use ncio 
    use yelmo 
    use yelmo_tools, only : stagger_aa_acx, stagger_aa_acy
    use yelmo_dynamics, only: check_vel_convergence
    use deformation 

    use mismip3D 

    implicit none 

    type(yelmo_class)     :: yelmo1

    character(len=56)  :: domain, grid_name  
    character(len=256) :: outfldr, file2D, file1D, file_restart
    character(len=512) :: path_par, path_const 
    character(len=56)  :: experiment, exp_now, res  
    real(prec) :: time_init, time_end, time, dtt, dt2D_out   
    integer    :: n  

    real(prec) :: xmax, ymin, ymax 
    real(prec) :: dx 
    integer    :: i, j, nx, ny 
    real(prec) :: x_gl, x_gl_stnd
    real(prec) :: time_mod_1, time_mod_2 

    integer :: n_att, n_att_tot, q_att, q
    real(prec), allocatable :: ATT_values(:)
    real(prec) :: ATT_time, ATT_dt 
    logical    :: is_converged, exit_loop 
    real(prec) :: err  
    real(prec), allocatable :: ux_bar_prev(:,:), uy_bar_prev(:,:) 

    real(4) :: start, finish

    ! Start timing 
    call cpu_time(start)

    ! Assume program is running from the output folder
    outfldr = "./"

    ! Define input and output locations 
    path_const = trim(outfldr)//"yelmo_const_MISMIP3D.nml"
    path_par   = trim(outfldr)//"yelmo_MISMIP3D.nml" 
    file2D     = trim(outfldr)//"yelmo2D.nc"
    file1D     = trim(outfldr)//"yelmo1D.nc"
    
    ! Define the domain, grid and experiment from parameter file
    call nml_read(path_par,"mismip","domain",       domain)        ! MISMIP3D
    call nml_read(path_par,"mismip","experiment",   experiment)    ! "Std", "RF"
    call nml_read(path_par,"mismip","dx",           dx)            ! [km] Grid resolution ! must be multiple of xmax and ymax!!

    ! Timing parameters 
    call nml_read(path_par,"mismip","time_init",    time_init)     ! [yr] Starting time
    call nml_read(path_par,"mismip","time_end",     time_end)      ! [yr] Ending time
    call nml_read(path_par,"mismip","dtt",          dtt)           ! [yr] Main loop time step 
    call nml_read(path_par,"mismip","dt2D_out",     dt2D_out)      ! [yr] Frequency of 2D output 

    ! Define default grid name for completeness 
    grid_name = "MISMIP3D" 

    ! General initialization of yelmo constants (used globally)
    call yelmo_global_init(path_const)

    ! Set up timing conditions
    if (trim(experiment) .eq. "Stnd") then 

        ! Shorter domain required 
        xmax = 800.0

        ! When to apply time modifications
        time_mod_1 = 15000.0   ! Switch from Stnd => P75S
        time_mod_2 = 15100.0   ! Switch back from P75S => Stnd 
        time_end   = time_mod_2 + 200.0 

        !time_mod_1 = 20000.0
        !time_mod_2 = time_mod_1
        !time_end   = time_mod_1

!         time_end =  1000.0 
!         dt2D_out =   100.0 

    else if (trim(experiment) .eq. "RF") then 
        ! For experiment=="RF" (ratefac)

        ! Longer domain required 
        xmax = 2000.0

        n_att     = 17
        allocate(ATT_values(n_att))
        ATT_values = sec_year*1e-26*[464.16,215.44,100.00,46.416,21.544,10.0,4.6416,2.1544,1.0, &
                                     2.1544,4.6416,10.0,21.544,46.416,100.00,215.44,464.16]
        ATT_time   = 20e3
        ATT_dt     =  2e3 
        time_end   = ATT_time + n_att*ATT_dt + 100e3
        dt2D_out   = 500.0 
        
        write(*,*) "time_init = ", time_init 
        write(*,*) "time_end  = ", time_end
        do q_att = 1, n_att 
            write(*,"(i3,f10.3,g15.5)") q_att, ATT_time+(q_att-1)*ATT_dt, ATT_values(q_att)
        end do
         
    else 

        write(*,*) "experiment not recognized: ", trim(experiment)
        stop 

    end if 
    
    ! Define the domain and grid
    ymax =  50.0
    ymin = -50.0
    call yelmo_init_grid(yelmo1%grd,grid_name,units="km",x0=0.0,dx=dx,nx=int(xmax/dx)+1,y0=ymin,dy=dx,ny=int((ymax-ymin)/dx)+1)

!     call grid_init(grid1,name="MISMIP3D",mtype="cartesian",units="kilometers", &
!                    x0=0.d0,dx=dx,nx=int(xmax/dx)+1,y0=dble(ymin),dy=dx,ny=int((ymax-ymin)/dx)+1)

    ! === Initialize ice sheet model =====

    ! Initialize data objects
    call yelmo_init(yelmo1,filename=path_par,grid_def="none",time=time_init,load_topo=.FALSE., &
                        domain=domain,grid_name=grid_name)

    if (trim(experiment) .eq. "RF") then
        ! Ensure rate factor (and counter) is set for first value if performing RF experiment 
        q_att = 1 
        yelmo1%mat%par%rf_const = ATT_values(q_att)
    end if 

    ! Allocate extra grid vars 
    allocate(ux_bar_prev(yelmo1%grd%nx,yelmo1%grd%ny))
    allocate(uy_bar_prev(yelmo1%grd%nx,yelmo1%grd%ny))

    ! Load boundary values

    yelmo1%bnd%z_sl     = 0.0
    yelmo1%bnd%bmb_shlf = 0.0 
    yelmo1%bnd%T_shlf   = T0  
    yelmo1%bnd%H_sed    = 0.0 
    yelmo1%bnd%H_w      = 0.0

    ! Make sure beta is defined externally and law is compatible 
    yelmo1%dyn%par%beta_method = 1 
    yelmo1%dyn%par%m_drag      = 3 

    call mismip3D_boundaries(yelmo1%bnd%T_srf,yelmo1%bnd%smb,yelmo1%bnd%Q_geo,yelmo1%dyn%now%C_bed, &
                yelmo1%grd%x*1e-3,yelmo1%grd%y*1e-3,x_gl=800.0,experiment=experiment)

    ! Check boundary values 
    call yelmo_print_bound(yelmo1%bnd)

    ! Initialize output file 
    call yelmo_write_init(yelmo1,file2D,time_init=time_init,units="years")
    
    ! Intialize topography 
    call mismip3D_topo_init(yelmo1%bnd%z_bed,yelmo1%tpo%now%H_ice,yelmo1%tpo%now%z_srf, &
                            yelmo1%grd%xc*1e-3,yelmo1%grd%yc*1e-3,experiment)
    
    time     = time_init 
    yelmo1%dyn%par%use_ssa = .TRUE. 

    ! Initialize the yelmo state (dyn,therm,mat)
    call yelmo_init_state(yelmo1,path_par,time=time_init,thrm_method="robin")

    ! Write initial state 
    x_gl      = find_x_gl(yelmo1%grd%x*1e-3,yelmo1%grd%y*1e-3,yelmo1%tpo%now%H_grnd)
    call write_step_2D(yelmo1,file2D,time=time,x_gl=x_gl) 

    ! Assign unrealistic previous values 
    ux_bar_prev = 1e8
    uy_bar_prev = 1e8
    exit_loop   = .FALSE. 

    ! Advance timesteps
    do n = 1, ceiling((time_end-time_init)/dtt)

        ! Get current time 
        time = time_init + n*dtt

        ! == Yelmo ice sheet ===================================================
        call yelmo_update(yelmo1,time)

        x_gl      = find_x_gl(yelmo1%grd%x*1e-3,yelmo1%grd%y*1e-3,yelmo1%tpo%now%H_grnd)
        
        ! == Update boundaries 
            
        if (trim(experiment) .eq. "RF") then 

            exp_now = "RF" 

            is_converged = .FALSE. 
            err = sqrt(sum(yelmo1%tpo%now%dHicedt**2)/yelmo1%grd%npts)
            if (err .lt. 1e-2) is_converged =.TRUE. 

            if (time .gt. ATT_time+ATT_dt) then 
                ! Ensure minimum time per step has been reached before checking convergence

                write(*,*) "err: ", time, ATT_time, err, yelmo1%mat%par%rf_const, q_att 
            
                if (is_converged .and. q_att == n_att) then 
                    ! If output timestep also reached,
                    ! then time to kill simulation 
                    if (mod(time,dt2D_out)==0) exit_loop = .TRUE. 
                else if (is_converged) then
                    ! Time to step ATT_value 
                    q_att = min(q_att+1,n_att)
                    yelmo1%mat%par%rf_const = ATT_values(q_att)
                    ATT_time = time 
                    dt2D_out = 500.0
                end if   

            end if 

        else ! Stnd + P75S 
            exp_now = trim(experiment)
            if (time .gt. time_mod_1) exp_now = "P75S"
            if (time .gt. time_mod_2) exp_now = trim(experiment)
            if (time .gt. time_mod_1-50.0) dt2D_out = 10.0 
            if (time .eq. time_mod_1) x_gl_stnd = x_gl 

        end if  

        call mismip3D_boundaries(yelmo1%bnd%T_srf,yelmo1%bnd%smb,yelmo1%bnd%Q_geo,yelmo1%dyn%now%C_bed, &
                yelmo1%grd%x*1e-3,yelmo1%grd%y*1e-3,x_gl=x_gl_stnd,experiment=exp_now)

        ! == MODEL OUTPUT =======================================================
        if (mod(time,dt2D_out)==0) then  
            call write_step_2D(yelmo1,file2D,time=time,x_gl=x_gl)    
        end if 

        if (mod(time,5.0*dtt)==0) then
            write(*,"(a,2f14.4,a10,g14.3,f10.2)") "time = ",  &
                time, maxval(yelmo1%tpo%now%H_ice), trim(exp_now), yelmo1%mat%par%rf_const, x_gl 
        end if 

        if (exit_loop) exit 

    end do 

    ! Write summary 
    write(*,*) "====== "//trim(domain)//"-"//trim(experiment)//" ======="
    write(*,*) "nz, H0 = ", yelmo1%par%nz_aa, maxval(yelmo1%tpo%now%H_ice)

    ! Finalize program
    call yelmo_end(yelmo1,time=time)

    ! Stop timing 
    call cpu_time(finish)

    print '("Time = ",f12.3," min.")', (finish-start)/60.0 

contains

    subroutine write_step_2D_small(ylmo,filename,time,x_gl)

        implicit none 
        
        type(yelmo_class), intent(IN) :: ylmo
        character(len=*),  intent(IN) :: filename
        real(prec), intent(IN) :: time
        real(prec), intent(IN) :: x_gl 

        ! Local variables
        integer    :: ncid, n
        real(prec) :: time_prev 

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_size(filename,"time",ncid)
        call nc_read(filename,"time",time_prev,start=[n],count=[1],ncid=ncid) 
        if (abs(time-time_prev).gt.1e-5) n = n+1 

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        ! 1D variables of interest 
        call nc_write(filename,"x_rf",ylmo%mat%par%rf_const,units="", &
            long_name="Rate factor",dim1="time",start=[n],count=[1],ncid=ncid)
        call nc_write(filename,"x_gl",x_gl,units="", &
            long_name="Grounding line position",dim1="time",start=[n],count=[1],ncid=ncid)

        ! == yelmo_topography ==
        call nc_write(filename,"H_ice",ylmo%tpo%now%H_ice,units="m",long_name="Ice thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf",ylmo%tpo%now%z_srf,units="m",long_name="Surface elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_bed",ylmo%tpo%now%mask_bed,units="",long_name="Bed mask", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"z_bed",ylmo%bnd%z_bed,units="m",long_name="Bedrock elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_sl",ylmo%bnd%z_sl,units="m",long_name="Sea level rel. to present", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine write_step_2D_small

    subroutine write_step_2D(ylmo,filename,time,x_gl)

        implicit none 
        
        type(yelmo_class), intent(IN) :: ylmo
        character(len=*),  intent(IN) :: filename
        real(prec), intent(IN) :: time
        real(prec), intent(IN) :: x_gl

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

        ! 1D variables of interest 
        call nc_write(filename,"x_rf",ylmo%mat%par%rf_const,units="", &
            long_name="Rate factor",dim1="time",start=[n],count=[1],ncid=ncid)
        call nc_write(filename,"x_gl",x_gl,units="", &
            long_name="Grounding line position",dim1="time",start=[n],count=[1],ncid=ncid)

        ! == yelmo_topography ==
        call nc_write(filename,"H_ice",ylmo%tpo%now%H_ice,units="m",long_name="Ice thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf",ylmo%tpo%now%z_srf,units="m",long_name="Surface elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_bed",ylmo%tpo%now%mask_bed,units="",long_name="Bed mask", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"N_eff",ylmo%tpo%now%N_eff,units="Pa",long_name="Effective pressure", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"dzsrfdt",ylmo%tpo%now%dzsrfdt,units="m/a",long_name="Surface elevation change", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dHicedt",ylmo%tpo%now%dHicedt,units="m/a",long_name="Ice thickness change", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"H_grnd",ylmo%tpo%now%H_grnd,units="m",long_name="Ice thickness overburden", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"dist_grline",ylmo%tpo%now%dist_grline,units="km", &
                      long_name="Distance to nearest grounding-line point", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"is_grline",ylmo%tpo%now%is_grline,units="--", &
                      long_name="Grounding-line point", &
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

        ! == yelmo_material ==
        call nc_write(filename,"visc_bar",ylmo%mat%now%visc_bar,units="Pa a^-1",long_name="Vertically averaged viscosity", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"visc",ylmo%mat%now%visc,units="Pa a^-1",long_name="Viscosity", &
!                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"ATT_bar",ylmo%mat%now%ATT_bar,units="Pa^-3 a^-1",long_name="Vertically averaged rate factor", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"ATT",ylmo%mat%now%ATT,units="Pa^-3 a^-1",long_name="Rate factor", &
!                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        ! == yelmo_dynamics ==

        call nc_write(filename,"ssa_mask_acx",ylmo%dyn%now%ssa_mask_acx,units="1",long_name="SSA mask (acx)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"ssa_mask_acy",ylmo%dyn%now%ssa_mask_acy,units="1",long_name="SSA mask (acy)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"C_bed",ylmo%dyn%now%C_bed,units="m a^-1 Pa^-2",long_name="Dragging constant", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"beta",ylmo%dyn%now%beta,units="Pa m^3 a^-3",long_name="Dragging coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"beta_acx",ylmo%dyn%now%beta_acx,units="Pa m^3 a^-3",long_name="Dragging coefficient (acx)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"beta_acy",ylmo%dyn%now%beta_acy,units="Pa m^3 a^-3",long_name="Dragging coefficient (acy)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dyn_visc_eff",ylmo%dyn%now%visc_eff,units="Pa a^-1",long_name="Vertically averaged viscosity", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"taud_acx",ylmo%dyn%now%taud_acx,units="Pa",long_name="Driving stress, x-direction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taud_acy",ylmo%dyn%now%taud_acy,units="Pa",long_name="Driving stress, y-direction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"sigma_horiz_sq",ylmo%dyn%now%sigma_horiz_sq,units="1",long_name="Horizontal stress components squared", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"lhs_x",ylmo%dyn%now%lhs_x,units="Pa",long_name="Shear reduction (x)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"lhs_y",ylmo%dyn%now%lhs_y,units="Pa",long_name="Shear reduction (y)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"lhs_xy",ylmo%dyn%now%lhs_xy,units="Pa",long_name="Shear reduction magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"duxdz",ylmo%dyn%now%duxdz,units="1/a",long_name="Vertical shear (x)", &
                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"duydz",ylmo%dyn%now%duydz,units="1/a",long_name="Vertical shear (y)", &
                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        
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
        
        call nc_write(filename,"ux_bar",ylmo%dyn%now%ux_bar,units="m/a",long_name="Vertically integrated velocity (x)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uy_bar",ylmo%dyn%now%uy_bar,units="m/a",long_name="Vertically integrated velocity (y)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_bar",ylmo%dyn%now%uxy_bar,units="m/a",long_name="Vertically integrated velocity magnitude", &
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

        call nc_write(filename,"de",ylmo%mat%now%strn%de,units="",long_name="Vertically averaged strain rate", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        ! == yelmo_bound ==

        call nc_write(filename,"z_bed",ylmo%bnd%z_bed,units="m",long_name="Bedrock elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
!         call nc_write(filename,"z_sl",ylmo%bnd%z_sl,units="m",long_name="Sea level rel. to present", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
!         call nc_write(filename,"H_sed",ylmo%bnd%H_sed,units="m",long_name="Sediment thickness", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
!         call nc_write(filename,"H_w",ylmo%bnd%H_w,units="m",long_name="Basal water pressure", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"smb",ylmo%bnd%smb,units="m/a",long_name="Surface mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"T_srf",ylmo%bnd%T_srf,units="K",long_name="Surface temperature", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
     
!         call nc_write(filename,"bmb",ylmo%tpo%now%bmb,units="m/a",long_name="Basal mass balance", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"bmb_shlf",ylmo%bnd%bmb_shlf,units="m/a",long_name="Basal mass balance (shelf)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"Q_geo",ylmo%bnd%Q_geo,units="mW/m^2",long_name="Geothermal heat flux", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine write_step_2D
    
end program yelmo_mismip 



