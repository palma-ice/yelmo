program yelmo_mask_ice
    ! Synthetic stability test for bnd%mask_ice border treatment.
    ! Sets up a small box domain with sinusoidal + noisy bedrock and
    ! a slab of ice with noise, then forces all four border rings to a
    ! single mask_ice value (mask_bnd = -1, 0, or +1) read from the
    ! namelist, after yelmo_init has run. Drives the model forward
    ! and writes 2D/1D output for inspection.

    use nml
    use ncio
    use yelmo
    use timestepping

    implicit none

    type(tstep_class) :: ts
    type(yelmo_class) :: yelmo1

    character(len=56)  :: domain, grid_name
    character(len=256) :: outfldr, file2D, file1D, file_restart
    character(len=512) :: path_par

    ! Timing
    real(wp) :: time_init, time_end, dtt, dt1D_out, dt2D_out

    ! Geometry
    real(wp) :: lx, ly, dx
    integer  :: nx, ny

    ! Bedrock construction
    real(wp) :: zb_mean
    real(wp) :: zb_slope
    real(wp) :: zb_amp_sin
    real(wp) :: zb_wave_x, zb_wave_y
    real(wp) :: zb_amp_noise

    ! Ice and boundary forcing
    real(wp) :: H_ice_init
    real(wp) :: H_ice_noise
    real(wp) :: z_sl_const
    real(wp) :: T_srf_const
    real(wp) :: smb_const
    real(wp) :: Q_geo_const

    ! Mask choice per border ring (west/east/south/north)
    integer  :: mask_bnd_w, mask_bnd_e, mask_bnd_s, mask_bnd_n

    ! Local
    integer  :: i, j
    real(wp), allocatable :: rnd(:,:)
    real(8)  :: cpu_start_time, cpu_end_time, cpu_dtime
    real(wp) :: two_pi

    two_pi = 2.0_wp * pi

    call yelmo_cpu_time(cpu_start_time)

    outfldr = "./"
    call yelmo_load_command_line_args(path_par)

    file2D       = trim(outfldr)//"yelmo2D.nc"
    file1D       = trim(outfldr)//"yelmo1D.nc"
    file_restart = trim(outfldr)//"yelmo_restart.nc"

    ! === Load control parameters ===
    call nml_read(path_par,"ctrl","domain",       domain)

    call nml_read(path_par,"ctrl","time_init",    time_init)
    call nml_read(path_par,"ctrl","time_end",     time_end)
    call nml_read(path_par,"ctrl","dtt",          dtt)
    call nml_read(path_par,"ctrl","dt1D_out",     dt1D_out)
    call nml_read(path_par,"ctrl","dt2D_out",     dt2D_out)

    call nml_read(path_par,"ctrl","dx",           dx)
    call nml_read(path_par,"ctrl","lx",           lx)
    call nml_read(path_par,"ctrl","ly",           ly)

    call nml_read(path_par,"ctrl","zb_mean",      zb_mean)
    call nml_read(path_par,"ctrl","zb_slope",     zb_slope)
    call nml_read(path_par,"ctrl","zb_amp_sin",   zb_amp_sin)
    call nml_read(path_par,"ctrl","zb_wave_x",    zb_wave_x)
    call nml_read(path_par,"ctrl","zb_wave_y",    zb_wave_y)
    call nml_read(path_par,"ctrl","zb_amp_noise", zb_amp_noise)

    call nml_read(path_par,"ctrl","H_ice_init",   H_ice_init)
    call nml_read(path_par,"ctrl","H_ice_noise",  H_ice_noise)
    call nml_read(path_par,"ctrl","z_sl_const",   z_sl_const)
    call nml_read(path_par,"ctrl","T_srf_const",  T_srf_const)
    call nml_read(path_par,"ctrl","smb_const",    smb_const)
    call nml_read(path_par,"ctrl","Q_geo_const",  Q_geo_const)

    call nml_read(path_par,"ctrl","mask_bnd_w",   mask_bnd_w)
    call nml_read(path_par,"ctrl","mask_bnd_e",   mask_bnd_e)
    call nml_read(path_par,"ctrl","mask_bnd_s",   mask_bnd_s)
    call nml_read(path_par,"ctrl","mask_bnd_n",   mask_bnd_n)

    call check_mask_bnd(mask_bnd_w,"mask_bnd_w")
    call check_mask_bnd(mask_bnd_e,"mask_bnd_e")
    call check_mask_bnd(mask_bnd_s,"mask_bnd_s")
    call check_mask_bnd(mask_bnd_n,"mask_bnd_n")

    ! === Grid ===
    grid_name = "MASK_ICE"
    nx = int(lx/dx) + 1
    ny = int(ly/dx) + 1

    call yelmo_init_grid(yelmo1%grd,grid_name,units="km", &
                            x0=-0.5_wp*lx,dx=dx,nx=nx, &
                            y0=-0.5_wp*ly,dy=dx,ny=ny)

    ! === Initialize ice sheet model (no topo loading) ===
    call yelmo_init(yelmo1,filename=path_par,grid_def="none",time=time_init, &
                        load_topo=.FALSE.,domain=domain,grid_name=grid_name)

    nx = yelmo1%grd%nx
    ny = yelmo1%grd%ny

    ! === Define bedrock topography ===
    ! z_bed = zb_mean + slope*x + amp_sin * sin(2*pi*x/Lx_wave) * cos(2*pi*y/Ly_wave)
    !               + amp_noise * rand([-1,+1])

    allocate(rnd(nx,ny))
    call random_seed_init()
    call random_number(rnd)
    rnd = 2.0_wp*rnd - 1.0_wp     ! map [0,1] -> [-1,+1]

    do j = 1, ny
    do i = 1, nx
        yelmo1%bnd%z_bed(i,j) = zb_mean &
            + zb_slope * (yelmo1%grd%x(i,j)*1e-3_wp) &
            + zb_amp_sin * sin(two_pi*yelmo1%grd%x(i,j)*1e-3_wp/zb_wave_x) &
                         * cos(two_pi*yelmo1%grd%y(i,j)*1e-3_wp/zb_wave_y) &
            + zb_amp_noise * rnd(i,j)
    end do
    end do

    ! === Define initial ice thickness ===
    ! H_ice = H_ice_init + H_ice_noise * rand([-1,+1])
    call random_number(rnd)
    rnd = 2.0_wp*rnd - 1.0_wp
    yelmo1%tpo%now%H_ice = max(0.0_wp, H_ice_init + H_ice_noise*rnd)

    ! Surface elevation
    yelmo1%tpo%now%z_srf = yelmo1%bnd%z_bed + yelmo1%tpo%now%H_ice

    ! Reference ice thickness (used when mask_ice == 0)
    yelmo1%bnd%H_ice_ref = yelmo1%tpo%now%H_ice

    ! === Boundary fields ===
    yelmo1%bnd%z_sl     = z_sl_const
    yelmo1%bnd%bmb_shlf = 0.0_wp
    yelmo1%bnd%T_shlf   = yelmo1%bnd%c%T0
    yelmo1%bnd%H_sed    = 0.0_wp
    yelmo1%bnd%T_srf    = yelmo1%bnd%c%T0 + T_srf_const
    yelmo1%bnd%smb      = smb_const
    yelmo1%bnd%Q_geo    = Q_geo_const

    ! === Set mask_ice on the four border rings ===
    ! Interior is left as initialized by ybound_define_mask_ice (= 1 for
    ! the MASK_ICE experiment, which routes through the DEFAULT case).
    yelmo1%bnd%mask_ice(1,:)  = mask_bnd_w      ! west  (x = -lx/2)
    yelmo1%bnd%mask_ice(nx,:) = mask_bnd_e      ! east  (x = +lx/2)
    yelmo1%bnd%mask_ice(:,1)  = mask_bnd_s      ! south (y = -ly/2)
    yelmo1%bnd%mask_ice(:,ny) = mask_bnd_n      ! north (y = +ly/2)

    write(*,*) "yelmo_mask_ice:: domain    = ", trim(domain)
    write(*,*) "yelmo_mask_ice:: nx, ny    = ", nx, ny
    write(*,*) "yelmo_mask_ice:: dx [km]   = ", dx
    write(*,"(a,4i4)") " yelmo_mask_ice:: mask_bnd w/e/s/n = ", &
        mask_bnd_w, mask_bnd_e, mask_bnd_s, mask_bnd_n
    write(*,*) "yelmo_mask_ice:: H_ice min/max = ", minval(yelmo1%tpo%now%H_ice), maxval(yelmo1%tpo%now%H_ice)
    write(*,*) "yelmo_mask_ice:: z_bed min/max = ", minval(yelmo1%bnd%z_bed),     maxval(yelmo1%bnd%z_bed)

    call yelmo_print_bound(yelmo1%bnd)

    ! === Initialize state (dyn, therm, mat) ===
    call yelmo_init_state(yelmo1,time=time_init,thrm_method="robin")

    ! === Initialize timestepping ===
    call tstep_init(ts,time_init,time_end,method="const",units="year", &
                                            time_ref=0.0_wp,const_rel=0.0_wp,const_cal=0.0_wp)

    ! === Write initial state ===
    call yelmo_write_init(yelmo1,file2D,time_init=ts%time,units="years")
    call write_step_2D(yelmo1,file2D,time=ts%time)

    call yelmo_write_reg_init(yelmo1,file1D,time_init=ts%time,units="years", &
                              mask=(yelmo1%bnd%mask_ice /= -1))
    call yelmo_write_reg_step(yelmo1,file1D,time=ts%time)

    ! === Advance timesteps ===
    do while (.not. ts%is_finished)

        call tstep_update(ts,dtt,verbose=.FALSE.)

        call yelmo_update(yelmo1,ts%time)

        if (mod(nint(ts%time_elapsed*100),nint(dt2D_out*100))==0) then
            call write_step_2D(yelmo1,file2D,time=ts%time)
        end if

        if (mod(nint(ts%time_elapsed*100),nint(dt1D_out*100))==0) then
            call yelmo_write_reg_step(yelmo1,file1D,time=ts%time)
        end if

        if (mod(ts%time,100.0)==0 .and. (.not. yelmo_log)) then
            write(*,"(a,f14.4,a,2g14.4)") "yelmo:: time = ", ts%time, &
                "  H_ice min/max = ", minval(yelmo1%tpo%now%H_ice), maxval(yelmo1%tpo%now%H_ice)
        end if

    end do

    ! === Summary ===
    write(*,"(a,a,4i4,a)") " ====== "//trim(domain)//"  mask_bnd w/e/s/n = ", "", &
        mask_bnd_w, mask_bnd_e, mask_bnd_s, mask_bnd_n, "  ======="
    write(*,*) "H_ice min/max = ", minval(yelmo1%tpo%now%H_ice), maxval(yelmo1%tpo%now%H_ice)

    call yelmo_restart_write(yelmo1,file_restart,time=ts%time)
    call yelmo_end(yelmo1,time=ts%time)

    call yelmo_cpu_time(cpu_end_time,cpu_start_time,cpu_dtime)
    write(*,"(a,f12.3,a)") "Time  = ",cpu_dtime/60.0 ," min"
    write(*,"(a,f12.1,a)") "Speed = ",(1e-3*(ts%time_end-ts%time_init))/(cpu_dtime/3600.0), " kiloyears / hr"

contains

    subroutine check_mask_bnd(val,name)
        integer,          intent(IN) :: val
        character(len=*), intent(IN) :: name
        if (val .lt. -1 .or. val .gt. 1) then
            write(*,*) "yelmo_mask_ice:: Error: "//trim(name)//" must be -1, 0, or 1. Got: ", val
            stop
        end if
    end subroutine check_mask_bnd

    subroutine random_seed_init()
        ! Deterministic seed so the noise field is reproducible across runs.
        integer :: n
        integer, allocatable :: seed(:)
        call random_seed(size=n)
        allocate(seed(n))
        seed = 20260513
        call random_seed(put=seed)
    end subroutine random_seed_init

    subroutine write_step_2D(ylmo,filename,time)

        implicit none

        type(yelmo_class), intent(IN) :: ylmo
        character(len=*),  intent(IN) :: filename
        real(wp),          intent(IN) :: time

        integer :: ncid, n

        call nc_open(filename,ncid,writable=.TRUE.)
        n = nc_time_index(filename,"time",time,ncid)
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        call yelmo_write_step_model_metrics(filename,ylmo,n,ncid)

        ! Topography
        call nc_write(filename,"H_ice",ylmo%tpo%now%H_ice,units="m",long_name="Ice thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf",ylmo%tpo%now%z_srf,units="m",long_name="Surface elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_bed",ylmo%bnd%z_bed,units="m",long_name="Bedrock elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"H_grnd",ylmo%tpo%now%H_grnd,units="m",long_name="Grounded ice thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_grnd",ylmo%tpo%now%f_grnd,units="1",long_name="Grounded fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_ice",ylmo%tpo%now%f_ice,units="1",long_name="Ice-covered fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_bed",ylmo%tpo%now%mask_bed,units="",long_name="Bed mask", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! Mass balance components
        call nc_write(filename,"mb_net",ylmo%tpo%now%mb_net,units="m/yr",long_name="Net mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mb_resid",ylmo%tpo%now%mb_resid,units="m/yr",long_name="Residual mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mb_relax",ylmo%tpo%now%mb_relax,units="m/yr",long_name="Relaxation mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dHidt",ylmo%tpo%now%dHidt,units="m/yr",long_name="dH_ice/dt", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"tau_relax",ylmo%tpo%now%tau_relax,units="yr",long_name="Relaxation timescale", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! Boundary fields
        call nc_write(filename,"mask_ice",ylmo%bnd%mask_ice,units="",long_name="Boundary ice mask (-1,0,1)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"H_ice_ref",ylmo%bnd%H_ice_ref,units="m",long_name="Reference ice thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"smb",ylmo%tpo%now%smb,units="m/yr",long_name="Surface mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"T_srf",ylmo%bnd%T_srf,units="K",long_name="Surface temperature", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! Dynamics
        call nc_write(filename,"uxy_bar",ylmo%dyn%now%uxy_bar,units="m/yr",long_name="Vertically averaged velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_b",ylmo%dyn%now%uxy_b,units="m/yr",long_name="Basal sliding velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taud",ylmo%dyn%now%taud,units="Pa",long_name="Driving stress magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_close(ncid)

        return

    end subroutine write_step_2D

end program yelmo_mask_ice
