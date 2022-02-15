

module yelmo_io
    
    use ncio 
    
    use yelmo_defs 
    use yelmo_grid 

    use interp2D
    use coordinates_mapping_scrip, only : map_scrip_class, map_scrip_init, map_scrip_field, &
                                            gen_map_filename

    implicit none

    interface nc_read_interp
        module procedure    nc_read_interp_int_2D
        module procedure    nc_read_interp_int_3D
        module procedure    nc_read_interp_wp_2D
        module procedure    nc_read_interp_wp_3D
        module procedure    nc_read_interp_logical_2D
    end interface

    private 
    public :: yelmo_write_init
    public :: yelmo_write_step
    public :: yelmo_write_step_model_metrics
    public :: yelmo_write_step_pd_metrics
    public :: yelmo_restart_write
    public :: yelmo_restart_read_topo_bnd
    public :: yelmo_restart_read


contains

    subroutine yelmo_write_init(ylmo,filename,time_init,units)

        implicit none 

        type(yelmo_class), intent(IN) :: ylmo 
        character(len=*),  intent(IN) :: filename, units 
        real(wp),          intent(IN) :: time_init
        
        ! Initialize file by writing grid info
        call yelmo_grid_write(ylmo%grd,filename,ylmo%par%domain,ylmo%par%grid_name,create=.TRUE.)


        ! Initialize netcdf file and dimensions
        call nc_write_dim(filename,"month",     x=1,dx=1,nx=12,         units="month")
        call nc_write_dim(filename,"zeta",      x=ylmo%par%zeta_aa,     units="1")
        call nc_write_dim(filename,"zeta_ac",   x=ylmo%par%zeta_ac,     units="1")
        call nc_write_dim(filename,"zeta_rock", x=ylmo%thrm%par%zr%zeta_aa,units="1")
        call nc_write_dim(filename,"age_iso",   x=ylmo%mat%par%age_iso, units="kyr")
        call nc_write_dim(filename,"pd_age_iso",x=ylmo%dta%pd%age_iso,  units="kyr")
        
        call nc_write_dim(filename,"time",      x=time_init,dx=1.0_prec,nx=1,units=trim(units),unlimited=.TRUE.)

        ! Static information
        call nc_write(filename,"basins",  ylmo%bnd%basins, dim1="xc",dim2="yc",units="(0 - 8)",long_name="Hydrological basins")
        call nc_write(filename,"regions", ylmo%bnd%regions,dim1="xc",dim2="yc",units="(0 - 8)",long_name="Domain regions")
        
        ! Additional optional static information 
        call nc_write(filename,"z_bed_sd", ylmo%bnd%z_bed_sd,dim1="xc",dim2="yc",units="m",long_name="Stdev(z_bed)")
        
        return

    end subroutine yelmo_write_init

    subroutine yelmo_write_step(ylmo,filename,time,nms,compare_pd)

        implicit none 
        
        type(yelmo_class), intent(IN) :: ylmo        
        character(len=*),  intent(IN) :: filename
        real(prec),        intent(IN) :: time
        character(len=*),  optional, intent(IN) :: nms(:)
        logical, optional, intent(IN) :: compare_pd

        ! Local variables
        integer    :: ncid, n, q, qtot
        real(prec) :: time_prev 
        character(len=56), allocatable :: names(:) 
        logical ::  write_pd_metrics 

        if (present(nms)) then 
            qtot = size(nms,1)
            allocate(names(qtot))
            do q = 1, qtot 
                names(q) = trim(nms(q))
            end do 
        else 
            qtot = 13 
            allocate(names(qtot))
            names(1)  = "H_ice"
            names(2)  = "z_srf"
            names(3)  = "z_bed"
            names(4)  = "mask_bed"
            names(5)  = "uxy_b"
            names(6)  = "uxy_s"
            names(7)  = "beta"
            names(8)  = "visc_bar"
            names(9)  = "T_prime_b"
            names(10) = "H_w"
            names(11) = "smb"
            names(12) = "bmb"
            names(13) = "z_sl"

        end if 

        write_pd_metrics = .TRUE. 
        if (present(compare_pd)) write_pd_metrics = compare_pd

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
 
        if (write_pd_metrics) then 
            ! Write present-day data metrics (rmse[H],etc)
            call yelmo_write_step_pd_metrics(filename,ylmo,n,ncid)
        end if  
        
        ! Determine number of variables to write 
        qtot = size(names) 

        ! Loop over variables and apply the appropriate write statement
        do q = 1, qtot 

            select case(trim(names(q)))

                case("H_ice")
                    call nc_write(filename,"H_ice",ylmo%tpo%now%H_ice,units="m",long_name="Ice thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
                case("z_srf")
                    call nc_write(filename,"z_srf",ylmo%tpo%now%z_srf,units="m",long_name="Surface elevation", &
                                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
                case("z_bed")
                    call nc_write(filename,"z_bed",ylmo%bnd%z_bed,units="m",long_name="Bedrock elevation", &
                                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid) 
                case("mask_bed")
                    call nc_write(filename,"mask_bed",ylmo%tpo%now%mask_bed,units="",long_name="Bed mask", &
                                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
                case("uxy_b")
                    call nc_write(filename,"uxy_b",ylmo%dyn%now%uxy_b,units="m/yr",long_name="Basal sliding velocity magnitude", &
                                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
                case("uxy_s")
                    call nc_write(filename,"uxy_s",ylmo%dyn%now%uxy_s,units="m/yr",long_name="Surface velocity magnitude", &
                                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
                case("beta")
                    call nc_write(filename,"beta",ylmo%dyn%now%beta,units="Pa yr m-1",long_name="Basal friction coefficient", &
                                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
                case("visc_bar")
                    call nc_write(filename,"visc_bar",ylmo%mat%now%visc_bar,units="Pa yr",long_name="Vertically averaged viscosity", &
                                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
                case("T_prime_b")
                    call nc_write(filename,"T_prime_b",ylmo%thrm%now%T_prime_b,units="K",long_name="Basal homologous ice temperature", &
                                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
                case("H_w")
                    call nc_write(filename,"H_w",ylmo%thrm%now%H_w,units="m water equiv.",long_name="Basal water layer thickness", &
                                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
                case("smb")
                    call nc_write(filename,"smb",ylmo%bnd%smb,units="m/yr ice equiv.",long_name="Surface mass balance", &
                                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
                case("bmb")
                    call nc_write(filename,"bmb",ylmo%tpo%now%bmb,units="m/yr ice equiv.",long_name="Basal mass balance", &
                                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
                case("z_sl")
                    call nc_write(filename,"z_sl",ylmo%bnd%z_sl,units="m",long_name="Sea level rel. to present", &
                                    dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
                
                case DEFAULT 

                    write(io_unit_err,*) 
                    write(io_unit_err,*) "yelmo_write_step:: Error: variable not yet supported."
                    write(io_unit_err,*) "variable = ", trim(names(q))
                    write(io_unit_err,*) "filename = ", trim(filename)
                    stop 
                    
            end select


        end do 
            
        ! External data
        ! call nc_write(filename,"dzbdt",isos%now%dzbdt,units="m/a",long_name="Bedrock uplift rate", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! call nc_write(filename,"Ta_ann",snp%now%ta_ann,units="K",long_name="Near-surface air temperature (ann)", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        ! call nc_write(filename,"Ta_sum",snp%now%ta_sum,units="K",long_name="Near-surface air temperature (sum)", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        ! call nc_write(filename,"Pr_ann",snp%now%pr_ann*1e-3,units="m/a water equiv.",long_name="Precipitation (ann)", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! call nc_write(filename,"dTa_ann",snp%now%ta_ann-snp%clim0%ta_ann,units="K",long_name="Near-surface air temperature anomaly (ann)", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        ! call nc_write(filename,"dTa_sum",snp%now%ta_sum-snp%clim0%ta_sum,units="K",long_name="Near-surface air temperature anomaly (sum)", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        ! call nc_write(filename,"dPr_ann",(snp%now%pr_ann-snp%clim0%pr_ann)*1e-3,units="m/a water equiv.",long_name="Precipitation anomaly (ann)", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        ! call nc_write(filename,"dT_shlf",mshlf%now%dT_shlf,units="K",long_name="Shelf temperature anomaly", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine yelmo_write_step

    subroutine yelmo_write_step_model_metrics(filename,ylmo,n,ncid)
        ! Write model performance metrics (speed, dt, eta) 

        implicit none 

        character(len=*),  intent(IN) :: filename 
        type(yelmo_class), intent(IN) :: ylmo 
        integer                       :: n 
        integer, optional             :: ncid 
        
        ! Write model speed 
        call nc_write(filename,"speed",ylmo%time%model_speed,units="kyr/hr",long_name="Model speed (Yelmo only)", &
                      dim1="time",start=[n],count=[1],missing_value=mv,ncid=ncid)
        call nc_write(filename,"dt_avg",ylmo%time%dt_avg,units="yr",long_name="Average timestep", &
                      dim1="time",start=[n],count=[1],missing_value=mv,ncid=ncid)
        call nc_write(filename,"eta_avg",ylmo%time%eta_avg,units="m a**-1",long_name="Average eta (maximum PC truncation error)", &
                      dim1="time",start=[n],count=[1],missing_value=mv,ncid=ncid)
        call nc_write(filename,"ssa_iter_avg",ylmo%time%ssa_iter_avg,units="",long_name="Average Picard iterations for SSA convergence", &
                      dim1="time",start=[n],count=[1],missing_value=mv,ncid=ncid)

        call nc_write(filename,"pc_tau_max",abs(ylmo%time%pc_tau_max),units="m a**-1", &
                        long_name="Maximum truncation error over last N timestep (magnitude)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],count=[ylmo%grd%nx,ylmo%grd%ny,1], &
                      missing_value=mv,ncid=ncid)
        
        return 

    end subroutine yelmo_write_step_model_metrics

    subroutine yelmo_write_step_pd_metrics(filename,ylmo,n,ncid)
        ! Write present-day data comparison metrics (rmse[H],etc)

        implicit none 

        character(len=*),  intent(IN) :: filename 
        type(yelmo_class), intent(IN) :: ylmo 
        integer                       :: n 
        integer, optional             :: ncid 
        
        call nc_write(filename,"rmse_H",ylmo%dta%pd%rmse_H,units="m",long_name="RMSE - Ice thickness", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"rmse_zsrf",ylmo%dta%pd%rmse_zsrf,units="m",long_name="RMSE - Surface elevation", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"rmse_uxy",ylmo%dta%pd%rmse_uxy,units="m/a",long_name="RMSE - Surface velocity", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"rmse_uxy_log",ylmo%dta%pd%rmse_loguxy,units="log(m/a)",long_name="RMSE - Log surface velocity", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"rmse_iso",ylmo%dta%pd%rmse_iso,units="m",long_name="RMSE - isochronal layer depths", &
                  dim1="pd_age_iso",dim2="time",start=[1,n],missing_value=mv,ncid=ncid)

        return 

    end subroutine yelmo_write_step_pd_metrics

    subroutine yelmo_restart_write(dom,filename,time,init)

        implicit none 

        type(yelmo_class), intent(IN) :: dom
        character(len=*),  intent(IN) :: filename 
        real(wp),          intent(IN) :: time 
        logical, optional, intent(IN) :: init 

        ! Local variables
        integer  :: ncid, n, nx, ny, nz, nz_ac, nz_r   
        logical  :: initialize_file 
        real(wp) :: time_prev 
        character(len=2) :: xnm, ynm 

        initialize_file = .TRUE. 
        if (present(init)) initialize_file = init

        xnm = "xc"
        ynm = "yc" 

        nx    = dom%grd%nx
        ny    = dom%grd%ny 
        nz    = size(dom%par%zeta_aa,1)
        nz_ac = size(dom%par%zeta_ac,1)
        nz_r  = size(dom%thrm%par%zr%zeta_aa,1)

        ! Write all yelmo data to file, so that it can be
        ! read later to restart a simulation.
        
        ! == Initialize netcdf file ==============================================

        if (initialize_file) then

            ! Initialize file by writing grid info
            call yelmo_grid_write(dom%grd,filename,dom%par%domain,dom%par%grid_name,create=.TRUE.)

            call nc_write_dim(filename,"zeta",     x=dom%par%zeta_aa,       units="1")
            call nc_write_dim(filename,"zeta_ac",  x=dom%par%zeta_ac,       units="1")
            call nc_write_dim(filename,"zeta_rock",x=dom%thrm%par%zr%zeta_aa,units="1")
            call nc_write_dim(filename,"age_iso",  x=dom%mat%par%age_iso,   units="kyr")
            call nc_write_dim(filename,"time",     x=time,dx=1.0_prec,nx=1, units="years")
            
            call nc_write_dim(filename,"pc_steps", x=[1,2,3],units="1")
            
            if (dom%grd%is_projection) then 
                call nc_write_attr(filename,xnm,"standard_name","projection_x_coordinate")
                call nc_write_attr(filename,ynm,"standard_name","projection_y_coordinate")
            end if 

        end if 
        
        ! == Begin writing data ==============================================
        
        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)
        
        if (initialize_file) then 

            ! Current time index to write will be the first and only one 
            n = 1 

        else 

            ! Determine current writing time step 
            n = nc_size(filename,"time",ncid)
            call nc_read(filename,"time",time_prev,start=[n],count=[1],ncid=ncid) 
            if (abs(time-time_prev).gt.1e-5) n = n+1 

            ! Update the time step
            call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        end if 

        ! == Predictor-corrector (pc) variables ===
        ! (these will not be read in by yelmo_restart_read, but can be useful to output for diagnostics)

        call nc_write(filename,"pc_tau",       dom%time%pc_tau,        units="m/yr",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"pc_tau_masked",dom%time%pc_tau_masked, units="m/yr",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"pc_tau_max",   dom%time%pc_tau_max,    units="m/yr",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        
        ! == time variables ===

        ! Note: these variables are not defined on the 2D grid, so they 
        ! will give interpolation errors to `cdo remapcon`. They 
        ! are currently not used by a restarted simulation, so they 
        ! are commented out here. However, it could be useful in the future
        ! to have these variables in the restart file, and potentially 
        ! others that do not depend on the grid.
        
        ! call nc_write(filename,"pc_dt",        dom%time%pc_dt,         units="yr",  dim1="pc_steps",dim2="time",ncid=ncid,start=[1,n],count=[3,1],grid_mapping="")
        ! call nc_write(filename,"pc_eta",       dom%time%pc_eta,        units="m/yr",dim1="pc_steps",dim2="time",ncid=ncid,start=[1,n],count=[3,1],grid_mapping="")
        
        ! == ytopo variables ===
        call nc_write(filename,"H_ice",       dom%tpo%now%H_ice,       units="m",  dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"z_srf",       dom%tpo%now%z_srf,       units="m",  dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])

        call nc_write(filename,"dzsrfdt",     dom%tpo%now%dzsrfdt,     units="m/yr",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"dHicedt",     dom%tpo%now%dHicedt,     units="m/yr",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"fmb",         dom%tpo%now%fmb,         units="m/yr",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"bmb",         dom%tpo%now%bmb,         units="m/yr",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"mb_applied",  dom%tpo%now%mb_applied,  units="m/yr",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"mb_resid",    dom%tpo%now%mb_resid,    units="m/yr",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"calv_flt",    dom%tpo%now%calv_flt,    units="m/yr",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"calv_grnd",   dom%tpo%now%calv_grnd,   units="m/yr",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"calv",        dom%tpo%now%calv,        units="m/yr",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"eps_eff",     dom%tpo%now%eps_eff,     units="1/yr",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"tau_eff",     dom%tpo%now%tau_eff,     units="Pa",  dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"dzsdx",       dom%tpo%now%dzsdx,       units="m/m", dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])  
        call nc_write(filename,"dzsdy",       dom%tpo%now%dzsdy,       units="m/m", dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])  
        call nc_write(filename,"dHicedx",     dom%tpo%now%dHicedx,     units="m/m", dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"dHicedy",     dom%tpo%now%dHicedy,     units="m/m", dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"H_eff",       dom%tpo%now%H_eff,       units="m",   dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"H_grnd",      dom%tpo%now%H_grnd,      units="b",   dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"f_grnd",      dom%tpo%now%f_grnd,      units="1",   dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])  
        call nc_write(filename,"f_grnd_acx",  dom%tpo%now%f_grnd_acx,  units="1",   dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"f_grnd_acy",  dom%tpo%now%f_grnd_acy,  units="1",   dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"f_ice",       dom%tpo%now%f_ice,       units="1",   dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])  
        call nc_write(filename,"dist_margin", dom%tpo%now%dist_margin, units="m",   dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"dist_grline", dom%tpo%now%dist_grline, units="m",   dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"mask_bed",    dom%tpo%now%mask_bed,    units="1",   dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])  
        call nc_write(filename,"mask_grz",    dom%tpo%now%mask_grz,    units="1",   dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])  
            
        call nc_write(filename,"dHdt_n",      dom%tpo%now%dHdt_n,      units="m/yr",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"H_ice_n",     dom%tpo%now%H_ice_n,     units="m",   dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"H_ice_pred",  dom%tpo%now%H_ice_pred,  units="m",   dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"H_ice_corr",  dom%tpo%now%H_ice_corr,  units="m",   dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        
        call nc_write(filename,"z_srf_n",     dom%tpo%now%z_srf_n,     units="m",  dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        
        ! == ydyn variables ===

        call nc_write(filename,"ux",            dom%dyn%now%ux,     units="m/a",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1]) 
        call nc_write(filename,"uy",            dom%dyn%now%uy,     units="m/a",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1]) 
        call nc_write(filename,"uxy",           dom%dyn%now%uxy,    units="m/a",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1]) 
        call nc_write(filename,"uz",            dom%dyn%now%uz,     units="m/a",dim1="xc",dim2="yc",dim3="zeta_ac",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz_ac,1]) 
      
        call nc_write(filename,"ux_bar",        dom%dyn%now%ux_bar, units="m/a",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"uy_bar",        dom%dyn%now%uy_bar, units="m/a",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"uxy_bar",       dom%dyn%now%uxy_bar,units="m/a",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 

        call nc_write(filename,"ux_bar_prev",   dom%dyn%now%ux_bar_prev, units="m/a",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"uy_bar_prev",   dom%dyn%now%uy_bar_prev, units="m/a",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        
        call nc_write(filename,"ux_b",          dom%dyn%now%ux_b,   units="m/a",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"uy_b",          dom%dyn%now%uy_b,   units="m/a",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"uxy_b",         dom%dyn%now%uxy_b,  units="m/a",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 

        call nc_write(filename,"ux_s",          dom%dyn%now%ux_s,     units="m/a",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"uy_s",          dom%dyn%now%uy_s,     units="m/a",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"uxy_s",         dom%dyn%now%uxy_s,    units="m/a",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 

        call nc_write(filename,"ux_i",          dom%dyn%now%ux_i,     units="m/a",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1]) 
        call nc_write(filename,"uy_i",          dom%dyn%now%uy_i,     units="m/a",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1]) 
        call nc_write(filename,"ux_i_bar",      dom%dyn%now%ux_i_bar, units="m/a",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"uy_i_bar",      dom%dyn%now%uy_i_bar, units="m/a",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"uxy_i_bar",     dom%dyn%now%uxy_i_bar,units="m/a",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 

        call nc_write(filename,"duxydt",        dom%dyn%now%duxydt,   units="m/a^2",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 

        call nc_write(filename,"duxdz",         dom%dyn%now%duxdz,     units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1]) 
        call nc_write(filename,"duydz",         dom%dyn%now%duydz,     units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1]) 
        call nc_write(filename,"duxdz_bar",     dom%dyn%now%duxdz_bar, units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"duydz_bar",     dom%dyn%now%duydz_bar, units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 

        call nc_write(filename,"taud_acx",      dom%dyn%now%taud_acx, units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"taud_acy",      dom%dyn%now%taud_acy, units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"taud",          dom%dyn%now%taud,     units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 

        call nc_write(filename,"taub_acx",      dom%dyn%now%taub_acx, units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"taub_acy",      dom%dyn%now%taub_acy, units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"taub",          dom%dyn%now%taub,     units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 

        call nc_write(filename,"qq_gl_acx",     dom%dyn%now%qq_gl_acx,units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"qq_gl_acy",     dom%dyn%now%qq_gl_acy,units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        
        call nc_write(filename,"qq_acx",        dom%dyn%now%qq_acx,   units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"qq_acy",        dom%dyn%now%qq_acy,   units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"qq",            dom%dyn%now%qq,       units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 

        call nc_write(filename,"visc_eff",      dom%dyn%now%visc_eff,    units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1]) 
        call nc_write(filename,"visc_eff_int",  dom%dyn%now%visc_eff_int,units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 

        call nc_write(filename,"N_eff",         dom%dyn%now%N_eff,    units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])        
        call nc_write(filename,"cb_ref",        dom%dyn%now%cb_ref,   units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"c_bed",         dom%dyn%now%c_bed,    units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"beta_acx",      dom%dyn%now%beta_acx, units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"beta_acy",      dom%dyn%now%beta_acy, units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"beta",          dom%dyn%now%beta,     units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 

        call nc_write(filename,"beta_eff",      dom%dyn%now%beta_eff, units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        
        call nc_write(filename,"f_vbvs",        dom%dyn%now%f_vbvs,   units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 

        call nc_write(filename,"ssa_mask_acx",  dom%dyn%now%ssa_mask_acx, units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"ssa_mask_acy",  dom%dyn%now%ssa_mask_acy, units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"ssa_err_acx",   dom%dyn%now%ssa_err_acx,  units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"ssa_err_acy",   dom%dyn%now%ssa_err_acy,  units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 

        ! == ymat variables === 

        call nc_write(filename,"strn2D_dxx", dom%mat%now%strn2D%dxx,  units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"strn2D_dyy", dom%mat%now%strn2D%dyy,  units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"strn2D_dxy", dom%mat%now%strn2D%dxy,  units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"strn2D_dxz", dom%mat%now%strn2D%dxz,  units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"strn2D_dyz", dom%mat%now%strn2D%dyz,  units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"strn2D_de",  dom%mat%now%strn2D%de,   units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"strn2D_div", dom%mat%now%strn2D%div,  units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"strn2D_f_shear",dom%mat%now%strn2D%f_shear,units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        
        call nc_write(filename,"strn_dxx",     dom%mat%now%strn%dxx,        units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1]) 
        call nc_write(filename,"strn_dyy",     dom%mat%now%strn%dyy,        units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1]) 
        call nc_write(filename,"strn_dxy",     dom%mat%now%strn%dxy,        units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1]) 
        call nc_write(filename,"strn_dxz",     dom%mat%now%strn%dxz,        units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1]) 
        call nc_write(filename,"strn_dyz",     dom%mat%now%strn%dyz,        units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1]) 
        call nc_write(filename,"strn_de",      dom%mat%now%strn%de,         units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1]) 
        call nc_write(filename,"strn_div",     dom%mat%now%strn%div,        units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1]) 
        call nc_write(filename,"strn_f_shear", dom%mat%now%strn%f_shear,    units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1]) 

        call nc_write(filename,"strs2D_txx", dom%mat%now%strs2D%txx,  units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"strs2D_tyy", dom%mat%now%strs2D%tyy,  units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"strs2D_txy", dom%mat%now%strs2D%txy,  units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"strs2D_txz", dom%mat%now%strs2D%txz,  units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"strs2D_tyz", dom%mat%now%strs2D%tyz,  units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"strs2D_te",  dom%mat%now%strs2D%te,   units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"strs2D_tau_eig_1",dom%mat%now%strs2D%tau_eig_1,units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"strs2D_tau_eig_2",dom%mat%now%strs2D%tau_eig_2,units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        
        call nc_write(filename,"strs_txx",     dom%mat%now%strs%txx,        units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1]) 
        call nc_write(filename,"strs_tyy",     dom%mat%now%strs%tyy,        units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1]) 
        call nc_write(filename,"strs_txy",     dom%mat%now%strs%txy,        units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1]) 
        call nc_write(filename,"strs_txz",     dom%mat%now%strs%txz,        units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1]) 
        call nc_write(filename,"strs_tyz",     dom%mat%now%strs%tyz,        units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1]) 
        call nc_write(filename,"strs_te",      dom%mat%now%strs%te,         units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1]) 

        call nc_write(filename,"enh",         dom%mat%now%enh,           units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1]) 
        call nc_write(filename,"enh_bnd",     dom%mat%now%enh_bnd,       units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1]) 
        call nc_write(filename,"enh_bar",     dom%mat%now%enh_bar,       units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"ATT",         dom%mat%now%ATT,           units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1]) 
        call nc_write(filename,"ATT_bar",     dom%mat%now%ATT_bar,       units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 
        call nc_write(filename,"visc",        dom%mat%now%visc,          units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1]) 
        call nc_write(filename,"visc_int",    dom%mat%now%visc_int,      units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 

        call nc_write(filename,"f_shear_bar", dom%mat%now%f_shear_bar,   units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1]) 

        call nc_write(filename,"dep_time",    dom%mat%now%dep_time,      units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1])     
        call nc_write(filename,"depth_iso",   dom%mat%now%depth_iso,     units="",dim1="xc",dim2="yc",dim3="age_iso",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,size(dom%mat%now%depth_iso,3),1])     
        
        ! == ytherm variables ===

        call nc_write(filename,"enth",        dom%thrm%now%enth,       units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1])      
        call nc_write(filename,"T_ice",       dom%thrm%now%T_ice,      units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1])      
        call nc_write(filename,"omega",       dom%thrm%now%omega,      units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1])      
        call nc_write(filename,"T_pmp",       dom%thrm%now%T_pmp,      units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1])      
        
        call nc_write(filename,"f_pmp",       dom%thrm%now%f_pmp,      units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])        
        call nc_write(filename,"bmb_grnd",    dom%thrm%now%bmb_grnd,   units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])     
        call nc_write(filename,"Q_strn",      dom%thrm%now%Q_strn,     units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1])     
        call nc_write(filename,"Q_b",         dom%thrm%now%Q_b,        units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])          
        call nc_write(filename,"Q_ice_b",     dom%thrm%now%Q_ice_b,    units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])          
        call nc_write(filename,"T_prime_b",   dom%thrm%now%T_prime_b,  units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])    
        call nc_write(filename,"H_w",         dom%thrm%now%H_w,        units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"dHwdt",       dom%thrm%now%dHwdt,      units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        
        call nc_write(filename,"cp",          dom%thrm%now%cp,         units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1])         
        call nc_write(filename,"kt",          dom%thrm%now%kt,         units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1])         
        call nc_write(filename,"H_cts",       dom%thrm%now%H_cts,      units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])      
        
        call nc_write(filename,"advecxy",     dom%thrm%now%advecxy,    units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1])      
        
        call nc_write(filename,"Q_rock",      dom%thrm%now%Q_rock,     units="mW m-2",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])        
        call nc_write(filename,"enth_rock",   dom%thrm%now%enth_rock,  units="J m-3", dim1="xc",dim2="yc",dim3="zeta_rock",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz_r,1])      
        call nc_write(filename,"T_rock",      dom%thrm%now%T_rock,     units="K",     dim1="xc",dim2="yc",dim3="zeta_rock",dim4="time",ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz_r,1])      

        ! == ybound variables ===

        call nc_write(filename,"z_bed",       dom%bnd%z_bed,       units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"z_bed_sd",    dom%bnd%z_bed_sd,    units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"z_sl",        dom%bnd%z_sl,        units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"H_sed",       dom%bnd%H_sed,       units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"smb",         dom%bnd%smb,         units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"T_srf",       dom%bnd%T_srf,       units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"bmb_shlf",    dom%bnd%bmb_shlf,    units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"T_shlf",      dom%bnd%T_shlf,      units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"Q_geo",       dom%bnd%Q_geo,       units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])

        call nc_write(filename,"enh_srf",     dom%bnd%enh_srf,     units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])

        call nc_write(filename,"basins",      dom%bnd%basins,      units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"basin_mask",  dom%bnd%basin_mask,  units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"regions",     dom%bnd%regions,     units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"region_mask", dom%bnd%region_mask, units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])

        call nc_write(filename,"ice_allowed", dom%bnd%ice_allowed, units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"calv_mask",   dom%bnd%calv_mask,   units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"H_ice_ref",   dom%bnd%H_ice_ref,   units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])
        call nc_write(filename,"z_bed_ref",   dom%bnd%z_bed_ref,   units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n],count=[nx,ny,1])

        ! == ydata variables ===

        ! TO DO (not necessary for restart, but let's see...)

        ! Close the netcdf file
        call nc_close(ncid)

        ! Write summary 
        write(*,*) 
        write(*,*) "time = ", time, " : saved restart file: ", trim(filename)
        write(*,*) 

        return 

    end subroutine yelmo_restart_write

    subroutine yelmo_read_interp_2D(var2D,filename,vname,domain,grid_name)  
        ! Load a variable from a file.
        ! Interpolate to current grid as needed. 
        
        implicit none 

        real(wp), intent(OUT) :: var2D(:,:)  
        character(len=*),  intent(IN)    :: filename 
        character(len=*),  intent(IN)    :: vname
        character(len=*),  intent(IN)    :: domain 
        character(len=*),  intent(IN)    :: grid_name  
        
        ! Local variables
        integer :: nx, ny, n 
        character(len=56) :: file_domain 
        character(len=56) :: file_grid_name 
        type(map_scrip_class) :: mps 

        nx = size(var2D,1)
        ny = size(var2D,2) 

        ! Load restart file grid attributes 
        if (nc_exists_attr(filename,"domain")) then 
            call nc_read_attr(filename,"domain",    file_domain)
        else 
            file_domain = trim(domain)
        end if 

        if (nc_exists_attr(filename,"grid_name")) then 
            call nc_read_attr(filename,"grid_name", file_grid_name)
        else 
            file_grid_name = trim(grid_name)
        end if 

        ! Determine which slice to get (last one)
        n = nc_size(filename,"time")

        if (trim(file_grid_name) .eq. trim(grid_name) ) then 
            ! File's grid and yelmo grid are the same

            ! Load the data without interpolation (by not specifying mps argument)
            call nc_read(filename,vname,var2D,start=[1,1,n],count=[nx,ny,1])

        else 
            ! Restart grid is different than Yelmo grid 

            ! ! Load the scrip map from file (should already have been generated via cdo externally)
            ! call map_scrip_init(mps,restart_grid_name,dom%par%grid_name, &
            !                         method="con",fldr="maps",load=.TRUE.)

            ! call yelmo_read_interp_internal(dom,filename,time,mps) 
            
        end if 

        return 

    end subroutine yelmo_read_interp_2D

    subroutine yelmo_restart_read_topo_bnd(dom,filename,time)  
        ! Load yelmo variables from restart file: [tpo] 
        ! [dyn,therm,mat] variables loaded using yelmo_restart_read
        
        implicit none 

        type(yelmo_class), intent(INOUT) :: dom 
        character(len=*),  intent(IN)    :: filename 
        real(wp),          intent(IN)    :: time 
        
        ! Local variables
        character(len=56) :: restart_domain 
        character(len=56) :: restart_grid_name 
        type(map_scrip_class) :: mps 

        ! Load restart file grid attributes 
        if (nc_exists_attr(filename,"domain")) then 
            call nc_read_attr(filename,"domain",    restart_domain)
        else 
            restart_domain = trim(dom%par%domain)
        end if 

        if (nc_exists_attr(filename,"grid_name")) then 
            call nc_read_attr(filename,"grid_name", restart_grid_name)
        else 
            restart_grid_name = trim(dom%par%grid_name)
        end if 


        if (trim(restart_grid_name) .eq. trim(dom%par%grid_name) ) then 
            ! Restart file grid and yelmo grid are the same

            ! Load the data without interpolation (by not specifying mps argument)
            call yelmo_restart_read_topo_bnd_internal(dom,filename,time)

            ! Set yelmo flag too
            dom%par%restart_interpolated = .FALSE. 

        else 
            ! Restart grid is different than Yelmo grid 

            ! Load the scrip map from file (should already have been generated via cdo externally)
            call map_scrip_init(mps,restart_grid_name,dom%par%grid_name, &
                                    method="con",fldr="maps",load=.TRUE.)

            call yelmo_restart_read_topo_bnd_internal(dom,filename,time,mps) 

            ! Set yelmo flag too
            dom%par%restart_interpolated = .TRUE. 
            
        end if 

        return 

    end subroutine yelmo_restart_read_topo_bnd

    subroutine yelmo_restart_read(dom,filename,time)
        ! Load yelmo variables from restart file: [dyn,therm,mat] 
        ! [tpo] variables loaded using yelmo_restart_read_topo

        implicit none 

        type(yelmo_class), intent(INOUT) :: dom 
        character(len=*),  intent(IN)    :: filename 
        real(wp),          intent(IN)    :: time 
        
        ! Local variables 
        character(len=56) :: restart_domain 
        character(len=56) :: restart_grid_name 
        type(map_scrip_class) :: mps 
        
        ! Load restart file grid attributes 
        if (nc_exists_attr(filename,"domain")) then 
            call nc_read_attr(filename,"domain",    restart_domain)
        else 
            restart_domain = trim(dom%par%domain)
        end if 

        if (nc_exists_attr(filename,"grid_name")) then 
            call nc_read_attr(filename,"grid_name", restart_grid_name)
        else 
            restart_grid_name = trim(dom%par%grid_name)
        end if 
        
        
        if (trim(restart_grid_name) .eq. trim(dom%par%grid_name) ) then 
            ! Restart file grid and yelmo grid are the same

            ! Load the data without interpolation (by not specifying mps argument)
            call yelmo_restart_read_internal(dom,filename,time)

        else 
            ! Restart grid is different than Yelmo grid 

            ! Load the scrip map from file (should already have been generated via cdo externally)
            call map_scrip_init(mps,restart_grid_name,dom%par%grid_name, &
                                    method="con",fldr="maps",load=.TRUE.)

            call yelmo_restart_read_internal(dom,filename,time,mps) 

        end if 
        
        ! ajr: testing
        call yelmo_restart_write(dom,"yelmo_restart_init.nc",time)
        
        return 

    end subroutine yelmo_restart_read
    
    subroutine yelmo_restart_read_topo_bnd_internal(dom,filename,time,mps)  
        ! Load yelmo variables from restart file: [tpo] 
        ! [dyn,therm,mat] variables loaded using yelmo_restart_read
        
        implicit none 

        type(yelmo_class), intent(INOUT) :: dom 
        character(len=*),  intent(IN)    :: filename 
        real(wp),          intent(IN)    :: time 
        type(map_scrip_class), optional, intent(IN) :: mps 

        ! Local variables
        integer  :: ncid, n, nx, ny, nz, nz_ac 
        real(wp) :: time_of_restart_file 

        ! Read all yelmo data from file,
        ! in order to restart a simulation.
        
        ! Open the file for reading

        call nc_open(filename,ncid,writable=.FALSE.)

        if (.FALSE.) then 
            ! No need to read in the dimension information,
            ! this will be initialized by Yelmo itself

            call nc_read(filename,"xc",dom%grd%xc,ncid=ncid)
            dom%grd%xc = dom%grd%xc*1e3      
            call nc_read(filename,"yc",dom%grd%yc,ncid=ncid)   
            dom%grd%yc = dom%grd%yc*1e3   
            call nc_read(filename,"zeta",dom%par%zeta_aa,ncid=ncid)   
            call nc_read(filename,"zeta_ac",dom%par%zeta_ac,ncid=ncid) 
            call nc_read(filename,"time",time_of_restart_file,ncid=ncid)    

        end if 
        
        ! Define dimensions of variables 
        nx    = size(dom%grd%xc,1)
        ny    = size(dom%grd%yc,1)
        nz    = size(dom%par%zeta_aa,1) 
        nz_ac = size(dom%par%zeta_ac,1) 
        
        ! Assume that first time dimension value is to be read in
        n = 1 

        ! == time variables ===

        ! call nc_read(filename,"pc_dt",       dom%time%pc_dt, ncid=ncid)
        ! call nc_read(filename,"pc_eta",      dom%time%pc_eta,ncid=ncid)
        
        ! == ytopo variables ===

        call nc_read_interp(filename,"H_ice",       dom%tpo%now%H_ice,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"z_srf",       dom%tpo%now%z_srf,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"dzsrfdt",     dom%tpo%now%dzsrfdt,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"dHicedt",     dom%tpo%now%dHicedt,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"fmb",         dom%tpo%now%fmb,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"bmb",         dom%tpo%now%bmb,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"mb_applied",  dom%tpo%now%mb_applied,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"mb_resid",    dom%tpo%now%mb_resid,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"calv_flt",    dom%tpo%now%calv_flt,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"calv_grnd",   dom%tpo%now%calv_grnd,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"calv",        dom%tpo%now%calv,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"eps_eff",     dom%tpo%now%eps_eff,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"tau_eff",     dom%tpo%now%tau_eff,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"dzsdx",       dom%tpo%now%dzsdx,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)  
        call nc_read_interp(filename,"dzsdy",       dom%tpo%now%dzsdy,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)  
        call nc_read_interp(filename,"dHicedx",     dom%tpo%now%dHicedx,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"dHicedy",     dom%tpo%now%dHicedy,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"H_eff",       dom%tpo%now%H_eff,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"H_grnd",      dom%tpo%now%H_grnd,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"f_grnd",      dom%tpo%now%f_grnd,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)  
        call nc_read_interp(filename,"f_grnd_acx",  dom%tpo%now%f_grnd_acx,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"f_grnd_acy",  dom%tpo%now%f_grnd_acy,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"f_ice",       dom%tpo%now%f_ice,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)  
        call nc_read_interp(filename,"dist_margin", dom%tpo%now%dist_margin,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"dist_grline", dom%tpo%now%dist_grline,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"mask_bed",    dom%tpo%now%mask_bed,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"mask_grz",    dom%tpo%now%mask_grz,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)

        call nc_read_interp(filename,"dHdt_n",      dom%tpo%now%dHdt_n,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"H_ice_n",     dom%tpo%now%H_ice_n,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"H_ice_pred",  dom%tpo%now%H_ice_pred,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"H_ice_corr",  dom%tpo%now%H_ice_corr,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        
        call nc_read_interp(filename,"z_srf_n",     dom%tpo%now%z_srf_n,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        
        ! == ybound variables ===

        call nc_read_interp(filename,"z_bed",       dom%bnd%z_bed,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"z_bed_sd",    dom%bnd%z_bed_sd,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"z_sl",        dom%bnd%z_sl,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"H_sed",       dom%bnd%H_sed,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"ice_allowed", dom%bnd%ice_allowed,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"calv_mask",   dom%bnd%calv_mask,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"H_ice_ref",   dom%bnd%H_ice_ref,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"z_bed_ref",   dom%bnd%z_bed_ref,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 

        call nc_read_interp(filename,"smb",         dom%bnd%smb,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"T_srf",       dom%bnd%T_srf,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"bmb_shlf",    dom%bnd%bmb_shlf,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"T_shlf",      dom%bnd%T_shlf,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"Q_geo",       dom%bnd%Q_geo,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 

        call nc_read_interp(filename,"enh_srf",     dom%bnd%enh_srf,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 

        call nc_read_interp(filename,"basins",      dom%bnd%basins,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"basin_mask",  dom%bnd%basin_mask,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"regions",     dom%bnd%regions,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"region_mask", dom%bnd%region_mask,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        
        ! Close the netcdf file
        call nc_close(ncid)

        dom%tpo%par%time = time
        dom%dyn%par%time = time

        ! Write summary 
        write(*,*) 
        write(*,*) "time = ", time, " : loaded restart file: ", trim(filename)
        write(*,*) 
        
        return 

    end subroutine yelmo_restart_read_topo_bnd_internal


    subroutine yelmo_restart_read_internal(dom,filename,time,mps)
        ! Load yelmo variables from restart file: [dyn,therm,mat] 
        ! [tpo] variables loaded using yelmo_restart_read_topo_bnd

        implicit none 

        type(yelmo_class), intent(INOUT) :: dom 
        character(len=*),  intent(IN)    :: filename 
        real(wp),          intent(IN)    :: time 
        type(map_scrip_class), optional, intent(IN) :: mps

        ! Local variables
        integer :: ncid, n, nx, ny, nz, nz_ac, nz_r, n_iso 
        
        ! Read all yelmo data from file,
        ! in order to restart a simulation.
        
        ! Define dimensions of variables 
        nx    = size(dom%grd%xc,1)
        ny    = size(dom%grd%yc,1)
        nz    = size(dom%par%zeta_aa,1) 
        nz_ac = size(dom%par%zeta_ac,1) 
        
        nz_r  = size(dom%thrm%now%enth_rock,3)
        n_iso = size(dom%mat%now%depth_iso,3) 

        ! Assume that first time dimension value is to be read in
        n = 1 

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.FALSE.)
        
        ! == ydyn variables ===

        call nc_read_interp(filename,"ux",            dom%dyn%now%ux,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 
        call nc_read_interp(filename,"uy",            dom%dyn%now%uy,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 
        call nc_read_interp(filename,"uxy",           dom%dyn%now%uxy,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 
        call nc_read_interp(filename,"uz",            dom%dyn%now%uz,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz_ac,1],mps=mps) 
      
        call nc_read_interp(filename,"ux_bar",        dom%dyn%now%ux_bar,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"uy_bar",        dom%dyn%now%uy_bar,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"uxy_bar",       dom%dyn%now%uxy_bar,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 

        call nc_read_interp(filename,"ux_bar_prev",   dom%dyn%now%ux_bar_prev,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"uy_bar_prev",   dom%dyn%now%uy_bar_prev,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        
        call nc_read_interp(filename,"ux_b",          dom%dyn%now%ux_b,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"uy_b",          dom%dyn%now%uy_b,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"uxy_b",         dom%dyn%now%uxy_b,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 

        call nc_read_interp(filename,"ux_s",          dom%dyn%now%ux_s,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"uy_s",          dom%dyn%now%uy_s,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"uxy_s",         dom%dyn%now%uxy_s,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 

        call nc_read_interp(filename,"ux_i",          dom%dyn%now%ux_i,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 
        call nc_read_interp(filename,"uy_i",          dom%dyn%now%uy_i,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 
        call nc_read_interp(filename,"ux_i_bar",      dom%dyn%now%ux_i_bar,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"uy_i_bar",      dom%dyn%now%uy_i_bar,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"uxy_i_bar",     dom%dyn%now%uxy_i_bar,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 

        call nc_read_interp(filename,"duxydt",        dom%dyn%now%duxydt,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 

        call nc_read_interp(filename,"duxdz",         dom%dyn%now%duxdz,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 
        call nc_read_interp(filename,"duydz",         dom%dyn%now%duydz,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 
        call nc_read_interp(filename,"duxdz_bar",     dom%dyn%now%duxdz_bar,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"duydz_bar",     dom%dyn%now%duydz_bar,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 

        call nc_read_interp(filename,"taud_acx",      dom%dyn%now%taud_acx,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"taud_acy",      dom%dyn%now%taud_acy,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"taud",          dom%dyn%now%taud,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 

        call nc_read_interp(filename,"taub_acx",      dom%dyn%now%taub_acx,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"taub_acy",      dom%dyn%now%taub_acy,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"taub",          dom%dyn%now%taub,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 

        call nc_read_interp(filename,"qq_gl_acx",     dom%dyn%now%qq_gl_acx,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"qq_gl_acy",     dom%dyn%now%qq_gl_acy,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        
        call nc_read_interp(filename,"qq_acx",        dom%dyn%now%qq_acx,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"qq_acy",        dom%dyn%now%qq_acy,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"qq",            dom%dyn%now%qq,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 

        call nc_read_interp(filename,"visc_eff",      dom%dyn%now%visc_eff,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 
        call nc_read_interp(filename,"visc_eff_int",  dom%dyn%now%visc_eff_int,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 

        call nc_read_interp(filename,"N_eff",         dom%dyn%now%N_eff,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)        
        call nc_read_interp(filename,"cb_ref",        dom%dyn%now%cb_ref,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"c_bed",         dom%dyn%now%c_bed,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"beta_acx",      dom%dyn%now%beta_acx,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"beta_acy",      dom%dyn%now%beta_acy,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"beta",          dom%dyn%now%beta,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 

        call nc_read_interp(filename,"beta_eff",      dom%dyn%now%beta_eff,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 

        call nc_read_interp(filename,"f_vbvs",        dom%dyn%now%f_vbvs,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 

        call nc_read_interp(filename,"ssa_mask_acx",  dom%dyn%now%ssa_mask_acx,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"ssa_mask_acy",  dom%dyn%now%ssa_mask_acy,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"ssa_err_acx",   dom%dyn%now%ssa_err_acx,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"ssa_err_acy",   dom%dyn%now%ssa_err_acy,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 

        ! == ymat variables ===

        call nc_read_interp(filename,"strn2D_dxx", dom%mat%now%strn2D%dxx,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"strn2D_dyy", dom%mat%now%strn2D%dyy,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"strn2D_dxy", dom%mat%now%strn2D%dxy,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"strn2D_dxz", dom%mat%now%strn2D%dxz,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"strn2D_dyz", dom%mat%now%strn2D%dyz,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"strn2D_de",  dom%mat%now%strn2D%de, ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"strn2D_div", dom%mat%now%strn2D%div,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"strn2D_f_shear",dom%mat%now%strn2D%f_shear,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        
        call nc_read_interp(filename,"strn_dxx",     dom%mat%now%strn%dxx,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 
        call nc_read_interp(filename,"strn_dyy",     dom%mat%now%strn%dyy,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 
        call nc_read_interp(filename,"strn_dxy",     dom%mat%now%strn%dxy,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 
        call nc_read_interp(filename,"strn_dxz",     dom%mat%now%strn%dxz,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 
        call nc_read_interp(filename,"strn_dyz",     dom%mat%now%strn%dyz,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 
        call nc_read_interp(filename,"strn_de",      dom%mat%now%strn%de,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 
        call nc_read_interp(filename,"strn_div",     dom%mat%now%strn%div,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 
        call nc_read_interp(filename,"strn_f_shear", dom%mat%now%strn%f_shear,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 

        call nc_read_interp(filename,"strs2D_txx", dom%mat%now%strs2D%txx,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"strs2D_tyy", dom%mat%now%strs2D%tyy,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"strs2D_txy", dom%mat%now%strs2D%txy,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"strs2D_txz", dom%mat%now%strs2D%txz,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"strs2D_tyz", dom%mat%now%strs2D%tyz,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"strs2D_te",  dom%mat%now%strs2D%te,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"strs2D_tau_eig_1",dom%mat%now%strs2D%tau_eig_1,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"strs2D_tau_eig_2",dom%mat%now%strs2D%tau_eig_2,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        
        call nc_read_interp(filename,"strs_txx", dom%mat%now%strs%txx,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 
        call nc_read_interp(filename,"strs_tyy", dom%mat%now%strs%tyy,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 
        call nc_read_interp(filename,"strs_txy", dom%mat%now%strs%txy,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 
        call nc_read_interp(filename,"strs_txz", dom%mat%now%strs%txz,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 
        call nc_read_interp(filename,"strs_tyz", dom%mat%now%strs%tyz,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 
        call nc_read_interp(filename,"strs_te",  dom%mat%now%strs%te,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 

        call nc_read_interp(filename,"enh",         dom%mat%now%enh,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 
        call nc_read_interp(filename,"enh_bnd",     dom%mat%now%enh_bnd,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"enh_bar",     dom%mat%now%enh_bar,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"ATT",         dom%mat%now%ATT,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 
        call nc_read_interp(filename,"ATT_bar",     dom%mat%now%ATT_bar,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"visc",        dom%mat%now%visc,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 
        call nc_read_interp(filename,"visc_int",    dom%mat%now%visc_int,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 

        call nc_read_interp(filename,"f_shear_bar", dom%mat%now%f_shear_bar,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 

        call nc_read_interp(filename,"dep_time",    dom%mat%now%dep_time,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 
        call nc_read_interp(filename,"depth_iso",   dom%mat%now%depth_iso,ncid=ncid,start=[1,1,1,n],count=[nx,ny,n_iso,1],mps=mps) 
        
        ! == ytherm variables ===

        call nc_read_interp(filename,"enth",        dom%thrm%now%enth,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps)   
        call nc_read_interp(filename,"T_ice",       dom%thrm%now%T_ice,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps)    
        call nc_read_interp(filename,"omega",       dom%thrm%now%omega,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 
        call nc_read_interp(filename,"T_pmp",       dom%thrm%now%T_pmp,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 
        
        call nc_read_interp(filename,"f_pmp",       dom%thrm%now%f_pmp,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"bmb_grnd",    dom%thrm%now%bmb_grnd,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)    
        call nc_read_interp(filename,"Q_strn",      dom%thrm%now%Q_strn,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps)      
        call nc_read_interp(filename,"Q_b",         dom%thrm%now%Q_b,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)         
        call nc_read_interp(filename,"Q_ice_b",     dom%thrm%now%Q_ice_b,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)         
        call nc_read_interp(filename,"T_prime_b",   dom%thrm%now%T_prime_b,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)  
        call nc_read_interp(filename,"H_w",         dom%thrm%now%H_w,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"dHwdt",       dom%thrm%now%dHwdt,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        
        call nc_read_interp(filename,"cp",          dom%thrm%now%cp,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 
        call nc_read_interp(filename,"kt",          dom%thrm%now%kt,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps)      
        call nc_read_interp(filename,"H_cts",       dom%thrm%now%H_cts,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)       
        
        call nc_read_interp(filename,"advecxy",     dom%thrm%now%advecxy,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps)   
        
        call nc_read_interp(filename,"Q_rock",      dom%thrm%now%Q_rock,     ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)        
        call nc_read_interp(filename,"enth_rock",   dom%thrm%now%enth_rock,  ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz_r,1],mps=mps)      
        call nc_read_interp(filename,"T_rock",      dom%thrm%now%T_rock,     ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz_r,1],mps=mps)      

        ! Close the netcdf file
        call nc_close(ncid)

        ! Write summary 

        dom%thrm%par%time = time
        dom%mat%par%time  = time 

        write(*,*) 
        write(*,*) "time = ", time, " : loaded restart file: ", trim(filename)
        write(*,*) 
        
        return 

    end subroutine yelmo_restart_read_internal
    
    subroutine nc_read_interp_wp_2D(filename,vnm,var2D,var2D_in,ncid,start,count,mps,method)
        ! Read in a 2D field and interpolate it using scrip map
        ! as needed. 

        implicit none

        character(len=*),   intent(IN)  :: filename 
        character(len=*),   intent(IN)  :: vnm 
        real(wp),           intent(OUT) :: var2D(:,:) 
        real(wp),optional,  intent(IN)  :: var2D_in(:,:) 
        integer, optional,  intent(IN)  :: ncid 
        integer, optional,  intent(IN)  :: start(:) 
        integer, optional,  intent(IN)  :: count(:) 
        type(map_scrip_class), optional, intent(IN) :: mps
        character(len=*),      optional, intent(IN) :: method 

        ! Local variables 
        integer :: nx, ny 
        integer,  allocatable :: dims(:) 
        real(wp), allocatable :: var2D_src(:,:) 
        character(len=56) :: mapping_method 

        if (present(var2D_in)) then 
            ! Get source array from argument (useful for handling 3D arrays with this routine)

            nx = size(var2D_in,1)
            ny = size(var2D_in,2) 

            allocate(var2D_src(nx,ny))

            var2D_src = var2D_in 

        else 
            ! Load 2D array from file 

            ! Determine dimensions of current variable in the source file
            call nc_dims(filename,vnm,dims=dims)

            nx = dims(1)
            ny = dims(2) 

            allocate(var2D_src(nx,ny))

            ! Load the variable from the file to the local 2D array
            call nc_read(filename,vnm,var2D_src,ncid=ncid,start=start,count=[nx,ny,1],missing_value=MV)
        
        end if 


        if (nx .eq. size(var2D,1) .and. ny .eq. size(var2D,2)) then 
            ! Assume no interpolation needed, copy variable for output directly

            var2D = var2D_src 

        else 
            ! Map local source array it to our target array 

            ! Determine mapping method for this variable 
            mapping_method = "mean"
            if (present(method)) mapping_method = trim(method) 

            ! Safety check 
            if (.not. present(mps)) then 
                write(io_unit_err,*) ""
                write(io_unit_err,*) "nc_read_interp:: Error: map_scrip_class object must &
                        &be provided as an argument since array read from file does not &
                        &match the Yelmo array size."
                write(io_unit_err,*) "filename: ", trim(filename)
                write(io_unit_err,*) "variable: ", trim(vnm)
                write(io_unit_err,*) "dims in file:         ", nx, ny 
                write(io_unit_err,*) "dims in yelmo object: ", size(var2D,1), size(var2D,2)
                stop 
            end if 

            ! Perform conservative interpolation 
            var2D = MV 
            call map_scrip_field(mps,vnm,var2D_src,var2D,method=mapping_method, &
                                        missing_value=MV,fill_method="nn")

        end if 

        return

    end subroutine nc_read_interp_wp_2D

    subroutine nc_read_interp_wp_3D(filename,vnm,var3D,ncid,start,count,mps,method)
        ! Read in a 2D field and interpolate it using scrip map
        ! as needed. 

        implicit none

        character(len=*),   intent(IN)  :: filename 
        character(len=*),   intent(IN)  :: vnm 
        real(wp),           intent(OUT) :: var3D(:,:,:) 
        integer, optional,  intent(IN)  :: ncid 
        integer, optional,  intent(IN)  :: start(:) 
        integer, optional,  intent(IN)  :: count(:) 
        type(map_scrip_class), optional, intent(IN) :: mps
        character(len=*),      optional, intent(IN) :: method 

        ! Local variables 
        integer :: nx, ny, nz, k  
        integer,  allocatable :: dims(:) 
        real(wp), allocatable :: var3D_src(:,:,:) 

        ! Determine dimensions of current variable in the source file
        call nc_dims(filename,vnm,dims=dims)

        nx = dims(1)
        ny = dims(2) 
        nz = dims(3) 

        allocate(var3D_src(nx,ny,nz))

        ! Safety check 
        if (nz .ne. size(var3D,3)) then 

            write(io_unit_err,*) ""
            write(io_unit_err,*) "nc_read_interp_wp_3D:: Error: vertical dimension of variable in &
                    &input file does not match vertical dimension of yelmo object. Vertical &
                    & interpolation is not yet supported."
            write(io_unit_err,*) "filename  = ", trim(filename)
            write(io_unit_err,*) "variable  = ", trim(vnm)
            write(io_unit_err,*) "nz[file]  = ", nz
            write(io_unit_err,*) "nz[yelmo] = ", size(var3D,3) 
            stop 
        end if 
            
        ! Read in full 3D variable of interest 
        call nc_read(filename,vnm,var3D_src,ncid=ncid,start=start,count=[nx,ny,nz,1],missing_value=MV)
        
        ! Loop over vertical dimension and apply interpolation 
        do k = 1, nz 
            
            call nc_read_interp_wp_2D(filename,vnm,var3D(:,:,k),var3D_src(:,:,k),ncid, &
                                                                start,count,mps,method)

        end do 

        return

    end subroutine nc_read_interp_wp_3D

    subroutine nc_read_interp_int_2D(filename,vnm,var2D,var2D_in,ncid,start,count,mps,method)
        ! Read in a 2D field and interpolate it using scrip map
        ! as needed. 

        implicit none

        character(len=*),   intent(IN)  :: filename 
        character(len=*),   intent(IN)  :: vnm 
        integer,            intent(OUT) :: var2D(:,:) 
        integer, optional,  intent(IN)  :: var2D_in(:,:) 
        integer, optional,  intent(IN)  :: ncid 
        integer, optional,  intent(IN)  :: start(:) 
        integer, optional,  intent(IN)  :: count(:) 
        type(map_scrip_class), optional, intent(IN) :: mps
        character(len=*),      optional, intent(IN) :: method 

        ! Local variables 
        integer :: nx, ny 
        integer,  allocatable :: dims(:) 
        integer,  allocatable :: var2D_src(:,:) 
        character(len=56) :: mapping_method 

        if (present(var2D_in)) then 
            ! Get source array from argument (useful for handling 3D arrays with this routine)

            nx = size(var2D_in,1)
            ny = size(var2D_in,2) 

            allocate(var2D_src(nx,ny))

            var2D_src = var2D_in 

        else 
            ! Load 2D array from file 

            ! Determine dimensions of current variable in the source file
            call nc_dims(filename,vnm,dims=dims)

            nx = dims(1)
            ny = dims(2) 

            allocate(var2D_src(nx,ny))

            ! Load the variable from the file to the local 2D array
            call nc_read(filename,vnm,var2D_src,ncid=ncid,start=start,count=[nx,ny,1],missing_value=int(MV))
        
        end if 


        if (nx .eq. size(var2D,1) .and. ny .eq. size(var2D,2)) then 
            ! Assume no interpolation needed, copy variable for output directly

            var2D = var2D_src 

        else 
            ! Map local source array it to our target array 

            ! Determine mapping method for this variable 
            mapping_method = "mean"
            if (present(method)) mapping_method = trim(method) 

            ! Perform conservative interpolation 
            var2D = MV 
            call map_scrip_field(mps,vnm,var2D_src,var2D,method=mapping_method, &
                                        missing_value=int(MV),fill_method="nn")

        end if 

        return

    end subroutine nc_read_interp_int_2D

    subroutine nc_read_interp_int_3D(filename,vnm,var3D,ncid,start,count,mps,method)
        ! Read in a 2D field and interpolate it using scrip map
        ! as needed. 

        implicit none

        character(len=*),   intent(IN)  :: filename 
        character(len=*),   intent(IN)  :: vnm 
        integer,            intent(OUT) :: var3D(:,:,:) 
        integer, optional,  intent(IN)  :: ncid 
        integer, optional,  intent(IN)  :: start(:) 
        integer, optional,  intent(IN)  :: count(:) 
        type(map_scrip_class), optional, intent(IN) :: mps
        character(len=*),      optional, intent(IN) :: method 

        ! Local variables 
        integer :: nx, ny, nz, k  
        integer,  allocatable :: dims(:) 
        integer,  allocatable :: var3D_src(:,:,:) 

        ! Determine dimensions of current variable in the source file
        call nc_dims(filename,vnm,dims=dims)

        nx = dims(1)
        ny = dims(2) 
        nz = dims(3) 

        allocate(var3D_src(nx,ny,nz))

        ! Safety check 
        if (nz .ne. size(var3D,3)) then 

            write(io_unit_err,*) ""
            write(io_unit_err,*) "nc_read_interp_int_3D:: Error: vertical dimension of variable in &
                    &input file does not match vertical dimension of yelmo object. Vertical &
                    & interpolation is not yet supported."
            write(io_unit_err,*) "filename  = ", trim(filename)
            write(io_unit_err,*) "variable  = ", trim(vnm)
            write(io_unit_err,*) "nz[file]  = ", nz
            write(io_unit_err,*) "nz[yelmo] = ", size(var3D,3) 
            stop 
        end if 
            
        ! Read in full 3D variable of interest 
        call nc_read(filename,vnm,var3D_src,ncid=ncid,start=start,count=[nx,ny,nz,1],missing_value=int(MV))
        
        ! Loop over vertical dimension and apply interpolation 
        do k = 1, nz 
            
            call nc_read_interp_int_2D(filename,vnm,var3D(:,:,k),var3D_src(:,:,k),ncid, &
                                                                    start,count,mps,method)

        end do 

        return

    end subroutine nc_read_interp_int_3D

    subroutine nc_read_interp_logical_2D(filename,vnm,var2D,var2D_in,ncid,start,count,mps,method)
        ! Read in a 2D field and interpolate it using scrip map
        ! as needed. 

        implicit none

        character(len=*),   intent(IN)  :: filename 
        character(len=*),   intent(IN)  :: vnm 
        logical,            intent(OUT) :: var2D(:,:) 
        logical, optional,  intent(IN)  :: var2D_in(:,:) 
        integer, optional,  intent(IN)  :: ncid 
        integer, optional,  intent(IN)  :: start(:) 
        integer, optional,  intent(IN)  :: count(:) 
        type(map_scrip_class), optional, intent(IN) :: mps
        character(len=*),      optional, intent(IN) :: method 

        ! Local variables 
        integer :: nx, ny 
        integer,  allocatable :: dims(:) 
        integer,  allocatable :: var2D_src(:,:) 
        integer,  allocatable :: var2D_int(:,:) 
        character(len=56) :: mapping_method 

        if (present(var2D_in)) then 
            ! Get source array from argument (useful for handling 3D arrays with this routine)

            nx = size(var2D_in,1)
            ny = size(var2D_in,2) 

            allocate(var2D_src(nx,ny))

            var2D_src = 0 
            where (var2D_in) var2D_src = 1 

        else 
            ! Load 2D array from file 

            ! Determine dimensions of current variable in the source file
            call nc_dims(filename,vnm,dims=dims)

            nx = dims(1)
            ny = dims(2) 

            allocate(var2D_src(nx,ny))

            ! Load the variable from the file to the local 2D array
            call nc_read(filename,vnm,var2D_src,ncid=ncid,start=start,count=[nx,ny,1],missing_value=int(MV))
        
        end if 


        if (nx .eq. size(var2D,1) .and. ny .eq. size(var2D,2)) then 
            ! Assume no interpolation needed, copy variable for output directly

            var2D = .FALSE. 
            where(var2D_src .eq. 1) var2D = .TRUE.  

        else 
            ! Map local source array it to our target array 

            ! Determine mapping method for this variable 
            mapping_method = "count"
            if (present(method)) mapping_method = trim(method) 

            ! Perform conservative interpolation 
            allocate(var2D_int(size(var2D,1),size(var2D,2)))
            var2D_int = int(MV)
            call map_scrip_field(mps,vnm,var2D_src,var2D_int,method=mapping_method, &
                                        missing_value=int(MV),fill_method="nn")

            var2D = .FALSE. 
            where(var2D_int .eq. 1) var2D = .TRUE.  

        end if 

        return

    end subroutine nc_read_interp_logical_2D

end module yelmo_io


