

module yelmo_io
    
    use ncio 
    
    use yelmo_defs 
    use yelmo_tools, only : get_region_indices
    use yelmo_grid 
    
    use variable_io
    use interp2D
    use coordinates_mapping_scrip, only : map_scrip_class, map_scrip_init, map_scrip_field, &
                                            gen_map_filename, nc_read_interp
    
    implicit none

    type yelmo_io_tables
        type(var_io_type) :: v
        type(var_io_type), allocatable :: tpo(:)
        type(var_io_type), allocatable :: dyn(:)
        type(var_io_type), allocatable :: mat(:)
        type(var_io_type), allocatable :: thrm(:)
        type(var_io_type), allocatable :: bnd(:)
        type(var_io_type), allocatable :: dta(:)
    end type
    
    private 
    public :: yelmo_write_init
    public :: yelmo_write_step
    public :: yelmo_write_var
    public :: yelmo_write_step_model_metrics
    public :: yelmo_write_step_pd_metrics
    public :: yelmo_restart_write
    public :: yelmo_restart_read_topo_bnd
    public :: yelmo_restart_read

contains

    subroutine yelmo_write_init(ylmo,filename,time_init,units,irange,jrange)
        ! Initialize a NetCDF file for Yelmo output.
        ! Produces a file with all possible dimension information 
        ! to be able to plot different variables of the user's choice. 
        ! Also, irange=[i1,i2] and jrange=[j1,j2] can be used to limit 
        ! the dimensions to a specific horizontal region of the domain.
        
        implicit none 

        type(yelmo_class), intent(IN) :: ylmo 
        character(len=*),  intent(IN) :: filename
        real(wp),          intent(IN) :: time_init
        character(len=*),  intent(IN) :: units 
        integer,           intent(IN), optional :: irange(2)
        integer,           intent(IN), optional :: jrange(2)
        
        ! Initialize file by writing grid info
        call yelmo_grid_write(ylmo%grd, filename, ylmo%par%domain, ylmo%par%grid_name, create=.TRUE.,irange=irange,jrange=jrange)

        ! Initialize netcdf file and dimensions
        call nc_write_dim(filename,"month",     x=1,dx=1,nx=12,         units="month")
        call nc_write_dim(filename,"zeta",      x=ylmo%par%zeta_aa,     units="1")
        call nc_write_dim(filename,"zeta_ac",   x=ylmo%par%zeta_ac,     units="1")
        call nc_write_dim(filename,"zeta_rock", x=ylmo%thrm%par%zr%zeta_aa,units="1")
        call nc_write_dim(filename,"age_iso",   x=ylmo%mat%par%age_iso, units="kyr")
        call nc_write_dim(filename,"pd_age_iso",x=ylmo%dta%pd%age_iso,  units="kyr")
        call nc_write_dim(filename,"pc_steps",  x=1,dx=1,nx=3,          units="1")
        
        call nc_write_dim(filename,"time",      x=time_init,dx=1.0_wp,nx=1,units=trim(units),unlimited=.TRUE.)

        if (ylmo%grd%is_projection) then 
            call nc_write_attr(filename,"xc","standard_name","projection_x_coordinate")
            call nc_write_attr(filename,"yc","standard_name","projection_y_coordinate")
        end if 

        return

    end subroutine yelmo_write_init

    subroutine yelmo_write_step(ylmo,filename,time,nms_tpo,nms_dyn,nms_mat,nms_therm, &
                                                nms_bnd,nms_dta,compare_pd,irange,jrange)

        implicit none 
        
        type(yelmo_class), intent(IN) :: ylmo        
        character(len=*),  intent(IN) :: filename
        real(prec),        intent(IN) :: time
        character(len=*),  intent(IN), optional :: nms_tpo(:)
        character(len=*),  intent(IN), optional :: nms_dyn(:)
        character(len=*),  intent(IN), optional :: nms_mat(:)
        character(len=*),  intent(IN), optional :: nms_therm(:)
        character(len=*),  intent(IN), optional :: nms_bnd(:)
        character(len=*),  intent(IN), optional :: nms_dta(:)
        logical,           intent(IN), optional :: compare_pd
        integer,           intent(IN), optional :: irange(2)
        integer,           intent(IN), optional :: jrange(2)
        
        ! Local variables
        integer    :: ncid, n, q, qtot
        character(len=56), allocatable :: names(:) 
        logical ::  write_pd_metrics 
        integer :: i1, i2, j1, j2

        type(yelmo_io_tables) :: io

        ! Get indices for current domain of interest
        call get_region_indices(i1,i2,j1,j2,ylmo%grd%nx,ylmo%grd%ny,irange,jrange)


        if (present(nms_tpo)) then 
            qtot = size(nms_tpo,1)
            allocate(names(qtot))
            do q = 1, qtot 
                names(q) = trim(nms_tpo(q))
            end do 
        else 
            qtot = 18 
            allocate(names(qtot))
            names(1)  = "H_ice"
            names(2)  = "z_srf"
            names(3)  = "z_bed"
            names(4)  = "mask_bed"
            names(5)  = "uxy_b"
            names(6)  = "uxy_s"
            names(7)  = "uxy_bar"
            names(8)  = "ux_bar"
            names(9)  = "uy_bar"
            names(10) = "beta"
            names(11) = "visc_bar"
            names(12) = "T_prime_b"
            names(13) = "H_w"
            names(14) = "mb_net"
            names(15) = "smb"
            names(16) = "bmb"
            names(17) = "cmb"
            names(18) = "z_sl"

        end if 

        write_pd_metrics = .TRUE. 
        if (present(compare_pd)) write_pd_metrics = compare_pd

        ! Load io tables
        call load_var_io_table(io%tpo,"input/yelmo-variables-ytopo.md")
        call load_var_io_table(io%dyn,"input/yelmo-variables-ydyn.md")
        call load_var_io_table(io%mat,"input/yelmo-variables-ymat.md")
        call load_var_io_table(io%thrm,"input/yelmo-variables-ytherm.md")
        call load_var_io_table(io%bnd,"input/yelmo-variables-ybound.md")
        call load_var_io_table(io%dta,"input/yelmo-variables-ydata.md")

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_time_index(filename,"time",time,ncid)

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        if (n .eq. 1) then
            ! First time index, write some static fields

            call nc_write(filename,"basins",   ylmo%bnd%basins(i1:i2,j1:j2),  dim1="xc",dim2="yc",units="(0 - 8)",long_name="Hydrological basins")
            call nc_write(filename,"regions",  ylmo%bnd%regions(i1:i2,j1:j2), dim1="xc",dim2="yc",units="(0 - 8)",long_name="Domain regions") 
            call nc_write(filename,"z_bed_sd", ylmo%bnd%z_bed_sd(i1:i2,j1:j2),dim1="xc",dim2="yc",units="m",long_name="Stdev(z_bed)")
        
        end if

        ! Write model metrics (model speed, dt, eta)
        call yelmo_write_step_model_metrics(filename,ylmo,n,ncid)
 
        if (write_pd_metrics) then 
            ! Write present-day data metrics (rmse[H],etc)
            call yelmo_write_step_pd_metrics(filename,ylmo,n,ncid)
        end if  
        
        ! Determine number of variables to write 
        qtot = size(names) 

        ! Loop over variables and write each variable
if (.FALSE.) then
        ! Use variable_io

        if (present(nms_tpo)) then
            do q = 1, size(nms_tpo) 
                call find_var_io_in_table(io%v,nms_tpo(q),io%tpo)
                call yelmo_write_var_io_ytopo(filename,io%v,ylmo,n,ncid)
            end do
        end if

else
        ! Use hard-coded variable writing
        do q = 1, qtot 
               call yelmo_write_var(filename,names(q),ylmo,n,ncid)
        end do 
end if

        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine yelmo_write_step

    subroutine yelmo_write_var_io_ytopo(filename,v,ylmo,n,ncid,irange,jrange)

        implicit none

        character(len=*),  intent(IN) :: filename 
        type(var_io_type), intent(IN) :: v
        type(yelmo_class), intent(IN) :: ylmo 
        integer,           intent(IN) :: n 
        integer,           intent(IN), optional :: ncid 
        integer,           intent(IN), optional :: irange(2)
        integer,           intent(IN), optional :: jrange(2)
        
        ! Local variables
        integer :: i1, i2, j1, j2
        character(len=32), allocatable :: dims(:)

        ! Get indices for current domain of interest
        call get_region_indices(i1,i2,j1,j2,ylmo%grd%nx,ylmo%grd%ny,irange,jrange)

        ! Allocate local representation of dims to be able to add "time" as last dimension
        allocate(dims(v%ndims+1))
        dims(1:v%ndims) = v%dims
        dims(v%ndims+1) = "time"

        select case(trim(v%varname))

            case("H_ice")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%H_ice(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("dHidt")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%dHidt(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("dHidt_dyn")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%dHidt_dyn(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("mb_net")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%mb_net(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("mb_relax")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%mb_relax(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("mb_resid")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%mb_resid(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("mb_err")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%mb_err(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("smb")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%smb(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("bmb")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%bmb(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("fmb")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%fmb(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("dmb")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%dmb(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("cmb")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%cmb(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("bmb_ref")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%bmb_ref(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("fmb_ref")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%fmb_ref(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("dmb_ref")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%dmb_ref(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("cmb_flt")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%cmb_flt(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("cmb_grnd")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%cmb_grnd(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("z_srf")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%z_srf(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("dzsdt")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%dzsdt(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("mask_adv")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%mask_adv(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("eps_eff")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%eps_eff(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("tau_eff")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%tau_eff(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("z_base")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%z_base(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("dzsdx")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%dzsdx(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("dzsdy")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%dzsdy(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("dHidx")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%dHidx(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("dHidy")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%dHidy(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("dzbdx")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%dzbdx(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("dzbdy")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%dzbdy(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("H_eff")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%H_eff(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("H_grnd")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%H_grnd(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("H_calv")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%H_calv(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("kt_calv")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%kt(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("z_bed_filt")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%z_bed_filt(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("f_grnd")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%f_grnd(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("f_grnd_acx")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%f_grnd_acx(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("f_grnd_acy")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%f_grnd_acy(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("f_grnd_ab")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%f_grnd_ab(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("f_ice")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%f_ice(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("f_grnd_bmb")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%f_grnd_bmb(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("f_grnd_pin")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%f_grnd_pin(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("dist_margin")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%dist_margin(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("dist_grline")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%dist_grline(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("mask_bed")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%mask_bed(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("mask_grz")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%mask_grz(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("mask_frnt")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%mask_frnt(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("dHidt_dyn_n")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%dHidt_dyn_n(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("H_ice_n")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%H_ice_n(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("z_srf_n")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%z_srf_n(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("H_ice_dyn")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%H_ice_dyn(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("f_ice_dyn")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%f_ice_dyn(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            
            case("pc_pred_H_ice")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%pred%H_ice(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("pc_pred_dHidt_dyn")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%pred%dHidt_dyn(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("pc_pred_mb_net")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%pred%mb_net(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("pc_pred_smb")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%pred%smb(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("pc_pred_bmb")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%pred%bmb(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("pc_pred_fmb")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%pred%fmb(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("pc_pred_dmb")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%pred%dmb(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("pc_pred_cmb")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%pred%cmb(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("pc_corr_H_ice")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%corr%H_ice(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("pc_corr_dHidt_dyn")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%corr%dHidt_dyn(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("pc_corr_mb_net")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%corr%mb_net(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("pc_corr_smb")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%corr%smb(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("pc_corr_bmb")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%corr%bmb(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("pc_corr_fmb")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%corr%fmb(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("pc_corr_dmb")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%corr%dmb(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("pc_corr_cmb")
                call nc_write(filename,trim(v%varname),ylmo%tpo%now%corr%cmb(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            
            case DEFAULT 

                write(io_unit_err,*) 
                write(io_unit_err,*) "yelmo_write_var_io_ytopo:: Error: variable not yet supported."
                write(io_unit_err,*) "variable = ", trim(v%varname)
                write(io_unit_err,*) "filename = ", trim(filename)
                stop 
                
        end select

        return

    end subroutine yelmo_write_var_io_ytopo

    subroutine yelmo_write_var_io_ydyn(filename,v,ylmo,n,ncid,irange,jrange)

        implicit none

        character(len=*),  intent(IN) :: filename 
        type(var_io_type), intent(IN) :: v
        type(yelmo_class), intent(IN) :: ylmo 
        integer,           intent(IN) :: n 
        integer,           intent(IN), optional :: ncid 
        integer,           intent(IN), optional :: irange(2)
        integer,           intent(IN), optional :: jrange(2)
        
        ! Local variables
        integer :: i1, i2, j1, j2
        character(len=32), allocatable :: dims(:)

        ! Get indices for current domain of interest
        call get_region_indices(i1,i2,j1,j2,ylmo%grd%nx,ylmo%grd%ny,irange,jrange)

        ! Allocate local representation of dims to be able to add "time" as last dimension
        allocate(dims(v%ndims+1))
        dims(1:v%ndims) = v%dims
        dims(v%ndims+1) = "time"
        
        select case(trim(v%varname))
            
            case("ux") ! 3D
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%ux(i1:i2,j1:j2,:), &
                            start=[1,1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("uy") ! 3D
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%uy(i1:i2,j1:j2,:), &
                            start=[1,1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("uxy") ! 3D
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%uxy(i1:i2,j1:j2,:), &
                            start=[1,1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("uz") ! 3D
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%uz(i1:i2,j1:j2,:), &
                            start=[1,1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("uz_star") ! 3D
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%uz_star(i1:i2,j1:j2,:), &
                            start=[1,1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("ux_bar")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%ux_bar(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("uy_bar")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%uy_bar(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("uxy_bar")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%uxy_bar(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("ux_bar_prev")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%ux_bar_prev(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("uy_bar_prev")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%uy_bar_prev(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("ux_b")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%ux_b(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("uy_b")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%uy_b(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("uz_b")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%uz_b(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("uxy_b")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%uxy_b(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("ux_s")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%ux_s(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("uy_s")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%uy_s(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("uz_s")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%uz_s(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("uxy_s")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%uxy_s(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("ux_i") ! 3D
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%ux_i(i1:i2,j1:j2,:), &
                            start=[1,1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("uy_i") ! 3D
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%uy_i(i1:i2,j1:j2,:), &
                            start=[1,1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("ux_i_bar")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%ux_i_bar(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("uy_i_bar")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%uy_i_bar(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("uxy_i_bar")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%uxy_i_bar(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("duxydt")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%duxydt(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("duxdz") ! 3D
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%duxdz(i1:i2,j1:j2,:), &
                            start=[1,1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("duydz") ! 3D
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%duydz(i1:i2,j1:j2,:), &
                            start=[1,1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("duxdz_bar")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%duxdz_bar(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("duydz_bar")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%duydz_bar(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("taud_acx")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%taud_acx(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("taud_acy")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%taud_acy(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("taud")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%taud(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("taub_acx")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%taub_acx(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("taub_acy")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%taub_acy(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("taub")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%taub(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("taul_int_acx")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%taul_int_acx(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("taul_int_acy")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%taul_int_acy(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("qq_gl_acx")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%qq_gl_acx(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("qq_gl_acy")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%qq_gl_acy(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("qq_acx")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%qq_acx(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("qq_acy")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%qq_acy(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("qq")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%qq(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("de_eff") ! 3D
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%de_eff(i1:i2,j1:j2,:), &
                            start=[1,1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("visc_eff") ! 3D
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%visc_eff(i1:i2,j1:j2,:), &
                            start=[1,1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("visc_eff_int")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%visc_eff_int(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("N_eff")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%N_eff(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("cb_tgt")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%cb_tgt(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("cb_ref")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%cb_ref(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("c_bed")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%c_bed(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("beta_acx")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%beta_acx(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("beta_acy")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%beta_acy(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("beta")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%beta(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("beta_eff")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%beta_eff(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("f_vbvs")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%f_vbvs(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("ssa_mask_acx")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%ssa_mask_acx(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("ssa_mask_acy")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%ssa_mask_acy(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("ssa_err_acx")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%ssa_err_acx(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("ssa_err_acy")
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%ssa_err_acy(i1:i2,j1:j2), &
                            start=[1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("jvel_dxx") ! 3D
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%jvel%dxx(i1:i2,j1:j2,:), &
                            start=[1,1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("jvel_dxy") ! 3D
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%jvel%dxy(i1:i2,j1:j2,:), &
                            start=[1,1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("jvel_dxz") ! 3D
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%jvel%dxz(i1:i2,j1:j2,:), &
                            start=[1,1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("jvel_dyx") ! 3D
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%jvel%dyx(i1:i2,j1:j2,:), &
                            start=[1,1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("jvel_dyy") ! 3D
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%jvel%dyy(i1:i2,j1:j2,:), &
                            start=[1,1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("jvel_dyz") ! 3D
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%jvel%dyz(i1:i2,j1:j2,:), &
                            start=[1,1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("jvel_dzx") ! 3D
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%jvel%dzx(i1:i2,j1:j2,:), &
                            start=[1,1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("jvel_dzy") ! 3D
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%jvel%dzy(i1:i2,j1:j2,:), &
                            start=[1,1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("jvel_dzz") ! 3D
                call nc_write(filename,trim(v%varname),ylmo%dyn%now%jvel%dzz(i1:i2,j1:j2,:), &
                            start=[1,1,1,n],units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            
            case DEFAULT 

                write(io_unit_err,*) 
                write(io_unit_err,*) "yelmo_write_var_io_ydyn:: Error: variable not yet supported."
                write(io_unit_err,*) "variable = ", trim(v%varname)
                write(io_unit_err,*) "filename = ", trim(filename)
                stop 
                
        end select

        return

    end subroutine yelmo_write_var_io_ydyn

    subroutine yelmo_write_var(filename,varname,ylmo,n,ncid)

        implicit none

        character(len=*),  intent(IN) :: filename 
        character(len=*),  intent(IN) :: varname
        type(yelmo_class), intent(IN) :: ylmo 
        integer                       :: n 
        integer, optional             :: ncid 
        
        select case(trim(varname))

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
            case("uxy_bar")
                call nc_write(filename,"uxy_bar",ylmo%dyn%now%uxy_bar,units="m/yr",long_name="Depth-averaged velocity magnitude", &
                                dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            case("ux_bar")
                call nc_write(filename,"ux_bar",ylmo%dyn%now%ux_bar,units="m/yr",long_name="Depth-averaged x-velocity", &
                                dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            case("uy_bar")
                call nc_write(filename,"uy_bar",ylmo%dyn%now%uy_bar,units="m/yr",long_name="Depth-averaged y-velocity", &
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
            case("mb_net")
                call nc_write(filename,"mb_net",ylmo%tpo%now%mb_net,units="m/yr ice equiv.",long_name="Net applied mass balance", &
                                dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            case("smb")
                call nc_write(filename,"smb",ylmo%tpo%now%smb,units="m/yr ice equiv.",long_name="Net surface mass balance", &
                                dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            case("bmb")
                call nc_write(filename,"bmb",ylmo%tpo%now%bmb,units="m/yr ice equiv.",long_name="Basal mass balance", &
                                dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            case("cmb")
                call nc_write(filename,"cmb",ylmo%tpo%now%bmb,units="m/yr ice equiv.",long_name="Calving mass balance", &
                                dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            case("z_sl")
                call nc_write(filename,"z_sl",ylmo%bnd%z_sl,units="m",long_name="Sea level rel. to present", &
                                dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
            
            
            case("ux")
                call nc_write(filename,"ux",ylmo%dyn%now%ux,units="m/yr",long_name="Velocity, x-direction", &
                                dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
            case("uy")
                call nc_write(filename,"uy",ylmo%dyn%now%ux,units="m/yr",long_name="Velocity, y-direction", &
                                dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
            case("uz")
                call nc_write(filename,"uz",ylmo%dyn%now%ux,units="m/yr",long_name="Velocity, z-direction", &
                                dim1="xc",dim2="yc",dim3="zeta_ac",dim4="time",start=[1,1,1,n],ncid=ncid)
                                
            case DEFAULT 

                write(io_unit_err,*) 
                write(io_unit_err,*) "yelmo_write_var:: Error: variable not yet supported."
                write(io_unit_err,*) "variable = ", trim(varname)
                write(io_unit_err,*) "filename = ", trim(filename)
                stop 
                
        end select

        return

    end subroutine yelmo_write_var

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
        call nc_write(filename,"rmse_uxy",ylmo%dta%pd%rmse_uxy,units="m/yr",long_name="RMSE - Surface velocity", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"rmse_uxy_log",ylmo%dta%pd%rmse_loguxy,units="log(m/yr)",long_name="RMSE - Log surface velocity", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"rmse_iso",ylmo%dta%pd%rmse_iso,units="m",long_name="RMSE - isochronal layer depths", &
                  dim1="pd_age_iso",dim2="time",start=[1,n],missing_value=mv,ncid=ncid)

        return 

    end subroutine yelmo_write_step_pd_metrics

    subroutine yelmo_restart_write(dom,filename,time,init,irange,jrange)
        ! Write all yelmo data to file, so that it can be
        ! read later to restart a simulation.
        
        implicit none 

        type(yelmo_class), intent(IN) :: dom
        character(len=*),  intent(IN) :: filename 
        real(wp),          intent(IN) :: time 
        logical,           intent(IN), optional :: init 
        integer,           intent(IN), optional :: irange(2)
        integer,           intent(IN), optional :: jrange(2)
        
        ! Local variables
        integer  :: ncid, n, q, nt
        integer :: i1, i2, j1, j2 
        logical  :: initialize_file  
        
        type(yelmo_io_tables) :: io

        initialize_file = .TRUE. 
        if (present(init)) initialize_file = init

        ! Get indices for current domain of interest
        call get_region_indices(i1,i2,j1,j2,dom%grd%nx,dom%grd%ny,irange,jrange)

        ! Load variable io tables
        call load_var_io_table(io%tpo,"input/yelmo-variables-ytopo.md")
        call load_var_io_table(io%dyn,"input/yelmo-variables-ydyn.md")
        call load_var_io_table(io%mat,"input/yelmo-variables-ymat.md")
        call load_var_io_table(io%thrm,"input/yelmo-variables-ytherm.md")
        call load_var_io_table(io%bnd,"input/yelmo-variables-ybound.md")
        call load_var_io_table(io%dta,"input/yelmo-variables-ydata.md")

        ! == Initialize netcdf file ==============================================

        if (initialize_file) then
            call yelmo_write_init(dom,filename,time,"years",irange,jrange)
        end if 
        
        ! == Begin writing data ==============================================
        
        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)
        
        if (initialize_file) then 

            ! Current time index to write will be the first and only one 
            n = 1 

        else 

            ! Determine current writing time step 
            n = nc_time_index(filename,"time",time,ncid)

            ! Update the time step
            call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        end if 

        ! == Predictor-corrector (pc) variables ===
        ! (these will not be read in by yelmo_restart_read, but can be useful to output for diagnostics)

        call nc_write(filename,"pc_tau",       dom%time%pc_tau,        units="m/yr",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])
        call nc_write(filename,"pc_tau_masked",dom%time%pc_tau_masked, units="m/yr",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])
        call nc_write(filename,"pc_tau_max",   dom%time%pc_tau_max,    units="m/yr",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])
        
        ! == time variables ===

        ! Note: these variables are not defined on the 2D grid, so they 
        ! will give interpolation errors to `cdo remapcon`. They 
        ! are currently not used by a restarted simulation, so they 
        ! are commented out here. However, it could be useful in the future
        ! to have these variables in the restart file, and potentially 
        ! others that do not depend on the grid.
        
        ! ajr: writing these values is reactivated to see if it improves restart file performance

        call nc_write(filename,"pc_dt",        dom%time%pc_dt,         units="yr",  dim1="pc_steps",dim2="time",ncid=ncid,start=[1,n],count=[3,1],grid_mapping="")
        call nc_write(filename,"pc_eta",       dom%time%pc_eta,        units="m/yr",dim1="pc_steps",dim2="time",ncid=ncid,start=[1,n],count=[3,1],grid_mapping="")

        ! == ytopo variables ===
        do q = 1, size(io%tpo)
            call yelmo_write_var_io_ytopo(filename,io%tpo(q),dom,n,ncid,irange,jrange)
        end do

        ! = ytopo_pc variables (just for diagnostic output) ===
        
        ! call nc_write(filename,"pc_pred_H_ice",     dom%tpo%now%pred%H_ice,     units="m",    dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])
        ! call nc_write(filename,"pc_pred_dHidt_dyn", dom%tpo%now%pred%dHidt_dyn, units="m",    dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])
        ! call nc_write(filename,"pc_pred_mb_net",    dom%tpo%now%pred%mb_net,    units="m/yr", dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])
        ! call nc_write(filename,"pc_pred_smb",       dom%tpo%now%pred%smb,       units="m/yr", dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])
        ! call nc_write(filename,"pc_pred_bmb",       dom%tpo%now%pred%bmb,       units="m/yr", dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])
        ! call nc_write(filename,"pc_pred_fmb",       dom%tpo%now%pred%fmb,       units="m/yr", dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])
        ! call nc_write(filename,"pc_pred_dmb",       dom%tpo%now%pred%dmb,       units="m/yr", dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])
        ! call nc_write(filename,"pc_pred_cmb",       dom%tpo%now%pred%cmb,       units="m/yr", dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])
        ! call nc_write(filename,"pc_corr_H_ice",     dom%tpo%now%corr%H_ice,     units="m",    dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])
        ! call nc_write(filename,"pc_corr_dHidt_dyn", dom%tpo%now%corr%dHidt_dyn, units="m",    dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])
        ! call nc_write(filename,"pc_corr_mb_net",    dom%tpo%now%corr%mb_net,    units="m/yr", dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])
        ! call nc_write(filename,"pc_corr_smb",       dom%tpo%now%corr%smb,       units="m/yr", dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])
        ! call nc_write(filename,"pc_corr_bmb",       dom%tpo%now%corr%bmb,       units="m/yr", dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])
        ! call nc_write(filename,"pc_corr_fmb",       dom%tpo%now%corr%fmb,       units="m/yr", dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])
        ! call nc_write(filename,"pc_corr_dmb",       dom%tpo%now%corr%dmb,       units="m/yr", dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])
        ! call nc_write(filename,"pc_corr_cmb",       dom%tpo%now%corr%cmb,       units="m/yr", dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])

        ! == ydyn variables ===
        do q = 1, size(io%dyn)
            call yelmo_write_var_io_ydyn(filename,io%dyn(q),dom,n,ncid,irange,jrange)
        end do

        ! == ymat variables === 

        call nc_write(filename,"strn2D_dxx", dom%mat%now%strn2D%dxx,  units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n]) 
        call nc_write(filename,"strn2D_dyy", dom%mat%now%strn2D%dyy,  units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n]) 
        call nc_write(filename,"strn2D_dxy", dom%mat%now%strn2D%dxy,  units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n]) 
        call nc_write(filename,"strn2D_dxz", dom%mat%now%strn2D%dxz,  units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n]) 
        call nc_write(filename,"strn2D_dyz", dom%mat%now%strn2D%dyz,  units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n]) 
        call nc_write(filename,"strn2D_de",  dom%mat%now%strn2D%de,   units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n]) 
        call nc_write(filename,"strn2D_div", dom%mat%now%strn2D%div,  units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n]) 
        call nc_write(filename,"strn2D_f_shear",dom%mat%now%strn2D%f_shear,units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n]) 
        
        call nc_write(filename,"strn_dxx",     dom%mat%now%strn%dxx,        units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n]) 
        call nc_write(filename,"strn_dyy",     dom%mat%now%strn%dyy,        units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n]) 
        call nc_write(filename,"strn_dxy",     dom%mat%now%strn%dxy,        units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n]) 
        call nc_write(filename,"strn_dxz",     dom%mat%now%strn%dxz,        units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n]) 
        call nc_write(filename,"strn_dyz",     dom%mat%now%strn%dyz,        units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n]) 
        call nc_write(filename,"strn_de",      dom%mat%now%strn%de,         units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n]) 
        call nc_write(filename,"strn_div",     dom%mat%now%strn%div,        units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n]) 
        call nc_write(filename,"strn_f_shear", dom%mat%now%strn%f_shear,    units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n]) 

        call nc_write(filename,"strs2D_txx", dom%mat%now%strs2D%txx,  units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n]) 
        call nc_write(filename,"strs2D_tyy", dom%mat%now%strs2D%tyy,  units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n]) 
        call nc_write(filename,"strs2D_txy", dom%mat%now%strs2D%txy,  units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n]) 
        call nc_write(filename,"strs2D_txz", dom%mat%now%strs2D%txz,  units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n]) 
        call nc_write(filename,"strs2D_tyz", dom%mat%now%strs2D%tyz,  units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n]) 
        call nc_write(filename,"strs2D_te",  dom%mat%now%strs2D%te,   units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n]) 
        call nc_write(filename,"strs2D_tau_eig_1",dom%mat%now%strs2D%tau_eig_1,units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n]) 
        call nc_write(filename,"strs2D_tau_eig_2",dom%mat%now%strs2D%tau_eig_2,units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n]) 
        
        call nc_write(filename,"strs_txx",     dom%mat%now%strs%txx,        units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n]) 
        call nc_write(filename,"strs_tyy",     dom%mat%now%strs%tyy,        units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n]) 
        call nc_write(filename,"strs_txy",     dom%mat%now%strs%txy,        units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n]) 
        call nc_write(filename,"strs_txz",     dom%mat%now%strs%txz,        units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n]) 
        call nc_write(filename,"strs_tyz",     dom%mat%now%strs%tyz,        units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n]) 
        call nc_write(filename,"strs_te",      dom%mat%now%strs%te,         units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n]) 

        call nc_write(filename,"enh",         dom%mat%now%enh,           units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n]) 
        call nc_write(filename,"enh_bnd",     dom%mat%now%enh_bnd,       units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n]) 
        call nc_write(filename,"enh_bar",     dom%mat%now%enh_bar,       units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n]) 
        call nc_write(filename,"ATT",         dom%mat%now%ATT,           units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n]) 
        call nc_write(filename,"ATT_bar",     dom%mat%now%ATT_bar,       units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n]) 
        call nc_write(filename,"visc",        dom%mat%now%visc,          units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n]) 
        call nc_write(filename,"visc_int",    dom%mat%now%visc_int,      units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n]) 

        call nc_write(filename,"f_shear_bar", dom%mat%now%f_shear_bar,   units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n]) 

        call nc_write(filename,"dep_time",    dom%mat%now%dep_time,      units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n])     
        call nc_write(filename,"depth_iso",   dom%mat%now%depth_iso,     units="",dim1="xc",dim2="yc",dim3="age_iso",dim4="time",ncid=ncid,start=[1,1,1,n])     
        
        ! == ytherm variables ===

        call nc_write(filename,"enth",        dom%thrm%now%enth,       units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n])      
        call nc_write(filename,"T_ice",       dom%thrm%now%T_ice,      units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n])      
        call nc_write(filename,"omega",       dom%thrm%now%omega,      units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n])      
        call nc_write(filename,"T_pmp",       dom%thrm%now%T_pmp,      units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n])      

        call nc_write(filename,"f_pmp",       dom%thrm%now%f_pmp,      units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])        
        call nc_write(filename,"bmb_grnd",    dom%thrm%now%bmb_grnd,   units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])     
        call nc_write(filename,"Q_strn",      dom%thrm%now%Q_strn,     units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n])     
        call nc_write(filename,"dQsdt",       dom%thrm%now%dQsdt,      units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n])     
        call nc_write(filename,"Q_b",         dom%thrm%now%Q_b,        units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])          
        call nc_write(filename,"Q_ice_b",     dom%thrm%now%Q_ice_b,    units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])          
        call nc_write(filename,"T_prime_b",   dom%thrm%now%T_prime_b,  units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])    
        call nc_write(filename,"H_w",         dom%thrm%now%H_w,        units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])
        call nc_write(filename,"dHwdt",       dom%thrm%now%dHwdt,      units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])
       
        call nc_write(filename,"cp",          dom%thrm%now%cp,         units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n])         
        call nc_write(filename,"kt",          dom%thrm%now%kt,         units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n])         
        call nc_write(filename,"H_cts",       dom%thrm%now%H_cts,      units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])      
        
        call nc_write(filename,"advecxy",     dom%thrm%now%advecxy,    units="",dim1="xc",dim2="yc",dim3="zeta",dim4="time",ncid=ncid,start=[1,1,1,n])      

        call nc_write(filename,"Q_rock",      dom%thrm%now%Q_rock,     units="mW m-2",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])        
        call nc_write(filename,"enth_rock",   dom%thrm%now%enth_rock,  units="J m-3", dim1="xc",dim2="yc",dim3="zeta_rock",dim4="time",ncid=ncid,start=[1,1,1,n])      
        call nc_write(filename,"T_rock",      dom%thrm%now%T_rock,     units="K",     dim1="xc",dim2="yc",dim3="zeta_rock",dim4="time",ncid=ncid,start=[1,1,1,n])      

        ! == ybound variables ===

        call nc_write(filename,"z_bed",       dom%bnd%z_bed,       units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])
        call nc_write(filename,"z_bed_sd",    dom%bnd%z_bed_sd,    units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])
        call nc_write(filename,"z_sl",        dom%bnd%z_sl,        units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])
        call nc_write(filename,"H_sed",       dom%bnd%H_sed,       units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])
        call nc_write(filename,"smb_ref",     dom%bnd%smb,         units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])
        call nc_write(filename,"T_srf",       dom%bnd%T_srf,       units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])
        call nc_write(filename,"bmb_shlf",    dom%bnd%bmb_shlf,    units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])
        call nc_write(filename,"fmb_shlf",    dom%bnd%fmb_shlf,    units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])
        call nc_write(filename,"T_shlf",      dom%bnd%T_shlf,      units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])
        call nc_write(filename,"Q_geo",       dom%bnd%Q_geo,       units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])
        call nc_write(filename,"enh_srf",     dom%bnd%enh_srf,     units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])

        call nc_write(filename,"basins",      dom%bnd%basins,      units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])
        call nc_write(filename,"basin_mask",  dom%bnd%basin_mask,  units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])
        call nc_write(filename,"regions",     dom%bnd%regions,     units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])
        call nc_write(filename,"region_mask", dom%bnd%region_mask, units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])

        call nc_write(filename,"ice_allowed", dom%bnd%ice_allowed, units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])
        call nc_write(filename,"calv_mask",   dom%bnd%calv_mask,   units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])
        
        call nc_write(filename,"H_ice_ref",   dom%bnd%H_ice_ref,   units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])
        call nc_write(filename,"z_bed_ref",   dom%bnd%z_bed_ref,   units="",dim1="xc",dim2="yc",dim3="time",ncid=ncid,start=[1,1,n])

        nt = nc_size(filename,"time",ncid)
        write(*,*) "nt = ", nt
        
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

    subroutine yelmo_restart_read_topo_bnd(tpo,bnd,tme,restart_interpolated,grd,domain,grid_name,filename,time)  
        ! Load yelmo variables from restart file: [tpo] 
        ! [dyn,therm,mat] variables loaded using yelmo_restart_read
        
        implicit none 

        type(ytopo_class),  intent(INOUT) :: tpo 
        type(ybound_class), intent(INOUT) :: bnd 
        type(ytime_class),  intent(INOUT) :: tme
        integer,            intent(OUT)   :: restart_interpolated
        type(ygrid_class),  intent(IN)    :: grd
        character(len=*),   intent(IN)    :: domain
        character(len=*),   intent(IN)    :: grid_name
        character(len=*),   intent(IN)    :: filename 
        real(wp),           intent(IN)    :: time  

        ! Local variables
        character(len=56) :: restart_domain 
        character(len=56) :: restart_grid_name 
        type(map_scrip_class) :: mps 
        integer  :: nx_restart
        real(wp) :: dx_restart 
        real(wp), allocatable :: xc_restart(:)
        
        ! Load restart file grid attributes 
        if (nc_exists_attr(filename,"domain")) then 
            call nc_read_attr(filename,"domain",    restart_domain)
        else 
            restart_domain = trim(domain)
        end if 

        if (nc_exists_attr(filename,"grid_name")) then 
            call nc_read_attr(filename,"grid_name", restart_grid_name)
        else 
            restart_grid_name = trim(grid_name)
        end if 


        if (trim(restart_grid_name) .eq. trim(grid_name) ) then 
            ! Restart file grid and yelmo grid are the same

            ! Load the data without interpolation (by not specifying mps argument)
            call yelmo_restart_read_topo_bnd_internal(tpo,bnd,tme,filename,time)

            ! Set yelmo flag too
            restart_interpolated = 0

        else 
            ! Restart grid is different than Yelmo grid 

            ! Load the scrip map from file (should already have been generated via cdo externally)
            call map_scrip_init(mps,restart_grid_name,grid_name,method="con",fldr="maps",load=.TRUE.)

            ! Load the data with interpolation
            call yelmo_restart_read_topo_bnd_internal(tpo,bnd,tme,filename,time,mps) 

            ! Determine whether interpolation is from low to high resolution (1)
            ! or from high to low resolution (-1)

            ! Load the x-axis from the restart file 
            ! and determine grid resolution from first points
            nx_restart = nc_size(filename,"xc")
            allocate(xc_restart(nx_restart))
            call nc_read(filename,"xc",xc_restart)
            dx_restart = xc_restart(2) - xc_restart(1) 

            if (dx_restart .lt. grd%dx) then 
                ! Low to high resolution
                restart_interpolated = 1
            else 
                ! High to low resolution
                restart_interpolated = -1 
            end if 

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
        call yelmo_restart_write(dom,"./yelmo_restart_init.nc",time)
        
        return 

    end subroutine yelmo_restart_read
    
    subroutine yelmo_restart_read_topo_bnd_internal(tpo,bnd,tme,filename,time,mps)  
        ! Load yelmo variables from restart file: [tpo] 
        ! [dyn,therm,mat] variables loaded using yelmo_restart_read
        
        implicit none 

        type(ytopo_class),  intent(INOUT) :: tpo 
        type(ybound_class), intent(INOUT) :: bnd 
        type(ytime_class),  intent(INOUT) :: tme
        character(len=*),  intent(IN)    :: filename 
        real(wp),          intent(IN)    :: time 
        type(map_scrip_class), optional, intent(IN) :: mps 

        ! Local variables
        integer  :: ncid, n, nx, ny
        real(wp) :: time_of_restart_file 

        ! Read all yelmo data from file,
        ! in order to restart a simulation.
        
        ! Open the file for reading

        call nc_open(filename,ncid,writable=.FALSE.)

        ! Note: no need to read in the dimension information,
        ! this will be initialized by Yelmo itself

        ! Define dimensions of variables 
        nx    = tpo%par%nx
        ny    = tpo%par%ny

        ! Assume that first time dimension value is to be read in
        n = 1 

        ! == time variables ===

        ! ajr: testing reading these variables too to improve restart file performance
        call nc_read(filename,"pc_dt",       tme%pc_dt, start=[1,n],count=[3,1],ncid=ncid)
        call nc_read(filename,"pc_eta",      tme%pc_eta,start=[1,n],count=[3,1],ncid=ncid)
        
        call nc_read_interp(filename,"pc_tau",       tme%pc_tau,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"pc_tau_masked",tme%pc_tau_masked,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"pc_tau_max",   tme%pc_tau_max,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        
        ! == ytopo variables ===

        call nc_read_interp(filename,"H_ice",       tpo%now%H_ice,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"z_srf",       tpo%now%z_srf,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"z_base",      tpo%now%z_base,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"dzsdt",       tpo%now%dzsdt,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"dHidt",       tpo%now%dHidt,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"dHidt_dyn",   tpo%now%dHidt_dyn,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"mb_net",      tpo%now%mb_net,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"mb_relax",    tpo%now%mb_relax,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"mb_resid",    tpo%now%mb_resid,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"mb_err",      tpo%now%mb_err,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"smb",         tpo%now%smb,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"bmb",         tpo%now%bmb,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"fmb",         tpo%now%fmb,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"dmb",         tpo%now%dmb,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"cmb",         tpo%now%cmb,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"cmb_flt",     tpo%now%cmb_flt,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"cmb_grnd",    tpo%now%cmb_grnd,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"eps_eff",     tpo%now%eps_eff,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"tau_eff",     tpo%now%tau_eff,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"dzsdx",       tpo%now%dzsdx,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)  
        call nc_read_interp(filename,"dzsdy",       tpo%now%dzsdy,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)  
        call nc_read_interp(filename,"dHidx",       tpo%now%dHidx,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"dHidy",       tpo%now%dHidy,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"dzbdx",       tpo%now%dzbdx,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)  
        call nc_read_interp(filename,"dzbdy",       tpo%now%dzbdy,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)  
        call nc_read_interp(filename,"H_eff",       tpo%now%H_eff,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"H_grnd",      tpo%now%H_grnd,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"f_grnd",      tpo%now%f_grnd,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)  
        call nc_read_interp(filename,"f_grnd_acx",  tpo%now%f_grnd_acx,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"f_grnd_acy",  tpo%now%f_grnd_acy,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"f_ice",       tpo%now%f_ice,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)  
        call nc_read_interp(filename,"dist_margin", tpo%now%dist_margin,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"dist_grline", tpo%now%dist_grline,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"mask_bed",    tpo%now%mask_bed,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"mask_grz",    tpo%now%mask_grz,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)

        call nc_read_interp(filename,"dHidt_dyn_n", tpo%now%dHidt_dyn_n,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"H_ice_n",     tpo%now%H_ice_n,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"z_srf_n",     tpo%now%z_srf_n,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        
        call nc_read_interp(filename,"H_ice_dyn",   tpo%now%H_ice_dyn,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"f_ice_dyn",   tpo%now%f_ice_dyn,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)  
        
        ! = ytopo_pc variables ===
        
        call nc_read_interp(filename,"pc_pred_H_ice",    tpo%now%pred%H_ice,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"pc_pred_dHidt_dyn",tpo%now%pred%dHidt_dyn,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"pc_pred_mb_net",   tpo%now%pred%mb_net,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"pc_pred_smb",      tpo%now%pred%smb,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"pc_pred_bmb",      tpo%now%pred%bmb,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"pc_pred_fmb",      tpo%now%pred%fmb,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"pc_pred_dmb",      tpo%now%pred%dmb,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"pc_pred_cmb",      tpo%now%pred%cmb,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"pc_corr_H_ice",    tpo%now%corr%H_ice,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"pc_corr_dHidt_dyn",tpo%now%corr%dHidt_dyn,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"pc_corr_mb_net",   tpo%now%corr%mb_net,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"pc_corr_smb",      tpo%now%corr%smb,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"pc_corr_bmb",      tpo%now%corr%bmb,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"pc_corr_fmb",      tpo%now%corr%fmb,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"pc_corr_dmb",      tpo%now%corr%dmb,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        call nc_read_interp(filename,"pc_corr_cmb",      tpo%now%corr%cmb,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)

        ! == ybound variables ===

        call nc_read_interp(filename,"z_bed",       bnd%z_bed,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"z_bed_sd",    bnd%z_bed_sd,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"z_sl",        bnd%z_sl,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"H_sed",       bnd%H_sed,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"smb_ref",     bnd%smb,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"T_srf",       bnd%T_srf,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"bmb_shlf",    bnd%bmb_shlf,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"fmb_shlf",    bnd%fmb_shlf,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"T_shlf",      bnd%T_shlf,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"Q_geo",       bnd%Q_geo,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"enh_srf",     bnd%enh_srf,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 

        call nc_read_interp(filename,"basins",      bnd%basins,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"basin_mask",  bnd%basin_mask,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"regions",     bnd%regions,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"region_mask", bnd%region_mask,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        
        call nc_read_interp(filename,"ice_allowed", bnd%ice_allowed,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"calv_mask",   bnd%calv_mask,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        
        call nc_read_interp(filename,"H_ice_ref",   bnd%H_ice_ref,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"z_bed_ref",   bnd%z_bed_ref,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 

        ! Close the netcdf file
        call nc_close(ncid)
        
        tpo%par%time = time
        
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
        
        ! == ytopo variables ===

        ! Reload mask_bed since it contains thermodynamic information too
        call nc_read_interp(filename,"mask_bed",      dom%tpo%now%mask_bed,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)
        
        ! == ydyn variables ===

        call nc_read_interp(filename,"ux",            dom%dyn%now%ux,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 
        call nc_read_interp(filename,"uy",            dom%dyn%now%uy,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 
        call nc_read_interp(filename,"uxy",           dom%dyn%now%uxy,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 
        call nc_read_interp(filename,"uz",            dom%dyn%now%uz,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz_ac,1],mps=mps) 
        call nc_read_interp(filename,"uz_star",       dom%dyn%now%uz_star,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz_ac,1],mps=mps) 
      
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

        call nc_read_interp(filename,"taul_int_acx",  dom%dyn%now%taul_int_acx,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"taul_int_acy",  dom%dyn%now%taul_int_acy,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 

        call nc_read_interp(filename,"qq_gl_acx",     dom%dyn%now%qq_gl_acx,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"qq_gl_acy",     dom%dyn%now%qq_gl_acy,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        
        call nc_read_interp(filename,"qq_acx",        dom%dyn%now%qq_acx,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"qq_acy",        dom%dyn%now%qq_acy,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
        call nc_read_interp(filename,"qq",            dom%dyn%now%qq,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 

        call nc_read_interp(filename,"de_eff",        dom%dyn%now%de_eff,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 
        call nc_read_interp(filename,"visc_eff",      dom%dyn%now%visc_eff,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 
        call nc_read_interp(filename,"visc_eff_int",  dom%dyn%now%visc_eff_int,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 

        call nc_read_interp(filename,"N_eff",         dom%dyn%now%N_eff,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps)        
        call nc_read_interp(filename,"cb_tgt",        dom%dyn%now%cb_tgt,ncid=ncid,start=[1,1,n],count=[nx,ny,1],mps=mps) 
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
        
        call nc_read_interp(filename,"jvel_dxx",      dom%dyn%now%jvel%dxx, ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 
        call nc_read_interp(filename,"jvel_dxy",      dom%dyn%now%jvel%dxy, ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 
        call nc_read_interp(filename,"jvel_dxz",      dom%dyn%now%jvel%dxz, ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 
        call nc_read_interp(filename,"jvel_dyx",      dom%dyn%now%jvel%dyx, ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 
        call nc_read_interp(filename,"jvel_dyy",      dom%dyn%now%jvel%dyy, ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 
        call nc_read_interp(filename,"jvel_dyz",      dom%dyn%now%jvel%dyz, ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps) 
        call nc_read_interp(filename,"jvel_dzx",      dom%dyn%now%jvel%dzx, ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz_ac,1],mps=mps) 
        call nc_read_interp(filename,"jvel_dzy",      dom%dyn%now%jvel%dzy, ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz_ac,1],mps=mps)
        call nc_read_interp(filename,"jvel_dzz",      dom%dyn%now%jvel%dzz, ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz_ac,1],mps=mps) 

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
        call nc_read_interp(filename,"dQsdt",       dom%thrm%now%dQsdt,ncid=ncid,start=[1,1,1,n],count=[nx,ny,nz,1],mps=mps)      
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
        dom%dyn%par%time  = time
        
        write(*,*) 
        write(*,*) "time = ", time, " : loaded restart file: ", trim(filename)
        write(*,*) 
        
        return 

    end subroutine yelmo_restart_read_internal
    
end module yelmo_io


