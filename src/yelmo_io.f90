

module yelmo_io
    
    use ncio 
    
    use yelmo_defs 
    use yelmo_grid 

    implicit none

    private 
    public :: yelmo_write_init
    public :: yelmo_restart_write
    public :: yelmo_restart_read_1
    public :: yelmo_restart_read_2


contains

    subroutine yelmo_write_init(ylmo,filename,time_init,units)

        implicit none 

        type(yelmo_class), intent(IN) :: ylmo 
        character(len=*),  intent(IN) :: filename, units 
        real(prec),        intent(IN) :: time_init

        ! Initialize file by writing grid info
        call yelmo_grid_write(ylmo%grd,filename,create=.TRUE.)

        ! Initialize netcdf file and dimensions
        call nc_write_dim(filename,"month",  x=1,dx=1,nx=12,         units="month")
        call nc_write_dim(filename,"zeta",   x=ylmo%par%zeta_aa,     units="1")
        call nc_write_dim(filename,"zeta_ac",x=ylmo%par%zeta_ac,     units="1")
        call nc_write_dim(filename,"time",   x=time_init,dx=1.0_prec,nx=1,units=trim(units),unlimited=.TRUE.)

        ! Static information
        call nc_write(filename,"basins",  ylmo%bnd%basins, dim1="xc",dim2="yc",units="(0 - 8)",long_name="Hydrological basins")
        call nc_write(filename,"regions", ylmo%bnd%regions,dim1="xc",dim2="yc",units="(0 - 8)",long_name="Domain regions")
        
        ! Additional optional static information 
        call nc_write(filename,"z_bed_sd", ylmo%bnd%z_bed_sd,dim1="xc",dim2="yc",units="m",long_name="Stdev(z_bed)")
        
        return

    end subroutine yelmo_write_init 

    subroutine yelmo_restart_write(dom,filename,time)

        implicit none 

        type(yelmo_class), intent(IN) :: dom 
        character(len=*),  intent(IN) :: filename 
        real(prec),        intent(IN) :: time 
        
        ! Local variables
        integer :: ncid 
        
        ! Write all yelmo data to file, so that it can be
        ! read later to restart a simulation.
        
        ! == Initialize netcdf file ==============================================

        call nc_create(filename)
        call nc_write_dim(filename,"xc",       x=dom%grd%xc*1e-3,       units="kilometers")
        call nc_write_dim(filename,"yc",       x=dom%grd%yc*1e-3,       units="kilometers")
        call nc_write_dim(filename,"zeta",     x=dom%par%zeta_aa,       units="1")
        call nc_write_dim(filename,"zeta_ac",  x=dom%par%zeta_ac,       units="1")
        call nc_write_dim(filename,"time",     x=time,dx=1.0_prec,nx=1, units="years ago")
        
        ! == Begin writing data ==============================================
        
        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)
        
        ! == ytopo variables ===
        call nc_write(filename,"H_ice",   dom%tpo%now%H_ice,   units="m",  dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"z_srf",   dom%tpo%now%z_srf,   units="m",  dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"dzsrfdt", dom%tpo%now%dzsrfdt, units="m/a",dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"dHicedt", dom%tpo%now%dHicedt, units="m/a",dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"bmb",     dom%tpo%now%bmb,     units="m/a",dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"mb_applied",dom%tpo%now%mb_applied,units="m/a",dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"calv",    dom%tpo%now%calv,    units="m/a",dim1="xc",dim2="yc",ncid=ncid)
        
        call nc_write(filename,"dzsdx",   dom%tpo%now%dzsdx,   units="m/m",  dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"dzsdy",   dom%tpo%now%dzsdy,   units="m/m",  dim1="xc",dim2="yc",ncid=ncid)
        
        call nc_write(filename,"f_grnd",   dom%tpo%now%f_grnd,   units="1",dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"f_ice",    dom%tpo%now%f_ice,    units="1",dim1="xc",dim2="yc",ncid=ncid)

        call nc_write(filename,"mask_bed", dom%tpo%now%mask_bed, units="1",dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"is_float", dom%tpo%now%is_float, units="1",dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"is_grline",dom%tpo%now%is_grline,units="1",dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"is_grz",   dom%tpo%now%is_grz,   units="1",dim1="xc",dim2="yc",ncid=ncid)
        
        ! == ydyn variables ===

        ! ajr: to do !!!
!         call nc_write(filename,"ux_i",    dom%dyn%now%ux_i,units="m/a",dim1="xc",dim2="yc",dim3="zeta",ncid=ncid)
!         call nc_write(filename,"uy_i",    dom%dyn%now%uy_i,units="m/a",dim1="xc",dim2="yc",dim3="zeta",ncid=ncid)
!         call nc_write(filename,"ux_i_bar",dom%dyn%now%ux_i_bar,units="m/a",dim1="xc",dim2="yc",ncid=ncid)
!         call nc_write(filename,"uy_i_bar",dom%dyn%now%uy_i_bar,units="m/a",dim1="xc",dim2="yc",ncid=ncid)
!         call nc_write(filename,"ux_b",    dom%dyn%now%ux_b,units="m/a",dim1="xc",dim2="yc",ncid=ncid)
!         call nc_write(filename,"uy_b",    dom%dyn%now%uy_b,units="m/a",dim1="xc",dim2="yc",ncid=ncid)
!         call nc_write(filename,"ux_bar",  dom%dyn%now%ux_bar,units="m/a",dim1="xc",dim2="yc",ncid=ncid)
!         call nc_write(filename,"uy_bar",  dom%dyn%now%uy_bar,units="m/a",dim1="xc",dim2="yc",ncid=ncid)

!         call nc_write(filename,"ux",    dom%dyn%now%ux,units="m/a",dim1="xc",dim2="yc",dim3="zeta",ncid=ncid)
!         call nc_write(filename,"uy",    dom%dyn%now%uy,units="m/a",dim1="xc",dim2="yc",dim3="zeta",ncid=ncid)
!         call nc_write(filename,"uz",    dom%dyn%now%uz,units="m/a",dim1="xc",dim2="yc",dim3="zeta_ac",ncid=ncid)
        
!         call nc_write(filename,"f_vbvs",    dom%dyn%now%f_vbvs,units="1",dim1="xc",dim2="yc",ncid=ncid)
        
!         call nc_write(filename,"C_bed",  dom%dyn%now%C_bed,  units="?",dim1="xc",dim2="yc",ncid=ncid)
!         call nc_write(filename,"N_eff",    dom%dyn%now%N_eff,    units="bar",dim1="xc",dim2="yc",ncid=ncid)
!         call nc_write(filename,"beta_c",dom%dyn%now%beta_c,units="?",dim1="xc",dim2="yc",ncid=ncid)
        
!         call nc_write(filename,"tau_b",  dom%dyn%now%tau_b,  units="?",dim1="xc",dim2="yc",ncid=ncid)

        ! == ymat variables ===
        call nc_write(filename,"strn2D_dxx",dom%mat%now%strn2D%dxx,units="m/m",dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"strn2D_dyy",dom%mat%now%strn2D%dyy,units="m/m",dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"strn2D_dxy",dom%mat%now%strn2D%dxy,units="m/m",dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"strn2D_de", dom%mat%now%strn2D%de, units="m/m",dim1="xc",dim2="yc",ncid=ncid)
        
        call nc_write(filename,"strn_dxx",dom%mat%now%strn%dxx,units="m/m",dim1="xc",dim2="yc",dim3="zeta",ncid=ncid)
        call nc_write(filename,"strn_dyy",dom%mat%now%strn%dyy,units="m/m",dim1="xc",dim2="yc",dim3="zeta",ncid=ncid)
        call nc_write(filename,"strn_dxy",dom%mat%now%strn%dxy,units="m/m",dim1="xc",dim2="yc",dim3="zeta",ncid=ncid)
        call nc_write(filename,"strn_dxz",dom%mat%now%strn%dxz,units="m/m",dim1="xc",dim2="yc",dim3="zeta",ncid=ncid)
        call nc_write(filename,"strn_dyz",dom%mat%now%strn%dyz,units="m/m",dim1="xc",dim2="yc",dim3="zeta",ncid=ncid)
        call nc_write(filename,"strn_de", dom%mat%now%strn%de, units="m/m",dim1="xc",dim2="yc",dim3="zeta",ncid=ncid)
        call nc_write(filename,"strn_f_shear",dom%mat%now%strn%f_shear,units="1",dim1="xc",dim2="yc",dim3="zeta",ncid=ncid)

        call nc_write(filename,"f_shear_bar",dom%mat%now%f_shear_bar,units="",dim1="xc",dim2="yc",ncid=ncid)

        call nc_write(filename,"enh",      dom%mat%now%enh,         units="",dim1="xc",dim2="yc",dim3="zeta",ncid=ncid)
        call nc_write(filename,"enh_bar",  dom%mat%now%enh_bar,     units="",dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"ATT",      dom%mat%now%ATT,         units="",dim1="xc",dim2="yc",dim3="zeta",ncid=ncid)
        call nc_write(filename,"ATT_bar",  dom%mat%now%ATT_bar,     units="",dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"visc",     dom%mat%now%visc,        units="",dim1="xc",dim2="yc",dim3="zeta",ncid=ncid)
        call nc_write(filename,"visc_int", dom%mat%now%visc_int,    units="",dim1="xc",dim2="yc",ncid=ncid)
        
        ! == ytherm variables ===
        call nc_write(filename,"T_ice",     dom%thrm%now%T_ice,     units="",dim1="xc",dim2="yc",dim3="zeta",ncid=ncid)
        call nc_write(filename,"T_pmp",     dom%thrm%now%T_pmp,     units="",dim1="xc",dim2="yc",dim3="zeta",ncid=ncid)
        call nc_write(filename,"phid",      dom%thrm%now%phid,      units="",dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"f_pmp",     dom%thrm%now%f_pmp,     units="",dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"bmb_grnd",  dom%thrm%now%bmb_grnd,  units="",dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"Q_b",       dom%thrm%now%Q_b,       units="",dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"Q_strn",    dom%thrm%now%Q_strn,    units="",dim1="xc",dim2="yc",dim3="zeta", ncid=ncid)
        call nc_write(filename,"T_rock",    dom%thrm%now%T_rock,    units="",dim1="xc",dim2="yc",dim3="sr",ncid=ncid)
        
        ! == ybound variables ===
        call nc_write(filename,"z_bed",   dom%bnd%z_bed,   units="m",    dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"z_sl",    dom%bnd%z_sl,    units="m",    dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"H_sed",   dom%bnd%H_sed,   units="m",    dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"H_w",     dom%bnd%H_w,     units="m",    dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"smb",     dom%bnd%smb,     units="m/a",  dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"T_srf",   dom%bnd%T_srf,   units="K",    dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"bmb_shlf",dom%bnd%bmb_shlf,units="m/a",  dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"T_shlf",  dom%bnd%T_shlf,  units="K",    dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"Q_geo",   dom%bnd%Q_geo,   units="mW/m2",dim1="xc",dim2="yc",ncid=ncid)

        call nc_write(filename,"basins",     dom%bnd%basins,     units="",dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"basin_mask", dom%bnd%basin_mask, units="",dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"regions",    dom%bnd%regions,    units="",dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"region_mask",dom%bnd%region_mask,units="",dim1="xc",dim2="yc",ncid=ncid)


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

    subroutine yelmo_restart_read_1(dom,filename,time)
        ! Load yelmo variables from restart file: [tpo] 
        ! [dyn,therm,mat] variables loaded using yelmo_restart_read_2
        
        implicit none 

        type(yelmo_class), intent(INOUT) :: dom 
        character(len=*),  intent(IN)    :: filename 
        real(prec),        intent(IN)    :: time 
        
        ! Local variables
        integer :: ncid
        
        ! Read all yelmo data from file,
        ! in order to restart a simulation.
        
        ! Open the file for writing
        call nc_open(filename,ncid,writable=.FALSE.)
        
                ! == ytopo variables ===
        call nc_read(filename,"H_ice",   dom%tpo%now%H_ice,   ncid=ncid)
        call nc_read(filename,"z_srf",   dom%tpo%now%z_srf,   ncid=ncid)
        call nc_read(filename,"dzsrfdt", dom%tpo%now%dzsrfdt, ncid=ncid)
        call nc_read(filename,"dHicedt", dom%tpo%now%dHicedt, ncid=ncid)
        call nc_read(filename,"bmb",     dom%tpo%now%bmb,     ncid=ncid)
        call nc_read(filename,"mb_applied",dom%tpo%now%mb_applied,ncid=ncid)
        call nc_read(filename,"calv",    dom%tpo%now%calv,    ncid=ncid)
        
        call nc_read(filename,"dzsdx",   dom%tpo%now%dzsdx,   ncid=ncid)
        call nc_read(filename,"dzsdy",   dom%tpo%now%dzsdy,   ncid=ncid)
        
        call nc_read(filename,"f_grnd",   dom%tpo%now%f_grnd,   ncid=ncid)
        call nc_read(filename,"f_ice",    dom%tpo%now%f_ice,    ncid=ncid)
        
        call nc_read(filename,"mask_bed", dom%tpo%now%mask_bed, ncid=ncid)
        call nc_read(filename,"is_float", dom%tpo%now%is_float, ncid=ncid)
        call nc_read(filename,"is_grline",dom%tpo%now%is_grline,ncid=ncid)
        call nc_read(filename,"is_grz",   dom%tpo%now%is_grz,   ncid=ncid)
        
        ! Close the netcdf file
        call nc_close(ncid)

        ! Write summary 
        write(*,*) 
        write(*,*) "time = ", time, " : loaded restart file: ", trim(filename)
        write(*,*) 
        
        return 

    end subroutine yelmo_restart_read_1 


    subroutine yelmo_restart_read_2(dom,filename,time)
        ! Load yelmo variables from restart file: [dyn,therm,mat] 
        ! [tpo] variables loaded using yelmo_restart_read_1

        implicit none 

        type(yelmo_class), intent(INOUT) :: dom 
        character(len=*),  intent(IN)    :: filename 
        real(prec),        intent(IN)    :: time 
        
        ! Local variables
        integer :: ncid
        
        ! Read all yelmo data from file,
        ! in order to restart a simulation.
        
        ! Open the file for writing
        call nc_open(filename,ncid,writable=.FALSE.)
        
        ! == ydyn variables ===

        ! ajr: to do !!!
!         call nc_read(filename,"ux_sia",dom%dyn%now%ux_sia,ncid=ncid)
!         call nc_read(filename,"uy_sia",dom%dyn%now%uy_sia,ncid=ncid)
!         call nc_read(filename,"ux_ssa",dom%dyn%now%ux_ssa,ncid=ncid)
!         call nc_read(filename,"uy_ssa",dom%dyn%now%uy_ssa,ncid=ncid)
!         call nc_read(filename,"ux_bar",dom%dyn%now%ux_bar,ncid=ncid)
!         call nc_read(filename,"uy_bar",dom%dyn%now%uy_bar,ncid=ncid)
!         call nc_read(filename,"diff",  dom%dyn%now%diff,  ncid=ncid)
        
!         call nc_read(filename,"ux",    dom%dyn%now%ux,ncid=ncid)
!         call nc_read(filename,"uy",    dom%dyn%now%uy,ncid=ncid)
!         call nc_read(filename,"uz",    dom%dyn%now%uz,ncid=ncid)
        
!         call nc_read(filename,"f_ssa",     dom%dyn%now%f_ssa,     ncid=ncid)
!         call nc_read(filename,"f_ssa_mx",  dom%dyn%now%f_ssa_mx,  ncid=ncid)
!         call nc_read(filename,"f_ssa_my",  dom%dyn%now%f_ssa_my,  ncid=ncid)
!         call nc_read(filename,"ssa_active",dom%dyn%now%ssa_active,ncid=ncid)
        
!         call nc_read(filename,"f_vbvs",    dom%dyn%now%f_vbvs,    ncid=ncid)

!         call nc_read(filename,"N_eff",    dom%dyn%now%N_eff,    ncid=ncid)

!         call nc_read(filename,"beta",      dom%dyn%now%beta,      ncid=ncid)
!         call nc_read(filename,"beta_c",    dom%dyn%now%beta_c,    ncid=ncid)
        
!         call nc_read(filename,"tau_mx",  dom%dyn%now%tau_mx,      ncid=ncid)
!         call nc_read(filename,"tau_my",  dom%dyn%now%tau_my,      ncid=ncid)
        
        ! == ymat variables ===
        call nc_read(filename,"strn2D_dxx",dom%mat%now%strn2D%dxx,ncid=ncid)
        call nc_read(filename,"strn2D_dyy",dom%mat%now%strn2D%dyy,ncid=ncid)
        call nc_read(filename,"strn2D_dxy",dom%mat%now%strn2D%dxy,ncid=ncid)
        call nc_read(filename,"strn2D_de", dom%mat%now%strn2D%de, ncid=ncid)
        
        call nc_read(filename,"strn_dxx",dom%mat%now%strn%dxx,        ncid=ncid)
        call nc_read(filename,"strn_dyy",dom%mat%now%strn%dyy,        ncid=ncid)
        call nc_read(filename,"strn_dxy",dom%mat%now%strn%dxy,        ncid=ncid)
        call nc_read(filename,"strn_dxz",dom%mat%now%strn%dxz,        ncid=ncid)
        call nc_read(filename,"strn_dyz",dom%mat%now%strn%dyz,        ncid=ncid)
        call nc_read(filename,"strn_de", dom%mat%now%strn%de,         ncid=ncid)
        call nc_read(filename,"strn_f_shear",dom%mat%now%strn%f_shear,ncid=ncid)

        call nc_read(filename,"f_shear_bar",dom%mat%now%f_shear_bar,  ncid=ncid)

        call nc_read(filename,"enh",      dom%mat%now%enh,            ncid=ncid)
        call nc_read(filename,"enh_bar",  dom%mat%now%enh_bar,        ncid=ncid)
        call nc_read(filename,"ATT",      dom%mat%now%ATT,            ncid=ncid)
        call nc_read(filename,"ATT_bar",  dom%mat%now%ATT_bar,        ncid=ncid)
        call nc_read(filename,"visc",     dom%mat%now%visc,           ncid=ncid)
        call nc_read(filename,"visc_int", dom%mat%now%visc_int,       ncid=ncid)
        
        ! == ytherm variables ===
        call nc_read(filename,"T_ice",     dom%thrm%now%T_ice,        ncid=ncid)
        call nc_read(filename,"T_pmp",     dom%thrm%now%T_pmp,        ncid=ncid)
        call nc_read(filename,"phid",      dom%thrm%now%phid,         ncid=ncid)
        call nc_read(filename,"f_pmp",     dom%thrm%now%f_pmp,        ncid=ncid)
        call nc_read(filename,"bmb_grnd",  dom%thrm%now%bmb_grnd,     ncid=ncid)
        call nc_read(filename,"Q_b",       dom%thrm%now%Q_b,          ncid=ncid)
        call nc_read(filename,"Q_strn",    dom%thrm%now%Q_strn,       ncid=ncid)
        call nc_read(filename,"T_rock",    dom%thrm%now%T_rock,       ncid=ncid)
        
        ! == ybound variables ===
        call nc_read(filename,"z_bed",   dom%bnd%z_bed,               ncid=ncid)
        call nc_read(filename,"z_sl",    dom%bnd%z_sl,                ncid=ncid)
        call nc_read(filename,"H_sed",   dom%bnd%H_sed,               ncid=ncid)
        call nc_read(filename,"H_w",     dom%bnd%H_w,                 ncid=ncid)
        call nc_read(filename,"smb",     dom%bnd%smb,                 ncid=ncid)
        call nc_read(filename,"T_srf",   dom%bnd%T_srf,               ncid=ncid)
        call nc_read(filename,"bmb_shlf",dom%bnd%bmb_shlf,            ncid=ncid)
        call nc_read(filename,"T_shlf",  dom%bnd%T_shlf,              ncid=ncid)
        call nc_read(filename,"Q_geo",   dom%bnd%Q_geo,               ncid=ncid)

        call nc_read(filename,"basins",     dom%bnd%basins,           ncid=ncid)
        call nc_read(filename,"basin_mask", dom%bnd%basin_mask,       ncid=ncid)
        call nc_read(filename,"regions",    dom%bnd%regions,          ncid=ncid)
        call nc_read(filename,"region_mask",dom%bnd%region_mask,      ncid=ncid)

        ! Close the netcdf file
        call nc_close(ncid)

        ! Write summary 
        write(*,*) 
        write(*,*) "time = ", time, " : loaded restart file: ", trim(filename)
        write(*,*) 
        
        return 

    end subroutine yelmo_restart_read_2 
    

end module yelmo_io 


