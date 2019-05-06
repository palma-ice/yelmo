
module yelmo_regions

    use ncio 
    use yelmo_defs 

    implicit none
    
    private
    public :: calc_yregions
    public :: write_yreg_init
    public :: write_yreg_step 

contains

    subroutine yregions_init()

        implicit none 




        return 

    end subroutine yregions_init 

    subroutine calc_yregions(reg,tpo,dyn,thrm,mat,bnd,mask) 
        ! Calculate a set of regional variables (averages, totals)
        ! for a given domain defined by `mask`.

        implicit none 

        type(yregions_class), intent(INOUT) :: reg 
        type(ytopo_class),    intent(IN)    :: tpo 
        type(ydyn_class),     intent(IN)    :: dyn
        type(ytherm_class),   intent(IN)    :: thrm 
        type(ymat_class),     intent(IN)    :: mat
        type(ybound_class),   intent(IN)    :: bnd  
        logical,              intent(IN)    :: mask(:,:) 

        ! Local variables
        integer    :: nx, ny
        real(prec) :: npts_tot, npts_grnd, npts_flt 

        real(prec) :: m3_km3 = 1e-9 
        real(prec) :: m2_km2 = 1e-6 
        real(prec) :: conv_km3a_Sv

        logical, allocatable :: mask_tot(:,:) 
        logical, allocatable :: mask_grnd(:,:) 
        logical, allocatable :: mask_flt(:,:) 
        
        ! Conversion parameter 
        conv_km3a_Sv = 1e-6*(1e9*rho_w/rho_ice)/sec_year

        ! Grid size 
        nx = size(mask,1)
        ny = size(mask,2)

        ! Allocate masks 
        allocate(mask_tot(nx,ny))
        allocate(mask_grnd(nx,ny))
        allocate(mask_flt(nx,ny))
        
        ! Define masks 
        mask_tot  = (mask .and. tpo%now%H_ice .gt. 0.0) 
        mask_grnd = (mask .and. tpo%now%H_ice .gt. 0.0 .and. tpo%now%f_grnd .gt. 0.0)
        mask_flt  = (mask .and. tpo%now%H_ice .gt. 0.0 .and. tpo%now%f_grnd .eq. 0.0)
         
        
        npts_tot  = real(count(mask_tot),prec)
        npts_grnd = real(count(mask_grnd),prec)
        npts_flt  = real(count(mask_flt),prec)
        
        ! ===== Total ice variables =====

        if (npts_tot .gt. 0) then 

            ! ytopo variables 
            reg%H_ice      = sum(tpo%now%H_ice,mask=mask_tot)/npts_tot     ! [m]
            reg%z_srf      = sum(tpo%now%z_srf,mask=mask_tot)/npts_tot     ! [m]
            reg%dHicedt    = sum(tpo%now%dHicedt,mask=mask_tot)/npts_tot   ! [m/a]
            reg%H_ice_max  = maxval(tpo%now%H_ice,mask=mask_tot)       ! [m]
            reg%dzsrfdt    = sum(tpo%now%dzsrfdt,mask=mask_tot)/npts_tot   ! [m/a]
            
            reg%V_ice      = sum(tpo%now%H_ice,mask=mask_tot)*tpo%par%dx*tpo%par%dy*m3_km3             ! [km^3]
            reg%A_ice      = count(tpo%now%H_ice .gt. 0.0 .and. mask_tot)*tpo%par%dx*tpo%par%dy*m2_km2 ! [km^2]
            reg%dVicedt    = sum(tpo%now%dHicedt,mask=mask_tot)*tpo%par%dx*tpo%par%dy*m3_km3           ! [km^3/a]
            reg%fwf        = -reg%dVicedt*conv_km3a_Sv   ! [Sv]

            ! ydyn variables 
            reg%uxy_bar    = sum(dyn%now%uxy_bar,mask=mask_tot)/npts_tot      ! [m/a]
            reg%uxy_s      = sum(dyn%now%uxy_s,mask=mask_tot)/npts_tot        ! [m/a]
            reg%uxy_b      = sum(dyn%now%uxy_b,mask=mask_tot)/npts_tot        ! [m/a]
            
            ! Boundary variables
            reg%z_bed      = sum(bnd%z_bed,mask=mask_tot)/npts_tot
            reg%smb        = sum(bnd%smb,mask=mask_tot)/npts_tot
            reg%T_srf      = sum(bnd%T_srf,mask=mask_tot)/npts_tot
            reg%bmb        = sum(tpo%now%bmb,mask=mask_tot)/npts_tot
            
        else 

            ! ytopo variables 
            reg%H_ice      = 0.0_prec 
            reg%z_srf      = 0.0_prec 
            reg%dHicedt    = 0.0_prec 
            reg%H_ice_max  = 0.0_prec 
            reg%dzsrfdt    = 0.0_prec 
            
            reg%V_ice      = 0.0_prec 
            reg%A_ice      = 0.0_prec 
            reg%dVicedt    = 0.0_prec 
            reg%fwf        = 0.0_prec 
            
            ! ydyn variables 
            reg%uxy_bar    = 0.0_prec 
            reg%uxy_s      = 0.0_prec 
            reg%uxy_b      = 0.0_prec 
            
            ! Boundary variables
            reg%z_bed      = 0.0_prec 
            reg%smb        = 0.0_prec 
            reg%T_srf      = 0.0_prec 
            reg%bmb        = 0.0_prec 
            
        end if 

        ! ===== Grounded ice variables =====

        if (npts_grnd .gt. 0) then 

            ! ytopo variables 
            reg%H_ice_g      = sum(tpo%now%H_ice,mask=mask_grnd)/npts_grnd     ! [m]
            reg%z_srf_g      = sum(tpo%now%z_srf,mask=mask_grnd)/npts_grnd     ! [m]
            
            reg%V_ice_g      = sum(tpo%now%H_ice,mask=mask_grnd)*tpo%par%dx*tpo%par%dy*m3_km3             ! [km^3]
            reg%A_ice_g      = count(tpo%now%H_ice .gt. 0.0 .and. mask_grnd)*tpo%par%dx*tpo%par%dy*m2_km2 ! [km^2]
            
            ! ydyn variables 
            reg%uxy_bar_g    = sum(dyn%now%uxy_bar,mask=mask_grnd)/npts_grnd      ! [m/a]
            reg%uxy_s_g      = sum(dyn%now%uxy_s,mask=mask_grnd)/npts_grnd        ! [m/a]
            reg%uxy_b_g      = sum(dyn%now%uxy_b,mask=mask_grnd)/npts_grnd        ! [m/a]
            
            ! ythrm variables 
            reg%f_pmp        = sum(thrm%now%f_pmp,mask=mask_grnd)/npts_grnd       ! [fraction]

            ! Boundary variables
            reg%H_w          = sum(bnd%H_w,mask=mask_grnd)/npts_grnd
            reg%bmb_g        = sum(tpo%now%bmb,mask=mask_grnd)/npts_grnd

        else 

            ! ytopo variables 
            reg%H_ice_g      = 0.0_prec 
            reg%z_srf_g      = 0.0_prec 

            reg%V_ice_g      = 0.0_prec 
            reg%A_ice_g      = 0.0_prec 

            ! ydyn variables 
            reg%uxy_bar_g    = 0.0_prec 
            reg%uxy_s_g      = 0.0_prec 
            reg%uxy_b_g      = 0.0_prec 
            
            ! ythrm variables 
            reg%f_pmp        = 0.0_prec 

            ! Boundary variables
            reg%H_w          = 0.0_prec 
            reg%bmb_g        = 0.0_prec 
            
        end if 
        
        ! ===== Floating ice variables =====

        if (npts_flt .gt. 0) then 

            ! ytopo variables 
            reg%H_ice_f      = sum(tpo%now%H_ice,mask=mask_flt)/npts_flt     ! [m]
            reg%V_ice_f      = sum(tpo%now%H_ice,mask=mask_flt)*tpo%par%dx*tpo%par%dy*m3_km3             ! [km^3]
            reg%A_ice_f      = count(tpo%now%H_ice .gt. 0.0 .and. mask_flt)*tpo%par%dx*tpo%par%dy*m2_km2 ! [km^2]

            ! ydyn variables 
            reg%uxy_bar_f    = sum(dyn%now%uxy_bar,mask=mask_flt)/npts_flt      ! [m/a]
            reg%uxy_s_f      = sum(dyn%now%uxy_s,mask=mask_flt)/npts_flt        ! [m/a]
            reg%uxy_b_f      = sum(dyn%now%uxy_b,mask=mask_flt)/npts_flt        ! [m/a]

            ! Boundary variables
            reg%z_sl         = sum(bnd%z_sl,mask=mask_grnd)/npts_grnd             ! [m]
            reg%bmb_shlf     = sum(tpo%now%bmb,mask=mask_flt)/npts_flt
            reg%T_shlf       = sum(bnd%T_shlf,mask=mask_flt)/npts_flt
            
        else 

            ! ytopo variables 
            reg%H_ice_f      = 0.0_prec 

            reg%V_ice_f      = 0.0_prec 
            reg%A_ice_f      = 0.0_prec 

            ! ydyn variables 
            reg%uxy_bar_f    = 0.0_prec 
            reg%uxy_s_f      = 0.0_prec 
            reg%uxy_b_f      = 0.0_prec 
            
            ! Boundary variables
            reg%z_sl         = 0.0_prec 
            reg%bmb_shlf     = 0.0_prec 
            reg%T_shlf       = 0.0_prec 
            
        end if 
        
        return 

    end subroutine calc_yregions 

    subroutine write_yreg_init(dom,filename,time_init,units,mask)

        implicit none 

        type(yelmo_class), intent(IN) :: dom 
        character(len=*),  intent(IN) :: filename, units 
        real(prec),        intent(IN) :: time_init
        logical,           intent(IN) :: mask(:,:) 

        ! Initialize netcdf file and dimensions
        call nc_create(filename)
        call nc_write_dim(filename,"xc",        x=dom%grd%xc*1e-3,      units="kilometers")
        call nc_write_dim(filename,"yc",        x=dom%grd%yc*1e-3,      units="kilometers")
        call nc_write_dim(filename,"zeta",      x=dom%par%zeta_aa,      units="1")
        call nc_write_dim(filename,"time",      x=time_init,dx=1.0_prec,nx=1,units=trim(units),unlimited=.TRUE.)
        
        ! Static information
        call nc_write(filename,"mask", mask,  units="1",long_name="Region mask",dim1="xc",dim2="yc")
        
        return

    end subroutine write_yreg_init 

    subroutine write_yreg_step(reg,filename,time)

        implicit none 
        
        type(yregions_class), intent(IN) :: reg
        character(len=*),     intent(IN) :: filename
        real(prec),           intent(IN) :: time

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

        ! ===== Total ice variables =====

        call nc_write(filename,"H_ice",reg%H_ice,units="m",long_name="Mean ice thickness", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"z_srf",reg%z_srf,units="m",long_name="Mean surface elevation", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"dHicedt",reg%dHicedt,units="m/a",long_name="Mean rate ice thickness change", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"H_ice_max",reg%H_ice_max,units="m/a",long_name="Max ice thickness", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"dzsrfdt",reg%dzsrfdt,units="m/a",long_name="Mean rate surface elevation change", &
                      dim1="time",start=[n],ncid=ncid)
        
        call nc_write(filename,"V_ice",reg%V_ice*1e-6,units="1e6 km^3",long_name="Ice volume", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"A_ice",reg%A_ice*1e-6,units="1e6 km^2",long_name="Ice area", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"dVicedt",reg%dVicedt,units="km^3/a",long_name="Rate volume change", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"fwf",reg%fwf,units="Sv",long_name="Rate volume change", &
                      dim1="time",start=[n],ncid=ncid)

        call nc_write(filename,"uxy_bar",reg%uxy_bar,units="m/a",long_name="Mean depth-ave velocity", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"uxy_s",reg%uxy_s,units="m/a",long_name="Mean surface velocity", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"uxy_b",reg%uxy_b,units="m/a",long_name="Mean basal velocity", &
                      dim1="time",start=[n],ncid=ncid)

        call nc_write(filename,"z_bed",reg%z_bed,units="m",long_name="Mean bedrock elevation", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"smb",reg%smb,units="m/a",long_name="Mean surface mass balance", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"T_srf",reg%T_srf,units="K",long_name="Mean surface temperature", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"bmb",reg%bmb,units="m/a",long_name="Mean total basal mass balance", &
                      dim1="time",start=[n],ncid=ncid)
        


        ! ===== Grounded ice variables =====

        call nc_write(filename,"H_ice_g",reg%H_ice_g,units="m",long_name="Mean ice thickness (grounded)", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"z_srf_g",reg%z_srf_g,units="m",long_name="Mean surface elevation (grounded)", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"V_ice_g",reg%V_ice_g*1e-6,units="1e6 km^3",long_name="Ice volume (grounded)", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"A_ice_g",reg%A_ice_g*1e-6,units="1e6 km^2",long_name="Ice area (grounded)", &
                      dim1="time",start=[n],ncid=ncid)

        call nc_write(filename,"uxy_bar_g",reg%uxy_bar_g,units="m/a",long_name="Mean depth-ave velocity (grounded)", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"uxy_s_g",reg%uxy_s_g,units="m/a",long_name="Mean surface velocity (grounded)", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"uxy_b_g",reg%uxy_b_g,units="m/a",long_name="Mean basal velocity (grounded)", &
                      dim1="time",start=[n],ncid=ncid)
        
        call nc_write(filename,"f_pmp",reg%f_pmp,units="1",long_name="Temperate fraction (grounded)", &
                      dim1="time",start=[n],ncid=ncid)
        
        call nc_write(filename,"H_w",reg%H_w,units="m",long_name="Mean basal water thickness (grounded)", &
                      dim1="time",start=[n],ncid=ncid)
        
        call nc_write(filename,"bmb_g",reg%bmb_g,units="m/a",long_name="Mean basal mass balance (grounded)", &
                      dim1="time",start=[n],ncid=ncid)
        
        ! ===== Floating ice variables =====

        call nc_write(filename,"H_ice_f",reg%H_ice_f,units="m",long_name="Mean ice thickness (floating)", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"V_ice_f",reg%V_ice_f*1e-6,units="1e6 km^3",long_name="Ice volume (floating)", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"A_ice_f",reg%A_ice_f*1e-6,units="1e6 km^2",long_name="Ice area (floating)", &
                      dim1="time",start=[n],ncid=ncid)

        call nc_write(filename,"uxy_bar_f",reg%uxy_bar_f,units="m/a",long_name="Mean depth-ave velocity (floating)", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"uxy_s_f",reg%uxy_s_f,units="m/a",long_name="Mean surface velocity (floating)", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"uxy_b_f",reg%uxy_b_f,units="m/a",long_name="Mean basal velocity (floating)", &
                      dim1="time",start=[n],ncid=ncid)
        
        call nc_write(filename,"z_sl",reg%z_sl,units="m",long_name="Mean sea level (floating)", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"bmb_shlf",reg%bmb_shlf,units="m/a",long_name="Mean basal mass balance (floating)", &
                      dim1="time",start=[n],ncid=ncid)
        call nc_write(filename,"T_shlf",reg%T_shlf,units="K",long_name="Mean marine shelf temperature (floating)", &
                      dim1="time",start=[n],ncid=ncid)

        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine write_yreg_step

    subroutine calc_write_yreg_step(dom,filename,time,mask)

        implicit none 
        
        type(yelmo_class), intent(IN) :: dom 
        character(len=*),  intent(IN) :: filename
        real(prec),        intent(IN) :: time
        logical,           intent(IN) :: mask(:,:) 

        ! Local variables
        integer    :: ncid, n
        real(prec) :: time_prev 

        type(yregions_class) :: reg
        
        ! 1. Calculate the regional values 
        call calc_yregions(reg,dom%tpo,dom%dyn,dom%thrm,dom%mat,dom%bnd,mask) 

        ! 2. Write the region to the specified file 
        call write_yreg_step(reg,filename,time)

        return 

    end subroutine calc_write_yreg_step

end module yelmo_regions
