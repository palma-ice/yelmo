
module yelmo_data

    use nml 
    use ncio 
    use yelmo_defs 
    use yelmo_tools, only : adjust_topography_gradients
    use topography 
    
    implicit none
    
    private
    public :: ydata_alloc, ydata_dealloc
    public :: ydata_par_load, ydata_load 
    public :: ydata_compare

contains

    subroutine ydata_compare(dta,tpo,dyn,mat,thrm,bnd,domain)

        implicit none 

        type(ydata_class),  intent(INOUT) :: dta
        type(ytopo_class),  intent(IN)    :: tpo 
        type(ydyn_class),   intent(IN)    :: dyn 
        type(ymat_class),   intent(IN)    :: mat
        type(ytherm_class), intent(IN)    :: thrm
        type(ybound_class), intent(IN)    :: bnd
        character(len=*),   intent(IN)    :: domain 

        ! Local variables 
        integer :: q, q1, nx, ny 
        logical,  allocatable :: mask_region(:,:)
        logical,  allocatable :: mask(:,:)
        real(dp), allocatable :: tmp(:,:) 
        real(dp), allocatable :: tmp1(:,:) 
        
        real(dp), parameter :: tol = 1e-3 

        nx = size(tpo%now%H_ice,1)
        ny = size(tpo%now%H_ice,2)
        
        allocate(mask(nx,ny))
        allocate(mask_region(nx,ny))

        ! ======================================================
        ! Calculate errors

        dta%pd%err_H_ice   = tpo%now%H_ice - dta%pd%H_ice 
        dta%pd%err_z_srf   = tpo%now%z_srf - dta%pd%z_srf 
        dta%pd%err_uxy_s   = dyn%now%uxy_s - dta%pd%uxy_s 
        
        ! Isochronal layer error 
        dta%pd%err_depth_iso = mv

        do q = 1, dta%par%pd_age_n_iso
            ! Loop over observed isochronal layer depths 

            do q1 = 1, mat%par%n_iso
                ! Loop over isochronal layer depths in Yelmo 

                if (abs(mat%par%age_iso(q1)-dta%pd%age_iso(q)) .lt. tol) then 
                    ! Isochronal layer in data matches this one 

                    where(dta%pd%depth_iso(:,:,q) .ne. mv) 
                        dta%pd%err_depth_iso(:,:,q) = mat%now%depth_iso(:,:,q1) - dta%pd%depth_iso(:,:,q)
                    elsewhere 
                        dta%pd%err_depth_iso(:,:,q) = mv 
                    end where 

                end if

            end do 

        end do 

        ! ======================================================
        ! Whole ice sheet error metrics (rmse)

        ! Calculate region over which to calculate metrics 

        ! By default, calculate everywhere 
        mask_region = .TRUE. 

        ! For Greenland, limit to continental Greenland where data is defined currently
        if (trim(domain) .eq. "Greenland") then 
            where (bnd%regions .ne. 1.3) mask_region = .FALSE. 
        end if 

        allocate(tmp(nx,ny))
        allocate(tmp1(nx,ny))

        ! == rmse[Ice thickness] ===================
        
        tmp = tpo%now%H_ice-dta%pd%H_ice
        
        ! Define mask over which to perform comparison with data 
        mask = (tpo%now%H_ice .ne. 0.0 .or. dta%pd%H_ice .ne. 0.0) .and. mask_region 
        
        if (count(mask) .gt. 0) then 
            dta%pd%rmse_H = sqrt(sum(tmp**2,mask=mask)/count(mask))
        else 
            dta%pd%rmse_H = mv 
        end if 

        if (dta%pd%rmse_H .eq. 0.0_wp) dta%pd%rmse_H = mv 

        ! == rmse[Surface elevation] ===================

        tmp = dta%pd%err_z_srf
        
        mask = tmp .ne. 0.0 .and. mask_region 
         
        if (count(mask) .gt. 0) then 
            dta%pd%rmse_zsrf = sqrt(sum(tmp**2,mask=mask)/count(mask))
        else 
            dta%pd%rmse_zsrf = mv 
        end if 

        if (dta%pd%rmse_zsrf .eq. 0.0_wp) dta%pd%rmse_zsrf = mv 
        
        ! == rmse[Surface velocity] ===================
        tmp = dta%pd%err_uxy_s

        mask = tmp .ne. 0.0 .and. mask_region 
         
        if (count(mask) .gt. 0) then
            dta%pd%rmse_uxy = sqrt(sum(tmp**2,mask=mask)/count(mask))
        else
            dta%pd%rmse_uxy = mv
        end if 

        if (dta%pd%rmse_uxy .eq. 0.0_wp) dta%pd%rmse_uxy = mv 
        
        ! == rmse[log(Surface velocity)] ===================
        tmp = dta%pd%uxy_s 
        where(dta%pd%uxy_s .gt. 0.0) tmp = log(tmp)
        tmp1 = dyn%now%uxy_s 
        where(dyn%now%uxy_s .gt. 0.0) tmp1 = log(tmp1)
        
        mask = (tmp .ne. 0.0 .or. tmp1 .ne. 0.0) .and. mask_region 
        
        if (count(mask) .gt. 0) then 
            dta%pd%rmse_loguxy = sqrt(sum((tmp1-tmp)**2,mask=mask)/count(mask))
        else
            dta%pd%rmse_loguxy = mv
        end if 
        
        if (dta%pd%rmse_loguxy .eq. 0.0_wp) dta%pd%rmse_loguxy = mv 
        
        ! == rmse[isochronal layer depth] ============

        do q1 = 1, dta%par%pd_age_n_iso
            mask = dta%pd%err_depth_iso(:,:,q1) .ne. mv
            if (count(mask) .gt. 0) then 
                dta%pd%rmse_iso(q1) = sqrt( sum(dta%pd%err_depth_iso(:,:,q1)**2,mask=mask) / count(mask) )
            else 
                dta%pd%rmse_iso(q1) = mv 
            end if 
        end do 

        return 

    end subroutine ydata_compare

    subroutine ydata_load(dta,bnd,par_path,grad_lim_zb,dx,boundaries)

        implicit none 

        type(ydata_class),  intent(INOUT) :: dta 
        type(ybound_class), intent(IN)    :: bnd 
        character(len=*),   intent(IN)    :: par_path 
        real(wp),           intent(IN)    :: grad_lim_zb
        real(wp),           intent(IN)    :: dx 
        character(len=*),   intent(IN)    :: boundaries 

        ! Local variables 
        character(len=1028) :: filename 
        character(len=56)   :: nms(5) 
        real(wp)            :: z_bed_f_sd
        real(wp), allocatable :: z_bed_sd(:,:) 
        real(wp), allocatable :: tmp(:,:,:) 

        real(wp), parameter :: z_sl_pd = 0.0_wp     ! [m] Define present day relative sea level as zero

        ! Allocate temporary array for loading monthly data 
        allocate(z_bed_sd(size(dta%pd%H_ice,1),size(dta%pd%H_ice,2)))
        allocate(tmp(size(dta%pd%H_ice,1),size(dta%pd%H_ice,2),12))

        if (dta%par%pd_topo_load) then 
            ! Load present-day data from specified files and fields

            ! =========================================
            ! Load topography data from netcdf file 
            filename = dta%par%pd_topo_path
            nms(1:5) = dta%par%pd_topo_names 

            write(*, *) "Reading Hice at PD..."
            call nc_read(filename,nms(1), dta%pd%H_ice, missing_value=mv)
            write(*, *) "Reading zbed at PD..."
            call nc_read(filename,nms(2), dta%pd%z_bed, missing_value=mv) 
            write(*, *) "Reading mask at PD..."
            call nc_read(filename,nms(5), dta%pd%mask, missing_value=mv) 
            write(*, *) "Done reading mask at PD!"
            
            ! If available read in bedrock standard deviation field
            if (trim(nms(3)) .ne. ""     .and. &
                trim(nms(3)) .ne. "none" .and. &
                trim(nms(3)) .ne. "None") then 

                ! Read in z_bed_sd
                call nc_read(filename,nms(3),z_bed_sd)

                ! Determine scaling factor from yelmo_init_topo parameter choice
                ! Note: reading from "yelmo_init_topo" section is not optimal,
                ! but it is important to keep the choice of z_bed_f_sd consistent
                ! between data loaded for initializing the model, and data loaded
                ! as the reference present day topography. This ensures that the
                ! parameter choice is the same for both cases. Perhaps this can
                ! be improved though.

                call nml_read(par_path,"yelmo_init_topo","z_bed_f_sd",z_bed_f_sd)

                ! Apply scaling to adjust z_bed depending on standard deviation
                ! Use scaling suppied as input argument
                dta%pd%z_bed = dta%pd%z_bed + z_bed_f_sd*z_bed_sd 

            else 
                ! No scaling loaded or applied

                z_bed_f_sd = 0.0_wp
                z_bed_sd   = 0.0_wp 

            end if 
            
            call nc_read(filename,nms(4), dta%pd%z_srf, missing_value=mv)
            
            ! Remove englacial lakes for better comparison with model
            ! Assume sea level is present day level of 0.
            call remove_englacial_lakes(dta%pd%H_ice,dta%pd%z_bed, &
                                        dta%pd%z_srf,dta%pd%z_srf*0.0_wp, &
                                        bnd%c%rho_ice,bnd%c%rho_sw)

            write(*,*) "ydata_load:: removed englacial lakes from PD reference ice thickness."

            ! Clean up field 
            where(dta%pd%H_ice  .lt. 1.0) dta%pd%H_ice = 0.0 

            ! Adjust bedrock and ice thickness for smoothness
            call adjust_topography_gradients(dta%pd%z_bed,dta%pd%H_ice,grad_lim_zb,dx,boundaries)

            ! Artificially delete ice from locations that are not allowed
            where (.not. bnd%ice_allowed) 
                dta%pd%H_ice = 0.0 
                dta%pd%z_srf = max(dta%pd%z_bed,0.0)
            end where 
            
            ! Calculate H_grnd (ice thickness overburden)
            ! (extracted from `calc_H_grnd` in topography)
            dta%pd%H_grnd = dta%pd%H_ice - (bnd%c%rho_sw/bnd%c%rho_ice)*max(z_sl_pd-dta%pd%z_bed,0.0_wp)

            ! Define the mask to be consistent with internal mask_bed calculations
            dta%pd%mask_bed = mask_bed_ocean
            where(dta%pd%H_ice .eq. 0.0 .and. dta%pd%z_srf .gt. 0.0)  dta%pd%mask_bed = mask_bed_land
            where(dta%pd%H_ice .gt. 0.0 .and. dta%pd%H_grnd .lt. 0.0) dta%pd%mask_bed = mask_bed_float
            where(dta%pd%H_ice .gt. 0.0 .and. dta%pd%H_grnd .ge. 0.0) dta%pd%mask_bed = mask_bed_frozen

        end if 

        if (dta%par%pd_tsrf_load) then 
            ! Load present-day data for surface temperature (or near-surface temperature)

            ! =========================================
            ! Load climate data from netcdf file 
            filename = dta%par%pd_tsrf_path
            nms(1)   = dta%par%pd_tsrf_name 

            if (dta%par%pd_tsrf_monthly) then
                ! Monthly data => annual mean 
                call nc_read(filename,nms(1), tmp, missing_value=mv)
                dta%pd%T_srf = sum(tmp,dim=3) / 12.0
            else 
                ! Annual mean 
                call nc_read(filename,nms(1), dta%pd%T_srf, missing_value=mv)
            end if 

            ! Make sure temperatures are in Kelvin 
            if (minval(dta%pd%T_srf,mask=dta%pd%T_srf .ne. mv) .lt. 100.0) then
                ! Probably in Celcius, convert...
                dta%pd%T_srf = dta%pd%T_srf + 273.15 
            end if 

        end if 

        if (dta%par%pd_smb_load) then 
            ! Load present-day data for surface mass balance 

            ! =========================================
            ! Load smb data from netcdf file 
            filename = dta%par%pd_smb_path
            nms(1)   = dta%par%pd_smb_name

            if (dta%par%pd_smb_monthly) then 
                call nc_read(filename,nms(1), tmp, missing_value=mv)
                dta%pd%smb = sum(tmp,dim=3) / 12.0

                ! Convert from mm we / day to m ie / a 
                dta%pd%smb = dta%pd%smb * bnd%c%conv_mmdwe_maie

            else 
                call nc_read(filename,nms(1), dta%pd%smb, missing_value=mv)

                ! Convert from mm we / a to m ie / a 
                dta%pd%smb = dta%pd%smb * bnd%c%conv_mmawe_maie

            end if 

            ! Clean smb to avoid tiny values 
            where (abs(dta%pd%smb) .lt. 1e-3) dta%pd%smb = 0.0
            
        end if 

        if (dta%par%pd_vel_load) then 
            ! Load present-day data for surface velocity 

            ! =========================================
            ! Load vel data from netcdf file 
            filename = dta%par%pd_vel_path
            nms(1:2) = dta%par%pd_vel_names

            call nc_read(filename,nms(1), dta%pd%ux_s, missing_value=mv)
            call nc_read(filename,nms(2), dta%pd%uy_s, missing_value=mv)
            where (dta%pd%ux_s .ne. mv .and. dta%pd%uy_s .ne. mv) &
                dta%pd%uxy_s = sqrt(dta%pd%ux_s**2 + dta%pd%uy_s**2)

            ! Make sure that velocity is zero where no ice exists 
            ! (if data was loaded)
            if (dta%par%pd_topo_load) then 

                where (dta%pd%H_ice .lt. 1.0) 
                    dta%pd%ux_s  = 0.0 
                    dta%pd%uy_s  = 0.0 
                    dta%pd%uxy_s = 0.0 
                end where 

            end if 

        end if 

        if (dta%par%pd_age_load) then 
            ! Load present-day data for isochrones (ice ages) 

            ! =========================================
            ! Load vel data from netcdf file 
            filename = dta%par%pd_age_path
            nms(1:2) = dta%par%pd_age_names

            call nc_read(filename,nms(1), dta%pd%age_iso,   missing_value=mv)
            call nc_read(filename,nms(2), dta%pd%depth_iso, missing_value=mv)

        end if 

        ! Summarize data loading 
        write(*,*) "ydata_load:: range(H_ice):     ",   minval(dta%pd%H_ice),   maxval(dta%pd%H_ice)
        write(*,*) "ydata_load:: range(z_srf):     ",   minval(dta%pd%z_srf),   maxval(dta%pd%z_srf)
        write(*,*) "ydata_load:: range(z_bed):     ",   minval(dta%pd%z_bed),   maxval(dta%pd%z_bed)
        write(*,*) "ydata_load:: range(z_bed_sd):  ",   minval(z_bed_sd),       maxval(z_bed_sd)
        write(*,*) "ydata_load:: range(mask_bed):  ",   minval(dta%pd%mask_bed),maxval(dta%pd%mask_bed)
        write(*,*) "ydata_load:: scaling fac z_bed_f_sd: ", z_bed_f_sd  
        write(*,*) "ydata_load:: range(T_srf):     ",   minval(dta%pd%T_srf),   maxval(dta%pd%T_srf)
        write(*,*) "ydata_load:: range(smb):       ",   minval(dta%pd%smb,dta%pd%smb .ne. mv), &
                                                        maxval(dta%pd%smb,dta%pd%smb .ne. mv)
        write(*,*) "ydata_load:: range(uxy_s):     ",   minval(dta%pd%uxy_s,dta%pd%uxy_s .ne. mv), &
                                                        maxval(dta%pd%uxy_s,dta%pd%uxy_s .ne. mv)
        write(*,*) "ydata_load:: range(age_iso):   ",   minval(dta%pd%age_iso), maxval(dta%pd%age_iso)
        write(*,*) "ydata_load:: range(depth_iso): ",   minval(dta%pd%depth_iso,dta%pd%depth_iso .ne. mv), &
                                                        maxval(dta%pd%depth_iso,dta%pd%depth_iso .ne. mv)
                

        return 

    end subroutine ydata_load

    subroutine ydata_par_load(par,filename,domain,grid_name,init)

        type(ydata_param_class), intent(OUT) :: par
        character(len=*),        intent(IN)  :: filename, domain, grid_name   
        logical, optional,       intent(IN)  :: init 

        ! Local variables
        logical :: init_pars 

        init_pars = .FALSE.
        if (present(init)) init_pars = .TRUE. 
 
        ! Store parameter values in output object
        call nml_read(filename,"yelmo_data","pd_topo_load",    par%pd_topo_load,    init=init_pars)
        call nml_read(filename,"yelmo_data","pd_topo_path",    par%pd_topo_path,    init=init_pars)
        call nml_read(filename,"yelmo_data","pd_topo_names",   par%pd_topo_names,   init=init_pars)
        call nml_read(filename,"yelmo_data","pd_tsrf_load",    par%pd_tsrf_load,    init=init_pars)
        call nml_read(filename,"yelmo_data","pd_tsrf_path",    par%pd_tsrf_path,    init=init_pars)
        call nml_read(filename,"yelmo_data","pd_tsrf_name",    par%pd_tsrf_name,    init=init_pars)
        call nml_read(filename,"yelmo_data","pd_tsrf_monthly", par%pd_tsrf_monthly, init=init_pars)
        call nml_read(filename,"yelmo_data","pd_smb_load",     par%pd_smb_load,     init=init_pars)
        call nml_read(filename,"yelmo_data","pd_smb_path",     par%pd_smb_path,     init=init_pars)
        call nml_read(filename,"yelmo_data","pd_smb_name",     par%pd_smb_name,     init=init_pars)
        call nml_read(filename,"yelmo_data","pd_smb_monthly",  par%pd_smb_monthly,  init=init_pars)
        call nml_read(filename,"yelmo_data","pd_vel_load",     par%pd_vel_load,     init=init_pars)
        call nml_read(filename,"yelmo_data","pd_vel_path",     par%pd_vel_path,     init=init_pars)
        call nml_read(filename,"yelmo_data","pd_vel_names",    par%pd_vel_names,    init=init_pars)
        call nml_read(filename,"yelmo_data","pd_age_load",     par%pd_age_load,     init=init_pars)
        call nml_read(filename,"yelmo_data","pd_age_path",     par%pd_age_path,     init=init_pars)
        call nml_read(filename,"yelmo_data","pd_age_names",    par%pd_age_names,    init=init_pars)
        
        ! Subsitute domain/grid_name
        call yelmo_parse_path(par%pd_topo_path,domain,grid_name)
        call yelmo_parse_path(par%pd_tsrf_path,domain,grid_name)
        call yelmo_parse_path(par%pd_smb_path, domain,grid_name)
        call yelmo_parse_path(par%pd_vel_path, domain,grid_name)
        call yelmo_parse_path(par%pd_age_path, domain,grid_name)

        ! Internal parameters 
        par%domain = trim(domain) 
        
        ! Get number of isochrone layers to load 
        if (par%pd_age_load) then 
            par%pd_age_n_iso = nc_size(par%pd_age_path,par%pd_age_names(1))
        else 
            par%pd_age_n_iso = 1        ! To avoid allocation errors  
        end if 

        return

    end subroutine ydata_par_load

    subroutine ydata_alloc(pd,nx,ny,nz,n_iso)

        implicit none 

        type(ydata_pd_class) :: pd 
        integer :: nx, ny, nz, n_iso  

        call ydata_dealloc(pd)
        
        allocate(pd%H_ice(nx,ny))
        allocate(pd%z_srf(nx,ny))
        allocate(pd%z_bed(nx,ny))
        allocate(pd%H_grnd(nx,ny))
        allocate(pd%mask(nx,ny))       
        allocate(pd%mask_bed(nx,ny))
        
        allocate(pd%T_srf(nx,ny))
        allocate(pd%smb(nx,ny))
        
        allocate(pd%age_iso(n_iso))
        allocate(pd%depth_iso(nx,ny,n_iso))

        allocate(pd%ux_s(nx,ny))
        allocate(pd%uy_s(nx,ny))
        allocate(pd%uxy_s(nx,ny))
        
        allocate(pd%err_H_ice(nx,ny))
        allocate(pd%err_z_srf(nx,ny))
        allocate(pd%err_z_bed(nx,ny))
        allocate(pd%err_uxy_s(nx,ny))
        allocate(pd%err_depth_iso(nx,ny,n_iso))

        allocate(pd%rmse_iso(n_iso))
        
        pd%H_ice         = 0.0 
        pd%z_srf         = 0.0 
        pd%z_bed         = 0.0 
        pd%H_grnd        = 0.0 
        pd%mask_bed      = 0

        pd%T_srf         = 0.0 
        pd%smb           = 0.0 

        pd%age_iso       = 0.0 
        pd%depth_iso     = 0.0 

        pd%ux_s          = 0.0 
        pd%uy_s          = 0.0 
        pd%uxy_s         = 0.0 
        
        pd%err_H_ice     = 0.0 
        pd%err_z_srf     = 0.0 
        pd%err_z_bed     = 0.0 
        pd%err_uxy_s     = 0.0 
        pd%err_depth_iso = 0.0 
        
        pd%rmse_iso      = mv 

        ! Assign some other initial values here too
        pd%rmse_H        = mv  
        pd%rmse_zsrf     = mv 
        pd%rmse_uxy      = mv 
        pd%rmse_loguxy   = mv 

        return 

    end subroutine ydata_alloc

    subroutine ydata_dealloc(pd)

        implicit none 

        type(ydata_pd_class) :: pd

        if (allocated(pd%H_ice))            deallocate(pd%H_ice)
        if (allocated(pd%z_srf))            deallocate(pd%z_srf)
        if (allocated(pd%z_bed))            deallocate(pd%z_bed)
        if (allocated(pd%H_grnd))           deallocate(pd%H_grnd)
        if (allocated(pd%mask))             deallocate(pd%mask)
        if (allocated(pd%mask_bed))         deallocate(pd%mask_bed)
        
        if (allocated(pd%T_srf))            deallocate(pd%T_srf)
        if (allocated(pd%smb))              deallocate(pd%smb)
        
        if (allocated(pd%depth_iso))        deallocate(pd%depth_iso)
        
        if (allocated(pd%ux_s))             deallocate(pd%ux_s)
        if (allocated(pd%uy_s))             deallocate(pd%uy_s)
        if (allocated(pd%uxy_s))            deallocate(pd%uxy_s)
        
        if (allocated(pd%err_H_ice))        deallocate(pd%err_H_ice)
        if (allocated(pd%err_z_srf))        deallocate(pd%err_z_srf)
        if (allocated(pd%err_z_bed))        deallocate(pd%err_z_bed)
        if (allocated(pd%err_uxy_s))        deallocate(pd%err_uxy_s)
        if (allocated(pd%err_depth_iso))    deallocate(pd%err_depth_iso)
        
        if (allocated(pd%rmse_iso))         deallocate(pd%rmse_iso)
        
        return 

    end subroutine ydata_dealloc

end module yelmo_data
