
module yelmo_data

    use nml 
    use ncio 
    use yelmo_defs 

    implicit none
    
    private
    public :: ydata_alloc, ydata_dealloc
    public :: ydata_par_load, ydata_load 
    public :: ydata_compare

contains

    subroutine ydata_compare(dta,tpo,dyn,thrm,bnd)

        implicit none 

        type(ydata_class),  intent(INOUT) :: dta
        type(ytopo_class),  intent(IN)    :: tpo 
        type(ydyn_class),   intent(IN)    :: dyn 
        type(ytherm_class), intent(IN)    :: thrm
        type(ybound_class), intent(IN)    :: bnd
        
        ! Calculate errors
        dta%pd%err_H_ice   = tpo%now%H_ice - dta%pd%H_ice 
        dta%pd%err_z_srf   = tpo%now%z_srf - dta%pd%z_srf 
        dta%pd%err_uxy_s   = dyn%now%uxy_s - dta%pd%uxy_s 
        
        return 

    end subroutine ydata_compare 

    subroutine ydata_load(dta,ice_allowed)

        implicit none 

        type(ydata_class), intent(INOUT) :: dta 
        logical,           intent(IN)    :: ice_allowed(:,:) 

        ! Local variables 
        character(len=1028) :: filename 
        character(len=56)   :: nms(4) 
        real(prec), allocatable :: tmp(:,:,:) 

        ! Allocate temporary array for loading monthly data 
        allocate(tmp(size(dta%pd%H_ice,1),size(dta%pd%H_ice,2),12))

        if (dta%par%pd_topo_load) then 
            ! Load present-day data from specified files and fields

            ! =========================================
            ! Load topography data from netcdf file 
            filename = dta%par%pd_topo_path
            nms(1:3) = dta%par%pd_topo_names 

            call nc_read(filename,nms(1), dta%pd%H_ice, missing_value=mv)
            call nc_read(filename,nms(2), dta%pd%z_srf, missing_value=mv)
            call nc_read(filename,nms(3), dta%pd%z_bed, missing_value=mv) 

            ! Clean up field 
            where(dta%pd%H_ice  .lt. 1.0) dta%pd%H_ice = 0.0 

            ! Artificially delete ice from locations that are not allowed
            where (.not. ice_allowed) 
                dta%pd%H_ice = 0.0 
                dta%pd%z_srf = max(dta%pd%z_bed,0.0)
            end where 
            
            ! Calculate H_grnd (ice thickness overburden)
            ! (extracted from `calc_H_grnd` in yelmo_topography)
            dta%pd%H_grnd = dta%pd%H_ice - (rho_sw/rho_ice)*(0.0_prec-dta%pd%z_bed)

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
                dta%pd%smb = dta%pd%smb * conv_mmdwe_maie

            else 
                call nc_read(filename,nms(1), dta%pd%smb, missing_value=mv)

                ! Convert from mm we / a to m ie / a 
                dta%pd%smb = dta%pd%smb * conv_mmawe_maie

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

        ! Summarize data loading 
        write(*,*) "ydata_load:: range(H_ice): ", minval(dta%pd%H_ice),   maxval(dta%pd%H_ice)
        write(*,*) "ydata_load:: range(z_srf): ", minval(dta%pd%z_srf),   maxval(dta%pd%z_srf)
        write(*,*) "ydata_load:: range(z_bed): ", minval(dta%pd%z_bed),   maxval(dta%pd%z_bed)
        write(*,*) "ydata_load:: range(T_srf): ", minval(dta%pd%T_srf), maxval(dta%pd%T_srf)
        write(*,*) "ydata_load:: range(smb):   ", minval(dta%pd%smb,dta%pd%smb .ne. mv), &
                                                  maxval(dta%pd%smb,dta%pd%smb .ne. mv)
        write(*,*) "ydata_load:: range(uxy_s): ", minval(dta%pd%uxy_s,dta%pd%uxy_s .ne. mv), &
                                                  maxval(dta%pd%uxy_s,dta%pd%uxy_s .ne. mv)
            

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
        
        ! Subsitute domain/grid_name
        call yelmo_parse_path(par%pd_topo_path,domain,grid_name)
        call yelmo_parse_path(par%pd_tsrf_path,domain,grid_name)
        call yelmo_parse_path(par%pd_smb_path, domain,grid_name)
        call yelmo_parse_path(par%pd_vel_path, domain,grid_name)

        ! Internal parameters 
        par%domain = trim(domain) 
         
        return

    end subroutine ydata_par_load

    subroutine ydata_alloc(pd,nx,ny,nz)

        implicit none 

        type(ydata_pd_class) :: pd 
        integer :: nx, ny, nz  

        call ydata_dealloc(pd)
        
        allocate(pd%H_ice(nx,ny))
        allocate(pd%z_srf(nx,ny))
        allocate(pd%z_bed(nx,ny))
        allocate(pd%H_grnd(nx,ny))
        
        allocate(pd%T_srf(nx,ny))
        allocate(pd%smb(nx,ny))
        
        allocate(pd%ux_s(nx,ny))
        allocate(pd%uy_s(nx,ny))
        allocate(pd%uxy_s(nx,ny))
        
        allocate(pd%err_H_ice(nx,ny))
        allocate(pd%err_z_srf(nx,ny))
        allocate(pd%err_z_bed(nx,ny))
        
        allocate(pd%err_uxy_s(nx,ny))
        
        pd%H_ice        = 0.0 
        pd%z_srf        = 0.0 
        pd%z_bed        = 0.0 
        pd%H_grnd       = 0.0 
        
        pd%T_srf        = 0.0 
        pd%smb          = 0.0 

        pd%ux_s         = 0.0 
        pd%uy_s         = 0.0 
        pd%uxy_s        = 0.0 
        
        pd%err_H_ice    = 0.0 
        pd%err_z_srf    = 0.0 
        pd%err_z_bed    = 0.0 
        
        pd%err_uxy_s = 0.0 
        
        return 
    end subroutine ydata_alloc 

    subroutine ydata_dealloc(pd)

        implicit none 

        type(ydata_pd_class) :: pd

        if (allocated(pd%H_ice))        deallocate(pd%H_ice)
        if (allocated(pd%z_srf))        deallocate(pd%z_srf)
        if (allocated(pd%z_bed))        deallocate(pd%z_bed)
        if (allocated(pd%H_grnd))       deallocate(pd%H_grnd)
        
        if (allocated(pd%T_srf))        deallocate(pd%T_srf)
        if (allocated(pd%smb))          deallocate(pd%smb)
        
        if (allocated(pd%ux_s))         deallocate(pd%ux_s)
        if (allocated(pd%uy_s))         deallocate(pd%uy_s)
        if (allocated(pd%uxy_s))        deallocate(pd%uxy_s)
        
        if (allocated(pd%err_H_ice))    deallocate(pd%err_H_ice)
        if (allocated(pd%err_z_srf))    deallocate(pd%err_z_srf)
        if (allocated(pd%err_z_bed))    deallocate(pd%err_z_bed)
        
        if (allocated(pd%err_uxy_s))    deallocate(pd%err_uxy_s)
        
        return 

    end subroutine ydata_dealloc 

end module yelmo_data
