
module yelmo_data

    use nml 
    use ncio 
    use yelmo_defs 

    implicit none
    
    private
    public :: ydata_alloc, ydata_dealloc
    public :: ydata_par_load, ydata_load 

contains

    subroutine ydata_load(dta)

        implicit none 

        type(ydata_class), intent(INOUT) :: dta 

        ! Local variables 
        character(len=1028) :: filename 
        character(len=56)   :: nms(4) 
        real(prec), allocatable :: tmp(:,:,:) 

        ! Allocate temporary array for loading monthly data 
        allocate(tmp(size(dta%pd%H_ice,1),size(dta%pd%H_ice,2),12))

        if (dta%par%pd_load_data) then 
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

            ! =========================================
            ! Load climate data from netcdf file 
            filename = dta%par%pd_clim_path
            nms      = dta%par%pd_clim_names 

            ! == t2m, pr ==
            
            if (dta%par%pd_clim_monthly) then
                
                call nc_read(filename,nms(2), tmp, missing_value=mv)
                dta%pd%t2m_ann = sum(tmp,dim=3) / 12.0

                if (trim(dta%par%domain) == "Antarctica") then 
                    ! Southern Hemisphere summer
                    dta%pd%t2m_sum = (tmp(:,:,12) + tmp(:,:,1) + tmp(:,:,2)) / 3.0 
                else 
                    ! Northern Hemisphere summer 
                    dta%pd%t2m_sum = sum(tmp(:,:,6:8),dim=3) / 3.0
                end if

                call nc_read(filename,nms(3), tmp, missing_value=mv)
                dta%pd%pr_ann = sum(tmp,dim=3) / 12.0 
                
                if (trim(nms(3)) == "sf" .and. trim(nms(4)) == "rf") then 
                    ! Load rainfall too, to get sum 
                    call nc_read(filename,nms(4), tmp, missing_value=mv)
                    dta%pd%pr_ann = dta%pd%pr_ann + sum(tmp,dim=3) / 12.0
                end if 

                ! Convert from mm we / day to m ie / a 
                dta%pd%pr_ann = dta%pd%pr_ann * conv_mmdwe_maie

            else 
                call nc_read(filename,nms(2), dta%pd%t2m_ann, missing_value=mv)
                dta%pd%t2m_sum = mv

                call nc_read(filename,nms(3), dta%pd%pr_ann, missing_value=mv)
                if (trim(nms(3)) == "sf" .and. trim(nms(4)) == "rf") then
                    call nc_read(filename,nms(3), tmp(:,:,1), missing_value=mv)
                    dta%pd%pr_ann = dta%pd%pr_ann + tmp(:,:,1)
                end if

                ! Convert from mm we / a to m ie / a 
                dta%pd%pr_ann = dta%pd%pr_ann * conv_mmdwe_maie

            end if 

            ! Make sure temperatures are in Kelvin 
            if (minval(dta%pd%t2m_ann,mask=dta%pd%t2m_ann .ne. mv) .lt. 100.0) then
                ! Probably in Celcius, convert...
                dta%pd%t2m_ann = dta%pd%t2m_ann + 273.15 
                dta%pd%t2m_sum = dta%pd%t2m_sum + 273.15
            end if 

            
            ! =========================================
            ! Load smb data from netcdf file 
            filename = dta%par%pd_smb_path
            nms(1)   = dta%par%pd_smb_name

            ! == smb == 

            if (dta%par%pd_smb_monthly) then 
                call nc_read(filename,nms(1), tmp, missing_value=mv)
                dta%pd%smb_ann = sum(tmp,dim=3) / 12.0

                ! Convert from mm we / day to m ie / a 
                dta%pd%smb_ann = dta%pd%smb_ann * conv_mmdwe_maie

            else 
                call nc_read(filename,nms(1), dta%pd%smb_ann, missing_value=mv)

                ! Convert from mm we / a to m ie / a 
                dta%pd%smb_ann = dta%pd%smb_ann * conv_mmawe_maie

            end if 

            ! Clean smb to avoid tiny values 
            where (abs(dta%pd%smb_ann) .lt. 1e-3) dta%pd%smb_ann = 0.0
            
            ! Limit smb to ice and land points 
            !where (dta%pd%z_srf .eq. 0.0) dta%pd%smb_ann = mv 

            ! =========================================
            ! Load vel data from netcdf file 
            filename = dta%par%pd_vel_path
            nms(1:2) = dta%par%pd_vel_names

            call nc_read(filename,nms(1), dta%pd%ux_s, missing_value=mv)
            call nc_read(filename,nms(2), dta%pd%uy_s, missing_value=mv)
            where (dta%pd%ux_s .ne. mv .and. dta%pd%uy_s .ne. mv) &
                dta%pd%uxy_s = sqrt(dta%pd%ux_s**2 + dta%pd%uy_s**2)

            ! Make sure that velocity is zero where no ice exists 
            where (dta%pd%H_ice .lt. 1.0) 
                dta%pd%ux_s  = 0.0 
                dta%pd%uy_s  = 0.0 
                dta%pd%uxy_s = 0.0 
            end where 

            ! Summarize data loading 
            write(*,*) "ydata_load:: range(H_ice):   ", minval(dta%pd%H_ice),   maxval(dta%pd%H_ice)
            write(*,*) "ydata_load:: range(z_srf):   ", minval(dta%pd%z_srf),   maxval(dta%pd%z_srf)
            write(*,*) "ydata_load:: range(z_bed):   ", minval(dta%pd%z_bed),   maxval(dta%pd%z_bed)
            write(*,*) "ydata_load:: range(t2m_ann): ", minval(dta%pd%t2m_ann), maxval(dta%pd%t2m_ann)
            write(*,*) "ydata_load:: range(t2m_sum): ", minval(dta%pd%t2m_sum), maxval(dta%pd%t2m_sum)
            write(*,*) "ydata_load:: range(pr_ann):  ", minval(dta%pd%pr_ann),  maxval(dta%pd%pr_ann)
            write(*,*) "ydata_load:: range(smb_ann): ", minval(dta%pd%smb_ann,dta%pd%smb_ann .ne. mv), &
                                                        maxval(dta%pd%smb_ann,dta%pd%smb_ann .ne. mv)
            write(*,*) "ydata_load:: range(uxy_s):   ", minval(dta%pd%uxy_s,dta%pd%uxy_s .ne. mv), &
                                                        maxval(dta%pd%uxy_s,dta%pd%uxy_s .ne. mv)
            
        end if 

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
        call nml_read(filename,"yelmo_data","pd_load_data",    par%pd_load_data,    init=init_pars)
        call nml_read(filename,"yelmo_data","pd_topo_path",    par%pd_topo_path,    init=init_pars)
        call nml_read(filename,"yelmo_data","pd_topo_names",   par%pd_topo_names,   init=init_pars)
        call nml_read(filename,"yelmo_data","pd_clim_path",    par%pd_clim_path,    init=init_pars)
        call nml_read(filename,"yelmo_data","pd_clim_names",   par%pd_clim_names,   init=init_pars)
        call nml_read(filename,"yelmo_data","pd_clim_monthly", par%pd_clim_monthly, init=init_pars)
        call nml_read(filename,"yelmo_data","pd_smb_path",     par%pd_smb_path,     init=init_pars)
        call nml_read(filename,"yelmo_data","pd_smb_name",     par%pd_smb_name,     init=init_pars)
        call nml_read(filename,"yelmo_data","pd_smb_monthly",  par%pd_smb_monthly,  init=init_pars)
        call nml_read(filename,"yelmo_data","pd_vel_path",     par%pd_vel_path,     init=init_pars)
        call nml_read(filename,"yelmo_data","pd_vel_names",    par%pd_vel_names,    init=init_pars)
        
        ! Subsitute domain/grid_name
        call yelmo_parse_path(par%pd_topo_path,domain,grid_name)
        call yelmo_parse_path(par%pd_clim_path,domain,grid_name)
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
        
        allocate(pd%t2m_ann(nx,ny))
        allocate(pd%t2m_sum(nx,ny))
        allocate(pd%pr_ann(nx,ny))
        allocate(pd%smb_ann(nx,ny))
        
        allocate(pd%ux_s(nx,ny))
        allocate(pd%uy_s(nx,ny))
        allocate(pd%uxy_s(nx,ny))
        
        allocate(pd%err_H_ice(nx,ny))
        allocate(pd%err_z_srf(nx,ny))
        allocate(pd%err_z_bed(nx,ny))
        
        allocate(pd%err_uxy_s(nx,ny))
        
        pd%H_ice       = 0.0 
        pd%z_srf       = 0.0 
        pd%z_bed       = 0.0 
        
        pd%t2m_ann     = 0.0 
        pd%t2m_sum     = 0.0 
        pd%pr_ann      = 0.0 
        pd%smb_ann     = 0.0

        pd%ux_s      = 0.0 
        pd%uy_s      = 0.0 
        pd%uxy_s     = 0.0 
        
        pd%err_H_ice   = 0.0 
        pd%err_z_srf   = 0.0 
        pd%err_z_bed   = 0.0 
        
        pd%err_uxy_s = 0.0 
        
        return 
    end subroutine ydata_alloc 

    subroutine ydata_dealloc(pd)

        implicit none 

        type(ydata_pd_class) :: pd

        if (allocated(pd%H_ice)) deallocate(pd%H_ice)
        if (allocated(pd%z_srf)) deallocate(pd%z_srf)
        if (allocated(pd%z_bed)) deallocate(pd%z_bed)
        
        if (allocated(pd%t2m_ann)) deallocate(pd%t2m_ann)
        if (allocated(pd%t2m_sum)) deallocate(pd%t2m_sum)
        if (allocated(pd%pr_ann))  deallocate(pd%pr_ann)
        if (allocated(pd%smb_ann)) deallocate(pd%smb_ann)
        
        if (allocated(pd%ux_s))  deallocate(pd%ux_s)
        if (allocated(pd%uy_s))  deallocate(pd%uy_s)
        if (allocated(pd%uxy_s)) deallocate(pd%uxy_s)
        
        if (allocated(pd%err_H_ice)) deallocate(pd%err_H_ice)
        if (allocated(pd%err_z_srf)) deallocate(pd%err_z_srf)
        if (allocated(pd%err_z_bed)) deallocate(pd%err_z_bed)
        
        if (allocated(pd%err_uxy_s)) deallocate(pd%err_uxy_s)
        
        return 

    end subroutine ydata_dealloc 

end module yelmo_data
