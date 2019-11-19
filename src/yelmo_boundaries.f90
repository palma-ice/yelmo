
module yelmo_boundaries

    use nml 
    use ncio 
    use yelmo_defs 

    implicit none
    
    private
    public :: ybound_load_masks
    public :: ybound_define_ice_allowed
    public :: ybound_load_z_bed 
    public :: ybound_alloc, ybound_dealloc
    
    public :: ybound_load_pd 
    
contains

    subroutine ybound_load_masks(bnd,nml_path,nml_group,domain,grid_name)
        ! Load masks for managing regions and basins, etc. 

        implicit none 

        type(ybound_class), intent(INOUT) :: bnd 
        character(len=*), intent(IN)      :: nml_path, nml_group
        character(len=*), intent(IN)      :: domain, grid_name 

        ! Local variables
        logical            :: load_var
        character(len=512) :: filename 
        character(len=56)  :: vnames(2)  

        ! ====================================
        !
        ! basins 
        !
        ! ====================================

        ! Specify default values 
        bnd%basin_mask = 1.0
        bnd%basins     = 1.0

        call nml_read(nml_path,nml_group,"basins_load",load_var)

        if (load_var) then

            call nml_read(nml_path,nml_group, "basins_path",filename)
            call yelmo_parse_path(filename,domain,grid_name)
            
            call nml_read(nml_path,nml_group,"basins_nms",vnames) 
            ! Load basin information from a file 
            call nc_read(filename,vnames(1),bnd%basins)

            if (trim(vnames(2)) .ne. "None") then 
                ! If basins have been extrapolated, also load the original basin extent mask
                call nc_read(filename,vnames(2),bnd%basin_mask)
            end if 

        end if 
        
        ! ====================================
        !
        ! Read in the regions
        !
        ! ====================================

        ! First assign default region values
        bnd%region_mask = 1.0 
        select case(trim(domain))
            case("North")
                bnd%regions = bnd%index_north
            case("Antarctica")
                bnd%regions = bnd%index_south
            case("Greenland")
                bnd%regions = bnd%index_grl
            case DEFAULT 
                ! Assign a default value everywhere 
                bnd%regions = 1.0  

        end select

        call nml_read(nml_path,nml_group,"regions_load",load_var)
          
        if (load_var) then
            
            call nml_read(nml_path,nml_group, "regions_path",filename)
            call yelmo_parse_path(filename,domain,grid_name)
            
            ! Load region information from a file 
            call nml_read(nml_path,nml_group,"regions_nms",vnames)
            call nc_read(filename,vnames(1),bnd%regions)

            if (trim(vnames(2)) .ne. "None") then 
                ! If regions have been extrapolated, also load the original region extent mask
                call nc_read(filename,vnames(2),bnd%region_mask)
            end if 
            
        end if 
        
        write(*,*) "ybound_load_masks:: range(basins):  ", minval(bnd%basins),  maxval(bnd%basins)
        write(*,*) "ybound_load_masks:: range(regions): ", minval(bnd%regions), maxval(bnd%regions)

        return 

    end subroutine ybound_load_masks

    subroutine ybound_define_ice_allowed(bnd,domain)
        ! Update mask of where ice is allowed to be greater than zero 

        implicit none 

        type(ybound_class), intent(INOUT) :: bnd
        character(len=*),   intent(IN)    :: domain    

        ! Local variables 
        integer :: i, nx, ny 

        nx = size(bnd%ice_allowed,1)
        ny = size(bnd%ice_allowed,2)

        ! Initially set ice_allowed to true everywhere 
        bnd%ice_allowed = .TRUE. 

        ! Also set calv_mask false everywhere (no imposed calving front)
        bnd%calv_mask   = .FALSE. 
        
        ! Determine allowed regions based on domain
        select case(trim(domain))

            case ("North") 
                ! Allow ice everywhere except the open ocean 

                where (bnd%regions .eq. 1.0) bnd%ice_allowed = .FALSE. 
                bnd%ice_allowed(1,:)  = .FALSE. 
                bnd%ice_allowed(nx,:) = .FALSE. 
                bnd%ice_allowed(:,1)  = .FALSE. 
                bnd%ice_allowed(:,ny) = .FALSE. 

            case ("Eurasia") 
                ! Allow ice only in the Eurasia domain (1.2*)

                where (bnd%regions .lt. 1.2 .or. bnd%regions .gt. 1.29) bnd%ice_allowed = .FALSE. 
                bnd%ice_allowed(1,:)  = .FALSE. 
                bnd%ice_allowed(nx,:) = .FALSE. 
                bnd%ice_allowed(:,1)  = .FALSE. 
                bnd%ice_allowed(:,ny) = .FALSE. 

            case ("Greenland") 

                where (bnd%regions .ne. 1.3) bnd%ice_allowed = .FALSE. 
                bnd%ice_allowed(1,:)  = .FALSE. 
                bnd%ice_allowed(nx,:) = .FALSE. 
                bnd%ice_allowed(:,1)  = .FALSE. 
                bnd%ice_allowed(:,ny) = .FALSE. 

            case ("Antarctica") 

                where (bnd%regions .eq. 2.0) bnd%ice_allowed = .FALSE. 
                bnd%ice_allowed(1,:)  = .FALSE. 
                bnd%ice_allowed(nx,:) = .FALSE. 
                bnd%ice_allowed(:,1)  = .FALSE. 
                bnd%ice_allowed(:,ny) = .FALSE. 

            
            case ("EISMINT")

                ! Ice can grow everywhere, except borders 
                bnd%ice_allowed       = .TRUE. 
                bnd%ice_allowed(1,:)  = .FALSE. 
                bnd%ice_allowed(nx,:) = .FALSE. 
                bnd%ice_allowed(:,1)  = .FALSE. 
                bnd%ice_allowed(:,ny) = .FALSE. 

            case ("MISMIP") 

                ! Ice can grow everywhere, except farthest x-border  
                bnd%ice_allowed       = .TRUE. 
                bnd%ice_allowed(nx,:) = .FALSE. 
                
            case DEFAULT 
                ! Let ice grow everywhere for unknown domain (ice_allowed can always be modified later)

                bnd%ice_allowed = .TRUE. 

        end select

                
        return 

    end subroutine ybound_define_ice_allowed

    subroutine ybound_load_z_bed(bnd,nml_path,nml_group,domain,grid_name)
        ! Load the bedrock elevation boundary field

        implicit none

        type(ybound_class), intent(INOUT) :: bnd
        character(len=*),   intent(IN)    :: nml_path, nml_group
        character(len=*),   intent(IN)    :: domain, grid_name 

        ! Local variables
        logical            :: load_var
        character(len=512) :: filename
        character(len=56)  :: vname

        ! == z_bed =====
        call nml_read(nml_path,nml_group,"z_bed_load",load_var)

        if (load_var) then 
            call nml_read(nml_path,nml_group,"z_bed_path",filename)
            call yelmo_parse_path(filename,domain,grid_name)
            call nml_read(nml_path,nml_group,"z_bed_nm",  vname)
            call nc_read(filename,vname,bnd%z_bed)

            vname = ""
            call nml_read(nml_path,nml_group,"z_bed_sd_nm",  vname)
            
            if (trim(vname) .ne. "" .and. trim(vname) .ne. "none" .and. trim(vname) .ne. "None") then 
                call nc_read(filename,vname,bnd%z_bed_sd)
            else
                bnd%z_bed_sd = 0.0 
            end if 

        else 
            bnd%z_bed    = 0.0 
            bnd%z_bed_sd = 0.0  
        end if 
        
        write(*,*) "ybound_load_z_bed:: range(z_bed):     ", minval(bnd%z_bed),    maxval(bnd%z_bed)
        write(*,*) "ybound_load_z_bed:: range(z_bed_sd):  ", minval(bnd%z_bed_sd), maxval(bnd%z_bed_sd)

        return 

    end subroutine ybound_load_z_bed 
    
    subroutine ybound_alloc(now,nx,ny)

        implicit none 

        type(ybound_class) :: now 
        integer :: nx, ny 

        call ybound_dealloc(now)

        allocate(now%z_bed(nx,ny))
        allocate(now%z_bed_sd(nx,ny))
        allocate(now%z_sl(nx,ny))
        allocate(now%H_sed(nx,ny))
        allocate(now%smb(nx,ny))
        allocate(now%T_srf(nx,ny))
        allocate(now%bmb_shlf(nx,ny))
        allocate(now%T_shlf(nx,ny))
        allocate(now%Q_geo(nx,ny))

        allocate(now%basins(nx,ny))
        allocate(now%basin_mask(nx,ny))
        allocate(now%regions(nx,ny))
        allocate(now%region_mask(nx,ny))
        
        allocate(now%ice_allowed(nx,ny))
        allocate(now%calv_mask(nx,ny))
        
        allocate(now%H_ice_ref(nx,ny))
!mmr
        allocate(now%z_bed_ref(nx,ny)) !mmr 
!mmr
        
        now%z_bed       = 0.0_prec 
        now%z_bed_sd    = 0.0_prec
        now%z_sl        = 0.0_prec 
        now%H_sed       = 0.0_prec 
        now%smb         = 0.0_prec 
        now%T_srf       = 0.0_prec 
        now%bmb_shlf    = 0.0_prec 
        now%T_shlf      = 0.0_prec 
        now%Q_geo       = 0.0_prec 

        now%basins      = 0.0_prec 
        now%basin_mask  = 0.0_prec 
        now%regions     = 0.0_prec 
        now%region_mask = 0.0_prec 
        
        now%ice_allowed = .TRUE.  ! By default allow ice everywhere 
        now%calv_mask   = .FALSE. ! By default no, no calving mask 

        now%H_ice_ref   = 0.0_prec 
!mmr
        now%z_bed_ref = 0.0_prec !mmr
!mr
        return 
    end subroutine ybound_alloc 

    subroutine ybound_dealloc(now)

        implicit none 

        type(ybound_class) :: now

        if (allocated(now%z_bed))       deallocate(now%z_bed)
        if (allocated(now%z_bed_sd))    deallocate(now%z_bed_sd)
        if (allocated(now%z_sl))        deallocate(now%z_sl)
        if (allocated(now%H_sed))       deallocate(now%H_sed)
        if (allocated(now%smb))         deallocate(now%smb)
        if (allocated(now%T_srf))       deallocate(now%T_srf)
        if (allocated(now%bmb_shlf))    deallocate(now%bmb_shlf)
        if (allocated(now%T_shlf))      deallocate(now%T_shlf)
        if (allocated(now%Q_geo))       deallocate(now%Q_geo)
        
        if (allocated(now%basins))      deallocate(now%basins)
        if (allocated(now%basin_mask))  deallocate(now%basin_mask)
        if (allocated(now%regions))     deallocate(now%regions)
        if (allocated(now%region_mask)) deallocate(now%region_mask)
        
        if (allocated(now%ice_allowed)) deallocate(now%ice_allowed)
        if (allocated(now%calv_mask))   deallocate(now%calv_mask)
        
        if (allocated(now%H_ice_ref))   deallocate(now%H_ice_ref)
!mmr 
       if (allocated(now%z_bed_ref))   deallocate(now%z_bed_ref) !mmr
!mmr        
        return 

    end subroutine ybound_dealloc 

! ajr: unused...
! Note: some of this functionality is now in yelmo_data. Careful thought
! needed about which routines are relevant and in which module they should appear...
    
    subroutine ybound_load_pd(bnd,nml_path,nml_group,domain,grid_name)
        ! Load all boundary fields for the present day (except z_bed, handled via ybound_load_zbed) 

        ! ajr: not used yet, needs testing...

        implicit none

        type(ybound_class), intent(INOUT) :: bnd 
        character(len=*),   intent(IN)    :: nml_path, nml_group 
        character(len=*),   intent(IN)    :: domain, grid_name

        ! Local variables
        logical            :: load_var
        character(len=512) :: filename
        character(len=56)  :: vname 

        ! == z_sl =====

        bnd%z_sl = 0.0 

        ! == H_sed =====
        call nml_read(nml_path,nml_group,"H_sed_load",load_var)
        if (load_var) then 
            call nml_read(nml_path,nml_group,"H_sed_path",filename)
            call yelmo_parse_path(filename,domain,grid_name)
            call nml_read(nml_path,nml_group,"H_sed_nm",  vname)
            call nc_read(filename,vname,bnd%H_sed)
        else 
            bnd%H_sed = mv 
        end if 
        
        ! == smb =====
        call nml_read(nml_path,nml_group,"smb_load",load_var)
        if (load_var) then 
            call nml_read(nml_path,nml_group,"smb_path",filename)
            call yelmo_parse_path(filename,domain,grid_name)
            call nml_read(nml_path,nml_group,"smb_nm",  vname)
            call nc_read(filename,vname,bnd%smb)
        else 
            bnd%smb = mv 
        end if 

        ! == T_srf =====
        call nml_read(nml_path,nml_group,"T_srf_load",load_var)
        if (load_var) then 
            call nml_read(nml_path,nml_group,"T_srf_path",filename)
            call yelmo_parse_path(filename,domain,grid_name)
            call nml_read(nml_path,nml_group,"T_srf_nm",  vname)
            call nc_read(filename,vname,bnd%T_srf)
        else 
            bnd%T_srf = mv 
        end if 

        ! == bmb_shlf =====
        call nml_read(nml_path,nml_group,"bmb_shlf_load",load_var)
        if (load_var) then 
            call nml_read(nml_path,nml_group,"bmb_shlf_path",filename)
            call yelmo_parse_path(filename,domain,grid_name)
            call nml_read(nml_path,nml_group,"bmb_shlf_nm",  vname)
            call nc_read(filename,vname,bnd%bmb_shlf)
        else 
            bnd%bmb_shlf = mv 
        end if 

        ! == T_shlf =====
        call nml_read(nml_path,nml_group,"T_shlf_load",load_var)
        if (load_var) then 
            call nml_read(nml_path,nml_group,"T_shlf_path",filename)
            call yelmo_parse_path(filename,domain,grid_name)
            call nml_read(nml_path,nml_group,"T_shlf_nm",  vname)
            call nc_read(filename,vname,bnd%T_shlf)
        else 
            bnd%T_shlf = mv 
        end if 

        ! == Q_geo =====
        call nml_read(nml_path,nml_group,"Q_geo_load",load_var)
        if (load_var) then 
            call nml_read(nml_path,nml_group,"Q_geo_path",filename)
            call yelmo_parse_path(filename,domain,grid_name)
            call nml_read(nml_path,nml_group,"Q_geo_nm",  vname)
            call nc_read(filename,vname,bnd%Q_geo)
        else 
            bnd%Q_geo = mv 
        end if 

        return 

    end subroutine ybound_load_pd 

end module yelmo_boundaries
