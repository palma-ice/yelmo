
module yelmo_boundaries

    use nml 
    use ncio 
    use yelmo_defs 

    implicit none
    
    private
    public :: ybound_define_physical_constants
    public :: ybound_load_masks
    public :: ybound_define_ice_allowed
    public :: ybound_alloc, ybound_dealloc
    
contains

    subroutine ybound_define_physical_constants(c,phys_const,domain,grid_name)

        implicit none

        type(ybound_const_class), intent(OUT) :: c 
        character(len=*), intent(IN) :: phys_const 
        character(len=*), intent(IN) :: domain 
        character(len=*), intent(IN) :: grid_name 
        
        ! Local variables 
        logical :: init_pars 

        select case(trim(phys_const))

            case("Earth")
                c%sec_year   = 31536000.0       ! [s/a]  365*24*3600
                c%g          = 9.81             ! [m/s^2] Gravitational accel.
                c%T0         = 273.15           ! [K] Reference freezing temperature 
                c%rho_ice    =  910.0           ! [kg/m^3] Density ice            
                c%rho_w      = 1000.0           ! [kg/m^3] Density water          
                c%rho_sw     = 1028.0           ! [kg/m^3] Density seawater       
                c%rho_a      = 3300.0           ! [kg m-3] Density asthenosphere
                c%rho_rock   = 2000.0           ! [kg m-3] Density bedrock (mantle/lithosphere)
                c%L_ice      = 333500.0         ! [J kg-1] Latent heat of fusion for ice/water
                c%T_pmp_beta = 9.8e-8           ! [K Pa^-1] Greve and Blatter (2009)   

            case("EISMINT","EISMINT1","EISMINT2")

                c%sec_year   = 31556926.0       ! [s/a] EISMINT value
                c%g          = 9.81             ! [m s-2] Gravitational accel.
                c%T0         = 273.15           ! [K] Reference freezing temperature 
                c%rho_ice    =  910.0           ! [kg m-3] Density ice            
                c%rho_w      = 1000.0           ! [kg m-3] Density water          
                c%rho_sw     = 1028.0           ! [kg m-3] Density seawater       
                c%rho_a      = 3300.0           ! [kg m-3] Density asthenosphere
                c%rho_rock   = 2000.0           ! [kg m-3] Density bedrock (mantle/lithosphere)
                c%L_ice      = 333500.0         ! [J kg-1] Latent heat of fusion for ice/water
                c%T_pmp_beta = 9.7e-8           ! [K Pa^-1] EISMINT2 value (beta1 = 8.66e-4 [K m^-1])
            
            case("MISMIP3D")

                c%sec_year   = 31556926.0       ! [s/a] EISMINT value
                c%g          = 9.81             ! [m s-2] Gravitational accel.
                c%T0         = 273.15           ! [K] Reference freezing temperature 
                c%rho_ice    =  900.0           ! [kg m-3] Density ice            
                c%rho_w      = 1000.0           ! [kg m-3] Density water          
                c%rho_sw     = 1028.0           ! [kg m-3] Density seawater       
                c%rho_a      = 3300.0           ! [kg m-3] Density asthenosphere
                c%rho_rock   = 2000.0           ! [kg m-3] Density bedrock (mantle/lithosphere)
                c%L_ice      = 333500.0         ! [J kg-1] Latent heat of fusion for ice/water
                c%T_pmp_beta = 9.7e-8           ! [K Pa^-1] EISMINT2 value (beta1 = 8.66e-4 [K m^-1])

            case("TROUGH")

                c%sec_year   = 31556926.0       ! [s/a] EISMINT value
                c%g          = 9.81             ! [m s-2] Gravitational accel.
                c%T0         = 273.15           ! [K] Reference freezing temperature 
                c%rho_ice    =  918.0           ! [kg m-3] Density ice            
                c%rho_w      = 1000.0           ! [kg m-3] Density water          
                c%rho_sw     = 1028.0           ! [kg m-3] Density seawater       
                c%rho_a      = 3300.0           ! [kg m-3] Density asthenosphere
                c%rho_rock   = 2000.0           ! [kg m-3] Density bedrock (mantle/lithosphere)
                c%L_ice      = 333500.0         ! [J kg-1] Latent heat of fusion for ice/water
                c%T_pmp_beta = 9.7e-8           ! [K Pa^-1] EISMINT2 value (beta1 = 8.66e-4 [K m^-1])

            case DEFAULT
                ! Load parameter values from parameter file

                init_pars = .TRUE. 
                
                call nml_read(phys_const,"yconst","sec_year",    c%sec_year,   init=init_pars)
                call nml_read(phys_const,"yconst","g",           c%g,          init=init_pars)
                call nml_read(phys_const,"yconst","T0",          c%T0,         init=init_pars)
                call nml_read(phys_const,"yconst","rho_ice",     c%rho_ice,    init=init_pars)
                call nml_read(phys_const,"yconst","rho_w",       c%rho_w,      init=init_pars)
                call nml_read(phys_const,"yconst","rho_sw",      c%rho_sw,     init=init_pars)
                call nml_read(phys_const,"yconst","rho_a",       c%rho_a,      init=init_pars)
                call nml_read(phys_const,"yconst","rho_rock",    c%rho_rock,   init=init_pars)
                call nml_read(phys_const,"yconst","L_ice",       c%L_ice,      init=init_pars)
                call nml_read(phys_const,"yconst","T_pmp_beta",  c%T_pmp_beta, init=init_pars)

        end select


        ! Define conversion factors too

        c%conv_we_ie          = c%rho_w/c%rho_ice
        c%conv_mmdwe_maie     = 1e-3*365*c%conv_we_ie
        c%conv_mmawe_maie     = 1e-3*c%conv_we_ie
        
        c%conv_m3_Gt          = c%rho_ice *1e-12                ! [kg/m3] * [Gigaton/1e12kg]
        c%conv_km3_Gt         = (1e9) * c%conv_m3_Gt            ! [1e9m^3/km^3]
        c%conv_millionkm3_Gt  = (1e6) * (1e9) *c%conv_m3_Gt     ! [1e6km3/1] * [1e9m^3/km^3] * conv
        
        c%area_seasurf        = 3.618e8                         ! [km^2]
        c%conv_km3_sle        = (1e-3) / 394.7                  ! [m/mm] / [km^3 to raise ocean by 1mm] => m sle, see https://sealevel.info/conversion_factors.html

        if (yelmo_log) then
            write(*,*) "yelmo:: loaded physical constants: "
            write(*,*) "domain:    ", trim(domain)
            write(*,*) "grid_name: ", trim(grid_name)
            write(*,*) "    sec_year   = ", c%sec_year 
            write(*,*) "    g          = ", c%g 
            write(*,*) "    T0         = ", c%T0 
            write(*,*) "    rho_ice    = ", c%rho_ice 
            write(*,*) "    rho_w      = ", c%rho_w 
            write(*,*) "    rho_sw     = ", c%rho_sw 
            write(*,*) "    rho_a      = ", c%rho_a 
            write(*,*) "    rho_rock   = ", c%rho_rock 
            write(*,*) "    L_ice      = ", c%L_ice 
            write(*,*) "    T_pmp_beta = ", c%T_pmp_beta 
            
        end if 

        return

    end subroutine ybound_define_physical_constants
    
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

                bnd%ice_allowed = .FALSE. 
                where (bnd%regions .eq. 1.3)  bnd%ice_allowed = .TRUE.      ! Main Greenland region
                where (bnd%regions .eq. 1.11) bnd%ice_allowed = .TRUE.      ! Ellesmere Island
                where (bnd%regions .eq. 1.0)  bnd%ice_allowed = .TRUE.      ! Open ocean (included some connections between 1.3 and 1.11)
                
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
        allocate(now%fmb_shlf(nx,ny))
        allocate(now%T_shlf(nx,ny))
        allocate(now%Q_geo(nx,ny))

        allocate(now%enh_srf(nx,ny)) 

        allocate(now%basins(nx,ny))
        allocate(now%basin_mask(nx,ny))
        allocate(now%regions(nx,ny))
        allocate(now%region_mask(nx,ny))
        
        allocate(now%ice_allowed(nx,ny))
        allocate(now%calv_mask(nx,ny))
        
        allocate(now%H_ice_ref(nx,ny))
        allocate(now%z_bed_ref(nx,ny))
        
        allocate(now%z_bed_corr(nx,ny))

        now%z_bed       = 0.0_prec 
        now%z_bed_sd    = 0.0_prec
        now%z_sl        = 0.0_prec 
        now%H_sed       = 0.0_prec 
        now%smb         = 0.0_prec 
        now%T_srf       = 0.0_prec 
        now%bmb_shlf    = 0.0_prec 
        now%fmb_shlf    = 0.0_prec 
        now%T_shlf      = 0.0_prec 
        now%Q_geo       = 0.0_prec 

        now%enh_srf     = 1.0_prec 

        now%basins      = 0.0_prec 
        now%basin_mask  = 0.0_prec 
        now%regions     = 0.0_prec 
        now%region_mask = 0.0_prec 
        
        now%ice_allowed = .TRUE.  ! By default allow ice everywhere 
        now%calv_mask   = .FALSE. ! By default no, no calving mask 

        now%H_ice_ref   = 0.0_prec 
        now%z_bed_ref   = 0.0_prec

        now%z_bed_corr  = 0.0_prec

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
        if (allocated(now%fmb_shlf))    deallocate(now%fmb_shlf)
        if (allocated(now%T_shlf))      deallocate(now%T_shlf)
        if (allocated(now%Q_geo))       deallocate(now%Q_geo)
        
        if (allocated(now%enh_srf))     deallocate(now%enh_srf)
        
        if (allocated(now%basins))      deallocate(now%basins)
        if (allocated(now%basin_mask))  deallocate(now%basin_mask)
        if (allocated(now%regions))     deallocate(now%regions)
        if (allocated(now%region_mask)) deallocate(now%region_mask)
        
        if (allocated(now%ice_allowed)) deallocate(now%ice_allowed)
        if (allocated(now%calv_mask))   deallocate(now%calv_mask)
        
        if (allocated(now%H_ice_ref))   deallocate(now%H_ice_ref) 
        if (allocated(now%z_bed_ref))   deallocate(now%z_bed_ref)

        if (allocated(now%z_bed_corr))  deallocate(now%z_bed_corr)

        return 

    end subroutine ybound_dealloc

end module yelmo_boundaries
