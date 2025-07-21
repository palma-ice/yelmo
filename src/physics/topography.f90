module topography 

    use yelmo_defs, only : wp, dp, io_unit_err, pi
    use yelmo_tools, only : get_neighbor_indices
    use subgrid, only : calc_subgrid_array, calc_subgrid_array_cell

    implicit none 

    ! Key for matching bed types given by mask_bed 
    integer, parameter :: mask_bed_ocean   = 0 
    integer, parameter :: mask_bed_land    = 1
    integer, parameter :: mask_bed_frozen  = 2
    integer, parameter :: mask_bed_stream  = 3
    integer, parameter :: mask_bed_grline  = 4
    integer, parameter :: mask_bed_float   = 5
    integer, parameter :: mask_bed_island  = 6
    integer, parameter :: mask_bed_partial = 7

    private  

    public :: gen_mask_bed

    public :: calc_ice_fraction
    public :: calc_ice_fraction_new
    public :: calc_ice_front

    public :: calc_z_srf
    public :: calc_z_srf_max
    public :: calc_z_srf_gl_subgrid_area
    public :: calc_H_eff
    public :: calc_H_grnd
    public :: calc_H_af
    public :: calc_f_grnd_subgrid_area_aa
    public :: calc_f_grnd_subgrid_area
    public :: calc_f_grnd_subgrid_linear
    public :: calc_f_grnd_pinning_points
    public :: remove_englacial_lakes
    public :: calc_distance_to_ice_margin
    public :: calc_distance_to_grounding_line
    public :: calc_grounding_line_zone
    public :: calc_bmb_total
    public :: calc_fmb_total
    
    ! ajr: these routines are slow, do not use...
    !public :: distance_to_grline
    !public :: distance_to_margin
    
    public :: determine_grounded_fractions

    ! Integers
    public :: mask_bed_ocean  
    public :: mask_bed_land  
    public :: mask_bed_frozen
    public :: mask_bed_stream
    public :: mask_bed_grline
    public :: mask_bed_float 
    public :: mask_bed_island
    
contains 

    elemental subroutine gen_mask_bed(mask,f_ice,f_pmp,f_grnd,mask_grz)
        ! Generate an output mask for model conditions at bed
        ! based on input masks 
        ! 0: ocean, 1: land, 2: sia, 3: streams, grline: 4, floating: 5, islands: 6
        ! 7: partially-covered ice cell.

        implicit none 

        integer,  intent(OUT) :: mask 
        real(wp), intent(IN)  :: f_ice, f_pmp, f_grnd
        integer,  intent(IN)  :: mask_grz

        if (mask_grz .eq. 0) then
            ! Grounding line

            mask = mask_bed_grline

        else if (f_ice .eq. 0.0) then 
            ! Ice-free points 

            if (f_grnd .gt. 0.0) then
                ! Ice-free land

                mask = mask_bed_land

            else
                ! Ice-free ocean

                mask = mask_bed_ocean

            end if 

        else if (f_ice .gt. 0.0 .and. f_ice .lt. 1.0) then 
            ! Partially ice-covered points 

            mask = mask_bed_partial

        else
            ! Fully ice-covered points 

            if (f_grnd .gt. 0.0) then
                ! Grounded ice-covered points 

                if (f_pmp .gt. 0.5) then 
                    ! Temperate points

                    mask = mask_bed_stream 

                else
                    ! Frozen points 

                    mask = mask_bed_frozen 

                end if 

            else
                ! Floating ice-covered points 

                mask = mask_bed_float

            end if 

        end if 

        return 

    end subroutine gen_mask_bed

    subroutine find_connected_mask(mask,mask_ref,mask_now,boundaries)
        ! Brute-force routine to find all points 
        ! that touch or are connected with other points in a mask.
        ! Here use to find any floating ice points that are
        ! not connected to grounded ice or ice-free land. 

        ! AJR: TO DO, this routine is not ready!!!

        implicit none 

        logical, intent(INOUT) :: mask(:,:)         ! Connected points of interest
        logical, intent(IN)    :: mask_ref(:,:)     ! Points to be connected to
        logical, intent(IN)    :: mask_now(:,:)     ! Points of interest
        character(len=*), intent(IN) :: boundaries
    
        ! Local variables 
        integer :: i, j, q, nx, ny
        integer :: im1, ip1, jm1, jp1 
        integer :: n_unfilled 

        logical, allocatable :: mask0(:,:) 

        integer, parameter :: qmax = 1000 

        nx = size(mask,1)
        ny = size(mask,2) 

        ! Allocate local mask object for diagnosing points of interest
        allocate(mask0(nx,ny))

        ! Initially assume all points are unconnected
        mask = .FALSE. 

        ! Set purely land points to zero 
        where (mask_ref .and. mask_now) mask = .TRUE. 

        ! Iteratively fill in open-ocean points that are found
        ! next to known open-ocean points
        do q = 1, qmax 

            n_unfilled = 0 

            do j = 1, ny 
            do i = 1, nx 

                ! Get neighbor indices
                call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

                if (mask_now(i,j)) then 
                    ! This is a point of interest
                    ! Define any neighbor points of interest as connected

                    if (mask_now(im1,j)) then 
                        mask(im1,j) = .TRUE.
                        n_unfilled = n_unfilled + 1
                    end if 

                    if (mask_now(ip1,j)) then 
                        mask(ip1,j) = .TRUE. 
                        n_unfilled = n_unfilled + 1
                    end if

                    if (mask_now(i,jm1)) then 
                        mask(i,jm1) = .TRUE.
                        n_unfilled = n_unfilled + 1
                    end if 

                    if (mask_now(i,jp1)) then 
                        mask(i,jp1) = .TRUE.
                        n_unfilled = n_unfilled + 1
                    end if 

                end if
                    
            end do 
            end do  

            !write(*,*) q, n_unfilled, count(mask .eq. -1) 

            ! Exit loop if no more open-ocean points are found 
            if (n_unfilled .eq. 0) exit 

        end do 

        return 

    end subroutine find_connected_mask
    
    subroutine calc_ice_fraction_new(f_ice,H_ice,z_bed,z_sl,rho_ice,rho_sw,boundaries,flt_subgrid)
        ! Determine the area fraction of a grid cell
        ! that is ice-covered. Assume that marginal points
        ! have equal thickness to inland neighbors 

        ! Note: routine works well, but apparently not as stable as routine below
        ! for GRL-16KM initmip tests. Reasons are not really clear (ajr, 2022-02-24)

        implicit none 

        real(wp), intent(OUT) :: f_ice(:,:)             ! [--] Ice covered fraction (aa-nodes)
        real(wp), intent(IN)  :: H_ice(:,:)             ! [m] Ice thickness on standard grid (aa-nodes)
        real(wp), intent(IN)  :: z_bed(:,:)             ! [m] Bedrock elevation
        real(wp), intent(IN)  :: z_sl(:,:)              ! [m] Sea-level elevation
        real(wp), intent(IN)  :: rho_ice 
        real(wp), intent(IN)  :: rho_sw 
        character(len=*), intent(IN) :: boundaries
        logical, optional     :: flt_subgrid            ! Option to allow fractions for floating ice margins             
        
        ! Local variables 
        integer  :: i, j, nx, ny
        integer  :: im1, ip1, jm1, jp1 
        real(wp) :: H_eff  
        logical  :: get_fractional_cover 

        real(wp) :: H_neighb(4)
        logical  :: cf_neighb(4) 
        logical  :: mask(4) 
        real(wp), allocatable :: H_grnd(:,:)            ! [m] Thickness until flotation - floating if H_grnd<=0 (aa-nodes)
        logical, allocatable :: mask_cf(:,:) 
                
        nx = size(H_ice,1)
        ny = size(H_ice,2)

        allocate(H_grnd(nx,ny))
        allocate(mask_cf(nx,ny)) 

        ! By default, fractional cover will be determined
        get_fractional_cover = .TRUE. 
        if (present(flt_subgrid)) get_fractional_cover = flt_subgrid 

        if (get_fractional_cover) then 
            ! Calculate H_grnd, without accounting for fractional ice cover
            call calc_H_grnd(H_grnd,H_ice,f_ice,z_bed,z_sl,rho_ice,rho_sw,use_f_ice=.FALSE.)
        end if 

        ! Initialize mask_cf to False 
        mask_cf = .FALSE. 

        ! Initialize f_ice to binary values first 
        where(H_ice .gt. 0.0_wp)
            f_ice = 1.0_wp 
        elsewhere
            f_ice = 0.0_wp 
        end where

        if (get_fractional_cover) then 
            ! For floating ice-covered points with ice-free neighbors (ie, at the floating calving margin),
            ! determine the fraction of grid point that should be ice covered. 

            ! First determine which points are at the floating calving front.

            do j = 1, ny
            do i = 1, nx 

                ! Get neighbor indices
                call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
    
                ! Determine calving front mask 
                if (H_ice(i,j) .gt. 0.0 .and. H_grnd(i,j) .lt. 0.0) then 
                    ! Current point is floating 

                    ! If any neighbor is ice-free ocean (ie, is grounded below sea level),
                    ! then this point is a calving front point.
                    mask_cf(i,j) = .FALSE. 
                    if (H_ice(im1,j) .eq. 0.0 .and. H_grnd(im1,j) .lt. 0.0) mask_cf(i,j) = .TRUE.
                    if (H_ice(ip1,j) .eq. 0.0 .and. H_grnd(ip1,j) .lt. 0.0) mask_cf(i,j) = .TRUE.
                    if (H_ice(i,jm1) .eq. 0.0 .and. H_grnd(i,jm1) .lt. 0.0) mask_cf(i,j) = .TRUE.
                    if (H_ice(i,jp1) .eq. 0.0 .and. H_grnd(i,jp1) .lt. 0.0) mask_cf(i,j) = .TRUE.

                end if 

            end do 
            end do

            ! Next calculate the effective ice thickness of points at the 
            ! calving front and determine f_ice. 

            do j = 1, ny
            do i = 1, nx 

                ! Get neighbor indices
                call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

                if (mask_cf(i,j)) then 
                    ! This is a calving front point 

                    H_neighb  = [H_ice(im1,j),H_ice(ip1,j),H_ice(i,jm1),H_ice(i,jp1)]
                    cf_neighb = [mask_cf(im1,j),mask_cf(ip1,j),mask_cf(i,jm1),mask_cf(i,jp1)]
                    
                    ! Get a mask of available ice-covered neighbors that are not
                    ! at the calving front. 
                    mask = H_neighb .gt. 0.0 .and. (.not. cf_neighb)
                    
                    if ( count(mask) .eq. 0 ) then 
                        ! This point does not have any upstream ice-covered neighbors,
                        ! simply assign a fraction of 1.0.

                        f_ice(i,j) = 1.0_wp 

                    else
                        ! This point does have upstream ice-covered neighbors,
                        ! determine fraction.

                        H_eff = minval(H_neighb,mask=mask)

                        ! If thickest upstream ice is thinner than current
                        ! point, then set effective ice thickness equal
                        ! to that of current point. This ensures
                        ! f_ice below will be bounded 0 < f_ice <= 1
                        if (H_eff .lt. H_ice(i,j)) H_eff = H_ice(i,j) 

                        ! Determine the cell ice fraction
                        ! Note: fraction is determined as a ratio of 
                        ! thicknesses, derived from volume conservation 
                        ! vol = H_ice*dx*dy = H_eff*area_frac 
                        ! f_ice = area_frac / (dx*dy)
                        ! f_ice = H_ice/H_eff 
                        
                        f_ice(i,j) = H_ice(i,j) / H_eff
                        
                    end if 

                end if 

            end do 
            end do 

        end if 

        return 

    end subroutine calc_ice_fraction_new
    
    subroutine calc_ice_fraction(f_ice,H_ice,z_bed,z_sl,rho_ice,rho_sw,boundaries,flt_subgrid)
        ! Determine the area fraction of a grid cell
        ! that is ice-covered. Assume that marginal points
        ! have equal thickness to inland neighbors 

        implicit none 

        real(wp), intent(OUT) :: f_ice(:,:)             ! [--] Ice covered fraction (aa-nodes)
        real(wp), intent(IN)  :: H_ice(:,:)             ! [m] Ice thickness on standard grid (aa-nodes)
        real(wp), intent(IN)  :: z_bed(:,:)             ! [m] Bedrock elevation
        real(wp), intent(IN)  :: z_sl(:,:)              ! [m] Sea-level elevation
        real(wp), intent(IN)  :: rho_ice
        real(wp), intent(IN)  :: rho_sw
        character(len=*), intent(IN) :: boundaries
        logical,  intent(IN), optional :: flt_subgrid   ! Option to allow fractions for floating ice margins             
        
        ! Local variables 
        integer  :: i, j, nx, ny
        integer  :: im1, ip1, jm1, jp1 
        real(wp) :: H_eff  
        logical  :: get_fractional_cover 

        real(wp) :: H_neighb(4)
        real(wp) :: Hg_neighb(4)
        logical  :: float_neighb(4)
        integer  :: n_neighb(4)
        logical  :: mask(4) 
        logical  :: mask_grnd(4)
        integer  :: n_now
        integer, allocatable  :: n_ice(:,:) 
        real(wp), allocatable :: H_grnd(:,:)            ! [m] Thickness until flotation - floating if H_grnd<=0 (aa-nodes)
        
        real(wp), parameter :: H_lim        = 100.0_wp 

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        allocate(n_ice(nx,ny))
        allocate(H_grnd(nx,ny))

        ! By default, fractional cover will be determined
        get_fractional_cover = .TRUE. 
        if (present(flt_subgrid)) get_fractional_cover = flt_subgrid 

        if (get_fractional_cover) then 
            ! Calculate H_grnd, without accounting for fractional ice cover
            call calc_H_grnd(H_grnd,H_ice,f_ice,z_bed,z_sl,rho_ice,rho_sw,use_f_ice=.FALSE.)
        end if 

        ! Initialize f_ice to binary values first 
        where(H_ice .gt. 0.0_wp)
            f_ice = 1.0_wp 
        elsewhere
            f_ice = 0.0_wp 
        end where

        if (get_fractional_cover) then
            ! For ice-covered points with ice-free neighbors (ie, at the floating or grounded margin),
            ! determine the fraction of grid point that should be ice covered. 

            !!$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1, H_neighb,mask)
            do j = 1, ny
            do i = 1, nx 

                ! Get neighbor indices
                call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
                
                ! Count how many neighbors are ice covered  
                H_neighb   = [H_ice(im1,j),H_ice(ip1,j),H_ice(i,jm1),H_ice(i,jp1)]
                mask       = H_neighb .gt. 0.0_wp 
                n_ice(i,j) = count(mask) 

                ! ajr: test.
                ! Count how many neighbors are ice covered and floating
                ! This appears to perform worse.
                ! H_neighb   = [H_ice(im1,j),H_ice(ip1,j),H_ice(i,jm1),H_ice(i,jp1)]
                ! Hg_neighb  = [H_grnd(im1,j),H_grnd(ip1,j),H_grnd(i,jm1),H_grnd(i,jp1)]
                ! mask       = H_neighb .gt. 0.0_wp .and. Hg_neighb .le. 0.0_wp 
                ! n_ice(i,j) = count(mask) 

            end do 
            end do
            !!$omp end parallel do

            ! Determine ice fractional cover for margin points 

            !!$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1, H_neighb,mask,n_now,H_eff)
            do j = 1, ny
            do i = 1, nx 

                ! Get neighbor indices
                call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
                
                if (H_ice(i,j) .gt. 0.0_wp .and. n_ice(i,j) .eq. 0) then 
                    ! First, treat a special case:
                    ! Island point, assume the cell is fully covered
                    ! to ensure it is dynamically active.

                    f_ice(i,j) = 1.0_wp

                else if (H_ice(i,j) .gt. 0.0_wp .and. n_ice(i,j) .lt. 4) then
                    ! This point is ice-covered, but is at the ice margin 

                    ! Get neighbor ice thicknesses and border-ice counts
                    H_neighb = [H_ice(im1,j),H_ice(ip1,j),H_ice(i,jm1),H_ice(i,jp1)]
                    n_neighb = [n_ice(im1,j),n_ice(ip1,j),n_ice(i,jm1),n_ice(i,jp1)]
                    
                    ! Get mask indicating neighbors that both have ice and are surrounded by ice
                    mask     = H_neighb .gt. 0.0_wp .and. n_neighb .eq. 4
                    n_now    = count(mask)

                    if (H_grnd(i,j) .le. 0.0) then
                        ! Floating point 

                        if (n_now .gt. 0) then 
                            ! Get minimum value of upstream neighbors
                            H_eff = minval(H_neighb,mask=mask)
                        else 
                            ! No upstream neighbors available, compare
                            ! ice thickness against minimum allowed
                            ! 'full' ice thickness value H_lim.
                            H_eff = max(H_ice(i,j),H_lim)
                        end if

                    else 
                        ! Grounded point, set H_eff = H_ice following CISM
                        ! (do not allow partially filled cells for grounded ice)

                        H_eff = H_ice(i,j) 

                    end if
                
                    ! Determine the cell ice fraction
                    ! Note: fraction is determined as a ratio of 
                    ! thicknesses, derived from volume conservation 
                    ! vol = H_ice*dx*dy = H_eff*area_frac 
                    ! f_ice = area_frac / (dx*dy)
                    ! f_ice = H_ice/H_eff 
                    ! Note: H_eff == 0.0 probably won't happen, but keep if-statement 
                    ! for safety 

                    if (H_eff .gt. 0.0_wp) then 
                        f_ice(i,j) = min( H_ice(i,j) / H_eff, 1.0_wp ) 
                    else 
                        f_ice(i,j) = 1.0_wp 
                    end if 
                
                end if 

            end do 
            end do 
            !!$omp end parallel do

        end if 

        return 

    end subroutine calc_ice_fraction
    
    subroutine calc_ice_front(mask_frnt,f_ice,f_grnd,z_bed,z_sl,boundaries)
        ! Calculate a mask of ice front points that 
        ! demarcates both the ice front (last ice-covered point)
        ! and the ice-free front (first ice-free point), 
        ! and distinguishes between floating, marine and grounded fronts.

        implicit none 

        integer,  intent(OUT) :: mask_frnt(:,:) 
        real(wp), intent(IN)  :: f_ice(:,:) 
        real(wp), intent(IN)  :: f_grnd(:,:) 
        real(wp), intent(IN)  :: z_bed(:,:)
        real(wp), intent(IN)  :: z_sl(:,:)
        character(len=*), intent(IN) :: boundaries 
        
        ! Local variables 
        integer  :: i, j, nx, ny
        integer  :: im1, ip1, jm1, jp1 
        integer  :: n 
        real(wp) :: f_neighb(4) 

        integer, parameter :: val_ice_free  = -1 
        integer, parameter :: val_flt       = 1
        integer, parameter :: val_marine    = 2
        integer, parameter :: val_grnd      = 3
        
        nx = size(mask_frnt,1) 
        ny = size(mask_frnt,2) 

        ! Initialize mask to zero everywhere to start 
        mask_frnt = 0

        !!$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,n,f_neighb)
        do j = 1, ny
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

            f_neighb = [f_ice(im1,j),f_ice(ip1,j),f_ice(i,jm1),f_ice(i,jp1)]
            n = count(f_neighb .lt. 1.0)

            if (f_ice(i,j) .eq. 1.0 .and. n .gt. 0) then 
                ! This point is an ice front. 

                if (f_grnd(i,j) .gt. 0.0 .and. (z_sl(i,j) .le. z_bed(i,j)) ) then 
                    ! Ice front grounded above sea level 

                    mask_frnt(i,j) = val_grnd 

                else if (f_grnd(i,j) .gt. 0.0) then 
                    ! Ice front grounded below sea level 

                    mask_frnt(i,j) = val_marine 

                else
                    ! Floating ice front 

                    mask_frnt(i,j) = val_flt 

                end if 

                ! Ensure adjacent ice-free points are marked too
                if (f_ice(im1,j) .lt. 1.0) mask_frnt(im1,j) = val_ice_free
                if (f_ice(ip1,j) .lt. 1.0) mask_frnt(ip1,j) = val_ice_free
                if (f_ice(i,jm1) .lt. 1.0) mask_frnt(i,jm1) = val_ice_free
                if (f_ice(i,jp1) .lt. 1.0) mask_frnt(i,jp1) = val_ice_free

            end if 

        end do
        end do 
        !!$omp end parallel do

        return 

    end subroutine calc_ice_front

    elemental subroutine calc_z_srf(z_srf,H_ice,f_ice,H_grnd,z_bed,z_sl,rho_ice,rho_sw)
        ! Calculate surface elevation

        implicit none 

        real(wp), intent(INOUT) :: z_srf
        real(wp), intent(IN)    :: H_ice
        real(wp), intent(IN)    :: f_ice
        real(wp), intent(IN)    :: H_grnd
        real(wp), intent(IN)    :: z_bed
        real(wp), intent(IN)    :: z_sl
        real(wp), intent(IN)    :: rho_ice 
        real(wp), intent(IN)    :: rho_sw
        
        ! Local variables 
        real(wp) :: rho_ice_sw
        real(wp)   :: H_eff 

        rho_ice_sw = rho_ice/rho_sw ! Ratio of density of ice to seawater [--]
        
        ! Get effective ice thickness
        call calc_H_eff(H_eff,H_ice,f_ice,set_frac_zero=.TRUE.)

        ! Update the surface elevation based on z_bed, H_ice and overburden ice thickness 
        if (H_grnd .gt. 0.0) then 
            ! Grounded ice or ice-free land

            z_srf = z_bed + H_eff 

        else
            ! Floating ice or open ocean

            z_srf = z_sl + (1.0-rho_ice_sw)*H_eff

        end if 
        
        return 

    end subroutine calc_z_srf

    elemental subroutine calc_z_srf_max(z_srf,H_ice,f_ice,z_bed,z_sl,rho_ice,rho_sw)
        ! Calculate surface elevation
        ! Adapted from Pattyn (2017), Eq. 1
        
        implicit none 

        real(wp), intent(INOUT) :: z_srf 
        real(wp), intent(IN)    :: H_ice
        real(wp), intent(IN)    :: f_ice
        real(wp), intent(IN)    :: z_bed
        real(wp), intent(IN)    :: z_sl
        real(wp), intent(IN)    :: rho_ice 
        real(wp), intent(IN)    :: rho_sw
        
        ! Local variables
        integer :: i, j, nx, ny 
        real(wp) :: rho_ice_sw
        real(wp) :: H_eff

        rho_ice_sw = rho_ice/rho_sw ! Ratio of density of ice to seawater [--]
        
        ! Get effective ice thickness
        call calc_H_eff(H_eff,H_ice,f_ice,set_frac_zero=.TRUE.)

        ! Initially calculate surface elevation everywhere 
        z_srf = max(z_bed + H_eff, z_sl + (1.0-rho_ice_sw)*H_eff)
        
        return 

    end subroutine calc_z_srf_max

    subroutine calc_z_srf_gl_subgrid_area(z_srf,f_grnd,H_ice,f_ice,z_bed,z_sl,gz_nx,rho_ice,rho_sw,boundaries)
        ! Interpolate variables at grounding line to subgrid level to 
        ! calculate the average z_srf value for the aa-node cell

        implicit none
        
        real(wp), intent(OUT) :: z_srf(:,:)       ! aa-nodes 
        real(wp), intent(IN)  :: f_grnd(:,:)      ! aa-nodes
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: z_bed(:,:)
        real(wp), intent(IN)  :: z_sl(:,:)
        integer,  intent(IN)  :: gz_nx            ! Number of interpolation points per side (nx*nx)
        real(wp), intent(IN)  :: rho_ice 
        real(wp), intent(IN)  :: rho_sw
        character(len=*), intent(IN) :: boundaries
    
        ! Local variables
        integer  :: i, j, nx, ny 
        integer  :: im1, ip1, jm1, jp1
        real(wp) :: H_eff 
        real(wp) :: f_grnd_neighb(4) 
        logical  :: is_grline 

        real(wp), allocatable :: z_srf_int(:,:) 
        real(wp), allocatable :: H_ice_int(:,:)
        real(wp), allocatable :: f_ice_int(:,:)  
        real(wp), allocatable :: z_bed_int(:,:) 
        real(wp), allocatable :: z_sl_int(:,:) 
        
        nx = size(z_srf,1)
        ny = size(z_srf,2) 

        ! Allocate the subgrid arrays 
        allocate(z_srf_int(gz_nx,gz_nx))
        allocate(H_ice_int(gz_nx,gz_nx))
        allocate(f_ice_int(gz_nx,gz_nx))
        allocate(z_bed_int(gz_nx,gz_nx))
        allocate(z_sl_int(gz_nx,gz_nx))
        
        ! ajr: assume f_ice_int=1 everywhere this is used for now. 
        ! Needs to be fixed in the future potentially. 
        f_ice_int = 1.0_wp 

        write(*,*) "calc_z_srf_gl_subgrid_area:: routine not ready for f_ice values. Fix!"
        stop 
        
        ! Calculate the surface elevation based on whole grid values,
        ! except at the grounding line which is treated with subgrid interpolations. 
        !!$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,f_grnd_neighb,is_grline,v1,v2,v3,v4,H_ice_int,z_bed_int,z_sl_int,z_srf_int)
        do j = 1, ny 
        do i = 1, nx

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

            f_grnd_neighb = [f_grnd(im1,j),f_grnd(ip1,j),f_grnd(i,jm1),f_grnd(i,jp1)]

            if (f_grnd(i,j) .eq. 0.0 .and. count(f_grnd_neighb .gt. 0.0).gt.0) then
                is_grline = .TRUE. 
            else if (f_grnd(i,j) .gt. 0.0 .and. count(f_grnd_neighb .eq. 0.0).gt.0) then
                is_grline = .TRUE. 
            else 
                is_grline = .FALSE. 
            end if 

            if (is_grline .and. f_ice(i,j) .eq. 1.0) then 
                ! Only treat grounding line points that are fully ice-covered:  
                ! Perform subgrid calculations 

                ! Calculate subgrid values for this cell
                call calc_subgrid_array(H_ice_int,H_ice,gz_nx,i,j,im1,ip1,jm1,jp1)
                call calc_subgrid_array(z_bed_int,z_bed,gz_nx,i,j,im1,ip1,jm1,jp1)
                call calc_subgrid_array(z_sl_int, z_sl, gz_nx,i,j,im1,ip1,jm1,jp1)
                
                ! Calculate subgrid surface elevations
                call calc_z_srf_max(z_srf_int,H_ice_int,f_ice_int,z_bed_int,z_sl_int,rho_ice,rho_sw)

                ! Calculate full grid z_srf value as the mean of subgrid values 
                z_srf(i,j) = sum(z_srf_int) / real(gz_nx*gz_nx,wp)

            end if 

        end do 
        end do
        !!$omp end parallel do

        return
        
    end subroutine calc_z_srf_gl_subgrid_area

    elemental subroutine calc_H_eff(H_eff,H_ice,f_ice,set_frac_zero)
        ! Calculate ice-thickness, scaled at margins to actual thickness
        ! but as if it covered the whole grid cell.
        
        implicit none

        real(wp), intent(OUT) :: H_eff 
        real(wp), intent(IN)  :: H_ice 
        real(wp), intent(IN)  :: f_ice 
        logical,    intent(IN), optional :: set_frac_zero 

        if (f_ice .gt. 0.0) then 
            H_eff = H_ice / f_ice
        else 
            H_eff = H_ice 
        end if

        if (present(set_frac_zero)) then 
            if (set_frac_zero .and. f_ice .lt. 1.0) H_eff = 0.0_wp 
        end if 

        return

    end subroutine calc_H_eff

    elemental subroutine calc_H_grnd(H_grnd,H_ice,f_ice,z_bed,z_sl,rho_ice,rho_sw,use_f_ice)
        ! Calculate ice thickness overburden, H_grnd
        ! When H_grnd >= 0, grounded, when H_grnd < 0, floating 
        ! Also calculate rate of change for diagnostic related to grounding line 

        ! Note that this will not be identical to thickness above flotation,
        ! since it accounts for height above the water via z_bed too. 

        implicit none 

        real(wp), intent(INOUT) :: H_grnd
        real(wp), intent(IN)    :: H_ice
        real(wp), intent(IN)    :: f_ice
        real(wp), intent(IN)    :: z_bed
        real(wp), intent(IN)    :: z_sl 
        real(wp), intent(IN)    :: rho_ice 
        real(wp), intent(IN)    :: rho_sw
        logical,  intent(IN), optional :: use_f_ice

        ! Local variables   
        real(wp) :: rho_sw_ice 
        real(wp) :: H_eff 
        logical  :: use_f_ice_now 

        ! By default use f_ice to set points with f_ice < 1 to H_eff=0.0 
        use_f_ice_now = .TRUE. 
        if (present(use_f_ice)) use_f_ice_now = use_f_ice 

        rho_sw_ice = rho_sw/rho_ice ! Ratio of density of seawater to ice [--]
        
        ! Get effective ice thickness
        if (use_f_ice_now) then 
            call calc_H_eff(H_eff,H_ice,f_ice,set_frac_zero=.TRUE.)
        else 
            H_eff = H_ice 
        end if 

        ! Calculate new H_grnd (ice thickness overburden)
        !H_grnd = H_eff - rho_sw_ice*max(z_sl-z_bed,0.0_wp)

        ! ajr: testing. This ensures that ice-free ground above sea level
        ! also has H_grnd > 0.
        if (z_sl-z_bed .gt. 0.0) then 
            ! Grounded below sea level, diagnose overburden minus water thickness
            H_grnd = H_eff - rho_sw_ice*(z_sl-z_bed)
        else
            ! Grounded above sea level, simply sum elevation above sea level and ice thickness
            H_grnd = H_eff + (z_bed-z_sl)
        end if 

        ! ajr: to test somewhere eventually, more closely follows Gladstone et al (2010), Leguy et al (2021)
        !H_grnd = -( (z_sl-z_bed) - rho_ice_sw*H_eff )

        return 

    end subroutine calc_H_grnd

    elemental subroutine calc_H_af(H_af,H_ice,f_ice,z_bed,z_sl,rho_ice,rho_sw,use_f_ice)
        ! Calculate ice thickness above flotation, H_af

        implicit none 

        real(wp), intent(INOUT) :: H_af
        real(wp), intent(IN)    :: H_ice
        real(wp), intent(IN)    :: f_ice
        real(wp), intent(IN)    :: z_bed
        real(wp), intent(IN)    :: z_sl
        real(wp), intent(IN)    :: rho_ice 
        real(wp), intent(IN)    :: rho_sw
        
        logical,  intent(IN), optional :: use_f_ice

        ! Local variables   
        real(wp) :: rho_sw_ice 
        real(wp) :: H_eff 
        logical  :: use_f_ice_now 

        ! By default use f_ice to set points with f_ice < 1 to H_eff=0.0 
        use_f_ice_now = .TRUE. 
        if (present(use_f_ice)) use_f_ice_now = use_f_ice 

        rho_sw_ice = rho_sw/rho_ice ! Ratio of density of seawater to ice [--]
        
        ! Get effective ice thickness
        if (use_f_ice_now) then 
            call calc_H_eff(H_eff,H_ice,f_ice,set_frac_zero=.TRUE.)
        else 
            H_eff = H_ice 
        end if 

        ! Calculate ice thickness above flotation 
        ! ie, total ice column thickness minus the corresponding depth of seawater
        H_af = H_eff - rho_sw_ice*max(z_sl-z_bed,0.0_wp)

        ! Limit to positive values, since when it is negative, it is completely floating
        H_af = max(H_af,0.0_wp)

        return 

    end subroutine calc_H_af

    subroutine calc_f_grnd_subgrid_area_aa(f_grnd,H_grnd,gz_nx,boundaries)
        ! Use H_grnd to determined grounded area fraction of grid point.

        implicit none
        
        real(wp), intent(OUT) :: f_grnd(:,:)        ! aa-nodes 
        real(wp), intent(IN)  :: H_grnd(:,:)        ! aa-nodes
        integer,  intent(IN)  :: gz_nx              ! Number of interpolation points per side (nx*nx)
        character(len=*), intent(IN) :: boundaries

        ! Local variables
        integer  :: i, j, nx, ny
        integer  :: im1, ip1, jm1, jp1 
        real(wp) :: Hg_int(gz_nx,gz_nx)

        !integer, parameter :: nx_interp = 15

        nx = size(H_grnd,1)
        ny = size(H_grnd,2) 

        ! First binary estimate of f_grnd based on aa-nodes
        f_grnd = 1.0
        where (H_grnd .lt. 0.0) f_grnd = 0.0
        
        ! Find grounding line cells and determine fraction 
        do j = 1, ny 
        do i = 1, nx

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

            if (maxval(H_grnd(im1:ip1,jm1:jp1)) .ge. 0.0 .and. minval(H_grnd(im1:ip1,jm1:jp1)) .lt. 0.0) then 
                ! Point contains grounding line, get grounded area  
                
                call calc_subgrid_array(Hg_int, H_grnd,gz_nx,i,j,im1,ip1,jm1,jp1)
                
                ! Calculate weighted fraction (assume all points have equal weight)
                f_grnd(i,j) = real(count(Hg_int .ge. 0.0),wp) / real(gz_nx*gz_nx,wp)

            end if 

        end do 
        end do 

        return
        
    end subroutine calc_f_grnd_subgrid_area_aa
    
    subroutine calc_f_grnd_subgrid_area(f_grnd,f_grnd_acx,f_grnd_acy,H_grnd,gz_nx,boundaries)
        ! Use H_grnd to determined grounded area fraction of grid point.

        implicit none
        
        real(wp), intent(OUT) :: f_grnd(:,:)        ! aa-nodes 
        real(wp), intent(OUT) :: f_grnd_acx(:,:)    ! ac-nodes
        real(wp), intent(OUT) :: f_grnd_acy(:,:)    ! ac-nodes
        real(wp), intent(IN)  :: H_grnd(:,:)        ! aa-nodes
        integer,  intent(IN)  :: gz_nx          ! Number of interpolation points per side (nx*nx)
        character(len=*), intent(IN) :: boundaries
        
        ! Local variables
        integer  :: i, j, nx, ny
        real(wp) :: Hg_1, Hg_2, Hg_3, Hg_4
        real(wp) :: Hg_min, Hg_max  
        integer  :: im1, ip1, jm1, jp1 
        real(wp) :: Hg_int(gz_nx,gz_nx)

        !integer, parameter :: nx_interp = 15

        nx = size(H_grnd,1)
        ny = size(H_grnd,2) 

        ! Initialize all masks to zero (fully floating) to start
        f_grnd     = 0.0_wp 
        f_grnd_acx = 0.0_wp 
        f_grnd_acy = 0.0_wp 

        ! Find grounding line cells and determine fraction 
        !!$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,Hg_1,Hg_2,Hg_3,Hg_4,Hg_max,Hg_min)
        do j = 1, ny 
        do i = 1, nx

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

            ! === f_grnd at aa-nodes ===

            ! Calculate Hg at corners (ab-nodes)
            Hg_1 = 0.25_wp*(H_grnd(i,j) + H_grnd(ip1,j) + H_grnd(ip1,jp1) + H_grnd(i,jp1))
            Hg_2 = 0.25_wp*(H_grnd(i,j) + H_grnd(im1,j) + H_grnd(im1,jp1) + H_grnd(i,jp1))
            Hg_3 = 0.25_wp*(H_grnd(i,j) + H_grnd(im1,j) + H_grnd(im1,jm1) + H_grnd(i,jm1))
            Hg_4 = 0.25_wp*(H_grnd(i,j) + H_grnd(ip1,j) + H_grnd(ip1,jm1) + H_grnd(i,jm1))
            
            Hg_min = min(Hg_1,Hg_2,Hg_3,Hg_4)
            Hg_max = max(Hg_1,Hg_2,Hg_3,Hg_4)

            if (Hg_max .ge. 0.0 .and. Hg_min .lt. 0.0) then 
                ! Point contains grounding line, get grounded area  
                
                call calc_subgrid_array_cell(Hg_int,Hg_1,Hg_2,Hg_3,Hg_4,gz_nx)

                ! Calculate weighted fraction (assume all points have equal weight)
                f_grnd(i,j) = real(count(Hg_int .ge. 0.0),wp) / real(gz_nx*gz_nx,wp)

            else if (Hg_min .ge. 0.0) then 
                ! Fully grounded point

                f_grnd(i,j) = 1.0_wp 

            end if 

            ! === f_grnd at acx nodes === 

            ! First, calculate Hg at corners (acy-nodes)
            Hg_1 = 0.5_wp*(H_grnd(ip1,j) + H_grnd(ip1,jp1))
            Hg_2 = 0.5_wp*(H_grnd(i,j)   + H_grnd(i,jp1))
            Hg_3 = 0.5_wp*(H_grnd(i,j)   + H_grnd(i,jm1))
            Hg_4 = 0.5_wp*(H_grnd(ip1,j) + H_grnd(ip1,jm1))
            
            Hg_min = min(Hg_1,Hg_2,Hg_3,Hg_4)
            Hg_max = max(Hg_1,Hg_2,Hg_3,Hg_4)

            if (Hg_max .ge. 0.0 .and. Hg_min .lt. 0.0) then 
                ! Point contains grounding line, get grounded area  
                
                call calc_subgrid_array_cell(Hg_int,Hg_1,Hg_2,Hg_3,Hg_4,gz_nx)

                ! Calculate weighted fraction (assume all points have equal weight)
                f_grnd_acx(i,j) = real(count(Hg_int .ge. 0.0),wp) / real(gz_nx*gz_nx,wp)

            else if (Hg_min .ge. 0.0) then
                ! Purely grounded point 

                f_grnd_acx(i,j) = 1.0_wp 

            end if 

            ! === f_grnd at acy-nodes ===
        
            ! First, calculate Hg at corners (acx-nodes)
            Hg_1 = 0.5_wp*(H_grnd(i,jp1)   + H_grnd(ip1,jp1))
            Hg_2 = 0.5_wp*(H_grnd(im1,jp1) + H_grnd(i,jp1))
            Hg_3 = 0.5_wp*(H_grnd(im1,j)   + H_grnd(i,j))
            Hg_4 = 0.5_wp*(H_grnd(ip1,j)   + H_grnd(i,j))
            
            Hg_min = min(Hg_1,Hg_2,Hg_3,Hg_4)
            Hg_max = max(Hg_1,Hg_2,Hg_3,Hg_4)

            if (Hg_max .ge. 0.0 .and. Hg_min .lt. 0.0) then 
                ! Point contains grounding line, get grounded area  
                
                call calc_subgrid_array_cell(Hg_int,Hg_1,Hg_2,Hg_3,Hg_4,gz_nx)

                ! Calculate weighted fraction (assume all points have equal weight)
                f_grnd_acy(i,j) = real(count(Hg_int .ge. 0.0),wp) / real(gz_nx*gz_nx,wp)

            else if (Hg_min .ge. 0.0) then 
                ! Purely grounded point 
                    
                f_grnd_acy(i,j) = 1.0_wp 
                
            end if 

        end do 
        end do 
        !!$omp end parallel do


if (.TRUE.) then 
    ! Replace subgrid acx/acy estimates with linear average to ac-nodes 

        ! acx-nodes 
        do j = 1, ny 
        do i = 1, nx-1
            f_grnd_acx(i,j) = 0.5_wp*(f_grnd(i,j) + f_grnd(i+1,j))
        end do 
        end do
        f_grnd_acx(nx,:) = f_grnd_acx(nx-1,:) 

        ! acy-nodes 
        do j = 1, ny-1 
        do i = 1, nx
            f_grnd_acy(i,j) = 0.5_wp*(f_grnd(i,j) + f_grnd(i,j+1))
        end do 
        end do
        f_grnd_acy(:,ny) = f_grnd_acy(:,ny-1) 

end if 

        return
        
    end subroutine calc_f_grnd_subgrid_area
    
    subroutine calc_f_grnd_subgrid_linear(f_grnd,f_grnd_x,f_grnd_y,H_grnd)
        ! Calculate the grounded fraction of a cell in the x- and y-directions
        ! at the ac nodes
        !
        ! Given point 1 is grounded and point 2 is floating, we are looking for
        ! when H_grnd=0. The equation along an x-axis between point 1 (grounded point) 
        ! and point 2 (floating point) for H_grnd is:
        ! H_grnd = H_grnd_1 + f*(H_grnd_2-H_grnd_1)

        ! To find fraction f for when H_grnd == 0, rearrange equation:
        ! 0 = H_grnd_1 + f*(H_grnd_2-H_grnd_1)
        ! f = -H_grnd_1 / (H_grnd_2-H_grnd_1)
        ! where f is the distance along the 0:1 axis between the first and second points. 
        
        implicit none 

        real(wp), intent(OUT) :: f_grnd(:,:)
        real(wp), intent(OUT) :: f_grnd_x(:,:)
        real(wp), intent(OUT) :: f_grnd_y(:,:)
        real(wp), intent(IN)  :: H_grnd(:,:)

        ! Local variables  
        integer :: i, j, nx, ny 
        real(wp) :: H_grnd_1, H_grnd_2

        nx = size(f_grnd,1)
        ny = size(f_grnd,2)

        ! Central aa-node
        f_grnd = 1.0
        where (H_grnd < 0.0) f_grnd = 0.0
        
        ! x-direction, ac-node
        f_grnd_x = 1.0
        !!$omp parallel do collapse(2) private(i,j,H_grnd_1,H_grnd_2)
        do j = 1, ny 
        do i = 1, nx-1 

            if (H_grnd(i,j) .gt. 0.0 .and. H_grnd(i+1,j) .le. 0.0) then 
                ! Point is grounded, neighbor is floating 

                H_grnd_1 = H_grnd(i,j) 
                H_grnd_2 = H_grnd(i+1,j) 

                ! Calculate fraction 
                f_grnd_x(i,j) = -H_grnd_1 / (H_grnd_2 - H_grnd_1)

            else if (H_grnd(i,j) .le. 0.0 .and. H_grnd(i+1,j) .gt. 0.0) then 
                ! Point is floating, neighbor is grounded 

                H_grnd_1 = H_grnd(i+1,j) 
                H_grnd_2 = H_grnd(i,j) 

                ! Calculate fraction 
                f_grnd_x(i,j) = -H_grnd_1 / (H_grnd_2 - H_grnd_1)

            else if (H_grnd(i,j) .le. 0.0 .and. H_grnd(i+1,j) .le. 0.0) then 
                ! Point is floating, neighbor is floating
                f_grnd_x(i,j) = 0.0 

            else 
                ! Point is grounded, neighbor is grounded
                f_grnd_x(i,j) = 1.0 

            end if 

        end do 
        end do 
        !!$omp end parallel do

        ! y-direction, ac-node
        f_grnd_y = 1.0
        !!$omp parallel do collapse(2) private(i,j,H_grnd_1,H_grnd_2)
        do j = 1, ny-1 
        do i = 1, nx 

            if (H_grnd(i,j) .gt. 0.0 .and. H_grnd(i,j+1) .le. 0.0) then 
                ! Point is grounded, neighbor is floating 

                H_grnd_1 = H_grnd(i,j) 
                H_grnd_2 = H_grnd(i,j+1) 

                ! Calculate fraction 
                f_grnd_y(i,j) = -H_grnd_1 / (H_grnd_2 - H_grnd_1)

            else if (H_grnd(i,j) .le. 0.0 .and. H_grnd(i,j+1) .gt. 0.0) then 
                ! Point is floating, neighbor is grounded 

                H_grnd_1 = H_grnd(i,j+1) 
                H_grnd_2 = H_grnd(i,j) 

                ! Calculate fraction 
                f_grnd_y(i,j) = -H_grnd_1 / (H_grnd_2 - H_grnd_1)
                
            else if (H_grnd(i,j) .le. 0.0 .and. H_grnd(i,j+1) .le. 0.0) then 
                ! Point is floating, neighbor is floating
                f_grnd_y(i,j) = 0.0 

            else 
                ! Point is grounded, neighbor is grounded
                f_grnd_y(i,j) = 1.0 

            end if 

        end do 
        end do 
        !!$omp end parallel do

        ! Set boundary points equal to neighbor for aesthetics 
        f_grnd_x(nx,:) = f_grnd_x(nx-1,:) 
        f_grnd_y(:,ny) = f_grnd_y(:,ny-1) 
        
        return 

    end subroutine calc_f_grnd_subgrid_linear
    
    subroutine calc_f_grnd_pinning_points(f_grnd,H_ice,f_ice,z_bed,z_bed_sd,z_sl,rho_ice,rho_sw)
        ! For floating points, determine how much of bed could be
        ! touching the base of the ice shelf due to subgrid pinning
        ! points following the distribution z = N(z_bed,z_bed_sd)

        implicit none

        real(wp), intent(OUT) :: f_grnd(:,:) 
        real(wp), intent(IN)  :: H_ice(:,:) 
        real(wp), intent(IN)  :: f_ice(:,:) 
        real(wp), intent(IN)  :: z_bed(:,:) 
        real(wp), intent(IN)  :: z_bed_sd(:,:) 
        real(wp), intent(IN)  :: z_sl(:,:) 
        real(wp), intent(IN)  :: rho_ice 
        real(wp), intent(IN)  :: rho_sw 

        ! Local variables 
        integer :: i, j, nx, ny 
        real(wp) :: x, mu, sigma
        real(wp) :: rho_ice_sw
        real(wp) :: z_base_now
        real(wp) :: H_eff_now

        nx = size(f_grnd,1) 
        ny = size(f_grnd,2) 

        rho_ice_sw = rho_ice/rho_sw ! Ratio of density of ice to seawater [--]
        
        ! Initialize f_grnd to zero everywhere
        f_grnd = 0.0_wp 

        !!$omp parallel do collapse(2) private(i,j,H_eff_now,z_base_now,x,mu,sigma)
        do j = 1, ny 
        do i = 1, nx 

            ! Get effective ice thickness
            call calc_H_eff(H_eff_now,H_ice(i,j),f_ice(i,j))

            ! Determine depth of base of ice shelf, assuming
            ! for now that it is floating. 
            z_base_now = z_sl(i,j) - rho_ice_sw*H_eff_now

            ! Check if it is actually floating: 

            if (z_base_now .gt. z_bed(i,j)) then 
                ! Floating point, calculate pinning fraction 

                ! Define parameters to calculate inverse CDF
                x     = z_base_now 
                mu    = z_bed(i,j) 
                sigma = z_bed_sd(i,j) 

                if (sigma .eq. 0.0_wp) then 
                    ! sigma not available, set f_grnd to zero 

                    f_grnd(i,j) = 0.0_wp 

                else

                    ! Calculate inverse cdf at z = ice_base, ie,
                    ! probability to have points at or higher than z,
                    ! which would be grounded. 

                    !f_grnd(i,j) = cdf(x,mu,sigma,inv=.TRUE.)

                    ! Alternatively, use approximation given by
                    ! Pollard and DeConto (2012), Eq. 13:

                    f_grnd(i,j) = 0.5_wp*max(0.0_wp,(1.0_wp-(x-mu)/sigma))

                end if 

            end if
            
        end do 
        end do 
        !!$omp end parallel do

        return

    end subroutine calc_f_grnd_pinning_points

    subroutine remove_englacial_lakes(H_ice,z_bed,z_srf,z_sl,rho_ice,rho_sw)
        ! Diagnose where ice should be grounded, but a gap exists.
        ! In these locations, increase ice thickness in order
        ! to fill in the englacial lake. 
        ! Also check if ice thickness is below bedrock, in these
        ! places reduce ice thickness. 
        ! Used when loading datasets. 

        implicit none

        real(wp), intent(INOUT) :: H_ice(:,:) 
        real(wp), intent(IN)    :: z_bed(:,:) 
        real(wp), intent(IN)    :: z_srf(:,:) 
        real(wp), intent(IN)    :: z_sl(:,:) 
        real(wp), intent(IN)    :: rho_ice 
        real(wp), intent(IN)    :: rho_sw 

        ! Local variables 
        integer :: i, j, nx , ny 
        real(wp) :: H_grnd_now 
        real(wp) :: rho_sw_ice
        
        rho_sw_ice = rho_sw/rho_ice ! Ratio of density of seawater to ice [--]
        
        nx = size(H_ice,1)
        ny = size(H_ice,2)

        ! Surface elevation is available, use to remove englacial lakes
        ! from ice thickness field. 
        do j = 1, ny 
        do i = 1, nx 

            ! Calculate H_grnd (ice thickness overburden)
            ! (ie, diagnose whether ice is grounded or floating)
            H_grnd_now = H_ice(i,j) - rho_sw_ice*max(z_sl(i,j)-z_bed(i,j),0.0_wp)

            if (H_grnd_now .gt. 0.0_wp .and. z_srf(i,j) - H_ice(i,j) .gt. z_bed(i,j)) then 
                ! Ice should be grounded, but a gap exists, 
                ! so a subglacial lake exists, increase thickness.

                H_ice(i,j) = z_srf(i,j) - z_bed(i,j) 

            else if (z_srf(i,j) - H_ice(i,j) .lt. z_bed(i,j)) then 
                ! Ice is erroneously thick, reduce it

                H_ice(i,j) = z_srf(i,j) - z_bed(i,j) 

            end if 

        end do 
        end do

        return

    end subroutine remove_englacial_lakes

    subroutine calc_distance_to_ice_margin(dist_mrgn,f_ice,dx,boundaries)
        ! Calculate distance to the ice margin
        
        ! Note: this subroutine is a wrapper that calls the
        ! grounding-line distance routine, since the algorithm
        ! works the same way. Simply substitute f_ice for f_grnd. 

        implicit none 

        real(wp), intent(OUT) :: dist_mrgn(:,:) ! [km] Distance to grounding line
        real(wp), intent(IN)  :: f_ice(:,:)     ! [1]  Fraction of grid-cell ice coverage 
        real(wp), intent(IN)  :: dx             ! [m]  Grid resolution (assume dy=dx)
        character(len=*), intent(IN) :: boundaries
    
        call calc_distance_to_grounding_line(dist_mrgn,f_ice,dx,boundaries)

        return 

    end subroutine calc_distance_to_ice_margin
    
    subroutine calc_distance_to_grounding_line(dist_gl,f_grnd,dx,boundaries)
        ! Calculate distance to the grounding line 
        
        implicit none 

        real(wp), intent(OUT) :: dist_gl(:,:)   ! [km] Distance to grounding line
        real(wp), intent(IN)  :: f_grnd(:,:)    ! [1]  Grounded grid-cell fraction 
        real(wp), intent(IN)  :: dx             ! [m]  Grid resolution (assume dy=dx)
        character(len=*), intent(IN) :: boundaries

        ! Local variables 
        integer  :: i, j, nx, ny, q
        integer  :: im1, ip1, jm1, jp1
        real(wp) :: dist_direct_min 
        real(wp) :: dist_corners_min
        real(wp) :: dx_km 
        real(wp) :: dists(8) 

        real(wp), parameter :: dist_max = 1e10          ! [km]
        real(wp), parameter :: sqrt_2   = sqrt(2.0_wp) 
        integer,  parameter :: iter_max = 1000

        real(wp), allocatable :: dist_gl_ref(:,:) 

        nx = size(dist_gl,1)
        ny = size(dist_gl,2)

        allocate(dist_gl_ref(nx,ny)) 

        ! Units
        dx_km = dx*1e-3 

        ! 0: Assign a big distance to all points to start  ======================

        dist_gl  = dist_max

        ! 1. Next, determine grounding line =====================================

        !!$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1)
        do j = 1, ny 
        do i = 1, nx

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

            ! Grounded point or partially floating point with floating neighbors
            if (f_grnd(i,j) .gt. 0.0 .and. &
                (f_grnd(im1,j) .eq. 0.0 .or. f_grnd(ip1,j) .eq. 0.0 .or. &
                 f_grnd(i,jm1) .eq. 0.0 .or. f_grnd(i,jp1) .eq. 0.0) ) then 
                
                dist_gl(i,j)  = 0.0_wp 

            end if 

        end do 
        end do 
        !!$omp end parallel do

        ! 2. Next, determine distances to grounding line ======================
        
        do q = 1, iter_max  
            ! Iterate distance of one neighbor at a time until grid is filled in

            dist_gl_ref = dist_gl

            !!$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,dists,dist_direct_min,dist_corners_min)
            do j = 1, ny 
            do i = 1, nx

                if (dist_gl(i,j) .eq. dist_max) then 
                    ! Distance needs to be determined for this point 

                    ! Get neighbor indices
                    call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

                    ! Get distances to direct and corner neighbors
                    dists = [dist_gl_ref(im1,j),dist_gl_ref(ip1,j), &       ! Direct neighbors
                             dist_gl_ref(i,jm1),dist_gl_ref(i,jp1), &       ! Direct neighbors
                             dist_gl_ref(im1,jp1),dist_gl_ref(ip1,jp1), &   ! Corner neighbors
                             dist_gl_ref(im1,jm1),dist_gl_ref(ip1,jm1)]     ! Corner neighbors

                    if (count(dists .lt. dist_max) .ge. 2) then 
                        ! Some neighbors have been calculated already 
                        ! Note: check for at least 2 neighbors already
                        ! calculated to reduce errors in points with 
                        ! neighbors on both sides with similar distances. This 
                        ! Waiting for even more neighbors would reduce errors
                        ! further but greatly increases iterations. 

                        ! Determine minimum distance to grounding line for 
                        ! direct and diagonal neighbors, separately.
                        dist_direct_min  = minval(dists(1:4))
                        dist_corners_min = minval(dists(5:8))

                        if (dist_direct_min .le. dist_corners_min) then 
                            ! Assume nearest path to grounding line is along
                            ! direct neighbor path - add dx to dist for this point 

                            dist_gl(i,j) = dist_direct_min + dx_km 

                        else
                            ! Assume nearest path to grounding line is via 
                            ! a diagonal neighbor - add sqrt(2) to dist for this point

                            dist_gl(i,j) = dist_corners_min + sqrt_2*dx_km

                        end if 

                    end if 
                end if

            end do 
            end do
            !!$omp end parallel do
            
            if (count(dist_gl .eq. dist_max) .eq. 0) then 
                ! No more points to check 
                exit 
            end if 

        end do

        ! Set all floating-point distances to negative values 
        where (f_grnd .eq. 0.0_wp) 
            dist_gl  = -dist_gl
        end where

        
        return 

    end subroutine calc_distance_to_grounding_line
    
    subroutine calc_grounding_line_zone(mask_grz,dist_gl,dist_grz)
        ! Define a mask of the grounding-zone region 
        ! mask_grz == -2: floating-point out of grounding zone
        ! mask_grz == -1: floating-point in grounding zone
        ! mask_grz ==  0: grounding line 
        ! mask_grz ==  1: grounded-point in grounding zone
        ! mask_grz ==  2: grounded-point out of grounding zone

        implicit none

        integer,  intent(OUT) :: mask_grz(:,:)  ! [1]  Grounding-zone mask 
        real(wp), intent(IN)  :: dist_gl(:,:)   ! [km] Distance to grounding line
        real(wp), intent(IN)  :: dist_grz       ! [km] Distance to define as grounding zone

        ! Finally define the grounding-zone mask too
        where(dist_gl .eq. 0.0_wp)
            mask_grz =  0 
        else where(dist_gl .lt. 0.0_wp .and. abs(dist_gl) .le. dist_grz)
            mask_grz = -1
        else where(dist_gl .gt. 0.0_wp .and.     dist_gl  .le. dist_grz)
            mask_grz =  1
        else where(dist_gl .lt. 0.0_wp)
            mask_grz = -2
        elsewhere
            mask_grz =  2
        end where

        return

    end subroutine calc_grounding_line_zone

    subroutine calc_bmb_total(bmb,bmb_grnd,bmb_shlf,H_ice,H_grnd,f_grnd,gz_Hg0,gz_Hg1, &
                                                            gz_nx,bmb_gl_method,grounded_melt,mask_pd,boundaries)

        implicit none 

        real(wp),         intent(OUT) :: bmb(:,:) 
        real(wp),         intent(IN)  :: bmb_grnd(:,:) 
        real(wp),         intent(IN)  :: bmb_shlf(:,:) 
        real(wp),         intent(IN)  :: H_ice(:,:)
        real(wp),         intent(IN)  :: H_grnd(:,:)
        real(wp),         intent(IN)  :: f_grnd(:,:) 
        real(wp),         intent(IN)  :: gz_Hg0
        real(wp),         intent(IN)  :: gz_Hg1
        integer,          intent(IN)  :: gz_nx
        character(len=*), intent(IN)  :: bmb_gl_method 
        ! for optimization (optional, check!)
        logical,  intent(IN), optional  :: grounded_melt
        integer,  intent(IN), optional  :: mask_pd(:,:)
        character(len=*), intent(IN)  :: boundaries 

        ! Local variables
        integer    :: i, j, nx, ny
        integer    :: n_float  
        real(wp) :: bmb_shlf_now 

        nx = size(bmb,1)
        ny = size(bmb,2) 

        ! Combine floating and grounded parts into one field =========================
        ! Apply the floating basal mass balance according 
        ! to different subgridding options at the grounding line
        ! (following the notation of Leguy et al., 2021 - see Fig. 3)
        
        select case(bmb_gl_method)

            case("fcmp")
                ! Flotation criterion melt parameterization 
                ! Apply full bmb_shlf value where flotation criterion is met 

                where(H_grnd .le. 0.0_wp)
                    
                    bmb = bmb_shlf

                elsewhere

                    bmb = bmb_grnd 

                end where 

            case("fmp")
                ! Full melt parameterization
                ! Apply full bmb_shlf value to any cell that is at least
                ! partially floating. Perhaps unrealistic.

                where(f_grnd .lt. 1.0_wp)

                    bmb = bmb_shlf 

                elsewhere

                    bmb = bmb_grnd 

                end where 

            case("pmp")
                ! Partial melt parameterization
                ! Apply bmb_shlf to floating fraction of cell 

                where(f_grnd .lt. 1.0_wp)

                    bmb = f_grnd*bmb_grnd + (1.0_wp-f_grnd)*bmb_shlf 

                elsewhere

                    bmb = bmb_grnd 

                end where 

            case("pmpt")
                ! Partial melt parameterization with tidal grounding zone

                call calc_bmb_gl_pmpt(bmb,bmb_grnd,bmb_shlf,H_grnd,gz_Hg0,gz_Hg1,gz_nx,boundaries)

            case("nmp")
                ! No melt parameterization
                ! Apply bmb_shlf only where fully floating diagnosed

                where(f_grnd .eq. 0.0_wp)

                    bmb = bmb_shlf 

                elsewhere

                    bmb = bmb_grnd 

                end where 

        end select

        ! Melting in grounded zone for optimization
        if (present(mask_pd) .and. grounded_melt) then
            where(mask_pd .eq. mask_bed_float) bmb = bmb_shlf
            where(mask_pd .eq. mask_bed_ocean) bmb = bmb_shlf
        end if

        ! For aesthetics, also make sure that bmb is zero on ice-free land
        where (H_grnd .gt. 0.0_wp .and. H_ice .eq. 0.0_wp) bmb = 0.0_wp 

        return 

    end subroutine calc_bmb_total

    subroutine calc_fmb_total(fmb,fmb_shlf,bmb_shlf,H_ice,H_grnd,f_ice, &
                                fmb_method,fmb_scale,rho_ice,rho_sw,dx,boundaries)

        implicit none 

        real(wp), intent(OUT) :: fmb(:,:) 
        real(wp), intent(IN)  :: fmb_shlf(:,:) 
        real(wp), intent(IN)  :: bmb_shlf(:,:) 
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: H_grnd(:,:) 
        real(wp), intent(IN)  :: f_ice(:,:)
        integer,  intent(IN)  :: fmb_method 
        real(wp), intent(IN)  :: fmb_scale
        real(wp), intent(IN)  :: rho_ice
        real(wp), intent(IN)  :: rho_sw
        real(wp), intent(IN)  :: dx
        character(len=*), intent(IN) :: boundaries
    
        ! Local variables
        integer    :: i, j, nx, ny, n_margin
        integer    :: im1, ip1, jm1, jp1
        real(wp) :: H_eff, dz  
        real(wp) :: area_flt
        real(wp) :: area_tot
        real(wp) :: bmb_eff 
        logical  :: mask(4) 

        real(wp) :: rho_ice_sw 
        
        rho_ice_sw = rho_ice/rho_sw ! Ratio of density of ice to seawater [--]
        
        ! Total cell area 
        area_tot = dx*dx 

        nx = size(fmb,1)
        ny = size(fmb,2) 

        select case(fmb_method)

            case(0)
                ! fmb provided by boundary field fmb_shlf 

                fmb = fmb_shlf

            case(1,2) 
                ! Calculate fmb (1) as proportional to local bmb_shlf value, 
                ! or (2) from a reference field, but scaled to the area
                ! of the grid cell itself where it will be applied

                !!$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,mask,n_margin,dz,area_flt,bmb_eff)
                do j = 1, ny 
                do i = 1, nx 

                    ! Get neighbor indices
                    call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

                    ! Get mask of neighbors that are ice free 
                    mask = ( [H_ice(im1,j),H_ice(ip1,j),H_ice(i,jm1),H_ice(i,jp1)].eq.0.0 )

                    ! Count how many neighbors are ice free 
                    n_margin = count(mask) 

                    if (H_ice(i,j) .gt. 0.0_wp .and. &
                        H_grnd(i,j) .lt. H_ice(i,j) .and. &
                                            n_margin .gt. 0) then 
                        ! Cell is ice-covered, [grounded below sea level or floating] and at the ice margin

                        ! Get effective ice thickness
                        call calc_H_eff(H_eff,H_ice(i,j),f_ice(i,j))

                        ! Determine depth of adjacent water using centered cell information
                        if (H_grnd(i,j) .lt. 0.0_wp) then 
                            ! Cell is floating, calculate submerged ice thickness 
                            
                            dz = (H_eff*rho_ice_sw)
                            
                        else 
                            ! Cell is grounded, recover depth of seawater

                            dz = max( (H_eff - H_grnd(i,j)) * rho_ice_sw, 0.0_wp)

                        end if 

                        ! Get area of ice submerged and adjacent to seawater
                        area_flt = real(n_margin,wp)*dz*dx 

                        if (fmb_method .eq. 1) then
                            ! Also calculate the mean bmb_shlf value for the ice-free neighbors 
                            bmb_eff = sum([bmb_shlf(im1,j),bmb_shlf(ip1,j),bmb_shlf(i,jm1),bmb_shlf(i,jp1)], &
                                            mask=mask) / real(n_margin,wp)

                            ! Finally calculate the effective front mass balance rate 

                            fmb(i,j) = bmb_eff*(area_flt/area_tot)*fmb_scale

                        else if (fmb_method .eq. 2) then

                            ! Finally calculate the scaled front mass balance rate 

                            fmb(i,j) = fmb_shlf(i,j)*(area_flt/area_tot)
                        end if

                    else 
                        ! Set front mass balance equal to zero 

                        fmb(i,j) = 0.0_wp 

                    end if 

                end do 
                end do
                !!$omp end parallel do

            case DEFAULT 

                write(*,*) "calc_fmb_total:: Error: fmb_method not recognized."
                write(*,*) "fmb_method = ", fmb_method 
                stop 

        end select

        return 

    end subroutine calc_fmb_total

    subroutine calc_bmb_gl_pmpt(bmb,bmb_grnd,bmb_shlf,H_grnd,gz_Hg0,gz_Hg1,nxi,boundaries)
        ! Calculate basal mass balance, with bmb at the grounding line
        ! determined via subgrid calculation of flotation and parameterization
        ! for tidal-induced grounded melt. 

        implicit none
        
        real(wp), intent(OUT) :: bmb(:,:)           ! aa-nodes 
        real(wp), intent(IN)  :: bmb_grnd(:,:)      ! aa-nodes 
        real(wp), intent(IN)  :: bmb_shlf(:,:)      ! aa-nodes 
        real(wp), intent(IN)  :: H_grnd(:,:)        ! aa-nodes
        real(wp), intent(IN)  :: gz_Hg0             ! Lower limit in H_grnd for grounding zone
        real(wp), intent(IN)  :: gz_Hg1             ! Upper limit in H_grnd for grounding zone
        integer,  intent(IN)  :: nxi                ! Number of interpolation points per side (nxi*nxi)
        character(len=*), intent(IN) :: boundaries 

        ! Local variables
        integer  :: i, j, i1, j1, nx, ny
        integer  :: im1, ip1, jm1, jp1 
        real(wp) :: Hg_1, Hg_2, Hg_3, Hg_4, Hg_mid  
        real(wp) :: wt 

        real(wp), allocatable :: Hg_int(:,:)
        real(wp), allocatable :: bmb_int(:,:)

        ! Consistency check
        if (gz_Hg0 .gt. 0.0) then 
            write(io_unit_err,*) "calc_bmb_gl_pmpt:: Error: lower limit on grounding zone must be <= 0.0."
            write(io_unit_err,*) "gz_Hg0 = ", gz_Hg0
            stop 
        end if 
        if (gz_Hg1 .lt. 0.0) then 
            write(io_unit_err,*) "calc_bmb_gl_pmpt:: Error: upper limit on grounding zone must be >= 0.0."
            write(io_unit_err,*) "gz_Hg1 = ", gz_Hg1
            stop 
        end if 

        nx = size(H_grnd,1)
        ny = size(H_grnd,2) 

        ! Allocate subgrid arrays
        allocate(Hg_int(nxi,nxi))
        allocate(bmb_int(nxi,nxi))

        !!$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,Hg_1,Hg_2,Hg_3,Hg_4,Hg_int,i1,j1,wt)
        do j = 1, ny 
        do i = 1, nx

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

            if (minval(H_grnd(im1:ip1,jm1:jp1)) .ge. gz_Hg1) then 
                ! Entire cell is grounded

                bmb(i,j) = bmb_grnd(i,j)
            
            else if (maxval(H_grnd(im1:ip1,jm1:jp1)) .lt. gz_Hg0) then 
                ! Entire cell is floating

                bmb(i,j) = bmb_shlf(i,j) 

            else
                ! Point contains the grounding zone
                
                ! Calculate subgrid values of H_grnd
                call calc_subgrid_array(Hg_int,H_grnd,nxi,i,j,im1,ip1,jm1,jp1)

                ! Calculate individual bmb values for each subgrid point
                do j1 = 1, nxi
                do i1 = 1, nxi

                    if (Hg_int(i1,j1) .lt. gz_Hg0) then
                        ! Floating, outside of grounding zone 
                        wt = 0.0
                    else if (Hg_int(i1,j1) .ge. gz_Hg1) then
                        ! Grounded, outside of grounding zone
                        wt = 1.0 
                    else
                        ! Within grounding zone
                        wt = (Hg_int(i1,j1)-gz_Hg0) / (gz_Hg1 - gz_Hg0)
                    end if 

                    ! Get subgrid bmb weighted between floating and grounded contributions
                    bmb_int(i1,j1) = wt*bmb_grnd(i,j) + (1.0-wt)*bmb_shlf(i,j) 

                end do
                end do

                ! Get the mean bmb rate for the entire cell
                bmb(i,j) = sum(bmb_int) / real(nxi*nxi,wp)

            end if 

        end do 
        end do
        !!$omp end parallel do

        return
        
    end subroutine calc_bmb_gl_pmpt
    
!! f_grnd calculations from IMAU-ICE / CISM 

! == Routines for determining the grounded fraction on all four grids
  
  subroutine determine_grounded_fractions(f_grnd,f_grnd_acx,f_grnd_acy,f_grnd_ab,H_grnd,boundaries)
    ! Determine the grounded fraction of centered and staggered grid points
    ! Uses the bilinear interpolation scheme (with analytical solutions) 
    ! from CISM (Leguy et al., 2021), as adapted from IMAU-ICE v2.0 code (rev. 4776833b)
    
    implicit none
    
    real(wp), intent(OUT) :: f_grnd(:,:) 
    real(wp), intent(OUT), optional :: f_grnd_acx(:,:) 
    real(wp), intent(OUT), optional :: f_grnd_acy(:,:) 
    real(wp), intent(OUT), optional :: f_grnd_ab(:,:) 
    real(wp), intent(IN)  :: H_grnd(:,:) 
    character(len=*), intent(IN) :: boundaries
    
    ! Local variables
    integer :: i, j, nx, ny 
    integer :: im1, ip1, jm1, jp1

    real(wp), allocatable :: f_grnd_NW(:,:)
    real(wp), allocatable :: f_grnd_NE(:,:)
    real(wp), allocatable :: f_grnd_SW(:,:)
    real(wp), allocatable :: f_grnd_SE(:,:)
    real(wp), allocatable :: f_flt(:,:) 

    nx = size(f_grnd,1)
    ny = size(f_grnd,2) 

    allocate(f_grnd_NW(nx,ny))
    allocate(f_grnd_NE(nx,ny))
    allocate(f_grnd_SW(nx,ny))
    allocate(f_grnd_SE(nx,ny))
    
    allocate(f_flt(nx,ny))

    ! Define aa-node variable f_flt as the flotation function
    ! following Leguy et al. (2021), Eq. 6. 
    ! Note: -H_grnd is not exactly the same, since it is in ice thickness,
    ! whereas the L21 equation is in water equivalent thickness, and it includes
    ! bedrock above sea level. But test this as it is first. 
    f_flt = -H_grnd
    
    ! Calculate grounded fractions of all four quadrants of each a-grid cell
    call determine_grounded_fractions_CISM_quads(f_grnd_NW,f_grnd_NE,f_grnd_SW,f_grnd_SE,f_flt,boundaries)
    
    ! Get grounded fractions on all four grids by averaging over the quadrants
    !!$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1)
    do j = 1, ny
    do i = 1, nx 
        
        ! Get neighbor indices
        call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

      ! aa-nodes
      f_grnd(i,j)     = 0.25_wp * (f_grnd_NW(i,j) + f_grnd_NE(i,j) + f_grnd_SW(i,j) + f_grnd_SE(i,j))
      
      if (present(f_grnd_acx)) then
        ! acx-nodes
        f_grnd_acx(i,j) = 0.25_wp * (f_grnd_NE(i,j) + f_grnd_SE(i,j) + f_grnd_NW(ip1,j) + f_grnd_SW(ip1,j))
      end if 

      if (present(f_grnd_acy)) then
        ! acy-nodes
        f_grnd_acy(i,j) = 0.25_wp * (f_grnd_NE(i,j) + f_grnd_NW(i,j) + f_grnd_SE(i,jp1) + f_grnd_SW(i,jp1))
      end if 

      if (present(f_grnd_ab)) then
        ! ab-nodes
        f_grnd_ab(i,j)  = 0.25_wp * (f_grnd_NE(i,j) + f_grnd_NW(ip1,j) + f_grnd_SE(i,jp1) + f_grnd_SW(ip1,jp1))
      end if

    end do
    end do
    !!$omp end parallel do
    
    return 

  end subroutine determine_grounded_fractions

  subroutine determine_grounded_fractions_CISM_quads(f_grnd_NW,f_grnd_NE,f_grnd_SW,f_grnd_SE,f_flt,boundaries)
    ! Calculate grounded fractions of all four quadrants of each a-grid cell
    ! (using the approach from CISM, where grounded fractions are calculated
    !  based on analytical solutions to the bilinear interpolation)
    
    implicit none
    
    real(wp), intent(OUT) :: f_grnd_NW(:,:)
    real(wp), intent(OUT) :: f_grnd_NE(:,:)
    real(wp), intent(OUT) :: f_grnd_SW(:,:)
    real(wp), intent(OUT) :: f_grnd_SE(:,:)
    real(wp), intent(IN)  :: f_flt(:,:) 
    character(len=*), intent(IN) :: boundaries

    ! Local variables:
    integer  :: i, j, ii, jj, nx, ny
    integer  :: im1, ip1, jm1, jp1  
    real(wp) :: f_NW, f_N, f_NE, f_W, f_m, f_E, f_SW, f_S, f_SE
    real(wp) :: fq_NW, fq_NE, fq_SW, fq_SE
    
    nx = size(f_flt,1)
    ny = size(f_flt,2)

    ! Calculate grounded fractions of all four quadrants of each a-grid cell
    !!$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1, f_NW,f_N,f_NE,f_W,f_m,f_E,f_SW,f_S,f_SE, fq_NW,fq_NE,fq_SW,fq_SE)
    do j = 1, ny
    do i = 1, nx
        
        ! Get neighbor indices
        call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

      f_NW = 0.25_wp * (f_flt(im1,jp1) + f_flt(i,jp1)   + f_flt(im1,j)   + f_flt(i,j))
      f_N  = 0.50_wp * (f_flt(i,jp1)   + f_flt(i,j))
      f_NE = 0.25_wp * (f_flt(i,jp1)   + f_flt(ip1,jp1) + f_flt(i,j)     + f_flt(ip1,j))
      f_W  = 0.50_wp * (f_flt(im1,j)   + f_flt(i,j))
      f_m  = f_flt(i,j) 
      f_E  = 0.50_wp * (f_flt(i,j)     + f_flt(ip1,j))
      f_SW = 0.25_wp * (f_flt(im1,j)   + f_flt(i,j)     + f_flt(im1,jm1) + f_flt(i,jm1))
      f_S  = 0.50_wp * (f_flt(i,j)     + f_flt(i,jm1))
      f_SE = 0.25_wp * (f_flt(i,j)     + f_flt(ip1,j)   + f_flt(i,jm1)   + f_flt(ip1,jm1))

      ! NW
      fq_NW = f_NW
      fq_NE = f_N
      fq_SW = f_W
      fq_SE = f_m
      call calc_fraction_above_zero( fq_NW, fq_NE,  fq_SW,  fq_SE,  f_grnd_NW(i,j) )
      
      ! NE
      fq_NW = f_N
      fq_NE = f_NE
      fq_SW = f_m
      fq_SE = f_E
      call calc_fraction_above_zero( fq_NW, fq_NE,  fq_SW,  fq_SE,  f_grnd_NE(i,j) )
      
      ! SW
      fq_NW = f_W
      fq_NE = f_m
      fq_SW = f_SW
      fq_SE = f_S
      call calc_fraction_above_zero( fq_NW, fq_NE,  fq_SW,  fq_SE,  f_grnd_SW(i,j) )
      
      ! SE
      fq_NW = f_m
      fq_NE = f_E
      fq_SW = f_S
      fq_SE = f_SE
      call calc_fraction_above_zero( fq_NW, fq_NE,  fq_SW,  fq_SE,  f_grnd_SE(i,j) )
      
    end do
    end do
    !!$omp end parallel do

    return 

  end subroutine determine_grounded_fractions_CISM_quads

  subroutine calc_fraction_above_zero( f_NW, f_NE, f_SW, f_SE, phi)
    ! Given a square with function values at the four corners,
    ! calculate the fraction phi of the square where the function is larger than zero.
    
    ! Note: calculations below require double precision!! 
    
    implicit none
    
    ! In/output variables:
    real(wp), intent(IN)    :: f_NW, f_NE, f_SW, f_SE
    real(wp), intent(OUT)   :: phi
    
    ! Local variables:
    real(dp) :: f_NWp, f_NEp, f_SWp, f_SEp
    real(dp) :: aa,bb,cc,dd,x,f1,f2
    integer  :: scen

    real(wp), parameter :: ftol = 1e-4_dp
    
    ! The analytical solutions sometime give problems when one or more of the corner
    ! values is VERY close to zero; avoid this.
    if (f_NW == 0.0_dp) then
      f_NWp = ftol
    else if (f_NW > 0.0_dp) then
      f_NWp = MAX(  ftol, f_NW)
    else if (f_NW < 0.0_dp) then
      f_NWp = MIN( -ftol, f_NW)
    else
      f_NWp = f_NW
    end if
    if (f_NE == 0.0_dp) then
      f_NEp = ftol
    else if (f_NE > 0.0_dp) then
      f_NEp = MAX(  ftol, f_NE)
    else if (f_NE < 0.0_dp) then
      f_NEp = MIN( -ftol, f_NE)
    else
      f_NEp = f_NE
    end if
    if (f_SW == 0.0_dp) then
      f_SWp = ftol
    else if (f_SW > 0.0_dp) then
      f_SWp = MAX(  ftol, f_SW)
    else if (f_SW < 0.0_dp) then
      f_SWp = MIN( -ftol, f_SW)
    else
      f_SWp = f_SW
    end if
    if (f_SE == 0.0_dp) then
      f_SEp = ftol
    else if (f_SE > 0.0_dp) then
      f_SEp = MAX(  ftol, f_SE)
    else if (f_SE < 0.0_dp) then
      f_SEp = MIN( -ftol, f_SE)
    else
      f_SEp = f_SE
    end if
    
    if (f_NWp <= 0.0_dp .AND. f_NEp <= 0.0_dp .AND. f_SWp <= 0.0_dp .AND. f_SEp <= 0.0_dp) then
      ! All four corners are grounded.
      
      phi = 1.0_wp
      
    else if (f_NWp >= 0.0_dp .AND. f_NEp >= 0.0_dp .AND. f_SWp >= 0.0_dp .AND. f_SEp >= 0.0_dp) then
      ! All four corners are floating
      
      phi = 0.0_wp
      
    else
      ! At least one corner is grounded and at least one is floating;
      ! the grounding line must pass through this square!
      
      ! Only four "scenarios" exist (with rotational symmetries):
      ! 1: SW grounded, rest floating
      ! 2: SW floating, rest grounded
      ! 3: south grounded, north floating
      ! 4: SW & NE grounded, SE & NW floating
      ! Rotate the four-corner world until it matches one of these scenarios.
      call rotate_quad_until_match( f_NWp, f_NEp, f_SWp, f_SEp, scen)
    
      ! Calculate initial values of coefficients, and make correction
      ! for when d=0 (to avoid problems)

      aa  = f_SWp
      bb  = f_SEp - f_SWp
      cc  = f_NWp - f_SWp
      dd  = f_NEp + f_SWp - f_NWp - f_SEp

      ! Exception for when d=0
      if (ABS(dd) < ftol) then
        if (f_SWp > 0.0_dp) then
          f_SWp = f_SWp + 0.1_dp
        else
          f_SWp = f_SWp - 0.1_dp
        end if
        aa  = f_SWp
        bb  = f_SEp - f_SWp
        cc  = f_NWp - f_SWp
        dd  = f_NEp + f_SWp - f_NWp - f_SEp
      end if
        
      if (scen == 1) then
        ! 1: SW grounded, rest floating
        
        phi = ((bb*cc - aa*dd) * log(abs(1.0_dp - (aa*dd)/(bb*cc))) + aa*dd) / (dd**2)
         
      else if (scen == 2) then
        ! 2: SW floating, rest grounded
        ! Assign negative coefficients to calculate floating fraction,
        ! then get complement to obtain grounded fraction. 

        aa  = -(f_SWp)
        bb  = -(f_SEp - f_SWp)
        cc  = -(f_NWp - f_SWp)
        dd  = -(f_NEp + f_SWp - f_NWp - f_SEp)

        ! Exception for when d=0
        if (ABS(dd) < 1e-4_dp) then
          if (f_SWp > 0.0_dp) then
            f_SWp = f_SWp + 0.1_dp
          else
            f_SWp = f_SWp - 0.1_dp
          end if
          aa  = -(f_SWp)
          bb  = -(f_SEp - f_SWp)
          cc  = -(f_NWp - f_SWp)
          dd  = -(f_NEp + f_SWp - f_NWp - f_SEp)
        end if
        
        phi = 1.0_dp - ((bb*cc - aa*dd) * log(abs(1.0_dp - (aa*dd)/(bb*cc))) + aa*dd) / (dd**2)
        
      else if (scen == 3) then
        ! 3: south grounded, north floating
        
        ! Exception for when the GL runs parallel to the x-axis
        if (abs( 1.0_dp - f_NWp/f_NEp) < 1e-6_dp .and. abs( 1.0_dp - f_SWp/f_SEp) < 1e-6_dp) then
          
          phi = f_SWp / (f_SWp - f_NWp)
          
        else
            
          x   = 0.0_dp
          f1  = ((bb*cc - aa*dd) * log(abs(cc+dd*x)) - bb*dd*x) / (dd**2)
          x   = 1.0_dp
          f2  = ((bb*cc - aa*dd) * log(abs(cc+dd*x)) - bb*dd*x) / (dd**2)
          phi = f2-f1
                  
        end if
        
      else if (scen == 4) then
        ! 4: SW & NE grounded, SE & NW floating
        ! (recalculate coefficients here explicitly for two cases)

        ! SW corner
        aa  = f_SWp
        bb  = f_SEp - f_SWp
        cc  = f_NWp - f_SWp
        dd  = f_NEp + f_SWp - f_NWp - f_SEp
        phi = ((bb*cc - aa*dd) * log(abs(1.0_dp - (aa*dd)/(bb*cc))) + aa*dd) / (dd**2)
        
        ! NE corner
        call rotate_quad( f_NWp, f_NEp, f_SWp, f_SEp)
        call rotate_quad( f_NWp, f_NEp, f_SWp, f_SEp)
        aa  = f_SWp
        bb  = f_SEp - f_SWp
        cc  = f_NWp - f_SWp
        dd  = f_NEp + f_SWp - f_NWp - f_SEp
        phi = phi + ((bb*cc - aa*dd) * log(abs(1.0_dp - (aa*dd)/(bb*cc))) + aa*dd) / (dd**2)
        
      else
        write(io_unit_err,*) 'determine_grounded_fractions_CISM_quads - calc_fraction_above_zero - ERROR: unknown scenario [', scen, ']!'
        stop
      end if
      
    end if
    
    if (pi < -0.01_wp .OR. phi > 1.01_wp .OR. phi /= phi) then
      write(io_unit_err,*) 'calc_fraction_above_zero - ERROR: phi = ', phi
      write(io_unit_err,*) 'scen = ', scen
      write(io_unit_err,*) 'f = [', f_NWp, ',', f_NEp, ',', f_SWp, ',', f_SEp, ']'
      write(io_unit_err,*) 'aa = ', aa, ', bb = ', bb, ', cc = ', cc, ', dd = ', dd, ', f1 = ', f1, ',f2 = ', f2
      stop
    end if
    
    phi = MAX( 0.0_wp, MIN( 1.0_wp, phi))
    
    return

  end subroutine calc_fraction_above_zero

  subroutine rotate_quad_until_match( f_NW, f_NE, f_SW, f_SE, scen)
    ! Rotate the four corners until one of the four possible scenarios is found.
    ! 1: SW grounded, rest floating
    ! 2: SW floating, rest grounded
    ! 3: south grounded, north floating
    ! 4: SW & NE grounded, SE & NW floating
    
    implicit none
    
    ! In/output variables:
    real(dp), intent(INOUT) :: f_NW, f_NE, f_SW, f_SE
    integer,  intent(OUT)   :: scen
    
    ! Local variables:
    logical :: found_match
    integer :: nit
    
    found_match = .FALSE.
    scen        = 0
    nit         = 0
    
    do while (.not. found_match)
      
      nit = nit+1
      
      call rotate_quad( f_NW, f_NE, f_SW, f_SE)
      
      if     (f_SW < 0.0_wp .AND. f_SE > 0.0_wp .AND. f_NE > 0.0_wp .AND. f_NW > 0.0_wp) then
        ! 1: SW grounded, rest floating
        scen = 1
        found_match = .TRUE.
      else if (f_SW > 0.0_wp .AND. f_SE < 0.0_wp .AND. f_NE < 0.0_wp .AND. f_NW < 0.0_wp) then
        ! 2: SW floating, rest grounded
        scen = 2
        found_match = .TRUE.
      else if (f_SW < 0.0_wp .AND. f_SE < 0.0_wp .AND. f_NE > 0.0_wp .AND. f_NW > 0.0_wp) then
        ! 3: south grounded, north floating
        scen = 3
        found_match = .TRUE.
      else if (f_SW < 0.0_wp .AND. f_SE > 0.0_wp .AND. f_NE < 0.0_wp .AND. f_NW > 0.0_wp) then
        ! 4: SW & NE grounded, SE & NW floating
        scen = 4
        found_match = .TRUE.
      end if
      
      if (nit > 4) then
        write(io_unit_err,*) 
        write(io_unit_err,*) 'determine_grounded_fractions_CISM_quads - rotate_quad_until_match - ERROR: couldnt find matching scenario!'
        write(io_unit_err,*) 'f_SW, f_SE, f_NE, f_NW: ', f_SW, f_SE, f_NE, f_NW
        stop 
      end if
      
    end do
    
    return 

  end subroutine rotate_quad_until_match

  subroutine rotate_quad( f_NW, f_NE, f_SW, f_SE)
    ! Rotate the four corners anticlockwise by 90 degrees
    
    implicit none
    
    ! In/output variables:
    real(dp), intent(INOUT) :: f_NW, f_NE, f_SW, f_SE
    
    ! Local variables:
    real(dp) :: fvals(4)
    
    fvals = [f_NW,f_NE,f_SE,f_SW]
    f_NW = fvals( 2)
    f_NE = fvals( 3)
    f_SE = fvals( 4)
    f_SW = fvals( 1)
    
    return 

  end subroutine rotate_quad


    elemental function cdf(x,mu,sigma,inv) result(F)
        ! Solve for cumulative probability below (cdf)
        ! or above (cdf(inv=TRUE)) the value x
        ! given a normal distribution N(mu,sigma)

        ! See equation for F(x) and Q(x) in, e.g.:
        ! https://en.wikipedia.org/wiki/Normal_distribution

        implicit none

        real(wp), intent(IN)  :: x
        real(wp), intent(IN)  :: mu
        real(wp), intent(IN)  :: sigma
        logical,  intent(IN), optional :: inv
        real(wp) :: F
        
        real(wp), parameter :: sqrt2 = sqrt(2.0_wp)

        ! Calculate CDF at value x
        F = 0.5_wp*( 1.0_wp + error_function((x-mu)/(sqrt2*sigma)) )

        if (present(inv)) then 
            if (inv) then 
                ! Calculate inverse CDF at value x 
                F = 1.0_wp - F 
            end if 
        end if 

        return

    end function cdf

    elemental function error_function(X) result(ERR)
        ! Purpose: Compute error function erf(x)
        ! Input:   x   --- Argument of erf(x)
        ! Output:  ERR --- erf(x)
        
        ! Note: also separately defined in thermodynamics.f90

        implicit none 

        real(wp), intent(IN)  :: X
        real(wp) :: ERR
        
        ! Local variables:
        real(wp)              :: EPS
        real(wp)              :: X2
        real(wp)              :: ER
        real(wp)              :: R
        real(wp)              :: C0
        integer                 :: k
        
        EPS = 1.0e-15
        X2  = X * X
        if (abs(X) < 3.5) then
            ER = 1.0
            R  = 1.0
            do k = 1, 50
                R  = R * X2 / (real(k, wp) + 0.5)
                ER = ER+R
                if(abs(R) < abs(ER) * EPS) then
                    C0  = 2.0 / sqrt(pi) * X * exp(-X2)
                    ERR = C0 * ER
                    EXIT
                end if
            end do
        else
            ER = 1.0
            R  = 1.0
            do k = 1, 12
                R  = -R * (real(k, wp) - 0.5) / X2
                ER = ER + R
                C0  = EXP(-X2) / (abs(X) * sqrt(pi))
                ERR = 1.0 - C0 * ER
                if(X < 0.0) ERR = -ERR
            end do
        end if

        return

    end function error_function

end module topography


