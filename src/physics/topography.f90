module topography 

    use yelmo_defs, only : wp, dp, prec, io_unit_err, sec_year, pi, T0, g, rho_ice, rho_sw, rho_w 

    implicit none 

    private  

    public :: extend_floating_slab
    public :: remove_fractional_ice
    public :: remove_icebergs 

    public :: calc_ice_fraction

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

contains 

    subroutine extend_floating_slab(H_ice,f_grnd,H_slab,n_ext)
        ! Extend ice field so that there is always 
        ! floating ice next to grounded marine margins
        ! Extended ice should be very thin, will 
        ! be assigned value H_slab. Slab will be extended
        ! n_ext points away from marine margin

        implicit none 

        real(wp), intent(INOUT) :: H_ice(:,:) 
        real(wp), intent(IN)    :: f_grnd(:,:) 
        real(wp), intent(IN)    :: H_slab       ! Typically 1 or 0.1 m. 
        integer,  intent(IN)    :: n_ext        ! Number of points to extend slab
        
        ! Local variables 
        integer :: i, j, nx, ny, iter 
        integer :: im1, ip1, jm1, jp1
        logical :: is_marine 

        logical  :: ms4(4)
        real(wp) :: Hi4(4) 
        real(wp) :: fg4(4)
        
        logical,  allocatable :: mask_slab(:,:)
        real(wp), allocatable :: H_new(:,:) 

        nx = size(H_ice,1) 
        ny = size(H_ice,2) 

        allocate(mask_slab(nx,ny)) 
        allocate(H_new(nx,ny)) 

        mask_slab = .FALSE. 
        H_new     = H_ice 

        do iter = 1, n_ext

            do j = 1, ny 
            do i = 1, nx 

                ! Get neighbor indices
                im1 = max(i-1,1) 
                ip1 = min(i+1,nx) 
                jm1 = max(j-1,1) 
                jp1 = min(j+1,ny) 
                
                if ( f_grnd(i,j) .eq. 0.0 .and. H_ice(i,j) .eq. 0.0 ) then 
                    ! Floating ice-free ocean point
                    
                    ! Get neighbor values in convenient arrays
                    fg4 = [f_grnd(im1,j),f_grnd(ip1,j),f_grnd(i,jm1),f_grnd(i,jp1)]
                    Hi4 = [H_ice(im1,j),H_ice(ip1,j),H_ice(i,jm1),H_ice(i,jp1)]
                    ms4 = [mask_slab(im1,j),mask_slab(ip1,j),mask_slab(i,jm1),mask_slab(i,jp1)]

                    if ( (count(fg4 .gt. 0.0 .and. Hi4 .gt. 0.0) .gt. 0) .or. &
                         (count(ms4) .gt. 0) ) then 
                        ! At least one neighbors is either a grounded point
                        ! or an extended slab point - make this point extended slab.

                        H_new(i,j)     = H_slab 
                        
                    end if

                end if 

            end do 
            end do 

            ! Update H_ice to current array 
            H_ice = H_new 

            ! Update mask_slab
            where(H_ice .eq. H_slab) 
                mask_slab = .TRUE. 
            elsewhere
                mask_slab = .FALSE.
            end where

        end do 

        return

    end subroutine extend_floating_slab

    subroutine remove_fractional_ice(H_ice,f_ice)
        ! Eliminate fractional ice covered points that only 
        ! have fractional ice neighbors. 

        implicit none 

        real(wp), intent(INOUT) :: H_ice(:,:) 
        real(wp), intent(INOUT) :: f_ice(:,:) 
        
        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: im1, ip1, jm1, jp1 
        real(wp), allocatable :: H_new(:,:) 

        nx = size(H_ice,1) 
        ny = size(H_ice,2) 

        allocate(H_new(nx,ny)) 

        H_new = H_ice 

        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices
            im1 = max(i-1,1) 
            ip1 = min(i+1,nx) 
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny) 
            
            if (f_ice(i,j) .gt. 0.0 .and. f_ice(i,j) .lt. 1.0) then 
                ! Fractional ice-covered point 

                if ( count([f_ice(im1,j),f_ice(ip1,j), &
                        f_ice(i,jm1),f_ice(i,jp1)] .eq. 1.0) .eq. 0) then 
                    ! No fully ice-covered neighbors available.
                    ! Point should be removed. 

                    H_new(i,j) = 0.0_wp 

                end if

            end if

        end do 
        end do

        ! Update ice thickness 
        H_ice = H_new 

        return

    end subroutine remove_fractional_ice

    subroutine remove_icebergs(H_ice)

        implicit none 

        real(wp), intent(INOUT) :: H_ice(:,:) 

        return

    end subroutine remove_icebergs

    subroutine find_connected_mask(mask,mask_ref,mask_now)
        ! Brute-force routine to find all points 
        ! that touch or are connected with other points in a mask.
        ! Here use to find any floating ice points that are
        ! not connected to grounded ice or ice-free land. 

        ! AJR: TO DO, this routine is not ready!!!

        implicit none 

        logical, intent(INOUT) :: mask(:,:)         ! Connected points of interest
        logical, intent(IN)    :: mask_ref(:,:)     ! Points to be connected to
        logical, intent(IN)    :: mask_now(:,:)     ! Points of interest

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
                im1 = max(i-1,1)
                ip1 = min(i+1,nx)
                jm1 = max(j-1,1)
                jp1 = min(j+1,ny)
                
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
    
    subroutine calc_ice_fraction_new(f_ice,H_ice,z_bed,z_sl,flt_subgrid)
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
            call calc_H_grnd(H_grnd,H_ice,f_ice,z_bed,z_sl,use_f_ice=.FALSE.)
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
            im1 = max(i-1,1) 
            ip1 = min(i+1,nx) 
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny) 
            
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
                im1 = max(i-1,1) 
                ip1 = min(i+1,nx) 
                jm1 = max(j-1,1) 
                jp1 = min(j+1,ny) 
                
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
    
        subroutine calc_ice_fraction(f_ice,H_ice,z_bed,z_sl,flt_subgrid)
        ! Determine the area fraction of a grid cell
        ! that is ice-covered. Assume that marginal points
        ! have equal thickness to inland neighbors 

        implicit none 

        real(wp), intent(OUT) :: f_ice(:,:)             ! [--] Ice covered fraction (aa-nodes)
        real(wp), intent(IN)  :: H_ice(:,:)             ! [m] Ice thickness on standard grid (aa-nodes)
        real(wp), intent(IN)  :: z_bed(:,:)             ! [m] Bedrock elevation
        real(wp), intent(IN)  :: z_sl(:,:)              ! [m] Sea-level elevation
        logical, optional     :: flt_subgrid            ! Option to allow fractions for floating ice margins             
        
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
            call calc_H_grnd(H_grnd,H_ice,f_ice,z_bed,z_sl,use_f_ice=.FALSE.)
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

            do j = 1, ny
            do i = 1, nx 

                ! Get neighbor indices
                im1 = max(i-1,1) 
                ip1 = min(i+1,nx) 
                jm1 = max(j-1,1) 
                jp1 = min(j+1,ny) 
                
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

            ! Determine ice fractional cover for margin points 

            do j = 1, ny
            do i = 1, nx 

                ! Get neighbor indices
                im1 = max(i-1,1) 
                ip1 = min(i+1,nx) 
                jm1 = max(j-1,1) 
                jp1 = min(j+1,ny) 
                
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

        end if 

        return 

    end subroutine calc_ice_fraction
    
    elemental subroutine calc_z_srf(z_srf,H_ice,f_ice,H_grnd,z_bed,z_sl)
        ! Calculate surface elevation

        implicit none 

        real(prec), intent(INOUT) :: z_srf
        real(prec), intent(IN)    :: H_ice
        real(prec), intent(IN)    :: f_ice
        real(prec), intent(IN)    :: H_grnd
        real(prec), intent(IN)    :: z_bed
        real(prec), intent(IN)    :: z_sl

        ! Local variables 
        real(prec) :: rho_ice_sw
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

    elemental subroutine calc_z_srf_max(z_srf,H_ice,f_ice,z_bed,z_sl)
        ! Calculate surface elevation
        ! Adapted from Pattyn (2017), Eq. 1
        
        implicit none 

        real(prec), intent(INOUT) :: z_srf 
        real(prec), intent(IN)    :: H_ice
        real(prec), intent(IN)    :: f_ice
        real(prec), intent(IN)    :: z_bed
        real(prec), intent(IN)    :: z_sl

        ! Local variables
        integer :: i, j, nx, ny 
        real(prec) :: rho_ice_sw
        real(prec) :: H_eff

        rho_ice_sw = rho_ice/rho_sw ! Ratio of density of ice to seawater [--]
        
        ! Get effective ice thickness
        call calc_H_eff(H_eff,H_ice,f_ice,set_frac_zero=.TRUE.)

        ! Initially calculate surface elevation everywhere 
        z_srf = max(z_bed + H_eff, z_sl + (1.0-rho_ice_sw)*H_eff)
        
        return 

    end subroutine calc_z_srf_max

    subroutine calc_z_srf_gl_subgrid_area(z_srf,f_grnd,H_ice,f_ice,z_bed,z_sl,gl_sep_nx)
        ! Interpolate variables at grounding line to subgrid level to 
        ! calculate the average z_srf value for the aa-node cell

        implicit none
        
        real(prec), intent(OUT) :: z_srf(:,:)       ! aa-nodes 
        real(prec), intent(IN)  :: f_grnd(:,:)      ! aa-nodes
        real(prec), intent(IN)  :: H_ice(:,:)
        real(prec), intent(IN)  :: f_ice(:,:)
        real(prec), intent(IN)  :: z_bed(:,:)
        real(prec), intent(IN)  :: z_sl(:,:)
        integer,    intent(IN)  :: gl_sep_nx        ! Number of interpolation points per side (nx*nx)

        ! Local variables
        integer  :: i, j, nx, ny
        real(wp) :: v1, v2, v3, v4 
        integer  :: im1, ip1, jm1, jp1
        real(wp) :: H_eff 
        real(wp) :: f_grnd_neighb(4) 
        logical  :: is_grline 

        real(prec), allocatable :: z_srf_int(:,:) 
        real(prec), allocatable :: H_ice_int(:,:)
        real(prec), allocatable :: f_ice_int(:,:)  
        real(prec), allocatable :: z_bed_int(:,:) 
        real(prec), allocatable :: z_sl_int(:,:) 
        
        nx = size(z_srf,1)
        ny = size(z_srf,2) 

        ! Allocate the subgrid arrays 
        allocate(z_srf_int(gl_sep_nx,gl_sep_nx))
        allocate(H_ice_int(gl_sep_nx,gl_sep_nx))
        allocate(f_ice_int(gl_sep_nx,gl_sep_nx))
        allocate(z_bed_int(gl_sep_nx,gl_sep_nx))
        allocate(z_sl_int(gl_sep_nx,gl_sep_nx))
        
        ! ajr: assume f_ice_int=1 everywhere this is used for now. 
        ! Needs to be fixed in the future potentially. 
        f_ice_int = 1.0_wp 

        write(*,*) "calc_z_srf_gl_subgrid_area:: routine not ready for f_ice values. Fix!"
        stop 
        
        ! Calculate the surface elevation based on whole grid values,
        ! except at the grounding line which is treated with subgrid interpolations. 
        do j = 1, ny 
        do i = 1, nx

            ! Get neighbor indices
            im1 = max(i-1,1) 
            ip1 = min(i+1,nx) 
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny) 
            
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

                ! Calculate values at corners (ab-nodes) and interpolate
                
                ! == H_ice == 
                v1 = 0.25_prec*(H_ice(i,j) + H_ice(ip1,j) + H_ice(ip1,jp1) + H_ice(i,jp1))
                v2 = 0.25_prec*(H_ice(i,j) + H_ice(im1,j) + H_ice(im1,jp1) + H_ice(i,jp1))
                v3 = 0.25_prec*(H_ice(i,j) + H_ice(im1,j) + H_ice(im1,jm1) + H_ice(i,jm1))
                v4 = 0.25_prec*(H_ice(i,j) + H_ice(ip1,j) + H_ice(ip1,jm1) + H_ice(i,jm1))
                call calc_subgrid_array(H_ice_int,v1,v2,v3,v4,gl_sep_nx)
                
                ! == z_bed == 
                v1 = 0.25_prec*(z_bed(i,j) + z_bed(ip1,j) + z_bed(ip1,jp1) + z_bed(i,jp1))
                v2 = 0.25_prec*(z_bed(i,j) + z_bed(im1,j) + z_bed(im1,jp1) + z_bed(i,jp1))
                v3 = 0.25_prec*(z_bed(i,j) + z_bed(im1,j) + z_bed(im1,jm1) + z_bed(i,jm1))
                v4 = 0.25_prec*(z_bed(i,j) + z_bed(ip1,j) + z_bed(ip1,jm1) + z_bed(i,jm1))
                call calc_subgrid_array(z_bed_int,v1,v2,v3,v4,gl_sep_nx)
                
                ! == z_sl == 
                v1 = 0.25_prec*(z_sl(i,j) + z_sl(ip1,j) + z_sl(ip1,jp1) + z_sl(i,jp1))
                v2 = 0.25_prec*(z_sl(i,j) + z_sl(im1,j) + z_sl(im1,jp1) + z_sl(i,jp1))
                v3 = 0.25_prec*(z_sl(i,j) + z_sl(im1,j) + z_sl(im1,jm1) + z_sl(i,jm1))
                v4 = 0.25_prec*(z_sl(i,j) + z_sl(ip1,j) + z_sl(ip1,jm1) + z_sl(i,jm1))
                call calc_subgrid_array(z_sl_int,v1,v2,v3,v4,gl_sep_nx)
                
                ! Calculate subgrid surface elevations
                call calc_z_srf_max(z_srf_int,H_ice_int,f_ice_int,z_bed_int,z_sl_int)

                ! Calculate full grid z_srf value as the mean of subgrid values 
                z_srf(i,j) = sum(z_srf_int) / real(gl_sep_nx*gl_sep_nx,prec)

            end if 

        end do 
        end do 

        return
        
    end subroutine calc_z_srf_gl_subgrid_area

    elemental subroutine calc_H_eff(H_eff,H_ice,f_ice,set_frac_zero)
        ! Calculate ice-thickness, scaled at margins to actual thickness
        ! but as if it covered the whole grid cell.
        
        implicit none

        real(prec), intent(OUT) :: H_eff 
        real(prec), intent(IN)  :: H_ice 
        real(prec), intent(IN)  :: f_ice 
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

    elemental subroutine calc_H_grnd(H_grnd,H_ice,f_ice,z_bed,z_sl,use_f_ice)
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
        !H_grnd = H_eff - rho_sw_ice*max(z_sl-z_bed,0.0_prec)

        ! ajr: testing. This ensures that ice-free ground above sea level
        ! also has H_grnd > 0.
        if (z_sl-z_bed .gt. 0.0) then 
            ! Grounded below sea level, diagnose overburden minus water thickness
            H_grnd = H_eff - rho_sw_ice*(z_sl-z_bed)
        else
            ! Grounded above sea level, simply sum elevation above sea level and ice thickness
            H_grnd = H_eff + (z_bed-z_sl)
        end if 

        ! ajr: to test somewhere eventually, more closely follows Gladstone et al (2010), Leguy etl (2021)
        !H_grnd = -( (z_sl-z_bed) - rho_ice_sw*H_eff )

        return 

    end subroutine calc_H_grnd

    elemental subroutine calc_H_af(H_af,H_ice,f_ice,z_bed,z_sl,use_f_ice)
        ! Calculate ice thickness above flotation, H_af

        implicit none 

        real(wp), intent(INOUT) :: H_af
        real(wp), intent(IN)    :: H_ice
        real(wp), intent(IN)    :: f_ice
        real(wp), intent(IN)    :: z_bed
        real(wp), intent(IN)    :: z_sl 
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

    subroutine calc_f_grnd_subgrid_area_aa(f_grnd,H_grnd,gl_sep_nx)
        ! Use H_grnd to determined grounded area fraction of grid point.

        implicit none
        
        real(prec), intent(OUT) :: f_grnd(:,:)      ! aa-nodes 
        real(prec), intent(IN)  :: H_grnd(:,:)      ! aa-nodes
        integer,    intent(IN)  :: gl_sep_nx        ! Number of interpolation points per side (nx*nx)

        ! Local variables
        integer    :: i, j, nx, ny
        real(prec) :: Hg_1, Hg_2, Hg_3, Hg_4, Hg_mid  
        integer    :: im1, ip1, jm1, jp1 

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
            im1 = max(i-1,1) 
            ip1 = min(i+1,nx) 
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny) 
                
            ! Calculate Hg at corners (ab-nodes)
            Hg_1 = 0.25_prec*(H_grnd(i,j) + H_grnd(ip1,j) + H_grnd(ip1,jp1) + H_grnd(i,jp1))
            Hg_2 = 0.25_prec*(H_grnd(i,j) + H_grnd(im1,j) + H_grnd(im1,jp1) + H_grnd(i,jp1))
            Hg_3 = 0.25_prec*(H_grnd(i,j) + H_grnd(im1,j) + H_grnd(im1,jm1) + H_grnd(i,jm1))
            Hg_4 = 0.25_prec*(H_grnd(i,j) + H_grnd(ip1,j) + H_grnd(ip1,jm1) + H_grnd(i,jm1))
            
            if (max(Hg_1,Hg_2,Hg_3,Hg_4) .ge. 0.0 .and. min(Hg_1,Hg_2,Hg_3,Hg_4) .lt. 0.0) then 
                ! Point contains grounding line, get grounded area  
                
                call calc_grounded_fraction_cell(f_grnd(i,j),Hg_1,Hg_2,Hg_3,Hg_4,gl_sep_nx)

            end if 

        end do 
        end do 

        return
        
    end subroutine calc_f_grnd_subgrid_area_aa
    
    subroutine calc_f_grnd_subgrid_area(f_grnd,f_grnd_acx,f_grnd_acy,H_grnd,gl_sep_nx)
        ! Use H_grnd to determined grounded area fraction of grid point.

        implicit none
        
        real(wp), intent(OUT) :: f_grnd(:,:)        ! aa-nodes 
        real(wp), intent(OUT) :: f_grnd_acx(:,:)    ! ac-nodes
        real(wp), intent(OUT) :: f_grnd_acy(:,:)    ! ac-nodes
        real(wp), intent(IN)  :: H_grnd(:,:)        ! aa-nodes
        integer,  intent(IN)  :: gl_sep_nx          ! Number of interpolation points per side (nx*nx)

        ! Local variables
        integer  :: i, j, nx, ny
        real(wp) :: Hg_1, Hg_2, Hg_3, Hg_4
        real(wp) :: Hg_min, Hg_max  
        integer  :: im1, ip1, jm1, jp1 

        !integer, parameter :: nx_interp = 15

        nx = size(H_grnd,1)
        ny = size(H_grnd,2) 

        ! Initialize all masks to zero (fully floating) to start
        f_grnd     = 0.0_wp 
        f_grnd_acx = 0.0_wp 
        f_grnd_acy = 0.0_wp 

        ! Find grounding line cells and determine fraction 
        do j = 1, ny 
        do i = 1, nx

            ! Get neighbor indices
            im1 = max(i-1,1) 
            ip1 = min(i+1,nx) 
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny) 
            
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
                
                call calc_grounded_fraction_cell(f_grnd(i,j),Hg_1,Hg_2,Hg_3,Hg_4,gl_sep_nx)

            else if (Hg_max .ge. 0.0 .and. Hg_min .ge. 0.0) then 
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
                
                call calc_grounded_fraction_cell(f_grnd_acx(i,j),Hg_1,Hg_2,Hg_3,Hg_4,gl_sep_nx)

            else if (Hg_max .ge. 0.0 .and. Hg_min .ge. 0.0) then
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
                
                call calc_grounded_fraction_cell(f_grnd_acy(i,j),Hg_1,Hg_2,Hg_3,Hg_4,gl_sep_nx)

            else if (Hg_max .ge. 0.0 .and. Hg_min .ge. 0.0) then 
                ! Purely grounded point 
                    
                f_grnd_acy(i,j) = 1.0_wp 
                
            end if 

        end do 
        end do 


if (.TRUE.) then 
    ! Replace subgrid acx/acy estimates with linear average to ac-nodes 

        ! acx-nodes 
        do j = 1, ny 
        do i = 1, nx-1
            f_grnd_acx(i,j) = 0.5_prec*(f_grnd(i,j) + f_grnd(i+1,j))
        end do 
        end do
        f_grnd_acx(nx,:) = f_grnd_acx(nx-1,:) 

        ! acy-nodes 
        do j = 1, ny-1 
        do i = 1, nx
            f_grnd_acy(i,j) = 0.5_prec*(f_grnd(i,j) + f_grnd(i,j+1))
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

        real(prec), intent(OUT) :: f_grnd(:,:)
        real(prec), intent(OUT) :: f_grnd_x(:,:)
        real(prec), intent(OUT) :: f_grnd_y(:,:)
        real(prec), intent(IN)  :: H_grnd(:,:)

        ! Local variables  
        integer :: i, j, nx, ny 
        real(prec) :: H_grnd_1, H_grnd_2

        nx = size(f_grnd,1)
        ny = size(f_grnd,2)

        ! Central aa-node
        f_grnd = 1.0
        where (H_grnd < 0.0) f_grnd = 0.0
        
        ! x-direction, ac-node
        f_grnd_x = 1.0
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

        ! y-direction, ac-node
        f_grnd_y = 1.0
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

        ! Set boundary points equal to neighbor for aesthetics 
        f_grnd_x(nx,:) = f_grnd_x(nx-1,:) 
        f_grnd_y(:,ny) = f_grnd_y(:,ny-1) 
        
        return 

    end subroutine calc_f_grnd_subgrid_linear
    
    subroutine calc_f_grnd_pinning_points(f_grnd,H_ice,f_ice,z_bed,z_bed_sd,z_sl)
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

        return

    end subroutine calc_f_grnd_pinning_points

    subroutine remove_englacial_lakes(H_ice,z_bed,z_srf,z_sl)
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

    subroutine calc_distance_to_ice_margin(dist_mrgn,f_ice,dx)
        ! Calculate distance to the ice margin
        
        ! Note: this subroutine is a wrapper that calls the
        ! grounding-line distance routine, since the algorithm
        ! works the same way. Simply substitute f_ice for f_grnd. 

        implicit none 

        real(wp), intent(OUT) :: dist_mrgn(:,:) ! [km] Distance to grounding line
        real(wp), intent(IN)  :: f_ice(:,:)     ! [1]  Fraction of grid-cell ice coverage 
        real(wp), intent(IN)  :: dx             ! [m]  Grid resolution (assume dy=dx)

        call calc_distance_to_grounding_line(dist_mrgn,f_ice,dx)

        return 

    end subroutine calc_distance_to_ice_margin
    
    subroutine calc_distance_to_grounding_line(dist_gl,f_grnd,dx)
        ! Calculate distance to the grounding line 
        
        implicit none 

        real(wp), intent(OUT) :: dist_gl(:,:)   ! [km] Distance to grounding line
        real(wp), intent(IN)  :: f_grnd(:,:)    ! [1]  Grounded grid-cell fraction 
        real(wp), intent(IN)  :: dx             ! [m]  Grid resolution (assume dy=dx)

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

        do j = 1, ny 
        do i = 1, nx

            ! Get neighbor indices
            im1 = max(1,i-1)
            ip1 = min(nx,i+1)
            jm1 = max(1,j-1)
            jp1 = min(ny,j+1)

            ! Grounded point or partially floating point with floating neighbors
            if (f_grnd(i,j) .gt. 0.0 .and. &
                (f_grnd(im1,j) .eq. 0.0 .or. f_grnd(ip1,j) .eq. 0.0 .or. &
                 f_grnd(i,jm1) .eq. 0.0 .or. f_grnd(i,jp1) .eq. 0.0) ) then 
                
                dist_gl(i,j)  = 0.0_wp 

            end if 

        end do 
        end do 

        ! 2. Next, determine distances to grounding line ======================
        
        do q = 1, iter_max  
            ! Iterate distance of one neighbor at a time until grid is filled in

            dist_gl_ref = dist_gl

            do j = 1, ny 
            do i = 1, nx

                if (dist_gl(i,j) .eq. dist_max) then 
                    ! Distance needs to be determined for this point 

                    ! Get neighbor indices
                    im1 = max(1, i-1)
                    ip1 = min(nx,i+1)
                    jm1 = max(1, j-1)
                    jp1 = min(ny,j+1)

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

    subroutine calc_bmb_total(bmb,bmb_grnd,bmb_shlf,H_ice,H_grnd,f_grnd,bmb_gl_method,diffuse_bmb_shlf)

        implicit none 

        real(prec),       intent(OUT) :: bmb(:,:) 
        real(prec),       intent(IN)  :: bmb_grnd(:,:) 
        real(prec),       intent(IN)  :: bmb_shlf(:,:) 
        real(prec),       intent(IN)  :: H_ice(:,:)
        real(prec),       intent(IN)  :: H_grnd(:,:)
        real(prec),       intent(IN)  :: f_grnd(:,:) 
        character(len=*), intent(IN)  :: bmb_gl_method 
        logical,          intent(IN)  :: diffuse_bmb_shlf 

        ! Local variables
        integer    :: i, j, nx, ny
        integer    :: n_float  
        real(prec) :: bmb_shlf_now 

        nx = size(bmb,1)
        ny = size(bmb,2) 

        ! Combine floating and grounded parts into one field =========================

        ! Initialize bmb to zero everywhere to start and apply bmb_grnd to grounded ice 
        bmb = 0.0_wp 
        where(f_grnd .eq. 1.0) bmb = bmb_grnd 

        ! Apply the floating basal mass balance according 
        ! to different subgridding options at the grounding line
        ! (following the notation of Leguy et al., 2021 - see Fig. 3)
        select case(bmb_gl_method)

            case("fcmp")
                ! Flotation criterion melt parameterization 
                ! Apply full bmb_shlf value where flotation criterion is met 

                where(H_grnd .le. 0.0_wp)
                    
                    bmb = bmb_shlf

                end where 

            case("fmp")
                ! Full melt parameterization
                ! Apply full bmb_shlf value to any cell that is at least
                ! partially floating. Perhaps unrealistic.

                where(f_grnd .lt. 1.0_wp)

                    bmb = bmb_shlf 

                end where 

            case("pmp")
                ! Partial melt parameterization
                ! Apply bmb_shlf to floating fraction of cell 

                where(f_grnd .lt. 1.0_wp)

                    bmb = f_grnd*bmb_grnd + (1.0_wp-f_grnd)*bmb_shlf 

                end where 

            case("nmp")
                ! No melt parameterization
                ! Apply bmb_shlf only where fully floating diagnosed

                where(f_grnd .eq. 0.0_wp)

                    bmb = bmb_shlf 

                end where 

        end select

        ! For aesthetics, also make sure that bmb is zero on ice-free land
        where (H_grnd .gt. 0.0_wp .and. H_ice .eq. 0.0_wp) bmb = 0.0_wp 

        if (diffuse_bmb_shlf) then 
            ! Allow marine melt (bmb_shlf) to permeate inland at the grounding line,
            ! to induce more effective retreat in warm periods 
            ! Note: this method is not recommended, especially
            ! when a subgrid melting parameterization is used above! 
            
            do j = 1, ny
            do i = 1, nx

                if (f_grnd(i,j) .eq. 1.0) then 
                    ! Grounded point, look for floating neighbors 

                    if (.FALSE.) then
                        ! 9-neighbor method

                        n_float = count(f_grnd(i-1:i+1,j-1:j+1) .lt. 1.0)

                        if (n_float .gt. 0) then 
                            ! bmb_shelf is the mean of the neighbours
                            bmb_shlf_now = sum(bmb_shlf(i-1:i+1,j-1:j+1),mask=f_grnd(i-1:i+1,j-1:j+1) .lt. 1.0) / real(n_float,prec)
                            bmb(i,j)     = (1.0-n_float/9.0)*bmb_grnd(i,j) + (n_float/9.0)*bmb_shlf_now

                        end if

                    else 
                        ! 5-neighbor method 

                        n_float = count([f_grnd(i-1,j),f_grnd(i+1,j),f_grnd(i,j-1),f_grnd(i,j+1)].lt. 1.0)

                        if (n_float .gt. 0) then
                            ! Floating points exist 
                            bmb_shlf_now = sum([bmb_shlf(i-1,j),bmb_shlf(i+1,j),bmb_shlf(i,j-1),bmb_shlf(i,j+1)], &
                                            mask=[f_grnd(i-1,j),f_grnd(i+1,j),f_grnd(i,j-1),f_grnd(i,j+1)] .lt. 1.0) / real(n_float,prec)
                            bmb(i,j)     = (1.0-n_float/5.0)*bmb_grnd(i,j) + (n_float/5.0)*bmb_shlf_now
                        end if

                    end if


                end if

            end do
            end do

        end if 

        return 

    end subroutine calc_bmb_total

    subroutine calc_fmb_total(fmb,fmb_shlf,bmb_shlf,H_ice,H_grnd,f_ice,fmb_method,fmb_scale,dx)

        implicit none 

        real(wp), intent(OUT) :: fmb(:,:) 
        real(wp), intent(IN)  :: fmb_shlf(:,:) 
        real(wp), intent(IN)  :: bmb_shlf(:,:) 
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: H_grnd(:,:) 
        real(wp), intent(IN)  :: f_ice(:,:)
        integer,  intent(IN)  :: fmb_method 
        real(wp), intent(IN)  :: fmb_scale
        real(wp), intent(IN)  :: dx

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

            case(1) 
                ! Calculate fmb as proportional to local bmb_shlf value and 
                ! scaled to the area of the grid cell itself where it will be applied

                ! Note: fmb_scale=10 suggested by DeConto and Pollard (2016, nat) based
                ! on plume modeling work of Slater et al. (2015, grl)

                do j = 1, ny 
                do i = 1, nx 

                    ! Get neighbor indices
                    im1 = max(i-1,1) 
                    ip1 = min(i+1,nx) 
                    jm1 = max(j-1,1) 
                    jp1 = min(j+1,ny) 
                    
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

                        ! Also calculate the mean bmb_shlf value for the ice-free neighbors 
                        bmb_eff = sum([bmb_shlf(im1,j),bmb_shlf(ip1,j),bmb_shlf(i,jm1),bmb_shlf(i,jp1)], &
                                        mask=mask) / real(n_margin,wp)

                        ! Finally calculate the effective front mass balance rate 

                        fmb(i,j) = bmb_eff*(area_flt/area_tot)*fmb_scale

                    else 
                        ! Set front mass balance equal to zero 

                        fmb(i,j) = 0.0_wp 

                    end if 

                end do 
                end do 

            case DEFAULT 

                write(*,*) "calc_fmb_total:: Error: fmb_method not recongized."
                write(*,*) "fmb_method = ", fmb_method 
                stop 

        end select

        return 

    end subroutine calc_fmb_total

    subroutine calc_subgrid_array(vint,v1,v2,v3,v4,nx)
        ! Given the four corners of a cell in quadrants 1,2,3,4,
        ! calculate the subgrid values via linear interpolation

        implicit none 

        real(prec), intent(OUT) :: vint(:,:)  
        real(prec), intent(IN)  :: v1,v2,v3,v4
        integer,    intent(IN)  :: nx                    ! Number of interpolation points 

        ! Local variables 
        integer :: i, j 
        real(prec) :: x(nx), y(nx) 

        ! Populate x,y axes for interpolation points (between 0 and 1)
        do i = 1, nx 
            x(i) = 0.0 + real(i-1)/real(nx-1)
        end do 
        y = x 
        
        ! Calculate interpolated value      
        vint = 0.0 
        do i = 1, nx 
        do j = 1, nx 

            vint(i,j) = interp_bilin_pt(v1,v2,v3,v4,x(i),y(j))

        end do 
        end do 

        return 

    end subroutine calc_subgrid_array

    subroutine calc_grounded_fraction_cell(f_g,Hg_1,Hg_2,Hg_3,Hg_4,nx)
        ! Given the four corners of a cell in quadrants 1,2,3,4,
        ! calculate the grounded fraction (ie area with Hg>0)
        !
        ! Convention:
        !   
        !  Hg_2---Hg_1
        !    |     |
        !    |     |
        !  Hg_3---Hg_4

        implicit none 

        real(wp), intent(OUT) :: f_g 
        real(wp), intent(IN)  :: Hg_1, Hg_2, Hg_3, Hg_4
        integer,    intent(IN)  :: nx                    ! Number of interpolation points 

        ! Local variables 
        integer  :: i, j 
        real(wp) :: x(nx), y(nx) 
        real(wp) :: Hg_int(nx,nx)  

        ! Populate x,y axes for interpolation points (between 0 and 1)
        do i = 1, nx 
            x(i) = 0.0 + real(i-1)/real(nx-1)
        end do 
        y = x 

        ! Perform interpolation of Hg onto fine grid
        Hg_int = 0.0 
        do i = 1, nx 
        do j = 1, nx 

            Hg_int(i,j) = interp_bilin_pt(Hg_1,Hg_2,Hg_3,Hg_4,x(i),y(j))

        end do 
        end do 

        ! Calculate weighted fraction (assume all points have equal weight)
        f_g = real(count(Hg_int .ge. 0.0),wp) / real(nx*nx,wp)

        return 

    end subroutine calc_grounded_fraction_cell

    function interp_bilin_pt(z1,z2,z3,z4,xout,yout) result(zout)
        ! Interpolate a point given four neighbors at corners of square (0:1,0:1)
        ! z2    z1
        !    x,y
        ! z3    z4 
        ! 

        implicit none 

        real(prec), intent(IN) :: z1, z2, z3, z4 
        real(prec), intent(IN) :: xout, yout 
        real(prec) :: zout 

        ! Local variables 
        real(prec) :: x0, x1, y0, y1 
        real(prec) :: alpha1, alpha2, p0, p1 

        x0 = 0.0 
        x1 = 1.0 
        y0 = 0.0 
        y1 = 1.0 

        alpha1  = (xout - x0) / (x1-x0)
        p0      = z3 + alpha1*(z4-z3)
        p1      = z2 + alpha1*(z1-z2)
            
        alpha2  = (yout - y0) / (y1-y0)
        zout    = p0 + alpha2*(p1-p0)

        return 

    end function interp_bilin_pt

    function distance_to_margin(H_ice,dx) result(dist)
        ! ajr: very slow and obsolete routine - do not use! 
        ! instead, use calc_distance_to_margin above

        implicit none 

        real(prec), intent(IN) :: H_ice(:,:) 
        real(prec), intent(IN) :: dx 
        real(prec) :: dist(size(H_ice,1),size(H_ice,2))

        ! Local variables 
        integer :: i, j, nx, ny
        integer :: i1, j1  
        real(prec) :: dx_km
        real(prec), allocatable :: dists(:,:) 

        real(prec), parameter :: dist_max = 1e10 

        nx = size(H_ice,1)
        ny = size(H_ice,2) 

        allocate(dists(nx,ny))

        dx_km = dx*1e-3 

        ! Initially set everything to zero distance to margin 
        dist = 0.0 

        do j = 1, ny 
        do i = 1, nx 

            if (H_ice(i,j) .gt. 0.0) then 
                ! Ice-covered point, check min distance to margin 

                dists = dist_max 

                do j1 = 1, ny 
                do i1 = 1, nx 
                    if (H_ice(i1,j1) .eq. 0.0) then 
                        ! Check distance 
                        dists(i1,j1) = sqrt(((i1-i)*dx_km)**2 + ((j1-j)*dx_km)**2) 
                    end if 
                end do 
                end do 

                ! Get minimum distance 
                dist(i,j) = minval(dists)
                
            end if 

        end do 
        end do 


        return 

    end function distance_to_margin
    
    function distance_to_grline(is_grline,f_grnd,dx) result(dist)
        ! ajr: very slow and obsolete routine - do not use! 
        ! instead, use calc_distance_to_grounding_line above
        
        implicit none 
         
        logical,    intent(IN) :: is_grline(:,:)
        real(prec), intent(IN) :: f_grnd(:,:) 
        real(prec), intent(IN) :: dx 
        real(prec) :: dist(size(is_grline,1),size(is_grline,2))

        ! Local variables 
        integer :: i, j, nx, ny
        integer :: i1, j1  
        real(prec) :: dx_km
        real(prec), allocatable :: dists(:,:) 

        real(prec), parameter :: dist_max = 1e10 

        nx = size(is_grline,1)
        ny = size(is_grline,2) 

        allocate(dists(nx,ny))

        dx_km = dx*1e-3 

        ! Initially set everything to zero distance to margin 
        dist = 0.0 

        if (count(is_grline) .gt. 0) then 
            ! If grounding-line points exist, loop over points checking distances

            do j = 1, ny 
            do i = 1, nx 

                if (.not. is_grline(i,j)) then 
                    ! Not at the grounding line, check min distance to grounding line 

                    dists = dist_max 

                    do j1 = 1, ny 
                    do i1 = 1, nx 
                        if (is_grline(i1,j1)) then 
                            ! Check distance 
                            dists(i1,j1) = sqrt(((i1-i)*dx_km)**2 + ((j1-j)*dx_km)**2) 
                        end if 
                    end do 
                    end do 

                    ! Get minimum distance 
                    dist(i,j) = minval(dists)
                    
                    ! For floating points set distance to negative value 
                    if (f_grnd(i,j) .eq. 0.0) dist(i,j) = -dist(i,j) 

                end if 

            end do 
            end do 

        end if 

        return 

    end function distance_to_grline


!! f_grnd calculations from IMAU-ICE / CISM 

! == Routines for determining the grounded fraction on all four grids
  
  subroutine determine_grounded_fractions(f_grnd,f_grnd_acx,f_grnd_acy,f_grnd_ab,H_grnd)
    ! Determine the grounded fraction of centered and staggered grid points
    ! Uses the bilinear interpolation scheme (with analytical solutions) 
    ! from CISM (Leguy et al., 2021), as adapted from IMAU-ICE v2.0 code (rev. 4776833b)
    
    implicit none
    
    real(wp), intent(OUT) :: f_grnd(:,:) 
    real(wp), intent(OUT), optional :: f_grnd_acx(:,:) 
    real(wp), intent(OUT), optional :: f_grnd_acy(:,:) 
    real(wp), intent(OUT), optional :: f_grnd_ab(:,:) 
    real(wp), intent(IN)  :: H_grnd(:,:) 
    
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
    call determine_grounded_fractions_CISM_quads(f_grnd_NW,f_grnd_NE,f_grnd_SW,f_grnd_SE,f_flt)
    
    ! Get grounded fractions on all four grids by averaging over the quadrants
    do j = 1, ny
    do i = 1, nx 
        
      ! Get neighbor indices
      im1 = max(i-1,1) 
      ip1 = min(i+1,nx) 
      jm1 = max(j-1,1) 
      jp1 = min(j+1,ny) 

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
    
    return 

  end subroutine determine_grounded_fractions

  subroutine determine_grounded_fractions_CISM_quads(f_grnd_NW,f_grnd_NE,f_grnd_SW,f_grnd_SE,f_flt)
    ! Calculate grounded fractions of all four quadrants of each a-grid cell
    ! (using the approach from CISM, where grounded fractions are calculated
    !  based on analytical solutions to the bilinear interpolation)
    
    implicit none
    
    real(wp), intent(OUT) :: f_grnd_NW(:,:)
    real(wp), intent(OUT) :: f_grnd_NE(:,:)
    real(wp), intent(OUT) :: f_grnd_SW(:,:)
    real(wp), intent(OUT) :: f_grnd_SE(:,:)
    real(wp), intent(IN)  :: f_flt(:,:) 

    ! Local variables:
    integer  :: i, j, ii, jj, nx, ny
    integer  :: im1, ip1, jm1, jp1  
    real(wp) :: f_NW, f_N, f_NE, f_W, f_m, f_E, f_SW, f_S, f_SE
    real(wp) :: fq_NW, fq_NE, fq_SW, fq_SE
    
    nx = size(f_flt,1)
    ny = size(f_flt,2)

    ! Calculate grounded fractions of all four quadrants of each a-grid cell
    do j = 1, ny
    do i = 1, nx
        
      ! Get neighbor indices
      im1 = max(i-1,1) 
      ip1 = min(i+1,nx) 
      jm1 = max(j-1,1) 
      jp1 = min(j+1,ny) 

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
      if (ABS(dd) < 1e-4_dp) then
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
        write(io_unit_err,*) 'determine_grounded_fractions_CISM_quads - rotate_quad_until_match - ERROR: couldnt find matching scenario!'
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

        real(prec), intent(IN)  :: X
        real(prec) :: ERR
        
        ! Local variables:
        real(prec)              :: EPS
        real(prec)              :: X2
        real(prec)              :: ER
        real(prec)              :: R
        real(prec)              :: C0
        integer                 :: k
        
        EPS = 1.0e-15
        X2  = X * X
        if (abs(X) < 3.5) then
            ER = 1.0
            R  = 1.0
            do k = 1, 50
                R  = R * X2 / (real(k, prec) + 0.5)
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
                R  = -R * (real(k, prec) - 0.5) / X2
                ER = ER + R
                C0  = EXP(-X2) / (abs(X) * sqrt(pi))
                ERR = 1.0 - C0 * ER
                if(X < 0.0) ERR = -ERR
            end do
        end if

        return

    end function error_function

end module topography


