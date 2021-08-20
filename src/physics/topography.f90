module topography 

    use yelmo_defs, only : wp, prec, sec_year, pi, T0, g, rho_ice, rho_sw, rho_w 

    implicit none 

    private  

    public :: calc_ice_fraction
    public :: calc_z_srf
    public :: calc_z_srf_max
    public :: calc_z_srf_subgrid_area
    public :: calc_H_eff
    public :: calc_H_grnd
    public :: calc_f_grnd_subgrid_area
    public :: calc_f_grnd_subgrid_linear
    public :: filter_f_grnd
    public :: calc_grline
    public :: calc_bmb_total
    public :: calc_fmb_total
    public :: distance_to_grline
    public :: distance_to_margin
    

contains 



    ! ============================================================
    !
    ! Calculations
    !
    ! ============================================================

    subroutine calc_ice_fraction(f_ice,H_ice,f_grnd)
        ! Determine the area fraction of a grid cell
        ! that is ice-covered. Assume that marginal points
        ! have equal thickness to inland neighbors 

        implicit none 

        real(wp), intent(OUT) :: f_ice(:,:)             ! [--] Ice covered fraction (aa-nodes)
        real(wp), intent(IN)  :: H_ice(:,:)             ! [m] Ice thickness on standard grid (aa-nodes)
        real(wp), intent(IN)  :: f_grnd(:,:)            ! [--] Grounded fraction (aa-nodes)

        ! Local variables 
        integer  :: i, j, nx, ny
        integer  :: im1, ip1, jm1, jp1
        real(wp) :: H_ice_x, H_ice_y 
        real(wp) :: H_eff  
        logical  :: get_fractional_cover 

        real(wp) :: H_neighb(4)
        integer  :: n_neighb(4)
        logical  :: mask(4) 
        integer  :: n_now
        integer, allocatable  :: n_ice(:,:) 

        real(wp), parameter :: f_ice_island = 0.1 
        
        nx = size(H_ice,1)
        ny = size(H_ice,2)

        allocate(n_ice(nx,ny))
        n_ice = 0
        
        ! By default, fractional cover will be determined
        get_fractional_cover = .TRUE. 

        ! ajr: see if we can avoid this without relying on "boundaries"
        ! if (trim(boundaries) .eq. "MISMIP3D") then
        !     ! Do not use f_ice treatment for MISMIP3D, it is problematic
        !     ! at the domain boundaries.
        !     get_fractional_cover = .FALSE. 
        ! end if 

        ! For ice-covered points with ice-free neighbors (ie, at the floating or grounded margin),
        ! determine the fraction of grid point that should be ice covered. 

        do j = 1, ny
        do i = 1, nx 

            ! Get neighbor indices
            im1 = max(i-1,1) 
            ip1 = min(i+1,nx) 
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny) 
            
            ! Initially set any point with ice to f_ice=1
            if (H_ice(i,j) .gt. 0.0) then 
                f_ice(i,j) = 1.0 
            else 
                f_ice(i,j) = 0.0 
            end if 

            ! Also count how many neighbors are ice covered 
            if (get_fractional_cover) then 
                H_neighb   = [H_ice(im1,j),H_ice(ip1,j),H_ice(i,jm1),H_ice(i,jp1)]
                mask       = H_neighb .gt. 0.0_wp 
                n_ice(i,j) = count(mask) 
            end if 
        end do 
        end do

        if (get_fractional_cover) then 
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
                    ! Island point, assume the cell is not full to ensure it
                    ! is not dynamically active.

                    f_ice(i,j) = f_ice_island

                else if (H_ice(i,j) .gt. 0.0_wp .and. n_ice(i,j) .lt. 4) then
                    ! This point is ice-covered, but is at the ice margin 

                    ! Get neighbor ice thicknesses and border-ice counts
                    H_neighb = [H_ice(im1,j),H_ice(ip1,j),H_ice(i,jm1),H_ice(i,jp1)]
                    n_neighb = [n_ice(im1,j),n_ice(ip1,j),n_ice(i,jm1),n_ice(i,jp1)]
                    
                    ! Get mask indicating neighbors that both have ice and are surrounded by ice
                    mask     = H_neighb .gt. 0.0_wp .and. n_neighb .eq. 4
                    n_now    = count(mask)

                    if (n_now .ge. 1) then 
                        ! Some neighbors have full ice coverage

                        ! Determine height to give to potentially partially filled cell
                        if (f_grnd(i,j) .eq. 0.0) then 
                            ! Floating point, set H_eff = minimum of neighbors

                            H_eff = minval(H_neighb,mask=mask)

                        else 
                            ! Grounded point, set H_eff < H_min arbitrarily (0.5 works well)

                            ! ajr: following CISM, do not allow partially filled cells
                            ! for grounded ice 
                            H_eff = H_ice(i,j) 

                            ! Alternative approach:
                            ! H_eff = 0.5*minval(H_neighb,mask=mask)

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

                    else 
                        ! No neighbors have full ice coverage, 
                        ! set point to partial coveragae like an island

                        f_ice(i,j) = f_ice_island 

                    end if 

                end if 

            end do 
            end do 

        end if 

        return 

    end subroutine calc_ice_fraction
    
    elemental subroutine calc_z_srf(z_srf,H_ice,H_grnd,z_bed,z_sl)
        ! Calculate surface elevation

        implicit none 

        real(prec), intent(INOUT) :: z_srf
        real(prec), intent(IN)    :: H_ice
        real(prec), intent(IN)    :: H_grnd
        real(prec), intent(IN)    :: z_bed
        real(prec), intent(IN)    :: z_sl

        ! Local variables 
        real(prec) :: rho_ice_sw

        rho_ice_sw = rho_ice/rho_sw ! Ratio of density of ice to seawater [--]
        
        ! Update the surface elevation based on z_bed, H_ice and overburden ice thickness 
        if (H_grnd .gt. 0.0) then 
            ! Grounded ice or ice-free land

            z_srf = z_bed + H_ice 

        else
            ! Floating ice or open ocean

            z_srf = z_sl + (1.0-rho_ice_sw)*H_ice

        end if 
        
        return 

    end subroutine calc_z_srf

    subroutine calc_z_srf_max(z_srf,H_ice,z_bed,z_sl)
        ! Calculate surface elevation
        ! Adapted from Pattyn (2017), Eq. 1
        
        implicit none 

        real(prec), intent(INOUT) :: z_srf(:,:) 
        real(prec), intent(IN)    :: H_ice(:,:)
        real(prec), intent(IN)    :: z_bed(:,:)
        real(prec), intent(IN)    :: z_sl(:,:)

        ! Local variables
        real(prec) :: rho_ice_sw
        real(prec) :: H_mrgn 

        rho_ice_sw = rho_ice/rho_sw ! Ratio of density of ice to seawater [--]
        
        ! Initially calculate surface elevation everywhere 
        z_srf = max(z_bed + H_ice, z_sl + (1.0-rho_ice_sw)*H_ice)
        
        return 

    end subroutine calc_z_srf_max

    subroutine calc_z_srf_subgrid_area(z_srf,f_grnd,H_ice,z_bed,z_sl,gl_sep_nx)
        ! Interpolate variables at grounding line to subgrid level to 
        ! calculate the average z_srf value for the aa-node cell

        implicit none
        
        real(prec), intent(OUT) :: z_srf(:,:)       ! aa-nodes 
        real(prec), intent(IN)  :: f_grnd(:,:)      ! aa-nodes
        real(prec), intent(IN)  :: H_ice(:,:)
        real(prec), intent(IN)  :: z_bed(:,:)
        real(prec), intent(IN)  :: z_sl(:,:)
        integer,    intent(IN)  :: gl_sep_nx        ! Number of interpolation points per side (nx*nx)

        ! Local variables
        integer    :: i, j, nx, ny
        real(prec) :: v1, v2, v3, v4 
        integer    :: im1, ip1, jm1, jp1
        real(prec) :: rho_ice_sw

        real(prec), allocatable :: z_srf_int(:,:) 
        real(prec), allocatable :: H_ice_int(:,:) 
        real(prec), allocatable :: z_bed_int(:,:) 
        real(prec), allocatable :: z_sl_int(:,:) 
        real(prec), allocatable :: H_grnd_int(:,:) 

        rho_ice_sw = rho_ice/rho_sw ! Ratio of density of ice to seawater [--]

        nx = size(z_srf,1)
        ny = size(z_srf,2) 

        ! Allocate the subgrid arrays 
        allocate(z_srf_int(gl_sep_nx,gl_sep_nx))
        allocate(H_ice_int(gl_sep_nx,gl_sep_nx))
        allocate(z_bed_int(gl_sep_nx,gl_sep_nx))
        allocate(z_sl_int(gl_sep_nx,gl_sep_nx))
        allocate(H_grnd_int(gl_sep_nx,gl_sep_nx))
        
        ! Calculate the surface elevation based on whole grid values,
        ! except at the grounding line which is treated with subgrid interpolations. 
        do j = 1, ny 
        do i = 1, nx

            if (f_grnd(i,j) .eq. 1.0) then 
                ! Fully grounded grid point 

                z_srf(i,j) = z_bed(i,j) + H_ice(i,j)

            else if (f_grnd(i,j) .eq. 0.0) then 
                ! Fully floating grid point 

                z_srf(i,j) = z_sl(i,j) + (1.0-rho_ice_sw)*H_ice(i,j)

            else 
                ! Partially grounded and floating (ie, grounding line) 
                ! Perform subgrid calculations 

                ! Get neighbor indices
                im1 = max(i-1,1) 
                ip1 = min(i+1,nx) 
                jm1 = max(j-1,1) 
                jp1 = min(j+1,ny) 
                
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
                
                
                ! Now calculate H_grnd and z_srf for each subgrid point 
                call calc_H_grnd(H_grnd_int,H_ice_int,z_bed_int,z_sl_int)

                where (H_grnd_int .gt. 0.0)
                    ! Fully grounded ice or ice-free land 
                    z_srf_int = z_bed_int + H_ice_int 

                elsewhere 
                    ! Fully floating ice or open ocean
                    z_srf_int = z_sl_int + (1.0-rho_ice_sw)*H_ice_int 

                end where 

                ! Calculate full grid z_srf value as the mean of subgrid values 
                z_srf(i,j) = sum(z_srf_int) / real(gl_sep_nx*gl_sep_nx,prec)

            end if 

        end do 
        end do 

        return
        
    end subroutine calc_z_srf_subgrid_area

    elemental subroutine calc_H_eff(H_corr,H_ice,f_ice)
        ! Calculate ice-thickness, scaled at margins to actual thickness
        ! but as if it covered the whole grid cell.
        
        implicit none

        real(prec), intent(OUT) :: H_corr 
        real(prec), intent(IN)  :: H_ice 
        real(prec), intent(IN)  :: f_ice 

        
        if (f_ice .gt. 0.0) then 
            H_corr = H_ice / f_ice
        else 
            H_corr = H_ice 
        end if

        return

    end subroutine calc_H_eff

    elemental subroutine calc_H_grnd(H_grnd,H_ice,z_bed,z_sl)
        ! Calculate ice thickness overburden, H_grnd
        ! When H_grnd >= 0, grounded, when H_grnd < 0, floating 
        ! Also calculate rate of change for diagnostic related to grounding line 

        ! Note that thickness above flotation H_sl is therefore max(H_grnd,0)
        
        implicit none 

        real(prec), intent(INOUT) :: H_grnd
        real(prec), intent(IN)    :: H_ice
        real(prec), intent(IN)    :: z_bed
        real(prec), intent(IN)    :: z_sl 

        ! Local variables   
        real(prec) :: rho_sw_ice 
        
        rho_sw_ice = rho_sw/rho_ice ! Ratio of density of seawater to ice [--]
        
        ! Calculate new H_grnd (ice thickness overburden)
        H_grnd = H_ice - rho_sw_ice*max(z_sl-z_bed,0.0_prec)

        return 

    end subroutine calc_H_grnd

    subroutine calc_f_grnd_subgrid_area(f_grnd,f_grnd_acx,f_grnd_acy,H_grnd,gl_sep_nx)
        ! Use H_grnd to determined grounded area fraction of grid point.

        implicit none
        
        real(prec), intent(OUT) :: f_grnd(:,:)      ! aa-nodes 
        real(prec), intent(OUT) :: f_grnd_acx(:,:)  ! ac-nodes
        real(prec), intent(OUT) :: f_grnd_acy(:,:)  ! ac-nodes
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

if (.TRUE.) then 
        ! Linear average to ac-nodes 

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

else 
        ! Subgrid area calcs centered on ac-nodes directly

        ! acx-nodes 
        do j = 1, ny 
        do i = 1, nx

            ! Get neighbor indices
            im1 = max(i-1,1) 
            ip1 = min(i+1,nx) 
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny) 
            
            ! === acy nodes === 

            ! First, calculate Hg at corners (acy-nodes)
            Hg_1 = 0.5_prec*(H_grnd(ip1,j) + H_grnd(ip1,jp1))
            Hg_2 = 0.5_prec*(H_grnd(i,j)   + H_grnd(i,jp1))
            Hg_3 = 0.5_prec*(H_grnd(i,j)   + H_grnd(i,jm1))
            Hg_4 = 0.5_prec*(H_grnd(ip1,j) + H_grnd(ip1,jm1))
            
            if (max(Hg_1,Hg_2,Hg_3,Hg_4) .ge. 0.0 .and. min(Hg_1,Hg_2,Hg_3,Hg_4) .lt. 0.0) then 
                ! Point contains grounding line, get grounded area  
                
                call calc_grounded_fraction_cell(f_grnd_acx(i,j),Hg_1,Hg_2,Hg_3,Hg_4,gl_sep_nx)

            else 
                ! Purely grounded or floating point 

                Hg_mid = 0.5_prec*(H_grnd(i,j)  + H_grnd(ip1,j))

                if (Hg_mid .lt. 0.0) then 
                    f_grnd_acx(i,j) = 0.0_prec 
                else 
                    f_grnd_acx(i,j) = 1.0_prec 
                end if 

            end if 

            ! === acy nodes === 

            ! First, calculate Hg at corners (acx-nodes)
            Hg_1 = 0.5_prec*(H_grnd(i,jp1)   + H_grnd(ip1,jp1))
            Hg_2 = 0.5_prec*(H_grnd(im1,jp1) + H_grnd(i,jp1))
            Hg_3 = 0.5_prec*(H_grnd(im1,j)   + H_grnd(i,j))
            Hg_4 = 0.5_prec*(H_grnd(ip1,j)   + H_grnd(i,j))
            
            if (max(Hg_1,Hg_2,Hg_3,Hg_4) .ge. 0.0 .and. min(Hg_1,Hg_2,Hg_3,Hg_4) .lt. 0.0) then 
                ! Point contains grounding line, get grounded area  
                
                call calc_grounded_fraction_cell(f_grnd_acy(i,j),Hg_1,Hg_2,Hg_3,Hg_4,gl_sep_nx)

            else 
                ! Purely grounded or floating point 
                
                Hg_mid = 0.5_prec*(H_grnd(i,j)  + H_grnd(i,jp1))
                
                if (Hg_mid .lt. 0.0) then 
                    f_grnd_acy(i,j) = 0.0_prec 
                else 
                    f_grnd_acy(i,j) = 1.0_prec 
                end if 
                
            end if 

        end do 
        end do

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
    
    subroutine filter_f_grnd(f_grnd)
        ! Remove isolated floating points inside of grounded ice 

        implicit none 

        real(prec), intent(INOUT) :: f_grnd(:,:)
        
        ! Local variables 
        integer :: i, j, nx, ny 
        real(prec), allocatable :: f_grnd_old(:,:) 
        logical, allocatable :: is_grline(:,:) 
        logical, allocatable :: is_grz(:,:) 
        
        nx = size(f_grnd,1)
        ny = size(f_grnd,2)

        allocate(f_grnd_old(nx,ny))
        allocate(is_grline(nx,ny))
        allocate(is_grz(nx,ny))
        
        ! == LAKES ==
        f_grnd_old = f_grnd
        do i = 2, nx-1 
        do j = 2, ny-1 

            if (f_grnd_old(i,j) .lt. 1.0) then 
                ! Check partially floating points 

                if (f_grnd_old(i-1,j)+f_grnd_old(i+1,j)+f_grnd_old(i,j-1)+f_grnd_old(i,j+1) .eq. 4.0) then 
                    ! Point is surrounded by fully grounded ice (ie, a lake),
                    ! so treat it as grounded ice for now 
                    f_grnd(i,j) = 1.0 

                end if 

            end if 

        end do 
        end do 

        ! == GRLINE ISLANDS == 
        f_grnd_old = f_grnd 

        ! Calculate intermediate grounding line estimate 
        call calc_grline(is_grline,is_grz,f_grnd)

        do i = 2, nx-1 
        do j = 2, ny-1 

            if (is_grline(i,j)) then 
                ! Grounding line defined point 

                if ( (f_grnd_old(i-1,j) .eq. 0.0 .or. is_grline(i-1,j)) .and. &
                     (f_grnd_old(i+1,j) .eq. 0.0 .or. is_grline(i+1,j)) .and. &
                     (f_grnd_old(i,j-1) .eq. 0.0 .or. is_grline(i,j-1)) .and. &
                     (f_grnd_old(i,j+1) .eq. 0.0 .or. is_grline(i,j+1)) ) then   
                    ! Point is surrounded by floating or grounding-line ice,
                    ! so treat it as floating ice for now 
                    f_grnd(i,j) = 0.0 

                end if 

            end if 

        end do 
        end do 

        return 

    end subroutine filter_f_grnd

    subroutine calc_grline(is_grline,is_grz,f_grnd)

        implicit none 

        logical,    intent(OUT) :: is_grline(:,:)   ! grounding line 
        logical,    intent(OUT) :: is_grz(:,:)      ! grounding zone 
        real(prec), intent(IN)  :: f_grnd(:,:)

        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: im1, ip1, jm1, jp1 

        nx = size(is_grline,1)
        ny = size(is_grline,2)


        ! 1. Determine grounding line ===========================================

        is_grline = .FALSE. 
 
        do j = 1, ny 
        do i = 1, nx

            im1 = max(1, i-1)
            ip1 = min(nx,i+1)
            
            jm1 = max(1, j-1)
            jp1 = min(ny,j+1)

            ! Grounded point or partially floating point with floating neighbors
            if (f_grnd(i,j) .gt. 0.0 .and. &
                (f_grnd(im1,j) .eq. 0.0 .or. f_grnd(ip1,j) .eq. 0.0 .or. &
                 f_grnd(i,jm1) .eq. 0.0 .or. f_grnd(i,jp1) .eq. 0.0) ) then 
                
                is_grline(i,j) = .TRUE. 

            end if 

        end do 
        end do 

        ! 2. Determine grounding line zone ===========================================
        
        is_grz = .FALSE. 
 
        do j = 1, ny 
        do i = 1, nx

            im1 = max(1, i-1)
            ip1 = min(nx,i+1)
            
            jm1 = max(1, j-1)
            jp1 = min(ny,j+1)

            ! Neighbor with the grounding line, or grounding line 
            is_grz(i,j) = ( count(is_grline(im1:ip1,jm1:jp1)) .gt. 0 )

        end do 
        end do 
        
        return 

    end subroutine calc_grline
    
    subroutine calc_bmb_total(bmb,bmb_grnd,bmb_shlf,f_grnd,diffuse_bmb_shlf)

        implicit none 

        real(prec), intent(OUT) :: bmb(:,:) 
        real(prec), intent(IN)  :: bmb_grnd(:,:) 
        real(prec), intent(IN)  :: bmb_shlf(:,:) 
        real(prec), intent(IN)  :: f_grnd(:,:) 
        logical,    intent(IN)  :: diffuse_bmb_shlf 

        ! Local variables
        integer    :: i, j, nx, ny
        integer    :: n_float  
        real(prec) :: bmb_shlf_now 

        nx = size(bmb,1)
        ny = size(bmb,2) 

        ! Combine floating and grounded parts into one field =========================
        ! (Note: use of where-statement allows bmb_shlf to be set negative everywhere,
        ! but only be applied where it is relevant.)
        where (f_grnd .eq. 1.0) 
            ! Purely grounded, use only bmb_grnd 
            bmb = bmb_grnd 
        elsewhere
            ! Comine bmb_grnd and bmb_shlf 

            ! Add the two fields for now to avoid complications with f_grnd
            bmb = bmb_grnd + bmb_shlf 

            ! Weighted average
            !bmb = f_grnd*bmb_grnd + (1.0-f_grnd)*bmb_shlf

        end where 


        if (diffuse_bmb_shlf) then 
            ! Allow marine melt (bmb_shlf) to permeate inland at the grounding line,
            ! to induce more effective retreat in warm periods 

            do j = 1, ny
            do i = 1, nx

                if (f_grnd(i,j) .eq. 1.0) then 
                    ! Grounded point, look for floating neighbors 

                    if (.FALSE.) then
                        ! 9-neighbours method

                        n_float = count(f_grnd(i-1:i+1,j-1:j+1) .lt. 1.0)

                        if (n_float .gt. 0) then 
                            ! bmb_shelf is the mean of the neighbours
                            bmb_shlf_now = sum(bmb_shlf(i-1:i+1,j-1:j+1),mask=f_grnd(i-1:i+1,j-1:j+1) .lt. 1.0) / real(n_float,prec)
                            bmb(i,j)     = (1.0-n_float/9.0)*bmb_grnd(i,j) + (n_float/9.0)*bmb_shlf_now

                        end if

                    else 
                        ! 5-neighbors method 

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
        real(wp) :: H_ref, dz  
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

                        ! Get margin-scaled ice thickness (f_ice > 0 since H_ice > 0, except on first timestep, be safe!)
                        if (f_ice(i,j) .eq. 0.0_wp) then 
                            ! f_ice not yet defined, assume point covers cell 
                            H_ref = H_ice(i,j) 
                        else 
                            H_ref = H_ice(i,j) / f_ice(i,j) 
                        end if 
                        
                        ! Determine depth of adjacent water using centered cell information
                        if (H_grnd(i,j) .lt. 0.0_wp) then 
                            ! Cell is floating, calculate submerged ice thickness 
                            
                            dz = (H_ref*rho_ice_sw)
                            
                        else 
                            ! Cell is grounded, recover depth of seawater

                            dz = max( (H_ref - H_grnd(i,j)) * rho_ice_sw, 0.0_wp)

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

        implicit none 

        real(prec), intent(OUT) :: f_g 
        real(prec), intent(IN)  :: Hg_1,Hg_2,Hg_3,Hg_4
        integer,    intent(IN)  :: nx                    ! Number of interpolation points 

        ! Local variables 
        integer :: i, j 
        real(prec) :: x(nx), y(nx) 
        real(prec) :: Hg_int(nx,nx)  

        ! Populate x,y axes for interpolation points (between 0 and 1)
        do i = 1, nx 
            x(i) = 0.0 + real(i-1)/real(nx-1)
        end do 
        y = x 


        ! Calculate interpolation heights
        Hg_int = 0.0 
        do i = 1, nx 
        do j = 1, nx 

            Hg_int(i,j) = interp_bilin_pt(Hg_1,Hg_2,Hg_3,Hg_4,x(i),y(j))

        end do 
        end do 

        ! Calculate weighted fraction (assume all points have equal weight)
        f_g = real(count(Hg_int .ge. 0.0),prec) / real(nx*nx,prec)

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

end module topography 


