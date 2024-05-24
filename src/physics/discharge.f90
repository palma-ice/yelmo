module discharge

    use yelmo_defs,  only : sp, dp, wp, TOL_UNDERFLOW, degrees_to_radians

    implicit none

    private
    public :: calc_mb_discharge

contains

    subroutine calc_mb_discharge(mb_discharge,H_ice,z_srf,z_bed_sd,dist_grline, &
                    dist_margin,f_ice,dx,alpha_max,tau_mbd,sigma_ref,m_d,m_r)
        ! Calculate implicit subgrid calving discharge rate
        ! following Calov et al. (2015)

        implicit none

        real(wp), intent(OUT) :: mb_discharge(:,:)
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: z_srf(:,:)
        real(wp), intent(IN)  :: z_bed_sd(:,:)
        real(wp), intent(IN)  :: dist_grline(:,:)       ! [km]
        real(wp), intent(IN)  :: dist_margin(:,:)       ! [km]
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: dx                     ! [m]
        real(wp), intent(IN)  :: alpha_max
        real(wp), intent(IN)  :: tau_mbd
        real(wp), intent(IN)  :: sigma_ref
        real(wp), intent(IN)  :: m_d
        real(wp), intent(IN)  :: m_r

        ! Local variables
        integer :: i, j, nx, ny
        real(wp) :: cosalpha_max
        real(wp) :: cosalpha
        real(wp) :: dx_km
        real(wp) :: f_sd
        real(wp) :: f_l
        real(wp) :: f_r  

        real(wp), parameter :: dl = 32.0            ! [km] Length scale 
        real(wp), parameter :: dist_max = 500.0     ! [km] Maximum distance from coast to calculate discharge
    
        dx_km = dx*1e-3

        nx = size(mb_discharge,1)
        ny = size(mb_discharge,2)

        ! Get cosine of alpha_max, alpha=60deg --> cosalpha_max=0.5
        cosalpha_max = cos(alpha_max*degrees_to_radians)

        ! Loop over domain, calculate discharge at each relevant point
        do j = 1, ny 
        do i = 1, nx 

            ! First assume discharge is zero
            mb_discharge(i,j) = 0.0

            if (H_ice(i,j) .gt. 0.0 .and. dist_grline(i,j) .ge. 0.0 .and. &
                                                dist_margin(i,j) .le. dist_max) then 
            !if (H_ice(i,j) .gt. 0.0 .and. dist_grline(i,j) .ge. 0.0) then
                ! Ice exists, is grounded and lies within dist_max distance
                ! to the ice margin.

                ! Calculate the angle of descent relative to the direction of the coast
                call calc_coastal_cosine_angle(cosalpha,z_srf,dist_grline,dx,i,j)

                if (cosalpha .ge. cosalpha_max) then
                    ! Ice is flowing towards the coast, more or less, apply parameterization

                    ! Calculate scaling factors (roughness, distance to coast and resolution)
                    f_sd = tanh(z_bed_sd(i,j)/sigma_ref)
                    f_l  = ( dl / (dl + dist_grline(i,j)) )**m_d
                    f_r  = ( dx_km/dl )**m_r 

                    ! Calculate subgrid discharge mass balance rate (negative, for mass loss)
                    mb_discharge(i,j) = -( f_sd * f_l * f_r * (H_ice(i,j) / tau_mbd) )

                end if

            end if

        end do
        end do

        return

    end subroutine calc_mb_discharge

    subroutine calc_coastal_cosine_angle(cosalpha,z_srf,dist,dx,i,j)

        implicit none

        real(wp), intent(OUT) :: cosalpha
        real(wp), intent(IN)  :: z_srf(:,:)
        real(wp), intent(IN)  :: dist(:,:)
        real(wp), intent(IN)  :: dx
        integer,  intent(IN)  :: i 
        integer,  intent(IN)  :: j 

        ! Local variables
        integer :: nx, ny 
        real(wp) :: twodx
        real(wp) :: dzsdx 
        real(wp) :: dzsdy
        real(wp) :: dzsdxy
        real(wp) :: dldx 
        real(wp) :: dldy
        real(wp) :: dldxy

        nx = size(z_srf,1)
        ny = size(z_srf,2) 

        twodx = 2.0 * dx 

        if (i .eq. 1) then 
            dzsdx = (z_srf(i+1,j) - z_srf(i,j)) / dx 
            dldx  = (dist(i+1,j)  - dist(i,j))  / dx
        else if (i .eq. nx) then 
            dzsdx = (z_srf(i,j) - z_srf(i-1,j)) / dx 
            dldx  = (dist(i,j)  - dist(i-1,j))  / dx
        else
            dzsdx = (z_srf(i+1,j) - z_srf(i-1,j)) / twodx 
            dldx  = (dist(i+1,j)  - dist(i-1,j))  / twodx
        end if 

        if (j .eq. 1) then 
            dzsdy = (z_srf(i,j+1) - z_srf(i,j)) / dx 
            dldy  = (dist(i,j+1)  - dist(i,j))  / dx
        else if (j .eq. ny) then 
            dzsdy = (z_srf(i,j) - z_srf(i,j-1)) / dx 
            dldy  = (dist(i,j)  - dist(i,j-1))  / dx
        else
            dzsdy = (z_srf(i,j+1) - z_srf(i,j-1)) / twodx 
            dldy  = (dist(i,j+1)  - dist(i,j-1))  / twodx
        end if 

        ! Get magnitudes
        dzsdxy = sqrt( dzsdx**2 + dzsdy**2)
        dldxy  = sqrt( dldx**2 + dldy**2)
        
        if (dzsdxy .gt. 0.0 .and. dldxy .gt. 0.0) then 
            cosalpha = (dzsdx*dldx + dzsdy*dldy) / (dzsdxy*dldxy)
        else
            cosalpha = 0.0
        end if
        
        return

    end subroutine calc_coastal_cosine_angle

end module discharge