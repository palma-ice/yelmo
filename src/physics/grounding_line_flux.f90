module grounding_line_flux
    ! Routines to prescribe grounding-line fluxes (eg, Schoof or Tsai laws)

    use yelmo_defs, only : sp, dp, prec, rho_ice, rho_sw, g  

    implicit none 


    private 

    public :: calc_grounding_line_flux 
    


contains 

    subroutine calc_grounding_line_flux(qq_gl_acx,qq_gl_acy,H_ice,ATT_bar,f_grnd,f_grnd_acx,f_grnd_acy,n_glen,Q0,f_drag)

        implicit none 

        real(prec), intent(OUT) :: qq_gl_acx(:,:)       ! [m^2 a^-1] Diagnosed grounding line flux (acx-nodes)
        real(prec), intent(OUT) :: qq_gl_acy(:,:)       ! [m^2 a^-1] Diagnosed grounding line flux (acy-nodes)
        real(prec), intent(IN)  :: H_ice(:,:)           ! [m] Ice thickness (aa-nodes)
        real(prec), intent(IN)  :: ATT_bar(:,:)         ! [a^-1 Pa^-3] Depth-averaged rate factor
        real(prec), intent(OUT) :: f_grnd(:,:)          ! [--] Grounding line fraction (aa-nodes) - binary variable
        real(prec), intent(OUT) :: f_grnd_acx(:,:)      ! [--] Grounding line fraction (aa-nodes) - binary variable
        real(prec), intent(OUT) :: f_grnd_acy(:,:)      ! [--] Grounding line fraction (aa-nodes) - binary variable
        real(prec), intent(IN)  :: n_glen               ! [--] Glen's flow law exponent 
        real(prec), intent(IN)  :: Q0                   ! [?] Scaling coefficient, in the range of 0.60-0.65
        real(prec), intent(IN)  :: f_drag               ! [--] Dragging coefficient, f_drag ~ 0.6 

        ! Local variables
        integer    :: i, j, nx, ny
        logical    :: is_gl  
        real(prec) :: H_gl, A_gl, qq_gl    

        nx = size(H_ice,1)
        ny = size(H_ice,2) 

        ! Initially set gl-flux to zero everywhere 
        qq_gl_acx = 0.0
        qq_gl_acy = 0.0 

        ! acx-nodes 
        do j = 1, ny 
        do i = 1, nx
            
            ! Determine if grounding-line sits between aa-nodes
            is_gl = (f_grnd(i,j) .gt. 0.0 .and. f_grnd(i+1,j) .eq. 0.0) .or. & 
                    (f_grnd(i,j) .eq. 0.0 .and. f_grnd(i+1,j) .gt. 0.0)

            if (is_gl) then 
                ! For grounding-line points, diagnose flux 

                ! Determine ice thicknes at the grounding line 
                if (f_grnd(i,j) .gt. 0.0 .and. f_grnd(i+1,j) .eq. 0.0) then 
                    ! Floating to the right

                    H_gl = H_ice(i,j)*f_grnd_acx(i,j)   + (1.0-f_grnd_acx(i,j))*H_ice(i+1,j)
                    A_gl = ATT_bar(i,j)*f_grnd_acx(i,j) + (1.0-f_grnd_acx(i,j))*ATT_bar(i+1,j)
                    
                    qq_gl = calc_gl_flux_tsai(H_gl,A_gl,n_glen,Q0,f_drag)

                else if (f_grnd(i,j) .eq. 0.0 .and. f_grnd(i+1,j) .gt. 0.0) then
                    ! Floating to the left

                    H_gl = H_ice(i,j)*(1.0-f_grnd_acx(i,j))   + f_grnd_acx(i,j)*H_ice(i+1,j)
                    A_gl = ATT_bar(i,j)*(1.0-f_grnd_acx(i,j)) + f_grnd_acx(i,j)*ATT_bar(i+1,j)

                    qq_gl = calc_gl_flux_tsai(H_gl,A_gl,n_glen,Q0,f_drag)
                    
                end if 

                
            end if 

        end do 
        end do 

        ! acy-nodes 
        do j = 1, ny 
        do i = 1, nx

        end do 
        end do 



        return 

    end subroutine calc_grounding_line_flux

    elemental function calc_gl_flux_tsai(H_gl,A_gl,n_glen,Q0,f_drag) result(qq_gl)

        implicit none 

        real(prec), intent(IN)  :: H_gl                 ! [m] Grounding-line ice thickness
        real(prec), intent(IN)  :: A_gl                 ! [a-1 Pa-3] Grounding-line rate factor
        real(prec), intent(IN)  :: n_glen               ! [--] Glen's flow law exponent 
        real(prec), intent(IN)  :: Q0                   ! [?] Scaling coefficient, in the range of 0.60-0.65
        real(prec), intent(IN)  :: f_drag               ! [--] Dragging coefficient, f_drag ~ 0.6 
        real(prec) :: qq_gl                             ! [m2 / a] Grounding-line flux 
        
        ! Local variables 
        real(prec) :: density_factor

        ! Calculate constant density factor 
        density_factor = (1.0 - rho_ice/rho_sw)**(n_glen-1.0)

        ! Calculate grounding line flux following Tsai et al. (2015), Eq. 38
        qq_gl = Q0 * 8.0 * A_gl * (rho_ice*g)**n_glen / ((4.0**n_glen)*f_drag) & 
                                * density_factor * H_gl**(n_glen+2.0)

        return 

    end function calc_gl_flux_tsai



end module grounding_line_flux 
