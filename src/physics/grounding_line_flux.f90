module grounding_line_flux
    ! Routines to prescribe grounding-line fluxes (eg, Schoof or Tsai laws)

    use yelmo_defs, only : sp, dp, prec, rho_ice, rho_sw, g  

    implicit none 


    private 

    public :: calc_grounding_line_flux 
    


contains 

    subroutine calc_grounding_line_flux(qq_gl_acx,qq_gl_acy,H_ice,ATT_bar,C_bed,ux_bar,uy_bar, &
                                            f_grnd,f_grnd_acx,f_grnd_acy,n_glen,m_drag,Q0,f_drag,glf_method)

        implicit none 

        real(prec), intent(OUT) :: qq_gl_acx(:,:)       ! [m^2 a^-1] Diagnosed grounding line flux (acx-nodes)
        real(prec), intent(OUT) :: qq_gl_acy(:,:)       ! [m^2 a^-1] Diagnosed grounding line flux (acy-nodes)
        real(prec), intent(IN)  :: H_ice(:,:)           ! [m] Ice thickness (aa-nodes)
        real(prec), intent(IN)  :: ATT_bar(:,:)         ! [a^-1 Pa^-n] Depth-averaged rate factor
        real(prec), intent(IN)  :: C_bed(:,:)           ! [Pa a m^-1] Basal friction coefficient 
        real(prec), intent(IN)  :: ux_bar(:,:)          ! [m a^-1] Depth-averaged horizontal veloctiy, x-direction (acx-nodes)
        real(prec), intent(IN)  :: uy_bar(:,:)          ! [m a^-1] Depth-averaged horizontal veloctiy, y-direction (acy-nodes)
        real(prec), intent(IN)  :: f_grnd(:,:)          ! [--] Grounding line fraction (aa-nodes) - binary variable
        real(prec), intent(IN)  :: f_grnd_acx(:,:)      ! [--] Grounding line fraction (acx-nodes)
        real(prec), intent(IN)  :: f_grnd_acy(:,:)      ! [--] Grounding line fraction (acy-nodes)
        real(prec), intent(IN)  :: n_glen               ! [--] Glen's flow law exponent 
        real(prec), intent(IN)  :: m_drag               ! [--] Power law exponent
        real(prec), intent(IN)  :: Q0                   ! [?] Scaling coefficient, in the range of 0.60-0.65
        real(prec), intent(IN)  :: f_drag               ! [--] Dragging coefficient, f_drag ~ 0.6 
        character(len=*), intent(IN) :: glf_method      ! "power" or "coulomb" method 

        ! Local variables
        integer    :: i, j, nx, ny
        logical    :: is_gl, float_right, float_left   
        real(prec) :: H_gl, A_gl, C_bed_gl, qq_gl    
        real(prec) :: qq_left, qq_right, H_now, f_lin

        nx = size(H_ice,1)
        ny = size(H_ice,2) 

        ! Initially set gl-flux to zero everywhere 
        qq_gl_acx = 0.0
        qq_gl_acy = 0.0 

        ! acx-nodes 
        do j = 1, ny 
        do i = 3, nx-3
            
            ! Determine if grounding-line sits between aa-nodes
            float_right = (f_grnd(i,j) .gt. 0.0 .and. f_grnd(i+1,j) .eq. 0.0)
            float_left  = (f_grnd(i,j) .eq. 0.0 .and. f_grnd(i+1,j) .gt. 0.0)
            is_gl       = (float_right .or. float_left)
                    

            if (is_gl) then 
                ! For grounding-line points, diagnose flux 

                ! 1. Determine quantities at the grounding line ==============================

                if (float_right) then 
                    ! Floating to the right

                    H_gl = H_ice(i,j)  *f_grnd_acx(i,j) + (1.0-f_grnd_acx(i,j))*H_ice(i+1,j)
                    A_gl = ATT_bar(i,j)*f_grnd_acx(i,j) + (1.0-f_grnd_acx(i,j))*ATT_bar(i+1,j)
                    
                    C_bed_gl = C_bed(i,j)   ! Upstream (non-zero) C_bed 

                else
                    ! Floating to the left

                    H_gl = H_ice(i,j)  *(1.0-f_grnd_acx(i,j)) + f_grnd_acx(i,j)*H_ice(i+1,j)
                    A_gl = ATT_bar(i,j)*(1.0-f_grnd_acx(i,j)) + f_grnd_acx(i,j)*ATT_bar(i+1,j)

                    C_bed_gl = C_bed(i+1,j)   ! Upstream (non-zero) C_bed 
                    
                end if 

                ! 2. Calculate flux directly at the grounding line ============================

                select case(trim(glf_method))

                    case("power")
                        ! Calculate magnitude of flux given a Coulomb friction relation, following Schoof (2007)
                        
                        qq_gl = calc_gl_flux_power(H_gl,A_gl,C_bed_gl,n_glen,m_drag)

                    case("coulomb")
                        ! Calculate magnitude of flux given a Coulomb friction relation, following Tsai et al. (2015)
                        
                        qq_gl = calc_gl_flux_coulomb(H_gl,A_gl,n_glen,Q0,f_drag)

                    case DEFAULT 

                        write(*,*) "calc_grounding_line_flux:: Error: glf_method not recognized, &
                                   &must be 'coulomb' (Tsai) or 'power' (Schoof)"
                        write(*,*) "glf_method: ", trim(glf_method)
                        stop 

                end select 

                ! Assign a direction to the flux as well (direction of flow)
                qq_gl = sign(qq_gl,ux_bar(i,j))

                ! 3. Interpolate prescribed gl-flux to ac-nodes to the left and right of gl. ==
                
                if (f_grnd_acx(i,j) .lt. 0.5) then 
                    ! Grounding line is between acx(i,j) and acx(i-1,j), so 
                    ! it should be interpolated to those two points

                    ! Left side 

                    H_now   = 0.5*(H_ice(i-2,j) + H_ice(i-1,j)) 
                    qq_left = H_now * ux_bar(i-2,j) 

                    ! Get the ratio of left grid cell width to the 
                    ! total distance between qq_left and qq_gl 
                    f_lin = 1.0 / ( 1.0 + (0.5 + f_grnd_acx(i,j)) )

                    qq_gl_acx(i-1,j) = (1.0-f_lin)*qq_left + f_lin*qq_gl

                    ! Right side 

                    H_now   = 0.5*(H_ice(i+1,j) + H_ice(i+2,j)) 
                    qq_right = H_now * ux_bar(i+1,j) 

                    ! Get the ratio of distance from gl to ac(i,j) to the 
                    ! total distance between qq_gl and qq_right 
                    f_lin = (0.5-f_grnd_acx(i,j)) / ( 1.0 + (0.5-f_grnd_acx(i,j)) )

                    qq_gl_acx(i,j) = (1.0-f_lin)*qq_gl + f_lin*qq_right

                else
                    ! Grounding line is between acx(i,j) and acx(i+1,j)

                    ! Left side 

                    H_now   = 0.5*(H_ice(i-1,j) + H_ice(i,j)) 
                    qq_left = H_now * ux_bar(i-1,j) 

                    ! Get the ratio of left grid cell width to the 
                    ! total distance between qq_left and qq_gl 
                    f_lin = 1.0 / ( 1.0 + (f_grnd_acx(i,j)-0.5) )

                    qq_gl_acx(i,j) = (1.0-f_lin)*qq_left + f_lin*qq_gl

                    ! Right side 

                    H_now   = 0.5*(H_ice(i+2,j) + H_ice(i+3,j)) 
                    qq_right = H_now * ux_bar(i+2,j) 

                    ! Get the ratio of distance from gl to ac(i,j) to the 
                    ! total distance between qq_gl and qq_right  
                    f_lin = (1.5-f_grnd_acx(i,j)) / ( 1.0 + (1.5-f_grnd_acx(i,j)) )

                    qq_gl_acx(i+1,j) = (1.0-f_lin)*qq_gl + f_lin*qq_right

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

    elemental function calc_gl_flux_power(H_gl,A_gl,C_bed,n_glen,m_drag) result(qq_gl)
        ! Calculate the analytical grounding-line flux solution
        ! following Schoof (2007) 
        
        implicit none 

        real(prec), intent(IN)  :: H_gl                 ! [m] Grounding-line ice thickness
        real(prec), intent(IN)  :: A_gl                 ! [a-1 Pa-3] Grounding-line rate factor
        real(prec), intent(IN)  :: C_bed                ! [Pa a/m] Basal friction coefficient
        real(prec), intent(IN)  :: n_glen               ! [--] Glen's flow law exponent 
        real(prec), intent(IN)  :: m_drag               ! [--] Power-law dragging exponent
        real(prec) :: qq_gl                             ! [m2 / a] Grounding-line flux 
        
        ! Calculate grounding line flux following Schoof (2007), Eq. 20
        ! Note: m_drag = 1/3 in the original Schoof formula
        qq_gl = ( A_gl * (rho_ice*g)**(n_glen+1.0) * (1.0 - rho_ice/rho_sw)**n_glen / ((4.0**n_glen)*C_bed) ) **(1.0/(m_drag+1.0)) & 
                                 * H_gl**( (m_drag+n_glen+3.0)/(m_drag+1.0) )

        return 

    end function calc_gl_flux_power

    elemental function calc_gl_flux_coulomb(H_gl,A_gl,n_glen,Q0,f_drag) result(qq_gl)
        ! Calculate the analytical grounding-line flux solution
        ! following Tsai et al. (2015)

        implicit none 

        real(prec), intent(IN)  :: H_gl                 ! [m] Grounding-line ice thickness
        real(prec), intent(IN)  :: A_gl                 ! [a-1 Pa-3] Grounding-line rate factor
        real(prec), intent(IN)  :: n_glen               ! [--] Glen's flow law exponent 
        real(prec), intent(IN)  :: Q0                   ! [?] Scaling coefficient, in the range of 0.60-0.65
        real(prec), intent(IN)  :: f_drag               ! [--] Dragging coefficient, f_drag ~ 0.6 
        real(prec) :: qq_gl                             ! [m2 / a] Grounding-line flux 
        
        ! Calculate grounding line flux following Tsai et al. (2015), Eq. 38
        qq_gl = Q0 * 8.0 * A_gl * (rho_ice*g)**n_glen / ((4.0**n_glen)*f_drag) & 
                                * (1.0 - rho_ice/rho_sw)**(n_glen-1.0) * H_gl**(n_glen+2.0)

        return 

    end function calc_gl_flux_coulomb


end module grounding_line_flux 
