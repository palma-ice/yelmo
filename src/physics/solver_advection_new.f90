module solver_advection_new 

    use yelmo_defs, only : sp, dp, wp, tol_underflow

    implicit none 




contains 


    subroutine calc_advec2D(dXdt,X,ux,uy,dx,dy,dt)
        ! Notation following documentation of MITGCM:
        ! https://mitgcm.readthedocs.io/en/latest/algorithm/adv-schemes.html#centered-second-order-advection-diffusion

        implicit none 


        real(wp), intent(OUT)   :: dXdt(:,:) 
        real(wp), intent(IN)    :: X(:,:) 
        real(wp), intent(IN)    :: ux(:,:) 
        real(wp), intent(IN)    :: uy(:,:) 
        real(wp), intent(IN)    :: dx 
        real(wp), intent(IN)    :: dy 
        real(wp), intent(IN)    :: dt 
        
        ! Local variables 
        real(wp), allocatable   :: Fx(:,:) 
        real(wp), allocatable   :: Fy(:,:) 





        return 

    end subroutine calc_advec2D

    subroutine calc_flux2D_expl_2nd(Fx,Fy,X,ux,uy,dx,dy)

        implicit none

        real(wp), intent(OUT)   :: Fx(:,:) 
        real(wp), intent(OUT)   :: Fy(:,:) 
        real(wp), intent(IN   ) :: X(:,:) 
        real(wp), intent(IN)    :: ux(:,:) 
        real(wp), intent(IN)    :: uy(:,:) 
        real(wp), intent(IN)    :: dx 
        real(wp), intent(IN)    :: dy

        
        return 

    end subroutine calc_flux2D_expl_2nd

end module solver_advection_new 