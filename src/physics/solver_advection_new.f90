module solver_advection_new 

    use yelmo_defs, only : sp, dp, prec, tol_underflow
    use yelmo_tools, only : stagger_aa_ab
    
    implicit none 




contains 


    subroutine calc_advec2D(dXdt,X,ux,uy,dx,dy,dt)
        ! Notation following documentation of MITGCM:
        ! https://mitgcm.readthedocs.io/en/latest/algorithm/adv-schemes.html#centered-second-order-advection-diffusion

        implicit none 


        real(prec), intent(OUT)   :: dXdt(:,:) 
        real(prec), intent(INOUT) :: X(:,:) 
        real(prec), intent(IN)    :: ux(:,:) 
        real(prec), intent(IN)    :: uy(:,:) 
        real(prec), intent(IN)    :: dx 
        real(prec), intent(IN)    :: dy 
        real(prec), intent(IN)    :: dt 
        
        ! Local variables 
        real(prec), allocatable   :: Fx(:,:) 
        real(prec), allocatable   :: Fy(:,:) 





        return 

    end subroutine calc_advec2D 

    subroutine calc_flux2D_expl_2nd(Fx,Fy,X,ux,uy,dx,dy)

        implicit none

        real(prec), intent(OUT)   :: Fx(:,:) 
        real(prec), intent(OUT)   :: Fy(:,:) 
        real(prec), intent(INOUT) :: X(:,:) 
        real(prec), intent(IN)    :: ux(:,:) 
        real(prec), intent(IN)    :: uy(:,:) 
        real(prec), intent(IN)    :: dx 
        real(prec), intent(IN)    :: dy

        
        return 

    end subroutine calc_flux2D_expl_2nd

end module solver_advection_new 