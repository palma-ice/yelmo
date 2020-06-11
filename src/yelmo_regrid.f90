module yelmo_regridding 

    use yelmo_defs
    
    implicit none 


contains 

    subroutine calc_regrid_weights(wts,x0,y0,x1,y1)

        implicit none 

        real(prec), intent(OUT) :: wts(:,:) 
        real(prec), intent(IN)  :: x0(:) 
        real(prec), intent(IN)  :: y0(:)
        real(prec), intent(IN)  :: x1(:) 
        real(prec), intent(IN)  :: y1(:)


        return 

    end subroutine calc_regrid_weights 

end module yelmo_regridding
