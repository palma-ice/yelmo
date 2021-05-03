module yelmo_regridding 

    use yelmo_defs
    
    implicit none 

    private 
    public :: yelmo_regrid 

contains 
    
    subroutine yelmo_regrid(dom1,dom0,dx1)

        implicit none 

        type(yelmo_class), intent(INOUT) :: dom1        ! Target Yelmo object 
        type(yelmo_class), intent(IN)    :: dom0        ! Initial Yelmo object 
        real(prec),        intent(IN)    :: dx1         ! [km] Target grid resolution

        ! Local variables 
        integer :: i, j, k 
        real(prec) :: f_grid      ! Ratio between initial and target grid
        real(prec), allocatable :: f_grid_allowed(:) 

        ! First, make sure that the f_grid ratio is appropriate
        ! (target grid must be a multiple of the source grid) 
        f_grid_allowed = [0.0625_prec,0.125_prec,0.25_prec,0.5_prec,1.0_prec,2.0_prec,4.0_prec,8.0_prec,16.0_prec]
        
        f_grid = dom0%grd%dx*1e-3 / dx1 

        if (.not. any(f_grid .eq. f_grid_allowed,1)) then 

            write(*,*) "yelmo_regrid:: Error: target resolution must be a multiple of the original resolution."
            write(*,*) "f_grid = ", f_grid 
            stop 

        end if 


        
        return 

    end subroutine yelmo_regrid

    
    subroutine calc_regrid_weights(wts,x0,y0,x1,y1)

        implicit none 

        real(prec), intent(OUT) :: wts(:,:) 
        real(prec), intent(IN)  :: x0(:) 
        real(prec), intent(IN)  :: y0(:)
        real(prec), intent(IN)  :: x1(:) 
        real(prec), intent(IN)  :: y1(:)


        return 

    end subroutine calc_regrid_weights

    subroutine regrid_downsample(dst,src,dx_dst,dx_src)

        implicit none 

        real(prec), intent(OUT) :: dst(:,:) 
        real(prec), intent(IN)  :: src(:,:) 
        real(prec), intent(IN)  :: dx_dst 
        real(prec), intent(IN)  :: dx_src 

        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: i1, j1, n 
        integer :: scale 

        nx = size(dst,1)
        ny = size(dst,2)

        ! Scale should be a positive integer 
        scale = int(dx_src/dx_dst) 

        ! Loop over target grid 
        do j = 1, ny 
        do i = 1, nx 

            ! Determine point on source grid that coincides 
            ! with this point 
            !i1 = 


        end do 
        end do  

        return 

    end subroutine regrid_downsample

end module yelmo_regridding
