module climate_adjustments
    
    use yelmo_defs, only : sp, dp, prec, pi, missing_value, mv, tol_underflow, rho_ice, g 
    
    implicit none 

    private 


contains 

    subroutine downscale_smb_gradient(smb,z_srf,smb_in,z_srf_in,xx,yy,radius,method)

        implicit none 

        real(prec), intent(OUT) :: smb(:,:)
        real(prec), intent(IN)  :: z_srf(:,:) 
        real(prec), intent(IN)  :: smb_in(:,:) 
        real(prec), intent(IN)  :: z_srf_in(:,:) 
        real(prec), intent(IN)  :: xx(:,:) 
        real(prec), intent(IN)  :: yy(:,:) 
        real(prec), intent(IN)  :: radius 
        character(len=*), intent(IN) :: method 

        ! Local variables 
        integer :: i, j, nx, ny, npts 
        integer :: i0, i1, j0, j1 
        real(prec) :: x0, x1, y0, y1 

        real(prec) :: xpts(1000), ypts(1000) 

        nx = size(smb,1)
        ny = size(smb,2) 

        smb = MV 

        ! Loop over each target smb point, determine neighbors,
        ! build relationship of smb_in vs z_srf_in in neighborhood (separate for positive and negative smb),
        ! Determine best-guess smb for this point given current z_srf.
        ! Note: instead of a circle, just use a square to avoid calculating distances.
        do j = 1, ny 
        do i = 1, nx 

            ! Step 1: Find neighbors to current point, and populate input data points xpts and ypts 

            select case(trim(method))

                case("square") 

                    x0 = xx(i,j) - radius 
                    x1 = xx(i,j) + radius 
                    y0 = yy(i,j) - radius 
                    y1 = yy(i,j) + radius 

                    i0 = minloc(abs(xx(:,j)-x0),dim=1)
                    i1 = minloc(abs(xx(:,j)-x1),dim=1)
                    j0 = minloc(abs(yy(i,:)-y0),dim=1)
                    j1 = minloc(abs(yy(i,:)-y1),dim=1)

                    npts = (i1-i0+1)*(j1-j0+1)

                    ! Store all neighbors in box 
                    xpts = mv 
                    ypts = mv 
                    xpts(1:npts) = reshape(z_srf_in(i0:i1,j0:j1),[npts])
                    ypts(1:npts) = reshape(smb_in(i0:i1,j0:j1),[npts])

                case("circle")
                    ! To do 

                    write(*,*) "downscale_smb_gradient:: Error: method not yet implemented: "//trim(method)

                case("basin")
                    ! To do 

                    write(*,*) "downscale_smb_gradient:: Error: method not yet implemented: "//trim(method)

                case DEFAULT 
                    write(*,*) "downscale_smb_gradient:: Error: method not recognized: "//trim(method)
                    stop 

            end select 

            ! Step 2: Use neighbors to calculate smb vs z_srf relationship 


            
        end do 
        end do 


        return 

    end subroutine downscale_smb_gradient

    subroutine calc_distances(dist,x,y,xx,yy,mask)

        implicit none 

        real(prec), intent(OUT) :: dist(:,:)  
        real(prec), intent(IN)  :: x 
        real(prec), intent(IN)  :: y 
        real(prec), intent(IN)  :: xx(:,:) 
        real(prec), intent(IN)  :: yy(:,:) 
        logical,    intent(IN)  :: mask(:,:) 

        ! Local variables 
        integer :: i, j, nx, ny 

        nx = size(dist,1)
        ny = size(dist,2) 


        ! Populate distance matrix where possible 
        dist = MV 

        do j = 1, ny 
        do i = 1, nx 
            if (mask(i,j)) then 
                dist(i,j) = sqrt( (xx(i,j)-x)**2 + (yy(i,j)-y)**2 ) 
            end if 
        end do 
        end do 

        return 

    end subroutine calc_distances

end module climate_adjustments 
