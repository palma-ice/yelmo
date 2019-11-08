module yelmo_tools
    ! Generic functions and subroutines that could be used in many contexts:
    ! math, vectors, sorting, etc. 

    use yelmo_defs, only : sp, dp, prec, missing_value, tol_underflow, pi
    
    implicit none 

    private 
    public :: calc_magnitude 
    public :: calc_magnitude_from_staggered
    public :: calc_magnitude_from_staggered_ice
    public :: stagger_ac_aa
    public :: stagger_aa_ab
    public :: stagger_aa_ab_ice 
    public :: stagger_ab_aa 
    public :: stagger_aa_acx
    public :: stagger_aa_acy
    public :: stagger_acx_aa
    public :: stagger_acy_aa
    public :: stagger_ab_acx
    public :: stagger_ab_acy 
    public :: calc_gradient_ac
    public :: calc_gradient_ac_ice
    public :: calc_gradient_ac_gl

    public :: mean_mask
    public :: minmax
    public :: fill_borders_2D
    public :: fill_borders_3D 
    
    public :: smooth_gauss_2D
    public :: smooth_gauss_3D
    public :: gauss_values

    public :: regularize2D 

    ! Integration functions
    public :: test_integration
    public :: integrate_trapezoid1D_pt
    public :: integrate_trapezoid1D_1D
    public :: calc_vertical_integrated_2D
    public :: calc_vertical_integrated_3D
    
contains 

    elemental function calc_magnitude(u,v) result(umag)
        ! Get the vector magnitude from two components at the same grid location

        implicit none 

        real(prec), intent(IN)  :: u, v 
        real(prec) :: umag 

        umag = sqrt(u*u+v*v)

        return

    end function calc_magnitude 
    
    function calc_magnitude_from_staggered(u,v) result(umag)
        ! Calculate the centered (aa-nodes) magnitude of a vector 
        ! from the staggered (ac-nodes) components

        implicit none 
        
        real(prec), intent(IN)  :: u(:,:), v(:,:)  
        real(prec) :: umag(size(u,1),size(u,2)) 

        ! Local variables 
        integer :: i, j, nx, ny 
        real(prec) :: unow, vnow 

        nx = size(u,1)
        ny = size(u,2) 

        umag = 0.0_prec 

        do j = 2, ny-1 
        do i = 2, nx-1 
            unow = 0.5_prec*(u(i,j)+u(i-1,j))
            vnow = 0.5_prec*(v(i,j)+v(i,j-1))
            umag(i,j) = sqrt(unow*unow+vnow*vnow)
        end do 
        end do 

        return

    end function calc_magnitude_from_staggered 
    
    function calc_magnitude_from_staggered_ice(u,v,H) result(umag)
        ! Calculate the centered (aa-nodes) magnitude of a vector 
        ! from the staggered (ac-nodes) components

        implicit none 
        
        real(prec), intent(IN)  :: u(:,:), v(:,:), H(:,:) 
        real(prec) :: umag(size(u,1),size(u,2)) 

        ! Local variables 
        integer :: i, j, nx, ny 
        real(prec) :: unow, vnow 
        real(prec) :: f1, f2, H1, H2 
        nx = size(u,1)
        ny = size(u,2) 

        umag = 0.0_prec 

        do j = 2, ny-1 
        do i = 2, nx-1 

            H1 = 0.5_prec*(H(i-1,j)+H(i,j))
            H2 = 0.5_prec*(H(i,j)+H(i+1,j))

            f1 = 0.5_prec 
            f2 = 0.5_prec 
            if (H1 .eq. 0.0) f1 = 0.0_prec  
            if (H2 .eq. 0.0) f2 = 0.0_prec   

            if (f1+f2 .gt. 0.0) then 
                unow = (f1*u(i-1,j) + f2*u(i,j)) / (f1+f2)
                if (abs(unow) .lt. tol_underflow) unow = 0.0_prec 
            else 
                unow = 0.0 
            end if 

            H1 = 0.5_prec*(H(i,j-1)+H(i,j))
            H2 = 0.5_prec*(H(i,j)+H(i,j+1))

            f1 = 0.5_prec 
            f2 = 0.5_prec 
            if (H1 .eq. 0.0) f1 = 0.0_prec  
            if (H2 .eq. 0.0) f2 = 0.0_prec   

            if (f1+f2 .gt. 0.0) then 
                vnow = (f1*v(i,j-1) + f2*v(i,j)) / (f1+f2)
                if (abs(vnow) .lt. tol_underflow) vnow = 0.0_prec 
            else 
                vnow = 0.0 
            end if 

            umag(i,j) = sqrt(unow*unow+vnow*vnow)
        end do 
        end do 

        return

    end function calc_magnitude_from_staggered_ice 
    
    function stagger_ac_aa(u,v) result(umag)
        ! Calculate the centered (aa-node) magnitude of a scalar 
        ! from the staggered (ac-node) components

        implicit none 
        
        real(prec), intent(IN)  :: u(:,:), v(:,:)    ! acx-, acy-nodes 
        real(prec) :: umag(size(u,1),size(u,2))      ! aa-nodes 

        ! Local variables 
        integer :: i, j, nx, ny 

        nx = size(u,1)
        ny = size(u,2) 

        umag = 0.0_prec 

        do j = 2, ny-1 
        do i = 2, nx-1 
            umag(i,j) = 0.25_prec*(u(i,j)+u(i-1,j)+v(i,j)+v(i,j-1))
        end do 
        end do 

        return

    end function stagger_ac_aa 
    
    function stagger_aa_ab(u) result(ustag)
        ! Stagger from Aa => Ab
        ! Four point average from corner Aa nodes to central Ab node 

        implicit none 

        real(prec), intent(IN)  :: u(:,:) 
        real(prec) :: ustag(size(u,1),size(u,2)) 

        ! Local variables 
        integer :: i, j, nx, ny  

        nx = size(u,1)
        ny = size(u,2) 

        ustag = 0.0_prec 

        do j = 1, ny-1 
        do i = 1, nx-1
            ustag(i,j) = 0.25_prec*(u(i+1,j+1)+u(i+1,j)+u(i,j+1)+u(i,j))
        end do 
        end do 

        return

    end function stagger_aa_ab 
    
    function stagger_aa_ab_ice(u,H_ice) result(ustag)
        ! Stagger from Aa => Ab
        ! Four point average from corner Aa nodes to central Ab node 

        implicit none 

        real(prec), intent(IN)  :: u(:,:) 
        real(prec), intent(IN)  :: H_ice(:,:) 
        real(prec) :: ustag(size(u,1),size(u,2)) 

        ! Local variables 
        integer :: i, j, nx, ny, k   

        nx = size(u,1)
        ny = size(u,2) 

        ustag = 0.0_prec 

        do j = 1, ny-1 
        do i = 1, nx-1
            k = 0 
            ustag(i,j) = 0.0 
            if (H_ice(i,j) .gt. 0.0) then 
                ustag(i,j) = ustag(i,j) + u(i,j) 
                k = k+1
            end if 

            if (H_ice(i+1,j) .gt. 0.0) then 
                ustag(i,j) = ustag(i,j) + u(i+1,j) 
                k = k+1 
            end if 
            
            if (H_ice(i,j+1) .gt. 0.0) then 
                ustag(i,j) = ustag(i,j) + u(i,j+1) 
                k = k+1 
            end if 
            
            if (H_ice(i+1,j+1) .gt. 0.0) then 
                ustag(i,j) = ustag(i,j) + u(i+1,j+1) 
                k = k+1 
            end if 
            
            if (k .gt. 0) then 
                ustag(i,j) = ustag(i,j) / real(k,prec)
            end if 

            !ustag(i,j) = 0.25_prec*(u(i+1,j+1)+u(i+1,j)+u(i,j+1)+u(i,j))
        end do 
        end do 

        return

    end function stagger_aa_ab_ice 
    
    function stagger_ab_aa(u) result(ustag)
        ! Stagger from Ab => Aa
        ! Four point average from corner Ab nodes to central Aa node 

        implicit none 

        real(prec), intent(IN)  :: u(:,:) 
        real(prec) :: ustag(size(u,1),size(u,2)) 

        ! Local variables 
        integer :: i, j, nx, ny  

        nx = size(u,1)
        ny = size(u,2) 

        ustag = 0.0_prec 

        do j = 2, ny 
        do i = 2, nx
            ustag(i,j) = 0.25_prec*(u(i,j)+u(i-1,j)+u(i,j-1)+u(i-1,j-1))
        end do 
        end do 

        return

    end function stagger_ab_aa 
    
    function stagger_aa_acx(u) result(ustag)
        ! Stagger from Aa => Ac, x-direction 

        implicit none

        real(prec), intent(IN)  :: u(:,:) 
        real(prec) :: ustag(size(u,1),size(u,2)) 

        ! Local variables 
        integer :: i, j, nx, ny  

        nx = size(u,1)
        ny = size(u,2) 

        ustag = 0.0_prec 

        do j = 1, ny 
        do i = 1, nx-1
            ustag(i,j) = 0.5_prec*(u(i,j)+u(i+1,j))
        end do 
        end do 

        return

    end function stagger_aa_acx 
    
    function stagger_aa_acy(u) result(ustag)
        ! Stagger from Aa => Ac 

        implicit none 

        real(prec), intent(IN)  :: u(:,:) 
        real(prec) :: ustag(size(u,1),size(u,2)) 

        ! Local variables 
        integer :: i, j, nx, ny  

        nx = size(u,1)
        ny = size(u,2) 

        ustag = 0.0_prec 

        do j = 1, ny-1 
        do i = 1, nx
            ustag(i,j) = 0.5_prec*(u(i,j)+u(i,j+1))
        end do 
        end do 

        return

    end function stagger_aa_acy 
    
    function stagger_acx_aa(u) result(ustag)
        ! Stagger from Aa => Ac, x-direction 

        implicit none

        real(prec), intent(IN)  :: u(:,:) 
        real(prec) :: ustag(size(u,1),size(u,2)) 

        ! Local variables 
        integer :: i, j, nx, ny  

        nx = size(u,1)
        ny = size(u,2) 

        ustag = 0.0_prec 

        do j = 1, ny 
        do i = 2, nx
            ustag(i,j) = 0.5_prec*(u(i-1,j)+u(i,j))
        end do 
        end do 

        return

    end function stagger_acx_aa 
    
    function stagger_acy_aa(u) result(ustag)
        ! Stagger from Aa => Ac 

        implicit none 

        real(prec), intent(IN)  :: u(:,:) 
        real(prec) :: ustag(size(u,1),size(u,2)) 

        ! Local variables 
        integer :: i, j, nx, ny  

        nx = size(u,1)
        ny = size(u,2) 

        ustag = 0.0_prec 

        do j = 2, ny
        do i = 1, nx
            ustag(i,j) = 0.5_prec*(u(i,j-1)+u(i,j))
        end do 
        end do 

        return

    end function stagger_acy_aa 
    
    function stagger_ab_acx(u) result(ustag)
        ! Stagger from Ab => Ac, x-direction 

        implicit none

        real(prec), intent(IN)  :: u(:,:) 
        real(prec) :: ustag(size(u,1),size(u,2)) 

        ! Local variables 
        integer :: i, j, nx, ny  

        nx = size(u,1)
        ny = size(u,2) 

        ustag = 0.0_prec 

        do j = 2, ny 
        do i = 1, nx
            ustag(i,j) = 0.5_prec*(u(i,j)+u(i,j-1))
        end do 
        end do 

        return

    end function stagger_ab_acx 
    
    function stagger_ab_acy(u) result(ustag)
        ! Stagger from Ab => Ac 

        implicit none 

        real(prec), intent(IN)  :: u(:,:) 
        real(prec) :: ustag(size(u,1),size(u,2)) 

        ! Local variables 
        integer :: i, j, nx, ny  

        nx = size(u,1)
        ny = size(u,2) 

        ustag = 0.0_prec 

        do j = 1, ny 
        do i = 2, nx
            ustag(i,j) = 0.5_prec*(u(i,j)+u(i-1,j))
        end do 
        end do 

        return

    end function stagger_ab_acy 
    
    subroutine calc_gradient_ac(dvardx,dvardy,var,dx)
        ! Calculate gradient on ac nodes 

        implicit none 

        real(prec), intent(OUT) :: dvardx(:,:) 
        real(prec), intent(OUT) :: dvardy(:,:) 
        real(prec), intent(IN)  :: var(:,:) 
        real(prec), intent(IN)  :: dx 

        ! Local variables 
        integer :: i, j, nx, ny 
        real(prec) :: dy 

        nx = size(var,1)
        ny = size(var,2)

        ! Assume y-resolution is identical to x-resolution 
        dy = dx 

        ! Slope in x-direction
        do j = 1, ny 
        do i = 1, nx-1 
            dvardx(i,j) = (var(i+1,j)-var(i,j))/dx 
        end do 
        end do 

        ! Slope in y-direction
        do j = 1, ny-1 
        do i = 1, nx 
            dvardy(i,j) = (var(i,j+1)-var(i,j))/dy
        end do 
        end do 


        return 

    end subroutine calc_gradient_ac
    
    subroutine calc_gradient_ac_ice(dvardx,dvardy,var,H_ice,dx,margin2nd,grad_lim)
        ! Calculate gradient on ac nodes 
        ! for an ice sheet, only using non-zero thickness points

        implicit none 

        real(prec), intent(OUT) :: dvardx(:,:) 
        real(prec), intent(OUT) :: dvardy(:,:) 
        real(prec), intent(IN)  :: var(:,:) 
        real(prec), intent(IN)  :: H_ice(:,:)
        real(prec), intent(IN)  :: dx 
        logical,    intent(IN)  :: margin2nd 
        real(prec), intent(IN)  :: grad_lim 

        ! Local variables 
        integer :: i, j, nx, ny 
        real(prec) :: dy 
        real(prec) :: H0, H1, H2 

        nx = size(var,1)
        ny = size(var,2)

        ! Assume y-resolution is identical to x-resolution 
        dy = dx 

        ! Slope in x-direction
        do j = 1, ny 
        do i = 1, nx-1 
            dvardx(i,j) = (var(i+1,j)-var(i,j))/dx 
        end do 
        end do 

        ! Slope in y-direction
        do j = 1, ny-1 
        do i = 1, nx 
            dvardy(i,j) = (var(i,j+1)-var(i,j))/dy
        end do 
        end do 

        ! === Modify margin gradients =========================
        ! Following Saito et al (2007) by applying a second-order, upwind gradient

        if (margin2nd) then

            ! Slope in x-direction
            do j = 1, ny 
            do i = 3, nx-3

                if (H_ice(i,j) .gt. 0.0 .and. H_ice(i+1,j) .eq. 0.0) then 
                    ! Margin point (ice-free to the right)

                    H0 = var(i+1,j) 
                    H1 = var(i,j)
                    H2 = var(i-1,j)

                    if (H_ice(i-1,j) .gt. 0.0 .and. H_ice(i-2,j) .gt. 0.0) then 
                        ! Second-order upwind if possible 
                        dvardx(i,j) = (3.0*H0 - 4.0*H1 + H2) / (2.0*dx)
                    end if 

                else if (H_ice(i+1,j) .gt. 0.0 .and. H_ice(i,j) .eq. 0.0) then
                    ! Margin point (ice-free to the left)

                    H0 = var(i,j) 
                    H1 = var(i+1,j)
                    H2 = var(i+2,j)

                    if (H_ice(i+2,j) .gt. 0.0 .and. H_ice(i+3,j) .gt. 0.0) then 
                        ! Second-order upwind if possible
                        dvardx(i,j) = -(3.0*H0 - 4.0*H1 + H2) / (2.0*dx) 
                    end if 


                end if 

            end do 
            end do 
            
            ! Slope in y-direction
            do j = 3, ny-3 
            do i = 1, nx

                if (H_ice(i,j) .gt. 0.0 .and. H_ice(i,j+1) .eq. 0.0) then 
                    ! Margin point (ice-free to the top)

                    H0 = var(i,j+1) 
                    H1 = var(i,j)
                    H2 = var(i,j-1)

                    if (H_ice(i,j-1) .gt. 0.0 .and. H_ice(i,j-2) .gt. 0.0) then 
                        ! Second-order upwind if possible
                        dvardy(i,j) = (3.0*H0 - 4.0*H1 + H2) / (2.0*dy)
                    end if 

                else if (H_ice(i,j+1) .gt. 0.0 .and. H_ice(i,j) .eq. 0.0) then
                    ! Margin point (ice-free to the bottom)

                    H0 = var(i,j) 
                    H1 = var(i,j+1)
                    H2 = var(i,j+2)

                    if (H_ice(i,j+2) .gt. 0.0 .and. H_ice(i,j+3) .gt. 0.0) then 
                        ! Second-order upwind if possible
                        dvardy(i,j) = -(3.0*H0 - 4.0*H1 + H2) / (2.0*dy)
                    end if 
                    
                end if 

            end do 
            end do

        end if 

        ! Finally, ensure that gradient is beneath desired limit 
        call minmax(dvardx,grad_lim)
        call minmax(dvardy,grad_lim)

        return 

    end subroutine calc_gradient_ac_ice
    
    subroutine calc_gradient_ac_gl(dvardx,dvardy,var,H_ice, &
                                      f_grnd_acx,f_grnd_acy,dx,method,grad_lim)

        implicit none 

        real(prec), intent(OUT) :: dvardx(:,:)
        real(prec), intent(OUT) :: dvardy(:,:) 
        real(prec), intent(IN)  :: var(:,:) 
        real(prec), intent(IN)  :: H_ice(:,:)
        real(prec), intent(IN)  :: f_grnd_acx(:,:)
        real(prec), intent(IN)  :: f_grnd_acy(:,:)
        real(prec), intent(IN)  :: dx 
        integer,    intent(IN)  :: method           ! Which gl gradient calculation to use
        real(prec), intent(IN)  :: grad_lim         ! Very high limit == 0.05, low limit < 0.01 

        ! Local variables 
        integer :: i, j, nx, ny 
        real(prec) :: dy
        real(prec) :: dvardx_1, dvardx_2 

        nx = size(H_ice,1)
        ny = size(H_ice,2) 

        dy = dx 

        select case(method)

            case(0)  
                ! Do nothing, use the standard no-subgrid treatment 

            case(1)
                ! Weighted average using the grounded fraction (ac-nodes)
                ! or one-sided choice
                ! between surface slope and virtual slope of 
                ! floating ice (using ice thickness)

                ! x-direction 
                do j = 1, ny 
                do i = 1, nx-1 

                    if ( f_grnd_acx(i,j) .gt. 0.0 .and. f_grnd_acx(i,j) .lt. 1.0) then 
                        ! Grounding line point (ac-node)

                        ! Get slope of grounded point and virtual floating point (using H_ice),
                        ! then assume slope is the weighted average of the two 
                        dvardx_1    = (var(i+1,j)-var(i,j)) / dx 
                        dvardx_2    = 0.0 !(H_ice(i+1,j)-H_ice(i,j)) / dx 
                        dvardx(i,j) = f_grnd_acx(i,j)*dvardx_1 + (1.0-f_grnd_acx(i,j))*dvardx_2  
                        
                        ! Limit the slope 
                        call minmax(dvardx(i,j),grad_lim)  
                                   
                    end if 

                end do 
                end do 

                ! y-direction 
                do j = 1, ny-1 
                do i = 1, nx 

                    if ( f_grnd_acy(i,j) .gt. 0.0 .and. f_grnd_acy(i,j) .lt. 1.0) then 
                        ! Grounding line point (ac-node)

                        ! Get slope of grounded point and virtual floating point (using H_ice),
                        ! then assume slope is the weighted average of the two 
                        dvardx_1    = (var(i,j+1)-var(i,j)) / dx 
                        dvardx_2    = 0.0 !(H_ice(i,j+1)-H_ice(i,j)) / dx 
                        dvardy(i,j) = f_grnd_acy(i,j)*dvardx_1 + (1.0-f_grnd_acy(i,j))*dvardx_2  
                        
                        ! Limit the slope 
                        call minmax(dvardy(i,j),grad_lim)  
                         
                    end if 

                end do 
                end do 

            case(2)
                ! One-sided differences upstream and downstream of the grounding line
                ! analgous to Feldmann et al. (2014, JG)

                ! x-direction 
                do j = 1, ny 
                do i = 1, nx-1 

                    if ( f_grnd_acx(i,j) .gt. 0.0 .and. f_grnd_acx(i,j) .lt. 1.0) then 
                        ! Grounding line point (ac-node)

                        if (f_grnd_acx(i,j) .gt. 0.5) then 
                            ! Consider grounded 
                            dvardx(i,j) = (var(i+1,j)-var(i,j)) / dx 
                        else 
                            ! Consider floating 
                            !dvardx(i,j) = (H_ice(i+1,j)-H_ice(i,j)) / dx
                            dvardx(i,j) = 0.0 
                        end if 

                        ! Limit the slope 
                        call minmax(dvardx(i,j),grad_lim)  

                    end if 

                end do 
                end do 

                ! y-direction 
                do j = 1, ny-1 
                do i = 1, nx 

                    if ( f_grnd_acy(i,j) .gt. 0.0 .and. f_grnd_acy(i,j) .lt. 1.0) then 
                        ! Grounding line point (ac-node)

                        if (f_grnd_acy(i,j) .gt. 0.5) then 
                            ! Consider grounded 
                            dvardy(i,j) = (var(i,j+1)-var(i,j)) / dy 
                        else 
                            ! Consider floating 
!                             dvardy(i,j) = (H_ice(i,j+1)-H_ice(i,j)) / dy
                            dvardy(i,j) = 0.0 
                        end if 
                        
                        ! Limit the slope 
                        call minmax(dvardy(i,j),grad_lim)  

                    end if 

                end do 
                end do 

            case DEFAULT  
                
                write(*,*) "calc_gradient_ac_gl:: Error: grad_gl_method not recognized."
                write(*,*) "grad_gl_method = ", method 
                stop 

        end select

        return 

    end subroutine calc_gradient_ac_gl

    function mean_mask(var,mask) result(ave)

        implicit none 

        real(prec), intent(IN) :: var(:,:) 
        logical,    intent(IN) :: mask(:,:) 
        real(prec) :: ave 
        integer :: n 

        n = count(mask)
        
        if (n .gt. 0) then 
            ave = sum(var,mask=mask) / real(n,prec)
        else 
            ave = 0.0 
        end if 

        return 

    end function mean_mask 
    
    elemental subroutine minmax(var,var_lim)

        implicit none 

        real(prec), intent(INOUT) :: var 
        real(prec), intent(IN)    :: var_lim 

        if (var .lt. -var_lim) then 
            var = -var_lim 
        else if (var .gt. var_lim) then 
            var =  var_lim 
        end if 

        return 

    end subroutine minmax 

    subroutine fill_borders_2D(var,nfill)

        implicit none 

        real(prec), intent(INOUT) :: var(:,:) 
        integer,    intent(IN)    :: nfill        ! How many neighbors to fill in 

        ! Local variables 
        integer :: i, j, nx, ny, q 
        
        nx = size(var,1)
        ny = size(var,2)

        do q = 1, nfill 
            var(q,:)      = var(nfill+1,:)      
            var(nx-q+1,:) = var(nx-nfill,:)   
            
            var(:,q)      = var(:,nfill+1)     
            var(:,ny-q+1) = var(:,ny-nfill)  
        end do 

        return 

    end subroutine fill_borders_2D

    subroutine fill_borders_3D(var,nfill)
        ! 3rd dimension is not filled (should be vertical dimension)

        implicit none 

        real(prec), intent(INOUT) :: var(:,:,:) 
        integer,    intent(IN)    :: nfill        ! How many neighbors to fill in 

        ! Local variables 
        integer :: i, j, nx, ny, q 
        
        nx = size(var,1)
        ny = size(var,2)

        do q = 1, nfill 
            var(q,:,:)      = var(nfill+1,:,:)      
            var(nx-q+1,:,:) = var(nx-nfill,:,:)   
            
            var(:,q,:)      = var(:,nfill+1,:)     
            var(:,ny-q+1,:) = var(:,ny-nfill,:)  
        end do 

        return 

    end subroutine fill_borders_3D

    subroutine smooth_gauss_3D(var,mask_apply,dx,n_smooth,mask_use)

        ! Smooth out strain heating to avoid noise 

        implicit none

        real(prec), intent(INOUT) :: var(:,:,:)      ! nx,ny,nz_aa: 3D variable
        logical,    intent(IN)    :: mask_apply(:,:) 
        real(prec), intent(IN)    :: dx 
        integer,    intent(IN)    :: n_smooth  
        logical,    intent(IN), optional :: mask_use(:,:) 

        ! Local variables
        integer    :: k, nz_aa  

        nz_aa = size(var,3)

        do k = 1, nz_aa 
             call smooth_gauss_2D(var(:,:,k),mask_apply,dx,n_smooth,mask_use)
        end do 

        return 

    end subroutine smooth_gauss_3D
    
    subroutine smooth_gauss_2D(var,mask_apply,dx,n_smooth,mask_use)
        ! Smooth out a field to avoid noise 
        ! mask_apply designates where smoothing should be applied 
        ! mask_use   designates which points can be considered in the smoothing filter 

        implicit none

        real(prec), intent(INOUT) :: var(:,:)      ! [nx,ny] 2D variable
        logical,    intent(IN)    :: mask_apply(:,:) 
        real(prec), intent(IN)    :: dx 
        integer,    intent(IN)    :: n_smooth  
        logical,    intent(IN), optional :: mask_use(:,:) 

        ! Local variables
        integer :: i, j, nx, ny, n, n2
        real(prec) :: sigma    
        real(prec), allocatable :: filter0(:,:), filter(:,:) 
        real(prec), allocatable :: var_old(:,:) 
        logical,    allocatable :: mask_use_local(:,:) 

        nx    = size(var,1)
        ny    = size(var,2)
        n     = 5 
        n2    = (n-1)/2 

        sigma = dx*n_smooth 

        allocate(var_old(nx,ny))
        allocate(mask_use_local(nx,ny))
        allocate(filter0(n,n))
        allocate(filter(n,n))

        ! Check whether mask_use is available 
        if (present(mask_use)) then 
            ! use mask_use to define neighborhood points
            
            mask_use_local = mask_use 

        else
            ! Assume that mask_apply also gives the points to use for smoothing 

            mask_use_local = mask_apply
        
        end if

        ! Calculate default 2D Gaussian smoothing kernel
        filter0 = gauss_values(dx,dx,sigma=sigma,n=n)

        var_old = var 

        do j = n2+1, ny-n2
        do i = n2+1, nx-n2

            if (mask_apply(i,j)) then 
                ! Apply smoothing to this point 

                filter = filter0 
                where(.not. mask_use_local(i-n2:i+n2,j-n2:j+n2)) filter = 0.0
                if (sum(filter) .gt. 0.0) then
                    ! If neighbors are available, perform smoothing   
                    filter = filter/sum(filter)
                    var(i,j) = sum(var_old(i-n2:i+n2,j-n2:j+n2)*filter) 
                end if  

            end if 

        end do 
        end do 

        return 

    end subroutine smooth_gauss_2D

    function gauss_values(dx,dy,sigma,n) result(filt)
        ! Calculate 2D Gaussian smoothing kernel
        ! https://en.wikipedia.org/wiki/Gaussian_blur

        implicit none 

        real(prec), intent(IN) :: dx 
        real(prec), intent(IN) :: dy 
        real(prec), intent(IN) :: sigma 
        integer,    intent(IN) :: n 
        real(prec) :: filt(n,n) 

        ! Local variables 
        real(prec) :: x, y  
        integer    :: n2, i, j, i1, j1  

        if (mod(n,2) .ne. 1) then 
            write(*,*) "gauss_values:: error: n can only be odd."
            write(*,*) "n = ", n 
        end if 

        n2 = (n-1)/2 

        do j = -n2, n2 
        do i = -n2, n2 
            x = i*dx 
            y = j*dy 

            i1 = i+1+n2 
            j1 = j+1+n2 
            filt(i1,j1) = 1.0/(2.0*pi*sigma**2)*exp(-(x**2+y**2)/(2*sigma**2))

        end do 
        end do 
        
        ! Normalize to ensure sum to 1
        filt = filt / sum(filt)

        return 

    end function gauss_values

    ! ================================================================================
    !
    ! Regularizing/smoothing functions 
    !
    ! ================================================================================

    subroutine regularize2D(var,H_ice)
        ! Ensure smoothness in 2D fields (ie, no checkerboard patterns)

        implicit none 

        real(prec), intent(INOUT) :: var(:,:)      ! aa-nodes
        real(prec), intent(IN)    :: H_ice(:,:)     ! aa-nodes
        
        ! Local variables
        integer    :: i, j, nx, ny, nlow,  nhi, n   
        integer    :: im1, ip1, jm1, jp1 
        real(prec), allocatable :: var0(:,:) 
        real(prec) :: varx(2), vary(2), var9(3,3)
        logical    :: check_x, check_y 
        
        nx = size(var,1)
        ny = size(var,2) 

        allocate(var0(nx,ny))
        var0 = var 

        do j = 2, ny-1 
        do i = 2, nx-1

            if (H_ice(i,j) .gt. 0.0) then 
                ! Only apply to ice-covered points 

                im1 = max(1, i-1)
                ip1 = min(nx,i+1)
                
                jm1 = max(1, j-1)
                jp1 = min(ny,j+1)

                varx = [var0(im1,j),var0(ip1,j)]
                where([H_ice(im1,j),H_ice(ip1,j)] .eq. 0.0_prec) varx = missing_value 

                vary = [var0(i,jm1),var0(i,jp1)]
                where([H_ice(i,jm1),H_ice(i,jp1)] .eq. 0.0_prec) vary = missing_value 
                
                ! Check if checkerboard exists in each direction 
                check_x = (count(varx .gt. var0(i,j) .and. varx.ne.missing_value) .eq. 2 .or. &
                           count(varx .lt. var0(i,j) .and. varx.ne.missing_value) .eq. 2) 

                check_y = (count(vary .gt. var0(i,j) .and. vary.ne.missing_value) .eq. 2 .or. &
                           count(vary .lt. var0(i,j) .and. vary.ne.missing_value) .eq. 2) 
                
                if (check_x .or. check_y) then 
                    ! Checkerboard exists, apply 9-point neighborhood average

                    var9 = var0(i-1:i+1,j-1:j+1)
                    where(H_ice(i-1:i+1,j-1:j+1) .eq. 0.0_prec) var9 = missing_value 

                    n = count(var9 .ne. missing_value) 

                    var(i,j) = sum(var9,mask=var9.ne.missing_value) / real(n,prec)

                end if 

            end if 

        end do 
        end do 

        return 

    end subroutine regularize2D 

    ! === Generic integration functions ============

    function calc_vertical_integrated_3D(var,zeta) result(var_int)
        ! Vertically integrate a field 3D field (nx,ny,nz)
        ! layer by layer (in the z-direction), return a 3D array

        implicit none

        real(prec), intent(IN) :: var(:,:,:)
        real(prec), intent(IN) :: zeta(:)
        real(prec) :: var_int(size(var,1),size(var,2),size(var,3))

        ! Local variables 
        integer :: i, j, nx, ny

        nx = size(var,1)
        ny = size(var,2)

        do j = 1, ny
        do i = 1, nx
            var_int(i,j,:) = integrate_trapezoid1D_1D(var(i,j,:),zeta)
        end do
        end do

        return

    end function calc_vertical_integrated_3D

    function calc_vertical_integrated_2D(var,zeta) result(var_int)
        ! Vertically integrate a field 3D field (nx,ny,nz) 
        ! to the surface, return a 2D array (nx,ny)
        
        implicit none

        real(prec), intent(IN) :: var(:,:,:)
        real(prec), intent(IN) :: zeta(:)
        real(prec) :: var_int(size(var,1),size(var,2))

        ! Local variables 
        integer :: i, j, nx, ny

        nx = size(var,1)
        ny = size(var,2)

        do j = 1, ny
        do i = 1, nx
            var_int(i,j) = integrate_trapezoid1D_pt(var(i,j,:),zeta)
        end do
        end do

        return

    end function calc_vertical_integrated_2D
    
    function integrate_trapezoid1D_pt(var,zeta) result(var_int)
        ! Integrate a variable from the base to height zeta(nk) in the ice column.
        ! The value of the integral using the trapezium rule can be found using
        ! integral = (b - a)*((f(a) +f(b))/2 + Σ_1_n-1(f(k)) )/n 
        ! Returns a point of integrated value of var at level zeta(nk).

        implicit none

        real(prec), intent(IN) :: var(:)
        real(prec), intent(IN) :: zeta(:)
        real(prec) :: var_int

        ! Local variables 
        integer :: k, nk

        nk = size(var,1)

        ! Initial value is zero
        var_int = 0.0_prec 

        ! Intermediate values include sum of all previous values 
        ! Take current value as average between points
        do k = 2, nk
             var_int = var_int + 0.5_prec*(var(k)+var(k-1))*(zeta(k) - zeta(k-1))
        end do

        return

    end function integrate_trapezoid1D_pt

    function integrate_trapezoid1D_1D(var,zeta) result(var_int)
        ! Integrate a variable from the base to each layer zeta of the ice column.
        ! Note this is designed assuming indices 1 = base, nk = surface 
        ! The value of the integral using the trapezium rule can be found using
        ! integral = (b - a)*((f(a) +f(b))/2 + Σ_1_n-1(f(k)) )/n 
        ! Returns a 1D array with integrated value at each level 

        implicit none

        real(prec), intent(IN) :: var(:)
        real(prec), intent(IN) :: zeta(:)
        real(prec) :: var_int(size(var,1))

        ! Local variables 
        integer :: k, nk

        nk = size(var,1)

        ! Initial value is zero
        var_int(1:nk) = 0.0_prec 

        ! Intermediate values include sum of all previous values 
        ! Take current value as average between points
        do k = 2, nk
             var_int(k:nk) = var_int(k:nk) + 0.5_prec*(var(k)+var(k-1))*(zeta(k) - zeta(k-1))
        end do
        
        return

    end function integrate_trapezoid1D_1D

    subroutine simpne(x,y,result)
        !*****************************************************************************80
        !
        !! SIMPNE approximates the integral of unevenly spaced data.
        !
        !  Discussion:
        !
        !    The routine repeatedly interpolates a 3-point Lagrangian polynomial 
        !    to the data and integrates that exactly.
        !
        !  Modified:
        !
        !    10 February 2006
        !
        !  Reference:
        !
        !    Philip Davis, Philip Rabinowitz,
        !    Methods of Numerical Integration,
        !    Second Edition,
        !    Dover, 2007,
        !    ISBN: 0486453391,
        !    LC: QA299.3.D28.
        !
        !  Parameters:
        !
        !    Input, integer ( kind = 4 ) NTAB, number of data points.  
        !    NTAB must be at least 3.
        !
        !    Input, real ( kind = 8 ) X(NTAB), contains the X values of the data,
        !    in order.
        !
        !    Input, real ( kind = 8 ) Y(NTAB), contains the Y values of the data.
        !
        !    Output, real ( kind = 8 ) RESULT.
        !    RESULT is the approximate value of the integral.
        
        implicit none

        real(prec) :: x(:)
        real(prec) :: y(:)
        real(prec) :: result

        integer :: ntab

        real(prec) :: del(3)
        real(prec) :: e
        real(prec) :: f
        real(prec) :: feints
        real(prec) :: g(3)
        integer    :: i
        integer    :: n
        real(prec) :: pi(3)
        real(prec) :: sum1

        real(prec) :: x1
        real(prec) :: x2
        real(prec) :: x3

        ntab = size(x,1) 

        result = 0.0D+00

        if ( ntab <= 2 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'SIMPNE - Fatal error!'
            write ( *, '(a)' ) '  NTAB <= 2.'
            stop 1
        end if
     
        n = 1
     
        do
     
            x1 = x(n)
            x2 = x(n+1)
            x3 = x(n+2)
            e = x3 * x3- x1 * x1
            f = x3 * x3 * x3 - x1 * x1 * x1
            feints = x3 - x1

            del(1) = x3 - x2
            del(2) = x1 - x3
            del(3) = x2 - x1

            g(1) = x2 + x3
            g(2) = x1 + x3
            g(3) = x1 + x2

            pi(1) = x2 * x3
            pi(2) = x1 * x3
            pi(3) = x1 * x2
     
            sum1 = 0.0D+00
            do i = 1, 3
                sum1 = sum1 + y(n-1+i) * del(i) &
                    * ( f / 3.0D+00 - g(i) * 0.5D+00 * e + pi(i) * feints )
            end do
            result = result - sum1 / ( del(1) * del(2) * del(3) )
     
            n = n + 2

            if ( ntab <= n + 1 ) then
            exit
            end if

        end do
     
        if ( mod ( ntab, 2 ) /= 0 ) then
            return
        end if

        n = ntab - 2
        x3 = x(ntab)
        x2 = x(ntab-1)
        x1 = x(ntab-2)
        e = x3 * x3 - x2 * x2
        f = x3 * x3 * x3 - x2 * x2 * x2
        feints = x3 - x2

        del(1) = x3 - x2
        del(2) = x1 - x3
        del(3) = x2 - x1

        g(1) = x2 + x3
        g(2) = x1 + x3
        g(3) = x1 + x2

        pi(1) = x2 * x3
        pi(2) = x1 * x3
        pi(3) = x1 * x2
     
        sum1 = 0.0D+00
        do i = 1, 3
            sum1 = sum1 + y(n-1+i) * del(i) * &
                ( f / 3.0D+00 - g(i) * 0.5D+00 * e + pi(i) * feints )
        end do
     
        result = result - sum1 / ( del(1) * del(2) * del(3) )
     
        return

    end subroutine simpne 

    subroutine test_integration()

        implicit none 

        ! Local variables
        integer :: i, j, n, k, t 
        integer :: nn(11) 
        real(prec), allocatable :: zeta0(:) 
        real(prec), allocatable :: zeta(:)
        real(prec), allocatable :: var0(:) 
        real(prec), allocatable :: var(:) 
        real(prec), allocatable :: var_ints(:)
        real(prec) :: var_int 
        real(prec) :: var_int_00

        write(*,*) "=== test_integration ======"
        
        nn = [11,21,31,41,51,61,71,81,91,101,1001]

        do k = 1, size(nn)

            n = nn(k)

            allocate(zeta0(n))
            allocate(zeta(n))
            allocate(var0(n))
            allocate(var(n))
            allocate(var_ints(n))

            do i = 1, n 
                zeta0(i) = real(i-1)/real(n-1)
!                 var0(i) = real(i-1)
                var0(i)  = (n-1)-real(i-1)
            end do 

            ! Linear zeta 
            zeta = zeta0
            var = var0 

            ! Non-linear zeta 
            zeta = zeta0*zeta0 
            do i = 1, n 
                do j = 1, n 
                    if (zeta0(j) .ge. zeta(i)) exit 
                end do 

                if (zeta0(j) .eq. zeta(i)) then 
                    var(i) = var0(j) 
                else 
                    var(i) = var0(j-1) + (var0(j)-var0(j-1))*(zeta(i)-zeta0(j-1))/(zeta0(j)-zeta0(j-1))
                end if 
            end do 

!             do i = 1, n 
!                 write(*,*) zeta0(i), var0(i), zeta(i), var(i) 
!             end do 
!             stop 
            
            ! Analytical average value 
!             var_int_00 = real(n-1)/2.0

            do t = 1, 10000
                
                ! Test trapezoid1D solver 
                var_int  = integrate_trapezoid1D_pt(var,zeta)

                ! Determine "analytical" value from simpson approximation solver 
                call simpne(zeta,var,var_int_00)
            end do 

            ! Test trapezoid 1D_1D solver, check last value for full average over column
!             var_ints = integrate_trapezoid1D_1D(var,zeta)
!             var_int  = var_ints(n) 

            write(*,*) "mean (0:",n ,") = ", var_int_00, var_int, var_int-var_int_00, 100.0*(var_int-var_int_00)/var_int_00 

            deallocate(zeta0,var0,zeta,var,var_ints)

        end do 

        stop 

        return 

    end subroutine test_integration
    
end module yelmo_tools 