module yelmo_tools
    ! Generic functions and subroutines that could be used in many contexts:
    ! math, vectors, sorting, etc. 

    use yelmo_defs, only : sp, dp, wp, missing_value, TOL_UNDERFLOW, pi, &
                            io_unit_err

    !$ use omp_lib
    
    implicit none 

    private 
    public :: get_region_indices
    public :: get_neighbor_indices
    public :: calc_magnitude 
    public :: calc_magnitude_from_staggered
    public :: stagger_ac_aa
    public :: stagger_aa_ab
    public :: stagger_aa_ab_ice 
    public :: stagger_ab_aa 
    public :: stagger_ab_aa_ice
    public :: stagger_aa_acx
    public :: stagger_aa_acy
    public :: stagger_acx_aa
    public :: stagger_acy_aa
    public :: stagger_ab_acx
    public :: stagger_ab_acy 

    public :: calc_gradient_acx
    public :: calc_gradient_acy

    public :: acx_to_nodes_3D
    public :: acy_to_nodes_3D
    public :: aa_to_nodes_3D
    public :: acz_to_nodes_3D
    public :: acx_to_nodes
    public :: acy_to_nodes
    public :: aa_to_nodes
    
    public :: calc_gradient_ac
    public :: calc_gradient_ac_ice
    public :: calc_gradient_ac_gl

    public :: mean_mask
    public :: minmax

    public :: set_boundaries_2D_aa
    public :: set_boundaries_3D_aa
    public :: set_boundaries_2D_acx
    public :: set_boundaries_3D_acx
    public :: set_boundaries_2D_acy 
    public :: set_boundaries_3D_acy 

    public :: fill_borders_2D
    public :: fill_borders_3D 
    
    public :: smooth_gauss_2D
    public :: smooth_gauss_3D
    public :: gauss_values

    public :: adjust_topography_gradients 
    
    public :: regularize2D 

    ! Integration functions
    public :: test_integration
    public :: integrate_trapezoid1D_pt
    public :: integrate_trapezoid1D_1D
    public :: calc_vertical_integrated_2D
    public :: calc_vertical_integrated_3D
    
contains 

    subroutine get_region_indices(i1,i2,j1,j2,nx,ny,irange,jrange)
        ! Get indices for a region based on bounds. 
        ! If no bounds provided, use whole domain.

        implicit none

        integer, intent(OUT) :: i1
        integer, intent(OUT) :: i2
        integer, intent(OUT) :: j1
        integer, intent(OUT) :: j2
        integer, intent(IN)  :: nx
        integer, intent(IN)  :: ny
        integer, intent(IN), optional :: irange(2)
        integer, intent(IN), optional :: jrange(2)

        
        if (present(irange)) then
            i1 = irange(1)
            i2 = irange(2)
        else
            i1 = 1
            i2 = nx
        end if

        if (present(jrange)) then
            j1 = jrange(1)
            j2 = jrange(2)
        else
            j1 = 1
            j2 = ny
        end if

        return

    end subroutine get_region_indices

    subroutine get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

        implicit none

        integer, intent(OUT) :: im1 
        integer, intent(OUT) :: ip1 
        integer, intent(OUT) :: jm1 
        integer, intent(OUT) :: jp1 
        integer, intent(IN)  :: i 
        integer, intent(IN)  :: j
        integer, intent(IN)  :: nx 
        integer, intent(IN)  :: ny
        
        character(len=*), intent(IN) :: boundaries

        select case(trim(boundaries))

            case("infinite")
                im1 = max(i-1,1)
                ip1 = min(i+1,nx) 
                jm1 = max(j-1,1)
                jp1 = min(j+1,ny) 

            case("MISMIP3D","TROUGH")
                im1 = max(i-1,1)
                ip1 = min(i+1,nx) 
                jm1 = j-1
                if (jm1 .eq. 0)    jm1 = ny
                jp1 = j+1
                if (jp1 .eq. ny+1) jp1 = 1 
                
            case DEFAULT 
                ! Periodic

                im1 = i-1
                if (im1 .eq. 0)    im1 = nx 
                ip1 = i+1
                if (ip1 .eq. nx+1) ip1 = 1 

                jm1 = j-1
                if (jm1 .eq. 0)    jm1 = ny
                jp1 = j+1
                if (jp1 .eq. ny+1) jp1 = 1 

        end select 

        return

    end subroutine get_neighbor_indices

    elemental function calc_magnitude(u,v) result(umag)
        ! Get the vector magnitude from two components at the same grid location

        implicit none 

        real(wp), intent(IN)  :: u, v 
        real(wp) :: umag 

        umag = sqrt(u*u+v*v)

        return

    end function calc_magnitude
    
    function calc_magnitude_from_staggered(u,v,f_ice,boundaries) result(umag)
        ! Calculate the centered (aa-nodes) magnitude of a vector 
        ! from the staggered (ac-nodes) components

        implicit none 
        
        real(wp), intent(IN)  :: u(:,:), v(:,:)
        real(wp), intent(IN)  :: f_ice(:,:) 
        real(wp) :: umag(size(u,1),size(u,2)) 
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: im1, ip1, jm1, jp1
        real(wp) :: unow, vnow 
        real(wp) :: f1, f2, H1, H2 
        
        nx = size(u,1)
        ny = size(u,2) 

        umag = 0.0_wp 

        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
            
            if (f_ice(i,j) .eq. 1.0) then 
                unow = 0.5*(u(im1,j)+u(i,j))
                vnow = 0.5*(v(i,jm1)+v(i,j))
                
                if (abs(unow) .lt. TOL_UNDERFLOW) unow = 0.0_wp 
                if (abs(vnow) .lt. TOL_UNDERFLOW) vnow = 0.0_wp 
            else 
                unow = 0.0 
                vnow = 0.0 
            end if 

            umag(i,j) = sqrt(unow*unow+vnow*vnow)

            if (abs(umag(i,j)) .lt. TOL_UNDERFLOW) umag(i,j) = 0.0_wp

        end do 
        end do 

        return

    end function calc_magnitude_from_staggered

    function stagger_ac_aa(u,v) result(umag)
        ! Calculate the centered (aa-node) magnitude of a scalar 
        ! from the staggered (ac-node) components

        implicit none 
        
        real(wp), intent(IN)  :: u(:,:), v(:,:)    ! acx-, acy-nodes 
        real(wp) :: umag(size(u,1),size(u,2))      ! aa-nodes 

        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: im1, ip1, jm1, jp1 

        nx = size(u,1)
        ny = size(u,2) 

        umag = 0.0_wp 

        do j = 1, ny
        do i = 1, nx 
            ! BC: Periodic boundary conditions
            im1 = i-1
            if (im1 == 0) then
                im1 = nx
            end if
            ip1 = i+1
            if (ip1 == nx+1) then
                ip1 = 1
            end if

            jm1 = j-1
            if (jm1 == 0) then
                jm1 = ny
            end if
            jp1 = j+1
            if (jp1 == ny+1) then
                jp1 = 1
            end if

            umag(i,j) = 0.25_wp*(u(i,j)+u(im1,j)+v(i,j)+v(i,jm1))
        end do 
        end do 

        return

    end function stagger_ac_aa
    
    function stagger_aa_ab(u) result(ustag)
        ! Stagger from Aa => Ab
        ! Four point average from corner Aa nodes to central Ab node 

        implicit none 

        real(wp), intent(IN)  :: u(:,:) 
        real(wp) :: ustag(size(u,1),size(u,2)) 

        ! Local variables 
        integer :: i, j, nx, ny  
        integer :: im1, ip1, jm1, jp1 
        
        nx = size(u,1)
        ny = size(u,2) 

        ustag = 0.0_wp 

        do j = 1, ny 
        do i = 1, nx
            ! BC: Periodic boundary conditions
            ip1 = i+1
            if (ip1 == nx+1) then
                ip1 = 1
            end if
            jp1 = j+1
            if (jp1 == ny+1) then
                jp1 = 1
            end if

            ustag(i,j) = 0.25_wp*(u(ip1,jp1)+u(ip1,j)+u(i,jp1)+u(i,j))
        end do 
        end do 

        return

    end function stagger_aa_ab
    
    function stagger_aa_ab_ice(u,H_ice,f_ice) result(ustag)
        ! Stagger from Aa => Ab
        ! Four point average from corner Aa nodes to central Ab node 

        implicit none 

        real(wp), intent(IN)  :: u(:,:) 
        real(wp), intent(IN)  :: H_ice(:,:) 
        real(wp), intent(IN)  :: f_ice(:,:) 
        real(wp) :: ustag(size(u,1),size(u,2)) 

        ! Local variables 
        integer :: i, j, nx, ny, k   
        integer :: im1, ip1, jm1, jp1 

        nx = size(u,1)
        ny = size(u,2) 

        ustag = 0.0_wp 

        do j = 1, ny 
        do i = 1, nx

            ! BC: Periodic boundary conditions
            ip1 = i+1
            if (ip1 == nx+1) then
                ip1 = 1
            end if
            jp1 = j+1
            if (jp1 == ny+1) then
                jp1 = 1
            end if

            k = 0 
            ustag(i,j) = 0.0 
            if (f_ice(i,j) .eq. 0.0) then 
                ustag(i,j) = ustag(i,j) + u(i,j) 
                k = k+1
            end if 

            if (f_ice(ip1,j) .eq. 0.0) then 
                ustag(i,j) = ustag(i,j) + u(ip1,j) 
                k = k+1 
            end if 
            
            if (f_ice(i,jp1) .eq. 0.0) then 
                ustag(i,j) = ustag(i,j) + u(i,jp1) 
                k = k+1 
            end if 
            
            if (f_ice(ip1,jp1) .eq. 0.0) then 
                ustag(i,j) = ustag(i,j) + u(ip1,jp1) 
                k = k+1 
            end if 
            
            if (k .gt. 0) then 
                ustag(i,j) = ustag(i,j) / real(k,wp)
            else 
                ustag(i,j) = 0.25_wp*(u(ip1,jp1)+u(ip1,j)+u(i,jp1)+u(i,j))
            end if 

        end do 
        end do 

        return

    end function stagger_aa_ab_ice
    
    function stagger_ab_aa(u) result(ustag)
        ! Stagger from Ab => Aa
        ! Four point average from corner Ab nodes to central Aa node 

        implicit none 

        real(wp), intent(IN)  :: u(:,:) 
        real(wp) :: ustag(size(u,1),size(u,2)) 

        ! Local variables 
        integer :: i, j, nx, ny  
        integer :: im1, jm1, ip1, jp1

        nx = size(u,1)
        ny = size(u,2) 

        ustag = 0.0_wp 

        do j = 1, ny 
        do i = 1, nx
            ! BC: Periodic boundary conditions
            im1 = i-1
            if (im1 == 0) then
                im1 = nx
            end if
            jm1 = j-1
            if (jm1 == 0) then
                jm1 = ny
            end if

            ustag(i,j) = 0.25_wp*(u(i,j)+u(im1,j)+u(i,jm1)+u(im1,jm1))
        end do 
        end do 

        return

    end function stagger_ab_aa
    
    function stagger_ab_aa_ice(u,H_ice,f_ice) result(ustag)
        ! Stagger from ab => aa
        ! Four point average from corner ab-nodes to central aa-node 

        implicit none 

        real(wp), intent(IN)  :: u(:,:) 
        real(wp), intent(IN)  :: H_ice(:,:) 
        real(wp), intent(IN)  :: f_ice(:,:) 
        real(wp) :: ustag(size(u,1),size(u,2)) 

        ! Local variables 
        integer :: i, j, nx, ny, k   
        integer :: im1, jm1, ip1, jp1 
        real(wp), allocatable ::H_ice_ab(:,:) 

        nx = size(u,1)
        ny = size(u,2) 

        allocate(H_ice_ab(nx,ny))

        H_ice_ab = stagger_aa_ab_ice(H_ice,H_ice,f_ice) 

        ustag = 0.0_wp 

        do j = 1, ny 
        do i = 1, nx

            ! BC: Periodic boundary conditions
            im1 = i-1
            if (im1 == 0) then
                im1 = nx
            end if

            jm1 = j-1
            if (jm1 == 0) then
                jm1 = ny
            end if

            k = 0 
            ustag(i,j) = 0.0 
            if (H_ice_ab(i,j) .gt. 0.0) then 
                ustag(i,j) = ustag(i,j) + u(i,j) 
                k = k+1
            end if 

            if (H_ice_ab(im1,j) .gt. 0.0) then 
                ustag(i,j) = ustag(i,j) + u(im1,j) 
                k = k+1 
            end if 
            
            if (H_ice_ab(i,jm1) .gt. 0.0) then 
                ustag(i,j) = ustag(i,j) + u(i,jm1) 
                k = k+1 
            end if 
            
            if (H_ice_ab(im1,jm1) .gt. 0.0) then 
                ustag(i,j) = ustag(i,j) + u(im1,jm1) 
                k = k+1 
            end if 
            
            if (k .gt. 0) then 
                ustag(i,j) = ustag(i,j) / real(k,wp)
            else 
                ustag(i,j) = 0.25_wp*(u(im1,jm1)+u(im1,j)+u(i,jm1)+u(i,j))
            end if 

        end do 
        end do 

        return

    end function stagger_ab_aa_ice

! ===== NEW =======================

    subroutine acx_to_nodes_3D(varn,varx,i,j,k,xn,yn,zn,im1,ip1,jm1,jp1,mask)
        ! 3D !!
        ! Variable defined on acx-nodes horizontally, aa-nodes vertically. 

        real(wp), intent(OUT) :: varn(:) 
        real(wp), intent(IN)  :: varx(:,:,:) 
        integer,  intent(IN)  :: i 
        integer,  intent(IN)  :: j
        integer,  intent(IN)  :: k 
        real(wp), intent(IN)  :: xn(:)      ! Four points in horizontal plane
        real(wp), intent(IN)  :: yn(:)      ! Four points in horizontal plane
        real(wp), intent(IN)  :: zn         ! Distance to horizontal planes in vertical direction
        integer,  intent(IN)  :: im1, ip1, jm1, jp1 
        logical,  intent(IN), optional :: mask(:,:)

        ! Local variables 
        integer  :: nx, ny, nz_aa
        real(wp) :: wt
        real(wp) :: vx_up(4)
        real(wp) :: vx_md(4)
        real(wp) :: vx_dn(4) 

        nx    = size(varx,1)
        ny    = size(varx,2)
        nz_aa = size(varx,3) 

        ! Get four nodes in horizontal aa-node (middle) plane 
        call acx_to_nodes(vx_md,varx(:,:,k),i,j,xn,yn,im1,ip1,jm1,jp1,mask)

        ! Get four nodes in plane below
        if (k .gt. 1) then 
            call acx_to_nodes(vx_dn,varx(:,:,k-1),i,j,xn,yn,im1,ip1,jm1,jp1,mask)
            wt = zn/2.0
            varn(1:4) = wt*vx_dn + (1.0-wt)*vx_md
        else 
            varn(1:4) = vx_md 
        end if 

        ! Get four nodes in plane above
        if (k .lt. nz_aa) then 
            call acx_to_nodes(vx_up,varx(:,:,k+1),i,j,xn,yn,im1,ip1,jm1,jp1,mask)
            wt = zn/2.0
            varn(5:8) = wt*vx_up + (1.0-wt)*vx_md
        else 
            varn(5:8) = vx_md 
        end if 

        return

    end subroutine acx_to_nodes_3D

    subroutine acy_to_nodes_3D(varn,vary,i,j,k,xn,yn,zn,im1,ip1,jm1,jp1,mask)
        ! 3D !!
        ! Variable defined on acx-nodes horizontally, aa-nodes vertically. 

        real(wp), intent(OUT) :: varn(:) 
        real(wp), intent(IN)  :: vary(:,:,:) 
        integer,  intent(IN)  :: i 
        integer,  intent(IN)  :: j
        integer,  intent(IN)  :: k 
        real(wp), intent(IN)  :: xn(:)      ! Four points in horizontal plane
        real(wp), intent(IN)  :: yn(:)      ! Four points in horizontal plane
        real(wp), intent(IN)  :: zn         ! Distance to horizontal planes in vertical direction
        integer,  intent(IN)  :: im1, ip1, jm1, jp1 
        logical,  intent(IN), optional :: mask(:,:)

        ! Local variables 
        integer  :: nx, ny, nz_aa
        real(wp) :: wt
        real(wp) :: vy_up(4)
        real(wp) :: vy_md(4)
        real(wp) :: vy_dn(4) 

        nx    = size(vary,1)
        ny    = size(vary,2)
        nz_aa = size(vary,3)

        ! Get four nodes in horizontal aa-node (middle) plane 
        call acy_to_nodes(vy_md,vary(:,:,k),i,j,xn,yn,im1,ip1,jm1,jp1,mask)

        ! Get four nodes in plane below
        if (k .gt. 1) then 
            call acy_to_nodes(vy_dn,vary(:,:,k-1),i,j,xn,yn,im1,ip1,jm1,jp1,mask)
            wt = zn/2.0
            varn(1:4) = wt*vy_dn + (1.0-wt)*vy_md
        else 
            varn(1:4) = vy_md 
        end if 

        ! Get four nodes in plane above
        if (k .lt. nz_aa) then 
            call acy_to_nodes(vy_up,vary(:,:,k+1),i,j,xn,yn,im1,ip1,jm1,jp1,mask)
            wt = zn/2.0
            varn(5:8) = wt*vy_up + (1.0-wt)*vy_md
        else 
            varn(5:8) = vy_md 
        end if 

        return

    end subroutine acy_to_nodes_3D
    
    subroutine aa_to_nodes_3D(varn,var,i,j,k,xn,yn,zn,im1,ip1,jm1,jp1,mask)
        ! 3D !!
        ! Variable defined on aa-nodes horizontally, aa-nodes vertically. 

        real(wp), intent(OUT) :: varn(:) 
        real(wp), intent(IN)  :: var(:,:,:) 
        integer,  intent(IN)  :: i 
        integer,  intent(IN)  :: j
        integer,  intent(IN)  :: k 
        real(wp), intent(IN)  :: xn(:)      ! Four points in horizontal plane
        real(wp), intent(IN)  :: yn(:)      ! Four points in horizontal plane
        real(wp), intent(IN)  :: zn         ! Distance to horizontal planes in vertical direction
        integer,  intent(IN)  :: im1, ip1, jm1, jp1 
        logical,  intent(IN), optional :: mask(:,:)

        ! Local variables 
        integer  :: nx, ny, nz_aa
        real(wp) :: wt
        real(wp) :: vv_up(4)
        real(wp) :: vv_md(4)
        real(wp) :: vv_dn(4) 

        nx    = size(var,1)
        ny    = size(var,2)
        nz_aa = size(var,3) 

        ! Get four nodes in horizontal aa-node (middle) plane 
        call aa_to_nodes(vv_md,var(:,:,k),i,j,xn,yn,im1,ip1,jm1,jp1,mask)

        ! Get four nodes in plane below
        if (k .gt. 1) then 
            call aa_to_nodes(vv_dn,var(:,:,k-1),i,j,xn,yn,im1,ip1,jm1,jp1,mask)
            wt = zn/2.0
            varn(1:4) = wt*vv_dn + (1.0-wt)*vv_md
        else 
            varn(1:4) = vv_md 
        end if 

        ! Get four nodes in plane above
        if (k .lt. nz_aa) then 
            call aa_to_nodes(vv_up,var(:,:,k+1),i,j,xn,yn,im1,ip1,jm1,jp1,mask)
            wt = zn/2.0
            varn(5:8) = wt*vv_up + (1.0-wt)*vv_md
        else 
            varn(5:8) = vv_md 
        end if 

        return

    end subroutine aa_to_nodes_3D
    
    subroutine acz_to_nodes_3D(varn,var,i,j,k,xn,yn,zn,im1,ip1,jm1,jp1,mask)
        ! 3D !!
        ! Variable defined on aa-nodes horizontally, acz-nodes vertically. 

        real(wp), intent(OUT) :: varn(:) 
        real(wp), intent(IN)  :: var(:,:,:) 
        integer,  intent(IN)  :: i 
        integer,  intent(IN)  :: j
        integer,  intent(IN)  :: k 
        real(wp), intent(IN)  :: xn(:)      ! Four points in horizontal plane
        real(wp), intent(IN)  :: yn(:)      ! Four points in horizontal plane
        real(wp), intent(IN)  :: zn         ! Distance to horizontal planes in vertical direction
        integer,  intent(IN)  :: im1, ip1, jm1, jp1 
        logical,  intent(IN), optional :: mask(:,:)

        ! Local variables 
        integer  :: nx, ny, nz_ac
        real(wp) :: wt
        real(wp) :: vv_up(4)
        real(wp) :: vv_md(4)
        real(wp) :: vv_dn(4) 

        nx    = size(var,1)
        ny    = size(var,2)
        nz_ac = size(var,3)

        ! Get four nodes in horizontal plane of top cell face
        call aa_to_nodes(vv_up,var(:,:,k),i,j,xn,yn,im1,ip1,jm1,jp1,mask)
        
        ! Get four nodes in horizontal plane of bottom cell face
        if (k .gt. 1) then 
            call aa_to_nodes(vv_dn,var(:,:,k-1),i,j,xn,yn,im1,ip1,jm1,jp1,mask)
        else 
            vv_dn = vv_up 
        end if 

        ! Interpolate to upper and lower nodes
        wt = (1.0-zn) / 2.0
        varn(1:4) = wt*vv_dn + (1.0-wt)*vv_up 
        varn(5:8) = (1.0-wt)*vv_dn + wt*vv_up 

        return

    end subroutine acz_to_nodes_3D
    
    subroutine acx_to_nodes(varn,varx,i,j,xn,yn,im1,ip1,jm1,jp1,mask)
        ! 2D !!
        ! Variable defined on acx-nodes horizontally 
        ! Assumed no vertical dimension. 

        real(wp), intent(OUT) :: varn(:) 
        real(wp), intent(IN)  :: varx(:,:) 
        integer,  intent(IN)  :: i 
        integer,  intent(IN)  :: j
        real(wp), intent(IN)  :: xn(:)
        real(wp), intent(IN)  :: yn(:) 
        integer,  intent(IN)  :: im1, ip1, jm1, jp1 
        logical,  intent(IN), optional :: mask(:,:)

        ! Local variables 
        integer  :: k, nx, ny, n 
        integer  :: i0, i1, j0, j1 
        real(wp) :: v0, v1, wt 
        real(wp) :: va, vb, vc, vd

        n = size(xn,1)
        
        nx = size(varx,1)
        ny = size(varx,2)

        ! Initialize interpolated variable at nodes of interest
        varn = 0.0

        do k = 1, n

            if (yn(k) > 0) then
                j0 = j
                j1 = jp1
                wt = yn(k) / 2.0 
            else
                j0 = jm1
                j1 = j
                wt = (1.0 + (1.0+yn(k)) ) / 2.0 
            end if

            if (present(mask)) then
                ! Limit neighbors to within mask

                write(io_unit_err,*) "acx_to_nodes:: Error: mask is not yet working properly. &
                &Do not include a mask argument."
                stop

                ! Get left and right-side 
                va = varx(im1,j0)
                if (.not. mask(im1,j0)) va = varx(i,j)
                vb = varx(im1,j1)
                if (.not. mask(im1,j1)) vb = varx(i,j)
                vc = varx(i,j0)
                if (.not. mask(i,j0)) vc = varx(i,j)
                vd = varx(i,j1)
                if (.not. mask(i,j1)) vd = varx(i,j)
                
                ! Get left and right-side 
                v0 = (1.0-wt)*va + wt*vb
                v1 = (1.0-wt)*vc + wt*vd
            else
                ! Use all values

                ! Get left and right-side 
                v0 = (1.0-wt)*varx(im1,j0) + wt*varx(im1,j1)
                v1 = (1.0-wt)*varx(i,j0)   + wt*varx(i,j1)
            end if
            
            ! Interpolate horizontally to the node location 
            wt = (1.0 + xn(k)) / 2.0
            varn(k) = (1.0-wt)*v0 + wt*v1

        end do
        
        return

    end subroutine acx_to_nodes

    subroutine acy_to_nodes(varn,vary,i,j,xn,yn,im1,ip1,jm1,jp1,mask)
        ! 2D !!
        ! Variable defined on acy-nodes horizontally 
        ! Assumed no vertical dimension. 
        
        real(wp), intent(OUT) :: varn(:) 
        real(wp), intent(IN)  :: vary(:,:) 
        integer,  intent(IN)  :: i 
        integer,  intent(IN)  :: j
        real(wp), intent(IN)  :: xn(:)
        real(wp), intent(IN)  :: yn(:) 
        integer,  intent(IN)  :: im1, ip1, jm1, jp1
        logical,  intent(IN), optional :: mask(:,:)

        ! Local variables 
        integer  :: k, nx, ny, n 
        integer  :: i0, i1, j0, j1  
        real(wp) :: v0, v1, wt 
        real(wp) :: va, vb, vc, vd

        n = size(xn,1)
        
        nx = size(vary,1)
        ny = size(vary,2)

        ! Initialize interpolated variable at nodes of interest
        varn = 0.0

        do k = 1, n

            if (xn(k) > 0) then
                i0 = i
                i1 = ip1
                wt = xn(k) / 2.0 
            else
                i0 = im1
                i1 = i
                wt = (1.0 + (1.0+xn(k)) ) / 2.0 
            end if

            if (present(mask)) then
                ! Limit neighbors to within mask

                write(io_unit_err,*) "acx_to_nodes:: Error: mask is not yet working properly. &
                &Do not include a mask argument."
                stop

                ! Get left and right-side 
                va = vary(i0,jm1)
                if (.not. mask(i0,jm1)) va = vary(i,j)
                vb = vary(i1,jm1)
                if (.not. mask(i1,jm1)) vb = vary(i,j)
                vc = vary(i0,j)
                if (.not. mask(i0,j)) vc = vary(i,j)
                vd = vary(i1,j)
                if (.not. mask(i1,j)) vd = vary(i,j)
                
                ! Get left and right-side 
                v0 = (1.0-wt)*va + wt*vb
                v1 = (1.0-wt)*vc + wt*vd
            else
                ! Use all values

                ! Get top and bottom-side 
                v0 = (1.0-wt)*vary(i0,jm1) + wt*vary(i1,jm1);
                v1 = (1.0-wt)*vary(i0,j)   + wt*vary(i1,j);
            end if
                
            ! Interpolate vertically to the node location 
            wt = (1.0 + yn(k)) / 2.0;
            varn(k) = (1.0-wt)*v0 + wt*v1;

        end do
        
        return

    end subroutine acy_to_nodes
    
    subroutine aa_to_nodes(varn,var,i,j,xn,yn,im1,ip1,jm1,jp1,mask)
        ! 2D !!
        ! Variable defined on aa-nodes horizontally 
        ! Assumed no vertical dimension. 

        real(wp), intent(OUT) :: varn(:) 
        real(wp), intent(IN)  :: var(:,:) 
        integer,  intent(IN)  :: i 
        integer,  intent(IN)  :: j
        real(wp), intent(IN)  :: xn(:)
        real(wp), intent(IN)  :: yn(:) 
        integer,  intent(IN)  :: im1, ip1, jm1, jp1 
        logical,  intent(IN), optional :: mask(:,:)

        ! Local variables 
        integer  :: k, nx, ny, n 
        integer  :: i0, i1, j0, j1 
        real(wp) :: v0, v1, wt 
        real(wp) :: va, vb, vc, vd

        n = size(xn,1)
        
        nx = size(var,1)
        ny = size(var,2)

        ! Initialize interpolated variable at nodes of interest
        varn = 0.0

        do k = 1, n

            if (yn(k) > 0) then
                j0 = j
                j1 = jp1
                wt = yn(k) / 2.0 
            else
                j0 = jm1
                j1 = j
                wt = (1.0 + (1.0+yn(k)) ) / 2.0 
            end if

            if (xn(k) > 0) then
                i0 = i
                i1 = ip1 
            else
                i0 = im1
                i1 = i
            end if

            if (present(mask)) then
                ! Limit neighbors to within mask

                write(io_unit_err,*) "acx_to_nodes:: Error: mask is not yet working properly. &
                &Do not include a mask argument."
                stop
                
                ! Get left and right-side 
                va = var(i0,j0)
                if (.not. mask(i0,j0)) va = var(i,j)
                vb = var(i0,j1)
                if (.not. mask(i0,j1)) vb = var(i,j)
                vc = var(i1,j0)
                if (.not. mask(i1,j0)) vc = var(i,j)
                vd = var(i1,j1)
                if (.not. mask(i1,j1)) vd = var(i,j)
                
                ! Get left and right-side 
                v0 = (1.0-wt)*va + wt*vb
                v1 = (1.0-wt)*vc + wt*vd

            else
                ! Use all values

                ! Get left and right-side 
                v0 = (1.0-wt)*var(i0,j0) + wt*var(i0,j1)
                v1 = (1.0-wt)*var(i1,j0) + wt*var(i1,j1)
            end if

            ! Interpolate horizontally to the node location 
            wt = abs(xn(k)) / 2.0
            varn(k) = (1.0-wt)*v0 + wt*v1

        end do
        
        return

    end subroutine aa_to_nodes

! ===== end NEW =======================



    function stagger_aa_acx(u) result(ustag)
        ! Stagger from Aa => Ac, x-direction 

        implicit none

        real(wp), intent(IN)  :: u(:,:) 
        real(wp) :: ustag(size(u,1),size(u,2)) 

        ! Local variables 
        integer :: i, j, nx, ny  

        nx = size(u,1)
        ny = size(u,2) 

        ustag = 0.0_wp 

        do j = 1, ny 
        do i = 1, nx-1
            ustag(i,j) = 0.5_wp*(u(i,j)+u(i+1,j))
        end do 
        end do 

        return

    end function stagger_aa_acx
    
    function stagger_aa_acy(u) result(ustag)
        ! Stagger from Aa => Ac 

        implicit none 

        real(wp), intent(IN)  :: u(:,:) 
        real(wp) :: ustag(size(u,1),size(u,2)) 

        ! Local variables 
        integer :: i, j, nx, ny  

        nx = size(u,1)
        ny = size(u,2) 

        ustag = 0.0_wp 

        do j = 1, ny-1 
        do i = 1, nx
            ustag(i,j) = 0.5_wp*(u(i,j)+u(i,j+1))
        end do 
        end do 

        return

    end function stagger_aa_acy
    
    function stagger_acx_aa(u) result(ustag)
        ! Stagger from Aa => Ac, x-direction 

        implicit none

        real(wp), intent(IN)  :: u(:,:) 
        real(wp) :: ustag(size(u,1),size(u,2)) 

        ! Local variables 
        integer :: i, j, nx, ny  

        nx = size(u,1)
        ny = size(u,2) 

        ustag = 0.0_wp 

        do j = 1, ny 
        do i = 2, nx
            ustag(i,j) = 0.5_wp*(u(i-1,j)+u(i,j))
        end do 
        end do 

        return

    end function stagger_acx_aa
    
    function stagger_acy_aa(u) result(ustag)
        ! Stagger from Aa => Ac 

        implicit none 

        real(wp), intent(IN)  :: u(:,:) 
        real(wp) :: ustag(size(u,1),size(u,2)) 

        ! Local variables 
        integer :: i, j, nx, ny  

        nx = size(u,1)
        ny = size(u,2) 

        ustag = 0.0_wp 

        do j = 2, ny
        do i = 1, nx
            ustag(i,j) = 0.5_wp*(u(i,j-1)+u(i,j))
        end do 
        end do 

        return

    end function stagger_acy_aa
    
    function stagger_ab_acx(u) result(ustag)
        ! Stagger from Ab => Ac, x-direction 

        implicit none

        real(wp), intent(IN)  :: u(:,:) 
        real(wp) :: ustag(size(u,1),size(u,2)) 

        ! Local variables 
        integer :: i, j, nx, ny  

        nx = size(u,1)
        ny = size(u,2) 

        ustag = 0.0_wp 

        do j = 2, ny 
        do i = 1, nx
            ustag(i,j) = 0.5_wp*(u(i,j)+u(i,j-1))
        end do 
        end do 

        return

    end function stagger_ab_acx
    
    function stagger_ab_acy(u) result(ustag)
        ! Stagger from Ab => Ac 

        implicit none 

        real(wp), intent(IN)  :: u(:,:) 
        real(wp) :: ustag(size(u,1),size(u,2)) 

        ! Local variables 
        integer :: i, j, nx, ny  

        nx = size(u,1)
        ny = size(u,2) 

        ustag = 0.0_wp 

        do j = 1, ny 
        do i = 2, nx
            ustag(i,j) = 0.5_wp*(u(i,j)+u(i-1,j))
        end do 
        end do 

        return

    end function stagger_ab_acy
    
    subroutine calc_gradient_ac(dvardx,dvardy,var,dx)
        ! Calculate gradient on ac nodes 

        implicit none 

        real(wp), intent(OUT) :: dvardx(:,:) 
        real(wp), intent(OUT) :: dvardy(:,:) 
        real(wp), intent(IN)  :: var(:,:) 
        real(wp), intent(IN)  :: dx 

        ! Local variables 
        integer    :: i, j, nx, ny
        integer    :: im1, ip1, jm1, jp1  
        real(wp) :: dy 

        nx = size(var,1)
        ny = size(var,2)

        ! Assume y-resolution is identical to x-resolution 
        dy = dx 

        do j = 1, ny 
        do i = 1, nx 

            im1 = max(1, i-1)
            ip1 = min(nx,i+1)
            
            jm1 = max(1, j-1)
            jp1 = min(ny,j+1)

            ! Slope in x-direction
            dvardx(i,j) = (var(ip1,j)-var(i,j))/dx 

            ! Slope in y-direction
            dvardy(i,j) = (var(i,jp1)-var(i,j))/dy
        end do 
        end do 

        return 

    end subroutine calc_gradient_ac
    
    subroutine calc_gradient_acx(dvardx,var,f_ice,dx,grad_lim,margin2nd,zero_outside,boundaries)
        ! Calculate gradient on ac-nodes, accounting for ice margin if needed

        implicit none 

        real(wp), intent(OUT) :: dvardx(:,:) 
        real(wp), intent(IN)  :: var(:,:) 
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: dx 
        real(wp), intent(IN)  :: grad_lim 
        logical,  intent(IN)  :: margin2nd 
        logical,  intent(IN)  :: zero_outside 
        character(len=*), intent(IN) :: boundaries  ! Boundary conditions to apply 
        
        ! Local variables 
        integer  :: i, j, nx, ny 
        integer  :: im1, ip1, jm1, jp1
        integer  :: im2, ip2, jm2, jp2
        real(wp) :: V0, V1, V2 

        nx = size(var,1)
        ny = size(var,2)

        !!$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,ip2,V0,V1,V2)
        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
            
            V0 = var(i,j) 
            V1 = var(ip1,j) 

            if (zero_outside) then 

                if (f_ice(i,j)   .lt. 1.0) V0 = 0.0 
                if (f_ice(ip1,j) .lt. 1.0) V1 = 0.0 
                
            end if 

            dvardx(i,j) = (V1-V0)/dx 

if (margin2nd) then 
            ! === Modify margin gradients =========================
            ! Following Saito et al (2007) by applying a second-order, upwind gradient

            if (f_ice(i,j) .eq. 1.0 .and. f_ice(ip1,j) .lt. 1.0) then 
                ! Ice-free to the right

                if (f_ice(im1,j) .eq. 1.0) then 
                    V0 = var(ip1,j)
                    if (zero_outside) V0 = 0.0 
                    V1 = var(i,j)
                    V2 = var(im1,j)
                    dvardx(i,j) = (1.0*V2-4.0*V1+3.0*V0)/dx
                else 
                    dvardx(i,j) = 0.0
                end if 

            else if (f_ice(i,j) .lt. 1.0 .and. f_ice(ip1,j) .eq. 1.0) then
                ! Ice-free to the left

                if (ip1 .lt. nx) then 
                    ip2 = ip1+1 
                    if (f_ice(ip2,j) .eq. 1.0) then
                        V0 = var(i,j)
                        if (zero_outside) V0 = 0.0 
                        V1 = var(ip1,j)
                        V2 = var(ip2,j)
                        dvardx(i,j) = -(1.0*V2-4.0*V1+3.0*V0)/dx
                    else 
                        dvardx(i,j) = 0.0
                    end if
                else
                    dvardx(i,j) = 0.0
                end if

            end if 

end if

        end do 
        end do
        !!$omp end parallel do

        ! Special case for infinite boundary conditions - ensure that slope 
        ! is the same, not the variable itself.
        select case(trim(boundaries))

            case("infinite","MISMIP3D","TROUGH") 
                dvardx(1,:)  = dvardx(2,:)
                dvardx(nx,:) = dvardx(nx-1,:)
            
        end select

        ! Finally, ensure that gradient is beneath desired limit 
        call minmax(dvardx,grad_lim)

        return 

    end subroutine calc_gradient_acx
    
subroutine calc_gradient_acy(dvardy,var,f_ice,dy,grad_lim,margin2nd,zero_outside,boundaries)
        ! Calculate gradient on ac-nodes, accounting for ice margin if needed

        implicit none 

        real(wp), intent(OUT) :: dvardy(:,:) 
        real(wp), intent(IN)  :: var(:,:) 
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: dy 
        real(wp), intent(IN)  :: grad_lim 
        logical,  intent(IN)  :: margin2nd 
        logical,  intent(IN)  :: zero_outside 
        character(len=*), intent(IN) :: boundaries  ! Boundary conditions to apply 
        
        ! Local variables 
        integer  :: i, j, nx, ny 
        integer  :: im1, ip1, jm1, jp1
        integer  :: im2, ip2, jm2, jp2
        real(wp) :: V0, V1, V2 

        nx = size(var,1)
        ny = size(var,2)

        !!$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,jp2,V0,V1,V2)
        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
             
            V0 = var(i,j) 
            V1 = var(i,jp1) 

            if (zero_outside) then 

                if (f_ice(i,j)   .lt. 1.0) V0 = 0.0 
                if (f_ice(i,jp1) .lt. 1.0) V1 = 0.0 
                
            end if 

            dvardy(i,j) = (V1-V0)/dy

if (margin2nd) then 
            ! === Modify margin gradients =========================
            ! Following Saito et al (2007) by applying a second-order, upwind gradient

            if (f_ice(i,j) .eq. 1.0 .and. f_ice(i,jp1) .lt. 1.0) then 
                ! Ice-free to the top

                if (f_ice(i,jm1) .eq. 1.0) then 
                    V0 = var(i,jp1)
                    if (zero_outside) V0 = 0.0 
                    V1 = var(i,j)
                    V2 = var(i,jm1)
                    dvardy(i,j) = (1.0*V2-4.0*V1+3.0*V0)/dy
                else 
                    dvardy(i,j) = 0.0
                end if 

            else if (f_ice(i,j) .lt. 1.0 .and. f_ice(i,jp1) .eq. 1.0) then
                ! Ice-free to the bottom

                if (jp1 .lt. ny) then 
                    jp2 = jp1+1 
                    if (f_ice(i,jp2) .eq. 1.0) then
                        V0 = var(i,j)
                        if (zero_outside) V0 = 0.0 
                        V1 = var(i,jp1)
                        V2 = var(i,jp2)
                        dvardy(i,j) = -(1.0*V2-4.0*V1+3.0*V0)/dy
                    else 
                        dvardy(i,j) = 0.0
                    end if
                else
                    dvardy(i,j) = 0.0
                end if

            end if 

end if

        end do 
        end do
        !!$omp end parallel do
        
        ! Special case for infinite boundary conditions - ensure that slope 
        ! is the same, not the variable itself.
        select case(trim(boundaries))

            case("infinite") 
                dvardy(:,1)  = dvardy(:,2)
                dvardy(:,ny) = dvardy(:,ny-1)

        end select

        ! Finally, ensure that gradient is beneath desired limit 
        call minmax(dvardy,grad_lim)

        return 

    end subroutine calc_gradient_acy
    
    subroutine calc_gradient_ac_ice(dvardx,dvardy,var,f_ice,dx,margin2nd,grad_lim,boundaries,zero_outside)
        ! Calculate gradient on ac nodes, accounting for ice margin if needed

        implicit none 

        real(wp), intent(OUT) :: dvardx(:,:) 
        real(wp), intent(OUT) :: dvardy(:,:) 
        real(wp), intent(IN)  :: var(:,:) 
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: dx 
        logical,    intent(IN)  :: margin2nd 
        real(wp), intent(IN)  :: grad_lim 
        character(len=*), intent(IN) :: boundaries  ! Boundary conditions to apply 
        logical,    intent(IN), optional :: zero_outside 

        ! Local variables 
        integer    :: i, j, nx, ny 
        integer    :: im1, ip1, jm1, jp1 
        real(wp) :: dy 
        real(wp) :: H0, H1, H2 
        logical    :: use_zeros_outside 


        use_zeros_outside = .FALSE. 
        if (present(zero_outside)) use_zeros_outside = zero_outside 

        nx = size(var,1)
        ny = size(var,2)

        ! Assume y-resolution is identical to x-resolution 
        dy = dx 

        do j = 1, ny 
        do i = 1, nx 

            ip1 = min(nx,i+1)
            jp1 = min(ny,j+1)

            ! x-direction 
            H0 = var(i,j) 
            H1 = var(ip1,j) 

            if (use_zeros_outside) then 

                if (f_ice(i,j)   .lt. 1.0) H0 = 0.0 
                if (f_ice(ip1,j) .lt. 1.0) H1 = 0.0 
                
            end if 

            dvardx(i,j) = (H1-H0)/dx 


            ! y-direction 
            H0 = var(i,j) 
            H1 = var(i,jp1) 

            if (use_zeros_outside) then 

                if (f_ice(i,j)   .lt. 1.0) H0 = 0.0 
                if (f_ice(i,jp1) .lt. 1.0) H1 = 0.0 
                
            end if 

            dvardy(i,j) = (H1-H0)/dy

        end do 
        end do 

        ! === Modify margin gradients =========================
        ! Following Saito et al (2007) by applying a second-order, upwind gradient

        if (margin2nd) then

            ! Slope in x-direction
            do j = 1, ny 
            do i = 2, nx-3

                if (f_ice(i,j) .eq. 1.0 .and. f_ice(i+1,j) .lt. 1.0) then 
                    ! Margin point (ice-free to the right)

                    H0 = var(i,j) 
                    H1 = var(i-1,j)

                    ! Second-order upwind
                    dvardx(i,j) = -(1.0/(3.0*dx))*(H0+H1)

                else if (f_ice(i,j) .lt. 1.0 .and. f_ice(i+1,j) .eq. 1.0) then
                    ! Margin point (ice-free to the left)

                    H0 = var(i+1,j) 
                    H1 = var(i+2,j)

                    ! Second-order upwind
                    dvardx(i,j) = (1.0/(3.0*dx))*(H0+H1)

                end if 

            end do 
            end do 
            
            ! Slope in y-direction
            do j = 2, ny-3 
            do i = 1, nx

                if (f_ice(i,j) .eq. 1.0 .and. f_ice(i,j+1) .lt. 1.0) then 
                    ! Margin point (ice-free above)

                    H0 = var(i,j) 
                    H1 = var(i,j-1)

                    ! Second-order upwind
                    dvardy(i,j) = -(1.0/(3.0*dy))*(H0+H1)

                else if (f_ice(i,j) .lt. 1.0 .and. f_ice(i,j+1) .eq. 1.0) then
                    ! Margin point (ice-free below)

                    H0 = var(i,j+1) 
                    H1 = var(i,j+2)

                    ! Second-order upwind
                    dvardy(i,j) = (1.0/(3.0*dy))*(H0+H1)

                end if 

            end do 
            end do

        end if 

        ! Finally, ensure that gradient is beneath desired limit 
        call minmax(dvardx,grad_lim)
        call minmax(dvardy,grad_lim)

        if (trim(boundaries) .eq. "periodic") then 

            dvardx(1,:)    = dvardx(nx-2,:) 
            dvardx(nx-1,:) = dvardx(2,:) 
            dvardx(nx,:)   = dvardx(3,:) 
            dvardx(:,1)    = dvardx(:,ny-1)
            dvardx(:,ny)   = dvardx(:,2) 

            dvardy(1,:)    = dvardy(nx-1,:) 
            dvardy(nx,:)   = dvardy(2,:) 
            dvardy(:,1)    = dvardy(:,ny-2)
            dvardy(:,ny-1) = dvardy(:,2) 
            dvardy(:,ny)   = dvardy(:,3)

        else if (trim(boundaries) .eq. "infinite") then 

            dvardx(1,:)    = dvardx(2,:) 
            dvardx(nx-1,:) = dvardx(nx-2,:) 
            dvardx(nx,:)   = dvardx(nx-2,:) 
            dvardx(:,1)    = dvardx(:,2)
            dvardx(:,ny)   = dvardx(:,ny-1) 

            dvardy(1,:)    = dvardy(2,:) 
            dvardy(nx,:)   = dvardy(nx-1,:) 
            dvardy(:,1)    = dvardy(:,2)
            dvardy(:,ny-1) = dvardy(:,ny-2) 
            dvardy(:,ny)   = dvardy(:,ny-2)

        end if 

        return 

    end subroutine calc_gradient_ac_ice
    
    subroutine calc_gradient_ac_gl(dvardx,dvardy,var,H_ice, &
                                      f_grnd_acx,f_grnd_acy,dx,method,grad_lim)

        implicit none 

        real(wp), intent(OUT) :: dvardx(:,:)
        real(wp), intent(OUT) :: dvardy(:,:) 
        real(wp), intent(IN)  :: var(:,:) 
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: f_grnd_acx(:,:)
        real(wp), intent(IN)  :: f_grnd_acy(:,:)
        real(wp), intent(IN)  :: dx 
        integer,    intent(IN)  :: method           ! Which gl gradient calculation to use
        real(wp), intent(IN)  :: grad_lim         ! Very high limit == 0.05, low limit < 0.01 

        ! Local variables 
        integer :: i, j, nx, ny 
        real(wp) :: dy
        real(wp) :: dvardx_1, dvardx_2 

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

    subroutine find_upstream_neighbor(i0,j0,i,j,ux,uy)
        ! From point [i,j], determine the indices
        ! of the best defined upstream point, as 
        ! determined from the velocity components ux and uy. 

        implicit none 

        integer, intent(OUT) :: i0 
        integer, intent(OUT) :: j0
        integer, intent(IN)  :: i 
        integer, intent(IN)  :: j
        integer, intent(IN)  :: ux(:,:) 
        integer, intent(IN)  :: uy(:,:) 

        ! Local variables
        integer :: im1, ip1, jm1, jp1 
        integer :: nx, ny 
        real(wp) :: ux_aa
        real(wp) :: uy_aa
        
        nx = size(ux,1) 
        ny = size(ux,2) 

        ! Define neighbor indices
        im1 = max(i-1,1)
        ip1 = min(i+1,nx)
        jm1 = max(j-1,1)
        jp1 = min(j+1,ny)
        
        ! Determine upstream node(s) 

        ux_aa = 0.5*(ux(i,j)+ux(im1,j))
        uy_aa = 0.5*(uy(i,j)+uy(i,jm1))
        
        if (ux_aa .ge. 0.0) then 
            i0 = im1
        else 
            i0 = ip1 
        end if 

        if (uy_aa .ge. 0.0) then 
            j0 = jm1
        else 
            j0 = jp1  
        end if 
        
        return 

    end subroutine find_upstream_neighbor

    function mean_mask(var,mask) result(ave)

        implicit none 

        real(wp), intent(IN) :: var(:,:) 
        logical,    intent(IN) :: mask(:,:) 
        real(wp) :: ave 
        integer :: n 

        n = count(mask)
        
        if (n .gt. 0) then 
            ave = sum(var,mask=mask) / real(n,wp)
        else 
            ave = 0.0 
        end if 

        return 

    end function mean_mask
    
    elemental subroutine minmax(var,var_lim)

        implicit none 

        real(wp), intent(INOUT) :: var 
        real(wp), intent(IN)    :: var_lim 

        if (var .lt. -var_lim) then 
            var = -var_lim 
        else if (var .gt. var_lim) then 
            var =  var_lim 
        end if 

        return 

    end subroutine minmax

    subroutine set_boundaries_2D_aa(var,boundaries,var_ref)

        implicit none 

        real(wp), intent(INOUT) :: var(:,:) 
        character(len=*), intent(IN) :: boundaries 
        real(wp), intent(IN), optional :: var_ref(:,:) 

        ! Local variables 
        integer :: nx, ny  

        nx = size(var,1) 
        ny = size(var,2) 

        select case(trim(boundaries))

            case("zeros","EISMINT")

                ! Set border values to zero
                var(1,:)  = 0.0
                var(nx,:) = 0.0

                var(:,1)  = 0.0
                var(:,ny) = 0.0

            case("periodic","periodic-xy") 

                var(1:2,:)     = var(nx-3:nx-2,:) 
                var(nx-1:nx,:) = var(2:3,:) 

                var(:,1:2)     = var(:,ny-3:ny-2) 
                var(:,ny-1:ny) = var(:,2:3) 
            
            case("periodic-x") 

                ! Periodic x 
                var(1:2,:)     = var(nx-3:nx-2,:) 
                var(nx-1:nx,:) = var(2:3,:) 
                
                ! Infinite (free-slip too)
                var(:,1)  = var(:,2)
                var(:,ny) = var(:,ny-1)

            case("MISMIP3D")

                ! === MISMIP3D =====
                var(1,:)    = var(2,:)          ! x=0, Symmetry 
                var(nx,:)   = 0.0               ! x=800km, no ice
                
!                var(:,1)    = var(:,2)          ! y=-50km, Free-slip condition
!                var(:,ny)   = var(:,ny-1)       ! y= 50km, Free-slip condition

            case("TROUGH")

                ! === MISMIP3D =====
                var(1,:)    = var(2,:)          ! x=0, Symmetry 
                var(nx,:)   = 0.0               ! x=800km, no ice
                
            case("infinite")
                ! Set border points equal to inner neighbors 

                call fill_borders_2D(var,nfill=1)

            case("fixed") 
                ! Set border points equal to prescribed values from array

                call fill_borders_2D(var,nfill=1,fill=var_ref)

            case DEFAULT 

                write(io_unit_err,*) "set_boundaries_2D_aa:: error: boundary method not recognized."
                write(io_unit_err,*) "boundaries = ", trim(boundaries)
                stop 

        end select 

        return

    end subroutine set_boundaries_2D_aa

    subroutine set_boundaries_3D_aa(var,boundaries,var_ref)

        implicit none 

        real(wp), intent(INOUT) :: var(:,:,:) 
        character(len=*), intent(IN) :: boundaries 
        real(wp), intent(IN), optional :: var_ref(:,:,:) 

        ! Local variables 
        integer :: nx, ny, nz  
        integer :: k 

        nx = size(var,1) 
        ny = size(var,2) 
        nz = size(var,3) 

        if (present(var_ref)) then

            do k = 1, nz 
                call set_boundaries_2D_aa(var(:,:,k),boundaries,var_ref(:,:,k))
            end do 
        
        else 
        
            do k = 1, nz 
                call set_boundaries_2D_aa(var(:,:,k),boundaries)
            end do 
        
        end if 
        
        return

    end subroutine set_boundaries_3D_aa

    subroutine set_boundaries_2D_acx(var_acx,boundaries)

        implicit none 

        real(wp), intent(INOUT) :: var_acx(:,:) 
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer :: nx, ny  

        nx = size(var_acx,1) 
        ny = size(var_acx,2) 

        select case(trim(boundaries))

            case("periodic") 

                var_acx(1,:)    = var_acx(nx-2,:) 
                var_acx(nx-1,:) = var_acx(2,:) 
                var_acx(nx,:)   = var_acx(3,:) 
                var_acx(:,1)    = var_acx(:,ny-1)
                var_acx(:,ny)   = var_acx(:,2) 
                
            case("periodic-x") 
                
                var_acx(1,:)    = var_acx(nx-2,:) 
                var_acx(nx-1,:) = var_acx(2,:) 
                var_acx(nx,:)   = var_acx(3,:) 
                var_acx(:,1)    = var_acx(:,2)
                var_acx(:,ny)   = var_acx(:,ny-1) 

            case("infinite") 
                
                var_acx(1,:)    = var_acx(2,:) 
                var_acx(nx-1,:) = var_acx(nx-2,:) 
                var_acx(nx,:)   = var_acx(nx-1,:) 
                var_acx(:,1)    = var_acx(:,2)
                var_acx(:,ny)   = var_acx(:,ny-1) 

            case("MISMIP3D","TROUGH") 
                
                !var_acx(1,:)    = var_acx(2,:) 
                var_acx(nx-1,:) = var_acx(nx-2,:) 
                var_acx(nx,:)   = var_acx(nx-1,:) 
                var_acx(:,1)    = var_acx(:,2)
                var_acx(:,ny)   = var_acx(:,ny-1) 

        end select 

        return 

    end subroutine set_boundaries_2D_acx

    subroutine set_boundaries_3D_acx(var_acx,boundaries)

        implicit none 

        real(wp), intent(INOUT) :: var_acx(:,:,:) 
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer :: nx, ny  

        nx = size(var_acx,1) 
        ny = size(var_acx,2) 

        select case(trim(boundaries))

            case("periodic") 

                var_acx(1,:,:)    = var_acx(nx-2,:,:) 
                var_acx(nx-1,:,:) = var_acx(2,:,:) 
                var_acx(nx,:,:)   = var_acx(3,:,:) 
                var_acx(:,1,:)    = var_acx(:,ny-1,:)
                var_acx(:,ny,:)   = var_acx(:,2,:) 
                
            case("periodic-x") 
                
                var_acx(1,:,:)    = var_acx(nx-2,:,:) 
                var_acx(nx-1,:,:) = var_acx(2,:,:) 
                var_acx(nx,:,:)   = var_acx(3,:,:) 
                var_acx(:,1,:)    = var_acx(:,2,:)
                var_acx(:,ny,:)   = var_acx(:,ny-1,:) 

            case("infinite") 
                
                var_acx(1,:,:)    = var_acx(2,:,:) 
                var_acx(nx-1,:,:) = var_acx(nx-2,:,:) 
                var_acx(nx,:,:)   = var_acx(nx-1,:,:) 
                var_acx(:,1,:)    = var_acx(:,2,:)
                var_acx(:,ny,:)   = var_acx(:,ny-1,:) 

            case("MISMIP3D","TROUGH") 
                
                var_acx(1,:,:)    = var_acx(2,:,:) 
                var_acx(nx-1,:,:) = var_acx(nx-2,:,:) 
                var_acx(nx,:,:)   = var_acx(nx-1,:,:) 
                var_acx(:,1,:)    = var_acx(:,2,:)
                var_acx(:,ny,:)   = var_acx(:,ny-1,:) 

        end select 

        return 

    end subroutine set_boundaries_3D_acx
    
    subroutine set_boundaries_2D_acy(var_acy,boundaries)

        implicit none 

        real(wp), intent(INOUT) :: var_acy(:,:) 
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer :: nx, ny  

        nx = size(var_acy,1) 
        ny = size(var_acy,2) 

        select case(trim(boundaries))

            case("periodic") 

                var_acy(1,:)    = var_acy(nx-1,:) 
                var_acy(nx,:)   = var_acy(2,:) 
                var_acy(:,1)    = var_acy(:,ny-2)
                var_acy(:,ny-1) = var_acy(:,2) 
                var_acy(:,ny)   = var_acy(:,3)

            case("periodic-x") 
                
                var_acy(1,:)    = var_acy(nx-1,:) 
                var_acy(nx,:)   = var_acy(2,:) 
                var_acy(:,1)    = var_acy(:,2)
                var_acy(:,ny-1) = var_acy(:,ny-2) 
                var_acy(:,ny)   = var_acy(:,ny-1)

            case("infinite") 
                
                var_acy(1,:)    = var_acy(2,:) 
                var_acy(nx,:)   = var_acy(nx-1,:) 
                var_acy(:,1)    = var_acy(:,2)
                var_acy(:,ny-1) = var_acy(:,ny-2) 
                var_acy(:,ny)   = var_acy(:,ny-1)

            case("MISMIP3D","TROUGH") 
                
                var_acy(1,:)    = var_acy(2,:) 
                var_acy(nx,:)   = var_acy(nx-1,:) 
                var_acy(:,1)    = var_acy(:,2)
                var_acy(:,ny-1) = var_acy(:,ny-2) 
                var_acy(:,ny)   = var_acy(:,ny-1)

        end select 

        return 

    end subroutine set_boundaries_2D_acy

    subroutine set_boundaries_3D_acy(var_acy,boundaries)

        implicit none 

        real(wp), intent(INOUT) :: var_acy(:,:,:) 
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer :: nx, ny  

        nx = size(var_acy,1) 
        ny = size(var_acy,2) 

        select case(trim(boundaries))

            case("periodic") 

                var_acy(1,:,:)    = var_acy(nx-1,:,:) 
                var_acy(nx,:,:)   = var_acy(2,:,:) 
                var_acy(:,1,:)    = var_acy(:,ny-2,:)
                var_acy(:,ny-1,:) = var_acy(:,2,:) 
                var_acy(:,ny,:)   = var_acy(:,3,:)

            case("periodic-x") 
                
                var_acy(1,:,:)    = var_acy(nx-1,:,:) 
                var_acy(nx,:,:)   = var_acy(2,:,:) 
                var_acy(:,1,:)    = var_acy(:,2,:)
                var_acy(:,ny-1,:) = var_acy(:,ny-2,:) 
                var_acy(:,ny,:)   = var_acy(:,ny-1,:)

            case("infinite") 
                
                var_acy(1,:,:)    = var_acy(2,:,:) 
                var_acy(nx,:,:)   = var_acy(nx-1,:,:) 
                var_acy(:,1,:)    = var_acy(:,2,:)
                var_acy(:,ny-1,:) = var_acy(:,ny-2,:) 
                var_acy(:,ny,:)   = var_acy(:,ny-1,:)

            case("MISMIP3D","TROUGH") 
                
                var_acy(1,:,:)    = var_acy(2,:,:) 
                var_acy(nx,:,:)   = var_acy(nx-1,:,:) 
                var_acy(:,1,:)    = var_acy(:,2,:)
                var_acy(:,ny-1,:) = var_acy(:,ny-2,:) 
                var_acy(:,ny,:)   = var_acy(:,ny-1,:)

        end select 

        return 

    end subroutine set_boundaries_3D_acy

    subroutine fill_borders_2D(var,nfill,fill)

        implicit none 

        real(wp), intent(INOUT) :: var(:,:) 
        integer,    intent(IN)    :: nfill        ! How many neighbors to fill in 
        real(wp), intent(IN), optional :: fill(:,:) ! Values to impose 

        ! Local variables 
        integer :: i, j, nx, ny, q 
        
        nx = size(var,1)
        ny = size(var,2)

        if (present(fill)) then 
            ! Fill with prescribed values from array 'fill' 

            do q = 1, nfill 
                var(q,:)      = fill(nfill+1,:)      
                var(nx-q+1,:) = fill(nx-nfill,:)   
                
                var(:,q)      = fill(:,nfill+1)     
                var(:,ny-q+1) = fill(:,ny-nfill)  
            end do 

        else 
            ! Fill with interior neighbor values 

            do q = 1, nfill 
                var(q,:)      = var(nfill+1,:)      
                var(nx-q+1,:) = var(nx-nfill,:)   
                
                var(:,q)      = var(:,nfill+1)     
                var(:,ny-q+1) = var(:,ny-nfill)  
            end do 

        end if 

        return 

    end subroutine fill_borders_2D

    subroutine fill_borders_3D(var,nfill)
        ! 3rd dimension is not filled (should be vertical dimension)

        implicit none 

        real(wp), intent(INOUT) :: var(:,:,:) 
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

    subroutine smooth_gauss_3D(var,dx,f_sigma,mask_apply,mask_use)

        ! Smooth out strain heating to avoid noise 

        implicit none

        real(wp),   intent(INOUT) :: var(:,:,:)      ! nx,ny,nz_aa: 3D variable
        real(wp),   intent(IN)    :: dx 
        real(wp),   intent(IN)    :: f_sigma  
        logical,    intent(IN), optional :: mask_apply(:,:) 
        logical,    intent(IN), optional :: mask_use(:,:) 

        ! Local variables
        integer    :: k, nz_aa  

        nz_aa = size(var,3)

        do k = 1, nz_aa 
             call smooth_gauss_2D(var(:,:,k),dx,f_sigma,mask_apply,mask_use)
        end do 

        return 

    end subroutine smooth_gauss_3D
    
    subroutine smooth_gauss_2D(var,dx,f_sigma,mask_apply,mask_use)
        ! Smooth out a field to avoid noise 
        ! mask_apply designates where smoothing should be applied 
        ! mask_use   designates which points can be considered in the smoothing filter 

        implicit none

        real(wp),   intent(INOUT) :: var(:,:)      ! [nx,ny] 2D variable
        real(wp),   intent(IN)    :: dx 
        real(wp),   intent(IN)    :: f_sigma  
        logical,    intent(IN), optional :: mask_apply(:,:) 
        logical,    intent(IN), optional :: mask_use(:,:) 

        ! Local variables
        integer  :: i, j, nx, ny, n, n2
        real(wp) :: sigma    
        real(wp), allocatable :: filter0(:,:), filter(:,:) 
        real(wp), allocatable :: var_old(:,:) 
        logical,  allocatable :: mask_apply_local(:,:) 
        logical,  allocatable :: mask_use_local(:,:)

        nx    = size(var,1)
        ny    = size(var,2)

        ! Safety check
        if (f_sigma .lt. 1.0_wp) then 
            write(io_unit_err,*) ""
            write(io_unit_err,*) "smooth_gauss_2D:: Error: f_sigma must be >= 1."
            write(io_unit_err,*) "f_sigma: ", f_sigma 
            write(io_unit_err,*) "dx:      ", dx 
            stop 
        end if 

        ! Get smoothing radius as standard devation of Gaussian function
        sigma = dx*f_sigma 

        ! Determine half-width of filter as 3-sigma
        n2 = 3*ceiling(f_sigma)

        ! Get total number of points for filter window in each direction
        n = 2*n2+1
        
        allocate(var_old(nx+2*n2,ny+2*n2))
        allocate(mask_apply_local(nx+2*n2,ny+2*n2))
        allocate(mask_use_local(nx+2*n2,ny+2*n2))
        allocate(filter0(n,n))
        allocate(filter(n,n))

        ! Check whether mask_apply is available 
        if (present(mask_apply)) then 
            ! use mask_use to define neighborhood points
            
            mask_apply_local = .FALSE. 
            mask_apply_local(n2+1:n2+nx,n2+1:n2+ny) = mask_apply 

        else
            ! Assume that everywhere should be smoothed

            mask_apply_local = .FALSE. 
            mask_apply_local(n2+1:n2+nx,n2+1:n2+ny) = .TRUE.
        
        end if

        ! Check whether mask_use is available 
        if (present(mask_use)) then 
            ! use mask_use to define neighborhood points
            
            mask_use_local = .TRUE.
            mask_use_local(n2+1:n2+nx,n2+1:n2+ny) = mask_use 

        else
            ! Assume that mask_apply also gives the points to use for smoothing 

            mask_use_local = mask_apply_local
        
        end if

        ! Calculate default 2D Gaussian smoothing kernel
        filter0 = gauss_values(dx,dx,sigma=sigma,n=n)

        var_old = 0.0 
        var_old(n2+1:n2+nx,n2+1:n2+ny) = var 
        var_old(1:n2,n2+1:n2+ny)       = var(n2:1:-1,:)
        var_old(nx+1:nx+n2,n2+1:n2+ny) = var((nx-n2+1):nx,:)
        var_old(n2+1:n2+nx,1:n2)       = var(:,n2:1:-1)
        var_old(n2+1:n2+nx,ny+1:ny+n2) = var(:,(ny-n2+1):ny)
        var_old(1:n2,n2+1:n2+ny)       = var(n2:1:-1,:)
        var_old(nx+1:nx+n2,n2+1:n2+ny) = var((nx-n2+1):nx,:)
        var_old(n2+1:n2+nx,1:n2)       = var(:,n2:1:-1)
        var_old(n2+1:n2+nx,ny+1:ny+n2) = var(:,(ny-n2+1):ny)
        
        !!$omp parallel do collapse(2) private(i,j,filter)
        do j = n2+1, n2+ny 
        do i = n2+1, n2+nx 

            if (mask_apply_local(i,j)) then 
                ! Apply smoothing to this point 

                ! Limit filter input to neighbors of interest
                filter = filter0 
                where(.not. mask_use_local(i-n2:i+n2,j-n2:j+n2)) filter = 0.0

                ! If neighbors are available, normalize and perform smoothing  
                if (sum(filter) .gt. 0.0) then 
                    filter = filter/sum(filter)
                    var(i-n2,j-n2) = sum(var_old(i-n2:i+n2,j-n2:j+n2)*filter) 
                end if  

            end if 

        end do 
        end do 
        !!$omp end parallel do

        return 

    end subroutine smooth_gauss_2D

    function gauss_values(dx,dy,sigma,n) result(filt)
        ! Calculate 2D Gaussian smoothing kernel
        ! https://en.wikipedia.org/wiki/Gaussian_blur

        implicit none 

        real(wp), intent(IN) :: dx 
        real(wp), intent(IN) :: dy 
        real(wp), intent(IN) :: sigma 
        integer,    intent(IN) :: n 
        real(wp) :: filt(n,n) 

        ! Local variables 
        real(wp) :: x, y  
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

    subroutine adjust_topography_gradients(z_bed,H_ice,grad_lim,dx,boundaries)
        ! Smooth the bedrock topography and corresponding ice thickness,
        ! so that specified limit on gradients is not exceeded. Only
        ! apply smoothing directly in places where gradient is too large.

        implicit none 

        real(wp), intent(INOUT) :: z_bed(:,:) 
        real(wp), intent(INOUT) :: H_ice(:,:) 
        real(wp), intent(IN)    :: grad_lim 
        real(wp), intent(IN)    :: dx 
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer  :: i, j, q, nx, ny 
        integer  :: im1, ip1, jm1, jp1 
        real(wp) :: dy
        real(wp), allocatable :: dzbdx(:,:)
        real(wp), allocatable :: dzbdy(:,:)
        real(wp), allocatable :: f_ice(:,:)
        logical,  allocatable :: mask_apply(:,:) 
        logical,  allocatable :: mask_use(:,:) 

        integer, parameter :: iter_max = 50 

        nx = size(z_bed,1)
        ny = size(z_bed,2) 

        dy = dx 

        allocate(dzbdx(nx,ny))
        allocate(dzbdy(nx,ny))
        allocate(f_ice(nx,ny))
        allocate(mask_apply(nx,ny))
        allocate(mask_use(nx,ny))

        ! Smooth z_bed in specific locations if gradients are exceeded.

        mask_use = .TRUE. 
        f_ice    = 1.0 

        do q = 1, iter_max

            ! Calculate bedrock gradients (f_ice and grad_lim are not used)
            call calc_gradient_acx(dzbdx,z_bed,f_ice,dx,grad_lim=100.0_wp, &
                                        margin2nd=.FALSE.,zero_outside=.FALSE.,boundaries=boundaries)
            call calc_gradient_acy(dzbdy,z_bed,f_ice,dy,grad_lim=100.0_wp, &
                                        margin2nd=.FALSE.,zero_outside=.FALSE.,boundaries=boundaries)

            ! Determine where gradients are too large
            mask_apply = .FALSE.
            do j = 3, ny-3
            do i = 3, nx-3
                
                ! Get neighbor indices
                call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

                if (abs(dzbdx(i,j)) .ge. grad_lim) then 
                    mask_apply(i,j)   = .TRUE. 
                    mask_apply(ip1,j) = .TRUE. 
                end if

                if (abs(dzbdy(i,j)) .ge. grad_lim) then 
                    mask_apply(i,j)   = .TRUE. 
                    mask_apply(i,jp1) = .TRUE. 
                end if
            end do 
            end do 

            write(*,*) "z_bed smoothing: ", q, count(mask_apply),  &
                                        maxval(abs(dzbdx(3:nx-3,3:ny-3))), &
                                        maxval(abs(dzbdy(3:nx-3,3:ny-3)))

            if (count(mask_apply) .eq. 0) exit 

            ! Smooth z_bed at desired locations, and H_ice so that H_ice avoids spurious patterns
            call smooth_gauss_2D(z_bed,dx=dx,f_sigma=2.0_wp,mask_apply=mask_apply,mask_use=mask_use)
            call smooth_gauss_2D(H_ice,dx=dx,f_sigma=2.0_wp,mask_apply=mask_apply,mask_use=mask_use)
            
        end do

        return

    end subroutine adjust_topography_gradients

    subroutine regularize2D_gauss(var,H_ice,dx)
        ! Ensure smoothness in 2D fields (ie, no checkerboard patterns)
        ! ajr: doesnt work!

        implicit none 

        real(wp), intent(INOUT) :: var(:,:)       ! aa-nodes
        real(wp), intent(IN)    :: H_ice(:,:)     ! aa-nodes
        real(wp), intent(IN)    :: dx 

        ! Local variables
        integer    :: i, j, nx, ny  
        integer    :: im1, ip1, jm1, jp1 
        real(wp) :: varx(2), vary(2)
        logical    :: check_x, check_y 
        
        logical, allocatable :: bad_pts(:,:) 

        nx = size(var,1)
        ny = size(var,2) 

        allocate(bad_pts(nx,ny)) 

        ! All points are good initially 
        bad_pts = .FALSE. 

        do j = 1, ny 
        do i = 1, nx

            if (H_ice(i,j) .gt. 0.0) then 
                ! Only check ice-covered points 

                im1 = max(1, i-1)
                ip1 = min(nx,i+1)
                
                jm1 = max(1, j-1)
                jp1 = min(ny,j+1)

                varx = [var(im1,j),var(ip1,j)]
                where([H_ice(im1,j),H_ice(ip1,j)] .eq. 0.0_wp) varx = missing_value 

                vary = [var(i,jm1),var(i,jp1)]
                where([H_ice(i,jm1),H_ice(i,jp1)] .eq. 0.0_wp) vary = missing_value 
                
                ! Check if checkerboard exists in each direction 
                check_x = (count(varx .gt. var(i,j) .and. varx.ne.missing_value) .eq. 2 .or. &
                           count(varx .lt. var(i,j) .and. varx.ne.missing_value) .eq. 2) 

                check_y = (count(vary .gt. var(i,j) .and. vary.ne.missing_value) .eq. 2 .or. &
                           count(vary .lt. var(i,j) .and. vary.ne.missing_value) .eq. 2) 
                
                ! If check is true, mark point for later treatment 
                if (check_x .or. check_y) bad_pts(i,j) = .TRUE.  

            end if 

        end do 
        end do 

        ! Now apply Gaussian smoothing to bad points, with a wide radius 
        call smooth_gauss_2D(var,mask_apply=bad_pts,dx=dx,f_sigma=5.0_wp, &
                                mask_use=H_ice.gt.0.0_wp .and. (.not. bad_pts))

        return 

    end subroutine regularize2D_gauss

    subroutine regularize2D(var,H_ice,dx)
        ! Ensure smoothness in 2D fields (ie, no checkerboard patterns)

        implicit none 

        real(wp), intent(INOUT) :: var(:,:)       ! aa-nodes
        real(wp), intent(IN)    :: H_ice(:,:)     ! aa-nodes
        real(wp), intent(IN)    :: dx 

        ! Local variables
        integer    :: i, j, nx, ny, n   
        integer    :: im1, ip1, jm1, jp1 
        real(wp), allocatable :: var0(:,:) 
        real(wp) :: varx(2), vary(2), var9(3,3)
        logical    :: check_x, check_y 
        
        integer    :: q, qmax, npts

        qmax = 10

        nx = size(var,1)
        ny = size(var,2) 

        allocate(var0(nx,ny))

        do q = 1, qmax

            var0 = var 
            npts = 0

        do j = 2, ny-1 
        do i = 2, nx-1

            if (H_ice(i,j) .gt. 0.0) then 
                ! Only apply to ice-covered points 

                im1 = max(1, i-1)
                ip1 = min(nx,i+1)
                
                jm1 = max(1, j-1)
                jp1 = min(ny,j+1)

                varx = [var0(im1,j),var0(ip1,j)]
                where([H_ice(im1,j),H_ice(ip1,j)] .eq. 0.0_wp) varx = missing_value 

                vary = [var0(i,jm1),var0(i,jp1)]
                where([H_ice(i,jm1),H_ice(i,jp1)] .eq. 0.0_wp) vary = missing_value 
                
                ! Check if checkerboard exists in each direction 
                check_x = (count(varx .gt. var0(i,j) .and. varx.ne.missing_value) .eq. 2 .or. &
                           count(varx .lt. var0(i,j) .and. varx.ne.missing_value) .eq. 2) 

                check_y = (count(vary .gt. var0(i,j) .and. vary.ne.missing_value) .eq. 2 .or. &
                           count(vary .lt. var0(i,j) .and. vary.ne.missing_value) .eq. 2) 
                
                if (check_x .or. check_y) then 
                    ! Checkerboard exists, apply 9-point neighborhood average

                    var9 = var0(i-1:i+1,j-1:j+1)
                    where(H_ice(i-1:i+1,j-1:j+1) .eq. 0.0_wp) var9 = missing_value 

                    n = count(var9 .ne. missing_value) 

                    var(i,j) = sum(var9,mask=var9.ne.missing_value) / real(n,wp)
                    npts     = npts + 1
                end if 

            end if 

        end do 
        end do 

        if (npts .eq. 0) exit 

        end do 

        return 

    end subroutine regularize2D

    ! === Generic integration functions ============

    function calc_vertical_integrated_3D(var,zeta) result(var_int)
        ! Vertically integrate a field 3D field (nx,ny,nz)
        ! layer by layer (in the z-direction), return a 3D array

        implicit none

        real(wp), intent(IN) :: var(:,:,:)
        real(wp), intent(IN) :: zeta(:)
        real(wp) :: var_int(size(var,1),size(var,2),size(var,3))

        ! Local variables 
        integer :: i, j, nx, ny

        nx = size(var,1)
        ny = size(var,2)

        !!$omp parallel do collapse(2) private(i,j)
        do j = 1, ny
        do i = 1, nx
            var_int(i,j,:) = integrate_trapezoid1D_1D(var(i,j,:),zeta)
        end do
        end do
        !!$omp end parallel do

        return

    end function calc_vertical_integrated_3D

    function calc_vertical_integrated_2D(var,zeta) result(var_int)
        ! Vertically integrate a field 3D field (nx,ny,nz) 
        ! to the surface, return a 2D array (nx,ny)
        
        implicit none

        real(wp), intent(IN) :: var(:,:,:)
        real(wp), intent(IN) :: zeta(:)
        real(wp) :: var_int(size(var,1),size(var,2))

        ! Local variables 
        integer :: i, j, nx, ny

        nx = size(var,1)
        ny = size(var,2)

        !!$omp parallel do collapse(2) private(i,j)
        do j = 1, ny
        do i = 1, nx
            var_int(i,j) = integrate_trapezoid1D_pt(var(i,j,:),zeta)
        end do
        end do
        !!$omp end parallel do 

        return

    end function calc_vertical_integrated_2D
    
    function integrate_trapezoid1D_pt(var,zeta) result(var_int)
        ! Integrate a variable from the base to height zeta(nk) in the ice column.
        ! The value of the integral using the trapezium rule can be found using
        ! integral = (b - a)*((f(a) +f(b))/2 + Σ_1_n-1(f(k)) )/n 
        ! Returns a point of integrated value of var at level zeta(nk).

        implicit none

        real(wp), intent(IN) :: var(:)
        real(wp), intent(IN) :: zeta(:)
        real(wp) :: var_int

        ! Local variables 
        integer :: k, nk
        real(wp) :: var_mid 
        
        nk = size(var,1)

        ! Initial value is zero
        var_int = 0.0_wp 

        ! Intermediate values include sum of all previous values 
        ! Take current value as average between points
        do k = 2, nk
            var_mid = 0.5_wp*(var(k)+var(k-1))
            if (abs(var_mid) .lt. TOL_UNDERFLOW) var_mid = 0.0_wp 
            var_int = var_int + var_mid*(zeta(k) - zeta(k-1))
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

        real(wp), intent(IN) :: var(:)
        real(wp), intent(IN) :: zeta(:)
        real(wp) :: var_int(size(var,1))

        ! Local variables 
        integer    :: k, nk
        real(wp) :: var_mid 

        nk = size(var,1)

        ! Initial value is zero
        var_int(1:nk) = 0.0_wp 

        ! Intermediate values include sum of all previous values 
        ! Take current value as average between points
        do k = 2, nk
            var_mid = 0.5_wp*(var(k)+var(k-1))
            if (abs(var_mid) .lt. TOL_UNDERFLOW) var_mid = 0.0_wp 
            var_int(k:nk) = var_int(k:nk) + 0.5_wp*(var(k)+var(k-1))*(zeta(k) - zeta(k-1))
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

        real(wp) :: x(:)
        real(wp) :: y(:)
        real(wp) :: result

        integer :: ntab

        real(wp) :: del(3)
        real(wp) :: e
        real(wp) :: f
        real(wp) :: feints
        real(wp) :: g(3)
        integer    :: i
        integer    :: n
        real(wp) :: pi(3)
        real(wp) :: sum1

        real(wp) :: x1
        real(wp) :: x2
        real(wp) :: x3

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
        real(wp), allocatable :: zeta0(:) 
        real(wp), allocatable :: zeta(:)
        real(wp), allocatable :: var0(:) 
        real(wp), allocatable :: var(:) 
        real(wp), allocatable :: var_ints(:)
        real(wp) :: var_int 
        real(wp) :: var_int_00

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
