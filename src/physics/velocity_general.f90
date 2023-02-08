module velocity_general 
    ! This module contains general routines that are used by several solvers. 
    
    use yelmo_defs ,only  : sp, dp, wp, tol_underflow, io_unit_err, rho_ice, rho_sw, rho_w, g, &
                            jacobian_3D_class
    use yelmo_tools, only : get_neighbor_indices, stagger_aa_ab, stagger_aa_ab_ice, &
                    acx_to_nodes, acy_to_nodes, acx_to_nodes_3D, acy_to_nodes_3D, &
                    integrate_trapezoid1D_1D, integrate_trapezoid1D_pt, minmax

    use deformation, only : calc_strain_rate_horizontal_2D

    use solver_ssa_ac, only : ssa_diagnostics_write_init, ssa_diagnostics_write_step

    implicit none 

    private 
    public :: calc_uz_3D_jac
    public :: calc_uz_3D
    !public :: calc_uz_3D_aa
    public :: calc_driving_stress
    public :: calc_driving_stress_gl
    public :: calc_lateral_bc_stress_2D
    public :: set_inactive_margins
    public :: calc_ice_flux
    public :: calc_vel_ratio
    public :: limit_vel

    public :: picard_calc_error 
    public :: picard_calc_error_angle 
    public :: picard_calc_convergence_l1rel_matrix
    public :: picard_calc_convergence_l2
    public :: picard_relax 
    
contains 
    
    subroutine calc_uz_3D_jac(uz,uz_star,ux,uy,jvel,H_ice,f_ice,f_grnd,smb,bmb,dHdt,dzsdt, &
                                    dzsdx,dzsdy,dzbdx,dzbdy,zeta_aa,zeta_ac,dx,dy,use_bmb,boundaries)
        ! Following algorithm outlined by the Glimmer ice sheet model:
        ! https://www.geos.ed.ac.uk/~mhagdorn/glide/glide-doc/glimmer_htmlse9.html#x17-660003.1.5

        ! Note: rate of ice-base elevation change (dzbdt) is deduced from dzbdt = dzsdt - dhdt. 
        ! This formulation does not depend on rate of bedrock uplift (which is implicit in dzsdt),
        ! and is valid for both grounded and floating ice. 

        implicit none 

        real(wp), intent(OUT) :: uz(:,:,:)          ! nx,ny,nz_ac
        real(wp), intent(OUT) :: uz_star(:,:,:)     ! nx,ny,nz_ac
        real(wp), intent(IN)  :: ux(:,:,:)          ! nx,ny,nz_aa
        real(wp), intent(IN)  :: uy(:,:,:)          ! nx,ny,nz_aa
        type(jacobian_3D_class), intent(IN) :: jvel 
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: f_grnd(:,:)
        real(wp), intent(IN)  :: smb(:,:) 
        real(wp), intent(IN)  :: bmb(:,:) 
        real(wp), intent(IN)  :: dHdt(:,:) 
        real(wp), intent(IN)  :: dzsdt(:,:) 
        real(wp), intent(IN)  :: dzsdx(:,:) 
        real(wp), intent(IN)  :: dzsdy(:,:) 
        real(wp), intent(IN)  :: dzbdx(:,:) 
        real(wp), intent(IN)  :: dzbdy(:,:) 
        real(wp), intent(IN)  :: zeta_aa(:)    ! z-coordinate, aa-nodes 
        real(wp), intent(IN)  :: zeta_ac(:)    ! z-coordinate, ac-nodes  
        real(wp), intent(IN)  :: dx 
        real(wp), intent(IN)  :: dy
        logical,  intent(IN)  :: use_bmb
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer :: i, j, k, nx, ny, nz_aa, nz_ac
        integer :: im1, ip1, jm1, jp1  
        integer  :: im1m, ip1m, jm1m, jp1m 
        integer :: kup, kdn, kmid
        real(wp) :: f_bmb
        real(wp) :: H_now
        real(wp) :: H_inv
        real(wp) :: dzbdx_aa
        real(wp) :: dzbdy_aa
        real(wp) :: dzsdx_aa
        real(wp) :: dzsdy_aa 
        real(wp) :: dudx_aa
        real(wp) :: dvdy_aa
        real(wp) :: dudz_aa 
        real(wp) :: dvdz_aa 
        real(wp) :: ux_aa 
        real(wp) :: uy_aa 
        real(wp) :: uz_grid 
        real(wp) :: uz_srf 
        real(wp) :: zeta_now 
        real(wp) :: c_x 
        real(wp) :: c_y 
        real(wp) :: c_t 

        real(wp) :: dzsdtn(4)
        real(wp) :: dhdtn(4)
        real(wp) :: dzbdxn(4)
        real(wp) :: dzbdyn(4)
        real(wp) :: dzsdxn(4)
        real(wp) :: dzsdyn(4)
        real(wp) :: dudxn(4) 
        real(wp) :: dvdyn(4) 
        real(wp) :: uxn_up(4) 
        real(wp) :: uxn_dn(4) 
        real(wp) :: uxn(4) 
        real(wp) :: uyn_up(4) 
        real(wp) :: uyn_dn(4)
        real(wp) :: uyn(4)  
        real(wp) :: dudzn(4) 
        real(wp) :: dvdzn(4) 
        
        real(wp) :: wt0
        real(wp) :: xn(4) 
        real(wp) :: yn(4) 
        real(wp) :: wtn(4)
        real(wp) :: wt2D 

        real(wp) :: dudxn8(8) 
        real(wp) :: dvdyn8(8) 
        real(wp) :: zn 
        real(wp) :: wtn8(8)
        real(wp) :: wt3D 

        real(wp) :: dzsdt_now
        real(wp) :: dhdt_now
        real(wp) :: dzbdt_now 

        real(wp), parameter :: uz_min = -10.0     ! [m/yr] Minimum allowed vertical velocity downwards for stability
        
        nx    = size(ux,1)
        ny    = size(ux,2)
        nz_aa = size(zeta_aa,1)
        nz_ac = size(zeta_ac,1) 
        
        ! Get nodes and weighting 
        wt0  = 1.0/sqrt(3.0)
        xn   = [wt0,-wt0,-wt0, wt0]
        yn   = [wt0, wt0,-wt0,-wt0]
        wtn  = [1.0,1.0,1.0,1.0]
        wt2D = 4.0   ! Surface area of square [-1:1,-1:1]=> 2x2 => 4 

        ! Get nodes and weighting 
        zn   = wt0
        wtn8 = [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0]
        wt3D = 8.0   ! Volume of square [-1:1,-1:1, -1:1]=> 2x2x2 => 8

        ! Initialize vertical velocity to zero 
        uz = 0.0 

        ! Define switch for bmb
        if (use_bmb) then 
            f_bmb = 1.0 
        else 
            f_bmb = 0.0 
        end if 

        ! Next, calculate vertical velocity at each point through the column

        !$omp parallel do 
        do j = 1, ny
        do i = 1, nx

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

            ! Get neighbor indices limited to ice-covered points
            im1m = im1
            if (f_ice(im1,j) .lt. 1.0) im1m = i  
            ip1m = ip1
            if (f_ice(ip1,j) .lt. 1.0) ip1m = i  
            jm1m = jm1 
            if (f_ice(i,jm1) .lt. 1.0) jm1m = j 
            jp1m = jp1 
            if (f_ice(i,jp1) .lt. 1.0) jp1m = j

            ! Diagnose rate of ice-base elevation change (needed for all points)
            dzsdt_now = dzsdt(i,j) 
            dhdt_now  = dhdt(i,j) 
            dzbdt_now = dzsdt_now - dhdt_now

            if (f_ice(i,j) .eq. 1.0) then

                H_now  = H_ice(i,j) 
                H_inv = 1.0/H_now 

                ! Get the centered ice-base gradient
                call acx_to_nodes(dzbdxn,dzbdx,i,j,xn,yn,im1m,ip1m,jm1m,jp1m)
                dzbdx_aa = sum(dzbdxn*wtn)/wt2D
                
                call acy_to_nodes(dzbdyn,dzbdy,i,j,xn,yn,im1m,ip1m,jm1m,jp1m)
                dzbdy_aa = sum(dzbdyn*wtn)/wt2D
                
                ! Get the centered surface gradient
                call acx_to_nodes(dzsdxn,dzsdx,i,j,xn,yn,im1m,ip1m,jm1m,jp1m)
                dzsdx_aa = sum(dzsdxn*wtn)/wt2D
                
                call acy_to_nodes(dzsdyn,dzsdy,i,j,xn,yn,im1m,ip1m,jm1m,jp1m)
                dzsdy_aa = sum(dzsdyn*wtn)/wt2D
                
                ! Get the aa-node centered horizontal velocity at the base
                call acx_to_nodes(uxn,ux(:,:,1),i,j,xn,yn,im1m,ip1m,jm1m,jp1m)
                ux_aa = sum(uxn*wtn)/wt2D
                
                call acy_to_nodes(uyn,uy(:,:,1),i,j,xn,yn,im1m,ip1m,jm1m,jp1m)
                uy_aa = sum(uyn*wtn)/wt2D
                
                ! Determine grid vertical velocity at the base due to sigma-coordinates 
                ! Glimmer, Eq. 3.35 
                ! ajr, 2020-01-27, untested:::
!                 uz_grid = dzsdt(i,j) + (ux_aa*dzsdx_aa + uy_aa*dzsdy_aa) &
!                             - ( (1.0_wp-zeta_ac(1))*dHdt(i,j) + ux_aa*dHdx_aa + uy_aa*dHdy_aa )
                uz_grid = 0.0_wp 

                ! ===================================================================
                ! Greve and Blatter (2009) style:

                ! Determine basal vertical velocity for this grid point 
                ! Following Eq. 5.31 of Greve and Blatter (2009)
                uz(i,j,1) = dzbdt_now + uz_grid + f_bmb*bmb(i,j) + ux_aa*dzbdx_aa + uy_aa*dzbdy_aa
                if (abs(uz(i,j,1)) .lt. TOL_UNDERFLOW) uz(i,j,1) = 0.0_wp 
                
                ! Set stability limit on basal uz value.
                ! This only gets applied in rare cases when something
                ! is going wrong in the model. 
                if (uz(i,j,1) .lt. uz_min) uz(i,j,1) = uz_min 

                ! Determine surface vertical velocity following kinematic boundary condition 
                ! Glimmer, Eq. 3.10 [or Folwer, Chpt 10, Eq. 10.8]
                !uz_srf = dzsdt(i,j) + ux_aa*dzsdx_aa + uy_aa*dzsdy_aa - smb(i,j) 
                
                ! Integrate upward to each point above base until just below surface is reached 
                ! Integrate on vertical ac-nodes (ie, vertical cell borders between aa-node centers)
                do k = 2, nz_ac 

                    ! Note: center of cell below this one is zeta_aa(k-1)
                    kmid = k-1

                    ! Get dudz/dvdz values at vertical aa-nodes, in order
                    ! to vertically integrate each cell up to ac-node border.

if (.FALSE.) then
    ! 2D QUADRATURE
                    call acx_to_nodes(dudxn,jvel%dxx(:,:,kmid),i,j,xn,yn,im1m,ip1m,jm1m,jp1m)
                    dudx_aa = sum(dudxn*wtn)/wt2D

                    call acy_to_nodes(dvdyn,jvel%dyy(:,:,kmid),i,j,xn,yn,im1m,ip1m,jm1m,jp1m)
                    dvdy_aa = sum(dvdyn*wtn)/wt2D

else 
    ! 3D QUADRATURE
                    call acx_to_nodes_3D(dudxn8,jvel%dxx,i,j,kmid,xn,yn,zn,im1m,ip1m,jm1m,jp1m)
                    dudx_aa = sum(dudxn8*wtn8)/wt3D

                    call acy_to_nodes_3D(dvdyn8,jvel%dyy,i,j,kmid,xn,yn,zn,im1m,ip1m,jm1m,jp1m)
                    dvdy_aa = sum(dvdyn8*wtn8)/wt3D

end if 

                    ! Calculate vertical velocity of this layer
                    ! (Greve and Blatter, 2009, Eq. 5.95)
                    uz(i,j,k) = uz(i,j,k-1) - H_now*(zeta_ac(k)-zeta_ac(k-1))*(dudx_aa+dvdy_aa)

                    ! Apply correction to match kinematic boundary condition at surface 
                    !uz(i,j,k) = uz(i,j,k) - zeta_ac(k)*(uz(i,j,k)-uz_srf)

                    if (abs(uz(i,j,k)) .lt. TOL_UNDERFLOW) uz(i,j,k) = 0.0_wp 
                    
                end do 


                ! === Also calculate adjusted vertical velocity to be used for temperature advection
                
                do k = 1, nz_ac 

                    ! Get the centered horizontal velocity of box on vertical ac-nodes at the level k 
                    ! ajr: Given that the correction is applied to uz, which is defined on 
                    ! ac-nodes, it seems the correction should also be calculated on ac-nodes.
                    ! Note: nz_ac = nz_aa + 1
    
                    if (k .eq. 1) then 
                        kup = k 
                        kdn = k 
                    else if (k .eq. nz_ac) then 
                        kup = k-1 
                        kdn = k-1 
                    else
                        kup = k 
                        kdn = k-1 
                    end if
                    
                    call acx_to_nodes(uxn_up,ux(:,:,kup),i,j,xn,yn,im1m,ip1m,jm1m,jp1m)
                    call acx_to_nodes(uxn_dn,ux(:,:,kdn),i,j,xn,yn,im1m,ip1m,jm1m,jp1m)
                    uxn = 0.5_wp*(uxn_up+uxn_dn)
                    ux_aa = sum(uxn*wtn)/wt2D
                    
                    call acy_to_nodes(uyn_up,uy(:,:,kup),i,j,xn,yn,im1m,ip1m,jm1m,jp1m)
                    call acy_to_nodes(uyn_dn,uy(:,:,kdn),i,j,xn,yn,im1m,ip1m,jm1m,jp1m)
                    uyn = 0.5_wp*(uyn_up+uyn_dn)
                    uy_aa = sum(uyn*wtn)/wt2D
                    
                    ! Take zeta directly at vertical cell edge where uz is calculated
                    ! (this is also where ux_aa and uy_aa are calculated above)
                    zeta_now = zeta_ac(k)

                    ! Calculate sigma-coordinate derivative correction factors
                    ! (Greve and Blatter, 2009, Eqs. 5.131 and 5.132, 
                    !  also shown in 1D with Eq. 5.145)

                    ! Note: not dividing by H here, since this is done in the thermodynamics advection step
                    c_x = -( (1.0-zeta_now)*dzbdx_aa  + zeta_now*dzsdx_aa )
                    c_y = -( (1.0-zeta_now)*dzbdy_aa  + zeta_now*dzsdy_aa )
                    c_t = -( (1.0-zeta_now)*dzbdt_now + zeta_now*dzsdt_now )
                    
                    ! Calculate adjusted vertical velocity for advection 
                    ! of this layer
                    ! (e.g., Greve and Blatter, 2009, Eq. 5.148)
                    uz_star(i,j,k) = uz(i,j,k) + ux_aa*c_x + uy_aa*c_y + c_t 

                    if (abs(uz_star(i,j,k)) .lt. TOL_UNDERFLOW) uz_star(i,j,k) = 0.0_wp
                    
                end do 
                
            else 
                ! No ice here, set vertical velocity equal to negative accum and bedrock change 

                do k = 1, nz_ac

                    uz(i,j,k) = dzbdt_now - max(smb(i,j),0.0)
                    if (abs(uz(i,j,k)) .lt. TOL_UNDERFLOW) uz(i,j,k) = 0.0_wp 

                    uz_star(i,j,k) = uz(i,j,k)

               end do 

            end if 

        end do 
        end do 
        !$omp end parallel do 

        return 

    end subroutine calc_uz_3D_jac

    subroutine calc_uz_3D(uz,uz_star,ux,uy,H_ice,f_ice,f_grnd,smb,bmb,dHdt,dzsdt, &
                                    dzsdx,dzsdy,dzbdx,dzbdy,zeta_aa,zeta_ac,dx,dy,use_bmb,boundaries)
        ! Following algorithm outlined by the Glimmer ice sheet model:
        ! https://www.geos.ed.ac.uk/~mhagdorn/glide/glide-doc/glimmer_htmlse9.html#x17-660003.1.5

        ! Note: rate of ice-base elevation change (dzbdt) is deduced from dzbdt = dzsdt - dhdt. 
        ! This formulation does not depend on rate of bedrock uplift (which is implicit in dzsdt),
        ! and is valid for both grounded and floating ice. 

        implicit none 

        real(wp), intent(OUT) :: uz(:,:,:)          ! nx,ny,nz_ac
        real(wp), intent(OUT) :: uz_star(:,:,:)     ! nx,ny,nz_ac
        real(wp), intent(IN)  :: ux(:,:,:)          ! nx,ny,nz_aa
        real(wp), intent(IN)  :: uy(:,:,:)          ! nx,ny,nz_aa
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: f_grnd(:,:)
        real(wp), intent(IN)  :: smb(:,:) 
        real(wp), intent(IN)  :: bmb(:,:) 
        real(wp), intent(IN)  :: dHdt(:,:) 
        real(wp), intent(IN)  :: dzsdt(:,:) 
        real(wp), intent(IN)  :: dzsdx(:,:) 
        real(wp), intent(IN)  :: dzsdy(:,:) 
        real(wp), intent(IN)  :: dzbdx(:,:) 
        real(wp), intent(IN)  :: dzbdy(:,:) 
        real(wp), intent(IN)  :: zeta_aa(:)    ! z-coordinate, aa-nodes 
        real(wp), intent(IN)  :: zeta_ac(:)    ! z-coordinate, ac-nodes  
        real(wp), intent(IN)  :: dx 
        real(wp), intent(IN)  :: dy
        logical,  intent(IN)  :: use_bmb
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer :: i, j, k, nx, ny, nz_aa, nz_ac
        integer :: im1, ip1, jm1, jp1   
        integer :: kup, kdn
        real(wp) :: f_bmb
        real(wp) :: H_now
        real(wp) :: H_inv
        real(wp) :: dzbdx_aa
        real(wp) :: dzbdy_aa
        real(wp) :: dzsdx_aa
        real(wp) :: dzsdy_aa 
        real(wp) :: dudx_aa
        real(wp) :: dvdy_aa
        real(wp) :: dudz_aa 
        real(wp) :: dvdz_aa 
        real(wp) :: ux_aa 
        real(wp) :: uy_aa 
        real(wp) :: uz_grid 
        real(wp) :: uz_srf 
        real(wp) :: zeta_now 
        real(wp) :: c_x 
        real(wp) :: c_y 
        real(wp) :: c_t 

        real(wp) :: dzsdtn(4)
        real(wp) :: dhdtn(4)
        real(wp) :: dzbdxn(4)
        real(wp) :: dzbdyn(4)
        real(wp) :: dzsdxn(4)
        real(wp) :: dzsdyn(4)
        real(wp) :: dudxn(4) 
        real(wp) :: dvdyn(4) 
        real(wp) :: uxn_up(4) 
        real(wp) :: uxn_dn(4) 
        real(wp) :: uxn(4) 
        real(wp) :: uyn_up(4) 
        real(wp) :: uyn_dn(4)
        real(wp) :: uyn(4)  
        real(wp) :: dudzn(4) 
        real(wp) :: dvdzn(4) 
        
        real(wp) :: wt0
        real(wp) :: xn(4) 
        real(wp) :: yn(4) 
        real(wp) :: wtn(4)
        real(wp) :: wt2D 

        real(wp) :: dzsdt_now
        real(wp) :: dhdt_now
        real(wp) :: dzbdt_now 

        real(wp), allocatable :: dudx(:,:,:)
        real(wp), allocatable :: dvdy(:,:,:)
        real(wp), allocatable :: dudy(:,:)
        real(wp), allocatable :: dvdx(:,:)
        
        real(wp), parameter :: uz_min = -10.0     ! [m/yr] Minimum allowed vertical velocity downwards for stability
        
        nx    = size(ux,1)
        ny    = size(ux,2)
        nz_aa = size(zeta_aa,1)
        nz_ac = size(zeta_ac,1) 
        
        ! Allocate arrays 
        allocate(dudx(nx,ny,nz_aa))
        allocate(dvdy(nx,ny,nz_aa))
        allocate(dudy(nx,ny))
        allocate(dvdx(nx,ny))
        
        ! Get nodes and weighting 
        wt0  = 1.0/sqrt(3.0)
        xn   = [wt0,-wt0,-wt0, wt0]
        yn   = [wt0, wt0,-wt0,-wt0]
        wtn  = [1.0,1.0,1.0,1.0]
        wt2D = 4.0   ! Surface area of square [-1:1,-1:1]=> 2x2 => 4 

        ! Initialize vertical velocity to zero 
        uz = 0.0 

        ! Define switch for bmb
        if (use_bmb) then 
            f_bmb = 1.0 
        else 
            f_bmb = 0.0 
        end if 

        ! First calculate horizontal strain rates at each layer for later use,
        ! with no correction factor for sigma-transformation.
        ! Note: we only need dudx and dvdy, but routine also calculate cross terms, which will not be used.

        !$omp parallel do 
        do k = 1, nz_aa
            call calc_strain_rate_horizontal_2D(dudx(:,:,k),dudy,dvdx,dvdy(:,:,k),ux(:,:,k),uy(:,:,k),f_ice,dx,dy,boundaries)
        end do

        ! Next, calculate vertical velocity at each point through the column

        !$omp parallel do 
        do j = 1, ny
        do i = 1, nx

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

            ! Diagnose rate of ice-base elevation change (needed for all points)
            dzsdt_now = dzsdt(i,j) 
            dhdt_now  = dhdt(i,j) 
            dzbdt_now = dzsdt_now - dhdt_now

            if (f_ice(i,j) .eq. 1.0) then

                H_now  = H_ice(i,j) 
                H_inv = 1.0/H_now 

                ! Get the centered ice-base gradient
                call acx_to_nodes(dzbdxn,dzbdx,i,j,xn,yn,im1,ip1,jm1,jp1)
                dzbdx_aa = sum(dzbdxn*wtn)/wt2D
                
                call acy_to_nodes(dzbdyn,dzbdy,i,j,xn,yn,im1,ip1,jm1,jp1)
                dzbdy_aa = sum(dzbdyn*wtn)/wt2D
                
                ! Get the centered surface gradient
                call acx_to_nodes(dzsdxn,dzsdx,i,j,xn,yn,im1,ip1,jm1,jp1)
                dzsdx_aa = sum(dzsdxn*wtn)/wt2D
                
                call acy_to_nodes(dzsdyn,dzsdy,i,j,xn,yn,im1,ip1,jm1,jp1)
                dzsdy_aa = sum(dzsdyn*wtn)/wt2D
                
                ! Get the aa-node centered horizontal velocity at the base
                call acx_to_nodes(uxn,ux(:,:,1),i,j,xn,yn,im1,ip1,jm1,jp1)
                ux_aa = sum(uxn*wtn)/wt2D
                
                call acy_to_nodes(uyn,uy(:,:,1),i,j,xn,yn,im1,ip1,jm1,jp1)
                uy_aa = sum(uyn*wtn)/wt2D
                
                ! Determine grid vertical velocity at the base due to sigma-coordinates 
                ! Glimmer, Eq. 3.35 
                ! ajr, 2020-01-27, untested:::
!                 uz_grid = dzsdt(i,j) + (ux_aa*dzsdx_aa + uy_aa*dzsdy_aa) &
!                             - ( (1.0_wp-zeta_ac(1))*dHdt(i,j) + ux_aa*dHdx_aa + uy_aa*dHdy_aa )
                uz_grid = 0.0_wp 

                ! ===================================================================
                ! Greve and Blatter (2009) style:

                ! Determine basal vertical velocity for this grid point 
                ! Following Eq. 5.31 of Greve and Blatter (2009)
                uz(i,j,1) = dzbdt_now + uz_grid + f_bmb*bmb(i,j) + ux_aa*dzbdx_aa + uy_aa*dzbdy_aa
                if (abs(uz(i,j,1)) .lt. TOL_UNDERFLOW) uz(i,j,1) = 0.0_wp 
                
                ! Set stability limit on basal uz value.
                ! This only gets applied in rare cases when something
                ! is going wrong in the model. 
                if (uz(i,j,1) .lt. uz_min) uz(i,j,1) = uz_min 

                ! Determine surface vertical velocity following kinematic boundary condition 
                ! Glimmer, Eq. 3.10 [or Folwer, Chpt 10, Eq. 10.8]
                !uz_srf = dzsdt(i,j) + ux_aa*dzsdx_aa + uy_aa*dzsdy_aa - smb(i,j) 
                
                ! Integrate upward to each point above base until just below surface is reached 
                ! Integrate on vertical ac-nodes (ie, vertical cell borders between aa-node centers)
                do k = 2, nz_ac 

                    ! Calculate sigma-coordinate derivative correction factors for this layer
                    ! (Greve and Blatter, 2009, Eqs. 5.131 and 5.132, 
                    !  also shown in 1D with Eq. 5.145)
                    ! See src/physics/deformation.f90::calc_jacobian_vel_3D()

                    ! Take zeta at the center of the cell below the current vertical ac boundary
                    ! (this is also where dudz_aa and dvz_aa are calculated above)
                    zeta_now = zeta_aa(k-1)

                    c_x = -H_inv * ( (1.0-zeta_now)*dzbdx_aa + zeta_now*dzsdx_aa )
                    c_y = -H_inv * ( (1.0-zeta_now)*dzbdy_aa + zeta_now*dzsdy_aa )

                    ! Get dudz/dvdz values at vertical aa-nodes, in order
                    ! to vertically integrate each cell up to ac-node border.
                    ! Note: nz_ac = nz_aa + 1
                    if (k .eq. 2) then
                        kup = k 
                        kdn = k-1 
                    else if (k .eq. nz_ac) then 
                        kup = k-1 
                        kdn = k-2 
                    else 
                        ! Centered on k-1 
                        kup = k 
                        kdn = k-2
                    end if

                    call acx_to_nodes(uxn_up,ux(:,:,kup),i,j,xn,yn,im1,ip1,jm1,jp1)
                    call acx_to_nodes(uxn_dn,ux(:,:,kdn),i,j,xn,yn,im1,ip1,jm1,jp1)
                    dudzn = (uxn_up - uxn_dn) / (zeta_aa(kup)-zeta_aa(kdn))
                    dudz_aa = sum(dudzn*wtn)/wt2D

                    call acy_to_nodes(uyn_up,uy(:,:,kup),i,j,xn,yn,im1,ip1,jm1,jp1)
                    call acy_to_nodes(uyn_dn,uy(:,:,kdn),i,j,xn,yn,im1,ip1,jm1,jp1)
                    dvdzn = (uyn_up - uyn_dn) / (zeta_aa(kup)-zeta_aa(kdn))
                    dvdz_aa = sum(dvdzn*wtn)/wt2D

                    ! Calculate sigma-corrected derivatives
                    call acx_to_nodes(dudxn,dudx(:,:,k-1),i,j,xn,yn,im1,ip1,jm1,jp1)
                    dudx_aa = sum(dudxn*wtn)/wt2D  +  c_x*dudz_aa 

                    call acy_to_nodes(dvdyn,dvdy(:,:,k-1),i,j,xn,yn,im1,ip1,jm1,jp1)
                    dvdy_aa = sum(dvdyn*wtn)/wt2D  +  c_y*dvdz_aa 

                    ! Calculate vertical velocity of this layer
                    ! (Greve and Blatter, 2009, Eq. 5.95)
                    uz(i,j,k) = uz(i,j,k-1) - H_now*(zeta_ac(k)-zeta_ac(k-1))*(dudx_aa+dvdy_aa)

                    ! Apply correction to match kinematic boundary condition at surface 
                    !uz(i,j,k) = uz(i,j,k) - zeta_ac(k)*(uz(i,j,k)-uz_srf)

                    if (abs(uz(i,j,k)) .lt. TOL_UNDERFLOW) uz(i,j,k) = 0.0_wp 
                    
                end do 


                ! === Also calculate adjusted vertical velocity to be used for temperature advection
                
                do k = 1, nz_ac 

                    ! Get the centered horizontal velocity of box on vertical ac-nodes at the level k 
                    ! ajr: Given that the correction is applied to uz, which is defined on 
                    ! ac-nodes, it seems the correction should also be calculated on ac-nodes.
                    ! Note: nz_ac = nz_aa + 1
    
                    if (k .eq. 1) then 
                        kup = k 
                        kdn = k 
                    else if (k .eq. nz_ac) then 
                        kup = k-1 
                        kdn = k-1 
                    else
                        kup = k 
                        kdn = k-1 
                    end if
                    
                    call acx_to_nodes(uxn_up,ux(:,:,kup),i,j,xn,yn,im1,ip1,jm1,jp1)
                    call acx_to_nodes(uxn_dn,ux(:,:,kdn),i,j,xn,yn,im1,ip1,jm1,jp1)
                    uxn = 0.5_wp*(uxn_up+uxn_dn)
                    ux_aa = sum(uxn*wtn)/wt2D
                    
                    call acy_to_nodes(uyn_up,uy(:,:,kup),i,j,xn,yn,im1,ip1,jm1,jp1)
                    call acy_to_nodes(uyn_dn,uy(:,:,kdn),i,j,xn,yn,im1,ip1,jm1,jp1)
                    uyn = 0.5_wp*(uyn_up+uyn_dn)
                    uy_aa = sum(uyn*wtn)/wt2D
                    
                    ! Take zeta directly at vertical cell edge where uz is calculated
                    ! (this is also where ux_aa and uy_aa are calculated above)
                    zeta_now = zeta_ac(k)

                    ! Calculate sigma-coordinate derivative correction factors
                    ! (Greve and Blatter, 2009, Eqs. 5.131 and 5.132, 
                    !  also shown in 1D with Eq. 5.145)

                    ! Note: not dividing by H here, since this is done in the thermodynamics advection step
                    c_x = -( (1.0-zeta_now)*dzbdx_aa  + zeta_now*dzsdx_aa )
                    c_y = -( (1.0-zeta_now)*dzbdy_aa  + zeta_now*dzsdy_aa )
                    c_t = -( (1.0-zeta_now)*dzbdt_now + zeta_now*dzsdt_now )
                    
                    ! Calculate adjusted vertical velocity for advection 
                    ! of this layer
                    ! (e.g., Greve and Blatter, 2009, Eq. 5.148)
                    uz_star(i,j,k) = uz(i,j,k) + ux_aa*c_x + uy_aa*c_y + c_t 

                    if (abs(uz_star(i,j,k)) .lt. TOL_UNDERFLOW) uz_star(i,j,k) = 0.0_wp
                    
                end do 
                
            else 
                ! No ice here, set vertical velocity equal to negative accum and bedrock change 

                do k = 1, nz_ac

                    uz(i,j,k) = dzbdt_now - max(smb(i,j),0.0)
                    if (abs(uz(i,j,k)) .lt. TOL_UNDERFLOW) uz(i,j,k) = 0.0_wp 

                    uz_star(i,j,k) = uz(i,j,k)

               end do 

            end if 

        end do 
        end do 
        !$omp end parallel do 

        return 

    end subroutine calc_uz_3D

    subroutine calc_uz_3D_aa(uz,uz_star,ux,uy,H_ice,f_ice,f_grnd,z_bed,z_srf,smb,bmb,dHdt,dzsdt, &
                                            dzsdx,dzsdy,dzbdx,dzbdy,zeta_aa,zeta_ac,dx,dy,use_bmb,boundaries)
        ! Following algorithm outlined by the Glimmer ice sheet model:
        ! https://www.geos.ed.ac.uk/~mhagdorn/glide/glide-doc/glimmer_htmlse9.html#x17-660003.1.5

        ! Note: rate of bedrock uplift (dzbdt) no longer considered, since the rate is 
        ! very small and now z_bed is updated externally (ie, now assume dzbdt = 0.0 here)

        implicit none 

        real(wp), intent(OUT) :: uz(:,:,:)          ! nx,ny,nz_ac
        real(wp), intent(OUT) :: uz_star(:,:,:)     ! nx,ny,nz_ac
        real(wp), intent(IN)  :: ux(:,:,:)          ! nx,ny,nz_aa
        real(wp), intent(IN)  :: uy(:,:,:)          ! nx,ny,nz_aa
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: f_grnd(:,:)
        real(wp), intent(IN)  :: z_bed(:,:) 
        real(wp), intent(IN)  :: z_srf(:,:) 
        real(wp), intent(IN)  :: smb(:,:) 
        real(wp), intent(IN)  :: bmb(:,:) 
        real(wp), intent(IN)  :: dHdt(:,:) 
        real(wp), intent(IN)  :: dzsdt(:,:)
        real(wp), intent(IN)  :: dzsdx(:,:) 
        real(wp), intent(IN)  :: dzsdy(:,:) 
        real(wp), intent(IN)  :: dzbdx(:,:) 
        real(wp), intent(IN)  :: dzbdy(:,:) 
        real(wp), intent(IN)  :: zeta_aa(:)    ! z-coordinate, aa-nodes 
        real(wp), intent(IN)  :: zeta_ac(:)    ! z-coordinate, ac-nodes  
        real(wp), intent(IN)  :: dx 
        real(wp), intent(IN)  :: dy
        logical,  intent(IN)  :: use_bmb
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer  :: i, j, k, nx, ny, nz_aa, nz_ac
        integer  :: im1, ip1, jm1, jp1   
        real(wp) :: f_bmb 
        real(wp) :: H_now
        real(wp) :: H_inv
        real(wp) :: dzbdx_aa
        real(wp) :: dzbdy_aa
        real(wp) :: dzsdx_aa
        real(wp) :: dzsdy_aa
        real(wp) :: duxdx_aa
        real(wp) :: duydy_aa
        real(wp) :: duxdz_aa 
        real(wp) :: duydz_aa
        real(wp) :: duxdx_now 
        real(wp) :: duydy_now 
        real(wp) :: ux_aa 
        real(wp) :: uy_aa 
        real(wp) :: uz_grid 
        real(wp) :: uz_srf 
        real(wp) :: zeta_now 
        real(wp) :: c_x 
        real(wp) :: c_y 
        real(wp) :: c_t 

        real(wp) :: dzsdt_now
        real(wp) :: dhdt_now
        real(wp) :: dzbdt_now 

        real(wp), parameter :: dzbdt        = 0.0     ! For posterity, keep dzbdt variable, but set to zero 
        real(wp), parameter :: uz_min       = -10.0   ! [m/yr] Minimum allowed vertical velocity downwards for stability
        
        nx    = size(ux,1)
        ny    = size(ux,2)
        nz_aa = size(zeta_aa,1)
        nz_ac = size(zeta_ac,1) 

        ! Initialize vertical velocity to zero 
        uz = 0.0 

        ! Define switch for bmb
        if (use_bmb) then 
            f_bmb = 1.0 
        else 
            f_bmb = 0.0 
        end if 
        
        ! Next, calculate velocity 

        !$omp parallel do 
        do j = 1, ny
        do i = 1, nx

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

            ! Diagnose rate of ice-base elevation change (needed for all points)
            dzsdt_now = dzsdt(i,j) 
            dhdt_now  = dhdt(i,j) 
            dzbdt_now = dzsdt_now - dhdt_now

            if (f_ice(i,j) .eq. 1.0) then

                ! Get weighted ice thickness for stability
!                 H_now = (4.0*H_ice(i,j) + 2.0*(H_ice(im1,j)+H_ice(ip1,j)+H_ice(i,jm1)+H_ice(i,jp1))) / 16.0 &
!                       + (H_ice(im1,jm1)+H_ice(ip1,jm1)+H_ice(ip1,jp1)+H_ice(im1,jp1)) / 16.0 

                H_now  = H_ice(i,j) 
                H_inv = 1.0/H_now 

                ! Get the centered bedrock gradient 
                dzbdx_aa = 0.5*(dzbdx(i,j)+dzbdx(im1,j))
                dzbdy_aa = 0.5*(dzbdy(i,j)+dzbdy(i,jm1))

                ! Get the centered surface gradient 
                dzsdx_aa = 0.5*(dzsdx(i,j)+dzsdx(im1,j))
                dzsdy_aa = 0.5*(dzsdy(i,j)+dzsdy(i,jm1))

                ! Get the centered horizontal velocity at the base
                ux_aa = 0.5_wp* (ux(im1,j,1) + ux(i,j,1))
                uy_aa = 0.5_wp* (uy(i,jm1,1) + uy(i,j,1))
                
                ! Determine grid vertical velocity at the base due to sigma-coordinates 
                ! Glimmer, Eq. 3.35 
                ! ajr, 2020-01-27, untested:::
!                 uz_grid = dzsdt(i,j) + (ux_aa*dzsdx_aa + uy_aa*dzsdy_aa) &
!                             - ( (1.0_wp-zeta_ac(1))*dHdt(i,j) + ux_aa*dHdx_aa + uy_aa*dHdy_aa )
                uz_grid = 0.0_wp 

                ! ===================================================================
                ! Greve and Blatter (2009) style:

                ! Determine basal vertical velocity for this grid point 
                ! Following Eq. 5.31 of Greve and Blatter (2009)
                uz(i,j,1) = dzbdt + uz_grid + f_bmb*bmb(i,j) + ux_aa*dzbdx_aa + uy_aa*dzbdy_aa
                if (abs(uz(i,j,1)) .lt. TOL_UNDERFLOW) uz(i,j,1) = 0.0_wp 
                
                ! Set stability limit on basal uz value for grounded ice.
                ! This only gets applied in rare cases when something
                ! is going wrong in the model. 
                if (f_grnd(i,j) .eq. 1.0 .and. uz(i,j,1) .lt. uz_min) uz(i,j,1) = uz_min 

                ! Determine surface vertical velocity following kinematic boundary condition 
                ! Glimmer, Eq. 3.10 [or Folwer, Chpt 10, Eq. 10.8]
                !uz_srf = dzsdt(i,j) + ux_aa*dzsdx_aa + uy_aa*dzsdy_aa - smb(i,j) 
                
                ! Integrate upward to each point above base until just below surface is reached 
                do k = 2, nz_ac 

                    ! Greve and Blatter (2009), Eq. 5.72
                    ! Bueler and Brown  (2009), Eq. 4
                    
                    duxdx_aa  = (ux(i,j,k-1)   - ux(im1,j,k-1)  )/dx
                    duydy_aa  = (uy(i,j,k-1)   - uy(i,jm1,k-1)  )/dy

                    ! Note: nz_ac = nz_aa + 1
                    if (k .eq. 2) then
                        duxdz_aa  = ( 0.5*(ux(i,j,k)+ux(im1,j,k)) &
                            - 0.5*(ux(i,j,k-1)+ux(im1,j,k-1)) ) / (zeta_aa(k)-zeta_aa(k-1))
                        duydz_aa  = ( 0.5*(uy(i,j,k)+uy(i,jm1,k)) &
                            - 0.5*(uy(i,j,k-1)+uy(i,jm1,k-1)) ) / (zeta_aa(k)-zeta_aa(k-1))
                    else if (k .eq. nz_ac) then
                        duxdz_aa  = ( 0.5*(ux(i,j,k-1)+ux(im1,j,k-1)) &
                            - 0.5*(ux(i,j,k-2)+ux(im1,j,k-2)) ) / (zeta_aa(k-1)-zeta_aa(k-2))
                        duydz_aa  = ( 0.5*(uy(i,j,k-1)+uy(i,jm1,k-1)) &
                            - 0.5*(uy(i,j,k-2)+uy(i,jm1,k-2)) ) / (zeta_aa(k-1)-zeta_aa(k-2))
                    else
                        ! Centered difference vertically, centered on k-1
                        duxdz_aa  = ( 0.5*(ux(i,j,k)+ux(im1,j,k)) &
                            - 0.5*(ux(i,j,k-2)+ux(im1,j,k-2)) ) / (zeta_aa(k)-zeta_aa(k-2))
                        duydz_aa  = ( 0.5*(uy(i,j,k)+uy(i,jm1,k)) &
                            - 0.5*(uy(i,j,k-2)+uy(i,jm1,k-2)) ) / (zeta_aa(k)-zeta_aa(k-2))
                    end if 

                    ! Calculate sigma-coordinate derivative correction factors
                    ! (Greve and Blatter, 2009, Eqs. 5.131 and 5.132)
                    c_x = -H_inv * ( (1.0-zeta_ac(k))*dzbdx_aa + zeta_ac(k)*dzsdx_aa )
                    c_y = -H_inv * ( (1.0-zeta_ac(k))*dzbdy_aa + zeta_ac(k)*dzsdy_aa )

                    ! Calculate derivatives
                    duxdx_now = duxdx_aa + c_x*duxdz_aa 
                    duydy_now = duydy_aa + c_y*duydz_aa 
                    
                    ! Calculate vertical velocity of this layer
                    ! (Greve and Blatter, 2009, Eq. 5.95)
                    uz(i,j,k) = uz(i,j,k-1) & 
                        - H_now*(zeta_ac(k)-zeta_ac(k-1))*(duxdx_now+duydy_now)

                    ! Apply correction to match kinematic boundary condition at surface 
                    !uz(i,j,k) = uz(i,j,k) - zeta_ac(k)*(uz(i,j,k)-uz_srf)

                    if (abs(uz(i,j,k)) .lt. TOL_UNDERFLOW) uz(i,j,k) = 0.0_wp 
                    
                end do 
                
                ! === Also calculate adjusted vertical velocity to be used for temperature advection
                
                do k = 1, nz_ac 

                    ! Get the centered horizontal velocity of box
                    ! on vertical ac-nodes at the level k 
                    ! ajr: Given that the correction is 
                    ! applied to uz, which is defined on 
                    ! ac-nodes, it seems the correction
                    ! should also be calculated on ac-nodes.
                    ! Note: nz_ac = nz_aa + 1
    
                    if (k .eq. 1) then 
                        ux_aa = 0.5_wp* (ux(im1,j,k) + ux(i,j,k))
                        uy_aa = 0.5_wp* (uy(i,jm1,k) + uy(i,j,k))
                    else if (k .eq. nz_ac) then  
                        ux_aa = 0.5_wp* (ux(im1,j,k-1) + ux(i,j,k-1))
                        uy_aa = 0.5_wp* (uy(i,jm1,k-1) + uy(i,j,k-1))
                    else
                        ux_aa = 0.25_wp* (ux(im1,j,k) + ux(i,j,k) + ux(im1,j,k-1) + ux(i,j,k-1))
                        uy_aa = 0.25_wp* (uy(i,jm1,k) + uy(i,j,k) + uy(i,jm1,k-1) + uy(i,j,k-1))
                    end if 

                    ! Take zeta directly at vertical cell edge where uz is calculated
                    ! (this is also where ux_aa and uy_aa are calculated above)
                    zeta_now = zeta_ac(k)

                    ! Calculate sigma-coordinate derivative correction factors
                    ! (Greve and Blatter, 2009, Eqs. 5.131 and 5.132, 
                    !  also shown in 1D with Eq. 5.145)

                    ! Note: not dividing by H here, since this is done in the thermodynamics advection step
                    c_x = -ux_aa * ( (1.0_wp-zeta_now)*dzbdx_aa  + zeta_now*dzsdx_aa )
                    c_y = -uy_aa * ( (1.0_wp-zeta_now)*dzbdy_aa  + zeta_now*dzsdy_aa )
                    c_t =         -( (1.0_wp-zeta_now)*dzbdt_now + zeta_now*dzsdt_now )

                    ! Calculate adjusted vertical velocity for advection 
                    ! of this layer
                    ! (e.g., Greve and Blatter, 2009, Eq. 5.148)
                    uz_star(i,j,k) = uz(i,j,k) + c_x + c_y + c_t 

                    if (abs(uz_star(i,j,k)) .lt. TOL_UNDERFLOW) uz_star(i,j,k) = 0.0_wp 
                    
                end do 
                
            else 
                ! No ice here, set vertical velocity equal to negative accum and bedrock change 

                do k = 1, nz_ac 

                    uz(i,j,k) = dzbdt - max(smb(i,j),0.0)
                    if (abs(uz(i,j,k)) .lt. TOL_UNDERFLOW) uz(i,j,k) = 0.0_wp 

                    uz_star(i,j,k) = uz(i,j,k)

               end do 

            end if 

        end do 
        end do 
        !$omp end parallel do 

        return 

    end subroutine calc_uz_3D_aa

    subroutine calc_driving_stress(taud_acx,taud_acy,H_ice,f_ice,dzsdx,dzsdy,dx,taud_lim,boundaries)
        ! Calculate driving stress on staggered grid points
        ! Units: taud [Pa] == [kg m-1 s-2]
        ! taud = rho_ice*g*H_ice*dzs/dx

        ! Note: interpolation to Ab nodes no longer used here.

        ! Note: use parameter taud_lim to limit maximum allowed driving stress magnitude applied in the model.
        ! Should be an extreme value. eg, if dzdx = taud_lim / (rho*g*H), 
        ! then for taud_lim=5e5 and H=3000m, dzdx = 2e5 / (910*9.81*3000) = 0.02, 
        ! which is a rather steep slope for a shallow-ice model.

        implicit none 

        real(wp), intent(OUT) :: taud_acx(:,:)    ! [Pa]
        real(wp), intent(OUT) :: taud_acy(:,:)    ! [Pa]
        real(wp), intent(IN)  :: H_ice(:,:)       ! [m]
        real(wp), intent(IN)  :: f_ice(:,:)       ! [--]
        real(wp), intent(IN)  :: dzsdx(:,:)       ! [--]
        real(wp), intent(IN)  :: dzsdy(:,:)       ! [--]
        real(wp), intent(IN)  :: dx               ! [m] 
        real(wp), intent(IN)  :: taud_lim         ! [Pa]
        character(len=*), intent(IN) :: boundaries  ! Boundary conditions to apply 

        ! Local variables 
        integer :: i, j, nx, ny 
        integer    :: im1, ip1, jm1, jp1 
        real(wp) :: dy, rhog 
        real(wp) :: H_mid
        real(wp) :: taud_mean 
        real(wp) :: taud_eps 

        nx = size(H_ice,1)
        ny = size(H_ice,2) 

        ! Define shortcut parameter 
        rhog = rho_ice * g 

        ! Assume grid resolution is symmetrical 
        dy = dx 

        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

            ! x-direction
            if (f_ice(i,j) .eq. 1.0_wp .and. f_ice(ip1,j) .lt. 1.0_wp) then 
                H_mid = 0.5_wp*(H_ice(i,j)+0.0_wp)
            else if (f_ice(i,j) .lt. 1.0_wp .and. f_ice(ip1,j) .eq. 1.0_wp) then 
                H_mid = 0.5_wp*(0.0_wp+H_ice(ip1,j))
            else  
                H_mid = 0.5_wp*(H_ice(i,j)+H_ice(ip1,j)) 
            end if
            taud_acx(i,j) = rhog * H_mid * dzsdx(i,j) 

            ! Apply limit
            if (taud_acx(i,j) .gt.  taud_lim) taud_acx(i,j) =  taud_lim 
            if (taud_acx(i,j) .lt. -taud_lim) taud_acx(i,j) = -taud_lim 
            
            ! y-direction 
            if (f_ice(i,j) .eq. 1.0_wp .and. f_ice(i,jp1) .lt. 1.0_wp) then 
                H_mid = 0.5_wp*(H_ice(i,j)+0.0_wp)
            else if (f_ice(i,j) .lt. 1.0_wp .and. f_ice(i,jp1) .eq. 1.0_wp) then 
                H_mid = 0.5_wp*(0.0_wp+H_ice(i,jp1))
            else  
                H_mid = 0.5_wp*(H_ice(i,j)+H_ice(i,jp1)) 
            end if
            taud_acy(i,j) = rhog * H_mid * dzsdy(i,j) 

            ! Apply limit
            if (taud_acy(i,j) .gt.  taud_lim) taud_acy(i,j) =  taud_lim 
            if (taud_acy(i,j) .lt. -taud_lim) taud_acy(i,j) = -taud_lim 
            
        end do
        end do 

        ! SPECIAL CASES: edge cases for stability...

        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

            ! === x-direction ===

            if (f_ice(i,j) .eq. 1.0_wp .and. f_ice(ip1,j) .eq. 1.0_wp) then 

                taud_mean = 0.5*(taud_acx(im1,j)+taud_acx(ip1,j))
                taud_eps  = (taud_acx(i,j) - taud_mean) / taud_lim

                if (abs(taud_acx(i,j)) .eq. taud_lim .and. abs(taud_eps) .gt. 0.4) then 
                    ! Set driving stress equal to weighted average of neighbors
                    taud_acx(i,j) = 0.3*taud_acx(i,j) + 0.35*taud_acx(im1,j)+ 0.35*taud_acx(ip1,j)
                end if 


            end if

            ! === y-direction ===

            if (f_ice(i,j) .eq. 1.0_wp .and. f_ice(i,jp1) .eq. 1.0_wp) then 

                taud_mean = 0.5*(taud_acy(i,jm1)+taud_acy(i,jp1))
                taud_eps  = (taud_acy(i,j) - taud_mean) / taud_lim

                if (abs(taud_acy(i,j)) .eq. taud_lim .and. abs(taud_eps) .gt. 0.4) then 
                    ! Set driving stress equal to weighted average of neighbors
                    taud_acy(i,j) = 0.30*taud_acy(i,j) + 0.35*taud_acy(i,jm1)+ 0.35*taud_acy(i,jp1)
                end if 


            end if
            
        end do
        end do 

        return 

    end subroutine calc_driving_stress

    subroutine calc_driving_stress_gl(taud_acx,taud_acy,H_ice,z_srf,z_bed,z_sl,H_grnd, &
                                      f_grnd,f_grnd_acx,f_grnd_acy,dx,method,beta_gl_stag)
        ! taud = rho_ice*g*H_ice
        ! Calculate driving stress on staggered grid points, with 
        ! special treatment of the grounding line 
        ! Units: taud [Pa] == [kg m-1 s-2]
        
        ! Note: interpolation to Ab nodes no longer used here.

        implicit none 

        real(wp), intent(OUT) :: taud_acx(:,:)
        real(wp), intent(OUT) :: taud_acy(:,:) 
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: z_srf(:,:)
        real(wp), intent(IN)  :: z_bed(:,:)
        real(wp), intent(IN)  :: z_sl(:,:)
        real(wp), intent(IN)  :: H_grnd(:,:)
        real(wp), intent(IN)  :: f_grnd(:,:)
        real(wp), intent(IN)  :: f_grnd_acx(:,:)
        real(wp), intent(IN)  :: f_grnd_acy(:,:)
        real(wp), intent(IN)  :: dx 
        integer,    intent(IN)  :: method        ! Which driving stress calculation to use
        integer,    intent(IN)  :: beta_gl_stag  ! Method of grounding line staggering of beta 

        ! Local variables 
        integer :: i, j, nx, ny
        integer :: im1, ip1, jm1, jp1  
        real(wp) :: dy, rhog 
        real(wp) :: taud_grnd, taud_flt, taud_now 
        real(wp) :: H_mid, H_gl, z_gl, H_grnd_mid 
        real(wp) :: dzsdx, dzsdy
        real(wp) :: dzsdx_1, dzsdx_2
        real(wp) :: H_1, H_2  
        real(wp) :: taud_old, fac_gl   

        real(wp), parameter :: slope_max = 0.05   ! Very high limit == 0.05, low limit < 0.01 

        nx = size(H_ice,1)
        ny = size(H_ice,2) 

        ! Define shortcut parameter 
        rhog = rho_ice * g 

        ! Assume grid resolution is symmetrical 
        dy = dx 

        ! === Subgrid treatment === 
        
        ! Next, refine the driving stress at the grounding line
        ! as desired 

        select case(method)

            case(-1)
                ! One-sided choice
                ! between upstream or downstream driving stress 
                ! at grounding line

                if (beta_gl_stag .eq. 1) then 
                    ! Upstream beta assigned at gl (ie, beta=beta_upstream)

if (.FALSE.) then 
                    ! x-direction 
                    do j = 1, ny 
                    do i = 2, nx-1 

                        if (f_grnd(i,j) .eq. 0.0 .and. f_grnd(i+1,j) .gt. 0.0) then 
                            taud_acx(i,j) = taud_acx(i+1,j) 
                        else if (f_grnd(i,j) .gt. 0.0 .and. f_grnd(i+1,j) .eq. 0.0) then  
                            taud_acx(i,j) = taud_acx(i-1,j)
                        end if 

                    end do 
                    end do 

                    ! y-direction 
                    do j = 2, ny-1 
                    do i = 1, nx 

                        if (f_grnd(i,j) .eq. 0.0 .and. f_grnd(i,j+1) .gt. 0.0) then 
                            taud_acy(i,j) = taud_acy(i,j+1) 
                        else if (f_grnd(i,j) .gt. 0.0 .and. f_grnd(i,j+1) .eq. 0.0) then  
                            taud_acy(i,j) = taud_acy(i,j-1)
                        end if 

                    end do 
                    end do 
end if 
                else if (beta_gl_stag .eq. 2) then 
                    ! Downstream beta assigned at gl (ie, beta=0)

                    ! x-direction 
                    do j = 1, ny 
                    do i = 2, nx-1 

                        if (f_grnd(i,j) .eq. 0.0 .and. f_grnd(i+1,j) .gt. 0.0) then 
                            taud_acx(i,j) = taud_acx(i-1,j) 
                        else if (f_grnd(i,j) .gt. 0.0 .and. f_grnd(i+1,j) .eq. 0.0) then  
                            taud_acx(i,j) = taud_acx(i+1,j)
                        end if 

                    end do 
                    end do 

                    ! y-direction 
                    do j = 2, ny-1 
                    do i = 1, nx 

                        if (f_grnd(i,j) .eq. 0.0 .and. f_grnd(i,j+1) .gt. 0.0) then 
                            taud_acy(i,j) = taud_acy(i,j-1) 
                        else if (f_grnd(i,j) .gt. 0.0 .and. f_grnd(i,j+1) .eq. 0.0) then  
                            taud_acy(i,j) = taud_acy(i,j+1)
                        end if 

                    end do 
                    end do 

                else 

                    write(*,*) "calc_driving_stress_gl:: Error: Wrong choice of beta_gl_stag for this method."
                    stop 

                end if 

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

                        ! Get the ice thickness at the ac-node as the average of two neighbors
                        H_gl    = 0.5_wp*(H_ice(i,j)+H_ice(i+1,j))

                        ! Get slope of grounded point and virtual floating point (using H_ice),
                        ! then assume slope is the weighted average of the two 
                        dzsdx_1 = (z_srf(i+1,j)-z_srf(i,j)) / dx 
                        dzsdx_2 = 0.0 !(H_ice(i+1,j)-H_ice(i,j)) / dx 
                        dzsdx   = f_grnd_acx(i,j)*dzsdx_1 + (1.0-f_grnd_acx(i,j))*dzsdx_2  
                        
                        ! Limit the slope
                        call minmax(dzsdx,slope_max)  
                                                    
                        ! Get the driving stress
                        taud_old = taud_acx(i,j) 
                        taud_acx(i,j) = rhog * H_gl * dzsdx
                        
                        if (j .eq. 6) then 
                            write(*,"(a,i3,12g12.3)") "taud: ", i, f_grnd_acx(i,j), taud_old, taud_acx(i,j)
                        end if 

                    end if 

                end do 
                end do 

                ! y-direction 
                do j = 1, ny-1 
                do i = 1, nx 

                    if ( f_grnd_acy(i,j) .gt. 0.0 .and. f_grnd_acy(i,j) .lt. 1.0) then 
                        ! Grounding line point (ac-node)

                        ! Get the ice thickness at the ac-node as the average of two neighbors
                        H_gl    = 0.5_wp*(H_ice(i,j)+H_ice(i,j+1))

                        ! Get slope of grounded point and virtual floating point (using H_ice),
                        ! then assume slope is the weighted average of the two 
                        dzsdx_1 = (z_srf(i,j+1)-z_srf(i,j)) / dx 
                        dzsdx_2 = 0.0 !(H_ice(i,j+1)-H_ice(i,j)) / dx 
                        dzsdy   = f_grnd_acy(i,j)*dzsdx_1 + (1.0-f_grnd_acy(i,j))*dzsdx_2  
                        
                        call minmax(dzsdy,slope_max)  

                        ! Get the driving stress
                        taud_acy(i,j) = rhog * H_gl * dzsdy
                        
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

                        H_grnd_mid = 0.5_wp*(H_grnd(i,j) + H_grnd(i+1,j))

                        if (H_grnd_mid .gt. 0.0) then 
                            ! Consider grounded 
                            dzsdx = (z_srf(i+1,j)-z_srf(i,j)) / dx 
                        else 
                            ! Consider floating 
                            dzsdx = (H_ice(i+1,j)-H_ice(i,j)) / dx
                        end if 
                        call minmax(dzsdx,slope_max)  

                        ! Get the ice thickness at the ac-node
                        H_gl    = 0.5_wp*(H_ice(i,j)+H_ice(i+1,j))

                        ! Get the driving stress
                        taud_old = taud_acx(i,j) 
                        taud_acx(i,j) = rhog * H_gl * dzsdx
                        
                        if (j .eq. 6) then 
                            write(*,"(a,i3,12g12.3)") "taud: ", i, f_grnd_acx(i,j), taud_old, taud_acx(i,j)
                        end if 

                    end if 

                end do 
                end do 

                ! y-direction 
                do j = 1, ny-1 
                do i = 1, nx 

                    if ( f_grnd_acy(i,j) .gt. 0.0 .and. f_grnd_acy(i,j) .lt. 1.0) then 
                        ! Grounding line point (ac-node)

                        H_grnd_mid = 0.5_wp*(H_grnd(i,j) + H_grnd(i,j+1))

                        if (H_grnd_mid .gt. 0.0) then 
                            ! Consider grounded 
                            dzsdx = (z_srf(i,j+1)-z_srf(i,j)) / dx 
                        else 
                            ! Consider floating 
                            dzsdx = (H_ice(i,j+1)-H_ice(i,j)) / dx
                        end if 
                        call minmax(dzsdx,slope_max)  

                        ! Get the ice thickness at the ac-node
                        H_gl    = 0.5_wp*(H_ice(i,j)+H_ice(i,j+1))

                        ! Get the driving stress
                        taud_acy(i,j) = rhog * H_gl * dzsdx
                        
                    end if 

                end do 
                end do 

            case(3) 
                ! Linear interpolation following Gladstone et al. (2010, TC) Eq. 27

                 
                do j = 1, ny 
                do i = 1, nx

                    im1 = max(1, i-1)
                    ip1 = min(nx,i+1)
                    
                    jm1 = max(1, j-1)
                    jp1 = min(ny,j+1)

                    ! x-direction
                    if (H_grnd(i,j) .gt. 0.0 .and. H_grnd(ip1,j) .le. 0.0) then 
                        ! Grounding line point 

                        ! (i,j) grounded; (ip1,j) floating
                        call integrate_gl_driving_stress_linear(taud_acx(i,j),H_ice(i,j),H_ice(ip1,j), &
                                                    z_bed(i,j),z_bed(ip1,j),z_sl(i,j), z_sl(ip1,j),dx)
                    
                    else if (H_grnd(i,j) .le. 0.0 .and. H_grnd(ip1,j) .gt. 0.0) then 
                        ! (i,j) floating; (ip1,j) grounded 

                        call integrate_gl_driving_stress_linear(taud_acx(i,j),H_ice(ip1,j),H_ice(i,j), &
                                                    z_bed(ip1,j),z_bed(i,j),z_sl(ip1,j),z_sl(i,j),dx)
                        
                        ! Set negative for direction 
                        taud_acx(i,j) = -taud_acx(i,j)

                    end if 

                    ! y-direction 
                    if (H_grnd(i,j) .gt. 0.0 .and. H_grnd(i,jp1) .le. 0.0) then
                        ! Grounding line point
 
                        ! (i,j) grounded; (i,jp1) floating 
                        call integrate_gl_driving_stress_linear(taud_acy(i,j),H_ice(i,j),H_ice(i,jp1), &
                                                    z_bed(i,j),z_bed(i,jp1),z_sl(i,j),z_sl(i,jp1),dx)

                    else if (H_grnd(i,j) .le. 0.0 .and. H_grnd(i,jp1) .gt. 0.0) then
                        ! (i,j) floating; (i,jp1) grounded

                        call integrate_gl_driving_stress_linear(taud_acy(i,j),H_ice(i,jp1),H_ice(i,j), &
                                                    z_bed(i,jp1),z_bed(i,j),z_sl(i,jp1),z_sl(i,j),dx)
                        
                        ! Set negative for direction 
                        taud_acy(i,j) = -taud_acy(i,j) 
                    end if

                end do 
                end do 

            case DEFAULT  
                ! Do nothing, use the standard no-subgrid treatment 
        end select 

        return 

    end subroutine calc_driving_stress_gl

    subroutine calc_lateral_bc_stress_2D(tau_bc_int_acx,tau_bc_int_acy,mask_frnt,H_ice,f_ice,z_srf,z_sl,rho_ice,rho_sw,boundaries)
            ! Calculate the vertically integrated lateral stress [Pa m] boundary condition
            ! at the ice front. 

        implicit none 

        real(wp), intent(OUT) :: tau_bc_int_acx(:,:) 
        real(wp), intent(OUT) :: tau_bc_int_acy(:,:) 
        integer,  intent(IN)  :: mask_frnt(:,:) 
        real(wp), intent(IN)  :: H_ice(:,:) 
        real(wp), intent(IN)  :: f_ice(:,:) 
        real(wp), intent(IN)  :: z_srf(:,:) 
        real(wp), intent(IN)  :: z_sl(:,:) 
        real(wp), intent(IN)  :: rho_ice 
        real(wp), intent(IN)  :: rho_sw 
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer  :: i, j, nx, ny
        integer  :: im1, ip1, jm1, jp1 
        integer  :: i1, j1  
        real(wp) :: H_ice_now 
        real(wp) :: z_srf_now 
        real(wp) :: z_sl_now 
        
        nx = size(tau_bc_int_acx,1) 
        ny = size(tau_bc_int_acx,2) 

        ! Intialize boundary fields to zero 
        tau_bc_int_acx = 0.0 
        tau_bc_int_acy = 0.0 

        do j = 1, ny
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

            ! == acx nodes == 

            ! if ( (f_ice(i,j) .eq. 1.0 .and. f_ice(ip1,j) .lt. 1.0) .or. &
            !      (f_ice(i,j) .lt. 1.0 .and. f_ice(ip1,j) .eq. 1.0) ) then 
            !     ! Ice margin detected, proceed to calculating the boundary stress

            if ( (mask_frnt(i,j) .gt. 0 .and. mask_frnt(ip1,j) .lt. 0) .or. & 
                 (mask_frnt(i,j) .lt. 0 .and. mask_frnt(ip1,j) .gt. 0) ) then 
                ! Ice margin detected, proceed to calculating the boundary stress

                ! Determine index of the ice covered cell.
                if (mask_frnt(i,j) .lt. 0.0) then 
                    i1 = ip1 
                else
                    i1 = i 
                end if 

                ! Get current ice thickness, surface elevation 
                ! and sea level at ice front from aa-node values.
                ! (No f_ice scaling since f_ice=1)
                H_ice_now = H_ice(i1,j)     
                z_srf_now = z_srf(i1,j) 
                z_sl_now  = z_sl(i1,j) 

                ! Calculate the lateral stress bc for this point
                call calc_lateral_bc_stress(tau_bc_int_acx(i,j),H_ice_now, &
                                                z_srf_now,z_sl_now,rho_ice,rho_sw)

            end if 


            ! == acy nodes == 

            ! if ( (f_ice(i,j) .eq. 1.0 .and. f_ice(i,jp1) .lt. 1.0) .or. &
            !      (f_ice(i,j) .lt. 1.0 .and. f_ice(i,jp1) .eq. 1.0) ) then 
            !     ! Ice margin detected, proceed to calculating the boundary stress

            if ( (mask_frnt(i,j) .gt. 0 .and. mask_frnt(i,jp1) .lt. 0) .or. & 
                 (mask_frnt(i,j) .lt. 0 .and. mask_frnt(i,jp1) .gt. 0) ) then 
                ! Ice margin detected, proceed to calculating the boundary stress
                
                ! Determine index of the ice covered cell.
                if (mask_frnt(i,j) .lt. 0.0) then 
                    j1 = jp1 
                else
                    j1 = j
                end if 

                ! Get current ice thickness, surface elevation 
                ! and sea level at ice front from aa-node values.
                ! (No f_ice scaling since f_ice=1)
                H_ice_now = H_ice(i,j1)     
                z_srf_now = z_srf(i,j1) 
                z_sl_now  = z_sl(i,j1) 

                ! Calculate the lateral stress bc for this point
                call calc_lateral_bc_stress(tau_bc_int_acy(i,j),H_ice_now, &
                                                z_srf_now,z_sl_now,rho_ice,rho_sw)

            end if 

        end do 
        end do

        return

    end subroutine calc_lateral_bc_stress_2D
    
    subroutine calc_lateral_bc_stress(tau_bc_int,H_ice,z_srf,z_sl,rho_ice,rho_sw)
            ! Calculate the vertically integrated lateral stress [Pa m] boundary condition
            ! at the ice front for given conditions.

            ! =========================================================
            ! Generalized solution for all ice fronts (floating and grounded)
            ! See Lipscomb et al. (2019), Eqs. 11 & 12, and 
            ! Winkelmann et al. (2011), Eq. 27 

        implicit none 

        real(wp), intent(OUT) :: tau_bc_int
        real(wp), intent(IN)  :: H_ice 
        real(wp), intent(IN)  :: z_srf
        real(wp), intent(IN)  :: z_sl
        real(wp), intent(IN)  :: rho_ice 
        real(wp), intent(IN)  :: rho_sw 

        ! Local variables 
        real(wp) :: f_submerged 
        real(wp) :: H_ocn 

        ! Get current ocean thickness bordering ice sheet
        ! (for bedrock above sea level, this will give zero)
        f_submerged = 1.d0 - min((z_srf-z_sl)/H_ice,1.d0)
        H_ocn       = H_ice*f_submerged
        
        tau_bc_int = 0.5d0*rho_ice*g*H_ice**2 &         ! tau_out_int                                                ! p_out
                   - 0.5d0*rho_sw *g*H_ocn**2           ! tau_in_int

        return

    end subroutine calc_lateral_bc_stress

    subroutine integrate_gl_driving_stress_linear(taud,H_a,H_b,zb_a,zb_b,z_sl_a,z_sl_b,dx)
        ! Compute the driving stress for the grounding line more precisely (subgrid)
        ! following Gladstone et al. (2010, TC), Eq. 27 
        ! Note: here cell i is grounded and cell i+1 is floating 
        ! Units: taud [Pa] 

        implicit none 

        real(wp), intent(OUT) :: taud
        real(wp), intent(IN) :: H_a, H_b          ! Ice thickness cell i and i+1, resp.
        real(wp), intent(IN) :: zb_a, zb_b        ! Bedrock elevation cell i and i+1, resp.
        real(wp), intent(IN) :: z_sl_a,  z_sl_b   ! Sea level cell i and i+1, resp.
        real(wp), intent(IN) :: dx  

        ! Local variables 
        real(wp) :: Ha, Hb, Sa, Sb, Ba, Bb, sla, slb 
        real(wp) :: dl, dx_ab 
        real(wp) :: H_mid, dzsdx 
        real(wp) :: rho_sw_ice, rho_ice_sw 
        integer  :: n 
        integer, parameter :: ntot = 100 

        ! Parameters 
        rho_sw_ice = rho_sw / rho_ice 
        rho_ice_sw = rho_ice / rho_sw 

        ! Get step size (dimensionless) and step resolution
        dl    = 1.0_wp / real(ntot-1.0_wp,wp)
        dx_ab  = dx*dl 

        ! Initialize driving stress to zero 
        taud = 0.0_wp 

        ! Step through the grid cell and calculate
        ! each piecewise value of driving stress
        do n = 1, ntot 

            Ha  = H_a + (H_b-H_a)*dl*(n-1)
            Hb  = H_a + (H_b-H_a)*dl*(n)

            Ba  = zb_a + (zb_b-zb_a)*dl*(n-1)
            Bb  = zb_a + (zb_b-zb_a)*dl*(n)
            
            sla = z_sl_a + (z_sl_b-z_sl_a)*dl*(n-1)
            slb = z_sl_a + (z_sl_b-z_sl_a)*dl*(n)
            
            if (Ha < rho_sw_ice*(sla-Ba)) then 
                Sa = (1.0-rho_ice_sw)*Ha
            else 
                Sa = Ba + Ha 
            end if 

            if (Hb < rho_sw_ice*(slb-Bb)) then 
                Sb = (1.0-rho_ice_sw)*Hb
            else 
                Sb = Bb + Hb 
            end if 
            
            H_mid = 0.5_wp * (Ha+Hb)
            dzsdx = (Sb-Sa) / dx_ab 

            taud  = taud + (H_mid * dzsdx)*dl 

        end do 

        ! Finally multiply with rho_ice*g 
        taud = rho_ice*g *taud

        return 

    end subroutine integrate_gl_driving_stress_linear
    
    subroutine set_inactive_margins(ux,uy,f_ice,boundaries)

        implicit none

        real(wp), intent(INOUT) :: ux(:,:) 
        real(wp), intent(INOUT) :: uy(:,:) 
        real(wp), intent(IN)    :: f_ice(:,:) 
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: im1, ip1, jm1, jp1

        nx = size(f_ice,1) 
        ny = size(f_ice,2) 

        ! Find partially-filled outer margins and set velocity to zero
        ! (this will also treat all other ice-free points too) 
        
        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

            if (f_ice(i,j) .lt. 1.0_wp .and. f_ice(ip1,j) .eq. 0.0_wp) then 
                ux(i,j) = 0.0_wp 
            end if
            if (f_ice(i,j) .eq. 0.0_wp .and. f_ice(ip1,j) .lt. 1.0_wp) then 
                ux(i,j) = 0.0_wp
            end if 

            if (f_ice(i,j) .lt. 1.0_wp .and. f_ice(i,jp1) .eq. 0.0_wp) then 
                uy(i,j) = 0.0_wp 
            end if
            if (f_ice(i,j) .eq. 0.0_wp .and. f_ice(i,jp1) .lt. 1.0_wp) then 
                uy(i,j) = 0.0_wp 
            end if

        end do
        end do

        return

    end subroutine set_inactive_margins

    subroutine calc_ice_flux(qq_acx,qq_acy,ux_bar,uy_bar,H_ice,dx,dy)
        ! Calculate the ice flux at a given point.
        ! Note: calculated on ac-nodes.
        ! qq      [m3 a-1] 
        ! ux,uy   [m a-1]
        ! H_ice   [m] 

        implicit none 

        real(wp), intent(OUT) :: qq_acx(:,:)     ! [m3 a-1] Ice flux (acx nodes)
        real(wp), intent(OUT) :: qq_acy(:,:)     ! [m3 a-1] Ice flux (acy nodes)
        real(wp), intent(IN)  :: ux_bar(:,:)     ! [m a-1]  Vertically averaged velocity (acx nodes)
        real(wp), intent(IN)  :: uy_bar(:,:)     ! [m a-1]  Vertically averaged velocity (acy nodes)
        real(wp), intent(IN)  :: H_ice(:,:)      ! [m]      Ice thickness, aa-nodes
        real(wp), intent(IN)  :: dx              ! [m]      Horizontal resolution, x-dir
        real(wp), intent(IN)  :: dy              ! [m]      Horizontal resolution, y-dir 

        ! Local variables 
        integer :: i, j, nx, ny 
        real(wp) :: area_ac 

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        ! Reset fluxes to zero 
        qq_acx = 0.0 
        qq_acy = 0.0 

        ! acx-nodes 
        do j = 1, ny 
        do i = 1, nx-1 
            area_ac     = (0.5_wp*(H_ice(i,j)+H_ice(i+1,j))) * dx 
            qq_acx(i,j) = area_ac*ux_bar(i,j)
        end do 
        end do 

        ! acy-nodes 
        do j = 1, ny-1 
        do i = 1, nx 
            area_ac     = (0.5_wp*(H_ice(i,j)+H_ice(i,j+1))) * dy 
            qq_acy(i,j) = area_ac*uy_bar(i,j)
        end do 
        end do 
        
        return 

    end subroutine calc_ice_flux

    elemental function calc_vel_ratio(uxy_base,uxy_srf) result(f_vbvs)

        implicit none 

        real(wp), intent(IN) :: uxy_base, uxy_srf  
        real(wp) :: f_vbvs 

        ! Calculate the basal to surface velocity ratio, f_vbvs
        if ( uxy_srf .gt. 0.0) then
            f_vbvs = min(1.0, uxy_base / uxy_srf) 
        else 
            ! No ice (or no velocity)
            f_vbvs = 1.0 
        end if 

        return 

    end function calc_vel_ratio

    elemental subroutine limit_vel(u,u_lim)
        ! Apply a velocity limit (for stability)

        implicit none 

        real(wp), intent(INOUT) :: u 
        real(wp), intent(IN)    :: u_lim

        real(wp), parameter :: tol = 1e-10
        
        u = min(u, u_lim)
        u = max(u,-u_lim)

        ! Also avoid underflow errors 
        if (abs(u) .lt. tol) u = 0.0 

        return 

    end subroutine limit_vel
    
    subroutine picard_calc_error(corr,ux,uy,ux_prev,uy_prev)
        ! Calculate the current error, ie, the 'correction vector'
        ! as defined by De Smedt et al. (2010):
        ! corr = U_now - U_prev.

        implicit none

        real(wp), intent(OUT) :: corr(:) 
        real(wp), intent(IN)  :: ux(:,:) 
        real(wp), intent(IN)  :: uy(:,:) 
        real(wp), intent(IN)  :: ux_prev(:,:)
        real(wp), intent(IN)  :: uy_prev(:,:) 
        
        ! Local variables 
        integer :: i, j, nx, ny, k  

        nx = size(ux,1)
        ny = size(ux,2) 

        ! Consistency check 
        if (size(corr,1) .ne. 2*nx*ny) then 
            write(*,*) "calc_convergence_angle:: Error: corr(N) must have N=2*nx*ny."
            stop 
        end if 

        k = 0

        do j = 1, ny 
        do i = 1, nx  

            k = k+1 
            corr(k) = ux(i,j) - ux_prev(i,j) 

            k = k+1 
            corr(k) = uy(i,j) - uy_prev(i,j) 

        end do 
        end do 

        return

    end subroutine picard_calc_error

    subroutine picard_calc_error_angle(theta,corr_nm1,corr_nm2)
        ! Calculate the current error, ie, the 'correction vector'
        ! as defined by De Smedt et al. (2010):
        ! corr = U_now - U_prev.

        implicit none

        real(wp), intent(OUT) :: theta
        real(wp), intent(IN)  :: corr_nm1(:) 
        real(wp), intent(IN)  :: corr_nm2(:) 
        
        
        ! Local variables   
        real(dp) :: val 

        real(dp), parameter :: tol = 1e-5 

        val = sum(corr_nm1*corr_nm2) / &
                (sqrt(sum(corr_nm1**2))*sqrt(sum(corr_nm2**2))+tol)

        theta = acos(val) 

        return

    end subroutine picard_calc_error_angle

    subroutine picard_calc_convergence_l1rel_matrix(err_x,err_y,ux,uy,ux_prev,uy_prev)

        implicit none 

        real(wp), intent(OUT) :: err_x(:,:)
        real(wp), intent(OUT) :: err_y(:,:)
        real(wp), intent(IN)  :: ux(:,:) 
        real(wp), intent(IN)  :: uy(:,:) 
        real(wp), intent(IN)  :: ux_prev(:,:) 
        real(wp), intent(IN)  :: uy_prev(:,:)  

        ! Local variables

        real(wp), parameter :: ssa_vel_tolerance = 1e-2   ! [m/a] only consider points with velocity above this tolerance limit
        real(wp), parameter :: tol = 1e-5 

        ! Error in x-direction
        where (abs(ux) .gt. ssa_vel_tolerance) 
            err_x = 2.0_wp * abs(ux - ux_prev) / abs(ux + ux_prev + tol)
        elsewhere 
            err_x = 0.0_wp
        end where 

        ! Error in y-direction 
        where (abs(uy) .gt. ssa_vel_tolerance) 
            err_y = 2.0_wp * abs(uy - uy_prev) / abs(uy + uy_prev + tol)
        elsewhere 
            err_y = 0.0_wp
        end where 

        return 

    end subroutine picard_calc_convergence_l1rel_matrix

    subroutine picard_calc_convergence_l2(is_converged,resid,ux,uy,ux_prev,uy_prev, &
                                                mask_acx,mask_acy,resid_tol,iter,iter_max,log)

        ! Calculate the current level of convergence using the 
        ! L2 relative error norm between the current and previous
        ! velocity solutions.

        ! Note for parameter norm_method defined below: 

        ! norm_method=1: as defined by De Smedt et al. (2010):
        ! conv = ||U_now - U_prev||/||U_prev||.

        ! norm_method=2: as defined by Gagliardini et al., GMD, 2013, Eq. 65:
        ! conv = 2*||U_now - U_prev||/||U_now+U_prev||.
        
        implicit none 

        logical,  intent(OUT) :: is_converged
        real(wp), intent(OUT) :: resid 
        real(wp), intent(IN) :: ux(:,:) 
        real(wp), intent(IN) :: uy(:,:) 
        real(wp), intent(IN) :: ux_prev(:,:) 
        real(wp), intent(IN) :: uy_prev(:,:)  
        logical,  intent(IN) :: mask_acx(:,:) 
        logical,  intent(IN) :: mask_acy(:,:) 
        real(wp), intent(IN) :: resid_tol 
        integer,  intent(IN) :: iter 
        integer,  intent(IN) :: iter_max 
        logical,  intent(IN) :: log 

        ! Local variables
        integer :: i, j, nx, ny, k  
        real(dp) :: res1, res2
        
        real(wp) :: ux_resid_max 
        real(wp) :: uy_resid_max 
        integer  :: nx_check, ny_check  
        character(len=1) :: converged_txt 

        real(dp), parameter :: du_reg  = 1e-10              ! [m/yr] Small regularization factor to avoid divide-by-zero
        real(wp), parameter :: vel_tol = 1e-5               ! [m/yr] only consider points with velocity above this tolerance limit
        integer,  parameter :: norm_method = 1              ! See note above.
        
        ! Count how many points should be checked for convergence
        nx_check = count(abs(ux).gt.vel_tol .and. mask_acx)
        ny_check = count(abs(uy).gt.vel_tol .and. mask_acy)

        if ( (nx_check+ny_check) .gt. 0 ) then

            select case(norm_method)

                case(1)
                    
                    res1 = sqrt( sum((ux-ux_prev)*(ux-ux_prev),mask=abs(ux).gt.vel_tol .and. mask_acx) &
                               + sum((uy-uy_prev)*(uy-uy_prev),mask=abs(uy).gt.vel_tol .and. mask_acy) )

                    res2 = sqrt( sum((ux_prev)*(ux_prev),mask=abs(ux).gt.vel_tol .and. mask_acx) &
                               + sum((uy_prev)*(uy_prev),mask=abs(uy).gt.vel_tol .and. mask_acy) )

                    resid = res1/(res2+du_reg)

                case(2)

                    res1 = sqrt( sum((ux-ux_prev)*(ux-ux_prev),mask=abs(ux).gt.vel_tol .and. mask_acx) &
                               + sum((uy-uy_prev)*(uy-uy_prev),mask=abs(uy).gt.vel_tol .and. mask_acy) )

                    res2 = sqrt( sum((ux+ux_prev)*(ux+ux_prev),mask=abs(ux).gt.vel_tol .and. mask_acx) &
                               + sum((uy+uy_prev)*(uy+uy_prev),mask=abs(uy).gt.vel_tol .and. mask_acy) )

                    resid = 2.0_wp*res1/(res2+du_reg)

            end select 

             

        else 
            ! No points available for comparison, set residual equal to zero 

            resid = 0.0_wp 

        end if

        ! Check for convergence
        if (resid .lt. resid_tol) then 
            is_converged = .TRUE. 
            converged_txt = "C"
        else if (iter .eq. iter_max) then 
            is_converged = .TRUE. 
            converged_txt = "X" 
        else 
            is_converged = .FALSE. 
            converged_txt = ""
        end if 

        ! Also calculate maximum error magnitude for perspective
        if (nx_check .gt. 0) then 
            ux_resid_max = maxval(abs(ux-ux_prev),mask=abs(ux).gt.vel_tol .and. mask_acx)
        else 
            ux_resid_max = 0.0 
        end if 

        if (ny_check .gt. 0) then 
            uy_resid_max = maxval(abs(uy-uy_prev),mask=abs(uy).gt.vel_tol .and. mask_acy)
        else 
            uy_resid_max = 0.0 
        end if 

        !if (log .and. is_converged) then
        if (log) then
            ! Write summary to log if desired and iterations have completed

            ! Write summary to log
            write(*,"(a,a2,i4,g12.4,a3,2i8,2g12.4)") &
                "ssa: ", trim(converged_txt), iter, resid, " | ", nx_check, ny_check, ux_resid_max, uy_resid_max 

        end if 


if (.TRUE.) then
        if (ux_resid_max .ge. 9999.0_wp .or. uy_resid_max .ge. 9999.0_wp) then 
            ! Strange case is occurring. Poor convergence with high error, investigate

            write(io_unit_err,*) "ssa: Error: strange case occurring."
            write(io_unit_err,"(a,a2,i4,g12.4,a3,2i8,2g12.4)") &
            "ssa: ", trim(converged_txt), iter, resid, " | ", nx_check, ny_check, ux_resid_max, uy_resid_max 

            write(io_unit_err,*) "Writing diagnostic file: ssa_check.nc."

            ! Initialize output file 
            call ssa_diagnostics_write_init("ssa_check.nc",nx=size(ux,1),ny=size(ux,2),time_init=real(iter,wp))

            ! Write file with dummy zero values for variables we don't have
            call ssa_diagnostics_write_step("ssa_check.nc",ux,uy,resid,ux*0.0_wp,ux*0.0_wp,ux*0.0_wp, &
                                    int(ux*0.0_wp),int(ux*0.0_wp),ux-ux_prev,uy-uy_prev,ux*0.0_wp,ux*0.0_wp,ux*0.0_wp,ux*0.0_wp, &
                                    ux*0.0_wp,ux*0.0_wp,ux*0.0_wp,ux*0.0_wp,ux*0.0_wp,ux*0.0_wp,ux_prev,uy_prev,time=real(iter,wp))

            stop 

        end if 
end if 
        
        return 

    end subroutine picard_calc_convergence_l2

    elemental subroutine picard_relax(ux,uy,ux_prev,uy_prev,rel)
        ! Relax velocity solution with previous iteration 

        implicit none 

        real(wp), intent(INOUT) :: ux
        real(wp), intent(INOUT) :: uy
        real(wp), intent(IN)    :: ux_prev
        real(wp), intent(IN)    :: uy_prev
        real(wp), intent(IN)    :: rel

        ! Apply relaxation 
        ux = ux_prev + rel*(ux-ux_prev)
        uy = uy_prev + rel*(uy-uy_prev)
        
        return 

    end subroutine picard_relax




    

end module velocity_general
