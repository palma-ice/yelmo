module velocity_general 
    ! This module contains general routines that are used by several solvers. 
    
    use yelmo_defs ,only  : sp, dp, wp, prec, tol_underflow, rho_ice, rho_sw, rho_w, g
    use yelmo_tools, only : stagger_aa_ab, stagger_aa_ab_ice, &
                    integrate_trapezoid1D_1D, integrate_trapezoid1D_pt, minmax

    implicit none 

    private 
    public :: calc_uz_3D
    public :: calc_driving_stress
    public :: calc_driving_stress_gl
    public :: calc_ice_flux
    public :: calc_vel_ratio
    public :: limit_vel

contains 

    subroutine calc_uz_3D(uz,ux,uy,H_ice,f_ice,z_bed,z_srf,smb,bmb,dHdt,dzsdt,zeta_aa,zeta_ac,dx,dy)
        ! Following algorithm outlined by the Glimmer ice sheet model:
        ! https://www.geos.ed.ac.uk/~mhagdorn/glide/glide-doc/glimmer_htmlse9.html#x17-660003.1.5

        ! Note: rate of bedrock uplift (dzbdt) no longer considered, since the rate is 
        ! very small and now z_bed is updated externally (ie, now assume dzbdt = 0.0 here)

        implicit none 

        real(prec), intent(OUT) :: uz(:,:,:)        ! nx,ny,nz_ac
        real(prec), intent(IN)  :: ux(:,:,:)        ! nx,ny,nz_aa
        real(prec), intent(IN)  :: uy(:,:,:)        ! nx,ny,nz_aa
        real(prec), intent(IN)  :: H_ice(:,:)
        real(prec), intent(IN)  :: f_ice(:,:)
        real(prec), intent(IN)  :: z_bed(:,:) 
        real(prec), intent(IN)  :: z_srf(:,:) 
        real(prec), intent(IN)  :: smb(:,:) 
        real(prec), intent(IN)  :: bmb(:,:) 
        real(prec), intent(IN)  :: dHdt(:,:) 
        real(prec), intent(IN)  :: dzsdt(:,:) 
        real(prec), intent(IN)  :: zeta_aa(:)    ! z-coordinate, aa-nodes 
        real(prec), intent(IN)  :: zeta_ac(:)    ! z-coordinate, ac-nodes  
        real(prec), intent(IN)  :: dx 
        real(prec), intent(IN)  :: dy

        ! Local variables 
        integer :: i, j, k, nx, ny, nz_aa, nz_ac
        integer :: im1, ip1, jm1, jp1   
        real(prec) :: H_ij
        !real(prec) :: dHdx_aa, dHdy_aa, dzsdx_aa, dzsdy_aa 
        real(prec) :: dzbdx_aa
        real(prec) :: dzbdy_aa
        real(prec) :: duxdx_aa
        real(prec) :: duydy_aa
        real(prec) :: ux_aa 
        real(prec) :: uy_aa 
        real(prec) :: uz_grid 
        real(prec) :: uz_srf 
        real(prec) :: corr 

        real(prec), parameter :: dzbdt = 0.0   ! For posterity, keep dzbdt variable, but set to zero 
        real(prec), parameter :: tol   = 1e-4 

        nx    = size(ux,1)
        ny    = size(ux,2)
        nz_aa = size(zeta_aa,1)
        nz_ac = size(zeta_ac,1) 

        ! Initialize vertical velocity to zero 
        uz = 0.0 

        ! Next, calculate velocity 

        !$omp parallel do 
        do j = 1, ny
        do i = 1, nx

            ! Define neighbor indices
            im1 = max(1,i-1)
            ip1 = min(nx,i+1)
            jm1 = max(1,j-1)
            jp1 = min(ny,j+1)
            
            if (f_ice(i,j) .eq. 1.0) then

                ! Get weighted ice thickness for stability
!                 H_ij = (4.0*H_ice(i,j) + 2.0*(H_ice(im1,j)+H_ice(ip1,j)+H_ice(i,jm1)+H_ice(i,jp1))) / 16.0 &
!                       + (H_ice(im1,jm1)+H_ice(ip1,jm1)+H_ice(ip1,jp1)+H_ice(im1,jp1)) / 16.0 

                H_ij = H_ice(i,j) 

                ! Get the centered bedrock gradient 
                dzbdx_aa = (z_bed(ip1,j)-z_bed(im1,j))/(2.0_prec*dx)
                dzbdy_aa = (z_bed(i,jp1)-z_bed(i,jm1))/(2.0_prec*dy)
                
                ! Get the centered horizontal velocity at the base
                ux_aa = 0.5_prec* (ux(im1,j,1) + ux(i,j,1))
                uy_aa = 0.5_prec* (uy(i,jm1,1) + uy(i,j,1))
                
!                 ! Get the centered surface gradient 
!                 dzsdx_aa = (z_srf(ip1,j)-z_srf(im1,j))/(2.0_prec*dx)
!                 dzsdy_aa = (z_srf(i,jp1)-z_srf(i,jm1))/(2.0_prec*dy)
                
!                 ! Get the centered ice thickness gradient 
!                 dHdx_aa = (H_ice(ip1,j)-H_ice(im1,j))/(2.0_prec*dx)
!                 dHdy_aa = (H_ice(i,jp1)-H_ice(i,jm1))/(2.0_prec*dy)
                
                ! Determine grid vertical velocity at the base due to sigma-coordinates 
                ! Glimmer, Eq. 3.35 
                ! ajr, 2020-01-27, untested:::
!                 uz_grid = dzsdt(i,j) + (ux_aa*dzsdx_aa + uy_aa*dzsdy_aa) &
!                             - ( (1.0_prec-zeta_ac(1))*dHdt(i,j) + ux_aa*dHdx_aa + uy_aa*dHdy_aa )
                uz_grid = 0.0_prec 

                ! ===================================================================
                ! Greve and Blatter (2009) style:

                ! Determine basal vertical velocity for this grid point 
                ! Following Eq. 5.31 of Greve and Blatter (2009)
                uz(i,j,1) = dzbdt + uz_grid + bmb(i,j) + ux_aa*dzbdx_aa + uy_aa*dzbdy_aa

                if (abs(uz(i,j,1)) .le. tol) uz(i,j,1) = 0.0_prec 
                
                ! Determine surface vertical velocity following kinematic boundary condition 
                ! Glimmer, Eq. 3.10 [or Folwer, Chpt 10, Eq. 10.8]
                !uz_srf = dzsdt(i,j) + ux_aa*dzsdx_aa + uy_aa*dzsdy_aa - smb(i,j) 
                
                ! Integrate upward to each point above base until surface is reached 
                do k = 2, nz_ac 

                    ! Greve and Blatter (2009), Eq. 5.72
                    ! Bueler and Brown  (2009), Eq. 4
                    duxdx_aa  = (ux(i,j,k-1)   - ux(im1,j,k-1)  )/dx
                    duydy_aa  = (uy(i,j,k-1)   - uy(i,jm1,k-1)  )/dy
                    
                    ! Testing wider stencil for stability (no effect so far)
!                     duxdx_aa  = 0.5*((ux(i,jp1,k) - ux(im1,jp1,k))/dx + (ux(i,jm1,k) - ux(im1,jm1,k))/dx)
!                     duydy_aa  = 0.5*((uy(ip1,j,k) - uy(i+1,jm1,k))/dy + (uy(im1,j,k) - uy(im1,jm1,k))/dy)

                    uz(i,j,k) = uz(i,j,k-1) & 
                        - H_ij*(zeta_ac(k)-zeta_ac(k-1))*(duxdx_aa+duydy_aa)

                    ! Apply correction to match kinematic boundary condition at surface 
                    !uz(i,j,k) = uz(i,j,k) - zeta_ac(k)*(uz(i,j,k)-uz_srf)

                    if (abs(uz(i,j,k)) .le. tol) uz(i,j,k) = 0.0_prec 
                    
                end do 
                
            else 
                ! No ice here, set vertical velocity equal to negative accum and bedrock change 

                do k = 1, nz_ac 
                    uz(i,j,k) = dzbdt - max(smb(i,j),0.0)
                    if (abs(uz(i,j,k)) .le. tol) uz(i,j,k) = 0.0_prec 
               end do 

            end if 

        end do 
        end do 
        !$omp end parallel do 

        ! Calculate and apply correction for sigma-coordinate stretching 
        call calc_advec_vertical_column_correction(uz,ux,uy,H_ice,z_srf,dHdt,dzsdt,zeta_ac,dx)
        
        return 

    end subroutine calc_uz_3D

    subroutine calc_advec_vertical_column_correction(uz,ux,uy,H_ice,z_srf,dHdt,dzsdt,zeta_ac,dx)
        ! Calculate the corrected vertical velocity, accounting for stretching of 
        ! the vertical axis between grid cells due to the use of sigma-coordinates. 

        ! Note: parameter max_corr may be necessary for very steep topography that violates 
        ! shallow-model assumptions. Imposing this limit ensures the model can continue. 
        
        implicit none 

        real(prec), intent(INOUT) :: uz(:,:,:)          ! nx,ny,nz_ac
        real(prec), intent(IN)    :: ux(:,:,:)          ! nx,ny,nz_aa
        real(prec), intent(IN)    :: uy(:,:,:)          ! nx,ny,nz_aa
        real(prec), intent(IN)    :: H_ice(:,:)         ! nx,ny 
        real(prec), intent(IN)    :: z_srf(:,:)         ! nx,ny 
        real(prec), intent(IN)    :: dHdt(:,:)          ! nx,ny 
        real(prec), intent(IN)    :: dzsdt(:,:)         ! nx,ny 
        real(prec), intent(IN)    :: zeta_ac(:)         ! nz_ac
        real(prec), intent(IN)    :: dx   

        ! Local variables 
        integer :: i, j, k, nx, ny, nz_ac
        integer :: im1, ip1, jm1, jp1 
        real(prec) :: ux_aa, uy_aa 
        real(prec) :: dx_inv, dx_inv2
        real(prec) :: c_x, c_y, c_t 
        real(prec) :: corr 
        real(prec) :: uz_corr                           ! [m/a] nz_ac 
        
        real(prec), parameter :: tol = 1e-4 
        real(prec), parameter :: max_corr = 1.0_prec    ! Maximum allowed deviation from original uz (eg 200%)

        nx    = size(H_ice,1)
        ny    = size(H_ice,2)
        nz_ac = size(zeta_ac,1) 

        ! Define some constants 
        dx_inv  = 1.0_prec / dx 
        dx_inv2 = 1.0_prec / (2.0_prec*dx)
        
        !$omp parallel do 
        do j = 1, ny 
        do i = 1, nx 

            ! Define neighbor indices
            im1 = max(1,i-1)
            ip1 = min(nx,i+1)
            jm1 = max(1,j-1)
            jp1 = min(ny,j+1)
            
            do k = 1, nz_ac 

                ! Estimate direction of current flow into cell (x and y), centered horizontally in grid point
                ! and averaged to staggered cell edges where uz is defined.
                if (k .eq. 1) then 
                    ux_aa = 0.5_prec*(ux(i,j,k)+ux(im1,j,k))
                    uy_aa = 0.5_prec*(uy(i,j,k)+uy(i,jm1,k))
                else if (k .eq. nz_ac) then 
                    ux_aa = 0.5_prec*(ux(i,j,k-1))
                    uy_aa = 0.5_prec*(uy(i,j,k-1))
                else 
                    ux_aa = 0.25_prec*(ux(i,j,k-1)+ux(im1,j,k-1) + ux(i,j,k)+ux(im1,j,k))
                    uy_aa = 0.25_prec*(uy(i,j,k-1)+uy(i,jm1,k-1) + uy(i,j,k)+uy(i,jm1,k))
                end if 

                if (abs(ux_aa) .lt. TOL_UNDERFLOW) ux_aa = 0.0_wp 
                if (abs(uy_aa) .lt. TOL_UNDERFLOW) uy_aa = 0.0_wp 

                ! Get horizontal scaling correction terms 
                c_x = (1.0_prec-zeta_ac(k))*(H_ice(ip1,j)-H_ice(im1,j))*dx_inv2 - (z_srf(ip1,j)-z_srf(im1,j))*dx_inv2
                c_y = (1.0_prec-zeta_ac(k))*(H_ice(i,jp1)-H_ice(i,jm1))*dx_inv2 - (z_srf(i,jp1)-z_srf(i,jm1))*dx_inv2
                
                ! Get grid velocity term 
                c_t = (1.0_prec-zeta_ac(k))*dHdt(i,j) - dzsdt(i,j) 
                
                ! Calculate total correction term, and limit it to within max_corr 
                corr = ux_aa*c_x + uy_aa*c_y + c_t  
                corr = sign(min(abs(corr),abs(max_corr*uz(i,j,k))),corr)

                ! Apply correction 
                uz_corr = uz(i,j,k) + corr 

                ! Limit new velocity to avoid underflow errors 
                if (abs(uz_corr) .le. tol) uz_corr = 0.0_prec 

                ! Set uz equal to new corrected uz 
                uz(i,j,k) = uz_corr  

            end do         

        end do 
        end do 
        !$omp end parallel do 

        return 

    end subroutine calc_advec_vertical_column_correction

    subroutine calc_driving_stress(taud_acx,taud_acy,H_ice,dzsdx,dzsdy,dx,taud_lim,boundaries)
        ! Calculate driving stress on staggered grid points
        ! Units: taud [Pa] == [kg m-1 s-2]
        ! taud = rho_ice*g*H_ice*dzs/dx

        ! Note: interpolation to Ab nodes no longer used here.

        ! Note: use parameter taud_lim to limit maximum allowed driving stress magnitude applied in the model.
        ! Should be an extreme value. eg, if dzdx = taud_lim / (rho*g*H), 
        ! then for taud_lim=5e5 and H=3000m, dzdx = 2e5 / (910*9.81*3000) = 0.02, 
        ! which is a rather steep slope for a shallow-ice model.

        implicit none 

        real(prec), intent(OUT) :: taud_acx(:,:)    ! [Pa]
        real(prec), intent(OUT) :: taud_acy(:,:)    ! [Pa]
        real(prec), intent(IN)  :: H_ice(:,:)       ! [m]
        real(prec), intent(IN)  :: dzsdx(:,:)       ! [--]
        real(prec), intent(IN)  :: dzsdy(:,:)       ! [--]
        real(prec), intent(IN)  :: dx               ! [m] 
        real(prec), intent(IN)  :: taud_lim         ! [Pa]
        character(len=*), intent(IN) :: boundaries  ! Boundary conditions to apply 

        ! Local variables 
        integer :: i, j, nx, ny 
        integer    :: im1, ip1, jm1, jp1 
        real(prec) :: dy, rhog 
        real(prec) :: H_mid
        
        nx = size(H_ice,1)
        ny = size(H_ice,2) 

        ! Define shortcut parameter 
        rhog = rho_ice * g 

        ! Assume grid resolution is symmetrical 
        dy = dx 

        do j = 1, ny 
        do i = 1, nx 

            im1 = max(1, i-1)
            ip1 = min(nx,i+1)
            
            jm1 = max(1, j-1)
            jp1 = min(ny,j+1)

            ! x-direction
            H_mid         = 0.5_prec*(H_ice(i,j)+H_ice(ip1,j)) 
            taud_acx(i,j) = rhog * H_mid * dzsdx(i,j) 

            ! y-direction 
            H_mid         = 0.5_prec*(H_ice(i,j)+H_ice(i,jp1))
            taud_acy(i,j) = rhog * H_mid * dzsdy(i,j) 

        end do
        end do 

        ! Apply limit 
        where(abs(taud_acx) .gt. taud_lim) taud_acx = sign(taud_lim,taud_acx)
        where(abs(taud_acy) .gt. taud_lim) taud_acy = sign(taud_lim,taud_acy)
        
        if (trim(boundaries) .eq. "periodic") then 

            taud_acx(1,:)    = taud_acx(nx-2,:) 
            taud_acx(nx-1,:) = taud_acx(2,:) 
            taud_acx(nx,:)   = taud_acx(3,:) 
            taud_acx(:,1)    = taud_acx(:,ny-1)
            taud_acx(:,ny)   = taud_acx(:,2) 

            taud_acy(1,:)    = taud_acy(nx-1,:) 
            taud_acy(nx,:)   = taud_acy(2,:) 
            taud_acy(:,1)    = taud_acy(:,ny-2)
            taud_acy(:,ny-1) = taud_acy(:,2) 
            taud_acy(:,ny)   = taud_acy(:,3)

        else if (trim(boundaries) .eq. "infinite") then 

            taud_acx(1,:)    = taud_acx(2,:) 
            taud_acx(nx-1,:) = taud_acx(nx-2,:) 
            taud_acx(nx,:)   = taud_acx(nx-2,:) 
            taud_acx(:,1)    = taud_acx(:,2)
            taud_acx(:,ny)   = taud_acx(:,ny-1) 

            taud_acy(1,:)    = taud_acy(2,:) 
            taud_acy(nx,:)   = taud_acy(nx-1,:) 
            taud_acy(:,1)    = taud_acy(:,2)
            taud_acy(:,ny-1) = taud_acy(:,ny-2) 
            taud_acy(:,ny)   = taud_acy(:,ny-2)

        end if 

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

        real(prec), intent(OUT) :: taud_acx(:,:)
        real(prec), intent(OUT) :: taud_acy(:,:) 
        real(prec), intent(IN)  :: H_ice(:,:)
        real(prec), intent(IN)  :: z_srf(:,:)
        real(prec), intent(IN)  :: z_bed(:,:)
        real(prec), intent(IN)  :: z_sl(:,:)
        real(prec), intent(IN)  :: H_grnd(:,:)
        real(prec), intent(IN)  :: f_grnd(:,:)
        real(prec), intent(IN)  :: f_grnd_acx(:,:)
        real(prec), intent(IN)  :: f_grnd_acy(:,:)
        real(prec), intent(IN)  :: dx 
        integer,    intent(IN)  :: method        ! Which driving stress calculation to use
        integer,    intent(IN)  :: beta_gl_stag  ! Method of grounding line staggering of beta 

        ! Local variables 
        integer :: i, j, nx, ny
        integer :: im1, ip1, jm1, jp1  
        real(prec) :: dy, rhog 
        real(prec) :: taud_grnd, taud_flt, taud_now 
        real(prec) :: H_mid, H_gl, z_gl, H_grnd_mid 
        real(prec) :: dzsdx, dzsdy
        real(prec) :: dzsdx_1, dzsdx_2
        real(prec) :: H_1, H_2  
        real(prec) :: taud_old, fac_gl   

        real(prec), parameter :: slope_max = 0.05   ! Very high limit == 0.05, low limit < 0.01 

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
                        H_gl    = 0.5_prec*(H_ice(i,j)+H_ice(i+1,j))

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
                        H_gl    = 0.5_prec*(H_ice(i,j)+H_ice(i,j+1))

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

                        H_grnd_mid = 0.5_prec*(H_grnd(i,j) + H_grnd(i+1,j))

                        if (H_grnd_mid .gt. 0.0) then 
                            ! Consider grounded 
                            dzsdx = (z_srf(i+1,j)-z_srf(i,j)) / dx 
                        else 
                            ! Consider floating 
                            dzsdx = (H_ice(i+1,j)-H_ice(i,j)) / dx
                        end if 
                        call minmax(dzsdx,slope_max)  

                        ! Get the ice thickness at the ac-node
                        H_gl    = 0.5_prec*(H_ice(i,j)+H_ice(i+1,j))

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

                        H_grnd_mid = 0.5_prec*(H_grnd(i,j) + H_grnd(i,j+1))

                        if (H_grnd_mid .gt. 0.0) then 
                            ! Consider grounded 
                            dzsdx = (z_srf(i,j+1)-z_srf(i,j)) / dx 
                        else 
                            ! Consider floating 
                            dzsdx = (H_ice(i,j+1)-H_ice(i,j)) / dx
                        end if 
                        call minmax(dzsdx,slope_max)  

                        ! Get the ice thickness at the ac-node
                        H_gl    = 0.5_prec*(H_ice(i,j)+H_ice(i,j+1))

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
    
    subroutine calc_ice_flux(qq_acx,qq_acy,ux_bar,uy_bar,H_ice,dx,dy)
        ! Calculate the ice flux at a given point.
        ! Note: calculated on ac-nodes.
        ! qq      [m3 a-1] 
        ! ux,uy   [m a-1]
        ! H_ice   [m] 

        implicit none 

        real(prec), intent(OUT) :: qq_acx(:,:)     ! [m3 a-1] Ice flux (acx nodes)
        real(prec), intent(OUT) :: qq_acy(:,:)     ! [m3 a-1] Ice flux (acy nodes)
        real(prec), intent(IN)  :: ux_bar(:,:)     ! [m a-1]  Vertically averaged velocity (acx nodes)
        real(prec), intent(IN)  :: uy_bar(:,:)     ! [m a-1]  Vertically averaged velocity (acy nodes)
        real(prec), intent(IN)  :: H_ice(:,:)      ! [m]      Ice thickness, aa-nodes
        real(prec), intent(IN)  :: dx              ! [m]      Horizontal resolution, x-dir
        real(prec), intent(IN)  :: dy              ! [m]      Horizontal resolution, y-dir 

        ! Local variables 
        integer :: i, j, nx, ny 
        real(prec) :: area_ac 

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        ! Reset fluxes to zero 
        qq_acx = 0.0 
        qq_acy = 0.0 

        ! acx-nodes 
        do j = 1, ny 
        do i = 1, nx-1 
            area_ac     = (0.5_prec*(H_ice(i,j)+H_ice(i+1,j))) * dx 
            qq_acx(i,j) = area_ac*ux_bar(i,j)
        end do 
        end do 

        ! acy-nodes 
        do j = 1, ny-1 
        do i = 1, nx 
            area_ac     = (0.5_prec*(H_ice(i,j)+H_ice(i,j+1))) * dy 
            qq_acy(i,j) = area_ac*uy_bar(i,j)
        end do 
        end do 
        
        return 

    end subroutine calc_ice_flux

    elemental function calc_vel_ratio(uxy_base,uxy_srf) result(f_vbvs)

        implicit none 

        real(prec), intent(IN) :: uxy_base, uxy_srf  
        real(prec) :: f_vbvs 

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

        real(prec), intent(INOUT) :: u 
        real(prec), intent(IN)    :: u_lim

        real(prec), parameter :: tol = 1e-10
        
        u = min(u, u_lim)
        u = max(u,-u_lim)

        ! Also avoid underflow errors 
        if (abs(u) .lt. tol) u = 0.0 

        return 

    end subroutine limit_vel
    
end module velocity_general
