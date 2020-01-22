module velocity_hybrid_pd12

    use nml 
    use yelmo_defs ,only  : sp, dp, prec, tol_underflow, rho_ice, rho_sw, rho_w, g
    use yelmo_tools, only : stagger_aa_ab, stagger_aa_ab_ice, &
                    integrate_trapezoid1D_1D, integrate_trapezoid1D_pt, minmax

    implicit none 

    private 
    public :: set_ssa_masks 
    public :: calc_vel_basal
    public :: calc_uz_3D
    public :: calc_vertical_integrated_3D_ice
    public :: calc_shear_reduction
    public :: calc_shear_3D
    public :: calc_ice_flux
    public :: calc_basal_stress
    public :: calc_driving_stress_ac 
    public :: calc_driving_stress_gl_ac
    public :: calc_visc_eff
    public :: calc_stress_eff_horizontal_squared
    public :: calc_vel_ratio

contains 

    subroutine set_ssa_masks(ssa_mask_acx,ssa_mask_acy,beta_acx,beta_acy,H_ice,f_grnd_acx,f_grnd_acy,beta_max,use_ssa)
        ! Define where ssa calculations should be performed
        ! Note: could be binary, but perhaps also distinguish 
        ! grounding line/zone to use this mask for later gl flux corrections
        ! mask = 0: no ssa calculated
        ! mask = 1: shelfy-stream ssa calculated 
        ! mask = 2: shelf ssa calculated 

        implicit none 
        
        integer,    intent(OUT) :: ssa_mask_acx(:,:) 
        integer,    intent(OUT) :: ssa_mask_acy(:,:)
        real(prec), intent(IN)  :: beta_acx(:,:)
        real(prec), intent(IN)  :: beta_acy(:,:)
        real(prec), intent(IN)  :: H_ice(:,:)
        real(prec), intent(IN)  :: f_grnd_acx(:,:)
        real(prec), intent(IN)  :: f_grnd_acy(:,:)
        real(prec), intent(IN)  :: beta_max
        logical,    intent(IN)  :: use_ssa       ! SSA is actually active now? 

        ! Local variables
        integer    :: i, j, nx, ny
        real(prec) :: H_acx, H_acy
        
        nx = size(H_ice,1)
        ny = size(H_ice,2)
        
        ! Initially no active ssa points
        ssa_mask_acx = 0
        ssa_mask_acy = 0
        
        if (use_ssa) then 

            ! x-direction
            do j = 1, ny
            do i = 1, nx-1

                if (H_ice(i,j) .gt. 0.0 .or. H_ice(i+1,j) .gt. 0.0) then 
                    ! Ice is present on ac-node
                    
                    if (f_grnd_acx(i,j) .gt. 0.0) then 
                        ! Grounded ice or grounding line (ie, shelfy-stream)
                        ssa_mask_acx(i,j) = 1
                    else 
                        ! Shelf ice 
                        ssa_mask_acx(i,j) = 2
                    end if 

                    ! Deactivate if dragging is to high and away from grounding line
                    if ( beta_acx(i,j) .ge. beta_max .and. f_grnd_acx(i,j) .eq. 1.0 ) ssa_mask_acx(i,j) = 0 
                    
                end if

            end do 
            end do

            ! y-direction
            do j = 1, ny-1
            do i = 1, nx

                if (H_ice(i,j) .gt. 0.0 .or. H_ice(i,j+1) .gt. 0.0) then 
                    ! Ice is present on ac-node
                    
                    if (f_grnd_acy(i,j) .gt. 0.0) then 
                        ! Grounded ice or grounding line (ie, shelfy-stream)
                        ssa_mask_acy(i,j) = 1
                    else 
                        ! Shelf ice 
                        ssa_mask_acy(i,j) = 2
                    end if 

                    ! Deactivate if dragging is to high and away from grounding line
                    if ( beta_acy(i,j) .ge. beta_max .and. f_grnd_acy(i,j) .eq. 1.0 ) ssa_mask_acy(i,j) = 0 
                    
                end if
                 
            end do 
            end do

            ! Final check on both masks to avoid isolated non-ssa points
            do j = 2, ny-1
            do i = 2, nx-1

                ! acx-nodes 
                if (  ssa_mask_acx(i,j) .eq. 0 .and. &
                    ssa_mask_acx(i+1,j) .gt. 0 .and. ssa_mask_acx(i-1,j) .gt. 0 .and.  &
                    ssa_mask_acx(i,j+1) .gt. 0 .and. ssa_mask_acx(i,j-1) .gt. 0 ) then 

                    if (f_grnd_acx(i,j) .gt. 0.0) then 
                        ! Grounded ice or grounding line (ie, shelfy-stream)
                        ssa_mask_acx(i,j) = 1
                    else 
                        ! Shelf ice 
                        ssa_mask_acx(i,j) = 2
                    end if 

                end if 

                ! acy-nodes 
                if (  ssa_mask_acy(i,j) .eq. 0 .and. &
                    ssa_mask_acy(i+1,j) .gt. 0 .and. ssa_mask_acy(i-1,j) .gt. 0 .and.  &
                    ssa_mask_acy(i,j+1) .gt. 0 .and. ssa_mask_acy(i,j-1) .gt. 0 ) then 

                    if (f_grnd_acy(i,j) .gt. 0.0) then 
                        ! Shelfy stream 
                        ssa_mask_acy(i,j) = 1
                    else 
                        ! Shelf 
                        ssa_mask_acy(i,j) = 2
                    end if 

                end if 

            end do 
            end do

        end if 

        return
        
    end subroutine set_ssa_masks 
    
    subroutine calc_vel_basal(ux_b,uy_b,ux_bar,uy_bar,ux_i,uy_i)
        ! Calculate basal sliding from the difference of
        ! vertical mean velocity and vertical mean internal shear velocity
        ! Pollard and de Conto (2012), after Eqs. 2a & 2b
        
        implicit none
        
        real(prec), intent(OUT) :: ux_b(:,:) 
        real(prec), intent(OUT) :: uy_b(:,:)
        real(prec), intent(IN)  :: ux_bar(:,:) 
        real(prec), intent(IN)  :: uy_bar(:,:)
        real(prec), intent(IN)  :: ux_i(:,:) 
        real(prec), intent(IN)  :: uy_i(:,:)
        
        real(prec), parameter :: du_tol = 0.0_prec

        where(abs(ux_bar-ux_i) .gt. du_tol)
            ux_b = sign(ux_bar-ux_i,ux_bar)
        elsewhere
            ux_b = 0.0_prec
        end where
        
        where(abs(uy_bar-uy_i) .gt. du_tol)
            uy_b = sign(uy_bar-uy_i,uy_bar)
        elsewhere
            uy_b = 0.0_prec
        end where
        
        return
        
    end subroutine calc_vel_basal

    subroutine calc_uz_3D(uz,ux,uy,H_ice,z_bed,smb,bmb,zeta_aa,zeta_ac,dx,dy)
        ! Following algorithm outlined by the Glimmer ice sheet model:
        ! https://www.geos.ed.ac.uk/~mhagdorn/glide/glide-doc/glimmer_htmlse9.html#x17-660003.1.5

        ! Note: rate of bedrock uplift (dzbdt) no longer considered, since the rate is 
        ! very small and now z_bed is updated externally (ie, now assume dzbdt = 0.0 here)

        implicit none 

        real(prec), intent(OUT) :: uz(:,:,:)        ! nx,ny,nz_ac
        real(prec), intent(IN)  :: ux(:,:,:)        ! nx,ny,nz_aa
        real(prec), intent(IN)  :: uy(:,:,:)        ! nx,ny,nz_aa
        real(prec), intent(IN)  :: H_ice(:,:)
        real(prec), intent(IN)  :: z_bed(:,:) 
        real(prec), intent(IN)  :: smb(:,:) 
        real(prec), intent(IN)  :: bmb(:,:) 
        real(prec), intent(IN)  :: zeta_aa(:)    ! z-coordinate, aa-nodes 
        real(prec), intent(IN)  :: zeta_ac(:)    ! z-coordinate, ac-nodes  
        real(prec), intent(IN)  :: dx 
        real(prec), intent(IN)  :: dy

        ! Local variables 
        integer :: i, j, k, nx, ny, nz_aa, nz_ac   
        real(prec) :: H_ij 
        real(prec) :: dzbdx_ac
        real(prec) :: dzbdy_ac
        real(prec) :: duxdx_aa
        real(prec) :: duydy_aa

        real(prec), parameter :: dzbdt = 0.0   ! For posterity, keep dzbdt variable, but set to zero 

        nx    = size(ux,1)
        ny    = size(ux,2)
        nz_aa = size(zeta_aa,1)
        nz_ac = size(zeta_ac,1) 

        ! Initialize vertical velocity to zero 
        uz = 0.0 

        ! Next, calculate velocity 
        do j = 2, ny-1
        do i = 2, nx-1

            if (H_ice(i,j) .gt. 0.0) then

                ! Get weighted ice thickness for stability
!                 H_ij = (4.0*H_ice(i,j) + 2.0*(H_ice(i-1,j)+H_ice(i+1,j)+H_ice(i,j-1)+H_ice(i,j+1))) / 16.0 &
!                       + (H_ice(i-1,j-1)+H_ice(i+1,j-1)+H_ice(i+1,j+1)+H_ice(i-1,j+1)) / 16.0 

                H_ij = H_ice(i,j) 

                ! Get the staggered bedrock gradient 
                dzbdx_Ac = (z_bed(i+1,j)-z_bed(i,j))/dx
                dzbdy_Ac = (z_bed(i,j+1)-z_bed(i,j))/dy
                
                ! ===================================================================
                ! Greve and Blatter (2009) style:

                ! Determine basal vertical velocity for this grid point 
                ! Following Eq. 5.31 of Greve and Blatter (2009)
                uz(i,j,1) = dzbdt + bmb(i,j) + ux(i,j,1)*dzbdx_Ac + uy(i,j,1)*dzbdy_Ac

                ! Integrate upward to each point above base until surface is reached 
                do k = 2, nz_ac 

                    ! Greve and Blatter (2009), Eq. 5.72
                    ! Bueler and Brown  (2009), Eq. 4
                    duxdx_aa  = (ux(i,j,k)   - ux(i-1,j,k)  )/dx
                    duydy_aa  = (uy(i,j,k)   - uy(i,j-1,k)  )/dy
                    
                    ! Testing wider stencil for stability (no effect so far)
!                     duxdx_aa  = 0.5*((ux(i,j+1,k) - ux(i-1,j+1,k))/dx + (ux(i,j-1,k) - ux(i-1,j-1,k))/dx)
!                     duydy_aa  = 0.5*((uy(i+1,j,k)   - uy(i+1,j-1,k)  )/dy + (uy(i-1,j,k)   - uy(i-1,j-1,k)  )/dy)

                    uz(i,j,k) = uz(i,j,k-1) & 
                        - H_ij*(zeta_ac(k)-zeta_ac(k-1))*(duxdx_aa+duydy_aa)

                end do 
                
            else 
                ! No ice here, set vertical velocity equal to negative accum and bedrock change 

                uz(i,j,:) = dzbdt - max(smb(i,j),0.0)

            end if 

        end do 
        end do 

        return 

    end subroutine calc_uz_3D 

    function calc_vertical_integrated_3D_ice(var,H_ice,sigma) result(var_int)
        ! Vertically integrate a field 3D field (nx,ny,nz)
        ! layer by layer (in the z-direction), return a 3D array
        
        implicit none

        real(prec), intent(IN) :: var(:,:,:)
        real(prec), intent(IN) :: H_ice(:,:) 
        real(prec), intent(IN) :: sigma(:)
        real(prec) :: var_int(size(var,1),size(var,2),size(var,3))

        ! Local variables 
        integer :: i, j, nx, ny

        nx = size(var,1)
        ny = size(var,2)

        do j = 1, ny
        do i = 1, nx
            var_int(i,j,:) = integrate_trapezoid1D_1D(var(i,j,:),sigma*H_ice(i,j))
        end do
        end do

        return

    end function calc_vertical_integrated_3D_ice

    subroutine calc_shear_reduction(lhs_x,lhs_y,ux_b,uy_b,visc_eff,dx)
        ! Calculate reduction in driving stress expected from
        ! basal sliding
        ! Pollard and de Conto (2012), obtained from Eqs. 2a/2b
        
        implicit none
        
        real(prec), intent(OUT) :: lhs_x(:,:)    ! acx node
        real(prec), intent(OUT) :: lhs_y(:,:)    ! acy node
        real(prec), intent(IN)  :: ux_b(:,:)     ! acx node
        real(prec), intent(IN)  :: uy_b(:,:)     ! acx node 
        real(prec), intent(IN)  :: visc_eff(:,:) ! aa node 
        real(prec), intent(IN)  :: dx 

        ! Local variables 
        integer :: i, j, nx, ny 
        real(prec) :: dy 
        real(prec) :: duxdx, duydy, duydx, duxdy 
        real(prec), allocatable :: varx1_aa(:,:), vary1_aa(:,:) 
        real(prec), allocatable :: varx2_aa(:,:), vary2_aa(:,:) 

        ! Assume dx and dy are the same 
        dy = dx 

        ! Get matrix sizes 
        nx = size(lhs_x,1)
        ny = size(lhs_x,2) 

        allocate(varx1_aa(nx,ny))
        allocate(vary1_aa(nx,ny)) 
        allocate(varx2_aa(nx,ny))
        allocate(vary2_aa(nx,ny)) 

        ! Initially set lhs equal to zero
        lhs_x = 0.0
        lhs_y = 0.0
        
        varx1_aa = 0.0 
        vary1_aa = 0.0 
        varx2_aa = 0.0 
        vary2_aa = 0.0 

        ! First get arrays of intermediate values
        do j = 2, ny-1 
        do i = 2, nx-1 

            ! Calculate intermediate terms of
            ! lhs_x: 2*visc*(2*dudx+dvdy), visc*(dudy+dvdx)
            ! lhs_y: 2*visc*(2*dvdy+dudx), visc*(dudy+dvdx)
            ! on aa nodes 

            ! Gradients on aa nodes 
            duxdx = (ux_b(i,j)-ux_b(i-1,j))/dx 
            duydy = (uy_b(i,j)-uy_b(i,j-1))/dy 

            ! Cross terms on aa nodes 
            duxdy  = ( 0.25_prec*(ux_b(i,j)+ux_b(i-1,j)+ux_b(i,j+1)+ux_b(i-1,j+1)) &
                      -0.25_prec*(ux_b(i,j)+ux_b(i-1,j)+ux_b(i,j-1)+ux_b(i-1,j-1))) /dy 

            duydx  = ( 0.25_prec*(uy_b(i,j)+uy_b(i,j-1)+uy_b(i+1,j)+uy_b(i+1,j-1)) &
                      -0.25_prec*(uy_b(i,j)+uy_b(i,j-1)+uy_b(i-1,j)+uy_b(i-1,j-1))) /dx  

            ! Intermediate terms on aa nodes - x-direction  
            varx1_aa(i,j) = 2.0_prec*visc_eff(i,j)*(2.0_prec*duxdx+duydy)
            vary1_aa(i,j) = visc_eff(i,j)*(duxdy+duydx)

            ! Intermediate terms on aa nodes - y-direction  
            vary2_aa(i,j) = 2.0_prec*visc_eff(i,j)*(2.0_prec*duydy+duxdx)
            varx2_aa(i,j) = visc_eff(i,j)*(duxdy+duydx)

        end do 
        end do 

        ! Next calculate the total lhs_x term, on acx nodes 
        do j = 2, ny-1 
        do i = 1, nx-1 

            ! On acx nodes 
            duxdx = (varx1_aa(i+1,j)-varx1_aa(i,j))/dx 
            duxdy = ( 0.25_prec*(vary1_aa(i,j)+vary1_aa(i+1,j)+vary1_aa(i,j+1)+vary1_aa(i+1,j+1)) &
                     -0.25_prec*(vary1_aa(i,j)+vary1_aa(i+1,j)+vary1_aa(i,j-1)+vary1_aa(i+1,j-1))) /dy

            lhs_x(i,j) = duxdx + duydy 

        end do 
        end do 

        ! Next calculate the total lhs_y term, on acy nodes 
        do j = 1, ny-1 
        do i = 2, nx-1 

            ! On acy nodes 
            duydy = (vary2_aa(i,j+1)-vary2_aa(i,j))/dy 

            duydx = ( 0.25_prec*(varx2_aa(i,j)+varx2_aa(i+1,j)+varx2_aa(i,j+1)+varx2_aa(i+1,j+1)) &
                     -0.25_prec*(varx2_aa(i,j)+varx2_aa(i,j+1)+varx2_aa(i-1,j)+varx2_aa(i-1,j+1))) /dx

            lhs_y(i,j) = duydy + duydx 

        end do 
        end do 

        return
        
    end subroutine calc_shear_reduction
    
    subroutine calc_shear_3D(duxdz,duydz,dd_ab,taud_acx,taud_acy,ATT,lhs_x,lhs_y, &
                             sigma_horiz_sq,zeta_aa,n_glen,boundaries)
        ! Calculate internal shear via SIA minus stretching
        ! Pollard and de Conto (2012), Eq. 1
        ! Calculate coefficients on Ab nodes, then
        ! stagger to ac nodes for each component 
        
        implicit none
        
        real(prec), intent(OUT) :: duxdz(:,:,:)         ! nx,ny,nz_aa [1/a]
        real(prec), intent(OUT) :: duydz(:,:,:)         ! nx,ny,nz_aa [1/a]
        real(prec), intent(OUT) :: dd_ab(:,:,:)         ! nx,ny,nz_aa [m2/a]
        real(prec), intent(IN)  :: taud_acx(:,:)        ! dzsdx(:,:)
        real(prec), intent(IN)  :: taud_acy(:,:)        ! dzsdy(:,:) 
        real(prec), intent(IN)  :: ATT(:,:,:)           ! nx,ny,nz_aa 
        real(prec), intent(IN)  :: lhs_x(:,:)           ! ac-nodes
        real(prec), intent(IN)  :: lhs_y(:,:)           ! ac-nodes
        real(prec), intent(IN)  :: sigma_horiz_sq(:,:)  ! [Pa], aa-nodes
        real(prec), intent(IN)  :: zeta_aa(:)           ! Vertical axis
        real(prec), intent(IN)  :: n_glen               ! Glen law exponent
        character(len=*), intent(IN) :: boundaries      ! Boundary conditions to apply 

        ! Local variables
        integer    :: i, j, k, nx, ny, nz_aa 
        real(prec) :: exp1, depth    
        real(prec), allocatable :: sigma_ab(:,:)
        real(prec), allocatable :: sigma_xz_ab(:,:), sigma_yz_ab(:,:)   ! [Pa]
        real(prec), allocatable :: sigma_xz(:,:), sigma_yz(:,:)         ! [Pa]
        real(prec), allocatable :: sigma_horiz_sq_ab(:,:)               ! [Pa]
        real(prec), allocatable :: sigma_tot_sq_ab(:,:)                 ! [Pa]
        real(prec) :: ATT_ab, ddx, ddy 
        real(prec) :: ATT_ac, sigma_tot_sq_ac
        logical :: is_mismip 

        is_mismip = .FALSE. 
        if (trim(boundaries) .eq. "MISMIP3D") is_mismip = .TRUE. 

        ! Define exponent for later use 
        exp1 = (n_glen-1.0_prec)/2.0_prec 

        ! Get array dimensions
        nx    = size(duxdz,1)
        ny    = size(duxdz,2)
        nz_aa = size(zeta_aa,1)

        ! Allocate local variables  
        allocate(sigma_ab(nx,ny)) 
        allocate(sigma_xz_ab(nx,ny)) 
        allocate(sigma_yz_ab(nx,ny)) 
        allocate(sigma_xz(nx,ny)) 
        allocate(sigma_yz(nx,ny))
        allocate(sigma_horiz_sq_ab(nx,ny))
        allocate(sigma_tot_sq_ab(nx,ny))

        ! Initially set shear to zero everywhere 
        ! (and diagnostic diffusion variable)
        duxdz   = 0.0
        duydz   = 0.0
        dd_ab   = 0.0  

        ! Stagger horizontal stress contribution to ab-nodes 
        sigma_horiz_sq_ab = stagger_aa_ab(sigma_horiz_sq)

        ! Calculate shear for each layer of the ice sheet 
        ! with diffusivity on Ab nodes (ie, EISMINT type I)
        do k = 1, nz_aa 

            ! Determine current depth fraction
            depth = (1.0_prec-zeta_aa(k))

            ! Determine the vertical shear for current layer (Ac nodes)
            ! Pollard and de Conto (2012), Eq. 3 
            ! ajr: note, initial negative sign moved to be explicit
            ! in the final calculation of dux/dz duy/dz, so that 
            ! the definition of diffusivity dd_ab is consistent with
            ! other SIA derivations  
!             sigma_xz = (rhog*H_ice_acx*dzsdx - lhs_x)*depth
!             sigma_yz = (rhog*H_ice_acy*dzsdy - lhs_y)*depth
            sigma_xz = (taud_acx - lhs_x)*depth
            sigma_yz = (taud_acy - lhs_y)*depth
            
            ! Acx => Ab
            sigma_xz_ab = 0.0 
            do j = 1, ny-1
            do i = 1, nx 
                sigma_xz_ab(i,j) = 0.5_prec*(sigma_xz(i,j)+sigma_xz(i,j+1))
            end do 
            end do 

            ! Acy => Ab
            sigma_yz_ab = 0.0
            do j = 1, ny
            do i = 1, nx-1 
                sigma_yz_ab(i,j) = 0.5_prec*(sigma_yz(i,j)+sigma_yz(i+1,j))
            end do 
            end do 

            ! Get total sigma terms squared (input to Eq. 1a/1b) on ab-nodes
            ! with stress from horizontal stretching `sigma_horiz_sq`
            sigma_tot_sq_ab = sigma_xz_ab**2 + sigma_yz_ab**2 + sigma_horiz_sq_ab
            
            ! Determine the diffusivity for current layer on Ab nodes
            ! Pollard and de Conto (2012), Eq. 1a/dd 
            do j = 1, ny-1
            do i = 1, nx-1 

                ! Stagger rate factor: Aa => Ab nodes
                ATT_ab       = 0.25_prec * (ATT(i,j,k)+ATT(i+1,j,k)+ATT(i,j+1,k)+ATT(i+1,j+1,k))

                ! Calculate diffusivity on Ab nodes
                dd_ab(i,j,k) = 2.0_prec*ATT_ab*(sigma_tot_sq_ab(i,j)**exp1)
            end do 
            end do

            ! Calculate the shear for current layer on Ac nodes
            ! Pollard and de Conto (2012), Eq. 1a/1b
            ! Ac-x nodes
            do j = 2, ny
            do i = 1, nx 

                ! Get sigma_tot on ac-nodes 
                !sigma_tot_sq_ac = 0.5_prec * (sigma_tot_sq_ab(i,j-1)+sigma_tot_sq_ab(i,j))
                !ATT_ac          = 0.5_prec * (ATT(i,j,k) + ATT(i+1,j,k))
                !ddx             = 2.0_prec * ATT_ac * (sigma_tot_sq_ac**exp1)

                ! Get coefficient on Ac nodes
                ddx = 0.5_prec*(dd_ab(i,j,k)+dd_ab(i,j-1,k))

                ! Calculate shear on ac nodes 
                duxdz(i,j,k) = -ddx*sigma_xz(i,j)

            end do 
            end do

            ! Ac-y nodes
            do j = 1, ny
            do i = 2, nx 

                ! Get sigma_tot on ac-nodes 
                !sigma_tot_sq_ac = 0.5_prec * (sigma_tot_sq_ab(i-1,j)+sigma_tot_sq_ab(i,j))
                !ATT_ac          = 0.5_prec * (ATT(i,j,k) + ATT(i,j+1,k))
                !ddy             = 2.0_prec * ATT_ac * (sigma_tot_sq_ac**exp1)
                
                ! Get coefficient on Ac nodes
                ddy = 0.5_prec*(dd_ab(i,j,k)+dd_ab(i-1,j,k))

                ! Calculate shear on ac nodes 
                duydz(i,j,k) = -ddy*sigma_yz(i,j)

            end do 
            end do
            
        end do 
        
        if (is_mismip) then 
            ! Fill in boundaries 
            ! ajr, 2018-11-12: mismip is still showing slight assymmetry in the y-direction when 
            ! shear is included in solution (ie, mix_method=1). Could be related to boundary conditions here.

            duxdz(:,1,:)  = duxdz(:,2,:) 
            duxdz(:,ny,:) = duxdz(:,ny-1,:) 
            duydz(:,1,:)  = duydz(:,2,:) 
            duydz(:,ny,:) = duydz(:,ny-1,:) 
            
            duxdz(1,:,:)  = duxdz(2,:,:) 
            duxdz(nx,:,:) = duxdz(nx-1,:,:) 
            duydz(1,:,:)  = duydz(2,:,:) 
            duydz(nx,:,:) = duydz(nx-1,:,:) 
            
        end if 

        return
        
    end subroutine calc_shear_3D
    
    function calc_visc_eff(ux,uy,duxdz,duydz,H_ice,ATT,zeta_aa,dx,dy,n) result(visc)
        ! Calculate effective viscosity eta to be used in SSA solver
        ! Pollard and de Conto (2012), Eqs. 2a/b and Eq. 4 (`visc=mu*H_ice*A**(-1/n)`)
        ! Note: calculated on same nodes as eps_sq (aa-nodes by default)
        ! Note: this is equivalent to the vertically-integrated viscosity, 
        ! since it is multiplied with H_ice 

        implicit none 
        
        real(prec), intent(IN) :: ux(:,:)      ! Vertically averaged horizontal velocity, x-component
        real(prec), intent(IN) :: uy(:,:)      ! Vertically averaged horizontal velocity, y-component
        real(prec), intent(IN) :: duxdz(:,:)   ! Only from internal shear, vertical mean (ac-nodes)
        real(prec), intent(IN) :: duydz(:,:)   ! Only from internal shear, vertical mean (ac-nodes)
        real(prec), intent(IN) :: H_ice(:,:)   ! Ice thickness
        real(prec), intent(IN) :: ATT(:,:,:)   ! nx,ny,nz_aa Rate factor
        real(prec), intent(IN) :: zeta_aa(:)  ! Vertical axis (sigma-coordinates from 0 to 1)
        real(prec), intent(IN) :: dx, dy
        real(prec), intent(IN) :: n
        real(prec) :: visc(size(ux,1),size(ux,2))   ! [Pa a m]
        
        ! Local variables 
        integer :: i, j, nx, ny, k, nz_aa 
        integer :: im1, ip1, jm1, jp1  
        real(prec) :: inv_4dx, inv_4dy
        real(prec) :: dudx, dudy
        real(prec) :: dvdx, dvdy 
        real(prec) :: duxdz_aa, duydz_aa
        real(prec) :: eps_sq, mu  

        real(prec), allocatable :: visc1D(:) 
        
        real(prec), parameter :: visc_min     = 1e3_prec 
        real(prec), parameter :: epsilon_sq_0 = 1e-6_prec   ! [a^-1] Bueler and Brown (2009), Eq. 26
        
        nx    = size(ux,1)
        ny    = size(ux,2)
        nz_aa = size(zeta_aa,1) 

        allocate(visc1D(nz_aa))

        inv_4dx = 1.0_prec / (4.0_prec*dx) 
        inv_4dy = 1.0_prec / (4.0_prec*dy) 

        ! Initialize viscosity to minimum value 
        visc = visc_min 
        
        ! Loop over domain to calculate viscosity at each aa-node

        do j = 1, ny
        do i = 1, nx

            im1 = max(i-1,1) 
            ip1 = min(i+1,nx) 
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny) 
            
            ! 1. Calculate effective strain components from horizontal stretching

            ! Aa node
            dudx = (ux(i,j) - ux(im1,j))/dx
            dvdy = (uy(i,j) - uy(i,jm1))/dy

            ! Calculation of cross terms on central aa-nodes (symmetrical results)
            dudy = ((ux(i,jp1)   - ux(i,jm1))    &
                  + (ux(im1,jp1) - ux(im1,jm1))) * inv_4dx 
            dvdx = ((uy(ip1,j)   - uy(im1,j))    &
                  + (uy(ip1,jm1) - uy(im1,jm1))) * inv_4dy 

            ! 2. Un-stagger shear terms to central aa-nodes

            duxdz_aa = 0.5_prec*(duxdz(i,j) + duxdz(im1,j))
            duydz_aa = 0.5_prec*(duydz(i,j) + duydz(im1,j))
            
            ! Avoid underflows 
            if (abs(dudx) .lt. tol_underflow) dudx = 0.0 
            if (abs(dvdy) .lt. tol_underflow) dvdy = 0.0 
            if (abs(dudy) .lt. tol_underflow) dudy = 0.0 
            if (abs(dvdx) .lt. tol_underflow) dvdx = 0.0 

            if (abs(duxdz_aa) .lt. tol_underflow) duxdz_aa = 0.0 
            if (abs(duydz_aa) .lt. tol_underflow) duydz_aa = 0.0 
            
            ! 3. Calculate the total effective strain rate
            ! from Pollard and de Conto (2012), Eq. 6
            ! (Note: equation in text seems to have typo concerning cross terms)

            eps_sq = dudx**2 + dvdy**2 + dudx*dvdy + 0.25*(dudy+dvdx)**2 &
                          + 0.25_prec*duxdz_aa**2 + 0.25_prec*duydz_aa**2 &
                          + epsilon_sq_0
            
            ! 4. Calculate the effective visocity (`eta` in Greve and Blatter, 2009)
            ! Pollard and de Conto (2012), Eqs. 2a/b and Eq. 4 (`visc=A**(-1/n)*mu*H_ice`)

            mu = 0.5_prec*(eps_sq)**((1.0_prec - n)/(2.0_prec*n))

            do k = 1, nz_aa 
                visc1D(k) = ATT(i,j,k)**(-1.0_prec/n) * mu
                
            end do 

            visc(i,j) = integrate_trapezoid1D_pt(visc1D,zeta_aa) 
            if (H_ice(i,j) .gt. 1.0_prec) visc(i,j) = visc(i,j) * H_ice(i,j) 
            
        end do 
        end do  

        return
        
    end function calc_visc_eff

    function calc_stress_eff_horizontal_squared(ux,uy,ATT_bar,dx,dy,n) result(sigma_sq)
        ! Calculate squared effective stress of horizontal stretching terms
        ! (from vertically averaged quantities)
        ! as defined in Pollard and de Conto (2012), Eq. 7
        ! Note: calculated on aa-nodes 

        implicit none 
        
        real(prec), intent(IN) :: ux(:,:)      ! Vertically averaged horizontal velocity, x-component
        real(prec), intent(IN) :: uy(:,:)      ! Vertically averaged horizontal velocity, y-component
        real(prec), intent(IN) :: ATT_bar(:,:)   ! nx,ny,nz_aa Rate factor
        real(prec), intent(IN) :: dx, dy
        real(prec), intent(IN) :: n
        real(prec) :: sigma_sq(size(ux,1),size(ux,2))
        
        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: im1, ip1, jm1, jp1  
        real(prec) :: inv_4dx, inv_4dy
        real(prec) :: dudx, dudy
        real(prec) :: dvdx, dvdy 
        real(prec) :: duxdz_aa, duydz_aa
        real(prec) :: eps_sq, mu  

        real(prec), parameter :: epsilon_sq_0 = 1e-6_prec   ! [a^-1] Bueler and Brown (2009), Eq. 26
        
        nx    = size(ux,1)
        ny    = size(ux,2)

        inv_4dx = 1.0_prec / (4.0_prec*dx) 
        inv_4dy = 1.0_prec / (4.0_prec*dy) 

        ! Initialize viscosity to zero
        sigma_sq = 0.0 
        
        ! Loop over domain to calculate viscosity at each Aa node

        do j = 1, ny
        do i = 1, nx

            im1 = max(i-1,1) 
            ip1 = min(i+1,nx) 
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny) 
            
            ! 1. Calculate effective strain components from horizontal stretching

            ! Aa node
            dudx = (ux(i,j) - ux(im1,j))/dx
            dvdy = (uy(i,j) - uy(i,jm1))/dy

            ! Calculation of cross terms on central Aa nodes (symmetrical results)
            dudy = ((ux(i,jp1)   - ux(i,jm1))    &
                  + (ux(im1,jp1) - ux(im1,jm1))) * inv_4dx 
            dvdx = ((uy(ip1,j)   - uy(im1,j))    &
                  + (uy(ip1,jm1) - uy(im1,jm1))) * inv_4dy 

            ! 2. Calculate the total effective strain rate
            ! from Pollard and de Conto (2012), Eq. 6
            ! (Note: equation in text seems to have typo concerning cross terms)

            ! Avoid underflows 
            if (abs(dudx) .lt. tol_underflow) dudx = 0.0 
            if (abs(dvdy) .lt. tol_underflow) dvdy = 0.0 
            if (abs(dudy) .lt. tol_underflow) dudy = 0.0 
            if (abs(dvdx) .lt. tol_underflow) dvdx = 0.0 

            eps_sq = dudx**2 + dvdy**2 + dudx*dvdy + 0.25*(dudy+dvdx)**2 + epsilon_sq_0
            
            ! 3. Calculate the effective visocity (`eta` in Greve and Blatter, 2009)
            ! Pollard and de Conto (2012), Eqs. 2a/b and Eq. 4 (`visc=A**(-1/n)*mu*H_ice`)

            mu = 0.5_prec*(eps_sq)**((1.0_prec - n)/(2.0_prec*n))

            ! 4. Pollard and de Conto (2012), Eq. 7 on aa-nodes:
            ! Note: equation in text seems to have typo concerning cross terms
            sigma_sq(i,j) = (2.0_prec*mu * ATT_bar(i,j)**(-1.0_prec/n) )**2 * eps_sq

        end do 
        end do  

        return
        
    end function calc_stress_eff_horizontal_squared

    subroutine calc_basal_stress(taub_acx,taub_acy,beta_acx,beta_acy,ux_b,uy_b)
        ! Calculate the basal stress resulting from sliding (friction times velocity)
        ! Note: calculated on ac-nodes.
        ! taub [Pa] 
        ! beta [Pa a m-1]
        ! u    [m a-1]
        ! taub = -beta*u 

        implicit none 

        real(prec), intent(OUT) :: taub_acx(:,:)   ! [Pa] Basal stress (acx nodes)
        real(prec), intent(OUT) :: taub_acy(:,:)   ! [Pa] Basal stress (acy nodes)
        real(prec), intent(IN)  :: beta_acx(:,:)   ! [Pa a m-1] Basal friction (acx nodes)
        real(prec), intent(IN)  :: beta_acy(:,:)   ! [Pa a m-1] Basal friction (acy nodes)
        real(prec), intent(IN)  :: ux_b(:,:)       ! [m a-1] Basal velocity (acx nodes)
        real(prec), intent(IN)  :: uy_b(:,:)       ! [m a-1] Basal velocity (acy nodes)
        
        real(prec), parameter :: tol = 1e-3_prec 

        ! Calculate basal stress 
        taub_acx = -beta_acx * ux_b 
        taub_acy = -beta_acy * uy_b 

        ! Avoid underflows
        where(abs(taub_acx) .lt. tol) taub_acx = 0.0_prec 
        where(abs(taub_acy) .lt. tol) taub_acy = 0.0_prec 
        
        return 

    end subroutine calc_basal_stress

    subroutine calc_driving_stress_ac(taud_acx,taud_acy,H_ice,dzsdx,dzsdy,dx)
        ! taud = rho_ice*g*H_ice
        ! Calculate driving stress on staggered grid points, with 
        ! special treatment of the grounding line 
        ! Units: taud [Pa] == [kg m-1 s-2]
        
        ! Note: interpolation to Ab nodes no longer used here.

        implicit none 

        real(prec), intent(OUT) :: taud_acx(:,:)
        real(prec), intent(OUT) :: taud_acy(:,:) 
        real(prec), intent(IN)  :: H_ice(:,:)
        real(prec), intent(IN)  :: dzsdx(:,:)
        real(prec), intent(IN)  :: dzsdy(:,:)
        real(prec), intent(IN)  :: dx 

        ! Local variables 
        integer :: i, j, nx, ny 
        real(prec) :: dy, rhog 
        real(prec) :: H_mid

        real(prec), allocatable :: Hi_ab(:,:) 

        nx = size(H_ice,1)
        ny = size(H_ice,2) 

        ! Allocate Hi_ab
!         allocate(Hi_ab(nx,ny))

        ! Stagger H_ice to Ab nodes:
        ! This will be used to calculate H_mid on the acx/acy nodes,
        ! but it should come from ab-nodes instead of ac-nodes for stability 
        ! Note: this is disabled, as it seemed not to affect results
!         Hi_ab = stagger_aa_ab_ice(H_ice,H_ice)
        
        ! Define shortcut parameter 
        rhog = rho_ice * g 

        ! Assume grid resolution is symmetrical 
        dy = dx 

        ! x-direction
        taud_acx = 0.0_prec  
        do j = 2, ny 
        do i = 1, nx-1 
!             H_mid         = 0.5_prec * (Hi_ab(i,j)+Hi_ab(i,j-1))
            H_mid         = 0.5_prec*(H_ice(i,j)+H_ice(i+1,j)) 
            taud_acx(i,j) = rhog * H_mid * dzsdx(i,j) 
        end do 
        end do 
        taud_acx(nx,:) = taud_acx(nx-1,:) 
        taud_acx(:,1)  = taud_acx(:,2) 

        ! y-direction
        taud_acy = 0.0_prec  
        do j = 1, ny-1 
        do i = 2, nx 
!             H_mid         = 0.5_prec * (Hi_ab(i,j)+Hi_ab(i-1,j))
            H_mid         = 0.5_prec*(H_ice(i,j)+H_ice(i,j+1))
            taud_acy(i,j) = rhog * H_mid * dzsdy(i,j) 
        end do 
        end do   
        taud_acy(:,ny) = taud_acy(:,ny-1)  
        taud_acy(1,:)  = taud_acy(2,:)

        return 

    end subroutine calc_driving_stress_ac 

    subroutine calc_driving_stress_gl_ac(taud_acx,taud_acy,H_ice,z_srf,z_bed,z_sl,H_grnd, &
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

                    write(*,*) "calc_driving_stress_gl_ac:: Error: Wrong choice of beta_gl_stag for this method."
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

                ! Note: this does not handle 'floating on the left' case correctly - fix! 
                stop 

                ! x-direction 
                do j = 1, ny 
                do i = 1, nx-1 
                    if (H_grnd(i,j) .gt. 0.0 .and. H_grnd(i+1,j) .le. 0.0) then 
                        ! Grounding line point 

                        ! (i,j) grounded; (i+1,j) floating
                        taud_acx(i,j) = integrate_gl_driving_stress_linear(H_ice(i,j),H_ice(i+1,j), &
                                                    z_bed(i,j),z_bed(i+1,j),z_sl(i,j), z_sl(i+1,j),dx)
                    
                    else if (H_grnd(i,j) .le. 0.0 .and. H_grnd(i+1,j) .gt. 0.0) then 
                        ! (i,j) floating; (i+1,j) grounded 

                        taud_acx(i,j) = integrate_gl_driving_stress_linear(H_ice(i+1,j),H_ice(i,j), &
                                                    z_bed(i+1,j),z_bed(i,j),z_sl(i+1,j),z_sl(i,j),dx)
                    end if 

                end do 
                end do 
                taud_acx(nx,:) = taud_acx(nx-1,:) 
        
                ! y-direction 
                do j = 1, ny-1 
                do i = 1, nx 
                    if (H_grnd(i,j) .gt. 0.0 .and. H_grnd(i,j+1) .le. 0.0) then
                        ! Grounding line point
 
                        ! (i,j) grounded; (i,j+1) floating 
                        taud_acy(i,j) = integrate_gl_driving_stress_linear(H_ice(i,j),H_ice(i,j+1), &
                                                    z_bed(i,j),z_bed(i,j+1),z_sl(i,j),z_sl(i,j+1),dx)

                    else if (H_grnd(i,j) .le. 0.0 .and. H_grnd(i,j+1) .gt. 0.0) then
                        ! (i,j) floating; (i,j+1) grounded

                        taud_acy(i,j) = integrate_gl_driving_stress_linear(H_ice(i,j+1),H_ice(i,j), &
                                                    z_bed(i,j+1),z_bed(i,j),z_sl(i,j+1),z_sl(i,j),dx)
                         
                    end if 
                end do 
                end do 
                taud_acy(:,ny) = taud_acy(:,ny-1) 
                
            case DEFAULT  
                ! Do nothing, use the standard no-subgrid treatment 
        end select 

        return 

    end subroutine calc_driving_stress_gl_ac 

    function integrate_gl_driving_stress_linear(H_ice,H_ice1,z_bed,z_bed1,z_sl,z_sl1,dx) result(taud)
        ! Compute the driving stress for the grounding line more precisely (subgrid)
        ! following Gladstone et al. (2010, TC), Eq. 27 
        ! Note: here cell i is grounded and cell i+1 is floating 
        ! Units: taud [Pa] 

        implicit none 

        real(prec), intent(IN) :: H_ice, H_ice1     ! Ice thickness cell i and i+1, resp.
        real(prec), intent(IN) :: z_bed, z_bed1     ! Bedrock elevation cell i and i+1, resp.
        real(prec), intent(IN) :: z_sl,  z_sl1      ! Sea level cell i and i+1, resp.
        real(prec), intent(IN) :: dx 
        real(prec) :: taud 

        ! Local variables 
        real(prec) :: Ha, Hb, Sa, Sb, Ba, Bb, sla, slb 
        real(prec) :: dl, dxab 
        real(prec) :: H_mid, dzsdx 
        real(prec) :: rho_sw_ice, rho_ice_sw 
        integer :: n 
        integer, parameter :: ntot = 100 

        ! Parameters 
        rho_sw_ice = rho_sw / rho_ice 
        rho_ice_sw = rho_ice / rho_sw 

        ! Get step size (dimensionless) and step resolution
        dl    = 1.0_prec / real(ntot,prec)
        dxab  = dx*dl 

        ! Initialize driving stress to zero 
        taud = 0.0_prec 

        ! Step through the grid cell and calculate
        ! each piecewise value of driving stress
        do n = 1, ntot 

            Ha  = H_ice + (H_ice1-H_ice)*dl*(n-1)
            Hb  = H_ice + (H_ice1-H_ice)*dl*(n)

            Ba  = z_bed + (z_bed1-z_bed)*dl*(n-1)
            Bb  = z_bed + (z_bed1-z_bed)*dl*(n)
            
            sla = z_sl + (z_sl1-z_sl)*dl*(n-1)
            slb = z_sl + (z_sl1-z_sl)*dl*(n)
            
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
            
            H_mid = 0.5_prec * (Ha+Hb)
            dzsdx = (Sb-Sa) / dxab 

            taud  = taud + (H_mid * dzsdx)*dl 

        end do 

        ! Finally multiply with rho_ice*g 
        taud = rho_ice*g *taud

        return 

    end function integrate_gl_driving_stress_linear 
    
    subroutine calc_ice_flux(qq_acx,qq_acy,ux_bar,uy_bar,H_ice,dx,dy)
        ! Calculate the basal stress resulting from sliding (friction times velocity)
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

end module velocity_hybrid_pd12
