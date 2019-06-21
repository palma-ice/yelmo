module mass_conservation

    use yelmo_defs !, only :: sp, dp, prec 
    use yelmo_tools, only : stagger_aa_ab
    use mass_conservation_impl_sico, only : calc_adv2D_expl_sico, calc_adv2D_impl_sico

    implicit none 

    private
    public :: calc_ice_thickness

contains 

    subroutine calc_ice_thickness(H_ice,mb_applied,f_grnd,H_ocn,ux,uy,mbal,calv,z_bed_sd,dx,dt, &
                                    solver,boundaries,ice_allowed,H_min,sd_min,sd_max,calv_max)
        ! Interface subroutine to update ice thickness through application
        ! of advection, vertical mass balance terms and calving 

        implicit none 

        real(prec),       intent(INOUT) :: H_ice(:,:)           ! [m]   Ice thickness 
        real(prec),       intent(OUT)   :: mb_applied(:,:)      ! [m/a] Actual mass balance applied to real ice points
        real(prec),       intent(IN)    :: f_grnd(:,:)          ! [--]  Grounded fraction 
        real(prec),       intent(IN)    :: H_ocn(:,:)           ! [m]   Ocean thickness (ie, depth)
        real(prec),       intent(IN)    :: ux(:,:)              ! [m/a] Depth-averaged velocity, x-direction (ac-nodes)
        real(prec),       intent(IN)    :: uy(:,:)              ! [m/a] Depth-averaged velocity, y-direction (ac-nodes)
        real(prec),       intent(IN)    :: mbal(:,:)            ! [m/a] Net mass balance; mbal = smb+bmb  !-calv 
        real(prec),       intent(IN)    :: calv(:,:)            ! [m/a] Calving rate 
        real(prec),       intent(IN)    :: z_bed_sd(:,:)        ! [m]   Standard deviation of bed topography
        real(prec),       intent(IN)    :: dx                   ! [m]   Horizontal resolution
        real(prec),       intent(IN)    :: dt                   ! [a]   Timestep 
        character(len=*), intent(IN)    :: solver               ! Solver to use for the ice thickness advection equation
        character(len=*), intent(IN)    :: boundaries           ! What boundary conditions should apply?
        logical,          intent(IN)    :: ice_allowed(:,:)     ! Mask of where ice is allowed to be greater than zero 
        real(prec),       intent(IN)    :: H_min                ! [m]   Minimum allowed ice thickness parameter
        real(prec),       intent(IN)    :: sd_min               ! [m]   Minimum stdev(z_bed) parameter
        real(prec),       intent(IN)    :: sd_max               ! [m]   Maximum stdev(z_bed) parameter
        real(prec),       intent(IN)    :: calv_max             ! [m]   Maximum grounded calving rate parameter

        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: n 
        real(prec), allocatable :: calv_grnd(:,:) 

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        allocate(calv_grnd(nx,ny))
        calv_grnd = 0.0 

        ! 1. Apply mass conservation =================

        ! First, only resolve the dynamic part (ice advection)
        select case(trim(solver))
            ! Choose solver to use 

            case("expl")

                call calc_adv2D_expl(H_ice,ux,uy,mbal*0.0,dx,dx,dt)

            case("impl-upwind") 

                call calc_adv2D_impl_upwind(H_ice,ux,uy,mbal*0.0,dx,dx,dt,f_upwind=1.0)
            
            ! Other solvers below...
            case("expl-sico")

                call calc_adv2D_expl_sico(H_ice,ux,uy,mbal*0.0,dx,dx,dt)

            case("impl-sico")

                call calc_adv2D_impl_sico(H_ice,ux,uy,mbal*0.0,dx,dx,dt,use_lis=.FALSE.)

            case("impl-sico-lis")

                call calc_adv2D_impl_sico(H_ice,ux,uy,mbal*0.0,dx,dx,dt,use_lis=.TRUE.)

            case DEFAULT 

                write(*,*) "calc_ice_thickness:: Error: solver not recognized."
                write(*,*) "solver = ", trim(solver)
                stop 

        end select 
        
        ! Determine grounded calving 
        calv_grnd = calc_calving_rate_grounded(H_ice,f_grnd,z_bed_sd,sd_min,sd_max,calv_max,dt)

        ! Next, handle mass balance in order to be able to diagnose
        ! precisely how much mass was lost/gained 
        mb_applied = mbal - calv_grnd !- calv 

        ! Ensure ice cannot form in open ocean 
        where(f_grnd .eq. 0.0 .and. H_ice .eq. 0.0)  mb_applied = 0.0  

        ! Ensure melt is limited to amount of available ice to melt  
        where((H_ice+dt*mb_applied) .lt. 0.0) mb_applied = -H_ice/dt

        ! Apply modified mass balance to update the ice thickness 
        H_ice = H_ice + dt*mb_applied

        ! Limit grounded ice thickness to below maximum threshold value
        ! based on shear stress limit 
        call limit_grounded_margin_thickness_stress(H_ice,mb_applied,f_grnd,H_ocn,dt)

        ! Limit grounded ice thickess to above minimum and below inland neighbor at the margin
!         call limit_grounded_margin_thickness(H_ice,mb_applied,f_grnd,H_min,dt) 
        call limit_grounded_margin_thickness_flux(H_ice,mb_applied,f_grnd,mbal,ux,uy,dx,dt,H_min)

        ! Also ensure tiny numeric ice thicknesses are removed
        where (H_ice .lt. 1e-5) H_ice = 0.0 

        ! Finally, treat calving at the floating ice margin 
!         where (calv .gt. 0.0) H_ice = max(0.0,H_ice-calv*dt)
!         where (calv .gt. 0.0) H_ice = 0.0
        
        ! Artificially delete ice from locations that are not allowed
        ! according to boundary mask (ie, EISMINT, BUELER-A, open ocean)
        where (.not. ice_allowed) H_ice = 0.0 

        ! 2. Post processing of H_ice ================

        select case(trim(boundaries))

            case("zeros")

                ! Set border values to zero
                H_ice(1,:)  = 0.0
                H_ice(nx,:) = 0.0
                H_ice(:,1)  = 0.0
                H_ice(:,ny) = 0.0

            case("EISMINT")

                ! Set border values to zero
                H_ice(1,:)  = 0.0
                H_ice(nx,:) = 0.0
                H_ice(:,1)  = 0.0
                H_ice(:,ny) = 0.0
                
            case("MISMIP3D")

                ! === MISMIP3D =====
                H_ice(1,:)    = H_ice(2,:)       ! x=0, Symmetry 
                H_ice(nx,:)   = 0.0              ! x=800km, no ice
                H_ice(:,1)    = H_ice(:,2)       ! y=-50km, Free-slip condition
                H_ice(:,ny)   = H_ice(:,ny-1)    ! y= 50km, Free-slip condition

            case("infinite")
                ! ajr: we should check setting border H values equal to inner neighbors
                
                write(*,*) "calc_ice_thickness:: error: boundary method not implemented yet: "//trim(boundaries)
                write(*,*) "TO DO!"
                stop 

            case DEFAULT 

                write(*,*) "calc_ice_thickness:: error: boundary method not recognized: "//trim(boundaries)
                stop 

        end select 
        
        return 

    end subroutine calc_ice_thickness

    subroutine calc_adv2D_expl(H_ice, ux, uy, mdot, dx, dy, dt)
        ! Solve 2D advection equation for ice sheet thickness via explicit flux divergence:
        ! d[H]/dt = -grad[H*(ux,uy)] + mdot 
        !
        ! ajr: adapted from IMAU-ICE code from Heiko Goelzer (h.goelzer@uu.nl) 2016
        ! Note: original algorithm called for interpolation of H_ice to ab-nodes explicitly. 
        ! It also works using the ac-nodes directly, but is less stable. This can be chosen via the parameter
        ! use_ab_expl. 

        implicit none 

        real(prec), intent(INOUT) :: H_ice(:,:)             ! [m] aa-nodes, Ice thickness 
        real(prec), intent(IN)    :: ux(:,:)                ! [m a^-1] ac-nodes, Horizontal velocity, x-direction
        real(prec), intent(IN)    :: uy(:,:)                ! [m a^-1] ac-nodes, Horizontal velocity, y-direction
        real(prec), intent(IN)    :: mdot(:,:)              ! [m a^-1] aa-nodes, Source term, net rate of change at top and bottom of cell (no vertical advection) 
        real(prec), intent(IN)    :: dx                     ! [m] Horizontal grid spacing, x-direction
        real(prec), intent(IN)    :: dy                     ! [m] Horizontal grid spacing, y-direction
        real(prec), intent(IN)    :: dt                     ! [a] Timestep 

        ! Local variables:
        integer                 :: i, j, nx, ny 

        real(prec), allocatable :: dHdt(:,:)                ! [m] aa-nodes, Total change this timestep due to fluxes divergence and mdot 
        real(prec), allocatable :: flux_xr(:,:)             ! [m^2 a^-1] ac-nodes, Flux in the x-direction to the right
        real(prec), allocatable :: flux_xl(:,:)             ! [m^2 a^-1] ac-nodes, Flux in the x-direction to the left
        real(prec), allocatable :: flux_yu(:,:)             ! [m^2 a^-1] ac-nodes, Flux in the y-direction upwards
        real(prec), allocatable :: flux_yd(:,:)             ! [m^2 a^-1] ac-nodes, Flux in the y-direction downwards
        real(prec), allocatable :: H_ice_ab(:,:)            ! [m] ab-nodes, Ice thickness on staggered grid ab

        logical, parameter :: use_ab_expl = .TRUE.          ! Explicitly calculate ice thickness on ab nodes?
                                                            ! (ajr: this may be more stable, follows Macayeal staggering) 

        real(prec), parameter :: dHdt_lim = 1e3             ! [m a-1] Maximum rate of change allowed (high value for extreme changes)

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        allocate(dHdt(nx,ny))
        allocate(flux_xr(nx,ny))
        allocate(flux_xl(nx,ny))
        allocate(flux_yu(nx,ny))
        allocate(flux_yd(nx,ny))
        allocate(H_ice_ab(nx,ny))

        if (use_ab_expl) then 
            ! Use explicit calculation of H_ice on ab-nodes 

            ! Stagger ice thickness to ab-nodes 
            H_ice_ab = stagger_aa_ab(H_ice)

            ! Calculate the flux across each boundary [m^2 a^-1]
            do i = 2,nx-1
            do j = 2,ny-1
                flux_xr(i,j) = ux(i  ,j  ) * 0.5 * (H_ice_ab(i  ,j  ) + H_ice_ab(i  ,j-1))
                flux_xl(i,j) = ux(i-1,j  ) * 0.5 * (H_ice_ab(i-1,j  ) + H_ice_ab(i-1,j-1))
                flux_yu(i,j) = uy(i  ,j  ) * 0.5 * (H_ice_ab(i  ,j  ) + H_ice_ab(i-1,j  ))
                flux_yd(i,j) = uy(i  ,j-1) * 0.5 * (H_ice_ab(i  ,j-1) + H_ice_ab(i-1,j-1))
            end do
            end do

        else 
            ! Calculate fluxes on ac-nodes directly 
            
            ! Calculate the flux across each boundary [m^2 a^-1]
            do i = 2,nx-1
            do j = 2,ny-1
                flux_xr(i,j) = ux(i  ,j  ) * 0.5 * (H_ice(i  ,j  ) + H_ice(i+1,j  ))
                flux_xl(i,j) = ux(i-1,j  ) * 0.5 * (H_ice(i-1,j  ) + H_ice(i  ,j  ))
                flux_yu(i,j) = uy(i  ,j  ) * 0.5 * (H_ice(i  ,j  ) + H_ice(i  ,j+1))
                flux_yd(i,j) = uy(i  ,j-1) * 0.5 * (H_ice(i  ,j-1) + H_ice(i  ,j  ))
            end do
            end do

        end if 

        ! Calculate flux divergence on aa-nodes 
        dHdt = 0.0
        do i = 2,nx-1
        do j = 2,ny-1
            dHdt(i,j) = (1.0 / dx) * (flux_xl(i,j) - flux_xr(i,j)) + (1.0 / dy) * (flux_yd(i,j) - flux_yu(i,j))
        end do
        end do
        
        ! Limit dHdt for stability 
        where (dHdt .lt. -dHdt_lim) dHdt = -dHdt_lim 
        where (dHdt .gt.  dHdt_lim) dHdt =  dHdt_lim 
         
        ! Update H_ice:
        H_ice = H_ice + dt*dHdt

        ! Corners
        ! BL
        H_ice(1,1)   = H_ice(2,1)+H_ice(1,2)-H_ice(2,2)
        ! BR
        H_ice(nx,1)  = H_ice(nx-1,1)+H_ice(nx,2)-H_ice(nx-1,2)
        ! TL
        H_ice(1,ny)  = H_ice(2,ny)+H_ice(1,ny-1)-H_ice(2,ny-1)
        ! BL
        H_ice(nx,ny) = H_ice(nx-1,ny)+H_ice(nx,ny-1)-H_ice(nx-1,ny-1)

        return 

    end subroutine calc_adv2D_expl

    subroutine  calc_adv2D_impl_upwind(H,ux,uy,mdot,dx,dy,dt,f_upwind)
        ! To solve the 2D adevection equation:
        ! dH/dt =
        ! M H = Frelax
        ! ajr: adapted from GRISLI (Ritz et al., 1997)

        implicit none

        real(prec), intent(INOUT) :: H(:,:)         ! Ice thickness (aa-node)
        real(prec), intent(IN)    :: ux(:,:)        ! Depth-averaged velocity - x direction (ac-node)
        real(prec), intent(IN)    :: uy(:,:)        ! Depth-averaged velocity - y direction (ac-node)
        real(prec), intent(IN)    :: mdot(:,:)      ! Total column mass balance (aa-node)
        real(prec), intent(IN)    :: dx             ! [m] x-resolution
        real(prec), intent(IN)    :: dy             ! [m] y-resolution
        real(prec), intent(IN)    :: dt             ! [a] Timestep (assumes dx=dy)
        real(prec), intent(IN)    :: f_upwind       ! [-] Fraction of "upwind-ness" to apply (ajr: experimental!) - between 0.5 and 1.0, default f_upwind=1.0
        ! Local variables
        integer    :: i, j, nx, ny
        integer    :: iter, ierr 
        real(prec) :: dtdx, dtdx2
        real(prec) :: reste, delh
        real(prec), allocatable :: crelax(:,:)      ! diagnonale de M
        real(prec), allocatable :: arelax(:,:)      ! sous diagonale selon x
        real(prec), allocatable :: brelax(:,:)      ! sur  diagonale selon x
        real(prec), allocatable :: drelax(:,:)      ! sous diagonale selon y
        real(prec), allocatable :: erelax(:,:)      ! sur  diagonale selon y
        real(prec), allocatable :: frelax(:,:)      ! vecteur
        real(prec), allocatable :: c_west(:,:)      ! sur demi mailles Ux
        real(prec), allocatable :: c_east(:,:)      ! sur demi mailles Ux
        real(prec), allocatable :: c_north(:,:)     ! sur demi mailles Uy
        real(prec), allocatable :: c_south(:,:)     ! sur demi mailles Uy
        real(prec), allocatable :: deltaH(:,:)      ! Change in H

        logical,    parameter :: use_upwind = .TRUE.  ! Apply upwind advection scheme or central scheme?   
                                                      ! (now this is redundant with f_upwind parameter) 
        
        ! Note: f_upwind=0.6 gives good agreement with EISMINT1 summit elevation
        ! for the moving and fixed margin experiments, when using the calc_shear_3D approach.
        ! (f_upwind=0.5, ie central method, works well when using the velocity_sia approach)

        ! Note: it may be that f_upwind>0.5 is more stable for real dynamic simulations

        ! Determine array size
        nx = size(H,1)
        ny = size(H,2)

        ! Allocate local arrays
        allocate(crelax(nx,ny))
        allocate(arelax(nx,ny))
        allocate(brelax(nx,ny))
        allocate(drelax(nx,ny))
        allocate(erelax(nx,ny))
        allocate(frelax(nx,ny))
        allocate(c_west(nx,ny))
        allocate(c_east(nx,ny))
        allocate(c_north(nx,ny))
        allocate(c_south(nx,ny))

        allocate(deltaH(nx,ny))

        ! Define some helpful values
        dtdx2 = dt/(dx**2)
        dtdx  = dt/dx

        ! Initialize relaxation arrays
        arelax = 0.0
        brelax = 0.0
        drelax = 0.0
        erelax = 0.0
        crelax = 1.0
        frelax = 0.0
        deltaH = 0.0

        ! Modify coefficients depending on method (upwind, central)
        if (use_upwind) then
            ! Upwind method

            if (f_upwind .eq. 1.0) then
                ! Apply normal upwind scheme

                where (ux.ge.0.0)
                    c_west = 1.0
                    c_east = 0.0
                elsewhere
                    c_west = 0.0
                    c_east = 1.0
                end where

                where (uy.ge.0.0)
                    c_south = 1.0
                    c_north = 0.0
                elsewhere
                    c_south = 0.0
                    c_north = 1.0
                end where

            else 
                ! Apply fractional upwind scheme to reduce numerical diffusion,
                ! but benefit from upwind stability (ajr: experimental!)
                ! f_upwind = 0.5 => central difference, f_upwind = 1.0 => full upwind 

                where (ux.ge.0.0)
                    c_west = f_upwind
                    c_east = 1.0 - f_upwind
                elsewhere
                    c_west = 1.0 - f_upwind
                    c_east = f_upwind
                end where

                where (uy.ge.0.0)
                    c_south = f_upwind
                    c_north = 1.0 - f_upwind
                elsewhere
                    c_south = 1.0 - f_upwind
                    c_north = f_upwind
                end where

            end if 

        else
            ! Central method

            c_west  = 0.5
            c_east  = 0.5
            c_south = 0.5
            c_north = 0.5

        end if

        ! Populate diagonals
        do j=2,ny-1
        do i=2,nx-1

            !  sous diagonale en x
            arelax(i,j) = -dtdx*c_west(i-1,j)*ux(i-1,j)    ! partie advective en x

            !  sur diagonale en x
            brelax(i,j) = +dtdx*c_east(i,j)*ux(i,j)        ! partie advective

            !  sous diagonale en y
            drelax(i,j) = -dtdx*c_south(i,j-1)*uy(i,j-1)   ! partie advective en y

            !  sur diagonale en y
            erelax(i,j) = +dtdx*c_north(i,j)*uy(i,j)       ! partie advective


            ! diagonale
            crelax(i,j) = 1.0 + dtdx* &
                       ((c_west(i,j)*ux(i,j) - c_east(i-1,j)*ux(i-1,j)) &
                      +(c_south(i,j)*uy(i,j) - c_north(i,j-1)*uy(i,j-1)))

            ! Combine all terms
            frelax(i,j) = H(i,j) + dt*mdot(i,j)

        end do
        end do

        ! Initialize new H solution to zero (to get zeros at boundaries)
        H  = 0.0

        ! Initially assume convergence criterion is not satisfied 
        ierr = -1   ! convergence criterion not fulfilled
        
        do iter = 1, 1000 
            ! Relaxation loop 

            ! Calculate change in H
            do j=2,ny-1
            do i=2,nx-1

                reste = (((arelax(i,j)*H(i-1,j) + drelax(i,j)*H(i,j-1)) &
                        + (brelax(i,j)*H(i+1,j) + erelax(i,j)*H(i,j+1))) &
                        +  crelax(i,j)*H(i,j))  - frelax(i,j)

                deltaH(i,j) = reste/crelax(i,j)

            end do
            end do

            ! Adjust H to new value
            H = H - deltaH

            ! Check stopping criterion (something like rmse of remaining change in H)
            delh = 0
            do j=2,ny-1
            do i=2,nx-1
                delh = delh + deltaH(i,j)**2
            end do
            end do
            delh = sqrt(delh)/((nx-2)*(ny-2))
            
            ! Use simple stopping criterion: maximum remaining change in H
            ! Note: this is less likely to converge given the same stopping
            ! criterion.
!             delh = maxval(abs(deltaH))
            
            if ( delh .lt. 1e-4) then
                ! Solution has converged, exit  
                ierr = 0 
                exit 
            end if 

        end do ! End of relaxation loop

        !write(6,'(10x,a,i0,5x,i2)') 'calc_adv2D_impl_upwind: iter = ', iter, ierr
        
        return

    end subroutine calc_adv2D_impl_upwind

    subroutine limit_grounded_margin_thickness(H_ice,mb_applied,f_grnd,H_min,dt)

        implicit none 

        real(prec), intent(INOUT) :: H_ice(:,:) 
        real(prec), intent(INOUT) :: mb_applied(:,:) 
        real(prec), intent(IN)    :: f_grnd(:,:) 
        real(prec), intent(IN)    :: H_min
        real(prec), intent(IN)    :: dt
        
        ! Local variables 
        integer :: i, j, nx, ny 
        logical :: is_margin 
        real(prec) :: H_min_neighbor, H_diff 

        real(prec), allocatable :: H_ice_0(:,:) 

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        allocate(H_ice_0(nx,ny))
        H_ice_0 = H_ice 

        do j = 2, ny-1 
        do i = 2, nx-1 

            ! Determine if ice-covered point has an ice-free neighbor (ie, at the ice margin)
            is_margin = (H_ice_0(i,j) .gt. 0.0 .and. f_grnd(i,j) .gt. 0.0 &
                .and. minval([H_ice_0(i-1,j),H_ice_0(i+1,j),H_ice_0(i,j-1),H_ice_0(i,j+1)]) .eq. 0.0)

            if (is_margin .and. H_ice_0(i,j) .lt. H_min .and. mb_applied(i,j) .le. 0.0) then 
                ! If grounded margin point is too thin and mass balance not positive, 
                ! impose ablation to set ice thickness to zero

                mb_applied(i,j) = mb_applied(i,j) - H_ice_0(i,j)/dt
                H_ice(i,j)      = 0.0 

            else if (is_margin) then 
                ! Check if margin is thicker than inland neighbors

                H_min_neighbor = minval([H_ice_0(i-1,j),H_ice_0(i+1,j),H_ice_0(i,j-1),H_ice_0(i,j+1)], &
                                    mask=[H_ice_0(i-1,j),H_ice_0(i+1,j),H_ice_0(i,j-1),H_ice_0(i,j+1)] .gt. 0.0)

                if(H_ice_0(i,j) .gt. H_min_neighbor) then 
                    H_diff = H_ice_0(i,j)-0.5*H_min_neighbor
                    mb_applied(i,j) = mb_applied(i,j) - H_diff/dt 
                    H_ice(i,j)      = H_ice_0(i,j)-H_diff
                end if 

            end if
            
        end do 
        end do 

        return 

    end subroutine limit_grounded_margin_thickness

    subroutine limit_grounded_margin_thickness_flux(H_ice,mb_applied,f_grnd,mbal,ux,uy,dx,dt,H_min)
        ! Remove marginal ice that is too thin (H_ice < H_min) with mass balance
        ! and cannot be replenished from upstream flux

        implicit none 

        real(prec), intent(INOUT) :: H_ice(:,:)         ! [m] Ice thickness 
        real(prec), intent(INOUT) :: mb_applied(:,:) 
        real(prec), intent(IN) :: f_grnd(:,:)           ! [-] Grounded fraction
        real(prec), intent(IN) :: mbal(:,:)             ! [m/a] Net mass balance 
        real(prec), intent(IN) :: ux(:,:)               ! [m/a] velocity, x-direction (ac-nodes)
        real(prec), intent(IN) :: uy(:,:)               ! [m/a] velocity, y-direction (ac-nodes)
        real(prec), intent(IN) :: dx, dt 
        real(prec), intent(IN) :: H_min                 ! [m] Threshold for calving

        ! Local variables 
        integer :: i, j, nx, ny
        real(prec) :: eps_xx, eps_yy  
        logical :: test_mij, test_pij, test_imj, test_ipj
        logical :: is_margin, positive_mb 
        real(prec), allocatable :: dHdt(:,:), Hdiff(:,:), Hfrac(:,:) 
        real(prec), allocatable :: f_ice(:,:)                ! [-] Ice area fraction
        real(prec), allocatable :: H_ice_0(:,:) 
        integer :: n_free

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        allocate(dHdt(nx,ny))
        allocate(Hdiff(nx,ny))
        allocate(Hfrac(nx,ny))
        allocate(f_ice(nx,ny))
        allocate(H_ice_0(nx,ny))
        
        ! Local representation of f_ice (binary) - improve in future
        f_ice = 0.0
        where(H_ice .gt. 0.0) f_ice = 1.0 

        ! Determine ice thickness accounting for fraction  
        where (H_ice .gt. 1.0 .and. f_ice .gt. 0.0) 
            Hfrac = H_ice / f_ice 
        elsewhere 
            Hfrac = 0.0
        end where 

        ! Ice thickness above threshold
        Hdiff = Hfrac - H_min

        ! Diagnosed lagrangian rate of change
        dHdt = 0.0 

        do j = 2, ny
        do i = 2, nx
        
                ! Calculate strain rate locally (aa-node)
                ! Note: dx should probably account for f_ice somehow,  
                ! but would maybe only be a minor adjustment
                eps_xx = (ux(i,j) - ux(i-1,j))/dx
                eps_yy = (uy(i,j) - uy(i,j-1))/dx

                ! Calculate thickness change via conservation
                dHdt(i,j) = mbal(i,j) - Hfrac(i,j)*(eps_xx+eps_yy)

        end do 
        end do
        
        ! Store original ice thickness field 
        H_ice_0 = H_ice 
        
        do j = 2, ny-1
        do i = 2, nx-1

            ! Determine how many ice-free points are bordering grounded ice point

            ! Determine if grounded, ice-covered point has an ice-free neighbor (ie, at the grounded ice margin)
            is_margin = (H_ice_0(i,j) .gt. 0.0 .and. f_grnd(i,j) .gt. 0.0 &
                .and. minval([H_ice_0(i-1,j),H_ice_0(i+1,j),H_ice_0(i,j-1),H_ice_0(i,j+1)]) .eq. 0.0)

            if (f_grnd(i,j) .gt. 0.0 .and. f_ice(i,j) .gt. 0.0) then
                ! Grounded ice-covered point

                n_free = count([f_ice(i-1,j).eq.0.0, &
                                f_ice(i+1,j).eq.0.0, &
                                f_ice(i,j-1).eq.0.0, &
                                f_ice(i,j+1).eq.0.0])

            else
                ! No ice-free points bordering point 

                n_free = 0 

            end if

            if (n_free .gt. 0 .and. Hdiff(i,j).le.0.0) then 
                ! Check if current point is at the margin,
                ! and has thickness less than threshold, or if
                ! ice below H_min limit, accounting for mass flux from inland

                positive_mb = (mbal(i,j).gt.0.0)

                test_mij = ( ((Hdiff(i-1,j).gt.0.0).and.(ux(i,j).ge.0.0)  &  ! neighbor (i-1,j) total > H_min
                    .and.  (dHdt(i-1,j).gt.(-Hdiff(i-1,j)*abs(ux(i-1,j)/dx)))) & 
                    .or.(f_grnd(i-1,j).gt.0.0.and.positive_mb) ) !

                test_pij = ( ((Hdiff(i+1,j).gt.0.0).and.(ux(i,j).le.0.0) & ! neighbor (i+1,j) total > H_min
                    .and.(dHdt(i+1,j).gt.(-Hdiff(i+1,j)*abs(ux(i,j)/dx)))) &
                    .or.(f_grnd(i+1,j).gt.0.0.and.positive_mb) ) !

                test_imj = ( ((Hdiff(i,j-1).gt.0.0).and.(uy(i,j).ge.0.0)  &  ! neighbor (i,j-1) total > H_min
                    .and.(dHdt(i,j-1).gt.(-Hdiff(i,j-1)*abs(uy(i,j-1)/dx))))&
                    .or.(f_grnd(i,j-1).gt.0.0.and.positive_mb) ) !

                test_ipj = ( ((Hdiff(i,j+1).gt.0.0).and.(uy(i,j).le.0.0) & ! neighbor (i,j+1) total > H_min
                    .and.(dHdt(i,j+1).gt.(-Hdiff(i,j+1)*abs(uy(i,j)/dx))))&
                    .or.(f_grnd(i,j+1).gt.0.0.and.positive_mb) ) !

                if ((.not.(test_mij.or.test_pij.or.test_imj.or.test_ipj))) then
                    mb_applied(i,j) = mb_applied(i,j) - H_ice_0(i,j) / dt    
                    H_ice(i,j)      = 0.0          
                end if  

            end if

        end do
        end do

        return 

    end subroutine limit_grounded_margin_thickness_flux
    
    subroutine limit_grounded_margin_thickness_stress(H_ice,mb_applied,f_grnd,H_ocn,dt)
        ! Remove marginal ice that exceeds a stress threshold following
        ! Bassis and Walker (2012), Eq. 2.12 

        implicit none 

        real(prec), intent(INOUT) :: H_ice(:,:)             ! [m] Ice thickness 
        real(prec), intent(INOUT) :: mb_applied(:,:)        ! [m/a] Applied mass balance
        real(prec), intent(IN)    :: f_grnd(:,:)            ! [-] Grounded fraction
        real(prec), intent(IN)    :: H_ocn(:,:)             ! [m] Ocean thickness (depth)
        real(prec), intent(IN)    :: dt 

        ! Local variables 
        integer    :: i, j, nx, ny
        real(prec) :: tau_c, H_max, H_ocn_now 
        logical    :: is_margin  
        real(prec) :: rho_ice_g, rho_sw_ice, rho_ice_sw  
        real(prec), allocatable :: H_ice_0(:,:) 

        real(prec), parameter :: C0    = 1e6                ! [Pa] Depth-averaged shear stress in ice 
        real(prec), parameter :: alpha = 0.0                ! [--] Friction coefficient for Bassis and Walker (2012), Eq. 2.13
        real(prec), parameter :: r     = 0.0                ! [--] Crevasse fraction 
        
        rho_ice_g  = rho_ice * g 
        rho_sw_ice = rho_sw / rho_ice 
        rho_ice_sw = rho_ice / rho_sw 

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        allocate(H_ice_0(nx,ny))
        H_ice_0 = H_ice 
        
        do j = 2, ny-1
        do i = 2, nx-1 

            ! Determine if grounded, ice-covered point has an ice-free neighbor (ie, at the grounded ice margin)
            is_margin = (H_ice_0(i,j) .gt. 0.0 .and. f_grnd(i,j) .gt. 0.0 &
                .and. minval([H_ice_0(i-1,j),H_ice_0(i+1,j),H_ice_0(i,j-1),H_ice_0(i,j+1)]) .eq. 0.0)

            if (is_margin) then 
                ! Determine if this margin point should fail 

                ! Calculate depth of seawater (limited by ice thickness and flotation criterion)
                if (H_ocn(i,j) .gt. 0.0) then 
                    H_ocn_now = min(rho_ice_sw*H_ice_0(i,j),H_ocn(i,j))
                else 
                    H_ocn_now = 0.0 
                end if 

                ! Get depth-averaged shear-stress in ice, Bassis and Walker (2012), Eq. 2.13 vertically integrated
                ! alpha = 0.65: model S1 validated for cold ice 
                ! alpha = 0.4 : model S2 for warmer ice 
                ! alpha = 0.0 : model S3 for purely plastic yielding (default)
                tau_c = C0 + 0.5*alpha*rho_ice_g*H_ice_0(i,j)

                ! Get critical ice thickness to cause stress failure
                H_max = (1.0-r)*tau_c/rho_ice_g + sqrt(((1.0-r)*tau_c/rho_ice_g)**2 + rho_sw_ice*H_ocn_now**2)

                if (H_ice_0(i,j) .gt. H_max) then 
                    ! Critical stress exceeded, calve this ice - pass to mb_applied for now...

                    mb_applied(i,j) = mb_applied(i,j) - H_ice_0(i,j) / dt    
                    H_ice(i,j)      = 0.0  

                end if 

            end if 

        end do 
        end do 

        return 

    end subroutine limit_grounded_margin_thickness_stress
    
    function calc_calving_rate_grounded(H_ice,f_grnd,z_bed_sd,sd_min,sd_max,calv_max,dt) result(calv)
        ! Parameterize grounded ice-margin calving as a function of 
        ! standard deviation of bedrock at each grid point.
        ! Assumes that higher variability in subgrid implies cliffs
        ! that are not represented at low resolution. 

        implicit none 

        real(prec), intent(IN) :: H_ice(:,:)                ! [m] Ice thickness 
        real(prec), intent(IN) :: f_grnd(:,:)               ! [-] Grounded fraction
        real(prec), intent(IN) :: z_bed_sd(:,:)             ! [m] Standard deviation of bedrock topography
        real(prec), intent(IN) :: sd_min                    ! [m] stdev(z_bed) at/below which calv=0
        real(prec), intent(IN) :: sd_max                    ! [m] stdev(z_bed) at/above which calv=calv_max 
        real(prec), intent(IN) :: calv_max                  ! [m/a] Maximum allowed calving rate
        real(prec), intent(IN) :: dt 
        real(prec) :: calv(size(H_ice,1),size(H_ice,2))     ! [m/a] Calculated calving rate 

        ! Local variables
        integer :: i, j, nx, ny  
        real(prec) :: f_scale 
        logical    :: is_grnd_margin 

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        calv = 0.0 
        
        if (calv_max .gt. 0.0) then 
            ! Determine grounded calving rate 

            do j = 2, ny-1
            do i = 2, nx-1 

                ! Determine if grounded, ice-covered point has an ice-free neighbor (ie, at the grounded ice margin)
                is_grnd_margin = (H_ice(i,j) .gt. 0.0 .and. f_grnd(i,j) .gt. 0.0 &
                    .and. minval([H_ice(i-1,j),H_ice(i+1,j),H_ice(i,j-1),H_ice(i,j+1)]) .eq. 0.0)

                if (is_grnd_margin) then
                    ! Grounded ice-covered point

                    f_scale = (z_bed_sd(i,j) - sd_min)/(sd_max-sd_min)
                    if (f_scale .lt. 0.0) f_scale = 0.0 
                    if (f_scale .gt. 1.0) f_scale = 1.0 

                    ! Calculate calving rate from linear function, limited
                    ! to available ice thickness 
                    calv(i,j) = min(f_scale*calv_max, H_ice(i,j)/dt) 
                    
                end if 

            end do 
            end do 

        end if 
        
        return 

    end function calc_calving_rate_grounded

    subroutine calc_ice_margin(H_ice,H_ref,f_ice,f_grnd)
        ! Determine the area fraction of a grid cell
        ! that is ice-covered. Assume that marginal points
        ! have equal thickness to inland neighbors 

        implicit none 

        real(prec), intent(INOUT) :: H_ice(:,:)                ! [m] Ice thickness on standard grid (aa-nodes)
        real(prec), intent(INOUT) :: H_ref(:,:)                ! [m] Margin ice thickness for partially filled cells, V_cell = H_ref*f_ice
        real(prec), intent(INOUT) :: f_ice(:,:)                ! [--] Ice covered fraction (aa-nodes)
        real(prec), intent(IN) :: f_grnd(:,:)               ! [--] Grounded fraction (aa-nodes)

        ! Local variables 
        integer :: i, j, nx, ny 
        real(prec) :: H_neighb(4)
        logical :: mask_neighb(4)
        real(prec) :: H_mrgn 

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        ! Initially set fraction to one everywhere there is ice 
        ! and zero everywhere there is no ice
        f_ice = 0.0  
        where (H_ice .gt. 0.0) f_ice = 1.0

        ! For ice-covered points with ice-free neighbors (ie, at the floating or grounded margin),
        ! determine the fraction of grid point that should be ice covered. 

        do j = 2, ny-1
        do i = 2, nx-1 

            if (f_ice(i,j) .gt. 0.0 .and. &
                count([f_ice(i-1,j),f_ice(i+1,j),f_ice(i,j-1),f_ice(i,j+1)].eq.0) .gt. 0) then 
                ! This point is at the ice margin

                H_neighb    = [H_ice(i-1,j),H_ice(i+1,j),H_ice(i,j-1),H_ice(i,j+1)]
                mask_neighb = (H_neighb .gt. 0.0)

                if (count(mask_neighb) .gt. 0) then 
                    ! Neighbors with ice should generally be found, but put this check just in case

                    ! Determine height to give to partially filled cell as average of neighbors
                    H_mrgn = sum(H_neighb,mask=mask_neighb)/real(count(mask_neighb))

                    ! If margin point is grounded, then assign it with 
                    ! a thickness of half of neighbor-average
                    if (f_grnd(i,j) .eq. 1.0) H_mrgn = 0.5 * H_mrgn

                    ! Determine the cell ice fraction
                    ! Note: fraction is determined as a ratio of 
                    ! thicknesses, derived from volume conservation 
                    ! vol = H_ice*dx*dy = H_mrgn*area_frac 
                    ! f_ice = area_frac / (dx*dy)
                    ! f_ice = H_ice/H_mrgn 
                    f_ice(i,j) = min( H_ice(i,j) / H_mrgn, 1.0 ) 

                end if

            end if  

        end do 
        end do 

        return 

    end subroutine calc_ice_margin
    
end module mass_conservation
