module solver_advection
    
    use yelmo_defs, only : sp, dp, prec, tol_underflow
    use yelmo_tools, only : stagger_aa_ab
    use solver_advection_sico, only : calc_adv2D_expl_sico, calc_adv2D_impl_sico

    implicit none 



    private 
    public :: calc_advec2D
    public :: calc_adv2D_expl
    public :: calc_adv2D_impl_upwind

contains 

    subroutine calc_advec2D(var,ux,uy,var_dot,dx,dy,dt,solver)
        ! General routine to apply 2D advection equation to variable `var` 
        ! with source term `var_dot`. Various solvers are possible

        real(prec),       intent(INOUT) :: var(:,:)             ! [var] Variable to be advected
        real(prec),       intent(IN)    :: ux(:,:)              ! [m/a] 2D velocity, x-direction (ac-nodes)
        real(prec),       intent(IN)    :: uy(:,:)              ! [m/a] 2D velocity, y-direction (ac-nodes)
        real(prec),       intent(IN)    :: var_dot(:,:)         ! [dvar/dt] Source term for variable
        real(prec),       intent(IN)    :: dx                   ! [m]   Horizontal resolution, x-direction
        real(prec),       intent(IN)    :: dy                   ! [m]   Horizontal resolution, y-direction
        real(prec),       intent(IN)    :: dt                   ! [a]   Timestep 
        character(len=*), intent(IN)    :: solver               ! Solver to use for the ice thickness advection equation

        select case(trim(solver))
            ! Choose solver to use 

            case("expl")

                call calc_adv2D_expl(var,ux,uy,var_dot,dx,dy,dt)

            case("impl-upwind") 

                call calc_adv2D_impl_upwind(var,ux,uy,var_dot,dx,dy,dt,f_upwind=1.0_prec)
            
            ! Other solvers below...
            case("expl-sico")

                call calc_adv2D_expl_sico(var,ux,uy,var_dot,dx,dy,dt)

            case("impl-sico")

                call calc_adv2D_impl_sico(var,ux,uy,var_dot,dx,dy,dt,use_lis=.FALSE.)

            case("impl-sico-lis")

                call calc_adv2D_impl_sico(var,ux,uy,var_dot,dx,dy,dt,use_lis=.TRUE.)

            case DEFAULT 

                write(*,*) "calc_advec2D:: Error: solver not recognized."
                write(*,*) "solver = ", trim(solver)
                stop 

        end select 
        
        return 

    end subroutine calc_advec2D

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
        real(prec) :: reste, delh, tmp 
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

        ! Avoid underflows 
        where (abs(arelax) .lt. tol_underflow) arelax = 0.0_prec 
        where (abs(brelax) .lt. tol_underflow) brelax = 0.0_prec 
        where (abs(drelax) .lt. tol_underflow) drelax = 0.0_prec 
        where (abs(erelax) .lt. tol_underflow) erelax = 0.0_prec 
        
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
            where(abs(deltaH) .lt. tol_underflow) deltaH = 0.0_prec      ! Avoid underflows
            delh = sqrt(sum(deltaH**2)) / ((nx-2)*(ny-2))
            
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

end module solver_advection 