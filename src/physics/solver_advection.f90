module solver_advection
    
    use yelmo_defs, only : sp, dp, wp, tol_underflow
    use yelmo_tools, only : get_neighbor_indices, stagger_aa_ab, stagger_nodes_aa_ab_ice
    use solver_advection_sico, only : calc_adv2D_expl_sico, calc_adv2D_impl_sico

    implicit none 



    private 
    public :: calc_advec2D
    public :: calc_adv2D_expl
    public :: calc_adv2D_impl_upwind

contains 

    subroutine calc_advec2D(dvdt,var,f_ice,ux,uy,var_dot,dx,dy,dt,solver,boundaries)
        ! General routine to apply 2D advection equation to variable `var` 
        ! with source term `var_dot`. Various solvers are possible

        real(wp),       intent(OUT)   :: dvdt(:,:)            ! [dvdt] Variable rate of change
        real(wp),       intent(IN)    :: var(:,:)             ! [var]  Variable to be advected
        real(wp),       intent(IN)    :: f_ice(:,:)             ! [var]  Variable to be advected
        real(wp),       intent(IN)    :: ux(:,:)              ! [m/a] 2D velocity, x-direction (ac-nodes)
        real(wp),       intent(IN)    :: uy(:,:)              ! [m/a] 2D velocity, y-direction (ac-nodes)
        real(wp),       intent(IN)    :: var_dot(:,:)         ! [dvar/dt] Source term for variable
        real(wp),       intent(IN)    :: dx                   ! [m]   Horizontal resolution, x-direction
        real(wp),       intent(IN)    :: dy                   ! [m]   Horizontal resolution, y-direction
        real(wp),       intent(IN)    :: dt                   ! [a]   Timestep 
        character(len=*), intent(IN)    :: solver               ! Solver to use for the ice thickness advection equation
        character(len=*), intent(IN)    :: boundaries           ! Boundary conditions to impose

        ! Local variables 
        real(wp), allocatable :: var_now(:,:) 

        allocate(var_now(size(var,1),size(var,2)))

        ! Assign local variable to be modified 
        var_now = var 

        select case(trim(solver))
            ! Choose solver to use 

            case("none") 

                ! Do nothing: no advection considered. 

            case("expl")

                call calc_adv2D_expl(var_now,f_ice,ux,uy,var_dot,dx,dy,dt,boundaries)

            case("expl-new")

                call calc_adv2D_expl_new(var_now,f_ice,ux,uy,var_dot,dx,dy,dt,boundaries)

            case("expl-upwind")

                call calc_adv2D_expl_upwind(var_now,ux,uy,var_dot,dx,dy,dt,boundaries)

            case("impl-upwind") 

                call calc_adv2D_impl_upwind(var_now,ux,uy,var_dot,dx,dy,dt,boundaries,f_upwind=1.0_wp)
            
            ! Other solvers below...
            case("expl-sico")

                call calc_adv2D_expl_sico(var_now,ux,uy,var_dot,dx,dy,dt)

            case("impl-sico")

                call calc_adv2D_impl_sico(var_now,ux,uy,var_dot,dx,dy,dt,use_lis=.FALSE.)

            case("impl-sico-lis")

                call calc_adv2D_impl_sico(var_now,ux,uy,var_dot,dx,dy,dt,use_lis=.TRUE.)

            case DEFAULT 

                write(*,*) "calc_advec2D:: Error: solver not recognized."
                write(*,*) "solver = ", trim(solver)
                stop 

        end select 
        
        ! Determine rate of change 
        dvdt = (var_now - var) / dt 

        return 

    end subroutine calc_advec2D

    subroutine calc_adv2D_expl_new(H_ice, f_ice, ux, uy, mdot, dx, dy, dt, boundaries)
        ! Solve 2D advection equation for ice sheet thickness via explicit flux divergence:
        ! d[H]/dt = -grad[H*(ux,uy)] + mdot 
        
        implicit none 

        real(wp), intent(INOUT) :: H_ice(:,:)             ! [m] aa-nodes, Ice thickness 
        real(wp), intent(IN)    :: f_ice(:,:)
        real(wp), intent(IN)    :: ux(:,:)                ! [m a^-1] ac-nodes, Horizontal velocity, x-direction
        real(wp), intent(IN)    :: uy(:,:)                ! [m a^-1] ac-nodes, Horizontal velocity, y-direction
        real(wp), intent(IN)    :: mdot(:,:)              ! [m a^-1] aa-nodes, Source term, net rate of change at top and bottom of cell (no vertical advection) 
        real(wp), intent(IN)    :: dx                     ! [m] Horizontal grid spacing, x-direction
        real(wp), intent(IN)    :: dy                     ! [m] Horizontal grid spacing, y-direction
        real(wp), intent(IN)    :: dt                     ! [a] Timestep 
        character(len=*), intent(IN) :: boundaries 

        ! Local variables:
        integer    :: i, j, nx, ny 
        integer    :: im1, ip1, jm1, jp1
        integer    :: im2, ip2, jm2, jp2
        real(wp)   :: H_ac(2)
        real(wp)   :: wt_ab(4)
        real(wp)   :: H_ab(4)
        real(wp)   :: ux_ab(4)
        real(wp)   :: uy_ab(4)
        real(wp)   :: fx_ab(4)
        real(wp)   :: fy_ab(4) 
        real(wp)   :: wt
        real(wp)   :: ux_now, uy_now 
        real(wp)   :: dfxdx, dfydy 
        real(wp), allocatable :: dHdt(:,:)                ! [m] aa-nodes, Total change this timestep due to fluxes divergence and mdot 
        real(wp), allocatable :: fx(:,:) 
        real(wp), allocatable :: fy(:,:) 

        real(wp), parameter :: dHdt_lim = 1e3             ! [m a-1] Maximum rate of change allowed (high value for extreme changes)

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        ! Initialize dHdt 
        allocate(dHdt(nx,ny))
        dHdt = 0.0_wp 

        allocate(fx(nx,ny))
        allocate(fy(nx,ny))

        fx = 0.0_wp 
        fy = 0.0_wp 

        ! First calculate flux on aa-nodes over entire domain 
        do i = 1, nx
        do j = 1, ny

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

            ! Get ab-node weighting based on whether ice is present 
            wt_ab = 0.0_wp 
            if (count([f_ice(i,j),f_ice(ip1,j),f_ice(i,jp1),f_ice(ip1,jp1)].lt.1.0_wp) .eq. 0) then 
                wt_ab(1) = 1.0_wp
            end if
            if (count([f_ice(i,j),f_ice(im1,j),f_ice(i,jp1),f_ice(im1,jp1)].lt.1.0_wp) .eq. 0) then 
                wt_ab(2) = 1.0_wp 
            end if 
            if (count([f_ice(i,j),f_ice(im1,j),f_ice(i,jm1),f_ice(im1,jm1)].lt.1.0_wp) .eq. 0) then 
                wt_ab(3) = 1.0_wp
            end if 
            if (count([f_ice(i,j),f_ice(ip1,j),f_ice(i,jm1),f_ice(ip1,jm1)].lt.1.0_wp) .eq. 0) then 
                wt_ab(4) = 1.0_wp 
            end if 
            
            wt = sum(wt_ab)
            
            H_ab(1) = 0.25_wp*(H_ice(i,j)+H_ice(ip1,j)+H_ice(i,jp1)+H_ice(ip1,jp1)) 
            H_ab(2) = 0.25_wp*(H_ice(i,j)+H_ice(im1,j)+H_ice(i,jp1)+H_ice(im1,jp1)) 
            H_ab(3) = 0.25_wp*(H_ice(i,j)+H_ice(im1,j)+H_ice(i,jm1)+H_ice(im1,jm1)) 
            H_ab(4) = 0.25_wp*(H_ice(i,j)+H_ice(ip1,j)+H_ice(i,jm1)+H_ice(ip1,jm1)) 

            ! === ux =========

            ux_ab(1) = 0.5_wp*(ux(i,j)+ux(i,jp1))
            ux_ab(2) = 0.5_wp*(ux(im1,j)+ux(im1,jp1))
            ux_ab(3) = 0.5_wp*(ux(im1,jm1)+ux(im1,j))
            ux_ab(4) = 0.5_wp*(ux(i,jm1)+ux(i,j))
            where (abs(ux_ab) .lt. TOL_UNDERFLOW) ux_ab = 0.0_wp 

            ! === uy =========

            uy_ab(1) = 0.5_wp*(uy(i,j)+uy(ip1,j))
            uy_ab(2) = 0.5_wp*(uy(i,j)+uy(im1,j))
            uy_ab(3) = 0.5_wp*(uy(i,jm1)+uy(im1,jm1))
            uy_ab(4) = 0.5_wp*(uy(i,jm1)+uy(ip1,jm1))
            where (abs(uy_ab) .lt. TOL_UNDERFLOW) uy_ab = 0.0_wp 

            fx_ab = H_ab*ux_ab
            fy_ab = H_ab*uy_ab 

            dfxdx = 0.5_wp*( (fx_ab(1)-fx_ab(2))/dx + (fx_ab(4)-fx_ab(3))/dx )
            dfydy = 0.5_wp*( (fy_ab(1)-fy_ab(4))/dy + (fy_ab(2)-fy_ab(3))/dy )
            
            dHdt(i,j) = -dfxdx - dfydy 

            ! Limit dHdt for stability 
            if (dHdt(i,j) .lt. -dHdt_lim) dHdt(i,j) = -dHdt_lim 
            if (dHdt(i,j) .gt.  dHdt_lim) dHdt(i,j) =  dHdt_lim 
            


            ! ! Get centered fluxes (aa-nodes)
            ! fx(i,j) = 0.25_wp*(sum(H_ab*ux_ab))
            ! fy(i,j) = 0.25_wp*(sum(H_ab*uy_ab))
            
            ! ! Get centered fluxes (aa-nodes)
            ! H_ac(1) = 0.5_wp*(H_ice(i,j)+H_ice(im1,j))
            ! H_ac(2) = 0.5_wp*(H_ice(i,j)+H_ice(ip1,j))
            ! fx(i,j) = 0.5_wp*(H_ac(1)*ux(im1,j) + H_ac(2)*ux(i,j))

            ! H_ac(1) = 0.5_wp*(H_ice(i,j)+H_ice(i,jm1))
            ! H_ac(2) = 0.5_wp*(H_ice(i,j)+H_ice(i,jp1))
            ! fy(i,j) = 0.5_wp*(H_ac(1)*uy(i,jm1) + H_ac(2)*uy(i,j))
            
        end do
        end do 

if (.FALSE.) then
        ! Now calculate spatial derivatives dfx/dx / dfy/dy
        ! to obtain dh/dt
        do i = 1, nx
        do j = 1, ny

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

            !im2 = max(i-2,1)
            !ip2 = min(i+2,nx)
            !jm2 = max(j-2,1)
            !jp2 = min(j+2,ny)
            
            ! ajr: works for simple (EISMINT tests), but not for non-linear
            ! realistic ice sheets
            ! Second-order spatial derivative following
            ! https://web.media.mit.edu/~crtaylor/calculator.html
            dfxdx = (1.0_wp/(2.0_wp*dx))*(-1.0_wp*fx(im1,j)+1.0_wp*fx(ip1,j))
            dfydy = (1.0_wp/(2.0_wp*dy))*(-1.0_wp*fy(i,jm1)+1.0_wp*fy(i,jp1))

            ! ajr: works for simple (EISMINT tests), but not for non-linear
            ! realistic ice sheets
            ! Fourth order spatial derivative following
            ! Man and Tsai (2008).
            ! dfxdx = (1.0_wp/(12.0_wp*dx))*(fx(im2,j)-8.0_wp*fx(im1,j)+8.0_wp*fx(ip1,j)-fx(ip2,j))
            ! dfydy = (1.0_wp/(12.0_wp*dy))*(fy(i,jm2)-8.0_wp*fy(i,jm1)+8.0_wp*fy(i,jp1)-fy(i,jp1))


            dHdt(i,j) = -dfxdx - dfydy 

            ! Limit dHdt for stability 
            if (dHdt(i,j) .lt. -dHdt_lim) dHdt(i,j) = -dHdt_lim 
            if (dHdt(i,j) .gt.  dHdt_lim) dHdt(i,j) =  dHdt_lim 
            
        end do
        end do

end if 

        ! Update H_ice:
        H_ice = H_ice + dt*dHdt + dt*mdot 

        return 

    end subroutine calc_adv2D_expl_new

    subroutine calc_adv2D_expl(H_ice, f_ice, ux, uy, mdot, dx, dy, dt, boundaries)
        ! Solve 2D advection equation for ice sheet thickness via explicit flux divergence:
        ! d[H]/dt = -grad[H*(ux,uy)] + mdot 
        !
        ! ajr: adapted from IMAU-ICE code from Heiko Goelzer (h.goelzer@uu.nl) 2016
        ! Note: It also works using the ac-nodes directly, but is less stable, particularly
        ! for EISMINT2-expf. 

        implicit none 

        real(wp), intent(INOUT) :: H_ice(:,:)               ! [m] aa-nodes, Ice thickness 
        real(wp), intent(IN)    :: f_ice(:,:)               ! Ice fraction
        real(wp), intent(IN)    :: ux(:,:)                  ! [m a^-1] ac-nodes, Horizontal velocity, x-direction
        real(wp), intent(IN)    :: uy(:,:)                  ! [m a^-1] ac-nodes, Horizontal velocity, y-direction
        real(wp), intent(IN)    :: mdot(:,:)                ! [m a^-1] aa-nodes, Source term, net rate of change at top and bottom of cell (no vertical advection) 
        real(wp), intent(IN)    :: dx                       ! [m] Horizontal grid spacing, x-direction
        real(wp), intent(IN)    :: dy                       ! [m] Horizontal grid spacing, y-direction
        real(wp), intent(IN)    :: dt                       ! [a] Timestep 
        character(len=*), intent(IN) :: boundaries 
        
        ! Local variables:
        integer  :: i, j, nx, ny 
        integer  :: im1, ip1, jm1, jp1
        real(wp) :: H_ab(4)
        real(wp) :: flux_xr                                 ! [m^2 a^-1] ac-nodes, Flux in the x-direction to the right
        real(wp) :: flux_xl                                 ! [m^2 a^-1] ac-nodes, Flux in the x-direction to the left
        real(wp) :: flux_yu                                 ! [m^2 a^-1] ac-nodes, Flux in the y-direction upwards
        real(wp) :: flux_yd                                 ! [m^2 a^-1] ac-nodes, Flux in the y-direction downwards
        real(wp), allocatable :: dHdt(:,:)                  ! [m] aa-nodes, Total change this timestep due to fluxes divergence and mdot 
        real(wp), allocatable :: f_ice_now(:,:)
        real(wp), parameter :: dHdt_lim = 1e3               ! [m a-1] Maximum rate of change allowed (high value for extreme changes)

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        ! Initialize dHdt 
        allocate(dHdt(nx,ny))
        dHdt = 0.0_wp 

        allocate(f_ice_now(nx,ny))
        f_ice_now = 1.0 

        do j = 1, ny
        do i = 1, nx
        
            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

            ! Get the ice thickness on ab-nodes
            call stagger_nodes_aa_ab_ice(H_ab,H_ice,f_ice_now,i,j)

            ! Calculate the flux across each boundary [m^2 a^-1]
            flux_xr = ux(i,j)   * 0.5*(H_ab(1)+H_ab(4))
            flux_xl = ux(im1,j) * 0.5*(H_ab(2)+H_ab(3))
            flux_yu = uy(i,j)   * 0.5*(H_ab(1)+H_ab(2))
            flux_yd = uy(i,jm1) * 0.5*(H_ab(3)+H_ab(4))

            ! Calculate flux divergence on aa-nodes 
            dHdt(i,j) = (1.0 / dx) * (flux_xl - flux_xr) + (1.0 / dy) * (flux_yd - flux_yu)

            ! Limit dHdt for stability 
            if (dHdt(i,j) .lt. -dHdt_lim) dHdt(i,j) = -dHdt_lim 
            if (dHdt(i,j) .gt.  dHdt_lim) dHdt(i,j) =  dHdt_lim 
            
        end do
        end do

        ! Update H_ice:
        H_ice = H_ice + dt*dHdt + dt*mdot 

        return 

    end subroutine calc_adv2D_expl

    subroutine calc_adv2D_expl_upwind(H_ice, ux, uy, mdot, dx, dy, dt, boundaries)
        ! Solve 2D advection equation for ice sheet thickness via explicit flux divergence:
        ! d[H]/dt = -grad[H*(ux,uy)] + mdot 

        ! Second-order upwind routine, see eg, Winkelmann et al (2011), Eq. 18 

        implicit none 

        real(wp), intent(INOUT) :: H_ice(:,:)             ! [m] aa-nodes, Ice thickness 
        real(wp), intent(IN)    :: ux(:,:)                ! [m a^-1] ac-nodes, Horizontal velocity, x-direction
        real(wp), intent(IN)    :: uy(:,:)                ! [m a^-1] ac-nodes, Horizontal velocity, y-direction
        real(wp), intent(IN)    :: mdot(:,:)              ! [m a^-1] aa-nodes, Source term, net rate of change at top and bottom of cell (no vertical advection) 
        real(wp), intent(IN)    :: dx                     ! [m] Horizontal grid spacing, x-direction
        real(wp), intent(IN)    :: dy                     ! [m] Horizontal grid spacing, y-direction
        real(wp), intent(IN)    :: dt                     ! [a] Timestep 
        character(len=*), intent(IN) :: boundaries 
        
        ! Local variables:
        integer    :: i, j, nx, ny
        integer    :: im1, ip1, jm1, jp1 
        real(wp) :: H_now 
        real(wp) :: flux_xr                               ! [m^2 a^-1] ac-nodes, Flux in the x-direction to the right
        real(wp) :: flux_xl                               ! [m^2 a^-1] ac-nodes, Flux in the x-direction to the left
        real(wp) :: flux_yu                               ! [m^2 a^-1] ac-nodes, Flux in the y-direction upwards
        real(wp) :: flux_yd                               ! [m^2 a^-1] ac-nodes, Flux in the y-direction downwards
        real(wp), allocatable :: dHdt(:,:)                ! [m] aa-nodes, Total change this timestep due to fluxes divergence and mdot 

        real(wp), parameter :: dHdt_lim = 1e3             ! [m a-1] Maximum rate of change allowed (high value for extreme changes)

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        ! Initialize dHdt 
        allocate(dHdt(nx,ny))
        dHdt = 0.0_wp 

        do j = 1, ny
        do i = 1, nx
        
            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

            ! Calculate the flux across each boundary [m^2 a^-1]
            if (ux(i,j) .gt. 0.0) then 
                flux_xr = ux(i,j)*H_ice(i,j)
            else
                flux_xr = ux(i,j)*H_ice(ip1,j) 
            end if 

            if (ux(im1,j) .gt. 0.0) then 
                flux_xl = ux(im1,j)*H_ice(im1,j)
            else
                flux_xl = ux(im1,j)*H_ice(i,j) 
            end if 

            if (uy(i,j) .gt. 0.0) then 
                flux_yu = uy(i,j)*H_ice(i,j)
            else
                flux_yu = uy(i,j)*H_ice(i,jp1) 
            end if 

            if (uy(i,jm1) .gt. 0.0) then 
                flux_yd = uy(i,jm1)*H_ice(i,jm1)
            else
                flux_yd = uy(i,jm1)*H_ice(i,j) 
            end if 

            ! Calculate flux divergence on aa-nodes 
            dHdt(i,j) = (1.0 / dx) * (flux_xl - flux_xr) + (1.0 / dy) * (flux_yd - flux_yu)

            ! Limit dHdt for stability 
            if (dHdt(i,j) .lt. -dHdt_lim) dHdt(i,j) = -dHdt_lim 
            if (dHdt(i,j) .gt.  dHdt_lim) dHdt(i,j) =  dHdt_lim 
            
        end do
        end do

        ! Update H_ice:
        H_ice = H_ice + dt*dHdt + dt*mdot 

        return 

    end subroutine calc_adv2D_expl_upwind

    subroutine  calc_adv2D_impl_upwind(H,ux,uy,mdot,dx,dy,dt,boundaries,f_upwind)
        ! To solve the 2D adevection equation:
        ! dH/dt =
        ! M H = Frelax
        ! ajr: adapted from GRISLI (Ritz et al., 1997)

        implicit none

        real(wp), intent(INOUT) :: H(:,:)         ! Ice thickness (aa-node)
        real(wp), intent(IN)    :: ux(:,:)        ! Depth-averaged velocity - x direction (ac-node)
        real(wp), intent(IN)    :: uy(:,:)        ! Depth-averaged velocity - y direction (ac-node)
        real(wp), intent(IN)    :: mdot(:,:)      ! Total column mass balance (aa-node)
        real(wp), intent(IN)    :: dx             ! [m] x-resolution
        real(wp), intent(IN)    :: dy             ! [m] y-resolution
        real(wp), intent(IN)    :: dt             ! [a] Timestep (assumes dx=dy)
        character(len=*), intent(IN) :: boundaries 
        real(wp), intent(IN)    :: f_upwind       ! [-] Fraction of "upwind-ness" to apply (ajr: experimental!) - between 0.5 and 1.0, default f_upwind=1.0
        
        ! Local variables
        integer    :: i, j, nx, ny
        integer    :: im1, ip1, jm1, jp1 
        integer    :: iter, ierr 
        real(wp) :: dtdx, dtdx2
        real(wp) :: reste, delh, tmp 
        real(wp) :: ux_i, ux_im1, uy_j, uy_jm1
        real(wp), allocatable :: crelax(:,:)      ! diagnonale de M
        real(wp), allocatable :: arelax(:,:)      ! sous diagonale selon x
        real(wp), allocatable :: brelax(:,:)      ! sur  diagonale selon x
        real(wp), allocatable :: drelax(:,:)      ! sous diagonale selon y
        real(wp), allocatable :: erelax(:,:)      ! sur  diagonale selon y
        real(wp), allocatable :: frelax(:,:)      ! vecteur
        real(wp), allocatable :: c_west(:,:)      ! sur demi mailles Ux
        real(wp), allocatable :: c_east(:,:)      ! sur demi mailles Ux
        real(wp), allocatable :: c_north(:,:)     ! sur demi mailles Uy
        real(wp), allocatable :: c_south(:,:)     ! sur demi mailles Uy
        real(wp), allocatable :: deltaH(:,:)      ! Change in H

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
        do j = 1, ny
        do i = 1, nx

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

            !  sous diagonale en x
            ux_im1 = ux(im1,j)
            if (abs(ux_im1) .lt. TOL_UNDERFLOW) ux_im1 = 0.0
            arelax(i,j) = -dtdx*c_west(im1,j)*ux_im1    ! partie advective en x

            !  sur diagonale en x
            ux_i = ux(i,j)
            if (abs(ux_i) .lt. TOL_UNDERFLOW) ux_i = 0.0
            brelax(i,j) = +dtdx*c_east(i,j)*ux_i        ! partie advective

            !  sous diagonale en y
            uy_jm1 = uy(i,jm1)
            if (abs(uy_jm1) .lt. TOL_UNDERFLOW) uy_jm1 = 0.0
            drelax(i,j) = -dtdx*c_south(i,jm1)*uy_jm1   ! partie advective en y

            !  sur diagonale en y
            uy_j = uy(i,j)
            if (abs(uy_j) .lt. TOL_UNDERFLOW) uy_j = 0.0
            erelax(i,j) = +dtdx*c_north(i,j)*uy_j       ! partie advective


            ! diagonale
            crelax(i,j) = 1.0 + dtdx* &
                       ((c_west(i,j)*ux_i - c_east(im1,j)*ux_im1) &
                      +(c_south(i,j)*uy_j - c_north(i,jm1)*uy_jm1))

            ! Combine all terms
            frelax(i,j) = H(i,j) + dt*mdot(i,j)

        end do
        end do

        ! Avoid underflows 
        where (abs(arelax) .lt. tol_underflow) arelax = 0.0_wp 
        where (abs(brelax) .lt. tol_underflow) brelax = 0.0_wp 
        where (abs(drelax) .lt. tol_underflow) drelax = 0.0_wp 
        where (abs(erelax) .lt. tol_underflow) erelax = 0.0_wp 
        
        ! Initialize new H solution to zero (to get zeros at boundaries)
        H  = 0.0

        ! Initially assume convergence criterion is not satisfied 
        ierr = -1   ! convergence criterion not fulfilled
        
        do iter = 1, 1000 
            ! Relaxation loop 

            ! Calculate change in H
            do j = 1, ny
            do i = 1, nx

                ! Get neighbor indices
                call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

                reste = (((arelax(i,j)*H(im1,j) + drelax(i,j)*H(i,jm1)) &
                        + (brelax(i,j)*H(ip1,j) + erelax(i,j)*H(i,jp1))) &
                        +  crelax(i,j)*H(i,j))  - frelax(i,j)

                deltaH(i,j) = reste/crelax(i,j)

            end do
            end do

            ! Avoid underflows
            where(abs(deltaH) .lt. tol_underflow) deltaH = 0.0_wp
            
            ! Adjust H to new value
            H = H - deltaH

            ! Avoid underflows
            where(abs(H) .lt. tol_underflow) H = 0.0_wp
            
            ! Check stopping criterion (something like rmse of remaining change in H)
            delh = sqrt(sum(deltaH**2)) / ((nx-2)*(ny-2))
            
            ! Use simple stopping criterion: maximum remaining change in H
            ! Note: this is less likely to converge given the same stopping
            ! criterion.
!             delh = maxval(abs(deltaH))
            
            if ( delh .lt. 1e-6) then
                ! Solution has converged, exit  
                ierr = 0 
                exit 
            end if 

        end do ! End of relaxation loop

        !write(6,'(10x,a,i0,5x,i2)') 'calc_adv2D_impl_upwind: iter = ', iter, ierr
        
        return

    end subroutine calc_adv2D_impl_upwind

end module solver_advection
