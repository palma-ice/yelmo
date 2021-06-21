module mass_conservation

    use yelmo_defs, only : sp, dp, wp, TOL_UNDERFLOW, g, rho_ice, rho_sw  
    use yelmo_tools, only : fill_borders_2D

    use solver_advection, only : calc_advec2D  
    use velocity_general, only : set_inactive_margins 

    implicit none 

    private
    public :: calc_ice_thickness_dyn
    public :: calc_ice_thickness_mbal
    public :: apply_ice_thickness_boundaries
    public :: relax_ice_thickness

contains 

    subroutine calc_ice_thickness_dyn(H_ice,dHdt_n,H_ice_n,H_ice_pred,f_ice,f_grnd,ux,uy, &
                                      solver,dx,dt,beta,pc_step)
        ! Interface subroutine to update ice thickness through application
        ! of advection, vertical mass balance terms and calving 

        implicit none 

        real(wp),         intent(INOUT) :: H_ice(:,:)           ! [m]   Ice thickness 
        real(wp),         intent(INOUT) :: dHdt_n(:,:)          ! [m/a] Advective rate of ice thickness change from previous=>current timestep 
        real(wp),         intent(INOUT) :: H_ice_n(:,:)         ! [m]   Ice thickness from previous=>current timestep 
        real(wp),         intent(INOUT) :: H_ice_pred(:,:)      ! [m]   Ice thickness from predicted timestep 
        real(wp),         intent(IN)    :: f_ice(:,:)           ! [--]  Ice area fraction 
        real(wp),         intent(IN)    :: f_grnd(:,:)          ! [--]  Fraction of grounded ice 
        real(wp),         intent(IN)    :: ux(:,:)              ! [m/a] Depth-averaged velocity, x-direction (ac-nodes)
        real(wp),         intent(IN)    :: uy(:,:)              ! [m/a] Depth-averaged velocity, y-direction (ac-nodes)
        character(len=*), intent(IN)    :: solver               ! Solver to use for the ice thickness advection equation
        real(wp),         intent(IN)    :: dx                   ! [m]   Horizontal resolution
        real(wp),         intent(IN)    :: dt                   ! [a]   Timestep 
        real(wp),         intent(IN)    :: beta(4)              ! Timestep weighting parameters
        character(len=*), intent(IN)    :: pc_step              ! Current predictor-corrector step ('predictor' or 'corrector')

        ! Local variables 
        integer :: i, j, nx, ny
        integer :: im1, ip1, jm1, jp1  
        real(wp), allocatable :: mbal_zero(:,:) 
        real(wp), allocatable :: dHdt_advec(:,:) 
        real(wp), allocatable :: ux_tmp(:,:) 
        real(wp), allocatable :: uy_tmp(:,:) 

        real(wp), parameter :: dHdt_advec_lim = 10.0_wp     ! [m/a] Hard limit on advection rate

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        allocate(mbal_zero(nx,ny))
        mbal_zero = 0.0_wp 

        allocate(ux_tmp(nx,ny))
        allocate(uy_tmp(nx,ny))
        allocate(dHdt_advec(nx,ny))

        dHdt_advec = 0.0_wp 

        ! Set local velocity fields with no margin treatment intially
        ux_tmp = ux
        uy_tmp = uy
        
        ! Ensure that no velocity is defined for outer boundaries of partially-filled margin points
        call set_inactive_margins(ux_tmp,uy_tmp,f_ice)

        ! ===================================================================================
        ! Resolve the dynamic part (ice advection) using multistep method

        select case(trim(pc_step))
        
            case("predictor") 
                
                ! Store ice thickness from time=n
                H_ice_n   = H_ice 

                ! Store advective rate of change from saved from previous timestep (now represents time=n-1)
                dHdt_advec = dHdt_n 

                ! Determine current advective rate of change (time=n)
                call calc_advec2D(dHdt_n,H_ice,ux_tmp,uy_tmp,mbal_zero,dx,dx,dt,solver)

                ! ajr: testing stability fix for spin-up, limit advection rate!
                ! where(dHdt_n .gt.  dHdt_advec_lim) dHdt_n = dHdt_advec_lim
                ! where(dHdt_n .lt. -dHdt_advec_lim) dHdt_n = -dHdt_advec_lim
                
                ! Calculate rate of change using weighted advective rates of change 
                dHdt_advec = beta(1)*dHdt_n + beta(2)*dHdt_advec 
                
                ! Calculate predicted ice thickness (time=n+1,pred)
                H_ice = H_ice_n + dt*dHdt_advec 

            case("corrector") ! corrector 

                ! Determine advective rate of change based on predicted H,ux/y fields (time=n+1,pred)
                call calc_advec2D(dHdt_advec,H_ice_pred,ux_tmp,uy_tmp,mbal_zero,dx,dx,dt,solver)

                ! ajr: testing stability fix for spin-up, limit advection rate!
                ! where(dHdt_advec .gt.  dHdt_advec_lim) dHdt_advec = dHdt_advec_lim
                ! where(dHdt_advec .lt. -dHdt_advec_lim) dHdt_advec = -dHdt_advec_lim
                
                ! Calculate rate of change using weighted advective rates of change 
                dHdt_advec = beta(3)*dHdt_advec + beta(4)*dHdt_n 
                
                ! Calculate corrected ice thickness (time=n+1)
                H_ice = H_ice_n + dt*dHdt_advec 

                ! Finally, update dHdt_n with correct term to use as n-1 on next iteration
                dHdt_n = dHdt_advec 

        end select

        ! Post processing of H_ice ================

        ! Also ensure tiny numeric ice thicknesses are removed
        where (H_ice .lt. 1e-5) H_ice = 0.0 

        return 

    end subroutine calc_ice_thickness_dyn

    subroutine calc_ice_thickness_mbal(H_ice,mb_applied,calv,f_ice,f_grnd,H_ocn, &
                                       mbal,calv_flt,calv_grnd,dx,dt,reset)
        ! Interface subroutine to update ice thickness through application
        ! of advection, vertical mass balance terms and calving 

        implicit none 

        real(wp),       intent(INOUT) :: H_ice(:,:)             ! [m]   Ice thickness 
        real(wp),       intent(INOUT) :: mb_applied(:,:)        ! [m/a] Actual mass balance applied to real ice points
        real(wp),       intent(INOUT) :: calv(:,:)              ! [m/a] Actual calving rate 
        real(wp),       intent(IN)    :: f_ice(:,:)             ! [--]  Ice area fraction 
        real(wp),       intent(IN)    :: f_grnd(:,:)            ! [--]  Grounded fraction 
        real(wp),       intent(IN)    :: H_ocn(:,:)             ! [m]   Ocean thickness (ie, depth)
        real(wp),       intent(IN)    :: mbal(:,:)              ! [m/a] Net mass balance; mbal = smb+bmb (calving separate) 
        real(wp),       intent(IN)    :: calv_flt(:,:)          ! [m/a] Potential calving rate (floating)
        real(wp),       intent(IN)    :: calv_grnd(:,:)         ! [m/a] Potential calving rate (grounded)
        real(wp),       intent(IN)    :: dx                     ! [m]   Horizontal resolution
        real(wp),       intent(IN)    :: dt                     ! [a]   Timestep 
        logical,        intent(IN)    :: reset 

        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: im1, ip1, jm1, jp1 
        real(wp) :: calv_applied
        real(wp), allocatable :: calv_resid(:,:) 
        real(wp) :: wts(4)

        nx = size(H_ice,1)
        ny = size(H_ice,2) 

        allocate(calv_resid(nx,ny))
        calv_resid = 0.0 

        if (reset) then 
            ! Make sure to first initialize calv and mb_applied to 
            ! zero 

            calv        = 0.0_wp 
            mb_applied  = 0.0_wp 

        end if 

        ! ==== MASS BALANCE =====

        ! First apply mass balance 
        mb_applied = mb_applied + mbal

        ! Ensure ice cannot form in open ocean 
        where(f_grnd .eq. 0.0 .and. H_ice .eq. 0.0)  mb_applied = 0.0  

        ! Ensure melt is limited to amount of available ice to melt  
        where((H_ice+dt*mb_applied) .lt. 0.0) mb_applied = -H_ice/dt


        ! ===== CALVING ======
        ! Limit calving contributions to margin points 
        ! (for now assume this was done well externally)

        ! Combine grounded and floating calving into one field for output
        calv = calv + (calv_flt + calv_grnd) 


        ! Diagnose residual calving 
        where((H_ice-dt*calv) .lt. 0.0_wp) 
            calv_resid = -(H_ice-dt*calv)/dt
        elsewhere
            calv_resid = 0.0_wp
        end where

        ! Eliminate residual calving from calving field 
        calv = calv - calv_resid 
        where (abs(calv) .lt. TOL_UNDERFLOW) calv = 0.0_wp 

if (.TRUE.) then
    ! Ensure that excess calving gets applied to upstream neighbors with 
    ! full ice cover. 

        ! Determine if any residual calving should migrate to neighboring inland
        ! cells if all ice in the margin cell is deleted.
        do j = 1, ny 
        do i = 1, nx 

            ! Define neighbor indices
            im1 = max(i-1,1)
            ip1 = min(i+1,nx)
            jm1 = max(j-1,1)
            jp1 = min(j+1,ny)
            
            if (calv_resid(i,j) .gt. 0.0_wp) then 
                ! Calving diagnosed for this point 

                if (f_ice(ip1,j) .eq. 1.0) wts(1) = 1.0_wp 
                if (f_ice(i,jp1) .eq. 1.0) wts(2) = 1.0_wp 
                if (f_ice(im1,j) .eq. 1.0) wts(3) = 1.0_wp 
                if (f_ice(i,jm1) .eq. 1.0) wts(4) = 1.0_wp 
                
                if (sum(wts) .eq. 0.0) then 
                    ! This shouldn't happen, something went wrong! 
                    write(*,*) "calc_ice_thickness_mbal:: Error: &
                    &calving point found with no fully ice-covered neighbors. Check!"
                    write(*,*) "i, j: ", i, j 
                else
                    wts = wts / sum(wts) 
                end if 

                calv(ip1,j) = calv(ip1,j) + calv_resid(i,j)*wts(1)
                calv(i,jp1) = calv(i,jp1) + calv_resid(i,j)*wts(2)
                calv(im1,j) = calv(im1,j) + calv_resid(i,j)*wts(3)
                calv(i,jm1) = calv(i,jm1) + calv_resid(i,j)*wts(4)

                calv_resid(i,j) = calv_resid(i,j) - sum(wts)*calv_resid(i,j)
            
            end if 
        end do 
        end do 

        if (maxval(abs(calv_resid)) .gt. 1e-3) then 
            write(*,*) "calc_ice_thickness_mbal:: Error: residual calving not &
            & properly accounted for."
            write(*,*) "calv_resid: ", minval(calv_resid), maxval(calv_resid)
            stop 
        end if 
end if 

        ! Subtract calving from mb_applied 
        mb_applied = mb_applied - calv 

        ! Apply modified mass balance to update the ice thickness 
        H_ice = H_ice + dt*mb_applied

        ! Ensure tiny numeric ice thicknesses are removed
        where (H_ice .lt. 1e-5) H_ice = 0.0 

        return 

    end subroutine calc_ice_thickness_mbal

    subroutine apply_ice_thickness_boundaries(H_ice,mb_resid,f_ice,f_grnd,uxy_b,ice_allowed,boundaries,H_ice_ref, &
                                                H_min_flt,H_min_grnd,dt)

        implicit none

        real(wp),           intent(INOUT)   :: H_ice(:,:)               ! [m] Ice thickness 
        real(wp),           intent(OUT)     :: mb_resid(:,:)            ! [m/yr] Residual mass balance
        real(wp),           intent(IN)      :: f_ice(:,:)               ! [--] Fraction of ice cover
        real(wp),           intent(IN)      :: f_grnd(:,:)              ! [--] Grounded ice fraction
        real(wp),           intent(IN)      :: uxy_b(:,:)               ! [m/a] Basal sliding speed, aa-nodes
        logical,            intent(IN)      :: ice_allowed(:,:)         ! Mask of where ice is allowed to be greater than zero 
        character(len=*),   intent(IN)      :: boundaries               ! Boundary condition choice
        real(wp),           intent(IN)      :: H_ice_ref(:,:)           ! [m]  Reference ice thickness to fill with for boundaries=="fixed"
        real(wp),           intent(IN)      :: H_min_flt                ! [m] Minimum allowed floating ice thickness 
        real(wp),           intent(IN)      :: H_min_grnd               ! [m] Minimum allowed grounded ice thickness 
        real(wp),           intent(IN)      :: dt                       ! [yr] Timestep

        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: im1, ip1, jm1, jp1 
        real(wp), allocatable :: H_ice_new(:,:)
        real(wp) :: H_eff 
        logical  :: is_margin 
        logical  :: is_island 
        logical  :: is_isthmus_x 
        logical  :: is_isthmus_y 

        nx = size(H_ice,1)
        ny = size(H_ice,2) 

        allocate(H_ice_new(nx,ny)) 
        H_ice_new = H_ice 

        ! Apply special case for symmetric EISMINT domain when basal sliding is active
        ! (ensure summit thickness does not grow disproportionately)
        if (trim(boundaries) .eq. "EISMINT" .and. maxval(uxy_b) .gt. 0.0) then 
            i = (nx-1)/2 
            j = (ny-1)/2
            H_ice_new(i,j) = (H_ice(i-1,j)+H_ice(i+1,j) &
                                    +H_ice(i,j-1)+H_ice(i,j+1)) / 4.0 
        end if  
        
        ! Artificially delete ice from locations that are not allowed
        ! according to boundary mask (ie, EISMINT, BUELER-A, open ocean)
        where (.not. ice_allowed) H_ice_new = 0.0 

        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices 
            im1 = max(i-1,1)
            ip1 = min(i+1,nx)
            jm1 = max(j-1,1)
            jp1 = min(j+1,ny)

            is_margin = f_ice(i,j) .gt. 0.0 .and. &
                count([f_ice(im1,j),f_ice(ip1,j),f_ice(i,jm1),f_ice(i,jp1)].eq.0.0) .gt. 0

            if (is_margin) then
                ! Ice covered point at the margin

                ! Calculate current ice thickness 
                if (f_ice(i,j) .gt. 0.0) then 
                    H_eff = H_ice_new(i,j) / f_ice(i,j) 
                else 
                    H_eff = H_ice_new(i,j) 
                end if 

                ! Remove ice that is too thin 
                if (f_grnd(i,j) .lt. 1.0_wp .and. H_eff .lt. H_min_flt)  H_ice_new(i,j) = 0.0_wp 
                if (f_grnd(i,j) .eq. 1.0_wp .and. H_eff .lt. H_min_grnd) H_ice_new(i,j) = 0.0_wp 

            end if 

            ! Check for ice islands
            is_island = H_ice(i,j) .gt. 0.0 .and. &
                count([H_ice(im1,j),H_ice(ip1,j),H_ice(i,jm1),H_ice(i,jp1)].gt.0.0) .eq. 0

! ajr: Treatment of isthmuses makes the model 
! more unstable. Disable for now - maybe it is not needed.
if (.FALSE.) then
            is_isthmus_x = H_ice(i,j) .gt. 0.0 .and. &
                count([H_ice(im1,j),H_ice(ip1,j)].gt.0.0) .eq. 0

            is_isthmus_y = H_ice(i,j) .gt. 0.0 .and. &
                count([H_ice(i,jm1),H_ice(i,jp1)].gt.0.0) .eq. 0
else       
            is_isthmus_x = .FALSE. 
            is_isthmus_y = .FALSE. 
end if 

            if (is_island) then 
                ! Ice-covered island
                ! Redistribute ice to neighbors
                ! (slowly island will disappear or merge with bigger ice body)

                H_ice_new(im1,j) = H_ice_new(im1,j) + 0.125_wp*H_ice(i,j)
                H_ice_new(ip1,j) = H_ice_new(ip1,j) + 0.125_wp*H_ice(i,j)
                H_ice_new(i,jm1) = H_ice_new(i,jm1) + 0.125_wp*H_ice(i,j)
                H_ice_new(i,jp1) = H_ice_new(i,jp1) + 0.125_wp*H_ice(i,j)
                H_ice_new(i,j)   = 0.500_wp*H_ice(i,j)
            
            else if (is_isthmus_x) then 

                H_ice_new(im1,j) = H_ice_new(im1,j) + 0.250_wp*H_ice(i,j)
                H_ice_new(ip1,j) = H_ice_new(ip1,j) + 0.250_wp*H_ice(i,j)
                H_ice_new(i,j)   = 0.500_wp*H_ice(i,j)

            else if (is_isthmus_y) then 

                H_ice_new(i,jm1) = H_ice_new(i,jm1) + 0.250_wp*H_ice(i,j)
                H_ice_new(i,jp1) = H_ice_new(i,jp1) + 0.250_wp*H_ice(i,j)
                H_ice_new(i,j)   = 0.500_wp*H_ice(i,j)

            end if 

        end do 
        end do 

        select case(trim(boundaries))

            case("zeros","EISMINT")

                ! Set border values to zero
                H_ice_new(1,:)  = 0.0
                H_ice_new(nx,:) = 0.0

                H_ice_new(:,1)  = 0.0
                H_ice_new(:,ny) = 0.0

            case("periodic","periodic-xy") 

                H_ice_new(1:2,:)     = H_ice_new(nx-3:nx-2,:) 
                H_ice_new(nx-1:nx,:) = H_ice_new(2:3,:) 

                H_ice_new(:,1:2)     = H_ice_new(:,ny-3:ny-2) 
                H_ice_new(:,ny-1:ny) = H_ice_new(:,2:3) 
            
            case("periodic-x") 

                ! Periodic x 
                H_ice_new(1:2,:)     = H_ice_new(nx-3:nx-2,:) 
                H_ice_new(nx-1:nx,:) = H_ice_new(2:3,:) 
                
                ! Infinite (free-slip too)
                H_ice_new(:,1)  = H_ice_new(:,2)
                H_ice_new(:,ny) = H_ice_new(:,ny-1)

            case("MISMIP3D")

                ! === MISMIP3D =====
                H_ice_new(1,:)    = H_ice_new(2,:)       ! x=0, Symmetry 
                H_ice_new(nx,:)   = 0.0              ! x=800km, no ice
                
                H_ice_new(:,1)    = H_ice_new(:,2)       ! y=-50km, Free-slip condition
                H_ice_new(:,ny)   = H_ice_new(:,ny-1)    ! y= 50km, Free-slip condition

            case("infinite")
                ! Set border points equal to inner neighbors 

                call fill_borders_2D(H_ice_new,nfill=1)

            case("fixed") 
                ! Set border points equal to prescribed values from array

                call fill_borders_2D(H_ice_new,nfill=1,fill=H_ice_ref)

            case DEFAULT 

                write(*,*) "apply_ice_thickness_boundaries:: error: boundary method not recognized: "//trim(boundaries)
                stop 

        end select 
        
        ! Determine mass balance related to changes applied here 
        mb_resid = (H_ice_new - H_ice) / dt 

        ! Reset actual ice thickness to new values 
        H_ice = H_ice_new 

        return

    end subroutine apply_ice_thickness_boundaries

    subroutine relax_ice_thickness(H_ice,f_grnd,H_ref,topo_rel,tau,dt)
        ! This routines allows ice within a given mask to be
        ! relaxed to a reference state with certain timescale tau 
        ! (if tau=0), then H_ice = H_ice_ref directly 

        implicit none 

        real(wp), intent(INOUT) :: H_ice(:,:) 
        real(wp), intent(IN)    :: f_grnd(:,:)  
        real(wp), intent(IN)    :: H_ref(:,:) 
        integer,    intent(IN)    :: topo_rel 
        real(wp), intent(IN)    :: tau
        real(wp), intent(IN)    :: dt 

        ! Local variables 
        integer    :: i, j, nx, ny 
        logical    :: apply_relax 
        real(wp) :: dHdt 

        nx = size(H_ice,1)
        ny = size(H_ice,2) 

        do j = 2, ny-1 
            do i = 2, nx-1 

                ! No relaxation to start
                apply_relax = .FALSE.

                select case(topo_rel)

                    case(1) 
                        ! Relax the shelf (floating) ice and ice-free points
                    
                        if (f_grnd(i,j) .eq. 0.0 .or. H_ref(i,j) .eq. 0.0) apply_relax = .TRUE. 
                
                    case(2) 
                        ! Relax the shelf (floating) ice and ice-free points
                        ! and the grounding-line ice too
                        
                        if (f_grnd(i,j) .eq. 0.0 .or. H_ref(i,j) .eq. 0.0) apply_relax = .TRUE. 
                        
                        if (f_grnd(i,j) .gt. 0.0 .and. &
                         (f_grnd(i-1,j) .eq. 0.0 .or. f_grnd(i+1,j) .eq. 0.0 &
                            .or. f_grnd(i,j-1) .eq. 0.0 .or. f_grnd(i,j+1) .eq. 0.0)) apply_relax = .TRUE. 
                
                    case(3)
                        ! Relax all points
                        
                        apply_relax = .TRUE. 
                
                    case DEFAULT
                        ! No relaxation

                        apply_relax = .FALSE.

                end select
                

                if (apply_relax) then 

                    if (tau .eq. 0.0) then
                        ! Impose ice thickness 

                        H_ice(i,j) = H_ref(i,j) 

                    else
                        ! Apply relaxation to reference state 

                        dHdt = (H_ref(i,j) - H_ice(i,j)) / tau 

                        H_ice(i,j) = H_ice(i,j) + dHdt*dt 

                    end if 
                end if 

            end do 
        end do 


        return 

    end subroutine relax_ice_thickness


end module mass_conservation
