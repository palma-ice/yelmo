module calving
    ! Definitions for various calving laws 

    use yelmo_defs, only : sp, dp, prec 

    implicit none 


    private 

    public :: apply_calving
    public :: calc_calving_rate_simple
    public :: calc_calving_rate_flux 
    public :: calc_calving_rate_eigen
    public :: calc_calving_rate_kill 
    
contains 

    subroutine apply_calving(H_ice,calv,f_grnd,H_min_flt,dt)

        implicit none 

        real(prec), intent(INOUT) :: H_ice(:,:) 
        real(prec), intent(INOUT) :: calv(:,:) 
        real(prec), intent(IN)    :: f_grnd(:,:) 
        real(prec), intent(IN)    :: H_min_flt
        real(prec), intent(IN)    :: dt  
        
        ! Ensure calving is limited to amount of available ice to calve  
        where(f_grnd .eq. 0.0 .and. (H_ice-dt*calv) .lt. 0.0) calv = H_ice/dt

        ! Additionally modify calving to remove any ice (margin or not)
        ! less than H_min_flt 
        where (f_grnd .eq. 0.0 .and. H_ice .lt. H_min_flt) calv = H_ice/dt
        
        ! Apply modified mass balance to update the ice thickness 
        H_ice = H_ice - dt*calv
        
        ! Also ensure tiny numeric ice thicknesses are removed
        where (f_grnd .eq. 0.0 .and. H_ice .lt. 1e-5) H_ice = 0.0 
        
        return 
        
    end subroutine apply_calving

    subroutine calc_calving_rate_simple(calv,H_ice,f_grnd,H_calv,tau)
        ! Calculate the calving rate [m/a] based on a simple threshold rule
        ! H_ice < H_calv

        implicit none 

        real(prec), intent(OUT) :: calv(:,:)
        real(prec), intent(IN)  :: H_ice(:,:)                ! [m] Ice thickness 
        real(prec), intent(IN)  :: f_grnd(:,:)               ! [-] Grounded fraction
        real(prec), intent(IN)  :: H_calv 
        real(prec), intent(IN)  :: tau                       ! [a] Calving timescale, ~ 1yr

        ! Local variables 
        integer :: i, j, nx, ny
        real(prec) :: H_mrgn
        logical :: is_front 
        integer :: n_ocean 

        nx = size(H_ice,1)
        ny = size(H_ice,2)
            
            
        ! Initially set calving rate to zero 
        calv = 0.0 

        do j=2,ny-1
        do i=2,nx-1

            if ( (f_grnd(i,j) .eq. 0.0 .and. H_ice(i,j) .gt. 0.0) .and. &
                   ( (f_grnd(i-1,j) .eq. 0.0 .and. H_ice(i-1,j).eq.0.0) .or. &
                     (f_grnd(i+1,j) .eq. 0.0 .and. H_ice(i+1,j).eq.0.0) .or. &
                     (f_grnd(i,j-1) .eq. 0.0 .and. H_ice(i,j-1).eq.0.0) .or. &
                     (f_grnd(i,j+1) .eq. 0.0 .and. H_ice(i,j+1).eq.0.0) ) ) then 
                ! Ice-shelf floating margin: floating ice point with open ocean neighbor 
                ! If this point is an ice front, check for calving

                if (H_ice(i,j) .lt. H_calv) then 
                    ! Apply calving at front, delete all ice in point (H_ice) 

                    calv(i,j) = (H_calv-H_ice(i,j)) / tau 

                end if 

            end if

        end do
        end do

        return 

    end subroutine calc_calving_rate_simple
    
    subroutine calc_calving_rate_flux(calv,H_ice,f_grnd,mbal,ux,uy,dx,H_calv,tau)
        ! Calculate the calving rate [m/a] based on a simple threshold rule
        ! H_ice < H_calv

        implicit none 

        real(prec), intent(OUT) :: calv(:,:)
        real(prec), intent(IN)  :: H_ice(:,:)                ! [m] Ice thickness 
        real(prec), intent(IN)  :: f_grnd(:,:)               ! [-] Grounded fraction
        real(prec), intent(IN)  :: mbal(:,:)                 ! [m/a] Net mass balance 
        real(prec), intent(IN)  :: ux(:,:)                   ! [m/a] velocity, x-direction (ac-nodes)
        real(prec), intent(IN)  :: uy(:,:)                   ! [m/a] velocity, y-direction (ac-nodes)
        real(prec), intent(IN)  :: dx 
        real(prec), intent(IN)  :: H_calv                    ! [m] Threshold for calving
        real(prec), intent(IN)  :: tau                       ! [a] Calving timescale, ~ 1yr

        ! Local variables 
        integer :: i, j, nx, ny
        real(prec) :: eps_xx, eps_yy  
        logical :: test_mij, test_pij, test_imj, test_ipj
        logical :: positive_mb 
        real(prec), allocatable :: dHdt(:,:), H_diff(:,:)  

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        allocate(dHdt(nx,ny))
        allocate(H_diff(nx,ny))

        ! Ice thickness above threshold
        H_diff = H_ice - H_calv

        ! Diagnosed lagrangian rate of change
        dHdt = 0.0 

        do j = 2, ny
        do i = 2, nx
        
                ! Calculate strain rate locally (aa-node)
                eps_xx = (ux(i,j) - ux(i-1,j))/dx
                eps_yy = (uy(i,j) - uy(i,j-1))/dx

                ! Calculate thickness change via conservation
                dHdt(i,j) = mbal(i,j) - H_ice(i,j)*(eps_xx+eps_yy)

        end do 
        end do
        

        ! Initially set calving rate to zero 
        calv = 0.0 

        do j = 2, ny-1
        do i = 2, nx-1

            if ( (f_grnd(i,j) .eq. 0.0 .and. H_ice(i,j) .gt. 0.0 .and. H_diff(i,j).lt.0.0) .and. &
                   ( (f_grnd(i-1,j) .eq. 0.0 .and. H_ice(i-1,j).eq.0.0) .or. &
                     (f_grnd(i+1,j) .eq. 0.0 .and. H_ice(i+1,j).eq.0.0) .or. &
                     (f_grnd(i,j-1) .eq. 0.0 .and. H_ice(i,j-1).eq.0.0) .or. &
                     (f_grnd(i,j+1) .eq. 0.0 .and. H_ice(i,j+1).eq.0.0) ) ) then 
                ! Ice-shelf floating margin: floating ice point with open ocean neighbor 
                 
                ! Check if current point is at the floating ice front,
                ! and has thickness less than threshold, or if
                ! ice below H_calv limit, accounting for mass flux from inland

                positive_mb = (mbal(i,j).gt.0.0)

                test_mij = ( ((H_diff(i-1,j).gt.0.0).and.(ux(i-1,j).gt.0.0)  &  ! neighbor (i-1,j) total > H_calv
                    .and.  (dHdt(i-1,j).gt.(-H_diff(i-1,j)*abs(ux(i-1,j)/dx)))) & 
                    .or.(f_grnd(i-1,j).gt.0.0.and.positive_mb ))

                test_pij = ( ((H_diff(i+1,j).gt.0.0).and.(ux(i,j).lt.0.0) & ! neighbor (i+1,j) total > H_calv
                    .and.(dHdt(i+1,j).gt.(-H_diff(i+1,j)*abs(ux(i,j)/dx)))) &
                    .or.(f_grnd(i+1,j).gt.0.0.and.positive_mb ))

                test_imj = ( ((H_diff(i,j-1).gt.0.0).and.(uy(i,j-1).gt.0.0)  &  ! neighbor (i,j-1) total > H_calv
                    .and.(dHdt(i,j-1).gt.(-H_diff(i,j-1)*abs(uy(i,j-1)/dx))))&
                    .or.(f_grnd(i,j-1).gt.0.0.and.positive_mb ))

                test_ipj = ( ((H_diff(i,j+1).gt.0.0).and.(uy(i,j).lt.0.0) & ! neighbor (i,j+1) total > H_calv
                    .and.(dHdt(i,j+1).gt.(-H_diff(i,j+1)*abs(uy(i,j)/dx))))&
                    .or.(f_grnd(i,j+1).gt.0.0.and.positive_mb ))

                if ((.not.(test_mij.or.test_pij.or.test_imj.or.test_ipj))) then
                    calv(i,j) = (H_calv - H_ice(i,j)) / tau             
                end if  

            end if

        end do
        end do

        return 

    end subroutine calc_calving_rate_flux
    
    subroutine calc_calving_rate_eigen(calv,H_ice,f_grnd,ux_bar,uy_bar,dx,dy,H_calv,k_calv)
        ! Calculate the calving rate [m/a] based on the "eigencalving" law
        ! from Levermann et al. (2012)

        implicit none 

        real(prec), intent(OUT) :: calv(:,:)
        real(prec), intent(IN)  :: H_ice(:,:)
        real(prec), intent(IN)  :: f_grnd(:,:)  
        real(prec), intent(IN)  :: ux_bar(:,:), uy_bar(:,:)
        real(prec), intent(IN)  :: dx, dy 
        real(prec), intent(IN)  :: H_calv 
        real(prec), intent(IN)  :: k_calv

        ! Local variables 
        integer :: i, j, nx, ny
        real(prec), allocatable :: eps_xx(:,:), eps_yy(:,:) 
        real(prec) :: eps_xx_max, eps_yy_max 
        integer :: imax_xx, jmax_xx, imax_yy, jmax_yy                                          
        real(prec), allocatable :: spr(:,:)   ! total spreading rate in the two directions epsxx * epsyyy 

        logical, allocatable :: is_front(:,:)   
        integer :: n_ocean 

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        allocate(spr(nx,ny),eps_xx(nx,ny),eps_yy(nx,ny))
        allocate(is_front(nx,ny))

        do j = 1, ny
        do i = 1, nx
            ! Calculate strain rate locally (Aa node)
            eps_xx(i,j) = (ux_bar(i,j) - ux_bar(i-1,j))/dx
            eps_yy(i,j) = (uy_bar(i,j) - uy_bar(i,j-1))/dy            
        end do
        end do
        

        calv = 0.0

        do j = 2, ny-1
        do i = 2, nx-1  
            
            if ( (f_grnd(i,j) .eq. 0.0 .and. H_ice(i,j) .gt. 0.0) .and. &
                   ( (f_grnd(i-1,j) .eq. 0.0 .and. H_ice(i-1,j).eq.0.0) .or. &
                     (f_grnd(i+1,j) .eq. 0.0 .and. H_ice(i+1,j).eq.0.0) .or. &
                     (f_grnd(i,j-1) .eq. 0.0 .and. H_ice(i,j-1).eq.0.0) .or. &
                     (f_grnd(i,j+1) .eq. 0.0 .and. H_ice(i,j+1).eq.0.0) ) ) then 
                ! Ice-shelf floating margin: floating ice point with open ocean neighbor 
                ! If this point is an ice front, check for calving

                if ((eps_xx(i,j).gt.0.0).and.(eps_yy(i,j).gt.0.0)) then                   
                    ! Divergence in both directions, apply calving law 
                    ! Flux condition + calving rate with spreading:       

                    calv(i,j) = k_calv * eps_xx(i,j)*eps_yy(i,j)                                       
                
                end if

            end if

        end do
        end do

        return 

    end subroutine calc_calving_rate_eigen

    subroutine calc_calving_rate_kill(calv,H_ice,f_grnd,tau)

        implicit none 

        real(prec), intent(OUT) :: calv(:,:)
        real(prec), intent(IN)  :: H_ice(:,:)
        real(prec), intent(IN)  :: f_grnd(:,:) 
        real(prec), intent(IN)  :: tau 

        ! Kill all floating ice using a characteristic timescale tau
        where (f_grnd .eq. 0.0) calv = H_ice / tau 

        return 

    end subroutine calc_calving_rate_kill

end module calving 
