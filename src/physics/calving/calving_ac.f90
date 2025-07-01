module calving_ac
    ! Definitions for various calving laws on ac nodes. 
    ! This will be used for flux calving (lsf) 

    use yelmo_defs, only : sp, dp, wp, prec, TOL_UNDERFLOW
    use yelmo_tools, only : get_neighbor_indices
    use topography, only : calc_H_eff 

    implicit none 
    private 
    
    ! === Calving related stress/strain routines === 
    !public :: calc_eps_eff
    !public :: calc_tau_eff_ac

    ! === Floating calving routines === 
    public :: calc_calving_threshold_lsf
    
    ! === Grounded calving routines === 

    ! === CalvMIP calving rates ===
    public :: calvmip_exp1
    public :: calvmip_exp2
    public :: calvmip_exp5_ac
    public :: calvmip_exp5_aa

contains 

    ! ===================================================================
    !
    ! Stress and strain rates
    !
    ! ===================================================================

    ! TO DO
    subroutine calc_eps_eff_ac(eps_eff,eps_eig_1,eps_eig_2,f_ice,boundaries)
        ! Effective strain rate based. Levermann et al. (2012)

        implicit none 

        real(wp), intent(OUT) :: eps_eff(:,:)
        real(wp), intent(IN)  :: eps_eig_1(:,:)
        real(wp), intent(IN)  :: eps_eig_2(:,:) 
        real(wp), intent(IN)  :: f_ice(:,:) 
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer  :: i, j, nx, ny, n  
        integer  :: im1, jm1, ip1, jp1 
        real(wp) :: eps_eff_neighb(4)

        nx = size(eps_eff,1)
        ny = size(eps_eff,2) 

        !!$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,n,eps_eff_neighb)
        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
            
            if (f_ice(i,j) .eq. 0.0_wp) then 
                ! Ice-free point, no strain

                eps_eff(i,j) = 0.0_wp 

            else if (eps_eig_1(i,j) .eq. 0.0 .and. eps_eig_2(i,j) .eq. 0.0) then 
                ! Margin point was likely just advected, no stresses available, 
                ! use maximum value of eps_eff from upstream neighbors.

                eps_eff_neighb = 0.0_wp 

                if (f_ice(im1,j).gt.0.0) eps_eff_neighb(1) = eps_eig_1(im1,j) * eps_eig_2(im1,j)
                if (f_ice(ip1,j).gt.0.0) eps_eff_neighb(2) = eps_eig_1(ip1,j) * eps_eig_2(ip1,j)
                if (f_ice(i,jm1).gt.0.0) eps_eff_neighb(3) = eps_eig_1(i,jm1) * eps_eig_2(i,jm1)
                if (f_ice(i,jp1).gt.0.0) eps_eff_neighb(4) = eps_eig_1(i,jp1) * eps_eig_2(i,jp1)

                n = count(eps_eff_neighb.ne.0.0_wp)

                if (n .gt. 0) then 
                    eps_eff(i,j) = sum(eps_eff_neighb,mask=eps_eff_neighb.ne.0.0_wp) / real(n,wp)
                else 
                    eps_eff(i,j) = 0.0_wp
                end if 

            else 
                ! Stresses are available at this margin point. 
                ! Calculate the effective strain rate directly.
                eps_eff(i,j) = eps_eig_1(i,j) * eps_eig_2(i,j)
                
            end if 
            
        end do 
        end do 
        !!$omp end parallel do

        return 

    end subroutine calc_eps_eff_ac
    
    subroutine calc_tau_eff_ac(tau_eff,tau_eig_1,tau_eig_2,f_ice,w2,boundaries)
        ! Effective stress rates. Based on additional of principal stresses.

        implicit none 

        real(wp), intent(OUT) :: tau_eff(:,:)
        real(wp), intent(IN)  :: tau_eig_1(:,:)
        real(wp), intent(IN)  :: tau_eig_2(:,:) 
        real(wp), intent(IN)  :: f_ice(:,:) 
        real(wp), intent(IN)  :: w2 
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer  :: i, j, nx, ny, n  
        integer  :: im1, jm1, ip1, jp1 
        real(wp) :: tau_eff_neighb(4) 

        nx = size(tau_eff,1)
        ny = size(tau_eff,2) 
        
        !!$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,n,tau_eff_neighb)
        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
            
            if (f_ice(i,j) .eq. 0.0_wp) then 
                ! Ice-free point, no stress

                tau_eff(i,j) = 0.0_wp 

            else if (tau_eig_1(i,j) .eq. 0.0 .and. tau_eig_2(i,j) .eq. 0.0) then 
                ! Margin point was likely just advected, no stresses available, 
                ! use maximum value of tau_eff from upstream neighbors.

                tau_eff_neighb = 0.0_wp 

                if (f_ice(im1,j).gt.0.0) tau_eff_neighb(1) = calc_tau_eff_now_ac(tau_eig_1(im1,j),tau_eig_2(im1,j),w2)
                if (f_ice(ip1,j).gt.0.0) tau_eff_neighb(2) = calc_tau_eff_now_ac(tau_eig_1(ip1,j),tau_eig_2(ip1,j),w2)
                if (f_ice(i,jm1).gt.0.0) tau_eff_neighb(3) = calc_tau_eff_now_ac(tau_eig_1(i,jm1),tau_eig_2(i,jm1),w2)
                if (f_ice(i,jp1).gt.0.0) tau_eff_neighb(4) = calc_tau_eff_now_ac(tau_eig_1(i,jp1),tau_eig_2(i,jp1),w2)

                n = count(tau_eff_neighb.ne.0.0_wp)

                if (n .gt. 0) then 
                    tau_eff(i,j) = sum(tau_eff_neighb,mask=tau_eff_neighb.ne.0.0_wp) / real(n,wp)
                else 
                    tau_eff(i,j) = 0.0_wp 
                end if 

            else 
                ! Stresses are available at this margin point. 
                ! Calculate the effective strain rate directly.

                tau_eff(i,j) = calc_tau_eff_now_ac(tau_eig_1(i,j),tau_eig_2(i,j),w2)
            
            end if 

        end do 
        end do
        !!$omp end parallel do

        return 

    end subroutine calc_tau_eff_ac
    
    elemental function calc_tau_eff_now_ac(teig1,teig2,w2) result(tau_eff) 

        implicit none 

        real(wp), intent(IN) :: teig1 
        real(wp), intent(IN) :: teig2
        real(wp), intent(IN) :: w2
        real(wp) :: tau_eff

        ! Local variables 
        real(wp) :: tau1, tau2

        tau1    = max(teig1,0.0_wp)
        tau2    = max(teig2,0.0_wp)
        tau_eff = sqrt(tau1**2 + (w2 * tau2)**2)

        return 

    end function calc_tau_eff_now_ac

    ! ===================================================================
    !
    ! Calving - floating ice 
    !
    ! ===================================================================
    
    subroutine calc_calving_threshold_lsf(cr_acx,cr_acy,u_acx,v_acy,H_ice,H_ice_c,boundaries)
        ! Threshold calving rate flux based on CalvingMIP experiment 5.
        ! Valid for floating and grounded ice.
            
        implicit none
            
        real(wp), intent(OUT) :: cr_acx(:,:), cr_acy(:,:)
        real(wp), intent(IN)  :: u_acx(:,:),  v_acy(:,:)
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: H_ice_c
        character(len=*), intent(IN)  :: boundaries             ! Boundary conditions to impose
                
        ! local variables
        integer  :: i, j, ip1, im1, jp1, jm1, nx, ny
        real(wp) :: uxy_aa,uxy_acx,uxy_acy,u_acy,v_acx
        real(wp), allocatable :: H_ice_fill(:,:), wv_aa(:,:) 
            
        nx = size(u_acx,1)
        ny = size(u_acx,2) 
        allocate(H_ice_fill(nx,ny))
        allocate(wv_aa(nx,ny))
    
        ! Initialize    
        uxy_aa     = 0.0_wp
        H_ice_fill = H_ice
        wv_aa      = 0.0_wp
                
        ! since we compute on ac-nodes and ice thickness are on aa-nodes
        ! we need to extrapolate ice thickness to the ocean
        call extrapolate_ocn_laplace_simple(H_ice_fill,H_ice,H_ice)
    
        do j = 1, ny
            do i = 1, nx
                call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
                ! velocity on aa-node
                uxy_aa     = ((0.5*(u_acx(i,j)+u_acx(im1,j)))**2 + (0.5*(v_acy(i,j)+v_acy(i,jm1)))**2)**0.5
                wv_aa(i,j) = MAX(0.0_wp,1.0_wp+(H_ice_c-H_ice_fill(i,j))/H_ice_c)*uxy_aa                
            end do
        end do

        do j = 1, ny
            do i = 1, nx
                ! Stagger velocities x/y ac-velocities into y/x ac-nodes
                call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
                u_acy = 0.25_wp*(u_acx(i,j)+u_acx(im1,j)+u_acx(im1,jp1)+u_acx(i,jp1))
                v_acx = 0.25_wp*(v_acy(i,j)+v_acy(i,jm1)+v_acy(ip1,jm1)+v_acy(ip1,j))
                ! x-direction
                uxy_acx     = MAX(1e-8,(u_acx(i,j)**2 + v_acx**2)**0.5)
                cr_acx(i,j) = -(u_acx(i,j)/uxy_acx)*0.5*(wv_aa(i,j)+wv_aa(ip1,j))
                ! y-direction
                uxy_acy     = MAX(1e-8,(v_acy(i,j)**2 + u_acy**2)**0.5)
                cr_acy(i,j) = -(v_acy(i,j)/uxy_acy)*0.5*(wv_aa(i,j)+wv_aa(i,jp1))
            end do
        end do

        deallocate(H_ice_fill)
        deallocate(wv_aa)
    
        return
        
    end subroutine calc_calving_threshold_lsf

    subroutine calc_calving_rate_vonmises_l19(mb_calv,H_ice,f_ice,f_grnd,tau_eff,dx,kt,boundaries)
        ! Calculate the 'horizontal' calving rate [m/yr] based on the 
        ! von Mises stress approach, as outlined by Lipscomb et al. (2019)
        ! Eqs. 73-75.
        ! L19: kt = 0.0025 m yr-1 Pa-1, w2=25
    
        implicit none 
    
        real(wp), intent(OUT) :: mb_calv(:,:)
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: f_grnd(:,:)  
        real(wp), intent(IN)  :: tau_eff(:,:)
        real(wp), intent(IN)  :: dx
        real(wp), intent(IN)  :: kt(:,:)
        character(len=*), intent(IN) :: boundaries 
    
        ! Local variables 
        integer  :: i, j
        integer  :: im1, jm1, ip1, jp1
        integer  :: nx, ny 
        real(wp) :: dy   
        real(wp) :: calv_ref
        real(wp) :: calv_now
        real(wp) :: H_eff 
        real(wp) :: wt
            
        real(wp), parameter :: calv_lim = 2000.0_wp     ! To avoid really high calving values
    
        nx = size(H_ice,1)
        ny = size(H_ice,2)
    
        ! Assume square grid cells 
        dy = dx 
    
        mb_calv = 0.0_wp
    
            !!$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,wt,calv_ref,H_eff,calv_now)
        do j = 1, ny
        do i = 1, nx  

            ! Calculate lateral calving rate 
            calv_ref = kt(i,j)*tau_eff(i,j) 
    
            calv_now = min(calv_now,calv_lim)
                        
            ! Get calving mass balance rate
            mb_calv(i,j) = -calv_now

        end do
        end do
        !!$omp end parallel do
    
        return 
    
    end subroutine calc_calving_rate_vonmises_l19

    subroutine calc_calving_rate_vonmises_v23(mb_calv,tau_eff,tau_ice,u_acx,v_acy)
        ! Calculate the 'horizontal' calving rate [m/yr] based on the 
        ! von Mises stress approach, as outlined by Wilner et al. (2023)
        ! Eq. 2.
        ! c = v*tau_eff/tau_ice

        implicit none 

        real(wp), intent(OUT) :: mb_calv(:,:) 
        real(wp), intent(IN)  :: tau_eff(:,:)
        real(wp), intent(IN)  :: tau_ice
        real(wp), intent(IN)  :: u_acx(:,:)
        real(wp), intent(IN)  :: v_acy(:,:)

        ! Local variables 
        integer  :: i, j
        real(wp) :: calv_ref
        real(wp) :: calv_now
        
        real(wp), parameter :: calv_lim = 2000.0_wp     ! To avoid really high calving values

        ! Initializa calving 
        mb_calv = 0.0_wp

        !!$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,wt,calv_ref,H_eff,calv_now)
        do j = 1, size(tau_eff,2)
        do i = 1, size(tau_eff,1)  
            
            ! Calculate lateral calving rate
            ! jablasco: to do! 
            calv_ref = u_acx(i,j)*tau_eff(i,j)/tau_ice 

            ! Apply calving limit
            calv_now = min(calv_ref,calv_lim)
                    
            ! Get calving mass balance rate
            mb_calv(i,j) = -calv_now

        end do
        end do
        !!$omp end parallel do

        return 

    end subroutine calc_calving_rate_vonmises_v23
       
    subroutine calc_calving_rate_eigen(mb_calv,H_ice,f_ice,f_grnd,eps_eff,dx,k2,boundaries)
        ! Calculate the 'horizontal' calving rate [m/yr] based on the 
        ! von Mises stress approach, as outlined by Lipscomb et al. (2019)
        ! Eqs. 73-75.
        ! L19: kt = 0.0025 m yr-1 Pa-1, w2=25

        implicit none 

        real(wp), intent(OUT) :: mb_calv(:,:)
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: f_ice(:,:)
        real(wp), intent(IN)  :: f_grnd(:,:)  
        real(wp), intent(IN)  :: eps_eff(:,:)
        real(wp), intent(IN)  :: dx
        real(wp), intent(IN)  :: k2
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer  :: i, j
        integer  :: im1, jm1, ip1, jp1
        integer  :: nx, ny
        integer  :: n_ocean 
        logical  :: is_margin 
        real(wp) :: dy   
        real(wp) :: calv_ref
        real(wp) :: calv_now
        real(wp) :: H_eff 

        real(wp) :: dxx, dyy, dxy 
        real(wp) :: eps_eig_1_now, eps_eig_2_now
        real(wp) :: eps_eff_neighb(4)
        real(wp) :: wt
        
        real(wp), parameter :: calv_lim = 2000.0_wp     ! To avoid really high calving values

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        ! Assume square grid cells 
        dy = dx 

        mb_calv = 0.0_wp

        !!$omp parallel do collapse(2) private(i,j,im1,ip1,jm1,jp1,wt,calv_ref,H_eff,calv_now)
        do j = 1, ny
        do i = 1, nx  
            ! Calculate lateral calving rate 
            calv_ref = max( k2*eps_eff(i,j), 0.0_wp )

            ! Apply calving limit
            calv_now = min(calv_now,calv_lim)

            ! Get calving mass balance rate
            mb_calv(i,j) = -calv_now

        end do
        end do
        !!$omp end parallel do

        return 

    end subroutine calc_calving_rate_eigen
     
    ! ===================================================================
    !
    ! Calving - grounded ice 
    !
    ! ===================================================================

    
    ! ===================================================================
    !
    !                      CalvMIP experiments
    !
    ! ===================================================================

    subroutine calvmip_exp1(cr_acx,cr_acy,u_acx,v_acy,lsf_aa,dx,boundaries)
        ! Experiment 1 & 3 of CalvMIP
        implicit none
    
        real(wp), intent(OUT) :: cr_acx(:,:),cr_acy(:,:)   ! Calving rates on ac-nodes
        real(wp), intent(IN)  :: u_acx(:,:),v_acy(:,:)     ! Velocities on ac-nodes
        real(wp), intent(IN)  :: lsf_aa(:,:)               ! LSF mask on aa-nodes
        real(wp), intent(IN)  :: dx                        ! Ice resolution
        character(len=*), intent(IN)  :: boundaries        ! Boundary conditions to impose
    
        ! Local variables
        integer  :: i, j, im1, ip1, jm1, jp1, nx, ny
        real(wp) :: r, rip1, rim1, rjp1, rjm1
    
        nx = size(u_acx,1)
        ny = size(u_acx,2)
    
        ! Initialize calving rates to opposite as velocity
        cr_acx = -u_acx 
        cr_acy = -v_acy
    
        r    = 0.0_wp
        rip1 = 0.0_wp
        rim1 = 0.0_wp
        rjp1 = 0.0_wp
        rjm1 = 0.0_wp
        
        do j = 1, ny
        do i = 1, nx
    
            ! Below radius, no calving.
            ! aa-nodes indices
            r = sqrt((0.5*(nx+1)-i)*(0.5*(nx+1)-i) + (0.5*(ny+1)-j)*(0.5*(ny+1)-j))*dx

            ! Below radius
            if (r .lt. 750e3) then
                ! Now treat border points
                call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
                rip1 = sqrt((0.5*(nx+1)-ip1)*(0.5*(nx+1)-ip1) + (0.5*(ny+1)-j)*(0.5*(ny+1)-j))*dx
                rim1 = sqrt((0.5*(nx+1)-im1)*(0.5*(nx+1)-im1) + (0.5*(ny+1)-j)*(0.5*(ny+1)-j))*dx
                rjp1 = sqrt((0.5*(nx+1)-i)*(0.5*(nx+1)-i) + (0.5*(ny+1)-jp1)*(0.5*(ny+1)-jp1))*dx
                rjm1 = sqrt((0.5*(nx+1)-i)*(0.5*(nx+1)-i) + (0.5*(ny+1)-jm1)*(0.5*(ny+1)-jm1))*dx
        
                ! === Check direction ===
                ! x-direction
                if (rip1 .ge. 750e3) then
                    ! border point right
                    cr_acx(i,j)   = -u_acx(i,j)
                else
                    cr_acx(i,j)   = 0.0_wp
                end if

                if (rim1 .ge. 750e3) then
                    ! border point left
                    cr_acx(im1,j) = -u_acx(im1,j)
                else
                    cr_acx(im1,j) = 0.0_wp
                end if

                ! y-direction
                if (rjp1 .ge. 750e3) then
                    ! border point top
                    cr_acy(i,j)   = -v_acy(i,j)
                else
                    !border point top
                    cr_acy(i,j)   = 0.0_wp
                end if

                if (rjm1 .ge. 750e3) then
                    ! border point bottom
                    cr_acy(i,jm1) = -v_acy(i,jm1)
                else
                    ! border point bottom
                    cr_acy(i,jm1) = 0.0_wp
                end if
                
            end if
    
        end do
        end do
    
        return
    
    end subroutine calvmip_exp1
    
    subroutine calvmip_exp2(cr_acx,cr_acy,u_acx,v_acy,time,boundaries)
        ! Experiment 2 & 4 of CalvMIP
    
        implicit none
    
        real(wp), intent(OUT) :: cr_acx(:,:), cr_acy(:,:)
        real(wp), intent(IN)  :: u_acx(:,:),  v_acy(:,:)
        real(wp), intent(IN)  :: time
        character(len=*), intent(IN)  :: boundaries             ! Boundary conditions to impose
    
        ! local variables
        integer  :: i, j, ip1, im1, jp1, jm1, nx, ny
        real(wp) :: wv,uxy_acx,uxy_acy,u_acy,v_acx   
        real(wp), parameter :: pi = acos(-1.0)  ! Calculate pi intrinsically

        nx = size(u_acx,1)
        ny = size(u_acx,2) 

        ! Initialize    
        wv      = -300.0 * sin(2.0 * pi * time / 1000.0) 
        uxy_acx = 0.0_wp
        uxy_acy = 0.0_wp
        u_acy   = 0.0_wp
        v_acx   = 0.0_wp
        cr_acx  = 0.0_wp
        cr_acy  = 0.0_wp
    
        do j = 1, ny
        do i = 1, nx
            ! Stagger velocities x/y ac-velocities into y/x ac-nodes
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
            u_acy = 0.25_wp*(u_acx(i,j)+u_acx(im1,j)+u_acx(im1,jp1)+u_acx(i,jp1))
            v_acx = 0.25_wp*(v_acy(i,j)+v_acy(i,jm1)+v_acy(ip1,jm1)+v_acy(ip1,j))
            ! x-direction
            uxy_acx     = MAX(1e-8,(u_acx(i,j)**2 + v_acx**2)**0.5)
            cr_acx(i,j) = -u_acx(i,j)+(u_acx(i,j)/uxy_acx)*wv
            ! y-direction
            uxy_acy     = MAX(1e-8,(v_acy(i,j)**2 + u_acy**2)**0.5)
            cr_acy(i,j) = -v_acy(i,j)+(v_acy(i,j)/uxy_acy)*wv
        end do
        end do

        return
    
    end subroutine calvmip_exp2

    subroutine calvmip_exp5_ac(cr_acx,cr_acy,u_acx,v_acy,H_ice,H_ice_c,boundaries)
        ! Experiment 5 of CalvMIP
        
        implicit none
        
        real(wp), intent(OUT) :: cr_acx(:,:), cr_acy(:,:)
        real(wp), intent(IN)  :: u_acx(:,:),  v_acy(:,:)
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: H_ice_c
        character(len=*), intent(IN)  :: boundaries             ! Boundary conditions to impose
            
        ! local variables
        integer  :: i, j, ip1, im1, jp1, jm1, nx, ny
        real(wp) :: wv_acx,wv_acy,wvr,uxy_acx,uxy_acy,u_acy,v_acx 
        real(wp), allocatable :: H_ice_fill(:,:)  
        
        nx = size(u_acx,1)
        ny = size(u_acx,2) 
        allocate(H_ice_fill(nx,ny))

        ! Initialize    
        wv_acx  = 0.0_wp
        wv_acy  = 0.0_wp
        uxy_acx = 0.0_wp
        uxy_acy = 0.0_wp
        u_acy   = 0.0_wp
        v_acx   = 0.0_wp
        cr_acx  = 0.0_wp
        cr_acy  = 0.0_wp
        H_ice_fill = H_ice
            
        ! since we compute on ac-nodes and ice thickness are on aa-nodes
        ! we need to extrapolate ice thickness to the ocean
        call extrapolate_ocn_laplace_simple(H_ice_fill,H_ice,H_ice)

        do j = 1, ny
            do i = 1, nx
                call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
                ! x-direction
                wv_acx = MAX(0.0_wp,1.0_wp+(H_ice_c-(0.5*(H_ice_fill(i,j)+H_ice_fill(ip1,j))))/H_ice_c)
                cr_acx(i,j) = -u_acx(i,j)*wv_acx
                ! y-direction
                wv_acy = MAX(0.0_wp,1.0_wp+(H_ice_c-(0.5*(H_ice_fill(i,j)+H_ice_fill(i,jp1))))/H_ice_c)
                cr_acy(i,j) = -v_acy(i,j)*wv_acy
            end do
        end do

        deallocate(H_ice_fill)

        return
    
    end subroutine calvmip_exp5_ac

    subroutine calvmip_exp5_aa(cr_acx,cr_acy,u_acx,v_acy,H_ice,H_ice_c,boundaries)
        ! Experiment 5 of CalvMIP
            
        implicit none
            
        real(wp), intent(OUT) :: cr_acx(:,:), cr_acy(:,:)
        real(wp), intent(IN)  :: u_acx(:,:),  v_acy(:,:)
        real(wp), intent(IN)  :: H_ice(:,:)
        real(wp), intent(IN)  :: H_ice_c
        character(len=*), intent(IN)  :: boundaries             ! Boundary conditions to impose
                
        ! local variables
        integer  :: i, j, ip1, im1, jp1, jm1, nx, ny
        real(wp) :: uxy_aa,uxy_acx,uxy_acy,u_acy,v_acx
        real(wp), allocatable :: H_ice_fill(:,:), wv_aa(:,:) 
            
        nx = size(u_acx,1)
        ny = size(u_acx,2) 
        allocate(H_ice_fill(nx,ny))
        allocate(wv_aa(nx,ny))
    
        ! Initialize    
        uxy_aa     = 0.0_wp
        H_ice_fill = H_ice
        wv_aa      = 0.0_wp
                
        ! since we compute on ac-nodes and ice thickness are on aa-nodes
        ! we need to extrapolate ice thickness to the ocean
        call extrapolate_ocn_laplace_simple(H_ice_fill,H_ice,H_ice)
    
        do j = 1, ny
            do i = 1, nx
                call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
                ! velocity on aa-node
                uxy_aa     = ((0.5*(u_acx(i,j)+u_acx(im1,j)))**2 + (0.5*(v_acy(i,j)+v_acy(i,jm1)))**2)**0.5
                wv_aa(i,j) = MAX(0.0_wp,1.0_wp+(H_ice_c-H_ice_fill(i,j))/H_ice_c)*uxy_aa                
            end do
        end do

        do j = 1, ny
            do i = 1, nx
                ! Stagger velocities x/y ac-velocities into y/x ac-nodes
                call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
                u_acy = 0.25_wp*(u_acx(i,j)+u_acx(im1,j)+u_acx(im1,jp1)+u_acx(i,jp1))
                v_acx = 0.25_wp*(v_acy(i,j)+v_acy(i,jm1)+v_acy(ip1,jm1)+v_acy(ip1,j))
                ! x-direction
                uxy_acx     = MAX(1e-8,(u_acx(i,j)**2 + v_acx**2)**0.5)
                cr_acx(i,j) = -(u_acx(i,j)/uxy_acx)*0.5*(wv_aa(i,j)+wv_aa(ip1,j))
                ! y-direction
                uxy_acy     = MAX(1e-8,(v_acy(i,j)**2 + u_acy**2)**0.5)
                cr_acy(i,j) = -(v_acy(i,j)/uxy_acy)*0.5*(wv_aa(i,j)+wv_aa(i,jp1))
            end do
        end do

        deallocate(H_ice_fill)
        deallocate(wv_aa)
    
        return
        
    end subroutine calvmip_exp5_aa
    
    ! ===================================================================
    !
    !                 Ocean extrapolation routines
    !
    ! ===================================================================

    subroutine extrapolate_ocn_laplace_simple(mask_fill, mask_orig,mask)
        ! Routine to extrapolate values using the Laplace equation.
        ! Assumes that value 0 in mask represents ice-free points
                
        implicit none
            
        real(wp), intent(INOUT) :: mask_fill(:,:)
        real(wp), intent(IN) :: mask_orig(:,:)
        real(wp), intent(IN) :: mask(:,:)
                
        ! Local variables
        integer :: i, j, iter
        real(wp) :: error, tol
        real(wp), allocatable :: mask_new(:,:)
                
        ! Allocate memory for the temporary array
        allocate(mask_new(size(mask_orig,1), size(mask_orig,2)))
                
        ! Initialize variables
        mask_fill = mask_orig
        mask_new  = mask_orig
        tol       = 1e-2_wp      ! Tolerance for convergence
        error     = tol + 1.0_wp
        iter      = 0
                
        ! Jacobi iteration
        do while (error > tol)
            error = 0.0_wp
            iter = iter + 1
                
            do i = 2, size(mask_orig,1)-1
                do j = 2, size(mask_orig,2)-1
                    if (mask(i,j) .eq. 0.0_wp) then
                        mask_new(i,j) = 0.25_wp * (mask_fill(i+1,j) + mask_fill(i-1,j) + mask_fill(i,j+1) + mask_fill(i,j-1))
                        error = error + abs(mask_new(i,j) - mask_fill(i,j))
                    end if
                end do
            end do
            mask_fill = mask_new
        end do
                
        deallocate(mask_new)
        
        return
            
    end subroutine extrapolate_ocn_laplace_simple

end module calving_ac