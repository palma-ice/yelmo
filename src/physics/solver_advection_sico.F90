module solver_advection_sico

    use yelmo_defs, only : sp, dp, wp 
    use yelmo_tools, only : get_neighbor_indices
    use solver_linear

    implicit none 
    
    private
    public :: linear_solver_save_advection
    public :: linear_solver_matrix_advection_csr_2D

    ! To be deprecated...
    public :: calc_adv2D_expl_sico
    public :: calc_adv2D_impl_sico 
    public :: calc_advec_horizontal_point_expl

contains 

    subroutine linear_solver_save_advection(H,lgs)
        ! Extract ice thickness solution from lgs object. 

        implicit none 

        real(wp), intent(OUT) :: H(:,:)                 ! [X]
        type(linear_solver_class), intent(IN) :: lgs 

        ! Local variables 
        integer :: i, j, n, nr 

        do n = 1, lgs%nmax

            i = lgs%n2i((n+1)/2)
            j = lgs%n2j((n+1)/2)

            nr = n
            H(i,j) = lgs%x_value(nr)

        end do

        return

    end subroutine linear_solver_save_advection

    subroutine linear_solver_matrix_advection_csr_2D(lgs,H,ux,uy,F,mask,dx,dy,dt,boundaries)
        ! Define sparse matrices A*x=b in format 'compressed sparse row' (csr)
        ! for 2D advection equations with velocity components
        ! ux and uy defined on ac-nodes (right and top borders of i,j grid cell)
        ! and variable to be advected H defined on aa nodes.
        ! Store sparse matrices in linear_solver_class object 'lgs' for later use.

        implicit none 

        type(linear_solver_class), intent(INOUT) :: lgs
        real(wp), intent(INOUT)   :: H(:,:)         ! [X] Variable of interest (aa nodes)
        real(wp), intent(IN)      :: ux(:,:)        ! [m a-1] Horizontal velocity x-direction (ac nodes)
        real(wp), intent(IN)      :: uy(:,:)        ! [m a-1] Horizontal velocity y-direction (ac nodes)
        real(wp), intent(IN)      :: F(:,:)         ! [m a-1] Net source/sink terms (aa nodes)
        integer,  intent(IN)      :: mask(:,:)      ! Advection mask
        real(wp), intent(IN)      :: dx             ! [m] Horizontal step x-direction
        real(wp), intent(IN)      :: dy             ! [m] Horizontal step y-direction 
        real(wp), intent(IN)      :: dt             ! [a] Time step 
        character(len=*), intent(IN) :: boundaries 

        ! Local variables  
        integer  :: i, j, k, nx, ny
        integer  :: im1, ip1, jm1, jp1 
        integer  :: n, nr, nc
        real(wp) :: dt_darea
        logical  :: is_periodic

        real(wp), allocatable  :: ux_1(:,:), ux_2(:,:)
        real(wp), allocatable  :: uy_1(:,:), uy_2(:,:)
        real(wp), allocatable  :: Hx_1(:,:), Hx_2(:,:)
        real(wp), allocatable  :: Hy_1(:,:), Hy_2(:,:)

        real(wp), parameter :: WOVI = 1.0     ! Weighing parameter for the over-implicit scheme 

        nx = size(H,1)
        ny = size(H,2) 

        dt_darea = dt/(dx*dy)

        ! Check if domain is periodic
        if (trim(boundaries) .eq. "periodic") then 
            is_periodic = .TRUE. 
        else 
            is_periodic = .FALSE. 
        end if 

        ! Safety check for initialization
        if (.not. allocated(lgs%x_value)) then 
            ! Object 'lgs' has not been initialized yet, do so now.

            call linear_solver_init(lgs,nx,ny,nvar=1,n_terms=5)

        end if

        ! Allocate local variables
        allocate(ux_1(nx,ny))
        allocate(ux_2(nx,ny))
        allocate(uy_1(nx,ny))
        allocate(uy_2(nx,ny))
        
        allocate(Hx_1(nx,ny))
        allocate(Hx_2(nx,ny))
        allocate(Hy_1(nx,ny))
        allocate(Hy_2(nx,ny))
        
        ! ====================================================================================
        ! Step 1: populate variables representing velocity components (ux,uy) at each
        ! cell face (left,right,bottom,top) and upstream variable (Hx,Hy) at each cell face. 

        ux_1 = 0.0
        ux_2 = 0.0
        uy_1 = 0.0
        uy_2 = 0.0
        
        Hx_1 = 0.0
        Hx_2 = 0.0
        Hy_1 = 0.0
        Hy_2 = 0.0

        do j = 1, ny
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

            ! Get velocity components on each cell face 
            ux_1(i,j) = ux(im1,j)
            ux_2(i,j) = ux(i,j)
            uy_1(i,j) = uy(i,jm1)
            uy_2(i,j) = uy(i,j)

            if (ux_1(i,j) >= 0.0) then
                Hx_1(i,j) = H(im1,j)
            else
                Hx_1(i,j) = H(i,j)
            end if

            if (ux_2(i,j) >= 0.0) then
                Hx_2(i,j) = H(i,j)
            else
                Hx_2(i,j) = H(ip1,j)
            end if

            if (uy_1(i,j) >= 0.0) then
                Hy_1(i,j) = H(i,jm1)
            else
                Hy_1(i,j) = H(i,j)
            end if

            if (uy_2(i,j) >= 0.0) then
                Hy_2(i,j) = H(i,j)
            else
                Hy_2(i,j) = H(i,jp1)
            end if

        end do
        end do

        !-------- Assembly of the system of linear equations
        !             (matrix storage: compressed sparse row CSR) --------

        ! Initialize values to zero 
        lgs%a_index = 0.0
        lgs%a_value = 0.0 
        lgs%b_value = 0.0 
        lgs%x_value = 0.0 

        lgs%a_ptr(1) = 1

        k = 0

        do n = 1, lgs%nmax

            ! Get i,j indices of current point
            i = lgs%n2i(n)
            j = lgs%n2j(n)

            nr = n   ! row counter

            ! Get neighbor indices assuming periodic domain
            ! (all other boundary conditions are treated with special cases below)
            im1 = i-1
            if (im1 .eq. 0)    im1 = nx 
            ip1 = i+1
            if (ip1 .eq. nx+1) ip1 = 1 

            jm1 = j-1
            if (jm1 .eq. 0)    jm1 = ny
            jp1 = j+1
            if (jp1 .eq. ny+1) jp1 = 1 
            

            ! Handle special cases first, otherwise populate with normal inner discretization

            if (mask(i,j) .eq. 0) then 
                ! Zero thickness imposed 

                k = k+1
                lgs%a_index(k) = nr
                lgs%a_value(k) = 1.0_wp   ! diagonal element only
                
                lgs%b_value(nr) = 0.0_wp
                lgs%x_value(nr) = 0.0_wp

            else if (mask(i,j) .eq. -1) then 
                ! Prescribed ice thickness imposed 

                k = k+1
                lgs%a_index(k) = nr
                lgs%a_value(k) = 1.0_wp   ! diagonal element only
                
                lgs%b_value(nr) = 0.0_wp
                lgs%x_value(nr) = H(i,j)
            
            else if ( (.not. is_periodic) .and. &
                (i .eq. 1 .or. i .eq. nx .or. j .eq. 1 .or. j .eq. ny) ) then
                ! At outer boundary, and not periodic,
                ! for now set to zero boundary condition
                ! but should be handled in the future via 'boundaries' argument
                ! and/or mask

                k = k+1
                lgs%a_index(k) = nr
                lgs%a_value(k) = 1.0_wp   ! diagonal element only
                
                lgs%b_value(nr) = 0.0_wp
                lgs%x_value(nr) = 0.0_wp

            else
                ! Inner point

                k = k+1
                lgs%a_index(k) = lgs%ij2n(i,jm1)                    ! for H(i,jm1)
                if (uy_1(i,j) > 0.0) &
                    lgs%a_value(k) = -dt_darea*uy_1(i,j)*dx*WOVI

                k = k+1
                lgs%a_index(k) = lgs%ij2n(im1,j)                    ! for H(im1,j)
                if (ux_1(i,j) > 0.0) &
                    lgs%a_value(k) = -dt_darea*ux_1(i,j)*dy*WOVI

                k = k+1
                lgs%a_index(k) = nr                                 ! for H(i,j)
                lgs%a_value(k) = 1.0                                ! (diagonal element)
                if (uy_1(i,j) < 0.0) &
                    lgs%a_value(k) = lgs%a_value(k) &
                                    - dt_darea*uy_1(i,j)*dx*WOVI
                if (ux_1(i,j) < 0.0) &
                    lgs%a_value(k) = lgs%a_value(k) &
                                    - dt_darea*ux_1(i,j)*dy*WOVI
                if (ux_2(i,j) > 0.0) &
                    lgs%a_value(k) = lgs%a_value(k) &
                                    + dt_darea*ux_2(i,j)*dy*WOVI
                if (uy_2(i,j) > 0.0) &
                    lgs%a_value(k) = lgs%a_value(k) &
                                    + dt_darea*uy_2(i,j)*dx*WOVI

                k = k+1
                lgs%a_index(k) = lgs%ij2n(ip1,j)                    ! for H(ip1,j)
                if (ux_2(i,j) < 0.0) &
                    lgs%a_value(k) = dt_darea*ux_2(i,j)*dy*WOVI

                k = k+1
                lgs%a_index(k) = lgs%ij2n(i,jp1)                    ! for H(i,jp1)
                if (uy_2(i,j) < 0.0) &
                    lgs%a_value(k) = dt_darea*uy_2(i,j)*dx*WOVI


                ! Right-hand side 

                lgs%b_value(nr) = H(i,j) + dt*F(i,j) &
                                -(1.0-WOVI) * dt_darea &
                                     * (  ( ux_2(i,j)*Hx_2(i,j)*dy      &
                                           -ux_1(i,j)*Hx_1(i,j)*dy )    &
                                        + ( uy_2(i,j)*Hy_2(i,j)*dx      &
                                           -uy_1(i,j)*Hy_1(i,j)*dx  ) )

                ! Initial guess == previous H

                ! lgs%x_value(nr) = H(i,j)                            
                
            end if

            lgs%x_value(nr) = H(i,j)                            
                
            lgs%a_ptr(nr+1) = k+1   ! row is completed, store index to next row
            
        end do

        ! Done: A, x and b matrices in Ax=b have been populated 
        ! and stored in lgs object. 

        return

    end subroutine linear_solver_matrix_advection_csr_2D

    ! ==== SICOPOLIS5-dev explicit solver =====

    subroutine calc_adv2D_expl_sico(uu,ux,uy,F,dx,dy,dt)
        ! Second-order upwind advection scheme
        ! If input units are uu==[X] and vx/vy==[m/a], returns [X/a]

        implicit none 

        real(wp), intent(OUT)   :: uu(:,:)  ! [X] Variable of interest (aa nodes)
        real(wp), intent(IN)    :: ux(:,:)  ! [m a-1] Horizontal velocity x-direction (ac nodes)
        real(wp), intent(IN)    :: uy(:,:)  ! [m a-1] Horizontal velocity y-direction (ac nodes)
        real(wp), intent(IN)    :: F(:,:)   ! [m a-1] Source term (aa nodes)
        real(wp), intent(IN)    :: dx       ! [m] Horizontal step x-direction
        real(wp), intent(IN)    :: dy       ! [m] Horizontal step y-direction 
        real(wp), intent(IN)    :: dt       ! [a] Time step 

        ! Local variables  
        integer    :: i, j, nx, ny 
        real(wp), allocatable :: uu_adv(:,:)
        
        nx = size(uu,1)
        ny = size(uu,2) 

        allocate(uu_adv(nx,ny))

        uu_adv = 0.0 

        do i = 2, nx-1 
        do j = 2, ny-1 

            call calc_advec_horizontal_point_expl(uu_adv(i,j),uu,ux,uy,dx,dy,i,j)

        end do 
        end do 

        ! Apply advection and source term
        ! [X] = [X] - [time]*[X/time]+[time]*[X/time]
        uu = uu + dt*F - dt*uu_adv
        
        return

    end subroutine calc_adv2D_expl_sico

    subroutine calc_advec_horizontal_point_expl(advecxy,var,ux,uy,dx,dy,i,j)
        ! Newly implemented advection algorithms (ajr)
        ! First-order and second-order upwind implementations   
        ! [m-1] * [m a-1] * [X] = [X a-1]
        ! Output: [X a-1]

        ! Note: Adapted from sicopolis v5-dev 

        implicit none

        real(wp), intent(OUT) :: advecxy
        real(wp), intent(IN)  :: var(:,:)   ! Enth, T, age, H_ice, etc...
        real(wp), intent(IN)  :: ux(:,:) 
        real(wp), intent(IN)  :: uy(:,:)
        real(wp), intent(IN)  :: dx, dy   
        integer,  intent(IN)  :: i, j 

        ! Local variables 
        integer  :: k, nx, ny, nz 
        real(wp) :: ux_1, ux_2, uy_1, uy_2
        real(wp) :: up_x_1, up_x_2, up_y_1, up_y_2 

        nx = size(var,1)
        ny = size(var,2)

        if (i .eq. 1 .or. i .eq. nx .or. j .eq. 1 .or. j .eq. ny) then 
            write(*,*) "Cannot advect at border!"
            stop 
        end if 

        ux_1 = ux(i-1,j)
        ux_2 = ux(i,j)
        uy_1 = uy(i,j-1)
        uy_2 = uy(i,j)

        up_x_1 = 0.0 
        up_x_2 = 0.0 
        up_y_1 = 0.0 
        up_y_2 = 0.0 

        if (ux_1 >= 0.0) then
            up_x_1 = var(i-1,j)
        else
            up_x_1 = var(i,j)
        end if

        if (ux_2 >= 0.0) then
            up_x_2 = var(i,j)
        else
            up_x_2 = var(i+1,j)
        end if

        if (uy_1 >= 0.0) then
            up_y_1 = var(i,j-1)
        else
            up_y_1 = var(i,j)
        end if

        if (uy_2 >= 0.0) then
            up_y_2 = var(i,j)
        else
            up_y_2 = var(i,j+1)
        end if
 
        ! Combine advection terms for total contribution 
        advecxy = (1.0/(dx*dy)) &
                * ( ( ux_2*up_x_2*dx - ux_1*up_x_1*dx ) &
                  + ( uy_2*up_y_2*dy - uy_1*up_y_1*dy ) )
        
        return 

    end subroutine calc_advec_horizontal_point_expl

    ! ==== SICOPOLIS5-dev over-implicit solver =====
    
    subroutine calc_adv2D_impl_sico(uu,ux,uy,F,dx,dy,dt,use_lis)
        !------------------------------------------------------------------------------
        ! Over-implicit solver for the general ice thickness equation.
        !------------------------------------------------------------------------------

        implicit none

        real(wp), intent(INOUT)   :: uu(:,:)  ! [X] Variable of interest (aa nodes)
        real(wp), intent(IN)      :: ux(:,:)  ! [m a-1] Horizontal velocity x-direction (ac nodes)
        real(wp), intent(IN)      :: uy(:,:)  ! [m a-1] Horizontal velocity y-direction (ac nodes)
        real(wp), intent(IN)      :: F(:,:)   ! [m a-1] Net source/sink terms (aa nodes)
        real(wp), intent(IN)      :: dx       ! [m] Horizontal step x-direction
        real(wp), intent(IN)      :: dy       ! [m] Horizontal step y-direction 
        real(wp), intent(IN)      :: dt       ! [a] Time step 
        logical,  intent(IN)      :: use_lis  ! use_lis or sor? 

        ! Local variables  
        integer :: i, j, nx, ny 
        integer :: IMAX, JMAX
        integer :: n, m, k, nnz
        real(wp), allocatable  :: ux_1(:,:), ux_2(:,:)
        real(wp), allocatable  :: uy_1(:,:), uy_2(:,:)
        real(wp), allocatable  :: up_x_1(:,:), up_x_2(:,:)
        real(wp), allocatable  :: up_y_1(:,:), up_y_2(:,:)
        integer,  allocatable  :: ii(:), jj(:), nn(:,:) 
        real(wp) :: dt_darea
        real(wp), parameter    :: OVI_WEIGHT = 1.0  ! Weighing parameter for the over-implicit scheme 
        real(wp), parameter    :: OMEGA_SOR  = 1.0  ! Relaxation parameter for the iterative SOR solver (0 < OMEGA_SOR < 2)
        real(wp), parameter    :: EPS_SOR    = 1e-3 ! [m] Error tolerance
        
        type(linear_solver_class) :: lgs 

! Include header for lis solver fortran interface
#include "lisf.h"
        
        LIS_INTEGER              :: ierr
        LIS_INTEGER              :: iter
        LIS_INTEGER              :: nc, nr
        LIS_INTEGER              :: nmax, n_sprs
        LIS_INTEGER, allocatable :: lgs_a_ptr(:), lgs_a_index(:)
        LIS_INTEGER, allocatable :: lgs_a_diag_index(:)
        LIS_MATRIX               :: lgs_a
        LIS_VECTOR               :: lgs_b, lgs_x
        LIS_SCALAR,  allocatable :: lgs_a_value(:), lgs_b_value(:), lgs_x_value(:)
        LIS_SOLVER               :: solver
        character(len=256)       :: ch_solver_set_option

        nx = size(uu,1)
        ny = size(uu,2)
        
        nmax   =   nx*ny 
        n_sprs = 5*nx*ny 

        allocate(ii(nx*ny),jj(nx*ny))
        allocate(nn(nx,ny))

        allocate(ux_1(nx,ny))
        allocate(ux_2(nx,ny))
        allocate(uy_1(nx,ny))
        allocate(uy_2(nx,ny))
        
        allocate(up_x_1(nx,ny))
        allocate(up_x_2(nx,ny))
        allocate(up_y_1(nx,ny))
        allocate(up_y_2(nx,ny))

        ! ajr:        
        call linear_solver_init(lgs,nx,ny,nvar=1,n_terms=5)

        ! =======================================================================
        !-------- Construction of a vector (with index n) from a 2-d array
        !         (with indices i, j) by diagonal numbering --------
        ! ajr: note ii, jj, and nn are built with 0-indexing in mind, 
        ! later when i,j values are extracted, make sure to add a 1.
        ! ajr: note, this can be done once outside this routine, but for now
        ! do it here.

        ! For keeping consistency with sico (0:nx-1,0:ny-1) indexing        
        IMAX = nx-1
        JMAX = ny-1 

        n=1

        do m=0, IMAX+JMAX
           do i=m, 0, -1
              j = m-i
              if ((i <= IMAX).and.(j <= JMAX)) then
                 ii(n)   = i+1
                 jj(n)   = j+1
                 nn(i+1,j+1) = n
                 n=n+1
              end if
           end do
        end do
        ! =======================================================================


        !-------- Abbreviations --------

        dt_darea = dt/(dx*dy)

        ux_1 = 0.0
        ux_2 = 0.0
        uy_1 = 0.0
        uy_2 = 0.0
        
        up_x_1 = 0.0
        up_x_2 = 0.0
        up_y_1 = 0.0
        up_y_2 = 0.0

        do i = 2, nx-1 
        do j = 2, ny-1

            ux_1(i,j) = ux(i-1,j)
            ux_2(i,j) = ux(i,j)
            uy_1(i,j) = uy(i,j-1)
            uy_2(i,j) = uy(i,j)

            if (ux_1(i,j) >= 0.0) then
                up_x_1(i,j) = uu(i-1,j)
            else
                up_x_1(i,j) = uu(i,j)
            end if

            if (ux_2(i,j) >= 0.0) then
                up_x_2(i,j) = uu(i,j)
            else
                up_x_2(i,j) = uu(i+1,j)
            end if

            if (uy_1(i,j) >= 0.0) then
                up_y_1(i,j) = uu(i,j-1)
            else
                up_y_1(i,j) = uu(i,j)
            end if

            if (uy_2(i,j) >= 0.0) then
                up_y_2(i,j) = uu(i,j)
            else
                up_y_2(i,j) = uu(i,j+1)
            end if

        end do
        end do

        !-------- Assembly of the system of linear equations
        !                     (matrix storage: compressed sparse row CSR) --------

        allocate(lgs_a_value(n_sprs), lgs_a_index(n_sprs), lgs_a_ptr(nmax+1))
        allocate(lgs_a_diag_index(nmax), lgs_b_value(nmax), lgs_x_value(nmax))

        lgs_a_value = 0.0
        lgs_a_index = 0
        lgs_a_ptr   = 0
        lgs_b_value = 0.0
        lgs_x_value = 0.0

        lgs_a_ptr(1) = 1

        k = 0

        do nr=1, nmax   ! loop over rows

            i = ii(nr)
            j = jj(nr)

            if (i .gt. 1 .and. i .lt. nx .and. j .gt. 1 .and. j .lt. ny) then
                ! Inner point 

                k=k+1 ; nc=nn(i,j-1) ; lgs_a_index(k)=nc   ! for uu(i,j-1)
                if (uy_1(i,j) > 0.0) &
                 lgs_a_value(k) = -dt_darea*uy_1(i,j)*dx*OVI_WEIGHT

                k=k+1 ; nc=nn(i-1,j) ; lgs_a_index(k)=nc   ! for uu(i-1,j)
                if (ux_1(i,j) > 0.0) &
                 lgs_a_value(k) = -dt_darea*ux_1(i,j)*dy*OVI_WEIGHT

                k=k+1 ; lgs_a_index(k)=nr ; lgs_a_diag_index(nr)=k  ! for uu(i,j)
                lgs_a_value(k) = 1.0                             ! (diagonal element)
                if (uy_1(i,j) < 0.0) &
                 lgs_a_value(k) = lgs_a_value(k) &
                                  - dt_darea*uy_1(i,j)*dx*OVI_WEIGHT
                if (ux_1(i,j) < 0.0) &
                 lgs_a_value(k) = lgs_a_value(k) &
                                  - dt_darea*ux_1(i,j)*dy*OVI_WEIGHT
                if (ux_2(i,j) > 0.0) &
                 lgs_a_value(k) = lgs_a_value(k) &
                                  + dt_darea*ux_2(i,j)*dy*OVI_WEIGHT
                if (uy_2(i,j) > 0.0) &
                 lgs_a_value(k) = lgs_a_value(k) &
                                  + dt_darea*uy_2(i,j)*dx*OVI_WEIGHT

                k=k+1 ; nc=nn(i+1,j) ; lgs_a_index(k)=nc   ! for uu(i+1,j)
                if (ux_2(i,j) < 0.0) &
                 lgs_a_value(k) = dt_darea*ux_2(i,j)*dy*OVI_WEIGHT

                k=k+1 ; nc=nn(i,j+1) ; lgs_a_index(k)=nc   ! for uu(i,j+1)
                if (uy_2(i,j) < 0.0) &
                 lgs_a_value(k) = dt_darea*uy_2(i,j)*dx*OVI_WEIGHT

                lgs_b_value(nr) = uu(i,j) &
                                +dt*F(i,j) &
                                -(1.0-OVI_WEIGHT) &
                                   * dt_darea &
                                     * (  ( ux_2(i,j)*up_x_2(i,j)*dy   &
                                           -ux_1(i,j)*up_x_1(i,j)*dy ) &
                                        + ( uy_2(i,j)*up_y_2(i,j)*dx    &
                                           -uy_1(i,j)*up_y_1(i,j)*dx  ) )
                ! right-hand side

            else   ! zero-thickness boundary condition

                k = k+1
                lgs_a_value(k)       = 1.0   ! diagonal element only
                lgs_a_diag_index(nr) = k
                lgs_a_index(k)       = nr
                lgs_b_value(nr)      = 0.0

            end if

            lgs_x_value(nr) = uu(i,j)   ! old variable value,
            ! initial guess for solution vector

            lgs_a_ptr(nr+1) = k+1   ! row is completed, store index to next row

        end do

        nnz = k   ! number of non-zero elements of the matrix

        !-------- Solution of the system of linear equations --------

        if (use_lis) then 

            !  ------ Settings for Lis
                    
            call lis_initialize(ierr)           ! Important for parallel computing environments

            call lis_matrix_create(LIS_COMM_WORLD, lgs_a, ierr)
            call lis_vector_create(LIS_COMM_WORLD, lgs_b, ierr)
            call lis_vector_create(LIS_COMM_WORLD, lgs_x, ierr)

            call lis_matrix_set_size(lgs_a, 0, nmax, ierr)
            call lis_vector_set_size(lgs_b, 0, nmax, ierr)
            call lis_vector_set_size(lgs_x, 0, nmax, ierr)

            do nr=1, nmax

                do nc=lgs_a_ptr(nr), lgs_a_ptr(nr+1)-1
                    call lis_matrix_set_value(LIS_INS_VALUE, nr, lgs_a_index(nc), &
                                                        lgs_a_value(nc), lgs_a, ierr)
                end do

                call lis_vector_set_value(LIS_INS_VALUE, nr, lgs_b_value(nr), lgs_b, ierr)
                call lis_vector_set_value(LIS_INS_VALUE, nr, lgs_x_value(nr), lgs_x, ierr)

            end do

            call lis_matrix_set_type(lgs_a, LIS_MATRIX_CSR, ierr)
            call lis_matrix_assemble(lgs_a, ierr)

            !  ------ Solution with Lis

            call lis_solver_create(solver, ierr)

            ch_solver_set_option = '-i bicg -p ilu -maxiter 1000 -tol 1.0e-12 -initx_zeros false'

            call lis_solver_set_option(trim(ch_solver_set_option), solver, ierr)
            call CHKERR(ierr)

            call lis_solve(lgs_a, lgs_b, lgs_x, solver, ierr)
            call CHKERR(ierr)

            call lis_solver_get_iter(solver, iter, ierr)
            !write(6,'(10x,a,i0,5x,i2)') 'calc_adv2D_impl_sico [lis]: iter = ', iter, ierr

            lgs_x_value = 0.0
            call lis_vector_gather(lgs_x, lgs_x_value, ierr)
            call lis_matrix_destroy(lgs_a, ierr)
            call lis_vector_destroy(lgs_b, ierr)
            call lis_vector_destroy(lgs_x, ierr)
            call lis_solver_destroy(solver, ierr)
                   
            call lis_finalize(ierr)     ! Important for parallel computing environments 

        else 
            ! Use internal SOR solver 

            call sor_sprs(lgs_a_value(1:nnz),lgs_a_index(1:nnz),  &
                          lgs_a_diag_index,lgs_a_ptr,lgs_b_value, &
                          OMEGA_SOR, EPS_SOR, lgs_x_value, iter, ierr)

            !write(6,'(10x,a,i0,5x,i2)') 'calc_adv2D_impl_sico [sor]: iter = ', iter, ierr
            
        end if 

        do nr = 1, nmax
            i       = ii(nr)
            j       = jj(nr)
            uu(i,j) = lgs_x_value(nr)
        end do

        deallocate(lgs_a_value, lgs_a_index, lgs_a_ptr)
        deallocate(lgs_a_diag_index, lgs_b_value, lgs_x_value)
        
        return 

    end subroutine calc_adv2D_impl_sico
    
    subroutine calc_adv2D_impl_sico_0(uu,ux,uy,F,dx,dy,dt,use_lis)
        !------------------------------------------------------------------------------
        ! Over-implicit solver for the general ice thickness equation.
        !------------------------------------------------------------------------------

        implicit none

        real(wp), intent(INOUT)   :: uu(:,:)  ! [X] Variable of interest (aa nodes)
        real(wp), intent(IN)      :: ux(:,:)  ! [m a-1] Horizontal velocity x-direction (ac nodes)
        real(wp), intent(IN)      :: uy(:,:)  ! [m a-1] Horizontal velocity y-direction (ac nodes)
        real(wp), intent(IN)      :: F(:,:)   ! [m a-1] Net source/sink terms (aa nodes)
        real(wp), intent(IN)      :: dx       ! [m] Horizontal step x-direction
        real(wp), intent(IN)      :: dy       ! [m] Horizontal step y-direction 
        real(wp), intent(IN)      :: dt       ! [a] Time step 
        logical,  intent(IN)      :: use_lis  ! use_lis or sor? 

        ! Local variables  
        integer :: i, j, nx, ny 
        integer :: IMAX, JMAX
        integer :: n, m, k, nnz
        real(wp), allocatable  :: ux_1(:,:), ux_2(:,:)
        real(wp), allocatable  :: uy_1(:,:), uy_2(:,:)
        real(wp), allocatable  :: up_x_1(:,:), up_x_2(:,:)
        real(wp), allocatable  :: up_y_1(:,:), up_y_2(:,:)
        integer,  allocatable  :: ii(:), jj(:), nn(:,:) 
        real(wp) :: dt_darea
        real(wp), parameter    :: OVI_WEIGHT = 1.0  ! Weighing parameter for the over-implicit scheme 
        real(wp), parameter    :: OMEGA_SOR  = 1.0  ! Relaxation parameter for the iterative SOR solver (0 < OMEGA_SOR < 2)
        real(wp), parameter    :: EPS_SOR    = 1e-3 ! [m] Error tolerance

! Include header for lis solver fortran interface
#include "lisf.h"
        
        LIS_INTEGER              :: ierr
        LIS_INTEGER              :: iter
        LIS_INTEGER              :: nc, nr
        LIS_INTEGER              :: nmax, n_sprs
        LIS_INTEGER, allocatable :: lgs_a_ptr(:), lgs_a_index(:)
        LIS_INTEGER, allocatable :: lgs_a_diag_index(:)
        LIS_MATRIX               :: lgs_a
        LIS_VECTOR               :: lgs_b, lgs_x
        LIS_SCALAR,  allocatable :: lgs_a_value(:), lgs_b_value(:), lgs_x_value(:)
        LIS_SOLVER               :: solver
        character(len=256)       :: ch_solver_set_option

        nx = size(uu,1)
        ny = size(uu,2)
        
        nmax   =   nx*ny 
        n_sprs = 5*nx*ny 

        allocate(ii(nx*ny),jj(nx*ny))
        allocate(nn(nx,ny))

        allocate(ux_1(nx,ny))
        allocate(ux_2(nx,ny))
        allocate(uy_1(nx,ny))
        allocate(uy_2(nx,ny))
        
        allocate(up_x_1(nx,ny))
        allocate(up_x_2(nx,ny))
        allocate(up_y_1(nx,ny))
        allocate(up_y_2(nx,ny))
        
        ! =======================================================================
        !-------- Construction of a vector (with index n) from a 2-d array
        !         (with indices i, j) by diagonal numbering --------
        ! ajr: note ii, jj, and nn are built with 0-indexing in mind, 
        ! later when i,j values are extracted, make sure to add a 1.
        ! ajr: note, this can be done once outside this routine, but for now
        ! do it here.

        ! For keeping consistency with sico (0:nx-1,0:ny-1) indexing        
        IMAX = nx-1
        JMAX = ny-1 

        n=1

        do m=0, IMAX+JMAX
           do i=m, 0, -1
              j = m-i
              if ((i <= IMAX).and.(j <= JMAX)) then
                 ii(n)   = i+1
                 jj(n)   = j+1
                 nn(i+1,j+1) = n
                 n=n+1
              end if
           end do
        end do
        ! =======================================================================


        !-------- Abbreviations --------

        dt_darea = dt/(dx*dy)

        ux_1 = 0.0
        ux_2 = 0.0
        uy_1 = 0.0
        uy_2 = 0.0
        
        up_x_1 = 0.0
        up_x_2 = 0.0
        up_y_1 = 0.0
        up_y_2 = 0.0

        do i = 2, nx-1 
        do j = 2, ny-1

            ux_1(i,j) = ux(i-1,j)
            ux_2(i,j) = ux(i,j)
            uy_1(i,j) = uy(i,j-1)
            uy_2(i,j) = uy(i,j)

            if (ux_1(i,j) >= 0.0) then
                up_x_1(i,j) = uu(i-1,j)
            else
                up_x_1(i,j) = uu(i,j)
            end if

            if (ux_2(i,j) >= 0.0) then
                up_x_2(i,j) = uu(i,j)
            else
                up_x_2(i,j) = uu(i+1,j)
            end if

            if (uy_1(i,j) >= 0.0) then
                up_y_1(i,j) = uu(i,j-1)
            else
                up_y_1(i,j) = uu(i,j)
            end if

            if (uy_2(i,j) >= 0.0) then
                up_y_2(i,j) = uu(i,j)
            else
                up_y_2(i,j) = uu(i,j+1)
            end if

        end do
        end do

        !-------- Assembly of the system of linear equations
        !                     (matrix storage: compressed sparse row CSR) --------

        allocate(lgs_a_value(n_sprs), lgs_a_index(n_sprs), lgs_a_ptr(nmax+1))
        allocate(lgs_a_diag_index(nmax), lgs_b_value(nmax), lgs_x_value(nmax))

        lgs_a_value = 0.0
        lgs_a_index = 0
        lgs_a_ptr   = 0
        lgs_b_value = 0.0
        lgs_x_value = 0.0

        lgs_a_ptr(1) = 1

        k = 0

        do nr=1, nmax   ! loop over rows

            i = ii(nr)
            j = jj(nr)

            if (i .gt. 1 .and. i .lt. nx .and. j .gt. 1 .and. j .lt. ny) then
                ! Inner point 

                k=k+1 ; nc=nn(i,j-1) ; lgs_a_index(k)=nc   ! for uu(i,j-1)
                if (uy_1(i,j) > 0.0) &
                 lgs_a_value(k) = -dt_darea*uy_1(i,j)*dx*OVI_WEIGHT

                k=k+1 ; nc=nn(i-1,j) ; lgs_a_index(k)=nc   ! for uu(i-1,j)
                if (ux_1(i,j) > 0.0) &
                 lgs_a_value(k) = -dt_darea*ux_1(i,j)*dy*OVI_WEIGHT

                k=k+1 ; lgs_a_index(k)=nr ; lgs_a_diag_index(nr)=k  ! for uu(i,j)
                lgs_a_value(k) = 1.0                             ! (diagonal element)
                if (uy_1(i,j) < 0.0) &
                 lgs_a_value(k) = lgs_a_value(k) &
                                  - dt_darea*uy_1(i,j)*dx*OVI_WEIGHT
                if (ux_1(i,j) < 0.0) &
                 lgs_a_value(k) = lgs_a_value(k) &
                                  - dt_darea*ux_1(i,j)*dy*OVI_WEIGHT
                if (ux_2(i,j) > 0.0) &
                 lgs_a_value(k) = lgs_a_value(k) &
                                  + dt_darea*ux_2(i,j)*dy*OVI_WEIGHT
                if (uy_2(i,j) > 0.0) &
                 lgs_a_value(k) = lgs_a_value(k) &
                                  + dt_darea*uy_2(i,j)*dx*OVI_WEIGHT

                k=k+1 ; nc=nn(i+1,j) ; lgs_a_index(k)=nc   ! for uu(i+1,j)
                if (ux_2(i,j) < 0.0) &
                 lgs_a_value(k) = dt_darea*ux_2(i,j)*dy*OVI_WEIGHT

                k=k+1 ; nc=nn(i,j+1) ; lgs_a_index(k)=nc   ! for uu(i,j+1)
                if (uy_2(i,j) < 0.0) &
                 lgs_a_value(k) = dt_darea*uy_2(i,j)*dx*OVI_WEIGHT

                lgs_b_value(nr) = uu(i,j) &
                                +dt*F(i,j) &
                                -(1.0-OVI_WEIGHT) &
                                   * dt_darea &
                                     * (  ( ux_2(i,j)*up_x_2(i,j)*dy   &
                                           -ux_1(i,j)*up_x_1(i,j)*dy ) &
                                        + ( uy_2(i,j)*up_y_2(i,j)*dx    &
                                           -uy_1(i,j)*up_y_1(i,j)*dx  ) )
                ! right-hand side

            else   ! zero-thickness boundary condition

                k = k+1
                lgs_a_value(k)       = 1.0   ! diagonal element only
                lgs_a_diag_index(nr) = k
                lgs_a_index(k)       = nr
                lgs_b_value(nr)      = 0.0

            end if

            lgs_x_value(nr) = uu(i,j)   ! old variable value,
            ! initial guess for solution vector

            lgs_a_ptr(nr+1) = k+1   ! row is completed, store index to next row

        end do

        nnz = k   ! number of non-zero elements of the matrix

        !-------- Solution of the system of linear equations --------

        if (use_lis) then 

            !  ------ Settings for Lis
                    
            call lis_initialize(ierr)           ! Important for parallel computing environments

            call lis_matrix_create(LIS_COMM_WORLD, lgs_a, ierr)
            call lis_vector_create(LIS_COMM_WORLD, lgs_b, ierr)
            call lis_vector_create(LIS_COMM_WORLD, lgs_x, ierr)

            call lis_matrix_set_size(lgs_a, 0, nmax, ierr)
            call lis_vector_set_size(lgs_b, 0, nmax, ierr)
            call lis_vector_set_size(lgs_x, 0, nmax, ierr)

            do nr=1, nmax

                do nc=lgs_a_ptr(nr), lgs_a_ptr(nr+1)-1
                    call lis_matrix_set_value(LIS_INS_VALUE, nr, lgs_a_index(nc), &
                                                        lgs_a_value(nc), lgs_a, ierr)
                end do

                call lis_vector_set_value(LIS_INS_VALUE, nr, lgs_b_value(nr), lgs_b, ierr)
                call lis_vector_set_value(LIS_INS_VALUE, nr, lgs_x_value(nr), lgs_x, ierr)

            end do

            call lis_matrix_set_type(lgs_a, LIS_MATRIX_CSR, ierr)
            call lis_matrix_assemble(lgs_a, ierr)

            !  ------ Solution with Lis

            call lis_solver_create(solver, ierr)

            ch_solver_set_option = '-i bicg -p ilu -maxiter 1000 -tol 1.0e-12 -initx_zeros false'

            call lis_solver_set_option(trim(ch_solver_set_option), solver, ierr)
            call CHKERR(ierr)

            call lis_solve(lgs_a, lgs_b, lgs_x, solver, ierr)
            call CHKERR(ierr)

            call lis_solver_get_iter(solver, iter, ierr)
            !write(6,'(10x,a,i0,5x,i2)') 'calc_adv2D_impl_sico_0 [lis]: iter = ', iter, ierr

            lgs_x_value = 0.0
            call lis_vector_gather(lgs_x, lgs_x_value, ierr)
            call lis_matrix_destroy(lgs_a, ierr)
            call lis_vector_destroy(lgs_b, ierr)
            call lis_vector_destroy(lgs_x, ierr)
            call lis_solver_destroy(solver, ierr)
                   
            call lis_finalize(ierr)     ! Important for parallel computing environments 

        else 
            ! Use internal SOR solver 

            call sor_sprs(lgs_a_value(1:nnz),lgs_a_index(1:nnz),  &
                          lgs_a_diag_index,lgs_a_ptr,lgs_b_value, &
                          OMEGA_SOR, EPS_SOR, lgs_x_value, iter, ierr)

            !write(6,'(10x,a,i0,5x,i2)') 'calc_adv2D_impl_sico [sor]: iter = ', iter, ierr
            
        end if 

        do nr = 1, nmax
            i       = ii(nr)
            j       = jj(nr)
            uu(i,j) = lgs_x_value(nr)
        end do

        deallocate(lgs_a_value, lgs_a_index, lgs_a_ptr)
        deallocate(lgs_a_diag_index, lgs_b_value, lgs_x_value)
        
        return 

    end subroutine calc_adv2D_impl_sico_0
    
    subroutine sor_sprs(lgs_a_value, lgs_a_index, lgs_a_diag_index, lgs_a_ptr, &
                        lgs_b_value, omega, eps_sor, lgs_x_value, iter, ierr)
        !-------------------------------------------------------------------------------
        !> SOR solver for a system of linear equations lgs_a*lgs_x=lgs_b
        !! [matrix storage: compressed sparse row CSR,
        !! represented by arrays lgs_a_value(values), lgs_a_index (indices)
        !! and lgs_a_ptr (pointers)].
        !<------------------------------------------------------------------------------

        implicit none

        real(dp),   intent(in)    :: lgs_a_value(:)         ! size=nnz
        integer,    intent(in)    :: lgs_a_index(:)         ! size=nnz
        integer,    intent(in)    :: lgs_a_diag_index(:)    ! size=nmax
        integer,    intent(in)    :: lgs_a_ptr(:)           ! size=nmax+1
        real(dp),   intent(in)    :: lgs_b_value(:)         ! size=nmax
        real(wp),   intent(in)    :: omega, eps_sor
        real(dp),   intent(inout) :: lgs_x_value(:)         ! size=nmax
        integer,    intent(out)   :: iter
        integer,    intent(out)   :: ierr
        
        ! Local variables 
        integer    :: nnz, nmax
        integer    :: nr, k
        real(wp) :: b_nr
        logical    :: flag_convergence
        real(wp), allocatable :: lgs_x_value_prev(:)
        
        nnz  = size(lgs_a_value)
        nmax = size(lgs_x_value) 

        allocate(lgs_x_value_prev(nmax))

        ! Initially assume convergence criterion is not satisfied 
        ierr = -1   ! convergence criterion not fulfilled
        
        do iter = 1, 1000   ! iter_loop 

            lgs_x_value_prev = lgs_x_value

            do nr = 1, nmax

                b_nr = 0.0_dp 

                do k=lgs_a_ptr(nr), lgs_a_ptr(nr+1)-1
                    b_nr = b_nr + lgs_a_value(k)*lgs_x_value(lgs_a_index(k))
                end do

                lgs_x_value(nr) = lgs_x_value(nr) &
                                -omega*(b_nr-lgs_b_value(nr)) &
                                /lgs_a_value(lgs_a_diag_index(nr))

            end do

            flag_convergence = .true.
            do nr = 1, nmax
                if (abs(lgs_x_value(nr)-lgs_x_value_prev(nr)) > eps_sor) then
                    flag_convergence = .false.
                    exit
                end if
            end do

            if (flag_convergence) then
!                 write(6,'(11x,a,i5)') 'sor_sprs: iter = ', iter
                ierr = 0   ! convergence criterion fulfilled
!                 deallocate(lgs_x_value_prev)
!                 return
                exit 
            end if

        end do ! End iter_loop

        deallocate(lgs_x_value_prev)

        return 

    end subroutine sor_sprs

end module solver_advection_sico

