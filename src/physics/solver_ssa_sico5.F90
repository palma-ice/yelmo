module solver_ssa_sico5
    ! This ssa solver code was adapted from SICOPOLIS (v5-dev, svn revision 1421) module calc_vxy_m.F90. 
    
    use yelmo_defs, only : sp, dp, prec, rho_ice, rho_sw, g

    implicit none 

    private 
    public :: calc_vxy_ssa_matrix 

contains 

    subroutine calc_vxy_ssa_matrix(vx_m,vy_m,vx_bnd,vy_bnd,beta_acx,beta_acy,visc_eff, &
                    ssa_mask_acx,ssa_mask_acy,H_ice,taud_acx,taud_acy,H_grnd,z_sl,z_bed, &
                    dx,dy,ulim,boundaries,gfa1,gfa2,gfb1,gfb2)
        ! Solution of the system of linear equations for the horizontal velocities
        ! vx_m, vy_m in the shallow shelf approximation.
        ! Adapted from sicopolis version 5-dev (svn revision 1421)
        ! Uses the LIS library 
    
        implicit none

        real(prec), intent(INOUT) :: vx_m(:,:)
        real(prec), intent(INOUT) :: vy_m(:,:)
        real(prec), intent(IN)    :: vx_bnd(:,:) 
        real(prec), intent(IN)    :: vy_bnd(:,:) 
        real(prec), intent(IN)    :: beta_acx(:,:)
        real(prec), intent(IN)    :: beta_acy(:,:)
        real(prec), intent(IN)    :: visc_eff(:,:)
        integer,    intent(IN)    :: ssa_mask_acx(:,:) 
        integer,    intent(IN)    :: ssa_mask_acy(:,:) 
        real(prec), intent(IN)    :: H_ice(:,:)
        real(prec), intent(IN)    :: taud_acx(:,:)
        real(prec), intent(IN)    :: taud_acy(:,:)   
        real(prec), intent(IN)    :: H_grnd(:,:)  
        real(prec), intent(IN)    :: z_sl(:,:) 
        real(prec), intent(IN)    :: z_bed(:,:) 
        real(prec), intent(IN)    :: dx, dy
        real(prec), intent(IN)    :: ulim 
        character(len=*), intent(IN) :: boundaries 
        integer,    intent(INOUT) :: gfa1(:,:) 
        logical,    intent(INOUT) :: gfa2(:,:) 
        logical,    intent(INOUT) :: gfb1(:,:) 
        logical,    intent(INOUT) :: gfb2(:,:) 
        
        ! Local variables
        integer    :: nx, ny
        real(prec) :: dxi, deta
        integer    :: i, j, k, n, m 
        integer    :: i1, j1, i00, j00
        real(prec) :: inv_dxi, inv_deta, inv_dxi_deta, inv_dxi2, inv_deta2
        real(prec) :: factor_rhs_2, factor_rhs_3a, factor_rhs_3b
        real(prec) :: rho_sw_ice, H_mid, beta_now, taud_now, H_ocn_now   
        character(len=256) :: ch_solver_set_option
        integer    :: IMAX, JMAX 

        integer, allocatable    :: n2i(:), n2j(:)
        integer, allocatable    :: ij2n(:,:)
        integer, allocatable    :: maske(:,:)
        !logical, allocatable    :: flag_grounded_front_a_1(:,:) 
        !logical, allocatable    :: flag_grounded_front_a_2(:,:) 
        logical, allocatable    :: flag_grounded_front_b_1(:,:) 
        logical, allocatable    :: flag_grounded_front_b_2(:,:) 
        logical, allocatable    :: flag_grounding_line_1(:,:) 
        logical, allocatable    :: flag_grounding_line_2(:,:) 
        logical, allocatable    :: flag_ice_front_1(:,:)
        logical, allocatable    :: flag_ice_front_2(:,:)  
        real(prec), allocatable :: vis_int_g(:,:) 
        real(prec), allocatable :: vis_int_sgxy(:,:) 
        logical :: is_mismip 

! Include header for lis solver fortran interface
#include "lisf.h"
        
        LIS_INTEGER :: ierr
        LIS_INTEGER :: nc, nr
        ! LIS_INTEGER :: iter
        LIS_MATRIX  :: lgs_a
        LIS_VECTOR  :: lgs_b, lgs_x
        LIS_SOLVER  :: solver

        LIS_INTEGER :: nmax, n_sprs 
!         LIS_INTEGER, parameter                 :: nmax   =  2*size(H_ice,1)*size(H_ice,2)
!         LIS_INTEGER, parameter                 :: n_sprs = 20*size(H_ice,1)*size(H_ice,2)
        LIS_INTEGER, allocatable, dimension(:) :: lgs_a_ptr, lgs_a_index
        LIS_SCALAR,  allocatable, dimension(:) :: lgs_a_value, lgs_b_value, lgs_x_value

        is_mismip = .FALSE. 
        if (trim(boundaries) .eq. "MISMIP3D") is_mismip = .TRUE. 

        nx = size(H_ice,1)
        ny = size(H_ice,2)
        
        nmax   =  2*nx*ny 
        n_sprs = 20*nx*ny 

        allocate(n2i(nx*ny),n2j(nx*ny))
        allocate(ij2n(nx,ny))

        allocate(maske(nx,ny))
        allocate(flag_grounding_line_1(nx,ny))
        allocate(flag_grounding_line_2(nx,ny))
        allocate(flag_ice_front_1(nx,ny))
        allocate(flag_ice_front_2(nx,ny))
        
        allocate(vis_int_g(nx,ny))
        allocate(vis_int_sgxy(nx,ny))

        !--- External yelmo arguments => local sicopolis variable names ---
        dxi          = dx 
        deta         = dy 

        vis_int_g    = visc_eff 

        rho_sw_ice   = rho_sw/rho_ice ! Ratio of density of seawater to ice [--]

        !-------- Abbreviations --------

        inv_dxi       = 1.0_prec/dxi
        inv_deta      = 1.0_prec/deta
        inv_dxi_deta  = 1.0_prec/(dxi*deta)
        inv_dxi2      = 1.0_prec/(dxi*dxi)
        inv_deta2     = 1.0_prec/(deta*deta)

        factor_rhs_2  = 0.5_prec*rho_ice*G*(rho_sw-rho_ice)/rho_sw
        factor_rhs_3a = 0.5_prec*rho_ice*G
        factor_rhs_3b = 0.5_prec*rho_sw*G

        ! Set maske and grounding line / calving front flags

!         call set_sico_masks(maske,flag_ice_front_1,flag_ice_front_2, &
!                     flag_grounding_line_1,flag_grounding_line_2, H_ice, H_grnd)
        call set_sico_masks_old(maske,flag_ice_front_1,flag_ice_front_2, &
                    flag_grounding_line_1,flag_grounding_line_2, H_ice, H_grnd)

        gfa1 = maske
        gfb1 = flag_ice_front_1
        gfb2 = flag_ice_front_2
        
        !-------- Depth-integrated viscosity on the staggered grid
        !                                       [at (i+1/2,j+1/2)] --------

        call stagger_visc_aa_ab(vis_int_sgxy,vis_int_g,H_ice)

        !-------- Basal drag parameter (for shelfy stream) --------

        !  (now provided as input matrix)

        ! =======================================================================
        !-------- Reshaping of a 2-d array (with indices i, j)
        !                                  to a vector (with index n) --------
        ! ajr: note, this can be done once outside this routine, but for now
        ! do it here to keep solver portable.

        n=1

        do i=1, nx
        do j=1, ny
           n2i(n)    = i
           n2j(n)    = j
           ij2n(i,j) = n
           n=n+1
        end do
        end do

        ! =======================================================================

        !-------- Assembly of the system of linear equations
        !                         (matrix storage: compressed sparse row CSR) --------

        allocate(lgs_a_value(n_sprs), lgs_a_index(n_sprs), lgs_a_ptr(nmax+1))
        allocate(lgs_b_value(nmax), lgs_x_value(nmax))

        lgs_a_value = 0.0
        lgs_a_index = 0
        lgs_a_ptr   = 0

        lgs_b_value = 0.0
        lgs_x_value = 0.0

        lgs_a_ptr(1) = 1

        k = 0

do n=1, nmax-1, 2

   i = n2i((n+1)/2)
   j = n2j((n+1)/2)

!  ------ Equations for vx_m (at (i+1/2,j))

   nr = n   ! row counter

    if ( (i /= nx).and.(j /= 1).and.(j /= ny) ) then
      ! inner point on the staggered grid in x-direction

      H_mid     = 0.50_prec*(H_ice(i,j)+H_ice(i+1,j))
      beta_now  = beta_acx(i,j) 
      taud_now  = taud_acx(i,j) 

      if (is_mismip .and. i == 1) then
            ! MISMIP3D: free-slip border, vx = neighbor vx

            k  = k+1
            ! if (k > n_sprs) stop ' calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k)  = 1.0   ! diagonal element only
            lgs_a_index(k)  = nr

            lgs_b_value(nr) = 0.0
            lgs_x_value(nr) = 0.0
      
      else if (ssa_mask_acx(i,j) .eq. -1) then 
        ! Assign prescribed boundary velocity to this point
        ! (eg for prescribed velocity corresponding to prescribed grounding line flux)

        ! TO DO 
        write(*,*) "solver_ssa_sico5:: Error: prescribed boundary conditions not yet tested!"
        stop "solver_ssa_sico5 error, see log."
        
        k  = k+1
        ! if (k > n_sprs) stop ' calc_vxy_ssa_matrix: n_sprs too small!'
        lgs_a_value(k)  = 1.0   ! diagonal element only
        lgs_a_index(k)  = nr

        lgs_b_value(nr) = vx_bnd(i,j)
        lgs_x_value(nr) = vx_bnd(i,j)
      
      ! === Proceed with normal ssa checks =================
      
      else if (ssa_mask_acx(i,j) .gt. 0) then 

        if ( &
              ( flag_ice_front_1(i,j).and.flag_ice_front_2(i+1,j) ) &
              .or. &
              ( flag_ice_front_2(i,j).and.flag_ice_front_1(i+1,j) ) &
            ) then
            ! one neighbour is floating ice and the other is ocean
            ! (calving front)

           if (flag_ice_front_1(i,j)) then
              i1 = i     ! floating ice marker
           else   ! flag_ice_front_1(i+1,j)==.true.
              i1 = i+1   ! floating ice marker
           end if

           if (.not.( flag_ice_front_2(i1-1,j) &
                      .and. &
                      flag_ice_front_2(i1+1,j) ) ) then
              ! discretization of the x-component of the BC

              nc = 2*ij2n(i1-1,j)-1
                       ! smallest nc (column counter), for vx_m(i1-1,j)
              k  = k+1
              lgs_a_value(k) = -4.0_prec*inv_dxi*vis_int_g(i1,j)
              lgs_a_index(k) = nc

              nc = 2*ij2n(i1,j-1)
                       ! next nc (column counter), for vy_m(i1,j-1)
              k  = k+1
              lgs_a_value(k) = -2.0_prec*inv_deta*vis_int_g(i1,j)
              lgs_a_index(k) = nc

              nc = 2*ij2n(i1,j)-1
                       ! next nc (column counter), for vx_m(i1,j)
              k  = k+1
              lgs_a_value(k) = 4.0_prec*inv_dxi*vis_int_g(i1,j)
              lgs_a_index(k) = nc

              nc = 2*ij2n(i1,j)
                       ! largest nc (column counter), for vy_m(i1,j)
              k  = k+1
              lgs_a_value(k) = 2.0_prec*inv_deta*vis_int_g(i1,j)
              lgs_a_index(k) = nc

              !lgs_b_value(nr) = factor_rhs_2*H_ice(i1,j)*H_ice(i1,j)

              ! =========================================================
              ! Generalized solution for all ice fronts (floating and grounded)

              if (z_sl(i1,j)-z_bed(i1,j) .gt. 0.0) then 
                ! Bed below sea level 
                  H_ocn_now = min(rho_ice/rho_sw*H_ice(i1,j), &    ! Flotation depth 
                                  z_sl(i1,j)-z_bed(i1,j))          ! Grounded depth 

              else 
                ! Bed above sea level 
                H_ocn_now = 0.0 

              end if 

               lgs_b_value(nr) = factor_rhs_3a*H_ice(i1,j)*H_ice(i1,j) &
                               - factor_rhs_3b*H_ocn_now*H_ocn_now

              !write(*,*) "ssaxc:", H_ice(i1,j), H_ocn_now, H_ocn_now/H_ice(i1,j)

              ! =========================================================
              

              lgs_x_value(nr) = vx_m(i,j)

           else   !      (flag_ice_front_2(i1-1,j)==.true.)
                  ! .and.(flag_ice_front_2(i1+1,j)==.true.);
                  ! velocity assumed to be zero

              k  = k+1
              lgs_a_value(k) = 1.0_prec   ! diagonal element only
              lgs_a_index(k) = nr

              lgs_b_value(nr) = 0.0_prec

              lgs_x_value(nr) = 0.0_prec

           end if 

        else if ( &
                  ( (maske(i,j)==3).and.(maske(i+1,j)==1) ) &
                  .or. &
                  ( (maske(i,j)==1).and.(maske(i+1,j)==3) ) &
                ) then
                ! one neighbour is floating ice and the other is ice-free land;
                ! velocity assumed to be zero

           k  = k+1
           ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
           lgs_a_value(k) = 1.0_prec   ! diagonal element only
           lgs_a_index(k) = nr

           lgs_b_value(nr) = 0.0_prec

           lgs_x_value(nr) = 0.0_prec

        else
            ! inner shelfy stream or floating ice 

            if (i .eq. 1) then  ! ajr: filler to avoid LIS errors 
              k  = k+1
              lgs_a_value(k) = 1.0 
              lgs_a_index(k) = nr 
            else 
              nc = 2*ij2n(i-1,j)-1
                     ! smallest nc (column counter), for vx_m(i-1,j)
              k  = k+1
              ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
              lgs_a_value(k) = 4.0_prec*inv_dxi2*vis_int_g(i,j)
              lgs_a_index(k) = nc
            end if 

            nc = 2*ij2n(i,j-1)-1
                     ! next nc (column counter), for vx_m(i,j-1)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = inv_deta2*vis_int_sgxy(i,j-1)
            lgs_a_index(k) = nc

            nc = 2*ij2n(i,j-1)
                     ! next nc (column counter), for vy_m(i,j-1)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = inv_dxi_deta &
                                    *(2.0_prec*vis_int_g(i,j)+vis_int_sgxy(i,j-1))
            lgs_a_index(k) = nc

            nc = 2*ij2n(i,j)-1
                     ! next nc (column counter), for vx_m(i,j)
!             if (nc /= nr) then   ! (diagonal element)
!                errormsg = ' >>> calc_vxy_ssa_matrix: ' &
!                              //'Check for diagonal element failed!'
!                call error(errormsg)
!             end if
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = -4.0_prec*inv_dxi2 &
                                    *(vis_int_g(i+1,j)+vis_int_g(i,j)) &
                             -inv_deta2 &
                                    *(vis_int_sgxy(i,j)+vis_int_sgxy(i,j-1)) &
                             -beta_now
            lgs_a_index(k) = nc

            nc = 2*ij2n(i,j)
                     ! next nc (column counter), for vy_m(i,j)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = -inv_dxi_deta &
                                    *(2.0_prec*vis_int_g(i,j)+vis_int_sgxy(i,j))
            lgs_a_index(k) = nc

            nc = 2*ij2n(i,j+1)-1
                     ! next nc (column counter), for vx_m(i,j+1)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = inv_deta2*vis_int_sgxy(i,j)
            lgs_a_index(k) = nc

            nc = 2*ij2n(i+1,j-1)
                     ! next nc (column counter), for vy_m(i+1,j-1)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = -inv_dxi_deta &
                                  *(2.0_prec*vis_int_g(i+1,j)+vis_int_sgxy(i,j-1))
            lgs_a_index(k) = nc

            nc = 2*ij2n(i+1,j)-1
                     ! next nc (column counter), for vx_m(i+1,j)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = 4.0_prec*inv_dxi2*vis_int_g(i+1,j)
            lgs_a_index(k) = nc

            nc = 2*ij2n(i+1,j)
                     ! largest nc (column counter), for vy_m(i+1,j)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = inv_dxi_deta &
                                    *(2.0_prec*vis_int_g(i+1,j)+vis_int_sgxy(i,j))
            lgs_a_index(k) = nc

            lgs_b_value(nr) = taud_now

            lgs_x_value(nr) = vx_m(i,j)

        end if

      else   ! neither neighbour is floating or grounded ice,
             ! velocity assumed to be zero

         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = 1.0_prec   ! diagonal element only
         lgs_a_index(k) = nr

         lgs_b_value(nr) = 0.0_prec

         lgs_x_value(nr) = 0.0_prec

      end if

   else
      ! boundary condition, velocity assumed to be zero (unless mismsip)

        if (is_mismip .and. j == 1) then
            ! MISMIP3D: free-slip border, vx = neighbor vx

            k  = k+1
            ! if (k > n_sprs) stop ' calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = 1.0   ! diagonal element only
            lgs_a_index(k) = nr

            lgs_b_value(nr) = vx_m(i,j+1)
            lgs_x_value(nr) = vx_m(i,j+1)
            
        else if (is_mismip .and. j == ny) then
            ! MISMIP3D: free-slip border, vx = neighbor vx
            k  = k+1
            ! if (k > n_sprs) stop ' calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = 1.0   ! diagonal element only
            lgs_a_index(k) = nr

            lgs_b_value(nr) = vx_m(i,j-1)
            lgs_x_value(nr) = vx_m(i,j-1)
        
        else 
            ! velocity assumed to be zero

            k  = k+1
            ! if (k > n_sprs) stop ' calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = 1.0   ! diagonal element only
            lgs_a_index(k) = nr

            lgs_b_value(nr) = 0.0

            lgs_x_value(nr) = 0.0

        end if 

   end if

   lgs_a_ptr(nr+1) = k+1   ! row is completed, store index to next row

!  ------ Equations for vy_m (at (i,j+1/2))

   nr = n+1   ! row counter


   if ( (j /= ny).and.(i /= 1).and.(i /= nx) ) then
      ! inner point on the staggered grid in y-direction

      H_mid     = 0.50_prec*(H_ice(i,j)+H_ice(i,j+1))
      beta_now  = beta_acy(i,j) 
      taud_now  = taud_acy(i,j) 

      if (is_mismip .and. j == 1) then
            ! MISMIP3D: free-slip border, vy=0

            k  = k+1
            ! if (k > n_sprs) stop ' calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = 1.0   ! diagonal element only
            lgs_a_index(k) = nr

            lgs_b_value(nr) = 0.0
            lgs_x_value(nr) = 0.0
            
      else if (ssa_mask_acy(i,j) .eq. -1) then 
        ! Assign prescribed boundary velocity to this point
        ! (eg for prescribed velocity corresponding to prescribed grounding line flux)

        ! TO DO
        write(*,*) "solver_ssa_sico5:: Error: prescribed boundary conditions not yet tested!"
        stop "solver_ssa_sico5 error, see log."

        k  = k+1
        ! if (k > n_sprs) stop ' calc_vxy_ssa_matrix: n_sprs too small!'
        lgs_a_value(k)  = 1.0   ! diagonal element only
        lgs_a_index(k)  = nr

        lgs_b_value(nr) = vy_bnd(i,j)
        lgs_x_value(nr) = vy_bnd(i,j)
      
      ! === Proceed with normal ssa checks =================

      else if ( ssa_mask_acy(i,j) .gt. 0 ) then 
          
        if ( &
              ( flag_ice_front_1(i,j).and.flag_ice_front_2(i,j+1) ) &
              .or. &
              ( flag_ice_front_2(i,j).and.flag_ice_front_1(i,j+1) ) &
            ) then
            ! one neighbour is floating ice and the other is ocean
            ! (calving front)

         if (flag_ice_front_1(i,j)) then
            j1 = j     ! floating ice marker
         else   ! flag_ice_front_1(i,j+1)==.true.
            j1 = j+1   ! floating ice marker
         end if

         if (.not.( flag_ice_front_2(i,j1-1) &
                    .and. &
                    flag_ice_front_2(i,j1+1) ) ) then
            ! discretization of the y-component of the BC

            nc = 2*ij2n(i-1,j1)-1
                     ! smallest nc (column counter), for vx_m(i-1,j1)
            k  = k+1
            lgs_a_value(k) = -2.0_prec*inv_dxi*vis_int_g(i,j1)
            lgs_a_index(k) = nc

            nc = 2*ij2n(i,j1-1)
                     ! next nc (column counter), for vy_m(i,j1-1)
            k  = k+1
            lgs_a_value(k) = -4.0_prec*inv_deta*vis_int_g(i,j1)
            lgs_a_index(k) = nc

            nc = 2*ij2n(i,j1)-1
                     ! next nc (column counter), for vx_m(i,j1)
            k  = k+1
            lgs_a_value(k) = 2.0_prec*inv_dxi*vis_int_g(i,j1)
            lgs_a_index(k) = nc

            nc = 2*ij2n(i,j1)
                     ! largest nc (column counter), for vy_m(i,j1)
            k  = k+1
            lgs_a_value(k) = 4.0_prec*inv_deta*vis_int_g(i,j1)
            lgs_a_index(k) = nc

!             lgs_b_value(nr) = factor_rhs_2*H_ice(i,j1)*H_ice(i,j1)

              ! =========================================================
              ! Generalized solution for all ice fronts (floating and grounded)
              
              if (z_sl(i,j1)-z_bed(i,j1) .gt. 0.0) then 
                ! Bed below sea level 
                  H_ocn_now = min(rho_ice/rho_sw*H_ice(i,j1), &    ! Flotation depth 
                                  z_sl(i,j1)-z_bed(i,j1))          ! Grounded depth 

              else 
                ! Bed above sea level 
                H_ocn_now = 0.0 

              end if 
              
               lgs_b_value(nr) = factor_rhs_3a*H_ice(i,j1)*H_ice(i,j1) &
                               - factor_rhs_3b*H_ocn_now*H_ocn_now


              ! =========================================================
              
            lgs_x_value(nr) = vy_m(i,j)

         else   !      (flag_ice_front_2(i,j1-1)==.true.)
                ! .and.(flag_ice_front_2(i,j1+1)==.true.);
                ! velocity assumed to be zero

            k  = k+1
            lgs_a_value(k) = 1.0_prec   ! diagonal element only
            lgs_a_index(k) = nr

            lgs_b_value(nr) = 0.0_prec

            lgs_x_value(nr) = 0.0_prec

         end if

        else if ( &
                ( (maske(i,j)==3).and.(maske(i,j+1)==1) ) &
                .or. &
                ( (maske(i,j)==1).and.(maske(i,j+1)==3) ) &
              ) then
           ! one neighbour is floating ice and the other is ice-free land;
           ! velocity assumed to be zero

           k  = k+1
           ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
           lgs_a_value(k) = 1.0_prec   ! diagonal element only
           lgs_a_index(k) = nr

           lgs_b_value(nr) = 0.0_prec

           lgs_x_value(nr) = 0.0_prec

        else
            ! inner shelfy stream or floating ice 

            nc = 2*ij2n(i-1,j)-1
                     ! smallest nc (column counter), for vx_m(i-1,j)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = inv_dxi_deta &
                                    *(2.0_prec*vis_int_g(i,j)+vis_int_sgxy(i-1,j))
            lgs_a_index(k) = nc

            nc = 2*ij2n(i-1,j)
                     ! next nc (column counter), for vy_m(i-1,j)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = inv_dxi2*vis_int_sgxy(i-1,j)
            lgs_a_index(k) = nc

            nc = 2*ij2n(i-1,j+1)-1
                     ! next nc (column counter), for vx_m(i-1,j+1)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = -inv_dxi_deta &
                                  *(2.0_prec*vis_int_g(i,j+1)+vis_int_sgxy(i-1,j))
            lgs_a_index(k) = nc

            if (j .eq. 1) then  ! ajr: filler to avoid LIS errors 
              k  = k+1
              lgs_a_value(k) = 1.0   ! diagonal element only
              lgs_a_index(k) = nr
            else
              nc = 2*ij2n(i,j-1)
                       ! next nc (column counter), for vy_m(i,j-1)
              k  = k+1
              ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
              lgs_a_value(k) = 4.0_prec*inv_deta2*vis_int_g(i,j)
              lgs_a_index(k) = nc
            end if 

            nc = 2*ij2n(i,j)-1
                     ! next nc (column counter), for vx_m(i,j)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = -inv_dxi_deta &
                                    *(2.0_prec*vis_int_g(i,j)+vis_int_sgxy(i,j))
            lgs_a_index(k) = nc

            nc = 2*ij2n(i,j)
                     ! next nc (column counter), for vy_m(i,j)
!             if (nc /= nr) then   ! (diagonal element)
!                errormsg = ' >>> calc_vxy_ssa_matrix: ' &
!                              //'Check for diagonal element failed!'
!                call error(errormsg)
!             end if
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = -4.0_prec*inv_deta2 &
                                    *(vis_int_g(i,j+1)+vis_int_g(i,j)) &
                             -inv_dxi2 &
                                    *(vis_int_sgxy(i,j)+vis_int_sgxy(i-1,j)) &
                             -beta_now
            lgs_a_index(k) = nc

            nc = 2*ij2n(i,j+1)-1
                     ! next nc (column counter), for vx_m(i,j+1)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = inv_dxi_deta &
                                    *(2.0_prec*vis_int_g(i,j+1)+vis_int_sgxy(i,j))
            lgs_a_index(k) = nc

            nc = 2*ij2n(i,j+1)
                     ! next nc (column counter), for vy_m(i,j+1)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = 4.0_prec*inv_deta2*vis_int_g(i,j+1)
            lgs_a_index(k) = nc

            nc = 2*ij2n(i+1,j)
                     ! largest nc (column counter), for vy_m(i+1,j)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = inv_dxi2*vis_int_sgxy(i,j)
            lgs_a_index(k) = nc

            lgs_b_value(nr) = taud_now 

            lgs_x_value(nr) = vy_m(i,j)

        end if

      else   ! neither neighbour is floating or grounded ice,
             ! velocity assumed to be zero

         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = 1.0_prec   ! diagonal element only
         lgs_a_index(k) = nr

         lgs_b_value(nr) = 0.0_prec
         lgs_x_value(nr) = 0.0_prec

      end if

   else   ! boundary condition, velocity assumed to be zero

      k  = k+1
      ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
      lgs_a_value(k) = 1.0_prec   ! diagonal element only
      lgs_a_index(k) = nr

      lgs_b_value(nr) = 0.0_prec
      lgs_x_value(nr) = 0.0_prec

   end if

   lgs_a_ptr(nr+1) = k+1   ! row is completed, store index to next row

end do

!-------- Settings for Lis --------

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

!-------- Solution of the system of linear equations with Lis --------

call lis_solver_create(solver, ierr)

ch_solver_set_option = '-i bicgsafe -p jacobi '// &
                        '-maxiter 100 -tol 1.0e-4 -initx_zeros false'

call lis_solver_set_option(trim(ch_solver_set_option), solver, ierr)
call CHKERR(ierr)

call lis_solve(lgs_a, lgs_b, lgs_x, solver, ierr)
call CHKERR(ierr)

!call lis_solver_get_iter(solver, lin_iter, ierr)
!write(6,'(a,i0,a)', advance='no') 'lin_iter = ', lin_iter, ', '

!!! call lis_solver_get_time(solver,solver_time,ierr)
!!! print *, 'calc_vxy_ssa_matrix: time (s) = ', solver_time

lgs_x_value = 0.0_prec
call lis_vector_gather(lgs_x, lgs_x_value, ierr)
call lis_matrix_destroy(lgs_a, ierr)
call lis_vector_destroy(lgs_b, ierr)
call lis_vector_destroy(lgs_x, ierr)
call lis_solver_destroy(solver, ierr)

do n=1, nmax-1, 2

   i = n2i((n+1)/2)
   j = n2j((n+1)/2)

   nr = n
   vx_m(i,j) = lgs_x_value(nr)

   nr = n+1
   vy_m(i,j) = lgs_x_value(nr)

end do

deallocate(lgs_a_value, lgs_a_index, lgs_a_ptr)
deallocate(lgs_b_value, lgs_x_value)

! Limit velocity for thin ice thicknesses ==========
call limit_vel_thickness(vx_m,vy_m,H_ice,u_lim_thk=1000.0)

! Limit the velocity generally =====================
call limit_vel(vx_m,ulim)
call limit_vel(vy_m,ulim)

return 

end subroutine calc_vxy_ssa_matrix

    subroutine stagger_visc_aa_ab(visc_ab,visc,H_ice)

        implicit none 

        real(prec), intent(OUT) :: visc_ab(:,:) 
        real(prec), intent(IN)  :: visc(:,:) 
        real(prec), intent(IN)  :: H_ice(:,:) 

        ! Local variables 
        integer :: i, j, k
        integer :: nx, ny 

        nx = size(visc,1)
        ny = size(visc,2)

        ! Initialisation
        visc_ab = 0.0_prec 

        ! Stagger viscosity only using contributions from neighbors that have ice  
        do i = 1, nx-1 
        do j = 1, ny-1 

            k=0

            if (H_ice(i,j) > 0) then
                k = k+1                              ! floating or grounded ice
                visc_ab(i,j) = visc_ab(i,j) + visc(i,j)
            end if

            if (H_ice(i+1,j) > 0) then
                k = k+1                                  ! floating or grounded ice
                visc_ab(i,j) = visc_ab(i,j) + visc(i+1,j)
            end if

            if (H_ice(i,j+1) > 0) then
                k = k+1                                  ! floating or grounded ice
                visc_ab(i,j) = visc_ab(i,j) + visc(i,j+1)
            end if

            if (H_ice(i+1,j+1) > 0) then
                k = k+1                                      ! floating or grounded ice
                visc_ab(i,j) = visc_ab(i,j) + visc(i+1,j+1)
            end if

            if (k>0) visc_ab(i,j) = visc_ab(i,j)/real(k,prec)

        end do
        end do

        return 

    end subroutine stagger_visc_aa_ab

    subroutine set_sico_masks(maske,front1,front2,gl1,gl2,H_ice,H_grnd)
        ! Define where ssa calculations should be performed
        ! Note: could be binary, but perhaps also distinguish 
        ! grounding line/zone to use this mask for later gl flux corrections
        ! mask = 0: Grounded ice 
        ! mask = 1: Ice-free land 
        ! mask = 2: Open ocean  
        ! mask = 3: Ice shelf 

        ! Note: this mask is defined on central Aa nodes 
        
        implicit none 
        
        integer,    intent(OUT) :: maske(:,:) 
        logical,    intent(OUT) :: front1(:,:)
        logical,    intent(OUT) :: front2(:,:)  
        logical,    intent(OUT) :: gl1(:,:)
        logical,    intent(OUT) :: gl2(:,:) 
        real(prec), intent(IN)  :: H_ice(:,:)
        real(prec), intent(IN)  :: H_grnd(:,:)
        
        ! Local variables
        integer    :: i, j, nx, ny
        integer    :: i1, i2, j1, j2 
        logical    :: is_float 
        
        nx = size(maske,1)
        ny = size(maske,2)
        
        gl1 = .FALSE. 
        gl2 = .FALSE. 

        ! First determine general ice coverage mask 
        ! (land==1,ocean==2,floating_ice==3,grounded_ice==0)

        do j = 1, ny
        do i = 1, nx
            
            ! Check if this point would be floating
            is_float = H_grnd(i,j) .le. 0.0 

            if (H_ice(i,j) .eq. 0.0) then 
                ! Ice-free point 

                if (is_float) then 
                    ! Open ocean 
                    maske(i,j) = 2
                else 
                    ! Ice-free land 
                    maske(i,j) = 1 
                end if 

            else 
                ! Ice-covered point

                if (is_float) then 
                    ! Ice shelf 
                    maske(i,j) = 3
                else
                    ! Grounded ice 
                    maske(i,j) = 0 
                end if 
                
            end if 

        end do 
        end do 
        
        !-------- Detection of the grounding line and the calving front --------

        front1  = .false. 
        front2  = .false. 

        do j = 1, ny
        do i = 1, nx
         
          i1 = max(i-1,1)
          i2 = min(i+1,nx)
          j1 = max(j-1,1)
          j2 = min(j+1,ny)
          
          if (H_ice(i,j) .gt. 0.0 .and. &
              (H_ice(i1,j) .eq. 0.0 .or. &
               H_ice(i2,j) .eq. 0.0 .or. &
               H_ice(i,j1) .eq. 0.0 .or. &
               H_ice(i,j2) .eq. 0.0)) then 

            front1(i,j) = .TRUE. 

          end if 

          if (H_ice(i,j) .eq. 0.0 .and. &
              (H_ice(i1,j) .gt. 0.0 .or. &
               H_ice(i2,j) .gt. 0.0 .or. &
               H_ice(i,j1) .gt. 0.0 .or. &
               H_ice(i,j2) .gt. 0.0)) then 

            front2(i,j) = .TRUE. 

          end if 

        end do
        end do
        
        return
        
    end subroutine set_sico_masks 
    
    subroutine set_sico_masks_old(maske,front1,front2,gl1,gl2,H_ice,H_grnd)
        ! Define where ssa calculations should be performed
        ! Note: could be binary, but perhaps also distinguish 
        ! grounding line/zone to use this mask for later gl flux corrections
        ! mask = 0: Grounded ice 
        ! mask = 1: Ice-free land 
        ! mask = 2: Open ocean  
        ! mask = 3: Ice shelf 

        ! Note: this mask is defined on central Aa nodes 
        
        implicit none 
        
        integer,    intent(OUT) :: maske(:,:) 
        logical,    intent(OUT) :: front1(:,:)
        logical,    intent(OUT) :: front2(:,:)
        logical,    intent(OUT) :: gl1(:,:) 
        logical,    intent(OUT) :: gl2(:,:)  
        real(prec), intent(IN)  :: H_ice(:,:)
        real(prec), intent(IN)  :: H_grnd(:,:)
        
        ! Local variables
        integer    :: i, j, nx, ny
        logical    :: is_float 
        
        nx = size(maske,1)
        ny = size(maske,2)
        
        ! First determine general ice coverage mask 
        ! (land==1,ocean==2,floating_ice==3,grounded_ice==0)

        do j = 1, ny
        do i = 1, nx
            
            ! Check if this point would be floating
            !is_float = H_ice(i,j) < rho_sw_ice*(z_sl(i,j)-z_bed(i,j))
            is_float = H_grnd(i,j) .le. 0.0 
            !is_float = f_grnd(i,j) .eq. 0.0

            if (H_ice(i,j) .eq. 0.0) then 
                ! Ice-free point 

                if (is_float) then 
                    ! Open ocean 
                    maske(i,j) = 2
                else 
                    ! Ice-free land 
                    maske(i,j) = 1 
                end if 

            else 
                ! Ice-covered point

                if (is_float) then 
                    ! Ice shelf 
                    maske(i,j) = 3
                else
                    ! Grounded ice 
                    maske(i,j) = 0 
                end if 
                
            end if 

        end do 
        end do 
        
        !-------- Detection of the grounding line and the calving front --------

        front1  = .false. 
        front2  = .false. 

        ! Grid interior points

        do i=2, nx-1
        do j=2, ny-1

          if ( (maske(i,j)==3) &               ! floating ice
                .and. &
                  (    (maske(i+1,j)==2)   &    ! with
                   .or.(maske(i-1,j)==2)   &    ! one
                   .or.(maske(i,j+1)==2)   &    ! neighbouring
                   .or.(maske(i,j-1)==2) ) &    ! ocean point
              ) &
           front1(i,j) = .true.

           if ( (maske(i,j)==2) &               ! ocean
                .and. &
                  (    (maske(i+1,j)==3)   &    ! with
                   .or.(maske(i-1,j)==3)   &    ! one
                   .or.(maske(i,j+1)==3)   &    ! neighbouring
                   .or.(maske(i,j-1)==3) ) &    ! floating ice point
              ) &
           front2(i,j) = .true.

if (.FALSE.) then 
             if ( (maske(i,j)==0) &              ! grounded ice point
                  .and. &
                    (    (maske(i+1,j)==1)   &   ! with
                     .or.(maske(i-1,j)==1)   &   ! one
                     .or.(maske(i,j+1)==1)   &   ! neighbouring
                     .or.(maske(i,j-1)==1) ) &   ! ice-free land point
                ) &
             front1(i,j) = .true.

             if ( (maske(i,j)==1) &              ! ice-free land point
                  .and. &
                    (    (maske(i+1,j)==0)   &   ! with
                     .or.(maske(i-1,j)==0)   &   ! one
                     .or.(maske(i,j+1)==0)   &   ! neighbouring
                     .or.(maske(i,j-1)==0) ) &   ! grounded ice point
                ) &
             front2(i,j) = .true.
end if 
if (.FALSE.) then 
             if ( (maske(i,j)==0) &              ! grounded ice point
                  .and. &
                    (    (maske(i+1,j)==2)   &   ! with
                     .or.(maske(i-1,j)==2)   &   ! one
                     .or.(maske(i,j+1)==2)   &   ! neighbouring
                     .or.(maske(i,j-1)==2) ) &   ! ocean point
                ) &
             front1(i,j) = .true.

             if ( (maske(i,j)==2) &              ! ocean point
                  .and. &
                    (    (maske(i+1,j)==0)   &   ! with
                     .or.(maske(i-1,j)==0)   &   ! one
                     .or.(maske(i,j+1)==0)   &   ! neighbouring
                     .or.(maske(i,j-1)==0) ) &   ! grounded ice point
                ) &
             front2(i,j) = .true.
end if 

           if ( (maske(i,j)==0) &               ! grounded ice
                .and. &
                  (    (maske(i+1,j)==3)   &    ! with
                   .or.(maske(i-1,j)==3)   &    ! one
                   .or.(maske(i,j+1)==3)   &    ! neighbouring
                   .or.(maske(i,j-1)==3) ) &    ! floating ice point
              ) &
           gl1(i,j) = .true.

           if ( (maske(i,j)==3) &               ! floating ice
                .and. &
                  (    (maske(i+1,j)==0)   &    ! with
                   .or.(maske(i-1,j)==0)   &    ! one
                   .or.(maske(i,j+1)==0)   &    ! neighbouring
                   .or.(maske(i,j-1)==0) ) &    ! grounded ice point
              ) &
           gl2(i,j) = .true.

        end do
        end do

        ! Grid boundary points 

        do i=2, nx-1

           j=1

           if ( (maske(i,j)==2) &               ! ocean
                .and. (maske(i,j+1)==3) &       ! with one neighbouring
              ) &                               ! floating ice point
           front2(i,j) = .true.

           if ( (maske(i,j)==1) &               ! ice-free land point
                .and. (maske(i,j+1)==0) &       ! with one neighbouring
              ) &                               ! grounded ice point
           front2(i,j) = .true.

           if ( (maske(i,j)==2) &               ! ocean point
                .and. (maske(i,j+1)==0) &       ! with one neighbouring
              ) &                               ! grounded ice point
           front2(i,j) = .true.

           j=ny

           if ( (maske(i,j)==2) &               ! ocean
                .and. (maske(i,j-1)==3) &       ! with one neighbouring
              ) &                               ! floating ice point
           front2(i,j) = .true.

           if ( (maske(i,j)==1) &               ! ice-free land point
                .and. (maske(i,j-1)==0) &       ! with one neighbouring
              ) &                               ! grounded ice point
           front2(i,j) = .true.

           if ( (maske(i,j)==2) &               ! ocean point
                .and. (maske(i,j-1)==0) &       ! with one neighbouring
              ) &                               ! grounded ice point
           front2(i,j) = .true.

        end do

        do j=2, ny-1

           i=1

           if ( (maske(i,j)==2) &               ! ocean
                .and. (maske(i+1,j)==3) &       ! with one neighbouring
              ) &                               ! floating ice point
           front2(i,j) = .true.

           if ( (maske(i,j)==1) &               ! ice-free land point
                .and. (maske(i+1,j)==0) &       ! with one neighbouring
              ) &                               ! grounded ice point
           front2(i,j) = .true.

           if ( (maske(i,j)==2) &               ! ocean point
                .and. (maske(i+1,j)==0) &       ! with one neighbouring
              ) &                               ! grounded ice point
           front2(i,j) = .true.

           i=nx

           if ( (maske(i,j)==2) &               ! ocean
                .and. (maske(i-1,j)==3) &       ! with one neighbouring
              ) &                               ! floating ice point
           front2(i,j) = .true.

           if ( (maske(i,j)==1) &               ! ice-free land point
                .and. (maske(i-1,j)==0) &       ! with one neighbouring
              ) &                               ! grounded ice point
           front2(i,j) = .true.

           if ( (maske(i,j)==2) &               ! ocean point
                .and. (maske(i-1,j)==0) &       ! with one neighbouring
              ) &                               ! grounded ice point
           front2(i,j) = .true.

        end do
        
        return
        
    end subroutine set_sico_masks_old 
    
    subroutine limit_vel_thickness(ux,uy,H_ice,u_lim_thk)

        implicit none 

        real(prec), intent(INOUT) :: ux(:,:) 
        real(prec), intent(INOUT) :: uy(:,:) 
        real(prec), intent(IN)    :: H_ice(:,:) 
        real(prec), intent(IN)    :: u_lim_thk 

        ! Local variables 
        integer    :: i, j, nx, ny 
        real(prec) :: H_mid 
        real(prec), parameter :: H_lim = 2.0 

        nx = size(H_ice,1)
        ny = size(H_ice,2) 

        ! acx
        do j = 1, ny 
        do i = 1, nx-1

          H_mid = 0.5*(H_ice(i,j)+H_ice(i+1,j))
          if (H_mid .le. H_lim) call limit_vel(ux(i,j),u_lim_thk)

        end do 
        end do  

        ! acy
        do j = 1, ny-1 
        do i = 1, nx

          H_mid = 0.5*(H_ice(i,j)+H_ice(i,j+1))
          if (H_mid .le. H_lim) call limit_vel(uy(i,j),u_lim_thk)

        end do 
        end do  
        
        return 

    end subroutine limit_vel_thickness

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

end module solver_ssa_sico5
