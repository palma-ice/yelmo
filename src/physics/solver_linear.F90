module solver_linear

! Uncomment next line if petsc should be used
!#define USEPETSC 

    use yelmo_defs, only : sp, dp, wp, io_unit_err

    implicit none 

    type linear_solver_class

        integer :: nmax
        integer :: n_terms
        integer :: n_sprs
        
        integer,  allocatable :: n2i(:)
        integer,  allocatable :: n2j(:)
        integer,  allocatable :: ij2n(:,:)
        
        integer,  allocatable :: a_ptr(:)
        integer,  allocatable :: a_index(:)
        real(wp), allocatable :: a_value(:)
        real(wp), allocatable :: b_value(:)
        real(wp), allocatable :: x_value(:)

        real(wp), allocatable :: resid(:)
        real(wp) :: L1_norm, L2_norm, L2_rel_norm

    end type

    private 

    public :: linear_solver_class
    public :: linear_solver_init
    public :: linear_solver_matrix_solve

contains 
    
    subroutine linear_solver_init(lgs,nx,ny,nvar,n_terms)
        ! Initialize the lgs object that will hold
        ! arrays A, x and b needed to solve linear
        ! equation Ax=b. Arrays are populated in 
        ! another subroutine below - see
        ! linear_solver_matrix_ssa_ac_csr_2D. 

        implicit none

        type(linear_solver_class), intent(INOUT) :: lgs
        integer, intent(IN) :: nx
        integer, intent(IN) :: ny 
        integer, intent(IN) :: nvar 
        integer, intent(IN) :: n_terms 

        ! Local variables
        integer :: i, j, n 

        ! Define array sizes
        lgs%nmax    = nvar*nx*ny 
        lgs%n_terms = n_terms
        lgs%n_sprs  = nvar*n_terms*nx*ny 

        ! Ensure all object arrays are deallocated first
        if (allocated(lgs%n2i))     deallocate(lgs%n2i)
        if (allocated(lgs%n2j))     deallocate(lgs%n2j)
        if (allocated(lgs%ij2n))    deallocate(lgs%ij2n)
        if (allocated(lgs%a_ptr))   deallocate(lgs%a_ptr)
        if (allocated(lgs%a_index)) deallocate(lgs%a_index)
        if (allocated(lgs%a_value)) deallocate(lgs%a_value)
        if (allocated(lgs%b_value)) deallocate(lgs%b_value)
        if (allocated(lgs%x_value)) deallocate(lgs%x_value)
        
        ! Allocate arrays to proper size
        allocate(lgs%n2i(nx*ny))
        allocate(lgs%n2j(nx*ny))
        allocate(lgs%ij2n(nx,ny))
        allocate(lgs%a_value(lgs%n_sprs))
        allocate(lgs%a_index(lgs%n_sprs))
        allocate(lgs%a_ptr(lgs%nmax+1))
        allocate(lgs%b_value(lgs%nmax))
        allocate(lgs%x_value(lgs%nmax))
        
        ! Initialize array values to zero
        lgs%a_value = 0.0
        lgs%a_index = 0
        lgs%a_ptr   = 0

        lgs%b_value = 0.0
        lgs%x_value = 0.0
        
        ! Define indices for reshaping of a 2-d array (with indices i, j)
        ! to a vector (with index n)

        n = 1

        do i = 1, nx
        do j = 1, ny
            lgs%n2i(n)    = i
            lgs%n2j(n)    = j
            lgs%ij2n(i,j) = n
            n = n+1
        end do
        end do

        return

    end subroutine linear_solver_init

! ==== GENERAL INTERFACE TO MATRIX SOLVE ROUTINES =====
    
    subroutine linear_solver_matrix_solve(lgs,lis_settings)
        ! Use LIS to solve matrix equation Ax=b. Take
        ! predefined A, x and b arrays from the lgs object,
        ! and pass them to lis-specific variables, solve
        ! for new x and return to lgs object.

        implicit none 

        type(linear_solver_class), intent(INOUT) :: lgs 
        character(len=*), intent(IN) :: lis_settings        ! LIS solver settings

! ==== PETSC SPECIFIC CODE =====
#ifdef USEPETSC
    
        ! Solve matrix equation with PETSc
        call linear_solver_matrix_solve_petsc(lgs,rtol=1e-6_wp,atol=1e-6_wp,maxits=1000)

#else
        
        ! Solve matrix equation with LIS
        call linear_solver_matrix_solve_lis(lgs,lis_settings)
#endif

        return

    end subroutine linear_solver_matrix_solve

! ==== LIS SPECIFIC CODE =====

    subroutine linear_solver_matrix_solve_lis(lgs,lis_settings)
        ! Use LIS to solve matrix equation Ax=b. Take
        ! predefined A, x and b arrays from the lgs object,
        ! and pass them to lis-specific variables, solve
        ! for new x and return to lgs object.

        implicit none 

        type(linear_solver_class), intent(INOUT) :: lgs 
        character(len=*), intent(IN) :: lis_settings        ! LIS solver settings

        ! =========================================================
        ! LIS-specific variables 

        ! Include header for lis solver fortran interface
#include "lisf.h"
        
        LIS_INTEGER :: ierr
        LIS_INTEGER :: nr
        LIS_INTEGER :: nc
        LIS_INTEGER :: nmax
        LIS_INTEGER :: lin_iter
        LIS_REAL    :: residual 
        LIS_REAL    :: solver_time  
        LIS_MATRIX  :: lgs_a
        LIS_VECTOR  :: lgs_b, lgs_x
        LIS_SOLVER  :: solver

        LIS_INTEGER :: lgs_a_index_now
        LIS_SCALAR  :: lgs_a_value_now
        LIS_SCALAR  :: lgs_b_value_now
        LIS_SCALAR  :: lgs_x_value_now
        LIS_SCALAR, allocatable :: lgs_x_value_out(:)

        ! =========================================================
        
        ! Store nmax in local LIS-specific variable
        nmax = lgs%nmax

        !-------- Settings for Lis --------
               
        call lis_initialize(ierr)           ! Important for parallel computing environments   
        call CHKERR(ierr)

        call lis_matrix_create(LIS_COMM_WORLD, lgs_a, ierr)
        call CHKERR(ierr)
        call lis_vector_create(LIS_COMM_WORLD, lgs_b, ierr)
        call lis_vector_create(LIS_COMM_WORLD, lgs_x, ierr)

        call lis_matrix_set_size(lgs_a, 0, nmax, ierr)
        call CHKERR(ierr)
        call lis_vector_set_size(lgs_b, 0, nmax, ierr)
        call lis_vector_set_size(lgs_x, 0, nmax, ierr)

        ! === Storage order: compressed sparse row (CSR) ===

        do nr=1, nmax

            do nc=lgs%a_ptr(nr), lgs%a_ptr(nr+1)-1

                ! Use temporary values with LIS data types for use with lis routines
                lgs_a_index_now = lgs%a_index(nc)
                lgs_a_value_now = lgs%a_value(nc)
                
                call lis_matrix_set_value(LIS_INS_VALUE, nr, lgs_a_index_now, &
                                                        lgs_a_value_now, lgs_a, ierr)
            end do

            ! Use temporary values with LIS data types for use with lis routines
            lgs_b_value_now = lgs%b_value(nr)
            lgs_x_value_now = lgs%x_value(nr) 

            call lis_vector_set_value(LIS_INS_VALUE, nr, lgs_b_value_now, lgs_b, ierr)
            call lis_vector_set_value(LIS_INS_VALUE, nr, lgs_x_value_now, lgs_x, ierr)

        end do 


        call lis_matrix_set_type(lgs_a, LIS_MATRIX_CSR, ierr)
        call lis_matrix_assemble(lgs_a, ierr)

        !-------- Solution of the system of linear equations with Lis --------

        call lis_solver_create(solver, ierr)

        ! ch_solver_set_option = '-i bicgsafe -p jacobi '// &
        !                         '-maxiter 100 -tol 1.0e-4 -initx_zeros false'
        call lis_solver_set_option(trim(lis_settings), solver, ierr)
        call CHKERR(ierr)

        call lis_solve(lgs_a, lgs_b, lgs_x, solver, ierr)
        call CHKERR(ierr)

        ! Get solver solution information
        call lis_solver_get_iter(solver, lin_iter, ierr)
        call lis_solver_get_time(solver,solver_time,ierr)
        
        ! Obtain the relative L2_norm == ||b-Ax|| / ||b||
        call lis_solver_get_residualnorm(solver,residual,ierr)

        ! Store in lgs object too
        lgs%L2_rel_norm = residual

        ! Print a summary
        !write(*,*) "solve_lis: [time (s), iter, L2_rel_norm] = ", solver_time, lin_iter, residual

        ! Gather x values in local array of lis-type
        allocate(lgs_x_value_out(nmax))
        lgs_x_value_out = 0.0_wp
        call lis_vector_gather(lgs_x, lgs_x_value_out, ierr)
        call CHKERR(ierr)
        
        ! Save to lgs object
        lgs%x_value = lgs_x_value_out
        
        ! Destroy all lis variables
        call lis_matrix_destroy(lgs_a, ierr)
        call CHKERR(ierr)

        call lis_vector_destroy(lgs_b, ierr)
        call lis_vector_destroy(lgs_x, ierr)
        call lis_solver_destroy(solver, ierr)
        call CHKERR(ierr)

        ! Finalize lis.
        call lis_finalize(ierr)           ! Important for parallel computing environments
        call CHKERR(ierr)

        return

    end subroutine linear_solver_matrix_solve_lis


! ==== PETSC SPECIFIC CODE =====
#ifdef USEPETSC

    subroutine linear_solver_matrix_solve_petsc(lgs,rtol,atol,maxits)
        ! Use PETSC Krylov solver to solve matrix equation Ax=b. Take
        ! predefined A, x and b arrays from the lgs object,
        ! and pass them to lis-specific variables, solve
        ! for new x and return to lgs object.

#include <petsc/finclude/petscksp.h>

        use petscksp
        use mpi

        implicit none 

        type(linear_solver_class), intent(INOUT) :: lgs 
        real(wp) :: rtol                    ! Input tolerance parameter
        real(wp) :: atol                    ! Input tolerance parameter
        integer  :: maxits                  ! Input iteration parameter 

        ! Local variables
        
        type(PetscErrorCode) :: perr
        type(PetscInt)       :: lin_iter
        type(tMat)           :: A
        type(tVec)           :: b
        type(tVec)           :: x
        type(tKSP)           :: KSP_solver
        type(PetscReal)      :: PETSc_rtol
        type(PetscReal)      :: PETSc_atol
        type(PetscInt)       :: PETSc_maxits
        
        ! Initialise PETSc MPI stuff
        call initialise_petsc()
        
        ! Convert the CSR matrix from the lgs storage structure to a PETSc matrix
        call convert_lgs_to_petsc( A, b, x, lgs)
        
        ! Set up the KSP solver
        call KSPcreate( PETSC_COMM_WORLD, KSP_solver, perr)
        call handle_error( 'KSPcreate', perr)
      
        ! Set operators. Here the matrix that defines the linear system
        ! also serves as the preconditioning matrix.
        call KSPSetOperators( KSP_solver, A, A, perr)
        call handle_error( 'KSPSetOperators', perr)

        ! Iterative solver tolerances
        PETSc_rtol   = rtol         ! The relative convergence tolerance; relative decrease in the (possibly preconditioned) residual norm
        PETSc_atol   = atol         ! The absolute convergence tolerance; absolute size of the (possibly preconditioned) residual norm 
        PETSc_maxits = maxits       ! Maximum allowed iterations
        call KSPSetTolerances( KSP_solver, PETSc_rtol, PETSc_atol, PETSC_DEFAULT_REAL, PETSc_maxits, perr)
        call handle_error( 'KSPSetTolerances', perr)
        
        ! To start from nonzero guess
        !call KSPSetInitialGuessNonzero( KSP_solver,PetscBool flg);

        ! Set runtime options, e.g.,
        !     -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
        ! These options will override those specified above as long as
        ! KSPSetFromOptions() is called _after_ any other customization routines.
        call KSPSetFromOptions( KSP_solver, perr)
        call handle_error( 'KSPSetFromOptions', perr)
        
        ! Solve the linear system
        call KSPSolve( KSP_solver, b, x, perr)
        call handle_error( 'KSPSolve', perr)
        
        ! Check convergence information
        !call KSPSetConvergenceTest( KSP_solver, PetscErrorCode (*test)(KSP ksp,PetscInt it,PetscReal rnorm, KSPConvergedReason *reason,void *ctx),void *ctx,PetscErrorCode (*destroy)(void *ctx));

        ! Find out how many iterations it took
        call KSPGetIterationNumber( KSP_solver, lin_iter, perr)
        call handle_error( 'KSPGetIterationNumber', perr)
        write(*,*) "  PETSc solved Ax=b in ", lin_iter, " iterations"
        
        ! Get the solution back to the lgs storage structure
        call convert_petsc_solution_to_lgs( lgs, x)
        call sync
        
        ! Clean up after yourself
        call KSPDestroy( KSP_solver, perr)
        call handle_error( 'KSPDestroy', perr)
        call VecDestroy( x, perr)
        call handle_error( 'VecDestroy', perr)
        call VecDestroy( b, perr)
        call handle_error( 'VecDestroy', perr)
        call MatDestroy( A, perr)
        call handle_error( 'MatDestroy', perr)
        
        ! Finalise PETSc MPI stuff
        call finalise_petsc()

        return

        contains
            ! Petsc specific routines 
            
            ! subroutine parse_solver_settings_petsc(rtol,atol,dtol,maxits,settings)
            !     -ksp_rtol <rtol> -ksp_atol <atol> -ksp_divtol <dtol> -ksp_max_it <its>

            !     implicit none

            !     real(wp), intent(OUT) :: rtol 
            !     real(wp), intent(OUT) :: atol 
            !     real(wp), intent(OUT) :: dtol 
            !     integer,  intent(OUT) :: matits 
            !     character(len=*), intent(IN) :: settings 

            !     return

            ! end subroutine parse_solver_settings_petsc

            subroutine convert_lgs_to_petsc(A,b,x,lgs)
                ! Convert lgs storage object to PETSc matrices A, b and x.
                ! Adapted from IMAU-ICE petsc_module.F90, which was
                ! mostly copied from: https://petsc.org/release/documentation/manual/getting_started/#parallel-programming
                  
                implicit none
                
                ! In- and output variables:
                type(tMat), intent(INOUT) :: A
                type(tVec), intent(INOUT) :: b
                type(tVec), intent(INOUT) :: x
                type(linear_solver_class), intent(IN) :: lgs
                
                ! Local variables
                type(PetscErrorCode) :: perr
                integer              :: k, i, j, istart, iend
                type(PetscReal)      :: v
            
                ! == Matrix A ==
                ! ==============
            
                ! Initialise the matrix object
                call MatCreate( PETSC_COMM_WORLD, A, perr)
                call handle_error( 'MatCreate', perr)
                
                ! Set the matrix type to parallel (MPI) Aij
                call MatSetType( A, 'mpiaij', perr)
                call handle_error( 'MatSetType', perr)
                
                ! Set the size, let PETSc automatically determine parallelisation domains
                !call MatSetSizes( A, PETSC_DECIDE, PETSC_DECIDE, CSR%m, CSR%n, perr)
                call MatSetSizes( A, PETSC_DECIDE, PETSC_DECIDE, lgs%nmax, lgs%nmax, perr)
                call handle_error( 'MatSetSizes', perr)
                
                ! Not entirely sure what this one does, but apparently it's really important
                call MatSetFromOptions( A, perr)
                call handle_error( 'MatSetFromOptions', perr)
                
                ! Tell PETSc how much memory needs to be allocated
                !call MatMPIAIJSetPreallocation( A, CSR%nnz_per_row_max+1, PETSC_NULL_INTEGER, CSR%nnz_per_row_max+1, PETSC_NULL_INTEGER, perr)
                call MatMPIAIJSetPreallocation( A, lgs%n_terms+1, PETSC_NULL_INTEGER, lgs%n_terms+1, PETSC_NULL_INTEGER, perr)
                call handle_error( 'MatMPIAIJSetPreallocation', perr)
                
                ! Get parallelisation domains ("ownership ranges")
                call MatGetOwnershipRange( A, istart, iend, perr)
                call handle_error( 'MatGetOwnershipRange', perr)
            
                ! Fill in matrix values
                do i = istart+1, iend ! +1 because PETSc indexes from 0
                    do k = lgs%a_ptr(i), lgs%a_ptr(i+1)-1

                        j = lgs%a_index(k)
                        v = lgs%a_value(k)

                        call MatSetValues( A, 1, i-1, 1, j-1, v, INSERT_VALUES, perr)
                        call handle_error( 'MatSetValues', perr)

                    end do
                end do
                call sync
            
                ! == Vectors b,x ==
                ! =================

                ! Create parallel vectors.
                call VecCreate( PETSC_COMM_WORLD, x, perr)
                call handle_error( 'VecCreate', perr)
                !call VecSetSizes( x, PETSC_DECIDE, CSR%n, perr)
                call VecSetSizes( x, PETSC_DECIDE, lgs%nmax, perr)
                call handle_error( 'VecSetSizes', perr)
                call VecSetFromOptions( x, perr)
                call handle_error( 'VecSetFromOptions', perr)
                
                call VecCreate( PETSC_COMM_WORLD, b, perr)
                call handle_error( 'VecCreate', perr)
                !call VecSetSizes( b, PETSC_DECIDE, CSR%m, perr)
                call VecSetSizes( b, PETSC_DECIDE, lgs%nmax, perr)
                call handle_error( 'VecSetSizes', perr)
                call VecSetFromOptions( b, perr)
                call handle_error( 'VecSetFromOptions', perr)
            
                ! Fill in vector values
                do i = istart+1, iend ! +1 because PETSc indexes from 0

                    v = lgs%b_value( i)
                    call VecSetValues( b, 1, i-1, v, INSERT_VALUES, perr)
                    call handle_error( 'VecSetValues', perr)

                    v = lgs%x_value( i)
                    call VecSetValues( x, 1, i-1, v, INSERT_VALUES, perr)
                    call handle_error( 'VecSetValues', perr)

                end do
                call sync
            
                ! Assemble matrix and vectors, using the 2-step process:
                !   MatAssemblyBegin(), MatAssemblyEnd()
                ! Computations can be done while messages are in transition
                ! by placing code between these two statements.
                
                call MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY, perr)
                call handle_error( 'MatAssemblyBegin', perr)
                call VecAssemblyBegin( b, perr)
                call handle_error( 'VecAssemblyBegin', perr)
                call VecAssemblyBegin( x, perr)
                call handle_error( 'VecAssemblyBegin', perr)
                
                call MatAssemblyEnd(   A, MAT_FINAL_ASSEMBLY, perr)
                call handle_error( 'MatAssemblyEnd', perr)
                call VecAssemblyEnd(   b, perr)
                call handle_error( 'VecAssemblyEnd', perr)
                call VecAssemblyEnd(   x, perr)
                call handle_error( 'VecAssemblyEnd', perr)
                
                return

            end subroutine convert_lgs_to_petsc
        
            subroutine convert_petsc_solution_to_lgs(lgs,x)
      
                implicit none
                
                ! In- and output variables:
                type(linear_solver_class), intent(INOUT) :: lgs
                type(tVec), intent(IN) :: x
                
                ! Local variables
                type(PetscErrorCode) :: perr
                type(PetscInt)       :: istart, iend, i
                type(PetscInt)       :: ix(1)
                type(PetscScalar)    :: v(1)
                
                ! Get parallelisation domains ("ownership ranges")
                call VecGetOwnershipRange( x, istart, iend, perr)
                call handle_error( 'VecGetOwnershipRange', perr)
                
                ! Get values
                do i = istart+1,iend
                    ix(1) = i-1
                    call VecGetValues( x, 1, ix, v, perr)
                    call handle_error( 'VecGetValues', perr)
                    lgs%x_value(i) = v(1)
                end do
                call sync
                
                return

            end subroutine convert_petsc_solution_to_lgs
  
            ! == General PETSc initialisation and finalisation
            subroutine initialise_petsc()
                ! Initialise PETSc

                implicit none

                type(PetscErrorCode) :: perr

                ! Initialise PETSc MPI stuff
                call PetscInitialize( PETSC_NULL_CHARACTER, perr)
                call handle_error( 'PetscInitialize', perr)
                
                return

            end subroutine initialise_petsc

            subroutine finalise_petsc()
                ! Finalise PETSc

                implicit none

                type(PetscErrorCode) :: perr

                ! Finalise PETSc MPI stuff
                call PetscFinalize( perr)
                call handle_error( 'PetscFinalize', perr)
                
                return

            end subroutine finalise_petsc
  
            ! == Error handler
            subroutine handle_error( routine_name, perr)

                implicit none 

                ! In/output variables:
                character(LEN=*),     intent(IN) :: routine_Name
                type(PetscErrorCode), intent(IN) :: perr

                ! Local variables
                integer :: cerr, ierr 

                if (perr /= 0) THEN
                  write(io_unit_err,*) '    PETSc routine "', routine_name, '" returned error flag ', perr
                  call MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
                  stop
                end if

                return 

            end subroutine handle_error

            ! Synchronise the different processes
            subroutine sync
                ! Use MPI_BARRIER to synchronise all the processes         

                implicit none

                integer :: ierr

                call MPI_BARRIER( MPI_COMM_WORLD, ierr)

                return

            end subroutine sync

    end subroutine linear_solver_matrix_solve_petsc

#endif

end module solver_linear
