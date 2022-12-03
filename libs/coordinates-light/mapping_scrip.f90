module coordinates_mapping_scrip
    ! Define and perform mapping using the 
    ! SCRIP file format for storing mapping
    ! weights and neighbors.

    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit

    use ncio 
    use index 
    use interp2D
    use gaussian_filter, only : filter_gaussian, filter_gaussian_fast 
    use grid_to_cdo, only : call_system_cdo

    use, intrinsic :: iso_fortran_env, only: error_unit

    implicit none 


    ! === coord_constants ===================================

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the precision of the coord library (sp,dp)
    integer,  parameter :: wp = sp 


    ! Missing value and aliases
    real(dp), parameter :: MISSING_VALUE_DEFAULT = -9999.0_dp 
    real(dp), parameter :: MV = MISSING_VALUE_DEFAULT
    
    ! Error distance (very large) and error index 
    real(dp), parameter :: ERR_DIST = 1E8_dp 
    integer,  parameter :: ERR_IND  = -1 

    ! Mathematical constants
    real(dp), parameter  :: pi  = 2._dp*acos(0._dp)
    real(dp), parameter  :: degrees_to_radians = pi / 180._dp  ! Conversion factor between radians and degrees
    real(dp), parameter  :: radians_to_degrees = 180._dp / pi  ! Conversion factor between degrees and radians
     
    ! ======================================================

    type map_scrip_class

        character(len=256) :: src_name 
        character(len=256) :: dst_name 
        character(len=512) :: map_fname 

        ! ========================================================
        ! Variables below are defined to be consistent with 
        ! a SCRIP format netcdf file
        ! ========================================================
        
        integer :: src_grid_size
        integer :: dst_grid_size 
        integer :: dst_grid_corners
        integer :: src_grid_rank
        integer :: dst_grid_rank
        integer :: num_links
        integer :: num_wgts

        integer, allocatable :: src_grid_dims(:) 
        integer, allocatable :: dst_grid_dims(:) 
        real(8), allocatable :: src_grid_center_lat(:) 
        real(8), allocatable :: dst_grid_center_lat(:) 
        real(8), allocatable :: src_grid_center_lon(:) 
        real(8), allocatable :: dst_grid_center_lon(:) 
        real(8), allocatable :: dst_grid_corner_lat(:,:)
        real(8), allocatable :: dst_grid_corner_lon(:,:)
        integer, allocatable :: src_grid_imask(:) 
        integer, allocatable :: dst_grid_imask(:) 
        real(8), allocatable :: src_grid_area(:) 
        real(8), allocatable :: dst_grid_area(:)
        real(8), allocatable :: src_grid_frac(:) 
        real(8), allocatable :: dst_grid_frac(:)
        integer, allocatable :: src_address(:)
        integer, allocatable :: dst_address(:)
        real(8), allocatable :: remap_matrix(:,:) 

    end type 

    interface map_scrip_field
        module procedure map_scrip_field_double
        module procedure map_scrip_field_float
        module procedure map_scrip_field_integer 
        module procedure map_scrip_field_logical
    end interface 

    interface nc_read_interp
        module procedure    nc_read_interp_dp_2D
        module procedure    nc_read_interp_dp_3D
        module procedure    nc_read_interp_sp_2D
        module procedure    nc_read_interp_sp_3D
        module procedure    nc_read_interp_int_2D
        module procedure    nc_read_interp_int_3D
        module procedure    nc_read_interp_logical_2D
    end interface

    private 
    public :: map_scrip_class 
    public :: map_scrip_field
    public :: map_scrip_init
    public :: map_scrip_init_from_griddesc
    public :: map_scrip_load 
    public :: map_scrip_end

    public :: nc_read_interp 

    public :: gen_map_filename

contains 
    
    subroutine map_scrip_field_logical(map,var_name,var1,var2,mask2,method,reset,missing_value, &
                                        mask_pack,fill_method,filt_method,filt_par,verbose)
        ! Map a variable field var1 from a src_grid to variable field var2 on dst_grid 

        ! Note: method='mean' is analogous to the method normalize_opt='fracarea' 
        ! desribed in the SCRIP documention (Fig. 2.4 in scripusers.pdf). The 
        ! other methods normalize_opt=['destarea','none'] have not been implemented.

        implicit none 

        type(map_scrip_class),  intent(IN), target :: map 
        character(len=*),       intent(IN)    :: var_name 
        logical,                intent(IN)    :: var1(:,:) 
        logical,                intent(INOUT) :: var2(:,:) 
        integer,                intent(OUT), optional :: mask2(:,:) 
        character(len=*),       intent(IN),  optional :: method
        logical,                intent(IN),  optional :: reset           ! Fill cells with no available values?
        double precision,       intent(IN),  optional :: missing_value   ! Points not included in mapping
        logical,                intent(IN),  optional :: mask_pack(:,:)  ! Mask for where to interpolate
        character(len=*),       intent(IN),  optional :: fill_method     ! Method to fill in remaining missing values
        character(len=*),       intent(IN),  optional :: filt_method     ! Method to use for filtering
        double precision,       intent(IN),  optional :: filt_par(:)     ! gaussian=[sigma,dx]; poisson=[tol]
        logical,                intent(IN),  optional :: verbose         ! Print information
        
        ! Local variables 
        real(dp) :: missing_val
        real(dp) :: filt_pars(2)
        real(dp), allocatable :: var1dp(:,:) 
        real(dp), allocatable :: var2dp(:,:) 
        
        ! By defualt missing value is the coordinates package default value
        missing_val = mv 
        if (present(missing_value)) missing_val = dble(missing_value)

        filt_pars = [1.0_dp,1.0_dp] 
        if (present(filt_par)) filt_pars = dble(filt_par) 

        allocate(var1dp(size(var1,1),size(var1,2)))
        allocate(var2dp(size(var2,1),size(var2,2)))
        
        var1dp = 0.0_dp 
        where(var1) var1dp = 1.0_dp 

        var2dp = 0.0_dp 
        where(var2) var2dp = 1.0_dp 
        
        call map_scrip_field_double(map,var_name,var1dp,var2dp,mask2,method,reset,missing_val, &
                                            mask_pack,fill_method,filt_method,filt_pars)

        var2 = .FALSE. 
        where(int(var2dp) .eq. 1.0_dp) var2 = .TRUE. 

        return 

    end subroutine map_scrip_field_logical
    
    subroutine map_scrip_field_integer(map,var_name,var1,var2,mask2,method,reset,missing_value, &
                                        mask_pack,fill_method,filt_method,filt_par,verbose)
        ! Map a variable field var1 from a src_grid to variable field var2 on dst_grid 

        ! Note: method='mean' is analogous to the method normalize_opt='fracarea' 
        ! desribed in the SCRIP documention (Fig. 2.4 in scripusers.pdf). The 
        ! other methods normalize_opt=['destarea','none'] have not been implemented.

        implicit none 

        type(map_scrip_class),  intent(IN), target :: map 
        character(len=*),       intent(IN)    :: var_name 
        integer,                intent(IN)    :: var1(:,:) 
        integer,                intent(INOUT) :: var2(:,:) 
        integer,                intent(OUT), optional :: mask2(:,:) 
        character(len=*),       intent(IN),  optional :: method
        logical,                intent(IN),  optional :: reset           ! Fill cells with no available values?
        integer,                intent(IN),  optional :: missing_value   ! Points not included in mapping
        logical,                intent(IN),  optional :: mask_pack(:,:)  ! Mask for where to interpolate
        character(len=*),       intent(IN),  optional :: fill_method     ! Method to fill in remaining missing values
        character(len=*),       intent(IN),  optional :: filt_method     ! Method to use for filtering
        integer,                intent(IN),  optional :: filt_par(:)     ! gaussian=[sigma,dx]; poisson=[tol]
        logical,                intent(IN),  optional :: verbose         ! Print information
        
        ! Local variables 
        real(dp) :: missing_val
        real(dp) :: filt_pars(2)
        real(dp), allocatable :: var1dp(:,:) 
        real(dp), allocatable :: var2dp(:,:) 
        
        ! By defualt missing value is the coordinates package default value
        missing_val = mv 
        if (present(missing_value)) missing_val = dble(missing_value)

        filt_pars = [1.0_dp,1.0_dp] 
        if (present(filt_par)) filt_pars = dble(filt_par) 

        allocate(var1dp(size(var1,1),size(var1,2)))
        allocate(var2dp(size(var2,1),size(var2,2)))
        
        var1dp = real(var1,dp)
        var2dp = real(var2,dp)
        
        call map_scrip_field_double(map,var_name,var1dp,var2dp,mask2,method,reset,missing_val, &
                                            mask_pack,fill_method,filt_method,filt_pars)

        var2 = int(var2dp) 

        return 

    end subroutine map_scrip_field_integer

    subroutine map_scrip_field_float(map,var_name,var1,var2,mask2,method,reset,missing_value, &
                                            mask_pack,fill_method,filt_method,filt_par,verbose)
        ! Map a variable field var1 from a src_grid to variable field var2 on dst_grid 

        ! Note: method='mean' is analogous to the method normalize_opt='fracarea' 
        ! desribed in the SCRIP documention (Fig. 2.4 in scripusers.pdf). The 
        ! other methods normalize_opt=['destarea','none'] have not been implemented.

        implicit none 

        type(map_scrip_class),  intent(IN), target :: map 
        character(len=*),       intent(IN)    :: var_name 
        real(sp),               intent(IN)    :: var1(:,:) 
        real(sp),               intent(INOUT) :: var2(:,:) 
        integer,                intent(OUT), optional :: mask2(:,:) 
        character(len=*),       intent(IN),  optional :: method
        logical,                intent(IN),  optional :: reset           ! Reset var2 initially to missing_value?
        real(sp),               intent(IN),  optional :: missing_value   ! Points not included in mapping
        logical,                intent(IN),  optional :: mask_pack(:,:)  ! Mask for where to interpolate
        character(len=*),       intent(IN),  optional :: fill_method     ! Method to fill in remaining missing values
        character(len=*),       intent(IN),  optional :: filt_method     ! Method to use for filtering
        real(sp),               intent(IN),  optional :: filt_par(:)     ! gaussian=[sigma,dx]; poisson=[tol]
        logical,                intent(IN),  optional :: verbose         ! Print information
        
        ! Local variables 
        real(dp) :: missing_val
        real(dp) :: filt_pars(2) 
        real(dp), allocatable :: var1dp(:,:) 
        real(dp), allocatable :: var2dp(:,:) 
        
        ! By defualt missing value is the coordinates package default value
        missing_val = mv 
        if (present(missing_value)) missing_val = dble(missing_value)

        filt_pars = [1.0_dp,1.0_dp] 
        if (present(filt_par)) filt_pars = dble(filt_par) 

        allocate(var1dp(size(var1,1),size(var1,2)))
        allocate(var2dp(size(var2,1),size(var2,2)))
        
        var1dp = real(var1,dp)
        var2dp = real(var2,dp)
        
        call map_scrip_field_double(map,var_name,var1dp,var2dp,mask2,method,reset,missing_val, &
                                            mask_pack,fill_method,filt_method,filt_pars)

        var2 = real(var2dp,sp) 

        return 

    end subroutine map_scrip_field_float

    subroutine map_scrip_field_double(map,var_name,var1,var2,mask2,method,reset,missing_value, &
                                                        mask_pack,fill_method,filt_method,filt_par,verbose)
        ! Map a variable field var1 from a src_grid to variable field var2 on dst_grid 

        ! Note: method='mean' is analogous to the method normalize_opt='fracarea' 
        ! desribed in the SCRIP documention (Fig. 2.4 in scripusers.pdf). The 
        ! other methods normalize_opt=['destarea','none'] have not been implemented.

        implicit none 

        type(map_scrip_class),  intent(IN), target :: map 
        character(len=*),       intent(IN)    :: var_name 
        real(8),                intent(IN)    :: var1(:,:) 
        real(8),                intent(INOUT) :: var2(:,:) 
        integer,                intent(OUT), optional :: mask2(:,:) 
        character(len=*),       intent(IN),  optional :: method
        logical,                intent(IN),  optional :: reset           ! Reset var2 initially to missing_value?
        double precision,       intent(IN),  optional :: missing_value   ! Points not included in mapping
        logical,                intent(IN),  optional :: mask_pack(:,:)  ! Mask for where to interpolate
        character(len=*),       intent(IN),  optional :: fill_method     ! Method to fill in remaining missing values
        character(len=*),       intent(IN),  optional :: filt_method     ! Method to use for filtering
        double precision,       intent(IN),  optional :: filt_par(:)     ! gaussian=[sigma,dx]; poisson=[tol]
        logical,                intent(IN),  optional :: verbose         ! Print information
        
        ! Local variables 
        integer :: n, k, npts1, npts2
        character(len=56) :: method_interp         
        logical           :: reset_pts
        double precision  :: missing_val 
        logical           :: verbose_out
        logical, allocatable  :: maskp(:)
        real(dp), allocatable :: area(:)
        integer :: i, j, j1, j2  

        real(dp), allocatable :: var1_vec(:)
        real(dp), allocatable :: var2_vec(:) 
        integer, allocatable  :: mask2_vec(:) 
        real(dp) :: area_tot, pt_ave, pt_var   
        integer  :: npt_now, num_links_now 
        integer, parameter :: max_num_links_now = 10000
        real(dp), dimension(max_num_links_now) :: var1_now 
        real(dp), dimension(max_num_links_now) :: wts1_now 
        real(dp) :: wts1_tot 

        logical, allocatable  :: maskp2d(:,:) 

        integer :: npts_apply 
        real(dp) :: mean2, mean2b           ! Check mean before/after filtering

        npts1 = size(var1,1)*size(var1,2)
        npts2 = size(var2,1)*size(var2,2)

        ! Confirm that source (var1) and destination (var2)
        ! arrays match map. 
        if (npts1 .ne. map%src_grid_size) then 
            write(*,*) "map_scrip_field:: Error: source array and map size do not match."
            write(*,*) "size(var1): ", size(var1,1), size(var1,2), " = ", npts1
            write(*,*) "map%src_grid_size = ", map%src_grid_size 
            stop 
        end if 
        if (npts2 .ne. map%dst_grid_size) then 
            write(*,*) "map_scrip_field:: Error: dest. array and map size do not match."
            write(*,*) "size(var2): ", size(var2,1), size(var2,2), " = ", npts2
            write(*,*) "map%dst_grid_size = ", map%dst_grid_size 
            stop 
        end if 
            
        ! By default, method is 'mean' (conservative)
        method_interp = "mean"
        if (present(method)) method_interp = trim(method) 

        ! By default, reset target grid points to missing values initially
        reset_pts = .TRUE. 
        if (present(reset)) reset_pts = reset 

        ! By defualt missing value is the coordinates package default value
        missing_val = mv 
        if (present(missing_value)) missing_val = missing_value

        ! By default, not verbose output 
        verbose_out = .FALSE. 
        if (present(verbose)) verbose_out = verbose 

        ! By default, all var2 points are interpolated
        allocate(maskp(npts2))
        maskp = .TRUE. 
        if (present(mask_pack)) maskp = reshape(mask_pack,[npts2])

        ! Count total points to be applied 
        npts_apply = count(maskp)
        
        ! Store var1 in vector format
        allocate(var1_vec(npts1)) 
        var1_vec = reshape(var1,[npts1])

        ! Store var2 in vector format
        allocate(var2_vec(npts2))
        var2_vec = reshape(var2,[npts2])

        ! Allocate mask to keep track of which points have been interpolated 
        allocate(mask2_vec(npts2))
        mask2_vec = 0 

        ! Reset output points to missing values
        if (reset_pts) var2_vec = missing_val 


        j1 = 0 
        j2 = 0 

        ! Loop over target points
        do k = 1, npts2 

            ! Find the range of link indices that correspond 
            ! to the current point k, ie, such that:
            ! map%dst_address(j1:j2) == k 
            ! Note: dst_address can be expected to be sorted 
            ! in ascending order.
            j1 = j2+1 

            ! Check index associated with this address. If 
            ! it is greater than the current index k, it 
            ! means this point has no interpolation links,
            ! so skip this iteration of the main loop.
            if (map%dst_address(j1) .gt. k) then 
                j1 = j1-1 
                cycle 
            end if 

            ! Given j1 is the start of the addresses associated 
            ! with the current index k, find the upper range 
            ! such that map%dst_address(j1:j2) == k and it 
            ! covers all addresses equal to k.
            do j = j1, map%num_links
                if (map%dst_address(j) .eq. map%dst_address(j1) ) then 
                    j2 = j 
                else 
                    exit 
                end if 
            end do 

            ! Determine the number of links 
            num_links_now = j2-j1+1

            if (num_links_now>max_num_links_now) then
              write(*,*) "map_scrip_field:: Error: num_links_now>max_num_links_now: ", &
                                                            num_links_now, max_num_links_now
              write(*,*) " To avoid this error, increase hard-coded variable 'max_num_links_now' &
                         &in coordinates_mapping_scrip.f90."
              write(*,*) 
              stop 
            endif

            if (maskp(k)) then 
                ! Only interpolate for desired target points 
            
                ! Assign data and weights to pointers
                var1_now(1:num_links_now) = var1_vec(map%src_address(j1:j2))
                wts1_now(1:num_links_now) = map%remap_matrix(1,j1:j2)

                ! Calculate the total weight associated with this point,
                ! accounting for missing values in the source array.
                wts1_tot = sum(wts1_now(1:num_links_now),mask=var1_now(1:num_links_now) .ne. missing_val)

                if (wts1_tot .gt. 0.0d0) then 
                    ! Interpolation data found, proceed to interpolate this point
                    
                    var2_vec(k)  = 0.0d0 
                    mask2_vec(k) = 1 

                    select case(trim(method_interp))

                        case("mean")
                            ! Calculate the area-weighted mean 

                            var2_vec(k) = sum((wts1_now(1:num_links_now)/wts1_tot)*var1_now(1:num_links_now), &
                                                                    mask=var1_now(1:num_links_now) .ne. missing_val)

                        case("count")
                            ! Choose the most frequently occurring value, weighted by area

                            var2_vec(k) = maxcount(var1_now(1:num_links_now),wts1_now(1:num_links_now),missing_val)

                        case("stdev")
                            ! Calculate the weighted standard deviation 
                            ! using unbiased estimator correction 

                            npt_now = count(var1_now(1:num_links_now) .ne. missing_val)

                            if (npt_now .gt. 2) then
                                ! Only calculate stdev for 2 or more input points

                                pt_ave      = sum((wts1_now(1:num_links_now)/wts1_tot)*var1_now(1:num_links_now), &
                                                                    mask=var1_now(1:num_links_now) .ne. missing_val)
                                var2_vec(k) = (npt_now/(npt_now - 1.0)) &
                                               * sum((wts1_now(1:num_links_now)/wts1_tot)*(var1_now(1:num_links_now)-pt_ave)**2, & 
                                                                            mask=var1_now(1:num_links_now) .ne. missing_val)
                                var2_vec(k) = sqrt(var2_vec(k))
                                
                            else
                                ! Otherwise assume standard deviation is zero 
                                var2_vec(k) = 0.0d0 

                            end if 

                        case DEFAULT 

                            write(*,*) "map_scrip_field:: Error: interpolation method not recognized."
                            write(*,*) "method = ", trim(method_interp) 
                            stop 

                    end select 

                end if 

            end if 

        end do 

        ! Send back to 2D array 
        var2 = reshape(var2_vec,[size(var2,1),size(var2,2)])

        ! Get interpolation mask too if desired 
        if (present(mask2)) then 
            mask2 = reshape(mask2_vec,[size(var2,1),size(var2,2)])
        end if 

        ! Allocate mask if needed 
        if (present(fill_method) .or. present(filt_method)) then

            allocate(maskp2d(size(var2,1),size(var2,2)))
            maskp2d = reshape(maskp,[size(var2,1),size(var2,2)])

        end if 

        ! === Filling ===
        ! Fill in remaining missing values 

        if (present(fill_method)) then 

            select case(trim(fill_method))

                case("weighted")

                    call fill_weighted(var2,missing_val,n=6,mask=maskp2d)

                case("nn")
                    
                    call fill_nearest(var2,missing_val,mask=maskp2d)

                case("none") ! eg "none"

                    ! Pass - no filling applied 

                case DEFAULT 

                    write(*,*) "map_scrip_field:: Error: fill method not recognized: "//trim(fill_method)
                    write(*,*) "  fill_method = [weighted,nn,none]."
                    write(*,*) 
                    stop 
                    
            end select

        end if 

        ! === Filtering ===
        ! Now perform filtering (smoothing) steps if
        ! the right arguments have been provided. 
        
        if (present(filt_method)) then 

            ! Update maskp2d to also reflect missing values 
            maskp2d = reshape(maskp,[size(var2,1),size(var2,2)])
            maskp2d = (maskp2d .and. var2 .ne. missing_val)

            ! Calculate grid average before filtering 
            if (verbose_out .and. npts_apply .gt. 0) then 
                mean2 = sum(var2,mask=maskp2d) / real(npts_apply,dp)
            end if 

            select case(trim(filt_method))

                case("gaussian")

                    call filter_gaussian(var2,sigma=filt_par(1),dx=filt_par(2),mask=maskp2d)
                    
                case("gaussian-fast")

                    call filter_gaussian_fast(var2,sigma=filt_par(1),dx=filt_par(2),mask=maskp2d)
        

                case("poisson")

                    call filter_poisson(var2,mask=maskp2d,tol=filt_par(1), &
                                    missing_value=missing_val,wrapx=.FALSE.,verbose=.FALSE.)

                case("none")

                    ! Pass - no filtering applied 

                case DEFAULT

                    write(*,*) "map_scrip_field:: Error: filtering method not recognized: "//trim(filt_method)
                    write(*,*) "  filt_method = [gaussian,gaussian-fast,poisson,none]."
                    write(*,*) 
                    stop 
                    
            end select 

            ! Calculate grid average after filtering 
            if (verbose_out .and. npts_apply .gt. 0) then 
                mean2b = sum(var2,mask=maskp2d) / real(npts_apply,dp)
            end if 

            if (verbose_out) then 
                ! Print summary of filtering 
                write(*,"(4a,2g14.5)") var_name, " - ",filt_method, ": mean[orig,filtered]: ", mean2, mean2b
            end if 

        end if 


        return 

    end subroutine map_scrip_field_double

    subroutine map_scrip_init(mps,grid_name_src,grid_name_tgt,method,fldr,load,clean)
        ! Generate mapping weights from grid1 to grid2

        implicit none 

        type(map_scrip_class), intent(INOUT) :: mps
        character(len=*), intent(IN) :: grid_name_src
        character(len=*), intent(IN) :: grid_name_tgt
        character(len=*), intent(IN), optional :: method
        character(len=*), intent(IN), optional :: fldr      ! Directory in which to save/load map
        logical,          intent(IN), optional :: load      ! Whether loading is desired if map exists already
        logical,          intent(IN), optional :: clean     ! Whether to delete intermediate grid desc / grid files

        ! Local variables 
        logical :: load_file, fldr_exists, file_exists 
        character(len=256) :: mapmethod
        character(len=256) :: mapfldr 
        character(len=512) :: filename

        ! Load file if it exists by default
        load_file = .TRUE. 
        if (present(load)) load_file = load 

        mapfldr = "maps"
        if (present(fldr)) mapfldr = trim(fldr)

        mapmethod = "con"
        if (present(method)) mapmethod = trim(method)

        ! Note: do not assign max distance here, save all distances
        ! up until the maximum number of neighbors
        ! Later, when loading map, let use choose max_distance
        ! In this way, less recalculation of maps will be needed
        ! when the max_distance changes.

        ! Determine if file matching these characteristics exists
        filename = gen_map_filename(grid_name_src,grid_name_tgt,mapfldr,mapmethod)
        inquire(file=filename,exist=file_exists)

        !! Now load map information from file if exists and is desired
        !! or else calculate weights and store in file. 
        if ( load_file .and. file_exists ) then 

            ! Read map from file
            call map_scrip_load(mps,grid_name_src,grid_name_tgt,mapfldr,mapmethod)

        else

            ! == Generate the SCRIP map via a cdo call:

            ! To do: see coordinates_mapping_scrip.f90 in coordinates package.
            ! The approach depends on the coordinates object, which has not been
            ! included in this standalone module. For now, it is only possible
            ! to load a scrip map that has already been generated externally. 

            write(error_unit,*) ""
            write(error_unit,*) "map_scrip_init:: Error: scrip map file not found. &
                        &This file should be pregenerated using cdo before running Yelmo. &
                        &See maps/readme.md in the main Yelmo directory for details."
            write(error_unit,*) "scrip map filename: ", trim(filename)
        end if 

        return 

    end subroutine map_scrip_init

    subroutine map_scrip_init_from_gridnc(map,src_name,dst_name,fldr,method,load)
        ! Use cdo to generate scrip map based on grid 
        ! definitions. 

        ! 1. Assume that grid description text files already exist
        !    for each grid. 

        implicit none 

        type(map_scrip_class), intent(INOUT) :: map     ! map object to be initialized
        character(len=*), intent(IN) :: src_name        ! Source grid name
        character(len=*), intent(IN) :: dst_name        ! Dest./target grid name
        character(len=*), intent(IN) :: fldr            ! Folder where grid desciptions can be found
        character(len=*), intent(IN) :: method          ! Interpolation method (con,bil,bic)
        logical,          intent(IN), optional :: load  ! Load map from file if available? 

        ! Local variables 
        character(len=512)  :: fnm1
        character(len=512)  :: fnm2
        character(len=512)  :: fnm_map 
        character(len=2048) :: cdo_cmd
        logical :: load_map 
        logical :: map_exists  
        logical :: cdo_success 

        ! Determine whether map file should be loaded if available 
        load_map = .TRUE. 
        if (present(load)) load_map = load 

        ! Step 1: call cdo to generate mapping weights in a scrip file 

        ! Generate grid description filenames 
        fnm1 = trim(fldr)//"/"//"grid_"//trim(src_name)//".nc"
        fnm2 = trim(fldr)//"/"//"grid_"//trim(dst_name)//".nc"

        ! Determine map filename from grid names and folder 
        fnm_map = gen_map_filename(src_name,dst_name,fldr,method)
        
        ! Check if scrip weights file already exists  
        inquire(file=trim(fnm_map),exist=map_exists)

        if ( (.not. map_exists) .or. (.not. load_map) ) then 
            ! If no map exists yet, or loading is not desired, 
            ! then call cdo to generate a new map file. 

            ! Define cdo command to generate mapping weights from 
            ! src grid (fnm1) to dest grid (fnm2) using example netcdf 
            ! grid file (src_nc) and storing weights in map file (fnm_map).
            cdo_cmd = "cdo gen"//trim(method)//","//trim(fnm2)// &
                            " "//trim(fnm1)//" "//trim(fnm_map)

            ! Call cdo command via system call
            call call_system_cdo(cdo_cmd)

        end if 

        ! Step 2: load map weights and initialize map_scrip_class object 
        call map_scrip_load(map,src_name,dst_name,fldr,method)

        return 

    end subroutine map_scrip_init_from_gridnc

    subroutine map_scrip_init_from_griddesc(map,src_name,dst_name,fldr,src_nc,method,load)
        ! Use cdo to generate scrip map based on grid 
        ! definitions. 

        ! 1. Assume that grid description text files already exist
        !    for each grid. 

        implicit none 

        type(map_scrip_class), intent(INOUT) :: map     ! map object to be initialized
        character(len=*), intent(IN) :: src_name        ! Source grid name
        character(len=*), intent(IN) :: dst_name        ! Dest./target grid name
        character(len=*), intent(IN) :: fldr            ! Folder where grid desciptions can be found
        character(len=*), intent(IN) :: src_nc          ! Path to source netcdf file containing grid/variables (needed by cdo) 
        character(len=*), intent(IN) :: method          ! Interpolation method (con,bil,bic)
        logical,          intent(IN), optional :: load  ! Load map from file if available? 

        ! Local variables 
        character(len=512)  :: fnm1
        character(len=512)  :: fnm2
        character(len=512)  :: fnm_map 
        character(len=2048) :: cdo_cmd
        logical :: load_map 
        logical :: map_exists  
        logical :: cdo_success 

        ! Determine whether map file should be loaded if available 
        load_map = .TRUE. 
        if (present(load)) load_map = load 

        ! Step 1: call cdo to generate mapping weights in a scrip file 

        ! Generate grid description filenames 
        fnm1 = trim(fldr)//"/"//"grid_"//trim(src_name)//".txt"
        fnm2 = trim(fldr)//"/"//"grid_"//trim(dst_name)//".txt"

        ! Determine map filename from grid names and folder 
        fnm_map = gen_map_filename(src_name,dst_name,fldr,method)
        
        ! Check if scrip weights file already exists  
        inquire(file=trim(fnm_map),exist=map_exists)

        if ( (.not. map_exists) .or. (.not. load_map) ) then 
            ! If no map exists yet, or loading is not desired, 
            ! then call cdo to generate a new map file. 

            ! Define cdo command to generate mapping weights from 
            ! src grid (fnm1) to dest grid (fnm2) using example netcdf 
            ! grid file (src_nc) and storing weights in map file (fnm_map).
            cdo_cmd = "cdo gen"//trim(method)//","//trim(fnm2)//" -setgrid,"//trim(fnm1)// &
                    " "//trim(src_nc)//" "//trim(fnm_map)

            ! Call cdo command via system call
            call call_system_cdo(cdo_cmd)

        end if 

        ! Step 2: load map weights and initialize map_scrip_class object 
        call map_scrip_load(map,src_name,dst_name,fldr,method)

        return 

    end subroutine map_scrip_init_from_griddesc

    subroutine map_scrip_load(map,src_name,dst_name,fldr,method)
        ! Load a map_scrip_class object into memory
        ! from a netcdf file. 

        implicit none 

        type(map_scrip_class), intent(INOUT) :: map
        character(len=*), intent(IN) :: src_name
        character(len=*), intent(IN) :: dst_name 
        character(len=*), intent(IN) :: fldr  
        character(len=*), intent(IN) :: method
        
        ! Local variables 
        integer, allocatable :: dims(:) 
        character(len=56), allocatable :: dim_names(:) 

        ! Define map names 
        map%src_name = trim(src_name) 
        map%dst_name = trim(dst_name) 

        ! Determine filename from grid names and folder 
        map%map_fname = gen_map_filename(src_name,dst_name,fldr,method)
        
        ! write(*,*) "Loading SCRIP map from file: "//trim(map%map_fname) 
        ! write(*,*) "" 

        call nc_dims(map%map_fname,"src_grid_center_lat",dim_names,dims)
        map%src_grid_size = dims(1) 
        call nc_dims(map%map_fname,"dst_grid_center_lat",dim_names,dims)
        map%dst_grid_size = dims(1) 
        if (nc_exists_var(map%map_fname,"dst_grid_corner_lat")) then 
          call nc_dims(map%map_fname,"dst_grid_corner_lat",dim_names,dims)
          map%dst_grid_corners = dims(1) 
        else
          map%dst_grid_corners = 0
        endif
        call nc_dims(map%map_fname,"src_grid_dims",dim_names,dims)
        map%src_grid_rank = dims(1) 
        call nc_dims(map%map_fname,"dst_grid_dims",dim_names,dims)
        map%dst_grid_rank = dims(1) 
        call nc_dims(map%map_fname,"remap_matrix",dim_names,dims)
        map%num_wgts  = dims(1)
        map%num_links = dims(2) 

        ! write(*,*) "src_grid_size:    ", map%src_grid_size 
        ! write(*,*) "dst_grid_size:    ", map%dst_grid_size 
        ! write(*,*) "dst_grid_corners: ", map%dst_grid_corners 
        ! write(*,*) "src_grid_rank:    ", map%src_grid_rank 
        ! write(*,*) "dst_grid_rank:    ", map%dst_grid_rank 
        ! write(*,*) "num_links:        ", map%num_links 
        ! write(*,*) "num_wgts:         ", map%num_wgts 
        
        ! Allocate map_scrip to match dimensions 
        call map_scrip_alloc(map)

        ! Load map from file 
        ! Note: it seems dst_grid_corner_lat and dst_grid_corner_lon
        ! do not exist in all map files generated by cdo. This may be 
        ! related to the type of grid, but it is unclear. Nonetheless,
        ! these variables are so far not used in the map_field routine above.
        ! Below they are only read-in if available. 
        call nc_read(map%map_fname,"src_grid_dims",map%src_grid_dims)
        call nc_read(map%map_fname,"dst_grid_dims",map%dst_grid_dims)
        call nc_read(map%map_fname,"src_grid_center_lat",map%src_grid_center_lat)
        call nc_read(map%map_fname,"dst_grid_center_lat",map%dst_grid_center_lat)
        call nc_read(map%map_fname,"src_grid_center_lon",map%src_grid_center_lon)
        call nc_read(map%map_fname,"dst_grid_center_lon",map%dst_grid_center_lon)
        if (nc_exists_var(map%map_fname,"dst_grid_corner_lat")) then 
            call nc_read(map%map_fname,"dst_grid_corner_lat",map%dst_grid_corner_lat)
        end if 
        if (nc_exists_var(map%map_fname,"dst_grid_corner_lon")) then 
            call nc_read(map%map_fname,"dst_grid_corner_lon",map%dst_grid_corner_lon)
        end if
        call nc_read(map%map_fname,"src_grid_imask",map%src_grid_imask)
        call nc_read(map%map_fname,"dst_grid_imask",map%dst_grid_imask)
        if (nc_exists_var(map%map_fname,"src_grid_area")) then 
          call nc_read(map%map_fname,"src_grid_area",map%src_grid_area)
        endif
        if (nc_exists_var(map%map_fname,"dst_grid_area")) then 
          call nc_read(map%map_fname,"dst_grid_area",map%dst_grid_area)
        endif
        call nc_read(map%map_fname,"src_grid_frac",map%src_grid_frac)
        call nc_read(map%map_fname,"dst_grid_frac",map%dst_grid_frac)
        call nc_read(map%map_fname,"src_address",map%src_address)
        call nc_read(map%map_fname,"dst_address",map%dst_address)
        call nc_read(map%map_fname,"remap_matrix",map%remap_matrix)

        !write(*,*) "range(remap_matrix): ", minval(map%remap_matrix), maxval(map%remap_matrix)
        !write(*,*) "Loaded SCRIP weights."

        ! Summary print line
        !write(*,*) "Loaded SCRIP map from file: "//trim(filename) 
        write(*,*) "Loaded "//trim(method)//" SCRIP map: "//trim(map%src_name)//" => "//trim(map%dst_name) 

        return 

    end subroutine map_scrip_load

    subroutine map_scrip_end(map)
        ! Allocate arrays in map_scrip_class. This
        ! routine assumes that size parameters within
        ! the object are already specified. 

        implicit none 

        type(map_scrip_class), intent(INOUT) :: map 
        
        ! Simply deallocate the map object 
        call map_scrip_dealloc(map) 

        return 

    end subroutine map_scrip_end

    subroutine map_scrip_alloc(map)
        ! Allocate arrays in map_scrip_class. This
        ! routine assumes that size parameters within
        ! the object are already specified. 

        implicit none 

        type(map_scrip_class), intent(INOUT) :: map 

        ! First deallocate everything for safety 

        call map_scrip_dealloc(map) 

        ! Proceed to allocation 

        allocate(map%src_grid_dims(map%src_grid_rank))
        allocate(map%dst_grid_dims(map%dst_grid_rank))
        allocate(map%src_grid_center_lat(map%src_grid_size))
        allocate(map%dst_grid_center_lat(map%dst_grid_size))
        allocate(map%src_grid_center_lon(map%src_grid_size))
        allocate(map%dst_grid_center_lon(map%dst_grid_size))
        if (map%dst_grid_corners>0) then
          allocate(map%dst_grid_corner_lat(map%dst_grid_corners,map%dst_grid_size))
          allocate(map%dst_grid_corner_lon(map%dst_grid_corners,map%dst_grid_size))
        endif
        allocate(map%src_grid_imask(map%src_grid_size))
        allocate(map%dst_grid_imask(map%dst_grid_size))
        allocate(map%src_grid_area(map%src_grid_size))
        allocate(map%dst_grid_area(map%dst_grid_size))
        allocate(map%src_grid_frac(map%src_grid_size))
        allocate(map%dst_grid_frac(map%dst_grid_size))
        allocate(map%src_address(map%num_links))
        allocate(map%dst_address(map%num_links))
        allocate(map%remap_matrix(map%num_wgts,map%num_links))

        return 

    end subroutine map_scrip_alloc

    subroutine map_scrip_dealloc(map)
        ! Allocate arrays in map_scrip_class. This
        ! routine assumes that size parameters within
        ! the object are already specified. 

        implicit none 

        type(map_scrip_class), intent(INOUT) :: map 
        
        if (allocated(map%src_grid_dims))       deallocate(map%src_grid_dims)
        if (allocated(map%dst_grid_dims))       deallocate(map%dst_grid_dims)
        if (allocated(map%src_grid_center_lat)) deallocate(map%src_grid_center_lat)
        if (allocated(map%dst_grid_center_lat)) deallocate(map%dst_grid_center_lat)
        if (allocated(map%src_grid_center_lon)) deallocate(map%src_grid_center_lon)
        if (allocated(map%dst_grid_center_lon)) deallocate(map%dst_grid_center_lon)
        if (allocated(map%dst_grid_corner_lat)) deallocate(map%dst_grid_corner_lat)
        if (allocated(map%dst_grid_corner_lon)) deallocate(map%dst_grid_corner_lon)
        if (allocated(map%src_grid_imask))      deallocate(map%src_grid_imask)
        if (allocated(map%dst_grid_imask))      deallocate(map%dst_grid_imask)
        if (allocated(map%src_grid_area))       deallocate(map%src_grid_area)
        if (allocated(map%dst_grid_area))       deallocate(map%dst_grid_area)
        if (allocated(map%src_grid_frac))       deallocate(map%src_grid_frac)
        if (allocated(map%dst_grid_frac))       deallocate(map%dst_grid_frac)
        if (allocated(map%src_address))         deallocate(map%src_address)
        if (allocated(map%dst_address))         deallocate(map%dst_address)
        if (allocated(map%remap_matrix))        deallocate(map%remap_matrix)
        
        return 

    end subroutine map_scrip_dealloc

    function gen_map_filename(src_name,dst_name,fldr,method) result(filename)
        ! Output the standard map filename with input folder name
        implicit none 

        character(len=*), intent(IN) :: src_name
        character(len=*), intent(IN) :: dst_name 
        character(len=*), intent(IN) :: fldr 
        character(len=*), intent(IN) :: method
        character(len=256) :: filename

        filename = trim(fldr)//"/scrip-"//trim(method)//"_"//trim(src_name)//"_"//trim(dst_name)//".nc"

        return

    end function gen_map_filename


! === NCIO Extension functions ===

    subroutine nc_read_interp_dp_2D(filename,vnm,var2D,var2D_in,ncid,start,count,mps,method, &
                                        fill_method,filt_method,filt_par)
        ! Read in a 2D field and interpolate it using scrip map
        ! as needed. 

        implicit none

        character(len=*),   intent(IN)  :: filename 
        character(len=*),   intent(IN)  :: vnm 
        real(dp),           intent(OUT) :: var2D(:,:) 
        real(dp),optional,  intent(IN)  :: var2D_in(:,:) 
        integer, optional,  intent(IN)  :: ncid 
        integer, optional,  intent(IN)  :: start(:) 
        integer, optional,  intent(IN)  :: count(:) 
        type(map_scrip_class), intent(IN),  optional :: mps
        character(len=*),      intent(IN),  optional :: method 
        character(len=*),      intent(IN),  optional :: fill_method     ! Method to fill in remaining missing values
        character(len=*),      intent(IN),  optional :: filt_method     ! Method to use for filtering
        real(dp),              intent(IN),  optional :: filt_par(:)     ! gaussian=[sigma,dx]; poisson=[tol]
        
        ! Local variables 
        integer :: nx, ny 
        integer,  allocatable :: dims(:) 
        real(dp), allocatable :: var2D_src(:,:) 
        character(len=56) :: mapping_method 

        if (present(var2D_in)) then 
            ! Get source array from argument (useful for handling 3D arrays with this routine)

            nx = size(var2D_in,1)
            ny = size(var2D_in,2) 

            allocate(var2D_src(nx,ny))

            var2D_src = var2D_in 

        else 
            ! Load 2D array from file 

            ! Determine dimensions of current variable in the source file
            call nc_dims(filename,vnm,dims=dims)

            nx = dims(1)
            ny = dims(2) 

            allocate(var2D_src(nx,ny))

            ! Load the variable from the file to the local 2D array
            call nc_read(filename,vnm,var2D_src,ncid=ncid,start=start,count=[nx,ny,1],missing_value=real(MV,dp))
        
        end if 


        if (nx .eq. size(var2D,1) .and. ny .eq. size(var2D,2)) then 
            ! Assume no interpolation needed, copy variable for output directly

            var2D = var2D_src 

        else 
            ! Map local source array it to our target array 

            ! Determine mapping method for this variable 
            mapping_method = "mean"
            if (present(method)) mapping_method = trim(method) 

            ! Safety check 
            if (.not. present(mps)) then 
                write(error_unit,*) ""
                write(error_unit,*) "nc_read_interp:: Error: map_scrip_class object must &
                        &be provided as an argument since array read from file does not &
                        &match the Yelmo array size."
                write(error_unit,*) "filename: ", trim(filename)
                write(error_unit,*) "variable: ", trim(vnm)
                write(error_unit,*) "dims in file:         ", nx, ny 
                write(error_unit,*) "dims in yelmo object: ", size(var2D,1), size(var2D,2)
                stop 
            end if 

            ! Perform conservative interpolation 
            var2D = real(MV,dp) 
            call map_scrip_field(mps,vnm,var2D_src,var2D,method=mapping_method,missing_value=real(MV,dp), &
                        fill_method=fill_method,filt_method=filt_method,filt_par=filt_par)

        end if 

        return

    end subroutine nc_read_interp_dp_2D

    subroutine nc_read_interp_dp_3D(filename,vnm,var3D,ncid,start,count,mps,method, &
                                        fill_method,filt_method,filt_par)
        ! Read in a 3D field and interpolate it using scrip map
        ! as needed. 

        implicit none

        character(len=*),   intent(IN)  :: filename 
        character(len=*),   intent(IN)  :: vnm 
        real(dp),           intent(OUT) :: var3D(:,:,:) 
        integer, optional,  intent(IN)  :: ncid 
        integer, optional,  intent(IN)  :: start(:) 
        integer, optional,  intent(IN)  :: count(:) 
        type(map_scrip_class), intent(IN),  optional :: mps
        character(len=*),      intent(IN),  optional :: method 
        character(len=*),      intent(IN),  optional :: fill_method     ! Method to fill in remaining missing values
        character(len=*),      intent(IN),  optional :: filt_method     ! Method to use for filtering
        real(dp),              intent(IN),  optional :: filt_par(:)     ! gaussian=[sigma,dx]; poisson=[tol]
        
        ! Local variables 
        integer :: nx, ny, nz, k  
        integer,  allocatable :: dims(:) 
        real(dp), allocatable :: var3D_src(:,:,:) 

        ! Determine dimensions of current variable in the source file
        call nc_dims(filename,vnm,dims=dims)

        nx = dims(1)
        ny = dims(2) 
        nz = dims(3) 

        allocate(var3D_src(nx,ny,nz))

        ! Safety check 
        if (nz .ne. size(var3D,3)) then 

            write(error_unit,*) ""
            write(error_unit,*) "nc_read_interp_dp_3D:: Error: third dimension of variable in &
                    &input file does not match third dimension of array provided as an argument. &
                    &Interpolation of the third dimension is not yet supported."
            write(error_unit,*) "filename  = ", trim(filename)
            write(error_unit,*) "variable  = ", trim(vnm)
            write(error_unit,*) "nz[file]  = ", nz
            write(error_unit,*) "nz[array] = ", size(var3D,3) 
            stop 
        end if 
            
        ! Read in full 3D variable of interest 
        call nc_read(filename,vnm,var3D_src,ncid=ncid,start=start,count=[nx,ny,nz,1],missing_value=real(MV,dp))
        
        ! Loop over third dimension and apply interpolation 
        do k = 1, nz 
            
            call nc_read_interp_dp_2D(filename,vnm,var3D(:,:,k),var3D_src(:,:,k),ncid,start,count, &
                                mps,method,fill_method=fill_method,filt_method=filt_method,filt_par=filt_par)

        end do 

        return

    end subroutine nc_read_interp_dp_3D

    subroutine nc_read_interp_sp_2D(filename,vnm,var2D,var2D_in,ncid,start,count,mps,method, &
                                        fill_method,filt_method,filt_par)
        ! Read in a 2D field and interpolate it using scrip map
        ! as needed. 

        implicit none

        character(len=*),   intent(IN)  :: filename 
        character(len=*),   intent(IN)  :: vnm 
        real(sp),           intent(OUT) :: var2D(:,:) 
        real(sp),optional,  intent(IN)  :: var2D_in(:,:) 
        integer, optional,  intent(IN)  :: ncid 
        integer, optional,  intent(IN)  :: start(:) 
        integer, optional,  intent(IN)  :: count(:) 
        type(map_scrip_class), intent(IN),  optional :: mps
        character(len=*),      intent(IN),  optional :: method 
        character(len=*),      intent(IN),  optional :: fill_method     ! Method to fill in remaining missing values
        character(len=*),      intent(IN),  optional :: filt_method     ! Method to use for filtering
        real(sp),              intent(IN),  optional :: filt_par(:)     ! gaussian=[sigma,dx]; poisson=[tol]
        
        ! Local variables 
        integer :: nx, ny 
        integer,  allocatable :: dims(:) 
        real(sp), allocatable :: var2D_src(:,:) 
        character(len=56) :: mapping_method 

        if (present(var2D_in)) then 
            ! Get source array from argument (useful for handling 3D arrays with this routine)

            nx = size(var2D_in,1)
            ny = size(var2D_in,2) 

            allocate(var2D_src(nx,ny))

            var2D_src = var2D_in 

        else 
            ! Load 2D array from file 

            ! Determine dimensions of current variable in the source file
            call nc_dims(filename,vnm,dims=dims)

            nx = dims(1)
            ny = dims(2) 

            allocate(var2D_src(nx,ny))

            ! Load the variable from the file to the local 2D array
            call nc_read(filename,vnm,var2D_src,ncid=ncid,start=start,count=[nx,ny,1],missing_value=real(MV,wp))
        
        end if 


        if (nx .eq. size(var2D,1) .and. ny .eq. size(var2D,2)) then 
            ! Assume no interpolation needed, copy variable for output directly

            var2D = var2D_src 

        else 
            ! Map local source array it to our target array 

            ! Determine mapping method for this variable 
            mapping_method = "mean"
            if (present(method)) mapping_method = trim(method) 

            ! Safety check 
            if (.not. present(mps)) then 
                write(error_unit,*) ""
                write(error_unit,*) "nc_read_interp:: Error: map_scrip_class object must &
                        &be provided as an argument since array read from file does not &
                        &match the Yelmo array size."
                write(error_unit,*) "filename: ", trim(filename)
                write(error_unit,*) "variable: ", trim(vnm)
                write(error_unit,*) "dims in file:         ", nx, ny 
                write(error_unit,*) "dims in yelmo object: ", size(var2D,1), size(var2D,2)
                stop 
            end if 

            ! Perform conservative interpolation 
            var2D = real(MV,sp) 
            call map_scrip_field(mps,vnm,var2D_src,var2D,method=mapping_method,missing_value=real(MV,sp), &
                            fill_method=fill_method,filt_method=filt_method,filt_par=filt_par)

        end if 

        return

    end subroutine nc_read_interp_sp_2D

    subroutine nc_read_interp_sp_3D(filename,vnm,var3D,ncid,start,count,mps,method, &
                                        fill_method,filt_method,filt_par)
        ! Read in a 2D field and interpolate it using scrip map
        ! as needed. 

        implicit none

        character(len=*),   intent(IN)  :: filename 
        character(len=*),   intent(IN)  :: vnm 
        real(sp),           intent(OUT) :: var3D(:,:,:) 
        integer, optional,  intent(IN)  :: ncid 
        integer, optional,  intent(IN)  :: start(:) 
        integer, optional,  intent(IN)  :: count(:) 
        type(map_scrip_class), intent(IN),  optional :: mps
        character(len=*),      intent(IN),  optional :: method 
        character(len=*),      intent(IN),  optional :: fill_method     ! Method to fill in remaining missing values
        character(len=*),      intent(IN),  optional :: filt_method     ! Method to use for filtering
        real(sp),              intent(IN),  optional :: filt_par(:)     ! gaussian=[sigma,dx]; poisson=[tol]
        
        ! Local variables 
        integer :: nx, ny, nz, k  
        integer,  allocatable :: dims(:) 
        real(sp), allocatable :: var3D_src(:,:,:) 

        ! Determine dimensions of current variable in the source file
        call nc_dims(filename,vnm,dims=dims)

        nx = dims(1)
        ny = dims(2) 
        nz = dims(3) 

        allocate(var3D_src(nx,ny,nz))

        ! Safety check 
        if (nz .ne. size(var3D,3)) then 

            write(error_unit,*) ""
            write(error_unit,*) "nc_read_interp_sp_3D:: Error: vertical dimension of variable in &
                    &input file does not match vertical dimension of yelmo object. Vertical &
                    & interpolation is not yet supported."
            write(error_unit,*) "filename  = ", trim(filename)
            write(error_unit,*) "variable  = ", trim(vnm)
            write(error_unit,*) "nz[file]  = ", nz
            write(error_unit,*) "nz[yelmo] = ", size(var3D,3) 
            stop 
        end if 
            
        ! Read in full 3D variable of interest 
        call nc_read(filename,vnm,var3D_src,ncid=ncid,start=start,count=[nx,ny,nz,1],missing_value=real(MV,sp))
        
        ! Loop over vertical dimension and apply interpolation 
        do k = 1, nz 
            
            call nc_read_interp_sp_2D(filename,vnm,var3D(:,:,k),var3D_src(:,:,k),ncid,start,count,mps,method, &
                                fill_method=fill_method,filt_method=filt_method,filt_par=filt_par)

        end do 

        return

    end subroutine nc_read_interp_sp_3D

    subroutine nc_read_interp_int_2D(filename,vnm,var2D,var2D_in,ncid,start,count,mps,method, &
                                        fill_method)
        ! Read in a 2D field and interpolate it using scrip map
        ! as needed. 

        implicit none

        character(len=*),   intent(IN)  :: filename 
        character(len=*),   intent(IN)  :: vnm 
        integer,            intent(OUT) :: var2D(:,:) 
        integer, optional,  intent(IN)  :: var2D_in(:,:) 
        integer, optional,  intent(IN)  :: ncid 
        integer, optional,  intent(IN)  :: start(:) 
        integer, optional,  intent(IN)  :: count(:) 
        type(map_scrip_class), intent(IN),  optional :: mps
        character(len=*),      intent(IN),  optional :: method 
        character(len=*),      intent(IN),  optional :: fill_method     ! Method to fill in remaining missing values

        ! Local variables 
        integer :: nx, ny 
        integer,  allocatable :: dims(:) 
        integer,  allocatable :: var2D_src(:,:) 
        character(len=56) :: mapping_method 

        if (present(var2D_in)) then 
            ! Get source array from argument (useful for handling 3D arrays with this routine)

            nx = size(var2D_in,1)
            ny = size(var2D_in,2) 

            allocate(var2D_src(nx,ny))

            var2D_src = var2D_in 

        else 
            ! Load 2D array from file 

            ! Determine dimensions of current variable in the source file
            call nc_dims(filename,vnm,dims=dims)

            nx = dims(1)
            ny = dims(2) 

            allocate(var2D_src(nx,ny))

            ! Load the variable from the file to the local 2D array
            call nc_read(filename,vnm,var2D_src,ncid=ncid,start=start,count=[nx,ny,1],missing_value=int(MV))
        
        end if 


        if (nx .eq. size(var2D,1) .and. ny .eq. size(var2D,2)) then 
            ! Assume no interpolation needed, copy variable for output directly

            var2D = var2D_src 

        else 
            ! Map local source array it to our target array 

            ! Determine mapping method for this variable 
            mapping_method = "mean"
            if (present(method)) mapping_method = trim(method) 

            ! Perform conservative interpolation 
            var2D = int(MV)
            call map_scrip_field(mps,vnm,var2D_src,var2D,method=mapping_method, &
                                    missing_value=int(MV),fill_method=fill_method)

        end if 

        return

    end subroutine nc_read_interp_int_2D

    subroutine nc_read_interp_int_3D(filename,vnm,var3D,ncid,start,count,mps,method,fill_method)
        ! Read in a 2D field and interpolate it using scrip map
        ! as needed. 

        implicit none

        character(len=*),   intent(IN)  :: filename 
        character(len=*),   intent(IN)  :: vnm 
        integer,            intent(OUT) :: var3D(:,:,:) 
        integer, optional,  intent(IN)  :: ncid 
        integer, optional,  intent(IN)  :: start(:) 
        integer, optional,  intent(IN)  :: count(:) 
        type(map_scrip_class), intent(IN),  optional :: mps
        character(len=*),      intent(IN),  optional :: method 
        character(len=*),      intent(IN),  optional :: fill_method     ! Method to fill in remaining missing values

        ! Local variables 
        integer :: nx, ny, nz, k  
        integer,  allocatable :: dims(:) 
        integer,  allocatable :: var3D_src(:,:,:) 

        ! Determine dimensions of current variable in the source file
        call nc_dims(filename,vnm,dims=dims)

        nx = dims(1)
        ny = dims(2) 
        nz = dims(3) 

        allocate(var3D_src(nx,ny,nz))

        ! Safety check 
        if (nz .ne. size(var3D,3)) then 

            write(error_unit,*) ""
            write(error_unit,*) "nc_read_interp_int_3D:: Error: vertical dimension of variable in &
                    &input file does not match vertical dimension of yelmo object. Vertical &
                    & interpolation is not yet supported."
            write(error_unit,*) "filename  = ", trim(filename)
            write(error_unit,*) "variable  = ", trim(vnm)
            write(error_unit,*) "nz[file]  = ", nz
            write(error_unit,*) "nz[yelmo] = ", size(var3D,3) 
            stop 
        end if 
            
        ! Read in full 3D variable of interest 
        call nc_read(filename,vnm,var3D_src,ncid=ncid,start=start,count=[nx,ny,nz,1],missing_value=int(MV))
        
        ! Loop over vertical dimension and apply interpolation 
        do k = 1, nz 
            
            call nc_read_interp_int_2D(filename,vnm,var3D(:,:,k),var3D_src(:,:,k),ncid,start,count,mps,method, &
                                fill_method=fill_method)

        end do 

        return

    end subroutine nc_read_interp_int_3D

    subroutine nc_read_interp_logical_2D(filename,vnm,var2D,var2D_in,ncid,start,count,mps,method)
        ! Read in a 2D field and interpolate it using scrip map
        ! as needed. 

        implicit none

        character(len=*),   intent(IN)  :: filename 
        character(len=*),   intent(IN)  :: vnm 
        logical,            intent(OUT) :: var2D(:,:) 
        logical, optional,  intent(IN)  :: var2D_in(:,:) 
        integer, optional,  intent(IN)  :: ncid 
        integer, optional,  intent(IN)  :: start(:) 
        integer, optional,  intent(IN)  :: count(:) 
        type(map_scrip_class), intent(IN),  optional :: mps
        character(len=*),      intent(IN),  optional :: method 

        ! Local variables 
        integer :: nx, ny 
        integer,  allocatable :: dims(:) 
        integer,  allocatable :: var2D_src(:,:) 
        integer,  allocatable :: var2D_int(:,:) 
        character(len=56) :: mapping_method 

        if (present(var2D_in)) then 
            ! Get source array from argument (useful for handling 3D arrays with this routine)

            nx = size(var2D_in,1)
            ny = size(var2D_in,2) 

            allocate(var2D_src(nx,ny))

            var2D_src = 0 
            where (var2D_in) var2D_src = 1 

        else 
            ! Load 2D array from file 

            ! Determine dimensions of current variable in the source file
            call nc_dims(filename,vnm,dims=dims)

            nx = dims(1)
            ny = dims(2) 

            allocate(var2D_src(nx,ny))

            ! Load the variable from the file to the local 2D array
            call nc_read(filename,vnm,var2D_src,ncid=ncid,start=start,count=[nx,ny,1],missing_value=int(MV))
        
        end if 


        if (nx .eq. size(var2D,1) .and. ny .eq. size(var2D,2)) then 
            ! Assume no interpolation needed, copy variable for output directly

            var2D = .FALSE. 
            where(var2D_src .eq. 1) var2D = .TRUE.  

        else 
            ! Map local source array it to our target array 

            ! Determine mapping method for this variable 
            mapping_method = "count"
            if (present(method)) mapping_method = trim(method) 

            ! Perform conservative interpolation 
            allocate(var2D_int(size(var2D,1),size(var2D,2)))
            var2D_int = int(MV)
            call map_scrip_field(mps,vnm,var2D_src,var2D_int,method=mapping_method, &
                                                missing_value=int(MV),fill_method="nn")

            var2D = .FALSE. 
            where(var2D_int .eq. 1) var2D = .TRUE.  

        end if 

        return

    end subroutine nc_read_interp_logical_2D

end module coordinates_mapping_scrip
