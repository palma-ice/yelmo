
module yelmo_defs
    
    use nml 

    implicit none 

    ! =========================================================================
    !
    ! CONSTANTS (program precision, global constants)
    !
    ! =========================================================================

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the precision of the library (sp,dp)
    integer,  parameter :: prec = sp 

    ! Define more conventional wp = working precision 
    ! ajr: slowly transition from 'prec' to 'wp' 
    integer,  parameter :: wp = prec 

    ! Missing value and aliases
    real(wp), parameter :: MISSING_VALUE_DEFAULT = real(-9999.0,wp)
    real(wp), parameter :: MISSING_VALUE = MISSING_VALUE_DEFAULT
    real(wp), parameter :: MV = MISSING_VALUE_DEFAULT
    
    ! Error distance (very large), error index, and smallest number epsilon 
    real(wp), parameter :: ERR_DIST = real(1E8,wp) 
    integer,    parameter :: ERR_IND  = -1 
    real(wp), parameter :: tol_underflow = real(1E-15,wp)

    ! Mathematical constants
    real(wp), parameter :: pi  = real(2._dp*acos(0.0_dp),wp)
    real(wp), parameter :: degrees_to_radians = real(pi / 180._dp,wp)  ! Conversion factor between radians and degrees
    real(wp), parameter :: radians_to_degrees = real(180._dp / pi,wp)  ! Conversion factor between degrees and radians
    
    ! The constants below should be loaded using the global subroutine
    ! defined below `yelmo_constants_load`.
    ! Note: The key limitation imposed by defining the parameters defined 
    ! globally is that these constants must be the same for all domains 
    ! being run in the same program. 

    ! Yelmo configuration options 
    logical :: yelmo_log
    logical :: yelmo_use_omp 

    ! Physical constants 
    real(wp) :: sec_year       ! [s] seconds per year 
    real(wp) :: g              ! [m s-2] Gravitational accel.  
    real(wp) :: T0             ! [K] Reference freezing temperature  
    real(wp) :: rho_ice        ! [kg m-3] Density ice           
    real(wp) :: rho_w          ! [kg m-3] Density water          
    real(wp) :: rho_sw         ! [kg m-3] Density seawater      
    real(wp) :: rho_a          ! [kg m-3] Density asthenosphere  
    real(wp) :: rho_m          ! [kg m-3] Density mantle (lith) 
    real(wp) :: L_ice          ! [J kg-1] Latent heat           
    real(wp) :: T_pmp_beta     ! [K Pa-1] Melt point pressure slope

    ! Internal parameters 
    real(wp) :: conv_we_ie        ! Conversion water equiv. => m/a ice equiv. 
    real(wp) :: conv_mmdwe_maie   ! Conversion mm/d water equiv. => m/a ice equiv.
    real(wp) :: conv_mmawe_maie   ! Conversion mm/a water equiv. => m/a ice equiv. 
    
    ! =========================================================================
    !
    ! YELMO objects: ytopo 
    !
    ! =========================================================================
    
    ! ytopo parameters
    type ytopo_param_class
        character(len=256) :: solver
        character(len=256) :: calv_method  
        integer            :: surf_gl_method 
        logical            :: margin2nd 
        logical            :: use_bmb  
        logical            :: topo_fixed
        integer            :: topo_rel
        real(wp)         :: topo_rel_tau 
        real(wp)         :: calv_H_lim
        real(wp)         :: calv_tau  
        real(wp)         :: H_min_grnd
        real(wp)         :: H_min_flt 
        real(wp)         :: sd_min 
        real(wp)         :: sd_max 
        real(wp)         :: calv_max  
        real(wp)         :: grad_lim 
        integer            :: gl_sep 
        integer            :: gl_sep_nx 
        logical            :: diffuse_bmb_shlf 

        ! Internal parameters 
        real(dp)           :: time 
        integer            :: nx, ny
        real(wp)         :: dx, dy
        character(len=256) :: boundaries 
        
        character(len=256) :: pc_step 
        real(wp) :: dt_zeta, dt_beta(4)

        real(wp) :: speed, speed_pred, speed_corr 

    end type

    ! ytopo state variables
    type ytopo_state_class
        ! Model variables that the define the state of the domain 

        real(wp), allocatable :: H_ice(:,:)      ! Ice thickness [m] 
        real(wp), allocatable :: z_srf(:,:)      ! Surface elevation [m]
        real(wp), allocatable :: dzsrfdt(:,:)    ! Surface elevation rate of change [m/a] 
        real(wp), allocatable :: dHicedt(:,:)    ! Ice thickness rate of change [m/a] 
        real(wp), allocatable :: bmb(:,:)        ! Combined field of bmb_grnd and bmb_shlf 
        real(wp), allocatable :: mb_applied(:,:) ! Actual mass balance applied [m/a], for mass balance accounting
        real(wp), allocatable :: calv_grnd(:,:)  ! Grounded calving [m/a]
        real(wp), allocatable :: calv(:,:)       ! Calving [m/a]
        
        real(wp), allocatable :: H_margin(:,:)   ! [m] Margin ice thickness in partially filled cells 

        real(wp), allocatable :: dzsdx(:,:)      ! Surface elevation slope [m m-1], Ac x nodes
        real(wp), allocatable :: dzsdy(:,:)      ! Surface elevation slope [m m-1], Ac y nodes
        real(wp), allocatable :: dHicedx(:,:)    ! Ice thickness gradient slope [m m-1], Ac x nodes
        real(wp), allocatable :: dHicedy(:,:)    ! Ice thickness gradient slope [m m-1], Ac y nodes
        
        real(wp), allocatable :: H_grnd(:,:)       ! Ice thickness overburden [m]
        
        ! Masks 
        real(wp), allocatable :: f_grnd(:,:)       ! Grounded fraction (grounding line fraction between 0 and 1)
        real(wp), allocatable :: f_grnd_acx(:,:)   ! Grounded fraction (acx nodes)
        real(wp), allocatable :: f_grnd_acy(:,:)   ! Grounded fraction (acy nodes)
        real(wp), allocatable :: f_ice(:,:)        ! Ice-covered fraction 

        real(wp), allocatable :: dist_margin(:,:)  ! Distance to nearest margin point 
        real(wp), allocatable :: dist_grline(:,:)  ! Distance to nearest grounding-line point 
        
        ! Additional masks 
        integer,    allocatable :: mask_bed(:,:)    ! Multi-valued bed mask
        logical,    allocatable :: is_grline(:,:)   ! Grounding line points
        logical,    allocatable :: is_grz(:,:)      ! Grounding line plus grounded neighbors
        
        real(wp), allocatable :: dHdt_n(:,:)      ! [m/a] Ice thickness change due to advection only
        real(wp), allocatable :: H_ice_n(:,:)     ! [m] Ice thickness from the previous timestep 
        real(wp), allocatable :: H_ice_pred(:,:)  ! [m] Ice thickness, predicted, for time=n+1
        real(wp), allocatable :: H_ice_corr(:,:)  ! [m] Ice thickness, corrected, for time=n+1 
        
        real(wp), allocatable :: z_srf_n(:,:)     ! [m] Surface elevation from the previous timestep 
        
    end type

    ! ytopo class
    type ytopo_class

        type(ytopo_param_class) :: par          ! Parameters
        type(ytopo_state_class) :: now          ! Variables

    end type

    ! =========================================================================
    !
    ! YELMO objects: ydyn 
    !
    ! =========================================================================
    
    ! ydyn parameters
    type ydyn_param_class

        character(len=256) :: solver 
        integer    :: visc_method 
        real(wp) :: visc_const 
        integer    :: beta_method
        real(wp) :: beta_const
        real(wp) :: beta_q                ! Friction law exponent
        real(wp) :: beta_u0               ! [m/a] Friction law velocity threshold 
        integer    :: beta_gl_scale         ! Beta grounding-line scaling method (beta => 0 at gl?)
        integer    :: beta_gl_sep           ! Beta grounding-line sub-element (subgrid) parameterization
        integer    :: beta_gl_stag          ! Beta grounding-line staggering method 
        real(wp) :: beta_gl_f             ! Fraction of beta at gl 
        integer    :: taud_gl_method        ! Driving stress grounding line treatment 
        real(wp) :: H_grnd_lim 
        real(wp) :: H_sed_sat
        integer    :: cb_method
        logical    :: cb_with_pmp           ! Scale friction coefficient between frozen and streaming values?
        logical    :: cb_margin_pmp         ! Ensure margin and grline are considered streaming?
        character(len=256) :: cb_scale      ! Method for scaling cb with elevation
        real(wp) :: cb_z0  
        real(wp) :: cb_z1
        real(wp) :: cb_min      
        real(wp) :: cf_frozen
        real(wp) :: cf_stream
        integer    :: n_sm_beta 
        real(wp) :: beta_min              ! Minimum allowed value of beta
        real(wp) :: eps_0                 ! Minimum assumed strain rate for effective viscosity regularization
        character(len=256) :: ssa_lis_opt 
        real(wp) :: ssa_beta_max          ! Maximum value of beta for which ssa should be calculated
        real(wp) :: ssa_vel_max
        integer    :: ssa_iter_max 
        real(wp) :: ssa_iter_rel 
        real(wp) :: ssa_iter_conv 
        real(wp) :: taud_lim 
        real(wp) :: cb_sia
        
        integer    :: neff_method
        real(wp) :: neff_const
        real(wp) :: neff_p 
        logical    :: neff_set_water 
        real(wp) :: neff_w_max
        real(wp) :: neff_N0
        real(wp) :: neff_delta 
        real(wp) :: neff_e0 
        real(wp) :: neff_Cc 

        real(wp) :: till_phi_const 
        real(wp) :: till_phi_min 
        real(wp) :: till_phi_max 
        real(wp) :: till_phi_zmin 
        real(wp) :: till_phi_zmax 
        
        ! Internal parameters 
        character(len=256) :: boundaries 
        logical    :: use_ssa                   ! Should ssa be used? 
        logical    :: use_bmb                   ! Set to match `use_bmb` in ytopo_param_class 
        integer    :: nx, ny, nz_aa, nz_ac 
        real(wp) :: dx, dy
        real(wp), allocatable :: zeta_aa(:)   ! Layer centers (aa-nodes), plus base and surface: nz_aa points 
        real(wp), allocatable :: zeta_ac(:)   ! Layer borders (ac-nodes), plus base and surface: nz_ac == nz_aa-1 points
        real(dp)   :: time

        integer    :: ssa_iter_now              ! Number of iterations used for Picard iteration to solve ssa this timestep
        real(wp) :: speed 

    end type

    ! ydyn state variables
    type ydyn_state_class
        ! Model variables that the define the state of the domain 

        real(wp), allocatable :: ux(:,:,:) 
        real(wp), allocatable :: uy(:,:,:) 
        real(wp), allocatable :: uxy(:,:,:)
        real(wp), allocatable :: uz(:,:,:)  

        real(wp), allocatable :: ux_bar(:,:) 
        real(wp), allocatable :: uy_bar(:,:)
        real(wp), allocatable :: uxy_bar(:,:)

        real(wp), allocatable :: ux_bar_prev(:,:) 
        real(wp), allocatable :: uy_bar_prev(:,:)

        real(wp), allocatable :: ux_b(:,:) 
        real(wp), allocatable :: uy_b(:,:)
        real(wp), allocatable :: uxy_b(:,:)

        ! Surface velocity: eventually these could be pointers since it is simply
        ! the top layer in ux(:,:,:), etc. and only used, not calculated.
        real(wp), allocatable :: ux_s(:,:) 
        real(wp), allocatable :: uy_s(:,:)
        real(wp), allocatable :: uxy_s(:,:)
        
        real(wp), allocatable :: ux_i(:,:,:) 
        real(wp), allocatable :: uy_i(:,:,:)
        real(wp), allocatable :: ux_i_bar(:,:) 
        real(wp), allocatable :: uy_i_bar(:,:)
        real(wp), allocatable :: uxy_i_bar(:,:) 
        
        real(wp), allocatable :: duxydt(:,:) 

        real(wp), allocatable :: dd_ab(:,:,:)  
        real(wp), allocatable :: dd_ab_bar(:,:)  
        
        real(wp), allocatable :: sigma_horiz_sq(:,:)
        real(wp), allocatable :: lhs_x(:,:) 
        real(wp), allocatable :: lhs_y(:,:) 
        real(wp), allocatable :: lhs_xy(:,:) 

        real(wp), allocatable :: duxdz(:,:,:) 
        real(wp), allocatable :: duydz(:,:,:)
        real(wp), allocatable :: duxdz_bar(:,:) 
        real(wp), allocatable :: duydz_bar(:,:)

        real(wp), allocatable :: taud_acx(:,:) 
        real(wp), allocatable :: taud_acy(:,:) 
        real(wp), allocatable :: taud(:,:) 
        
        real(wp), allocatable :: taub_acx(:,:) 
        real(wp), allocatable :: taub_acy(:,:) 
        real(wp), allocatable :: taub(:,:)
        
        real(wp), allocatable :: qq_gl_acx(:,:) 
        real(wp), allocatable :: qq_gl_acy(:,:) 
        
        real(wp), allocatable :: qq_acx(:,:) 
        real(wp), allocatable :: qq_acy(:,:) 
        real(wp), allocatable :: qq(:,:)
        
        real(wp), allocatable :: visc_eff(:,:,:)
        real(wp), allocatable :: visc_eff_int(:,:)

        real(wp), allocatable :: N_eff(:,:)       ! Effective pressure
        real(wp), allocatable :: cf_ref(:,:)
        real(wp), allocatable :: c_bed(:,:)  
        real(wp), allocatable :: beta_acx(:,:) 
        real(wp), allocatable :: beta_acy(:,:) 
        real(wp), allocatable :: beta(:,:)         
        real(wp), allocatable :: beta_eff(:,:) 

        real(wp), allocatable :: f_vbvs(:,:) 

        integer,    allocatable :: ssa_mask_acx(:,:) 
        integer,    allocatable :: ssa_mask_acy(:,:) 
        real(wp), allocatable :: ssa_err_acx(:,:) 
        real(wp), allocatable :: ssa_err_acy(:,:) 

    end type

    ! ydyn class
    type ydyn_class

        type(ydyn_param_class)    :: par        ! physical parameters
        type(ydyn_state_class)    :: now

    end type

    ! =========================================================================
    !
    ! YELMO objects: ymat 
    !
    ! =========================================================================
    
    type strain_2D_class 
        real(wp), allocatable :: dxx(:,:) 
        real(wp), allocatable :: dyy(:,:) 
        real(wp), allocatable :: dxy(:,:) 
        real(wp), allocatable :: de(:,:) 
    end type 

    type strain_3D_class 
        real(wp), allocatable :: dxx(:,:,:) 
        real(wp), allocatable :: dyy(:,:,:) 
        real(wp), allocatable :: dzz(:,:,:)
        real(wp), allocatable :: dxy(:,:,:) 
        real(wp), allocatable :: dxz(:,:,:) 
        real(wp), allocatable :: dyz(:,:,:) 
        real(wp), allocatable :: de(:,:,:) 
        real(wp), allocatable :: f_shear(:,:,:) 
    end type 
    
    type ymat_param_class
        
        character(len=56)   :: flow_law
        integer             :: rf_method 
        real(wp)          :: rf_const
        logical             :: rf_use_eismint2
        logical             :: rf_with_water 
        real(wp)          :: n_glen                       ! Flow law exponent (n_glen=3)
        real(wp)          :: visc_min  
        real(wp)          :: de_max 
        character(len=56)   :: enh_method  
        real(wp)          :: enh_shear
        real(wp)          :: enh_stream
        real(wp)          :: enh_shlf
        real(wp)          :: enh_umin 
        real(wp)          :: enh_umax
        logical             :: calc_age
        real(wp), allocatable :: age_iso(:)
        character(len=56)   :: tracer_method  
        real(wp)          :: tracer_impl_kappa
        
        ! Internal parameters
        real(dp)   :: time 
        real(wp) :: dx, dy  
        integer    :: nx, ny, nz_aa, nz_ac  
        integer    :: n_iso 

        real(wp) :: speed 

        real(wp), allocatable :: zeta_aa(:)   ! Layer centers (aa-nodes), plus base and surface: nz_aa points 
        real(wp), allocatable :: zeta_ac(:)   ! Layer borders (ac-nodes), plus base and surface: nz_ac == nz_aa-1 points
        
    end type 

    type ymat_state_class 

        type(strain_2D_class)   :: strn2D
        type(strain_3D_class)   :: strn 

        real(wp), allocatable :: enh(:,:,:)
        real(wp), allocatable :: enh_bnd(:,:,:)
        real(wp), allocatable :: enh_bar(:,:)
        real(wp), allocatable :: ATT(:,:,:) 
        real(wp), allocatable :: ATT_bar(:,:)
        real(wp), allocatable :: visc(:,:,:) 
        real(wp), allocatable :: visc_int(:,:) 

        real(wp), allocatable :: f_shear_bar(:,:) 
        
        real(wp), allocatable :: dep_time(:,:,:)      ! Ice deposition time (for online age tracing)
        real(wp), allocatable :: depth_iso(:,:,:)     ! Depth of specific isochronal layers

    end type 

    type ymat_class
        type(ymat_param_class) :: par 
        type(ymat_state_class) :: now 
    end type

    ! =========================================================================
    !
    ! YELMO objects: ytherm 
    !
    ! =========================================================================
    
    !ytherm parameters 
    type ytherm_param_class
        character(len=256)  :: method  
        character(len=256)  :: dt_method  
        character(len=256)  :: solver_advec 
        integer             :: nx, ny 
        real(wp)          :: dx, dy  
        integer             :: nz_aa     ! Number of vertical points in ice (layer centers, plus base and surface)
        integer             :: nz_ac     ! Number of vertical points in ice (layer boundaries)
        integer             :: nzr       ! Number of vertical points in bedrock 
        real(wp)          :: gamma  
        logical             :: use_strain_sia 
        integer             :: n_sm_qstrn    ! Standard deviation (in points) for Gaussian smoothing of strain heating
        integer             :: n_sm_qb       ! Standard deviation (in points) for Gaussian smoothing of basal heating
        logical             :: use_const_cp 
        real(wp)          :: const_cp 
        logical             :: use_const_kt 
        real(wp)          :: const_kt 
        real(wp)          :: enth_cr  
        real(wp)          :: omega_max 
        real(wp)          :: till_rate 
        real(wp)          :: H_w_max 
        
        real(wp), allocatable :: zeta_aa(:)   ! Layer centers (aa-nodes), plus base and surface: nz_aa points 
        real(wp), allocatable :: zeta_ac(:)   ! Layer borders (ac-nodes), plus base and surface: nz_ac == nz_aa-1 points

        real(wp), allocatable :: dzeta_a(:)
        real(wp), allocatable :: dzeta_b(:)
        
        real(dp)   :: time
        real(wp) :: dt_zeta, dt_beta(2)

        real(wp) :: speed 

    end type

    ! ytherm state variables
    type ytherm_state_class
        real(wp), allocatable :: enth(:,:,:)      ! [J m-3] Ice enthalpy 
        real(wp), allocatable :: T_ice(:,:,:)     ! [K]     Ice temp. 
        real(wp), allocatable :: omega(:,:,:)     ! [--]    Ice water content
        real(wp), allocatable :: T_pmp(:,:,:)     ! Pressure-corrected melting point
        
        real(wp), allocatable :: f_pmp(:,:)       ! fraction of cell at pressure melting point
        real(wp), allocatable :: bmb_grnd(:,:)    ! Grounded basal mass balance 
        real(wp), allocatable :: Q_strn(:,:,:)    ! Internal heat production 
        real(wp), allocatable :: Q_b(:,:)         ! Basal friction heat production
        real(wp), allocatable :: Q_ice_b(:,:)     ! Basal ice heat flux 
        real(wp), allocatable :: T_prime_b(:,:)   ! Homologous temperature at the base 
        real(wp), allocatable :: H_w(:,:)         ! [m] Basal water layer thickness 
        real(wp), allocatable :: dHwdt(:,:)       ! [m/a] Basal water layer thickness rate of change
        
        real(wp), allocatable :: cp(:,:,:)        ! Specific heat capacity  
        real(wp), allocatable :: kt(:,:,:)        ! Heat conductivity  
        real(wp), allocatable :: H_cts(:,:)       ! Height of the cts
        
        real(wp), allocatable :: advecxy(:,:,:)

    end type

    ! ytherm class
    type ytherm_class

        type(ytherm_param_class)   :: par        ! physical parameters
        type(ytherm_state_class)   :: now

    end type

    ! =========================================================================
    !
    ! YELMO objects: ybound
    !
    ! =========================================================================
    
    ! ybnd variables (intent IN)
    type ybound_class

        ! Region constants
        real(wp) :: index_north = 1.0   ! Northern Hemisphere region number
        real(wp) :: index_south = 2.0   ! Antarctica region number
        real(wp) :: index_grl   = 1.3   ! Greenland region number

        ! Variables that save the current boundary conditions
        real(wp), allocatable :: z_bed(:,:)
        real(wp), allocatable :: z_bed_sd(:,:) 
        real(wp), allocatable :: z_sl(:,:)
        real(wp), allocatable :: H_sed(:,:)
        real(wp), allocatable :: smb(:,:)
        real(wp), allocatable :: T_srf(:,:)
        real(wp), allocatable :: bmb_shlf(:,:)
        real(wp), allocatable :: T_shlf(:,:)
        real(wp), allocatable :: Q_geo(:,:)

        real(wp), allocatable :: enh_srf(:,:)

        ! Useful masks
        real(wp), allocatable :: basins(:,:) 
        real(wp), allocatable :: basin_mask(:,:)
        real(wp), allocatable :: regions(:,:) 
        real(wp), allocatable :: region_mask(:,:) 

        logical,    allocatable :: ice_allowed(:,:)     ! Locations where ice thickness can be greater than zero 
        logical,    allocatable :: calv_mask(:,:)       ! for calv_method="kill-loc", where calv_mask==False, calv.
        
        real(wp), allocatable :: H_ice_ref(:,:)       ! Reference ice thickness, may be used for relaxation routines
        real(wp), allocatable :: z_bed_ref(:,:)       ! Reference bedrock elevation, may be used for relaxation routines

        ! Other external variables that can be useful, ie maybe with tracers
        ! to do 

    end type

    ! =========================================================================
    !
    ! YELMO objects: ydata
    !
    ! =========================================================================
    
    type ydata_param_class 
        logical             :: pd_topo_load 
        character(len=1028) :: pd_topo_path 
        character(len=56)   :: pd_topo_names(3)
        logical             :: pd_tsrf_load  
        character(len=1028) :: pd_tsrf_path 
        character(len=56)   :: pd_tsrf_name
        logical             :: pd_tsrf_monthly
        logical             :: pd_smb_load 
        character(len=1028) :: pd_smb_path 
        character(len=56)   :: pd_smb_name
        logical             :: pd_smb_monthly 
        logical             :: pd_vel_load  
        character(len=1028) :: pd_vel_path 
        character(len=56)   :: pd_vel_names(2) 
        logical             :: pd_age_load 
        character(len=1028) :: pd_age_path 
        character(len=56)   :: pd_age_names(2) 
        integer             :: pd_age_n_iso 

        character(len=56)   :: domain 
    end type 

    type ydata_pd_class   ! pd = present-day
        ! Variables that contain observations / reconstructions for comparison/inversion
        real(wp), allocatable :: H_ice(:,:), z_srf(:,:), z_bed(:,:), H_grnd(:,:)
        real(wp), allocatable :: ux_s(:,:), uy_s(:,:), uxy_s(:,:) 
        real(wp), allocatable :: T_srf(:,:), smb(:,:)
        real(wp), allocatable :: depth_iso(:,:,:)  
        
        ! Comparison metrics 
        real(wp), allocatable :: err_H_ice(:,:), err_z_srf(:,:), err_z_bed(:,:)
        real(wp), allocatable :: err_uxy_s(:,:)
        real(wp), allocatable :: err_depth_iso(:,:,:) 

        ! Axis 
        real(wp), allocatable :: age_iso(:) 

        real(wp) :: rmse_H 
        real(wp) :: rmse_zsrf
        real(wp) :: rmse_uxy 
        real(wp) :: rmse_loguxy 
        real(wp), allocatable :: rmse_iso(:) 
    
    end type

    type ydata_class 
        type(ydata_param_class) :: par 
        type(ydata_pd_class)    :: pd 
    end type 

    ! =========================================================================
    !
    ! YELMO objects: yregions
    !
    ! =========================================================================
    
    ! yregions variables
    type yregions_class
        ! Individual values of interest to output from a Yelmo domain 

        ! ===== Total ice variables =====
        real(wp) :: H_ice, z_srf, dHicedt, H_ice_max, dzsrfdt
        real(wp) :: V_ice, A_ice, dVicedt, fwf, V_sl 
        real(wp) :: uxy_bar, uxy_s, uxy_b, z_bed, smb, T_srf, bmb

        ! ===== Grounded ice variables =====
        real(wp) :: H_ice_g, z_srf_g, V_ice_g, A_ice_g, uxy_bar_g, uxy_s_g, uxy_b_g
        real(wp) :: f_pmp, H_w, bmb_g 

        ! ===== Floating ice variables =====
        real(wp) :: H_ice_f, V_ice_f, A_ice_f, uxy_bar_f, uxy_s_f, uxy_b_f, z_sl, bmb_shlf, T_shlf
        
    end type

    ! =========================================================================
    !
    ! YELMO objects: ygrid 
    !
    ! =========================================================================
    
    type ygrid_class 

        ! Grid name 
        character(len=256) :: name 
        
        ! Parameters
        integer    :: nx, ny, npts
        real(wp) :: dx, dy

        ! Projection parameters (optional)
        character(len=256) :: mtype 
        real(wp) :: lambda
        real(wp) :: phi
        real(wp) :: alpha
        real(wp) :: scale
        real(wp) :: x_e
        real(wp) :: y_n
        logical    :: is_projection 

        ! Axes
        real(wp), allocatable :: xc(:)    
        real(wp), allocatable :: yc(:) 

        ! Grid arrays 
        real(wp), allocatable :: x(:,:)
        real(wp), allocatable :: y(:,:)
        real(wp), allocatable :: lon(:,:)
        real(wp), allocatable :: lat(:,:)
        real(wp), allocatable :: area(:,:)
        
    end type 

    ! =========================================================================
    !
    ! YELMO objects: yelmo 
    !
    ! =========================================================================
    
    ! Define all parameters needed to represent a given domain
    type yelmo_param_class

        ! Domain and experiment definition
        character (len=256) :: domain
        character (len=256) :: grid_name
        character (len=512) :: grid_path
        character (len=256) :: experiment
        character (len=512) :: restart

        ! Data logging 
        logical             :: log_timestep 

        ! Vertical dimension definition
        character (len=56)  :: zeta_scale 
        real(wp)          :: zeta_exp 
        integer             :: nz_ac
        integer             :: nz_aa 

        ! Yelmo timesteps 
        integer             :: dt_method 
        real(wp)          :: dt_min
        real(wp)          :: cfl_max 
        real(wp)          :: cfl_diff_max 
        character (len=56)  :: pc_method
        character (len=56)  :: pc_controller
        logical             :: pc_filter_vel 
        logical             :: pc_use_H_pred 
        integer             :: pc_n_redo 
        real(wp)          :: pc_tol 
        real(wp)          :: pc_eps  

        ! Sigma coordinates (internal parameter)
        real(wp), allocatable :: zeta_aa(:)   ! Layer centers (aa-nodes), plus base and surface: nz_aa points 
        real(wp), allocatable :: zeta_ac(:)   ! Layer borders (ac-nodes), plus base and surface: nz_ac == nz_aa-1 points
        
        ! Other internal parameters
        logical :: use_restart 
        
    end type

    type ytime_class 

        ! Time step parameters for predictor-corrector (PC) method (Cheng et al, 2017)
        real(wp) :: pc_dt(3)
        real(wp) :: pc_eta(3)
        real(wp), allocatable :: pc_tau(:,:)
        real(wp), allocatable :: pc_tau_masked(:,:)
        
        real(wp), allocatable :: dt_adv(:,:) 
        real(wp), allocatable :: dt_diff(:,:) 
        real(wp), allocatable :: dt_adv3D(:,:,:)
        
        ! Timing information
        real(wp) :: model_speed 
        real(wp) :: model_speeds(100)         ! Eg, 100 timesteps for running mean 
        real(wp) :: dt_avg 
        real(wp) :: dts(100)                  ! Eg, 100 timesteps for running mean
        real(wp) :: eta_avg 
        real(wp) :: etas(100)                 ! Eg, 100 timesteps for running mean
        real(wp) :: ssa_iter_avg 
        real(wp) :: ssa_iters(100)            ! Eg, 100 timesteps for running mean
        
        ! Truncation error information over several timesteps
        real(wp), allocatable :: pc_taus(:,:,:)
        real(wp), allocatable :: pc_tau_max(:,:)

        character(len=512)   :: log_timestep_file 

    end type 

    ! Define the overall yelmo_class, which is a container for
    ! all information needed to model a given domain (eg, Greenland, Antarctica, NH)
    type yelmo_class
        type(yelmo_param_class) :: par      ! General domain parameters
        type(ygrid_class)       :: grd      ! Grid definition
        type(ytime_class)       :: time     ! Timestep and timing variables
        type(ytopo_class)       :: tpo      ! Topography variables
        type(ydyn_class)        :: dyn      ! Dynamics variables
        type(ymat_class)        :: mat      ! Material variables
        type(ytherm_class)      :: thrm     ! Thermodynamics variables
        type(ybound_class)      :: bnd      ! Boundary variables to drive model
        type(ydata_class)       :: dta      ! Data variables for comparison
        type(yregions_class)    :: reg      ! Regionally aggregated variables for whole domain 
    end type

    public   ! All yelmo defs are public

contains 

    function yelmo_get_precision() result(yelmo_prec)

        implicit none 

        integer :: yelmo_prec 

        yelmo_prec = kind(wp)

        return 

    end function yelmo_get_precision

        
    subroutine yelmo_parse_path(path,domain,grid_name)

        implicit none 

        character(len=*), intent(INOUT) :: path 
        character(len=*), intent(IN)    :: domain, grid_name 

        call nml_replace(path,"{domain}",   trim(domain))
        call nml_replace(path,"{grid_name}",trim(grid_name))
        
        return 

    end subroutine yelmo_parse_path

    subroutine yelmo_global_init(filename)

        !$ use omp_lib 

        implicit none 

        character(len=*), intent(IN)  :: filename
        
        ! Local variables
        logical :: init_pars
        integer :: n_threads 
        character(len=10) :: n_threads_str 

        init_pars = .TRUE. 
        
        ! Check openmp status - set global variable to use as a switch 
        yelmo_use_omp = .FALSE. 
        !$ yelmo_use_omp = .TRUE.

        ! Output some information about openmp status 
        if (yelmo_use_omp) then 
            
            n_threads = 1
            !$ n_threads = omp_get_max_threads() 

            write(n_threads_str,"(i10)") n_threads 
            n_threads_str = adjustl(n_threads_str)

            write(*,*) "yelmo_global_init:: openmp is active, Yelmo will run on "//trim(n_threads_str)//" thread(s)."
            
        else 
            
            n_threads = 1
            write(*,*) "yelmo_global_init:: openmp is not active, Yelmo will run on 1 thread."

        end if 

        ! Load parameter values 

        call nml_read(filename,"yelmo_config","yelmo_log",yelmo_log,init=init_pars)

        call nml_read(filename,"yelmo_constants","sec_year",    sec_year,   init=init_pars)
        call nml_read(filename,"yelmo_constants","g",           g,          init=init_pars)
        call nml_read(filename,"yelmo_constants","T0",          T0,         init=init_pars)
        
        call nml_read(filename,"yelmo_constants","rho_ice",     rho_ice,    init=init_pars)
        call nml_read(filename,"yelmo_constants","rho_w",       rho_w,      init=init_pars)
        call nml_read(filename,"yelmo_constants","rho_sw",      rho_sw,     init=init_pars)
        call nml_read(filename,"yelmo_constants","rho_a",       rho_a,      init=init_pars)
        call nml_read(filename,"yelmo_constants","rho_m",       rho_m,      init=init_pars)
        call nml_read(filename,"yelmo_constants","L_ice",       L_ice,      init=init_pars)
        call nml_read(filename,"yelmo_constants","T_pmp_beta",  T_pmp_beta, init=init_pars)
        
        if (yelmo_log) then 
            write(*,*) "yelmo:: configuration:"
            write(*,*) "    yelmo_log          = ", yelmo_log
            
            write(*,*) "yelmo:: loaded global constants:"
            write(*,*) "    sec_year   = ", sec_year 
            write(*,*) "    g          = ", g 
            write(*,*) "    T0         = ", T0 
            write(*,*) "    rho_ice    = ", rho_ice 
            write(*,*) "    rho_w      = ", rho_w 
            write(*,*) "    rho_sw     = ", rho_sw 
            write(*,*) "    rho_a      = ", rho_a 
            write(*,*) "    rho_m      = ", rho_m 
            write(*,*) "    L_ice      = ", L_ice 
            write(*,*) "    T_pmp_beta = ", T_pmp_beta 
            
        end if 

        ! Define conversion factors too

        conv_we_ie  = rho_w/rho_ice
        conv_mmdwe_maie = 1e-3*365*conv_we_ie
        conv_mmawe_maie = 1e-3*conv_we_ie
        
        return

    end subroutine yelmo_global_init
    
    subroutine yelmo_load_command_line_args(path_par)
        ! Load the parameter filename from the command line 
        ! call eg: ./yelmo_test.x yelmo_Greenland.nml 

        implicit none 

        character(len=*), intent(OUT) :: path_par 

        ! Local variables 
        integer :: narg 

        narg = command_argument_count()

        if (narg .ne. 1) then 
            write(*,*) "yelmo_load_command_line_args:: Error: The following &
            &argument must be provided: path_par"
            stop 
        end if 

        call get_command_argument(1,path_par)

        return 

    end subroutine yelmo_load_command_line_args

    subroutine yelmo_cpu_time(time,time0,dtime)
        ! Calculate time intervals using system_clock.

        ! Note: for mulithreading, cpu_time() won't work properly.
        ! Instead, system_clock() should be used as it is here, 
        ! unless use_cpu_time=.TRUE. 

        !$ use omp_lib

        implicit none 

        real(8), intent(OUT) :: time 
        real(8), intent(IN),  optional :: time0 
        real(8), intent(OUT), optional :: dtime 

        ! Local variables
        logical    :: using_omp 
        integer(4) :: clock 
        integer(4) :: clock_rate
        integer(4) :: clock_max 
        real(8)    :: wtime 

        ! Check openmp status - do not use global switch, since it may not have been initialized yet
        using_omp = .FALSE. 
        !$ using_omp = .TRUE.

        if (using_omp) then 
            ! --------------------------------------
            ! omp_get_wtime must be used for multithread openmp execution to get timing on master thread 
            ! The following lines will overwrite time with the result from omp_get_wtime on the master thread 

            !$ time = omp_get_wtime()

            ! --------------------------------------
            
        else 

            ! cpu_time can be used for serial execution to get timing on 1 processor
            call cpu_time(time)

        end if 

        if (present(dtime)) then 
            ! Calculate time interval 

            if (.not. present(time0)) then  
                write(*,*) "yelmo_cpu_time:: Error: time0 argument is missing, but necessary."
                stop
            end if 
            
            ! Calculate the difference between current time and time0 in [s]
            dtime = time - time0

            ! Limit dtime to non-zero number 
            if (dtime .eq. 0.0d0) then
                write(*,*) "yelmo_cpu_time:: Error: dtime cannot equal zero - check precision of timing variables, &
                            &which should be real(kind=8) to maintain precision."
                write(*,*) "clock", time, time0, dtime  
                stop  
            end if 

        end if 

        return 

    end subroutine yelmo_cpu_time 

    subroutine yelmo_calc_speed(speed,model_time0,model_time1,cpu_time0,cpu_time1)
        ! Calculate the model computational speed [model-kyr / hr]
        ! given model times for start and end of window in [yr]
        ! and cpu processing times for start and end of window in [sec]

        implicit none 

        real(wp), intent(OUT) :: speed            ! [kyr / hr]
        real(wp), intent(IN)  :: model_time0      ! [yr]
        real(wp), intent(IN)  :: model_time1      ! [yr]
        real(8),    intent(IN)  :: cpu_time0        ! [sec]
        real(8),    intent(IN)  :: cpu_time1        ! [sec]

        ! Local variables 
        real(wp) :: cpu_dtime 
        real(wp) :: model_dtime 

        cpu_dtime   = (cpu_time1 - cpu_time0)/3600.d0       ! [sec] => [hr] 
        model_dtime = (model_time1 - model_time0)*1d-3      ! [yr] => [kyr] 

        if (cpu_dtime .gt. 0.0) then 
            speed = model_dtime / cpu_dtime                 ! [kyr / hr]
        else 
            speed = 0.0 
        end if 

        return 

    end subroutine yelmo_calc_speed 
    
    subroutine yelmo_calc_running_stats(val_out,vals,val_now,stat)

        implicit none 

        real(wp), intent(OUT)   :: val_out 
        real(wp), intent(INOUT) :: vals(:) 
        real(wp), intent(IN)    :: val_now 
        character(len=*), intent(IN) :: stat          ! 'mean', 'min', 'max', 'stdev'
        
        ! Local variables 
        integer    :: n 
        real(wp) :: val_mean 

        ! Shift rates vector to eliminate oldest entry, and add current entry in the first position
        vals    = cshift(vals,-1)
        vals(1) = val_now  

        ! Calculate running stats value 
        n    = count(vals .ne. 0.0_wp)
        if (n .gt. 0) then 

            select case(trim(stat))

                case("mean")
                    val_out = sum(vals,mask=vals .ne. 0.0_wp) / real(n,wp)
                case("min")
                    val_out = minval(vals,mask=vals .ne. 0.0_wp)
                case("max")
                    val_out = maxval(vals,mask=vals .ne. 0.0_wp)
                case("stdev")
                    val_mean = sum(vals,mask=vals .ne. 0.0_wp) / real(n,wp)
                    val_out  = sqrt(sum((vals-val_mean)**2,mask=vals .ne. 0.0_wp) / real(n,wp))

                case DEFAULT 
                    write(*,*) "yelmo_calc_running_stats:: Error: stat not found."
                    write(*,*) "stat = ", trim(stat) 
                    stop 

            end select 
                   
        else

            val_out = 0.0_wp
        
        end if 

        return 

    end subroutine yelmo_calc_running_stats

    subroutine yelmo_calc_running_stats_2D(val_out,vals,val_now,stat)

        implicit none 

        real(wp), intent(OUT)   :: val_out(:,:) 
        real(wp), intent(INOUT) :: vals(:,:,:) 
        real(wp), intent(IN)    :: val_now(:,:) 
        character(len=*), intent(IN) :: stat            ! 'mean', 'min', 'max', 'stdev'
        
        ! Local variables
        integer :: i, j, nx, ny  
        integer :: n, k, n_now  

        nx = size(val_out,1)
        ny = size(val_out,2) 

        do j = 1, ny 
        do i = 1, nx 

            call yelmo_calc_running_stats(val_out(i,j),vals(i,j,:),val_now(i,j),stat)

        end do 
        end do  

        return 

    end subroutine yelmo_calc_running_stats_2D

end module yelmo_defs

