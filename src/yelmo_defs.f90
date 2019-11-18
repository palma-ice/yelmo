
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

    ! Missing value and aliases
    real(prec), parameter :: MISSING_VALUE_DEFAULT = real(-9999.0,prec)
    real(prec), parameter :: MISSING_VALUE = MISSING_VALUE_DEFAULT
    real(prec), parameter :: MV = MISSING_VALUE_DEFAULT
    
    ! Error distance (very large), error index, and smallest number epsilon 
    real(prec), parameter :: ERR_DIST = real(1E8,prec) 
    integer,    parameter :: ERR_IND  = -1 
    real(prec), parameter :: tol_underflow = real(1E-15,prec)

    ! Mathematical constants
    real(prec), parameter :: pi  = real(2._dp*acos(0.0_dp),prec)
    real(prec), parameter :: degrees_to_radians = real(pi / 180._dp,prec)  ! Conversion factor between radians and degrees
    real(prec), parameter :: radians_to_degrees = real(180._dp / pi,prec)  ! Conversion factor between degrees and radians
    
    ! The constants below should be loaded using the global subroutine
    ! defined below `yelmo_constants_load`.
    ! Note: The key limitation imposed by defining the parameters defined 
    ! globally is that these constants must be the same for all domains 
    ! being run in the same program. 

    ! Yelmo configuration options 
    logical :: yelmo_log

    ! Physical constants 
    real(prec) :: sec_year       ! [s] seconds per year 
    real(prec) :: g              ! [m s-2] Gravitational accel.  
    real(prec) :: T0             ! [K] Reference freezing temperature  
    real(prec) :: rho_ice        ! [kg m-3] Density ice           
    real(prec) :: rho_w          ! [kg m-3] Density water          
    real(prec) :: rho_sw         ! [kg m-3] Density seawater      
    real(prec) :: rho_a          ! [kg m-3] Density asthenosphere  
    real(prec) :: rho_m          ! [kg m-3] Density mantle (lith) 
    real(prec) :: L_ice          ! [J kg-1] Latent heat           
    real(prec) :: T_pmp_beta     ! [K Pa-1] Melt point pressure slope

    ! Internal parameters 
    real(prec) :: conv_we_ie        ! Conversion water equiv. => m/a ice equiv. 
    real(prec) :: conv_mmdwe_maie   ! Conversion mm/d water equiv. => m/a ice equiv.
    real(prec) :: conv_mmawe_maie   ! Conversion mm/a water equiv. => m/a ice equiv. 
    
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
        real(prec)         :: topo_rel_tau 
        real(prec)         :: calv_H_lim
        real(prec)         :: calv_tau  
        real(prec)         :: H_min_grnd
        real(prec)         :: H_min_flt 
        real(prec)         :: sd_min 
        real(prec)         :: sd_max 
        real(prec)         :: calv_max  
        real(prec)         :: grad_lim 
        integer            :: gl_sep 
        integer            :: gl_sep_nx 
        logical            :: diffuse_bmb_shlf 

        ! Internal parameters 
        real(prec)         :: time 
        integer            :: nx, ny
        real(prec)         :: dx, dy
        character(len=256) :: boundaries 
        
    end type

    ! ytopo state variables
    type ytopo_state_class
        ! Model variables that the define the state of the domain 

        real(prec), allocatable :: H_ice(:,:)      ! Ice thickness [m] 
        real(prec), allocatable :: z_srf(:,:)      ! Surface elevation [m]
        real(prec), allocatable :: dzsrfdt(:,:)    ! Surface elevation rate of change [m/a] 
        real(prec), allocatable :: dHicedt(:,:)    ! Ice thickness rate of change [m/a] 
        real(prec), allocatable :: bmb(:,:)        ! Combined field of bmb_grnd and bmb_shlf 
        real(prec), allocatable :: mb_applied(:,:) ! Actual mass balance applied [m/a], for mass balance accounting
        real(prec), allocatable :: calv_grnd(:,:)  ! Grounded calving [m/a]
        real(prec), allocatable :: calv(:,:)       ! Calving [m/a]
        real(prec), allocatable :: calv_mean(:,:)  ! Calving [m/a]
        real(prec), allocatable :: calv_times(:)   ! Calving

        real(prec), allocatable :: H_margin(:,:)   ! [m] Margin ice thickness in partially filled cells 

        real(prec), allocatable :: dzsdx(:,:)      ! Surface elevation slope [m m-1], Ac x nodes
        real(prec), allocatable :: dzsdy(:,:)      ! Surface elevation slope [m m-1], Ac y nodes
        real(prec), allocatable :: dHicedx(:,:)    ! Ice thickness gradient slope [m m-1], Ac x nodes
        real(prec), allocatable :: dHicedy(:,:)    ! Ice thickness gradient slope [m m-1], Ac y nodes
        
        real(prec), allocatable :: H_grnd(:,:)       ! Ice thickness overburden [m]
        
        ! Masks 
        real(prec), allocatable :: f_grnd(:,:)       ! Grounded fraction (grounding line fraction between 0 and 1)
        real(prec), allocatable :: f_grnd_acx(:,:)   ! Grounded fraction (acx nodes)
        real(prec), allocatable :: f_grnd_acy(:,:)   ! Grounded fraction (acy nodes)
        real(prec), allocatable :: f_ice(:,:)        ! Ice-covered fraction 

        real(prec), allocatable :: dist_margin(:,:)  ! Distance to nearest margin point 
        real(prec), allocatable :: dist_grline(:,:)  ! Distance to nearest grounding-line point 
        
        ! Additional masks 
        integer,    allocatable :: mask_bed(:,:)    ! Multi-valued bed mask
        logical,    allocatable :: is_float(:,:)    ! Fully floating grid points
        logical,    allocatable :: is_grline(:,:)   ! Grounding line points
        logical,    allocatable :: is_grz(:,:)      ! Grounding line plus grounded neighbors
        
    end type

    ! ytopo class
    type ytopo_class

        type(ytopo_param_class) :: par        ! Parameters
        type(ytopo_state_class) :: now        ! Variables

    end type

    ! =========================================================================
    !
    ! YELMO objects: ydyn 
    !
    ! =========================================================================
    
    ! ydyn parameters
    type ydyn_param_class

        character(len=256) :: solver 
        character(len=256) :: sia_solver 
        integer    :: mix_method            ! Method for mixing sia and ssa velocity solutions
        logical    :: calc_diffusivity      ! Calculate diagnostic diffusivity field
        integer    :: beta_method
        real(prec) :: beta_const
        real(prec) :: beta_q                ! Friction law exponent
        real(prec) :: beta_u0               ! [m/a] Friction law velocity threshold 
        integer    :: beta_gl_scale         ! Beta grounding-line scaling method (beta => 0 at gl?)
        integer    :: beta_gl_sep           ! Beta grounding-line sub-element (subgrid) parameterization
        integer    :: beta_gl_stag          ! Beta grounding-line staggering method 
        real(prec) :: beta_gl_f             ! Fraction of beta at gl 
        integer    :: taud_gl_method        ! Driving stress grounding line treatment 
        real(prec) :: H_grnd_lim 
        real(prec) :: H_sed_sat
        integer    :: cb_method
        logical    :: cb_with_pmp           ! Scale friction coefficient between frozen and streaming values?
        logical    :: cb_margin_pmp         ! Ensure margin and grline are considered streaming?
        character(len=256) :: cb_scale      ! Method for scaling cb with elevation
        real(prec) :: cb_z0  
        real(prec) :: cb_z1
        real(prec) :: cb_min      
        real(prec) :: cf_frozen
        real(prec) :: cf_stream
        integer    :: n_sm_beta 
        real(prec) :: beta_min              ! Minimum allowed value of beta
        real(prec) :: ssa_beta_max          ! Maximum value of beta for which ssa should be calculated
        real(prec) :: ssa_vel_max
        integer    :: ssa_iter_max 
        real(prec) :: ssa_iter_rel 
        real(prec) :: ssa_iter_conv 
        real(prec) :: cb_sia
        
        integer    :: neff_method
        real(prec) :: neff_const
        real(prec) :: neff_p 
        logical    :: neff_set_water 
        real(prec) :: neff_w_max
        real(prec) :: neff_N0
        real(prec) :: neff_delta 
        real(prec) :: neff_e0 
        real(prec) :: neff_Cc 

        real(prec) :: till_phi_const 
        real(prec) :: till_phi_min 
        real(prec) :: till_phi_max 
        real(prec) :: till_phi_zmin 
        real(prec) :: till_phi_zmax 
        
        ! Internal parameters 
        character(len=256) :: boundaries 
        logical    :: use_ssa                    ! Should ssa be used? 
        logical    :: use_bmb                    ! Set to match `use_bmb` in ytopo_param_class 
        integer    :: nx, ny, nz_aa, nz_ac 
        real(prec) :: dx, dy
        real(prec), allocatable :: zeta_aa(:)   ! Layer centers (aa-nodes), plus base and surface: nz_aa points 
        real(prec), allocatable :: zeta_ac(:)   ! Layer borders (ac-nodes), plus base and surface: nz_ac == nz_aa-1 points
        real(prec) :: time

    end type

    ! ydyn state variables
    type ydyn_state_class
        ! Model variables that the define the state of the domain 

        real(prec), allocatable :: ux(:,:,:) 
        real(prec), allocatable :: uy(:,:,:) 
        real(prec), allocatable :: uxy(:,:,:)
        real(prec), allocatable :: uz(:,:,:)  

        real(prec), allocatable :: ux_bar(:,:) 
        real(prec), allocatable :: uy_bar(:,:)
        real(prec), allocatable :: uxy_bar(:,:)

        real(prec), allocatable :: ux_b(:,:) 
        real(prec), allocatable :: uy_b(:,:)
        real(prec), allocatable :: uxy_b(:,:)

        ! Surface velocity: eventually these could be pointers since it is simply
        ! the top layer in ux(:,:,:), etc. and only used, not calculated.
        real(prec), allocatable :: ux_s(:,:) 
        real(prec), allocatable :: uy_s(:,:)
        real(prec), allocatable :: uxy_s(:,:)
        
        real(prec), allocatable :: ux_i(:,:,:) 
        real(prec), allocatable :: uy_i(:,:,:)
        real(prec), allocatable :: ux_i_bar(:,:) 
        real(prec), allocatable :: uy_i_bar(:,:)
        real(prec), allocatable :: uxy_i_bar(:,:) 
        
        real(prec), allocatable :: dd_ab(:,:,:)  
        real(prec), allocatable :: dd_ab_bar(:,:)  
        
        real(prec), allocatable :: sigma_horiz_sq(:,:)
        real(prec), allocatable :: lhs_x(:,:) 
        real(prec), allocatable :: lhs_y(:,:) 
        real(prec), allocatable :: lhs_xy(:,:) 

        real(prec), allocatable :: duxdz(:,:,:) 
        real(prec), allocatable :: duydz(:,:,:)
        real(prec), allocatable :: duxdz_bar(:,:) 
        real(prec), allocatable :: duydz_bar(:,:)

        real(prec), allocatable :: taud_acx(:,:) 
        real(prec), allocatable :: taud_acy(:,:) 
        real(prec), allocatable :: taud(:,:) 
        
        real(prec), allocatable :: taub_acx(:,:) 
        real(prec), allocatable :: taub_acy(:,:) 
        real(prec), allocatable :: taub(:,:)
        
        real(prec), allocatable :: qq_gl_acx(:,:) 
        real(prec), allocatable :: qq_gl_acy(:,:) 
        
        real(prec), allocatable :: qq_acx(:,:) 
        real(prec), allocatable :: qq_acy(:,:) 
        real(prec), allocatable :: qq(:,:)
        
        real(prec), allocatable :: visc_eff(:,:)

        real(prec), allocatable :: N_eff(:,:)       ! Effective pressure
        real(prec), allocatable :: cf_ref(:,:)
        real(prec), allocatable :: c_bed(:,:)  
        real(prec), allocatable :: beta_acx(:,:) 
        real(prec), allocatable :: beta_acy(:,:) 
        real(prec), allocatable :: beta(:,:) 
        
        real(prec), allocatable :: f_vbvs(:,:) 

        integer,    allocatable :: ssa_mask_acx(:,:) 
        integer,    allocatable :: ssa_mask_acy(:,:) 

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
        real(prec), allocatable :: dxx(:,:) 
        real(prec), allocatable :: dyy(:,:) 
        real(prec), allocatable :: dxy(:,:) 
        real(prec), allocatable :: de(:,:) 
    end type 

    type strain_3D_class 
        real(prec), allocatable :: dxx(:,:,:) 
        real(prec), allocatable :: dyy(:,:,:) 
        real(prec), allocatable :: dzz(:,:,:)
        real(prec), allocatable :: dxy(:,:,:) 
        real(prec), allocatable :: dxz(:,:,:) 
        real(prec), allocatable :: dyz(:,:,:) 
        real(prec), allocatable :: de(:,:,:) 
        real(prec), allocatable :: f_shear(:,:,:) 
    end type 
    
    type ymat_param_class
        
        character(len=56) :: flow_law
        integer    :: rf_method 
        real(prec) :: rf_const
        logical    :: rf_use_eismint2
        logical    :: rf_with_water 
        real(prec) :: n_glen                     ! Flow law exponent (n_glen=3)
        real(prec) :: visc_min  
        logical    :: use_2D_enh
        real(prec) :: enh_shear
        real(prec) :: enh_stream
        real(prec) :: enh_shlf
        
        
        character(len=56) :: age_method  
        real(prec)        :: age_impl_kappa

        ! Internal parameters
        real(prec) :: time 
        logical    :: calc_age
        real(prec) :: dx, dy  
        integer    :: nx, ny, nz_aa, nz_ac  

        real(prec), allocatable :: zeta_aa(:)   ! Layer centers (aa-nodes), plus base and surface: nz_aa points 
        real(prec), allocatable :: zeta_ac(:)   ! Layer borders (ac-nodes), plus base and surface: nz_ac == nz_aa-1 points
        
    end type 

    type ymat_state_class 

        type(strain_2D_class)   :: strn2D
        type(strain_3D_class)   :: strn 

        real(prec), allocatable :: enh(:,:,:)
        real(prec), allocatable :: enh_bar(:,:)
        real(prec), allocatable :: ATT(:,:,:) 
        real(prec), allocatable :: ATT_bar(:,:)
        real(prec), allocatable :: visc(:,:,:) 
        real(prec), allocatable :: visc_int(:,:) 

        real(prec), allocatable :: f_shear_bar(:,:) 
        
        real(prec), allocatable :: dep_time(:,:,:)    ! Ice deposition time (for online age tracing)

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
        character(len=256)  :: solver_advec 
        integer             :: nx, ny 
        real(prec)          :: dx, dy  
        integer             :: nz_aa     ! Number of vertical points in ice (layer centers, plus base and surface)
        integer             :: nz_ac     ! Number of vertical points in ice (layer boundaries)
        integer             :: nzr       ! Number of vertical points in bedrock 
        real(prec)          :: gamma  
        logical             :: use_strain_sia 
        integer             :: n_sm_qstrn    ! Standard deviation (in points) for Gaussian smoothing of strain heating
        integer             :: n_sm_qb       ! Standard deviation (in points) for Gaussian smoothing of basal heating
        logical             :: use_const_cp 
        real(prec)          :: const_cp 
        logical             :: use_const_kt 
        real(prec)          :: const_kt 
        real(prec)          :: enth_cr  
        real(prec)          :: omega_max 
        
        real(prec), allocatable :: zeta_aa(:)   ! Layer centers (aa-nodes), plus base and surface: nz_aa points 
        real(prec), allocatable :: zeta_ac(:)   ! Layer borders (ac-nodes), plus base and surface: nz_ac == nz_aa-1 points

        real(prec), allocatable :: dzeta_a(:)
        real(prec), allocatable :: dzeta_b(:)
        
        real(prec) :: time

    end type

    ! ytherm state variables
    type ytherm_state_class
        real(prec), allocatable :: enth(:,:,:)      ! [J m-3] Ice enthalpy 
        real(prec), allocatable :: T_ice(:,:,:)     ! [K]     Ice temp. 
        real(prec), allocatable :: omega(:,:,:)     ! [--]    Ice water content
        real(prec), allocatable :: T_pmp(:,:,:)     ! Pressure-corrected melting point
        real(prec), allocatable :: f_pmp(:,:)       ! fraction of cell at pressure melting point
        real(prec), allocatable :: bmb_grnd(:,:)    ! Grounded basal mass balance 
        real(prec), allocatable :: Q_strn(:,:,:)    ! Internal heat production 
        real(prec), allocatable :: Q_b(:,:)         ! Basal friction heat production
        real(prec), allocatable :: Q_ice_b(:,:)     ! Basal ice heat flux 
        real(prec), allocatable :: cp(:,:,:)        ! Specific heat capacity  
        real(prec), allocatable :: kt(:,:,:)        ! Heat conductivity  
        real(prec), allocatable :: H_cts(:,:)       ! Height of the cts
        real(prec), allocatable :: T_prime_b(:,:)   ! Homologous temperature at the base 
        real(prec), allocatable :: H_w(:,:)         ! [m] Basal water layer thickness 
        real(prec), allocatable :: dHwdt(:,:)       ! [m/a] Basal water layer thickness rate of change
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
        real(prec) :: index_north = 1.0   ! Northern Hemisphere region number
        real(prec) :: index_south = 2.0   ! Antarctica region number
        real(prec) :: index_grl   = 1.3   ! Greenland region number

        ! Variables that save the current boundary conditions
        real(prec), allocatable :: z_bed(:,:)
        real(prec), allocatable :: z_bed_sd(:,:) 
        real(prec), allocatable :: z_sl(:,:)
        real(prec), allocatable :: H_sed(:,:)
        real(prec), allocatable :: smb(:,:)
        real(prec), allocatable :: T_srf(:,:)
        real(prec), allocatable :: bmb_shlf(:,:)
        real(prec), allocatable :: T_shlf(:,:)
        real(prec), allocatable :: Q_geo(:,:)

        ! Useful masks
        real(prec), allocatable :: basins(:,:) 
        real(prec), allocatable :: basin_mask(:,:)
        real(prec), allocatable :: regions(:,:) 
        real(prec), allocatable :: region_mask(:,:) 

        logical,    allocatable :: ice_allowed(:,:)     ! Locations where ice thickness can be greater than zero 
        logical,    allocatable :: calv_mask(:,:)       ! for calv_method="kill-loc", where calv_mask==False, calv.
        real(prec), allocatable :: H_ice_ref(:,:)       ! Reference ice thickness, may be used for relaxation routines

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
        
        character(len=56)   :: domain 
    end type 

    type ydata_pd_class   ! pd = present-day
        ! Variables that contain observations / reconstructions for comparison/inversion
        real(prec), allocatable :: H_ice(:,:), z_srf(:,:), z_bed(:,:), H_grnd(:,:)
        real(prec), allocatable :: ux_s(:,:), uy_s(:,:), uxy_s(:,:) 
        real(prec), allocatable :: T_srf(:,:), smb(:,:) 
        ! Comparison metrics 
        real(prec), allocatable :: err_H_ice(:,:), err_z_srf(:,:), err_z_bed(:,:)
        real(prec), allocatable :: err_uxy_s(:,:)
        
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
        real(prec) :: H_ice, z_srf, dHicedt, H_ice_max, dzsrfdt
        real(prec) :: V_ice, A_ice, dVicedt, fwf
        real(prec) :: uxy_bar, uxy_s, uxy_b, z_bed, smb, T_srf, bmb

        ! ===== Grounded ice variables =====
        real(prec) :: H_ice_g, z_srf_g, V_ice_g, A_ice_g, uxy_bar_g, uxy_s_g, uxy_b_g
        real(prec) :: f_pmp, H_w, bmb_g 

        ! ===== Floating ice variables =====
        real(prec) :: H_ice_f, V_ice_f, A_ice_f, uxy_bar_f, uxy_s_f, uxy_b_f, z_sl, bmb_shlf, T_shlf
        
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
        real(prec) :: dx, dy

        ! Projection parameters (optional)
        character(len=256) :: mtype 
        real(prec) :: lambda
        real(prec) :: phi
        real(prec) :: alpha
        real(prec) :: scale
        real(prec) :: x_e
        real(prec) :: y_n
        logical    :: is_projection 

        ! Axes
        real(prec), allocatable :: xc(:)    
        real(prec), allocatable :: yc(:) 

        ! Grid arrays 
        real(prec), allocatable :: x(:,:)
        real(prec), allocatable :: y(:,:)
        real(prec), allocatable :: lon(:,:)
        real(prec), allocatable :: lat(:,:)
        real(prec), allocatable :: area(:,:)
        
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
        real(prec)          :: zeta_exp 
        integer             :: nz_ac
        integer             :: nz_aa 

        ! Yelmo timesteps 
        integer             :: dt_method 
        real(prec)          :: dt_min
        real(prec)          :: dt_ref 
        real(prec)          :: cfl_max 
        real(prec)          :: cfl_diff_max 
        real(prec)          :: pc_ebs 
        real(prec)          :: pc1_ebs 

        ! Sigma coordinates (internal parameter)
        real(prec), allocatable :: zeta_aa(:)   ! Layer centers (aa-nodes), plus base and surface: nz_aa points 
        real(prec), allocatable :: zeta_ac(:)   ! Layer borders (ac-nodes), plus base and surface: nz_ac == nz_aa-1 points
        
        ! Other internal parameters
        real(prec), allocatable :: dt_adv(:,:) 
        real(prec), allocatable :: dt_diff(:,:) 
        real(prec), allocatable :: dt_adv3D(:,:,:)
        logical :: use_restart 
        
        ! Time step parameters for predictor-corrector (PC) method (Cheng et al, 2017)
        real(prec) :: pc_dt 
        real(prec) :: pc_eta 
        real(prec) :: pc1_dt
        real(prec) :: pc1_eta 
        real(prec), allocatable :: pc_tau(:,:)
        real(prec), allocatable :: pc1_tau(:,:)

        ! Timing information 
        real(prec) :: model_speed 
        real(prec) :: model_speeds(10)    ! Use 10 timesteps for running mean  

        character(len=512)   :: log_timestep_file 

    end type

    ! Define the overall yelmo_class, which is a container for
    ! all information needed to model a given domain (eg, Greenland, Antarctica, NH)
    type yelmo_class
        type(yelmo_param_class) :: par      ! General domain parameters
        type(ygrid_class)       :: grd      ! Grid definition
        type(ytopo_class)       :: tpo      ! Topography variables
        type(ydyn_class)        :: dyn      ! Dynamics variables
        type(ymat_class)        :: mat      ! Material variables
        type(ytherm_class)      :: thrm     ! Thermodynamics variables
        type(ybound_class)      :: bnd      ! Boundary variables to drive model
        type(ydata_class)       :: dta      ! Data variables for comparison
        type(yregions_class)    :: reg      ! Regionally aggregated variables  
    end type

    public   ! All yelmo defs are public

contains 

    function yelmo_get_precision() result(yelmo_prec)

        implicit none 

        integer :: yelmo_prec 

        yelmo_prec = kind(prec)

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

        character(len=*), intent(IN)  :: filename
        
        ! Local variables
        logical :: init_pars 

        init_pars = .TRUE. 
        
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

     subroutine yelmo_calc_speed(rate,rates,model_time0,model_time1,cpu_time0)
        ! Calculate the model computational speed [model-kyr / hr]
        ! Note: uses a running mean of rates over the last X steps, in order
        ! to provide a smoother estimate of the rate 

        implicit none 

        real(prec), intent(OUT)   :: rate        ! [kyr / hr]
        real(prec), intent(INOUT) :: rates(:)    ! [kyr / hr]

        real(prec), intent(IN) :: model_time0    ! [yr]
        real(prec), intent(IN) :: model_time1    ! [yr]
        real(prec), intent(IN) :: cpu_time0      ! [sec]
        
        ! Local variables
        integer    :: ntot, n 
        real(4)    :: cpu_time1      ! [sec]
        real(prec) :: rate_now 

        ! Get current time 
        call cpu_time(cpu_time1)

        if (model_time1 .gt. model_time0) then 
            ! Model has advanced in time, calculate rate 

            ! Calculate the model speed [model-yr / sec]
            rate_now = (model_time1-model_time0) / (cpu_time1-cpu_time0)

            ! Convert to more useful rate [model-kyr / hr]
            rate_now = rate_now*1e-3*3600.0 

        else 
            rate_now = 0.0 
        end if 

        ! Shift rates vector to eliminate oldest entry, and add current entry
        n = size(rates) 
        rates    = cshift(rates,1)
        rates(n) = rate_now 

        ! Avoid underflows
        where(abs(rates) .lt. tol_underflow) rates = 0.0_prec
        
        ! Calculate running average rate 
        n    = count(rates .gt. 0.0)
        if (n .gt. 0) then 
            rate = sum(rates,mask=rates .gt. 0.0) / real(n,prec)
        else 
            rate = 0.0 
        end if 

        return 

    end subroutine yelmo_calc_speed 
    
end module yelmo_defs

