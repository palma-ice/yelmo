
module yelmo_topography

    use nml  
    use ncio
    
    use yelmo_defs
    use yelmo_tools 
    
    use mass_conservation
    use calving
    use topography 
    use deformation 

    implicit none
    
    ! Key for matching bed types given by mask_bed 
    integer, parameter :: mask_bed_ocean   = 0 
    integer, parameter :: mask_bed_land    = 1
    integer, parameter :: mask_bed_frozen  = 2
    integer, parameter :: mask_bed_stream  = 3
    integer, parameter :: mask_bed_grline  = 4
    integer, parameter :: mask_bed_float   = 5
    integer, parameter :: mask_bed_island  = 6
    integer, parameter :: mask_bed_partial = 7

    private
    public :: calc_ytopo
    public :: ytopo_par_load, ytopo_alloc, ytopo_dealloc
     
    ! Integers
    public :: mask_bed_ocean  
    public :: mask_bed_land  
    public :: mask_bed_frozen
    public :: mask_bed_stream
    public :: mask_bed_grline
    public :: mask_bed_float 
    public :: mask_bed_island
    
contains

    subroutine calc_ytopo(tpo,dyn,mat,thrm,bnd,time,topo_fixed)
        ! Calculate adjustments to surface elevation, bedrock elevation
        ! and ice thickness 

        implicit none 

        type(ytopo_class),  intent(INOUT) :: tpo
        type(ydyn_class),   intent(IN)    :: dyn
        type(ymat_class),   intent(IN)    :: mat
        type(ytherm_class), intent(IN)    :: thrm  
        type(ybound_class), intent(IN)    :: bnd 
        real(prec),         intent(IN)    :: time  
        logical,            intent(IN)    :: topo_fixed 

        ! Local variables 
        real(prec) :: dx, dt   
        integer :: i, j, nx, ny  
        real(prec), allocatable :: mbal(:,:) 
        real(prec), allocatable :: H_ref(:,:) 

        type(strain_2D_class) :: strn2D 
        type(stress_2D_class) :: strs2D

        real(8)    :: cpu_time0, cpu_time1
        real(prec) :: model_time0, model_time1 
        real(prec) :: speed 

        nx = size(tpo%now%H_ice,1)
        ny = size(tpo%now%H_ice,2)

        allocate(mbal(nx,ny))
        allocate(H_ref(nx,ny))

        ! Strain and stress 
        allocate(strn2D%dxx(nx,ny))
        allocate(strn2D%dyy(nx,ny))
        allocate(strn2D%dxy(nx,ny))
        allocate(strn2D%dxz(nx,ny))
        allocate(strn2D%dyz(nx,ny))
        allocate(strn2D%de(nx,ny))
        allocate(strn2D%f_shear(nx,ny))

        allocate(strs2D%txx(nx,ny))
        allocate(strs2D%tyy(nx,ny))
        allocate(strs2D%txy(nx,ny))
        allocate(strs2D%txz(nx,ny))
        allocate(strs2D%tyz(nx,ny))
        allocate(strs2D%te(nx,ny))
        allocate(strs2D%teig1(nx,ny))
        allocate(strs2D%teig2(nx,ny))
        
        ! Initialize time if necessary 
        if (tpo%par%time .gt. dble(time)) then 
            tpo%par%time = dble(time) 
        end if 

        ! Get time step
        dt = dble(time) - tpo%par%time 

        ! Store initial cpu time and model time for metrics later
        call yelmo_cpu_time(cpu_time0)
        model_time0 = tpo%par%time 

        ! Combine basal mass balance into one field accounting for 
        ! grounded/floating fraction of grid cells 
        call calc_bmb_total(tpo%now%bmb,thrm%now%bmb_grnd,bnd%bmb_shlf,tpo%now%f_grnd,tpo%par%diffuse_bmb_shlf)
        
        ! Combine frontal mass balance into one field, and 
        ! calculate as needed 
        call calc_fmb_total(tpo%now%fmb,bnd%fmb_shlf,bnd%bmb_shlf,tpo%now%H_ice, &
                        tpo%now%H_grnd,tpo%now%f_ice,tpo%par%fmb_method,tpo%par%fmb_scale,tpo%par%dx)

        ! 1. Perform topography calculations ------------------

        if ( .not. topo_fixed .and. dt .gt. 0.0 ) then 

            ! === Step 1: ice thickness evolution from dynamics alone ===

            ! Calculate the ice thickness conservation from dynamics only -----
            call calc_ice_thickness_dyn(tpo%now%H_ice,tpo%now%dHdt_n,tpo%now%H_ice_n,tpo%now%H_ice_pred, &
                                        tpo%now%f_ice,tpo%now%f_grnd,dyn%now%ux_bar,dyn%now%uy_bar, &
                                        solver=tpo%par%solver,dx=tpo%par%dx,dt=dt,beta=tpo%par%dt_beta,pc_step=tpo%par%pc_step)

            ! Update ice fraction mask 
            call calc_ice_fraction(tpo%now%f_ice,tpo%now%H_ice,tpo%now%f_grnd)
            
            ! If desired, relax solution to reference state
            ! ajr: why is this here? Shouldn't it go after all dyn+calv+mb steps applied?
            ! Test two options to see.
            if (tpo%par%topo_rel .ne. 0) then 

                call relax_ice_thickness(tpo%now%H_ice,tpo%now%f_grnd,bnd%H_ice_ref, &
                                            tpo%par%topo_rel,tpo%par%topo_rel_tau,dt)
                
            end if 

            ! === Step 2: ice thickness evolution from vertical column mass balance ===

            tpo%now%mb_applied  = 0.0_wp  
            tpo%now%calv        = 0.0_wp 

            ! Also, define temporary variable for total column mass balance (without calving)
           
            mbal = bnd%smb + tpo%now%bmb + tpo%now%fmb 
            
            if (.not. tpo%par%use_bmb) then
                ! WHEN RUNNING EISMINT1 ensure bmb is not accounted for here !!!
                mbal = bnd%smb + tpo%now%fmb  
            end if 
            
if (.FALSE.) then 

            ! Now, apply mass-conservation step 2: mass balance (no calving terms)
            call calc_ice_thickness_mbal(tpo%now%H_ice,tpo%now%mb_applied,tpo%now%calv,tpo%now%f_ice, &
                                         tpo%now%f_grnd,bnd%z_sl-bnd%z_bed,dyn%now%ux_bar,dyn%now%uy_bar, &
                                         mbal,tpo%now%calv_flt*0.0_wp,tpo%now%calv_grnd*0.0_wp,bnd%z_bed_sd,tpo%par%dx,dt)

            ! Update ice fraction mask 
            call calc_ice_fraction(tpo%now%f_ice,tpo%now%H_ice,tpo%now%f_grnd)

end if 

            ! Next, diagnose CALVING

            ! Diagnose potential floating-ice calving rate [m/a]
            select case(trim(tpo%par%calv_flt_method))

                case("zero")

                    tpo%now%calv_flt = 0.0 

                case("simple") 
                    ! Use simple threshold method

                    call calc_calving_rate_simple(tpo%now%calv_flt,tpo%now%H_ice,tpo%now%f_grnd,tpo%now%f_ice, &
                                                    tpo%par%calv_H_lim,tpo%par%calv_tau)
                    
                case("flux") 
                    ! Use threshold+flux method from Peyaud et al. (2007), ie, GRISLI 

                    call calc_calving_rate_flux(tpo%now%calv_flt,tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd,mbal,dyn%now%ux_bar, &
                                                dyn%now%uy_bar,tpo%par%dx,tpo%par%calv_H_lim,tpo%par%calv_tau)
                
                case("vm-l19")
                    ! Use von Mises calving as defined by Lipscomb et al. (2019)

                    call calc_calving_rate_vonmises_l19(tpo%now%calv_flt,tpo%now%H_ice,tpo%now%f_ice,tpo%now%f_grnd, &
                                                        mat%now%strs2D%teig1,mat%now%strs2D%teig2,mat%now%ATT_bar, &
                                                        mat%now%visc_bar,tpo%par%dx,tpo%par%dx,tpo%par%kt,tpo%par%w2,mat%par%n_glen)

                case("kill") 
                    ! Delete all floating ice (using characteristic time parameter)
                    call calc_calving_rate_kill(tpo%now%calv_flt,tpo%now%H_ice,tpo%now%f_grnd.eq.0.0_prec,tpo%par%calv_tau,dt)

                case("kill-pos")
                    ! Delete all floating ice beyond a given location (using characteristic time parameter)

                    call calc_calving_rate_kill(tpo%now%calv_flt,tpo%now%H_ice, &
                                                    ( tpo%now%f_grnd .eq. 0.0_prec .and. &
                                                      tpo%now%H_ice  .gt. 0.0_prec .and. &
                                                      bnd%calv_mask ), tau=0.0_prec, dt=dt )

                case DEFAULT 

                    write(*,*) "calc_ytopo:: Error: calving method not recognized."
                    write(*,*) "calv_flt_method = ", trim(tpo%par%calv_flt_method)
                    stop 

            end select
            

            ! Diagnose potential grounded-ice calving rate [m/a]

            ! For now, set it to zero 
            tpo%now%calv_grnd = 0.0

if (.FALSE.) then 
            ! Now, apply mass-conservation step 3: calving (no mbal term)
            call calc_ice_thickness_mbal(tpo%now%H_ice,tpo%now%mb_applied,tpo%now%calv,tpo%now%f_ice, &
                                         tpo%now%f_grnd,bnd%z_sl-bnd%z_bed,dyn%now%ux_bar,dyn%now%uy_bar, &
                                         mbal*0.0,tpo%now%calv_flt,tpo%now%calv_grnd,bnd%z_bed_sd,tpo%par%dx,dt)

else

            ! Apply mass-conservation step (mbal and calving together)
            call calc_ice_thickness_mbal(tpo%now%H_ice,tpo%now%mb_applied,tpo%now%calv,tpo%now%f_ice, &
                                         tpo%now%f_grnd,bnd%z_sl-bnd%z_bed,dyn%now%ux_bar,dyn%now%uy_bar, &
                                         mbal,tpo%now%calv_flt,tpo%now%calv_grnd,bnd%z_bed_sd,tpo%par%dx,dt)


end if 

            ! Update ice fraction mask 
            call calc_ice_fraction(tpo%now%f_ice,tpo%now%H_ice,tpo%now%f_grnd)
            
            ! Finally, apply all additional (generally artificial) ice thickness adjustments 
            ! and store changes in residual mass balance field. 
            call apply_ice_thickness_boundaries(tpo%now%H_ice,tpo%now%mb_resid,tpo%now%f_ice,tpo%now%f_grnd, &
                                                dyn%now%uxy_b,bnd%ice_allowed,tpo%par%boundaries,bnd%H_ice_ref, &
                                                tpo%par%H_min_flt,tpo%par%H_min_grnd,dt)

            ! Update ice fraction mask 
            call calc_ice_fraction(tpo%now%f_ice,tpo%now%H_ice,tpo%now%f_grnd)
            
            ! Save the rate of change of ice thickness in output variable [m/a]
            tpo%now%dHicedt = (tpo%now%H_ice - tpo%now%H_ice_n) / dt 

        end if 

        ! 2. Calculate additional topographic properties ------------------

        ! Store previous surface elevation on predictor step for calculating
        ! rate of change of surface elevation.
        if (trim(tpo%par%pc_step) .eq. "predictor") then 
            tpo%now%z_srf_n = tpo%now%z_srf 
        end if 

        ! Calculate the surface elevation 
        select case(tpo%par%surf_gl_method)
            ! Choose method to treat grounding line points when calculating surface elevation

            case(0)
                ! Binary (grounded elevation or floating elevation via archemedes)
                ! Note: two functions that should give the same results
                
                !call calc_z_srf(tpo%now%z_srf,tpo%now%H_ice,tpo%now%H_grnd,bnd%z_bed,bnd%z_sl)
                call calc_z_srf_max(tpo%now%z_srf,tpo%now%H_ice,bnd%z_bed,bnd%z_sl)
            
            case(1)
                ! Subgrid z_srf calculations at the grounding line 

                call calc_z_srf_subgrid_area(tpo%now%z_srf,tpo%now%f_grnd,tpo%now%H_ice,bnd%z_bed,bnd%z_sl,tpo%par%gl_sep_nx)

        end select 

        ! Determine rate of surface elevation change 
        if (dt .gt. 0.0) then 
            tpo%now%dzsrfdt = (tpo%now%z_srf-tpo%now%z_srf_n) / dt 
        end if 
        
        ! Calculate the surface slope (on staggered Ac x/y nodes)
        call calc_gradient_ac_ice(tpo%now%dzsdx,tpo%now%dzsdy,tpo%now%z_srf,tpo%now%f_ice,tpo%par%dx, &
                                                tpo%par%margin2nd,tpo%par%grad_lim,tpo%par%boundaries)
        call calc_gradient_ac_ice(tpo%now%dHicedx,tpo%now%dHicedy,tpo%now%H_ice,tpo%now%f_ice,tpo%par%dx, &
                                                tpo%par%margin2nd,tpo%par%grad_lim,tpo%par%boundaries)
        
        ! ajr: experimental, doesn't seem to work properly yet! ===>
        ! Modify surface slope gradient at the grounding line if desired 
!         call calc_gradient_ac_gl(tpo%now%dzsdx,tpo%now%dzsdy,tpo%now%z_srf,tpo%now%H_ice, &
!                                       tpo%now%f_grnd_acx,tpo%now%f_grnd_acy,tpo%par%dx,method=2,grad_lim=tpo%par%grad_lim)
        
        ! 3. Calculate new masks ------------------------------

        ! Calculate grounding overburden ice thickness 
        call calc_H_grnd(tpo%now%H_grnd,tpo%now%H_ice,bnd%z_bed,bnd%z_sl)


        ! Calculate the grounded fraction and grounding line mask of each grid cell
        select case(tpo%par%gl_sep)

            case(1) 
                ! Binary f_grnd, linear f_grnd_acx/acy based on H_grnd

                call calc_f_grnd_subgrid_linear(tpo%now%f_grnd,tpo%now%f_grnd_acx,tpo%now%f_grnd_acy,tpo%now%H_grnd)

            case(2)
                ! Grounded area f_gnrd, average to f_grnd_acx/acy 

                call calc_f_grnd_subgrid_area(tpo%now%f_grnd,tpo%now%f_grnd_acx,tpo%now%f_grnd_acy,tpo%now%H_grnd,tpo%par%gl_sep_nx)
            
        end select
        
!         ! Filter f_grnd to avoid lakes of one grid point inside of grounded ice 
!         ! ajr: note, this should be improved to treat ac-nodes too 
!         call filter_f_grnd(tpo%now%f_grnd)

        ! Calculate the grounding line mask 
        call calc_grline(tpo%now%is_grline,tpo%now%is_grz,tpo%now%f_grnd)

        ! Calculate the bed mask
        call gen_mask_bed(tpo%now%mask_bed,tpo%now%f_ice,thrm%now%f_pmp,tpo%now%f_grnd,tpo%now%is_grline)

        ! Calculate distance to ice margin (really slow if always on)
        !tpo%now%dist_margin = distance_to_margin(tpo%now%H_ice,tpo%par%dx)

        ! Calculate distance to grounding line (really slow if always on)
        !tpo%now%dist_grline = distance_to_grline(tpo%now%is_grline,tpo%now%f_grnd,tpo%par%dx)

        ! ================================

        ! Calculate computational performance (model speed in kyr/hr)
        call yelmo_cpu_time(cpu_time1)
        model_time1 = tpo%par%time 
        call yelmo_calc_speed(speed,model_time0,model_time1,cpu_time0,cpu_time1)

        ! Store the speed variable in predictor or corrector speed variable
        if (trim(tpo%par%pc_step) .eq. "predictor") then 
            tpo%par%speed_pred = speed 
        else 
            tpo%par%speed_corr = speed 

            ! If corrector step, then also calculate the speed of both 
            ! predictor+corrector calls: mean of the predictor and corrector speeds
            ! divided by two, since two calls were made. 
            tpo%par%speed = 0.5 * (0.5*(tpo%par%speed_pred+tpo%par%speed_corr))
            
        end if 

        if (trim(tpo%par%pc_step) .eq. "corrector") then 
            ! Advance ytopo timestep on corrector step 

            tpo%par%time = dble(time)
            
        end if 

!         if (yelmo_log) then 

!             if (count(tpo%now%H_ice.gt.0.0) .gt. 0) then 
!                 write(*,"(a,f14.4,f10.4,f10.2)") "calc_ytopo::  time = ", tpo%par%time, dt, &
!                     sum(tpo%now%H_ice,mask=tpo%now%H_ice.gt.0.0)/real(count(tpo%now%H_ice.gt.0.0))
!             else 
!                 write(*,"(a,f14.4,f10.4,f10.2)") "calc_ytopo::  time = ", tpo%par%time, dt, 0.0
!             end if 

!         end if 

        return 

    end subroutine calc_ytopo

    elemental subroutine gen_mask_bed(mask,f_ice,f_pmp,f_grnd,is_grline)
        ! Generate an output mask for model conditions at bed
        ! based on input masks 
        ! 0: ocean, 1: land, 2: sia, 3: streams, grline: 4, floating: 5, islands: 6
        ! 7: partially-covered ice cell.

        implicit none 

        integer,    intent(OUT) :: mask 
        real(prec), intent(IN)  :: f_ice, f_pmp, f_grnd
        logical,    intent(IN)  :: is_grline

        if (is_grline) then
            ! Grounding line

            mask = mask_bed_grline

        else if (f_ice .eq. 0.0) then 
            ! Ice-free points 

            if (f_grnd .gt. 0.0) then
                ! Ice-free land

                mask = mask_bed_land

            else
                ! Ice-free ocean

                mask = mask_bed_ocean

            end if 

        else if (f_ice .gt. 0.0 .and. f_ice .lt. 1.0) then 
            ! Partially ice-covered points 

            mask = mask_bed_partial

        else
            ! Fully ice-covered points 

            if (f_grnd .gt. 0.0) then
                ! Grounded ice-covered points 

                if (f_pmp .gt. 0.5) then 
                    ! Temperate points

                    mask = mask_bed_stream 

                else
                    ! Frozen points 

                    mask = mask_bed_frozen 

                end if 

            else
                ! Floating ice-covered points 

                mask = mask_bed_float

            end if 

        end if 

        return 

    end subroutine gen_mask_bed

    subroutine ytopo_par_load(par,filename,nx,ny,dx,init)

        type(ytopo_param_class), intent(OUT) :: par
        character(len=*),        intent(IN)  :: filename
        integer,                  intent(IN)  :: nx, ny 
        real(prec),               intent(IN)  :: dx  
        logical, optional,       intent(IN)  :: init 

        ! Local variables
        logical :: init_pars 

        init_pars = .FALSE.
        if (present(init)) init_pars = .TRUE. 
 
        ! Store parameter values in output object
        call nml_read(filename,"ytopo","solver",            par%solver,           init=init_pars)
        call nml_read(filename,"ytopo","surf_gl_method",    par%surf_gl_method,   init=init_pars)
        call nml_read(filename,"ytopo","calv_flt_method",   par%calv_flt_method,  init=init_pars)
        call nml_read(filename,"ytopo","calv_grnd_method",  par%calv_grnd_method, init=init_pars)
        call nml_read(filename,"ytopo","fmb_method",        par%fmb_method,       init=init_pars)
        call nml_read(filename,"ytopo","margin2nd",         par%margin2nd,        init=init_pars)
        call nml_read(filename,"ytopo","use_bmb",           par%use_bmb,          init=init_pars)
        call nml_read(filename,"ytopo","topo_fixed",        par%topo_fixed,       init=init_pars)
        call nml_read(filename,"ytopo","topo_rel",          par%topo_rel,         init=init_pars)
        call nml_read(filename,"ytopo","topo_rel_tau",      par%topo_rel_tau,     init=init_pars)
        call nml_read(filename,"ytopo","calv_H_lim",        par%calv_H_lim,       init=init_pars)
        call nml_read(filename,"ytopo","calv_tau",          par%calv_tau,         init=init_pars)
        call nml_read(filename,"ytopo","H_min_grnd",        par%H_min_grnd,       init=init_pars)
        call nml_read(filename,"ytopo","H_min_flt",         par%H_min_flt,        init=init_pars)
        call nml_read(filename,"ytopo","sd_min",            par%sd_min,           init=init_pars)
        call nml_read(filename,"ytopo","sd_max",            par%sd_max,           init=init_pars)
        call nml_read(filename,"ytopo","calv_max",          par%calv_max,         init=init_pars)
        call nml_read(filename,"ytopo","grad_lim",          par%grad_lim,         init=init_pars)
        call nml_read(filename,"ytopo","gl_sep",            par%gl_sep,           init=init_pars)
        call nml_read(filename,"ytopo","gl_sep_nx",         par%gl_sep_nx,        init=init_pars)
        call nml_read(filename,"ytopo","diffuse_bmb_shlf",  par%diffuse_bmb_shlf, init=init_pars)
        call nml_read(filename,"ytopo","fmb_scale",         par%fmb_scale,        init=init_pars)
        call nml_read(filename,"ytopo","kt",                par%kt,               init=init_pars)
        call nml_read(filename,"ytopo","w2",                par%w2,               init=init_pars)
        
        ! === Set internal parameters =====

        par%nx  = nx 
        par%ny  = ny 
        par%dx  = dx 
        par%dy  = dx 

        ! Define how boundaries of grid should be treated 
        ! This should only be modified by the dom%par%experiment variable
        ! in yelmo_init. By default set boundaries to zero 
        par%boundaries = "zeros" 
        
        ! Define current time as unrealistic value
        par%time = 1000000000   ! [a] 1 billion years in the future 

        ! Intialize timestepping parameters to Forward Euler (beta2=beta4=0: no contribution from previous timestep)
        par%dt_zeta     = 1.0 
        par%dt_beta(1)  = 1.0 
        par%dt_beta(2)  = 0.0 
        par%dt_beta(3)  = 1.0 
        par%dt_beta(4)  = 0.0 

        return

    end subroutine ytopo_par_load
    
    subroutine ytopo_alloc(now,nx,ny)

        implicit none 

        type(ytopo_state_class), intent(INOUT) :: now 
        integer, intent(IN) :: nx, ny  

        call ytopo_dealloc(now)

        allocate(now%H_ice(nx,ny))
        allocate(now%z_srf(nx,ny))
        allocate(now%dzsrfdt(nx,ny))
        allocate(now%dHicedt(nx,ny))
        allocate(now%bmb(nx,ny))
        allocate(now%fmb(nx,ny))
        allocate(now%mb_applied(nx,ny))
        allocate(now%mb_resid(nx,ny))
        
        allocate(now%calv(nx,ny))
        allocate(now%calv_flt(nx,ny))
        allocate(now%calv_grnd(nx,ny))
        
        allocate(now%dzsdx(nx,ny))
        allocate(now%dzsdy(nx,ny))

        allocate(now%dHicedx(nx,ny))
        allocate(now%dHicedy(nx,ny))
        
        allocate(now%H_eff(nx,ny))
        allocate(now%H_grnd(nx,ny))

        ! Masks 
        allocate(now%f_grnd(nx,ny))
        allocate(now%f_grnd_acx(nx,ny))
        allocate(now%f_grnd_acy(nx,ny))
        allocate(now%f_ice(nx,ny))

        allocate(now%dist_margin(nx,ny))
        allocate(now%dist_grline(nx,ny))
        
        allocate(now%mask_bed(nx,ny))
        allocate(now%is_grline(nx,ny))
        allocate(now%is_grz(nx,ny))

        allocate(now%dHdt_n(nx,ny))
        allocate(now%H_ice_n(nx,ny))
        allocate(now%H_ice_pred(nx,ny))
        allocate(now%H_ice_corr(nx,ny))
        
        allocate(now%z_srf_n(nx,ny))
        
        now%H_ice       = 0.0 
        now%z_srf       = 0.0  
        now%dzsrfdt     = 0.0 
        now%dHicedt     = 0.0
        now%bmb         = 0.0  
        now%fmb         = 0.0
        now%mb_applied  = 0.0 
        now%mb_resid    = 0.0
        now%calv        = 0.0
        now%calv_flt    = 0.0
        now%calv_grnd   = 0.0
        now%dzsdx       = 0.0 
        now%dzsdy       = 0.0 
        now%dHicedx     = 0.0 
        now%dHicedy     = 0.0
        now%H_eff       = 0.0 
        now%H_grnd      = 0.0  
        now%f_grnd      = 0.0  
        now%f_grnd_acx  = 0.0  
        now%f_grnd_acy  = 0.0  
        now%f_ice       = 0.0  
        now%dist_margin = 0.0
        now%dist_grline = 0.0 
        
        now%mask_bed    = 0.0 
        now%is_grline   = .FALSE. 
        now%is_grz      = .FALSE. 
         
        now%dHdt_n      = 0.0  
        now%H_ice_n     = 0.0 
        now%H_ice_pred  = 0.0 
        now%H_ice_corr  = 0.0 
        
        now%z_srf_n     = 0.0 

        return 

    end subroutine ytopo_alloc

    subroutine ytopo_dealloc(now)

        implicit none 

        type(ytopo_state_class), intent(INOUT) :: now

        if (allocated(now%H_ice))       deallocate(now%H_ice)
        if (allocated(now%z_srf))       deallocate(now%z_srf)
        
        if (allocated(now%dzsrfdt))     deallocate(now%dzsrfdt)
        if (allocated(now%dHicedt))     deallocate(now%dHicedt)
        if (allocated(now%bmb))         deallocate(now%bmb)
        if (allocated(now%fmb))         deallocate(now%fmb)
        if (allocated(now%mb_applied))  deallocate(now%mb_applied)
        if (allocated(now%mb_resid))    deallocate(now%mb_resid)
        
        if (allocated(now%calv))        deallocate(now%calv)
        if (allocated(now%calv_flt))    deallocate(now%calv_flt)
        if (allocated(now%calv_grnd))   deallocate(now%calv_grnd)
            
        if (allocated(now%dzsdx))       deallocate(now%dzsdx)
        if (allocated(now%dzsdy))       deallocate(now%dzsdy)
        if (allocated(now%dHicedx))     deallocate(now%dHicedx)
        if (allocated(now%dHicedy))     deallocate(now%dHicedy)
        
        if (allocated(now%H_eff))       deallocate(now%H_eff)
        if (allocated(now%H_grnd))      deallocate(now%H_grnd)

        if (allocated(now%f_grnd))      deallocate(now%f_grnd)
        if (allocated(now%f_grnd_acx))  deallocate(now%f_grnd_acx)
        if (allocated(now%f_grnd_acy))  deallocate(now%f_grnd_acy)
        
        if (allocated(now%f_ice))       deallocate(now%f_ice)

        if (allocated(now%dist_margin)) deallocate(now%dist_margin)
        if (allocated(now%dist_grline)) deallocate(now%dist_grline)
        
        if (allocated(now%mask_bed))    deallocate(now%mask_bed)
        if (allocated(now%is_grline))   deallocate(now%is_grline)
        if (allocated(now%is_grz))      deallocate(now%is_grz)

        if (allocated(now%dHdt_n))      deallocate(now%dHdt_n)
        if (allocated(now%H_ice_n))     deallocate(now%H_ice_n)
        if (allocated(now%H_ice_pred))  deallocate(now%H_ice_pred)
        if (allocated(now%H_ice_corr))  deallocate(now%H_ice_corr)
        
        if (allocated(now%z_srf_n))     deallocate(now%z_srf_n)
        
        return 

    end subroutine ytopo_dealloc
    
end module yelmo_topography
