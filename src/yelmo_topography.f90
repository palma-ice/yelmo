
module yelmo_topography

    use nml  
    use ncio
    
    use yelmo_defs
    use yelmo_tools 
    
    use mass_conservation
    use calving
    use topography 
    
    implicit none
    
    ! Key for matching bed types given by mask_bed 
    integer, parameter :: mask_bed_ocean  = 0 
    integer, parameter :: mask_bed_land   = 1
    integer, parameter :: mask_bed_frozen = 2
    integer, parameter :: mask_bed_stream = 3
    integer, parameter :: mask_bed_grline = 4
    integer, parameter :: mask_bed_float  = 5
    integer, parameter :: mask_bed_island = 6
    
    private
    public :: calc_ytopo, calc_ytopo_masks
    public :: ytopo_load_H_ice
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

    subroutine calc_ytopo(tpo,dyn,thrm,bnd,time,topo_fixed)
        ! Calculate adjustments to surface elevation, bedrock elevation
        ! and ice thickness 

        implicit none 

        type(ytopo_class),  intent(INOUT) :: tpo
        type(ydyn_class),   intent(IN)    :: dyn
        type(ytherm_class), intent(IN)    :: thrm  
        type(ybound_class), intent(IN)    :: bnd 
        real(prec),         intent(IN)    :: time  
        logical,            intent(IN)    :: topo_fixed 

        ! Local variables 
        real(prec) :: dx, dt   
        integer :: i, j, nx, ny  
        real(prec), allocatable :: mbal(:,:) 
        
        nx = size(tpo%now%H_ice,1)
        ny = size(tpo%now%H_ice,2)

        allocate(mbal(nx,ny))

        ! Initialize time if necessary 
        if (tpo%par%time .gt. time) then 
            tpo%par%time = time 
        end if 

        ! Get time step
        dt = time - tpo%par%time 

        
        ! Combine basal mass balance into one field accounting for 
        ! grounded/floating fraction of grid cells 
        call calc_bmb_total(tpo%now%bmb,thrm%now%bmb_grnd,bnd%bmb_shlf,tpo%now%f_grnd,tpo%par%diffuse_bmb_shlf)
        
        
        ! Perform topography calculations 
        if ( .not. topo_fixed .and. dt .gt. 0.0 ) then 
            
            ! Define temporary variable for total column mass balance 
           
            mbal = bnd%smb + tpo%now%bmb           
            
            if (.not. tpo%par%use_bmb) then
                ! WHEN RUNNING EISMINT1 ensure bmb is not accounted for here !!!
                mbal = bnd%smb 
            end if 
            
            ! 1. Calculate the ice thickness conservation -----
            call calc_ice_thickness(tpo%now%H_ice,tpo%now%H_margin,tpo%now%f_ice,tpo%now%mb_applied, &
                                    tpo%now%dHdt_n,tpo%now%H_ice_n,tpo%now%H_ice_pred, &
                                    tpo%now%f_grnd,bnd%z_sl-bnd%z_bed,dyn%now%ux_bar,dyn%now%uy_bar, &
                                    mbal=mbal,calv=tpo%now%calv,z_bed_sd=bnd%z_bed_sd,dx=tpo%par%dx,dt=dt, &
                                    solver=trim(tpo%par%solver),boundaries=trim(tpo%par%boundaries), &
                                    ice_allowed=bnd%ice_allowed,H_min=tpo%par%H_min_grnd, &
                                    sd_min=tpo%par%sd_min,sd_max=tpo%par%sd_max,calv_max=tpo%par%calv_max, &
                                    beta1=tpo%par%dt_beta1,beta2=tpo%par%dt_beta2,pc_step=tpo%par%pc_step)
            
            ! If desired, relax solution to reference state
            if (tpo%par%topo_rel .ne. 0) then 

                call relax_ice_thickness(tpo%now%H_ice,tpo%now%f_grnd,bnd%H_ice_ref, &
                                            tpo%par%topo_rel,tpo%par%topo_rel_tau,dt)
                
            end if 


            ! ====== CALVING ======

            ! Diagnose calving rate [m/a]
            select case(trim(tpo%par%calv_method))

                case("zero")

                    tpo%now%calv = 0.0 

                case("simple") 
                    ! Use simple threshold method

                    call calc_calving_rate_simple(tpo%now%calv,tpo%now%H_ice,tpo%now%f_grnd,tpo%now%f_ice, &
                                                    tpo%par%calv_H_lim,tpo%par%calv_tau)
                    
                case("flux") 
                    ! Use threshold+flux method from GRISLI 

                    call calc_calving_rate_flux(tpo%now%calv,tpo%now%H_ice,tpo%now%f_grnd,tpo%now%f_ice,mbal,dyn%now%ux_bar, &
                                                dyn%now%uy_bar,tpo%par%dx,tpo%par%calv_H_lim,tpo%par%calv_tau)
                    
                case("kill") 
                    ! Delete all floating ice (using characteristic time parameter)
                    call calc_calving_rate_kill(tpo%now%calv,tpo%now%H_ice,tpo%now%f_grnd.eq.0.0_prec,tpo%par%calv_tau)

                case("kill-pos")
                    ! Delete all floating ice beyond a given location (using characteristic time parameter)

                    call calc_calving_rate_kill(tpo%now%calv,tpo%now%H_ice, &
                                                    ( tpo%now%f_grnd .eq. 0.0_prec .and. &
                                                      tpo%now%H_ice  .gt. 0.0_prec .and. &
                                                      bnd%calv_mask ), tpo%par%calv_tau)

                case DEFAULT 

                    write(*,*) "calc_ytopo:: Error: calving method not recognized."
                    write(*,*) "calv_method = ", trim(tpo%par%calv_method)
                    stop 

            end select
 
            ! Apply calving
            call apply_calving(tpo%now%H_ice,tpo%now%calv,tpo%now%f_grnd,tpo%par%H_min_flt,dt)

            ! Apply special case for symmetric EISMINT domain when basal sliding is active
            ! (ensure summit thickness does not grow disproportionately)
            if (trim(tpo%par%boundaries) .eq. "EISMINT" .and. maxval(dyn%now%uxy_b) .gt. 0.0) then 
                i = (tpo%par%nx-1)/2 
                j = (tpo%par%ny-1)/2
                tpo%now%H_ice(i,j) = (tpo%now%H_ice(i-1,j)+tpo%now%H_ice(i+1,j) &
                                        +tpo%now%H_ice(i,j-1)+tpo%now%H_ice(i,j+1)) / 4.0 
            end if  
             
            ! Determine the rate of change of ice thickness [m/a]
            tpo%now%dHicedt = (tpo%now%H_ice - tpo%now%H_ice_n)/dt

        else 

            ! Set rate of change to zero 
            tpo%now%dHicedt = 0.0 

        end if 

!         ! 2. Calculate new masks ------------------------------

!         ! Calculate grounding overburden ice thickness 
!         call calc_H_grnd(tpo%now%H_grnd,tpo%now%H_ice,bnd%z_bed,bnd%z_sl)


!         ! Calculate the grounded fraction and grounding line mask of each grid cell
!         select case(tpo%par%gl_sep)

!             case(1) 
!                 ! Binary f_grnd, linear f_grnd_acx/acy based on H_grnd

!                 call calc_f_grnd_subgrid_linear(tpo%now%f_grnd,tpo%now%f_grnd_acx,tpo%now%f_grnd_acy,tpo%now%H_grnd)

!             case(2)
!                 ! Grounded area f_gnrd, average to f_grnd_acx/acy 

!                 call calc_f_grnd_subgrid_area(tpo%now%f_grnd,tpo%now%f_grnd_acx,tpo%now%f_grnd_acy,tpo%now%H_grnd,tpo%par%gl_sep_nx)
            
!         end select
        
! !         ! Filter f_grnd to avoid lakes of one grid point inside of grounded ice 
! !         ! ajr: note, this should be improved to treat ac-nodes too 
! !         call filter_f_grnd(tpo%now%f_grnd)

!         ! Calculate the grounding line mask 
!         call calc_grline(tpo%now%is_grline,tpo%now%is_grz,tpo%now%f_grnd)

!         ! Calculate the bed mask
!         tpo%now%mask_bed = gen_mask_bed(tpo%now%H_ice,thrm%now%f_pmp,tpo%now%f_grnd,tpo%now%is_grline)


        ! 3. Calculate additional topographic properties ------------------

        ! Store previous surface elevation 
        tpo%now%dzsrfdt = tpo%now%z_srf
        
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
            tpo%now%dzsrfdt = (tpo%now%z_srf-tpo%now%dzsrfdt) / dt 
        else 
            tpo%now%dzsrfdt = 0.0 
        end if 
        
        ! Calculate the surface slope (on staggered Ac x/y nodes)
        call calc_gradient_ac_ice(tpo%now%dzsdx,tpo%now%dzsdy,tpo%now%z_srf,tpo%now%H_ice,tpo%par%dx,tpo%par%margin2nd,tpo%par%grad_lim)
        call calc_gradient_ac_ice(tpo%now%dHicedx,tpo%now%dHicedy,tpo%now%H_ice,tpo%now%H_ice,tpo%par%dx,tpo%par%margin2nd,tpo%par%grad_lim)
        
        ! ajr: experimental, doesn't seem to work properly yet! ===>
        ! Modify surface slope gradient at the grounding line if desired 
!         call calc_gradient_ac_gl(tpo%now%dzsdx,tpo%now%dzsdy,tpo%now%z_srf,tpo%now%H_ice, &
!                                       tpo%now%f_grnd_acx,tpo%now%f_grnd_acy,tpo%par%dx,method=2,grad_lim=tpo%par%grad_lim)

        ! Calculate distance to ice margin (really slow if always on)
        !tpo%now%dist_margin = distance_to_margin(tpo%now%H_ice,tpo%par%dx)

        ! Calculate distance to grounding line (really slow if always on)
        !tpo%now%dist_grline = distance_to_grline(tpo%now%is_grline,tpo%now%f_grnd,tpo%par%dx)

        ! Advance ytopo timestep 
        tpo%par%time = time

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

    subroutine calc_ytopo_masks(tpo,dyn,thrm,bnd)
        ! Calculate adjustments to surface elevation, bedrock elevation
        ! and ice thickness 

        implicit none 

        type(ytopo_class),  intent(INOUT) :: tpo
        type(ydyn_class),   intent(IN)    :: dyn
        type(ytherm_class), intent(IN)    :: thrm  
        type(ybound_class), intent(IN)    :: bnd 

        ! Local variables 
        real(prec) :: dx   

        ! 2. Calculate new masks ------------------------------

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
        tpo%now%mask_bed = gen_mask_bed(tpo%now%H_ice,thrm%now%f_pmp,tpo%now%f_grnd,tpo%now%is_grline)


        ! Calculate distance to ice margin (really slow if always on)
        !tpo%now%dist_margin = distance_to_margin(tpo%now%H_ice,tpo%par%dx)

        ! Calculate distance to grounding line (really slow if always on)
        !tpo%now%dist_grline = distance_to_grline(tpo%now%is_grline,tpo%now%f_grnd,tpo%par%dx)

        return 

    end subroutine calc_ytopo_masks

    subroutine ytopo_load_H_ice(tpo,nml_path,nml_group,domain,grid_name,ice_allowed)
        ! Update the topography of the yelmo domain 

        implicit none

        type(ytopo_class), intent(INOUT) :: tpo
        character(len=*), intent(IN)     :: nml_path, nml_group
        character(len=*), intent(IN)     :: domain, grid_name 
        logical,          intent(IN)     :: ice_allowed(:,:) 

        ! Local variables
        logical            :: load_var
        character(len=512) :: filename
        character(len=56)  :: vname 

        ! == H_ice =====
        call nml_read(nml_path,nml_group,"H_ice_load",load_var)

        if (load_var) then 
            call nml_read(nml_path,nml_group,"H_ice_path",filename)
            call yelmo_parse_path(filename, domain,grid_name)
            call nml_read(nml_path,nml_group,"H_ice_nm",  vname)
            call nc_read(filename,vname,tpo%now%H_ice)
        else 
            tpo%now%H_ice = 0.0  
        end if 
        
        ! Clean up field 
        where(tpo%now%H_ice  .lt. 1.0) tpo%now%H_ice = 0.0 

        ! Artificially delete ice from locations that are not allowed
        where (.not. ice_allowed) tpo%now%H_ice = 0.0 
        
        write(*,*) "ytopo_load_H_ice:: range(H_ice):  ", minval(tpo%now%H_ice),  maxval(tpo%now%H_ice)

        return 

    end subroutine ytopo_load_H_ice 

    elemental function gen_mask_bed(H_ice,f_pmp,f_grnd,is_grline) result(mask)
        ! Generate an output mask for model conditions at bed
        ! based on input masks 
        ! 0: ocean, 1: land, 2: sia, 3: streams, grline: 4, floating: 5, islands: 6

        implicit none 

        real(prec), intent(IN) :: H_ice, f_pmp, f_grnd
        logical, intent(IN) :: is_grline
        integer :: mask

        if (is_grline) then
            mask = mask_bed_grline        ! Grounding line

        else if (f_grnd .gt. 0.0 .and. H_ice .eq. 0.0) then 
            mask = mask_bed_land        ! Ice-free land

        else if (f_grnd .gt. 0.0 .and. f_pmp .lt. 0.5) then 
            mask = mask_bed_frozen        ! Inland frozen bed

        else if (f_grnd .gt. 0.0) then 
            mask = mask_bed_stream        ! Inland stream

        else if (f_grnd .eq. 0.0 .and. H_ice .gt. 0.0) then 
            mask = mask_bed_float        ! Floating ice shelves

        else 
            mask = mask_bed_ocean        ! Ocean 

        end if 

        return 

    end function gen_mask_bed 

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
        call nml_read(filename,"ytopo","calv_method",       par%calv_method,      init=init_pars)
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

        ! Intialize timestepping parameters to Forward Euler (beta2=0: no contribution from previous timestep)
        par%dt_zeta     = 1.0 
        par%dt_beta1    = 1.0 
        par%dt_beta2    = 0.0 

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
        allocate(now%mb_applied(nx,ny))
        allocate(now%calv_grnd(nx,ny))
        allocate(now%calv(nx,ny))

        allocate(now%H_margin(nx,ny))
        
        allocate(now%dzsdx(nx,ny))
        allocate(now%dzsdy(nx,ny))

        allocate(now%dHicedx(nx,ny))
        allocate(now%dHicedy(nx,ny))
        
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
        
        now%H_ice       = 0.0 
        now%z_srf       = 0.0  
        now%dzsrfdt     = 0.0 
        now%dHicedt     = 0.0
        now%bmb         = 0.0  
        now%mb_applied  = 0.0 
        now%calv_grnd   = 0.0
        now%calv        = 0.0
        now%H_margin    = 0.0 
        now%dzsdx       = 0.0 
        now%dzsdy       = 0.0 
        now%dHicedx     = 0.0 
        now%dHicedy     = 0.0
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
        if (allocated(now%mb_applied))  deallocate(now%mb_applied)
        if (allocated(now%calv_grnd))   deallocate(now%calv_grnd)
        if (allocated(now%calv))        deallocate(now%calv)

        if (allocated(now%H_margin))    deallocate(now%H_margin)
        
        if (allocated(now%dzsdx))       deallocate(now%dzsdx)
        if (allocated(now%dzsdy))       deallocate(now%dzsdy)
        if (allocated(now%dHicedx))     deallocate(now%dHicedx)
        if (allocated(now%dHicedy))     deallocate(now%dHicedy)
        
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
        
        return 

    end subroutine ytopo_dealloc 
    
end module yelmo_topography
