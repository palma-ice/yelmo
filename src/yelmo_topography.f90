
module yelmo_topography

    use nml  
    use ncio
    
    use yelmo_defs
    use yelmo_tools 
    
    use mass_conservation
    use calving
    use basal_dragging, only : calc_effective_pressure  

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
    public :: calc_ytopo
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
        real(prec) :: dx, dt, dt_calv    
        integer :: i, j, nx, ny  
        real(prec), allocatable :: mbal(:,:) 
        real(prec) :: N_eff_min 

        nx = size(tpo%now%H_ice,1)
        ny = size(tpo%now%H_ice,2)

        allocate(mbal(nx,ny))

        ! Initialize time if necessary 
        if (tpo%par%time .gt. time) then 
            tpo%par%time      = time
            tpo%par%time_calv = time 
        end if 

        ! Get time steps
        dt                = time - tpo%par%time 
        dt_calv           = time - tpo%par%time_calv
        

        
        ! Combine basal mass balance into one field accounting for 
        ! grounded/floating fraction of grid cells 
        call calc_bmb_total(tpo%now%bmb,thrm%now%bmb_grnd,bnd%bmb_shlf,tpo%now%f_grnd,tpo%par%diffuse_bmb_shlf)
        


        ! Perform topography calculations 
        if ( .not. topo_fixed .and. dt .gt. 0.0 ) then 

            ! Store previous ice thickness
            tpo%now%dHicedt = tpo%now%H_ice 
            
            ! Define temporary variable for total column mass balance 
           
            mbal = bnd%smb + tpo%now%bmb           
            
            if (.not. tpo%par%use_bmb) then
                ! WHEN RUNNING EISMINT1 ensure bmb is not accounted for here !!!
                mbal = bnd%smb 
            end if 

            ! 1. Calculate the ice thickness conservation and apply bedrock uplift -----
            call calc_ice_thickness(tpo%now%H_ice,tpo%now%mb_applied, &
                                    tpo%now%f_grnd,dyn%now%ux_bar,dyn%now%uy_bar, &
                                    mbal=mbal,calv=tpo%now%calv*0.0,dx=tpo%par%dx,dt=dt, &
                                    solver=trim(tpo%par%solver),boundaries=trim(tpo%par%boundaries), &
                                    ice_allowed=bnd%ice_allowed,H_min=tpo%par%H_min)
            
            ! ====== CALVING ======
            if (dt_calv .ge. tpo%par%calv_dt) then 
                ! Diagnose calving rate at desired timestep frequency [m/a]

                ! Update the ice-covered fraction to be more precise with calving 
                tpo%now%f_ice = calc_ice_fraction(tpo%now%H_ice,tpo%now%f_grnd)
                
                select case(trim(tpo%par%calv_method))

                    case("zero")

                        tpo%now%calv = 0.0 

                    case("simple") 
                        ! Use simple threshold method

                        tpo%now%calv  = calc_calving_rate_simple(tpo%now%H_ice,tpo%now%f_grnd,tpo%now%f_ice,dt_calv,tpo%par%H_calv)
                    
                    case("flux") 
                        ! Use threshold+flux method from GRISLI 

                        tpo%now%calv  = calc_calving_rate_flux(tpo%now%H_ice,tpo%now%f_grnd,tpo%now%f_ice, &
                                               mbal,dyn%now%ux_bar,dyn%now%uy_bar,tpo%par%dx,dt_calv,tpo%par%H_calv)
                    
                    case("kill") 
                        ! Delete all floating ice 
                        tpo%now%calv = calc_calving_rate_kill(tpo%now%H_ice,tpo%now%f_grnd,dt_calv)

                    case DEFAULT 

                        write(*,*) "calc_ytopo:: Error: calving method not recognized."
                        write(*,*) "calv_method = ", trim(tpo%par%calv_method)
                        stop 

                end select

                ! Apply calving rate, ensure only available ice is deleted 
                call apply_calving(tpo%now%H_ice,tpo%now%calv,dt_calv,H_min=tpo%par%H_min)
                
                ! Updating current calving time 
                tpo%par%time_calv = time 

                
            else 
                ! dt too small, set calving to zero 
                tpo%now%calv = 0.0 

            end if 

            ! Determine the rate of change of ice thickness [m/a]
            tpo%now%dHicedt = (tpo%now%H_ice - tpo%now%dHicedt)/dt

        else 

            ! Set rate of change to zero 
            tpo%now%dHicedt = 0.0 

        end if 

        ! 2. Calculate new masks ------------------------------

        ! Calculate grounding overburden ice thickness 
        call calc_H_grnd(tpo%now%H_grnd,tpo%now%H_ice,bnd%z_bed,bnd%z_sl)


        ! Calculate the grounded fraction and grounding line mask of each grid cell
        select case(tpo%par%gl_sep)

            case(1) 
                ! Binary f_grnd, linear f_grnd_acx/acy based on H_grnd

                call calc_f_grnd_subgrid_ac_linear(tpo%now%f_grnd,tpo%now%f_grnd_acx,tpo%now%f_grnd_acy,tpo%now%H_grnd)

            case(2)
                ! Grounded area f_gnrd, average to f_grnd_acx/acy 

                call calc_f_grnd_subgrid_area(tpo%now%f_grnd,tpo%now%f_grnd_acx,tpo%now%f_grnd_acy,tpo%now%H_grnd,tpo%par%gl_sep_nx)
            
        end select
        
!         ! Filter f_grnd to avoid lakes of one grid point inside of grounded ice 
!         ! ajr: note, this should be improved to treat ac-nodes too 
!         call filter_f_grnd(tpo%now%f_grnd)

        ! Calculate the grounding line mask 
        call calc_grline(tpo%now%is_grline,tpo%now%f_grnd)

        ! Calculate the ice-covered fraction of each grid cell 
        ! ajr: note, this should be improved to treat both floating and grounded ice points...
        tpo%now%f_ice = calc_ice_fraction(tpo%now%H_ice,tpo%now%f_grnd)
        
        ! Calculate the bed mask
        tpo%now%mask_bed = gen_mask_bed(tpo%now%H_ice,thrm%now%f_pmp,tpo%now%f_grnd,tpo%now%is_grline)


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
        call calc_gradient_ac_ice(tpo%now%dzsdx,tpo%now%dzsdy,tpo%now%z_srf,tpo%now%H_ice,tpo%par%dx,tpo%par%margin2nd)
        call calc_gradient_ac_ice(tpo%now%dHicedx,tpo%now%dHicedy,tpo%now%H_ice,tpo%now%H_ice,tpo%par%dx,tpo%par%margin2nd)
        
        ! Calculate effective pressure [bar = 1e-5 Pa == 1e-5 kg m^-1 s^-2]
        tpo%now%N_eff = calc_effective_pressure(tpo%now%H_ice,bnd%z_bed,bnd%z_sl,bnd%H_w,p=tpo%par%neff_p)

!         ! Ensure effective pressure is never zero 
!         N_eff_min = minval(tpo%now%N_eff,mask=tpo%now%N_eff .ne. 0.0)
!         where(tpo%now%N_eff .lt. N_eff_min) tpo%now%N_eff = N_eff_min
        
        ! Calculate distance to ice margin 
        !tpo%now%dist_margin = distance_to_margin(tpo%now%H_ice,tpo%par%dx)

        ! Calculate distance to grounding line 
        !tpo%now%dist_grline = distance_to_grline(tpo%now%is_grline,tpo%now%f_grnd,tpo%par%dx)

        ! Advance ytopo timestep 
        tpo%par%time = time

        if (yelmo_write_log) then 

            if (count(tpo%now%H_ice.gt.0.0) .gt. 0) then 
                write(*,"(a,f14.4,f10.4,f10.2)") "calc_ytopo::  time = ", tpo%par%time, dt, &
                    sum(tpo%now%H_ice,mask=tpo%now%H_ice.gt.0.0)/real(count(tpo%now%H_ice.gt.0.0))
            else 
                write(*,"(a,f14.4,f10.4,f10.2)") "calc_ytopo::  time = ", tpo%par%time, dt, 0.0
            end if 

        end if 

        return 

    end subroutine calc_ytopo

    subroutine ytopo_load_H_ice(tpo,nml_path,nml_group,domain,grid_name)
        ! Update the topography of the yelmo domain 

        implicit none

        type(ytopo_class), intent(INOUT) :: tpo
        character(len=*), intent(IN)     :: nml_path, nml_group
        character(len=*), intent(IN)     :: domain, grid_name 

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

        ! Reset H_ice to zero if initialization parameter is set to "icefree"
        if (trim(tpo%par%init) == "icefree") tpo%now%H_ice = 0.d0

        write(*,*) "ytopo_load_H_ice:: range(H_ice):  ", minval(tpo%now%H_ice),  maxval(tpo%now%H_ice)

        return 

    end subroutine ytopo_load_H_ice 

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
        call nml_read(filename,"ytopo","method",            par%method,            init=init_pars)
        call nml_read(filename,"ytopo","init",              par%init,              init=init_pars)
        call nml_read(filename,"ytopo","solver",            par%solver,            init=init_pars)
        call nml_read(filename,"ytopo","margin2nd",         par%margin2nd,         init=init_pars)
        call nml_read(filename,"ytopo","surf_gl_method",    par%surf_gl_method,    init=init_pars)
        call nml_read(filename,"ytopo","calv_method",       par%calv_method,       init=init_pars)
        
        call nml_read(filename,"ytopo","use_bmb",           par%use_bmb,          init=init_pars)
        call nml_read(filename,"ytopo","use_calv_subgrid",  par%use_calv_subgrid, init=init_pars)
        call nml_read(filename,"ytopo","ocean_kill",        par%ocean_kill,       init=init_pars)
        call nml_read(filename,"ytopo","grline_fixed",      par%grline_fixed,     init=init_pars)
        call nml_read(filename,"ytopo","topo_fixed",        par%topo_fixed,       init=init_pars)
        call nml_read(filename,"ytopo","topo_relax_dt",     par%topo_relax_dt,    init=init_pars)
        call nml_read(filename,"ytopo","topo_fixed_dt",     par%topo_fixed_dt,    init=init_pars)
        call nml_read(filename,"ytopo","calv_dt",           par%calv_dt,          init=init_pars)
        call nml_read(filename,"ytopo","H_calv",            par%H_calv,           init=init_pars)
        call nml_read(filename,"ytopo","H_min",             par%H_min,            init=init_pars)
        call nml_read(filename,"ytopo","gl_sep",            par%gl_sep,           init=init_pars)
        call nml_read(filename,"ytopo","gl_sep_nx",         par%gl_sep_nx,        init=init_pars)
        call nml_read(filename,"ytopo","diffuse_bmb_shlf",  par%diffuse_bmb_shlf, init=init_pars)
                
        call nml_read(filename,"ytopo","neff_p",            par%neff_p,           init=init_pars)
        
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
        par%time_calv = par%time 

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
        allocate(now%calv(nx,ny))

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

        allocate(now%N_eff(nx,ny))

        now%H_ice      = 0.0 
        now%z_srf      = 0.0  
        now%dzsrfdt    = 0.0 
        now%dHicedt    = 0.0
        now%bmb        = 0.0  
        now%mb_applied = 0.0 
        now%calv       = 0.0
        now%dzsdx      = 0.0 
        now%dzsdy      = 0.0 
        now%dHicedx    = 0.0 
        now%dHicedy    = 0.0
        now%H_grnd     = 0.0  
        now%f_grnd     = 0.0  
        now%f_grnd_acx = 0.0  
        now%f_grnd_acy = 0.0  
        now%f_ice      = 0.0  
        now%dist_margin = 0.0
        now%dist_grline = 0.0
        now%mask_bed   = 0.0 
        now%is_grline  = .FALSE. 
        now%is_grz     = .FALSE. 
        now%N_eff      = 0.0 

        return 
    end subroutine ytopo_alloc 

    subroutine ytopo_dealloc(now)

        implicit none 

        type(ytopo_state_class), intent(INOUT) :: now

        if (allocated(now%H_ice))      deallocate(now%H_ice)
        if (allocated(now%z_srf))      deallocate(now%z_srf)
        
        if (allocated(now%dzsrfdt))    deallocate(now%dzsrfdt)
        if (allocated(now%dHicedt))    deallocate(now%dHicedt)
        if (allocated(now%bmb))        deallocate(now%bmb)
        if (allocated(now%mb_applied)) deallocate(now%mb_applied)
        if (allocated(now%calv))       deallocate(now%calv)

        if (allocated(now%dzsdx))      deallocate(now%dzsdx)
        if (allocated(now%dzsdy))      deallocate(now%dzsdy)
        if (allocated(now%dHicedx))    deallocate(now%dHicedx)
        if (allocated(now%dHicedy))    deallocate(now%dHicedy)
        
        if (allocated(now%H_grnd))     deallocate(now%H_grnd)

        ! Masks
        if (allocated(now%f_grnd))     deallocate(now%f_grnd)
        if (allocated(now%f_grnd_acx)) deallocate(now%f_grnd_acx)
        if (allocated(now%f_grnd_acy)) deallocate(now%f_grnd_acy)
        
        if (allocated(now%f_ice))      deallocate(now%f_ice)

        if (allocated(now%dist_margin)) deallocate(now%dist_margin)
        if (allocated(now%dist_grline)) deallocate(now%dist_grline)
        
        if (allocated(now%mask_bed))   deallocate(now%mask_bed)
        if (allocated(now%is_grline))  deallocate(now%is_grline)
        if (allocated(now%is_grz))     deallocate(now%is_grz)

        if (allocated(now%N_eff))      deallocate(now%N_eff)
        
        return 

    end subroutine ytopo_dealloc 

    ! ============================================================
    !
    ! Calculations
    !
    ! ============================================================

    function calc_ice_fraction(H_ice,f_grnd) result(f_ice)

        implicit none 

        real(prec), intent(IN) :: H_ice(:,:) 
        real(prec), intent(IN) :: f_grnd(:,:)  
        real(prec) :: f_ice(size(H_ice,1),size(H_ice,2))

        ! Local variables 
        integer :: i, j, nx, ny 
        real(prec) :: dx, dy 
        real(prec) :: H_ref(4)
        logical :: mask_ref(4)
        real(prec) :: H_frac 

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        ! Set dx/dy to one, since the volume is only used in a relative sense 
        dx = 1.0 
        dy = 1.0 

        ! Initially set fraction to one everywhere there is ice 
        ! and zero everywhere there is no ice 
        where (f_grnd .eq. 0.0 .and. H_ice .gt. 0.0)
            ! Floating ice 
            f_ice = 1.0
        else where ( f_grnd .gt. 0.0 .and. H_ice .gt. 0.0 )
            ! Grounded ice 
            f_ice = 1.0 
        elsewhere
            ! No ice  
            f_ice = 0.0 
        end where

if (.FALSE.) then 
        ! For floating points at the border of ice and no ice,
        ! determine the fraction of ice in the cell.
        do j = 2, ny-1
        do i = 2, nx-1 

            if (f_grnd(i,j) .eq. 0.0 .and. f_ice(i,j) .gt. 0.0 .and. &
                count([f_ice(i-1,j),f_ice(i+1,j),f_ice(i,j-1),f_ice(i,j+1)].eq.0) .gt. 0) then 
                ! This point is at the calving front 

                H_ref    = [H_ice(i-1,j),H_ice(i+1,j),H_ice(i,j-1),H_ice(i,j+1)]
                mask_ref = ([f_grnd(i-1,j),f_grnd(i+1,j),f_grnd(i,j-1),f_grnd(i,j+1)] .eq. 0.0 &
                            .and. H_ref .gt. 0.0) .or. &
                           ([f_grnd(i-1,j),f_grnd(i+1,j),f_grnd(i,j-1),f_grnd(i,j+1)] .gt. 0.0 &
                            .and. H_ref .gt. 0.0)

                if (count(mask_ref) .gt. 0) then 
                    ! Neighbors with ice should generally be found, but put this check just in case

                    ! Determine height to give to partially filled cell as average of neighbors
                    H_frac = sum(H_ref,mask=mask_ref)/real(count(mask_ref))

                    ! Determine the cell ice fraction
                    ! Note: fraction is determined as a ratio of 
                    ! thicknesses, derived from volume conservation 
                    ! vol = H_ice*dx*dy = H_frac*area_frac 
                    ! f_ice = area_frac / (dx*dy)
                    ! f_ice = H_ice/H_frac 
                    f_ice(i,j) = min( H_ice(i,j) / H_frac, 1.0 ) 

                end if

            end if  

        end do 
        end do 

end if 
        ! Note: for now, do not treat fraction of grounded ice,
        ! however this should be considered in the future.

        return 

    end function calc_ice_fraction

    elemental subroutine calc_z_srf(z_srf,H_ice,H_grnd,z_bed,z_sl)
        ! Calculate surface elevation

        implicit none 

        real(prec), intent(INOUT) :: z_srf
        real(prec), intent(IN)    :: H_ice
        real(prec), intent(IN)    :: H_grnd
        real(prec), intent(IN)    :: z_bed
        real(prec), intent(IN)    :: z_sl

        ! Local variables 
        real(prec) :: rho_ice_sw

        rho_ice_sw = rho_ice/rho_sw ! Ratio of density of ice to seawater [--]
        
        ! Update the surface elevation based on z_bed, H_ice and overburden ice thickness 
        if (H_grnd .gt. 0.0) then 
            ! Grounded ice or ice-free land

            z_srf = z_bed + H_ice 

        else
            ! Floating ice or open ocean

            z_srf = z_sl + (1.0-rho_ice_sw)*H_ice

        end if 
        
        return 

    end subroutine calc_z_srf 

    elemental subroutine calc_z_srf_max(z_srf,H_ice,z_bed,z_sl)
        ! Calculate surface elevation
        ! Adapted from Pattyn (2017), Eq. 1
        
        implicit none 

        real(prec), intent(INOUT) :: z_srf
        real(prec), intent(IN)    :: H_ice
        real(prec), intent(IN)    :: z_bed
        real(prec), intent(IN)    :: z_sl

        ! Local variables
        real(prec) :: rho_ice_sw

        rho_ice_sw = rho_ice/rho_sw ! Ratio of density of ice to seawater [--]
        
        z_srf = max(z_bed + H_ice, z_sl + (1.0-rho_ice_sw)*H_ice)
        
        return 

    end subroutine calc_z_srf_max 

    subroutine calc_z_srf_subgrid_area(z_srf,f_grnd,H_ice,z_bed,z_sl,gl_sep_nx)
        ! Interpolate variables at grounding line to subgrid level to 
        ! calculate the average z_srf value for the aa-node cell

        implicit none
        
        real(prec), intent(OUT) :: z_srf(:,:)       ! aa-nodes 
        real(prec), intent(IN)  :: f_grnd(:,:)      ! aa-nodes
        real(prec), intent(IN)  :: H_ice(:,:)
        real(prec), intent(IN)  :: z_bed(:,:)
        real(prec), intent(IN)  :: z_sl(:,:)
        integer,    intent(IN)  :: gl_sep_nx        ! Number of interpolation points per side (nx*nx)

        ! Local variables
        integer    :: i, j, nx, ny
        real(prec) :: v1, v2, v3, v4 
        integer    :: i1, i2, j1, j2 
        real(prec) :: rho_ice_sw

        real(prec), allocatable :: z_srf_int(:,:) 
        real(prec), allocatable :: H_ice_int(:,:) 
        real(prec), allocatable :: z_bed_int(:,:) 
        real(prec), allocatable :: z_sl_int(:,:) 
        real(prec), allocatable :: H_grnd_int(:,:) 

        rho_ice_sw = rho_ice/rho_sw ! Ratio of density of ice to seawater [--]

        nx = size(z_srf,1)
        ny = size(z_srf,2) 

        ! Allocate the subgrid arrays 
        allocate(z_srf_int(gl_sep_nx,gl_sep_nx))
        allocate(H_ice_int(gl_sep_nx,gl_sep_nx))
        allocate(z_bed_int(gl_sep_nx,gl_sep_nx))
        allocate(z_sl_int(gl_sep_nx,gl_sep_nx))
        allocate(H_grnd_int(gl_sep_nx,gl_sep_nx))
        
        ! Calculate the surface elevation based on whole grid values,
        ! except at the grounding line which is treated with subgrid interpolations. 
        do j = 1, ny 
        do i = 1, nx

            if (f_grnd(i,j) .eq. 1.0) then 
                ! Fully grounded grid point 

                z_srf(i,j) = z_bed(i,j) + H_ice(i,j)

            else if (f_grnd(i,j) .eq. 0.0) then 
                ! Fully floating grid point 

                z_srf(i,j) = z_sl(i,j) + (1.0-rho_ice_sw)*H_ice(i,j)

            else 
                ! Partially grounded and floating (ie, grounding line) 
                ! Perform subgrid calculations 

                i1 = max(i-1,1) 
                i2 = min(i+1,nx) 
                j1 = max(j-1,1) 
                j2 = min(j+1,ny) 

                ! Calculate values at corners (ab-nodes) and interpolate
                
                ! == H_ice == 
                v1 = 0.25_prec*(H_ice(i,j) + H_ice(i2,j) + H_ice(i2,j2) + H_ice(i,j2))
                v2 = 0.25_prec*(H_ice(i,j) + H_ice(i1,j) + H_ice(i1,j2) + H_ice(i,j2))
                v3 = 0.25_prec*(H_ice(i,j) + H_ice(i1,j) + H_ice(i1,j1) + H_ice(i,j1))
                v4 = 0.25_prec*(H_ice(i,j) + H_ice(i2,j) + H_ice(i2,j1) + H_ice(i,j1))
                call calc_subgrid_array(H_ice_int,v1,v2,v3,v4,gl_sep_nx)
                
                ! == z_bed == 
                v1 = 0.25_prec*(z_bed(i,j) + z_bed(i2,j) + z_bed(i2,j2) + z_bed(i,j2))
                v2 = 0.25_prec*(z_bed(i,j) + z_bed(i1,j) + z_bed(i1,j2) + z_bed(i,j2))
                v3 = 0.25_prec*(z_bed(i,j) + z_bed(i1,j) + z_bed(i1,j1) + z_bed(i,j1))
                v4 = 0.25_prec*(z_bed(i,j) + z_bed(i2,j) + z_bed(i2,j1) + z_bed(i,j1))
                call calc_subgrid_array(z_bed_int,v1,v2,v3,v4,gl_sep_nx)
                
                ! == z_sl == 
                v1 = 0.25_prec*(z_sl(i,j) + z_sl(i2,j) + z_sl(i2,j2) + z_sl(i,j2))
                v2 = 0.25_prec*(z_sl(i,j) + z_sl(i1,j) + z_sl(i1,j2) + z_sl(i,j2))
                v3 = 0.25_prec*(z_sl(i,j) + z_sl(i1,j) + z_sl(i1,j1) + z_sl(i,j1))
                v4 = 0.25_prec*(z_sl(i,j) + z_sl(i2,j) + z_sl(i2,j1) + z_sl(i,j1))
                call calc_subgrid_array(z_sl_int,v1,v2,v3,v4,gl_sep_nx)
                
                
                ! Now calculate H_grnd and z_srf for each subgrid point 
                call calc_H_grnd(H_grnd_int,H_ice_int,z_bed_int,z_sl_int)

                where (H_grnd_int .gt. 0.0)
                    ! Fully grounded ice or ice-free land 
                    z_srf_int = z_bed_int + H_ice_int 

                elsewhere 
                    ! Fully floating ice or open ocean
                    z_srf_int = z_sl_int + (1.0-rho_ice_sw)*H_ice_int 

                end where 

                ! Calculate full grid z_srf value as the mean of subgrid values 
                z_srf(i,j) = sum(z_srf_int) / real(gl_sep_nx*gl_sep_nx,prec)

            end if 

        end do 
        end do 

        return
        
    end subroutine calc_z_srf_subgrid_area
    
    elemental subroutine calc_H_grnd(H_grnd,H_ice,z_bed,z_sl)
        ! Calculate ice thickness overburden, H_grnd
        ! When H_grnd >= 0, grounded, when H_grnd < 0, floating 
        ! Also calculate rate of change for diagnostic related to grounding line 

        implicit none 

        real(prec), intent(INOUT) :: H_grnd
        real(prec), intent(IN)    :: H_ice
        real(prec), intent(IN)    :: z_bed
        real(prec), intent(IN)    :: z_sl 

        ! Local variables   
        real(prec) :: rho_sw_ice 

        rho_sw_ice = rho_sw/rho_ice ! Ratio of density of seawater to ice [--]
        
        ! Calculate new H_grnd (ice thickness overburden)
        H_grnd = H_ice - rho_sw_ice*(z_sl-z_bed)

        return 

    end subroutine calc_H_grnd

    subroutine calc_f_grnd_subgrid_area(f_grnd,f_grnd_acx,f_grnd_acy,H_grnd,gl_sep_nx)
        ! Use H_grnd to determined grounded area fraction of grid point.

        implicit none
        
        real(prec), intent(OUT) :: f_grnd(:,:)      ! aa-nodes 
        real(prec), intent(OUT) :: f_grnd_acx(:,:)  ! ac-nodes
        real(prec), intent(OUT) :: f_grnd_acy(:,:)  ! ac-nodes
        real(prec), intent(IN)  :: H_grnd(:,:)      ! aa-nodes
        integer,    intent(IN)  :: gl_sep_nx        ! Number of interpolation points per side (nx*nx)

        ! Local variables
        integer    :: i, j, nx, ny
        real(prec) :: Hg_1, Hg_2, Hg_3, Hg_4, Hg_mid  
        integer    :: i1, i2, j1, j2 

        !integer, parameter :: nx_interp = 15

        nx = size(H_grnd,1)
        ny = size(H_grnd,2) 

        ! First binary estimate of f_grnd based on aa-nodes
        f_grnd = 1.0
        where (H_grnd < 0.0) f_grnd = 0.0
        
        ! Find grounding line cells and determine fraction 
        do j = 1, ny 
        do i = 1, nx

            i1 = max(i-1,1) 
            i2 = min(i+1,nx) 
            j1 = max(j-1,1) 
            j2 = min(j+1,ny) 

            ! Calculate Hg at corners (ab-nodes)
            Hg_1 = 0.25_prec*(H_grnd(i,j) + H_grnd(i2,j) + H_grnd(i2,j2) + H_grnd(i,j2))
            Hg_2 = 0.25_prec*(H_grnd(i,j) + H_grnd(i1,j) + H_grnd(i1,j2) + H_grnd(i,j2))
            Hg_3 = 0.25_prec*(H_grnd(i,j) + H_grnd(i1,j) + H_grnd(i1,j1) + H_grnd(i,j1))
            Hg_4 = 0.25_prec*(H_grnd(i,j) + H_grnd(i2,j) + H_grnd(i2,j1) + H_grnd(i,j1))
            
            if (max(Hg_1,Hg_2,Hg_3,Hg_4) .ge. 0.0 .and. min(Hg_1,Hg_2,Hg_3,Hg_4) .lt. 0.0) then 
                ! Point contains grounding line, get grounded area  
                
                call calc_grounded_fraction_cell(f_grnd(i,j),Hg_1,Hg_2,Hg_3,Hg_4,gl_sep_nx)

            end if 

        end do 
        end do 

if (.TRUE.) then 
        ! Linear average to ac-nodes 

        ! acx-nodes 
        do j = 1, ny 
        do i = 1, nx-1
            f_grnd_acx(i,j) = 0.5_prec*(f_grnd(i,j) + f_grnd(i+1,j))
        end do 
        end do
        f_grnd_acx(nx,:) = f_grnd_acx(nx-1,:) 

        ! acy-nodes 
        do j = 1, ny-1 
        do i = 1, nx
            f_grnd_acy(i,j) = 0.5_prec*(f_grnd(i,j) + f_grnd(i,j+1))
        end do 
        end do
        f_grnd_acy(:,ny) = f_grnd_acy(:,ny-1) 

else 
        ! Subgrid area calcs centered on ac-nodes directly

        ! acx-nodes 
        do j = 1, ny 
        do i = 1, nx
            i1 = max(i-1,1) 
            i2 = min(i+1,nx) 
            j1 = max(j-1,1) 
            j2 = min(j+1,ny) 

            ! Calculate Hg at corners (acy-nodes)
            Hg_1 = 0.5_prec*(H_grnd(i2,j) + H_grnd(i2,j2))
            Hg_2 = 0.5_prec*(H_grnd(i,j)  + H_grnd(i,j2))
            Hg_3 = 0.5_prec*(H_grnd(i,j)  + H_grnd(i,j1))
            Hg_4 = 0.5_prec*(H_grnd(i2,j) + H_grnd(i2,j1))
            
            if (max(Hg_1,Hg_2,Hg_3,Hg_4) .ge. 0.0 .and. min(Hg_1,Hg_2,Hg_3,Hg_4) .lt. 0.0) then 
                ! Point contains grounding line, get grounded area  
                
                call calc_grounded_fraction_cell(f_grnd_acx(i,j),Hg_1,Hg_2,Hg_3,Hg_4,gl_sep_nx)

            else 
                ! Purely grounded or floating point 

                Hg_mid = 0.5_prec*(H_grnd(i,j)  + H_grnd(i2,j))

                if (Hg_mid .lt. 0.0) then 
                    f_grnd_acx(i,j) = 0.0_prec 
                else 
                    f_grnd_acx(i,j) = 1.0_prec 
                end if 

            end if 

        end do 
        end do

        ! acy-nodes 
        do j = 1, ny 
        do i = 1, nx
            i1 = max(i-1,1) 
            i2 = min(i+1,nx) 
            j1 = max(j-1,1) 
            j2 = min(j+1,ny) 

            ! Calculate Hg at corners (acx-nodes)
            Hg_1 = 0.5_prec*(H_grnd(i,j2)  + H_grnd(i2,j2))
            Hg_2 = 0.5_prec*(H_grnd(i1,j2) + H_grnd(i,j2))
            Hg_3 = 0.5_prec*(H_grnd(i1,j)  + H_grnd(i,j))
            Hg_4 = 0.5_prec*(H_grnd(i2,j)  + H_grnd(i,j))
            
            if (max(Hg_1,Hg_2,Hg_3,Hg_4) .ge. 0.0 .and. min(Hg_1,Hg_2,Hg_3,Hg_4) .lt. 0.0) then 
                ! Point contains grounding line, get grounded area  
                
                call calc_grounded_fraction_cell(f_grnd_acy(i,j),Hg_1,Hg_2,Hg_3,Hg_4,gl_sep_nx)

            else 
                ! Purely grounded or floating point 
                
                Hg_mid = 0.5_prec*(H_grnd(i,j)  + H_grnd(i,j2))
                
                if (Hg_mid .lt. 0.0) then 
                    f_grnd_acy(i,j) = 0.0_prec 
                else 
                    f_grnd_acy(i,j) = 1.0_prec 
                end if 
                
            end if 

        end do 
        end do

end if 

        return
        
    end subroutine calc_f_grnd_subgrid_area
    
    subroutine calc_f_grnd_subgrid_ac_linear(f_grnd,f_grnd_x,f_grnd_y,H_grnd)
        ! Calculate the grounded fraction of a cell in the x- and y-directions
        ! at the ac nodes
        !
        ! Given point 1 is grounded and point 2 is floating, we are looking for
        ! when H_grnd=0. The equation along an x-axis between point 1 (grounded point) 
        ! and point 2 (floating point) for H_grnd is:
        ! H_grnd = H_grnd_1 + f*(H_grnd_2-H_grnd_1)

        ! To find fraction f for when H_grnd == 0, rearrange equation:
        ! 0 = H_grnd_1 + f*(H_grnd_2-H_grnd_1)
        ! f = -H_grnd_1 / (H_grnd_2-H_grnd_1)
        ! where f is the distance along the 0:1 axis between the first and second points. 
        
        implicit none 

        real(prec), intent(OUT) :: f_grnd(:,:)
        real(prec), intent(OUT) :: f_grnd_x(:,:)
        real(prec), intent(OUT) :: f_grnd_y(:,:)
        real(prec), intent(IN)  :: H_grnd(:,:)

        ! Local variables  
        integer :: i, j, nx, ny 
        real(prec) :: H_grnd_1, H_grnd_2

        nx = size(f_grnd,1)
        ny = size(f_grnd,2)

        ! Central Aa node
        f_grnd = 1.0
        where (H_grnd < 0.0) f_grnd = 0.0
        
        ! x-direction, Ac node
        f_grnd_x = 1.0
        do j = 1, ny 
        do i = 1, nx-1 

            if (H_grnd(i,j) .gt. 0.0 .and. H_grnd(i+1,j) .le. 0.0) then 
                ! Point is grounded, neighbor is floating 

                H_grnd_1 = H_grnd(i,j) 
                H_grnd_2 = H_grnd(i+1,j) 

                ! Calculate fraction 
                f_grnd_x(i,j) = -H_grnd_1 / (H_grnd_2 - H_grnd_1)

            else if (H_grnd(i,j) .le. 0.0 .and. H_grnd(i+1,j) .gt. 0.0) then 
                ! Point is floating, neighbor is grounded 

                H_grnd_1 = H_grnd(i+1,j) 
                H_grnd_2 = H_grnd(i,j) 

                ! Calculate fraction 
                f_grnd_x(i,j) = -H_grnd_1 / (H_grnd_2 - H_grnd_1)

            else if (H_grnd(i,j) .le. 0.0 .and. H_grnd(i+1,j) .le. 0.0) then 
                ! Point is floating, neighbor is floating
                f_grnd_x(i,j) = 0.0 

            else 
                ! Point is grounded, neighbor is grounded
                f_grnd_x(i,j) = 1.0 

            end if 

        end do 
        end do 

        ! y-direction, Ac node
        f_grnd_y = 1.0
        do j = 1, ny-1 
        do i = 1, nx 

            if (H_grnd(i,j) .gt. 0.0 .and. H_grnd(i,j+1) .le. 0.0) then 
                ! Point is grounded, neighbor is floating 

                H_grnd_1 = H_grnd(i,j) 
                H_grnd_2 = H_grnd(i,j+1) 

                ! Calculate fraction 
                f_grnd_y(i,j) = -H_grnd_1 / (H_grnd_2 - H_grnd_1)

            else if (H_grnd(i,j) .le. 0.0 .and. H_grnd(i,j+1) .gt. 0.0) then 
                ! Point is floating, neighbor is grounded 

                H_grnd_1 = H_grnd(i,j+1) 
                H_grnd_2 = H_grnd(i,j) 

                ! Calculate fraction 
                f_grnd_y(i,j) = -H_grnd_1 / (H_grnd_2 - H_grnd_1)
                
            else if (H_grnd(i,j) .le. 0.0 .and. H_grnd(i,j+1) .le. 0.0) then 
                ! Point is floating, neighbor is floating
                f_grnd_y(i,j) = 0.0 

            else 
                ! Point is grounded, neighbor is grounded
                f_grnd_y(i,j) = 1.0 

            end if 

        end do 
        end do 

        ! Set boundary points equal to neighbor for aesthetics 
        f_grnd_x(nx,:) = f_grnd_x(nx-1,:) 
        f_grnd_y(:,ny) = f_grnd_y(:,ny-1) 
        
        return 

    end subroutine calc_f_grnd_subgrid_ac_linear
    
    subroutine filter_f_grnd(f_grnd)
        ! Remove isolated floating points inside of grounded ice 

        implicit none 

        real(prec), intent(INOUT) :: f_grnd(:,:)
        
        ! Local variables 
        integer :: i, j, nx, ny 
        real(prec), allocatable :: f_grnd_old(:,:) 
        logical, allocatable :: is_grline(:,:) 

        nx = size(f_grnd,1)
        ny = size(f_grnd,2)

        allocate(f_grnd_old(nx,ny))
        allocate(is_grline(nx,ny))

        ! == LAKES ==
        f_grnd_old = f_grnd
        do i = 2, nx-1 
        do j = 2, ny-1 

            if (f_grnd_old(i,j) .lt. 1.0) then 
                ! Check partially floating points 

                if (f_grnd_old(i-1,j)+f_grnd_old(i+1,j)+f_grnd_old(i,j-1)+f_grnd_old(i,j+1) .eq. 4.0) then 
                    ! Point is surrounded by fully grounded ice (ie, a lake),
                    ! so treat it as grounded ice for now 
                    f_grnd(i,j) = 1.0 

                end if 

            end if 

        end do 
        end do 

        ! == GRLINE ISLANDS == 
        f_grnd_old = f_grnd 

        ! Calculate intermediate grounding line estimate 
        call calc_grline(is_grline,f_grnd)

        do i = 2, nx-1 
        do j = 2, ny-1 

            if (is_grline(i,j)) then 
                ! Grounding line defined point 

                if ( (f_grnd_old(i-1,j) .eq. 0.0 .or. is_grline(i-1,j)) .and. &
                     (f_grnd_old(i+1,j) .eq. 0.0 .or. is_grline(i+1,j)) .and. &
                     (f_grnd_old(i,j-1) .eq. 0.0 .or. is_grline(i,j-1)) .and. &
                     (f_grnd_old(i,j+1) .eq. 0.0 .or. is_grline(i,j+1)) ) then   
                    ! Point is surrounded by floating or grounding-line ice,
                    ! so treat it as floating ice for now 
                    f_grnd(i,j) = 0.0 

                end if 

            end if 

        end do 
        end do 

        return 

    end subroutine filter_f_grnd

    subroutine calc_grline(is_grline,f_grnd)

        implicit none 

        logical, intent(OUT) :: is_grline(:,:) 
        real(prec), intent(IN)  :: f_grnd(:,:)

        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: im1, ip1, jm1, jp1 

        nx = size(is_grline,1)
        ny = size(is_grline,2)

        is_grline = .FALSE. 
 
        do j = 1, ny 
        do i = 1, nx

            im1 = max(1, i-1)
            ip1 = min(nx,i+1)
            
            jm1 = max(1, j-1)
            jp1 = min(ny,j+1)

            ! Grounded point or partially floating point with floating neighbors
            if (f_grnd(i,j) .gt. 0.0 .and. &
                (f_grnd(im1,j) .eq. 0.0 .or. f_grnd(ip1,j) .eq. 0.0 .or. &
                 f_grnd(i,jm1) .eq. 0.0 .or. f_grnd(i,jp1) .eq. 0.0) ) then 
                
                is_grline(i,j) = .TRUE. 

            end if 

        end do 
        end do 

        return 

    end subroutine calc_grline

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

    subroutine calc_bmb_total(bmb,bmb_grnd,bmb_shlf,f_grnd,diffuse_bmb_shlf)

        implicit none 

        real(prec), intent(OUT) :: bmb(:,:) 
        real(prec), intent(IN)  :: bmb_grnd(:,:) 
        real(prec), intent(IN)  :: bmb_shlf(:,:) 
        real(prec), intent(IN)  :: f_grnd(:,:) 
        logical,    intent(IN)  :: diffuse_bmb_shlf 

        ! Local variables
        integer    :: i, j, nx, ny
        integer    :: n_float  
        real(prec) :: bmb_shlf_now 

        nx = size(bmb,1)
        ny = size(bmb,2) 

        ! Combine floating and grounded parts into one field =========================
        
        ! Weighted average
        !bmb = f_grnd*bmb_grnd + (1.0-f_grnd)*bmb_shlf

        ! Simply add the two fields for now to avoid complications with f_grnd
        bmb = bmb_grnd + bmb_shlf 


        if (diffuse_bmb_shlf) then 
            ! Allow marine melt (bmb_shlf) to permeate inland at the grounding line,
            ! to induce more effective retreat in warm periods 

            do j = 1, ny
            do i = 1, nx

                if (f_grnd(i,j) .eq. 1.0) then 
                    ! Grounded point, look for floating neighbors 

                    if (.FALSE.) then
                        ! 9-neighbours method

                        n_float = count(f_grnd(i-1:i+1,j-1:j+1) .lt. 1.0)

                        if (n_float .gt. 0) then 
                            ! bmb_shelf is the mean of the neighbours
                            bmb_shlf_now = sum(bmb_shlf(i-1:i+1,j-1:j+1),mask=f_grnd(i-1:i+1,j-1:j+1) .lt. 1.0) / real(n_float,prec)
                            bmb(i,j)     = (1.0-n_float/9.0)*bmb_grnd(i,j) + (n_float/9.0)*bmb_shlf_now

                        end if

                    else 
                        ! 4-neighbors method 

                        n_float = count([f_grnd(i-1,j),f_grnd(i+1,j),f_grnd(i,j-1),f_grnd(i,j+1)].lt. 1.0)

                        if (n_float .gt. 0) then
                            ! Floating points exist 
                            bmb_shlf_now = sum([bmb_shlf(i-1,j),bmb_shlf(i+1,j),bmb_shlf(i,j-1),bmb_shlf(i,j+1)], &
                                            mask=[f_grnd(i-1,j),f_grnd(i+1,j),f_grnd(i,j-1),f_grnd(i,j+1)] .lt. 1.0) / real(n_float,prec)
                            bmb(i,j)     = (1.0-n_float/5.0)*bmb_grnd(i,j) + (n_float/5.0)*bmb_shlf_now
                        end if

                    end if


                end if

            end do
            end do

        end if 

        return 

    end subroutine calc_bmb_total 

    subroutine calc_subgrid_array(vint,v1,v2,v3,v4,nx)
        ! Given the four corners of a cell in quadrants 1,2,3,4,
        ! calculate the subgrid values via linear interpolation

        implicit none 

        real(prec), intent(OUT) :: vint(:,:)  
        real(prec), intent(IN)  :: v1,v2,v3,v4
        integer,    intent(IN)  :: nx                    ! Number of interpolation points 

        ! Local variables 
        integer :: i, j 
        real(prec) :: x(nx), y(nx) 

        ! Populate x,y axes for interpolation points (between 0 and 1)
        do i = 1, nx 
            x(i) = 0.0 + real(i-1)/real(nx-1)
        end do 
        y = x 
        
        ! Calculate interpolated value      
        vint = 0.0 
        do i = 1, nx 
        do j = 1, nx 

            vint(i,j) = interp_bilin_pt(v1,v2,v3,v4,x(i),y(j))

        end do 
        end do 

        return 

    end subroutine calc_subgrid_array

    subroutine calc_grounded_fraction_cell(f_g,Hg_1,Hg_2,Hg_3,Hg_4,nx)
        ! Given the four corners of a cell in quadrants 1,2,3,4,
        ! calculate the grounded fraction (ie area with Hg>0)

        implicit none 

        real(prec), intent(OUT) :: f_g 
        real(prec), intent(IN)  :: Hg_1,Hg_2,Hg_3,Hg_4
        integer,    intent(IN)  :: nx                    ! Number of interpolation points 

        ! Local variables 
        integer :: i, j 
        real(prec) :: x(nx), y(nx) 
        real(prec) :: Hg_int(nx,nx)  

        ! Populate x,y axes for interpolation points (between 0 and 1)
        do i = 1, nx 
            x(i) = 0.0 + real(i-1)/real(nx-1)
        end do 
        y = x 


        ! Calculate interpolation heights
        Hg_int = 0.0 
        do i = 1, nx 
        do j = 1, nx 

            Hg_int(i,j) = interp_bilin_pt(Hg_1,Hg_2,Hg_3,Hg_4,x(i),y(j))

        end do 
        end do 

        ! Calculate weighted fraction (assume all points have equal weight)
        f_g = real(count(Hg_int .ge. 0.0),prec) / real(nx*nx,prec)

        return 

    end subroutine calc_grounded_fraction_cell

    function interp_bilin_pt(z1,z2,z3,z4,xout,yout) result(zout)
        ! Interpolate a point given four neighbors at corners of square (0:1,0:1)
        ! z2    z1
        !    x,y
        ! z3    z4 
        ! 

        implicit none 

        real(prec), intent(IN) :: z1, z2, z3, z4 
        real(prec), intent(IN) :: xout, yout 
        real(prec) :: zout 

        ! Local variables 
        real(prec) :: x0, x1, y0, y1 
        real(prec) :: alpha1, alpha2, p0, p1 

        x0 = 0.0 
        x1 = 1.0 
        y0 = 0.0 
        y1 = 1.0 

        alpha1  = (xout - x0) / (x1-x0)
        p0      = z3 + alpha1*(z4-z3)
        p1      = z2 + alpha1*(z1-z2)
            
        alpha2  = (yout - y0) / (y1-y0)
        zout    = p0 + alpha2*(p1-p0)

        return 

    end function interp_bilin_pt

    function distance_to_margin(H_ice,dx) result(dist)

        implicit none 

        real(prec), intent(IN) :: H_ice(:,:) 
        real(prec), intent(IN) :: dx 
        real(prec) :: dist(size(H_ice,1),size(H_ice,2))

        ! Local variables 
        integer :: i, j, nx, ny
        integer :: i1, j1  
        real(prec) :: dx_km
        real(prec), allocatable :: dists(:,:) 

        real(prec), parameter :: dist_max = 1e10 

        nx = size(H_ice,1)
        ny = size(H_ice,2) 

        allocate(dists(nx,ny))

        dx_km = dx*1e-3 

        ! Initially set everything to zero distance to margin 
        dist = 0.0 

        do j = 1, ny 
        do i = 1, nx 

            if (H_ice(i,j) .gt. 0.0) then 
                ! Ice-covered point, check min distance to margin 

                dists = dist_max 

                do j1 = 1, ny 
                do i1 = 1, nx 
                    if (H_ice(i1,j1) .eq. 0.0) then 
                        ! Check distance 
                        dists(i1,j1) = sqrt(((i1-i)*dx_km)**2 + ((j1-j)*dx_km)**2) 
                    end if 
                end do 
                end do 

                ! Get minimum distance 
                dist(i,j) = minval(dists)
                
            end if 

        end do 
        end do 


        return 

    end function distance_to_margin
    
    function distance_to_grline(is_grline,f_grnd,dx) result(dist)

        implicit none 
         
        logical,    intent(IN) :: is_grline(:,:)
        real(prec), intent(IN) :: f_grnd(:,:) 
        real(prec), intent(IN) :: dx 
        real(prec) :: dist(size(is_grline,1),size(is_grline,2))

        ! Local variables 
        integer :: i, j, nx, ny
        integer :: i1, j1  
        real(prec) :: dx_km
        real(prec), allocatable :: dists(:,:) 

        real(prec), parameter :: dist_max = 1e10 

        nx = size(is_grline,1)
        ny = size(is_grline,2) 

        allocate(dists(nx,ny))

        dx_km = dx*1e-3 

        ! Initially set everything to zero distance to margin 
        dist = 0.0 

        if (count(is_grline) .gt. 0) then 
            ! If grounding-line points exist, loop over points checking distances

            do j = 1, ny 
            do i = 1, nx 

                if (.not. is_grline(i,j)) then 
                    ! Not at the grounding line, check min distance to grounding line 

                    dists = dist_max 

                    do j1 = 1, ny 
                    do i1 = 1, nx 
                        if (is_grline(i1,j1)) then 
                            ! Check distance 
                            dists(i1,j1) = sqrt(((i1-i)*dx_km)**2 + ((j1-j)*dx_km)**2) 
                        end if 
                    end do 
                    end do 

                    ! Get minimum distance 
                    dist(i,j) = minval(dists)
                    
                    ! For floating points set distance to negative value 
                    if (f_grnd(i,j) .eq. 0.0) dist(i,j) = -dist(i,j) 

                end if 

            end do 
            end do 

        end if 

        return 

    end function distance_to_grline
    
end module yelmo_topography
