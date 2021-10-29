
module yelmo_material

    use nml 

    use yelmo_defs
    use yelmo_tools, only : calc_vertical_integrated_2D, calc_vertical_integrated_3D, regularize2D
    
    use deformation
    use ice_tracer  

    ! Note: 3D arrays defined such that first index (k=1) == base, and max index (k=nk) == surface 
    
    implicit none
    
    private
    public :: ymat_par_load, ymat_alloc, ymat_dealloc
    public :: calc_ymat 

contains

    subroutine calc_ymat(mat,tpo,dyn,thrm,bnd,time)
        ! Calculate material properties given dynamic and thermodynamic state

        implicit none
        
        type(ymat_class),   intent(INOUT) :: mat
        type(ytopo_class),  intent(IN)    :: tpo 
        type(ydyn_class),   intent(IN)    :: dyn 
        type(ytherm_class), intent(IN)    :: thrm
        type(ybound_class), intent(IN)    :: bnd     
        real(prec),         intent(IN)    :: time    ! Current time (for age tracing)

        ! Local variables
        integer    :: k, nz_aa
        real(prec) :: dt

        real(prec), allocatable :: X_srf(:,:) 
        logical,    allocatable :: mask_tracers(:,:) 

        real(prec), parameter   :: enh_min = 0.1_prec       ! Minimum allowed enhancement factor value (for enh_method="paleo-shear")
        real(prec), parameter   :: enh_max = 10.0_prec      ! Maximum allowed enhancement factor value (for enh_method="paleo-shear")

        nz_aa = mat%par%nz_aa

        ! Allocate temporary arrays 
        allocate(X_srf(mat%par%nx,mat%par%ny))
        allocate(mask_tracers(mat%par%nx,mat%par%ny))

        ! Initialize time if necessary 
        if (mat%par%time .gt. time) then 
            mat%par%time = time
        end if 

        ! Get time step and advance current time 
        dt            = time - mat%par%time 
        mat%par%time  = time 
        
        ! 00. First update ice age if possible
        if (mat%par%calc_age .and. dt .gt. 0.0) then 
            ! Perform calculations of age tracer: dep_time (deposition time)

            ! Set surface boundary condition to current time 
            X_srf = time 

            ! Define limits to where to calculate age tracers
            ! (avoid very fast-flowing ice, as they are not interesting here)
            ! Surface value will be imposed in these places 
            mask_tracers = .TRUE. 
            where (dyn%now%uxy_bar .gt. 500.0_prec) mask_tracers = .FALSE. 

            call calc_tracer_3D(mat%now%dep_time,X_srf,dyn%now%ux,dyn%now%uy,dyn%now%uz, &
                tpo%now%H_ice,tpo%now%bmb,mat%par%zeta_aa,mat%par%zeta_ac,mat%par%tracer_method, &
                mat%par%tracer_impl_kappa,dt,thrm%par%dx,time,mask=mask_tracers)

            ! Calculate isochrones too
            call calc_isochrones(mat%now%depth_iso,mat%now%dep_time,tpo%now%H_ice,mat%par%age_iso, &
                                                                                mat%par%zeta_aa,time)

        end if 

        ! 0. Update strain rate 
        ! Note: this calculation of strain rate is particular to the material module, and 
        ! may differ from the effective strain rate calculated locally in the dynamics module

        ! Calculate the strain rate from the full 3D tensor

        if (.TRUE.) then 
            ! 3D strain rate - standard 

            call calc_strain_rate_3D(mat%now%strn,dyn%now%ux,dyn%now%uy,dyn%now%uz,tpo%now%H_ice, &
                                     tpo%now%f_grnd,mat%par%zeta_aa,mat%par%zeta_ac,mat%par%dx,mat%par%de_max)


            ! And get the vertical average of shear
            mat%now%f_shear_bar = calc_vertical_integrated_2D(mat%now%strn%f_shear,mat%par%zeta_aa) 
            
            ! Get the 2D average of strain rate in case it is needed 
            mat%now%strn2D%dxx = calc_vertical_integrated_2D(mat%now%strn%dxx,mat%par%zeta_aa)
            mat%now%strn2D%dyy = calc_vertical_integrated_2D(mat%now%strn%dyy,mat%par%zeta_aa)
            mat%now%strn2D%dxy = calc_vertical_integrated_2D(mat%now%strn%dxy,mat%par%zeta_aa)
            mat%now%strn2D%de  = calc_vertical_integrated_2D(mat%now%strn%de, mat%par%zeta_aa)
        
        else 
            ! 2D strain rate - for testing only, since it is less accurate than the 3D strain rate

            call calc_strain_rate_2D(mat%now%strn2D,dyn%now%ux_bar,dyn%now%uy_bar,mat%par%dx,mat%par%dx)
            
            do k = 1, nz_aa 
                mat%now%strn%dxx(:,:,k) = mat%now%strn2D%dxx 
                mat%now%strn%dyy(:,:,k) = mat%now%strn2D%dyy 
                mat%now%strn%dxy(:,:,k) = mat%now%strn2D%dxy 
                mat%now%strn%de(:,:,k)  = mat%now%strn2D%de 
            end do 

            mat%now%f_shear_bar = 1.0 

        end if 

        ! 1. Update enhancement factor ======================

        select case(trim(mat%par%enh_method))

            case("simple","simple-tracer")
                ! Grounded ice: enh = enh_shear 
                ! Floating ice: enh = enh_shlf

                ! First define spatially varying enhancement factor (2D only),
                ! for lowest layer of 3D enh field
                ! Specify enh_stream = enh_shear 
                mat%now%enh(:,:,1) = define_enhancement_factor_2D(tpo%now%f_grnd,mat%now%f_shear_bar,dyn%now%uxy(:,:,nz_aa), &
                                                               mat%par%enh_shear,mat%par%enh_shear,mat%par%enh_shlf)
            
                ! Fill in the remaining 3D enh field layers too 
                do k = 2, nz_aa
                    mat%now%enh(:,:,k) = mat%now%enh(:,:,1)
                end do 

            case("shear2D","shear2D-tracer")
                ! Calculate 2D enhancement factor based on depth-averaged
                ! shear fraction (f_shear_bar)
                ! enh = enh_shear*f_shear_bar + enh_stream*(1-f_shear_bar)
                ! Note, floating ice always has enh = enh_shelf 

                ! First define spatially varying enhancement factor (2D only),
                ! for lowest layer of 3D enh field
                mat%now%enh(:,:,1) = define_enhancement_factor_2D(tpo%now%f_grnd,mat%now%f_shear_bar,dyn%now%uxy(:,:,nz_aa), &
                                                               mat%par%enh_shear,mat%par%enh_stream,mat%par%enh_shlf)
            
                ! Fill in the remaining 3D enh field layers too 
                do k = 2, nz_aa
                    mat%now%enh(:,:,k) = mat%now%enh(:,:,1)
                end do 

            case("shear3D","shear3D-tracer") 
                ! Calculate 3D enhancement factor based on 3D
                ! shear fraction field (f_shear)
                ! enh = enh_shear*f_shear_bar + enh_stream*(1-f_shear_bar)
                ! Note, floating ice always has enh = enh_shelf 
                
                ! Define spatially varying enhancement factor
                mat%now%enh = define_enhancement_factor_3D(mat%now%strn%f_shear,tpo%now%f_grnd,dyn%now%uxy(:,:,nz_aa), &
                                                           mat%par%enh_shear,mat%par%enh_stream,mat%par%enh_shlf)

            case DEFAULT 

                write(*,*) "calc_ymat:: Error: enhancement method not recognized: "//trim(mat%par%enh_method)
                stop 

        end select 

        ! If enh_method is one of the "*-tracer" methods, then
        ! additionally scale enh field by enh_bnd tracer field. 
        select case(trim(mat%par%enh_method))

            case("simple-tracer","shear2D-tracer","shear3D-tracer")
                ! Calculate 3D enhancement factor multiplier enh_bnd as the evolution 
                ! of an imposed enhancement factor at the surface propogating
                ! as a tracer inside of the ice sheet. Assume that propogation 
                ! is only valid for slow-flowing (ie, shearing ice), and impose 
                ! value of enh_bnd=1.0 for the fast-flowing and floating ice, respectively. 

                if (dt .gt. 0.0) then 
                    ! Update anisotropic enhancement factor tracer field if advancing timestep 
                    ! (if not, do nothing) 

                    ! Set surface boundary condition to boundary enh field
                    X_srf = bnd%enh_srf  

                    ! Define limits to where to calculate tracers, 
                    ! surface value will be imposed in fast regions
                    mask_tracers = .TRUE. 
                    where (dyn%now%uxy_bar .gt. mat%par%enh_umax) mask_tracers = .FALSE. 

                    call calc_tracer_3D(mat%now%enh_bnd,X_srf,dyn%now%ux,dyn%now%uy,dyn%now%uz,tpo%now%H_ice, &
                                        tpo%now%bmb,mat%par%zeta_aa,mat%par%zeta_ac,mat%par%tracer_method, &
                                        mat%par%tracer_impl_kappa,dt,thrm%par%dx,time,mask=mask_tracers)

                end if 

                ! Ensure enh_bnd is always non-zero and positive,
                ! but also not extremely high (eg 0.1 <= enh <= 10)
                where (mat%now%enh_bnd .lt. enh_min) mat%now%enh_bnd = enh_min
                where (mat%now%enh_bnd .gt. enh_max) mat%now%enh_bnd = enh_max
                
                ! Additionally update field to impose a value of one in streaming/floating regimes 
                call modify_enhancement_factor_bnd(mat%now%enh_bnd,tpo%now%f_grnd,dyn%now%uxy_bar,enh_stream=1.0_prec, &
                                enh_shlf=1.0_prec,umin=mat%par%enh_umin,umax=mat%par%enh_umax)

        
                ! Finally scale enh by enh_bnd 
                mat%now%enh = mat%now%enh * mat%now%enh_bnd 
                
        end select 

        ! Finally get the vertical average
        mat%now%enh_bar = calc_vertical_integrated_2D(mat%now%enh,mat%par%zeta_aa)


        ! 2. Update rate factor ==========================

        select case(mat%par%rf_method)

            case(-1)
                ! rater factor has been defined externally - do nothing

            case(0)
                ! Use constant parameter value for rate factor 

                mat%now%ATT     = mat%par%rf_const 
                mat%now%ATT_bar = mat%par%rf_const 

            case (1) 
                ! Calculate rate factor from ice temp., enhancement factor and water content 

                if (mat%par%rf_use_eismint2) then 
                    ! Use EISMINT2 (Payne et al, 2000) constants
                    mat%now%ATT = calc_rate_factor_eismint(thrm%now%T_ice,thrm%now%T_pmp,mat%now%enh)
                else 
                    ! Use Greve and Blatter (2009) constants 
                    mat%now%ATT = calc_rate_factor(thrm%now%T_ice,thrm%now%T_pmp,mat%now%enh)
                end if 

                ! Scale rate factor by water content if desired 
                if (mat%par%rf_with_water) then 
                    call scale_rate_factor_water(mat%now%ATT,thrm%now%omega)
                end if 

                ! Ensure rate factor is relatively smooth
!                 do k = 1, nz_aa
!                     call regularize2D(mat%now%ATT(:,:,k),tpo%now%H_ice,tpo%par%dx)
!                 end do 

                ! Get vertically averaged value 
                mat%now%ATT_bar = calc_vertical_integrated_2D(mat%now%ATT,mat%par%zeta_aa)
            
            case(2)
                ! Calculate rate factor from viscosity and effective strain rate 
                ! (only works when dyn%par%visc_method=0 with visc_const prescribed,
                ! and n_glen=1)

                if (dyn%par%visc_method .ne. 0 .or. mat%par%n_glen .ne. 1.0_wp) then 
                    write(*,*) "calc_ymat:: Error: rf_method=2 only works when viscosity &
                    &is prescribed (ydyn.visc_method=0) and ymat.n_glen=1.0."
                    write(*,*) "ydyn.visc_method = ", dyn%par%visc_method 
                    write(*,*) "ymat.n_glen      = ", mat%par%n_glen 
                    stop 
                end if 

                ! ATT = (2.0*visc_eff)^(-1) 
                
                mat%now%ATT     = 1.0_wp / (2.0_wp*dyn%par%visc_const)
                mat%now%ATT_bar = 1.0_wp / (2.0_wp*dyn%par%visc_const)
            
            case DEFAULT 
                ! Not recognized 

                write(*,*) "calc_ymat:: Error: rf_method not recognized."
                write(*,*) "rf_method = ", mat%par%rf_method
                stop 

        end select 

        ! 2. Calculate the updated visocity =====
        ! Note: this viscosity which is used for material and thermodynamic properties may be
        ! different than the visocity calculated locally in the dynamics module
        
        mat%now%visc = calc_viscosity_glen(mat%now%strn%de,mat%now%ATT,mat%par%n_glen,mat%par%visc_min)
        
        ! Ensure viscosity is relatively smooth
!         do k = 1, nz_aa
!             call regularize2D(mat%now%visc(:,:,k),tpo%now%H_ice,tpo%par%dx)
!         end do 

        ! Calculate visc_int (vertically integrated visc) as diagnostic quantity
        mat%now%visc_int = calc_vertical_integrated_2D(mat%now%visc,mat%par%zeta_aa)
        where(tpo%now%H_ice .gt. 0.0) mat%now%visc_int = mat%now%visc_int*tpo%now%H_ice 
        
        return
        
    end subroutine calc_ymat
    
    subroutine ymat_par_load(par,filename,zeta_aa,zeta_ac,nx,ny,dx,init)

        implicit none 

        type(ymat_param_class), intent(OUT) :: par
        character(len=*),       intent(IN)  :: filename
        real(prec),             intent(IN)  :: zeta_aa(:)
        real(prec),             intent(IN)  :: zeta_ac(:)
        integer,                intent(IN)  :: nx, ny 
        real(prec),             intent(IN)  :: dx  
        logical, optional,      intent(IN)  :: init 

        ! Local variables 
        logical :: init_pars 
        real(prec) :: age_iso(100) 

        age_iso = 0.0 

        init_pars = .FALSE.
        if (present(init)) init_pars = .TRUE. 
 
        ! Store local parameter values in output object
        call nml_read(filename,"ymat","flow_law",               par%flow_law,               init=init_pars)
        call nml_read(filename,"ymat","rf_method",              par%rf_method,              init=init_pars)
        call nml_read(filename,"ymat","rf_const",               par%rf_const,               init=init_pars)
        call nml_read(filename,"ymat","rf_use_eismint2",        par%rf_use_eismint2,        init=init_pars)
        call nml_read(filename,"ymat","rf_with_water",          par%rf_with_water,          init=init_pars)
        call nml_read(filename,"ymat","n_glen",                 par%n_glen,                 init=init_pars)
        call nml_read(filename,"ymat","visc_min",               par%visc_min,               init=init_pars)
        call nml_read(filename,"ymat","de_max",                 par%de_max,                 init=init_pars)
        call nml_read(filename,"ymat","enh_method",             par%enh_method,             init=init_pars)
        call nml_read(filename,"ymat","enh_shear",              par%enh_shear,              init=init_pars)
        call nml_read(filename,"ymat","enh_stream",             par%enh_stream,             init=init_pars)
        call nml_read(filename,"ymat","enh_shlf",               par%enh_shlf,               init=init_pars)
        call nml_read(filename,"ymat","enh_umin",               par%enh_umin,               init=init_pars)
        call nml_read(filename,"ymat","enh_umax",               par%enh_umax,               init=init_pars)
        call nml_read(filename,"ymat","calc_age",               par%calc_age,               init=init_pars)
        call nml_read(filename,"ymat","age_iso",                age_iso,                    init=init_pars)
        call nml_read(filename,"ymat","tracer_method",          par%tracer_method,          init=init_pars)
        call nml_read(filename,"ymat","tracer_impl_kappa",      par%tracer_impl_kappa,      init=init_pars)
        
        ! Set internal parameters
        par%nx    = nx 
        par%ny    = ny 
        par%dx    = dx
        par%dy    = dx  
        par%nz_aa = size(zeta_aa,1)  
        par%nz_ac = size(zeta_ac,1)
        
        if (allocated(par%zeta_aa)) deallocate(par%zeta_aa)
        allocate(par%zeta_aa(par%nz_aa))
        par%zeta_aa = zeta_aa 
        
        if (allocated(par%zeta_ac)) deallocate(par%zeta_ac)
        allocate(par%zeta_ac(par%nz_ac))
        par%zeta_ac = zeta_ac 
        
        if (count(age_iso .eq. 0.0) .eq. size(age_iso)) then 
            ! No isochrones to be calculated, fill with one layer for present day (age=0)
            par%n_iso = 1 
        else 
            ! Assume all isochrone ages are not equal to present day
            par%n_iso = count(age_iso .ne. 0.0) 
        end if 

        if (allocated(par%age_iso)) deallocate(par%age_iso)
        allocate(par%age_iso(par%n_iso))
        par%age_iso = age_iso(1:par%n_iso)

        ! Define current time as unrealistic value
        par%time = 1000000000   ! [a] 1 billion years in the future

        return 

    end subroutine ymat_par_load

    subroutine ymat_alloc(now,nx,ny,nz_aa,nz_ac,n_iso)

        implicit none 

        type(ymat_state_class), intent(INOUT) :: now 
        integer :: nx, ny, nz_aa, nz_ac, n_iso   

        ! First make sure fields are deallocated
        call ymat_dealloc(now)

        ! Allocate fields to desired dimensions 
        allocate(now%strn2D%dxx(nx,ny))
        allocate(now%strn2D%dyy(nx,ny))
        allocate(now%strn2D%dxy(nx,ny))
        allocate(now%strn2D%de(nx,ny))
        
        allocate(now%strn%dxx(nx,ny,nz_aa))
        allocate(now%strn%dyy(nx,ny,nz_aa))
        allocate(now%strn%dzz(nx,ny,nz_aa))
        allocate(now%strn%dxy(nx,ny,nz_aa))
        allocate(now%strn%dxz(nx,ny,nz_aa))
        allocate(now%strn%dyz(nx,ny,nz_aa))
        allocate(now%strn%de(nx,ny,nz_aa))
        allocate(now%strn%f_shear(nx,ny,nz_aa))
        
        allocate(now%enh(nx,ny,nz_aa))
        allocate(now%enh_bnd(nx,ny,nz_aa))
        allocate(now%enh_bar(nx,ny))
        allocate(now%ATT(nx,ny,nz_aa))
        allocate(now%ATT_bar(nx,ny))
        
        allocate(now%visc(nx,ny,nz_aa))
        allocate(now%visc_int(nx,ny))

        allocate(now%f_shear_bar(nx,ny)) 

        allocate(now%dep_time(nx,ny,nz_aa)) 
        allocate(now%depth_iso(nx,ny,n_iso)) 

        now%strn2D%dxx   = 0.0 
        now%strn2D%dyy   = 0.0 
        now%strn2D%dxy   = 0.0
        now%strn2D%de    = 0.0 
        
        now%strn%dxx     = 0.0 
        now%strn%dyy     = 0.0 
        now%strn%dzz     = 0.0
        now%strn%dxy     = 0.0 
        now%strn%dxz     = 0.0
        now%strn%dyz     = 0.0
        now%strn%de      = 0.0
        now%strn%f_shear = 0.0 
     
        now%enh          = 1.0 
        now%enh_bnd      = 1.0 
        now%enh_bar      = 0.0 
        now%ATT          = 0.0 
        now%ATT_bar      = 0.0    
        now%visc         = 0.0 
        now%visc_int     = 0.0 

        now%f_shear_bar  = 0.0 

        now%dep_time     = 0.0 
        now%depth_iso    = 0.0 
        
        return 

    end subroutine ymat_alloc 

    subroutine ymat_dealloc(now)

        implicit none 

        type(ymat_state_class), intent(INOUT) :: now 

        if (allocated(now%strn2D%dxx))      deallocate(now%strn2D%dxx)
        if (allocated(now%strn2D%dyy))      deallocate(now%strn2D%dyy)
        if (allocated(now%strn2D%dxy))      deallocate(now%strn2D%dxy)
        if (allocated(now%strn2D%de))       deallocate(now%strn2D%de)
        
        if (allocated(now%strn%dxx))        deallocate(now%strn%dxx)
        if (allocated(now%strn%dyy))        deallocate(now%strn%dyy)
        if (allocated(now%strn%dxy))        deallocate(now%strn%dxy)
        if (allocated(now%strn%dxz))        deallocate(now%strn%dxz)
        if (allocated(now%strn%dyz))        deallocate(now%strn%dyz)
        if (allocated(now%strn%de))         deallocate(now%strn%de)
        if (allocated(now%strn%f_shear))    deallocate(now%strn%f_shear)
        
        if (allocated(now%enh))             deallocate(now%enh)
        if (allocated(now%enh_bnd))         deallocate(now%enh_bnd)
        if (allocated(now%enh_bar))         deallocate(now%enh_bar)
        if (allocated(now%ATT))             deallocate(now%ATT)
        if (allocated(now%ATT_bar))         deallocate(now%ATT_bar)

        if (allocated(now%visc))            deallocate(now%visc)
        if (allocated(now%visc_int))        deallocate(now%visc_int)
        
        if (allocated(now%f_shear_bar))     deallocate(now%f_shear_bar)

        if (allocated(now%dep_time))        deallocate(now%dep_time)
        if (allocated(now%depth_iso))       deallocate(now%depth_iso)
        
        return 

    end subroutine ymat_dealloc 

end module yelmo_material
