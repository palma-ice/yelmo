
module yelmo_material

    use nml 

    use yelmo_defs
    use yelmo_tools, only : calc_vertical_integrated_2D, calc_vertical_integrated_3D
    
    use deformation
    use iceage 

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

        nz_aa = mat%par%nz_aa

        ! Initialize time if necessary 
        if (mat%par%time .gt. time) then 
            mat%par%time = time
        end if 

        ! Get time step and advance current time 
        dt            = time - mat%par%time 
        mat%par%time  = time 
        
        ! 0. Update strain rate 
        ! Note: this calculation of strain rate is particular to the material module, and 
        ! may differ from the strain rate calculated locally in the dynamics module

        ! Calculate the strain rate from the full 3D tensor

        if (.TRUE.) then 
            ! 3D strain rate - standard 

            call calc_strain_rate_3D(mat%now%strn,dyn%now%ux,dyn%now%uy,dyn%now%uz,tpo%now%H_ice, &
                                     tpo%now%f_grnd,mat%par%zeta_aa,mat%par%zeta_ac,mat%par%dx)


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

        if (mat%par%use_2D_enh) then 
            ! Calculate 2D enhancement factor 

            ! Next, define spatially varying enhancement factor (2D only)
            mat%now%enh_bar = define_enhancement_factor_2D(tpo%now%f_grnd,mat%now%f_shear_bar,dyn%now%uxy(:,:,nz_aa), &
                                                           mat%par%enh_shear,mat%par%enh_stream,mat%par%enh_shlf)
            
            ! Fill in 3D enh field too 
            do k = 1, nz_aa
                mat%now%enh(:,:,k) = mat%now%enh_bar
            end do 

        else 
            ! Calculate 3D enhancement factor 

            ! Next, define spatially varying enhancement factor
            mat%now%enh = define_enhancement_factor(mat%now%strn%f_shear,tpo%now%f_grnd,dyn%now%uxy(:,:,nz_aa), &
                                                    mat%par%enh_shear,mat%par%enh_stream,mat%par%enh_shlf)

            ! And get the vertical average
            mat%now%enh_bar     = calc_vertical_integrated_2D(mat%now%enh,mat%par%zeta_aa)
 
        end if 

        ! 2. Update rate factor ==========================

        select case(mat%par%rf_method)

            case(-1)
                ! rater factor has been defined externally - do nothing

            case(0)
                ! Use constant parameter value for rate factor 

                mat%now%ATT     = mat%par%rf_const 
                mat%now%ATT_bar = mat%par%rf_const 

            case (1) 
                ! Calculate rate factor from ice temp. and enhancement factor 

                if (mat%par%rf_use_eismint2) then 
                    ! Use EISMINT2 (Payne et al, 2000) constants
                    mat%now%ATT = calc_rate_factor_eismint(thrm%now%T_ice,thrm%now%T_pmp,mat%now%enh)
                else 
                    ! Use Greve and Blatter (2009) constants 
                    mat%now%ATT = calc_rate_factor(thrm%now%T_ice,thrm%now%T_pmp,mat%now%enh)
                end if 

                ! Get vertically averaged value 
                mat%now%ATT_bar = calc_vertical_integrated_2D(mat%now%ATT,mat%par%zeta_aa)
            
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
        
        ! Calculate visc_bar (vertically averaged visc) as diagnostic quantity
        mat%now%visc_bar = calc_vertical_integrated_2D(mat%now%visc,mat%par%zeta_aa)
        where(tpo%now%H_ice .gt. 0.0) mat%now%visc_bar = mat%now%visc_bar*tpo%now%H_ice 
        
        if (mat%par%calc_age .and. dt .gt. 0.0) then 
            ! Perform calculations of age tracer: dep_time (deposition time)

            call calc_tracer_3D(mat%now%dep_time,dyn%now%ux,dyn%now%uy,dyn%now%uz,tpo%now%H_ice, &
                tpo%now%bmb,mat%par%zeta_aa,mat%par%zeta_ac,thrm%par%dzeta_a,thrm%par%dzeta_b, &
                mat%par%age_method,mat%par%age_impl_kappa,dt,thrm%par%dx,time)

        end if 

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

        init_pars = .FALSE.
        if (present(init)) init_pars = .TRUE. 
 
        ! Store local parameter values in output object
        call nml_read(filename,"ymat","flow_law",               par%flow_law,               init=init_pars)
        call nml_read(filename,"ymat","rf_method",              par%rf_method,              init=init_pars)
        call nml_read(filename,"ymat","rf_const",               par%rf_const,               init=init_pars)
        call nml_read(filename,"ymat","rf_use_eismint2",        par%rf_use_eismint2,        init=init_pars)
        call nml_read(filename,"ymat","n_glen",                 par%n_glen,                 init=init_pars)
        call nml_read(filename,"ymat","visc_min",               par%visc_min,               init=init_pars)
        call nml_read(filename,"ymat","use_2D_enh",             par%use_2D_enh,             init=init_pars)
        call nml_read(filename,"ymat","enh_shear",              par%enh_shear,              init=init_pars)
        call nml_read(filename,"ymat","enh_stream",             par%enh_stream,             init=init_pars)
        call nml_read(filename,"ymat","enh_shlf",               par%enh_shlf,               init=init_pars)
        
        call nml_read(filename,"ymat","age_method",             par%age_method,             init=init_pars)
        call nml_read(filename,"ymat","age_impl_kappa",         par%age_impl_kappa,         init=init_pars)
        
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
        
        ! Define current time as unrealistic value
        par%time = 1000000000   ! [a] 1 billion years in the future

        ! Determine whether age is actually being calculated or not 
        par%calc_age = .TRUE. 
        if (trim(par%age_method) .eq. "None") par%calc_age = .FALSE. 

        return 

    end subroutine ymat_par_load

    subroutine ymat_alloc(now,nx,ny,nz_aa,nz_ac)

        implicit none 

        type(ymat_state_class), intent(INOUT) :: now 
        integer :: nx, ny, nz_aa, nz_ac  

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
        allocate(now%enh_bar(nx,ny))
        allocate(now%ATT(nx,ny,nz_aa))
        allocate(now%ATT_bar(nx,ny))
        
        allocate(now%visc(nx,ny,nz_aa))
        allocate(now%visc_bar(nx,ny))

        allocate(now%f_shear_bar(nx,ny)) 

        allocate(now%dep_time(nx,ny,nz_aa)) 

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
     
        now%enh          = 0.0 
        now%enh_bar      = 0.0 
        now%ATT          = 0.0 
        now%ATT_bar      = 0.0    
        now%visc         = 0.0 
        now%visc_bar     = 0.0 

        now%f_shear_bar  = 0.0 

        now%dep_time     = 0.0 

        return 

    end subroutine ymat_alloc 

    subroutine ymat_dealloc(now)

        implicit none 

        type(ymat_state_class), intent(INOUT) :: now 

        if (allocated(now%strn2D%dxx))  deallocate(now%strn2D%dxx)
        if (allocated(now%strn2D%dyy))  deallocate(now%strn2D%dyy)
        if (allocated(now%strn2D%dxy))  deallocate(now%strn2D%dxy)
        if (allocated(now%strn2D%de))   deallocate(now%strn2D%de)
        
        if (allocated(now%strn%dxx))      deallocate(now%strn%dxx)
        if (allocated(now%strn%dyy))      deallocate(now%strn%dyy)
        if (allocated(now%strn%dxy))      deallocate(now%strn%dxy)
        if (allocated(now%strn%dxz))      deallocate(now%strn%dxz)
        if (allocated(now%strn%dyz))      deallocate(now%strn%dyz)
        if (allocated(now%strn%de))       deallocate(now%strn%de)
        if (allocated(now%strn%f_shear))  deallocate(now%strn%f_shear)
        
        if (allocated(now%enh))         deallocate(now%enh)
        if (allocated(now%enh_bar))     deallocate(now%enh_bar)
        if (allocated(now%ATT))         deallocate(now%ATT)
        if (allocated(now%ATT_bar))     deallocate(now%ATT_bar)

        if (allocated(now%visc))        deallocate(now%visc)
        if (allocated(now%visc_bar))    deallocate(now%visc_bar)
        
        if (allocated(now%f_shear_bar)) deallocate(now%f_shear_bar)

        if (allocated(now%dep_time))    deallocate(now%dep_time)
        
        return 

    end subroutine ymat_dealloc 

end module yelmo_material
