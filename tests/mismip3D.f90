module mismip3D

    use ncio 
    use yelmo_defs, only : sp, dp, prec, sec_year, T0  
    
    implicit none 

    
    private 
    public :: mismip3D_topo_init
    public :: mismip3D_boundaries
    public :: perturb_friction
    public :: find_x_gl
    public :: find_x_gl_2D 

contains 

    subroutine mismip3D_topo_init(z_bed,H_ice,z_srf,x,y,experiment)

        implicit none 

        real(prec), intent(OUT) :: z_bed(:,:) 
        real(prec), intent(OUT) :: H_ice(:,:) 
        real(prec), intent(OUT) :: z_srf(:,:) 
        real(prec), intent(IN)  :: x(:), y(:)  
        character(len=*), intent(IN) :: experiment 

        ! Local variables 
        integer :: i, j, nx, ny 

        nx = size(z_bed,1)
        ny = size(z_bed,2)
        
        ! Define bedrock as a slope                                                                
        select case(trim(experiment))

            case("Stnd","P75S")
                ! MISMIP3D - Pattyn (2012)
                
                do j = 1, ny
                    z_bed(:,j) = -100.0 - x
                end do 

            case("RF")
                ! Schoof (2007)

                do j = 1, ny
                    z_bed(:,j) = 720.0 - 778.5*(x/750.0)
                end do 

            case DEFAULT 

                write(*,*) "mismip3D_topo_init:: Error: experiment not recognized."
                write(*,*) "experiment = ", trim(experiment)
                stop 

        end select 

        ! Set ice thickness to 10 m everywhere initially 
        H_ice = 10.0

!         ! Set ice thickness to 1000 m everywhere initially 
!         ! But linearly decreasing thickness with bed depth 
!         H_ice = 1000.0
!         where(z_bed .lt. 0.0) H_ice = max(0.0,1000.0-z_bed*0.9)

        ! Remove ice from deep bed to ensure ice is floating for most of 
        ! region of domain of interest, and initial grounding line advances
        ! to equilibrium, since here we are testing retreat 
        where(z_bed .lt. -500.0) H_ice = 0.0 

        ! Adjust for floating ice later, for now assume fully grounded
        z_srf = z_bed + H_ice 

        return 

    end subroutine mismip3D_topo_init

    subroutine mismip3D_boundaries(T_srf,smb,ghf,experiment)

        implicit none 

        real(prec), intent(OUT) :: T_srf(:,:) 
        real(prec), intent(OUT) :: smb(:,:) 
        real(prec), intent(OUT) :: ghf(:,:) 
        
        character(len=*), intent(IN) :: experiment 

        select case(trim(experiment))

            case("Stnd","Stnd-0.3","P75S")

                ! Surface temperature 
                T_srf = T0 

                ! Surface mass balance 
                smb = 0.5    ! [m/a]

                if (trim(experiment) .eq. "Stnd-0.3") smb = 0.3   ! [m/a] 

                ! Geothermal heat flux 
                ghf = 42.0   ! [mW/m2]

            case("RF")

                ! Surface temperature 
                T_srf = T0 

                ! Surface mass balance 
                smb = 0.3    ! [m/a]

                ! Geothermal heat flux 
                ghf = 42.0   ! [mW/m2]

            case DEFAULT 

                write(*,*) "Experiment not recognized: "//trim(experiment)
                stop 

        end select 

        return 

    end subroutine mismip3D_boundaries

    subroutine perturb_friction(cf_bed,xx,yy,x_gl)

        implicit none 

        real(prec), intent(INOUT) :: cf_bed(:,:) 
        real(prec), intent(IN)    :: xx(:,:) 
        real(prec), intent(IN)    :: yy(:,:) 
        real(prec), intent(IN)    :: x_gl

        ! Local variables 
        real(prec) :: xc 
        real(prec) :: yc 
        real(prec) :: y_gl, a_ref  
        
        ! For Stnd experiment: cf_bed_ref = 1e7     [Pa m^-1/3 s^1/3] / (1/sec_year)^1/3 = 3.165176e4
        ! For RF experiment:   cf_bed_ref = 7.624e6 [Pa m^-1/3 s^1/3] / (1/sec_year)^1/3 = 2.412596e4
        real(prec), parameter :: cf_bed_ref = 3.165176e4
        !real(prec), parameter :: cf_bed_ref = 2.412596e4

        ! Set bed friction coefficient
        xc    = 150.0 
        yc    =  10.0 
        y_gl  =   0.0  
        a_ref =  0.75 
        cf_bed = cf_bed_ref * &
            ( 1.0 - a_ref*exp(-(xx-x_gl)**2/(2.0*xc**2)-(yy-y_gl)**2/(2.0*yc**2)) )

        return 

    end subroutine perturb_friction

    function find_x_gl(xx,yy,H_grnd) result(x_gl)
        ! Find the position of the grounding line in the x-direction 
        implicit none 

        real(prec), intent(IN) :: xx(:,:) 
        real(prec), intent(IN) :: yy(:,:) 
        real(prec), intent(IN) :: H_grnd(:,:) 
        real(prec) :: x_gl 

        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: inow, jnow 
        real(prec) :: f_grnd 

        nx = size(xx,1)
        ny = size(xx,2) 

        ! Determine y index of y=0km
        jnow = minloc(yy(1,:),dim=1) 

        ! Determine x index of grounding line
        if (H_grnd(1,jnow) .le. 0.0) then 
            ! No grounded ice 

            x_gl = 0.0 

        else 
            ! Grounded ice exists, find grounding line
            do i = 1, nx-1 
                if (H_grnd(i,jnow) .gt. 0.0 .and. H_grnd(i+1,jnow) .le. 0.0) then 
                    inow = i 
                    exit 
                end if 
            end do 

            ! Determine x-coordinate of grounding line,
            ! where H_grnd==0.0 
            ! Interpolation: H = H1 + f*(H2-H1), so when H = 0
            ! f = -H1 / (H2-H1)  
            ! f_grnd of index inow should be a fraction of the cell 
            f_grnd = -H_grnd(inow,jnow) / (H_grnd(inow+1,jnow)-H_grnd(inow,jnow))
            x_gl   = xx(inow,jnow) + f_grnd*(xx(inow+1,jnow)-xx(inow,jnow))

        end if 

        return 
        
    end function find_x_gl

    subroutine find_x_gl_2D(x_gl,x_gl_std,xx,yy,f_grnd)
        ! Find the position of the grounding line in the x-direction 
        implicit none 

        real(prec), intent(OUT) :: x_gl
        real(prec), intent(OUT) :: x_gl_std 
        real(prec), intent(IN)  :: xx(:,:) 
        real(prec), intent(IN)  :: yy(:,:) 
        real(prec), intent(IN)  :: f_grnd(:,:) 
        
        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: n_gl  

        real(prec), allocatable :: x_gl_2D(:,:) 

        nx = size(xx,1)
        ny = size(xx,2) 

        allocate(x_gl_2D(nx,ny))
        x_gl_2D = 0.0 

        do j = 1, ny 
        do i = 1, nx 

            if (f_grnd(i,j) .gt. 0.0 .and. &
                    count(f_grnd(i-1:i+1,j-1:j+1) .eq. 0.0) .gt. 0) then 
                ! Grounding line point 

                x_gl_2D(i,j) = sqrt(xx(i,j)**2+yy(i,j)**2)

            end if 

        end do 
        end do 

        n_gl = count(x_gl_2D .gt. 0.0) 

        if (n_gl .gt. 0.0) then 

            x_gl     = sum(x_gl_2D,mask=x_gl_2D.gt.0.0) / real(n_gl,prec)
            x_gl_std = sqrt(sum((x_gl_2D-x_gl)**2,mask=x_gl_2D.gt.0.0)) / real(n_gl,prec)

        else 

            x_gl     = 0.0 
            x_gl_std = 0.0 

        end if 

        return 
        
    end subroutine find_x_gl_2D

end module mismip3D

