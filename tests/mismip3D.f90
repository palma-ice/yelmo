module mismip3D

    use ncio 
    use yelmo_defs, only : sp, dp, prec, sec_year, T0  
    
    implicit none 

    
    private 
    public :: mismip3D_topo_init
    public :: mismip3D_boundaries
    public :: find_x_gl

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

    subroutine mismip3D_boundaries(T_srf,smb,ghf,C_bed,xx,yy,x_gl,experiment)

        implicit none 

        real(prec), intent(OUT) :: T_srf(:,:) 
        real(prec), intent(OUT) :: smb(:,:) 
        real(prec), intent(OUT) :: ghf(:,:) 
        real(prec), intent(OUT) :: C_bed(:,:) 
        real(prec), intent(IN)  :: xx(:,:) 
        real(prec), intent(IN)  :: yy(:,:) 
        real(prec), intent(IN)  :: x_gl
        character(len=*), intent(IN) :: experiment 

        ! Local variables 
        integer    :: i, j, nx, ny 
        real(prec) :: y_gl, xc, yc, a_ref

        ! 1e7 [Pa^1/3 m^-1/3 s^1/3] *[a/sec_year]^1/3 = 31651.76 [Pa^1/3 m^-1/3 a^1/3]**3 = 3.170981e+13 [Pa^1 m^-1 a^1]
!         real(prec), parameter :: C_bed_ref = 31651.76**3    ! [Pa (m/a)^-1] == 1e7 [Pa^1/3 (m/s)^-1/3]
        real(prec), parameter :: C_bed_ref = 3.170981e13    ! [Pa (m/a)^-1] == 1e7 [Pa^1/3 (m/s)^-1/3]
        
        nx = size(T_srf,1)
        ny = size(T_srf,2)

        select case(trim(experiment))

            case("Stnd")

                ! Surface temperature 
                T_srf = T0 

                ! Surface mass balance 
                smb = 0.5    ! [m/a]

                ! Geothermal heat flux 
                ghf = 42.0   ! [mW/m2]

                ! Set bed friction coefficient 
                C_bed = C_bed_ref 

            case("RF")

                ! Surface temperature 
                T_srf = T0 

                ! Surface mass balance 
                smb = 0.3    ! [m/a]

                ! Geothermal heat flux 
                ghf = 42.0   ! [mW/m2]

                ! Set bed friction coefficient 
                C_bed = C_bed_ref 

            case("P75S")

                ! Surface temperature 
                T_srf = T0 

                ! Surface mass balance 
                smb = 0.5    ! [m/a]

                ! Geothermal heat flux 
                ghf = 42.0   ! [mW/m2]

                ! Set bed friction coefficient
                xc    = 150.0 
                yc    =  10.0 
                y_gl  =   0.0  
                a_ref = 0.75 
                C_bed = C_bed_ref * &
                    ( 1.0 - a_ref*exp(-(xx-x_gl)**2/(2.0*xc**2)-(yy-y_gl)**2/(2.0*yc**2)) )

            case DEFAULT 

                write(*,*) "Experiment not recognized: "//trim(experiment)
                stop 

        end select 

        return 

    end subroutine mismip3D_boundaries  

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

end module mismip3D 

