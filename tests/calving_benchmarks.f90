module calving_benchmarks

    use yelmo_defs, only : wp, pi

    implicit none

    private
    public :: calvmip_init
    public :: set_calving_mask

contains

    subroutine calvmip_init(z_bed,xx,yy,domain)

        implicit none

        real(wp), intent(OUT) :: z_bed(:,:)

        real(wp), intent(IN)  :: xx(:,:)
        real(wp), intent(IN)  :: yy(:,:)
        character(len=*), intent(IN) :: domain
        
        ! Local variables
        real(wp) :: R0, Bc, Bl, Ba, rc

        select case(trim(domain))

            case("circular")
                
                R0 = 800e3;
                Bc = 900; 
                Bl = -2000;   
                rc = 0;
                
                call define_bedrock_circular(z_bed,xx,yy,R0,Bc,Bl,rc)

            case("thule")

                R0 = 800e3;
                Bc = 900; 
                Bl = -2000; 
                Ba = 1100;  
                rc = 0;

                call define_bedrock_thule(z_bed,xx,yy,R0,Bc,Bl,Ba,rc)

            case DEFAULT

                write(*,*) "calvmip_init:: Error: domain not recognized."
                write(*,*) "domain = ", trim(domain)
                stop 

        end select

        return

    end subroutine calvmip_init

    elemental subroutine define_bedrock_circular(z_bed,xx,yy,R0,Bc,Bl,rc)
        ! Details here:
        ! https://github.com/JRowanJordan/CalvingMIP/wiki/circular-domain

        implicit none

        real(wp), intent(OUT) :: z_bed
        real(wp), intent(IN)  :: xx
        real(wp), intent(IN)  :: yy
        real(wp), intent(IN)  :: R0
        real(wp), intent(IN)  :: Bc
        real(wp), intent(IN)  :: Bl
        real(wp), intent(IN)  :: rc
        
        ! Local variables
        real(wp) :: r
        real(wp) :: theta

        r = sqrt(xx*xx+yy*yy);
        theta = atan2(yy,xx);
        z_bed = Bc-(Bc-Bl)*(r-rc)**2 / (R0-rc)**2;

        return

    end subroutine define_bedrock_circular

    elemental subroutine define_bedrock_thule(z_bed,xx,yy,R0,Bc,Bl,Ba,rc)
        ! Details here: 
        ! https://github.com/JRowanJordan/CalvingMIP/wiki/Thule-domain

        implicit none

        real(wp), intent(OUT) :: z_bed
        real(wp), intent(IN)  :: xx
        real(wp), intent(IN)  :: yy
        real(wp), intent(IN)  :: R0
        real(wp), intent(IN)  :: Bc
        real(wp), intent(IN)  :: Bl
        real(wp), intent(IN)  :: Ba
        real(wp), intent(IN)  :: rc
        
        ! Local variables
        real(wp) :: r
        real(wp) :: theta
        real(wp) :: l, a 

        r = sqrt(xx*xx+yy*yy);
        theta = atan2(yy,xx);
        l = R0-cos(2*theta)*R0/2;
        a = Bc-(Bc-Bl)*(r-rc)**2 / (R0-rc)**2;
        z_bed = Ba*cos(3*pi*r / l)+a ;

        return

    end subroutine define_bedrock_thule

    elemental subroutine set_calving_mask(mask,xx,yy,r_lim)

        implicit none

        logical, intent(OUT) :: mask
        real(wp), intent(IN)  :: xx
        real(wp), intent(IN)  :: yy
        real(wp), intent(IN)  :: r_lim 

        real(wp) :: r
        
        r = sqrt(xx*xx+yy*yy);

        if (r .gt. r_lim) then 
            mask = .TRUE.
        else
            mask = .FALSE.
        end if

        return

    end subroutine set_calving_mask

    subroutine calving_boundaries(T_srf,smb,ghf,xx,yy,H_ice,experiment,time,period,rad_el,dT_test,dsmb_test)

        implicit none

        real(wp), intent(OUT) :: T_srf(:,:)
        real(wp), intent(OUT) :: smb(:,:)
        real(wp), intent(OUT) :: ghf(:,:)
        real(wp), intent(IN)  :: xx(:,:)      ! [m]
        real(wp), intent(IN)  :: yy(:,:)      ! [m]
        real(wp), intent(IN)  :: H_ice(:,:)
        character(len=*), intent(IN) :: experiment
        real(wp), intent(IN) :: time
        real(wp), intent(IN) :: period
        real(wp), intent(IN), optional :: rad_el
        real(wp), intent(IN), optional :: dT_test
        real(wp), intent(IN), optional :: dsmb_test

        real(wp), parameter :: x_summit = 0.0
        real(wp), parameter :: y_summit = 0.0

    end subroutine calving_boundaries

end module calving_benchmarks
