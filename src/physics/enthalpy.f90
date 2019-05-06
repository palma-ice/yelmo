module enthalpy
    ! This module contains functions for computing the conversion of [T,w] => [enth] 
    ! and the inverse [enth,w] => [T] 
    ! Functions ported from sicopolis5-dev 

    use yelmo_defs, only : sp, dp, prec, rho_ice  
!     use yelmo_tools, only : integrate_trapezoid1D_1D  ! Not used yet, to do...

    implicit none 

    ! Private variables that remain populated for general use (tables)

    ! Temperature integral of the specific heat of ice.
    ! Index is temperature in deg C.
    real(prec), dimension(-256:255)       :: c_int_table

    ! Inverse of the temperature integral of the specific heat
    ! of ice. Index is enthalpy in J/kg (zero for 0 deg C).
    real(prec), dimension(-524288:524287) :: c_int_inv_table

    ! Lower index limit of properly defined values in c_int_table
    ! (n_temp_min >= -256).
    integer :: n_temp_min

    ! Upper index limit of properly defined values in c_int_table
    ! (n_temp_max <= 255).
    integer :: n_temp_max

    ! Lower index limit of properly defined values in c_int_inv_table
    ! (n_enth_min >= -524288).
    integer :: n_enth_min

    ! Upper index limit of properly defined values in c_int_inv_table
    ! (n_enth_max <= 524287).
    integer :: n_enth_max

    ! Latent heat of ice and inverse
    real(prec), parameter :: L     = 3.35e+05_prec    ! in J/kg
    real(prec), parameter :: L_inv = 1.0_prec/L
    
    private 
    public :: enthalpy_init_tables
    public :: enth_fct_temp_omega
    public :: temp_fct_enth
    public :: omega_fct_enth

contains 

    function enth_fct_temp_omega(temp_val, omega_val)
        !-------------------------------------------------------------------------------
        !> Enthalpy as a function of temperature and water content.
        !<------------------------------------------------------------------------------

        implicit none

        real(prec)              :: enth_fct_temp_omega

        real(prec), intent(in)  :: temp_val, omega_val

        enth_fct_temp_omega = c_int_val(temp_val) + L*omega_val

        return 

    end function enth_fct_temp_omega
    
    function temp_fct_enth(enth_val, temp_m_val)
        !-------------------------------------------------------------------------------
        !> Temperature as a function of enthalpy.
        !<------------------------------------------------------------------------------

        implicit none

        real(prec)             :: temp_fct_enth

        real(prec), intent(in) :: enth_val
        real(prec), intent(in) :: temp_m_val

        real(prec) :: enth_i

        enth_i = c_int_val(temp_m_val)   ! Enthalpy of pure ice at the melting point

        if (enth_val < enth_i) then   ! cold ice
           temp_fct_enth = c_int_inv_val(enth_val)
        else   ! temperate ice
           temp_fct_enth = temp_m_val
        end if

        return 

    end function temp_fct_enth
    
    function omega_fct_enth(enth_val, temp_m_val)
        !-------------------------------------------------------------------------------
        !> Water content as a function of enthalpy.
        !<------------------------------------------------------------------------------

        implicit none

        real(prec)             :: omega_fct_enth

        real(prec), intent(in) :: enth_val
        real(prec), intent(in) :: temp_m_val

        real(prec) :: enth_i

        enth_i = c_int_val(temp_m_val)   ! Enthalpy of pure ice at the melting point

        omega_fct_enth = max((enth_val-enth_i)*L_inv, 0.0_prec)

        return 

    end function omega_fct_enth


    subroutine enthalpy_init_tables()

        implicit none 

        real(prec) :: c_table(-190:10)
        integer :: i 

        ! Define c_table values 
        do i = -190, 10
            c_table(i) = calc_specific_heat_capacity(real(i,prec)+273.15)
        end do 

        call calc_c_int_table(c_table, -190, 10)

        call calc_c_int_inv_table()

!         write(*,*) "c_int_table: ", minval(c_int_table), maxval(c_int_table)
!         write(*,*) "c_int_inv_table: ", minval(c_int_inv_table), maxval(c_int_inv_table)
!         stop 
        
        return 

    end subroutine enthalpy_init_tables 

    elemental function calc_specific_heat_capacity(T_ice) result(cp)

        implicit none 

        real(prec), intent(IN) :: T_ice  
        real(prec) :: cp 

        ! Specific heat capacity (Greve and Blatter, 2009, Eq. 4.39; Ritz, 1987)
        cp = (146.3 +7.253*T_ice)    ! [J kg-1 K-1]

        return 

    end function calc_specific_heat_capacity

    ! ========================================================================
    !
    ! Conversion from temperature (temp) and water content (omega) to enthalpy
    ! (enth) and vice versa.
    !
    ! ========================================================================

    function c_int_val(temp_val)
        !-------------------------------------------------------------------------------
        !> Temperature integral of the specific heat of ice
        !! (enthalpy as function of temperature).
        !<------------------------------------------------------------------------------
        implicit none

        real(prec)             :: c_int_val

        real(prec), intent(in) :: temp_val

        integer :: n_temp_1, n_temp_2

        ! character(len=256) :: errormsgg

        n_temp_1 = floor(temp_val)

        n_temp_1 = max(min(n_temp_1, n_temp_max-1), n_temp_min)
        n_temp_2 = n_temp_1 + 1

        ! if ((n_temp_1 < n_temp_min-1).or.(n_temp_2 > n_temp_max+1)) then
        !    errormsgg = ' >>> c_int_val: Temperature argument out of allowed range!'
        !    call error(errormsgg)
        ! end if
        !    *** Commented out after some testing in order to save computing time. ***

        c_int_val = c_int_table(n_temp_1) &
                    + (c_int_table(n_temp_2)-c_int_table(n_temp_1)) &
                      * (temp_val-real(n_temp_1,prec))   ! Linear interpolation

        return 

    end function c_int_val

        
    function c_int_inv_val(enth_val)
        !-------------------------------------------------------------------------------
        !> Inverse function of c_int_val (temperature as function of enthalpy).
        !<------------------------------------------------------------------------------

        implicit none

        real(prec)             :: c_int_inv_val

        real(prec), intent(in) :: enth_val

        integer :: n_enth_1, n_enth_2

        ! character(len=256) :: errormsgg

        n_enth_1 = floor(enth_val)
        n_enth_1 = max(min(n_enth_1, n_enth_max-1), n_enth_min)
        n_enth_2 = n_enth_1 + 1

        ! if ((n_enth_1 < n_enth_min-1).or.(n_enth_2 > n_enth_max+1)) then
        !    errormsgg = ' >>> c_int_inv_val: Enthalpy argument out of allowed range!'
        !    call error(errormsgg)
        ! end if
        !    *** Commented out after some testing in order to save computing time. ***

        c_int_inv_val = c_int_inv_table(n_enth_1) &
                        + (c_int_inv_table(n_enth_2)-c_int_inv_table(n_enth_1)) &
                          * (enth_val-real(n_enth_1,prec))   ! Linear interpolation

        return 

    end function c_int_inv_val
    
    subroutine calc_c_int_table(c_table, n_tmp_min, n_tmp_max)
        !-------------------------------------------------------------------------------
        !> Computation of the temperature integral of the specific heat of ice as a
        !! table (c_int_table).
        !<------------------------------------------------------------------------------
        implicit none

        integer,                             intent(in) :: n_tmp_min, n_tmp_max
        real(prec), dimension(n_tmp_min:n_tmp_max), intent(in) :: c_table

        integer            :: n
        real(prec)         :: c_int_zero
        character(len=256) :: errormsgg

        !-------- Initialisation --------

        c_int_table = 0.0_prec

        n_temp_min = n_tmp_min
        n_temp_max = n_tmp_max

        if ((n_temp_min <= -256).or.(n_temp_max >= 255)) then
           write(*,*) ' >>> calc_c_int_table: ' &
                          //'Temperature indices out of allowed range!'
           stop 
        end if

        !-------- Numerical integration with the trapezoidal rule (spacing
        !         of data in c_table and c_int_table assumed to be 1 deg C) --------

        do n=n_temp_min+1, n_temp_max
           c_int_table(n) = c_int_table(n-1) + 0.5_prec*(c_table(n-1)+c_table(n))
                            ! that's the real stuff
        end do

        do n=n_temp_max+1, 255
           c_int_table(n) = c_int_table(n_temp_max)   ! dummy values
        end do

        !-------- Shift of the zero level to 0 deg C --------

        c_int_zero = c_int_table(0)

        do n=-256, 255
           c_int_table(n) = c_int_table(n) - c_int_zero
        end do

        return 

    end subroutine calc_c_int_table

    
    subroutine calc_c_int_inv_table()
        !-------------------------------------------------------------------------------
        !> Computation of the inverse of the temperature integral of the specific heat
        !! of ice as a table (c_int_inv_table).
        !<------------------------------------------------------------------------------

        implicit none

        integer       :: n
        integer       :: n_temp_1, n_temp_2
        real(prec)           :: enth_min, enth_max
        real(prec)           :: enth_val, enth_1, enth_2
        character(len=256) :: errormsgg

        !-------- Initialisation --------

        c_int_inv_table = 0.0_prec

        enth_min = c_int_val(real(n_temp_min,prec))
        enth_max = c_int_val(real(n_temp_max,prec))

        n_enth_min = ceiling(enth_min)
        n_enth_max = floor(enth_max)

        if ((n_enth_min <= -524288).or.(n_enth_max >= 524287)) then
           write(*,*) ' >>> calc_c_int_inv_table: ' &
                          //'Enthalpy indices out of allowed range!'
           stop 
        end if

        !-------- Linear interpolation between adjacent enthalpy values --------

        n_temp_1 = n_temp_min
        n_temp_2 = n_temp_min+1

        do n=n_enth_min, n_enth_max

           enth_val = real(n,prec)

           do

              if ((n_temp_1 > n_temp_max).or.(n_temp_2 > n_temp_max)) then
                 write(*,*) ' >>> calc_c_int_inv_table: ' &
                                //'Temperature indices out of allowed range!'
                 stop 
              end if

              enth_1 = c_int_val(real(n_temp_1,prec))
              enth_2 = c_int_val(real(n_temp_2,prec))

              if ( (enth_1 <= enth_val).and.(enth_2 >= enth_val) ) exit

              n_temp_1 = n_temp_1+1
              n_temp_2 = n_temp_2+1

           end do

           c_int_inv_table(n) = real(n_temp_1,prec) &
                                + (real(n_temp_2,prec)-real(n_temp_1,prec)) &
                                  * (enth_val-enth_1)/(enth_2-enth_1)
                                                      ! Linear interpolation
        end do

        do n=-524288, n_enth_min-1
           c_int_inv_table(n) = c_int_inv_table(n_enth_min)   ! dummy values
        end do

        do n=n_enth_max+1, 524287
           c_int_inv_table(n) = c_int_inv_table(n_enth_max)   ! dummy values
        end do

        return 

    end subroutine calc_c_int_inv_table
    
end module enthalpy 

