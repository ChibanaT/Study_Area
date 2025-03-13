program capacitor

    ! Purpose:
    !   To calculate the behavior of a capacitor as follows:
    !   1. If capacitance and voltage are known, calculate
    !       charge, number of electrons and energy stored.
    !   2. If charge and voltage are known, calculate capa-
    !       citance, number of electrons and energy stored.

    ! Record of Revisions:
    !       Date                Programmer              Description of Change
    !   ===========         =================       =============================
    !   2025/03/11               Chibana                Original code

    implicit none

    ! Data Dictionary: declare constants
    real, parameter :: electrons_per_coulomb = 6.241461e18

    ! Data Dictionary: declare variable types, definitions, & units
    real :: c               ! Capacitance of the capacitor (farads)
    real :: charge          ! Charge on the capacitor (coulombs)
    real :: electrons       ! Number of electrons on the plates of the capacitor
    real :: energy          ! Energy stored in the eletric field (joules)
    integer :: type         ! Type of input data available for the calculation:
                            !   1: C and V
                            !   2: Charge and V
    real :: v               !   Voltage on the capacitor (volts)

    ! Prompt user for the type of input data available
    write (*,100)
    100 format ('This program calculates information about a ' &
                'capacitor.',/, 'Please specify the type of information' &
                ' available from the following list: ',/, &
                '   1 -- capacitance and voltage ',/, &
                '   2 -- charge and voltage ',/, &
                ' Select options 1 or 2: ')

    ! Get response and validate it
    do
        read (*,*) type
        if ((type == 1) .or. (type == 2)) exit
        write (*,110) type
        110 format ('Invalid response: ', i6, '. Please enter 1 or 2: ')
    end do

    ! Get addicitional data based upon the type of calculation
    input: if (type == 1) then 

        ! Get capacitance
        write (*,*) 'Enter capacitance in farads: '
        read (*,*) c

        ! Get voltage
        write (*,*) 'Enter voltage in volts: '
        read (*,*) v
    else
        
        ! Get charge
        write (*,*) 'Enter charge in coulombs: '
        read (*,*) charge

        ! Get voltage
        write (*,*) 'Enter voltage in volts: '
        read (*,*) v
    end if input

    ! Calculate the unknown quantities
    calculate: if (type == 1) then
        charge = c * v                              ! Charge
    else
        c = charge / v                              ! Capacitance
    end if calculate
    electrons = charge * electrons_per_coulomb      ! Electrons
    energy = 0.5 * c * v**2                         ! Energy

    ! Write out answers
    write (*,120) v, c, charge, electrons, energy
    120 format ('For this capacitor: ' ,/, &
                '  Voltage              = ', f10.2, ' V' ,/, &
                '  Capacitance          = ', es10.3, ' F' ,/, &
                '  Total Charge         = ', es10.3, 'C' ,/, &
                '  Number of electrons  = ', es10.3 ,/, &
                '  Total energy         = ', f10.4, ' joules')

    !

end program capacitor