program free_fall

    ! Purpose:
    !   To calculate the velocity of impact of a free fall body

    ! Record of Revisions:
    !          Date                Programmer              Description of Change
    !       ===========         ================        ============================
    !       2025/03/04              Chibana              Original code

    implicit none

    ! Data Dictionary: declare constants
    real, parameter :: g = 9.80665          ! Gravity acceleration (m/s^2)

    ! Data Dictionary: declare variable types, definitions, & units
    real :: velocity                ! Impact velocity (m/s)
    real :: height                  ! Body height from the impact surface

    ! Prompt user for variables
    write (*,*) 'Enter the height of the body from the impact surface: '
    read (*,*) height

    ! Echo the variables
    write (*,*) 'The body is ', height, ' meters from the surface.'

    ! Perform Calculations
    velocity = (2 * g * height) ** 0.5

    ! Display Results
    write (*,*) '====================================================================='
    write (*,*) 'The impact velocity is: ', velocity, ' m/s'

    ! Finish

end program free_fall