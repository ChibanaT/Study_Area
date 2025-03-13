program ball

    ! Purpose:
    !   To calculate distance traveled by a ball thrown at a specified
    !   angle THETA and at a specified velocity V0 from a point on the
    !   surface of th earth, ignoring the effects of air friction and 
    !   the earth' curvature.

    ! Record of Revisions:
    !       Date                Programmer              Description of Change
    !   ============        =================       ============================
    !    2025/03/07              Chibana                Original code

    implicit none

    ! Data Dictionary: declare constants
    real, parameter :: degrees_2_rad = 0.01745329       ! Deg ==> rad conv.
    real, parameter :: gravity = -9.81                  ! Accel. due to gracity (m/s)

    ! Data Dictionary: declare variable types, definitions, & units
    integer :: max_degrees          ! Angle at which the max rng occurs (degrees)
    real :: max_range               ! Maximum range for the ball at vel V0 (meters)
    real :: range                   ! Range of the ball at a particular angle (meters)
    real :: radian                  ! Angle at which the ball was thrown (in radians)
    integer :: theta                ! Angle at which the ball was thrown (in degrees)
    real ::v0                       ! Velocity of the ball in (m/s)

    ! Initialize variables
    max_range = 0.
    max_degrees = 0
    v0 = 20.

    ! Loop over all specified angles

    loop: do theta = 0, 90

        ! Get angle in radians
        radian = real(theta) * degrees_2_rad

        ! Calculate range in meters
        range = (-2. * v0**2 / gravity) * sin(radian) * cos(radian)

        ! Write out the range for this angle
        Write (*,*) 'Theta  = ', theta, ' degrees; Range = ', range, &
                    ' meters'

        ! Compare the range to the previous maxium range. If this
        ! range is larger, save it and the angle at which i occurred.
        if (range > max_range) then
            max_range = range
            max_degrees = theta
        
        end if
    end do loop

    ! Skip a line, and then write out the maximum range and the angle
    ! at which it occurred
    write (*,*) ' '
    write (*,*) 'Max range = ', max_range, ' at ', max_degrees, ' degrees'

    ! Finish

end program ball