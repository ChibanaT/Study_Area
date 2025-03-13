program power

    ! Purpose:
    !   To calculate the current, real, reactive, and apparent power,
    !   and the power factor supplied to a load.

    ! Record of Revisions:
    !       Date                Programmer              Description of Change
    !   ============         ================        ============================
    !    2025/03/04              Chibana                Original code

    implicit none

    ! Data Dictionary: declare constants
    real, parameter :: deg_2_rad = 0.01745329       ! Deg to radians factor

    ! Data Dictionary: declare variable types, definitions, & units
    real :: amps        ! Current in the load (A)
    real :: p           ! Real power of the load (W)
    real :: pf          ! Power factor of load (no units)
    real :: q           ! Reactive power of the load (VAR)
    real :: s           ! Apparent power of the load (VA)
    real :: theta       ! Impedance angle of the load (deg)
    real :: volts       ! Rms voltage of the power source (V)
    real :: z           ! Magnitude of the load impedance (ohms)

    ! Prompt the user for the rms voltage
    write (*,*) 'Enter the rms voltage of the source: '
    read (*,*) volts

    ! Prompt the user for the magnitude and angle of the impedance
    write (*,*) 'Enter the magnitude and angle of impedance '
    write (*,*) 'in ohms and degrees: '
    read (*,*) z, theta

    ! Perform calculations
    amps = volts / z
    p = volts * amps * cos (theta * deg_2_rad)
    q = volts * amps * sin (theta * deg_2_rad)
    s = volts * amps
    pf = cos (theta * deg_2_rad)

    ! Write out the results
    write (*,*) 'Voltage        = ', volts, ' volts'
    write (*,*) 'Impedance      = ', z, ' ohms at ', theta, ' degrees'
    write (*,*) 'Current        = ', amps, ' amps'
    write (*,*) 'Real Power     = ', p, ' watts'
    write (*,*) 'Reactive Power = ', q, ' VAR'
    write (*,*) 'Apparent Power = ', s, ' VA'
    write (*,*) 'Power Factor   = ', pf

    ! Finish up
    
end program power