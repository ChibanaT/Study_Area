program roots

    ! Purpose:
    !   This program solves the roots of a quadratic equation of the
    !   form a*x**2 + b*x + c = 0. It calculates the answers regardless
    !   of the type of roots that the equation possesses.

    ! Record of Revisions:
    !       Date            Programmer              Description of Change
    !   ===========     =================       =============================
    !   2025/03/04           Chibana                Original code

    implicit none

    ! Data dictionary: declare variable types, definition, & units
    real :: a                   ! Coefficient of x**2 term of equation
    real :: b                   ! Coefficient of x term of equation
    real :: c                   ! Constant term of equation
    real :: discriminant        ! Discriminant of the equation
    real :: imag_part           ! Imaginary part of equation (for complex roots)
    real :: real_part           ! Real part of equation (for complex roots)
    real :: x1                  ! First solution of equation (for real roots)
    real :: x2                  ! Second solution of equation (for real roots)

    ! Prompt the user for the coefficients of the equation
    write (*,*) 'This program solves for the roots of a quadratic '
    write (*,*) ' equation of the form A * X ** 2 + B * X + C = 0.'
    write (*,*) 'Enter the coefficients A, B and C: '
    read (*,*) a, b, c

    ! Echo back coefficients
    write (*,*) 'The coefficients A, B and C are: ', a, b, c

    ! Calculate discriminant
    discriminant = b**2 - 4. * a * c

    ! Solve for the roots, depending upon the value of the discriminant
    if ( discriminant > 0. ) then        ! There are two real roots, so...
        x1 = (-b + sqrt(discriminant)) / (2. * a)
        x2 = (-b - sqrt(discriminant)) / (2. * a)

        write (*,*) 'This equation has two real roots: '
        write (*,*) 'X1 = ', x1
        write (*,*) 'X2 = ', x2

    else if (discriminant < 0.) then
        real_part = (-b) / (2. * a)
        imag_part = sqrt(abs(discriminant)) / (2. * a)

        write (*,*) 'This equation has complex roots: '
        write (*,*) 'X1 = ', real_part, ' +i ', imag_part
        write (*,*) 'X2 = ', real_part, ' -i ', imag_part

    else if (discriminant == 0.) then       ! There is one repeated root, so...
        x1 = (-b) / (2. * a)

        write (*,*) 'This equation has two identical real roots: '
        write (*,*) 'X1 = X2 = ', x1


    end if
    
    
end program roots