program fun_xy
    ! Purpose:
    !   This program solves the function f(x,y) for a user specified x and y,
    !   where f(x,y) is defined as:

    !                _
    !               |
    !               | X + Y             X >= 0 and Y >= 0
    !               | X + Y**2          X >= 0 and Y < 0
    !   F(X,Y) =    | X**2 + Y          X < 0 and Y >= 0
    !               | X**2 + Y**2       X < 0 and Y < 0
    !               |_
    !            
    ! Record of Revisions:
    !      Date                Programmer              Description of Change
    !   ===========         ================        ============================
    !   2025/03/06              Chibana                 Original code

    implicit none

    ! Data Dictionary: declare variable types, definitions, & units
    real :: x           ! First independent variable
    real :: y           ! Second independent variable
    real :: fun         ! Resulting function

    ! Prompt the user for the values x and y
    write (*,*) 'Enter the coefficients x and y: '
    read (*,*) x, y

    ! Write the coefficients of x and y
    write (*,*) 'The coefficients x and y are: ', x, y
    
    ! Calculate the function f(x,y) based upon the signs of x and y
    if ((x >= 0.) .and. (y >= 0.)) then
        fun = x + y
    else if ((x >= 0.) .and. (y < 0.)) then
        fun = x + y**2
    else if ((x < 0.) .and. (y >= 0.)) then
        fun = x**2 + y
    else
        fun = x**2 + y**2
    end if

    ! Write the value of the function
    write (*,*) 'The value of the function is: ', fun

    ! Finish
    
end program fun_xy