program doy

    ! Purpose
    !   This program calculates the day of year corresponding to a
    !   specified date. It illustrates the use of counting loops
    !   and the SELECT CASE construct

    ! Record of Revisions:
    !       Date                Programmer              Description of Change
    !   ===========         ================        ============================
    !   2025/03/06               Chibana                Original code
    
    implicit none

    ! Data Dictionary: declare variable types, definitions, & units
    integer :: day                  ! Day(dd)
    integer :: day_of_year          ! Day of year
    integer :: i                    ! Index variable
    integer :: leap_day             ! Extra day for leap year
    integer :: month                ! Month (mm)
    integer :: year                 ! Year (yyyy)

    ! Get day, month and year to convert
    write (*,*) 'This program calculates the day of year given the '
    write (*,*) 'current date. Enter current month (1-12), day (1-31),'
    write (*,*) 'and year in that order: '
    read (*,*) month, day, year

    ! Check for leap year, and add extra day if necessary
    if (mod(year, 400) == 0) then
        leap_day = 1                ! Years divisible by 400 are leap years
    else if (mod(year, 100) == 0) then
        leap_day = 0                ! Other centuries are not leap years
    else if (mod(year, 4) == 0) then
        leap_day = 1                ! Otherwise every 4th year is a leap year
    else 
        leap_day = 0                ! Other years are not leap years
    end if

    ! Calculate day of year
    day_of_year = day
    do i = 1, month - 1

        ! Add days in months from January to last month
        select case (i)
        case (1, 3, 5, 7, 8, 10, 12)
            day_of_year = day_of_year + 31
        case (4, 6, 9, 11)
            day_of_year = day_of_year + 30
        case (2)
            day_of_year = day_of_year + 28 + leap_day
        end select

    end do

    ! Tell user 
    write (*,*) 'Day                = ', day
    write (*,*) 'Month              = ', month
    write (*,*) 'Year               = ', year
    write (*,*) 'Day of year        = ', day_of_year

    ! Finish
    
end program doy