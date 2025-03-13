program stats_3

    ! Purpose:
    !   To calculate mean and the standard deviation of an input
    !   data set, where each input value can be positive, negative,
    !   or zero

    ! Record of Revisions:
    !       Date                Programmer              Description of Change
    !   ===========         ================        ============================
    !   2025/03/06              Chibana                 Original code

    implicit none

    ! Data Dictionary: declare variable types, definitions, & units
    integer :: i                ! Loop index
    integer :: n = 0            ! The number of input samples
    real :: std_dev             ! The standard deviation of the input samples
    real :: sum_x = 0.          ! The sum of the input values
    real :: sum_x2  = 0.        ! The sum of the squares of the input values
    real :: x = 0.              ! An input data value
    real :: x_bar               ! The average of the input samples

    ! Get the number of points to input
    write (*,*) 'Enter number of points: '
    read (*,*) n

    ! Check to see if we have enough input data
    if (n < 2) then ! Insufficient data
        write (*,*) 'At least 2 values must be entered.'

    else ! We will have enought data, so let's get it
        
        ! Loop to read input values
        do i = 1, n

            ! Read values
            write (*,*) 'Enter number: '
            read (*,*) x
            write (*,*) 'The number is ', x

            ! Accumulate sums
            sum_x = sum_x + x
            sum_x2 = sum_x2 + x**2

        end do

        ! Now calculate statistics
        x_bar = sum_x / real(n)
        std_dev = sqrt((real(n) * sum_x2 - sum_x**2) / (real(n) * real(n-1)))

        ! Tell user
        write (*,*) 'The mean of this data set is: ', x_bar
        write (*,*) 'The standard deviation is: ', std_dev
        write (*,*) 'The number of data points is: ', n

    end if

    ! Finish
    
end program stats_3