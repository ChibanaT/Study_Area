program stats_1

    ! Purpose:
    !   To calculate mean and standard deviation of an input
    !   data set containing an arbitrary number of input values

    ! Record of Revisions:
    !       Date                Programmer              Description of Change
    !   ============        =================       ============================
    !    2025/03/06              Chibana                Original code

    implicit none

    ! Data dictionary: declare variables types, definitions, & units
    integer :: n = 0                ! The number of input samples
    real :: std_dev = 0.            ! The standad deviation of the input samples
    real :: sum_x = 0.              ! The sum of the input values
    real :: sum_x2 = 0.             ! The sum of the squares of the unput values
    real :: x = 0.                  ! An input data value
    real :: x_bar                   ! The average of the input samples

    ! While loop to read input values.
    do
        ! Read in next value
        write (*,*) 'Enter number: '
        read (*,*) x
        write (*,*) 'The number is ', x

        ! Test for loop exit
        if (x < 0) exit

        ! Otherwise, accumulate sums
        n       = n + 1
        sum_x   = sum_x + x
        sum_x2  = sum_x2 + x**2
    end do

    ! Calculate the mean and standard deviation
    x_bar = sum_x / real(n)
    std_dev = sqrt((real(n) * sum_x2 - sum_x**2) / (real(n) * real(n - 1)))

    ! Tell user
    write (*,*) 'The mean of this data set is: ', x_bar
    write (*,*) 'The standard deviation is: ', std_dev
    write (*,*) 'The number of data points is: ', n

    ! Finish
    
end program stats_1