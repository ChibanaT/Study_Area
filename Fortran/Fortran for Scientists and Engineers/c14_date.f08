program c14_date

    ! Purpose:
    !   To calculate the age of an organic sample from the percentage
    !   of the original carbon 14 remaining in the sample

    ! Record of Revisions:
    !       Date               Programmer              Description of Change
    !  ==============       ================        ===========================
    !   2025/03/04              Chibana                 Original code

    implicit none

    ! Data Dictionary: declare constants
    real, parameter :: lambda = 0.00012097      ! The radioactive decay
                                                ! constant of carbon 14,
                                                ! in units of 1/year.

    ! Data Dictionary: declare variables types, definitions, & units
    real :: age             ! The age of the sample (years)
    real :: percent         ! The percentage of carbon 14 remaining at the time
                            ! of the measurement (%)
    real :: ratio           ! The ratio of carbon 14 remaining at the time
                            ! of the measurement to the original amount of
                            ! carbon 14 (no units)

    ! Prompt the user for percentage of C-14 remaining
    write (*,*) 'Enter the percentage of carbon 14 remaining: '
    read (*,*) percent

    ! Echo the user's input value
    write (*,*) 'The remaining carbon 14 = ', percent, ' %.'

    ! Perform calculation
    ratio = percent / 100.                  ! Convert to fractional ratio
    age = (-1.0 / lambda) * log(ratio)      ! Get age in years

    ! Tell the user about the age of the sample
    write (*,*) 'The age of the sample is ', age, ' years.'

    !Finish up

end program c14_date