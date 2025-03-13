program mass_consumption

    ! Purpose:
    !   To calculate the mass necessary to produce energy for a year

    ! Record of Revisions:
    !       Date            Programmer          Description fo Change
    !   ============     ===============     ===========================
    !    2025/03/04          Chibana            Original code

    implicit none

    ! Data Dictionary: declare constants
    real, parameter :: c = 2.9979e8                           ! Speed of Light (m/s)
    real, parameter :: seconds_per_year = 3600 * 24 * 365     ! Convert 1 year to seconds
    ! Data Dictionary: declare variable types, definition, & units
    real :: mass                ! Necessary mass to produce energy     
    real :: energy_consumption  ! Energy consumption in a year (MJ/s)
    
    ! Prompt user for energy consumption in year
    write (*,*) 'Enter the energy consumption in year: '
    read (*,*) energy_consumption

    ! Echo the energy consumption
    write (*,*) 'The Energy Consumption per second is: ', energy_consumption, 'MJ/s'

    ! Perform Calculation
    mass = ((energy_consumption * (10 ** 6)) * seconds_per_year) / (c ** 2)
    
    ! Display Results
    write (*,*) '======================================================================='
    write (*,*) 'The necessary mass per year is: ', mass, ' kg'

    ! Finish
end program mass_consumption