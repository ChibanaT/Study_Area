program total_energy
    
    ! Purpose:
    !   To calculate total energy possessed by an object
    !   in the Earth's gravitational field.

    ! Record of Revisions:
    !           Date               Programmer              Description of Change
    !       ===========         ================         ===========================
    !       2025/03/04              Chibana                 Original code

    implicit none
    
    ! Data Dictionary: declare constants
    real, parameter :: g = 9.80665         ! Gravity acceleration

    ! Data Dictionary: declare variable types, definitions, & units
    real :: energy                      ! Total energy (J)
    real :: potential_energy            ! Potential energy (J)
    real :: kinectic_energy             ! Kinectic energy (J)
    real :: mass                        ! Mass of the body (kg)
    real :: height                      ! Height of the body above the surface of the Earth (m)
    real :: velocity                    ! Velocity of the body (m/s)

    ! Prompt user for variables
    write (*,*) 'Enter the object mass: '
    read (*,*) mass
    write (*,*) 'Enter the object velocity: '
    read (*,*) velocity
    write (*,*) 'Enter the object height above surface of the Earth: '
    read (*,*) height

    ! Echo user the variables
    write (*,*) 'The object mass is:                              ', mass, ' kg'
    write (*,*) 'The object velocity is:                          ', velocity, ' m/s'
    write (*,*) 'The object height above surface of the Earth is: ', height, ' m'

    ! Perform calculation
    potential_energy = mass * height * g
    kinectic_energy = 0.5 * mass * (velocity ** 2)
    energy = potential_energy + kinectic_energy

    ! Display the results
    write(*,*) '======================================================================'
    write (*,*) 'The Potential Energy is:                   ', potential_energy, ' J'
    write (*,*) 'The Kinematic Energy is:                   ', kinectic_energy, ' J'
    write (*,*) 'The Total Energy is:                       ', energy, ' J'

    ! Finish

end program total_energy
