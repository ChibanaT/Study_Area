program table

    ! Purpose:
    !   To illustrate the use of formatted WRITE statements. This
    !   program generates a table containing the square roots, squares,
    !   and cubes of all integers betwenn 1 and 10. The table includes
    !   a title and column headings.

    ! Record of Revisions:
    !       Date                Programmer              Description of Change
    !   ============         ================        ===========================
    !    2025/03/07              Chibana                Original code

    implicit none

    integer :: cube         ! The cube of i
    integer :: i            ! Index variable
    integer :: square       ! The square of i
    real :: square_root     ! The square root of i

    ! Print the title of the table on a new page
    write (*,100)
    100 format (t3, 'Table of Square Roots, Squares, and Cubes')

    ! Print the column headings after skipping one line
    write (*,110)
    110 format (t4, 'Number', t13, 'Square Root', t29, 'Square', t39, 'Cube')
    write (*,120)
    120 format (t4, '========', t13, '============', t29, '========', t39, '========')

    ! Generate the required values, and print them out
    do i = 1, 10
        square_root = sqrt(real(i))
        square = i**2
        cube = i**3
        write (*,130) i, square_root, square, cube
        130 format (1x, t4, i4, t13, f10.6, t27, i6, t37, i6)
    end do

    ! Finish
    
end program table