program read_file

    ! Purpose:
    !   To illustrate how to read an unknown number of values from
    !   an input data file, detecting both any formatting erros and
    !   the end of the file.

    ! Record of Revisions:
    !   Date       Programmer      Description of Change 
    !   ====       ==========      =====================
    ! 2025/03/25    Chibana         Original code

    implicit none

    ! Data Dictionary: declare variable types, definitions, & units
    character(len=20) :: filename       ! Nme of file to open
    character(len=80) :: msg            ! Error Message
    integer :: nvals = 0                ! Number of values read in
    integer :: status                   ! I/O status
    real :: value                       ! The real value read in

    ! Get the file name, and echo it back to the user
    write(*,*) 'Please enter input file name: '
    read(*,*) filename
    write (*,1000) filename
    1000 format('Input file name is: ',A)

    ! Open the file, and check for errors on open.
    open(unit = 3, file = filename, status = 'old', action = 'read', &
         iostat = status, iomsg = msg)
    openif : if (status /= 0) then
        
        ! Open was ok. Read values.
        readloop : do 
            read(3,*, iostat = status) value    ! Get next value
            if (status /= 0) exit               ! Exit ifnot valid
            nvals = nvals + 1                   ! Valid: increase count
            write
    
end program read_file