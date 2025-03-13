program compare

    ! Purpose:
    !   To compare two strings to see if they are equivalent,
    !   ignoring case.

    ! Record of Revisions:
    !       Date               Programmer              Description of Change
    !   ===========         ================        =============================
    !   2025/03/07              Chibana                 Original code

    implicit none

    ! Data dictionary: declare variable types, definitions, & units
    integer :: i                        ! Loop index
    character(len=20) :: str1           ! First string to compare
    character(len=20) :: str1a          ! Copy of first string to compare
    character(len=20) :: str2           ! Second string to compare
    character(len=20) :: str2a          ! Copy of second string to compare
    
    ! Prompt for the strings
    write (*,*) 'Enter first string to compare: '
    read (*,*) str1
    write (*,*) 'Enter second string to compare '
    read (*,*) str2

    ! Make copies so that the original strings are not modified
    str1a = str1
    str2a = str2

    ! Now shift lowercase letters to uppercase
    do i = 1, len(str1a)
        if (str1a(i : i) >= 'a' .and. str1a(i : i) <= 'z') then
            str1a(i : i) = achar (iachar (str1a(i : i)) - 32)
        end if
    end do

    do i = 1, len(str2a)
        if (str2a(i : i) >= 'a' .and. str2a(i : i) <= 'z') then
            str2a(i : i) = achar (iachar (str2a(i : i)) -32)
        end if
    end do

    ! Compare strings and write result
    if (str1a == str2a) then
        write (*,*) " ' ", str1, "' = '", str2, "' ignoring case."
    else
        write (*,*) " ' ", str1, "' /= ", str2, "' ignoring case."
    end if

    ! Finish

end program compare