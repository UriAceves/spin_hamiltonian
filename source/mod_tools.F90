!> This module contains useful tools
module mod_tools
  implicit none

contains

  !> This function converts an integer to a string (character)
  character(len=100) function ItoS(i)
    implicit none
    integer, intent(in) :: i

    write(ItoS, "(i0)") i
  end function ItoS
end module tools