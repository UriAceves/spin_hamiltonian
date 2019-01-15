!> This module contains useful constants
module mod_constants
  use mod_types, only: dp
  implicit none
  complex(kind=dp),parameter   :: ci = (0.d0,1.d0)
  !! Imaginary unit
  complex(kind=dp),parameter   :: c0 = (0.d0,0.d0)
  !! Imaginary zero
  complex(kind=dp),parameter   :: c1 = (1.d0,0.d0)
  !! Imaginary one
  real(kind=dp), parameter     :: mub = 5.7883818d-2
  !! Bohr magneton
  real(kind=dp), parameter     :: g = 2.11d0
  !! g-factor
end module mod_constants