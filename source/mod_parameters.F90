!> This module contains parameters for the system
module mod_parameters
  use mod_types, only: sp, dp
  implicit none
  real(kind=sp)      :: S
  !! Spin quantum number
  integer            :: dim
  !! Dimension of the Hamiltonian: 2S+1
  real(kind=dp)      :: D, E, B(3) = 0.d0
  !! Hamiltonian parameters
  integer            :: npoints = 1
  !! Number of points in the magnetig field loop
  integer            :: comp = 3
  !! Component that will be varied
  real(kind=dp)      :: Bmax
  !! Maximum value of the magnetic field
  namelist /parameters/ S, D, E, B 
  namelist /Bloop/ Bmax, npoints, comp
end module mod_parameters