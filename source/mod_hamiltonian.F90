!> This module contains variables and subroutine to build tthe spin hamiltonian
module mod_hamiltonian
  use mod_types, only: dp
  implicit none
  complex(kind=dp),dimension(:,:), allocatable :: hamilt
  !! Hamiltonian matrix
  real(kind=dp),   dimension(:)  , allocatable :: eigenvalues
  !! Real eigenvectors 
contains

  subroutine allocate_hamiltonian(dim)
    implicit none
    integer, intent(in) :: dim

    allocate( hamilt(dim,dim), eigenvalues(dim) )
  end subroutine allocate_hamiltonian

  !> This subroutine builds the spin hamiltonian
  subroutine build_hamiltonian()
    use mod_constants,  only: g, mub
    use mod_parameters, only: dim, D, E, B
    use mod_spins,      only: Sx, Sy, Sz, Sx2, Sy2, Sz2
    use mod_tools,      only: ItoS
    implicit none
    character(len=30)  :: format_var
    !! Format variable for the Hamiltonian
    integer            :: i,j
    !! Counter for loops

    ! Building the Hamiltonian
    hamilt = D*Sz2 + E*(Sx2-Sy2) + g*mub*( B(1)*Sx + B(2)*Sy + B(3)*Sz)

    ! Checking if Hamiltonian is hermitian
    if( sum(abs(conjg(transpose(hamilt))-hamilt)) > 1.d-12 ) then
      write(*,"('Hamiltonian is not hermitian!')")
      stop
    end if

    ! Writing the hamiltonian to file
    write(unit=11,fmt="('Hamiltonian:')")
    format_var = "(" // trim(itos(dim)) // "(es9.2,'+i',es9.2,3x))"
    do i = 1,dim
      write(unit=11,fmt=trim(format_var)) (real(hamilt(i,j)),aimag(hamilt(i,j)), j=1,dim)
    end do
  end subroutine build_hamiltonian

  subroutine deallocate_hamiltonian()
    implicit none

    deallocate( hamilt,eigenvalues )
  end subroutine deallocate_hamiltonian

end module mod_hamiltonian