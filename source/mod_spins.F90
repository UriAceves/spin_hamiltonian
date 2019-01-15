!> Module for the spin operator matrices
!> The general definition is used from:
!> https://en.wikipedia.org/wiki/Spin_(physics)#Higher_spins
module mod_spins
  use mod_types,     only: dp
  use mod_constants, only: ci, c0
  implicit none
  complex(kind=dp), dimension(:,:), allocatable :: Sx, Sy, Sz
  !! x, y and z components of the spin operators
  complex(kind=dp), dimension(:,:), allocatable :: Sx2, Sy2, Sz2, S2
  !! Square of the x, y and z components of the spin operators and total S^2
  real(kind=dp),     dimension(:),  allocatable :: Sx_exp,Sy_exp,Sz_exp,S2_exp
  !! Expectation values of Sx, Sy, Sz and S^2 for each eigenstate
contains

  subroutine allocate_spins(dim)
    implicit none
    integer, intent(in) :: dim

    allocate( Sx(dim,dim),Sy(dim,dim),Sz(dim,dim) )
    allocate( Sx2(dim,dim),Sy2(dim,dim),Sz2(dim,dim),S2(dim,dim) )
    allocate( Sx_exp(dim),Sy_exp(dim),Sz_exp(dim),S2_exp(dim) )
  end subroutine allocate_spins

  subroutine build_spin_operators(dim)
    use mod_constants,  only: ci, c0
    use mod_parameters, only: S
    implicit none
    integer, intent(in) :: dim
    integer :: i

    Sx = c0
    Sy = c0
    Sz = c0
    do i=1,dim
      Sz(i,i) = S + 1 - i
      if (i<dim) then
        Sx(i+1,i) = 0.5*sqrt( (s+1)*2*i-i*(i+1) )
        Sx(i,i+1) = Sx(i+1,i)
        Sy(i+1,i) = ci*Sx(i+1,i)
        Sy(i,i+1) = conjg( Sy(i+1,i) )
      end if
    end do
    Sx2 = matmul(Sx,Sx)
    Sy2 = matmul(Sy,Sy)
    Sz2 = matmul(Sz,Sz)
    S2  = Sx2 + Sy2 + Sz2

  end subroutine build_spin_operators

  subroutine calculate_expectation_values(hamilt)
    implicit none
    integer :: i,j,k
    !! Loop integers
    complex(kind=dp),dimension(:,:), allocatable :: hamilt
    !! Hamiltonian matrix


    write(unit=11,fmt="('Expectation values:')")
    Sx_exp = 0.d0
    Sy_exp = 0.d0
    Sz_exp = 0.d0
    S2_exp = 0.d0
    different_eigenvalues_exp: do j = 1,size(hamilt,2)
      do i=1,size(hamilt,1)
        do k=1,size(hamilt,1)
          Sx_exp(j) = Sx_exp(j) + real( conjg(hamilt(k,j))*Sx(k,i)*hamilt(i,j) )
          Sy_exp(j) = Sy_exp(j) + real( conjg(hamilt(k,j))*Sy(k,i)*hamilt(i,j) )
          Sz_exp(j) = Sz_exp(j) + real( conjg(hamilt(k,j))*Sz(k,i)*hamilt(i,j) )
          S2_exp(j) = S2_exp(j) + real( conjg(hamilt(k,j))*S2(k,i)*hamilt(i,j) )
        end do
      end do
      write(unit=11,fmt="('n = ',i0,2x,'Sx = ',es9.2,2x'Sy = ',es9.2,2x'Sz = ',es9.2,2x,'S2 = ',es9.2,2x)") &
                                 j,        real(Sx_exp(j)), real(Sy_exp(j)), real(Sz_exp(j)), real(S2_exp(j))
    end do different_eigenvalues_exp


  end subroutine calculate_expectation_values

  subroutine deallocate_spins()
    implicit none
    deallocate( Sx,Sy,Sz,Sx2,Sy2,Sz2,S2 )
    deallocate( Sx_exp,Sy_exp,Sz_exp,S2_exp )
  end subroutine deallocate_spins

end module mod_spins