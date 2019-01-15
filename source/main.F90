!> This program mounts the spin Hamiltonian:
!> H = DS_z^2 + E(S_x^2 - S_y^2) + \vec{B}\cdot\vec{S}
!> In the basis of the eigenvectors of S^2 and S_z, |S M>
!> and diagonalizes it.
program spin_hamiltonian
  use mod_types,       only: sp, dp
  use mod_constants,   only: ci, c0
  use mod_parameters,  only: S, dim, parameters, Bloop, comp, Bmax, npoints, B
  use mod_spins,       only: allocate_spins, build_spin_operators, calculate_expectation_values, deallocate_spins
  use mod_hamiltonian, only: hamilt, eigenvalues, allocate_hamiltonian, build_hamiltonian, deallocate_hamiltonian
  use mod_tools,       only: ItoS
  implicit none
  integer            :: i,j,k,l
  !! Counter for loops
  integer            :: ios
  !! Error variable (for I/O)
  character(len=10)  :: input_file = "input", output_file = "output", results_file = "results"
  !! Input and output filenames
  character(len=30)  :: format_var
  !! Format variable for the hamiltonian
  real(kind=dp)      :: Bmin
  !! Minimum value for the magnetic field (from B(comp))
  real(kind=dp)      :: Bstep
  !! Step size between different fields
  
  ! Reading spin and Hamiltonian parameters from input file
  open(unit=10,file=input_file)
  read(unit=10,nml=parameters,iostat=ios)
  if(ios /= 0) then
    write(*,"('Error reading ""parameter"" list from ',a)") input_file
    stop
  end if
  read(unit=10,nml=Bloop,iostat=ios)
  if(ios /= 0) then
    write(*,"('Error reading ""Bloop"" list from ',a)") input_file
    stop
  end if
  close(unit=10)

  ! Checking if spin is integer or half-integer and non-zero
  if(( dble(nint(2*S))-2*S > 1.d-12 ).or.( S == 0.0 )) then
    write(*,"('S is not an integer or half-integer. S = ',es9.2)") S
    stop
  end if

  ! Dimension of the Hamiltonian
  dim = 2*S+1

  ! Allocating spin variables and building the spin operators
  call allocate_spins(dim)
  call build_spin_operators(dim)

  ! Allocating hamiltonian and eigenvalues
  call allocate_hamiltonian(dim)

  ! Writing the parameters to output file
  open(unit=11,file=output_file)
  write(unit=11,fmt="('Input parameters:')")
  write(unit=11,nml=parameters)

  ! Opening results file and writing header
  open(unit=12,file=results_file)
  ! Print the Git version (VERSION is defined via CMake macro and defined with compiler flag -DVERSION='')
#if defined(VERSION)
  write(unit=11,fmt="('# Git version: ',a)") VERSION
  write(unit=12,fmt="('# Git version: ',a)") VERSION
#else
  write(unit=11,fmt="('# Git version: unknown')")
  write(unit=12,fmt="('# Git version: unknown')")
#endif

  write(unit=12,fmt="('#       B(',i0,')     ,  Transition energies(2:dim) ')") comp

  ! Looping over magnetic fields (starting with value given in 'parameters')
  Bmin = B(comp)
  Bstep = (Bmax - Bmin)/npoints
  mag_field: do k=1,npoints
    B(comp) = Bmin + (k-1)*Bstep

    call build_hamiltonian()

    ! Diagonalizing the hamiltonian
    call eigensolver( dim, hamilt, eigenvalues )

    ! Writing eigenvectors and eigenvalues
    format_var = "(" // trim(itos(dim)) // "(es9.2,14x))"
    write(unit=11,fmt="('Eigenvalues:')")
    write(unit=11,fmt=trim(format_var)) (eigenvalues(j), j=1,dim)

    format_var = "(" // trim(itos(dim)) // "(es9.2,'+i',es9.2,3x))"
    write(unit=11,fmt="('Eigenvectors:')")
    do i = 1,dim
      write(unit=11,fmt=trim(format_var)) (real(hamilt(i,j)),aimag(hamilt(i,j)), j=1,dim)
    end do

    ! Calculating expectation values of S^2 and S_z
    call calculate_expectation_values(hamilt)

    ! Calculating transition energies from the ground state
    write(unit=11,fmt="('Transition energies from the ground state:')")
    do j = 2,dim
      write(unit=11,fmt="('0 -> ',i0,3x,'\Delta E = ',es9.2)") j-1,eigenvalues(j)-eigenvalues(1)
    end do
    format_var = "(" // trim(itos(dim)) // "(es16.9,2x))"
    write(unit=12,fmt=format_var) B(comp),(eigenvalues(j)-eigenvalues(1), j=2,dim)
  end do mag_field


  ! Closing files
  close(unit=11)
  close(unit=12)

  ! Deallocating variables
  call deallocate_spins()
  call deallocate_hamiltonian()

end program spin_hamiltonian