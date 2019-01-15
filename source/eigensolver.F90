!> This is an interface to the LAPACK eigensolver
subroutine eigensolver(n,a,w)
  implicit none
  integer,          intent(in)    :: n
  double complex,   intent(inout) :: a(n,n)
  double precision, intent(out)   :: w(n)
  ! Workspace variables
  integer                       :: lwork, info
  double complex, allocatable   :: work(:)
  double precision, allocatable :: rwork(:)

  lwork = 3*n
  allocate(work(lwork),rwork(lwork))
  call zheev('V','U',n,a,n,w,work,lwork,rwork,info)
  if (info /= 0) write(*,'("eigensolver: failure in zheev")')
  deallocate(work,rwork)

end subroutine eigensolver