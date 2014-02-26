program nw_vector
  implicit none
  integer :: i, j, im, iset
  integer :: unitno
  integer :: lentit
  integer :: lenbas
  integer :: nsets
  integer :: nbf
  integer :: nmo(2)
  real*8, dimension(:),   allocatable  :: occ, evals
  real*8, dimension(:,:), allocatable  :: kvecs
  character(len=20) :: scftype20
  character(len=80) :: title
  character(len=20) :: basis_name

  unitno = 15
  open (unitno,status='old',form='unformatted', &
       file='h2o.movecs')
  read (unitno) ! skip convergence info
  read (unitno) scftype20
  print *, 'scftype20 ', scftype20
  read (unitno) lentit
  print *, 'lentit', lentit

  title = ' '
  if (lentit .le. 80) then
     read (unitno) title(1:lentit)
     print *, 'title ', title
  else
     print *, 'title name too short', lentit
  end if
  
  read (unitno) lenbas
  print *, 'lenbas ', lenbas
  basis_name = ' '
  read (unitno) basis_name (1:lenbas)
  print *, 'basis_name ', basis_name
  
  read (unitno) nsets
  print *, 'nsets ', nsets
  read (unitno) nbf
  print *, 'nbf ', nbf
  
  nmo = 0
  read (unitno) (nmo(i), i = 1, nsets)
  print *, 'nmo(i) ', nmo(1:nsets)

  allocate (occ(nbf))
  allocate (evals(nbf))
  allocate (kvecs(nbf,nbf))
  !---
  do iset = 1, nsets
     read (unitno) (occ(j), j = 1, nbf)
     print *, 'occ(j) ', occ(1:nbf)
     read (unitno) (evals(j), j = 1, nbf)
     print *, 'evals(j) ', evals(1:nbf)
     do im = 1, nmo(iset)
        read (unitno) (kvecs(j,im), j = 1, nbf)
        print *, 'kvecs(i,j) ', im, kvecs(1:nbf,im)
     end do
  end do

  do j = 1, nbf
     if (occ(j) .gt. 0.0d0) then
       print *, j
     end if
  end do
  deallocate (occ, evals, kvecs)

  close (unitno)

end program nw_vector
