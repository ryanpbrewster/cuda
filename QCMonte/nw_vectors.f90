module nw_vectors
  use qc_kind
  use qc_geom
  use qc_basis
  use qc_step
  implicit none
  integer :: nw_nsets
  integer :: nw_nbf
  integer :: nw_nmo(2)
  integer :: nw_iocc
  integer :: nw_icore
  real(dp), dimension(:), allocatable   :: nw_en
  real(dp), dimension(:,:), allocatable :: nw_co

contains

  subroutine nw_vectors_read
    integer :: i, j, im, jb
    integer :: unitno
    integer :: lentit
    integer :: lenbas
    real(dp), dimension(:),   allocatable  :: occ
    character(len=20) :: scftype20
    character(len=80) :: title
    character(len=20) :: basis_name

    if (.TRUE.) then
       unitno = 15
       open (unitno,status='old',form='unformatted', &
            file='nwchem.movecs')
       read (unitno) ! skip convergence info
       read (unitno) scftype20
       read (unitno) lentit

       title = ' '
       if (lentit .le. 80) then
          read (unitno) title(1:lentit)
       else
          print *, 'title name too short', lentit
       end if

       read (unitno) lenbas
       basis_name = ' '
       read (unitno) basis_name (1:lenbas)

       read (unitno) nw_nsets
       read (unitno) nw_nbf
       print *, 'nw_vectors: nbf ', nw_nbf

       nw_nmo = 0
       read (unitno) (nw_nmo(i), i = 1, nw_nsets)
       print *, 'nw_vectors: nmo ', nw_nmo(1:nw_nsets)
    
       allocate (occ(nw_nbf))
    end if

    
    allocate (nw_en(nw_nbf))
    allocate (nw_co(nw_nbf,nw_nbf))
    !---
    ! only close system
    !do iset = 1, nsets
    !
    if (.TRUE.) then
       read (unitno) (occ(jb), jb = 1, nw_nbf)
       read (unitno) (nw_en(jb), jb = 1, nw_nbf)
       do im = 1, nw_nmo(1)
!          print *, 'mo i', im, nw_en(im)
          read (unitno) (nw_co(jb,im), jb = 1, nw_nbf)
       end do
       close (unitno)
       !--
       nw_iocc = 0
       do j = 1, nw_nbf
          if (occ(j) .gt. 0.0_dp) then
             nw_iocc = j
          end if
       end do

       deallocate (occ)! evals, kvecs)

       call orbital_check

    end if

    if (qc_ngfs .ne. nw_nbf) then
       call qc_abort ('You might use the different basis sets or geometry')
    end if

  end subroutine nw_vectors_read


  subroutine orbital_check
    integer :: i, j
    integer :: znum, iocc

    j = 0
    do i = 1, natom
       j = j + atom_znum(i)
    end do

    iocc = j/2

    !if (iocc .ne. nw_iocc) then
    !   call qc_abort('iocc error')
    !end if


    nw_icore = 0
    do i = 1, natom
       znum = atom_znum(i)
       if (znum .ge. 3 .and. znum .le. 10) nw_icore = nw_icore + 1
    end do

  end subroutine orbital_check


end module nw_vectors
