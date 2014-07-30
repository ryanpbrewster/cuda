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
       print *, 'K = ', nw_nbf

       nw_nmo = 0
       read (unitno) (nw_nmo(i), i = 1, nw_nsets)
       print *, 'N = ', nw_nmo(1:nw_nsets)
    
       allocate (occ(nw_nbf))
    end if

    
    allocate (nw_en(nw_nbf))
    allocate (nw_co(nw_nbf,nw_nbf))
    !---
    ! only close system
    !do iset = 1, nsets
    !
    read (unitno) (occ(jb), jb = 1, nw_nbf)
    read (unitno) (nw_en(jb), jb = 1, nw_nbf)
    do im = 1, nw_nmo(1)
       read (unitno) (nw_co(jb,im), jb = 1, nw_nbf)
    end do
    close (unitno)

    do i = 1, nw_nbf
        print *, 'e(', i, ') = ', nw_en(i)
    end do

    do i = 1, nw_nbf
        do j = 1, nw_nbf
            print *, 'c(', i, ',', j, ') = ', nw_co(j,i)
        end do
    end do

    if (qc_ngfs .ne. nw_nbf) then
       call qc_abort ('You might use the different basis sets or geometry')
    end if

  end subroutine nw_vectors_read


end module nw_vectors
