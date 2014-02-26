module qc_geom
  use qc_kind
  use qc_constant
  implicit none

  integer :: natom

  real(dp), dimension(:,:), allocatable :: atom_pos
  integer,  dimension(:),   allocatable :: atom_znum

contains

  subroutine qc_geom_read
    integer :: i, j
    character(len=2) :: atname
    real(dp) :: xi, yi, zi
    real(dp), dimension(:), allocatable :: pos_local
    
    if (.TRUE.) then
       open (unit=15, file='geom.xyz', status='old')
       read (15, *) natom
       read (15, *)
    end if

    call qc_geom_mem_init

    allocate (pos_local(natom*3))

    if (.TRUE.) then
       j = 0
       do i = 1, natom
          read (15, *) atname, xi, yi, zi
          j = j + 1 
          pos_local(j) = xi*ang_to_bohr
          j = j + 1 
          pos_local(j) = yi*ang_to_bohr
          j = j + 1 
          pos_local(j) = zi*ang_to_bohr
          atom_znum(i) = atomic_number_get (atname)
       end do
       close (15)
    end if

    j = 0
    do i = 1, natom
       j = j+1
       atom_pos(1,i) = pos_local(j)
       j = j+1
       atom_pos(2,i) = pos_local(j)
       j = j+1
       atom_pos(3,i) = pos_local(j)
    end do

    deallocate (pos_local)

  end subroutine qc_geom_read


  subroutine qc_geom_mem_init

    allocate (atom_pos(3,natom))
    allocate (atom_znum(natom))

  end subroutine qc_geom_mem_init


  subroutine qc_geom_mem_free
    deallocate (atom_pos)
    deallocate (atom_znum)
  end subroutine qc_geom_mem_free

end module qc_geom
