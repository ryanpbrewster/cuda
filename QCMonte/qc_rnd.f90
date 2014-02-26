module qc_rnd
  use qc_kind
  implicit none
  
contains

  subroutine rnd_init
    integer, allocatable :: seed(:)
    integer :: my_seed
    integer :: the_size
    integer :: j

    ! seed should be set to a large odd integer
    my_seed = 7654321
    call random_seed (size = the_size)
    allocate (seed(the_size))
    do j = 1, the_size
       seed(j) = my_seed + (j-1)
    end do
    call random_seed (put = seed)

  end subroutine rnd_init


  function rnd_val () result (ran1)
    real(dp) :: ran1
    
    call random_number (ran1)

  end function rnd_val


  function rnd_val3 () result (ran3)
    real(dp),dimension(3) :: ran3
    
    call random_number (ran3)

  end function rnd_val3
  
end module qc_rnd


