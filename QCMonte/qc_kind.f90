module qc_kind
  implicit none

  integer, parameter :: dp = selected_real_kind(15)
  
contains

  subroutine qc_kind_test
    real(dp) :: f

    f = 3.145_dp

    print *, precision(f), f
    
    f = 3.145E-4_dp
    print *, precision(f), f

  end subroutine qc_kind_test

end module qc_kind


