program qc_main
  use qc_step
  use qc_input
  use qc_monte
  use nw_vectors
  use iso_c_binding

  implicit none

  print *, 'Starting execution'

  call qc_start

  call qc_input_read

  call nw_vectors_read

  call qc_basis_print

  call qc_input_mem_free

  call qc_end

  print *, 'Done executing'

end program qc_main
