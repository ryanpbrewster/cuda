program qc_main
  use qc_step
  use qc_input
  use qc_monte
  use nw_vectors

  implicit none

  print *, 'Starting execution'
  
  call qc_start

  call qc_input_read

  call nw_vectors_read

  print *, 'Done reading inputs'
  
  call qc_monte_energy 

  print *, 'Done with MC'
  
  call qc_input_mem_free
 
  call qc_end

  print *, 'Done executing'

end program qc_main
