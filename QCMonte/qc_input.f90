module qc_input
  use qc_kind
  use qc_constant, only:ang_to_bohr
  use qc_geom
  use qc_basis
!  use qc_kvc, only: kvc, qc_kvc_mem_init, qc_kvc_mem_free
  
  implicit none

  character(len=20) :: job_name
  integer :: job_type
  integer :: job_theory
  integer :: charge
!  integer :: ncel_overlap
!  integer :: ncel_coulomb
!  integer :: scf_cycle, scf_diis
!  real(dp):: period
!  real(dp):: eri4_tol
!  real(dp):: scf_conv, scf_relax
  !--- 
  integer, parameter :: job_type_ener  = 0
  !---
  integer, parameter :: job_theory_hf   = 0
  integer, parameter :: job_theory_mp2  = 1
  integer, parameter :: job_theory_monte = 2
  integer, parameter :: job_theory_quad1 = 3
  integer, parameter :: job_theory_quad2 = 4
  integer, parameter :: job_theory_mcmp2 = 5
  integer, parameter :: job_theory_mcmp3 = 6
  !--
  integer :: mc_ntrial, mc_pair_num
  real(dp):: mc_delx
  logical :: mc_lband
contains

  subroutine qc_input_read

    call input_read
    call qc_geom_read
    call qc_basis_read

  end subroutine qc_input_read


  subroutine input_read
    character(len=12) :: key, cval
   
    ! default values
    if (.TRUE.) then
       open (15, file='input', status='old')
       read (15, *) key, job_name
       !--- JOBTYPE ---
       read (15, *) key, cval
       if (cval(1:4) .eq. 'ENER')  job_type= job_type_ener
       !--- THEORY---
       read (15, *) key, cval
       job_theory = -1
       if (cval(1:2) .eq. 'HF')    job_theory = job_theory_hf
       if (cval(1:3) .eq. 'MP2')   job_theory = job_theory_mp2
       if (cval(1:5) .eq. 'QUAD1')  job_theory = job_theory_quad1
       if (cval(1:5) .eq. 'QUAD2')  job_theory = job_theory_quad2
       if (cval(1:5) .eq. 'MONTE') job_theory = job_theory_monte
       if (cval(1:5) .eq. 'MCMP2') job_theory = job_theory_mcmp2
       if (cval(1:5) .eq. 'MCMP3') job_theory = job_theory_mcmp3
       if (job_theory .eq. -1) then
          print *, 'JOB_THEORY IS NOT ASSIGNED'
          stop
       end if
       read (15, *) key, charge
       ! Spherical Basis set
       read (15, *) key, lspherical
       !---
       if (job_theory .eq. job_theory_monte .or. &
            job_theory .eq. job_theory_mcmp2) then
          read (15, *) key, mc_ntrial
          read (15, *) key, mc_pair_num
          read (15, *) key, mc_delx
          read (15, *) key, mc_lband
       end if
       
       close (15)
    end if

    return

  end subroutine input_read



  subroutine qc_input_mem_free

    call qc_basis_mem_free
    call qc_geom_mem_free

  end subroutine qc_input_mem_free


end module qc_input
