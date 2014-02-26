module qc_monte
  use qc_kind
  use qc_constant
  use qc_step
  use qc_geom
  use nw_vectors
  use qc_rnd
  use qc_input, only: mc_ntrial, mc_delx, mc_pair_num, mc_lband
  use qc_psi
  use qc_mcmp2

  implicit none

  real(dp), parameter :: rnd_min = 1.0e-3_dp
  real(dp) :: delx

  !----(mc_basis)
  type mc_basis_typ
     integer :: znum
     real(dp):: alpha(2), norm(2)
  end type mc_basis_typ

  type(mc_basis_typ), dimension(:), allocatable, target :: mc_basis_list
  !integer :: iocc1, iocc2, ivir1, ivir2
  !--
  integer :: mc_nbas
  integer, dimension(:), allocatable :: atom_ibas
  !-----

contains

  subroutine qc_monte_energy 

    iocc1 = nw_icore + 1
    iocc2 = nw_iocc
    ivir1 = nw_iocc + 1
    ivir2 = nw_nmo(1) 

    call rnd_init
    call qc_psi_mem_init
    call mc_basis_read

    call monte_energy

    call mc_basis_mem_free
    call qc_psi_mem_free

  end subroutine qc_monte_energy


  subroutine monte_energy
    integer  :: ip, it, nsucc, nfail, ntrial0
    real(dp) :: emp2, qeps(2), qeps_ave(2)
    real(dp) :: diff, emp2_ave, emp2_var, qeps_var(2)
    real(dp) :: scs_emp2, scs_emp2_ave, scs_emp2_var
    real(dp) :: dcpu, cpu_start, cpu_end
    real(dp) :: ratio, emp2_tmp, var_tmp, qeps_tmp(2), qvar_tmp(2)
    real(dp) :: scs_emp2_tmp, scs_var_tmp
    real(dp) :: sys_nproc
    !-- te_val --

    call mc_te_val_init 

    !--- (Get Normalized Factor of importance sampling P func)---
    call mc_eri2v 

    el_pair_num = mc_pair_num

    call mc_mem_init
    ! --- Initial Positions
    call mc_elec_pos_init

    do ip = 1, el_pair_num
       call mc_weight_func_get (ip) !el_pos1, el_pos2, wgt12)
    end do

    ! Perform Metropolis Monte Carlo

    nsucc = 0
    nfail = 0
    ntrial0 = 100000
    delx    = mc_delx

    do it = 1, ntrial0

       call mc_move_scheme (nsucc, nfail)

       if (mod(it, 1000) .eq. 0) then
          ratio = dble(nfail)/dble(nfail + nsucc)
          if (ratio .lt. 0.5_dp) then
             ratio = max(1.0_dp/(2.0_dp*ratio),1.1_dp)
          else
             ratio = max(0.9_dp, 1.0_dp/(2.0_dp*ratio))
          end if
          delx = delx * ratio !1.1_dp
          print *, 'Iteration', it

          nsucc = 0
          nfail = 0
       end if
    end do

    call qc_step_start ("MONTE_MP2")

    print *, 'Back in monte_energy'

    cpu_start = 0.0_dp
    if (.TRUE.) cpu_start = cputime ()

    nsucc = 0
    nfail = 0

    emp2_ave = 0.0_dp
    emp2_var = 0.0_dp
    qeps_ave = 0.0_dp
    qeps_var = 0.0_dp

    scs_emp2_ave = 0.0_dp
    scs_emp2_var = 0.0_dp
    
    ! --- initialize
    do ip = 1, el_pair_num
       el_pair_list(ip)%is_new = .true.
    end do

    print *, 'Calling mc_local_energy'
    call mc_local_energy (emp2, scs_emp2, qeps)

    print *, 'Starting mc_ntrials'
    do it = 1, mc_ntrial
        if (mod(it, 10000) == 0) then
            print *, 'MC Iteration', it
        endif
            

       call mc_move_scheme (nsucc, nfail)

       call mc_local_energy (emp2, scs_emp2, qeps)
       
       diff = emp2 - emp2_ave
       emp2_ave = emp2_ave + diff/dble(it) 
       emp2_var = emp2_var + diff*diff*dble(it-1)/dble(it)

       diff = scs_emp2 - scs_emp2_ave
       scs_emp2_ave = scs_emp2_ave + diff/dble(it) 
       scs_emp2_var = scs_emp2_var + diff*diff*dble(it-1)/dble(it)

       diff = qeps(1) - qeps_ave(1)
       qeps_ave(1) = qeps_ave(1) + diff/dble(it)
       qeps_var(1) = qeps_var(1) + diff*diff*dble(it-1)/dble(it)
       diff = qeps(2) - qeps_ave(2)
       qeps_ave(2) = qeps_ave(2) + diff/dble(it)
       qeps_var(2) = qeps_var(2) + diff*diff*dble(it-1)/dble(it)


       if (mod(it, 1000) .eq. 0) then
          ratio = dble(nfail)/dble(nfail + nsucc)
          if (ratio .lt. 0.5_dp) then
             ratio = max(1.0_dp/(2.0_dp*ratio),1.1_dp)
          else
             ratio = max(0.9_dp, 1.0_dp/(2.0_dp*ratio))
          end if
          delx = delx * ratio !1.1_dp

          nsucc = 0
          nfail = 0
       end if

       if (mod (it,100) .eq. 0) then

          emp2_tmp = emp2_ave
          scs_emp2_tmp = scs_emp2_ave
          var_tmp = emp2_var
          scs_var_tmp = scs_emp2_var
          qeps_tmp = qeps_ave
          qvar_tmp = qeps_var
          
          sys_nproc = 1.0_dp
          if (.TRUE.) then
             emp2_tmp = emp2_tmp/dble(sys_nproc)
             var_tmp = sqrt(var_tmp/dble(sys_nproc*it*(sys_nproc*it - 1.0_dp)) )

             scs_emp2_tmp = scs_emp2_tmp/dble(sys_nproc)
             scs_var_tmp = sqrt(scs_var_tmp/dble(sys_nproc*it*(sys_nproc*it - 1.0_dp)) )

             cpu_end = cputime ()
             dcpu = cpu_end - cpu_start
             cpu_start = cpu_end
             write(25, '(f10.2,4f14.6,f10.2)') dble(it)/1000_dp, emp2_ave, &
                  sqrt(emp2_var/dble(it*(it-1.0_dp))), emp2_tmp, var_tmp, dcpu
             write(26, '(f10.2,4f14.6,f10.2)') dble(it)/1000_dp, scs_emp2_ave, &
                  sqrt(scs_emp2_var/dble(it*(it-1.0_dp))), &
                  scs_emp2_tmp, scs_var_tmp, dcpu

             if (mc_lband) then
                qeps_tmp = qeps_tmp/dble(sys_nproc)
                qvar_tmp(1) = sqrt(qvar_tmp(1)/dble(sys_nproc*it*(sys_nproc*it - 1.0_dp)) )
                qvar_tmp(2) = sqrt(qvar_tmp(2)/dble(sys_nproc*it*(sys_nproc*it - 1.0_dp)) )
                write(27, '(f10.2,4f14.6)') dble(it)/1000_dp, qeps_ave(1:2), &
                     qeps_tmp(1:2)
                write(28, '(f10.2,4f14.6)') dble(it)/1000_dp, qeps_var(1:2), &
                     qvar_tmp(1:2)
                call flush(27)
                call flush(28)
             end if
             call flush(26)
             call flush(25)
             !call flush(24)
             call flush(23)
             !call flush(22)
          end if
       end if

    end do

    call mc_mem_free
    call qc_step_end

  end subroutine monte_energy


  subroutine mc_mem_init

    call mcmp2_mem_init

  end subroutine mc_mem_init


  subroutine mc_mem_free

    call mcmp2_mem_free

  end subroutine mc_mem_free


  subroutine mc_elec_pos_init
    integer :: ia, ip
    real(dp):: twopi, rval, amp1, amp2, theta1, theta2
    real(dp):: pos_i(3)

    twopi = 2.0_dp * pi

    do ip = 1, el_pair_num
       ia = int(natom*rnd_val()) + 1
       pos_i= atom_pos(1:3,ia)

       rval = rnd_val() * 0.2_dp
       amp1 = sqrt(-0.5_dp*log(rval))
       rval = rnd_val() * 0.5_dp
       amp2 = sqrt(-0.5_dp*log(rval))
       theta1 = twopi* rnd_val()
       theta2 = pi* rnd_val()
       el_pair_list(ip)%pos1(1) = pos_i(1) + amp1*cos(theta1)
       el_pair_list(ip)%pos1(2) = pos_i(2) + amp1*sin(theta1)
       el_pair_list(ip)%pos1(3) = pos_i(3) + amp2*cos(theta2)
       
       !-- elec position 2
       ia = int(natom*rnd_val()) + 1
       pos_i= atom_pos(1:3,ia)

       rval = rnd_val() * 0.2_dp
       amp1 = sqrt(-0.5_dp*log(rval))
       rval = rnd_val() * 0.5_dp
       amp2 = sqrt(-0.5_dp*log(rval))
       theta1 = twopi* rnd_val()
       theta2 = pi* rnd_val()
       el_pair_list(ip)%pos2(1) = pos_i(1) + amp1*cos(theta1)
       el_pair_list(ip)%pos2(2) = pos_i(2) + amp1*sin(theta1)
       el_pair_list(ip)%pos2(3) = pos_i(3) + amp2*cos(theta2)
    end do

  end subroutine mc_elec_pos_init


  subroutine mc_weight_func_get (ip)
    integer, intent(in)  :: ip
    type(el_pair_typ), pointer :: el_pair
    !---
    integer :: ia, ib, ie
    real(dp) :: r1, r2, r12, dr(3)
    real(dp) :: gf1, gf2, azi, normi

    el_pair => el_pair_list(ip)

    dr = el_pair%pos1 - el_pair%pos2
    r12 = sqrt(dot_product(dr, dr))

    !if (r12 .lt. 1.0E-3_dp) print *, 'r12 is too close', r12

    gf1 = 0.0_dp
    gf2 = 0.0_dp
    do ia = 1, natom 
       ib = atom_ibas(ia)

       dr = el_pair%pos1 - atom_pos(1:3,ia)
       r1 = dot_product(dr, dr)

       dr = el_pair%pos2 - atom_pos(1:3,ia)
       r2 = dot_product(dr, dr)
          
       do ie = 1, 2
          azi = mc_basis_list(ib)%alpha(ie)
          normi = mc_basis_list(ib)%norm(ie)

          gf1 = gf1 + exp(-azi*r1)*normi
          gf2 = gf2 + exp(-azi*r2)*normi
       end do
    end do

    el_pair%r12 = r12
    el_pair%wgt = gf1*gf2/r12
    
  end subroutine mc_weight_func_get


  subroutine mc_move_scheme (nsucc, nfail)
    integer, intent(inout) :: nsucc, nfail
    !---
    integer :: ip
    real(dp) :: pos1_old(3), pos2_old(3), wgt_old, r12_old
    real(dp) :: ratio, rnd3(3), rval
    type (el_pair_typ), pointer :: el_pair

    do ip = 1, el_pair_num
       el_pair => el_pair_list(ip)

       pos1_old = el_pair%pos1
       pos2_old = el_pair%pos2
       wgt_old  = el_pair%wgt
       r12_old  = el_pair%r12
       
       rnd3     = rnd_val3()
       el_pair%pos1 = pos1_old + (rnd3 - 0.5_dp)*delx

       rnd3     = rnd_val3()
       el_pair%pos2 = pos2_old + (rnd3 - 0.5_dp)*delx

       call mc_weight_func_get (ip)

       ratio = el_pair%wgt/wgt_old
       rval = rnd_val ()
       if (rval .lt. rnd_min) rval = rnd_min

       el_pair%is_new = .true.
       
       if (ratio .gt. rval) then
          nsucc = nsucc + 1
       else
          el_pair%pos1 = pos1_old
          el_pair%pos2 = pos2_old
          el_pair%wgt  = wgt_old
          el_pair%r12  = r12_old
          el_pair%is_new = .false.
          nfail = nfail + 1
       end if

    end do

  end subroutine mc_move_scheme


  subroutine mc_local_energy (emp2, scs_emp2, qeps)
    real(dp), intent(out):: emp2, scs_emp2, qeps(2)

    integer  :: ip
    type (el_pair_typ), pointer :: el_pair1

    ! (1) Update PSI
    do ip = 1, el_pair_num
       el_pair1 => el_pair_list(ip)

       if (el_pair1%is_new) then
          call qc_psi_get (el_pair1%pos1, el_pair1%psi1)
          call qc_psi_get (el_pair1%pos2, el_pair1%psi2)
       end if
    end do

    ! (2) Calc. MP2
    call mcmp2_local_energy (emp2, scs_emp2, qeps, mc_lband)

  end subroutine mc_local_energy



  subroutine mc_basis_read
    integer :: ib, ia, znum
    real(dp) :: alpha(2), coef(2)
    type(mc_basis_typ), pointer :: ibas
    character(len=2) :: atname

    if (.TRUE.) then
       open(unit=15, file='mc_basis.dat', status='old')
       read (15, *) mc_nbas
    end if

    
    call mc_basis_mem_init (mc_nbas)

    do ib = 1, mc_nbas
       ibas => mc_basis_list(ib)
       if (.TRUE.) then
          read (15, *) atname
          read (15, *) alpha(1), coef(1)
          read (15, *) alpha(2), coef(2)
          znum     = atomic_number_get(atname)
       end if

       ibas%alpha(1) = alpha(1)
       !ibas%norm(1)  = exp(1.5_dp*log(alpha/pi)) * coef
       ibas%norm(1)  = exp(0.75_dp*log(2.0_dp*alpha(1)/pi)) * coef(1)
       ibas%alpha(2) = alpha(2)
       !ibas%norm(2)  = exp(1.5_dp*log(alpha/pi)) * coef
       ibas%norm(2)  = exp(0.75_dp*log(2.0_dp*alpha(2)/pi)) * coef(2)
       ibas%znum     = znum
    end do
       
    do ia = 1, natom 
       atom_ibas(ia) = atom_to_mc_basis(ia)
    end do

    if (.TRUE.) then
       close (15)
    end if

  contains
    
    function atom_to_mc_basis (ia) result (ib)
      integer, intent(in) :: ia
      integer :: ib, jb, znum
      
      znum = atom_znum(ia)
      
      ib = 0
      do jb = 1, mc_nbas
         if (znum .eq. mc_basis_list(jb)%znum) ib = jb
      end do
      
      if (ib .eq. 0) then
         stop
      end if
      
    end function atom_to_mc_basis

  end subroutine mc_basis_read
  

  
  subroutine mc_basis_mem_init (mc_nbas)
    integer, intent(in) :: mc_nbas

    allocate (atom_ibas(natom))
    allocate (mc_basis_list(mc_nbas))

  end subroutine mc_basis_mem_init



  subroutine mc_basis_mem_free

    deallocate (atom_ibas)
    deallocate (mc_basis_list)

  end subroutine mc_basis_mem_free



  subroutine mc_eri2v !(g_wgt)
    real(dp) :: g_wgt
    integer  :: ia, ja, ib, jb, m, ts, ie, je
    real(dp) :: pisub2, eri
    real(dp) :: azi, azj, comz, gzi, normi, normj
    real(dp) :: rr, tt, cc, h
    real(dp) :: dum1, dum2, dum3, dum4, dum5, dum6
    real(dp), dimension(3) :: posi, posj, dr
    type(mc_basis_typ), pointer :: ibas, jbas

    pisub2 = 2.0_dp*sqrt_pi**5

    g_wgt = 0.0_dp
    
    ! --- electron (1) ----
    do ia = 1, natom 

       posi = atom_pos(1:3,ia)

       ib   = atom_ibas(ia) 
       ibas => mc_basis_list(ib)
       ! --- electron (2) ----
       do ja = 1, natom 
          posj = atom_pos(1:3,ja)

          jb   = atom_ibas(ja) 
          jbas => mc_basis_list(jb)

          dr = posi - posj
          rr = dot_product (dr, dr) 

          do ie = 1, 2
             azi   = ibas%alpha(ie)
             normi = ibas%norm (ie)

             do je = 1, 2
                azj   = jbas%alpha(je)
                normj = jbas%norm (je)

                comz = azi + azj
                gzi  = azi*azj/comz 

                tt = gzi*rr
                cc = pisub2/(sqrt(comz)*azi*azj)
         
                ! for only S
                m = 0
                if (tt .lt. tf(m)) then
                   ts = nint (tt*20.0_dp)
                   h  = 0.05_dp*dble(ts) - tt
                   dum1 = igamma(ts,m+6)*h*0.166666666666667_dp + igamma(ts,m+5)
                   dum2 = dum1*h*0.2_dp  + igamma(ts,m+4)
                   dum3 = dum2*h*0.25_dp + igamma(ts,m+3)
                   dum4 = dum3*h/3.0_dp  + igamma(ts,m+2)
                   dum5 = dum4*h*0.5_dp  + igamma(ts,m+1)
                   dum6 = dum5*h         + igamma(ts,m+0)
                   eri  = cc* dum6
                else
                   eri = cc*f1(m)*exp(f2(m)*log(tt))
                end if
                eri = eri*normi*normj
                !print *, ia, ib, eri
                g_wgt = g_wgt + eri
             end do
          end do
       end do
    end do

    if (.TRUE.) print *, 'g_wgt', g_wgt

    do ib = 1, mc_nbas
       ibas => mc_basis_list(ib)
       ibas%norm = ibas%norm/sqrt(g_wgt)
       !mc_basis_list(ib)%norm(2) = mc_basis_list(ib)%norm(2)/sqrt(g_wgt)
    end do

  end subroutine mc_eri2v
  

end module qc_monte
