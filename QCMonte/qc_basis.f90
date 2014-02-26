module qc_basis
  use qc_kind
  use qc_step
  use qc_constant
  use qc_geom
  implicit none

  type qc_shl_typ
     integer :: am
     integer :: iat
     integer :: isgs  ! if lspherical = true, isph is the index to spherical cgs
     integer :: icgs  ! cgs_list(21)
     integer :: nprim
     integer :: ncgs 
     integer :: nsgs 
     !integer,  dimension(:),   allocatable :: csh   ! cartesian shell info
     !integer,  dimension(:),   allocatable :: ssh   ! spherical shell info
     real(dp), dimension(:),   allocatable :: alpha !
     real(dp), dimension(:,:), allocatable :: norm  !coef
  end type qc_shl_typ

  integer :: qc_ncgs
  integer :: qc_ngfs
  integer :: qc_nshl
  logical :: lspherical
  type(qc_shl_typ), dimension(:), allocatable, target :: qc_shl_list

contains

  subroutine qc_basis_read
    integer :: ia, znum, nshell
    integer :: ncgs0, nsgs0, npgs0
    integer :: ncgs, nsgs, npgs, nshl
    integer :: eof
    integer :: is, ip
    character(len=2) :: atname, sym
    integer, dimension(20) :: am, nprim
    real(dp), dimension(20,20) :: alpha, coef, coef2
    logical :: lread, lok

    lread = .false.
    lok   = .false.
    if (.TRUE.) then
       inquire (file='basis.dat', exist=lok)
       
       if (lok) then

          open(unit=15, file='basis.dat', status='old')
          
          lread = .true.
          ! Gaussian94 format
          ! S, SP, P, D
          ncgs = 0
          nsgs = 0
          npgs = 0
          nshl = 0
          do ia = 1, natom
             rewind (15)
             do 
                read (15, '(A2,I5)', iostat=eof) atname, nshell
                if (eof /= 0) then
                   lread = .false.
                end if
                znum = atomic_number_get (atname)
                if (znum == atom_znum(ia)) then
                   call basis_read1 (nshell, ncgs0, nsgs0, npgs0)
                   nshl = nshl + nshell
                   ncgs = ncgs + ncgs0
                   nsgs = nsgs + nsgs0
                   npgs = npgs + npgs0
                   exit
                end if
             end do
          end do
       end if
    end if
    
    if (.not. lok) then
       call qc_abort ('QC_BASIS_READ: FAIL TO OPEN [basis.dat]')
    end if
    if (.not. lread) then
       call qc_abort ('QC_BASIS_READ: FAIL TO READ [basis.dat]')
    end if

    qc_ncgs = ncgs
    qc_ngfs = ncgs
    if (lspherical) qc_ngfs = nsgs
    qc_nshl = nshl

    call qc_basis_mem_init

    nshl = 0
    ncgs = 0
    nsgs = 0
       
    do ia = 1, natom
       !--- MPI : READ
       if (.TRUE.) then
          rewind (15)
       
          do 
             read (15, '(A2,I5)', iostat=eof) atname, nshell
             if (eof /= 0) then
                stop
             end if

             znum = atomic_number_get (atname)
             if (znum == atom_znum(ia)) then
                am = 0
                nprim = 0
                alpha = 0.0_dp
                coef  = 0.0_dp
                coef2 = 0.0_dp

                do is = 1, nshell
                   read (15, *) sym, nprim(is)
                   if (sym == 'SP') then
                      am(is) = -1
                   else if (sym == 'S ') then
                      am(is) = 0
                   else if (sym == 'P ') then
                      am(is) = 1
                   else if (sym == 'D ') then
                      am(is) = 2
                   else if (sym == 'F ') then
                      am(is) = 3
                   else if (sym == 'G ') then
                      am(is) = 4
                   end if
                   
                   if (am(is) .eq. -1) then
                      do ip = 1, nprim(is)
                         read (15, *) alpha(ip,is), coef(ip,is), coef2(ip,is)
                      end do
                   else
                      do ip = 1, nprim(is)
                         read (15, *) alpha(ip,is), coef(ip,is)
                      end do
                   end if

                end do

                exit
             end if
          end do
       end if

       call basis_read2 (ia, nshell, nshl, ncgs, nsgs, &
            nprim, am, alpha, coef, coef2)
    end do

    if (.TRUE.) then
       print *, 'NSHL ', qc_nshl
       print *, 'NGFS ', qc_ngfs
       print *, 'NCGS ', qc_ncgs
       close (15)
    end if
    
  end subroutine qc_basis_read


  subroutine basis_read1 (nshell, ncgs, nsgs, npgs)
    integer, intent(in) :: nshell
    integer, intent(out) :: ncgs, npgs, nsgs

    integer  :: is, ip, nprim
    real(dp) :: alpha
    character(len=2) :: sym

    ncgs = 0
    nsgs = 0
    npgs = 0

    do is = 1, nshell
       read (15, *) sym, nprim
         
       if (sym == 'S ') then
          ncgs = ncgs + 1
          nsgs = nsgs + 1
          npgs = npgs + nprim
          do ip = 1, nprim
             read (15, *) alpha
          end do
       else if (sym == 'SP') then
          ncgs = ncgs + 4
          nsgs = nsgs + 4
          npgs = npgs + 4*nprim
          do ip = 1, nprim
             read (15, *) alpha
          end do
       else if (sym == 'P ') then
          ncgs = ncgs + 3
          nsgs = nsgs + 3
          npgs = npgs + 3*nprim
          do ip = 1, nprim
             read (15, *) alpha
          end do
       else if (sym == 'D ') then
          ncgs = ncgs + 6
          nsgs = nsgs + 5
          npgs = npgs + 6*nprim
          do ip = 1, nprim
             read (15, *) alpha
          end do
       else if (sym == 'F ') then
          ncgs = ncgs + 10
          nsgs = nsgs + 7
          npgs = npgs + 10*nprim
          do ip = 1, nprim
             read (15, *) alpha
          end do
       else if (sym == 'G ') then
          ncgs = ncgs + 15
          nsgs = nsgs + 9
          npgs = npgs + 15*nprim
          do ip = 1, nprim
             read (15, *) alpha
          end do
       end if

    end do
    
  end subroutine basis_read1


  
  subroutine basis_read2 (ia, nshell, nshl, ncgs, nsgs, &
                          nprim, am, alpha, coef, coef2)
    integer, intent(in)    :: ia, nshell
    integer, intent(inout) :: ncgs, nsgs, nshl
    integer, intent(in)    :: nprim(20), am(20)
    real(dp),intent(in)    :: alpha(20,20),coef(20,20), coef2(20,20)
    integer  :: is, ip, jp
    real(dp) :: cnorm, aa, dum, fac, facs, pi32
    type(qc_shl_typ), pointer :: ishl

    do is = 1, nshell

       nshl = nshl + 1
       ishl => qc_shl_list(nshl)
       ishl%iat  = ia
       ishl%isgs = nsgs + 1
       ishl%icgs = ncgs + 1
       ishl%nprim= nprim(is)

       ishl%am   = am(is)  
       allocate (ishl%alpha(nprim(is) ))

       if (am(is) == -1) then
          !-- (S)--

          ishl%ncgs = 4
          ishl%nsgs = 4

          allocate (ishl%norm(2,nprim(is) ))

          nsgs = nsgs + 4
          ncgs = ncgs + 4

          do ip = 1, nprim(is)
             cnorm = exp(0.75_dp*log(2.0_dp*alpha(ip,is)/pi))
             ishl%alpha(ip)  = alpha(ip,is)
             ishl%norm(1,ip) = coef(ip,is)*cnorm
             cnorm = cnorm*sqrt(4.0_dp*alpha(ip,is))
             ishl%norm(2,ip)  = coef2(ip,is)*cnorm
          end do

       else if (am(is) == 0) then
          ishl%ncgs = 1
          ishl%nsgs = 1

          allocate (ishl%norm(1,nprim(is) ))

          do ip = 1, nprim(is)
             cnorm = exp(0.75_dp*log(2.0_dp*alpha(ip,is)/pi))
             ishl%alpha(ip) = alpha(ip,is)
             ishl%norm(1,ip)= coef(ip,is)*cnorm
          end do

          ncgs = ncgs + 1
          nsgs = nsgs + 1
          
          !--- normalize basis functions

          facs = 0.0_dp
          do ip = 1, nprim(is)
             do jp = 1, ip
                aa = ishl%alpha(ip) + ishl%alpha(jp)
                fac = aa * sqrt(aa)
                dum = ishl%norm(1,ip)*ishl%norm(1,jp)/fac
                if (ip .ne. jp) dum = dum + dum
                
                facs = facs + dum
             end do
          end do
          pi32 = 5.56832799683170_dp
          facs = 1.0_dp/sqrt (facs*pi32)

          do ip = 1, nprim(is)
             ishl%norm(1,ip) = ishl%norm(1,ip)*facs
          end do

       else if (am(is) == 1) then
          ishl%ncgs = 3
          ishl%nsgs = 3

          allocate (ishl%norm(1,nprim(is) ))

          do ip = 1, nprim(is)
             cnorm = exp(0.75_dp*log(2.0_dp*alpha(ip,is)/pi))
             cnorm = cnorm*sqrt(4.0_dp*alpha(ip,is))
             ishl%alpha(ip)   = alpha(ip,is)
             ishl%norm (1,ip) = coef(ip,is)*cnorm
          end do

          ncgs = ncgs + 3
          nsgs = nsgs + 3

          facs = 0.0_dp
          do ip = 1, nprim(is)
             do jp = 1, ip
                aa = ishl%alpha(ip) + ishl%alpha(jp)
                fac = aa * sqrt(aa)
                dum = 0.5_dp*ishl%norm(1,ip)*ishl%norm(1,jp)/(aa*fac)
                if (ip .ne. jp) dum = dum + dum
                
                facs = facs + dum
             end do
          end do
          pi32 = 5.56832799683170_dp
          facs = 1.0_dp/sqrt (facs*pi32)

          do ip = 1, nprim(is)
             ishl%norm(1,ip) = ishl%norm(1,ip)*facs
          end do

       else if (am(is) == 2) then
          ishl%ncgs = 6
          ishl%nsgs = 5
          allocate (ishl%norm(1,nprim(is) ))

          do ip = 1, nprim(is)
             cnorm = exp(0.75_dp*log(2.0_dp*alpha(ip,is)/pi))*sqrt(4.0_dp*alpha(ip,is))**2
             cnorm = cnorm/sqrt(3.0_dp)

             ishl%alpha(ip) = alpha(ip,is)
             ishl%norm (1,ip) = coef(ip,is)*cnorm  ! dxx
          end do

          ncgs = ncgs + 6
          nsgs = nsgs + 5

       else if (am(is) == 3) then
          ishl%ncgs = 10
          ishl%nsgs = 7

          allocate (ishl%norm(1,nprim(is) ))

          do ip = 1, nprim(is)
             cnorm = exp(0.75_dp*log(2.0_dp*alpha(ip,is) /pi))*sqrt(4.0_dp*alpha(ip,is))**3
             cnorm = cnorm/sqrt(15.0_dp)

             ishl%alpha(ip)   = alpha(ip,is)
             ishl%norm (1,ip) = coef(ip,is)*cnorm
          end do

          ncgs = ncgs + 10
          nsgs = nsgs + 7

       else if (am(is) == 4) then
          ishl%ncgs = 15
          ishl%nsgs = 9

          allocate (ishl%norm(1,nprim(is) ))

          do ip = 1, nprim(is)
             cnorm = exp(0.75_dp*log(2.0_dp*alpha(ip,is) /pi))*sqrt(4.0_dp*alpha(ip,is))**4
             cnorm = cnorm/sqrt(7.0_dp*15.0_dp)

             ishl%alpha(ip)   = alpha(ip,is)
             ishl%norm (1,ip) = coef(ip,is)*cnorm
          end do

          ncgs = ncgs + 15
          nsgs = nsgs + 9  
       end if

    end do
    
    
  end subroutine basis_read2



  subroutine qc_basis_mem_init

    allocate (qc_shl_list(qc_nshl))

  end subroutine qc_basis_mem_init


  subroutine qc_basis_mem_free
    integer :: is

    do is = 1, qc_nshl
       deallocate (qc_shl_list(is)%alpha)
       deallocate (qc_shl_list(is)%norm)
    end do

    deallocate (qc_shl_list)

  end subroutine qc_basis_mem_free


  subroutine qc_basis_print
    integer :: is
    
    do is = 1, qc_nshl
       print *, 'am  ', is, qc_shl_list(is)%am
       print *, 'iat ', is, qc_shl_list(is)%iat
       print *, 'isgs', is, qc_shl_list(is)%isgs
       print *, 'icgs', is, qc_shl_list(is)%icgs
       print *, 'nprim', is, qc_shl_list(is)%nprim
       print *, 'norm', is, qc_shl_list(is)%norm(1,1)
    end do

  end subroutine qc_basis_print

end module qc_basis
