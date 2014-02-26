module qc_step
  use qc_kind
  use qc_constant
  implicit none
  
  private
  real(dp):: cpu_start, step_cpu_start
  
  public  :: qc_start, qc_end, qc_abort
  public  :: qc_step_start, qc_step_end
  public  :: cputime

contains

  subroutine qc_start
    integer :: ierr, hostnm, pid, getpid
    character(len=30) :: host_name, user_name
    character(len=28) :: datx

    call qc_constant_init

    cpu_start = cputime ()

    if (.TRUE.) then
       ierr = hostnm (host_name)
       call getlog   (user_name)
       pid = getpid ()
       
       call datum (datx)

       write (*,'("|",68(1H=),"|")')
       write (*,"(A,T35,A30)") "  PROGRAM STARTED AT ", ADJUSTR(datx)
       write (*,"(A,T31,A30)") "  PROGRAM STARTED ON ", ADJUSTR(host_name)
       write (*,"(A,T31,A30)") "  PROGRAM STARTED BY ", ADJUSTR(user_name)
       write (*,"(A,T51,I10)") "  PROGRAM PROCESS ID ", pid
    end if

  end subroutine qc_start


  subroutine qc_end
    real(dp) :: cpu_end
    real(dp) :: dcpu
    character (len=70) :: buffer1, buffer2
    character (len=50) :: cpu_info
    integer            :: len, leni
    buffer1= &
    "|==>                                                              <==|"
    buffer2= &
    "|==>                                                              <==|"

    if (.TRUE.) then
       cpu_end = cputime ()

       dcpu = cpu_end - cpu_start
       if (dcpu .lt. 86400) then
          write (cpu_info,'("TOTAL CPU TIME  : ", F12.3, " sec. (", F8.2, " hr )")') &
               dcpu, (dcpu)/3600.0d0
          len = len_trim (cpu_info)
          leni = 35 - len/2

          buffer1(leni:leni+len) = trim(cpu_info)

          write (cpu_info,'("TOTAL CPU TIME  : ", F12.3, " min. (", F8.2, " hr )")') &
               dcpu/60.0d0, (dcpu)/3600.0d0
          len = len_trim (cpu_info)
          leni = 35 - len/2

          buffer2(leni:leni+len) = trim(cpu_info)
       else
          write (cpu_info,'("TOTAL CPU TIME  : ", F12.3, " min. (", F8.2, " hr )")') &
               dcpu/60.0d0, (dcpu)/3600.0d0
          len = len_trim (cpu_info)
          leni = 35 - len/2

          buffer1(leni:leni+len) = trim(cpu_info)

          write (cpu_info,'("TOTAL CPU TIME  : ", F12.3, " hr.  (", F8.2, " day)")') &
               dcpu/3600.0d0, (dcpu)/86400.0d0
          len = len_trim (cpu_info)
          leni = 35 - len/2

          buffer2(leni:leni+len) = trim(cpu_info)
       end if

       write (*,*)
       write (*,'("|",68(1H=),"|")')
       write (*,'(A70)') buffer1
       write (*,'(A70)') buffer2
       write (*,'("|",68(1H=),"|")')

       call flush (6)
    end if

  end subroutine qc_end


  subroutine qc_abort(comment)
    real(dp) :: cpu_end
    real(dp) :: dcpu
    character (len=*)  :: comment
    character (len=70) :: buffer1, buffer2, buffer
    character (len=50) :: cpu_info
    integer            :: len, leni
    buffer1= &
    "|==>                                                              <==|"
    buffer2= &
    "|==>                                                              <==|"
    buffer = &
    "|@--                                                              ---|"

    if (.TRUE.) then
       len = len_trim (comment)

       leni = 35 - len/2
       buffer (leni:leni+len) = trim(comment)

       cpu_end = cputime ()

       dcpu = cpu_end - cpu_start
       if (dcpu .lt. 86400) then
          write (cpu_info,'("TOTAL CPU TIME  : ", F12.3, " sec. (", F8.2, " hr )")') &
               dcpu, (dcpu)/3600.0d0
          len = len_trim (cpu_info)
          leni = 35 - len/2

          buffer1(leni:leni+len) = trim(cpu_info)

          write (cpu_info,'("TOTAL CPU TIME  : ", F12.3, " min. (", F8.2, " hr )")') &
               dcpu/60.0d0, (dcpu)/3600.0d0
          len = len_trim (cpu_info)
          leni = 35 - len/2

          buffer2(leni:leni+len) = trim(cpu_info)
       else
          write (cpu_info,'("TOTAL CPU TIME  : ", F12.3, " min. (", F8.2, " hr )")') &
               dcpu/60.0d0, (dcpu)/3600.0d0
          len = len_trim (cpu_info)
          leni = 35 - len/2

          buffer1(leni:leni+len) = trim(cpu_info)

          write (cpu_info,'("TOTAL CPU TIME  : ", F12.3, " hr.  (", F8.2, " day)")') &
               dcpu/3600.0d0, (dcpu)/86400.0d0
          len = len_trim (cpu_info)
          leni = 35 - len/2

          buffer2(leni:leni+len) = trim(cpu_info)
       end if

       write (*,*)
       write (*,'("|",68(1H=),"|")')
       write (*,'(A70)') buffer
       write (*,'(A70)') buffer1
       write (*,'(A70)') buffer2
       write (*,'("|",68(1H=),"|")')

       call flush (6)
    end if

    stop

  end subroutine qc_abort


  subroutine qc_step_start (jobs)
    character (len=*)  :: jobs
    character (len=70) :: buffer
    integer            :: len, leni
    buffer = "|@--                                                              ---|"

    len = len_trim (jobs)

    leni = 35 - len/2
    buffer (leni:leni+len) = trim(jobs)

    if (.TRUE.) then
       write (*,*)
       write (*,'("|",68(1H-),"|")')
       write (*,'(A70)') buffer
       write (*,'("|",68(1H-),"|")')
       write (*,*)


       call flush(6)
    end if

    step_cpu_start = cputime ()

    print *, 'End of qc_step_start'

  end subroutine qc_step_start



  subroutine qc_step_end
    real(dp) :: cpu_end
    real(dp) :: dcpu
    character (len=70) :: buffer
    character (len=50) :: cpu_info
    integer            :: len, leni
    buffer = &
    "|>>>                                                              <<<|"

    !return
    cpu_end = cputime ()

    dcpu = cpu_end - step_cpu_start

    if (.TRUE.) then
       if (dcpu .lt. 3600.0d0) then
          write (cpu_info,'(" STEP CPU TIME  : ", F12.3, " sec. (", F8.2, " min)")') &
               dcpu, (dcpu)/60.0d0
       else
          write (cpu_info,'(" STEP CPU TIME  : ", F12.3, " min. (", F8.2, " hr )")') &
               dcpu/60.0d0, (dcpu)/3600.0d0
       end if
       
       len = len_trim (cpu_info)
       
       leni = 35 - len/2
       
       buffer (leni:leni+len) = trim(cpu_info)
       
       write (*,'("|",68(1H-),"|")')
       write (*,'(A70)') buffer
       write (*,'("|",68(1H-),"|")')
       write (*,*)

       call flush (6)
    end if

  end subroutine qc_step_end




  function cputime () result (second)
    character(len=10)     :: date,time,zone
    integer, dimension(8) :: value
    real(dp)              :: second

    call date_and_time (date, time, zone, value)
    ! value(3) : day
    ! vaule(5) : hr
    ! value(6) : min
    ! vaule(7) : sec
    ! value(8) : millisec.
    second = value(8)*0.001d0 + value(7) + value(6)*60.0d0 + &
         (value(5) + 24.0d0 * value(3))* 3600.0d0

  end function cputime


  subroutine datum (cal_date)
    character(len=*), intent(out) :: cal_date
    character(len=10) :: time
    character(len=5 ) :: zone
    character(len=8 ) :: date
    integer, dimension(8) :: values
    
    call date_and_time (date=date, time=time, zone=zone, values=values)

    cal_date=date//" "//time

  end subroutine datum

end module qc_step
