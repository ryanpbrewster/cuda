module qc_psi
  use qc_kind
  use qc_basis
  use qc_constant
  use nw_vectors
  implicit none

  public :: qc_psi_mem_init
  public :: qc_psi_mem_free
  public :: qc_psi_get

  private
  real(dp), dimension(:), allocatable :: icgs
  integer :: iocc1, iocc2, ivir1, ivir2
  real(dp) :: ang(15), cd(2), cf (7), cg(11)

contains

  subroutine qc_psi_mem_init 
    
    iocc1 = nw_icore + 1
    iocc2 = nw_iocc
    ivir1 = nw_iocc+1
    ivir2 = nw_nmo(1)

    !print *, 'iocc1 : iocc2 ', iocc1, iocc2
    !print *, 'ivir1 : ivir2 ', ivir1, ivir2
    allocate (icgs(nw_nbf))
    
    cd(1) = sqrt(3.0_dp)
    cd(2) = cd(1)*0.5_dp
    !--
    cf(1) = sqrt(2.5_dp)*0.5_dp
    cf(2) = cf(1)*3.0_dp
    cf(3) = sqrt(15.0_dp)
    cf(4) = sqrt(1.5_dp)*0.5_dp
    cf(5) = sqrt(6.0_dp)
    cf(6) = 1.5_dp
    cf(7) = cf(3)*0.5_dp
    !--
    cg(1) = 2.9580398915498085_dp ! (3,1,0) (1,3,0)
    cg(2) = 6.2749501990055672_dp ! 
    cg(3) = 2.0916500663351894_dp ! 
    cg(4) = 1.1180339887498949_dp !
    cg(5) = 6.7082039324993694_dp
    cg(6) = 2.3717082451262845_dp !
    cg(7) = 3.1622776601683795_dp
    cg(8) = 0.55901699437494745_dp !
    cg(9) = 3.3541019662496847_dp !
    cg(10)= 0.73950997288745213_dp !
    cg(11)= 4.4370598373247132_dp !

  end subroutine qc_psi_mem_init


  subroutine qc_psi_mem_free

    deallocate (icgs)

  end subroutine qc_psi_mem_free



  subroutine qc_psi_get (pos, psi)
    real(dp), intent(in)  :: pos(3)
    real(dp), intent(out) :: psi(iocc1:ivir2)
    integer :: im, ic
    real(dp) :: psii

    !print *, 'pos ', pos
    call cgs_get (pos)

    !do ic = 1, nw_nbf
    !   print *, ic, icgs(ic)
    !end do
    
    mo_loop : do im = iocc1, ivir2
       psii = 0.0_dp !cmplx(0.0_dp, 0.0_dp, dp)

       do ic = 1, nw_nbf
          psii = psii + icgs(ic)*nw_co(ic,im)
          if (abs(nw_co(ic,im)) .gt. 1.0e-3_dp ) then
             !print *, 'co ic im', ic, im, nw_co(ic,im)
          end if
       end do

       psi(im) = psii
       !print *, 'psi im', im, psi(im)
    end do mo_loop

!    stop

  end subroutine qc_psi_get




  subroutine cgs_get (pos)
    real(dp), dimension(3), intent(in) :: pos
    integer :: is, iat, iam, ic, ip
    real(dp):: r2, zti, rad(2)
    real(dp):: pos_i(3), dr(3), x, y, z
    type (qc_shl_typ), pointer :: ishl

    do is = 1, qc_nshl
       ishl => qc_shl_list(is)
       iat  =  ishl%iat
       pos_i = atom_pos(1:3,iat)  

       iam  =  ishl%am
       
       dr = pos - pos_i
       r2 = dot_product (dr, dr)
       
       x = dr(1)
       y = dr(2)
       z = dr(3)

       !---
       rad = 0.0_dp

       do ip = 1, ishl%nprim
          zti = ishl%alpha(ip)
          rad(1) = rad(1) + exp(-zti*r2)*ishl%norm(1,ip) ! S
       end do

       if (ishl%am .eq. -1) then  !SP
          do ip = 1, ishl%nprim
             zti = ishl%alpha(ip)
             rad(2) = rad(2) + exp(-zti*r2)*ishl%norm(2,ip) ! Px
          end do
       end if
       !---
       ang = 0.0_dp

       if (lspherical) then
          
          ic = ishl%isgs - 1
          if (ishl%am .eq. 0) then
             
             icgs(ic+1) = rad(1)

          else if (ishl%am .eq. -1) then
             
             icgs(ic+1) = rad(1) 
             icgs(ic+2) = rad(2)*x
             icgs(ic+3) = rad(2)*y
             icgs(ic+4) = rad(2)*z
          else if (ishl%am .eq. 1) then
             
             icgs(ic+1) = rad(1)*x
             icgs(ic+2) = rad(1)*y
             icgs(ic+3) = rad(1)*z

          else if (ishl%am .eq. 2) then
             ! l = 2, m = -2
             ang(1) = cd(1)*x*y ! dxy
             ! l = 2, m = -1
             ang(2) = cd(1)*y*z ! dyz
             ! l = 2, m = 0
             ang(3) = 0.5_dp*(2.0_dp*z*z - x*x - y*y) ! dxx, dyy, dzz
             ang(4) = -cd(1)*x*z ! dxz
             ang(5) =  cd(2)*(x*x - y*y) ! dxx, dyy
             !--(test) --
             !ang(1) = x*y ! dxy
             !ang(2) = y*z ! dyz
             !ang(3) = 0.5_dp*(2.0_dp*z*z - x*x - y*y) ! dxx, dyy, dzz
             !ang(4) = -x*z ! dxz
             !ang(5) =  cd(2)*(x*x - y*y) ! dxx, dyy
             !--
             icgs(ic+1) = rad(1)*ang(1)
             icgs(ic+2) = rad(1)*ang(2)
             icgs(ic+3) = rad(1)*ang(3)
             icgs(ic+4) = rad(1)*ang(4)
             icgs(ic+5) = rad(1)*ang(5)

          else if (ishl%am .eq. 3) then

             ang(1) = y*(cf(2)*x*x-cf(1)*y*y) ! xxy, yyy
             ang(2) = cf(3)*x*y*z
             ang(3) = y*(cf(5)*z*z-cf(4)*(x*x+y*y))
             ang(4) = z*(z*z-cf(6)*(x*x+y*y))
             ang(5) =-x*(cf(5)*z*z-cf(4)*(x*x+y*y))
             ang(6) = z*cf(7)*(x*x-y*y)
             ang(7) = x*(cf(2)*y*y-cf(1)*x*x)

             icgs(ic+1) = rad(1)*ang(1)
             icgs(ic+2) = rad(1)*ang(2)
             icgs(ic+3) = rad(1)*ang(3)
             icgs(ic+4) = rad(1)*ang(4)
             icgs(ic+5) = rad(1)*ang(5)
             icgs(ic+6) = rad(1)*ang(6)
             icgs(ic+7) = rad(1)*ang(7)

          else if (ishl%am .eq. 4) then
             ! m = -4
             ang(1) = cg(1)*(x*x*x*y-x*y*y*y) ! xxxy, xyyy
             ! m = -3
             ang(2) = y*z*(cg(2)*x*x - cg(3)*y*y) ! (2,1,1) (0,3,1)
             ! m = -2
             ang(3) = x*y*cg(4)*(-x*x - y*y) + cg(5)*x*y*z*z
             ! m = -1
             ang(4) = -cg(6)*x*x*y*z - cg(6)*y*y*y*z + cg(7)*y*z*z*z
             ! m = 0
             ang(5) = 0.375_dp*(x*x*x*x + y*y*y*y + 2.0_dp*x*x*y*y) + &
                  z*z*z*z - 3.0_dp*z*z*(x*x + y*y) 
             ! m = 1
             ang(6) = cg(6)*x*x*x*z + cg(6)*x*y*y*z - cg(7)*x*z*z*z 
             ! m = 2
             ang(7) = cg(8)*(y*y*y*y - x*x*x*x) + &
                      cg(9)*z*z*(x*x - y*y)
             ! m = 3
             ang(8) = x*z*(cg(2)*y*y - cg(3)*x*x) ! (1,2,1) (3,0,1)
             ! m = 4
             ang(9) = cg(10)*(x*x*x*x + y*y*y*y) - cg(11)*x*x*y*y

             icgs(ic+1) = rad(1)*ang(1)
             icgs(ic+2) = rad(1)*ang(2)
             icgs(ic+3) = rad(1)*ang(3)
             icgs(ic+4) = rad(1)*ang(4)
             icgs(ic+5) = rad(1)*ang(5)
             icgs(ic+6) = rad(1)*ang(6)
             icgs(ic+7) = rad(1)*ang(7)
             icgs(ic+8) = rad(1)*ang(8)
             icgs(ic+9) = rad(1)*ang(9)
          end if

       else
          ! cartesian GTO
          !print *, 'ish am ic ', is, ishl%am, ishl%icgs

          ic = ishl%icgs - 1
          if (ishl%am .eq. 0) then
             
             icgs(ic+1) = rad(1)

          else if (ishl%am .eq. -1) then
             
             icgs(ic+1) = rad(1) 
             icgs(ic+2) = rad(2)*x
             icgs(ic+3) = rad(2)*y
             icgs(ic+4) = rad(2)*z
          else if (ishl%am .eq. 1) then
             
             icgs(ic+1) = rad(1)*x
             icgs(ic+2) = rad(1)*y
             icgs(ic+3) = rad(1)*z

          else if (ishl%am .eq. 2) then

             ang(1) = x*x ! dxx
             ang(2) = x*y !cd(1)*x*y ! dxy
             ang(3) = x*z !cd(1)*x*z ! dxz
             ang(4) = y*y ! dyy
             ang(5) = y*z ! cd(1)*y*z ! dyz
             ang(6) = z*z ! dzz
             !--
             icgs(ic+1) = rad(1)*ang(1)
             icgs(ic+2) = rad(1)*ang(2)
             icgs(ic+3) = rad(1)*ang(3)
             icgs(ic+4) = rad(1)*ang(4)
             icgs(ic+5) = rad(1)*ang(5)
             icgs(ic+6) = rad(1)*ang(6)

          else if (ishl%am .eq. 3) then

             ang(1) = x*x*x       ! fxxx
             ang(2) = x*x*y       !sqrt_5*x*x*y ! fxxy
             ang(3) = x*x*z       !sqrt_5*x*x*z ! fxxz
             ang(4) = x*y*y       !sqrt_5*x*y*y ! fxyy
             ang(5) = x*y*z       !cf(3)*x*y*z ! fxyz
             ang(6) = x*z*z       !sqrt_5*x*z*z ! fxzz
             ang(7) = y*y*y ! fyyy
             ang(8) = y*y*z !sqrt_5*y*y*z ! fyyz
             ang(9) = y*z*z !sqrt_5*y*z*z ! fyzz
             ang(10) = z*z*z ! fzzz

             icgs(ic+1) = rad(1)*ang(1)
             icgs(ic+2) = rad(1)*ang(2)
             icgs(ic+3) = rad(1)*ang(3)
             icgs(ic+4) = rad(1)*ang(4)
             icgs(ic+5) = rad(1)*ang(5)
             icgs(ic+6) = rad(1)*ang(6)
             icgs(ic+7) = rad(1)*ang(7)
             icgs(ic+8) = rad(1)*ang(8)
             icgs(ic+9) = rad(1)*ang(9)
             icgs(ic+10) = rad(1)*ang(10)
             
          else if (ishl%am .eq. 4) then

             ang(1) = x*x*x*x    ! (4,0,0)
             ang(2) = x*x*x*y    ! (3,1,0)  
             ang(3) = x*x*x*z    ! (3,0,1)   
             ang(4) = x*x*y*y    ! (2,2,0)   
             ang(5) = x*x*y*z    ! (2,1,1)   
             ang(6) = x*x*z*z    ! (2,0,2)   
             ang(7) = x*y*y*y    ! (1,3,0)
             ang(8) = x*y*y*z    ! (1,2,1)
             ang(9) = x*y*z*z    ! (1,1,2)
             ang(10)= x*z*z*z    ! (1,0,3)
             !
             ang(11)= y*y*y*y    ! (0,4,0)  
             ang(12)= y*y*y*z    ! (0,3,1)
             ang(13)= y*y*z*z    ! (0,2,2)
             ang(14)= y*z*z*z    ! (0,1,3)
             ang(15)= z*z*z*z    ! (0,0,4)

             icgs(ic+1) = rad(1)*ang(1)
             icgs(ic+2) = rad(1)*ang(2)
             icgs(ic+3) = rad(1)*ang(3)
             icgs(ic+4) = rad(1)*ang(4)
             icgs(ic+5) = rad(1)*ang(5)
             icgs(ic+6) = rad(1)*ang(6)
             icgs(ic+7) = rad(1)*ang(7)
             icgs(ic+8) = rad(1)*ang(8)
             icgs(ic+9) = rad(1)*ang(9)
             icgs(ic+10) = rad(1)*ang(10)
             icgs(ic+11) = rad(1)*ang(11)
             icgs(ic+12) = rad(1)*ang(12)
             icgs(ic+13) = rad(1)*ang(13)
             icgs(ic+14) = rad(1)*ang(14)
             icgs(ic+15) = rad(1)*ang(15)
          end if

       end if
    end do

  end subroutine cgs_get


end module qc_psi
