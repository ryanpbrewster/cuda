module qc_constant
  use qc_kind
  implicit none

  real(dp), parameter :: pi = 3.14159265358979_dp
  real(dp), parameter :: rttwo = 1.41421356237309504880_dp
  real(dp), parameter :: bohr_to_ang = 0.529177249_dp
  real(dp), parameter :: ev = 27.2113957_dp
  real(dp), parameter :: ang_to_bohr = 1.0_dp/bohr_to_ang !0.52917720859_dp
  real(dp) :: sqrt_pi, sqrt_pi3, sqrt_5
  real(dp) :: pisub
  character(len=2), parameter :: atom_name(56) = &
        (/'H ','He', & 
          'Li','Be','B ','C ','N ','O ','F ','Ne', & 
          'Na','Mg','Al','Si','P ','S ','Cl','Ar', & 
          'K ','Ca', & 
          'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn', & 
          'Ga','Ge','As','Se','Br','Kr', & 
          'Rb','Sr', & 
          'Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd', & 
          'In','Sn','Sb','Te','I ','Xe', & 
          'Cs','Ba'/) !, &
!          'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg'/)
          
  ! fmt
  real(dp), parameter :: tf(0:16) = &
       (/33.0_dp,37.0_dp,41.0_dp,43.0_dp,46.0_dp,49.0_dp,51.0_dp,54.0_dp,56.0_dp, &
       58.0_dp,61.0_dp,63.0_dp,66.0_dp,68.0_dp,70.0_dp,72.0_dp,74.0_dp/)
  real(dp) :: igamma(0:1500,0:16)

  !--- (for int1 and eri)---
  real(dp) :: f1(0:9), f2(0:9)

contains

  subroutine qc_constant_init
    integer :: i, j, k
    real(dp)  :: a, b, c

    sqrt_pi = sqrt(pi)
    sqrt_pi3 = sqrt_pi**3
    sqrt_5  = sqrt(5.0_dp)

    pisub = 2.0d0/sqrt_pi
    f1(0) = 1.0d0/pisub
    f1(1) = 1.0d0/2.0d0/pisub
    f1(2) = 3.0d0/4.0d0/pisub
    f1(3) = 15.0d0/8.0d0/pisub
    f1(4) = 105.0_dp/16.0_dp/pisub
    f1(5) = 945.0_dp/32.0_dp/pisub
    f1(6) = 10395.0_dp/64.0_dp/pisub
    f1(7) = 135135.0_dp/128.0_dp/pisub
    f1(8) = 2027025.0_dp/256.0_dp/pisub
    f1(9) = 34459425.0_dp/512.0_dp/pisub
    
    f2(0) = -0.5_dp
    f2(1) = -1.5_dp
    f2(2) = -2.5_dp
    f2(3) = -3.5_dp
    f2(4) = -4.5_dp
    f2(5) = -5.5_dp
    f2(6) = -6.5_dp
    f2(7) = -7.5_dp
    f2(8) = -8.5_dp
    f2(9) = -9.5_dp

    ! preliminary for fmt evalulation
    igamma = 0.0_dp

    do i = 0, 1500
       a = dble(i)/20.0_dp
       b = exp(-a)
       c = 0.0_dp
       k = 50 + int(i/3.0_dp)
       do j = k, 0, -1
          c = (2.0_dp*a*c + b)/dble(2*j + 1)
          if (j <= 16) igamma(i,j) = c
       end do
    end do
   
  end subroutine qc_constant_init


  function atomic_number_get (atnm) result (znum)
    character(len=2), intent(in) :: atnm
    integer :: i, nsize, znum
    
    znum = 0

    nsize = size (atom_name)
    do i = 1, nsize
       if (atom_name(i) == atnm) then
          znum = i
          exit
       end if
    end do

  end function atomic_number_get

end module qc_constant
