module qc_mcmp2
  use qc_kind
  use qc_constant
  use nw_vectors

  implicit none

  type el_pair_typ
     real(dp), dimension(3) :: pos1, pos2
     real(dp), dimension(:), allocatable :: psi1, psi2
     real(dp)  :: wgt, r12
     logical   :: is_new
  end type el_pair_typ

  type ij_pair_typ
     real(dp), dimension(21) :: o_13, o_14, o_23, o_24
     real(dp), dimension(21) :: v_13, v_14, v_23, v_24
     real(dp) :: emp2, qeps(2), scs_emp2
  end type ij_pair_typ

  integer :: el_pair_num
  type (el_pair_typ), dimension(:),   allocatable, target :: el_pair_list
  type (ij_pair_typ), dimension(:,:), allocatable, target :: ij_pair_list

  integer :: iocc1, iocc2, ivir1, ivir2
  !--
  real(dp), dimension(:,:), allocatable :: te_val

contains


  subroutine mcmp2_mem_init
    integer :: ip

    allocate (el_pair_list(el_pair_num))
    allocate (ij_pair_list(el_pair_num, el_pair_num))

    do ip = 1, el_pair_num
       allocate (el_pair_list(ip)%psi1(iocc1:ivir2))
       allocate (el_pair_list(ip)%psi2(iocc1:ivir2))
    end do

  end subroutine mcmp2_mem_init


  subroutine mcmp2_mem_free
    integer :: ip

    do ip = 1, el_pair_num
       deallocate (el_pair_list(ip)%psi1)
       deallocate (el_pair_list(ip)%psi2)
    end do

    deallocate (el_pair_list)
    deallocate (ij_pair_list)

  end subroutine mcmp2_mem_free




  subroutine mc_te_val_init 
    real(dp), PARAMETER  :: xgk(21) &
         = (/-0.995657163025808080735527280689003_dp, &
         0.995657163025808080735527280689003_dp, &
         -0.973906528517171720077964012084452_dp, &
         0.973906528517171720077964012084452_dp, &
         -0.930157491355708226001207180059508_dp, &
         0.930157491355708226001207180059508_dp, &
         -0.865063366688984510732096688423493_dp, &
         0.865063366688984510732096688423493_dp, &
         -0.780817726586416897063717578345042_dp, &
         0.780817726586416897063717578345042_dp, &
         -0.679409568299024406234327365114874_dp, &
         0.679409568299024406234327365114874_dp, &
         -0.562757134668604683339000099272694_dp, &
         0.562757134668604683339000099272694_dp, &
         -0.433395394129247190799265943165784_dp, &
         0.433395394129247190799265943165784_dp, &
         -0.294392862701460198131126603103866_dp, &
         0.294392862701460198131126603103866_dp, &
         -0.148874338981631210884826001129720_dp, &
         0.148874338981631210884826001129720_dp, &
         0.000000000000000000000000000000000_dp/)

    integer  :: itau, im, am
    real(dp) :: center, hlgth, tau, xx, en_i, en_a

    iocc1 = nw_icore + 1
    iocc2 = nw_iocc
    ivir1 = nw_iocc + 1
    ivir2 = nw_nmo(1) 

    allocate (te_val(iocc1:ivir2, 21))

    center = 0.5_dp
    hlgth  = 0.5_dp

    do itau = 1, 21
       tau  = center + hlgth * xgk(itau)

       xx = (1.0_dp - tau)/tau

       do im = iocc1, iocc2
          en_i = nw_en(im)
          te_val(im,itau) = exp(en_i*xx)
       end do
       do am = ivir1, ivir2
          en_a = nw_en(am)
          te_val(am,itau) = exp(-en_a*xx)
       end do
    end do

  end subroutine mc_te_val_init


  subroutine mcmp2_local_energy (emp2, scs_emp2, qeps, lband)
    real(dp), intent(out):: emp2, scs_emp2, qeps(2)
    logical, intent(in)  :: lband
    real(dp), PARAMETER  :: xgk(21) &
    = (/-0.995657163025808080735527280689003_dp, &
         0.995657163025808080735527280689003_dp, &
        -0.973906528517171720077964012084452_dp, &
         0.973906528517171720077964012084452_dp, &
        -0.930157491355708226001207180059508_dp, &
         0.930157491355708226001207180059508_dp, &
        -0.865063366688984510732096688423493_dp, &
         0.865063366688984510732096688423493_dp, &
        -0.780817726586416897063717578345042_dp, &
         0.780817726586416897063717578345042_dp, &
        -0.679409568299024406234327365114874_dp, &
         0.679409568299024406234327365114874_dp, &
        -0.562757134668604683339000099272694_dp, &
         0.562757134668604683339000099272694_dp, &
        -0.433395394129247190799265943165784_dp, &
         0.433395394129247190799265943165784_dp, &
        -0.294392862701460198131126603103866_dp, &
         0.294392862701460198131126603103866_dp, &
        -0.148874338981631210884826001129720_dp, &
         0.148874338981631210884826001129720_dp, &
         0.000000000000000000000000000000000_dp/)
    real(dp), parameter :: wgk(21) &
     = (/0.011694638867371874278064396062192_dp, &
         0.011694638867371874278064396062192_dp, &
         0.032558162307964727478818972459390_dp, &
         0.032558162307964727478818972459390_dp, &
         0.054755896574351996031381300244580_dp, &
         0.054755896574351996031381300244580_dp, &
         0.075039674810919952767043140916190_dp, &
         0.075039674810919952767043140916190_dp, &
         0.093125454583697605535065465083366_dp, &
         0.093125454583697605535065465083366_dp, &
         0.109387158802297641899210590325805_dp, &
         0.109387158802297641899210590325805_dp, &
         0.123491976262065851077958109831074_dp, &
         0.123491976262065851077958109831074_dp, &
         0.134709217311473325928054001771707_dp, &
         0.134709217311473325928054001771707_dp, &
         0.142775938577060080797094273138717_dp, &
         0.142775938577060080797094273138717_dp, &
         0.147739104901338491374841515972068_dp, &
         0.147739104901338491374841515972068_dp, &
         0.149445554002916905664936468389821_dp/)

    integer  :: itau, ip, jp
    integer  :: icount2
    real(dp) :: center, hlgth, tau
    real(dp) :: a_resk, a_val, emp2a
    real(dp) :: b_resk, b_val, emp2b
    real(dp) :: q_resk(2), q_val(2)
    real(dp) :: r12, r34
    real(dp) :: wgt12, wgt34
    type (el_pair_typ), pointer :: el_pair1, el_pair2

    center = 0.5_dp
    hlgth  = 0.5_dp

    emp2   = 0.0_dp
    scs_emp2   = 0.0_dp
    qeps   = 0.0_dp

    icount2 = 0
    do ip = 1, el_pair_num - 1
       el_pair1 => el_pair_list(ip)

       do jp = ip + 1, el_pair_num
          el_pair2 => el_pair_list(jp)

          if (el_pair1%is_new .or. el_pair2%is_new) then
             a_resk = 0.0_dp
             b_resk = 0.0_dp
             q_resk = 0.0_dp

             do itau = 1, 21
                tau = center + hlgth*xgk(itau)
                call mc_emp2_func (ip, jp, tau, itau, a_val, b_val, q_val,lband)
                a_resk = a_resk + wgk(itau)*a_val
                b_resk = b_resk + wgk(itau)*b_val
                q_resk = q_resk + wgk(itau)*q_val
             end do

             r12   = el_pair1%r12
             r34   = el_pair2%r12
             wgt12 = el_pair1%wgt
             wgt34 = el_pair2%wgt

             emp2a = hlgth*a_resk/(r12*r34)/(wgt12*wgt34)
             emp2b = hlgth*b_resk/(r12*r34)/(wgt12*wgt34)

             ij_pair_list(ip,jp)%emp2 = -2.0_dp*emp2a + emp2b
             ij_pair_list(ip,jp)%scs_emp2 = -5.0_dp/6.0_dp*(2.0_dp*emp2a) + 1.0_dp/3.0_dp* emp2b
             ij_pair_list(ip,jp)%qeps = &
                  hlgth*q_resk/(r12*r34)/(wgt12*wgt34)
          end if
          icount2 = icount2 + 1
          emp2 = emp2 + ij_pair_list(ip,jp)%emp2
          scs_emp2 = scs_emp2 + ij_pair_list(ip,jp)%scs_emp2
          qeps = qeps + ij_pair_list(ip,jp)%qeps
       end do
    end do

    !icount2 = el_pair_num * (el_pair_num - 1)/2
    emp2 = emp2/dble(icount2)
    scs_emp2 = scs_emp2/dble(icount2)
    qeps = qeps/dble(icount2)

  end subroutine mcmp2_local_energy



  subroutine mc_emp2_func (ip, jp, tau, itau, aval, bval, qeps, lband)
    integer,  intent(in)  :: ip, jp
    real(dp), intent(in)  :: tau
    integer,  intent(in)  :: itau
    real(dp), intent(out) :: aval, bval, qeps(2)
    logical,  intent(in)  :: lband
    !---
    type(el_pair_typ), pointer :: el_pair1, el_pair2
    integer :: im, am, gm
    real(dp) :: o_13, o_14, o_23, o_24
    real(dp) :: v_13, v_14, v_23, v_24
    real(dp) :: p_13, p_23, p_14, p_24
    real(dp) :: xx
    real(dp) :: t_i, t_a, t_g

    xx    = (1.0_dp - tau)/tau

    el_pair1 => el_pair_list(ip)
    el_pair2 => el_pair_list(jp)

    o_13 = 0.0_dp !cmplx (0.0_dp, 0.0_dp, dp)
    o_14 = 0.0_dp !cmplx (0.0_dp, 0.0_dp, dp)
    o_23 = 0.0_dp !cmplx (0.0_dp, 0.0_dp, dp)
    o_24 = 0.0_dp !cmplx (0.0_dp, 0.0_dp, dp)
    v_13 = 0.0_dp !cmplx (0.0_dp, 0.0_dp, dp)
    v_14 = 0.0_dp !cmplx (0.0_dp, 0.0_dp, dp)
    v_23 = 0.0_dp !cmplx (0.0_dp, 0.0_dp, dp)
    v_24 = 0.0_dp !cmplx (0.0_dp, 0.0_dp, dp)

    qeps = 0.0_dp

    if (xx .le. 25.0d0) then 

       do im = iocc1, iocc2
          t_i  = te_val(im,itau) 
          o_13  = o_13 + (el_pair1%psi1(im))*el_pair2%psi1(im)*t_i 
          o_14  = o_14 + (el_pair1%psi1(im))*el_pair2%psi2(im)*t_i 
          o_23  = o_23 + (el_pair1%psi2(im))*el_pair2%psi1(im)*t_i 
          o_24  = o_24 + (el_pair1%psi2(im))*el_pair2%psi2(im)*t_i 
       end do

       do am = ivir1, ivir2
          t_a  = te_val(am,itau) 
          v_13 = v_13 + el_pair1%psi1(am)*(el_pair2%psi1(am))*t_a 
          v_14 = v_14 + el_pair1%psi1(am)*(el_pair2%psi2(am))*t_a 
          v_23 = v_23 + el_pair1%psi2(am)*(el_pair2%psi1(am))*t_a 
          v_24 = v_24 + el_pair1%psi2(am)*(el_pair2%psi2(am))*t_a 
       end do
    end if

    ij_pair_list(ip,jp)%o_13(itau) = o_13
    ij_pair_list(ip,jp)%o_14(itau) = o_14
    ij_pair_list(ip,jp)%o_23(itau) = o_23
    ij_pair_list(ip,jp)%o_24(itau) = o_24
    ij_pair_list(ip,jp)%v_13(itau) = v_13
    ij_pair_list(ip,jp)%v_14(itau) = v_14
    ij_pair_list(ip,jp)%v_23(itau) = v_23
    ij_pair_list(ip,jp)%v_24(itau) = v_24
    !---
    ij_pair_list(jp,ip)%o_13(itau) = o_13
    ij_pair_list(jp,ip)%o_14(itau) = o_23
    ij_pair_list(jp,ip)%o_23(itau) = o_14
    ij_pair_list(jp,ip)%o_24(itau) = o_24
    ij_pair_list(jp,ip)%v_13(itau) = v_13
    ij_pair_list(jp,ip)%v_14(itau) = v_23
    ij_pair_list(jp,ip)%v_23(itau) = v_14
    ij_pair_list(jp,ip)%v_24(itau) = v_24
    
    !---
    aval = (o_13*o_24*v_13*v_24)/(tau*tau)
    bval = (o_14*o_23*v_13*v_24)/(tau*tau)

    aval = aval + (o_14*o_23*v_14*v_23)/(tau*tau)
    bval = bval + (o_13*o_24*v_14*v_23)/(tau*tau)

    if (lband) then
       qeps = 0.0_dp
       do im = 1, 2  ! homo and lumo
          gm = (im - 1) + iocc2
          if (im .eq. 1) t_g  = te_val(gm,itau) 
          if (im .eq. 2) t_g  = 1.0_dp/te_val(gm,itau) 
          p_13 = (el_pair1%psi1(gm))*el_pair2%psi1(gm)*t_g
          p_23 = (el_pair1%psi2(gm))*el_pair2%psi1(gm)*t_g
          p_14 = (el_pair1%psi1(gm))*el_pair2%psi2(gm)*t_g
          p_24 = (el_pair1%psi2(gm))*el_pair2%psi2(gm)*t_g
          qeps(im) = qeps(im) - 2.0_dp*dble(p_13*o_24*v_13*v_24)/(tau*tau)
          qeps(im) = qeps(im) +        dble(o_14*p_23*v_13*v_24)/(tau*tau)
          qeps(im) = qeps(im) - 2.0_dp*dble(p_14*o_23*v_14*v_23)/(tau*tau)
          qeps(im) = qeps(im) +        dble(o_13*p_24*v_14*v_23)/(tau*tau)
          !---
          p_13 = el_pair1%psi1(gm)*(el_pair2%psi1(gm))/t_g
          p_14 = el_pair1%psi1(gm)*(el_pair2%psi2(gm))/t_g
          qeps(im) = qeps(im) + 2.0_dp*dble(o_13*o_24*p_13*v_24)/(tau*tau)
          qeps(im) = qeps(im) -        dble(o_14*o_23*p_13*v_24)/(tau*tau)
          qeps(im) = qeps(im) + 2.0_dp*dble(o_14*o_23*p_14*v_23)/(tau*tau)
          qeps(im) = qeps(im) -        dble(o_13*o_24*p_14*v_23)/(tau*tau)
       end do
    end if

    aval = 0.5_dp*aval
    bval = 0.5_dp*bval
    qeps = 0.5_dp*qeps

  end subroutine mc_emp2_func


end module qc_mcmp2
