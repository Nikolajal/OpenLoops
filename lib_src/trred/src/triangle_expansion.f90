!******************************************************************************!
!                                                                              !
!    triangle_expansion.f90                                                    !
!    is part of trred & OpenLoops2                                             !
!    Copyright (C) 2017-2018 Federico Buccioni, Jean-Nicolas Lang,             !
!                            Stefano Pozzorini, Hantian Zhang and Max Zoller   !
!                                                                              !
!    trred has been developed by J.-N. Lang, H. Zhang and F. Buccioni          !
!    trred is licenced under the GNU GPL version 3,                            !
!    see COPYING for details.                                                  !
!                                                                              !
!******************************************************************************!


module triangle_expansion_DP
  use c0_0mm_DP, only: C0_n_0mm
  use c0_m00_DP, only: C0_n_m00 
  use c0_mmm_DP, only: C0_n_mmm, C0_n_mmm_init, C0_n_mmm_update
  use c0_000_DP, only: C0_n_000, C0_n_000_EP1
  use b0_mm_DP, only: B0_n_mm, B0_n_mm_init, B0_n_mm_update
  use b0_DP, only: B0_0_n, B0_n
  use b0_m0m1_DP, only: B0_n_m0m1
  use triangle_aux_DP, only: dp, target_precision, cone, cnul, rnul
  
  implicit none

  abstract interface
    function coeff_func(p2,m2,muUV2,muIR2,n)
      use triangle_aux_DP, only: dp
      implicit none
      complex(dp), intent(in) :: p2,m2(:),muUV2,muIR2
      integer,       intent(in) :: n
      complex(dp)             :: coeff_func
    end function coeff_func
  end interface

  contains

  function coeff_template(p2,m2,muUV2,muIR2,d,offset,coeff_init,coeff_update) result(coeff)
    ! Computes the coefficient function func(p2,d,m2,muUV2,muIR2,offset) := 
    ! \sum_{n=offset}^\infty d^{n-offset} func^(n)
    complex(dp),                  intent(in) :: m2(:)
    real(dp),                     intent(in) :: p2,muUV2,muIR2,d
    integer, optional,              intent(in) :: offset
    procedure(coeff_func), pointer, intent(in) :: coeff_init,coeff_update
    complex(dp) :: coeff,coeffn
    integer    :: n,m

    if (d .lt. 0 .or. d .gt. 1) then
      write (*,*) 'ERROR: called coeff_template with invalid d=', d
      stop
    end if

    if (present(offset)) then
      if (offset .lt. 0) then
        write (*,*) 'ERROR: called coeff_template with offset<0'
        stop
      end if
      m = offset
    else
      m = 0
    end if

    n = m
    coeff = coeff_init(cmplx(p2,kind=dp),m2,cmplx(muUV2,kind=dp),cmplx(muIR2,kind=dp),n)

    n = n + 1
    coeffn = d**(n-m)*coeff_update(cmplx(p2,kind=dp),m2,cmplx(muUV2,kind=dp),cmplx(muIR2,kind=dp),n)

    do while (abs(coeffn/coeff) .gt. target_precision)
      coeff = coeff + coeffn
      n = n + 1
      coeffn = d**(n-m)*coeff_update(cmplx(p2,kind=dp),m2,cmplx(muUV2,kind=dp),cmplx(muIR2,kind=dp),n)
    end do
    coeff = coeff + coeffn
    
  end function coeff_template

  function list_template(p2,m2,muUV2,muIR2,d,clist,list_init,list_update) result(list)
    ! Computes \sum_i clist(i)*func(p2,d,m2,muUV2,muIR2,i) in an efficient way,
    ! computing equal coefficients only once.
    complex(dp),                  intent(in) :: m2(:)
    real(dp),                     intent(in) :: p2,muUV2,muIR2,d
    complex(dp),                  intent(in) :: clist(0:)
    procedure(coeff_func), pointer, intent(in) :: list_init,list_update
    complex(dp) :: list,listn
    integer       :: n,m,k

    if (d .lt. 0 .or. d .gt. 1) then
      write (*,*) 'ERROR: called list_template with invalid d=',d
      stop
    end if

    do m = 1, size(clist)
      if (clist(m-1) .ne. 0._dp) then
        exit
      end if
    end do
    if (m .eq. size(clist) .and. clist(m-1) .eq. 0._dp  ) then
      write (*,*) 'ERROR: called list_template with zero clist.'
      stop 
    end if
    m = m - 1

    n = m
    list = clist(m)*list_init(cmplx(p2,kind=dp),m2,cmplx(muUV2,kind=dp),cmplx(muIR2,kind=dp),n)

    n = n + 1
    listn = 0._dp
    do k = m, min(n,size(clist)-1)
      if (clist(k) .ne. 0._dp) then
        listn = listn + d**(n-k)*clist(k)
      end if
    end do
    listn = listn*list_update(cmplx(p2,kind=dp),m2,cmplx(muUV2,kind=dp),cmplx(muIR2,kind=dp),n)

    ! the condition n .lt. size(clist)-1 ensures that for d=0 all coefficients
    ! are taken into account. For instance for a list with a zero
    ! clist=[1.,1.,0.,1.] this would otherwise lead to a mistake.
    do while (abs(listn/list) .gt. target_precision .or. n .lt. size(clist)-1)
      list = list + listn
      n = n + 1
      listn = cnul
      do k = m, min(n,size(clist)-1)
        if (clist(k) .ne. 0._dp) then
          listn = listn + d**(n-k)*clist(k)
        end if
      end do
      listn = listn*list_update(cmplx(p2,kind=dp),m2,cmplx(muUV2,kind=dp),cmplx(muIR2,kind=dp),n)
    end do
    list = list + listn
    
  end function list_template


  !!!!!!!!
  !  B0  !
  !!!!!!!!

  function B0d_coeff(p2,m2,muUV2,d,offset) result(B0)
    complex(dp), intent(in)           :: m2
    real(dp),    intent(in)           :: p2,d,muUV2
    integer,       intent(in), optional :: offset
    complex(dp) :: B0
    procedure(coeff_func), pointer :: coeff_init

    coeff_init => B0_n
    B0 = coeff_template(p2,[m2],muUV2,rnul,d,offset,coeff_init,coeff_init)

  end function B0d_coeff

  function B0d_list(p2,m2,muUV2,d,clist) result(B0)
    complex(dp), intent(in) :: m2
    real(dp),    intent(in) :: p2,d,muUV2
    complex(dp), intent(in) :: clist(0:)
    complex(dp) :: B0
    procedure(coeff_func), pointer :: coeff_init

    coeff_init => B0_n
    B0 = list_template(p2,[m2],muUV2,rnul,d,clist,coeff_init,coeff_init)

  end function B0d_list

  function B0d_mm_opt(p2,m2,muUV2,d,offset) result(B0)
    complex(dp), intent(in)        :: m2
    real(dp),    intent(in)        :: p2,d,muUV2
    integer,    intent(in), optional :: offset
    complex(dp)                    :: B0
    procedure(coeff_func), pointer   :: coeff_init,coeff_update

    coeff_init => B0_n_mm_init
    coeff_update => B0_n_mm_update
    B0 = coeff_template(p2,[m2],muUV2,rnul,d,offset,coeff_init,coeff_init)

  end function B0d_mm_opt

  function B0d_mm_opt_list(p2,m2,muUV2,d,clist) result(B0)
    complex(dp), intent(in)      :: m2
    real(dp),    intent(in)      :: p2,d,muUV2
    complex(dp), intent(in)      :: clist(0:)
    complex(dp)                  :: B0
    procedure(coeff_func), pointer :: coeff_init,coeff_update

    coeff_init => B0_n_mm_init
    coeff_update => B0_n_mm_update
    B0 = list_template(p2,[m2],muUV2,rnul,d,clist,coeff_init,coeff_update)

  end function B0d_mm_opt_list

  function B0d_0_coeff(p2,muUV2,d,offset) result(B0)
    real(dp),    intent(in)           :: p2,d,muUV2
    integer,       intent(in), optional :: offset
    complex(dp)                       :: B0
    procedure(coeff_func),      pointer :: coeff_init

    coeff_init => B0_0_n
    B0 = coeff_template(p2,[cnul],muUV2,rnul,d,offset,coeff_init,coeff_init)
  end function B0d_0_coeff

  function B0d_0_list(p2,muUV2,d,clist) result(B0)
    real(dp),    intent(in)      :: p2,d,muUV2
    complex(dp), intent(in)      :: clist(0:)
    complex(dp)                  :: B0
    procedure(coeff_func), pointer :: coeff_init,coeff_update

    coeff_init => B0_0_n
    B0 = list_template(p2,[cnul],muUV2,rnul,d,clist,coeff_init,coeff_init)

  end function B0d_0_list

  function B0d_m0m1_coeff(p2,m02,m12,muUV2,d,offset) result(B0)
    use b0_m0m1_DP, only: B0_n_m0m1
    real(dp),    intent(in)           :: p2,d,muUV2
    complex(dp), intent(in)           :: m02,m12
    integer,       intent(in), optional :: offset
    complex(dp)                       :: B0
    procedure(coeff_func), pointer      :: coeff_init

    coeff_init => B0_n_m0m1
    B0 = coeff_template(p2,[m02,m12],muUV2,rnul,d,offset,coeff_init,coeff_init)

  end function B0d_m0m1_coeff

  function B0d_m0m1_list(p2,m02,m12,muUV2,d,clist) result(B0)
    use b0_m0m1_DP, only: B0_n_m0m1
    real(dp),    intent(in)      :: p2,d,muUV2
    complex(dp), intent(in)      :: m02,m12
    complex(dp), intent(in)      :: clist(0:)
    complex(dp)                  :: B0
    procedure(coeff_func), pointer :: coeff_init

    coeff_init => B0_n_m0m1
    B0 = list_template(p2,[m02,m12],muUV2,rnul,d,clist,coeff_init,coeff_init)

  end function b0d_m0m1_list

  !!!!!!!!
  !  C0  !
  !!!!!!!!

  function C0d_0mm_coeff(p2,m2,d,offset) result(C0)
    complex(dp), intent(in)           :: m2
    real(dp),    intent(in)           :: p2,d
    integer,       intent(in), optional :: offset
    complex(dp)                       :: C0
    procedure(coeff_func),     pointer  :: coeff_init

    coeff_init => C0_n_0mm
    C0 = coeff_template(p2,[m2],rnul,rnul,d,offset,coeff_init,coeff_init)

  end function C0d_0mm_coeff

  function C0d_m00_coeff(p2,m2,muIR2,d,offset) result(C0)
    complex(dp), intent(in)           :: m2
    real(dp),    intent(in)           :: p2,d,muIR2
    integer,       intent(in), optional :: offset
    complex(dp)                       :: C0
    procedure(coeff_func), pointer   :: coeff_init

    coeff_init => C0_n_m00
    C0 = coeff_template(p2,[m2],rnul,muIR2,d,offset,coeff_init,coeff_init)

  end function C0d_m00_coeff

  function C0d_mmm_table(p2,m2,d,offset) result(C0)
    complex(dp), intent(in)           :: m2
    real(dp),    intent(in)           :: p2,d
    integer,       intent(in), optional :: offset
    complex(dp) :: C0
    procedure(coeff_func), pointer   :: coeff_init,coeff_update

    coeff_init => C0_n_mmm_init
    coeff_update => C0_n_mmm_update
    C0 = coeff_template(p2,[m2],rnul,rnul,d,offset,coeff_init,coeff_init)

  end function C0d_mmm_table

  function C0d_000_coeff(p2,muIR2,d,offset) result(C0)
    real(dp),    intent(in)           :: p2,d,muIR2
    integer,       intent(in), optional :: offset
    complex(dp) :: C0
    procedure(coeff_func), pointer   :: coeff_init

    coeff_init => C0_n_000
    C0 = coeff_template(p2,[cnul],rnul,muIR2,d,offset,coeff_init,coeff_init)

  end function C0d_000_coeff

  function C0d_000_EP1_coeff(p2,muIR2,d,offset) result(C0)
    real(dp),    intent(in)           :: p2,d,muIR2
    integer,       intent(in), optional :: offset
    complex(dp) :: C0
    procedure(coeff_func), pointer   :: coeff_init

    coeff_init => C0_n_000_EP1
    C0 = coeff_template(p2,[cnul],rnul,muIR2,d,offset,coeff_init,coeff_init)

  end function C0d_000_EP1_coeff

  function C0d_m0m1m1_coeff(p2,m02,m12,d,offset) result(C0)
    use c0_m0m1m1_DP, only: C0_n_m0m1m1
    real(dp),    intent(in)           :: p2,d
    complex(dp), intent(in)           :: m02,m12
    integer,       intent(in), optional :: offset
    complex(dp)                       :: C0
    procedure(coeff_func), pointer      :: coeff_init

    coeff_init => C0_n_m0m1m1
    C0 = coeff_template(p2,[m02,m12],rnul,rnul,d,offset,coeff_init,coeff_init)

  end function C0d_m0m1m1_coeff

  function C0d_m0m1m1_list(p2,m02,m12,d,clist) result(C0)
    use c0_m0m1m1_DP, only: C0_n_m0m1m1
    real(dp),    intent(in)           :: p2,d
    complex(dp), intent(in)           :: m02,m12
    complex(dp), intent(in)           :: clist(0:)
    complex(dp)                       :: C0
    procedure(coeff_func), pointer      :: coeff_init

    coeff_init => C0_n_m0m1m1
    C0 = list_template(p2,[m02,m12],rnul,rnul,d,clist,coeff_init,coeff_init)

  end function C0d_m0m1m1_list

end module triangle_expansion_DP

module triangle_expansion_QP
  use c0_0mm_QP, only: C0_n_0mm
  use c0_m00_QP, only: C0_n_m00 
  use c0_mmm_QP, only: C0_n_mmm, C0_n_mmm_init, C0_n_mmm_update
  use c0_000_QP, only: C0_n_000, C0_n_000_EP1
  use b0_mm_QP, only: B0_n_mm, B0_n_mm_init, B0_n_mm_update
  use b0_QP, only: B0_0_n, B0_n
  use b0_m0m1_QP, only: B0_n_m0m1
  use triangle_aux_QP, only: qp, target_precision, cone, cnul, rnul
  
  implicit none

  abstract interface
    function coeff_func(p2,m2,muUV2,muIR2,n)
      use triangle_aux_QP, only: qp
      implicit none
      complex(qp), intent(in) :: p2,m2(:),muUV2,muIR2
      integer,       intent(in) :: n
      complex(qp)             :: coeff_func
    end function coeff_func
  end interface

  contains

  function coeff_template(p2,m2,muUV2,muIR2,d,offset,coeff_init,coeff_update) result(coeff)
    ! Computes the coefficient function func(p2,d,m2,muUV2,muIR2,offset) := 
    ! \sum_{n=offset}^\infty d^{n-offset} func^(n)
    complex(qp),                  intent(in) :: m2(:)
    real(qp),                     intent(in) :: p2,muUV2,muIR2,d
    integer, optional,              intent(in) :: offset
    procedure(coeff_func), pointer, intent(in) :: coeff_init,coeff_update
    complex(qp) :: coeff,coeffn
    integer    :: n,m

    if (d .lt. 0 .or. d .gt. 1) then
      write (*,*) 'ERROR: called coeff_template with invalid d=', d
      stop
    end if

    if (present(offset)) then
      if (offset .lt. 0) then
        write (*,*) 'ERROR: called coeff_template with offset<0'
        stop
      end if
      m = offset
    else
      m = 0
    end if

    n = m
    coeff = coeff_init(cmplx(p2,kind=qp),m2,cmplx(muUV2,kind=qp),cmplx(muIR2,kind=qp),n)

    n = n + 1
    coeffn = d**(n-m)*coeff_update(cmplx(p2,kind=qp),m2,cmplx(muUV2,kind=qp),cmplx(muIR2,kind=qp),n)

    do while (abs(coeffn/coeff) .gt. target_precision)
      coeff = coeff + coeffn
      n = n + 1
      coeffn = d**(n-m)*coeff_update(cmplx(p2,kind=qp),m2,cmplx(muUV2,kind=qp),cmplx(muIR2,kind=qp),n)
    end do
    coeff = coeff + coeffn
    
  end function coeff_template

  function list_template(p2,m2,muUV2,muIR2,d,clist,list_init,list_update) result(list)
    ! Computes \sum_i clist(i)*func(p2,d,m2,muUV2,muIR2,i) in an efficient way,
    ! computing equal coefficients only once.
    complex(qp),                  intent(in) :: m2(:)
    real(qp),                     intent(in) :: p2,muUV2,muIR2,d
    complex(qp),                  intent(in) :: clist(0:)
    procedure(coeff_func), pointer, intent(in) :: list_init,list_update
    complex(qp) :: list,listn
    integer       :: n,m,k

    if (d .lt. 0 .or. d .gt. 1) then
      write (*,*) 'ERROR: called list_template with invalid d=',d
      stop
    end if

    do m = 1, size(clist)
      if (clist(m-1) .ne. 0._qp) then
        exit
      end if
    end do
    if (m .eq. size(clist) .and. clist(m-1) .eq. 0._qp  ) then
      write (*,*) 'ERROR: called list_template with zero clist.'
      stop 
    end if
    m = m - 1

    n = m
    list = clist(m)*list_init(cmplx(p2,kind=qp),m2,cmplx(muUV2,kind=qp),cmplx(muIR2,kind=qp),n)

    n = n + 1
    listn = 0._qp
    do k = m, min(n,size(clist)-1)
      if (clist(k) .ne. 0._qp) then
        listn = listn + d**(n-k)*clist(k)
      end if
    end do
    listn = listn*list_update(cmplx(p2,kind=qp),m2,cmplx(muUV2,kind=qp),cmplx(muIR2,kind=qp),n)

    ! the condition n .lt. size(clist)-1 ensures that for d=0 all coefficients
    ! are taken into account. For instance for a list with a zero
    ! clist=[1.,1.,0.,1.] this would otherwise lead to a mistake.
    do while (abs(listn/list) .gt. target_precision .or. n .lt. size(clist)-1)
      list = list + listn
      n = n + 1
      listn = cnul
      do k = m, min(n,size(clist)-1)
        if (clist(k) .ne. 0._qp) then
          listn = listn + d**(n-k)*clist(k)
        end if
      end do
      listn = listn*list_update(cmplx(p2,kind=qp),m2,cmplx(muUV2,kind=qp),cmplx(muIR2,kind=qp),n)
    end do
    list = list + listn
    
  end function list_template


  !!!!!!!!
  !  B0  !
  !!!!!!!!

  function B0d_coeff(p2,m2,muUV2,d,offset) result(B0)
    complex(qp), intent(in)           :: m2
    real(qp),    intent(in)           :: p2,d,muUV2
    integer,       intent(in), optional :: offset
    complex(qp) :: B0
    procedure(coeff_func), pointer :: coeff_init

    coeff_init => B0_n
    B0 = coeff_template(p2,[m2],muUV2,rnul,d,offset,coeff_init,coeff_init)

  end function B0d_coeff

  function B0d_list(p2,m2,muUV2,d,clist) result(B0)
    complex(qp), intent(in) :: m2
    real(qp),    intent(in) :: p2,d,muUV2
    complex(qp), intent(in) :: clist(0:)
    complex(qp) :: B0
    procedure(coeff_func), pointer :: coeff_init

    coeff_init => B0_n
    B0 = list_template(p2,[m2],muUV2,rnul,d,clist,coeff_init,coeff_init)

  end function B0d_list

  function B0d_mm_opt(p2,m2,muUV2,d,offset) result(B0)
    complex(qp), intent(in)        :: m2
    real(qp),    intent(in)        :: p2,d,muUV2
    integer,    intent(in), optional :: offset
    complex(qp)                    :: B0
    procedure(coeff_func), pointer   :: coeff_init,coeff_update

    coeff_init => B0_n_mm_init
    coeff_update => B0_n_mm_update
    B0 = coeff_template(p2,[m2],muUV2,rnul,d,offset,coeff_init,coeff_init)

  end function B0d_mm_opt

  function B0d_mm_opt_list(p2,m2,muUV2,d,clist) result(B0)
    complex(qp), intent(in)      :: m2
    real(qp),    intent(in)      :: p2,d,muUV2
    complex(qp), intent(in)      :: clist(0:)
    complex(qp)                  :: B0
    procedure(coeff_func), pointer :: coeff_init,coeff_update

    coeff_init => B0_n_mm_init
    coeff_update => B0_n_mm_update
    B0 = list_template(p2,[m2],muUV2,rnul,d,clist,coeff_init,coeff_update)

  end function B0d_mm_opt_list

  function B0d_0_coeff(p2,muUV2,d,offset) result(B0)
    real(qp),    intent(in)           :: p2,d,muUV2
    integer,       intent(in), optional :: offset
    complex(qp)                       :: B0
    procedure(coeff_func),      pointer :: coeff_init

    coeff_init => B0_0_n
    B0 = coeff_template(p2,[cnul],muUV2,rnul,d,offset,coeff_init,coeff_init)
  end function B0d_0_coeff

  function B0d_0_list(p2,muUV2,d,clist) result(B0)
    real(qp),    intent(in)      :: p2,d,muUV2
    complex(qp), intent(in)      :: clist(0:)
    complex(qp)                  :: B0
    procedure(coeff_func), pointer :: coeff_init,coeff_update

    coeff_init => B0_0_n
    B0 = list_template(p2,[cnul],muUV2,rnul,d,clist,coeff_init,coeff_init)

  end function B0d_0_list

  function B0d_m0m1_coeff(p2,m02,m12,muUV2,d,offset) result(B0)
    use b0_m0m1_QP, only: B0_n_m0m1
    real(qp),    intent(in)           :: p2,d,muUV2
    complex(qp), intent(in)           :: m02,m12
    integer,       intent(in), optional :: offset
    complex(qp)                       :: B0
    procedure(coeff_func), pointer      :: coeff_init

    coeff_init => B0_n_m0m1
    B0 = coeff_template(p2,[m02,m12],muUV2,rnul,d,offset,coeff_init,coeff_init)

  end function B0d_m0m1_coeff

  function B0d_m0m1_list(p2,m02,m12,muUV2,d,clist) result(B0)
    use b0_m0m1_QP, only: B0_n_m0m1
    real(qp),    intent(in)      :: p2,d,muUV2
    complex(qp), intent(in)      :: m02,m12
    complex(qp), intent(in)      :: clist(0:)
    complex(qp)                  :: B0
    procedure(coeff_func), pointer :: coeff_init

    coeff_init => B0_n_m0m1
    B0 = list_template(p2,[m02,m12],muUV2,rnul,d,clist,coeff_init,coeff_init)

  end function b0d_m0m1_list

  !!!!!!!!
  !  C0  !
  !!!!!!!!

  function C0d_0mm_coeff(p2,m2,d,offset) result(C0)
    complex(qp), intent(in)           :: m2
    real(qp),    intent(in)           :: p2,d
    integer,       intent(in), optional :: offset
    complex(qp)                       :: C0
    procedure(coeff_func),     pointer  :: coeff_init

    coeff_init => C0_n_0mm
    C0 = coeff_template(p2,[m2],rnul,rnul,d,offset,coeff_init,coeff_init)

  end function C0d_0mm_coeff

  function C0d_m00_coeff(p2,m2,muIR2,d,offset) result(C0)
    complex(qp), intent(in)           :: m2
    real(qp),    intent(in)           :: p2,d,muIR2
    integer,       intent(in), optional :: offset
    complex(qp)                       :: C0
    procedure(coeff_func), pointer   :: coeff_init

    coeff_init => C0_n_m00
    C0 = coeff_template(p2,[m2],rnul,muIR2,d,offset,coeff_init,coeff_init)

  end function C0d_m00_coeff

  function C0d_mmm_table(p2,m2,d,offset) result(C0)
    complex(qp), intent(in)           :: m2
    real(qp),    intent(in)           :: p2,d
    integer,       intent(in), optional :: offset
    complex(qp) :: C0
    procedure(coeff_func), pointer   :: coeff_init,coeff_update

    coeff_init => C0_n_mmm_init
    coeff_update => C0_n_mmm_update
    C0 = coeff_template(p2,[m2],rnul,rnul,d,offset,coeff_init,coeff_init)

  end function C0d_mmm_table

  function C0d_000_coeff(p2,muIR2,d,offset) result(C0)
    real(qp),    intent(in)           :: p2,d,muIR2
    integer,       intent(in), optional :: offset
    complex(qp) :: C0
    procedure(coeff_func), pointer   :: coeff_init

    coeff_init => C0_n_000
    C0 = coeff_template(p2,[cnul],rnul,muIR2,d,offset,coeff_init,coeff_init)

  end function C0d_000_coeff

  function C0d_000_EP1_coeff(p2,muIR2,d,offset) result(C0)
    real(qp),    intent(in)           :: p2,d,muIR2
    integer,       intent(in), optional :: offset
    complex(qp) :: C0
    procedure(coeff_func), pointer   :: coeff_init

    coeff_init => C0_n_000_EP1
    C0 = coeff_template(p2,[cnul],rnul,muIR2,d,offset,coeff_init,coeff_init)

  end function C0d_000_EP1_coeff

  function C0d_m0m1m1_coeff(p2,m02,m12,d,offset) result(C0)
    use c0_m0m1m1_QP, only: C0_n_m0m1m1
    real(qp),    intent(in)           :: p2,d
    complex(qp), intent(in)           :: m02,m12
    integer,       intent(in), optional :: offset
    complex(qp)                       :: C0
    procedure(coeff_func), pointer      :: coeff_init

    coeff_init => C0_n_m0m1m1
    C0 = coeff_template(p2,[m02,m12],rnul,rnul,d,offset,coeff_init,coeff_init)

  end function C0d_m0m1m1_coeff

  function C0d_m0m1m1_list(p2,m02,m12,d,clist) result(C0)
    use c0_m0m1m1_QP, only: C0_n_m0m1m1
    real(qp),    intent(in)           :: p2,d
    complex(qp), intent(in)           :: m02,m12
    complex(qp), intent(in)           :: clist(0:)
    complex(qp)                       :: C0
    procedure(coeff_func), pointer      :: coeff_init

    coeff_init => C0_n_m0m1m1
    C0 = list_template(p2,[m02,m12],rnul,rnul,d,clist,coeff_init,coeff_init)

  end function C0d_m0m1m1_list

end module triangle_expansion_QP
