!   MyModule
!
!   contains functions:
!            subroutines: Growth, nearest_interp_1d
!
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!

  module psModule

     implicit none

!   setup some typical mesh parameters for mesh points
     integer, parameter :: nexp = 12
     integer, parameter :: n = 2**nexp
     integer, parameter :: NMAX=2**nexp
     integer, parameter :: dp = selected_real_kind(15)

!   fixed max number of bessel transforms
!     integer, parameter :: my_nmax = 30
     integer :: my_nmax

     real(kind=dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
     real(kind=dp), parameter :: Omega_M = 0.2648
     real(kind=dp), parameter :: Omega_L = 0.7352

  contains

!================================================================
!                        l_max function
!================================================================

integer function l_max (fixedk)  ! returns the integer as a linear response to input k

	implicit none
	
	real(kind=dp), intent(in):: fixedk
	real(kind=dp) :: l0,le,k0,ke

        l0 = 15
        le = 30
        k0 = 0.0
        ke = 3.0

        l_max = INT(((le-l0)*fixedk+(l0*ke-le*k0))/(ke-k0))

end function l_max

!================================================================
!                        prectent function
!================================================================

subroutine percent(p,t,s,pp)  
! writes p-th precent value out of t total with step of s, given that pp was the last output

	implicit none
	
	integer, intent(in):: p,t,s
	integer, intent(inout):: pp
	integer :: pp_local

        pp_local = (100*p)/t

        if (mod(pp_local,s)==0.and.pp_local>pp)then
          write(6,*), pp_local,"%"
        endif

        pp = pp_local

    return
end subroutine percent

!================================================================
!                        Growth Rate
!================================================================

subroutine Growth(parz,Omega_M,Omega_L,Dz,fz)
  implicit none

  real (kind=dp), intent(in) :: parz, Omega_M, Omega_L
  real (kind=dp), intent(out) :: Dz, fz
  real (kind=dp) :: a, q, p
  real (kind=dp) :: Hub, Hub0
  real (kind=dp) :: beta, beta0

  real (kind=dp) :: alogam
  real (kind=dp) :: beta_log
  real (kind=dp) :: betain
  integer (kind=4) :: ifault

  a    = 1 / (1+parz)
  Hub  = Sqrt(Omega_M*a**(-3) + Omega_L)
  Hub0 = Sqrt(Omega_M + Omega_L)

  q = 5.0_dp/6.0_dp
  p = 2.0_dp/3.0_dp

! -- Note: ifort does not support lgamma function. alogam is explicitly given in asa063.f90. -- !
!  beta_log = lgamma(q)+lgamma(p)-lgamma(q+p)
  beta_log = alogam(q,ifault)+alogam(p,ifault)-alogam(q+p,ifault)
  beta  = 2.0/sqrt(pi)*gamma(q)*gamma(p)*betain(Omega_L*a**3/(Omega_M+Omega_L*a**3),q,p,beta_log,ifault)
  beta0 = 2.0/sqrt(pi)*gamma(q)*gamma(p)*betain(Omega_L/(Omega_M+Omega_L),q,p,beta_log,ifault)

  Dz = Hub/Hub0*beta/beta0
  fz = 3.0_dp*Omega_M**p*Omega_L**q/(a**2*beta*Hub**3) - 3.0_dp/2.0_dp*Omega_M/(Omega_M + a**3*Omega_L)

  return
end subroutine Growth

!====================================================================================================================================!
!            Integrals
!====================================================================================================================================!

subroutine PS_dm(qqin,XXin,YYin,sgin,XXinloop,YYinloop,sginloop,VVin,TTin,fixedk,resPS)

      implicit none

!==== fftlog transform variables ====!

      integer :: dir,kropt
      real(kind=dp) :: mu,q,kr,k,dlnr, dlogr,logkc,logrc,logrmax,logrmin,nc,rk
      real(kind=dp), dimension(NMAX) :: a,kk
      real(kind=dp), dimension(2*NMAX+3*(NMAX/2)+19) :: wsave
      logical :: ok

!==== my variables ====!
      integer :: nindex,i
      real(kind=dp), dimension(n), intent(in) :: qqin,XXin,YYin,sgin
      real(kind=dp), dimension(n), intent(in) :: XXinloop,YYinloop,sginloop
      real(kind=dp), dimension(n), intent(in) :: VVin(1:n),TTin(1:n)
      real(kind=dp), intent(in) :: fixedk
      real(kind=dp), dimension(3), intent(out) :: resPS
      real(kind=dp), dimension(1) :: kout, numout

      logrmin = log10(qqin(1))
      logrmax = log10(qqin(n))
      logrc = (logrmin+logrmax)/2.0_dp
      dlogr = (logrmax - logrmin)/n
      dlnr = dlogr*log(10.0_dp)

      kout(1)=fixedk
      numout(1)=0.0_dp
      resPS=0.0_dp

!======================!
!   fftlog setup

!    order of Bessel function
      mu=0.5_dp
!    bias exponent: q = 0 is unbiased
      q=0.0_dp
!    sensible approximate choice of k_c r_c
      kr=1.0_dp
!    tell fhti to change kr to low-ringing value
!         kropt = 0 to use input kr as is;
!                 1 to change kr to nearest low-ringing kr, quietly;
!                 2 to change kr to nearest low-ringing kr, verbosely;
!                 3 for option to change kr interactively.
      kropt=1
!    forward transform
      dir=1
!     central index (1/2 integral if n is even)
      nc=dble(n+1)/2.0_dp

      my_nmax = l_max(fixedk)

! == Compute PS == !
      do nindex=0,my_nmax,1
         mu = 0.5_dp + dble(nindex)     ! Bessel transform order
!  initialize FFTLog transform - note fhti resets kr
         kr=1.0_dp
         call fhti(n,mu,q,dlnr,kr,kropt,wsave,ok)
!        write (6,*) 'Parameter ok=',ok
!  Input paramaters for the transforms
         logkc=log10(kr)-logrc
         rk=10.0_dp**(logrc-logkc)
!        write (6,*), "kc=",10.0_dp**(logkc),", dlogr=", dlogr
!  set up the k space meshes
         do i = 1,n
           k = 10.0_dp**(logkc+(i-nc)*dlogr)
           kk(i)=k
        enddo
!  call subrutines: Z, A, W

        call int_Z(a,qqin,XXin,YYin,sgin,fixedk,nindex)
        call fht(n,a,dir,wsave)
        call nearest_interp_1d(n,kk,a,1,kout,numout)
        resPS(1) = resPS(1) + 4.0_dp*pi*sqrt(pi/2.d0)*kout(1)**(-3.d0/2.d0)*numout(1)

        call int_A(a,qqin,XXin,YYin,sgin,XXinloop,YYinloop,sginloop,fixedk,nindex)
        call fht(n,a,dir,wsave)
        call nearest_interp_1d(n,kk,a,1,kout,numout)
        resPS(2) = resPS(2) - 2.0_dp*pi*sqrt(pi/2.d0)*kout(1)**(-3.d0/2.d0)*numout(1)

        if (nindex >= 1) then
        call int_W(a,qqin,XXin,YYin,VVin,TTin,fixedk,nindex)
        call fht(n,a,dir,wsave)
        call nearest_interp_1d(n,kk,a,1,kout,numout)
        resPS(3) = resPS(3) + 2.0_dp*pi*sqrt(pi/2.d0)*kout(1)**(-3.d0/2.d0)*numout(1)
        endif

       enddo

end subroutine PS_dm

!====================================================================================================================================!

subroutine PS_dm_eft(qqin,XXin,YYin,sgin,XXin_eft_1,YYin_eft_1,sgin_eft_1,VVin_eft_2,VVin_eft_3,TTin_eft_3,fixedk,resPS)

      implicit none

!==== fftlog transform variables ====!

      integer :: dir,kropt
      real(kind=dp) :: mu,q,kr,k,dlnr, dlogr,logkc,logrc,logrmax,logrmin,nc,rk
      real(kind=dp), dimension(NMAX) :: a,kk
      real(kind=dp), dimension(2*NMAX+3*(NMAX/2)+19) :: wsave
      logical :: ok

!==== my variables ====!
      integer :: nindex,i
      real(kind=dp), dimension(n), intent(in) :: qqin,XXin,YYin,sgin
      real(kind=dp), dimension(n), intent(in) :: XXin_eft_1,YYin_eft_1,sgin_eft_1
      real(kind=dp), dimension(n), intent(in) :: VVin_eft_2,VVin_eft_3,TTin_eft_3
      real(kind=dp), intent(in) :: fixedk
      real(kind=dp), dimension(3), intent(out) :: resPS
      real(kind=dp), dimension(1) :: kout, numout

      logrmin = log10(qqin(1))
      logrmax = log10(qqin(n))
      logrc = (logrmin+logrmax)/2.0_dp
      dlogr = (logrmax - logrmin)/n
      dlnr = dlogr*log(10.0_dp)

      kout(1)=fixedk
      numout(1)=0.0_dp
      resPS=0.0_dp

!======================!
!   fftlog setup

!    order of Bessel function
      mu=0.5_dp
!    bias exponent: q = 0 is unbiased
      q=0.0_dp
!    sensible approximate choice of k_c r_c
      kr=1.0_dp
!    tell fhti to change kr to low-ringing value
!         kropt = 0 to use input kr as is;
!                 1 to change kr to nearest low-ringing kr, quietly;
!                 2 to change kr to nearest low-ringing kr, verbosely;
!                 3 for option to change kr interactively.
      kropt=1
!    forward transform
      dir=1
!     central index (1/2 integral if n is even)
      nc=dble(n+1)/2.0_dp

      my_nmax = l_max(fixedk)

! == Compute PS == !
      do nindex=0,my_nmax,1
         mu = 0.5_dp + dble(nindex)     ! Bessel transform order
!  initialize FFTLog transform - note fhti resets kr
         kr=1.0_dp
         call fhti(n,mu,q,dlnr,kr,kropt,wsave,ok)
!        write (6,*) 'Parameter ok=',ok
!  Input paramaters for the transforms
         logkc=log10(kr)-logrc
         rk=10.0_dp**(logrc-logkc)
!        write (6,*), "kc=",10.0_dp**(logkc),", dlogr=", dlogr
!  set up the k space meshes
         do i = 1,n
           k = 10.0_dp**(logkc+(i-nc)*dlogr)
           kk(i)=k
        enddo
!  call subrutines:

        call int_eft_1(a,qqin,XXin,YYin,sgin,XXin_eft_1,YYin_eft_1,sgin_eft_1,fixedk,nindex)
        call fht(n,a,dir,wsave)
        call nearest_interp_1d(n,kk,a,1,kout,numout)
        resPS(1) = resPS(1) - 2.0_dp*pi*sqrt(pi/2.d0)*kout(1)**(-3.d0/2.d0)*numout(1)

        if (nindex >= 1) then
        call int_eft_2(a,qqin,XXin,YYin,VVin_eft_2,fixedk,nindex)
        call fht(n,a,dir,wsave)
        call nearest_interp_1d(n,kk,a,1,kout,numout)
        resPS(2) = resPS(2) + 2.0_dp*pi*sqrt(pi/2.d0)*kout(1)**(-3.d0/2.d0)*numout(1)

        call int_eft_3(a,qqin,XXin,YYin,VVin_eft_3,TTin_eft_3,fixedk,nindex)
        call fht(n,a,dir,wsave)
        call nearest_interp_1d(n,kk,a,1,kout,numout)
        resPS(3) = resPS(3) + 2.0_dp*pi*sqrt(pi/2.d0)*kout(1)**(-3.d0/2.d0)*numout(1)
        endif

       enddo

end subroutine PS_dm_eft

!====================================================================================================================================!

subroutine PS_bias_df(qqin,XXin,YYin,sgin,Xiinlin,U10inlin,U10insum,U20inloop,U11inloop,XX10inloop,YY10inloop,fixedk,resPS)

      implicit none

!==== fftlog transform variables ====!

      integer :: dir,kropt
      real(kind=dp) :: mu,q,kr,k,dlnr, dlogr,logkc,logrc,logrmax,logrmin,nc,rk
      real(kind=dp), dimension(NMAX) :: a,kk
      real(kind=dp), dimension(2*NMAX+3*(NMAX/2)+19) :: wsave
      logical :: ok

!==== my variables ====!
      integer :: nindex,i
      real(kind=dp), dimension(n), intent(in) :: qqin,XXin,YYin,sgin
      real(kind=dp), dimension(n), intent(in) :: Xiinlin,U10inlin,U10insum
      real(kind=dp), dimension(n), intent(in) :: U20inloop,U11inloop,XX10inloop,YY10inloop
      real(kind=dp), intent(in) :: fixedk
      real(kind=dp), dimension(5), intent(out) :: resPS
      real(kind=dp), dimension(1) :: kout, numout

      logrmin = log10(qqin(1))
      logrmax = log10(qqin(n))
      logrc = (logrmin+logrmax)/2.0_dp
      dlogr = (logrmax - logrmin)/n
      dlnr = dlogr*log(10.0_dp)

      kout(1)=fixedk
      numout(1)=0.0_dp
      resPS=0.0_dp

!======================!
!   fftlog setup

!    order of Bessel function
      mu=0.5_dp
!    bias exponent: q = 0 is unbiased
      q=0.0_dp
!    sensible approximate choice of k_c r_c
      kr=1.0_dp
!    tell fhti to change kr to low-ringing value
!         kropt = 0 to use input kr as is;
!                 1 to change kr to nearest low-ringing kr, quietly;
!                 2 to change kr to nearest low-ringing kr, verbosely;
!                 3 for option to change kr interactively.
      kropt=1
!    forward transform
      dir=1
!     central index (1/2 integral if n is even)
      nc=dble(n+1)/2.0_dp

      my_nmax = l_max(fixedk)

! == Compute PS == !
      do nindex=0,my_nmax,1
         mu = 0.5_dp + dble(nindex)     ! Bessel transform order
!  initialize FFTLog transform - note fhti resets kr
         kr=1.0_dp
         call fhti(n,mu,q,dlnr,kr,kropt,wsave,ok)
!        write (6,*) 'Parameter ok=',ok
!  Input paramaters for the transforms
         logkc=log10(kr)-logrc
         rk=10.0_dp**(logrc-logkc)
!        write (6,*), "kc=",10.0_dp**(logkc),", dlogr=", dlogr
!  set up the k space meshes
         do i = 1,n
           k = 10.0_dp**(logkc+(i-nc)*dlogr)
           kk(i)=k
        enddo

!  call subrutines: \df, \df\df, \df^2, \df^2\df^2, \df\df^2

        call int_d(a,qqin,XXin,YYin,sgin,U10insum,XX10inloop,YY10inloop,fixedk,nindex)
        call fht(n,a,dir,wsave)
        call nearest_interp_1d(n,kk,a,1,kout,numout)
        resPS(1) = resPS(1) + 4.0_dp*pi*sqrt(pi/2.d0)*kout(1)**(-3.d0/2.d0)*numout(1)

        call int_dd(a,qqin,XXin,YYin,Xiinlin,U10inlin,U11inloop,fixedk,nindex)
        call fht(n,a,dir,wsave)
        call nearest_interp_1d(n,kk,a,1,kout,numout)
        resPS(2) = resPS(2) + 4.0_dp*pi*sqrt(pi/2.d0)*kout(1)**(-3.d0/2.d0)*numout(1)

        call int_d2(a,qqin,XXin,YYin,U10inlin,U20inloop,fixedk,nindex)
        call fht(n,a,dir,wsave)
        call nearest_interp_1d(n,kk,a,1,kout,numout)
        resPS(3) = resPS(3) + 4.0_dp*pi*sqrt(pi/2.d0)*kout(1)**(-3.d0/2.d0)*numout(1)

        call int_d2d2(a,qqin,XXin,YYin,Xiinlin,fixedk,nindex)
        call fht(n,a,dir,wsave)
        call nearest_interp_1d(n,kk,a,1,kout,numout)
        resPS(4) = resPS(4) + 4.0_dp*pi*sqrt(pi/2.d0)*kout(1)**(-3.d0/2.d0)*numout(1)

        if (nindex >= 1) then
        call int_dd2(a,qqin,XXin,YYin,Xiinlin,U10inlin,fixedk,nindex)
        call fht(n,a,dir,wsave)
        call nearest_interp_1d(n,kk,a,1,kout,numout)
        resPS(5) = resPS(5) + 4.0_dp*pi*sqrt(pi/2.d0)*kout(1)**(-3.d0/2.d0)*numout(1)
        endif

       enddo

end subroutine PS_bias_df

!====================================================================================================================================!

subroutine PS_bias_s2(qqin,XXin,YYin,sgin,Zetainloop,Chiinloop,V10inloop,V12inloop,XX20inloop,YY20inloop,fixedk,resPS)

      implicit none

!==== fftlog transform variables ====!

      integer :: dir,kropt
      real(kind=dp) :: mu,q,kr,k,dlnr, dlogr,logkc,logrc,logrmax,logrmin,nc,rk
      real(kind=dp), dimension(NMAX) :: a,kk
      real(kind=dp), dimension(2*NMAX+3*(NMAX/2)+19) :: wsave
      logical :: ok

!==== my variables ====!
      integer :: nindex,i
      real(kind=dp), dimension(n), intent(in) :: qqin,XXin,YYin,sgin
      real(kind=dp), dimension(n), intent(in) :: Zetainloop,Chiinloop
      real(kind=dp), dimension(n), intent(in) :: V10inloop,V12inloop,XX20inloop,YY20inloop
      real(kind=dp), intent(in) :: fixedk
      real(kind=dp), dimension(4), intent(out) :: resPS
      real(kind=dp), dimension(1) :: kout, numout

      logrmin = log10(qqin(1))
      logrmax = log10(qqin(n))
      logrc = (logrmin+logrmax)/2.0_dp
      dlogr = (logrmax - logrmin)/n
      dlnr = dlogr*log(10.0_dp)

      kout(1)=fixedk
      numout(1)=0.0_dp
      resPS=0.0_dp

!======================!
!   fftlog setup

!    order of Bessel function
      mu=0.5_dp
!    bias exponent: q = 0 is unbiased
      q=0.0_dp
!    sensible approximate choice of k_c r_c
      kr=1.0_dp
!    tell fhti to change kr to low-ringing value
!         kropt = 0 to use input kr as is;
!                 1 to change kr to nearest low-ringing kr, quietly;
!                 2 to change kr to nearest low-ringing kr, verbosely;
!                 3 for option to change kr interactively.
      kropt=1
!    forward transform
      dir=1
!     central index (1/2 integral if n is even)
      nc=dble(n+1)/2.0_dp

      my_nmax = l_max(fixedk)

! == Compute PS == !
      do nindex=0,my_nmax,1
         mu = 0.5_dp + dble(nindex)     ! Bessel transform order
!  initialize FFTLog transform - note fhti resets kr
         kr=1.0_dp
         call fhti(n,mu,q,dlnr,kr,kropt,wsave,ok)
!        write (6,*) 'Parameter ok=',ok
!  Input paramaters for the transforms
         logkc=log10(kr)-logrc
         rk=10.0_dp**(logrc-logkc)
!        write (6,*), "kc=",10.0_dp**(logkc),", dlogr=", dlogr
!  set up the k space meshes
         do i = 1,n
           k = 10.0_dp**(logkc+(i-nc)*dlogr)
           kk(i)=k
        enddo

!  call subrutines: s^2, \df s^2, \df^2 s^2, s^2s^2

        call int_d(a,qqin,XXin,YYin,sgin,V10inloop,XX20inloop,YY20inloop,fixedk,nindex)
        call fht(n,a,dir,wsave)
        call nearest_interp_1d(n,kk,a,1,kout,numout)
        resPS(1) = resPS(1) + 4.0_dp*pi*sqrt(pi/2.d0)*kout(1)**(-3.d0/2.d0)*numout(1)

        if (nindex >= 1) then
        call int_ds2(a,qqin,XXin,YYin,V12inloop,fixedk,nindex)
        call fht(n,a,dir,wsave)
        call nearest_interp_1d(n,kk,a,1,kout,numout)
        resPS(2) = resPS(2) + 4.0_dp*pi*sqrt(pi/2.d0)*kout(1)**(-3.d0/2.d0)*numout(1)
        endif

        call int_s2s2(a,qqin,XXin,YYin,Chiinloop,fixedk,nindex)
        call fht(n,a,dir,wsave)
        call nearest_interp_1d(n,kk,a,1,kout,numout)
        resPS(3) = resPS(3) + 4.0_dp*pi*sqrt(pi/2.d0)*kout(1)**(-3.d0/2.d0)*numout(1)

        call int_s2s2(a,qqin,XXin,YYin,Zetainloop,fixedk,nindex)
        call fht(n,a,dir,wsave)
        call nearest_interp_1d(n,kk,a,1,kout,numout)
        resPS(4) = resPS(4) + 4.0_dp*pi*sqrt(pi/2.d0)*kout(1)**(-3.d0/2.d0)*numout(1)

       enddo

end subroutine PS_bias_s2

!====================================================================================================================================!

subroutine PS_bias_hd(qqin,XXin,YYin,sgin,dXiinlin,dU10inlin,fixedk,resPS)

      implicit none

!==== fftlog transform variables ====!

      integer :: dir,kropt
      real(kind=dp) :: mu,q,kr,k,dlnr, dlogr,logkc,logrc,logrmax,logrmin,nc,rk
      real(kind=dp), dimension(NMAX) :: a,kk
      real(kind=dp), dimension(2*NMAX+3*(NMAX/2)+19) :: wsave
      logical :: ok

!==== my variables ====!
      integer :: nindex,i
      real(kind=dp), dimension(n), intent(in) :: qqin,XXin,YYin,sgin
      real(kind=dp), dimension(n), intent(in) :: dXiinlin,dU10inlin
      real(kind=dp), intent(in) :: fixedk
      real(kind=dp), dimension(2), intent(out) :: resPS
      real(kind=dp), dimension(1) :: kout, numout

      logrmin = log10(qqin(1))
      logrmax = log10(qqin(n))
      logrc = (logrmin+logrmax)/2.0_dp
      dlogr = (logrmax - logrmin)/n
      dlnr = dlogr*log(10.0_dp)

      kout(1)=fixedk
      numout(1)=0.0_dp
      resPS=0.0_dp

!======================!
!   fftlog setup

!    order of Bessel function
      mu=0.5_dp
!    bias exponent: q = 0 is unbiased
      q=0.0_dp
!    sensible approximate choice of k_c r_c
      kr=1.0_dp
!    tell fhti to change kr to low-ringing value
!         kropt = 0 to use input kr as is;
!                 1 to change kr to nearest low-ringing kr, quietly;
!                 2 to change kr to nearest low-ringing kr, verbosely;
!                 3 for option to change kr interactively.
      kropt=1
!    forward transform
      dir=1
!     central index (1/2 integral if n is even)
      nc=dble(n+1)/2.0_dp

      my_nmax = l_max(fixedk)

! == Compute PS == !
      do nindex=0,my_nmax,1
         mu = 0.5_dp + dble(nindex)     ! Bessel transform order
!  initialize FFTLog transform - note fhti resets kr
         kr=1.0_dp
         call fhti(n,mu,q,dlnr,kr,kropt,wsave,ok)
!        write (6,*) 'Parameter ok=',ok
!  Input paramaters for the transforms
         logkc=log10(kr)-logrc
         rk=10.0_dp**(logrc-logkc)
!        write (6,*), "kc=",10.0_dp**(logkc),", dlogr=", dlogr
!  set up the k space meshes
         do i = 1,n
           k = 10.0_dp**(logkc+(i-nc)*dlogr)
           kk(i)=k
        enddo

!  call subrutines: d^2\df, \dfd^2\df

        if (nindex >= 1) then
        call int_ds2(a,qqin,XXin,YYin,dU10inlin,fixedk,nindex)
        call fht(n,a,dir,wsave)
        call nearest_interp_1d(n,kk,a,1,kout,numout)
        resPS(1) = resPS(1) + 4.0_dp*pi*sqrt(pi/2.d0)*kout(1)**(-3.d0/2.d0)*numout(1)
        endif

        call int_s2s2(a,qqin,XXin,YYin,dXiinlin,fixedk,nindex)
        call fht(n,a,dir,wsave)
        call nearest_interp_1d(n,kk,a,1,kout,numout)
        resPS(2) = resPS(2) + 4.0_dp*pi*sqrt(pi/2.d0)*kout(1)**(-3.d0/2.d0)*numout(1)

       enddo

end subroutine PS_bias_hd

!====================================================================================================================================!
!            Integrands for dark matter
!====================================================================================================================================!

subroutine int_Z(ya,rr,XXin,YYin,sgin,fixedk,nin)

      implicit none

      integer, intent(in) :: nin
      real(kind=dp), intent(in) :: fixedk
      real(kind=dp), intent(in) :: rr(1:n),XXin(1:n),YYin(1:n),sgin(1:n)
      real(kind=dp), intent(out) :: ya(1:n)
      real(kind=dp) :: ntemp1,ntemp2

! Two cases when n=0 and n>0
     if (nin==0) then
        ya(1:n) = rr(1:n)**(3.d0/2.d0)*(exp(-1.d0/2.d0*fixedk**2*(XXin(1:n)+YYin(1:n))) - exp(-fixedk**2*sgin(1:n)))
     else
        ntemp1 = dble(nin)
        ntemp2 = 3.d0/2.d0 - ntemp1 
        ya(1:n) = rr(1:n)**ntemp2*(fixedk*YYin(1:n))**ntemp1*exp(-1.d0/2.d0*fixedk**2.d0*(XXin(1:n)+YYin(1:n)))
    endif   

end subroutine int_Z

!============================================

subroutine int_A(ya,rr,XXin,YYin,sgin,XXinloop,YYinloop,sginloop,fixedk,nin)

      implicit none

      integer, intent(in) :: nin
      real(kind=dp), intent(in) :: fixedk
      real(kind=dp), intent(in) :: rr(1:n),XXin(1:n),YYin(1:n),sgin(1:n)
      real(kind=dp), intent(in) :: XXinloop(1:n),YYinloop(1:n),sginloop(1:n)
      real(kind=dp), intent(out) :: ya(1:n)
      real(kind=dp) :: ntemp1,ntemp2

! Two cases when n=0 and n>0
      if (nin==0) then
        ya(1:n) = rr(1:n)**(3.d0/2.d0)*(fixedk**2*(XXinloop(1:n)+YYinloop(1:n))&
                  *exp(-1.d0/2.d0*fixedk**2*(XXin(1:n)+YYin(1:n)))&
                  - 2.d0*fixedk**2*sginloop(1:n)*exp(-fixedk**2*sgin(1:n)))
      else
        ntemp1 = dble(nin)
        ntemp2 = 3.d0/2.d0 - ntemp1
        ya(1:n) = rr(1:n)**ntemp2*(fixedk*YYin(1:n))**ntemp1*(fixedk**2*XXinloop(1:n)+fixedk**2*YYinloop(1:n)&
                  -2.d0*ntemp1*YYinloop(1:n)/YYin(1:n))*exp(-1.d0/2.d0*fixedk**2.d0*(XXin(1:n)+YYin(1:n)))
      endif

end subroutine int_A

!============================================

subroutine int_W(ya,rr,XXin,YYin,VVin,TTin,fixedk,nin)

       implicit none

       integer, intent(in) :: nin
       real(kind=dp), intent(in) :: fixedk
       real(kind=dp), intent(in) :: rr(1:n),XXin(1:n),YYin(1:n)
       real(kind=dp), intent(in) :: VVin(1:n),TTin(1:n)
       real(kind=dp), intent(out) :: ya(1:n)
       real(kind=dp) :: ntemp1,ntemp2

! Two cases when n=1 and n>1
       if (nin==1) then
          ya(1:n) = rr(1:n)**(3.d0/2.d0)*fixedk**3*(VVin(1:n)+(1.d0/3.d0)*TTin(1:n))*exp(-1.d0/2.d0*fixedk**2*(XXin(1:n)+YYin(1:n)))
       else
          ntemp1 = dble(nin) - 1.d0
          ntemp2 = 3.d0/2.d0 - ntemp1
          ya(1:n) = rr(1:n)**ntemp2*(fixedk*YYin(1:n))**ntemp1*(fixedk**3*VVin(1:n)+fixedk**3*(1.d0/3.d0)*TTin(1:n)&
                   -(2.d0/3.d0)*ntemp1*fixedk*TTin(1:n)/YYin(1:n))*exp(-1.d0/2.d0*fixedk**2.d0*(XXin(1:n)+YYin(1:n)))
       endif

end subroutine int_W

!====================================================================================================================================!
!            Integrands for eft
!====================================================================================================================================!

subroutine int_eft_1(ya,rr,XXin,YYin,sgin,XXin_eft_1,YYin_eft_1,sgin_eft_1,fixedk,nin)

      implicit none

      integer, intent(in) :: nin
      real(kind=dp), intent(in) :: fixedk
      real(kind=dp), intent(in) :: rr(1:n),XXin(1:n),YYin(1:n),sgin(1:n)
      real(kind=dp), intent(in) :: XXin_eft_1(1:n),YYin_eft_1(1:n),sgin_eft_1(1:n)
      real(kind=dp), intent(out) :: ya(1:n)
      real(kind=dp) :: ntemp1,ntemp2

! Two cases when n=0 and n>0
      if (nin==0) then
        ya(1:n) = rr(1:n)**(3.d0/2.d0)*(fixedk**2*(XXin_eft_1(1:n)+YYin_eft_1(1:n))&
                  *exp(-1.d0/2.d0*fixedk**2*(XXin(1:n)+YYin(1:n)))&
                  - 2.d0*fixedk**2*sgin_eft_1(1:n)*exp(-fixedk**2*sgin(1:n)))
      else
        ntemp1 = dble(nin)
        ntemp2 = 3.d0/2.d0 - ntemp1
        ya(1:n) = rr(1:n)**ntemp2*(fixedk*YYin(1:n))**ntemp1*(fixedk**2*XXin_eft_1(1:n)+fixedk**2*YYin_eft_1(1:n)&
                  -2.d0*ntemp1*YYin_eft_1(1:n)/YYin(1:n))*exp(-1.d0/2.d0*fixedk**2.d0*(XXin(1:n)+YYin(1:n)))
      endif

end subroutine int_eft_1

!============================================

subroutine int_eft_2(ya,rr,XXin,YYin,VVin_eft_2,fixedk,nin)

       implicit none

       integer, intent(in) :: nin
       real(kind=dp), intent(in) :: fixedk
       real(kind=dp), intent(in) :: rr(1:n),XXin(1:n),YYin(1:n)
       real(kind=dp), intent(in) :: VVin_eft_2(1:n)
       real(kind=dp), intent(out) :: ya(1:n)
       real(kind=dp) :: ntemp1,ntemp2

! Two cases when n=1 and n>1
       if (nin==1) then
          ya(1:n) = rr(1:n)**(3.d0/2.d0)*fixedk**3*VVin_eft_2(1:n)*exp(-1.d0/2.d0*fixedk**2*(XXin(1:n)+YYin(1:n)))
       else
          ntemp1 = dble(nin) - 1.d0
          ntemp2 = 3.d0/2.d0 - ntemp1
          ya(1:n) = rr(1:n)**ntemp2*(fixedk*YYin(1:n))**ntemp1*fixedk**3*VVin_eft_2(1:n)*exp(-1.d0/2.d0*fixedk**2.d0*(XXin(1:n)+YYin(1:n)))
       endif

end subroutine int_eft_2

!============================================

subroutine int_eft_3(ya,rr,XXin,YYin,VVin_eft_3,TTin_eft_3,fixedk,nin)

       implicit none

       integer, intent(in) :: nin
       real(kind=dp), intent(in) :: fixedk
       real(kind=dp), intent(in) :: rr(1:n),XXin(1:n),YYin(1:n)
       real(kind=dp), intent(in) :: VVin_eft_3(1:n),TTin_eft_3(1:n)
       real(kind=dp), intent(out) :: ya(1:n)
       real(kind=dp) :: ntemp1,ntemp2

! Two cases when n=1 and n>1
       if (nin==1) then
          ya(1:n) = rr(1:n)**(3.d0/2.d0)*fixedk**3*(VVin_eft_3(1:n)+(1.d0/3.d0)*TTin_eft_3(1:n))*exp(-1.d0/2.d0*fixedk**2*(XXin(1:n)+YYin(1:n)))
       else
          ntemp1 = dble(nin) - 1.d0
          ntemp2 = 3.d0/2.d0 - ntemp1
          ya(1:n) = rr(1:n)**ntemp2*(fixedk*YYin(1:n))**ntemp1*(fixedk**3*VVin_eft_3(1:n)+fixedk**3*(1.d0/3.d0)*TTin_eft_3(1:n)&
                   -(2.d0/3.d0)*ntemp1*fixedk*TTin_eft_3(1:n)/YYin(1:n))*exp(-1.d0/2.d0*fixedk**2.d0*(XXin(1:n)+YYin(1:n)))
       endif

end subroutine int_eft_3

!====================================================================================================================================!
!            Integrands for biasing
!====================================================================================================================================!

subroutine int_d(ya,rr,XXin,YYin,sgin,U10insum,XX10inloop,YY10inloop,fixedk,nin)

      implicit none

      integer, intent(in) :: nin
      integer :: ncut
      real(kind=dp), intent(in) :: fixedk
      real(kind=dp), intent(in) :: rr(1:n),XXin(1:n),YYin(1:n),sgin(1:n)
      real(kind=dp), intent(in) :: U10insum(1:n),XX10inloop(1:n),YY10inloop(1:n)
      real(kind=dp), intent(out) :: ya(1:n)
      real(kind=dp) :: ntemp1,ntemp2
      real(kind=dp) :: sig10(1:n)

      ncut = 50
      sig10(1:n) = 1.d0/2.d0*sum(XX10inloop((n-ncut):n))/dble(ncut+1)*rr(1:n)/rr(1:n)

! Two cases when n=0 and n>0
      if (nin==0) then
        ya(1:n) = rr(1:n)**(3.d0/2.d0)*(-1)*(fixedk**2*(XX10inloop(1:n)+YY10inloop(1:n))&
                  *exp(-1.d0/2.d0*fixedk**2*(XXin(1:n)+YYin(1:n)))&
                  - 2.d0*fixedk**2*sig10(1:n)*exp(-fixedk**2*sgin(1:n)))
      else
        ntemp1 = dble(nin)
        ntemp2 = 3.d0/2.d0 - ntemp1
        ya(1:n) = rr(1:n)**ntemp2*(fixedk*YYin(1:n))**ntemp1*(&
                  -fixedk**2*(XX10inloop(1:n)+YY10inloop(1:n))&
                  +2.d0*(ntemp1*YY10inloop(1:n)-rr(1:n)*U10insum(1:n))/YYin(1:n))&
                  *exp(-1.d0/2.d0*fixedk**2.d0*(XXin(1:n)+YYin(1:n)))
      endif

end subroutine int_d

!====================================================================================================================================!

subroutine int_dd(ya,rr,XXin,YYin,Xiinlin,U10inlin,U11inloop,fixedk,nin)

      implicit none

      integer, intent(in) :: nin
      real(kind=dp), intent(in) :: fixedk
      real(kind=dp), intent(in) :: rr(1:n),XXin(1:n),YYin(1:n)
      real(kind=dp), intent(in) :: Xiinlin(1:n),U10inlin(1:n),U11inloop(1:n)
      real(kind=dp), intent(out) :: ya(1:n)
      real(kind=dp) :: ntemp1,ntemp2

! Two cases when n=0 and n>0
      if (nin==0) then
        ya(1:n) = rr(1:n)**(3.d0/2.d0)*(Xiinlin(1:n)-fixedk**2*(U10inlin(1:n))**2)&
                  *exp(-1.d0/2.d0*fixedk**2*(XXin(1:n)+YYin(1:n)))
      else
        ntemp1 = dble(nin)
        ntemp2 = 3.d0/2.d0 - ntemp1
        ya(1:n) = rr(1:n)**ntemp2*(fixedk*YYin(1:n))**ntemp1*(Xiinlin(1:n)-fixedk**2*(U10inlin(1:n))**2&
                    +(2.d0*ntemp1*(U10inlin(1:n))**2 - rr(1:n)*U11inloop(1:n))/YYin(1:n))&
                    *exp(-1.d0/2.d0*fixedk**2.d0*(XXin(1:n)+YYin(1:n)))
      endif

end subroutine int_dd

!====================================================================================================================================!

subroutine int_d2(ya,rr,XXin,YYin,U10inlin,U20inloop,fixedk,nin)

      implicit none

      integer, intent(in) :: nin
      real(kind=dp), intent(in) :: fixedk
      real(kind=dp), intent(in) :: rr(1:n),XXin(1:n),YYin(1:n)
      real(kind=dp), intent(in) :: U10inlin(1:n),U20inloop(1:n)
      real(kind=dp), intent(out) :: ya(1:n)
      real(kind=dp) :: ntemp1,ntemp2

! Two cases when n=0 and n>0
      if (nin==0) then
        ya(1:n) = rr(1:n)**(3.d0/2.d0)*(-fixedk**2*(U10inlin(1:n))**2)*exp(-1.d0/2.d0*fixedk**2*(XXin(1:n)+YYin(1:n)))
      else
        ntemp1 = dble(nin)
        ntemp2 = 3.d0/2.d0 - ntemp1
        ya(1:n) = rr(1:n)**ntemp2*(fixedk*YYin(1:n))**ntemp1*(-fixedk**2*(U10inlin(1:n))**2&
                 +(2.d0*ntemp1*(U10inlin(1:n))**2-rr(1:n)*U20inloop(1:n))/YYin(1:n))&
                 *exp(-1.d0/2.d0*fixedk**2.d0*(XXin(1:n)+YYin(1:n)))
      endif

end subroutine int_d2

!====================================================================================================================================!

subroutine int_d2d2(ya,rr,XXin,YYin,Xiinlin,fixedk,nin)

      implicit none

      integer, intent(in) :: nin
      real(kind=dp), intent(in) :: fixedk
      real(kind=dp), intent(in) :: rr(1:n),XXin(1:n),YYin(1:n)
      real(kind=dp), intent(in) :: Xiinlin(1:n)
      real(kind=dp), intent(out) :: ya(1:n)
      real(kind=dp) :: ntemp1,ntemp2

! Two cases when n=0 and n>0
      if (nin==0) then
        ya(1:n) = rr(1:n)**(3.d0/2.d0)*2.d0*(Xiinlin(1:n))**2*exp(-1.d0/2.d0*fixedk**2*(XXin(1:n)+YYin(1:n)))
      else
        ntemp1 = dble(nin)
        ntemp2 = 3.d0/2.d0 - ntemp1
        ya(1:n) = rr(1:n)**ntemp2*(fixedk*YYin(1:n))**ntemp1*2.d0*(Xiinlin(1:n))**2&
                  *exp(-1.d0/2.d0*fixedk**2.d0*(XXin(1:n)+YYin(1:n)))
      endif

end subroutine int_d2d2

!====================================================================================================================================!

subroutine int_dd2(ya,rr,XXin,YYin,Xiinlin,U10inlin,fixedk,nin)

      implicit none

      integer, intent(in) :: nin
      real(kind=dp), intent(in) :: fixedk
      real(kind=dp), intent(in) :: rr(1:n),XXin(1:n),YYin(1:n)
      real(kind=dp), intent(in) :: Xiinlin(1:n),U10inlin(1:n)
      real(kind=dp), intent(out) :: ya(1:n)
      real(kind=dp) :: ntemp1,ntemp2

! Two cases when n>=1

        ntemp1 = dble(nin)
        ntemp2 = 3.d0/2.d0 - ntemp1
        ya(1:n) = rr(1:n)**ntemp2*(fixedk*YYin(1:n))**ntemp1*(-2.d0*Xiinlin(1:n)*rr(1:n)*U10inlin(1:n)/YYin(1:n))&
                  *exp(-1.d0/2.d0*fixedk**2.d0*(XXin(1:n)+YYin(1:n)))

end subroutine int_dd2

!====================================================================================================================================!

subroutine int_ds2(ya,rr,XXin,YYin,V12inloop,fixedk,nin)

      implicit none

      integer, intent(in) :: nin
      real(kind=dp), intent(in) :: fixedk
      real(kind=dp), intent(in) :: rr(1:n),XXin(1:n),YYin(1:n)
      real(kind=dp), intent(in) :: V12inloop(1:n)
      real(kind=dp), intent(out) :: ya(1:n)
      real(kind=dp) :: ntemp1,ntemp2

! Two cases when n>=1

        ntemp1 = dble(nin)
        ntemp2 = 3.d0/2.d0 - ntemp1
        ya(1:n) = rr(1:n)**ntemp2*(fixedk*YYin(1:n))**ntemp1*(-rr(1:n)*V12inloop(1:n)/YYin(1:n))&
                  *exp(-1.d0/2.d0*fixedk**2.d0*(XXin(1:n)+YYin(1:n)))

end subroutine int_ds2

!====================================================================================================================================!

subroutine int_s2s2(ya,rr,XXin,YYin,Chiinloop,fixedk,nin)

      implicit none

      integer, intent(in) :: nin
      real(kind=dp), intent(in) :: fixedk
      real(kind=dp), intent(in) :: rr(1:n),XXin(1:n),YYin(1:n)
      real(kind=dp), intent(in) :: Chiinloop(1:n)
      real(kind=dp), intent(out) :: ya(1:n)
      real(kind=dp) :: ntemp1,ntemp2

! Two cases when n=0 and n>0
      if (nin==0) then
        ya(1:n) = rr(1:n)**(3.d0/2.d0)*Chiinloop(1:n)*exp(-1.d0/2.d0*fixedk**2*(XXin(1:n)+YYin(1:n)))

      else
        ntemp1 = dble(nin)
        ntemp2 = 3.d0/2.d0 - ntemp1
        ya(1:n) = rr(1:n)**ntemp2*(fixedk*YYin(1:n))**ntemp1*Chiinloop(1:n)&
                  *exp(-1.d0/2.d0*fixedk**2.d0*(XXin(1:n)+YYin(1:n)))
      endif

end subroutine int_s2s2

!==================================================================
!    NEAREST_INTERP_1D evaluates the nearest neighbor interpolant.
!==================================================================

subroutine nearest_interp_1d ( nd, xd, yd, ni, xi, yi )

!  Discussion:
!    The nearest neighbor interpolant L(ND,XD,YD)(X) is the piecewise
!    constant function which interpolates the data (XD(I),YD(I)) for I = 1
!    to ND.
!  Licensing:  This code is distributed under the GNU LGPL license.
!  Modified:   04 September 2012
!  Author:     John Burkardt
!  Modified:   Zvonimir Vlah 18.04.2013.
!  Parameters:
!    Input, integer ( kind = 4 ) ND, the number of data points.
!    ND must be at least 1.
!    Input, real ( kind = 8 ) XD(ND), the data points.
!    Input, real ( kind = 8 ) YD(ND), the data values.
!    Input, integer ( kind = 4 ) NI, the number of interpolation points.
!    Input, real ( kind = 8 ) XI(NI), the interpolation points.
!    Output, real ( kind = 8 ) YI(NI), the interpolated values.

  implicit none
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  real ( kind = 8 ) d
  real ( kind = 8 ) d2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
!  real ( kind = 8 ) xdin(nd)
!  real ( kind = 8 ) xiin(ni)
!  real ( kind = 8 ) ydin(nd)
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) xi(ni)
  real ( kind = 8 ) yd(nd)
  real ( kind = 8 ) yi(ni)

!  xd(1:nd) = log10(xdin(1:nd))
!  xi(1:ni) = log10(xiin(1:ni))
!  yd(1:nd) = log10(ydin(1:nd))

  do i = 1, ni

    k = 1
    d = abs ( xi(i) - xd(k) )

    do j = 2, nd

      d2 = abs ( xi(i) - xd(j) )

      if ( d2 < d ) then
        k = j
        d = d2
      end if
    end do
       
    if ( xi(i) > xd(k) ) then
    l=k+1
    yi(i) = yd(l) + ( yd(k) - yd(l) )/( xd(k) - xd(l) )*(xi(i)-xd(l))
    else if ( xi(i) < xd(k) ) then
    l=k-1
    yi(i) = yd(l) + ( yd(k) - yd(l) )/( xd(k) - xd(l) )*(xi(i)-xd(l))
    else
    yi(i) = yd(k)
    end if

   end do

  return
end subroutine nearest_interp_1d

end module psModule
