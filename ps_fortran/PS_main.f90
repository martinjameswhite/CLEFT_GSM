!===================================================!

program PS_main

!       Program computes the power spectrum terms in LEFT formalism.
!       Program uses the spherical Hankel transform routine.

  use psModule
  implicit none

!==== My variables ====!

  integer, parameter :: numkk=400
  integer :: i,j,l,p,nz,pp
  real(kind=dp) :: dlnr, dlogr,logkc,logrc,logrmax,logrmin,nc,r
  real(kind=dp) :: Dz, fz
  real(kind=dp) :: kkoutMin,kkoutMax,dkk,cfkkout,fixedk
  real(kind=dp) :: k_val
  real(kind=dp), dimension(3) :: P_dm_val
  real(kind=dp), dimension(5) :: P_bias_df_val
  real(kind=dp), dimension(4) :: P_bias_s2_val
  real(kind=dp), dimension(2) :: P_bias_hd_val
  real(kind=dp), dimension(n,24) :: in
  real(kind=dp), dimension(n) :: rr
  real(kind=dp), dimension(n) :: qqin,XXin,YYin,X0in,sgin
  real(kind=dp), dimension(n) :: XXinloop,YYinloop,X0inloop,sginloop
  real(kind=dp), dimension(n) :: VVin,TTin
  real(kind=dp), dimension(n) :: Xiinlin,U10inlin,U10insum
!  real(kind=dp), dimension(n) :: Xiinlin,U10inlin,U10inloop
  real(kind=dp), dimension(n) :: U20inloop,U11inloop,XX10inloop,YY10inloop
  real(kind=dp), dimension(n) :: Zetainloop,Chiinloop
  real(kind=dp), dimension(n) :: V10inloop,V12inloop,XX20inloop,YY20inloop
  real(kind=dp), dimension(n) :: dXiinlin,dU10inlin
  real(kind=dp), dimension(numkk) :: kkout

  real(kind=dp), dimension(9) :: parz = (/0.00,0.25,0.55,0.60,0.65,0.70,0.75,0.80,1.00/)

  character(len=1024) :: filename
  character(len=3), dimension(9) :: zstr = (/ "000","025","055","060","065","070","075","080","100"/)


!-----------------------------------------------------!
!   setup some typical mesh parameters for mesh points

  nc=dble(n+1)/2.0_dp
  logrmin = log10(0.0001_dp)
  logrmax = log10(1000000.0_dp)           
  logrc = (logrmin+logrmax)/2.0_dp
  dlogr = (logrmax - logrmin)/n
  dlnr = dlogr*log(10.0_dp)

!   set up the r space meshes
  do i = 1,n
     r = 10.0_dp**(logrc+(i-nc)*dlogr)
     rr(i) = r
  enddo

!---     input-data    ---!
! Read input data --- Note that data was prepared to fit the initialized rr(i) mash!

  open(20,file='Ens_ANL15_FFT_2p12.dat',status='old')
  read(20,*)
  read(20,*)
  read(20,*)
    do j=1,n
     read(20,*) in(j,:)
    enddo
  close(20)

!---   redshift loop   ---!
! do nz=1,1,1
 do nz=1,9,1

!---   growth  rate    ---!
  Dz = 0.0_dp
  fz = 0.0_dp
  call Growth(parz(nz),Omega_M,Omega_L,Dz,fz)

  write (6,*), "Redshifr  z = ", parz(nz)
  write (6,*), "Growth D(z) = ", Dz
  write (6,*), "Log growth f(z) = ", fz

!---     input-data    ---!

  qqin = in(:,1)
  sgin = Dz**2.0_dp * in(:,2)/3.0_dp
  X0in = Dz**2.0_dp * 2.0_dp/3.0_dp * (        - in(:,3) - in(:,4))
  XXin = Dz**2.0_dp * 2.0_dp/3.0_dp * (in(:,2) - in(:,3) - in(:,4))
  YYin = Dz**2.0_dp * 2.0_dp * in(:,4)

  sginloop = Dz**4.0_dp * in(:,5)/3.0_dp
  X0inloop = Dz**4.0_dp * 2.0_dp/3.0_dp * (        - in(:,6) - in(:,7))
  XXinloop = Dz**4.0_dp * 2.0_dp/3.0_dp * (in(:,5) - in(:,6) - in(:,7))
  YYinloop = Dz**4.0_dp * 2.0_dp * in(:,7)

  VVin = Dz**4.0_dp * 1.0_dp/5.0_dp* ( 2.0_dp*in(:,8) - 3.0_dp*in(:,9) )
  TTin = Dz**4.0_dp * 3.0_dp* in(:,9)

  Xiinlin = Dz**2.0_dp * in(:,10)
  Zetainloop = Dz**4.0_dp * in(:,11)
  Chiinloop = Dz**4.0_dp * in(:,12)

  U10inlin  = Dz**2.0_dp * in(:,13)
!  U10inloop = Dz**4.0_dp * in(:,12)
  U10insum  = Dz**2.0_dp * in(:,13) + Dz**4.0_dp * in(:,14)

  U11inloop  = Dz**4.0_dp * in(:,15)
  U20inloop  = Dz**4.0_dp * in(:,16)
  V10inloop  = Dz**4.0_dp * in(:,17)
  V12inloop  = Dz**4.0_dp * in(:,18)
  XX10inloop = Dz**4.0_dp * in(:,19)
  YY10inloop = Dz**4.0_dp * in(:,20)
  XX20inloop = Dz**4.0_dp * in(:,21)
  YY20inloop = Dz**4.0_dp * in(:,22)

  dXiinlin = Dz**4.0_dp * in(:,23)
  dU10inlin = Dz**4.0_dp * in(:,24)

!---------------------------------------------------------!
!   q-data check if input data and the mash are consistent

!    do j=1,n,15
!     write (6,*), "d=",qqin(j),rr(j),abs(qqin(j)-rr(j))/rr(j)
!    enddo
!    stop

!---     out k-data grid    ---!

  kkoutMin=0.001_dp
  kkoutMax=3.0_dp
  kkout(1) = kkoutMin
! log k-mash
  dkk = (log10(kkoutMax)-log10(kkoutMin))/(numkk-1)
  cfkkout = 10.0_dp**(dkk)
  do l = 2,numkk
     kkout(l) = cfkkout*kkout(l-1)
  enddo
! lin k-mash
!  dkk = (kkoutMax-kkoutMin)/(numkk-1)
!  do l = 2,numkk
!     kkout(l) = kkout(l-1)+dkk
!  enddo

!---  write into the file  ---!

    write(*,*) 'Power Spectrum Terms'
    write (filename,"(A10,A3,A4)") "PS_hh_46_z",zstr(nz),".dat"
!    open(unit=10, file="Data/"//trim(filename), status='replace',access='sequential')
    open(unit=10, file=trim(filename), status='replace',access='sequential')

    write (10,"(A53)")"#---------------------------------------------------#"
    write (10,"(A, F8.5, A, F8.5, A, F8.5, A)") "#--- z =",parz(nz),", D(z) =",Dz,", f(z) =",fz," ---#"
    write (10,"(A301)") "#----------------------------------------------------------------------------&
                         -----------------------------------------------------------------------------&
                         -----------------------------------------------------------------------------&
                         ---------------------------------------------------------------------#"
    write (10,"(A301)") "#      k [h/Mpc]             P_Zel                P_A                 P_W                 &
                   P_d                P_dd                P_d^2              P_d^2d^2             P_dd^2              &
                   P_s^2               P_ds^2             P_d^2s^2            P_s^2s^2             P_D2d              P_dD2d      #"
    write (10,"(A301)") "#----------------------------------------------------------------------------&
                         -----------------------------------------------------------------------------&
                         -----------------------------------------------------------------------------&
                         ---------------------------------------------------------------------#"

!---  compute spectra  ---!

    pp=0 ! starting percent
    write (6,*), pp,"%"

 do p=1,numkk,1

    call percent(p,numkk,10,pp)
    k_val = kkout(p)
!    write (6,*), k_val
 
    P_dm_val = 0.0_dp
    P_bias_df_val = 0.0_dp
    P_bias_s2_val = 0.0_dp
    P_bias_hd_val = 0.0_dp

    call PS_dm(qqin,XXin,YYin,sgin,XXinloop,YYinloop,sginloop,VVin,TTin,k_val,P_dm_val)
    call PS_bias_df(qqin,XXin,YYin,sgin,Xiinlin,U10inlin,U10insum,U20inloop,U11inloop,XX10inloop,YY10inloop,k_val,P_bias_df_val)
    call PS_bias_s2(qqin,XXin,YYin,sgin,Zetainloop,Chiinloop,V10inloop,V12inloop,XX20inloop,YY20inloop,k_val,P_bias_s2_val)
    call PS_bias_hd(qqin,XXin,YYin,sgin,dXiinlin,dU10inlin,k_val,P_bias_hd_val)

!---  output spectra data ---!

    write(10,'(15ES20.10)') k_val, (P_dm_val(i), i=1,3), (P_bias_df_val(i), i=1,5), (P_bias_s2_val(i), i=1,4), (P_bias_hd_val(i), i=1,2)

 enddo ! k-loop
   close(unit=10)
 enddo ! z-loop

!======================!
  stop
end program PS_main



