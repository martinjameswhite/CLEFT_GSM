!===================================================!

program ps02_main

!       Program computes the power spectrum terms in LEFT formalism.
!       Program uses the spherical Hankel transform routine.

  use ps02Module
  implicit none

!==== My variables ====!

  integer, parameter :: numkk=400
  integer :: i,j,l,p,nz,pp
  real(kind=dp) :: dlnr, dlogr,logkc,logrc,logrmax,logrmin,nc,r
  real(kind=dp) :: Dz, fz
  real(kind=dp) :: kkoutMin,kkoutMax,dkk,cfkkout,fixedk
  real(kind=dp) :: k_val
  real(kind=dp), dimension(3,2) :: P02_dm_val
  real(kind=dp), dimension(3,2) :: P02_bias_val
  real(kind=dp), dimension(n,16) :: in
  real(kind=dp), dimension(n,5) :: in2
  real(kind=dp), dimension(n) :: rr
  real(kind=dp), dimension(n) :: qqin,XXin,YYin,X0in,sgin
  real(kind=dp), dimension(n) :: XXinloop,YYinloop,X0inloop,sginloop
  real(kind=dp), dimension(n) :: XX2inloop,YY2inloop,sg2inloop
  real(kind=dp), dimension(n) :: XX22inloop,YY22inloop
  real(kind=dp), dimension(n) :: XX13inloop,YY13inloop
  real(kind=dp), dimension(n) :: VVin, TTin, VV112, TT112, VVintot, TTintot
!  real(kind=dp), dimension(n) :: Xiinlin,U10insum
  real(kind=dp), dimension(n) :: Xiinlin,U10inlin,U10inloop
  real(kind=dp), dimension(n) :: U20inloop,U11inloop,XX10inloop,YY10inloop
  real(kind=dp), dimension(numkk) :: kkout

  real(kind=dp), dimension(9) :: parz = (/0.00,0.25,0.55,0.60,0.65,0.70,0.75,0.80,1.00/)

  character(len=1024) :: filename1,filename2
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

  open(20,file='../Data/Ens_ANL15_FFT_4096_46_lin&loop.dat',status='old')
  read(20,*)
  read(20,*)
  read(20,*)
    do j=1,n
     read(20,*) in(j,:)
    enddo
  close(20)

  open(21,file='../Data/Ens_ANL15_FFT_4096_46_APP_lin&loop_add.dat',status='old')
  read(21,*)
  read(21,*)
  read(21,*)
    do j=1,n
     read(21,*) in2(j,:)
    enddo
  close(21)

!---   redshift loop   ---!
! do nz=8,8,1
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
  sgin = Dz**2.0_dp * in(:,3)/3.0_dp
  X0in = Dz**2.0_dp * 2.0_dp/3.0_dp * (        - in(:,4) - in(:,5))
  XXin = Dz**2.0_dp * 2.0_dp/3.0_dp * (in(:,3) - in(:,4) - in(:,5))
  YYin = Dz**2.0_dp * 2.0_dp * in(:,5)

  sginloop = Dz**4.0_dp * in(:,6)/3.0_dp
  X0inloop = Dz**4.0_dp * 2.0_dp/3.0_dp * (        - in(:,7) - in(:,8))
  XXinloop = Dz**4.0_dp * 2.0_dp/3.0_dp * (in(:,6) - in(:,7) - in(:,8))
  YYinloop = Dz**4.0_dp * 2.0_dp * in(:,8)

  XX22inloop = Dz**4.0_dp * 2.0_dp/3.0_dp * (in2(1,2)*in2(:,1)/in2(:,1) - in2(:,2) - in2(:,3))
  YY22inloop = Dz**4.0_dp * 2.0_dp * in2(:,3)
  XX13inloop = Dz**4.0_dp * 2.0_dp/3.0_dp * (in2(1,4)*in2(:,1)/in2(:,1) - in2(:,4) - in2(:,5))
  YY13inloop = Dz**4.0_dp * 2.0_dp * in2(:,5)

  XX2inloop = XX22inloop + 3.d0/4.d0*XX13inloop
  YY2inloop = YY22inloop + 3.d0/4.d0*YY13inloop
  sg2inloop = 1.d0/2.d0*XX2inloop(n)*qqin/qqin

!  XX2inloop = XX22inloop + XX13inloop
!  YY2inloop = YY22inloop + YY13inloop
!  sg2inloop = 1.d0/2.d0*XX2inloop(n)*qqin/qqin
  
  VVin = Dz**4.0_dp * 1.0_dp/5.0_dp* ( 2.0_dp*in(:,9) - 3.0_dp*in(:,10) )
  TTin = Dz**4.0_dp * 3.0_dp* in(:,10)
  VV112 = 1.d0/3.d0*VVin;
  TT112 = 1.d0/3.d0*TTin;
  VVintot = 2.d0*VVin - VV112
  TTintot = 2.d0*TTin - TT112

  Xiinlin    = Dz**2.0_dp * in(:,2)

  U10inlin  = Dz**2.0_dp * in(:,11)
  U10inloop = Dz**4.0_dp * in(:,12)
!  U10insum  = Dz**2.0_dp * in(:,11) + Dz**4.0_dp * in(:,12)

  U20inloop  = Dz**4.0_dp * in(:,13)
  U11inloop  = Dz**4.0_dp * in(:,14)
  XX10inloop = Dz**4.0_dp * in(:,15)
  YY10inloop = Dz**4.0_dp * in(:,16)



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

    write(*,*) 'Power Spectrum Terms: spin 0 & 2'
!--- spin=0 ---!
    write (filename1,"(A19,A3,A4)") "mathP2_12_0_hh_46_z",zstr(nz),".dat"
    open(unit=10, file="Data/"//trim(filename1), status='replace',access='sequential')

    write (10,*) "#---------------------------------------------------#"
    write (10,"(A, F8.5, A, F8.5, A, F8.5, A)") " #--- z =",parz(nz),", D(z) =",Dz,", f(z) =",fz," ---#"
    write (10,"(A163)") "#---------------------------------------------------------------------------------------&
                   --------------------------------------------------------------------------#"
    write (10,"(A143)") "#     k [h/Mpc]           P01_Zel              P01_A              P01_W                 &
                 P01_F'            P01_(F')^2            P01_F''   #"
    write (10,"(A163)") "#---------------------------------------------------------------------------------------&
                   --------------------------------------------------------------------------#"
!--- spin=2 ---!
    write (filename2,"(A19,A3,A4)") "mathP2_12_2_hh_46_z",zstr(nz),".dat"
    open(unit=11, file="Data/"//trim(filename2), status='replace',access='sequential')

    write (11,*) "#---------------------------------------------------#"
    write (11,"(A, F8.5, A, F8.5, A, F8.5, A)") " #--- z =",parz(nz),", D(z) =",Dz,", f(z) =",fz," ---#"
    write (11,"(A163)") "#---------------------------------------------------------------------------------------&
                   --------------------------------------------------------------------------#"
    write (11,"(A143)") "#     k [h/Mpc]           P01_Zel              P01_A              P01_W                 &
                 P01_F'            P01_(F')^2            P01_F''   #"
    write (11,"(A163)") "#---------------------------------------------------------------------------------------&
                   --------------------------------------------------------------------------#"

!---  compute spectra  ---!

    pp=0 ! starting percent
    write (6,*), pp,"%"

 do p=1,numkk,1

    call percent(p,numkk,10,pp)
    k_val = kkout(p)
!    write (6,*), k_val
 
    P02_dm_val(:,:) =  0.0_dp
    P02_bias_val(:,:)  =  0.0_dp

    call ps02_dm(qqin,XXin,YYin,sgin,XX2inloop,YY2inloop,sg2inloop,VVintot,TTintot,k_val,P02_dm_val)
    call ps02_bias(qqin,XXin,YYin,sgin,Xiinlin,U10inlin,XX10inloop,YY10inloop,k_val,P02_bias_val)

!---  output spectra data ---!

    write(10,'(7ES20.10)') k_val, (P02_dm_val(i,1), i=1,3), (P02_bias_val(i,1), i=1,3)
    write(11,'(7ES20.10)') k_val, (P02_dm_val(i,2), i=1,3), (P02_bias_val(i,2), i=1,3)

 enddo ! k-loop
   close(unit=10)
   close(unit=11)
 enddo ! z-loop

!======================!
  stop
end program ps02_main



