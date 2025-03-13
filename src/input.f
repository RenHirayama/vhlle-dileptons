c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE INPUT
c..
c.. This is the main coarse-graining subroutine
c.. for UrQMD file14 (with multiple timesteps) read-in,
c.. calling the other calculation routines
c*****|****************************************************************|X

      subroutine input(grd,grd_z,dx,vol,vac)
      implicit none

      include 'defs.f'

c. general
      real*8 dx,dt,vol,vol4
      integer grd,grd_z,vac

c. event header, one for every collision

      integer filetype
      integer UrQMDver
c.. number of events
      integer noe,noefo,evts,nouse,noecount
      logical emptyev
c.. impact parameter
      real*8 b(evtmax),bmin,bmax
c.. projectile & target
      integer proA,proZ,tarA,tarZ
c.. cross-section
      real*8 sigma
c.. transformation betas
      real*8 betaNN,betaLAB,betaPRO
c..  beam & c.o.m. energy / beam momentum
      real*8 elab,ecm,plab,pcm
      real*8 snn
c.. time & timestep
      real*8 deltat
      integer tottime

c. functions
      real*8 temp,chem,cs,mupik
      real*8 ranff
      external ranff
      integer ityp
      external ityp
      integer charge
      external charge
c.. temperature & chemical potential
      real*8 t,mub,muq,lambda

c. timestep header, one for every timestep
      integer nexit(evtmax,timemax)
      integer nexitcount

c. event body

c.. read OUTgoing particles
      real*8 OUTr0(timemax,outmax),OUTrx(timemax,outmax)
      real*8 OUTry(timemax,outmax),OUTrz(timemax,outmax)
      real*8 OUTp0(timemax,outmax),OUTpx(timemax,outmax)
      real*8 OUTpy(timemax,outmax),OUTpz(timemax,outmax)
      real*8 OUTmass(timemax,outmax)
      integer OUTityp(timemax,outmax),OUTiso3(timemax,outmax)
      integer OUTch(timemax,outmax),OUTlcoll(timemax,outmax)
      integer OUTcoll(timemax,outmax),OUTpaproc(timemax,outmax)
      integer OUTorigin(timemax,outmax)
      integer tformflag

c. for subroutines
c.. GRID
c... cell properties
      real*8 nopart(timemax,grd,grd,grd_z)

      real*8 norho0(timemax,grd,grd,grd_z)
      real*8 mrho0(timemax,grd,grd,grd_z)
      real*8 p0rho0(timemax,grd,grd,grd_z)
      real*8 pxrho0(timemax,grd,grd,grd_z)
      real*8 pyrho0(timemax,grd,grd,grd_z)
      real*8 pzrho0(timemax,grd,grd,grd_z)

      real*8 nodelta(timemax,grd,grd,grd_z)
      real*8 mdelta(timemax,grd,grd,grd_z)
      real*8 p0delta(timemax,grd,grd,grd_z)
      real*8 pxdelta(timemax,grd,grd,grd_z)
      real*8 pydelta(timemax,grd,grd,grd_z)
      real*8 pzdelta(timemax,grd,grd,grd_z)

      real*8 noomega(timemax,grd,grd,grd_z)
      real*8 momega(timemax,grd,grd,grd_z)
      real*8 p0omega(timemax,grd,grd,grd_z)
      real*8 pxomega(timemax,grd,grd,grd_z)
      real*8 pyomega(timemax,grd,grd,grd_z)
      real*8 pzomega(timemax,grd,grd,grd_z)

      real*8 nophi(timemax,grd,grd,grd_z)
      real*8 mphi(timemax,grd,grd,grd_z)
      real*8 p0phi(timemax,grd,grd,grd_z)
      real*8 pxphi(timemax,grd,grd,grd_z)
      real*8 pyphi(timemax,grd,grd,grd_z)
      real*8 pzphi(timemax,grd,grd,grd_z)

      real*8 p0(timemax,grd,grd,grd_z)
      real*8 px(timemax,grd,grd,grd_z)
      real*8 py(timemax,grd,grd,grd_z)
      real*8 pz(timemax,grd,grd,grd_z)
      real*8 n0(timemax,grd,grd,grd_z)
      real*8 npi(timemax,grd,grd,grd_z)
      real*8 nk(timemax,grd,grd,grd_z)
      real*8 jx(timemax,grd,grd,grd_z)
      real*8 jy(timemax,grd,grd,grd_z)
      real*8 jz(timemax,grd,grd,grd_z)
      real*8 rhoeff(timemax,grd,grd,grd_z)

      real*8 t00(timemax,grd,grd,grd_z)
      real*8 t01(timemax,grd,grd,grd_z)
      real*8 t02(timemax,grd,grd,grd_z)
      real*8 t03(timemax,grd,grd,grd_z)
      real*8 t11(timemax,grd,grd,grd_z)
      real*8 t22(timemax,grd,grd,grd_z)
      real*8 t33(timemax,grd,grd,grd_z)
      real*8 t12(timemax,grd,grd,grd_z)
      real*8 t13(timemax,grd,grd,grd_z)
      real*8 t23(timemax,grd,grd,grd_z)


c.. TRANSFORM
c... local rest-frame properties
      real*8 tau(timemax,grd,grd,grd_z)
      real*8 p0_lrf(timemax,grd,grd,grd_z)
      real*8 n0_lrf(timemax,grd,grd,grd_z)
      real*8 npi_lrf(timemax,grd,grd,grd_z)
      real*8 nk_lrf(timemax,grd,grd,grd_z)
      real*8 jx_lrf(timemax,grd,grd,grd_z)
      real*8 jy_lrf(timemax,grd,grd,grd_z)
      real*8 jz_lrf(timemax,grd,grd,grd_z)
      real*8 edens_lrf(timemax,grd,grd,grd_z)
      real*8 rhoe_lrf(timemax,grd,grd,grd_z)
      real*8 gam(timemax,grd,grd,grd_z)
      real*8 volce(timemax,grd,grd,grd_z)
c.. cell properties
      real*8 edens_cell,n0_cell,vxce,vyce,vzce,gce,volumce,rhonuc
      real*8 pidens_cell,kdens_cell,mupion,mukaon,maxmuk,maxmupi
      real*8 temp_cell(timemax,grd,grd,grd_z)
      real*8 mub_cell(timemax,grd,grd,grd_z)
      real*8 mupi_cell(timemax,grd,grd,grd_z)
      real*8 muk_cell(timemax,grd,grd,grd_z)

      real*8 qgp_vol(timemax)
      real*8 qgp_4vl(timemax)
      real*8 qgp_tot4vl
      real*8 qgp_avtemp(timemax)
      real*8 qgp_tottemp
      real*8 bar_vol(timemax)
      real*8 bar_totnum(timemax)

      real*8 mupi_av(timemax)
      real*8 muk_av(timemax)
      integer cells_t_gt05(timemax)

      integer collnoel(timemax)
      integer collnoin(timemax)

      integer elcount(timemax,grd,grd,grd_z)
      integer inelcount(timemax,grd,grd,grd_z)

      integer fo_cell(timemax,grd,grd,grd_z)
      integer fo_marker(grd,grd,grd_z)
c... counter for filled / empty cells
      integer cells_bar,cells_mes_nobar,cells_lt_minpart

c.. other
      integer e,f,g,h,i,ii,iii,iiii,j,jj,jjj,jjjj,k,z,zz,zzz,zzzz,m
      integer n,o,p,q,r,rr,rrr,w,ww,www,wwww,a1,a2,a3,a4
      integer timesteps, timestepcount, multi, mark1, mark2, freestr
      real*8 bev(evtmax)
      real*8 lack, percent, coarsefac, crse
      integer lack_, counter, counterx, count_fo, evt_modulo
      character dummy, dot
      real*8 dNchdeta(evtmax),dNchdetatot,dNchmin,dNchmax,part_y
      real*8 avtempfo
      real*8 avcoll
      real*8 plateau
      logical delemit

      logical lhcinput

      real*8 nucmass
      PARAMETER(nucmass=0.938d0)

      integer urqmdityp,urqmdcharge

      character*32 betalab_,ecm_
      character*32 control

c.. set all to initial values
      cells_bar=0
      cells_mes_nobar=0
      cells_lt_minpart=0

      qgp_tot4vl=0.0d0
      qgp_tottemp=0.0d0
      avtempfo=0.0d0
      noe=0
      noefo=0
      counter=0
      counterx=1
      count_fo=0
      dot='.'
      multi=1

      maxmuk=0.0d0
      maxmupi=0.0d0

      dNchdetatot=0.0d0
      dNchmin=1.0d3
      dNchmax=0.0d0
      nouse=0
      do z=1,evtmax
        dNchdeta(z)=0.0d0
      end do

c.. initialization of cell variables
      write(0,*)'Initialization'
c      do e=1,timemax
c       write(0,'(a1)',advance='NO'),dot
c       write(0,'(a1)'),dot
c       counter=counter+1
c       if(counter.eq.10) then
c        percent=100d0*counterx*counter/timemax
c        counterx=counterx+1
c        counter=0
c        write(0,'(f4.0)',advance='NO'),percent
c        write(0,'(f4.0)'),percent
c        write(0,*),'%'
c       end if
       qgp_vol(1:timemax)=0.0d0
       qgp_4vl(1:timemax)=0.0d0
       qgp_avtemp(1:timemax)=0.0d0
       bar_vol(1:timemax)=0.0d0
       bar_totnum(1:timemax)=0.0d0
       mupi_av(1:timemax)=0.0d0
       muk_av(1:timemax)=0.0d0
       cells_t_gt05(1:timemax)=0
       collnoel(1:timemax)=0
       collnoin(1:timemax)=0
c       do f=1,grd
c        do g=1,grd
c         do h=1,grd_z
          nopart(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          norho0(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          noomega(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          nophi(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          nodelta(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          mrho0(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          momega(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          mphi(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          mdelta(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0

          p0(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          px(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          py(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          pz(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          n0(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          npi(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          nk(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          jx(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          jy(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          jz(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          rhoeff(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          p0_lrf(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          n0_lrf(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          npi_lrf(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          nk_lrf(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          jx_lrf(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          jy_lrf(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          jz_lrf(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          edens_lrf(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          rhoe_lrf(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          volce(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          temp_cell(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          mub_cell(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          mupi_cell(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          muk_cell(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          gam(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          tau(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0

          t01(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          t02(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          t03(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0

          t00(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          t11(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          t22(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          t33(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0

          t12(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          t13(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0
          t23(1:timemax,1:grd,1:grd,1:grd_z)=0.0d0

          fo_cell(1:timemax,1:grd,1:grd,1:grd_z)=0
          fo_marker(1:grd,1:grd,1:grd_z)=0

          elcount(1:timemax,1:grd,1:grd,1:grd_z)=0
          inelcount(1:timemax,1:grd,1:grd,1:grd_z)=0

c         end do
c        end do
c       end do
c      end do
      write(0,*)'completed'
      write(0,*)'****************************************************'
c------------------------------------------------------------------------
c. ****DECIDE WHICH FORMAT TO READ IN (UrQMD or SMASH)****

       goto(991,992)model

c------------------------------------------------------------------------
c. **** START READ-IN UrQMD FILE 14****
 991  continue
c. ---- start loop over EVENTS ----
c      write(0,*)'START READING INPUT FILE 14...'

      do 112 m=1,evtmax
      noe=noe+1
      emptyev=.TRUE.

      if (m.eq.1) then
c. read collision/event header (first event)
      read(*,131,end=99) dummy,lack_,lack_,UrQMDver,dummy,filetype
      read(*,132,end=99) dummy,dummy,proA,proZ,dummy,
     &                   dummy,tarA,tarZ,dummy
      read(*,133,end=99) dummy,betaNN,betaLAB,betaPRO
      read(*,134,end=99) dummy,b(noe),bmin,bmax,dummy,sigma
      read(*,135,end=99) dummy,lack_,dummy,elab,dummy,ecm,dummy,plab
      read(*,136,end=99) dummy,noe,dummy,lack_,dummy,dummy,tottime,
     &                   dummy,deltat

      write(0,'(A20,2x,I5)')'UrQMD version used:',UrQMDver

c. Determine whether calculation is for LHC energies or not. For LHC
c. the f14 file format is different.

      if(ecm.gt.1000.0d0) then
       lhcinput=.TRUE.
      else
       lhcinput=.FALSE.
      endif

      if(lhcinput) write(0,*)'Calculation for LHC energies!'

c. For high energies the precision of transformation betas given in
c. file14 is not sufficient. In this case calculate more precise values
c. manually.

      if(ecm.gt.50.0d0) then
       snn=ecm**2
       pcm=sqrt((snn-(2*nucmass)**2)*(snn)/(4.0*snn))
       betaLAB=pcm/sqrt(pcm**2+nucmass**2)
      endif

c. The determination of UrQMDver, the label for the version of UrQMD
c. which is used for the calculation, is important as the format of
c. the file14 output format differs slightly.

      if(UrQMDver.lt.20030) then
       stop 'Use UrQMD version 2.3 or newer!!!'
      else if(UrQMDver.ge.20030.AND.UrQMDver.lt.30400) then
       read(*,*,end=99)
       read(*,*,end=99)
       read(*,*,end=99)
       read(*,*,end=99)
       read(*,*,end=99)
       read(*,*,end=99)
       read(*,*,end=99)
       read(*,*,end=99)
      else if(UrQMDver.ge.30400) then
       read(*,*,end=99)
       read(*,*,end=99)
       read(*,*,end=99)
       read(*,*,end=99)
       read(*,*,end=99)
       read(*,*,end=99)
       read(*,*,end=99)
       read(*,*,end=99)
       read(*,*,end=99)
       read(*,*,end=99)
       read(*,*,end=99)
      end if

      evt_modulo=mod(noe,10)
c      write(0,*)'modulo',evt_modulo !Debug only
      if(evt_modulo.eq.0) write(0,*)'event no.',noe

      write(0,'(A26,F9.2,2x,F10.2)')'-->Input for ECM / E_lab:',ecm,elab
      write(0,'(A25,F12.9,2x,F12.9)')'-->beta(NN) / beta(LAB):',betaNN,betaLAB
      write(0,'(A21,I5,2x,I5)')'-->Projectile A / Z:',proA,proZ
      write(0,'(A17,I5,2x,I5)')'-->Target A / Z:',tarA,tarZ
      write(0,'(A35,F6.2,2x,F6.2)')'-->Impact Parameter b_min / b_max:',bmin,bmax
      write(0,'(A27,I5,2x,F6.3)')'-->Timesteps / Stepleghth:',int(tottime/deltat+0.01d0),deltat
      write(0,*)'****************************************************'
      if(outputs.ne.1) write(0,*)'***OUTPUT SUPPRESSED!***'
      write(0,*)'START READING INPUT FILE...'

c.. Write header line for standard output
      write(*,761)deltat,int(tottime/deltat+0.01d0),dx,proA,proZ,
     &         ecm,betaLAB,bmin,bmax

c.. Check number of timesteps < timemax
      timesteps=int(tottime/deltat+0.01d0)
      dt=deltat
      if(int(tottime/deltat+0.01d0).gt.timemax) then
        write(0,*)'TOO MANY TIMESTEPS -> INCREASE TIMEMAX!'
        stop
      endif
c      write(0,*)'TIMESTEPS, DT:',timesteps,dt !debug only

      end if ! reading of first event header

      if (m.gt.1) then
c. read collision/event header (first event)
      do z=1,3
      read(*,*,end=99)
      end do
      read(*,134,end=99) dummy,b(noe),bmin,bmax,dummy,sigma
      read(*,135,end=99) dummy,lack_,dummy,elab,dummy,ecm,dummy,plab
      read(*,136,end=99) dummy,noe,dummy,lack_,dummy,dummy,tottime,
     &                   dummy,deltat

      if(UrQMDver.ge.20030.AND.UrQMDver.lt.30400) then
       do zz=1,8
        read(*,*,end=99)
       end do
      else if(UrQMDver.ge.30400) then
         do zz=1,11
        read(*,*,end=99)
       end do
      end if

c.. Number of timesteps
      timesteps=int(tottime/deltat+0.01d0)
      dt=deltat

      evt_modulo=mod(noe,10)
c      write(0,*)'modulo',evt_modulo !Debug only
      if(evt_modulo.eq.0) write(0,*)'event no.',noe

      end if ! reading of event header (noe.ne.1)

c ---- start loop over TIMESTEPS ----
c      write(0,*)'loop over timesteps' !debug only
c.now go through for every timestep

      do 111  i=1,timesteps
c      write(0,*)'timestep',i !debug only
       read(*,301,end=99) nexit(noe,i),lack_
       read(*,302,end=99) lack_,lack_,lack_,lack_,lack_,lack_,
     &                    lack_,lack_

c. check for correct filetype
c..STANDARD OUTPUT FILE 14
      if (filetype.eq.14) then

c. read outgoing particles
c       write(0,*)'nexit',nexit(noe,i) !debug only
       do 21 ii=1,nexit(noe,i)
c       write(0,*)'particle',ii !debug only
        if (ii.ge.outmax) stop 'increase outmax!'
        if(lhcinput) then
         read(*,244,end=99) OUTr0(i,ii),OUTrx(i,ii),
     &    OUTry(i,ii),OUTrz(i,ii),
     &    OUTp0(i,ii),OUTpx(i,ii),OUTpy(i,ii),OUTpz(i,ii),
     &    OUTmass(i,ii),OUTityp(i,ii),OUTiso3(i,ii),
     &    OUTch(i,ii),OUTlcoll(i,ii),OUTcoll(i,ii),
     &    OUTpaproc(i,ii)
        else
         read(*,214,end=99) OUTr0(i,ii),OUTrx(i,ii),
     &    OUTry(i,ii),OUTrz(i,ii),
     &    OUTp0(i,ii),OUTpx(i,ii),OUTpy(i,ii),OUTpz(i,ii),
     &    OUTmass(i,ii),OUTityp(i,ii),OUTiso3(i,ii),
     &    OUTch(i,ii),OUTlcoll(i,ii),OUTcoll(i,ii),
     &    OUTpaproc(i,ii)
        endif
c        write(0,*)'part,OUTp0,OUTpx,OUTpy,OUTpz',ii,OUTp0(i,ii),
c     &             OUTpx(i,ii),OUTpy(i,ii),OUTpz(i,ii),
c     &             OUTrx(i,ii),OUTry(i,ii),OUTrz(i,ii)!Debug only
c.. Check, whether it is an empty event
        if(OUTpaproc(i,ii).ne.0) emptyev=.FALSE.
c.. Count charged particles for dN_ch/deta (only for NA60)
        if(na60mode.ge.1) then
          part_y=0.5d0*log((OUTp0(i,ii)+OUTpz(i,ii))/
     &                     (OUTp0(i,ii)-OUTpz(i,ii)))
          if(abs(part_y).lt.0.5d0.AND.OUTch(i,ii).ne.0
     &       .AND.i.eq.timesteps) then
             dNchdetatot=dNchdetatot+1.0d0
             dNchdeta(m)=dNchdeta(m)+1.0d0
          endif
        endif
 21     continue

c..EXPERIMENTAL!:
c..MODIFIED FILE 44 FOR HYDRO+COARSE USE
c..In this output file an additional flag is used to indicate whether
c..a particle from hydro is already produced (t>tform) or does not yet
c..exist
      else if (filetype.eq.44) then

c. read outgoing particles
c       write(0,*)'nexit',nexit(noe,i) !debug only
       do 22 ii=1,nexit(noe,i)
c       write(0,*)'particle',ii !debug only
        if (ii.ge.outmax) stop 'increase outmax!'
        if(lhcinput) then
         read(*,2441,end=99) tformflag,OUTr0(i,ii),
     &    OUTrx(i,ii),OUTry(i,ii),OUTrz(i,ii),
     &    OUTp0(i,ii),OUTpx(i,ii),OUTpy(i,ii),OUTpz(i,ii),
     &    OUTmass(i,ii),OUTityp(i,ii),OUTiso3(i,ii),
     &    OUTch(i,ii),OUTlcoll(i,ii),OUTcoll(i,ii),
     &    OUTpaproc(i,ii)
        else
         read(*,2141,end=99) tformflag,OUTr0(i,ii),
     &    OUTrx(i,ii),OUTry(i,ii),OUTrz(i,ii),
     &    OUTp0(i,ii),OUTpx(i,ii),OUTpy(i,ii),OUTpz(i,ii),
     &    OUTmass(i,ii),OUTityp(i,ii),OUTiso3(i,ii),
     &    OUTch(i,ii),OUTlcoll(i,ii),OUTcoll(i,ii),
     &    OUTpaproc(i,ii)
        endif
c        write(0,*)'part,OUTp0,OUTpx,OUTpy,OUTpz',ii,OUTp0(i,ii),
c     &             OUTpx(i,ii),OUTpy(i,ii),OUTpz(i,ii),
c     &             OUTrx(i,ii),OUTry(i,ii),OUTrz(i,ii)!Debug only
c        write(0,*)'tformflag',tformflag !Debug only
        if(tformflag.eq.0.AND.OUTpaproc(i,ii).eq.96) OUTityp(i,ii)=0
c.. Check, whether it is an empty event
        if(emptyev.AND.OUTpaproc(i,ii).ne.0) then
           emptyev=.FALSE.
c           write(0,*)'noe,paproc,emptyev',noe,OUTpaproc(i,ii),emptyev !Debug only
        endif
c        write(0,*)'part,OUTp0,OUTpx,OUTpy,OUTpz',ii,OUTp0(i,ii),
c     &             OUTpx(i,ii),OUTpy(i,ii),OUTpz(i,ii),
c     &             OUTrx(i,ii),OUTry(i,ii),OUTrz(i,ii)!Debug only
c.. Count charged particles for dN_ch/deta (only for NA60)
        if(na60mode.ge.1) then
          part_y=0.5d0*log((OUTp0(i,ii)+OUTpz(i,ii))/
     &                     (OUTp0(i,ii)-OUTpz(i,ii)))
          if(abs(part_y).lt.0.5d0.AND.OUTch(i,ii).ne.0
     &       .AND.i.eq.timesteps) then
             dNchdetatot=dNchdetatot+1.0d0
             dNchdeta(m)=dNchdeta(m)+1.0d0
          endif
        endif
 22   continue

      else
       write(0,*)'wrong input file format'
       stop
      end if

 111  continue
c. now put particles into the grid
c       if(timesteps.gt.timemax) timesteps=timemax
c       write(0,*)'EVENT, DNchDETA(EVENT)',m,noe,dNchdeta(m) !Debug only
       if((dNchdeta(m).ge.30.0d0.OR.na60mode.ne.1).AND.(.NOT.emptyev)) then
c       write(0,*)'CALL GRID' !Debug only
       call grid(nexit,noe,timesteps,grd,grd_z,dx,dt,vol,OUTr0,
     &            OUTrx,OUTry,OUTrz,OUTp0,OUTpx,OUTpy,OUTpz,OUTmass,
     &            OUTityp,OUTiso3,OUTch,OUTlcoll,OUTcoll,OUTpaproc,
     &            OUTorigin,nopart,norho0,mrho0,p0rho0,pxrho0,pyrho0,
     &            pzrho0,nodelta,mdelta,p0delta,pxdelta,pydelta,pzdelta,
     &            noomega,momega,p0omega,pxomega,pyomega,pzomega,
     &            nophi,mphi,p0phi,pxphi,pyphi,pzphi,
     &            p0,px,py,pz,n0,npi,nk,jx,jy,jz,t00,t01,t02,t03,t11,
     &            t22,t33,t12,t13,t23,rhoeff)
       else if(dNchdeta(m).lt.30.0d0.AND.na60mode.eq.1.AND.(.NOT.emptyev)) then
        nouse=nouse+1
        dNchdetatot=dNchdetatot-dNchdeta(m)
       else if(emptyev) then
        nouse=nouse+1
       end if

c. New Minimum / Maximum for dN_ch/deta? (for NA60 cut only)
      if(na60mode.ge.1) then
        if(dNchdeta(m).gt.dNchmax) dNchmax=dNchdeta(m)
        if(dNchdeta(m).lt.dNchmin) dNchmin=dNchdeta(m)
      end if
c. start all over with next event
      if ((noe.ge.evtmax).OR.(noe.ge.100.AND.focalc)) then
c     if (noe.ge.100) then !FOR TESTPURPOSE ONLY!
        evts=noe
       goto 89
      endif
 112  continue


c. **** START READ-IN SMASH OSCAR FILE****
 992  continue
c. ---- start loop over EVENTS ----
      write(0,*)'START READING SMASH OSCAR FILE...'

      call getenv('ecm',ecm_)
      read(ecm_,*)  ecm

      snn=ecm**2
      pcm=sqrt((snn-(2*nucmass)**2)*(snn)/(4.0*snn))
      betaLAB=pcm/sqrt(pcm**2+nucmass**2)

      emptyev=.FALSE.

c.. Read file header
c     write(0,*)'Read file header...' !Debug only
      read(*,*,end=99)
      read(*,*,end=99)
      read(*,*,end=99)

      noecount=0
      timestepcount=0
      nexitcount=0

 9112 continue
c      write(0,*)'Read timestep / event header' !Debug only
      read(*,*,end=98) dummy,dummy,noe,dummy,dummy,control,nexitcount
c      write(0,*)'noe,nexitcount,control',noe,nexitcount,control !Debug only

c.. **OLD EVENT** finished - Put particles from last event on the grid
      if(control.eq.'end') then
        call grid(nexit,noecount,timesteps,grd,grd_z,dx,dt,vol,
     &   OUTr0,OUTrx,OUTry,OUTrz,OUTp0,OUTpx,OUTpy,OUTpz,OUTmass,
     &   OUTityp,OUTiso3,OUTch,OUTlcoll,OUTcoll,OUTpaproc,
     &   OUTorigin,nopart,norho0,mrho0,p0rho0,pxrho0,pyrho0,
     &   pzrho0,nodelta,mdelta,p0delta,pxdelta,pydelta,pzdelta,
     &   noomega,momega,p0omega,pxomega,pyomega,pzomega,
     &   nophi,mphi,p0phi,pxphi,pyphi,pzphi,
     &   p0,px,py,pz,n0,npi,nk,jx,jy,jz,t00,t01,t02,t03,t11,
     &   t22,t33,t12,t13,t23,rhoeff)
        goto 9112
      endif
c.. **NEW EVENT** starts
      if(noe.gt.noecount) then
c.. After Event 1, define timesteps and duration
        if(noecount.eq.1) then
          timesteps=timestepcount
          tottime=OUTr0(i,1)-OUTr0(1,1)
          deltat=OUTr0(2,1)-OUTr0(1,1)
          write(0,*)OUTr0(2,1),OUTr0(1,1)
          dt=deltat
          write(0,*)'****************************************************'
          write(0,'(A27,F9.2)')'--> Center-of-mass energy:',ecm
          write(0,'(A27,F9.6)')'-->   Transformation beta:',betaLAB
c          write(0,'(A27,I5,2x,F6.3)')'-->Timesteps / Stepleghth:',timesteps,dt
          write(0,*)'-->Timesteps / Stepleghth:',timesteps,dt
          write(0,*)'****************************************************'
        endif
c..Start with new event, except if mamimum number of events is reached
c... Event counter
        noecount=noe
c... Reset timestep counter
        i=1
        if(noecount.eq.1) timestepcount=1
        if(noecount.gt.evtmax) then
          noecount=evtmax
          noe=evtmax
          goto 89
        endif
c... Print out event number
        evt_modulo=mod(noe,10)
        if(evt_modulo.eq.0) write(0,*)'event no.',noe
c.. A **NEW TIMESTEP** starts
      else !(.NOT.(noe.gt.noecount))
        i=i+1
        if(noecount.eq.1) timestepcount=timestepcount+1
c.. Stop routine if  mamimum number of timesteps is reached
        if(i.gt.timemax) then
          write(0,*)'TOO MANY TIMESTEPS -> INCREASE TIMEMAX!'
          stop
        endif
c... Ignore last timestep if negative
        if(OUTr0(i-1,1).le.0.0000000001d0) then
          i=1
          if(noecount.eq.1) timestepcount=1
        endif
      endif

      nexit(noe,i)=nexitcount

c..Now read in the particle output
      do ii=1,nexit(noe,i)
         read(*,*,end=99) OUTr0(i,ii),OUTrx(i,ii),
     &    OUTry(i,ii),OUTrz(i,ii),OUTmass(i,ii),
     &    OUTp0(i,ii),OUTpx(i,ii),OUTpy(i,ii),OUTpz(i,ii),
     &    OUTityp(i,ii),OUTpaproc(i,ii)
c         write(0,*) OUTr0(i,ii),OUTrx(i,ii),
c     &    OUTry(i,ii),OUTrz(i,ii),OUTmass(i,ii),
c     &    OUTp0(i,ii),OUTpx(i,ii),OUTpy(i,ii),OUTpz(i,ii),
c     &    OUTityp(i,ii),OUTpaproc(i,ii) !Debug only
c...NOTE: Particle ID in SMASH-OSCAR file is PDG-ID. However,
c...      for further calculation the program uses the UrQMD-ID.
          urqmdcharge=charge(OUTityp(i,ii))
          OUTch(i,ii)=urqmdcharge
          urqmdityp=ityp(OUTityp(i,ii))
          OUTityp(i,ii)=urqmdityp
c          write(0,*)'ITYP,CHARGE',OUTityp(i,ii),OUTch(i,ii)
      enddo

      goto 9112

c ---- all read in ----
 98   continue
      call grid(nexit,noecount,timesteps,grd,grd_z,dx,dt,vol,
     & OUTr0,OUTrx,OUTry,OUTrz,OUTp0,OUTpx,OUTpy,OUTpz,OUTmass,
     & OUTityp,OUTiso3,OUTch,OUTlcoll,OUTcoll,OUTpaproc,
     & OUTorigin,nopart,norho0,mrho0,p0rho0,pxrho0,pyrho0,
     & pzrho0,nodelta,mdelta,p0delta,pxdelta,pydelta,pzdelta,
     & noomega,momega,p0omega,pxomega,pyomega,pzomega,
     & nophi,mphi,p0phi,pxphi,pyphi,pzphi,
     & p0,px,py,pz,n0,npi,nk,jx,jy,jz,t00,t01,t02,t03,t11,
     & t22,t33,t12,t13,t23,rhoeff)
      evts=noe
      goto 89
 99   continue
      evts=noe-1
 89   continue
      if(timesteps.gt.timemax) timesteps=timemax
c...SUBTRACT EXTRACTED EVENTS FROM TOTAL NUMBER (FOR NA60 CUT ONLY)
      if(na60mode.eq.1) then
        evts=evts-nouse
      endif
c.. Write out dN_ch/d_eta for |eta|<0.5
      if(na60mode.ge.1.AND.model.eq.1) then
        write(0,*)'****************************************************'
        write(0,*)'dN_ch/d_eta for |eta|<0.5 =',dNchdetatot/evts
        write(0,'(A25,F7.2,F7.2)')'Minimum & Maximum Value:',dNchmin,dNchmax
        write(0,*)'****************************************************'
      end if !na60mode
      write(0,*)'READ-IN COMPLETED'
      write(0,*)'****************************************************'
c------------------------------------------------------------------------
c..EXPERIMENTAL (!) PART FOR FREEZE-OUT CALCULATIONS. NEEDS FURTHER WORK, NOT
c..TO BE USED FOR ANY PHYSICAL CALCULATIONS IN THE PRESENT FORM!
c..Activate by setting focalc=.true. in main program coarse

c. First,check whether FREEZE-OUT CELLS shall be determined
      if(focalc) then
c. **** CALL "FREEZEOUT" TO READ-IN FILE 15****
       write(0,*)'Now read in f15 to determine FREEZE-OUT CELLS...'
       CALL freezeout(grd,grd_z,dx,dt,vol,timemax,fo_cell,collnoel,
     &                collnoin,elcount,inelcount,noefo)
       write(0,*)'READ-IN OF FILE 15 COMPLETED'
       write(0,*)'****************************************************'
       write(74,*)'#ev_urqmd, delta_t',noe,dt
       write(74,599)'#timestep','coll_in,coll_el'
c. Global counter of collisions per timestep
       do w=1,timemax
       write(74,511)w,dble(collnoel(w))/dble(noefo),dble(collnoin(w))/dble(noefo)
       enddo
c. Determine the freeze-out time for each cell
       do a2=1,grd
        do a3=1,grd
         do a4=1,grd_z
          avcoll=0
          plateau=0
          mark1=0
          freestr=0
          write(0,*)'space-cell (x,y,z):',a2,a3,a4 ! Debug only
          do a1=1,timemax
           write(0,*)'timestep, N_coll',a1,inelcount(a1,a2,a3,a4) ! Debug only
           if(a1.gt.7) then
            avcoll=(dble(inelcount(a1-2,a2,a3,a4))+dble(inelcount(a1-3,a2,a3,a4))
     &             +dble(inelcount(a1-4,a2,a3,a4))+dble(inelcount(a1-5,a2,a3,a4))
     &             +dble(inelcount(a1-6,a2,a3,a4))+dble(inelcount(a1-7,a2,a3,a4)))/6.0d0
           else
           avcoll=0.0d0
           endif
           if(fo_marker(a2,a3,a4).gt.-6) then
            if(dble(inelcount(a1,a2,a3,a4)).gt.0.85d0*avcoll.AND.dble(inelcount(a1,a2,a3,a4))
     &         .lt.1.15d0*avcoll.AND.dble(inelcount(a1,a2,a3,a4)).gt.0) then
             fo_marker(a2,a3,a4)=fo_marker(a2,a3,a4)-1
             plateau=plateau+dble(inelcount(a1,a2,a3,a4))/6.0d0
            endif
            if(dble(inelcount(a1,a2,a3,a4)).lt.0.85d0*avcoll.AND.dble(inelcount(a1,a2,a3,a4))
     &         .gt.1.15d0*avcoll.AND.dble(inelcount(a1,a2,a3,a4)).gt.0) then
             fo_marker(a2,a3,a4)=0
             plateau=0.0d0
            endif
           endif
           if(fo_marker(a2,a3,a4).eq.-6.AND.inelcount(a1,a2,a3,a4).lt.0.8d0*plateau.AND.freestr.eq.0) then
             mark1=a1
             freestr=1
           endif
           if(fo_marker(a2,a3,a4).eq.-6.AND.inelcount(a1,a2,a3,a4).gt.0.9d0*plateau.AND.freestr.eq.1) then
             mark1=0
             freestr=0
           endif
           if(fo_marker(a2,a3,a4).eq.-6.AND.inelcount(a1,a2,a3,a4).lt.0.5d0*plateau) then
             fo_marker(a2,a3,a4)=mark1
             write(0,*)'Chemical freezeout at timestep',fo_marker(a2,a3,a4) ! Debug only
           endif
          enddo !a1
         enddo !a4
        enddo !a3
       enddo !a2
      endif !focalc
c------------------------------------------------------------------------
c. *** START ACTUAL CALCULATIONS ****
      write(0,*)' CONTINUE WITH TRANSFORMATIONS'
c.. LORENTZ TRANSFORMATION of lab-frame quantities
c.. into rest-frame quantities

      call transform(evts,grd,grd_z,dx,timesteps,dt,nopart,p0,n0,npi,nk,
     &              jx,jy,jz,px,py,pz,tau,p0_lrf,n0_lrf,npi_lrf,nk_lrf,
     &              jx_lrf,jy_lrf,jz_lrf,cells_bar,cells_mes_nobar,gam,
     &              vol,volce,t00,t01,t02,t03,t11,t22,t33,t12,t13,t23,
     &              edens_lrf,rhoeff,rhoe_lrf,cells_lt_minpart)

      write(0,*)'Events, not-used',evts,nouse
      write(0,*)'Cells with baryon content = ',cells_bar
      write(0,*)'Cells with meson and without baryon content = ',
     &           cells_mes_nobar
      if(rates.eq.-1) then
       write(0,*)'Cells with less than 9 part. / 7 bar.(ignored) = ',
     &           cells_lt_minpart
      end if

c.. do calculations using the transformed quantites
      if(rates.eq.-1) then
         write(0,*)'START PHOTON EMISSION...'
      else
         write(0,*)'START DILEPTON EMISSION...'
      end if
c... loop over cells
      counter=0
      counterx=0
      do 804 h=1,timesteps
c      do 804 h=1,20
c      do 804 h=21,timesteps

c      call timestepdil(h,dt)

       do 803 i=1,grd

         write(0,'(a1)',advance='NO'),dot
c         write(0,'(a1)'),dot

        do 802 j=1,grd
         do 801 k=1,grd_z

         counter=counter+1
         counterx=counterx+1
         if(counter.eq.100000) then
          percent=100d0*counterx/(timesteps*grd*grd*grd_z)
          counter=0
          write(0,'(f7.2)',advance='NO'),percent
c          write(0,'(f4.0)'),percent
          write(0,*),'%'
         end if


c... get cell properties from arrays

          n0_cell=n0_lrf(h,i,j,k)
          edens_cell=edens_lrf(h,i,j,k)
          pidens_cell=npi_lrf(h,i,j,k)
          kdens_cell=nk_lrf(h,i,j,k)

          if(edens_cell.eq.0.0d0.OR.n0_cell.eq.0.0d0) cycle

c          write(0,*),h,i,j,k !debug only

          vxce=jx(h,i,j,k)/n0(h,i,j,k)
          vyce=jy(h,i,j,k)/n0(h,i,j,k)
          vzce=jz(h,i,j,k)/n0(h,i,j,k)
          gce=gam(h,i,j,k)
          volumce=(volce(h,i,j,k)) ! in fm**[3]
          vol4=tau(h,i,j,k)*volumce

c          write(0,*)'vol4=',vol4 !Debug only
c          write(0,*)'n0_cell',n0_cell !debug only
c          write(0,*)'edens_cell',edens_cell !debug only

c... Conversion of energy and baryon density in units of n_0 and
c... e_0 (necessary for funcitions temp and chem)

          n0_cell=n0_cell/n_0
          edens_cell=(edens_cell*1.0d3)/e_0

c          if(edens_cell.lt.4.0d-4) cycle

c          write(0,*)'n0_cell',n0_cell !debug only
c          write(0,*)'edens_cell',edens_cell !debug only

c... determine TEMPERATURE of cell

c          write(0,*)'Determine TEMPERATURE of cell' !debug only

          t=temp(edens_cell,n0_cell)


c          write(0,*)'temp',t !Debug only

c... determine BARYON CHEMICAL POTENTIAL of cell

c          write(0,*)'Determine CHEMICAL POTENTIAL of cell' !debug only

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!
c..        ATTENTION!!! Note the difference between mu_q and mu_b!!! !
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!

           muq=chem(edens_cell,n0_cell)
           mub=3.0d0*muq

           if(baryons.eq.0) then
            mub=0.0d0
            muq=0.0d0
           endif

c          write(0,*)'mub',mub ! Debug only


c... determine FRACTION OF QGP in cell

           lambda=0.0d0
           if(eos.eq.3.OR.eos.eq.5) then
             lambda=cs(edens_cell,n0_cell)
             if(eos.eq.5.AND.t.gt.175.0d0) lambda=1.0d0
             if(eos.eq.3.AND.t.gt.200.0d0) lambda=1.0d0
             if(eos.eq.3.AND.t.lt.100.0d0) lambda=0.0d0
             if(lambda.gt.0.0d0) then
              qgp_vol(h)=qgp_vol(h)+lambda*volumce
              qgp_4vl(h)=qgp_4vl(h)+lambda*vol4
              qgp_tot4vl=qgp_tot4vl+lambda*vol4
              qgp_avtemp(h)=qgp_avtemp(h)+t*lambda*vol4
              qgp_tottemp=qgp_tottemp+t*lambda*vol4
c              write(0,*)qgp_vol(h),qgp_4vl(h) !Debug only
c              write(0,*)lambda*volumce,lambda*vol4 !Debug only
             endif
           endif
           if(eos.eq.4) then
             if(t.gt.170.0d0) lambda=1.0d0
              if(lambda.gt.0.0d0) then
              qgp_vol(h)=qgp_vol(h)+lambda*volumce
              qgp_4vl(h)=qgp_4vl(h)+lambda*vol4
              qgp_tot4vl=qgp_tot4vl+lambda*vol4
              qgp_avtemp(h)=qgp_avtemp(h)+t*lambda*vol4
              qgp_tottemp=qgp_tottemp+t*lambda*vol4
c              write(0,*)qgp_vol(h),qgp_4vl(h) !Debug only
c              write(0,*)lambda*volumce,lambda*vol4 !Debug only
             endif
           endif
           if(eos.eq.6) then
             if(t.gt.170.0d0) lambda=1.0d0
              if(lambda.gt.0.0d0) then
              qgp_vol(h)=qgp_vol(h)+lambda*volumce
              qgp_4vl(h)=qgp_4vl(h)+lambda*vol4
              qgp_tot4vl=qgp_tot4vl+lambda*vol4
              qgp_avtemp(h)=qgp_avtemp(h)+t*lambda*vol4
              qgp_tottemp=qgp_tottemp+t*lambda*vol4
c              write(0,*)qgp_vol(h),qgp_4vl(h) !Debug only
c              write(0,*)lambda*volumce,lambda*vol4 !Debug only
             endif
            endif

c... Conversion back to *GeV*

          t=abs(t)*1.0d-3
          mub=abs(mub)*1.0d-3
          muq=abs(muq)*1.0d-3

c... determine RHOEFF / PION and KAON CHEMICAL POTENTIAL of cell
           mupion=0.0d0
           mukaon=0.0d0

          if(pikchem.eq.1) then
           mupion=mupik(t,edens_cell,pidens_cell,n0_cell,lambda,0)
           if(rates.eq.-1.OR.rates.eq.1.OR.rates.eq.3) then
            mukaon=mupik(t,edens_cell,kdens_cell,n0_cell,lambda,1)
           endif

           if(.NOT.(mupion.ge.-1.0d0)) mupion=-1.0d0
           if(.NOT.(mukaon.ge.-1.0d0)) mukaon=-1.0d0

           if(mupion.gt.0.141d0) then
c            write(0,*)'large mu_pion:',mupion ! Debug only
            mupion=0.1409d0
           endif
           if(mukaon.gt.0.450d0) then
c            write(0,*)'large mu_kaon',mukaon ! Debug only
            mukaon=0.4499d0
           endif
c          write(0,*)'MUPI in cell:',h,i,j,k,mupion ! Debug only
          endif !(pikchem.eq.1)

          if(t.gt.0.1d0.AND.t.lt.0.175d0) then
             maxmuk=max(maxmuk,mukaon)
             maxmupi=max(maxmupi,mupion)

             mupi_av(h)=mupi_av(h)+mupion
             muk_av(h)=muk_av(h)+mukaon

             cells_t_gt05(h)=cells_t_gt05(h)+1
          endif


          if (t.gt.0.01d0) then
           bar_vol(h)=bar_vol(h)+volce(h,i,j,k)
           bar_totnum(h)=bar_totnum(h)+n0_lrf(h,i,j,k)*volce(h,i,j,k)
          endif

c... get EFFECTIVE BARYON DENSITY (in units of rho_0)
          rhonuc=rhoe_lrf(h,i,j,k)/n_0

          if(baryons.eq.0) rhonuc=0.0d0

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c.. Set cutoff for rho_nuc                                                 !
          if(rhonuc.gt.rhoeffmax) then                                     !
c          write(0,*)'rhonuc set to rhomax',h,i,j,k,rhonuc ! Debug only    !
           rhonuc=rhoeffmax                                                !
          endif                                                            !
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c---------------------------------------------------------------------------

          temp_cell(h,i,j,k)=t
          mub_cell(h,i,j,k)=mub
          mupi_cell(h,i,j,k)=mupion
          muk_cell(h,i,j,k)=mukaon

c---------------------------------------------------------------------------
c.. WRITE OUPUT 68, 69 & 70
c...For central cells
c          write(68,*)'FOR CENTRAL CELLS (x=y=z=0):'
          if(i.eq.int((grd+1)/2).AND.j.eq.int((grd+1)/2).AND.
     &       k.eq.int((grd_z+1)/2)) then
          write(68,721)h,i,j,k,deltat,vxce,vyce,vzce,edens_cell,
     &                  n0_cell,t,mub,mupion,mukaon,lambda,gce
          endif
c...For z-axis (longitudinal)
c          write(69,*)'FOR LONGITUDINAL AXIS (x=y=0):'
          if(i.eq.int((grd+1)/2).AND.j.eq.int((grd+1)/2)) then
           write(69,721)h,i,j,k,deltat,vxce,vyce,vzce,edens_cell,
     &                  n0_cell,t,mub,mupion,mukaon,lambda,gce
          endif
c...For x-axis (transversal)
c          write(70,*)'FOR TRANSVERSAL AXIS (y=z=0):'
          if(j.eq.int((grd+1)/2).AND.k.eq.int((grd_z+1)/2)) then
           write(70,721)h,i,j,k,deltat,vxce,vyce,vzce,edens_cell,
     &                  n0_cell,t,mub,mupion,mukaon,lambda,gce
          endif
c---------------------------------------------------------------------------
c.. DETERMINE AVERAGE FREEZE-OUT TEMPERATURE
          if(focalc) then
           if(fo_marker(i,j,k).eq.h) then
             count_fo=count_fo+1
             avtempfo=avtempfo+((temp_cell(h,i,j,k)+temp_cell(h-1,i,j,k)+temp_cell(h-2,i,j,k))/3.0d0)
           endif !marker
          endif !focalc
c---------------------------------------------------------------------------
c.. Coarsegraining of low-temp. cells to speed up calculation for collisions
c.. at high energies
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          multi=1                                                         !
          if(ecm.gt.3.0d0.AND.ecm.lt.20.0d0.AND.t.le.0.065d0) then        !
           coarsefac=ranff(seedinit)                                      !
           multi=8                                                        !
           crse=1.0d0/dble(multi)                                         !
           if(coarsefac.gt.crse) cycle                                    !
          else if(ecm.gt.3.0d0.AND.ecm.lt.20.0d0.AND.t.gt.0.065d0.AND.    !
     &    t.le.0.80d0) then                                               !
           coarsefac=ranff(seedinit)                                      !
           multi=4                                                        !
           crse=1.0d0/dble(multi)                                         !
           if(coarsefac.gt.crse) cycle                                    !
          else if(ecm.ge.20.0d0.AND.t.le.0.080d0) then                    !
           coarsefac=ranff(seedinit)                                      !
           multi=20                                                       !
           crse=1.0d0/dble(multi)                                         !
           if(coarsefac.gt.crse) cycle                                    !
          else if(ecm.ge.20.0d0.AND.t.gt.0.080d0.AND.t.le.0.100d0) then   !
           coarsefac=ranff(seedinit)                                      !
           multi=10                                                       !
           crse=1.0d0/dble(multi)                                         !
           if(coarsefac.gt.crse) cycle                                    !                                                    !
          else if(ecm.ge.1000.0d0.AND.t.gt.0.10d0.AND.t.le.0.135d0) then  !
           coarsefac=ranff(seedinit)                                      !
           multi=5                                                        !
           crse=1.0d0/dble(multi)                                         !
           if(coarsefac.gt.crse) cycle                                    !
          else if(ecm.ge.1000.0d0.AND.t.gt.0.135d0.AND.t.le.0.250d0) then !
           coarsefac=ranff(seedinit)                                      !
           multi=3                                                        !
           crse=1.0d0/dble(multi)                                         !
           if(coarsefac.gt.crse) cycle                                    !
          end if                                                          !
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c... Calculate DILEPTON / PHOTON  EMISSION
c         write(0,*)'outputs =',outputs ! Debug only
         if(outputs.eq.1) then
          if(t.le.0.95d0) cycle
          if(t.le.0.1d0.AND.rates.eq.-1) cycle
c          write(0,*)'h,i,j,k,temp,mu_b',h,i,j,k,t,mub !Debug only

c... HADRONIC VECTOR-MESON CONTRIBUTION
c DILEPTON
          if(lambda.lt.0.99999d0.AND.rates.ge.0) then
           if(rates.eq.1) then ! RAPP SPECTRAL FUNCTION
            call dilemit_rapp(t,rhonuc,mukaon,mupion,gce,vxce,vyce,vzce,
     &                        vol4,multi,betaLAB,dt,h,lambda)
           elseif(rates.ge.2) then ! RAPP SPECTRAL FUNCTION - RHO & OMEGA / PHI ONLY
             call dilemit_rapp_hr(t,rhonuc,mukaon,mupion,gce,vxce,vyce,vzce,
     &                        vol4,multi,betaLAB,dt,h,lambda)
           elseif(rates.eq.0) then ! ELETSKY SPECTRAL FUNCTION
            call dilemit(vac,t,mub,mupion,gce,vxce,vyce,vzce,vol4,multi,betaLAB,
     &                   dt,h,lambda)
           endif
c PHOTON
          elseif(lambda.lt.0.99999d0.AND.rates.eq.-1) then
            call rapp_photon(t,mub,mupion,gce,vxce,vyce,vzce,
     &                        vol4,multi,betaLAB,dt,h,lambda)
          endif

c... QGP CONTRIBUTION
c DILEPTON
          if(lambda.gt.0.1d-5.AND.rates.ne.3.AND.rates.ge.0) then
           if(latqgp.eq.0) then
             call qgpemit(lambda,t,muq,gce,vxce,vyce,vzce,multi,vol4,
     &                    betaLAB,dt,h)
           elseif(latqgp.eq.1) then
             call qgpemit_lat(lambda,t,muq,gce,vxce,vyce,vzce,multi,vol4,
     &                        betaLAB,dt,h)
           endif
c PHOTON
          elseif(lambda.gt.0.1d-5.AND.rates.eq.-1) then
           call qgp_photon(t,mub,mupion,gce,vxce,vyce,vzce,
     &                  vol4,multi,betaLAB,dt,h,lambda)
          endif

c... 4-PION CONTRIBUTION
c DILEPTON ONLY
          if(lambda.lt.0.99999d0.AND.rates.ne.3.AND.rates.ge.0) then
           if(fourpimix.eq.0) then
            call fopiemit(t,mub,mupion,gce,vxce,vyce,vzce,vol4,multi,
     &                    betaLAB,dt,h,lambda)
           elseif(fourpimix.eq.1) then
            call fopiemit_mix(t,mub,rhonuc,mupion,mukaon,gce,vxce,vyce,
     &                        vzce,vol4,multi,betaLAB,dt,h,lambda)
           endif
          endif

c... BREMSSTRAHLUNG & PI-RHO-OMEGA CONTRIBUTION
C PHOTON ONLY
          if(lambda.lt.0.99999d0.AND.rates.eq.-1) then
           call brems_photon(t,mub,mupion,mukaon,gce,vxce,vyce,vzce,
     &                        vol4,multi,betaLAB,dt,h,lambda)

           call prom_photon(t,mub,mupion,gce,vxce,vyce,vzce,
     &                        vol4,multi,betaLAB,dt,h,lambda)
          endif
         endif

 801     continue
 802    continue
 803   continue
 804  continue !loop over cells/timesteps

      if(eos.eq.3.OR.eos.eq.4.OR.eos.eq.5.OR.eos.eq.6) then
c        write(0,*)'***QGP VOLUME, 4-VOLUME & av. TEMP. in timesteps***'
        do 701 r=1,timesteps
          qgp_avtemp(r)=qgp_avtemp(r)/qgp_4vl(r)
c          write(0,*)r,qgp_vol(r),qgp_4vl(r),qgp_avtemp(r)
 701    continue
        qgp_tottemp=qgp_tottemp/qgp_tot4vl
c        write(0,*)'TOTAL QGP 4-VOLUME & AVERAGE TEMPERATURE:',qgp_tot4vl,qgp_tottemp
      endif
c      write(0,*)'***VOLUME with BARYON CONTENT and BARYON NUMBER***'
      do 702 rr=1,timesteps
c        write(0,*)rr,bar_vol(rr),bar_totnum(rr)
 702  continue
      if(rates.gt.0.AND.pikchem.eq.1) then
c        write(0,*)'***Average PION & KAON CH. POT. in timesteps***'
        do 703 rrr=1,timesteps
         mupi_av(rrr)=mupi_av(rrr)/cells_t_gt05(rrr)
         muk_av(rrr)=muk_av(rrr)/cells_t_gt05(rrr)
c         write(0,*)rrr,mupi_av(rrr),muk_av(rrr)
 703    continue
      endif
      if(focalc) then
        write(0,*)'***Average FREEZE-OUT temperature***'
        avtempfo=avtempfo/dble(count_fo)
        write(0,*)'T_fo =',avtempfo
      endif

      write(0,*),'100%'

c-----|------------------------------------------------------------------
c. **** SHINING FOR LOW-TEMPERATURE CELLS ****
c... By default no dileptons from the Delta(1232) are calculated. To change, set
c... DELEMIT to .TRUE.
c........................!
      delemit=.TRUE.    !
c........................!

c..ONLY FOR DILEPTONS:
      if(rates.ge.0) then
      write(0,*)'CALCULATE DILEPTONS FROM LOW-TEMP. CELLS...'
c... loop over cells
      counter=0
      counterx=0
      do 904 n=1,timesteps
       do 903 o=1,grd

         write(0,'(a1)',advance='NO'),dot
c         write(0,'(a1)'),dot

        do 902 p=1,grd
         do 901 q=1,grd_z

         counter=counter+1
         counterx=counterx+1
         if(counter.eq.100000) then
          percent=100d0*counterx/(timesteps*grd*grd*grd_z)
          counter=0
          write(0,'(f7.2)',advance='NO'),percent
c          write(0,'(f4.0)'),percent
          write(0,*),'%'
         end if

c... Check whether temperature < 0.050 GeV (i.e. cell was not considered for
c... thermal emission)
c         if(outputs.eq.1) then
          if(temp_cell(n,o,p,q).lt.0.95d0) then
           if(norho0(n,o,p,q).gt.0.0d0) then
            call dirrho(dt,betaLAB,evts,n,mrho0(n,o,p,q),p0rho0(n,o,p,q),
     &              pxrho0(n,o,p,q),pyrho0(n,o,p,q),pzrho0(n,o,p,q),
     &              norho0(n,o,p,q))
           endif
           if(nodelta(n,o,p,q).gt.0.0d0.AND.delemit) then
c            write(0,*)'No. of DELTA in cell',n,o,p,q,nodelta(n,o,p,q) !Debug only
            call shdelta(dt,betaLAB,evts,n,mdelta(n,o,p,q),
     &              p0delta(n,o,p,q),pxdelta(n,o,p,q),pydelta(n,o,p,q),
     &              pzdelta(n,o,p,q),nodelta(n,o,p,q))
           endif
           if(noomega(n,o,p,q).gt.0.0d0) then
c            write(0,*)'No. of OMEGAS in cell',n,o,p,q,noomega(n,o,p,q) !Debug only
            call diromega(dt,betaLAB,evts,n,momega(n,o,p,q),
     &              p0omega(n,o,p,q),pxomega(n,o,p,q),pyomega(n,o,p,q),
     &              pzomega(n,o,p,q),noomega(n,o,p,q))
            call dalomega(dt,betaLAB,evts,n,momega(n,o,p,q),
     &              p0omega(n,o,p,q),pxomega(n,o,p,q),pyomega(n,o,p,q),
     &              pzomega(n,o,p,q),noomega(n,o,p,q))
           endif
           if(nophi(n,o,p,q).gt.0.0d0) then
c            write(0,*)'No. of PHIS in cell',n,o,p,q,noomega(n,o,p,q) !Debug only
            call dirphi(dt,betaLAB,evts,n,mphi(n,o,p,q),
     &              p0phi(n,o,p,q),pxphi(n,o,p,q),pyphi(n,o,p,q),
     &              pzphi(n,o,p,q),nophi(n,o,p,q))
           endif
          endif
c         endif

 901     continue
 902    continue
 903   continue
 904  continue !loop over cells/timesteps

      endif ! (rates.ge.0)
c      call setarrayzero
c 888  continue

      write(0,*),'100%'

c-----|------------------------------------------------------------------
c. **** OUTPUT ****
c..This writes the values of T and mu for each cell into a seperate file
c      if(outputs.eq.1) then
      call output(grd,grd_z,dx,deltat,timesteps,tau,p0,n0,edens_lrf,
     &            n0_lrf,jx,jy,jz,mub_cell,temp_cell,gam)
c      endif
c-----|------------------------------------------------------------------

c. **** FORMAT DEFINITIONS ****
c format: file13, UrQMD 2.3
 213  format(9e16.8,i11,2i3,i9,i5,i4,8e16.8)
c    LHC
c 213  format(9e24.16,i11,2i3,i9,i5,i4,8e24.16)

c formats: file14, UrQMD 2.3
 214  format(9e16.8,i11,2i3,i9,i5,i4)
 2141 format(i4,9e16.8,i11,2i3,i9,i5,i4)
c    LHC
 244  format(9e24.16,i11,2i3,i9,i5,i4)
 2441 format(i4,9e24.16,i11,2i3,i9,i5,i4)

c header-line for each collision in file14

 131  format(a20,3i7,a15,i2)
 132  format(a13,a13,i4,i4,a12,a13,i4,i4,a1)
 133  format(a36,3f11.7)
 134  format(a36,3f6.2,a31,f9.2)
 135  format(a20,i3,a15,e11.4,a15,e11.4,a15,e11.4)
 136  format(a7,i9,a13,i12,a9,a20,i7,a20,f11.3)
 137  format(a2,15(i3,a2))
 138  format(a2,12(e11.4,a2))
 139  format(a171)

 301  format(2i12)
 302  format(8i8)

c format for output 72
 555  format(e14.7,4f12.7,i9,2f12.7)
 556  format(I5,1X,13(E16.8E3,1X),I4,1X,2(E16.8E3,1X),I9,1X,
     &F12.8,1X,2(E16.8,1x))

c format for output 74
 511  format(i5,2x,e14.8,2x,e14.8)
 599  format(10(a12,4x))

c format for output 68, 69 & 70
 721  format(4(i4,2x),10(e14.7,2x),2(e14.4,2x))

c format for standard output header
 761  format(F6.3,2x,I5,2x,F6.3,2x,2(I5,2x),F9.2,2x,F12.9,2(2x,F6.2))

      return
      end

c*****|****************************************************************|X
