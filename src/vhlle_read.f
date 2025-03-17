c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE INPUT
c..
c.. This is the main coarse-graining subroutine
c.. for UrQMD file14 (with multiple timesteps) read-in,
c.. calling the other calculation routines
c*****|****************************************************************|X

      subroutine vhlle_read(grd,grd_z,dx,vol,vac)
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
      character*32 control,dummylong
      
      write(0,*)'CODE HIJACKED'
      write(0,*)'****************************************************'
      
      t=0.15
      rhonuc=0.2
      muq=0.3
      mukaon=0
      mupion=0
      gce=0.1
      vxce=0
      vyce=0
      vzce=0.3
      vol4=1
      multi=100
      betaLAB=0.3
      dt=1
      h=1
      lambda=0.15
      do h=1,500
       if(mod(h,100).eq.0) write(0,*)h
       call qgpemit_lat(lambda,t,muq,gce,vxce,vyce,vzce,multi,vol4,
     &                  betaLAB,dt,h)
       call dilemit_rapp_hr(t,rhonuc,mukaon,mupion,gce,vxce,vyce,vzce,
     &                     vol4,multi,betaLAB,dt,h,lambda)
       call fopiemit_mix(t,mub,rhonuc,mupion,mukaon,gce,vxce,vyce,
     &                   vzce,vol4,multi,betaLAB,dt,h,lambda)
      end do

      write(0,*)'CALLED DILETPON EMISSIONS'
      write(0,*)'****************************************************'
c... HADRONIC VECTOR-MESON CONTRIBUTION
c DILEPTON
          if(lambda.lt.0.99999d0.AND.rates.ge.0) then
           if(rates.eq.1) then ! RAPP SPECTRAL FUNCTION
            call dilemit_rapp(t,rhonuc,mukaon,mupion,gce,vxce,vyce,vzce,
     &                        vol4,multi,betaLAB,dt,h,lambda)
           elseif(rates.ge.2) then ! RAPP SPECTRAL FUNCTION - RHO & OMEGA / PHI ONLY or all 3
             call dilemit_rapp_hr(t,rhonuc,mukaon,mupion,gce,vxce,vyce,vzce,
     &                        vol4,multi,betaLAB,dt,h,lambda)
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

c-----|------------------------------------------------------------------
c. **** OUTPUT ****
c..This writes the values of T and mu for each cell into a seperate file
c      if(outputs.eq.1) then
c      call output(grd,grd_z,dx,deltat,timesteps,tau,p0,n0,edens_lrf,
c     &            n0_lrf,jx,jy,jz,mub_cell,temp_cell,gam)
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
