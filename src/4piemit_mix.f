c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE FOPIEMIT_MIX
c..
c..   This program returns 4-pion dilepton rates rates dR/dM 
c..
c..         mu: chemical potential [GeV]
c..          T: temperature [GeV]
c..          M: invariant mass [GeV]
c..      dR/dM: in-medium dilepton rate [1/GeV]
c..   
c..   dN/dM=4*alpha_em**2/pi**2*T*BessK1(M/T)*Im(PI_em)
c..
c..   with
c..   Pi_em=1/3*M**2*WF
c..
c..   WF is the electromagnetic current current correlator
c..
c..   To obtain dR/dM^2 (often plotted in papers) 
c..   one has to further divide for 2M,i.e.     
c..   dR/dM^2=(dR/dM)/(2M)
c..
c*****|****************************************************************|X

      SUBROUTINE fopiemit_mix(temp,mub,rhonuc,pipot,kpot,gce,vxce,vyce,
     &                     vzce,vol4,multi,beta_lab,dt,time,lambda)

      implicit none

      include 'defs.f'

c-----|----------------------------------------------------------------|X

      real*8 lambda,temp,mub,rhonuc,pipot,kpot,gce,vxce,vyce,vzce
      real*8 mass,p0l,pxl,pyl,pzl
      integer multi,stat
      real*8 contr,vol4
      real*8 factor,rate,result
      real*8 mmin4pi,mmax4pi
c      parameter(mmin4pi=0.7875d0,mmax4pi=2.5625d0)
      parameter(mmin4pi=0.560d0,mmax4pi=2.75d0)
 
      integer flag4pi
      parameter (flag4pi=4)

      real*8 dt,time

      integer noe,ityp
      real*8 bev,acce,accp
      real*8 effvol4

      real*8 p0_po_lab,px_po_lab,py_po_lab,pz_po_lab
      real*8 p0_el_lab,px_el_lab,py_el_lab,pz_el_lab
      real*8 beta_lab

c     functions
      real*8 distm4pi_mix
      external distm4pi_mix

c     common variables
      real*8 T,rhn,ppt,kpt
      common /th_dyn_r/T,rhn,ppt,kpt

      integer i,l

c-----|----------------------------------------------------------------|X

c     rename variables and put in common
      T=temp
      rhn=rhonuc
      ppt=pipot
      kpt=kpot

c     minimum mass for integration: mmin4pi

c     maximum mass for integration:  mmax4pi

c     calculate mass and momentum integrated rate
c     "rate" has units [fm^{-4}]

c     integration over mass dependent function--Simpson integration
      rate=0d0
      call qsimp_had(distm4pi_mix,mmin4pi,mmax4pi,result)
c     constant factor
      factor=4.0d0*alpha_em**2/(3.0d0*pi**2)
      rate=result*factor
c     conversion to fm^{-4}
      rate=rate/(hqc**4)

c      write(0,*) 'In FOPIEMIT_MIX:' !Debug only
c      write(0,*) 'rate',rate !Debug only

c. **** Loop n-times for better mass and momentum statsitics ****
c.. Number of Loops
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
      if(t.gt.0.1d0.AND.dexp(pipot/temp).gt.2.0d0) then     !
        stat=3                                              !
      elseif(t.lt.0.1d0.AND.dexp(pipot/temp).gt.3.0d0) then !
        stat=3                                              ! 
      else                                                  !
        stat=10                                             !
      endif                                                 !                              
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc! 
      do 678 l=1,stat

c. **** Generate mass ****
      CALL fourpi_mdist_mix(mass,mmin4pi,mmax4pi)

c      write(0,*) 'mass',mass !Debug only

c     in case of too many iterations 
c     mass is set to zero in massdist
      if(mass.eq.0) return

c. **** Generate 3-momenta ****
      CALL fourpi_momdist_mix(mass,gce,vxce,vyce,
     &        vzce,p0l,pxl,pyl,pzl)

c     in case of too many interations
c     p0l is set to zero in momdist
      if(p0l.eq.0) return

c      write(0,*)'p0l,pxl,pyl,pzl',p0l,pxl,pyl,pzl !Debug only

      CALL pel_ppo_dist(beta_lab,mass,p0l,pxl,pyl,pzl,p0_el_lab,
     &px_el_lab,py_el_lab,pz_el_lab,p0_po_lab,px_po_lab,py_po_lab,
     &pz_po_lab)     

c      write(0,*)'p0e,pxe,pye,pze',p0_el_lab,px_el_lab,py_el_lab,
c     &                            pz_el_lab ! Debug only 
c      write(0,*)'p0e,pxe,pye,pze',p0_po_lab,px_po_lab,py_po_lab,
c     &                            pz_po_lab ! Debug only 

c.. Determine cell contribution
      contr=rate*(1d0-lambda)*dble(multi)/dble(stat)*vol4

       if(contr.lt.0.0d0.OR.contr.gt.(5.0d0*vol4)) then
c        write(0,*)'Suspicious weight!'
c        write(0,*)'Weight,vol4,T,rhnuc,ppt',contr,vol4,T,rhn,pipot
        contr=0.0d0
        return
       endif

c. **** Write into output file f71 ****
      ityp=444
      noe=1
      bev=0.0d0
      acce=1.0d0
      accp=1.0d0
      effvol4=vol4*(1d0-lambda)*dble(multi)/dble(stat)
c........................................................................
c. NOTE: The extended output format writes out ***lab-system momenta*** !
c.       for e+ and e-                                                  !
c.......................................................................!
      if(vHLLE_out) then
       write(71,*)time,ityp,contr/multi,mass,p0l,pxl,pyl,pzl,effvol4/multi,t,mub
c       write(0,*)time,ityp,contr/multi,mass,p0l,pxl,pyl,pzl,effvol4/multi,t,mub
      elseif(ext_out) then !extended output format
       write(71,556)ityp,contr,mass,p0_el_lab,px_el_lab,py_el_lab,
     & pz_el_lab,p0_po_lab,px_po_lab,py_po_lab,pz_po_lab,dt,time,effvol4,
     & flag4pi,mub,temp,noe,lambda,acce,accp
      else !standard output format
       write(71,555)contr,mass,pxl,pyl,pzl,flag4pi,mub,temp
      endif

 555  format(e14.7,4f12.7,i9)	 
 556  format(I5,1X,13(E16.8E3,1X),I4,1X,2(E16.8E3,1X),I9,1X,
     &F12.8,1X,2(E16.8,1x))
 557  format(I3,I3,e14.7,5f12.7,3f6.3)

 678  continue

      return
      end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE FOURPI_MDIST_MIX
c..
c..   This subroutine generates the dilepton mass for 4pi emission.
c..
c..   Input : mmin4pi,mmax4pi (max and min values for the mass)
c..   Output: m (the dilepton mass)
c..
c*****|****************************************************************|X

      SUBROUTINE fourpi_mdist_mix(m,mmin4pi,mmax4pi)
      implicit none

c-----|----------------------------------------------------------------|X

      real*8 mmin4pi,mmax4pi,h,hmax,hm
      real*8 m,mtest
      real*8 ranff
      integer i,j
      real*8 aux
      real*8 zbar
      real*8 hqc
      real*8 truef
      parameter (hqc  = 0.197327) ! value of $\hbar c$

c     functions
      real*8 distm4pi_mix,bessk1

c     common variables
      real*8 T,rhn,ppt,kpt
      common /th_dyn_r/T,rhn,ppt,kpt
      integer seedinit
      common /randomf/ seedinit
      
c-----|----------------------------------------------------------------|X      
      
c      write(0,*)'mmin,mmax,T,rhn,ppt,kpt',mmin4pi,mmax4pi,T,rhn,ppt,kpt !Debug only

      truef=0.0d0
      hmax=0.0d0
      mtest=mmin4pi
      do 100 i=1,30
        mtest=mtest+(i*0.025d0)
        hm=distm4pi_mix(mtest)
c       write(0,*)'m,distm',mtest,hm !Debug only
        hmax=max(hmax,hm)
c       write(0,*)'HMAX =',hmax !Debug only
 100  continue

      hmax=1.5d0*hmax
 
      j=0
 2    continue
      j=j+1
      m=mmin4pi+ranff(seedinit)*(mmax4pi-mmin4pi)
c     write(*,*) 'generated mass:', m
      h=ranff(seedinit)*hmax
      
      truef=distm4pi_mix(m)

      if(truef.gt.hmax) then
         aux=truef
         write(0,*)'hmax too small'
         write(0,*)'hmax',hmax
         write(0,*)'distm4pi_mix(m)',aux
         write(0,*)'maximum violated by', (aux-hmax)/hmax
c         if((aux-hmax)/hmax.gt.1.d-3) stop
c      stop
      endif
      
c     mass could not be generated
      if(j.gt.10000000) then
         m=0.d0
         write(0,*)'from fourpi_mdist_mix:'
         write(0,*)'too many trials to generate mass'
         write(0,*) 'number of trials',j
         return
      endif
      
      if(h.gt.truef) goto 2
      
      return
      end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE FOURPI_MOMDIST_MIX
c..
c..   This subroutine generates the dilepton momenta for 4pi emission
c..
c*****|****************************************************************|X

      SUBROUTINE fourpi_momdist_mix(m,g,vx,vy,vz,p0,px,py,pz) 
      implicit none

c-----|----------------------------------------------------------------|X

      real*8 m,g,vx,vy,vz,p0,px,py,pz
      integer i
      real*8 ranff
      real*8 ftrial
      real*8 q0,qx,qy,qz,q,phi,cost,sint,v2
      real*8 xmin,xmax,x,x2,c,compf,truef
      real*8 pi
      parameter (pi = 3.1415926535) ! useful constant
      real*8 q0up,q0down

c     common variables
      real*8 qup,qdown
      common /qupdow/ qup,qdown
      real*8 T,rhn,ppt,kpt
      common /th_dyn_r/T,rhn,ppt,kpt
      integer seedinit
      common /randomf/ seedinit

c-----|----------------------------------------------------------------|X      

c.. maximum energy allowed from tabulated momenta
      q0up=sqrt(qup**2+m**2)
c.. minimum energy
      q0down=m

c     apply transformation+rejection method in order to generate
c     the energy q0

c     first, the transformation method is applied to
c     generate a random deviate q0 according to 
c     the distribution  exp(-q0/T)

      xmin=0.d0
      c=-T*exp(-q0down/T)   
      xmax=-T*exp(-q0up/T)-c

c     start random generation

c     generate random phi
      phi=ranff(seedinit)*2.d0*pi
c     generate random cos($\theta$)
      cost=-1d0+ranff(seedinit)*(1d0+1.d0)

      i=0
 1    continue
      i=i+1
c..   generate q0 according to the distribution
c     exp(-q0/T)
      x=ranff(seedinit)*xmax
      q0=-(T*Log((T*Exp(-m/T)-x)/T))
      if((T*Exp(-m/T)-x)/T.le.0d0) then
         write(*,*) 'numerical rumor in fourpi_momdist_mix'
         write(*,*) 'extracted log of non positive number'
         write(*,*)'logarg',(T*Exp(-m/T)-x)/T
         write(*,*) 'x,T,m',x,T,m
         stop
      endif
c..    in principle(mathematically) the obtained q0 is 
c..    always q0.ge.m... but numerics plays bad games
c..    when q0 is about m. This is the reason of this
c..    conditional line
       if(q0**2.lt.m**2) q0=m
       q=sqrt(q0**2-m**2)

c     now apply rejection method in order to generate the
c     distribution p(x)=q*exp(-q0/T). The comparison function 
c     is f(x)=qup*exp(-q0/T). The distribution f(x) has been
c     already generated by transformation method (up to the 
c     constant factor qup).

      compf=qup*exp(-q0/T)
      truef=q*exp(-q0/T)
c..   added for stability
      if(truef.le.0.d0) then
         p0=0.d0
         return
      endif

      if(compf.lt.truef) then
         write(*,*) 'from fourpi_momdist_mix:'
         write(*,*)'comparison function too small'
         write(*,*) 'ERROR'
         stop
      endif
      x2=ranff(seedinit)
      if(x2.gt.truef/compf.and.i.lt.1000000) then
         goto 1
      endif
      if(i.ge.1000000) then
         write(*,*) 'too many iteration in fourpi_momdist_mix' 
         p0=0.d0
         return
      endif

c     find cartesiane coordinates
       sint=sqrt(1.d0-cost**2)
       qx=q*cos(phi)*sint
       qy=q*sin(phi)*sint
       qz=q*cost
c     boost to cf
       if(g.gt.1.0d0) then
       v2=vx**2+vy**2+vz**2
              
       p0=g*q0+vx*g*qx+vy*g*qy+vz*g*qz
       px=vx*g*q0
     &    +(1.0d0+(g-1.0d0)*vx**2/v2)*qx
     &    +(g-1.0d0)*vx*vy/v2*qy
     &    +(g-1.0d0)*vx*vz/v2*qz
       py=vy*g*q0
     &    +(g-1.0d0)*vy*vx/v2*qx
     &    +(1.0d0+(g-1.0d0)*vy**2/v2)*qy
     &    +(g-1.0d0)*vy*vz/v2*qz
       pz=vz*g*q0
     &    +(g-1.0d0)*vz*vx/v2*qx
     &    +(g-1.0d0)*vz*vy/v2*qy
     &    +(1.0d0+(g-1.0d0)*vz**2/v2)*qz
       else
       p0=q0
       px=qx
       py=qy
       pz=qz
       endif
       return 
       end

c-----|----------------------------------------------------------------|X      
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. FUNCTION DISTM4PI_MIX
c..
c..   This program returns the dR/dM distribution function 
c..   (without constant factors).
c..
c..   mass : dilepton mass (in GeV)
c..
c*****|****************************************************************|X

      real*8 function distm4pi_mix(mass)
      implicit none

c-----|----------------------------------------------------------------|X

      real*8 mass
      real*8 WF
      real*8 hqc
      parameter (hqc  = 0.197327) ! value of $\hbar c$

c     functions
      real*8 bessk1
      complex*16 PiVmix

c     common variables
      real*8 T,rhn,ppt,kpt
      common /th_dyn_r/T,rhn,ppt,kpt

c-----|----------------------------------------------------------------|X

      WF=3.0d0*dimag(PiVmix(dcmplx(mass**2)))

      distm4pi_mix=mass**2*T*bessk1((mass/T))*WF
      end

c-----|----------------------------------------------------------------|X
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c.. THE FOLLOWING ROUTINES ARE PROVIDED BY HENDRIK VAN HEES
c..
c.. They are used to calculate the multi-pion contribution of the 
c.. thermal dilepton radiation acoording to vector/axial-vector
c.. correlators from tau-decay data provided by ALEPH [Eur.Phys.J. C4, 
c.. 409-431 (1998)].
c..
c.. The mixing here is according to Dey, Eletsky and Ioffe [Phys.Lett. 
c.. B252, 620 (1990)] and takes the following form [H. van Hees, 
c.. R. Rapp, Nucl. Phys. A 806, 339 (2008), Equation (12)]:
c..
c..   Pi_V,mix=[(1-epsilon)*z_pi**4*Pi_V,4pi]+
c..            [(epsilon/2)*z_pi**3*Pi_A,3pi]+
c..            [(epsilon/2)*(z_pi**4+z_pi**5)*Pi_A,5pi]
c..
c.. where epsilon=(1/2)*e(T,mu_pi)/e(T_c,0)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function PiV(s)

c     Good only for imaginary part!!!

      implicit none

      double complex s,FV,i
      double precision eA,delA

      double precision gapi

      real*8 pi,mass_pion
      parameter (pi=3.141592653589793d0)
      parameter (mass_pion=0.139d0)

      real*8 T,rhn,ppt,kpt
      common /th_dyn_r/T,rhn,ppt,kpt

      gapi=dexp(ppt/T)

c      eA=1.235d0 ! Ralf's value
c      delA=0.14d0 !Ralf's value

c      eA=1.26d0
c      delA=0.12d0

      eA=1.15d0
      delA=0.11d0

      i=(0.0d0,1.0d0)

      PiV=(0.0d0,0.0d0)

      if (dble(s) .gt. 16.0d0*mass_pion**2) then
      PiV=PiV+i/(8.0d0*pi)*(1.0d0+0.22d0/log(1.0d0+sqrt(s)/0.2d0))
     $        /(1.0d0+exp((eA-sqrt(s))/delA))*
     $        (1.0d0-(4.0d0*mass_pion)**2/s)**2
     $        *gapi**4
      end if

      end


c-----|----------------------------------------------------------------|X      
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. FUNCTION PIVMIX
c..
c.. This routine calculates the mixing effect on the (isovector) 
c.. four-pion part with the mixing paramter. The pion chemical 
c.. potentials are implemented as an overall fugacity factor, 
c.. z_pi = exp(mu_pi/T) in Boltzmann approximation [see SUBROUTINE 
c.. MUPIK in temp.f] (they do not here, but already in FUNCTIONS PIV,
c.. PIA3PI and PiA5PI)
c.. 
c.. The two-pion piece, as well as the three-pion piece corresponding to
c.. a1 decay, a_1 --> pion + rho, have been removed as they are
c.. included via the rho-spectral function.
c..
c.. The expression for the mixing effect on the vector-isovector 
c.. current correlation function is according to H. van Hees, 
c.. R. Rapp, Nucl. Phys. A 806, 339 (2008), Equation (12)
c..
c*****|****************************************************************|X

      double complex function PiVmix(s)

      implicit none

c-----|----------------------------------------------------------------|X 

      double complex s,PiV,PiA3pi,PiA5pi
      double precision eps,tadpii

c..   tadpic=tadpii(T_C)=tadpii(0.175GeV)=0.49513467
      real*8 tadpic
      parameter (tadpic=0.49513467d0)

      double precision gapi

      real*8 T,rhn,ppt,kpt
      common /th_dyn_r/T,rhn,ppt,kpt

c-----|----------------------------------------------------------------|X 

      gapi=dexp(ppt/T)

      eps=tadpii(T)/(2.0d0*tadpic)

c      Tc=0.175d0

      PiVmix=(1.0d0-eps)*PiV(s)+
     $     0.5d0*eps*(PiA3pi(s)+PiA5pi(s)+PiA5pi(s)/gapi)

      end

c-----|----------------------------------------------------------------|X 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function PiA3pi(s)

      implicit none


      double precision gapi

      double complex s,sthr,i
      double complex a0,a1,a3,a4
      double precision a2

      real*8 pi
      parameter (pi=3.141592653589793d0)

      real*8 T,rhn,ppt,kpt
      common /th_dyn_r/T,rhn,ppt,kpt

      gapi=dexp(ppt/T)

c     Imaginary part own fit

      sthr=(0.55d0,0.0d0)
      i=(0.0d0,1.0d0)
      a0=0.0305d0
      a1=0.1764d0
      a2=47.692d0
      a3=1.1981d0
      a4=0.620959d0

c      eA=1.55d0 !my values without threshold
c      delA=0.16d0

      if (dble(s) .gt. dble(a1)) then
        PiA3pi=i*pi*a0/s*((1.0d0,0.0d0)-a1**2/s**2)**a2/
     $         ((s-a3)**2+a4**2)*gapi**3
      else
         piA3pi=(0.0d0,0.0d0)
      end if

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function PiA5pi(s)

      implicit none


      double precision gapi

      double complex s,sthr,i
      double complex eA,delA

      real*8 pi,mass_pion
      parameter (pi=3.141592653589793d0)
      parameter (mass_pion=0.139d0)

      real*8 T,rhn,ppt,kpt
      common /th_dyn_r/T,rhn,ppt,kpt

      gapi=dexp(ppt/T)

c     Imaginary part own fit

      sthr=(0.55d0,0.0d0)
      i=(0.0d0,1.0d0)

      eA=1.6d0
      delA=0.14d0

      if (dble(s) .gt. 25.0d0*mass_pion**2) then
       PiA5pi=i/(8.0d0*pi)*(1.0d0+0.22d0/log(1.0d0+sqrt(s)/0.2d0))
     $        /(1.0d0+exp((eA-sqrt(s))/delA))
     $        *(1.0d0-(5.0d0*mass_pion)**2/s)**2.5d0
     $        *gapi**5
      else
         piA5pi=(0.0d0,0.0d0)
      end if

      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. FUNCTION PITAD_INTERPOL
c..
c.. Interpolates the tabulated pion-tadpole diagram at given temperature
c..
c*****|****************************************************************|X

      real*8 function tadpii(temp)
      implicit none

c-----|----------------------------------------------------------------|X

      integer k
      real*8 temp

c-----|----------------------------------------------------------------|X
c.. TEMP BINS
      real*8 tbin(151)
      data tbin /0.050,0.051,0.052,0.053,0.054,0.055,0.056,0.057,0.058,
     &0.059,0.060,0.061,0.062,0.063,0.064,0.065,0.066,0.067,0.068,0.069,
     &      0.070,0.071,0.072,0.073,0.074,0.075,0.076,0.077,0.078,0.079,
     &      0.080,0.081,0.082,0.083,0.084,0.085,0.086,0.087,0.088,0.089,
     &      0.090,0.091,0.092,0.093,0.094,0.095,0.096,0.097,0.098,0.099,
     &      0.100,0.101,0.102,0.103,0.104,0.105,0.106,0.107,0.108,0.109,
     &      0.110,0.111,0.112,0.113,0.114,0.115,0.116,0.117,0.118,0.119,
     &      0.120,0.121,0.122,0.123,0.124,0.125,0.126,0.127,0.128,0.129,
     &      0.130,0.131,0.132,0.133,0.134,0.135,0.136,0.137,0.138,0.139,
     &      0.140,0.141,0.142,0.143,0.144,0.145,0.146,0.147,0.148,0.149,
     &      0.150,0.151,0.152,0.153,0.154,0.155,0.156,0.157,0.158,0.159,
     &      0.160,0.161,0.162,0.163,0.164,0.165,0.166,0.167,0.168,0.169,
     &      0.170,0.171,0.172,0.173,0.174,0.175,0.176,0.177,0.178,0.179,
     &      0.180,0.181,0.182,0.183,0.184,0.185,0.186,0.187,0.188,0.189,
     &      0.190,0.191,0.192,0.193,0.194,0.195,0.196,0.197,0.198,0.199,
     &      0.200/

c.. PI-TADPOLE BINS
      real*8 pitad(151)
      data pitad/0.052484015,0.055966439,0.0595307621,0.0631723997,
     &0.0668868373,0.0706695309,0.074516008,0.0784218749,0.0823827521,
     &0.0863944054,0.0904526758,0.0945535052,0.0986929449,0.1028671600,
     &0.1070724190,0.1113051630,0.1155619050,0.1198393050,0.1241341670,
     &0.1284433720,0.1327639820,0.1370931680,0.1414282320,0.1457665990,
     &0.1501058190,0.1544435520,0.1587776230,0.1631059280,0.1674264950,
     &0.1717374690,0.1760371100,0.1803237510,0.1845958700,0.1888520320,
     &0.1930909010,0.1973112360,0.2015118890,0.2056917890,0.2098499840,
     &0.2139855730,0.2180977480,0.2221857740,0.2262489930,0.2302868290,
     &0.2342987420,0.2382842890,0.2422430800,0.2461747870,0.2500791290,
     &0.2539559060,0.2578049480,0.2616261410,0.2654194210,0.2691847700,
     &0.2729222130,0.2766318180,0.2803136990,0.2839679880,0.2875948710,
     &0.2911945660,0.2947673140,0.2983134060,0.3018331440,0.3053268640,
     &0.3087949290,0.3122377250,0.3156556630,0.3190491750,0.3224187150,
     &0.3257647570,0.3290877980,0.3323883370,0.3356669000,0.3389240400,
     &0.3421603080,0.3453762740,0.3485725200,0.3517496430,0.3549082490,
     &0.3580490380,0.3611724700,0.3642792650,0.3673700680,0.3704455290,
     &0.3735063200,0.3765530970,0.3795865380,0.3826073250,0.3856161430,
     &0.3886136840,0.3916005660,0.3945776490,0.3975455580,0.4005050040,
     &0.4034567000,0.4064013600,0.4093397110,0.4122724720,0.4152003700,
     &0.4181241380,0.4210445010,0.4239621970,0.4268779620,0.4297925360,
     &0.4327066590,0.4356210740,0.4385365290,0.4414537700,0.4443735480,
     &0.4472966130,0.4502237200,0.4531556240,0.4560930830,0.4590368570,
     &0.4619877120,0.4649464070,0.4679137090,0.4708903920,0.4738772250,
     &0.4768749820,0.4798844410,0.4829063830,0.4859415890,0.4889908450,
     &0.4920549410,0.4951346700,0.5027521360,0.5104316360,0.5181732020,
     &0.5259768660,0.5338426600,0.5417706160,0.5497607660,0.5578131390,
     &0.5659277680,0.5741046810,0.5823439070,0.5906454740,0.5990094140,
     &0.6074357520,0.6159245160,0.6244757340,0.6330894320,0.6417656360,
     &0.6505043720,0.6593056650,0.6681695420,0.6770960260,0.6860851410,
     &0.6951369140,0.7042513650/
c-----|----------------------------------------------------------------|X
      
      if(temp.lt.0.051d0) then
        tadpii=pitad(1)
      elseif(temp.gt.0.200d0) then
        tadpii=pitad(151)
      else
        k=0
 110    k=k+1
        if(temp.ge.tbin(k)) then
           goto 110
        end if
        tadpii=pitad(k-1)+(pitad(k)-pitad(k-1))/
     &         (tbin(k)-tbin(k-1))*(temp-tbin(k-1))
      endif
      end

c-----|----------------------------------------------------------------|X
