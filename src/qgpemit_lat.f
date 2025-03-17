c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE QGPEMIT_lat
c..
c.. This program returns dilepton rates dR/dM for the QUARK-GLUON
c.. PLASMA according to lattice calculations. The according fit was
c.. provided by HENDRIK VAN HEES.
c..
c..         mu: chemical potential [GeV]
c..          T: temperature [GeV]
c..          M: invariant mass [GeV]
c..      dR/dM: in-medium dilepton rate [1/GeV]
c..   
c.. The QGP dilepton rate is calculated according to the expression
c..
c.. dN/(d^4xd^4q) = -alpha_EM**2*dilphsp(M)/(pi**3*M**2)*f_B(q_0;T)
c..                     *ImPi_EM(M,q;mu_B,T)
c..
c.. which corresponds to eq.1 in arXiv:1304.2309.
c..
c.. As d^4q=dE*d^3q and M*dM=E*dE, we can transform this into
c..
c.. dN/(d^4x*d^3q*dM) = dN/(d^4xd^4q) * M/E
c..
c.. and in consequence get (with \int d^3q=4pi \int q**2 dq):
c..
c.. dN/(d^4x*dM)=M/E*(\int d^3q dR/d^4q) = 4*pi*(\int dq q**2 dR/d^4q)
c..
c*****|****************************************************************|X

      SUBROUTINE qgpemit_lat(lam,temp,muquark,gce,vxce,vyce,vzce,
     &  multi,vol4,beta_lab,dt,timestep) 

      implicit none
      include 'defs.f'

c-----|----------------------------------------------------------------|X

      real*8 lam,temp,muquark,gce,vxce,vyce,vzce
      real*8 mass,p0l,pxl,pyl,pzl
      integer pospt,posm
      real*8 pt
      integer multi,stat

c     common variables
      real*8 T,muq
      common /th_dyn/T,muq
      common /leptype/m_lept 
 
      integer flagqgp
      parameter (flagqgp=2)

      real*8 tmin_qgplat,tmax_qgplat
      parameter (tmin_qgplat=0.100d0)
      parameter (tmax_qgplat=0.490d0)

      integer save_posmaxpt 
      real*8 factor,rate !,res
c     functions
      real*8 distmqgp_lat
      real*8 fqgp_lat
      external distmqgp_lat
      external fqgp_lat
	    
      integer noe,ityp
      real*8 bev,acce,accp
      real*8 effvol4

      real*8 dt,time
      integer timestep

      real*8 ammax,ammin
      real*8 result,dam,am(100),s,hp,fam1
      integer nsimp,nsimp1,ic
      integer i1

      real*8 beta_lab,p0_el_lab,px_el_lab,py_el_lab,pz_el_lab
      real*8 p0_po_lab,px_po_lab,py_po_lab,pz_po_lab

      real*8 contr,vol4

      integer l

      logical zwstep
      real*8 zwmin(1:4),zwmax(1:4) 
      integer j

c-----|----------------------------------------------------------------|X
c      write(0,*)'***IN QGPEMIT_LAT***' !Debug only

c.. Rename variables and put in common
      T=temp
      muq=muquark

c.. First check whether T is in temperature range
      if(T.lt.tmin_qgplat) return
      if(T.gt.tmax_qgplat) T=tmax_qgplat

c.. Determine the mass of the lepton-type and put it in common

      if(dimuon) then
         m_lept=mass_muon
      else !dielectron
         m_lept=mass_electron
      endif

c.. Invariant mass threshold (not put in common, as possibly different from mmin)
       ammin=2.0d0*m_lept*1.000001d0
       if(.NOT.dimuon) ammin=0.01d0
       ammax=2.75d0

c.. MASS STEPS
      zwmin(1)=ammin
      zwmax(1)=0.45d0   
      zwmin(2)=0.4500001d0
      zwmax(2)=0.9d0   
      zwmin(3)=0.9000001d0
      zwmax(3)=1.8d0
      zwmin(4)=1.8000001d0
      zwmax(4)=ammax

c-----|----------------------------------------------------------------|X
c. **** If: zwstep=.TRUE. : Loop to increase statistics at high mass ****
c*********************!   
      zwstep=.true.   !
c     zwstep=.false.  !
c*********************!

c.. Set counter to zero
      j=0

 123  continue

      if(zwstep) then
        j=j+1
        ammin=zwmin(j)
        ammax=zwmax(j)
c        write(0,*)'j,ammin,ammax',j,ammin,ammax !Debug only
      endif

c-----|----------------------------------------------------------------|X

c.. calculate mass and momentum integrated rate
c.. NOTE: "rate" has units [fm^{-4}]!

c.. Integration over mass dependent function - SIMPSON INTEGRATION
         rate=0d0

      call qsimp_had(distmqgp_lat,ammin,ammax,result)

c.................................................................

      factor=alpha_em**2/(pi**3)

      rate=result*factor
c.. Conversion to fm^{-4}
      rate=rate/(hqc**4)

c      write(0,*)'temp,rate',temp,rate !Debug only

c-----|----------------------------------------------------------------|X
   
c. **** Loop n-times for better mass and momentum statsitics ****

c.. Number of Loops (for better stistics, especially at
c.. low beam energies at FAIR, where only few cells emit
c.. dileptons from the QGP)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
      if(na60mode.eq.1.OR.beta_lab.gt.0.99996d0) then    !
         stat=1                                          !
      else if(na60mode.eq.0.AND.beta_lab.lt.0.94d0) then !
         stat=15                                         !
      else                                               !
         stat=10                                         !
      end if                                             !         
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc! 

      do 678 l=1,stat

c.. ****GENERATE MASS****
      call massdistqgp_lat(mass,ammin,ammax)

c      write(0,*)'mass',mass !Debug only

c.. In case of too many iterations 
c.. mass is set to zero in massdistqgp_lat
      if(mass.eq.0) then
       write(0,*)"mass is zero"
       return
      endif

c.. ****GENERATE 3-MOMENTA****
      call qgpmomdist_lat(gce,vxce,vyce,vzce,p0l,pxl,pyl,pzl)

c.. In case of too many interations
c.. p0l is set to zero in qgpmomdist_lat
      if(p0l.eq.0) then
       write(0,*)"p0L is zero"
       return
      endif
c      write(0,*)'p0l,pxl,pyl,pzl',p0l,pxl,pyl,pzl ! Debug only 

c.. Determine momenta of e+ and e-

      CALL pel_ppo_dist(beta_lab,mass,p0l,pxl,pyl,pzl,p0_el_lab,
     &px_el_lab,py_el_lab,pz_el_lab,p0_po_lab,px_po_lab,py_po_lab,
     &pz_po_lab)     

c      write(0,*)'p0e,pxe,pye,pze',p0_el_lab,px_el_lab,py_el_lab,
c     &                            pz_el_lab ! Debug only 
c      write(0,*)'p0p,pxp,pyp,pzp',p0_po_lab,px_po_lab,py_po_lab,
c     &                            pz_po_lab ! Debug only 

c-----|----------------------------------------------------------------|X

c.. ****WRITE CELL CONTRIBUTION****
      contr=rate*lam*dble(multi)/dble(stat)*vol4
c      write(0,*)'contr=',contr

       if(contr.lt.0.0d0.OR.contr.gt.(5.0d0*vol4)) then
        write(0,*)'Suspicious weight!'
        write(0,*)'Weight,vol4,T,mu_q',contr,vol4,T,muquark
        contr=0.0d0
        return
       endif

c      write(0,*)'contr, mass, p',contr,mass,p0l,pxl,pyl,pzl

c. **** Write into output file f71 ****
      ityp=222
      noe=1
      bev=0.0d0
      acce=1.0d0
      accp=1.0d0
      time=dble(timestep*dt)
      effvol4=vol4*lam*dble(multi)/dble(stat)
c........................................................................
c. NOTE: The extended output format writes out ***lab-system momenta*** !
c.       for e+ and e-                                                  !
c.......................................................................!
      if(vHLLE_out) then
       write(71,557)ityp,contr/multi,mass,p0l,pxl,pyl,pzl,effvol4/multi,t,3*muquark
      elseif(ext_out) then !extended output format
c       write(71,556)ityp,contr,mass
       write(71,*)ityp,contr,mass,p0_el_lab,px_el_lab,py_el_lab,
     & pz_el_lab,p0_po_lab,px_po_lab,py_po_lab,pz_po_lab,dt,time,effvol4,
     & flagqgp,muquark,temp,noe,lam,acce,accp
      else !standard output format
       write(71,555)contr,mass,pxl,pyl,pzl,flagqgp,muquark,temp
      endif

 555  format(e14.7,4f12.7,I2,2E14.7)
 556  format(I5,1X,13(E16.8E3,1X),I4,1X,2(E16.8E3,1X),I9,1X,
     &F12.8,1X,2(E16.8,1x))
 557  format(I3,e14.7,5f12.7,3f6.3)

 678  continue

c-----|----------------------------------------------------------------|X

      if(zwstep.AND.j.lt.4) goto 123

      return
      end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE MASSDISTQGP_LAT
c..
c.. This subroutine generates the dilepton mass for QGP emission.
c..
c*****|****************************************************************|X

      SUBROUTINE massdistqgp_lat(mass,ammin,ammax)
      implicit none
      
      include 'defs.f'

c-----|----------------------------------------------------------------|X

      real*8 temp,muquark,mass
      real*8 fmax,f,truef,ranff,distmqgp_lat
      integer i     
      real*8 ammin,ammax
      
c.. Common variables
      real*8 T,muq
      common /th_dyn/T,muq

      real*8 qgpmommaxtab(0:80,0:1499),qgpmassmaxtab(0:400)
      common /qgpmaxmom/ qgpmommaxtab,qgpmassmaxtab

c-----|----------------------------------------------------------------|X

c.. Rename variables put in common
      temp=T
      muquark=muq

c.. Maximum of the distribution (approximated)
c      fmax=0.95d0*temp**3
       if(ammin.lt.0.5d0) then
        fmax=qgpmassmaxtab(min(int((T-0.100d0)/0.001d0)+1,400))
       else
        fmax=2.0d0*distmqgp_lat(ammin)
       endif
c      write(0,*) 'maximum of the distribution', fmax !Debug only

c-----|----------------------------------------------------------------|X

      i=0
      truef=0.0d0

 1    continue
      i=i+1

      mass=ammin+ranff(seedinit)*(ammax-ammin)
      f=ranff(seedinit)*fmax

      truef=distmqgp_lat(mass)

      if(truef.gt.fmax)then
         write(0,*)'in massdistqgp_lat:'
         write(0,*)'fmax too small'
         write(0,*) 'fmax',fmax
         write(0,*) 'fp',truef
      endif 

      if(f.gt.truef.and.i.le.10000000) then
        go to 1
      end if
      if(i.gt.10000000)then
         mass=0.0d0
         stop 'too many iterations in massdistiqgp'
      end if 
c       write(0,*) 'number of trials in massdistqgp_lat:', i  !Debug only

       return 
       end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c..  FUNCTION DISTMQGP_LAT
c..
c..  dR/dM distribution function (without constant factors)
c..
c..  mass : dilepton mass (in GeV) 
c..
c*****|****************************************************************|X

      real*8 FUNCTION distmqgp_lat(mass)
      implicit none

c-----|----------------------------------------------------------------|X

      real*8 res,res2 
      real*8 dilphsp
      real*8 mass 

c     momentum dependent function
      real*8 fqgp_lat
      real*8 bessk1
      external fqgp_lat
      external bessk1

c     common variables
      real*8 am,T,muq
      common /th_dyn/T,muq
      common /dilmass/am

      real*8 m_lept
      common /leptype/m_lept 
     
c-----|----------------------------------------------------------------|X

      distmqgp_lat=0.d0 !initialization

c     rename variable and put in common
      am=mass

c     dilepton phase-space
      dilphsp=sqrt(1-4.d0*m_lept**2/am**2)*(1+2.d0*m_lept**2/am**2)
     
      if(1-4.d0*m_lept**2/am**2.lt.0) return

c     integrate momentum dependent funtion with routine
c     qromb(func,a,b,ss)

c      call qromb(fqgp_lat,1.d-25,sqrt((30d0*T+am)**2-am**2),res)
      call qromb(fqgp_lat,0.001d0,5.0d0,res)

      distmqgp_lat=dilphsp*res/(am**2)   
   
      return
      end

c234567**1*********2*********3*********4*********5*********6*********7**X
c.. FUNCTION FQGP_LAT
c..
c.. THIS FUNCTION IS PROVIDED BY HENDRIK VAN HEES 
c.. (has been modified for use in Coarse-Graining routine)
c..
c.. Momentum dependent function (to be integrated)
c.. Input variables:  
c..
c..  q : modulus of the dilepton 3-momentum (in GeV)
c..
c.. It calculates a 'formfactor' fitted to lattice QCD correlators for
c.. 3-momentum q=0 (Kaczmarek et al. 2010) and with a lightlike limit
c.. consistent with the LO alpha_s photon production rate (augmented 
c.. by K-factor)
c..
c*****|****************************************************************|X

      real*8 FUNCTION fqgp_lat(q)
      implicit none

c-----|----------------------------------------------------------------|X

c     VM energy, VM 3-momentum (modulus)
      real*8 q0,q     

c     boson thermal distrib.
      real*8 fb

c     common variables
      real*8 am,T,muq
      common /th_dyn/T,muq
      common /dilmass/am
      
      real*8 CEM,fhat2,xpl,xmi,K,La,F,als,W
      real*8 pi
      parameter (pi=3.141592653589793d0)

c-----|----------------------------------------------------------------|X
      
      if(T.eq.0.) then
         fqgp_lat=0.0d0
         return
      endif
      
      q0=sqrt(q*q+am*am)
      fb=1.d0/(exp(q0/T)-1.d0)


      CEM=5.0d0/9.0d0
      xpl=exp(-(q0+q)/(2.0d0*T))
      xmi=exp(-(q0-q)/(2.0d0*T))
      fhat2=1.0d0+2.0d0*T/q*dlog((1.0d0+xpl)/(1.0d0+xmi))
      La=2.0d0*T
      F=La**2/(La**2+am**2)
      K=2.0d0
      als=6.0d0*pi/(28.0d0*dlog(T/0.022d0))  

c.. W corresponds to $\im \Pi_{\text{EM}}$ in Ralf Rapp's paper,
c.. eq.(1), arXiv:1304.2309v2
   
      W=CEM*3.0d0/(12.0d0*pi)*am**2*
     $     (fhat2+2.0d0*pi*als*T**2/am**2*K*F*
     $     (2.0d0+am**2/q0**2)/3.0d0*
     $     dlog(1.0d0+2.912d0*q0/(4.0d0*pi*als*T)))

c.. dN/(dx^4dM)=M/q0 \int d^3q dR/d^4q
c.. \int d^3q=4pi \int p**2 dp

      fqgp_lat=4.0d0*pi*W*fb*am/q0*q**2

      return
      end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE QGPMOMDIST_LAT
c..
c.. This subroutine generates the dilepton momenta for QGP emission.
c..
c*****|****************************************************************|X

      SUBROUTINE qgpmomdist_lat(g,vx,vy,vz,p0,px,py,pz) 
      implicit none

      include 'defs.f'

c-----|----------------------------------------------------------------|X

      real*8 g,vx,vy,vz,v2,p0,px,py,pz,mu
      real*8 qx,qy,qz,q0,qel
      integer a,i
      real*8 fmax,f,fp,distpqgp_lat,ptmax,pzmax
      real*8 ranff,den 		
      real*8 qlocal,pmunu

c     common variables
      real*8 am,T,muq
      common /th_dyn/T,muq
      common /dilmass/am
      real*8  time,  acttime, bdist, ebeam, bimp,bmin,ecm

      real*8 pxmax,pxmin,pymax,pymin,pzmin,lhs,rsph
      real*8 pel,pmax,pmin,ranx,phi_dil,theta_dil,sin_theta_dil

      real*8 truef

      real*8 qgpmommaxtab(0:80,0:1499),qgpmassmaxtab(0:400)
      common /qgpmaxmom/ qgpmommaxtab,qgpmassmaxtab

c-----|----------------------------------------------------------------|X

c.. maximum pt and pz range for particle generation       
c       ptmax=ecm
c       pzmax=2.d0*ecm

c.. bose distribution       
c       a=-1
c.. no chemical potential for dil-pair (boson-like)
c       mu=0.d0

c.. maximum of the distribution is approximated with the maximum 
c.. of the boson distribution
c       den=vx**2+vy**2+vz**2-1.0d0
       
        px=0.0d0
        py=0.0d0
        pz=0.0d0
        p0=0.0d0

c       px=abs(vx)/vx*sqrt(-vx**2*am**2/den)
c       py=abs(vy)/vy*sqrt(-vy**2*am**2/den)
c       pz=abs(vz)/vz*sqrt(-vz**2*am**2/den)
c       p0=sqrt(am**2+px**2+py**2+pz**2)

c       fmax=fp(p0,px,py,pz,vx,vy,vz,g,a,T,mu)


        if(am.gt.1.7d0.AND.T.le.0.280d0.AND.T.gt.0.220d0) then
          fmax=qgpmommaxtab(min(int((T-0.100d0)/0.005d0),80),
     &                      min(int((am-0.001d0)/0.002d0),1499))*1.175d0!*1.65d0
        else if(am.gt.1.1d0.AND.T.le.0.220d0) then
          fmax=qgpmommaxtab(min(int((T-0.100d0)/0.005d0),80),
     &                      min(int((am-0.001d0)/0.002d0),1499))*1.25d0!*1.65d0
        else
          fmax=qgpmommaxtab(min(int((T-0.100d0)/0.005d0),80),
     &                      min(int((am-0.001d0)/0.002d0),1499))*1.1d0!*1.65d0
        endif

c       write(0,*)'IN QGPMOMDIST_LAT: ecm,pzmax,fmax',ecm,pzmax,fmax !Debug only

c       write(0,*) 'momenta at max:',px,py,pz,p0 !Debug only
       i=0

c-----|----------------------------------------------------------------|X

c       define "radius" for cut at exp(-30)*fmax
c        rsph=sqrt((30.d0*T+am)**2-am**2)
c        lhs=exp(30.d0)*(exp(am/T)+a)-a
c        if(lhs.lt.0.d0) then
c       	  write(0,*)'ERROR:negative argument for log!'
c          stop
c        endif
c        if((T*log(lhs))**2-am**2.lt.0) then
c	  write(0,*)'ERROR:sqrt of negative nuber!'
c          stop
c        endif

c..	radius
c        rsph=sqrt((T*log(lhs))**2-am**2)  
       
c        write(0,*)'T,mu,radius',T,mu,rsph
 
c        call lorsphere(rsph,g,vx,vy,vz,am,pxmax,pxmin,
c     & pymax,pymin,pzmax,pzmin)
         
        pmin=0.01d0
        pmax=3.5d0
         
ccccccccccccccccccccccccccc       check     cccccccccccccccccccccccccccX
c        if(px.gt.pxmax.or.px.lt.pxmin) then
c            write(0,*)'smt wrong with max/min values of px'
c            write(0,*)'px at max:',px,'pxmax,pxmin:',pxmax,pxmin
c            stop
c        endif
c        if(py.gt.pymax.or.py.lt.pymin) then
c            write(0,*)'smt wrong with max/min values of py'
c            write(0,*)'py at max:',py,'pymax,pymin:',pymax,pymin
c            stop
c        endif
c        if(pz.gt.pzmax.or.pz.lt.pzmin) then
c            write(0,*)'smt wrong with max/min values of pz'
c            write(0,*)'pz at max:',pz,'pzmax,pzmin:',pzmax,pzmin
c            stop
c        endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccX
   
 1     continue
       i=i+1
 
c       qel=(ranff(seedinit)*(pmax-pmin))+pmin

c       phi_dil=2d0*pi*ranff(seedinit)
c       ranx=ranff(seedinit)     

c       sin_theta_dil=(ranx-0.5d0)*2d0    
c       theta_dil=asin(sin_theta_dil)+pi/2d0

c...in the CELL REST FRAME
c       qx=sin(theta_dil)*cos(phi_dil)*qel               
c       qy=sin(theta_dil)*sin(phi_dil)*qel
c       qz=cos(theta_dil)*qel

c...determine *absolute* values for q
        qx=(ranff(seedinit)*(pmax-pmin))+pmin
        qy=(ranff(seedinit)*(pmax-pmin))+pmin
        qz=(ranff(seedinit)*(pmax-pmin))+pmin

c...determign sign of momenta
        if(ranff(seedinit).gt.0.5d0) qx=-1.0d0*qx
        if(ranff(seedinit).gt.0.5d0) qy=-1.0d0*qy
        if(ranff(seedinit).gt.0.5d0) qz=-1.0d0*qz

        qel=sqrt(qx**2+qy**2+qz**2)

        q0=sqrt(am**2+qx**2+qy**2+qz**2)

c...check
       if(q0**2.le.am**2) goto 1 

c        write(0,*) 'generated momenta:',qx,qy,qz,q0,am,i  !Debug only
        f=ranff(seedinit)*fmax
c        write(0,*)'fmax,f',fmax,f !Debug only
       
c       pmunu=g*(p0-px*vx-py*vy-pz*vz)
 
c       if(pmunu**2.le.am**2) goto 1 

c       qlocal=sqrt(pmunu**2-am**2)
        qlocal=qel

c        if(qlocal.gt.pmax.OR.qlocal.lt.pmin) then
c         write(0,*)'IN QGPEMIT_LAT: qlocal out of range!'
c         write(0,*)'qx,qy,qz,q0,qlocal',qx,qy,qz,q0,qlocal
c         stop
c        endif

c        write(0,*)'qlocal,q0',qlocal,q0 ! Debug only
        truef=distpqgp_lat(qlocal)
c        write(0,*)'distpqgp_lat=',truef ! Debug only

c       write(0,*) 'pmunu,qlocal',pmunu,qlocal
      
        if(truef.gt.fmax)then
        write(0,*)'in qgpmomdist_lat:'
        write(0,*)'fmax too small'
        write(0,*)'fmax',fmax
        write(0,*)'fp',truef
        write(0,*)'mass,T',am,T
       endif 

       
       if(f.gt.truef.and.i.le.1500000) then
        go to 1
       end if
       if(i.gt.1500000)then
        q0=0.0d0
        write(0,*)'too many iterations in qgpmomdist_lat'
        write(0,*)'mass,T,fmax',am,T,fmax
        write(0,*)'g,vx,vy,vz',g,vx,vy,vz
        write(0,*)'pxmin,pxmax,pymin,pymax,pzmin,pzmax',pxmin,pxmax,pymin,pymax,pzmin,pzmax
c        stop 'too many iterations in qgpmomdist_lat'
        return
       end if 
c       write(0,*) 'number of trials in qgpmomdist_lat:', i ! Debug only

c     boost to cf
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

c       write(0,*) 'COM momenta:',px,py,pz,p0  !Debug only
      

       return 
       end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE DISTPQGP_LAT
c..
c.. dR/(dMd^3q) distribution function (without constant factors)
c..
c.. Input variables:  
c.. q : modulus of the dilepton 3-momentum in the thermal 
c..     rest frame  --i.e. p_/mu u^/mu -- (in GeV)
c.. p0: energy of the dilepton in the frame the mom. distr. is 
c..     generated, e.g. typically c.m. frame
c..
c*****|****************************************************************|X

      real*8 FUNCTION distpqgp_lat(q)
      implicit none

c-----|----------------------------------------------------------------|X

c     VM energy, VM 3-momentum (modulus)
      real*8 q0,q,p0   
c     boson thermal distrib.
      real*8 fb
      real*8 dilphsp

c     common variables
      real*8 am,T,muq
      common /th_dyn/T,muq
      common /dilmass/am
      
      real*8 m_lept
      common /leptype/m_lept 

      real*8 ImD_lat
      real*8 CEM,fhat2,xpl,xmi,K,La,F,als
      real*8 pi
      parameter (pi=3.141592653589793d0)

      
c-----|----------------------------------------------------------------|X

      if(T.eq.0.0d0) then
         distpqgp_lat=0.
         return
      endif

c     dilepton phase-space
      dilphsp=sqrt(1-4.d0*m_lept**2/am**2)*(1+2.d0*m_lept**2/am**2)     
      if(1-4.d0*m_lept**2/am**2.lt.0) return      
      
      q0=sqrt(q*q+am*am)
      fb=1.d0/(exp(q0/T)-1.d0)


      CEM=5.0d0/9.0d0
      xpl=exp(-(q0+q)/(2.0d0*T))
      xmi=exp(-(q0-q)/(2.0d0*T))
      fhat2=1.0d0+2.0d0*T/q*dlog((1.0d0+xpl)/(1.0d0+xmi))
      La=2.0d0*T
      F=La**2/(La**2+am**2)
      K=2.0d0
      als=6.0d0*pi/(28.0d0*dlog(T/0.022d0))  
   
      ImD_lat=CEM*fb*3.0d0/(12.0d0*pi)*am**2*
     $     (fhat2+2.0d0*pi*als*T**2/am**2*K*F*
     $     (2.0d0+am**2/q0**2)/3.0d0*
     $     dlog(1.0d0+2.912d0*q0/(4.0d0*pi*als*T)))

c      distpqgp_lat=ImD_lat/(am**2)*dilphsp*am/p0
      distpqgp_lat=ImD_lat/(am**2)*dilphsp*am/q0!*(q0/p0)

c      write(0,*)'IN DISTPQGP_LAT: q,q0,p0',q,q0,p0 !Debug only
 
      return
      end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE QGPTABLES_LAT
c..
c.. Reads maxima for 
c..
c*****|****************************************************************|X

      SUBROUTINE qgptables_lat
      implicit none

      include 'defs.f'

c-----|----------------------------------------------------------------|X

      real*8 Targ,amarg
      real*8 dummy

      real*8 qgpmommaxtab(0:80,0:1499),qgpmassmaxtab(0:400)
      common /qgpmaxmom/ qgpmommaxtab,qgpmassmaxtab

      integer k,kk
      integer l

c-----|----------------------------------------------------------------|X

      open(unit=8,file='diltables/qgp/table_distmqgpmax_ext.dat')
      open(unit=9,file='diltables/qgp/table_distpqgpmax_ext.dat')

      do 881 l=0,400,1
        if(dimuon) then
          read(8,99) targ,dummy,qgpmassmaxtab(l)
        else !dielectron
          read(8,99) targ,qgpmassmaxtab(l),dummy
        endif
c        write(0,*)'T,l,mmax',targ,l,qgpmassmaxtab(l) ! Debug only
 881  continue

      close(8)

      do 991 k=0,80,1
       do 992 kk=0,1499,1
        read(9,99) targ,amarg,qgpmommaxtab(k,kk)
c        write(0,*)'T,m,k,kk,qmax',targ,amarg,k,kk,qgpmommaxtab(k,kk) ! Debug only
 992   continue
 991  continue

      close(9)

 99   format(2x,e14.8,2x,e14.8,2x,e14.8)

      end 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
c..                                                             !  
c..     **** INTEGRATION ROUTINES ****                          !
c..      --> SEE FILE ***qgpemit.f***                           !
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!

