c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE QGPEMIT
c..
c.. This program returns dilepton rates dR/dM for the QUARK-GLUON
c.. PLASMA.
c..
c*****|****************************************************************|X

      SUBROUTINE qgpemit(lam,temp,muquark,gce,vxce,vyce,vzce,
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
      logical htl_corr
      common /th_dyn/T,muq
      common /leptype/m_lept 
      common /htl/htl_corr
 
      integer flagqgp
      parameter (flagqgp=2)

      integer save_posmaxpt 
      real*8 factor,rate !,res
c     functions
      real*8 distmqgp
	    
      integer noe,ityp
      real*8 bev,acce,accp
      real*8 effvol4

      real*8 dt,time
      integer timestep

      real*8 ammin,ammax
      real*8 result,dam,am(100),s,hp,fam1
      integer nsimp,nsimp1,ic
      integer i1

      real*8 beta_lab,p0_el_lab,px_el_lab,py_el_lab,pz_el_lab
      real*8 p0_po_lab,px_po_lab,py_po_lab,pz_po_lab

      real*8 contr,vol4

      integer l

c.. Steps to increase statistics
      logical zwstep
      common /steps/zwstep
      real*8 zwmin(1:4),zwmax(1:4) 
      integer j

c-----|----------------------------------------------------------------|X

c.. Rename variables and put in common
      T=temp
      muq=muquark
      htl_corr=htlcorr
   
c     write(0,*)'IN QGPEMIT: muq=',muq ! Debug only

c.. Determine the mass of the lepton-type and put it in common

      if(dimuon) then
         m_lept=mass_muon
      else !dielectron
         m_lept=mass_electron
      endif

c.. Invariant mass threshold (not put in common)
      ammin=2.0d0*m_lept*1.000001d0
      if(.NOT.dimuon) ammin=0.01d0
      ammax=3.0d0

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
      endif

c-----|----------------------------------------------------------------|X

c.. calculate mass and momentum integrated rate
c.. NOTE: "rate" has units [fm^{-4}]!

c.. Integration over mass dependent function - SIMPSON INTEGRATION
           rate=0d0
           nsimp=40
          
           dam=(ammax-ammin)/nsimp
           nsimp1=nsimp+1
           do i1=1,nsimp1
              am(i1)=ammin+dam*(i1-1)
           enddo

           fam1=distmqgp(am(1))

           s=0d0
           do i1=2,nsimp
              ic=2*(i1-2*((i1-1)/2))
              s=s+ic*distmqgp(am(i1))
           enddo
           result=dam/3.d0*(s+fam1+distmqgp(am(nsimp1)))

c.. Fractional charge for 2 quark flavors (u,d): (1d0/9d0+4d0/9d0)
      factor=alpha_em**2/4.d0/pi**4*(1d0/9d0+4d0/9d0)

      rate=result*factor*4.d0*pi
c.. Conversion to fm^{-4}
      rate=rate/(hqc**4)

c-----|----------------------------------------------------------------|X
c. **** Loop again n-times for better mass and momentum statsitics ****


c.. Number of Loops (for better stistics, especially at
c.. low beam energies at FAIR, where only few cells emit
c.. dileptons from the QGP)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
      if(na60mode.eq.1.OR.eos.eq.3.OR.eos.eq.5) then     !
         stat=1                                          !
      else if(na60mode.eq.0.AND.beta_lab.lt.0.94d0.AND.  !
     &        (.NOT.(eos.eq.3.OR.eos.eq.5))) then        !
         stat=10                                         !
      else                                               !
         stat=2                                          !
      end if                                             !         
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc! 

      do 678 l=1,stat

c.. ****GENERATE MASS****
      call massdistqgp(mass,ammin,ammax)

c      write(0,*)'mass',mass !Debug only

c.. In case of too many iterations 
c.. mass is set to zero in massdistqgp
      if(mass.eq.0) return

c.. ****GENERATE 3-MOMENTA****
      call qgpmomdist(gce,vxce,vyce,vzce,p0l,pxl,pyl,pzl)

c.. In case of too many interations
c.. p0l is set to zero in qgpmomdist
      if(p0l.eq.0) return

c.. Determine momenta of e+ and e-

      CALL pel_ppo_dist(beta_lab,mass,p0l,pxl,pyl,pzl,p0_el_lab,
     &px_el_lab,py_el_lab,pz_el_lab,p0_po_lab,px_po_lab,py_po_lab,
     &pz_po_lab)     

c      write(0,*)'p0e,pxe,pye,pze',p0_el_lab,px_el_lab,py_el_lab,
c     &                            pz_el_lab ! Debug only 
c      write(0,*)'p0e,pxe,pye,pze',p0_po_lab,px_po_lab,py_po_lab,
c     &                            pz_po_lab ! Debug only 

c-----|----------------------------------------------------------------|X

c.. ****WRITE CELL CONTRIBUTION****
      contr=rate*lam*float(multi)/dble(stat)*vol4

       if(contr.lt.0.0d0.OR.contr.gt.(5.0d0*vol4)) then
        write(0,*)'Suspicious weight!'
        write(0,*)'Weight,vol4,T,mu_q',contr,vol4,T,muquark
        contr=0.0d0
        return
       endif

c. **** Write into output file f71 ****
      ityp=222
      noe=1
      bev=0.0d0
      acce=1.0d0
      accp=1.0d0
      time=dble(timestep*dt)
      effvol4=vol4*lam*float(multi)/dble(stat)
c........................................................................
c. NOTE: The extended output format writes out ***lab-system momenta*** !
c.       for e+ and e-                                                  !
c.......................................................................!
      if(vHLLE_out) then
       write(71,557)ityp,contr,mass,p0l,pxl,pyl,pzl,t,muquark
      elseif(ext_out) then !extended output format
       write(71,556)ityp,contr,mass,p0_el_lab,px_el_lab,py_el_lab,
     & pz_el_lab,p0_po_lab,px_po_lab,py_po_lab,pz_po_lab,dt,time,effvol4,
     & flagqgp,muquark,temp,noe,lam,acce,accp
      else !standard output format
       write(71,555)contr,mass,pxl,pyl,pzl,flagqgp,muquark,temp
      endif

 555  format(e14.7,4f12.7)
 556  format(I5,1X,13(E16.8E3,1X),I4,1X,2(E16.8E3,1X),I9,1X,
     &F12.8,1X,2(E16.8,1x))
 557  format(I3,e14.7,5f12.7,2f6.3)

 678  continue

c-----|----------------------------------------------------------------|X

      if(zwstep.AND.j.lt.4) goto 123

      return
      end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE MASSDISTQGP
c..
c.. This subroutine generates the dilepton mass for QGP emission.
c..
c*****|****************************************************************|X

      SUBROUTINE massdistqgp(mass,ammin,ammax)
      implicit none
      
      include 'defs.f'

c-----|----------------------------------------------------------------|X

      real*8 temp,muquark,mass
      real*8 fmax,f,ranff,distmqgp
      integer i     
      real*8 ammin,ammax
      
c.. Common variables
      real*8 T,muq
      common /th_dyn/T,muq
c.. Minimum mass (see file "comdil.f")

      real*8 mtest
      integer count

c-----|----------------------------------------------------------------|X

c.. Rename variables put in common
      temp=T
      muquark=muq

c.. Maximum of the distribution (approximated)
c      fmax=0.95d0*temp**3
       fmax=0.0d0
       if(ammin.lt.0.5d0) then
        count=0
 7011   continue
        count=count+1
        mtest=max(0.0d0,(ammin-0.002d0)+count*0.02d0)
        fmax=max(fmax,3.0d0*distmqgp(mtest))
        if(mtest.lt.ammax) goto 7011
       else
        fmax=3.0d0*distmqgp(ammin)
       endif
c      write(0,*) 'maximum of the distribution', fmax !Debug only

c-----|----------------------------------------------------------------|X

      i=0
 1    continue
      i=i+1

      mass=ammin+ranff(seedinit)*(ammax-ammin)
      f=ranff(seedinit)*fmax

      if(distmqgp(mass).gt.fmax)then
         write(0,*) 'IN MASSDISTQGP: fmax too small'
         write(0,*) 'fmax',fmax
         write(0,*) 'f',distmqgp(mass)
      endif 

      if(f.gt.distmqgp(mass).and.i.le.10000000) then
        go to 1
      end if
      if(i.gt.10000000)then
         mass=0.0d0
         stop 'too many iterations in massdistiqgp'
      end if 
c       write(0,*) 'number of trials in massdistqgp:', i  !Debug only

       return 
       end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c..  FUNCTION DISTMQGP
c..
c..  dR/dM distribution function (without constant factors)
c..
c..  mass : dilepton mass (in GeV) 
c..
c*****|****************************************************************|X

      real*8 FUNCTION distmqgp(mass)
      implicit none

c-----|----------------------------------------------------------------|X

      real*8 res,res2 
      real*8 dilphsp
      real*8 mass 

c     momentum dependent function
      real*8 fqgp 
      external fqgp 

c     common variables
      real*8 am,T,muq
      common /th_dyn/T,muq
      common /dilmass/am

      real*8 m_lept
      common /leptype/m_lept 
     
c-----|----------------------------------------------------------------|X

      distmqgp=0.d0 !initialization

c     rename variable and put in common
      am=mass

c     dilepton phase-space
      dilphsp=sqrt(1-4.d0*m_lept**2/am**2)*(1+2.d0*m_lept**2/am**2)
     
      if(1-4.d0*m_lept**2/am**2.lt.0) return

c     integrate momentum dependent funtion with routine
c     qromb(func,a,b,ss) and  midinf(funk,aa,bb,s,n)

c      call  qromb(fqgp,0.d0,3.d0,res)
c      call midinf(fqgp,3.d0,1.d2,res2,10)
      
c      if(res.lt.0d0.or.res2.lt.0.d0) return
c      
c      res=res+res2

c..	changed
	call qromb(fqgp,0.d0,sqrt((30d0*T+am)**2-am**2),res)


c..   dR/dM=2*M*dR/dM^2
      distmqgp=res*dilphsp*2.d0*mass
     
      return
      end

c234567**1*********2*********3*********4*********5*********6*********7**X
c.. FUNCTION FQGP
c..
c.. Momentum dependent function (to be integrated)
c.. Input variables:  
c..
c..  q : modulus of the dilepton 3-momentum (in GeV)
c..
c*****|****************************************************************|X

      real*8 FUNCTION fqgp(q)
      implicit none

c-----|----------------------------------------------------------------|X

c     VM energy, VM 3-momentum (modulus)
      real*8 q0,q     

c     boson thermal distrib.
      real*8 fb,f
      real*8 xmin,xmax

c     hard thermal loop
      real*8 htlcor
      external htlcor

c     common variables
      real*8 am,T,muq
      logical htl_corr
      common /th_dyn/T,muq
      common /dilmass/am
      common /htl/htl_corr
      
c-----|----------------------------------------------------------------|X
      
      if(T.eq.0.) then
         fqgp=0.
         return
      endif
      
      q0=sqrt(q*q+am*am)
      fb=1.d0/(exp(q0/T)-1.d0)

      xmin=exp(-(q0-q)/2.d0/T)
      xmax=exp(-(q0+q)/2.d0/T)

      f=(xmin+exp(-(q0+muq)/T))/(xmax+exp(-(q0+muq)/T))*
     &  (xmax+exp(-muq/T))/(xmin+exp(-muq/T))
     
c      f=log(f)*fb*T/q

      f=log(f)*fb*T

c     transformation to invariant mass distributions
c
c     dN/(dx^4dM^2)=1/(2q0)/int d^3q dR/d^4q

c      fqgp=f*q*q/2.d0/q0

      if(htl_corr) fqgp=f*q/2.d0/q0*htlcor(q0,T)   
      if(.NOT.htl_corr) fqgp=f*q/2.d0/q0   

c      write(*,*) '# q, fqgp(q)'
c      write(*,*) q,fqgp

      return
      end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE QGPMOMDIST
c..
c.. This subroutine generates the dilepton momenta for QGP emission.
c..
c*****|****************************************************************|X

      SUBROUTINE qgpmomdist(g,vx,vy,vz,p0,px,py,pz) 
      implicit none

      include 'defs.f'

c-----|----------------------------------------------------------------|X

      real*8 g,vx,vy,vz,p0,px,py,pz,mu
      integer a,i
      real*8 fmax,f,fp,distpqgp,ptmax,pzmax
      real*8 ranff,den 		
      real*8 qlocal,pmunu

c     common variables
      real*8 am,T,muq
      common /th_dyn/T,muq
      common /dilmass/am
      real*8  time,  acttime, bdist, ebeam, bimp,bmin,ecm
      common /rsys/ time,acttime,bdist,bimp,bmin,ebeam,ecm

c     hard thermal loop
      real*8 htlcor
      external htlcor

c     TEST!!!!
c      real*8 ecm
c      parameter (ecm=5d0)

c      new
       real*8 pxmax,pxmin,pymax,pymin,pzmin,lhs,rsph
       logical new_cooper 

       real*8 minmom,minp0

c-----|----------------------------------------------------------------|X

       new_cooper=.true.

c.. maximum pt and pz range for particle generation       
       ptmax=ecm
       pzmax=2.d0*ecm
c.. bose distribution       
       a=-1
c.. no chemical potential for dil-pair (boson-like)
       mu=0.d0

c.. maximum of the distribution is approximated with the maximum 
c.. of the boson distribution
       den=vx**2+vy**2+vz**2-1.0d0
       
       px=abs(vx)/vx*sqrt(-vx**2*am**2/den)
       py=abs(vy)/vy*sqrt(-vy**2*am**2/den)
       pz=abs(vz)/vz*sqrt(-vz**2*am**2/den)
       p0=sqrt(am**2+px**2+py**2+pz**2)


       fmax=fp(p0,px,py,pz,vx,vy,vz,g,a,T,mu)
      
c       write(0,*) 'momenta at max:',px,py,pz,p0 !Debug only
       i=0

c-----|----------------------------------------------------------------|X
c     ////////// added  ////////////////////////
        if(new_cooper) then

c       define "radius" for cut at exp(-30)*fmax
c        rsph=sqrt((30.d0*T+am)**2-am**2)
        lhs=exp(30.d0)*(exp(am/T)+a)-a
        if(lhs.lt.0.d0) then
       	  write(0,*)'ERROR:negative argument for log!'
          stop
        endif
        if((T*log(lhs))**2-am**2.lt.0) then
	  write(0,*)'ERROR:sqrt of negative nuber!'
          stop
        endif

c..	radius
        rsph=sqrt((T*log(lhs))**2-am**2)  
       
c        write(*,*)'T,mu,radius',T,mu,rsph
 
        call lorsphere(rsph,g,vx,vy,vz,am,pxmax,pxmin,
     & pymax,pymin,pzmax,pzmin)

c----|----------------------------------------------------------------|X

c        minmom=min(dabs(pxmin),dabs(pymin))
c        minmom=min(dabs(minmom),dabs(pzmin))
        minmom=0.01d0
        minp0=dsqrt(am**2+minmom**2)

c        fmax=1.0d2*fmax*htlcor(minp0,T)
        fmax=distpqgp(minmom,minp0)

ccccccccccccccccccccccccccc       check     cccccccccccccccccccccccccccX
        if(px.gt.pxmax.or.px.lt.pxmin) then
            write(0,*)'smt wrong with max/min values of px'
            write(0,*)'px at max:',px,'pxmax,pxmin:',pxmax,pxmin
            stop
        endif
        if(py.gt.pymax.or.py.lt.pymin) then
            write(0,*)'smt wrong with max/min values of py'
            write(0,*)'py at max:',py,'pymax,pymin:',pymax,pymin
            stop
        endif
        if(pz.gt.pzmax.or.pz.lt.pzmin) then
            write(0,*)'smt wrong with max/min values of pz'
            write(0,*)'pz at max:',pz,'pzmax,pzmin:',pzmax,pzmin
            stop
        endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccX
c----|----------------------------------------------------------------|X

        endif
c    ////////////////////////////////////////////////////////////

   
 1     continue
       i=i+1


       if(new_cooper) then
c       //////////////  changed  ///////////////////////////////
          px=pxmin+ranff(seedinit)*(pxmax-pxmin)
          py=pymin+ranff(seedinit)*(pymax-pymin)
          pz=pzmin+ranff(seedinit)*(pzmax-pzmin)
c       /////////////////////////////////////////////////////////
       else
          px=-ptmax+2.*ranff(seedinit)*ptmax
          py=-ptmax+2.*ranff(seedinit)*ptmax
          pz=-pzmax+2.*ranff(seedinit)*pzmax
       endif

       p0=sqrt(am**2+px**2+py**2+pz**2)

c       write(*,*) 'generated momenta:',px,py,pz,p0,m,i 
       f=ranff(seedinit)*fmax
       
       pmunu=g*(p0-px*vx-py*vy-pz*vz)

       if(pmunu**2.le.am**2) goto 1 

       qlocal=sqrt(pmunu**2-am**2)

c       write(*,*) 'pmunu,qlocal',pmunu,qlocal
      
       if(distpqgp(qlocal,pmunu).gt.fmax)then
        write(0,*) 'IN QGPMOMDIST: fmax too small'
        write(0,*) 'fmax',fmax
        write(0,*) 'fp',distpqgp(qlocal,pmunu)
        write(0,*) 'mass,T,muq',am,T,muq
       endif 
       
       if(f.gt.distpqgp(qlocal,pmunu).and.i.le.1000000) then
        go to 1
       end if
       if(i.gt.1000000)then
        p0=0.0d0
        write(0,*) 'too many iterations in qgpmomdist'
c        stop 'too many iterations in qgpmomdist'
       end if 
c       write(*,*) 'number of trials in qgpmomdist:', i

       return 
       end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE DISTQGP
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

      real*8 FUNCTION distpqgp(q,p0)
      implicit none

c-----|----------------------------------------------------------------|X

c     VM energy, VM 3-momentum (modulus)
      real*8 q0,q,p0   
c     boson thermal distrib.
      real*8 fb,f
      real*8 xmin,xmax
      real*8 dilphsp

c     hard thermal loop
      real*8 htlcor
      external htlcor

c     common variables
      real*8 am,T,muq
      logical htl_corr
      common /th_dyn/T,muq
      common /dilmass/am
      common /htl/htl_corr
      
      real*8 m_lept
      common /leptype/m_lept 
      
c-----|----------------------------------------------------------------|X

      if(T.eq.0.) then
         distpqgp=0.
         return
      endif

c     dilepton phase-space
      dilphsp=sqrt(1-4.d0*m_lept**2/am**2)*(1+2.d0*m_lept**2/am**2)     
      if(1-4.d0*m_lept**2/am**2.lt.0) return      
      
      q0=sqrt(q*q+am*am)
      fb=1.d0/(exp(q0/T)-1.d0)

      xmin=exp(-(q0-q)/2.d0/T)
      xmax=exp(-(q0+q)/2.d0/T)

      f=(xmin+exp(-(q0+muq)/T))/(xmax+exp(-(q0+muq)/T))*
     &  (xmax+exp(-muq/T))/(xmin+exp(-muq/T))
     
      f=log(f)*fb*T/q*dilphsp

c     transformation to invariant mass distributions
c
c     dN/(d^4x d^3p dM)=M/(p0) dR/d^4q(q0=pmunu)

c      if(htl_corr) distpqgp=f*am/p0*htlcor(q0,T) 
c      if(.NOT.htl_corr)  distpqgp=f*am/p0
 
      distpqgp=f*am/p0

      return
      end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c..  FUNCTION HTLCOR
c..
c..  This function calculates the hard thermal loop (HTL) correction for
c..  the QGP dilepton rates
c..
c..  (This subroutine was developed by Hendrik van Hees)
c..
c*****|****************************************************************|X

      real*8 FUNCTION htlcor(q0,temarg)
      IMPLICIT NONE

c-----|----------------------------------------------------------------|X

      real*8 q0,temarg,temp,alphsT,gsT,emqth,x
      
      real*8 pi  
      parameter (pi=3.141592653589793)

c-----|----------------------------------------------------------------|X
      temp=temarg*1.0d3 ! temp is in MeV

      alphsT=6.0d0*pi/(28.0d0*log(temp/22.0d0))
      gsT   =sqrt(alphsT*4.0d0*pi)
      emqth =gsT*temarg/sqrt(6.0d0)
      x=q0/emqth
      htlcor=0.5263d0*
     $     (1.9d0*(1.0d0-exp(-x**(1.28d0)/7.7d0))+
     $     1.6d0/x+50.0d0/x**2+60.0d0/x**3)
      end

c-----|----------------------------------------------------------------|X
c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$
c$$$c	subroutine outqgp(dx,dtime,xstep,ystep,it_step)
c$$$        subroutine outqgp(dx,dtime,xstep,ystep)
c$$$c       This routine writes the output 
c$$$c       The rate R [fm^-4] stored in "weightqgp" is multiplied
c$$$c       by the cell volume and the time-step
c$$$	implicit none
c$$$
c$$$	integer i,j
c$$$	real*8 sum
c$$$        real*8 dx,dtime
c$$$        integer xstep,ystep !,it_step
c$$$        real*8 vol4
c$$$        real*8 summ
c$$$
c$$$	include 'comdil.f'
c$$$        include 'comdilout.f'
c$$$
c$$$	sum=0.d0
c$$$c     4-volume in fm^4
c$$$c        vol4=dx**3*dtime*float(xstep*ystep*it_step)
c$$$
c$$$c     since it_step can assume different values according to the application 
c$$$c     or not of course graining in time for the emission, 
c$$$c     multiplication for it_step already performed in MC_dilqgp,
c$$$c     
c$$$         vol4=dx**3*dtime*float(xstep*ystep)
c$$$	do i=0,posmaxm
c$$$           summ=0d0
c$$$           do j=0,posmaxpt
c$$$	      if(weightqgp(i,j).gt.0d0)
c$$$     &         write(44,'(2f12.7,e14.7)') mmin+0.5d0*dm+i*dm,
c$$$     &             0.5d0*dpt+j*dpt,weightqgp(i,j)*vol4
c$$$              sum=sum+weightqgp(i,j)*vol4
c$$$              summ=summ+weightqgp(i,j)*vol4
c$$$           enddo
c$$$c           write(44,'(f12.7,e14.7)') mmin+0.5d0*dm+i*dm,
c$$$c     &      summ
c$$$	enddo
c$$$	    
c$$$        write(*,*)'#sumqgp',sum
c$$$      
c$$$        end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
c..                                                             !  
c..     **** INTEGRATION ROUTINES ****                          !
c..                                                             !
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!

c///////////////////////////////////////////////////////////////!

c234567**1*********2*********3*********4*********5*********6*********7**X
c..  SUBROUTINE QROMB
c..
c..  Performs a Romberg integration
c..
c*****|****************************************************************|X

      SUBROUTINE qromb(func,a,b,ss)

c-----|----------------------------------------------------------------|X

      INTEGER JMAX,JMAXP,K,KM
      REAL*8 a,b,func,ss,EPS
      EXTERNAL func
      PARAMETER (EPS=1.e-6, JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
CU    USES polint,trapzd
      INTEGER j
      REAL*8 dss,h(JMAXP),s(JMAXP)

c-----|----------------------------------------------------------------|X

      h(1)=1.
      do 11 j=1,JMAX
c        write(0,*)'in qromb:a,b,ss',a,b,ss ! Debug only
c        write(0,*)'func(a),func(b)',func(a),func(b) !Debug only
        call trapzd(func,a,b,s(j),j)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=0.25*h(j)
11    continue
      stop 'too many steps in qromb'
      END

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c..  SUBROUTINE MIDINF
c..
c..  Performs an integration of improper integral 
c..
c*****|****************************************************************|X     
         
      SUBROUTINE midinf(funk,aa,bb,s,n)

c-----|----------------------------------------------------------------|X

      INTEGER n
      REAL*8 aa,bb,s,funk
      EXTERNAL funk
      INTEGER it,j
      REAL*8 a,b,ddel,del,sum,tnm,func,x

c-----|----------------------------------------------------------------|X

c     hier change of variables
      func(x)=funk(1./x-1)/x**2
      b=1.d0/(1.d0+aa)
      a=1.d0/(1.d0+bb)
      if (n.eq.1) then
        s=(b-a)*func(0.5d0*(a+b))
      else
        it=3**(n-2)
        tnm=it
        del=(b-a)/(3.d0*tnm)
        ddel=del+del
        x=a+0.5d0*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+ddel
          sum=sum+func(x)
          x=x+del
11      continue
        s=(s+(b-a)*sum/tnm)/3.d0
      endif
      return
      END

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c..  SUBROUTINE LORSPHERE
c..
c*****|****************************************************************|X

      SUBROUTINE lorsphere(r,g,vx,vy,vz,m,pxmax,pxmin,
     &     pymax,pymin,pzmax,pzmin)
      implicit none

c-----|----------------------------------------------------------------|X

      real*8 r,g,vx,vy,vz,m,pxmax,pxmin,pymax,pymin,pzmax,pzmin
      real*8 p1x,p2x,p1y,p2y,p1z,p2z
      real*8 bqx,bqy,bqz
      real*8 pcx,pcy,pcz
      real*8 k,v2,qxbar,qybar,qzbar

c-----|----------------------------------------------------------------|X      

      v2=vx**2+vy**2+vz**2
      
c     boost point (qx,qy,qz)=(0,0,r)
      
      p1x=bqx(0.d0,0.d0,r,g,vx,vy,vz,m)
      p1y=bqy(0.d0,0.d0,r,g,vx,vy,vz,m)
      p1z=bqz(0.d0,0.d0,r,g,vx,vy,vz,m)
      
c     boost point (qx,qy,qz)=(0,0,-r)
      p2x=bqx(0.d0,0.d0,-r,g,vx,vy,vz,m)
      p2y=bqy(0.d0,0.d0,-r,g,vx,vy,vz,m)
      p2z=bqz(0.d0,0.d0,-r,g,vx,vy,vz,m)
      
c     write(*,*)'p1z,p2z',p1z,p2z
      
      pzmin=p2z
      pzmax=p1z
      
c     find if bqz has an intrinsic local minimum
c     or maximum along the (0,0,qz) axis
      k=g*vz*v2/(v2+vz**2*(g-1.d0))
      if(abs(k).gt.1) then
         if(vz.lt.0) then
            qzbar=m/sqrt(k**2-1)
            if(qzbar.lt.r) then
               p1z=bqz(0.d0,0.d0,qzbar,g,vx,vy,vz,m)
               pzmax=p1z
            endif
         elseif(vz.gt.0) then
            qzbar=-m/sqrt(k**2-1)
            if(qzbar.gt.-r) then
               p2z=bqz(0.d0,0.d0,qzbar,g,vx,vy,vz,m)
               pzmin=p2z
            endif
         endif
      endif
      
c     boost point (qx,qy,qz)=(0,r,0)
      p1y=bqy(0.d0,r,0.d0,g,vx,vy,vz,m)


c     boost point (qx,qy,qz)=(0,-r,0)
      p2y=bqy(0.d0,-r,0.d0,g,vx,vy,vz,m)
      
      
      pymin=p2y                 !min(p2y,pymin)
      pymax=p1y                 !max(p1y,pymax)
      
c     find if bqy has an intrinsic local minimum
c     or maximum along the (0,qy,0) axis
      k=g*vy*v2/(v2+vy**2*(g-1.d0))
      if(abs(k).gt.1) then
         if(vy.lt.0) then
            qybar=m/sqrt(k**2-1)
            if(qybar.lt.r) then
               p1y=bqy(0.d0,qybar,0.d0,g,vx,vy,vz,m)
               pymax=max(pymax,p1y)
            endif
         elseif(vy.gt.0) then
            qybar=-m/sqrt(k**2-1)
            if(qybar.gt.-r) then
               p2y=bqy(0.d0,qybar,0.d0,g,vx,vy,vz,m)
               pymin=min(pymin,p2y)
            endif
         endif
      endif
      
c     boost point (qx,qy,qz)=(r,0,0)
      
      p1x=bqx(r,0.d0,0.d0,g,vx,vy,vz,m)
c     !      p1y=bqy(r,0.d0,0.d0,g,vx,vy,vz,m)
c     !      p1z=bqz(r,0.d0,0.d0,g,vx,vy,vz,m)
      
c     boost point (qx,qy,qz)=(-r,0,0)
      p2x=bqx(-r,0.d0,0.d0,g,vx,vy,vz,m)
c     !      p2y=bqy(-r,0.d0,0.d0,g,vx,vy,vz,m)
c     !      p2z=bqz(-r,0.d0,0.d0,g,vx,vy,vz,m)
      
      
      pxmin=p2x                 !min(p2x,pxmin)
c     !      pymin=min(min(p1y,p2y),pymin)
c     !      pzmin=min(min(p1z,p2z),pzmin)
      pxmax=p1x                 !max(p1x,pxmax)
c     !      pymax=max(max(p1y,p2y),pymax)
c     !      pzmax=max(max(p1z,p2z),pzmax)
      
c     find if bqx has an intrinsic local minimum
c     or maximum along the (qx,0,0) axis
      k=g*vx*v2/(v2+vx**2*(g-1.d0))
      if(abs(k).gt.1) then
         if(vx.lt.0) then
            qxbar=m/sqrt(k**2-1)
            if(qxbar.lt.r) then
               p1x=bqx(qxbar,0.d0,0.d0,g,vx,vy,vz,m)
               pxmax=max(pxmax,p1x)
            endif
         elseif(vx.gt.0) then
            qxbar=-m/sqrt(k**2-1)
            if(qxbar.gt.-r) then
               p2x=bqx(qxbar,0.d0,0.d0,g,vx,vy,vz,m)
               pxmin=min(pxmin,p2x)
            endif
         endif
      endif
      
      
      return
      end
      
c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c..  FUNCTION FP
c..
c*****|****************************************************************|X

       real*8 FUNCTION fp(p0,px,py,pz,vx,vy,vz,g,a,T,mu)
       implicit none

c-----|----------------------------------------------------------------|X

       real*8 p0,px,py,pz,vx,vy,vz,g,T,mu
       integer a
       real*8 pmum

c-----|----------------------------------------------------------------|X  
    
       pmum=g*(p0-px*vx-py*vy-pz*vz)
       fp=1./(exp((pmum-mu)/T)+a)
       
       return
       end

c-----|----------------------------------------------------------------|X
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     returns boosted qx
      real*8 function bqx(qx,qy,qz,g,vx,vy,vz,m)
      implicit none
      real*8 qx,qy,qz,g,vx,vy,vz,m
      real*8 q0,v2,vsq
      
      q0=sqrt(m**2+qx**2+qy**2+qz**2)
      v2=vx**2+vy**2+vz**2
      vsq=vx*qx+vy*qy+vz*qz

      bqx=qx+(g-1.d0)/v2*vsq*vx+g*vx*q0
      
      
      return
      end
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       returns boosted qy
      real*8 function bqy(qx,qy,qz,g,vx,vy,vz,m)
      implicit none
      real*8 qx,qy,qz,g,vx,vy,vz,m
      real*8 q0,v2,vsq

      q0=sqrt(m**2+qx**2+qy**2+qz**2)
      v2=vx**2+vy**2+vz**2
      vsq=vx*qx+vy*qy+vz*qz

      bqy=qy+(g-1.d0)/v2*vsq*vy+g*vy*q0
      
      return
      end
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       returns boosted qz
      real*8 function bqz(qx,qy,qz,g,vx,vy,vz,m)
      implicit none
      real*8 qx,qy,qz,g,vx,vy,vz,m
      real*8 q0,v2,vsq
      
      q0=sqrt(m**2+qx**2+qy**2+qz**2)
      v2=vx**2+vy**2+vz**2
      vsq=vx*qx+vy*qy+vz*qz
      
      bqz=qz+(g-1.d0)/v2*vsq*vz+g*vz*q0
      
      return
      end
      
