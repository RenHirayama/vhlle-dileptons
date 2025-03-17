c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE DILEMIT
c..
c..   This program returns in-medium resp. vacuum dilepton rates dR/dM 
c..   (up to known constant factors). It uses tables of the kind
c..
c..               mu, T, M, dR/dM
c..
c..   where mu: chemical potential [GeV]
c..          T: temperature [GeV]
c..          M: invariant mass [GeV]
c..      dR/dM: in-medium dilepton rate [a.u.]
c..   
c..   To obtain the rate in physical units the following
c..   factor has to be multiplied:
c..   
c..   factor=alpha_em**2*4.d0/pi**2*mrho**4/grho**2
c..
c..   to obtain the rate in units of [GeV^-2*fm^-4]
c..   one has to further divide for (hqc**4).
c..
c..   In summary:
c..   -----------------------------------------------
c..   |  dR/dM[phys.u.]=dR/dM[a.u.]*factor/(hqc**4) |
c..   -----------------------------------------------
c..
c..   To obtain dR/dM^2 (often plotted in papers) 
c..   one has to further divide for 2M,i.e.     
c..   dR/dM^2[phys.u.]=dR/dM[phys.u.]/2./M
c..
c*****|****************************************************************|X

      SUBROUTINE dilemit(vac,temp,mub,pipot,gce,vxce,vyce,vzce,vol4,multi,
     &                   beta_lab,dt,timestep,lambda) 
      implicit none

      include 'defs.f'

c-----|----------------------------------------------------------------|X

      real*8 lambda,temp,mub,pipot,gce,vxce,vyce,vzce
      real*8 mass,p0l,pxl,pyl,pzl
      integer it_step
      real*8 contr,vol4
      real*8 factor,rate,result
      real*8 ammax,ammin  

      real*8 fugacity

      real*8 grho
      parameter (grho=5.03)
      real*8 mminit(1:200)    
      real*8 massit(1:200)

      integer multi,stat,l
      integer flagrho,vac
      parameter (flagrho=1)

      integer noe,ityp
      real*8 bev,acce,accp
      real*8 effvol4

      real*8 p0_po_lab,px_po_lab,py_po_lab,pz_po_lab
      real*8 p0_el_lab,px_el_lab,py_el_lab,pz_el_lab
      real*8 beta_lab

      real*8 dt,time
      integer timestep
      
c.. functions
      real*8 distmhad
      external distmhad

c. Common Blocks
      real*8 T,mu
      common /th_dyn/T,mu
      logical vacuum
      common /vflag/vacuum

c.. Steps to increase statistics
      logical zwstep
      real*8 zwmin(1:12),zwmax(1:12)
      integer j

c-----|----------------------------------------------------------------|

c.. Set common variables
      T=temp
      mu=mub      

c      write(0,*)'T,temp,mu,mub',T,temp,mu,mub ! Debug only
c      write(0,*)'vol4',vol4 ! Debug only   
   
      vacuum=.TRUE.
      if(vac.eq.0) vacuum=.FALSE.

c.. Parameters
c      write(0,*)'Dilepton emission' ! Debug only

      massit(104)=0.77d0
      mminit(104)=0.277d0

c.. Determine the mass of the lepton-type 
      if(dimuon) then
         m_lept=mass_muon
      else !dielectron
         m_lept=mass_electron
      endif

c      write(0,*)'lepton mass',m_lept ! Debug only 

c..Define masses for common block
      mmin=2d0*m_lept

c.. Threshold is fixed at 2*m_lept
      ammin=mmin

c.. For vacuum calculation threshold coincides with 
c.. the rho MASS THRESHOLD

c      write(0,*)'vac =',vac ! Debug only
      if(vac.eq.1) ammin=mminit(104)

c.. MAXIMUM MASS for integration
      ammax=mmax-0.00001d0

c.. MASS STEPS
      zwmin(1)=mmin
      zwmax(1)=0.05d0
      zwmin(2)=0.05001d0
      zwmax(2)=0.1d0
      zwmin(3)=0.10001d0
      zwmax(3)=0.15d0
      zwmin(4)=0.15001d0
      zwmax(4)=0.2d0
      zwmin(5)=0.20001d0
      zwmax(5)=0.3d0
      zwmin(6)=0.30001d0
      zwmax(6)=0.4d0
      zwmin(7)=0.40001d0
      zwmax(7)=0.5d0
      zwmin(8)=0.50001d0
      zwmax(8)=0.6d0
      zwmin(9)=0.60001d0
      zwmax(9)=0.7d0
      zwmin(10)=0.70001d0
      zwmax(10)=0.8d0
      zwmin(11)=0.80001d0
      zwmax(11)=1.0d0
      zwmin(12)=1.00001d0
      zwmax(12)=mmax-0.00001d0

c.. Time

      time=dt*dble(timestep)

c-----|----------------------------------------------------------------|X
c. **** If: zwstep=.TRUE. : Loop to increase statistics at high masses ****
c*********************!   
      zwstep=.true.   !
c     zwstep=.false.  !
c*********************!

c.. Set counter to zero
      j=0

 123  continue

      if(zwstep.AND.(.NOT.dimuon).AND.vac.ne.1) then
        j=j+1
        ammin=zwmin(j)
        ammax=zwmax(j)
      endif

c-----|----------------------------------------------------------------|X

c     calculate mass and momentum integrated rate
c     "rate" has units [fm^{-4}]

c. **** Integration over mass dependent function - SIMPSON INTEGRATION ****
      rate=0d0

c      write(0,*)'Simpson Integration' ! Debug only 
c      write(0,*)'ammin,ammax,result',ammin,ammax,result ! Debug only
      CALL qsimp_had(distmhad,ammin,ammax,result)
c      write(0,*)'ammin,ammax,result',ammin,ammax,result ! Debug only
c      write(0,*)'Integration finished' ! Debug only 
c     constant factor
           factor=alpha_em**2*4.d0/pi**2*massit(104)**4/grho**2
           rate=result*factor
c     conversion to fm^{-4}
      rate=rate/(hqc**4)

c      write(0,*) 'rate',rate ! Debug only

c-----|----------------------------------------------------------------|X
c. **** Loop n-times for better mass and momentum statsitics ****

c.. Number of Loops
cccccccccccccccccccccccccccc!
      stat=1                !                            
cccccccccccccccccccccccccccc! 

      do 678 l=1,stat

c. **** Generate mass ****
      mass=0.0d0
c      write(0,*)'mass',mass ! Debug only
      CALL massdist(vac,mass,ammin,ammax)
c      write(0,*)'mass',mass ! Debug only

c.. in case of too many iterations 
c.. mass is set to zero in massdist
      if(mass.eq.0) return

c. **** Generate 3-momenta ****   
c      write(0,*)'vac,mass,gce,vxce,vyce,vzce',vac,mass,gce,vxce,vyce,
c     &                                        vzce ! Debug only 
      CALL momdist(vac,mass,gce,vxce,vyce,vzce,p0l,pxl,pyl,pzl)

c.. In case of too many interations
c.. p0l is set to zero in momdist
 
c      write(0,*)'p0l,pxl,pyl,pzl',p0l,pxl,pyl,pzl ! Debug only 
      if(p0l.eq.0d0) then
       write(0,*)'p0l = 0, return'
       return
      endif            
 
c.. Determine momenta of e+ and e-

      CALL pel_ppo_dist(beta_lab,mass,p0l,pxl,pyl,pzl,p0_el_lab,
     &px_el_lab,py_el_lab,pz_el_lab,p0_po_lab,px_po_lab,py_po_lab,
     &pz_po_lab)     

c      write(0,*)'p0e,pxe,pye,pze',p0_el_lab,px_el_lab,py_el_lab,
c     &                            pz_el_lab ! Debug only 
c      write(0,*)'p0e,pxe,pye,pze',p0_po_lab,px_po_lab,py_po_lab,
c     &                            pz_po_lab ! Debug only 

c.. Additional FUGACITY factor
      fugacity=dexp(2.0d0*pipot/T)

c.. Determine cell contribution
      contr=rate*vol4*(1.0d0-lambda)*dble(multi)/dble(stat)*fugacity

c. **** Write into output file f71 ****
      ityp=104
      noe=1
      bev=0.0d0
      acce=1.0d0
      accp=1.0d0
      effvol4=vol4*(1.0d0-lambda)*dble(multi)/dble(stat)
c........................................................................
c. NOTE: The extended output format writes out ***lab-system momenta*** !
c.       for e+ and e-                                                  !
c.......................................................................!
      if(ext_out) then !extended output format
       write(71,556)ityp,contr,mass,p0_el_lab,px_el_lab,py_el_lab,
     & pz_el_lab,p0_po_lab,px_po_lab,py_po_lab,pz_po_lab,dt,time,effvol4,
     & flagrho,mub,temp,noe,lambda,acce,accp
      else !standard output format
       write(71,555)contr,mass,pxl,pyl,pzl,flagrho,mub,temp
      endif

 555  format(e14.7,4f12.7,i9,2f12.7)
 556  format(I5,1X,13(E16.8E3,1X),I4,1X,2(E16.8E3,1X),I9,1X,
     &F12.8,1X,2(E16.8,1x))

 678  continue

c-----|----------------------------------------------------------------|X

      if(zwstep.AND.(.NOT.dimuon).AND.vac.ne.1.AND.j.lt.12) goto 123

      return
      end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE MASSDIST
c..
c.. This subroutine generates the dilepton mass for rho-meson emission
c*****|****************************************************************|X

      subroutine massdist(vac,m,ammin,ammax)
      implicit none

      include 'defs.f'

c-----|----------------------------------------------------------------|X

      real*8 ammax,ammin,h,distmhad,hmax,bessk1
      real*8 mminit(1:200)
      real*8 m,mtest,htest
      real*8 ranff
      external ranff
      integer j,b,vac
      real*8 drdmmax
      real*8 aux
c      real*8 rand

c. Common Blocks
      real*8 T,mu
      common /th_dyn/T,mu
  
c-----|----------------------------------------------------------------|X
    
      mminit(104)=0.277d0

c      ammin=mmin
c.. For vacuum calculation threshold coincides with 
c.. the rho mass threshold
c      if(vac.eq.1) ammin=mminit(104)
c      ammax=mmax
c      write(0,*)'ammin,ammax',ammin,ammax ! Debug only
c...VACUUM CASE
      if(vac.eq.1) then
c         hmax=bessk1(ammin/T)*9.d0*T ![GeV^{-2}]
         hmax=0.d0
         do b=1,10
           hmax=hmax+(1.0d0/b)*bessk1(b*(ammin/T))
         end do 
         hmax=hmax*9.d0*T       ![GeV^{-2}]
c.	added 01.05.2010
c         write(0,*)'in massdist: hmax',hmax !Debug only
         hmax=min(hmax,drdmmax(mu,T))
c...IN-MEDIUM CASE
      else
c      maximum of the distribution 
       mtest=ammin-0.0025d0
       if(mtest.le.0.0d0) mtest=0.0001d0
       hmax=0.0d0

362    continue
       htest=distmhad(mtest)
       if(htest.gt.hmax) hmax=htest
       if(mtest.lt.(ammax).AND.mtest.lt.1.485d0) then
         mtest=mtest+0.0025d0 
         goto 362
       endif
      endif

       hmax=1.2*hmax
c      hmax=drdmmax(mu,T)      

c     write(0,*) 'hmax',hmax ! Debug only

      j=0
 2    continue
      j=j+1

c      write(0,*)'in massdist: ammin,ammax',ammin,ammax ! Debug only
c      m=ammin+dble(rand(0))*(ammax-ammin)
      m=ammin+(ranff(seedinit))*(ammax-ammin)
c      write(0,*) 'generated mass:', m ! Debug only
c      h=dble(rand())*hmax
      h=ranff(seedinit)*hmax
      
      if(distmhad(m).gt.hmax) then
         aux=distmhad(m)
         write(0,*)'hmax too small'
         write(0,*)'hmax',hmax
         write(0,*)'distmhad(m)',aux
         write(0,*)'maximum violated by', (aux-hmax)/hmax
         if((aux-hmax)/hmax.gt.5.d-2) stop
c      stop
      endif
      
c     mass could not be generated
      if(j.gt.10000000) then
         m=0.d0
         write(0,*)'too many trials to generate mass'
         write(0,*) 'number of trials',j
         return
      endif
      
      if(h.gt.distmhad(m)) goto 2
      
      end


c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE MOMDIST
c..
c..  This subroutine generates the dilepton momenta for rho-meson
c..  emission
c*****|****************************************************************|X

      SUBROUTINE momdist(vac,m,g,vx,vy,vz,p0,px,py,pz) 
      implicit none

c-----|----------------------------------------------------------------|X

      real*8 m,g,vx,vy,vz,p0,px,py,pz
      integer i,vac
c      real*8 ranf
      real*8 ranff
c      external ranff
      real*8 ftrial,massit
      real*8 drdmd3qmax
      real*8 q0,qx,qy,qz,q,phi,cost,sint,v2
      real*8 xmin,xmax,x,x2,c,Fint,fb,compf,truef
      real*8 imDrho
      real*8 pi
      parameter (pi=3.1415926535) ! useful constant
      real*8 s9,s10,eps
      parameter(eps=1.d-6)
c      real*8 rand
      real*8 test1, test2

c     common variables
      real*8 qup,qdown,q0up,q0down
      common /qupdow/ qup,qdown
      real*8 T,mub
      common /th_dyn/T,mub
      integer seedinit
      common /randomf/ seedinit

c-----|----------------------------------------------------------------|X

c.. maximum energy allowed from tabulated momenta
       q0up=sqrt(qup**2+m**2)
c       write(0,*) 'q0up',q0up ! Debug only
c..   minimum energy
       q0down=m

      if(vac.eq.1) then
         ftrial=1.d0
      else
         ftrial=drdmd3qmax(mub,T,m)
      endif
c      write(0,*)'ftrial',ftrial ! Debug only
      xmin=0.d0
      c=Fint(m,q0down,T)
      xmax=Fint(m,q0up,T)-c
c      write(0,*)'xmax',xmax

c..generate random phi
c       phi=ranf(0)*2.d0*pi
       phi=ranff(seedinit)*2.d0*pi
c       phi=dble(rand(1))*2.d0*pi
c..generate random cos($\theta$)
c       cost=-1d0+ranf(0)*(1d0+1.d0)
       cost=-1d0+ranff(seedinit)*(1d0+1.d0)
c       cost=-1d0+dble(rand())*(1d0+1.d0)
       
       i=0
 1     continue
       i=i+1
c..generate q0 according to the distribution
c     1/-1+exp(-q0/T)
c..the comparison functions is qup*imDrhomax
c       x=ranf(0)*xmax
       x=ranff(seedinit)*xmax
c       x=dble(rand())*xmax
       
       if(m/T.lt.10d0) then
             q0=-T*Log(1d0-Exp(x/T)*(1.d0-Exp(-m/T)))
c             write(0,*)'T,x,q0',T,x,q0 ! Debug only
             if(1d0-Exp(x/T)*(1.d0-Exp(-m/T)).le.0d0) then
                write(0,*) 'numerical rumor'
                write(0,*) 'extracted log of non positive number'
                write(0,*) 'logarg',1d0-Exp(x/T)*(1.d0-Exp(-m/T))
                write(0,*) 'x,T,m',x,T,m
                stop
             endif
       else
c        Boltzmann approximation for m/T>10
          q0=-(T*Log((T*Exp(-m/T)-x)/T))
c          write(0,*)'T,x,Boltzm. approx. for q0',T,x,q0 ! Debug only
          if((T*Exp(-m/T)-x)/T.le.0d0) then
             write(0,*) 'numerical rumor in Bolt. appr.'
             write(0,*) 'extracted log of non positive number'
             write(0,*) 'logarg',(T*Exp(-m/T)-x)/T
             write(0,*) 'x,T,m',x,T,m
             stop
          endif
       endif

c.. In principle(mathematically) the obtained q0 is 
c.. always greater m, but numerics plays bad games
c.. when q0 is about m. This is the reason of this
c.. conditional line
       if(q0**2.lt.m**2) q0=m
       q=sqrt(q0**2-m**2)
c       write(0,*)'q =',q ! Debug only 
c.. boson distribution function
       fb=1.d0/(exp(q0/T)-1.d0)
c       write(0,*)'boson distr. func.',fb ! Debug only 
       compf=qup*ftrial*fb
c       write(0,*)'compf =',compf ! Debug only 
       if(vac.eq.1) then
           truef=q*fb
c           write(0,*)'in vacuum truef=',truef ! Debug only
       else
          if(q.ge.qdown) then
             truef=q*fb*imDrho(mub,T,q,m)
c             write(0,*)'imDrho(mub,T,q,m)',imDrho(mub,T,q,m) ! Debug only 
c             write(0,*)'in-medium truef=',truef ! Debug only 
          else
c.. for q<qdown imDrho is approx with its value at qdown
             truef=q*fb*imDrho(mub,T,qdown,m)
c             write(0,*)'in-medium (q<qdown) truef=',truef ! Debug only 
          endif  
       endif
c.. added for stability
c       write(0,*) 'truef',truef ! Debug only
       if(truef.le.0.d0) then
         p0=0.d0
       return
       endif
 
       if(compf.lt.truef) then
          write(0,*)'comparison function too small'
          write(0,*)'T,mu',T,mub
          write(0,*)'m,q,q0',m,q,q0
          write(0,*)'qup,qdown,q0up,q0down',qup,qdown,q0up,q0down
          write(0,*)'ftrial,fb',ftrial,fb
          write(0,*)'imDrho,drdmd3qmax',imDrho(mub,T,q,m),drdmd3qmax(mub,T,m)
          write(0,*)'compf,truef',compf,truef
          write(0,*) 'ERROR'
          stop
       endif
c       x2=ranf(0)
       x2=ranff(seedinit)
c       x2=dble(rand())
c       if(x2.gt.truef/compf) goto 1 
       if(x2.gt.truef/compf.and.i.lt.1000000) then
          goto 1
       endif
       if(i.ge.1000000) then
         write(0,*) 'too many iteration in had_mom' 
         p0=0.d0
          return
       endif

c.. Find cartesian coordinates
       sint=sqrt(1.d0-cost**2)
       qx=q*cos(phi)*sint
       qy=q*sin(phi)*sint
       qz=q*cost
c.. Boost to cf
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
c       write(0,*) 'number of trials in had_momdist:', i ! Debug only

c.. Check for correct values of momenta

c       test1=sqrt(m**2+px**2+py**2+pz**2)
c       write(0,*)'CHECK:Are momenta correct?',test1,p0 ! Debug only
c       test2=1.0d0/sqrt(1.0d0-v2)-g
c       write(0,*)'g comparison =',test2 ! Debug only
c       p0=test1

       return 
       end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. FUNCTION FINT
c..
c*****|****************************************************************|X

      real*8 FUNCTION Fint(m,q0,temp)
      implicit none

c-----|----------------------------------------------------------------|X

      real*8 m,temp,q0

c-----|----------------------------------------------------------------|X
      
      if(m/temp.lt.10d0) then
         Fint=temp*log(-1d0+exp(q0/temp))-q0
      else
c.. Boltzmann approximation for m/T>10
	Fint=-temp*exp(-q0/temp)
      endif

      return
      end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. FUNCTION DRDMD3QMAX
c..
c..  Given values of (mub,temp,mass) defing a point P in the
c..  mub-temp-mass space, this routine returns the
c..  maximum of the function imDrho within the 2x2x2 cube 
c..  centered around the point P
c*****|****************************************************************|X

      real*8 FUNCTION drdmd3qmax(mub,temp,mass)
      implicit none
  
      include 'defs.f'

c-----|----------------------------------------------------------------|X

      real*8 mass
      real*8 temp,mub
      integer m
      parameter(m=2)
      
      integer i,j,k,ii,jj,kk
c      real*8 dmu
c      parameter (dmu=1.d-2)

      real*8 masstab(kmmax),dmd3qmax(imumax,jtmax,kmmax)
      real*8 masstabmin,masstabmax
      common /vmmass/ masstab,masstabmax,masstabmin
      common /maxdmd3q/dmd3qmax
      
      real*8 mutab(imumax),ttab(jtmax)
      real*8 mutabmax,mutabmin,ttabmax,ttabmin
      common /tbvarel/mutab,ttab,mutabmax,mutabmin,
     &     ttabmax,ttabmin

c..   auxiliary variables
      real*8 maxaux(imumax),dy,res

c-----|----------------------------------------------------------------|X

c      write(0,*)'temp,mub',temp,mub ! Debug only

      if(mub.gt.mutabmax) then
         write(0,*) 'chemical potential too high! need bigger table'
         mub=mutabmax
      endif

c..   determine i (chemical potential index)
      i=int((mub-mutabmin)/dmu)+1
      ii=min(max(i-(m-1)/2,1),imumax+1-m)
      
      if(mass.gt.masstabmax) then
          write(0,*) 'invariant mass too high! need bigger mass table'
         stop
      endif

c..   determine k (mass index)
      k=int((mass-masstabmin)/dm)+1
      kk=min(max(k-(m-1)/2,1),kmmax+1-m)

c      write(0,*) 'i,ii,k,kk',i,ii,k,kk ! Debug only

      if(temp.gt.ttabmax) then
         write(0,*) 'from drdmd3qmax'
         write(0,*) 'temperature too high! need bigger table'
         write(0,*) 'temp,ttabmax',temp,ttabmax
         temp=ttabmax
      endif

c..   determine j (temperature index)
      j=int((temp-ttabmin)/(ttabmax-ttabmin)
     &     *dble(jtmax-1))+1
      jj=min(max(j-(m-1)/2,1),jtmax+1-m)
      
c      write(0,*) 'j,jj',j,jj

      res=0.d0
      do i=ii,ii+m-1
c     find maximum on the  (m,m) matrix
c     (1,1)=(jj,kk)  (1,2)=(jj,kk+1) ... (1,m)=(jj,kk+m-1)
c     (2,1)=(jj+1,kk)  ...
c     ...
c     (m,1)=(jj+m-1,kk)   ....            (m,m)=(jj+m-1,kk+m-1)
         maxaux(i)=0.d0
            do j=1,m
               do k=1,m
                  maxaux(i)=max(maxaux(i),dmd3qmax(i,jj+j-1,kk+k-1))
               enddo
            enddo
c     find maximum on the mu grid
            res=max(res,maxaux(i))
      enddo
      

      drdmd3qmax=res      
c     write(0,*)'res',res      

      return
      end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. FUNCTION IMDRHO
c..
c..   Imaginary part of the in-medium rho-meson propagator
c*****|****************************************************************|X

      real*8 FUNCTION imDrho(mub,temp,q,mass)!(M,prho)
      implicit none

      include 'defs.f'

c-----|----------------------------------------------------------------|X

      real*8 mub,temp,q,mass
      complex*16 selfaux,self
      real*8 reself,imself
      real*8 dilphsp
      real*8 massit(1:200)

c-----|----------------------------------------------------------------|X
 
c.. Determine the mass of the lepton-type 
      if(dimuon) then
         m_lept=mass_muon
      else !dielectron
         m_lept=mass_electron
      endif
   
c      write(0,*)'in imDrho: mub,temp,q,mass',mub,temp,q,mass ! Debug only

      massit(104)=0.77d0

      imDrho=0.d0

c.. Dilepton phase-space
      dilphsp=sqrt(1.d0-4.d0*m_lept**2/mass**2)*
     &   (1.d0+2.d0*m_lept**2/mass**2)

      if(1.d0-4.d0*m_lept**2/mass**2.lt.0.d0) return

      selfaux=self(mub,temp,q,mass)
      reself=dreal(selfaux)
      imself=dimag(selfaux)
c      write(0,*)'selfaux,reself,imself',selfaux,reself,imself ! Debug only
      imDrho=-imself/((mass**2-massit(104)**2-reself)**2+imself**2)

      return
      end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. FUNCTION DISTMHAD
c..
c..  dR/dM distribution function (without constant factors)
c*****|****************************************************************|X

      real*8 FUNCTION distmhad(mass)
      implicit none

      include 'defs.f'

c-----|----------------------------------------------------------------|X

      real*8 mass
      real*8 temp,mub
      integer m
      parameter(m=4)
      
      integer i,j,k,ii,jj,kk
c      real*8 dmu
c      parameter (dmu=1.d-2)

      integer b
      real*8  bosesum,bessk1
      complex*16 vacselfrho
      real*8 reself,imself,imDrho
      real*8 massit(1:200)

      real*8 masstab(kmmax),ratetab(imumax,jtmax,kmmax)
      real*8 masstabmin,masstabmax
      common /vmmass/ masstab,masstabmax,masstabmin
      common /drdm/ratetab
      real*8 mutab(imumax),ttab(jtmax)
      real*8 mutabmax,mutabmin,ttabmax,ttabmin
      common /tbvarel/mutab,ttab,mutabmax,mutabmin,
     &     ttabmax,ttabmin

c.. Auxiliary variables
      real*8 rateaux(m,m),yre(imumax),dy,res

c. Common Blocks
      real*8 T,mu
      common /th_dyn/T,mu
      logical vacuum
      common /vflag/vacuum
     
c-----|----------------------------------------------------------------|X
      
      massit(104)=0.77d0

c. **** VACUUM CASE ****
      if(vacuum) then
         bosesum=0d0
         do b=1,10
            bosesum=bosesum+(1.0d0/b)*bessk1(b*(mass/T)) 
         enddo         
         reself=dreal(vacselfrho(mass))
         imself=dimag(vacselfrho(mass))
         imDrho=-imself/((mass**2-massit(104)**2-reself)**2+imself**2)
         distmhad=bosesum*T*imDrho
c         write(0,*)'distmhad,reself,imself,bosesum,T',distmhad,reself,
c     &              imself,bosesum,T ! Debug only
         return
      endif

c-----|----------------------------------------------------------------|X

c. **** IN-MEDIUM CASE ****
c.. Rename variables
      temp=T
      mub=mu

c      write(0,*)'in distmhad:temp,mub',temp,mub ! Debug only
      
      if(mub.gt.mutabmax) then
         write(0,*) 'chemical potential too high! need bigger table'
         mub=mutabmax
      endif

c.. Determine i (chemical potential index)
      i=int((mub-mutabmin)/dmu)+1
      ii=min(max(i-(m-1)/2,1),imumax+1-m)
      
      if(mass.gt.masstabmax) then
          write(0,*) 'from distmhad'
          write(0,*) 'invariant mass too high! need bigger mass table'
          write(0,*) 'mass,masstabmax',mass,masstabmax
         stop
      endif

c.. Determine k (mass index)
      k=int((mass-masstabmin)/dm)+1
      kk=min(max(k-(m-1)/2,1),kmmax+1-m)

      if(temp.gt.ttabmax) then
         write(0,*) 'from distmhad'
         write(0,*) 'temperature too high! need bigger table'
         write(0,*) 'temp,ttabmax',temp,ttabmax
         temp=ttabmax
      endif

c.. Determine j (temperature index)
      j=int((temp-ttabmin)/(ttabmax-ttabmin)
     &     *dble(jtmax-1))+1
      jj=min(max(j-(m-1)/2,1),jtmax+1-m)
c      write(0,*) 'ttabmin,ttabmax',ttabmin,ttabmax ! Debug only
c      write(0,*) 'j,jj',j,jj ! Debug only

      do i=ii,ii+m-1

c-----|----------------------------------------------------------------|X
c.. Map indeces on a (m,m) matrix:                                     |
c.. (1,1)=(jj,kk)  (1,2)=(jj,kk+1) ... (1,m)=(jj,kk+m-1))              |
c.. (2,1)=(jj+1,kk)  ...)                                              |
c..  ...                                                               |
c.. (m,1)=(jj+m-1,kk)   ....            (m,m)=(jj+m-1,kk+m-1)          |
c-----|----------------------------------------------------------------|X
     
       do j=1,m
        do k=1,m
           rateaux(j,k)=ratetab(i,jj+j-1,kk+k-1)   
        enddo
       enddo

c.. Interpolate the (T,m) grid
c      write(0,*)'CALL polint2' ! Debug only
      CALL polin2(ttab(jj),masstab(kk),rateaux,m,m,
     &        temp,mass,yre(i),dy)
c      write(0,*)'rateaux(m,m),m',rateaux(m,m),m !Debug only
      enddo

c.. Now interpolate results on the mu raw
c      write(0,*)'CALL polint' ! Debug only
c      write(0,*)'mutab(ii),yre(ii),res,dy',mutab(ii),yre(ii),res,dy ! Debug only
      CALL polint(mutab(ii),yre(ii),m,mub,res,dy)
c      write(0,*)'mutab(ii),yre(ii),res,dy',mutab(ii),yre(ii),res,dy ! Debug only
 
      distmhad=0.0d0
      distmhad=res      

c      write(0,*)'distmhad',distmhad ! Debug only
      res=0

      return
      end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. FUNCTION DRDMMAX
c..
c..  Given values of (mub,temp) defing a point P in the
c..  mub-temp plane, this routine returns the
c..  maximum of the invariant mass distribution dR/dM 
c..  (up to known constant factors)
c*****|****************************************************************|X

      real*8 function drdmmax(mub,temp)
      implicit none

      include 'defs.f'
      
c-----|----------------------------------------------------------------|X

      real*8 temp,mub
      integer m
      parameter(m=4)
      integer i,j,ii,jj
c      real*8 dmu
c      parameter (dmu=1.d-2)
     
c     common variables
c      integer imumax,jtmax
c      parameter (imumax=95,jtmax=101) 

      real*8 ratemax(imumax,jtmax)
      common/maxdr/ratemax
      real*8 mutab(imumax),ttab(jtmax)
      real*8 mutabmax,mutabmin,ttabmax,ttabmin
      common /tbvarel/mutab,ttab,mutabmax,mutabmin,
     &     ttabmax,ttabmin

c.. auxiliary variables
      real*8 maxaux(imumax),dy,res
           
c-----|----------------------------------------------------------------|X

       if(mub.gt.mutabmax) then
         write(0,*) 'chemical potential too high! need bigger table'
         mub=mutabmax
      endif

c.. determine i (chemical potential index)
      i=int((mub-mutabmin)/dmu)+1
      ii=min(max(i-(m-1)/2,1),imumax+1-m)
c      write(0,*)'drdmmx: i, ii ',i,ii ! Debug only

      if(temp.gt.ttabmax) then
         write(0,*) 'from drdmmax'
         write(0,*) 'temperature too high! need bigger table'
         write(0,*) 'temp,ttabmax',temp,ttabmax
         temp=ttabmax
      endif

c.. determine j (temperature index)
      j=int((temp-ttabmin)/(ttabmax-ttabmin)
     &     *dble(jtmax-1))+1
      jj=min(max(j-(m-1)/2,1),jtmax+1-m)
c      write(0,*)'j-(m-1)/2,jtmax+1-m',j-(m-1)/2,jtmax+1-m ! Debug only
c      write(0,*)'drdmmx: temp, ttabmin, ttabmax',temp,ttabmin,ttabmax ! Debug only      
c      write(0,*)'drdmmx: j, jj ',j,jj ! Debug only      

      res=0.d0
      do i=ii,ii+m-1
c.. find maximum on the T grid
         maxaux(i)=0.d0
         do j=jj,jj+m-1
            maxaux(i)=max(maxaux(i),ratemax(i,j))
         enddo
c.. find maximum on the mu grid
      res=max(res,maxaux(i))
      enddo
      
      drdmmax=res
c      write(0,*)'drdmmax',drdmmax ! Debug only

      end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. FUNCTION SELF
c..
c..  VM self energy
c*****|****************************************************************|X

      complex*16 FUNCTION self(mub,temp,q,mass)
      implicit none

      include 'defs.f'

c-----|----------------------------------------------------------------|X

      real*8 mub,temp,q,mass
      integer i,j,k,ii,jj,kk
      integer m
      parameter (m=4)
      real*8 dq
      parameter (dq=1.d-2)
      complex*16 vacselfrho

c     common variables
c      integer imumax,jtmax,kqmax
c      parameter (imumax=95,jtmax=101,kqmax=300)

      real*8 mutab(imumax),ttab(jtmax)
      real*8 qtab(kqmax),reself(imumax,jtmax,kqmax),
     &     imself(imumax,jtmax,kqmax)
      real*8 mutabmax,mutabmin,ttabmax,ttabmin
      real*8 qtabmax,qtabmin
      common /tbvarel/mutab,ttab,mutabmax,mutabmin,
     &     ttabmax,ttabmin
      common /eltbself/qtab,reself,imself,qtabmax,qtabmin

c..   auxiliary variables
      real*8 reaux(m,m),imaux(m,m),
     &     yre(imumax),yim(imumax),dy,res,ims

c-----|----------------------------------------------------------------|X

      if(mub.gt.mutabmax) then
         write(0,*) 'chemical potential too high! need bigger table'
         mub=mutabmax
      endif

c.. Determine i (chemical potential index)
      i=int((mub-mutabmin)/dmu)+1
      ii=min(max(i-(m-1)/2,1),imumax+1-m)


      if(q.gt.qtabmax) then
         write(0,*) 'momentum too high!need bigger table'
         stop
      endif
c.. Determine k (momentum index)
      k=int((q-qtabmin)/dq)+1
      kk=min(max(k-(m-1)/2,1),kqmax+1-m)

      if(temp.gt.ttabmax) then
         write(0,*) 'from self'
         write(0,*) 'temperature too high! need bigger table'
         write(0,*) 'temp,ttabmax',temp,ttabmax
         temp=ttabmax
      endif

c.. Determine j (temperature index)
      j=int((temp-ttabmin)/(ttabmax-ttabmin)
     &     *dble(jtmax-1))+1
      jj=min(max(j-(m-1)/2,1),jtmax+1-m)

      do i=ii,ii+m-1

c-----|----------------------------------------------------------------|X
c.. Map indeces on a (m,m) matrix:                                     |
c.. (1,1)=(jj,kk)  (1,2)=(jj,kk+1) ... (1,m)=(jj,kk+m-1)               |
c.. (2,1)=(jj+1,kk)  ...)                                              |
c..  ...                                                               |
c.. (m,1)=(jj+m-1,kk)   ....            (m,m)=(jj+m-1,kk+m-1)          |
c-----|----------------------------------------------------------------|X

         do j=1,m
            do k=1,m
               reaux(j,k)=reself(i,jj+j-1,kk+k-1)
               imaux(j,k)=imself(i,jj+j-1,kk+k-1)
            enddo
         enddo

c.. Interpolate the (T,q) grid
         call polin2(ttab(jj),qtab(kk),reaux,m,m,
     &        temp,q,yre(i),dy)
         call polin2(ttab(jj),qtab(kk),imaux,m,m,
     &        temp,q,yim(i),dy)
      enddo

c.. Now interpolate results on the mu raw
      CALL polint(mutab(ii),yre(ii),m,mub,res,dy)
      CALL polint(mutab(ii),yim(ii),m,mub,ims,dy)
 
c      write(0,*)'dcmplx,vacselfrho',dcmplx(res,ims),vacselfrho(mass) ! Debug only
      self=dcmplx(res,ims)+vacselfrho(mass)

      end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. FUNCTION VACSELFRHO
c..
c*****|****************************************************************|X

      complex*16 FUNCTION vacselfrho(M)
      implicit none

c-----|----------------------------------------------------------------|X

      real*8 M
      real*8 reself,imself,grho2,pi,ome0,p0
      real*8 mrho,mpi
      real*8 vacwidth
      parameter(vacwidth=.15d0)
      real*8 massit(1:200)

c-----|----------------------------------------------------------------|X
      massit(101)=0.14d0     
      massit(104)=0.77d0

      mrho=massit(104)
      mpi=massit(101)

      vacselfrho=0.d0
      if(M.le.2.d0*mpi) return

      pi=4.d0*datan(1.d0)

      ome0=mrho/2.d0
      p0=sqrt(ome0**2-mpi**2)
      
      grho2=vacwidth/mrho/(p0/ome0)**3*48.d0*pi

      reself=grho2*M**2/48d0/pi**2*(
     & (1d0-4d0*mpi**2/M**2)**1.5*
     &     dlog(abs((1d0+sqrt(1d0-4.d0*mpi**2/M**2))/
     &     (1d0-sqrt(1-4.d0*mpi**2/M**2))))+
     &     8.d0*mpi**2*(1d0/M**2-1d0/mrho**2)-2d0*(p0/ome0)**3*
     &     dlog((ome0+p0)/mpi))
      
      imself=-grho2*M**2/48d0/pi*(1.d0-4d0*mpi**2/M**2)**1.5

      vacselfrho=dcmplx(reself,imself)
      
      return
      end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE SETARRAYZERO
c..
c*****|****************************************************************|X

      SUBROUTINE setarrayzero
      implicit none

      integer separ
      parameter (separ=999)

c..   write event separator
c        write(44,'(5i12)') separ, separ, separ, separ, separ
        write(71,'(6i12)') separ, separ, separ, separ, separ, separ

      end

c-----|----------------------------------------------------------------|X

c//////////////////////////////////////////////////////////////////////|X
c.. ROUTINES FOR POLYNOMIAL INTERPOLATION                               X
c//////////////////////////////////////////////////////////////////////|X

c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE POLINT
c..
c*****|****************************************************************|X      

      SUBROUTINE polint(xa,ya,n,x,y,dy)

      INTEGER n,NMAX
      REAL*8 dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      REAL*8 den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)

      ns=1
      dif=abs(x-xa(1))
c      write(0,*)'CALL POLINT: x=',x,dif
      do 11 i=1,n
c      write(0,*)'in polint: i, xa(i), ya(i)',i,xa(i),ya(i)
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
c          write(0,*)'in polint: m, i',m,i
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)stop 'failure in polint'
          den=w/den
c          write(0,*)'in polint: ho, hp, w, den',ho,hp,w,den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
c        write(0,*)'in polint: y =',y
c        write(0,*)'in polint: dy =',dy
        y=y+dy

13    continue
      return
      END

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE POLIN2
c..
c*****|****************************************************************|X  

      SUBROUTINE polin2(x1a,x2a,ya,m,n,x1,x2,y,dy)
      INTEGER m,n,NMAX,MMAX
      REAL*8 dy,x1,x2,y,x1a(m),x2a(n),ya(m,n)
      PARAMETER (NMAX=20,MMAX=20)
CU    USES polint
      INTEGER j,k
      REAL*8 ymtmp(MMAX),yntmp(NMAX)
      do 12 j=1,m
        do 11 k=1,n
          yntmp(k)=ya(j,k)
c           write(0,*) 'x2a,yntmp',x2a(k),yntmp(k) ! Debug only
11      continue
c          write(0,*) 'n,x2,x1a',n,x2,x1a(j) ! Debug only
        call polint(x2a,yntmp,n,x2,ymtmp(j),dy)
c        call ratint(x2a,yntmp,n,x2,ymtmp(j),dy)
c         write(0,*) 'j,ymtmp(j)',j,ymtmp(j) ! Debug only
12    continue
c        write(0,*) 'm,x1',m,x1 ! Debug only
      call polint(x1a,ymtmp,m,x1,y,dy)
c      call ratint(x1a,ymtmp,m,x1,y,dy)
      return
      END

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE QSIMP_HAD
c..
CU. ***USES trapzd***
c..
c.. Returns as 's' the integral of the function 'func' from 'a' to 'b'. 
c.. The parameters EPS can be setto the desired fractional accuracy 
c.. and JMAX so that 2 to the power JMAX-1 is the maximum allowed 
c.. number of steps. Integration is performed by Simpson's rule.
c..
c*****|****************************************************************|X 
     
      SUBROUTINE qsimp_had(func,a,b,s)
      implicit none
      
      INTEGER JMAX
      REAL*8 a,b,func,s,EPS
      EXTERNAL func
c      PARAMETER (EPS=1.d-6, JMAX=40)
c      PARAMETER (EPS=1.d-5, JMAX=35)
      PARAMETER (EPS=1.d-5, JMAX=35)
      INTEGER j
      REAL*8 os,ost,st
c      write(0,*)'***SUBROUTINE qsimp_had***'
c      write(0,*)'a,b,s',a,b,s ! Debug only
      ost=-1.d30
      os= -1.d30
      do 11 j=1,JMAX
c      write(0,*)'***CALL trapzd***'
        call trapzd(func,a,b,st,j)
c      write(0,*)'***END trapzd***'
        s=(4.d0*st-ost)/3.d0
c        write(0,*)'s',s ! Debug only
c        write(0,*) 'j,sum',j,st ! Debug only
        if (abs(s-os).lt.EPS*abs(os)) return
        if (s.eq.0.d0.and.os.eq.0.d0.and.j.gt.6) return
        os=s
        ost=st
11    continue
      write(0,*)'too many steps in qsimp_had'
      s=0.0d0
c      stop 'too many steps in qsimp_had'
c      write(0,*)'***END qsimp_had***'
      return
      END

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE TRAPZD
c..
c.. This routine computes the nth stage of refinement of an extended 
c.. trapezoidal rule. 'func' is input as the name of the function to be 
c.. integrated between limits 'a' and 'b', also input. When called with n=1, 
c.. the routine returns as s the crudest estimate of integral of f(x) 
c.. from a to b. Subsequent calls with n=2,3,... (in that sequential order) 
c.. will improve the accuracy of s by adding 2n-2 additional interior 
c.. points. 's' should not be modified between sequential calls!
c..
c*****|****************************************************************|X 
     
      SUBROUTINE trapzd(func,a,b,s,n)

      INTEGER n
      REAL*8 a,b,s,func
      EXTERNAL func
      INTEGER it,j
      REAL*8 del,sum,tnm,x

c      write(0,*)'trapzd:a,b,s,n',a,b,s,n ! Debug only
      if (n.eq.1) then
        s=0.5d0*(b-a)*(func(a)+func(b))
c        write(0,*)'func(a),func(b),s',func(a),func(b),s ! Debug only
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5d0*del
        sum=0.d0
c        write(0,*) 'it',it ! Debug only
        do 11 j=1,it
          sum=sum+func(x)
c          write(0,*)j,x,func(x),sum ! Debug only
          x=x+del
11      continue
        s=0.5d0*(s+(b-a)*sum/tnm)
      endif
      return
      END

c-----|----------------------------------------------------------------|X

      SUBROUTINE RATINT(XA,YA,N,X,Y,DY)
      INTEGER NMAX,N,NS,I,M
      REAL*8 TINY
      PARAMETER (NMAX=10,TINY=1.D-25)
      REAL*8 XA(N),YA(N),C(NMAX),D(NMAX),X,Y,H,HH,DY,T,DD
c      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      HH=ABS(X-XA(1))
      DO 11 I=1,N
        H=ABS(X-XA(I))
        IF (H.EQ.0.)THEN
          Y=YA(I)
          DY=0.0
          RETURN
        ELSE IF (H.LT.HH) THEN
          NS=I
          HH=H
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)+TINY
 11   CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          W=C(I+1)-D(I)
          H=XA(I+M)-X
          T=(XA(I)-X)*D(I)/H
          DD=T-C(I+1)
          IF(DD.EQ.0.) STOP 'ERROR IN RATINT'
          DD=W/DD
          D(I)=C(I+1)*DD
          C(I)=T*DD
 12     CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
 13   CONTINUE
      RETURN
      END

c-----|----------------------------------------------------------------|X

      SUBROUTINE qtrap(func,a,b,s)  
      INTEGER JMAX  
      REAL*8 a,b,func,s,EPS  
      EXTERNAL func  
c      PARAMETER (EPS=1.e-6, JMAX=20)  
      PARAMETER (EPS=1.d-5, JMAX=35)  
CU    USES trapzd  
      INTEGER j  
      REAL*8 olds  
      olds=-1.d30  
      do 11 j=1,JMAX  
        call trapzd(func,a,b,s,j)  
        if (abs(s-olds).lt.EPS*abs(olds)) return  
        if (s.eq.0.0d0.AND.olds.eq.0.0d0.AND.j.gt.6) return  
        olds=s  
 11   continue  
      write(0,*)'too many steps in qtrap'  
      s=0.0d0
      return
      END  

c-----|----------------------------------------------------------------|X

