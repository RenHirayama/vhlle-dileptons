c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE FOPIEMIT
c..
c..   This program returns 4-pion dilepton rates rates dR/dM 
c..   (up to known constant factors).
c..
c..         mu: chemical potential [GeV]
c..          T: temperature [GeV]
c..          M: invariant mass [GeV]
c..      dR/dM: in-medium dilepton rate [a.u.]
c..   
c..   -----------------------------------------------
c..   |  dR/dM[phys.u.]=dR/dM[a.u.]*factor/(hqc**4) |
c..   -----------------------------------------------
c..
c..   To obtain dR/dM^2 (often plotted in papers) 
c..   one has to further divide for 2M,i.e.     
c..   dR/dM^2[phys.u.]=dR/dM[phys.u.]/2./M
c..
c*****|****************************************************************|X

      SUBROUTINE fopiemit(temp,mub,pipot,gce,vxce,vyce,vzce,vol4,multi,
     &                    beta_lab,dt,timestep,lambda)
      implicit none

      include 'defs.f'

c-----|----------------------------------------------------------------|X

      real*8 lambda,temp,mub,pipot,gce,vxce,vyce,vzce
      real*8 mass,p0l,pxl,pyl,pzl
      integer multi,stat
      real*8 contr,vol4
      real*8 factor,rate,result
      real*8 mmin4pi,mmax4pi
      parameter(mmin4pi=0.7876d0,mmax4pi=2.5624d0)
c      parameter(mmin4pi=0.560d0,mmax4pi=2.75d0)
 
      integer flag4pi
      parameter (flag4pi=4)

      real*8 dt,time
      integer timestep

      integer noe,ityp
      real*8 bev,acce,accp

      real*8 p0_po_lab,px_po_lab,py_po_lab,pz_po_lab
      real*8 p0_el_lab,px_el_lab,py_el_lab,pz_el_lab
      real*8 beta_lab

      real*8 fugacity

c     functions
      real*8 distm4pi
      external distm4pi

c     common variables
      real*8 T,mu
      common /th_dyn/T,mu

      integer i,l

c-----|----------------------------------------------------------------|X

c     rename variables and put in common
      T=temp
      mu=mub

c     minimum mass for integration: mmin4pi

c     maximum mass for integration:  mmax4pi

c     calculate mass and momentum integrated rate
c     "rate" has units [fm^{-4}]

c     integration over mass dependent function--Simpson integration
      rate=0.0d0
      result=0.0d0
      mass=0.0d0
      call qsimp_had(distm4pi,mmin4pi,mmax4pi,result)
      if(result.eq.0.0d0) call qtrap(distm4pi,0.790d0,2.55d0,result)            
      if(result.eq.0.0d0) return
c     constant factor
      factor=temp/4.d0/pi**4
      rate=result*factor
c     conversion to fm^{-4}
      rate=rate/(hqc**4)

c      write(0,*) 'rate',rate !Debug only

c. **** Loop n-times for better mass and momentum statsitics ****

c.. Number of Loops
cccccccccccccccccccccccccccc!
      stat=1                !                            
cccccccccccccccccccccccccccc! 

      do 678 l=1,stat

c. **** Generate mass ****
      CALL fourpi_mdist(mass,mmin4pi,mmax4pi)

c      write(0,*) 'mass',mass !Debug only

c     in case of too many iterations 
c     mass is set to zero in massdist
      if(mass.eq.0) return
      if(mass.gt.mmax4pi) return

c. **** Generate 3-momenta ****
      CALL fourpi_momdist(mass,gce,vxce,vyce,
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

c.. Additional FUGACITY factor
  
      fugacity=dexp(4.0d0*pipot/T)

c.. Determine cell contribution
      contr=rate*(1d0-lambda)*dble(multi)/dble(stat)*vol4*fugacity

c. **** Write into output file f71 ****
      ityp=444
      noe=1
      bev=0.0d0
      acce=1.0d0
      accp=1.0d0
      time=dble(timestep*dt)
c........................................................................
c. NOTE: The extended output format writes out ***lab-system momenta*** !
c.       for e+ and e-                                                  !
c.......................................................................!
      if(vHLLE_out) then
       write(71,557)ityp,contr,mass,p0l,pxl,pyl,pzl,t,mub
      elseif(ext_out) then !extended output format
       write(71,556)ityp,contr,mass,p0_el_lab,px_el_lab,py_el_lab,
     & pz_el_lab,p0_po_lab,px_po_lab,py_po_lab,pz_po_lab,dt,time,vol4,
     & flag4pi,mub,temp,noe,lambda,acce,accp
      else !standard output format
       write(71,555)contr,mass,pxl,pyl,pzl,flag4pi,mub,temp
      endif

 555  format(e14.7,4f12.7,i9)	 
 556  format(I5,1X,13(E16.8E3,1X),I4,1X,2(E16.8E3,1X),I9,1X,
     &F12.8,1X,2(E16.8,1x))
 557  format(I3,e14.7,5f12.7,2f6.3)

 678  continue

      return
      end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE FOURPI_MDIST
c..
c..   This subroutine generates the dilepton mass for 4pi emission.
c..
c..   Input : mmin4pi,mmax4pi (max and min values for the mass)
c..   Output: m (the dilepton mass)
c..
c*****|****************************************************************|X

      SUBROUTINE fourpi_mdist(m,mmin4pi,mmax4pi)
      implicit none

c-----|----------------------------------------------------------------|X

      real*8 mmin4pi,mmax4pi,h,hmax,maxsigma4pi,hm
      real*8 m,mtest
      real*8 ranff
      integer i,j
      real*8 aux
      real*8 zbar
      real*8 hqc
      parameter (hqc  = 0.197327) ! value of $\hbar c$

c     functions
      real*8 distm4pi,bessk1,sigma4pi

c     common variables
      real*8 T,mu
      common /th_dyn/T,mu
      integer seedinit
      common /randomf/ seedinit
      
c-----|----------------------------------------------------------------|X      

c...  root of equation 5*K1(z)-z*K2(z)=0
c      zbar=3.4138d0
c...  majorant of the distribution
c      maxsigma4pi=sigma4pi(1.4785d0)*1.d-7 /(hqc**2)
c      hmax=maxsigma4pi*T**4*zbar**4*bessk1(zbar)
 
      hmax=0.0d0
      mtest=mmin4pi
      do 100 i=1,30
        mtest=mmin4pi+(i*0.025d0)
        hm=distm4pi(mtest)
c        write(0,*)'m,distm',mtest,hm !Debug only
        hmax=max(hmax,hm)
c       write(0,*)'HMAX =',hmax !Debug only
 100  continue

      hmax=1.3d0*hmax


      j=0
 2    continue
      j=j+1
      m=mmin4pi+ranff(seedinit)*(mmax4pi-mmin4pi)
c     write(*,*) 'generated mass:', m
      h=ranff(seedinit)*hmax
      
      if(distm4pi(m).gt.hmax) then
         aux=distm4pi(m)
         write(*,*)'hmax too small'
         write(*,*)'hmax',hmax
         write(*,*)'distm4pi(m)',aux
         write(*,*)'maximum violated by', (aux-hmax)/hmax
         if((aux-hmax)/hmax.gt.1.d-3) stop
c      stop
      endif
      
c     mass could not be generated
      if(j.gt.10000000) then
         m=0.d0
         write(*,*)'from fourpi_mdist:'
         write(*,*)'too many trials to generate mass'
         write(*,*) 'number of trials',j
         return
      endif
      
      if(h.gt.distm4pi(m)) goto 2
      
      return
      end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE FOURPI_MOMDIST
c..
c..   This subroutine generates the dilepton momenta for 4pi emission
c..
c*****|****************************************************************|X

      SUBROUTINE fourpi_momdist(m,g,vx,vy,vz,p0,px,py,pz) 
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
      real*8 T,mub
      common /th_dyn/T,mub
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
         write(*,*) 'numerical rumor in fourpi_momdist'
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
         write(*,*) 'from fourpi_momdist:'
         write(*,*)'comparison function too small'
         write(*,*) 'ERROR'
         stop
      endif
      x2=ranff(seedinit)
      if(x2.gt.truef/compf.and.i.lt.1000000) then
         goto 1
      endif
      if(i.ge.1000000) then
         write(*,*) 'too many iteration in fourpi_momdist' 
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
c.. FUNCTION DISTM4PI
c..
c..   This program returns the dR/dM distribution function 
c..   (without constant factors).
c..
c..   mass : dilepton mass (in GeV)
c..
c*****|****************************************************************|X

      real*8 function distm4pi(mass)
      implicit none

c-----|----------------------------------------------------------------|X

      real*8 mass
      real*8 bessk1,sigma4pi,sigma
      real*8 hqc
      parameter (hqc  = 0.197327) ! value of $\hbar c$

c     common variables
      real*8 T,mu
      common /th_dyn/T,mu

c-----|----------------------------------------------------------------|X

c..   conversion from nb to fm**2
      sigma=sigma4pi(mass)*1.d-7 
c..   conversion from fm**2 to Gev^{-2}
      sigma=sigma/(hqc**2)

      distm4pi=mass**4*bessk1((mass/T))*sigma
      end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. FUNCTION SIGMA4PI
c..
c..   This program returns the cross section for e+e- -> 4pi [nb].
c..
c..   The total cross section is the sum of the
c..   cross section for the reaction  e+ e- --> pi+ pi- pi+ pi-
c..   and the one for the reaction    e+ e- --> pi+ pi- pi0 pi0
c..
c..   The function is obtained as linear interpolation of exp. data
c..   measured by the BaBar collaboration (restricted to ecm<2.6 GeV)
c..   [Phys. Rev. D 71, 052001 (2005) & arXiv:0710.3455 [hep-ex]]
c..
c*****|****************************************************************|X
      
      real*8 FUNCTION sigma4pi(x)
      implicit none

c-----|----------------------------------------------------------------|X

      real*8 x,dx,dy,sbabar,sbabar2
      integer i

c-----|----------------------------------------------------------------|X

c..   ecm bins [GeV]
      real*8 xb(72)
      data xb/0.7875,0.8125,0.8375,0.8625,0.8875,0.9125, 
     &        0.9375,0.9625,0.9875,1.0125,1.0375,1.0625, 
     &        1.0875,1.1125,1.1375,1.1625,1.1875,1.2125,
     &        1.2375,1.2625,1.2875,1.3125,1.3375,1.3625,
     &        1.3875,1.4125,1.4375,1.4625,1.4875,1.5125,
     &        1.5375,1.5625,1.5875,1.6125,1.6375,1.6625, 
     &        1.6875,1.7125,1.7375,1.7625,1.7875,1.8125,
     &        1.8375,1.8625,1.8875,1.9125,1.9375,1.9625,
     &        1.9875,2.0125,2.0375,2.0625,2.0875,2.1125, 
     &        2.1375,2.1625,2.1875,2.2125,2.2375,2.2625,
     &        2.2875,2.3125,2.3375,2.3625,2.3875,2.4125, 
     &        2.4375,2.4625,2.4875,2.5125,2.5375,2.5625 / 
c..   sigma for e+ e- --> pi+ pi- pi+ pi- [nb]
      real*8 yb(72)
      data yb/0.06,0.08,0.11,0.28,0.37,0.44,
     &	      0.36,0.78,0.94,1.14,1.76,2.65,
     &        3.07,3.82,5.02,7.1,7.97,10.5,
     &        12.3,13.4,16.0,18.2,20.2,21.7,
     &        24.9,27.0,28.3,29.3,30.2,29.8,
     &        28.7,26.4,26.0,22.9,22.0,19.8,
     &        17.7,16.2,14.9,12.9,10.7,10.0,
     &        8.29,6.99,6.86,6.23,6.55,6.29,
     &        5.92,5.48,5.72,5.38,5.5,4.6,
     &        4.78,4.73,3.82,3.49,3.55,3.43,
     &        3.11,2.69,3.13,2.51,2.11,2.3,
     &        1.94,2.18,1.76,1.73,1.62,1.69/
c..   sigma for e+ e- --> pi+ pi- pi0 pi0 [nb] (clicked)
      real*8 yb2(72)
      data yb2/0. ,0. ,0. ,0. ,0. ,0.,
     &         1.171,2.957,4.051,5.723,8.49,8.719,
     &         12.01,13.33,14.94,16.73,17.31,19.79, 
     &         20.19,21.51,23.3,25.78,26.3,26.87,
     &         28.37,30.16,31.94,32.06,32.23,31.65,
     &         29.51,29.8,27.89,27.2,25.87,24.83,
     &         23.56,22.58,20.32,19.46,16.86,15.07,
     &         12.47,11.89,10.39,9.235,9.233,8.539,
     &         8.595,9.17,9.11,8.359,7.578, 7.405, 
     &         6.594,6.245,6.201,5.983,5.284,4.891,
     &         4.716,4.41,3.886,3.581,3.537,3.231,
     &         3.188,3.013,2.751,2.795,2.533,2.271/ 
  
c-----|----------------------------------------------------------------|X

      if(x.gt.xb(72)) then
         write(0,*) 'from sbabar: x too high'
         write(0,*) 's=',x      
         sigma4pi=0.0d0
         return
c         stop
      endif
      
      if(x.eq.xb(72)) then
            sbabar=yb(72) 
            sbabar2=yb2(72)
      else
         do i=1,72
            if(x.ge.xb(i).and.x.lt.xb(i+1)) then
               dy=yb(i+1)-yb(i)
               dx=xb(i+1)-xb(i)
               sbabar=yb(i)+dy/dx*(x-xb(i))
               dy=yb2(i+1)-yb2(i)
               sbabar2=yb2(i)+dy/dx*(x-xb(i))
               goto 11
            endif
         enddo
      endif

 11   continue

      sigma4pi=sbabar+sbabar2

      return
      end
      
c-----|----------------------------------------------------------------|X
