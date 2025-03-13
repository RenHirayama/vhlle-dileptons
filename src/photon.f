c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE RAPP_PHOTON
c..
c..   This program returns in-medium resp. vacuum photon rates dR/dM 
c..   according to the Rapp Spectral Functions. We use a parametrization
c..   which is given in Heffernan et al., PRC 91, 027902 (2015), Eqs.
c..   (2) - (7) 
c..
c..  where mu_B: baryon chemical potential [in units GeV]
c..           T: temperature [GeV]
c..       mu_pi: pion chemical potenital [GeV]
c..        mu_k: kaon chemical potential [GeV]
c..           q: invariant mass [GeV]
c..       dR/dq: in-medium photon rate [1/(GeV**3)]  
c..    
c..  Note that usually the photon rates are plotted as q0*dR/dq3 with
c..  [dR/dq]=(1/q0)*pi*(q**2)*[q0*dR/dq3]
c..
c*****|****************************************************************|X

      SUBROUTINE rapp_photon(temp,mub,pipot,gce,vxce,vyce,
     &                     vzce,vol4,multi,beta_lab,dt,timestep,lambda) 
      implicit none

      include 'defs.f'

c-----|----------------------------------------------------------------|X

      real*8 lambda,temp,mub,pipot,kpot,gce,vxce,vyce,vzce
      real*8 mass,p0l,pxl,pyl,pzl
      integer it_step
      real*8 contr,vol4
      real*8 factor,rate,result
      real*8 ammax,ammin  

      integer type(1:10)

      integer multi,loops,l
      integer flagrho,vac

      integer noe,ityp
      real*8 bev,acce,accp

      real*8 p0_po_lab,px_po_lab,py_po_lab,pz_po_lab
      real*8 p0_el_lab,px_el_lab,py_el_lab,pz_el_lab
      real*8 beta_lab

      real*8 dt,time
      integer timestep
      
      integer pc
      real*8 qup,qdown
      real*8 minq

      real*8 fugacity
      
c.. functions
      real*8 rate_rho
      external rate_rho

c. Common Blocks
      real*8 T,mu
      common /th_dyn/T,mu
      logical vacuum
      common /vflag/vacuum
      integer meson
      common /typedef/meson

c.. Steps to increase statistics
      logical zwstep
      common /steps/zwstep
      real*8 zwmin(1:18),zwmax(1:18)
      integer j,jj

c-----|----------------------------------------------------------------|

      T=0.0d0
      mu=0.0d0

c.. Set common variables
      T=temp
c.. NOTE: Parametrization of spectral function only works up to a mu_B
c.. of approx. 400 MeV. If higher values are obtained, set mu_B manually
c.. to 450 MeV (realistic results down to 0.2 GeV).
c.. This only for the rho contribution!
       if(mub.le.0.45d0) then
         mu=mub
       elseif(mub.gt.0.45d0) then
         mu=0.45d0
       endif

c.. Parametrization at low momenta becomes inaccurate for high temperatures
       if(t.le.0.175d0) then
         minq=1.0d-1
       elseif(t.gt.0.175d0) then
         minq=4.0d-1
       endif

c.. Restrict calculation to temperatures below 200 MeV
       if(t.gt.0.200d0) t=0.200d0

c      write(0,*)'vol4',vol4 ! Debug only
   
      vacuum=.FALSE.
      vac=0

      meson=0
     

c.. FOR PHOTONS: MASS IS ZERO!

      mass=0.0d0

c.. Parameters
c... Time
      time=dt*dble(timestep)

c-----|----------------------------------------------------------------|X
c. **** Loop n-times (to increase statistics) ****
      
c*********************!
      loops=1         !
c*********************!                             

      do 678 l=1,loops

       ityp=104 ! rho meson SF

c       write(0,*)'ityp=',ityp ! Debug only
c       write(0,*)'T,mub,pipot,kpt',T,mub,pipot,kpt ! Debug only
       flagrho=1

c    calculate mass and momentum integrated rate

c. **** Loop over p-slices and perform a SIMPSON INTEGRATION for each slice****
      rate=0d0
      result=0d0

      do 707 pc=1,15
       qdown=minq+(pc-1)*0.35d0
       qup=qdown+0.35d0
      

c      write(0,*)'Simpson Integration' ! Debug only 
c      write(0,*)'ammin,ammax,result',ammin,ammax,result ! Debug only
      CALL qsimp_had(rate_rho,qdown,qup,result)
c      write(0,*)'ammin,ammax,result',ammin,ammax,result ! Debug only
c      write(0,*)'Integration finished' ! Debug only 

c     constant factor
         factor=1.0d0
c         factor=alpha_em**2*4.d0/pi**2*massit(104)**4/grho**2
         rate=result*factor
         if(rate.eq.0.0d0) cycle
c     conversion to fm^{-4}
c      rate=rate/(hqc**4)
 
c      write(0,*) 'rate',rate ! Debug only

c. **** Generate 3-momenta ****   
c      write(0,*)'vac,mass,gce,vxce,vyce,vzce',vac,mass,gce,vxce,vyce,
c     &                                        vzce ! Debug only 
      CALL rho_ph_momdist(qup,qdown,gce,vxce,vyce,vzce,p0l,pxl,pyl,pzl)
c      write(0,*)'p0l,pxl,pyl,pzl',p0l,pxl,pyl,pzl ! Debug only


c.. In case of too many interations
c.. p0l is set to zero in r_momdist
 
c      write(0,*)'p0l,pxl,pyl,pzl',p0l,pxl,pyl,pzl ! Debug only 
      if(p0l.eq.0d0) then
       write(0,*)'p0l = 0, return'
       write(0,*)'meson,p0,temp,rate',meson,p0l,temp,rate ! Debug only
       cycle
      endif            

c.. **** Additional FUGACITY factor ****
c... The parametrization pertains to hadronic matter in
c... chemical equilibrium (CE), i.e., for mu_B=mu_\bar{B}
c... without any meson chemical potentials (and therefore
c... mu_B_i=mu_N (for each B_i=N,Delta,N*,...). An extension 
c... of the rate parametrizations to fully incorporate the
c... mu_i dependencies is not practical. However, their 
c... leading effect can be rather accurately captured by 
c... fugacity factors. As mu_rho=2*mu_pi, an extra over-
c... all factor of (z_pi)^2 has to be added.
  
      fugacity=dexp(2.0d0*pipot/T)

c... In Heffernan et al., PRC 91, 027902 (2015) it is further 
c... argued, that enhanced resonance abundances in URHICs, can be
c... accounted for by including an extra factor of z_pi, together
c... with together with using mu_B=mu_N (which is the case for the
c... EoS used for the present study).

      if(beta_lab.gt.0.95d0) fugacity=dexp(3.0d0*pipot/T)

c.. Determine cell contribution
       contr=rate*vol4*(1.0d0-lambda)*dble(multi)*fugacity/dble(loops)

c       contr=rate*(1.0d0-lambda)*dble(multi)/dble(loops)
c     &       /(4.0d0/3.0d0*3.14d0*p0l)

       if(contr.lt.0.0d0.OR.contr.gt.(5.0d0*vol4)) then
        write(0,*)'Suspicious weight!'
        write(0,*)'Weight,vol4,T,mub,pipot,fug',contr,vol4,T,mub,pipot,fugacity
        contr=0.0d0
        cycle
       endif

c. **** Write into output file f71 ****
      noe=1
      bev=0.0d0

c      acce=sqrt(vxce**2+vyce**2)
c      accp=vzce

      acce=pipot
      accp=0.0d0
c........................................................................
c. NOTE: The extended output format writes out ***lab-system momenta*** !
c.       for e+ and e-                                                  !
c.......................................................................!
      if(ext_out) then !extended output format
       write(71,556)ityp,contr,p0l,pxl,pyl,pzl,dt,time,vol4,
     & flagrho,mub,temp,noe,lambda,acce,accp
      else !standard output format
       write(71,555)contr,p0l,pxl,pyl,pzl,flagrho,mub,temp
      endif

 555  format(e14.7,4f12.7,i9,2f12.7)
 556  format(I5,1X,8(E16.8E3,1X),I4,1X,2(E16.8E3,1X),I9,1X,
     &F12.8,1X,2(E16.8,1x))

 707  continue

 678  continue


c-----|----------------------------------------------------------------|X

      return
      end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE RHO_PH_MOMDIST
c..
c..  This subroutine generates the dilepton momenta for rho-meson
c..  emission
c*****|****************************************************************|X

      SUBROUTINE rho_ph_momdist(qup,qdown,g,vx,vy,vz,p0,px,py,pz) 
      implicit none

c-----|----------------------------------------------------------------|X

      real*8 temp,mub,pichem,kchem

      real*8 m,g,vx,vy,vz,p,px,py,pz
      integer i,vac

      real*8 ranff
      external ranff

      real*8 rate_rho
      external rate_rho

      real*8 qup,qdown

      real*8 ftrial,massit
      real*8 r_drdmd3qmax_hr
      real*8 q0,p0,qx,qy,qz,q,phi,cost,sint,v2
      real*8 xmin,xmax,x,x2,c,r_Fint_hr,fb,compf,truef
      real*8 r_imDrho_hr

      real*8 aux,h,hmax,htest,qtest,j
      real*8 test1, test2
      real*8 factor

      real*8 pi
      parameter (pi=3.1415926535) ! useful constant
      real*8 alpha_em
      parameter (alpha_em=1.0d0/137.0d0)
      real*8 hqc
      parameter (hqc=0.197327d0)

c     common variables
      real*8 T,mu
      common /th_dyn/T,mu
      integer seedinit
      common /randomf/ seedinit

c-----|----------------------------------------------------------------|X

      vac=0

c.. Rename variables
      temp=T
      mub=mu

c..generate random phi
      phi=ranff(seedinit)*2.d0*pi
c..generate random cos($\theta$)
      cost=-1d0+ranff(seedinit)*(1d0+1.d0)
     

c.. maximum of the distribution......................
      qtest=qdown-1.0d-4
      if(qdown.le.0.00d0) qtest=1.0d-5
      hmax=0.0d0

362   continue
        htest=rate_rho(qtest)
        if(htest.gt.hmax) hmax=htest
        if(qtest.lt.(qup).AND.qtest.lt.5.5d0) then
         qtest=qtest+0.00125d0 
         goto 362
      endif

c....................................................

      hmax=1.1d0*hmax

      j=0
 2    continue
      j=j+1

      q=qdown+(ranff(seedinit))*(qup-qdown)
      h=ranff(seedinit)*hmax
      
      if(rate_rho(q).gt.hmax) then
       aux=rate_rho(q)
       write(0,*)'hmax too small'
       write(0,*)'hmax,q',hmax,q
       write(0,*)'rate_rho(q)',aux
       write(0,*)'maximum violated by', (aux-hmax)/hmax
       if((aux-hmax)/hmax.gt.0.5d0) stop
c      stop
      endif
      
c     mass could not be generated
      if(j.gt.10000) then
         q=0.d0
         write(0,*)'too many trials to generate momentum'
         write(0,*) 'number of trials',j
         return
      endif
      
      if(h.gt.rate_rho(q)) goto 2


c.. Find cartesian coordinates
       sint=sqrt(1.d0-cost**2)
       q0=q
       qx=q*cos(phi)*sint
       qy=q*sin(phi)*sint
       qz=q*cost
c.. Boost to cf
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


       return 
       end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. FUNCTION RATE_RHO
c..
c..   Imaginary part of the in-medium meson propagator
c*****|****************************************************************|X

      real*8 FUNCTION rate_rho(q)
      implicit none

      include 'defs.f'

c-----|----------------------------------------------------------------|X

      real*8 q,temp,mub

      real*8 rate_nomub
      external rate_nomub
  
      real*8 func_rho
      external func_rho

      real*8 T,mu
      common /th_dyn/T,mu
c-----|----------------------------------------------------------------|X
c.... Rename variables
      temp=T
      mub=mu

      rate_rho=rate_nomub(q,temp)*func_rho(q,mub,temp)*4.0d0*pi*q

      return
      end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. FUNCTION RATE_NOMUB
c..
c..  Interpolates In-medium vector meson propagator from table
c*****|****************************************************************|X

      real*8 FUNCTION rate_nomub(q,temp)
      implicit none

      include 'defs.f'

c-----|----------------------------------------------------------------|X

      real*8 q,temp

      real*8 a,b,c

c-----|----------------------------------------------------------------|X

      a=-31.21d0+353.61d0*temp-1739.4d0*temp**2+3105.0d0*temp**3
      b=-5.513d0-42.2d0*temp+333.0d0*temp**2-570.0d0*temp**3
      c=-6.153d0+57.0d0*temp-134.61d0*temp**2+8.31d0*temp**3

      rate_nomub=dexp(a*q+b+(c/(q+0.2d0)))

      return
      end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. FUNCTION FUNC_RHO
c..
c*****|****************************************************************|X

      real*8 FUNCTION func_rho(q,mub,temp)
      implicit none

      include 'defs.f'

c-----|----------------------------------------------------------------|X

      real*8 q,mub,temp

      real*8 n,p,r,s,v,w,alpha,beta,eta,d,k,m
 
      real*8 rho_dRd3q

c-----|----------------------------------------------------------------|X
      
      n=-0.04d0+2.3d0*temp-12.8d0*temp**2
      p=23.66d0-354.0d0*temp+1175.0d0*temp**2
      r=-54.3d0+742.6d0*temp-2350.0d0*temp**2

      s=-22.11d0+808.7d0*temp-11604.4d0*temp**2+81700.0d0*temp**3
     &  -282480.0d0*temp**4+384116.0d0*temp**5
      v=-1.6d0-121.7d0*temp+1775.0d0*temp**2-5516.0d0*temp**3
      w=-9.874d0+469.0d0*temp-4371.5d0*temp**2+11000.0d0*temp**3

      alpha=84.784d0-3028.6d0*temp+42434.0d0*temp**2-291390.0d0*temp**3
     &      +981000.0d0*temp**4-1295400.0d0*temp**5
      beta=59.64d0-726.46d0*temp+1093.4*temp**2+4256.0d0*temp**3
      eta=-73.9d0+458.66d0*temp+2450.0d0*temp**2-12348.0d0*temp**3
     
      
      d=n*mub+p*mub**2+r*mub**3
      k=s*mub+v*mub**2+w*mub**3
      m=alpha*mub+beta*mub**2+eta*mub**3

      rho_dRd3q=dexp(d-(k/(q**2))-(m/(q)))
 
      func_rho=rho_dRd3q
 
      return
      end

c-----|----------------------------------------------------------------|X
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccc BREMSSTRAHLUNG cccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE BREMS_PHOTON
c..
c..   This program returns the bremsstrahlung contributions for photon
c..   spectra.
c..
c*****|****************************************************************|X
c*****|****************************************************************|X

      SUBROUTINE brems_photon(temp,mub,pipot,kpot,gce,vxce,vyce,
     &                     vzce,vol4,multi,beta_lab,dt,timestep,lambda) 
      implicit none

      include 'defs.f'

c-----|----------------------------------------------------------------|X

      real*8 lambda,temp,mub,pipot,kpot,gce,vxce,vyce,vzce
      real*8 mass,p0l,pxl,pyl,pzl
      integer it_step
      real*8 contr,vol4
      real*8 factor,rate,result
      real*8 ammax,ammin  

      integer type(1:10)

      integer multi,loops,l
      integer flag,vac

      integer noe,ityp
      real*8 bev,acce,accp

      real*8 p0_po_lab,px_po_lab,py_po_lab,pz_po_lab
      real*8 p0_el_lab,px_el_lab,py_el_lab,pz_el_lab
      real*8 beta_lab

      real*8 dt,time
      integer timestep
      
      integer pc
      real*8 qup,qdown

      real*8 fugacity

      real*8 minq
      
c.. functions
      real*8 rate_brems
      external rate_brems

c. Common Blocks
      real*8 T,mu
      common /th_dyn/T,mu
      logical vacuum
      common /vflag/vacuum
      integer meson
      common /typedef/meson

c.. Steps to increase statistics
      logical zwstep
      common /steps/zwstep
      real*8 zwmin(1:18),zwmax(1:18)
      integer j,jj

c.. Bremsstrahlung contribution
      integer brcont
      common /bremsstr/brcont

c-----|----------------------------------------------------------------|

      T=0.0d0
      mu=0.0d0

c.. Set common variables
      T=temp
      mu=mub

c.. Parametrization at low momenta becomes inaccurate for high temperatures
      if(t.le.0.175d0) then
        minq=1.0d-1
      elseif(t.gt.0.175d0) then
        minq=4.0d-1
      endif

c.. Restrict calculation to temperatures below 200 MeV
      if(t.gt.0.200d0) t=0.200d0

c      write(0,*)'vol4',vol4 ! Debug only
   
      vacuum=.FALSE.
      vac=0

      meson=0

c.. FOR PHOTONS: MASS IS ZERO!

      mass=0.0d0

c.. Parameters
c... Time
      time=dt*dble(timestep)

c-----|----------------------------------------------------------------|X
c. **** Loop n-times for different bremsstrahlung contributions ****
      
      loops=6  !DO NOT CHANGE THIS!!! 
      brcont=0                            

      do 678 l=1,loops

       brcont=brcont+1
       ityp=330+l
       flag=1+l

c       write(0,*)'ityp=',ityp ! Debug only
c       write(0,*)'T,mub,pipot,kpt',T,mub,pipot,kpt ! Debug only

c    calculate mass and momentum integrated rate

c. **** Loop over p-slices and perform a SIMPSON INTEGRATION for each slice****
      rate=0d0
      result=0d0

      do 707 pc=1,15
       qdown=minq+(pc-1)*0.35d0
       qup=qdown+0.35d0
      

c      write(0,*)'Simpson Integration' ! Debug only 
c      write(0,*)'ammin,ammax,result',ammin,ammax,result ! Debug only
      CALL qsimp_had(rate_brems,qdown,qup,result)
c      write(0,*)'ammin,ammax,result',ammin,ammax,result ! Debug only
c      write(0,*)'Integration finished' ! Debug only 

c     constant factor
         factor=1.0d0
c         factor=alpha_em**2*4.d0/pi**2*massit(104)**4/grho**2
         rate=result*factor
         if(rate.eq.0.0d0) cycle
c     conversion to fm^{-4}
c      rate=rate/(hqc**4)
 
c      write(0,*) 'rate',rate ! Debug only

c. **** Generate 3-momenta ****   
c      write(0,*)'vac,mass,gce,vxce,vyce,vzce',vac,mass,gce,vxce,vyce,
c     &                                        vzce ! Debug only 
      CALL br_ph_momdist(qup,qdown,gce,vxce,vyce,vzce,p0l,pxl,pyl,pzl)
c      write(0,*)'p0l,pxl,pyl,pzl',p0l,pxl,pyl,pzl ! Debug only


c.. In case of too many interations
c.. p0l is set to zero in r_momdist
 
c      write(0,*)'p0l,pxl,pyl,pzl',p0l,pxl,pyl,pzl ! Debug only 
      if(p0l.eq.0d0) then
       write(0,*)'p0l = 0, return'
       write(0,*)'meson,p0,temp,rate',meson,p0l,temp,rate ! Debug only
       cycle
      endif            
 
c.. Additional FUGACITY factor

      if(brcont.eq.1) then  
       fugacity=dexp(2.0d0*pipot/T)+(0.2d0*dexp(pipot/T)*dexp(kpot/T))
      elseif(brcont.eq.2) then    
       fugacity=dexp(3.0d0*pipot/T)
      elseif(brcont.eq.3) then
       fugacity=dexp(2.0d0*pipot/T)*dexp(kpot/T)
      elseif(brcont.eq.4) then
       fugacity=dexp(pipot/T)*dexp(kpot/T)
      elseif(brcont.eq.5) then
       fugacity=dexp(pipot/T)
      elseif(brcont.eq.6) then
       fugacity=dexp(kpot/T)
      endif

c.. Determine cell contribution
       contr=rate*vol4*(1.0d0-lambda)*dble(multi)*fugacity

c       contr=rate*(1.0d0-lambda)*dble(multi)
c     &       /(4.0d0/3.0d0*3.14d0*p0l)

       if(contr.lt.0.0d0.OR.contr.gt.(5.0d0*vol4)) then
        write(0,*)'Suspicious weight!'
        write(0,*)'Weight,vol4,T,mub,pipot,fug',contr,vol4,T,mub,pipot,fugacity
        contr=0.0d0
        cycle
       endif

c. **** Write into output file f71 ****
      noe=1
      bev=0.0d0

c      acce=sqrt(vxce**2+vyce**2)
c      accp=vzce

       acce=pipot
       accp=kpot
c........................................................................
c. NOTE: The extended output format writes out ***lab-system momenta*** !
c.       for e+ and e-                                                  !
c.......................................................................!
      if(ext_out) then !extended output format
       write(71,556)ityp,contr,p0l,pxl,pyl,pzl,dt,time,vol4,
     & flag,mub,temp,noe,lambda,acce,accp
      else !standard output format
       write(71,555)contr,p0l,pxl,pyl,pzl,flag,mub,temp
      endif

 555  format(e14.7,4f12.7,i9,2f12.7)
 556  format(I5,1X,8(E16.8E3,1X),I4,1X,2(E16.8E3,1X),I9,1X,
     &F12.8,1X,2(E16.8,1x))

 707  continue

 678  continue


c-----|----------------------------------------------------------------|X

      return
      end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE BR_PH_MOMDIST
c..
c..  This subroutine generates the dilepton momenta for rho-meson
c..  emission
c*****|****************************************************************|X

      SUBROUTINE br_ph_momdist(qup,qdown,g,vx,vy,vz,p0,px,py,pz) 
      implicit none

c-----|----------------------------------------------------------------|X

      real*8 temp,mub,pichem,kchem

      real*8 m,g,vx,vy,vz,p,px,py,pz
      integer i,vac

      real*8 ranff
      external ranff

      real*8 rate_brems
      external rate_brems

      real*8 qup,qdown

      real*8 ftrial,massit
      real*8 r_drdmd3qmax_hr
      real*8 q0,p0,qx,qy,qz,q,phi,cost,sint,v2
      real*8 xmin,xmax,x,x2,c,r_Fint_hr,fb,compf,truef

      real*8 aux,h,hmax,htest,qtest,j
      real*8 test1, test2
      real*8 factor

      real*8 pi
      parameter (pi=3.1415926535) ! useful constant
      real*8 alpha_em
      parameter (alpha_em=1.0d0/137.0d0)
      real*8 hqc
      parameter (hqc=0.197327d0)

c     common variables
      real*8 T,mu
      common /th_dyn/T,mu
      integer seedinit
      common /randomf/ seedinit

c-----|----------------------------------------------------------------|X

      vac=0

c.. Rename variables
      temp=T
      mub=mu

c..generate random phi
      phi=ranff(seedinit)*2.d0*pi
c..generate random cos($\theta$)
      cost=-1d0+ranff(seedinit)*(1d0+1.d0)
     

c.. maximum of the distribution......................
      qtest=qdown-1.0d-4
      if(qtest.le.0.0d0) qtest=1.0d-5
      hmax=0.0d0

362   continue
        htest=rate_brems(qtest)
        if(htest.gt.hmax) hmax=htest
        if(qtest.lt.(qup).AND.qtest.lt.5.5d0) then
         qtest=qtest+0.00125d0 
         goto 362
      endif

c....................................................

      hmax=1.1d0*hmax

      j=0
 2    continue
      j=j+1

      q=qdown+(ranff(seedinit))*(qup-qdown)
      h=ranff(seedinit)*hmax
      
      if(rate_brems(q).gt.hmax) then
       aux=rate_brems(q)
       write(0,*)'hmax too small'
       write(0,*)'hmax,q',hmax,q
       write(0,*)'rate_brems(q)',aux
       write(0,*)'maximum violated by', (aux-hmax)/hmax
       if((aux-hmax)/hmax.gt.0.5d0) stop
c      stop
      endif
      
c     mass could not be generated
      if(j.gt.10000) then
         q=0.d0
         write(0,*)'too many trials to generate momentum'
         write(0,*) 'number of trials',j
         return
      endif
      
      if(h.gt.rate_brems(q)) goto 2


c.. Find cartesian coordinates
       sint=sqrt(1.d0-cost**2)
       q0=q
       qx=q*cos(phi)*sint
       qy=q*sin(phi)*sint
       qz=q*cost
c.. Boost to cf
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


       return 
       end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. FUNCTION RATE_BREMS
c..
c..   Photon emission rate dR/d3q from pi-pi bremsstrahlung according
c..   to Eqs. (8) and (9) in Heffernan et al., PRC 91, 027902 (2015)
c*****|****************************************************************|X

      real*8 FUNCTION rate_brems(q)
      implicit none

c-----|----------------------------------------------------------------|X
c.. Useful constants
      real*8 pi
      parameter (pi=3.1415926535) 
      real*8 alpha_em
      parameter (alpha_em=1/137.0d0)
      real*8 hqc
      parameter (hqc=0.197327d0)
      real*8 sq2
      parameter (sq2=dsqrt(2.0d0))
      real*8 alpiq
      parameter (alpiq=2.0d0)

c.. Parametrizations
      real*8 alpha_B,beta_B,gamma_B,delta_B
      real*8 a0,a1,a2,a3
      real*8 b0,b1,b2,b3
      real*8 c0,c1,c2,c3
      real*8 d0,d1,d2,d3

c.. Form factor
      real*8 ff
c.. Ext. function
      real*8 ff_brems
      external ff_brems

c.. Variables
      real*8 q,temp

c.. Common variables
      real*8 T,mu
      common /th_dyn/T,mu
c.. Bremsstrahlung contribution identifier
      integer brcont
      common /bremsstr/brcont

c-----|----------------------------------------------------------------|X
      temp=T  
   
      rate_brems=0.0d0 
       
      ff=ff_brems(q)

      if(brcont.eq.1) then !pi-pi + pi-K bremsstrahlung
        alpha_B=-16.28d0+62.45d0*temp-93.4d0*temp**2-7.5*temp**3
        beta_B=-35.54d0+414.8d0*temp-2054d0*temp**2+3718.8d0*temp**3
        gamma_B=0.7364d0-10.7d0*temp+56.32*temp**2-103.5*temp**3
        delta_B=-2.51d0+58.152d0*temp-318.24d0*temp**2+610.7d0*temp**3

        rate_brems=dexp(alpha_B+beta_B*q+gamma_B*q**2+delta_B/(q+0.2d0))
     &             *(4.0d0*pi*q**2)/q

        return
        
      elseif(brcont.eq.2) then ! pi+rho->pi+gamma (omega t-channel exchg.)
c.....OLD PARAMETRIZATION:
c        rate_brems=ff*temp**(2.8d0)*dexp(-(1.461d0*temp**(2.3094d0)
c     &             +0.727d0)/((2.0d0*temp*q)**(0.86d0))
c     &             +(0.566d0*temp**(1.4094d0)-0.9957d0)*q/temp)
c     &             *(4.0d0*pi*q**2)/q

c.....NEW PARAMETRIZATION (H. v. Hees, priv. comm., 03.08.2015):       
        a0=1.0d-31*(-9.589d-09*temp**(-8.072d0)-0.0784d0*temp)
        a1=0.0438d0*temp**(0.493d0)+0.918d0
        a2=0.00506d0*temp**(-1.338d0)-0.00296d0
        a3=-1.493d0*temp**(-1.169d0)-0.667d0

        rate_brems=ff*a0*(q**a1-1.12d0*q)*exp(77.0d0*q**a2+a3*q**0.75d0)
     &             *(4.0d0*pi*q**2)/q

        return
          
      elseif(brcont.eq.3) then ! pi+K*->pi+gamma
c.....OLD PARAMETRIZATION:
c        rate_brems=ff*temp**(3.75d0)*dexp(-(0.35d0)/((2.0d0*temp*q)**(1.05d0))
c     &             +(2.3894d0*temp**(0.03435d0)-3.222d0)*q/temp)
c     &             *(4.0d0*pi*q**2)/q

c.....NEW PARAMETRIZATION (H. v. Hees, priv. comm., 03.08.2015):
        a0=2.2030d0*temp**(4.968d0)-1.026d-05
        a1=-0.1039d0*temp**(-1.374d0)+0.2253d0
        a2=-0.1382d0*temp**(-1.124d0)-0.1382d0
        a3=-0.9063d0*temp**(-0.952d0)+0.349d0

        rate_brems=ff*a0/sq2*q**a1*dexp(a2*q**(-0.97d0)+a3*q**1.08d0)
     &             *(4.0d0*pi*q**2)/q

        return

      elseif(brcont.eq.4) then ! pi+K->K*+gamma
c.....OLD PARAMETRIZATION:
c        rate_brems=ff*(1.0d0/temp**3)*dexp(-(5.4018d0
c     &             *temp**(-0.6864d0)-1.51d0)*(2.0d0*temp*q)**(0.07d0)
c     &             -0.91d0*q/temp)*(4.0d0*pi*q**2)/q

c.....NEW PARAMETRIZATION (H. v. Hees, priv. comm., 03.08.2015):
        b0=0.00946d0*temp**5.72d0-1.017d-08
        b1=-temp**0.439d0+1.87d0
        b2=-1.225d0*temp**(-0.949d0)+0.901d0

        rate_brems=ff*b0/sq2*q**b1*dexp(b2*q**(0.99d0)+5.0d0*q**(-0.29d0))
     &             *(4.0d0*pi*q**2)/q

        return

      elseif(brcont.eq.5) then ! rho+K->K+gamma
c.....OLD PARAMETRIZATION:
c        rate_brems=ff*temp**(3.5d0)*dexp(-(0.9386d0
c     &             *temp**(1.551d0)+0.634d0)/((2.0d0*temp*q)**(1.01d0))
c     &             +(0.568d0*temp**(0.5397d0)-1.164d0)*q/temp)
c     &             *(4.0d0*pi*q**2)/q

c.....NEW PARAMETRIZATION (H. v. Hees, priv. comm., 03.08.2015):
        c0=0.588d0*temp**(3.03d0)+9.899d-05
        c1=-0.442d0*temp**(-1.138d0)+0.879d0
        c2=-0.6953d0*temp**(-0.992d0)+0.0389d0
        c3=-0.6348d0*temp**(-1.084d0)-0.157d0

        rate_brems=ff*c0/sq2*q**c1*dexp(c2*q**(1.12d0)+c3*q**(-0.78d0))
     &             *(4.0d0*pi*q**2)/q
        return

      elseif(brcont.eq.6) then ! K*+K->pi+gamma
c.....OLD PARAMETRIZATION:
c        rate_brems=ff*temp**(3.7d0)*dexp(-(6.096d0*temp**(1.889d0)
c     &   +1.0299d0)/((2.0d0*temp*q)**(-1.613d0*temp**(2.162d0)+0.975d0))
c     &   -0.96d0*q/temp)*(4.0d0*pi*q**2)/q

c.....NEW PARAMETRIZATION (H. v. Hees, priv. comm., 03.08.2015):
        d0=1.396d0*temp**(4.74d0)+7.62d-04
        d1=-4.567d0*temp**(-0.433d0)+6.978d0
        d2=-0.432d0*temp**(-1.174d0)-0.388d0
        d3=-1.245d0*temp**(-0.887d0)+1.503d0

        rate_brems=ff*d0/sq2*q**d1*dexp(d2*q**(1.1d0)+d3*q**(-0.8d0))
     &             *(4.0d0*pi*q**2)/q
        return
      end if

      end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. FUNCTION FF_BREMS
c..
c..   Dipole form factors for meson-meson bremsstrahlung according
c..   to Eqs. (10) and (11) in Turbide et al., PRC 69, 014903 (2004)
c*****|****************************************************************|X

      real*8 FUNCTION ff_brems(q)
      implicit none

c-----|----------------------------------------------------------------|X

      real*8 q,tav

      real*8 alpiq2
      parameter (alpiq2=2.0d0)

c.. Bremsstrahlung contribution identifier
      integer brcont
      common /bremsstr/brcont

c-----|----------------------------------------------------------------|X
c     From Ralf Rapp via Hendrik van Hees (private comm.)

      if(brcont.ge.2.AND.brcont.le.4.AND.q.gt.0.2d0) then

      tav=34.5096d0*q**0.737d0-67.557d0*q**(0.7584d0)+
     $           32.858d0*q**0.7806d0
      ff_brems=(alpiq2/(alpiq2-tav))**8

      elseif(brcont.ge.5.AND.brcont.le.6.AND.q.gt.0.2d0) then

      tav=-76.424d0*q**0.6236d0+36.944d0*q**0.6604d0+39.0448d0*
     $           q**0.5873d0
      ff_brems=(alpiq2/(alpiq2-tav))**8
        
      else

      ff_brems=1.0d0

      endif

      return
      end

c-----|----------------------------------------------------------------|X
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccc QUARK-GLUON-PLASMA cccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE QGP_PHOTON
c..
c..   This program returns in-medium resp. vacuum photon rates dR/dM 
c..   according to the Rapp Spectral Functions. We use a parametrization
c..   which is given in Heffernan et al., PRC 91, 027902 (2015), Eqs.
c..   (2) - (7) 
c..
c..  where mu_B: baryon chemical potential [in units GeV]
c..           T: temperature [GeV]
c..       mu_pi: pion chemical potenital [GeV]
c..        mu_k: kaon chemical potential [GeV]
c..           q: invariant mass [GeV]
c..       dR/dq: in-medium photon rate [1/GeV]  
c..    
c..  Note that usually the photon rates are plotted as q0*dR/dq3 with
c..  [dR/dq]=(1/q0)*4*pi*(q**2)*[q0*dR/dq3]
c..
c*****|****************************************************************|X

      SUBROUTINE qgp_photon(temp,mub,pipot,gce,vxce,vyce,
     &                     vzce,vol4,multi,beta_lab,dt,timestep,lambda) 
      implicit none

      include 'defs.f'

c-----|----------------------------------------------------------------|X

      real*8 lambda,temp,mub,pipot,kpot,gce,vxce,vyce,vzce
      real*8 mass,p0l,pxl,pyl,pzl
      integer it_step
      real*8 contr,vol4
      real*8 factor,rate,result
      real*8 ammax,ammin  

      integer type(1:10)

      integer multi,loops,l
      integer flag,vac

      integer noe,ityp
      real*8 bev,acce,accp

      real*8 p0_po_lab,px_po_lab,py_po_lab,pz_po_lab
      real*8 p0_el_lab,px_el_lab,py_el_lab,pz_el_lab
      real*8 beta_lab

      real*8 dt,time
      integer timestep
      
      integer pc
      real*8 qup,qdown

c.. functions
      real*8 rate_qgp
      external rate_qgp

c. Common Blocks
      real*8 T,mu
      common /th_dyn/T,mu
      logical vacuum
      common /vflag/vacuum
      integer meson
      common /typedef/meson

c-----|----------------------------------------------------------------|

      T=0.0d0
      mu=0.0d0

c.. Set common variables
      T=temp
      mu=mub

c      write(0,*)'vol4',vol4 ! Debug only
   
      vacuum=.FALSE.
      vac=0

      meson=0


c.. FOR PHOTONS: MASS IS ZERO!

      mass=0.0d0

c.. Parameters
c... Time
      time=dt*dble(timestep)

c-----|----------------------------------------------------------------|X
c. **** Loop n-times (to increase statistics) ****
      
c*********************!
      loops=1         !
c*********************!                             

      do 678 l=1,loops

       ityp=222 ! QGP

c       write(0,*)'ityp=',ityp ! Debug only
c       write(0,*)'T,mub',T,mub ! Debug only
       flag=8

c    calculate mass and momentum integrated rate

c. **** Loop over p-slices and perform a SIMPSON INTEGRATION for each slice****
      rate=0d0
      result=0d0

      do 707 pc=1,15
       qdown=1.0d-1+(pc-1)*0.35d0
       qup=qdown+0.35d0
      

c      write(0,*)'Simpson Integration' ! Debug only 
      CALL qsimp_had(rate_qgp,qdown,qup,result)
c      write(0,*)'qdown,qup,result',qdown,qup,result ! Debug only
c      write(0,*)'Integration finished' ! Debug only 

c     constant factor
         factor=1.0d0
c         factor=alpha_em**2*4.d0/pi**2*massit(104)**4/grho**2
         rate=result*factor
         if(rate.eq.0.0d0) cycle
c     conversion to fm^{-4}
c      rate=rate/(hqc**4)
 
c      write(0,*) 'rate',rate ! Debug only

c. **** Generate 3-momenta ****   
c      write(0,*)'vac,mass,gce,vxce,vyce,vzce',vac,mass,gce,vxce,vyce,
c     &                                        vzce ! Debug only 
      CALL qgp_ph_momdist(qup,qdown,gce,vxce,vyce,vzce,p0l,pxl,pyl,pzl)
c      write(0,*)'p0l,pxl,pyl,pzl',p0l,pxl,pyl,pzl ! Debug only


c.. In case of too many interations
c.. p0l is set to zero in r_momdist
 
c      write(0,*)'p0l,pxl,pyl,pzl',p0l,pxl,pyl,pzl ! Debug only 
      if(p0l.eq.0d0) then
       write(0,*)'p0l = 0, return'
       write(0,*)'meson,p0,temp,rate',meson,p0l,temp,rate ! Debug only
       cycle
      endif            
      
c.. Determine cell contribution
       contr=rate*vol4*lambda*dble(multi)/dble(loops)

c       contr=rate*(1.0d0-lambda)*dble(multi)/dble(loops)
c     &       /(4.0d0/3.0d0*3.14d0*p0l)

       if(contr.lt.0.0d0.OR.contr.gt.(5.0d0*vol4)) then
        write(0,*)'Suspicious weight!'
        write(0,*)'Weight,vol4,T,mub,pipot',contr,vol4,T,mub,pipot
        contr=0.0d0
        cycle
       endif

c. **** Write into output file f71 ****
      noe=1
      bev=0.0d0
c      acce=sqrt(vxce**2+vyce**2)
c      accp=vzce

       acce=0.0d0
       accp=0.0d0
c........................................................................
c. NOTE: The extended output format writes out ***lab-system momenta*** !
c.       for e+ and e-                                                  !
c.......................................................................!
      if(ext_out) then !extended output format
       write(71,556)ityp,contr,p0l,pxl,pyl,pzl,dt,time,vol4,
     & flag,mub,temp,noe,lambda,acce,accp
      else !standard output format
       write(71,555)contr,p0l,pxl,pyl,pzl,flag,mub,temp
      endif

 555  format(e14.7,4f12.7,i9,2f12.7)
 556  format(I5,1X,8(E16.8E3,1X),I4,1X,2(E16.8E3,1X),I9,1X,
     &F12.8,1X,2(E16.8,1x))

 707  continue

 678  continue


c-----|----------------------------------------------------------------|X

      return
      end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE QGP_PH_MOMDIST
c..
c..  This subroutine generates the dilepton momenta for rho-meson
c..  emission
c*****|****************************************************************|X

      SUBROUTINE qgp_ph_momdist(qup,qdown,g,vx,vy,vz,p0,px,py,pz) 
      implicit none

c-----|----------------------------------------------------------------|X

      real*8 temp,mub,pichem,kchem

      real*8 m,g,vx,vy,vz,p,px,py,pz
      integer i,vac

      real*8 ranff
      external ranff

      real*8 rate_qgp
      external rate_qgp

      real*8 qup,qdown

      real*8 ftrial,massit
      real*8 r_drdmd3qmax_hr
      real*8 q0,p0,qx,qy,qz,q,phi,cost,sint,v2
      real*8 xmin,xmax,x,x2,c,r_Fint_hr,fb,compf,truef
      real*8 r_imDrho_hr

      real*8 aux,h,hmax,htest,qtest,j
      real*8 test1, test2
      real*8 factor

      real*8 pi
      parameter (pi=3.1415926535d0) ! useful constant
      real*8 alpha_em
      parameter (alpha_em=1.0d0/137.0d0)
      real*8 hqc
      parameter (hqc=0.197327d0)

c     common variables
      real*8 T,mu
      common /th_dyn/T,mu
      integer seedinit
      common /randomf/ seedinit

c-----|----------------------------------------------------------------|X

      vac=0

c.. Rename variables
      temp=T
      mub=mu

c..generate random phi
      phi=ranff(seedinit)*2.d0*pi
c..generate random cos($\theta$)
      cost=-1d0+ranff(seedinit)*(1d0+1.d0)
     

c.. maximum of the distribution......................
      qtest=qdown-1.0d-4
      if(qtest.le.0.0d0) qtest=0.1d-5
      hmax=0.0d0

362   continue
        htest=rate_qgp(qtest)
        if(htest.gt.hmax) hmax=htest
        if(qtest.lt.(qup).AND.qtest.lt.5.5d0) then
         qtest=qtest+0.00125d0 
         goto 362
      endif

c....................................................

      hmax=1.1d0*hmax

      j=0
 2    continue
      j=j+1

      q=qdown+(ranff(seedinit))*(qup-qdown)
      h=ranff(seedinit)*hmax
      
      if(rate_qgp(q).gt.hmax) then
       aux=rate_qgp(q)
       write(0,*)'hmax too small'
       write(0,*)'hmax,q',hmax,q
       write(0,*)'rate_qgp(q)',aux
       write(0,*)'maximum violated by', (aux-hmax)/hmax
       if((aux-hmax)/hmax.gt.0.5d0) stop
c      stop
      endif
      
c     mass could not be generated
      if(j.gt.10000) then
         q=0.d0
         write(0,*)'too many trials to generate momentum'
         write(0,*) 'number of trials',j
         return
      endif
      
      if(h.gt.rate_qgp(q)) goto 2


c.. Find cartesian coordinates
       sint=sqrt(1.d0-cost**2)
       q0=q
       qx=q*cos(phi)*sint
       qy=q*sin(phi)*sint
       qz=q*cost
c.. Boost to cf
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
  
c      write(0,*)'q0,qx,qy,qz,p0,px,py,pz',q0,qx,qy,qz,p0,px,py,pz !Debug only

       return 
       end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. FUNCTION RATE_QGP
c..
c..   Returns the QGP photon emission rate dR/d3q
c*****|****************************************************************|X

      real*8 FUNCTION rate_qgp(q)
      implicit none

c-----|----------------------------------------------------------------|X

      real*8 als,prfph,gsT,xRR,C22,Cab,Ctot,aux,dRd3q

      real*8 q,temp

      real*8 pi
      parameter (pi=3.141592653d0) ! useful constant
      real*8 al
      parameter (al=1.0d0/137.0d0)
      real*8 hqc
      parameter (hqc=0.197327d0)

      real*8 T,mu
      common /th_dyn/T,mu

c.. QGP type
      integer phqgp,hgpar
      common /inputph/ hgpar,phqgp

c-----|----------------------------------------------------------------|X
      temp=T

c.. From Ralf Rapp via Hendrik van Hees (private comm.)
     
      als=6.0d0*pi/(27.0d0*dlog(temp/0.022d0))
      prfph=als*temp**2*(5.d0/9.d0)
      gsT=sqrt(als*4.0d0*pi)
      xRR=q/temp
      C22=0.041d0/xRR-0.3615d0+1.01d0*dexp(-1.35d0*xRR)
      Cab=sqrt(1.5d0)*(0.548d0/xRR**(1.5d0)*dlog(12.28d0+1.0d0/xRR)
     $     +0.133d0*xRR/dsqrt(1.d0+xRR/16.27d0))
      Ctot=0.5d0*log(2.0d0*xRR)+C22+Cab

      aux=dexp(-xRR)

      if(phqgp.eq.1) then
c.. FULL PARAMETRIZATION

      dRd3q=al/pi**2*prfph*aux/(1.0d0+aux)*(dlog(sqrt(3.d0)/gsT)+Ctot)
 
      rate_qgp=dRd3q*(4.0d0*pi*q**2)/q/(hqc**4)

      elseif(phqgp.eq.0) then
c.. pQCD RATES

      dRd3q=al/(2.0d0*pi**2)*prfph*aux*(dlog(1.0d0+(2.912d0/gsT**2)*xRR))

      rate_qgp=dRd3q*(4.0d0*pi*q**2)/q/(hqc**4)

      endif

      return
      end

c-----|----------------------------------------------------------------|X
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccc pi-rho-omega EMISS. ccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE PROM_PHOTON
c..
c..   This program returns the bremsstrahlung contributions for photon
c..   spectra.
c..
c*****|****************************************************************|X
c*****|****************************************************************|X

      SUBROUTINE prom_photon(temp,mub,pipot,gce,vxce,vyce,
     &                     vzce,vol4,multi,beta_lab,dt,timestep,lambda) 
      implicit none

      include 'defs.f'

c-----|----------------------------------------------------------------|X

      real*8 lambda,temp,mub,pipot,kpot,gce,vxce,vyce,vzce
      real*8 mass,p0l,pxl,pyl,pzl
      integer it_step
      real*8 contr,vol4
      real*8 factor,rate,result
      real*8 ammax,ammin  

      integer type(1:10)

      integer multi,loops,l
      integer flag,vac

      integer noe,ityp
      real*8 bev,acce,accp

      real*8 p0_po_lab,px_po_lab,py_po_lab,pz_po_lab
      real*8 p0_el_lab,px_el_lab,py_el_lab,pz_el_lab
      real*8 beta_lab

      real*8 dt,time
      integer timestep
      
      integer pc
      real*8 qup,qdown

      real*8 fugacity

      real*8 minq
      
c.. functions
      real*8 rate_prom
      external rate_prom

c. Common Blocks
      real*8 T,mu
      common /th_dyn/T,mu
      logical vacuum
      common /vflag/vacuum
      integer meson
      common /typedef/meson

c.. Steps to increase statistics
      logical zwstep
      common /steps/zwstep
      real*8 zwmin(1:18),zwmax(1:18)
      integer j,jj

c.. pi-rho-omega contribution
      integer promcont
      common /prom/promcont

c-----|----------------------------------------------------------------|

      T=0.0d0
      mu=0.0d0

c.. Set common variables
      T=temp
      mu=mub

c.. Parametrization at low momenta becomes inaccurate for 
c.. high temperatures
      if(t.le.0.175d0) then
        minq=1.0d-1
      elseif(t.gt.0.175d0) then
        minq=4.0d-1
      endif

c.. Restrict calculation to temperatures below 200 MeV
      if(t.gt.0.200d0) t=0.200d0

c      write(0,*)'vol4',vol4 ! Debug only
   
      vacuum=.FALSE.
      vac=0

      meson=0

c.. FOR PHOTONS: MASS IS ZERO!

      mass=0.0d0

c.. Parameters
c... Time
      time=dt*dble(timestep)

c-----|----------------------------------------------------------------|X
c. **** Loop n-times for different pi-rho-omega contributions ****
      
      loops=3  !DO NOT CHANGE THIS!!! 
      promcont=0                            

      do 678 l=1,loops

       promcont=promcont+1
       ityp=550+l
       flag=10+l

c       write(0,*)'ityp=',ityp ! Debug only
c       write(0,*)'T,mub,pipot,kpt',T,mub,pipot,kpt ! Debug only

c    calculate mass and momentum integrated rate

c. **** Loop over p-slices and perform a SIMPSON INTEGRATION for each slice****
      rate=0d0
      result=0d0

      do 707 pc=1,15
       qdown=minq+(pc-1)*0.35d0
       qup=qdown+0.35d0
      

c      write(0,*)'Simpson Integration' ! Debug only 
c      write(0,*)'ammin,ammax,result',ammin,ammax,result ! Debug only
      CALL qsimp_had(rate_prom,qdown,qup,result)
c      write(0,*)'ammin,ammax,result',ammin,ammax,result ! Debug only
c      write(0,*)'Integration finished' ! Debug only 

c     constant factor
         factor=1.0d0
c         factor=alpha_em**2*4.d0/pi**2*massit(104)**4/grho**2
         rate=result*factor
         if(rate.eq.0.0d0) cycle
c     conversion to fm^{-4}
c      rate=rate/(hqc**4)
 
c      write(0,*) 'rate',rate ! Debug only

c. **** Generate 3-momenta ****   
c      write(0,*)'vac,mass,gce,vxce,vyce,vzce',vac,mass,gce,vxce,vyce,
c     &                                        vzce ! Debug only 
      CALL prom_ph_momdist(qup,qdown,gce,vxce,vyce,vzce,p0l,pxl,pyl,pzl)
c      write(0,*)'p0l,pxl,pyl,pzl',p0l,pxl,pyl,pzl ! Debug only


c.. In case of too many interations
c.. p0l is set to zero in r_momdist
 
c      write(0,*)'p0l,pxl,pyl,pzl',p0l,pxl,pyl,pzl ! Debug only 
      if(p0l.eq.0d0) then
       write(0,*)'p0l = 0, return'
       write(0,*)'meson,p0,temp,rate',meson,p0l,temp,rate ! Debug only
       cycle
      endif            
 
c.. Additional FUGACITY factor
      if(promcont.eq.1) then
       fugacity=1.0d0
      elseif(promcont.eq.2) then
       fugacity=dexp(4.0d0*pipot/T)
      elseif(promcont.eq.3) then
       fugacity=dexp(5.0d0*pipot/T) 
      endif

c.. Determine cell contribution
       contr=rate*vol4*(1.0d0-lambda)*dble(multi)*fugacity

c       contr=rate*(1.0d0-lambda)*dble(multi)
c     &       /(4.0d0/3.0d0*3.14d0*p0l)

       if(contr.lt.0.0d0.OR.contr.gt.(5.0d0*vol4)) then
        write(0,*)'Suspicious weight!'
        write(0,*)'Weight,vol4,T,mub,pipot,fug',contr,vol4,T,mub,pipot,fugacity
        contr=0.0d0
        cycle
       endif

c. **** Write into output file f71 ****
      noe=1
      bev=0.0d0

c      acce=sqrt(vxce**2+vyce**2)
c      accp=vzce

       acce=pipot
       accp=0.0d0
c........................................................................
c. NOTE: The extended output format writes out ***lab-system momenta*** !
c.       for e+ and e-                                                  !
c.......................................................................!
      if(ext_out) then !extended output format
       write(71,556)ityp,contr,p0l,pxl,pyl,pzl,dt,time,vol4,
     & flag,mub,temp,noe,lambda,acce,accp
      else !standard output format
       write(71,555)contr,p0l,pxl,pyl,pzl,flag,mub,temp
      endif

 555  format(e14.7,4f12.7,i9,2f12.7)
 556  format(I5,1X,8(E16.8E3,1X),I4,1X,2(E16.8E3,1X),I9,1X,
     &F12.8,1X,2(E16.8,1x))

 707  continue

 678  continue


c-----|----------------------------------------------------------------|X

      return
      end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE PROM_PH_MOMDIST
c..
c..  This subroutine generates the dilepton momenta for rho-meson
c..  emission
c*****|****************************************************************|X

      SUBROUTINE prom_ph_momdist(qup,qdown,g,vx,vy,vz,p0,px,py,pz) 
      implicit none

c-----|----------------------------------------------------------------|X

      real*8 temp,mub,pichem,kchem

      real*8 m,g,vx,vy,vz,p,px,py,pz
      integer i,vac

      real*8 ranff
      external ranff

      real*8 rate_prom
      external rate_prom

      real*8 qup,qdown

      real*8 ftrial,massit
      real*8 q0,p0,qx,qy,qz,q,phi,cost,sint,v2
      real*8 xmin,xmax,x,x2,c,r_Fint_hr,fb,compf,truef

      real*8 aux,h,hmax,htest,qtest,j
      real*8 test1, test2
      real*8 factor

      real*8 pi
      parameter (pi=3.1415926535) ! useful constant
      real*8 alpha_em
      parameter (alpha_em=1.0d0/137.0d0)
      real*8 hqc
      parameter (hqc=0.197327d0)

c     common variables
      real*8 T,mu
      common /th_dyn/T,mu
      integer seedinit
      common /randomf/ seedinit

c-----|----------------------------------------------------------------|X

      vac=0

c.. Rename variables
      temp=T
      mub=mu

c..generate random phi
      phi=ranff(seedinit)*2.d0*pi
c..generate random cos($\theta$)
      cost=-1d0+ranff(seedinit)*(1d0+1.d0)
     

c.. maximum of the distribution......................
      qtest=qdown-1.0d-4
      if(qtest.le.0.0d0) qtest=1.0d-5
      hmax=0.0d0

362   continue
        htest=rate_prom(qtest)
        if(htest.gt.hmax) hmax=htest
        if(qtest.lt.(qup).AND.qtest.lt.5.5d0) then
         qtest=qtest+0.00125d0 
         goto 362
      endif

c....................................................

      hmax=1.1d0*hmax

      j=0
 2    continue
      j=j+1

      q=qdown+(ranff(seedinit))*(qup-qdown)
      h=ranff(seedinit)*hmax
      
      if(rate_prom(q).gt.hmax) then
       aux=rate_prom(q)
       write(0,*)'hmax too small'
       write(0,*)'hmax,q',hmax,q
       write(0,*)'rate_prom(q)',aux
       write(0,*)'maximum violated by', (aux-hmax)/hmax
       if((aux-hmax)/hmax.gt.0.5d0) stop
c      stop
      endif
      
c     mass could not be generated
      if(j.gt.10000) then
         q=0.d0
         write(0,*)'too many trials to generate momentum'
         write(0,*) 'number of trials',j
         return
      endif
      
      if(h.gt.rate_prom(q)) goto 2


c.. Find cartesian coordinates
       sint=sqrt(1.d0-cost**2)
       q0=q
       qx=q*cos(phi)*sint
       qy=q*sin(phi)*sint
       qz=q*cost
c.. Boost to cf
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


       return 
       end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. FUNCTION RATE_PROM
c..
c..   Photon emission rate dR/d3q from pi-rho-omega complex according
c..   to Eqs. (B1) - (B6) in Holt et al., arXiv:1506.09205 [hep-ph] 
c*****|****************************************************************|X

      real*8 FUNCTION rate_prom(q)
      implicit none

c-----|----------------------------------------------------------------|X

      real*8 pi
      parameter (pi=3.1415926535) ! useful constant
      real*8 alpha_em
      parameter (alpha_em=1/137.0d0)
      real*8 hqc
      parameter (hqc=0.197327d0)
      real*8 a_1,a_2,a_3,a_4,a_5,a_6,a_7

      real*8 q,temp

c.. Common variables
      real*8 T,mu
      common /th_dyn/T,mu
c.. pi-rho-omega contribution identifier
      integer promcont
      common /prom/promcont

c-----|----------------------------------------------------------------|X
      temp=T  
   
      rate_prom=0.0d0 
       
      if(promcont.eq.1) then ! pi+rho->gamma+omega
        a_1=-35.8991d0+460.425d0*temp-2592.04d0*temp**2+5342.32d0*temp**3
        a_2=-41.9725d0+601.952d0*temp-3587.8d0*temp**2+7604.97d0*temp**3
        a_3=0.740436d0-16.7159d0*temp+133.526d0*temp**2-347.589d0*temp**3
        a_4=2.00611d0-3.79343d0*temp+29.3101d0*temp**2-72.8725d0*temp**3
        a_5=-8.33046d0+121.091d0*temp-801.676d0*temp**2+1712.16d0*temp**3
        a_6=17.9029d0-388.5d0*temp+2779.03d0*temp**2-6448.4d0*temp**3
        a_7=-15.622d0+340.651d0*temp-2483.18d0*temp**2+5870.61d0*temp**3

        rate_prom=dexp(a_1*q+a_2+a_3*q**(a_4)+a_5*(q+a_6)**(a_7))
     &            *(4.0d0*pi*q**2)/q 
        return
        
      elseif(promcont.eq.2) then ! pi+omega->gamma+rho
        a_1=-29.4663d0+291.356d0*temp-1301.27d0*temp**2+2102.12d0*temp**3
        a_2=-45.081+688.929d0*temp-4150.15d0*temp**2+8890.76d0*temp**3
        a_3=-0.260076d0+8.92875d0*temp-60.868*temp**2+136.57d0*temp**3
        a_4=2.2663d0-8.30596d0*temp+49.3342d0*temp**2-90.8501d0*temp**3
        a_5=10.2955d0-317.077d0*temp+2412.15d0*temp**2-6020.9d0*temp**3
        a_6=3.12251d0-47.5277d0*temp+222.61d0*temp**2-241.9d0*temp**3
        a_7=-3.39045d0+56.5927d0*temp-336.97d0*temp**2+622.756d0*temp**3

        rate_prom=dexp(a_1*q+a_2+a_3*q**(a_4)+a_5*(q+a_6)**(a_7))
     &            *(4.0d0*pi*q**2)/q   
        return
          
      elseif(promcont.eq.3) then ! rho+omega->gamma+pi
        a_1=-29.6866d0+331.769d0*temp-1618.66d0*temp**2+2918.53d0*temp**3
        a_2=-15.3332d0+90.2225d0*temp-300.185d0*temp**2+428.386d0*temp**3
        a_3=-7.35061d0+109.288d0*temp-630.396d0*temp**2+1227.69d0 *temp**3
        a_4=-10.6044d0+109.1d0*temp-500.718d0*temp**2+872.951d0 *temp**3
        a_5=0.0d0
        a_6=0.0d0
        a_7=0.0d0

        rate_prom=dexp(a_1*q+a_2+(a_3/(q+0.2d0))+(a_4/(q+0.2d0)**2))
     &            *(4.0d0*pi*q**2)/q       

        return
      end if

      end

c-----|----------------------------------------------------------------|X
