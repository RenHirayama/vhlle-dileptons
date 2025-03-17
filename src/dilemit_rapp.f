c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE DILEMIT_RAPP
c..
c..   This program returns in-medium resp. vacuum dilepton rates dR/dM 
c..   according to the Rapp Spectral Functions. It uses tables of the 
c..   kind:
c..
c..               T, rho, mu_pi, mu_k, M, dR/dM
c..
c..  where rho:  effective daryon density [in units of rho_0]
c..          T:  temperature [GeV]
c..      mu_pi:  pion chemical potenital [GeV]
c..       mu_k:  kaon chemical potential [GeV]
c..          M:  invariant mass [GeV]
c..      dR/dM:  in-medium dilepton rate [1/GeV]
c..   
c..   To obtain dR/dM^2 (often plotted in papers) 
c..   one has to further divide for 2M,i.e.     
c..   dR/dM^2[phys.u.]=dR/dM[phys.u.]/2./M
c..
c*****|****************************************************************|X

      SUBROUTINE dilemit_rapp(temp,rhonuc,kpot,pipot,gce,vxce,vyce,
     &                     vzce,vol4,multi,beta_lab,dt,timestep,lambda) 
      implicit none

      include 'defs.f'

c-----|----------------------------------------------------------------|X

      real*8 lambda,temp,rhonuc,pipot,kpot,gce,vxce,vyce,vzce
      real*8 mass,p0l,pxl,pyl,pzl
      integer it_step
      real*8 contr,vol4
      real*8 factor,rate,result
      real*8 ammax,ammin  

      real*8 grho
      parameter (grho=5.03)
      real*8 mminit(1:200)    
      real*8 massit(1:200)

      integer type(1:10)

      integer multi,mesons,l
      integer flagrho,vac

      integer noe,ityp
      real*8 bev,acce,accp

      real*8 p0_po_lab,px_po_lab,py_po_lab,pz_po_lab
      real*8 p0_el_lab,px_el_lab,py_el_lab,pz_el_lab
      real*8 beta_lab

      real*8 dt,time
      integer timestep
      
c.. functions
      real*8 r_distmhad
      external r_distmhad

c. Common Blocks
      real*8 T,rhn,ppt,kpt
      common /th_dyn_r/T,rhn,ppt,kpt
      logical vacuum
      common /vflag/vacuum
      integer meson
      common /typedef/meson

c-----|----------------------------------------------------------------|

c.. Set common variables
      T=temp
      rhn=rhonuc
      ppt=pipot
      kpt=kpot      

c      write(0,*)'vol4',vol4 ! Debug only
   
      vacuum=.FALSE.
      vac=0

      meson=0

c.. Parameters
c      write(0,*)'Dilepton emission' ! Debug only

c.. OMEGA properties
      massit(103)=0.783d0
      mminit(103)=0.277d0
c.. RHO properties
      massit(104)=0.77d0
      mminit(104)=0.277d0
c.. PHI properties
      massit(109)=1.019d0
      mminit(109)=0.277d0

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

c.. MAXIMUM MASS for integration
      ammax=mmax

c.. Time
      time=dt*dble(timestep)

c-----|----------------------------------------------------------------|X
c. **** Loop n-times, once for each vector meson ****

c.. Number of Loops
      if(rates.eq.1) mesons=3 ! RHO, OMEGA AND PHI    
      if(rates.eq.2) mesons=1 ! RHO ONLY   
       
c*********************!
      type(1)=104     !
      type(2)=103     !
      type(3)=109     !
c*********************!                             

      do 678 l=1,mesons
       meson=l
       ityp=type(l)
c       write(0,*)'ityp=',ityp ! Debug only
c       write(0,*)'T,rhn,ppt,kpt',T,rhn,ppt,kpt ! Debug only
       flagrho=l+5

c    calculate mass and momentum integrated rate

c. **** Integration over mass dependent function - SIMPSON INTEGRATION ****
      rate=0d0
      result=0d0

c      write(0,*)'Simpson Integration' ! Debug only 
c      write(0,*)'ammin,ammax,result',ammin,ammax,result ! Debug only
      CALL qsimp_had(r_distmhad,ammin,ammax,result)
c      write(0,*)'ammin,ammax,result',ammin,ammax,result ! Debug only
c      write(0,*)'Integration finished' ! Debug only 

c     constant factor
         factor=1.0d0
c         factor=alpha_em**2*4.d0/pi**2*massit(104)**4/grho**2
         rate=result*factor
c     conversion to fm^{-4}
c      rate=rate/(hqc**4)
 
c      write(0,*) 'rate',rate ! Debug only

c. **** Generate mass ****
      mass=0.0d0
c      write(0,*)'mass',mass ! Debug only
      CALL r_massdist(mass)
c      write(0,*)'mass',mass ! Debug only
   
c.. in case of too many iterations 
c.. mass is set to zero in r_massdist
      if(mass.eq.0) return

c. **** Generate 3-momenta ****   
c      write(0,*)'vac,mass,gce,vxce,vyce,vzce',vac,mass,gce,vxce,vyce,
c     &                                        vzce ! Debug only 
      CALL r_momdist(mass,gce,vxce,vyce,vzce,p0l,pxl,pyl,pzl)
c      write(0,*)'p0l,pxl,pyl,pzl',p0l,pxl,pyl,pzl ! Debug only


c.. In case of too many interations
c.. p0l is set to zero in r_momdist
 
c      write(0,*)'p0l,pxl,pyl,pzl',p0l,pxl,pyl,pzl ! Debug only 
      if(p0l.eq.0d0) then
       write(0,*)'p0l = 0, return'
       write(0,*)'mass,rate',mass,rate ! Debug only
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

c.. Determine cell contribution
       contr=rate*vol4*(1.0d0-lambda)*dble(multi)

c. **** Write into output file f71 ****
      noe=1
      bev=0.0d0
      acce=1.0d0
      accp=1.0d0
c........................................................................
c. NOTE: The extended output format writes out ***lab-system momenta*** !
c.       for e+ and e-                                                  !
c.......................................................................!
      if(ext_out) then !extended output format
       write(71,556)ityp,contr,mass,p0_el_lab,px_el_lab,py_el_lab,
     & pz_el_lab,p0_po_lab,px_po_lab,py_po_lab,pz_po_lab,dt,time,vol4,
     & flagrho,rhn,temp,noe,lambda,acce,accp
      else !standard output format
       write(71,555)contr,mass,pxl,pyl,pzl,flagrho,rhn,temp
      endif

 555  format(e14.7,4f12.7,i9,2f12.7)
 556  format(I5,1X,13(E16.8E3,1X),I4,1X,2(E16.8E3,1X),I9,1X,
     &F12.8,1X,2(E16.8,1x))

 678  continue


c      if(meson.eq.3) stop ! Debug only

      return
      end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE R_MASSDIST
c..
c.. This subroutine generates the dilepton mass for rho-meson emission
c*****|****************************************************************|X

      subroutine r_massdist(m)
      implicit none

      include 'defs.f'

c-----|----------------------------------------------------------------|X

      real*8 temp,rhonuc,pichem,kchem

      real*8 ammax,ammin,h,r_distmhad,hmax,bessk1
      real*8 mminit(1:200)
      real*8 m
      real*8 ranff
      external ranff
      integer j,b
      real*8 r_drdmmax
      real*8 aux
c      real*8 rand

c. Common Blocks
      real*8 T,rhn,ppt,kpt
      common /th_dyn_r/T,rhn,ppt,kpt
  
c-----|----------------------------------------------------------------|X

c.. Rename variables
      temp=T
      rhonuc=rhn
      pichem=ppt
      kchem=kpt
    

      mminit(104)=0.277d0

      ammin=mmin
      ammax=mmax

c      write(0,*)'ammin,ammax',ammin,ammax ! Debug only

c     maximum of the distribution 
      hmax=r_drdmmax(rhonuc,pichem,kchem,temp)

c      write(0,*) 'hmax',hmax ! Debug only

      j=0
 2    continue
      j=j+1

c      write(0,*)'in r_massdist: ammin,ammax',ammin,ammax ! Debug only
c      m=ammin+dble(rand(0))*(ammax-ammin)
      m=ammin+(ranff(seedinit))*(ammax-ammin)
c      write(0,*)'generated mass:', m ! Debug only
c      write(0,*)'r_distmhad(mass):',r_distmhad(m) ! Debug only
c      h=dble(rand())*hmax
      h=ranff(seedinit)*hmax
c      write(0,*)'ranf*hmax=',h ! Debug only
      
      if(r_distmhad(m).gt.hmax) then
         aux=r_distmhad(m)
         write(0,*)'hmax too small'
         write(0,*)'hmax',hmax
         write(0,*)'r_distmhad(m)',aux
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
      
      if(h.gt.r_distmhad(m)) goto 2
      
      end


c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE R_MOMDIST
c..
c..  This subroutine generates the dilepton momenta for rho-meson
c..  emission
c*****|****************************************************************|X

      SUBROUTINE r_momdist(m,g,vx,vy,vz,p0,px,py,pz) 
      implicit none

c-----|----------------------------------------------------------------|X

      real*8 temp,rhonuc,pichem,kchem

      real*8 m,g,vx,vy,vz,p0,px,py,pz
      integer i,vac
c      real*8 ranf
      real*8 ranff
c      external ranff
      real*8 ftrial,massit
      real*8 r_drdmd3qmax
      real*8 q0,qx,qy,qz,q,phi,cost,sint,v2
      real*8 xmin,xmax,x,x2,c,r_Fint,fb,compf,truef
      real*8 r_imDrho
      real*8 pi
      parameter (pi=3.1415926535) ! useful constant
      real*8 s9,s10,eps
      parameter(eps=1.d-6)
c      real*8 rand
      real*8 test1, test2
      real*8 factor

      real*8 grho
      parameter (grho=5.03)
      real*8 alpha_em
      parameter (alpha_em=1/137.0d0)
      real*8 hqc
      parameter (hqc=0.197327d0)

c     common variables
      real*8 qup,qdown,q0up,q0down
      common /qupdow/ qup,qdown
      real*8 T,rhn,ppt,kpt
      common /th_dyn_r/T,rhn,ppt,kpt
      integer seedinit
      common /randomf/ seedinit

c-----|----------------------------------------------------------------|X

      vac=0

c.. Rename variables
      temp=T
      rhonuc=rhn
      pichem=ppt
      kchem=kpt

c.. maximum energy allowed from tabulated momenta
       q0up=sqrt(qup**2+m**2)
c       write(0,*) 'q0up',q0up ! Debug only
c..   minimum energy
       q0down=m

      if(vac.eq.1) then
         ftrial=1.d0
      else
         ftrial=r_drdmd3qmax(rhonuc,pichem,kchem,temp,m)
      endif
c      write(0,*)'ftrial',ftrial ! Debug only
      xmin=0.d0
      c=r_Fint(m,q0down,temp)
      xmax=r_Fint(m,q0up,temp)-c
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
c     1/-1+exp(-q0/temp)
c..the comparison functions is qup*r_imDrhomax
c       x=ranf(0)*xmax
       x=ranff(seedinit)*xmax
c       x=dble(rand())*xmax
       
       if(m/T.lt.10d0) then
             q0=-T*Log(1d0-Exp(x/temp)*(1.d0-Exp(-m/temp)))
c             write(0,*)'temp,x,q0',temp,x,q0 ! Debug only
             if(1d0-Exp(x/temp)*(1.d0-Exp(-m/temp)).le.0d0) then
                write(0,*) 'numerical rumor'
                write(0,*) 'extracted log of non positive number'
                write(0,*) 'logarg',1d0-Exp(x/temp)*(1.d0-Exp(-m/temp))
                write(0,*) 'x,temp,m',x,temp,m
                stop
             endif
       else
c        Boltzmann approximation for m/T>10
          q0=-(temp*Log((temp*Exp(-m/temp)-x)/temp))
c          write(0,*)'temp,x,Boltzm. approx. for q0',temp,x,q0 ! Debug only
          if((temp*Exp(-m/temp)-x)/temp.le.0d0) then
             write(0,*) 'numerical rumor in Bolt. appr.'
             write(0,*) 'extracted log of non positive number'
             write(0,*) 'logarg',(temp*Exp(-m/temp)-x)/temp
             write(0,*) 'x,temp,m',x,temp,m
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
       fb=1.d0/(exp(q0/temp)-1.d0)
c       write(0,*)'boson distr. func.',fb ! Debug only 
       compf=qup*ftrial*fb
c       write(0,*)'compf =',compf ! Debug only 
       if(vac.eq.1) then
           truef=q*fb
c           write(0,*)'in vacuum truef=',truef ! Debug only
       else
          if(q.ge.qdown) then
             truef=q*fb*r_imDrho(q,m)
c             write(0,*)'r_imDrho(q,m)',r_imDrho(q,m) ! Debug only 
c             write(0,*)'in-medium truef=',truef ! Debug only 
          else
c.. for q<qdown r_imDrho is approx with its value at qdown
             truef=q*fb*r_imDrho(qdown,m)
c             write(0,*)'in-medium (q<qdown) truef=',truef ! Debug only 
          endif  
       endif
c.. added for stability
c       write(0,*) 'truef',truef ! Debug only
       if(truef.le.0.d0) then
         write(0,*)'truef < 0!'
         p0=0.d0
       return
       endif
 
       if(compf.lt.truef) then
          write(0,*)'comparison function too small'
          write(0,*)'temp,rhonuc',temp,rhonuc
          write(0,*)'m,q,q0',m,q,q0
          write(0,*)'qup,qdown,q0up,q0down',qup,qdown,q0up,q0down
          write(0,*)'ftrial,fb',ftrial,fb
          write(0,*)'r_imDrho,r_drdmd3qmax',r_imDrho(q,m),
     &                       r_drdmd3qmax(rhonuc,pichem,kchem,temp,m)
          write(0,*)'compf,truef',compf,truef
          write(0,*) 'ERROR'
          stop
       endif
c       x2=ranf(0)
       x2=ranff(seedinit)
c       x2=dble(rand())
c       if(x2.gt.truef/compf) goto 1 
       if(x2.gt.truef/compf.and.i.lt.10000000) then
          goto 1
       endif
       if(i.ge.10000000) then
         write(0,*) 'too many iteration in r_momdist' 
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

c       write(0,*) 'number of trials in omdist:', i ! Debug only

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
c.. FUNCTION R_FINT
c..
c*****|****************************************************************|X

      real*8 FUNCTION r_Fint(m,q0,temp)
      implicit none

c-----|----------------------------------------------------------------|X

      real*8 m,temp,q0

c-----|----------------------------------------------------------------|X
      
      if(m/temp.lt.10d0) then
         r_Fint=temp*log(-1d0+exp(q0/temp))-q0
      else
c.. Boltzmann approximation for m/T>10
	r_Fint=-temp*exp(-q0/temp)
      endif

      return
      end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. FUNCTION R_DRDMD3QMAX
c..
c..  Given values of (mub,temp,mass) defing a point P in the
c..  mub-temp-mass space, this routine returns the
c..  maximum of the function r_imDrho within the 2x2x2 cube 
c..  centered around the point P
c*****|****************************************************************|X

      real*8 FUNCTION r_drdmd3qmax(rhn,ppt,kpt,temp,mass)
      implicit none
  
      include 'defs.f'

c-----|----------------------------------------------------------------|X

      real*8 mass
      real*8 temp,rhn,ppt,kpt
      integer m
      parameter(m=2)
      
      integer i,j,k,o,p,ii,jj,kk,oo,pp

      real*8 dmd3qmx_rh(itmax,jrhomax,kpimax,lkmax,mmmaxrhom)
      real*8 dmd3qmx_om(itmax,jrhomax,kpimax,lkmax,mmmaxrhom)
      real*8 dmd3qmx_ph(itmax,jrhomax,kpimax,lkmax,mmmaxrhom)

      real*8 masstab(mmmaxrhom)
      real*8 masstabmin,masstabmax
      common /vmmass/ masstab,masstabmax,masstabmin
      common /maxdmd3q/dmd3qmx_rh,dmd3qmx_om,dmd3qmx_ph
      
      real*8 mutab(jrhomax),ttab(itmax),pitab(kpimax),ktab(lkmax)
      real*8 mutabmax,mutabmin,ttabmax,ttabmin,pitabmin,pitabmax
      real*8 ktabmin,ktabmax

      common /tbvar/mutab,ttab,pitab,ktab,mutabmax,mutabmin,
     &     ttabmax,ttabmin,pitabmin,pitabmax,ktabmin,ktabmax
 
      integer meson
      common /typedef/meson

c..   auxiliary variables
      real*8 maxaux(jrhomax),dy,res1(lkmax),res2(kpimax),res

c-----|----------------------------------------------------------------|X

c      write(0,*)'temp,rhn,ppt,kpt,meson',temp,rhn,ppt,kpt,meson ! Debug only

      if(rhn.gt.mutabmax) then
         write(0,*) 'chemical potential too high! need bigger table'
         rhn=mutabmax
      endif

c..   determine i (chemical potential index)
      i=int((rhn-mutabmin)/dmur)+1
      ii=min(max(i-(m-1)/2,1),jrhomax+1-m)
      
      if(mass.gt.masstabmax) then
          write(0,*) 'invariant mass too high! need bigger mass table'
         stop
      endif

c..   determine k (mass index)
      k=int((mass-masstabmin)/dmr)+1
      kk=min(max(k-(m-1)/2,1),mmmax+1-m)

c      write(0,*) 'i,ii,k,kk',i,ii,k,kk ! Debug only

      if(temp.gt.ttabmax) then
         write(0,*) 'from r_drdmd3qmax'
         write(0,*) 'temperature too high! need bigger table'
         write(0,*) 'temp,ttabmax',temp,ttabmax
         temp=ttabmax
      endif

c..   determine j (temperature index)
      j=int((temp-ttabmin)/(ttabmax-ttabmin)
     &     *dble(itmax-1))+1
      jj=min(max(j-(m-1)/2,1),itmax+1-m)

c.. Determine o (pion chemical potential index)
      o=int((ppt/0.141d0*3.0d0)+0.5d0)+1
      oo=min(max(o-(m-1)/2,1),kpimax+1-m)

c.. Determine p (pion chemical potential index)
      p=int((kpt/0.450d0*3.0d0)+0.5d0)+1
      pp=min(max(o-(m-1)/2,1),lkmax+1-m)
      

c      write(0,*) 'i,ii',i,ii
c      write(0,*) 'j,jj',j,jj
c      write(0,*) 'o,oo',o,oo
c      write(0,*) 'p,pp',p,pp
c      write(0,*) 'k,kk',k,kk

      res=0.d0
      do o=oo,oo+m-1
       res2(o)=0.0d0
       do p=oo,oo+m-1
         res1(p)=0.0d0
         do i=ii,ii+m-1
          maxaux(i)=0.d0
          do j=1,m
           do k=1,m
                if(meson.eq.1) then
                 maxaux(i)=max(maxaux(i),dmd3qmx_rh(jj+j-1,i,o,p,kk+k-1))
c                 write(0,*)'i,maxaux(i)',i,maxaux(i)
                endif
                if(meson.eq.2) then
                 maxaux(i)=max(maxaux(i),dmd3qmx_om(jj+j-1,i,o,p,kk+k-1))
                endif
                if(meson.eq.3) then
                 maxaux(i)=max(maxaux(i),dmd3qmx_ph(jj+j-1,i,o,p,kk+k-1))
                endif
             enddo
          enddo
c     find maximum on the mu grid
          res1(p)=max(res1(p),maxaux(i))
         enddo
c     find maximum on the pichem grid
         res2(o)=max(res2(o),res1(p))
        enddo
c     find maximum on the kchem grid
        res=max(res,res2(o))

      enddo
      

      r_drdmd3qmax=res      
c     write(0,*)'res',res      

      return
      end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. FUNCTION R_IMDRHO
c..
c..   Imaginary part of the in-medium meson propagator
c*****|****************************************************************|X

      real*8 FUNCTION r_imDrho(q,mass)
      implicit none

      include 'defs.f'

c-----|----------------------------------------------------------------|X

      real*8 rhn,temp,ppt,kpt,q,mass
      real*8 improp
      real*8 dilphsp

c-----|----------------------------------------------------------------|X
 
c.. Determine the mass of the lepton-type 
      if(dimuon) then
         m_lept=mass_muon
      else !dielectron
         m_lept=mass_electron
      endif
   
c      write(0,*)'in r_imDrho: mub,temp,q,mass',mub,temp,q,mass ! Debug only

      r_imDrho=0.d0

c.. Dilepton phase-space
      dilphsp=sqrt(1.d0-4.d0*m_lept**2/mass**2)*
     &   (1.d0+2.d0*m_lept**2/mass**2)

      if(1.d0-4.d0*m_lept**2/mass**2.lt.0.d0) return

      r_imDrho=dabs(improp(q,mass))

      return
      end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. FUNCTION IMPROP
c..
c..  Interpolates In-medium vector meson propagator from table
c*****|****************************************************************|X

      real*8 FUNCTION improp(q,mass)
      implicit none

      include 'defs.f'

c-----|----------------------------------------------------------------|X

      real*8 temp,rhonuc,pichem,kchem,q,mass
      integer i,j,k,n,o,p,ii,jj,kk,nn,oo,pp
      integer m
      parameter (m=2)

      real*8 improp_rho(itmax,jrhomax,kpimax,lkmax,mmmaxrhom,nqmax)
      real*8 improp_ome(itmax,jrhomax,kpimax,lkmax,mmmaxrhom,nqmax)
      real*8 improp_phi(itmax,jrhomax,kpimax,lkmax,mmmaxrhom,nqmax)

      real*8 mutab(jrhomax),ttab(itmax),pitab(kpimax),ktab(lkmax)
      real*8 qtab(nqmax)
      real*8 mutabmax,mutabmin,ttabmax,ttabmin,pitabmin,pitabmax
      real*8 qtabmax,qtabmin,ktabmin,ktabmax
      real*8 masstab(mmmaxrhom)
      real*8 masstabmin,masstabmax

      common /vmmass/ masstab,masstabmax,masstabmin
      common /tbvar/mutab,ttab,pitab,ktab,mutabmax,mutabmin,
     &     ttabmax,ttabmin,pitabmin,pitabmax,ktabmin,ktabmax
      common /tbself/qtab,improp_rho,improp_ome,improp_phi,
     &     qtabmax,qtabmin

      real*8 T,rhn,ppt,kpt
      common /th_dyn_r/T,rhn,ppt,kpt
      integer meson
      common /typedef/meson

c..   auxiliary variables
      real*8 imaux(m,m),dy,ims
      real*8 yim1(jrhomax),yim2(mmmaxrhom),yim3(lkmax),yim4(kpimax)

c-----|----------------------------------------------------------------|X

c.. Rename variables
      temp=T
      rhonuc=rhn
      pichem=ppt
      kchem=kpt

      if(rhonuc.gt.mutabmax) then
         write(0,*) 'chemical potential too high! need bigger table'
         rhonuc=mutabmax
      endif

c.. Determine i (chemical potential index)
      i=int((rhonuc-mutabmin)/dmur)+1
      ii=min(max(i-(m-1)/2,1),jrhomax+1-m)


      if(q.gt.qtabmax) then
         write(0,*) 'momentum too high!need bigger table'
         stop
      endif

c.. Determine k (mass index)
      k=int((mass-masstabmin)/dmr)+1
      kk=min(max(k-(m-1)/2,1),mmmax+1-m)

c.. Determine n (momentum index)
      n=int((q-qtabmin)/dqr)+1
      nn=min(max(n-(m-1)/2,1),nqmax+1-m)

c.. Determine o (pion chemical potential index)
      o=int((pichem/0.141d0*3.0d0)+0.5d0)+1
      oo=min(max(o-(m-1)/2,1),4+1-m)

c.. Determine p (pion chemical potential index)
      p=int((kchem/0.450d0*3.0d0)+0.5d0)+1
      pp=min(max(o-(m-1)/2,1),4+1-m)

      if(temp.gt.ttabmax) then
         write(0,*) 'from self'
         write(0,*) 'temperature too high! need bigger table'
         write(0,*) 'temp,ttabmax',temp,ttabmax
         temp=ttabmax
      endif

c.. Determine j (temperature index)
      j=int((temp-ttabmin)/(ttabmax-ttabmin)
     &     *dble(itmax-1))+1
      jj=min(max(j-(m-1)/2,1),itmax+1-m)
 
c      write(0,*)'indices:temp,rhonuc,pichem,kchem,mass,q',
c     &                j,i,o,p,k,n ! Debug only

      do o=oo,oo+m-1
       do p=pp,pp+m-1
        do k=kk,kk+m-1
         do i=ii,ii+m-1

c-----|----------------------------------------------------------------|X
c.. Map indeces on a (m,m) matrix:                                     |
c.. (1,1)=(jj,kk)  (1,2)=(jj,kk+1) ... (1,m)=(jj,kk+m-1)               |
c.. (2,1)=(jj+1,kk)  ...)                                              |
c..  ...                                                               |
c.. (m,1)=(jj+m-1,kk)   ....            (m,m)=(jj+m-1,kk+m-1)          |
c-----|----------------------------------------------------------------|X

           do j=1,m 
            do n=1,m
             if(meson.eq.1)imaux(j,n)=improp_rho(jj+j-1,i,o,p,k,nn+n-1)
             if(meson.eq.2)imaux(j,n)=improp_ome(jj+j-1,i,o,p,k,nn+n-1)
             if(meson.eq.3)imaux(j,n)=improp_phi(jj+j-1,i,o,p,k,nn+n-1)       
c             if(imaux(j,n).le.0.0d0) write(0,*)'improp < 0:',imaux(j,n)
            enddo
           enddo

c.. Interpolate the (T,q) grid
           call polin2(ttab(jj),qtab(nn),imaux,m,m,
     &          temp,q,yim1(i),dy)
c           write(0,*)'ttab(jj),qtab(nn),yim1(i)',ttab(jj),qtab(nn),yim1(i) !Debug onl
          enddo

c.. Now interpolate results on the mu raw
          CALL polint(mutab(ii),yim1(ii),m,rhonuc,yim2(k),dy)
c          write(0,*)'mutab(ii),yim1(ii),yim2(k),dy',mutab(ii),yim1(ii),yim2(k),dy ! Debug only
         enddo

c.. Now interpolate results on the mass raw
         CALL polint(masstab(kk),yim2(kk),m,mass,yim3(p),dy)
c         write(0,*)'masstab(kk),yim2(kk),yim3(p),dy',masstab(kk),yim2(kk),yim3(p),dy ! Debug only
        enddo

c.. Now interpolate results on the kchem raw
        CALL polint(ktab(pp),yim3(pp),m,kchem,yim4(o),dy)
c         write(0,*)'ktab(pp),yim3(pp),yim4(o),dy',ktab(pp),yim3(pp),yim4(o),dy ! Debug only
       enddo

c.. Now interpolate results on the pichem raw
       CALL polint(pitab(oo),yim4(oo),m,pichem,ims,dy)
c       write(0,*)'pitab(oo),yim4(oo),ims,dy',pitab(oo),yim4(oo),ims,dy ! Debug only
c       write(0,*)'in improp:temp,rhonuc,pichem,kchem,mass,q',
c     &                temp,rhonuc,pichem,kchem,mass,q ! Debug only

      improp=ims
c      if(improp.le.0.0d0) write(0,*)'meson,improp < 0: ',meson,improp ! Debug only
      end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. FUNCTION R_DISTMHAD
c..
c..  dR/dM distribution function (without constant factors)
c*****|****************************************************************|X

      real*8 FUNCTION r_distmhad(mass)
      implicit none

      include 'defs.f'

c-----|----------------------------------------------------------------|X

      real*8 mass
      real*8 temp,rhonuc,pichem,kchem
      integer m
      parameter(m=2)
      real*8 dilphsp
      
      integer i,j,k,o,p,ii,jj,kk,oo,pp

      real*8 masstab(mmmaxrhom)
      real*8 rtab_rho(itmax,jrhomax,kpimax,lkmax,mmmaxrhom)
      real*8 rtab_ome(itmax,jrhomax,kpimax,lkmax,mmmaxrhom)
      real*8 rtab_phi(itmax,jrhomax,kpimax,lkmax,mmmaxrhom)
      real*8 masstabmin,masstabmax
      common /vmmass/ masstab,masstabmax,masstabmin
      common /drdm/rtab_rho,rtab_ome,rtab_phi
      real*8 mutab(jrhomax),ttab(itmax),pitab(kpimax),ktab(lkmax)
      real*8 mutabmax,mutabmin,ttabmax,ttabmin,pitabmax,pitabmin
      real*8 ktabmax,ktabmin

      common /tbvar/mutab,ttab,pitab,ktab,mutabmax,mutabmin,
     &     ttabmax,ttabmin,pitabmin,pitabmax,ktabmin,ktabmax

c.. Auxiliary variables
      real*8 rateaux(m,m),yre(jrhomax),dy,res1(1:4),res2(1:4),res

c. Common Blocks
      real*8 T,rhn,ppt,kpt
      common /th_dyn_r/T,rhn,ppt,kpt
      logical vacuum
      common /vflag/vacuum
      integer meson
      common /typedef/meson
     
c-----|----------------------------------------------------------------|X

c.. Rename variables
      temp=T
      rhonuc=rhn
      pichem=ppt
      kchem=kpt

c.. Determine the mass of the lepton-type 
      if(dimuon) then
         m_lept=mass_muon
      else !dielectron
         m_lept=mass_electron
      endif

      dilphsp=1.0d0

c      write(0,*)'in r_distmhad:temp,rhonuc,pichem,kchem',temp,rhonuc,
c     &                                                 pichem,kchem ! Debug only
      
      if(rhonuc.gt.mutabmax) then
         write(0,*) 'chemical potential too high! need bigger table'
         rhn=mutabmax
      endif

c.. Determine i (chemical potential index)
      i=int((rhonuc-mutabmin)/dmur)+1
      ii=min(max(i-(m-1)/2,1),jrhomax+1-m)
      
      if(mass.gt.masstabmax) then
          write(0,*) 'from r_distmhad'
          write(0,*) 'invariant mass too high! need bigger mass table'
         stop
      endif

c.. Determine m (mass index)
      k=int((mass-masstabmin)/dmr)+1
      kk=min(max(k-(m-1)/2,1),mmmax+1-m)

      if(temp.gt.ttabmax) then
         write(0,*) 'from r_distmhad'
         write(0,*) 'temperature too high! need bigger table'
         write(0,*) 'temp,ttabmax',temp,ttabmax
         temp=ttabmax
      endif

c.. Determine j (temperature index)
      j=int((temp-ttabmin)/(ttabmax-ttabmin)
     &     *dble(itmax-1))+1
      jj=min(max(j-(m-1)/2,1),itmax+1-m)
c      write(0,*) 'ttabmin,ttabmax',ttabmin,ttabmax ! Debug only
c      write(0,*) 'j,jj',j,jj ! Debug only

c.. Determine o (pion chemical potential index)
      o=int((pichem/0.141d0*3.0d0)+0.5d0)+1
      oo=min(max(o-(m-1)/2,1),4+1-m)

c.. Determine p (pion chemical potential index)
      p=int((kchem/0.450d0*3.0d0)+0.5d0)+1
      pp=min(max(o-(m-1)/2,1),4+1-m)
     
      do o=oo,oo+m-1
       do p=pp,pp+m-1 
        do i=ii,ii+m-1
         do j=1,m
          do k=1,m
             if(meson.eq.1) rateaux(j,k)=rtab_rho(jj+j-1,i,o,p,kk+k-1)
             if(meson.eq.2) rateaux(j,k)=rtab_ome(jj+j-1,i,o,p,kk+k-1)
             if(meson.eq.3) rateaux(j,k)=rtab_phi(jj+j-1,i,o,p,kk+k-1) 
          enddo
         enddo

c.. Interpolate the (T,m) grid
c         write(0,*)'CALL polint2' ! Debug only
         CALL polin2(ttab(jj),masstab(kk),rateaux,m,m,
     &        temp,mass,yre(i),dy)
c         write(0,*)'rateaux(m,m),m',rateaux(m,m),m !Debug only
         enddo

c.. Now interpolate results on the mu raw
c         write(0,*)'CALL polint' ! Debug only
c         write(0,*)'mutab(ii),yre(ii),res1(p),dy',mutab(ii),yre(ii),res1(p),dy ! Debug only
         CALL polint(mutab(ii),yre(ii),m,rhonuc,res1(p),dy)
c         write(0,*)'mutab(ii),yre(ii),res1(p),dy',mutab(ii),yre(ii),res1(p),dy ! Debug only
c        write(0,*)'o,p',o,p ! Debug only
        enddo
c.. Now interpolate results on the kchem raw        
        CALL polint(ktab,res1,m,kchem,res2(o),dy)
c        write(0,*)'ktab(pp),kchem,res1(pp),res2(o),dy',ktab(pp),kchem,res1(pp),res2(o),dy ! Debug only
       enddo
c.. Now interpolate results on the kchem raw        
       CALL polint(pitab(oo),res2(oo),m,pichem,res,dy)
c       write(0,*)'pitab(oo),pichem,res2(oo),res,dy',pitab(oo),pichem,res2(oo),res,dy ! Debug only

       if(rates.eq.2) then
        dilphsp=sqrt(1.d0-4.d0*m_lept**2/mass**2)*
     &          (1.d0+2.d0*m_lept**2/mass**2)
        if(mass.le.2*m_lept) dilphsp=0.0d0
       endif

      r_distmhad=0.0d0
      r_distmhad=res*dilphsp      

c      write(0,*)'r_distmhad',r_distmhad ! Debug only
      res=0
 
c      if(.NOT.(pichem.ge.0.0d0)) stop ! Debug only

      return
      end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. FUNCTION R_DRDMMAX
c..
c..  Given values of (mub,temp) defing a point P in the
c..  mub-temp plane, this routine returns the
c..  maximum of the invariant mass distribution dR/dM 
c..  (up to known constant factors)
c*****|****************************************************************|X

      real*8 function r_drdmmax(rhn,ppt,kpt,temp)
      implicit none

      include 'defs.f'
      
c-----|----------------------------------------------------------------|X

      real*8 temp,rhn,ppt,kpt
      integer m
      parameter(m=2)
      integer i,j,o,p,ii,jj,oo,pp

      real*8 rmax_rho(itmax,jrhomax,kpimax,lkmaxrhom)
      real*8 rmax_ome(itmax,jrhomax,kpimax,lkmaxrhom)
      real*8 rmax_phi(itmax,jrhomax,kpimax,lkmaxrhom)
      common /maxdr/rmax_rho,rmax_ome,rmax_phi

      real*8 mutab(jrhomax),ttab(itmax),pitab(kpimax),ktab(lkmax)
      real*8 mutabmax,mutabmin,ttabmax,ttabmin
      real*8 pitabmin,pitabmax,ktabmin,ktabmax

      common /tbvar/mutab,ttab,pitab,ktab,mutabmax,mutabmin,
     &     ttabmax,ttabmin,pitabmin,pitabmax,ktabmin,ktabmax

      integer meson
      common /typedef/meson

c.. auxiliary variables
      real*8 maxaux(jrhomax),dy,res,res1(1:4),res2(1:4)
           
c-----|----------------------------------------------------------------|X

c       write(0,*)'mutabmax,mutabmin,ttabmax,ttabmin,pitabmin,pitabmax,ktabmin,ktabmax'
c       write(0,*)mutabmax,mutabmin,ttabmax,ttabmin,pitabmin,pitabmax,ktabmin,ktabmax !Debug only

       if(rhn.gt.mutabmax) then
         write(0,*) 'chemical potential too high! need bigger table'
         write(0,*) 'mub,mutabmax',rhn,mutabmax
         rhn=mutabmax
      endif

c.. determine i (chemical potential index)
      i=int((rhn-mutabmin)/dmur)+1
      ii=min(max(i-(m-1)/2,1),jrhomax+1-m)
c      write(0,*)'drdmmx: i, ii ',i,ii ! Debug only

      if(temp.gt.ttabmax) then
         write(0,*) 'from r_drdmmax'
         write(0,*) 'temperature too high! need bigger table'
         write(0,*) 'temp,ttabmax',temp,ttabmax
         temp=ttabmax
      endif

c.. determine j (temperature index)
      j=int((temp-ttabmin)/(ttabmax-ttabmin)
     &     *dble(itmax-1))+1
      jj=min(max(j-(m-1)/2,1),itmax+1-m)
c      write(0,*)'j-(m-1)/2,itmax+1-m',j-(m-1)/2,itmax+1-m ! Debug only
c      write(0,*)'drdmmx: temp, ttabmin, ttabmax',temp,ttabmin,ttabmax ! Debug only      
c      write(0,*)'drdmmx: j, jj ',j,jj ! Debug only      

c.. Determine o (pion chemical potential index)
      o=int(ppt/0.141d0*3.0d0)+1
      oo=min(max(o-(m-1)/2,1),4+1-m)

c.. Determine p (pion chemical potential index)
      p=int(kpt/0.450d0*3.0d0)+1
      pp=min(max(p-(m-1)/2,1),4+1-m)

c      write(0,*)'i,ii,j,jj,o,oo,p,pp',i,ii,j,jj,o,oo,p,pp !Debug only

      res=0.d0

      do o=oo,oo+m-1
       res2(o)=0.d0
       do p=pp,pp+m-1
        res1(p)=0.d0
        do i=ii,ii+m-1
c.. find maximum on the T grid
         maxaux(i)=0.d0
         do j=jj,jj+m-1
            if(meson.eq.1) maxaux(i)=max(maxaux(i),rmax_rho(j,i,o,p))
            if(meson.eq.2) maxaux(i)=max(maxaux(i),rmax_ome(j,i,o,p))
            if(meson.eq.3) maxaux(i)=max(maxaux(i),rmax_phi(j,i,o,p))
         enddo
c.. find maximum on the mu grid
        res1(p)=max(res1(p),maxaux(i))
        enddo
       res2(o)=max(res2(o),res1(p))
       enddo
      res=max(res,res2(o))
      enddo
      
      r_drdmmax=res
c      write(0,*)'r_drdmmax',r_drdmmax ! Debug only

      end

c-----|----------------------------------------------------------------|X

