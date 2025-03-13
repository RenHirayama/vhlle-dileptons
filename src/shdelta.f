c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE SHDELTA
c..
c.. Dalitz decay of the Delta(1232) baryon for "freeze-out" cells (i.e.
c.. with T < 50 MeV or no baryon content 
c*****|****************************************************************|X

      SUBROUTINE shdelta(tau,beta,noe,timestep,mres,p0res,pxres,
     &                  pyres,pzres,factor)
      IMPLICIT NONE

      INCLUDE 'defs.f'

c-----|----------------------------------------------------------------|X 

      integer multi
      integer i,n,timestep,count
      integer ityp,flagdelta

      real*8 mass_lepton,mass_nucleon
      real*8 dwidth,weight,tau,gev,br
      real*8 vacmass_delta,gamma_tot,cv,ph

      real*8 smin_delta,smax_delta,tmax_delta
      real*8 dgamma,mgstar
      real*8 s,t
      real*8 tgamma
      real*8 xgstar,ygstar,zgstar,tgstar,xgamma,ygamma,zgamma
      real*8 p0_gstar,p0_nucleon,p0_pion,p0_gamma
      real*8 p0elab,p0plab  

c... pair & e+/e- Properties
      real*8 mres,p0res,pxres,pyres,pzres
      real*8 p0_el_lab,px_el_lab,py_el_lab,pz_el_lab
      real*8 p0_po_lab,px_po_lab,py_po_lab,pz_po_lab

c... Input/Output
      integer noe
      real*8 beta,gamma,time,factor
      real*8 temp,mub,vol4,acce,accp,lambda

c-----|----------------------------------------------------------------|X 
 
      mass_nucleon=0.938d0
      gev=0.197d0
 
      gamma=1d0/sqrt(1d0-beta**2)

c***************!
      multi=25  !
c***************!

c... DIMUON or DIELECTRON output
      mass_lepton=mass_electron
      if(dimuon) mass_lepton=mass_muon

c... OUTPUT PARAMETERS
      time=timestep*tau

      ityp=17
      flagdelta=7
      vol4=0.0d0
      lambda=0.0d0
      acce=1.0d0
      accp=1.0d0

      mub=0.0d0
      temp=0.0d0
      
c-----|----------------------------------------------------------------|X 

c... Calculation of Delta(1232) Decay Width
      call dgamma_sum_delta(tau,mres,multi,weight,smin_delta,
     &smax_delta,tmax_delta)

      if(weight.eq.0.0d0) return
      
      weight=weight*factor/dble(multi)

      do 999 i=1,multi
      
 111   continue 

       call gamma_star(smin_delta,smax_delta,tmax_delta,s,t,xgstar,
     & ygstar,zgstar,tgstar,xgamma,ygamma,zgamma,tgamma)

c       write(0,*)'smin_delta,smax_delta,tmax_delta',smin_delta,smax_delta,tmax_delta !Debug only
c       write(0,*)'s, t',s,t !Debug only
c       write(0,*)'xgstar,ygstar,zgstar,tgstar',xgstar,ygstar,zgstar,tgstar !Debug only
c       write(0,*)'xgamma,ygamma,zgamma,tgamma',xgamma,ygamma,zgamma,tgamma !Debug only

       call daldelta(tau,s,mres,dgamma)
 
c       write(0,*)'dgamma',dgamma

       if(t.gt.dgamma) goto 111      
     
       mgstar=s
       p0_gstar=mres/2d0-(mass_nucleon**2-mgstar**2)/(2d0*mres)
       p0_nucleon=mres/2d0+(mass_nucleon**2-mgstar**2)/(2d0*mres)

c       write(0,*)'mgstar,p0_gstar,p0_nucleon',mgstar,p0_gstar,p0_nucleon !Debug only
c       write(0,*)'mres,p0res,pxres,pyres,pzres',mres,p0res,pxres,pyres,pzres !Debug only
c       write(0,*)'beta,gamma',beta,gamma !Debug only

       call lobo_dal(p0_gstar,p0_nucleon,mass_nucleon,mgstar,beta,gamma,
     & xgstar,ygstar,zgstar,tgstar,xgamma,ygamma,zgamma,mres,
     & p0res,pxres,pyres,pzres,p0_el_lab,px_el_lab,py_el_lab,pz_el_lab,
     & p0_po_lab,px_po_lab,py_po_lab,pz_po_lab)    
                          
c       write(0,*)'p0_el_lab,px_el_lab,py_el_lab,pz_el_lab',p0_el_lab,px_el_lab,py_el_lab,pz_el_lab !Debug only
c       write(0,*)'p0_po_lab,px_po_lab,py_po_lab,pz_po_lab',p0_po_lab,px_po_lab,py_po_lab,pz_po_lab !Debug only
     
       if(weight.lt.1.0d0.and.weight.gt.0.0d0) then
c... Write into Output File f72
        if(ext_out) then !extended output format
         write(72,556)ityp,weight,mgstar,p0_el_lab,px_el_lab,py_el_lab,
     &   pz_el_lab,p0_po_lab,px_po_lab,py_po_lab,pz_po_lab,tau,time,
     &   vol4,flagdelta,mub,temp,noe,lambda,acce,accp
        else !standard output format
         write(71,555)weight,mgstar,xgstar,ygstar,zgstar,flagdelta,
     &   mub,temp
        endif
       endif
                                                                                                                                 
 999  continue
   
c-----|----------------------------------------------------------------|X 
c format for output 72 
 555  format(e14.7,4f12.7,i9,2f12.7)
 556  format(I5,1X,13(E16.8E3,1X),I4,1X,2(E16.8E3,1X),I9,1X,
     &F12.8,1X,2(E16.8,1x))

      end

c-----|----------------------------------------------------------------|X 
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE DGAMMA_SUM_DELTA
c..
c..
c*****|****************************************************************|X

      SUBROUTINE dgamma_sum_delta(tau,mres,multi,weight_delta,smin_del,
     &smax_del,tmax_del)

      IMPLICIT NONE

      include 'defs.f'

c******************** DELTA(1232) DALITZ DECAY **************************!         

      integer i,multi

c... Mass Resolution
      real*8 dx
      parameter(dx=0.001d0)

c... (Vacuum) Masses
      real*8 vacmass_delta,mass_nucleon 
      real*8 mass_lepton
      parameter(vacmass_delta=1.232d0)
      parameter(mass_nucleon=0.938d0)
          
c... Constants               
      real*8 gamma_tot_const,g
      parameter(gamma_tot_const=0.12d0,g=1.98d0)

      integer z_del
      real*8 mres,e,mt,ml,lambda,gamma0,gamma1,dgamma,gamma_tot
      real*8 weight_delta,dgamma_sum,dgamma_0,dgamma_m  
      real*8 smin_del,smax_del,tmax_del,mx,gev,tau,ph
               
c************************************************************************!
c... This routine uses the description by M.Zetenyi and Gy.Wolf, Heavy   !
c... Ion Physics 17/1 (2003) 27-39 [nucl-th/0202047], which is identical !
c... to Krivoruchenko's approach.                                        ! 
c************************************************************************!

      gev=0.197d0
 
c... DIMUON or DIELECTRON output
      mass_lepton=mass_electron
      if(dimuon) mass_lepton=mass_muon
 
c************************************************************************!
 
      smin_del=2*mass_lepton+0.0000001d0
      smax_del=mres-mass_nucleon

c... For the case of dimuon emission it might happen, that the Delta mass
c... is not sufficient to decay into a muon pair
      if(smax_del.lt.smin_del) then
       weight_delta=0.0d0
       goto 334
      endif

      z_del=int((smax_del-smin_del)/dx)

      dgamma=0.0d0
      dgamma_sum=0.0d0
      weight_delta=0.0d0

c************************************************************************!
      
      do 333 i=1, z_del     
      
      mx = smin_del+dx*i

c-----------------------------------------------------!
c... Phase-Space Factor                               !
      ph=sqrt(1.0d0-(4.0d0*mass_lepton**2/(mx**2)))*  !
     &(1.0d0+(2.0d0*mass_lepton**2/(mx**2)))          !
c-----------------------------------------------------!
                
      e=sqrt(4d0*pi*alpha_em)
      
      mt=(e*g)**2/(12.0d0*mres**2*mass_nucleon**2)*
     &((mres-mass_nucleon)**2-mx**2)*(3.0d0*mres**4+
     &6.0d0*mres**3*mass_nucleon+4.0d0*mres**2*mass_nucleon**2+
     &2.0d0*mres*mass_nucleon**3+mass_nucleon**4-
     &2.0d0*mres*mass_nucleon*mx**2-2.0d0*mass_nucleon**2*mx**2+mx**4)
      
      ml=(e*g)**2*mx**2/(3.0d0*mass_nucleon**2)*
     &((mres-mass_nucleon)**2-mx**2)
      
      lambda=mx**4+mass_nucleon**4+mres**4-2d0*
     & (mx**2*mass_nucleon**2+mx**2*mres**2+mass_nucleon**2*mres**2)


c... To avoid some numerical instabilites (can happen with very rare 
c... combinations of values)

      if(abs(lambda).le.10d-15) then
      lambda=0.0d0
      endif
            
      gamma0=sqrt(lambda)/(16d0*pi*mres**3)*(2d0*mt+ml)

      dgamma=(2d0*alpha_em)/(3d0*pi*mx)*gamma0*(tau/gev)*ph

c... ALTERNATIVE PARAMETRIZATION (KRIVORUCHENKO)..........................
c      gamma0=alpha_em*(mres+mass_nucleon)**2/(16*mres**3
c     &       *mass_nucleon**2)
c     &       *sqrt((mres+mass_nucleon)**2-mx**2)*((mres-mass_nucleon)**2
c     &       -mx**2)**(3/2)*(2.97**2+3*0.03**2)
c
c      gamma1=alpha_em*(mx**2+2*mass_lepton**2)*
c     &       sqrt(1-((4*mass_lepton**2)/(mx**2)))
c
c      dgamma=gamma0*gamma1/(pi*mx**3)*(tau/gev)
c.........................................................................

      if(i.eq.1) then 
      tmax_del=dgamma
      cycle
      endif  
 
      if(tmax_del.lt.dgamma) tmax_del=dgamma
      
      dgamma_sum = dgamma_sum + dgamma
      dgamma=0.0d0
 
  333 continue

c************************************************************************!              
 
      weight_delta=dgamma_sum*dx/multi
c      write(0,*)'weight_delta',weight_delta

  334 continue

      end

c-----|----------------------------------------------------------------|X 
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE GAMMA_STAR
c..
c..
c*****|****************************************************************|X

      SUBROUTINE gamma_star(smin,smax,tmax,s,t,xgstar,ygstar,
     &zgstar,tgstar,xparticle,yparticle,zparticle,tparticle)
      
      IMPLICIT NONE

c*****************************************************************!
      
      real*8 pi,y,phi_gstar,phi_particle
      real*8 theta_particle,sin_theta_gstar,theta_gstar
      real*8 s,smin,smax,t,tmax,x,yy,xparticle,yparticle,zparticle
      real*8 xgstar,ygstar,zgstar,tgstar
      real*8 tparticle

c... Random Function
      real*8 ranff
      external ranff

      integer seedinit
      common /randomf/ seedinit

      parameter(pi=3.141592654d0) 

c*****************************************************************!

      phi_gstar=2d0*pi*ranff(seedinit)
        
      y=ranff(seedinit)
      sin_theta_gstar=(y-0.5)*2d0
      theta_gstar=asin(sin_theta_gstar)+pi/2d0
                              
                              
      phi_particle=phi_gstar + pi
      if (phi_particle.gt.2*pi) then
      phi_particle=phi_particle-2*pi
      endif
                                                      
      if (theta_gstar.le.pi) then
      theta_particle=pi-theta_gstar
      else 
      theta_particle=0
      endif     
      
      xgstar=sin(theta_gstar)*cos(phi_gstar)
      ygstar=sin(theta_gstar)*sin(phi_gstar)
      zgstar=cos(theta_gstar)
      tgstar=sin(theta_gstar)
      
      xparticle=sin(theta_particle)*cos(phi_particle)
      yparticle=sin(theta_particle)*sin(phi_particle)
      zparticle=cos(theta_particle)
      tparticle=sin(theta_particle)
      
        
      t=ranff(seedinit)*tmax
      s=smin+ranff(seedinit)*(smax-smin)

c*****************************************************************!

      end

c-----|----------------------------------------------------------------|X 
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE DALDELTA
c..
c..
c*****|****************************************************************|X

      SUBROUTINE daldelta(tau,mx,mres,dgamma)

      IMPLICIT NONE
          
      include 'defs.f'      

c******************** DELTA(1232) DALITZ DECAY **************************!
c... This routine uses the description by M.Zetenyi and Gy.Wolf, Heavy   !
c... Ion Physics 17/1 (2003) 27-39 [nucl-th/0202047], which is identical !
c... to Krivoruchenko's approach.                                        ! 
c************************************************************************!

      real*8 mass_lepton
      real*8 vacmass_delta,mass_nucleon 
      real*8 mx,gamma_tot_const,g,ph
      real*8 mres,e,mt,ml,lambda,gamma0,gamma1,dgamma
      real*8 dgamma_sum,tau,gev

c... Constants & Masses
      parameter(vacmass_delta=1.232d0)
      parameter(mass_nucleon=0.938d0)

c... Coupling  
      parameter(g=1.98d0)

c... Gamma (Delta --> N+gamma)
      parameter(gamma_tot_const=0.12d0)
      
c************************************************************************!

      gev=0.197d0

c... DIMUON or DIELECTRON output
      mass_lepton=mass_electron
      if(dimuon) mass_lepton=mass_muon

c-----------------------------------------------------!
c... Phase-Space Factor                               !
      ph=sqrt(1.0d0-(4.0d0*mass_lepton**2/(mx**2)))*  !
     &(1.0d0+(2.0d0*mass_lepton**2/(mx**2)))          !
c-----------------------------------------------------!
        
      e=sqrt(4d0*pi*alpha_em)
      
      mt=(e*g)**2/(12.0d0*mres**2*mass_nucleon**2)*
     &((mres-mass_nucleon)**2-mx**2)*(3.0d0*mres**4+
     &6.0d0*mres**3*mass_nucleon+4.0d0*mres**2*mass_nucleon**2+
     &2.0d0*mres*mass_nucleon**3+mass_nucleon**4-
     &2.0d0*mres*mass_nucleon*mx**2-2.0d0*mass_nucleon**2*mx**2+mx**4)
      
      ml=(e*g)**2*mx**2/(3.0d0*mass_nucleon**2)*
     &((mres-mass_nucleon)**2-mx**2)
      
      lambda=mx**4+mass_nucleon**4+mres**4-2d0*
     &(mx**2*mass_nucleon**2+mx**2*mres**2+mass_nucleon**2*mres**2)
      
      gamma0=sqrt(lambda)/(16d0*pi*mres**3)*(2d0*mt+ml)

      dgamma=(2d0*alpha_em)/(3d0*pi*mx)*gamma0*(tau/gev)*ph

c... ALTERNATIVE PARAMETRIZATION (KRIVORUCHENKO)..........................
c      gamma0=alpha_em*(mres+mass_nucleon)**2/(16*mres**3
c     &       *mass_nucleon**2)
c     &       *sqrt((mres+mass_nucleon)**2-mx**2)*((mres-mass_nucleon)**2
c     &       -mx**2)**(3/2)*(2.97**2+3*0.03**2)
c
c      gamma1=alpha_em*(mx**2+2*mass_lepton**2)*
c     &       sqrt(1-((4*mass_lepton**2)/(mx**2)))
c
c      dgamma=gamma0*gamma1/(pi*mx**3)*(tau/gev)
c.........................................................................
c************************************************************************!
           
      end

c-----|----------------------------------------------------------------|X 
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE LOBO_DAL
c..
c..
c*****|****************************************************************|X

      SUBROUTINE lobo_dal(p0_gstar,p0_particle,mass_particle,m_gstar,
     &beta,gamma,x_gstar,y_gstar,z_gstar,t_gstar,x_particle,y_particle,
     &z_particle,m_res,
     &p0_res,px_res,py_res,pz_res,p0_el_lab,px_el_lab,
     &py_el_lab,pz_el_lab,p0_po_lab,px_po_lab,py_po_lab,pz_po_lab)
                                
      IMPLICIT NONE

      include 'defs.f'
 
c**********************************************************************!
     
      real*8 m_gstar,p0_gstar,p0_particle,p_gstar,px_gstar,py_gstar
      real*8 pz_gstar,pt_gstar,p_particle,px_particle,py_particle
      real*8 pz_particle
      real*8 phi_el,x,sin_theta_el,theta_el,phi_po,theta_po
      real*8 p0_el,p0_po,p_el,p_po,px_el,py_el,pz_el,px_po,py_po,pz_po
      real*8 beta_x_gstar,beta_y_gstar,beta_z_gstar,beta_gstar
      real*8 gamma_gstar,p0_el_res,px_el_res,py_el_res,pz_el_res
      real*8 p0_po_res,px_po_res,py_po_res,pz_po_res,beta_x_res
      real*8 beta_y_res,beta_z_res,beta_res,gamma_res,p0_el_cm,px_el_cm
      real*8 py_el_cm,pz_el_cm,p0_po_cm,px_po_cm,py_po_cm,pz_po_cm
      real*8 beta_x_urqmd,beta_y_urqmd,beta_z_urqmd,beta_urqmd
      real*8 gamma_urqmd,p0_el_lab,px_el_lab,py_el_lab,pz_el_lab
      real*8 p0_po_lab,px_po_lab,py_po_lab,pz_po_lab,x_gstar,y_gstar
      real*8 z_gstar,t_gstar,x_particle,y_particle,z_particle,gamma,beta
      real*8 m_res,mass_particle,x_el,y_el,z_el,x_po,y_po,z_po,s
      real*8 p0_res,px_res,py_res,pz_res
      real*8 pt_el_lab,pt_po_lab,pt_di_lab
      real*8 m_total,m_partgstar,m_inv_gstar,m_elpo

c... Random Function
      real*8 ranff
      external ranff

c**********************************************************************!
c... Energies of the gamma* and the particle follow from energy 
c... and momentum conservation in the decay process
c... p0_gstar=m_res/2-(mass_particle**2-m_gstar**2)/(2*m_res)
c... p0_particle=m_res/2+(mass_particle**2-m_gstar**2)/(2*m_res)
c**********************************************************************!

      p_gstar = sqrt(p0_gstar**2 - m_gstar**2)
      px_gstar=x_gstar*p_gstar
      py_gstar=y_gstar*p_gstar
      pz_gstar=z_gstar*p_gstar
      pt_gstar=t_gstar*p_gstar

      p_particle=sqrt(p0_particle**2-mass_particle**2)
      px_particle=x_particle*p_particle
      py_particle=y_particle*p_particle
      pz_particle=z_particle*p_particle

      phi_el=2d0*pi*ranff(seedinit)
      x=ranff(seedinit)                                      
      sin_theta_el=(x-0.5)*2d0                        
      theta_el=asin(sin_theta_el)+pi/2d0

      phi_po=phi_el + pi                            

      if(phi_po.gt.2*pi) then                               
      phi_po=phi_po-2*pi
      end if

      if(theta_el.le.pi) then
      theta_po=pi-theta_el
      else
      theta_po=0
      end if

      x_el=sin(theta_el)*cos(phi_el)               
      y_el=sin(theta_el)*sin(phi_el)
      z_el=cos(theta_el)

      x_po=sin(theta_po)*cos(phi_po)               
      y_po=sin(theta_po)*sin(phi_po) 
      z_po=cos(theta_po)


      p0_el=m_gstar/2d0                               
      p0_po=m_gstar/2d0

c-------------------------------------------!
      if(.not.(dimuon)) then                !
c... DIELECTRONS                            !
      p_el=sqrt(p0_el**2-mass_electron**2)  !      
      p_po=sqrt(p0_po**2-mass_electron**2)  !
      endif                                 ! 
                                            !
      if(dimuon) then                       !
c... DIMUONS                                !
      p_el=sqrt(p0_el**2-mass_muon**2)      !  
      p_po=sqrt(p0_po**2-mass_muon**2)      !
      endif                                 !
c-------------------------------------------!

      px_el=x_el*p_el                              
      py_el=y_el*p_el 
      pz_el=z_el*p_el

      px_po=x_po*p_po                              
      py_po=y_po*p_po
      pz_po=z_po*p_po



      beta_x_gstar=-px_gstar/p0_gstar
      beta_y_gstar=-py_gstar/p0_gstar
      beta_z_gstar=-pz_gstar/p0_gstar

      beta_gstar=sqrt(beta_x_gstar**2+beta_y_gstar**2+beta_z_gstar**2)
      gamma_gstar=1d0/sqrt(1-beta_gstar**2)

      p0_el_res=gamma_gstar*p0_el - beta_x_gstar*gamma_gstar*px_el - 
     &beta_y_gstar*gamma_gstar*py_el -beta_z_gstar*gamma_gstar*pz_el
                                        
      px_el_res=(-beta_x_gstar*gamma_gstar)*p0_el + 
     &(1+(gamma_gstar-1)*((beta_x_gstar*beta_x_gstar)/
     &(beta_gstar*beta_gstar)))*px_el + 
     &(gamma_gstar-1)*((beta_x_gstar*beta_y_gstar)/
     &(beta_gstar*beta_gstar))*py_el + 
     &(gamma_gstar-1)*((beta_x_gstar*beta_z_gstar)/
     &(beta_gstar*beta_gstar))*pz_el

      py_el_res=(-beta_y_gstar*gamma_gstar)*p0_el + 
     &(gamma_gstar-1)*((beta_x_gstar*beta_y_gstar)/
     &(beta_gstar*beta_gstar))*px_el + 
     &(1+(gamma_gstar-1)*((beta_y_gstar*beta_y_gstar)/
     &(beta_gstar*beta_gstar)))*py_el + 
     &(gamma_gstar-1)*((beta_y_gstar*beta_z_gstar)/
     &(beta_gstar*beta_gstar))*pz_el

      pz_el_res=(-beta_z_gstar*gamma_gstar)*p0_el + 
     &(gamma_gstar-1)*((beta_x_gstar*beta_z_gstar)/
     &(beta_gstar*beta_gstar))*px_el +
     &(gamma_gstar-1)*((beta_y_gstar*beta_z_gstar)/
     &(beta_gstar*beta_gstar))*py_el +
     &(1+(gamma_gstar-1)*((beta_z_gstar*beta_z_gstar)/
     &(beta_gstar*beta_gstar)))*pz_el 
     

      p0_po_res=gamma_gstar*p0_po - beta_x_gstar*gamma_gstar*px_po - 
     &beta_y_gstar*gamma_gstar*py_po -beta_z_gstar*gamma_gstar*pz_po

      px_po_res=(-beta_x_gstar*gamma_gstar)*p0_po + 
     &(1+(gamma_gstar-1)*((beta_x_gstar*beta_x_gstar)/
     &(beta_gstar*beta_gstar)))*px_po + 
     &(gamma_gstar-1)*((beta_x_gstar*beta_y_gstar)/
     &(beta_gstar*beta_gstar))*py_po +
     &(gamma_gstar-1)*((beta_x_gstar*beta_z_gstar)/
     &(beta_gstar*beta_gstar))*pz_po

      py_po_res=(-beta_y_gstar*gamma_gstar)*p0_po + 
     &(gamma_gstar-1)*((beta_x_gstar*beta_y_gstar)/
     &(beta_gstar*beta_gstar))*px_po +
     &(1+(gamma_gstar-1)*((beta_y_gstar*beta_y_gstar)/
     &(beta_gstar*beta_gstar)))*py_po + 
     &(gamma_gstar-1)*((beta_y_gstar*beta_z_gstar)/
     &(beta_gstar*beta_gstar))*pz_po 

      pz_po_res=(-beta_z_gstar*gamma_gstar)*p0_po +
     &(gamma_gstar-1)*((beta_x_gstar*beta_z_gstar)/
     &(beta_gstar*beta_gstar))*px_po +
     &(gamma_gstar-1)*((beta_y_gstar*beta_z_gstar)/
     &(beta_gstar*beta_gstar))*py_po +
     &(1+(gamma_gstar-1)*((beta_z_gstar*beta_z_gstar)/
     &(beta_gstar*beta_gstar)))*pz_po 



      beta_x_res=-px_res/p0_res
      beta_y_res=-py_res/p0_res
      beta_z_res=-pz_res/p0_res

      beta_res=sqrt(beta_x_res**2+beta_y_res**2+beta_z_res**2)

      gamma_res=1d0/sqrt(1-beta_res**2)

      p0_el_cm=gamma_res*p0_el_res-beta_x_res*gamma_res*px_el_res-
     &beta_y_res*gamma_res*py_el_res-beta_z_res*gamma_res*pz_el_res

      px_el_cm=(-beta_x_res*gamma_res)*p0_el_res+
     &(1+(gamma_res-1)*((beta_x_res*beta_x_res)/
     &(beta_res*beta_res)))*px_el_res +
     &(gamma_res-1)*((beta_x_res*beta_y_res)/
     &(beta_res*beta_res))*py_el_res +
     &(gamma_res-1)*((beta_x_res*beta_z_res)/
     &(beta_res*beta_res))*pz_el_res

      py_el_cm=(-beta_y_res*gamma_res)*p0_el_res+
     &(gamma_res-1)*((beta_x_res*beta_y_res)/
     &(beta_res*beta_res))*px_el_res+
     &(1+(gamma_res-1)*((beta_y_res*beta_y_res)/
     &(beta_res*beta_res)))*py_el_res+
     &(gamma_res-1)*((beta_y_res*beta_z_res)/
     &(beta_res*beta_res))*pz_el_res 

      pz_el_cm=(-beta_z_res*gamma_res)*p0_el_res+
     &(gamma_res-1)*((beta_x_res*beta_z_res)/
     &(beta_res*beta_res))*px_el_res+
     &(gamma_res-1)*((beta_y_res*beta_z_res)/
     &(beta_res*beta_res))*py_el_res+
     &(1+(gamma_res-1)*((beta_z_res*beta_z_res)/
     &(beta_res*beta_res)))*pz_el_res 


      p0_po_cm=gamma_res*p0_po_res-beta_x_res*gamma_res*px_po_res- 
     &beta_y_res*gamma_res*py_po_res -beta_z_res*gamma_res*pz_po_res

      px_po_cm=(-beta_x_res*gamma_res)*p0_po_res+
     &(1+(gamma_res-1)*((beta_x_res*beta_x_res)/
     &(beta_res*beta_res)))*px_po_res+
     &(gamma_res-1)*((beta_x_res*beta_y_res)/
     &(beta_res*beta_res))*py_po_res+
     &(gamma_res-1)*((beta_x_res*beta_z_res)/
     &(beta_res*beta_res))*pz_po_res

      py_po_cm=(-beta_y_res*gamma_res)*p0_po_res+
     &(gamma_res-1)*((beta_x_res*beta_y_res)/
     &(beta_res*beta_res))*px_po_res+
     &(1+(gamma_res-1)*((beta_y_res*beta_y_res)/
     &(beta_res*beta_res)))*py_po_res+
     &(gamma_res-1)*((beta_y_res*beta_z_res)/
     &(beta_res*beta_res))*pz_po_res 

      pz_po_cm=(-beta_z_res*gamma_res)*p0_po_res+
     &(gamma_res-1)*((beta_x_res*beta_z_res)/
     &(beta_res*beta_res))*px_po_res+
     &(gamma_res-1)*((beta_y_res*beta_z_res)/
     &(beta_res*beta_res))*py_po_res+
     &(1+(gamma_res-1)*((beta_z_res*beta_z_res)/
     &(beta_res*beta_res)))*pz_po_res 


      beta_x_urqmd=0
      beta_y_urqmd=0
      beta_z_urqmd=-beta
      
      
      beta_urqmd=beta
      gamma_urqmd=gamma

      p0_el_lab=gamma_urqmd*p0_el_cm-beta_x_urqmd*gamma_urqmd*px_el_cm- 
     &beta_y_urqmd*gamma_urqmd*py_el_cm-beta_z_urqmd*gamma_urqmd*
     &pz_el_cm

      px_el_lab=(-beta_x_urqmd*gamma)*p0_el_cm + 
     &(1+(gamma_urqmd-1)*((beta_x_urqmd*beta_x_urqmd)/
     &(beta_urqmd*beta_urqmd)))*px_el_cm + 
     &(gamma_urqmd-1)*((beta_x_urqmd*beta_y_urqmd)/
     &(beta_urqmd*beta_urqmd))*py_el_cm + 
     &(gamma_urqmd-1)*((beta_x_urqmd*beta_z_urqmd)/
     &(beta_urqmd*beta_urqmd))*pz_el_cm

      py_el_lab=(-beta_y_urqmd*gamma)*p0_el_cm+ 
     &(gamma_urqmd-1)*((beta_x_urqmd*beta_y_urqmd)/
     &(beta_urqmd*beta_urqmd))*px_el_cm+ 
     &(1+(gamma_urqmd-1)*((beta_y_urqmd*beta_y_urqmd)/
     &(beta_urqmd*beta_urqmd)))*py_el_cm+
     &(gamma_urqmd-1)*((beta_y_urqmd*beta_z_urqmd)/
     &(beta_urqmd*beta_urqmd))*pz_el_cm 

      pz_el_lab=(-beta_z_urqmd*gamma)*p0_el_cm+
     &(gamma_urqmd-1)*((beta_x_urqmd*beta_z_urqmd)/
     &(beta_urqmd*beta_urqmd))*px_el_cm+
     &(gamma_urqmd-1)*((beta_y_urqmd*beta_z_urqmd)/
     &(beta_urqmd*beta_urqmd))*py_el_cm+
     &(1+(gamma_urqmd-1)*((beta_z_urqmd*beta_z_urqmd)/
     &(beta_urqmd*beta_urqmd)))*pz_el_cm 
  
      p0_po_lab=gamma_urqmd*p0_po_cm-beta_x_urqmd*gamma_urqmd*px_po_cm- 
     &beta_y_urqmd*gamma_urqmd*py_po_cm -beta_z_urqmd*gamma_urqmd*
     &pz_po_cm

      px_po_lab=(-beta_x_urqmd*gamma)*p0_po_cm+
     &(1+(gamma_urqmd-1)*((beta_x_urqmd*beta_x_urqmd)/
     &(beta_urqmd*beta_urqmd)))*px_po_cm+ 
     &(gamma_urqmd-1)*((beta_x_urqmd*beta_y_urqmd)/
     &(beta_urqmd*beta_urqmd))*py_po_cm+
     &(gamma_urqmd-1)*((beta_x_urqmd*beta_z_urqmd)/
     &(beta_urqmd*beta_urqmd))*pz_po_cm
     
      py_po_lab=(-beta_y_urqmd*gamma)*p0_po_cm+
     &(gamma_urqmd-1)*((beta_x_urqmd*beta_y_urqmd)/
     &(beta_urqmd*beta_urqmd))*px_po_cm+
     &(1+(gamma_urqmd-1)*((beta_y_urqmd*beta_y_urqmd)/
     &(beta_urqmd*beta_urqmd)))*py_po_cm+
     &(gamma_urqmd-1)*((beta_y_urqmd*beta_z_urqmd)/
     &(beta_urqmd*beta_urqmd))*pz_po_cm 

      pz_po_lab=(-beta_z_urqmd*gamma)*p0_po_cm+
     &(gamma_urqmd-1)*((beta_x_urqmd*beta_z_urqmd)/
     &(beta_urqmd*beta_urqmd))*px_po_cm+
     &(gamma_urqmd-1)*((beta_y_urqmd*beta_z_urqmd)/
     &(beta_urqmd*beta_urqmd))*py_po_cm+
     &(1+(gamma_urqmd-1)*((beta_z_urqmd*beta_z_urqmd)/
     &(beta_urqmd*beta_urqmd)))*pz_po_cm 
     
      m_total=sqrt(p0_res**2-px_res**2-py_res**2-pz_res**2)
      
      m_partgstar=sqrt((p0_gstar+p0_particle)**2-(px_gstar+
     &px_particle)**2-(py_gstar+
     &py_particle)**2-(pz_gstar+pz_particle)**2)
      
      m_inv_gstar=sqrt(p0_gstar**2-px_gstar**2-py_gstar**2-pz_gstar**2)
      
      m_elpo=sqrt(p0_el+p0_po)**2-(px_el+px_po)**2-(py_el+py_po)**2-
     &(pz_el+pz_po)**2        
 
c**********************************************************************!
      
      end
