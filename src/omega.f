      SUBROUTINE diromega(tau,betaLAB,noe,timestep,m_pair,p0_pair,px_pair,
     &                  py_pair,pz_pair,factor)

      IMPLICIT NONE

      include 'defs.f'

c************************ OMEGA DECAYS ***************************!

      integer ityp,flagomega
      integer n,timestep

      real*8 mass_lepton
      real*8 mres,dwidth,br,weight,gev,tau
      real*8 vacmass_omega,gamma_tot,cv,ph
      real*8 br_omega_mumu,br_omega_elpo

      real*8 beta_om, gamma_om, tau_om

c... pair & e+/e- Properties
      real*8 m_pair,p0_pair,px_pair,py_pair,pz_pair
      real*8 p0_el_lab,px_el_lab,py_el_lab,pz_el_lab
      real*8 p0_po_lab,px_po_lab,py_po_lab,pz_po_lab

c... Constants & Masses
      parameter(vacmass_omega=0.78265d0)
      parameter(gamma_tot=0.00849d0)
      parameter(br_omega_mumu=9.00d-5)
      parameter(br_omega_elpo=7.28d-5)

c... Input/Output
      integer noe
      real*8 betaLAB,time,factor
      real*8 temp,mub,vol4,acce,accp,lambda

c*****************************************************************!

      gev=0.197d0 

      mres=m_pair

c... BETA & GAMMA FACTOR & TIME IN THE REST FRAME

      beta_om=sqrt(px_pair**2+py_pair**2+pz_pair**2)/p0_pair
      gamma_om=1.0d0/(sqrt(1.0d0-beta_om**2))
      tau_om=tau/gamma_om

c... DIMUON or DIELECTRON output
      mass_lepton=mass_electron
      if(dimuon) mass_lepton=mass_muon

c... For DIMUONS
      if(dimuon) then
      cv=br_omega_mumu*gamma_tot/vacmass_omega
c... For DIELECTRONS
      else
      cv=br_omega_elpo*gamma_tot/vacmass_omega
      endif

c-----------------------------------------------------!
c... Phase-Space Factor                               !
      ph=sqrt(1.0d0-(4.0d0*mass_lepton**2/mres**2))*  !
     &(1.0d0+(2.0d0*mass_lepton**2/mres**2))          !
c-----------------------------------------------------!

c... Calculation of Gamma(V --> l+l-)  
      dwidth=cv/(mres)**3*(vacmass_omega)**4*ph
    
c... Deermine Weight       
c      br=dwidth*(tau/gev)
      br=dwidth*(tau_om/gev)

      weight=br*factor

c-----------------------------------------------------------------!
c... Determine e+/ei momenta

      call pel_ppo_dist(betaLAB,m_pair,p0_pair,px_pair,py_pair,
     &          pz_pair,p0_el_lab,px_el_lab,py_el_lab,pz_el_lab,
     &          p0_po_lab,px_po_lab,py_po_lab,pz_po_lab)   

c-----------------------------------------------------------------!
c... Write into Output File f72
 
      time=timestep*tau

      ityp=103
      flagomega=7
      vol4=0.0d0
      lambda=0.0d0
      acce=1.0d0
      accp=1.0d0

      mub=0.0d0
      temp=0.0d0
      
      if(ext_out) then !extended output format
       write(72,556)ityp,weight,m_pair,p0_el_lab,px_el_lab,py_el_lab,
     &  pz_el_lab,p0_po_lab,px_po_lab,py_po_lab,pz_po_lab,tau_om,time,vol4,
     &  flagomega,mub,temp,noe,lambda,acce,accp
       else !standard output format
        write(72,555)weight,m_pair,px_pair,py_pair,pz_pair,flagomega,
     &  mub,temp
       endif

c*****************************************************************!
c format for output 72 
 555  format(e14.7,4f12.7,i9,2f12.7)
 556  format(I5,1X,13(E16.8E3,1X),I4,1X,2(E16.8E3,1X),I9,1X,
     &F12.8,1X,2(E16.8,1x))


c*****************************************************************!

      end
                              
c*****************************************************************!
c*****************************************************************!
c*****************************************************************!      
      
      SUBROUTINE dalomega(tau,betaLAB,noe,timestep,m_pair,p0_pair,px_pair,
     &                  py_pair,pz_pair,factor)

      IMPLICIT NONE

      include 'defs.f'

c************************ OMEGA DALITZ DECAYS ********************!

      integer ityp,flagomega,multi
      integer n,timestep,count

      real*8 beta_om, gamma_om, tau_om
     
      real*8 smin,smax,tmax     
      real*8 dgamma,mgstar,tau,time_start,time_end
      real*8 mass_nucleon,mass_gamma,gev
      real*8 p0_gstar,p0_nucleon,p0_pion,p0_gamma

      real*8 gamma_photon,vacmass_omega,ph 
      real*8 t1,t2,f1,lambda_omega,gamma_omega,br
      real*8 mass_lepton

c... pair & e+/e- Properties
      real*8 p0res,pxres,pyres,pzres,mres,ptres,weight 
      real*8 s,t
      real*8 tgamma
      real*8 xgstar,ygstar,zgstar,tgstar,xgamma,ygamma,zgamma

      real*8 m_pair,p0_pair,px_pair,py_pair,pz_pair
      real*8 p0_el_lab,px_el_lab,py_el_lab,pz_el_lab
      real*8 p0_po_lab,px_po_lab,py_po_lab,pz_po_lab

c... Constants & Masses
      parameter(vacmass_omega=0.7826d0)
      parameter(gamma_photon=0.000702972d0)
                               
c... For the form factor
      parameter(lambda_omega=0.65d0)
      parameter(gamma_omega=0.04d0)

c... Branching Ration omega --> pi0+gamma
      parameter(br=0.0828d0)

c... Input/Output
      integer noe
      real*8 betaLAB,gammaLAB,time,factor
      real*8 temp,mub,vol4,acce,accp,lambda

c******************************************************************!
      count=0

      gev=hqc 

      mres=m_pair 
      p0res=p0_pair
      pxres=px_pair
      pyres=py_pair
      pzres=pz_pair

c... BETA & GAMMA FACTOR & TIME IN THE REST FRAME

      beta_om=sqrt(px_pair**2+py_pair**2+pz_pair**2)/p0_pair
      gamma_om=1.0d0/(sqrt(1.0d0-beta_om**2))
      tau_om=tau/gamma_om

c... DIMUON or DIELECTRON output
      mass_lepton=mass_electron
      if(dimuon) mass_lepton=mass_muon

c******************************************************************X
c... MULTI: Multiplier for better statistics for rare resonances
c... Note: Changing the multiplier does not alter the final dilepton
c... yield. it rather serves as a tool to improve statistics for
c... filters and momentum smearing, etc.
c******************!
      multi=25     !
c******************!     
c******************************************************************X

      CALL dgamma_omega(tau_om,multi,weight,smin,smax,tmax,mres,factor)

 103  continue

      CALL gamma_star(smin,smax,tmax,s,t,xgstar,ygstar,zgstar,tgstar,
     &                xgamma,ygamma,zgamma,tgamma)

c-----------------------------------------------------------------!

c-----------------------------------------------------!
c... Phase-Space Factor                               !
      ph=sqrt(1.0d0-(4.0d0*mass_lepton**2/(s**2)))*   !
     &(1.0d0+(2.0d0*mass_lepton**2/(s**2)))           !
c-----------------------------------------------------!

c... Gamma (gamma* --> l+l-)
      t1=((2d0*alpha_em)/(3d0*pi*s))*ph
     
c... Gamma (omega --> pi0+gamma*) [without form factor]
      t2=gamma_photon*(sqrt((1d0+(s**2)/
     &(vacmass_omega**2-mass_pion**2))**2-((2d0*vacmass_omega*s)/
     &(vacmass_omega**2-mass_pion**2))**2))**3 
    
c... Form factor 
      f1=(lambda_omega**2*(lambda_omega**2+gamma_omega**2))/
     &((lambda_omega**2-s**2)**2+(lambda_omega**2*gamma_omega**2))

c... Determine Width       
      dgamma=t1*t2*f1*(tau_om/gev)
c      write(0,*)'dgamma',dgamma !Debug only
c-----------------------------------------------------------------!

      if (t.le.dgamma) then
      mgstar=s

c... Energy of the gamma* and the pion due to energy and  
c... momentum conservation 

      p0_gstar=mres/2d0-(mass_pion**2-mgstar**2)/(2d0*mres)
      p0_pion=mres/2d0+(mass_pion**2-mgstar**2)/(2d0*mres)
c      write(0,*)'p0_gstar,p0_pion',p0_gstar,p0_pion ! Debug only
                  
c... Check if energy is larger than mass (only for omega, for
c... low mass omegas it might happen that p0_gstar < m_gstar       
    
      if(p0_gstar.gt.mgstar) then

c-----------------------------------------------------------------!
c... Determine e+/e- momenta    
      gammaLAB=1.0d0/(sqrt(1.0d0-betaLAB**2))
  
      call lobo_dal(p0_gstar,p0_pion,mass_pion,mgstar,betaLAB,gammaLAB,
     &xgstar,ygstar,zgstar,tgstar,xgamma,ygamma,zgamma,mres,
     &p0res,pxres,pyres,pzres,p0_el_lab,px_el_lab,py_el_lab,pz_el_lab,
     &p0_po_lab,px_po_lab,py_po_lab,pz_po_lab)
     
c      write(0,*)p0_gstar,p0_pion,mass_pion,mgstar,beta_om,gamma_om,
c     &xgstar,ygstar,zgstar,tgstar,xgamma,ygamma,zgamma,mres,
c     &p0res,pxres,pyres,pzres,p0_el_lab,px_el_lab,py_el_lab,pz_el_lab,
c     &p0_po_lab,px_po_lab,py_po_lab,pz_po_lab ! Debug only
c-----------------------------------------------------------------!

      if(weight.gt.0.0d0) then 
      count=count+1         
c-----------------------------------------------------------------!
c... Write into Output File f72
 
      time=timestep*tau

      ityp=103
      flagomega=77
      vol4=0.0d0
      lambda=0.0d0
      acce=1.0d0
      accp=1.0d0

      mub=0.0d0
      temp=0.0d0
      
      if(ext_out) then !extended output format
       write(72,556)ityp,weight,mgstar,p0_el_lab,px_el_lab,py_el_lab,
     &  pz_el_lab,p0_po_lab,px_po_lab,py_po_lab,pz_po_lab,tau_om,time,vol4,
     &  flagomega,mub,temp,noe,lambda,acce,accp
       else !standard output format
        write(72,555)weight,m_pair,px_pair,py_pair,pz_pair,flagomega,
     &  mub,temp
       endif            
c-----------------------------------------------------------------!
      endif

      endif !(p0_gstar.gt.mgstar)
                                                                       
      endif !(t.le.dgamma)

      if(count.lt.multi) goto 103       

      return

c*****************************************************************!
c format for output 72 
 555  format(e14.7,4f12.7,i9,2f12.7)
 556  format(I5,1X,13(E16.8E3,1X),I4,1X,2(E16.8E3,1X),I9,1X,
     &F12.8,1X,2(E16.8,1x))


c*****************************************************************!

      end

c*****************************************************************!
c*****************************************************************!
c*****************************************************************!
      SUBROUTINE dgamma_omega(tau,multi,weight,smin_omega,
     &                      smax_omega,tmax_omega,mres,factor)

      IMPLICIT NONE

      include 'defs.f'
          
c************************************************************************! 

c... Mass Resolution
      real*8 dx
      parameter(dx=0.001d0)

c************************************************************************!
c... The following constants are updated according to the values given   !
c... by J. Beringer et al. (Particle Data Group), Phys. Rev. D86, 010001 !
c... (2012)                                                              !
c************************************************************************! 
      
c...Constants & Masses
      real*8 mass_pi0,mass_lepton
      parameter(mass_pi0=0.1349766d0)

c... omega Dalitz
      real*8 vacmass_omega,dgamma_om,dgamma_sum_omega
      real*8 mx_omega,t1_omega,t2_omega,f1_omega,lambda_omega
      real*8 gamma_omega,br_omega
      integer z_omega
      real*8 tau,gev,gamma_photon,gamma_tot
      real*8 weight
 
      parameter (vacmass_omega=0.7826d0)

c     For the form factor
      parameter(lambda_omega=0.65d0,gamma_omega=0.04d0)

c     Gamma (omega --> pi0 + photon)
      parameter(gamma_photon=0.000702972d0)

c     Gamma (omega --> pi0 + photon) / Gamma_tot
      parameter(br_omega=0.0828d0)
            
      real*8 smin_omega,smax_omega,tmax_omega
           
c... Other      
      integer ityp,i,multi
      real*8 mres,factor

c************************************************************************!
c... This routine is based the descriptions / form factors by            !  
c... L.G. Landsberg, Phys. Rep. 128, 301 (1985) and                      ! 
c... G.Q. Li, C.M. Ko et al., Nucl. Phys. A610, 324C (1996)              !                           
c************************************************************************! 

      dgamma_sum_omega=0.0d0
     
      gev=0.197d0

c... DIMUON or DIELECTRON output
      mass_lepton=mass_electron
      if(dimuon) mass_lepton=mass_muon

c      write(0,*)'logical dimuon =',dimuon
c      write(0,*)'ityp,tau,mres',ityp,tau,mres

c************************************************************************!                      
c... omega

      smin_omega=(2.0d0*mass_lepton)+0.0000001d0
      smax_omega=mres-mass_pi0
      tmax_omega=0.086086d0
              
      z_omega=int((smax_omega-smin_omega)/dx)
      
      do 444 i=1,z_omega

      mx_omega=smin_omega+dx*i

c... gamma* --> mu+ mu-
      t1_omega=((2.0d0*alpha_em)/(3.0d0*pi*mx_omega))*
     &sqrt(1.0d0-(4.0d0*mass_lepton**2)/
     &(mx_omega**2))*(1.0d0+(2.0d0*mass_lepton**2)/(mx_omega**2))
      
c... omega --> pi0 + gamma*
      t2_omega=(sqrt((1.0d0+(mx_omega**2)/(mres**2-
     &mass_pi0**2))**2 
     &-((2.0d0*mres*mx_omega)/(mres**2-mass_pi0**2))**2))**3

c... form factor
      f1_omega=(lambda_omega**2*(lambda_omega**2+gamma_omega**2))/
     &((lambda_omega**2-mx_omega**2)**2+
     &(lambda_omega**2*gamma_omega**2))
      
c... dN(l+l-)/dM
      dgamma_om=t1_omega*t2_omega*f1_omega*gamma_photon*(tau/gev)

c---------------------------------------------------------------------!
c... if no shining shall be applied:
c      dgamma_om=t1_omega*t2_omega*f1_omega*br_omega
c---------------------------------------------------------------------!

c... Summation over all invariant masses      
      dgamma_sum_omega = dgamma_sum_omega + dgamma_om

c... Determine tmax
      if(i.eq.1) then 
      tmax_omega=dgamma_om
      cycle
      endif  
 
      if(tmax_omega.lt.dgamma_om) tmax_omega=dgamma_om

  444 continue      
      
      weight=dgamma_sum_omega*dx/multi*factor
      dgamma_sum_omega=0.0d0      


c************************************************************************!

      end                        
