      SUBROUTINE dirrho(tau,betaLAB,noe,timestep,m_pair,p0_pair,px_pair,
     &                  py_pair,pz_pair,factor)
      
      IMPLICIT NONE

      include 'defs.f'

c********************* RHO 0 DIRECT DECAYS ***********************! 

      integer n,timestep,count
      integer ityp,flagrho

      real*8 mass_lepton
      real*8 mres,dwidth,weight,tau,gev,br
      real*8 vacmass_rho,gamma_tot,cv,ph
      real*8 br_rho_mumu,br_rho_elpo

      real*8 beta_rho, gamma_rho, tau_rho

c... pair & e+/e- Properties
      real*8 m_pair,p0_pair,px_pair,py_pair,pz_pair
      real*8 p0_el_lab,px_el_lab,py_el_lab,pz_el_lab
      real*8 p0_po_lab,px_po_lab,py_po_lab,pz_po_lab

c... Constants & Masses
      parameter(vacmass_rho=0.7755d0)
      parameter(gamma_tot=0.1494d0)
      parameter(br_rho_mumu=4.55d-5)
      parameter(br_rho_elpo=4.72d-5)

c... Input/Output
      integer noe
      real*8 betaLAB,time,factor
      real*8 temp,mub,vol4,acce,accp,lambda

c*****************************************************************!

      gev=0.197d0 

      mres=m_pair

c... BETA & GAMMA FACTOR & TIME IN THE REST FRAME

      beta_rho=sqrt(px_pair**2+py_pair**2+pz_pair**2)/p0_pair
      gamma_rho=1.0d0/(sqrt(1.0d0-beta_rho**2))
      tau_rho=tau/gamma_rho

c... DIMUON or DIELECTRON output
      mass_lepton=mass_electron
      if(dimuon) mass_lepton=mass_muon

c... For DIMUONS
      if(dimuon) then
      cv=br_rho_mumu*gamma_tot/vacmass_rho
c... For DIELECTRONS
      else
      cv=br_rho_elpo*gamma_tot/vacmass_rho
      endif          
                        
c     count=count+1

c-----------------------------------------------------!
c... Phase-Space Factor                               !
      ph=sqrt(1.0d0-(4.0d0*mass_lepton**2/mres**2))*  !
     &(1.0d0+(2.0d0*mass_lepton**2/mres**2))          !
c-----------------------------------------------------!

c... Calculation of Gamma(V --> l+l-) 
      dwidth=(cv/(mres)**3)*(vacmass_rho)**4*ph  

c... Deermine Weight      
c      br=dwidth*(tau/gev)
       br=dwidth*(tau_rho/gev)
c      br=1-exp(-(tau/gev)*dwidth)
   
      weight=br*factor

c-----------------------------------------------------------------!
c... Determine e+/ei momenta

      call pel_ppo_dist(betaLAB,m_pair,p0_pair,px_pair,py_pair,
     &          pz_pair,p0_el_lab,px_el_lab,py_el_lab,pz_el_lab,
     &          p0_po_lab,px_po_lab,py_po_lab,pz_po_lab)   

c-----------------------------------------------------------------!
c... Write into Output File f72
 
      time=timestep*tau

      ityp=104
      flagrho=1
      vol4=0.0d0
      lambda=0.0d0
      acce=1.0d0
      accp=1.0d0

      mub=0.0d0
      temp=0.0d0
      
      if(ext_out) then !extended output format
       write(72,556)ityp,weight,m_pair,p0_el_lab,px_el_lab,py_el_lab,
     &  pz_el_lab,p0_po_lab,px_po_lab,py_po_lab,pz_po_lab,tau_rho,time,vol4,
     &  flagrho,mub,temp,noe,lambda,acce,accp
       else !standard output format
        write(72,555)weight,m_pair,px_pair,py_pair,pz_pair,flagrho,
     &  mub,temp
       endif

c*****************************************************************!
c format for output 72 
 555  format(e14.7,4f12.7,i9,2f12.7)
 556  format(I5,1X,13(E16.8E3,1X),I4,1X,2(E16.8E3,1X),I9,1X,
     &F12.8,1X,2(E16.8,1x))

      end
                                                           
