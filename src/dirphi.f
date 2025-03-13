      SUBROUTINE dirphi(tau,betaLAB,noe,timestep,m_pair,p0_pair,px_pair,
     &                  py_pair,pz_pair,factor)

      IMPLICIT NONE   
       
      include 'defs.f'   
      
c********************* PHI DIRECT DECAYS *************************!

      integer ityp,flagphi
      integer n,timestep
            
      real*8 mass_lepton
      real*8 vacmass_phi,gamma_tot,cv,tau,gev
      real*8 mres,dwidth,br,weight,ph
      real*8 br_phi_mumu,br_phi_elpo

      real*8 beta_phi, gamma_phi, tau_phi

c... pair & e+/e- Properties
      real*8 m_pair,p0_pair,px_pair,py_pair,pz_pair
      real*8 p0_el_lab,px_el_lab,py_el_lab,pz_el_lab
      real*8 p0_po_lab,px_po_lab,py_po_lab,pz_po_lab

c... Constants & Masses
      parameter(vacmass_phi=1.01946d0)
      parameter(gamma_tot=0.00426d0)            
      parameter(br_phi_mumu=2.87d-4)
      parameter(br_phi_elpo=2.954d-4)

c... Input/Output
      integer noe
      real*8 betaLAB,time,factor
      real*8 temp,mub,vol4,acce,accp,lambda

c*****************************************************************!

      gev=0.197d0 

      mres=m_pair

c... BETA & GAMMA FACTOR & TIME IN THE REST FRAME

      beta_phi=sqrt(px_pair**2+py_pair**2+pz_pair**2)/p0_pair
      gamma_phi=1.0d0/(sqrt(1.0d0-beta_phi**2))
      tau_phi=tau/gamma_phi

c... DIMUON or DIELECTRON output
      mass_lepton=mass_electron
      if(dimuon) mass_lepton=mass_muon

c... For DIMUONS
      if(dimuon) then
      cv=br_phi_mumu*gamma_tot/vacmass_phi
c... For DIELECTRONS
      else
      cv=br_phi_elpo*gamma_tot/vacmass_phi
      endif

c-----------------------------------------------------!
c... Phase-Space Factor                               !
      ph=sqrt(1.0d0-(4.0d0*mass_lepton**2/mres**2))*  !
     &(1.0d0+(2.0d0*mass_lepton**2/mres**2))          !
c-----------------------------------------------------!

c... Calculation of Gamma(V --> l+l-)                                                            
      dwidth=cv/(mres)**3*(vacmass_phi)**4*ph      

c... Determine Weight          
c      br=dwidth*(tau/gev)	
      br=dwidth*(tau_phi/gev)

      weight=br*factor         

c-----------------------------------------------------------------!
c... Determine e+/ei momenta

      call pel_ppo_dist(betaLAB,m_pair,p0_pair,px_pair,py_pair,
     &          pz_pair,p0_el_lab,px_el_lab,py_el_lab,pz_el_lab,
     &          p0_po_lab,px_po_lab,py_po_lab,pz_po_lab)   

c-----------------------------------------------------------------!
c... Write into Output File f72
 
      time=timestep*tau

      ityp=109
      flagphi=8
      vol4=0.0d0
      lambda=0.0d0
      acce=1.0d0
      accp=1.0d0

      mub=0.0d0
      temp=0.0d0
      
      if(ext_out) then !extended output format
       write(72,556)ityp,weight,mres,p0_el_lab,px_el_lab,py_el_lab,
     &  pz_el_lab,p0_po_lab,px_po_lab,py_po_lab,pz_po_lab,tau_phi,time,vol4,
     &  flagphi,mub,temp,noe,lambda,acce,accp
       else !standard output format
        write(72,555)weight,mres,px_pair,py_pair,pz_pair,flagphi,
     &  mub,temp
       endif

c*****************************************************************!
c format for output 72 
 555  format(e14.7,4f12.7,i9,2f12.7)
 556  format(I5,1X,13(E16.8E3,1X),I4,1X,2(E16.8E3,1X),I9,1X,
     &F12.8,1X,2(E16.8,1x))                          
         
c*****************************************************************!
                                                     
      end
