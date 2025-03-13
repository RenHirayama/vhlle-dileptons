      SUBROUTINE dgamma_omega(tau,multi,weight_omega,smin_omega,
     &                      smax_omega,tmax_omega,mres)

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
      real*8 vacmass_omega,dgamma_omega,dgamma_sum_omega
      real*8 mx_omega,t1_omega,t2_omega,f1_omega,lambda_omega
      real*8 gamma_omega,br_omega
      integer z_omega
      real*8 tau,gev,gamma_photon,gamma_tot

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
      real*8 mres

c************************************************************************!
c... This routine is based the descriptions / form factors by            !  
c... L.G. Landsberg, Phys. Rep. 128, 301 (1985) and                      ! 
c... G.Q. Li, C.M. Ko et al., Nucl. Phys. A610, 324C (1996)              !                           
c************************************************************************! 

      dgamma_sum_omega=0.0d0
     
      gev=hqc

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
      dgamma_omega=t1_omega*t2_omega*f1_omega*gamma_photon*(tau/gev)

c---------------------------------------------------------------------!
c... if no shining shall be applied:
c      dgamma_omega=t1_omega*t2_omega*f1_omega*br_omega
c---------------------------------------------------------------------!

c... Summation over all invariant masses      
      dgamma_sum_omega = dgamma_sum_omega + dgamma_omega

c... Determine tmax
      if(i.eq.1) then 
      tmax_omega=dgamma_omega
      cycle
      endif  
 
      if(tmax_omega.lt.dgamma_omega) tmax_omega=dgamma_omega

  444 continue      
      
      weight_omega=dgamma_sum_omega*dx/multi
      dgamma_sum_omega=0.0d0      


c************************************************************************!

      end                                          
                                                                          
