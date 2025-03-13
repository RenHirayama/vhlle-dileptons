c************************************************************8
c.. Authors: Sascha Vogel, Elvira Santini, Katharina Schmidt,
c..          Marcus Bleicher
c.. Dilepton routine for UrQMD output
c************************************************************8
 
      subroutine dilep(ityp,chg,p0res,pxres,pyres,pzres,mres,tau,
     &time_start,time_end,flag,dens_cre,dens_abs,beta,noe,bev)
      implicit none
      
cyyyyy     cccccc resonance identifiers (as from UrQMD) ccccc
c     

c     omega = 103
c     pion=101
c     rho=104
c     phi=109
c     eta =102
c     etaprime=107
c     delta1232=17

                                                  

      integer omega,rho,phi,pion,etaprime,eta,delta,omegadir,omegadal
      integer ityp,chg,multi,n,validm,flag
   
c     ccccc variables ccccc
      
      real*8 dens_cre,dens_abs
      real*8 p0res,pxres,pyres,pzres,mres,ptres,weight
      real*8 pxelab,pyelab,pzelab,pxplab,pyplab,pzplab
      real*8 p0elab,p0plab     
      real*8 s,t
      real*8 tgamma
      real*8 xgstar,ygstar,zgstar,tgstar,xgamma,ygamma,zgamma
      real*8 weight_pion,weight_omega,weight_etaprime,weight_eta
      real*8 weight_delta
      
      real*8 smin_pion,smax_pion,tmax_pion
      real*8 smin_omega,smax_omega,tmax_omega
      real*8 smin_etaprime,smax_etaprime,tmax_etaprime
      real*8 smin_eta,smax_eta,tmax_eta
      real*8 smin_delta,smax_delta,tmax_delta
      
      real*8 dgamma,mgstar,tau,time_start,time_end
      real*8 mass_nucleon,mass_pion,mass_gamma
      real*8 p0_gstar,p0_nucleon,p0_pion,p0_gamma

c     ccccc data input ccccc

      character*9 name1
      character*9 name2(3)
      character*1 help

      integer event
      integer noe
      real*8 bev
      real*8 beta,gamma
      real*8 list(15)

      integer anpion

c     multi: multiplier for better statistics for rare resonances
c     note: changing the multiplier does not alter the final dilepton
c     yield. it rather serves as a tool to improve statistics for filters
c     and momentum smearing etc.
c
      multi=1
     

c     ccccc UrQMD columns and particle identifiers ccccc
          
c      ityp=10
c      chg=12
c      omega direct and omega dalitz will are labeled 1031 and 1032      
      
      omega=103
      rho=104
      phi=109
      pion=101
      etaprime=107
      eta=102
      delta=17
      omegadir=1031
      omegadal=1032
      
                                
      gamma=1d0/sqrt(1d0-beta**2)

c     ccccc calculation of decay widths ccccc
      call dgamma_sum(ityp,tau,multi,weight_pion,weight_omega,
     &weight_etaprime,
     &weight_eta,smin_pion,smax_pion,tmax_pion,smin_omega,smax_omega,
     &tmax_omega,smin_etaprime,smax_etaprime,tmax_etaprime,smin_eta,
     &smax_eta,tmax_eta)


      validm=0

c     ccccc direct and dalitz decay of the omega meson ccccc

      if (ityp.eq.omega) then

      do 222 n=1,multi   
      call diromega(tau,mres,multi,weight) 
       
      call lobo_dir(beta,gamma,mres,p0res,pxres,pyres,pzres,ptres,
     &p0elab,pxelab,pyelab,pzelab,p0plab,pxplab,pyplab,pzplab)

      call output (omegadir,weight,mres,pxres,pyres,pzres,
     &p0elab,pxelab,pyelab,pzelab,p0plab,pxplab,pyplab,pzplab,multi,
     &tau,time_start,time_end,flag,dens_cre,dens_abs,noe,bev)


  222 continue
  
  666 if(validm.lt.multi) then
      call t_omega(tau,multi,
     &smin_omega,smax_omega,tmax_omega)
           
      call gamma_star(smin_omega,smax_omega,tmax_omega,s,t,xgstar,
     &ygstar,zgstar,tgstar,xgamma,ygamma,zgamma,tgamma)

      call dalomega(tau,s,dgamma,mass_pion)
                     
      if (t.le.dgamma) then
      mass_pion=0.135d0
      mgstar=s

c     energy of the gamma* and the pion due to energy and  
c     momentum conservation 

      p0_gstar=mres/2d0-(mass_pion**2-mgstar**2)/(2d0*mres)
      p0_pion=mres/2d0+(mass_pion**2-mgstar**2)/(2d0*mres)
                  
c     check if energy is larger than mass (only for omega, for
c     low mass omegas it might happen that p0_gstar < m_gstar       
    
      if(p0_gstar.gt.mgstar) then
      validm=validm+1
    
      call lobo_dal(p0_gstar,p0_pion,mass_pion,mgstar,beta,gamma,
     &xgstar,ygstar,zgstar,tgstar,xgamma,ygamma,zgamma,mres,
     &p0res,pxres,pyres,pzres,p0elab,pxelab,pyelab,pzelab,p0plab,
     &pxplab,pyplab,pzplab)
                                                      
      call output (omegadal,weight_omega,mgstar,pxres,pyres,pzres,
     &p0elab,pxelab,pyelab,pzelab,p0plab,pxplab,pyplab,pzplab,
     &multi,tau,time_start,time_end,flag,dens_cre,dens_abs,noe,bev)
      endif
      endif
                                                                       
      goto 666
      endif
      validm=0
                                                                                               
c     ccccc direct decay of the rho meson ccccc
      
      elseif (ityp.eq. rho .and. chg .eq. 0) then
      do 333 n=1,multi
      call dirrho(tau,mres,multi,weight)
      
      call lobo_dir(beta,gamma,mres,p0res,pxres,pyres,pzres,ptres,
     &p0elab,pxelab,pyelab,pzelab,p0plab,pxplab,pyplab,pzplab)
  
      call output (rho,weight,mres,pxres,pyres,pzres,p0elab,
     &pxelab,pyelab,pzelab,p0plab,pxplab,pyplab,pzplab,multi,
     &tau,time_start,time_end,flag,dens_cre,dens_abs,noe,bev)

     
  333 continue
c      goto 10
     

c     ccccc direct decay of the phi meson ccccc
      elseif (ityp.eq.phi) then
      do 444 n=1,multi
      call dirphi(tau,mres,multi,weight)
            
      call lobo_dir(beta,gamma,mres,p0res,pxres,pyres,pzres,ptres,
     &p0elab,pxelab,pyelab,pzelab,p0plab,pxplab,pyplab,pzplab)

      call output (phi,weight,mres,pxres,pyres,pzres,
     &p0elab,pxelab,pyelab,pzelab,p0plab,pxplab,pyplab,pzplab,multi,
     &tau,time_start,time_end,flag,dens_cre,dens_abs,noe,bev)
                                      
  444 continue
c      goto 10


c     ccccc dalitz decay of the pi0 meson ccccc
      elseif (ityp.eq.pion.and.chg.eq.0) then
      anpion=anpion+1        
  555 if(validm.lt.multi)then 
      call gamma_star(smin_pion,smax_pion,tmax_pion,s,t,xgstar,ygstar,
     &zgstar,tgstar,xgamma,ygamma,zgamma,tgamma)
          
      call dalpi(s,dgamma)

      if (t.le.dgamma) then
      mass_gamma=0
      mgstar=s

      p0_gstar=mres/2-(mass_gamma**2-mgstar**2)/(2*mres)
      p0_gamma=mres/2+(mass_gamma**2-mgstar**2)/(2*mres)

      validm=validm+1
      call lobo_dal(p0_gstar,p0_gamma,mass_gamma,mgstar,beta,gamma,
     &xgstar,ygstar,zgstar,tgstar,xgamma,ygamma,zgamma,mres,
     &p0res,pxres,pyres,pzres,p0elab,pxelab,pyelab,pzelab,p0plab,
     &pxplab,pyplab,pzplab)
                            
      call output (pion,weight_pion,mgstar,pxres,pyres,pzres,
     &p0elab,pxelab,pyelab,pzelab,p0plab,pxplab,pyplab,pzplab,multi,
     &tau,time_start,time_end,flag,dens_cre,dens_abs,noe,bev)                                 
      endif
                                                        
      goto 555
      endif
      validm=0
c      goto 10
      
                                                               
c     ccccc dalitz decay of the eta prime meson ccccc
      elseif (ityp.eq.etaprime) then
               
                                
  777 if (validm.lt.multi) then
      call gamma_star(smin_etaprime,smax_etaprime,tmax_etaprime,s,t,
     &xgstar,ygstar,zgstar,tgstar,xgamma,ygamma,zgamma,tgamma)
                                 
      call daletaprime(s,dgamma)
             
      if (t.le.dgamma) then
      mass_gamma=0
      mgstar=s

      p0_gstar=mres/2d0-(mass_gamma**2-mgstar**2)/(2d0*mres)
      p0_gamma=mres/2d0+(mass_gamma**2-mgstar**2)/(2d0*mres)

      validm=validm+1                                                             
      call lobo_dal(p0_gstar,p0_gamma,mass_gamma,mgstar,beta,gamma,
     &xgstar,ygstar,zgstar,tgstar,xgamma,ygamma,zgamma,mres,
     &p0res,pxres,pyres,pzres,p0elab,pxelab,pyelab,pzelab,p0plab,
     &pxplab,pyplab,pzplab)
                                                                                    
      call output (etaprime,weight_etaprime,mgstar,pxres,pyres,pzres,
     &p0elab,pxelab,pyelab,pzelab,p0plab,pxplab,pyplab,pzplab,multi,
     &tau,time_start,time_end,flag,dens_cre,dens_abs,noe,bev)
      endif
                                                            
      goto 777
      endif
      validm=0                                                                                                   
c      goto 10
                                                                                                                           
c     ccccc dalitz decay of the eta meson ccccc

      elseif (ityp.eq.eta) then
              
  888 if (validm.lt.multi) then
      call gamma_star(smin_eta,smax_eta,tmax_eta,s,t,xgstar,ygstar,
     &zgstar,tgstar,xgamma,ygamma,zgamma,tgamma)
                                                   
      call daleta(s,dgamma)
      
      if (t.le.dgamma) then
      mass_gamma=0
      mgstar=s
      p0_gstar=mres/2d0-(mass_gamma**2-mgstar**2)/(2d0*mres)
      p0_gamma=mres/2d0+(mass_gamma**2-mgstar**2)/(2d0*mres)

      validm=validm+1                                                       
      call lobo_dal(p0_gstar,p0_gamma,mass_gamma,mgstar,beta,gamma,
     &xgstar,ygstar,zgstar,tgstar,xgamma,ygamma,zgamma,mres,
     &p0res,pxres,pyres,pzres,p0elab,pxelab,pyelab,pzelab,p0plab,
     &pxplab,pyplab,pzplab)
      
                                                                            
      call output (eta,weight_eta,mgstar,pxres,pyres,pzres,
     &p0elab,pxelab,pyelab,pzelab,p0plab,pxplab,pyplab,pzplab,multi,
     &tau,time_start,time_end,flag,dens_cre,dens_abs,noe,bev)
      endif
                                                                                                   
      goto 888
      endif  
      
      validm=0
      


c     ccccc dalitz decay of the Delta(1232) baryon ccccc
      elseif((ityp.eq.delta.and.chg.eq.0).or.
     &(ityp.eq.delta.and.chg.eq.1)) then
c      write(0,*)'mres_delta',mres
      call dgamma_sum_delta(tau,mres,multi,weight_delta,smin_delta,
     &smax_delta,tmax_delta)
      
  999 if (validm.lt.multi) then

      call t_delta(tau,mres,multi,tmax_delta)      
      call gamma_star(smin_delta,smax_delta,tmax_delta,s,t,xgstar,
     &ygstar,zgstar,tgstar,xgamma,ygamma,zgamma,tgamma)

      call daldelta(tau,s,mres,dgamma,mass_nucleon)
           
      if (t.le.dgamma) then
      mgstar=s
      p0_gstar=mres/2d0-(mass_nucleon**2-mgstar**2)/(2d0*mres)
      p0_nucleon=mres/2d0+(mass_nucleon**2-mgstar**2)/(2d0*mres)
      
      validm=validm+1
      call lobo_dal(p0_gstar,p0_nucleon,mass_nucleon,mgstar,beta,gamma,
     &xgstar,ygstar,zgstar,tgstar,xgamma,ygamma,zgamma,mres,
     &p0res,pxres,pyres,pzres,p0elab,pxelab,pyelab,pzelab,p0plab,
     &pxplab,pyplab,pzplab)    
                               
              
      if(weight_delta.lt.1.0d0) then
      call output (delta,weight_delta,mgstar,pxres,pyres,pzres,
     &p0elab,pxelab,pyelab,pzelab,p0plab,pxplab,pyplab,pzplab,multi,
     &tau,time_start,time_end,flag,dens_cre,dens_abs,noe,bev)
      endif
      if(weight_delta.gt.1.0d0) then
      write(0,*) 'suspicious Delta weight',weight_delta
      endif      

      endif
                                                         
      goto 999
      endif
                                                                           
      endif
      validm=0
      weight_pion=0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end       
