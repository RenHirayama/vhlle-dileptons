c234567**1*********2*********3*********4*********5*********6*********7**X
c..   ***COARSE-GRAINING PROGRAM FOR URQMD OUTPUT***                    X
c..                                                                     X
c..                                                                     X
c..   The coordinate system                                             X
c..   x = orthogonal to beam, in reaction plane                         X
c..   y = orthogonal to beam, out of reaction plane                     X
c..   z = beam axis                                                     X
c..                                                                     X
c..   Maximal grid size is grid size is (grd*grd*grd_z) cells           X
c..                                                                     X
c..   This program operates in GeV                                      X
c..                                                                     X
c****&|****************************************************************|X
c..                                                                     X
c.. USE RUNFILE TO RUN THIS ROUTINE                                     X
c.. (not all necessary definitions listed, see example runfile)         X
c..                                                                     X  
c.. export vacuum=*spectral function: 1 for vacuum, 0 for in-medium*    X
c.. export gridsize=*gridsize (number of cells per direction)*          X
c.. export cellsize=*length of cells (in fm)*                           X
c.. export eos=*equation of state*                                      X
c.. export leptype=*0 for dielectrons, 1 for dimuons*                   X
c.. ...                                                                 X
c.. ...                                                                 X
c.. export ftn71=*file for regular coarse-graining output*              X
c.. export ftn72=*file for additional rho & Delta shinig output*        X
c..                                                                     X
c..                                                                     X
c..                                                                     X
c.. cat *infile* | ./coarse > *outfile*                                 X
c..                                                                     X
c****&|****************************************************************|X
      
      PROGRAM vhlle_dil

c..ONLY FOR INTEL COMPILER!
cccccccccccccccccccccccccc!
c      use ifport         !
cccccccccccccccccccccccccc!

      implicit none

c.......................................................................
c.***NOTE***: To switch between dimuon/dielectron output and to        . 
c.            determine the equation of state, use the file defs.f     .
c.......................................................................
 
      include 'defs.f'

c. general varibles & grid parameter
      
      character*32 dx_,grd_,vacuum_,eos_,leptype_,rates_,na60mode_
      character*32 output_,pikchem_,fourpimix_,baryons_,latqgp_,random_
      character*32 phqgp_,hgpar_,model_
      character*96 file68,file69,file70,file71,file72
      real*8 dx,vol
      integer grd,grd_z,vac
      integer leptype
      
      real*8 time_used,mtest,distmhad
      external distmhad
c      integer time
      integer timen
      integer now(1:3),day(1:3)

      real*8 T,mu
      common /th_dyn/T,mu

      real*8 rhonuc,muq,mub,mukaon,mupion,gce,betaLAB,dt,vol4,lambda
      real*8 vxce,vyce,vzce
      integer h, multi
      
c.. define general values
 
c..FREEZE-OUT.CALCULATION?..!
c..USE NOT RECOMMENDED......!
c      focalc=.true.        !
      focalc=.false.        !
cccccccccccccccccccccccccccc!

      vac=0   
      eos=6 

      dx=0.0d0
      vol=0.0d0
      grd=0 
      na60mode=0

      outputs=1
      
c.. SPECIFY MODEL (1=UrQMD, 2=SMASH)
      call getenv('model',model_)
      read(model_,*)  model

c.. Vacuum / In-medium calculation
      call getenv('vacuum',vacuum_)
      read(vacuum_,*)  vac
c.. get grid size
      call getenv('gridsize',grd_)
      read(grd_,*)  grd
c.. get value of dx
      call getenv('cellsize',dx_)
      read(dx_,*)  dx

c.. Determination of Equation of State:
c. *eos=1* for ULTRARELATIVISTIC GAS
c. *eos=2* for HADRON GAS EoS
c. *eos=3* for CHIRAL EoS
c. *eos=4* for HG + QGP above T_c
c. *eos=5* for CHIRAL EoS + QGP above T_c
c. *eos=6* for HG EoS + QGP from LATTICE EoS above T_c

      call getenv('eos',eos_)
      read(eos_,*)  eos

c.. Determine Rates (DILEPTON: 0-3, PHOTON: -1)
      call getenv('rates',rates_)
      read(rates_,*)  rates

c.. Pi/K Chemical Potentials (Yes: *pikchem=1*, No (mu=0): *pikchem=0*)
      call getenv('pikchem',pikchem_)
      read(pikchem_,*)  pikchem

c.. Baryon Chemical Potential (Yes: *baryons=1*, No (mu=0): *baryons=0*)
      call getenv('baryons',baryons_)
      read(baryons_,*)  baryons

c.. Suppress Output in f71 an f72 for *output=0*
      call getenv('output',output_)
      read(output_,*)  outputs

c.. Random Seed Input
      call getenv('random',random_)
      read(random_,*)  seedinit

c...FOR DILEPTON CALCULATION ONLY.....................................
      if (rates.ge.0) then
c.. Dimuon *OR* Dielectron Emission
      call getenv('leptype',leptype_) 
      read(leptype_,*)  leptype

c.. 4-Pion Rates (V-A Mixing: *fourpimix=1*, Pure V: *fourpimix=0*)
      call getenv('fourpimix',fourpimix_)
      read(fourpimix_,*)  fourpimix

c.. Lattice QGP rates (Yes: *latqgp=1*, No: *latqgp=0*)
      call getenv('latqgp',latqgp_)
      read(latqgp_,*)  latqgp

c.. dNch/deta > 30 selection for NA60
      call getenv('na60mode',na60mode_)
      read(na60mode_,*)  na60mode

c...FOR PHOTON CALCULATION ONLY.....................................
      elseif (rates.lt.0) then

c.. Full QGP parametrization (Yes: *phqgp=1*, No: *phqgp=0*)
      call getenv('phqgp',phqgp_)
      read(phqgp_,*)  phqgp       

c.. Parametrization for HG-Contribution (1: Rapp/Heffernan/Turbide, 2: Kapusta/Lichard)
      call getenv('hgpar',hgpar_)
      read(hgpar_,*)  hgpar       

      endif
c...................................................................

      if (hgpar.ne.1.AND.rates.lt.0) then
      stop 'Alternative HG-Rates not yet implemented!'
      endif

      if (vac.eq.1.AND.rates.ge.1) rates=0

      if(dx.eq.0.0d0) stop 'Must define cellsize!'
      if(grd.eq.0) stop 'Must define gridsize!'
      if(eos.eq.0) stop 'Must define equation of state!'
      if(abs(latqgp).gt.1) stop 'Wrong parameter for QGP rates!'
      if(abs(fourpimix).gt.1) stop 'Wrong parameter for 4-pi-mixing!'

c.. define grid size in beam direction
      grd_z=3*grd

c.. volume of cell
      vol=dx*dx*dx 

c.. Dielectron / Dimuon emission
      dimuon=.FALSE.
      if(leptype.eq.1) dimuon=.TRUE. 
 
c..Set Output Units
      call getenv('ftn68',file68)
      call getenv('ftn69',file69)
      call getenv('ftn70',file70)
      call getenv('ftn71',file71)
      call getenv('ftn72',file72)

      if (file68(1:4).ne.'    ') then
         OPEN(UNIT=68,FILE=file68,STATUS='replace',FORM='FORMATTED')
      endif
      if (file69(1:4).ne.'    ') then
         OPEN(UNIT=69,FILE=file69,STATUS='replace',FORM='FORMATTED')
      endif
      if (file70(1:4).ne.'    ') then
         OPEN(UNIT=70,FILE=file70,STATUS='replace',FORM='FORMATTED')
      endif
      if (file71(1:4).ne.'    ') then
         OPEN(UNIT=71,FILE=file71,STATUS='replace',FORM='FORMATTED')
      endif
      if (file72(1:4).ne.'    ') then
         OPEN(UNIT=72,FILE=file72,STATUS='replace',FORM='FORMATTED')
      endif

c.. write the defined quantities
      write(0,*)'****************************************************'
      write(0,*)'*** COARSE-GRAINING PROGRAM FOR TRANSPORT INPUT  ***'  
      write(0,*)'****************************************************'
      write(0,*)'***  This program reads UrQMD or SMASH output    ***'      
      write(0,*)'***  and puts the particles on a grid of small   ***'
      write(0,*)'***  space-time cells to determine temperature   ***'  
      write(0,*)'***  and chem. potential, and calculates the     ***'  
      write(0,*)'***  corresponding thermal dilepton or photon    ***'
      write(0,*)'***  emission.                                   ***'
      write(0,*)'***                                              ***'
      write(0,*)'***  When using for calculations, please cite:   ***'
      write(0,*)'***  -> Phys.Rev. C 91 (2015) 054911 -> NA60     ***'
      write(0,*)'***  -> Phys.Rev. C 92 (2015) 014911 -> HADES    ***'
      write(0,*)'***  -> Phys.Rev. C 93 (2016) 054901 -> FAIR     ***'
      write(0,*)'***  -> Phys.Rev. C 94 (2016) 024912 -> RHIC/LHC ***'
      write(0,*)'****************************************************'
      if (vac.eq.1) write(0,*)'---VACUUUM CALCULATION---'
      if (vac.eq.0) write(0,*)'---IN-MEDIUM CALCULATION---'
c... Determine the the lepton-type 
      if(rates.ge.0) then
       if(dimuon) then
       write(0,*)'---of DIMUON output---'
       else !dielectron
       write(0,*)'---of DIELECTRON output---'
       endif
      else if(rates.eq.-1) then
       write(0,*)'---of PHOTON output---'
      endif
      if(eos.eq.1) write(0,*)'---with ULTRAS-REL. GAS EoS---'
      if(eos.eq.2) write(0,*)'---with HADRON GAS EoS---'
      if(eos.eq.3) write(0,*)'---with CHIRAL EoS---'
      if(eos.eq.4) write(0,*)'---with HG Eos + QGP above T_c---'
      if(eos.eq.5) write(0,*)'---with CHIRAL EoS + QGP above T_c---'
      if(eos.eq.6) write(0,*)'---with HADRON EoS (LATTICE EoS around T_c)---'
      if(rates.eq.-1) then
      write(0,*)'**IM-SPECTRAL FUNCTION FOR PHOTONS**'
       if(eos.ge.3) then
       write(0,*)'***  (+ BREMSSTRAHLUNG AND QGP)  ***'
        if(phqgp.eq.1) then
        write(0,*)'**** FULL QGP PARAMETRIZATION  ****'
        elseif(phqgp.eq.0) then
        write(0,*)'********  pQCD QGP-Rates  *********'
        endif
       else
       write(0,*)'******  (+ BREMSSTRAHLUNG)  ********'
       endif
      endif
      if(rates.eq.0) write(0,*)'***ELETSKY SPECTRAL FUNCTION***'
      if(rates.eq.1) then
      write(0,*)'***OLD RAPP/WAMBACH SPECTRAL FUNCTION***'
      write(0,*)'***FOR RHO, OMEGA AND PHI RATES***'
      write(0,*)'***DO NOT USE ANY LONGER !!! ***' 
      endif
      if(rates.eq.2) then
      write(0,*)'***OLD RAPP/WAMBACH SPECTRAL FUNCTION***'
      write(0,*)'***RHO/OMEGA ONLY & HIGH MASS RESOLUTION***'
      write(0,*)'***DO NOT USE ANY LONGER !!! ***' 
      endif
      if(rates.eq.3) then
      write(0,*)'***RAPP/WAMBACH SPECTRAL FUNCTION***'
      write(0,*)'***PHI IN HIGH MASS RESOLUTION ONLY***' 
      write(0,*)'***NO QGP / MULTI_PION EMISSION!!!***'
      endif
      if(rates.eq.4) then
      write(0,*)'***UPDATED RAPP/WAMBACH SF***'
      write(0,*)'***RHO/OMEGA ONLY & HIGH MASS RESOLUTION***' 
      endif
      if(rates.ne.3.AND.rates.ne.-1) then
       if(htlcorr.AND.(eos.eq.3.OR.eos.eq.4.OR.eos.eq.5.OR.eos.eq.6).AND.latqgp.eq.0) then
       write(0,*)'***HTL-improved QGP rates***'
       endif
       if((eos.eq.3.OR.eos.eq.4.OR.eos.eq.5.OR.eos.eq.6).AND.latqgp.eq.1) then
       write(0,*)'***Lattice QGP rates***'
       endif
       if(fourpimix.eq.1) then
       write(0,*)'***V-A Mixing for 4-pion Rates***'
       endif
      endif
      if(baryons.eq.0) write(0,*)'***NO BARYONIC EFFECTS IN CALCULATION***'
      write(0,*)'*************************************************'
      write(0,*)'started with:'
      write(0,'(A23,F5.2)')'--> Cell length (fm) = ',dx
      write(0,'(A26,F12.4)')'--> Cell volume (fm**3) = ',vol
      write(0,'(A22,I3)')'--> Grid size (x,y) = ',grd
      write(0,'(A20,I3)')'--> Grid size (z) = ',grd_z
      write(0,*)'*************************************************'
      if(na60mode.eq.1) then
      write(0,*)'******CALCULATION FOR NA60: dN_ch/deta > 30******'
      write(0,*)'*************************************************'
      endif
c-----|----------------------------------------------------------------|X

c.. Initialize the RANDOM SEED
c      call IDATE(day)
c      call ITIME(now)
c      call SYSTEM_CLOCK(timen)
c      call SYSTEM_CLOCK(count=timen)

c      seedinit=abs((now(1)+1)*(now(2)+1)*(now(3)+1)+(now(3)+1)+
c     &         day(1)*day(2)*day(3)*timen*(now(3)+55))

c      write(0,*)'time, day',now(1),now(2),now(3),day(1),day(2),day(3) ! Debug only

      write(0,*)'RANDOM SEED',seedinit
  
c-----|----------------------------------------------------------------|X
c.  Set values for TABLE RESOLUTION

        itmaxor=0
        jrhomaxor=0
        kpimaxor=0
        lkmaxor=0
        nqmaxor=0
        mmmaxor=0

      if(rates.eq.2.OR.rates.eq.4) then
        itmaxor=itmaxrhom
        jrhomaxor=jrhomaxrhom
        kpimaxor=kpimaxrhom
        lkmaxor=lkmaxrhom
        nqmaxor=nqmaxrhom
        mmmaxor=mmmaxrhom
      else if(rates.eq.3) then
        itmaxor=itmaxphi
        jrhomaxor=jrhomaxphi
        kpimaxor=kpimaxphi
        lkmaxor=lkmaxphi
        nqmaxor=nqmaxphi
        mmmaxor=mmmaxphi
      endif 

c-----|----------------------------------------------------------------|X
c.  Read in EoS file and tabulated dilepton rates
      
      call readeos()
      

      if(rates.eq.1) then
        call readdil_rapp()
      elseif(rates.eq.2.OR.rates.eq.4) then
        call readdil_rapp_rho()
      elseif(rates.eq.3) then
        call readdil_rapp_phi()
      elseif(rates.eq.0) then
        call readdil_elet(vac)
      endif
   
      if(rates.ne.-1.AND.latqgp.eq.1) call qgptables_lat

c-----|----------------------------------------------------------------|X
c
c.. FOR DEBUGGING: This is just to check the results for distmhad(T,mu,m) 
c      mtest=0.001022d0
c      T=0.050d0
c      mu=0.050d0
c      do while (mtest.lt.1.1d0)
c      write (72,*),mtest,distmhad(mtest)
c      mtest = mtest+5.0d-3
c      end do
c  
c-----|----------------------------------------------------------------|X
      
c.  now read in input from UrQMD file14 and do what has to be done
      t=0.15
      rhonuc=0.2
      muq=0.3
      mukaon=0
      mupion=0
      gce=0.1
      vxce=0
      vyce=0
      vzce=0.3
      vol4=1
      multi=100
      betaLAB=0.3
      dt=1
      h=1
      lambda=0.15
      do h=1,1000
       if(mod(h,100).eq.0) write(0,*)h
       call qgpemit_lat(lambda,t,muq,gce,vxce,vyce,vzce,multi,vol4,
     &                  betaLAB,dt,h)
       call dilemit_rapp_hr(t,rhonuc,mukaon,mupion,gce,vxce,vyce,vzce,
     &                     vol4,multi,betaLAB,dt,h,lambda)
       call fopiemit_mix(t,mub,rhonuc,mupion,mukaon,gce,vxce,vyce,
     &                   vzce,vol4,multi,betaLAB,dt,h,lambda)
      end do

c      call vhlle_read(grd,grd_z,dx,vol,vac)
      
c-----|----------------------------------------------------------------|X
c. **** CALCULATION FINISHED ****
      
      write(0,*)'COARSE-GRAINING completed successfully'

c-----|----------------------------------------------------------------|X      
c. **** TIME ELAPSED ****

      call cpu_time(time_used)
      write(0,*)'Elapsed CPU time =',time_used  
 
      close(68)
      close(69)
      close(70)
      close(71)
      close(72)

      stop '****Calculation finished****'

c-----|----------------------------------------------------------------|X 
      
      end       
      
