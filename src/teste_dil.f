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

      PROGRAM teste_dil

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

      character*32 dx_,grd_,vacuum_,eos_,leptype_,rates_,na60mode_,rhonuc_
      character*32 output_,pikchem_,fourpimix_,baryons_,latqgp_,random_
      character*32 phqgp_,hgpar_,model_
      character*96 file71
      character*96 inputFile, inputDir
      real*8 dx,vol
      integer grd,grd_z,vac
      integer leptype

      real*8 time_used,mtest,r_distmhad_hr,dummy,bessk1,mupik
      external r_distmhad_hr,bessk1,mupik
c      integer time
      integer timen
      integer now(1:3),day(1:3)

      real*8 T,mu
      common /th_dyn/T,mu

      real*8 rhonuc,muq,mub,mukaon,mupion,gce,betaLAB,dt,vol4,lambda
      real*8 vxce,vyce,vzce
      integer h, multi

      logical vhlle_finished

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

c.. Random Seed Input
      call getenv('random',random_)
      read(random_,*)  seedinit

c...FOR DILEPTON CALCULATION ONLY.....................................
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

c.. Dielectron / Dimuon emission
      dimuon=.FALSE.
      if(leptype.eq.1) dimuon=.TRUE.

c..Set Input File
      call getenv('inputFile',inputFile)
c..Set Output Units
      call getenv('ftn71',file71)

      if (file71(1:4).ne.'    ') then
         OPEN(UNIT=71,FILE=file71,STATUS='replace',FORM='FORMATTED')
      endif

      write(0,*)'RANDOM SEED',seedinit

c-----|----------------------------------------------------------------|X
c.  Set values for TABLE RESOLUTION


      t=0.15
      mub=0.1
      mupion=0
      mukaon=0
      lambda=0.5
      rhonuc=0.02
      vxce=0
      vyce=0
      vzce=0.2
      multi=1
      betaLAB=0.1
      dt=0.02
      vol4=1
      gce=1/sqrt(1+vxce**2+vyce**2+vzce**2)

      call readeos()
      if(rates.ne.-1.AND.latqgp.eq.1) call qgptables_lat

      itmaxor=itmaxrhom
      jrhomaxor=jrhomaxrhom
      kpimaxor=kpimaxrhom
      lkmaxor=lkmaxrhom
      nqmaxor=nqmaxrhom
      mmmaxor=mmmaxrhom
      rates=4
      call readdil_rapp_rho()
      muq=mub/3
c      write(0,413)h,t,mub,mupion,mukaon,vxce,vyce,vzce,lambda
c      call sleep(1)
       if(mod(h,100).eq.0) write(0,*)h
       if (lambda.gt.0.00001) then
         call qgpemit_lat(lambda,t,muq,gce,vxce,vyce,vzce,multi,vol4,
     &                    betaLAB,dt,h)
       endif
       if (lambda.lt.0.99999) then
       call dilemit_rapp_hr(t,rhonuc,mukaon,mupion,gce,vxce,vyce,vzce,
     &                   vol4,multi,betaLAB,dt,h,lambda)
       call fopiemit_mix(t,mub,rhonuc,mupion,mukaon,gce,vxce,vyce,
     &                 vzce,vol4,multi,betaLAB,dt,h,lambda)
       endif
101   write(0,*)"Did rho, omega, multi-pi, QGP"
        
      itmaxor=itmaxphi
      jrhomaxor=jrhomaxphi
      kpimaxor=kpimaxphi
      lkmaxor=lkmaxphi
      nqmaxor=nqmaxphi
      mmmaxor=mmmaxphi
      rates=3
      call readdil_rapp_phi()
      if(lambda.lt.0.99999) then
        if(mod(h,100).eq.0) write(0,*)h
          call dilemit_rapp_hr(t,rhonuc,mukaon,mupion,gce,vxce,vyce,vzce,
     &                         vol4,multi,betaLAB,dt,h,lambda)
        endif
102   write(0,*)"Did phi"

413   format(I3,7f7.4,1f2.1)
c-----|----------------------------------------------------------------|X
c. **** CALCULATION FINISHED ****

      write(0,*)'COARSE-GRAINING completed successfully',rhonuc

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

