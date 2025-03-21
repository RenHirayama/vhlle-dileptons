c....&|...Version 1.0...last change...25.09.2012...Stephan.Endres......|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE READEOS
c..
c.. Subroutine which reads in EoS matrices. Readeos reads in pure 
c.. hadronic EoS.
c*****|****************************************************************|X
 
      SUBROUTINE readeos()
      implicit none
      
      include 'defs.f'

c-----|----------------------------------------------------------------|X

c. Variables
      real*8 t,mu,e,p,n,s,mstar,mus,lam

c.. arrays read in
      real*8 ptab(0:2000,0:400),ttab(0:2000,0:400),lamtab(0:200,0:239)
      real*8 mutab(0:2000,0:400),stab(0:2000,0:400),msttab(0:2000,0:400)
      real*8 ptab2(0:200,0:200),ttab2(0:200,0:200)
      real*8 mutab2(0:200,0:200),stab2(0:200,0:200),msttab2(0:200,0:200)
      real*8 mustab(0:2000,0:400),mustab2(0:200,0:200)
      real*8 cstab(0:2000,0:400),cstab2(0:200,0:200)
      real*8 ptab3(0:200,0:200),ttab3(0:200,0:200)
      real*8 mutab3(0:200,0:200),stab3(0:200,0:200),msttab3(0:200,0:200)
      real*8 mustab3(0:200,0:200),cstab3(0:200,0:200)

      real*8 latttab(0:43440),latemin,latemax
      real*8 latttab_ht(0:13549),latemin_ht,latemax_ht

c.. general
      integer in, ie, j
     
c. Common-Blocks.
     
      common /tab/ ptab,ttab,mutab,stab,lamtab,ptab2,ttab2,mutab2,stab2,
     &     mustab,mustab2,cstab,cstab2,cstab3,ptab3,ttab3,mutab3
     &     ,stab3,mustab3

      common /lattab/ latttab,latemin,latemax
      common /lattab_ht/ latttab_ht,latemin_ht,latemax_ht
      
c-----|----------------------------------------------------------------|X
      if(eos.eq.2.OR.eos.eq.4.or.eos.eq.6) then
     
      write(0,*)'Read EoS tables for HG...'

c. Open the table files

      open(unit=53,
     $     file='eostables/hadgas_eos.dat')
      open(unit=54,
     $     file='eostables/hg_eos_small.dat')     
      open(unit=55,
     $     file='eostables/hg_eos_mini.dat') 


c. Read in EoS.
c.. first read hadgas_eos.dat
      do 1169 in = 0,400,1
         j = 53
         do 1116 ie = 0,2000,1
            read(j,7787) t,mu,e,p,n,s,mus,lam
            ptab(ie,in) = p
            ttab(ie,in) = t
            mutab(ie,in) = mu
            stab(ie,in) = s
            mustab(ie,in) = mus
            cstab(ie,in) = lam

 1116    continue
 1169 continue
      close(53)

c.. secondly read hg_eos_small.dat
      do 1269 in = 0,200,1
         j = 54
         do 1216 ie = 0,200,1
            read(j,7787) t,mu,e,p,n,s,mus,lam
            ptab2(ie,in) = p
            ttab2(ie,in) = t
            mutab2(ie,in) = mu
            stab2(ie,in) = s
            mustab2(ie,in) = mus
            cstab2(ie,in) = lam

 1216    continue
 1269 continue
      close(54)      

c.. finally read hg_eos_mini.dat
      do 1569 in = 0,200,1
         j = 55
         do 1516 ie = 0,200,1
            read(j,7787) t,mu,e,p,n,s,mus,lam
            ptab3(ie,in) = p
            ttab3(ie,in) = t
            mutab3(ie,in) = mu
            stab3(ie,in) = s
            mustab3(ie,in) = mus
            cstab3(ie,in) = lam
 1516    continue
 1569 continue
      close(55)
 
      else if(eos.eq.3.OR.eos.eq.5) then

      write(0,*)'Read chiral EoS tables...' 

c     Open files of EoS table.
c     
      open(unit=51,
     $     file='eostables/chiraleos.dat')
      open(unit=52,
     $     file='eostables/chiralsmall.dat')     
      open(unit=56,
     $     file='eostables/chiralmini.dat')  



c  Read in EoS.
c
         do 1109 in = 0,400,1
            j = 51
            do 1110 ie = 0,2000,1
            read(j,7777) t,mu,e,p,n,s,mstar,lam
            ptab(ie,in) = p
            ttab(ie,in) = t
            mutab(ie,in) = mu
            stab(ie,in) = s
            msttab(ie,in) = mstar
            cstab(ie,in) = lam
 1110    continue
 1109 continue
      close(51)

      do 1209 in = 0,200,1
         j = 52
         do 1210 ie = 0,200,1
            read(j,7777) t,mu,e,p,n,s,mstar,lam
            ptab2(ie,in) = p
            ttab2(ie,in) = t
            mutab2(ie,in) = mu
            stab2(ie,in) = s
            msttab2(ie,in) = mstar
            cstab2(ie,in) = lam
 1210    continue
 1209 continue
      close(52)

      do 1509 in = 0,200,1
         j = 56
         do 1510 ie = 0,200,1
            read(j,7777) t,mu,e,p,n,s,mstar,lam
            ptab3(ie,in) = p
            ttab3(ie,in) = t
            mutab3(ie,in) = mu
            stab3(ie,in) = s
            msttab3(ie,in) = mstar
            cstab3(ie,in) = lam
 1510    continue
 1509 continue
      close(56)
      
        endif

c..Read in addtional LATTICE EOS for eos=6........

      if(eos.eq.6) then 

       open(unit=59,file='eostables/table_lattice_eos_ext.dat')  

        write(0,*)'Read additional LATTICE EoS tables...'

        do 1601 in = 0,33611,1
          read (59,7797) e,t
c          write(0,*)e,t
          if(in.eq.0) latemin=e
          if(in.eq.33611) latemax=e
          latttab(in) = t
 1601   continue
       close(59)

       open(unit=60,file='eostables/table_lattice_eos_HT.dat')  

        do 1602 in = 0,13549,1
          read (60,7797) e,t
c          write(0,*)e,t
          if(in.eq.0) latemin_ht=e
          if(in.eq.13549) latemax_ht=e
          latttab_ht(in) = t
 1602   continue
       close(60)

      endif

c.................................................

      return

 7777 format(2(1x,f8.3),6(1x,e15.7))
 7787 format(2(1x,f8.3),6(1x,e15.7))
 7797 format(2x,e14.8,2x,e14.8)

      end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE READDIL
c..
c.. This subroutine reads, in following order, these tables:
c..    
c.. 1) the table where the values of 
c..    dR/dM(T,mu,M) (up to factors) are stored.
c.. 2) the table where - for each T and mu -
c..    the values of the maxima of  dR/dM are stored.
c.. 3) the table where the VM in-medium self-energies
c..    are stored
c.. 4) the table where - for each T, mu and M -
c..    the values of the maxima of the function 
c..    ImDrho(M,q;T,mu), with respect to the variable q, 
c..    are stored.
c*****|****************************************************************|X

      SUBROUTINE readdil_elet(vac)
      implicit none
      include 'defs.f'

c-----|----------------------------------------------------------------|X

c.. counter variables
      integer i,j,k
      real*8 massaux,muaux,taux,qaux

c.. Flag for vacuum calculation:
c.. If *vac=1*, vacuum calculation is performed,
c.. if *vac=0*, in-medium calculation is performed.

      integer vac
     
c..   common variables
      real*8 masstab(kmmax)
      real*8 masstabmin,masstabmax
      real*8 ratetab(imumax,jtmax,kmmax)
      real*8 dmd3qmax(imumax,jtmax,kmmax)
      real*8 ratemax(imumax,jtmax)

      common /vmmass/ masstab,masstabmax,masstabmin
      common /drdm/ratetab
      common /maxdr/ratemax
      common /maxdmd3q/dmd3qmax
        
      real*8 mutab(imumax),ttab(jtmax)
      real*8 qtab(kqmax),reself(imumax,jtmax,kqmax),
     &     imself(imumax,jtmax,kqmax)
      real*8 mutabmax,mutabmin,ttabmax,ttabmin
      real*8 qtabmax,qtabmin

      common /tbvarel/mutab,ttab,mutabmax,mutabmin,
     &     ttabmax,ttabmin
      common /eltbself/qtab,reself,imself,qtabmax,qtabmin
      real*8 qup,qdown
      common /qupdow/ qup,qdown

c-----|----------------------------------------------------------------|X

c. 1) READ TABULATED VACUUM/IN-MEDIUM RATES (up to known factors)

      if(dimuon) then
***DI-MUONS***
c.. VACUUM case
      if(vac.eq.1) then
      write(0,*)'Read vacuum dilmuon tables...'
         open(1,file=
     &        'diltables/elka/tab_muon_vacrate.dat', 
     &        status='old')

c.. IN-MEDIUM case     
      else
      write(0,*)'Read in-medium dimuon tables...'
         open(1,file='diltables/elka/tab_muon_imrate.all.dat', 
     &        status='old')
      endif

      else ! di-electron
 
***DI-ELECTRONS***
c.. VACUUM case
      if(vac.eq.1) then
      write(0,*)'Read vacuum dielectron tables...'
         open(1,file=
     &        'diltables/elka/tab_vacrate.dat', 
     &        status='old')

c.. IN-MEDIUM case     
      else
      write(0,*)'Read in-medium dielectron tables...'
         open(1,file='diltables/elka/tab_imrate.all.dat', 
     &        status='old')
      endif
      endif
c.. Here only rates are stored. Values of mu and T will be read later

c.. Loop over baryon chemical potential
      do i=1,imumax
         read(1,*) !skip counter
c.. Loop over temperatures
         do j=1,jtmax
c.. Loop over masses
            do k=1,kmmax
               read(1,12)taux,massaux,ratetab(i,j,k)
               if(i.eq.1.and.j.eq.1) masstab(k)=massaux
            enddo
         enddo
      enddo

      close(1)

      masstabmin=masstab(1)
      masstabmax=masstab(kmmax)
      mmin=masstab(1)
      mmax=masstabmax

c      write(0,*)'masstabmin,masstabmax,mmin,mmax',masstabmin,
c     &            masstabmax,mmin,mmax !Debug only

c-----|----------------------------------------------------------------|X

c. 2)  READ CORRESPONDING MAXIMA (up to known factors)

      if(dimuon) then

c ***DI-MUONS***
       write(0,*)'Read corresponding dimuon maxima...'
c.. VACUUM case
      if(vac.eq.1) then
          open(1,file=
     &        'diltables/elka/max_muon_vacrate.dat', 
     &        status='old')
c.. IN-MEDIUM case   
      else
         open(1,file='diltables/elka/max_muon_imrate.dat', 
     &        status='old')
      endif

      else ! di-electrons

c ***DI-ELECTRONS***
       write(0,*)'Read corresponding dielectron maxima...'

c.. VACUUM case
      if(vac.eq.1) then
          open(1,file=
     &        'diltables/elka/max_vacrate.dat', 
     &        status='old')
c.. IN-MEDIUM case   
      else
         open(1,file='diltables/elka/max_imrate.dat', 
     &        status='old')
      endif
      endif
c.. Loop over baryon chemical potential
      do i=1,imumax
         read(1,*) !skip counter
c.. Loop over temperatures
         do j=1,jtmax
           read(1,13) taux,ratemax(i,j)
c           write(0,*),i,j,taux,ratemax(i,j) ! Debug only
	 enddo
      enddo
      close(1)

c-----|----------------------------------------------------------------|X
c     NOTE:
c     the following two tables are needed only for the in-medium
c     calculation. However, we read them also in the case of the 
c     vacuum calculation, so that we do not have to alter the 
c     structures of the common blocks.
c-----|----------------------------------------------------------------|X

c. 3)  READ TABULATED IN-MEDIUM SELF ENERGIES (vacuum excluded)
       write(0,*)'Read self energies...'

      open(1,file='diltables/elka/tab_self.all.dat',
     &     status='old')

c..   loop over baryon chemical potential
      do i=1,imumax
         read(1,*) !skip counter
c..   loop over temperatures
         do j=1,jtmax
c..   loop over momenta
            do k=1,kqmax
               read(1,14)muaux,taux,qaux,reself(i,j,k),imself(i,j,k)
               if(i.eq.1.and.j.eq.1) qtab(k)=qaux
               if(k.eq.1) ttab(j)=taux
               if(j.eq.1.and.k.eq.1) mutab(i)=muaux
            enddo 
         enddo
      enddo

      close(1)
      
      mutabmax=mutab(imumax)
      mutabmin=mutab(1)
      ttabmax=ttab(jtmax)
      ttabmin=ttab(1)
      qtabmax=qtab(kqmax)
      qtabmin=qtab(1)
      qup=qtabmax
      qdown=qtabmin
 
      write(0,*)'mumax,mumin,tmax,tmin,qmax,qmin,qup,qdown',mutabmax,
     &           mutabmin,ttabmax,ttabmin,qtabmax,qtabmin,qup,qdown ! Debug only

c-----|----------------------------------------------------------------|X

c. 4) READ TABULATED MAXIMA OF THE FUNCTION ImDrho(M,q;T,mu)
      write(0,*)'Read maxima of ImDrho(M,q;T,mu)...'

      if(dimuon) then
c ***DI-MUONS***
      open(1,file='diltables/elka/max_imrho_muon.dat', 
     &	status='old')
 
      else !dielectron
c ***DI-ELECTRONS***     
      open(1,file='diltables/elka/max_imrho.dat', 
     &	status='old')
      endif

c..   here store only values of the maxima. 
c     Values of mu, T and M have been aready stored.

c..   loop over baryon chemical potential
      do i=1,imumax
         read(1,*) !skip counter
c..   loop over temperatures
         do j=1,jtmax
c..   loop over masses
            do k=1,kmmax
               read(1,12)taux,massaux,dmd3qmax(i,j,k)
            enddo
         enddo
      enddo

      close(1)

 12   format(2f12.5,1x,e14.5) 
 13   format(f12.5,1x,e14.5)
 14   format(3f12.5,1x,2e14.5)

      return
      end

c*****|****************************************************************|X
c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE READDIL_RAPP
c..
c.. This subroutine reads, in following order, these tables:
c..    
c.. 1) the table where the values of 
c..    dR/dM are stored.
c.. 2) the table where - for each T and ch. potentials -
c..    the values of the maxima of  dR/dM are stored.
c.. 3) the table where the VM in-medium propagators
c..    are stored
c.. 4) the table where - for each T, ch. potentials and M -
c..    the values of the maxima of the function 
c..    ImDrho(M,q;T,mu), with respect to the variable q, 
c..    are stored.
c*****|****************************************************************|X

      SUBROUTINE readdil_rapp()
      implicit none
      include 'defs.f'

c-----|----------------------------------------------------------------|X

c.. counter variables
      integer i,j,k,l,m,n
      real*8 massaux,rhoaux,taux,piaux,kaux,qaux

c.. Flag for vacuum calculation:
c.. If *vac=1*, vacuum calculation is performed,
c.. if *vac=0*, in-medium calculation is performed.

      integer vac
     
c..   common variables
      real*8 masstab(mmmaxrhom)
      real*8 masstabmin,masstabmax
      real*8 rtab_rho(itmax,jrhomax,kpimax,lkmax,mmmaxrhom)
      real*8 rtab_ome(itmax,jrhomax,kpimax,lkmax,mmmaxrhom)
      real*8 rtab_phi(itmax,jrhomax,kpimax,lkmax,mmmaxrhom)
      real*8 dmd3qmx_rh(itmax,jrhomax,kpimax,lkmax,mmmaxrhom)
      real*8 dmd3qmx_om(itmax,jrhomax,kpimax,lkmax,mmmaxrhom)
      real*8 dmd3qmx_ph(itmax,jrhomax,kpimax,lkmax,mmmaxrhom)
      real*8 rmax_rho(itmax,jrhomax,kpimax,lkmaxrhom)
      real*8 rmax_ome(itmax,jrhomax,kpimax,lkmaxrhom)
      real*8 rmax_phi(itmax,jrhomax,kpimax,lkmaxrhom)

      common /vmmass/ masstab,masstabmax,masstabmin
      common /drdm/rtab_rho,rtab_ome,rtab_phi
      common /maxdr/rmax_rho,rmax_ome,rmax_phi
      common /maxdmd3q/dmd3qmx_rh,dmd3qmx_om,dmd3qmx_ph
        
      real*8 mutab(jrhomax),ttab(itmax),pitab(kpimax),ktab(lkmax)
      real*8 qtab(nqmax)

      real*8 improp_rho(itmax,jrhomax,kpimax,lkmax,mmmaxrhom,nqmax)
      real*8 improp_ome(itmax,jrhomax,kpimax,lkmax,mmmaxrhom,nqmax)
      real*8 improp_phi(itmax,jrhomax,kpimax,lkmax,mmmaxrhom,nqmax)

      real*8 mutabmax,mutabmin,ttabmax,ttabmin,pitabmin,pitabmax
      real*8 qtabmax,qtabmin,ktabmin,ktabmax

      common /tbvar/mutab,ttab,pitab,ktab,mutabmax,mutabmin,
     &     ttabmax,ttabmin,pitabmin,pitabmax,ktabmin,ktabmax
      common /tbself/qtab,improp_rho,improp_ome,improp_phi,
     &     qtabmax,qtabmin
      real*8 qup,qdown
      common /qupdow/ qup,qdown

c-----|----------------------------------------------------------------|X
c.. SET VALUE OF MASS_RESOLUTION

      dmr=dmr_lr

c-----|----------------------------------------------------------------|X

c. 1) READ TABULATED VACUUM/IN-MEDIUM RATES (up to known factors)

      vac=0 ! No vacuum calculation available for this Routine!!!

      if(dimuon) then
***DI-MUONS***
c.. VACUUM case
      if(vac.eq.1) then
       stop'No vacuum calculation possible'
c.. IN-MEDIUM case     
      else
      write(0,*)'Read Rapp spectral function dimuon tables...'
c       open(1,file='diltables/rawaold/Rates-tot.dat', 
c     &      status='old')
       open(1,file='diltables/rawaold/Rates-ext.dat', 
     &      status='old')
      endif

      else ! di-electron
 
***DI-ELECTRONS***
c.. VACUUM case
      if(vac.eq.1) then
       stop'No vacuum calculation possible'
c.. IN-MEDIUM case     
      else
      write(0,*)'Read Rapp spectral function dielectron tables...'
       open(1,file='diltables/rawaold/Rates-el.dat', 
     &      status='old')
      endif
      endif
c.. Here only rates are stored. Values of mu and T will be read later

c.. Loop over temperature
      do i=1,itmax!32
         read(1,*) !skip counter
c.. Loop over densities
         do j=1,jrhomax!9
          do k=1,kpimax!4
           do l=1,lkmax!4
c.. Loop over masses
            do m=1,mmmax!31
               read(1,12)rhoaux,piaux,kaux,massaux,rtab_rho(i,j,k,l,m),
     &                   rtab_ome(i,j,k,l,m),rtab_phi(i,j,k,l,m)
c               write(0,*)i,j,k,l,m,massaux ! Debug only
               if(i.eq.1.AND.j.eq.1.AND.k.eq.1.AND.l.eq.1) then
                masstab(m)=massaux
               endif
            enddo
           enddo
          enddo
         enddo
      enddo

      close(1)

      masstabmin=masstab(1)
      masstabmax=masstab(mmmax)
      mmin=masstab(1)
      mmax=masstabmax

      write(0,*)'masstabmin,masstabmax,mmin,mmax'
      write(0,21)masstabmin,masstabmax,mmin,mmax !Debug only

c-----|----------------------------------------------------------------|X

c. 2)  READ CORRESPONDING MAXIMA (up to known factors)

      if(dimuon) then

c ***DI-MUONS***
       write(0,*)'Read corresponding dimuon maxima...'
c        open(1,file='diltables/rawaold/max-Rates.dat', 
c     &       status='old')
       open(1,file='diltables/rawaold/max_ext-Rates.dat', 
     &       status='old') 
      else ! di-electrons

c ***DI-ELECTRONS***
       write(0,*)'Read corresponding dielectron maxima...'
       open(1,file='diltables/rawaold/max_el-Rates.dat', 
     &        status='old')
       endif
c.. Loop over baryon chemical potential
      do i=1,itmax!32
         read(1,*) !skip counter
c.. Loop over temperatures
       do j=1,jrhomax!9
        do k=1,kpimax!4
         do l=1,lkmax!4
          read(1,13) rhoaux,piaux,kaux,rmax_rho(i,j,k,l),
     &               rmax_ome(i,j,k,l),rmax_phi(i,j,k,l)
c           write(0,*),i,j,k,l,rhoaux,piaux,kaux,rmax_phi(i,j,k,l) ! Debug only
	 enddo
        enddo
       enddo
      enddo

      close(1)

c-----|----------------------------------------------------------------|X
c. 3)  READ TABULATED IN-MEDIUM PROPAGATORS 
       write(0,*)'Read self energies...'

c      open(1,file='diltables/rawaold/ImDrho-total.dat',
c     &     status='old')
      open(1,file='diltables/rawaold/ImDrho-ext.dat',
     &     status='old')

c..   loop over temperature
      do i=1,itmax
         read(1,*) !skip counter
         ttab(i)=0.050d0+((i-1)*0.010d0)
c..   loop over densities
         do j=1,jrhomax!9
          do k=1,kpimax!4
           do l=1,lkmax!4
c..   loop over mass and momenta
            do n=1,nqmax!21
             do m=1,mmmax!31
              read(1,14)rhoaux,piaux,kaux,massaux,qaux,
     &                improp_rho(i,j,k,l,m,n),
     &                improp_ome(i,j,k,l,m,n),improp_phi(i,j,k,l,m,n)
              if(i.eq.1.AND.j.eq.1.AND.k.eq.1.AND.l.eq.1.AND.m.eq.1)then
               qtab(n)=qaux*1d-3
              endif
              if(i.eq.1.AND.k.eq.1.AND.l.eq.1.AND.m.eq.1.AND.n.eq.1)then
               mutab(j)=rhoaux
              endif
              if(i.eq.1.AND.j.eq.1.AND.l.eq.1.AND.m.eq.1.AND.n.eq.1)then
               pitab(k)=piaux*1d-3
              endif
              if(i.eq.1.AND.j.eq.1.AND.k.eq.1.AND.m.eq.1.AND.n.eq.1)then
              ktab(l)=kaux*1d-3
              endif
            enddo
           enddo
          enddo 
         enddo
        enddo
      enddo

      close(1)
      
      mutabmax=mutab(jrhomax)
      mutabmin=mutab(1)
      ttabmax=ttab(itmax)
      ttabmin=ttab(1)
      pitabmax=pitab(kpimax)
      pitabmin=pitab(1)
      ktabmax=ktab(lkmax)
      ktabmin=ktab(1)
      qtabmax=qtab(nqmax)
      qtabmin=qtab(1)
      qup=qtabmax
      qdown=qtabmin
      
      write(0,*)itmax
      write(0,*)'mumax,mumin,tmax,tmin'
      write(0,21)mutabmax,mutabmin,ttabmax,ttabmin ! Debug only
      write(0,*)'pimax,pimin,kmax,kmin,qup,qdown'
      write(0,22)pitabmax,pitabmin,ktabmax,ktabmin,qup,qdown ! Debug only

      rhoeffmax=mutabmax

c-----|----------------------------------------------------------------|X

c. 4) READ TABULATED MAXIMA OF THE FUNCTION ImDrho(M,q;T,mu)
      write(0,*)'Read maxima of ImDrho(M,q;T,mu)...'

c      open(1,file='diltables/rawaold/max_ImDrho.dat', 
c     &	status='old')
      open(1,file='diltables/rawaold/max_ext_ImDrho.dat', 
     &	status='old')
 
c..   here store only values of the maxima. 
c     Values of mu, T and M have been aready stored.

c..   loop over temperature
      do i=1,itmax!32
         read(1,*) !skip counter
c..   loop over densitites
         do j=1,jrhomax!9
          do k=1,kpimax!4
           do l=1,lkmax!4
c..   loop over masses
            do m=1,mmmax
             read(1,12)rhoaux,piaux,kaux,massaux,dmd3qmx_rh(i,j,k,l,m),
     &                 dmd3qmx_om(i,j,k,l,m),dmd3qmx_ph(i,j,k,l,m)
c             write(0,*)'dmd3qmax:rho,omega,phi',dmd3qmx_rh(i,j,k,l,m),
c     &                  dmd3qmx_om(i,j,k,l,m),dmd3qmx_ph(i,j,k,l,m) ! Debug only
            enddo
           enddo
          enddo
         enddo 
      enddo

      close(1)

 12   format(1X,9(E10.4,2X))
 13   format(1X,9(E10.4,2X))
 14   format(1X,9(E10.4,2X))

 21   format(1X,4(f8.4,2X))
 22   format(1X,4(f8.4,2X))

      return
      end

c*****|****************************************************************|X
c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE READDIL_RAPP_RHO
c..
c.. This subroutine reads, in following order, these tables:
c..    
c.. 1) the table where the values of 
c..    dR/dM are stored.
c.. 2) the table where - for each T and ch. potentials -
c..    the values of the maxima of  dR/dM are stored.
c.. 3) the table where the VM in-medium propagators
c..    are stored
c.. 4) the table where - for each T, ch. potentials and M -
c..    the values of the maxima of the function 
c..    ImDrho(M,q;T,mu), with respect to the variable q, 
c..    are stored.
c*****|****************************************************************|X

      SUBROUTINE readdil_rapp_rho()
      implicit none
      include 'defs.f'

c-----|----------------------------------------------------------------|X

c.. counter variables
      integer i,j,k,l,m,n
      real*8 massaux,rhoaux,taux,piaux,kaux,qaux

c.. Flag for vacuum calculation:
c.. If *vac=1*, vacuum calculation is performed,
c.. if *vac=0*, in-medium calculation is performed.

      integer vac
     
c..   common variables
      real*8 masstab(mmmaxrhom)
      real*8 masstabmin,masstabmax
      real*8 rtab_rho(itmax,jrhomax,kpimax,lkmax,mmmaxrhom)
      real*8 rtab_ome(itmax,jrhomax,kpimax,lkmax,mmmaxrhom)
      real*8 rtab_phi(itmax,jrhomax,kpimax,lkmax,mmmaxrhom)
      real*8 dmd3qmx_rh(itmax,jrhomax,kpimax,lkmax,mmmaxrhom)
      real*8 dmd3qmx_om(itmax,jrhomax,kpimax,lkmax,mmmaxrhom)
      real*8 dmd3qmx_ph(itmax,jrhomax,kpimax,lkmax,mmmaxrhom)
      real*8 rmax_rho(itmax,jrhomax,kpimax,lkmax)
      real*8 rmax_ome(itmax,jrhomax,kpimax,lkmax)
      real*8 rmax_phi(itmax,jrhomax,kpimax,lkmax)

      common /vmmass/ masstab,masstabmax,masstabmin
      common /drdm/rtab_rho,rtab_ome,rtab_phi
      common /maxdr/rmax_rho,rmax_ome,rmax_phi
      common /maxdmd3q/dmd3qmx_rh,dmd3qmx_om,dmd3qmx_ph
        
      real*8 mutab(jrhomax),ttab(itmax),pitab(kpimax),ktab(lkmax)
      real*8 qtab(nqmax)

      real*8 improp_rho(itmax,jrhomax,kpimax,lkmax,mmmaxrhom,nqmax)
      real*8 improp_ome(itmax,jrhomax,kpimax,lkmax,mmmaxrhom,nqmax)
      real*8 improp_phi(itmax,jrhomax,kpimax,lkmax,mmmaxrhom,nqmax)

      real*8 mutabmax,mutabmin,ttabmax,ttabmin,pitabmin,pitabmax
      real*8 qtabmax,qtabmin,ktabmin,ktabmax

      common /tbvar/mutab,ttab,pitab,ktab,mutabmax,mutabmin,
     &     ttabmax,ttabmin,pitabmin,pitabmax,ktabmin,ktabmax
      common /tbself/qtab,improp_rho,improp_ome,improp_phi,
     &     qtabmax,qtabmin
      real*8 qup,qdown
      common /qupdow/ qup,qdown

c-----|----------------------------------------------------------------|X
c.. SET VALUE OF MASS_RESOLUTION

      dmr=dmr_hr

c-----|----------------------------------------------------------------|X

c. 1) READ TABULATED VACUUM/IN-MEDIUM RATES (up to known factors)

      vac=0 ! No vacuum calculation available for this Routine!!!

      if(dimuon) then
***DI-MUONS***
c.. VACUUM case
      if(vac.eq.1) then
       stop'No vacuum calculation possible'
c.. IN-MEDIUM case     
      else
      write(0,*)'Read Rapp spectral function dimuon tables...'
       open(1,file='diltables/rawa/Rates-or-mu.dat', 
     &      status='old')
      endif
      else ! di-electron
 
***DI-ELECTRONS***
c.. VACUUM case
      if(vac.eq.1) then
       stop'No vacuum calculation possible'
c.. IN-MEDIUM case     
      else
      write(0,*)'Read Rapp spectral function dielectron tables...'
       open(1,file='diltables/rawa/Rates-or-f.dat', 
     &      status='old')
      endif
      endif
c.. Here only rates are stored. Values of mu and T will be read later

c.. Loop over temperature
      do i=1,itmaxor!32
         read(1,*) !skip counter
c.. Loop over densities
         do j=1,jrhomaxor!13
          do k=1,kpimaxor!4
           do l=1,lkmaxor!1
c.. Loop over masses
            do m=1,mmmaxor!76
               read(1,12)rhoaux,piaux,massaux,rtab_rho(i,j,k,l,m),
     &                   rtab_ome(i,j,k,l,m)
c               write(0,*)i,j,k,l,m,massaux ! Debug only
               if(i.eq.1.AND.j.eq.1.AND.k.eq.1.AND.l.eq.1) then
                masstab(m)=massaux
               endif
            enddo
           enddo
          enddo
         enddo
      enddo

      close(1)

      masstabmin=masstab(1)
      masstabmax=masstab(mmmaxor)
      mmin=masstab(1)
      mmax=masstabmax

      write(0,*)'masstabmin,masstabmax,mmin,mmax'
      write(0,21)masstabmin,masstabmax,mmin,mmax !Debug only

c-----|----------------------------------------------------------------|X

c. 2)  READ CORRESPONDING MAXIMA (up to known factors)

      if(dimuon) then

c ***DI-MUONS***
       write(0,*)'Read corresponding dimuon maxima...'
       open(1,file='diltables/rawa/max_Rates-or-mu.dat', 
     &       status='old') 
      else ! di-electrons

c ***DI-ELECTRONS***
       write(0,*)'Read corresponding dielectron maxima...'
       open(1,file='diltables/rawa/max_Rates-or-f.dat', 
     &        status='old')
       endif
c.. Loop over baryon chemical potential
      do i=1,itmaxor!32
         read(1,*) !skip counter
c.. Loop over temperatures
       do j=1,jrhomaxor!13
        do k=1,kpimaxor!4
         do l=1,lkmaxor!1
          read(1,13) rhoaux,piaux,rmax_rho(i,j,k,l),rmax_ome(i,j,k,l)
c           write(0,*),i,j,k,l,rhoaux,piaux,kaux,rmax_rho(i,j,k,l),
c     &               rmax_ome(i,j,k,l) ! Debug only
	 enddo
        enddo
       enddo
      enddo

      close(1)

c-----|----------------------------------------------------------------|X
c. 3) READ TABULATED IN-MEDIUM PROPAGATORS 
      write(0,*)'Read self energies...'

      open(1,file='diltables/rawa/ImDrho-or-f.dat',
     &     status='old')

c..   loop over temperature
      do i=1,itmaxor
         read(1,*) !skip counter
         ttab(i)=0.050d0+((i-1)*0.010d0)
c..   loop over densities
         do j=1,jrhomaxor!13
          do k=1,kpimaxor!4
           do l=1,lkmaxor!1
c..   loop over mass and momenta
            do n=1,nqmaxor!21
             do m=1,mmmaxor!76
              read(1,14)rhoaux,piaux,massaux,qaux,improp_rho(i,j,k,l,m,n),
     &                  improp_ome(i,j,k,l,m,n)
              if(i.eq.1.AND.j.eq.1.AND.k.eq.1.AND.l.eq.1.AND.m.eq.1)then
               qtab(n)=qaux*1d-3
              endif
              if(i.eq.1.AND.k.eq.1.AND.l.eq.1.AND.m.eq.1.AND.n.eq.1)then
               mutab(j)=rhoaux
              endif
              if(i.eq.1.AND.j.eq.1.AND.l.eq.1.AND.m.eq.1.AND.n.eq.1)then
               pitab(k)=piaux*1d-3
              endif
              if(i.eq.1.AND.j.eq.1.AND.k.eq.1.AND.m.eq.1.AND.n.eq.1)then
              ktab(l)=0.0d0
              endif
            enddo
           enddo
          enddo 
         enddo
        enddo
      enddo

      close(1)
      
      mutabmax=mutab(jrhomaxor)
      mutabmin=mutab(1)
      ttabmax=ttab(itmaxor)
      ttabmin=ttab(1)
      pitabmax=pitab(kpimaxor)
      pitabmin=pitab(1)
      ktabmax=ktab(lkmaxor)
      ktabmin=ktab(1)
      qtabmax=qtab(nqmaxor)
      qtabmin=qtab(1)
      qup=qtabmax
      qdown=qtabmin
 
      write(0,*)'mumax,mumin,tmax,tmin'
      write(0,21)mutabmax,mutabmin,ttabmax,ttabmin ! Debug only
      write(0,*)'pimax,pimin,kmax,kmin,qup,qdown'
      write(0,22)pitabmax,pitabmin,ktabmax,ktabmin,qup,qdown ! Debug only

      rhoeffmax=mutabmax

c-----|----------------------------------------------------------------|X

c. 4) READ TABULATED MAXIMA OF THE FUNCTION ImDrho(M,q;T,mu)
      write(0,*)'Read maxima of ImDrho(M,q;T,mu)...'

      open(1,file='diltables/rawa/max_ImDrho-or-f.dat', 
     &	status='old')

c..   here store only values of the maxima. 
c     Values of mu, T and M have been aready stored.

c..   loop over temperature
      do i=1,itmaxor!32
         read(1,*) !skip counter
c..   loop over densitites
         do j=1,jrhomaxor!13
          do k=1,kpimaxor!4
           do l=1,lkmaxor!1
c..   loop over masses
            do m=1,mmmaxor!76
             read(1,12)rhoaux,piaux,massaux,dmd3qmx_rh(i,j,k,l,m),
     &                 dmd3qmx_om(i,j,k,l,m)
c             write(0,*)i,j,k,l,m
c             write(0,*)'dmd3qmax:rho,omega,phi',dmd3qmx_rh(i,j,k,l,m),
c     &                  dmd3qmx_om(i,j,k,l,m),dmd3qmx_ph(i,j,k,l,m) ! Debug only
            enddo
           enddo
          enddo
         enddo 
      enddo
      close(1)

 12   format(1X,9(E10.4,2X))
 13   format(1X,9(E10.4,2X))
 14   format(1X,9(E10.4,2X))

 21   format(1X,4(f8.4,2X))
 22   format(1X,4(f8.4,2X))

      return
      end

c*****|****************************************************************|X
c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE READDIL_RAPP_PHI
c..
c.. This subroutine reads, in following order, these tables:
c..    
c.. 1) the table where the values of 
c..    dR/dM are stored.
c.. 2) the table where - for each T and ch. potentials -
c..    the values of the maxima of  dR/dM are stored.
c.. 3) the table where the VM in-medium propagators
c..    are stored
c.. 4) the table where - for each T, ch. potentials and M -
c..    the values of the maxima of the function 
c..    ImDrho(M,q;T,mu), with respect to the variable q, 
c..    are stored.
c*****|****************************************************************|X

      SUBROUTINE readdil_rapp_phi()
      implicit none
      include 'defs.f'

c-----|----------------------------------------------------------------|X

c.. counter variables
      integer i,j,k,l,m,n
      real*8 massaux,rhoaux,taux,piaux,kaux,qaux

c.. Flag for vacuum calculation:
c.. If *vac=1*, vacuum calculation is performed,
c.. if *vac=0*, in-medium calculation is performed.

      integer vac
     
c..   common variables
      real*8 masstab(mmmaxrhom)
      real*8 masstabmin,masstabmax
      real*8 rtab_rho(itmax,jrhomax,kpimax,lkmax,mmmaxrhom)
      real*8 rtab_ome(itmax,jrhomax,kpimax,lkmax,mmmaxrhom)
      real*8 rtab_phi(itmax,jrhomax,kpimax,lkmax,mmmaxrhom)
      real*8 dmd3qmx_rh(itmax,jrhomax,kpimax,lkmax,mmmaxrhom)
      real*8 dmd3qmx_om(itmax,jrhomax,kpimax,lkmax,mmmaxrhom)
      real*8 dmd3qmx_ph(itmax,jrhomax,kpimax,lkmax,mmmaxrhom)
      real*8 rmax_rho(itmax,jrhomax,kpimax,lkmax)
      real*8 rmax_ome(itmax,jrhomax,kpimax,lkmax)
      real*8 rmax_phi(itmax,jrhomax,kpimax,lkmax)

      common /vmmass/ masstab,masstabmax,masstabmin
      common /drdm/rtab_rho,rtab_ome,rtab_phi
      common /maxdr/rmax_rho,rmax_ome,rmax_phi
      common /maxdmd3q/dmd3qmx_rh,dmd3qmx_om,dmd3qmx_ph
        
      real*8 mutab(jrhomax),ttab(itmax),pitab(kpimax),ktab(lkmax)
      real*8 qtab(nqmax)

      real*8 improp_rho(itmax,jrhomax,kpimax,lkmax,mmmaxrhom,nqmax)
      real*8 improp_ome(itmax,jrhomax,kpimax,lkmax,mmmaxrhom,nqmax)
      real*8 improp_phi(itmax,jrhomax,kpimax,lkmax,mmmaxrhom,nqmax)

      real*8 mutabmax,mutabmin,ttabmax,ttabmin,pitabmin,pitabmax
      real*8 qtabmax,qtabmin,ktabmin,ktabmax

      common /tbvar/mutab,ttab,pitab,ktab,mutabmax,mutabmin,
     &     ttabmax,ttabmin,pitabmin,pitabmax,ktabmin,ktabmax
      common /tbself/qtab,improp_rho,improp_ome,improp_phi,
     &     qtabmax,qtabmin
      real*8 qup,qdown
      common /qupdow/ qup,qdown

c-----|----------------------------------------------------------------|X
c.. SET VALUE OF MASS_RESOLUTION

      dmr=dmr_hr

c-----|----------------------------------------------------------------|X

c. 1) READ TABULATED VACUUM/IN-MEDIUM RATES (up to known factors)

      vac=0 ! No vacuum calculation available for this Routine!!!

      if(dimuon) then
***DI-MUONS***
c.. VACUUM case
      if(vac.eq.1) then
       stop'No vacuum calculation possible'
c.. IN-MEDIUM case     
      else
      write(0,*)'Read Rapp spectral function dimuon tables...'
c       open(1,file='diltables/rawa/Rates-tot.dat', 
c     &      status='old')
       open(1,file='diltables/rawa/Rates-phi-el.dat', 
     &      status='old')
      endif

      else ! di-electron
 
***DI-ELECTRONS***
c.. VACUUM case
      if(vac.eq.1) then
       stop'No vacuum calculation possible'
c.. IN-MEDIUM case     
      else
      write(0,*)'Read Rapp spectral function dielectron tables...'
       open(1,file='diltables/rawa/Rates-phi-el.dat', 
     &      status='old')
      endif
      endif
c.. Here only rates are stored. Values of mu and T will be read later

c.. Loop over temperature
      do i=1,itmaxphi!32
         read(1,*) !skip counter
c.. Loop over densities
         do j=1,jrhomaxphi!13
          do k=1,kpimaxphi!4
           do l=1,lkmaxphi!4
c.. Loop over masses
            do m=1,mmmaxphi!101
               read(1,12)rhoaux,piaux,kaux,massaux,rtab_phi(i,j,k,l,m)
c               write(0,*)i,j,k,l,m,massaux ! Debug only
               if(i.eq.1.AND.j.eq.1.AND.k.eq.1.AND.l.eq.1) then
                masstab(m)=massaux
               endif
            enddo
           enddo
          enddo
         enddo
      enddo

      close(1)

      masstabmin=masstab(1)
      masstabmax=masstab(mmmaxphi)
      mmin=masstab(1)
      mmax=masstabmax

      write(0,*)'masstabmin,masstabmax,mmin,mmax'
      write(0,21)masstabmin,masstabmax,mmin,mmax !Debug only

c-----|----------------------------------------------------------------|X

c. 2)  READ CORRESPONDING MAXIMA (up to known factors)

      if(dimuon) then

c ***DI-MUONS***
       write(0,*)'Read corresponding dimuon maxima...'
       open(1,file='diltables/rawa/max_Rates-phi-mu.dat', 
     &       status='old') 
      else ! di-electrons

c ***DI-ELECTRONS***
       write(0,*)'Read corresponding dielectron maxima...'
       open(1,file='diltables/rawa/max_Rates-phi-el.dat', 
     &        status='old')
       endif
c.. Loop over baryon chemical potential
      do i=1,itmaxphi!32
         read(1,*) !skip counter
c.. Loop over temperatures
       do j=1,jrhomaxphi!13
        do k=1,kpimaxphi!4
         do l=1,lkmaxphi!4
          read(1,13) rhoaux,piaux,kaux,rmax_phi(i,j,k,l)
c           write(0,*),i,j,k,l,rhoaux,piaux,kaux,rmax_rho(i,j,k,l),
c     &               rmax_ome(i,j,k,l) ! Debug only
	 enddo
        enddo
       enddo
      enddo

      close(1)

c-----|----------------------------------------------------------------|X
c. 3)  READ TABULATED IN-MEDIUM PROPAGATORS 
       write(0,*)'Read self energies...'

      open(1,file='diltables/rawa/ImDrho-phi.dat',
     &     status='old')

c..   loop over temperature
      do i=1,itmaxphi
         read(1,*) !skip counter
         ttab(i)=0.050d0+((i-1)*0.010d0)
c..   loop over densities
         do j=1,jrhomaxphi!13
          do k=1,kpimaxphi!4
           do l=1,lkmaxphi!4
c..   loop over mass and momenta
            do n=1,nqmaxphi!21
             do m=1,mmmaxphi!101
              read(1,14)rhoaux,piaux,kaux,massaux,qaux,improp_phi(i,j,k,l,m,n)
              if(i.eq.1.AND.j.eq.1.AND.k.eq.1.AND.l.eq.1.AND.m.eq.1)then
               qtab(n)=qaux*1d-3
              endif
              if(i.eq.1.AND.k.eq.1.AND.l.eq.1.AND.m.eq.1.AND.n.eq.1)then
               mutab(j)=rhoaux
              endif
              if(i.eq.1.AND.j.eq.1.AND.l.eq.1.AND.m.eq.1.AND.n.eq.1)then
               pitab(k)=piaux*1d-3
              endif
              if(i.eq.1.AND.j.eq.1.AND.k.eq.1.AND.m.eq.1.AND.n.eq.1)then
              ktab(l)=kaux*1d-3
              endif
            enddo
           enddo
          enddo 
         enddo
        enddo
      enddo

      close(1)
      
      mutabmax=mutab(jrhomaxphi)
      mutabmin=mutab(1)
      ttabmax=ttab(itmaxphi)
      ttabmin=ttab(1)
      pitabmax=pitab(kpimaxphi)
      pitabmin=pitab(1)
      ktabmax=ktab(lkmaxphi)
      ktabmin=ktab(1)
      qtabmax=qtab(nqmaxphi)
      qtabmin=qtab(1)
      qup=qtabmax
      qdown=qtabmin
 
      write(0,*)'mumax,mumin,tmax,tmin'
      write(0,21)mutabmax,mutabmin,ttabmax,ttabmin ! Debug only
      write(0,*)'pimax,pimin,kmax,kmin,qup,qdown'
      write(0,22)pitabmax,pitabmin,ktabmax,ktabmin,qup,qdown ! Debug only

      rhoeffmax=mutabmax

c-----|----------------------------------------------------------------|X

c. 4) READ TABULATED MAXIMA OF THE FUNCTION ImDrho(M,q;T,mu)
      write(0,*)'Read maxima of ImDrho(M,q;T,mu)...'


      open(1,file='diltables/rawa/max_ImDrho-phi.dat', 
     &	status='old')
 
c..   here store only values of the maxima. 
c     Values of mu, T and M have been aready stored.

c..   loop over temperature
      do i=1,itmaxphi!32
         read(1,*) !skip counter
c..   loop over densitites
         do j=1,jrhomaxphi!13
          do k=1,kpimaxphi!4
           do l=1,lkmaxphi!4
c..   loop over masses
            do m=1,mmmaxphi!101
             read(1,12)rhoaux,piaux,kaux,massaux,dmd3qmx_ph(i,j,k,l,m)
c             write(0,*)'dmd3qmax:rho,omega,phi',dmd3qmx_rh(i,j,k,l,m),
c     &                  dmd3qmx_om(i,j,k,l,m),dmd3qmx_ph(i,j,k,l,m) ! Debug only
            enddo
           enddo
          enddo
         enddo 
      enddo

      close(1)

 12   format(1X,9(E10.4,2X))
 13   format(1X,9(E10.4,2X))
 14   format(1X,9(E10.4,2X))

 21   format(1X,4(f8.4,2X))
 22   format(1X,4(f8.4,2X))

      return
      end

c*****|****************************************************************|X
