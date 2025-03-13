c234567**1*********2*********3*********4*********5*********6*********7**X
c.. FUNCTION TEMP
c..
c.. This function determines the temperature of the underlying EoS.
c*****|****************************************************************|X

      FUNCTION temp(e,n)
      implicit none

      INCLUDE 'defs.f'

c-----|----------------------------------------------------------------|X

c.. arrays from tables
      real*8 ptab(0:2000,0:400),ttab(0:2000,0:400),lamtab(0:200,0:239)
      real*8 mutab(0:2000,0:400),stab(0:2000,0:400)
      real*8 ptab2(0:200,0:200),ttab2(0:200,0:200)
      real*8 mutab2(0:200,0:200),stab2(0:200,0:200),msttab2(0:200,0:200)
      real*8 mustab(0:2000,0:400),mustab2(0:200,0:200)
      real*8 cstab(0:2000,0:400),cstab2(0:200,0:200)
      real*8 ptab3(0:200,0:200),ttab3(0:200,0:200)
      real*8 mutab3(0:200,0:200),stab3(0:200,0:200),msttab3(0:200,0:200)
      real*8 mustab3(0:200,0:200),cstab3(0:200,0:200)

      real*8 latttab(0:43440),latemin,latemax
      real*8 latttab_ht(0:13549),latemin_ht,latemax_ht

c.. Common-block variables
      integer stabil,antit

c. function variables     
      real*8 e,n,temp,et0,mu
      real*8 de,dn,p1,p2,p3,p4,p13,p24
      real*8 p12,p34

      real*8 lattemp,mix
     
c.     Common-blocks.
      common /tab/ ptab,ttab,mutab,stab,lamtab,ptab2,ttab2,mutab2,stab2,
     &     mustab,mustab2,cstab,cstab2,cstab3,ptab3,ttab3,mutab3
     &     ,stab3,mustab3

      common /lattab/ latttab,latemin,latemax
      common /lattab_ht/ latttab_ht,latemin_ht,latemax_ht

c-----|----------------------------------------------------------------|X

      temp=0.0d0
      lattemp=0.0d0

c      write(0,*)'in function temp: e,n',e,n ! Debug only

c. Case of negative baryon density      
      stabil= 0
      antit = 0
      if (n.lt.0)then
         n=-n
         antit = 1
      end if

c. Now CALCULATE TEMPERATURE (in MeV).


c. *eos=1*: Temperature is that of an ULTRARELATIVISTIC GAS
c. of particles with one degree of freedom.
            
      if (eos.eq.1) then
       temp = (e*e_0*0.33d0*hc3*30d0/pi2)**0.25d0
c       write(0,*)'temp for eos=1:',temp ! Debug only

c. *eos=2*: Use tabellized equation of state for a HADRON GAS.
     
c     If e and n are both smaller than e_0, n_0, the nuclei are in the
c     ground state.      

      else if (eos.eq.2.OR.eos.eq.4.OR.eos.eq.6) then      
         if (e.le.1000.0d0) then
c... from hg_eos_mini.dat
            if((e.lt.0.1d0).and.(n.lt.0.02d0)) then
c               write(0,*)'read hg_eos_mini' !Debug only
               de = 0.0005d0
               dn = 0.0001d0
               p1 = ttab3(idint(e/de),idint(n/dn))
               p2 = ttab3(idint(e/de)+1,idint(n/dn))
               p3 = ttab3(idint(e/de),idint(n/dn)+1)
               p4 = ttab3(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               temp = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
c... from hg_eos_small.dat
            if((e.lt.10.0d0).and.(n.lt.2.0d0).and.
     &              ((e.ge.0.1d0).or.(n.ge.0.02d0))) then
c               write(0,*)'read hg_eos_small' !Debug only
               de = 0.05d0
               dn = 0.01d0
               p1 = ttab2(idint(e/de),idint(n/dn))
               p2 = ttab2(idint(e/de)+1,idint(n/dn))
               p3 = ttab2(idint(e/de),idint(n/dn)+1)
               p4 = ttab2(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               temp = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
c... from hadgas_eos.dat
            if (((e.ge.10.0d0).or.(n.ge.2.0d0)).and.
     &          (e.lt.1000.0d0)) then
c               write(0,*)'read hadgas_eos' !Debug only
               de = 0.5d0
               dn = 0.1d0
               if(n.ge.40.0d0) n=39.99d0
               p1 = ttab(idint(e/de),idint(n/dn))
               p2 = ttab(idint(e/de)+1,idint(n/dn))
               p3 = ttab(idint(e/de),idint(n/dn)+1)
               p4 = ttab(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               temp = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if 

         else if (e.gt.1000d0) then
            temp=350.0d0
         else
            temp = 0d0
         end if
            

      else if (eos.eq.3.OR.eos.eq.5) then      
         if (e.le.400.0d0) then

            if((e.lt.0.1d0).and.(n.lt.0.02d0)) then
               de = 0.0005d0
               dn = 0.0001d0
               p1 = ttab3(idint(e/de),idint(n/dn))
               p2 = ttab3(idint(e/de)+1,idint(n/dn))
               p3 = ttab3(idint(e/de),idint(n/dn)+1)
               p4 = ttab3(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               temp = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            if((e.lt.10.0d0).and.(n.lt.2.0d0).and.
     $              ((e.ge.0.1d0).or.(n.ge.0.02d0))) then
               de = 0.05d0
               dn = 0.01d0
               p1 = ttab2(idint(e/de),idint(n/dn))
               p2 = ttab2(idint(e/de)+1,idint(n/dn))
               p3 = ttab2(idint(e/de),idint(n/dn)+1)
               p4 = ttab2(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               temp = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            if ((e.ge.10.0d0).or.(n.ge.2.0d0)) then
               de = 0.5d0
               dn = 0.1d0

               p1 = ttab(idint(e/de),idint(n/dn))
               p2 = ttab(idint(e/de)+1,idint(n/dn))
               p3 = ttab(idint(e/de),idint(n/dn)+1)
               p4 = ttab(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               temp = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if 
            
         else if (e.ge.400.0d0) then
            temp = 350.0d0
c     else
            
c     temp = 0d0
            
         end if
         
      end if

c..TEMPERATURE from LATTICE EOS for eos=6..............................|X

      if(eos.eq.6) then
        
       if(e.gt.latemax_ht) then
        lattemp=1000.0d0
        temp=1000.0d0
        goto 350
       elseif(e.gt.latemax.AND.e.lt.latemax_ht) then

        de = 1.0d0

        p1 = latttab_ht(idint((e-latemin_ht)/de))
        p2 = latttab_ht(idint((e-latemin_ht)/de)+1)
         
        lattemp = (p2-p1)/de*
     &            ((e-latemin_ht)-dfloat(idint((e-latemin_ht)/de))*de) + p1
        temp=lattemp

       elseif(e.gt.latemin.AND.e.lt.latemax) then 

        de = 0.01d0

        p1 = latttab(idint((e-latemin)/de))
        p2 = latttab(idint((e-latemin)/de)+1)
         
        lattemp = (p2-p1)/de*
     &            ((e-latemin)-dfloat(idint((e-latemin)/de))*de) + p1

c        write(0,*)'temp,lattemp',temp,lattemp !Debug only
c...OLD MATCHING!:
c        if(lattemp.gt.170.0d0.AND.lattemp.gt.temp) temp=lattemp
c        if(lattemp.gt.152.0d0.AND.lattemp.lt.170.0d0.AND.lattemp.gt.temp) then
c          mix=((lattemp-152.0d0)/18.0d0)**(1.0d0/7.0d0)
c          temp=mix*lattemp+(1.0d0-mix)*temp
c          write(0,*)'mixtemp =',temp !Debug only
c        endif
c...NEW MATCHING!:
        if(lattemp.gt.180.0d0) temp=lattemp
        if(lattemp.gt.170.0d0.AND.lattemp.le.180.0d0
     &     .AND.lattemp.gt.(1.2d0*temp)) temp=0.2d0*temp+0.8*lattemp
        if(lattemp.gt.170.0d0.AND.lattemp.le.180.0d0
     &     .AND.lattemp.le.(1.2d0*temp)) temp=lattemp
        if(lattemp.gt.160.0d0.AND.lattemp.lt.170.0d0.AND.lattemp.gt.temp) then
          mix=0.8d0*((lattemp-150.0d0)/20.0d0)**(10.0d0)
          temp=mix*lattemp+(1.0d0-mix)*temp
c          write(0,*)'mixtemp =',temp !Debug only
        endif
       else
        lattemp=0.0d0
       endif
      
      endif

 350  continue

c.....|................................................................|X

      if ((e.le.1d0).and.(n.le.1d0).and.(stabil.eq.1).and.(n.ge.-1d0))
     &     then
         temp = 0d0
      end if
      
      
      if (antit.eq.1)then
         n=-n
      end if    
      return
      end     

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. FUNCTION CHEM
c..
c.. This function determines the chemical potential of the underlying 
c.. EoS.
c*****|****************************************************************|X

      FUNCTION chem(e,n)
      implicit none

      INCLUDE 'defs.f'

c-----|----------------------------------------------------------------|X

c.. Read in from tables
      real*8 ptab(0:2000,0:400),ttab(0:2000,0:400),lamtab(0:200,0:239)
      real*8 mutab(0:2000,0:400),stab(0:2000,0:400)
      real*8 ptab2(0:200,0:200),ttab2(0:200,0:200)
      real*8 mutab2(0:200,0:200),stab2(0:200,0:200),msttab2(0:200,0:200)
      real*8 mustab(0:2000,0:400),mustab2(0:200,0:200)
      real*8 cstab(0:2000,0:400),cstab2(0:200,0:200)
      real*8 ptab3(0:200,0:200),ttab3(0:200,0:200)
      real*8 mutab3(0:200,0:200),stab3(0:200,0:200),msttab3(0:200,0:200)
      real*8 mustab3(0:200,0:200),cstab3(0:200,0:200)
c.. Common-block variables
      integer stabil,antic

c.. Function variables
      real*8 e,n,chem,et0,temp
      real*8 de,dn,p1,p2,p3,p4,p13,p24
      real*8 p12,p34

c. Common-blocks
      common /tab/ ptab,ttab,mutab,stab,lamtab,ptab2,ttab2,mutab2,stab2,
     &     mustab,mustab2,cstab,cstab2,cstab3,ptab3,ttab3,mutab3
     &     ,stab3,mustab3

c-----|----------------------------------------------------------------|X

c. Case of negative baryon density
      stabil= 0
      antic = 0
      if (n.lt.0)then
         n=-n
         antic = 1
      end if

c. Now calculate CHEMICAL POTENTIAL (in MeV)

c.. *eos=1*: Chemical potential is *zero* for ULTRARELATIVISTIC GAS.
         
      if (eos.eq.1) then
         chem = 0d0
     
c.. *eos=2*: Use tabellized equation of state for HADRON GAS.  

c... If e and n are both smaller than e_0, n_0, the nuclei are in the
c... ground state.      

      else if (eos.eq.2.OR.eos.eq.4.OR.eos.eq.6) then   
         if (eos.eq.6.AND.e.gt.4.97d0) then
           chem=0.0d0
           return
         endif     
         if (e.lt.1000.0d0) then
c... from hg_eos_mini.dat
            if((e.lt.0.1d0).and.(n.lt.0.02d0)) then
               de = 0.0005d0
               dn = 0.0001d0
               p1 = mutab3(idint(e/de),idint(n/dn))
               p2 = mutab3(idint(e/de)+1,idint(n/dn))
               p3 = mutab3(idint(e/de),idint(n/dn)+1)
               p4 = mutab3(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3

               if(p4.ne.0.0d0.AND.p2.eq.0.0d0) then
                  p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p4
               endif

               if(p3.ne.0.0d0.AND.p1.eq.0.0d0) then
                  p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p3
               endif

               chem = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
 
            end if
c... from hg_eos_small.dat
            if((e.lt.10.0d0).and.(n.lt.2.0d0).and.
     &          ((e.ge.0.1d0).or.(n.ge.0.02d0))) then
               de = 0.05d0
               dn = 0.01d0
               p1 = mutab2(idint(e/de),idint(n/dn))
               p2 = mutab2(idint(e/de)+1,idint(n/dn))
               p3 = mutab2(idint(e/de),idint(n/dn)+1)
               p4 = mutab2(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3

               if(p4.ne.0.0d0.AND.p2.eq.0.0d0) then
                  p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p4
               endif

               if(p3.ne.0.0d0.AND.p1.eq.0.0d0) then
                  p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p3
               endif
                        
               chem = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13

            end if
c... from hadgas_eos.dat
            if (((e.ge.10.0d0).or.(n.ge.2.0d0)).and.(e.lt.1000.0d0)) then       
               de = 0.5d0
               dn = 0.1d0
               if(n.ge.40.0d0) n=39.99d0
               p1 = mutab(idint(e/de),idint(n/dn))
               p2 = mutab(idint(e/de)+1,idint(n/dn))
               p3 = mutab(idint(e/de),idint(n/dn)+1)
               p4 = mutab(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
 
               if(p4.ne.0.0d0.AND.p2.eq.0.0d0) then
                  p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p4
               endif

               if(p3.ne.0.0d0.AND.p1.eq.0.0d0) then
                  p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p3
               endif
           
               chem = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
               
            end if
         else if (e.gt.1000.0d0) then
            chem = 1d0
!     else
!     chem = 0d0
         end if
         
      else if (eos.eq.3.OR.eos.eq.5) then        
         if (e.lt.400.0d0) then

            if((e.lt.0.1d0).and.(n.lt.0.02d0)) then
               de = 0.0005d0
               dn = 0.0001d0
               p1 = mutab3(idint(e/de),idint(n/dn))
               p2 = mutab3(idint(e/de)+1,idint(n/dn))
               p3 = mutab3(idint(e/de),idint(n/dn)+1)
               p4 = mutab3(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
            
               chem = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            if((e.lt.10.0d0).and.(n.lt.2.0d0).and.
     $           ((e.ge.0.1d0).or.(n.ge.0.02d0))) then
               de = 0.05d0
               dn = 0.01d0
               p1 = mutab2(idint(e/de),idint(n/dn))
               p2 = mutab2(idint(e/de)+1,idint(n/dn))
               p3 = mutab2(idint(e/de),idint(n/dn)+1)
               p4 = mutab2(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
            
               chem = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            if ((e.ge.10.0d0).or.(n.ge.2.0d0)) then       
               de = 0.5d0
               dn = 0.1d0

               p1 = mutab(idint(e/de),idint(n/dn))
               p2 = mutab(idint(e/de)+1,idint(n/dn))
               p3 = mutab(idint(e/de),idint(n/dn)+1)
               p4 = mutab(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
            
               chem = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
               
            end if
         else if (e.gt.400.0d0) then
c     JS muss bestimmt werden
            chem = 1d0
c     else
c     chem = 0d0
         end if

      end if
      
      if ((e.le.1d0).and.(n.le.1d0).and.(stabil.eq.1)) then
         chem = e_0/n_0/3
      end if
      
            
c. Case of negative baryon density      
      if (antic.eq.1)then
         n=-n
         chem = -chem
      end if

      return
      end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c..  FUNCTION CS
c..
c..  This function-subprogram determines the fraction of QGP
c*****|****************************************************************|X

      FUNCTION cs(e,n)
      INCLUDE 'defs.f'
c
c  This function-subprogram determines the chemical potential of the
c  underlying EoS.
c  It is used in subroutine fileo.
c
c  Type declarations for variables in common-blocks.
c
  
      real*8 ptab(0:2000,0:400),ttab(0:2000,0:400),lamtab(0:200,0:239)
      real*8 mutab(0:2000,0:400),stab(0:2000,0:400)
      real*8 ptab2(0:200,0:200),ttab2(0:200,0:200)
      real*8 mutab2(0:200,0:200),stab2(0:200,0:200),msttab2(0:200,0:200)
      real*8 mustab(0:2000,0:400),mustab2(0:200,0:200)
      real*8 cstab(0:2000,0:400),cstab2(0:200,0:200)
      real*8 ptab3(0:200,0:200),ttab3(0:200,0:200)
      real*8 mutab3(0:200,0:200),stab3(0:200,0:200),msttab3(0:200,0:200)
      real*8 mustab3(0:200,0:200),cstab3(0:200,0:200)
      integer stabil,antic
c
c  Type declarations for variables used in function schem.
c
      real*8 e,n,chem,et0,temp,schem,cs
      real*8 de,dn,p1,p2,p3,p4,p13,p24
c      real*8 p12,p34

c  Common-block

      common /tab/ ptab,ttab,mutab,stab,lamtab,ptab2,ttab2,mutab2,stab2,
     $     mustab,mustab2,cstab,cstab2,cstab3,ptab3,ttab3,mutab3
     $     ,stab3,mustab3
c
c  Fermi energy density in the QGP for given baryon density
c  (in units of e0).
c
      et0 = 1d0/54d0/pi2*dabs(40.5d0*pi2*n*n0*hc3)**(4d0/3d0)+B
      et0 = et0/e0/hc3
      

      antic = 0
      if (n.lt.0)then
         n=-n
         antic = 1
      end if
c
c  Calculate chemical potential (in MeV).
c  If eos=1, chemical potential is zero.
c

            
      if (eos.eq.1.or.eos.eq.2.or.eos.eq.4) then
         cs = 0.0d0
c     
c  If flag eos=[else], use tabellized equation of state.
c     

c     If e and n are both smaller than e0, n0, the nuclei are in the
c     ground state. The pressure is set to zero by hand, if flag
c     stabil = 1.
c     

      else if (eos.eq.3.or.eos.eq.5) then        
         if (e.lt.400.0d0) then

            if((e.lt.0.1d0).and.(n.lt.0.02d0)) then
               de = 0.0005d0
               dn = 0.0001d0
               p1 = cstab3(idint(e/de),idint(n/dn))
               p2 = cstab3(idint(e/de)+1,idint(n/dn))
               p3 = cstab3(idint(e/de),idint(n/dn)+1)
               p4 = cstab3(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3


               cs = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            if((e.lt.10.0d0).and.(n.lt.2.0d0).and.
     $           ((e.ge.0.1d0).or.(n.ge.0.02d0))) then
               de = 0.05d0
               dn = 0.01d0
               p1 = cstab2(idint(e/de),idint(n/dn))
               p2 = cstab2(idint(e/de)+1,idint(n/dn))
               p3 = cstab2(idint(e/de),idint(n/dn)+1)
               p4 = cstab2(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
            


               cs = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            if ((e.ge.10.0d0).or.(n.ge.2.0d0)) then 
               de = 0.5d0
               dn = 0.1d0
               
               p1 = cstab(idint(e/de),idint(n/dn))
               p2 = cstab(idint(e/de)+1,idint(n/dn))
               p3 = cstab(idint(e/de),idint(n/dn)+1)
               p4 = cstab(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
            
               cs = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
c     chem = 0.25d0*(p13+p24+p12+p34)    
            end if
         else if (e.gt.400.0d0) then
c     JS muss bestimmt werden
            cs = 1/sqrt(3d0)
c     else
c     chem = 0d0
         end if
         
      end if
      
      if ((e.le.1d0).and.(n.le.1d0).and.(stabil.eq.1)) then
         cs = 0
      end if
      
      if (antic.eq.1)then
         n=-n
         cs = cs

      end if
      return
      end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c..  FUNCTION MUPIK
c..
c..  This function-subprogram determines the pion and kaon chemical
c..  potential (for a weakly interacting gas of relativistic bosons)
c*****|****************************************************************|X

      FUNCTION mupik(temp,e,n,nB,lam,pik)
      implicit none
      INCLUDE 'defs.f'

c-----|----------------------------------------------------------------|X

      integer pik
      real*8 mupik,temp,T,n,nB,e,m,z,g,lam
      real*8 mass_pi,mass_k
      real*8 neff

      real*8 bessk
      external bessk

      PARAMETER (mass_pi=0.139d0)
      PARAMETER (mass_k=0.493d0)

c-----|----------------------------------------------------------------|X

      T=temp
      m=0.0d0
      g=1.0d0
      if(pik.eq.0) m=mass_pi
      if(pik.eq.1) m=mass_k
  
      if(pik.eq.0) g=3.0d0
      if(pik.eq.1) g=2.0d0

      if (T.lt.0.05) then
        mupik=0.0d0
        return    
      endif 

      z=m/T

      neff=n
      if (eos.eq.3.OR.eos.eq.5) neff=lam*n

c... Ideal gas in Boltzmann limit [J. Sollfrank, P. Koch, U. Heinz, 
c... Z.Phys. C52 (1991) 593-609, Equation (33); See also B. Kämpfer,
c... P.Koch, O. P. Pavlenko, Phys.Rev. C49 (1994) 1132-1138, 
c... Equation (4)]

c      mupik=T*dlog(2.0d0*pi**2*(neff*hqc**3)/
c     &               (g*T*(m**2)*bessk(2,z)))

c... It is nearly identical to the realtivistic ideal gas ansatz from 
c... [W. Greiner, N. Neise, H. Stoecker,"Thermodynamik und 
c... Statistische Mechanik", Frankfurt 1993; Equation (12.18)]:

      mupik=-1.0d0*T*dlog(T*m**2/(2*pi**2*(neff/g)*hqc**3)
     &                       *bessk(2,z))

c... Alternative Formula (do not use)

c      mupik=-3.0d0/2.0d0*T*
c     &      dlog(T*m/(2.0d0*pi*(neff/g*hqc**3)**(2.0d0/3.0d0)))
  
c      write(0,*)'in MUPIK: temp,n,pik',temp,n,pik ! Debug only
c      write(0,*)'z, bessk_2(z)',z,bessk(2,z) ! Debug only
c      write(0,*)'mu_pi/k =',mupik ! Debug only 
 
      return

      END

c*****|****************************************************************|X
