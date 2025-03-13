c....&|...Version 1.0...last change...14.08.2012...Stephan Endres......|X
c------------------------------------------------------------------------
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE TRANSFORM
c..
c.. The subroutine calculates local rest frame quantities from
c.. laboratory frame quantities. We hereby use the definition of 
c.. Eckart, where the rest frame is defined via particle (baryon)
c.. 4-flow N^mu: N^i=0 (i=1,2,3)
c****&|****************************************************************|X

      SUBROUTINE transform(noe,grd,grd_z,dx,timesteps,time,nopart,p0,
     &    n0,npi,nk,jx,jy,jz,px,py,pz,tau,p0_lrf,n0_lrf,npi_lrf,nk_lrf,
     &    jx_lrf,jy_lrf,jz_lrf,cells_bar,cells_mes_nobar,gam,vol,volce,
     &    taa,tab,tac,tad,tbb,tcc,tdd,tbc,tbd,tcd,edens_lrf,rhoeff,
     &    rhoe_lrf,cells_lt_minpart)
   
      implicit none
      include 'defs.f'


c. variables(from subroutine input)
      integer noe,grd,grd_z,timesteps
      real*8 time,timegev,vol,dx,dt

c. variables(from subroutine grid)
      real*8 nopart(timemax,grd,grd,grd_z)
      real*8 p0(timemax,grd,grd,grd_z),n0(timemax,grd,grd,grd_z)
      real*8 npi(timemax,grd,grd,grd_z),nk(timemax,grd,grd,grd_z)
      real*8 jx(timemax,grd,grd,grd_z),jy(timemax,grd,grd,grd_z)
      real*8 jz(timemax,grd,grd,grd_z)  
      real*8 px(timemax,grd,grd,grd_z),py(timemax,grd,grd,grd_z)
      real*8 pz(timemax,grd,grd,grd_z)  
      real*8 rhoeff(timemax,grd,grd,grd_z)

c.. quantities in local rest frame       
      real*8 tau(timemax,grd,grd,grd_z)
      real*8 p0_lrf(timemax,grd,grd,grd_z)
      real*8 n0_lrf(timemax,grd,grd,grd_z)
      real*8 npi_lrf(timemax,grd,grd,grd_z)
      real*8 nk_lrf(timemax,grd,grd,grd_z)
      real*8 jx_lrf(timemax,grd,grd,grd_z)
      real*8 jy_lrf(timemax,grd,grd,grd_z)
      real*8 jz_lrf(timemax,grd,grd,grd_z) 

      real*8 volce(timemax,grd,grd,grd_z)
      real*8 dt_lrf(timemax,grd,grd,grd_z)
      real*8 dx_lrf(timemax,grd,grd,grd_z)
      real*8 dy_lrf(timemax,grd,grd,grd_z)
      real*8 dz_lrf(timemax,grd,grd,grd_z) 
      
      real*8 edens_lrf(timemax,grd,grd,grd_z)    
      real*8 rhoe_lrf(timemax,grd,grd,grd_z) 

c.. general
      real*8 betax,betay,betaz,beta
      real*8 gamma 
      real*8 gam(timemax,grd,grd,grd_z)
      real*8 etrans1,etrans2,etrans3,edens
      real*8 tab_tr,tab_tr2,tac_tr,tad_tr
      real*8 taa_tr,tbb_tr,tcc_tr,tdd_tr
      real*8 taa_tr2,tbb_tr2,tcc_tr2,tdd_tr2
      real*8 tbc_tr,tbd_tr,tcd_tr

      real*8 taa(timemax,grd,grd,grd_z)
      real*8 tab(timemax,grd,grd,grd_z)
      real*8 tac(timemax,grd,grd,grd_z)
      real*8 tad(timemax,grd,grd,grd_z)
      real*8 tbb(timemax,grd,grd,grd_z)
      real*8 tcc(timemax,grd,grd,grd_z)
      real*8 tdd(timemax,grd,grd,grd_z)
      real*8 tbc(timemax,grd,grd,grd_z)
      real*8 tbd(timemax,grd,grd,grd_z)
      real*8 tcd(timemax,grd,grd,grd_z)

      real*8 px_lrf,py_lrf,pz_lrf

      real*8 p0_total(timemax)
      real*8 n0_total(timemax)
      real*8 p0_lor(timemax)
      real*8 n0_lor(timemax)

c.. counter for filled / empty cells
      integer cells_bar,cells_mes_nobar,cells_lt_minpart
      integer cmesnobar_time,cbar_time

c.. anisotropy
      real*8 P_trvse,P_par
      real*8 aniso,x_aniso,dr_dx,r_aniso

c.. external function
      real*8 ranff
      external ranff

c.. other
      integer h,i,j,k,hh
      integer counter,counterx
      real*8 percent
      real*8 ddx,ddy,ddz
      character dot
      real*8 minpart,minbar
      real*8 beta_t

c-----------------------------------------------------------------------X
c. Set everything to zero
      counter=0
      counterx=1
      dot='.'      

      cells_bar=0
      cells_mes_nobar=0

      do 976 hh=1,timemax
      n0_total(hh)=0.0d0
      p0_total(hh)=0.0d0     
      n0_lor(hh)=0.0d0
      p0_lor(hh)=0.0d0  
 976  continue
 
      write(0,*)'Total events',noe
      write(0,*)'Subroutine TRANSFORM started'   
 
c-----------------------------------------------------------------------X
c. LOOP OVER CELLS

      do 314 h=1,timesteps
       cbar_time=0
       cmesnobar_time=0
c.. Printing Progress of Subroutine ..............!
       write(0,'(a1)',advance='NO'),dot           !
c       write(0,'(a1)'),dot                       !
       counter=counter+1                          !
       if(counter.eq.10) then                     ! 
        percent=100d0*counterx*counter/timesteps  !
        counterx=counterx+1                       !
        counter=0                                 !
        write(0,'(f4.0)',advance='NO'),percent    !
c        write(0,'(f4.0)'),percent                !
        write(0,*),'%'                            !
       end if                                     !
c.................................................!
       do 313 i=1,grd
        do 312 j=1,grd
         do 311 k=1,grd_z
c          write(0,*),h,i,j,k ! Debug only
     
 
          if(n0(h,i,j,k).eq.0d0.AND.nopart(h,i,j,k).gt.0d0) then
           cells_mes_nobar=cells_mes_nobar+1
           cmesnobar_time=cmesnobar_time+1
          endif

c... skip the rest, if baryon current is zero in cell 
          if(n0(h,i,j,k).eq.0.0d0) cycle 

c... For PHOTON emission
c... Special treatment for cells with too little content
          if(rates.eq.-1) then
           minpart=9.0d0
           minbar=7.0d0/vol
           beta_t=sqrt((jx(h,i,j,k)/n0(h,i,j,k))**2
     &                +(jy(h,i,j,k)/n0(h,i,j,k))**2)
           if(nopart(h,i,j,k).lt.minpart.OR.
     &       (n0(h,i,j,k).lt.minbar.AND.beta_t.gt.0.4d0)) then
            if(i.gt.2.AND.i.lt.(grd-1).AND.j.gt.2.AND.j.lt.(grd-1).AND.
     &         K.gt.2.AND.K.lt.(grd_z-1)) then
              nopart(h,i,j,k)=(nopart(h,i,j,k)+nopart(h,i+1,j,k)+nopart(h,i-1,j,k)+nopart(h,i,j+1,k)
     &                        +nopart(h,i,j-1,k)+nopart(h,i,j,k-1)+nopart(h,i,j,k-1))/7.0d0                             
              n0(h,i,j,k)=(n0(h,i,j,k)+n0(h,i+1,j,k)+n0(h,i-1,j,k)+n0(h,i,j+1,k)+n0(h,i,j-1,k)
     &                    +n0(h,i,j,k-1)+n0(h,i,j,k-1))/7.0d0    
            endif
            if(nopart(h,i,j,k).lt.minpart.OR.
     &        (n0(h,i,j,k).lt.minbar.AND.beta_t.gt.0.4d0)) then  
              cells_lt_minpart=cells_lt_minpart+1
              cycle
            else
            jx(h,i,j,k)=(jx(h,i,j,k)+jx(h,i+1,j,k)+jx(h,i-1,j,k)+jx(h,i,j+1,k)
     &                  +jx(h,i,j-1,k)+jx(h,i,j,k-1)+jx(h,i,j,k-1))/7.0d0
            jy(h,i,j,k)=(jy(h,i,j,k)+jy(h,i+1,j,k)+jy(h,i-1,j,k)+jy(h,i,j+1,k)
     &                  +jy(h,i,j-1,k)+jy(h,i,j,k-1)+jy(h,i,j,k-1))/7.0d0
            jz(h,i,j,k)=(jz(h,i,j,k)+jz(h,i+1,j,k)+jz(h,i-1,j,k)+jz(h,i,j+1,k)
     &                  +jz(h,i,j-1,k)+jz(h,i,j,k-1)+jz(h,i,j,k-1))/7.0d0 

            p0(h,i,j,k)=(p0(h,i,j,k)+p0(h,i+1,j,k)+p0(h,i-1,j,k)+p0(h,i,j+1,k)
     &                  +p0(h,i,j-1,k)+p0(h,i,j,k-1)+p0(h,i,j,k-1))/7.0d0
            px(h,i,j,k)=(px(h,i,j,k)+px(h,i+1,j,k)+px(h,i-1,j,k)+px(h,i,j+1,k)
     &                  +px(h,i,j-1,k)+px(h,i,j,k-1)+px(h,i,j,k-1))/7.0d0
            py(h,i,j,k)=(py(h,i,j,k)+py(h,i+1,j,k)+py(h,i-1,j,k)+py(h,i,j+1,k)
     &                  +py(h,i,j-1,k)+py(h,i,j,k-1)+py(h,i,j,k-1))/7.0d0
            pz(h,i,j,k)=(pz(h,i,j,k)+pz(h,i+1,j,k)+pz(h,i-1,j,k)+pz(h,i,j+1,k)
     &                  +pz(h,i,j-1,k)+pz(h,i,j,k-1)+pz(h,i,j,k-1))/7.0d0

            taa(h,i,j,k)=(taa(h,i,j,k)+taa(h,i+1,j,k)+taa(h,i-1,j,k)+taa(h,i,j+1,k)
     &                   +taa(h,i,j-1,k)+taa(h,i,j,k-1)+taa(h,i,j,k-1))/7.0d0
            tab(h,i,j,k)=(tab(h,i,j,k)+tab(h,i+1,j,k)+tab(h,i-1,j,k)+tab(h,i,j+1,k)
     &                   +tab(h,i,j-1,k)+tab(h,i,j,k-1)+tab(h,i,j,k-1))/7.0d0
            tac(h,i,j,k)=(tac(h,i,j,k)+tac(h,i+1,j,k)+tac(h,i-1,j,k)+tac(h,i,j+1,k)
     &                   +tac(h,i,j-1,k)+tac(h,i,j,k-1)+tac(h,i,j,k-1))/7.0d0
            tad(h,i,j,k)=(tad(h,i,j,k)+tad(h,i+1,j,k)+tad(h,i-1,j,k)+tad(h,i,j+1,k)
     &                   +tad(h,i,j-1,k)+tad(h,i,j,k-1)+tad(h,i,j,k-1))/7.0d0

            tbb(h,i,j,k)=(tbb(h,i,j,k)+tbb(h,i+1,j,k)+tbb(h,i-1,j,k)+tbb(h,i,j+1,k)
     &                   +tbb(h,i,j-1,k)+tbb(h,i,j,k-1)+tbb(h,i,j,k-1))/7.0d0
            tcc(h,i,j,k)=(tcc(h,i,j,k)+tcc(h,i+1,j,k)+tcc(h,i-1,j,k)+tcc(h,i,j+1,k)
     &                   +tcc(h,i,j-1,k)+tcc(h,i,j,k-1)+tcc(h,i,j,k-1))/7.0d0
            tdd(h,i,j,k)=(tdd(h,i,j,k)+tdd(h,i+1,j,k)+tdd(h,i-1,j,k)+tdd(h,i,j+1,k)
     &                   +tdd(h,i,j-1,k)+tdd(h,i,j,k-1)+tdd(h,i,j,k-1))/7.0d0
          
            tbc(h,i,j,k)=(tbc(h,i,j,k)+tbc(h,i+1,j,k)+tbc(h,i-1,j,k)+tbc(h,i,j+1,k)
     &                   +tbc(h,i,j-1,k)+tbc(h,i,j,k-1)+tbc(h,i,j,k-1))/7.0d0
            tbd(h,i,j,k)=(tbd(h,i,j,k)+tbd(h,i+1,j,k)+tbd(h,i-1,j,k)+tbd(h,i,j+1,k)
     &                   +tbd(h,i,j-1,k)+tbd(h,i,j,k-1)+tbd(h,i,j,k-1))/7.0d0
            tcd(h,i,j,k)=(tcd(h,i,j,k)+tcd(h,i+1,j,k)+tcd(h,i-1,j,k)+tcd(h,i,j+1,k)
     &                   +tcd(h,i,j-1,k)+tcd(h,i,j,k-1)+tcd(h,i,j,k-1))/7.0d0

            rhoeff(h,i,j,k)=(rhoeff(h,i,j,k)+rhoeff(h,i+1,j,k)+rhoeff(h,i-1,j,k)+rhoeff(h,i,j+1,k)
     &                      +rhoeff(h,i,j-1,k)+rhoeff(h,i,j,k-1)+rhoeff(h,i,j,k-1))/7.0d0  
            npi(h,i,j,k)=(npi(h,i,j,k)+npi(h,i+1,j,k)+npi(h,i-1,j,k)+npi(h,i,j+1,k)+npi(h,i,j-1,k)
     &                   +npi(h,i,j,k-1)+npi(h,i,j,k-1))/7.0d0
            nk(h,i,j,k)=(nk(h,i,j,k)+nk(h,i+1,j,k)+nk(h,i-1,j,k)+nk(h,i,j+1,k)+nk(h,i,j-1,k)
     &                  +nk(h,i,j,k-1)+nk(h,i,j,k-1))/7.0d0    
            endif
           endif
          endif !rates.eq.-1

          cells_bar=cells_bar+1
          cbar_time=cbar_time+1
          
c... REMARK ON UNITS: [p_mu]=GeV
c... Normalize cells to number of events

          nopart(h,i,j,k)=nopart(h,i,j,k)/dble(noe)
          n0(h,i,j,k)=n0(h,i,j,k)/dble(noe)
          rhoeff(h,i,j,k)=rhoeff(h,i,j,k)/dble(noe)
          npi(h,i,j,k)=npi(h,i,j,k)/dble(noe)
          nk(h,i,j,k)=nk(h,i,j,k)/dble(noe)

          p0(h,i,j,k)=p0(h,i,j,k)/dble(noe)

          px(h,i,j,k)=px(h,i,j,k)/dble(noe)
          py(h,i,j,k)=py(h,i,j,k)/dble(noe)
          pz(h,i,j,k)=pz(h,i,j,k)/dble(noe)
          
          jx(h,i,j,k)=jx(h,i,j,k)/dble(noe)
          jy(h,i,j,k)=jy(h,i,j,k)/dble(noe)
          jz(h,i,j,k)=jz(h,i,j,k)/dble(noe)

c... energy momentum tensor

c          write(0,*)'t00,t01,t02,t03',taa(h,i,j,k),tab(h,i,j,k),
c     &               tac(h,i,j,k),tad(h,i,j,k) ! Debug only
c          write(0,*)'t11,t12,t13',tbb(h,i,j,k),tbc(h,i,j,k),
c     &               tbd(h,i,j,k) ! Debug only
c          write(0,*)'t22,t23',tcc(h,i,j,k),tcd(h,i,j,k) ! Debug only
c          write(0,*)'t33',tdd(h,i,j,k) ! Debug only

          taa(h,i,j,k)=taa(h,i,j,k)/dble(noe)
          tab(h,i,j,k)=tab(h,i,j,k)/dble(noe)
          tac(h,i,j,k)=tac(h,i,j,k)/dble(noe)
          tad(h,i,j,k)=tad(h,i,j,k)/dble(noe)

          tbb(h,i,j,k)=tbb(h,i,j,k)/dble(noe)
          tcc(h,i,j,k)=tcc(h,i,j,k)/dble(noe)
          tdd(h,i,j,k)=tdd(h,i,j,k)/dble(noe)

          tbc(h,i,j,k)=tbc(h,i,j,k)/dble(noe)
          tbd(h,i,j,k)=tbd(h,i,j,k)/dble(noe)
          tcd(h,i,j,k)=tcd(h,i,j,k)/dble(noe)


c          write(0,*)'t00,t01,t02,t03',taa(h,i,j,k),tab(h,i,j,k),
c     &               tac(h,i,j,k),tad(h,i,j,k) ! Debug only
c          write(0,*)'t11,t12,t13',tbb(h,i,j,k),tbc(h,i,j,k),
c     &               tbd(h,i,j,k) ! Debug only
c          write(0,*)'t22,t23',tcc(h,i,j,k),tcd(h,i,j,k) ! Debug only
c          write(0,*)'t33',tdd(h,i,j,k) ! Debug only

c-----------------------------------------------------------------------X

c.. calculate BETA, GAMMA and TAU 

          betax=jx(h,i,j,k)/n0(h,i,j,k)
          betay=jy(h,i,j,k)/n0(h,i,j,k)
          betaz=jz(h,i,j,k)/n0(h,i,j,k)
c..In rare cases with only few particles and in the presence of anti-...!
c..baryons "betaz" might get unphysically high values...................!
          if(dabs(betax).ge.1.0d0) then                                 !
            betax=ranff(seedinit)*0.1d0*dsign(1.0d0,betax)              !
            jx(h,i,j,k)=n0(h,i,j,k)*betax                               !
          end if                                                        !
          if(dabs(betay).ge.1.0d0) then                                 !
            betay=ranff(seedinit)*0.1d0*dsign(1.0d0,betay)              !
            jy(h,i,j,k)=n0(h,i,j,k)*betay                               !
          end if                                                        !
          if(dabs(betaz).ge.1.0d0) then                                 !
            betaz=ranff(seedinit)*0.9d0*dsign(1.0d0,betaz)              !
            jz(h,i,j,k)=n0(h,i,j,k)*betaz                               !
          end if                                                        !
c.......................................................................!
c.FOR PHOTONS: High trasvese flow of single cells strongly blue shifts..!
c............. the photon results.......................................!
          if(rates.eq.-1.AND.n0(h,i,j,k).lt.(30.0d0/vol/dble(noe))) then!
c         if(n0(h,i,j,k).lt.(30.0d0/vol/dble(noe))) then                !
           if(dabs(betax).ge.0.45d0) then                               !
            betax=0.45d0*dsign(1.0d0,betax)                             !
            jx(h,i,j,k)=n0(h,i,j,k)*betax                               !
           end if                                                       !
           if(dabs(betay).ge.0.45d0) then                               !
            betay=0.45d0*dsign(1.0d0,betay)                             !
            jy(h,i,j,k)=n0(h,i,j,k)*betay                               !
           end if                                                       !
          end if                                                        !
c.......................................................................!



 777      continue
          beta=sqrt((betax**2)+(betay**2)+(betaz**2))
          if(beta.gt.1.0d0) then
            betax=betax*0.5d0
            jx(h,i,j,k)=jx(h,i,j,k)*0.5d0
            betay=betay*0.5d0
            jy(h,i,j,k)=jy(h,i,j,k)*0.5d0
            goto 777 
          endif

          gamma=1d0/sqrt(1d0-beta**2)
          gam(h,i,j,k)=gamma

c-----------------------------------------------------------------------X

c.. calculate local rest frame quantities
         p0_lrf(h,i,j,k)=gamma*(p0(h,i,j,k)-betax*px(h,i,j,k)-
     &          betay*py(h,i,j,k)-betaz*pz(h,i,j,k))

         n0_lrf(h,i,j,k)=gamma*(n0(h,i,j,k)-betax*jx(h,i,j,k)-
     &          betay*jy(h,i,j,k)-betaz*jz(h,i,j,k))


         npi_lrf(h,i,j,k)=npi(h,i,j,k)/gamma
         nk_lrf(h,i,j,k)=nk(h,i,j,k)/gamma
         rhoe_lrf(h,i,j,k)=rhoeff(h,i,j,k)/gamma


         jx_lrf(h,i,j,k)=(-gamma*betax*n0(h,i,j,k))+
     &          (1+(gamma-1)*((betax**2)/(beta**2)))*jx(h,i,j,k)+
     &          ((gamma-1)*((betax*betay)/(beta**2)))*jy(h,i,j,k)+
     &          ((gamma-1)*((betax*betaz)/(beta**2)))*jz(h,i,j,k)

         jy_lrf(h,i,j,k)=(-gamma*betay*n0(h,i,j,k))+
     &          ((gamma-1)*((betax*betay)/(beta**2)))*jx(h,i,j,k)+
     &          (1+(gamma-1)*((betay**2)/(beta**2)))*jy(h,i,j,k)+
     &          ((gamma-1)*((betay*betaz)/(beta**2)))*jz(h,i,j,k)

         jz_lrf(h,i,j,k)=(-gamma*betaz*n0(h,i,j,k))+
     &          ((gamma-1)*((betax*betaz)/(beta**2)))*jx(h,i,j,k)+
     &          ((gamma-1)*((betay*betaz)/(beta**2)))*jy(h,i,j,k)+
     &          (1+(gamma-1)*((betaz**2)/(beta**2)))*jz(h,i,j,k)


         px_lrf=(-gamma*betax*p0(h,i,j,k))+
     &          (1+(gamma-1)*((betax**2)/(beta**2)))*px(h,i,j,k)+
     &          ((gamma-1)*((betax*betay)/(beta**2)))*py(h,i,j,k)+
     &          ((gamma-1)*((betax*betaz)/(beta**2)))*pz(h,i,j,k)

         py_lrf=(-gamma*betay*p0(h,i,j,k))+
     &          ((gamma-1)*((betax*betay)/(beta**2)))*px(h,i,j,k)+
     &          (1+(gamma-1)*((betay**2)/(beta**2)))*py(h,i,j,k)+
     &          ((gamma-1)*((betay*betaz)/(beta**2)))*pz(h,i,j,k)

         pz_lrf=(-gamma*betaz*p0(h,i,j,k))+
     &          ((gamma-1)*((betax*betaz)/(beta**2)))*px(h,i,j,k)+
     &          ((gamma-1)*((betay*betaz)/(beta**2)))*py(h,i,j,k)+
     &          (1+(gamma-1)*((betaz**2)/(beta**2)))*pz(h,i,j,k)


         tau(h,i,j,k)=time/gamma
         volce(h,i,j,k)=((dx)**3)*gamma 

c-----------------------------------------------------------------------X
c... Transformation of components of T_mu,nu 

      taa_tr=gamma*gamma*(taa(h,i,j,k)-betax*tab(h,i,j,k)-
     &         betay*tac(h,i,j,k)-betaz*tad(h,i,j,k))-
     &        betax*gamma*gamma*(tab(h,i,j,k)-betax*tbb(h,i,j,k)
     &         -betay*tbc(h,i,j,k)-betaz*tbd(h,i,j,k))-
     &        betay*gamma*gamma*(tac(h,i,j,k)-betax*tbc(h,i,j,k)
     &         -betay*tcc(h,i,j,k)-betaz*tcd(h,i,j,k))-
     &        betaz*gamma*gamma*(tad(h,i,j,k)-betax*tbd(h,i,j,k)
     &         -betay*tcd(h,i,j,k)-betaz*tdd(h,i,j,k))

      tab_tr=-betax*gamma*gamma*(taa(h,i,j,k)-betax*tab(h,i,j,k)
     &          -betay*tac(h,i,j,k)-betaz*tad(h,i,j,k))+
     &        (1.0d0+(gamma-1.0d0)*(betax**2)/(beta**2))
     &          *gamma*(tab(h,i,j,k)-betax*tbb(h,i,j,k)-
     &          betay*tbc(h,i,j,k)-betaz*tbd(h,i,j,k))+
     &        ((gamma-1.0d0)*betax*betay/(beta**2))*gamma*
     &          (tac(h,i,j,k)-betax*tbc(h,i,j,k)-betay*tcc(h,i,j,k)
     &          -betaz*tcd(h,i,j,k))+
     &        ((gamma-1.0d0)*betax*betaz/(beta**2))*gamma*
     &          (tad(h,i,j,k)-betax*tbd(h,i,j,k)-betay*tcd(h,i,j,k)
     &          -betaz*tdd(h,i,j,k))


c      tab_tr2=gamma*(-betax*gamma*taa(h,i,j,k)+(1.0d0+(gamma-1.0d0)*(betax**2)/(beta**2))*tab(h,i,j,k)+
c     &                  ((gamma-1.0d0)*betax*betay/(beta**2))*tac(h,i,j,k)+((gamma-1.0d0)*betax*betaz/(beta**2))*tad(h,i,j,k))-
c     &        betax*gamma*(-betax*gamma*tab(h,i,j,k)+(1.0d0+(gamma-1.0d0)*(betax**2)/(beta**2))*tbb(h,i,j,k)+
c     &                  ((gamma-1.0d0)*betax*betay/(beta**2))*tbc(h,i,j,k)+((gamma-1.0d0)*betax*betaz/(beta**2))*tbd(h,i,j,k))-
c     &        betay*gamma*(-betax*gamma*tac(h,i,j,k)+(1.0d0+(gamma-1.0d0)*(betax**2)/(beta**2))*tbc(h,i,j,k)+
c     &                  ((gamma-1.0d0)*betax*betay/(beta**2))*tcc(h,i,j,k)+((gamma-1.0d0)*betax*betaz/(beta**2))*tcd(h,i,j,k))-
c     &        betaz*gamma*(-betax*gamma*tad(h,i,j,k)+(1.0d0+(gamma-1.0d0)*(betax**2)/(beta**2))*tbd(h,i,j,k)+
c     &                  ((gamma-1.0d0)*betax*betay/(beta**2))*tcd(h,i,j,k)+((gamma-1.0d0)*betax*betaz/(beta**2))*tdd(h,i,j,k))


       tac_tr=-betay*gamma*gamma*(taa(h,i,j,k)-betax*tab(h,i,j,k)-betay*tac(h,i,j,k)-betaz*tad(h,i,j,k))+
     &        ((gamma-1.0d0)*(betax*betay)/(beta**2))*gamma*(tab(h,i,j,k)-betax*tbb(h,i,j,k)-
     &           betay*tbc(h,i,j,k)-betaz*tbd(h,i,j,k))+
     &        (1.0d0+(gamma-1.0d0)*betay**2/(beta**2))*gamma*(tac(h,i,j,k)-betax*tbc(h,i,j,k)-
     &           betay*tcc(h,i,j,k)-betaz*tcd(h,i,j,k))+
     &        ((gamma-1.0d0)*betay*betaz/(beta**2))*gamma*(tad(h,i,j,k)-betax*tbd(h,i,j,k)-betay*tcd(h,i,j,k)-betaz*tdd(h,i,j,k))


       tad_tr=-betaz*gamma*gamma*(taa(h,i,j,k)-betax*tab(h,i,j,k)-betay*tac(h,i,j,k)-betaz*tad(h,i,j,k))+
     &        ((gamma-1.0d0)*(betax*betaz)/(beta**2))*gamma*(tab(h,i,j,k)-betax*tbb(h,i,j,k)-
     &           betay*tbc(h,i,j,k)-betaz*tbd(h,i,j,k))+
     &        ((gamma-1.0d0)*betay*betaz/(beta**2))*gamma*(tac(h,i,j,k)-betax*tbc(h,i,j,k)-
     &           betay*tcc(h,i,j,k)-betaz*tcd(h,i,j,k))+
     &        (1.0d0+(gamma-1.0d0)*betaz**2/(beta**2))*gamma*(tad(h,i,j,k)-betax*tbd(h,i,j,k)-betay*tcd(h,i,j,k)-betaz*tdd(h,i,j,k))


      tbb_tr=-betax*gamma*(-betax*gamma*taa(h,i,j,k)+(1.0d0+(gamma-1.0d0)*betax**2/(beta**2))*tab(h,i,j,k)+
     &       ((gamma-1.0d0)*betax*betay/(beta**2))*tac(h,i,j,k)+((gamma-1.0d0)*betax*betaz/(beta**2))*tad(h,i,j,k))+
     &        (1.0d0+(gamma-1.0d0)*betax**2/(beta**2))*(-betax*gamma*tab(h,i,j,k)+
     &        (1.0d0+(gamma-1.0d0)*betax**2/(beta**2))*tbb(h,i,j,k)+
     &       ((gamma-1.0d0)*betax*betay/(beta**2))*tbc(h,i,j,k)+((gamma-1.0d0)*betax*betaz/(beta**2))*tbd(h,i,j,k))+
     &        ((gamma-1.0d0)*betax*betay/(beta**2))*(-betax*gamma*tac(h,i,j,k)+
     &        (1.0d0+(gamma-1.0d0)*betax**2/(beta**2))*tbc(h,i,j,k)+
     &       ((gamma-1.0d0)*betax*betay/(beta**2))*tcc(h,i,j,k)+((gamma-1.0d0)*betax*betaz/(beta**2))*tcd(h,i,j,k))+
     &        ((gamma-1.0d0)*betax*betaz/(beta**2))*(-betax*gamma*tad(h,i,j,k)+
     &        (1.0d0+(gamma-1.0d0)*betax**2/(beta**2))*tbd(h,i,j,k)+
     &       ((gamma-1.0d0)*betax*betay/(beta**2))*tcd(h,i,j,k)+((gamma-1.0d0)*betax*betaz/(beta**2))*tdd(h,i,j,k))
  

      tcc_tr=-betay*gamma*(-betay*gamma*taa(h,i,j,k)+
     &        ((gamma-1.0d0)*betax*betay/(beta**2))*tab(h,i,j,k)+
     &        (1.0d0+(gamma-1.0d0)*betay**2/(beta**2))*tac(h,i,j,k)+
     &        ((gamma-1.0d0)*betay*betaz/(beta**2))*tad(h,i,j,k))+
     &       ((gamma-1.0d0)*betax*betay/(beta**2))*
     &        (-betay*gamma*tab(h,i,j,k)+
     &        ((gamma-1.0d0)*betax*betay/(beta**2))*tbb(h,i,j,k)+
     &        (1.0d0+(gamma-1.0d0)*betay**2/(beta**2))*tbc(h,i,j,k)+
     &        ((gamma-1.0d0)*betay*betaz/(beta**2))*tbd(h,i,j,k))+
     &       (1.0d0+(gamma-1.0d0)*betay**2/(beta**2))*
     &        (-betay*gamma*tac(h,i,j,k)+
     &        ((gamma-1.0d0)*betax*betay/(beta**2))*tbc(h,i,j,k)+
     &        (1.0d0+(gamma-1.0d0)*betay**2/(beta**2))*tcc(h,i,j,k)+
     &        ((gamma-1.0d0)*betay*betaz/(beta**2))*tcd(h,i,j,k))+
     &       ((gamma-1.0d0)*betay*betaz/(beta**2))*
     &        (-betay*gamma*tad(h,i,j,k)+
     &        ((gamma-1.0d0)*betax*betay/(beta**2))*tbd(h,i,j,k)+
     &        (1.0d0+(gamma-1.0d0)*betay**2/(beta**2))*tcd(h,i,j,k)+
     &        ((gamma-1.0d0)*betay*betaz/(beta**2))*tdd(h,i,j,k))

      tdd_tr=-betaz*gamma*(-betaz*gamma*taa(h,i,j,k)+
     &        ((gamma-1.0d0)*betax*betaz/(beta**2))*tab(h,i,j,k)+
     &        ((gamma-1.0d0)*betay*betaz/(beta**2))*tac(h,i,j,k)+
     &        (1.0d0+(gamma-1.0d0)*betaz**2/(beta**2))*tad(h,i,j,k))+
     &       ((gamma-1.0d0)*betax*betaz/(beta**2))*
     &        (-betaz*gamma*tab(h,i,j,k)+
     &        ((gamma-1.0d0)*betax*betaz/(beta**2))*tbb(h,i,j,k)+
     &        ((gamma-1.0d0)*betay*betaz/(beta**2))*tbc(h,i,j,k)+
     &        (1.0d0+(gamma-1.0d0)*betaz**2/(beta**2))*tbd(h,i,j,k))+
     &       ((gamma-1.0d0)*betay*betaz/(beta**2))*
     &        (-betaz*gamma*tac(h,i,j,k)+
     &        ((gamma-1.0d0)*betax*betaz/(beta**2))*tbc(h,i,j,k)+
     &        ((gamma-1.0d0)*betay*betaz/(beta**2))*tcc(h,i,j,k)+
     &        (1.0d0+(gamma-1.0d0)*betaz**2/(beta**2))*tcd(h,i,j,k))+
     &       (1.0d0+(gamma-1.0d0)*betaz**2/(beta**2))*
     &        (-betaz*gamma*tad(h,i,j,k)+
     &        ((gamma-1.0d0)*betax*betaz/(beta**2))*tbd(h,i,j,k)+
     &        ((gamma-1.0d0)*betay*betaz/(beta**2))*tcd(h,i,j,k)+
     &        (1.0d0+(gamma-1.0d0)*betaz**2/(beta**2))*tdd(h,i,j,k))

      tbc_tr=-betax*gamma*(-betay*gamma*taa(h,i,j,k)+
     &        ((gamma-1.0d0)*betax*betay/(beta**2))*tab(h,i,j,k)+
     &        (1.0d0+(gamma-1.0d0)*betay**2/(beta**2))*tac(h,i,j,k)+
     &        ((gamma-1.0d0)*betay*betaz/(beta**2))*tad(h,i,j,k))+
     &       (1.0d0+(gamma-1.0d0)*betax**2/(beta**2))*
     &        (-betay*gamma*tab(h,i,j,k)+
     &        ((gamma-1.0d0)*betax*betay/(beta**2))*tbb(h,i,j,k)+
     &        (1.0d0+(gamma-1.0d0)*betay**2/(beta**2))*tbc(h,i,j,k)+
     &        ((gamma-1.0d0)*betay*betaz/(beta**2))*tbd(h,i,j,k))+
     &       ((gamma-1.0d0)*betax*betay/(beta**2))*
     &        (-betay*gamma*tac(h,i,j,k)+
     &        ((gamma-1.0d0)*betax*betay/(beta**2))*tbc(h,i,j,k)+
     &        (1.0d0+(gamma-1.0d0)*betay**2/(beta**2))*tcc(h,i,j,k)+
     &        ((gamma-1.0d0)*betay*betaz/(beta**2))*tcd(h,i,j,k))+
     &       ((gamma-1.0d0)*betax*betaz/(beta**2))*
     &        (-betay*gamma*tad(h,i,j,k)+
     &        ((gamma-1.0d0)*betax*betay/(beta**2))*tbd(h,i,j,k)+
     &        (1.0d0+(gamma-1.0d0)*betay**2/(beta**2))*tcd(h,i,j,k)+
     &        ((gamma-1.0d0)*betay*betaz/(beta**2))*tdd(h,i,j,k))

      tbd_tr=-betax*gamma*(-betaz*gamma*taa(h,i,j,k)+
     &        ((gamma-1.0d0)*betax*betaz/(beta**2))*tab(h,i,j,k)+
     &        ((gamma-1.0d0)*betay*betaz/(beta**2))*tac(h,i,j,k)+
     &        (1.0d0+(gamma-1.0d0)*betaz**2/(beta**2))*tad(h,i,j,k))+
     &       (1.0d0+(gamma-1.0d0)*betax**2/(beta**2))*
     &        (-betaz*gamma*tab(h,i,j,k)+
     &        ((gamma-1.0d0)*betax*betaz/(beta**2))*tbb(h,i,j,k)+
     &        ((gamma-1.0d0)*betay*betaz/(beta**2))*tbc(h,i,j,k)+
     &        (1.0d0+(gamma-1.0d0)*betaz**2/(beta**2))*tbd(h,i,j,k))+
     &       ((gamma-1.0d0)*betax*betay/(beta**2))*
     &        (-betaz*gamma*tac(h,i,j,k)+
     &        ((gamma-1.0d0)*betax*betaz/(beta**2))*tbc(h,i,j,k)+
     &        ((gamma-1.0d0)*betay*betaz/(beta**2))*tcc(h,i,j,k)+
     &        (1.0d0+(gamma-1.0d0)*betaz**2/(beta**2))*tcd(h,i,j,k))+
     &       ((gamma-1.0d0)*betax*betaz/(beta**2))*
     &        (-betaz*gamma*tad(h,i,j,k)+
     &        ((gamma-1.0d0)*betax*betaz/(beta**2))*tbd(h,i,j,k)+
     &        ((gamma-1.0d0)*betay*betaz/(beta**2))*tcd(h,i,j,k)+
     &        (1.0d0+(gamma-1.0d0)*betaz**2/(beta**2))*tdd(h,i,j,k))

      tcd_tr=-betay*gamma*(-betaz*gamma*taa(h,i,j,k)+
     &        ((gamma-1.0d0)*betax*betaz/(beta**2))*tab(h,i,j,k)+
     &        ((gamma-1.0d0)*betay*betaz/(beta**2))*tac(h,i,j,k)+
     &        (1.0d0+(gamma-1.0d0)*betaz**2/(beta**2))*tad(h,i,j,k))+
     &       ((gamma-1.0d0)*betax*betay/(beta**2))*
     &        (-betaz*gamma*tab(h,i,j,k)+
     &        ((gamma-1.0d0)*betax*betaz/(beta**2))*tbb(h,i,j,k)+
     &        ((gamma-1.0d0)*betay*betaz/(beta**2))*tbc(h,i,j,k)+
     &        (1.0d0+(gamma-1.0d0)*betaz**2/(beta**2))*tbd(h,i,j,k))+
     &       (1.0d0+(gamma-1.0d0)*betay**2/(beta**2))*
     &        (-betaz*gamma*tac(h,i,j,k)+
     &        ((gamma-1.0d0)*betax*betaz/(beta**2))*tbc(h,i,j,k)+
     &        ((gamma-1.0d0)*betay*betaz/(beta**2))*tcc(h,i,j,k)+
     &        (1.0d0+(gamma-1.0d0)*betaz**2/(beta**2))*tcd(h,i,j,k))+
     &       ((gamma-1.0d0)*betay*betaz/(beta**2))*
     &        (-betaz*gamma*tad(h,i,j,k)+
     &        ((gamma-1.0d0)*betax*betaz/(beta**2))*tbd(h,i,j,k)+
     &        ((gamma-1.0d0)*betay*betaz/(beta**2))*tcd(h,i,j,k)+
     &        (1.0d0+(gamma-1.0d0)*betaz**2/(beta**2))*tdd(h,i,j,k))

c      write(0,*),taa_tr,tbb_tr,tcc_tr,tdd_tr ! Debug only


c-----------------------------------------------------------------------X
c... ANISOTROPIC CASE ...

      P_trvse=(tbb_tr+tcc_tr)/2
      P_par=tdd_tr
      aniso=P_par/P_trvse

      taa_tr2=0.0d0
      tbb_tr2=0.0d0
      tcc_tr2=0.0d0
      tdd_tr2=0.0d0

      x_aniso=P_trvse**(1.333333d0)/P_par**(1.333333d0)
c        write(0,*)'anisotr. param. x =',x_aniso !Debug only
      if(x_aniso.lt.1.0d0) then
        r_aniso=0.5d0*x_aniso**(-0.33333)*(1.0d0+(x_aniso*
     &          atanh(sqrt(abs(x_aniso-1.0d0))))/
     &          (sqrt(abs(x_aniso-1.0d0))))
      else if(x_aniso.ge.1.0d0) then
        r_aniso=0.5d0*x_aniso**(-0.33333)*(1.0d0+(x_aniso*
     &          datan(sqrt(abs(x_aniso-1.0d0))))/
     &          (sqrt(abs(x_aniso-1.0d0))))
      end if
c        write(0,*)'anisotr. fuc. r =',r_aniso !Debug only
      dr_dx=-1.0d0*((2.0d0+x_aniso)*sqrt(abs(x_aniso-1.0d0))+
     &        (-4.0d0+x_aniso)*x_aniso*
     &        atanh(sqrt(abs(x_aniso-1.0d0)))/
     &        (12.0d0*abs(x_aniso-1.0d0)**(1.5d0)*x_aniso**(1.33333d0)))
c        write(0,*)'dr/dx =',dr_dx !Debug only

cc...............................DEBUGGING.ONLY.................................!
c      if(i.eq.int(grd/2+1).AND.j.eq.int(grd/2+1).AND.k.eq.int(grd_z/2+1)) then !
c        write(72,*)'cell',h,i,j,k                                              !
c        write(72,*)'P_trv,P_par,aniso,edens =',P_trvse,P_par,aniso,taa_tr      !
c         write(72,'(I5,1x,4(E14.6,1x))')h,P_trvse,P_par,x_aniso,r_aniso        ! 
c        write(72,*)'LRF T_0_1/2/3:',tab_tr,tac_tr,tad_tr                       !
c        write(72,*)'LRF T_1_2/3:',tbc_tr,tbd_tr                                ! 
c        write(72,*)'LRF T_2_3:',tcd_tr                                         !    
c        write(72,*)'LAB T_0_1/2/3:',tab(h,i,j,k),                              !
c     &                             tac(h,i,j,k),tad(h,i,j,k)                   !
c        write(72,*)'LAB T_1_2/3:',tbc(h,i,j,k),tbd(h,i,j,k)                    ! 
c        write(72,*)'LAB T_2_3:',tcd(h,i,j,k)                                   !
c      endif                                                                    !
cc..............................................................................!

      

      if(aniso.gt.2.5d0.AND.tdd_tr.gt.1d-10) then

c...Save "old" values of matrix entries (anisotropic)
       taa_tr2=taa_tr
       tbb_tr2=tbb_tr
       tcc_tr2=tcc_tr
       tdd_tr2=tdd_tr

c...Set new isotropic values
       taa_tr=taa_tr/r_aniso
       tbb_tr=P_trvse/(r_aniso+3*x_aniso*dr_dx)
       tcc_tr=P_trvse/(r_aniso+3*x_aniso*dr_dx)
       tdd_tr=P_par/(r_aniso-6*x_aniso*dr_dx)

c       write(0,*),taa_tr,tbb_tr,tcc_tr,tdd_tr !Debug only

      end if


c       if(aniso.lt.0.5d0) then
c        x_aniso=P_trvse**(1.333333d0)/P_par**(1.333333d0)
c        write(0,*)'anisotr. param. x =',x_aniso
c        r_aniso=0.5d0*x_aniso**(-0.33333)*(1.0d0+(x_aniso*
c     &          atan(sqrt(x_aniso-1.0d0)))/
c     &          (sqrt(x_aniso-1.0d0)))
c        write(0,*)'anisotr. fuc. r =',r_aniso
c        dr_dx=((2.0d0+x_aniso)*sqrt(x_aniso-1.0d0)+
c     &        (-4.0d0+x_aniso)*x_aniso*
c     &        atan(sqrt(x_aniso-1.0d0)))/
c     &        (12.0d0*(x_aniso-1.0d0)**(1.5d0)*x_aniso**(1.33333d0))
c        write(0,*)'dr/dx =',dr_dx
c       endif

c-----------------------------------------------------------------------X

c.. Set local rest frame energy density
          edens_lrf(h,i,j,k)=taa_tr

c.. Check if particle-flow vanishes in local rest frame
c         if(jx_lrf(h,i,j,k).gt.1d-11.OR.jy_lrf(h,i,j,k).gt.1d-11.OR.
c     &      jz_lrf(h,i,j,k).gt.1d-11) then
c          write(0,*)h,i,j,k
c          write(0,*)jx(h,i,j,k),jy(h,i,j,k),jz(h,i,j,k)
c          write(0,*)jx_lrf(h,i,j,k),jy_lrf(h,i,j,k),jz_lrf(h,i,j,k)
c          stop '***ERROR!*** Particle-flow non-vanishing in LRF!'
c          end if

c.. baryon number and total energy     
         n0_total(h)=n0_total(h)+n0(h,i,j,k)*vol
         p0_total(h)=p0_total(h)+p0(h,i,j,k)
         n0_lor(h)=n0_lor(h)+n0_lrf(h,i,j,k)*volce(h,i,j,k)
         p0_lor(h)=p0_lor(h)+p0_lrf(h,i,j,k)*gamma

c... FOR DEBUGGING ONLY .............................................!
c         write(72,*)'***********************************************'!
c         write(72,*)'CELL',h,i,j,k                                  !
c         write(72,*)'nopart(t,x,y,z)',nopart(h,i,j,k)               !
c         write(72,*)'vol',vol                                       !
c         write(72,*)'time',time                                     !
c         write(72,*)'n0(t,x,y,z)',n0(h,i,j,k)                       !
c         write(72,*)'p0(t,x,y,z)',p0(h,i,j,k)                       !
c         write(72,*)'px(t,x,y,z)',px(h,i,j,k)                       !
c         write(72,*)'py(t,x,y,z)',py(h,i,j,k)                       !
c         write(72,*)'pz(t,x,y,z)',pz(h,i,j,k)                       !
c         write(72,*)'betax',betax                                   !
c         write(72,*)'betay',betay                                   !
c         write(72,*)'betaz',betaz                                   !
c         write(72,*)'BETA',beta                                     !
c         write(72,*)'GAMMA',gamma                                   !
c         write(72,*)'TAU',tau(h,i,j,k)                              !
c         write(72,*)'LRF:p0,n0,jx,jy,jz'                            !
c         write(72,*)p0_lrf(h,i,j,k),n0_lrf(h,i,j,k),                !
c     &             jx_lrf(h,i,j,k),jy_lrf(h,i,j,k),                 !
c     &             jz_lrf(h,i,j,k)                                  !
c         write(72,*)'LRF:px,py,pz'                                  !
c         write(72,*)px_lrf,py_lrf,pz_lrf                            !
c         write(72,*)'LAB T_mu_mu:',taa(h,i,j,k),tbb(h,i,j,k),       !
c     &                             tcc(h,i,j,k),tdd(h,i,j,k)        !
c         write(72,*)'LAB T_0_1/2/3:',tab(h,i,j,k),                  !
c     &                               tac(h,i,j,k),tad(h,i,j,k)      !
c         write(72,*)'LAB T_1_2/3:',tbc(h,i,j,k),tbd(h,i,j,k)        ! 
c         write(72,*)'LAB T_2_3:',tcd(h,i,j,k)                       ! 
c         write(72,*)'LRF T_mu_mu:',taa_tr,tbb_tr,tcc_tr,tdd_tr      !
c         write(72,*)'LRF ANISO:',taa_tr2,tbb_tr2,tcc_tr2,tdd_tr2    !
c         write(72,*)'LRF T_0_1/2/3:',tab_tr,tac_tr,tad_tr           !
c         write(72,*)'LRF T_1_2/3:',tbc_tr,tbd_tr                    ! 
c         write(72,*)'LRF T_2_3:',tcd_tr                             ! 
c         write(72,*)'LRF cell volume:',volce(h,i,j,k)               !
c....................................................................!


 311     continue
 312    continue
 313   continue

c....................................................................!
c       write(0,*)'timestep',h !debug only                           !
c       write(0,*)'Cells with baryon content = ',cbar_time           !
c       write(0,*)'Cells with meson and without baryon content = ',  !
c     &           cmesnobar_time                                     !
c....................................................................!

 314  continue !loop over timesteps


      write(0,*)'completed'
      write(0,*)'LORENTZ-TRANSFORMATION finished'

c..............FOR.DEBUGGING.ONLY..................! 
c.. Check if BARYON NUMBER and ENERGY are conserved!
c      do 977 hh=1,timesteps                       !
c      write(0,*)'timestep =',hh*time              !
c      write(0,*)'total baryons =',n0_total(hh)    !
c      write(0,*)'total energy =',p0_total(hh)     !
c      write(0,*)'baryons after LTF =',n0_lor(hh)  !
c      write(0,*)'energy after LTF =',p0_lor(hh)   !
c 977  continue                                    !
c..................................................!
      return
      end

c****&|****************************************************************|X

