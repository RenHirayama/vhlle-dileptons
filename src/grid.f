c------------------------------------------------------------------------
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE GRID
c..
c.. This subroutine defines the cell quantities
c..
c.. CAUTION: All quantities are just summed up for all events ,they will
c..          be normalized later in TRANSFORM to number of events
c****&|****************************************************************|X

      SUBROUTINE grid(nexit,noe,timesteps,grd,grd_z,dx,dt,vol,
     &              OUTr0,OUTrx,OUTry,OUTrz,OUTp0,OUTpx,OUTpy,OUTpz,
     &              OUTmass,OUTityp,OUTiso3,OUTch,OUTlcoll,OUTcoll,
     &              OUTpaproc,OUTorigin,nopart,norho0,mrho0,p0rho0,
     &              pxrho0,pyrho0,pzrho0,nodelta,mdelta,p0delta,
     &              pxdelta,pydelta,pzdelta,noomega,momega,p0omega,
     &              pxomega,pyomega,pzomega,nophi,mphi,p0phi,
     &              pxphi,pyphi,pzphi,p0,px,py,pz,n0,npi,nk,jx,
     &              jy,jz,t00,t01,t02,t03,t11,t22,t33,t12,t13,t23,
     &              rhoeff)


      implicit none
      include 'defs.f'


c. general    
      real*8 dx,dt,vol,rapmax
      integer grd,grd_z

c.. event number
      integer noe

c. event information
      integer nexit(evtmax,timemax)
      integer timesteps
      
c.. OUTgoing particles
      real*8 OUTr0(timemax,outmax),OUTrx(timemax,outmax)
      real*8 OUTry(timemax,outmax),OUTrz(timemax,outmax) 
      real*8 OUTp0(timemax,outmax),OUTpx(timemax,outmax)
      real*8 OUTpy(timemax,outmax),OUTpz(timemax,outmax)
      real*8 OUTmass(timemax,outmax)
      integer OUTityp(timemax,outmax),OUTiso3(timemax,outmax)
      integer OUTch(timemax,outmax),OUTlcoll(timemax,outmax)
      integer OUTcoll(timemax,outmax),OUTpaproc(timemax,outmax)
      integer OUTorigin(timemax,outmax)

c. grid
c.. cell numbers in x/y/z-direction
      integer cx(timemax,outmax)
      integer cy(timemax,outmax)
      integer cz(timemax,outmax)

c.. cell position  
      real*8 xmin,xmax,ymin,ymax,zmin,zmax
      integer x,y,z

c.. cell properties
      real*8 nopart(timemax,grd,grd,grd_z)

      real*8 norho0(timemax,grd,grd,grd_z)
      real*8 mrho0(timemax,grd,grd,grd_z)
      real*8 p0rho0(timemax,grd,grd,grd_z) 
      real*8 pxrho0(timemax,grd,grd,grd_z)
      real*8 pyrho0(timemax,grd,grd,grd_z)
      real*8 pzrho0(timemax,grd,grd,grd_z)

      real*8 nodelta(timemax,grd,grd,grd_z)
      real*8 mdelta(timemax,grd,grd,grd_z)
      real*8 p0delta(timemax,grd,grd,grd_z) 
      real*8 pxdelta(timemax,grd,grd,grd_z)
      real*8 pydelta(timemax,grd,grd,grd_z)
      real*8 pzdelta(timemax,grd,grd,grd_z)

      real*8 noomega(timemax,grd,grd,grd_z)
      real*8 momega(timemax,grd,grd,grd_z)
      real*8 p0omega(timemax,grd,grd,grd_z) 
      real*8 pxomega(timemax,grd,grd,grd_z)
      real*8 pyomega(timemax,grd,grd,grd_z)
      real*8 pzomega(timemax,grd,grd,grd_z)

      real*8 nophi(timemax,grd,grd,grd_z)
      real*8 mphi(timemax,grd,grd,grd_z)
      real*8 p0phi(timemax,grd,grd,grd_z) 
      real*8 pxphi(timemax,grd,grd,grd_z)
      real*8 pyphi(timemax,grd,grd,grd_z)
      real*8 pzphi(timemax,grd,grd,grd_z)

      real*8 p0(timemax,grd,grd,grd_z) 
      real*8 px(timemax,grd,grd,grd_z)
      real*8 py(timemax,grd,grd,grd_z)
      real*8 pz(timemax,grd,grd,grd_z)
      real*8 n0(timemax,grd,grd,grd_z) 
      real*8 npi(timemax,grd,grd,grd_z) 
      real*8 nk(timemax,grd,grd,grd_z) 
      real*8 jx(timemax,grd,grd,grd_z)
      real*8 jy(timemax,grd,grd,grd_z)
      real*8 jz(timemax,grd,grd,grd_z)
      real*8 rhoeff(timemax,grd,grd,grd_z)

      real*8 t00(timemax,grd,grd,grd_z)
      real*8 t01(timemax,grd,grd,grd_z)
      real*8 t02(timemax,grd,grd,grd_z)
      real*8 t03(timemax,grd,grd,grd_z)
      real*8 t11(timemax,grd,grd,grd_z)
      real*8 t22(timemax,grd,grd,grd_z)
      real*8 t33(timemax,grd,grd,grd_z)
      real*8 t12(timemax,grd,grd,grd_z)
      real*8 t13(timemax,grd,grd,grd_z)
      real*8 t23(timemax,grd,grd,grd_z)
   
      real*8 dn0,djx,djy,djz,dp0,dpx,dpy,dpz
   
      real*8 rapidity
 
c.. other
      integer e,f,g,h,i,j,k,l,m,o
    
c. start analysis
c.. define cell properties / positions of cells

c.. loop over all timesteps & particles
c      write(0,*)'noe',noe !debug only
      do 212 i=1,timesteps
c       write(0,*)'timestep',i,timesteps !debug only
       do 211 j=1,nexit(noe,i)
c        write(0,*)'nexit',j,nexit(noe,i) !debug only

        if(OUTityp(i,j).eq.0) cycle     
 
c.. loop over lattice cell coordinates

        cx(i,j)=0
        cy(i,j)=0
        cz(i,j)=0

c.. x-coordinate
        do 221 k=1,grd
                      
         xmin=(-0.5d0*dx*grd)+((k-1)*dx)
         xmax=(-0.5d0*dx*grd)+(k*dx)

            if (OUTrx(i,j).gt.xmin.AND.OUTrx(i,j).le.xmax) then
             cx(i,j)=k
            end if
 221    continue 
    
c.. y-coordinate
        do 222 l=1,grd

         ymin=(-0.5d0*dx*grd)+((l-1)*dx)
         ymax=(-0.5d0*dx*grd)+(l*dx)

            if (OUTry(i,j).gt.ymin.AND.OUTry(i,j).le.ymax) then
             cy(i,j)=l
            end if
 222    continue
   
c.. z-coordinate
        do 223 m=1,grd_z

         zmin=(-0.5d0*dx*grd_z)+((m-1)*dx)
         zmax=(-0.5d0*dx*grd_z)+(m*dx)

            if (OUTrz(i,j).gt.zmin.AND.OUTrz(i,j).le.zmax) then
             cz(i,j)=m
            end if
 223    continue

        x=cx(i,j)
        y=cy(i,j)
        z=cz(i,j)
 
c.. check if all particles are in the grid
        if((x.eq.0).OR.(y.eq.0).OR.(z.eq.0)) then
c         write(0,*)'out of grid:',OUTrx(i,j),OUTry(i,j),OUTrz(i,j)
c         stop 'grid size too small'
          cycle
        end if
 
 
c        if((OUTityp(i,j).lt.100d0.AND.OUTityp(i,j).gt.0d0).OR.
c     &      OUTityp(i,j).gt.-100d0.AND.OUTityp(i,j).lt.0d0) then

c        if(OUTityp(i,j).gt.100d0.OR.OUTityp(i,j).lt.-100d0) then

        nopart(i,x,y,z)=nopart(i,x,y,z)+1d0

        dp0=OUTp0(i,j)
        p0(i,x,y,z)=p0(i,x,y,z)+dp0
        if(.NOT.(abs(p0(i,x,y,z)).gt.0.0d0)) then
         write(0,*)'WRONG P0:',i,x,y,z,p0(i,x,y,z),dp0
         stop
        endif

        dpx=OUTpx(i,j)
        px(i,x,y,z)=px(i,x,y,z)+dpx

        dpy=OUTpy(i,j)
        py(i,x,y,z)=py(i,x,y,z)+dpy

        dpz=OUTpz(i,j)
        pz(i,x,y,z)=pz(i,x,y,z)+dpz

        t00(i,x,y,z)=t00(i,x,y,z)+OUTp0(i,j)/vol
        t01(i,x,y,z)=t01(i,x,y,z)+OUTpx(i,j)/vol
        t02(i,x,y,z)=t02(i,x,y,z)+OUTpy(i,j)/vol
        t03(i,x,y,z)=t03(i,x,y,z)+OUTpz(i,j)/vol

        t11(i,x,y,z)=t11(i,x,y,z)+((OUTpx(i,j)**2)/OUTp0(i,j))/vol
        t22(i,x,y,z)=t22(i,x,y,z)+((OUTpy(i,j)**2)/OUTp0(i,j))/vol
        t33(i,x,y,z)=t33(i,x,y,z)+((OUTpz(i,j)**2)/OUTp0(i,j))/vol

        t12(i,x,y,z)=t12(i,x,y,z)+((OUTpx(i,j)*OUTpy(i,j))/OUTp0(i,j))/vol
        t13(i,x,y,z)=t13(i,x,y,z)+((OUTpx(i,j)*OUTpz(i,j))/OUTp0(i,j))/vol
        t23(i,x,y,z)=t23(i,x,y,z)+((OUTpy(i,j)*OUTpz(i,j))/OUTp0(i,j))/vol

c        endif
   
c........for debugging only.........................................!     
        if(i.eq.10.AND.x.eq.8.AND.y.eq.8.AND.z.eq.23) then         !
        write(0,*)'noe',noe                                        !
        write(0,*),i,x,y,z                                         !
        write(0,*),OUTrx(i,j),OUTry(i,j),OUTrz(i,j)                !
        write(0,*),p0(i,x,y,z),px(i,x,y,z),py(i,x,y,z),pz(i,x,y,z) !
        write(0,*)'vol',vol                                        !
        end if                                                     !
c...................................................................!

        dn0=1d0/vol
c        write(0,*),dn0 !debug only
        djx=OUTpx(i,j)/OUTp0(i,j)/vol
c        write(0,*),djx !debug only
        djy=OUTpy(i,j)/OUTp0(i,j)/vol
c        write(0,*),djy !debug only
        djz=OUTpz(i,j)/OUTp0(i,j)/vol
c        write(0,*),djz !debug only
        
c.. define baryon density [n(x_mu)=(1/Cell-Volume)*(Sum over baryons(1))]
c.. and baryon current [j(x_mu)=(1/Cell-Volume)*(Sum over baryons(p_mu/p_0)]
c... BARYONS
        if(OUTityp(i,j).lt.100d0.AND.OUTityp(i,j).gt.0d0) then
        n0(i,x,y,z)=n0(i,x,y,z)+dn0
c          if(i.eq.10.AND.x.eq.8.AND.y.eq.8.AND.z.eq.23) then 
c          write(0,*),'dn added to n0',i,x,y,z,dn0,n0(i,x,y,z) !debug only
c          end if 
        jx(i,x,y,z)=jx(i,x,y,z)+djx
        jy(i,x,y,z)=jy(i,x,y,z)+djy
        jz(i,x,y,z)=jz(i,x,y,z)+djz
        end if
        if(OUTityp(i,j).eq.1) rhoeff(i,x,y,z)=rhoeff(i,x,y,z)+dn0
        if(OUTityp(i,j).lt.100d0.AND.OUTityp(i,j).ge.2) then
         rhoeff(i,x,y,z)=rhoeff(i,x,y,z)+(0.5d0*dn0)
        endif
c... ANTI-BARYONS
        if(OUTityp(i,j).gt.-100d0.AND.OUTityp(i,j).lt.0d0) then
        n0(i,x,y,z)=n0(i,x,y,z)-dn0
c        write(0,*),'dn substracted from n0',i,x,y,z,dn0,n0(i,x,y,z) !debug only
        jx(i,x,y,z)=jx(i,x,y,z)-djx   
        jy(i,x,y,z)=jy(i,x,y,z)-djy
        jz(i,x,y,z)=jz(i,x,y,z)-djz
        end if
        if(OUTityp(i,j).eq.-1) rhoeff(i,x,y,z)=rhoeff(i,x,y,z)+dn0
        if(OUTityp(i,j).gt.-100d0.AND.OUTityp(i,j).le.-2) then
         rhoeff(i,x,y,z)=rhoeff(i,x,y,z)+(0.5d0*dn0)
        endif

c.. Define kaon and pion density [n(x_pi/k)=(1/Cell-Volume)*(Sum over pions/kaons))]
c... PIONS
        if(OUTityp(i,j).eq.101) then
        npi(i,x,y,z)=npi(i,x,y,z)+dn0
c       write(0,*),'dn added to npi',i,x,y,z,dn0,npi(i,x,y,z) !debug only
c... Also include resonances, that decay into pions according to their branching ratios
c... and the number of pions they decay into:
c...Mesons
c        elseif(OUTityp(i,j).eq.104) then
c        npi(i,x,y,z)=npi(i,x,y,z)+2.0d0*dn0
c        elseif(OUTityp(i,j).eq.103) then
c        npi(i,x,y,z)=npi(i,x,y,z)+2.6d0*dn0   
c....N*s
c        elseif(abs(OUTityp(i,j)).eq.2) then
c        npi(i,x,y,z)=npi(i,x,y,z)+1.45d0*dn0          
c        elseif(abs(OUTityp(i,j)).eq.3) then
c        npi(i,x,y,z)=npi(i,x,y,z)+1.45d0*dn0    
c        elseif(abs(OUTityp(i,j)).eq.4) then
c        npi(i,x,y,z)=npi(i,x,y,z)+0.75d0*dn0   
c        elseif(abs(OUTityp(i,j)).eq.5) then
c        npi(i,x,y,z)=npi(i,x,y,z)+0.8d0*dn0   
c        elseif(abs(OUTityp(i,j)).ge.6.AND.abs(OUTityp(i,j)).le.16) then
c        npi(i,x,y,z)=npi(i,x,y,z)+dn0 
c....Delta / Delta*s
c        elseif(abs(OUTityp(i,j)).eq.17) then
c        npi(i,x,y,z)=npi(i,x,y,z)+dn0   
c        elseif(abs(OUTityp(i,j)).eq.18) then
c        npi(i,x,y,z)=npi(i,x,y,z)+1.7d0*dn0   
c        elseif(abs(OUTityp(i,j)).eq.19) then
c        npi(i,x,y,z)=npi(i,x,y,z)+1.85d0*dn0   
c        elseif(abs(OUTityp(i,j)).eq.20) then
c        npi(i,x,y,z)=npi(i,x,y,z)+1.8d0*dn0   
c        elseif(abs(OUTityp(i,j)).ge.21.AND.abs(OUTityp(i,j)).le.26) then
c        npi(i,x,y,z)=npi(i,x,y,z)+1.5d0*dn0  
c........
        endif
    
c... KAONS
c        if(OUTityp(i,j).eq.106.OR.OUTityp(i,j).eq.-106) then
        if(OUTityp(i,j).eq.106) then
        nk(i,x,y,z)=nk(i,x,y,z)+dn0
c       write(0,*),'dn added to nk',i,x,y,z,dn0,nk(i,x,y,z) !debug only
        endif

c.. Finally count the number of rho_0 / omega / phi / Delta per cell, for if the temperature in 
c.. the cell will be too low for thermal emission a shining procedure will be
c.. applied
        if(OUTityp(i,j).eq.104.AND.OUTch(i,j).eq.0) then
         norho0(i,x,y,z)=norho0(i,x,y,z)+1.0d0
         mrho0(i,x,y,z)=OUTmass(i,j)
         p0rho0(i,x,y,z)=OUTp0(i,j)
         pxrho0(i,x,y,z)=OUTpx(i,j)
         pyrho0(i,x,y,z)=OUTpy(i,j)
         pzrho0(i,x,y,z)=OUTpz(i,j)
        end if
        if((OUTityp(i,j).eq.17.AND.OUTch(i,j).eq.0).OR.
     &     (OUTityp(i,j).eq.17.AND.OUTch(i,j).eq.1)) then
         nodelta(i,x,y,z)=nodelta(i,x,y,z)+1.0d0
         mdelta(i,x,y,z)=OUTmass(i,j)
         p0delta(i,x,y,z)=OUTp0(i,j)
         pxdelta(i,x,y,z)=OUTpx(i,j)
         pydelta(i,x,y,z)=OUTpy(i,j)
         pzdelta(i,x,y,z)=OUTpz(i,j)
        end if        
        if(OUTityp(i,j).eq.103.AND.OUTch(i,j).eq.0) then
         noomega(i,x,y,z)=noomega(i,x,y,z)+1.0d0
         momega(i,x,y,z)=OUTmass(i,j)
         p0omega(i,x,y,z)=OUTp0(i,j)
         pxomega(i,x,y,z)=OUTpx(i,j)
         pyomega(i,x,y,z)=OUTpy(i,j)
         pzomega(i,x,y,z)=OUTpz(i,j)
        end if
        if(OUTityp(i,j).eq.109.AND.OUTch(i,j).eq.0) then
         nophi(i,x,y,z)=nophi(i,x,y,z)+1.0d0
         mphi(i,x,y,z)=OUTmass(i,j)
         p0phi(i,x,y,z)=OUTp0(i,j)
         pxphi(i,x,y,z)=OUTpx(i,j)
         pyphi(i,x,y,z)=OUTpy(i,j)
         pzphi(i,x,y,z)=OUTpz(i,j)
        end if



 211   continue ! loop over particles
 212  continue ! loop over timesteps
     
      return
      end

c****&|****************************************************************|X

