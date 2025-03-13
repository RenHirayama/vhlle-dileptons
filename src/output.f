c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE OUTPUT
c..
c.. Routine writes the output into a textfile
c*****|****************************************************************|X

      SUBROUTINE output(grd,grd_z,dx,dt,timesteps,tau,p0,n0,edens_lrf,
     &                  n0_lrf,jx,jy,jz,mub_cell,temp_cell,
     &                  gce)

      implicit none
      include 'defs.f'

c. variables 
      integer grd,grd_z,timesteps
      real*8 tau(timemax,grd,grd,grd_z)
      real*8 p0(timemax,grd,grd,grd_z)
      real*8 n0(timemax,grd,grd,grd_z)
      real*8 edens_lrf(timemax,grd,grd,grd_z)
      real*8 jx(timemax,grd,grd,grd_z)
      real*8 jy(timemax,grd,grd,grd_z)
      real*8 jz(timemax,grd,grd,grd_z)
      real*8 n0_lrf(timemax,grd,grd,grd_z)
      real*8 mub_cell(timemax,grd,grd,grd_z) 
      real*8 temp_cell(timemax,grd,grd,grd_z) 
      real*8 gce(timemax,grd,grd,grd_z)
      real*8 xmin,xmax,ymin,ymax,zmin,zmax,dx,dt 
      real*8 vxce,vyce,vzce

c. others
      integer h,i,j,k
      integer counter,counterx,separ
      real*8 percent,t,mub
      character dot

      parameter(separ=999)

c-----------------------------------------------------------------------X
      
      counter=0
      counterx=1  
      dot='.'    

c. write the output into textfile
      write(0,*)'Write CELL PROPERTIES from CG in standard output'         
c      write(0,*)'dt =',dt ! Debug only
      do 804 h=1,timesteps
c       call timestep(h,dt) ! Write timestep header in output
       write(0,'(a1)',advance='NO'),dot
c       write(0,'(a1)'),dot
       counter=counter+1
       if(counter.eq.10) then
        percent=100d0*counterx*counter/timesteps
        counterx=counterx+1
        counter=0
        write(0,'(f4.0)',advance='NO'),percent
c        write(0,'(f4.0)'),percent
        write(0,*),'%'
       end if

       do 803 i=1,grd
        do 802 j=1,grd
         do 801 k=1,grd_z

          xmin=(-0.5d0*dx*grd)+((i-1)*dx)
          xmax=(-0.5d0*dx*grd)+(i*dx)

          ymin=(-0.5d0*dx*grd)+((j-1)*dx)
          ymax=(-0.5d0*dx*grd)+(j*dx)

          zmin=(-0.5d0*dx*grd_z)+((k-1)*dx)
          zmax=(-0.5d0*dx*grd_z)+(k*dx)

          if(temp_cell(h,i,j,k).lt.0.050d0) cycle

          if(n0_lrf(h,i,j,k).ne.0.0d0) then

           vxce=jx(h,i,j,k)/n0(h,i,j,k)
           vyce=jy(h,i,j,k)/n0(h,i,j,k)
           vzce=jz(h,i,j,k)/n0(h,i,j,k)

C...OLD OUTPUT FORMAT
c           write(*,111)h,i,j,k,(xmax+xmin)/2.0d0,(ymax+ymin)/2.0d0,
c     &                (zmax+zmin)/2.0d0,vxce,vyce,vzce,
c     &                p0(h,i,j,k),n0(h,i,j,k),
c     &                edens_lrf(h,i,j,k),n0_lrf(h,i,j,k),
c     &                mub_cell(h,i,j,k),temp_cell(h,i,j,k),
c     &                gce(h,i,j,k)

C...NEW OUTPUT FORMAT optimized for charm calculations
           write(*,121)h*dt,(xmax+xmin)/2.0d0,(ymax+ymin)/2.0d0,
     &                (zmax+zmin)/2.0d0,vxce,vyce,vzce,
     &                mub_cell(h,i,j,k),temp_cell(h,i,j,k)

          end if 

 801     continue
 802    continue
 803   continue
 804  continue

c      write(*,'((16i5,2x))'),separ,separ,separ,separ,
c     &        separ,separ,separ,separ,separ,separ,separ 

      write(0,*)'completed'  
      write(0,*)'*****************************************'      

c------------------------------------------------------------------------
c. Output format definitions 

 111  format(4(i4,2x),6(e11.4,2x),6(e14.7,2x),1(e14.4,2x))
 121  format(7(e11.4,2x),2(e14.7,2x))   
 
      return                                               
      end  

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE TIMESTEP
c..
c*****|****************************************************************|X

      SUBROUTINE timestep(step,dt)
      implicit none

      integer separ,step
      real*8 dt,time
      parameter (separ=888)

c      write(0,*)'dt =',dt! Debug only

      time=step*dt

c..   write event separator
c        write(44,'(5i12)') separ, separ, separ, separ, separ
        write(*,'((i4,2x),(2e12.4,2x),(8i5,2x))'),separ,time,dt,step,
     &        step,step,step,step,step,step,step 

      end

c-----|----------------------------------------------------------------|X
c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE TIMESTEPDIL
c..
c*****|****************************************************************|X

      SUBROUTINE timestepdil(step,dt)
      implicit none

      integer separ,step
      real*8 dt,time
      parameter (separ=888)

c      write(0,*)'dt =',dt! Debug only

      time=step*dt

c..   write event separator
c        write(44,'(5i12)') separ, separ, separ, separ, separ
        write(71,'((i4,2x),(2e12.4,2x),(8i5,2x))'),separ,time,dt,step,
     &        step,step,step,step,step,step,step 

      end

c-----|----------------------------------------------------------------|X
