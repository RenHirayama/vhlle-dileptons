c************************************************************8
c.. Authors: Marcus Bleicher, modified by Stephan Endres
c.. file 15 analysis...
c************************************************************8
      subroutine freezeout(grd,grd_z,dx,dt,vol,timemax,fo_cell,
     &                     collnoel,collnoin,elcount,inelcount,
     &                     noefo)

      implicit none
      
c. array sizes
      integer wmax,inmax,outmax
      parameter (wmax=20000)
      parameter (inmax=10)
      parameter (outmax=100)

c. general
c. interaction number ww (1..wmax)
      integer ww
c.. number of events
      integer noe,noefo

c. interaction header, one for every interaction ww
      integer nin(wmax),nexit(wmax),iline(wmax),ctag(wmax)
      real*8 acttime(wmax),sqrts(wmax),stot(wmax)
      real*8 sigpart(wmax),cdens(wmax)

c.. read INgoing particles
      integer INind(wmax,inmax)
      real*8 INr0(wmax,inmax),INrx(wmax,inmax)
      real*8 INry(wmax,inmax),INrz(wmax,inmax) 
      real*8 INp0(wmax,inmax),INpx(wmax,inmax)
      real*8 INpy(wmax,inmax),INpz(wmax,inmax)
      real*8 INmass(wmax,inmax)
      integer INityp(wmax,inmax),INiso3(wmax,inmax)
      integer INch(wmax,inmax),INlcoll(wmax,inmax)
      integer INcoll(wmax,inmax),INistr(wmax,inmax)
      integer INorigin(wmax,inmax)

c.. read OUTgoing particles
      integer OUTind(wmax,outmax)
      real*8 OUTr0(wmax,outmax),OUTrx(wmax,outmax)
      real*8 OUTry(wmax,outmax),OUTrz(wmax,outmax) 
      real*8 OUTp0(wmax,outmax),OUTpx(wmax,outmax)
      real*8 OUTpy(wmax,outmax),OUTpz(wmax,outmax)
      real*8 OUTmass(wmax,outmax)
      integer OUTityp(wmax,outmax),OUTiso3(wmax,outmax)
      integer OUTch(wmax,outmax),OUTlcoll(wmax,outmax)
      integer OUTcoll(wmax,outmax),OUTistr(wmax,outmax)
      integer OUTorigin(wmax,outmax)

c.. coarse-graining
      integer grd,grd_z,timemax
      real*8 dt,dx,vol
      integer fo_cell(timemax,grd,grd,grd_z)
      integer collnoel(timemax)  
      integer collnoin(timemax)  
      integer elcount(timemax,grd,grd,grd_z)  
      integer inelcount(timemax,grd,grd,grd_z)  

c.. input file
      character*96 file15,file73,file74
      
c.. other
      integer i,k,l,m
      
c................................................................!
      
      call getenv('ftn15',file15)
      if (file15(1:4).ne.'    ') then
        OPEN(UNIT=15,FILE=file15,STATUS='old',FORM='FORMATTED')
      else 
       stop 'File15 not defined'
      endif   

      call getenv('ftn73',file73)
      if (file73(1:4).ne.'    ') then
        OPEN(UNIT=73,FILE=file73,STATUS='replace',FORM='FORMATTED')
      else 
       stop 'File73 not defined'
      endif   
      
      call getenv('ftn74',file74)
      if (file74(1:4).ne.'    ') then
        OPEN(UNIT=74,FILE=file74,STATUS='replace',FORM='FORMATTED')
      else 
       stop 'File74 not defined'
      endif   

c.. set all to zero
      ww=0
      noe=0

      call setzero(nin,nexit,iline,ctag,
     & acttime,sqrts,stot,
     & sigpart,cdens,
     & INind,
     & INr0,INrx,INry,INrz, 
     & INp0,INpx,INpy,INpz,
     & INmass,INityp,INiso3,INch,INlcoll,INcoll,INistr,
     & INorigin,
     & OUTind,
     & OUTr0,OUTrx,OUTry,OUTrz, 
     & OUTp0,OUTpx,OUTpy,OUTpz,
     & OUTmass,OUTityp,OUTiso3,OUTch,OUTlcoll,OUTcoll,OUTistr,
     & OUTorigin)

c ---- start loop ----

 1    continue
      ww=ww+1
      if (ww.ge.wmax) stop 'increase wmax!'

c. read collision/event header
      read(15,502,end=99) nin(ww),nexit(ww),iline(ww),ctag(ww),
     &  acttime(ww),sqrts(ww),stot(ww),sigpart(ww),cdens(ww)

c..      the event header (nin(ww)=0)
      if (nin(ww).eq.-1) then
       ww=ww-1
       if (noe.gt.0) then

c... do what ever should be done
       call ana(ww,noe,nin,nexit,iline,ctag,
     & acttime,sqrts,stot,
     & sigpart,cdens,
     & INind,
     & INr0,INrx,INry,INrz, 
     & INp0,INpx,INpy,INpz,
     & INmass,INityp,INiso3,INch,INlcoll,INcoll,INistr,
     & INorigin,
     & OUTind,
     & OUTr0,OUTrx,OUTry,OUTrz, 
     & OUTp0,OUTpx,OUTpy,OUTpz,
     & OUTmass,OUTityp,OUTiso3,OUTch,OUTlcoll,OUTcoll,OUTistr,
     & OUTorigin,grd,grd_z,dx,dt,vol,timemax,fo_cell,
     & collnoel,collnoin,elcount,inelcount)

       endif
c.. a new event starts. set all to zero and noe+1
       noe=noe+1
c       write(0,*)'event: ',noe
       ww=0
c... not necessary to set things to zero again
c      call setzero(nin,nexit,iline,ctag,
c     & acttime,sqrts,stot,
c     & sigpart,cdens,
c     & INind,
c     & INr0,INrx,INry,INrz, 
c     & INp0,INpx,INpy,INpz,
c     & INmass,INityp,INiso3,INch,INlcoll,INcoll,INistr,
c     & INorigin,
c     & OUTind,
c     & OUTr0,OUTrx,OUTry,OUTrz, 
c     & OUTp0,OUTpx,OUTpy,OUTpz,
c     & OUTmass,OUTityp,OUTiso3,OUTch,OUTlcoll,OUTcoll,OUTistr,
c     & OUTorigin)

c.. start all over
       goto 1
      endif
      
c. store all the ingoing particle of this interaction (i=1..nin(ww))
      do 11 i=1,nin(ww)
       if (i.ge.inmax) stop 'increase inmax!'
c.. get particle information
      read(15,501) INind(ww,i),INr0(ww,i),INrx(ww,i),
     &    INry(ww,i),INrz(ww,i), 
     &    INp0(ww,i),INpx(ww,i),INpy(ww,i),INpz(ww,i),
     &    INmass(ww,i),INityp(ww,i),INiso3(ww,i),INch(ww,i),
     &    INlcoll(ww,i),INcoll(ww,i),INistr(ww,i),INorigin(ww,i)
 11   continue
         
c. read outgoing particles
      do 20 i=1,nexit(ww)
       if (i.ge.outmax) stop 'increase outmax!'
      read(15,501) OUTind(ww,i),OUTr0(ww,i),OUTrx(ww,i),
     &    OUTry(ww,i),OUTrz(ww,i), 
     &    OUTp0(ww,i),OUTpx(ww,i),OUTpy(ww,i),OUTpz(ww,i),
     &    OUTmass(ww,i),OUTityp(ww,i),OUTiso3(ww,i),OUTch(ww,i),
     &    OUTlcoll(ww,i),OUTcoll(ww,i),OUTistr(ww,i),OUTorigin(ww,i)
 20   continue

c. start all over with next event    
      goto 1

c.. all read in
 99   continue
      noefo=noe
      write(0,*)'Read in of file 15 finished'
      return
c....
c formats: file15
c 501  format(i5,9e16.8,i7,2i3,i6,i5,i3,i15)
c header-line for each collision in file15
c 502  format(i1,i8,i4,i7,f8.3,4e12.4)
c same with index for file15
 501  format(i5,9e16.8,i11,2i3,i9,i5,i3,i15)
c header-line for each collision in file15
c 502  format(i1,i8,i4,i7,f8.3,4e12.4)
 502  format(i8,i8,i4,i7,f8.3,4e12.4)
      end

c--------------------------------------------------------------------

c... set all arrays to zero...

      subroutine setzero(nin,nexit,iline,ctag,
     & acttime,sqrts,stot,
     & sigpart,cdens,
     & INind,
     & INr0,INrx,INry,INrz, 
     & INp0,INpx,INpy,INpz,
     & INmass,INityp,INiso3,INch,INlcoll,INcoll,INistr,
     & INorigin,
     & OUTind,
     & OUTr0,OUTrx,OUTry,OUTrz, 
     & OUTp0,OUTpx,OUTpy,OUTpz,
     & OUTmass,OUTityp,OUTiso3,OUTch,OUTlcoll,OUTcoll,OUTistr,
     & OUTorigin)

      implicit none
c. array sizes
      integer wmax,inmax,outmax
      parameter (wmax=20000)
      parameter (inmax=10)
      parameter (outmax=50)

c. general
c. interaction number ww (1..wmax)
      integer ww

c. interaction header, one for every interaction ww
      integer nin(wmax),nexit(wmax),iline(wmax),ctag(wmax)
      real*8 acttime(wmax),sqrts(wmax),stot(wmax)
      real*8 sigpart(wmax),cdens(wmax)

c.. read INgoing particles
      integer INind(wmax,inmax)
      real*8 INr0(wmax,inmax),INrx(wmax,inmax)
      real*8 INry(wmax,inmax),INrz(wmax,inmax) 
      real*8 INp0(wmax,inmax),INpx(wmax,inmax)
      real*8 INpy(wmax,inmax),INpz(wmax,inmax)
      real*8 INmass(wmax,inmax)
      integer INityp(wmax,inmax),INiso3(wmax,inmax)
      integer INch(wmax,inmax),INlcoll(wmax,inmax)
      integer INcoll(wmax,inmax),INistr(wmax,inmax)
      integer INorigin(wmax,inmax)

c.. read OUTgoing particles
      integer OUTind(wmax,outmax)
      real*8 OUTr0(wmax,outmax),OUTrx(wmax,outmax)
      real*8 OUTry(wmax,outmax),OUTrz(wmax,outmax) 
      real*8 OUTp0(wmax,outmax),OUTpx(wmax,outmax)
      real*8 OUTpy(wmax,outmax),OUTpz(wmax,outmax)
      real*8 OUTmass(wmax,outmax)
      integer OUTityp(wmax,outmax),OUTiso3(wmax,outmax)
      integer OUTch(wmax,outmax),OUTlcoll(wmax,outmax)
      integer OUTcoll(wmax,outmax),OUTistr(wmax,outmax)
      integer OUTorigin(wmax,outmax)

c.. other
      integer i

      do 1 ww=1,wmax
       nin(ww)=0
       nexit(ww)=0
       iline(ww)=0
       ctag(ww)=0
       acttime(ww)=0
       sqrts(ww)=0
       stot(ww)=0
       sigpart(ww)=0
       cdens(ww)=0

       do 2 i=1,inmax
        INind(ww,i)=0
        INr0(ww,i)=0
        INrx(ww,i)=0
        INry(ww,i)=0
        INrz(ww,i)=0 
        INp0(ww,i)=0
        INpx(ww,i)=0
        INpy(ww,i)=0
        INpz(ww,i)=0
        INmass(ww,i)=0
        INityp(ww,i)=0
        INiso3(ww,i)=0
        INch(ww,i)=0
        INlcoll(ww,i)=0
        INcoll(ww,i)=0
        INistr(ww,i)=0
        INorigin(ww,i)=0
 2     continue
       
       do 3 i=1,outmax
        OUTind(ww,i)=0
        OUTr0(ww,i)=0
        OUTrx(ww,i)=0
        OUTry(ww,i)=0
        OUTrz(ww,i)=0 
        OUTp0(ww,i)=0
        OUTpx(ww,i)=0
        OUTpy(ww,i)=0
        OUTpz(ww,i)=0
        OUTmass(ww,i)=0
        OUTityp(ww,i)=0
        OUTiso3(ww,i)=0
        OUTch(ww,i)=0
        OUTlcoll(ww,i)=0
        OUTcoll(ww,i)=0
        OUTistr(ww,i)=0
        OUTorigin(ww,i)=0
 3     continue
 1    continue
      return
      end
c---------------------------------------------------------------
c... analyze the event...

      subroutine ana(ww,noe,nin,nexit,iline,ctag,
     & acttime,sqrts,stot,
     & sigpart,cdens,
     & INind,
     & INr0,INrx,INry,INrz, 
     & INp0,INpx,INpy,INpz,
     & INmass,INityp,INiso3,INch,INlcoll,INcoll,INistr,
     & INorigin,
     & OUTind,
     & OUTr0,OUTrx,OUTry,OUTrz, 
     & OUTp0,OUTpx,OUTpy,OUTpz,
     & OUTmass,OUTityp,OUTiso3,OUTch,OUTlcoll,OUTcoll,OUTistr,
     & OUTorigin,grd,grd_z,dx,dt,vol,timemax,fo_cell,
     & collnoel,collnoin,elcount,inelcount)

      implicit none
c. array sizes
      integer wmax,inmax,outmax
      parameter (wmax=20000)
      parameter (inmax=10)
      parameter (outmax=100)

c. general
c. interaction number ww (1..wmax)
      integer ww
c.. event number
      integer noe

c. interaction header, one for every interaction ww
      integer nin(wmax),nexit(wmax),iline(wmax),ctag(wmax)
      real*8 acttime(wmax),sqrts(wmax),stot(wmax)
      real*8 sigpart(wmax),cdens(wmax)

c.. read INgoing particles
      integer INind(wmax,inmax)
      real*8 INr0(wmax,inmax),INrx(wmax,inmax)
      real*8 INry(wmax,inmax),INrz(wmax,inmax) 
      real*8 INp0(wmax,inmax),INpx(wmax,inmax)
      real*8 INpy(wmax,inmax),INpz(wmax,inmax)
      real*8 INmass(wmax,inmax)
      integer INityp(wmax,inmax),INiso3(wmax,inmax)
      integer INch(wmax,inmax),INlcoll(wmax,inmax)
      integer INcoll(wmax,inmax),INistr(wmax,inmax)
      integer INorigin(wmax,inmax)

c.. read OUTgoing particles
      integer OUTind(wmax,outmax)
      real*8 OUTr0(wmax,outmax),OUTrx(wmax,outmax)
      real*8 OUTry(wmax,outmax),OUTrz(wmax,outmax) 
      real*8 OUTp0(wmax,outmax),OUTpx(wmax,outmax)
      real*8 OUTpy(wmax,outmax),OUTpz(wmax,outmax)
      real*8 OUTmass(wmax,outmax)
      integer OUTityp(wmax,outmax),OUTiso3(wmax,outmax)
      integer OUTch(wmax,outmax),OUTlcoll(wmax,outmax)
      integer OUTcoll(wmax,outmax),OUTistr(wmax,outmax)
      integer OUTorigin(wmax,outmax)

c.. other
      integer j,i,ii,iii,iiii,iww,iww2,iwwf
      integer i1,i2,i3,i4,k,l,m
      integer iproc, iproc2,foproc
      integer IITYP,iclass
      integer ix,check1,check2,obs,check

      real*8 xmin,xmax,ymin,ymax,zmin,zmax
      integer h,x,y,z
      integer cx(wmax),cy(wmax),cz(wmax)
   
c.. coarse-graining
      integer grd,grd_z,timemax
      real*8 dt,dx,vol
      integer fo_cell(timemax,grd,grd,grd_z)  
      integer marker(timemax,grd,grd,grd_z) 
      integer elcount(timemax,grd,grd,grd_z)  
      integer inelcount(timemax,grd,grd,grd_z)  
      integer collnoel(timemax)  
      integer collnoin(timemax)  

c.. start analysis -------------BELOW--------------
c.. some examples (one complete event is here and can be analysed)
c      write(0,*)'Event #, WW #',noe,ww !Debug only
c      write(0,*)'grd,grd_z,dx,dt,vol,timemax',grd,grd_z,dx,dt,vol,timemax !Debug only
      marker(1:timemax,1:grd,1:grd,1:grd_z)=0

c.. loop over all interactions ww
      do  iww=1,ww
c       write(0,*)'iww =',iww ! Debug only
c.. loop over ingoing particles
	if (nin(iww).eq.1) then ! resonance decay
 	 iproc=1
	elseif (nin(iww).eq.2 .and. nexit(iww).eq.1) then ! resonance formation
	 iproc=2
	elseif (nin(iww).eq.2 .and. nexit(iww).eq.2) then ! 2->2 collision
	 iproc=3
	elseif (nin(iww).eq.2 .and. nexit(iww).gt.2) then !2->n inelastic coll.
	 iproc=4
	elseif (nexit(iww).eq.0) then ! Paul blocked
	 iproc=5
	else 
	 write (0,*) 'undefined iproc'
	 stop
	endif

c        write(0,*)'iproc =',iproc

        foproc=0
	iproc2=0
	ix=1
	if (iproc.eq.2) then
c.. test if resonance decays in same particles as initial particles?
	  iproc2=1 ! inelastic, check if it was elastic,i.e. iproc2->0
          foproc=1
	  do iww2=iww+1,ww
	   if (nin(iww2).eq.1) then ! only decays are interesting
            if (INind(iww2,1).eq.OUTind(iww,1))then ! same resonance?
	     if (nexit(iww2).eq.2) then ! only 2 particle decays 
	      if ( (INityp(iww,1).eq.OUTityp(iww2,1)) 
     &	       .and. (INityp(iww,2).eq.OUTityp(iww2,2)) ) then
	        iproc2=0 ! elastic
                foproc=0
	        ix=iww2
	      endif
	      if ( (INityp(iww,1).eq.OUTityp(iww2,2)) 
     &	       .and. (INityp(iww,2).eq.OUTityp(iww2,1)) ) then
	        iproc2=0 ! elastic
                foproc=0
	        ix=iww2
	      endif
             endif
	    endif
	   endif
	  enddo
	endif

      	if (iproc.eq.3) then ! elastic or not?
	 iproc2=1
         foproc=1
	 if ( (INityp(iww,1).eq.OUTityp(iww,1)) 
     &	  .and. (INityp(iww,2).eq.OUTityp(iww,2)) ) then
	   iproc2=0 ! elastic
           foproc=0
	 endif
	 if ( (INityp(iww,1).eq.OUTityp(iww,2)) 
     &	  .and. (INityp(iww,2).eq.OUTityp(iww,1)) ) then
	   iproc2=0 ! elastic
           foproc=0
	 endif
	end if

	do iwwf=iww+1,ww
	 if (nin(iwwf).eq.2) then ! only scattering processes are interesting
	  do j=1,nexit(iww)
	   if (INind(iwwf,1).eq.OUTind(iww,j).OR.INind(iwwf,2).eq.OUTind(iww,j))then ! same particle scattering again?
            if (.not.(iline(iwwf).eq.13.OR.iline(iwwf).eq.17.OR.iline(iwwf).eq.19
     &          .OR.iline(iwwf).eq.22.OR.iline(iwwf).eq.26.OR.iline(iwwf).eq.38))then
	     foproc=0
            endif
	   endif
	  enddo
	 endif
	enddo

c output
	if (iproc2.eq.0) then ! elastic
	 if (iproc.eq.2 .or. iproc.eq. 3) then
           write(73,*) "0 ",INr0(iww,1),INrx(iww,1),INry(iww,1),INrz(iww,1) 
           collnoin(int(acttime(iww)/dt))=collnoin(int(acttime(iww)/dt))+1
	 endif
	endif
	if (iproc2.eq.1) then ! inelastic
	 if (iproc.eq. 2 .or. iproc.eq. 3 .or.iproc.eq. 4) then
           write(73,*) "1 ",INr0(iww,1),INrx(iww,1),INry(iww,1),INrz(iww,1)
           collnoel(int(acttime(iww)/dt))=collnoel(int(acttime(iww)/dt))+1
	 endif
	endif

	
c.. IF PROCESS IS FREEZEOUT, MARK CORRESPONDING CELL AS FO-CELL
        if(foproc.eq.1.AND.(iproc.eq.2.OR.iproc.eq.3.or.iproc.eq.4))then
c        write(0,*)'FO conditions fulfilled' !Debug only
        
c.. loop over grid cell coordinates
        cx(iww)=0
        cy(iww)=0
        cz(iww)=0

c.. x-coordinate
        do 221 k=1,grd
                      
         xmin=(-0.5d0*dx*grd)+((k-1)*dx)
         xmax=(-0.5d0*dx*grd)+(k*dx)

            if (INrx(iww,1).gt.xmin.AND.INrx(iww,1).le.xmax) then
             cx(iww)=k
c             write(0,*)'cx =',k !Debug only
            end if
 221    continue 
    
c.. y-coordinate
        do 222 l=1,grd

         ymin=(-0.5d0*dx*grd)+((l-1)*dx)
         ymax=(-0.5d0*dx*grd)+(l*dx)

            if (INry(iww,1).gt.ymin.AND.INry(iww,1).le.ymax) then
             cy(iww)=l
c             write(0,*)'cy =',l !Debug only
            end if
 222    continue
   
c.. z-coordinate
        do 223 m=1,grd_z

         zmin=(-0.5d0*dx*grd_z)+((m-1)*dx)
         zmax=(-0.5d0*dx*grd_z)+(m*dx)

            if (INrz(iww,1).gt.zmin.AND.INrz(iww,1).le.zmax) then
             cz(iww)=m
c             write(0,*)'cz =',m !Debug only
            end if
 223    continue

        x=cx(iww)
        y=cy(iww)
        z=cz(iww)
        
c.. timestep
        
        h=int(acttime(iww)/dt)
c        write(0,*)'timestep =',h !Debug only

       if(x.gt.0.AND.y.gt.0.AND.z.gt.0)then
        if(h.le.timemax.AND.x.le.grd.AND.y.le.grd.AND.z.le.grd_z)then
c.. collision counter for cell       
        if(iproc2.eq.1) inelcount(h,x,y,z)=inelcount(h,x,y,z)+1
        if(iproc2.eq.0) elcount(h,x,y,z)=elcount(h,x,y,z)+1
c.. freeze-out cell or not?
         if(marker(h,x,y,z).eq.0) then
          fo_cell(h,x,y,z)=fo_cell(h,x,y,z)+1
          marker(h,x,y,z)=1
c          write(0,*)'freezeout: (h,x,y,z)=',h,x,y,z !Debug only
         endif
        endif
       endif

      end if ! iproc2.eq.1
      end do
      
      return
c formats: file15
 999  format(10e16.8,i7)
      end
