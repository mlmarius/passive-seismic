      program raytrace
* ANU: 29/03/2001
*           Alex Gorby
* P/S-wave and Bulk/Shear version for region/global inversions

      implicit real*8 (a-h, o-z)
      integer, parameter :: npcom=440000,nlayer=100,nlayeri=100
      integer, parameter :: nunkn=500000
      integer, parameter :: P=1, S=2, Pdiff=1, SKS=3, PP=4
      integer iraygeometry
c array of 2 elements each with 1 chars
      character*1 cmodel(2)
      parameter (msg = 256)
c     parameter for regional study, all values in colat and colon
      real scaleloc,scaleglob
      common /refmodel/ np15
      dimension  w(3,msg+1)
      real*8 xor(nunkn)
      dimension icellnum(npcom)
      character*54 head1,head2*154
      character*256 fpar,fvel,fpergl,fperlc,fdat,fsynt,frel,fray,fdel
      common /param/dx,dy,rmilo,rmalo,rmila,
     &rmala,dxi,dyi,numblock,nx,ny,nz,nxi,nyi,nzi,icellnum
      dimension r(npcom),v(npcom),d(npcom)
      dimension r1(npcom),v1(npcom),d1(npcom)
      dimension xlayer(0:nlayer),dz(nlayer),dzsum(nlayer)
      dimension xlayeri(0:nlayeri)
      common /value/ pi,con,rcon,deg,r,v,d,
     >r1,v1,d1,dz,dzsum
      common /model/ xlayer,xlayeri
      common /perturb/ xan(nunkn)
      integer (kind=8) :: JEV, NBLOCK, NST

C'source_block', 'station_block',
C                'residual', 'event_number',
C                'source_longitude', 'source_latitude',
C                'source_depth', 'station_longitude', 'station_latitude',
C                'observed_tt', 'locations2degrees', 'P_or_S'

      TYPE :: ARRIVAL
          INTEGER(KIND=8) :: SOURCE_BLOCK
          INTEGER(KIND=8) :: STATION_BLOCK
          REAL :: RESIDUAL
          INTEGER(KIND=8) :: EVENT_NUMBER
          REAL :: SOURCE_LONGITUDE
          REAL :: SOURCE_LATITUDE
          REAL :: SOURCE_DEPTH
          REAL :: STATION_LONGITUDE
          REAL :: STATTION_LATITUDE
          REAL :: OBSERVED_TT
          REAL :: DELTA
          INTEGER(KIND=1) :: P_OR_S
      END TYPE ARRIVAL

      TYPE(ARRIVAL) :: NEXT_ARRIVAL

*-------------------- FILES AND PARAMETER DEFINITIONS ------------------

c I have this local and global parmetrisations, example param2x2
      data fpar
     >/'/home/sudipta/repos/passive-seismic/raytracer/params/param2x2'/

c reference 1d velocity model for P and S
      data fvel
     > /'/home/sudipta/repos/passive-seismic/raytracer/params/ak135.15.SKS'/

c perturbation of 1D velocity model in global grid, example is perglob2x2.dat
      data fpergl
     >/'/home/sudipta/repos/passive-seismic/raytracer/params/perglob2x2.dat'/
c perturbation of 1D velocity model in local grid , example is perloc2x2.dat
      data fperlc 
     >/'/home/sudipta/repos/passive-seismic/raytracer/params/perloc2x2.dat'/
c This is P or S arrivals
      data fdat
     > /'/home/sudipta/repos/passive-seismic/raytracer/region_p.csv'/

c OUTPUTS
c relocation terms
      data frel /'Preloc.chk.it2'/
c traced rays
      data fray
     > /'Praylength.chk.it2'/
c time delays/residuals
      data fdel
     > /'Presiduals.chk.it2'/
c generated synthetic arrivals if requested by parameters below
      data fsynt /'P5-asia-restestchk.dat'/

*   ICHCKR  -> 1 - normal inversion, 2 - checkerboard generation,
*   results are writen into fsynt
      ichckr =1  

*    P or S model ?

      icmodel = P
      iraygeometry = P

* Scale parameter for perturbations, 0 is 1D model and 1 is 3D velocity model
      scaleloc =1.0d0
      scaleglob=0.0d0

*-------------------------------END---------------------------------------------
      pi=4.0*atan(1.0)
      deg=111.195d0
      con=pi/180.
      rcon=180./pi
      dpi = asin(1.)/ 90.
      r2d = 90./asin(1.)
      ieq=0
      rms=0.0e10
      nmdat=0
      cmodel(1)='P'
      cmodel(2)='S'


      if(ichckr.EQ.2)then
      write(*,*)'Checkerboard data generation!!!'
      open(16,file=fsynt)
      endif

      if(scaleglob.EQ.0.0)write(*,*)'Global model is ZERO'
      if(scaleloc.EQ.0.0)write(*,*)'Local model is ZERO'

* The data file:
      open(2,file= fdat, form='formatted')
*    To write calculated time delay
      open(3,file=fdel,form='formatted')
      if(icmodel.EQ.P)then
        write(3,'(a1)') cmodel(P)
      else
        write(3,'(a1)') cmodel(S)
      endif
*    To write  relocation coefficients
      open(11,file=frel)
* Read velocity model
      open(10,file=fvel)
* To save the calculated raypath data:
      open(14,file=fray)
* Read the layer division from an array outside:
      open(19,file=fpar)
* To provide anomalies:
* Local model
      open(20,file=fperlc)
* Global surrounding model
      open(21,file=fpergl)


* NOTE AK135 contains: radius,depth,Vp,Vs,density
* MODIFIED FOR AK135 & S-WAVES : (use nint & vs1 & vs !!!)
* NOTE: v1(i)=vs1 & v(i)=vs are a MUST for S !!!
 
      i=1
56    continue
      read(10,*,end=57)r(i),d(i),vp,vs,rho
      zzz=d(i)
      d(i)=float(nint(d(i)))
      r(i)=float(nint(r(i)))
      if(vs.LE.0.0000000000)vs=vp
      if(icmodel.EQ.P)then
        v(i)=vp
      else
        v(i)=vs
      endif
      i=i+1
      goto 56
57    np15=i-1
 
***************************************************
*     READ parameterization
      read(19,*)nx,ny,nz,dx,dy
* Read global layer division data:
      do i=0,nz
      read(19,*)xlayer(i)
      enddo
* Read local layer division data
      read(19,*)rmilo,rmalo,rmila,rmala,nxi,nyi,nzi
      do i=0,nzi
      read(19,*)xlayeri(i)
      enddo

 
      rmila=90.-rmila
      rmala=90.-rmala
      call minmax(rmila,rmala)
      call minmax(rmilo,rmalo)
      dxi=((rmalo-rmilo)/float(nxi))
      dyi=((rmala-rmila)/float(nyi))
      nmax1=nxi*nyi*nzi
      nmax2=nx*ny*nz+nmax1
      read(19,*)numcel,numblock 
      if(numcel.NE.nmax1)then
        write(*,*)'Wrong parameterization!'
        stop
      endif
* Read the target block number:
      read(19,*)(icellnum(i),i=1,nmax1)

***************************************************
* Velocity perturbations for each cell
* local model
      read(20,911)head1
      read(20,912)head2
      read(20,*)(xor(i),i=1,nmax1)
      xor(nmax1+1)=0.0

* global model
      read(21,911)head1
      read(21,912)head2
      read(21,*)(xor(i),i=nmax1+1,nmax2)
      xor(nmax2+1)=0.0
c777   format(10f)

      do j=1,nmax1
         xan(j)=xor(j)*0.01*scaleloc
      enddo
      do j=nmax1+1,nmax2
         xan(j)=xor(j)*0.01*scaleglob
      enddo

*     Begining for the ray tracing
      deltamin=1000.
      deltamax=0.
      depthmin=1000.
      depthmax=0.

10    read(2, *) NEXT_ARRIVAL
      nblock = NEXT_ARRIVAL%SOURCE_BLOCK
      nst = NEXT_ARRIVAL%STATION_BLOCK
      RES = NEXT_ARRIVAL%RESIDUAL
      JEV = NEXT_ARRIVAL%EVENT_NUMBER
      VLON = NEXT_ARRIVAL%SOURCE_LONGITUDE
      VLAT = NEXT_ARRIVAL%SOURCE_LATITUDE
      DEPTH = NEXT_ARRIVAL%SOURCE_DEPTH/1000.0
      SLON = NEXT_ARRIVAL%STATION_LONGITUDE
      SLAT = NEXT_ARRIVAL%STATTION_LATITUDE
      DELAY = NEXT_ARRIVAL%OBSERVED_TT
      DELTA = NEXT_ARRIVAL%DELTA


C QUESTION FOR ALEXEI: DELAY == OBSERVERD_TT?

!    nblock,nst,res,jev
!    >,vlon,vlat,
!    >depth,
!    >slon,slat,
!    >delay,delta
!     print *, NEXT_ARRIVAL
c     if(delta.GT.95.0)goto 10

      if(deltamin.GT.delta)deltamin=delta
      if(deltamax.LT.delta)deltamax=delta
      if(depthmin.GT.depth) depthmin=depth
      if(depthmax.LT.depth) depthmax=depth

150   format(2i8,f6.1,i8,
     >f9.3,f8.3,
     >f6.1,
     >f9.3,f8.3,
     >f9.3,f9.3,i4)


c  16092   17524   5.4  587486  125.588  78.790  23.1  134.351 -19.819 822.100  99.000

151   format(2i8,f6.1,i8,
     >f9.3,f8.3,
     >f6.1,
     >f9.3,f8.3,
     >2f8.3)

!     if(nmdat.LE.15)then
      write(*, *) 'ARRIVAL', NMDAT, ':', nblock, nst, res,
     > NEXT_ARRIVAL%EVENT_NUMBER, JEV, vlon, vlat,
     > depth, slon, slat, delay, delta, NEXT_ARRIVAL%P_or_S
!     endif


      nmdat=nmdat+1


c         write(*,*)vlon,vlat

        if(slon.GT.180.)then
           slon=slon-360.
        endif
        if(vlon.GT.180.)then
           vlon=vlon-360.
        endif

        hr=0.0e10 
        aar=slat
        bbr=slon
        hs=depth
        aas=vlat
        bbs=vlon
        dth2=0.0
        dth3=0.0
        dth4=0.0
c       write(*,*)'Event:Lat,Lon,Dp, Station: Lat, Lon, Dp'
c       write(*,*)aas, bbs, hs, aar, bbr, hr

        call pbr(aas, bbs, hs, aar, bbr, hr,w,ni,iraygeometry)
        ieq=ieq+1

c       write(*,*)w(2,1),w(3,1),w(1,1)
c       write(*,*)
c    &w(2,int((ni+1)/2)),w(3,int((ni+1)/2)),w(1,int((ni+1)/2))
c       write(*,*)w(2,ni+1),w(3,ni+1),w(1,ni+1)
      
c       write(*,*)'rp1= ',rp

        call cderiv( w, ni, dth2, dth3, dth4)
        dth1=1.0
        write(11, *) nblock, dth1, dth2, dth3, dth4

        if(ichckr.EQ.2)delay=0.0e10

        ichr=14
        ichd=3
        call rayl(ichr,ichd,w,ni,delay)


      if(ichckr.EQ.2)then

*     Rewind the output file in order to save the space on the disk
      rewind(3)
      rewind(11)
      rewind(14)

*     Write data on the disk

      delay=-1.*delay
      write(16,150)nblock,nst,res,jev
     >,vlon,vlat,
     >depth,
     >slon,slat,
     >delay,delta
      endif


        rms=rms+delay**2
      goto 10

200   continue
* ---------------------------------------------------------------

*New:
110   format(2i6,f7.2,i4,i6,f8.3,f8.3,
     >f6.1,
     >f8.4,
     >f8.3,f8.3,2f9.3,
     >2i2,f7.2,f7.2,f7.2,f9.3,f7.2,f7.2,i4,a8,f8.3,i5,i3)

120   format(i8,4f8.2)
555   format(i7,1x,f12.2)
666   format(f7.2)

911   format(a54)
912   format(a154)

      if((ieq-1).GT.0)then
      rms=sqrt(rms/(ieq-1))
      else
      rms=0.
      endif

      write(*,*)'Total used rays: ',ieq
      write(*,*)'Total RMS: ',rms
      write(*,*)'Deltamin= ',deltamin
      write(*,*)'Deltamax= ',deltamax
      write(*,*)'DepthMin: ',depthmin
      write(*,*)'DepthMax: ',depthmax
 
      stop 'No worries mate ...'
      end
 
* ---------------------------------------------------------------
      subroutine minmax(rmin,rmax)
      implicit real*8 (a-h, o-z)
      if(rmin.LE.rmax)return
      dest=rmin
      rmin=rmax
      rmax=dest
      return      
      end

      subroutine rayl(ichr,ichd,w,ni,delay)
      implicit real*8 (a-h, o-z)
      real*8 la,lo
      parameter (msg = 256)
      dimension  w(3,msg+1)

c     dinter - is the step for interpolation
c     finding the boundary of cell

      dinter=1.000
      rlength=0.0e10
      rdelay =0.0e10
      r2d = 90./asin(1.)
      dpi = asin(1.)/ 90.
* ibl is the flag that is writen in order to mark
* rays penetraiting lower 2900 km.
      ibl=0

c     write(*,*)'rayl'

c     Begin here for each point
      i=2
      do while(i.LE.ni+1)

c     Find location of two points

      call cellvel(w(2,i-1),w(3,i-1),w(1,i-1),Vel1,nbl1)
      call cellvel(w(2,i),w(3,i),w(1,i),Vel2,nbl2)

c       is it in the same block ? 
c       if Yes - find the parameters
c       and begin again for next two points (see go to 99)

         if(nbl1.EQ.nbl2)then
         
c           find length between points
            call dist(w(2,i-1),w(3,i-1),w(1,i-1),
     &                w(2,i),w(3,i),w(1,i),dst)

            rlength=rlength+dst
            rdelay =rdelay+dst/((Vel1+Vel2)/2.) 
            iwr=0
            go to 99

         endif
          
c    if they are not in the same block then
c    find the last point in the same block

         if(nbl1.NE.nbl2)then

         Vel2=Vel1

c        find last point within the block
c        firstly find the distance between points

            call dist(w(2,i-1),w(3,i-1),w(1,i-1),
     &                w(2,i),w(3,i),w(1,i),dst)

c        find number of steps for fixed length dinter

            nms=int(dst/dinter)
            if(nms.LE.0)nms=1

c        find dx, dy and dz for each step

                  ar=w(2,i-1)*dpi
                  as=w(2,i)*dpi
                  br=w(3,i-1)*dpi
                  bs=w(3,i)*dpi
                  rr=w(1,i-1)
                  rs=w(1,i)
                 x1 = rr*sin(ar)*cos(br)
                 y1 = rr*sin(ar)*sin(br)
                 z1 = rr*cos(ar)
                 x2 = rs*sin(as)*cos(bs)
                 y2 = rs*sin(as)*sin(bs)
                 z2 = rs*cos(as)
                 dx = (x2-x1)/nms
                 dy = (y2-y1)/nms
                 dz = (z2-z1)/nms

c         begin looking for the last point

          do j=1,nms

c         gradually increase distance toward the second point
          x = x1 + dx*j
          y = y1 + dy*j
          z = z1 + dz*j
          r = sqrt(x**2 + y**2 + z**2)
          
          acosa=z/r
          if(acosa.LT.-1.)acosa=-1.
          if(acosa.GT.1)acosa=1.
          la = acos(acosa)*r2d
          acosa=(x/r)/sin(la*dpi)
          if(acosa.LT.-1.)acosa=-1.
          if(acosa.GT.1.)acosa=1.
          lo = acos(acosa)*r2d
          if(y.LT.0.00000)lo=360.00000-lo



c         write(*,*)'tut?',i

c          where is it ?
            call cellvel(la,lo,r,Velc,nbl)      

c          if it is in the next cell begin to find parameters
c          for the last distance, redefine first point and
c          begin for new two points ( see go to 66)
c          if they still  are in the same block
c          increase distance more and performe approximation again

            if(nbl.NE.nbl1)then 

               call dist(w(2,i-1),w(3,i-1),w(1,i-1),
     &                la,lo,r,dst)

                   rlength=rlength+dst
                   rdelay =rdelay+dst/((Vel1+Vel2)/2.)

c       Write results into the file if ray is within grid net

                   if(nbl1.GT.0)then
                     write(ichr,*)nbl1,rlength
                     iwr=1
                   else
                     ibl=1
                   endif

                   rlength=0.0e10
                   w(2,i-1)=la
                   w(3,i-1)=lo
                   w(1,i-1)=r
                              go to 66
            endif

          Vel2=Velc

          enddo


         endif
99    continue
      i=i+1
66    continue

      enddo

      if((nbl1.GT.0).AND.(iwr.EQ.0))then
         write(ichr,*)nbl1,rlength
         iwr=1
      else
         ibl=1
      endif


c     The ray is finished
      write(ichr,*)'0 -1'
      write(ichd,*)(delay-rdelay),' ',ibl
      delay=delay-rdelay
      return
      end



      subroutine dist(rla1,rlo1,r1,rla2,rlo2,r2,dst)
      implicit real*8 (a-h, o-z)
      dpi = asin(1.)/ 90.

      as=rla1*dpi
      ar=rla2*dpi
      bs=rlo1*dpi
      br=rlo2*dpi
      rs=r1
      rr=r2
                 x1 = rr*sin(ar)*cos(br)
                 y1 = rr*sin(ar)*sin(br)
                 z1 = rr*cos(ar)
                 x2 = rs*sin(as)*cos(bs)
                 y2 = rs*sin(as)*sin(bs)
                 z2 = rs*cos(as)
                 dx = (x2-x1) 
                 dy = (y2-y1) 
                 dz = (z2-z1) 
                 dst=sqrt(dx**2+dy**2+dz**2)
      return
      end



      subroutine km2deg(ala,alo,adp,dx,dy,bla,blo,bdp)
      implicit real*8 (a-h, o-z)
c     This subroutine calculate position of new point
c     in polar coordinates basing on the coordinates
c     of main point in radians ( la is colatitude) and dx and dy in kilometers
      dpi = asin(1.)/ 90.
      dps=adp*SIN(ala)
      blo=alo+atan2(dx,dps)
      bla=ala+atan2(dy,adp)
        if(bla.gt.(180.*dpi))then
           bla=360.*dpi-bla
           blo=blo+180.*dpi
        endif
        if(bla.lt.0.)then
           bla=abs(bla)
           blo=blo+180.*dpi
        endif
        if(blo.lt.0.)blo=360.*dpi+blo
        if(blo.gt.(360.*dpi))blo=blo-(360.*dpi)
      bdp=sqrt(adp**2+dx**2+dy**2)
      return
      end

c                                                                       
c 3d pseudo-bending for the continuous, spherical earth
c   based on kazuki koketsu (eri, univ. tokyo) algorithm
c
c  ray is tracing from receiver to source
c  input coordinates are in degrees and km
c  latitude from -90 to 90 and longitude from -180 to 180
c  depth is radius from the earth's center

      subroutine pbr(aas, bbs, hs, aar, bbr, hr,w,ni,iraygeometry)    
      implicit real*8 (a-h, o-z)                                        
      parameter (msg = 256)
      dimension  r(msg+1), a(msg+1), b(msg+1)
      dimension  w(3,msg+1)
      integer ni,i
      external   vel,rtim,rlen
      real*8 RNULL
      common /coord/ shiftlo
      data      ro / 6371.00 /
     
      data    RNULL /0.0e10/
c     write(*,*)aas, bbs, hs, aar, bbr, hr
c                                                                       
c       parameters for calculation                                      
c         xfac   = enhancement factor (see um and thurber, 1987)        
c         nloop  = number of bending iterations
c         n1, n2 = min & max of ray segments
c         mins  = min. length of segment (km)
c       initialization                                                  

c      flim   = 1.e-2/100. = 0.0001 sec
c      flim   = 1.e-2/10. = 0.001 sec

      data Ddepth /3481./
      ni     = msg+1
      xfac   = 1.7
      n1     = 2
      n2     = 256
      nloop  = 200
      flim   = 1.e-2/100.
      mins   = 33.3


      dpi = asin(1.)/ 90.
      r2d = 90./asin(1.)

c      Check coordinates
      if(aas.LT.-90.OR.aas.GT.90.)then
        write(*,*)'Latitude of source is out of range'
        stop
      endif
      if(aar.LT.-90.OR.aar.GT.90.)then
        write(*,*)'Latitude of station is out of range'
        stop
      endif
      if(bbs.LT.-180.OR.bbs.GT.180.)then
        write(*,*)'Longitude of source is out of range'
        stop
      endif
      if(bbr.LT.-180.OR.bbr.GT.180.)then
        write(*,*)'Longitude of station is out of range'
        stop
      endif


c     Rotate coordinates in order to have
c     longitude and latitude range from 0 to 180. 
c     This program does not work with angles
c     greater than 180.       

c     Pass from latitude to colatitude

      as = (90.00-aas) * dpi
      ar = (90.00-aar) * dpi

      if(bbr.LT.0.0)then
         bre=360.+bbr
      else
         bre=bbr
      endif

      if(bbs.LT.0.0)then
         bso=360.+bbs
      else
         bso=bbs
      endif
      dlo=abs(bso-bre)

      if(dlo.LT.180.)then


c     write(*,*)dlo

      shiftlo=0.0e10
      if(bso.LT.bre)then
         shiftlo=bso-(180.-dlo)/2.
         bbs=(180.-dlo)/2.
         bbr=bbs+dlo
      else
         shiftlo=bre-(180.-dlo)/2.
         bbr=(180.-dlo)/2.
         bbs=bbr+dlo
      endif

      else

      dlo=360.0000-dlo

      shiftlo=0.0e10
      if(bso.LT.bre)then
         shiftlo=bso-(dlo+(180.-dlo)/2.)
         bbs=(180.-dlo)/2.+dlo 
         bbr=bbs-dlo
      else    
         shiftlo=bre-(dlo+(180.-dlo)/2.)
         bbr=(180.-dlo)/2.+dlo
         bbs=bbr-dlo
      endif
 

      endif
c        write(*,*)'bbr,bbs,shift'
c        write(*,*)'aar,aas'
c        write(*,*)aar,aas
c        write(*,*)'bbr,bbs,shiftlo'
c        write(*,*)bbr,bbs,shiftlo

         bs = bbs * dpi
         br = bbr * dpi  

      ad = (as + ar) / 2.                                               
      rs = ro - hs
      rr = ro + hr                                                      
c
c *** initial straight ray ***                                           
c       ni : number of ray segments
      ni = n1
      x1 = rr*sin(ar)*cos(br)                                       
      y1 = rr*sin(ar)*sin(br)                                       
      z1 = rr*cos(ar)                                               
      x2 = rs*sin(as)*cos(bs)                                       
      y2 = rs*sin(as)*sin(bs)                                       
      z2 = rs*cos(as)       
      dx = (x2-x1) / ni
      dy = (y2-y1) / ni
      dz = (z2-z1) / ni
c     write(*,*)(x2-x1),(y2-y1),(z2-z1)
c     write(*,*)sqrt((x2-x1)**2+(y2-y1)**2)
      do 155  j=1,ni+1
          x = x1 + dx*(j-1)
          y = y1 + dy*(j-1)
          z = z1 + dz*(j-1)                                             
          r(j) = sqrt(x**2 + y**2 + z**2)                         
          acosa=z/r(j)
          if(acosa.LT.-1.)acosa=-1.
          if(acosa.GT.1)acosa=1.
          a(j) = acos(acosa)                                   
          acosa=x/r(j)/sin(a(j))
          if(acosa.LT.-1.)acosa=-1.
          if(acosa.GT.1)acosa=1.
          b(j) = acos(acosa)
          if(y.LT.0.00000)b(j)=360.00000*dpi-b(j)

cc    write(*,*)a(j)*r2d,b(j)*r2d,r(j)
 155  continue
c
      to = rtim(ni+1,r,a,b)
      tp = to
c      write(6,*) to
      do 112  i=1,ni+1
          w(1,i) = r(i)
          w(2,i) = a(i)
          w(3,i) = b(i)
 112  continue
c *** number of points loop ***
      loops = 0
      do while(ni .le. n2)
c                                                                       
c *** interation loop ***                                               
          do 3000  l=1,nloop                                                 
              loops = loops + 1
              do 2250  kk=2,ni
c               see um & thurber (1987) p.974.
c              write(*,*)'Tut 1'
                  if(mod(kk,2) .eq. 0) then
                      k = kk/2 + 1
                  else
                      k = ni+1 - (kk-1)/2
                  endif
                  r1 = r(k-1)                                 
                  a1 = a(k-1)                                 
                  b1 = b(k-1)                                 
                  x1 = r1*sin(a1)*cos(b1)                           
                  y1 = r1*sin(a1)*sin(b1)                           
                  z1 = r1*cos(a1)                                   
                  r3 = r(k+1)                                 
                  a3 = a(k+1)                                 
                  b3 = b(k+1)                                 
                  x3 = r3*sin(a3)*cos(b3)                           
                  y3 = r3*sin(a3)*sin(b3)                           
                  z3 = r3*cos(a3)                                   
                  dx = x3 - x1                                      
                  dy = y3 - y1                                      
                  dz = z3 - z1                                      
                  x2 = x1 + dx/2                                    
                  y2 = y1 + dy/2                                    
                  z2 = z1 + dz/2                                    
                  r2 = sqrt(x2**2 + y2**2 + z2**2)                  
                  acosa=z2/r2
                  if(acosa.LT.-1.)acosa=-1.
                  if(acosa.GT.1)acosa=1.
                  a2 = acos(acosa)
                  sina = sin(a2)                                    
                  cosa = cos(a2)                                    
                  acosa=x2/r2/sina
                  if(acosa.LT.-1.)acosa=-1.
                  if(acosa.GT.1)acosa=1.
                  b2 = acos(acosa)                          
                  if(y.LT.0.00000)b2=360.00000*dpi-b2
c
c
                  dn = dx**2 + dy**2 + dz**2                        
                  ddn = sqrt(dn)                                    
                  dr = (r3-r1) / ddn                                
                  da = (a3-a1) / ddn                                
                  db = (b3-b1) / ddn

c  Begin find the gradients and velocities

c                 first find the length of segment
                  dseg=sqrt((dx/2)**2+(dy/2)**2+(dz/2)**2)

                  ddseg=dseg/2.
c                 write(*,*)'dseg ',dseg 

c                 Now ddseg will be a distance to find dV
c                 along the coordinates 

c                 Determine velocity at 3 points
c                 write(*,*)'before',ni
c                 write(*,*)
c           write(*,*)'Pered opredeleniem glubini ------------------'
c           do j=1,ni+1
c           write(*,*)j,ni,r(j),a(j)*r2d,b(j)*r2d
c           enddo
c                 write(*,*)'grohnulas pri',ni,' segmentov'
c                 write(*,*)'l,k ',l,k
c                 write(*,*)r1,a1*r2d,b1*r2d


                  v1 = vel(r1,a1,b1,1)                            

c                 write(*,*)'V1= ',v1
c                 write(*,*)r2,a2*r2d,b2*r2d

                  v2 = vel(r2,a2,b2,2)                            

c                 write(*,*)'V2= ',v2
c                 write(*,*)r3,a3*r2d,b3*r2d

                  v3 = vel(r3,a3,b3,3)                            

c                 write(*,*)'V3= ',v3

c                 Begin to determine coordinates
c                 of pints surroundibg point a2,b2,r2
c                 at the distance ddseg

                  upz = r2+ddseg
                  dwz = r2-ddseg

                   if(upz.gt.ro)then
                      upz=ro
                      dwz=upz-dseg
                   endif

                   if(dwz.le.0.)then
                      dwz=0.00000001
                      upz=ro
                   endif

C Comment out the following if-endif for SKS & P: It's just for S and Pdiff
                   if(iraygeometry.LE.2)then
                   if(dwz.le.Ddepth)then
                      dwz=Ddepth
                      upz=dwz+dseg
                   endif
                   endif
C -------------------------------------------------


c                 find dV along depth
                  vr1 = vel(upz,a2,b2,4)
                  vr2 = vel(dwz,a2,b2,5)
                  vr=(vr1-vr2)/dseg

c                 write(*,*)'vr= ',vr

c                 find dV along longitude

                  call  km2deg(a2,b2,r2,ddseg,RNULL,adV,bdV,rdV) 
                  vb2 = vel(rdV,adV,bdV,6)
                  call  km2deg(a2,b2,r2,-1.*ddseg,RNULL,adV,bdV,rdV)
                  vb1 = vel(rdV,adV,bdV,7)
                  vb=-1.*(vb1-vb2)/dseg

c                 write(*,*)'vb= ',vb

c                 find dV along latitude

                  call  km2deg(a2,b2,r2,RNULL,ddseg,adV,bdV,rdV)
                  va2 = vel(rdV,adV,bdV,8)
                  call  km2deg(a2,b2,r2,RNULL,-1.*ddseg,adV,bdV,rdV)
                  va1 = vel(rdV,adV,bdV,9)
                  va=-1.*(va1-va2)/dseg

c                 write(*,*)'va= ',va
c                 write(*,*)'**************'
               
c               write(*,*)'after',ni


c           do j=1,ni+1
c           write(*,*)j,ni,r(j),a(j)*r2d,b(j)*r2d
c           enddo
      
			             
c               spherical
c                   velocity gradient
c                 va = va / r2
c                 vb = vb / r2 / sina
c                   (tangential vector) = (slowness vector) / s
c                 write(*,*)'perturbations !!!'
                  pr = dr
                  pa = r2 * da
                  pb = r2 * sina * db
c                 write(*,*)pr,pa,pb
                  vrd = pr*vr + pa*va + pb*vb
c                 write(*,*)vrd
                  rvr = vr - vrd*pr
                  rva = va - vrd*pa
                  rvb = vb - vrd*pb
c                 write(*,*)rvr,rva,rvb
                  rvs = sqrt(rvr*rvr + rva*rva + rvb*rvb)               
c                 write(*,*)rvs
c                 write(*,*)'posle rvs'
                  if(rvs .eq. 0.) then                              
                      r(k) = r2                                   
                      a(k) = a2                                   
                      b(k) = b2                                   
                  else                                              
                      rvr = rvr / rvs                                 
                      rva = rva / rvs                                 
                      rvb = rvb / rvs                                 
c                     write(*,*)rvr,rva,rvb
                      cc   = (1./v1+1./v3)/2.                          
c                     write(*,*)'cc ',cc
                      rcur = vr*rvr + va*rva + vb*rvb               
c   Tut esli rcur < 0.0 proishodit hernia
c   poetomu postavlen abs. Ne yasno mozhno li eto delat
c   ili net no rabotaet. Obichno oshibka poyavliaetsia
c   ochen redko v nekotorih tochkah
c v etom sluchae abs prosto ne daet oshibki y posledniaya iteraciya
c  uzhe ne imeet rcur negativnim y podgoniaet normalno reshenie
c  ( mozhet bit)
c                     write(*,*)'rcur ', rcur
                      if(rcur.LE.0.0)then
c                     write(*,*)'Negative'
                      rcur=abs(rcur)
                      endif
                      rcur = (cc*v2+1.) / (4.*cc*rcur)                
c                     write(*,*)'rcur ',rcur
                      rcur = -rcur + sqrt(rcur**2+dn/(8.*cc*v2))     
c                     write(*,*)'rcur ',rcur
                      rdr = rvr * rcur
                      rda = rva * rcur
                      rdb = rvb * rcur
c                     write(*,*)'rdr,rda,rdb ',rdr,rda,rdb
                      rp  = r2 + rdr
                      ap  = a2 + rda/r2
                      bp  = b2 + rdb/(r2*sina)
c                     write(*,*)rp,ap,bp
                      r(k) = (rp-r(k))*xfac + r(k)
                      a(k) = (ap-a(k))*xfac + a(k)
                      b(k) = (bp-b(k))*xfac + b(k)
              
* Fixing ray geometry

*                    Prevent ray penetration beyond D''
                     if(iraygeometry.LE.2)then
                        if(r(k).le.Ddepth)then
                        r(k)=Ddepth
                        endif
                     endif

*                    Fix middle point for PP at the surface

                     if(iraygeometry.EQ.4)then
                       if(k.EQ.(ni/2+1))then
                          r(k)=ro
                       endif
                     endif

                  endif
c             write(*,*)'after vector',ni
c             write(*,*)r(k),a(k)*r2d,b(k)*r2d
c             write(*,*)'++++++++++++++++++++++++++++++++++++++'
 2250         continue                                              
c                                                                       
c             write(*,*)'pered ciclom', ni
              idstn=ni
              do 212  j=1,ni+1
                  w(1,j) = r(j)
                  w(2,j) = a(j)
                  w(3,j) = b(j)
c             write(*,*)j,ni,w(1,j),w(2,j)*r2d,w(3,j)*r2d
 212          continue
              ni=idstn
c              here trace each iteration
c             write(*,*)'pered rtim',ni
              tk = rtim(ni+1,r,a,b)
           
c             write(*,*)to,tk,ni
c              write(6,*) ni, l, (tk-atim)/atim*100.
c               tk > to $b$n;~$bdd;_$5$;$k$?$a(b abs $b$oiu$1$j$$!%(b
              if(abs(to-tk) .le. to*flim)  go to 310                        
              to = tk                                                       
 3000     continue                                                          
 310      continue   
              to=tk
c          skip increasing of segment number if minimum length
c          of segment is exceed or maximum number of segments
c          was reached

           if(dseg.lt.mins.or.ni.ge.n2) go to 66666
c
c           double the number of points.
c           write(*,*)'double number'

          ni = ni * 2
         
c         write(*,*)ni
          do i=1,ni/2+1
              r(i*2-1) = w(1,i)
              a(i*2-1) = w(2,i)
              b(i*2-1) = w(3,i)
          enddo
          do k=2,ni,2
              r1 = r(k-1)
              a1 = a(k-1)                                 
              b1 = b(k-1)                                 
              x1 = r1*sin(a1)*cos(b1)                           
              y1 = r1*sin(a1)*sin(b1)                           
              z1 = r1*cos(a1)                                   
              r3 = r(k+1)                                 
              a3 = a(k+1)                                 
              b3 = b(k+1)                                 
              x3 = r3*sin(a3)*cos(b3)                           
              y3 = r3*sin(a3)*sin(b3)                           
              z3 = r3*cos(a3)                                   
              dx = x3 - x1                                      
              dy = y3 - y1                                      
              dz = z3 - z1                                      
              x2 = x1 + dx/2                                    
              y2 = y1 + dy/2                                    
              z2 = z1 + dz/2                                    
              r2 = sqrt(x2**2 + y2**2 + z2**2)                  
              acosa=z2/r2
              if(acosa.LT.-1.)acosa=-1.
              if(acosa.GT.1)acosa=1.
              a2 = acos(acosa)
              sina = sin(a2)                                    
              acosa=x2/r2/sina
              if(acosa.LT.-1.)acosa=-1.
              if(acosa.GT.1)acosa=1.
              b2 = acos(acosa)
              if(y.LT.0.00000)b2=360.00000*dpi-b2
              r(k) = r2
              a(k) = a2
              b(k) = b2
          enddo
          tk = rtim(ni+1,r,a,b)
c           tk > tp $b$n;~$bdd;_$5$;$k$?$a(b abs $b$oiu$1$j$$!%(b
c here i change tp and put to
          if(abs(to-tk) .le. to*flim)  go to 99999
          to = tk 
      enddo
c                                                                       
99999 continue
c     write(*,*)'poslednii',ni
      idstn=ni
      do i=1,ni+1
          w(1,i) = r(i)
          w(2,i) = a(i)
          w(3,i) = b(i)
      enddo
      ni=idstn
66666 continue
c Return coordinates to the origin
      idstn=ni
      do k=1,ni+1
          w(2,k) = w(2,k)*r2d
          w(3,k) = w(3,k)*r2d+shiftlo
       if(w(3,k).lt.0.)w(3,k)=360.+w(3,k)
c     write(*,*)k,ni,w(1,k),w(2,k),w(3,k)
      enddo 
      ni=idstn
c     do i=1,ni+1
c     write(*,*)w(1,i),w(2,i),w(3,i)
c     enddo
c     write(*,*)'Ray is done with ',ni,' segments'
      return
      end                                                               
c                                                                       
c
      function rlen(n, r, a, b)                                         
      implicit real*8 (a-h, o-z)                                        
      parameter (msg = 256)
      integer n
      dimension  r(msg+1), a(msg+1), b(msg+1)
      rlen = 0.                                                         
      do 100  i=1,n-1                                                   
          j  = i + 1                                                    
          xi = r(i)*sin(a(i))*cos(b(i))                                 
          yi = r(i)*sin(a(i))*sin(b(i))                                 
          zi = r(i)*cos(a(i))                                           
          xj = r(j)*sin(a(j))*cos(b(j))                                 
          yj = r(j)*sin(a(j))*sin(b(j))                                 
          zj = r(j)*cos(a(j))                                           
          dl = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2                     
          rlen = rlen + sqrt(dl)                                        
  100 continue                                                          
      return                                                            
      end                                                               
c
c                                                                       
      function rtim(m, r, a, b)                      
      implicit real*8 (a-h, o-z)                                        
      parameter (msg = 256)
      dimension  r(msg+1), a(msg+1), b(msg+1)
      integer m
      if(m.GT.(msg+1))write(*,*)'*'
      rtim = 0.
      rv1 = 1./vel(r(1),a(1),b(1),-1)
      do 100  j=1,m-1
          x1 = r(j)*sin(a(j))*cos(b(j))
          y1 = r(j)*sin(a(j))*sin(b(j))
          z1 = r(j)*cos(a(j))                                               
          x2 = r(j+1)*sin(a(j+1))*cos(b(j+1))
          y2 = r(j+1)*sin(a(j+1))*sin(b(j+1))
          z2 = r(j+1)*cos(a(j+1))
          dl = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
          rv2 = 1./vel(r(j+1),a(j+1),b(j+1),-2)
          sm = (rv1 + rv2) / 2.         
          if(dl.LE.0.00000000000)dl=1.00000000000
          rtim = rtim + sqrt(dl)*sm
          rv1 = rv2
  100 continue                                                          
      return                                                            
      end                                                               
c
c
      function vel(r, pa, ra, k)                              
      implicit real*8 (a-h, o-z)                                    
      common /coord/ shiftlo
      r2d = 90./asin(1.)                                                        

c     convert to degree and rotate coordinates 
      clat=pa*r2d
      clon=ra*r2d+shiftlo
      if(clon.lt.0.)clon=360.+clon

c     write(*,*)k
c     if(k.EQ.-1)write(*,*)clat,clon,r
      call cellvel(clat,clon,r,V,nbl)
c     write(*,*)'exit vel'
      if(V.LE.0.000000)then
      write(*,*)V,' at ',nbl,clon,clat,r
      endif
      vel=V

      return
      end


      subroutine cellvel(rla,rlo,era,Vl,nbl)
      implicit real*8 (a-h, o-z)
      real*8 la,lo,era,dp
c This subroutine determine the velocity in the
c block using co Lat,Lon and R.
c   Lat and Lon in degrees
c   era - is the radius
c   dp is the depth
      parameter (nunkn=500000)
      parameter(npcom=440000,nlayer=100)
c     Note, the coordinates are colongitude and colatitude
      parameter (nlayeri=100)
      dimension icellnum(npcom)
      common /param/dx,dy,rmilo,rmalo,rmila,
     &rmala,dxi,dyi,numblock,nx,ny,nz,nxi,nyi,nzi,icellnum


      integer  num
      dimension r(npcom),v(npcom),d(npcom)
      dimension r1(npcom),v1(npcom),d1(npcom)
      dimension xlayer(0:nlayer),dz(nlayer),dzsum(nlayer)
      dimension xlayeri(0:nlayeri)
      common /refmodel/ np15

      common /value/ pi,con,rcon,deg,r,v,d,
     >r1,v1,d1,dz,dzsum

      common/model/xlayer,xlayeri
      common /perturb/ xan(nunkn)

      la=rla
      lo=rlo

      Re=6371.00
      dp=Re-era
     
      nlr=0
      nla=0
      nlo=0

c     Additional parameters for zone of study

      dxi=((rmalo-rmilo)/float(nxi))
      dyi=((rmala-rmila)/float(nyi))

      if(dp.LT.0.00)then
c       dp=Re+dp
        dp=0.0
      endif

      do while(lo.GE.360.)
         lo=lo-360.
      enddo

      do while(lo.LT.0.0e10)
         lo=lo+360.00
      enddo

      if(la.LT.0.000)then
         la=abs(la)
      endif

      if(la.GT.180.)then
         la=360.-la
      endif

      if(la.EQ.180.)then
         nla=nla-0.0001000
      endif

c     Firstly determine the number of layer within 1D velocity model
      if(dp.LT.0.0)then
      ilr=1
      else
      do i=2,np15
 
       if(dp.GE.d(i-1).AND.dp.LT.d(i))then
           ilr=i-1
           go to 5
       endif    

      enddo 
5     continue
      endif

****  IF   *******************************************************
c     Where is the point - inside/outside of the region
c     if(dp.GT.1600.)write(*,*)'------'

      if(lo.GE.rmilo.AND.lo.LE.rmalo.AND.la.GE.rmila
     &.AND.la.LE.rmala.AND.dp.LE.xlayeri(nzi))then

c      inside of the study region

      if(dp.GT.xlayeri(nzi))then
      write(*,*)'Velocity determination
     & is deeper than regional nodes, go to the global scale'
c     go to  24
      endif

c     number of the layer in the zone of study
      do i=1,nzi

       if(dp.GE.xlayeri(i-1).AND.dp.LT.xlayeri(i))then
           nlr=i
           go to 21
       endif

      enddo
21     continue

c     Then Column (lon) number

      nlo=int((lo-rmilo)/dxi)+1

c     and then Row (lat) number

      nla=int((la-rmila)/dyi)+1
      num=(nla-1)*nxi+nlo+(nlr-1)*nxi*nyi
      Vl=v(ilr)*(1.+xan(num))
      num=icellnum(num)
      nbl=num 


** ELSE *************************************** ELSE
      else
c           off the region of study

c24     continue


      if(dp.GT.xlayer(nz))then
c     write(*,*)'Velocity determination
c    &                               is deeper than nodes'
      go to 4
      endif

*     What layer in the global net
      do i=1,nz

       if(dp.GE.xlayer(i-1).AND.dp.LT.xlayer(i))then
           nlr=i
           go to 1
       endif

      enddo 
1     continue

c     Fix boundaries

      if(dp.ge.xlayer(nz))nlr=nz
      if(dp.le.0)nlr=1

c     Then Column (lon) number

      nlo=int(lo/dx)+1

c     and then Row (lat) number

      nla=int(la/dy)+1

      if(nlr.LE.0.OR.nlr.GT.nz)then
         write(*,*)dp,'NLR: ',nlr
         stop 'error of depth'
      endif
      if(nlo.LE.0.OR.nlo.GT.nx)then
         write(*,*)lo,'NLO: ',nlo
         stop 'error of LON'
      endif
      if(nla.LE.0.OR.nla.GT.ny)then
         write(*,*)la,'NLA: ',nla
         stop 'error of LAT'
      endif

c     The number of cell in global scale will be
    
      num=(nla-1)*nx+nlo+(nlr-1)*nx*ny+numblock

c       write(*,*)'global',num
      Vl=v(ilr)*(1.+xan(num))
      nbl=num

      endif
******************* ENDIF ******************

c     and velocity in this block

      return
4     continue
      Vl=v(ilr)
      nbl=-1
      return
      end

      subroutine cderiv(w,ni,dth2,dth3,dth4)
      implicit real*8 (a-h, o-z)
      parameter (msg = 256)
      dimension  w(3,msg+1)
      ve=6.2
      nbl=0
      dpi = asin(1.)/ 90.
      r2d = 90./asin(1.)
      r0        = 6371.00
      epc       = 1.0e-5
      nrp       = ni+1
      pe        = w(2,nrp)*dpi
      re        = w(3,nrp)*dpi
      he        = w(1,nrp)
      nrp1      = nrp-1
      pe1       = w(2,nrp1)*dpi
      re1       = w(3,nrp1)*dpi
      he1       = w(1,nrp1)
      call conaz(re,pe,re1,pe1,des,az)
      des       = des*(he/r0)
      dhe       = he1-he
      adh       = abs(dhe)
      call cellvel(pe*r2d,re*r2d,he,ve,nbl)
      if(des.lt.epc.and.adh.lt.epc) dhe = 1.0
      th        = atan2(des,dhe)
      dtddel    = (sin(th)*he/ve)
c     write(*,*)'rp= ',dtddel
      dtdr      = abs(cos(th)/ve)
c     write(*,*)'az= ',az*r2d
      if(az*r2d.GT.360.)az=az-360.*dpi
      dth3 =-dtddel*cos(az)
      dth4 = dtddel*sin(az)*sin(pe)
      dth2 =-sqrt((he/ve)**2-dtddel**2)/he
      return
      end

      subroutine conaz(re,pe,rs,ps,del,az)
      implicit real*8 (a-h, o-z)
      data r0,eps,pi15/6371.00,1.0e-7,4.712389/
      sps  = sin(ps)
      cps  = cos(ps)
      spe  = sin(pe)
      cpe  = cos(pe)
      ses  = sin(re-rs)
      ces  = cos(re-rs)
      x    = sps*ses
      y    = cpe*sps*ces - spe*cps
      s    = sqrt(x*x + y*y)
      c    = spe*sps*ces + cpe*cps
      del  = atan2(s,c)*r0
      az   = 0.0
      ax   = abs(x)
      ay   = abs(y)
      if(ax.lt.eps.and.ay.lt.eps) return
      az   = pi15-atan2(y,x)
      return
      end
