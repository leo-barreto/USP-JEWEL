C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C++          JEWEL Add-On to Read v-USPhydro 2+1D Profiles          ++
C++                                                                 ++
C++ The program is part of the developed interface between the      ++
C++ parton propagation of JEWEL and an external hydrodynamic 2+1D   ++
C++ medium profile, intended for v-USPhydro.                        ++
C++                                                                 ++
C++ This code implements multiple auxiliary functions for reading   ++
C++ and interpolating data of a medium profile.                     ++
C++                                                                 ++
C++                                                                 ++
C++ Created by:                                                     ++
C++  - Fabio M. Canedo [fabio.canedo@usp.br]                        ++
C++  - Leonardo Barreto [leonardo.barreto.campos@usp.br]            ++
C++  Instituto de Fisica, Universidade de Sao Paulo, Brazil         ++
C++  2019                                                           ++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine reader(filename,np,nt,timesteps,tprofile,u,theta)
      implicit none
      integer i,j,k
      integer np,nt
      integer geti
      integer ios
      character*100 filename
      double precision timesteps(60)
      double precision tprofile(np,np,60)
      double precision u(np,np,60)
      double precision theta(np,np,60)
      double precision t,x,y,temp,vx,vy
      double precision pi
      DATA PI/3.141592653589793d0/
      integer linecounter

      open(unit=1,file=filename,iostat=ios)
      write(*,*) "Opening file: ",filename
      
      do k=2,60
      timesteps(k)=1.d20
      end do
      k=1
      timesteps(k)=0.d0

      linecounter=0
      do while (ios.eq.0)
      
      read(1,*,iostat=ios) t,x,y,temp,vx,vy
      i=geti(x,np)
      j=geti(y,np)

      linecounter=linecounter+1

      if(t.ne.timesteps(k)) then
      k=k+1
      timesteps(k)=t
      end if


      tprofile(i,j,k)=temp
      u(i,j,k)=sqrt(vx**2+vy**2)
      if(vx.ne.0.d0) then
            theta(i,j,k)=atan(vy/vx)
      else
            if(vy.gt.0.d0) then
                  theta(i,j,k)=pi/2.d0
            else
                  theta(i,j,k)=-pi/2.d0
            end if
      end if

      end do
      
      write(*,*) "File has ",linecounter, " lines"

      !Stop simulations if file has not enough lines (arbitrary)
      if (linecounter.lt.10) then
           STOP
      end if

      end subroutine



      ! Create matrix of probabilities given
      ! temperature at tau0 and initvtxmap.dat
      subroutine read_initvtx(np, tprofile, vtxprofile)
      implicit none
      COMMON/logfile/logfid
      INTEGER logfid
      integer np, ios, i, j, nlines, k
      logical f_exist
      double precision tprofile(np, np, 60)
      double precision vtxprofile(np, np), tempmap(2000), vtxmap(2000)
      double precision temp, dT, dsup, dinf
C--hydro auxiliary files
      COMMON /HYDROF/ INITVTXF, IDEALN0, CUSTOMN0F
      CHARACTER*200 INITVTXF, CUSTOMN0F
      LOGICAL IDEALN0


      ! Load initial vertex map from INITVTXF
      tempmap = 0.0
      vtxmap = 0.0
      nlines = 1

      write(*,*) 'Loading initial vertex transformation from ', INITVTXF
      write(logfid,*) 'read initial vertex table from ', INITVTXF
      inquire(file=INITVTXF, exist=f_exist)
      if (f_exist .eqv. .false.) then
          write(*,*) 'INITVTXF does not exist. Killing simulation.'
          stop
      end if      

      open(1, file=INITVTXF, iostat=ios)
      do while (ios .eq. 0)
          read(1, *, iostat=ios) tempmap(nlines), vtxmap(nlines)
          nlines = nlines + 1
      end do

      close(1)


      if (nlines .gt. 2001 .or. nlines .lt. 2) then
          write(*,*) 'WARNING: init vertex table size is wrong.
     &Check read_initvtx'
      end if

      ! Transform tprofile IC (tau0) into hard scattering
      ! probability from linear interpolation of vertex map

      ! Check boundaries
      if ((maxval(tprofile) .gt. maxval(tempmap)) .or.
     &(minval(tprofile) .lt. minval(tempmap))) then
          write(*,*) 'Temperature (', temp, ' GeV) out of bounds
     &in initvtx. Check initvertex table.'
          stop
      end if


      do i = 1, np
          do j = 1, np
              temp = tprofile(i, j, 2)

              ! No temperature => no entropy
              if (temp .eq. 0) then
                  vtxprofile(i, j) = 0.d0

              else

                  ! Check where temp is found in tempmap
                  do k = 2, nlines
                      if (tempmap(k) .gt. temp) then
                          dT = tempmap(k) - tempmap(k - 1)
                          dsup = vtxmap(k - 1) * (tempmap(k) - temp)
                          dinf = vtxmap(k) * (temp - tempmap(k - 1))
                          ! Interpolate
                          vtxprofile(i, j) = (dsup + dinf) / dT

                          ! Only do this for the first k, since temp is
                          ! increasing with k
                          exit
                      end if
                  end do
              end if
          end do
      end do

      ! Normalize probabilities
      vtxprofile = vtxprofile / sum(vtxprofile)

      end subroutine



      ! Read array for n0 function interpolation 
      subroutine read_customn0(tempmap, n0array)
      implicit none
      COMMON/logfile/logfid
      INTEGER logfid
      integer np, ios, i, j, nlines, k
      logical f_exist
      double precision tempmap(2000), n0array(2000)
C--hydro auxiliary files
      COMMON /HYDROF/ INITVTXF, IDEALN0, CUSTOMN0F
      CHARACTER*200 INITVTXF, CUSTOMN0F
      LOGICAL IDEALN0


      ! Load n0array from CUSTOMN0F
      tempmap = 0.0
      n0array = 0.0
      nlines = 1

      if (IDEALN0) then
          write(logfid,*) 
     &'Assuming ideal number density (propto T ** 3).'
      else 
          write(*,*) 'Loading n0 custom function from ', CUSTOMN0F
          write(logfid,*) 'read n0 interpolation table from ', CUSTOMN0F
          inquire(file=CUSTOMN0F, exist=f_exist)
          if (f_exist .eqv. .false.) then
              write(*,*) 'CUSTOMN0F does not exist. Killing simulation.'
              stop
          end if      

          open(1, file=CUSTOMN0F, iostat=ios)
          do while (ios .eq. 0)
              read(1, *, iostat=ios) tempmap(nlines), n0array(nlines)
              nlines = nlines + 1
          end do

          close(1)


          if (nlines .gt. 2001 .or. nlines .lt. 2) then
              write(*,*) 'WARNING: custom n0 table size is wrong.
     &Check read_customn0'
          end if
      end if

      end subroutine



      integer function geti(x,np)
      implicit none
      integer np
      double precision x,xmin,xmax,dx
      
      xmax=25.d0
      xmin=-25.d0
      dx=(xmax-xmin)/(np-1)
      geti=1+(x-xmin)/dx
      end function



      integer function getk(t,timesteps)
      implicit none
      integer k
      double precision t
      double precision timesteps(60)
      
      getk=1
      do k=1,60

      if(timesteps(k).le.t) then
      getk=k
      endif

      enddo
      end function



      double precision function interpol(t,x,y,np,timesteps,tgrid,norm)
      implicit none
      integer i,j,ii,jj,iii,jjj,np
      integer k,kk
      logical norm
      double precision timesteps(60)
      double precision tgrid(np,np,60),igrid(4,4),xgrid(4,4),ygrid(4,4)
      double precision xmax,xmin,dx,dt
      double precision t,x,y,xa,ya
      double precision f(2)
      double precision bicubic
      integer getk

      COMMON/GAMMAMAX/GAMMAMAXIMUM,VELMAXIMUM
      DOUBLE PRECISION GAMMAMAXIMUM,VELMAXIMUM

      k=getk(t,timesteps)
      if(timesteps(k).eq.0.d0.and.timesteps(k+1).eq.0.d0) then
      dt=1e30
      else
      dt=timesteps(k+1)-timesteps(k)
      end if
      
      xmax=25.d0
      xmin=-25.d0
      dx=(xmax-xmin)/(np-1)

      i=1+floor((x-xmin)/dx)
      j=1+floor((y-xmin)/dx)
      xa=mod(x-xmin,dx)
      ya=mod(y-xmin,dx)

      do kk=1,2
      !write(*,*) "k:",k-1+kk
      do ii=1,4
            do jj=1,4
                  iii=i+ii-2
                  jjj=j+jj-2
                  if (iii.gt.np.or.iii.lt.1) then
                        igrid(ii,jj)=0.d0
                  else if (jjj.gt.np.or.jjj.lt.1) then
                        igrid(ii,jj)=0.d0
                  else if (k-1+kk.gt.59.or.k-1+kk.lt.1) then
                        igrid(ii,jj)=0.d0
                  else
                        igrid(ii,jj)=tgrid(iii,jjj,k-1+kk)
                  end if
                  xgrid(ii,jj)=xa+(ii-1)*dx
                  ygrid(ii,jj)=ya+(jj-1)*dx
            enddo
      enddo

      f(kk)=bicubic(xa,ya,dx,xmin,xmax,igrid,xgrid,ygrid,norm)
      enddo


      interpol=f(1)+(t-timesteps(k))*(f(2)-f(1))/dt
      if (norm .and. interpol.gt.VELMAXIMUM) then
        !write(*,*) "V > VELMAXIMUM", interpol, VELMAXIMUM
        interpol=VELMAXIMUM
      end if
      if (norm .and. interpol.lt.0) then
        interpol=0.d0
        !write(*,*) "Negative velocity norm"
      end if
      !interpol=0.1
      end function
     

 
      double precision function 
     &       bicubic(xa,ya,dx,xmin,xmax,igrid,xc,yc,norm)
      implicit none
      integer i,j
      logical norm
      double precision y(2,2),y1(2,2),y2(2,2),y12(2,2)
      double precision igrid(4,4)
      double precision xc(4,4),yc(4,4),ansy
      double precision xa,ya,dx,xmin,xmax
      double precision dertwospline
      
      do i=1,2
            do j=1,2
                  y(i,j)=igrid(i+1,j+1)
                  y1(i,j)=(igrid(i+2,j+1)-igrid(i,j+1))/(xc(i+2,j+1)
     &-xc(i,j+1))
                  y2(i,j)=(igrid(i+1,j+2)-igrid(i+1,j))/(yc(i+1,j+2)
     &-yc(i+1,j))
                  y12(i,j)=(igrid(i+2,j+2)-igrid(i+2,j)-igrid(i,j+2)
     &+igrid(i,j))/(yc(i+1,j+2)-yc(i+1,j))*(xc(i+2,j+1)-xc(i,j+1))
            end do
      end do
      
      bicubic=dertwospline(xc(2:3,2:3),yc(2:3,2:3),y,y1,y2,y12,
     &xa,ya,norm)
      end function



      double precision function
     &   dertwospline(x1,x2,y,y1,y2,y12,xa,xb,norm)
      implicit none
      integer i
      logical norm
      double precision x1(2,2)
      double precision x2(2,2)
      double precision y(2,2)
      double precision y1(2,2)
      double precision y2(2,2)
      double precision y12(2,2)
      double precision w(2),w1(2)
      double precision xa,xb
      double precision derspline

      w(1)=derspline(x1(:,1),y(:,1),y1(:,1),xa,norm)
      w(2)=derspline(x1(:,2),y(:,2),y1(:,1),xa,norm)
      w1(1)=derspline(x1(:,1),y2(:,1),y12(:,1),xa,norm)
      w1(2)=derspline(x1(:,2),y2(:,2),y12(:,1),xa,norm)
      
      dertwospline=derspline(x2(1,:),w(:),w1(:),xb,norm)

      end function



      double precision function derspline(x,y,yprime,xval,norm)
      implicit none
      integer i
      logical norm
      double precision x(2),y(2),yprime(2),c(4),yvec(4)
      double precision dx,xval,t
      double precision a(4,4)
      data a/1.d0,0.d0,-3.d0,2.d0,0.d0,0.d0,3.d0,-2.d0,0.d0,1.d0,
     &-2.d0,1.d0,0.d0,0.d0,-1.d0,1.d0/



      dx=x(2)-x(1)
      do i=1,2
          yvec(i)=y(i)
          yvec(i+2)=dx*yprime(i)
      end do

      c=matmul(transpose(a),yvec) 

      t=(xval-x(1))/dx

      derspline=0.d0
      do i=1,4
          derspline=derspline+c(i)*t**(i-1)
      end do

      end function
      


      DOUBLE PRECISION FUNCTION INTERPOLATEN0(T)
      IMPLICIT NONE
      COMMON/logfile/logfid
      INTEGER logfid
C--number density parameters
      common/n0par/ densconst, n0array(2000), tempn0array(2000)
      double precision densconst, n0array, tempn0array
C--local variables
      double precision T, Tmin, Tmax, n0min, n0max
      integer counter

      interpolaten0 = 0.d0

      ! Flag if out of bounds
      if ((T .lt. minval(tempn0array)) .or.
     &(T .gt. maxval(tempn0array))) then
          write(logfid,*) 
     &"Temperature out of CUSTOMN0F bounds in INTERPOLATEN0"
          write(logfid,*) "Will continue with n0 = 0 (no density)"
          return
      end if

      ! Find Tmin and Tmax
      do counter = 1, 1800
          if ((T .gt. tempn0array(counter)) .and.
     &(T .le. tempn0array(counter + 1))) then
          
              Tmin = tempn0array(counter)
              Tmax = tempn0array(counter + 1)
              n0min = n0array(counter)
              n0max = n0array(counter + 1)

              ! write(*, *) "Tmin = ", Tmin 
              ! write(*, *) "Tmax = ", Tmax 

              ! write(*, *) "n0 = ",  
      ! &n0min + (T - Tmin) * (n0max - n0min) / (Tmax - Tmin)
              ! write(*, *) "n0min = ", n0min 
              ! write(*, *) "n0max = ", n0max 

              ! write(*, *) "counter = ", counter

              ! Linear interpolation
              interpolaten0 = n0min + (T - Tmin) * (n0max - n0min) 
     &/ (Tmax - Tmin)
              return
         end if 
      end do 
          
      write(logfid,*) "INTERPOLATEN0: T bounded but no inner value"
      write(logfid,*) "This should not happen, check function"
      write(logfid,*) "Will continue with n0 = 0 (no density)"
      
      return
      
      END
