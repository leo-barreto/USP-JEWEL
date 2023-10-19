C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C++                            USP-JEWEL                            ++
C++                                                                 ++
C++ The program is part of the developted interface between the     ++
C++ parton propagation of JEWEL and an external hydrodynamic 2+1D   ++
C++ medium profile, intended for v-USPhydro.                        ++
C++                                                                 ++
C++ This code implements an interpolation function estimate a       ++
C++ value based on input coordinate (X4,Y4,T4) and a grid tprofile. ++
C++                                                                 ++
C++                                                                 ++
C++ Created by:                                                     ++
C++  - Fabio M. Canedo [fabio.canedo@usp.br]                        ++
C++  Instituto de Fisica, Universidade de Sao Paulo, Brazil         ++
C++  2019                                                           ++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DOUBLE PRECISION FUNCTION INTERPOLATE(X4,Y4,T4)
      IMPLICIT NONE
      INTEGER NX,NY,NT,NF
      DOUBLE PRECISION DX,DT,XMAX,XMIN,TMAX
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU

      common/grid/timesteps(60),tprofile(834,834,60)
      double precision timesteps,tprofile
      double precision tempx(2,2),tempy(2,2),tempxy(2,2)
      double precision p(4,4),a(4,4),au(4,4),ptemp
      double precision vu(4),vt(4)
      DOUBLE PRECISION X4,Y4,T4
      INTEGER I,J,N,K,ii,jj,iii,jjj,bigi
      DOUBLE PRECISION T,U,AX,BX,AY,BY,STEP,FRA
      DOUBLE PRECISION ITEMP,FTEMP
    
      nx=834
      ny=834
      nt=60
      xmax=25.d0
      xmin=-25.d0
      dx=(xmax-xmin)/(834.d0-1.d0)

      ITEMP=0.0d0
      FTEMP=0.0d0
      N=INT(MAX((XMAX-XMIN+dx)/dx,0.0d0))
      AX=(X4-xmin)/(xmax-xmin)
      AY=(Y4-xmin)/(xmax-xmin)
      BX=AX*(N-1)+1
      BY=AY*(N-1)+1
      I=FLOOR(BX)
      J=FLOOR(BY)
      T=BX-I
      U=BY-J

      k=1
      do ii=1,60
      if(timesteps(ii).lt.t4) then
      k=ii
      fra=(t4-timesteps(ii))/(timesteps(ii+1)-timesteps(ii))
      end if
      end do


      IF(I>nx.OR.J>nx.OR.K>nt.OR.I<1.OR.J<1.OR.K<1) THEN
      INTERPOLATE=0.0d0
      RETURN
      ENDIF

      do bigi=0,1
      do ii=0,1
            do jj=0,1
                  IF(i>nx-3.OR.j>nx-3.OR.i<2.OR.j<2) THEN
                        tempx(ii,jj)=0.d0
                        tempy(ii,jj)=0.d0
                        tempxy(ii,jj)=0.d0
                  else
                        tempx(ii+1,jj+1)=(tprofile(i+ii+1,j+jj,k+bigi)- &
     &                  tprofile(i+ii-1,j+jj,k+bigi))/2.
                        tempy(ii+1,jj+1)=(tprofile(i+ii,j+jj+1,k+bigi)- &
     &                  tprofile(i+ii,j+jj-1,k))/2.
                        tempxy(ii+1,jj+1)=(tprofile(i+ii-1,j+jj-1,k+bigi&
     &                  )-tprofile(i+ii-1,j+jj+1,k+bigi)-               &
     &                  tprofile(i+ii+1,j+jj-1,k+bigi)+                 &
     &                  tprofile(i+ii+1,j+jj+1,k+bigi))/4.
                  endif                  
            enddo
      enddo
        
      do ii=0,1
            do jj=0,1
      iii=1
      jjj=1
      p(iii+ii,jjj+jj)=tprofile(i+ii,j+jj,k+bigi)

      iii=3
      jjj=1
      p(ii+iii,jj+jjj)=tempx(1+ii,1+jj)
    
      iii=1
      jjj=3
      p(ii+iii,jj+jjj)=tempy(1+ii,1+jj)
    
      iii=3
      jjj=3
      p(ii+iii,jj+jjj)=tempxy(1+ii,1+jj)
            enddo
      enddo

      au(1,1)=1.d0
      au(1,2)=0.d0
      au(1,3)=0.d0
      au(1,4)=0.d0
      au(2,1)=0.d0
      au(2,2)=0.d0
      au(2,3)=1.d0
      au(2,4)=0.d0
      au(3,1)=-3.d0
      au(3,2)=3.d0
      au(3,3)=-2.d0
      au(3,4)=-1.d0
      au(4,1)=2.d0 
      au(4,2)=-2.d0
      au(4,3)=1.d0
      au(4,4)=1.d0
       
      a=au*p*transpose(au)
       
      do ii=1,4
            vu(ii)=u**(ii-1)
            vt(ii)=t**(ii-1)
      enddo     
      
      ptemp=0.0d0
      do ii=1,4
            do jj=1,4
                  ptemp=ptemp+vu(ii)*a(ii,jj)*vt(jj)
            enddo
      enddo
      
      if(bigi.eq.0d0) then
      itemp = ptemp
      else
      ftemp = ptemp
      endif
      
      enddo
 
C     ITEMP=(1.0-T)*(1.0-U)*TEMP(I,J,K)+(T)*(1.0-U)*TEMP(I+1,J,K)
C    &+(1.0-T)*(U)*TEMP(I,J+1,K)+(T)*(U)*TEMP(I+1,J+1,K)


C     FTEMP=(1-T)*(1-U)*TEMP(I,J,K+1)+(T)*(1-U)*TEMP(I+1,J,K+1)
C    &+(1-T)*(U)*TEMP(I,J+1,K+1)+(T)*(U)*TEMP(I+1,J+1,K+1)
      !write(*,*) "itemp=",itemp
      !write(*,*) "ftemp=",ftemp
      !write(*,*) "fra=",fra
      INTERPOLATE=ITEMP+(FTEMP-ITEMP)*FRA
      if(interpolate.lt.0.0d0) interpolate=0.0d0
      !write(*,*) "interpolate=",interpolate
      

      IF(INTERPOLATE.GT.1) THEN
      !WRITE(*,*) "Interpolate: ",FTEMP, ITEMP, x4, y4, t4
      !write(*,*) u, t, fra
      ENDIF

      END
