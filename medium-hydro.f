C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C++ Copyright (C) 2017 Korinna C. Zapp [Korinna.Zapp@cern.ch]       ++
C++                                                                 ++
C++ This file is part of JEWEL 2.2.0                                ++
C++                                                                 ++
C++ The JEWEL homepage is jewel.hepforge.org                        ++
C++                                                                 ++
C++ The medium model was partly implemented by Jochen Klein.        ++
C++ Raghav Kunnawalkam Elayavalli helped with the implementation    ++
C++ of the V+jet processes.                                         ++
C++                                                                 ++
C++ Please follow the MCnet GUIDELINES and cite Eur.Phys.J. C74     ++
C++ (2014) no.2, 2762 [arXiv:1311.0048] for the code and            ++
C++ JHEP 1303 (2013) 080 [arXiv:1212.1599] and                      ++
C++ optionally EPJC 60 (2009) 617 [arXiv:0804.3568] for the         ++
C++ physics. The reference for V+jet processes is EPJC 76 (2016)    ++
C++ no.12 695 [arXiv:1608.03099] and for recoil effects it is       ++
C++ arXiv:1707.01539.
C++                                                                 ++
C++ JEWEL relies heavily on PYTHIA 6 for the event generation. The  ++
C++ modified version of PYTHIA 6.4.25 that is distributed with      ++
C++ JEWEL is, however, not an official PYTHIA release and must not  ++
C++ be used for anything else. Please refer to results as           ++
C++ "JEWEL+PYTHIA".                                                 ++
C++                                                                 ++
C++ JEWEL also uses code provided by S. Zhang and J. M. Jing        ++
C++ (Computation of Special Functions, John Wiley & Sons, New York, ++
C++ 1996 and http://jin.ece.illinois.edu) for computing the         ++
C++ exponential integral Ei(x).                                     ++
C++                                                                 ++
C++                                                                 ++
C++ JEWEL  is free software; you can redistribute it and/or         ++
C++ modify it under the terms of the GNU General Public License     ++
C++ as published by the Free Software Foundation; either version 2  ++
C++ of the License, or (at your option) any later version.          ++
C++                                                                 ++
C++ JEWEL is distributed in the hope that it will be useful,        ++
C++ but WITHOUT ANY WARRANTY; without even the implied warranty of  ++
C++ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the    ++
C++ GNU General Public License for more details.                    ++
C++                                                                 ++
C++ You should have received a copy of the GNU General Public       ++  
C++ License along with this program; if not, write to the Free      ++
C++ Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, ++
C++ MA 02110-1301 USA                                               ++
C++                                                                 ++
C++ Linking JEWEL statically or dynamically with other modules is   ++
C++ making a combined work based on JEWEL. Thus, the terms and      ++
C++ conditions of the GNU General Public License cover the whole    ++
C++ combination.                                                    ++
C++                                                                 ++
C++ In addition, as a special exception, I give you permission to   ++
C++ combine JEWEL with the code for the computation of special      ++
C++ functions provided by S. Zhang and J. M. Jing. You may copy and ++
C++ distribute such a system following the terms of the GNU GPL for ++
C++ JEWEL and the licenses of the other code concerned, provided    ++
C++ that you include the source code of that other code when and as ++
C++ the GNU GPL requires distribution of source code.               ++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE MEDINIT(FILE,id,etam,mass)
      IMPLICIT NONE
C--medium parameters
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     & NX,NY,NT,NF,DX,DT,XMAX,XMIN,TMAX
      INTEGER NX,NY,NT,NF
      DOUBLE PRECISION DX,DT,XMAX,XMIN,TMAX
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      common/grid/timesteps(60),tprofile(834,834,60),vtxmap(834,834)
      double precision timesteps,tprofile,vtxmap
      common/gridvel/ux(834,834,60),uy(834,834,60)
      double precision ux,uy
      COMMON/MEDFILEC/MEDFILE,NLIST,endoff
      CHARACTER*200 MEDFILE
      INTEGER NLIST
      logical endoff
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,MODMED,MEDFILELIST
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON,MODMED,MEDFILELIST
C--max rapidity
      common/rapmax2/etamax2
      double precision etamax2
C--longitudinal boost of momentum distribution
      common/boostmed/boost
      logical boost
C--factor to vary Debye mass
      COMMON/MDFAC/MDFACTOR,MDSCALEFAC
      DOUBLE PRECISION MDFACTOR,MDSCALEFAC
      common/temperature/tempfac
      double precision tempfac
C--nuclear thickness function
      COMMON /THICKFNC/ RMAX,TA(100,2)
      DOUBLE PRECISION RMAX,TA
C--geometrical cross section
      COMMON /CROSSSEC/ IMPMAX,CROSS(200,3)
      DOUBLE PRECISION IMPMAX,CROSS
C--hydrodynamic quantities
      COMMON /HYDROLIM/ MIDRAPLIM, TMAXLIM, TVELMAXLIM, BOOSTTR, 
     &GLOBALLIMS, PRETAUHYDRO
      DOUBLE PRECISION MIDRAPLIM, TMAXLIM, TVELMAXLIM
      LOGICAL BOOSTTR, GLOBALLIMS, PRETAUHYDRO
C--hydro auxiliary files
      COMMON /HYDROF/ INITVTXF
      CHARACTER*200 INITVTXF
C--identifier of log file
      common/logfile/logfid
      integer logfid
C--grid parameters
      common/gridpar/ gdt,gdx,gxmax,gxmin,gnx,gny,gnt
      double precision gdt,gdx,gxmax,gxmin
      integer gnx,gny,gnt

      DATA RAU/10./
      DATA D3/0.9d0/
      DATA ZETA3/1.2d0/
C--local variables
      INTEGER I,LUN,POS,IOS,id,mass
      double precision etam
      CHARACTER*100 BUFFER,LABEL,tempbuf
      CHARACTER*100 FILE
      character firstchar
      logical fileexist

      ! Use MIDRAPLIM for limits and
      ! etamax2 for the simulation
      etamax2 = etam
      logfid = id

      IOS=0
      LUN=77

      NLIST=1
      endoff=.false.

C--default settings
      TAUI=0.6d0          ! Not used for hydro
      TI=0.36d0           ! Not used for hydro
      TC=0.17d0
      WOODSSAXON=.TRUE.   ! Not used for hydro
      CENTRMIN=0.d0       ! Not used for hydro
      CENTRMAX=10.d0      ! Not used for hydro
      NF=3
      A=mass
      N0=0.17d0
      D=0.54d0
      SIGMANN=6.2
      MDFACTOR=0.45d0
      MDSCALEFAC=0.9d0
      tempfac=1.0d0
      boost =.true.
      breal=0.d0
C--hydro settings
      MODMED=.TRUE.       ! Deprecated
      MEDFILELIST=.FALSE. ! Deprecated
      GLOBALLIMS = .true. ! Global limits vs JEWEL default impl.
      ! Extrapolation of soft mid rap lim, +0.1 then tuning 
      ! (see nucl-ex/1612.08966)
      MIDRAPLIM=3.3d0
      TMAXLIM=0.675d0   ! Approx vUSP+TRENTo temp lim (PbPb 5TeV)
      BOOSTTR= .true. ! Logical for transverse velocity (F => u=0)
      TVELMAXLIM=0.927d0   ! Approx vUSP+TRENTo uT lim (PbPb 5TeV)
      ! Temperature assumption before TAU0, i.e. F => T(tau < tau0) = 0
      PRETAUHYDRO=.false.
      !GRIDN=834     ! Number of points in grid (per dimension)
      ! Initial vertex map file
      INITVTXF='/sampa/leonardo/USP-JEWEL/initvertexmap.dat' 


C--read settings from file
	write(logfid,*)
	inquire(file=FILE,exist=fileexist)
	if(fileexist)then
        write(logfid,*)'Reading medium parameters from ',FILE
        OPEN(unit=LUN,file=FILE,status='old',err=10)
	  do 20 i=1,1000
          READ(LUN, '(A)', iostat=ios) BUFFER
	     if (ios.ne.0) goto 30
	     firstchar = buffer(1:1)
	     if (firstchar.eq.'#') goto 20
          POS=SCAN(BUFFER,' ')
          LABEL=BUFFER(1:POS)
          BUFFER=BUFFER(POS+1:)
          IF (LABEL=="TAUI")THEN
            READ(BUFFER,*,IOSTAT=IOS) TAUI
          ELSE IF (LABEL=="TI") THEN
            READ(BUFFER,*,IOSTAT=IOS) TI
          ELSE IF (LABEL=="TC") THEN
            READ(BUFFER,*,IOSTAT=IOS) TC
          ELSE IF (LABEL=="WOODSSAXON") THEN
            READ(BUFFER,*,IOSTAT=IOS) WOODSSAXON
          ELSE IF (LABEL=="MODMED") THEN
            READ(BUFFER,*,IOSTAT=IOS) MODMED
          ELSE IF (LABEL=="MEDFILE") THEN
            READ(BUFFER,'(50A)',IOSTAT=IOS) MEDFILE
          ELSE IF (LABEL=="CENTRMIN") THEN
            READ(BUFFER,*,IOSTAT=IOS) CENTRMIN
          ELSE IF (LABEL=="CENTRMAX") THEN
            READ(BUFFER,*,IOSTAT=IOS) CENTRMAX
          ELSE IF (LABEL=="NF") THEN
            READ(BUFFER,*,IOSTAT=IOS) NF
          ELSE IF (LABEL=="N0") THEN
            READ(BUFFER,*,IOSTAT=IOS) N0
          ELSE IF (LABEL=="D") THEN
            READ(BUFFER,*,IOSTAT=IOS) D
          ELSE IF (LABEL=="SIGMANN") THEN
            READ(BUFFER,*,IOSTAT=IOS) SIGMANN
          ELSE IF (LABEL=="MDFACTOR") THEN
            READ(BUFFER,*,IOSTAT=IOS) MDFACTOR
          ELSE IF (LABEL=="MDSCALEFAC") THEN
            READ(BUFFER,*,IOSTAT=IOS) MDSCALEFAC
          ELSE IF (LABEL=="GLOBALLIMS") THEN
            READ(BUFFER,*,IOSTAT=IOS) GLOBALLIMS
          ELSE IF (LABEL=="MIDRAPLIM") THEN
            READ(BUFFER,*,IOSTAT=IOS) MIDRAPLIM
          ELSE IF (LABEL=="TMAXLIM") THEN
            READ(BUFFER,*,IOSTAT=IOS) TMAXLIM
          ELSE IF (LABEL=="BOOSTTR") THEN
            READ(BUFFER,*,IOSTAT=IOS) BOOSTTR
          ELSE IF (LABEL=="TVELMAXLIM") THEN
            READ(BUFFER,*,IOSTAT=IOS) TVELMAXLIM
          ELSE IF (LABEL=="PRETAUHYDRO") THEN
            READ(BUFFER,*,IOSTAT=IOS) PRETAUHYDRO
          !ELSE IF (LABEL=="GRIDN") THEN
            !READ(BUFFER,*,IOSTAT=IOS) GRIDN
          ELSE IF (LABEL=="INITVTXF") THEN
            READ(BUFFER,*,IOSTAT=IOS) INITVTXF
          ! JEWEL original boost (longitudinal) as user option
          ELSE IF (LABEL=="BOOSTZ") THEN
            READ(BUFFER,*,IOSTAT=IOS) boost
          ELSE IF (LABEL=="BREAL") THEN
            READ(BUFFER,*,IOSTAT=IOS) breal

	  else
       write(logfid,*)'unknown label ',label
	     endif
 20	  continue

 30	  close(LUN,status='keep')
	  write(logfid,*)'...done'
	  goto 40

 10     write(logfid,*)'Could not open medium parameter file, '//
     &	'will run with default settings.'

	else
	  write(logfid,*)'No medium parameter file found, '//
     &	'will run with default settings.'
	endif

 40   write(logfid,*)'using parameters:'
      write(logfid,*)'TAUI        = ',TAUI
      write(logfid,*)'TI          = ',TI
      write(logfid,*)'TC          = ',TC
      write(logfid,*)'WOODSSAXON  = ',WOODSSAXON
      write(logfid,*)'MODMED      = ',MODMED
      write(logfid,*)'MEDFILELIST = ',MEDFILELIST
      write(logfid,*)'CENTRMIN    = ',CENTRMIN
      write(logfid,*)'CENTRMAX    = ',CENTRMAX
      write(logfid,*)'NF          = ',NF
      write(logfid,*)'A           = ',A
      write(logfid,*)'N0          = ',N0
      write(logfid,*)'D           = ',D
      write(logfid,*)'SIGMANN     = ',SIGMANN
      write(logfid,*)'MDFACTOR    = ',MDFACTOR
      write(logfid,*)'MDSCALEFAC  = ',MDSCALEFAC
      write(logfid,*)'BREAL       = ',breal
      write(logfid,*)'GLOBALLIMS  = ',GLOBALLIMS
      write(logfid,*)'MIDRAPLIM   = ',MIDRAPLIM
      write(logfid,*)'TMAXLIM     = ',TMAXLIM
      write(logfid,*)'BOOSTTR    = ',BOOSTTR
      write(logfid,*)'TVELMAXLIM  = ',TVELMAXLIM
      write(logfid,*)'PRETAUHYDRO = ',PRETAUHYDRO
      !write(logfid,*)'GRIDN       = ',GRIDN
      write(logfid,*)'INITVTXF    = ',INITVTXF
      write(logfid,*)'BOOSTZ      = ',boost
      write(logfid,*)
      write(logfid,*)
      write(logfid,*)


      if (.not. boost) then
        write(logfid,*) 'No longitudinal boost => MIDRAPLIM = 0'
        MIDRAPLIM = 0.0
      else if (MIDRAPLIM.lt.etamax2) then
        MIDRAPLIM = etamax2
        write(logfid,*) 'ETAMAX > MIDRAPLIM'
        write(logfid,*) 'Extrapolating mid-rapidity limit' 
        write(logfid,*)
      end if

      if (.not. BOOSTTR) then
        TVELMAXLIM = 0.d0 
        write(logfid,*) 'No transverse u => TVELMAXLIM = 0'
        write(logfid,*)
      end if

C--Call the modified medium setup
      CALL MYMED()
      write(logfid,*) 'Hydrodynamic profile loaded: ', MEDFILE
      write(logfid,*)
     
      END

      SUBROUTINE MEDNEXTEVT
      IMPLICIT NONE
C--medium parameters
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     & NX,NY,NT,NF,DX,DT,XMAX,XMIN,TMAX
      INTEGER NX,NY,NT,NF
      DOUBLE PRECISION DX,DT,XMAX,XMIN,TMAX
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      common/grid/timesteps(60),tprofile(834,834,60),vtxmap(834,834)
      double precision timesteps,tprofile,vtxmap
      common/gridvel/ux(834,834,60),uy(834,834,60)
      double precision ux,uy
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,MODMED
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON,MODMED
C--geometrical cross section
      COMMON /CROSSSEC/ IMPMAX,CROSS(200,3)
      DOUBLE PRECISION IMPMAX,CROSS
C--local variables
      integer i,j
      DOUBLE PRECISION PYR,R,b1,b2,gettemp

      ! Dummy function
      END



      SUBROUTINE PICKVTX(x, y)
      IMPLICIT NONE
      DOUBLE PRECISION x, y
C--medium parameters
      common/grid/timesteps(60), tprofile(834,834,60),vtxmap(834,834)
      double precision timesteps,tprofile,vtxmap
      integer tries, maxpos, minpos, irand, jrand, counter
      double precision pyr, zval
      common/logfile/logfid
      integer logfid

      do counter = 1, 5
          zval = pyr(0) * maxval(vtxmap)

          ! Consider only a 10 fm by 10 fm (k = 250, 583)
          ! since anywhere else won't have entropy (PbPb)
          maxpos = 583
          minpos = 250

          do tries = 1, 1000000
              irand = int(pyr(0) * (maxpos - minpos) + minpos)
              jrand = int(pyr(0) * (maxpos - minpos) + minpos)

              if (vtxmap(irand, jrand) .gt. zval) then
                  x = -25.d0 + (irand - 1) * (50.d0 / 833.d0)
                  y = -25.d0 + (jrand - 1) * (50.d0 / 833.d0)
                  return
              end if
          end do

          write(*,*) 'Failed to find vtx, restarting selection.'
          write(*,*) 'This should not happen often.'
      end do

      ! Only do process 5 times, kill simulation otherwise
      write(*,*) 'No initial vertex found. Check initialvtx table.'
      write(logfid,*) 'No initial vertex found. Check initialvtx table.'
      stop

      END SUBROUTINE



      SUBROUTINE GETSCATTERER(X,Y,Z,T,TYPE,PX,PY,PZ,E,MS)
      IMPLICIT NONE
C--medium parameters
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     & NX,NY,NT,NF,DX,DT,XMAX,XMIN,TMAX
      INTEGER NX,NY,NT,NF
      DOUBLE PRECISION DX,DT,XMAX,XMIN,TMAX
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      common/grid/timesteps(60),tprofile(834,834,60),vtxmap(834,834)
      double precision timesteps,tprofile,vtxmap
C--internal medium parameters
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,MODMED
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON,MODMED
C--longitudinal boost of momentum distribution
      common/boostmed/boost
      logical boost
C--function calls
      DOUBLE PRECISION GETTEMP,GETMD,GETMOM,GETMS
C--identifier of log file
      common/logfile/logfid
      integer logfid
C--local variables
      DOUBLE PRECISION X,Y,Z,T,MS,PX,PY,PZ,E,MD,TEMP
      double precision u,ux,uy,px2,py2,px3,py3,e3
      double precision getu,getutheta,eta,vel
      INTEGER TYPE
      DOUBLE PRECISION R,PYR,pmax,wt,tau,theta,phi,pi,p,ys,pz2,e2
      double precision e4
      DATA PI/3.141592653589793d0/

      R=PYR(0)
      IF(R.LT.(2.*12.*NF*D3/3.)/(2.*12.*NF*D3/3.+3.*16.*ZETA3/2.))THEN
         TYPE=2
      ELSE
         TYPE=21
      ENDIF

      ! JEWEL original
      !MS=GETMS(X,Y,Z,T)
      !MD=GETMD(X,Y,Z,T)
      
      MD=GETMD(X,Y,Z,T)
      MS = MD / sqrt(2.)

      TEMP=GETTEMP(X,Y,Z,T)
      tau=sqrt(t**2-z**2)
      if (boost) then
        ys = 0.5d0*log((t+z)/(t-z))
      else
        ys = 0.d0
      endif

      pmax = 10.*temp

      IF(TEMP.LT.1.D-2) THEN
        write(logfid,*)'asking for a scattering centre without medium:'
        write(logfid,*)'at (x,y,z,t)=',X,Y,Z,T
        write(logfid,*)'making one up to continue but '//
     &'something is wrong!'
        TYPE=21
        PX=0.d0
        PY=0.d0
        PZ=0.d0

       ! JEWEL original
       !MS=GETMS(0.d0,0.d0,0.d0,0.d0)
       !MD=GETMD(0.d0,0.d0,0.d0,0.d0)

        MD=GETMD(0.d0,0.d0,0.d0,0.d0)
        MS=MD / sqrt(2.)
        E=SQRT(PX**2+PY**2+PZ**2+MS**2)
        RETURN
      ENDIF

 10	p = pyr(0)**0.3333333*pmax
	E2 = sqrt(p**2+ms**2)
	if (type.eq.2) then
	  wt = (exp(ms/temp)-1.)/(exp(E2/temp)-1.)
	else
	  wt = (exp(ms/temp)+1.)/(exp(E2/temp)+1.)
	endif
	if (wt.gt.1.) write(logfid,*)'Error in getscatterer: weight = ',wt
	if (wt.lt.0.) write(logfid,*)'Error in getscatterer: weight = ',wt
	if (pyr(0).gt.wt) goto 10
	phi = pyr(0)*2.*pi
	theta = -acos(2.*pyr(0)-1.)+pi
	px  = p*sin(theta)*cos(phi)
	py  = p*sin(theta)*sin(phi)
	pz2 = p*cos(theta)

      !JEWEL original longitudinal boost
      !E   = cosh(ys)*E2 + sinh(ys)*pz2
      !pz  = sinh(ys)*E2 + cosh(ys)*pz2
 
      !Perform boost
      E = E2
      pz = pz2
      call LorentzLocalBoost(E,px,py,pz,x,y,z,t)

      IF (E.lt.0.d0) THEN
        write(logfid,*) 'Negative energy (', E, ') in GETSCATTERER'
      END IF
      END




      SUBROUTINE AVSCATCEN(X,Y,Z,T,PX,PY,PZ,E,m)
      IMPLICIT NONE
C--longitudinal boost of momentum distribution
       common/boostmed/boost
       logical boost
C--max rapidity
      common/rapmax2/etamax2
      double precision etamax2
      double precision getu,getutheta,theta,u,gettemp
C--local variables
      double precision x,y,z,t,px,py,pz,e,getms,m,ys,temp

      !Original JEWEL implementation 
      !if (boost) then
      !  ys = 0.5*log((t+z)/(t-z))
      !  if ((z.eq.0.d0).and.(t.eq.0.d0)) ys =0.d0
      !  if (ys.gt.etamax2) ys=etamax2
      !  if (ys.lt.-etamax2) ys=-etamax2

      !else
      !  ys = 0.d0
      !endif
     
      !m  = getms(x,y,z,t)
      !e = m*cosh(ys)
      !px = 0.d0
      !py = 0.d0
      !pz = m*sinh(ys)
     
      
      m  = getms(x,y,z,t)
      e = m
      px = 0.d0
      py = 0.d0
      pz = 0.d0
     
      !Apply boost
      call LorentzLocalBoost(e,px,py,pz,x,y,z,t)

      end


      SUBROUTINE maxscatcen(PX,PY,PZ,E,m)
      IMPLICIT NONE
C--longitudinal boost of momentum distribution
      common/boostmed/boost
      logical boost
C--max rapidity
      common/rapmax2/etamax2
      double precision etamax2
      double precision getu,getutheta,eta,u,theta
      COMMON/GAMMAMAX/GAMMAMAXIMUM,VELMAXIMUM,DENSITYMAXIMUM,
     &DENSITYMINIMUMDENSITYMINIMUM, TAUMIN
      DOUBLE PRECISION GAMMAMAXIMUM,VELMAXIMUM,DENSITYMAXIMUM,
     &DENSITYMINIMUMDENSITYMINIMUM, TAUMIN
      double precision pi
      DATA PI/3.141592653589793d0/
      double precision gettempmax
C--hydrodynamic limits
      COMMON /HYDROLIM/ MIDRAPLIM, TMAXLIM, TVELMAXLIM, BOOSTTR, 
     &GLOBALLIMS, PRETAUHYDRO
      DOUBLE PRECISION MIDRAPLIM, TMAXLIM, TVELMAXLIM
      LOGICAL BOOSTTR, GLOBALLIMS, PRETAUHYDRO
C--local variables
      double precision px,py,pz,e,getmsmax,m,ys
      
      !Original JEWEL implementation
      !if (boost) then
      !  ys = etamax2
      !else
      !  ys = 0.d0
      !endif

      ! GETSSCAT limits are boost invariant
      ! thus we keep the original code without
      ! a transverse boost of the 4-momentum
      ! This must be checked.

      if (boost) then
        ys = MIDRAPLIM
      else
        ys = 0.d0
      endif
      
      m = getmsmax()
      e = m * cosh(ys)
      px = 0
      py = 0
      pz = m * sinh(ys)
      end


      DOUBLE PRECISION FUNCTION GETMD(X1,Y1,Z1,T1)
      IMPLICIT NONE
C--factor to vary Debye mass
      COMMON/MDFAC/MDFACTOR,MDSCALEFAC
      DOUBLE PRECISION MDFACTOR,MDSCALEFAC
      DOUBLE PRECISION X1,Y1,Z1,T1,GETTEMP
      double precision getmdmin
      GETMD=MDSCALEFAC*3.*GETTEMP(X1,Y1,Z1,T1)
      !GETMD=MAX(GETMD,MDFACTOR)
      GETMD=MAX(GETMD,GETMDMIN())
      END



      DOUBLE PRECISION FUNCTION GETMS(X2,Y2,Z2,T2)
      IMPLICIT NONE
      DOUBLE PRECISION X2,Y2,Z2,T2,GETMD
      ! Deprecated
      GETMS=GETMD(X2,Y2,Z2,T2)/SQRT(2.)
      END



      DOUBLE PRECISION FUNCTION GETNEFF(X3,Y3,Z3,T3,P0,P1,P2,P3)
      IMPLICIT NONE
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     & NX,NY,NT,NF,DX,DT,XMAX,XMIN,TMAX
      INTEGER NX,NY,NT,NF
      DOUBLE PRECISION DX,DT,XMAX,XMIN,TMAX
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      common/grid/timesteps(60),tprofile(834,834,60),vtxmap(834,834)
      double precision timesteps,tprofile,vtxmap
      common/gridvel/ux(834,834,60),uy(834,834,60)
      double precision ux,uy
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,MODMED
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON,MODMED
      common/temperature/tempfac
      double precision tempfac
      common/boostmed/boost
      logical boost
C--max rapidity
      common/rapmax2/etamax2
      double precision etamax2
C--   local variables
      DOUBLE PRECISION X3,Y3,Z3,T3,PI,GETTEMP,tau,cosheta
      double precision getu,getutheta
      double precision umx, umy, umz, umr, umtheta, gam, vp, gamp, 
     &localtemp
      double precision J0, J1, J2, J3, P0, P1, P2, P3, ys
      DATA PI/3.141592653589793d0/

      getneff=0.d0
      localtemp = GETTEMP(X3,Y3,Z3,T3)
      IF ((ABS(Z3).gt.T3) .OR. (localtemp.le.TC)) THEN 
        RETURN
      END IF

      tau = sqrt(t3**2-z3**2)
      if (boost) then 
        umz = z3 / tau
      else 
        umz = 0.d0 
      end if
      
      !Medium 4-velocity
      umr = getu(x3, y3, z3, t3, localtemp)
      umtheta = getutheta(x3, y3, z3, t3, localtemp) 
      umx = umr * cos(umtheta)
      umy = umr * sin(umtheta)
      gam = sqrt(1 + umr ** 2 + umz ** 2)
      

      !cosheta = t3/tau
      GETNEFF=(2.*6.*NF*D3*2./3. + 16.*ZETA3*3./2.)
     &     *localtemp**3/PI**2
      !write(*,*) "Original neff = ", GETNEFF/cosheta

      !Scattering center 4-current
      J0 = getneff * gam
      J1 = getneff * umx
      J2 = getneff * umy
      J3 = getneff * umz

      !Effective density transform as 
      !J0 = p_mu J^mu / p0, check nucl-th/0612068
      getneff = (P0 * J0 - P1 * J1 - P2 * J2 - P3 * J3) / P0

      END
      
      

      DOUBLE PRECISION FUNCTION GETTEMP(X4,Y4,Z4,T4)
      IMPLICIT NONE
C--medium parameters
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     & NX,NY,NT,NF,DX,DT,XMAX,XMIN,TMAX
      INTEGER NX,NY,NT,NF
      DOUBLE PRECISION DX,DT,XMAX,XMIN,TMAX
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      common/grid/timesteps(60),tprofile(834,834,60),vtxmap(834,834)
      double precision timesteps,tprofile,vtxmap
      common/gridvel/u(834,834,60),utheta(834,834,60)
      double precision u,utheta
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,MODMED
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON,MODMED
      COMMON/GAMMAMAX/GAMMAMAXIMUM,VELMAXIMUM,DENSITYMAXIMUM,
     &DENSITYMINIMUMDENSITYMINIMUM, TAUMIN
      DOUBLE PRECISION GAMMAMAXIMUM,VELMAXIMUM,DENSITYMAXIMUM,
     &DENSITYMINIMUMDENSITYMINIMUM, TAUMIN
C--max rapidity
      common/rapmax2/etamax2
      double precision etamax2
      common/temperature/tempfac
      double precision tempfac
C--hydrodynamic quantities
            COMMON /HYDROLIM/ MIDRAPLIM, TMAXLIM, TVELMAXLIM, BOOSTTR, 
     &      GLOBALLIMS, PRETAUHYDRO
            DOUBLE PRECISION MIDRAPLIM, TMAXLIM, TVELMAXLIM
            LOGICAL BOOSTTR, GLOBALLIMS, PRETAUHYDRO
C--local variables
      DOUBLE PRECISION X4,Y4,Z4,T4,TAU,NPART,EPS0,EPSIN,TEMPIN,PI,
     &NTHICK,ys,MEDPART,interpolate
      double precision gettempmax
      DATA PI/3.141592653589793d0/

      GETTEMP=0.D0

      IF (ABS(Z4).gt.T4) RETURN

      TAU=SQRT(T4**2-Z4**2)
      ! NO ASSUMPTION BEFORE TAU0
      if ((tau.lt.TAUMIN-0.001) .and. (.not.PRETAUHYDRO)) return

      ! Consider only relevant regions for calculation
      ! otherwise temp = 0 => prob of interaction = 0
      ys = 0.5d0 * log((T4 + Z4) / (T4 - Z4))
      if (ys.gt.etamax2) return

      ! Grid values are between -25 and 25
      if ((abs(X4).gt.25.d0) .or. (abs(Y4).gt.25.d0)) return


      GETTEMP=tempfac*interpolate(X4,Y4,tau,1)
      if(gettemp.lt.tc) gettemp=0.0d0
      if(gettemp.ge.gettempmax()) gettemp=gettempmax()

      RETURN

      END

      double precision function getu(x,y,z,t,localtemperature)
            implicit none
            integer np
            double precision x,y,z,t,tau,localtemperature
            double precision interpol 
            common/grid/timesteps(60),tprofile(834,834,60),
     &      vtxmap(834,834)
            double precision timesteps,tprofile,vtxmap
            common/gridvel/u(834,834,60),utheta(834,834,60)
            double precision u,utheta           
            COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &      N0,SIGMANN,A,WOODSSAXON,MODMED,MEDFILELIST
            DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &      SIGMANN
            INTEGER A
            LOGICAL WOODSSAXON,MODMED,MEDFILELIST
C--hydrodynamic quantities
            COMMON /HYDROLIM/ MIDRAPLIM, TMAXLIM, TVELMAXLIM, BOOSTTR, 
     &      GLOBALLIMS, PRETAUHYDRO
            DOUBLE PRECISION MIDRAPLIM, TMAXLIM, TVELMAXLIM
            LOGICAL BOOSTTR, GLOBALLIMS, PRETAUHYDRO

            getu = 0.d0
            IF (BOOSTTR .eqv. .false.) then
              RETURN
            END IF
            tau = sqrt(t**2 - z**2)
            IF ((tau.le.TAUI) .OR. (localtemperature.le.TC)) THEN
              RETURN
            END IF
            getu=interpol(tau,x,y,834,timesteps,u,.true.)
            return
      end function

      double precision function getutheta(x,y,z,t,localtemperature)
            implicit none
            integer np
            double precision x,y,z,t,tau,localtemperature
            double precision interpol 
            common/grid/timesteps(60),tprofile(834,834,60),
     &      vtxmap(834,834)
            double precision timesteps,tprofile,vtxmap
            common/gridvel/u(834,834,60),utheta(834,834,60)
            double precision u,utheta           
            COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &      N0,SIGMANN,A,WOODSSAXON,MODMED,MEDFILELIST
            DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &      SIGMANN
            INTEGER A
            LOGICAL WOODSSAXON,MODMED,MEDFILELIST
C--hydrodynamic quantities
            COMMON /HYDROLIM/ MIDRAPLIM, TMAXLIM, TVELMAXLIM, BOOSTTR, 
     &      GLOBALLIMS, PRETAUHYDRO
            DOUBLE PRECISION MIDRAPLIM, TMAXLIM, TVELMAXLIM
            LOGICAL BOOSTTR, GLOBALLIMS, PRETAUHYDRO

            getutheta = 0.d0
            IF (BOOSTTR .eqv. .false.) then
              RETURN
            END IF

            tau = sqrt(t**2 - z**2)
            IF ((tau.le.TAUI) .OR. (localtemperature.le.TC)) THEN
              RETURN
            END IF
            getutheta=interpol(tau,x,y,834,timesteps,utheta,.false.)
            return
      end function

      
      subroutine LorentzLocalBoost(e, px, py, pz, x, y, z, t)
            implicit none
C--longitudinal boost of momentum distribution
            common/boostmed/boost
            logical boost
C--max rapidity
            common/rapmax2/etamax2
            double precision etamax2
            double precision fourmom(4,1), fourmomboost(4,1), 
     &      boostm(4,4)
            double precision ux, uy, uz, gam, u2, ux2, uy2, uz2, 
     &      unorm, uangle, tau
            double precision getu,getutheta
C--local variables
            double precision e, px, py, pz, x, y, z, t

            ! Use temp as 10.d0 to always calculate (T_C < 10)
            unorm = getu(x, y, z, t, 10.d0)
            uangle = getutheta(x, y, z, t, 10.d0)
            tau = sqrt(t ** 2 - z ** 2)

            ! Note that all 4-velocity components must be flipped
            ! since they are defined in the lab frame, thus the frame
            ! must be boosted as -u so the scattering centers are 
            ! boosted as u.
            if (tau.gt.0.d0 .and. boost) then
              uz = -z / tau
            else
              uz = 0.d0
            end if

            ! Respect simulation limit in uz
            if (abs(uz).gt.sinh(etamax2)) then
              if (z.gt.0d0) then
                uz = -sinh(etamax2)
              else 
                uz = sinh(etamax2)
              end if
            end if
            
            ! Only transform if there is velocity
            if (unorm.gt.0.d0 .or. abs(uz).gt.0.d0) then
              ux = -unorm * cos(uangle)
              uy = -unorm * sin(uangle)
              uz2 = uz ** 2
              u2 = unorm ** 2 + uz2
              ux2 = ux ** 2
              uy2 = uy ** 2
              gam = sqrt(1 + u2)

              ! Define initial four-momentum
              fourmom(1,1) = e
              fourmom(2,1) = px
              fourmom(3,1) = py
              fourmom(4,1) = pz

              ! Define boost matrix
              boostm(1,1) = gam
              boostm(1,2) = -ux
              boostm(1,3) = -uy
              boostm(1,4) = -uz

              boostm(2,1) = -ux
              boostm(2,2) = 1 + (gam - 1) * ux2 / u2
              boostm(2,3) = (gam - 1) * ux * uy / u2
              boostm(2,4) = (gam - 1) * ux * uz / u2

              boostm(3,1) = -uy
              boostm(3,2) = (gam - 1) * uy * ux / u2
              boostm(3,3) = 1 + (gam - 1) * uy2 / u2
              boostm(3,4) = (gam - 1) * uy * uz / u2
            
              boostm(4,1) = -uz
              boostm(4,2) = (gam - 1) * ux * uz / u2
              boostm(4,3) = (gam - 1) * uy * uz / u2
              boostm(4,4) = 1 + (gam - 1) * uz2 / u2

              ! Calculate new four-momentum
              fourmomboost = matmul(boostm, fourmom)
              e = fourmomboost(1,1) 
              px = fourmomboost(2,1) 
              py = fourmomboost(3,1) 
              pz = fourmomboost(4,1) 
            end if

      end subroutine 




      DOUBLE PRECISION FUNCTION GETTEMPMAX()
      IMPLICIT NONE
C--medium parameters
      COMMON/TEMPMAX/TEMPMAXIMUM
      DOUBLE PRECISION TEMPMAXIMUM
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     & NX,NY,NT,NF,DX,DT,XMAX,XMIN,TMAX
      INTEGER NX,NY,NT,NF
      DOUBLE PRECISION DX,DT,XMAX,XMIN,TMAX
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      common/grid/timesteps(60),tprofile(834,834,60),vtxmap(834,834)
      double precision timesteps,tprofile,vtxmap
      common/gridvel/u(834,834,60),utheta(834,834,60)
      double precision u,utheta
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,MODMED,MEDFILELIST
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON,MODMED,MEDFILELIST
C--function call
      DOUBLE PRECISION GETTEMP
C--hydrodynamic limits
      COMMON /HYDROLIM/ MIDRAPLIM, TMAXLIM, TVELMAXLIM, BOOSTTR, 
     &GLOBALLIMS, PRETAUHYDRO
      DOUBLE PRECISION MIDRAPLIM, TMAXLIM, TVELMAXLIM
      LOGICAL BOOSTTR, GLOBALLIMS, PRETAUHYDRO
      
      !GETTEMPMAX=GETTEMP(0.D0,0.D0,0.D0,TAUI)
      !write(*,*) "Max temp:", tempmaximum
      GETTEMPMAX=TMAXLIM
      END



      DOUBLE PRECISION FUNCTION GETMDMAX()
      IMPLICIT NONE
C--factor to vary Debye mass
      COMMON/MDFAC/MDFACTOR,MDSCALEFAC
      DOUBLE PRECISION MDFACTOR,MDSCALEFAC
      DOUBLE PRECISION GETTEMPMAX
      GETMDMAX=MDSCALEFAC*3.*GETTEMPMAX()
      GETMDMAX=MAX(GETMDMAX,MDFACTOR)
      END



      DOUBLE PRECISION FUNCTION GETMDMIN()
      IMPLICIT NONE
C--medium parameters
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     & NX,NY,NT,NF,DX,DT,XMAX,XMIN,TMAX
      INTEGER NX,NY,NT,NF
      DOUBLE PRECISION DX,DT,XMAX,XMIN,TMAX
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      common/grid/timesteps(60),tprofile(834,834,60),vtxmap(834,834)
      double precision timesteps,tprofile,vtxmap
      common/gridvel/u(834,834,60),utheta(834,834,60)
      double precision u,utheta
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,MODMED,MEDFILELIST
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON,MODMED,MEDFILELIST
C--factor to vary Debye mass
      COMMON/MDFAC/MDFACTOR,MDSCALEFAC
      DOUBLE PRECISION MDFACTOR,MDSCALEFAC
      DOUBLE PRECISION GETTEMPMAX
      GETMDMIN=MDSCALEFAC*3.*TC
      GETMDMIN=MAX(GETMDMIN,MDFACTOR)
      END



      DOUBLE PRECISION FUNCTION GETMSMAX()
      IMPLICIT NONE
      DOUBLE PRECISION GETMDMAX,SQRT
      GETMSMAX=GETMDMAX()/SQRT(2.D0)
      END



      DOUBLE PRECISION FUNCTION GETNATMDMIN()
      IMPLICIT NONE
C--medium parameters
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     & NX,NY,NT,NF,DX,DT,XMAX,XMIN,TMAX
      INTEGER NX,NY,NT,NF
      DOUBLE PRECISION DX,DT,XMAX,XMIN,TMAX
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      common/grid/timesteps(60),tprofile(834,834,60),vtxmap(834,834)
      double precision timesteps,tprofile,vtxmap
      common/gridvel/u(834,834,60),utheta(834,834,60)
      double precision u,utheta
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,MODMED,MEDFILELIST
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON,MODMED,MEDFILELIST
      COMMON/GAMMAMAX/GAMMAMAXIMUM,VELMAXIMUM,DENSITYMAXIMUM,
     &DENSITYMINIMUM, TAUMIN
      DOUBLE PRECISION GAMMAMAXIMUM,VELMAXIMUM,DENSITYMAXIMUM,
     &DENSITYMINIMUM, TAUMIN
C--max rapidity
      common/rapmax2/etamax2
      double precision etamax2
C--factor to vary Debye mass
      COMMON/MDFAC/MDFACTOR,MDSCALEFAC
      DOUBLE PRECISION MDFACTOR,MDSCALEFAC,PI
      DATA PI/3.141592653589793d0/
C--local variables
      DOUBLE PRECISION T,GETMDMIN
      T=GETMDMIN()/(MDSCALEFAC*3.)
      GETNATMDMIN=(2.*6.*NF*D3*2./3. + 16.*ZETA3*3./2.)
     &     *T**3/PI**2
      GETNATMDMIN=min(GETNATMDMIN,DENSITYMINIMUM)
      END



      DOUBLE PRECISION FUNCTION GETLTIMEMAX()
      IMPLICIT NONE
C--medium parameters
      COMMON/LTIME/MODLTIME
      DOUBLE PRECISION MODLTIME
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     & NX,NY,NT,NF,DX,DT,XMAX,XMIN,TMAX
      INTEGER NX,NY,NT,NF
      DOUBLE PRECISION DX,DT,XMAX,XMIN,TMAX
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      common/grid/timesteps(60),tprofile(834,834,60),vtxmap(834,834)
      double precision timesteps,tprofile,vtxmap
      common/gridvel/u(834,834,60),utheta(834,834,60)
      double precision u,utheta
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,MODMED,MEDFILELIST
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON,MODMED,MEDFILELIST
      COMMON/GAMMAMAX/GAMMAMAXIMUM,VELMAXIMUM,DENSITYMAXIMUM,
     &DENSITYMINIMUMDENSITYMINIMUM, TAUMIN
      DOUBLE PRECISION GAMMAMAXIMUM,VELMAXIMUM,DENSITYMAXIMUM,
     &DENSITYMINIMUMDENSITYMINIMUM, TAUMIN
C--max rapidity
      common/rapmax2/etamax2
      double precision etamax2
C--function call
      DOUBLE PRECISION GETTEMPMAX
      !if(medfilelist.eqv..false.) then
      if(.false.) then
        GETLTIMEMAX=TAUI*(GETTEMPMAX()/TC)**3*cosh(etamax2)
      else
      !Fabio: Putting my LTIME
      !write(*,*) "Lifetime whithout boost:",modltime
      GETLTIMEMAX=MODLTIME*cosh(etamax2)
      endif
      END


      DOUBLE PRECISION FUNCTION GETNEFFMAX()
      IMPLICIT NONE
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     & NX,NY,NT,NF,DX,DT,XMAX,XMIN,TMAX
      INTEGER NX,NY,NT,NF
      DOUBLE PRECISION DX,DT,XMAX,XMIN,TMAX
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      common/grid/timesteps(60),tprofile(834,834,60),vtxmap(834,834)
      double precision timesteps,tprofile,vtxmap
      !common/gridvel/ux(834,834,60),uy(834,834,60)
      double precision ux,uy
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,MODMED,MEDFILELIST
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON,MODMED,MEDFILELIST
      common/temperature/tempfac
      double precision tempfac

      COMMON/GAMMAMAX/GAMMAMAXIMUM,VELMAXIMUM,DENSITYMAXIMUM,
     &DENSITYMINIMUMDENSITYMINIMUM, TAUMIN
      DOUBLE PRECISION GAMMAMAXIMUM,VELMAXIMUM,DENSITYMAXIMUM,
     &DENSITYMINIMUMDENSITYMINIMUM, TAUMIN
C--max rapidity
      common/rapmax2/etamax2
      double precision etamax2
C--local variables
      DOUBLE PRECISION PI,GETTEMPMAX
      double precision J0, JR, J3, P0, P1, P2, P3, PR, gamtot
      DATA PI/3.141592653589793d0/
     
      GETNEFFMAX = DENSITYMAXIMUM
      END

      

      SUBROUTINE MYMED()
      IMPLICIT NONE
      COMMON/LTIME/MODLTIME
      DOUBLE PRECISION MODLTIME
      COMMON/TEMPMAX/TEMPMAXIMUM
      DOUBLE PRECISION TEMPMAXIMUM
      COMMON/GAMMAMAX/GAMMAMAXIMUM,VELMAXIMUM,DENSITYMAXIMUM,
     &DENSITYMINIMUM, TAUMIN
      DOUBLE PRECISION GAMMAMAXIMUM,VELMAXIMUM,DENSITYMAXIMUM,
     &DENSITYMINIMUM, TAUMIN
      DOUBLE PRECISION PTEMPERATURE,PVELOCITY,PDENSITY
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     & NX,NY,NT,NF,DX,DT,XMAX,XMIN,TMAX
      INTEGER NX,NY,NT,NF
      COMMON/logfile/logfid
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,MODMED,MEDFILELIST
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON,MODMED,MEDFILELIST
      common/boostmed/boost
      logical boost
      common/rapmax2/etamax2
      double precision etamax2
      COMMON/MEDFILEC/MEDFILE,NLIST,endoff
      CHARACTER*200 MEDFILE
      INTEGER NLIST
      LOGICAL ENDOFF
      INTEGER logfid
      DOUBLE PRECISION DX,DT,XMAX,XMIN,TMAX
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      common/grid/timesteps(60),tprofile(834,834,60),vtxmap(834,834)
      double precision timesteps,tprofile,vtxmap
      !common/gridvel/ux(834,834,60),uy(834,834,60)
      !double precision ux,uy
      common/gridvel/u(834,834,60),utheta(834,834,60)
      double precision u,utheta           
      double precision gridx(834,834),gridy(834,834)
      CHARACTER DUMMIE,CONTENT*100
      DOUBLE PRECISION tempsum,entropy
      INTEGER I,J,K,POS,II,kk,kkk,length
      logical ltime
      double precision hightemp, cdensmin, cdensmax, cvel, cveltot,
     &densconst, randuz, randpz, randpr, cdens, pyr, maxu
      integer nrandpoints, irand
      double precision densityarray(834,834,60)
      double precision PI
      DATA PI/3.141592653589793d0/
C--hydrodynamic limits
      COMMON /HYDROLIM/ MIDRAPLIM, TMAXLIM, TVELMAXLIM, BOOSTTR, 
     &GLOBALLIMS, PRETAUHYDRO
      DOUBLE PRECISION MIDRAPLIM, TMAXLIM, TVELMAXLIM
      LOGICAL BOOSTTR, GLOBALLIMS, PRETAUHYDRO
C--grid parameters
      common/gridpar/ gdt,gdx,gxmax,gxmin,gnx,gny,gnt
      double precision gdt,gdx,gxmax,gxmin,s
      integer gnx,gny,gnt,probcounter

      NX=834
      dx=50.d0/834.d0

      hightemp=0.d0


      do i=1,834
        do j=1,834
          do k=1,60
            tprofile(i,j,k)=0.d0

          end do

          gridx(i,j)=-25.d0+dx*(i-1)
          gridy(i,j)=-25.d0+dx*(j-1)

        end do
      end do
      
      call reader(medfile,834,60,timesteps,tprofile,u,utheta)

      call read_initvtx(834, tprofile, vtxmap)



      velmaximum = 0.d0
      hightemp = 0.d0

C--Loop that finds the medium lifetime and also
C--its evolution in entropy and temperature
C--as well as its highest temperature

      ! TAU0 is k = 2, k = 1 is a matrix of zeros
      TAUMIN = timesteps(2)
      do k=2,60      
        ltime=.true.
      
        do kk=1,NX
          do kkk=1,NX 
            if(tprofile(kk,kkk,k).ge.tc) then
              ltime=.false.

              ! Calculate other limts
              if (BOOSTTR .eqv. .true.) then
                cvel = u(kk, kkk, k)
              else 
                cvel = 0
              end if 

              if (cvel.gt.velmaximum) then
                velmaximum = cvel
              end if

            end if

            if (tprofile(kk,kkk,k) .gt. hightemp) then
              hightemp=tprofile(kk,kkk,k)
            end if

          enddo
        enddo
     
        write(*,*) "Ltime:",ltime

        if(.not.ltime.and.timesteps(k+1).gt.timesteps(k)) then
          write(*,*) "Lifetime not reached:",timesteps(k)
          nt=k+1
          modltime=timesteps(k+1)

        else
          write(*,*) "Lifetime reached:",timesteps(k)
        endif

      end do

      if (velmaximum.gt.TVELMAXLIM) then
        TVELMAXLIM = velmaximum
        write(logfid,*) 'Maximum trans u for hydro profile > TVELMAXLIM'
        write(logfid,*) 'Extrapolating trans u limit'
        write(logfid,*)
      end if 


      write(*,*)

      ! Define limits
      if (.not. GLOBALLIMS) then
        write(logfid,*) 'Using "medium profile-global" limits' 
        write(logfid,*) 'Similar to JEWEL original impl. (vs global)' 
        write(logfid,*) 'This could result in GETDELTAT weight errors'
        
        TMAXLIM = hightemp
        TVELMAXLIM = velmaximum

        if (boost) then
          MIDRAPLIM = ETAMAX2
        end if

        write(logfid,*) 'TMAXLIM = ', TMAXLIM
        write(logfid,*) 'TVELMAXLIM = ', TVELMAXLIM
        write(logfid,*) 'MIDRAPLIM = ', MIDRAPLIM
        write(logfid,*) 
      end if
        
      ! Setup global density limits
      ! neff = n * p_mu u^mu / p0, check GETNEFF
      ! neff_max is highest at tau0 (ux = uy = 0)
      ! neff_min is lowest when u_T is max
      ! due to temperature dependence
      gammamaximum=sqrt(1 + TVELMAXLIM ** 2)
      maxu = sqrt(TVELMAXLIM ** 2 + sinh(MIDRAPLIM) ** 2)
      densconst = (2.*6.*NF*D3*2./3. + 16.*ZETA3*3./2.)/PI**2

      densitymaximum = densconst * TMAXLIM ** 3 *
     &(sqrt(1 + sinh(MIDRAPLIM) ** 2) + sinh(MIDRAPLIM)) 

      densityminimum = densconst * TC ** 3 *
     &(sqrt(1 + maxu ** 2) - maxu) 

      write(*,*)
      write(*,*) "TAU max (for limits)= ", timesteps(nt - 1)
      write(*,*) "Highest temp = ", hightemp
      write(*,*) "Highest u_r = ", velmaximum
      write(*,*) "Highest gamma_r = ", gammamaximum
      write(*,*) "Highest u (with Bjorken) = ", maxu
      write(*,*) "Highest gamma (with Bjorken) = ", sqrt(1 + 
     &maxu ** 2)
      write(*,*) "Highest effective density = ", densitymaximum
      write(*,*) "Lowest effective density = ", densityminimum

      WRITE(*,*) "Temperature profile read succesfully :)"
      WRITE(*,*) 

      END


      DOUBLE PRECISION FUNCTION MEDPART(X4,Y4,Z4,T4)
      IMPLICIT NONE
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,MODMED,MEDFILELIST
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON,MODMED,MEDFILELIST
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     & NX,NY,NT,NF,DX,DT,XMAX,XMIN,TMAX
      INTEGER NX,NY,NT,NF
      DOUBLE PRECISION DX,DT,XMAX,XMIN,TMAX
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      common/grid/timesteps(60),tprofile(834,834,60),vtxmap(834,834)
      double precision timesteps,tprofile,vtxmap
      !common/gridvel/ux(834,834,60),uy(834,834,60)
      double precision ux,uy
      DOUBLE PRECISION X4,Y4,Z4,T4
      DOUBLE PRECISION STEP
      DOUBLE PRECISION TAU,interpolate
      STEP=(XMAX-XMIN)/(NX-1)
      TAU=SQRT(T4**2-Z4)
      TAU=0.0d0
      MEDPART=interpolate(X4,Y4,tau,1)
      END

      DOUBLE PRECISION FUNCTION MEDDERIV(XVAL,W)
      IMPLICIT NONE
      DOUBLE PRECISION XVAL
      INTEGER W
C--medium parameters
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,MODMED,MEDFILELIST
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON,MODMED,MEDFILELIST
C--variables for integration
      COMMON/INTEG/B,R
      DOUBLE PRECISION B,R

      IF (W.EQ.1) THEN
C--XVAL corresponds to z-coordinate
       MEDDERIV=N0/(1+EXP((SQRT(B**2+XVAL**2)-R)/D))
      ELSE 
       MEDDERIV=0.D0
      ENDIF
      END

      DOUBLE PRECISION FUNCTION INTEGRATE(TPROF)
      IMPLICIT NONE
      DOUBLE PRECISION TPROF,XMAX,XMIN,YMAX,YMIN,STEP
      DOUBLE PRECISION TERM
      DOUBLE PRECISION INTERMED(100)
      INTEGER N,I,J
      XMAX=10.d0
      XMIN=-10.d0
      YMAX=10.d0
      YMIN=-10.d0
      N=FLOOR(MAX((XMAX-XMIN+2*STEP)/(2*STEP),0.d0))
      INTEGRATE=0.d0
      DO 10 I=1,N
          TERM=0.d0
          INTERMED(I)=0.d0
          DO 20 J=1,N
              TERM=(STEP/3)*(TPROF(I,J)+TPROF(I,J+2)
     &        +4*TPROF(I,J+1))
              INTERMED(I)=INTERMED(I)+TERM
20        CONTINUE
10    CONTINUE
      TERM=0.d0
      DO 30 I=1,N
          TERM=(STEP/3)*(INTERMED(I)+INTERMED(I+2)
     &    +4*INTERMED(I+1))
          INTEGRATE=INTEGRATE+TERM
30    CONTINUE
      END


