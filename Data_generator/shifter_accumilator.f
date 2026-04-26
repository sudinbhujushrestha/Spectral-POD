C====================================================================
      PROGRAM STATISTIC

C     1.AVERAGE EACH VARIABLE W.R.T. PSI
C     2.SHIFT ALL THE VARIABLES TO THE SAME PHASE
C     3.CALCULATE PHASE AVERAGE
C     4.CALCULATE PHASE AVERAGE OF FLUCTUATION

      IMPLICIT NONE

      INCLUDE "mpif.h"

      INTEGER MYID,NUMPROCS,IERR,STATUS(MPI_STATUS_SIZE)
      
      INTEGER NXMAX,NYMAX,NZMAX,NZMAX2,NXMOD,NYMOD,NZ,NZ1,NZ2
      INTEGER NF,NTF,NWORK,NTRIGS,NTIME,NCPU
      INTEGER NCREQ,JJ1,JJ2,NYEND
      INTEGER I,J,K,IT,ITT,KF,ICPU,ID,IL,IOUT_STAG
      INTEGER IWAVY,NWAVE,IOUT1,IOUT2,IOUT3,IIN1,IIN2,IDUM

      REAL PEX,PEY,HBAR,TWOPI,DUM,XL,X,DX,Z
      REAL TIME,TIME0,TIMEWAVY,DTS,DT
      REAL FSIN,FCOS,HKA,HA,HK,HOMEG,TCOEF,VPHASE,CRAT
      REAL WW1,WW2
      REAL WW,WWF,QS
      REAL RE,RESBOT,RESTOP,USBOT,USTOP,CURV
      REAL CLBETA,CLGAMA
      REAL DLZ
      REAL, PARAMETER :: PI = 2.0 * ASIN(1.0)
      REAL theta
      REAL theta_rad


      PARAMETER(NXMAX=195,NYMAX=195,NZMAX=66,NZMAX2=NZMAX+1)
      PARAMETER(NXMOD=192,NYMOD=192,NZ=65)


C Ask professor what this NF and NTF stands for, it is not clear  

C CAREFUL !!! if all your files are in 1 single folder, and your 1 set of files has all the time sets then set NF=1
C Then set NTF to the number of time steps in that single file


C=======================================================================
C --- DESCRIPTION OF TIME STEP PARAMETERS ---
C=======================================================================
C
C  NF: Number of File sets (or snapshot folders).
C      In the original data structure, this represents the number of
C      separate directories (e.g., run_3225/, run_3275/) that the
C      program loops through. For this code, NF = 39. 
C
C  NTF: Number of Time steps per File set.
C       This represents the number of individual time step records
C       contained within each of the NF sets. For this code, NTF = 10. 
C
C  NTIME: Total Number of Time Steps.
C         This is the total number of snapshots that will be processed.
C         NTIME = NF * NTF = 39 * 10 = 390. [cite: 324]
C
C  ITT: Continuous Time Step Counter.
C       This variable serves as a single "master index" that converts
C       the two nested loop counters (KF and IT) into one continuous
C       count that runs from 1 to NTIME.
C       It is calculated as: ITT = (KF - 1) * NTF + IT.
C
C=======================================================================
      PARAMETER(NF=1,NTF=500)
      PARAMETER(NWORK=2*NXMAX*NYMAX)
      PARAMETER(NTRIGS=4*MAX(NXMAX,NYMAX))
      PARAMETER(NTIME=NF*NTF)
      PARAMETER(PEX=0.6666667,PEY=0.6666667,HBAR=1.57079)
      PARAMETER(NCPU=39)
      PARAMETER(IWAVY=3,NWAVE=6)
      PARAMETER(HKA=0.1,CRAT=12.0)
      PARAMETER(TIMEWAVY=1500.,DT=0.0075)
      PARAMETER(RESBOT=180.)
      PARAMETER(CLBETA=1.5,CLGAMA=20.0)
      PARAMETER(theta=0.0)

      CHARACTER*10, CFILE(NCPU),CFILE2(NCPU)
      CHARACTER*12, CPATH(NF)
 

      REAL U(NXMAX,NYMAX,NZMAX)


      REAL V(NXMAX,NYMAX,NZMAX)
      
      REAL W(NXMAX,NYMAX,NZMAX)

      REAL W_orig(NXMAX,NYMAX,NZMAX)

      REAL USHIFT(NXMAX,NYMAX,NZMAX)
      REAL VSHIFT(NXMAX,NYMAX,NZMAX)
      real WSHIFT(NXMAX,NYMAX,NZMAX)
      REAL PP(NXMAX,NYMAX,NZMAX)
      REAL PT(NXMAX,NYMAX,NZMAX)
     
      REAL PSHIFT(NXMAX,NYMAX,NZMAX)


      REAL HH(NXMAX,NYMAX),HT(NXMAX,NYMAX),HX(NXMAX,NYMAX)


      REAL ZZ1(NZMAX),ZZ2(NZMAX),ZZ(NZMAX),DZ(NZMAX),ZW(NZMAX)

      REAL FTMP(NXMAX)
      REAL WORK(NWORK)
      REAL TRIGSX(NTRIGS),TRIGSY(NTRIGS)
      INTEGER IFAX(19)

      DATA IOUT1,IOUT2,IOUT3/41,42,320/

      CALL MPI_INIT(IERR)     
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)


      theta_rad = theta * PI / 180.0


C-------------------------------
C     TURBULENCE PARAMETERS
C-------------------------------
      RESTOP=0.0
      IF(ABS(RESTOP).LT.1.E-6)THEN
         USBOT=1./(2.5*ALOG(HBAR*RESBOT)+5.0)
         USTOP=0.
         RE=RESBOT*(2.5*ALOG(HBAR*RESBOT)+5.0)
         GOTO 10
      ENDIF

      IF(ABS(RESBOT).LT.1.E-6)THEN
         USBOT=0.
         USTOP=1./(2.5*ALOG(HBAR*RESTOP)+5.0)
         RE=RESTOP*(2.5*ALOG(HBAR*RESTOP)+5.0)
         GOTO 10
      ENDIF

      USBOT=1./((2.5*ALOG(HBAR*RESBOT**2/(RESBOT+RESTOP))+5.0)
     1     +(2.5*ALOG(HBAR*RESTOP**2/(RESBOT+RESTOP))+5.0)*RESTOP
     1     /RESBOT)
      USTOP=1./((2.5*ALOG(HBAR*RESBOT**2/(RESBOT+RESTOP))+5.0)*RESBOT
     1     /RESTOP+(2.5*ALOG(HBAR*RESTOP**2/(RESBOT+RESTOP))+5.0))
      RE=RESBOT*(2.5*ALOG(HBAR*RESBOT**2/(RESBOT+RESTOP))+5.0)
     1  +RESTOP*(2.5*ALOG(HBAR*RESTOP**2/(RESBOT+RESTOP))+5.0)          

 10   CONTINUE

C-----END HERE

C-----------------------------------
C     COEF FOR BOTTOM WAVY WALL
C-----------------------------------
C     HK -- WAVE NUMBER
C     HA -- WAVE AMPLITUDE
C     HOMEG -- WAVE ANGLE FREQUENCY
C     NWAVE -- NUMBER OF WAVES IN X DIRECTION
C     HKA -- WAVE SLOPE
C     VPHASE -- WAVE PHASE SPEED
C     CRAT -- RATIO OF C TO USTOP
C     TIMEWAVY -- TIME FOR BOTTOM TO START MOVING

      HK=PEX*NWAVE
      HA=HKA/HK
      VPHASE=CRAT*USBOT
      HOMEG=HK*VPHASE

C-----END HERE

C-------------------------------
C     INITIALIZE VARIABLES
C-------------------------------
c
      NZ1=NZ+1
      NZ2=NZ+2

      TWOPI=2.*ACOS(-1.)
      XL=TWOPI/PEX
      DX=XL/FLOAT(NXMOD)

c Read OUTALL u, v, w, pp, WRITE(92+MYID*1000,*), based on the number of processors used, CFILES needs to set accordingly, the length of the name of the string is set at the beginning of the code
      CFILE(1)='fort.92'
      CFILE(2)='fort.1092'
      CFILE(3)='fort.2092'
      CFILE(4)='fort.3092'
      CFILE(5)='fort.4092'
      CFILE(6)='fort.5092'
      CFILE(7)='fort.6092'
      CFILE(8)='fort.7092'
      CFILE(9)='fort.8092'
      CFILE(10)='fort.9092'
      CFILE(11)='fort.10092'
      CFILE(12)='fort.11092'
      CFILE(13)='fort.12092'
      CFILE(14)='fort.13092'
      CFILE(15)='fort.14092'
      CFILE(16)='fort.15092'
      CFILE(17)='fort.16092'
      CFILE(18)='fort.17092'
      CFILE(19)='fort.18092'
      CFILE(20)='fort.19092'
      CFILE(21)='fort.20092'
      CFILE(22)='fort.21092'
      CFILE(23)='fort.22092'
      CFILE(24)='fort.23092'
      CFILE(25)='fort.24092'
      CFILE(26)='fort.25092'
      CFILE(27)='fort.26092'
      CFILE(28)='fort.27092'
      CFILE(29)='fort.28092'
      CFILE(30)='fort.29092'
      CFILE(31)='fort.30092'
      CFILE(32)='fort.31092'
      CFILE(33)='fort.32092'
      CFILE(34)='fort.33092'
      CFILE(35)='fort.34092'
      CFILE(36)='fort.35092'
      CFILE(37)='fort.36092'
      CFILE(38)='fort.37092'
      CFILE(39)='fort.38092'





c Read OUTALL ETA0, ETA, HH, WRITE(95+MYID*1000,*), again based on the number of processors used, CFILE2 needs to be set accordingly
c but remember in the main code 
c            ETA0(I,J)=0.
c            ETA(I,J)=ETA0(I,J)+HH(I,J)
c so there is no point of ETA0 and ETA, because of which ETA0 and ETA willb e read as dummy variables later on below
      CFILE2(1)='fort.95'
      CFILE2(2)='fort.1095'
      CFILE2(3)='fort.2095'
      CFILE2(4)='fort.3095'
      CFILE2(5)='fort.4095'
      CFILE2(6)='fort.5095'
      CFILE2(7)='fort.6095'
      CFILE2(8)='fort.7095'
      CFILE2(9)='fort.8095'
      CFILE2(10)='fort.9095'
      CFILE2(11)='fort.10095'
      CFILE2(12)='fort.11095'
      CFILE2(13)='fort.12095'
      CFILE2(14)='fort.13095'
      CFILE2(15)='fort.14095'
      CFILE2(16)='fort.15095'
      CFILE2(17)='fort.16095'
      CFILE2(18)='fort.17095'
      CFILE2(19)='fort.18095'
      CFILE2(20)='fort.19095'
      CFILE2(21)='fort.20095'
      CFILE2(22)='fort.21095'
      CFILE2(23)='fort.22095'
      CFILE2(24)='fort.23095' 
      CFILE2(25)='fort.24095'
      CFILE2(26)='fort.25095'
      CFILE2(27)='fort.26095'
      CFILE2(28)='fort.27095'
      CFILE2(29)='fort.28095'
      CFILE2(30)='fort.29095'
      CFILE2(31)='fort.30095'
      CFILE2(32)='fort.31095'
      CFILE2(33)='fort.32095'
      CFILE2(34)='fort.33095'
      CFILE2(35)='fort.34095'
      CFILE2(36)='fort.35095'
      CFILE2(37)='fort.36095'
      CFILE2(38)='fort.37095'
      CFILE2(39)='fort.38095'





C If the data files are inside a folder, we need to give the path as a string that is the purpose of CPATH, the length of CPATH is set at the beinning of the code
      CPATH(1)='../work_dir/'
c      CPATH(2)='../../0theta_11505-11880_waveon_CRAT6_ISTART1_192grid/'
c      CPATH(3)='../../0theta_11880-12255_waveon_CRAT6_ISTART1_192grid/'
c      PATH(4)='../run_3300/'
c      PATH(5)='../run_3325/'
c      PATH(6)='../run_3350/'
c      PATH(7)='../run_3375/'
c      PATH(8)='../run_3400/'
c      PATH(9)='../run_3425/'
c      PATH(10)='../run_3450/'
c      PATH(11)='../run_3475/'
c      PATH(12)='../run_3500/'
c      PATH(13)='../run_3525/'
c      PATH(14)='../run_3550/'
c      PATH(15)='../run_3575/'
c      PATH(16)='../run_3600/'
c      PATH(17)='../run_3625/'
c      PATH(18)='../run_3650/'
c      PATH(19)='../run_3675/'
c      PATH(20)='../run_3700/'
c      PATH(21)='../run_3725/'
c      PATH(22)='../run_3750/'
c      PATH(23)='../run_3775/'
c      PATH(24)='../run_3800/'
c      PATH(25)='../run_3825/'
c      PATH(26)='../run_3850/'
c      PATH(27)='../run_3875/'
c      PATH(28)='../run_3900/'
c      PATH(29)='../run_3925/'
c      PATH(30)='../run_3950/'
c      PATH(31)='../run_3975/'
c      PATH(32)='../run_4000/'
c      PATH(33)='../run_4025/'
c      PATH(34)='../run_4050/'
c      PATH(35)='../run_4075/'
c      PATH(36)='../run_4100/'
c      PATH(37)='../run_4125/'
c      PATH(38)='../run_4150/'
c      PATH(39)='../run_4175/'


C--------------------------------
C     GENERAGE GRID SYSTEM
C--------------------------------

      CALL CLUSTER_NL_1(NZ,NZMAX,1.0,CLBETA,CLGAMA,ZZ1,ZZ2,ZZ,ZW,DZ)

C-----END HERE

C---------------------------------------------------------
C     WRITE Z COORDINATES FOR PYTHON SPOD WEIGHTS
C---------------------------------------------------------
      IF(MYID .EQ. 0) THEN
         OPEN(UNIT=77, FILE='z_coordinates.txt', STATUS='UNKNOWN',
     1        FORM='FORMATTED')
         DO K=1,NZ
            WRITE(77,'(ES24.16)') ZZ(K)
         ENDDO
         CLOSE(77)
         PRINT*, "z_coordinates.txt written with ", NZ, " points."
      ENDIF

c Each of the **092 and **095 files are opened and are given unique unit numbers per file
      DO KF=1,NF
         DO ICPU=1,NCPU
            ID=ICPU-1
            IIN1=(92+ID*1000)*100+KF
            OPEN(UNIT=IIN1,FILE=CPATH(KF)//CFILE(ICPU),STATUS='OLD')
            IIN2=(95+ID*1000)*100+KF
            OPEN(UNIT=IIN2,FILE=CPATH(KF)//CFILE2(ICPU),STATUS='OLD')
         ENDDO
      ENDDO



      PRINT*, "OPEN FILE SUCCESSFULLY!"

      IF(MYID .EQ. 0) THEN 
         OPEN(UNIT=50, FILE='flattened_data.bin', STATUS='UNKNOWN', 
     1        FORM='UNFORMATTED', ACCESS='STREAM')
         PRINT*, "Opened flattened_data.bin for writing..."
      ENDIF


C Sudin_edit_open*~*~**~*~**~*~**~*~**~*~**~*~**~*~**~*~**~*~**~*~**~*~**~*~**~*~**~*~**~*~**~*~**~*~**~*~**~*~**~*~**~*~**~*~**~*~**~*
C*~*~**~*~**~*~**~*~**~*~**~*~**~*~**~*~**~*~**~*~**~*~**~*~**~*~**~*~**~*~**~*~**~*~**~*~**~*~**~*~**~*~**~*~**~*~**~*~**~*~**~*~**~*~

      
C skipping the 1st loop of which produces UPM

      DO 200 KF=1,NF
         
         DO 100 IT=1,NTF

            ITT=(KF-1)*NTF+IT

C----------------------------
C     READ IN 3D DOMAIN
C----------------------------

            DO ICPU=1,NCPU

               ID=ICPU-1
               
               NCREQ=NCPU-(NYMAX-NYMOD)*NCPU/NYMAX
               
               JJ1=NYMAX/NCPU-MOD(NYMAX-NYMOD,NYMAX/NCPU)
               JJ2=NYMAX/NCPU
               
               IF(ID.EQ.NCREQ-1) THEN
                  NYEND=JJ1
               ELSE IF(ID.LT.NCREQ-1) THEN
                  NYEND=JJ2
               ELSE
                  NYEND=0
               ENDIF

               IIN1=(92+ID*1000)*100+KF
               IIN2=(95+ID*1000)*100+KF
c Reading the variables from file 92 and 95 
               READ(IIN1,*) TIME,IDUM,IDUM,IDUM,DUM,DUM,DUM,DUM
               READ(IIN1,*) (((U(I,ID*NYMAX/NCPU+J,K),I=1,NXMOD),
     1              J=1,NYEND),K=1,NZ1)
               READ(IIN1,*) (((V(I,ID*NYMAX/NCPU+J,K),I=1,NXMOD),
     1              J=1,NYEND),K=1,NZ1)
               READ(IIN1,*) (((W(I,ID*NYMAX/NCPU+J,K),I=1,NXMOD),
     1              J=1,NYEND),K=1,NZ1)
               READ(IIN1,*) (((PT(I,ID*NYMAX/NCPU+J,K),
     1              I=1,NXMOD),J=1,NYEND),K=1,NZ1)
               READ(IIN2,*) ((DUM,I=1,NXMOD),J=1,NYEND)
               READ(IIN2,*) ((DUM,I=1,NXMOD),J=1,NYEND)
               READ(IIN2,*) ((HH(I,ID*NYMAX/NCPU+J),I=1,NXMOD),
     1              J=1,NYEND)

            ENDDO

            WRITE(*,5051) TIME
 5051        FORMAT(' T=',F12.4)
            PRINT*, "3D Domain data readin 1st time completed !"



C-----END HERE

C-------------------------------
C     4 DIFFERENT WAVY WALL
C-------------------------------

            TIME0=TIME-TIMEWAVY


C PHASE_AVERAGE_1D and PHASE_AVERAGE_2D : The averaged data obtained from AVERAGE_PSI_2D and AVERAGE_PSI_3D spans over N waves, these function averaged them onto a single wave 
           
            IF(IWAVY.EQ.1) THEN
               DTS=0.
            ELSE
C DTS is the number of wave cycles (or domain fractions) the wave has propagated through since time zero, expressed in units of domain length.      
c If 1 time unit is defined as the time required for a wave to travel one wavelength, then TIME0 is the total number of wavelengths the wave 
c has traveled. Dividing TIME0 by NWAVE gives the number of domain-lengths the wave has propagated through — since the domain contains NWAVE wavelengths.”  
               DTS=-TIME0/FLOAT(NWAVE)

            ENDIF
            
C     IWAVY=1: SOLID WAVY WALL
            
            IF(IWAVY.EQ.1) THEN
               
               DO I=1,NXMOD
                  DO J=1,NYMOD
                     X=(I-NXMOD/2-1)*DX
                     FSIN=SIN(HK*X)
                     FCOS=COS(HK*X)
                     HH(I,J)=-HA*FSIN
                     HT(I,J)=0.
                     HX(I,J)=-HA*HK*FCOS
                  ENDDO
               ENDDO
               
               GOTO 999
               
            ENDIF

C     IWAVY=2: VERTICAL MOVING WAVY WALL
      
            IF(IWAVY.EQ.2) THEN
               
               DO I=1,NXMOD
                  DO J=1,NYMOD
                     X=(I-NXMOD/2-1)*DX
                     FSIN=SIN(HK*X-HOMEG*TIME0)
                     FCOS=COS(HK*X-HOMEG*TIME0)
                     HH(I,J)=-HA*FSIN
                     HT(I,J)=HA*HOMEG*FCOS
                     HX(I,J)=-HA*HK*FCOS
                  ENDDO
               ENDDO
               
               GOTO 999
               
            ENDIF

C     IWAVY=3: WATER WAVE SURFACE

            IF(IWAVY.EQ.3) THEN
               
               DO I=1,NXMOD
                  DO J=1,NYMOD
                     X=(I-NXMOD/2-1)*DX
                     FSIN=SIN(HK*X-HOMEG*TIME0)
                     FCOS=COS(HK*X-HOMEG*TIME0)
                     HH(I,J)=-HA*FSIN
                     HT(I,J)=HA*HOMEG*FCOS
                     HX(I,J)=-HA*HK*FCOS
                  ENDDO
               ENDDO
               
               GOTO 999
               
            ENDIF

C     IWAVY=4: WATER WAVE SURFACE WITH WIND DRIFT

            IF(IWAVY.EQ.4) THEN
               
               DO I=1,NXMOD
                  DO J=1,NYMOD
                     X=(I-NXMOD/2-1)*DX
                     FSIN=SIN(HK*X-HOMEG*TIME0)
                     FCOS=COS(HK*X-HOMEG*TIME0)
                     HH(I,J)=-HA*FSIN
                     HT(I,J)=HA*HOMEG*FCOS
                     HX(I,J)=-HA*HK*FCOS
                  ENDDO
               ENDDO
               
               GOTO 999
               
            ENDIF

 999        CONTINUE

C----------------------------------
C     CALCULATE REAL PRESSURE
C----------------------------------

c            CALL PRESSURE(NXMAX,NYMAX,NZMAX,NXMOD,NYMOD,NZ,PEX,PEY,DZ,
c     1           ZZ,WORK,TRIGSX,TRIGSY,IFAX,PT,PP,HBAR,HH,HX,RE,DT)
            
c            PRINT*, 'Real pressure calculated !'

c            CALL PHASE_SHIFT_3D(NXMOD,NYMOD,NZ,NXMAX,NYMAX,NZMAX,NTIME,PP,
c     1           PSHIFT,-DTS,HOMEG,WORK,TRIGSX,IFAX)

            CALL PHASE_SHIFT_3D(NXMOD,NYMOD,NZ,NXMAX,NYMAX,NZMAX,NTIME,U,
     1           USHIFT,-DTS,HOMEG,WORK,TRIGSX,IFAX)
            CALL PHASE_SHIFT_3D(NXMOD,NYMOD,NZ,NXMAX,NYMAX,NZMAX,NTIME,V,
     1           VSHIFT,-DTS,HOMEG,WORK,TRIGSX,IFAX)

c shifting the W from the staggered grid to the same grid as U and V for the convenience of calculating W^2 and UW, VW
            DO K=1,NZ
               DO I=1,NXMOD
                  DO J=1,NYMOD
                     IF(K.EQ.1) THEN  
                        W_orig(i,j,K)=W(i,j,1)
                     ELSEIF (K.EQ.NZ) THEN
                        W_orig(i,j,K)=W(i,j,K-1)  
                     ELSE
                        W_orig(i,j,K)=0.5*(W(I,J,K-1)+W(I,J,K))
                  ENDIF        
                  ENDDO
               ENDDO
            ENDDO

            CALL PHASE_SHIFT_3D(NXMOD,NYMOD,NZ,NXMAX,NYMAX,NZMAX,NTIME,W_orig,
     1           WSHIFT,-DTS,HOMEG,WORK,TRIGSX,IFAX)


C           ---------------------------------------------------------
C           WRITE 3D DOMAIN TO THE OPEN BINARY FILE
C           ---------------------------------------------------------
            IF(MYID .EQ. 0) THEN 
               

               WRITE(50) (((USHIFT(I,J,K), I=1,NXMOD), J=1,NYMOD), 
     1              K=1,NZ)
               WRITE(50) (((VSHIFT(I,J,K), I=1,NXMOD), J=1,NYMOD), 
     1              K=1,NZ)
               WRITE(50) (((WSHIFT(I,J,K), I=1,NXMOD), J=1,NYMOD), 
     1              K=1,NZ)

               PRINT*, "Appended snapshot ", ITT
               
            ENDIF

 100     CONTINUE
 200  CONTINUE

C-----END HERE

      IF(MYID .EQ. 0) THEN
         CLOSE(50)
         PRINT*, "Finished writing all snapshots to flattened_data.bin"
      ENDIF

      DO KF=1,NF
         DO ICPU=1,NCPU
            ID=ICPU-1
            IIN1=(92+ID*1000)*100+KF
            IIN2=(95+ID*1000)*100+KF
            CLOSE(IIN1)
            CLOSE(IIN2)
         ENDDO
      ENDDO

      CALL MPI_FINALIZE(IERR)
      STOP
      END

C=====MAIN PROGRAM STATISTIC END HERE







C==========================================================================
      SUBROUTINE CLUSTER_NL_1(NZ,NZMAX,ZL,CLBETA,CLGAMA,ZZ1,ZZ2,ZZ,ZW,
     1     DZ)

C     THIS VERSION ONLY CLUSTER GRID NEAR K=1

C
C     CLGAMA: LARGEST GRID RATION
C     CLBETA: CLUSTERING PAPAMETER: SMALLER VALUE -> MORE CLUSTERED
C
      REAL ZZ(*),ZW(*),DZ(*),ZZ1(*),ZZ2(*)
C
      NEVEN=4
CC--
C      DO 10 K=NEVEN,NZ-NEVEN
C       ETAK=(K-NEVEN)/(NZ-2.*NEVEN)
C       ZZ1(K)=((1+CLBETA)*((CLBETA+1)/(CLBETA-1))**(2*ETAK-1.)
C     1   +1.-CLBETA)/(2.*(1.+((CLBETA+1.)/(CLBETA-1.))**(2*ETAK-1.)))
C 10   CONTINUE
CC--#####################3
C
C
CC--
C      a=1./(2.*clbeta)*alog((1+(exp(clbeta)-1.)*0.5)
C     1  /(1+(exp(-clbeta)-1.)*0.5))
CC--####################################
      a = 0.5
cc--@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      do 11 k=NEVEN+1,nz-NEVEN-1
      ETAK=(K-NEVEN)/(NZ-2.*NEVEN)/2.
      ZZ1(K)=0.5*(1.+sinh(clbeta*(etak-a))/sinh(clbeta*a))
 11      continue
      zz1(NEVEN)=0.
      zz1(nz-NEVEN)=0.5
C
      TWOPI=2.*ACOS(-1.)
      DO 18 K=NEVEN,NZ-NEVEN
 18         DZ(K)=(1.-COS(ZZ1(K)*TWOPI))*(CLGAMA-1.)/2.+1.
C
      DO 20 K=1,NEVEN-1
 20         DZ(K)=DZ(NEVEN)
C
CC--
      DO 30 K=NZ-NEVEN+1,NZ
 30         DZ(K)=DZ(NZ-NEVEN)
CC--@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
      DZ(1)=DZ(1)/2.
      DZ(NZ-1)=DZ(NZ-1)/2.
      DZ(NZ)=DZ(NZ)/2.
C
      ZZ2(1)=0.
      DO 40 K=1,NZ
 40         ZZ2(K+1)=ZZ2(K)+DZ(K)
C
      DO 45 K=1,NZ+1
c 45   ZZ(K)=(ZZ2(K)-ZZ2(1))/(ZZ2(NZ)-ZZ2(1))*ZL-ZL
 45         ZZ(K)=(ZZ2(K)-ZZ2(1))/(ZZ2(NZ)-ZZ2(1))*ZL
C
      DO 50 K=1,NZ
 50         DZ(K)=ZZ(K+1)-ZZ(K)
C
      WRITE(29,299)
 299    FORMAT('VARIABLES=K,ZZ,DZ,DZRATIO,RDZ,RATIOMAX')
      WRITE(29,399)
 399    FORMAT('ZONE T=""')
C
      UNIT=1
      RMAX=DZ((NZ-1)/2)/DZ(1)
      DO 60 K=NZ,2,-1
       WRITE(29,99)K,ZZ(K),DZ(K),DZ(K)/DZ(K-1),1./DZ(K),RMAX
 99        FORMAT(I5,5E12.4)
 60           CONTINUE
      WRITE(29,99)1,ZZ(1),DZ(1),UNIT,1./DZ(1),RMAX

      DO K=1,NZ
         IF(K.EQ.1) THEN
            ZW(K)=ZZ(K)
         ELSE IF(K.EQ.NZ-1) THEN
            ZW(K)=ZZ(K+1)
         ELSE IF(K.EQ.NZ) THEN
            ZW(K)=ZZ(K+1)+DZ(K)
         ELSE
            ZW(K)=(ZZ(K)+ZZ(K+1))/2.
         ENDIF
      ENDDO
C
      RETURN
      END

C=====SUBROUTINE CLUSTER_NL_1 END HERE

C========================================================================
      SUBROUTINE CLUSTER_NL(NZ,NZMAX,ZL,CLBETA,CLGAMA,ZZ1,ZZ2,ZZ,ZW,DZ)
C
C     CLGAMA: LARGEST GRID RATION
C     CLBETA: CLUSTERING PAPAMETER: SMALLER VALUE -> MORE CLUSTERED
C
      REAL ZZ(*),ZW(*),DZ(*),ZZ1(*),ZZ2(*)
C
      NEVEN=4
      DO 10 K=NEVEN,NZ-NEVEN
       ETAK=(K-NEVEN)/(NZ-2.*NEVEN)
       ZZ1(K)=((1+CLBETA)*((CLBETA+1)/(CLBETA-1))**(2*ETAK-1.)
     1   +1.-CLBETA)/(2.*(1.+((CLBETA+1.)/(CLBETA-1.))**(2*ETAK-1.)))
 10   CONTINUE
C
C      
      a=1./(2.*clbeta)*alog((1+(exp(clbeta)-1.)*0.5)
     1  /(1+(exp(-clbeta)-1.)*0.5))
      do 11 k=NEVEN+1,nz-NEVEN-1
      ETAK=(K-NEVEN)/(NZ-2.*NEVEN)
      ZZ1(K)=0.5*(1.+sinh(clbeta*(etak-a))/sinh(clbeta*a))
 11   continue
      zz1(NEVEN)=0.
      zz1(nz-NEVEN)=1.
C
      TWOPI=2.*ACOS(-1.)
      DO 18 K=NEVEN,NZ-NEVEN
 18   DZ(K)=(1.-COS(ZZ1(K)*TWOPI))*(CLGAMA-1.)/2.+1.
C
      DO 20 K=1,NEVEN-1
 20   DZ(K)=DZ(NEVEN)
C
      DO 30 K=NZ-NEVEN+1,NZ
 30   DZ(K)=DZ(NEVEN)
C
      DZ(1)=DZ(1)/2.
      DZ(NZ-1)=DZ(NZ-1)/2.
      DZ(NZ)=DZ(NZ)/2.
C
      ZZ2(1)=0.
      DO 40 K=1,NZ
 40   ZZ2(K+1)=ZZ2(K)+DZ(K)
C
      DO 45 K=1,NZ+1
c 45   ZZ(K)=(ZZ2(K)-ZZ2(1))/(ZZ2(NZ)-ZZ2(1))*ZL-ZL
 45   ZZ(K)=(ZZ2(K)-ZZ2(1))/(ZZ2(NZ)-ZZ2(1))*ZL
C
      DO 50 K=1,NZ
 50   DZ(K)=ZZ(K+1)-ZZ(K)
C
      WRITE(29,299)
 299  FORMAT('VARIABLES=K,ZZ,DZ,DZRATIO,RDZ,RATIOMAX')
      WRITE(29,399)
 399  FORMAT('ZONE T=""')
C
      UNIT=1
      RMAX=DZ((NZ-1)/2)/DZ(1)
      DO 60 K=NZ,2,-1
       WRITE(29,99)K,ZZ(K),DZ(K),DZ(K)/DZ(K-1),1./DZ(K),RMAX
 99    FORMAT(I5,5E12.4)
 60   CONTINUE
      WRITE(29,99)1,ZZ(1),DZ(1),UNIT,1./DZ(1),RMAX
      
      DO K=1,NZ
         IF(K.EQ.1) THEN
            ZW(K)=ZZ(K)
         ELSE IF(K.EQ.NZ-1) THEN
            ZW(K)=ZZ(K+1)
         ELSE IF(K.EQ.NZ) THEN
            ZW(K)=ZZ(K+1)+DZ(K)
         ELSE
            ZW(K)=(ZZ(K)+ZZ(K+1))/2.
         ENDIF
      ENDDO
C
      RETURN
      END

C=====SUBROUTINE CLUSTER_NL END HERE





C============================================================================
      SUBROUTINE PRESSURE(NXMAX,NYMAX,NZMAX,NXMOD,NYMOD,NZ,PEX,PEY,DZ,
     1     ZZ,WORK,TRIGSX,TRIGSY,IFAX,PT,PP,HBAR,HH,HX,RE,DT)

      IMPLICIT NONE

      INTEGER I,J,K
      INTEGER NXMAX,NYMAX,NZMAX,NXMOD,NYMOD,NZ,NCPU

      REAL PEX,PEY,HBAR,RE,DT

      INCLUDE "mpif.h"
      
      INTEGER MYID,NUMPROCS,IERR

      REAL PP(NXMAX,NYMAX,*),PT(NXMAX,NYMAX,*)
      REAL HH(NXMAX,*),HX(NXMAX,*)
      REAL EXR(NXMAX,NYMAX),EYR(NXMAX,NYMAX)
      REAL HXR(NXMAX,NYMAX),HYR(NXMAX,NYMAX)
      REAL HER(NXMAX,NYMAX)
      REAL DZ(*),ZZ(*)
      REAL T1(NXMAX,NYMAX),T2(NXMAX,NYMAX)
      REAL T3(NXMAX,NYMAX),T4(NXMAX,NYMAX)
      REAL T5(NXMAX,NYMAX),T6(NXMAX,NYMAX)
      REAL LAP(NXMAX,NYMAX,NZMAX)
      REAL WORK(*),TRIGSX(*),TRIGSY(*)
      REAL FTMP(NXMAX,NYMAX)
      INTEGER IFAX(*)

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

C------------------------
C     NONLINEAR COEF
C------------------------

      DO I=1,NXMOD
         DO J=1,NYMAX
            HER(I,J)=1./(HH(I,J)+HBAR)
         ENDDO
      ENDDO

      DO I=1,NXMOD
         DO J=1,NYMAX
            EXR(I,J)=HX(I,J)*HER(I,J)
            EYR(I,J)=0.
            HXR(I,J)=HX(I,J)*HER(I,J)
            HYR(I,J)=0.
         ENDDO
      ENDDO

C-----------------
C     LAP_PHI
C-----------------

C~~~~~~~~~~~~~~
C     K=2
C~~~~~~~~~~~~~~

c sudin_edit
c      CALL PDFY(PT(1,1,3),T1,FTMP,WORK,TRIGSY,IFAX,PEY,NXMOD,NYMOD,
c     1     NXMAX,NYMAX)
c      CALL PDFY(PT(1,1,2),T2,FTMP,WORK,TRIGSY,IFAX,PEY,NXMOD,NYMOD,
c     1     NXMAX,NYMAX)
c sudin_edit

      CALL PDFY(PT(1,1,3),T1,FTMP,WORK,TRIGSY,IFAX,PEY,NXMOD,NYMOD,
     1     NXMAX)
      CALL PDFY(PT(1,1,2),T2,FTMP,WORK,TRIGSY,IFAX,PEY,NXMOD,NYMOD,
     1     NXMAX)

      DO I=1,NXMOD
         DO J=1,NYMAX
            T3(I,J)=(HYR(I,J)-ZZ(3)*EYR(I,J))*(PT(I,J,4)-PT(I,J,2))/2.
     1           /DZ(2)
            T4(I,J)=(HYR(I,J)-ZZ(2)*EYR(I,J))*(-PT(I,J,4)+4.*PT(I,J,3)
     1           -3.*PT(I,J,2))/2./DZ(2)
            T5(I,J)=(T1(I,J)+3.*T2(I,J)+T3(I,J)+3.*T4(I,J))/3./DZ(2)
         ENDDO
      ENDDO

      DO I=1,NXMOD
         DO J=1,NYMAX
            T5(I,J)=(HYR(I,J)-ZZ(2)*EYR(I,J))*T5(I,J)
         ENDDO
      ENDDO

C                                                                                                                  C
c Sudin_edit     
c      CALL PDFX(PT(1,1,3),T1,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
c     1     NXMAX,NYMAX)
c      CALL PDFX(PT(1,1,2),T2,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
c     1     NXMAX,NYMAX)
c Sudin_edit    

      CALL PDFX(PT(1,1,3),T1,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
     1     NYMOD,NXMAX)
      CALL PDFX(PT(1,1,2),T2,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
     1     NYMOD,NXMAX)
          
      
      DO I=1,NXMAX
         DO J=1,NYMAX
            T3(I,J)=(HXR(I,J)-ZZ(3)*EXR(I,J))*(PT(I,J,4)-PT(I,J,2))/2.
     1           /DZ(2)
            T4(I,J)=(HXR(I,J)-ZZ(2)*EXR(I,J))*(-PT(I,J,4)+4.*PT(I,J,3)
     1           -3.*PT(I,J,2))/2./DZ(2)
            T4(I,J)=(T1(I,J)+3.*T2(I,J)+T3(I,J)+3.*T4(I,J))/3./DZ(2)
         ENDDO
      ENDDO
      
      DO I=1,NXMOD
         DO J=1,NYMAX
            T4(I,J)=(HXR(I,J)-ZZ(2)*EXR(I,J))*T4(I,J)
         ENDDO
      ENDDO

CC

      DO I=1,NXMOD
         DO J=1,NYMAX
            T1(I,J)=(-PT(I,J,4)+4.*PT(I,J,3)-3.*PT(I,J,2))/2./DZ(2)
            T2(I,J)=(HXR(I,J)-ZZ(2)*EXR(I,J))*T1(I,J)
            T3(I,J)=(HYR(I,J)-ZZ(2)*EYR(I,J))*T1(I,J)
         ENDDO
      ENDDO

c Sudin_edit
c      CALL PDFX_(T2,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
c     1     NMAX,NYMAX)
c      CALL PDFY_(T3,FTMP,WORK,TRIGSY,IFAX,PEY,NXMOD,NYMOD,
c     1     NXMAX,NYMAX)
c Sudin_edit     

      CALL PDFX_(T2,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
     1     NYMOD,NXMAX)
      CALL PDFY_(T3,FTMP,WORK,TRIGSY,IFAX,PEY,NXMOD,NYMOD,
     1     NXMAX)



CC

      DO I=1,NXMOD
         DO J=1,NYMAX
            T1(I,J)=HER(I,J)*(PT(I,J,3)-PT(I,J,2))/DZ(2)**2
         ENDDO
      ENDDO

      DO I=1,NXMOD
         DO J=1,NYMAX
            T1(I,J)=HER(I,J)*T1(I,J)
            LAP(I,J,2)=T1(I,J)+T2(I,J)+T3(I,J)+T4(I,J)+T5(I,J)
         ENDDO
      ENDDO


C Sudin_edit
c      CALL PDFXX(PT(1,1,2),T1,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
c     1     NXMAX,NYMAX)
c      CALL PDFYY(PT(1,1,2),T2,FTMP,WORK,TRIGSY,IFAX,PEY,NXMOD,NYMOD,
c     1     NXMAX,NYMAX)
C Sudin_edit

      CALL PDFXX(PT(1,1,2),T1,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
     1     NYMOD,NXMAX)
      CALL PDFYY(PT(1,1,2),T2,FTMP,WORK,TRIGSY,IFAX,PEY,NXMOD,NYMOD,
     1     NXMAX)

      DO I=1,NXMOD
         DO J=1,NYMAX
            LAP(I,J,2)=LAP(I,J,2)+T1(I,J)+T2(I,J)
         ENDDO
      ENDDO

C~~~~~END HERE

C~~~~~~~~~~~~~~
C     K=3
C~~~~~~~~~~~~~~

      CALL PDFY(PT(1,1,4),T1,FTMP,WORK,TRIGSY,IFAX,PEY,NXMOD,NYMOD,
     1     NXMAX)
      CALL PDFY(PT(1,1,2),T2,FTMP,WORK,TRIGSY,IFAX,PEY,NXMOD,NYMOD,
     1     NXMAX)
      
      DO I=1,NXMOD
         DO J=1,NYMAX
            T3(I,J)=(HYR(I,J)-ZZ(4)*EYR(I,J))*(PT(I,J,5)-PT(I,J,3))
     1           /(DZ(3)+DZ(4))
            T4(I,J)=(HYR(I,J)-ZZ(2)*EYR(I,J))*(-PT(I,J,4)+4.*PT(I,J,3)
     1           -3.*PT(I,J,2))/2./DZ(2)
            T5(I,J)=(T1(I,J)-T2(I,J)+T3(I,J)-T4(I,J))/2./DZ(2)           
         ENDDO
      ENDDO
      
      DO I=1,NXMOD
         DO J=1,NYMAX
            T5(I,J)=(HYR(I,J)-ZZ(3)*EYR(I,J))*T5(I,J)
         ENDDO
      ENDDO

CC
c Sudin_edit
c      CALL PDFX(PT(1,1,4),T1,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
c     1     NXMAX,NYMAX)
c      CALL PDFX(PT(1,1,2),T2,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
c     1     NXMAX,NYMAX)
c Sudin_edit     

      CALL PDFX(PT(1,1,4),T1,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
     1     NYMOD,NXMAX)
      CALL PDFX(PT(1,1,2),T2,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
     1     NYMOD,NXMAX)
      
      DO I=1,NXMOD
         DO J=1,NYMAX
            T3(I,J)=(HXR(I,J)-ZZ(4)*EXR(I,J))*(PT(I,J,5)-PT(I,J,3))
     1           /(DZ(3)+DZ(4))
            T4(I,J)=(HXR(I,J)-ZZ(2)*EXR(I,J))*(-PT(I,J,4)+4.*PT(I,J,3)
     1           -3.*PT(I,J,2))/2./DZ(2)
            T4(I,J)=(T1(I,J)-T2(I,J)+T3(I,J)-T4(I,J))/2./DZ(2)           
         ENDDO
      ENDDO
      
      DO I=1,NXMOD
         DO J=1,NYMAX
            T4(I,J)=(HXR(I,J)-ZZ(3)*EXR(I,J))*T4(I,J)
         ENDDO
      ENDDO

CC

      DO I=1,NXMOD
         DO J=1,NYMAX
            T1(I,J)=(PT(I,J,4)-PT(I,J,2))/2./DZ(2)
            T2(I,J)=(HXR(I,J)-ZZ(3)*EXR(I,J))*T1(I,J)
            T3(I,J)=(HYR(I,J)-ZZ(3)*EYR(I,J))*T1(I,J)
         ENDDO
      ENDDO

C Sudin_edit
c      CALL PDFX_(T2,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
c     1     NXMAX,NYMAX)
c      CALL PDFY_(T3,FTMP,WORK,TRIGSY,IFAX,PEY,NXMOD,NYMOD,
c     1     NXMAX,NYMAX)
C Sudin_edit     
 
      CALL PDFX_(T2,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
     1     NYMOD,NXMAX)
      CALL PDFY_(T3,FTMP,WORK,TRIGSY,IFAX,PEY,NXMOD,NYMOD,
     1     NXMAX)
      
CC

      DO I=1,NXMOD
         DO J=1,NYMAX
            T1(I,J)=HER(I,J)*(PT(I,J,4)-2.*PT(I,J,3)+PT(I,J,2))
     1           /(DZ(2))**2
         ENDDO
      ENDDO
      
      DO I=1,NXMOD
         DO J=1,NYMAX
            T1(I,J)=HER(I,J)*T1(I,J)
            LAP(I,J,3)=T1(I,J)+T2(I,J)+T3(I,J)+T4(I,J)+T5(I,J)
         ENDDO
      ENDDO

c Sudin_edit
c      CALL PDFXX(PT(1,1,3),T1,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
c     1     NXMAX,NYMAX)
c      CALL PDFYY(PT(1,1,3),T2,FTMP,WORK,TRIGSY,IFAX,PEY,NXMOD,NYMOD,
c     1     NXMAX,NYMAX)
c Sudin_edit

      CALL PDFXX(PT(1,1,3),T1,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
     1     NYMOD,NXMAX)
      CALL PDFYY(PT(1,1,3),T2,FTMP,WORK,TRIGSY,IFAX,PEY,NXMOD,NYMOD,
     1     NXMAX)

      DO I=1,NXMOD
         DO J=1,NYMAX
            LAP(I,J,3)=LAP(I,J,3)+T1(I,J)+T2(I,J)
         ENDDO
      ENDDO

C~~~~~END HERE

C~~~~~~~~~~~~~~~~~~~~
C     K=4,NZ-3
C~~~~~~~~~~~~~~~~~~~~

      DO 100 K=4,NZ-3

         CALL PDFY(PT(1,1,K+1),T1,FTMP,WORK,TRIGSY,IFAX,PEY,NXMOD,
     1        NYMOD,NXMAX)
         CALL PDFY(PT(1,1,K-1),T2,FTMP,WORK,TRIGSY,IFAX,PEY,NXMOD,
     1        NYMOD,NXMAX)
         
         DO I=1,NXMOD
            DO J=1,NYMAX
               T3(I,J)=(HYR(I,J)-ZZ(K+1)*EYR(I,J))
     1              *(PT(I,J,K+2)-PT(I,J,K))/(DZ(K)+DZ(K+1))
               T4(I,J)=(HYR(I,J)-ZZ(K-1)*EYR(I,J))
     1              *(PT(I,J,K)-PT(I,J,K-2))/(DZ(K-2)+DZ(K-1))
               T5(I,J)=(T1(I,J)-T2(I,J)+T3(I,J)-T4(I,J))/(DZ(K-1)+DZ(K))
            ENDDO
         ENDDO
         
         DO I=1,NXMOD
            DO J=1,NYMAX
               T5(I,J)=(HYR(I,J)-ZZ(K)*EYR(I,J))*T5(I,J)
            ENDDO
         ENDDO
         
CC

C sudin_edit
c         CALL PDFX(PT(1,1,K+1),T1,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
c     1        NXMAX,NYMAX)
c         CALL PDFX(PT(1,1,K-1),T2,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
c     1        NXMAX,NYMAX)
C sudin_edit      
         CALL PDFX(PT(1,1,K+1),T1,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
     1        NYMOD,NXMAX)
         CALL PDFX(PT(1,1,K-1),T2,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
     1        NYMOD,NXMAX)
         
         DO I=1,NXMOD
            DO J=1,NYMAX
               T3(I,J)=(HXR(I,J)-ZZ(K+1)*EXR(I,J))
     1              *(PT(I,J,K+2)-PT(I,J,K))/(DZ(K)+DZ(K+1))
               T4(I,J)=(HXR(I,J)-ZZ(K-1)*EXR(I,J))
     1              *(PT(I,J,K)-PT(I,J,K-2))/(DZ(K-2)+DZ(K-1))
               T4(I,J)=(T1(I,J)-T2(I,J)+T3(I,J)-T4(I,J))/(DZ(K-1)+DZ(K)) 
            ENDDO
         ENDDO
         
         DO I=1,NXMOD
            DO J=1,NYMAX
               T4(I,J)=(HXR(I,J)-ZZ(K)*EXR(I,J))*T4(I,J)
            ENDDO
         ENDDO
         
CC
         
         DO I=1,NXMOD
            DO J=1,NYMAX
               T1(I,J)=(PT(I,J,K+1)-PT(I,J,K-1))/(DZ(K-1)+DZ(K))
               T2(I,J)=(HXR(I,J)-ZZ(K)*EXR(I,J))*T1(I,J)
               T3(I,J)=(HYR(I,J)-ZZ(K)*EYR(I,J))*T1(I,J)
            ENDDO
         ENDDO


C sudin_edit
c         CALL PDFX_(T2,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
c     1        NXMAX,NYMAX)
c         CALL PDFY_(T3,FTMP,WORK,TRIGSY,IFAX,PEY,NXMOD,NYMOD,
c     1        NXMAX,NYMAX)
C sudin_edit
         CALL PDFX_(T2,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
     1        NYMOD,NXMAX)
         CALL PDFY_(T3,FTMP,WORK,TRIGSY,IFAX,PEY,NXMOD,NYMOD,
     1        NXMAX)
         
CC

         DO I=1,NXMOD
            DO J=1,NYMAX
               T1(I,J)=HER(I,J)*((PT(I,J,K+1)-PT(I,J,K))/DZ(K)
     1              -(PT(I,J,K)-PT(I,J,K-1))/DZ(K-1))*2./(DZ(K-1)+DZ(K))
            ENDDO
         ENDDO
         
         DO I=1,NXMOD
            DO J=1,NYMAX
               T1(I,J)=HER(I,J)*T1(I,J)
               LAP(I,J,K)=T1(I,J)+T2(I,J)+T3(I,J)+T4(I,J)+T5(I,J)
            ENDDO
         ENDDO

c Sudin_edit
c         CALL PDFXX(PT(1,1,K),T1,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
c     1        NXMAX,NYMAX)
c         CALL PDFYY(PT(1,1,K),T2,FTMP,WORK,TRIGSY,IFAX,PEY,NXMOD,NYMOD,
c     1        NXMAX,NYMAX)
C Sudin_edit         
         CALL PDFXX(PT(1,1,K),T1,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
     1        NYMOD,NXMAX)
         CALL PDFYY(PT(1,1,K),T2,FTMP,WORK,TRIGSY,IFAX,PEY,NXMOD,NYMOD,
     1        NXMAX)
         
         DO I=1,NXMOD
            DO J=1,NYMAX
               LAP(I,J,K)=LAP(I,J,K)+T1(I,J)+T2(I,J)
            ENDDO
         ENDDO

 100  CONTINUE

C~~~~~END HERE

C~~~~~~~~~~~~~~~~~
C     K=NZ-2
C~~~~~~~~~~~~~~~~~

      CALL PDFY(PT(1,1,NZ-1),T1,FTMP,WORK,TRIGSY,IFAX,PEY,NXMOD,
     1     NYMOD,NXMAX)
      CALL PDFY(PT(1,1,NZ-3),T2,FTMP,WORK,TRIGSY,IFAX,PEY,NXMOD,
     1     NYMOD,NXMAX)
      
      DO I=1,NXMOD
         DO J=1,NYMAX
            T3(I,J)=(HYR(I,J)-ZZ(NZ-1)*EYR(I,J))*(4.*PT(I,J,NZ)
     1           -3.*PT(I,J,NZ-1)-PT(I,J,NZ-2))/3./DZ(NZ-2)
            T4(I,J)=(HYR(I,J)-ZZ(NZ-3)*EYR(I,J))*(PT(I,J,NZ-2)
     1           -PT(I,J,NZ-4))/(DZ(NZ-4)+DZ(NZ-3))
            T5(I,J)=(T1(I,J)-T2(I,J)+T3(I,J)-T4(I,J))/2./DZ(NZ-2)
         ENDDO
      ENDDO
      
      DO I=1,NXMOD
         DO J=1,NYMAX
            T5(I,J)=(HYR(I,J)-ZZ(NZ-2)*EYR(I,J))*T5(I,J)
         ENDDO
      ENDDO

CC

C sudin_edit
c      CALL PDFX(PT(1,1,NZ-1),T1,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
c     1     NXMAX,NYMAX)
c      CALL PDFX(PT(1,1,NZ-3),T2,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
c     1     NXMAX,NYMAX)
C sudin_edit

      
      CALL PDFX(PT(1,1,NZ-1),T1,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
     1     NYMOD,NXMAX)
      CALL PDFX(PT(1,1,NZ-3),T2,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
     1     NYMOD,NXMAX)
      
      DO I=1,NXMOD
         DO J=1,NYMAX
            T3(I,J)=(HXR(I,J)-ZZ(NZ-1)*EXR(I,J))*(4.*PT(I,J,NZ)
     1           -3.*PT(I,J,NZ-1)-PT(I,J,NZ-2))/3./DZ(NZ-2)
            T4(I,J)=(HXR(I,J)-ZZ(NZ-3)*EXR(I,J))*(PT(I,J,NZ-2)
     1           -PT(I,J,NZ-4))/(DZ(NZ-4)+DZ(NZ-3))
            T4(I,J)=(T1(I,J)-T2(I,J)+T3(I,J)-T4(I,J))/2./DZ(NZ-2)           
         ENDDO
      ENDDO
      
      DO I=1,NXMOD
         DO J=1,NYMAX
            T4(I,J)=(HXR(I,J)-ZZ(NZ-2)*EXR(I,J))*T4(I,J)
         ENDDO
      ENDDO

CC

      DO I=1,NXMOD
         DO J=1,NYMAX
            T1(I,J)=(PT(I,J,NZ-1)-PT(I,J,NZ-3))/2./DZ(NZ-2)
            T2(I,J)=(HXR(I,J)-ZZ(NZ-2)*EXR(I,J))*T1(I,J)
            T3(I,J)=(HYR(I,J)-ZZ(NZ-2)*EYR(I,J))*T1(I,J)
         ENDDO
      ENDDO

C Sudin_edit
c      CALL PDFX_(T2,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
c     1     NXMAX,NYMAX)
c      CALL PDFY_(T3,FTMP,WORK,TRIGSY,IFAX,PEY,NXMOD,NYMOD,
c     1     NXMAX,NYMAX)
C Sudin_edit

      CALL PDFX_(T2,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
     1     NYMOD,NXMAX)
      CALL PDFY_(T3,FTMP,WORK,TRIGSY,IFAX,PEY,NXMOD,NYMOD,
     1     NXMAX)
      
CC

      DO I=1,NXMOD
         DO J=1,NYMAX
            T1(I,J)=HER(I,J)*(PT(I,J,NZ-1)-2.*PT(I,J,NZ-2)+PT(I,J,NZ-3))
     1           /DZ(NZ-2)**2
         ENDDO
      ENDDO

      DO I=1,NXMOD
         DO J=1,NYMAX
            T1(I,J)=HER(I,J)*T1(I,J)
            LAP(I,J,NZ-2)=T1(I,J)+T2(I,J)+T3(I,J)+T4(I,J)+T5(I,J)
         ENDDO
      ENDDO

c Sudin_edit
c      CALL PDFXX(PT(1,1,NZ-2),T1,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
c    1     NXMAX,NYMAX)
c      CALL PDFYY(PT(1,1,NZ-2),T2,FTMP,WORK,TRIGSY,IFAX,PEY,NXMOD,NYMOD,
c     1     NXMAX,NYMAX)
C Sudin_edit      

      CALL PDFXX(PT(1,1,NZ-2),T1,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
     1     NYMOD,NXMAX)
      CALL PDFYY(PT(1,1,NZ-2),T2,FTMP,WORK,TRIGSY,IFAX,PEY,NXMOD,NYMOD,
     1     NXMAX)
      
      DO I=1,NXMOD
         DO J=1,NYMAX
            LAP(I,J,NZ-2)=LAP(I,J,NZ-2)+T1(I,J)+T2(I,J)
         ENDDO
      ENDDO

C~~~~~END HERE

C~~~~~~~~~~~~~~~~~~
C     K=NZ-1
C~~~~~~~~~~~~~~~~~~

      CALL PDFY(PT(1,1,NZ),T1,FTMP,WORK,TRIGSY,IFAX,PEY,NXMOD,
     1     NYMOD,NXMAX)
      CALL PDFY(PT(1,1,NZ-1),T2,FTMP,WORK,TRIGSY,IFAX,PEY,NXMOD,
     1     NYMOD,NXMAX)
      CALL PDFY(PT(1,1,NZ-2),T3,FTMP,WORK,TRIGSY,IFAX,PEY,NXMOD,
     1     NYMOD,NXMAX)
      
      DO I=1,NXMOD
         DO J=1,NYMAX
            T4(I,J)=4.*(HYR(I,J)-ZZ(NZ)*EYR(I,J))*(8.*PT(I,J,NZ)-9.
     1           *PT(I,J,NZ-1)+PT(I,J,NZ-2))/3./DZ(NZ-2)-3.*(HYR(I,J)
     1           -ZZ(NZ-1)*EYR(I,J))*(4.*PT(I,J,NZ)-3.*PT(I,J,NZ-1)
     1           -PT(I,J,NZ-2))/3./DZ(NZ-2)-(HYR(I,J)-ZZ(NZ-2)
     1           *EYR(I,J))*(PT(I,J,NZ-1)-PT(I,J,NZ-3))
     1           /(DZ(NZ-2)+DZ(NZ-3))
            T5(I,J)=(4.*T1(I,J)-3.*T2(I,J)-T3(I,J)+T4(I,J))/3./DZ(NZ-2)
         ENDDO
      ENDDO
      
      DO I=1,NXMOD
         DO J=1,NYMAX
            T5(I,J)=(HYR(I,J)-ZZ(NZ-1)*EYR(I,J))*T5(I,J)
         ENDDO
      ENDDO

CC

c Sudin_edit
c      CALL PDFX(PT(1,1,NZ),T1,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
c     1     NXMAX,NYMAX)
c      CALL PDFX(PT(1,1,NZ-1),T2,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
c     1     NXMAX,NYMAX)
c      CALL PDFX(PT(1,1,NZ-2),T3,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
c     1     NXMAX,NYMAX)
c Sudin_edit

      CALL PDFX(PT(1,1,NZ),T1,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
     1     NYMOD,NXMAX)
      CALL PDFX(PT(1,1,NZ-1),T2,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
     1     NYMOD,NXMAX)
      CALL PDFX(PT(1,1,NZ-2),T3,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
     1     NYMOD,NXMAX)
      
      DO I=1,NXMOD
         DO J=1,NYMAX
            T4(I,J)=4.*(HXR(I,J)-ZZ(NZ)*EXR(I,J))*(8.*PT(I,J,NZ)-9.
     1           *PT(I,J,NZ-1)+PT(I,J,NZ-2))/3./DZ(NZ-2)-3.*(HXR(I,J)
     1           -ZZ(NZ-1)*EXR(I,J))*(4.*PT(I,J,NZ)-3.*PT(I,J,NZ-1)
     1           -PT(I,J,NZ-2))/3./DZ(NZ-2)-(HXR(I,J)-ZZ(NZ-2)
     1           *EXR(I,J))*(PT(I,J,NZ-1)-PT(I,J,NZ-3))
     1           /(DZ(NZ-2)+DZ(NZ-3))
            T4(I,J)=(4.*T1(I,J)-3.*T2(I,J)-T3(I,J)+T4(I,J))/3./DZ(NZ-2)
         ENDDO
      ENDDO
      
      DO I=1,NXMOD
         DO J=1,NYMAX
            T4(I,J)=(HXR(I,J)-ZZ(NZ-1)*EXR(I,J))*T4(I,J)
         ENDDO
      ENDDO

CC

      DO I=1,NXMOD
         DO J=1,NYMAX
            T1(I,J)=(4.*PT(I,J,NZ)-3.*PT(I,J,NZ-1)-PT(I,J,NZ-2))/3.
     1           /DZ(NZ-2)
            T2(I,J)=(HXR(I,J)-ZZ(NZ-1)*EXR(I,J))*T1(I,J)
            T3(I,J)=(HYR(I,J)-ZZ(NZ-1)*EYR(I,J))*T1(I,J)
         ENDDO
      ENDDO

c Sudin_edit
c      CALL PDFX_(T2,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
c     1     NXMAX,NYMAX)
c      CALL PDFY_(T3,FTMP,WORK,TRIGSY,IFAX,PEY,NXMOD,NYMOD,
c     1     NXMAX,NYMAX)
c Sudin_edit

      CALL PDFX_(T2,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
     1     NYMOD,NXMAX)
      CALL PDFY_(T3,FTMP,WORK,TRIGSY,IFAX,PEY,NXMOD,NYMOD,
     1     NXMAX)

CC

      DO I=1,NXMOD
         DO J=1,NYMAX
            T1(I,J)=HER(I,J)*(8.*PT(I,J,NZ)-12.*PT(I,J,NZ-1)+4.
     1           *PT(I,J,NZ-2))/3./DZ(NZ-2)**2
         ENDDO
      ENDDO

      DO I=1,NXMOD
         DO J=1,NYMAX
            T1(I,J)=HER(I,J)*T1(I,J)
            LAP(I,J,NZ-1)=T1(I,J)+T2(I,J)+T3(I,J)+T4(I,J)+T5(I,J)
         ENDDO
      ENDDO

c Sudin_edit
c      CALL PDFXX(PT(1,1,NZ-1),T1,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
c     1     NXMAX,NYMAX)
c      CALL PDFYY(PT(1,1,NZ-1),T2,FTMP,WORK,TRIGSY,IFAX,PEY,NXMOD,NYMOD,
c     1     NXMAX,NYMAX)
c Sudin_edit

      CALL PDFXX(PT(1,1,NZ-1),T1,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,
     1     NYMOD,NXMAX)
      CALL PDFYY(PT(1,1,NZ-1),T2,FTMP,WORK,TRIGSY,IFAX,PEY,NXMOD,NYMOD,
     1     NXMAX)
      
      DO I=1,NXMOD
         DO J=1,NYMAX
            LAP(I,J,NZ-1)=LAP(I,J,NZ-1)+T1(I,J)+T2(I,J)
         ENDDO
      ENDDO

C~~~~~END HERE

C-----END HERE

C--------------------
C     PRESSURE
C--------------------

      DO I=1,NXMOD
         DO J=1,NYMOD
            LAP(I,J,1)=(35.*LAP(I,J,2)-35.*LAP(I,J,3)+21.*LAP(I,J,4)
     1           -5.*LAP(I,J,5))/16.
            LAP(I,J,NZ)=(35.*LAP(I,J,NZ-1)-35.*LAP(I,J,NZ-2)
     1           +21.*LAP(I,J,NZ-3)-5.*LAP(I,J,NZ-4))/16.
         ENDDO
      ENDDO

      DO K=1,NZ
         DO J=1,NYMOD
            DO I=1,NXMOD
               PP(I,J,K)=PT(I,J,K)-DT/(2.*RE)*LAP(I,J,K)
            ENDDO
         ENDDO
      ENDDO

C-----END HERE

      RETURN
      END

C=====SUBROUTINE PRESSURE END HERE





C=============================================================================
      SUBROUTINE AVERAGE_PSI_2D(NXMOD,NYMOD,NXMAX,NYMAX,F,FMY)

C     AVERAGE F(I,J) W.R.T. PSI
c This subroutine computes the average of a 2D field F(x, y) over the y-dimension, resulting in a 1D field FMY(x)

      IMPLICIT NONE
      
      INTEGER NXMOD,NYMOD,NXMAX,NYMAX,I,J

      REAL F(NXMAX,*),FMY(*)
      
      DO I=1,NXMOD
         FMY(I)=0.
         DO J=1,NYMOD
            FMY(I)=FMY(I)+F(I,J)
         ENDDO
         FMY(I)=FMY(I)/FLOAT(NYMOD)
      ENDDO
      
      RETURN
      END

C=====SUBROUTINE AVERAGE_PSI_2D END HERE





C=============================================================================
      SUBROUTINE AVERAGE_PSI_3D(NXMOD,NYMOD,NZ,NXMAX,NYMAX,NZMAX,F,FMY)

c Given a 3D field F(x, y, z), this subroutine computes the average over the y-dimension at each (x, z) location, resulting in a 2D field FMY(x, z)
c Converts 3D  domain into a vertical X-Z plane averagerd slice
C     AVERAGE F(I,J,K) W.R.T. PSI

      IMPLICIT NONE
      
      INTEGER NXMOD,NYMOD,NXMAX,NYMAX,NZMAX,NZ,NZ1,I,J,K

      REAL F(NXMAX,NYMAX,*),FMY(NXMAX,*)

      NZ1=NZ+1
      
      DO I=1,NXMOD
         DO K=1,NZ1
            FMY(I,K)=0.
            DO J=1,NYMOD
               FMY(I,K)=FMY(I,K)+F(I,J,K)
            ENDDO
            FMY(I,K)=FMY(I,K)/FLOAT(NYMOD)
         ENDDO
      ENDDO

      RETURN
      END

C=====SUBROUTINE AVERAGE_PSI_3D END HERE





C=============================================================================
      SUBROUTINE PHASE_SHIFT_1D(NXMOD,NXMAX,NTIME,FMY,FSHIFT,DTS,HOMEG,
     1     WORK,TRIGSX,IFAX)

      IMPLICIT NONE

      INTEGER NXMOD,NXMAX,NTIME,I,IMOD
      REAL HOMEG,DTS,THETA

      REAL FMY(*),FSHIFT(*)
      REAL TMP(NXMAX),ETMP(NXMAX)
      REAL WORK(*),TRIGSX(*)
      INTEGER IFAX(*)

C-----------------------------
C     SHIFT TO SAME PHASE
C-----------------------------
      
      THETA=-HOMEG*DTS

      DO I=1,NXMOD
         TMP(I)=FMY(I)
      ENDDO

      CALL FFTX(TMP,WORK,TRIGSX,IFAX,NXMOD,1,NXMAX,-1)

      DO I=1,NXMOD-1,2
         IMOD=(I-1)/2
         ETMP(I)=TMP(I)*COS(IMOD*THETA)-TMP(I+1)*SIN(IMOD*THETA)
         ETMP(I+1)=TMP(I)*SIN(IMOD*THETA)+TMP(I+1)*COS(IMOD*THETA)
      ENDDO
      ETMP(NXMOD+1)=TMP(NXMOD+1)
      ETMP(NXMOD+2)=TMP(NXMOD+2)

      CALL FFTX(ETMP,WORK,TRIGSX,IFAX,NXMOD,1,NXMAX,1)

C-----END HERE

C------------------------
C     PHASE AVERAGE
C------------------------

      DO I=1,NXMOD
         FSHIFT(I)=ETMP(I)
      ENDDO
      
C-----END HERE

      RETURN
      END

C=====SUBROUTINE PHASE_SHIFT_1D END HERE





C=============================================================================
      SUBROUTINE PHASE_SHIFT_2D(NXMOD,NZ,NXMAX,NZMAX,NTIME,FMY,FSHIFT,
     1     DTS,HOMEG,WORK,TRIGSX,IFAX)

      IMPLICIT NONE

      INTEGER NXMOD,NXMAX,NZMAX,NTIME,I,K,NZ,NZ1,IMOD
      REAL HOMEG,DTS,THETA

      REAL FMY(NXMAX,*),FSHIFT(NXMAX,*)
      REAL TMP(NXMAX),ETMP(NXMAX)
      REAL WORK(*),TRIGSX(*)
      INTEGER IFAX(*)

      NZ1=NZ+1
      
c you're multiplying a rate of phase rotation (in radians per second) by a dimensionless 
c time shift in domain lengths — but with the implicit understanding that TIME0 has been scaled to wavelength-units.    
c HOMEG is the rate of phase change   
c “If I know how many domain lengths the wave has traveled (DTS), how much total phase has accumulated?”
c θ=total accumulated phase shift (in radians) across the domain

C---------------------------------------------------------------------
C PHASE SHIFT EXPLANATION:
C - One wave cycle corresponds to 2π radians of phase.
C - If the bottom boundary contains NWAVE = 6 full wave cycles,
C   then the total spatial phase span of the domain is:
C       total_phase = 6 × 2π = 12π radians.
C - The wave travels during time0 units; normalized to the domain
C   this gives a fractional displacement:
C       DTS = -TIME0 / NWAVE
C   ⇒ negative indicates backward shift for phase alignment.
C - The total phase that must be traced back to align with the
C   reference phase is:
C       THETA = -OMEGA × DTS
C - For example: if DTS = -100.2, then the wave has moved forward
C   by 100.2 domain lengths, and the corresponding phase shift is:
C       THETA = 2π × 6 × 100.2 = 1204.8π radians
C   ⇒ Each Fourier mode m is then rotated by m × THETA radians
C      to restore the wave to reference alignment.
C---------------------------------------------------------------------

      THETA=-HOMEG*DTS


      DO K=1,NZ1

C-----------------------------
C     SHIFT TO SAME PHASE
C-----------------------------
      
         DO I=1,NXMOD
            TMP(I)=FMY(I,K)
         ENDDO

         CALL FFTX(TMP,WORK,TRIGSX,IFAX,NXMOD,1,NXMAX,-1)

         DO I=1,NXMOD-1,2
            IMOD=(I-1)/2
            ETMP(I)=TMP(I)*COS(IMOD*THETA)-TMP(I+1)*SIN(IMOD*THETA)
            ETMP(I+1)=TMP(I)*SIN(IMOD*THETA)+TMP(I+1)*COS(IMOD*THETA)
         ENDDO
         ETMP(NXMOD+1)=TMP(NXMOD+1)
         ETMP(NXMOD+2)=TMP(NXMOD+2)

         CALL FFTX(ETMP,WORK,TRIGSX,IFAX,NXMOD,1,NXMAX,1)

C-----END HERE

C------------------------
C     PHASE AVERAGE
C------------------------

         DO I=1,NXMOD
            FSHIFT(I,K)=ETMP(I)
         ENDDO
      
C-----END HERE

      ENDDO

      RETURN
      END

C=====SUBROUTINE PHASE_SHIFT_2D END HERE


C=============================================================================
      SUBROUTINE PHASE_SHIFT_3D(NXMOD,NYMOD,NZ,NXMAX,NYMAX,NZMAX,NTIME,FMY,FSHIFT,
     1     DTS,HOMEG,WORK,TRIGSX,IFAX)

      IMPLICIT NONE

      INTEGER NXMOD,NYMOD,NXMAX,NYMAX,NZMAX,NTIME,I,J,K,NZ,NZ1,IMOD
      REAL HOMEG,DTS,THETA

      REAL FMY(NXMAX,NYMAX,*),FSHIFT(NXMAX,NYMAX,*)
      REAL TMP(NXMAX),ETMP(NXMAX)
      REAL WORK(*),TRIGSX(*)
      INTEGER IFAX(*)

      NZ1=NZ+1
      
c you're multiplying a rate of phase rotation (in radians per second) by a dimensionless 
c time shift in domain lengths — but with the implicit understanding that TIME0 has been scaled to wavelength-units.    
c HOMEG is the rate of phase change   
c “If I know how many domain lengths the wave has traveled (DTS), how much total phase has accumulated?”
c θ=total accumulated phase shift (in radians) across the domain

C---------------------------------------------------------------------
C PHASE SHIFT EXPLANATION:
C - One wave cycle corresponds to 2π radians of phase.
C - If the bottom boundary contains NWAVE = 6 full wave cycles,
C   then the total spatial phase span of the domain is:
C       total_phase = 6 × 2π = 12π radians.
C - The wave travels during time0 units; normalized to the domain
C   this gives a fractional displacement:
C       DTS = -TIME0 / NWAVE
C   ⇒ negative indicates backward shift for phase alignment.
C - The total phase that must be traced back to align with the
C   reference phase is:
C       THETA = -OMEGA × DTS
C - For example: if DTS = -100.2, then the wave has moved forward
C   by 100.2 domain lengths, and the corresponding phase shift is:
C       THETA = 2π × 6 × 100.2 = 1204.8π radians
C   ⇒ Each Fourier mode m is then rotated by m × THETA radians
C      to restore the wave to reference alignment.
C---------------------------------------------------------------------

      THETA=-HOMEG*DTS
      DO I=1,NXMAX
         TMP(I)=0.0
         ETMP(I)=0.0
      ENDDO

C-----------------------------
C     INITIALIZE OUTPUT ARRAY
C-----------------------------
      DO K = 1, NZMAX
         DO J = 1, NYMAX
            DO I = 1, NXMAX
               FSHIFT(I,J,K) = 0.0
            ENDDO
         ENDDO
      ENDDO

      DO J = 1, NYMOD
         DO K=1,NZ1

C-----------------------------
C     SHIFT TO SAME PHASE
C-----------------------------
      
            DO I=1,NXMOD
               TMP(I)=FMY(I,J,K)
            ENDDO

            CALL FFTX(TMP,WORK,TRIGSX,IFAX,NXMOD,1,NXMAX,-1)

            DO I=1,NXMOD-1,2
               IMOD=(I-1)/2
               ETMP(I)=TMP(I)*COS(IMOD*THETA)-TMP(I+1)*SIN(IMOD*THETA)
               ETMP(I+1)=TMP(I)*SIN(IMOD*THETA)+TMP(I+1)*COS(IMOD*THETA)
            ENDDO
            ETMP(NXMOD+1)=TMP(NXMOD+1)
            ETMP(NXMOD+2)=TMP(NXMOD+2)

            CALL FFTX(ETMP,WORK,TRIGSX,IFAX,NXMOD,1,NXMAX,1)

C-----END HERE

C------------------------
C     PHASE AVERAGE
C------------------------

            DO I=1,NXMOD
               FSHIFT(I,j,K)=ETMP(I)
            ENDDO
         
         ENDDO
      
C-----END HERE

      ENDDO

      RETURN
      END

C=====SUBROUTINE PHASE_SHIFT_2D END HERE



C=============================================================================
      SUBROUTINE PHASE_AVERAGE_1D(NXMOD,NXMAX,NTIME,FMY,FPM,DTS,HOMEG,
     1     WORK,TRIGSX,IFAX)

      IMPLICIT NONE

      INTEGER NXMOD,NXMAX,NTIME,I,IMOD
      REAL HOMEG,DTS,THETA

      REAL FMY(*),FPM(*)
      REAL TMP(NXMAX),ETMP(NXMAX)
      REAL WORK(*),TRIGSX(*)
      INTEGER IFAX(*)

C-----------------------------
C     SHIFT TO SAME PHASE
C-----------------------------
      
      THETA=-HOMEG*DTS

      DO I=1,NXMOD
         TMP(I)=FMY(I)
      ENDDO

      CALL FFTX(TMP,WORK,TRIGSX,IFAX,NXMOD,1,NXMAX,-1)

      DO I=1,NXMOD-1,2
         IMOD=(I-1)/2
         ETMP(I)=TMP(I)*COS(IMOD*THETA)-TMP(I+1)*SIN(IMOD*THETA)
         ETMP(I+1)=TMP(I)*SIN(IMOD*THETA)+TMP(I+1)*COS(IMOD*THETA)
      ENDDO
      ETMP(NXMOD+1)=TMP(NXMOD+1)
      ETMP(NXMOD+2)=TMP(NXMOD+2)

      CALL FFTX(ETMP,WORK,TRIGSX,IFAX,NXMOD,1,NXMAX,1)

C-----END HERE

C------------------------
C     PHASE AVERAGE
C------------------------

      DO I=1,NXMOD
         FPM(I)=FPM(I)+ETMP(I)/FLOAT(NTIME)
      ENDDO
      
C-----END HERE

      RETURN
      END

C=====SUBROUTINE PHASE_AVERAGE_1D END HERE






C=============================================================================
      SUBROUTINE PHASE_AVERAGE_2D(NXMOD,NZ,NXMAX,NZMAX,NTIME,FMY,FPM,
     1     DTS,HOMEG,WORK,TRIGSX,IFAX)

      IMPLICIT NONE

      INTEGER NXMOD,NXMAX,NZMAX,NTIME,I,K,NZ,NZ1,IMOD
      REAL HOMEG,DTS,THETA

      REAL FMY(NXMAX,*),FPM(NXMAX,*)
      REAL TMP(NXMAX),ETMP(NXMAX)
      REAL WORK(*),TRIGSX(*)
      INTEGER IFAX(*)

      NZ1=NZ+1
      
      THETA=-HOMEG*DTS
      
      DO K=1,NZ1

C-----------------------------
C     SHIFT TO SAME PHASE
C-----------------------------
c Convert the 2D domain to       
         DO I=1,NXMOD
            TMP(I)=FMY(I,K)
         ENDDO

         CALL FFTX(TMP,WORK,TRIGSX,IFAX,NXMOD,1,NXMAX,-1)

C Loop goes from 1 to NXMOD-1 in steps of 2
C Rotating the complex number f = a + i·b by angle θ = m · THETA
C Using Euler's identity: f_rotated = (a + i·b) · (cos(mθ) - i·sin(mθ))
C Real part:      f_real = a · cos(mθ) - b · sin(mθ)
C Imaginary part: f_imag = a · sin(mθ) + b · cos(mθ)

         DO I=1,NXMOD-1,2 
            IMOD=(I-1)/2
            ETMP(I)=TMP(I)*COS(IMOD*THETA)-TMP(I+1)*SIN(IMOD*THETA)
            ETMP(I+1)=TMP(I)*SIN(IMOD*THETA)+TMP(I+1)*COS(IMOD*THETA)
         ENDDO
         ETMP(NXMOD+1)=TMP(NXMOD+1)
         ETMP(NXMOD+2)=TMP(NXMOD+2)

         CALL FFTX(ETMP,WORK,TRIGSX,IFAX,NXMOD,1,NXMAX,1)

C-----END HERE

C------------------------
C     PHASE AVERAGE
C------------------------
C-----------------------------------------------------------
C Phase averaging:
C FPM(I,K) stores the average of all phase-aligned snapshots
C Each snapshot ETMP(I) is divided by NTIME and added here,
C so FPM ends up as: (ETMP₁ + ETMP₂ + ... + ETMPₙ) / NTIME
C This avoids filling FPM with a single snapshot — instead,
C we accumulate fractional contributions from each time step.
C-----------------------------------------------------------

         DO I=1,NXMOD
            FPM(I,K)=FPM(I,K)+ETMP(I)/FLOAT(NTIME)
         ENDDO
      
C-----END HERE

      ENDDO

      RETURN
      END

C=====SUBROUTINE PHASE_AVERAGE_2D END HERE
      

 




C=============================================================================
      SUBROUTINE FFTXY(F,WORK,TRIGS,IFAX,NXMOD,NYMOD,NXMAX,ISGN)
C
      REAL F(NXMAX,*),WORK(*),TRIGS(*)
      INTEGER IFAX(*)
C
      IF(ISGN.EQ.1) GO TO 50
C
      DO 10 J=NYMOD+1,NYMOD+2
      DO 10 I=1,NXMOD
 10   F(I,J)=0.
      DO 20 J=1,NYMOD+2
      DO 20 I=NXMOD+1,NXMOD+2
 20   F(I,J)=0.
C
      CALL FFTFAX(NXMOD,IFAX,TRIGS)
      CALL FFT991(F,WORK,TRIGS,IFAX,1,NXMAX,NXMOD,NYMOD,ISGN)
C
      CALL FFTFAX(NYMOD,IFAX,TRIGS)
      CALL FFT991(F,WORK,TRIGS,IFAX,NXMAX,1,NYMOD,NXMOD+2,ISGN)
      GOTO 99
C
 50   CONTINUE
C
      CALL FFTFAX(NYMOD,IFAX,TRIGS)
      CALL FFT991(F,WORK,TRIGS,IFAX,NXMAX,1,NYMOD,NXMOD+2,ISGN)
C
      CALL FFTFAX(NXMOD,IFAX,TRIGS)
      CALL FFT991(F,WORK,TRIGS,IFAX,1,NXMAX,NXMOD,NYMOD,ISGN)
C
 99   CONTINUE
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE FFTX(F,WORK,TRIGS,IFAX,NXMOD,NYMOD,NXMAX,ISGN)
C
      REAL F(NXMAX,*),WORK(*),TRIGS(*)
      INTEGER IFAX(*)
C
      IF(ISGN.EQ.1) GO TO 50
C
      DO 20 J=1,NYMOD
      DO 20 I=NXMOD+1,NXMOD+2
 20   F(I,J)=0.
C
      CALL FFTFAX(NXMOD,IFAX,TRIGS)
      CALL FFT991(F,WORK,TRIGS,IFAX,1,NXMAX,NXMOD,NYMOD,ISGN)
      GOTO 99
C
 50   CONTINUE
C
      CALL FFTFAX(NXMOD,IFAX,TRIGS)
      CALL FFT991(F,WORK,TRIGS,IFAX,1,NXMAX,NXMOD,NYMOD,ISGN)
C
 99   CONTINUE
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE FFTY(F,WORK,TRIGS,IFAX,NXMOD,NYMOD,NXMAX,ISGN)
C
      REAL F(NXMAX,*),WORK(*),TRIGS(*)
      INTEGER IFAX(*)
C
      IF(ISGN.EQ.1) GO TO 50
C
      DO 10 I=1,NXMOD
      DO 10 J=NYMOD+1,NYMOD+2
 10   F(I,J)=0.
C
      CALL FFTFAX(NYMOD,IFAX,TRIGS)
      CALL FFT991(F,WORK,TRIGS,IFAX,NXMAX,1,NYMOD,NXMOD,ISGN)
      GOTO 99
C
 50   CONTINUE
C
      CALL FFTFAX(NYMOD,IFAX,TRIGS)
      CALL FFT991(F,WORK,TRIGS,IFAX,NXMAX,1,NYMOD,NXMOD,ISGN)
C
 99   CONTINUE
      RETURN
      END
C----------------------------------------------------------------------C

C=====SUBROUTINE GROUP FFT END HERE






C===========================================================================
      SUBROUTINE PDFX(F,FX,FTMP,WORK,TRIGS,IFAX,PEX,NXMOD,NYMOD,NXMAX)
C
      REAL F(NXMAX,*),FX(NXMAX,*),FTMP(NXMAX,*),WORK(*),TRIGS(*)
      INTEGER IFAX(*)
C
      DO 10 I=1,NXMOD
      DO 10 J=1,NYMOD
       FTMP(I,J)=F(I,J)
 10   CONTINUE
      CALL FFTX(FTMP,WORK,TRIGS,IFAX,NXMOD,NYMOD,NXMAX,-1)
      DO 20 L=1,NXMOD+1,2
       PMODX=PEX*(L-1)/2
       DO 20 M=1,NYMOD
        FX(L,M)=-PMODX*FTMP(L+1,M)
        FX(L+1,M)=PMODX*FTMP(L,M)
 20   CONTINUE
      CALL FFTX(FX,WORK,TRIGS,IFAX,NXMOD,NYMOD,NXMAX,1)
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE PDFY(F,FY,FTMP,WORK,TRIGS,IFAX,PEY,NXMOD,NYMOD,NXMAX)
C
      REAL F(NXMAX,*),FY(NXMAX,*),FTMP(NXMAX,*),WORK(*),TRIGS(*)
      INTEGER IFAX(*)
C
      DO 10 I=1,NXMOD
      DO 10 J=1,NYMOD
       FTMP(I,J)=F(I,J)
 10   CONTINUE
      CALL FFTY(FTMP,WORK,TRIGS,IFAX,NXMOD,NYMOD,NXMAX,-1)
      DO 20 M=1,NYMOD+1,2
       PMODY=PEY*(M-1)/2
       DO 20 L=1,NXMOD
        FY(L,M)=-PMODY*FTMP(L,M+1)
        FY(L,M+1)=PMODY*FTMP(L,M)
 20   CONTINUE
      CALL FFTY(FY,WORK,TRIGS,IFAX,NXMOD,NYMOD,NXMAX,1)
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE PDFXY(F,FXY,FTMP,WORK,TRIGS,IFAX,PEX,PEY,
     1           NXMOD,NYMOD,NXMAX)
C
      REAL F(NXMAX,*),FXY(NXMAX,*),FTMP(NXMAX,*),WORK(*),TRIGS(*)
      INTEGER IFAX(*)
C
      DO 10 I=1,NXMOD
      DO 10 J=1,NYMOD
       FTMP(I,J)=F(I,J)
 10   CONTINUE
      CALL FFTXY(FTMP,WORK,TRIGS,IFAX,NXMOD,NYMOD,NXMAX,-1)
      DO 20 M=1,NYMOD+1,2
       MODY=(M-1)/2
       MP1=M+1
       DO 20 L=1,NXMOD+1,2
        MODX=(L-1)/2
        LP1=L+1
        AA=PEX*PEY*MODX*MODY
        FXY(L,M)=AA*FTMP(LP1,MP1)
        FXY(L,MP1)=-AA*FTMP(LP1,M)
        FXY(LP1,M)=-AA*FTMP(L,MP1)
        FXY(LP1,MP1)=AA*FTMP(L,M)
 20   CONTINUE
      CALL FFTXY(FXY,WORK,TRIGS,IFAX,NXMOD,NYMOD,NXMAX,1)
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE PDFXX(F,FXX,FTMP,WORK,TRIGS,IFAX,PEX,NXMOD,NYMOD,NXMAX)
C
      REAL F(NXMAX,*),FXX(NXMAX,*),FTMP(NXMAX,*),WORK(*),TRIGS(*)
      INTEGER IFAX(*)
C
      DO 10 I=1,NXMOD
      DO 10 J=1,NYMOD
       FTMP(I,J)=F(I,J)
 10   CONTINUE
      CALL FFTX(FTMP,WORK,TRIGS,IFAX,NXMOD,NYMOD,NXMAX,-1)
      DO 20 L=1,NXMOD+1,2
       AA=-(PEX*(L-1)/2)**2
       DO 20 M=1,NYMOD
        FXX(L,M)=AA*FTMP(L,M)
        FXX(L+1,M)=AA*FTMP(L+1,M)
 20   CONTINUE
      CALL FFTX(FXX,WORK,TRIGS,IFAX,NXMOD,NYMOD,NXMAX,1)
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE PDFYY(F,FYY,FTMP,WORK,TRIGS,IFAX,PEY,NXMOD,NYMOD,NXMAX)
C
      REAL F(NXMAX,*),FYY(NXMAX,*),FTMP(NXMAX,*),WORK(*),TRIGS(*)
      INTEGER IFAX(*)
C
      DO 10 I=1,NXMOD
      DO 10 J=1,NYMOD
       FTMP(I,J)=F(I,J)
 10   CONTINUE
      CALL FFTY(FTMP,WORK,TRIGS,IFAX,NXMOD,NYMOD,NXMAX,-1)
      DO 20 M=1,NYMOD+1,2
       AA=-(PEY*(M-1)/2)**2
       DO 20 L=1,NXMOD
        FYY(L,M)=AA*FTMP(L,M)
        FYY(L,M+1)=AA*FTMP(L,M+1)
 20   CONTINUE
      CALL FFTY(FYY,WORK,TRIGS,IFAX,NXMOD,NYMOD,NXMAX,1)
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE PDFX_(F,FTMP,WORK,TRIGS,IFAX,PEX,NXMOD,NYMOD,NXMAX)
C
      REAL F(NXMAX,*),FTMP(NXMAX,*),WORK(*),TRIGS(*)
      INTEGER IFAX(*)
C
      DO 10 I=1,NXMOD
      DO 10 J=1,NYMOD
       FTMP(I,J)=F(I,J)
 10   CONTINUE
      CALL FFTX(FTMP,WORK,TRIGS,IFAX,NXMOD,NYMOD,NXMAX,-1)
      DO 20 L=1,NXMOD+1,2
       PMODX=PEX*(L-1)/2
       DO 20 M=1,NYMOD
        F(L,M)=-PMODX*FTMP(L+1,M)
        F(L+1,M)=PMODX*FTMP(L,M)
 20   CONTINUE
      CALL FFTX(F,WORK,TRIGS,IFAX,NXMOD,NYMOD,NXMAX,1)
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE PDFY_(F,FTMP,WORK,TRIGS,IFAX,PEY,NXMOD,NYMOD,NXMAX)
C
      REAL F(NXMAX,*),FTMP(NXMAX,*),WORK(*),TRIGS(*)
      INTEGER IFAX(*)
C
      DO 10 I=1,NXMOD
      DO 10 J=1,NYMOD
       FTMP(I,J)=F(I,J)
 10   CONTINUE
      CALL FFTY(FTMP,WORK,TRIGS,IFAX,NXMOD,NYMOD,NXMAX,-1)
      DO 20 M=1,NYMOD+1,2
       PMODY=PEY*(M-1)/2
       DO 20 L=1,NXMOD
        F(L,M)=-PMODY*FTMP(L,M+1)
        F(L,M+1)=PMODY*FTMP(L,M)
 20   CONTINUE
      CALL FFTY(F,WORK,TRIGS,IFAX,NXMOD,NYMOD,NXMAX,1)
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE PDFXY_(F,FTMP,WORK,TRIGS,IFAX,PEX,PEY,
     *           NXMOD,NYMOD,NXMAX)
C
      REAL F(NXMAX,*),FTMP(NXMAX,*),WORK(*),TRIGS(*)
      INTEGER IFAX(*)
C
      DO 10 I=1,NXMOD
      DO 10 J=1,NYMOD
       FTMP(I,J)=F(I,J)
 10   CONTINUE
      CALL FFTXY(FTMP,WORK,TRIGS,IFAX,NXMOD,NYMOD,NXMAX,-1)
      DO 20 M=1,NYMOD+1,2
       MODY=(M-1)/2
       MP1=M+1
       DO 20 L=1,NXMOD+1,2
        MODX=(L-1)/2
        LP1=L+1
        AA=PEX*PEY*MODX*MODY
        F(L,M)=AA*FTMP(LP1,MP1)
        F(L,MP1)=-AA*FTMP(LP1,M)
        F(LP1,M)=-AA*FTMP(L,MP1)
        F(LP1,MP1)=AA*FTMP(L,M)
 20   CONTINUE
      CALL FFTXY(F,WORK,TRIGS,IFAX,NXMOD,NYMOD,NXMAX,1)
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE PDFXX_(F,FTMP,WORK,TRIGS,IFAX,PEX,NXMOD,NYMOD,NXMAX)
C
      REAL F(NXMAX,*),FTMP(NXMAX,*),WORK(*),TRIGS(*)
      INTEGER IFAX(*)
C
      DO 10 I=1,NXMOD
      DO 10 J=1,NYMOD
       FTMP(I,J)=F(I,J)
 10   CONTINUE
      CALL FFTX(FTMP,WORK,TRIGS,IFAX,NXMOD,NYMOD,NXMAX,-1)
      DO 20 L=1,NXMOD+1,2
       AA=-(PEX*(L-1)/2)**2
       DO 20 M=1,NYMOD
        F(L,M)=AA*FTMP(L,M)
        F(L+1,M)=AA*FTMP(L+1,M)
 20   CONTINUE
      CALL FFTX(F,WORK,TRIGS,IFAX,NXMOD,NYMOD,NXMAX,1)
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE PDFYY_(F,FTMP,WORK,TRIGS,IFAX,PEY,NXMOD,NYMOD,NXMAX)
C
      REAL F(NXMAX,*),FTMP(NXMAX,*),WORK(*),TRIGS(*)
      INTEGER IFAX(*)
C
      DO 10 I=1,NXMOD
      DO 10 J=1,NYMOD
       FTMP(I,J)=F(I,J)
 10   CONTINUE
      CALL FFTY(FTMP,WORK,TRIGS,IFAX,NXMOD,NYMOD,NXMAX,-1)
      DO 20 M=1,NYMOD+1,2
       AA=-(PEY*(M-1)/2)**2
       DO 20 L=1,NXMOD
        F(L,M)=AA*FTMP(L,M)
        F(L,M+1)=AA*FTMP(L,M+1)
 20   CONTINUE
      CALL FFTY(F,WORK,TRIGS,IFAX,NXMOD,NYMOD,NXMAX,1)
      RETURN
      END
C----------------------------------------------------------------------C

C=====SUBROUTINE GROUP DERIV END HERE






CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     FFT PACKAGE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C===========================================================================
C  Real FFTPACKTULT.F  29 January 1991
C
      SUBROUTINE FFT99(A,WORK,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)
C
C PURPOSE      PERFORMS MULTIPLE FAST FOURIER TRANSFORMS.  THIS PACKAGE
C              WILL PERFORM A NUMBER OF SIMULTANEOUS REAL/HALF-COMPLEX
C              PERIODIC FOURIER TRANSFORMS OR CORRESPONDING INVERSE
C              TRANSFORMS, I.E.  GIVEN A SET OF REAL DATA VECTORS, THE
C              PACKAGE RETURNS A SET OF 'HALF-COMPLEX' FOURIER
C              COEFFICIENT VECTORS, OR VICE VERSA.  THE LENGTH OF THE
C              TRANSFORMS MUST BE AN EVEN NUMBER GREATER THAN 4 THAT HAS
C              NO OTHER FACTORS EXCEPT POSSIBLY POWERS OF 2, 3, AND 5.
C              THIS IS AN ALL FORTRAN VERSION OF THE CRAYLIB PACKAGE
C              THAT IS MOSTLY WRITTEN IN CAL.
C
C              THE PACKAGE FFT99F CONTAINS SEVERAL USER-LEVEL ROUTINES:
C
C            SUBROUTINE FFTFAX
C                AN INITIALIZATION ROUTINE THAT MUST BE CALLED ONCE
C                BEFORE A SEQUENCE OF CALLS TO THE FFT ROUTINES
C                (PROVIDED THAT N IS NOT CHANGED).
C
C            SUBROUTINES FFT99 AND FFT991
C                TWO FFT ROUTINES THAT RETURN SLIGHTLY DIFFERENT
C                ARRANGEMENTS OF THE DATA IN GRIDPOINT SPACE.
C
C
C ACCESS       THIS FORTRAN VERSION MAY BE ACCESSED WITH
C
C                   *FORTRAN,P=XLIB,SN=FFT99F
C
C              TO ACCESS THE CRAY OBJECT CODE, CALLING THE USER ENTRY
C              POINTS FROM A CRAY PROGRAM IS SUFFICIENT.  THE SOURCE
C              FORTRAN AND CAL CODE FOR THE CRAYLIB VERSION MAY BE
C              ACCESSED USING
C
C                   FETCH P=CRAYLIB,SN=FFT99
C                   FETCH P=CRAYLIB,SN=CAL99
C
C USAGE        LET N BE OF THE FORM 2**P * 3**Q * 5**R, WHERE P .GE. 1,
C              Q .GE. 0, AND R .GE. 0.  THEN A TYPICAL SEQUENCE OF
C              CALLS TO TRANSFORM A GIVEN SET OF REAL VECTORS OF LENGTH
C              N TO A SET OF 'HALF-COMPLEX' FOURIER COEFFICIENT VECTORS
C              OF LENGTH N IS
C
C                   DIMENSION IFAX(13),TRIGS(3*N/2+1),A(M*(N+2)),
C                  +          WORK(M*(N+1))
C
C                   CALL FFTFAX (N, IFAX, TRIGS)
C                   CALL FFT99 (A,WORK,TRIGS,IFAX,INC,JUMP,N,M,ISIGN)
C
C              SEE THE INDIVIDUAL WRITE-UPS FOR FFTFAX, FFT99, AND
C              FFT991 BELOW, FOR A DETAILED DESCRIPTION OF THE
C              ARGUMENTS.
C
C HISTORY      THE PACKAGE WAS WRITTEN BY CLIVE TEMPERTON AT ECMWF IN
C              NOVEMBER, 1978.  IT WAS MODIFIED, DOCUMENTED, AND TESTED
C              FOR NCAR BY RUSS REW IN SEPTEMBER, 1980.
C
C-----------------------------------------------------------------------
C
C SUBROUTINE FFTFAX (N,IFAX,TRIGS)
C
C PURPOSE      A SET-UP ROUTINE FOR FFT99 AND FFT991.  IT NEED ONLY BE
C              CALLED ONCE BEFORE A SEQUENCE OF CALLS TO THE FFT
C              ROUTINES (PROVIDED THAT N IS NOT CHANGED).
C
C ARGUMENT     IFAX(13),TRIGS(3*N/2+1)
C DIMENSIONS
C
C ARGUMENTS
C
C ON INPUT     N
C               AN EVEN NUMBER GREATER THAN 4 THAT HAS NO PRIME FACTOR
C               GREATER THAN 5.  N IS THE LENGTH OF THE TRANSFORMS (SEE
C               THE DOCUMENTATION FOR FFT99 AND FFT991 FOR THE
C               DEFINITIONS OF THE TRANSFORMS).
C
C              IFAX
C               AN INTEGER ARRAY.  THE NUMBER OF ELEMENTS ACTUALLY USED
C               WILL DEPEND ON THE FACTORIZATION OF N.  DIMENSIONING
C               IFAX FOR 13 SUFFICES FOR ALL N LESS THAN A MILLION.
C
C              TRIGS
C               A FLOATING POINT ARRAY OF DIMENSION 3*N/2 IF N/2 IS
C               EVEN, OR 3*N/2+1 IF N/2 IS ODD.
C
C ON OUTPUT    IFAX
C               CONTAINS THE FACTORIZATION OF N/2.  IFAX(1) IS THE
C               NUMBER OF FACTORS, AND THE FACTORS THEMSELVES ARE STORED
C               IN IFAX(2),IFAX(3),...  IF FFTFAX IS CALLED WITH N ODD,
C               OR IF N HAS ANY PRIME FACTORS GREATER THAN 5, IFAX(1)
C               IS SET TO -99.
C
C              TRIGS
C               AN ARRAY OF TRIGNOMENTRIC FUNCTION VALUES SUBSEQUENTLY
C               USED BY THE FFT ROUTINES.
C
C-----------------------------------------------------------------------
C
C SUBROUTINE FFT991 (A,WORK,TRIGS,IFAX,INC,JUMP,N,M,ISIGN)
C                       AND
C SUBROUTINE FFT99 (A,WORK,TRIGS,IFAX,INC,JUMP,N,M,ISIGN)
C
C PURPOSE      PERFORM A NUMBER OF SIMULTANEOUS REAL/HALF-COMPLEX
C              PERIODIC FOURIER TRANSFORMS OR CORRESPONDING INVERSE
C              TRANSFORMS, USING ORDINARY SPATIAL ORDER OF GRIDPOINT
C              VALUES (FFT991) OR EXPLICIT CYCLIC CONTINUITY IN THE
C              GRIDPOINT VALUES (FFT99).  GIVEN A SET
C              OF REAL DATA VECTORS, THE PACKAGE RETURNS A SET OF
C              'HALF-COMPLEX' FOURIER COEFFICIENT VECTORS, OR VICE
C              VERSA.  THE LENGTH OF THE TRANSFORMS MUST BE AN EVEN
C              NUMBER THAT HAS NO OTHER FACTORS EXCEPT POSSIBLY POWERS
C              OF 2, 3, AND 5.  THESE VERSION OF FFT991 AND FFT99 ARE
C              OPTIMIZED FOR USE ON THE CRAY-1.
C
C ARGUMENT     A(M*(N+2)), WORK(M*(N+1)), TRIGS(3*N/2+1), IFAX(13)
C DIMENSIONS
C
C ARGUMENTS
C
C ON INPUT     A
C               AN ARRAY OF LENGTH M*(N+2) CONTAINING THE INPUT DATA
C               OR COEFFICIENT VECTORS.  THIS ARRAY IS OVERWRITTEN BY
C               THE RESULTS.
C
C              WORK
C               A WORK ARRAY OF DIMENSION M*(N+1)
C
C              TRIGS
C               AN ARRAY SET UP BY FFTFAX, WHICH MUST BE CALLED FIRST.
C
C              IFAX
C               AN ARRAY SET UP BY FFTFAX, WHICH MUST BE CALLED FIRST.
C
C              INC
C               THE INCREMENT (IN WORDS) BETWEEN SUCCESSIVE ELEMENTS OF
C               EACH DATA OR COEFFICIENT VECTOR (E.G.  INC=1 FOR
C               CONSECUTIVELY STORED DATA).
C
C              JUMP
C               THE INCREMENT (IN WORDS) BETWEEN THE FIRST ELEMENTS OF
C               SUCCESSIVE DATA OR COEFFICIENT VECTORS.  ON THE CRAY-1,
C               TRY TO ARRANGE DATA SO THAT JUMP IS NOT A MULTIPLE OF 8
C               (TO AVOID MEMORY BANK CONFLICTS).  FOR CLARIFICATION OF
C               INC AND JUMP, SEE THE EXAMPLES BELOW.
C
C              N
C               THE LENGTH OF EACH TRANSFORM (SEE DEFINITION OF
C               TRANSFORMS, BELOW).
C
C              M
C               THE NUMBER OF TRANSFORMS TO BE DONE SIMULTANEOUSLY.
C
C              ISIGN
C               = +1 FOR A TRANSFORM FROM FOURIER COEFFICIENTS TO
C                    GRIDPOINT VALUES.
C               = -1 FOR A TRANSFORM FROM GRIDPOINT VALUES TO FOURIER
C                    COEFFICIENTS.
C
C ON OUTPUT    A
C               IF ISIGN = +1, AND M COEFFICIENT VECTORS ARE SUPPLIED
C               EACH CONTAINING THE SEQUENCE:
C
C               A(0),B(0),A(1),B(1),...,A(N/2),B(N/2)  (N+2 VALUES)
C
C               THEN THE RESULT CONSISTS OF M DATA VECTORS EACH
C               CONTAINING THE CORRESPONDING N+2 GRIDPOINT VALUES:
C
C               FOR FFT991, X(0), X(1), X(2),...,X(N-1),0,0.
C               FOR FFT99, X(N-1),X(0),X(1),X(2),...,X(N-1),X(0).
C                   (EXPLICIT CYCLIC CONTINUITY)
C
C               WHEN ISIGN = +1, THE TRANSFORM IS DEFINED BY:
C                 X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
C                 WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
C                 AND I=SQRT (-1)
C
C               IF ISIGN = -1, AND M DATA VECTORS ARE SUPPLIED EACH
C               CONTAINING A SEQUENCE OF GRIDPOINT VALUES X(J) AS
C               DEFINED ABOVE, THEN THE RESULT CONSISTS OF M VECTORS
C               EACH CONTAINING THE CORRESPONDING FOURIER COFFICIENTS
C               A(K), B(K), 0 .LE. K .LE N/2.
C
C               WHEN ISIGN = -1, THE INVERSE TRANSFORM IS DEFINED BY:
C                 C(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*EXP(-2*I*J*K*PI/N))
C                 WHERE C(K)=A(K)+I*B(K) AND I=SQRT(-1)
C
C               A CALL WITH ISIGN=+1 FOLLOWED BY A CALL WITH ISIGN=-1
C               (OR VICE VERSA) RETURNS THE ORIGINAL DATA.
C
C               NOTE: THE FACT THAT THE GRIDPOINT VALUES X(J) ARE REAL
C               IMPLIES THAT B(0)=B(N/2)=0.  FOR A CALL WITH ISIGN=+1,
C               IT IS NOT ACTUALLY NECESSARY TO SUPPLY THESE ZEROS.
C
C EXAMPLES      GIVEN 19 DATA VECTORS EACH OF LENGTH 64 (+2 FOR EXPLICIT
C               CYCLIC CONTINUITY), COMPUTE THE CORRESPONDING VECTORS OF
C               FOURIER COEFFICIENTS.  THE DATA MAY, FOR EXAMPLE, BE
C               ARRANGED LIKE THIS:
C
C FIRST DATA   A(1)=    . . .                A(66)=             A(70)
C VECTOR       X(63) X(0) X(1) X(2) ... X(63) X(0)  (4 EMPTY LOCATIONS)
C
C SECOND DATA  A(71)=   . . .                                  A(140)
C VECTOR       X(63) X(0) X(1) X(2) ... X(63) X(0)  (4 EMPTY LOCATIONS)
C
C               AND SO ON.  HERE INC=1, JUMP=70, N=64, M=19, ISIGN=-1,
C               AND FFT99 SHOULD BE USED (BECAUSE OF THE EXPLICIT CYCLIC
C               CONTINUITY).
C
C               ALTERNATIVELY THE DATA MAY BE ARRANGED LIKE THIS:
C
C                FIRST         SECOND                          LAST
C                DATA          DATA                            DATA
C                VECTOR        VECTOR                          VECTOR
C
C                 A(1)=         A(2)=                           A(19)=
C
C                 X(63)         X(63)       . . .               X(63)
C        A(20)=   X(0)          X(0)        . . .               X(0)
C        A(39)=   X(1)          X(1)        . . .               X(1)
C                  .             .                               .
C                  .             .                               .
C                  .             .                               .
C
C               IN WHICH CASE WE HAVE INC=19, JUMP=1, AND THE REMAINING
C               PARAMETERS ARE THE SAME AS BEFORE.  IN EITHER CASE, EACH
C               COEFFICIENT VECTOR OVERWRITES THE CORRESPONDING INPUT
C               DATA VECTOR.
C
C-----------------------------------------------------------------------
      DIMENSION A(N),WORK(N),TRIGS(N),IFAX(1)
C
C     SUBROUTINE "FFT99" - MULTIPLE FAST REAL PERIODIC TRANSFORM
C     CORRESPONDING TO OLD SCALAR ROUTINE FFT9
C     PROCEDURE USED TO CONVERT TO HALF-LENGTH COMPLEX TRANSFORM
C     IS GIVEN BY COOLEY, LEWIS AND WELCH (J. SOUND VIB., VOL. 12
C     (1970), 315-337)
C
C     A IS THE ARRAY CONTAINING INPUT AND OUTPUT DATA
C     WORK IS AN AREA OF SIZE (N+1)*LOT
C     TRIGS IS A PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES
C     IFAX IS A PREVIOUSLY PREPARED LIST OF FACTORS OF N/2
C     INC IS THE INCREMENT WITHIN EACH DATA 'VECTOR'
C         (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
C     JUMP IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
C     N IS THE LENGTH OF THE DATA VECTORS
C     LOT IS THE NUMBER OF DATA VECTORS
C     ISIGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT
C           = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL
C
C     ORDERING OF COEFFICIENTS:
C         A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
C         WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED
C
C     ORDERING OF DATA:
C         X(N-1),X(0),X(1),X(2),...,X(N),X(0)
C         I.E. EXPLICIT CYCLIC CONTINUITY; (N+2) LOCATIONS REQUIRED
C
C     VECTORIZATION IS ACHIEVED ON CRAY BY DOING THE TRANSFORMS IN
C     PARALLEL
C
C     *** N.B. N IS ASSUMED TO BE AN EVEN NUMBER
C
C     DEFINITION OF TRANSFORMS:
C     -------------------------
C
C     ISIGN=+1: X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
C         WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
C
C     ISIGN=-1: A(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*COS(2*J*K*PI/N))
C               B(K)=-(1/N)*SUM(J=0,...,N-1)(X(J)*SIN(2*J*K*PI/N))
C
C
      NFAX=IFAX(1)
      NX=N+1
      NH=N/2
      INK=INC+INC
      IF (ISIGN.EQ.+1) GO TO 30
C
C     IF NECESSARY, TRANSFER DATA TO WORK AREA
      IGO=50
      IF (MOD(NFAX,2).EQ.1) GOTO 40
      IBASE=INC+1
      JBASE=1
      DO 20 L=1,LOT
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 10 M=1,N
      WORK(J)=A(I)
      I=I+INC
      J=J+1
   10 CONTINUE
      IBASE=IBASE+JUMP
      JBASE=JBASE+NX
   20 CONTINUE
C
      IGO=60
      GO TO 40
C
C     PREPROCESSING (ISIGN=+1)
C     ------------------------
C
   30 CONTINUE
      CALL FFT99A(A,WORK,TRIGS,INC,JUMP,N,LOT)
      IGO=60
C
C     COMPLEX TRANSFORM
C     -----------------
C
   40 CONTINUE
      IA=INC+1
      LA=1
      DO 80 K=1,NFAX
      IF (IGO.EQ.60) GO TO 60
   50 CONTINUE
      CALL VPASSM(A(IA),A(IA+INC),WORK(1),WORK(2),TRIGS,
     *   INK,2,JUMP,NX,LOT,NH,IFAX(K+1),LA)
      IGO=60
      GO TO 70
   60 CONTINUE
      CALL VPASSM(WORK(1),WORK(2),A(IA),A(IA+INC),TRIGS,
     *    2,INK,NX,JUMP,LOT,NH,IFAX(K+1),LA)
      IGO=50
   70 CONTINUE
      LA=LA*IFAX(K+1)
   80 CONTINUE
C
      IF (ISIGN.EQ.-1) GO TO 130
C
C     IF NECESSARY, TRANSFER DATA FROM WORK AREA
      IF (MOD(NFAX,2).EQ.1) GO TO 110
      IBASE=1
      JBASE=IA
      DO 100 L=1,LOT
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 90 M=1,N
      A(J)=WORK(I)
      I=I+1
      J=J+INC
   90 CONTINUE
      IBASE=IBASE+NX
      JBASE=JBASE+JUMP
  100 CONTINUE
C
C     FILL IN CYCLIC BOUNDARY POINTS
  110 CONTINUE
      IA=1
      IB=N*INC+1
CDIR$ IVDEP
      DO 120 L=1,LOT
      A(IA)=A(IB)
      A(IB+INC)=A(IA+INC)
      IA=IA+JUMP
      IB=IB+JUMP
  120 CONTINUE
      GO TO 140
C
C     POSTPROCESSING (ISIGN=-1):
C     --------------------------
C
  130 CONTINUE
      CALL FFT99B(WORK,A,TRIGS,INC,JUMP,N,LOT)
C
  140 CONTINUE
      RETURN
      END
C---------------------------------------------------------------C
      SUBROUTINE FFT99A(A,WORK,TRIGS,INC,JUMP,N,LOT)
      DIMENSION A(N),WORK(N),TRIGS(N)
C
C     SUBROUTINE FFT99A - PREPROCESSING STEP FOR FFT99, ISIGN=+1
C     (SPECTRAL TO GRIDPOINT TRANSFORM)
C
      NH=N/2
      NX=N+1
      INK=INC+INC
C
C     A(0) AND A(N/2)
      IA=1
      IB=N*INC+1
      JA=1
      JB=2
CDIR$ IVDEP
      DO 10 L=1,LOT
      WORK(JA)=A(IA)+A(IB)
      WORK(JB)=A(IA)-A(IB)
      IA=IA+JUMP
      IB=IB+JUMP
      JA=JA+NX
      JB=JB+NX
   10 CONTINUE
C
C     REMAINING WAVENUMBERS
      IABASE=2*INC+1
      IBBASE=(N-2)*INC+1
      JABASE=3
      JBBASE=N-1
C
      DO 30 K=3,NH,2
      IA=IABASE
      IB=IBBASE
      JA=JABASE
      JB=JBBASE
      C=TRIGS(N+K)
      S=TRIGS(N+K+1)
CDIR$ IVDEP
      DO 20 L=1,LOT
      WORK(JA)=(A(IA)+A(IB))-
     *    (S*(A(IA)-A(IB))+C*(A(IA+INC)+A(IB+INC)))
      WORK(JB)=(A(IA)+A(IB))+
     *    (S*(A(IA)-A(IB))+C*(A(IA+INC)+A(IB+INC)))
      WORK(JA+1)=(C*(A(IA)-A(IB))-S*(A(IA+INC)+A(IB+INC)))+
     *    (A(IA+INC)-A(IB+INC))
      WORK(JB+1)=(C*(A(IA)-A(IB))-S*(A(IA+INC)+A(IB+INC)))-
     *    (A(IA+INC)-A(IB+INC))
      IA=IA+JUMP
      IB=IB+JUMP
      JA=JA+NX
      JB=JB+NX
   20 CONTINUE
      IABASE=IABASE+INK
      IBBASE=IBBASE-INK
      JABASE=JABASE+2
      JBBASE=JBBASE-2
   30 CONTINUE
C
      IF (IABASE.NE.IBBASE) GO TO 50
C     WAVENUMBER N/4 (IF IT EXISTS)
      IA=IABASE
      JA=JABASE
CDIR$ IVDEP
      DO 40 L=1,LOT
      WORK(JA)=2.0*A(IA)
      WORK(JA+1)=-2.0*A(IA+INC)
      IA=IA+JUMP
      JA=JA+NX
   40 CONTINUE
C
   50 CONTINUE
      RETURN
      END
C---------------------------------------------------------------C
      SUBROUTINE FFT99B(WORK,A,TRIGS,INC,JUMP,N,LOT)
      DIMENSION WORK(N),A(N),TRIGS(N)
C
C     SUBROUTINE FFT99B - POSTPROCESSING STEP FOR FFT99, ISIGN=-1
C     (GRIDPOINT TO SPECTRAL TRANSFORM)
C
      NH=N/2
      NX=N+1
      INK=INC+INC
C
C     A(0) AND A(N/2)
      SCALE=1.0/FLOAT(N)
      IA=1
      IB=2
      JA=1
      JB=N*INC+1
CDIR$ IVDEP
      DO 10 L=1,LOT
      A(JA)=SCALE*(WORK(IA)+WORK(IB))
      A(JB)=SCALE*(WORK(IA)-WORK(IB))
      A(JA+INC)=0.0
      A(JB+INC)=0.0
      IA=IA+NX
      IB=IB+NX
      JA=JA+JUMP
      JB=JB+JUMP
   10 CONTINUE
C
C     REMAINING WAVENUMBERS
      SCALE=0.5*SCALE
      IABASE=3
      IBBASE=N-1
      JABASE=2*INC+1
      JBBASE=(N-2)*INC+1
C
      DO 30 K=3,NH,2
      IA=IABASE
      IB=IBBASE
      JA=JABASE
      JB=JBBASE
      C=TRIGS(N+K)
      S=TRIGS(N+K+1)
CDIR$ IVDEP
      DO 20 L=1,LOT
      A(JA)=SCALE*((WORK(IA)+WORK(IB))
     *   +(C*(WORK(IA+1)+WORK(IB+1))+S*(WORK(IA)-WORK(IB))))
      A(JB)=SCALE*((WORK(IA)+WORK(IB))
     *   -(C*(WORK(IA+1)+WORK(IB+1))+S*(WORK(IA)-WORK(IB))))
      A(JA+INC)=SCALE*((C*(WORK(IA)-WORK(IB))-S*(WORK(IA+1)+WORK(IB+1)))
     *    +(WORK(IB+1)-WORK(IA+1)))
      A(JB+INC)=SCALE*((C*(WORK(IA)-WORK(IB))-S*(WORK(IA+1)+WORK(IB+1)))
     *    -(WORK(IB+1)-WORK(IA+1)))
      IA=IA+NX
      IB=IB+NX
      JA=JA+JUMP
      JB=JB+JUMP
   20 CONTINUE
      IABASE=IABASE+2
      IBBASE=IBBASE-2
      JABASE=JABASE+INK
      JBBASE=JBBASE-INK
   30 CONTINUE
C
      IF (IABASE.NE.IBBASE) GO TO 50
C     WAVENUMBER N/4 (IF IT EXISTS)
      IA=IABASE
      JA=JABASE
      SCALE=2.0*SCALE
CDIR$ IVDEP
      DO 40 L=1,LOT
      A(JA)=SCALE*WORK(IA)
      A(JA+INC)=-SCALE*WORK(IA+1)
      IA=IA+NX
      JA=JA+JUMP
   40 CONTINUE
C
   50 CONTINUE
      RETURN
      END
C---------------------------------------------------------------C
      SUBROUTINE FFT991(A,WORK,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)
      DIMENSION A(N),WORK(N),TRIGS(N),IFAX(1)
C
C     SUBROUTINE "FFT991" - MULTIPLE REAL/HALF-COMPLEX PERIODIC
C     FAST FOURIER TRANSFORM
C
C     SAME AS FFT99 EXCEPT THAT ORDERING OF DATA CORRESPONDS TO
C     THAT IN MRFFT2
C
C     PROCEDURE USED TO CONVERT TO HALF-LENGTH COMPLEX TRANSFORM
C     IS GIVEN BY COOLEY, LEWIS AND WELCH (J. SOUND VIB., VOL. 12
C     (1970), 315-337)
C
C     A IS THE ARRAY CONTAINING INPUT AND OUTPUT DATA
C     WORK IS AN AREA OF SIZE (N+1)*LOT
C     TRIGS IS A PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES
C     IFAX IS A PREVIOUSLY PREPARED LIST OF FACTORS OF N/2
C     INC IS THE INCREMENT WITHIN EACH DATA 'VECTOR'
C         (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
C     JUMP IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
C     N IS THE LENGTH OF THE DATA VECTORS
C     LOT IS THE NUMBER OF DATA VECTORS
C     ISIGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT
C           = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL
C
C     ORDERING OF COEFFICIENTS:
C         A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
C         WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED
C
C     ORDERING OF DATA:
C         X(0),X(1),X(2),...,X(N-1)
C
C     VECTORIZATION IS ACHIEVED ON CRAY BY DOING THE TRANSFORMS IN
C     PARALLEL
C
C     *** N.B. N IS ASSUMED TO BE AN EVEN NUMBER
C
C     DEFINITION OF TRANSFORMS:
C     -------------------------
C
C     ISIGN=+1: X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
C         WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
C
C     ISIGN=-1: A(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*COS(2*J*K*PI/N))
C               B(K)=-(1/N)*SUM(J=0,...,N-1)(X(J)*SIN(2*J*K*PI/N))
C
      NFAX=IFAX(1)
      NX=N+1
      NH=N/2
      INK=INC+INC
      IF (ISIGN.EQ.+1) GO TO 30
C
C     IF NECESSARY, TRANSFER DATA TO WORK AREA
      IGO=50
      IF (MOD(NFAX,2).EQ.1) GOTO 40
      IBASE=1
      JBASE=1
      DO 20 L=1,LOT
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 10 M=1,N
      WORK(J)=A(I)
      I=I+INC
      J=J+1
   10 CONTINUE
      IBASE=IBASE+JUMP
      JBASE=JBASE+NX
   20 CONTINUE
C
      IGO=60
      GO TO 40
C
C     PREPROCESSING (ISIGN=+1)
C     ------------------------
C
   30 CONTINUE
      CALL FFT99A(A,WORK,TRIGS,INC,JUMP,N,LOT)
      IGO=60
C
C     COMPLEX TRANSFORM
C     -----------------
C
   40 CONTINUE
      IA=1
      LA=1
      DO 80 K=1,NFAX
      IF (IGO.EQ.60) GO TO 60
   50 CONTINUE
      CALL VPASSM(A(IA),A(IA+INC),WORK(1),WORK(2),TRIGS,
     *   INK,2,JUMP,NX,LOT,NH,IFAX(K+1),LA)
      IGO=60
      GO TO 70
   60 CONTINUE
      CALL VPASSM(WORK(1),WORK(2),A(IA),A(IA+INC),TRIGS,
     *    2,INK,NX,JUMP,LOT,NH,IFAX(K+1),LA)
      IGO=50
   70 CONTINUE
      LA=LA*IFAX(K+1)
   80 CONTINUE
C
      IF (ISIGN.EQ.-1) GO TO 130
C
C     IF NECESSARY, TRANSFER DATA FROM WORK AREA
      IF (MOD(NFAX,2).EQ.1) GO TO 110
      IBASE=1
      JBASE=1
      DO 100 L=1,LOT
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 90 M=1,N
      A(J)=WORK(I)
      I=I+1
      J=J+INC
   90 CONTINUE
      IBASE=IBASE+NX
      JBASE=JBASE+JUMP
  100 CONTINUE
C
C     FILL IN ZEROS AT END
  110 CONTINUE
      IB=N*INC+1
CDIR$ IVDEP
      DO 120 L=1,LOT
      A(IB)=0.0
      A(IB+INC)=0.0
      IB=IB+JUMP
  120 CONTINUE
      GO TO 140
C
C     POSTPROCESSING (ISIGN=-1):
C     --------------------------
C
  130 CONTINUE
      CALL FFT99B(WORK,A,TRIGS,INC,JUMP,N,LOT)
C
  140 CONTINUE
      RETURN
      END
C---------------------------------------------------------------C
      SUBROUTINE FFTFAX(N,IFAX,TRIGS)
      DIMENSION IFAX(13),TRIGS(1)
C
C MODE 3 IS USED FOR REAL/HALF-COMPLEX TRANSFORMS.  IT IS POSSIBLE
C TO DO COMPLEX/COMPLEX TRANSFORMS WITH OTHER VALUES OF MODE, BUT
C DOCUMENTATION OF THE DETAILS WERE NOT AVAILABLE WHEN THIS ROUTINE
C WAS WRITTEN.
C
      DATA MODE /3/
      CALL FAX (IFAX, N, MODE)
      I = IFAX(1)
      IF (IFAX(I+1) .GT. 5 .OR. N .LE. 4) IFAX(1) = -99
      IF (IFAX(1) .LE. 0 )THEN
        WRITE(*,1900)N
1900    FORMAT(' FFTFAX - Invalid N=',I20)
        RETURN
        ENDIF
      CALL FFTRIG (TRIGS, N, MODE)
      RETURN
      END
C---------------------------------------------------------------C
      SUBROUTINE FAX(IFAX,N,MODE)
      DIMENSION IFAX(10)
      NN=N
      IF (IABS(MODE).EQ.1) GO TO 10
      IF (IABS(MODE).EQ.8) GO TO 10
      NN=N/2
      IF ((NN+NN).EQ.N) GO TO 10
      IFAX(1)=-99
      RETURN
   10 K=1
C     TEST FOR FACTORS OF 4
   20 IF (MOD(NN,4).NE.0) GO TO 30
      K=K+1
      IFAX(K)=4
      NN=NN/4
      IF (NN.EQ.1) GO TO 80
      GO TO 20
C     TEST FOR EXTRA FACTOR OF 2
   30 IF (MOD(NN,2).NE.0) GO TO 40
      K=K+1
      IFAX(K)=2
      NN=NN/2
      IF (NN.EQ.1) GO TO 80
C     TEST FOR FACTORS OF 3
   40 IF (MOD(NN,3).NE.0) GO TO 50
      K=K+1
      IFAX(K)=3
      NN=NN/3
      IF (NN.EQ.1) GO TO 80
      GO TO 40
C     NOW FIND REMAINING FACTORS
   50 L=5
      INC=2
C     INC ALTERNATELY TAKES ON VALUES 2 AND 4
   60 IF (MOD(NN,L).NE.0) GO TO 70
      K=K+1
      IFAX(K)=L
      NN=NN/L
      IF (NN.EQ.1) GO TO 80
      GO TO 60
   70 L=L+INC
      INC=6-INC
      GO TO 60
   80 IFAX(1)=K-1
C     IFAX(1) CONTAINS NUMBER OF FACTORS
      NFAX=IFAX(1)
C     SORT FACTORS INTO ASCENDING ORDER
      IF (NFAX.EQ.1) GO TO 110
      DO 100 II=2,NFAX
      ISTOP=NFAX+2-II
      DO 90 I=2,ISTOP
      IF (IFAX(I+1).GE.IFAX(I)) GO TO 90
      ITEM=IFAX(I)
      IFAX(I)=IFAX(I+1)
      IFAX(I+1)=ITEM
   90 CONTINUE
  100 CONTINUE
  110 CONTINUE
      RETURN
      END
C---------------------------------------------------------------C
      SUBROUTINE FFTRIG(TRIGS,N,MODE)
c Sudin_edit
c      DIMENSION TRIGS(1)
c Sudin_edit
      DIMENSION TRIGS(*)
      PI=2.0*ASIN(1.0)
      IMODE=IABS(MODE)
      NN=N
      IF (IMODE.GT.1.AND.IMODE.LT.6) NN=N/2
      DEL=(PI+PI)/FLOAT(NN)
      L=NN+NN
      DO 10 I=1,L,2
         ANGLE=0.5*FLOAT(I-1)*DEL
         TRIGS(I)=COS(ANGLE)
         TRIGS(I+1)=SIN(ANGLE)
   10 CONTINUE
      IF (IMODE.EQ.1) RETURN
      IF (IMODE.EQ.8) RETURN
      DEL=0.5*DEL
      NH=(NN+1)/2
      L=NH+NH
      LA=NN+NN
      DO 20 I=1,L,2
         ANGLE=0.5*FLOAT(I-1)*DEL
         TRIGS(LA+I)=COS(ANGLE)
         TRIGS(LA+I+1)=SIN(ANGLE)
   20 CONTINUE
      IF (IMODE.LE.3) RETURN
      DEL=0.5*DEL
      LA=LA+NN
      IF (MODE.EQ.5) GO TO 40
      DO 30 I=2,NN
         ANGLE=FLOAT(I-1)*DEL
         TRIGS(LA+I)=2.0*SIN(ANGLE)
   30 CONTINUE
      RETURN
   40 CONTINUE
      DEL=0.5*DEL
      DO 50 I=2,N
         ANGLE=FLOAT(I-1)*DEL
         TRIGS(LA+I)=SIN(ANGLE)
   50 CONTINUE
      RETURN
      END
      SUBROUTINE VPASSM(A,B,C,D,TRIGS,INC1,INC2,INC3,INC4,LOT,N,IFAC,LA)
      DIMENSION A(N),B(N),C(N),D(N),TRIGS(N)
C
C     SUBROUTINE "VPASSM" - MULTIPLE VERSION OF "VPASSA"
C     PERFORMS ONE PASS THROUGH DATA
C     AS PART OF MULTIPLE COMPLEX FFT ROUTINE
C     A IS FIRST REAL INPUT VECTOR
C     B IS FIRST IMAGINARY INPUT VECTOR
C     C IS FIRST REAL OUTPUT VECTOR
C     D IS FIRST IMAGINARY OUTPUT VECTOR
C     TRIGS IS PRECALCULATED TABLE OF SINES " COSINES
C     INC1 IS ADDRESSING INCREMENT FOR A AND B
C     INC2 IS ADDRESSING INCREMENT FOR C AND D
C     INC3 IS ADDRESSING INCREMENT BETWEEN A"S & B"S
C     INC4 IS ADDRESSING INCREMENT BETWEEN C"S & D"S
C     LOT IS THE NUMBER OF VECTORS
C     N IS LENGTH OF VECTORS
C     IFAC IS CURRENT FACTOR OF N
C     LA IS PRODUCT OF PREVIOUS FACTORS
C
      DATA SIN36/0.587785252292473/,COS36/0.809016994374947/,
     *     SIN72/0.951056516295154/,COS72/0.309016994374947/,
     *     SIN60/0.866025403784437/
C
      M=N/IFAC
      IINK=M*INC1
      JINK=LA*INC2
      JUMP=(IFAC-1)*JINK
      IBASE=0
      JBASE=0
      IGO=IFAC-1
      IF (IGO.GT.4) RETURN
      GO TO (10,50,90,130),IGO
C
C     CODING FOR FACTOR 2
C
   10 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      DO 20 L=1,LA
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 15 IJK=1,LOT
      C(JA+J)=A(IA+I)+A(IB+I)
      D(JA+J)=B(IA+I)+B(IB+I)
      C(JB+J)=A(IA+I)-A(IB+I)
      D(JB+J)=B(IA+I)-B(IB+I)
      I=I+INC3
      J=J+INC4
   15 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   20 CONTINUE
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 40 K=LA1,M,LA
      KB=K+K-2
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      DO 30 L=1,LA
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 25 IJK=1,LOT
      C(JA+J)=A(IA+I)+A(IB+I)
      D(JA+J)=B(IA+I)+B(IB+I)
      C(JB+J)=C1*(A(IA+I)-A(IB+I))-S1*(B(IA+I)-B(IB+I))
      D(JB+J)=S1*(A(IA+I)-A(IB+I))+C1*(B(IA+I)-B(IB+I))
      I=I+INC3
      J=J+INC4
   25 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   30 CONTINUE
      JBASE=JBASE+JUMP
   40 CONTINUE
      RETURN
C
C     CODING FOR FACTOR 3
C
   50 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      IC=IB+IINK
      JC=JB+JINK
      DO 60 L=1,LA
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 55 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IC+I))
      C(JB+J)=(A(IA+I)-0.5*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)-B(IC+I)))
      C(JC+J)=(A(IA+I)-0.5*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)-B(IC+I)))
      D(JB+J)=(B(IA+I)-0.5*(B(IB+I)+B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I)))
      D(JC+J)=(B(IA+I)-0.5*(B(IB+I)+B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I)))
      I=I+INC3
      J=J+INC4
   55 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   60 CONTINUE
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 80 K=LA1,M,LA
      KB=K+K-2
      KC=KB+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      DO 70 L=1,LA
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 65 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IC+I))
      C(JB+J)=
     *    C1*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)-B(IC+I))))
     *   -S1*((B(IA+I)-0.5*(B(IB+I)+B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I))))
      D(JB+J)=
     *    S1*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)-B(IC+I))))
     *   +C1*((B(IA+I)-0.5*(B(IB+I)+B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I))))
      C(JC+J)=
     *    C2*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)-B(IC+I))))
     *   -S2*((B(IA+I)-0.5*(B(IB+I)+B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I))))
      D(JC+J)=
     *    S2*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)-B(IC+I))))
     *   +C2*((B(IA+I)-0.5*(B(IB+I)+B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I))))
      I=I+INC3
      J=J+INC4
   65 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   70 CONTINUE
      JBASE=JBASE+JUMP
   80 CONTINUE
      RETURN
C
C     CODING FOR FACTOR 4
C
   90 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      IC=IB+IINK
      JC=JB+JINK
      ID=IC+IINK
      JD=JC+JINK
      DO 100 L=1,LA
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 95 IJK=1,LOT
      C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
      C(JC+J)=(A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I))
      D(JA+J)=(B(IA+I)+B(IC+I))+(B(IB+I)+B(ID+I))
      D(JC+J)=(B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I))
      C(JB+J)=(A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I))
      C(JD+J)=(A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I))
      D(JB+J)=(B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I))
      D(JD+J)=(B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I))
      I=I+INC3
      J=J+INC4
   95 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  100 CONTINUE
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 120 K=LA1,M,LA
      KB=K+K-2
      KC=KB+KB
      KD=KC+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      C3=TRIGS(KD+1)
      S3=TRIGS(KD+2)
      DO 110 L=1,LA
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 105 IJK=1,LOT
      C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
      D(JA+J)=(B(IA+I)+B(IC+I))+(B(IB+I)+B(ID+I))
      C(JC+J)=
     *    C2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))
     *   -S2*((B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I)))
      D(JC+J)=
     *    S2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))
     *   +C2*((B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I)))
      C(JB+J)=
     *    C1*((A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I)))
     *   -S1*((B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I)))
      D(JB+J)=
     *    S1*((A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I)))
     *   +C1*((B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I)))
      C(JD+J)=
     *    C3*((A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I)))
     *   -S3*((B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I)))
      D(JD+J)=
     *    S3*((A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I)))
     *   +C3*((B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I)))
      I=I+INC3
      J=J+INC4
  105 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  110 CONTINUE
      JBASE=JBASE+JUMP
  120 CONTINUE
      RETURN
C
C     CODING FOR FACTOR 5
C
  130 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      IC=IB+IINK
      JC=JB+JINK
      ID=IC+IINK
      JD=JC+JINK
      IE=ID+IINK
      JE=JD+JINK
      DO 140 L=1,LA
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 135 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IE+I))+(B(IC+I)+B(ID+I))
      C(JB+J)=(A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     *  -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I)))
      C(JE+J)=(A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     *  +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I)))
      D(JB+J)=(B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     *  +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I)))
      D(JE+J)=(B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     *  -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I)))
      C(JC+J)=(A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     *  -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I)))
      C(JD+J)=(A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     *  +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I)))
      D(JC+J)=(B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     *  +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I)))
      D(JD+J)=(B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     *  -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I)))
      I=I+INC3
      J=J+INC4
  135 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  140 CONTINUE
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 160 K=LA1,M,LA
      KB=K+K-2
      KC=KB+KB
      KD=KC+KB
      KE=KD+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      C3=TRIGS(KD+1)
      S3=TRIGS(KD+2)
      C4=TRIGS(KE+1)
      S4=TRIGS(KE+2)
      DO 150 L=1,LA
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 145 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IE+I))+(B(IC+I)+B(ID+I))
      C(JB+J)=
     *    C1*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     *      -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))
     *   -S1*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     *      +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      D(JB+J)=
     *    S1*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     *      -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))
     *   +C1*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     *      +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      C(JE+J)=
     *    C4*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     *      +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))
     *   -S4*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     *      -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      D(JE+J)=
     *    S4*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     *      +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))
     *   +C4*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     *      -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      C(JC+J)=
     *    C2*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     *      -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))
     *   -S2*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     *      +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      D(JC+J)=
     *    S2*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     *      -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))
     *   +C2*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     *      +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      C(JD+J)=
     *    C3*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     *      +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))
     *   -S3*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     *      -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      D(JD+J)=
     *    S3*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     *      +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))
     *   +C3*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     *      -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      I=I+INC3
      J=J+INC4
  145 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  150 CONTINUE
      JBASE=JBASE+JUMP
  160 CONTINUE
      RETURN
      END
C---------------------------------------------------------------C
      SUBROUTINE FACT(N,IFAX)
C     FACTORIZATION ROUTINE THAT FIRST EXTRACTS ALL FACTORS OF 4
      DIMENSION IFAX(13)
      IF (N.GT.1) GO TO 10
      IFAX(1) = 0
      IF (N.LT.1) IFAX(1) = -99
      RETURN
   10 NN=N
      K=1
C     TEST FOR FACTORS OF 4
   20 IF (MOD(NN,4).NE.0) GO TO 30
      K=K+1
      IFAX(K)=4
      NN=NN/4
      IF (NN.EQ.1) GO TO 80
      GO TO 20
C     TEST FOR EXTRA FACTOR OF 2
   30 IF (MOD(NN,2).NE.0) GO TO 40
      K=K+1
      IFAX(K)=2
      NN=NN/2
      IF (NN.EQ.1) GO TO 80
C     TEST FOR FACTORS OF 3
   40 IF (MOD(NN,3).NE.0) GO TO 50
      K=K+1
      IFAX(K)=3
      NN=NN/3
      IF (NN.EQ.1) GO TO 80
      GO TO 40
C     NOW FIND REMAINING FACTORS
   50 L=5
      MAX = SQRT(FLOAT(NN))
      INC=2
C     INC ALTERNATELY TAKES ON VALUES 2 AND 4
   60 IF (MOD(NN,L).NE.0) GO TO 70
      K=K+1
      IFAX(K)=L
      NN=NN/L
      IF (NN.EQ.1) GO TO 80
      GO TO 60
   70 IF (L.GT.MAX) GO TO 75
      L=L+INC
      INC=6-INC
      GO TO 60
   75 K = K+1
      IFAX(K) = NN
   80 IFAX(1)=K-1
C     IFAX(1) NOW CONTAINS NUMBER OF FACTORS
      RETURN
      END

C=====SUBROUTINES FOR FFT PACKAGE END HERE
