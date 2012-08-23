

int RipleyK()
{
	int i,j;
	double u=0;		// Distancia entre puntos i & j
	double t=0;		// Distancia a la cual se evalua la funcion
		
	for(i=0; i<noPoints; i++)
		for(j=0; j<noPoints; j++)
		{
			u= sqrt( pow(x(i)-x(j),2) + pow(y(i)-y(j),2) );
			
			
			
	
x(j)
y(j)
}


C----------
C  RIPK.FOR
C----------

//----------
//  Melinda Moeur 
//  Intermountain Research Station
//  1221 S. Main, Moscow ID  83843
//  (208) 883-2316
//  E-mail: moeur@forest.moscowfsl.wsu.edu 

//  Original version 4/12/90

//  Ripley's K (or second-order neighborhood) analysis program for computing 
//  proportion of total possible pairs of points whose pair members are within 
//  a specified distance of one another.  

//  For use with stem map data on rectangular plot, includes edge correction. 
//  Works for both univariate analysis (WITHIN-CLASS ANALYSIS, in which FROM 
//  and TO points are the same), and bivariate analysis (BETWEEN-CLASS 
//  ANALYSIS, FROM and TO points different.  

//  If the analysis is univariate (testing the spatial pattern of a single 
//  population): The Monte Carlo routine (UNIMC) draws from a Poisson 
//  distributions of points, and the observed K-distribution is compared 
//  to the confidence envelope to determine departures from randomness.  

//  For bivariate patterns (testing the interaction between points in 
//  population 1 vs. population 2): The Monte Carlo routine (BIMC) shifts
//  the 'TO' points a random amount in the X and Y directions, while holding 
//  them in their original relative locations.  Then the K-distribution is 
//  computed between the original 'FROM' points and the shifted 'TO' points 
//  to determine the confidence envelope boundaries.
//----------
      INCLUDE 'RIPKCOM.FOR'

      char DATE[11],TIME[9]

      char NAMSAV[81],FMTSAV[81]

      double S0,S1,SEED

// 10 CONTINUE

//  SEED THE RANDOM GENERATOR
      S0=117193557D0 


      BIFLAG = 1
//----------
// READ TITLE AND FILE NAMES FROM RESPONSE FILE
//----------
//
//	NAME2 : Control file
//	NAME3 : From File
//	FMT3: Formato de NAME3
//	NAME4 : TO file
//	FMT4: Formato de NAME4
//  NAME7 : Output
//  NAME8 : Output
      
//	MCFLAG : Monte Carlo (0=NO, 1=YES)
//  NTIMES : Number of simulations
//	ALPHAMC: Alpha level
      
//  COMPUTE INDEX OF LOWER % AND UPPER% SIGNIFICANCE VALUES OF MC ENVELOPE
	IALPHA = IFIX(100.-(ALPHAMC*200.)+.5)
	ILOW = IFIX(ALPHAMC*NTIMES+.5)
	IUP = NTIMES - ILOW + 1

//----------
//  READ PLOT DIMENSIONS AND DISTANCE LIMITS FROM CONTROL FILE.
//----------
      Cntrl()

//----------
//  READ X,Y DATA IN "FROM" AND "TO" GROUPS.
//----------
//  20 CONTINUE

      Datain()

//----------
//  CALL THE IO ROUTINE TO ECHO INPUT.
//----------
      Iout()

//----------
// CALCULATE RIPLEY'S K FOR OBSERVED DATA.                            
//----------
      Ripk()

//----------
//  BEGIN MONTE CARLO SEQUENCE.
//----------
//----------
//  UNIVARIATE CASE: RANDOM DRAWS
//----------
	if(NAME3==NAME4)
	{
		for(NITER=0; NITER<NTIMES; NITER++)
      		UNIMC(NITER)
	}
	else
//----------
//  BIVARIATE CASE: TOROID SHIFT
//----------
	{
		for(NITER=0; NITER<NTIMES; NITER++)
      		BIMC(NITER)
	}

//----------
//  IF THIS IS A BIVARIATE ANALYSIS, SET THE LOGICAL FLAG.
//  SWITCH THE 'FROM' AND 'TO' POINTS AND REPEAT THE CALCULATIONS 
//  FOR THE ORIGINAL POINTS.
//----------
      if(BIFLAG==)
      {
	      if(NAME3!=NAME4)
	      {
		      BIFLAG = 2
	    	  NAMSAV = NAME3
		      NAME3 = NAME4
		      NAME4 = NAMSAV
		      FMTSAV = FMT3
		      FMT3 = FMT4
		      FMT4 = FMTSAV
		      GO TO 20
	      }
      }



C----------
C CALCULATE LINEARIZED STATISTICS AND WRITE OUTPUT.
C----------
      CALL RKOUT

  110 CONTINUE

C DATE AND TIME STAMP
      CALL GRDTIM (DATE,TIME)
      WRITE (7,6020) DATE,TIME

C  CHECK TO SEE IF THERE IS ANOTHER RUN TO PROCESS.
      GO TO 10

 1000 CONTINUE

      WRITE(*,6030)
 6030 FORMAT(' **DONE**')
      
      STOP
      END
C**********************************************************************

      SUBROUTINE CNTRL
C----------
C  READ VALUES FROM CONTROL FILE:
C  UNITS, MAXIMUM DISTANCE, DISTANCE INCREMENT, AND X,Y COORDINATES 
C  OF PLOT CORNERS (in order, SW, NW, NE, SE).
C----------
      INCLUDE 'RIPKCOM.FOR'

      OPEN(2,FILE=NAME2,STATUS='UNKNOWN')
      READ (2,2000) UNITS,DISTMAX,DINC,(CX(I),CY(I),I=1,4) 
 2000 FORMAT (A6/2F10.0,4(/2F10.0))

      NINC=1 + DISTMAX/DINC + .5
      DMAX2 = DISTMAX*DISTMAX

C COMPUTE BOUNDARY LENGTHS FOR THE TOROID ADJUSTMENT FOR BIVARIATE 
C  MONTE CARLO ANALYSIS.  ALSO COMPUTE PLOT AREA.
      SOUTH = ((CX(4)-CX(1)) + (CX(3)-CX(2)))/2.0
      WEST =  ((CY(2)-CY(1)) + (CY(3)-CY(4)))/2.0
      AREA = SOUTH*WEST

      CLOSE (2)
      RETURN
      END

C**********************************************************************

      SUBROUTINE DATAIN
C----------
C  READ X,Y COORDINATES OF POINTS IN "FROM" AND "TO" GROUPS, AND COUNT 
C  THE NUMBER OF POINTS IN EACH GROUP.                
C----------
      INCLUDE 'RIPKCOM.FOR'

      OPEN(3,FILE=NAME3,STATUS='UNKNOWN')
      NFROM=0
      DO 10 I=1,MAXFROM
      READ (3,*,END=15) XFROM(I),YFROM(I)
      NFROM = NFROM + 1
   10 CONTINUE
   15 CONTINUE
      CLOSE (3)  

      OPEN(4,FILE=NAME4,STATUS='UNKNOWN')
      NTO=0
      DO 20 I=1,MAXTO
      READ (4,*,END=25) XTO(I),YTO(I)
      NTO = NTO + 1
   20 CONTINUE
   25 CONTINUE
      CLOSE (4)  

      IF (BIFLAG .EQ. 1) THEN
      N1 = NFROM
      N2 = NTO
      ENDIF

      RETURN
      END

C*********************************************************************

      SUBROUTINE IOUT
C----------
C  WRITE OUTPUT.
C----------
      INCLUDE 'RIPKCOM.FOR'

      DOUBLE PRECISION S0,S1
      COMMON /RANCOM/ S0,S1

      IF (BIFLAG .EQ. 2) GO TO 10

      WRITE (7,6010) TITLE
 6010 FORMAT (' *************** Ripley''s K (second-order',
     &    ' neighborhood) Analysis ***************' // 
     &    1X,79('-')/ 1X,A80)

      IF (NAME3 .EQ. NAME4) THEN
      WRITE (7,6049) 
 6049 FORMAT 
     &(' ***** Univariate (within-class) Analysis ******' / 1X,79('-'))
      ELSE
      WRITE (7,6050) 
 6050 FORMAT 
     &(' ***** Bivariate (between-class) Analysis ******' / 1X,79('-'))
      ENDIF

      WRITE(7,6020) NAME2,NAME3,FMT3,NAME4,FMT4,NAME7,NAME8
 6020 FORMAT(' INPUT FILES:'/
     & ' UNIT2, CONTROL   :',  T22,A80/
     & ' UNIT3, ''FROM'' PTS:',T22,A80/
     & ' FMT3,            :',  T22,A80/
     & ' UNIT4, ''TO'' PTS  :',T22,A80/
     & ' FMT4,            :',  T22,A80/
     & ' UNIT7, REPORT    :',  T22,A80/
     & ' UNIT8, OUTPUT    :',  T22,A80)

      IF (MCFLAG .EQ. 1) THEN
      WRITE (7,6028) 
     &     MCFLAG,NTIMES,ALPHAMC,IALPHA,ILOW,IUP
 6028 FORMAT (/ 
     & ' Monte Carlo flag (0=NO, 1=YES) = ',I2 /
     & ' Number of iterations = ',I5,
     & '     One-sided Alpha Level = ',F6.3 /
     & ' Indices of ',I3,'% significance level for Monte ',
     &        'Carlo confidence bounds',2i5 )

      ELSE
      WRITE (7,6029) MCFLAG
 6029 FORMAT (/ 
     & ' Monte Carlo flag (0=NO, 1=YES) = ',I2)
      ENDIF
      WRITE (7,6030) ISEED,S0
 6030 FORMAT (/ 
     & ' Re-seed Random Number Generator (0=NO, 1=YES) = ',I2 /
     & ' Random number seed = ',f20.1)

      WRITE (7,6031) AREA,UNITS,DISTMAX,UNITS,DINC,UNITS
 6031 FORMAT (/ 
     & ' Plot Area: ',F10.3,'  Sq ',A6 //
     & ' Max. d for F(d) = ',F5.1,'  ',A6,'       d Increment = ', 
     & F5.2,'  ',A6)  

      WRITE (7,6060) (CX(I),CY(I),I=1,4),SOUTH,WEST  
 6060 FORMAT (/' X,Y Coordinates of Plot Corners' /
     &        '           X         Y '/
     &        ' SW:',2F10.3/' NW:',2F10.3/
     &        ' NE:',2F10.3/' SE:',2F10.3//
     &        ' Average Boundary Lengths (South--X, West--Y):',
     &          2f10.3)    

   10 CONTINUE
      IF (BIFLAG .EQ. 2) WRITE (7,6065)
 6065 FORMAT (/' ********** Switching ''FROM'' and ''TO'' Points for ',
     &        'Bivariate Analysis **********')
      WRITE (7,6040) NFROM,NTO
 6040 FORMAT (/' Number of Points in the ''FROM'' and ''TO'' Datasets:'/
     &       18X,2I10)

      WRITE (7,6070) XFROM(1),YFROM(1),XFROM(NFROM),YFROM(NFROM),
     &               XTO(1),YTO(1),XTO(NTO),YTO(NTO)
 6070 FORMAT (/'      X         Y   Coordinates of First and Last ',
     & 'Point in ''FROM'' Dataset',2(/2F10.2)/
     &        /'      X         Y   Coordinates of First and Last ',
     & 'Point in ''TO'' Dataset',2(/2F10.2) )

      WRITE (7,6080)
 6080 FORMAT (/' Distance Matrix for First and Last Points'/ 
     & ' NFROM  -------------------Distance to 20 Closest ''TO'' Points',
     & '-----------------')  

      RETURN
      END

C***********************************************************************

      SUBROUTINE RIPK
C----------
C  COMPUTE RIPLEY'S K, SECOND ORDER SPATIAL STATISTIC FOR MAPPED POINT 
C  PATTERNS, FOR THE OBSERVED POINT LOCATIONS.  
C
C  THE CUMULATIVE DISTRIBUTION (RIPLEY'S K) IS CORRECTED FOR EDGE
C  EFFECTS (THE PROPORTION OF CIRCUMFERENCE OF A CIRCLE WHICH LIES IN THE
C  PLOT) USING DIGGLE'S (1983) EDGE CORRECTION FACTORS.
C----------
      INCLUDE 'RIPKCOM.FOR'

      MAXND = 0
      MINND = 1000
      AVGND = 0.0

      DO 20 K = 1,NINC
      RIPLYK(K) = 0.0
   20 CONTINUE

      DO 190 I = 1,NFROM
      NDIST = 0
C----------
C COMPUTE DISTANCE FROM CURRENT POINT TO FOUR BOUNDARIES, THEN FIND NEAREST 
C TWO BOUNDARIES. ORDER OF BOUNDARY LINES: K=1(W), K=2(N),K=3(E),K=4(S)
C----------
      DB(1) = XFROM(I) - CX(1)
      DB(2) = CY(2) - YFROM(I)
      DB(3) = CX(4) - XFROM(I)
      DB(4) = YFROM(I) -CY(1)

      DBMIN1(I)=AMIN1(DB(1),DB(3))
      DBMIN2(I)=AMIN1(DB(2),DB(4))
      DBSQ(I) = DBMIN1(i)*DBMIN1(i) + DBMIN2(i)*DBMIN2(i)

C----------
C  COMPUTE DISTANCES BETWEEN CURRENT 'FROM' POINT AND ALL 'TO' POINTS.  
C  SAVE ONLY THOSE LESS THAN THE 'DISTMAX', THE MAXIMUM DISTANCE OF 
C  INTEREST IN THE RUN.
C----------
      DO 90 J = 1,NTO
      DIST(J) = 0.0
      D = (XFROM(I)-XTO(J))**2 + (YFROM(I)-YTO(J))**2

      IF (D .GT. DMAX2) GO TO 90

      NDIST = NDIST + 1   
      DIST(NDIST) = SQRT(D)
   90 CONTINUE

      IF (NDIST .GT. MAXND) MAXND = NDIST
      IF (NDIST .LT. MINND) MINND = NDIST
      AVGND = AVGND + NDIST
      IF ((I .EQ. 1) .OR. (I .EQ. NFROM)) THEN 
      WRITE (7,6090)  I,(DIST(J),J=1,20)
 6090 FORMAT (/I6,1X,10F7.2 / 7X,10F7.2)
      ENDIF

      DLIM = 0.0-DINC
      DO 185 K = 1,NINC
      DLIM = DLIM + DINC 

      DO 180 J = 1,NDIST
      IF (NDIST .EQ. 0) GO TO 185
C  SKIP SELF (I=J) FOR WITHIN-CLASS ANALYSIS
      IF (DIST(J) .EQ. 0.0) GO TO 180

C  SKIP IF DIST > DLIM
      IF (DIST(J) .GT. DLIM) GO TO 180

      DLIM2 = DIST(J)*DIST(J)

C  COMPUTE DIGGLE'S EDGE CORRECTION FOR DISTANCE BETWEEN POINTS GREATER 
C  THAN DISTANCE TO NEAREST BOUNDARY (BUT NOT BOTH BOUNDARIES). 
C  ALSO ACCOUNTS FOR CASE WHEN NO EDGE CORRECTION REQUIRED, FOR DISTANCE 
C  BETWEEN POINTS <= DISTANCE TO NEAREST BOUNDARY. 
      IF (DLIM2 .LE. DBSQ(I)) THEN
      DB1 = AMIN1(DBMIN1(I),DIST(J))/DIST(J) 
      DB2 = AMIN1(DBMIN2(I),DIST(J))/DIST(J) 
      WEIGHT = 1.0/(1.0-((ACOS(DB1)+ACOS(DB2))/3.14159))

C  COMPUTE DIGGLE'S EDGE CORRECTION FOR DISTANCE BETWEEN POINTS GREATER THAN 
C  both DISTANCE TO VERTICAL BOUNDARY AND DISTANCE TO HORIZONTAL BOUNDARY.

      ELSE
      DB1 = DBMIN1(I)/DIST(J)
      DB2 = DBMIN2(I)/DIST(J)
      WEIGHT = 1.0/(.75-((ACOS(DB1)+ACOS(DB2))/6.283185308))
      ENDIF

      IF (WEIGHT .GT. 4.0) WEIGHT = 4.0
      RIPLYK(K) = RIPLYK(K) + WEIGHT

  180 CONTINUE
  185 CONTINUE

  190 CONTINUE
      AVGND = AVGND/NFROM
      WRITE (7,6040) MINND,MAXND,AVGND
 6040 FORMAT ( /' Number of points within ''DISTMAX'' of any ''FROM''',
     &         ' Point'/
     &'   Min= ',i5,'   Max= ',i5,'   Average= ',f10.1)

      DO 200 K=1,NINC
      IF (BIFLAG .EQ. 1) THEN
      RIPK1(K) = RIPLYK(K)
      ELSE
      RIPK2(K) = RIPLYK(K)
      ENDIF
  200 CONTINUE

      RETURN
      END

C**********************************************************************

      SUBROUTINE RKOUT
C----------
C  COMPUTE RIPLEY'S K, DECIDING IF UNIVARIATE OR BIVARIATE ANALYSIS.  
C  COMPUTE LINEARIZED STATISTICS, CONFIDENCE ENVELOPES, AND OUTPUT.
C----------
      INCLUDE 'RIPKCOM.FOR'

      WRITE (8,'(A80)' ) NAME8
      WRITE (8,'(A80)') TITLE
      WRITE (7,7001)
 7001 FORMAT (/ 1X,79('='))
      
      CRIT05 = 0.0
      RIPUP = 0.0
      RIPLOW = 0.0
      IF (BIFLAG .NE. 1) GO TO 10
C----------
C  COMPUTE THE CRITICAL VALUE FOR THE 5% ACCEPTANCE LEVEL (RIPLEY 1979)
C  FOR WITHIN-CLASS ANALYSES, WHEN 'FROM' DATASET='TO' DATASET.  
C  OTHERWISE, MUST GENERATE MONTE CARLO CONFIDENCE ENVELOPES.
C----------
      CRIT05 = 1.42*SQRT(AREA)/FLOAT(NTO)
      WRITE (7,6000) CRIT05
 6000 FORMAT (/'     Critical Value', 
     & ' = 1.42*SQRT(AREA)/N   (0.05 Significance Level)' /
     & '     CRIT(.05) =', F7.3,
     & '   *** Use for Univariate Analysis Only *** ')
      RIPUP = CRIT05
      RIPLOW = -CRIT05
   10 CONTINUE

      WRITE (7,6006) IALPHA
 6006 FORMAT (/ 
     & ' Linearized Ripley''s K Statistic -- L(d) = SQRT[K(d)/pi] - d'/  
     & '       Sign.= -1: Uniform            Sign.= 1: Clustered'//
     & '                      ',I3,'% Confidence envelope    '/
     & '      d        L(d)   ',
     & ' L-lower   L-upper  Sign. '/)

C----------                                                
C  COMPUTE RIPLEY'S K USING PLOT AREA AND POINT DENSITY CONSTANTS.
C----------
      DLIM = 0.0 - DINC
      DO 200 K = 1,NINC
      DLIM = DLIM + DINC 

C  IF MONTE CARLO ANALYSIS, COMPUTE UPPER AND LOWER ENVELOPES
      IF (MCFLAG .EQ. 1) THEN
        IF (BIFLAG .EQ. 1) THEN
        RIPLOW = SQRT((RKLOW1(K)*AREA/(N1*N2))/3.14159265) - DLIM 
        RIPUP  = SQRT((RKUP1(K) *AREA/(N1*N2))/3.14159265) - DLIM 
      ELSE
        RIPLOW = SQRT((AREA/(N1+N2))*
     &           (RKLOW1(K)/N1 + RKLOW2(K)/N2)/3.14159265)-DLIM
        RIPUP  = SQRT((AREA/(N1+N2))*
     &           (RKUP1 (K)/N1 + RKUP2 (K)/N2)/3.14159265)-DLIM
      ENDIF
      ENDIF

C  COMPUTE THE LINEARIZED VARIANT, L(d)
C  UNIVARIATE CASE
      IF (BIFLAG .EQ. 1) THEN
        RIPL = SQRT((RIPLYK(K)*AREA/(N1*N2))/3.14159265) - DLIM
        ISIGN = 0
        IF (RIPL .GT. RIPUP)  ISIGN = 1
        IF (RIPL .LT. RIPLOW) ISIGN =-1

C  BIVARIATE CASE
      ELSE 
        RIPLYK(K) = (AREA/(N1+N2))*(RIPK1(K)/N1 + RIPK2(K)/N2)
        RIPL = SQRT(RIPLYK(K)/3.14159265) - DLIM
        ISIGN = 0
        IF (MCFLAG .NE. 0) THEN
        IF (RIPL .GT. RIPUP)  ISIGN = 1
        IF (RIPL .LT. RIPLOW) ISIGN =-1
      ENDIF
      ENDIF

      WRITE (7,6004) DLIM,RIPL,RIPLOW,RIPUP,ISIGN
      WRITE (8,6004) DLIM,RIPL,RIPLOW,RIPUP,ISIGN
 6004 FORMAT (F10.3,3F10.3,I5)
  200 CONTINUE

      RETURN
      END

C***********************************************************************

      SUBROUTINE UNIMC(NITER)
C----------
C  COMPUTE CONFIDENCE ENVELOPE FOR RIPLEY'S K FROM MONTE CARLO 
C  SIMULATIONS.
C
C  FOR UNIVARIATE ANALYSIS:
C  This code draws the smallest rectangle completely surrounding 
C  the non-rectangular plot; draws NFROM random coordinate pairs on a 
C  uniform [0,1] distribution and distributes them over the rectangular
C  region; computes K(d) for the random coordinate pairs; sorts the NTIMES 
C  values of K(d) at each value of d; discards the upper and lower ALPHA 
C  percent of values, and reports the MC envelope boundaries. 
C----------

      INCLUDE 'RIPKCOM.FOR'

      DOUBLE PRECISION S0,S1
      COMMON /RANCOM/ S0,S1

C  INITIALIZE IF THIS IS THE FIRST ITERATION
      IF (NITER .GT. 1) GO TO 15

      DO 10 K = 1,NINC
      RKLOW1(K) =  100000.0
      RKUP1(K)  = -100000.0
   10 CONTINUE

C  COMPUTE RANGES OVER WHICH TO DISTRIBUTE THE RANDOM VALUES OF X AND Y
      XRNG=CX(4)-CX(1) 
      YRNG=CY(2)-CY(1) 
   15 CONTINUE

      DO 190 I = 1,NFROM
C  DRAW RANDOM COORDINATE PAIR
      CALL RANN(X)
      CALL RANN(Y)

      XFROM(I) = CX(1)+X*XRNG
      YFROM(I) = CY(1)+Y*YRNG

C COMPUTE DISTANCE FROM CURRENT POINT TO FOUR BOUNDARIES
      DB(1) = XFROM(I) - CX(1)
      DB(2) = CY(2) - YFROM(I)
      DB(3) = CX(4) - XFROM(I)
      DB(4) = YFROM(I) - CY(1)

      DBMIN1(I)=AMIN1(DB(1),DB(3))
      DBMIN2(I)=AMIN1(DB(2),DB(4))
      DBSQ(I) = DBMIN1(I)*DBMIN1(I) + DBMIN2(I)*DBMIN2(I)

  190 CONTINUE
C----------
C  POINTS ARE DRAWN, NOW DO COMPUTATIONS.
C  COMPUTE DISTANCES BETWEEN CURRENT RANDOM 'FROM' POINT AND ALL 
C  POINTS.  SAVE ONLY THOSE LESS THAN THE 'DISTMAX', THE MAXIMUM 
C  DISTANCE OF INTEREST IN THE RUN.
C----------

      DO 206 K = 1,NINC
      CARLOK(K,NITER) = 0.0
  206 CONTINUE

      MAXND = 0
      DO 95 I = 1,NFROM
      NDIST = 0
      DO 90 J = 1,NFROM
      DIST(J) = 0.0

      D = (XFROM(I)-XFROM(J))**2 + (YFROM(I)-YFROM(J))**2

      IF (D .GT. DMAX2) GO TO 90
      NDIST = NDIST + 1   
      DIST(NDIST) = SQRT(D)
   90 CONTINUE

      DLIM = 0.0 - DINC
      DO 185 K = 1,NINC
      DLIM = DLIM + DINC 

      DO 180 J = 1,NDIST
      IF (NDIST .EQ. 0) GO TO 185
      IF (DIST(J) .EQ. 0.0) GO TO 180
      IF (DIST(J) .GT. DLIM) GO TO 180

      DLIM2 = DIST(J)*DIST(J)

C  COMPUTE WEIGHTS USING DIGGLE'S EDGE CORRECTION SCHEME

      IF (DLIM2 .LE. DBSQ(I)) THEN
      DB1 = AMIN1(DBMIN1(I),DIST(J))/DIST(J) 
      DB2 = AMIN1(DBMIN2(I),DIST(J))/DIST(J) 
      WEIGHT = 1.0/(1.0-((ACOS(DB1)+ACOS(DB2))/3.14159265))

      ELSE
      DB1 = DBMIN1(I)/DIST(J)
      DB2 = DBMIN2(I)/DIST(J)
      WEIGHT = 1.0/(.75-((ACOS(DB1)+ACOS(DB2))/6.283185308))
      ENDIF

      IF (WEIGHT .GT. 4.0) WEIGHT = 4.0
      CARLOK(K,NITER) = CARLOK(K,NITER) + WEIGHT

  180 CONTINUE
  185 CONTINUE

      IF (NDIST .GT. MAXND) MAXND = NDIST
   95 CONTINUE
C  END OF POINT LOOP

C  IF THIS IS THE LAST MONTE CARLO ITERATION, THEN DO THE SORTING AND 
C    COMPUTATIONS
      IF (NITER .EQ. NTIMES) THEN
      DO 300 K = 1,NINC
      DO 290 I = 1,NTIMES
      PASS(I) = CARLOK(K,I)
  290 CONTINUE
      CALL RDPSRT (NTIMES,PASS,INDEX)

      DO 291 I=1,NTIMES
      II = NTIMES-I+1
      DSORT(I) = PASS(INDEX(II))
  291 CONTINUE

      RKLOW1(K) = AMIN1(RKLOW1(K),DSORT(ILOW))
      RKUP1(K)  = AMAX1(RKUP1(K) ,DSORT(IUP))
  300 CONTINUE
      ENDIF
  
      RETURN
      END

C***********************************************************************

      SUBROUTINE BIMC(NITER)
C----------
C  This subroutine contains the code for the TOROIDAL SHIFT version of 
C  bivariate Ripley's K analysis.  In subroutine UNIMC, the confidence 
C  envelope boundary is computed from Monte Carlo draws of Poisson 
C  distributions of points, the number of random points drawn equal to 
C  the number of actual points in the observed dataset.  Then, the 
C  univariate empirical (observed) K-distribution can be compared to the 
C  confidence envelope to determine departures from randomness.  

C  In BIMC, the idea is to preserve the background spatial patterns 
C  of the two observed datasets, and to test their independence from 
C  each other.  Thus, one population (the 'FROM' points) are always held 
C  constant, fixed in their original locations.  The other population 
C  (the 'TO' points) are held in their relative spatial distribution 
C  (the points maintain the same position relative to each other), but 
C  the entire population is shifted randomly around the 'FROM' points in 
C  each Monte Carlo iteration.  The range of the shift is +/- the length 
C  of the X boundary in the X direction, and +/- the length of the Y 
C  boundary in the Y direction.  After the shift, if a point falls 
C  outside the original plot boundaries, it is toroidally shifted, or 
C  "wrapped" around the opposite boundary, so that it is relocated inside 
C  the plot. 
C----------

      INCLUDE 'RIPKCOM.FOR'

      DOUBLE PRECISION S0,S1
      COMMON /RANCOM/ S0,S1

      DIMENSION XTOS(MAXTO),YTOS(MAXTO),
     &   DBTO1(MAXTO),DBTO2(MAXTO),DBTOSQ(MAXTO)

C  INITIALIZE IF THIS IS THE FIRST ITERATION
      IF (NITER .GT. 1) GO TO 15
           
      DO 10 K = 1,NINC
      RKLOW1(K) =  100000.0
      RKUP1(K)  = -100000.0
      RKLOW2(K) =  100000.0
      RKUP2(K)  = -100000.0
   10 CONTINUE

   15 CONTINUE

      ISET = 1
      XSIGN = 1.0
      YSIGN = 1.0
C  DRAW RANDOM COORDINATE PAIR FOR SHIFTING X AND Y
      CALL RANN(X)
      CALL RANN(Y)

      CALL RANN(SS)
      IF (SS .LT. 0.5) XSIGN = -1.0
      CALL RANN(SS)
      IF (SS .LT. 0.5) YSIGN = -1.0

      DELTAX = X*XSIGN*SOUTH
      DELTAY = Y*YSIGN*WEST

      DO 190 I = 1,NTO

      X = DELTAX + XTO(I)
      Y = DELTAY + YTO(I)

C  CHECK WHETHER POINT IS INSIDE PLOT BOUNDARIES.  IF NOT, WRAP IT 
C  AROUND THE TOROID.

      IF (X .LT. CX(1)) X = CX(4) - ABS(X)
      IF (X .GT. CX(4)) X = X - CX(4)

      IF (Y .LT. CY(1)) Y = CY(2) - ABS(Y)
      IF (Y .GT. CY(2)) Y = Y - CY(2)

C----------
C  THE NEW 'TO' POINT IS WITHIN THE PLOT.  STORE IT AS 'TOS(shifted)'.
C----------  
      XTOS(I) = X
      YTOS(I) = Y

C COMPUTE DISTANCE FROM CURRENT POINT TO FOUR BOUNDARIES
      DB(1) = XTOS(I) - CX(1)
      DB(2) = CY(2) - YTOS(I)
      DB(3) = CX(4) - XTOS(I)
      DB(4) = YTOS(I) - CY(1)

      DBTO1(I)=AMIN1(DB(1),DB(3))
      DBTO2(I)=AMIN1(DB(2),DB(4))
      DBTOSQ(I) = DBTO1(I)*DBTO1(I) + DBTO2(I)*DBTO2(I)

  190 CONTINUE
C----------
C  'TOS' POINTS ARE DRAWN, NOW DO COMPUTATIONS.
C  COMPUTE DISTANCES BETWEEN ORIGINAL 'FROM' POINT AND ALL 'TOS' POINTS.  
C  SAVE ONLY THOSE LESS THAN THE 'DISTMAX', THE MAXIMUM DISTANCE OF 
C  INTEREST IN THE RUN.
C----------

  205 CONTINUE
C  ENTER HERE AFTER 'FROM' and 'TOS' POINTS ARE SWITCHED

      DO 206 K = 1,NINC
      CARLOK(K,NITER) = 0.0
  206 CONTINUE

C  SET LIMITS
C  BIVARIATE, FIRST TIME THROUGH
      IF (ISET .EQ. 1) THEN
        NFROM1 = 1  
        NFROM2 = NFROM
        NTO1 = 1
        NTO2 = NTO
C  BIVARIATE, SWITCH POINTS, SECOND TIME THROUGH
      ELSE 
        NFROM1 =  1  
        NFROM2 = NTO
        NTO1 = 1
        NTO2 = NFROM
      ENDIF

      MAXND = 0
      DO 95 I = NFROM1,NFROM2
      NDIST = 0
      DO 90 J = NTO1,NTO2
      DIST(J) = 0.0

      IF (ISET .EQ. 1) THEN
      D = (XFROM(I)-XTOS(J))**2 + (YFROM(I)-YTOS(J))**2
      ELSE
      D = (XTOS(I)-XFROM(J))**2 + (YTOS(I)-YFROM(J))**2
      ENDIF

      IF (D .GT. DMAX2) GO TO 90
      NDIST = NDIST + 1   
      DIST(NDIST) = SQRT(D)
   90 CONTINUE

      DLIM = 0.0 - DINC
      DO 185 K = 1,NINC
      DLIM = DLIM + DINC 

      DO 180 J = 1,NDIST
      IF (NDIST .EQ. 0) GO TO 185
      IF (DIST(J) .EQ. 0.0) GO TO 180
      IF (DIST(J) .GT. DLIM) GO TO 180

      DLIM2 = DIST(J)*DIST(J)

C  COMPUTE WEIGHTS USING DIGGLE'S EDGE CORRECTION SCHEME

C  BIVARIATE, FIRST TIME THROUGH
      IF (ISET .EQ. 1) THEN
        IF (DLIM2 .LE. DBSQ(I)) THEN
        DB1 = AMIN1(DBMIN1(I),DIST(J))/DIST(J) 
        DB2 = AMIN1(DBMIN2(I),DIST(J))/DIST(J) 
        WEIGHT = 1.0/(1.0-((ACOS(DB1)+ACOS(DB2))/3.14159265))
       ELSE
        DB1 = DBMIN1(I)/DIST(J)
        DB2 = DBMIN2(I)/DIST(J)
        WEIGHT = 1.0/(.75-((ACOS(DB1)+ACOS(DB2))/(2.*3.14159265)))
       ENDIF

C  BIVARIATE, SECOND TIME THROUGH AFTER POINTS SWITCHED
      ELSE 
        IF (DLIM2 .LE. DBTOSQ(I)) THEN
        DB1 = AMIN1(DBTO1(I),DIST(J))/DIST(J) 
        DB2 = AMIN1(DBTO2(I),DIST(J))/DIST(J) 
        WEIGHT = 1.0/(1.0-((ACOS(DB1)+ACOS(DB2))/3.14159265))
       ELSE
        DB1 = DBTO1(I)/DIST(J)
        DB2 = DBTO2(I)/DIST(J)
        WEIGHT = 1.0/(.75-((ACOS(DB1)+ACOS(DB2))/(2.*3.14159265)))
       ENDIF

      ENDIF

        IF (WEIGHT .GT. 4.0) WEIGHT = 4.0
        CARLOK(K,NITER) = CARLOK(K,NITER) + WEIGHT

  180 CONTINUE
  185 CONTINUE

      IF (NDIST .GT. MAXND) MAXND = NDIST
   95 CONTINUE
C  END OF POINT LOOP

C  IF THIS IS THE LAST MONTE CARLO ITERATION, THEN DO THE SORTING AND 
C    COMPUTATIONS
      IF (NITER .EQ. NTIMES) THEN
      DO 300 K = 1,NINC
      DO 290 I = 1,NTIMES
       PASS(I) = CARLOK(K,I)
  290 CONTINUE
      CALL RDPSRT (NTIMES,PASS,INDEX)

      DO 291 I=1,NTIMES
      II = NTIMES-I+1
      DSORT(I) = PASS(INDEX(II))
  291 CONTINUE

       IF (ISET .EQ. 1) THEN
       RKLOW1(K) = AMIN1(RKLOW1(K),DSORT(ILOW))
       RKUP1(K)  = AMAX1(RKUP1(K) ,DSORT(IUP))
       ELSE
       RKLOW2(K) = AMIN1(RKLOW2(K),DSORT(ILOW))
       RKUP2(K)  = AMAX1(RKUP2(K),DSORT(IUP))
       ENDIF
  300 CONTINUE
      ENDIF
  
C  SWITCH POINTS AND DO AGAIN IF BIVARIATE CASE.
      IF (ISET .EQ. 2) GO TO 1000
      ISET = 2
      GO TO 205

 1000 CONTINUE
      RETURN
      END

C**********************************************************************

      SUBROUTINE RANN(SEL)

      DOUBLE PRECISION S0,S1
      COMMON /RANCOM/ S0,S1

      S1=DMOD(16807D0*S0,2147483647D0)
      SEL=S1/2147483648D0
      S0=S1
      RETURN
      END

C**********************************************************************

      SUBROUTINE RDPSRT(N,A,INDEX)
C----------
C  THE VECTOR INDEX IS INITIALLY LOADED WITH VALUES FROM 1 TO N
C  INCLUSIVE.  **RDPSRT** REARRANGES THE ELEMENTS OF INDEX SO THAT
C  INDEX(1) IS THE SUBSCRIPT OF THE LARGEST ELEMENT IN THE VECTOR A,
C  INDEX(2) IS THE SUBSCRIPT OF THE SECOND LARGEST ELEMENT IN A,...,
C  AND INDEX(N) IS THE SUBSCRIPT OF THE SMALLEST ELEMENT IN A.  THE
C  PHYSICAL ARRANGEMENT OF THE VECTOR A IS NOT ALTERED.  THIS
C  ALGORITHM IS AN ADAPTATION OF THE TECHNIQUE DESCRIBED IN:

C       SCOWEN, R.A. 1965. ALGORITHM 271; QUICKERSORT. COMM ACM.
C                    8(11) 669-670.
C----------
      INTEGER INDEX,IPUSH,IL,IP,IU,INDIL,INDIP,INDIU,INDKL,INDKU,
     &          ITOP,JL,JU,KL,KU
      DIMENSION INDEX(1),A(1),IPUSH(33)

      DO 10 I=1,N
   10 INDEX(I)=I

C  RETURN IF FEWER THAN TWO ELEMENTS IN ARRAY A.
      IF(N.LT.2) GOTO 9999

C  BEGIN THE SORT.
      ITOP=0
      IL=1
      IU=N
   30 CONTINUE
      IF(IU.LE.IL) GO TO 40
      INDIL=INDEX(IL)
      INDIU=INDEX(IU)
      IF(IU.GT.IL+1) GO TO 50
      IF(A(INDIL).GE.A(INDIU)) GO TO 40
      INDEX(IL)=INDIU
      INDEX(IU)=INDIL
   40 CONTINUE
      IF(ITOP.EQ.0) GOTO 9998
      IL=IPUSH(ITOP-1)
      IU=IPUSH(ITOP)
      ITOP=ITOP-2
      GO TO 30
   50 CONTINUE
      IP=(IL+IU)/2
      INDIP=INDEX(IP)
      T=A(INDIP)
      INDEX(IP)=INDIL
      KL=IL
      KU=IU
   60 CONTINUE
      KL=KL+1
      IF(KL.GT.KU) GO TO 90
      INDKL=INDEX(KL)
      IF(A(INDKL).GE.T) GO TO 60
   70 CONTINUE
      INDKU=INDEX(KU)
      IF(KU.LT.KL) GO TO 100
      IF(A(INDKU).GT.T) GO TO 80
      KU=KU-1
      GO TO 70
   80 CONTINUE
      INDEX(KL)=INDKU
      INDEX(KU)=INDKL
      KU=KU-1
      GO TO 60
   90 CONTINUE
      INDKU=INDEX(KU)
  100 CONTINUE
      INDEX(IL)=INDKU
      INDEX(KU)=INDIP
      IF(KU.LE.IP) GO TO 110
      JL=IL
      JU=KU-1
      IL=KU+1
      GO TO 120
  110 CONTINUE
      JL=KU+1
      JU=IU
      IU=KU-1
  120 CONTINUE
      ITOP=ITOP+2
      IPUSH(ITOP-1)=JL
      IPUSH(ITOP)=JU
      GO TO 30
 9998 CONTINUE
 9999 CONTINUE

      RETURN
      END
C**********************************************************************

