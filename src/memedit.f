	SUBROUTINE MEMEDIT( MODE, LEVEL, DEF, AIM, RMAX, ACC )
C
C  PROGRAM FOR INTERACTIVELY SELECTING MEM ITERATION PARAMETERS
C
C  JULY 83 BY KDH AT IOA
* Dec 85 KDH at IOA - revised for new release of MEMSYS

	CHARACTER*1 PARAM

    1	TYPE *,' '
	TYPE *,'MEM iteration parameters:'
	TYPE *,'        MODE',MODE
	TYPE *,'       LEVEL',LEVEL
	TYPE *,'         AIM',AIM
	TYPE *,'        RMAX',RMAX
	TYPE *,'     DEFAULT',DEF
	TYPE *,' '

* decode mode

	M4=MODE/10000
	M3=(MODE-10000*M4)/1000
	M2=(MODE-10000*M4-1000*M3)/100
	M1=(MODE-10000*M4-1000*M3-100*M2)/10
	M0=(MODE-10000*M4-1000*M3-100*M2-10*M1)

* report search directions

    5 IF(M0.EQ.0) THEN
	TYPE *,'      full recursive generation of search directions'
      ELSE IF(M0.EQ.3) THEN
	TYPE *,'      3 search directions (plus current map)'
      ELSE IF(M0.EQ.4) THEN
	TYPE *,'      4 search directions (plus current map)'
      ELSE IF(M0.EQ.6) THEN
	TYPE *,'      6 search directions (plus current map)'
      ELSE
	TYPE *,CHAR(7),'INVALID NUMBER OF SEARCH DIRECTIONS', M0
	M0 = 0
	GOTO 5
      END IF

* report type of default

   10 IF(M1.EQ.0) THEN
	TYPE *,'      default =', DEF
      ELSE IF(M1.EQ.1) THEN
	TYPE *,'      default = exp(total(plnf))'
      ELSE IF(M1.EQ.3) THEN
	TYPE *,'      default = user image'
      ELSE IF(M1.EQ.4) THEN
	TYPE *,'      default =', DEF, ' with 0<f<1'
      ELSE IF(M1.EQ.5) THEN
	TYPE *,'      default =',DEF,' with total(f) held fixed'
      ELSE IF(M1.EQ.6) THEN
	TYPE *,'      default = exp(total(plnf)) with total(f) held fixed'
      ELSE IF(M1.EQ.8) THEN
	TYPE *,'      default = user image with total(f) held fixed'
      ELSE
	TYPE *,CHAR(7),'INVALID DEFAULT OPTION', M1
	M1 = 1
	GOTO 10
      END IF

* report constraint option 

   15 IF(M2.EQ.0) THEN
	TYPE *,'      chi-squared =', AIM
      ELSE IF(M2.EQ.1) THEN
	TYPE *,'      chi-squared = ', AIM, ' with sigma =', ACC
      ELSE IF(M2.EQ.5) THEN
	TYPE *,'      chi-squared rejecting outliers = ', AIM
      ELSE IF(M2.EQ.6) THEN
	TYPE *,'      chi-squared rejecting outliers = ', AIM, ' sigma =', ACC
      ELSE
	TYPE *,CHAR(7),'** INVALID CONSTRAINT OPTION', M2
	M2 = 1
	GOTO 15
      END IF

* report weight option

   20 IF(M3.EQ.0) THEN
	TYPE *,'      equal pixel weights'
      ELSE IF(M3.EQ.1) THEN
	TYPE *,'      user supplied pixel weights'
      ELSE
	TYPE *,CHAR(7),'** INVALID WEIGHT OPTION', M3
	M3 = 0
	GOTO 20
      END IF

	MODE = 1000*M3 + 100*M2 + 10*M1 + M0

* command node

	TYPE *,' '
	TYPE '(A)', '$What parameter to change (<CR> to proceed) ? '
	READ(*,'(A)') PARAM

	IF( PARAM .EQ. ' ' ) RETURN

      IF( PARAM .EQ. 'M' ) THEN
	TYPE *,'old MODE : ', MODE, ' enter new MODE : '
	READ(*,*,ERR=777) NEW
	MODE=NEW

      ELSE IF( PARAM .EQ. 'L' ) THEN
	TYPE *,'old LEVEL : ', LEVEL,' enter new LEVEL : '
	READ(*,*,ERR=777) NEW
	LEVEL=NEW

      ELSE IF( PARAM .EQ. 'D' ) THEN
	TYPE *,'old DEFAULT : ', DEF,' enter new DEFAULT : '
	READ(*,*,ERR=777) VALUE
	IF(VALUE.LE.0.) GOTO 777
	DEF=VALUE

      ELSE IF( PARAM .EQ. 'A' ) THEN
	TYPE *, 'old AIM : ', AIM, ' enter new AIM : '
	READ(*,*,ERR=777) VALUE
	IF(VALUE.LE.0.) GOTO 777
	AIM=VALUE

      ELSE IF( PARAM .EQ. 'R' ) THEN
	TYPE *,'old RMAX : ', RMAX, ' enter new RMAX : '
	READ(*,*,ERR=777) VALUE
	IF(VALUE.LE.0.) GOTO 777
	RMAX=VALUE

      END IF
	GOTO   1

  777	TYPE *, CHAR(7), '** INVALID REPLY, OLD VALUE KEPT.'
	GOTO   1
      END
