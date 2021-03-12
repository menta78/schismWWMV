!/ ------------------------------------------------------------------- /
      MODULE W3UOSTMD
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |             Lorenzo Mentaschi     |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         08-Oct-2018 |
!/                  +-----------------------------------+
!/
!/    Aug-2018 : Origination.                        ( version 6.07 )
!/
!/    Copyright 2010 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS.
!/       No unauthorized use without permission.
!/
!  1. Purpose : Parmeterization of the unresoled obstacles
!
!  2. Subroutines and functions :
!
!      Name               Type  Scope    Description
!     ----------------------------------------------------------------
!    UOST_INITGRID        Subr. Public   allocates UOST variables on each grid, 
!                                        and loads the matrices of alpha and beta coefficient
!    UOST_SRCTRM   Subr. Public   computes the source term for a given input spectrum
!     ----------------------------------------------------------------
!
!  3. Switches :
!
!  4. Source code :
!/
!/ ------------------------------------------------------------------- /
!/

      USE DATAPOOL

      IMPLICIT NONE
  
      PUBLIC :: UOST_INIT
      PUBLIC :: UOST_SRCTRM
  
  
      PRIVATE


      TYPE GRID
        CHARACTER(LEN=256)      :: UOSTFILELOCAL, UOSTFILESHADOW
        LOGICAL, ALLOCATABLE    :: UOST_LCL_OBSTRUCTED(:,:), UOST_SHD_OBSTRUCTED(:,:)
        INTEGER*1, ALLOCATABLE  :: UOSTLOCALALPHA(:,:,:,:), UOSTLOCALBETA(:,:,:,:)
        INTEGER*1, ALLOCATABLE  :: UOSTSHADOWALPHA(:,:,:,:), UOSTSHADOWBETA(:,:,:,:)
        REAL*4, ALLOCATABLE     :: UOSTCELLSIZE(:,:,:)
        REAL                    :: UOSTABMULTFACTOR = 100
        REAL                    :: UOSTCELLSIZEFACTOR = 1000
        REAL                    :: UOSTLOCALFACTOR = 1
        REAL                    :: UOSTSHADOWFACTOR = 1
        LOGICAL                 :: UOSTENABLED = .true.
        INTEGER                 :: NXGLOB, NXLOC, NY
      END TYPE GRID



      TYPE SGRD
        INTEGER                  :: NTH, NK
        REAL(rkind), ALLOCATABLE :: TH(:), COSTH(:), SINTH(:)
      END TYPE SGRD

  
  
      TYPE UOST_SOURCETERM
        REAL, ALLOCATABLE :: COSTH(:), SINTH(:)
        REAL :: GAMMAUP = 10
        REAL :: GAMMADOWN = 20
        TYPE(GRID) :: GRD
        TYPE(SGRD) :: SGD
        CONTAINS
        !PROCEDURE, PASS, PRIVATE :: COMPUTE_PSI => UOST_SOURCETERM_COMPUTE_PSI

        !compute_ld: estimates the local dissipation (private method)
        PROCEDURE, PASS, PRIVATE :: COMPUTE_LD => UOST_SOURCETERM_COMPUTE_LD
        !compute_se: estimates the shadow effect (private method)
        PROCEDURE, PASS, PRIVATE :: COMPUTE_SE => UOST_SOURCETERM_COMPUTE_SE
        !compute: estimates the whole dissipation
        PROCEDURE, PASS :: COMPUTE => UOST_SOURCETERM_COMPUTE
      END TYPE UOST_SOURCETERM
    
      ! srctrm: global singleton source term
      CLASS(UOST_SOURCETERM), ALLOCATABLE, TARGET :: SRCTRM


      CONTAINS


!/ ------------------------------------------------------------------- /
      SUBROUTINE UOST_INIT
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |          Lorenzo Mentaschi        |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         01-Oct-2018 |
!/                  +-----------------------------------+
!/
!/    Aug-2018 : Origination.                        ( version 6.07 )
!/
!  1. Purpose : allocate the UOST variables for a given grid, and load
!               them from file
!  2. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!     FILELOCAL, string:     file from where the alpha/beta coefficient
!                            for the local dissipation are loaded 
!     FILESHADOW, string:    file from where the alpha/beta coefficient
!                            for the shadow effect are loaded 
!     LOCALFACTOR, double:   adjustment parameter for the local 
!                            dissipation alpha and beta
!     SHADOWFACTOR, double:  adjustment parameter for the shadow
!                            dissipation alpha and beta
!     ----------------------------------------------------------------
!
!  3. Called by :
!
!      Name      Type  Module     Description
!     ----------------------------------------------------------------
!      W3IOGR    Subr. W3IOGRMD   Initialization of grid objects
!     ----------------------------------------------------------------
!
!  4. Source code :
!
!/ ------------------------------------------------------------------- /
      IMPLICIT NONE
  
      TYPE(GRID), POINTER :: GRD
      TYPE(SGRD), POINTER :: SGD
      REAL :: CGMAX, MINSIZE
  
      ALLOCATE(SRCTRM)

      GRD => SRCTRM%GRD
     ! in the future, load form ad hoc namelist, as done in ww3
      GRD%UOSTFILELOCAL = 'obstructions_local.in'
      GRD%UOSTFILESHADOW = 'obstructions_shadow.in'
      GRD%UOSTLOCALFACTOR = 1
      GRD%UOSTSHADOWFACTOR = 1
      GRD%NXGLOB = NP_GLOBAL
      GRD%NXLOC = NP_RES
      GRD%NY = 1

      SGD => SRCTRM%SGD
      SGD%NK = NUMSIG
      SGD%NTH = NUMDIR
      ALLOCATE(SGD%TH(NUMDIR), SGD%COSTH(NUMDIR), SGD%SINTH(NUMDIR))
      SGD%TH = SPDIR
      SGD%COSTH = COSTH
      SGD%SINTH = SINTH
  
      ALLOCATE( GRD%UOST_LCL_OBSTRUCTED(GRD%NXLOC, GRD%NY) )
      GRD%UOST_LCL_OBSTRUCTED = .FALSE.
      ALLOCATE( GRD%UOST_SHD_OBSTRUCTED(GRD%NXLOC, GRD%NY) )
      GRD%UOST_SHD_OBSTRUCTED = .FALSE.
      ALLOCATE( GRD%UOSTCELLSIZE(GRD%NXLOC, GRD%NY, SGD%NTH) )
      GRD%UOSTCELLSIZE = 0
      ALLOCATE( GRD%UOSTLOCALALPHA(GRD%NXLOC, GRD%NY, SGD%NK, SGD%NTH) )
      GRD%UOSTLOCALALPHA = 100
      ALLOCATE( GRD%UOSTLOCALBETA(GRD%NXLOC, GRD%NY, SGD%NK, SGD%NTH) )
      GRD%UOSTLOCALBETA = 100
      ALLOCATE( GRD%UOSTSHADOWALPHA(GRD%NXLOC, GRD%NY, SGD%NK, SGD%NTH) )
      GRD%UOSTSHADOWALPHA = 100
      ALLOCATE( GRD%UOSTSHADOWBETA(GRD%NXLOC, GRD%NY, SGD%NK, SGD%NTH) )
      GRD%UOSTSHADOWBETA = 100

  
      !  loading local/shadow alpha/beta
      CALL LOAD_ALPHABETA(GRD, SGD, 100)
  
      !  warning the user that for cells too small UOST may be inaccurate
      CGMAX = 20 ! simply taking a high value for the max group velocity to give an indication of this threshold
      MINSIZE = CGMAX*MAIN%DELT/1000
      WRITE(STAT%FHNDL,*)'*** WAVEWATCH-III WARNING IN W3UOST/UOST_INITGRID'
      WRITE(STAT%FHNDL,*)'UOST:'
      WRITE(STAT%FHNDL,*)'      global time step == ', MAIN%DELT, ' s'
      WRITE(STAT%FHNDL,*)'      FOR CELLS SMALLER THAN ABOUT ', MINSIZE,    &
                ' KM UOST MAY UNDERESTIMATE THE DISSIPATION'
      WRITE(STAT%FHNDL,*)
      END SUBROUTINE UOST_INIT
  
  
!/ ------------------------------------------------------------------- /  
      SUBROUTINE UOST_SRCTRM(IP, SPEC, CG, S, D)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |          Lorenzo Mentaschi        |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         01-Oct-2018 |
!/                  +-----------------------------------+
!/
!/    Aug-2018 : Origination.                        ( version 6.07 )
!/
!  1. Purpose : estimates the UOST source term for a give spectrum
!  2. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!     IGRID, integer:        id of the actual grid
!     ----------------------------------------------------------------
!
!  3. Called by :
!
!      Name      Type  Module     Description
!     ----------------------------------------------------------------
!      W3SRCE    Subr. W3SRCEMD   Computation of the source terms
!     ----------------------------------------------------------------
!
!  4. Source code :
!
!/ ------------------------------------------------------------------- /
      IMPLICIT NONE
  
      INTEGER, INTENT(IN) :: IP ! this is the node id in the local partition
      REAL(rkind), INTENT(IN) :: SPEC(SRCTRM%SGD%NK, SRCTRM%SGD%NTH), CG(SRCTRM%SGD%NK)
      REAL(rkind), INTENT(INOUT) :: S(SRCTRM%SGD%NK, SRCTRM%SGD%NTH), D(SRCTRM%SGD%NK, SRCTRM%SGD%NTH)
      REAL(rkind)     :: U10ABS, U10DIR

      CALL SET_WIND(IP, U10ABS, U10DIR)

      CALL SRCTRM%COMPUTE(IP, 1, SPEC, CG, U10ABS, U10DIR, S, D)
      END SUBROUTINE UOST_SRCTRM  
  

!/ ------------------------------------------------------------------- /  
      SUBROUTINE LOAD_ALPHABETA(GRD, SGD, FILEUNIT)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |          Lorenzo Mentaschi        |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         01-Oct-2018 |
!/                  +-----------------------------------+
!/
!/    Aug-2018 : Origination.                        ( version 6.07 )
!/
!  1. Purpose : loads local and shadow alpha and beta from files
!  2. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!     GRD, GRID type:        object representing the spatial grid to 
!                            be loaded
!     SGD, SGRD type:        object representing the current spectral grid
!     FILEUNIT, Integer:     unit id of the input files
!     ----------------------------------------------------------------
!
!  3. Called by :
!
!      Name            Type  Module     Description
!     ----------------------------------------------------------------
!      UOST_INITGRID   Subr. W3UOSTMD   Initialization of the UOST grid 
!     ----------------------------------------------------------------
!
!  4. Source code :
!
!/ ------------------------------------------------------------------- /
      IMPLICIT NONE

      TYPE(GRID), INTENT(INOUT) :: GRD
      TYPE(SGRD), INTENT(IN) :: SGD
      INTEGER, INTENT(IN) :: FILEUNIT
      CHARACTER(256) :: FILENAME
      LOGICAL :: FILEEXISTS
      INTEGER :: JG, J, L, I, IX, IY
!/S      INTEGER, SAVE           :: IENT   = 0

!/S      CALL STRACE (IENT, 'LOAD_ALPHABETA') 
    

      ! LOADING LOCAL ALPHA/BETA
      FILENAME = GRD%UOSTFILELOCAL
      INQUIRE(FILE=FILENAME, EXIST=FILEEXISTS)
      
      J = LEN_TRIM(FILENAME)
      IF (.NOT. FILEEXISTS) THEN
        WRITE(STAT%FHNDL,*)'*** WAVEWATCH III ERROR IN W3UOST: '// & 
                'FILE '//FILENAME(:J)//' NOT FOUND. QUITTING'
        STOP 9999 
      ENDIF  
      WRITE(STAT%FHNDL,*)'FILE '//FILENAME(:J)//' FOUND.'// &  
             'LOADING UOST SETTINGS FOR GRID'
  
      CALL LOAD_ALPHABETA_FROMFILE(FILEUNIT, FILENAME(:J), GRD%NXLOC, GRD%NY, SGD%NK, SGD%NTH,&
        GRD%UOSTABMULTFACTOR, GRD%UOSTLOCALALPHA, GRD%UOSTLOCALBETA,&
        GRD%UOSTCELLSIZE, GRD%UOST_LCL_OBSTRUCTED)
  
  
      ! LOADING SHADOW ALPHA/BETA
      FILENAME = GRD%UOSTFILESHADOW
      INQUIRE(FILE=FILENAME, EXIST=FILEEXISTS)
      
      J = LEN_TRIM(FILENAME)
      IF (.NOT. FILEEXISTS) THEN
        WRITE(STAT%FHNDL,*)'*** WAVEWATCH III ERROR IN W3UOST: '// &
                'FILE '//FILENAME(:J)//' NOT FOUND. QUITTING'
        STOP 9999
      ENDIF  
      WRITE(STAT%FHNDL,*)'FILE '//FILENAME(:J)//' FOUND.'//& 
              'LOADING UOST SETTINGS FOR GRID'
  
      CALL LOAD_ALPHABETA_FROMFILE(FILEUNIT, FILENAME(:J), GRD%NXLOC, GRD%NY, SGD%NK, SGD%NTH,&
        GRD%UOSTABMULTFACTOR, GRD%UOSTSHADOWALPHA, GRD%UOSTSHADOWBETA,&
        GRD%UOSTCELLSIZE, GRD%UOST_SHD_OBSTRUCTED)
  

      END SUBROUTINE LOAD_ALPHABETA
  
!/ ------------------------------------------------------------------- /  
      SUBROUTINE LOAD_ALPHABETA_FROMFILE(FILEUNIT, FILENAME, NXLOC, NY, NK, NTH,& 
          MULTFACTOR, ALPHAMTX, BETAMTX, CELLSIZE, ISOBSTRUCTED)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |          Lorenzo Mentaschi        |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         01-Oct-2018 |
!/                  +-----------------------------------+
!/
!/    Aug-2018 : Origination.                        ( version 6.07 )
!/
!  1. Purpose : loads alpha and beta from a single obstructions file
!  2. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!     FILEUNIT, Integer:             unit of the file to be opened
!     FILENAME, string:              name of the file
!     NXLOC, NY, NK, NTH, Integer:      size of the spatial/spectral grid (local of the partition)
!     MULTFACTOR, REAL:              multiplication factor for alpha and beta:
!                                    alpha and beta should be real in [0,1]
!                                    but to save memory the are stored in Integer*1
!     ALPHAMTX, BETAMTX, Integer*1:  loaded alpha and beta spatial/spectral matrices
!     CELLSIZE, REAL, REAL:          cell size for each spectral direction,
!                                    also loaded from the file
!     ISOBSTRUCTED, LOGICAL:         matrix of logicals, indicating for each cell
!                                    if it is obstructed or not
!     ----------------------------------------------------------------
!
!  3. Called by :
!
!      Name            Type  Module     Description
!     ----------------------------------------------------------------
!      LOAD_ALPHABETA  Subr. W3UOSTMD   Initialization of the UOST grid 
!     ----------------------------------------------------------------
!
!  4. Source code :
!
!/ ------------------------------------------------------------------- /
      IMPLICIT NONE
  
      CHARACTER(*), INTENT(IN) :: FILENAME
      REAL, INTENT(IN) :: MULTFACTOR
      INTEGER, INTENT(IN) :: FILEUNIT, NXLOC, NY, NK, NTH
      INTEGER*1, INTENT(INOUT) :: ALPHAMTX(:,:,:,:), BETAMTX(:,:,:,:)
      REAL*4, INTENT(INOUT) ::  CELLSIZE(:,:,:)
      LOGICAL, INTENT(INOUT) ::  ISOBSTRUCTED(:,:)
      CHARACTER(LEN=600) :: LINE
      INTEGER :: FIOSTAT
      LOGICAL :: HEADER, FILESTART, READINGCELLSIZE, READINGALPHA, &
                &NDISLOC
      INTEGER :: IXGLOB, IXLOC, IY, SPGRDS_SIZE, IK
      REAL, ALLOCATABLE :: TRANS(:), CLSZ(:)
!/S      INTEGER, SAVE           :: IENT   = 0

!/S      CALL STRACE (IENT, 'LOAD_ALPHABETA_FROMFILE') 
    
      !  INITIALIZING LOGICALS REPRESENTING THE DIFFERENT PHASES OF THE LOAD
      FILESTART = .TRUE.
      HEADER = .TRUE.;
      READINGCELLSIZE = .FALSE.
      READINGALPHA = .FALSE.
      IK = 0
      
      ALLOCATE(TRANS(NTH))
      ALLOCATE(CLSZ(NTH))
      
      OPEN(FILEUNIT, FILE=FILENAME, STATUS='OLD', ACTION='READ')
      READ_LOOP: DO
        READ(FILEUNIT, '(A)', IOSTAT=FIOSTAT) LINE
  
        IF (FIOSTAT .NE. 0) EXIT READ_LOOP
          
        IF (LINE(1:1) .EQ. '$') CYCLE
        
        IF (FILESTART) THEN
          !  reading the first line
          READ(LINE, '(I5)') SPGRDS_SIZE
          FILESTART = .FALSE.
        ELSEIF (HEADER) THEN
          !  reading the position of an obstructed cell
          READ(LINE, *), IXGLOB, IY
          NDISLOC = ipgl(IXGLOB)%rank .EQ. MYRANK
          IXLOC = ipgl(IXGLOB)%id
          IF (NDISLOC) ISOBSTRUCTED(IXLOC, IY) = .TRUE.
        IF ((IXLOC .GT. NXLOC) .OR. (IY .GT. NY)) THEN
          WRITE(STAT%FHNDL,*) '*** WAVEWATCH III ERROR IN W3UOST: '// &
                'GRID INDICES OUT OF RANGE.'// &
                'CHECK FILE '//FILENAME
          STOP 9999
        ENDIF
          !  marking the end of the reading of the header
          HEADER = .FALSE.
          IK = 1
          READINGCELLSIZE = .TRUE.
        ELSEIF (READINGCELLSIZE) THEN
          !  reading the sizes of the cell
          READ(LINE, *) CLSZ 
          IF (NDISLOC) CELLSIZE(IXLOC, IY, :) = CLSZ
          READINGCELLSIZE = .FALSE.
          READINGALPHA = .TRUE.
        ELSE        
          READ(LINE, *) TRANS
          IF (NDISLOC) THEN
            IF (READINGALPHA) THEN
              !  reading alpha for frequency IK
              ALPHAMTX(IXLOC, IY, IK, :) = NINT(TRANS*MULTFACTOR)
            ELSE
              !  reading beta for frequency IK
              BETAMTX(IXLOC, IY, IK, :) = NINT(TRANS*MULTFACTOR)
            ENDIF
          ENDIF
          IF (IK .LT. NK) THEN
            IK = IK + 1
          ELSE IF (READINGALPHA) THEN
            !  preparing to read the next cell
            READINGALPHA = .FALSE.
            IK = 1
          ELSE
            HEADER = .TRUE.
            IK = 1
          ENDIF
        ENDIF
      ENDDO READ_LOOP
      CLOSE(FILEUNIT)
  
      DEALLOCATE(TRANS)
      DEALLOCATE(CLSZ)
    
      END SUBROUTINE LOAD_ALPHABETA_FROMFILE
  
!/ ------------------------------------------------------------------- /  

      SUBROUTINE COMPUTE_REDUCTION_PSI(U10ABS, U10DIR, CGABS, CGDIR, PSI)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |          Lorenzo Mentaschi        |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         01-Oct-2018 |
!/                  +-----------------------------------+
!/
!/    Aug-2018 : Origination.                        ( version 6.07 )
!/
!  1. Purpose : In conditions of wind sea, the effect
!               of the unresolved obstacles is reduced.
!               Here a reduction psi is computed, as a function of the
!               wave age.
!  2. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!     U10ABS, real:     absolute value of U10
!     U10DIR, real:     direction of U10
!     CGABS, real:      absolute value of the group velocity
!     CGDIR, real:      direction of the group velocity
!     PSI, real:        output psi factor
!     ----------------------------------------------------------------
!
!  3. Called by :
!
!      Name                   Type  Module     Description
!     ----------------------------------------------------------------
!      UOST_SOURCETERM_COMPUTE_LD  Subr. W3UOSTMD   Computing the local dissipation
!      UOST_SOURCETERM_COMPUTE_SE  Subr. W3UOSTMD   Computing the shadow effect
!     ----------------------------------------------------------------
!
!  4. Source code :
!
!/ ------------------------------------------------------------------- /
      IMPLICIT NONE
      
      REAL, PARAMETER :: TOLERANCE = 0.000001
      REAL(rkind), PARAMETER :: WHTHR1 = .5, WHTHR2 = 1.5
  
      REAL(rkind), INTENT(IN) :: U10ABS, U10DIR, CGABS, CGDIR
      REAL(rkind), INTENT(OUT) :: PSI
      REAL :: THDELTA, CP, WA
!/S      INTEGER, SAVE           :: IENT   = 0

!/S      CALL STRACE (IENT, 'COMPUTE_REDUCTION_PSI') 

      ! computing the wave age
      THDELTA = ABS(U10DIR - CGDIR)
      DO WHILE (THDELTA .GT. PI)
        THDELTA = THDELTA - 2*PI
      ENDDO
      THDELTA = ABS(THDELTA)
      IF (PI/2 - THDELTA .GT. TOLERANCE) THEN
        CP = CGABS*2 ! this is scrictly valid only in deep water
        WA = CP/U10ABS/COS(THDELTA)
      ELSE
        WA = 9999999 ! a very high number
      ENDIF
  
      IF (WA .LE. WHTHR1) THEN
        !  if the wave age is less that 0.5, psi = 0, i.e. 
        !  no unresolved obstacle is considered
        PSI = 0
      ELSEIF ((WA .GT. WHTHR1) .AND. (WA .LT. WHTHR2)) THEN
        !  if the wave age is between 0.5 and 1.5
        !  psi scales linearly with WA 
        PSI = (WA - WHTHR1)/(WHTHR2 - WHTHR1)
      ELSE
        !  if the wave age is greater than 1.5 psi = 1
        PSI = 1
      ENDIF

      END SUBROUTINE COMPUTE_REDUCTION_PSI

!/ ------------------------------------------------------------------- /  

      SUBROUTINE UOST_SOURCETERM_COMPUTE_LD(THIS, IX, IY, SPEC, CG, U10ABS, U10DIR, S, D)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |          Lorenzo Mentaschi        |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         01-Oct-2018 |
!/                  +-----------------------------------+
!/
!/    Aug-2018 : Origination.                        ( version 6.07 )
!/
!  1. Purpose : Method of the class UOST_SOURCETERM. 
!               Computation of the local dissipation of the spectrum
!  2. Parameters :
!   
!     Parameter list
!     ----------------------------------------------------------------
!     THIS: UOST_SOURCETERM      instance of UOST_SOURCETERM passed to the method 
!                           (compulsory in oo programming)
!     IX, IY: Integer       coordinates of the actual cell
!     SPEC: real            input spectrum
!     CG: real              group velocity
!     U10ABS: real          absolute value of U10
!     U10DIR: real          direction of U10
!     S: real               source term
!     D: real               differential of the source term over the spectrum
!     ----------------------------------------------------------------
!
!  3. Called by :
!
!      Name                   Type  Module     Description
!     ----------------------------------------------------------------
!      UOST_SOURCETERM_COMPUTE     Subr. W3UOSTMD   Computing the source term
!     ----------------------------------------------------------------
!
!  4. Source code :
!
!/ ------------------------------------------------------------------- /
      IMPLICIT NONE
      
      CLASS(UOST_SOURCETERM), INTENT(INOUT) :: THIS
      INTEGER, INTENT(IN) :: IX, IY
      REAL(rkind), INTENT(IN) :: SPEC(THIS%SGD%NK, THIS%SGD%NTH), CG(THIS%SGD%NK)
      REAL(rkind), INTENT(OUT) :: S(THIS%SGD%NK, THIS%SGD%NTH), D(THIS%SGD%NK, THIS%SGD%NTH)
      REAL(rkind), INTENT(IN) :: U10ABS, U10DIR
      
      INTEGER :: IK, ITH, NK, NTH
      REAL(rkind) :: ALPHA, BETA, CGI, CELLSIZE, SPECI, SFC
      LOGICAL :: CELLOBSTRUCTED
      REAL(rkind) :: TH, PSI
!/S      INTEGER, SAVE           :: IENT   = 0

!/S      CALL STRACE (IENT, 'UOST_SOURCETERM_COMPUTE_LD') 

      S = 0
      D = 0
  
      CELLOBSTRUCTED =  THIS%GRD%UOST_LCL_OBSTRUCTED(IX, IY)
      IF (.NOT. CELLOBSTRUCTED) RETURN    
  
      NK = THIS%SGD%NK
      NTH = THIS%SGD%NTH
      
      DO IK = 1,NK
        CGI = CG(IK)
        DO ITH = 1,NTH
  
          !  Getting alpha and beta for local dissipation
          ALPHA = THIS%GRD%UOSTLOCALALPHA(IX, IY, IK, ITH)/THIS%GRD%UOSTABMULTFACTOR
          ALPHA = MAX(MIN(ALPHA*THIS%GRD%UOSTLOCALFACTOR, 1.), 0.)
          BETA = THIS%GRD%UOSTLOCALBETA(IX, IY, IK, ITH)/THIS%GRD%UOSTABMULTFACTOR
          BETA = MAX(MIN(BETA*THIS%GRD%UOSTLOCALFACTOR, 1.), 0.)
  
          
          IF (ALPHA .EQ. 1) CYCLE
          
          !  Getting the size of the cell along direction ith
          CELLSIZE = THIS%GRD%UOSTCELLSIZE(IX, IY, ITH)*THIS%GRD%UOSTCELLSIZEFACTOR
  
          SPECI = SPEC(IK, ITH)
  
          TH = THIS%SGD%TH(ITH)
          CALL COMPUTE_REDUCTION_PSI(U10ABS, U10DIR, CG(IK), TH, PSI)
  
          IF (BETA > 0.09) THEN
            !  Computing the local dissipation for a partially obstructed cell
            SFC = - CGI/CELLSIZE * (1 - BETA)/BETA
          ELSE
            !  The cell is almost completely obstructed. 
            !  Dissipating the energy almost completely.
            SFC = - CGI/CELLSIZE * THIS%GAMMAUP
          ENDIF
          
          S(IK, ITH) = SFC * SPECI * PSI
          D(IK, ITH) = SFC * PSI
        ENDDO
      ENDDO
    
      END SUBROUTINE UOST_SOURCETERM_COMPUTE_LD

!/ ------------------------------------------------------------------- /  
      SUBROUTINE UOST_SOURCETERM_COMPUTE_SE(THIS, IX, IY, SPEC, CG, U10ABS, U10DIR, S, D)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |          Lorenzo Mentaschi        |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         01-Oct-2018 |
!/                  +-----------------------------------+
!/
!/    Aug-2018 : Origination.                        ( version 6.07 )
!/
!  1. Purpose : Method of the class UOST_SOURCETERM. 
!               Computation of the shadow dissipation of the spectrum
!  2. Parameters :
!   
!     Parameter list
!     ----------------------------------------------------------------
!     THIS: UOST_SOURCETERM      instance of UOST_SOURCETERM passed to the method 
!                           (compulsory in oo programming)
!     IX, IY: Integer       coordinates of the actual cell
!     SPEC: real            input spectrum
!     CG: real              group velocity
!     U10ABS: real          absolute value of U10
!     U10DIR: real          direction of U10
!     S: real               source term
!     D: real               differential of the source term over the spectrum
!     ----------------------------------------------------------------
!
!  3. Called by :
!
!      Name                   Type  Module     Description
!     ----------------------------------------------------------------
!      UOST_SOURCETERM_COMPUTE     Subr. W3UOSTMD   Computing the source term
!     ----------------------------------------------------------------
!
!  4. Source code :
!
!/ ------------------------------------------------------------------- /
      IMPLICIT NONE
      
      CLASS(UOST_SOURCETERM), INTENT(INOUT), TARGET :: THIS
      INTEGER, INTENT(IN) :: IX, IY
      REAL(rkind), INTENT(IN) :: SPEC(THIS%SGD%NK, THIS%SGD%NTH), CG(THIS%SGD%NK)
      REAL(rkind), INTENT(OUT) :: S(THIS%SGD%NK, THIS%SGD%NTH), D(THIS%SGD%NK, THIS%SGD%NTH)
      REAL(rkind), INTENT(IN) :: U10ABS, U10DIR
      
      INTEGER :: IK, ITH, IS
      REAL :: CGI, SPECI, SFC, CELLSIZE, &
              SFCLEFT, SFCRIGHT, SFCCENTER, THDIAG, CGDIAG, &
              ALPHASH, BETASH, GAMMMA, GG
      INTEGER :: N = 8, ITHDIAG, NK, NTH!, NX, NY
      LOGICAL :: CELLOBSTRUCTED
      REAL(rkind) :: TH, PSI
!/S      INTEGER, SAVE           :: IENT   = 0

!/S      CALL STRACE (IENT, 'UOST_SOURCETERM_COMPUTE_SE') 
 
      S = 0
      D = 0
  
      NK = THIS%SGD%NK
      NTH = THIS%SGD%NTH
     !NX = THIS%GRD%NXLOC
     !NY = THIS%GRD%NY
      
     ! IF ((IX .EQ. 1) .OR. (IX .EQ. NX) .OR. (IY .EQ. 1) .OR. (IY .EQ. NY)) RETURN
  
      CELLOBSTRUCTED = THIS%GRD%UOST_SHD_OBSTRUCTED(IX, IY)
      IF (.NOT. CELLOBSTRUCTED) RETURN
      
      DO IK=1,NK
        DO ITH=1,NTH
  
          !  Getting alpha and beta of the shadow
          ALPHASH = THIS%GRD%UOSTSHADOWALPHA(IX, IY, IK, ITH)/THIS%GRD%UOSTABMULTFACTOR
          ALPHASH = MAX(MIN(ALPHASH*THIS%GRD%UOSTSHADOWFACTOR, 1.), 0.)
          BETASH = THIS%GRD%UOSTSHADOWBETA(IX, IY, IK, ITH)/THIS%GRD%UOSTABMULTFACTOR
          BETASH = MAX(MIN(BETASH*THIS%GRD%UOSTSHADOWFACTOR, 1.), 0.)
   
          IF (ALPHASH .EQ. 1) CYCLE        
   
          !  Getting the size of the cell along direction ith
          CELLSIZE = THIS%GRD%UOSTCELLSIZE(IX, IY, ITH)*THIS%GRD%UOSTCELLSIZEFACTOR
  
          CGI = CG(IK)
  
          GG = CGI/CELLSIZE
  
          IF (ALPHASH > 0.2) THEN
            !  Computing the shadow gamma coefficient for a partially obstructed cell
            GAMMMA = (BETASH/ALPHASH - 1)
          ELSE
            !  Alpha is small. The shadow dissipates the energy almost completely
            GAMMMA = THIS%GAMMADOWN
          ENDIF
          
          TH = THIS%SGD%TH(ITH)
          !  Computing the reduction psi related with the wind component of the spectrum
          CALL COMPUTE_REDUCTION_PSI(U10ABS, U10DIR, CG(IK), TH, PSI)
          
          SFC = - GG*GAMMMA
  
          SPECI = SPEC(IK, ITH)
          S(IK, ITH) = SFC * SPECI * PSI
          D(IK, ITH) = SFC * PSI
        ENDDO
      ENDDO  
      END SUBROUTINE UOST_SOURCETERM_COMPUTE_SE
  
!/ ------------------------------------------------------------------- /  
  
      SUBROUTINE UOST_SOURCETERM_COMPUTE(THIS, IX, IY, SPEC, CG, U10ABS, U10DIR, S, D)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |          Lorenzo Mentaschi        |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         01-Oct-2018 |
!/                  +-----------------------------------+
!/
!/    Aug-2018 : Origination.                        ( version 6.07 )
!/
!  1. Purpose : Method of the class UOST_SOURCETERM. 
!               Computation of the source term
!  2. Parameters :
!   
!     Parameter list
!     ----------------------------------------------------------------
!     THIS: UOST_SOURCETERM      instance of UOST_SOURCETERM passed to the method 
!                           (compulsory in oo programming)
!     IX, IY: Integer       coordinates of the actual cell
!     SPEC: real            input spectrum
!     CG: real              group velocity
!     U10ABS: real          absolute value of U10
!     U10DIR: real          direction of U10
!     S: real               source term
!     D: real               differential of the source term over the spectrum
!     ----------------------------------------------------------------
!
!  3. Called by :
!
!      Name                   Type  Module     Description
!     ----------------------------------------------------------------
!      UOST_SRCTRM     Subr. W3UOSTMD   Computing the source term
!     ----------------------------------------------------------------
!
!  4. Source code :
!
!/ ------------------------------------------------------------------- /
      IMPLICIT NONE
      
      CLASS(UOST_SOURCETERM), INTENT(INOUT) :: THIS
      INTEGER, INTENT(IN) :: IX, IY
      REAL(rkind), INTENT(IN) :: SPEC(THIS%SGD%NK, THIS%SGD%NTH), CG(THIS%SGD%NK)
      REAL(rkind), INTENT(IN) :: U10ABS, U10DIR
  
      REAL(rkind) :: S(THIS%SGD%NK, THIS%SGD%NTH), D(THIS%SGD%NK, THIS%SGD%NTH)
      REAL(rkind) :: S_LD(THIS%SGD%NK, THIS%SGD%NTH), S_SE(THIS%SGD%NK, THIS%SGD%NTH)
      REAL(rkind) :: D_LD(THIS%SGD%NK, THIS%SGD%NTH), D_SE(THIS%SGD%NK, THIS%SGD%NTH)
 
!/S      INTEGER, SAVE           :: IENT   = 0

!/S      CALL STRACE (IENT, 'UOST_SOURCETERM_COMPUTE') 
      
      IF (.NOT. THIS%GRD%UOSTENABLED) THEN
        S = 0
        D = 0
        RETURN
      ENDIF
  
      !  Initializing the LD and SE components
      S_LD = 0
      S_SE = 0
      !  Local dissipation
      CALL THIS%COMPUTE_LD(IX, IY, SPEC, CG, U10ABS, U10DIR, S_LD, D_LD)
      !  Shadow effect
      CALL THIS%COMPUTE_SE(IX, IY, SPEC, CG, U10ABS, U10DIR, S_SE, D_SE)
      S = S_LD + S_SE
      D = D_LD + D_SE

      END SUBROUTINE UOST_SOURCETERM_COMPUTE
  
  
!/ ------------------------------------------------------------------- /  
  
      END MODULE W3UOSTMD  

!/ ------------------------------------------------------------------- /  
  
