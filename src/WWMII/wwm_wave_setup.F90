#include "wwm_functions.h"
#undef DEBUG
#define DEBUG
      SUBROUTINE COMPUTE_LH_STRESS(F_X, F_Y)
      USE DATAPOOL
      implicit none
      real(rkind), intent(out) :: F_X(MNP), F_Y(MNP)
      real(rkind) :: INPUT(MNP)
      real(rkind) :: U_X1(MNP), U_Y1(MNP)
      real(rkind) :: U_X2(MNP), U_Y2(MNP)
      integer IP, ID, IS
      REAL(rkind) :: COSE2, SINE2, COSI2, WN, ELOC
      REAL(rkind) :: ACLOC(MSC,MDC)
      DO IP = 1, MNP
        ACLOC = AC2(IP,:,:)
        DO ID = 1, MDC
          DO IS = 2, MSC
            ELOC  = DS_INCR(IS)*DDIR*(SPSIG(IS)*ACLOC(IS,ID)+SPSIG(IS-1)*ACLOC(IS-1,ID))
            COSE2 = COS(SPDIR(ID))**TWO
            SINE2 = SIN(SPDIR(ID))**TWO
            COSI2 = COS(SPDIR(ID)) * SIN(SPDIR(ID))
            WN    = CG(IP,IS) / ( SPSIG(IS)/WK(IP,IS) )
            RSXX(IP) = RSXX(IP)+( WN * COSE2 + WN - ONEHALF)*ELOC
            RSXY(IP) = RSXY(IP)+( WN * COSI2               )*ELOC
            RSYY(IP) = RSYY(IP)+( WN * SINE2 + WN - ONEHALF)*ELOC
          ENDDO
        ENDDO
      END DO
      CALL DIFFERENTIATE_XYDIR(RSXX, U_X1, U_Y1)
      CALL DIFFERENTIATE_XYDIR(RSXY, U_X2, U_Y2)
      F_X = -U_X1 - U_Y2
      !
      CALL DIFFERENTIATE_XYDIR(RSYY, U_X1, U_Y1)
!     CALL DIFFERENTIATE_XYDIR(RSXY, U_X2, U_Y2)
      F_Y = -U_Y1 - U_X2
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_DIFF(IE, I1, UGRAD, VGRAD)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, intent(in) :: IE, I1
      REAL(rkind), intent(inout) :: UGRAD, VGRAD
      REAL(rkind) :: h
#ifdef DEBUG
      REAL(rkind) :: F1, F2, F3
#endif
      integer I2, I3, IP1, IP2, IP3
      INTEGER :: POS_TRICK(3,2)
      POS_TRICK(1,1) = 2
      POS_TRICK(1,2) = 3
      POS_TRICK(2,1) = 3
      POS_TRICK(2,2) = 1
      POS_TRICK(3,1) = 1
      POS_TRICK(3,2) = 2
      I2=POS_TRICK(I1, 1)
      I3=POS_TRICK(I1, 2)
      IP1=INE(I1,IE)
      IP2=INE(I2,IE)
      IP3=INE(I3,IE)
      h=TWO*TRIA(IE)
      UGRAD=-(YP(IP3)-YP(IP2))/h
      VGRAD= (XP(IP3)-XP(IP2))/h
#ifdef DEBUG
!      F1=(XP(IP1) - XP(IP2))*UGRAD + (YP(IP1) - YP(IP2))*VGRAD
!      F2=ZERO
!      F3=(XP(IP3) - XP(IP2))*UGRAD + (YP(IP3) - YP(IP2))*VGRAD
!      WRITE(200+MyRankD,*) 'F123=', F1, F2, F3
!      FLUSH(200+MyRankD)
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE REV_IDX_IA_JA(J, IP, JP)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, intent(in) :: J
      INTEGER, intent(out) :: IP, JP
      JP=JA(J)
      DO IP=1,MNP
        IF ((J .ge. IA(IP)) .and. (J .le. IA(IP+1)-1)) THEN
          RETURN
        END IF
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAVE_SETUP_COMPUTE_SYSTEM(ASPAR, B, FX, FY)
      USE DATAPOOL
      IMPLICIT NONE
      real(rkind), intent(in)  :: FX(MNP), FY(MNP)
      real(rkind), intent(out) :: ASPAR(NNZ)
      real(rkind), intent(out) :: B(MNP)
      INTEGER :: POS_TRICK(3,2), POS_SHIFT(3,3)
      integer I1, I2, I3, IP1, IP2, IP3
      integer IDX, IDX1, IDX2, IDX3
      INTEGER IE, IP, I, J, K, IPp, JPp
      real(rkind) :: eDep, eFX, eFY, eScal, eFact, eArea
      real(rkind) :: UGRAD, VGRAD, UGRAD1, VGRAD1
      INTEGER LIDX(2), KIDX(2), jdx
      POS_TRICK(1,1) = 2
      POS_TRICK(1,2) = 3
      POS_TRICK(2,1) = 3
      POS_TRICK(2,2) = 1
      POS_TRICK(3,1) = 1
      POS_TRICK(3,2) = 2
      ASPAR=0
      B=0
      DO I=1,3
        DO J=1,3
          K= I-J+1
          IF (K .le. 0) THEN
            K=K+3
          END IF
          IF (K .ge. 4) THEN
            K=K-3
          END IF
          POS_SHIFT(I,J)=K
        END DO
      END DO
      DO I=1,3
        jdx=0
        DO IDX=1,3
          K=POS_SHIFT(I,IDX)
          IF (K .ne. I) THEN
            jdx=jdx+1
            LIDX(jdx)=IDX
            KIDX(jdx)=K
          END IF
        END DO
        POS_SHIFT(I,LIDX(1))=KIDX(2)
        POS_SHIFT(I,LIDX(2))=KIDX(1)
      END DO
      DO IE=1,MNE
        IP1=INE(1,IE)
        IP2=INE(2,IE)
        IP3=INE(3,IE)
        eFX =( FX(IP1) +  FX(IP2) +  FX(IP3))/3.0_rkind
        eFY =( FY(IP1) +  FY(IP2) +  FY(IP3))/3.0_rkind
        eDep=(DEP(IP1) + DEP(IP2) + DEP(IP3))/3.0_rkind
        eArea=TRIA(IE)
        eFact=G9*eDep*eArea
        DO I1=1,3
          I2=POS_TRICK(I1,1)
          I3=POS_TRICK(I1,2)
          IP1=INE(I1,IE)
          IP2=INE(I2,IE)
          IP3=INE(I3,IE)
          IDX1=JA_IE(I1,1,IE)
          IDX2=JA_IE(I1,2,IE)
          IDX3=JA_IE(I1,3,IE)
          CALL COMPUTE_DIFF(IE, I1, UGRAD1, VGRAD1)
          eScal=UGRAD1*eFX + VGRAD1*eFY
          B(IP1) = B(IP1) + eScal*eArea
          !
          DO IDX=1,3
            K=POS_SHIFT(I1, IDX)
            CALL COMPUTE_DIFF(IE, K, UGRAD, VGRAD)
            eScal=UGRAD*UGRAD1 + VGRAD*VGRAD1
            J=JA_IE(I1,IDX,IE)
#ifdef DEBUG
!            WRITE(200+MyRankD,*) 'UGRAD=', UGRAD, 'VGRAD=', VGRAD
!            WRITE(200+MyRankD,*) 'UGRAD1=', UGRAD1, 'VGRAD1=', VGRAD1
!            WRITE(200+MyRankD,*) 'I1=', I1, ' K=', K
!            WRITE(200+MyRankD,*) 'I1=', I1, ' IDX=', IDX, 'eScal=', eScal
!            CALL REV_IDX_IA_JA(J, IPp, JPp)
!            WRITE(200+MyRankD,*) 'IPp=', IPp, ' JPp=', JPp
!            WRITE(200+MyRankD,*) '            -  -  -  -  -'
#endif
            ASPAR(J)=ASPAR(J)+eFact*eScal
          END DO
        END DO
#ifdef DEBUG
!        WRITE(200+MyRankD,*) '--------------------------------------'
#endif
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAVE_SETUP_APPLY_PRECOND(ASPAR, TheIn, TheOut)
      USE DATAPOOL
#ifdef MPI_PARALL_GRID
      USE elfe_msgp
#endif
      IMPLICIT NONE
      REAL(rkind), intent(in) :: ASPAR(NNZ)
      REAL(rkind), intent(in) :: TheIn(MNP)
      REAL(rkind), intent(out) :: TheOut(MNP)
      integer IP, J1, J, JP, J2
      REAL(rkind) :: eCoeff
      INTEGER :: ThePrecond = 2
      IF (ThePrecond .eq. 0) THEN
        TheOut=TheIn
      END IF
      IF (ThePrecond .eq. 1) THEN
        TheOut=0
        DO IP=1,NP_RES
          J1=I_DIAG(IP)
          DO J=IA(IP),IA(IP+1)-1
            JP=JA(J)
            IF (J .eq. J1) THEN
              eCoeff=ONE/ASPAR(J)
            ELSE
              J2=I_DIAG(JP)
#ifdef DEBUG
!            WRITE(200+MyRankD,*) 'aspar(J1)=', ASPAR(J1)
!            WRITE(200+MyRankD,*) 'aspar(J2)=', ASPAR(J2)
#endif
            
              eCoeff=-ASPAR(J) /(ASPAR(J1)*ASPAR(J2))
            END IF
            TheOut(IP)=TheOut(IP) + eCoeff*TheIn(JP)
          END DO
        END DO
#ifdef MPI_PARALL_GRID
        CALL EXCHANGE_P2D(TheOut)
#endif
      END IF
      IF (ThePrecond .eq. 2) THEN
        DO IP=1,NP_RES
          J=I_DIAG(IP)
          TheOut(IP)=TheIn(IP)/ASPAR(J)
        END DO
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAVE_SETUP_SYMMETRY_DEFECT(ASPAR)
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind), intent(in) :: ASPAR(NNZ)
      REAL(rkind) :: eVal, fVal, eSum
      INTEGER IP, J, JP, J2, IPb, nbM
      eSum=ZERO
      DO IP=1,NP_RES
        DO J=IA(IP),IA(IP+1)-1
          eVal=ASPAR(J)
          JP=JA(J)
          nbM=0
          DO J2=IA(JP),IA(JP+1)-1
            IPb=JA(J2)
            IF (IPb .eq. IP) THEN
              fVal=ASPAR(J2)
              eSum = eSum + abs(eVal - fVal)
              nbM=nbM+1
            END IF
          END DO
          IF (nbM .ne. 1) THEN
            WRITE(*,*) 'IP=', IP, 'J=', J, ' nbM=', nbM
            CALL WWM_ABORT('More errors to solve')
          END IF
        END DO
      END DO
      WRITE(200 + MyRankD,*) 'Symmetry error=', eSum
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAVE_SETUP_APPLY_FCT(ASPAR, TheIn, TheOut)
      USE DATAPOOL
#ifdef MPI_PARALL_GRID
      USE elfe_msgp
#endif
      IMPLICIT NONE
      REAL(rkind), intent(in) :: ASPAR(NNZ)
      REAL(rkind), intent(in) :: TheIn(MNP)
      REAL(rkind), intent(out) :: TheOut(MNP)
      integer IP, J, JP
      REAL(rkind) :: eCoeff
      TheOut=0
      DO IP=1,NP_RES
        DO J=IA(IP),IA(IP+1)-1
          JP=JA(J)
          eCoeff=ASPAR(J)
          TheOut(IP)=TheOut(IP) + eCoeff*TheIn(JP)
        END DO
      END DO
#ifdef MPI_PARALL_GRID
      CALL EXCHANGE_P2D(TheOut)
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAVE_SETUP_SCALAR_PROD(V1, V2, eScal)
      USE DATAPOOL, only : rkind, MNP
#ifdef MPI_PARALL_GRID
      USE DATAPOOL, only : nwild_loc_res
#endif
      USE DATAPOOL, only : NP_RES
#ifdef MPI_PARALL_GRID
      USE elfe_msgp, only : myrank, comm, ierr, nproc, istatus, rtype
#endif
      implicit none
      real(rkind), intent(in) :: V1(MNP), V2(MNP)
      real(rkind), intent(inout) :: eScal
      integer IP
#ifdef MPI_PARALL_GRID
      real(rkind) :: rScal(1), lScal(1)
      integer iProc
      lScal=0
      DO IP=1,NP_RES
        lScal(1)=lScal(1) + nwild_loc_res(IP)*V1(IP)*V2(IP)
      END DO
      IF (myrank == 0) THEN
        DO iProc=2,nproc
          CALL MPI_RECV(rScal,1,rtype, iProc-1, 19, comm, istatus, ierr)
          lScal = lScal + rScal
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(lScal,1,rtype, iProc-1, 23, comm, ierr)
        END DO
      ELSE
        CALL MPI_SEND(lScal,1,rtype, 0, 19, comm, ierr)
        CALL MPI_RECV(lScal,1,rtype, 0, 23, comm, istatus, ierr)
        eScal=lScal(1)
      END IF
#else
      eScal=0
      DO IP=1,NP_RES
        eScal=eScal + V1(IP)*V2(IP)
      END DO
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAVE_SETUP_SOLVE_POISSON_NEUMANN_DIR(ASPAR, B, TheOut)
      USE DATAPOOL
      IMPLICIT NONE
      real(rkind), intent(in) :: ASPAR(NNZ)
      real(rkind), intent(in) :: B(MNP)
      real(rkind), intent(out) :: TheOut(MNP)
      real(rkind) :: V_X(MNP), V_R(MNP), V_Z(MNP), V_P(MNP), V_Y(MNP)
      real(rkind) :: uO, uN, alphaV, h1, h2
      real(rkind) :: eNorm, CritVal, beta
      integer IP, nbIter
      nbIter=0
      CritVal=SOLVERTHR
      V_X=ZERO
      V_R=B
      CALL WAVE_SETUP_APPLY_PRECOND(ASPAR, V_R, V_Z)
      V_P=V_Z
      CALL WAVE_SETUP_SCALAR_PROD(V_Z, V_R, uO)
#ifdef DEBUG
      CALL WAVE_SETUP_SCALAR_PROD(B, B, eNorm)
      WRITE(200+MyRankD,*) 'sum(V_R)=', sum(V_R)
      WRITE(200+MyRankD,*) 'sum(V_Z)=', sum(V_Z)
      WRITE(200+MyRankD,*) 'Before loop, |B|=', eNorm
      FLUSH(200+MyRankD)
#endif
      DO
        nbIter=nbIter + 1
        Print *, 'nbIter=', nbIter
#ifdef DEBUG
!        WRITE(200+MyRankD,*) 'nbIter=', nbIter
!        FLUSH(200+MyRankD)
#endif
        CALL WAVE_SETUP_APPLY_FCT(ASPAR, V_P, V_Y)
        CALL WAVE_SETUP_SCALAR_PROD(V_P, V_Y, h2)
        alphaV=uO/h2
#ifdef DEBUG
!        WRITE(200+MyRankD,*) 'sum(V_P)=', sum(V_P)
!        WRITE(200+MyRankD,*) 'sum(V_Y)=', sum(V_Y)
!        WRITE(200+MyRankD,*) 'h2=', h2
!        WRITE(200+MyRankD,*) 'alphaV=', alphaV
!        FLUSH(200+MyRankD)
#endif
        !
        DO IP=1,MNP
          V_X(IP) = V_X(IP) + alphaV * V_P(IP)
          V_R(IP) = V_R(IP) - alphaV * V_Y(IP)
        END DO
        !
        CALL WAVE_SETUP_SCALAR_PROD(V_R, V_R, eNorm)
#ifdef DEBUG
        WRITE(200+MyRankD,*) 'nbIter=', nbIter, 'eNorm=', eNorm
        FLUSH(200+MyRankD)
#endif
        IF (eNorm .le. SOLVERTHR) THEN
          EXIT
        END IF
        !
        CALL WAVE_SETUP_APPLY_PRECOND(ASPAR, V_R, V_Z)
        CALL WAVE_SETUP_SCALAR_PROD(V_Z, V_R, uN)
        !
        beta=uN/uO
        uO=uN
        !
        DO IP=1,MNP
          V_P(IP)=V_Z(IP) + beta * V_P(IP)
        END DO
      END DO
      TheOut=V_X
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef PETSC
      SUBROUTINE PETSC_SOLVE_POISSON_NEUMANN(TheInp, TheOut)
      USE DATAPOOL
      USE ELFE_GLBL, ONLY : iplg, np_global
      USE elfe_msgp, only : myrank, nproc, comm
      use elfe_glbl, only: ipgl1=> ipgl
      ! iplg1 points to elfe_glbl::ipgl because ipgl exist allreay as integer in this function
      use petscpool
      use petscsys
      use petscmat
      implicit none
      real(rkind), intent(in) :: TheInp(MNP)
      real(rkind), intent(out) :: TheOut(MNP)
      integer :: I, J
      integer :: IP, IPGL, IE, POS
      integer :: I1, I2, I3
      integer :: POS_TRICK(3,2)

      real(kind=8)  :: X(MNP)
      real(kind=8)  :: B(MNP)
      real(kind=8)  :: ASPAR(NNZ)

      REAL    ::  TIME2
      ! solver timings
      real    ::  startTime, endTime
      real, save :: solverTimeSum = 0
!
! Petsc stuff
!
      PetscInt :: ncols
      PetscInt :: eCol
      PetscScalar :: eEntry

      integer :: counter

      KSPConvergedReason reason;
      ! solver iteration
      PetscInt iteration
      integer, save  :: iterationSum = 0        

      call PetscLogStagePush(stageFill, petscErr);CHKERRQ(petscErr)
         
      iteration = 0
!
! code for ASPAR and B. Need to write it.
!


! fill the new matrix
      ASPAR_petsc = 0
      oASPAR_petsc = 0
      counter = 1
      ncols = 0
      do i = 1, NP_RES
        ncols = IA(i+1) - IA(i)
        ! this is a interface node (row). ignore it. just increase counter
        if(ALOold2ALO(i-1) .eq. -999) then
          counter = counter + ncols
          cycle
        end if
        ! insert col by col into matrix
        do j = 1, ncols
          if(CSR_App2PetscLUT(counter) == -999) then
            oASPAR_petsc(o_CSR_App2PetscLUT(counter)) =  ASPAR(counter)
          else
            ASPAR_petsc(CSR_App2PetscLUT(counter)) = ASPAR(counter)
          endif
          counter = counter + 1
        end do
      end do
      call MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY, petscErr)
      CHKERRQ(petscErr)
      call MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY, petscErr)
      CHKERRQ(petscErr)

!fill RHS vector
!iterate over all resident (and interface) nodes
!map it to petsc global ordering
!and insert the value from B into RHS vector
      eEntry = 0;
      call VecSet(myB, eEntry, petscErr);CHKERRQ(petscErr)
      do i= 1, np
        ! this is a interface node (row). ignore it. just increase counter
        if(ALOold2ALO(i-1) .eq. -999) then
          cycle
        end if
        ! map to petsc global order
        eCol = AGO2PGO(iplg(i) - 1 )
        eEntry = B(i)
        call VecSetValue(myB, eCol, eEntry, ADD_VALUES, petscErr)
        CHKERRQ(petscErr)
      end do

      call VecAssemblyBegin(myB, petscErr);CHKERRQ(petscErr);
      call VecAssemblyEnd(myB, petscErr);CHKERRQ(petscErr);

      ! Copy the old solution from AC2 to myX to make the solver faster
      do i = 1, np
        eCol = AGO2PGO(iplg(i)-1)
        eEntry = AC2(i, ISS, IDD)
        call VecSetValue(myX,eCol,eEntry,INSERT_VALUES,petscErr)
        CHKERRQ(petscErr)
      end do
      call VecAssemblyBegin(myX, petscErr);CHKERRQ(petscErr);
      call VecAssemblyEnd(myX, petscErr);CHKERRQ(petscErr);

      ! Solve
      ! To solve successive linear systems that have different preconditioner matrices (i.e., the matrix elements
      ! and/or the matrix data structure change), the user must call KSPSetOperators() and KSPSolve() for each
      ! solve.
      if(samePreconditioner .eqv. .true.) call KSPSetOperators(Solver, matrix, matrix, SAME_PRECONDITIONER, petscErr);CHKERRQ(petscErr)
      call PetscLogStagePop(petscErr);CHKERRQ(petscErr)
      call PetscLogStagePush(stageSolve, petscErr);CHKERRQ(petscErr)
      call CPU_TIME(startTime)
      ! Solve!
      call KSPSolve(Solver, myB, myX, petscErr);CHKERRQ(petscErr);
      call CPU_TIME(endTime)
      call PetscLogStagePop(petscErr);CHKERRQ(petscErr)
         
      call KSPGetConvergedReason(Solver, reason, petscErr);CHKERRQ(petscErr);
      if (reason .LT. 0) then
        !CALL WWM_ABORT('Failure to converge')
        !write(stat%fhndl,*) 'Failure to converge'
      endif

#ifdef PETSC_DEBUG
      if(rank == 0) then
        if(reason .LT. 0 ) then
          write(DBG%FHNDL,*) "Failure to converge\n"
        else
          call KSPGetIterationNumber(Solver, iteration, petscErr)
          CHKERRQ(petscErr)
          ! print only the mean number of iteration
          iterationSum = iterationSum + iteration
          solverTimeSum = solverTimeSum + (endTime - startTime)
          if(ISS == MSC .and. IDD == MDC) then
            write(DBG%FHNDL,*) "mean number of iterations", iterationSum / real((MSC*MDC))
            print '("solver Time for all MSD MDC= ",f6.3," sec")', solverTimeSum
          endif
        endif
      endif
#endif

      X = 0.0_rkind
      !get the solution back to fortran.
      !iterate over all resident nodes (without interface and ghost nodes)
      !map the solution from petsc local ordering back to app old local ordering
      !(the app old ordering contains interface nodes)
      call VecGetArrayF90(myX, myXtemp, petscErr); CHKERRQ(petscErr)
      do i = 1, nNodesWithoutInterfaceGhosts
        X(ipgl1((PGO2AGO(PLO2PGO(i-1)))+1)%id) = myXtemp(i)
      end do
      call VecRestoreArrayF90(myX, myXtemp, petscErr)
      CHKERRQ(petscErr);
      !IF (SUM(X) .NE. SUM(X)) CALL WWM_ABORT('NaN in X')
      ! we have to fill the ghost and interface nodes with the solution from the other threads
      ! at least subroutine SOURCETERMS() make calculations on interface/ghost nodes which are
      ! normally set to 0, because they do net exist in petsc
      !call exchange_p2d(X)
      TheOut = X
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_WAVE_SETUP
# ifdef MPI_PARALL_GRID
      IF (AMETHOD ne 4) THEN
        CALL PETSC_INIT_PARALLEL
      END IF
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE FINALIZE_WAVE_SETUP
# ifdef MPI_PARALL_GRID
      IF (AMETHOD ne 4) THEN
        CALL PETSC_FINALIZE_PARALLEL
      END IF
# endif
      END SUBROUTINE
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAVE_SETUP_COMPUTATION
      USE DATAPOOL
      implicit none
      REAL(rkind) :: F_X(MNP), F_Y(MNP)
      REAL(rkind) :: ASPAR(NNZ), B(MNP)
#ifdef DEBUG
      REAL(rkind) :: Xtest(MNP), Vimg(MNP)
      REAL(rkind) :: eResidual, eResidual2
#endif
#ifdef DEBUG
      WRITE(200 + MyRankD,*) 'WAVE_SETUP_COMPUTATION, step 1'
      FLUSH(200 + MyRankD)
#endif
      CALL COMPUTE_LH_STRESS(F_X, F_Y)
      FLUSH(200 + MyRankD)
#ifdef DEBUG
      WRITE(200 + MyRankD,*) 'WAVE_SETUP_COMPUTATION, step 2'
      FLUSH(200 + MyRankD)
#endif
      CALL WAVE_SETUP_COMPUTE_SYSTEM(ASPAR, B, F_X, F_Y)
#ifdef DEBUG
      Xtest=ONE
      CALL WAVE_SETUP_APPLY_FCT(ASPAR, Xtest, Vimg)
      CALL WAVE_SETUP_SCALAR_PROD(Vimg, Vimg, eResidual)
      CALL WAVE_SETUP_SCALAR_PROD(Xtest, B, eResidual2)
      WRITE(200 + MyRankD,*) 'sum(abs(ASPAR))=', sum(abs(ASPAR))
      WRITE(200 + MyRankD,*) 'eResidual=', eResidual
      WRITE(200 + MyRankD,*) 'eResidual2=', eResidual2
      WRITE(200 + MyRankD,*) 'WAVE_SETUP_COMPUTATION, step 3'
      CALL WAVE_SETUP_SYMMETRY_DEFECT(ASPAR)
      FLUSH(200 + MyRankD)
#endif
      CALL WAVE_SETUP_SOLVE_POISSON_NEUMANN_DIR(ASPAR, B, ZETA_SETUP)
#ifdef DEBUG
      WRITE(200 + MyRankD,*) 'WAVE_SETUP_COMPUTATION, step 4'
      FLUSH(200 + MyRankD)
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
