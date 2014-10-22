#include "wwm_functions.h"
MODULE WWM_PGMCL_COUPLING_WITH_ROMS
      LOGICAL :: L_FIRST_ORDER_ARDHUIN
      LOGICAL :: L_STOKES_DRIFT_USING_INTEGRAL
      integer NlevelVert
      integer NlevelPartial
      integer NlevelIntegral
      CONTAINS
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef ROMS_WWM_PGMCL_COUPLING
# ifdef WWM_MPI
      SUBROUTINE WWM_CreateMatrixPartition
      USE DATAPOOL
      USE mod_coupler
      implicit none
      integer, allocatable :: rbuf_int(:)
      integer, allocatable :: TheIndex(:)
      integer, allocatable :: NumberNode(:), NumberTrig(:)
      integer, allocatable :: All_LocalToGlobal(:,:)
      integer i, eIdx, iProc, MNPloc, MNEloc, idx
      integer IPc, IP
#  ifdef DEBUG_WWM
      integer MinValIndex, MinValIndexInv, eVal
#  endif
      allocate(MatrixBelongingWAV(np_global, NnodesWAV), NumberNode(NnodesWAV), NumberTrig(NnodesWAV), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_coupl_roms, allocate error 3')
      IF (myrank.ne.MyRankLocal) THEN
        CALL WWM_ABORT('die from ignominious death')
      END IF
      IF (MyRankLocal.eq.0) THEN
        allocate(All_LocalToGlobal(np_global, NnodesWAV), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_coupl_roms, allo error 4')
        All_LocalToGlobal=0
        MatrixBelongingWAV=0
        DO i=1,MNP
          eIdx=iplg(i)
          MatrixBelongingWAV(eIdx,1)=i
          All_LocalToGlobal(i,1)=eIdx
        ENDDO
        NumberNode(1)=MNP
        DO iProc=2,NnodesWAV
          allocate(rbuf_int(1), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_coupl_roms, alloc error 5')
          CALL MPI_RECV(rbuf_int,1,itype, iProc-1, 194, WAV_COMM_WORLD, istatus, ierr)
          MNPloc=rbuf_int(1)
          NumberNode(iProc)=MNPloc
          deallocate(rbuf_int)
          !
          allocate(rbuf_int(MNPloc), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_coupl_roms, alloc error 6')
          CALL MPI_RECV(rbuf_int,MNPloc,itype, iProc-1, 195, WAV_COMM_WORLD, istatus, ierr)
          DO IP=1,MNPloc
            eIdx=rbuf_int(IP)
            MatrixBelongingWAV(eIdx,iProc)=IP
            All_LocalToGlobal(IP,iProc)=eIdx
          END DO
          deallocate(rbuf_int)
        END DO
        !
        allocate(rbuf_int(np_global*NnodesWAV), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_coupl_roms, alloc error 7')
        idx=0
        DO iProc=1,NnodesWAV
          DO IP=1,np_global
            idx=idx+1
            rbuf_int(idx)=MatrixBelongingWAV(IP,iProc)
          END DO
        END DO
        DO iProc=2,NnodesWAV
          CALL MPI_SEND(rbuf_int,np_global*NnodesWAV,itype, iProc-1, 196, WAV_COMM_WORLD, ierr)
        END DO
        deallocate(rbuf_int)
      ELSE
        allocate(rbuf_int(1), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_coupl_roms, alloc error 8')
        rbuf_int(1)=MNP
        CALL MPI_SEND(rbuf_int,1,itype, 0, 194, WAV_COMM_WORLD, ierr)
        deallocate(rbuf_int)

        CALL MPI_SEND(iplg,MNP,itype, 0, 195, WAV_COMM_WORLD, ierr)
!
        allocate(rbuf_int(np_global*NnodesWAV), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_coupl_roms, alloc error 10')
        CALL MPI_RECV(rbuf_int,np_global*NnodesWAV,itype, 0, 196, WAV_COMM_WORLD, istatus, ierr)
        idx=0
        DO iProc=1,NnodesWAV
          DO i=1,np_global
            idx=idx+1
            MatrixBelongingWAV(i,iProc)=rbuf_int(idx)
          END DO
        END DO
        deallocate(rbuf_int)
      ENDIF
      deallocate(NumberNode, NumberTrig)
      END SUBROUTINE
# else
      SUBROUTINE WWM_CreateMatrixPartition
      USE DATAPOOL
      USE mod_coupler
      IMPLICIT NONE
      integer IP
      allocate(MatrixBelongingWAV(MNP, 1), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_coupl_roms, allocate error 12')
      DO IP=1,MNP
        MatrixBelongingWAV(IP,1)=IP
      ENDDO
      END SUBROUTINE
# endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ROMS_COUPL_INITIALIZE
      USE DATAPOOL
      USE mod_coupler
      USE PGMCL_LIBRARY
      USE pgmcl_interp
      implicit none
      logical DoNearest
      integer rbuf_int(1)
      integer IP, iNodeSel, idx, eRankRecv
      real(rkind) eDiff, AbsDiff, SumDep1, SumDep2, SumDiff
      real(rkind) minBathy, maxBathy
!        character(len=40) :: FileSave1, FileSave2
!        character(len=3) :: eStrFi
      real(rkind) SumDepReceive
      !
      ! First part: initializations of the code
      !
# ifdef DEBUG_WWM
      WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 1, rnk=', myrank
      FLUSH(DBG%FHNDL)
# endif
      CALL SetComputationalNodes(ArrLocal, NnodesWAV, OCNid)
# ifdef DEBUG_WWM
      WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 1.2, rnk=', myrank
      FLUSH(DBG%FHNDL)
# endif
      CALL WWM_CreateMatrixPartition
# ifdef DEBUG_WWM
      WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 1.3, rnk=', myrank
      FLUSH(DBG%FHNDL)
# endif
      CALL GET_GRID_ARRAY_FE_r8(NP_TOTAL, NE_TOTAL, XPtotal, YPtotal, INEtotal, eGrid_wav)
# ifdef DEBUG_WWM
      WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 2, rnk=', myrank
      FLUSH(DBG%FHNDL)
# endif
      !
      ! Second part: exchanging grids
      !
      IF (MyRankLocal.eq.0) THEN
        CALL M2M_send_grid(ArrLocal, OCNid, eGrid_wav)
        CALL M2M_send_node_partition(ArrLocal, OCNid,                 &
     &        np_total, NnodesWAV, MatrixBelongingWAV)
      ENDIF
# ifdef DEBUG_WWM
      WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 3, rnk=', myrank
      FLUSH(DBG%FHNDL)
# endif
      CALL M2M_recv_grid(ArrLocal, OCNid, eGrid_ocn_rho)
# ifdef DEBUG_WWM
      WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 4, rnk=', myrank
      FLUSH(DBG%FHNDL)
# endif
      CALL M2M_recv_grid(ArrLocal, OCNid, eGrid_ocn_u)
# ifdef DEBUG_WWM
      WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 5, rnk=', myrank
      FLUSH(DBG%FHNDL)
# endif
      CALL M2M_recv_grid(ArrLocal, OCNid, eGrid_ocn_v)
# ifdef DEBUG_WWM
      WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 6, rnk=', myrank
      FLUSH(DBG%FHNDL)
# endif
      CALL M2M_recv_node_partition(ArrLocal, OCNid,                     &
     &   NnodeRho, NnodesOCN, MatrixBelongingOCN_rho)
# ifdef DEBUG_WWM
      WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 7, rnk=', myrank
      FLUSH(DBG%FHNDL)
# endif
      CALL M2M_recv_node_partition(ArrLocal, OCNid,                     &
     &   NnodeU, NnodesOCN, MatrixBelongingOCN_u)
# ifdef DEBUG_WWM
      WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 8, rnk=', myrank
      FLUSH(DBG%FHNDL)
# endif
      CALL M2M_recv_node_partition(ArrLocal, OCNid,                     &
     &   NnodeV, NnodesOCN, MatrixBelongingOCN_v)
# ifdef DEBUG_WWM
      WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 9, rnk=', myrank
      FLUSH(DBG%FHNDL)
# endif
      eRankRecv=ArrLocal % ListFirstRank(OCNid)
      CALL MPI_RECV(rbuf_int,4,itype, eRankRecv, 103, MPI_COMM_WORLD, istatus, ierr)
      Nlevel=rbuf_int(1)
      NlevelVert=rbuf_int(2)
      IF (rbuf_int(3) .eq. 1) THEN
        L_FIRST_ORDER_ARDHUIN=.TRUE.
        NlevelPartial=1
      ELSE
        L_FIRST_ORDER_ARDHUIN=.FALSE.
        NlevelPartial=Nlevel
      END IF
      IF (rbuf_int(4) .eq. 1) THEN
        L_STOKES_DRIFT_USING_INTEGRAL=.TRUE.
        NlevelIntegral=Nlevel+1
      ELSE
        L_STOKES_DRIFT_USING_INTEGRAL=.FALSE.
        NlevelIntegral=1
      END IF
# ifdef DEBUG_WWM
      WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 10, rnk=', myrank
      FLUSH(DBG%FHNDL)
# endif
      DoNearest=.TRUE.
# ifdef DEBUG_WWM
      WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 11, rnk=', myrank
      FLUSH(DBG%FHNDL)
# endif
      !
      ! Third part: computing interpolation matrices
      !
      FileSave_OCNtoWAV_rho='InterpSave_OCNtoWAV_rho'
      FileSave_OCNtoWAV_u='InterpSave_OCNtoWAV_u'
      FileSave_OCNtoWAV_v='InterpSave_OCNtoWAV_v'
      CALL SAVE_CreateInterpolationSparseMatrix_Parall(                 &
     &    FileSave_OCNtoWAV_rho, mMat_OCNtoWAV_rho, DoNearest,          &
     &    eGrid_ocn_rho, eGrid_wav,                                     &
     &    WAV_COMM_WORLD, MatrixBelongingWAV)
      CALL SAVE_CreateInterpolationSparseMatrix_Parall(                 &
     &    FileSave_OCNtoWAV_u, mMat_OCNtoWAV_u, DoNearest,              &
     &    eGrid_ocn_u, eGrid_wav,                                       &
     &    WAV_COMM_WORLD, MatrixBelongingWAV)
      CALL SAVE_CreateInterpolationSparseMatrix_Parall(                 &
     &    FileSave_OCNtoWAV_v, mMat_OCNtoWAV_v, DoNearest,              &
     &    eGrid_ocn_v, eGrid_wav,                                       &
     &    WAV_COMM_WORLD, MatrixBelongingWAV)
      CALL M2M_recv_sparseMatrix(ArrLocal, OCNid, mMat_WAVtoOCN_rho)
      CALL M2M_recv_sparseMatrix(ArrLocal, OCNid, mMat_WAVtoOCN_u)
      CALL M2M_recv_sparseMatrix(ArrLocal, OCNid, mMat_WAVtoOCN_v)
      IF (MyRankLocal .eq. 0) THEN
        CALL M2M_send_sparseMatrix(ArrLocal, OCNid, mMat_OCNtoWAV_rho)
        CALL M2M_send_sparseMatrix(ArrLocal, OCNid, mMat_OCNtoWAV_u)
        CALL M2M_send_sparseMatrix(ArrLocal, OCNid, mMat_OCNtoWAV_v)
      END IF
      CALL DEALLOCATE_GRID_ARRAY(eGrid_wav)
      CALL DEALLOCATE_GRID_ARRAY(eGrid_ocn_rho)
      CALL DEALLOCATE_GRID_ARRAY(eGrid_ocn_u)
      CALL DEALLOCATE_GRID_ARRAY(eGrid_ocn_v)
      !
      ! Fourth part: Computing restricted interpolation matrices
      ! and asynchronous arrays
      !
# ifdef DEBUG_WWM
      WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 14', myrank
      FLUSH(DBG%FHNDL)
# endif
      CALL MPI_INTERP_GetSystemOutputSide(ArrLocal, OCNid, WAVid,       &
     &    MatrixBelongingOCN_rho, MatrixBelongingWAV,                   &
     &    mMat_OCNtoWAV_rho, TheArr_OCNtoWAV_rho)
      CALL MPI_INTERP_GetAsyncOutput_r8(TheArr_OCNtoWAV_rho,            &
     &    3, TheAsync_OCNtoWAV_uvz)
      CALL MPI_INTERP_GetAsyncOutput_r8(TheArr_OCNtoWAV_rho,            &
     &    Nlevel+1, TheAsync_OCNtoWAV_rho)
      CALL MPI_INTERP_GetSystemOutputSide(ArrLocal, OCNid, WAVid,       &
     &    MatrixBelongingOCN_u, MatrixBelongingWAV,                     &
     &    mMat_OCNtoWAV_u, TheArr_OCNtoWAV_u)
      CALL MPI_INTERP_GetAsyncOutput_r8(TheArr_OCNtoWAV_u,              &
     &    NlevelVert, TheAsync_OCNtoWAV_u)
      CALL MPI_INTERP_GetSystemOutputSide(ArrLocal, OCNid, WAVid,       &
     &    MatrixBelongingOCN_v, MatrixBelongingWAV,                     &
     &    mMat_OCNtoWAV_v, TheArr_OCNtoWAV_v)
      CALL MPI_INTERP_GetAsyncOutput_r8(TheArr_OCNtoWAV_v,              &
     &    NlevelVert, TheAsync_OCNtoWAV_v)
      CALL MPI_INTERP_GetSystemInputSide(ArrLocal, WAVid, OCNid,        &
     &    MatrixBelongingWAV, MatrixBelongingOCN_rho,                   &
     &    mMat_WAVtoOCN_rho, TheArr_WAVtoOCN_rho)
      CALL MPI_INTERP_GetAsyncInput_r8(TheArr_WAVtoOCN_rho,             &
     &    19, TheAsync_WAVtoOCN_stat)
      CALL MPI_INTERP_GetSystemInputSide(ArrLocal, WAVid, OCNid,        &
     &    MatrixBelongingWAV, MatrixBelongingOCN_u,                     &
     &    mMat_WAVtoOCN_u, TheArr_WAVtoOCN_u)
      CALL MPI_INTERP_GetAsyncInput_r8(TheArr_WAVtoOCN_u,               &
     &    NlevelIntegral, TheAsync_WAVtoOCN_u)
      CALL MPI_INTERP_GetSystemInputSide(ArrLocal, WAVid, OCNid,        &
     &    MatrixBelongingWAV, MatrixBelongingOCN_v,                     &
     &    mMat_WAVtoOCN_v, TheArr_WAVtoOCN_v)
      CALL MPI_INTERP_GetAsyncInput_r8(TheArr_WAVtoOCN_v,               &
     &    NlevelIntegral, TheAsync_WAVtoOCN_v)
      deallocate(MatrixBelongingWAV)
      deallocate(MatrixBelongingOCN_rho)
      deallocate(MatrixBelongingOCN_u)
      deallocate(MatrixBelongingOCN_v)
      CALL DeallocSparseMatrix(mMat_OCNtoWAV_rho)
      CALL DeallocSparseMatrix(mMat_OCNtoWAV_u)
      CALL DeallocSparseMatrix(mMat_OCNtoWAV_v)
      CALL DeallocSparseMatrix(mMat_WAVtoOCN_rho)
      CALL DeallocSparseMatrix(mMat_WAVtoOCN_u)
      CALL DeallocSparseMatrix(mMat_WAVtoOCN_v)
      !
      ! Fifth part: more allocations and exchanges
      !
      allocate(A_wav_ur_3D(NlevelVert,MNP), A_wav_vr_3D(NlevelVert,MNP), U_wav(NlevelVert,MNP), V_wav(NlevelVert,MNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_coupl_roms, allocate error 23.4')
      allocate(CosAng(MNP), SinAng(MNP), dep_rho(MNP), A_wav_rho_3D(Nlevel+1,MNP), A_wav_stat(19,MNP), A_wav_uvz(3,MNP), A_wav_rho(MNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_coupl_roms, allocate error 16.1')
      allocate(A_wav_u_3D(NlevelIntegral,MNP), A_wav_v_3D(NlevelIntegral,MNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_coupl_roms, allocate error 17')
# ifdef DEBUG_WWM
      WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 24, rnk=', myrank
      FLUSH(DBG%FHNDL)
# endif
      CALL MPI_INTERP_RECV_r8(TheArr_OCNtoWAV_rho, 23, A_wav_rho)
      DO IP=1,MNP
        CosAng(IP)=COS(A_wav_rho(IP))
        SinAng(IP)=SIN(A_wav_rho(IP))
      END DO
# ifdef DEBUG_WWM
      WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 25, rnk=', myrank
      WRITE(DBG%FHNDL,*) 'MyRankGlobal=', MyRankGlobal
      FLUSH(DBG%FHNDL)
# endif
      CALL MPI_INTERP_RECV_r8(TheArr_OCNtoWAV_rho, 217, A_wav_rho)
# ifdef DEBUG_WWM
      SumDepReceive=0
# endif
      DO IP=1,MNP
        dep_rho(IP)=A_wav_rho(IP)
# ifdef DEBUG_WWM
        SumDepReceive=SumDepReceive + abs(A_wav_rho(idx))
# endif
      END DO
# ifdef DEBUG_WWM
      WRITE(DBG%FHNDL,*) 'SumDepReceive=', SumDepReceive
      WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, WAV, step 33'
      WRITE(DBG%FHNDL,*) 'WAV, rnk=', myrank
      FLUSH(DBG%FHNDL)
      AbsDiff=0
      SumDep1=0
      SumDep2=0
      SumDiff=0
      minBathy=140000
      maxBathy=0
      DO IP=1,MNP
        IF (dep_rho(IP).lt.minBathy) THEN
          minBathy=dep_rho(IP)
        END IF
        IF (dep_rho(IP).gt.maxBathy) THEN
          maxBathy=dep_rho(IP)
        END IF
      END DO
      WRITE(DBG%FHNDL,*) 'dep_rho, min=', minBathy, ' max=', maxBathy
      FLUSH(DBG%FHNDL)
      minBathy=140000
      maxBathy=0
      DO IP=1,MNP
        IF (DEP(IP).lt.minBathy) THEN
          minBathy=DEP(IP)
        END IF
        IF (DEP(IP).gt.maxBathy) THEN
          maxBathy=DEP(IP)
        END IF
      END DO
      WRITE(DBG%FHNDL,*) 'DEP, min=', minBathy, ' max=', maxBathy
      FLUSH(DBG%FHNDL)


      iNodeSel=-1
!        CALL MyGetString(MyRankGlobal, eStrFi)
!        FileSave1='DEP_infos' // eStrFi
!        FileSave2='Lookup_infos' // eStrFi
!        open(745, FILE=TRIM(FileSave1))
!        open(746, FILE=TRIM(FileSave2))
      DO IP=1,MNP
        DEP(IP)=dep_rho(IP)
        eDiff=abs(dep_rho(IP) - DEP(IP))
        SumDiff=SumDiff + eDiff
        SumDep1=SumDep1 + dep_rho(IP)
        SumDep2=SumDep2 + DEP(IP)
        IF (eDiff.gt.AbsDiff) THEN
          AbsDiff=eDiff
          iNodeSel=IP
        END IF
        IF ((DEP(IP).ge.200).and.(eDiff.ge.10)) THEN
          WRITE(DBG%FHNDL,*) 'AD, IP=', IP, dep_rho(IP), DEP(IP)
          WRITE(DBG%FHNDL,*) 'AD, xp, yp=', XP(IP), YP(IP)
          FLUSH(DBG%FHNDL)
        END IF
      END DO
!        close(745)
!        close(746)
!        WRITE(DBG%FHNDL,*) 'AD, AbsDiff=', AbsDiff
!        WRITE(DBG%FHNDL,*) 'AD, IP=', iNodeSel, dep_rho(iNodeSel), DEP(iNodeSel)
!        WRITE(DBG%FHNDL,*) 'AD, xp, yp=', XP(iNodeSel), YP(iNodeSel)
      WRITE(DBG%FHNDL,*) 'AD, SumDep1=', SumDep1, ' SumDep2=', SumDep2
      WRITE(DBG%FHNDL,*) 'AD, SumDiff=', SumDiff
      FLUSH(DBG%FHNDL)
# endif
      allocate(z_w_wav(0:Nlevel, MNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_coupl_roms, allocate error 19')
      allocate(USTOKES_wav(Nlevel,MNP), VSTOKES_wav(Nlevel,MNP), ZETA_CORR(MNP), J_PRESSURE(MNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_coupl_roms, allocate error 20')
# ifdef DEBUG_WWM
      WRITE(DBG%FHNDL,*) 'End ROMS_COUPL_INITIALIZE'
      FLUSH(DBG%FHNDL)
# endif
      CALL DEALLOCATE_Arr(TheArr_OCNtoWAV_rho)
      CALL DEALLOCATE_Arr(TheArr_OCNtoWAV_u)
      CALL DEALLOCATE_Arr(TheArr_OCNtoWAV_v)
      CALL DEALLOCATE_Arr(TheArr_WAVtoOCN_rho)
      CALL DEALLOCATE_Arr(TheArr_WAVtoOCN_u)
      CALL DEALLOCATE_Arr(TheArr_WAVtoOCN_v)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ROMS_COUPL_DEALLOCATE
      USE DATAPOOL
      USE mod_coupler
      USE PGMCL_LIBRARY
      USE pgmcl_interp
      implicit none
      logical DoNearest
      integer IP, iNodeSel, idx, eRankRecv
      real(rkind) eDiff, AbsDiff, SumDep1, SumDep2, SumDiff
      real(rkind) minBathy, maxBathy
      real(rkind) SumDepReceive
      deallocate(CosAng, SinAng, dep_rho)
      deallocate(A_wav_rho_3D, A_wav_stat, A_wav_uvz, A_wav_rho)
      deallocate(A_wav_u_3D, A_wav_v_3D, A_wav_ur_3D, A_wav_vr_3D)
      deallocate(z_w_wav, U_wav, V_wav)
      deallocate(USTOKES_wav, VSTOKES_wav)
      deallocate(ZETA_CORR, J_PRESSURE)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE STOKES_STRESS_INTEGRAL_ROMS
      USE mod_coupler
      USE DATAPOOL
      implicit none
      integer IP, k, ID, IS
      real(rkind) eF1, eF2, eDelta, TheInt, eDep, eHeight
      real(rkind) eFrac, eFracB, eQuot, TheIntChk
      real(rkind) eFct, eQuot1, eQuot2, eQuot3, eScal, eZeta
      real(rkind) eOmega, eMult, MFACT, kD, eSinc
      real(rkind) USTOKES1, USTOKES2, USTOKES3
      real(rkind) VSTOKES1, VSTOKES2, VSTOKES3
      real(rkind) USTOKESpart, VSTOKESpart, eJPress
      real(rkind) ACLOC, eWk, eSigma, eLoc, eSinhkd, eSinh2kd, eSinhkd2
      real(rkind) zMid, HS, ETOT, MinVal_MFACT, MaxVal_MFACT, MaxVal_eQuot1
      real(rkind) MinVal_MFACT_gl, MaxVal_MFACT_gl, MaxHS, SumHS, AvgHS
      real(rkind) WLM, KLM, AvgStokesNormA, AvgStokesNormB
      real(rkind) PPTAIL, CETAIL, CKTAIL
      real(rkind) ETOT1, EKTOT
      real(rkind) eQuotDispersion, eMaxAC, TotSumAC, eQuotAC, eQuotK
      real(rkind) StokesNormA, StokesNormB, cPhase
      integer IDsel, ISsel, SelectedK
      logical DoTail
      real(rkind) SumNormStokesA(Nlevel), SumNormStokesB(Nlevel)
      real(rkind) SumZetaCorr, MaxZetaCorr, AvgZetaCorr
      real(rkind) eMinMfact, eMaxMfact, SelectedHS
      real(rkind) MaxStokesNorm, MaxValSinc, StokesNorm, SelectedDEP
      real(rkind) CritError, USTOKES_bar, VSTOKES_bar
      real(rkind) USTOKES_bar_int, VSTOKES_bar_int
      real(rkind) eSum_tot, eSum_tot_int, eWkReal
      real(rkind) eSum_totA, eSum_totA_int
      real(rkind) eSum_totB, eSum_totB_int
      real(rkind) TotalBarotropicErrorUstokes, TotalBarotropicErrorVstokes
      real(rkind) TotalSumUstokes, TotalSumVstokes
      real(rkind) SumHeight
      real(rkind) eJPress_loc, eZetaCorr_loc, eProd, eUint, eVint
      real(rkind) z_r(Nlevel)
      real(rkind) z_w_loc(0:Nlevel), eUSTOKES_loc(Nlevel), eVSTOKES_loc(Nlevel)
      real(rkind) PartialU1(NlevelPartial), PartialV1(NlevelPartial)
      real(rkind) PartialU2(NlevelPartial), PartialV2(NlevelPartial)
      IF (.NOT. L_FIRST_ORDER_ARDHUIN) THEN
        eMinMfact=-3
        eMaxMfact=5
      END IF
      DO IP=1,MNP
        DO k=1,Nlevel
          z_r(k)=(z_w_wav(k,IP)+z_w_wav(k-1,IP))/2
        END DO
        z_w_loc=z_w_wav(:,IP)
        eDep=z_w_loc(Nlevel)-z_w_loc(0)
        IF (L_FIRST_ORDER_ARDHUIN) THEN
          PartialU1(1)=(U_wav(2,IP) - U_wav(1,IP))/(z_r(Nlevel)-z_r(Nlevel-1))
          PartialV1(1)=(V_wav(2,IP) - V_wav(1,IP))/(z_r(Nlevel)-z_r(Nlevel-1))
        ELSE
          DO k=2,Nlevel
            PartialU1(k)=(U_wav(k,IP) - U_wav(k-1,IP))/(z_r(k)-z_r(k-1))
            PartialV1(k)=(V_wav(k,IP) - V_wav(k-1,IP))/(z_r(k)-z_r(k-1))
          END DO
          PartialU1(1)=PartialU1(2)
          PartialV1(1)=PartialV1(2)
          DO k=2,Nlevel-1
            !we compute second differential with three values.
            !We have classic formula
            ! d2f/dx2 = (f(x+h) + f(x-h) -2f(x))/h^2
            ! and this is extended to three arbitrary positions
            ! but only first order accuracy.
            eF1=(z_r(k)-z_r(k-1))/(z_r(k+1)-z_r(k-1))
            eF2=(z_r(k+1)-z_r(k))/(z_r(k+1)-z_r(k-1))
            eDelta=(z_r(k) - z_r(k+1))*(z_r(k-1) - z_r(k))
            PartialU2(k)=(U_wav(k+1,IP)*eF1 + U_wav(k-1,IP)*eF2 - U_wav(k,IP))/eDelta
            PartialV2(k)=(V_wav(k+1,IP)*eF1 + V_wav(k-1,IP)*eF2 - V_wav(k,IP))/eDelta
          END DO
          PartialU2(1)=PartialU2(2)
          PartialV2(1)=PartialV2(2)
          PartialU2(Nlevel)=PartialU2(Nlevel-1)
          PartialV2(Nlevel)=PartialV2(Nlevel-1)
        END IF
        eUSTOKES_loc=0
        eVSTOKES_loc=0
        eJpress_loc=0
        eZetaCorr_loc=0
        IF (L_FIRST_ORDER_ARDHUIN) THEN
          DO IS=1,MSC
            eMult=SPSIG(IS)*DDIR*DS_INCR(IS)
            eWk=WK(IS,IP)
            kD=MIN(KDMAX, eWk*eDep)
            eWkReal=kD/eDep
            eSinh2kd=MySINH(2*kD)
            eSinhkd=MySINH(kD)
            eSinhkd2=eSinhkd**2
            eSigma=SPSIG(IS)
            IF (L_STOKES_DRIFT_USING_INTEGRAL) THEN
              eUint=0
              eVint=0
            END IF
            DO ID=1,MDC
              eLoc=AC2(IS,ID,IP)*eMult
              eScal=COSTH(ID)*PartialU1(1)+SINTH(ID)*PartialV1(1)
              eZeta=eWk/eSinhkd + (eWk/eSigma)*eScal
              eZetaCorr_loc=eZetaCorr_loc + eLoc*eZeta
              eJPress=G9*(kD/eSinh2kd)*(1/eDep) * eLoc
              eJPress_loc=eJPress_loc + eJPress
              IF (L_STOKES_DRIFT_USING_INTEGRAL) THEN
                eUint=eUint + eLoc*COSTH(ID)
                eVint=eVint + eLoc*SINTH(ID)
              END IF
            END DO
            IF (L_STOKES_DRIFT_USING_INTEGRAL) THEN
              DO k=1,Nlevel
                eFrac=(z_r(k) - z_w_loc(0))/eDep
                eHeight=z_w_loc(k)-z_w_loc(k-1)
                eFracB=eHeight/eDep
                eSinc=SINH(kD*eFracB)/(kD*eFracB)
                eQuot1=eSinc*MyCOSH(2*kD*eFrac)/eSinhkd2
                eProd=eSigma*eWkReal*eQuot1
                eUSTOKES_loc(k)=eUSTOKES_loc(k) + eUint*eProd
                eVSTOKES_loc(k)=eVSTOKES_loc(k) + eVint*eProd
              END DO
            END IF
          END DO
        ELSE
          DO IS=1,MSC
            eMult=SPSIG(IS)*DDIR*DS_INCR(IS)
            eWk=WK(IS,IP)
            kD=MIN(KDMAX, eWk*eDep)
            eWkReal=kD/eDep
            eSinh2kd=MySINH(2*kD)
            eSinhkd=MySINH(kD)
            eSinhkd2=eSinhkd**2
            eSigma=SPSIG(IS)
            DO ID=1,MDC
              eLoc=AC2(IS,ID,IP)*eMult
              TheInt=0
              DO k=1,Nlevel
                eHeight=z_w_loc(k)-z_w_loc(k-1)
                zMid=0.5*(z_w_loc(k)+z_w_loc(k-1))
                eFrac=(zMid - z_w_loc(0))/eDep
                eFracB=eHeight/eDep
                eSinc=MySINH(eFracB*kD)/(eFracB*kD)
                eQuot=eWkReal*2*MyCOSH(2*kD*eFrac)/eSinh2kd
                eFct=U_wav(k,IP)*COSTH(ID)+V_wav(k,IP)*SINTH(ID)
                TheInt=TheInt+eHeight*eFct*eQuot*eSinc
              END DO
              eOmega=eSigma + TheInt*eWkReal
              IF (L_STOKES_DRIFT_USING_INTEGRAL) THEN
                DO k=1,Nlevel
                  MFACT=eSigma/(eOmega - (U_wav(k,IP)*COSTH(ID)+V_wav(k,IP)*SINTH(ID))*eWkReal)
                  MFACT=MAX(MFACT, eMinMfact)
                  MFACT=MIN(MFACT, eMaxMfact)
                  eFrac=(z_r(k) - z_w_loc(0))/eDep
                  eHeight=z_w_loc(k)-z_w_loc(k-1)
                  eFracB=eHeight/eDep
                  eSinc=SINH(kD*eFracB)/(kD*eFracB)
                  eQuot1=eSinc*MyCOSH(2*kD*eFrac)/eSinhkd2
                  USTOKES1=MFACT*eSigma*COSTH(ID)*eWkReal*eQuot1
                  VSTOKES1=MFACT*eSigma*SINTH(ID)*eWkReal*eQuot1
                  eQuot2=eSinc*MySINH(2*kD*eFrac)/eSinhkd2
                  eQuot3=eSinc*(MySINH(kD*eFrac)/eSinhkd)**2
                  eScal=PartialU1(k)*COSTH(ID) + PartialV1(k)*SINTH(ID)
                  USTOKES2=0.5*(MFACT**2)*COSTH(ID)*eWkReal*eQuot2*eScal
                  USTOKES3=0.5*MFACT*PartialU2(k)*eQuot3
                  VSTOKES2=0.5*(MFACT**2)*SINTH(ID)*eWkReal*eQuot2*eScal
                  VSTOKES3=0.5*MFACT*PartialV2(k)*eQuot3
                  USTOKESpart=eLoc*(USTOKES1+USTOKES2+USTOKES3)
                  VSTOKESpart=eLoc*(VSTOKES1+VSTOKES2+VSTOKES3)
                  eUSTOKES_loc(k)=eUSTOKES_loc(k) + USTOKESpart
                  eVSTOKES_loc(k)=eVSTOKES_loc(k) + VSTOKESpart
                END DO
              ELSE
                MFACT=eSigma/(eOmega - (U_wav(Nlevel,IP)*COSTH(ID)+V_wav(Nlevel,IP)*SINTH(ID))*eWkReal)
              END IF
              eScal=COSTH(ID)*PartialU1(Nlevel)+SINTH(ID)*PartialV1(Nlevel)
              eZeta=eWk/eSinhkd + (MFACT*eWk/eSigma)*eScal
              eZetaCorr_loc=eZetaCorr_loc + MFACT*eLoc*eZeta
              eJPress=G9*(kD/eSinh2kd)*(1/eDep) * eLoc
              eJPress_loc=eJPress_loc + eJPress
            END DO
          END DO
        END IF
        USTOKES_wav(:,IP)=eUSTOKES_loc
        VSTOKES_wav(:,IP)=eVSTOKES_loc
        ZETA_CORR(IP)=eZetaCorr_loc
        J_PRESSURE(IP)=eJPress_loc
      ENDDO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PGMCL_ROMS_IN(K,IFILE,IT)
      USE pgmcl_library
      USE datapool
      USE mod_coupler
      implicit none
      INTEGER, INTENT(IN) :: K,IFILE,IT
      integer IP, kLev, i, idx
      real(rkind) u1, v1, u2, v2, z1
# ifdef DEBUG_WWM
      real(rkind) :: MaxUwind, SumUwind, avgUwind
      real(rkind) :: MaxVwind, SumVwind, avgVwind
      real(rkind) :: NbPoint
# endif
# ifdef DEBUG_WWM
      WRITE(DBG%FHNDL,*) 'WWM: Begin PGMCL_ROMS_IN'
      FLUSH(DBG%FHNDL)
# endif
      CALL MPI_INTERP_ARECV_3D_r8(TheAsync_OCNtoWAV_uvz, 201, A_wav_uvz)
# ifdef DEBUG_WWM
      WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_IN, After Data receive'
      FLUSH(DBG%FHNDL)
# endif
# ifdef DEBUG_WWM
      MaxUwind=0.0_r8
      SumUwind=0.0_r8
      MaxVwind=0.0_r8
      SumVwind=0.0_r8
      NbPoint=0
# endif
      WATLEVOLD=WATLEV
      DELTAT_WATLEV = MAIN%DTCOUP
      LCALC=.TRUE.
      DO IP=1,MNP
        u1=A_wav_uvz(1,IP)
        v1=A_wav_uvz(2,IP)
# ifdef DEBUG_WWM
        IF (abs(u1).gt.MaxUwind) THEN
          MaxUwind=abs(u1)
        ENDIF
        IF (abs(v1).gt.MaxVwind) THEN
          MaxVwind=abs(v1)
        ENDIF
        SumUwind=SumUwind + abs(u1)
        SumVwind=SumVwind + abs(v1)
        NbPoint=NbPoint+1
# endif
        u2=u1*CosAng(IP)-v1*SinAng(IP)
        v2=v1*CosAng(IP)+u1*SinAng(IP)
        z1=A_wav_uvz(3,IP)
        WINDXY(IP,1)=u2
        WINDXY(IP,2)=v2
        WATLEV(IP)=z1
      END DO
# ifdef DEBUG_WWM
      avgUwind=SumUwind/NbPoint
      avgVwind=SumVwind/NbPoint
      WRITE(DBG%FHNDL,*) 'WAV, MaxUwind=', MaxUwind, ' avgUwind=', avgUwind
      WRITE(DBG%FHNDL,*) 'WAV, MaxVwind=', MaxVwind, ' avgVwind=', avgVwind
      WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_IN, step 2'
      FLUSH(DBG%FHNDL)
# endif
      CALL MPI_INTERP_ARECV_3D_r8(TheAsync_OCNtoWAV_rho, 203, A_wav_rho_3D)
# ifdef DEBUG_WWM
      WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_IN, step 3'
      FLUSH(DBG%FHNDL)
# endif
      DO IP=1,MNP
        DO kLev=0,Nlevel
          z_w_wav(kLev,IP)=A_wav_rho_3D(kLev+1,IP)
        END DO
      END DO
# ifdef DEBUG_WWM
      WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_IN, step 4'
      FLUSH(DBG%FHNDL)
# endif
      CALL MPI_INTERP_ARECV_3D_r8(TheAsync_OCNtoWAV_u, 204, A_wav_ur_3D)
      DO IP=1,MNP
        U_wav(:,IP)=A_wav_ur_3D(:,IP)
      END DO
# ifdef DEBUG_WWM
      WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_IN, step 5'
      FLUSH(DBG%FHNDL)
# endif
      CALL MPI_INTERP_ARECV_3D_r8(TheAsync_OCNtoWAV_v, 205, A_wav_vr_3D)
# ifdef DEBUG_WWM
      WRITE(DBG%FHNDL,*) 'After the receive'
      FLUSH(DBG%FHNDL)
# endif
      DO IP=1,MNP
        V_wav(:,IP)=A_wav_vr_3D(:,IP)
      END DO
# ifdef DEBUG_WWM
      WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_IN, step 6'
      FLUSH(DBG%FHNDL)
# endif
      DO IP=1,MNP
        DO kLev=1,NlevelVert
          u1=U_wav(kLev,IP)
          v1=V_wav(kLev,IP)
          u2=u1*CosAng(IP)-v1*SinAng(IP)
          v2=v1*CosAng(IP)+u1*SinAng(IP)
          U_wav(kLev,IP)=u2
          V_wav(kLev,IP)=v2
        END DO
        CURTXY(IP,1)=u2
        CURTXY(IP,2)=v2
      END DO
# ifdef DEBUG_WWM
      WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_IN, step 7'
      FLUSH(DBG%FHNDL)
# endif
      END SUBROUTINE PGMCL_ROMS_IN
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PGMCL_ROMS_OUT(K)
      USE DATAPOOL
      USE pgmcl_library
      USE mod_coupler
# ifdef ST41
      USE W3SRC4MD, only : UFRIC
# endif
# ifdef ST42
      USE W3SRC4MD, only : UFRIC
# endif
      implicit none
      INTEGER, INTENT(IN)  :: K
      integer IP, kLev, idx
      real(rkind) u1, v1, u2, v2
      real(rkind) HS, TM01, TM02, KLM, WLM, TM10
      real(rkind) UBOT, ORBITAL, BOTEXPER, TMBOT
      real(rkind) FPP, TPP, CPP, WNPP, CGPP, KPP, LPP, PEAKDSPR, PEAKDM, DPEAK
      real(rkind) ETOTS,ETOTC,DM,DSPR
      REAL(RKIND) :: ACLOC(MSC,MDC)
      real(rkind) cPhase, eStokesNorm
      real(rkind) kD
      real(rkind) :: TPPD, KPPD, CGPD, CPPD
# ifdef DEBUG_WWM
      real(rkind) SumNormTau, MaxNormTau, AvgNormTau, eNorm
      real(rkind) AvgUFRICsqr, SumUFRICsqr
      real(rkind) AvgCd, SumCd
      real(rkind) AvgStressCd, SumStressCd, eStressCd, eMag
      real(rkind) AvgAlpha, SumAlpha, eAlpha, NbAlpha
      real(rkind) SumWind, AvgWind
      real(rkind) :: MaxHwave, SumHwave, avgHwave, NbPoint
      real(rkind) :: MaxLwave, SumLwave, avgLwave
      real(rkind) :: MaxTM02, SumTM02, AvgTM02
      real(rkind) :: MaxStokesNorm, SumStokesNorm, avgStokesNorm
      WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_OUT, step 1'
      FLUSH(DBG%FHNDL)
      SumNormTau=0
      MaxNormTau=0
# endif
      CALL STOKES_STRESS_INTEGRAL_ROMS
      DO IP=1,MNP
        IF (L_STOKES_DRIFT_USING_INTEGRAL) THEN
          DO kLev=1,Nlevel
            u1=USTOKES_wav(kLev,IP)
            v1=VSTOKES_wav(kLev,IP)
            u2=u1*CosAng(IP)+v1*SinAng(IP)
            v2=v1*CosAng(IP)-u1*SinAng(IP)
            A_wav_u_3D(kLev,IP)=u2
            A_wav_v_3D(kLev,IP)=v2
          END DO
        END IF
        u1=TAUWX(IP)
        v1=TAUWY(IP)
        u2=u1*CosAng(IP)+v1*SinAng(IP)
        v2=v1*CosAng(IP)-u1*SinAng(IP)
        IF (L_STOKES_DRIFT_USING_INTEGRAL) THEN
          A_wav_u_3D(Nlevel+1,IP)=u2
          A_wav_v_3D(Nlevel+1,IP)=v2
        ELSE
          A_wav_u_3D(1,IP)=u2
          A_wav_v_3D(1,IP)=v2
        END IF
# ifdef DEBUG_WWM
        eNorm=SQRT(u2*u2 + v2*v2)
        IF (eNorm.gt.MaxNormTau) THEN
          MaxNormTau=eNorm
        END IF
        SumNormTau=SumNormTau + eNorm
# endif
      END DO
# ifdef DEBUG_WWM
      AvgNormTau=SumNormTau / MNP
      WRITE(DBG%FHNDL,*) 'AvgNormTau=', AvgNormTau, 'MaxNormTau=', MaxNormTau
      WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_OUT, step 5.1'
      FLUSH(DBG%FHNDL)
# endif
      IF (L_STOKES_DRIFT_USING_INTEGRAL) THEN
        CALL MPI_INTERP_ASEND_3D_r8(TheAsync_WAVtoOCN_u, 208, A_wav_u_3D)
        CALL MPI_INTERP_ASEND_3D_r8(TheAsync_WAVtoOCN_v, 210, A_wav_u_3D)
        CALL MPI_INTERP_ASEND_3D_r8(TheAsync_WAVtoOCN_u, 209, A_wav_v_3D)
        CALL MPI_INTERP_ASEND_3D_r8(TheAsync_WAVtoOCN_v, 211, A_wav_v_3D)
      ELSE
        CALL MPI_INTERP_ASEND_3D_r8(TheAsync_WAVtoOCN_u, 208, A_wav_u_3D)
        CALL MPI_INTERP_ASEND_3D_r8(TheAsync_WAVtoOCN_v, 211, A_wav_v_3D)
      END IF
# ifdef DEBUG_WWM
      WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_OUT, step 5.3'
      FLUSH(DBG%FHNDL)
# endif
# ifdef DEBUG_WWM
      MaxHwave=0.0
      SumHwave=0.0
      MaxTM02=0
      SumTM02=0
      MaxLwave=0.0
      SumLwave=0.0
      SumStokesNorm=0
      MaxStokesNorm=0
      NbPoint=0.0
      SumUFRICsqr=0
      SumCd=0
      SumWind=0
      SumStressCd=0
      SumAlpha=0
      NbAlpha=0
# endif
      DO IP = 1, MNP
        ACLOC = AC2(:,:,IP)
        CALL MEAN_PARAMETER(IP,ACLOC,MSC,HS,TM01,TM02,TM10,KLM,WLM)
        CALL WAVE_CURRENT_PARAMETER(IP,ACLOC,UBOT,ORBITAL,BOTEXPER,TMBOT,'PGMCL_ROMS_OUT')
        CALL MEAN_DIRECTION_AND_SPREAD(IP,ACLOC,MSC,ETOTS,ETOTC,DM,DSPR)
        CALL PEAK_PARAMETER(IP,ACLOC,MSC,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD)
# ifdef DEBUG_WWM
        SumUFRICsqr=SumUFRICsqr + UFRIC(IP)*UFRIC(IP)
        eMag=SQRT(WINDXY(IP,1)**2 + WINDXY(IP,2)**2)
        eStressCd=CD(IP)*eMag*eMag
        SumWind=SumWind + eMag
        SumStressCd=SumStressCd + eStressCd
        IF (UFRIC(IP).gt.0) THEN
          eAlpha=G9*Z0(IP)/(UFRIC(IP) * UFRIC(IP))
          SumAlpha=SumAlpha + eAlpha
          NbAlpha=NbAlpha+1
        END IF
        IF (HS.gt.MaxHwave) THEN
          MaxHwave=HS
        ENDIF
        SumHwave=SumHwave + HS
        IF (TM02.gt.MaxTM02) THEN
          MaxTM02=TM02
        ENDIF
        SumTM02=SumTM02 + TM02
        IF (WLM.gt.MaxLwave) THEN
          MaxLwave=WLM
        ENDIF
        SumLwave=SumLwave + WLM
        kD=MIN(KDMAX, KLM*DEP(IP))
        cPhase=SQRT((G9/KLM)*REAL(MySINH(kD)/MyCOSH(kD)) )
        eStokesNorm=(G9*HS*HS/REAL(16))*2*(KLM/cPhase)
        IF (eStokesNorm.ne.eStokesNorm) THEN
          WRITE(DBG%FHNDL,*) 'eStokesNorm=', eStokesNorm
          WRITE(DBG%FHNDL,*) 'KLM=', KLM, 'WLM=', WLM
          WRITE(DBG%FHNDL,*) 'cPhase=', cPhase, 'kD=', kD
          WRITE(DBG%FHNDL,*) 'HS=', HS, ' DEP=', DEP(IP)
          FLUSH(DBG%FHNDL)
        END IF
        IF (eStokesNorm.gt.MaxStokesNorm) THEN
          MaxStokesNorm=eStokesNorm
        END IF
        SumStokesNorm=SumStokesNorm + eStokesNorm
        NbPoint=NbPoint + 1
# endif
        A_wav_stat(1, IP)=HS
        A_wav_stat(2, IP)=TM01
        A_wav_stat(3, IP)=TM02
        A_wav_stat(4, IP)=KLM
        A_wav_stat(5, IP)=WLM
        A_wav_stat(6, IP)=ORBITAL
        A_wav_stat(7, IP)=TMBOT
        A_wav_stat(8, IP)=DISSIPATION(IP)
        A_wav_stat(9, IP)=QBLOCAL(IP)
        A_wav_stat(10,IP)=DM
        A_wav_stat(11,IP)=TPP
        A_wav_stat(12,IP)=DSPR
        A_wav_stat(13,IP)=PEAKDSPR
        A_wav_stat(14,IP)=PEAKDM
        A_wav_stat(15,IP)=UFRIC(IP)
        A_wav_stat(16,IP)=Z0(IP)
        A_wav_stat(17,IP)=CD(IP)
        A_wav_stat(18,IP)=J_PRESSURE(IP)
        A_wav_stat(19,IP)=ZETA_CORR(IP)
      END DO
# ifdef DEBUG_WWM
      avgHwave=SumHwave/NbPoint
      avgLwave=SumLwave/NbPoint
      AvgTM02=SumTM02/NbPoint
      avgStokesNorm=SumStokesNorm/NbPoint
      WRITE(DBG%FHNDL,*) 'WAV, MaxHwave=', MaxHwave, ' avgHwave=', avgHwave
      WRITE(DBG%FHNDL,*) 'WAV, MaxLwave=', MaxLwave, ' avgLwave=', avgLwave
      WRITE(DBG%FHNDL,*) 'WAV, MaxStokesNorm=', MaxStokesNorm
      WRITE(DBG%FHNDL,*) 'WAV, avgStokesNorm=', avgStokesNorm
      WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_OUT, step 6'
      FLUSH(DBG%FHNDL)
      AvgUFRICsqr=SumUFRICsqr/MNP
      AvgStressCd=SumStressCd/MNP
      AvgAlpha=SumAlpha/NbAlpha
      AvgWind=SumWind/MNP
      AvgCd=SumCd/MNP
      WRITE(DBG%FHNDL,*) 'AvgNormTau=', AvgNormTau, 'AvgNormFV2=', AvgUFRICsqr
      WRITE(DBG%FHNDL,*) 'AvgNormTau=', AvgNormTau, 'AvgCdU2=', AvgStressCd
      WRITE(DBG%FHNDL,*) 'AvgCd=', AvgCd, ' AvgAlpha=', AvgAlpha
      WRITE(DBG%FHNDL,*) 'AvgWind=', AvgWind
      FLUSH(DBG%FHNDL)
# endif
      CALL MPI_INTERP_ASEND_3D_r8(TheAsync_WAVtoOCN_stat, 212, A_wav_stat)
# ifdef DEBUG_WWM
      WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_OUT, step 11'
      FLUSH(DBG%FHNDL)
# endif
      END SUBROUTINE
#endif
END MODULE
