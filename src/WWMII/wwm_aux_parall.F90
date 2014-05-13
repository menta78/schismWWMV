#include "wwm_functions.h"
#ifdef MPI_PARALL_GRID
!**********************************************************************
!* Some MPI_BARRIER are just not reliable. This construction makes    *
!* that every process receive and send to every other process         *
!**********************************************************************
      SUBROUTINE MYOWN_MPI_BARRIER(iresult)
      USE DATAPOOL
      IMPLICIT NONE
      integer, intent(in) :: iresult
      integer eInt(1)
      integer iRank, jRank, eTag
      eInt(1)=4
      WRITE(STAT%FHNDL,*) 'Before the loop of send/recv stat=', iresult
      FLUSH(STAT%FHNDL)
      DO iRank=0,nproc-1
        eTag=137 + iRank
        IF (myrank .eq. iRank) THEN
          DO jRank=0,nproc-1
            IF (iRank.ne. jRank) THEN
              CALL MPI_SEND(eInt,1,itype,jRank,eTag,comm,ierr)
            END IF
          END DO
        ELSE
          CALL MPI_RECV(eInt, 1, itype,iRank,eTag,comm,istatus,ierr)
        END IF
      END DO
      WRITE(STAT%FHNDL,*) 'After the loop of send/recv stat=', iresult
      FLUSH(STAT%FHNDL)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE AC_COHERENCY(AC, string)
      USE DATAPOOL
      IMPLICIT NONE
      character(*), intent(in) :: string
      REAL(rkind), intent(in) :: AC(MNP,MSC,MDC)
      REAL(rkind) :: ACwork(MSC,MDC,MNP)
      REAL(rkind) :: Lerror
      INTEGER IP
      DO IP=1,MNP
        ACwork(:,:,IP)=AC(IP,:,:)
      END DO
      CALL I5B_TOTAL_COHERENCY_ERROR(MSC, ACwork, Lerror)
      WRITE(STAT%FHNDL,*) 'coherency error between domains'
      WRITE(STAT%FHNDL,*) 'Lerror=', Lerror, ' mesg=', TRIM(string)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_TOTAL_COHERENCY_ERROR(MSCeffect, ACw, Lerror)
      USE DATAPOOL, only : MNP, MDC, rkind
      USE DATAPOOL, only : ListIPLG, ListMNP
      USE datapool, only : istatus, ierr, comm, rtype, myrank, nproc, iplg, np_global
      implicit none
      integer, intent(in) :: MSCeffect
      real(rkind), intent(in) :: ACw(MSCeffect, MDC, MNP)
      real(rkind), intent(out) :: Lerror
      real(rkind), allocatable :: ACtotal(:,:,:)
      real(rkind), allocatable :: ACloc(:,:,:)
      real(rkind) :: rbuf_real(1)
      integer, allocatable :: ListFirstMNP(:)
      integer, allocatable :: eStatus(:)
      integer IP, iProc, IPglob, IS, ID
      integer MNPloc
      integer istat
      IF (myrank == 0) THEN
        Lerror=0
        allocate(ListFirstMNP(nproc), eStatus(np_global), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 69')
        ListFirstMNP=0
        eStatus=0
        DO iProc=2,nproc
          ListFirstMNP(iProc)=ListFirstMNP(iProc-1) + ListMNP(iProc-1)
        END DO
        allocate(ACtotal(MSCeffect, MDC, np_global), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 70')
        DO IP=1,MNP
          IPglob=iplg(IP)
          ACtotal(:,:,IPglob)=ACw(:,:,IP)
          eStatus(IPglob)=1
        END DO
        DO iProc=2,nproc
          MNPloc=ListMNP(iProc)
          allocate(ACloc(MSCeffect, MDC, MNPloc), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 71')
          CALL MPI_RECV(ACloc,MNPloc*MSCeffect*MDC,rtype, iProc-1, 53, comm, istatus, ierr)
          DO IP=1,MNPloc
            IPglob=ListIPLG(IP+ListFirstMNP(iProc))
            IF (eStatus(IPglob) == 1) THEN
              DO IS=1,MSCeffect
                DO ID=1,MDC
                  Lerror=Lerror+abs(ACtotal(IS,ID,IPglob)-ACloc(IS,ID,IP))
                END DO
              END DO
            ELSE
              eStatus(IPglob)=1
              ACtotal(:,:,IPglob)=ACloc(:,:,IP)
            END IF
          END DO
          deallocate(ACloc)
        END DO
        deallocate(ListFirstMNP, ACtotal, eStatus)
        rbuf_real(1)=Lerror
        DO iProc=2,nproc
          CALL MPI_SEND(rbuf_real,1,rtype, iProc-1, 23, comm, ierr)
        END DO
      ELSE
        CALL MPI_SEND(ACw,MNP*MSCeffect*MDC,rtype, 0, 53, comm, ierr)
        CALL MPI_RECV(rbuf_real,1,rtype, 0, 23, comm, istatus, ierr)
        Lerror=rbuf_real(1)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_TOTAL_COHERENCY_ERROR_NPRES(MSCeffect, ACw, Lerror)
      USE DATAPOOL, only : MNP, MDC, NP_RES, rkind
      USE DATAPOOL, only : ListIPLG, ListMNP, ListNP_RES
      USE datapool, only : istatus, ierr, comm, rtype, myrank, nproc, iplg, np_global
      implicit none
      integer, intent(in) :: MSCeffect
      real(rkind), intent(in) :: ACw(MSCeffect, MDC, MNP)
      real(rkind), intent(out) :: Lerror
      real(rkind), allocatable :: ACtotal(:,:,:)
      real(rkind), allocatable :: ACloc(:,:,:)
      real(rkind) :: rbuf_real(1)
      integer, allocatable :: ListFirstMNP(:)
      integer, allocatable :: eStatus(:)
      integer IP, iProc, IPglob, IS, ID
      integer NP_RESloc
      integer istat
      IF (myrank == 0) THEN
        Lerror=0
        allocate(ListFirstMNP(nproc), eStatus(np_global), ACtotal(MSCeffect, MDC, np_global), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 72')
        ListFirstMNP=0
        eStatus=0
        DO iProc=2,nproc
          ListFirstMNP(iProc)=ListFirstMNP(iProc-1) + ListMNP(iProc-1)
        END DO
        DO IP=1,NP_RES
          IPglob=iplg(IP)
          ACtotal(:,:,IPglob)=ACw(:,:,IP)
          eStatus(IPglob)=1
        END DO
        DO iProc=2,nproc
          NP_RESloc=ListNP_RES(iProc)
          allocate(ACloc(MSCeffect, MDC, NP_RESloc), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 73')
          CALL MPI_RECV(ACloc,MSCeffect*MDC*NP_RESloc,rtype, iProc-1, 53, comm, istatus, ierr)
          DO IP=1,NP_RESloc
            IPglob=ListIPLG(IP+ListFirstMNP(iProc))
            IF (eStatus(IPglob) == 1) THEN
              DO IS=1,MSCeffect
                DO ID=1,MDC
                  Lerror=Lerror+abs(ACtotal(IS,ID,IPglob)-ACloc(IS,ID,IP))
                END DO
              END DO
            ELSE
              eStatus(IPglob)=1
              ACtotal(:,:,IPglob)=ACloc(:,:,IP)
            END IF
          END DO
          deallocate(ACloc)
        END DO
        deallocate(ListFirstMNP, ACtotal, eStatus)
        rbuf_real(1)=Lerror
        DO iProc=2,nproc
          CALL MPI_SEND(rbuf_real,1,rtype, iProc-1, 23, comm, ierr)
        END DO
      ELSE
        allocate(ACloc(MSCeffect, MDC, NP_RES), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 74')
        DO IP=1,NP_RES
          ACloc(:,:,IP)=ACw(:,:,IP)
        END DO
        CALL MPI_SEND(ACloc,NP_RES*MSCeffect*MDC,rtype, 0, 53, comm, ierr)
        deallocate(ACloc)
        CALL MPI_RECV(rbuf_real,1,rtype, 0, 23, comm, istatus, ierr)
        Lerror=rbuf_real(1)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COLLECT_ALL_IPLG
      USE DATAPOOL
      implicit none
      integer, allocatable :: rbuf_int(:)
      integer len, iProc, IP, idx, sumMNP
      allocate(ListMNP(nproc), ListNP_RES(nproc), rbuf_int(2), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 15')
      IF (myrank == 0) THEN
        ListMNP(1)=MNP
        ListNP_RES(1)=NP_RES
        DO iProc=2,nproc
          CALL MPI_RECV(rbuf_int,2,itype, iProc-1, 257, comm, istatus, ierr)
          ListMNP(iProc)=rbuf_int(1)
          ListNP_RES(iProc)=rbuf_int(2)
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(ListMNP,nproc,itype, iProc-1, 263, comm, ierr)
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(ListNP_RES,nproc,itype, iProc-1, 571, comm, ierr)
        END DO
      ELSE
        rbuf_int(1)=MNP
        rbuf_int(2)=NP_RES
        CALL MPI_SEND(rbuf_int,2,itype, 0, 257, comm, ierr)
        CALL MPI_RECV(ListMNP,nproc,itype, 0, 263, comm, istatus, ierr)
        CALL MPI_RECV(ListNP_RES,nproc,itype, 0, 571, comm, istatus, ierr)
      END IF
      deallocate(rbuf_int)
      sumMNP=sum(ListMNP)
      allocate(ListIPLG(sumMNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 16')
      IF (myrank == 0) THEN
        idx=0
        DO IP=1,MNP
          idx=idx+1
          ListIPLG(idx)=iplg(IP)
        END DO
        DO iProc=2,nproc
          len=ListMNP(iProc)
          allocate(rbuf_int(len), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 17')
          CALL MPI_RECV(rbuf_int,len,itype, iProc-1, 269, comm, istatus, ierr)
          DO IP=1,len
            idx=idx+1
            ListIPLG(idx)=rbuf_int(IP)
          END DO
          deallocate(rbuf_int)
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(ListIPLG,sumMNP,itype, iProc-1, 271, comm, ierr)
        END DO
      ELSE
        CALL MPI_SEND(iplg,MNP,itype, 0, 269, comm, ierr)
        CALL MPI_RECV(ListIPLG,sumMNP,itype, 0, 271, comm, istatus, ierr)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COLLECT_ALL_IA_JA
      USE DATAPOOL
      implicit none
      integer, allocatable :: rbuf_int(:)
      integer len, iProc, IP, idx
      integer sumIAsiz, sumNNZ
      allocate(ListNNZ(nproc), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 20')
      !
      ! Collecting NNZ
      !
      IF (myrank == 0) THEN
        ListNNZ(1)=NNZ
        allocate(rbuf_int(1), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 21')
        DO iProc=2,nproc
          CALL MPI_RECV(rbuf_int,1,itype, iProc-1, 257, comm, istatus, ierr)
          ListNNZ(iProc)=rbuf_int(1)
        END DO
        deallocate(rbuf_int)
        DO iProc=2,nproc
          CALL MPI_SEND(ListNNZ,nproc,itype, iProc-1, 263, comm, ierr)
        END DO
      ELSE
        allocate(rbuf_int(1), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 22')
        rbuf_int(1)=NNZ
        CALL MPI_SEND(rbuf_int,1,itype, 0, 257, comm, ierr)
        deallocate(rbuf_int)
        CALL MPI_RECV(ListNNZ,nproc,itype, 0, 263, comm, istatus, ierr)
      END IF
      !
      ! Collecting IA
      !
      sumIAsiz=sum(ListMNP) + nproc
      allocate(ListIA(sumIAsiz), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 23')
      IF (myrank == 0) THEN
        idx=0
        DO IP=1,MNP+1
          idx=idx+1
          ListIA(idx)=IA(IP)
        END DO
        DO iProc=2,nproc
          len=ListMNP(iProc)+1
          allocate(rbuf_int(len), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 24')
          CALL MPI_RECV(rbuf_int,len,itype, iProc-1, 269, comm, istatus, ierr)
          DO IP=1,len
            idx=idx+1
            ListIA(idx)=rbuf_int(IP)
          END DO
          deallocate(rbuf_int)
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(ListIA,sumIAsiz,itype, iProc-1, 271, comm, ierr)
        END DO
      ELSE
        CALL MPI_SEND(IA,MNP+1,itype, 0, 269, comm, ierr)
        CALL MPI_RECV(ListIA,sumIAsiz,itype, 0, 271, comm, istatus, ierr)
      END IF
      !
      ! Collecting JA
      !
      sumNNZ=sum(ListNNZ)
      allocate(ListJA(sumNNZ), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 25')
      IF (myrank == 0) THEN
        idx=0
        DO IP=1,NNZ
          idx=idx+1
          ListJA(idx)=JA(IP)
        END DO
        DO iProc=2,nproc
          len=ListNNZ(iProc)
          allocate(rbuf_int(len), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 26')
          CALL MPI_RECV(rbuf_int,len,itype, iProc-1, 569, comm, istatus, ierr)
          DO IP=1,len
            idx=idx+1
            ListJA(idx)=rbuf_int(IP)
          END DO
          deallocate(rbuf_int)
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(ListJA,sumNNZ,itype, iProc-1, 467, comm, ierr)
        END DO
      ELSE
        CALL MPI_SEND(JA,NNZ,itype, 0, 569, comm, ierr)
        CALL MPI_RECV(ListJA,sumNNZ,itype, 0, 467, comm, istatus, ierr)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SYMM_GRAPH_BUILD_ADJ(AdjGraph)
      USE DATAPOOL, only : wwm_nnbr, wwm_ListNeigh, myrank, Graph
      implicit none
      type(Graph), intent(inout) :: AdjGraph
      CALL KERNEL_GRAPH_BUILD_ADJ(AdjGraph, wwm_nnbr, wwm_ListNeigh)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GRAPH_BUILD_ADJ(AdjGraph)
      USE datapool, only : nnbr_p, nbrrank_p, Graph
      implicit none
      type(Graph), intent(inout) :: AdjGraph
      integer :: ListNe(nnbr_p)
      integer I
      DO I=1,nnbr_p
        ListNe(I)=nbrrank_p(I)+1
      END DO
      CALL KERNEL_GRAPH_BUILD_ADJ(AdjGraph, nnbr_p, ListNe)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE KERNEL_GRAPH_BUILD_ADJ(AdjGraph, nb, ListNe)
      USE DATAPOOL
      implicit none
      integer, intent(in) :: nb
      integer, intent(in) :: ListNe(nb)
      type(Graph), intent(inout) :: AdjGraph
      integer, allocatable :: rbuf_int(:)
      integer I, iProc
      integer idx, eDeg, nbEdge, iEdge
      AdjGraph % nbVert=nproc
      IF (myrank.eq.0) THEN
        allocate(AdjGraph % ListDegree(nproc), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 1')
        AdjGraph % ListDegree(1)=nb
        allocate(rbuf_int(1), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 2')
        DO iProc=2,nproc
          CALL MPI_RECV(rbuf_int,1,itype, iProc-1, 19, comm, istatus, ierr)
          AdjGraph % ListDegree(iProc)=rbuf_int(1)
        END DO
        deallocate(rbuf_int)
        nbEdge=0
        DO iProc=1,nproc
          nbEdge=nbEdge + AdjGraph % ListDegree(iProc)
        END DO
        AdjGraph % nbEdge=nbEdge
        allocate(AdjGraph % ListEdge(nbEdge,2), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 3')
        idx=0
        eDeg=AdjGraph % ListDegree(1)
        DO I=1,eDeg
          idx=idx+1
          AdjGraph % ListEdge(idx,1)=1
          AdjGraph % ListEdge(idx,2)=ListNe(I)
        END DO
        DO iProc=2,nproc
          eDeg=AdjGraph % ListDegree(iProc)
          allocate(rbuf_int(eDeg), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 4')
          CALL MPI_RECV(rbuf_int,eDeg,itype, iProc-1, 24, comm, istatus, ierr)
          DO I=1,eDeg
            idx=idx+1
            AdjGraph % ListEdge(idx,1)=iProc
            AdjGraph % ListEdge(idx,2)=rbuf_int(I)
          END DO
          deallocate(rbuf_int)
        END DO
        allocate(rbuf_int(1), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 5')
        rbuf_int(1)=nbEdge
        DO iProc=2,nproc
          CALL MPI_SEND(rbuf_int,1,itype, iProc-1, 30, comm, ierr)
        END DO
        deallocate(rbuf_int)
        !
        allocate(rbuf_int(nproc), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 6')
        DO iProc=1,nproc
          rbuf_int(iProc)=AdjGraph % ListDegree(iProc)
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(rbuf_int,nproc,itype, iProc-1, 32, comm, ierr)
        END DO
        deallocate(rbuf_int)
        !
        allocate(rbuf_int(nbEdge*2), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 7')
        DO iEdge=1,nbEdge
          rbuf_int(2*(iEdge-1)+1)=AdjGraph % ListEdge(iEdge,1)
          rbuf_int(2*(iEdge-1)+2)=AdjGraph % ListEdge(iEdge,2)
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(rbuf_int,2*nbEdge,itype, iProc-1, 34, comm, ierr)
        END DO
      ELSE
        allocate(rbuf_int(1), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 8')
        rbuf_int(1)=nb
        CALL MPI_SEND(rbuf_int,1,itype, 0, 19, comm, ierr)
        deallocate(rbuf_int)
        !
        allocate(rbuf_int(nb), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 9')
        DO I=1,nb
          rbuf_int(I)=ListNe(I)
        END DO
        CALL MPI_SEND(rbuf_int,nb,itype, 0, 24, comm, ierr)
        deallocate(rbuf_int)
        !
        allocate(rbuf_int(1), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 10')
        CALL MPI_RECV(rbuf_int,1,itype, 0, 30, comm, istatus, ierr)
        nbEdge=rbuf_int(1)
        deallocate(rbuf_int)
        AdjGraph % nbEdge=nbEdge
        !
        allocate(rbuf_int(nproc), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 11')
        CALL MPI_RECV(rbuf_int,nproc,itype, 0, 32, comm, istatus, ierr)
        allocate(AdjGraph % ListDegree(nproc), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 12')
        AdjGraph % ListDegree=rbuf_int
        deallocate(rbuf_int)
        !
        allocate(rbuf_int(2*nbEdge), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 13')
        CALL MPI_RECV(rbuf_int,2*nbEdge,itype, 0, 34, comm, istatus, ierr)
        allocate(AdjGraph % ListEdge(nbEdge,2), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 14')
        DO iEdge=1,nbEdge
          AdjGraph % ListEdge(iEdge,1)=rbuf_int(2*(iEdge-1)+1)
          AdjGraph % ListEdge(iEdge,2)=rbuf_int(2*(iEdge-1)+2)
        END DO
        deallocate(rbuf_int)
      ENDIF
      AdjGraph % MaxDeg=maxval(AdjGraph % ListDegree)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GRAPH_TEST_CONNECT(AdjGraph, result)
      USE DATAPOOL, only : Graph
      implicit none
      type(Graph), intent(in) :: AdjGraph
      integer, intent(out) :: result
      integer :: ListStatus(AdjGraph%nbVert)
      integer :: ListPosFirst(AdjGraph%nbVert)
      integer idx, iVert, nbVert, eAdj
      integer eDeg, sizConn, I, IsFinished
      idx=0
      nbVert=AdjGraph%nbVert
      ListStatus=0
      DO iVert=1,nbVert
        ListPosFirst(iVert)=idx
        idx=idx+AdjGraph % ListDegree(iVert)
      END DO
      ListStatus(1)=1
      DO
        IsFinished=1
        DO iVert=1,nbVert
          IF (ListStatus(iVert) == 1) THEN
            eDeg=AdjGraph%ListDegree(iVert)
            idx=ListPosFirst(iVert)
            DO I=1,eDeg
              eAdj=AdjGraph%ListEdge(idx+I,1)
              IF (ListStatus(eAdj) == 0) THEN
                IsFinished=0
              END IF
              ListStatus(eAdj)=1
            END DO
          END IF
        END DO
        IF (IsFinished == 1) THEN
          EXIT
        END IF
      END DO
      sizConn=0
      DO iVert=1,nbVert
        IF (ListStatus(iVert) == 1) THEN
          sizConn=sizConn+1
        END IF
      END DO
      IF (sizConn == nbVert) THEN
        result=1
      END IF
      result=0
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SETUP_ONED_SCATTER_ARRAY
      USE DATAPOOL
      IMPLICIT NONE
      integer :: ListFirst(nproc)
      integer MNPloc, iProc, IP, IP_glob
      integer, allocatable :: dspl_send(:)
      ListFirst=0
      DO iProc=2,nproc
        ListFirst(iProc)=ListFirst(iProc-1) + ListMNP(iProc-1)
      END DO
      IF (myrank .eq. 0) THEN
        allocate(oned_send_rqst(nproc-1), oned_send_stat(MPI_STATUS_SIZE,nproc-1), oned_send_type(nproc-1), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('error in IOBPtotal allocate')
        DO iProc=2,nproc
          MNPloc=ListMNP(iProc)
          WRITE(STAT%FHNDL,*) 'iProc, MNPloc=', iProc, MNPloc
          allocate(dspl_send(MNPloc), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('error in IOBPtotal allocate')
          DO IP=1,MNPloc
            IP_glob=ListIPLG(IP+ListFirst(iProc))
            dspl_send(IP)=IP_glob-1
          END DO
          call mpi_type_create_indexed_block(MNPloc,1,dspl_send,rtype,oned_send_type(iProc-1), ierr)
          call mpi_type_commit(oned_send_type(iProc-1), ierr)
          deallocate(dspl_send)
        END DO
        FLUSH(STAT%FHNDL)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SCATTER_ONED_ARRAY(Vtotal, Vlocal)
      USE DATAPOOL
      IMPLICIT NONE
      real(rkind) :: Vtotal(np_total)
      real(rkind) :: Vlocal(MNP)
      integer iProc, IP
      IF (myrank .eq. 0) THEN
        DO iProc=2,nproc
          CALL mpi_isend(Vtotal, 1, oned_send_type(iProc-1), iProc-1, 2030, comm, oned_send_rqst(iProc-1), ierr)
        END DO
        DO IP=1,MNP
          Vlocal(IP)=Vtotal(iplg(IP))
        END DO
        IF (nproc > 1) THEN
          CALL MPI_WAITALL(nproc-1, oned_send_rqst, oned_send_stat, ierr)
        END IF
      ELSE
        CALL MPI_RECV(Vlocal, MNP, rtype, 0, 2030, comm, istatus, ierr)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SETUP_BOUNDARY_SCATTER_ARRAY
      USE DATAPOOL
      IMPLICIT NONE
      integer :: ListFirst(nproc)
      integer MNPloc, iProc, IP, IP_glob
      integer, allocatable :: dspl_send(:), Indexes(:)
      integer :: NbSend(nproc)
      integer irank, eSend, idx, idx_nbproc
      ListFirst=0
      DO iProc=2,nproc
        ListFirst(iProc)=ListFirst(iProc-1) + ListMNP(iProc-1)
      END DO
      IF (myrank .eq. rank_boundary) THEN
        spparm_nbproc=0
        DO irank=0,nproc-1
          eSend=0
          IF (irank .ne. rank_boundary) THEN
            iProc=irank+1
            MNPloc=ListMNP(iProc)
            DO IP=1,MNPloc
              IP_glob=ListIPLG(IP+ListFirst(iProc))
              IF ((IOBPtotal(IP_glob) .eq. 2).or.(IOBPtotal(IP_glob) .eq. 4)) THEN
                eSend=eSend+1
              END IF
            END DO
            IF (eSend .gt. 0) THEN
              spparm_nbproc=spparm_nbproc+1
            END IF
          END IF
          NbSend(iProc)=eSend
        END DO
        allocate(spparm_listproc(spparm_nbproc), spparm_send_rqst(spparm_nbproc), spparm_send_stat(MPI_STATUS_SIZE,spparm_nbproc), spparm_send_type(spparm_nbproc), Indexes(np_total), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('error in IOBPtotal allocate')
        idx=0
        DO IP=1,np_total
          IF ((IOBPtotal(IP_glob) .eq. 2).or.(IOBPtotal(IP_glob) .eq. 4)) THEN
            idx=0
            Indexes(IP)=idx
          END IF
        END DO
        IF (IWBMNP .gt. 0) THEN
          allocate(Indexes_boundary(IWBMNP), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('error in IOBPtotal allocate')
          DO IP=1,IWBMNP
            IP_glob=iplg(IWBNDLC(IP))
            Indexes_boundary(IP)=Indexes(IP_glob)
          END DO
        END IF
        idx_nbproc=0
        DO irank=0,nproc-1
          eSend=NbSend(iProc)
          IF ((irank .ne. rank_boundary).and.(eSend.gt.0)) THEN
            iProc=irank+1
            idx_nbproc=idx_nbproc+1
            spparm_listproc(idx_nbproc)=iProc
            MNPloc=ListMNP(iProc)
            allocate(dspl_send(eSend), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('error in IOBPtotal allocate')
            idx=0
            DO IP=1,MNPloc
              IP_glob=ListIPLG(IP+ListFirst(iProc))
              IF ((IOBPtotal(IP_glob) .eq. 2).or.(IOBPtotal(IP_glob) .eq. 4)) THEN
                idx=idx+1
                dspl_send(idx)=Indexes(IP_glob)-1
              END IF
            END DO
            call mpi_type_create_indexed_block(eSend,8,dspl_send,rtype,spparm_send_type(idx_nbproc), ierr)
            call mpi_type_commit(spparm_send_type(idx_nbproc), ierr)
            deallocate(dspl_send)
          END IF
        END DO
        deallocate(Indexes)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SCATTER_BOUNDARY_ARRAY(Vtotal, Vlocal)
      USE DATAPOOL
      IMPLICIT NONE
      real(rkind) :: Vtotal(IWBMNPGL)
      real(rkind) :: Vlocal(IWBMNP)
      integer iProc, IP, irank
      IF ((IWBMNP .eq. 0).and.(myrank.ne.rank_boundary)) THEN
        RETURN
      END IF
      IF (myrank .eq. rank_boundary) THEN
        DO irank=1,spparm_nbproc
          CALL mpi_isend(Vtotal, 1, spparm_send_type(irank), spparm_listproc(irank)-1, 2030, comm, spparm_send_rqst(irank), ierr)
        END DO
        DO IP=1,IWBMNP
          Vlocal(IP)=Vtotal(Indexes_boundary(IP))
        END DO
        IF (spparm_nbproc > 0) THEN
          CALL MPI_WAITALL(spparm_nbproc, spparm_send_rqst, spparm_send_stat, ierr)
        END IF
      ELSE
        CALL MPI_RECV(Vlocal, IWBMNP, rtype, 0, 2030, comm, istatus, ierr)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
# if defined NETCDF && defined DEBUG
      SUBROUTINE NETCDF_WRITE_MATRIX(LocalColor, ASPAR)
      USE DATAPOOL
      USE NETCDF
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      integer, SAVE :: iSystem = 1
      integer MSCeffect
      integer ired, ncid, var_id
      MSCeffect=LocalColor%MSCeffect
      WRITE (FILE_NAME,10) TRIM(PRE_FILE_NAME),nproc, iSystem, myrank
  10  FORMAT (a,'_np',i2.2,'_syst',i3.3,'_iproc',i4.4, '.nc')
      iret = nf90_create(TRIM(FILE_NAME), NF90_CLOBBER, ncid)
      iret = nf90_def_dim(ncid, 'iter', NF90_UNLIMITED, iter_dims)
      iret = nf90_def_dim(ncid, 'three', 3, three_dims)
      iret = nf90_def_dim(ncid, 'msc', MSCeffect, msc_dims)
      iret = nf90_def_dim(ncid, 'mdc', MDC, mdc_dims)
      iret = nf90_def_dim(ncid, 'bbz', MNP, mnp_dims)
      iret = nf90_def_dim(ncid, 'mnpp', MNP+1, mnpp_dims)
      iret = nf90_def_dim(ncid, 'np_global', np_global, npgl_dims)
      iret = nf90_def_dim(ncid, 'np_res', NP_RES, np_res_dims)
      iret = nf90_def_dim(ncid, 'mne', MNE, mne_dims)
      iret = nf90_def_dim(ncid, 'nnz', NNZ, nnz_dims)
      iret = nf90_def_var(ncid,'ASPAR',NF90_DOUBLE,(/msc_dims, mdc_dims,nnz_dims/),var_id)
      iret = nf90_def_var(ncid,'IA',NF90_INT,(/mnpp_dims/),var_id)
      iret = nf90_def_var(ncid,'JA',NF90_INT,(/nnz_dims/),var_id)
      iret = nf90_def_var(ncid,'iplg',NF90_INT,(/ mnp_dims/),var_id)
      iret = nf90_close(ncid)
      !
      iret = nf90_open(TRIM(FILE_NAME), NF90_WRITE, ncid)
      iret=nf90_inq_varid(ncid, 'iplg', var_id)
      iret=nf90_put_var(ncid,var_id,iplg,start=(/1/), count = (/ MNP /))
      iret=nf90_inq_varid(ncid, 'IA', var_id)
      iret=nf90_put_var(ncid,var_id,IA,start=(/1/), count = (/ MNP+1 /))
      !
      iret=nf90_inq_varid(ncid, 'JA', var_id)
      iret=nf90_put_var(ncid,var_id,JA,start=(/1/), count = (/ NNZ /))
      !
      iret=nf90_inq_varid(ncid, 'ine', var_id)
      iret=nf90_put_var(ncid,var_id,INE,start=(/1,1/), count = (/ 3, MNE /))
      !
      iret=nf90_inq_varid(ncid, 'ASPAR', var_id)
      iret=nf90_put_var(ncid,var_id,ASPAR,start=(/1,1,1/), count = (/ MSC, MDC, NNZ/))
      iret=nf90_inq_varid(ncid, 'IA', var_id)
      iret=nf90_put_var(ncid,var_id,IA,start=(/1/), count = (/ MNP+1/))
      iret=nf90_inq_varid(ncid, 'JA', var_id)
      iret=nf90_put_var(ncid,var_id,JA,start=(/1/), count = (/ NNZ/))
      iret=nf90_inq_varid(ncid, 'ListPos', var_id)
      iret=nf90_put_var(ncid,var_id,ListPos,start=(/1/), count = (/ np_global/))
      iret = nf90_close(ncid)
      !
      END SUBROUTINE
# endif
#endif
