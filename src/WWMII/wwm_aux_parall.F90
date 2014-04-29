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
        DO iProc=2,nproc
          MNPloc=ListMNP(iProc)
          WRITE(STAT%FHNDL,*) 'iProc, MNPloc=', iProc, MNPloc
          allocate(dspl_send(MNPloc))
          DO IP=1,MNPloc
            IP_glob=ListIPLG(IP+ListFirst(iProc))
            dspl_send(IP)=IP_glob
            WRITE(STAT%FHNDL,*) 'IP,IP_glob=', IP, IP_glob
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
      Print *, 'Bonjour'
      DO IP=1,np_total
        Vtotal(IP)=DBLE(IP)
      END DO
      Print *, 'After Vtotal assignation'
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
      WRITE(STAT%FHNDL,*) 'MNP=', MNP
      DO IP=1,MNP
        IF (ABS(Vlocal(IP) - iplg(IP)) > 1) THEN
          Print *, 'BUG IP,Vloc,iplg=', IP, Vlocal(IP), iplg(IP)
          WRITE(STAT%FHNDL,*) 'BUG IP,Vloc,iplg=', IP, Vlocal(IP), iplg(IP)
        END IF
        WRITE(STAT%FHNDL,*) 'IP,Vloc,iplg=', IP, Vlocal(IP), iplg(IP)
      END DO
      FLUSH(STAT%FHNDL)
      Vlocal = 0
      Print *, 'Au revoir'
      END SUBROUTINE
#endif
