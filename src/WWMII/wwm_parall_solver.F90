#include "wwm_functions.h"
! I5 is under the assumptions that we have a lot of memory
!    and tries to minimize the number of MPI exchanges and
!    to have the processors be as busy as possible
!    We use memory ordered as AC(MNP,MSC,MDC)
! I5B is the same as I5. We use memory ordered as AC(MSC,MDC,MNP)
!    so reordering at the beginning but less operations later on.
! I4 is like I5 but we split the 1,MSC into Nblocks
!    so, there are actually Nblocks times more exchanges.
#define DEBUG
#undef DEBUG

#define PLAN_I4
#undef PLAN_I4

! This is for the reordering of ASPAR_pc and hopefully higher speed
! in the application of the preconditioner.
#undef REORDER_ASPAR_PC
#define REORDER_ASPAR_PC
! This is for the computation of ASPAR_block by a block algorithm
! with hopefully higher speed.
#undef ASPAR_B_COMPUTE_BLOCK
#define ASPAR_B_COMPUTE_BLOCK
! An algorithm that should be slightly faster for norm computations
#undef FAST_NORM
#define FAST_NORM
! Either we use the SELFE exchange routine or ours that exchanges only
! the ghost nodes and not the interface nodes.
#undef SELFE_EXCH
#define SELFE_EXCH
! Repeated CX/CY computations but less memory used.
#undef NO_MEMORY_CX_CY
#define NO_MEMORY_CX_CY
! New less memory intensive method obtained by rewriting the BCGS
! Needs 7 times MSC*MDC*MNP versus 9 times.
#undef BCGS_REORG
#define BCGS_REORG
! Rewriting of some exchange routines for LU solve
#undef LU_SOLVE_RWRT
#define LU_SOLVE_RWRT
! Operation L2U which are smaller than BLK_
#undef L2U_OPER
#define L2U_OPER
!**********************************************************************
!* We have to think on how the system is solved. Many questions are   *
!* mixed: the ordering of the nodes, the ghost nodes, the aspar array *
!* Here is a repository of the conclusions that have been reached     *
!*                                                                    *
!* Ordering 1> should be that way: We have two global nodes i and j.  *
!* -- If i and j belong to a common local grid, then we select the    *
!*    grid of lowest color and decide whether ipgl(i) < ipgl(j)       *
!* -- If i and j belong to two different grid then                    *
!*     ---If Color(i) < Color(j) or reverse we decide by that         *
!*     ---If Color(i) = Color(j) we decide by i<j or not (but it      *
!*        does not matter to the solution)                            *
!* The functions WRITE_EXPLICIT_ORDERING does exactly that and        *
!* provides an ordering that can be used. That is we start with the   *
!* nodes of lowest color and index until we arrive at highest color   *
!*                                                                    *
!* The ASPAR is computed correctly only on 1:NP_RES but this can be   *
!* extended by exchange routines.                                     *
!* We Compute on the resident nodes only. This means loops over       *
!* IP=1,NP_RES and backwards. This means that we do not have to do    *
!* exchanges of ASPAR values. Only the resident nodes are sent.       *
!* This is smaller and this is all that we ever need.                 *
!*                                                                    *
!* WRONG APPROACHES:                                                  *
!* to use all nodes 1:MNP may look simpler but it forces to have the  *
!* following property of the NNZ, IA, JA arrays. If two vertices      *
!* i and j are adjacent in a grid G, then they are adjacent in ANY    *
!* of the grid in which they are both contained.                      *
!* This property is actually not satisfied in general.                *
!* We could extend the IA, JA arrays by                               *
!* adding some vertices but that looks quite hazardous idea and it    *
!* it is actually not needed by the ILU0 preconditioner and other     *
!*                                                                    *
!* PROBLEM:                                                           *
!* There is an asymmetry in the construction of the ordering.         *
!* We start from low colors and upwards. If we had started with       *
!* high colors and gone downwards, then we get a different ordering   *
!* (even if we take the opposite, because of the interface nodes)     *
!* This requires the construction of many mappings.                   *
!* Our approach is actually to rebuild separate node sets.            *
!*                                                                    *
!* CHECKS:                                                            *
!* ---The sum of number of non-zero entries in Jstatus_L over all     *
!*    nodes should be equal to the sum of number of non-zero entries  *
!*    of Jstatus_U over all nodes.                                    *
!*    This is because number of upper diagonal entries should be      *
!*    equal to number of lower diagonal entries.                      *
!* ---We CANNOT have Jstatus_L(J)=1 and Jstatus_U(J)=1, i.e. a matrix *
!*    entry cannot be both lower and upper.                           *
!* ---We have sum of Jstatus_L + sum J_status_U + np_global should    *
!*    be equal to NNZ_global                                          *
!*                                                                    *
!* So, procedure is as follows:                                       *
!* ---compute ASPAR on 1,NP_RES nodes and no synchronization          *
!* ---compute on IP=1,NP_RES for L solving                            *
!* ---export to grids of higher rank, the values on nodes 1,NP_RES    *
!*    only. The other ghost points have invalid values or are         *
!*    resident of other grids of lower rank. (at this stage, some     *
!*    ghost values are wrong but are not exported)                    *
!* ---export to grid of lower rank in order to correct their ghost    *
!*    values and get the value                                        *
!* ---compute on IP=NP_RES,1,-1 for U solving                         *
!* ---export to grid of lower rank.                                   *
!*    Do everything similarly to L solve.                             *
!*                                                                    *
!* The basic approach is that we compute at a node S if and only if   *
!* it is a resident node. The twist come because some nodes are       *
!* resident for TWO domains. This is why we have the CovLower         *
!* We need to create disjoint domains for each node, so that          *
!* we have sum   sum(CovLower) = np_global                            *
!*                                                                    *
!* At the end of the resolution of the system, we need to do the      *
!* synchronization with respect to the unused nodes.                  *
!* Mystery?: When we apply the function, we do it on 1:NP_RES and     *
!* then call synchronizer. So, this means we need to do the           *
!* synchronization after the call to the preconditioner.              *
!* But there may be space for improvements here.                      *
!*                                                                    *
!* Description of specific exchange arrays:                           *
!* ---wwm_p2dsend_type/wwm_p2drecv_type                               *
!*    The points of 1:NP_RES are sent to nodes that contained them    *
!*    length=1                                                        *
!* ---wwmtot_p2dsend_type/wwmtot_p2drecv_type                         *
!*    same as above but length=MSC*MDC                                *
!* ---blk_p2dsend_type/blk_p2drecv_type                               *
!*    same as above but length=maxBlockLength for matrix exchanges    *
!* ---wwmmat_p2dsend_type/wwmmat_p2drecv_type                         *
!*    exchange of correct matrix elements, i.e. elements A(I,J)       *
!*    with I<=NP_RES                                                  *
!*    length is 1.                                                    *
!* ---u2l_p2dsend_type/u2l_p2drecv_type                               *
!*    upper 2 lower exchange arrays, depends on CovLower, so on       *
!*    the coloring chosen. length=maxBlockLength                      *
!*    exchange are from upper to lower.                               *
!* ---sync_p2dsend_type/sync_p2drecv_type                             *
!*    synchronize value, i.e. the CovLower=1 values  are send to all  *
!*    nodes. Length is MSC*MDC                                        *
!**********************************************************************
!* Mystery of the mpi_exchange routines.                              *
!* One way to have derived data types is to do                        *
!* call mpi_type_create_indexed_block(nbCommon,1,dspl_send,           *
!*                                    rtype,eType,ierr)               *
!* call mpi_type_commit(eType, ierr)                                  *
!* here 1 refers to the block length.                                 *
!* The problem is how to make exchanges. The standard method          *
!* mpi_isend(eMes, 1, eType, ....)                                    *
!* and so to send only 1 vector. A priori, it seems impossible        *
!* to send several vectors together because when we created the       *
!* derived data type, the length was not precised in this creation    *
!* one way could be to use mpi_type_resized, but this seems quite     *
!* inoperative or create seg-fault. So, we do not know how to make    *
!* exchanges of the form mpi_isend(AC, MSC*MDC, eType, ....)          *
!* and this means that we cannot send AC(MNP, MSC, MDC) simply.       *
!* Instead, we have to use U(MSC,MDC,MNP) in order to use the above   *
!**********************************************************************
MODULE WWM_PARALL_SOLVER
#if defined WWM_SOLVER && defined MPI_PARALL_GRID
      TYPE Graph
         integer nbVert
         integer MaxDeg
         integer nbEdge
         integer, dimension(:), pointer :: ListDegree
         integer, dimension(:,:), pointer :: ListEdge
      END TYPE Graph
      CONTAINS
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE EXCHANGE_P4D_WWM_TR(AC)
      USE DATAPOOL, only : MNP, MSC, MDC, rkind
      USE elfe_msgp, only : exchange_p4d_wwm
      implicit none
      real(rkind), intent(inout) :: AC(MNP,MSC,MDC)
      real(rkind) :: U(MSC,MDC,MNP)
      INTEGER :: IP, IS, ID
      DO IP = 1, MNP
        DO IS = 1, MSC
          DO ID = 1, MDC
            U(IS,ID,IP) = AC(IP,IS,ID)
          END DO
        END DO
      END DO
      CALL EXCHANGE_P4D_WWM(U)
      DO IP = 1, MNP
        DO IS = 1, MSC
          DO ID = 1, MDC
            AC(IP,IS,ID) = U(IS,ID,IP)
          END DO
        END DO
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_EXCHANGE_P4D_WWM(LocalColor, AC)
      USE DATAPOOL, only : MNP, MSC, MDC, rkind, LocalColorInfo
      USE DATAPOOL, only : wwm_nnbr_send, wwm_nnbr_recv
      USE DATAPOOL, only : wwm_ListNeigh_send, wwm_ListNeigh_recv
      USE DATAPOOL, only : wwmtot_p2dsend_type, wwmtot_p2drecv_type
      USE DATAPOOL, only : wwm_p2dsend_rqst, wwm_p2drecv_rqst
      USE DATAPOOL, only : wwm_p2dsend_stat, wwm_p2drecv_stat
      USE elfe_msgp, only : comm, ierr, myrank
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      real(rkind), intent(inout) :: AC(LocalColor%MSCeffect,MDC,MNP)
      integer iSync, iRank
      DO iSync=1,wwm_nnbr_send
        iRank=wwm_ListNeigh_send(iSync)
        CALL mpi_isend(AC, 1, wwmtot_p2dsend_type(iSync), iRank-1, 1020, comm, wwm_p2dsend_rqst(iSync), ierr)
      END DO
      DO iSync=1,wwm_nnbr_recv
        iRank=wwm_ListNeigh_recv(iSync)
        call mpi_irecv(AC,1,wwmtot_p2drecv_type(iSync),iRank-1,1020,comm,wwm_p2drecv_rqst(iSync),ierr)
      END DO
      IF (wwm_nnbr_send > 0) THEN
        call mpi_waitall(wwm_nnbr_send, wwm_p2dsend_rqst, wwm_p2dsend_stat,ierr)
      END IF
      IF (wwm_nnbr_recv > 0) THEN
        call mpi_waitall(wwm_nnbr_recv, wwm_p2drecv_rqst, wwm_p2drecv_stat,ierr)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_EXCHANGE_ASPAR(LocalColor, ASPAR_bl)
      USE DATAPOOL
      USE elfe_msgp, only : comm, ierr, myrank
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      real(rkind), intent(inout) :: ASPAR_bl(LocalColor%MSCeffect,MDC,NNZ)
      integer I, iProc
      INTEGER :: IZ, IS, ID
      do I=1,wwm_nnbr_m_send
        iProc=wwm_ListNeigh_m_send(I)
        call mpi_isend(ASPAR_bl,1,wwmmat_p2dsend_type(I),iProc-1,991,comm,wwmmat_p2dsend_rqst(i),ierr)
      enddo
      do I=1,wwm_nnbr_m_recv
        iProc=wwm_ListNeigh_m_recv(I)
        call mpi_irecv(ASPAR_bl,1,wwmmat_p2drecv_type(I),iProc-1,991,comm,wwmmat_p2drecv_rqst(i),ierr)
      enddo
      IF (wwm_nnbr_m_recv .gt. 0) THEN
        call mpi_waitall(wwm_nnbr_m_recv,wwmmat_p2drecv_rqst,wwmmat_p2drecv_stat,ierr)
      END IF
      IF (wwm_nnbr_m_send .gt. 0) THEN
        call mpi_waitall(wwm_nnbr_m_send,wwmmat_p2dsend_rqst,wwmmat_p2dsend_stat,ierr)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SYMM_GRAPH_BUILD_PROCESSOR_ADJACENCY(AdjGraph)
      USE DATAPOOL, only : wwm_nnbr, wwm_ListNeigh
      USE elfe_msgp, only : myrank
      implicit none
      type(Graph), intent(inout) :: AdjGraph
# ifdef DEBUG
      WRITE(740+myrank,*) 'wwm_nnbr=', wwm_nnbr
# endif
      CALL KERNEL_GRAPH_BUILD_PROCESSOR_ADJACENCY(AdjGraph, wwm_nnbr, wwm_ListNeigh)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GRAPH_BUILD_PROCESSOR_ADJACENCY(AdjGraph)
      USE elfe_msgp, only : nnbr_p, nbrrank_p
      implicit none
      type(Graph), intent(inout) :: AdjGraph
      integer :: ListNe(nnbr_p)
      integer I
      DO I=1,nnbr_p
        ListNe(I)=nbrrank_p(I)+1
      END DO
      CALL KERNEL_GRAPH_BUILD_PROCESSOR_ADJACENCY(AdjGraph, nnbr_p, ListNe)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE KERNEL_GRAPH_BUILD_PROCESSOR_ADJACENCY(AdjGraph, nb, ListNe)
      USE elfe_msgp, only : myrank, nproc, ierr, comm, istatus, itype
      USE DATAPOOL
      implicit none
      integer, intent(in) :: nb
      integer, intent(in) :: ListNe(nb)
      type(Graph), intent(inout) :: AdjGraph
      integer, allocatable :: rbuf_int(:)
      integer ierror, I, iProc
      integer idx, eDeg, nbEdge, iEdge
      integer istat
      AdjGraph % nbVert=nproc
# ifdef DEBUG
      WRITE(740+myrank,*) 'myrank=', myrank, ' KERNEL_GRAPH_BUILD_PROC...'
# endif
      IF (myrank.eq.0) THEN
        allocate(AdjGraph % ListDegree(nproc), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 1')
        AdjGraph % ListDegree(1)=nb
        allocate(rbuf_int(1), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 2')
        DO iProc=2,nproc
          CALL MPI_RECV(rbuf_int,1,itype, iProc-1, 19, comm, istatus, ierror)
          AdjGraph % ListDegree(iProc)=rbuf_int(1)
        END DO
        deallocate(rbuf_int)
        nbEdge=0
        DO iProc=1,nproc
          nbEdge=nbEdge + AdjGraph % ListDegree(iProc)
        END DO
        AdjGraph % nbEdge=nbEdge
# ifdef DEBUG
        WRITE(740+myrank,*) 'nbEdge=', nbEdge
# endif
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
# ifdef DEBUG
        WRITE(740+myrank,*) 'ListEdge exported'
# endif
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
# ifdef DEBUG
        WRITE(740+myrank,*) 'nbEdge=', nbEdge
# endif
        !
        allocate(rbuf_int(nproc), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 11')
        CALL MPI_RECV(rbuf_int,nproc,itype, 0, 32, comm, istatus, ierr)
        allocate(AdjGraph % ListDegree(nproc), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 12')
        AdjGraph % ListDegree=rbuf_int
        deallocate(rbuf_int)
# ifdef DEBUG
        WRITE(740+myrank,*) 'ListDegree assigned'
# endif
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
# ifdef DEBUG
        WRITE(740+myrank,*) 'ListEdge assigned'
# endif
      ENDIF
      AdjGraph % MaxDeg=maxval(AdjGraph % ListDegree)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GRAPH_TEST_CONNECT(AdjGraph, result)
      implicit none
      type(Graph), intent(in) :: AdjGraph
      integer, intent(out) :: result
      integer :: ListStatus(AdjGraph%nbVert)
      integer :: ListPosFirst(AdjGraph%nbVert)
      integer idx, iVert, nbVert, nbVertIsFinished, eAdj
      integer eDeg, sizConn, I, IsFinished
      integer istat
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
# ifdef DEBUG
      SUBROUTINE GRAPH_TEST_UNIDIRECT(AdjGraph)
      implicit none
      type(Graph), intent(in) :: AdjGraph
      integer nbEdge, iEdge, jEdge, jEdgeF
      integer eVert1, eVert2, fVert1, fVert2
      nbEdge=AdjGraph%nbEdge
      DO iEdge=1,nbEdge
        eVert1=AdjGraph%ListEdge(iEdge,1)
        eVert2=AdjGraph%ListEdge(iEdge,2)
        jEdgeF=-1
        DO jEdge=1,nbEdge
          fVert1=AdjGraph%ListEdge(jEdge,1)
          fVert2=AdjGraph%ListEdge(jEdge,2)
          IF ((eVert1.eq.fVert2).and.(eVert2.eq.fVert1)) THEN
            jEdgeF=jEdge
          END IF
        END DO
        IF (jEdgeF.eq.-1) THEN
          CALL WWM_ABORT('Error at test of symmetry')
        END IF
      END DO
      END SUBROUTINE
# endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COLLECT_ALL_IPLG
      USE DATAPOOL
      USE elfe_msgp
      USE elfe_glbl
      implicit none
      integer, allocatable :: rbuf_int(:)
      integer len, iProc, IP, idx, sumMNP
      integer istat
      allocate(ListMNP(nproc), ListNP_RES(nproc), rbuf_int(2), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 16')
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
# ifdef DEBUG
      DO iProc=1,nproc
        WRITE(740+myrank,*) 'iProc=', iProc, ' np_res=', ListNP_RES(iProc)
      END DO
      WRITE(740+myrank,*) 'max(np_res)=', maxval(ListNP_RES)
      WRITE(740+myrank,*) 'nnproc=', nproc
      WRITE(740+myrank,*) 'sumMNP=', sumMNP
# endif
      allocate(ListIPLG(sumMNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 17')
      IF (myrank == 0) THEN
        idx=0
        DO IP=1,MNP
          idx=idx+1
          ListIPLG(idx)=iplg(IP)
        END DO
        DO iProc=2,nproc
          len=ListMNP(iProc)
          allocate(rbuf_int(len), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 18')
          CALL MPI_RECV(rbuf_int,len,itype, iProc-1, 269, comm, istatus, ierr)
          DO IP=1,len
            idx=idx+1
            ListIPLG(idx)=rbuf_int(IP)
          END DO
          deallocate(rbuf_int)
        END DO
# ifdef DEBUG
        IF (idx /= sumMNP) THEN
          CALL WWM_ABORT('Inconsistency in IPLG creation')
        END IF
# endif
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
      SUBROUTINE COLLECT_ALL_COVLOWER(LocalColor)
      USE DATAPOOL
      USE elfe_msgp
      USE elfe_glbl
      implicit none
      type(LocalColorInfo), intent(inout) :: LocalColor
      integer, allocatable :: rbuf_int(:)
      integer len, iProc, IP, idx, sumMNP
      integer istat
      sumMNP=sum(ListMNP)
      allocate(LocalColor % ListCovLower(sumMNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 19')
      IF (myrank == 0) THEN
        idx=0
        DO IP=1,MNP
          idx=idx+1
          LocalColor % ListCovLower(idx)=LocalColor % CovLower(IP)
        END DO
        DO iProc=2,nproc
          len=ListMNP(iProc)
          allocate(rbuf_int(len), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 20')
          CALL MPI_RECV(rbuf_int,len,itype, iProc-1, 809, comm, istatus, ierr)
          DO IP=1,len
            idx=idx+1
            LocalColor % ListCovLower(idx)=rbuf_int(IP)
          END DO
          deallocate(rbuf_int)
        END DO
# ifdef DEBUG
        IF (idx /= sumMNP) THEN
          CALL WWM_ABORT('Inconsistency in CovLower creation')
        END IF
# endif
        DO iProc=2,nproc
          CALL MPI_SEND(LocalColor % ListCovLower,sumMNP,itype, iProc-1, 811, comm, ierr)
        END DO
      ELSE
        CALL MPI_SEND(LocalColor % CovLower,MNP,itype, 0, 809, comm, ierr)
        CALL MPI_RECV(LocalColor % ListCovLower,sumMNP,itype, 0, 811, comm, istatus, ierr)
      END IF
# ifdef DEBUG
      WRITE(740+myrank,*) 'COLLECT_ALL_COVLOWER'
      WRITE(740+myrank,*) 'sum(ListCovLower)=', sum(LocalColor % ListCovLower)
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COLLECT_ALL_IA_JA
      USE DATAPOOL
      USE elfe_msgp
      USE elfe_glbl
      implicit none
      integer, allocatable :: rbuf_int(:)
      integer len, iProc, IP, idx, sumMNP
      integer sumIAsiz, sumNNZ
      integer istat
      allocate(ListNNZ(nproc), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 21')
      !
      ! Collecting NNZ
      !
      IF (myrank == 0) THEN
        ListNNZ(1)=NNZ
        allocate(rbuf_int(1), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 22')
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
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 23')
        rbuf_int(1)=NNZ
        CALL MPI_SEND(rbuf_int,1,itype, 0, 257, comm, ierr)
        deallocate(rbuf_int)
        CALL MPI_RECV(ListNNZ,nproc,itype, 0, 263, comm, istatus, ierr)
      END IF
      !
      ! Collecting IA
      !
      sumIAsiz=sum(ListMNP) + nproc
# ifdef DEBUG
      DO iProc=1,nproc
        WRITE(740+myrank,*) 'iProc=', iProc, 'nnz=', ListNNZ(iProc)
      END DO
      WRITE(740+myrank,*) 'sumIAsiz=', sumIAsiz
# endif
      allocate(ListIA(sumIAsiz), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 24')
      IF (myrank == 0) THEN
        idx=0
        DO IP=1,MNP+1
          idx=idx+1
          ListIA(idx)=IA(IP)
        END DO
        DO iProc=2,nproc
          len=ListMNP(iProc)+1
          allocate(rbuf_int(len), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 25')
          CALL MPI_RECV(rbuf_int,len,itype, iProc-1, 269, comm, istatus, ierr)
          DO IP=1,len
            idx=idx+1
            ListIA(idx)=rbuf_int(IP)
          END DO
          deallocate(rbuf_int)
        END DO
# ifdef DEBUG
        IF (idx /= sumIAsiz) THEN
          CALL WWM_ABORT('Inconsistency in sumIAsiz')
        END IF
# endif
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
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 26')
      IF (myrank == 0) THEN
        idx=0
        DO IP=1,NNZ
          idx=idx+1
          ListJA(idx)=JA(IP)
        END DO
        DO iProc=2,nproc
          len=ListNNZ(iProc)
          allocate(rbuf_int(len), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 27')
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
# ifdef DEBUG
      SUBROUTINE TEST_ASPAR_SYMMETRY
      USE DATAPOOL, only : MNP, IA, JA
      USE elfe_msgp, only : myrank
      implicit none
      integer :: IsSymm
      integer IP, JP, IPB, J, J2, JFOUND
      IsSymm=1
      DO IP=1,MNP
        DO J=IA(IP),IA(IP+1)-1
          JP=JA(J)
          JFOUND=-1
          DO J2=IA(JP),IA(JP+1)-1
            IPb=JA(J2)
            IF (IPb == IP) THEN
              JFOUND=J2
            END IF
          END DO
          IF (JFOUND == -1) THEN
            IsSymm=0
          END IF
        END DO
      END DO
      IF (IsSymm == 0) THEN
        WRITE(740+myrank,*) 'NNZ_IA_JA is NOT symmetric'
      ELSE
        WRITE(740+myrank,*) 'NNZ_IA_JA IS symmetric'
      END IF
      END SUBROUTINE
# endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CREATE_WWM_P2D_EXCH(MSCeffect)
      USE DATAPOOL
      USE elfe_msgp
      USE elfe_glbl
      implicit none
      include 'mpif.h'
      integer, intent(in) :: MSCeffect
      integer :: ListFirst(nproc)
      integer :: ListNeigh01(nproc)
      integer :: ListCommon_recv(nproc)
      integer :: ListCommon_send(nproc)
      integer :: ListMapped(np_global)
      integer :: ListMappedB(np_global)
      integer, allocatable :: dspl_send(:), dspl_recv(:)
      integer, allocatable :: dspl_send_tot(:), dspl_recv_tot(:)
      integer IP, IP_glob, iProc, WeMatch, MNPloc, idx, NP_RESloc
      integer iNeigh, IPmap, eSize, eSizeRed, nbCommon
      integer nbCommon_send, nbCommon_recv, idx_send, idx_recv
      integer sumNbCommon_send, sumNbCommon_recv
      integer idxDspl_send, idxDspl_recv
      integer eExtent, eExtentRed, NewExtent, eLB, sizRType
      integer eType1, eType2
      integer istat
      ListFirst=0
      DO iProc=2,nproc
        ListFirst(iProc)=ListFirst(iProc-1) + ListMNP(iProc-1)
      END DO
      ListMapped=0
      DO IP=1,MNP
        IP_glob=iplg(IP)
# ifdef DEBUG
        IF (ListMapped(IP_glob) .gt. 0) THEN
          CALL WWM_ABORT('Clear error in ListMapped');
        ENDIF
# endif
        ListMapped(IP_glob)=IP
      END DO
      wwm_nnbr_send=0
      wwm_nnbr_recv=0
      ListCommon_send=0
      ListCommon_recv=0
      DO iProc=1,nproc
        IF (iPROC .ne. myrank+1) THEN
          MNPloc=ListMNP(iProc)
          NP_RESloc=ListNP_RES(iProc)
          ListMappedB=0
          DO IP=1,MNPloc
            IP_glob=ListIPLG(IP+ListFirst(iProc))
# ifdef DEBUG
            IF (ListMappedB(IP_glob) .gt. 0) THEN
              CALL WWM_ABORT('Clear error in ListMappedB');
            ENDIF
# endif
            ListMappedB(IP_glob)=IP
          END DO
          !
          nbCommon_recv=0
          DO IP=1,NP_RESloc
            IP_glob=ListIPLG(IP+ListFirst(iProc))
            IF (ListMapped(IP_glob).gt.0) THEN
              nbCommon_recv=nbCommon_recv+1
            END IF
          END DO
          IF (nbCommon_recv .gt. 0) THEN
            wwm_nnbr_recv=wwm_nnbr_recv+1
            ListCommon_recv(iProc)=nbCommon_recv
          END IF
          !
          nbCommon_send=0
          DO IP=1,NP_RES
            IP_glob=iplg(IP)
            IF (ListMappedB(IP_glob).gt.0) THEN
              nbCommon_send=nbCommon_send+1
            END IF
          END DO
          IF (nbCommon_send .gt. 0) THEN
            wwm_nnbr_send=wwm_nnbr_send+1
            ListCommon_send(iProc)=nbCommon_send
          END IF
        END IF
      END DO
# ifdef DEBUG
      WRITE(740+myrank,*) 'WWM_P2D: wwm_nnbr_send=', wwm_nnbr_send
      WRITE(740+myrank,*) 'WWM_P2D: wwm_nnbr_recv=', wwm_nnbr_recv
# endif
      allocate(wwm_ListNbCommon_send(wwm_nnbr_send), wwm_ListNbCommon_recv(wwm_nnbr_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 28')
      allocate(wwm_ListNeigh_send(wwm_nnbr_send), wwm_ListNeigh_recv(wwm_nnbr_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 29')
      idx_send=0
      idx_recv=0
      sumNbCommon_send=0
      sumNbCommon_recv=0
      DO iProc=1,nproc
        IF (ListCommon_send(iProc) .gt. 0) THEN
          idx_send=idx_send+1
          wwm_ListNeigh_send(idx_send)=iProc
          nbCommon=ListCommon_send(iProc)
          wwm_ListNbCommon_send(idx_send)=nbCommon
          sumNbCommon_send=sumNbCommon_send+nbCommon
        END IF
        IF (ListCommon_recv(iProc) .gt. 0) THEN
          idx_recv=idx_recv+1
          wwm_ListNeigh_recv(idx_recv)=iProc
          nbCommon=ListCommon_recv(iProc)
          wwm_ListNbCommon_recv(idx_recv)=nbCommon
          sumNbCommon_recv=sumNbCommon_recv+nbCommon
        END IF
      END DO
      allocate(wwm_ListDspl_send(sumNbCommon_send), wwm_ListDspl_recv(sumNbCommon_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 30')
      !
      ! Now the symmetric exchanges for color computations
      ! 
      wwm_nnbr=0
      DO iProc=1,nproc
        IF ((ListCommon_send(iProc).gt.0).or.(ListCommon_recv(iProc).gt.0)) THEN
          wwm_nnbr=wwm_nnbr+1
        END IF
      END DO
      allocate(wwm_ListNeigh(wwm_nnbr), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 31')
      idx=0
      DO iProc=1,nproc
        IF ((ListCommon_send(iProc).gt.0).or.(ListCommon_recv(iProc).gt.0)) THEN
          idx=idx+1
          wwm_ListNeigh(idx)=iProc
        END IF
      END DO
# ifdef DEBUG
      WRITE(740+myrank,*) 'WWM_P2D: wwm_nnbr=', wwm_nnbr
      WRITE(740+myrank,*) 'WWM_P2D: wwm_ListNeigh built'
# endif
      allocate(wwm_p2dsend_rqst(wwm_nnbr_send), wwm_p2drecv_rqst(wwm_nnbr_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 32')
      allocate(wwm_p2dsend_stat(MPI_STATUS_SIZE,wwm_nnbr_send), wwm_p2drecv_stat(MPI_STATUS_SIZE,wwm_nnbr_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 33')
      allocate(wwm_p2dsend_type(wwm_nnbr_send), wwm_p2drecv_type(wwm_nnbr_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 34')
      allocate(wwmtot_p2dsend_type(wwm_nnbr_send), wwmtot_p2drecv_type(wwm_nnbr_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 35')
# ifdef DEBUG
      WRITE(740+myrank,*) 'WWM_P2D: alloc done'
# endif
      idxDspl_send=0
      DO iNeigh=1,wwm_nnbr_send
        iProc=wwm_ListNeigh_send(iNeigh)
        MNPloc=ListMNP(iProc)
        NP_RESloc=ListNP_RES(iProc)
        ListMappedB=0
        DO IP=1,MNPloc
          IP_glob=ListIPLG(IP+ListFirst(iProc))
# ifdef DEBUG
          IF (ListMappedB(IP_glob) .gt. 0) THEN
            CALL WWM_ABORT('Clear error in ListMappedB');
          ENDIF
# endif
          ListMappedB(IP_glob)=IP
        END DO
        nbCommon=wwm_ListNbCommon_send(iNeigh)
        allocate(dspl_send(nbCommon), dspl_send_tot(nbCommon), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 36')
        idx=0
        DO IP=1,NP_RES
          IP_glob=iplg(IP)
          IF (ListMappedB(IP_glob).gt.0) THEN
            IPmap=ListMappedB(IP_glob)
            idx=idx+1
            dspl_send(idx)=IP-1
            dspl_send_tot(idx)=MSCeffect*MDC*(IP-1)
            idxDspl_send=idxDspl_send+1
            wwm_ListDspl_send(idxDspl_send)=IP
          END IF
        END DO
        call mpi_type_create_indexed_block(nbCommon,1,dspl_send,rtype,wwm_p2dsend_type(iNeigh), ierr)
        call mpi_type_commit(wwm_p2dsend_type(iNeigh), ierr)
        call mpi_type_create_indexed_block(nbCommon,MSCeffect*MDC,dspl_send_tot,rtype,wwmtot_p2dsend_type(iNeigh), ierr)
        call mpi_type_commit(wwmtot_p2dsend_type(iNeigh), ierr)
        deallocate(dspl_send)
        deallocate(dspl_send_tot)
      END DO
      idxDspl_recv=0
      DO iNeigh=1,wwm_nnbr_recv
        iProc=wwm_ListNeigh_recv(iNeigh)
        MNPloc=ListMNP(iProc)
        NP_RESloc=ListNP_RES(iProc)
        nbCommon=wwm_ListNbCommon_recv(iNeigh)
        allocate(dspl_recv(nbCommon), dspl_recv_tot(nbCommon), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 37')
        idx=0
        DO IP=1,NP_RESloc
          IP_glob=ListIPLG(IP+ListFirst(iProc))
          IF (ListMapped(IP_glob).gt.0) THEN
            IPmap=ListMapped(IP_glob)
            idx=idx+1
            dspl_recv(idx)=IPmap-1
            dspl_recv_tot(idx)=MSCeffect*MDC*(IPmap-1)
            idxDspl_recv=idxDspl_recv+1
            wwm_ListDspl_recv(idxDspl_recv)=IPmap
          END IF
        END DO
        !
        call mpi_type_create_indexed_block(nbCommon,1,dspl_recv,rtype,wwm_p2drecv_type(iNeigh),ierr)
        call mpi_type_commit(wwm_p2drecv_type(iNeigh), ierr)
        call mpi_type_create_indexed_block(nbCommon,MSCeffect*MDC,dspl_recv_tot,rtype,wwmtot_p2drecv_type(iNeigh),ierr)
        call mpi_type_commit(wwmtot_p2drecv_type(iNeigh), ierr)
        deallocate(dspl_recv)
        deallocate(dspl_recv_tot)
      END DO
# ifdef DEBUG
      WRITE(740+myrank,*) 'Leaving CREATE_WWM_P2D_EXCH'
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CREATE_WWM_MAT_P2D_EXCH(MSCeffect)
      USE DATAPOOL
      USE elfe_msgp
      USE elfe_glbl
      implicit none
      include 'mpif.h'
      integer, intent(in) :: MSCeffect
      integer :: ListFirstMNP(nproc), ListFirstNNZ(nproc)
      integer :: ListCommon_send(nproc), ListCommon_recv(nproc)
      integer :: ListMapped(np_global)
      integer :: ListMappedB(np_global)
      integer, allocatable :: dspl_send(:), dspl_recv(:)
      integer nbCommon_send, nbCommon_recv
      integer IAfirst
      integer IP, JP, I, J, J2, IP_glob, JP_glob, iProc
      integer WeMatch, MNPloc, NP_RESloc, JP_j
      integer IPloc, JPloc, Jfound, idx
      integer iNeigh, IPmap, nbCommon, nbCommonB, eSize, eSizeRed
      integer sumNbCommon_send, sumNbCommon_recv
      integer idxDspl_send, idxDspl_recv
      integer istat
      ListFirstNNZ=0
      ListFirstMNP=0
      DO iProc=2,nproc
        ListFirstMNP(iProc)=ListFirstMNP(iProc-1) + ListMNP(iProc-1)
        ListFirstNNZ(iProc)=ListFirstNNZ(iProc-1) + ListNNZ(iProc-1)
      END DO
      ListMapped=0
      DO IP=1,MNP
        IP_glob=iplg(IP)
# ifdef DEBUG
        IF (ListMapped(IP_glob) .gt. 0) THEN
          CALL WWM_ABORT('Clear error in ListMapped');
        ENDIF
# endif
        ListMapped(IP_glob)=IP
      END DO
      wwm_nnbr_m_recv=0
      wwm_nnbr_m_send=0
      ListCommon_recv=0
      DO I=1,wwm_nnbr_recv
        iProc=wwm_ListNeigh_recv(I)
        NP_RESloc=ListNP_RES(iProc)
        MNPloc=ListMNP(iProc)
        ListMappedB=0
        DO IP=1,MNPloc
          IP_glob=ListIPLG(IP+ListFirstMNP(iProc))
# ifdef DEBUG
          IF (ListMappedB(IP_glob) .gt. 0) THEN
            CALL WWM_ABORT('Clear error in ListMappedB');
          ENDIF
# endif
          ListMappedB(IP_glob)=IP
        END DO
        nbCommon_recv=0
        DO IP=1,NP_RESloc
          IP_glob=ListIPLG(IP+ListFirstMNP(iProc))
          IPloc=ListMapped(IP_glob)
          IF (IPloc.gt.0) THEN
            IAfirst=ListFirstMNP(iProc) + iProc-1
            DO J=ListIA(IP+IAfirst),ListIA(IP+IAfirst+1)-1
              JP=ListJA(J+ListFirstNNZ(iProc))
              JP_glob=ListIPLG(JP+ListFirstMNP(iProc))
              JPloc=ListMapped(JP_glob)
              IF (JPloc.gt.0) THEN
                JFOUND=-1
                DO J2=IA(IPloc),IA(IPloc+1)-1
                  IF (JA(J2) == JPloc) THEN
                    JFOUND=J2
                  END IF
                END DO
                IF (JFOUND .gt. 0) THEN
                  nbCommon_recv=nbCommon_recv+1
                END IF
              END IF
            END DO
          END IF
        END DO
        IF (nbCommon_recv .gt. 0) THEN
          wwm_nnbr_m_recv=wwm_nnbr_m_recv+1
          ListCommon_recv(iProc)=nbCommon_recv
        END IF
      END DO
      ListCommon_send=0
      DO I=1,wwm_nnbr_send
        iProc=wwm_ListNeigh_send(I)
        NP_RESloc=ListNP_RES(iProc)
        MNPloc=ListMNP(iProc)
        ListMappedB=0
        DO IP=1,MNPloc
          IP_glob=ListIPLG(IP+ListFirstMNP(iProc))
# ifdef DEBUG
          IF (ListMappedB(IP_glob) .gt. 0) THEN
            CALL WWM_ABORT('Clear error in ListMappedB');
          ENDIF
# endif
          ListMappedB(IP_glob)=IP
        END DO
        nbCommon_send=0
        DO IP=1,NP_RES
          IP_glob=iplg(IP)
          IPloc=ListMappedB(IP_glob)
          IF (IPloc.gt.0) THEN
            DO J=IA(IP),IA(IP+1)-1
              JP=JA(J)
              JP_glob=iplg(JP)
              JPloc=ListMappedB(JP_glob)
              IF (JPloc.gt.0) THEN
                IAfirst=ListFirstMNP(iProc) + iProc-1
                JFOUND=-1
                DO J2=ListIA(IPloc+IAfirst),ListIA(IPloc+IAfirst+1)-1
                  JP_j=ListJA(J2+ListFirstNNZ(iProc))
                  IF (JP_j == JPloc) THEN
                    JFOUND=J2
                  END IF
                END DO
                IF (JFOUND .gt. 0) THEN
                  nbCommon_send=nbCommon_send+1
                END IF
              END IF
            END DO
          END IF
        END DO
        IF (nbCommon_send .gt. 0) THEN
          wwm_nnbr_m_send=wwm_nnbr_m_send+1
          ListCommon_send(iProc)=nbCommon_send
        END IF
      END DO
# ifdef DEBUG
      WRITE(740+myrank,*) 'WWM_MAT_P2D: wwm_nnbr_m_send=', wwm_nnbr_m_send
      WRITE(740+myrank,*) 'WWM_MAT_P2D: wwm_nnbr_m_recv=', wwm_nnbr_m_recv
# endif
      allocate(wwmmat_p2dsend_rqst(wwm_nnbr_m_send), wwmmat_p2drecv_rqst(wwm_nnbr_m_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 38')
      allocate(wwmmat_p2dsend_stat(MPI_STATUS_SIZE,wwm_nnbr_m_send), wwmmat_p2drecv_stat(MPI_STATUS_SIZE,wwm_nnbr_m_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 39')
      allocate(wwmmat_p2dsend_type(wwm_nnbr_m_send), wwmmat_p2drecv_type(wwm_nnbr_m_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 40')
      allocate(wwm_ListNbCommon_m_send(wwm_nnbr_m_send), wwm_ListNbCommon_m_recv(wwm_nnbr_m_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 41')
      allocate(wwm_ListNeigh_m_recv(wwm_nnbr_m_recv), wwm_ListNeigh_m_send(wwm_nnbr_m_send), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 42')
# ifdef DEBUG
      WRITE(740+myrank,*) 'WWM_MAT_P2D: alloc done'
# endif
      idx=0
      sumNbCommon_send=0
      DO iProc=1,nproc
        IF (ListCommon_send(iProc) .gt. 0) THEN
          idx=idx+1
          wwm_ListNeigh_m_send(idx)=iProc
          nbCommon=ListCommon_send(iProc)
          wwm_ListNbCommon_m_send(idx)=nbCommon
          sumNbCommon_send=sumNbCommon_send+nbCommon
        END IF
      END DO
# ifdef DEBUG
      WRITE(740+myrank,*) 'WWM_MAT_P2D: wwm_ListNeigh_m_send built'
# endif
      idx=0
      sumNbCommon_recv=0
      DO iProc=1,nproc
        IF (ListCommon_recv(iProc) .gt. 0) THEN
          idx=idx+1
          wwm_ListNeigh_m_recv(idx)=iProc
          nbCommon=ListCommon_recv(iProc)
          wwm_ListNbCommon_m_recv(idx)=nbCommon
          sumNbCommon_recv=sumNbCommon_recv+nbCommon
        END IF
      END DO
# ifdef DEBUG
      WRITE(740+myrank,*) 'WWM_MAT_P2D: wwm_ListNeigh_m_recv built'
# endif
      allocate(wwm_ListDspl_m_send(sumNbCommon_send), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 43')
      idxDspl_send=0
      DO iNeigh=1,wwm_nnbr_m_send
        iProc=wwm_ListNeigh_m_send(iNeigh)
        MNPloc=ListMNP(iProc)
        ListMappedB=0
        DO IP=1,MNPloc
          IP_glob=ListIPLG(IP+ListFirstMNP(iProc))
# ifdef DEBUG
          IF (ListMappedB(IP_glob) .gt. 0) THEN
            CALL WWM_ABORT('Clear error in ListMappedB');
          ENDIF
# endif
          ListMappedB(IP_glob)=IP
        END DO
        nbCommon=wwm_ListNbCommon_m_send(iNeigh)
        allocate(dspl_send(nbCommon), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 44')
        idx=0
        DO IP=1,NP_RES
          IP_glob=iplg(IP)
          IPloc=ListMappedB(IP_glob)
          IF (IPloc.gt.0) THEN
            DO J=IA(IP),IA(IP+1)-1
              JP=JA(J)
              JP_glob=iplg(JP)
              JPloc=ListMappedB(JP_glob)
              IF (JPloc .gt. 0) THEN
                IAfirst=ListFirstMNP(iProc) + iProc-1
                JFOUND=-1
                DO J2=ListIA(IPloc+IAfirst),ListIA(IPloc+IAfirst+1)-1
                  JP_j=ListJA(J2+ListFirstNNZ(iProc))
                  IF (JP_j == JPloc) THEN
                    JFOUND=J2
                  END IF
                END DO
                IF (JFOUND .gt. 0) THEN
                  idxDspl_send=idxDspl_send+1
                  wwm_ListDspl_m_send(idxDspl_send)=J
                  idx=idx+1
                  dspl_send(idx)=(J-1)*MSCeffect*MDC
                END IF
              END IF
            END DO
          END IF
        END DO
        call mpi_type_create_indexed_block(nbCommon,MSCeffect*MDC,dspl_send,rtype,wwmmat_p2dsend_type(iNeigh), ierr)
        call mpi_type_commit(wwmmat_p2dsend_type(iNeigh), ierr)
        deallocate(dspl_send)
      END DO
# ifdef DEBUG
      WRITE(740+myrank,*) 'sumNbCommon_recv=', sumNbCommon_recv
# endif
      allocate(wwm_ListDspl_m_recv(sumNbCommon_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 45')
      idxDspl_recv=0
      DO iNeigh=1,wwm_nnbr_m_recv
        iProc=wwm_ListNeigh_m_recv(iNeigh)
        nbCommon=wwm_ListNbCommon_m_recv(iNeigh)
# ifdef DEBUG
        WRITE(740+myrank,*) 'nbCommon=', nbCommon
# endif
        allocate(dspl_recv(nbCommon), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 46')
        NP_RESloc=ListNP_RES(iProc)
        idx=0
        DO IP=1,NP_RESloc
          IP_glob=ListIPLG(IP+ListFirstMNP(iProc))
          IPloc=ListMapped(IP_glob)
          IF (IPloc.gt.0) THEN
            IAfirst=ListFirstMNP(iProc) + iProc-1
            DO J=ListIA(IP+IAfirst),ListIA(IP+IAfirst+1)-1
              JP=ListJA(J+ListFirstNNZ(iProc))
              JP_glob=ListIPLG(JP+ListFirstMNP(iProc))
              JPloc=ListMapped(JP_glob)
              IF (JPloc.gt.0) THEN
                JFOUND=-1
                DO J2=IA(IPloc),IA(IPloc+1)-1
                  IF (JA(J2) == JPloc) THEN
                    JFOUND=J2
                  END IF
                END DO
                IF (JFOUND /= -1) THEN
                  idxDspl_recv=idxDspl_recv+1
                  wwm_ListDspl_m_recv(idxDspl_recv)=Jfound
                  idx=idx+1
                  dspl_recv(idx)=(Jfound-1)*MSCeffect*MDC
                END IF
              END IF
            END DO
          END IF
        END DO
        call mpi_type_create_indexed_block(nbCommon,MSCeffect*MDC,dspl_recv,rtype,wwmmat_p2drecv_type(iNeigh), ierr)
        call mpi_type_commit(wwmmat_p2drecv_type(iNeigh), ierr)
        deallocate(dspl_recv)
      END DO
# ifdef DEBUG
      WRITE(740+myrank,*) 'Leaving CREATE_WWM_MAT_P2D_EXCH'
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
# ifdef DEBUG
      SUBROUTINE CHECK_STANDARD_SELFE_EXCH
      USE elfe_msgp, only : exchange_p2d, myrank
      USE DATAPOOL
      REAL(rkind) :: XPcopy(MNP)
      REAL(rkind) :: SumErr
      XPcopy=XP
      CALL exchange_p2d(XP)
      SumErr=sum(abs(XPcopy - XP))
!      WRITE(740+myrank,*) 'SELFE_EXCH SumErr=', SumErr
      END SUBROUTINE
# endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
# ifdef DEBUG
      SUBROUTINE WRITE_EXPLICIT_ORDERING(ListPos, ListPosRev, ListColor)
      USE DATAPOOL, only : MNP, ListMNP, ListNP_RES, ListIPLG
      USE elfe_msgp, only : nproc, myrank
      USE elfe_glbl, only : np_global
      IMPLICIT NONE
      integer, intent(inout) :: ListPos(np_global)
      integer, intent(inout) :: ListPosRev(np_global)
      integer, intent(in) :: ListColor(nproc)
      integer :: ListFirst(nproc)
      integer :: ListTotal(np_global)
      integer minColor, maxColor, eColor
      integer iProc, idx, IP, IPglob
      minColor=minval(ListColor)
      maxColor=maxval(ListColor)
      ListTotal=1
      ListFirst=0
      DO iProc=2,nproc
        ListFirst(iProc)=ListFirst(iProc-1) + ListMNP(iProc-1)
      END DO
      ListPos=0
      ListPosRev=0
      idx=0
      DO eColor=minColor,maxColor
        DO iProc=1,nproc
          IF (ListColor(iProc) .eq. eColor) THEN
            DO IP=1,ListNP_RES(iProc)
              IPglob=ListIPLG(IP+ListFirst(iProc))
              IF (ListTotal(IPglob) .eq. 1) THEN
                idx=idx+1
                ListPos(idx)=IPglob
                ListPosRev(IPglob)=idx
                ListTotal(IPglob)=0
              END IF
            END DO
          END IF
        END DO
      END DO
      IF (minval(ListPosRev) .eq. 0) THEN
        CALL WWM_ABORT('Please correct ')
      END IF
      IF (idx .ne. np_global) THEN
        DO IP=1,np_global
          WRITE(myrank+540,*) 'IP/Pos/PosR=', IP, ListPos(IP), ListPosRev(IP)
        END DO
        CALL WWM_ABORT('One more bug to solve')
      END IF
      END SUBROUTINE
# endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE BUILD_MULTICOLORING(AdjGraph, ListColor)
      USE elfe_msgp, only : myrank
      implicit none
      type(Graph), intent(in) :: AdjGraph
      integer, intent(out) :: ListColor(AdjGraph%nbVert)
      integer, allocatable :: CurrColor(:)
      integer MaxDeg, iVert, eVert, eColor, eDeg
      integer idx, I, ChromaticNr, nbVert
      integer, allocatable :: ListPosFirst(:)
      integer, allocatable :: TheOrdering(:)
      integer eColorF, iVertFound, eAdjColor
      integer nbUndef, MinDeg, eAdj, MinUndef, PosMin
      integer istat
      MaxDeg=AdjGraph % MaxDeg
      nbVert=AdjGraph % nbVert
      allocate(CurrColor(MaxDeg+1), ListPosFirst(nbVert), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 53')
      ListColor=0
      idx=0
      DO iVert=1,nbVert
        ListPosFirst(iVert)=idx
        idx=idx+AdjGraph % ListDegree(iVert)
      END DO
      MinDeg=MaxDeg+3
      PosMin=-1
      DO iVert=1,nbVert
        eDeg=AdjGraph % ListDegree(iVert)
        IF (eDeg .lt. MinDeg) THEN
          MinDeg=eDeg
          PosMin=iVert
        END IF
      END DO
      idx=ListPosFirst(PosMin)
      DO I=0,MinDeg
        IF (I.eq.0) THEN
          eVert=PosMin
        ELSE
          eVert=AdjGraph % ListEdge(idx+I,2)
        END IF
        ListColor(eVert)=I+1
      END DO
      DO
        MinUndef=nbVert
        iVertFound=0
        DO iVert=1,nbVert
          IF (ListColor(iVert) == 0) THEN
            idx=ListPosFirst(iVert)
            eDeg=AdjGraph % ListDegree(iVert)
            nbUndef=0
            DO I=1,eDeg
              eAdj=AdjGraph % ListEdge(idx+I,2)
              eAdjColor=ListColor(eAdj)
              IF (eAdjColor == 0) THEN
                nbUndef=nbUndef+1
              END IF
            END DO
            IF (nbUndef .lt. MinUndef) THEN
              MinUndef=nbUndef
              iVertFound=iVert
            END IF
          END IF
        END DO
        IF (iVertFound == 0) THEN
          EXIT
        END IF
        eDeg=AdjGraph % ListDegree(iVertFound)
        idx=ListPosFirst(iVertFound)
        CurrColor=0
        DO I=1,eDeg
          eVert=AdjGraph % ListEdge(idx+I,2)
          eColor=ListColor(eVert)
          IF (eColor.gt.0) THEN
            CurrColor(eColor)=1
          END IF
        END DO
        eColorF=-1
        DO I=1,MaxDeg+1
          IF (eColorF == -1) THEN
            IF (CurrColor(I) == 0) THEN
              eColorF=I
            END IF
          END IF
        END DO
        ListColor(iVertFound)=eColorF
      END DO
      deallocate(ListPosFirst)
      deallocate(CurrColor)
      ChromaticNr=maxval(ListColor)
# ifdef DEBUG
      WRITE(740+myrank,*) 'ChromaticNr=', ChromaticNr
      DO iVert=1,nbVert
        WRITE(740+myrank,*) 'iVert=', iVert, 'eColor=', ListColor(iVert)
      END DO
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DeallocateGraph(TheGraph)
      implicit none
      type(Graph), intent(inout) :: TheGraph
      deallocate(TheGraph % ListDegree)
      deallocate(TheGraph % ListEdge)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_BLOCK_FREQDIR(LocalColor, Nblock, MSCeffect)
      USE DATAPOOL, only : MNP, MDC, LocalColorInfo, rkind
      USE DATAPOOL, only : wwm_nnbr_send, wwm_nnbr_recv
      USE DATAPOOL, only : wwm_ListNbCommon_send, wwm_ListNbCommon_recv
      USE elfe_msgp, only : rtype, ierr, myrank
      USE elfe_glbl, only : iplg
      USE DATAPOOL, only : XP, YP
      implicit none
      type(LocalColorInfo), intent(inout) :: LocalColor
      integer, intent(in) :: MSCeffect
      integer, intent(in) :: Nblock
      integer Ntot, Hlen, Delta, iBlock, idx, ID, IS
      integer lenBlock, maxBlockLength
      integer IC, nbCommon, eFirst, I, idxSend, idxRecv
      integer istat
      REAL(rkind) :: eFieldStackA(MNP,3)
      Ntot=MyREAL(MSCeffect*MDC)
      Hlen=INT(Ntot/Nblock)
      Delta=Ntot - Hlen*Nblock
      iBlock=1
      idx=1
      LocalColor % Nblock=Nblock
      IF (Delta == 0) THEN
        maxBlockLength=Hlen
      ELSE
        maxBlockLength=Hlen+1
      ENDIF
      allocate(LocalColor % ISindex(Nblock, maxBlockLength), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 55')
      allocate(LocalColor % IDindex(Nblock, maxBlockLength), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 56')
      DO IS=1,MSCeffect
        DO ID=1,MDC
          LocalColor % ISindex(iBlock, idx)=IS
          LocalColor % IDindex(iBlock, idx)=ID
          IF (iBlock <= Delta) THEN
            lenBlock=Hlen+1
          ELSE
            lenBlock=Hlen
          END IF
          idx=idx+1
          IF (idx > lenBlock) THEN
            iBlock=iBlock+1
            idx=1
          ENDIF
        END DO
      END DO
      allocate(LocalColor % BlockLength(Nblock), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 57')
      DO iBlock=1,Nblock
        IF (iBlock <= Delta) THEN
          lenBlock=Hlen+1
        ELSE
          lenBlock=Hlen
        END IF
        LocalColor % BlockLength(iBlock)=lenBlock
      END DO
      LocalColor % maxBlockLength = maxBlockLength
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_BLK_L2U_ARRAY(LocalColor)
      USE DATAPOOL, only : MNP, MSC, MDC, LocalColorInfo, rkind
      USE DATAPOOL, only : wwm_nnbr_send, wwm_nnbr_recv
      USE DATAPOOL, only : wwm_ListNbCommon_send, wwm_ListNbCommon_recv
      USE DATAPOOL, only : wwm_ListDspl_send, wwm_ListDspl_recv
      USE DATAPOOL, only : wwm_ListNeigh_send, wwm_ListNeigh_recv
      USE elfe_msgp, only : rtype, ierr, myrank
      USE elfe_glbl, only : iplg
      USE DATAPOOL, only : XP, YP
      implicit none
      type(LocalColorInfo), intent(inout) :: LocalColor
# ifdef LU_SOLVE_RWRT
      integer, allocatable :: ListNeed(:), IdxRev(:)
      integer nbNeedSend_blk, nbNeedRecv_blk
      integer idx, IP
# endif
      integer nbUpp_send, nbLow_recv, iUpp, iLow, iRank
      integer maxBlockLength, idxSend, idxRecv
      integer I, IC, nbCommon, eFirst
      integer ListFirstCommon_send(wwm_nnbr_send)
      integer ListFirstCommon_recv(wwm_nnbr_recv)
      integer, allocatable :: dspl_send(:), dspl_recv(:)
      integer istat
      maxBlockLength=LocalColor % maxBlockLength
      ListFirstCommon_send=0
      DO I=2,wwm_nnbr_send
        ListFirstCommon_send(I)=ListFirstCommon_send(I-1)+wwm_ListNbCommon_send(I-1)
      END DO
      ListFirstCommon_recv=0
      DO I=2,wwm_nnbr_recv
        ListFirstCommon_recv(I)=ListFirstCommon_recv(I-1)+wwm_ListNbCommon_recv(I-1)
      END DO
# ifdef L2U_OPER
      nbUpp_send=LocalColor % nbUpp_send
      nbLow_recv=LocalColor % nbLow_recv
      allocate(LocalColor % l2u_p2dsend_type(nbUpp_send), LocalColor % l2u_p2drecv_type(nbLow_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 59')
      allocate(LocalColor % l2u_ListNeigh_send(nbUpp_send), LocalColor % l2u_ListNeigh_recv(nbLow_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 60')
# else
      allocate(LocalColor % blk_p2dsend_type(wwm_nnbr_send), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 59')
      allocate(LocalColor % blk_p2drecv_type(wwm_nnbr_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 60')
# endif
# ifdef DEBUG
      WRITE(740+myrank,*) 'maxBlockLength=', maxBlockLength
# endif
# ifdef LU_SOLVE_RWRT
#  ifdef L2U_OPER
      allocate(ListNeed(MNP), IdxRev(MNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 61a')
      ListNeed=0
      IdxRev=0
      nbNeedSend_blk=0
      DO iUpp=1,nbUpp_send
        I=LocalColor % ListIdxUpper_send(iUpp)
        iRank=wwm_ListNeigh_send(I)-1
        LocalColor % l2u_ListNeigh_send(iUpp)=iRank
        nbCommon=wwm_ListNbCommon_send(I)
        eFirst=ListFirstCommon_send(I)
        DO IC=1,nbCommon
          idxSend=wwm_ListDspl_send(eFirst+IC)
          IF (ListNeed(idxSend) .eq. 0) THEN
            ListNeed(idxSend)=1
            nbNeedSend_blk=nbNeedSend_blk+1
          END IF
        END DO
      END DO
      LocalColor % nbNeedSend_blk=nbNeedSend_blk
      allocate(LocalColor % IdxSend_blk(nbNeedSend_blk), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 61b')
      idx=0
      DO IP=1,MNP
        IF (ListNeed(IP) .eq. 1) THEN
          idx=idx+1
          LocalColor % IdxSend_blk(idx)=IP
          IdxRev(IP)=idx
        END IF
      END DO
      DO iUpp=1,nbUpp_send
        I=LocalColor % ListIdxUpper_send(iUpp)
        nbCommon=wwm_ListNbCommon_send(I)
        eFirst=ListFirstCommon_send(I)
        ALLOCATE(dspl_send(nbCommon), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 61')
        DO IC=1,nbCommon
          idxSend=wwm_ListDspl_send(eFirst+IC)
          dspl_send(IC)=(IdxRev(idxSend)-1)*maxBlockLength
        END DO
        call mpi_type_create_indexed_block(nbCommon,maxBlockLength,dspl_send,rtype,LocalColor % l2u_p2dsend_type(iUpp),ierr)
        call mpi_type_commit(LocalColor % l2u_p2dsend_type(iUpp), ierr)
        DEALLOCATE(dspl_send)
      END DO
      !
      !
      ListNeed=0
      IdxRev=0
      nbNeedRecv_blk=0
      DO iLow=1,nbLow_recv
        I=LocalColor % ListIdxLower_recv(iLow)
        iRank=wwm_ListNeigh_recv(I)-1
        LocalColor % l2u_ListNeigh_recv(iLow)=iRank
        nbCommon=wwm_ListNbCommon_recv(I)
        eFirst=ListFirstCommon_recv(I)
        DO IC=1,nbCommon
          idxRecv=wwm_ListDspl_recv(eFirst+IC)
          IF (ListNeed(idxRecv) .eq. 0) THEN
            ListNeed(idxRecv)=1
            nbNeedRecv_blk=nbNeedRecv_blk+1
          END IF
        END DO
      END DO
      LocalColor % nbNeedRecv_blk=nbNeedRecv_blk
      allocate(LocalColor % IdxRecv_blk(nbNeedRecv_blk), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 61b')
      idx=0
      DO IP=1,MNP
        IF (ListNeed(IP) .eq. 1) THEN
          idx=idx+1
          LocalColor % IdxRecv_blk(idx)=IP
          IdxRev(IP)=idx
        END IF
      END DO
      DO iLow=1,nbLow_recv
        I=LocalColor % ListIdxLower_recv(iLow)
        nbCommon=wwm_ListNbCommon_recv(I)
        eFirst=ListFirstCommon_recv(I)
        ALLOCATE(dspl_recv(nbCommon), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 62')
        DO IC=1,nbCommon
          idxRecv=wwm_ListDspl_recv(eFirst+IC)
          dspl_recv(IC)=(IdxRev(idxRecv)-1)*maxBlockLength
        END DO
        call mpi_type_create_indexed_block(nbCommon,maxBlockLength,dspl_recv,rtype,LocalColor % l2u_p2drecv_type(iLow),ierr)
        call mpi_type_commit(LocalColor % l2u_p2drecv_type(iLow), ierr)
        DEALLOCATE(dspl_recv)
      END DO
      deallocate(ListNeed, IdxRev)
#  else
      allocate(ListNeed(MNP), IdxRev(MNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 61a')
      ListNeed=0
      IdxRev=0
      nbNeedSend_blk=0
      DO I=1,wwm_nnbr_send
        nbCommon=wwm_ListNbCommon_send(I)
        eFirst=ListFirstCommon_send(I)
        DO IC=1,nbCommon
          idxSend=wwm_ListDspl_send(eFirst+IC)
          IF (ListNeed(idxSend) .eq. 0) THEN
            ListNeed(idxSend)=1
            nbNeedSend_blk=nbNeedSend_blk+1
          END IF
        END DO
      END DO
      LocalColor % nbNeedSend_blk=nbNeedSend_blk
      allocate(LocalColor % IdxSend_blk(nbNeedSend_blk), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 61b')
      idx=0
      DO IP=1,MNP
        IF (ListNeed(IP) .eq. 1) THEN
          idx=idx+1
          LocalColor % IdxSend_blk(idx)=IP
          IdxRev(IP)=idx
        END IF
      END DO
      DO I=1,wwm_nnbr_send
        nbCommon=wwm_ListNbCommon_send(I)
        eFirst=ListFirstCommon_send(I)
        ALLOCATE(dspl_send(nbCommon), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 61')
        DO IC=1,nbCommon
          idxSend=wwm_ListDspl_send(eFirst+IC)
          dspl_send(IC)=(IdxRev(idxSend)-1)*maxBlockLength
        END DO
        call mpi_type_create_indexed_block(nbCommon,maxBlockLength,dspl_send,rtype,LocalColor % blk_p2dsend_type(I),ierr)
        call mpi_type_commit(LocalColor % blk_p2dsend_type(I), ierr)
        DEALLOCATE(dspl_send)
      END DO
      !
      !
      ListNeed=0
      IdxRev=0
      nbNeedRecv_blk=0
      DO I=1,wwm_nnbr_recv
        nbCommon=wwm_ListNbCommon_recv(I)
        eFirst=ListFirstCommon_recv(I)
        DO IC=1,nbCommon
          idxRecv=wwm_ListDspl_recv(eFirst+IC)
          IF (ListNeed(idxRecv) .eq. 0) THEN
            ListNeed(idxRecv)=1
            nbNeedRecv_blk=nbNeedRecv_blk+1
          END IF
        END DO
      END DO
      LocalColor % nbNeedRecv_blk=nbNeedRecv_blk
      allocate(LocalColor % IdxRecv_blk(nbNeedRecv_blk), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 61b')
      idx=0
      DO IP=1,MNP
        IF (ListNeed(IP) .eq. 1) THEN
          idx=idx+1
          LocalColor % IdxRecv_blk(idx)=IP
          IdxRev(IP)=idx
        END IF
      END DO
      DO I=1,wwm_nnbr_recv
        nbCommon=wwm_ListNbCommon_recv(I)
        eFirst=ListFirstCommon_recv(I)
        ALLOCATE(dspl_recv(nbCommon), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 62')
        DO IC=1,nbCommon
          idxRecv=wwm_ListDspl_recv(eFirst+IC)
          dspl_recv(IC)=(IdxRev(idxRecv)-1)*maxBlockLength
        END DO
        call mpi_type_create_indexed_block(nbCommon,maxBlockLength,dspl_recv,rtype,LocalColor % blk_p2drecv_type(I),ierr)
        call mpi_type_commit(LocalColor % blk_p2drecv_type(I), ierr)
        DEALLOCATE(dspl_recv)
      END DO
      deallocate(ListNeed, IdxRev)
#  endif
# else
      DO I=1,wwm_nnbr_send
        nbCommon=wwm_ListNbCommon_send(I)
        eFirst=ListFirstCommon_send(I)
        ALLOCATE(dspl_send(nbCommon), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 61')
        DO IC=1,nbCommon
          idxSend=wwm_ListDspl_send(eFirst+IC)
          dspl_send(IC)=(idxSend-1)*maxBlockLength
        END DO
        call mpi_type_create_indexed_block(nbCommon,maxBlockLength,dspl_send,rtype,LocalColor % blk_p2dsend_type(I),ierr)
        call mpi_type_commit(LocalColor % blk_p2dsend_type(I), ierr)
        DEALLOCATE(dspl_send)
      END DO
      DO I=1,wwm_nnbr_recv
        nbCommon=wwm_ListNbCommon_recv(I)
        eFirst=ListFirstCommon_recv(I)
        ALLOCATE(dspl_recv(nbCommon), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 62')
        DO IC=1,nbCommon
          idxRecv=wwm_ListDspl_recv(eFirst+IC)
          dspl_recv(IC)=(idxRecv-1)*maxBlockLength
        END DO
        call mpi_type_create_indexed_block(nbCommon,maxBlockLength,dspl_recv,rtype,LocalColor % blk_p2drecv_type(I),ierr)
        call mpi_type_commit(LocalColor % blk_p2drecv_type(I), ierr)
        DEALLOCATE(dspl_recv)
      END DO
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SYMM_INIT_COLORING(LocalColor, NbBlock, MSCeffect)
      USE DATAPOOL, only : LocalColorInfo, MNP, rkind, XP, YP
      USE DATAPOOL, only : DO_SOLVE_L, DO_SOLVE_U
      USE DATAPOOL, only : DO_SYNC_UPP_2_LOW, DO_SYNC_LOW_2_UPP, DO_SYNC_FINAL
      USE elfe_msgp, only : myrank, nproc
      USE elfe_glbl, only : iplg
      implicit none
      type(LocalColorInfo), intent(inout) :: LocalColor
      integer, intent(in) :: MSCeffect
      integer, intent(in) :: NbBlock
      type(Graph) :: AdjGraph
      integer :: ListColor(nproc)
      integer :: ListColorWork(nproc)
      real(rkind) :: eFieldStackA(MNP,3)
      real(rkind) :: eFieldStackB(MNP,2)
      real(rkind) :: eFieldStackC(MNP,1)
      real(rkind) :: eFieldStackRevA(3,MNP)
      real(rkind) :: eFieldStackRevB(2,MNP)
      real(rkind) :: eFieldStackRevC(1,MNP)
      integer TheRes, istat
# ifdef DEBUG
      CALL TEST_ASPAR_SYMMETRY
      WRITE(740+myrank,*) 'After TEST_ASPAR_SYMMETRY'
# endif
# ifdef DEBUG
      CALL COMPUTE_TOTAL_INDEX_SHIFT(TheRes)
      WRITE(740+myrank,*) 'Total residual shift=', TheRes
# endif
      CALL COLLECT_ALL_IPLG
# ifdef DEBUG
      WRITE(740+myrank,*) 'After COLLECT_ALL_IPLG'
# endif
      CALL COLLECT_ALL_IA_JA
# ifdef DEBUG
      WRITE(740+myrank,*) 'After COLLECT_ALL_IA_JA'
# endif
      CALL CREATE_WWM_P2D_EXCH(MSCeffect)
# ifdef DEBUG
      WRITE(740+myrank,*) 'After CREATE_WWM_P2D_EXCH'
# endif
      CALL CREATE_WWM_MAT_P2D_EXCH(MSCeffect)
# ifdef DEBUG
      WRITE(740+myrank,*) 'After CREATE_WWM_MAT_P2D_EXCH'
# endif
      CALL SYMM_GRAPH_BUILD_PROCESSOR_ADJACENCY(AdjGraph)
# ifdef DEBUG
      WRITE(740+myrank,*) 'After SYMM_GRAPH_BUILD_PROCESSOR_ADJACENCY'
      CALL GRAPH_TEST_UNIDIRECT(AdjGraph)
      WRITE(740+myrank,*) 'After GRAPH_TEST_UNIDIRECT'
# endif
      CALL BUILD_MULTICOLORING(AdjGraph, ListColor)
# ifdef DEBUG
      WRITE(740+myrank,*) 'After BUILD_MULTICOLORING'
# endif
      CALL DeallocateGraph(AdjGraph)
      ListColorWork=-ListColor
      allocate(LocalColor % ListColor(nproc), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 63')
      LocalColor % ListColor=ListColorWork
# ifdef DEBUG
      WRITE(740+myrank,*) 'Before INIT_LOW_2_UPP_ARRAYS'
# endif
      CALL INIT_LOW_2_UPP_ARRAYS(LocalColor, ListColorWork)
      !
# ifdef DEBUG
      WRITE(740+myrank,*) 'Before CALL_BLOCK_FREQDIR'
# endif
      CALL INIT_BLOCK_FREQDIR(LocalColor, NbBlock, MSCeffect)
      !
# ifdef DEBUG
      WRITE(740+myrank,*) 'Before INIT_BLK_L2U_ARRAY'
# endif
      CALL INIT_BLK_L2U_ARRAY(LocalColor)
      !
# ifdef DEBUG
      WRITE(740+myrank,*) 'Before COLLECT_ALL_COVLOWER'
# endif
      CALL COLLECT_ALL_COVLOWER(LocalColor)
      !
# ifdef DEBUG
      WRITE(740+myrank,*) 'Before INIT_COVLOWER_ARRAY'
# endif
      CALL INIT_COVLOWER_ARRAY(LocalColor)
      !
# ifdef DEBUG
      WRITE(740+myrank,*) 'Before DETERMINE_JSTATUS_L_U'
# endif
      CALL DETERMINE_JSTATUS_L_U(LocalColor)
      !
      DO_SOLVE_L=.TRUE.
      DO_SOLVE_U=.TRUE.
      DO_SYNC_UPP_2_LOW=.TRUE.
      DO_SYNC_LOW_2_UPP=.TRUE.
      DO_SYNC_FINAL=.TRUE.
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_LOW_2_UPP_ARRAYS(LocalColor, ListColor)
      USE DATAPOOL, only : LocalColorInfo, MNP, rkind
      USE DATAPOOL, only : NNZ, IA, JA, NP_RES
      USE DATAPOOL, only : wwm_ListNbCommon_send, wwm_ListNbCommon_recv
      USE DATAPOOL, only : wwm_nnbr_send, wwm_nnbr_recv
      USE DATAPOOL, only : wwm_p2drecv_type, wwm_p2dsend_type
      USE DATAPOOL, only : wwm_ListNeigh_send, wwm_ListNeigh_recv
      USE DATAPOOL, only : wwm_ListDspl_recv, wwm_ListDspl_send
      USE elfe_msgp, only : myrank, nproc, comm, ierr, nbrrank_p
      implicit none
      include 'mpif.h'
      type(LocalColorInfo), intent(inout) :: LocalColor
      integer, intent(in) :: ListColor(nproc)
      real(rkind) :: p2d_data_send(MNP)
      real(rkind) :: CovLower(MNP), CovLower_meth2(MNP), CovLower_meth3(MNP)
      real(rkind) :: SumErr, SumDiff
      integer eColor, fColor, I, iRank, J
      integer nbLow_send, nbUpp_send, nbLow_recv, nbUpp_recv
      integer iProc, stat, eSize, iLow, iUpp, DoOper
      integer IC, eFirst, nbCommon, IP, IPloc, JP
      integer ListFirstCommon_send(wwm_nnbr_send)
      integer ListFirstCommon_recv(wwm_nnbr_recv)
      integer istat
      eColor=ListColor(myrank+1)
      nbUpp_send=0
      DO I=1,wwm_nnbr_send
        iRank=wwm_ListNeigh_send(I)
        fColor=ListColor(iRank)
# ifdef DEBUG
        WRITE(740+myrank,*) 'I=', I, 'iRank=', iRank, 'fColor=', fColor
        IF (fColor.eq.eColor) THEN
          call wwm_abort('Major error in the code')
        END IF
# endif
        IF (fColor.gt.eColor) THEN
          nbUpp_send=nbUpp_send + 1
        ENDIF
      END DO
      nbLow_recv=0
      DO I=1,wwm_nnbr_recv
        iRank=wwm_ListNeigh_recv(I)
        fColor=ListColor(iRank)
# ifdef DEBUG
        IF (fColor.eq.eColor) THEN
          call wwm_abort('Major error in the code')
        END IF
# endif
        IF (fColor.lt.eColor) THEN
          nbLow_recv=nbLow_recv + 1
        ENDIF
      END DO
      ListFirstCommon_send=0
      DO I=2,wwm_nnbr_send
        ListFirstCommon_send(I)=ListFirstCommon_send(I-1)+wwm_ListNbCommon_send(I-1)
      END DO
      ListFirstCommon_recv=0
      DO I=2,wwm_nnbr_recv
        ListFirstCommon_recv(I)=ListFirstCommon_recv(I-1)+wwm_ListNbCommon_recv(I-1)
      END DO
# ifdef DEBUG
      WRITE(740+myrank,*) 'SIC: nbLow_recv=', nbLow_recv, ' nbUpp_send=', nbUpp_send
# endif
      LocalColor % nbUpp_send=nbUpp_send
      LocalColor % nbLow_recv=nbLow_recv
      allocate(LocalColor % ListIdxUpper_send(nbUpp_send), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 64')
      allocate(LocalColor % ListIdxLower_recv(nbLow_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 65')
      iUpp=0
      DO I=1,wwm_nnbr_send
        iRank=wwm_ListNeigh_send(I)
        fColor=ListColor(iRank)
        IF (fColor.gt.eColor) THEN
          iUpp=iUpp + 1
          LocalColor % ListIdxUpper_send(iUpp)=I
        ENDIF
      END DO
      iLow=0
      DO I=1,wwm_nnbr_recv
        iRank=wwm_ListNeigh_recv(I)
        fColor=ListColor(iRank)
        IF (fColor.lt.eColor) THEN
          iLow=iLow + 1
          LocalColor % ListIdxLower_recv(iLow)=I
        ENDIF
      END DO
      allocate(LocalColor % Upp_s_rq(nbUpp_send), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 66')
      allocate(LocalColor % Upp_s_stat(MPI_STATUS_SIZE, nbUpp_send), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 67')
      allocate(LocalColor % Low_r_rq(nbLow_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 68')
      allocate(LocalColor % Low_r_stat(MPI_STATUS_SIZE, nbLow_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 69')
      p2d_data_send=0
      CovLower=1
      CovLower_meth2=1
      DO iUpp=1,nbUpp_send
        I=LocalColor % ListIdxUpper_send(iUpp)
        iRank=wwm_ListNeigh_send(i)
# ifdef DEBUG
        WRITE(740+myrank,*) 'ISEND: iUpp=', iUpp, 'I=', I, 'iRank=', iRank
# endif
        call mpi_isend(p2d_data_send,1,wwm_p2dsend_type(i),iRank-1,13,comm,LocalColor % Upp_s_rq(iUpp),ierr)
      END DO
      DO iLow=1,nbLow_recv
        I=LocalColor % ListIdxLower_recv(iLow)
        iRank=wwm_ListNeigh_recv(i)
# ifdef DEBUG
        WRITE(740+myrank,*) 'IRECV: iLow=', iLow, 'I=', I, 'iRank=', iRank
# endif
        nbCommon=wwm_ListNbCommon_recv(I)
        eFirst=ListFirstCommon_recv(I)
        DO IC=1,nbCommon
          IPloc=wwm_ListDspl_recv(eFirst+IC)
          CovLower_meth2(IPloc)=0
        END DO
        call mpi_irecv(CovLower,1,wwm_p2drecv_type(i),iRank-1,13,comm,LocalColor % Low_r_rq(iLow),ierr)
      END DO
      IF (nbUpp_send > 0) THEN
        call mpi_waitall(nbUpp_send, LocalColor%Upp_s_rq, LocalColor%Upp_s_stat,ierr)
      END IF
      IF (nbLow_recv > 0) THEN
        call mpi_waitall(nbLow_recv, LocalColor%Low_r_rq, LocalColor%Low_r_stat,ierr)
      END IF
      allocate(LocalColor % CovLower(MNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 70')
      LocalColor % CovLower=INT(CovLower)
      SumErr=sum(abs(CovLower-CovLower_meth2))
# ifdef DEBUG
      WRITE(740+myrank,*) 'SumErr(meth1/meth2) CovLower=', SumErr
      WRITE(740+myrank,*) 'MNP=', MNP, ' sum(CovLower)=', sum(CovLower)
      IF (SumErr .gt. 0) THEN
        DO IP=1,MNP
          WRITE(740+myrank,*) IP, CovLower(IP), CovLower_meth2(IP)
        END DO
      ENDIF
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_COVLOWER_ARRAY(LocalColor)
      USE DATAPOOL
      USE elfe_msgp
      USE elfe_glbl
      implicit none
      include 'mpif.h'
      type(LocalColorInfo), intent(inout) :: LocalColor
      integer :: ListFirst(nproc)
      integer :: ListNeigh01(nproc)
      integer :: ListCommon_recv(nproc)
      integer :: ListCommon_send(nproc)
      integer :: ListMapped0(np_global)
      integer :: ListMapped1(np_global)
      integer :: ListMapped0_B(np_global)
      integer :: ListMapped1_B(np_global)
      integer, allocatable :: dspl_send(:), dspl_recv(:)
      integer IP, IP_glob, iProc, WeMatch, MNPloc, idx, NP_RESloc
      integer iNeigh, IPmap, eSize, eSizeRed, nbCommon
      integer nbCommon_send, nbCommon_recv, idx_send, idx_recv
      integer sumNbCommon_send, sumNbCommon_recv
      integer eExtent, eExtentRed, NewExtent, eLB, sizRType
      integer eType1, eType2
      integer u2l_nnbr_send, u2l_nnbr_recv
      integer sync_nnbr_send, sync_nnbr_recv
      integer eCov, eColor, fColor, iSync
      integer maxBlockLength
      integer nbCase1, nbCase2
      integer nbMap0, nbMap1
      integer DoOper
      integer istat
# ifdef LU_SOLVE_RWRT
      integer, allocatable :: ListNeedSend(:), ListNeedRecv(:), IdxRev(:)
      integer nbNeedSend_u2l, nbNeedRecv_u2l
      integer lenMNP
# endif
      ListFirst=0
      DO iProc=2,nproc
        ListFirst(iProc)=ListFirst(iProc-1) + ListMNP(iProc-1)
      END DO
      ListMapped0=0
      ListMapped1=0
      nbMap0=0
      nbMap1=0
      DO IP=1,MNP
        IP_glob=iplg(IP)
        eCov=LocalColor % CovLower(IP)
        IF (eCov == 0) THEN
          ListMapped0(IP_glob)=IP
          nbMap0=nbMap0+1
        END IF
        IF (eCov == 1) THEN
          ListMapped1(IP_glob)=IP
          nbMap1=nbMap1+1
        END IF
      END DO
# ifdef DEBUG
      WRITE(740+myrank,*) 'nbMap0=', nbMap0, ' nbMap1=', nbMap1
# endif
      !
      ! First the Upper to lower (u2l) block arrays 
      !
      u2l_nnbr_send=0
      u2l_nnbr_recv=0
      ListCommon_send=0
      ListCommon_recv=0
      eColor=LocalColor % ListColor(myrank+1)
# ifdef LU_SOLVE_RWRT
      allocate(ListNeedRecv(MNP), ListNeedSend(MNP), IdxRev(MNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 79a')
      ListNeedRecv=0
      ListNeedSend=0
      IdxRev=0
      nbNeedSend_u2l=0
      nbNeedRecv_u2l=0
# endif
# ifdef DEBUG
      WRITE(740+myrank,*) 'U2L eColor=', eColor
# endif
      DO iNeigh=1,wwm_nnbr
        iProc=wwm_ListNeigh(iNeigh)
        fColor=LocalColor % ListColor(iProc)
# ifdef DEBUG
        WRITE(740+myrank,*) 'U2L iNeigh=', iNeigh, ' iProc=', iProc, 'fColor=', fColor
# endif
        MNPloc=ListMNP(iProc)
        NP_RESloc=ListNP_RES(iProc)
        IF (fColor .ge. eColor) THEN
          nbCommon_recv=0
          DO IP=1,NP_RESloc
            IP_glob=ListIPLG(IP+ListFirst(iProc))
            eCov=LocalColor % ListCovLower(IP+ListFirst(iProc))
            IF (eCov .eq. 1) THEN
              IPmap=ListMapped0(IP_glob)
              WeMatch=0
              IF (IPmap .gt. 0) THEN
                WeMatch=1
              ELSE
                IPmap=ListMapped1(IP_glob)
                IF ((IPmap .gt. 0).and.(IPmap .gt. NP_RES)) THEN
                  WeMatch=1
                END IF
              END IF
              IF (WeMatch .eq. 1) THEN
                nbCommon_recv=nbCommon_recv+1
# ifdef LU_SOLVE_RWRT
                IF (ListNeedRecv(IPmap) .eq. 0) THEN
                  nbNeedRecv_u2l=nbNeedRecv_u2l+1
                  ListNeedRecv(IPmap)=1
                END IF
# endif
              END IF
            END IF
          END DO
          IF (nbCommon_recv .gt. 0) THEN
            u2l_nnbr_recv=u2l_nnbr_recv+1
            ListCommon_recv(iProc)=nbCommon_recv
          END IF
# ifdef DEBUG
          WRITE(740+myrank,*) '   U2L nbCommon_recv=', nbCommon_recv
# endif
        END IF
        IF (fColor .le. eColor) THEN
          nbCommon_send=0
          DO IP=1,MNPloc
            IP_glob=ListIPLG(IP+ListFirst(iProc))
            eCov=LocalColor % ListCovLower(IP+ListFirst(iProc))
            IPmap=ListMapped1(IP_glob)
            IF ((IPmap .gt. 0).and.(IPmap .le. NP_RES)) THEN
              IF ((eCov .eq. 0).or.(IP.gt.NP_RESloc)) THEN
                nbCommon_send=nbCommon_send+1
# ifdef LU_SOLVE_RWRT
                IF (ListNeedSend(IPmap) .eq. 0) THEN
                  nbNeedSend_u2l=nbNeedSend_u2l+1
                  ListNeedSend(IPmap)=1
                END IF
# endif
              END IF
            END IF
          END DO
          IF (nbCommon_send .gt. 0) THEN
            u2l_nnbr_send=u2l_nnbr_send+1
            ListCommon_send(iProc)=nbCommon_send
          END IF
# ifdef DEBUG
          WRITE(740+myrank,*) '   U2L nbCommon_send=', nbCommon_send
# endif
        END IF
      END DO
# ifdef DEBUG
      WRITE(740+myrank,*) 'WWM_P2D: u2l_nnbr_send=', u2l_nnbr_send
      WRITE(740+myrank,*) 'WWM_P2D: u2l_nnbr_recv=', u2l_nnbr_recv
# endif
      LocalColor % u2l_nnbr_send=u2l_nnbr_send
      LocalColor % u2l_nnbr_recv=u2l_nnbr_recv
      allocate(LocalColor % u2l_ListNbCommon_send(u2l_nnbr_send), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 79')
      allocate(LocalColor % u2l_ListNbCommon_recv(u2l_nnbr_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 80')
      allocate(LocalColor % u2l_ListNeigh_send(u2l_nnbr_send), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 81')
      allocate(LocalColor % u2l_ListNeigh_recv(u2l_nnbr_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 82')
      idx_send=0
      idx_recv=0
      DO iProc=1,nproc
        IF (ListCommon_send(iProc) .gt. 0) THEN
          idx_send=idx_send+1
          LocalColor % u2l_ListNeigh_send(idx_send)=iProc
          nbCommon=ListCommon_send(iProc)
          LocalColor % u2l_ListNbCommon_send(idx_send)=nbCommon
        END IF
        IF (ListCommon_recv(iProc) .gt. 0) THEN
          idx_recv=idx_recv+1
          LocalColor % u2l_ListNeigh_recv(idx_recv)=iProc
          nbCommon=ListCommon_recv(iProc)
          LocalColor % u2l_ListNbCommon_recv(idx_recv)=nbCommon
        END IF
      END DO
# ifdef LU_SOLVE_RWRT
      LocalColor % nbNeedSend_u2l=nbNeedSend_u2l
      LocalColor % nbNeedRecv_u2l=nbNeedRecv_u2l
      allocate(LocalColor % IdxSend_u2l(nbNeedSend_u2l), LocalColor % IdxRecv_u2l(nbNeedRecv_u2l), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 82a')
# endif
      !
      ! Now creating the u2l exchange
      ! 
# ifdef DEBUG
      WRITE(740+myrank,*) 'WWM_P2D: wwm_nnbr=', wwm_nnbr
      WRITE(740+myrank,*) 'WWM_P2D: wwm_ListNeigh built'
# endif
      allocate(LocalColor % u2l_p2dsend_rqst(u2l_nnbr_send), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 83')
      allocate(LocalColor % u2l_p2drecv_rqst(u2l_nnbr_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 84')
      allocate(LocalColor % u2l_p2dsend_stat(MPI_STATUS_SIZE,u2l_nnbr_send), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 85')
      allocate(LocalColor % u2l_p2drecv_stat(MPI_STATUS_SIZE,u2l_nnbr_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 86')
      allocate(LocalColor % u2l_p2dsend_type(u2l_nnbr_send), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 87')
      allocate(LocalColor % u2l_p2drecv_type(u2l_nnbr_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 88')
# ifdef DEBUG
      WRITE(740+myrank,*) 'WWM_P2D: alloc done'
# endif
      maxBlockLength=LocalColor % maxBlockLength
# ifdef LU_SOLVE_RWRT
      idx=0
      DO IP=1,MNP
        IF (ListNeedSend(IP) .eq. 1) THEN
          idx=idx+1
          LocalColor % IdxSend_u2l(idx)=IP
          IdxRev(IP)=idx
        END IF
      END DO
# endif
      DO iNeigh=1,u2l_nnbr_send
        iProc=LocalColor % u2l_ListNeigh_send(iNeigh)
        MNPloc=ListMNP(iProc)
        NP_RESloc=ListNP_RES(iProc)
        ListMapped0_B=0
        ListMapped1_B=0
        DO IP=1,MNPloc
          IP_glob=ListIPLG(IP+ListFirst(iProc))
          eCov=LocalColor % ListCovLower(IP+ListFirst(iProc))
          IF (eCov == 0) THEN
            ListMapped0_B(IP_glob)=IP
          END IF
          IF (eCov == 1) THEN
            ListMapped1_B(IP_glob)=IP
          END IF
        END DO
        nbCommon=LocalColor % u2l_ListNbCommon_send(iNeigh)
        allocate(dspl_send(nbCommon), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 89')
        idx=0
        DO IP=1,NP_RES
          IF (LocalColor % CovLower(IP) .eq. 1) THEN
            IP_glob=iplg(IP)
            IPmap=ListMapped0_B(IP_glob)
            WeMatch=0
            IF (IPmap .gt. 0) THEN
              WeMatch=1
            ELSE
              IPmap=ListMapped1_B(IP_glob)
              IF ((IPmap .gt. 0).and.(IPmap .gt. NP_RESloc)) THEN
                WeMatch=1
              END IF
            END IF
            IF (WeMatch .eq. 1) THEN
              idx=idx+1
# ifdef LU_SOLVE_RWRT
              dspl_send(idx)=maxBlockLength*(IdxRev(IP)-1)
# else
              dspl_send(idx)=maxBlockLength*(IP-1)
# endif
            END IF
          END IF
        END DO
# ifdef DEBUG
        IF (idx .ne. nbCommon) THEN
          CALL WWM_ABORT('error in u2l_p2dsend')
        END IF
        WRITE(740+myrank,*) '   U2L idx=', idx, ' nbCommon=', nbCommon
# endif
        call mpi_type_create_indexed_block(nbCommon,maxBlockLength,dspl_send,rtype,LocalColor % u2l_p2dsend_type(iNeigh), ierr)
        call mpi_type_commit(LocalColor % u2l_p2dsend_type(iNeigh), ierr)
        deallocate(dspl_send)
      END DO
# ifdef LU_SOLVE_RWRT
      IdxRev=0
      idx=0
      DO IP=1,MNP
        IF (ListNeedRecv(IP) .eq. 1) THEN
          idx=idx+1
          LocalColor % IdxRecv_u2l(idx)=IP
          IdxRev(IP)=idx
        END IF
      END DO
# endif
      DO iNeigh=1,u2l_nnbr_recv
        iProc=LocalColor % u2l_ListNeigh_recv(iNeigh)
        MNPloc=ListMNP(iProc)
        NP_RESloc=ListNP_RES(iProc)
        ListMapped0_B=0
        ListMapped1_B=0
        DO IP=1,MNPloc
          IP_glob=ListIPLG(IP+ListFirst(iProc))
          eCov=LocalColor % ListCovLower(IP+ListFirst(iProc))
          IF (eCov == 0) THEN
            ListMapped0_B(IP_glob)=IP
          END IF
          IF (eCov == 1) THEN
            ListMapped1_B(IP_glob)=IP
          END IF
        END DO
        nbCommon=LocalColor % u2l_ListNbCommon_recv(iNeigh)
        allocate(dspl_recv(nbCommon), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 90')
        idx=0
        DO IP=1,NP_RESloc
          IP_glob=ListIPLG(IP+ListFirst(iProc))
          eCov=LocalColor % ListCovLower(IP+ListFirst(iProc))
          IF (eCov == 1) THEN
            IPmap=ListMapped0(IP_glob)
            WeMatch=0
            IF (IPmap .gt. 0) THEN
              WeMatch=1
            ELSE
              IPmap=ListMapped1(IP_glob)
              IF ((IPmap .gt. 0).and.(IPmap .gt. NP_RES)) THEN
                WeMatch=1
              END IF
            END IF
            IF (WeMatch .eq. 1) THEN
              idx=idx+1
# ifdef LU_SOLVE_RWRT
              dspl_recv(idx)=maxBlockLength*(IdxRev(IPmap)-1)
# else
              dspl_recv(idx)=maxBlockLength*(IPmap-1)
# endif
            END IF
          END IF
        END DO
# ifdef DEBUG
        IF (idx .ne. nbCommon) THEN
          CALL WWM_ABORT('error in u2l_p2drecv')
        END IF
        WRITE(740+myrank,*) '   U2L idx=', idx, ' nbCommon=', nbCommon
# endif
        call mpi_type_create_indexed_block(nbCommon,maxBlockLength,dspl_recv,rtype,LocalColor % u2l_p2drecv_type(iNeigh),ierr)
        call mpi_type_commit(LocalColor % u2l_p2drecv_type(iNeigh), ierr)
        deallocate(dspl_recv)
      END DO
# ifdef LU_SOLVE_RWRT
      deallocate(ListNeedRecv, ListNeedSend, IdxRev)
# endif
      !
      ! Now the synchronization arrays
      !
      sync_nnbr_send=0
      sync_nnbr_recv=0
      ListCommon_send=0
      ListCommon_recv=0
# ifdef DEBUG
      WRITE(740+myrank,*) 'wwm_nnbr=', wwm_nnbr
      WRITE(740+myrank,*) 'sum(ListCovLower)=', sum(LocalColor % ListCovLower)
# endif
      DO iNeigh=1,wwm_nnbr
        iProc=wwm_ListNeigh(iNeigh)
# ifdef DEBUG
        WRITE(740+myrank,*) 'iNeigh=', iNeigh, 'iProc=', iProc
# endif
        fColor=LocalColor % ListColor(iProc)
        MNPloc=ListMNP(iProc)
        NP_RESloc=ListNP_RES(iProc)
        ListMapped0_B=0
        ListMapped1_B=0
        DO IP=1,MNPloc
          IP_glob=ListIPLG(IP+ListFirst(iProc))
          eCov=LocalColor % ListCovLower(IP+ListFirst(iProc))
          IF (eCov == 0) THEN
            ListMapped0_B(IP_glob)=IP
          END IF
          IF (eCov == 1) THEN
            ListMapped1_B(IP_glob)=IP
          END IF
        END DO
        nbCommon_recv=0
        nbCase1=0
        nbCase2=0
        DO IP=1,NP_RESloc
          IP_glob=ListIPLG(IP+ListFirst(iProc))
          eCov=LocalColor % ListCovLower(IP+ListFirst(iProc))
          IF (eCov .eq. 1) THEN
            IF (ListMapped0(IP_glob) .gt. 0) THEN
              nbCommon_recv=nbCommon_recv+1
              nbCase1=nbCase1+1
            ELSE
              IPmap=ListMapped1(IP_glob)
              IF ((IPmap .gt. 0).and.(IPmap .gt. NP_RES)) THEN
                nbCommon_recv=nbCommon_recv+1
                nbCase2=nbCase2+1
              ENDIF
            END IF
          END IF
        END DO
# ifdef DEBUG
        WRITE(740+myrank,*) 'i=', iProc-1,  ' RnbCase12=', nbCase1, nbCase2
# endif
        IF (nbCommon_recv .gt. 0) THEN
          sync_nnbr_recv=sync_nnbr_recv+1
          ListCommon_recv(iProc)=nbCommon_recv
        END IF
        nbCommon_send=0
        DO IP=1,MNPloc
          IP_glob=ListIPLG(IP+ListFirst(iProc))
          eCov=LocalColor % ListCovLower(IP+ListFirst(iProc))
          IPmap=ListMapped1(IP_glob)
          IF ((IPmap .gt. 0).and.(IPmap .le. NP_RES)) THEN
            IF ((eCov .eq. 0).or.(IP.gt.NP_RESloc)) THEN
              nbCommon_send=nbCommon_send+1
            END IF
          END IF
        END DO
        IF (nbCommon_send .gt. 0) THEN
          sync_nnbr_send=sync_nnbr_send+1
          ListCommon_send(iProc)=nbCommon_send
        END IF
# ifdef DEBUG
        WRITE(740+myrank,*) '   nbCommon(send/recv)=', nbCommon_send, nbCommon_recv
# endif
      END DO
# ifdef DEBUG
      WRITE(740+myrank,*) 'WWM_P2D: sync_nnbr_send=', sync_nnbr_send
      WRITE(740+myrank,*) 'WWM_P2D: sync_nnbr_recv=', sync_nnbr_recv
# endif
      LocalColor % sync_nnbr_send=sync_nnbr_send
      LocalColor % sync_nnbr_recv=sync_nnbr_recv
      allocate(LocalColor % sync_ListNbCommon_send(sync_nnbr_send), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 91')
      allocate(LocalColor % sync_ListNbCommon_recv(sync_nnbr_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 92')
      allocate(LocalColor % sync_ListNeigh_send(sync_nnbr_send), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 93')
      allocate(LocalColor % sync_ListNeigh_recv(sync_nnbr_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 94')
      idx_send=0
      idx_recv=0
      DO iProc=1,nproc
        IF (ListCommon_send(iProc) .gt. 0) THEN
          idx_send=idx_send+1
          LocalColor % sync_ListNeigh_send(idx_send)=iProc
          nbCommon=ListCommon_send(iProc)
          LocalColor % sync_ListNbCommon_send(idx_send)=nbCommon
        END IF
        IF (ListCommon_recv(iProc) .gt. 0) THEN
          idx_recv=idx_recv+1
          LocalColor % sync_ListNeigh_recv(idx_recv)=iProc
          nbCommon=ListCommon_recv(iProc)
          LocalColor % sync_ListNbCommon_recv(idx_recv)=nbCommon
        END IF
      END DO
      !
      ! Now creating the sync exchange
      !
      allocate(LocalColor % sync_p2dsend_rqst(sync_nnbr_send), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 95')
      allocate(LocalColor % sync_p2drecv_rqst(sync_nnbr_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 96')
      allocate(LocalColor % sync_p2dsend_stat(MPI_STATUS_SIZE,sync_nnbr_send), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 97')
      allocate(LocalColor % sync_p2drecv_stat(MPI_STATUS_SIZE,sync_nnbr_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 98')
      allocate(LocalColor % sync_p2dsend_type(sync_nnbr_send), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 99')
      allocate(LocalColor % sync_p2drecv_type(sync_nnbr_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 100')
# ifdef DEBUG
      WRITE(740+myrank,*) 'SYNC sync_nnbr_send=', sync_nnbr_send
# endif
      DO iNeigh=1,sync_nnbr_send
        iProc=LocalColor % sync_ListNeigh_send(iNeigh)
        MNPloc=ListMNP(iProc)
        NP_RESloc=ListNP_RES(iProc)
        ListMapped0_B=0
        ListMapped1_B=0
        DO IP=1,MNPloc
          IP_glob=ListIPLG(IP+ListFirst(iProc))
          eCov=LocalColor % ListCovLower(IP+ListFirst(iProc))
          IF (eCov == 0) THEN
            ListMapped0_B(IP_glob)=IP
          END IF
          IF (eCov == 1) THEN
            ListMapped1_B(IP_glob)=IP
          END IF
        END DO
        nbCommon=LocalColor % sync_ListNbCommon_send(iNeigh)
# ifdef DEBUG
        WRITE(740+myrank,*) '   SYNC iNeigh=', iNeigh, ' nbCommon=', nbCommon
# endif
        allocate(dspl_send(nbCommon), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 101')
        idx=0
        DO IP=1,NP_RES
          IF (LocalColor % CovLower(IP) .eq. 1) THEN
            DoOper=0
            IP_glob=iplg(IP)
            IPmap=ListMapped0_B(IP_glob)
            IF (IPmap .gt. 0) THEN
              DoOper=1
            ELSE
              IPmap=ListMapped1_B(IP_glob)
              IF ((IPmap .gt. 0).and.(IPmap.gt.NP_RESloc)) THEN
                DoOper=1
              END IF
            END IF
            IF (DoOper == 1) THEN
              idx=idx+1
              dspl_send(idx)=MSC*MDC*(IP-1)
# ifdef DEBUG
              WRITE(740+myrank,*) 'idx=', idx, 'IP=', IP
              WRITE(740+myrank,*) '  IP_glob=', IP_glob, 'IPmap=', IPmap
# endif
            END IF
          END IF
        END DO
# ifdef DEBUG
        IF (idx .ne. nbCommon) THEN
          CALL WWM_ABORT('error in SYNC_p2dsend')
        END IF
        WRITE(740+myrank,*) '   SYNC_p2dsend iProc=', iProc-1, ' nbCommon=', nbCommon
# endif
        call mpi_type_create_indexed_block(nbCommon,MSC*MDC,dspl_send,rtype,LocalColor % sync_p2dsend_type(iNeigh), ierr)
        call mpi_type_commit(LocalColor % sync_p2dsend_type(iNeigh), ierr)
        deallocate(dspl_send)
      END DO
# ifdef DEBUG
      WRITE(740+myrank,*) 'SYNC sync_nnbr_recv=', sync_nnbr_recv
# endif
      DO iNeigh=1,sync_nnbr_recv
        iProc=LocalColor % sync_ListNeigh_recv(iNeigh)
        MNPloc=ListMNP(iProc)
        NP_RESloc=ListNP_RES(iProc)
        ListMapped0_B=0
        ListMapped1_B=0
        DO IP=1,MNPloc
          IP_glob=ListIPLG(IP+ListFirst(iProc))
          eCov=LocalColor % ListCovLower(IP+ListFirst(iProc))
          IF (eCov == 0) THEN
            ListMapped0_B(IP_glob)=IP
          END IF
          IF (eCov == 1) THEN
            ListMapped1_B(IP_glob)=IP
          END IF
        END DO
        nbCommon=LocalColor%sync_ListNbCommon_recv(iNeigh)
# ifdef DEBUG
        WRITE(740+myrank,*) '   SYNC iNeigh=', iNeigh, ' nbCommon=', nbCommon
# endif
        allocate(dspl_recv(nbCommon), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 102')
        idx=0
        DO IP=1,NP_RESloc
          IP_glob=ListIPLG(IP+ListFirst(iProc))
          eCov=LocalColor % ListCovLower(IP+ListFirst(iProc))
          IF (eCov == 1) THEN
            IPmap=ListMapped0(IP_glob)
            DoOper=0
            IF (IPmap .gt. 0) THEN
              DoOper=1
            ELSE
              IPmap=ListMapped1(IP_glob)
              IF ((IPmap .gt. 0).and.(IPmap.gt.NP_RES)) THEN
                DoOper=1
              END IF
            END IF
            IF (DoOper == 1) THEN
              idx=idx+1
              dspl_recv(idx)=MSC*MDC*(IPmap-1)
# ifdef DEBUG
              WRITE(740+myrank,*) 'idx=', idx, 'IP=', IP
              WRITE(740+myrank,*) '  IP_glob=', IP_glob, 'IPmap=', IPmap
# endif
            END IF
          END IF
        END DO
# ifdef DEBUG
        IF (idx .ne. nbCommon) THEN
          CALL WWM_ABORT('error in SYNC_p2drecv')
        END IF
        WRITE(740+myrank,*) '   SYNC_p2drecv iProc=', iProc-1, ' nbCommon=', nbCommon
# endif
# ifdef DEBUG
        WRITE(740+myrank,*) '   nbCase1=', nbCase1, ' nbCase2=', nbCase2
# endif
        call mpi_type_create_indexed_block(nbCommon,MSC*MDC,dspl_recv,rtype,LocalColor % sync_p2drecv_type(iNeigh),ierr)
        call mpi_type_commit(LocalColor % sync_p2drecv_type(iNeigh), ierr)
        deallocate(dspl_recv)
      END DO
# ifdef LU_SOLVE_RWRT
      lenMNP=0
      IF (LocalColor % nbNeedSend_blk > lenMNP) THEN
        lenMNP=LocalColor % nbNeedSend_blk
      END IF
      IF (LocalColor % nbNeedRecv_blk > lenMNP) THEN
        lenMNP=LocalColor % nbNeedRecv_blk
      END IF
      IF (LocalColor % nbNeedSend_u2l > lenMNP) THEN
        lenMNP=LocalColor % nbNeedSend_u2l
      END IF
      IF (LocalColor % nbNeedRecv_u2l > lenMNP) THEN
        lenMNP=LocalColor % nbNeedRecv_u2l
      END IF
#  ifdef DEBUG
      WRITE(7000+myrank,*) 'nbNeedSend_blk=', LocalColor % nbNeedSend_blk
      WRITE(7000+myrank,*) 'nbNeedRecv_blk=', LocalColor % nbNeedRecv_blk
      WRITE(7000+myrank,*) 'nbNeedSend_u2l=', LocalColor % nbNeedSend_u2l
      WRITE(7000+myrank,*) 'nbNeedRecv_u2l=', LocalColor % nbNeedRecv_u2l
      WRITE(7000+myrank,*) 'lenMNP=', lenMNP
#  endif
      allocate(LocalColor % ACexch(maxBlockLength, lenMNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 58')
# else
      allocate(LocalColor % ACexch(maxBlockLength, MNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 58')
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DETERMINE_JSTATUS_L_U(LocalColor)
      USE DATAPOOL, only : LocalColorInfo, MNP, rkind
      USE DATAPOOL, only : NNZ, IA, JA, NP_RES, I_DIAG
      USE DATAPOOL, only : wwm_nnbr_send, wwm_nnbr_recv
      USE DATAPOOL, only : wwm_p2drecv_type, wwm_p2dsend_type
      USE DATAPOOL, only : wwm_ListNeigh_send, wwm_ListNeigh_recv
      USE elfe_msgp, only : myrank, nproc, comm, ierr, nbrrank_p
      implicit none
      include 'mpif.h'
      type(LocalColorInfo), intent(inout) :: LocalColor
      integer Jstatus_L(NNZ), Jstatus_U(NNZ)
      integer IP, J, JP, DoOper
      integer istat
      integer nb, idx
      Jstatus_L=0
      Jstatus_U=0
      DO IP=1,NP_RES
        IF (LocalColor%CovLower(IP) == 1) THEN
          DO J=IA(IP),IA(IP+1)-1
            JP=JA(J)
            IF (JP /= IP) THEN
              IF (JP .lt. IP) THEN
                DoOper=1
              ELSE
                IF (LocalColor % CovLower(JP) == 0) THEN
                  DoOper=1
                ELSE
                  DoOper=0
                END IF
              END IF
            ELSE
              DoOper=0
            END IF
            Jstatus_L(J)=DoOper
          END DO
          DO J=IA(IP),IA(IP+1)-1
            JP=JA(J)
            IF (JP /= IP) THEN
              IF (JP .lt. IP) THEN
                DoOper=0
              ELSE
                IF (LocalColor % CovLower(JP) == 0) THEN
                  DoOper=0
                ELSE
                  DoOper=1
                END IF
              END IF
            ELSE
              DoOper=0
            END IF
            Jstatus_U(J)=DoOper
          END DO
        END IF
      END DO
# ifdef DEBUG
      DO IP=1,NP_RES
        DO J=IA(IP),IA(IP+1)-1
          JP=JA(J)
          IF ((Jstatus_L(J).eq.1).and.(Jstatus_U(J).eq.1)) THEN
            WRITE(myrank+919,*) 'MNP=', MNP, 'NP_RES=', NP_RES
            WRITE(myrank+919,*) 'IP=', IP, ' JP=', JP
            WRITE(myrank+919,*) 'IPcovLower=', LocalColor%CovLower(IP)
            WRITE(myrank+919,*) 'JPcovLower=', LocalColor%CovLower(JP)
            WRITE(myrank+919,*) 'We have major error'
            CALL WWM_ABORT('Please panic and debug')
          END IF
        END DO
      END DO
      WRITE(740+myrank,*) 'sum(Jstatus_L)=', sum(Jstatus_L)
      WRITE(740+myrank,*) 'sum(Jstatus_U)=', sum(Jstatus_U)
# endif
# ifdef REORDER_ASPAR_PC
      allocate(LocalColor % IA_L(NP_RES+1), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('determine_jstatus_L_U, error 1')
      allocate(LocalColor % IA_U(NP_RES+1), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('determine_jstatus_L_U, error 2')
      allocate(LocalColor % JA_LU(NNZ), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('determine_jstatus_L_U, error 3')
      allocate(LocalColor % JmapR(NNZ), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('determine_jstatus_L_U, error 5')
      LocalColor % JmapR=-1
      LocalColor % IA_L(1)=1
      idx=0
      DO IP=1,NP_RES
        nb=0
        DO J=IA(IP),IA(IP+1)-1
          IF (Jstatus_L(J).eq.1) THEN
            JP=JA(J)
            nb=nb+1
            idx=idx+1
            LocalColor % JmapR(J)=idx
            LocalColor % JA_LU(idx)=JP
          END IF
        END DO
        LocalColor % IA_L(IP+1)=LocalColor % IA_L(IP)+nb
      END DO
      LocalColor % IA_U(1)=LocalColor % IA_L(NP_RES+1)
      DO IP=1,NP_RES
        nb=0
        DO J=IA(IP),IA(IP+1)-1
          IF (Jstatus_U(J).eq.1) THEN
            JP=JA(J)
            nb=nb+1
            idx=idx+1
            LocalColor % JmapR(J)=idx
            LocalColor % JA_LU(idx)=JP
          END IF
        END DO
        nb=nb+1
        idx=idx+1
        J=I_DIAG(IP)
        LocalColor % JmapR(J)=idx
        LocalColor % JA_LU(idx)=JP
        LocalColor % IA_U(IP+1)=LocalColor % IA_U(IP)+nb
      END DO
# endif
      allocate(LocalColor % Jstatus_L(NNZ), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 103')
      allocate(LocalColor % Jstatus_U(NNZ), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 104')
      LocalColor % Jstatus_L=Jstatus_L
      LocalColor % Jstatus_U=Jstatus_U
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5_RECV_ASPAR_PC(LocalColor, SolDat)
      USE DATAPOOL, only : LocalColorInfo, I5_SolutionData
      USE DATAPOOL, only : MNP, MSC, MDC, rkind
      USE elfe_msgp, only : ierr, comm, rtype, istatus, nbrrank_p
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      type(I5_SolutionData), intent(inout) :: SolDat
      real(rkind), allocatable :: ASPAR_rs(:)
      integer idx, iNNZ, jNNZ, IS, ID, NNZ_l, siz
      integer iProc, i, iRank
      integer istat
      DO iProc=1,LocalColor % nbLow_send
        i=LocalColor % ListIdxUpper_send(iProc)
        iRank=nbrrank_p(i)
        NNZ_l=LocalColor % NNZ_len_r(iProc)
        siz=NNZ_l*MSC*MDC
        allocate(ASPAR_rs(NNZ_l*MSC*MDC), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 105')
        CALL mpi_recv(ASPAR_rs,siz,rtype,iRank,45,comm,istatus,ierr)
        idx=0
        DO iNNZ=1,NNZ_l
          jNNZ=LocalColor % NNZ_index_r(iProc,iNNZ)
          DO IS=1,MSC
            DO ID=1,MDC
              idx=idx+1
              SolDat % ASPAR_pc(jNNZ,IS,ID)=ASPAR_rs(idx)
            END DO
          END DO
        END DO
        deallocate(ASPAR_rs)
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5_SEND_ASPAR_PC(LocalColor, SolDat)
      USE DATAPOOL, only : LocalColorInfo, I5_SolutionData
      USE DATAPOOL, only : MSC, MDC, rkind
      USE elfe_msgp, only : ierr, comm, rtype, nbrrank_p
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      type(I5_SolutionData), intent(inout) :: SolDat
      real(rkind), allocatable :: ASPAR_rs(:)
      integer idx, iNNZ, jNNZ, IS, ID, NNZ_l, siz
      integer iProc, iRank, i
      integer istat
      DO iProc=1,LocalColor % nbUpp_send
        i=LocalColor % ListIdxUpper_send(iProc)
        iRank=nbrrank_p(i)
        NNZ_l=LocalColor % NNZ_len_s(iProc)
        siz=NNZ_l*MSC*MDC
        allocate(ASPAR_rs(NNZ_l*MSC*MDC), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 106')
        idx=0
        DO iNNZ=1,NNZ_l
          jNNZ=LocalColor % NNZ_index_s(iProc,iNNZ)
          DO IS=1,MSC
            DO ID=1,MDC
              idx=idx+1
              ASPAR_rs(idx)=SolDat % ASPAR_pc(jNNZ,IS,ID)
            END DO
          END DO
        END DO
        CALL mpi_isend(ASPAR_rs,siz,rtype,iRank,45,comm,LocalColor%Upp_s_rq(iProc),ierr)
        deallocate(ASPAR_rs)
      END DO
      IF (LocalColor % nbUpp_send > 0) THEN
        call mpi_waitall(LocalColor %nbUpp_send, LocalColor % Upp_s_rq, LocalColor % Upp_s_stat,ierr)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5_CREATE_PRECOND_ILU0(LocalColor, SolDat)
      USE DATAPOOL, only : MNP, NP_RES, LocalColorInfo, I5_SolutionData
      USE DATAPOOL, only : MSC, MDC, IA, JA, I_DIAG, rkind
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      type(I5_SolutionData), intent(inout) :: SolDat
      integer IP, JP, J, JP2, J2, J_FOUND, IS, ID
      integer, allocatable :: ListJ(:)
      integer istat
      real(rkind) tl
      SolDat%ASPAR_pc=SolDat%ASPAR_block
      CALL I5_RECV_ASPAR_PC(LocalColor, SolDat)
      allocate(ListJ(MNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 107')
      DO IP=1,NP_RES
        IF (LocalColor % CovLower(IP) == 1) THEN
          DO J=IA(IP),IA(IP+1)-1
            JP=JA(J)
            ListJ(JP)=J
          END DO
          DO J=IA(IP),I_DIAG(IP)-1
            JP=JA(J)
            DO IS=1,MSC
              DO ID=1,MDC
                tl=SolDat%ASPAR_pc(J,IS,ID)*SolDat%ASPAR_pc(I_DIAG(JP),IS,ID)
                DO J2=IA(JP),IA(JP+1)-1
                  JP2=JA(J2)
                  J_FOUND=ListJ(JP2)
                  IF (J_FOUND.gt.0) THEN ! Here is ILU0 approximation
                    SolDat%ASPAR_pc(J_FOUND,IS,ID)=SolDat%ASPAR_pc(J_FOUND,IS,ID) - tl*SolDat%ASPAR_pc(J2,IS,ID)
                  END IF
                END DO
                SolDat%ASPAR_pc(J,IS,ID)=tl
              END DO
            END DO
          END DO
          DO J=IA(IP),IA(IP+1)-1
            JP=JA(J)
            ListJ(JP)=0
          END DO
          J=I_DIAG(IP)
          DO IS=1,MSC
            DO ID=1,MDC
              SolDat%ASPAR_pc(J,IS,ID)=1.0_rkind/SolDat%ASPAR_pc(J,IS,ID)
            END DO
          END DO
        END IF
      END DO
      deallocate(ListJ)
      CALL I5_SEND_ASPAR_PC(LocalColor, SolDat)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_CREATE_PRECOND_ILU0(LocalColor, SolDat)
      USE DATAPOOL, only : MNP, NP_RES, LocalColorInfo, I5_SolutionData
      USE DATAPOOL, only : MSC, MDC, IA, JA, I_DIAG, rkind
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      type(I5_SolutionData), intent(inout) :: SolDat
      integer IP, JP, J, JP2, J2, J_FOUND, IS, ID
      integer, allocatable :: ListJ(:)
      integer istat
      real(rkind) tl
      CALL WWM_ABORT('Please program it')
      END SUBROUTINE
!**********************************************************************
!* We assign the values only for CovLower(IP)=1                       *
!* We could with some effort assign values for all with some effort   *
!* but the values would not be used                                   *
!**********************************************************************
      SUBROUTINE I5B_CREATE_PRECOND_SOR(LocalColor, SolDat)
      USE DATAPOOL, only : MNP, MDC, IA, JA, I_DIAG, NP_RES, rkind, ONE
      USE DATAPOOL, only : LocalColorInfo, I5_SolutionData
      USE elfe_msgp, only : exchange_p4d_wwm
# ifdef DEBUG
      USE elfe_msgp, only : myrank
# endif
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      type(I5_SolutionData), intent(inout) :: SolDat
      integer IP, ID, IS, JP, JR, J1, J, IPglob, JPglob
      real(rkind) eVal
# if defined DEBUG
      real(rkind) :: eSum, eSumB
# endif
      DO IP=1,NP_RES
        J=I_DIAG(IP)
        SolDat%AC4(:,:,IP)=ONE/SolDat % ASPAR_block(:,:,J)
      END DO
# ifdef SELFE_EXCHANGE
      CALL I5B_EXCHANGE_P4D_WWM(LocalColor, SolDat%AC4)
# else
      CALL EXCHANGE_P4D_WWM(SolDat%AC4)
# endif
      DO IP=1,NP_RES
        IF (LocalColor%CovLower(IP) == 1) THEN
# if defined REORDER_ASPAR_PC
          DO J=IA(IP),IA(IP+1)-1
            IF (LocalColor%Jstatus_L(J) == 1) THEN
              JP=JA(J)
              JR=LocalColor%JmapR(J)
              SolDat % ASPAR_pc(:,:,JR)=SolDat % ASPAR_block(:,:,J)*SolDat%AC4(:,:,JP)
            END IF
          ENDDO
          J=LocalColor% IA_U(IP+1)-1
          SolDat % ASPAR_pc(:,:,J)=SolDat%AC4(:,:,IP)
          DO J=IA(IP),IA(IP+1)-1
            IF (LocalColor%Jstatus_U(J) == 1) THEN
              JP=JA(J)
              JR=LocalColor%JmapR(J)
              SolDat % ASPAR_pc(:,:,JR)=SolDat % ASPAR_block(:,:,J)
            END IF
          END DO
# else
          DO J=IA(IP),IA(IP+1)-1
            IF (LocalColor%Jstatus_L(J) == 1) THEN
              JP=JA(J)
              SolDat % ASPAR_pc(:,:,J)=SolDat % ASPAR_block(:,:,J)*SolDat%AC4(:,:,JP)
            END IF
          ENDDO
          J=I_DIAG(IP)
          SolDat % ASPAR_pc(:,:,J)=SolDat%AC4(:,:,IP)
          DO J=IA(IP),IA(IP+1)-1
            IF (LocalColor%Jstatus_U(J) == 1) THEN
              JP=JA(J)
              SolDat % ASPAR_pc(:,:,J)=SolDat % ASPAR_block(:,:,J)
            END IF
          END DO
# endif
        END IF
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_CREATE_PRECOND(LocalColor, SolDat, TheMethod)
      USE DATAPOOL, only : LocalColorInfo, I5_SolutionData
# ifdef DEBUG
      USE elfe_msgp, only : myrank
# endif
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      type(I5_SolutionData), intent(inout) :: SolDat
      integer, intent(in) :: TheMethod
      IF (TheMethod == 1) THEN ! SOR 
        CALL I5B_CREATE_PRECOND_SOR(LocalColor, SolDat)
      ELSE IF (TheMethod == 2) THEN ! ILU0
        CALL I5B_CREATE_PRECOND_ILU0(LocalColor, SolDat)
      ELSE
        CALL WWM_ABORT('Wrong choice of preconditioner')
      ENDIF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_EXCHANGE_P3_LOW_2_UPP_Send(LocalColor, AC, iBlock)
      USE DATAPOOL, only : LocalColorInfo, MNP, MSC, MDC, rkind
      USE DATAPOOL, only : wwm_ListNeigh_send, wwm_p2dsend_type
      USE elfe_msgp, only : comm, ierr, myrank
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      REAL(rkind), intent(in) :: AC(LocalColor%MSCeffect, MDC, MNP)
      INTEGER, intent(in) :: iBlock
      integer iUpp, i, iRank, idx, lenBlock, maxBlockLength, IS, ID, IP, nbUpp_send
# ifdef LU_SOLVE_RWRT
      integer idxIP
# endif
      lenBlock=LocalColor % BlockLength(iBlock)
      maxBlockLength=LocalColor % maxBlockLength
# ifdef LU_SOLVE_RWRT
      DO idxIP=1,LocalColor % nbNeedSend_blk
        IP = LocalColor % IdxSend_blk(idxIP)
        DO idx=1,lenBlock
          IS=LocalColor % ISindex(iBlock, idx)
          ID=LocalColor % IDindex(iBlock, idx)
          LocalColor % ACexch(idx,idxIP)=AC(IS,ID,IP)
        END DO
      END DO
# else
      DO IP=1,MNP
        DO idx=1,lenBlock
          IS=LocalColor % ISindex(iBlock, idx)
          ID=LocalColor % IDindex(iBlock, idx)
          LocalColor % ACexch(idx,IP)=AC(IS,ID,IP)
        END DO
      END DO
# endif
      nbUpp_send=LocalColor % nbUpp_send
# ifdef L2U_OPER
      DO iUpp=1,nbUpp_send
        iRank=LocalColor % l2u_ListNeigh_send(iUpp)
        CALL mpi_isend(LocalColor % ACexch, 1, LocalColor % l2u_p2dsend_type(iUpp), iRank, 7, comm, LocalColor%Upp_s_rq(iUpp), ierr)
      END DO
# else
      DO iUpp=1,nbUpp_send
        i=LocalColor % ListIdxUpper_send(iUpp)
        iRank=wwm_ListNeigh_send(i)-1
        CALL mpi_isend(LocalColor % ACexch, 1, LocalColor % blk_p2dsend_type(i), iRank, 7, comm, LocalColor%Upp_s_rq(iUpp), ierr)
      END DO
# endif
      IF (nbUpp_send > 0) THEN
        call mpi_waitall(nbUpp_send, LocalColor % Upp_s_rq, LocalColor % Upp_s_stat,ierr)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_EXCHANGE_P3_LOW_2_UPP_Recv(LocalColor, AC, iBlock)
      USE DATAPOOL, only : LocalColorInfo, MNP, MSC, MDC, rkind
      USE DATAPOOL, only : wwm_ListNeigh_recv, wwm_p2drecv_type
      USE elfe_msgp, only : comm, ierr, myrank
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      REAL(rkind), intent(inout) :: AC(LocalColor%MSCeffect, MDC, MNP)
      INTEGER, intent(in) :: iBlock
      integer iProc, i, iRank, idx, lenBlock, maxBlockLength, IS, ID, IP, nbLow_recv
# ifdef LU_SOLVE_RWRT
      integer idxIP
# endif
      lenBlock=LocalColor % BlockLength(iBlock)
      maxBlockLength=LocalColor % maxBlockLength
# ifndef LU_SOLVE_RWRT
#  ifndef L2U_OPER
      DO IP=1,MNP
        DO idx=1,lenBlock
          IS=LocalColor % ISindex(iBlock, idx)
          ID=LocalColor % IDindex(iBlock, idx)
          LocalColor % ACexch(idx,IP)=AC(IS,ID,IP)
        END DO
      END DO
#  endif
# else
      DO idxIP=1,LocalColor % nbNeedRecv_blk
        IP = LocalColor % IdxRecv_blk(idxIP)
        DO idx=1,lenBlock
          IS=LocalColor % ISindex(iBlock, idx)
          ID=LocalColor % IDindex(iBlock, idx)
          LocalColor % ACexch(idx,idxIP) = AC(IS,ID,IP)
        END DO
      END DO
# endif
      nbLow_recv=LocalColor % nbLow_recv
# ifdef L2U_OPER
      DO iProc=1,nbLow_recv
        iRank=LocalColor % l2u_ListNeigh_recv(iProc)
        call mpi_irecv(LocalColor % ACexch,1,LocalColor % l2u_p2drecv_type(iProc),iRank,7,comm,LocalColor % Low_r_rq(iProc),ierr)
      END DO
# else
      DO iProc=1,nbLow_recv
        i=LocalColor % ListIdxLower_recv(iProc)
        iRank=wwm_ListNeigh_recv(i)-1
        call mpi_irecv(LocalColor % ACexch,1,LocalColor % blk_p2drecv_type(i),iRank,7,comm,LocalColor % Low_r_rq(iProc),ierr)
      END DO
# endif
      IF (nbLow_recv > 0) THEN
        call mpi_waitall(nbLow_recv, LocalColor%Low_r_rq, LocalColor%Low_r_stat,ierr)
      END IF
# ifdef LU_SOLVE_RWRT
      DO idxIP=1,LocalColor % nbNeedRecv_blk
        IP = LocalColor % IdxRecv_blk(idxIP)
        DO idx=1,lenBlock
          IS=LocalColor % ISindex(iBlock, idx)
          ID=LocalColor % IDindex(iBlock, idx)
          AC(IS,ID,IP) = LocalColor % ACexch(idx,idxIP)
        END DO
      END DO
# else
      DO IP=1,MNP
        DO idx=1,lenBlock
          IS=LocalColor % ISindex(iBlock, idx)
          ID=LocalColor % IDindex(iBlock, idx)
          AC(IS,ID,IP)=LocalColor % ACexch(idx,IP)
        END DO
      END DO
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_EXCHANGE_P3_UPP_2_LOW_Send(LocalColor, AC, iBlock)
      USE DATAPOOL, only : LocalColorInfo, MNP, MSC, MDC, rkind
      USE DATAPOOL, only : wwm_ListNeigh_send, wwm_p2dsend_type
      USE elfe_msgp, only : comm, ierr, myrank
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      REAL(rkind), intent(in) :: AC(LocalColor%MSCeffect, MDC, MNP)
      INTEGER, intent(in) :: iBlock
      integer iProc, iRank, idx, lenBlock, maxBlockLength, IS, ID, nbLow_send, IP
# ifdef LU_SOLVE_RWRT
      integer idxIP
# endif
      lenBlock=LocalColor % BlockLength(iBlock)
# ifdef LU_SOLVE_RWRT
      DO idxIP=1,LocalColor % nbNeedSend_u2l
        IP = LocalColor % IdxSend_u2l(idxIP)
        DO idx=1,lenBlock
          IS=LocalColor % ISindex(iBlock, idx)
          ID=LocalColor % IDindex(iBlock, idx)
          LocalColor % ACexch(idx,idxIP)=AC(IS,ID,IP)
        END DO
      END DO
# else
      DO IP=1,MNP
        DO idx=1,lenBlock
          IS=LocalColor % ISindex(iBlock, idx)
          ID=LocalColor % IDindex(iBlock, idx)
          LocalColor % ACexch(idx,IP)=AC(IS,ID,IP)
        END DO
      END DO
# endif
      nbLow_send=LocalColor % u2l_nnbr_send
      DO iProc=1,nbLow_send
        iRank=LocalColor % u2l_ListNeigh_send(iProc)
        call mpi_isend(LocalColor % ACexch,1,LocalColor%u2l_p2dsend_type(iProc),iRank-1,1151,comm,LocalColor%u2l_p2dsend_rqst(iProc),ierr)
      END DO
      IF (nbLow_send > 0) THEN
        call mpi_waitall(nbLow_send, LocalColor%u2l_p2dsend_rqst, LocalColor%u2l_p2dsend_stat,ierr)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_EXCHANGE_P3_UPP_2_LOW_Recv(LocalColor, AC, iBlock)
      USE DATAPOOL, only : LocalColorInfo, MNP, MSC, MDC, rkind
      USE DATAPOOL, only : wwm_ListNeigh_recv, wwm_p2drecv_type
      USE elfe_msgp, only : comm, ierr, myrank
      implicit none
      include 'mpif.h'
      type(LocalColorInfo), intent(in) :: LocalColor
      REAL(rkind), intent(inout) :: AC(LocalColor % MSCeffect, MDC, MNP)
      INTEGER, intent(in) :: iBlock
      integer iProc, i, iRank, idx, lenBlock, maxBlockLength
      integer nbUpp_recv, IS, ID, IP
# ifdef LU_SOLVE_RWRT
      integer idxIP
# endif
      lenBlock=LocalColor % BlockLength(iBlock)
# ifndef LU_SOLVE_RWRT
      DO IP=1,MNP
        DO idx=1,lenBlock
          IS=LocalColor % ISindex(iBlock, idx)
          ID=LocalColor % IDindex(iBlock, idx)
          LocalColor % ACexch(idx,IP)=AC(IS,ID,IP)
        END DO
      END DO
# else
      DO idxIP=1,LocalColor % nbNeedRecv_u2l
        IP = LocalColor % IdxRecv_u2l(idxIP)
        DO idx=1,lenBlock
          IS=LocalColor % ISindex(iBlock, idx)
          ID=LocalColor % IDindex(iBlock, idx)
          LocalColor % ACexch(idx,idxIP) = AC(IS,ID,IP)
        END DO
      END DO
# endif
      nbUpp_recv=LocalColor % u2l_nnbr_recv
      DO iProc=1,nbUpp_recv
        iRank=LocalColor % u2l_ListNeigh_recv(iProc)
        call mpi_irecv(LocalColor % ACexch,1,LocalColor%u2l_p2drecv_type(iProc),iRank-1,1151,comm,LocalColor % u2l_p2drecv_rqst(iProc),ierr)
      END DO
      IF (nbUpp_recv > 0) THEN
        call mpi_waitall(nbUpp_recv, LocalColor%u2l_p2drecv_rqst, LocalColor%u2l_p2drecv_stat,ierr)
      END IF
# ifdef LU_SOLVE_RWRT
      DO idxIP=1,LocalColor % nbNeedRecv_u2l
        IP = LocalColor % IdxRecv_u2l(idxIP)
        DO idx=1,lenBlock
          IS=LocalColor % ISindex(iBlock, idx)
          ID=LocalColor % IDindex(iBlock, idx)
          AC(IS,ID,IP) = LocalColor % ACexch(idx,idxIP)
        END DO
      END DO
# else
      DO IP=1,MNP
        DO idx=1,lenBlock
          IS=LocalColor % ISindex(iBlock, idx)
          ID=LocalColor % IDindex(iBlock, idx)
          AC(IS,ID,IP) = LocalColor % ACexch(idx,IP)
        END DO
      END DO
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
# ifdef DEBUG
      SUBROUTINE COMPUTE_TOTAL_INDEX_SHIFT(TheRes)
      USE DATAPOOL, only : NP_RES, IA, JA
      implicit none
      integer, intent(out) :: TheRes
      integer :: IP, JP, J
      TheRes=0
      DO IP=1,NP_RES
        DO J=IA(IP),IA(IP+1)-1
          JP=JA(J)
          TheRes=TheRes+abs(IP-JP)
        END DO
      END DO
      END SUBROUTINE
# endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_PARTIAL_SOLVE_L(LocalColor, SolDat, iBlock, ACret)
      USE DATAPOOL, only : LocalColorInfo, I5_SolutionData
      USE DATAPOOL, only : IA, JA, MSC, MDC, MNP, rkind, NP_RES
      USE DATAPOOL, only : DO_SOLVE_L, DO_SOLVE_U
# ifdef DEBUG
      USE elfe_msgp, only : myrank
      USE elfe_glbl, only : iplg
# endif
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      type(I5_SolutionData), intent(in) :: SolDat
      integer, intent(in) :: iBlock
      real(rkind), intent(inout) :: ACret(LocalColor%MSCeffect,MDC,MNP)
      real(rkind) :: eCoeff
      integer IP, idx, ID, IS, J, IP_glob
      integer lenBlock, JP
# ifdef DEBUG
      real(rkind) :: eSum
# endif
      lenBlock=LocalColor % BlockLength(iBlock)
      DO IP=1,NP_RES
        IF (LocalColor % CovLower(IP) == 1) THEN
# if defined REORDER_ASPAR_PC
          DO J=LocalColor % IA_L(IP),LocalColor % IA_L(IP+1)-1
            JP=LocalColor % JA_LU(J)
            DO idx=1,lenBlock
              IS=LocalColor % ISindex(iBlock, idx)
              ID=LocalColor % IDindex(iBlock, idx)
              eCoeff=SolDat % ASPAR_pc(IS,ID,J)
#  ifdef DEBUG
              IF ((IS.eq.4).and.(ID.eq.4)) THEN
                WRITE(myrank+790,*) 'eCoeff, AC12=', eCoeff, ACret(IS,ID,IP), ACret(IS,ID,JP)
              END IF
#  endif
              IF (DO_SOLVE_L) THEN
                ACret(IS,ID,IP)=ACret(IS,ID,IP) - eCoeff*ACret(IS,ID,JP)
              END IF
            END DO
          END DO
# else
          DO J=IA(IP),IA(IP+1)-1
            IF (LocalColor % Jstatus_L(J) .eq. 1) THEN
              JP=JA(J)
              DO idx=1,lenBlock
                IS=LocalColor % ISindex(iBlock, idx)
                ID=LocalColor % IDindex(iBlock, idx)
                eCoeff=SolDat % ASPAR_pc(IS,ID,J)
#  ifdef DEBUG
                IF ((IS.eq.4).and.(ID.eq.4)) THEN
                  WRITE(myrank+1490,*) 'IP, JP, J=', IP, JP, J
                  WRITE(myrank+790,*) 'eCoeff, AC12=', eCoeff, ACret(IS,ID,IP), ACret(IS,ID,JP)
                END IF
#  endif
                IF (DO_SOLVE_L) THEN
                  ACret(IS,ID,IP)=ACret(IS,ID,IP) - eCoeff*ACret(IS,ID,JP)
                END IF
              END DO
            END IF
          END DO
# endif
        ENDIF
      END DO
# if defined DEBUG
      WRITE(3000+myrank,*) 'Partial_SOLVE_L Sums of ACret'
      DO IS=1,LocalColor%MSCeffect
        DO ID=1,MDC
          eSum=sum(ACret(IS,ID,:))
          WRITE(3000+myrank,*) 'IS=', IS, 'ID=', ID, 'eSum=', eSum
        END DO
      END DO
      WRITE(3000+myrank,*) 'End of sums'
      CALL FLUSH(3000+myrank)
      CALL FLUSH(790+myrank)
      CALL FLUSH(1490+myrank)
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_PARTIAL_SOLVE_U(LocalColor, SolDat, iBlock, ACret)
      USE DATAPOOL, only : LocalColorInfo, I5_SolutionData
      USE DATAPOOL, only : MNP, IA, JA, I_DIAG, MSC, MDC, rkind, NP_RES
      USE DATAPOOL, only : DO_SOLVE_L, DO_SOLVE_U
      USE elfe_msgp, only : myrank
      USE elfe_glbl, only : iplg
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      type(I5_SolutionData), intent(in) :: SolDat
      integer, intent(in) :: iBlock
      real(rkind), intent(inout) :: ACret(LocalColor%MSCeffect,MDC,MNP)
      real(rkind) :: eCoeff
      integer lenBlock, IP, JP, idx, J, IS, ID, IP_glob
      integer DoOper
      lenBlock=LocalColor % BlockLength(iBlock)
      DO IP=NP_RES,1,-1
        IF (LocalColor % CovLower(IP) == 1) THEN
# if defined REORDER_ASPAR_PC
          DO J=LocalColor % IA_U(IP),LocalColor % IA_U(IP+1)-2
            JP=LocalColor % JA_LU(J)
            DO idx=1,lenBlock
              IS=LocalColor % ISindex(iBlock, idx)
              ID=LocalColor % IDindex(iBlock, idx)
              eCoeff=SolDat % ASPAR_pc(IS,ID,J)
              IF (DO_SOLVE_U) THEN
                ACret(IS,ID,IP)=ACret(IS,ID,IP) - eCoeff*ACret(IS,ID,JP)
              END IF
            END DO
          END DO
          J=LocalColor % IA_U(IP+1)-1
          DO idx=1,lenBlock
            IS=LocalColor % ISindex(iBlock, idx)
            ID=LocalColor % IDindex(iBlock, idx)
            IF (DO_SOLVE_U) THEN
              ACret(IS,ID,IP)=ACret(IS,ID,IP)*SolDat % ASPAR_pc(IS,ID,J)
            END IF
          END DO
# else
          DO J=IA(IP),IA(IP+1)-1
            IF (LocalColor % Jstatus_U(J) .eq. 1) THEN
              JP=JA(J)
              DO idx=1,lenBlock
                IS=LocalColor % ISindex(iBlock, idx)
                ID=LocalColor % IDindex(iBlock, idx)
                eCoeff=SolDat % ASPAR_pc(IS,ID,J)
                IF (DO_SOLVE_U) THEN
                  ACret(IS,ID,IP)=ACret(IS,ID,IP) - eCoeff*ACret(IS,ID,JP)
                END IF
              END DO
            END IF
          END DO
          J=I_DIAG(IP)
          DO idx=1,lenBlock
            IS=LocalColor % ISindex(iBlock, idx)
            ID=LocalColor % IDindex(iBlock, idx)
            IF (DO_SOLVE_U) THEN
              ACret(IS,ID,IP)=ACret(IS,ID,IP)*SolDat % ASPAR_pc(IS,ID,J)
            END IF
          END DO
# endif
        ENDIF
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_SYNC_SENDRECV(LocalColor, AC)
      USE DATAPOOL, only : LocalColorInfo, MNP, MSC, MDC, rkind
      USE DATAPOOL, only : wwm_ListNeigh_send, wwm_ListNeigh_recv
      USE elfe_msgp, only : comm, ierr, myrank
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      real(rkind), intent(inout) :: AC(MSC, MDC, MNP)
      INTEGER :: IP, IS, ID, i, iRank
      integer iSync
      integer nbSync_send, nbSync_recv
      nbSync_recv=LocalColor%sync_nnbr_recv
      nbSync_send=LocalColor%sync_nnbr_send
      DO iSync=1,nbSync_send
        iRank=LocalColor % sync_ListNeigh_send(iSync)
        CALL mpi_isend(AC, 1, LocalColor%sync_p2dsend_type(iSync), iRank-1, 1009, comm, LocalColor%sync_p2dsend_rqst(iSync), ierr)
      END DO
      DO iSync=1,nbSync_recv
        iRank=LocalColor % sync_ListNeigh_recv(iSync)
        call mpi_irecv(AC,1,LocalColor%sync_p2drecv_type(iSync),iRank-1,1009,comm,LocalColor%sync_p2drecv_rqst(iSync),ierr)
      END DO
      IF (nbSync_send > 0) THEN
        call mpi_waitall(nbSync_send, LocalColor%sync_p2dsend_rqst, LocalColor%sync_p2dsend_stat,ierr)
      END IF
      IF (nbSync_recv > 0) THEN
        call mpi_waitall(nbSync_recv, LocalColor%sync_p2drecv_rqst, LocalColor%sync_p2drecv_stat,ierr)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_APPLY_PRECOND(LocalColor, SolDat, ACret)
      USE DATAPOOL, only : LocalColorInfo, I5_SolutionData
      USE DATAPOOL, only : MSC, MDC, MNP, rkind
      USE elfe_msgp, only : myrank
      USE DATAPOOL, only : DO_SYNC_UPP_2_LOW, DO_SYNC_LOW_2_UPP, DO_SYNC_FINAL
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      type(I5_SolutionData), intent(in) :: SolDat
      real(rkind), intent(inout) :: ACret(LocalColor%MSCeffect, MDC, MNP)
      integer iBlock, lenBlock, idx, IS, ID
      integer maxBlockLength
      maxBlockLength=LocalColor % maxBlockLength
# ifdef DEBUG
      WRITE(myrank+7000,*) 'Nblock=', LocalColor % Nblock 
# endif
      DO iBlock=1,LocalColor % Nblock
# ifdef DEBUG
        WRITE(myrank+7000,*) 'iBlock=', iBlock
        WRITE(myrank+7000,*) '1: sum(ACret)=', sum(ACret)
# endif
        IF (DO_SYNC_LOW_2_UPP) THEN
          CALL I5B_EXCHANGE_P3_LOW_2_UPP_Recv(LocalColor, ACret, iBlock)
        END IF
# ifdef DEBUG
        WRITE(myrank+7000,*) '2: sum(ACret)=', sum(ACret)
# endif
        CALL I5B_PARTIAL_SOLVE_L(LocalColor, SolDat, iBlock, ACret)
# ifdef DEBUG
        WRITE(myrank+7000,*) '3: sum(ACret)=', sum(ACret)
# endif
        IF (DO_SYNC_LOW_2_UPP) THEN
          CALL I5B_EXCHANGE_P3_LOW_2_UPP_Send(LocalColor, ACret, iBlock)
        END IF
# ifdef DEBUG
        WRITE(myrank+7000,*) '4: sum(ACret)=', sum(ACret)
# endif
      END DO
      DO iBlock=1,LocalColor%Nblock
# ifdef DEBUG
        WRITE(myrank+7000,*) 'iBlock=', iBlock
        WRITE(myrank+7000,*) '1: sum(ACret)=', sum(ACret)
# endif
        IF (DO_SYNC_UPP_2_LOW) THEN
          CALL I5B_EXCHANGE_P3_UPP_2_LOW_Recv(LocalColor, ACret, iBlock)
        END IF
# ifdef DEBUG
        WRITE(myrank+7000,*) '2: sum(ACret)=', sum(ACret)
# endif
        CALL I5B_PARTIAL_SOLVE_U(LocalColor, SolDat, iBlock, ACret)
# ifdef DEBUG
        WRITE(myrank+7000,*) '3: sum(ACret)=', sum(ACret)
# endif
        IF (DO_SYNC_UPP_2_LOW) THEN
          CALL I5B_EXCHANGE_P3_UPP_2_LOW_Send(LocalColor, ACret, iBlock)
        END IF
# ifdef DEBUG
        WRITE(myrank+7000,*) '4: sum(ACret)=', sum(ACret)
# endif
      END DO
      IF (DO_SYNC_FINAL) THEN
        CALL I5B_SYNC_SENDRECV(LocalColor, ACret)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_APPLY_FCT(MSCeffect, SolDat,  ACin, ACret)
      USE DATAPOOL, only : I5_SolutionData, IA, JA, NP_RES, MDC, MNP, rkind
      USE elfe_msgp, only : exchange_p4d_wwm
      implicit none
      integer IP, J, idx
      type(I5_SolutionData), intent(inout) :: SolDat
      integer, intent(in) :: MSCeffect
      REAL(rkind), intent(in) :: ACin(MSCeffect, MDC, MNP)
      REAL(rkind), intent(inout) :: ACret(MSCeffect, MDC, MNP)
      REAL(rkind) :: eSum(MSCeffect,MDC)
      DO IP=1,NP_RES
        eSum=0
        DO J=IA(IP),IA(IP+1)-1
          idx=JA(J)
          eSum=eSum + SolDat % ASPAR_block(:,:,J)*ACin(:,:,idx)
        END DO
        ACret(:,:,IP)=eSum
      END DO
# ifdef SELFE_EXCHANGE
      CALL I5B_EXCHANGE_P4D_WWM(LocalColor, ACret)
# else
      CALL EXCHANGE_P4D_WWM(ACret)
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE REPLACE_NAN_ZERO(LocalColor, LScal)
      USE DATAPOOL, only : rkind, MSC, MDC, LocalColorInfo
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      real(rkind), intent(inout) :: LScal(LocalColor%MSCeffect,MDC)
      integer IS, ID
      DO IS=1,LocalColor%MSCeffect
        DO ID=1,MDC
          IF (LScal(IS,ID) .ne. LScal(IS,ID)) THEN
            LScal(IS,ID)=0
          END IF
        END DO
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_L2_LINF(MSCeffect, ACw1, ACw2, Norm_L2, Norm_LINF)
      USE DATAPOOL, only : rkind, MNP, MDC
      USE DATAPOOL, only : nwild_loc_res, NP_RES
      USE elfe_msgp, only : myrank, comm, ierr, nproc, istatus, rtype
      implicit none
      integer, intent(in) :: MSCeffect
      real(rkind), intent(in) :: ACw1(MSCeffect, MDC, MNP)
      real(rkind), intent(in) :: ACw2(MSCeffect, MDC, MNP)
      real(rkind), intent(inout) :: Norm_L2(MSCeffect, MDC)
      real(rkind), intent(inout) :: Norm_LINF(MSCeffect, MDC)
      real(rkind) :: LScal(MSCeffect, MDC, 2)
      real(rkind) :: RScal(MSCeffect, MDC, 2)
      integer IP, iProc, IS, ID
      LScal=0
      DO IP=1,NP_RES
        LScal(:,:,1)=LScal(:,:,1) + nwild_loc_res(IP)*((ACw1(:,:,IP) - ACw2(:,:,IP))**2)
        DO IS=1,MSCeffect
          DO ID=1,MDC
            LScal(IS,ID,2)=max(LScal(IS,ID,2), abs(ACw1(IS,ID,IP) - ACw2(IS,ID,IP)))
          END DO
        END DO
      END DO
      IF (myrank == 0) THEN
        DO iProc=2,nproc
          CALL MPI_RECV(RScal,MSCeffect*MDC*2,rtype, iProc-1, 19, comm, istatus, ierr)
          LScal(:,:,1) = LScal(:,:,1) + RScal(:,:,1)
          DO IS=1,MSCeffect
            DO ID=1,MDC
              LScal(IS,ID,2)=max(LScal(IS,ID,2), RScal(IS,ID,2))
            END DO
          END DO
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(LScal,MSCeffect*MDC*2,rtype, iProc-1, 23, comm, ierr)
        END DO
      ELSE
        CALL MPI_SEND(LScal,MSCeffect*MDC*2,rtype, 0, 19, comm, ierr)
        CALL MPI_RECV(LScal,MSCeffect*MDC*2,rtype, 0, 23, comm, istatus, ierr)
      END IF
      Norm_L2=LScal(:,:,1)
      Norm_LINF=LScal(:,:,2)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_SCALAR(MSCeffect, ACw1, ACw2, LScal)
      USE DATAPOOL, only : rkind, MNP, MDC
      USE DATAPOOL, only : nwild_loc_res, NP_RES
      USE elfe_msgp, only : myrank, comm, ierr, nproc, istatus, rtype
      implicit none
      integer, intent(in) :: MSCeffect
      real(rkind), intent(in) :: ACw1(MSCeffect, MDC, MNP)
      real(rkind), intent(in) :: ACw2(MSCeffect, MDC, MNP)
      real(rkind), intent(inout) :: LScal(MSCeffect, MDC)
      real(rkind) :: RScal(MSCeffect, MDC)
      integer IP, iProc
      LScal=0
      DO IP=1,NP_RES
        LScal=LScal + nwild_loc_res(IP)*ACw1(:,:,IP)*ACw2(:,:,IP)
      END DO
      IF (myrank == 0) THEN
        DO iProc=2,nproc
          CALL MPI_RECV(RScal,MSCeffect*MDC,rtype, iProc-1, 19, comm, istatus, ierr)
          LScal = LScal + RScal
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(LScal,MSCeffect*MDC,rtype, iProc-1, 23, comm, ierr)
        END DO
      ELSE
        CALL MPI_SEND(LScal,MSCeffect*MDC,rtype, 0, 19, comm, ierr)
        CALL MPI_RECV(LScal,MSCeffect*MDC,rtype, 0, 23, comm, istatus, ierr)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_TOTAL_COHERENCY_ERROR(MSCeffect, ACw, Lerror)
      USE DATAPOOL, only : MNP, MDC, rkind
      USE DATAPOOL, only : ListIPLG, ListMNP
      USE elfe_msgp, only : istatus, ierr, comm, rtype, myrank, nproc
      USE elfe_glbl, only : iplg, np_global
      implicit none
      integer, intent(in) :: MSCeffect
      real(rkind), intent(in) :: ACw(MSCeffect, MDC, MNP)
      real(rkind), intent(out) :: Lerror
      real(rkind), allocatable :: ACtotal(:,:,:)
      real(rkind), allocatable :: ACloc(:,:,:)
      real(rkind), allocatable :: rbuf_real(:)
      integer, allocatable :: ListFirstMNP(:)
      integer, allocatable :: eStatus(:)
      integer IP, iProc, IPglob, IS, ID
      integer MNPloc
      integer istat
      IF (myrank == 0) THEN
        Lerror=0
        allocate(ListFirstMNP(nproc), eStatus(np_global), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 108')
        ListFirstMNP=0
        eStatus=0
        DO iProc=2,nproc
          ListFirstMNP(iProc)=ListFirstMNP(iProc-1) + ListMNP(iProc-1)
        END DO
        allocate(ACtotal(MSCeffect, MDC, np_global), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 110')
        DO IP=1,MNP
          IPglob=iplg(IP)
          ACtotal(:,:,IPglob)=ACw(:,:,IP)
          eStatus(IPglob)=1
        END DO
        DO iProc=2,nproc
          MNPloc=ListMNP(iProc)
          allocate(ACloc(MSCeffect, MDC, MNPloc), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 111')
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
        deallocate(ListFirstMNP)
        deallocate(ACtotal)
        deallocate(eStatus)
        allocate(rbuf_real(1), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 112')
        rbuf_real(1)=Lerror
        DO iProc=2,nproc
          CALL MPI_SEND(rbuf_real,1,rtype, iProc-1, 23, comm, ierr)
        END DO
        deallocate(rbuf_real)
      ELSE
        CALL MPI_SEND(ACw,MNP*MSCeffect*MDC,rtype, 0, 53, comm, ierr)
        allocate(rbuf_real(1), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 113')
        CALL MPI_RECV(rbuf_real,1,rtype, 0, 23, comm, istatus, ierr)
        Lerror=rbuf_real(1)
        deallocate(rbuf_real)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_SUM_MAX(MSCeffect, ACw, LSum, LMax)
      USE DATAPOOL, only : rkind, MNP, MDC
      USE DATAPOOL, only : nwild_loc_res, NP_RES
      USE elfe_msgp, only : myrank, comm, ierr, nproc, istatus, rtype
      implicit none
      integer, intent(in) :: MSCeffect
      real(rkind), intent(in) :: ACw(MSCeffect, MDC, MNP)
      real(rkind), intent(inout) :: LSum(MSCeffect, MDC)
      real(rkind), intent(inout) :: LMax(MSCeffect, MDC)
      real(rkind) :: RScal(MSCeffect, MDC)
      integer IP, iProc, IS, ID
      LSum=0
      DO IP=1,NP_RES
        LSum=LSum + nwild_loc_res(IP)*ACw(:,:, IP)
      END DO
      DO IS=1,MSCeffect
        DO ID=1,MDC
          LMax(IS,ID)=maxval(ACw(IS,ID,:))
        END DO
      END DO
      IF (myrank == 0) THEN
        DO iProc=2,nproc
          CALL MPI_RECV(RScal,MSCeffect*MDC,rtype, iProc-1, 53, comm, istatus, ierr)
          LSum = LSum + RScal
        END DO
        DO iProc=2,nproc
          CALL MPI_RECV(RScal,MSCeffect*MDC,rtype, iProc-1, 59, comm, istatus, ierr)
          DO IS=1,MSCeffect
            DO ID=1,MDC
              IF (RScal(IS,ID) .gt. LMax(IS,ID)) THEN
                LMax(IS,ID)=RScal(IS,ID)
              END IF
            END DO
          END DO
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(LSum,MSCeffect*MDC,rtype, iProc-1, 197, comm, ierr)
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(LMax,MSCeffect*MDC,rtype, iProc-1, 199, comm, ierr)
        END DO
      ELSE
        CALL MPI_SEND(LSum,MSCeffect*MDC,rtype, 0, 53, comm, ierr)
        CALL MPI_SEND(LMax,MSCeffect*MDC,rtype, 0, 59, comm, ierr)
        CALL MPI_RECV(LSum,MSCeffect*MDC,rtype, 0, 197, comm, istatus, ierr)
        CALL MPI_RECV(LMax,MSCeffect*MDC,rtype, 0, 199, comm, istatus, ierr)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      FUNCTION I5B_SUMTOT(MSCeffect, ACw)
      USE DATAPOOL, only : rkind, MNP, MDC
      implicit none
      integer, intent(in) :: MSCeffect
      real(rkind), intent(in) :: ACw(MSCeffect, MDC, MNP)
      real(rkind) :: LSum(MSCeffect, MDC)
      real(rkind) :: LMax(MSCeffect, MDC)
      real(rkind) :: I5B_SUMTOT
      CALL I5B_SUM_MAX(MSCeffect, ACw, LSum, LMax)
      I5B_SUMTOT=sum(LSum)
      END FUNCTION
!**********************************************************************
!*                                                                    *
!**********************************************************************
! We use the notations of
! http://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method
! and we use K1=Id
! In this algorithm, the use of v_{i-1}, v_i can be replace to just "v"
! The same for x, r
! 
      SUBROUTINE I5B_BCGS_SOLVER(LocalColor, SolDat)
      USE DATAPOOL, only : MDC, MNP, NP_RES, NNZ, AC2, SOLVERTHR
      USE DATAPOOL, only : LocalColorInfo, I5_SolutionData, rkind
      USE DATAPOOL, only : PCmethod, STAT
# ifdef DEBUG
      USE elfe_msgp, only : myrank
# endif
      implicit none
      type(LocalColorInfo), intent(inout) :: LocalColor
      type(I5_SolutionData), intent(inout) :: SolDat
      REAL(rkind) :: Rho(LocalColor%MSCeffect,MDC)
      REAL(rkind) :: Prov(LocalColor%MSCeffect,MDC)
      REAL(rkind) :: Alpha(LocalColor%MSCeffect,MDC)
      REAL(rkind) :: Beta(LocalColor%MSCeffect,MDC)
      REAL(rkind) :: Omega(LocalColor%MSCeffect,MDC)
      REAL(rkind) :: MaxError, CritVal
      REAL(rkind) :: eSum1, eSum2
      REAL(rkind) :: TheTol
# ifdef DEBUG
      integer IS1, IS2
      REAL(rkind) :: Lerror
# endif
# ifdef FAST_NORM
      real(rkind) :: Norm_L2(LocalColor%MSCeffect,MDC), Norm_LINF(LocalColor%MSCeffect,MDC)
# endif
      integer :: MaxIter = 30
      integer IP, IS, ID, nbIter, MSCeffect
      MaxError=SOLVERTHR
      MSCeffect=LocalColor % MSCeffect
      CALL I5B_APPLY_FCT(MSCeffect, SolDat,  SolDat % AC2, SolDat % AC3)
      SolDat % AC1=0                               ! y
      SolDat % AC3=SolDat % B_block - SolDat % AC3 ! r residual
      SolDat % AC4=SolDat % AC3                    ! hat{r_0} term
      SolDat % AC5=0                               ! v
      SolDat % AC6=0                               ! p
      SolDat % AC7=0                               ! s
      SolDat % AC8=0                               ! z
      SolDat % AC9=0                               ! t
# ifdef DEBUG
      write(2000+myrank,*) 'Before loop'
      write(2000+myrank,*) 'sumtot(AC1)=', I5B_SUMTOT(MSCeffect, SolDat%AC1)
      write(2000+myrank,*) 'sumtot(AC2)=', I5B_SUMTOT(MSCeffect, SolDat%AC2)
      write(2000+myrank,*) 'sumtot(AC3)=', I5B_SUMTOT(MSCeffect, SolDat%AC3)
      write(2000+myrank,*) 'sumtot(AC4)=', I5B_SUMTOT(MSCeffect, SolDat%AC4)
      write(2000+myrank,*) 'sumtot(AC5)=', I5B_SUMTOT(MSCeffect, SolDat%AC5)
      write(2000+myrank,*) 'sumtot(AC6)=', I5B_SUMTOT(MSCeffect, SolDat%AC6)
      write(2000+myrank,*) 'sumtot(AC7)=', I5B_SUMTOT(MSCeffect, SolDat%AC7)
      write(2000+myrank,*) 'sumtot(AC8)=', I5B_SUMTOT(MSCeffect, SolDat%AC8)
      write(2000+myrank,*) 'sumtot(AC9)=', I5B_SUMTOT(MSCeffect, SolDat%AC9)
      CALL FLUSH(2000+myrank)
      DO IS=1,LocalColor%MSCeffect
        WRITE(myrank+240,*) 'IS, sum(AC3)=', IS, sum(SolDat%AC3(IS,:,:))
      END DO
# endif
      Rho=1
      Alpha=1
      Omega=1
      nbIter=0
      DO
        nbIter=nbIter+1

        ! L1: Rhoi =(\hat{r}_0, r_{i-1}
        CALL I5B_SCALAR(MSCeffect, SolDat % AC4, SolDat % AC3, Prov)
# ifdef DEBUG
        CALL I5B_TOTAL_COHERENCY_ERROR(MSCeffect, SolDat%AC4, Lerror)
        WRITE(myrank+240,*) 'error(AC4)=', Lerror
        CALL I5B_TOTAL_COHERENCY_ERROR(MSCeffect, SolDat%AC3, Lerror)
        WRITE(myrank+240,*) 'error(AC3)=', Lerror
        WRITE(myrank+240,*) 'sum(abs(AC4))=', sum(abs(SolDat%AC4))
        WRITE(myrank+240,*) 'sum(abs(AC3))=', sum(abs(SolDat%AC3))
        DO IS=1,LocalColor%MSCeffect
          WRITE(myrank+240,*) 'IS, sum(AC4)=', IS, sum(SolDat%AC4(IS,:,:))
        END DO
        DO IS=1,LocalColor%MSCeffect
          WRITE(myrank+240,*) 'IS, sum(Prov)=', IS, sum(Prov(IS,:))
        END DO
        DO IS=1,LocalColor%MSCeffect
          WRITE(myrank+240,*) 'IS, sum(AC4*AC3)=', IS, sum(SolDat%AC4(IS,:,:)*SolDat%AC2(IS,:,:))
        END DO
# endif

        ! L2: Beta=(RhoI/Rho(I-1))  *  (Alpha/Omega(i-1))
        Beta=(Prov/Rho)*(Alpha/Omega)
        CALL REPLACE_NAN_ZERO(LocalColor, Beta)
        Rho=Prov
# ifdef DEBUG
        DO IS=1,LocalColor%MSCeffect
          WRITE(myrank+240,*) 'IS, sum(Rho)=', IS, sum(Rho(IS,:))
        END DO
# endif

        ! L3: Pi = r(i-1) + Beta*(p(i-1) -omega(i-1)*v(i-1))
        DO IP=1,MNP
          SolDat%AC6(:,:,IP)=SolDat%AC3(:,:,IP)                        &
     &      + Beta(:,:)*SolDat%AC6(:,:,IP)                            &
     &      - Beta(:,:)*Omega(:,:)*SolDat%AC5(:,:,IP)
        END DO
# ifdef DEBUG
        write(2000+myrank,*) 'nbIter=', nbIter
        write(2000+myrank,*) 'sumtot(AC6)=', I5B_SUMTOT(MSCeffect, SolDat%AC6)
        CALL I5B_TOTAL_COHERENCY_ERROR(MSCeffect, SolDat%AC6, Lerror)
        WRITE(myrank+240,*) 'error(AC6)=', Lerror
        DO IS=1,LocalColor%MSCeffect
          WRITE(myrank+240,*) 'IS, sum(AC6)=', IS, sum(SolDat%AC6(IS,:,:))
        END DO
# endif

        ! L4 y=K^(-1) Pi
        SolDat%AC1=SolDat%AC6
        IF (PCmethod .gt. 0) THEN
          CALL I5B_APPLY_PRECOND(LocalColor, SolDat, SolDat%AC1)
        ENDIF
# ifdef DEBUG
        write(2000+myrank,*) 'nbIter=', nbIter
        write(2000+myrank,*) 'sumtot(AC1)=', I5B_SUMTOT(MSCeffect, SolDat%AC1)
        CALL I5B_TOTAL_COHERENCY_ERROR(MSCeffect, SolDat%AC1, Lerror)
        WRITE(myrank+240,*) 'error(AC1)=', Lerror
        DO IS=1,LocalColor%MSCeffect
          WRITE(myrank+240,*) 'IS, sum(AC1)=', IS, sum(SolDat%AC1(IS,:,:))
        END DO
# endif

        ! L5 vi=Ay
        CALL I5B_APPLY_FCT(MSCeffect, SolDat,  SolDat%AC1, SolDat%AC5)
# ifdef DEBUG
        write(2000+myrank,*) 'nbIter=', nbIter
        write(2000+myrank,*) 'sumtot(AC5)=', I5B_SUMTOT(MSCeffect, SolDat%AC5)
        CALL I5B_TOTAL_COHERENCY_ERROR(MSCeffect, SolDat%AC5, Lerror)
        WRITE(myrank+240,*) 'error(AC5)=', Lerror
        DO IS=1,LocalColor%MSCeffect
          WRITE(myrank+240,*) 'IS, sum(AC5)=', IS, sum(SolDat%AC5(IS,:,:))
        END DO
# endif

        ! L6 Alpha=Rho/(hat(r)_0, v_i)
        CALL I5B_SCALAR(MSCeffect, SolDat % AC4, SolDat % AC5, Prov)
# ifdef DEBUG
        DO IS=1,LocalColor%MSCeffect
          WRITE(myrank+240,*) 'IS, sum(Prov)=', IS, sum(Prov(IS,:))
        END DO
# endif
        Alpha(:,:)=Rho(:,:)/Prov(:,:)
        CALL REPLACE_NAN_ZERO(LocalColor, Alpha)
# ifdef DEBUG
        write(2000+myrank,*) 'nbIter=', nbIter
        write(2000+myrank,*) 'sumtot(Alpha)=', sum(Alpha)
# endif
        ! L7 s=r(i-1) - alpha v(i)
        DO IP=1,MNP
          SolDat%AC7(:,:,IP)=SolDat%AC3(:,:,IP)                        &
     &      - Alpha(:,:)*SolDat%AC5(:,:,IP)
        END DO
# ifdef DEBUG
        write(2000+myrank,*) 'nbIter=', nbIter
        write(2000+myrank,*) 'sum(AC7)=', sum(SolDat%AC7)
        write(2000+myrank,*) 'sumtot(AC7)=', I5B_SUMTOT(MSCeffect, SolDat%AC7)
        CALL I5B_TOTAL_COHERENCY_ERROR(MSCeffect, SolDat%AC7, Lerror)
        WRITE(myrank+240,*) 'error(AC7)=', Lerror
        DO IS=1,LocalColor%MSCeffect
          WRITE(myrank+240,*) 'IS, sum(AC7)=', IS, sum(SolDat%AC7(IS,:,:))
        END DO
        DO IS=1,LocalColor%MSCeffect
          WRITE(myrank+240,*) 'IS, sum(Alpha)=', IS, sum(Alpha(IS,:))
        END DO
# endif

        ! L8 z=K^(-1) s
        SolDat%AC8=SolDat%AC7
        IF (PCmethod .gt. 0) THEN
          CALL I5B_APPLY_PRECOND(LocalColor, SolDat, SolDat%AC8)
        END IF
# ifdef DEBUG
        write(2000+myrank,*) 'nbIter=', nbIter
        write(2000+myrank,*) 'sumtot(AC8)=', I5B_SUMTOT(MSCeffect, SolDat%AC8)
        CALL I5B_TOTAL_COHERENCY_ERROR(MSCeffect, SolDat%AC8, Lerror)
        WRITE(myrank+240,*) 'error(AC8)=', Lerror
        DO IS=1,LocalColor%MSCeffect
          WRITE(myrank+240,*) 'IS, sum(AC8)=', IS, sum(SolDat%AC8(IS,:,:))
        END DO
# endif

        ! L9 t=Az
        CALL I5B_APPLY_FCT(MSCeffect, SolDat,  SolDat%AC8, SolDat%AC9)
# ifdef DEBUG
        write(2000+myrank,*) 'nbIter=', nbIter
        write(2000+myrank,*) 'sumtot(AC9)=', I5B_SUMTOT(MSCeffect, SolDat%AC9)
        CALL I5B_TOTAL_COHERENCY_ERROR(MSCeffect, SolDat%AC9, Lerror)
        WRITE(myrank+240,*) 'error(AC9)=', Lerror
        DO IS=1,LocalColor%MSCeffect
          WRITE(myrank+240,*) 'IS, sum(AC9)=', IS, sum(SolDat%AC9(IS,:,:))
        END DO
# endif

        ! L10 omega=(t,s)/(t,t)
        CALL I5B_SCALAR(MSCeffect, SolDat % AC9, SolDat % AC7, Omega)
        CALL I5B_SCALAR(MSCeffect, SolDat % AC9, SolDat % AC9, Prov)
        Omega(:,:)=Omega(:,:)/Prov(:,:)
        CALL REPLACE_NAN_ZERO(LocalColor, Omega)

        ! L11 x(i)=x(i-1) + Alpha y + Omega z
        DO IP=1,MNP
          SolDat%AC2(:,:,IP)=SolDat%AC2(:,:,IP)                        &
     &      + Alpha(:,:)*SolDat%AC1(:,:,IP)                            &
     &      + Omega(:,:)*SolDat%AC8(:,:,IP)
        END DO
# ifdef DEBUG
        write(2000+myrank,*) 'nbIter=', nbIter
        write(2000+myrank,*) 'sumtot(AC2)=', I5B_SUMTOT(MSCeffect, SolDat%AC2)
        CALL I5B_TOTAL_COHERENCY_ERROR(MSCeffect, SolDat%AC2, Lerror)
        WRITE(myrank+240,*) 'error(AC2)=', Lerror
        DO IS=1,LocalColor%MSCeffect
          WRITE(myrank+240,*) 'IS, sum(AC2)=', IS, sum(SolDat%AC2(IS,:,:))
        END DO
# endif

        ! L12 If x is accurate enough finish
        CALL I5B_APPLY_FCT(MSCeffect, SolDat,  SolDat%AC2, SolDat%AC1)
# ifdef DEBUG
        write(2000+myrank,*) 'nbIter=', nbIter
        write(2000+myrank,*) 'sumtot(AC1)=', I5B_SUMTOT(MSCeffect, SolDat%AC1)
        CALL I5B_TOTAL_COHERENCY_ERROR(MSCeffect, SolDat%AC1, Lerror)
        WRITE(myrank+240,*) 'error(AC1)=', Lerror
        DO IS=1,LocalColor%MSCeffect
          WRITE(myrank+240,*) 'IS, sum(AC1)=', IS, sum(SolDat%AC1(IS,:,:))
        END DO
# endif
# if defined FAST_NORM
        CALL I5B_L2_LINF(MSCeffect, SolDat%AC1, SolDat%B_block, Norm_L2, Norm_LINF)
        CritVal=maxval(Norm_L2)
# else
        SolDat%AC8=SolDat%AC1 - SolDat%B_block
        CALL I5B_SCALAR(MSCeffect, SolDat % AC8, SolDat % AC8, Prov)
        CritVal=maxval(Prov)
# endif
        IF (CritVal .lt. MaxError) THEN
          EXIT
        ENDIF
        IF (nbIter .gt. MaxIter) THEN
          EXIT
        ENDIF

        ! L13 r=s-omega t
        DO IP=1,MNP
          SolDat%AC3(:,:,IP)=SolDat%AC7(:,:,IP)                        &
     &      - Omega(:,:)*SolDat%AC9(:,:,IP)
        END DO
# ifdef DEBUG
        write(2000+myrank,*) 'nbIter=', nbIter
        write(2000+myrank,*) 'sumtot(AC3)=', I5B_SUMTOT(MSCeffect, SolDat%AC3)
        CALL I5B_TOTAL_COHERENCY_ERROR(MSCeffect, SolDat%AC3, Lerror)
        WRITE(myrank+240,*) 'error(AC3)=', Lerror
        DO IS=1,LocalColor%MSCeffect
          WRITE(myrank+240,*) 'IS, sum(AC3)=', IS, sum(SolDat%AC3(IS,:,:))
        END DO
# endif
      END DO
      WRITE(STAT%FHNDL, *) 'nbIter=', nbIter
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
!
! With another node ordering, maybe better performance
!
      SUBROUTINE I5B_BCGS_REORG_SOLVER(LocalColor, SolDat)
      USE DATAPOOL, only : MDC, MNP, NP_RES, NNZ, AC2, SOLVERTHR
      USE DATAPOOL, only : LocalColorInfo, I5_SolutionData, rkind
      USE DATAPOOL, only : PCmethod, STAT
# ifdef DEBUG
      USE elfe_msgp, only : myrank
# endif
      implicit none
      type(LocalColorInfo), intent(inout) :: LocalColor
      type(I5_SolutionData), intent(inout) :: SolDat
      REAL(rkind) :: Rho(LocalColor%MSCeffect,MDC)
      REAL(rkind) :: Prov(LocalColor%MSCeffect,MDC)
      REAL(rkind) :: Alpha(LocalColor%MSCeffect,MDC)
      REAL(rkind) :: Beta(LocalColor%MSCeffect,MDC)
      REAL(rkind) :: Omega(LocalColor%MSCeffect,MDC)
      REAL(rkind) :: MaxError, CritVal
      REAL(rkind) :: eSum1, eSum2
      REAL(rkind) :: TheTol
      real(rkind) :: Norm_L2(LocalColor%MSCeffect,MDC), Norm_LINF(LocalColor%MSCeffect,MDC)
      integer :: MaxIter = 30
      integer IP, IS, ID, nbIter, MSCeffect
      MaxError=SOLVERTHR
      MSCeffect=LocalColor % MSCeffect
      CALL I5B_APPLY_FCT(MSCeffect, SolDat,  SolDat % AC2, SolDat % AC3)
      SolDat % AC1=0                               ! y
      SolDat % AC3=SolDat % B_block - SolDat % AC3 ! r residual
      SolDat % AC4=SolDat % AC3                    ! hat{r_0} term
      SolDat % AC5=0                               ! v
      SolDat % AC6=0                               ! p
      SolDat % AC7=0                               ! t
      Rho=1
      Alpha=1
      Omega=1
      nbIter=0
      DO
        nbIter=nbIter+1

        ! L1: Rhoi =(\hat{r}_0, r_{i-1}
        CALL I5B_SCALAR(MSCeffect, SolDat % AC4, SolDat % AC3, Prov)

        ! L2: Beta=(RhoI/Rho(I-1))  *  (Alpha/Omega(i-1))
        Beta=(Prov/Rho)*(Alpha/Omega)
        CALL REPLACE_NAN_ZERO(LocalColor, Beta)
        Rho=Prov

        ! L3: Pi = r(i-1) + Beta*(p(i-1) -omega(i-1)*v(i-1))
        DO IP=1,MNP
          SolDat%AC6(:,:,IP)=SolDat%AC3(:,:,IP)                        &
     &      + Beta(:,:)*SolDat%AC6(:,:,IP)                            &
     &      - Beta(:,:)*Omega(:,:)*SolDat%AC5(:,:,IP)
        END DO

        ! L4 y=K^(-1) Pi
        SolDat%AC1=SolDat%AC6
        IF (PCmethod .gt. 0) THEN
          CALL I5B_APPLY_PRECOND(LocalColor, SolDat, SolDat%AC1)
        ENDIF

        ! L5 vi=Ay
        CALL I5B_APPLY_FCT(MSCeffect, SolDat,  SolDat%AC1, SolDat%AC5)

        ! L6 Alpha=Rho/(hat(r)_0, v_i)
        CALL I5B_SCALAR(MSCeffect, SolDat % AC4, SolDat % AC5, Prov)
        Alpha(:,:)=Rho(:,:)/Prov(:,:)
        CALL REPLACE_NAN_ZERO(LocalColor, Alpha)

        ! L6.1 x(i)=x(i-1) + Alpha y
        DO IP=1,MNP
          SolDat%AC2(:,:,IP)=SolDat%AC2(:,:,IP)                        &
     &      + Alpha(:,:)*SolDat%AC1(:,:,IP)
        END DO

        ! L7 s=r(i-1) - alpha v(i)
        DO IP=1,MNP
          SolDat%AC3(:,:,IP)=SolDat%AC3(:,:,IP)                        &
     &      - Alpha(:,:)*SolDat%AC5(:,:,IP)
        END DO

        ! L8 z=K^(-1) s
        SolDat%AC1=SolDat%AC3
        IF (PCmethod .gt. 0) THEN
          CALL I5B_APPLY_PRECOND(LocalColor, SolDat, SolDat%AC1)
        END IF

        ! L9 t=Az
        CALL I5B_APPLY_FCT(MSCeffect, SolDat,  SolDat%AC1, SolDat%AC7)

        ! L10 omega=(t,s)/(t,t)
        CALL I5B_SCALAR(MSCeffect, SolDat % AC7, SolDat % AC3, Omega)
        CALL I5B_SCALAR(MSCeffect, SolDat % AC7, SolDat % AC7, Prov)
        Omega(:,:)=Omega(:,:)/Prov(:,:)
        CALL REPLACE_NAN_ZERO(LocalColor, Omega)

        ! L11 x(i)=x(i-1) + Omega z
        DO IP=1,MNP
          SolDat%AC2(:,:,IP)=SolDat%AC2(:,:,IP)                        &
     &      + Omega(:,:)*SolDat%AC1(:,:,IP)
        END DO

        ! L12 If x is accurate enough finish
        CALL I5B_APPLY_FCT(MSCeffect, SolDat,  SolDat%AC2, SolDat%AC1)
        CALL I5B_L2_LINF(MSCeffect, SolDat%AC1, SolDat%B_block, Norm_L2, Norm_LINF)
        CritVal=maxval(Norm_L2)
        IF (CritVal .lt. MaxError) THEN
          EXIT
        ENDIF
        IF (nbIter .gt. MaxIter) THEN
          EXIT
        ENDIF

        ! L13 r=s-omega t
        DO IP=1,MNP
          SolDat%AC3(:,:,IP)=SolDat%AC3(:,:,IP)                        &
     &      - Omega(:,:)*SolDat%AC7(:,:,IP)
        END DO
      END DO
      WRITE(STAT%FHNDL, *) 'nbIter=', nbIter
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I4_SPLIT_MSC(LocalColor, NbMSCblock)
      USE DATAPOOL
      implicit none
      type(LocalColorInfo), intent(inout) :: LocalColor
      integer, intent(in) :: NbMSCblock
      integer Hlen, Delta, MSCeffect, istat, iMSCblock, len
      integer ISbegin, IS1, IS2
      Hlen=INT(MSC/NbMSCblock)
      Delta=MSC - Hlen*NbMSCblock
      IF (Delta == 0) THEN
        MSCeffect=Hlen
      ELSE
        MSCeffect=Hlen+1
      ENDIF
      LocalColor % NbMSCblock=NbMSCblock
      LocalColor % MSCeffect=MSCeffect
      allocate(LocalColor % ISbegin(NbMSCblock), LocalColor % ISend(NbMSCblock), LocalColor % ISlen(NbMSCblock), stat=istat)
      IF (istat /=0) CALL WWM_ABORT('allocation error')
      ISbegin=0
      DO iMSCblock=1,NbMSCblock
        IF (iMSCblock <= Delta) THEN
          len=Hlen+1
        ELSE
          len=Hlen
        END IF
        IS1=ISbegin+1
        IS2=ISbegin+len
        LocalColor % ISbegin(iMSCblock)=IS1
        LocalColor % ISend  (iMSCblock)=IS2
        LocalColor % ISlen  (iMSCblock)=len
        ISbegin=ISbegin+len
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_ALLOCATE(SolDat, MSCeffect)
      USE DATAPOOL, only : I5_SolutionData, MNP, MSC, MDC, NNZ
      implicit none
      type(I5_SolutionData), intent(inout) :: SolDat
      integer, intent(in) :: MSCeffect
      integer istat
      allocate(SolDat % AC1(MSCeffect,MDC,MNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 119')
      allocate(SolDat % AC2(MSCeffect,MDC,MNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 120')
      allocate(SolDat % AC3(MSCeffect,MDC,MNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 121')
      allocate(SolDat % AC4(MSCeffect,MDC,MNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 122')
      allocate(SolDat % AC5(MSCeffect,MDC,MNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 123')
      allocate(SolDat % AC6(MSCeffect,MDC,MNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 124')
      allocate(SolDat % AC7(MSCeffect,MDC,MNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 125')
# ifndef BCGS_REORG
      allocate(SolDat % AC8(MSCeffect,MDC,MNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 126')
      allocate(SolDat % AC9(MSCeffect,MDC,MNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 127')
# endif
      allocate(SolDat % ASPAR_block(MSCeffect,MDC,NNZ), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 128')
      allocate(SolDat % B_block(MSCeffect,MDC, MNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 129')
      allocate(SolDat % ASPAR_pc(MSCeffect,MDC,NNZ), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 130')
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WWM_SOLVER_INIT
      implicit none
# ifdef PLAN_I4
      CALL I4_SOLVER_INIT
# else
      CALL I5B_SOLVER_INIT
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_SOLVER_INIT
      USE DATAPOOL
# ifdef DEBUG
      USE elfe_msgp, only : myrank
# endif
      implicit none
      NblockFreqDir = NB_BLOCK
      MainLocalColor%MSCeffect=MSC
      CALL SYMM_INIT_COLORING(MainLocalColor, NblockFreqDir, MSC)
# ifdef DEBUG
      WRITE(myrank+740,*) 'After SYMM_INIT_COLORING'
      CALL FLUSH(myrank+740)
# endif
      CALL I5B_ALLOCATE(SolDat, MSC)
# ifdef DEBUG
      WRITE(myrank+740,*) 'After I5B_ALLOCATE'
      CALL FLUSH(myrank+740)
# endif
      IF (PCmethod .eq. 2) THEN
!        CALL CREATE_ASPAR_EXCHANGE_ARRAY(LocalColor)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I4_SOLVER_INIT
      USE DATAPOOL
      implicit none
      integer NbMSCblock
      NblockFreqDir = NB_BLOCK
      NbMSCblock = 7
      CALL I4_SPLIT_MSC(MainLocalColor, NbMSCblock)
      CALL SYMM_INIT_COLORING(MainLocalColor, NblockFreqDir, MainLocalColor % MSCeffect)
      CALL I5B_ALLOCATE(SolDat, MainLocalColor % MSCeffect)
      IF (PCmethod .eq. 2) THEN
!        CALL CREATE_ASPAR_EXCHANGE_ARRAY(LocalColor)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_FREE(SolDat)
      USE DATAPOOL, only : I5_SolutionData
      implicit none
      type(I5_SolutionData), intent(inout) :: SolDat
      deallocate(SolDat % AC1)
      deallocate(SolDat % AC2)
      deallocate(SolDat % AC3)
      deallocate(SolDat % AC4)
      deallocate(SolDat % AC5)
      deallocate(SolDat % AC6)
      deallocate(SolDat % AC7)
# ifndef BCGS_REORG
      deallocate(SolDat % AC8)
      deallocate(SolDat % AC9)
# endif
      deallocate(SolDat % ASPAR_block)
      deallocate(SolDat % B_block)
      deallocate(SolDat % ASPAR_pc)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5_SUM(AC, eSum)
      USE DATAPOOL, only : MNP, MSC, MDC, rkind, ZERO
      implicit none
      real(rkind), intent(in) :: AC(MNP,MSC,MDC)
      real(rkind), intent(out) :: eSum
      integer IP,IS,ID
      eSum=ZERO
      DO IP=1,MNP
        DO IS=1,MSC
          DO ID=1,MDC
            eSum=eSum + AC(IP,IS,ID)
          END DO
        END DO
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
# ifdef DEBUG
      SUBROUTINE I5_LOCATE_MAX(AC)
      USE DATAPOOL, only : MNP, MSC, MDC, rkind, ZERO
      implicit none
      real(rkind), intent(in) :: AC(MNP,MSC,MDC)
      real(rkind) :: TheMax
      integer IP,IS,ID
      integer IPf,ISf,IDf
      TheMax=ZERO
      IPf=0
      ISf=0
      IDf=0
      DO IP=1,MNP
        DO IS=1,MSC
          DO ID=1,MDC
            IF (AC(IP,IS,ID) .gt. TheMax) THEN
              TheMax=AC(IP,IS,ID)
              IPf=IP
              ISf=IS
              IDf=ID
            END IF
          END DO
        END DO
      END DO
      Print *, 'f IP=', IPf, 'IS=', ISf, 'ID=', IDf
      END SUBROUTINE
# endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WWM_SOLVER_EIMPS(LocalColor, SolDat)
      USE DATAPOOL, only : LocalColorInfo, I5_SolutionData
      implicit none
      type(LocalColorInfo), intent(inout) :: LocalColor
      type(I5_SolutionData), intent(inout) :: SolDat
# if defined PLAN_I4
      CALL I4_EIMPS(LocalColor, SolDat)
# else
      CALL I5B_EIMPS(LocalColor, SolDat)
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I4_CADVXY_VECTOR(LocalColor, CX,CY, iMSCblock)
      USE DATAPOOL
      IMPLICIT NONE
      type(LocalColorInfo), intent(inout) :: LocalColor
      REAL(rkind), INTENT(OUT)  :: CX(LocalColor%MSCeffect,MDC,MNP), CY(LocalColor%MSCeffect,MDC,MNP)
      integer, intent(in) :: iMSCblock
      INTEGER     :: IP, IS, ID, IS1, IS2, ISr
      REAL(rkind)      :: DIFRU, USOC, WVC
!
! Loop over the resident nodes only ... exchange is done in the calling routine
!
      IS1=LocalColor%ISbegin(iMSCblock)
      IS2=LocalColor%ISend  (iMSCblock)
      DO IP = 1, MNP
        DO IS = IS1, IS2
          ISr=IS+1-IS1
          DO ID = 1, MDC
            IF (LSECU .OR. LSTCU) THEN
              CX(ISr,ID,IP) = CG(IP,IS)*COSTH(ID)+CURTXY(IP,1)
              CY(ISr,ID,IP) = CG(IP,IS)*SINTH(ID)+CURTXY(IP,2)
            ELSE
              CX(ISr,ID,IP) = CG(IP,IS)*COSTH(ID)
              CY(ISr,ID,IP) = CG(IP,IS)*SINTH(ID)
            END IF
            IF (LSPHE) THEN
              CX(ISr,ID,IP) = CX(ISr,ID,IP)*INVSPHTRANS(IP,1)
              CY(ISr,ID,IP) = CY(ISr,ID,IP)*INVSPHTRANS(IP,2)
            END IF
            IF (LDIFR) THEN
              CX(ISr,ID,IP) = CX(ISr,ID,IP)*DIFRM(IP)
              CY(ISr,ID,IP) = CY(ISr,ID,IP)*DIFRM(IP)
              IF (LSECU .OR. LSTCU) THEN
                IF (IDIFFR .GT. 1) THEN
                  WVC = SPSIG(IS)/WK(IP,IS)
                  USOC = (COSTH(ID)*CURTXY(IP,1) + SINTH(ID)*CURTXY(IP,2))/WVC
                  DIFRU = ONE + USOC * (ONE - DIFRM(IP))
                ELSE
                  DIFRU = DIFRM(IP)
                END IF
                CX(ISr,ID,IP) = CX(ISr,ID,IP) + DIFRU*CURTXY(IP,1)
                CY(ISr,ID,IP) = CY(ISr,ID,IP) + DIFRU*CURTXY(IP,2)
              END IF
            END IF
          END DO
        END DO
      END DO
      END SUBROUTINE
!**********************************************************************
!*
!**********************************************************************
      SUBROUTINE I4_EIMPS_ASPAR_B_BLOCK(LocalColor, ASPAR, B, U, iMSCblock)
      USE DATAPOOL
# ifdef DEBUG
      USE elfe_msgp, only : myrank
# endif
      IMPLICIT NONE
      type(LocalColorInfo), intent(inout) :: LocalColor
      REAL(rkind), intent(inout) :: ASPAR(LocalColor%MSCeffect, MDC, NNZ)
      REAL(rkind), intent(inout) :: B(LocalColor%MSCeffect, MDC, MNP)
      REAL(rkind), intent(in) :: U(LocalColor%MSCeffect, MDC, MNP)
      integer, intent(in) :: iMSCblock
      INTEGER :: POS_TRICK(3,2)
      REAL(rkind) :: FL11(LocalColor%MSCeffect,MDC), FL12(LocalColor%MSCeffect,MDC), FL21(LocalColor%MSCeffect,MDC), FL22(LocalColor%MSCeffect,MDC), FL31(LocalColor%MSCeffect,MDC), FL32(LocalColor%MSCeffect,MDC)
      REAL(rkind):: CRFS(LocalColor%MSCeffect,MDC,3), K1(LocalColor%MSCeffect,MDC), KM(LocalColor%MSCeffect,MDC,3), K(LocalColor%MSCeffect,MDC,3), TRIA03
      REAL(rkind) :: DELTAL(LocalColor%MSCeffect,MDC,3,MNE)
      INTEGER :: I1, I2, I3
      INTEGER :: IP, ID, IS, IE, POS
      INTEGER :: I, J, IPGL, IPrel, ISr, IS1, IS2
      REAL(rkind) :: KP(LocalColor%MSCeffect,MDC,3,MNE), NM(LocalColor%MSCeffect,MDC,MNE)
      REAL(rkind) :: DTK(LocalColor%MSCeffect,MDC), TMP3(LocalColor%MSCeffect,MDC)
      REAL(rkind) :: LAMBDA(LocalColor%MSCeffect,MDC,2)
# ifdef NO_MEMORY_CX_CY
      REAL(rkind) :: CX(LocalColor%MSCeffect,MDC,3), CY(LocalColor%MSCeffect,MDC,3)
      REAL(rkind) :: USOC, WVC, DIFRU
# else
      REAL(rkind) :: CX(LocalColor%MSCeffect,MDC,MNP), CY(LocalColor%MSCeffect,MDC,MNP)
# endif
      IS1=LocalColor%ISbegin(iMSCblock)
      IS2=LocalColor%ISend  (iMSCblock)

      POS_TRICK(1,1) = 2
      POS_TRICK(1,2) = 3
      POS_TRICK(2,1) = 3
      POS_TRICK(2,2) = 1
      POS_TRICK(3,1) = 1
      POS_TRICK(3,2) = 2

# ifndef NO_MEMORY_CX_CY
      CALL I4_CADVXY_VECTOR(LocalColor, CX,CY, iMSCblock)
# endif
!
!        Calculate countour integral quantities ...
!
      DO IE = 1, MNE
# ifndef NO_MEMORY_CX_CY
        I1 = INE(1,IE)
        I2 = INE(2,IE)
        I3 = INE(3,IE)
        LAMBDA(:,:,1) = ONESIXTH * (CX(:,:,I1) + CX(:,:,I2) + CX(:,:,I3))
        LAMBDA(:,:,2) = ONESIXTH * (CY(:,:,I1) + CY(:,:,I2) + CY(:,:,I3))
        K(:,:,1)  = LAMBDA(:,:,1) * IEN(1,IE) + LAMBDA(:,:,2) * IEN(2,IE)
        K(:,:,2)  = LAMBDA(:,:,1) * IEN(3,IE) + LAMBDA(:,:,2) * IEN(4,IE)
        K(:,:,3)  = LAMBDA(:,:,1) * IEN(5,IE) + LAMBDA(:,:,2) * IEN(6,IE)
        KP(:,:,:,IE) = MAX(ZERO,K)
        KM = MIN(0.0_rkind,K)
        FL11(:,:) = CX(:,:,I2)*IEN(1,IE)+CY(:,:,I2)*IEN(2,IE)
        FL12(:,:) = CX(:,:,I3)*IEN(1,IE)+CY(:,:,I3)*IEN(2,IE)
        FL21(:,:) = CX(:,:,I3)*IEN(3,IE)+CY(:,:,I3)*IEN(4,IE)
        FL22(:,:) = CX(:,:,I1)*IEN(3,IE)+CY(:,:,I1)*IEN(4,IE)
        FL31(:,:) = CX(:,:,I1)*IEN(5,IE)+CY(:,:,I1)*IEN(6,IE)
        FL32(:,:) = CX(:,:,I2)*IEN(5,IE)+CY(:,:,I2)*IEN(6,IE)
# else
        DO I=1,3
          IP=INE(I,IE)
          DO IS = IS1, IS2
            ISr=IS+1-IS1
            DO ID = 1, MDC
              IF (LSECU .OR. LSTCU) THEN
                CX(ISr,ID,I) = CG(IP,IS)*COSTH(ID)+CURTXY(IP,1)
                CY(ISr,ID,I) = CG(IP,IS)*SINTH(ID)+CURTXY(IP,2)
              ELSE
                CX(ISr,ID,I) = CG(IP,IS)*COSTH(ID)
                CY(ISr,ID,I) = CG(IP,IS)*SINTH(ID)
              END IF
              IF (LSPHE) THEN
                CX(ISr,ID,I) = CX(ISr,ID,I)*INVSPHTRANS(IP,1)
                CY(ISr,ID,I) = CY(ISr,ID,I)*INVSPHTRANS(IP,2)
              END IF
              IF (LDIFR) THEN
                CX(ISr,ID,I) = CX(ISr,ID,I)*DIFRM(IP)
                CY(ISr,ID,I) = CY(ISr,ID,I)*DIFRM(IP)
                IF (LSECU .OR. LSTCU) THEN
                  IF (IDIFFR .GT. 1) THEN
                    WVC = SPSIG(IS)/WK(IP,IS)
                    USOC = (COSTH(ID)*CURTXY(IP,1) + SINTH(ID)*CURTXY(IP,2))/WVC
                    DIFRU = ONE + USOC * (ONE - DIFRM(IP))
                  ELSE
                    DIFRU = DIFRM(IP)
                  END IF
                  CX(ISr,ID,I) = CX(ISr,ID,I) + DIFRU*CURTXY(IP,1)
                  CY(ISr,ID,I) = CY(ISr,ID,I) + DIFRU*CURTXY(IP,2)
                END IF
              END IF
            END DO
          END DO
        END DO
        LAMBDA(:,:,1) = ONESIXTH * (CX(:,:,1) + CX(:,:,2) + CX(:,:,3))
        LAMBDA(:,:,2) = ONESIXTH * (CY(:,:,1) + CY(:,:,2) + CY(:,:,3))
        K(:,:,1)  = LAMBDA(:,:,1) * IEN(1,IE) + LAMBDA(:,:,2) * IEN(2,IE)
        K(:,:,2)  = LAMBDA(:,:,1) * IEN(3,IE) + LAMBDA(:,:,2) * IEN(4,IE)
        K(:,:,3)  = LAMBDA(:,:,1) * IEN(5,IE) + LAMBDA(:,:,2) * IEN(6,IE)
        KP(:,:,:,IE) = MAX(ZERO,K)
        KM = MIN(0.0_rkind,K)
        FL11(:,:) = CX(:,:,2)*IEN(1,IE)+CY(:,:,2)*IEN(2,IE)
        FL12(:,:) = CX(:,:,3)*IEN(1,IE)+CY(:,:,3)*IEN(2,IE)
        FL21(:,:) = CX(:,:,3)*IEN(3,IE)+CY(:,:,3)*IEN(4,IE)
        FL22(:,:) = CX(:,:,1)*IEN(3,IE)+CY(:,:,1)*IEN(4,IE)
        FL31(:,:) = CX(:,:,1)*IEN(5,IE)+CY(:,:,1)*IEN(6,IE)
        FL32(:,:) = CX(:,:,2)*IEN(5,IE)+CY(:,:,2)*IEN(6,IE)
# endif
        CRFS(:,:,1) = - ONESIXTH *  (TWO *FL31(:,:) + FL32(:,:) + FL21(:,:) + TWO * FL22(:,:) )
        CRFS(:,:,2) = - ONESIXTH *  (TWO *FL32(:,:) + TWO * FL11(:,:) + FL12(:,:) + FL31(:,:) )
        CRFS(:,:,3) = - ONESIXTH *  (TWO *FL12(:,:) + TWO * FL21(:,:) + FL22(:,:) + FL11(:,:) )
        DELTAL(:,:,:,IE) = CRFS(:,:,:)- KP(:,:,:,IE)
        NM(:,:,IE)=ONE/MIN(-THR,KM(:,:,1) + KM(:,:,2) + KM(:,:,3))
      END DO
# if defined DEBUG
      WRITE(3000+myrank,*)  'sum(LAMBDA)=', sum(LAMBDA)
      WRITE(3000+myrank,*)  'sum(K     )=', sum(K)
      WRITE(3000+myrank,*)  'sum(KP    )=', sum(KP)
      WRITE(3000+myrank,*)  'sum(KM    )=', sum(KM)
      WRITE(3000+myrank,*)  'sum(FL11  )=', sum(FL11)
      WRITE(3000+myrank,*)  'sum(FL12  )=', sum(FL11)
      WRITE(3000+myrank,*)  'sum(FL21  )=', sum(FL21)
      WRITE(3000+myrank,*)  'sum(FL22  )=', sum(FL22)
      WRITE(3000+myrank,*)  'sum(FL31  )=', sum(FL31)
      WRITE(3000+myrank,*)  'sum(FL32  )=', sum(FL32)
      WRITE(3000+myrank,*)  'sum(CRFS  )=', sum(CRFS)
      WRITE(3000+myrank,*)  'sum(DELTAL)=', sum(DELTAL)
      WRITE(3000+myrank,*)  'sum(NM    )=', sum(NM)
      WRITE(3000+myrank,*)  'maxval(NM    )=', maxval(NM)
      DO IS=1,LocalColor%MSCeffect
        DO ID=1,MDC
          DO IE=1,MNE
            IF (ABS(NM(IS,ID,IE)) > 1000000) THEN
              WRITE(4000+myrank,*) 'IS, ID, IE=', IS, ID, IE
              WRITE(4000+myrank,*) '   NM=', NM(IS,ID,IE)
            END IF
          END DO
        END DO
      END DO
# endif
      J     = 0    ! Counter ...
      ASPAR = 0.0_rkind ! Mass matrix ...
      B     = 0.0_rkind ! Right hand side ...
!
! ... assembling the linear equation system ....
!
      DO IP = 1, NP_RES
        IF (IOBWB(IP) .EQ. 1 .AND. DEP(IP) .GT. DMIN) THEN
          DO I = 1, CCON(IP)
            J = J + 1
            IE    =  IE_CELL(J)
            POS   =  POS_CELL(J)
            K1(:,:)    =  KP(:,:,POS,IE) ! Flux Jacobian
            TRIA03 = ONETHIRD * TRIA(IE)
            DO ID=1,MDC
              DTK(:,ID)   =  K1(:,ID) * DT4A * IOBPD(ID,IP)
            END DO
            TMP3(:,:)  =  DTK(:,:) * NM(:,:,IE)
            I1    =  POSI(1,J) ! Position of the recent entry in the ASPAR matrix ... ASPAR is shown in fig. 42, p.122
            I2    =  POSI(2,J)
            I3    =  POSI(3,J)
            ASPAR(:,:,I1) =  TRIA03+DTK(:,:)- TMP3(:,:) * DELTAL(:,:,POS             ,IE) + ASPAR(:,:,I1)  ! Diagonal entry
            ASPAR(:,:,I2) =                 - TMP3(:,:) * DELTAL(:,:,POS_TRICK(POS,1),IE) + ASPAR(:,:,I2)  ! off diagonal entries ...
            ASPAR(:,:,I3) =                 - TMP3(:,:) * DELTAL(:,:,POS_TRICK(POS,2),IE) + ASPAR(:,:,I3)
            DO ID=1,MDC
              IF (IOBPD(ID,IP) .eq. 1) THEN
                B(:,ID,IP)     =  B(:,ID,IP) + TRIA03 * U(:,ID,IP)
              END IF
            END DO
          END DO !I: loop over connected elements ...
        ELSE
          DO I = 1, CCON(IP)
            J = J + 1
            IE    =  IE_CELL(J)
            TRIA03 = ONETHIRD * TRIA(IE)
            I1    =  POSI(1,J) ! Position of the recent entry in the ASPAR matrix ... ASPAR is shown in fig. 42, p.122
            ASPAR(:,:,I1) =  TRIA03 + ASPAR(:,:,I1)  ! Diagonal entry
            B(:,:,IP)     =  ZERO
          END DO
        END IF
      END DO
# if defined DEBUG
      WRITE(3000+myrank,*) 'iMSCblock=', iMSCblock
      WRITE(3000+myrank,*) 'IS12=', IS1, IS2
      DO IS=1,LocalColor%MSCeffect
        WRITE(3000+myrank,*) 'A: IS, sum(ASPAR)=', IS, sum(ASPAR(IS,:,:))
      END DO
# endif
      IF (LBCWA .OR. LBCSP) THEN
        IF (LINHOM) THEN
          IPrel=IP
        ELSE
          IPrel=1
        ENDIF
        DO IP = 1, IWBMNP
          IPGL = IWBNDLC(IP)
          ASPAR(:,:,I_DIAG(IPGL)) = SI(IPGL) ! Set boundary on the diagonal
          DO IS=IS1,IS2
            ISr=IS+1-IS1
            B(ISr,:,IPGL)         = SI(IPGL) * WBAC(IS,:,IPrel)
          END DO
        END DO
      END IF
# if defined DEBUG
      WRITE(3000+myrank,*) 'iMSCblock=', iMSCblock
      WRITE(3000+myrank,*) 'IS12=', IS1, IS2
      DO IS=1,LocalColor%MSCeffect
        WRITE(3000+myrank,*) 'B: IS, sum(ASPAR)=', IS, sum(ASPAR(IS,:,:))
      END DO
# endif
      IF (ICOMP .GE. 2 .AND. SMETHOD .GT. 0) THEN
        DO IP = 1, NP_RES
          IF (IOBWB(IP) .EQ. 1) THEN
            DO IS=IS1,IS2
              ISr=IS+1-IS1
              ASPAR(ISr,:,I_DIAG(IP)) = ASPAR(ISr,:,I_DIAG(IP)) + IMATDAA(IP,IS,:) * DT4A * SI(IP) ! Add source term to the diagonal
              B(ISr,:,IP)             = B(ISr,:,IP) + IMATRAA(IP,IS,:) * DT4A * SI(IP) ! Add source term to the right hand side
            END DO
          ENDIF
        END DO
      ENDIF
# if defined DEBUG
      WRITE(3000+myrank,*) 'sum(ASPAR )=', sum(ASPAR)
      WRITE(3000+myrank,*) 'sum(B     )=', sum(B)
      WRITE(3000+myrank,*) 'iMSCblock=', iMSCblock
      WRITE(3000+myrank,*) 'IS12=', IS1, IS2
      DO IS=1,LocalColor%MSCeffect
        WRITE(3000+myrank,*) 'C: IS, sum(ASPAR)=', IS, sum(ASPAR(IS,:,:))
      END DO
# endif
      END SUBROUTINE
!**********************************************************************
!*
!**********************************************************************
      SUBROUTINE EIMPS_ASPAR_B_BLOCK(ASPAR, B, U)
      USE DATAPOOL
# ifdef DEBUG
      USE elfe_msgp, only : myrank
# endif
      IMPLICIT NONE
      REAL(rkind), intent(inout) :: ASPAR(MSC, MDC, NNZ)
      REAL(rkind), intent(inout) :: B(MSC, MDC, MNP)
      REAL(rkind), intent(in)    :: U(MSC, MDC, MNP)
      INTEGER :: POS_TRICK(3,2)
      REAL(rkind) :: FL11(MSC,MDC), FL12(MSC,MDC), FL21(MSC,MDC), FL22(MSC,MDC), FL31(MSC,MDC), FL32(MSC,MDC)
      REAL(rkind):: CRFS(MSC,MDC,3), K1(MSC,MDC), KM(MSC,MDC,3), K(MSC,MDC,3), TRIA03
# ifndef NO_MEMORY_CX_CY
      REAL(rkind) :: CX(MSC,MDC,MNP), CY(MSC,MDC,MNP)
# else
      REAL(rkind) :: CX(MSC,MDC,3), CY(MSC,MDC,3)
      REAL(rkind)      :: DIFRU, USOC, WVC
# endif
      REAL(rkind) :: DELTAL(MSC,MDC,3,MNE)
      INTEGER :: I1, I2, I3
      INTEGER :: IP, ID, IS, IE, POS
      INTEGER :: I, J, IPGL, IPrel
      REAL(rkind) :: KP(MSC,MDC,3,MNE), NM(MSC,MDC,MNE)
      REAL(rkind) :: DTK(MSC,MDC), TMP3(MSC,MDC)
      REAL(rkind) :: LAMBDA(MSC,MDC,2)
# ifdef DEBUG
      WRITE(740+myrank,*) 'Begin of EIMPS_ASPAR_B_BLOCK'
# endif
      POS_TRICK(1,1) = 2
      POS_TRICK(1,2) = 3
      POS_TRICK(2,1) = 3
      POS_TRICK(2,2) = 1
      POS_TRICK(3,1) = 1
      POS_TRICK(3,2) = 2

# ifndef NO_MEMORY_CX_CY
      CALL CADVXY_VECTOR(CX, CY)
# endif
!
!        Calculate countour integral quantities ...
!
# ifdef DEBUG
      WRITE(740+myrank,*) ' Before MNE loop'
# endif
      DO IE = 1, MNE
# ifndef NO_MEMORY_CX_CY
        I1 = INE(1,IE)
        I2 = INE(2,IE)
        I3 = INE(3,IE)
        LAMBDA(:,:,1) = ONESIXTH * (CX(:,:,I1) + CX(:,:,I2) + CX(:,:,I3))
        LAMBDA(:,:,2) = ONESIXTH * (CY(:,:,I1) + CY(:,:,I2) + CY(:,:,I3))
        K(:,:,1)  = LAMBDA(:,:,1) * IEN(1,IE) + LAMBDA(:,:,2) * IEN(2,IE)
        K(:,:,2)  = LAMBDA(:,:,1) * IEN(3,IE) + LAMBDA(:,:,2) * IEN(4,IE)
        K(:,:,3)  = LAMBDA(:,:,1) * IEN(5,IE) + LAMBDA(:,:,2) * IEN(6,IE)
        KP(:,:,:,IE) = MAX(ZERO,K)
        KM = MIN(0.0_rkind,K)
        FL11(:,:) = CX(:,:,I2)*IEN(1,IE)+CY(:,:,I2)*IEN(2,IE)
        FL12(:,:) = CX(:,:,I3)*IEN(1,IE)+CY(:,:,I3)*IEN(2,IE)
        FL21(:,:) = CX(:,:,I3)*IEN(3,IE)+CY(:,:,I3)*IEN(4,IE)
        FL22(:,:) = CX(:,:,I1)*IEN(3,IE)+CY(:,:,I1)*IEN(4,IE)
        FL31(:,:) = CX(:,:,I1)*IEN(5,IE)+CY(:,:,I1)*IEN(6,IE)
        FL32(:,:) = CX(:,:,I2)*IEN(5,IE)+CY(:,:,I2)*IEN(6,IE)
# else
        DO I=1,3
          IP = INE(I,IE)
          DO IS=1,MSC
            DO ID=1,MDC
              IF (LSECU .OR. LSTCU) THEN
                CX(IS,ID,I) = CG(IP,IS)*COSTH(ID)+CURTXY(IP,1)
                CY(IS,ID,I) = CG(IP,IS)*SINTH(ID)+CURTXY(IP,2)
              ELSE
                CX(IS,ID,I) = CG(IP,IS)*COSTH(ID)
                CY(IS,ID,I) = CG(IP,IS)*SINTH(ID)
              END IF
              IF (LSPHE) THEN
                CX(IS,ID,I) = CX(IS,ID,I)*INVSPHTRANS(IP,1)
                CY(IS,ID,I) = CY(IS,ID,I)*INVSPHTRANS(IP,2)
              END IF
              IF (LDIFR) THEN
                CX(IS,ID,I) = CX(IS,ID,I)*DIFRM(IP)
                CY(IS,ID,I) = CY(IS,ID,I)*DIFRM(IP)
                IF (LSECU .OR. LSTCU) THEN
                  IF (IDIFFR .GT. 1) THEN
                    WVC = SPSIG(IS)/WK(IP,IS)
                    USOC = (COSTH(ID)*CURTXY(IP,1) + SINTH(ID)*CURTXY(IP,2))/WVC
                    DIFRU = ONE + USOC * (ONE - DIFRM(IP))
                  ELSE
                    DIFRU = DIFRM(IP)
                  END IF
                  CX(IS,ID,I) = CX(IS,ID,I) + DIFRU*CURTXY(IP,1)
                  CY(IS,ID,I) = CY(IS,ID,I) + DIFRU*CURTXY(IP,2)
                END IF
              END IF
            END DO
          END DO
        END DO
        LAMBDA(:,:,1) = ONESIXTH * (CX(:,:,1) + CX(:,:,2) + CX(:,:,3))
        LAMBDA(:,:,2) = ONESIXTH * (CY(:,:,1) + CY(:,:,2) + CY(:,:,3))
        K(:,:,1)  = LAMBDA(:,:,1) * IEN(1,IE) + LAMBDA(:,:,2) * IEN(2,IE)
        K(:,:,2)  = LAMBDA(:,:,1) * IEN(3,IE) + LAMBDA(:,:,2) * IEN(4,IE)
        K(:,:,3)  = LAMBDA(:,:,1) * IEN(5,IE) + LAMBDA(:,:,2) * IEN(6,IE)
        KP(:,:,:,IE) = MAX(ZERO,K)
        KM = MIN(0.0_rkind,K)
        FL11(:,:) = CX(:,:,2)*IEN(1,IE)+CY(:,:,2)*IEN(2,IE)
        FL12(:,:) = CX(:,:,3)*IEN(1,IE)+CY(:,:,3)*IEN(2,IE)
        FL21(:,:) = CX(:,:,3)*IEN(3,IE)+CY(:,:,3)*IEN(4,IE)
        FL22(:,:) = CX(:,:,1)*IEN(3,IE)+CY(:,:,1)*IEN(4,IE)
        FL31(:,:) = CX(:,:,1)*IEN(5,IE)+CY(:,:,1)*IEN(6,IE)
        FL32(:,:) = CX(:,:,2)*IEN(5,IE)+CY(:,:,2)*IEN(6,IE)
# endif
        CRFS(:,:,1) = - ONESIXTH *  (TWO *FL31(:,:) + FL32(:,:) + FL21(:,:) + TWO * FL22(:,:) )
        CRFS(:,:,2) = - ONESIXTH *  (TWO *FL32(:,:) + TWO * FL11(:,:) + FL12(:,:) + FL31(:,:) )
        CRFS(:,:,3) = - ONESIXTH *  (TWO *FL12(:,:) + TWO * FL21(:,:) + FL22(:,:) + FL11(:,:) )
        DELTAL(:,:,:,IE) = CRFS(:,:,:)- KP(:,:,:,IE)
        NM(:,:,IE)=ONE/MIN(-THR,KM(:,:,1) + KM(:,:,2) + KM(:,:,3))
      END DO
# if defined DEBUG
      WRITE(740+myrank,*) ' After MNE loop'
      WRITE(3000+myrank,*)  'sum(LAMBDA)=', sum(LAMBDA)
      WRITE(3000+myrank,*)  'sum(K     )=', sum(K)
      WRITE(3000+myrank,*)  'sum(KP    )=', sum(KP)
      WRITE(3000+myrank,*)  'sum(KM    )=', sum(KM)
      WRITE(3000+myrank,*)  'sum(FL11  )=', sum(FL11)
      WRITE(3000+myrank,*)  'sum(FL12  )=', sum(FL11)
      WRITE(3000+myrank,*)  'sum(FL21  )=', sum(FL21)
      WRITE(3000+myrank,*)  'sum(FL22  )=', sum(FL22)
      WRITE(3000+myrank,*)  'sum(FL31  )=', sum(FL31)
      WRITE(3000+myrank,*)  'sum(FL32  )=', sum(FL32)
      WRITE(3000+myrank,*)  'sum(CRFS  )=', sum(CRFS)
      WRITE(3000+myrank,*)  'sum(DELTAL)=', sum(DELTAL)
      WRITE(3000+myrank,*)  'sum(NM    )=', sum(NM)
      WRITE(3000+myrank,*)  'maxval(NM    )=', maxval(NM)
      DO IS=1,MSC
        DO ID=1,MDC
          DO IE=1,MNE
            IF (ABS(NM(IS,ID,IE)) > 1000000) THEN
              WRITE(4000+myrank,*) 'IS, ID, IE=', IS, ID, IE
              WRITE(4000+myrank,*) '   NM=', NM(IS,ID,IE)
            END IF
          END DO
        END DO
      END DO
# endif
      J     = 0    ! Counter ...
      ASPAR = 0.0_rkind ! Mass matrix ...
      B     = 0.0_rkind ! Right hand side ...
!
! ... assembling the linear equation system ....
!
      DO IP = 1, NP_RES
        IF (IOBWB(IP) .EQ. 1 .AND. DEP(IP) .GT. DMIN) THEN
          DO I = 1, CCON(IP)
            J = J + 1
            IE    =  IE_CELL(J)
            POS   =  POS_CELL(J)
            K1(:,:)    =  KP(:,:,POS,IE) ! Flux Jacobian
            TRIA03 = ONETHIRD * TRIA(IE)
            DO ID=1,MDC
              DTK(:,ID)   =  K1(:,ID) * DT4A * IOBPD(ID,IP)
            END DO
            TMP3(:,:)  =  DTK(:,:) * NM(:,:,IE)
            I1    =  POSI(1,J) ! Position of the recent entry in the ASPAR matrix ... ASPAR is shown in fig. 42, p.122
            I2    =  POSI(2,J)
            I3    =  POSI(3,J)
            ASPAR(:,:,I1) =  TRIA03+DTK(:,:)- TMP3(:,:) * DELTAL(:,:,POS             ,IE) + ASPAR(:,:,I1)  ! Diagonal entry
            ASPAR(:,:,I2) =                 - TMP3(:,:) * DELTAL(:,:,POS_TRICK(POS,1),IE) + ASPAR(:,:,I2)  ! off diagonal entries ...
            ASPAR(:,:,I3) =                 - TMP3(:,:) * DELTAL(:,:,POS_TRICK(POS,2),IE) + ASPAR(:,:,I3)
            DO ID=1,MDC
              B(:,ID,IP)     =  B(:,ID,IP) + IOBPD(ID,IP)*TRIA03 * U(:,ID,IP)
            END DO
          END DO !I: loop over connected elements ...
        ELSE
          DO I = 1, CCON(IP)
            J = J + 1
            IE    =  IE_CELL(J)
            TRIA03 = ONETHIRD * TRIA(IE)
            I1    =  POSI(1,J) ! Position of the recent entry in the ASPAR matrix ... ASPAR is shown in fig. 42, p.122
            ASPAR(:,:,I1) =  TRIA03 + ASPAR(:,:,I1)  ! Diagonal entry
            B(:,:,IP)     =  ZERO
          END DO
        END IF
      END DO
      IF (LBCWA .OR. LBCSP) THEN
        IF (LINHOM) THEN
          IPrel=IP
        ELSE
          IPrel=1
        ENDIF
        DO IP = 1, IWBMNP
          IPGL = IWBNDLC(IP)
          ASPAR(:,:,I_DIAG(IPGL)) = SI(IPGL) ! Set boundary on the diagonal
          B(:,:,IPGL)             = SI(IPGL) * WBAC(:,:,IPrel)
        END DO
      END IF
      IF (ICOMP .GE. 2 .AND. SMETHOD .GT. 0) THEN
        DO IP = 1, NP_RES
          IF (IOBWB(IP) .EQ. 1) THEN
            ASPAR(:,:,I_DIAG(IP)) = ASPAR(:,:,I_DIAG(IP)) + IMATDAA(IP,:,:) * DT4A * SI(IP) ! Add source term to the diagonal
            B(:,:,IP)             = B(:,:,IP) + IMATRAA(IP,:,:) * DT4A * SI(IP) ! Add source term to the right hand side
          ENDIF
        END DO
      ENDIF
# if defined DEBUG
      WRITE(3000+myrank,*)  'sum(ASPAR )=', sum(ASPAR)
      WRITE(3000+myrank,*)  'sum(B     )=', sum(B)
      DO IS=1,MSC
        WRITE(3000+myrank,*) 'IS, sum(ASPAR)=', IS, sum(ASPAR(IS,:,:))
      END DO
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_EIMPS(LocalColor, SolDat)
      USE DATAPOOL, only : LocalColorInfo, I5_SolutionData
      USE DATAPOOL, only : rkind, MSC, MDC, AC2, MNP, NNZ
      USE DATAPOOL, only : PCmethod, IOBPD, ZERO
      USE elfe_msgp, only : myrank, exchange_p4d_wwm
      implicit none
      type(LocalColorInfo), intent(inout) :: LocalColor
      type(I5_SolutionData), intent(inout) :: SolDat
# ifndef ASPAR_B_COMPUTE_BLOCK
      real(rkind) :: U(MNP), ASPAR(NNZ), B(MNP)
# endif
      integer IS, ID, IP
# ifdef DEBUG
      WRITE(740+myrank,*) 'Begin I5B_EIMPS'
# endif
      DO IP=1,MNP
        SolDat % AC2(:,:,IP)=AC2(IP,:,:)
      END DO
# if defined ASPAR_B_COMPUTE_BLOCK
#  ifdef DEBUG
      WRITE(740+myrank,*) 'MSC=', MSC
      WRITE(740+myrank,*) 'MDC=', MDC
      WRITE(740+myrank,*) 'MNP=', MNP
      WRITE(740+myrank,*) 'NNZ=', NNZ
      WRITE(740+myrank,*) 'Before EIMPS_ASPAR_B_BLOCK'
#  endif
      CALL EIMPS_ASPAR_B_BLOCK(SolDat%ASPAR_block, SolDat%B_block, SolDat%AC2)
#  ifdef DEBUG
      WRITE(740+myrank,*) 'After EIMPS_ASPAR_B_BLOCK'
#  endif
!      DO IS=1,MSC
!        DO ID=1,MDC
!          WRITE(6000+myrank,*) 'IS=', IS, 'ID=', ID
!          WRITE(6000+myrank,*) 'sum, ASPAR=', sum(SolDat%ASPAR_block(IS,ID,:)), 'B=', sum(SolDat%B_block(ID,ID,:))
!        END DO
!      END DO
# else
      DO IS=1,MSC
        DO ID=1,MDC
          U=AC2(:,IS,ID)
          CALL EIMPS_ASPAR_B(IS, ID, ASPAR, B, U)
          SolDat % ASPAR_block(IS,ID,:)=ASPAR
          SolDat % B_block(IS,ID,:)=B
!          WRITE(6000+myrank,*) 'IS=', IS, 'ID=', ID
!          WRITE(6000+myrank,*) 'sum, ASPAR=', sum(SolDat%ASPAR_block(IS,ID,:)), 'B=', sum(SolDat%B_block(ID,ID,:))
        END DO
      END DO
# endif
# ifdef DEBUG
      WRITE(740+myrank,*) 'After ASPAR init'
# endif
# ifdef SELFE_EXCHANGE
      CALL I5B_EXCHANGE_P4D_WWM(LocalColor, SolDat % B_block)
# else
      CALL EXCHANGE_P4D_WWM(SolDat % B_block)
# endif
# ifdef DEBUG
      WRITE(740+myrank,*) 'After EXCHANGE_P4D_WWM'
# endif
      CALL I5B_EXCHANGE_ASPAR(LocalColor, SolDat%ASPAR_block)
# ifdef DEBUG
      WRITE(740+myrank,*) 'After I5B_EXCHANGE_ASPAR'
# endif
!      IF (myrank .eq. 0) THEN
!        Print *, 'Before CREATE_PRECOND'
!      END IF
      CALL I5B_CREATE_PRECOND(LocalColor, SolDat, PCmethod)
# ifdef DEBUG
      WRITE(740+myrank,*) 'After I5B_CREATE_PRECOND'
# endif
# ifdef BCGS_REORG
      CALL I5B_BCGS_REORG_SOLVER(LocalColor, SolDat)
# else
      CALL I5B_BCGS_SOLVER(LocalColor, SolDat)
# endif
      DO IP=1,MNP
        DO IS=1,MSC
          DO ID=1,MDC
            AC2(IP,IS,ID)=MAX(ZERO, SolDat%AC2(IS,ID,IP))*MyREAL(IOBPD(ID,IP))
          END DO
        END DO
      END DO
# ifdef DEBUG
      IF (myrank == 1) THEN
        Write(myrank+591,*) 'Clearing ENDING'
      END IF
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I4_EIMPS(LocalColor, SolDat)
      USE DATAPOOL, only : LocalColorInfo, I5_SolutionData
      USE DATAPOOL, only : rkind, MSC, MDC, AC2, MNP, NNZ
      USE DATAPOOL, only : PCmethod, IOBPD, ZERO
      USE elfe_msgp, only : myrank, exchange_p4d_wwm
      implicit none
      type(LocalColorInfo), intent(inout) :: LocalColor
      type(I5_SolutionData), intent(inout) :: SolDat
      integer IS, ID, IP
      integer iMSCblock, IS1, IS2
      DO iMSCblock=1,LocalColor % NbMSCblock
# ifdef DEBUG
        WRITE(240+myrank,*) 'iMSCblock=', iMSCblock
# endif
        IS1=LocalColor%ISbegin(iMSCblock)
        IS2=LocalColor%ISbegin(iMSCblock)
        DO IP=1,MNP
          SolDat % AC2(1:LocalColor%ISlen(iMSCblock),:,IP)=AC2(IP,IS1:IS2,:)
        END DO
        CALL I4_EIMPS_ASPAR_B_BLOCK(LocalColor, SolDat%ASPAR_block, SolDat%B_block, SolDat%AC2, iMSCblock)
        CALL I5B_EXCHANGE_P4D_WWM(LocalColor, SolDat%B_block)
        CALL I5B_EXCHANGE_ASPAR(LocalColor, SolDat%ASPAR_block)
        CALL I5B_CREATE_PRECOND(LocalColor, SolDat, PCmethod)
# ifdef BCGS_REORG
        CALL I5B_BCGS_REORG_SOLVER(LocalColor, SolDat)
# else
        CALL I5B_BCGS_SOLVER(LocalColor, SolDat)
# endif
        DO IP=1,MNP
          DO IS=IS1,IS2
            DO ID=1,MDC
              AC2(IP,IS,ID)=MAX(ZERO, SolDat%AC2(IS+1-IS1,ID,IP))*MyREAL(IOBPD(ID,IP))
            END DO
          END DO
        END DO
      END DO
      END SUBROUTINE
#endif
END MODULE WWM_PARALL_SOLVER
