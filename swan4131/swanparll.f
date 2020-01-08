!
!     SWAN - routines for distributed-memory approach based on MPI
!
!  Contents of this file
!
!     SWINITMPI
!     SWEXITMPI
!     SWSYNC
!     SWSENDNB
!     SWRECVNB
!     SWBROADC
!     SWGATHER
!     SWREDUCE
!     SWREDUCI
!     SWREDUCR
!     SWSTRIP
!JAC!     SWORB
!     SWPARTIT
!     SWBLADM
!     SWDECOMP
!     SWEXCHG
!     SWRECVAC
!     SWSENDAC
!     SWCOLLECT
!     SWCOLOUT
!     SWCOLTAB
!     SWCOLSPC
!     SWCOLBLK
!JAC!     SWBLKCOL
!
!****************************************************************
!
      SUBROUTINE SWINITMPI
!
!****************************************************************
!
!MPI      USE MPI
      USE OCPCOMM4                                                        40.41
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!COH!     41.61: Homayoon Komijani
!
!  1. Updates
!
!     40.30, Jan. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!COH!     41.61, Aug. 14: coupling with COHERENS
!
!  2. Purpose
!
!     Join parallel application
!
!  3. Method
!
!     Start MPI and initialize some variables
!
!  4. Argument variables
!
!     ---
!
!  5. Parameter variables
!
!     ---
!
!  6. Local variables
!
!     CHARS :     array to pass character info to MSGERR
!     IERR  :     error value of MPI call
!     IF1   :     first non-character in string1
!     IF2   :     first non-character in string2
!     IL1   :     last non-character in string1
!     IL2   :     last non-character in string2
!     MSGSTR:     string to pass message to call MSGERR
!
      INTEGER      IERR, IF1, IF2, IL1, IL2
      CHARACTER*20 INTSTR, CHARS(2)
      CHARACTER*80 MSGSTR
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!MPI!     MPI_COMM_RANK    Get rank of processes in MPI communication context
!MPI!     MPI_COMM_SIZE    Get number of processes in MPI communication context
!MPI!     MPI_INIT         Enroll in MPI
!     MSGERR           Writes error message
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     Main program SWAN
!
! 10. Error messages
!
!     ---
!
! 12. Structure
!
!     Start MPI and initialize some common variables in module M_PARALL
!
! 13. Source text
!
      LEVERR = 0
      MAXERR = 1
      ITRACE = 0

!MPI!     --- enroll in MPI
!MPI
!MPI!NCOH      CALL MPI_INIT ( IERR )
!MPI!COH      IERR = MPI_SUCCESS                                                  41.61
!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
!MPI         CHARS(1) = INTSTR(IERR)
!MPI         CALL TXPBLA(CHARS(1),IF1,IL1)
!MPI         MSGSTR = 'MPI produces some internal error - '//
!MPI     &            'return code is '//CHARS(1)(IF1:IL1)
!MPI         CALL MSGERR ( 4, MSGSTR )
!MPI         RETURN
!MPI      END IF

!     --- initialize common variables

      INODE = 0
      NPROC = 1

!MPI!     --- get node number INODE
!MPI
!MPI!NCOH      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, INODE, IERR )
!MPI!COH      CALL MPI_COMM_RANK ( COMM, INODE, IERR )
      INODE = INODE + 1
!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
!MPI         CHARS(1) = INTSTR(IERR)
!MPI         CALL TXPBLA(CHARS(1),IF1,IL1)
!MPI         CHARS(2) = INTSTR(INODE)
!MPI         CALL TXPBLA(CHARS(2),IF2,IL2)
!MPI         MSGSTR = 'MPI produces some internal error - '//
!MPI     &            'return code is '//CHARS(1)(IF1:IL1)//
!MPI     &            ' and node number is '//CHARS(2)(IF2:IL2)
!MPI         CALL MSGERR ( 4, MSGSTR )
!MPI         RETURN
!MPI      END IF
!MPI
!MPI!     --- determine total number of processes
!MPI
!MPI!NCOH      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, NPROC, IERR )
!MPI!COH      CALL MPI_COMM_SIZE ( COMM, NPROC, IERR )
!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
!MPI         CHARS(1) = INTSTR(IERR)
!MPI         CALL TXPBLA(CHARS(1),IF1,IL1)
!MPI         CHARS(2) = INTSTR(INODE)
!MPI         CALL TXPBLA(CHARS(2),IF2,IL2)
!MPI         MSGSTR = 'MPI produces some internal error - '//
!MPI     &            'return code is '//CHARS(1)(IF1:IL1)//
!MPI     &            ' and node number is '//CHARS(2)(IF2:IL2)
!MPI         CALL MSGERR ( 4, MSGSTR )
!MPI         RETURN
!MPI      END IF

!     --- determine whether this is a parallel run or not

      IF ( NPROC.GT.1 ) THEN
         PARLL = .TRUE.
      ELSE
         PARLL = .FALSE.
      END IF

!     --- am I master?

      IAMMASTER = INODE.EQ.MASTER

!MPI!     --- define MPI constants for communication within SWAN
!MPI
!MPI      SWINT  = MPI_INTEGER
!MPI      SWREAL = MPI_REAL
!MPI      SWMAX  = MPI_MAX
!MPI      SWMIN  = MPI_MIN
!MPI      SWSUM  = MPI_SUM

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWEXITMPI
!
!****************************************************************
!
!MPI      USE MPI
      USE OCPCOMM4                                                        40.41
!COH      USE M_PARALL, ONLY: COMM                                            41.61
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!COH!     41.61: Homayoon Komijani
!
!  1. Updates
!
!     40.30, Jan. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!COH!     41.61, Aug. 14: coupling with COHERENS
!
!  2. Purpose
!
!     Exit parallel application
!
!  3. Method
!
!     Wrapper for MPI_FINALIZE
!
!  4. Argument variables
!
!     ---
!
!  6. Local variables
!
!     IERR    :   error value of MPI call
!     PARALMPI:   if true, parallel process is carried out with MPI
!
      INTEGER IERR
      LOGICAL PARALMPI
!
!  8. Subroutines used
!
!MPI!     MPI_ABORT        Abort MPI if severe error occurs
!MPI!     MPI_BARRIER      Blocks until all nodes have called this routine
!MPI!     MPI_INITIALIZED  Indicates whether MPI_Init has been called
!MPI!     MPI_FINALIZE     Cleans up the MPI state and exits
!
!  9. Subroutines calling
!
!     Main program SWAN
!
! 10. Error messages
!
!     ---
!
! 12. Structure
!
!MPI!     if MPI has been initialized
!MPI!        synchronize nodes
!MPI!        if severe error
!MPI!           abort MPI
!MPI!        else
!MPI!           close MPI
!MPI!
! 13. Source text
!
!MPI      CALL MPI_INITIALIZED ( PARALMPI, IERR )
!MPI      IF ( PARALMPI ) THEN
!MPI
!MPI!NCOH         CALL MPI_BARRIER ( MPI_COMM_WORLD, IERR )
!MPI!COH         CALL MPI_BARRIER ( COMM, IERR )
!MPI
!MPI         IF ( LEVERR.GE.4 ) THEN
!MPI
!MPI!        --- in case of a severe error abort all MPI processes
!MPI
!MPI!NCOH            CALL MPI_ABORT ( MPI_COMM_WORLD, LEVERR, IERR )
!MPI!COH            CALL MPI_ABORT ( COMM, LEVERR, IERR )
!MPI
!MPI         ELSE
!MPI
!MPI!        --- otherwise stop MPI operations on this computer
!MPI
!MPI!NCOH            CALL MPI_FINALIZE ( IERR )
!MPI
!MPI         END IF
!MPI
!MPI      END IF
!MPI
      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWSYNC
!
!****************************************************************
!
!MPI      USE MPI
      USE OCPCOMM4                                                        40.41
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Jan. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Synchronize nodes
!
!  3. Method
!
!     Wrapper for MPI_BARRIER
!
!  4. Argument variables
!
!     ---
!
!  5. Parameter variables
!
!     ---
!
!  6. Local variables
!
!     CHARS :     array to pass character info to MSGERR
!     IENT  :     number of entries
!     IERR  :     error value of MPI call
!     IF1   :     first non-character in string1
!     IF2   :     first non-character in string2
!     IL1   :     last non-character in string1
!     IL2   :     last non-character in string2
!     MSGSTR:     string to pass message to call MSGERR
!
      INTEGER      IENT, IERR, IF1, IF2, IL1, IL2
      CHARACTER*20 INTSTR, CHARS(2)
      CHARACTER*80 MSGSTR
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!MPI!     MPI_BARRIER      Blocks until all nodes have called this routine
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SWMAIN
!
! 10. Error messages
!
!     ---
!
! 12. Structure
!
!     Blocks until all nodes have called MPI_BARRIER routine.
!     In this way, all nodes are synchronized
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWSYNC')

!MPI!     --- blocks until all nodes have called this routine
!MPI
!MPI!NCOH      CALL MPI_BARRIER ( MPI_COMM_WORLD, IERR )
!MPI!COH      CALL MPI_BARRIER ( COMM, IERR )
!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
!MPI         CHARS(1) = INTSTR(IERR)
!MPI         CALL TXPBLA(CHARS(1),IF1,IL1)
!MPI         CHARS(2) = INTSTR(INODE)
!MPI         CALL TXPBLA(CHARS(2),IF2,IL2)
!MPI         MSGSTR = 'MPI produces some internal error - '//
!MPI     &            'return code is '//CHARS(1)(IF1:IL1)//
!MPI     &            ' and node number is '//CHARS(2)(IF2:IL2)
!MPI         CALL MSGERR ( 4, MSGSTR )
!MPI         RETURN
!MPI      END IF

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWSENDNB ( IPTR, ILEN, ITYPE, IDEST, ITAG )
!
!****************************************************************
!
!MPI      USE MPI
      USE OCPCOMM4                                                        40.41
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Data is sent to a neighbour
!
!  3. Method
!
!     Wrapper for MPI_SEND
!
!  4. Argument variables
!
!     IDEST       rank of the destination process
!     ILEN        length of array to be sent
!     IPTR        pointer to first element of array to be sent
!     ITAG        message type
!     ITYPE       type of data
!
      INTEGER IPTR, ILEN, ITYPE, IDEST, ITAG
!
!  5. Parameter variables
!
!     ---
!
!  6. Local variables
!
!     CHARS :     array to pass character info to MSGERR
!     IENT  :     number of entries
!     IERR  :     error value of MPI call
!     IF1   :     first non-character in string1
!     IF2   :     first non-character in string2
!     IL1   :     last non-character in string1
!     IL2   :     last non-character in string2
!     MSGSTR:     string to pass message to call MSGERR
!
      INTEGER      IENT, IERR, IF1, IF2, IL1, IL2
      CHARACTER*20 INTSTR, CHARS(2)
      CHARACTER*80 MSGSTR
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!MPI!     MPI_SEND         Immediately sends the data in the active
!MPI!                      MPI message buffer
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SWEXCHG
!
! 10. Error messages
!
!     ---
!
! 12. Structure
!
!MPI!     Data is sent to a neighbour with command MPI_SEND
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWSENDNB')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

!MPI!NCOH      CALL MPI_SEND ( IPTR, ILEN, ITYPE, IDEST-1,
!MPI!NCOH     &                ITAG, MPI_COMM_WORLD, IERR )
!MPI!COH      CALL MPI_SEND ( IPTR, ILEN, ITYPE, IDEST-1,
!MPI!COH     &                ITAG, COMM, IERR )
!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
!MPI         CHARS(1) = INTSTR(IERR)
!MPI         CALL TXPBLA(CHARS(1),IF1,IL1)
!MPI         CHARS(2) = INTSTR(INODE)
!MPI         CALL TXPBLA(CHARS(2),IF2,IL2)
!MPI         MSGSTR = 'MPI produces some internal error - '//
!MPI     &            'return code is '//CHARS(1)(IF1:IL1)//
!MPI     &            ' and node number is '//CHARS(2)(IF2:IL2)
!MPI         CALL MSGERR ( 4, MSGSTR )
!MPI         RETURN
!MPI      END IF

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWRECVNB ( IPTR, ILEN, ITYPE, ISOURCE, ITAG )
!
!****************************************************************
!
!MPI      USE MPI
      USE OCPCOMM4                                                        40.41
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Data is received from a neighbour
!
!  3. Method
!
!     Wrapper for MPI_RECV
!
!  4. Argument variables
!
!     ILEN        length of array to be received
!     IPTR        pointer to first element of array to be received
!     ISOURCE     rank of the source process
!     ITAG        message type
!     ITYPE       type of data
!
      INTEGER IPTR, ILEN, ITYPE, ISOURCE, ITAG
!
!  5. Parameter variables
!
!     ---
!
!  6. Local variables
!
!     CHARS :     array to pass character info to MSGERR
!     IENT  :     number of entries
!     IERR  :     error value of MPI call
!     IF1   :     first non-character in string1
!     IF2   :     first non-character in string2
!     IL1   :     last non-character in string1
!     IL2   :     last non-character in string2
!MPI!     ISTAT :     MPI status array
!     MSGSTR:     string to pass message to call MSGERR
!
      INTEGER      IENT, IERR, IF1, IF2, IL1, IL2
!MPI      INTEGER      ISTAT(MPI_STATUS_SIZE)
      CHARACTER*20 INTSTR, CHARS(2)
      CHARACTER*80 MSGSTR
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!MPI!     MPI_RECV         Immediately receives the data in the active
!MPI!                      MPI message buffer
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SWEXCHG
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!MPI!     Data is received from a neighbour with command MPI_RECV
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWRECVNB')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

!MPI!NCOH      CALL MPI_RECV ( IPTR, ILEN, ITYPE, ISOURCE-1, ITAG,
!MPI!NCOH     &                MPI_COMM_WORLD, ISTAT, IERR )
!MPI!COH      CALL MPI_RECV ( IPTR, ILEN, ITYPE, ISOURCE-1, ITAG,
!MPI!COH     &                COMM, ISTAT, IERR )
!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
!MPI         CHARS(1) = INTSTR(IERR)
!MPI         CALL TXPBLA(CHARS(1),IF1,IL1)
!MPI         CHARS(2) = INTSTR(INODE)
!MPI         CALL TXPBLA(CHARS(2),IF2,IL2)
!MPI         MSGSTR = 'MPI produces some internal error - '//
!MPI     &            'return code is '//CHARS(1)(IF1:IL1)//
!MPI     &            ' and node number is '//CHARS(2)(IF2:IL2)
!MPI         CALL MSGERR ( 4, MSGSTR )
!MPI         RETURN
!MPI      END IF

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWBROADC ( IPTR, ILEN, ITYPE )
!
!****************************************************************
!
!MPI      USE MPI
      USE OCPCOMM4                                                        40.41
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Broadcasts data from the master to all other processes
!
!  3. Method
!
!     Wrapper for MPI_BCAST
!
!  4. Argument variables
!
!     ILEN        length of array to be sent
!     IPTR        pointer to first element of array to be sent
!     ITYPE       type of data
!
      INTEGER IPTR, ILEN, ITYPE
!
!  5. Parameter variables
!
!     ---
!
!  6. Local variables
!
!     CHARS :     character for passing info to MSGERR
!     IENT  :     number of entries
!     IERR  :     error value of MPI call
!     IF    :     first non-character in string
!     IL    :     last non-character in string
!     MSGSTR:     string to pass message to call MSGERR
!
      INTEGER      IENT, IERR, IF, IL
      CHARACTER*20 INTSTR, CHARS
      CHARACTER*80 MSGSTR
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!MPI!     MPI_BCAST        Broadcasts a message from the master
!MPI!                      to all other processes of the group
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!TIMG!     SWTSTA           Start timing for a section of code
!TIMG!     SWTSTO           Stop timing for a section of code
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SNEXTI
!     FLFILE
!
! 10. Error messages
!
!     ---
!
! 12. Structure
!
!MPI!     Broadcasts data from the master to all other nodes
!MPI!     with command MPI_BCAST
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWBROADC')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

!TIMG      CALL SWTSTA(201)
!MPI!NCOH      CALL MPI_BCAST ( IPTR, ILEN, ITYPE, MASTER-1,
!MPI!NCOH     &                 MPI_COMM_WORLD, IERR )
!MPI!COH      CALL MPI_BCAST ( IPTR, ILEN, ITYPE, MASTER-1,
!MPI!COH     &                 COMM, IERR )
!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
!MPI         CHARS = INTSTR(IERR)
!MPI         CALL TXPBLA(CHARS,IF,IL)
!MPI         MSGSTR = 'MPI produces some internal error - '//
!MPI     &            'return code is '//CHARS(IF:IL)
!MPI         CALL MSGERR ( 4, MSGSTR )
!MPI         RETURN
!MPI      END IF
!TIMG      CALL SWTSTO(201)

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWGATHER ( IOPTR, IOLEN, IIPTR, IILEN, ITYPE )
!
!****************************************************************
!
!MPI      USE MPI
      USE OCPCOMM4                                                        40.41
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Gathers different amounts of data from each processor
!     to the master
!
!  3. Method
!
!     Wrapper for MPI_GATHERV
!
!  4. Argument variables
!
!     IILEN       length of input array
!     IIPTR       pointer to first element of input array (local)
!     IOLEN       length of output array
!     IOPTR       pointer to first element of output array (global)
!     ITYPE       type of data
!
      INTEGER IILEN, IIPTR, IOLEN, IOPTR, ITYPE
!
!  5. Parameter variables
!
!     ---
!
!  6. Local variables
!
!     CHARS :     character for passing info to MSGERR
!     I     :     loop counter
!     ICOUNT:     array specifying array size of data received
!                 from each processor
!     IDSPLC:     array specifying the starting address of the
!                 incoming data from each processor, relative
!                 to the global array
!     IENT  :     number of entries
!     IERR  :     error value of MPI call
!     IF    :     first non-character in string
!     IL    :     last non-character in string
!     MSGSTR:     string to pass message to call MSGERR
!
      INTEGER                 I, IENT, IERR, IF, IL
      INTEGER, ALLOCATABLE :: ICOUNT(:), IDSPLC(:)
      CHARACTER*20            INTSTR, CHARS
      CHARACTER*80            MSGSTR
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!MPI!     MPI_GATHER       Gathers data from all nodes to the master
!MPI!     MPI_GATHERV      Gathers different amounts of data from
!MPI!                      all nodes to the master
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
! 10. Error messages
!
!     ---
!
! 12. Structure
!
!     if not parallel, return
!
!MPI!     gather the array sizes to the master
!MPI!
!     check whether enough space has been allocated
!     for gathered data
!
!     calculate starting address of each local array
!     with respect to the global array
!
!MPI!     gather different amounts of data from each processor
!MPI!     to the master
!MPI!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWGATHER')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

      IF (IAMMASTER) THEN
         ALLOCATE(ICOUNT(0:NPROC-1))
         ALLOCATE(IDSPLC(0:NPROC-1))
      END IF

!MPI!     --- gather the array sizes to the master
!MPI
!MPI!NCOH      CALL MPI_GATHER( IILEN, 1, SWINT, ICOUNT, 1, SWINT,
!MPI!NCOH     &                 MASTER-1, MPI_COMM_WORLD, IERR )
!MPI!COH      CALL MPI_GATHER( IILEN, 1, SWINT, ICOUNT, 1, SWINT,
!MPI!COH     &                 MASTER-1, COMM, IERR )
!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
!MPI         CHARS = INTSTR(IERR)
!MPI         CALL TXPBLA(CHARS,IF,IL)
!MPI         MSGSTR = 'MPI produces some internal error - '//
!MPI     &            'return code is '//CHARS(IF:IL)
!MPI         CALL MSGERR ( 4, MSGSTR )
!MPI         RETURN
!MPI      END IF

!     --- check whether enough space has been allocated
!         for gathered data

      IF (IAMMASTER) THEN
         IF ( SUM(ICOUNT).GT.IOLEN ) THEN
            CALL MSGERR(4,
     &                  'Not enough space allocated for gathered data')
            RETURN
         END IF
      END IF

!     --- calculate starting address of each local array
!         with respect to the global array

      IF (IAMMASTER) THEN
         IDSPLC(0) = 0
         DO I = 1, NPROC-1
            IDSPLC(I) = ICOUNT(I-1) + IDSPLC(I-1)
         END DO
      END IF

!MPI!     --- gather different amounts of data from each processor
!MPI!         to the master
!MPI
!MPI!NCOH      CALL MPI_GATHERV( IIPTR, IILEN, ITYPE, IOPTR, ICOUNT, IDSPLC,
!MPI!NCOH     &                  ITYPE, MASTER-1, MPI_COMM_WORLD, IERR )
!MPI!COH      CALL MPI_GATHERV( IIPTR, IILEN, ITYPE, IOPTR, ICOUNT, IDSPLC,
!MPI!COH     &                  ITYPE, MASTER-1, COMM, IERR )
!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
!MPI         CHARS = INTSTR(IERR)
!MPI         CALL TXPBLA(CHARS,IF,IL)
!MPI         MSGSTR = 'MPI produces some internal error - '//
!MPI     &            'return code is '//CHARS(IF:IL)
!MPI         CALL MSGERR ( 4, MSGSTR )
!MPI         RETURN
!MPI      END IF

      IF (IAMMASTER) DEALLOCATE(ICOUNT,IDSPLC)

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWREDUCE ( IPTR, ILEN, ITYPE, ITYPRD )
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.96: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.96, Dec. 08: call SWREDUCI/R instead of passing startaddress of the array
!
!  2. Purpose
!
!     Performs a global reduction of type ITYPRD on
!     array (IPTR) of type ITYPE to collect values from
!     all processes
!
!  4. Argument variables
!
!     ILEN        length of array to be collect
!     IPTR        pointer to first element of array to be collect
!     ITYPE       type of data
!     ITYPRD      type of reduction
!
      INTEGER IPTR, ILEN, ITYPE, ITYPRD
!
!  5. Parameter variables
!
!     ---
!
!  6. Local variables
!
!     IENT  :     number of entries
!
      INTEGER IENT
!
!  8. Subroutines used
!
!     STRACE           Tracing routine for debugging
!
!  9. Subroutines calling
!
!     SWCOMP
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     Performs a global reduction of data across all nodes
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWREDUCE')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN
!
!     --- actual reduction of field array based on its type
      IF ( ITYPE.EQ.SWINT ) THEN
         CALL SWREDUCI ( IPTR, ILEN, ITYPRD )
      ELSE IF ( ITYPE.EQ.SWREAL ) THEN
         CALL SWREDUCR ( IPTR, ILEN, ITYPRD )
      END IF

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWREDUCI ( IARR, ILEN, ITYPRD )
!
!****************************************************************
!
!MPI      USE MPI
      USE OCPCOMM4
      USE M_PARALL
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.96: Marcel Zijlema
!
!  1. Updates
!
!     40.96, Dec. 08: New subroutine
!
!  2. Purpose
!
!     Performs a global reduction of type ITYPRD on integer
!     array IARR to collect values from all processes
!
!  3. Method
!
!     Wrapper for MPI_ALLREDUCE
!
!  4. Argument variables
!
!     IARR        integer array
!     ILEN        length of array to be collect
!     ITYPRD      type of reduction
!
      INTEGER ILEN, ITYPRD
      INTEGER IARR(ILEN)
!
!  6. Local variables
!
!     CHARS :     character for passing info to MSGERR
!     IENT  :     number of entries
!     IERR  :     error value of MPI call
!     IF    :     first non-character in string
!     IL    :     last non-character in string
!     ITEMP :     temporary array to store collected data
!     MSGSTR:     string to pass message to call MSGERR
!
      INTEGER      IENT, IERR, IF, IL
      CHARACTER*20 INTSTR, CHARS
      CHARACTER*80 MSGSTR

      INTEGER, ALLOCATABLE :: ITEMP(:)
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!MPI!     MPI_ALLREDUCE    Combines values from all processes and
!MPI!                      distribute the result back to all processes
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!TIMG!     SWTSTA           Start timing for a section of code
!TIMG!     SWTSTO           Stop timing for a section of code
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SWREDUCE
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!MPI!     Performs a global reduction of integers across all nodes
!MPI!     with command MPI_ALLREDUCE
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWREDUCI')

      ALLOCATE(ITEMP(ILEN))

!TIMG      CALL SWTSTA(202)
!MPI!NCOH      CALL MPI_ALLREDUCE ( IARR, ITEMP, ILEN, SWINT,
!MPI!NCOH     &                     ITYPRD, MPI_COMM_WORLD, IERR )
!MPI!COH      CALL MPI_ALLREDUCE ( IARR, ITEMP, ILEN, SWINT,
!MPI!COH     &                     ITYPRD, COMM, IERR )
!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
!MPI         CHARS = INTSTR(IERR)
!MPI         CALL TXPBLA(CHARS,IF,IL)
!MPI         MSGSTR = 'MPI produces some internal error - '//
!MPI     &            'return code is '//CHARS(IF:IL)
!MPI         CALL MSGERR ( 4, MSGSTR )
!MPI         RETURN
!MPI      END IF
      IARR = ITEMP
!TIMG      CALL SWTSTO(202)

      DEALLOCATE(ITEMP)

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWREDUCR ( ARR, ILEN, ITYPRD )
!
!****************************************************************
!
!MPI      USE MPI
      USE OCPCOMM4
      USE M_PARALL
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.96: Marcel Zijlema
!
!  1. Updates
!
!     40.96, Dec. 08: New subroutine
!
!  2. Purpose
!
!     Performs a global reduction of type ITYPRD on real
!     array ARR to collect values from all processes
!
!  3. Method
!
!     Wrapper for MPI_ALLREDUCE
!
!  4. Argument variables
!
!     ARR         real array
!     ILEN        length of array to be collect
!     ITYPRD      type of reduction
!
      INTEGER ILEN, ITYPRD
      REAL    ARR(ILEN)
!
!  6. Local variables
!
!     CHARS :     character for passing info to MSGERR
!     IENT  :     number of entries
!     IERR  :     error value of MPI call
!     IF    :     first non-character in string
!     IL    :     last non-character in string
!     MSGSTR:     string to pass message to call MSGERR
!     TEMP  :     temporary array to store collected data
!
      INTEGER      IENT, IERR, IF, IL
      CHARACTER*20 INTSTR, CHARS
      CHARACTER*80 MSGSTR

      REAL, ALLOCATABLE :: TEMP(:)
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!MPI!     MPI_ALLREDUCE    Combines values from all processes and
!MPI!                      distribute the result back to all processes
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!TIMG!     SWTSTA           Start timing for a section of code
!TIMG!     SWTSTO           Stop timing for a section of code
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SWREDUCE
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!MPI!     Performs a global reduction of reals across all nodes
!MPI!     with command MPI_ALLREDUCE
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWREDUCR')

      ALLOCATE(TEMP(ILEN))

!TIMG      CALL SWTSTA(202)
!MPI!NCOH      CALL MPI_ALLREDUCE ( ARR, TEMP, ILEN, SWREAL,
!MPI!NCOH     &                     ITYPRD, MPI_COMM_WORLD, IERR )
!MPI!COH      CALL MPI_ALLREDUCE ( ARR, TEMP, ILEN, SWREAL,
!MPI!COH     &                     ITYPRD, COMM, IERR )
!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
!MPI         CHARS = INTSTR(IERR)
!MPI         CALL TXPBLA(CHARS,IF,IL)
!MPI         MSGSTR = 'MPI produces some internal error - '//
!MPI     &            'return code is '//CHARS(IF:IL)
!MPI         CALL MSGERR ( 4, MSGSTR )
!MPI         RETURN
!MPI      END IF
      ARR = TEMP
!TIMG      CALL SWTSTO(202)

      DEALLOCATE(TEMP)

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWSTRIP ( IPOWN, IDIR, NPART, IWORK, MXC, MYC )
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Performs a stripwise partitioning with straight interfaces
!
!  3. Method
!
!     Each active point in a row/column will be assign to a part
!     according to its number and size (stored in IWORK).
!     The remaining points in the row/column will be assign
!     to the same part.
!
!  4. Argument variables
!
!     IDIR        direction of cutting
!                 1 = row
!                 2 = column
!     IPOWN       array giving the subdomain number of each gridpoint
!     IWORK       work array with the following meaning:
!                    IWORK(1,i) = number of i-th part to be created
!                    IWORK(2,i) = size of i-th part to be created
!     MXC         maximum counter of gridpoints in x-direction
!     MYC         maximum counter of gridpoints in y-direction
!     NPART       number of parts to be created
!
      INTEGER   IDIR, MXC, MYC, NPART
      INTEGER   IPOWN(*)
      INTEGER*8 IWORK(2,*)
!
!  6. Local variables
!
!     IC    :     index of (IX,IY)-point
!     ICC   :     index of (IX,IY)-point
!     IENT  :     number of entries
!     INCX  :     increment for adressing: 1 for x-dir, MXC for y-dir
!     INCY  :     increment for adressing: MXC for x-dir, 1 for y-dir
!     IX    :     index in x-direction
!     IY    :     index in y-direction
!     IYY   :     index in y-direction
!     IPART :     a part counter
!     MXCI  :     maximum counter of gridpoints in x/y-direction
!     MYCI  :     maximum counter of gridpoints in y/x-direction
!     NCURPT:     number of currently assigned points to a created part
!     NPREM :     number of remaining points in a row/column
!
      INTEGER IC, ICC, IENT, INCX, INCY, IX, IY, IYY, IPART,
     &        MXCI, MYCI, NCURPT, NPREM
!
!  8. Subroutines used
!
!     STRACE           Tracing routine for debugging
!
!  9. Subroutines calling
!
!     SWPARTIT
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     depending on cutting direction, determine indirect addressing
!     create first empty part
!     for all active points do
!         assign this point to the created part
!         if size of created part has been reached
!            determine remaining active points in the current column
!            if no remaining points, create next empty part
!            else remaining points belong to the current part
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWSTRIP')

!     --- depending on cutting direction, determine indirect addressing
!         for array IPOWN

      IF ( IDIR.EQ.1 ) THEN
         MXCI = MYC
         MYCI = MXC
         INCX = MXC
         INCY = 1
      ELSE IF ( IDIR.EQ.2 ) THEN
         MXCI = MXC
         MYCI = MYC
         INCX = 1
         INCY = MXC
      END IF

!     --- create first empty part

      IPART  = 1
      NCURPT = 0

!     --- for all active points do

      DO IX = 1, MXCI
         DO IY = 1, MYCI

            IC = IX*INCX + IY*INCY - MXC

            IF ( IPOWN(IC).EQ.1 ) THEN

!              --- assign this point to the created part

               IPOWN(IC) = IWORK(1,IPART)
               NCURPT    = NCURPT + 1

!              --- if size of created part has been reached

               IF ( NCURPT.GE.IWORK(2,IPART) ) THEN

!                 --- determine remaining active points in the
!                     current column

                  NPREM = 0
                  DO IYY = IY+1, MYCI
                     ICC = IX*INCX + IYY*INCY - MXC
                     IF (IPOWN(ICC).EQ.1) NPREM = NPREM +1
                  END DO

                  IF ( NPREM.EQ.0 ) THEN

!                    --- if no remaining points, create next empty part

                     IPART  = IPART + 1
                     NCURPT = 0

                  ELSE

!                    --- else remaining points belong to the current part

                     IWORK(2,IPART  ) = IWORK(2,IPART  ) + NPREM
                     IWORK(2,IPART+1) = IWORK(2,IPART+1) - NPREM

                  END IF

               END IF

            END IF

         END DO
      END DO

      RETURN
      END
!JAC!****************************************************************
!JAC!
!JAC      SUBROUTINE SWSTRIP ( IPOWN, IDIR, IPART, NPART, LPARTS,
!JAC     &                     MXC  , MYC )
!JAC!
!JAC!****************************************************************
!JAC!
!JAC      USE OCPCOMM4                                                        40.41
!JAC!
!JAC      IMPLICIT NONE
!JAC!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!JAC!
!JAC!  0. Authors
!JAC!
!JAC!     40.30: Marcel Zijlema
!JAC!     40.41: Marcel Zijlema
!JAC!
!JAC!  1. Updates
!JAC!
!JAC!     40.30, Feb. 03: New subroutine
!JAC!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!JAC!
!JAC!  2. Purpose
!JAC!
!JAC!     Performs a stripwise partitioning with straight interfaces
!JAC!
!JAC!  3. Method
!JAC!
!JAC!     This method is described in Ph.D. Thesis of M. Roest
!JAC!     entitled:
!JAC!     Partitioning for parallel finite difference computations
!JAC!     in coastal water simulation, DUT, 1997
!JAC!
!JAC!  4. Argument variables
!JAC!
!JAC!     IDIR        direction of cutting
!JAC!                 1 = row
!JAC!                 2 = column
!JAC!     IPART       part number that must be partitioned
!JAC!     IPOWN       array giving the subdomain number of each gridpoint
!JAC!     LPARTS      list of parts to be created
!JAC!                    lparts(1,i) = number of i-th part to be created
!JAC!                    lparts(2,i) = size of i-th part to be created
!JAC!     MXC         maximum counter of gridpoints in x-direction
!JAC!     MYC         maximum counter of gridpoints in y-direction
!JAC!     NPART       number of parts to be created
!JAC!
!JAC      INTEGER   IDIR, IPART, MXC, MYC, NPART
!JAC      INTEGER   IPOWN(*)
!JAC      INTEGER*8 LPARTS(2,*)
!JAC!
!JAC!  6. Local variables
!JAC!
!JAC!     IC    :     index of (IX,IY)-point
!JAC!     ICC   :     index of (IX,IY)-point
!JAC!     IENT  :     number of entries
!JAC!     INCX  :     increment for adressing: 1 for x-dir, MXC for y-dir
!JAC!     INCY  :     increment for adressing: MXC for x-dir, 1 for y-dir
!JAC!     IX    :     index in x-direction
!JAC!     IY    :     index in y-direction
!JAC!     IYY   :     index in y-direction
!JAC!     JPART :     a part counter
!JAC!     MXCI  :     maximum counter of gridpoints in x/y-direction
!JAC!     MYCI  :     maximum counter of gridpoints in y/x-direction
!JAC!     NBACK :     number of points in a row already assigned to a part
!JAC!     NFORW :     number of points in a row remaining to be assigned
!JAC!     NINPRT:     number of points currently assigned to a new part
!JAC!
!JAC      INTEGER IC, ICC, IENT, INCX, INCY, IX, IY, IYY, JPART,
!JAC     &        MXCI, MYCI, NBACK, NFORW, NINPRT
!JAC!
!JAC!  8. Subroutines used
!JAC!
!JAC!     STRACE           Tracing routine for debugging
!JAC!
!JAC!  9. Subroutines calling
!JAC!
!JAC!     SWORB
!JAC!     SWPARTIT
!JAC!
!JAC! 10. Error messages
!JAC!
!JAC!     ---
!JAC!
!JAC! 11. Remarks
!JAC!
!JAC!     ---
!JAC!
!JAC! 12. Structure
!JAC!
!JAC!     depending on IDIR, determine indirect addressing in IPOWN
!JAC!     start by creating the first part, which is currently empty
!JAC!     for all points in IPOWN do
!JAC!        if point belongs to the part that must be partitioned then
!JAC!           when current part has reached its planned size
!JAC!              see how many points in this row remain to be assigned
!JAC!              if no more points to be assigned, go on to next part
!JAC!              else, if majority of row has been assigned take the rest
!JAC!              else, leave this row to next part
!JAC!           assign point IC to part that is currently being created
!JAC!
!JAC! 13. Source text
!JAC!
!JAC      SAVE IENT
!JAC      DATA IENT/0/
!JAC      IF (LTRACE) CALL STRACE (IENT,'SWSTRIP')
!JAC
!JAC!     --- depending on IDIR, determine indirect addressing in IPOWN
!JAC
!JAC      IF ( IDIR.EQ.1 ) THEN
!JAC         MXCI = MYC
!JAC         MYCI = MXC
!JAC         INCX = MXC
!JAC         INCY = 1
!JAC      ELSE IF ( IDIR.EQ.2 ) THEN
!JAC         MXCI = MXC
!JAC         MYCI = MYC
!JAC         INCX = 1
!JAC         INCY = MXC
!JAC      END IF
!JAC
!JAC!     --- start by creating the first part, which is currently empty
!JAC
!JAC      JPART  = 1
!JAC      NINPRT = 0
!JAC
!JAC!     --- for all points in IPOWN do
!JAC
!JAC      DO IX = 1, MXCI
!JAC         DO IY = 1, MYCI
!JAC
!JAC            IC = IX*INCX + IY*INCY - MXC
!JAC
!JAC!           --- if this point belongs to the part that must be partitioned then
!JAC
!JAC            IF ( IPOWN(IC).EQ.IPART ) THEN
!JAC
!JAC!              --- when current part has reached its planned size
!JAC
!JAC               IF ( NINPRT.GE.LPARTS(2,JPART) ) THEN
!JAC
!JAC!                 --- see how many points in this row have been assigned
!JAC
!JAC                  NBACK = 0
!JAC                  DO IYY = 1, IY-1
!JAC
!JAC                     ICC = IX*INCX + IYY*INCY - MXC
!JAC                     IF (IPOWN(ICC).EQ.LPARTS(1,JPART)) NBACK = NBACK +1
!JAC
!JAC                  END DO
!JAC
!JAC!                 --- see how many points in this row remain to be assigned
!JAC
!JAC                  NFORW = 0
!JAC                  DO IYY = IY, MYCI
!JAC
!JAC                     ICC = IX*INCX + IYY*INCY - MXC
!JAC                     IF (IPOWN(ICC).EQ.IPART) NFORW = NFORW +1
!JAC
!JAC                  END DO
!JAC
!JAC!                 --- if no more points to be assigned, go on to next part
!JAC
!JAC                  IF ( NFORW.EQ.0 ) THEN
!JAC
!JAC                     JPART  = JPART + 1
!JAC                     NINPRT = 0
!JAC
!JAC                  ELSE IF ( (NBACK-NFORW).GT.0 ) THEN
!JAC
!JAC!                    --- if majority of row has been assigned take the rest
!JAC
!JAC                     LPARTS(2,JPART  ) = LPARTS(2,JPART  ) + NFORW
!JAC                     LPARTS(2,JPART+1) = LPARTS(2,JPART+1) - NFORW
!JAC
!JAC                  ELSE
!JAC!                    --- else, leave this row to next part
!JAC
!JAC                     LPARTS(2,JPART  ) = LPARTS(2,JPART  ) - NBACK
!JAC                     LPARTS(2,JPART+1) = LPARTS(2,JPART+1) + NBACK
!JAC
!JAC                     DO IYY = 1, IY-1
!JAC
!JAC                        ICC = IX*INCX + IYY*INCY - MXC
!JAC                        IF ( IPOWN(ICC).EQ.LPARTS(1,JPART) ) THEN
!JAC                           IPOWN(ICC) = LPARTS(1,JPART+1)
!JAC                        END IF
!JAC
!JAC                     END DO
!JAC
!JAC                     JPART  = JPART + 1
!JAC                     NINPRT = NBACK
!JAC
!JAC                  END IF
!JAC
!JAC               END IF
!JAC
!JAC!              --- assign point IC to part that is currently being created
!JAC
!JAC               IPOWN(IC) = LPARTS(1,JPART)
!JAC               NINPRT    = NINPRT + 1
!JAC
!JAC            END IF
!JAC
!JAC         END DO
!JAC      END DO
!JAC
!JAC      RETURN
!JAC      END
!JAC!****************************************************************
!JAC!
!JAC      SUBROUTINE SWORB ( IPOWN, IDIR, IPART, NPART, LPARTS,
!JAC     &                   MXC  , MYC )
!JAC!
!JAC!****************************************************************
!JAC!
!JAC      USE OCPCOMM4                                                        40.41
!JAC!
!JAC      IMPLICIT NONE
!JAC!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!JAC!
!JAC!  0. Authors
!JAC!
!JAC!     40.30: Marcel Zijlema
!JAC!     40.41: Marcel Zijlema
!JAC!
!JAC!  1. Updates
!JAC!
!JAC!     40.30, Feb. 03: New subroutine
!JAC!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!JAC!
!JAC!  2. Purpose
!JAC!
!JAC!     Performs an Orthogonal Recursive Bisection partitioning
!JAC!
!JAC!  3. Method
!JAC!
!JAC!     Starting with a single part (the entire domain), each part is
!JAC!     recursively partitioned by bisecting it, until all parts have been
!JAC!     created. The bisection direction is swapped in each direction.
!JAC!
!JAC!     This method is described in Ph.D. Thesis of M. Roest
!JAC!     entitled:
!JAC!     Partitioning for parallel finite difference computations
!JAC!     in coastal water simulation, DUT, 1997
!JAC!
!JAC!  4. Argument variables
!JAC!
!JAC!     IDIR        direction of cutting
!JAC!                 1 = row
!JAC!                 2 = column
!JAC!     IPART       part number that must be partitioned
!JAC!     IPOWN       array giving the subdomain number of each gridpoint
!JAC!     LPARTS      list of parts to be created
!JAC!                    lparts(1,i) = number of i-th part to be created
!JAC!                    lparts(2,i) = size of i-th part to be created
!JAC!     MXC         maximum counter of gridpoints in x-direction
!JAC!     MYC         maximum counter of gridpoints in y-direction
!JAC!     NPART       number of parts to be created
!JAC!
!JAC      INTEGER   IDIR, IPART, MXC, MYC, NPART
!JAC      INTEGER   IPOWN(*)
!JAC      INTEGER*8 LPARTS(2,*)
!JAC!
!JAC!  6. Local variables
!JAC!
!JAC!     IDIFF :     the difference to be applied to a subdomain-size
!JAC!     IENT  :     number of entries
!JAC!     IP    :     counter of parts to be splitted
!JAC!     ISPLIT:     counter of parts to be created in splitting
!JAC!     ISSUCC:     flag indicating success in reducing a difference in size
!JAC!                 0=no
!JAC!                 1=yes
!JAC!     IWORK :     see description LPARTS
!JAC!     J     :     loop counter
!JAC!     JEND  :     number of last new part to be created by splitting
!JAC!     JPARTE:     number of last part in 1..npart belonging to jpart
!JAC!     JPARTS:     number of first part in 1..npart belonging to jpart
!JAC!     JSTART:     number of first new part to be created by splitting
!JAC!     KSPLIT:     number of parts to be created in a particular splitting
!JAC!     NP    :     number of parts to be created in a particular recursion
!JAC!     NSPLIT:     maximum number of parts to be created in one splitting
!JAC!                 (NB: 2 = bisection, 4 = quadrisection)
!JAC
!JAC      INTEGER IDIFF, IENT, IP, ISPLIT, ISSUCC, J, JEND, JPARTE, JPARTS,
!JAC     &        JSTART, KSPLIT, NP, NSPLIT
!JAC      INTEGER*8 IWORK(2,NPART)
!JAC      PARAMETER(NSPLIT = 2)
!JAC!
!JAC!  8. Subroutines used
!JAC!
!JAC!     MSGERR           Writes error message
!JAC!     STRACE           Tracing routine for debugging
!JAC!     SWSTRIP          Performs a stripwise partitioning with straight
!JAC!                      interfaces
!JAC!
!JAC!  9. Subroutines calling
!JAC!
!JAC!     SWPARTIT
!JAC!
!JAC! 10. Error messages
!JAC!
!JAC!     ---
!JAC!
!JAC! 11. Remarks
!JAC!
!JAC!     ---
!JAC!
!JAC! 12. Structure
!JAC!
!JAC!     while not enough parts have been created, do another recursion
!JAC!
!JAC!        for each part that currently exists
!JAC!
!JAC!          determine which final parts belong to this part
!JAC!          if the number of such parts > 1, do further splitting
!JAC!
!JAC!            determine into how many parts this part must be split
!JAC!            determine the sizes and numbers of parts to be created
!JAC!
!JAC!            do splitting
!JAC!
!JAC!            determine whether objective partsizes have been modified
!JAC!            and distribute the difference over the constituent parts
!JAC!
!JAC!        swap cutting direction
!JAC!
!JAC! 13. Source text
!JAC!
!JAC      SAVE IENT
!JAC      DATA IENT/0/
!JAC      IF (LTRACE) CALL STRACE (IENT,'SWORB')
!JAC
!JAC!     --- while not enough parts have been created, do another recursion
!JAC
!JAC      NP = 1
!JAC100   IF ( NP.LE.NPART ) THEN
!JAC
!JAC!        --- for each part that currently exists
!JAC
!JAC         DO IP = 1, NP
!JAC
!JAC!           --- determine which final parts belong to this part
!JAC
!JAC            JPARTS = (IP-1)*NPART/NP+1
!JAC            JPARTE = (IP  )*NPART/NP
!JAC
!JAC!           --- if the number of such parts > 1, do further splitting
!JAC
!JAC            IF ( (JPARTE-JPARTS+1).GT.1 ) THEN
!JAC
!JAC!              --- determine into how many parts this part must be split
!JAC
!JAC               KSPLIT = MIN(NSPLIT,JPARTE-JPARTS+1)
!JAC
!JAC!              --- determine the sizes and numbers of parts to be created
!JAC
!JAC               DO ISPLIT = 1, KSPLIT
!JAC
!JAC                  JSTART = JPARTS+(ISPLIT-1)*(JPARTE-JPARTS+1)/KSPLIT
!JAC                  JEND   = JPARTS+    ISPLIT*(JPARTE-JPARTS+1)/KSPLIT-1
!JAC
!JAC                  IWORK(1,ISPLIT) = JSTART
!JAC
!JAC                  IWORK(2,ISPLIT) = 0
!JAC                  DO J = JSTART, JEND
!JAC                     IWORK(2,ISPLIT) = IWORK(2,ISPLIT) + LPARTS(2,J)
!JAC                  END DO
!JAC
!JAC               END DO
!JAC
!JAC!              --- do splitting
!JAC
!JAC               CALL SWSTRIP (IPOWN,IDIR,JPARTS,KSPLIT,IWORK,MXC,MYC)
!JAC
!JAC!              --- determine whether objective partsizes have been modified
!JAC!                  in SWSTRIP in order to make straight interfaces
!JAC
!JAC               DO ISPLIT = 1, KSPLIT
!JAC
!JAC                  JSTART = JPARTS+(ISPLIT-1)*(JPARTE-JPARTS+1)/KSPLIT
!JAC                  JEND   = JPARTS+    ISPLIT*(JPARTE-JPARTS+1)/KSPLIT-1
!JAC
!JAC                  DO J = JSTART, JEND
!JAC                     IWORK(2,ISPLIT) = IWORK(2,ISPLIT) - LPARTS(2,J)
!JAC                  END DO
!JAC
!JAC!                 --- and distribute the difference over the contiguous
!JAC!                     parts making sure not to cause negative subdomain-sizes
!JAC
!JAC                  J      = JSTART
!JAC                  ISSUCC = 0
!JAC                  IF ( IWORK(2,ISPLIT).LT.0 ) THEN
!JAC                     IDIFF = -1
!JAC                  ELSE
!JAC                     IDIFF =  1
!JAC                  END IF
!JAC
!JAC!                 --- reduce the difference until nothing is left
!JAC
!JAC200               IF ( IWORK(2,ISPLIT).NE.0 ) THEN
!JAC
!JAC!                    --- only adjust parts if their size remains valid
!JAC
!JAC                     IF ( LPARTS(2,J).GT.0 .AND.
!JAC     &                    ((LPARTS(2,J)+IDIFF).GT.0) ) THEN
!JAC                        ISSUCC          = 1
!JAC                        LPARTS(2,J)     = LPARTS(2,J)     + IDIFF
!JAC                        IWORK(2,ISPLIT) = IWORK(2,ISPLIT) - IDIFF
!JAC                     END IF
!JAC
!JAC!                    --- go on to the next part
!JAC
!JAC                     J = J + 1
!JAC
!JAC!                    --- when all parts have been visited, go back to first
!JAC
!JAC                     IF ( J.GT.JEND ) THEN
!JAC
!JAC!                       --- check whether any reduction of the difference
!JAC!                           was done in last pass over all parts
!JAC
!JAC                        IF ( ISSUCC.EQ.0 ) THEN
!JAC                           CALL MSGERR (4,'Internal problem in SWORB')
!JAC                           RETURN
!JAC                        END IF
!JAC                        J = JSTART
!JAC                     END IF
!JAC
!JAC                     GOTO 200
!JAC                  END IF
!JAC
!JAC               END DO
!JAC
!JAC            END IF
!JAC
!JAC         END DO
!JAC
!JAC!        --- swap cutting direction
!JAC
!JAC         IDIR = MOD(IDIR,2) + 1
!JAC
!JAC         NP = NSPLIT*NP
!JAC         GOTO 100
!JAC
!JAC      END IF
!JAC
!JAC      RETURN
!JAC      END
!****************************************************************
!
      SUBROUTINE SWPARTIT ( IPOWN, MXC, MYC )
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!     40.41, Sep. 04: determines load per processor based on speed
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Carries out the partitioning of the SWAN computational grid
!
!  3. Method
!
!     Based on stripwise partitioning
!JAC!     Based on Orthogonal Recursive Bisection
!
!  4. Argument variables
!
!     IPOWN       array giving the subdomain number of each gridpoint
!     MXC         maximum counter of gridpoints in x-direction
!     MYC         maximum counter of gridpoints in y-direction
!
      INTEGER MXC, MYC
      INTEGER IPOWN(MXC,MYC)
!
!  6. Local variables
!
!     I     :     loop counter
!     ICNT  :     auxiliary integer to count weights
!     IDIR  :     direction of cutting
!                 1 = row
!                 2 = column
!     IENT  :     number of entries
!     IX    :     index in x-direction
!     IY    :     index in y-direction
!     IWORK :     work array with the following meaning:
!                    IWORK(1,i) = number of i-th part to be created
!                    IWORK(2,i) = size of i-th part to be created
!     NACTP :     total number of active gridpoints
!     NPCUM :     cumulative number of gridpoints
!
      INTEGER   I, IDIR, IENT, IX, IY
      INTEGER*8 ICNT, NACTP, NPCUM
      INTEGER*8 IWORK(2,NPROC)
!
!  8. Subroutines used
!
!     STPNOW           Logical indicating whether program must
!                      terminated or not
!     STRACE           Tracing routine for debugging
!JAC!     SWORB            Performs an Orthogonal Recursive Bisection partitioning
!     SWSTRIP          Performs a stripwise partitioning with straight
!                      interfaces
!
      LOGICAL STPNOW
!
!  9. Subroutines calling
!
!     SWDECOMP
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     if not parallel, return
!     determine direction of cutting
!     determine number of active points
!     determine numbers and sizes of parts to be created
!     partition grid
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWPARTIT')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

!     --- determine direction of cutting

      IF ( MXC.GT.MYC ) THEN
         IDIR = 2
      ELSE
         IDIR = 1
      END IF

!     --- determine number of active points and
!         set IPOWN to 1 in these points

      NACTP = 0
      DO IX = 1, MXC
         DO IY = 1, MYC
            IF ( KGRPGL(IX,IY).NE.1 ) THEN
               IPOWN(IX,IY) = 1
               NACTP        = NACTP + 1
            END IF
         END DO
      END DO

!     --- determine numbers and sizes of parts to be created

      NPCUM = 0
      ICNT  = 0
      DO I = 1, NPROC
         ICNT       = ICNT + IWEIG(I)
         IWORK(1,I) = I
         IWORK(2,I) = (NACTP*ICNT)/SUM(IWEIG) - NPCUM
         NPCUM      = (NACTP*ICNT)/SUM(IWEIG)
      END DO
      DEALLOCATE(IWEIG)

!     --- partition grid

      CALL SWSTRIP ( IPOWN, IDIR, NPROC, IWORK, MXC, MYC )
!JAC      CALL SWORB ( IPOWN, IDIR, 1, NPROC, IWORK, MXC, MYC )
!JAC      IF (STPNOW()) RETURN

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWBLADM ( IPOWN, MXC, MYC )
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!     40.41, Jul. 04: determine global bounds in subdomains
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     For the present node, carries out the block administration
!     and determines array bounds with respect to global grid
!
!  3. Method
!
!     Based on domain decomposition, the interface sizes are
!     determined that is needed for the setup of block
!     administration stored as IBLKAD
!
!  4. Argument variables
!
!     IPOWN       array giving the subdomain number of each gridpoint
!     MXC         maximum counter of gridpoints in x-direction
!     MYC         maximum counter of gridpoints in y-direction
!
      INTEGER MXC, MYC
      INTEGER IPOWN(MXC,MYC)
!
!  5. Parameter variables
!
!     ---
!
!  6. Local variables
!
!     CHARS :     character for passing info to MSGERR
!     I     :     loop counter
!     IC    :     index of (IX,IY)-point
!     ICOFF :     offset of IC-index
!     ICRECV:     array containing positions of unknowns
!                 to be received from neighbour
!     ICSEND:     array containing positions of unknowns
!                 to be sent to neighbour
!     IDOM  :     subdomain number
!     IENT  :     number of entries
!     IF    :     first non-character in string
!     IL    :     last non-character in string
!     INB   :     neighbour counter
!     ISTART:     startaddress for each size interface in array IBLKAD
!     IX    :     index in x-direction
!     IXOFF :     offset in x-direction
!     IY    :     index in y-direction
!     IYOFF :     offset in y-direction
!     IWORK :     array used to determine interface sizes
!                   IWORK(1,i) = number of the i-th neighbour
!                   IWORK(2,i) = position of the i-th neighbour with
!                                respect to present subdomain
!                                (resp. top, bottom, right, left)
!                   IWORK(3,i) = size of interface to i-th neighbour
!     JOFFS :     offsets at which a point of a neigbhour domain can be found
!     MSGSTR:     string to pass message to call MSGERR
!     MXSIZ :     size of present subdomain in x-direction
!     MYSIZ :     size of present subdomain in y-direction
!     NNEIGH:     number of neighbouring subdomains
!     NOVLU :     number of overlapping unknowns
!
      INTEGER      I, IC, ICOFF, IDOM, IENT, IF, IL, INB, ISTART,
     &             IX, IXOFF, IY, IYOFF, JOFFS(2,4), MXSIZ, MYSIZ,
     &             NNEIGH, NOVLU
      INTEGER      IWORK(3,NPROC),
     &             ICRECV(NPROC,MAX(MXC,MYC)),
     &             ICSEND(NPROC,MAX(MXC,MYC))
      CHARACTER*20 INTSTR, CHARS
      CHARACTER*80 MSGSTR
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SWDECOMP
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     intialize offsets to be used in searching for interfaces
!     determine enclosing box of present subdomain
!     if subdomain appears to be empty
!        give warning and set empty bounding box
!     else
!        extend enclosing box to include halo area
!     localize global bounds in present subdomain
!     determine size of enclosing box
!     determine interface sizes:
!
!        loop over global grid
!           if point belongs to this part
!              for each of the four sizes
!                  if a neighbouring subdomain is found there
!                     find it in the list of neighbours
!                     if not yet in the list, add it
!                     store position of neighbour
!                     update number of overlapping unknowns
!
!     store block administration
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWBLADM')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

!     --- intialize offsets to be used in searching for interfaces

      JOFFS = RESHAPE((/0,1,0,-1,1,0,-1,0/), (/2,4/))

!     --- determine enclosing box of present subdomain

      IF ( MXC.GT.MYC ) THEN
         MXF = MXC+1
         MXL = 0
      ELSE
         MYF = MYC+1
         MYL = 0
      END IF
!JAC      MXF = MXC+1
!JAC      MYF = MYC+1
!JAC      MXL = 0
!JAC      MYL = 0

      DO IX = 1, MXC
         DO IY = 1, MYC

            IF( IPOWN(IX,IY).EQ.INODE ) THEN

               MXF = MIN(IX,MXF)
               MYF = MIN(IY,MYF)
               MXL = MAX(IX,MXL)
               MYL = MAX(IY,MYL)

            END IF

         END DO
      END DO

!     --- if subdomain appears to be empty

      IF ( MXF.GT.MXL .OR. MYF.GT.MYL ) THEN

!        --- give warning and set empty bounding box

         CHARS = INTSTR(INODE)
         CALL TXPBLA(CHARS,IF,IL)
         MSGSTR = 'Empty subdomain is detected - '//
     &            ' node number is '//CHARS(IF:IL)
         CALL MSGERR ( 2, MSGSTR )

         MXF = 1
         MYF = 1
         MXL = 0
         MYL = 0

      ELSE

!        --- extend enclosing box to include halo area

         MXF = MAX(1  ,MXF-IHALOX)
         MYF = MAX(1  ,MYF-IHALOY)
         MXL = MIN(MXC,MXL+IHALOX)
         MYL = MIN(MYC,MYL+IHALOY)

      END IF

!     --- localize global bounds in present subdomain

      IF ( MXCGL.GT.MYCGL ) THEN
         LMXF = MXF.EQ.1     .AND. INODE.EQ.1
         LMXL = MXL.EQ.MXCGL .AND. INODE.EQ.NPROC
         LMYF = MYF.EQ.1
         LMYL = MYL.EQ.MYCGL
      ELSE
         LMXF = MXF.EQ.1
         LMXL = MXL.EQ.MXCGL
         LMYF = MYF.EQ.1     .AND. INODE.EQ.1
         LMYL = MYL.EQ.MYCGL .AND. INODE.EQ.NPROC
      END IF

!     --- determine size of enclosing box

      MXSIZ = MXL - MXF + 1
      MYSIZ = MYL - MYF + 1

      IWORK  = 0
      ICRECV = 0
      ICSEND = 0

!     --- determine interface sizes

      DO IX = 1, MXC
         DO IY = 1, MYC

!           --- if point belongs to this part

            IF ( IPOWN(IX,IY).EQ.INODE ) THEN

!              --- for each of the four sizes

               DO I = 1, 4

                  IXOFF = JOFFS(1,I)
                  IYOFF = JOFFS(2,I)

!                 --- if a neighbouring subdomain is found there

                  IF ( (IX+IXOFF).GT.0.AND.(IX+IXOFF).LE.MXC.AND.
     &                 (IY+IYOFF).GT.0.AND.(IY+IYOFF).LE.MYC ) THEN

                     IF ( IPOWN(IX+IXOFF,IY+IYOFF).NE.0.AND.
     &                    IPOWN(IX+IXOFF,IY+IYOFF).NE.INODE ) THEN

                        IC    = (IY      -MYF)*MXSIZ + (IX      -MXF+1)
                        ICOFF = (IY+IYOFF-MYF)*MXSIZ + (IX+IXOFF-MXF+1)

!                       --- find it in the list of neighbours

                        IDOM = IPOWN(IX+IXOFF,IY+IYOFF)

                        INB = 1
  100                   IF ( INB.LE.NPROC .AND.
     &                       IWORK(1,INB).NE.IDOM .AND.
     &                       IWORK(1,INB).NE.0 )  THEN
                           INB = INB + 1
                           GOTO 100
                        END IF

                        IF ( INB.GT.NPROC ) THEN
                          CALL MSGERR (4,'Found more neighbours than '//
     &                                 'subdomains in the partitioning')
                          RETURN
                        END IF

!                       --- if not yet in the list, add it

                        IF ( IWORK(1,INB).EQ.0 ) IWORK(1,INB) = IDOM

!                       --- store position of neighbour with respect to
!                           present subdomain

                        IWORK(2,INB) = I

!                       --- update number of overlapping unknowns

                        IWORK(3,INB) = IWORK(3,INB) + 1

                        ICSEND(INB,IWORK(3,INB)) = IC
                        ICRECV(INB,IWORK(3,INB)) = ICOFF

                     END IF

                  END IF

               END DO

            END IF

         END DO
      END DO

!     --- store block administration

      NNEIGH    = COUNT(IWORK(1,:)>0)
      IBLKAD(1) = NNEIGH
      ISTART    = 3*NNEIGH+2
      DO INB = 1, NNEIGH
         IBLKAD(3*INB-1) = IWORK(1,INB)
         IBLKAD(3*INB  ) = IWORK(2,INB)
         IBLKAD(3*INB+1) = ISTART
         NOVLU           = IWORK(3,INB)
         IBLKAD(ISTART)  = NOVLU
         DO I = 1, NOVLU
            IBLKAD(ISTART      +I) = ICSEND(INB,I)
            IBLKAD(ISTART+NOVLU+I) = ICRECV(INB,I)
         END DO
         ISTART = ISTART + 2*NOVLU+1
      END DO

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWDECOMP
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SWCOMM3                                                         40.41
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Carries out domain decomposition meant for
!     distributed-memory approach
!
!  3. Method
!
!     First, carry out the partitioning of the
!     SWAN computational grid and then do the
!     block administration
!
!  4. Argument variables
!
!     ---
!
!  6. Local variables
!
!     IENT  :     number of entries
!     IX    :     loop counter
!     IY    :     loop counter
!     IPOWN :     array giving the subdomain number of each gridpoint
!
      INTEGER IENT, IX, IY
      INTEGER, ALLOCATABLE :: IPOWN(:,:)
!
!  8. Subroutines used
!
!     STPNOW           Logical indicating whether program must
!                      terminated or not
!     STRACE           Tracing routine for debugging
!     SWBLADM          Carries out the block administration
!     SWPARTIT         Carries out the partitioning of the SWAN
!                      computational grid
!
      LOGICAL STPNOW
!
!  9. Subroutines calling
!
!     SWMAIN
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     allocate and initialize array for block administration
!     store the original values of MXC, MYC and MCGRD
!     if not parallel, return
!     carry out the partitioning of computational grid
!     carry out the block administration
!     compute MXC, MYC and MCGRD for each subdomain
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWDECOMP')

!     --- allocate and initialize array for block administration

      IF (.NOT.ALLOCATED(IBLKAD)) ALLOCATE(IBLKAD(41+20*MAX(MXC,MYC)))
      IBLKAD = 0

!     --- store the original values of MXC, MYC and MCGRD of
!         global computational grid

      MXCGL   = MXC
      MYCGL   = MYC
      MCGRDGL = MCGRD
      MXF     = 1
      MXL     = MXC
      MYF     = 1
      MYL     = MYC
      LMXF    = MXF.EQ.1                                                  40.41
      LMXL    = MXL.EQ.MXCGL                                              40.41
      LMYF    = MYF.EQ.1                                                  40.41
      LMYL    = MYL.EQ.MYCGL                                              40.41

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

!     --- carry out the partitioning of the SWAN
!         computational grid

      ALLOCATE(IPOWN(MXC,MYC))
      IPOWN = 0
      CALL SWPARTIT( IPOWN, MXC, MYC )
      IF (STPNOW()) RETURN

!     --- carry out the block administration

      CALL SWBLADM( IPOWN, MXC, MYC )
      IF (STPNOW()) RETURN

!     --- compute MXC, MYC and MCGRD for each subdomain

      MXC = MXL - MXF + 1
      MYC = MYL - MYF + 1

      MCGRD = 1
      DO IX = MXF, MXL
         DO IY = MYF, MYL
            IF ( KGRPGL(IX,IY).NE.1 ) MCGRD = MCGRD + 1
         END DO
      END DO

      DEALLOCATE(IPOWN)

      RETURN
      END
!JAC!****************************************************************
!JAC!
!JAC      SUBROUTINE SWEXCHG ( FIELD, SWPDIR, KGRPNT )
!JAC!
!JAC!****************************************************************
!JAC!
!JAC      USE OCPCOMM4                                                        40.41
!JAC      USE SWCOMM3                                                         40.41
!JAC      USE M_PARALL                                                        40.31
!JAC!
!JAC      IMPLICIT NONE
!JAC!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!JAC!
!JAC!  0. Authors
!JAC!
!JAC!     40.30: Marcel Zijlema
!JAC!     40.41: Marcel Zijlema
!JAC!
!JAC!  1. Updates
!JAC!
!JAC!     40.30, Feb. 03: New subroutine
!JAC!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!JAC!
!JAC!  2. Purpose
!JAC!
!JAC!     Updates geographical field array through exchanging
!JAC!     values between neighbouring subdomains depending on
!JAC!     sweep direction
!JAC!
!JAC!  3. Method
!JAC!
!JAC!     Made use of MPI by means of SWSENDNB and SWRECVNB
!JAC!     and also block administration (stored in IBLKAD)
!JAC!
!JAC!  4. Argument variables
!JAC!
!JAC!     FIELD       geographical field array for which 'halo' values must
!JAC!                 be copied from neighbouring subdomains
!JAC!     KGRPNT      indirect addressing for grid points
!JAC!     SWPDIR      sweep direction (0=all directions together)
!JAC!
!JAC      INTEGER SWPDIR
!JAC      INTEGER KGRPNT(MXC*MYC)
!JAC      REAL    FIELD(MCGRD)
!JAC!
!JAC!  6. Local variables
!JAC!
!JAC!     IDOM  :     subdomain number
!JAC!     IENT  :     number of entries
!JAC!     INB   :     neighbour counter
!JAC!     IPNB  :     position of neighbour (=top, bottom, right, left)
!JAC!     IPR   :     array containing positions of neighbours from
!JAC!                 which data is to be received
!JAC!     IPS   :     array containing positions of neighbours to
!JAC!                 which data is to be sent
!JAC!     ISTART:     pointer in array IBLKAD
!JAC!     ISWP  :     sweep direction
!JAC!     ITAG  :     message tag for sending and receiving
!JAC!     K     :     loop counter
!JAC!     NNEIGH:     number of neighbouring subdomains
!JAC!     NOVLU :     number of overlapping unknowns
!JAC!     WORK  :     work array to store data to be sent to or
!JAC!                 received from neighbour
!JAC!
!JAC      INTEGER IDOM, IENT, INB, IPNB, ISTART, ISWP, ITAG,
!JAC     &        K, NNEIGH, NOVLU
!JAC      INTEGER IPR(2,4), IPS(2,4)
!JAC      REAL    WORK(MAX(MXC,MYC))
!JAC!
!JAC!  8. Subroutines used
!JAC!
!JAC!     STPNOW           Logical indicating whether program must
!JAC!                      terminated or not
!JAC!     STRACE           Tracing routine for debugging
!JAC!     SWRECVNB         Data is received from a neighbour
!JAC!     SWSENDNB         Data is sent to a neighbour
!JAC!TIMG!     SWTSTA           Start timing for a section of code
!JAC!TIMG!     SWTSTO           Stop timing for a section of code
!JAC!
!JAC      LOGICAL STPNOW
!JAC!
!JAC!  9. Subroutines calling
!JAC!
!JAC!     SWCOMP
!JAC!
!JAC! 10. Error messages
!JAC!
!JAC!     ---
!JAC!
!JAC! 11. Remarks
!JAC!
!JAC!     ---
!JAC!
!JAC! 12. Structure
!JAC!
!JAC!     if not parallel, return
!JAC!
!JAC!     for all neighbouring subdomains do
!JAC!        get position
!JAC!        if position corresponds to sweep selection
!JAC!           get subdomain number, pointer and size
!JAC!           store data to be sent in array WORK
!JAC!           send array WORK
!JAC!
!JAC!     for all neighbouring subdomains do
!JAC!        get position
!JAC!        if position corresponds to sweep selection
!JAC!           get subdomain number, pointer and size
!JAC!           receive next array and store in WORK
!JAC!           store the received data
!JAC!
!JAC! 13. Source text
!JAC!
!JAC      SAVE IENT
!JAC      DATA IENT/0/
!JAC      IF (LTRACE) CALL STRACE (IENT,'SWEXCHG')
!JAC
!JAC!     --- if not parallel, return
!JAC      IF (.NOT.PARLL) RETURN
!JAC
!JAC      IPR = RESHAPE((/2,4,2,3,1,3,1,4/), (/2,4/))
!JAC      IPS = RESHAPE((/1,3,1,4,2,4,2,3/), (/2,4/))
!JAC
!JAC!TIMG      CALL SWTSTA(203)
!JAC
!JAC      ISWP = MAX(1,SWPDIR)
!JAC
!JAC      NNEIGH = IBLKAD(1)
!JAC
!JAC!     --- for all neighbouring subdomains do
!JAC
!JAC      DO INB = 1, NNEIGH
!JAC
!JAC!        --- get position
!JAC
!JAC         IPNB = IBLKAD(3*INB)
!JAC
!JAC!        --- if position corresponds to sweep selection
!JAC
!JAC         IF ( SWPDIR.EQ.0 .OR.
!JAC     &        IPNB.EQ.IPS(1,ISWP) .OR. IPNB.EQ.IPS(2,ISWP) ) THEN
!JAC
!JAC!           --- get subdomain number, pointer and size
!JAC
!JAC            IDOM   = IBLKAD(3*INB-1)
!JAC            ISTART = IBLKAD(3*INB+1)
!JAC            NOVLU  = IBLKAD(ISTART)
!JAC
!JAC!           --- store data to be sent in array WORK
!JAC
!JAC            DO K = 1, NOVLU
!JAC               WORK(K) = FIELD(KGRPNT(IBLKAD(ISTART+K)))
!JAC            END DO
!JAC
!JAC!           --- send array WORK
!JAC
!JAC            ITAG = 2
!JAC            CALL SWSENDNB ( WORK, NOVLU, SWREAL, IDOM, ITAG )
!JAC            IF (STPNOW()) RETURN
!JAC
!JAC         END IF
!JAC
!JAC      END DO
!JAC
!JAC!     --- for all neighbouring subdomains do
!JAC
!JAC      DO INB = 1, NNEIGH
!JAC
!JAC!        --- get position
!JAC
!JAC         IPNB = IBLKAD(3*INB)
!JAC
!JAC!        --- if position corresponds to sweep selection
!JAC
!JAC         IF ( SWPDIR.EQ.0 .OR.
!JAC     &        IPNB.EQ.IPR(1,ISWP) .OR. IPNB.EQ.IPR(2,ISWP) ) THEN
!JAC
!JAC!           --- get subdomain number, pointer and size
!JAC
!JAC            IDOM   = IBLKAD(3*INB-1)
!JAC            ISTART = IBLKAD(3*INB+1)
!JAC            NOVLU  = IBLKAD(ISTART)
!JAC
!JAC!           --- receive next array and store in WORK
!JAC
!JAC            ITAG  = 2
!JAC            CALL SWRECVNB ( WORK, NOVLU, SWREAL, IDOM, ITAG )
!JAC            IF (STPNOW()) RETURN
!JAC
!JAC!           --- store the received data
!JAC
!JAC            DO K = 1, NOVLU
!JAC               FIELD(KGRPNT(IBLKAD(ISTART+NOVLU+K))) = WORK(K)
!JAC            END DO
!JAC
!JAC         END IF
!JAC
!JAC      END DO
!JAC
!JAC!TIMG      CALL SWTSTO(203)
!JAC
!JAC      RETURN
!JAC      END
!****************************************************************
!
      SUBROUTINE SWEXCHG ( FIELD, KGRPNT )
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SWCOMM3                                                         40.41
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Updates geographical field array through exchanging
!     values between neighbouring subdomains
!
!  3. Method
!
!     Made use of MPI by means of SWSENDNB and SWRECVNB
!     and also block administration (stored in IBLKAD)
!
!  4. Argument variables
!
!     FIELD       geographical field array for which 'halo' values must
!                 be copied from neighbouring subdomains
!     KGRPNT      indirect addressing for grid points
!
      INTEGER KGRPNT(MXC*MYC)
      REAL    FIELD(MCGRD)
!
!  6. Local variables
!
!     IDOM  :     subdomain number
!     IENT  :     number of entries
!     INB   :     neighbour counter
!     ISTART:     pointer in array IBLKAD
!     ITAG  :     message tag for sending and receiving
!     K     :     loop counter
!     NNEIGH:     number of neighbouring subdomains
!     NOVLU :     number of overlapping unknowns
!     WORK  :     work array to store data to be sent to or
!                 received from neighbour
!
      INTEGER IDOM, IENT, INB, ISTART, ITAG, K, NNEIGH, NOVLU
      REAL    WORK(MAX(MXC,MYC))
!
!  8. Subroutines used
!
!     STPNOW           Logical indicating whether program must
!                      terminated or not
!     STRACE           Tracing routine for debugging
!     SWRECVNB         Data is received from a neighbour
!     SWSENDNB         Data is sent to a neighbour
!TIMG!     SWTSTA           Start timing for a section of code
!TIMG!     SWTSTO           Stop timing for a section of code
!
      LOGICAL STPNOW
!
!  9. Subroutines calling
!
!     ---
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     if not parallel, return
!
!     for all neighbouring subdomains do
!        get subdomain number, pointer and size
!        store data to be sent in array WORK
!        send array WORK
!
!     for all neighbouring subdomains do
!        get subdomain number, pointer and size
!        receive next array and store in WORK
!        store the received data
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWEXCHG')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

!TIMG      CALL SWTSTA(203)

      NNEIGH = IBLKAD(1)

!     --- for all neighbouring subdomains do

      DO INB = 1, NNEIGH

!        --- get subdomain number, pointer and size

         IDOM   = IBLKAD(3*INB-1)
         ISTART = IBLKAD(3*INB+1)
         NOVLU  = IBLKAD(ISTART)

!        --- store data to be sent in array WORK

         DO K = 1, NOVLU
            WORK(K) = FIELD(KGRPNT(IBLKAD(ISTART+K)))
         END DO

!        --- send array WORK

         ITAG = 2
         CALL SWSENDNB ( WORK, NOVLU, SWREAL, IDOM, ITAG )
         IF (STPNOW()) RETURN

      END DO

!     --- for all neighbouring subdomains do

      DO INB = 1, NNEIGH

!        --- get subdomain number, pointer and size

         IDOM   = IBLKAD(3*INB-1)
         ISTART = IBLKAD(3*INB+1)
         NOVLU  = IBLKAD(ISTART)

!        --- receive next array and store in WORK

         ITAG  = 2
         CALL SWRECVNB ( WORK, NOVLU, SWREAL, IDOM, ITAG )
         IF (STPNOW()) RETURN

!        --- store the received data

         DO K = 1, NOVLU
            FIELD(KGRPNT(IBLKAD(ISTART+NOVLU+K))) = WORK(K)
         END DO

      END DO

!TIMG      CALL SWTSTO(203)

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWRECVAC ( AC2, IS, J, SWPDIR, KGRPNT )
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SWCOMM3                                                         40.41
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Receives action density from neighbouring subdomains
!     depending on sweep direction
!
!  3. Method
!
!     Use of SWRECVNB and block administration (stored in IBLKAD)
!
!  4. Argument variables
!
!     AC2         action density
!     IS          start index of J-th row
!     J           J-th row
!     KGRPNT      indirect addressing for grid points
!     SWPDIR      sweep direction
!
      INTEGER IS, J, SWPDIR
      INTEGER KGRPNT(MXC,MYC)
      REAL    AC2(MDC,MSC,MCGRD)
!
!  6. Local variables
!
!     IDOM  :     subdomain number
!     IENT  :     number of entries
!     INB   :     neighbour counter
!     IPNB  :     position of neighbour (=top, bottom, right, left)
!     IPR   :     array containing positions of neighbours from
!                 which data is to be received
!     ITAG  :     message tag for sending and receiving
!     NNEIGH:     number of neighbouring subdomains
!     WORK  :     work array to store data to received from neighbour
!
      INTEGER IDOM, IENT, INB, IPNB, ITAG, NNEIGH
      INTEGER IPR(2,4)
      REAL    WORK(MDC,MSC)
!
!  8. Subroutines used
!
!     STPNOW           Logical indicating whether program must
!                      terminated or not
!     STRACE           Tracing routine for debugging
!     SWRECVNB         Data is received from a neighbour
!
      LOGICAL STPNOW
!
!  9. Subroutines calling
!
!     SWCOMP
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     if not parallel, return
!
!     for all neighbouring subdomains do
!        get position
!        if position corresponds to sweep selection
!           get subdomain number
!           receive next array and store in WORK
!           store the received data
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWRECVAC')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

      IPR = RESHAPE((/2,4,2,3,1,3,1,4/), (/2,4/))

      NNEIGH = IBLKAD(1)

!     --- for all neighbouring subdomains do

      DO INB = 1, NNEIGH

!        --- get position

         IPNB = IBLKAD(3*INB)

!        --- if position corresponds to sweep selection

         IF ( IPNB.EQ.IPR(1,SWPDIR) .OR. IPNB.EQ.IPR(2,SWPDIR) ) THEN

!           --- get subdomain number

            IDOM   = IBLKAD(3*INB-1)

!           --- receive next array and store in WORK

            ITAG  = 2
            CALL SWRECVNB ( WORK, MDC*MSC, SWREAL, IDOM, ITAG )
            IF (STPNOW()) RETURN

!           --- store the received data

            AC2(:,:,KGRPNT(IS,J)) = WORK(:,:)

         END IF

      END DO

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWSENDAC ( AC2, IE, J, SWPDIR, KGRPNT )
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SWCOMM3                                                         40.41
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Sends action density to neighbouring subdomains
!     depending on sweep direction
!
!  3. Method
!
!     Use of SWSENDNB and block administration (stored in IBLKAD)
!
!  4. Argument variables
!
!     AC2         action density
!     IE          end index of J-th row
!     J           J-th row
!     KGRPNT      indirect addressing for grid points
!     SWPDIR      sweep direction
!
      INTEGER IE, J, SWPDIR
      INTEGER KGRPNT(MXC,MYC)
      REAL    AC2(MDC,MSC,MCGRD)
!
!  6. Local variables
!
!     IDOM  :     subdomain number
!     IENT  :     number of entries
!     INB   :     neighbour counter
!     IPNB  :     position of neighbour (=top, bottom, right, left)
!     IPS   :     array containing positions of neighbours to
!                 which data is to be sent
!     ITAG  :     message tag for sending and receiving
!     NNEIGH:     number of neighbouring subdomains
!     WORK  :     work array to store data to be sent to neighbour
!
      INTEGER IDOM, IENT, INB, IPNB, ITAG, NNEIGH
      INTEGER IPS(2,4)
      REAL    WORK(MDC,MSC)
!
!  8. Subroutines used
!
!     STPNOW           Logical indicating whether program must
!                      terminated or not
!     STRACE           Tracing routine for debugging
!     SWSENDNB         Data is sent to a neighbour
!
      LOGICAL STPNOW
!
!  9. Subroutines calling
!
!     SWCOMP
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     if not parallel, return
!
!     for all neighbouring subdomains do
!        get position
!        if position corresponds to sweep selection
!           get subdomain number
!           store data to be sent in array WORK
!           send array WORK
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWSENDAC')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

      IPS = RESHAPE((/1,3,1,4,2,4,2,3/), (/2,4/))

      NNEIGH = IBLKAD(1)

!     --- for all neighbouring subdomains do

      DO INB = 1, NNEIGH

!        --- get position

         IPNB = IBLKAD(3*INB)

!        --- if position corresponds to sweep selection

         IF ( IPNB.EQ.IPS(1,SWPDIR) .OR. IPNB.EQ.IPS(2,SWPDIR) ) THEN

!           --- get subdomain number

            IDOM   = IBLKAD(3*INB-1)

!           --- store data to be sent in array WORK

            WORK(:,:) = AC2(:,:,KGRPNT(IE,J))

!           --- send array WORK

            ITAG = 2
            CALL SWSENDNB ( WORK, MDC*MSC, SWREAL, IDOM, ITAG )
            IF (STPNOW()) RETURN

         END IF

      END DO

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWCOLLECT ( FIELDGL, FIELD, FULL )
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SWCOMM3                                                         40.41
      USE M_GENARR                                                        40.31
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.51: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Mar. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.51, Feb. 05: extended to full arrays
!
!  2. Purpose
!
!     Collects geographical FIELD arrays from all nodes
!
!  3. Method
!
!     Made use of MPI by means of SWGATHER
!
!  4. Argument variables
!
!     FIELD       geographical field array in own subdomain
!     FIELDGL     global geographical field array gathered from all nodes
!     FULL        if true, full arrays are handled otherwise 1-D compact
!                 arrays are handled
!
      REAL    FIELD(*), FIELDGL(*)
      LOGICAL FULL
!
!  6. Local variables
!
!     FLDC  :     auxiliary array for collecting data
!     IARRC :     auxiliary array for collecting grid indices and counter
!     IARRL :     auxiliary array containing grid indices and counter
!     IENT  :     number of entries
!     ILEN  :     integer indicating length of an array
!     ILEN2 :     integer indicating length of another array
!     INDX  :     pointer in array
!     INDXC :     pointer in collected array
!     IOFF1 :     offset
!     IOFF2 :     another offset
!     IP    :     node number
!     IX    :     loop counter
!     IY    :     loop counter
!     KGRPTC:     auxiliary array for collecting indirect
!                 addressing for grid points
!     MXFGL :     first index w.r.t. global grid in x-direction
!     MXLGL :     last index w.r.t. global grid in x-direction
!     MYFGL :     first index w.r.t. global grid in y-direction
!     MYLGL :     last index w.r.t. global grid in y-direction
!
      INTEGER IENT, ILEN, ILEN2, INDX, INDXC, IOFF1, IOFF2, IP, IX, IY,
     &        MXFGL, MXLGL, MYFGL, MYLGL
      INTEGER IARRL(5), IARRC(5,0:NPROC-1)

      INTEGER, ALLOCATABLE :: KGRPTC(:)
      REAL,    ALLOCATABLE :: FLDC(:)
!
!  8. Subroutines used
!
!     STPNOW           Logical indicating whether program must
!                      terminated or not
!     STRACE           Tracing routine for debugging
!     SWGATHER         Gathers different amounts of data from all nodes
!
      LOGICAL STPNOW
!
!  9. Subroutines calling
!
!     ---
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     if sequential run
!        just make a copy
!     else
!        gather necessary arrays
!        copy gathered data to global array in appropriate manner
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWCOLLECT')

      IF (.NOT.PARLL) THEN

!        --- in case of sequential run, just make a copy

         IF (FULL) THEN
            DO INDX = 1, MXC*MYC
               FIELDGL(INDX) = FIELD(INDX)
            END DO
         ELSE
            DO INDX = 1, MCGRD
               FIELDGL(INDX) = FIELD(INDX)
            END DO
         END IF

      ELSE

!        --- gather necessary arrays

         IARRL(1) = MXF
         IARRL(2) = MXL
         IARRL(3) = MYF
         IARRL(4) = MYL
         IARRL(5) = MCGRD
         CALL SWGATHER (IARRC, 5*NPROC, IARRL, 5, SWINT )
         IF (STPNOW()) RETURN

         IF (.NOT.FULL) THEN
            IF (IAMMASTER) THEN
               ILEN = MXCGL*MYCGL +
     &                       4*MAX(IHALOX,IHALOY)*NPROC*MAX(MXCGL,MYCGL)
               ALLOCATE(KGRPTC(ILEN))
            END IF
            CALL SWGATHER ( KGRPTC, ILEN, KGRPNT, MXC*MYC, SWINT )
            IF (STPNOW()) RETURN
         ELSE
            IF (IAMMASTER) ALLOCATE(KGRPTC(0))
         END IF

         IF (IAMMASTER) THEN
            IF (FULL) THEN
               ILEN = MXCGL*MYCGL +
     &                       4*MAX(IHALOX,IHALOY)*NPROC*MAX(MXCGL,MYCGL)
            ELSE
               ILEN = MCGRDGL +
     &                       4*MAX(IHALOX,IHALOY)*NPROC*MAX(MXCGL,MYCGL)
            END IF
            ALLOCATE(FLDC(ILEN))
         END IF
         IF (FULL) THEN
            ILEN2 = MXC*MYC
         ELSE
            ILEN2 = MCGRD
         END IF
         CALL SWGATHER ( FLDC, ILEN, FIELD, ILEN2, SWREAL )
         IF (STPNOW()) RETURN

!        --- copy gathered data to global array in appropriate manner

         IF (IAMMASTER) THEN

            IOFF1 = 0
            IOFF2 = 0

            DO IP = 0, NPROC-1
               IF ( MXCGL.GT.MYCGL ) THEN
                  IF ( IARRC(1,IP).EQ.1 .AND. IP.EQ.0 ) THEN
                     MXFGL = 1
                  ELSE
                     MXFGL = IARRC(1,IP) + IHALOX
                  END IF
                  IF ( IARRC(3,IP).EQ.1 ) THEN
                     MYFGL = 1
                  ELSE
                     MYFGL = IARRC(3,IP) + IHALOY
                  END IF
                  IF ( IARRC(2,IP).EQ.MXCGL .AND. IP.EQ.NPROC-1 ) THEN
                     MXLGL = MXCGL
                  ELSE
                     MXLGL = IARRC(2,IP) - IHALOX
                  END IF
                  IF ( IARRC(4,IP).EQ.MYCGL ) THEN
                     MYLGL = MYCGL
                  ELSE
                     MYLGL = IARRC(4,IP) - IHALOY
                  END IF
               ELSE
                  IF ( IARRC(1,IP).EQ.1 ) THEN
                     MXFGL = 1
                  ELSE
                     MXFGL = IARRC(1,IP) + IHALOX
                  END IF
                  IF ( IARRC(3,IP).EQ.1 .AND. IP.EQ.0 ) THEN
                     MYFGL = 1
                  ELSE
                     MYFGL = IARRC(3,IP) + IHALOY
                  END IF
                  IF ( IARRC(2,IP).EQ.MXCGL ) THEN
                     MXLGL = MXCGL
                  ELSE
                     MXLGL = IARRC(2,IP) - IHALOX
                  END IF
                  IF ( IARRC(4,IP).EQ.MYCGL .AND. IP.EQ.NPROC-1 ) THEN
                     MYLGL = MYCGL
                  ELSE
                     MYLGL = IARRC(4,IP) - IHALOY
                  END IF
               END IF

               ILEN = IARRC(2,IP)-IARRC(1,IP)+1

               IF (FULL) THEN
                  DO IX = MXFGL, MXLGL
                     DO IY = MYFGL, MYLGL
                        INDX  = (IY-1)*MXCGL+IX
                        INDXC = (IY-IARRC(3,IP))*ILEN+IX
     &                             -IARRC(1,IP)+1+IOFF2
                        FIELDGL(INDX)= FLDC(INDXC)
                     END DO
                  END DO
               ELSE
                  DO IX = MXFGL, MXLGL
                     DO IY = MYFGL, MYLGL
                        INDX  = KGRPGL(IX,IY)
                        INDXC = KGRPTC((IY-IARRC(3,IP))*ILEN+IX
     &                                    -IARRC(1,IP)+1+IOFF2)
                        FIELDGL(INDX) = FLDC(INDXC+IOFF1)
                     END DO
                  END DO
               END IF

               IOFF1 = IOFF1 + IARRC(5,IP)
               IOFF2 = IOFF2 + (IARRC(2,IP)-IARRC(1,IP)+1)*
     &                         (IARRC(4,IP)-IARRC(3,IP)+1)

            END DO

         END IF

         IF (IAMMASTER) DEALLOCATE(KGRPTC,FLDC)

      END IF

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWCOLOUT ( OURQT, BLKND )                                40.51
!
!****************************************************************
!
      USE TIMECOMM                                                        40.41
      USE OCPCOMM2
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM2                                                         40.51
      USE SWCOMM3                                                         40.41
      USE OUTP_DATA
      USE M_PARALL                                                        40.31
      USE SwanGriddata, ONLY: nvertsg                                     41.36
      USE SIZES                                                           41.36
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.51: Marcel Zijlema
!     41.36: Marcel Zijlema
!
!  1. Updates
!
!     40.30, May  03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.41, Dec. 04: optimization output with respect to COMPGRID
!     40.51, Feb. 05: re-design output process in parallel mode
!     41.36, Jun. 06: collecting data for PunSWAN
!
!  2. Purpose
!
!     Collects output results
!
!  3. Method
!
!     Read individual process output files containing tables,
!     spectral outputs and block data and write them to
!     generic output files in appropriate manner
!
!  4. Argument variables
!
!     BLKND       collected array giving node number per subdomain        40.41
!     OURQT       array indicating at what time requested output          40.51
!                 is processed                                            40.51
!
      REAL    BLKND(MXCGL,MYCGL)                                          40.51
      REAL*8  OURQT(MAX_OUTP_REQ)                                         40.51
!
!  6. Local variables
!
!     CORQ  :     current item in list of request outputs
!     CROSS :     auxiliary logical array
!     CUOPS :     current item in list of point sets
!     EXIST :     logical whether a file exist or not
!     DIF   :     difference between end and actual times
!     DTTIWR:     to write time string
!     IC    :     loop variable
!     IENT  :     number of entries
!     ILPOS :     actual length of filename
!     IP    :     loop variable
!     IPROC :     loop counter
!     IRQ   :     request number
!     IT    :     time step counter
!     IT0   :     integer indicating first step of simulation
!     IT1   :     integer indicating last step of simulation
!     ITMP1 :     auxiliary integer
!     ITMP2 :     auxiliary integer
!     ITMP3 :     auxiliary integer
!     ITMP4 :     auxiliary integer
!     ITMP5 :     auxiliary integer
!     ITMP6 :     auxiliary integer
!     IUNIT :     counter for file unit numbers
!     MIP   :     total number of output points
!     MXK   :     number of points in x-direction of output frame
!     MYK   :     number of points in y-direction of output frame
!     OPENED:     logical whether a file is open or not
!     PSTYPE:     type of point set
!     RTMP1 :     auxiliary real
!     RTMP2 :     auxiliary real
!     RTMP3 :     auxiliary real
!     RTMP4 :     auxiliary real
!     RTYPE :     type of request
!     SNAMPF:     name of plot frame
!     TNEXT :     time of next requested output
!     XC    :     computational grid x-coordinate of output point
!     XP    :     user x-coordinate of output point
!     YC    :     computational grid y-coordinate of output point
!     YP    :     user y-coordinate of output point
!
      INTEGER   IC, IENT, IP, IRQ, IT, IT0, IT1, IUNIT, MIP, MXK, MYK
      INTEGER   ITMP1, ITMP2, ITMP3, ITMP4, ITMP5, ITMP6
      INTEGER   ILPOS, IPROC
      REAL*8    DIF, TNEXT
      REAL      RTMP1, RTMP2, RTMP3, RTMP4
      REAL, ALLOCATABLE :: XC(:), YC(:), XP(:), YP(:)                     40.51
      LOGICAL   OPENED
      LOGICAL   EXIST
      LOGICAL, ALLOCATABLE :: CROSS(:,:)                                  40.86
      CHARACTER PSTYPE*1, RTYPE*4, SNAMPF*8, DTTIWR*18
      TYPE(ORQDAT), POINTER :: CORQ
      TYPE(OPSDAT), POINTER :: CUOPS
!
!  8. Subroutines used
!
!     EQREAL           Logical comparing two reals
!     STPNOW           Logical indicating whether program must
!                      terminated or not
!     STRACE           Tracing routine for debugging
!     SWCOLBLK         Collects block output
!     SWCOLSPC         Collects spectral output
!     SWCOLTAB         Collects table ouput
!     SWOEXC           Computes coordinates of output points              40.51
!
      LOGICAL   EQREAL, STPNOW
!
!  9. Subroutines calling
!
!     SWMAIN
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     do for all COMPUTE commands
!        do for all time steps
!           do for all output requests
!              processing of output instructions necessary for collection
!              check time of output action
!              compute coordinates of output points
!              correct problem coordinates with offset values
!              rewrite table output by means of collection of output
!              rewrite spectral output by means of collection of output
!              rewrite block output by means of collection of output
!     close all files and delete process files
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWCOLOUT')

      IF ( NREOQ.EQ.0 ) RETURN

      NPROC = MNPROC                                                      41.36

!     --- do for all COMPUTE commands

      DO IC = 1, NCOMPT

         NSTATC = NINT(RCOMPT(IC,1))
         IF ( NSTATC.EQ.1 ) THEN
            IT0 = 0
         ELSE
            IT0 = 1
         END IF
         IT1   = NINT(RCOMPT(IC,2))
         TFINC = RCOMPT(IC,3)
         TINIC = RCOMPT(IC,4)
         DT    = RCOMPT(IC,5)
         TIMCO = TINIC

!        --- do for all time steps

         DO IT = IT0, IT1

            IF (NSTATM.GT.0) CHTIME = DTTIWR(ITMOPT, TIMCO)

!           --- do for all output requests

            CORQ => FORQ
            DO 100 IRQ = 1, NREOQ

!              --- processing of output instructions necessary
!                  for collection

!              --- check time of output action

               DIF = TFINC - TIMCO
               IF ( IT.EQ.IT0 .AND. IC.EQ.1 ) THEN
                  CORQ%OQR(1) = OURQT(IRQ)                                40.51
               END IF
               IF (CORQ%OQR(1).LT.TINIC) THEN
                  TNEXT = TINIC
               ELSE
                  TNEXT = CORQ%OQR(1)
               ENDIF
               IF ( ABS(DIF).LT.0.5*DT .AND. CORQ%OQR(2).LT.0. ) THEN
                  CORQ%OQR(1) = TIMCO
               ELSE IF ( CORQ%OQR(2).GT.0. .AND. TIMCO.GE.TNEXT ) THEN
                  CORQ%OQR(1) = TNEXT + CORQ%OQR(2)
               ELSE
                  GOTO 50
               END IF

               RTYPE  = CORQ%RQTYPE
               SNAMPF = CORQ%PSNAME

               IF (SNAMPF.EQ.'COMPGRID') THEN
                  LCOMPGRD=.TRUE.
               ELSE
                  LCOMPGRD=.FALSE.
               ENDIF

               CUOPS => FOPS
               DO
                 IF (CUOPS%PSNAME.EQ.SNAMPF) EXIT
                 IF (.NOT.ASSOCIATED(CUOPS%NEXTOPS)) GOTO 50
                 CUOPS => CUOPS%NEXTOPS
               END DO
               PSTYPE = CUOPS%PSTYPE

               IF ( PSTYPE.EQ.'F' .OR. PSTYPE.EQ.'H' ) THEN
                  MXK = CUOPS%OPI(1)
                  MYK = CUOPS%OPI(2)
                  MIP = MXK * MYK
               ELSE IF ( PSTYPE.EQ.'C' .OR. PSTYPE.EQ.'P' .OR.
     &                   PSTYPE.EQ.'N' ) THEN
                  MXK = 0
                  MYK = 0
                  MIP = CUOPS%MIP
               ELSE IF ( PSTYPE.EQ.'U' ) THEN
                  MIP = CUOPS%MIP
                  IF ( LCOMPGRD ) MIP = nvertsg
                  MXK = MIP
                  MYK = 1
               END IF

               IF (.NOT.ALLOCATED(XC)) ALLOCATE(XC(MIP))
               IF (.NOT.ALLOCATED(YC)) ALLOCATE(YC(MIP))
               IF (.NOT.ALLOCATED(XP)) ALLOCATE(XP(MIP))
               IF (.NOT.ALLOCATED(YP)) ALLOCATE(YP(MIP))
               IF (.NOT.ALLOCATED(CROSS)) ALLOCATE(CROSS(4,MIP))          40.86

!              --- compute coordinates of output points                   40.51

               ITMP1  = MXC
               ITMP2  = MYC
               ITMP3  = MCGRD
               ITMP4  = NGRBND
               ITMP5  = MXF
               ITMP6  = MYF
               RTMP1  = XCLMIN
               RTMP2  = XCLMAX
               RTMP3  = YCLMIN
               RTMP4  = YCLMAX
               MXC    = MXCGL
               MYC    = MYCGL
               MCGRD  = MCGRDGL
               NGRBND = NGRBGL
               MXF    = 1
               MYF    = 1
               XCLMIN = XCGMIN
               XCLMAX = XCGMAX
               YCLMIN = YCGMIN
               YCLMAX = YCGMAX
               CALL SWOEXC (PSTYPE              ,
     &                      CUOPS%OPI           ,CUOPS%OPR           ,
     &                      CUOPS%XP            ,CUOPS%YP            ,
     &                      MIP                 ,XP                  ,
     &                      YP                  ,XC                  ,
     &                      YC                  ,KGRPGL              ,
     &                      XGRDGL              ,YGRDGL              ,
     &                      CROSS               )                         40.86
               MXC    = ITMP1
               MYC    = ITMP2
               MCGRD  = ITMP3
               NGRBND = ITMP4
               MXF    = ITMP5
               MYF    = ITMP6
               XCLMIN = RTMP1
               XCLMAX = RTMP2
               YCLMIN = RTMP3
               YCLMAX = RTMP4

!              --- find global vertex index for output locations
!                  stored in array XC

               IF ( .NOT.LCOMPGRD ) THEN
                  XC = -1.
                  YC =  0.
                  FILENM = TRIM(LOCALDIR)//DIRCH2//'output.set'
                  ILPOS  = LEN(TRIM(LOCALDIR))
                  DO IPROC = 1, NPROC
                     WRITE(FILENM(ILPOS-3:ILPOS),'(I4.4)') IPROC-1
                     IUNIT = HIOPEN+NREOQ*(NPROC+1)+IPROC
                     INQUIRE ( FILE=FILENM, EXIST=EXIST, OPENED=OPENED )
                     IF ( .NOT.OPENED ) THEN
                        IF (EXIST) THEN
                           OPEN ( UNIT=IUNIT, FILE=FILENM,
     &                            FORM='UNFORMATTED' )
                        ELSE
                           CALL MSGERR( 4,
     &                                'file output.set does not exist' )
                           RETURN
                        ENDIF
                     ENDIF
                     READ(IUNIT) ITMP1, ITMP2
                     IF (ITMP1.EQ.IRQ) THEN
                        DO IP = 1, ITMP2
                           READ (IUNIT) ITMP3, ITMP4
                           XC(ITMP3) = REAL(ITMP4) - 1.
                        ENDDO
                     ELSE
                        CALL MSGERR( 4,
     &        'inconsistency found in SWCOLOUT: wrong request number ' )
                        RETURN
                     ENDIF
                  ENDDO
               ELSE
                  DO IP = 1, MIP
                     XC(IP) = REAL(IP) - 1.
                     YC(IP) = 0.
                  ENDDO
               ENDIF

!              --- correct problem coordinates with offset values

               DO IP = 1, MIP
                  RTMP1 = XP(IP)
                  RTMP2 = YP(IP)
                  IF (.NOT.EQREAL(RTMP1,OVEXCV(1))) XP(IP)=RTMP1+XOFFS
                  IF (.NOT.EQREAL(RTMP2,OVEXCV(2))) YP(IP)=RTMP2+YOFFS
               END DO

!              --- rewrite table output by means of collection of
!                  output locations

               IF ( RTYPE(1:3).EQ.'TAB' ) THEN
!NCF                 IF ( RTYPE.EQ.'TABC') THEN
!NCF!                 --- use "block" intermediate file facility to pass data between cores
!NCF                  CALL SWCOLBLK ( RTYPE, CORQ%OQI, CORQ%OQR, CORQ%IVTYP,
!NCF     &                            CORQ%FAC, SNAMPF, MIP, 1, IRQ,
!NCF     &                            BLKND, XC, YC )
!NCF                  IF (STPNOW()) RETURN
!NCF                 ELSE
                  CALL SWCOLTAB ( RTYPE, CORQ%OQI, CORQ%IVTYP, MIP, IRQ,
     &                            BLKND, XC, YC, XP, YP )                 40.51 40.41
                  IF (STPNOW()) RETURN
!NCF                 ENDIF
               END IF

!              --- rewrite spectral output by means of collection of
!                  output locations

               IF ( RTYPE(1:2).EQ.'SP' ) THEN
                  CALL SWCOLSPC ( RTYPE, CORQ%OQI, CORQ%OQR, MIP,         41.40 40.51 40.41
     &                            IRQ  , BLKND   , XC      , YC )         40.51
                  IF (STPNOW()) RETURN
               END IF

!              --- rewrite block output by means of collection of process
!                  output data

               IF ( RTYPE(1:3).EQ.'BLK' ) THEN
                  CALL SWCOLBLK ( RTYPE, CORQ%OQI, CORQ%OQR, CORQ%IVTYP,  41.40
     &                            CORQ%FAC, SNAMPF, MXK, MYK, IRQ,        40.51 40.41
     &                            BLKND, XC, YC )                         40.51 40.41
                  IF (STPNOW()) RETURN
               END IF

               IF (ALLOCATED(XP)) DEALLOCATE(XP)                          40.51
               IF (ALLOCATED(YP)) DEALLOCATE(YP)                          40.51
               IF (ALLOCATED(XC)) DEALLOCATE(XC)                          40.51
               IF (ALLOCATED(YC)) DEALLOCATE(YC)                          40.51

               IF (ALLOCATED(CROSS)) DEALLOCATE(CROSS)                    40.86

  50           CONTINUE
               CORQ => CORQ%NEXTORQ

 100        CONTINUE

            IF ( NSTATC.EQ.1.AND.IT.LT.IT1 ) TIMCO = TIMCO + DT

         END DO

      END DO

!     --- close all files and delete process files

      DO IUNIT = HIOPEN+1, HIOPEN+NREOQ
        INQUIRE ( UNIT=IUNIT, OPENED=OPENED )
        IF (OPENED) CLOSE(IUNIT)
      END DO
!NPUN      DO IUNIT = HIOPEN+NREOQ+1, HIOPEN+NREOQ*(NPROC+1)
      DO IUNIT = HIOPEN+NREOQ+1, HIOPEN+NREOQ*(NPROC+1)+NPROC
        INQUIRE ( UNIT=IUNIT, OPENED=OPENED )
        IF (OPENED) CLOSE ( UNIT=IUNIT, STATUS='delete' )
      END DO

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWCOLTAB ( RTYPE, OQI, IVTYP, MIP, IRQ, BLKND,           40.51
     &                      XC   , YC , XP   , YP )                       40.51
!
!****************************************************************
!
      USE OCPCOMM2, ONLY: LENFNM, DIRCH2                                  41.36
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE OUTP_DATA
      USE M_PARALL                                                        40.31
      USE SIZES, ONLY: GLOBALDIR, LOCALDIR                                41.36
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.51: Agnieszka Herman
!     40.51: Marcel Zijlema
!     41.36: Marcel Zijlema
!
!  1. Updates
!
!     40.30, May  03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.41, Dec. 04: optimization output with respect to COMPGRID
!     40.51, Feb. 05: further optimization
!     40.51, Feb. 05: re-design output process in parallel mode
!     41.36, Jun. 12: collecting data for PunSWAN
!
!  2. Purpose
!
!     Printing of table output based on point set by means of
!     collecting individual process output files
!
!  4. Argument variables
!
!     BLKND       collected array giving node number per subdomain
!     IRQ         request number
!     IVTYP       type of variable output
!     MIP         total number of output points
!     OQI         array containing output request data
!     RTYPE       type of output request
!     XC          computational grid x-coordinate of output point
!     XP          user x-coordinate of output point
!     YC          computational grid y-coordinate of output point
!     YP          user y-coordinate of output point
!
      INTEGER   IRQ, MIP, OQI(4), IVTYP(OQI(3))
      REAL      BLKND(MXCGL,MYCGL), XC(MIP), YC(MIP), XP(MIP), YP(MIP)    40.51
      CHARACTER RTYPE*4
!
!  6. Local variables
!
!     EXIST :     logical whether a file exist or not
!     FSTR  :     an auxiliary string
!     I     :     integer
!     IENT  :     number of entries
!     IF    :     first non-character in string
!     IL    :     last non-character in string
!     ILPOS :     actual length of filename
!     IP    :     loop counter
!     IPROC :     loop counter
!     IUNIT :     counter for file unit numbers
!     IUT   :     auxiliary integer representing reference number
!     IVTYPE:     type of output quantity
!     IXK   :     loop counter
!     IYK   :     loop counter
!     JVAR  :     loop counter
!     LFIELD:     actual length of a part of field OUTLIN
!     MSGSTR:     string to pass message to call MSGERR
!     NLINES:     number of lines in heading
!     NREF  :     unit reference number
!     NUMDEC:     number of decimals in the table
!     NVAR  :     number of output variables
!     OPENED:     logical whether a file is open or not
!     OUTLIN:     output line
!
      INTEGER       I, IENT, IF, IL, ILPOS, IP, IPROC, IUNIT,
     &              IUT, IVTYPE, IXK, IYK, JVAR, LFIELD, NLINES,
     &              NREF, NUMDEC, NVAR
      LOGICAL       EXIST, OPENED
      CHARACTER*80  MSGSTR
      CHARACTER*18  FSTR
      CHARACTER*512 OUTLIN
!
!  8. Subroutines used
!
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SWCOLOUT
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     if generic output file exist
!
!        if it is not open or already open thru another request
!
!           open generic output file or reset reference number
!
!           count lines of heading
!
!           open individual process output files and
!           write heading to generic output file
!
!     read output data from the proper process file and write             40.51
!     it in appropriate manner to generic output file                     40.51
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWCOLTAB')

      NREF   = OQI(1)
      NVAR   = OQI(3)
      NLINES = 0

!     --- if generic output file exist

      IF ( NREF.NE.0 .AND. NREF.NE.PRINTF ) THEN

         FILENM = OUTP_FILES(OQI(2))
         ILPOS  = INDEX(FILENM, '-0')
         IF (ILPOS.NE.0) FILENM = FILENM(1:ILPOS-1)
         ILPOS  = LEN(TRIM(LOCALDIR))+2
         FILENM = FILENM(ILPOS:LENFNM)
         FILENM = TRIM(GLOBALDIR)//DIRCH2//TRIM(FILENM)
         INQUIRE ( FILE=FILENM, OPENED=OPENED, NUMBER=IUT )

!        --- if it is not open or already open thru another request

         IF ( .NOT.OPENED .OR. NREF.LE.HIOPEN ) THEN

!           --- open generic output file or reset reference number

            IF ( .NOT.OPENED ) THEN
               NREF = HIOPEN + IRQ
               OPEN ( UNIT=NREF, FILE=FILENM )
            ELSE
               NREF = IUT
            END IF
            OQI(1) = NREF

!           --- count lines of heading

            IF ( RTYPE.NE.'TABD' ) THEN
               IF ( RTYPE.EQ.'TABP' .OR. RTYPE.EQ.'TABI' ) THEN
                  NLINES = NLINES + 7
               ELSE IF ( RTYPE.EQ.'TABT' .OR. RTYPE.EQ.'TABS' ) THEN
                  NLINES = NLINES + 5
                  IF ( RTYPE.EQ.'TABT' ) THEN
                     NLINES = NLINES + 1
                  ELSE
                     IF ( NSTATM.EQ.1 ) NLINES = NLINES + 2
                     NLINES = NLINES + 2 + MIP
                  END IF
                  DO JVAR = 1, NVAR
                     IVTYPE = IVTYP(JVAR)
                     IF ( OVSVTY(IVTYPE).LE.2 ) THEN
                        NLINES = NLINES + 3
                     ELSE
                        NLINES = NLINES + 6
                    END IF
                  END DO
               END IF
            END IF

!           --- open individual process output files and
!               write heading to generic output file

            FILENM = OUTP_FILES(OQI(2))
            ILPOS  = INDEX ( FILENM, ' ' )-1
            ILPOS  = LEN(TRIM(LOCALDIR))
            DO IPROC = 1, NPROC
               I = IPROC
               I = I - 1
               WRITE(FILENM(ILPOS-3:ILPOS),100) I
!NPUN 100           FORMAT('-',I3.3)
 100           FORMAT(I4.4)
               IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC
               INQUIRE ( FILE=FILENM, EXIST=EXIST, OPENED=OPENED )
               IF ( .NOT.OPENED ) THEN
                  IF (EXIST) THEN
                     OPEN ( UNIT=IUNIT, FILE=FILENM )
                  ELSE
!NPUN                     MSGSTR= 'file '//FILENM(1:ILPOS)//' does not exist'
                     MSGSTR= 'file '//TRIM(FILENM(ILPOS:LENFNM))//
     &                       ' does not exist'
                     CALL MSGERR( 4, MSGSTR )
                     RETURN
                  END IF
               END IF
               DO IP = 1, NLINES
                  READ (IUNIT,'(A)') OUTLIN
                  CALL TXPBLA(OUTLIN,IF,IL)
                  IF (IPROC.EQ.1) WRITE (NREF, '(A)') OUTLIN(1:IL)
               END DO
            END DO

         END IF

      END IF

!     --- read output data from the proper process file and write         40.51
!         it in appropriate manner to generic output file                 40.51

      IF ( NREF.NE.PRINTF ) THEN
         IF ( RTYPE.EQ.'TABS' .AND. NSTATM.EQ.1 ) THEN
            DO IPROC = 1, NPROC
               IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC
               READ (IUNIT,'(A)') OUTLIN
               CALL TXPBLA(OUTLIN,IF,IL)
               IF (IPROC.EQ.1) WRITE (NREF, '(A)') OUTLIN(1:IL)
            END DO
         END IF
         IPLOOP : DO IP = 1, MIP                                          40.51
            IXK = NINT(XC(IP)+100.)-99                                    41.07 40.51
            IYK = NINT(YC(IP)+100.)-99                                    41.07 40.51
            IF ( IXK.LT.1 .OR. IYK.LT.1 .OR. IXK.GT.MXCGL .OR.            41.07 40.51
     &           IYK.GT.MYCGL ) THEN                                      41.07 40.51
               CALL WREXCV                                                40.51
               CYCLE IPLOOP                                               40.51
            END IF                                                        40.51
            IPROC = 1                                                     40.51
            PROCLOOP : DO                                                 40.51
              IF ( NINT(BLKND(IXK,IYK)).EQ.IPROC ) THEN                   40.51
                 IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC         40.51
                 READ (IUNIT,'(A)') OUTLIN                                40.51
                 CALL TXPBLA(OUTLIN,IF,IL)                                40.51
                 WRITE (NREF, '(A)') OUTLIN(1:IL)                         40.51
                 EXIT PROCLOOP                                            40.51
              ELSE                                                        40.51
                 IPROC = IPROC + 1                                        40.51
                 IF ( IPROC.LE.NPROC ) THEN                               40.51
                    CYCLE PROCLOOP                                        40.51
                 ELSE                                                     40.51
                    CALL WREXCV                                           40.51
                    EXIT PROCLOOP                                         40.51
                 END IF                                                   40.51
              END IF                                                      40.51
            END DO PROCLOOP                                               40.51
         END DO IPLOOP                                                    40.51
      END IF                                                              40.51
      RETURN

      CONTAINS                                                            40.51
      SUBROUTINE WREXCV                                                   40.51
      IL = 1                                                              40.51
      OUTLIN = '    '                                                     40.51
      IF (RTYPE.EQ.'TABI') THEN                                           40.51
!        --- write point sequence number as first column                  40.51
         WRITE (OUTLIN(1:8), '(I8)') IP                                   40.51
         IL = 9                                                           40.51
         OUTLIN(IL:IL) = ' '                                              40.51
      END IF                                                              40.51
      DO JVAR = 1, NVAR                                                   40.51
         IVTYPE = IVTYP(JVAR)                                             40.51
         IF (IVTYPE.EQ.40) THEN                                           40.51
!           --- for time, 18 characters are needed                        40.51
            FSTR   = '(A18)'                                              40.51
            LFIELD = 18                                                   40.51
            OUTLIN(IL:IL+LFIELD-1) = CHTIME                               40.51
         ELSE                                                             40.51
           IF (RTYPE.EQ.'TABD') THEN                                      40.51
              FSTR   = FLT_TABLE                                          40.51
              LFIELD = FLD_TABLE                                          40.51
           ELSE                                                           40.51
              FSTR   = '(F11.X)'                                          40.51
              LFIELD = 11                                                 40.51
              NUMDEC = MAX (0,6-NINT(LOG10(ABS(OVHEXP(IVTYPE)))))         40.51
              IF (NUMDEC.GT.9) NUMDEC = 9                                 40.51
              WRITE (FSTR(6:6), '(I1)') NUMDEC                            40.51
           END IF                                                         40.51
!          --- write value into OUTLIN                                    40.51
           IF (IVTYPE.EQ.1) THEN                                          40.51
              WRITE (OUTLIN(IL:IL+LFIELD-1), FMT=FSTR) XP(IP)             40.51
           ELSE IF (IVTYPE.EQ.2) THEN                                     40.51
              WRITE (OUTLIN(IL:IL+LFIELD-1), FMT=FSTR) YP(IP)             40.51
           ELSE                                                           40.51
              WRITE (OUTLIN(IL:IL+LFIELD-1), FMT=FSTR) OVEXCV(IVTYPE)     40.51
           END IF                                                         40.51
           IF (OVSVTY(IVTYPE).EQ.3) THEN                                  40.51
              IL = IL + LFIELD + 1                                        40.51
!             --- write second component of a vectorial quantity          40.51
              WRITE (OUTLIN(IL:IL+LFIELD-1), FMT=FSTR) OVEXCV(IVTYPE)     40.51
           END IF                                                         40.51
         END IF                                                           40.51
         IL = IL + LFIELD + 1                                             40.51
         OUTLIN(IL-1:IL) = '  '                                           40.51
      END DO                                                              40.51
      CALL TXPBLA(OUTLIN,IF,IL)                                           40.51
      WRITE (NREF, '(A)') OUTLIN(1:IL)                                    40.51
      RETURN                                                              40.51
      END SUBROUTINE WREXCV                                               40.51

      END
!****************************************************************
!
      SUBROUTINE SWCOLSPC ( RTYPE, OQI, OQR, MIP, IRQ, BLKND, XC, YC )    41.40 40.51
!
!****************************************************************
!
      USE OCPCOMM2, ONLY: LENFNM, DIRCH2                                  41.36
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE OUTP_DATA
      USE M_PARALL                                                        40.31
!NCF      USE swn_outnc, ONLY: swn_outnc_colspc                               41.40
      USE SIZES, ONLY: GLOBALDIR, LOCALDIR                                41.36
      USE SwanGriddata, ONLY: xcugrdgl, ycugrdgl
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.51: Agnieszka Herman
!     40.51: Marcel Zijlema
!     41.36: Marcel Zijlema
!
!  1. Updates
!
!     40.30, May  03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.41, Dec. 04: optimization output with respect to COMPGRID
!     40.51, Feb. 05: further optimization
!     40.51, Feb. 05: re-design output process in parallel mode
!     41.36, Jun. 12: collecting data for PunSWAN
!
!  2. Purpose
!
!     Printing of spectral output based on point set by means of
!     collecting individual process output files
!
!  4. Argument variables
!
!     BLKND       collected array giving node number per subdomain
!     IRQ         request number
!     MIP         total number of output points
!     OQI         integer array containing output request data
!     OQR         real array containing output request data
!     RTYPE       type of output request
!     XC          computational grid x-coordinate of output point
!     YC          computational grid y-coordinate of output point
!
      INTEGER   IRQ, MIP, OQI(4)
      REAL*8    OQR(2)                                                    41.40
      REAL      BLKND(MXCGL,MYCGL), XC(MIP), YC(MIP)                      40.51
      CHARACTER RTYPE*4
!
!  6. Local variables
!
!     EMPTY :     logical whether a line is empty or not
!     EXIST :     logical whether a file exist or not
!     I     :     integer
!     IBLKN :     integer giving node number per subdomain
!     IENT  :     number of entries
!     IF    :     first non-character in string
!     IL    :     last non-character in string
!     ILPOS :     actual length of filename
!     IP    :     loop counter
!     IPROC :     loop counter
!     IS    :     loop counter
!     IUNIT :     counter for file unit numbers
!     IUT   :     auxiliary integer representing reference number
!     IXK   :     loop counter
!     IYK   :     loop counter
!     MSGSTR:     string to pass message to call MSGERR
!NCF!     NCF   :     true if netCDF file
!     NLINES:     number of lines in heading
!     NREF  :     unit reference number
!     OPENED:     logical whether a file is open or not
!     OTYPE :     integer indicating dimension of spectrum
!     OUTLIN:     output line
!
      INTEGER       I, IBLKN, IENT, IF, IL, ILPOS, IP, IPROC, IS,
     &              IUNIT, IUT, IXK, IYK, NLINES, NREF, OTYPE
      LOGICAL       EMPTY, EXIST, OPENED
!NCF      LOGICAL, SAVE :: NCF = .FALSE.                                      41.40
      CHARACTER*80  MSGSTR
      CHARACTER (LEN=LENSPO) OUTLIN
      CHARACTER (LEN=8) :: CRFORM = '(2F14.4)'
!
!  8. Subroutines used
!
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!NCF!     STPNOW           Logical indicating whether program must
!NCF!     swn_outnc_colspc Collect spectral output for netcdf
!     TXPBLA           Removes leading and trailing blanks in string
!
!NCF      LOGICAL STPNOW
!
!  9. Subroutines calling
!
!     SWCOLOUT
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     if generic output file exist
!
!        if it is not open or already open thru another request
!
!           open generic output file or reset reference number
!
!           count lines of heading
!
!           open individual process output files and
!           write heading to generic output file
!
!     read output data from the proper process file and write             40.51
!     it in appropriate manner to generic output file                     40.51
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWCOLSPC')

!NCF      FILENM = OUTP_FILES(OQI(2))
!NCF      NCF = INDEX( FILENM, '.NC' ).NE.0 .OR.                              41.40
!NCF     &      INDEX (FILENM, '.nc' ).NE.0                                   41.40
!NCF!
!NCF      IF ( NCF ) THEN                                                     41.40
!NCF         CALL swn_outnc_colspc( RTYPE, OQI, OQR, MIP, KGRPGL )
!NCF         RETURN
!NCF      ENDIF
!NCF
      NREF   = OQI(1)
      NLINES = 0

!     --- if generic output file exist

      IF ( NREF.NE.0 ) THEN

         FILENM = OUTP_FILES(OQI(2))
         ILPOS  = INDEX(FILENM, '-0')
         IF (ILPOS.NE.0) FILENM = FILENM(1:ILPOS-1)
         ILPOS  = LEN(TRIM(LOCALDIR))+2
         FILENM = FILENM(ILPOS:LENFNM)
         FILENM = TRIM(GLOBALDIR)//DIRCH2//TRIM(FILENM)
         INQUIRE ( FILE=FILENM, OPENED=OPENED, NUMBER=IUT )

!        --- if it is not open or already open thru another request

         IF ( .NOT.OPENED .OR. NREF.LE.HIOPEN ) THEN

!           --- open generic output file or reset reference number

            IF ( .NOT.OPENED ) THEN
               NREF = HIOPEN + IRQ
               OPEN ( UNIT=NREF, FILE=FILENM )
            ELSE
               NREF = IUT
            END IF
            OQI(1) = NREF

!           --- count lines of first part of heading

            NLINES = NLINES + 5
            IF ( LCOMPGRD    ) NLINES = NLINES - 1
            IF ( NSTATM.EQ.1 ) NLINES = NLINES + 2

!           --- open individual process output files and
!               write heading to generic output file

            FILENM = OUTP_FILES(OQI(2))
            ILPOS  = INDEX ( FILENM, ' ' )-1
            ILPOS  = LEN(TRIM(LOCALDIR))
            DO IPROC = 1, NPROC
               I = IPROC
               I = I - 1
               WRITE(FILENM(ILPOS-3:ILPOS),100) I
!NPUN 100           FORMAT('-',I3.3)
 100           FORMAT(I4.4)
               IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC
               INQUIRE ( FILE=FILENM, EXIST=EXIST, OPENED=OPENED )
               IF ( .NOT.OPENED ) THEN
                  IF (EXIST) THEN
                     OPEN ( UNIT=IUNIT, FILE=FILENM )
                  ELSE
!NPUN                     MSGSTR= 'file '//FILENM(1:ILPOS)//' does not exist'
                     MSGSTR= 'file '//TRIM(FILENM(ILPOS:LENFNM))//
     &                       ' does not exist'
                     CALL MSGERR( 4, MSGSTR )
                     RETURN
                  END IF
               END IF
               DO IP = 1, NLINES
                  READ (IUNIT,'(A)') OUTLIN
                  CALL TXPBLA(OUTLIN,IF,IL)
                  IF (IPROC.EQ.1) WRITE (NREF, '(A)') OUTLIN(1:IL)
               END DO
            END DO

!           --- write number of locations in case of whole grid
            IF ( LCOMPGRD ) THEN
               WRITE (NREF,'(I6,T41,A)') MIP, 'number of locations'
            ENDIF

!           --- write coordinates of output points to generic
!               output file

            IF ( .NOT.LCOMPGRD ) THEN
            DO IP = 1, MIP                                                40.51
               IF ( XC(IP).LT.-0.01 .OR. YC(IP).LT.-0.01 .OR.             40.51
     &              XC(IP).GT.REAL(MXCGL-1)+0.01 .OR.                     40.51
     &              YC(IP).GT.REAL(MYCGL-1)+0.01 ) THEN                   40.51
                  IBLKN = NPROC+1                                         40.51
               ELSE                                                       40.51
                  IXK   = NINT(XC(IP)+100.)-99                            41.07 40.51
                  IYK   = NINT(YC(IP)+100.)-99                            41.07 40.51
                  IBLKN = NINT(BLKND(IXK,IYK))                            40.51
               END IF                                                     40.51
               EMPTY = .TRUE.
               DO IPROC = 1, NPROC
                  IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC
                  READ (IUNIT,'(A)') OUTLIN
                  CALL TXPBLA(OUTLIN,IF,IL)
                  IF (EMPTY.AND.IBLKN.EQ.IPROC) THEN                      40.51
                     WRITE (NREF, '(A)') OUTLIN(1:IL)
                     EMPTY = .FALSE.
                  END IF
               END DO
               IF (EMPTY) WRITE (NREF, '(A)') OUTLIN(1:IL)
            END DO
            ELSE
            DO IP = 1, MIP
               WRITE (NREF, FMT=CRFORM) DBLE(xcugrdgl(IP)),
     &                                  DBLE(ycugrdgl(IP))
            ENDDO
            ENDIF

!           --- count lines of rest of heading and write heading
!               to generic output file

            NLINES = 2 + MSC
            IF (RTYPE(4:4).EQ.'C') THEN
               NLINES = NLINES + 7 + MDC
            ELSE
               NLINES = NLINES + 11
            END IF
            DO IPROC = 1, NPROC
               IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC
               DO IP = 1, NLINES
                  READ (IUNIT,'(A)') OUTLIN
                  CALL TXPBLA(OUTLIN,IF,IL)
                  IF (IPROC.EQ.1) WRITE (NREF, '(A)') OUTLIN(1:IL)
               END DO
            END DO

         END IF

      END IF

      IF (RTYPE(4:4).EQ.'C') THEN
         IF (RTYPE.EQ.'SPEC') THEN
            OTYPE = -2
         ELSE
            OTYPE =  2
         END IF
      ELSE
         IF (RTYPE.EQ.'SPE1') THEN
            OTYPE = -1
         ELSE
            OTYPE =  1
         END IF
      END IF

!     --- read output data from the proper process file and write         40.51
!         it in appropriate manner to generic output file                 40.51

      IF ( NSTATM.EQ.1 ) THEN
         DO IPROC = 1, NPROC
            IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC
            READ (IUNIT,'(A)') OUTLIN
            CALL TXPBLA(OUTLIN,IF,IL)
            IF (IPROC.EQ.1) WRITE (NREF, '(A)') OUTLIN(1:IL)
         END DO
      END IF

      IPLOOP : DO IP = 1, MIP                                             40.51
         IXK = NINT(XC(IP)+100.)-99                                       41.07 40.51
         IYK = NINT(YC(IP)+100.)-99                                       41.07 40.51
         IF ( IXK.LT.1 .OR. IYK.LT.1 .OR. IXK.GT.MXCGL .OR.               41.07 40.51
     &        IYK.GT.MYCGL ) THEN                                         41.07 40.51
            WRITE (NREF, '(A)') 'NODATA'                                  40.51
            CYCLE IPLOOP                                                  40.51
         END IF                                                           40.51
         IPROC = 1                                                        40.51
         PROCLOOP : DO                                                    40.51
           IF ( NINT(BLKND(IXK,IYK)).EQ.IPROC ) THEN                      40.51
              IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC            40.51
              READ (IUNIT,'(A)') OUTLIN                                   40.51
              CALL TXPBLA(OUTLIN,IF,IL)                                   40.51
              WRITE (NREF, '(A)') OUTLIN(1:IL)                            40.51
              IF ( OUTLIN(IF:IL).NE.'NODATA' ) THEN                       40.51
                 IF ( ABS(OTYPE).EQ.1 ) THEN                              40.51
                    DO IS = 1, MSC                                        40.51
                       READ (IUNIT,'(A)') OUTLIN                          40.51
                       CALL TXPBLA(OUTLIN,IF,IL)                          40.51
                       WRITE (NREF, '(A)') OUTLIN(1:IL)                   40.51
                    END DO                                                40.51
                 ELSE                                                     40.51
                    IF ( OUTLIN(IF:IL).NE.'ZERO' ) THEN                   40.51
                       DO IS = 1, 1+MSC                                   40.51
                          READ (IUNIT,'(A)') OUTLIN                       40.51
                          CALL TXPBLA(OUTLIN,IF,IL)                       40.51
                          WRITE (NREF, '(A)') OUTLIN(1:IL)                40.51
                       END DO                                             40.51
                    END IF                                                40.51
                 END IF                                                   40.51
              END IF                                                      40.51
              EXIT PROCLOOP                                               40.51
           ELSE                                                           40.51
              IPROC = IPROC + 1                                           40.51
              IF ( IPROC.LE.NPROC ) THEN                                  40.51
                 CYCLE PROCLOOP                                           40.51
              ELSE                                                        40.51
                 WRITE (NREF, '(A)') 'NODATA'                             40.51
                 EXIT PROCLOOP                                            40.51
              END IF                                                      40.51
           END IF                                                         40.51
         END DO PROCLOOP                                                  40.51
      END DO IPLOOP

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWCOLBLK ( RTYPE , OQI, OQR, IVTYP, FAC  ,               41.40
     &                      PSNAME, MXK, MYK, IRQ  , BLKND,               40.51
     &                      XC    , YC )                                  40.51
!
!****************************************************************
!
      USE OCPCOMM2, ONLY: LENFNM, DIRCH2                                  41.36
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM3, ONLY: NSTATM                                           41.62
      USE SWCOMM4, ONLY: KSPHER                                           41.62
      USE OUTP_DATA
      USE M_PARALL
      USE SIZES, ONLY: GLOBALDIR, LOCALDIR                                41.36
!NCF      USE SwanGridData, ONLY: XCUGRDGL, YCUGRDGL                          41.40
!NCF      USE swn_outnc                                                       41.40
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.51: Agnieszka Herman
!     40.51: Marcel Zijlema
!     41.36: Marcel Zijlema
!     41.62: Andre van der Westhuysen
!
!  1. Updates
!
!     40.31, Dec. 03: New subroutine
!     40.41, Jun. 04: some improvements with respect to MATLAB
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.41, Dec. 04: optimization output with respect to COMPGRID
!     40.51, Feb. 05: further optimization
!     40.51, Feb. 05: re-design output process in parallel mode
!     41.36, Jun. 12: collecting data for PunSWAN
!     41.62, Nov. 15: included output fields for wave partitioning
!
!  2. Purpose
!
!     Writing of block output by means of collecting
!     individual process output files
!
!  4. Argument variables
!
!     BLKND       collected array giving node number per subdomain
!     FAC         factors of multiplication of block output
!     IRQ         request number
!     IVTYP       type of variable output
!     MXK         number of points in x-direction of output frame
!     MYK         number of points in y-direction of output frame
!     OQI         integer array containing output request data
!     OQR         real array containing output request data
!     PSNAME      name of output locations
!     RTYPE       type of output request
!     XC          computational grid x-coordinate of output point
!     YC          computational grid y-coordinate of output point
!
      INTEGER   MXK, MYK, IRQ, OQI(4), IVTYP(OQI(3))
      REAL      BLKND(MXCGL,MYCGL), XC(MXK*MYK), YC(MXK*MYK)              40.51
      REAL*8    OQR(2)                                                    41.40
      REAL      FAC(OQI(3))                                               41.40
      CHARACTER RTYPE*4, PSNAME*8
!
!  6. Local variables
!
!     CTIM  :     string representing date of computation
!     DFAC  :     multiplication factor of block output
!     EXIST :     logical whether a file exist or not
!     FMAX  :     auxiliary real
!     FTIP  :     auxiliary real
!     FTIP1 :     auxiliary real
!     FTIP2 :     auxiliary real
!     I     :     integer
!     IBLKN :     integer giving node number per subdomain
!     IDLA  :     lay-out indicator
!     IF    :     first non-character in string
!     IFAC  :     auxiliary integer
!     IENT  :     number of entries
!     IL    :     last non-character in string
!     ILPOS :     actual length of filename
!     IP    :     loop counter
!     IPD   :     switch for printing on paper or writing to file
!     IPROC :     loop counter
!     IREC  :     direct access file record counter
!     IUNIT :     counter for file unit numbers
!     IUT   :     auxiliary integer representing reference number
!     IVTYPE:     type of output quantity
!     IXK   :     loop counter
!     IYK   :     loop counter
!     JVAR  :     loop counter
!     MATLAB:     indicates whether binary Matlab files are used
!     MSGSTR:     string to pass message to call MSGERR
!     NAMVAR:     name of MATLAB variable
!     NREF  :     unit reference number
!     NVAR  :     number of output variables
!     OPENED:     logical whether a file is open or not
!     VOQ   :     collected output variables
!
      INTEGER      I, IBLKN, IDLA, IENT, IF, IFAC, IL, ILPOS, IP, IPD,
     &             IPROC, IUNIT, IUT, IXK, IYK, IVTYPE, J, JVAR, NREF,
     &             NVAR
      REAL         DFAC, FMAX, FTIP, FTIP1, FTIP2
      LOGICAL      EXIST, OPENED
      INTEGER, SAVE :: IREC(MAX_OUTP_REQ)=0                               40.51 40.41
      LOGICAL, SAVE :: MATLAB=.FALSE.                                     40.41
!NCF      LOGICAL, SAVE :: NCF   =.FALSE.                                     41.40
      LOGICAL, SAVE :: RAWPRT=.FALSE.                                     41.62
      CHARACTER*80 MSGSTR
      CHARACTER (LEN=20) :: CTIM                                          40.41
      CHARACTER (LEN=30) :: NAMVAR                                        40.41
      CHARACTER (LEN=80) :: HTXT(3)                                       41.62
      REAL, ALLOCATABLE :: VOQ(:,:)
      INTEGER VOQR(NMOVAR)                                                41.62
!
!  8. Subroutines used
!
!     MSGERR           Writes error message
!     SBLKPT           Writes block output to an ASCII file
!     STRACE           Tracing routine for debugging
!     SWRMAT           Writes block output to a binary Matlab file
!     TABHED           Prints heading
!     TXPBLA           Removes leading and trailing blanks in string      40.41
!
!  9. Subroutines calling
!
!     SWCOLOUT
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     if generic output file exist
!
!        if it is not open or already open thru another request
!
!           open generic output file or reset reference number
!
!           open individual process output files
!
!     read data from the proper process file and write                    40.51
!     it in appropriate manner to generic output file                     40.51
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWCOLBLK')

      NREF = OQI(1)
      NVAR = OQI(3)
      IDLA = OQI(4)

      IF ( RTYPE.EQ.'BLKP' ) THEN
        IPD = 1
        IF (NREF.EQ.PRINTF) CALL TABHED ('SWAN', PRINTF)
      ELSE IF ( RTYPE.EQ.'BLKD' ) THEN
        IPD = 2
      ELSE
        IPD = 3
      ENDIF

!     --- if generic output file exist

      IF ( NREF.NE.0 .AND. NREF.NE.PRINTF ) THEN

         FILENM = OUTP_FILES(OQI(2))
         ILPOS  = INDEX(FILENM, '-0')
         IF (ILPOS.NE.0) FILENM = FILENM(1:ILPOS-1)
         ILPOS  = LEN(TRIM(LOCALDIR))+2
         FILENM = FILENM(ILPOS:LENFNM)
         FILENM = TRIM(GLOBALDIR)//DIRCH2//TRIM(FILENM)
         INQUIRE ( FILE=FILENM, OPENED=OPENED, NUMBER=IUT )
         MATLAB = INDEX( FILENM, '.MAT' ).NE.0 .OR.                       40.41
     &            INDEX (FILENM, '.mat' ).NE.0                            40.41
!NCF         NCF    = INDEX( FILENM, '.NC'  ).NE.0 .OR.                       41.40
!NCF     &            INDEX (FILENM, '.nc'  ).NE.0                            41.40
         RAWPRT = INDEX( FILENM, '.RAW' ).NE.0 .OR.                       41.62
     &            INDEX (FILENM, '.raw' ).NE.0                            41.62

!        --- if it is not open or already open thru another request

         IF ( .NOT.OPENED .OR. NREF.LE.HIOPEN ) THEN

!           --- open generic output file or reset reference number

!NCF            IF (NCF .AND. .NOT.OPENED) THEN                               41.40
!NPUN!NCF               CALL swn_outnc_openblockfile(FILENM, MYK, MXK,
!NPUN!NCF     &                                      OVLNAM, XGRDGL, YGRDGL,
!NPUN!NCF     &                                      OQI, OQR, IVTYP, IRQ)
!NCF               CALL swn_outnc_openblockfile(FILENM, MYK, MXK,
!NCF     &                                      OVLNAM, XCUGRDGL, YCUGRDGL,
!NCF     &                                      OQI, OQR, IVTYP, IRQ)
!NCF                IUT    = HIOPEN + IRQ
!NCF                OPENED = .TRUE.
!NCF            END IF
!NCF
            IF ( .NOT.OPENED ) THEN
               NREF = HIOPEN + IRQ
               OPEN ( UNIT=NREF, FILE=FILENM )
            ELSE
               NREF = IUT
            END IF
            OQI(1) = NREF

            IF (MATLAB .AND. .NOT.OPENED) THEN
               CLOSE(NREF)
               OPEN(UNIT=NREF, FILE=FILENM, FORM='UNFORMATTED',
     &              STATUS='REPLACE',
!MatL4     &              ACCESS='DIRECT', RECL=1)
     &              ACCESS='DIRECT', RECL=4)                              41.08
               IREC(IRQ) = 1                                              40.51
            END IF

            IF (RAWPRT .AND. IPD.EQ.1 .AND. .NOT.OPENED) THEN             41.62
              IF (NSTATM.EQ.1) THEN
                 WRITE (HTXT(1),'(a)') ' yyyymmdd hhmmss'
              ELSE
                 WRITE (HTXT(1),'(a)') ''
              ENDIF
              IF (KSPHER.EQ.0) THEN
                 WRITE (HTXT(2),'(a)') '         x             y'
              ELSE
                 WRITE (HTXT(2),'(a)') '       lat           lon'
              ENDIF
              WRITE (HTXT(3),'(a)')
     &              '       name       nprt depth uabs  udir cabs  cdir'
!
              WRITE(NREF,'(A26)') 'SWAN PARTITIONED DATA FILE'
              WRITE(NREF,'(A16,A24,A50)') TRIM(HTXT(1)), TRIM(HTXT(2)),
     &                                    TRIM(HTXT(3))
              WRITE(NREF,'(A31,A20)') '        hs     tp     lp       ',
     &                                'theta     sp      wf'
            END IF                                                        41.62

!           --- open individual process output files

            FILENM = OUTP_FILES(OQI(2))
            ILPOS  = INDEX ( FILENM, ' ' )-1
            ILPOS  = LEN(TRIM(LOCALDIR))
            DO IPROC = 1, NPROC
               I = IPROC
               I = I - 1
               WRITE(FILENM(ILPOS-3:ILPOS),100) I
!NPUN 100           FORMAT('-',I3.3)
 100           FORMAT(I4.4)
               IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC
               INQUIRE ( FILE=FILENM, EXIST=EXIST, OPENED=OPENED )
               IF ( .NOT.OPENED ) THEN
                  IF (EXIST) THEN
                     OPEN ( UNIT=IUNIT, FILE=FILENM, FORM='UNFORMATTED')
                  ELSE
!NPUN                     MSGSTR= 'file '//FILENM(1:ILPOS)//' does not exist'
                     MSGSTR= 'file '//TRIM(FILENM(ILPOS:LENFNM))//
     &                       ' does not exist'
                     CALL MSGERR( 4, MSGSTR )
                     RETURN
                  END IF
               END IF
            END DO

         END IF

      END IF

!     --- read data from the proper process file and write                40.51
!         it in appropriate manner to generic output file                 40.51

      CTIM = CHTIME                                                       40.41
      CALL TXPBLA(CTIM,IF,IL)                                             40.41
      CTIM(9:9)='_'                                                       40.41

      IF (.NOT.RAWPRT) THEN
         ALLOCATE(VOQ(MXK*MYK,2))
      ELSE
         ALLOCATE(VOQ(MXK*MYK,NMOVAR))
         VOQR=(/ (J, J=1, NMOVAR) /)
      ENDIF

      JLOOP: DO JVAR = 1, NVAR

         IVTYPE = IVTYP(JVAR)
         DFAC   = FAC(JVAR)

         IF (.NOT.RAWPRT) THEN                                            41.62
            VOQ = OVEXCV(IVTYPE)
            J   = 1
         ELSE
            J   = IVTYPE
         ENDIF

         DO IP = 1, MXK*MYK                                               40.51
            IXK = NINT(XC(IP)+100.)-99                                    41.07 40.51
            IYK = NINT(YC(IP)+100.)-99                                    41.07 40.51
            IF ( IXK.LT.1 .OR. IYK.LT.1 .OR. IXK.GT.MXCGL .OR.
     &           IYK.GT.MYCGL ) THEN
               IBLKN = 0
               IPROC = NPROC+1
            ELSE
               IBLKN = NINT(BLKND(IXK,IYK))
               IPROC = 1
            END IF
            PROCLOOP1 : DO                                                40.51
              IF ( IBLKN.EQ.IPROC ) THEN                                  40.51
                 IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC         40.51
                 READ (IUNIT) VOQ(IP,J)                                   40.51
                 EXIT PROCLOOP1                                           40.51
              ELSE                                                        40.51
                 IPROC = IPROC + 1                                        40.51
                 IF ( IPROC.LE.NPROC ) THEN                               40.51
                    CYCLE PROCLOOP1                                       40.51
                 ELSE                                                     40.51
                    EXIT PROCLOOP1                                        40.51
                 END IF                                                   40.51
              END IF                                                      40.51
            END DO PROCLOOP1                                              40.51
         END DO                                                           40.51

         IF ( OVSVTY(IVTYPE).GE.3 ) THEN

            DO IP = 1, MXK*MYK                                            40.51
               IXK = NINT(XC(IP)+100.)-99                                 41.07 40.51
               IYK = NINT(YC(IP)+100.)-99                                 41.07 40.51
               IF ( IXK.LT.1 .OR. IYK.LT.1 .OR. IXK.GT.MXCGL .OR.
     &              IYK.GT.MYCGL ) THEN
                  IBLKN = 0
                  IPROC = NPROC+1
               ELSE
                  IBLKN = NINT(BLKND(IXK,IYK))
                  IPROC = 1
               END IF
               PROCLOOP2 : DO                                             40.51
                 IF ( IBLKN.EQ.IPROC ) THEN                               40.51
                    IUNIT=HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC        40.51
                    READ (IUNIT) VOQ(IP,J+1)                              40.51
                    EXIT PROCLOOP2                                        40.51
                 ELSE                                                     40.51
                    IPROC = IPROC + 1                                     40.51
                    IF ( IPROC.LE.NPROC ) THEN                            40.51
                       CYCLE PROCLOOP2                                    40.51
                    ELSE                                                  40.51
                       EXIT PROCLOOP2                                     40.51
                    END IF                                                40.51
                 END IF                                                   40.51
               END DO PROCLOOP2                                           40.51
            END DO                                                        40.51
         END IF

         IF (RAWPRT) CYCLE JLOOP                                          41.62

         IF ( IPD.EQ.1 ) THEN
            IF ( DFAC.LE.0. ) THEN
               IF ( OVHEXP(IVTYPE).LT.0.5E10 ) THEN
                  IFAC = INT (10.+LOG10(OVHEXP(IVTYPE))) - 13
               ELSE
                  IF ( OVSVTY(IVTYPE).EQ.1 ) THEN
                     FMAX = 1.E-8
                     DO IP = 1, MXK*MYK
                        FTIP = ABS(VOQ(IP,1))
                        FMAX = MAX (FMAX, FTIP)
                     END DO
                  ELSE IF ( OVSVTY(IVTYPE).EQ.2 ) THEN
                     FMAX = 1000.
                  ELSE IF ( OVSVTY(IVTYPE).EQ.3 ) THEN
                     FMAX = 1.E-8
                     DO IP = 1, MXK*MYK
                        FTIP1 = ABS(VOQ(IP,1))
                        FTIP2 = ABS(VOQ(IP,2))
                        FMAX  = MAX (FMAX, FTIP1, FTIP2)
                     END DO
                  END IF
                  IFAC = INT (10.+LOG10(FMAX)) - 13
               END IF
               DFAC = 10.**IFAC
            END IF
         ELSE
           IF ( DFAC.LE.0. ) DFAC = 1.
         END IF

         IF (OVSVTY(IVTYPE) .LT. 3) THEN
            IF (MATLAB) THEN
               IF (IL.EQ.1 .OR. IVTYPE.LT.3 .OR. IVTYPE.EQ.52) THEN       40.94 40.41
                  NAMVAR = OVSNAM(IVTYPE)                                 40.41
               ELSE                                                       40.41
                  NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//   40.41
     &                     '_'//CTIM                                      40.41
               END IF                                                     40.41
               CALL SWRMAT( MYK, MXK, NAMVAR, VOQ(1,1), NREF,
     &                      IREC(IRQ), IDLA, OVEXCV(IVTYPE) )             40.51
!NCF            ELSE IF (NCF) THEN                                            41.40
!NCF               IF ( IVTYPE.GT.2.AND.IVTYPE.NE.40 ) THEN
!NCF                  CALL swn_outnc_appendblock(MYK, MXK, IVTYPE, OQI(1),
!NCF     &                                       IRQ, VOQ(1,1),
!NCF     &                                       OVEXCV(IVTYPE), 1)
!NCF               END IF
            ELSE
               CALL SBLKPT( IPD, NREF, DFAC, PSNAME, OVUNIT(IVTYPE),
     &                      MXK, MYK, IDLA, OVLNAM(IVTYPE), VOQ(1,1) )
            END IF
         ELSE
            IF (MATLAB) THEN
               IF (IL.EQ.1) THEN                                          40.41
                  NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//   40.41
     &                     '_x'                                           40.41
               ELSE                                                       40.41
                  NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//   40.41
     &                     '_x_'//CTIM                                    40.41
               END IF                                                     40.41
               CALL SWRMAT( MYK, MXK, NAMVAR,
     &                 VOQ(1,1), NREF, IREC(IRQ), IDLA, OVEXCV(IVTYPE) )  40.51
               IF (IL.EQ.1) THEN                                          40.41
                  NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//   40.41
     &                     '_y'                                           40.41
               ELSE                                                       40.41
                  NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//   40.41
     &                     '_y_'//CTIM                                    40.41
               END IF                                                     40.41
               CALL SWRMAT( MYK, MXK, NAMVAR,
     &                 VOQ(1,2), NREF, IREC(IRQ), IDLA, OVEXCV(IVTYPE) )  40.51
!NCF            ELSE IF (NCF) THEN                                            41.40
!NCF               IF ( IVTYPE.GT.3 ) THEN
!NCF                  CALL swn_outnc_appendblock(MYK, MXK, IVTYPE, OQI(1),
!NCF     &                                       IRQ, VOQ(1,1),
!NCF     &                                       OVEXCV(IVTYPE), 1)
!NCF                  CALL swn_outnc_appendblock(MYK, MXK, IVTYPE, OQI(1),
!NCF     &                                       IRQ, VOQ(1,2),
!NCF     &                                       OVEXCV(IVTYPE), 2)
!NCF               END IF
            ELSE
               CALL SBLKPT( IPD, NREF, DFAC, PSNAME, OVUNIT(IVTYPE),
     &                      MXK, MYK, IDLA, OVLNAM(IVTYPE)//'X-comp',
     &                      VOQ(1,1) )
               CALL SBLKPT( IPD, NREF, DFAC, PSNAME, OVUNIT(IVTYPE),
     &                      MXK, MYK, IDLA, OVLNAM(IVTYPE)//'Y-comp',
     &                      VOQ(1,2) )
            END IF
         END IF

      END DO JLOOP

!     generate a dump of the raw partition data, if appropriate
      IF (RAWPRT) CALL SRAWPT ( NREF, VOQR, VOQ, MXK, MYK )               41.62

!NCF      IF ( NCF ) CALL swn_outnc_close_on_end(OQI(1), IRQ)                 41.40
!NCF
      IF (IPD.EQ.1 .AND. NREF.EQ.PRINTF) WRITE (PRINTF, 111)

      DEALLOCATE(VOQ)

  111 FORMAT (///)

      RETURN
      END
!JAC!****************************************************************
!JAC!
!JAC      SUBROUTINE SWBLKCOL ( MCOLR, KGRPNT )
!JAC!
!JAC!****************************************************************
!JAC!
!JAC      USE OCPCOMM4                                                        40.41
!JAC      USE SWCOMM3                                                         40.41
!JAC      USE M_PARALL                                                        40.31
!JAC!
!JAC      IMPLICIT NONE
!JAC!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!JAC!
!JAC!  0. Authors
!JAC!
!JAC!     40.30: Marcel Zijlema
!JAC!     40.41: Marcel Zijlema
!JAC!
!JAC!  1. Updates
!JAC!
!JAC!     40.30, Mar. 03: New subroutine
!JAC!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!JAC!
!JAC!  2. Purpose
!JAC!
!JAC!     Colours the subdomains with red, yellow, green and black
!JAC!     in order to determine the sequence of sweeps during the
!JAC!     iteration process
!JAC!
!JAC!  3. Method
!JAC!
!JAC!     The four-colour ordering scheme is based on colouring
!JAC!     subdomains in each direction in alternating way with
!JAC!     red and black (similar to a chessboard colouring),
!JAC!     whereafter the product of resulting color directions
!JAC!     is taken.
!JAC!     Initially, the most left and under blocks are coloured
!JAC!     red, whereas other are uncoloured. Then, the colors
!JAC!     are propagated from left to right and from bottom to
!JAC!     top. Finally, when both directions of all subdomains
!JAC!     have been coloured, the subdomains are coloured by
!JAC!     taking the product of two color directions.
!JAC!
!JAC!  4. Argument variables
!JAC!
!JAC!     KGRPNT      indirect addressing for grid points
!JAC!     MCOLR       flag to indicate multi-colouring
!JAC!                 of subdomains (.TRUE.) or not (.FALSE.)
!JAC!
!JAC      INTEGER KGRPNT(MXC,MYC)
!JAC      LOGICAL MCOLR
!JAC!
!JAC!  5. Parameter variables
!JAC!
!JAC!     ITERMAX:    maximum number of iterations
!JAC!     IWHITE:     integer used to colour subdomains 'white'
!JAC!
!JAC      INTEGER IWHITE, ITERMAX
!JAC      PARAMETER (IWHITE=0, ITERMAX=100)
!JAC!
!JAC!  6. Local variables
!JAC!
!JAC!     CHARS :     character for passing info to MSGERR
!JAC!     ICOLNB:     color of neighbouring subdomain
!JAC!     ICONV :     indicator for convergence (0=yes, 1=no)
!JAC!     IENT  :     number of entries
!JAC!     IF    :     first non-character in string
!JAC!     IL    :     last non-character in string
!JAC!     ITER  :     iteration count
!JAC!     IXCOL :     color in x-direction of own subdomain
!JAC!     IYCOL :     color in y-direction of own subdomain
!JAC!     MSGSTR:     string to pass message to call MSGERR
!JAC!     XCOL  :     field array containing present color in x-direction
!JAC!     YCOL  :     field array containing present color in y-direction
!JAC!
!JAC      INTEGER      ICOLNB, ICONV, IENT, IF, IL, ITER, IXCOL, IYCOL
!JAC      CHARACTER*20 INTSTR, CHARS
!JAC      CHARACTER*80 MSGSTR
!JAC      REAL, ALLOCATABLE :: XCOL(:), YCOL(:)
!JAC!
!JAC!  8. Subroutines used
!JAC!
!JAC!     INTSTR           Converts integer to string
!JAC!     MSGERR           Writes error message
!JAC!     STRACE           Tracing routine for debugging
!JAC!     SWEXCHG          Updates geographical field array through
!JAC!                      exchanging values between subdomains
!JAC!     SWREDUCE         Performs a global reduction
!JAC!     TXPBLA           Removes leading and trailing blanks in string
!JAC!
!JAC!  9. Subroutines calling
!JAC!
!JAC!     SWMAIN
!JAC!
!JAC! 10. Error messages
!JAC!
!JAC!     ---
!JAC!
!JAC! 11. Remarks
!JAC!
!JAC!     ---
!JAC!
!JAC! 12. Structure
!JAC!
!JAC!     if not parallel or no colouring, return
!JAC!
!JAC!     initially, both x- and y-direction of the most left
!JAC!     and under subdomains are coloured red and that of all
!JAC!     other subdomains are marked white, i.e. not being
!JAC!     coloured
!JAC!
!JAC!     while not all subdomains are coloured do
!JAC!        exchange colors of both directions between subdomains
!JAC!        adjust color in x-direction of own subdomain
!JAC!        adjust color in y-direction of own subdomain
!JAC!        check whether all subdomains have been coloured
!JAC!
!JAC!     if not all subdomains have been coloured gives message and stops
!JAC!
!JAC!     finally, subdomains are coloured by taking the product of two
!JAC!     color directions
!JAC!
!JAC! 13. Source text
!JAC!
!JAC      SAVE IENT
!JAC      DATA IENT/0/
!JAC      IF (LTRACE) CALL STRACE (IENT,'SWBLKCOL')
!JAC
!JAC      IBCOL = IRED
!JAC
!JAC!     --- if not parallel or no colouring, return
!JAC      IF (.NOT.PARLL .OR. .NOT.MCOLR) RETURN
!JAC
!JAC      ALLOCATE(XCOL(MCGRD))
!JAC      ALLOCATE(YCOL(MCGRD))
!JAC
!JAC!     --- initially, both x- and y-direction of the most left
!JAC!         and under subdomains are coloured red and that of all
!JAC!         other subdomains are marked white, i.e. not being
!JAC!         coloured
!JAC
!JAC      IF ( MXF.EQ.1 ) THEN
!JAC         IXCOL = IRED
!JAC      ELSE
!JAC         IXCOL = IWHITE
!JAC      END IF
!JAC      IF ( MYF.EQ.1 ) THEN
!JAC         IYCOL = IRED
!JAC      ELSE
!JAC         IYCOL = IWHITE
!JAC      END IF
!JAC
!JAC!     --- while not all subdomains are coloured, adjust the
!JAC!         color of each direction of own subdomain depending
!JAC!         on the color directions of left and under neighbours
!JAC
!JAC      DO ITER = 1, ITERMAX
!JAC
!JAC         ICONV = 0
!JAC
!JAC!        --- exchange colors of both directions between subdomains
!JAC
!JAC         XCOL = REAL(IXCOL)
!JAC         YCOL = REAL(IYCOL)
!JAC         CALL SWEXCHG ( XCOL, 0, KGRPNT )
!JAC         CALL SWEXCHG ( YCOL, 0, KGRPNT )
!JAC
!JAC!        --- adjust color in x-direction of own subdomain
!JAC
!JAC         IF ( IXCOL.EQ.IWHITE ) THEN
!JAC
!JAC            ICOLNB = NINT(XCOL(KGRPNT(1,2)))
!JAC            IF ( ICOLNB.EQ.IRED ) THEN
!JAC               IXCOL = IBLACK
!JAC            ELSE IF ( ICOLNB.EQ.IBLACK ) THEN
!JAC               IXCOL = IRED
!JAC            ELSE
!JAC               ICONV = 1
!JAC            END IF
!JAC
!JAC         END IF
!JAC
!JAC!        --- adjust color in y-direction of own subdomain
!JAC
!JAC         IF ( IYCOL.EQ.IWHITE ) THEN
!JAC
!JAC            ICOLNB = NINT(YCOL(KGRPNT(2,1)))
!JAC            IF ( ICOLNB.EQ.IRED ) THEN
!JAC               IYCOL = IBLACK
!JAC            ELSE IF ( ICOLNB.EQ.IBLACK ) THEN
!JAC               IYCOL = IRED
!JAC            ELSE
!JAC               ICONV = 1
!JAC            END IF
!JAC
!JAC         END IF
!JAC
!JAC!        --- check whether all subdomains have been coloured
!JAC
!JAC         CALL SWREDUCE( ICONV, 1, SWINT, SWMAX )
!JAC         IF ( ICONV.EQ.0 ) GOTO 100
!JAC
!JAC      END DO
!JAC 100  CONTINUE
!JAC
!JAC!     --- if not all subdomains have been coloured
!JAC!         gives message and stops
!JAC
!JAC      IF (ICONV.NE.0 .AND. (IXCOL.EQ.IWHITE .OR. IYCOL.EQ.IWHITE)) THEN
!JAC         CHARS = INTSTR(INODE)
!JAC         CALL TXPBLA(CHARS,IF,IL)
!JAC         MSGSTR = 'Subdomain '//CHARS(IF:IL)//
!JAC     &            'has not been coloured'
!JAC         CALL MSGERR ( 4, MSGSTR )
!JAC         RETURN
!JAC      END IF
!JAC
!JAC!     --- finally, subdomains are coloured by taking the product
!JAC!         of two color directions
!JAC
!JAC      IF ( IXCOL.EQ.IRED .AND. IYCOL.EQ.IRED ) THEN
!JAC         IBCOL = IRED
!JAC      ELSE IF ( IXCOL.EQ.IBLACK .AND. IYCOL.EQ.IRED ) THEN
!JAC         IBCOL = IYELOW
!JAC      ELSE IF ( IXCOL.EQ.IBLACK .AND. IYCOL.EQ.IBLACK ) THEN
!JAC         IBCOL = IGREEN
!JAC      ELSE IF ( IXCOL.EQ.IRED .AND. IYCOL.EQ.IBLACK ) THEN
!JAC         IBCOL = IBLACK
!JAC      END IF
!JAC
!JAC      DEALLOCATE(XCOL,YCOL)
!JAC
!JAC      RETURN
!JAC      END
