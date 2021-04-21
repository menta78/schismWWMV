subroutine SwanSumOverNodes ( rval )
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
!   Authors
!
!   40.95: Marcel Zijlema
!
!   Updates
!
!   40.95, July 2008: New subroutine
!
!   Purpose
!
!   Performs a global sum of reals over all nodes
!
!   Modules used
!
    use ocpcomm4
    use SIZES, only: MNPROC
    use GLOBAL, only: COMM
    use MESSENGER, only: IERR
    use mpi
!
    implicit none
!
!   Argument variables
!
    real, intent(inout)     :: rval     ! input value
!
!   Local variables
!
    integer                 :: count    ! length of array to be collect
    integer, save           :: ient = 0 ! number of entries in this subroutine
    real                    :: sumval   ! sum total of all input values from all subdomains
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanSumOverNodes')
    !
    ! if not parallel, return
    !
    if ( MNPROC==1 ) return
    !
!TIMG    call SWTSTA(202)
    sumval = 0.
    count  = 1
    call MPI_ALLREDUCE ( rval, sumval, count, MPI_REAL4, MPI_SUM, COMM, IERR )
    !
    rval = sumval
!TIMG    call SWTSTO(202)
    !
end subroutine SwanSumOverNodes