      SUBROUTINE sed_write_debug()
!--------------------------------------------------------------------!
! This subroutine writes debug or additional information to a file   !
! e.g. mirror.out                                                    !
!                                                                    !
! Author: Ligia Pinto                                                !
! Date: xx/xx/xxxx                                                   !
!                                                                    !
! History: 2012/12 - F.Ganthy : form homogenisation of sediments     !
!          routines                                                  !
!                                                                    !
!--------------------------------------------------------------------!

      USE elfe_glbl, ONLY : nea,ntracers
      USE elfe_msgp, ONLY : myrank

      IMPLICIT NONE

!- Local variables --------------------------------------------------!

      INTEGER :: debug_file

!- Start Statement --------------------------------------------------!

!     16 means mirror.out
      debug_file = 16

      IF(myrank==0) THEN
        WRITE(debug_file,*)'SED: sed_write_debug'
        !WRITE(debug_file,*)'it: ',it
        !WRITE(debug_file,*)'dt: ',dt
        !WRITE(debug_file,*)'nstp: ',nstp
        !WRITE(debug_file,*)'nnew: ',nnew

!        DO i=1,nea
!!            WRITE(debug,'(A,I5,A,E,A,E,A,E)')'nea:',i,' bustr:',bustr(i),      &
!            WRITE(debug,'(A,I5,A,E11.8,A,E11.8,A,E11.8)')'nea:',i,' bustr:',bustr(i),      &
!     &            ' bvstr:',bvstr(i),' tau_w:',tau_w(i)
!          ENDDO !nea

!          DO i=1,Nbed
!            DO j=1,nea
!              DO k=1,ntracers
!!                WRITE(debug,'(A,I1,A,I5,A,I5,A,F11.8,A,F,A,F)')'Nbed:',i,      &
!                WRITE(debug,'(A,I1,A,I5,A,I5,A,F11.8,A,F11.8,A,F11.8)')'Nbed:',i,      &
!     &                ' nea:',j,' ntracers:',k,' bed_frac:',bed_frac(i,j,k),   &
!     &                ' bed_mass(1):',bed_mass(i,j,1,k),                       &
!     &                ' bed_mass(2):',bed_mass(i,j,2,k)
!              ENDDOo !ntracers
!            ENDDO !nea
!          ENDDO !Nbed

!          DO i=1,Nbed
!            DO j=1,nea
!              WRITE(debug,'(A,I1,A,I5,A,F11.8,A,F11.8,A,F11.8,A,F11.8)')       &
!     &              'Nbed:',i,' nea:',j,                                       &
!     &              ' bed_thick:',bed_thick(j),' bed(ithck):',bed(i,j,ithck),  &
!     &              ' bed(iaged):',bed(i,j,iaged),' bed(iporo):',bed(i,j,iporo)
!            ENDDO !nea
!          ENDDO !Nbed

!          DO i=1,nea
!!            WRITE(debug,'(A,I5,A,F11.8,A,F11.8,A,F11.8,A,F11.8,A,F)')          &
!            WRITE(debug,'(A,I5,A,F11.8,A,F11.8,A,F11.8,A,F11.8,A,F11.8)')          &
!     &            'nea:',i,' bottom(itauc):',bottom(i,itauc),                  &
!     &            ' bottom(isd50):',bottom(i,isd50),' bottom(iwsed):',bottom(i,iwsed),&
!     &            ' bottom(idens):',bottom(i,idens),' bottom(iactv):',bottom(i,iactv)
!          ENDDO !nea

!!           DO i=1,nea
!!             WRITE(debug,'(A,I5,A,E,A,E,A,E,A,E)')'nea:',i,                    &
!!     &              ' bottom(nea,itauc):',bottom(i,itauc),                     &
!!     &              ' bottom(nea,isd50):',bottom(i,isd50),                     &
!!     &              ' bottom(nea,iwsed):',bottom(i,iwsed),                     &
!!     &              ' bottom(nea,idens):',bottom(i,idens)
!!           ENDDO !nea

      ENDIF !myrank

!--------------------------------------------------------------------!
      END SUBROUTINE sed_write_debug
