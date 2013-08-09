#include "wwm_functions.h"
#ifdef SELFE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE STOKES_STRESS_INTEGRAL_SELFE
        USE DATAPOOL
        implicit none
        integer IP, k, ID, IS, IL
        real(rkind) eF1, eF2, eDelta, TheInt, eDep, eHeight
        real(rkind) eFrac, eFracB, eQuot
        real(rkind) eQuot1, eScal
        real(rkind) eOmega, eMult, kD, eSinc
        real(rkind) USTOKESpart, VSTOKESpart, eJPress
        real(rkind) ACLOC, eWk, eSigma, eLoc, eSinhkd, eSinh2kd, eSinhkd2
        real(rkind) PPTAIL, CETAIL, CKTAIL
        logical DoTail
        real(rkind) eWkReal
        real(rkind) SumHeight
        real(rkind) eJPress_loc, eProd, eUint, eVint
        real(rkind) eUSTOKES_loc(NVRT), eVSTOKES_loc(NVRT)
        real(rkind) ZZETA

        DO IP=1,MNP
          eDep=SHYFZETA(NLEV(IP),MNP)
          eUSTOKES_loc=0
          eVSTOKES_loc=0
          eJpress_loc=0
          DO IS=1,MSC
            eMult=SPSIG(IS)*DDIR*DS_INCR(IS)
            eWk=WK(IP,IS)
            kD=MIN(KDMAX, eWk*eDep)
            eWkReal=kD/eDep
            eSinh2kd=MySINH(2*kD)
            eSinhkd=MySINH(kD)
            eSinhkd2=eSinhkd**2
            eSigma=SPSIG(IS)
            eUint=0
            eVint=0
            DO ID=1,MDC
              eLoc=AC2(IP,IS,ID)*eMult
              eJPress=G9*(kD/eSinh2kd)*(1/eDep) * eLoc
              eJPress_loc=eJPress_loc + eJPress
              eUint=eUint + eLoc*COSTH(ID)
              eVint=eVint + eLoc*SINTH(ID)
            END DO
            DO IL = KBP(IP), NVRT
              ZZETA = ZETA(IL,IP)-ZETA(KBP(IP),IP) !from bottom; 'z+D'
              eFrac=ZZETA/eDep
! Need some better understanding of vertical levels in SELFE
! for putting those correction terms.
!              eHeight=z_w_loc(k)-z_w_loc(k-1)
!              eFracB=eHeight/eDep
!              eSinc=SINH(kD*eFracB)/(kD*eFracB)
!              eQuot1=eSinc*MyCOSH(2*kD*eFrac)/eSinhkd2
              eQuot1=MyCOSH(2*kD*eFrac)/eSinhkd2
              eProd=eSigma*eWkReal*eQuot1
              eUSTOKES_loc(k)=eUSTOKES_loc(k) + eUint*eProd
              eVSTOKES_loc(k)=eVSTOKES_loc(k) + eVint*eProd
            ENDDO
          END DO
          STOKES_X(:,IP)=eUSTOKES_loc
          STOKES_Y(:,IP)=eVSTOKES_loc
          JPRESS(IP)=eJPress_loc
        ENDDO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
       SUBROUTINE RADIATION_STRESS_SELFE

        use elfe_glbl, only: iplg,errmsg
        USE elfe_msgp !, only : myrank,parallel_abort

        USE DATAPOOL
        IMPLICIT NONE

        INTEGER :: IP,IL,IS,ID
        REAL(rkind)  :: ACLOC(MSC,MDC)
        REAL(rkind)  :: COSE2, SINE2, COSI2
        REAL(rkind)  :: EWK(MNP),EWS(MNP),EWN(MNP),ETOT(MNP),MDIR(MNP)
        REAL(rkind)  :: m0, m0d, tmp, EHFR, ELOC, EFTAIL, ZZETA, DVEC2RAD
        REAL(rkind)  :: DS, D, KW, KD, SINH2KD, SINHKW, COSH2KW, COSHKW, COSHKD, ETOTS, ETOTC, EWSIG, S11, S22
        REAL(rkind)  :: WNTMP,WKTMP,WCGTMP,WCTMP,WN,WKDEPTMP
        REAL(rkind)  :: WSTMP, DEPLOC

        INTEGER      :: ND1,ND2
        REAL(rkind)  :: SINHKD,FSS(NVRT,MNP),FCS(NVRT,MNP),FSC(NVRT,MNP),FCC(NVRT,MNP)
        REAL(rkind)  :: dr_dxy(2,NVRT,nsa),HTOT,SXX3D0(NVRT,MNP),SYY3D0(NVRT,MNP),SXY3D0(NVRT,MNP)
        REAL(rkind)  :: WILD1(NVRT,MNP),WILD2(NVRT,MNP),WILD3(2,NVRT,nsa),WILD4(3,NVRT,MNP),DSPX,DSPY, WILD5(10)

!GD: imet_dry allows to choose between 2 different methods to compute the
!derivative at the sides between wet and dry elements:
!
! imet_dry=1 : only the values at the 2 nodes of the side are used to
! compute the derivative (this older method showed to provide inconsistent
! wave force at the wet/dry interface).
!
! imet_dry=2 : a 4-point stencil (the 3 wet nodes and an artificial
! node at the center of the side) is used to compute the derivative.
! This method is similar to using shape functions to compute the
! derivative at the center of the element and assigning this value to the
! the side center.
     
        IMET_DRY = 2 

        SXX3D(:,:) = ZERO
        SYY3D(:,:) = ZERO
        SXY3D(:,:) = ZERO

        EFTAIL = ONE / (PTAIL(1) - ONE)

        ETOT = ZERO
        MDIR = ZERO

        IF (LETOT) THEN
!AR: Estimate zeroth moment m0, mean wave direction, dominant wave number, dominant sigma ...
          DO IP = 1, MNP
!            IF (ABS(IOBP(IP)) .GT. 0) CYCLE
            IF (DEP(IP) .LT. DMIN) CYCLE
            DEPLOC = DEP(IP)
            ACLOC = AC2(IP,:,:)
            m0    = ZERO
            EWSIG  = ZERO
            ETOTS  = ZERO
            ETOTC  = ZERO
            IF (MSC .GE. 2) THEN
              DO ID = 1, MDC
                m0d = ZERO
                DO IS = 2, MSC
                  tmp = 0.5_rkind*(SPSIG(IS)*ACLOC(IS,ID)+SPSIG(IS-1)*ACLOC(IS-1,ID))*DS_INCR(IS)*DDIR
                  m0 = m0 + tmp
                  EWSIG  = EWSIG  + SPSIG(IS) * tmp
                  m0d = m0d + tmp
                END DO
                IF (MSC > 3) THEN
                  EHFR = ACLOC(MSC,ID) * SPSIG(MSC)
                  m0 = m0 + DDIR * EHFR * SPSIG(MSC) * EFTAIL
                endif
                ETOTC  = ETOTC + m0d * COS(SPDIR(ID))
                ETOTS  = ETOTS + m0d * SIN(SPDIR(ID))
              END DO
            ELSE
              DS = SGHIGH - SGLOW
              DO ID = 1, MDC
                m0d = ACLOC(1,ID) * DS * DDIR
                m0 = m0 + m0d
              END DO
            END IF
            ETOT(IP) = m0
            IF (m0 .GT. small .and. .not. dep(ip) .lt. dmin) then
              EWS(IP) = EWSIG/m0
              WSTMP = EWSIG/m0
              CALL ALL_FROM_TABLE(WSTMP,DEPLOC,WKTMP,WCGTMP,WKDEPTMP,WNTMP,WCTMP)
              EWN(IP) = WNTMP
              EWK(IP) = WKTMP
              MDIR(IP) = DVEC2RAD (ETOTC, ETOTS)
            ELSE
              EWS(IP)  = ZERO 
              EWN(IP)  = ZERO 
              EWK(IP)  = 10. 
              MDIR(IP) = ZERO 
            END IF 
          END DO !IP
        END IF !LETOT

!AR: Here comes the whole story ... 
! Etot = 1/16 * Hs² = 1/8 * Hmono² => Hs² = 2 * Hmono² => Hs = sqrt(2) * Hmono => Hmono = Hs / SQRT(2) 
! Etot = 1/16 * Hs² = 1/16 * (4 * sqrt(m0))² = m0 
! Etot = 1/8 * Hmono² ... so the problem for the analytical solution evolved because we treat the Etot from Hs and Hmono there is a factor of 2 between this!
! Or in other words for the analytical solution we impose a Hs = X[m], we integrate m0 out of it and get Etot, since this Etot is a function of Hs and not Hmono^X^O
! it needs the factor of 2 between it! This should make now things clear forever. So the question is not how we calculate the total energy the question is 
! what is defined on the boundary that means we should always recalculate the boundary in terms of Hs =  SQRT(2) * Hmono !!!
! Or saying it again in other words our boundary conditions is wrong if we impose Hmono in wwminput.nml !!!

        SXX3D = ZERO
        SXY3D = ZERO
        SYY3D = ZERO
        WWAVE_FORCE=ZERO
        IF (RADFLAG .EQ. 'LON') THEN
          RSXX = ZERO
          RSXY = ZERO
          RSYY = ZERO
          DO IP = 1, MNP
!            IF (ABS(IOBP(IP)) .GT. 0) CYCLE
            IF (DEP(IP) .LT. DMIN) CYCLE
            IF (.NOT. LETOT) THEN
              ACLOC = AC2(IP,:,:)
              DO ID = 1, MDC
                DO IS = 2, MSC
                  ELOC  = 0.5_rkind * (SPSIG(IS)*ACLOC(IS,ID)+SPSIG(IS-1)*ACLOC(IS-1,ID))*DS_INCR(IS)*DDIR
                  COSE2 = COS(SPDIR(ID))**TWO
                  SINE2 = SIN(SPDIR(ID))**TWO
                  COSI2 = COS(SPDIR(ID)) * SIN(SPDIR(ID))
                  WN    = CG(IP,IS) / ( SPSIG(IS)/WK(IP,IS) )
                  RSXX(IP) = RSXX(IP) + ( WN * COSE2 + WN - 0.5_rkind) * ELOC   ! Units = [ 1/s + 1/s - 1/s ] * m²s = m²
                  RSXY(IP) = RSXY(IP) + ( WN * COSI2          ) * ELOC
                  RSYY(IP) = RSYY(IP) + ( WN * SINE2 + WN - 0.5_rkind) * ELOC
                ENDDO
              ENDDO
            ELSE IF (LETOT) THEN
              RSXX(IP) =  ETOT(IP) * (EWN(IP)*((EWK(IP)*SIN(MDIR(IP)))**TWO/EWK(IP)**TWO+ONE)-0.5_rkind)
              RSXY(IP) =  ETOT(IP) *  EWN(IP)* EWK(IP)*SIN(MDIR(IP))*EWK(IP)*COS(MDIR(IP))* ONE/EWK(IP)
              RSYY(IP) =  ETOT(IP) * (EWN(IP)*((EWK(IP)*COS(MDIR(IP)))**TWO/EWK(IP)**TWO+ONE)-0.5_rkind)
            END IF 
          END DO

          DO IP = 1, MNP
!            IF (ABS(IOBP(IP)) .GT. 0) CYCLE
            IF (DEP(IP) .LT. DMIN)  THEN
              SXX3D(:,IP) = ZERO
              SXY3D(:,IP) = ZERO
              SYY3D(:,IP) = ZERO
            ELSE
              SXX3D(:,IP) = RSXX(IP) / DEP(IP) * G9
              SXY3D(:,IP) = RSXY(IP) / DEP(IP) * G9
              SYY3D(:,IP) = RSYY(IP) / DEP(IP) * G9
            END IF
          END DO
          !Store as double for force later
          SXX3D0 = SXX3D
          SXY3D0 = SXY3D
          SYY3D0 = SYY3D
        ELSE IF (RADFLAG .EQ. 'XIA') THEN
          IF (LETOT) THEN
            DO IP = 1, MNP
!              IF (ABS(IOBP(IP)) .GT. 0) CYCLE
              IF (DEP(IP) .LT. DMIN .OR. EWK(IP) .LT. THR) CYCLE
              DO IL = KBP(IP), NVRT
                ZZETA = ZETA(IL,IP)-ZETA(KBP(IP),IP) !from bottom; 'z+D'
                KW           = EWK(IP) * ZZETA !k*(z+D)
                KD           = EWK(IP) * DEP(IP) !k*D
                SINH2KD      = MySINH(MIN(KDMAX,TWO*KD))
                COSHKD       = MyCOSH(MIN(KDMAX,KD))
                SINHKW       = MySINH(MIN(KDMAX,KW))
                COSH2KW      = MyCOSH(MIN(KDMAX,TWO*KW))
                COSHKW       = MyCOSH(MIN(KDMAX,KW))
                IF(ABS(SINH2KD) .LT. THR) THEN
                  write(errmsg,*)'R.S.: div by 0 (0);',iplg(IP),EWK(IP),DEP(IP)
                  call parallel_abort(errmsg)
                endif
                tmp         = -ETOT(IP) * EWK(IP) / SINH2KD * (COSH2KW - ONE)                      - &
     &                         ETOT(IP) * (ZETA(IL,IP)-ZETA(NVRT,IP)) / DEP(IP)**2                  + & 
     &                         ETOT(IP) * KW * SINHKW / ( DEP(IP) * COSHKD ) - &
     &                         ETOT(IP) / DEP(IP)  *  (ONE - COSHKW / COSHKD)
                SXX3D(IL,IP) = (ETOT(IP) * EWK(IP) / SINH2KD * (COSH2KW + ONE) * COS(MDIR(IP))**2+tmp)* G9
                SYY3D(IL,IP) = (ETOT(IP) * EWK(IP) / SINH2KD * (COSH2KW + ONE) * SIN(MDIR(IP))**2+tmp)* G9
                SXY3D(IL,IP) = ETOT(IP) * EWK(IP) / SINH2KD * (COSH2KW + ONE) * COS(MDIR(IP))*SIN(MDIR(IP))* G9
              END DO !IL
            END DO !IP
          ELSE !not mono ... random waves version ... treat eveery wave packet like a single monochr. wave with the height of Hs 
            DO IP = 1, MNP
!              IF (ABS(IOBP(IP)) .GT. 0) CYCLE
              IF (DEP(IP) .LT. DMIN) CYCLE
              ACLOC = AC2(IP,:,:)
              DO IL = KBP(IP), NVRT
                ZZETA = ZETA(IL,IP)-ZETA(KBP(IP),IP) !from bottom
                DO IS = 1, MSC !freq
                  KW           = WK(IP,IS) * ZZETA !k*(z+D)
                  KD           = WK(IP,IS) * DEP(IP) !k*D
                  SINH2KD      = MySINH(MIN(KDMAX,TWO*KD))
                  COSHKD       = MyCOSH(MIN(KDMAX,KD))
                  SINHKW       = MySINH(MIN(KDMAX,KW))
                  COSH2KW      = MyCOSH(MIN(KDMAX,TWO*KW))
                  COSHKW       = MyCOSH(MIN(KDMAX,KW))
                  IF(ABS(SINH2KD) .LT. THR) call parallel_abort('R.S.: div 0 (0)')
                  DO ID = 1, MDC !direction
                    !Dimension of ELOC = m^2
                    ELOC = ACLOC(IS,ID) * SIGPOW(IS,2) * DDIR * FRINTF
                    IF (ELOC .LT. 10E-8) CYCLE
                      tmp          =-ELOC * WK(IP,IS) / SINH2KD * (COSH2KW - ONE)            - &
     &                               ELOC * (ZETA(IL,IP)-ZETA(NVRT, IP)) / DEP(IP)**TWO        + &
     &                               ELOC * KW * SINHKW / ( DEP(IP) * COSHKD ) - &
     &                               ELOC / DEP(IP) *  (ONE - COSHKW / COSHKD)
                      tmp=tmp*G9

                      SXX3D(IL,IP) = SXX3D(IL,IP) + &
     &                               G9*ELOC * WK(IP,IS) / SINH2KD * (COSH2KW + ONE) * COS(SPDIR(ID))**TWO+tmp
                      SYY3D(IL,IP) = SYY3D(IL,IP) + &
     &                               G9*ELOC * WK(IP,IS) / SINH2KD * (COSH2KW + ONE) * SIN(SPDIR(ID))**TWO+tmp
                      SXY3D(IL,IP) = SXY3D(IL,IP)+G9*ELOC * WK(IP,IS) / SINH2KD * (COSH2KW + ONE) * SIN(SPDIR(ID))*COS(SPDIR(ID))
                  END DO
                END DO
              END DO !IL
            END DO !IP
          END IF !mono

          !Store as double for force later
          SXX3D0=SXX3D
          SXY3D0=SXY3D
          SYY3D0=SYY3D

!ZYL: Mellor 2003; random waves only; only good for traditional sigma coord.
!     did not include the last term in Warner et al. (2008) as it is not in Mellor 2005
        ELSE IF (RADFLAG .EQ. 'MEL') THEN
          IF(LETOT) call parallel_abort('R.S.: no mono for Mellor 88')
          IF(KZ/=1.or.THETA_F>1.e-4) call parallel_abort('R.S.: MEL must use sigma')
          !WILD1: \sum_{dir}{E} at nodes; WILD2: K*D at nodes
          WILD1=0; WILD2=0; WILD4=0 !for exceptions
          FSS=0; FCS=0; FSC=0; FCC=0
          DO IS = 1, MSC !freq
            DO IP = 1, MNP
!              IF (ABS(IOBP(IP)) .GT. 0) CYCLE
              IF (DEP(IP) .LT. DMIN) CYCLE
              KD           = WK(IP,IS) * DEP(IP) !k*D
              WILD2(:,IP)  = KD
              SINHKD       = MySINH(MIN(KDMAX,KD))
              COSHKD       = MyCOSH(MIN(KDMAX,KD))
              IF(ABS(SINHKD) .LT. THR) call parallel_abort('R.S.: div by 0 (1)')
              DO IL = KBP(IP), NVRT
                ZZETA = ZETA(IL,IP)-ZETA(KBP(IP),IP) !from bottom
                KW           = WK(IP,IS) * ZZETA !k*(z+D)
                SINHKW       = MySINH(MIN(KDMAX,KW))
                COSHKW       = MyCOSH(MIN(KDMAX,KW))
                FSS(IL,IP)   = SINHKW/SINHKD
                FCS(IL,IP)   = COSHKW/SINHKD
                FSC(IL,IP)   = SINHKW/COSHKD
                FCC(IL,IP)   = COSHKW/COSHKD
                DO ID = 1, MDC !direction
                  !Dimension of ELOC = m^2
                  ELOC = AC2(IP,IS,ID) * SIGPOW(IS,2) * DDIR * FRINTF * 1.
                  IF (ELOC .LT. 10E-8) CYCLE
                  WILD1(IL,IP)=WILD1(IL,IP)+ELOC
                  SXX3D(IL,IP) = SXX3D(IL,IP)+G9*ELOC*WK(IP,IS)*(FCS(IL,IP)*FCC(IL,IP)*(COS(SPDIR(ID))**2+1.)-FSS(IL,IP)*FCS(IL,IP)) 
                  SYY3D(IL,IP) = SYY3D(IL,IP)+G9*ELOC*WK(IP,IS)*(FCS(IL,IP)*FCC(IL,IP)*(SIN(SPDIR(ID))**2+1.)-FSS(IL,IP)*FCS(IL,IP))
                  SXY3D(IL,IP) = SXY3D(IL,IP)+G9*ELOC*WK(IP,IS)*FCS(IL,IP)*FCC(IL,IP)*SIN(SPDIR(ID))*COS(SPDIR(ID))
                END DO !ID
              END DO !IL
            END DO !IP

            !Compute horizontal derivatives of WILD[1,2], and temporarily store them in dr_dxy and WILD3
            !dr_dxy: d{sum{E}}/d{x,y} - not a function of sigma
            !WILD3: d{K*D}/d{x,y} - not a function of sigma
            !Can we sum up E or KD over freq. before differentiation to save exchange time?
            call hgrad_nodes(IMET_DRY,0,NVRT,MNP,nsa,WILD1,dr_dxy)
            call hgrad_nodes(IMET_DRY,0,NVRT,MNP,nsa,WILD2,WILD3)
            call exchange_s3d_2(dr_dxy)
            call exchange_s3d_2(WILD3)

            !Compute vertical derivatives of FCC etc. and store them in WILD4
            DO IP = 1, MNP
!              IF (ABS(IOBP(IP)) .GT. 0) CYCLE
              IF (DEP(IP) .LT. DMIN) CYCLE
              KD= WK(IP,IS) * DEP(IP) !k*D
              DO IL = KBP(IP), NVRT
                WILD4(1,IL,IP)=FSC(IL,IP)*KD !d{FCC}/d{sigma}
                WILD4(2,IL,IP)=FCS(IL,IP)*KD !d{FSS}/d{sigma}
                WILD4(3,IL,IP)=FSS(IL,IP)*KD !d{FCS}/d{sigma}
              END DO !IL
            ENDDO !IP

            !Vertical derivative parts of wave forces (more later for horizontal parts): D{SPX}/D{sigma}
            do IP=1,nsa !sides
              if(idry_s(IP)==1) CYCLE
              ND1=isidenode(IP,1) !1st node of the side
              ND2=isidenode(IP,2)
              HTOT=(eta2(ND1)+eta2(ND2)+DEP8(ND1)+DEP8(ND2))/2
              if(HTOT<=0) call parallel_abort('RADIATION_STRESS_SELFE: HTOT')
              DO IL=1,NVRT
                WILD5(1)=(WILD4(1,IL,ND1)+WILD4(1,IL,ND2))/2 !d{FCC}/d{sigma} at side
                WILD5(2)=(WILD4(2,IL,ND1)+WILD4(2,IL,ND2))/2 !d{FSS}/d{sigma} at side
                WILD5(3)=(WILD4(3,IL,ND1)+WILD4(3,IL,ND2))/2 !d{FCS}/d{sigma} at side
                WILD5(4)=(FCC(IL,ND1)+FCC(IL,ND2))/2 !FCC at side
                WILD5(5)=(FSS(IL,ND1)+FSS(IL,ND2))/2 !FSS at side
                WILD5(6)=(FCS(IL,ND1)+FCS(IL,ND2))/2 !FCS at side
                WILD5(7)=(WILD1(IL,ND1)+WILD1(IL,ND2))/2 !sum{E} at side
                WILD5(8)=(WILD2(IL,ND1)+WILD2(IL,ND2))/2 !K*D at side
                DSPX=dr_dxy(1,IL,IP)/2*((WILD5(1)-WILD5(2))*WILD5(5)+(WILD5(4)-WILD5(5))*WILD5(2))+&
     &WILD5(7)*WILD3(1,IL,IP)*((WILD5(1)-WILD5(2))*WILD5(6)*(1+SIGMACOR(IL))+ &
     &(WILD5(4)-WILD5(5))*WILD5(3)*(1+SIGMACOR(IL))+(WILD5(4)-WILD5(5))*WILD5(6))
                DSPY=dr_dxy(2,IL,IP)/2*((WILD5(1)-WILD5(2))*WILD5(5)+(WILD5(4)-WILD5(5))*WILD5(2))+&
     &WILD5(7)*WILD3(2,IL,IP)*((WILD5(1)-WILD5(2))*WILD5(6)*(1+SIGMACOR(IL))+ &
     &(WILD5(4)-WILD5(5))*WILD5(3)*(1+SIGMACOR(IL))+(WILD5(4)-WILD5(5))*WILD5(6))

                WWAVE_FORCE(IL,IP,1)=WWAVE_FORCE(IL,IP,1)+G9*DSPX/HTOT !m/s/s
                WWAVE_FORCE(IL,IP,2)=WWAVE_FORCE(IL,IP,2)+G9*DSPY/HTOT
              END DO !IL
            enddo !IP; sides
          END DO !IS

          !Store as double for force later
          SXX3D0=SXX3D
          SXY3D0=SXY3D
          SYY3D0=SYY3D

        ELSE
            call parallel_abort('R.S.: unknown R.S. model') 
        END IF !RADFLAG 

!       Integrate over depth for checking
        RSXX = ZERO
        DO IP = 1, MNP
          IF(IDRY(IP)==1) CYCLE
          DO IL = KBP(IP)+1, NVRT 
            RSXX(IP) = RSXX(IP) + 0.5_rkind*( SXX3D(IL,IP)+SXX3D(IL-1,IP) ) * ABS((ZETA(IL,IP) - ZETA(IL-1,IP)))/G9
          END DO !IL
        END DO !IP
!'
!       Computation in double precision here
!        IF (RADFLAG.EQ.'LON'.OR.RADFLAG.EQ.'XIA') THEN
!         SXX3D0() etc. should have dimension of m^2/s/s, defined at nodes and whole levels.
!         Use same arrays to temporarily store properly scaled Sxx etc
!         write(12,*)'Checking Sxx,Sxy,Syy:'
          do IP=1,MNP
            IF(IDRY(IP)==1) then
              SXX3D0(:,IP)=ZERO
              SYY3D0(:,IP)=ZERO
              SXY3D0(:,IP)=ZERO
              cycle
            endif

            do IL=KBP(IP),NVRT
!             D*(Sxx, Sxy, Syy)/rho in Xia et al. (2004) 
!             After this the dimension of sdbt should be m^3/s/s
              SXX3D0(IL,IP)=(DEP8(IP)+ETA2(IP))*SXX3D0(IL,IP) !D*Sxx/rho
              SXY3D0(IL,IP)=(DEP8(IP)+ETA2(IP))*SXY3D0(IL,IP) !D*Sxy/rho
              SYY3D0(IL,IP)=(DEP8(IP)+ETA2(IP))*SYY3D0(IL,IP) !D*Syy/rho
            enddo !k
          enddo !IP

!         Compute radiation stress force 
!         wwave_force(:,1:nsa,1:2) = Rsx, Rsy in my notes (the terms in momen. eq.)
!         and has a dimension of m/s/s
          call hgrad_nodes(IMET_DRY,0,NVRT,MNP,nsa,SXX3D0,dr_dxy)
          call exchange_s3d_2(dr_dxy)
          do IS=1,nsa
            if(idry_s(IS)==0) then
              HTOT=(eta2(isidenode(IS,1))+eta2(isidenode(IS,2))+DEP8(isidenode(IS,1))+DEP8(isidenode(IS,2)))/2
              if(HTOT<=0) call parallel_abort('RADIATION_STRESS: (999)')
              do il = 1, nvrt
                WWAVE_FORCE(il,IS,1)=WWAVE_FORCE(il,IS,1)-dr_dxy(1,il,IS)/HTOT
              end do
            endif
          enddo !IS

          call hgrad_nodes(IMET_DRY,0,NVRT,MNP,nsa,SYY3D0,dr_dxy)
          call exchange_s3d_2(dr_dxy)
          do IS=1,nsa
            if(idry_s(IS)==0) then
              HTOT=(eta2(isidenode(IS,1))+eta2(isidenode(IS,2))+DEP8(isidenode(IS,1))+DEP8(isidenode(IS,2)))/2
              if(HTOT<=0) call parallel_abort('RADIATION_STRESS: (998)')
              do il = 1, nvrt 
                WWAVE_FORCE(il,IS,2)=WWAVE_FORCE(il,IS,2)-dr_dxy(2,il,IS)/HTOT
              end do
            endif
          enddo !IS

          call hgrad_nodes(IMET_DRY,0,NVRT,MNP,nsa,SXY3D0,dr_dxy)
          call exchange_s3d_2(dr_dxy)
!
!          write(12,*)'Checking R.S.'
!
          do IS=1,nsa
            if(idry_s(IS)==0) then
              HTOT=(eta2(isidenode(IS,1))+eta2(isidenode(IS,2))+DEP8(isidenode(IS,1))+DEP8(isidenode(IS,2)))/2
              if(HTOT<=0) call parallel_abort('RADIATION_STRESS: (997)')
              WWAVE_FORCE(:,IS,1)=WWAVE_FORCE(:,IS,1)-dr_dxy(2,:,IS)/HTOT
              WWAVE_FORCE(:,IS,2)=WWAVE_FORCE(:,IS,2)-dr_dxy(1,:,IS)/HTOT
            endif
          enddo !IS

        END SUBROUTINE
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAVEFORCE

      USE DATAPOOL
      IMPLICIT NONE
      INTEGER  ID, IS, IP

      REAL(kind=rkind)  ::   COSE2, SINE2, COSI2
      REAL(kind=rkind)  ::   ELOC
      REAL(kind=rkind)  ::   CK
      REAL(kind=rkind)  ::   DSIGMA
      REAL(kind=rkind)  ::   D_RSXX(MNP,2)
      REAL(kind=rkind)  ::   D_RSXY(MNP,2)
      REAL(kind=rkind)  ::   D_RSYY(MNP,2)

      RSXX(:) = zero 
      RSXY(:) = zero 
      RSYY(:) = zero

      DO IP = 1, MNP

        DO IS = 1, MSC

        CK = CG(IP,IS) * WK(IP,IS) ! CG ~ Group Velocity ; K ~ Wave Number

!       WE HAVE TO TAKE CARE ABOUT THE GROUP VELOCITY IF CURRENTS ARE PRESENT !!!
!       THE GROUP VELOCITY WILL BE DOPPLER SHIFTET IN PRESENCE OF CURRENTS
!       IT WILL BECOME A FUNCTION OF THE WAVE DIRECTION

          DO ID = 1, MDC

          ELOC = AC2(IP,IS,ID)  * DDIR * FRINTF

!         AC2 ~ Local Action Density
!         Directional Property of the Action Spectrum

          COSE2 = COS(SPDIR(ID))**2
          SINE2 = SIN(SPDIR(ID))**2
          COSI2 = COS(SPDIR(ID))*SIN(SPDIR(ID))
!                    /
!     Sxx = rho grav | ((N cos^2(theta) + N - 1/2) sig Ac) d sig d theta
!                   /
!                     /
!     Sxy = rho grav | (N sin(theta) cos(theta) sig Ac) d sig d theta
!                   /
!                     /
!     Syy = rho grav | ((N sin^2(theta) + N - 1/2) sig Ac) d sig d theta
!
          RSXX(IP) = RSXX(IP) + (CK * COSE2 + CK - SPSIG(IS)/TWO) * ELOC
          RSXY(IP) = RSXY(IP) +  CK * COSI2 * ELOC
          RSYY(IP) = RSYY(IP) + (CK * SINE2 + CK - SPSIG(IS)/TWO) * ELOC
!
          ENDDO
        ENDDO
      END DO

      RSXX(:) = RSXX(:) * RHOW * G9
      RSXY(:) = RSXY(:) * RHOW * G9
      RSYY(:) = RSYY(:) * RHOW * G9
!
!     Arrays are called by Adress and not by Value that means "Array(1,1)" ....
!
      CALL DIFFERENTIATE_XYDIR(RSXX(1),D_RSXX(1,1),D_RSXX(1,2))
      CALL DIFFERENTIATE_XYDIR(RSXY(1),D_RSXY(1,1),D_RSXY(1,2))
      CALL DIFFERENTIATE_XYDIR(RSYY(1),D_RSYY(1,1),D_RSYY(1,2))

      IF (LSPHE) THEN
         DO IP = 1, MNP
            D_RSXX(IP,1) = D_RSXX(IP,1)/( DEGRAD*REARTH*COS(YP(IP)*DEGRAD) )
            D_RSXY(IP,1) = D_RSXY(IP,1)/( DEGRAD*REARTH*COS(YP(IP)*DEGRAD) )
            D_RSYY(IP,1) = D_RSYY(IP,1)/( DEGRAD*REARTH*COS(YP(IP)*DEGRAD) )
            D_RSXX(IP,2) = D_RSXX(IP,2)/( DEGRAD*REARTH )
            D_RSXY(IP,2) = D_RSXY(IP,2)/( DEGRAD*REARTH )
            D_RSYY(IP,2) = D_RSYY(IP,2)/( DEGRAD*REARTH )
         END DO
      END IF

!     Fx = - (@Sxx/@x + @Sxy/@y) = FORCE(IP,1) 
!     Fy = - (@Sxy/@x + @Syy/@y) = FORCE(IP,2) 

      FORCEXY(:,1) = - (D_RSXX(:,1)+D_RSXY(:,2))
      FORCEXY(:,2) = - (D_RSXY(:,1)+D_RSYY(:,2))

      !write(*,*) sum(forcexy)

     END SUBROUTINE WAVEFORCE
!**********************************************************************
!*                                                                    *
!**********************************************************************

