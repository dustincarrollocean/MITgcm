#include "SMOOTH_OPTIONS.h"

      SUBROUTINE SMOOTH_FILTERVAR3D( smoothOpNb, myThid )

C     *==========================================================*
C     | SUBROUTINE smooth_filtervar3D
C     | o Routine that computes the filter variance
C     |   field associated with a diffusion operator, as part
C     |   a 3D spatial correlation operator (smooth_correld3D.F)
C     |   See Weaver and Courtier 01 for details.
C     *==========================================================*

      IMPLICIT NONE
c     == global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SMOOTH.h"

      INTEGER smoothOpNb, myThid

c     == external functions ==
      REAL*8   port_rand, port_rand_norm
      EXTERNAL PORT_RAND, PORT_RAND_NORM

c     == local variables ==
      INTEGER i,j,k, bi, bj, ii, jj, kk
      INTEGER diLoc, djLoc,  dkLoc
      INTEGER nbRand, nbt_in
      CHARACTER*( 80) fnamegeneric
      _RL smoothTmpFld (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL smoothTmpMean(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL smoothTmpVar (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS smooth3Dmask (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
c     == end of interface ==

c--   allow a mask other than maskC
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO k = 1,Nr
         DO j = 1-OLy,sNy+OLy
          DO i = 1-OLx,sNx+OLx
           IF (smooth3DmaskName(smoothOpNb)(1:5).EQ.'maskC') THEN
            smooth3Dmask(i,j,k,bi,bj) = maskC(i,j,k,bi,bj)
           ELSEIF (smooth3DmaskName(smoothOpNb)(1:5).EQ.'maskW') THEN
            smooth3Dmask(i,j,k,bi,bj) = maskW(i,j,k,bi,bj)
           ELSEIF (smooth3DmaskName(smoothOpNb)(1:5).EQ.'maskS') THEN
            smooth3Dmask(i,j,k,bi,bj) = maskS(i,j,k,bi,bj)
           ENDIF
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

c if smooth3Dfilter(smoothOpNb)=0: the filter variance field
c has been computed earlier and is already in the run directory
c so this routine does not do anything

      IF (smooth3Dfilter(smoothOpNb).NE.0) THEN

       nbt_in = smooth3Dnbt(smoothOpNb)/2

c read smoothing [i.e diffusion] operator:
       WRITE(fnamegeneric(1:80),'(1A,I3.3)')
     &       'smooth3Doperator', smoothOpNb
       CALL READ_REC_3D_RL( fnamegeneric, smoothprec, Nr, smooth3D_Kwx,
     &                      1, 1, myThid )
       CALL READ_REC_3D_RL( fnamegeneric, smoothprec, Nr, smooth3D_Kwy,
     &                      2, 1, myThid )
       CALL READ_REC_3D_RL( fnamegeneric, smoothprec, Nr, smooth3D_Kwz,
     &                      3, 1, myThid )
       CALL READ_REC_3D_RL( fnamegeneric, smoothprec, Nr, smooth3D_Kux,
     &                      4, 1, myThid )
       CALL READ_REC_3D_RL( fnamegeneric, smoothprec, Nr, smooth3D_Kvy,
     &                      5, 1, myThid )
       CALL READ_REC_3D_RL( fnamegeneric, smoothprec, Nr, smooth3D_Kuz,
     &                      6, 1, myThid )
       CALL READ_REC_3D_RL( fnamegeneric, smoothprec, Nr, smooth3D_Kvz,
     &                      7, 1, myThid )
       CALL READ_REC_3D_RL( fnamegeneric, smoothprec, Nr, smooth3D_Kuy,
     &                      8, 1, myThid )
       CALL READ_REC_3D_RL( fnamegeneric, smoothprec, Nr, smooth3D_Kvx,
     &                      9, 1, myThid )
       CALL READ_REC_3D_RL( fnamegeneric, smoothprec, Nr,
     &                      smooth3D_kappaR, 10, 1, myThid )
       CALL EXCH_XYZ_RL( smooth3D_Kwx, myThid )
       CALL EXCH_XYZ_RL( smooth3D_Kwy, myThid )
       CALL EXCH_XYZ_RL( smooth3D_Kwz, myThid )
       CALL EXCH_XYZ_RL( smooth3D_Kux, myThid )
       CALL EXCH_XYZ_RL( smooth3D_Kvy, myThid )
       CALL EXCH_XYZ_RL( smooth3D_Kuz, myThid )
       CALL EXCH_XYZ_RL( smooth3D_Kvz, myThid )
       CALL EXCH_XYZ_RL( smooth3D_Kuy, myThid )
       CALL EXCH_XYZ_RL( smooth3D_Kvx, myThid )
       CALL EXCH_XYZ_RL( smooth3D_kappaR, myThid )

c initialize filter variance field:
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         DO k=1,Nr
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            smooth3Dnorm(i,j,k,bi,bj) = 0.
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO

       IF (smooth3Dfilter(smoothOpNb).EQ.2) THEN
c compute the normalization matrix using the approximate method
c
c This method can be quite expensive -- so that the approximate
c method (see below) is usually the prefered one.
c The exact method can be used to check the accuracy
c of the approximate method results (that can be predicted).
c
c note: the exact method requires the adjoint of smooth_diff2D.F (see below)

        diLoc = 15 !int(5*smooth_L/smooth_dx)
        djLoc = 20 !int(5*smooth_L/smooth_dx)
        dkLoc = 8

        DO kk=1,dkLoc
         DO ii=1,diLoc,2
          DO jj=1,djLoc,2

           DO bj=myByLo(myThid),myByHi(myThid)
            DO bi=myBxLo(myThid),myBxHi(myThid)
             DO k=1,Nr
              DO j=1-OLy,sNy+OLy
               DO i=1-OLx,sNx+OLx
                smoothTmpFld(i,j,k,bi,bj) = 0.
               ENDDO
              ENDDO
             ENDDO
             DO k=kk,Nr,dkLoc
              DO j=jj,sNy,djLoc
               DO i=ii,sNx,diLoc
                smoothTmpFld(i,j,k,bi,bj) = 1.
               ENDDO
              ENDDO
             ENDDO
            ENDDO
           ENDDO

c note: as we go to adjoint part, we need to have 0 in overlaps
c       so we must NOT have done an exchange for smoothTmpFld

c adjoint:
           WRITE(errorMessageUnit,'(A,/,A)' )
     &      'you need to have adsmooth_diff3D compiled and then:',
     &      'uncomment the line below and comment the STOP'
           CALL ALL_PROC_DIE( myThid )
           STOP 'ABNORMAL END: S/R smooth_filtervar3D'
c          CALL adsmooth_diff3D( smoothTmpFld, nbt_in, myThid )

c division by sqrt(volume)*sqrt(volume) [1 to end adj, 1 to begin fwd]
           DO bj=myByLo(myThid),myByHi(myThid)
            DO bi=myBxLo(myThid),myBxHi(myThid)
             DO k=1,Nr
              DO j=1,sNy
               DO i=1,sNx
c division by ~sqrt(volume):
                smoothTmpFld(i,j,k,bi,bj) = smoothTmpFld(i,j,k,bi,bj)
     &                     *(recip_rA(i,j,bi,bj)*recip_drF(k))
               ENDDO
              ENDDO
             ENDDO
            ENDDO
           ENDDO

c coming out of adjoint part: overlaps are 0
c going in fwd part: we need to fill them up
           CALL EXCH_XYZ_RL( smoothTmpFld , myThid )

c fwd:
           CALL SMOOTH_DIFF3D( smoothTmpFld, nbt_in, myThid )

c convert variance to normalization factor:
           DO bj=myByLo(myThid),myByHi(myThid)
            DO bi=myBxLo(myThid),myBxHi(myThid)
             DO k=1,Nr,dkLoc
              DO j=jj,sNy,djLoc
               DO i=ii,sNx,diLoc
                IF (smoothTmpFld(i,j,k,bi,bj).NE.0.) THEN
                 smooth3Dnorm(i,j,k,bi,bj) =
     &                    1/SQRT(smoothTmpFld(i,j,k,bi,bj))
                ENDIF
               ENDDO
              ENDDO
             ENDDO
            ENDDO
           ENDDO

          ENDDO      !DO ii=1,diLoc
         ENDDO      !DO jj=1,djLoc
        ENDDO      !DO kk=1,dkLoc

       ELSEIF (smooth3Dfilter(smoothOpNb).EQ.1) THEN
c compute the normalization matrix using the approximate method

        DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)
          DO k=1,Nr
           DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
             smoothTmpMean(i,j,k,bi,bj) = 0. _d 0
             smoothTmpVar(i,j,k,bi,bj)  = 0. _d 0
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO

c initialize random number generator
        smoothTmpFld(1,1,1,1,1) = port_rand(1.d0)
        nbRand = 1000

        DO ii=1,nbRand
         WRITE(standardMessageUnit,'(A,I4,A,I4)')
     &   'smooth_filtervar3D: ', ii, ' members done out of ', nbRand

c fill smoothTmpFld with random numbers:
         DO bj=myByLo(myThid),myByHi(myThid)
          DO bi=myBxLo(myThid),myBxHi(myThid)
           DO k=1,Nr
            DO j=1-OLy,sNy+OLy
             DO i=1-OLx,sNx+OLx
              smoothTmpFld(i,j,k,bi,bj) = 0. _d 0
              IF (smooth3dmask(i,j,k,bi,bj).NE.0) THEN
               smoothTmpFld(i,j,k,bi,bj)=port_rand_norm()
              ENDIF
c division by sqrt(volume):
              smoothTmpFld(i,j,k,bi,bj) = smoothTmpFld(i,j,k,bi,bj)
     &                    *SQRT(recip_rA(i,j,bi,bj)*recip_drF(k))
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO

         CALL EXCH_XYZ_RL( smoothTmpFld, myThid )

c smooth random number field
         CALL SMOOTH_DIFF3D( smoothTmpFld, nbt_in, myThid )

c accumulate statistics (to compute the variance later)
         DO bj=myByLo(myThid),myByHi(myThid)
          DO bi=myBxLo(myThid),myBxHi(myThid)
           DO k=1,Nr
            DO j=1-OLy,sNy+OLy
             DO i=1-OLx,sNx+OLx
              smoothTmpVar(i,j,k,bi,bj) = smoothTmpVar(i,j,k,bi,bj)
     &         + smoothTmpFld(i,j,k,bi,bj)
     &          *smoothTmpFld(i,j,k,bi,bj)/nbRand
              smoothTmpMean(i,j,k,bi,bj) = smoothTmpMean(i,j,k,bi,bj)
     &         + smoothTmpFld(i,j,k,bi,bj)/nbRand
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO

C-      ii loop end
        ENDDO

c compute variance and convert it to normalization factor:
        DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)
          DO k=1,Nr
           DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
             IF (smooth3dmask(i,j,k,bi,bj).NE.0) THEN
              smooth3Dnorm(i,j,k,bi,bj) = 1/SQRT(
     &           nbRand/(nbRand-1)
     &          *( smoothTmpVar(i,j,k,bi,bj)
     &           - smoothTmpMean(i,j,k,bi,bj)*smoothTmpMean(i,j,k,bi,bj)
     &           )                              )
             ENDIF
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO

C      end smooth3Dfilter() if block
       ENDIF

c write smooth3Dnorm_3D to file:
       WRITE(fnamegeneric(1:80),'(1A,I3.3)')
     &       'smooth3Dnorm', smoothOpNb
       CALL WRITE_REC_3D_RL( fnamegeneric, smoothprec,
     &                       Nr, smooth3Dnorm, 1, 1, myThid )
       CALL EXCH_XYZ_RL( smooth3Dnorm,  myThid )

C     end smooth3Dfilter() <> 0 if block
      ENDIF

      RETURN
      END
