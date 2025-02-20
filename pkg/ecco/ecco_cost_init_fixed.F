#include "ECCO_OPTIONS.h"
#include "AD_CONFIG.h"
#ifdef ALLOW_CTRL
# include "CTRL_OPTIONS.h"
#endif

      subroutine ecco_cost_init_fixed( myThid )

c     ==================================================================
c     SUBROUTINE ecco_cost_init_fixed
c     ==================================================================
c
c     o Set contributions to the cost function and the cost function
c       itself to zero. The cost function and the individual contribu-
c       tions are defined in the header file "ecco_cost.h".
c
c     started: Christian Eckert eckert@mit.edu 30-Jun-1999
c
c     changed: Christian Eckert eckert@mit.edu 25-Feb-2000
c
c              - Restructured the code in order to create a package
c                for the MITgcmUV.
c
c     changed: Ralf Giering 18-Jan-2001
c
c              - move namelist reading to cost_readparms.F
c
c     ==================================================================
c     SUBROUTINE ecco_cost_init_fixed
c     ==================================================================

      implicit none

c     == global variables ==

#include "EEPARAMS.h"
#include "SIZE.h"
#include "GRID.h"
#include "PARAMS.h"

#ifdef ALLOW_CAL
#include "cal.h"
#endif
#ifdef ALLOW_ECCO
# ifdef ECCO_CTRL_DEPRECATED
#  include "ecco_cost.h"
# else
#  include "ECCO_SIZE.h"
#  include "ECCO.h"
# endif
#endif
#ifdef ALLOW_CTRL
# include "optim.h"
#endif

c     == routine arguments ==

      integer myThid

c     == local variables ==

      integer k
      logical exst
      _RL     missingObsFlag
      PARAMETER ( missingObsFlag = 1. _d 23 )
#if ( defined ECCO_CTRL_DEPRECATED || defined ALLOW_GENCOST_1D )
      character*(MAX_LEN_MBUF) msgBuf
      integer ilo, ihi, irec, gwunit
#endif
#if ( defined ECCO_CTRL_DEPRECATED && defined ALLOW_TRANSPORT_COST_CONTRIBUTION )
      _RL     dummy
#endif
#ifdef ALLOW_GENCOST_CONTRIBUTION
      integer k2
#endif

c     == external functions ==

      integer  cal_IntYears
      external cal_IntYears
      integer  cal_IntMonths
      external cal_IntMonths
      integer  cal_IntDays
      external cal_IntDays
      integer  ifnblnk
      external ifnblnk
      integer  ilnblnk
      external ilnblnk

c     == end of interface ==

#ifdef ALLOW_CTRL
      eccoiter=optimcycle
#else
      eccoiter=0
#endif

#ifdef ALLOW_CAL

c--   The number of monthly and daily averages generated by the
c--   current model integration.
      nyearsrec = cal_IntYears( myThid )
      nmonsrec = cal_IntMonths( myThid )
      ndaysrec = cal_IntDays( myThid )

      _BEGIN_MASTER( myThid )

#ifdef ECCO_CTRL_DEPRECATED

c--     Get the complete dates of the ...
c--     ... TMI data.
        if ( tmidatfile .ne. ' ' )
     &   call cal_FullDate( tmistartdate1,   tmistartdate2,
     &                      tmistartdate,    myThid )
c--     ... SST data.
        if ( sstdatfile .ne. ' ' )
     &   call cal_FullDate( sststartdate1,   sststartdate2,
     &                      sststartdate,    myThid )
c--     ... SSS data.
        if ( sssdatfile .ne. ' ' )
     &   call cal_FullDate( sssstartdate1,   sssstartdate2,
     &                      sssstartdate,    myThid )
c--     ... BP data.
        if ( bpdatfile .ne. ' ' )
     &   call cal_FullDate( bpstartdate1,   bpstartdate2,
     &                      bpstartdate,    myThid )
c--     ... IES data.
        if ( iesdatfile .ne. ' ' )
     &   call cal_FullDate( iesstartdate1,   iesstartdate2,
     &                      iesstartdate,    myThid )
#ifdef ALLOW_SSH_MEAN_COST_CONTRIBUTION
c--     ... mdt data.
        if ( mdtdatfile .ne. ' ' )
     &  call cal_FullDate( mdtstartdate1, mdtstartdate2,
     &                     mdtstartdate,  myThid )
c--     ... mdt data.
        if ( mdtdatfile .ne. ' ' )
     &  call cal_FullDate( mdtenddate1, mdtenddate2,
     &                     mdtenddate,  myThid )
#endif /* ALLOW_SSH_MEAN_COST_CONTRIBUTION */
c--     ... T/P data.
        if ( topexfile .ne. ' ' )
     &  call cal_FullDate( topexstartdate1, topexstartdate2,
     &                     topexstartdate,  myThid )
c--     ... ERS data.
        if ( ersfile .ne. ' ' )
     &   call cal_FullDate( ersstartdate1,   ersstartdate2,
     &                      ersstartdate,    myThid )
c--     ... GFO data.
        if ( gfofile .ne. ' ' )
     &   call cal_FullDate( gfostartdate1,   gfostartdate2,
     &                      gfostartdate,    myThid )
c--     ... SCAT data.
        if ( scatxdatfile .ne. ' ' )
     &   call cal_FullDate( scatstartdate1,   scatstartdate2,
     &                      scatxstartdate,    myThid )
        if ( scatydatfile .ne. ' ' )
     &   call cal_FullDate( scatstartdate1,   scatstartdate2,
     &                      scatystartdate,    myThid )
c--     ... ARGO data.
        if ( argotfile .ne. ' ' )
     &  call cal_FullDate( argotstartdate1,   argotstartdate2,
     &                     argotstartdate,    myThid )
         if ( argosfile .ne. ' ' )
     &  call cal_FullDate( argosstartdate1,   argotstartdate2,
     &                     argosstartdate,    myThid )
#endif /* ECCO_CTRL_DEPRECATED */

#ifdef ALLOW_GENCOST_CONTRIBUTION
      do k = 1, NGENCOST

c--    skip averaging when several cost terms use the
c      same barfile or when barfile is undefined
       gencost_barskip(k)=.FALSE.
       if ( gencost_barfile(k).EQ.' ' )
     &      gencost_barskip(k)=.TRUE.
       do k2 = 1,k-1
         if ( gencost_barfile(k2).EQ.gencost_barfile(k) )
     &      gencost_barskip(k)=.TRUE.
       enddo

c--    set time averaging parameters
       if ( (using_gencost(k)).AND.( (gencost_flag(k).GE.1).OR.
     &         (gencost_avgperiod(k).NE.'     ') ) ) then
         if ( gencost_avgperiod(k) .EQ. 'day' .OR.
     &        gencost_avgperiod(k) .EQ. 'DAY' ) then
            gencost_nrec(k)   = ndaysrec
            gencost_period(k) = 86400.
         else if ( gencost_avgperiod(k) .EQ. 'month' .OR.
     &        gencost_avgperiod(k) .EQ. 'MONTH' ) then
            gencost_nrec(k)   =nmonsrec
            gencost_period(k) = 0.
         else if ( gencost_avgperiod(k) .EQ. 'step' .OR.
     &        gencost_avgperiod(k) .EQ. 'STEP' ) then
            gencost_nrec(k)   =nTimeSteps+1
            gencost_period(k) = dTtracerLev(1)
         else if ( gencost_avgperiod(k) .EQ. 'const' .OR.
     &        gencost_avgperiod(k) .EQ. 'CONST' ) then
            gencost_nrec(k)   =1
            gencost_period(k) = dTtracerLev(1)
         else if ( gencost_avgperiod(k) .EQ. 'year' .OR.
     &        gencost_avgperiod(k) .EQ. 'YEAR' ) then
           STOP
     &      'ecco_cost_init_fixed: yearly data not yet implemented'
         else
           STOP
     &      'ecco_cost_init_fixed: gencost_avgperiod wrongly specified'
         endif
       endif

c--    set observation start/enddate
       if (gencost_startdate1(k).GT.0) then
         call cal_FullDate(
     &     gencost_startdate1(k), gencost_startdate2(k),
     &     gencost_startdate(1,k), myThid )
       else
         call cal_CopyDate(modelStartDate,
     &     gencost_startdate(1,k),myThid)
         gencost_startdate1(k)=startdate_1
         gencost_startdate2(k)=startdate_2
       endif

       if (gencost_enddate1(k).GT.0) then
         call cal_FullDate(
     &     gencost_enddate1(k), gencost_enddate2(k),
     &     gencost_enddate(1,k), myThid )
       else
         call cal_CopyDate(modelEndDate,
     &     gencost_enddate(1,k),myThid)
       endif

#ifdef ALLOW_GENCOST_1D
      if ( (gencost_name(k).EQ.'gmbp') .OR.
     &     (gencost_name(k).EQ.'gmsl')) then
       if(gencost_nrec(k).GT.N1DDATA)then
        WRITE(msgBuf,'(2A,2i8)') 'ecco_cost_init_fixed: ',
     &   'Increase N1DDATA', N1DDATA, gencost_nrec(k)
        CALL PRINT_ERROR( msgBuf, myThid)
       endif

      do irec = 1, N1DDATA
       gencost_1DDATA(irec,k) = 0. _d 0
      enddo

      if(gencost_wei1d(k).NE.0. _d 0)then
        gencost_wei1d(k) = 1. _d 0 / gencost_wei1d(k)
     &      /gencost_wei1d(k)
      endif

      if ( gencost_datafile(k) .NE. ' ' ) then
      ilo = ifnblnk(gencost_datafile(k))
      ihi = ilnblnk(gencost_datafile(k))

      CALL OPEN_COPY_DATA_FILE(
     I                          gencost_datafile(k)(ilo:ihi),
     I                          'ECCO_COST_INIT_FIXED: ',
     O                          gwunit,
     I                          myThid )
      do irec = 1, gencost_nrec(k)
         read(gwunit,*) gencost_1DDATA(irec,k)
      enddo
      close(gwunit)
      _BARRIER
      endif

      endif
#endif /* ALLOW_GENCOST_1D */

      enddo !do k = 1, NGENCOST
#endif /* ALLOW_GENCOST_CONTRIBUTION */

      _END_MASTER( myThid )

#endif /* ALLOW_CAL */

      call ecco_check( myThid )

c--   Get the weights that are to be used for the individual cost
c--   function contributions.
      call ecco_cost_weights( myThid )

c--   Initialise adjoint of monthly mean files calculated
c--   in cost_averagesfields (and their ad...).
cph(
cph   The following init. shoud not be applied if in the middle
cph   of a divided adjoint run
cph)
#ifndef ALLOW_TANGENTLINEAR_RUN
cph!!! and I think it needs to be seen by TAF
cph!!! for repeated TLM runs
cph!!!
      inquire( file='costfinal', exist=exst )
      if ( .NOT. exst) then
         call ecco_cost_init_barfiles( myThid )
      endif
#endif

#ifdef ECCO_CTRL_DEPRECATED

#ifdef ALLOW_TRANSPORT_COST_CONTRIBUTION
      do irec = 1, ndaysrec
       wtransp(irec)   = 0. _d 0
       transpobs(irec) = 0. _d 0
      enddo

      if ( costTranspDataFile .NE. ' ' ) then
      _BEGIN_MASTER(myThid)
      ilo = ifnblnk(costTranspDataFile)
      ihi = ilnblnk(costTranspDataFile)
      CALL OPEN_COPY_DATA_FILE(
     I                          costTranspDataFile(ilo:ihi),
     I                          'ECCO_COST_INIT_FIXED',
     O                          gwunit,
     I                          myThid )
      do irec = 1, ndaysrec
c-- read daily transport time series
c-- 1st: transport in m/s
c-- 2nd: date in YYYYMMDD
c-- 3rd: uncertainty in m/s
         read(gwunit,*) transpobs(irec), dummy, wtransp(irec)
c-- convert std.dev. to weight
         if ( wtransp(irec) .NE. 0. )
     &        wtransp(irec) =1.0/(wtransp(irec)*wtransp(irec))
c-- set weight to zero for missing values
         if ( transpobs(irec) .EQ. missingObsFlag )
     &        wtransp(irec) = 0. _d 0
      enddo
      _END_MASTER(myThid)
      _BARRIER
      endif
#endif /* ALLOW_TRANSPORT_COST_CONTRIBUTION */

#ifdef ALLOW_NEW_SSH_COST

c--   Read flags for picking SSH time averages
      do irec = 1, ndaysrec
       tpTimeMask(irec)  = 1. _d 0
       ersTimeMask(irec) = 1. _d 0
       gfoTimeMask(irec) = 1. _d 0
      enddo
c
      _BEGIN_MASTER(myThid)
c
#ifdef ALLOW_SSH_TPANOM_COST_CONTRIBUTION
      if ( tpTimeMaskFile .NE. ' ' ) then
       ilo = ifnblnk(tpTimeMaskFile)
       ihi = ilnblnk(tpTimeMaskFile)
       CALL OPEN_COPY_DATA_FILE(
     I                          tpTimeMaskFile(ilo:ihi),
     I                          'cost_ssh tp',
     O                          gwunit,
     I                          myThid )
       do irec = 1, ndaysrec
        read(gwunit,*) tpTimeMask(irec)
       enddo
      endif
#endif
c
#ifdef ALLOW_SSH_ERSANOM_COST_CONTRIBUTION
      if ( ersTimeMaskFile .NE. ' ' ) then
       ilo = ifnblnk(ersTimeMaskFile)
       ihi = ilnblnk(ersTimeMaskFile)
       CALL OPEN_COPY_DATA_FILE(
     I                          ersTimeMaskFile(ilo:ihi),
     I                          'cost_ssh ers',
     O                          gwunit,
     I                          myThid )
       do irec = 1, ndaysrec
        read(gwunit,*) ersTimeMask(irec)
       enddo
      endif
#endif
c
#ifdef ALLOW_SSH_GFOANOM_COST_CONTRIBUTION
      if ( gfoTimeMaskFile .NE. ' ' ) then
       ilo = ifnblnk(gfoTimeMaskFile)
       ihi = ilnblnk(gfoTimeMaskFile)
       CALL OPEN_COPY_DATA_FILE(
     I                          gfoTimeMaskFile(ilo:ihi),
     I                          'cost_ssh gfo',
     O                          gwunit,
     I                          myThid )
       do irec = 1, ndaysrec
        read(gwunit,*) gfoTimeMask(irec)
       enddo
      endif
#endif
c
       do irec = 1, ndaysrec
        if (
     &   ( tpTimeMask(irec).NE.0. .AND. tpTimeMask(irec).NE.1. ) .OR.
     &   ( ersTimeMask(irec).NE.0. .AND. ersTimeMask(irec).NE.1. ) .OR.
     &   ( ersTimeMask(irec).NE.0. .AND. ersTimeMask(irec).NE.1. ) )
     &  then
           WRITE(msgBuf,'(2A,I10)')
     &        'ecco_cost_init_fixed: (SSH)TimeMask not 0. or 1. ',
     &        'for irec (=day) ', irec
           CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &          SQUEEZE_RIGHT , myThid )
           CALL PRINT_ERROR( msgBuf , myThid )
           STOP 'ABNORMAL END: S/R ECCO_COST_INIT_FIXED'
        endif
       enddo
c
      _END_MASTER(myThid)
      _BARRIER
#endif /* ALLOW_NEW_SSH_COST */

#endif /* ECCO_CTRL_DEPRECATED */

c--   Summarize the cost function setup.
      _BEGIN_MASTER( myThid )
      call ecco_summary( myThid )
      call ecco_cost_summary( myThid )
      _END_MASTER( myThid )

      _BARRIER

      RETURN
      END
