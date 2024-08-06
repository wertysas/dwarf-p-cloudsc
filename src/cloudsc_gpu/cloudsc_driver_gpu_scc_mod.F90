! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE CLOUDSC_DRIVER_GPU_SCC_MOD

  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMPHYDER, ONLY: STATE_TYPE
  USE YOECLDP, ONLY : NCLV, YRECLDP, TECLDP
  USE CLOUDSC_MPI_MOD, ONLY: NUMPROC, IRANK
  USE TIMER_MOD, ONLY : PERFORMANCE_TIMER, GET_THREAD_NUM

  USE CLOUDSC_GPU_SCC_MOD, ONLY: CLOUDSC_SCC

  IMPLICIT NONE

CONTAINS

  SUBROUTINE CLOUDSC_DRIVER_GPU_SCC( &
     & NUMOMP, NPROMA, NLEV, NGPTOT, NGPBLKS, NGPTOTG, KFLDX, PTSPHY, &
     & PT, PQ, &
     & BUFFER_CML, BUFFER_TMP, BUFFER_LOC, &
     & PVFA, PVFL, PVFI, PDYNA, PDYNL, PDYNI, &
     & PHRSW,    PHRLW, &
     & PVERVEL,  PAP,      PAPH, &
     & PLSM,     LDCUM,    KTYPE, &
     & PLU,      PLUDE,    PSNDE,    PMFU,     PMFD, &
     & PA, &
     & PCLV,     PSUPSAT,&
     & PLCRIT_AER,PICRIT_AER, PRE_ICE, &
     & PCCN,     PNICE,&
     & PCOVPTOT, PRAINFRAC_TOPRFZ, &
     & PFSQLF,   PFSQIF ,  PFCQNNG,  PFCQLNG, &
     & PFSQRF,   PFSQSF ,  PFCQRNG,  PFCQSNG, &
     & PFSQLTUR, PFSQITUR, &
     & PFPLSL,   PFPLSN,   PFHPSL,   PFHPSN &
     & )
    ! Driver routine that invokes the optimized CLAW-based CLOUDSC GPU kernel

    INTEGER(KIND=JPIM)                                    :: NUMOMP, NPROMA, NLEV, NGPTOT, NGPBLKS, NGPTOTG
    INTEGER(KIND=JPIM)                                    :: KFLDX
    REAL(KIND=JPRB)                                       :: PTSPHY       ! Physics timestep
    REAL(KIND=JPRB), INTENT(IN)    :: PT(NPROMA, NLEV, NGPBLKS) ! T at start of callpar
    REAL(KIND=JPRB), INTENT(IN)    :: PQ(NPROMA, NLEV, NGPBLKS) ! Q at start of callpar
    REAL(KIND=JPRB), INTENT(INOUT) :: BUFFER_CML(NPROMA,NLEV,3+NCLV,NGPBLKS) ! Storage buffer for TENDENCY_CML
    REAL(KIND=JPRB), INTENT(INOUT) :: BUFFER_TMP(NPROMA,NLEV,3+NCLV,NGPBLKS) ! Storage buffer for TENDENCY_TMP
    REAL(KIND=JPRB), INTENT(INOUT) :: BUFFER_LOC(NPROMA,NLEV,3+NCLV,NGPBLKS) ! Storage buffer for TENDENCY_LOC
    REAL(KIND=JPRB), INTENT(IN)    :: PVFA(NPROMA, NLEV, NGPBLKS)     ! CC from VDF scheme
    REAL(KIND=JPRB), INTENT(IN)    :: PVFL(NPROMA, NLEV, NGPBLKS)     ! Liq from VDF scheme
    REAL(KIND=JPRB), INTENT(IN)    :: PVFI(NPROMA, NLEV, NGPBLKS)     ! Ice from VDF scheme
    REAL(KIND=JPRB), INTENT(IN)    :: PDYNA(NPROMA, NLEV, NGPBLKS)    ! CC from Dynamics
    REAL(KIND=JPRB), INTENT(IN)    :: PDYNL(NPROMA, NLEV, NGPBLKS)    ! Liq from Dynamics
    REAL(KIND=JPRB), INTENT(IN)    :: PDYNI(NPROMA, NLEV, NGPBLKS)    ! Liq from Dynamics
    REAL(KIND=JPRB), INTENT(IN)    :: PHRSW(NPROMA, NLEV, NGPBLKS)    ! Short-wave heating rate
    REAL(KIND=JPRB), INTENT(IN)    :: PHRLW(NPROMA, NLEV, NGPBLKS)    ! Long-wave heating rate
    REAL(KIND=JPRB), INTENT(IN)    :: PVERVEL(NPROMA, NLEV, NGPBLKS)  !Vertical velocity
    REAL(KIND=JPRB), INTENT(IN)    :: PAP(NPROMA, NLEV, NGPBLKS)      ! Pressure on full levels
    REAL(KIND=JPRB), INTENT(IN)    :: PAPH(NPROMA, NLEV+1, NGPBLKS) ! Pressure on half levels
    REAL(KIND=JPRB), INTENT(IN)    :: PLSM(NPROMA, NGPBLKS)    ! Land fraction (0-1)
    LOGICAL, INTENT(IN)            :: LDCUM(NPROMA, NGPBLKS)    ! Convection active
    INTEGER(KIND=JPIM), INTENT(IN) :: KTYPE(NPROMA, NGPBLKS)    ! Convection type 0,1,2
    REAL(KIND=JPRB), INTENT(IN)    :: PLU(NPROMA, NLEV, NGPBLKS)      ! Conv. condensate
    REAL(KIND=JPRB), INTENT(INOUT) :: PLUDE(NPROMA, NLEV, NGPBLKS)    ! Conv. detrained water
    REAL(KIND=JPRB), INTENT(IN)    :: PSNDE(NPROMA, NLEV, NGPBLKS)    ! Conv. detrained snow
    REAL(KIND=JPRB), INTENT(IN)    :: PMFU(NPROMA, NLEV, NGPBLKS)     ! Conv. mass flux up
    REAL(KIND=JPRB), INTENT(IN)    :: PMFD(NPROMA, NLEV, NGPBLKS)     ! Conv. mass flux down
    REAL(KIND=JPRB), INTENT(IN)    :: PA(NPROMA, NLEV, NGPBLKS)       ! Original Cloud fraction (t)
    REAL(KIND=JPRB), INTENT(IN)    :: PCLV(NPROMA, NLEV, NCLV, NGPBLKS)
    REAL(KIND=JPRB), INTENT(IN)    :: PSUPSAT(NPROMA, NLEV, NGPBLKS)
    REAL(KIND=JPRB), INTENT(IN)    :: PLCRIT_AER(NPROMA, NLEV, NGPBLKS)
    REAL(KIND=JPRB), INTENT(IN)    :: PICRIT_AER(NPROMA, NLEV, NGPBLKS)
    REAL(KIND=JPRB), INTENT(IN)    :: PRE_ICE(NPROMA, NLEV, NGPBLKS)
    REAL(KIND=JPRB), INTENT(IN)    :: PCCN(NPROMA, NLEV, NGPBLKS)     ! liquid cloud condensation nuclei
    REAL(KIND=JPRB), INTENT(IN)    :: PNICE(NPROMA, NLEV, NGPBLKS)    ! ice number concentration (cf. CCN)

    REAL(KIND=JPRB), INTENT(INOUT) :: PCOVPTOT(NPROMA, NLEV, NGPBLKS)    ! Precip fraction
    REAL(KIND=JPRB), INTENT(OUT) :: PRAINFRAC_TOPRFZ(NPROMA, NGPBLKS)
    ! Flux diagnostics for DDH budget
    REAL(KIND=JPRB), INTENT(OUT) :: PFSQLF(NPROMA, NLEV+1, NGPBLKS)    ! Flux of liquid
    REAL(KIND=JPRB), INTENT(OUT) :: PFSQIF(NPROMA, NLEV+1, NGPBLKS)    ! Flux of ice
    REAL(KIND=JPRB), INTENT(OUT) :: PFCQLNG(NPROMA, NLEV+1, NGPBLKS)   ! -ve corr for liq
    REAL(KIND=JPRB), INTENT(OUT) :: PFCQNNG(NPROMA, NLEV+1, NGPBLKS)   ! -ve corr for ice
    REAL(KIND=JPRB), INTENT(OUT) :: PFSQRF(NPROMA, NLEV+1, NGPBLKS)    ! Flux diagnostics
    REAL(KIND=JPRB), INTENT(OUT) :: PFSQSF(NPROMA, NLEV+1, NGPBLKS)    !    for DDH, generic
    REAL(KIND=JPRB), INTENT(OUT) :: PFCQRNG(NPROMA, NLEV+1, NGPBLKS)   ! rain
    REAL(KIND=JPRB), INTENT(OUT) :: PFCQSNG(NPROMA, NLEV+1, NGPBLKS)   ! snow
    REAL(KIND=JPRB), INTENT(OUT) :: PFSQLTUR(NPROMA, NLEV+1, NGPBLKS)  ! liquid flux due to VDF
    REAL(KIND=JPRB), INTENT(OUT) :: PFSQITUR(NPROMA, NLEV+1, NGPBLKS)  ! ice flux due to VDF
    REAL(KIND=JPRB), INTENT(OUT) :: PFPLSL(NPROMA, NLEV+1, NGPBLKS)    ! liq+rain sedim flux
    REAL(KIND=JPRB), INTENT(OUT) :: PFPLSN(NPROMA, NLEV+1, NGPBLKS)    ! ice+snow sedim flux
    REAL(KIND=JPRB), INTENT(OUT) :: PFHPSL(NPROMA, NLEV+1, NGPBLKS)    ! Enthalpy flux for liq
    REAL(KIND=JPRB), INTENT(OUT) :: PFHPSN(NPROMA, NLEV+1, NGPBLKS)    ! ice number concentration (cf. CCN)

    INTEGER(KIND=JPIM) :: JKGLO,IBL,ICEND
    TYPE(PERFORMANCE_TIMER) :: TIMER
    INTEGER(KIND=JPIM) :: TID ! thread id from 0 .. NUMOMP - 1

    INTEGER(KIND=JPIM) :: BUFFER_BLOCK_SIZE     ! block size for blocks in outer loop /johan
    INTEGER(KIND=JPIM) :: BUFFER_COUNT          ! number of buffers
    INTEGER(KIND=JPIM) :: BUFFER_IDX            ! idx of current buffer
    INTEGER(KIND=JPIM) :: BLOCK_START            ! idx of current buffer
    INTEGER(KIND=JPIM) :: BLOCK_END            ! idx of current buffer

    ! Local copy of cloud parameters for offload
    TYPE(TECLDP) :: LOCAL_YRECLDP

    NGPBLKS = (NGPTOT / NPROMA) + MIN(MOD(NGPTOT,NPROMA), 1)
1003 format(5x,'NUMPROC=',i0,', NUMOMP=',i0,', NGPTOTG=',i0,', NPROMA=',i0,', NGPBLKS=',i0)
    if (irank == 0) then
      write(0,1003) NUMPROC,NUMOMP,NGPTOTG,NPROMA,NGPBLKS
    end if

    ! Global timer for the parallel region
    CALL TIMER%START(NUMOMP)

    ! Workaround for PGI / OpenACC oddities:
    ! Create a local copy of the parameter struct to ensure they get
    ! moved to the device the in ``acc data`` clause below
    LOCAL_YRECLDP = YRECLDP

    ! Local timer for each thread
    TID = GET_THREAD_NUM()
    CALL TIMER%THREAD_START(TID)

    BUFFER_BLOCK_SIZE=NGPTOT
    BUFFER_COUNT=(NGPTOT+BUFFER_BLOCK_SIZE-1)/BUFFER_BLOCK_SIZE

    DO BUFFER_IDX=0, BUFFER_COUNT-1
    BLOCK_START=BUFFER_IDX*BUFFER_BLOCK_SIZE+1
    BLOCK_END=MIN((BUFFER_IDX+1)*BUFFER_BLOCK_SIZE, NGPTOT)

!$acc data &
!$acc copyin( &
!$acc   pt(:,:, BLOCK_START:BLOCK_END), &
!$acc   pq(:,:,BLOCK_START:BLOCK_END),buffer_cml, &
!$acc   buffer_tmp(:,:,:,BLOCK_START:BLOCK_END), &
!$acc   pvfa(:,:,BLOCK_START:BLOCK_END), &
!$acc   pvfl(:,:,BLOCK_START:BLOCK_END), &
!$acc   pvfi(:,:,BLOCK_START:BLOCK_END), &
!$acc   pdyna(:,:,BLOCK_START:BLOCK_END), &
!$acc   pdynl(:,:,BLOCK_START:BLOCK_END), &
!$acc   pdyni(:,:,BLOCK_START:BLOCK_END),&
!$acc   phrsw(:,:,BLOCK_START:BLOCK_END), &
!$acc   phrlw(:,:,BLOCK_START:BLOCK_END), &
!$acc   pvervel(:,:,BLOCK_START:BLOCK_END), &
!$acc   pap(:,:,BLOCK_START:BLOCK_END), &
!$acc   paph(:,:,BLOCK_START:BLOCK_END), &
!$acc   plsm(:,BLOCK_START:BLOCK_END), &
!$acc   ldcum(:,BLOCK_START:BLOCK_END), &
!$acc   ktype(:,BLOCK_START:BLOCK_END), &
!$acc   plu(:,:,BLOCK_START:BLOCK_END), &
!$acc   psnde(:,:,BLOCK_START:BLOCK_END),
!$acc   pmfu(:,:,BLOCK_START:BLOCK_END), &
!$acc   pmfd(:,:,BLOCK_START:BLOCK_END), &
!$acc   pa(:,:,BLOCK_START:BLOCK_END), &
!$acc   pclv(:,:,:,BLOCK_START:BLOCK_END), &
!$acc   psupsat(:,:,BLOCK_START:BLOCK_END), &
!$acc   plcrit_aer(:,:,BLOCK_START:BLOCK_END), &
!$acc   picrit_aer(:,:,BLOCK_START:BLOCK_END), &
!$acc   pre_ice(:,:,BLOCK_START:BLOCK_END), &
!$acc   pccn(:,:,BLOCK_START:BLOCK_END), &
!$acc   pnice(:,:,BLOCK_START:BLOCK_END), &
!$acc   yrecldp) &
!$acc copy( &               ! initialized and copied to device then back to host after region is done
!$acc   buffer_loc(:,:,:,BLOCK_START:BLOCK_END), &
!$acc   plude(:,:,BLOCK_START:BLOCK_END), &
!$acc   pcovptot(:,:,BLOCK_START:BLOCK_END), &
!$acc   prainfrac_toprfz(:,BLOCK_START:BLOCK_END)) &
!$acc copyout( &
!$acc   pfsqlf(:,:,BLOCK_START:BLOCK_END), &
!$acc   pfsqif(:,:,BLOCK_START:BLOCK_END), &
!$acc   pfcqnng(:,:,BLOCK_START:BLOCK_END), &
!$acc   pfcqlng(:,:,BLOCK_START:BLOCK_END) , &
!$acc   pfsqrf(:,:,BLOCK_START:BLOCK_END), &
!$acc   pfsqsf(:,:,BLOCK_START:BLOCK_END), &
!$acc   pfcqrng(:,:,BLOCK_START:BLOCK_END), &
!$acc   pfcqsng(:,:,BLOCK_START:BLOCK_END), &
!$acc   pfsqltur(:,:,BLOCK_START:BLOCK_END), &
!$acc   pfsqitur(:,:,BLOCK_START:BLOCK_END), &
!$acc   pfplsl(:,:,BLOCK_START:BLOCK_END), &
!$acc   pfplsn(:,:,BLOCK_START:BLOCK_END), &
!$acc   pfhpsl(:,:,BLOCK_START:BLOCK_END), &
!$acc   pfhpsn(:,:,BLOCK_START:BLOCK_END))


!$acc parallel loop gang vector_length(NPROMA)
    DO JKGLO=BLOCK_START, BLOCK_END, NPROMA ! loops from 1 ... NGPTOT, with step size NPROMA
       IBL=(JKGLO-1)/NPROMA+1
       ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)

       CALL CLOUDSC_SCC &
        & (1, ICEND, NPROMA, NLEV, PTSPHY,&
        & PT(:,:,IBL), PQ(:,:,IBL), &
        & BUFFER_TMP(:,:,1,IBL), BUFFER_TMP(:,:,3,IBL), BUFFER_TMP(:,:,2,IBL), BUFFER_TMP(:,:,4:8,IBL), &
        & BUFFER_LOC(:,:,1,IBL), BUFFER_LOC(:,:,3,IBL), BUFFER_LOC(:,:,2,IBL), BUFFER_LOC(:,:,4:8,IBL), &
        & PVFA(:,:,IBL), PVFL(:,:,IBL), PVFI(:,:,IBL), PDYNA(:,:,IBL), PDYNL(:,:,IBL), PDYNI(:,:,IBL), &
        & PHRSW(:,:,IBL),    PHRLW(:,:,IBL),&
        & PVERVEL(:,:,IBL),  PAP(:,:,IBL),      PAPH(:,:,IBL),&
        & PLSM(:,IBL),       LDCUM(:,IBL),      KTYPE(:,IBL), &
        & PLU(:,:,IBL),      PLUDE(:,:,IBL),    PSNDE(:,:,IBL),    PMFU(:,:,IBL),     PMFD(:,:,IBL),&
        !---prognostic fields
        & PA(:,:,IBL),       PCLV(:,:,:,IBL),   PSUPSAT(:,:,IBL),&
        !-- arrays for aerosol-cloud interactions
        & PLCRIT_AER(:,:,IBL),PICRIT_AER(:,:,IBL),&
        & PRE_ICE(:,:,IBL),&
        & PCCN(:,:,IBL),     PNICE(:,:,IBL),&
        !---diagnostic output
        & PCOVPTOT(:,:,IBL), PRAINFRAC_TOPRFZ(:,IBL),&
        !---resulting fluxes
        & PFSQLF(:,:,IBL),   PFSQIF (:,:,IBL),  PFCQNNG(:,:,IBL),  PFCQLNG(:,:,IBL),&
        & PFSQRF(:,:,IBL),   PFSQSF (:,:,IBL),  PFCQRNG(:,:,IBL),  PFCQSNG(:,:,IBL),&
        & PFSQLTUR(:,:,IBL), PFSQITUR (:,:,IBL), &
        & PFPLSL(:,:,IBL),   PFPLSN(:,:,IBL),   PFHPSL(:,:,IBL),   PFHPSN(:,:,IBL),&
        & YRECLDP=LOCAL_YRECLDP)

    ENDDO
!$acc end parallel loop
!$acc end data

    ENDDO ! end of outer block loop

    CALL TIMER%THREAD_END(TID)


    CALL TIMER%END()

    ! On GPUs, adding block-level column totals is cumbersome and
    ! error prone, and of little value due to the large number of
    ! processing "thread teams". Instead we register the total here.
    CALL TIMER%THREAD_LOG(TID=TID, IGPC=NGPTOT)

    CALL TIMER%PRINT_PERFORMANCE(NPROMA, NGPBLKS, NGPTOT)

  END SUBROUTINE CLOUDSC_DRIVER_GPU_SCC

END MODULE CLOUDSC_DRIVER_GPU_SCC_MOD
