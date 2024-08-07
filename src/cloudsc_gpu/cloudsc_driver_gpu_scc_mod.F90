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
    
    ! Temporary buffers used for double blocked lopo todo: remove and explicitly transfer
    ! copyin
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pt_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pq_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: buffer_tmp_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pvfa_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pvfl_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pvfi_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pdyna_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pdynl_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pdyni_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: phrsw_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: phrlw_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pvervel_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pap_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: paph_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: plsm_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: ldcum_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: ktype_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: plu_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: psnde_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pmfu_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pmfd_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pa_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pclv_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: psupsat_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: plcrit_aer_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: picrit_aer_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pre_ice_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pccn_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pnice_block
    ! ! copy
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) ::buffer_loc_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: plude_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pcovptot_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: prainfrac_toprfz_block
    ! ! copyout
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pfsqlf_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pfsqif_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pfcqnng_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pfcqlng_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pfsqrf_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pfsqsf_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pfcqrng_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pfcqsng_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pfsqltur_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pfsqitur_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pfplsl_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pfplsn_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pfhpsl_block
    ! REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: pfhpsn_block

    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:) :: TEST_ARRAY
    REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:) :: TEST_ARRAY_BLOCK
    INTEGER(KIND=JPIM) :: J
    INTEGER(KIND=JPIM) :: I
    INTEGER(KIND=JPIM) :: BLK
    

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


!     ! buffer allocations
!     !copyin
!     ALLOCATE pt_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE) ! T at start of callpar
!     ALLOCATE pq_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE) ! Q at start of callpar
!     ALLOCATE buffer_tmp_block(NPROMA,NLEV,3+NCLV,BUFFER_BLOCK_SIZE) ! Storage buffer for TENDENCY_TMP
!     ! ALLOCATE BUFFER_CML_block(NPROMA,NLEV,3+NCLV,BUFFER_BLOCK_SIZE) ! Storage buffer for TENDENCY_CML
!     ALLOCATE pvfa_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)     ! CC from VDF scheme
!     ALLOCATE pvfl_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)     ! Liq from VDF scheme
!     ALLOCATE pvfi_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)     ! Ice from VDF scheme
!     ALLOCATE pdyna_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)    ! CC from Dynamics
!     ALLOCATE pdynl_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)    ! Liq from Dynamics
!     ALLOCATE pdyni_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)    ! Liq from Dynamics
!     ALLOCATE phrsw_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)    ! Short-wave heating rate
!     ALLOCATE phrlw_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)    ! Long-wave heating rate
!     ALLOCATE pvervel_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)  !Vertical velocity
!     ALLOCATE pap_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)      ! Pressure on full levels
!     ALLOCATE paph_block(NPROMA, NLEV+1, BUFFER_BLOCK_SIZE) ! Pressure on half levels
!     ALLOCATE plsm_block(NPROMA, BUFFER_BLOCK_SIZE)    ! Land fraction _block(0-1)
!     ALLOCATE ldcum_block(NPROMA, BUFFER_BLOCK_SIZE)    ! Convection active
!     ALLOCATE ktype_block(NPROMA, BUFFER_BLOCK_SIZE)    ! Convection type 0,1,2
!     ALLOCATE plu_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)      ! Conv. condensate
!     ALLOCATE psnde_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)    ! Conv. detrained snow
!     ALLOCATE pmfu_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)     ! Conv. mass flux up
!     ALLOCATE pmfd_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)     ! Conv. mass flux down
!     ALLOCATE pa_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)       ! Original Cloud fraction _block(t)
!     ALLOCATE pclv_block(NPROMA, NLEV, NCLV, BUFFER_BLOCK_SIZE)
!     ALLOCATE psupsat_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)
!     ALLOCATE plcrit_aer_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)
!     ALLOCATE picrit_aer_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)
!     ALLOCATE pre_ice_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)
!     ALLOCATE pccn_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)     ! liquid cloud condensation nuclei
!     ALLOCATE pnice_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)    ! ice number concentration _block(cf. CCN)
!     
!     ! copy
!     ALLOCATE buffer_loc_block(NPROMA,NLEV,3+NCLV,BUFFER_BLOCK_SIZE) ! Storage buffer for TENDENCY_LOC
!     ALLOCATE plude_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)    ! Conv. detrained water
!     ALLOCATE pcovptot_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)    ! Precip fraction
!     ALLOCATE prainfrac_toprfz_block(NPROMA, BUFFER_BLOCK_SIZE)
!     
!     !copyout
!     ALLOCATE pfsqlf_block(NPROMA, NLEV+1, BUFFER_BLOCK_SIZE)    ! Flux of liquid
!     ALLOCATE pfsqif_block(NPROMA, NLEV+1, BUFFER_BLOCK_SIZE)    ! Flux of ice
!     ALLOCATE pfcqnng_block(NPROMA, NLEV+1, BUFFER_BLOCK_SIZE)   ! -ve corr for ice
!     ALLOCATE pfcqlng_block(NPROMA, NLEV+1, BUFFER_BLOCK_SIZE)   ! -ve corr for liq
!     ALLOCATE pfsqrf_block(NPROMA, NLEV+1, BUFFER_BLOCK_SIZE)    ! Flux diagnostics
!     ALLOCATE pfsqsf_block(NPROMA, NLEV+1, BUFFER_BLOCK_SIZE)    !    for DDH, generic
!     ALLOCATE pfcqrng_block(NPROMA, NLEV+1, BUFFER_BLOCK_SIZE)   ! rain
!     ALLOCATE pfcqsng_block(NPROMA, NLEV+1, BUFFER_BLOCK_SIZE)   ! snow
!     ALLOCATE pfsqltur_block(NPROMA, NLEV+1, BUFFER_BLOCK_SIZE)  ! liquid flux due to VDF
!     ALLOCATE pfsqitur_block(NPROMA, NLEV+1, BUFFER_BLOCK_SIZE)  ! ice flux due to VDF
!     ALLOCATE pfplsl_block(NPROMA, NLEV+1, BUFFER_BLOCK_SIZE)    ! liq+rain sedim flux
!     ALLOCATE pfplsn_block(NPROMA, NLEV+1, BUFFER_BLOCK_SIZE)    ! ice+snow sedim flux
!     ALLOCATE pfhpsl_block(NPROMA, NLEV+1, BUFFER_BLOCK_SIZE)    ! Enthalpy flux for liq
!     ALLOCATE pfhpsn_block(NPROMA, NLEV+1, BUFFER_BLOCK_SIZE)    ! ice number concentration _block(cf. CCN)


    DO BUFFER_IDX=0, BUFFER_COUNT-1
    BLOCK_START=BUFFER_IDX*BUFFER_BLOCK_SIZE+1
    BLOCK_END=MIN((BUFFER_IDX+1)*BUFFER_BLOCK_SIZE, NGPTOT)

! pt=pt(:,:, BLOCK_START:BLOCK_END)
! pq=pq(:,:,BLOCK_START:BLOCK_END)
! buffer_tmp=buffer_tmp(:,:,:,BLOCK_START:BLOCK_END)
! pvfa=pvfa(:,:,BLOCK_START:BLOCK_END)
! pvfl=pvfl(:,:,BLOCK_START:BLOCK_END)
! pvfi=pvfi(:,:,BLOCK_START:BLOCK_END)
! pdyna=pdyna(:,:,BLOCK_START:BLOCK_END)
! pdynl=pdynl(:,:,BLOCK_START:BLOCK_END)
! pdyni=pdyni(:,:,BLOCK_START:BLOCK_END)
! phrsw=phrsw(:,:,BLOCK_START:BLOCK_END)
! phrlw=phrlw(:,:,BLOCK_START:BLOCK_END)
! pvervel=pvervel(:,:,BLOCK_START:BLOCK_END)
! pap=pap(:,:,BLOCK_START:BLOCK_END)
! paph=paph(:,:,BLOCK_START:BLOCK_END)
! plsm=plsm(:,BLOCK_START:BLOCK_END)
! ldcum=ldcum(:,BLOCK_START:BLOCK_END)
! ktype=ktype(:,BLOCK_START:BLOCK_END)
! plu=plu(:,:,BLOCK_START:BLOCK_END)
! psnde=psnde(:,:,BLOCK_START:BLOCK_END)
! pmfu=pmfu(:,:,BLOCK_START:BLOCK_END)
! pmfd=pmfd(:,:,BLOCK_START:BLOCK_END)
! pa=pa(:,:,BLOCK_START:BLOCK_END)
! pclv=pclv(:,:,:,BLOCK_START:BLOCK_END)
! psupsat=psupsat(:,:,BLOCK_START:BLOCK_END)
! plcrit_aer=plcrit_aer(:,:,BLOCK_START:BLOCK_END)
! picrit_aer=picrit_aer(:,:,BLOCK_START:BLOCK_END)
! pre_ice=pre_ice(:,:,BLOCK_START:BLOCK_END)
! pccn=pccn(:,:,BLOCK_START:BLOCK_END)
! pnice=pnice(:,:,BLOCK_START:BLOCK_END)


!$acc data &
!$acc copyin( &
!$acc   pt(:,:, BLOCK_START:BLOCK_END), &
!$acc   pq(:,:,BLOCK_START:BLOCK_END), &
!$acc   buffer_cml, &
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

    ALLOCATE TEST_ARRAY(1024,1024)
    TEST_ARRAY = 37 ! do I need TEST_ARRAY(:,:)=37 ?
    ALLOCATE TEST_ARRAY_BLOCK(1024,64)

    !$acc enter data create(TEST_ARRAY_BLOCK)

    DO BLK=1,1024,64

      !$acc host_data use_device(TEST_ARRAY_BLOCK)
      call acc_memcpy_to_device(TEST_ARRAY_BLOCK, TEST_ARRAY(:,BLK:BLK+63), 1024*64*SIZEOF(TEST_ARRAY(1,1)))
      !$acc end host_data

!!      !$acc serial present(TEST_ARRAY_BLOCK)  ! Inside serial region everythin is executed on device on 1 thread
!!      !$acc end serial

      !$acc parallel loop gang vector_length(64)
      DO J=BLK,BKL+63
        !$acc loop vector(64)
        DO I=1,1024
        TEST_ARRAY(I,J) = 5
        END DO
      END DO

      !$acc host_data use_device(TEST_ARRAY_BLOCK)
      call acc_memcpy_to_device(TEST_ARRAY(:,BLK:BLK+63), TEST_ARRAY_BLOCK, 1024*64*SIZEOF(TEST_ARRAY(1,1)))
      !$acc end host_data
      !$acc exit data delete(TEST_ARRAY_BLOCK)

    END DO

    ! CHECK OUTPUT
      DO J=1,1024
        DO I=1,1024
          IF (TEST_ARRAY(I,J) /= 5) print*, 'Incorect value in TEST_ARRAY (should be 5)'
        END DO
      END DO
    DO
    DEALLOCATE(TEST_ARRAY)
    DEALLOCATE(TEST_ARRAY_BLOCK)

    ! ! deallocate buffer arrays
    ! DEALLOCATE pt_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE) ! T at start of callpar
    ! DEALLOCATE pq_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE) ! Q at start of callpar
    ! DEALLOCATE buffer_tmp_block(NPROMA,NLEV,3+NCLV,BUFFER_BLOCK_SIZE) ! Storage buffer for TENDENCY_TMP
    ! ! DEALLOCATE BUFFER_CML_block(NPROMA,NLEV,3+NCLV,BUFFER_BLOCK_SIZE) ! Storage buffer for TENDENCY_CML
    ! DEALLOCATE pvfa_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)     ! CC from VDF scheme
    ! DEALLOCATE pvfl_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)     ! Liq from VDF scheme
    ! DEALLOCATE pvfi_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)     ! Ice from VDF scheme
    ! DEALLOCATE pdyna_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)    ! CC from Dynamics
    ! DEALLOCATE pdynl_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)    ! Liq from Dynamics
    ! DEALLOCATE pdyni_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)    ! Liq from Dynamics
    ! DEALLOCATE phrsw_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)    ! Short-wave heating rate
    ! DEALLOCATE phrlw_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)    ! Long-wave heating rate
    ! DEALLOCATE pvervel_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)  !Vertical velocity
    ! DEALLOCATE pap_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)      ! Pressure on full levels
    ! DEALLOCATE paph_block(NPROMA, NLEV+1, BUFFER_BLOCK_SIZE) ! Pressure on half levels
    ! DEALLOCATE plsm_block(NPROMA, BUFFER_BLOCK_SIZE)    ! Land fraction _block(0-1)
    ! DEALLOCATE ldcum_block(NPROMA, BUFFER_BLOCK_SIZE)    ! Convection active
    ! DEALLOCATE ktype_block(NPROMA, BUFFER_BLOCK_SIZE)    ! Convection type 0,1,2
    ! DEALLOCATE plu_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)      ! Conv. condensate
    ! DEALLOCATE psnde_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)    ! Conv. detrained snow
    ! DEALLOCATE pmfu_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)     ! Conv. mass flux up
    ! DEALLOCATE pmfd_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)     ! Conv. mass flux down
    ! DEALLOCATE pa_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)       ! Original Cloud fraction _block(t)
    ! DEALLOCATE pclv_block(NPROMA, NLEV, NCLV, BUFFER_BLOCK_SIZE)
    ! DEALLOCATE psupsat_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)
    ! DEALLOCATE plcrit_aer_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)
    ! DEALLOCATE picrit_aer_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)
    ! DEALLOCATE pre_ice_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)
    ! DEALLOCATE pccn_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)     ! liquid cloud condensation nuclei
    ! DEALLOCATE pnice_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)    ! ice number concentration _block(cf. CCN)
    !
    ! ! copy
    ! DEALLOCATE buffer_loc_block(NPROMA,NLEV,3+NCLV,BUFFER_BLOCK_SIZE) ! Storage buffer for TENDENCY_LOC
    ! DEALLOCATE plude_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)    ! Conv. detrained water
    ! DEALLOCATE pcovptot_block(NPROMA, NLEV, BUFFER_BLOCK_SIZE)    ! Precip fraction
    ! DEALLOCATE prainfrac_toprfz_block(NPROMA, BUFFER_BLOCK_SIZE)
    !
    ! !copyout
    ! DEALLOCATE pfsqlf_block(NPROMA, NLEV+1, BUFFER_BLOCK_SIZE)    ! Flux of liquid
    ! DEALLOCATE pfsqif_block(NPROMA, NLEV+1, BUFFER_BLOCK_SIZE)    ! Flux of ice
    ! DEALLOCATE pfcqnng_block(NPROMA, NLEV+1, BUFFER_BLOCK_SIZE)   ! -ve corr for ice
    ! DEALLOCATE pfcqlng_block(NPROMA, NLEV+1, BUFFER_BLOCK_SIZE)   ! -ve corr for liq
    ! DEALLOCATE pfsqrf_block(NPROMA, NLEV+1, BUFFER_BLOCK_SIZE)    ! Flux diagnostics
    ! DEALLOCATE pfsqsf_block(NPROMA, NLEV+1, BUFFER_BLOCK_SIZE)    !    for DDH, generic
    ! DEALLOCATE pfcqrng_block(NPROMA, NLEV+1, BUFFER_BLOCK_SIZE)   ! rain
    ! DEALLOCATE pfcqsng_block(NPROMA, NLEV+1, BUFFER_BLOCK_SIZE)   ! snow
    ! DEALLOCATE pfsqltur_block(NPROMA, NLEV+1, BUFFER_BLOCK_SIZE)  ! liquid flux due to VDF
    ! DEALLOCATE pfsqitur_block(NPROMA, NLEV+1, BUFFER_BLOCK_SIZE)  ! ice flux due to VDF
    ! DEALLOCATE pfplsl_block(NPROMA, NLEV+1, BUFFER_BLOCK_SIZE)    ! liq+rain sedim flux
    ! DEALLOCATE pfplsn_block(NPROMA, NLEV+1, BUFFER_BLOCK_SIZE)    ! ice+snow sedim flux
    ! DEALLOCATE pfhpsl_block(NPROMA, NLEV+1, BUFFER_BLOCK_SIZE)    ! Enthalpy flux for liq
    ! DEALLOCATE pfhpsn_block(NPROMA, NLEV+1, BUFFER_BLOCK_SIZE)    ! ice number concentration _block(cf. CCN)

    CALL TIMER%THREAD_END(TID)


    CALL TIMER%END()

    ! On GPUs, adding block-level column totals is cumbersome and
    ! error prone, and of little value due to the large number of
    ! processing "thread teams". Instead we register the total here.
    CALL TIMER%THREAD_LOG(TID=TID, IGPC=NGPTOT)

    CALL TIMER%PRINT_PERFORMANCE(NPROMA, NGPBLKS, NGPTOT)

  END SUBROUTINE CLOUDSC_DRIVER_GPU_SCC

END MODULE CLOUDSC_DRIVER_GPU_SCC_MOD
