!------------------------------------------------------------------------------
#include "MAPL_Generic.h"
!
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: AdvCore_GridCompMod
!
! !DESCRIPTION: 
!    This a MAPL component that can be used in
!    either with offline or online applications to advect an arbitrary set
!    of constituents.
!
! \paragraph{Scientific Description:}
!
!   The advection scheme used is that from the FVdycore grid-point
!   dynamical core.  It runs on a sphere and uses finite-volume
!   discretization techniques. The advection is time split into a
!   horizontal phase that is assumed to be vertically Lagrangian and a
!   vertical remap phase. A complete description of the core from
!   which this component is taken may be found in:
!
!   \begin{quote}
!   Lin, S.-J. 2004, A vertically Lagrangian Finite-Volume Dynamical 
!   Core for Global Models. {\em Mon. Wea. Rev.}, {\bf 132}, 2293-2307.
!   \end{quote}
!
!  \paragraph{Code Implementation:}
!
!    It code uses the MAPL (http://MAPLCode.org/maplwiki/) to
!    encapsulate the FV advection scheme as an ESMF gridded component
!    using the ESMF paradigm of initialize, run and finalize methods,
!    and their SetServices routine. As in all ESMF codes, only
!    SetServices is public and the interface consists of of a Clock
!    and Import and Export states.  The import state includes a
!    specialized description of the motion field in terms of C-grid
!    winds and mass fluxes. These are assumed to have been accumulated
!    over the time interval specified in the resource file. The
!    default of this interval is 1800 seconds. The layer pressure
!    thicknesses in the import state are assumed to be the
!    instantaneous values valid at the beginning of this interval.  If
!    these thicknesses are friendly they will be updated to values
!    valid at the end of the interval, consistent with the given
!    motion field.  Mixing ratios of the constituents to be advected
!    are placed ESMF Fields within an ESMF Bundle in the Import
!    state. Each Field in the Bundle is tested for ``Friendliness'' to
!    advection; if friendly it is advected and its values updated.
!
!    Currently no Export capability is implemented. 
!
! !INTERFACE:

module AdvCore_GridCompMod

!
! !USES:

      use ESMF
      use MAPL
      use m_set_eta,       only: set_eta
      use fv_arrays_mod,   only: fv_atmos_type, FVPRC, REAL4, REAL8
      use fms_mod,         only: fms_init, set_domain, nullify_domain
      use fv_control_mod,  only: fv_init1, fv_init2, fv_end
      use fv_tracer2d_mod, only: offline_tracer_advection
      use fv_mp_mod,       only: is,ie, js,je, is_master, tile
      use fv_grid_utils_mod, only: g_sum
#ifdef ADJOINT
      ! discrete_linear: TAF-generated exact adjoint of the step (fvdycore/model/adjoint)
      use adv_adj_discrete_mod, only: adv_adj_discrete_init, adv_adj_discrete_adj, adv_adj_discrete_fwd
      use adj_backend_mod, only: backend_set_halo_ad, backend_halo_test
#endif

      USE FV_StateMod,     only: AdvCoreTracers => T_TRACERS
      USE FV_StateMod,     only: FV_Atm
      use CubeGridPrototype, only: register_grid_and_regridders

      implicit none
      private

      integer     :: nx, ny
      integer     :: npes_x, npes_y
      real(FVPRC) :: dt
      logical     :: FV3_DynCoreIsRunning=.false.
      integer     :: AdvCore_Advection
      integer     :: Use_Total_Air_Pressure
      logical     :: chk_mass
      logical, parameter :: ADVCORE_ADJ_DEBUG = .false.
      ! .true. only for MODEL_PHASE=ADJOINT in an ADJOINT build
      logical                    :: isAdjoint=.false.
#ifdef ADJOINT
      character(len=ESMF_MAXSTR) :: modelPhase
      character(len=ESMF_MAXSTR) :: advAdjMethod='reverse_wind_pure'
      character(len=ESMF_MAXSTR) :: advAdjWeight='airmass'
      integer                    :: advAdjRescale=1
      ! ADV_ADJ_METHOD: discrete_linear (fvdycore/model/adjoint, gchp/adv_adj_discrete.F90)
      logical                    :: advAdjDiscrete=.false.
      character(len=ESMF_MAXSTR) :: advAdjDiscFixer='linearized'
      integer                    :: advAdjDiscCheck=0
      character(len=ESMF_MAXSTR) :: advAdjHaloAd='own'
      ! forward tracer whose pre-advection field the forward run exports as TRAJ_PREADV
      character(len=ESMF_MAXSTR) :: advAdjTrajTracer='SPC_CO2'
      logical                    :: advAdjDiscReady=.false.
      ! ADV_ADJ_DISCRETE_CHECK 2: forward trajectory of the previous adjoint call (q at T+dt)
      REAL(FVPRC), ALLOCATABLE, SAVE :: adjQprev(:,:,:)
#endif

      integer,  parameter :: ntiles_per_pe = 1

! Tracer I/O History stuff
! -------------------------------------
      integer, parameter         :: ntracers=11
      integer                    :: ntracer
      character(len=ESMF_MAXSTR) :: myTracer
      character(len=ESMF_MAXSTR) :: tMassStr
      real(FVPRC), SAVE          :: TMASS0(ntracers)
      real(REAL8), SAVE          :: MASS0
      logical    , SAVE          :: firstRun=.true.

! !PUBLIC MEMBER FUNCTIONS:

      public SetServices
      logical, allocatable, save :: grids_on_my_pe(:)

      private scale_tracers_by_pressure_ratio
      private global_integral
      private global_integral_trmass_from_p
      private global_integral_trmass_from_dp
      private global_integral_dp
      private global_integral_dp_rst
      private global_integral_vv
!EOP

!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------
!BOP
! !IROUTINE: SetServices - Externally visible registration routine
!
! !INTERFACE:
!
      subroutine SetServices(GC, rc)
!
! !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: GC
      integer, optional,   intent(  out) :: RC
!
! !DESCRIPTION:
!
!     User-supplied setservices routine.
!     The register routine sets the subroutines to be called
!     as the init, run, and finalize routines.  Note that those are
!     private to the module.
!
!EOP

      character(len=ESMF_MAXSTR)              :: IAm
      integer                                 :: STATUS
      character(len=ESMF_MAXSTR)              :: COMP_NAME
      type (MAPL_MetaComp),      pointer      :: MAPL
      character(len=ESMF_MAXSTR)              :: DYCORE
      type(ESMF_VM)                           :: VM
      integer                                 :: comm, ndt
      integer                                 :: p_split=1
      integer                                 :: Check_Mass_Conservation
!=============================================================================

! Begin...

      ! Get my name and set-up traceback handle
      ! ---------------------------------------
    
      call ESMF_GridCompGet( GC, NAME=COMP_NAME, vm=vm, RC=STATUS )
      VERIFY_(STATUS)
      Iam = trim(COMP_NAME) // 'SetServices'

!BOS

! !IMPORT STATE:
!
    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'MFX',                                       &
         LONG_NAME  = 'pressure_weighted_eastward_mass_flux',      &
         UNITS      = 'Pa m+2 s-1',                                &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'MFY',                                       &
         LONG_NAME  = 'pressure_weighted_northward_mass_flux',     &
         UNITS      = 'Pa m+2 s-1',                                &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'CX',                                        &
         LONG_NAME  = 'eastward_accumulated_courant_number',       &
         UNITS      = '',                                          &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'CY',                                        &
         LONG_NAME  = 'northward_accumulated_courant_number',      &
         UNITS      = '',                                          &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'PLE0',                                      &
         LONG_NAME  = 'pressure_at_layer_edges_before_advection',  &
         UNITS      = 'Pa',                                        &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'PLE1',                                      &
         LONG_NAME  = 'pressure_at_layer_edges_after_advection',   &
         UNITS      = 'Pa',                                        &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec( gc,                                   &
        SHORT_NAME = 'DELPDRY',                                    &
        LONG_NAME  = 'delta dry pressure across levels',           &
        UNITS      = 'Pa',                                         &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec( gc,                                   &
        SHORT_NAME = 'TRADV',                                      &
        LONG_NAME  = 'advected_quantities',                        &
        UNITS      = 'unknown',                                    &
        DATATYPE   = MAPL_BundleItem,                              &
        RC=STATUS  )
     VERIFY_(STATUS)

    ! For using dry pressure
    !-----------------------------------------------
    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'DryPLE0',                                   &
         LONG_NAME  = 'dry_pressure_at_layer_edges_before_advection',&
         UNITS      = 'Pa',                                        &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'DryPLE1',                                   &
         LONG_NAME  = 'dry_pressure_at_layer_edges_after_advection',&
         UNITS      = 'Pa',                                        &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,             RC=STATUS  )
     VERIFY_(STATUS)

    ! For using total pressure
    !-----------------------------------------------
    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'SPHU0',                                     &
         LONG_NAME  = 'specific_humidity_before_advection',        &
         UNITS      = 'kg kg-1',                                   &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,           RC=STATUS  )
    VERIFY_(STATUS)

#ifdef ADJOINT
       call MAPL_AddImportSpec ( gc,                                  &
          SHORT_NAME = 'ADV_AIRDEN',                                &
          LONG_NAME  = 'dry_air_density_for_advection',             &
          UNITS      = 'kg m-3',                                    &
          PRECISION  = ESMF_KIND_R8,                                &
          DIMS       = MAPL_DimsHorzVert,                           &
          VLOCATION  = MAPL_VLocationCenter,           RC=STATUS  )
       VERIFY_(STATUS)
#endif

  !EXPORT STATE:
     call MAPL_AddExportSpec ( gc,                                  &
          SHORT_NAME = 'AREA',                                      &
          LONG_NAME  = 'agrid_cell_area',                           &
          UNITS      = 'm+2'  ,                                     &
          PRECISION  = ESMF_KIND_R8,                                &
          DIMS       = MAPL_DimsHorzOnly,                           &
          VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

     ! GCHP: moist pressure export to pass to GEOS-Chem
     call MAPL_AddExportSpec ( gc,                                  &
          SHORT_NAME = 'PLE',                                       &
          LONG_NAME  = 'pressure_at_layer_edges',                   &
          UNITS      = 'Pa'   ,                                     &
          PRECISION  = ESMF_KIND_R8,                                &
          DIMS       = MAPL_DimsHorzVert,                           &
          VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

     ! GCHP: dry pressure export to pass to GEOS-Chem
     call MAPL_AddExportSpec ( gc,                                  &
          SHORT_NAME = 'DryPLE',                                    &
          LONG_NAME  = 'dry_pressure_at_layer_edges',               &
          UNITS      = 'Pa'   ,                                     &
          PRECISION  = ESMF_KIND_R8,                                &
          DIMS       = MAPL_DimsHorzVert,                           &
          VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

     ! GCHP: for diagnostics
     call MAPL_AddExportSpec ( gc,                                  &
          SHORT_NAME = 'PLEadv',                                    &
          LONG_NAME  = 'post_advection_pressure_at_layer_edges',    &
          UNITS      = 'Pa'   ,                                     &
          PRECISION  = ESMF_KIND_R8,                                &
          DIMS       = MAPL_DimsHorzVert,                           &
          VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

#ifdef ADJOINT
     ! Trajectory for the discrete_linear adjoint: the forward tracer ADV_ADJ_TRAJ_TRACER
     ! right before the forward advection of the step (after the start-up rescaling).
     ! GEOS-Chem imports it, so it is written to gcadj_import_checkpoint.<T> and restored
     ! with it in the adjoint run.
     call MAPL_AddExportSpec ( gc,                                  &
          SHORT_NAME = 'TRAJ_PREADV',                               &
          LONG_NAME  = 'tracer_before_advection_for_adjoint',       &
          UNITS      = 'mol mol-1',                                 &
          PRECISION  = ESMF_KIND_R8,                                &
          DIMS       = MAPL_DimsHorzVert,                           &
          VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)
#endif

     ! GCHP: The below three exports are the same as above  but REAL4
     ! This is a tempoary work-around for a bug in MAPL 2.55 preventing
     ! conversion from R8 to R4 in History (ewl, 5/27/25)
     call MAPL_AddExportSpec ( gc,                                  &
          SHORT_NAME = 'PLE_R4',                                    &
          LONG_NAME  = 'pressure_at_layer_edges',                   &
          UNITS      = 'Pa'   ,                                     &
          DIMS       = MAPL_DimsHorzVert,                           &
          VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)
     call MAPL_AddExportSpec ( gc,                                  &
          SHORT_NAME = 'DryPLE_R4',                                 &
          LONG_NAME  = 'dry_pressure_at_layer_edges',               &
          UNITS      = 'Pa'   ,                                     &
          DIMS       = MAPL_DimsHorzVert,                           &
          VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)
     call MAPL_AddExportSpec ( gc,                                  &
          SHORT_NAME = 'PLEadv_R4',                                 &
          LONG_NAME  = 'post_advection_pressure_at_layer_edges',    &
          UNITS      = 'Pa'   ,                                     &
          DIMS       = MAPL_DimsHorzVert,                           &
          VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

     ! 3D Tracers
     do ntracer=1,ntracers
        write(myTracer, "('TEST_TRACER',i5.5)") ntracer-1
        call MAPL_AddExportSpec ( gc,                             &
             SHORT_NAME = TRIM(myTracer),                         &
             LONG_NAME  = TRIM(myTracer),                         &
             UNITS      = '1',                                    &
             DIMS       = MAPL_DimsHorzVert,                      &
             VLOCATION  = MAPL_VLocationCenter,               RC=STATUS  )
        VERIFY_(STATUS)
     enddo

!EOS

      ! Set the Profiling timers
      !-------------------------
      call MAPL_TimerAdd(GC,    name="INITIALIZE"  ,RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_TimerAdd(GC,    name="RUN"         ,RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_TimerAdd(GC,    name="FINALIZE"    ,RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_TimerAdd(GC,    name="TOTAL"       ,RC=STATUS)
      VERIFY_(STATUS)


      ! Register methods with MAPL
      ! --------------------------
      call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE,  &
                                        Initialize, RC=status )
      VERIFY_(STATUS)
      call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,         &
                                        Run,        RC=status )
      VERIFY_(STATUS)
      call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_FINALIZE,    &
                                        Finalize,   RC=status)
      VERIFY_(STATUS)

      ! Check if AdvCore is running without FV3_DynCoreIsRunning. If yes 
      ! then setup the MAPL Grid 
      ! -----------------------------------------------------------------
      call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
      VERIFY_(STATUS)
      call MAPL_GetResource(MAPL,       &
                            DYCORE,     &
                            'DYCORE:',  &
                            default="", &
                            RC=STATUS )
      VERIFY_(STATUS)
      call MAPL_GetResource(MAPL,                       &
                            AdvCore_Advection ,         &
                            label='AdvCore_Advection:', &
                            default=1,                  &
                            RC=STATUS )
      VERIFY_(STATUS)
      if(adjustl(DYCORE)=="FV3") FV3_DynCoreIsRunning = .true.
      if(adjustl(DYCORE)=="FV3+ADV") FV3_DynCoreIsRunning = .true.

      ! Check if using total air pressure. If not then use dry air pressure.
      ! 1 = total air pressure; 0 = dry air pressure
      ! -----------------------------------------------------------------
      call MAPL_GetResource(MAPL,                                         &
                            Use_Total_Air_Pressure,                       &
                            label='USE_TOTAL_AIR_PRESSURE_IN_ADVECTION:', &
                            default=0,                                    &
                            RC=STATUS )
      VERIFY_(STATUS)

      ! Check if printing mass for mass conservation check
      ! 1 = true; 0 = false
      ! -----------------------------------------------------------------
      call MAPL_GetResource(MAPL,                                         &
                            Check_Mass_Conservation,                      &
                            label='PRINT_MASS_IN_ADVECTION:', &
                            default=0,                                    &
                            RC=STATUS )
      chk_mass=.FALSE.
      if (Check_Mass_Conservation>0) chk_mass=.TRUE.

      ! Start up FMS/MPP
      !-------------------------------------------
      call ESMF_VMGet(VM,mpiCommunicator=comm,rc=STATUS)
      VERIFY_(STATUS)
      call fms_init(comm)
      VERIFY_(STATUS)

      if (.NOT. FV3_DynCoreIsRunning) then

         ! Make sure FV3 is setup
         ! -----------------------
         call register_grid_and_regridders()
         call fv_init1(FV_Atm, dt, grids_on_my_pe, p_split)

         ! Get Domain decomposition
         !-------------------------
         call MAPL_GetResource( MAPL,      &
                                nx,        &
                                'NX:',     &
                                default=0, &
                                RC=STATUS )
         VERIFY_(STATUS)
         FV_Atm(1)%layout(1) = nx
         call MAPL_GetResource( MAPL,      &
                                ny,        &
                                'NY:',     &
                                default=0, &
                                RC=STATUS )
         VERIFY_(STATUS)
         if (FV_Atm(1)%flagstruct%grid_type == 4) then
            FV_Atm(1)%layout(2) = ny
         else
            FV_Atm(1)%layout(2) = ny / 6
         end if

         ! Get Resolution Information
         !---------------------------
         ! FV grid dimensions setup from MAPL
         call MAPL_GetResource( MAPL,                     &
                                FV_Atm(1)%flagstruct%npx, &
                                'IM:',                    &
                                default=32,               &
                                RC=STATUS )
         VERIFY_(STATUS)
         call MAPL_GetResource( MAPL,                     &
                                FV_Atm(1)%flagstruct%npy, &
                                'JM:',                    &
                                default=192,              &
                                RC=STATUS )
         VERIFY_(STATUS)
         call MAPL_GetResource( MAPL,                     &
                                FV_Atm(1)%flagstruct%npz, &
                                'LM:',                    &
                                default=72,               &
                                RC=STATUS )
         VERIFY_(STATUS)

         ! FV likes npx;npy in terms of cell vertices
         if (FV_Atm(1)%flagstruct%npy == 6*FV_Atm(1)%flagstruct%npx) then
            FV_Atm(1)%flagstruct%ntiles = 6
            FV_Atm(1)%flagstruct%npy    = FV_Atm(1)%flagstruct%npx+1
            FV_Atm(1)%flagstruct%npx    = FV_Atm(1)%flagstruct%npx+1
         else
            FV_Atm(1)%flagstruct%ntiles = 1
            FV_Atm(1)%flagstruct%npy    = FV_Atm(1)%flagstruct%npy+1
            FV_Atm(1)%flagstruct%npx    = FV_Atm(1)%flagstruct%npx+1
         endif
      endif

      call MAPL_GetResource( MAPL,       &
                             ndt,        &
                             'RUN_DT:',  &
                             default=0,  &
                             RC=STATUS )
      VERIFY_(STATUS)
      DT = ndt

#ifdef ADJOINT
      call MAPL_GetResource( MAPL,              &
                             modelPhase,        &
                             'MODEL_PHASE:',    &
                             default='FORWARD', &
                             RC=STATUS )
      _VERIFY(STATUS)
      isAdjoint = .false.
      if (trim(ModelPhase) == 'ADJOINT') &
           isAdjoint = .true.
      if (isAdjoint) dt = -dt

      ! The advection adjoint is the unmodified forward scheme run with the
      ! reversed winds (reverse_wind_pure). The variant with face-edge and
      ! corner corrections (reverse_wind) was removed on 2026-10-14: it was
      ! equal inside faces and worse at seams, corners and the polar crossing
      ! (tests_final P2a, P2c, T4b). ADV_ADJ_METHOD is still read so that run
      ! directories asking for the removed variant stop here.
      call MAPL_GetResource( MAPL,                         &
                             advAdjMethod,                 &
                             'ADV_ADJ_METHOD:',            &
                             default='reverse_wind_pure',  &
                             RC=STATUS )
      _VERIFY(STATUS)
      ! discrete_linear (added 2026-10-19): the exact transpose of the forward step, TAF-generated
      ! from fv_tracer2d / tp_core / fv_mapz (fvdycore/model/adjoint). Exact for a linear
      ! forward scheme (hord_tr = 2, kord_tr = 17, as in run_P2lin.pbs; the forward run must use
      ! the same scheme). It acts on the adjoint directly: no air-mass weighting, no rescaling.
      if ( trim(advAdjMethod) /= 'reverse_wind_pure' .and. &
           trim(advAdjMethod) /= 'discrete_linear' ) then
         _ASSERT(.false., 'ADV_ADJ_METHOD must be reverse_wind_pure or discrete_linear')
      endif
      advAdjDiscrete = ( trim(advAdjMethod) == 'discrete_linear' )
      ! mass fixer in the discrete adjoint:
      !   none       : factor 1 (the forward's s(q) is not differentiated)
      !   linearized : adjoint of q*s(q), with s(q) linearized about the forward q of the step.
      !                That q is the forward partner of the adjoint tracer (SPC_X for SPC_X_ADJ) in
      !                TRADV. gcadj_internal_checkpoint.<T> holds the field AFTER the forward
      !                advection of step T -> T+dt (it is written at the start of GEOS-Chem), and
      !                the column adjoints leave q(T+dt) in the internal state. The forward run
      !                therefore exports the field BEFORE advection as TRAJ_PREADV; GEOS-Chem
      !                imports it (import checkpoint) and, in the adjoint run, copies it into the
      !                forward partner at the end of its step, so this q is q(T). Only for the
      !                tracer ADV_ADJ_TRAJ_TRACER (default SPC_CO2). The limiter decisions of the
      !                discrete adjoint use the same q.
      ! linearized (default): adjoint of the mass fixer linearized about the forward CO2 (DL is
      ! then the adjoint of the model); none: fixer factor held at 1 (development only)
      call MAPL_GetResource( MAPL, advAdjDiscFixer, 'ADV_ADJ_DISCRETE_FIXER:', &
                             default='linearized', RC=STATUS )
      _VERIFY(STATUS)
      if ( advAdjDiscrete .and. trim(advAdjDiscFixer) /= 'none' .and. &
           trim(advAdjDiscFixer) /= 'linearized' ) then
         _ASSERT(.false., 'ADV_ADJ_DISCRETE_FIXER must be none or linearized')
      endif
      ! 1: every step, also run the forward step of a test field with offline_tracer_advection
      !    and with the discrete code on the same inputs, and print the difference
      call MAPL_GetResource( MAPL, advAdjDiscCheck, 'ADV_ADJ_DISCRETE_CHECK:', &
                             default=0, RC=STATUS )
      _VERIFY(STATUS)
      ! adjoint of the tracer halo update in the discrete adjoint: own (exact transpose of the
      ! forward exchange, any layout) or fms (mpp_update_domains_ad: wrong for non-square
      ! subdomains, e.g. 48 cores at C24; kept for comparison). CHECK 3: dot-product test of both
      ! at the first adjoint step.
      call MAPL_GetResource( MAPL, advAdjHaloAd, 'ADV_ADJ_HALO_AD:', default='own', RC=STATUS )
      _VERIFY(STATUS)
      ! forward tracer exported as TRAJ_PREADV (forward run) and used as the trajectory of the
      ! linearized discrete adjoint (adjoint run, restored by GEOS-Chem from the checkpoint)
      call MAPL_GetResource( MAPL, advAdjTrajTracer, 'ADV_ADJ_TRAJ_TRACER:', default='SPC_CO2', &
                             RC=STATUS )
      _VERIFY(STATUS)
      if ( trim(advAdjHaloAd) /= 'own' .and. trim(advAdjHaloAd) /= 'fms' ) then
         _ASSERT(.false., 'ADV_ADJ_HALO_AD must be own or fms')
      endif
      if (MAPL_Am_I_Root()) write(*,*) 'ADVCORE adv_adj_method=', trim(advAdjMethod), &
                                       ' discrete fixer=', trim(advAdjDiscFixer), ' check=', advAdjDiscCheck, &
                                       ' halo adjoint=', trim(advAdjHaloAd)

      ! Weight that turns the adjoint variable into a mixing-ratio-like
      ! quantity around the reversed-wind advection (GCHP.rc):
      !   airmass (default): dry air mass at t+dt before and at t after the advection
      !                    (all validated RW results; the adjoint gains or loses
      !                    div(mf)/dp of the cell every step)
      !   dpadv          : divide by the dry air mass at t+dt (pressure thickness x
      !                    area); multiply by the air mass this reversed step leaves
      !                    in the cell before the vertical remapping (PLEAdv), plus
      !                    (mass at t - mass at t+dt) x the divided adjoint. This is the
      !                    transpose of the forward update q' = (q dp1 + div F) / dpA,
      !                    which uses the post-advection thickness dpA = dp1 + div(mf),
      !                    not the met thickness (development option, not validated)
      ! (The air density is used only with Use_Total_Air_Pressure > 0; the former option
      ! airden, density on both sides with dry pressure, is removed.)
      call MAPL_GetResource( MAPL, advAdjWeight, 'ADV_ADJ_WEIGHT:', &
                             default='airmass', RC=STATUS )
      _VERIFY(STATUS)
      if ( trim(advAdjWeight) /= 'dpadv'   .and. &
           trim(advAdjWeight) /= 'airmass' ) then
         _ASSERT(.false., 'ADV_ADJ_WEIGHT must be airmass (default) or dpadv')
      endif
      ! 1 (default): restore the global sum of the adjoint after the advection
      call MAPL_GetResource( MAPL, advAdjRescale, 'ADV_ADJ_RESCALE:', &
                             default=1, RC=STATUS )
      _VERIFY(STATUS)
      if (MAPL_Am_I_Root()) write(*,*) 'ADVCORE adv_adj_weight=', trim(advAdjWeight), &
                                       ' adv_adj_rescale=', advAdjRescale

      if (MAPL_Am_I_Root()) then
          if (ADVCORE_ADJ_DEBUG) then
           write(*,*) 'ADVCORE_SETSERVICES_DT ndt=', ndt,               &
              ' model_phase=', trim(modelPhase),               &
              ' isAdjoint=', isAdjoint, ' final_dt=', dt
          endif
      endif
#else
      if (MAPL_Am_I_Root()) then
          if (ADVCORE_ADJ_DEBUG) then
           write(*,*) 'ADVCORE_SETSERVICES_DT ndt=', ndt,               &
              ' model_phase=FORWARD (no ADJOINT build)',       &
              ' isAdjoint=.false. final_dt=', dt
          endif
      endif
#endif

      ! Start up FV if AdvCore is running without FV3_DynCoreIsRunning
      !--------------------------------------------------
      if (.NOT. FV3_DynCoreIsRunning) then
         call fv_init2(FV_Atm, dt, grids_on_my_pe, p_split)
      end if

      ! Ending with a Generic SetServices call is a MAPL requirement 
      !-------------------------------------------------------------
      call MAPL_GenericSetServices    ( GC, rc=STATUS)
      VERIFY_(STATUS)

      RETURN_(ESMF_SUCCESS)

      end subroutine SetServices

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialize - initialization routine
!
! !INTERFACE:
!
  subroutine Initialize(GC, IMPORT, EXPORT, CLOCK, RC)
!
! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
      type(ESMF_State),    intent(inout) :: IMPORT ! Import state
      type(ESMF_State),    intent(inout) :: EXPORT ! Export state
      type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
!
! !OUTPUT PARAMETERS:
      integer, optional,   intent(  out) :: RC     ! Error code
!
! !DESCRIPTION:
!     This initialization routine creates the import and export states,
!     as well as the internal state, which is attached to the component.
!     It also determines the distribution (and therefore the grid) 
!     and performs allocations of persistent data, 
!
!EOP
!=============================================================================
!BOC

      character(len=ESMF_MAXSTR)         :: IAm
      integer                            :: STATUS
      character(len=ESMF_MAXSTR)         :: COMP_NAME
      type(ESMF_Config)                  :: CF
      type (MAPL_MetaComp),      pointer :: MAPL
      type (ESMF_VM)                     :: VM
      real(REAL8), pointer               :: temp2d_r8(:,:) => NULL()
      integer                            :: IS, IE, JS, JE
      logical                            :: gridCreated
      type(ESMF_Grid)                    :: grid

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

      Iam = "Initialize"
      call ESMF_GridCompGet ( GC, name=COMP_NAME, config=CF, vm=VM, RC=STATUS )
      VERIFY_(STATUS)
      Iam = trim(COMP_NAME) // trim(Iam)

      ! Retrieve the pointer to the state
      ! ---------------------------------
      call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
      VERIFY_(STATUS)

      call MAPL_TimerOn(MAPL,"TOTAL")
      call MAPL_TimerOn(MAPL,"INITIALIZE")

      gridCreated=.false.
      call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
      VERIFY_(STATUS)
      call ESMF_GridCompGet(GC,grid=grid,rc=status)
      if (status == ESMF_SUCCESS) then
         call ESMF_GridValidate(grid,rc=status)
         if (status==ESMF_SUCCESS) GridCreated = .true.
      end if

      if (.not.GridCreated) then
         call MAPL_GridCreate(GC, rc=status)
         VERIFY_(STATUS)
      endif

      call MAPL_GenericInitialize(GC, IMPORT, EXPORT, CLOCK, RC=STATUS)
      VERIFY_(STATUS)
      ! Compute Grid-Cell Area
      ! ----------------------
      if (.NOT. FV3_DynCoreIsRunning) then
         IS = FV_Atm(1)%bd%isc
         IE = FV_Atm(1)%bd%iec
         JS = FV_Atm(1)%bd%jsc
         JE = FV_Atm(1)%bd%jec
         call MAPL_GetPointer(EXPORT, temp2d_r8, 'AREA', ALLOC=.TRUE., rc=status)
         VERIFY_(STATUS)
         temp2d_r8 = FV_Atm(1)%gridstruct%area_64(IS:IE,JS:JE)
         temp2d_r8 => NULL()
      endif

      call MAPL_TimerOff(MAPL,"INITIALIZE")
      call MAPL_TimerOff(MAPL,"TOTAL")

      RETURN_(ESMF_SUCCESS)

      end subroutine Initialize
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Run - run routine
!
! !INTERFACE:
!
      subroutine Run(GC, IMPORT, EXPORT, CLOCK, RC)
!
! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
      type(ESMF_State),    intent(inout) :: IMPORT ! Import state
      type(ESMF_State),    intent(inout) :: EXPORT ! Export state
      type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
!
! !OUTPUT PARAMETERS:
      integer, optional,   intent(  out) :: RC     ! Error code
!
! !DESCRIPTION:
! 
! The Run method advanced the advection one long time step, as
! specified in the configuration.  This may be broken down int a
! number of internal, small steps, also configurable.
!
!EOP
!=============================================================================
!BOC
! !LOCAL VARIABLES:
      character(len=ESMF_MAXSTR)    :: IAm
      integer                       :: STATUS
      character(len=ESMF_MAXSTR)    :: COMP_NAME
      type (ESMF_Grid)              :: ESMFGRID
      type (MAPL_MetaComp), pointer :: MAPL
      type (ESMF_Alarm)             :: ALARM

! Imports
      REAL(REAL8), POINTER, DIMENSION(:,:,:)   :: iDELPDRY
      REAL(REAL8), POINTER, DIMENSION(:,:,:)   :: iCX
      REAL(REAL8), POINTER, DIMENSION(:,:,:)   :: iCY
      REAL(REAL8), POINTER, DIMENSION(:,:,:)   :: iMFX
      REAL(REAL8), POINTER, DIMENSION(:,:,:)   :: iMFY
      REAL(REAL8), POINTER, DIMENSION(:,:,:)   :: iPLE0
      REAL(REAL8), POINTER, DIMENSION(:,:,:)   :: iPLE1
      REAL(REAL8), POINTER, DIMENSION(:,:,:)   :: iDryPLE0 ! GCHP dry
      REAL(REAL8), POINTER, DIMENSION(:,:,:)   :: iDryPLE1 ! GCHP dry
      REAL(REAL8), POINTER, DIMENSION(:,:,:)   :: iSPHU0   ! GCHP total
#ifdef ADJOINT
      REAL(REAL8), POINTER, DIMENSION(:,:,:)   :: iAIRDEN  ! dry air density
#endif

! Exports
      REAL(REAL8), POINTER, DIMENSION(:,:,:)   :: ePLE       ! GCHP export to GEOS-Chem
      REAL(REAL8), POINTER, DIMENSION(:,:,:)   :: eDryPLE    ! GCHP dry export to GEOS-Chem
      REAL(REAL8), POINTER, DIMENSION(:,:,:)   :: ePLEadv    ! GCHP post-advection PLE
#ifdef ADJOINT
      REAL(REAL8), POINTER, DIMENSION(:,:,:)   :: eTRAJ => NULL() ! tracer before advection
#endif

      ! GCHP R4 exports as temporary work-around for MAPL 2.55 bug
      REAL(REAL4), POINTER, DIMENSION(:,:,:)   :: ePLE_R4
      REAL(REAL4), POINTER, DIMENSION(:,:,:)   :: eDryPLE_R4
      REAL(REAL4), POINTER, DIMENSION(:,:,:)   :: ePLEadv_R4

! Locals
      REAL(FVPRC), POINTER, DIMENSION(:,:,:)   :: CX
      REAL(FVPRC), POINTER, DIMENSION(:,:,:)   :: CY
      REAL(FVPRC), POINTER, DIMENSION(:,:,:)   :: MFX
      REAL(FVPRC), POINTER, DIMENSION(:,:,:)   :: MFY
      REAL(FVPRC), POINTER, DIMENSION(:,:,:)   :: PLE0
      REAL(FVPRC), POINTER, DIMENSION(:,:,:)   :: PLE1
      REAL(FVPRC), POINTER, DIMENSION(:,:,:)   :: DryPLE0 ! GCHP dry
      REAL(FVPRC), POINTER, DIMENSION(:,:,:)   :: DryPLE1 ! GCHP dry
      REAL(FVPRC), POINTER, DIMENSION(:,:,:)   :: PLEAdv  ! GCHP total
      REAL(FVPRC), POINTER, DIMENSION(:,:,:)   :: SPHU0   ! GCHP total
      REAL(FVPRC), POINTER, DIMENSION(:,:,:)   :: DELPDRY ! import
#ifdef ADJOINT
      REAL(FVPRC), POINTER, DIMENSION(:,:,:)   :: AIRDEN => NULL() ! selected air density
   logical, allocatable                      :: isAdjointTracer(:)
      ! discrete_linear: inputs of the forward step, adjoint of one tracer, check fields
      REAL(FVPRC), ALLOCATABLE :: adjMfxF(:,:,:), adjMfyF(:,:,:), adjCxF(:,:,:), adjCyF(:,:,:)
      REAL(FVPRC), ALLOCATABLE :: adjQtraj(:,:,:)
      integer                  :: NF, lenAdj
      REAL(FVPRC), ALLOCATABLE :: adjPle0F(:,:,:), adjPle1F(:,:,:), adjLam(:,:,:)
      REAL(FVPRC), ALLOCATABLE :: adjQt4(:,:,:,:), adjQt(:,:,:), adjPleT(:,:,:), adjPleAdvT(:,:,:)
      REAL(REAL8)              :: adjChkLoc, adjChkGlob, adjChkRef
   logical                                   :: isAdjField
   integer                                   :: nAdjointTracers
#endif
      REAL(FVPRC), POINTER, DIMENSION(:)       :: AK
      REAL(FVPRC), POINTER, DIMENSION(:)       :: BK
      REAL(REAL8), allocatable :: ak_r8(:),bk_r8(:)
      REAL(FVPRC), POINTER, DIMENSION(:,:,:,:) :: TRACERS
      REAL(FVPRC) :: MASS1, TMASS1(ntracers)
      REAL(FVPRC) :: spc1_vv, spc1_mass_dp_met, spc1_mass_dp_rst ! debug
      REAL(FVPRC) :: dp_met, dp_rst                              ! debug
      TYPE(AdvCoreTracers), POINTER :: advTracers(:)
      type(ESMF_FieldBundle) :: TRADV
      type(ESMF_Field)       :: field
      type(ESMF_Array)       :: array
      INTEGER :: IM, JM, LM, N, NQ, LS, NADV
      REAL(FVPRC) :: PTOP, PINT
      REAL(REAL8) :: ptop_r8,pint_r8
! Temporaries for exports/tracers
      REAL, POINTER :: temp3D(:,:,:)
      real(REAL4),        pointer     :: tracer_r4 (:,:,:)
      real(REAL8),        pointer     :: tracer_r8 (:,:,:)
      character(len=ESMF_MAXSTR)    :: fieldName
      type(ESMF_TypeKind_Flag)      :: kind
      character(len=ESMF_MAXSTR)    :: STRING
! Excluding tracers
      type(ESMF_FieldBundle), save        :: bundleAdv
      type (ESMF_Config)                  :: CF
      logical                             :: adjustTracers
      type(ESMF_Alarm)                    :: predictorAlarm
      type(ESMF_Grid)                     :: bgrid
      integer, save                       :: nq_saved = 0
      integer                             :: i,j
      integer                             :: nqt
      logical                             :: tend
      logical                             :: exclude
      real(FVPRC)                         :: mfxMin, mfxMax, mfyMin, mfyMax
      real(FVPRC)                         :: cxMin, cxMax, cyMin, cyMax
      character(len=ESMF_MAXSTR)          :: tmpstring
      character(len=ESMF_MAXSTR)          :: adjustTracerMode
      character(len=ESMF_MAXSTR), allocatable :: xlist(:)
      character(len=ESMF_MAXSTR), allocatable :: biggerlist(:)
      integer, parameter                  :: XLIST_MAX = 60

      logical :: advOn   ! advect in this call?
#ifdef ADJOINT
!     reverse time debug info
      integer, parameter                :: DI = 3, DJ = 4, DL = 5
      ! Debug variables
      INTEGER, parameter             :: I_DBG = 6, J_DBG = 5, L_DBG=1
   real(FVPRC)                    :: rhoDryMin, rhoDryMax
   real(FVPRC)                    :: rhoUseMin, rhoUseMax
   real(REAL8)                    :: adjLocalSum, adjGlobalSum
   real(REAL8), allocatable       :: adjPreAdvSums(:)
   logical                        :: adjUseAirMass
   logical                        :: adjUseDpAdv
   real(FVPRC), allocatable       :: adjMuIn(:,:,:,:)  ! adjoint / air mass before the advection
   type(ESMF_Time)                :: adjCurrTime, adjStopTime
   integer                        :: L
   integer                        :: nAdjScaled
   type(ESMF_VM)                  :: vmRun
#endif

! Get my name and set-up traceback handle
! ---------------------------------------

      Iam = 'Run'
      call ESMF_GridCompGet( GC, name=COMP_NAME, CONFIG=CF, grid=ESMFGRID, &
                             RC=STATUS )
      VERIFY_(STATUS)
      Iam = trim(COMP_NAME) // Iam

!WMP  if (AdvCore_Advection>0) then

! Get parameters from generic state.
!-----------------------------------
      call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_Get( MAPL, IM=IM, JM=JM, LM=LM,   &
                                RUNALARM = ALARM, &
                                      RC = STATUS )
      VERIFY_(STATUS)

      call MAPL_TimerOn(MAPL,"TOTAL")
      call MAPL_TimerOn(MAPL,"RUN")

! Get AKs and BKs for vertical grid
!----------------------------------
      AllOCATE( AK(LM+1) ,stat=STATUS )
      VERIFY_(STATUS)
      AllOCATE( BK(LM+1) ,stat=STATUS )
      VERIFY_(STATUS)
      AllOCATE( AK_r8(LM+1) ,stat=STATUS )
      VERIFY_(STATUS)
      AllOCATE( BK_r8(LM+1) ,stat=STATUS )
      VERIFY_(STATUS)
      call set_eta(LM,LS,ptop_r8,pint_r8,ak_r8,bk_r8)
      ptop=ptop_r8
      pint=pint_r8
      ak=ak_r8
      bk=bk_r8

      CALL MAPL_GetPointer(IMPORT, iPLE0, 'PLE0', ALLOC = .TRUE., RC=STATUS)
      VERIFY_(STATUS)
      CALL MAPL_GetPointer(IMPORT, iPLE1, 'PLE1', ALLOC = .TRUE., RC=STATUS)
      VERIFY_(STATUS)
      CALL MAPL_GetPointer(IMPORT, iMFX,   'MFX', ALLOC = .TRUE., RC=STATUS)
      VERIFY_(STATUS)
      CALL MAPL_GetPointer(IMPORT, iMFY,   'MFY', ALLOC = .TRUE., RC=STATUS)
      VERIFY_(STATUS)
      CALL MAPL_GetPointer(IMPORT, iCX,     'CX', ALLOC = .TRUE., RC=STATUS)
      VERIFY_(STATUS)
      CALL MAPL_GetPointer(IMPORT, iCY,     'CY', ALLOC = .TRUE., RC=STATUS)
      VERIFY_(STATUS)

      ! Using dry versus total air in GCHP
      IF ( Use_Total_Air_Pressure > 0 ) THEN
         CALL MAPL_GetPointer(IMPORT,iSPHU0,'SPHU0',ALLOC = .TRUE.,RC=STATUS)
         VERIFY_(STATUS)
         ALLOCATE( SPHU0(IM,JM,LM  ) )
         SPHU0   = iSPHU0
      ELSE
         CALL MAPL_GetPointer(IMPORT,iDryPLE0,'DryPLE0',ALLOC=.TRUE.,RC=STATUS)
         VERIFY_(STATUS)
         CALL MAPL_GetPointer(IMPORT,iDryPLE1,'DryPLE1',ALLOC=.TRUE.,RC=STATUS)
         VERIFY_(STATUS)
         ALLOCATE( DryPLE0(IM,JM,LM+1) )
         ALLOCATE( DryPLE1(IM,JM,LM+1) )
#ifdef ADJOINT
         if (isAdjoint) then
            DryPLE0 = iDryPLE1
            DryPLE1 = iDryPLE0
         else
#endif
            DryPLE0 = iDryPLE0
            DryPLE1 = iDryPLE1
#ifdef ADJOINT
         end if
#endif
      ENDIF

#ifdef ADJOINT
      if (isAdjoint) then
         CALL MAPL_GetPointer(IMPORT, iAIRDEN, 'ADV_AIRDEN', ALLOC=.TRUE., RC=STATUS)
         VERIFY_(STATUS)
         ALLOCATE( AIRDEN(IM,JM,LM) )
         AIRDEN = iAIRDEN

         rhoDryMin = minval(AIRDEN)
         rhoDryMax = maxval(AIRDEN)

         ! For total-air runs, convert dry-air density to moist-air density
         if ( Use_Total_Air_Pressure > 0 ) then
            where (1.0_FVPRC - SPHU0 > 0.0_FVPRC)
               AIRDEN = AIRDEN / (1.0_FVPRC - SPHU0)
            end where
         endif

         if (MAPL_Am_I_Root() .and. firstRun) then
            if (ADVCORE_ADJ_DEBUG) then
               rhoUseMin = minval(AIRDEN)
               rhoUseMax = maxval(AIRDEN)
               if (Use_Total_Air_Pressure > 0) then
                  write(*,*) 'ADVCORE_DENSITY first adjoint step: ',                 &
                             'dry_min=', rhoDryMin, ' dry_max=', rhoDryMax,         &
                             ' moist_used_min=', rhoUseMin, ' moist_used_max=', rhoUseMax
               else
                  write(*,*) 'ADVCORE_DENSITY first adjoint step: ',                 &
                             'dry_min=', rhoDryMin, ' dry_max=', rhoDryMax
               endif
            endif
         endif
      endif
#endif

      ! Advect in this call?  The adjoint calls at the end of the window (clock
      ! at the stop time) do not correspond to any forward step: their winds are
      ! those of the step after the window.  With physics on the adjoint is still
      ! zero there; in transport-only runs a built-in seed is already set, so
      ! advecting would add a step that the forward run does not have.
      advOn = (AdvCore_Advection > 0)
#ifdef ADJOINT
      if (isAdjoint) then
         call ESMF_ClockGet(CLOCK, currTime=adjCurrTime, stopTime=adjStopTime, rc=STATUS)
         VERIFY_(STATUS)
         if (adjCurrTime == adjStopTime) advOn = .false.
      endif
#endif

      ALLOCATE(  PLE0(IM,JM,LM+1) )
      ALLOCATE(  PLE1(IM,JM,LM+1) )
      ALLOCATE(PLEAdv(IM,JM,LM+1) )
      ALLOCATE(   MFX(IM,JM,LM  ) )
      ALLOCATE(   MFY(IM,JM,LM  ) )
      ALLOCATE(    CX(IM,JM,LM  ) )
      ALLOCATE(    CY(IM,JM,LM  ) )

      PLEAdv = 0.0d0
#ifdef ADJOINT
      if (isAdjoint) then
         ! Backward integration: swap pressure-edge fields
         ! MFX/MFY/CX/CY are already negated upstream by GCHPctmEnv (dt = -dt)
         PLE0 = iPLE1
         PLE1 = iPLE0
      else
#endif
         PLE0 = iPLE0
         PLE1 = iPLE1
#ifdef ADJOINT
      end if
#endif
      MFX = iMFX
      MFY = iMFY
      CX  = iCX
      CY  = iCY

      ! If first run, get import DELPDRY (restart field from GEOS-Chem internal state)
      if ( firstRun ) THEN
         call MAPL_GetPointer(IMPORT, iDELPDRY, 'DELPDRY', NotFoundOK=.TRUE., RC=STATUS)
         VERIFY_(STATUS)
      endif

      ! The quantities to be advected come as friendlies in a bundle
      !  in the import state.
      !--------------------------------------------------------------

      call ESMF_StateGet(IMPORT, "TRADV", TRADV, rc=STATUS)
      VERIFY_(STATUS)

      ! ALT: this section attempts to limit the amount of advected tracers
      !-------------------------------------------------------------------
      adjustTracers = .false.
      call MAPL_GetResource ( MAPL,                         &
                              adjustTracerMode,             &
                              'EXCLUDE_ADVECTION_TRACERS:', &
                              default='ALWAYS',             &
                              rc=status )
      VERIFY_(STATUS)
      if (adjustTracerMode == 'ALWAYS') then
         adjustTracers = .true.
      else if (adjustTracerMode == 'PREDICTOR') then
         !get PredictorAlarm from clock
         call ESMF_ClockGetAlarm(clock, alarmName='PredictorAlarm', &
              alarm=PredictorAlarm, rc=status)
         if (status == ESMF_SUCCESS) then
            !check if ringing
            if (ESMF_AlarmIsRinging(predictorAlarm)) then
               adjustTracers = .true.
            end if
         end if
      else if (adjustTracerMode == 'NO') then
         ! Proceed without warning
         adjustTracers = .false.
      else
         !call WRITE_PARALLEL('Invalid option, ignored')
         adjustTracers = .false.
      end if
      if (adjustTracers) then
         if (firstRun) then
            ! get the list of excluded tracers from resource
            n = 0
            call ESMF_ConfigFindLabel(CF,'EXCLUDE_ADVECTION_TRACERS_LIST:',&
                                      rc=STATUS )
            if(STATUS==ESMF_SUCCESS) then

               tend  = .false.
               allocate(xlist(XLIST_MAX), stat=status)
               VERIFY_(STATUS)
               do while (.not.tend)
                  call ESMF_ConfigGetAttribute (CF,value=tmpstring,&
                                                default='',rc=STATUS)
                                        !ALT: we don't check return status!!!
                  if (tmpstring /= '')  then
                     n = n + 1
                     if (n > size(xlist)) then
                        allocate( biggerlist(2*n), stat=status )
                        VERIFY_(STATUS)
                        biggerlist(1:n-1)=xlist
                        call move_alloc(from=biggerlist, to=xlist)
                     end if
                     xlist(n) = tmpstring
                  end if
                  call ESMF_ConfigNextLine(CF,tableEnd=tend,rc=STATUS )
                  VERIFY_(STATUS)
               enddo
            end if

            ! Count the number of tracers
            !---------------------
            call ESMF_FieldBundleGet(TRADV, grid=bgrid,fieldCount=nqt,  RC=STATUS)
            VERIFY_(STATUS)
            BundleAdv = ESMF_FieldBundleCreate ( name='xTRADV', rc=STATUS )
            VERIFY_(STATUS)
            call ESMF_FieldBundleSet ( BundleAdv, grid=bgrid, rc=STATUS )
            VERIFY_(STATUS)
            !loop over NQ in TRADV
            do i = 1, nqt
               !get field from TRADV and its name
               call ESMF_FieldBundleGet(TRADV, fieldIndex=i, field=field, &
                                        rc=status)
               VERIFY_(STATUS)
               call ESMF_FieldGet(FIELD, name=fieldname, RC=STATUS)
               VERIFY_(STATUS)
               !exclude everything that is not cloud/water species
               if ( (FV3_DynCoreIsRunning      ) .and. &
                   ( (TRIM(fieldname) == 'Q'       ) .or. &
                     (TRIM(fieldname) == 'QLCN'    ) .or. &
                     (TRIM(fieldname) == 'QLLS'    ) .or. &
                     (TRIM(fieldname) == 'QICN'    ) .or. &
                     (TRIM(fieldname) == 'QILS'    ) .or. &
                     (TRIM(fieldname) == 'CLCN'    ) .or. &
                     (TRIM(fieldname) == 'CLLS'    ) .or. &
                     (TRIM(fieldname) == 'NCPL'    ) .or. &
                     (TRIM(fieldname) == 'NCPI'    ) .or. &
                     (TRIM(fieldname) == 'QRAIN'   ) .or. &
                     (TRIM(fieldname) == 'QSNOW'   ) .or. &
                     (TRIM(fieldname) == 'QGRAUPEL') ) ) then
                   ! write(STRING,'(A,A)') "ADV is excluding ", TRIM(fieldname)
                   ! call WRITE_PARALLEL( trim(STRING)   )
                     n = n + 1
                     if (n > size(xlist)) then
                        allocate( biggerlist(2*n), stat=status )
                        VERIFY_(STATUS)
                        biggerlist(1:n-1)=xlist
                        call move_alloc(from=biggerlist, to=xlist)
                     end if
                     xlist(n) = TRIM(fieldname)
               end if
               !loop over exclude_list
               exclude = .false.
               do j = 1, n
                  if (fieldname == xlist(j)) then
                     exclude = .true.
                     exit
                  end if
               end do
               if (.not. exclude) then
                  call MAPL_FieldBundleAdd(BundleAdv, FIELD, RC=STATUS)
                  VERIFY_(STATUS)
               end if
            end do

            if (allocated(xlist)) then
           !   ! Just in case xlist was allocated, but nothing was in it,
           !   ! could have garbage
           !   if (n > 0) then
           !      call ESMF_FieldBundleRemove(TRADV, fieldNameList=xlist, &
           !         relaxedFlag=.true., rc=status)
           !      VERIFY_(STATUS)
           !   end if
               deallocate(xlist)
            end if

         end if ! firstRun
         TRADV = bundleAdv
      end if ! adjustTracers

      call ESMF_FieldBundleGet(TRADV, fieldCount=NQ,    rc=STATUS)
      VERIFY_(STATUS)

      ! GCHP: If using total air then add one tracer for moisture for use in
      ! post-advection MMR conversion from total to dry air
      IF ( Use_Total_Air_Pressure > 0 ) THEN
         NAdv = nq + 1
      ELSE
         NAdv = nq
      ENDIF

      if (NQ > 0) then

         ! Allocate list of tracers big enough to hold all items in the bundle
         !--------------------------------------------------------------------
         ALLOCATE( TRACERS(IM,JM,LM,NAdv),stat=STATUS ) ! Includes SPHU tracer
         VERIFY_(STATUS)
         ALLOCATE( advTracers(NQ),stat=STATUS ) ! Does not include SPHU tracer
         VERIFY_(STATUS)
#ifdef ADJOINT
         ALLOCATE( isAdjointTracer(NAdv), stat=STATUS )
         VERIFY_(STATUS)
         isAdjointTracer = .false.
         nAdjointTracers = 0
#endif

         if (NQ /= NQ_SAVED) then
            write(STRING,'(A,I5,A)') "AdvCore is Advecting the following ", nq, " tracers in FV3:"
            call WRITE_PARALLEL( trim(STRING)   )
         end if

         ! Go through the bundle copying the friendlies into the tracer list.
         !-------------------------------------------------------------------
         do N=1,NQ
            call ESMF_FieldBundleGet(TRADV, fieldIndex=N, field=FIELD,RC=STATUS)
            VERIFY_(STATUS)
            call ESMF_FieldGet(field, array=array, name=fieldName, RC=STATUS)
            VERIFY_(STATUS)
            call ESMF_ArrayGet(array,typekind=kind, rc=status )
            VERIFY_(STATUS)
            advTracers(N)%is_r4 = (kind == ESMF_TYPEKIND_R4)   ! Is real*4?
            advTracers(N)%tName = fieldName



#ifdef ADJOINT
            isAdjField = .false.
            if (len_trim(fieldName) >= 4) then
               isAdjField = (fieldName(len_trim(fieldName)-3:len_trim(fieldName)) == '_ADJ')
            endif
            isAdjointTracer(N) = isAdjField
            if (isAdjField) then
               nAdjointTracers = nAdjointTracers + 1
               if (isAdjoint .and. firstRun .and. MAPL_Am_I_Root()) then
                  if (ADVCORE_ADJ_DEBUG) then
                     write(*,*) 'ADVCORE_ADJ_TRACER_INDEX idx=', N, ' name=', trim(fieldName)
                  endif
               endif
            endif
#endif

            if (NQ /= NQ_SAVED) then
               call WRITE_PARALLEL( trim(fieldName) )
            endif

            if (advTracers(N)%is_r4) then
               call ESMF_ArrayGet(array,farrayptr=tracer_r4, rc=status )
               VERIFY_(STATUS)
               advTracers(N)%content_r4 => tracer_r4
               TRACERS(:,:,:,N) = advTracers(N)%content_r4
            else
               call ESMF_ArrayGet(array,farrayptr=tracer_r8, rc=status )
               VERIFY_(STATUS)
               advTracers(N)%content => tracer_r8
               TRACERS(:,:,:,N) = advTracers(N)%content
            end if

         end do

#ifdef ADJOINT
         if (isAdjoint .and. firstRun .and. MAPL_Am_I_Root()) then
            if (ADVCORE_ADJ_DEBUG) then
               write(*,*) 'ADVCORE_ADJ_TRACER_COUNT n_adj=', nAdjointTracers, ' n_total=', NQ
            endif
         endif

         if (isAdjoint .and. ADVCORE_ADJ_DEBUG) then
            call ESMF_VMGetCurrent(vmRun, rc=STATUS)
            VERIFY_(STATUS)
            do N=1,NQ
               if (.not. isAdjointTracer(N)) cycle
               adjLocalSum = SUM( REAL(TRACERS(:,:,:,N), REAL8) )
               call MAPL_CommsAllReduceSum(vmRun, sendbuf=adjLocalSum, recvbuf=adjGlobalSum, &
                                           cnt=1, rc=STATUS)
               VERIFY_(STATUS)
               if (MAPL_Am_I_Root()) then
                  write(*,*) 'ADVCORE_ADJ_GLOBAL_SUM stage=after_import idx=', N, &
                             ' name=', trim(advTracers(N)%tName), ' sum=', adjGlobalSum
               endif
            enddo
         endif
#endif

         !----------------------------------------------------------
         ! Optional pre-scaling debug : Sum area-weighted
         ! values and print (first run only)
         !----------------------------------------------------------
         !if (chk_mass .and. firstRun) then
         !   dP_rst=0.d0
         !   spc1_mass_dp_rst=0.d0
         !
         !   ! Delta pressure from the restart file
         !   if ( ASSOCIATED(iDELPDRY) ) &
         !        call global_integral_dp_rst(dp_rst, iDELPDRY, IM,JM,LM)
         !
         !   ! Delta pressure as computed from met
         !   call global_integral_dp(dp_met, DryPLE0, IM,JM,LM)
         !
         !   ! v/v for species 1
         !   call global_integral_vv(spc1_vv, Tracers(:,:,:,1), IM,JM,LM)
         !
         !   ! Species 1 v/v * dp (rst) (proxy for mass in restart file)
         !   if ( ASSOCIATED(iDELPDRY) ) &
         !        call global_integral_trmass_from_dp(spc1_mass_dp_rst, &
         !        Tracers(:,:,:,1), iDELPDRY, IM,JM,LM)
         !
         !   ! Species 1 v/v * dp (met) (dp computed in function)
         !   call global_integral_trmass_from_p(spc1_mass_dp_met, &
         !        Tracers(:,:,:,1), DryPLE0, IM,JM,LM)
         !
         !   if (is_master()) then
         !      write(6,*) "Global area-weighted sums: dp met, dp rst, spc1 vv, spc1 vv*dp_met, spc1 vv*dp_rst"
         !      write(6,100)  dp_met, dp_rst, spc1_vv, &
         !           spc1_mass_dp_met, spc1_mass_dp_rst
         !   endif
100       format('First run, before scaling: ',e21.14,' ',e21.14,' ',e21.14,' ',e21.14,' ',e21.14)
         !
         !endif ! chk_mass .and. firstRun

         ! If first timestep and delta pressure in the restart file is non-zero,
         ! then scale mixing ratios by ratio of restart file delta pressure to
         ! run-time met delta pressure in order to conserve restart file mass
#ifdef ADJOINT
         ! Not in adjoint runs: their first call is at the end of the window,
         ! where the forward run does not scale. The transpose of the forward
         ! scaling at the start of the window (multiply the adjoint by the same
         ! ratio after the last advection adjoint) is not implemented (TODO).
         if ( firstRun .and. associated(iDELPDRY) .and. .not. isAdjoint ) THEN
#else
         if ( firstRun .and. associated(iDELPDRY) ) THEN
#endif

            ! Only scale mixing ratios if non-zero delta pressures in the restart file
            if ( sum(iDELPDRY) > 0.d0 ) THEN
               ALLOCATE( DELPDRY(IM,JM,LM) )
               DELPDRY = iDELPDRY
               if (AdvCore_Advection>0) then
                  if (Use_Total_Air_Pressure > 0) then
                     call scale_tracers_by_pressure_ratio(tracers, PLE0, &
                          DELPDRY, IM, JM, LM, NAdv)
                  else
                     call scale_tracers_by_pressure_ratio(tracers, DryPLE0, &
                          DELPDRY, IM, JM, LM, NAdv)
                  endif
               else
                  if (Use_Total_Air_Pressure > 0) then
                     call scale_tracers_by_pressure_ratio(tracers, PLE1, &
                          DELPDRY, IM, JM, LM, NAdv)
                  else
                     call scale_tracers_by_pressure_ratio(tracers, DryPLE1, &
                          DELPDRY, IM, JM, LM, NAdv)
                  endif
               endif
               DEALLOCATE( DELPDRY )
            endif

         endif

#ifdef ADJOINT
         ! Forward run: export ADV_ADJ_TRAJ_TRACER before the advection of this step (after the
         ! start-up rescaling, in the units of the internal state) as TRAJ_PREADV. Not with
         ! Use_Total_Air_Pressure: the adjoint run would convert it with its own SPHU0.
         if ( .not. isAdjoint ) then
            call MAPL_GetPointer( EXPORT, eTRAJ, 'TRAJ_PREADV', NotFoundOK=.TRUE., RC=STATUS )
            VERIFY_(STATUS)
            if ( associated(eTRAJ) ) then
               eTRAJ = 0.0_REAL8
               do N=1,NQ
                  if ( trim(advTracers(N)%tName) == trim(advAdjTrajTracer) ) &
                     eTRAJ = REAL(TRACERS(:,:,:,N), REAL8)
               enddo
            endif
         endif
#endif

         ! If using total air then set extra tracer to specific humidity and
         ! convert all other tracers from kg/kg dry to kg/kg total air
         if ( Use_Total_Air_Pressure > 0 ) then
            tracers(:,:,:,nq+1) = sphu0(:,:,:)
            do N=1,NQ
               tracers(:,:,:,N) = tracers(:,:,:,N) * (1.0 - sphu0)
            enddo
          endif

         if (NQ /= NQ_SAVED) then
            NQ_SAVED = NQ
         end if

         ! Check Mass conservation
         if (chk_mass .and. firstRun) then

            ! Compute initial mass (proxy). Do this differently based on
            ! whether advection, on or off and whether using total or dry
            ! air pressure in advection. Mass of first run is saved for
            ! comparison in all subsequent timesteps.
            if (AdvCore_Advection>0) then
               if ( Use_Total_Air_Pressure > 0 ) then
                  MASS0 = g_sum( FV_Atm(1)%domain,            &
                       PLE0(:,:,LM),                &
                       is,                          &
                       ie,                          &
                       js,                          &
                       je,                          &
                       FV_Atm(1)%ng,                &
                       FV_Atm(1)%gridstruct%area_64,&
                       1,                           &
                       .true. )
                  call global_integral(TMASS0, TRACERS, PLE0, IM,JM,LM,NAdv)
               else
                  MASS0 = g_sum( FV_Atm(1)%domain,            &
                       DryPLE0(:,:,LM),             &
                       is,                          &
                       ie,                          &
                       js,                          &
                       je,                          &
                       FV_Atm(1)%ng,                &
                       FV_Atm(1)%gridstruct%area_64,&
                       1,                           &
                       .true. )
                  call global_integral(TMASS0, TRACERS, DryPLE0, IM,JM,LM,NAdv)
                  if (MASS0 /= 0.0) TMASS0=TMASS0/MASS0
               endif
            else
               if ( Use_Total_Air_Pressure > 0 ) then
                  MASS0 = g_sum( FV_Atm(1)%domain,            &
                       PLE1(:,:,LM),                &
                       is,                          &
                       ie,                          &
                       js,                          &
                       je,                          &
                       FV_Atm(1)%ng,                &
                       FV_Atm(1)%gridstruct%area_64,&
                       1,                           &
                       .true. )
                  call global_integral(TMASS0, TRACERS, PLE1, IM,JM,LM,NAdv)
               else
                  MASS0 = g_sum( FV_Atm(1)%domain,            &
                       DryPLE1(:,:,LM),             &
                       is,                          &
                       ie,                          &
                       js,                          &
                       je,                          &
                       FV_Atm(1)%ng,                &
                       FV_Atm(1)%gridstruct%area_64,&
                       1,                           &
                       .true.)
                  call global_integral(TMASS0, TRACERS, DryPLE1, IM,JM,LM,NQ)
               endif
               if (MASS0 /= 0.0) TMASS0=TMASS0/MASS0
            endif

            !----------------------------------------------------------
            ! Optional pre-advection debug : Sum additional
            ! area-weighted values and print (first run only)
            !----------------------------------------------------------
            !dP_rst=0.d0
            !spc1_mass_dp_rst=0.d0
            !
            !! Delta pressure from the restart file
            !if ( ASSOCIATED(iDELPDRY) ) &
            !     call global_integral_dp_rst(dp_rst, iDELPDRY, IM,JM,LM)
            !
            !! Delta pressure as computed from met
            !call global_integral_dp(dp_met, DryPLE0, IM,JM,LM)
            !
            !! v/v for species 1
            !call global_integral_vv(spc1_vv, Tracers(:,:,:,1), IM,JM,LM)
            !
            !! Species 1 v/v * dp (rst) (proxy for mass in restart file)
            !if ( ASSOCIATED(iDELPDRY) ) &
            !     call global_integral_trmass_from_dp(spc1_mass_dp_rst, &
            !     Tracers(:,:,:,1), iDELPDRY, IM,JM,LM)
            !
            !! Species 1 v/v * dp (met) (dp computed in function)
            !call global_integral_trmass_from_p(spc1_mass_dp_met, &
            !     Tracers(:,:,:,1), DryPLE0, IM,JM,LM)
            !
            !! Print the values
            !if (is_master()) then
            !   write(6,*) "Global area-weighted sums: sfc p met, dp met, dp rst, spc1 vv, spc1 vv*dp_met, spc1 vv*dp_rst"
            !   write(6,101)  MASS0, dp_met, dp_rst, spc1_vv, &
            !        spc1_mass_dp_met, spc1_mass_dp_rst
            !endif
101         format('First run, after scaling: ',e21.14,' ',e21.14,' ',e21.14,' ',e21.14,' ',e21.14,' ',e21.14)
            !----------------------------------------------------------

         endif ! chk_mass .and. firstRun

         ! Run FV3 advection
         !------------------
#ifdef ADJOINT
         if (advOn) then
         IF (MAPL_Am_I_Root()) THEN
            WRITE(*,546) dt
546         FORMAT(' calling offline_tracer_advection with timestep = ', f8.3)
         ENDIF
#else
         if (advOn) then
#endif

#ifdef ADJOINT
            if (isAdjoint .and. .not. advAdjDiscrete) then
               if (ADVCORE_ADJ_DEBUG) then
                  call ESMF_VMGetCurrent(vmRun, rc=STATUS)
                  VERIFY_(STATUS)
                  do N=1,NQ
                     if (.not. isAdjointTracer(N)) cycle
                     adjLocalSum = SUM( REAL(TRACERS(:,:,:,N), REAL8) )
                     call MAPL_CommsAllReduceSum(vmRun, sendbuf=adjLocalSum, recvbuf=adjGlobalSum, &
                                                 cnt=1, rc=STATUS)
                     VERIFY_(STATUS)
                     if (MAPL_Am_I_Root()) then
                        write(*,*) 'ADVCORE_ADJ_GLOBAL_SUM stage=before_div_airden idx=', N, &
                                   ' name=', trim(advTracers(N)%tName), ' sum=', adjGlobalSum
                     endif
                  enddo
               endif

               ! Save global sum before density scaling for conservation fix
               call ESMF_VMGetCurrent(vmRun, rc=STATUS)
               VERIFY_(STATUS)
               ALLOCATE( adjPreAdvSums(NQ), stat=STATUS )
               VERIFY_(STATUS)
               adjPreAdvSums = 0.0_REAL8
               do N=1,NQ
                  if (.not. isAdjointTracer(N)) cycle
                  adjLocalSum = SUM( REAL(TRACERS(:,:,:,N), REAL8) )
                  call MAPL_CommsAllReduceSum(vmRun, sendbuf=adjLocalSum, recvbuf=adjGlobalSum, &
                                             cnt=1, rc=STATUS)
                  VERIFY_(STATUS)
                  adjPreAdvSums(N) = adjGlobalSum
               enddo

               ! lambda / (air weight) is advected like a mixing ratio. In adjoint
               ! mode DryPLE0 holds the pressure edges at t+dt (start of the
               ! backward step) and DryPLE1 those at t.
               adjUseAirMass = ( ( trim(advAdjWeight) == 'airmass' .or. &
                                   trim(advAdjWeight) == 'dpadv' ) .and. &
                                 Use_Total_Air_Pressure < 1 )
               adjUseDpAdv   = ( trim(advAdjWeight) == 'dpadv' .and. &
                                 Use_Total_Air_Pressure < 1 )
               do N=1,NQ
                  if (.not. isAdjointTracer(N)) cycle
                  if (adjUseAirMass) then
                     do L=1,LM
                        TRACERS(:,:,L,N) = TRACERS(:,:,L,N) / &
                             ( (DryPLE0(:,:,L+1) - DryPLE0(:,:,L)) * &
                               FV_Atm(1)%gridstruct%area(IS:IE,JS:JE) )
                     enddo
                  else
                     where (AIRDEN > 0.0_FVPRC)
                        TRACERS(:,:,:,N) = TRACERS(:,:,:,N) / AIRDEN
                     end where
                  endif
               enddo
               if (adjUseDpAdv) then
                  ALLOCATE( adjMuIn(IM,JM,LM,NQ), stat=STATUS )
                  VERIFY_(STATUS)
                  adjMuIn = TRACERS(:,:,:,1:NQ)
               endif
            endif
#endif

            if (ADVCORE_ADJ_DEBUG) then
               if (MAPL_Am_I_Root()) then
                  mfxMin = minval(MFX)
                  mfxMax = maxval(MFX)
                  mfyMin = minval(MFY)
                  mfyMax = maxval(MFY)
                  cxMin  = minval(CX)
                  cxMax  = maxval(CX)
                  cyMin  = minval(CY)
                  cyMax  = maxval(CY)
#ifdef ADJOINT
                  if (isAdjoint .and. associated(AIRDEN)) then
                     rhoDryMin = minval(AIRDEN)
                     rhoDryMax = maxval(AIRDEN)
                     write(*,*) 'ADVCORE_INPUT dt=', dt, ' mfxmin=', mfxMin, ' mfxmax=', mfxMax, &
                                ' mfymin=', mfyMin, ' mfymax=', mfyMax, ' cxmin=', cxMin,        &
                                ' cxmax=', cxMax, ' cymin=', cyMin, ' cymax=', cyMax,            &
                                ' airden_min=', rhoDryMin, ' airden_max=', rhoDryMax
                  else
                     write(*,*) 'ADVCORE_INPUT dt=', dt, ' mfxmin=', mfxMin, ' mfxmax=', mfxMax, &
                                ' mfymin=', mfyMin, ' mfymax=', mfyMax, ' cxmin=', cxMin,        &
                                ' cxmax=', cxMax, ' cymin=', cyMin, ' cymax=', cyMax
                  endif
#else
                  write(*,*) 'ADVCORE_INPUT dt=', dt, ' mfxmin=', mfxMin, ' mfxmax=', mfxMax, &
                             ' mfymin=', mfyMin, ' mfymax=', mfyMax, ' cxmin=', cxMin,        &
                             ' cxmax=', cxMax, ' cymin=', cyMin, ' cymax=', cyMax
#endif
               endif
            endif

            ! Run offline advection
#ifdef ADJOINT
            ! Adjoint tracers can be negative (they are sensitivities, not
            ! concentrations). Disable fill so fillz does not zero them out.
            if (isAdjoint) FV_Atm(1)%flagstruct%fill = .false.
            if (isAdjoint .and. advAdjDiscrete) then
               ! ---- discrete_linear: adjoint of the forward step, tracer by tracer
               if (.not. advAdjDiscReady) then
                  call adv_adj_discrete_init(FV_Atm(1)%gridstruct, FV_Atm(1)%flagstruct, FV_Atm(1)%bd, &
                                             FV_Atm(1)%domain, FV_Atm(1)%tile, FV_Atm(1)%npx,         &
                                             FV_Atm(1)%npy, FV_Atm(1)%npz, AK, BK, PTOP, abs(dt))
                  call backend_set_halo_ad(trim(advAdjHaloAd))
                  if (advAdjDiscCheck == 3) call backend_halo_test(FV_Atm(1)%npz)
                  advAdjDiscReady = .true.
               endif
               ! Inputs of the FORWARD step that this adjoint step reverses. GCHPctmEnv gives the
               ! reverse run the negated forward fluxes (ADV_ADJ_MASSFLUX: forward), and the swap
               ! above makes DryPLE0 / PLE0 the edges at t+dt and DryPLE1 / PLE1 those at t.
               ! To be confirmed on Pleiades with tools/check_conventions.py and
               ! ADV_ADJ_DISCRETE_CHECK (HANDOFF_PLEIADES.md).
               ALLOCATE( adjMfxF(IM,JM,LM), adjMfyF(IM,JM,LM), adjCxF(IM,JM,LM), adjCyF(IM,JM,LM) )
               ALLOCATE( adjPle0F(IM,JM,LM+1), adjPle1F(IM,JM,LM+1), adjLam(IM,JM,LM) )
               adjMfxF = -MFX; adjMfyF = -MFY; adjCxF = -CX; adjCyF = -CY
               if ( Use_Total_Air_Pressure > 0 ) then
                  adjPle0F = PLE1; adjPle1F = PLE0
               else
                  adjPle0F = DryPLE1; adjPle1F = DryPLE0
               endif
               do N=1,NQ
                  if (.not. isAdjointTracer(N)) cycle
                  adjLam = TRACERS(:,:,:,N)
                  if ( trim(advAdjDiscFixer) == 'linearized' ) then
                     ! forward partner: same name without the trailing '_ADJ'
                     lenAdj = len_trim(advTracers(N)%tName)
                     NF = 0
                     do L=1,NQ
                        if ( trim(advTracers(L)%tName) == advTracers(N)%tName(1:lenAdj-4) ) NF = L
                     enddo
                     if (NF == 0) then
                        write(*,*) 'ADVCORE discrete_linear: no forward tracer for ', trim(advTracers(N)%tName)
                        STATUS = 1
                        VERIFY_(STATUS)
                     endif
                     if ( trim(advTracers(NF)%tName) /= trim(advAdjTrajTracer) ) then
                        write(*,*) 'ADVCORE discrete_linear: trajectory is restored only for ', &
                                   trim(advAdjTrajTracer), ' (ADV_ADJ_TRAJ_TRACER), not ', trim(advTracers(NF)%tName)
                        STATUS = 1
                        VERIFY_(STATUS)
                     endif
                     ALLOCATE( adjQtraj(IM,JM,LM) )
                     adjQtraj = TRACERS(:,:,:,NF)
                     call adv_adj_discrete_adj(adjLam, adjMfxF, adjMfyF, adjCxF, adjCyF, &
                                               adjPle0F, adjPle1F, 'linearized', q_traj=adjQtraj)
                     DEALLOCATE( adjQtraj )
                  else
                     call adv_adj_discrete_adj(adjLam, adjMfxF, adjMfyF, adjCxF, adjCyF, &
                                               adjPle0F, adjPle1F, trim(advAdjDiscFixer))
                  endif
                  TRACERS(:,:,:,N) = adjLam
               enddo
               if (advAdjDiscCheck > 0) then
                  ! the same forward step by offline_tracer_advection and by the discrete code;
                  ! they must agree to round-off (checks the grid copy and the build, not the
                  ! convention, which both get from here)
                  ALLOCATE( adjQt4(IM,JM,LM,1), adjQt(IM,JM,LM), adjPleT(IM,JM,LM+1), adjPleAdvT(IM,JM,LM+1) )
                  NF = 0
                  if (advAdjDiscCheck == 2) then
                     ! real field: the forward trajectory q(T) of the first adjoint tracer (its forward partner)
                     do N=1,NQ
                        if (.not. isAdjointTracer(N)) cycle
                        lenAdj = len_trim(advTracers(N)%tName)
                        do L=1,NQ
                           if ( trim(advTracers(L)%tName) == advTracers(N)%tName(1:lenAdj-4) ) NF = L
                        enddo
                        exit
                     enddo
                  endif
                  if (NF > 0) then
                     adjQt = TRACERS(:,:,:,NF)
                  else
                     do L=1,LM
                        adjQt(:,:,L) = 4.0e-4_FVPRC * (1.0_FVPRC + 0.01_FVPRC*real(L,FVPRC)/real(LM,FVPRC))
                     enddo
                     adjQt(:,:,LM/2) = adjQt(:,:,LM/2) * (1.0_FVPRC + 0.1_FVPRC*FV_Atm(1)%gridstruct%area(IS:IE,JS:JE) &
                                       / maxval(FV_Atm(1)%gridstruct%area(IS:IE,JS:JE)))
                  endif
                  adjQt4(:,:,:,1) = adjQt
                  adjPleT = adjPle1F
                  call offline_tracer_advection(adjQt4, adjPle0F, adjPleT, adjMfxF, adjMfyF, adjCxF, adjCyF, &
                                                FV_Atm(1)%gridstruct, FV_Atm(1)%flagstruct, FV_Atm(1)%bd,  &
                                                FV_Atm(1)%domain, AK, BK, PTOP, FV_Atm(1)%npx,             &
                                                FV_Atm(1)%npy, FV_Atm(1)%npz, 1, abs(dt), adjPleAdvT)
                  call adv_adj_discrete_fwd(adjQt, adjMfxF, adjMfyF, adjCxF, adjCyF, adjPle0F, adjPle1F, 'full')
                  call ESMF_VMGetCurrent(vmRun, rc=STATUS)
                  VERIFY_(STATUS)
                  adjChkLoc = maxval(abs(REAL(adjQt - adjQt4(:,:,:,1), REAL8)))
                  call MAPL_CommsAllReduceMax(vmRun, sendbuf=adjChkLoc, recvbuf=adjChkGlob, cnt=1, rc=STATUS)
                  VERIFY_(STATUS)
                  if (MAPL_Am_I_Root()) write(*,'(a,es10.2)') &
                     ' ADVCORE discrete_linear check: max |offline_tracer_advection - discrete forward| =', adjChkGlob
                  if (NF > 0) then
                     ! trajectory alignment: one forward step of q(T) must give the q(T+dt) of the previous adjoint
                     ! call (advection-only runs). Reference: |q(T) - q(T+dt)| (one step of change).
                     if (allocated(adjQprev)) then
                        adjChkLoc = maxval(abs(REAL(adjQt - adjQprev, REAL8)))
                        call MAPL_CommsAllReduceMax(vmRun, sendbuf=adjChkLoc, recvbuf=adjChkGlob, cnt=1, rc=STATUS)
                        VERIFY_(STATUS)
                        adjChkLoc = maxval(abs(REAL(TRACERS(:,:,:,NF) - adjQprev, REAL8)))
                        call MAPL_CommsAllReduceMax(vmRun, sendbuf=adjChkLoc, recvbuf=adjChkRef, cnt=1, rc=STATUS)
                        VERIFY_(STATUS)
                        if (MAPL_Am_I_Root()) write(*,'(a,es10.2,a,es10.2)') &
                           ' ADVCORE trajectory check: max |fwd(q_traj(T)) - q_traj(T+dt)| =', adjChkGlob, &
                           '   reference max |q_traj(T) - q_traj(T+dt)| =', adjChkRef
                     else
                        allocate(adjQprev(IM,JM,LM))
                     endif
                     adjQprev = TRACERS(:,:,:,NF)
                  endif
                  DEALLOCATE( adjQt4, adjQt, adjPleT, adjPleAdvT )
               endif
               DEALLOCATE( adjMfxF, adjMfyF, adjCxF, adjCyF, adjPle0F, adjPle1F, adjLam )
            else
#endif
            if ( Use_Total_Air_Pressure > 0 ) then
               call offline_tracer_advection( TRACERS,              &
                                              PLE0,                 &
                                              PLE1,                 &
                                              MFX,                  &
                                              MFY,                  &
                                              CX,                   &
                                              CY,                   &
                                              FV_Atm(1)%gridstruct, &
                                              FV_Atm(1)%flagstruct, &
                                              FV_Atm(1)%bd,         &
                                              FV_Atm(1)%domain,     &
                                              AK,                   &
                                              BK,                   &
                                              PTOP,                 &
                                              FV_Atm(1)%npx,        &
                                              FV_Atm(1)%npy,        &
                                              FV_Atm(1)%npz,        &
                                              NAdv,                 &
                                              dt,                   &
                                              PLEAdv )
            else
               call offline_tracer_advection(TRACERS,              &
                                             DryPLE0,              &
                                             DryPLE1,              &
                                             MFX,                  &
                                             MFY,                  &
                                             CX,                   &
                                             CY,                   &
                                             FV_Atm(1)%gridstruct, &
                                             FV_Atm(1)%flagstruct, &
                                             FV_Atm(1)%bd,         &
                                             FV_Atm(1)%domain,     &
                                             AK,                   &
                                             BK,                   &
                                             PTOP,                 &
                                             FV_Atm(1)%npx,        &
                                             FV_Atm(1)%npy,        &
                                             FV_Atm(1)%npz,        &
                                             NAdv,                 &
                                             dt,                   &
                                             PLEAdv )
            endif
#ifdef ADJOINT
            endif   ! discrete_linear
            if (isAdjoint) FV_Atm(1)%flagstruct%fill = .true.
#endif
         endif

         ! Compute tracer mass post-advection (proxy). Done for all timesteps.
         ! Note that TMASS0 and MASS0 refer to proxy mass of first timestep and are
         ! printed and diffed every timestep for comparison.
         !-------------------------------------------------------------------
         if (chk_mass) then
            if ( Use_Total_Air_Pressure > 0 ) then
               MASS1 = g_sum( FV_Atm(1)%domain, &
                              PLE1(:,:,LM),                &
                              is,                          &
                              ie,                          &
                              js,                          &
                              je,                          &
                              FV_Atm(1)%ng,                &
                              FV_Atm(1)%gridstruct%area_64,&
                              1,                           &
                              .true.)
               call global_integral(TMASS1, TRACERS, PLE1, IM, JM, LM, NAdv)
            else
               MASS1 = g_sum( FV_Atm(1)%domain,            &
                              DryPLE1(:,:,LM),             &
                              is,                          &
                              ie,                          &
                              js,                          &
                              je,                          &
                              FV_Atm(1)%ng,                &
                              FV_Atm(1)%gridstruct%area_64,&
                              1,                           &
                              .true.)
               call global_integral(TMASS1, TRACERS, DryPLE1, IM,JM,LM,NQ)
            endif
            if (MASS1 /= 0.0) TMASS1=TMASS1/MASS1

            if (is_master()) then
               write(6,102)  MASS0, TMASS0(1)
               write(6,103)  MASS1, TMASS1(1)
               write(6,104) ( MASS1   - MASS0   )/ MASS0   , &
                    (TMASS1(1)-TMASS0(1))/TMASS0(1)
102            format('Tracer M0  : ',e21.14,' ',e21.14)
103            format('Tracer M1  : ',e21.14,' ',e21.14)
104            format('Tracer Mdif: ',e21.14,' ',e21.14)
            endif

            ! Optional post-advection debug : Sum additional values and print
            !! Area-weighted delta pressure
            !call global_integral_dp(dp_met, DryPLE1, IM,JM,LM)
            !
            !! Area-weighted species v/v * dp using met (proxy mass)
            !call global_integral_trmass_from_p(spc1_mass_dp_met, &
            !     Tracers(:,:,:,1), DryPLE1, IM,JM,LM)
            !
            !if (is_master()) then
            !   write(6,*) "Global area-weighted sums: sfc p met, spc1 vv*dp_met"
            !   write(6,105) dp_met, spc1_mass_dp_met
            !endif
105         format('Post-advection: ',e21.14,' ',e21.14)

         endif ! chk_mass

#ifdef ADJOINT
            ! Only undo the pre-advection scaling when advection actually ran:
            ! adjPreAdvSums and the AIRDEN division live inside AdvCore_Advection>0.
            if (isAdjoint .and. advOn .and. .not. advAdjDiscrete) then
               do N=1,NQ
                  if (.not. isAdjointTracer(N)) cycle
                  if (adjUseDpAdv) then
                     do L=1,LM
                        TRACERS(:,:,L,N) = ( TRACERS(:,:,L,N) * &
                             (PLEAdv(:,:,L+1) - PLEAdv(:,:,L)) + adjMuIn(:,:,L,N) * &
                             ( (DryPLE1(:,:,L+1) - DryPLE1(:,:,L)) -            &
                               (DryPLE0(:,:,L+1) - DryPLE0(:,:,L)) ) ) *        &
                             FV_Atm(1)%gridstruct%area(IS:IE,JS:JE)
                     enddo
                  elseif (adjUseAirMass) then
                     do L=1,LM
                        TRACERS(:,:,L,N) = TRACERS(:,:,L,N) * &
                             ( (DryPLE1(:,:,L+1) - DryPLE1(:,:,L)) * &
                               FV_Atm(1)%gridstruct%area(IS:IE,JS:JE) )
                     enddo
                  else
                     where (AIRDEN > 0.0_FVPRC)
                        TRACERS(:,:,:,N) = TRACERS(:,:,:,N) * AIRDEN
                     end where
                  endif
               enddo
               if (allocated(adjMuIn)) DEALLOCATE( adjMuIn )

               ! Restore pre-scaling global sum to undo calcScalingFactor applied
               ! inside offline_tracer_advection to adjoint tracers
               call ESMF_VMGetCurrent(vmRun, rc=STATUS)
               VERIFY_(STATUS)
               do N=1,NQ
                  if (.not. isAdjointTracer(N)) cycle
                  if (adjPreAdvSums(N) == 0.0_REAL8) cycle
                  adjLocalSum = SUM( REAL(TRACERS(:,:,:,N), REAL8) )
                  call MAPL_CommsAllReduceSum(vmRun, sendbuf=adjLocalSum, recvbuf=adjGlobalSum, &
                                             cnt=1, rc=STATUS)
                  VERIFY_(STATUS)
                  if (MAPL_Am_I_Root() .and. adjGlobalSum /= 0.0_REAL8) &
                     write(*,'(a,i0,a,es22.14)') ' ADVCORE adjoint global sum before/after advection, tracer ', &
                           N, ': ', adjPreAdvSums(N) / adjGlobalSum
                  if (adjGlobalSum /= 0.0_REAL8 .and. advAdjRescale /= 0) then
                     TRACERS(:,:,:,N) = TRACERS(:,:,:,N) * REAL(adjPreAdvSums(N) / adjGlobalSum, FVPRC)
                  endif
               enddo
               DEALLOCATE( adjPreAdvSums, stat=STATUS )
               VERIFY_(STATUS)

               if (ADVCORE_ADJ_DEBUG) then
                  call ESMF_VMGetCurrent(vmRun, rc=STATUS)
                  VERIFY_(STATUS)
                  do N=1,NQ
                     if (.not. isAdjointTracer(N)) cycle
                     adjLocalSum = SUM( REAL(TRACERS(:,:,:,N), REAL8) )
                     call MAPL_CommsAllReduceSum(vmRun, sendbuf=adjLocalSum, recvbuf=adjGlobalSum, &
                                                 cnt=1, rc=STATUS)
                     VERIFY_(STATUS)
                     if (MAPL_Am_I_Root()) then
                        write(*,*) 'ADVCORE_ADJ_GLOBAL_SUM stage=after_mul_airden idx=', N, &
                                   ' name=', trim(advTracers(N)%tName), ' sum=', adjGlobalSum
                     endif
                  enddo
               endif
            endif
#endif

         ! If using total air pressure then convert all tracers from kg/kg total
         ! to kg/kg dry for use in GEOS-Chem. Use the post-advection specific
         ! humidity tracer for the conversion.
         if ( Use_Total_Air_Pressure > 0 ) then
            do N=1,NQ
               TRACERS(:,:,:,N) = TRACERS(:,:,:,N) &
                                  / (1.0d0 - TRACERS(:,:,:,NQ+1))
            enddo
         endif

         ! Go through the bundle copying tracers back to the bundle.
         !-----------------------------------------------------------------

         do N=1,NQ
            if (advTracers(N)%is_r4) then
               advTracers(N)%content_r4 = TRACERS(:,:,:,N)
            else
               advTracers(N)%content    = TRACERS(:,:,:,N)
            end if

            !-----------------------------------------------
            !--> Fill Export States
            !--> This section is used for diagnostics only.
            !--> It has no effect on CTM experiments.
            !-----------------------------------------------
            if (N<=ntracers) then
               write(myTracer, "('TEST_TRACER',i5.5)") N-1
               call MAPL_GetPointer(EXPORT, temp3D, TRIM(myTracer), rc=status)
               VERIFY_(STATUS)
               if (associated(temp3D)) temp3D = TRACERS(:,:,:,N)
            endif
         enddo

         ! Deallocate the list of tracers
         !----------------------------------------------------------------
         DEALLOCATE( TRACERS,stat=STATUS )
         VERIFY_(STATUS)
#ifdef ADJOINT
         DEALLOCATE( isAdjointTracer, stat=STATUS )
         VERIFY_(STATUS)
#endif

      end if ! NQ > 0

      ! GCHP: pressure edge exports
      if ( Use_Total_Air_Pressure < 1 ) then
         call MAPL_GetPointer ( EXPORT, eDryPLE, 'DryPLE', ALLOC=.TRUE., RC=STATUS )
         eDryPLE(:,:,:) = DryPLE1(:,:,:)

         ! Set R4 exports as temporary work-around for MAPL 2.55 bug
         call MAPL_GetPointer ( EXPORT, eDryPLE_R4, 'DryPLE_R4', &
                                NotFoundOK=.TRUE., _RC )
         IF ( ASSOCIATED(eDryPLE_R4) ) eDryPLE_R4(:,:,:) = eDryPLE(:,:,:)

      endif
      call MAPL_GetPointer ( EXPORT, ePLE, 'PLE', ALLOC=.TRUE., _RC )
      ePLE(:,:,:) = PLE1(:,:,:)
      call MAPL_GetPointer ( EXPORT, ePLEadv, 'PLEadv', ALLOC=.TRUE., _RC )
      ePLEadv(:,:,:) = PLEadv(:,:,:)

      ! GCHP: Set R4 exports as temporary work-around for MAPL 2.55 bug      
      call MAPL_GetPointer ( EXPORT, ePLE_R4, 'PLE_R4', NotFoundOK=.TRUE., _RC )
      IF ( ASSOCIATED(ePLE_R4) ) ePLE_R4(:,:,:) = ePLE(:,:,:)
      call MAPL_GetPointer ( EXPORT, ePLEadv_R4, 'PLEadv_R4', NotFoundOK=.TRUE., _RC )
      IF ( ASSOCIATED(ePLEadv_R4) ) ePLEadv_R4(:,:,:) = ePLEadv(:,:,:)

      deallocate( advTracers, stat=STATUS )
      VERIFY_(STATUS)
      DEALLOCATE( AK ,stat=STATUS )
      VERIFY_(STATUS)
      DEALLOCATE( BK ,stat=STATUS )
      VERIFY_(STATUS)

      DEALLOCATE( PLE0   )
      DEALLOCATE( PLE1   )
      DEALLOCATE( PLEAdv )
      DEALLOCATE(  MFX   )
      DEALLOCATE(  MFY   )
      DEALLOCATE(   CX   )
      DEALLOCATE(   CY   )

      if ( Use_Total_Air_Pressure > 0 ) then
         DEALLOCATE( SPHU0 )
      else
         DEALLOCATE( DryPLE0 )
         DEALLOCATE( DryPLE1 )
      endif

#ifdef ADJOINT
      if (associated(AIRDEN)) then
         DEALLOCATE( AIRDEN )
      endif
#endif

      call MAPL_TimerOff(MAPL,"RUN")
      call MAPL_TimerOff(MAPL,"TOTAL")

      !WMP  end if ! AdvCore_Advection

      firstRun=.false.

      RETURN_(ESMF_SUCCESS)

      end subroutine Run
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Finalize - user supplied finalize routine
!
! !INTERFACE:
!
  subroutine Finalize(GC, IMPORT, EXPORT, CLOCK, RC)
!
! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
      type(ESMF_State),    intent(inout) :: IMPORT ! Import state
      type(ESMF_State),    intent(inout) :: EXPORT ! Export state
      type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
!
! !OUTPUT PARAMETERS:
      integer, optional,   intent(  out) :: RC     ! Error code
!
! !DESCRIPTION:
!    Finalize merely destroys the FVadv object that was created in Initialize
!    and releases the space for the persistent data .
!
!EOP
!=============================================================================
!BOC
! !LOCAL VARIABLES:

      character(len=ESMF_MAXSTR)    :: IAm
      integer                       :: STATUS
      character(len=ESMF_MAXSTR)    :: COMP_NAME

      ! Get my name and set-up traceback handle
      ! ---------------------------------------

      Iam = 'Finalize'
      call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
      VERIFY_(STATUS)
      Iam = trim(COMP_NAME) // TRIM(Iam)

      ! Clean up FV if AdvCore is running without FV3_DynCoreIsRunning
      !--------------------------------------------------
      if (.NOT. FV3_DynCoreIsRunning) then
         call fv_end(FV_Atm, grids_on_my_pe, .false.)
      endif

      call MAPL_GenericFinalize(GC, IMPORT, EXPORT, CLOCK, RC)
      VERIFY_(STATUS)

      RETURN_(ESMF_SUCCESS)
      end subroutine Finalize


subroutine global_integral (QG,Q,PLE,IM,JM,KM,NQ)

      real(FVPRC), intent(OUT)   :: QG(NQ)
      real(FVPRC), intent(IN)    :: Q(IM,JM,KM,NQ)
      real(FVPRC), intent(IN)    :: PLE(IM,JM,KM+1)
      integer,     intent(IN)    :: IM,JM,KM,NQ

      ! Locals
      integer   :: k,n
      real(REAL8), allocatable ::    dp(:,:,:)
      real(FVPRC), allocatable :: qsum1(:,:)

      allocate(    dp(im,jm,km) )
      allocate( qsum1(im,jm)    )

      ! Compute Pressure Thickness
      ! --------------------------
      do k=1,KM
         dp(:,:,k) = PLE(:,:,k+1)-PLE(:,:,k)
      enddo

      ! Loop over Tracers
      ! -----------------
      do n=1,NQ
         qsum1(:,:) = 0.d0
         do k=1,KM
            qsum1(:,:) = qsum1(:,:) + Q(:,:,k,n)*dp(:,:,k)
         enddo
         qg(n) = g_sum( FV_Atm(1)%domain,            &
                        qsum1,                       &
                        is,                          &
                        ie,                          &
                        js,                          &
                        je,                          &
                        FV_Atm(1)%ng,                &
                        FV_Atm(1)%gridstruct%area_64,&
                        1,                           &
                        .true.)
      enddo

      deallocate( dp )
      deallocate( qsum1 )

end subroutine global_integral

subroutine global_integral_trmass_from_p (QG,Q,PLE,IM,JM,KM)

      real(FVPRC), intent(OUT)   :: QG
      real(FVPRC), intent(IN)    :: Q(IM,JM,KM)
      real(FVPRC), intent(IN)    :: PLE(IM,JM,KM+1)
      integer,     intent(IN)    :: IM,JM,KM

      ! Locals
      integer   :: k
      real(REAL8), allocatable ::    dp(:,:,:)
      real(FVPRC), allocatable :: qsum1(:,:)

      allocate(    dp(im,jm,km) )
      allocate( qsum1(im,jm)    )

      ! Pressure thickness [Pa]
      do k=1,KM
         dp(:,:,k) = PLE(:,:,k+1)-PLE(:,:,k)
      enddo

      ! Column sum
      qsum1(:,:) = 0.d0
      do k=1,KM
         qsum1(:,:) = qsum1(:,:) + Q(:,:,k)*dp(:,:,k)
      enddo

      ! Global sum
      qg = g_sum( FV_Atm(1)%domain,            &
           qsum1,                       &
           is,                          &
           ie,                          &
           js,                          &
           je,                          &
           FV_Atm(1)%ng,                &
           FV_Atm(1)%gridstruct%area_64,&
           0,                           & ! do not divide total sum by total aread
           .true.)

      deallocate( dp )
      deallocate( qsum1 )

end subroutine global_integral_trmass_from_p

subroutine global_integral_trmass_from_dp (QG,Q,DP,IM,JM,KM)

      real(FVPRC), intent(OUT)   :: QG
      real(FVPRC), intent(IN)    :: Q(IM,JM,KM)
      real(FVPRC), intent(IN)    :: DP(IM,JM,KM)
      integer,     intent(IN)    :: IM,JM,KM

      ! Locals
      integer   :: k
      real(FVPRC), allocatable :: qsum1(:,:)

      allocate( qsum1(im,jm)    )

      ! columm sum
      qsum1(:,:) = 0.d0
      do k=1,KM
         qsum1(:,:) = qsum1(:,:) + Q(:,:,k)*dp(:,:,k)
      enddo

      ! global sum
      qg = g_sum( FV_Atm(1)%domain,            &
           qsum1,                       &
           is,                          &
           ie,                          &
           js,                          &
           je,                          &
           FV_Atm(1)%ng,                &
           FV_Atm(1)%gridstruct%area_64,&
           0,                           & ! do not divide total sum by total area
           .true.)

      deallocate( qsum1 )

end subroutine global_integral_trmass_from_dp

subroutine global_integral_dp (QG,PLE,IM,JM,KM)

      real(REAL8), intent(OUT)   :: QG
      real(FVPRC), intent(IN)    :: PLE(IM,JM,KM+1)
      integer,     intent(IN)    :: IM,JM,KM

      ! Locals
      integer   :: k,n
      real(REAL8), allocatable ::    dp(:,:,:)
      real(FVPRC), allocatable :: qsum1(:,:)

      allocate(    dp(im,jm,km) )
      allocate( qsum1(im,jm)    )

      ! Pressure thickness [Pa]
      do k=1,KM
         dp(:,:,k) = PLE(:,:,k+1)-PLE(:,:,k)
      enddo

      ! Column sum
      qsum1(:,:) = 0.d0
      do k=1,KM
         qsum1(:,:) = qsum1(:,:) + dp(:,:,k)
      enddo

      ! Global sum
      qg = g_sum( FV_Atm(1)%domain,            &
           qsum1,                       &
           is,                          &
           ie,                          &
           js,                          &
           je,                          &
           FV_Atm(1)%ng,                &
           FV_Atm(1)%gridstruct%area_64,&
           0,                           & ! do not divide total sum by total area
           .true.)

      deallocate( dp )
      deallocate( qsum1 )

end subroutine global_integral_dp

subroutine global_integral_dp_rst (QG,DP,IM,JM,KM)

      real(REAL8), intent(OUT)   :: QG
      real(FVPRC), intent(IN)    :: DP(IM,JM,KM)
      integer,     intent(IN)    :: IM,JM,KM

      ! Locals
      integer   :: k
      real(FVPRC), allocatable :: qsum1(:,:)

      allocate( qsum1(im,jm)    )

      ! column sum [Pa]
      qsum1(:,:) = 0.d0
      do k=1,KM
         qsum1(:,:) = qsum1(:,:) + dp(:,:,k)
      enddo

      ! global sum
      qg = g_sum( FV_Atm(1)%domain,            &
           qsum1,                       &
           is,                          &
           ie,                          &
           js,                          &
           je,                          &
           FV_Atm(1)%ng,                &
           FV_Atm(1)%gridstruct%area_64,&
           0,                           & ! do not divide total sum by total area
           .true.)

      deallocate( qsum1 )

end subroutine global_integral_dp_rst

subroutine global_integral_vv (QG,Q,IM,JM,KM)

      real(FVPRC), intent(OUT)   :: QG
      real(FVPRC), intent(IN)    :: Q(IM,JM,KM)
      integer,     intent(IN)    :: IM,JM,KM

      ! Locals
      integer   :: k
      real(FVPRC), allocatable :: qsum1(:,:)

      allocate( qsum1(im,jm)    )

      ! column sum
      qsum1(:,:) = 0.d0
      do k=1,KM
         qsum1(:,:) = qsum1(:,:) + Q(:,:,k)
      enddo

      ! global sum
      qg = g_sum( FV_Atm(1)%domain,            &
           qsum1,                       &
           is,                          &
           ie,                          &
           js,                          &
           je,                          &
           FV_Atm(1)%ng,                &
           FV_Atm(1)%gridstruct%area_64,&
           0,                           & ! do not divide total sum by total area
           .true.)

      deallocate( qsum1 )

end subroutine global_integral_vv

subroutine scale_tracers_by_pressure_ratio (Q,PLE,DP,IM,JM,KM,NQ)

      real(FVPRC), intent(INOUT)       :: Q(IM,JM,KM,NQ)
      real(FVPRC), intent(IN)          :: PLE(IM,JM,KM+1)
      real(FVPRC), intent(IN)          :: DP(IM,JM,KM)
      integer,     intent(IN)          :: IM,JM,KM,NQ

      ! Locals
      integer   :: i,j,k,n

      ! Loop over levels and tracers and apply scaling
      do n=1,NQ
      do k=1,KM
      do j=1,JM
      do i=1,IM
         Q(i,j,k,n) = Q(i,j,k,n) * dp(i,j,k) / (PLE(i,j,k+1)-PLE(i,j,k))
      enddo
      enddo
      enddo
      enddo

end subroutine scale_tracers_by_pressure_ratio

end module AdvCore_GridCompMod
