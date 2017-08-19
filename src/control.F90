!-----------------------------------------------------------------------
!
! MODULE: control_module
!> @brief
!> This module contains the variables that control SNAP's solver
!> routines. This includes the time-dependent variables. It also
!> includes some KBA scheduling variables.
!
!-----------------------------------------------------------------------

MODULE control_module

  USE global_module, ONLY: i_knd, r_knd, l_knd, zero, one
#ifdef SHM
  USE plib_module
#endif
  IMPLICIT NONE

  PUBLIC

  SAVE
!_______________________________________________________________________
!
! Module Input Variables
!
! epsi     - convergence criterion
! iitm     - max inner iterations
! oitm     - max outer iterations
! timedep  - 0/1=no/yes perform a time-dependent calculation
! tf       - final time
! nsteps   - number of time steps to cover the ts -> tf range
!
! swp_typ  - 0/1=standard order/mini-KBA sweep
! cor_swp  - 0/1=no/yes concurrent octant mesh sweeps (corner starts)
!
! angcpy   - 1/2 copies of the time-edge angular flux
!
! it_det   - 0/1=no/yes full iteration details
! soloutp  - 0/1=no/yes print single k-plane solution to output file
! kplane   - 0/1+=default mid-plane/k-plane index to print with soloutp
! popout   - 0/1/2=no/final only/every cycle print population data to
!            output file
! fluxp    - 0/1/2=print none/scalar flux/all flux moments to file
!
! fixup    - 0/1=no/yes perform flux fixup
!_______________________________________________________________________

  INTEGER(i_knd) :: iitm=5, oitm=100, timedep=0, nsteps=1, swp_typ=0,  &
    cor_swp=0, angcpy=1, it_det=0, soloutp=0, kplane=0, popout=0,      &
    fluxp=0, fixup=1

  REAL(r_knd) :: epsi=1.0E-4_r_knd, tf=zero
!_______________________________________________________________________
!
! Run-time variables
!
! dt       - time-step size
!
! tolr      - parameter, small number used for determining how to
!             compute flux error
! dfmxi(ng) - max error of inner iteration
! dfmxo     - max error of outer iteration
!
! inrdone(ng)  - logical for inners being complete
! otrdone      - logical for outers being complete
! update_ptr   - true/false update the ptr_out array
!
! ncor             - number of corners from which sweeps begin
! last_oct         - last octant to be swept
! corner_sch(2,4)  - corner scheduling control array
! yzstg(ncor)      - KBA stage in yz plane per starting corner
!_______________________________________________________________________

#ifdef SHM
  LOGICAL(l_knd), DIMENSION(:), POINTER :: update_ptr
  LOGICAL(l_knd), DIMENSION(:), POINTER :: inrdone, otrdone
#else
  LOGICAL(l_knd), DIMENSION(1) :: update_ptr
  LOGICAL(l_knd), DIMENSION(1) :: otrdone
  LOGICAL(l_knd), ALLOCATABLE, DIMENSION(:) :: inrdone
#endif
  LOGICAL(l_knd) :: otrdone_l ! *PIP: copy used after deallocate.
  INTEGER(i_knd) :: ncor, last_oct, corner_sch(2,4)

#ifdef SHM
  INTEGER(i_knd), DIMENSION(:), POINTER :: gcy ! *PIP: shared copy of cy
  INTEGER(i_knd), DIMENSION(:), POINTER :: yzstg
  REAL(r_knd), DIMENSION(:), POINTER :: dfmxo
#else
  INTEGER(i_knd), ALLOCATABLE, DIMENSION(:) :: yzstg
  INTEGER(i_knd), DIMENSION(1) :: gcy ! DEBUG USE
  REAL(r_knd), DIMENSION(1) :: dfmxo
#endif
  REAL(r_knd) :: dt

  REAL(r_knd), PARAMETER :: tolr=1.0E-12_r_knd

#ifdef SHM
  REAL(r_knd), DIMENSION(:), POINTER :: dfmxi
#else
  REAL(r_knd), ALLOCATABLE, DIMENSION(:) :: dfmxi
#endif

  CONTAINS


  SUBROUTINE control_allocate ( ng, ierr )

!-----------------------------------------------------------------------
!
! Allocate control module variables.
!
!-----------------------------------------------------------------------

    INTEGER(i_knd), INTENT(IN) :: ng

    INTEGER(i_knd), INTENT(OUT) :: ierr
!_______________________________________________________________________

    ierr = 0
#ifdef SHM
    !* PIP: only updated at setup.
    ALLOCATE(yzstg(ncor))
    !* PIP: used in sweep
    CALL shm_allocate(ng, dfmxi, "dfmxi")
    CALL shm_allocate(ng, inrdone, "inrdone")
    !* PIP: shared single logical variable
    CALL shm_allocate(1, otrdone, "otrdone")
    CALL shm_allocate(1, dfmxo, "dfmxo")
    CALL shm_allocate(1, gcy, "gcy")
    CALL shm_allocate(1, update_ptr, "update_ptr")
#else
    ALLOCATE( yzstg(ncor), dfmxi(ng), inrdone(ng), STAT=ierr )
    IF ( ierr /= 0 ) RETURN
#endif

    yzstg = 0

    dfmxi = -one
    inrdone = .FALSE.
    dfmxo(1) = -one
    otrdone(1) = .FALSE.
    update_ptr(1) = .FALSE.
!_______________________________________________________________________
!_______________________________________________________________________

#ifdef SHM
  CALL shm_barrier
#endif
  CALL plib_dbg_rootprintf("control_allocate done")

  END SUBROUTINE control_allocate


  SUBROUTINE control_deallocate

!-----------------------------------------------------------------------
!
! Deallocate control module arrays.
!
!-----------------------------------------------------------------------
!_______________________________________________________________________

#ifdef SHM
    DEALLOCATE( yzstg )
    CALL shm_deallocate(dfmxi)
    CALL shm_deallocate(inrdone)
    CALL shm_deallocate(otrdone)
    CALL shm_deallocate(dfmxo)
    CALL shm_deallocate(gcy)
    CALL shm_deallocate(update_ptr)
#else
    DEALLOCATE( yzstg, dfmxi, inrdone )
#endif
!_______________________________________________________________________
!_______________________________________________________________________

  END SUBROUTINE control_deallocate


END MODULE control_module
