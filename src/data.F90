!-----------------------------------------------------------------------
!
! MODULE: data_module
!> @brief
!> This module contains the variables and setup subroutines for the mock
!> cross section data. It establishes the number of groups and
!> constructs the cross section arrays.
!
!-----------------------------------------------------------------------

MODULE data_module

  USE global_module, ONLY: i_knd, r_knd, zero

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
! ng       - number of groups
! mat_opt  - material layout, 0/1/2=homogeneous/center/corner, with
!            two materials, and material 2 nowhere/center/corner
! src_opt  - source layout, 0/1/2=homogenous/src-center/src-corner=
!            source everywhere/center of problem/corner, strength=10.0
! scatp    - 0/1=no/yes print the full scattering matrix to file 'slgg'
!_______________________________________________________________________

  INTEGER(i_knd) :: ng=1, mat_opt=0, src_opt=0, scatp=0
!_______________________________________________________________________
!
! Run-time variables
!
! v(ng)         - mock velocity array
! nmat          - number of materials
! mat(nx,ny,nz) - material identifier array
!
! qi(nx,ny,nz,ng)             - fixed source array for src_opt<3
! qim(nang,nx,ny,nz,noct,ng)  - fixed source array for src_opt>=3
!
! sigt(nmat,ng)          - total interaction
! siga(nmat,ng)          - absorption
! sigs(nmat,ng)          - scattering, total
! slgg(nmat,nmom,ng,ng)  - scattering matrix, all moments/groups
! vdelt(ng)              - time-absorption coefficient
!_______________________________________________________________________

  INTEGER(i_knd) :: nmat=1
#ifdef SHM
  INTEGER(i_knd), DIMENSION(:,:,:), POINTER :: mat

  REAL(r_knd), DIMENSION(:), POINTER :: v, vdelt

  REAL(r_knd), DIMENSION(:,:), POINTER :: sigt, siga, sigs

  REAL(r_knd), DIMENSION(:,:,:,:), POINTER :: qi, slgg

  REAL(r_knd), DIMENSION(:,:,:,:,:,:), POINTER :: qim
#else
  INTEGER(i_knd), ALLOCATABLE, DIMENSION(:,:,:) :: mat

  REAL(r_knd), ALLOCATABLE, DIMENSION(:) :: v, vdelt

  REAL(r_knd), ALLOCATABLE, DIMENSION(:,:) :: sigt, siga, sigs

  REAL(r_knd), ALLOCATABLE, DIMENSION(:,:,:,:) :: qi, slgg

  REAL(r_knd), ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: qim
#endif

  CONTAINS


  SUBROUTINE data_allocate ( nx, ny, nz, nmom, nang, noct, timedep,    &
    istat )

!-----------------------------------------------------------------------
!
! Allocate data_module arrays.
!
!-----------------------------------------------------------------------

    INTEGER(i_knd), INTENT(IN) :: nx, ny, nz, nmom, nang, noct, timedep

    INTEGER(i_knd), INTENT(INOUT) :: istat
!_______________________________________________________________________
!
!   Establish number of materials according to mat_opt
!_______________________________________________________________________

    IF ( mat_opt > 0 ) nmat = 2
!_______________________________________________________________________
!
!   Allocate velocities
!_______________________________________________________________________

    istat = 0

    IF ( timedep == 1 ) THEN
#ifdef SHM
      CALL shm_allocate(ng, v, "v")
#else
      ALLOCATE( v(ng), STAT=istat )
#endif
    ELSE
      ALLOCATE( v(0), STAT=istat )
    END IF
    IF ( istat /= 0 ) RETURN

    v = zero
!_______________________________________________________________________
!
!   Allocate the material identifier array. ny and nz are 1 if not
!   2-D/3-D.
!_______________________________________________________________________

#ifdef SHM
    CALL shm_allocate(nx, ny, nz, mat, "mat")
#else
    ALLOCATE( mat(nx,ny,nz), STAT=istat )
    IF ( istat /= 0 ) RETURN
#endif

    mat = 1
!_______________________________________________________________________
!
!   Allocate the fixed source array. If src_opt < 3, allocate the qi
!   array, not the qim. Do the opposite (store the full angular copy) of
!   the source, qim, if src_opt>=3 (MMS). Allocate array not used to 0.
!   ny and nz are 1 if not 2-D/3-D.
!_______________________________________________________________________

    IF ( src_opt < 3 ) THEN
#ifdef SHM
      CALL shm_allocate(nx,ny,nz,ng,qi,"qi")
#else
      ALLOCATE( qi(nx,ny,nz,ng), qim(0,0,0,0,0,0), STAT=istat )
      IF ( istat /= 0 ) RETURN
#endif
      qi = zero
    ELSE
#ifdef SHM
      CALL shm_allocate(nx,ny,nz,ng,qi,"qi")
      CALL shm_allocate(nang,nx,ny,nz,noct,ng,qim,"qim")
#else
      ALLOCATE( qi(nx,ny,nz,ng), qim(nang,nx,ny,nz,noct,ng),           &
        STAT=istat )
      IF ( istat /= 0 ) RETURN
#endif
      qi = zero
      qim = zero
    END IF
!_______________________________________________________________________
!
!   Allocate mock cross sections
!_______________________________________________________________________
#ifdef SHM
    CALL shm_allocate(nmat,ng,sigt,"sigt")
    CALL shm_allocate(nmat,ng,siga,"siga")
    CALL shm_allocate(nmat,ng,sigs,"sigs")
    CALL shm_allocate(nmat,nmom,ng,ng,slgg,"slgg")
#else
    ALLOCATE( sigt(nmat,ng), siga(nmat,ng), sigs(nmat,ng),             &
      slgg(nmat,nmom,ng,ng), STAT=istat )
    IF ( istat /= 0 ) RETURN
#endif

    sigt = zero
    siga = zero
    sigs = zero
    slgg = zero
!_______________________________________________________________________
!
!   Allocate the vdelt array
!_______________________________________________________________________

#ifdef SHM
    CALL shm_allocate(ng,vdelt,"vdelt")
#else
    ALLOCATE( vdelt(ng), STAT=istat )
    IF ( istat /= 0 ) RETURN
#endif

    vdelt = zero
!_______________________________________________________________________
!_______________________________________________________________________

#ifdef SHM
  CALL shm_barrier
#endif
  CALL plib_dbg_rootprintf("data_allocate done")

  END SUBROUTINE data_allocate


  SUBROUTINE data_deallocate

!-----------------------------------------------------------------------
!
! Deallocate data_module arrays.
!
!-----------------------------------------------------------------------
!_______________________________________________________________________

#ifdef SHM
    CALL shm_deallocate(v)
    CALL shm_deallocate(mat)
    CALL shm_deallocate(qi)
    CALL shm_deallocate(qim)
    CALL shm_deallocate(sigt)
    CALL shm_deallocate(siga)
    CALL shm_deallocate(sigs)
    CALL shm_deallocate(slgg)
    CALL shm_deallocate(vdelt)
#else
    DEALLOCATE( v )
    DEALLOCATE( mat )
    DEALLOCATE( qi, qim )
    DEALLOCATE( sigt, siga, sigs )
    DEALLOCATE( slgg )
    DEALLOCATE( vdelt )
#endif
!_______________________________________________________________________
!_______________________________________________________________________

  END SUBROUTINE data_deallocate


END MODULE data_module
