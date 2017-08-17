!-----------------------------------------------------------------------
!
! MODULE: plib_module
!> @brief
!> This module contains the variables that control parallel
!> decomposition and the subroutines for parallel environment setup.
!> Only module that requires MPI library interaction except for
!> time_module (MPI_WTIME).
!
!-----------------------------------------------------------------------

MODULE plib_module

  USE global_module, ONLY: i_knd, r_knd, l_knd

  USE time_module, ONLY: wtime

  USE iso_c_binding

#ifdef MPI
  USE mpi
#endif

  IMPLICIT NONE

  PUBLIC

  INTERFACE glmax
    MODULE PROCEDURE glmax_i, glmax_d, glmax_d_1d
  END INTERFACE glmax

  INTERFACE glmin
    MODULE PROCEDURE glmin_i, glmin_d
  END INTERFACE glmin

  INTERFACE glsum
    MODULE PROCEDURE glsum_d
  END INTERFACE glsum

  INTERFACE rtsum
    MODULE PROCEDURE rtsum_d_1d
  END INTERFACE rtsum

  INTERFACE bcast
    MODULE PROCEDURE bcast_i_scalar, bcast_i_1d, bcast_d_scalar,       &
      bcast_d_1d, bcast_l_scalar
  END INTERFACE bcast

  INTERFACE psend
    MODULE PROCEDURE psend_d_2d, psend_d_3d
  END INTERFACE psend

  INTERFACE isend
    MODULE PROCEDURE isend_d_3d
  END INTERFACE isend

  INTERFACE precv
    MODULE PROCEDURE precv_d_2d, precv_d_3d
  END INTERFACE precv

#ifdef SHM
  INTERFACE shm_allocate
    MODULE PROCEDURE shm_allocate_d_1d, shm_allocate_d_2d,             &
      shm_allocate_d_3d, shm_allocate_d_4d, shm_allocate_d_5d,         &
      shm_allocate_d_6d, shm_allocate_i_1d, shm_allocate_i_2d,         &
      shm_allocate_i_3d, shm_allocate_l_1d
  END INTERFACE shm_allocate

  INTERFACE shm_deallocate
    MODULE PROCEDURE shm_deallocate_d_1d, shm_deallocate_d_2d,         &
      shm_deallocate_d_3d, shm_deallocate_d_4d, shm_deallocate_d_5d,   &
      shm_deallocate_d_6d, shm_deallocate_i_1d, shm_deallocate_i_2d,   &
      shm_deallocate_i_3d, shm_deallocate_l_1d
  END INTERFACE shm_deallocate
#endif

  SAVE
!_______________________________________________________________________
!
! Module Input Variables
!
! npey     - Number of MPI processes in the y-direction
! npez     - Number of MPI processes in the z-direction
! ichunk   - Size of work chunks in the x-direction
! nthreads - Number of OpenMP threads for energy parallelism
! nnested  - number of nested threads
!_______________________________________________________________________

  INTEGER(i_knd) :: npey=1, npez=1, ichunk=4, nthreads=1, nnested=1
!_______________________________________________________________________
!
! Run-time variables
!
! Note: all ranks are zero based
!
! root     - root process for comm_snap, 0
!
! nproc    - Number of MPI processes
! iproc    - Rank of calling process in base communicator
!
! comm_snap  - base communicator, duplicated from MPI_COMM_WORLD
! comm_space - SDD communicator, ndimen-1 grid for 2-D (x-y) or
!              3-D (x-y-z) problems. Non-existent for 1-D (x) problems.
! sproc      - Rank of calling process in comm_space
!
! ycomm      - y-dimension process communicator
! zcomm      - z-dimension process communicator
! yproc      - PE column in SDD 2-D PE mesh (comm_space)
! zproc      - PE row in SDD 2-D PE mesh (comm_space)
! firsty     - logical determining if lowest yproc
! lasty      - logical determining if highest yproc
! firstz     - logical determining if lowest zproc
! lastz      - logical determining if highest zproc
! ylop       - rank of preceding yproc in ycomm
! yhip       - rank of succeeding yproc in ycomm
! zlop       - rank of preceding zproc in zcomm
! zhip       - rank of succeeding zproc in zcomm
!
! thread_level       - level of MPI thread support
! thread_single      - MPI_THREAD_SINGLE
! thread_funneled    - MPI_THREAD_FUNNELED
! thread_serialized  - MPI_THREAD_SERIALIZED
! thread_multiple    - MPI_THREAD_MULTIPLE
! lock(nthreads)     - OpenMP lock for each thread
!
! do_nested - true/false use nested threading
!
! use_lock  - true/false apply lock to threads MPI communications
!             during sweep
!
! pce    - Parallel computational efficiency of the run
!_______________________________________________________________________

  INTEGER(i_knd), PARAMETER :: root=0

  INTEGER(i_knd) :: nproc, iproc, comm_snap, comm_space, sproc, ycomm, &
    zcomm, yproc, zproc, ylop, yhip, zlop, zhip, thread_level,         &
    thread_single, thread_funneled, thread_serialized, thread_multiple,&
    max_threads, ynproc, znproc,snproc

  INTEGER(i_knd) :: wnproc, wiproc

#ifdef SHM
  INTEGER(i_knd) :: comm_maxshm, comm_shm, comm_shmgrp
  INTEGER(i_knd) :: maxshm_nproc, maxshm_iproc,                        &
    shm_nproc, shm_iproc
  LOGICAL(l_knd) :: is_shm_master
#endif
  LOGICAL(l_knd) :: firsty, lasty, firstz, lastz, do_nested,           &
    use_lock=.FALSE.

  REAL(r_knd) :: pce

#ifdef OPENMP
  INCLUDE 'omp_lib.h'

  INTEGER(OMP_LOCK_KIND), ALLOCATABLE, DIMENSION(:) :: lock
#endif


  CONTAINS


#ifdef MPI

  SUBROUTINE pinit ( t1 )

!-----------------------------------------------------------------------
!
! Initialize the MPI process environment. Replicate MPI_COMM_WORLD to
! local communicator. Get each process its rank and total size.
!
!-----------------------------------------------------------------------

    REAL(r_knd), INTENT(OUT) :: t1
!_______________________________________________________________________
!
!   Local variables
!_______________________________________________________________________

    INTEGER(i_knd) :: ierr
!_______________________________________________________________________
!
!   Set thread support levels
!_______________________________________________________________________

    thread_single     = MPI_THREAD_SINGLE
    thread_funneled   = MPI_THREAD_FUNNELED
    thread_serialized = MPI_THREAD_SERIALIZED
    thread_multiple   = MPI_THREAD_MULTIPLE

#ifdef SHM
    comm_shm = MPI_COMM_NULL
    comm_maxshm = MPI_COMM_NULL
    comm_shmgrp = MPI_COMM_NULL
#endif
    comm_snap = MPI_COMM_NULL
    comm_space = MPI_COMM_NULL
    ycomm = MPI_COMM_NULL
    zcomm = MPI_COMM_NULL
!_______________________________________________________________________
!
!   Initialize MPI and return available thread support level. Prefer
!   thread_multiple but require thread_serialized. Abort for insufficient
!   thread support. Do not want to use MPI_COMM_WORLD everywhere, so
!   duplicate to comm_snap. Start the timer.
!_______________________________________________________________________

#ifdef SHM
    CALL MPI_INIT_THREAD ( thread_single, thread_level, ierr )
#else
    CALL MPI_INIT_THREAD ( thread_multiple, thread_level, ierr )

#endif

    CALL wtime ( t1 )

    CALL MPI_COMM_DUP ( MPI_COMM_WORLD, comm_snap, ierr )
!_______________________________________________________________________
!
!   Get the communicator size and each process rank within the main
!   communicator
!_______________________________________________________________________

    CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, wnproc, ierr )
    CALL MPI_COMM_RANK ( MPI_COMM_WORLD, wiproc, ierr )

    CALL MPI_COMM_SIZE ( comm_snap, nproc, ierr )
    CALL MPI_COMM_RANK ( comm_snap, iproc, ierr )
!_______________________________________________________________________
!
!   Put a barrier for every process to reach this point.
!_______________________________________________________________________

#ifdef SHM
    CALL MPI_COMM_SPLIT_TYPE ( MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,      &
      MPI_INFO_NULL, comm_maxshm, ierr )
    CALL MPI_COMM_SIZE ( comm_maxshm, maxshm_nproc, ierr )
    CALL MPI_COMM_RANK ( comm_maxshm, maxshm_iproc, ierr )
#endif

    CALL barrier ( comm_snap )
!_______________________________________________________________________
!_______________________________________________________________________

  END SUBROUTINE pinit


  SUBROUTINE barrier ( comm )

!-----------------------------------------------------------------------
!
! MPI barrier
!
!-----------------------------------------------------------------------

    INTEGER(i_knd), INTENT(IN) :: comm
!_______________________________________________________________________
!
!   Local variables
!_______________________________________________________________________

    INTEGER(i_knd) :: ierr
!_______________________________________________________________________

    CALL MPI_BARRIER ( comm, ierr )
!_______________________________________________________________________
!_______________________________________________________________________

  END SUBROUTINE barrier


  SUBROUTINE pcomm_set

!-----------------------------------------------------------------------
!
! Setup the SDD communicator.
!
!-----------------------------------------------------------------------
!_______________________________________________________________________
!
!   Local variables
!_______________________________________________________________________

    INTEGER(i_knd) :: ierr

    INTEGER(i_knd), DIMENSION(2) :: dims

    LOGICAL(l_knd) :: reorder

    LOGICAL(l_knd), DIMENSION(2) :: periodic, remain
!_______________________________________________________________________
!
!   Use MPI routines to create a Cartesian topology
!_______________________________________________________________________

    dims(1) = npez
    dims(2) = npey

    periodic = .FALSE.
    reorder = .TRUE.

#ifdef SHM
    ! Create topology based on master communicator
    CALL MPI_CART_CREATE ( comm_shmgrp, 2, dims, periodic, reorder,      &
      comm_space, ierr )
#else
    CALL MPI_CART_CREATE ( comm_snap, 2, dims, periodic, reorder,      &
      comm_space, ierr )
#endif
    CALL MPI_COMM_RANK ( comm_space, sproc, ierr )
    CALL MPI_COMM_SIZE ( comm_space, snproc, ierr )

!_______________________________________________________________________
!
!   Set up the sub-communicators of the cartesian grid
!_______________________________________________________________________

    remain(1) = .FALSE.
    remain(2) = .TRUE.
    CALL MPI_CART_SUB ( comm_space, remain, ycomm, ierr )
    CALL MPI_COMM_RANK ( ycomm, yproc, ierr )
    CALL MPI_COMM_SIZE ( ycomm, ynproc, ierr )

    remain(1) = .TRUE.
    remain(2) = .FALSE.
    CALL MPI_CART_SUB ( comm_space, remain, zcomm, ierr )
    CALL MPI_COMM_RANK ( zcomm, zproc, ierr )
    CALL MPI_COMM_SIZE ( zcomm, znproc, ierr )

    write(*,*) 'iproc=', iproc , '/', nproc,             &
               'ycomm iproc=', yproc , '/', ynproc,      &
               'comm_space iproc=', sproc , '/', snproc, &
               'zcomm iproc=', zproc , '/', znproc

!_______________________________________________________________________
!
!   Set some variables used during solution
!_______________________________________________________________________

    IF ( yproc == 0 ) THEN
      firsty = .TRUE.
      ylop = yproc
    ELSE
      firsty = .FALSE.
      ylop = yproc - 1
    END IF

    IF ( yproc == npey-1 ) THEN
      lasty = .TRUE.
      yhip = yproc
    ELSE
      lasty = .FALSE.
      yhip = yproc + 1
    END IF

    IF ( zproc == 0 ) THEN
      firstz = .TRUE.
      zlop = zproc
    ELSE
      firstz = .FALSE.
      zlop = zproc - 1
    END IF

    IF ( zproc == npez-1 ) THEN
      lastz = .TRUE.
      zhip = zproc
    ELSE
      lastz = .FALSE.
      zhip = zproc + 1
    END IF


!_______________________________________________________________________
!_______________________________________________________________________

  END SUBROUTINE pcomm_set


  SUBROUTINE pend

!-----------------------------------------------------------------------
!
! Call to end MPI processes
!
!-----------------------------------------------------------------------
!_______________________________________________________________________
!
!   Local variables
!_______________________________________________________________________

    INTEGER(i_knd) :: ierr
!_______________________________________________________________________
#ifdef SHM
    CALL plib_shm_destroy
    IF (comm_shm /= MPI_COMM_NULL) THEN
        CALL MPI_COMM_FREE(comm_shm, ierr)
!        write(*,*) 'comm_shm free'
    END IF
    IF (comm_maxshm /= MPI_COMM_NULL) THEN
        CALL MPI_COMM_FREE(comm_maxshm, ierr)
!        write(*,*) 'comm_maxshm free'
    END IF
    IF (comm_shmgrp /= MPI_COMM_NULL) THEN
        CALL MPI_COMM_FREE(comm_shmgrp, ierr)
!        write(*,*) 'comm_shmgrp free'
    END IF
#endif

    IF (zcomm /= MPI_COMM_NULL) THEN
        CALL MPI_COMM_FREE(zcomm, ierr)
!        write(*,*) 'zcomm free'
    END IF
    IF (ycomm /= MPI_COMM_NULL) THEN
        CALL MPI_COMM_FREE(ycomm, ierr)
!        write(*,*) 'ycomm free'
    END IF
    IF (comm_space /= MPI_COMM_NULL) THEN
        CALL MPI_COMM_FREE(comm_space, ierr)
!        write(*,*) 'comm_space free'
    END IF
    IF (comm_snap /= MPI_COMM_NULL) THEN
        CALL MPI_COMM_FREE(comm_snap, ierr)
!        write(*,*) 'comm_snap free'
    END IF

    CALL MPI_FINALIZE ( ierr )
!_______________________________________________________________________
!_______________________________________________________________________

  END SUBROUTINE pend


  SUBROUTINE glmax_i ( value, comm )

!-----------------------------------------------------------------------
!
! All reduce global max value (integer). Use specified communicator.
!
!-----------------------------------------------------------------------

    INTEGER(i_knd), INTENT(IN) :: comm

    INTEGER(i_knd), INTENT(INOUT) :: value
!_______________________________________________________________________
!
!   Local variables
!_______________________________________________________________________

    INTEGER(i_knd) :: ierr, x
!_______________________________________________________________________

    IF ( comm == MPI_COMM_NULL ) RETURN
    CALL MPI_ALLREDUCE ( value, x, 1, MPI_INTEGER, MPI_MAX, comm, ierr )
    value = x
!_______________________________________________________________________
!_______________________________________________________________________
  END SUBROUTINE glmax_i


  SUBROUTINE glmax_d ( value, comm )

!-----------------------------------------------------------------------
!
! All reduce global max value (double precision float). Use specified
! communicator.
!
!-----------------------------------------------------------------------

    INTEGER(i_knd), INTENT(IN) :: comm

    REAL(r_knd), INTENT(INOUT) :: value
!_______________________________________________________________________
!
!   Local variables
!_______________________________________________________________________

    INTEGER(i_knd) :: ierr

    REAL(r_knd) :: x
!_______________________________________________________________________

    IF ( comm == MPI_COMM_NULL ) RETURN
    CALL MPI_ALLREDUCE ( value, x, 1, MPI_DOUBLE_PRECISION, MPI_MAX,   &
      comm, ierr )
    value = x

!_______________________________________________________________________
!_______________________________________________________________________

  END SUBROUTINE glmax_d


  SUBROUTINE glmax_d_1d ( value, dlen, comm )

!-----------------------------------------------------------------------
!
! All reduce global max value (double precision float) for 1-d array.
! Use specified communicator.
!
!-----------------------------------------------------------------------

    INTEGER(i_knd), INTENT(IN) :: dlen, comm

    REAL(r_knd), DIMENSION(dlen), INTENT(INOUT) :: value
!_______________________________________________________________________
!
!   Local variables
!_______________________________________________________________________

    INTEGER(i_knd) :: ierr

    REAL(r_knd), DIMENSION(dlen) :: x
!_______________________________________________________________________

    IF ( comm == MPI_COMM_NULL ) RETURN
    CALL MPI_ALLREDUCE ( value, x, dlen, MPI_DOUBLE_PRECISION, MPI_MAX,&
      comm, ierr )
    value = x
!_______________________________________________________________________
!_______________________________________________________________________

  END SUBROUTINE glmax_d_1d


  SUBROUTINE glmin_i ( value, comm )

!-----------------------------------------------------------------------
!
! All reduce global min value (integer). Use specified communicator.
!
!-----------------------------------------------------------------------

    INTEGER(i_knd), INTENT(IN) :: comm

    INTEGER(i_knd), INTENT(INOUT) :: value
!_______________________________________________________________________
!
!   Local variables
!_______________________________________________________________________

    INTEGER(i_knd) :: ierr, x
!_______________________________________________________________________

    IF ( comm == MPI_COMM_NULL ) RETURN
    CALL MPI_ALLREDUCE ( value, x, 1, MPI_INTEGER, MPI_MIN, comm, ierr )
    value = x
!_______________________________________________________________________
!_______________________________________________________________________
  END SUBROUTINE glmin_i


  SUBROUTINE glmin_d ( value, comm )

!-----------------------------------------------------------------------
!
! All reduce global min value (double precision float). Use specified
! communicator.
!
!-----------------------------------------------------------------------

    INTEGER(i_knd), INTENT(IN) :: comm

    REAL(r_knd), INTENT(INOUT) :: value
!_______________________________________________________________________
!
!   Local variables
!_______________________________________________________________________

    INTEGER(i_knd) :: ierr

    REAL(r_knd) :: x
!_______________________________________________________________________

    IF ( comm == MPI_COMM_NULL ) RETURN
    CALL MPI_ALLREDUCE ( value, x, 1, MPI_DOUBLE_PRECISION, MPI_MIN,   &
      comm, ierr )
    value = x
!_______________________________________________________________________
!_______________________________________________________________________
  END SUBROUTINE glmin_d


  SUBROUTINE glsum_d ( value, comm )

!-----------------------------------------------------------------------
!
! All reduce global sum value (double precision float). Use specified
! communicator.
!
!-----------------------------------------------------------------------

    INTEGER(i_knd), INTENT(IN) :: comm

    REAL(r_knd), INTENT(INOUT) :: value
!_______________________________________________________________________
!
!   Local variables
!_______________________________________________________________________

    INTEGER(i_knd) :: ierr

    REAL(r_knd) :: x
!_______________________________________________________________________

    IF ( comm == MPI_COMM_NULL ) RETURN
    CALL MPI_ALLREDUCE ( value, x, 1, MPI_DOUBLE_PRECISION, MPI_SUM,   &
      comm, ierr )
    value = x
!_______________________________________________________________________
!_______________________________________________________________________
  END SUBROUTINE glsum_d


  SUBROUTINE rtsum_d_1d ( value, dlen, comm, rtproc )

!-----------------------------------------------------------------------
!
! Normal reduce-to-root sum (double precision float) for 1-d array. Use
! specified communicator.
!
!-----------------------------------------------------------------------

    INTEGER(i_knd), INTENT(IN) :: dlen, comm, rtproc

    REAL(r_knd), DIMENSION(dlen), INTENT(INOUT) :: value
!_______________________________________________________________________
!
!   Local variables
!_______________________________________________________________________

    INTEGER(i_knd) :: ierr

    REAL(r_knd), DIMENSION(dlen) :: x
!_______________________________________________________________________

    IF ( comm == MPI_COMM_NULL ) RETURN
    CALL MPI_REDUCE ( value, x, dlen, MPI_DOUBLE_PRECISION, MPI_SUM,   &
      rtproc, comm, ierr )
    value = x
!_______________________________________________________________________
  END SUBROUTINE rtsum_d_1d


  SUBROUTINE bcast_i_scalar ( value, comm, bproc )

!-----------------------------------------------------------------------
!
! Broadcast (integer scalar). Use specified communicator and casting
! proc.
!
!-----------------------------------------------------------------------

    INTEGER(i_knd), INTENT(IN) :: comm, bproc

    INTEGER(i_knd), INTENT(INOUT) :: value
!_______________________________________________________________________
!
!   Local variables
!_______________________________________________________________________

    INTEGER(i_knd) :: ierr
!_______________________________________________________________________

    IF ( nproc == 1 ) RETURN
    IF ( comm == MPI_COMM_NULL ) RETURN

    CALL MPI_BCAST ( value, 1, MPI_INTEGER, bproc, comm, ierr )
!_______________________________________________________________________
!_______________________________________________________________________
  END SUBROUTINE bcast_i_scalar


  SUBROUTINE bcast_i_1d ( value, ilen, comm, bproc )

!-----------------------------------------------------------------------
!
! Broadcast (integer 1d array). Use specified communicator and casting
! proc.
!
!-----------------------------------------------------------------------

    INTEGER(i_knd), INTENT(IN) :: ilen, comm, bproc

    INTEGER(i_knd), DIMENSION(ilen), INTENT(INOUT) :: value
!_______________________________________________________________________
!
!   Local variables
!_______________________________________________________________________

    INTEGER(i_knd) :: ierr
!_______________________________________________________________________

    IF ( nproc == 1 ) RETURN
    IF ( comm == MPI_COMM_NULL ) RETURN

    CALL MPI_BCAST ( value, ilen, MPI_INTEGER, bproc, comm,   &
      ierr )
!_______________________________________________________________________
!_______________________________________________________________________
  END SUBROUTINE bcast_i_1d


  SUBROUTINE bcast_d_scalar ( value, comm, bproc )

!-----------------------------------------------------------------------
!
! Broadcast (double precision float scalar). Use specified communicator
! and casting proc.
!
!-----------------------------------------------------------------------

    INTEGER(i_knd), INTENT(IN) :: comm, bproc

    REAL(r_knd), INTENT(INOUT) :: value
!_______________________________________________________________________
!
!   Local variables
!_______________________________________________________________________

    INTEGER(i_knd) :: ierr
!_______________________________________________________________________

    IF ( nproc == 1 ) RETURN
    IF ( comm == MPI_COMM_NULL ) RETURN

    CALL MPI_BCAST ( value, 1, MPI_DOUBLE_PRECISION, bproc, comm, ierr )
!_______________________________________________________________________
!_______________________________________________________________________
  END SUBROUTINE bcast_d_scalar


  SUBROUTINE bcast_d_1d ( value, dlen, comm, bproc )

!-----------------------------------------------------------------------
!
! Broadcast (double precision float 1d array). Use specified
! communicator and casting proc.
!
!-----------------------------------------------------------------------

    INTEGER(i_knd), INTENT(IN) :: dlen, comm, bproc

    REAL(r_knd), DIMENSION(dlen), INTENT(INOUT) :: value
!_______________________________________________________________________
!
!   Local variables
!_______________________________________________________________________

    INTEGER(i_knd) :: ierr
!_______________________________________________________________________

    IF ( nproc == 1 ) RETURN
    IF ( comm == MPI_COMM_NULL ) RETURN

    CALL MPI_BCAST ( value, dlen, MPI_DOUBLE_PRECISION, bproc, comm,   &
      ierr )
!_______________________________________________________________________
!_______________________________________________________________________
  END SUBROUTINE bcast_d_1d

  SUBROUTINE bcast_l_scalar ( value, comm, bproc )

!-----------------------------------------------------------------------
!
! Broadcast (logical scalar). Use specified communicator
! and casting proc.
!
!-----------------------------------------------------------------------

    INTEGER(i_knd), INTENT(IN) :: comm, bproc

    LOGICAL(l_knd), INTENT(INOUT) :: value
!_______________________________________________________________________
!
!   Local variables
!_______________________________________________________________________

    INTEGER(i_knd) :: ierr
!_______________________________________________________________________

    IF ( nproc == 1 ) RETURN
    IF ( comm == MPI_COMM_NULL ) RETURN

    CALL MPI_BCAST ( value, 1, MPI_LOGICAL, bproc, comm, ierr )
!_______________________________________________________________________
!_______________________________________________________________________
  END SUBROUTINE bcast_l_scalar

  SUBROUTINE psend_d_2d ( proc, myproc, d1, d2, value, comm, mtag )

!-----------------------------------------------------------------------
!
! Send a rank-3 double precision array.
!
!-----------------------------------------------------------------------

    INTEGER(i_knd), INTENT(IN) :: proc, myproc, d1, d2, comm, mtag

    REAL(r_knd), DIMENSION(d1,d2), INTENT(IN) :: value
!_______________________________________________________________________
!
!   Local variables
!_______________________________________________________________________

    INTEGER(i_knd) :: dlen, ierr
!_______________________________________________________________________

    IF ( proc==myproc .OR. comm==MPI_COMM_NULL ) RETURN

    dlen = d1*d2

    CALL MPI_SEND ( value, dlen, MPI_DOUBLE_PRECISION, proc, mtag,     &
      comm, ierr )
!_______________________________________________________________________
!_______________________________________________________________________

  END SUBROUTINE psend_d_2d


  SUBROUTINE psend_d_3d ( proc, myproc, d1, d2, d3, value, comm, mtag )

!-----------------------------------------------------------------------
!
! Send a rank-3 double precision array.
!
!-----------------------------------------------------------------------

    INTEGER(i_knd), INTENT(IN) :: proc, myproc, d1, d2, d3, comm, mtag

    REAL(r_knd), DIMENSION(d1,d2,d3), INTENT(IN) :: value
!_______________________________________________________________________
!
!   Local variables
!_______________________________________________________________________

    INTEGER(i_knd) :: dlen, ierr
!_______________________________________________________________________

    IF ( proc==myproc .OR. comm==MPI_COMM_NULL ) RETURN

    dlen = d1*d2*d3

    CALL MPI_SEND ( value, dlen, MPI_DOUBLE_PRECISION, proc, mtag,     &
      comm, ierr )
!_______________________________________________________________________
!_______________________________________________________________________

  END SUBROUTINE psend_d_3d


  SUBROUTINE isend_d_3d ( proc, myproc, d1, d2, d3, value, comm, mtag, &
    req )

!-----------------------------------------------------------------------
!
! Non-blocking-send a rank-3 double precision array.
!
!-----------------------------------------------------------------------

    INTEGER(i_knd), INTENT(IN) :: proc, myproc, d1, d2, d3, comm, mtag

    INTEGER(i_knd), INTENT(INOUT) :: req

    REAL(r_knd), DIMENSION(d1,d2,d3), INTENT(IN) :: value
!_______________________________________________________________________
!
!   Local variables
!_______________________________________________________________________

    INTEGER(i_knd) :: dlen, ierr
!_______________________________________________________________________

    IF ( proc==myproc .OR. comm==MPI_COMM_NULL ) RETURN

    dlen = d1*d2*d3

    CALL MPI_ISEND ( value, dlen, MPI_DOUBLE_PRECISION, proc, mtag,     &
      comm, req, ierr )
!_______________________________________________________________________
!_______________________________________________________________________

  END SUBROUTINE isend_d_3d


  SUBROUTINE precv_d_2d ( proc, myproc, d1, d2, value, comm, mtag )

!-----------------------------------------------------------------------
!
! Send a rank-3 double precision array.
!
!-----------------------------------------------------------------------

    INTEGER(i_knd), INTENT(IN) :: proc, myproc, d1, d2, comm, mtag

    REAL(r_knd), DIMENSION(d1,d2), INTENT(INOUT) :: value
!_______________________________________________________________________
!
!   Local variables
!_______________________________________________________________________

    INTEGER(i_knd) :: dlen, ierr, istat(MPI_STATUS_SIZE)
!_______________________________________________________________________

    IF ( proc==myproc .OR. comm==MPI_COMM_NULL ) RETURN

    dlen = d1*d2

    CALL MPI_RECV ( value, dlen, MPI_DOUBLE_PRECISION, proc, mtag,     &
      comm, istat, ierr )
!_______________________________________________________________________
!_______________________________________________________________________

  END SUBROUTINE precv_d_2d


  SUBROUTINE precv_d_3d ( proc, myproc, d1, d2, d3, value, comm, mtag )

!-----------------------------------------------------------------------
!
! Send a rank-3 double precision array.
!
!-----------------------------------------------------------------------

    INTEGER(i_knd), INTENT(IN) :: proc, myproc, d1, d2, d3, comm, mtag

    REAL(r_knd), DIMENSION(d1,d2,d3), INTENT(INOUT) :: value
!_______________________________________________________________________
!
!   Local variables
!_______________________________________________________________________

    INTEGER(i_knd) :: dlen, ierr, istat(MPI_STATUS_SIZE)
!_______________________________________________________________________

    IF ( proc==myproc .OR. comm==MPI_COMM_NULL ) RETURN

    dlen = d1*d2*d3

    CALL MPI_RECV ( value, dlen, MPI_DOUBLE_PRECISION, proc, mtag,     &
      comm, istat, ierr )
!_______________________________________________________________________
!_______________________________________________________________________

  END SUBROUTINE precv_d_3d


  SUBROUTINE cartrank ( coord, rank, comm )

!-----------------------------------------------------------------------
!
! Return the rank of a proc defined by the coordinates of the Cartesian
! communicator.
!
!-----------------------------------------------------------------------

    INTEGER(i_knd), INTENT(OUT) :: rank

    INTEGER(i_knd), INTENT(IN) :: comm

    INTEGER(i_knd), DIMENSION(2), INTENT(IN) :: coord
!_______________________________________________________________________
!
!   Local variables
!_______________________________________________________________________

    INTEGER(i_knd) :: ierr
!_______________________________________________________________________

    IF ( comm == MPI_COMM_NULL ) RETURN

    CALL MPI_CART_RANK ( comm, coord, rank, ierr )
!_______________________________________________________________________
!_______________________________________________________________________

  END SUBROUTINE cartrank


  SUBROUTINE waitinit ( req, d1 )

!-----------------------------------------------------------------------
!
! Initialize asynchronous communication request array to null state.
!
!-----------------------------------------------------------------------

    INTEGER(i_knd), INTENT(IN) :: d1

    INTEGER(i_knd), DIMENSION(d1), INTENT(OUT) :: req
!_______________________________________________________________________

    req = MPI_REQUEST_NULL
!_______________________________________________________________________
!_______________________________________________________________________

  END SUBROUTINE waitinit


  SUBROUTINE waitall ( req, d1 )

!-----------------------------------------------------------------------
!
! Wait for all asynchronous communications encapsulated in the req array
! to finish.
!
!-----------------------------------------------------------------------

    INTEGER(i_knd), INTENT(IN) :: d1

    INTEGER(i_knd), DIMENSION(d1), INTENT(INOUT) :: req
!_______________________________________________________________________
!
!   Local variables
!_______________________________________________________________________

    INTEGER(i_knd) :: info

    INTEGER(i_knd), DIMENSION(MPI_STATUS_SIZE,d1) :: stat
!_______________________________________________________________________

    CALL MPI_WAITALL ( d1, req, stat, info )
!_______________________________________________________________________
!_______________________________________________________________________

  END SUBROUTINE waitall

#else

  SUBROUTINE pinit ( t1 )
    REAL(r_knd), INTENT(OUT) :: t1
    INTEGER(i_knd) :: ierr
    CALL wtime ( t1 )
    thread_single     = 0
    thread_funneled   = 0
    thread_serialized = 0
    thread_multiple   = 0
    comm_snap = 0
    nproc = 1
    iproc = 0
  END SUBROUTINE pinit


  SUBROUTINE barrier ( comm )
    INTEGER(i_knd), INTENT(IN) :: comm
  END SUBROUTINE barrier


  SUBROUTINE pcomm_set
    comm_space = 0
    ycomm = 0
    zcomm = 0
    sproc = 0
    yproc = 0
    zproc = 0
    firsty = .TRUE.
    ylop = 0
    lasty = .TRUE.
    yhip = 0
    firstz = .TRUE.
    zlop = 0
    lastz = .TRUE.
    zhip = 0
  END SUBROUTINE pcomm_set


  SUBROUTINE pend
  END SUBROUTINE pend


  SUBROUTINE glmax_i ( value, comm )
    INTEGER(i_knd), INTENT(IN) :: comm
    INTEGER(i_knd), INTENT(IN) :: value
  END SUBROUTINE glmax_i


  SUBROUTINE glmax_d ( value, comm )
    INTEGER(i_knd), INTENT(IN) :: comm
    REAL(r_knd), INTENT(IN) :: value
  END SUBROUTINE glmax_d


  SUBROUTINE glmax_d_1d ( value, dlen, comm )
    INTEGER(i_knd), INTENT(IN) :: dlen, comm
    REAL(r_knd), DIMENSION(dlen), INTENT(IN) :: value
  END SUBROUTINE glmax_d_1d


  SUBROUTINE glmin_i ( value, comm )
    INTEGER(i_knd), INTENT(IN) :: comm
    INTEGER(i_knd), INTENT(IN) :: value
  END SUBROUTINE glmin_i


  SUBROUTINE glmin_d ( value, comm )
    INTEGER(i_knd), INTENT(IN) :: comm
    REAL(r_knd), INTENT(IN) :: value
  END SUBROUTINE glmin_d


  SUBROUTINE glsum_d ( value, comm )
    INTEGER(i_knd), INTENT(IN) :: comm
    REAL(r_knd), INTENT(IN) :: value
  END SUBROUTINE glsum_d


  SUBROUTINE rtsum_d_1d ( value, dlen, comm, rtproc )
    INTEGER(i_knd), INTENT(IN) :: dlen, comm, rtproc
    REAL(r_knd), DIMENSION(dlen), INTENT(IN) :: value
  END SUBROUTINE rtsum_d_1d


  SUBROUTINE bcast_i_scalar ( value, comm, bproc )
    INTEGER(i_knd), INTENT(IN) :: comm, bproc
    INTEGER(i_knd), INTENT(IN) :: value
  END SUBROUTINE bcast_i_scalar


  SUBROUTINE bcast_i_1d ( value, ilen, comm, bproc )
    INTEGER(i_knd), INTENT(IN) :: ilen, comm, bproc
    INTEGER(i_knd), DIMENSION(ilen), INTENT(IN) :: value
  END SUBROUTINE bcast_i_1d


  SUBROUTINE bcast_d_scalar ( value, comm, bproc )
    INTEGER(i_knd), INTENT(IN) :: comm, bproc
    REAL(r_knd), INTENT(IN) :: value
  END SUBROUTINE bcast_d_scalar


  SUBROUTINE bcast_d_1d ( value, dlen, comm, bproc )
    INTEGER(i_knd), INTENT(IN) :: dlen, comm, bproc
    REAL(r_knd), DIMENSION(dlen), INTENT(IN) :: value
  END SUBROUTINE bcast_d_1d

  SUBROUTINE bcast_l_scalar ( value, comm, bproc )
    INTEGER(i_knd), INTENT(IN) :: comm, bproc
    LOGICAL(l_knd), INTENT(IN) :: value
  END SUBROUTINE bcast_l_scalar

  SUBROUTINE psend_d_2d ( proc, myproc, d1, d2, value, comm, mtag )
    INTEGER(i_knd), INTENT(IN) :: proc, myproc, d1, d2, comm, mtag
    REAL(r_knd), DIMENSION(d1,d2), INTENT(IN) :: value
  END SUBROUTINE psend_d_2d


  SUBROUTINE psend_d_3d ( proc, myproc, d1, d2, d3, value, comm, mtag )
    INTEGER(i_knd), INTENT(IN) :: proc, myproc, d1, d2, d3, comm, mtag
    REAL(r_knd), DIMENSION(d1,d2,d3), INTENT(IN) :: value
  END SUBROUTINE psend_d_3d


  SUBROUTINE isend_d_3d ( proc, myproc, d1, d2, d3, value, comm, mtag, &
    req )
    INTEGER(i_knd), INTENT(IN) :: proc, myproc, d1, d2, d3, comm, mtag
    INTEGER(i_knd), INTENT(IN) :: req
    REAL(r_knd), DIMENSION(d1,d2,d3), INTENT(IN) :: value
  END SUBROUTINE isend_d_3d


  SUBROUTINE precv_d_2d ( proc, myproc, d1, d2, value, comm, mtag )
    INTEGER(i_knd), INTENT(IN) :: proc, myproc, d1, d2, comm, mtag
    REAL(r_knd), DIMENSION(d1,d2), INTENT(IN) :: value
  END SUBROUTINE precv_d_2d


  SUBROUTINE precv_d_3d ( proc, myproc, d1, d2, d3, value, comm, mtag )
    INTEGER(i_knd), INTENT(IN) :: proc, myproc, d1, d2, d3, comm, mtag
    REAL(r_knd), DIMENSION(d1,d2,d3), INTENT(IN) :: value
  END SUBROUTINE precv_d_3d


  SUBROUTINE cartrank ( coord, rank, comm )
    INTEGER(i_knd), INTENT(OUT) :: rank
    INTEGER(i_knd), INTENT(IN) :: comm
    INTEGER(i_knd), DIMENSION(2), INTENT(IN) :: coord
    rank = 0
  END SUBROUTINE cartrank


  SUBROUTINE waitinit ( req, d1 )
    INTEGER(i_knd), INTENT(IN) :: d1
    INTEGER(i_knd), DIMENSION(d1), INTENT(IN) :: req
  END SUBROUTINE waitinit


  SUBROUTINE waitall ( req, d1 )
    INTEGER(i_knd), INTENT(IN) :: d1
    INTEGER(i_knd), DIMENSION(d1), INTENT(IN) :: req
  END SUBROUTINE waitall

#endif


  SUBROUTINE pinit_omp ( ierr, error )

!-----------------------------------------------------------------------
!
! Setup the number of OpenMP threads. Check if any proc is exceeding
! max threads. Reset and report if so.
!
!-----------------------------------------------------------------------

    CHARACTER(LEN=64), INTENT(OUT) :: error

    INTEGER(i_knd), INTENT(OUT) :: ierr
!_______________________________________________________________________
!
!   Local variables
!_______________________________________________________________________

    INTEGER(i_knd) :: max_threads
#if defined(SHM)
    INTEGER(i_knd) :: is_shm_master_int
#endif
!_______________________________________________________________________

    ierr = 0
    error = ' '

#ifdef OPENMP

    max_threads = OMP_GET_MAX_THREADS()

    IF ( nthreads > max_threads ) THEN
      ierr = 1
      nthreads = max_threads
    END IF

    CALL glmax ( ierr, comm_snap )
    IF ( ierr /= 0 ) THEN
      error = '*WARNING: PINIT_OMP: NTHREADS>MAX_THREADS; reset to MAX_THREADS'
    END IF
!_______________________________________________________________________
!
!   Setup for nested threading
!_______________________________________________________________________

    do_nested = nnested > 1
    CALL OMP_SET_NESTED ( do_nested )
!_______________________________________________________________________
!
!   Create an array of locks, one for each thread, that will be used
!   to control threaded communications
!_______________________________________________________________________

    CALL plock_omp ( 'init', nthreads )

#elif defined(SHM)
    max_threads = maxshm_nproc
    IF ( nthreads > max_threads ) THEN
      ierr = 1
      nthreads = max_threads
    END IF
    do_nested = .FALSE.

    ! Create local shm communicator (same as a thread team)
    CALL MPI_COMM_SPLIT ( comm_maxshm, maxshm_iproc/nthreads,  maxshm_iproc,  &
      comm_shm, ierr )
    CALL MPI_COMM_SIZE ( comm_shm, shm_nproc, ierr )
    CALL MPI_COMM_RANK ( comm_shm, shm_iproc, ierr )

    IF ( shm_iproc > 0 ) THEN
        is_shm_master = .FALSE.
    ELSE
        is_shm_master = .TRUE.
    END IF

    ! Create group PIP communicator, overwrite comm_snap.iproc/nproc
    ! Consist of PIPs with the same local id across nodes. Thus each
    ! communicator should be able to build consistant topology.
    CALL MPI_COMM_SPLIT ( MPI_COMM_WORLD, shm_iproc, iproc,                  &
      comm_shmgrp, ierr )
    CALL MPI_COMM_SIZE ( comm_shmgrp, nproc, ierr )
    CALL MPI_COMM_RANK ( comm_shmgrp, iproc, ierr )

    write(*,*) 'world iproc=', wiproc , '/', wnproc,   &
               'comm_shmgrp iproc=', iproc , '/', nproc,  &
               'comm_maxshm iproc=', maxshm_iproc , '/', maxshm_nproc, &
               'comm_shm iproc=', shm_iproc, '/', shm_nproc, &
               'is_shm_master', is_shm_master

#ifdef SHM
    CALL plib_shm_init(comm_shm)
#endif

#else

    max_threads = 1
    nthreads = 1
    do_nested = .FALSE.

#endif

!_______________________________________________________________________
!_______________________________________________________________________

  END SUBROUTINE pinit_omp


#ifdef OPENMP

  SUBROUTINE plock_omp ( dowhat, nlock )

!-----------------------------------------------------------------------
!
! Operate on an OpenMP lock
!
!-----------------------------------------------------------------------

    CHARACTER(LEN=*), INTENT(IN) :: dowhat

    INTEGER(i_knd), INTENT(IN) :: nlock
!_______________________________________________________________________
!
!   Local variables
!_______________________________________________________________________

    INTEGER(i_knd) :: i
!_______________________________________________________________________

    SELECT CASE ( dowhat )

      CASE ( 'init' )
        ALLOCATE( lock(nlock) )
        DO i = 1, nlock
          CALL OMP_INIT_LOCK ( lock(i) )
        END DO
        use_lock = nproc>1 .AND. nthreads>1 .AND.                      &
                   thread_level/=thread_multiple

      CASE ( 'set' )
        CALL OMP_SET_LOCK ( lock(nlock) )

      CASE ( 'unset' )
        CALL OMP_UNSET_LOCK ( lock(nlock) )

      CASE ( 'destroy' )
        DO i = 1, nlock
          CALL OMP_DESTROY_LOCK ( lock(i) )
        END DO
        DEALLOCATE( lock )

      CASE DEFAULT
        RETURN

    END SELECT
!_______________________________________________________________________
!_______________________________________________________________________

  END SUBROUTINE plock_omp


  FUNCTION thread_num ()

!-----------------------------------------------------------------------
!
! Return thread number of caller, [0, nthreads-1]. Maintains separation
! of main code and OpenMP by placing here.
!
!-----------------------------------------------------------------------

    INTEGER(i_knd) :: thread_num
!_______________________________________________________________________

    thread_num = OMP_GET_THREAD_NUM()
!_______________________________________________________________________
!_______________________________________________________________________

  END FUNCTION thread_num

#else

  SUBROUTINE plock_omp ( dowhat, nlock )
    CHARACTER(LEN=*), INTENT(IN) :: dowhat
    INTEGER(i_knd), INTENT(IN) :: nlock
  END SUBROUTINE plock_omp


  FUNCTION thread_num ()
    INTEGER(i_knd) :: thread_num
    thread_num = 0
  END FUNCTION thread_num



#endif

#ifdef SHM
  SUBROUTINE shm_barrier ()
    INTEGER(i_knd) :: ierr
    CALL plib_shm_barrier
  END SUBROUTINE shm_barrier

  SUBROUTINE shm_deallocate_l_1d(value)
    LOGICAL(l_knd), DIMENSION(:), POINTER, INTENT(IN) :: value

    ! Local variables
    TYPE(C_PTR) :: p

    IF (ASSOCIATED(value) .EQV. .TRUE.) THEN
      p = C_LOC(value(1))
      CALL plib_shm_deallocate (p)
    END IF
  END SUBROUTINE shm_deallocate_l_1d

  SUBROUTINE shm_deallocate_d_1d(value)
    REAL(r_knd), DIMENSION(:), POINTER, INTENT(IN) :: value

    ! Local variables
    TYPE(C_PTR) :: p

    IF (ASSOCIATED(value) .EQV. .TRUE.) THEN
      p = C_LOC(value(1))
      CALL plib_shm_deallocate (p)
    END IF
  END SUBROUTINE shm_deallocate_d_1d

  SUBROUTINE shm_deallocate_d_2d(value)
    REAL(r_knd), DIMENSION(:,:), POINTER, INTENT(IN) :: value

    ! Local variables
    TYPE(C_PTR) :: p

    IF (ASSOCIATED(value) .EQV. .TRUE.) THEN
      p = C_LOC(value(1,1))
      CALL plib_shm_deallocate (p)
    END IF
  END SUBROUTINE shm_deallocate_d_2d

  SUBROUTINE shm_deallocate_d_3d(value)
    REAL(r_knd), DIMENSION(:,:,:), POINTER, INTENT(IN) :: value

    ! Local variables
    TYPE(C_PTR) :: p

    IF (ASSOCIATED(value) .EQV. .TRUE.) THEN
      p = C_LOC(value(1,1,1))
      CALL plib_shm_deallocate (p)
    END IF
  END SUBROUTINE shm_deallocate_d_3d

  SUBROUTINE shm_deallocate_d_4d(value)
    REAL(r_knd), DIMENSION(:,:,:,:), POINTER, INTENT(IN) :: value

    ! Local variables
    TYPE(C_PTR) :: p

    IF (ASSOCIATED(value) .EQV. .TRUE.) THEN
      p = C_LOC(value(1,1,1,1))
      CALL plib_shm_deallocate (p)
    END IF
  END SUBROUTINE shm_deallocate_d_4d

  SUBROUTINE shm_deallocate_d_5d(value)
    REAL(r_knd), DIMENSION(:,:,:,:,:), POINTER, INTENT(IN) :: value

    ! Local variables
    TYPE(C_PTR) :: p

    IF (ASSOCIATED(value) .EQV. .TRUE.) THEN
      p = C_LOC(value(1,1,1,1,1))
      CALL plib_shm_deallocate (p)
    END IF
  END SUBROUTINE shm_deallocate_d_5d

  SUBROUTINE shm_deallocate_d_6d(value)
    REAL(r_knd), DIMENSION(:,:,:,:,:,:), POINTER, INTENT(IN) :: value

    ! Local variables
    TYPE(C_PTR) :: p

    IF (ASSOCIATED(value) .EQV. .TRUE.) THEN
      p = C_LOC(value(1,1,1,1,1,1))
      CALL plib_shm_deallocate (p)
    END IF
  END SUBROUTINE shm_deallocate_d_6d

  SUBROUTINE shm_deallocate_i_1d(value)
    INTEGER(i_knd), DIMENSION(:), POINTER, INTENT(IN) :: value

    ! Local variables
    TYPE(C_PTR) :: p

    IF (ASSOCIATED(value) .EQV. .TRUE.) THEN
      p = C_LOC(value(1))
      CALL plib_shm_deallocate (p)
    END IF
  END SUBROUTINE shm_deallocate_i_1d

  SUBROUTINE shm_deallocate_i_2d(value)
    INTEGER(i_knd), DIMENSION(:,:), POINTER, INTENT(IN) :: value

    ! Local variables
    TYPE(C_PTR) :: p

    IF (ASSOCIATED(value) .EQV. .TRUE.) THEN
      p = C_LOC(value(1,1))
      CALL plib_shm_deallocate (p)
    END IF
  END SUBROUTINE shm_deallocate_i_2d

  SUBROUTINE shm_deallocate_i_3d(value)
    INTEGER(i_knd), DIMENSION(:,:,:), POINTER, INTENT(IN) :: value

    ! Local variables
    TYPE(C_PTR) :: p

    IF (ASSOCIATED(value) .EQV. .TRUE.) THEN
      p = C_LOC(value(1,1,1))
      CALL plib_shm_deallocate (p)
    END IF
  END SUBROUTINE shm_deallocate_i_3d

  SUBROUTINE shm_allocate_cp(dlen, type_sz, cptr, str)
    INTEGER(i_knd), INTENT(IN) :: dlen
    INTEGER(C_SIZE_T), INTENT(IN) :: type_sz
    TYPE(C_PTR), INTENT(INOUT) :: cptr
    CHARACTER(*), INTENT(IN) :: str

    ! Allocate and bcast
    cptr = C_NULL_PTR
    CALL plib_shm_allocate (type_sz * dlen, cptr, str//C_NULL_CHAR)
  END SUBROUTINE shm_allocate_cp

  SUBROUTINE shm_allocate_l_1d(d1, value, str)
    INTEGER(i_knd), INTENT(IN) :: d1
    LOGICAL(l_knd), DIMENSION(:), POINTER, INTENT(INOUT) :: value
    CHARACTER(*), INTENT(IN) :: str

    ! Local variables
    TYPE(C_PTR) :: p
    LOGICAL(l_knd) :: type

    CALL shm_allocate_cp(d1, SIZEOF(type), p, str)
    CALL C_F_POINTER(p, value, [d1])
  END SUBROUTINE shm_allocate_l_1d

  SUBROUTINE shm_allocate_d_1d(d1, value, str)
    INTEGER(i_knd), INTENT(IN) :: d1
    REAL(r_knd), DIMENSION(:), POINTER, INTENT(INOUT) :: value
    CHARACTER(*), INTENT(IN) :: str

    ! Local variables
    TYPE(C_PTR) :: p
    REAL(r_knd) :: type

    CALL shm_allocate_cp(d1, SIZEOF(type), p, str)
    CALL C_F_POINTER(p, value, [d1])
!    write(*,*) 'wiproc', wiproc, 'shm_iproc', shm_iproc, 'allocate(d,1d)', d1
  END SUBROUTINE shm_allocate_d_1d

  SUBROUTINE shm_allocate_d_2d(d1, d2, value, str)
    INTEGER(i_knd), INTENT(IN) :: d1, d2
    REAL(r_knd), DIMENSION(:,:), POINTER, INTENT(INOUT) :: value
    CHARACTER(*), INTENT(IN) :: str

    ! Local variables
    TYPE(C_PTR) :: p
    REAL(r_knd) :: type

    CALL shm_allocate_cp(d1 * d2, SIZEOF(type), p, str)
    CALL C_F_POINTER(p, value, [d1, d2])
!    write(*,*) 'wiproc', wiproc, 'shm_iproc', shm_iproc, 'allocate(d,2d)',     &
!        d1, d2
  END SUBROUTINE shm_allocate_d_2d

  SUBROUTINE shm_allocate_d_3d(d1, d2, d3, value, str)
    INTEGER(i_knd), INTENT(IN) :: d1, d2, d3
    REAL(r_knd), DIMENSION(:,:,:), POINTER, INTENT(INOUT) :: value
    CHARACTER(*), INTENT(IN) :: str

    ! Local variables
    TYPE(C_PTR) :: p
    REAL(r_knd) :: type

    CALL shm_allocate_cp(d1 * d2 * d3, SIZEOF(type), p, str)
    CALL C_F_POINTER(p, value, [d1, d2, d3])
!    write(*,*) 'wiproc', wiproc, 'shm_iproc', shm_iproc, 'allocate(d,3d)',     &
!        d1, d2, d3
  END SUBROUTINE shm_allocate_d_3d

  SUBROUTINE shm_allocate_d_4d(d1, d2, d3, d4, value, str)
    INTEGER(i_knd), INTENT(IN) :: d1, d2, d3, d4
    REAL(r_knd), DIMENSION(:,:,:,:), POINTER, INTENT(INOUT) :: value
    CHARACTER(*), INTENT(IN) :: str

    ! Local variables
    TYPE(C_PTR) :: p
    REAL(r_knd) :: type

    CALL shm_allocate_cp(d1 * d2 * d3 * d4, SIZEOF(type), p, str)
    CALL C_F_POINTER(p, value, [d1, d2, d3, d4])
!    write(*,*) 'wiproc', wiproc, 'shm_iproc', shm_iproc, 'allocate(d,4d)',     &
!        d1, d2, d3, d4
  END SUBROUTINE shm_allocate_d_4d

  SUBROUTINE shm_allocate_d_5d(d1, d2, d3, d4, d5, value, str)
    INTEGER(i_knd), INTENT(IN) :: d1, d2, d3, d4, d5
    REAL(r_knd), DIMENSION(:,:,:,:,:), POINTER, INTENT(INOUT) :: value
    CHARACTER(*), INTENT(IN) :: str

    ! Local variables
    TYPE(C_PTR) :: p
    REAL(r_knd) :: type

    CALL shm_allocate_cp(d1 * d2 * d3 * d4 * d5, SIZEOF(type), p, str)
    CALL C_F_POINTER(p, value, [d1, d2, d3, d4, d5])
!    write(*,*) 'wiproc', wiproc, 'shm_iproc', shm_iproc, 'allocate(d,5d)',     &
!        d1, d2, d3, d4, d5
  END SUBROUTINE shm_allocate_d_5d

  SUBROUTINE shm_allocate_d_6d(d1, d2, d3, d4, d5, d6, value, str)
    INTEGER(i_knd), INTENT(IN) :: d1, d2, d3, d4, d5, d6
    REAL(r_knd), DIMENSION(:,:,:,:,:,:), POINTER, INTENT(INOUT) :: value
    CHARACTER(*), INTENT(IN) :: str

    ! Local variables
    TYPE(C_PTR) :: p
    REAL(r_knd) :: type

    CALL shm_allocate_cp(d1 * d2 * d3 * d4 * d5 * d6, SIZEOF(type), p, str)
    CALL C_F_POINTER(p, value, [d1, d2, d3, d4, d5, d6])
!    write(*,*) 'wiproc', wiproc, 'shm_iproc', shm_iproc, 'allocate(d,6d)',     &
!        d1, d2, d3, d4, d5, d6
  END SUBROUTINE shm_allocate_d_6d

  SUBROUTINE shm_allocate_i_1d(d1, value, str)
    INTEGER(i_knd), INTENT(IN) :: d1
    INTEGER(i_knd), DIMENSION(:), POINTER, INTENT(INOUT) :: value
    CHARACTER(*), INTENT(IN) :: str

    ! Local variables
    TYPE(C_PTR) :: p
    INTEGER(i_knd) :: type

    CALL shm_allocate_cp(d1, SIZEOF(type), p, str)
    CALL C_F_POINTER(p, value, [d1])
!    write(*,*) 'wiproc', wiproc, 'shm_iproc', shm_iproc, 'allocate(i,1d)', d1
  END SUBROUTINE shm_allocate_i_1d

  SUBROUTINE shm_allocate_i_2d(d1, d2, value, str)
    INTEGER(i_knd), INTENT(IN) :: d1, d2
    INTEGER(i_knd), DIMENSION(:,:), POINTER, INTENT(INOUT) :: value
    CHARACTER(*), INTENT(IN) :: str

    ! Local variables
    TYPE(C_PTR) :: p
    INTEGER(i_knd) :: type

    CALL shm_allocate_cp(d1 * d2, SIZEOF(type), p, str)
    CALL C_F_POINTER(p, value, [d1, d2])
!    write(*,*) 'wiproc', wiproc, 'shm_iproc', shm_iproc, 'allocate(i,1d)', d1
  END SUBROUTINE shm_allocate_i_2d

  SUBROUTINE shm_allocate_i_3d(d1, d2, d3, value, str)
    INTEGER(i_knd), INTENT(IN) :: d1, d2, d3
    INTEGER(i_knd), DIMENSION(:,:,:), POINTER, INTENT(INOUT) :: value
    CHARACTER(*), INTENT(IN) :: str

    ! Local variables
    TYPE(C_PTR) :: p
    INTEGER(i_knd) :: type

    CALL shm_allocate_cp(d1 * d2 * d3, SIZEOF(type), p, str)
    CALL C_F_POINTER(p, value, [d1, d2, d3])
!    write(*,*) 'wiproc', wiproc, 'shm_iproc', shm_iproc, 'allocate(i,3d)',     &
!        d1, d2, d3
  END SUBROUTINE shm_allocate_i_3d
#endif


END MODULE plib_module
