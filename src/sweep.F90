!-----------------------------------------------------------------------
!
! MODULE: sweep_module
!> @brief
!> This module controls the mesh sweep scheduling. It directs the flow
!> of KBA stages for different octants, groups, spatial work chunks.
!
!-----------------------------------------------------------------------

MODULE sweep_module

  USE global_module, ONLY: i_knd, zero, one

  USE data_module, ONLY: ng

  USE geom_module, ONLY: ndiag

  USE control_module, ONLY: inrdone, swp_typ, ncor, corner_sch

  USE octsweep_module, ONLY: octsweep

  USE solvar_module, ONLY: flkx, flky, flkz, fmin, fmax

  USE plib_module

  USE thrd_comm_module, ONLY: assign_thrd_set, sweep_wait_bdry

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: sweep

  SAVE


  CONTAINS


  SUBROUTINE sweep ( t, do_grp, ng_per_thrd, nnstd_used, grp_act )

!-----------------------------------------------------------------------
!
! Driver for the mesh sweeps. Manages the loops over octant pairs.
!
!-----------------------------------------------------------------------

    INTEGER(i_knd), INTENT(IN) :: t

    INTEGER(i_knd), INTENT(INOUT) :: ng_per_thrd, nnstd_used

    INTEGER(i_knd), DIMENSION(ng), INTENT(INOUT) :: do_grp

    INTEGER(i_knd), DIMENSION(ng,nthreads), INTENT(INOUT) :: grp_act

!_______________________________________________________________________
!
!   Local variables
!_______________________________________________________________________

    INTEGER(i_knd) :: corner, jd, kd, g, n

    INTEGER(i_knd), DIMENSION(2) :: reqs
!_______________________________________________________________________
!
!   Assign the work to threads. Main level threads always applied to
!   energy groups. Apply nested threads additionally to groups if
!   swp_typ is 0. Apply nested threads to mini-KBA if swp_typ is 1.
!_______________________________________________________________________

#ifdef SHM
  IF ( is_shm_master .EQV. .TRUE. ) THEN
#endif
  !$OMP MASTER

    do_grp = 1
    WHERE ( inrdone ) do_grp = 0

!    write(*,*) 'sweep: do_grp=', do_grp
    CALL assign_thrd_set ( do_grp, ng, ng_per_thrd, ndiag, nnstd_used, &
      grp_act )

  !$OMP END MASTER
  !$OMP BARRIER
#ifdef SHM
  END IF
  CALL shm_barrier
#endif

!_______________________________________________________________________
!
!   Each thread initializes the reqs (send request) array for
!   asynchronous sends.
!_______________________________________________________________________

    ! *PIP: thread/PIP local variable
    CALL waitinit ( reqs, SIZE( reqs ) )
!_______________________________________________________________________
!
!   Clean up and initialize some values.
!_______________________________________________________________________

    ! *PIP: each thread/PIP is taking a different t
    clean_loop: DO n = 1, ng_per_thrd

      g = grp_act(n,t)
      IF ( g == 0 ) EXIT clean_loop

      fmin(g) = HUGE( one )
      fmax(g) = zero

      flkx(:,:,:,g) = zero
      flky(:,:,:,g) = zero
      flkz(:,:,:,g) = zero

    END DO clean_loop
!_______________________________________________________________________
!
!   Loop over octant pairs; set the starting corner, i.e., the direction
!   in y and z. Spawn the nested threads if used.
!_______________________________________________________________________

  !$OMP PARALLEL NUM_THREADS(nnstd_used) IF(nnstd_used>1)              &
  !$OMP& DEFAULT(SHARED) PRIVATE(corner,jd,kd,g) PROC_BIND(CLOSE)
    corner_loop: DO corner = 1, ncor

      jd = corner_sch(1,corner)
      kd = corner_sch(2,corner)
!_______________________________________________________________________
!
!     Loop over the groups assigned to each thread.
!_______________________________________________________________________

      grp_loop: DO n = 1, ng_per_thrd
        g = grp_act(n,t)
        !write (*,*) 'sweep corner', corner, 'n', n, 't', t, 'start', jd, kd
!_______________________________________________________________________
!
!       Sweep all the chunks of an octant pair (+/- x-dir).
!_______________________________________________________________________

        CALL octsweep ( g, jd, kd, t, reqs, SIZE( reqs ) )
!_______________________________________________________________________
!
!       End the loops. Destroy the task set.
!_______________________________________________________________________

      END DO grp_loop

    END DO corner_loop
  !$OMP END PARALLEL

  !   Complete last asynchronous sends
  CALL sweep_wait_bdry ( jd, kd, t, reqs, SIZE( reqs ))
!_______________________________________________________________________
!_______________________________________________________________________

  END SUBROUTINE sweep


END MODULE sweep_module
