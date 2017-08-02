!-----------------------------------------------------------------------
!
! MODULE: solvar_module
!> @brief
!> This module contains several variables that are used in the solution
!> process, including their allocation and deallocation. Also includes
!> initialization of sweep parameters.
! 
!-----------------------------------------------------------------------

MODULE solvar_module

  USE global_module, ONLY: i_knd, r_knd, zero

  USE plib_module, ONLY: ichunk

  USE geom_module, ONLY: nx, ny, nz, nc

  USE sn_module, ONLY: nang, noct, nmom, cmom

  USE data_module, ONLY: ng

  USE control_module, ONLY: timedep, angcpy

#ifdef SHM
  USE plib_module
#endif

  IMPLICIT NONE

  PUBLIC

  SAVE
!_______________________________________________________________________
!
! Module variables
!
! ptr_in(nang,nx,ny,nz,noct,ng)   - Incoming time-edge flux pointer
! ptr_out(nang,nx,ny,nz,noct,ng)  - Outgoing time-edge flux pointer
!
! flux0(nx,ny,nz,ng)         - Scalar flux moments array
! flux0po(nx,ny,nz,ng)       - Previous outer copy of scalar flux array
! flux0pi(nx,ny,nz,ng)       - Previous inner copy of scalar flux array
! fluxm(cmom-1,nx,ny,nz,ng)  - Flux moments array
!
! q2grp0(nx,ny,nz,ng)        - Isotropic out-of-group + fixed sources
! q2grpm(cmom-1,nx,ny,nz,ng) - Anisotropic out-of-group + fixed sources
! qtot(cmom,ichunk,ny,nz,nc,ng) - Total source: q2grp0 + q2grpm +
!                                 within-group source
!
! t_xs(nx,ny,nz,ng)       - Total cross section on mesh
! a_xs(nx,ny,nz,ng)       - Absorption cross section on mesh
! s_xs(nx,ny,nz,nmom,ng)  - In-group scattering cross section on mesh
!
! psii(nang,ny,nz,ng)     - Working psi_x array
! psij(nang,ichunk,nz,ng) - Working psi_y array
! psik(nang,ichunk,ny,ng) - Working psi_z array
!
! jb_in(nang,ichunk,nz,ng)  - y-dir boundary flux in from comm
! jb_out(nang,ichunk,nz,ng) - y-dir boundary flux out to comm
! kb_in(nang,ichunk,ny,ng)  - z-dir boundary flux in from comm
! kb_out(nang,ichunk,ny,ng) - z-dir boundary flux out to comm
!
! flkx(nx+1,ny,nz,ng)     - x-dir leakage array
! flky(nx,ny+1,nz,ng)     - y-dir leakage array
! flkz(nx,ny,nz+1,ng)     - z-dir leakage array
!
! fmin(ng)       - dummy flux min
! fmax(ng)       - dummy flux max
!
! pop(ng)        - particle population spectrum
!
!_______________________________________________________________________

#ifdef SHM
  REAL(r_knd), DIMENSION(:), POINTER :: fmin, fmax, pop

  REAL(r_knd), DIMENSION(:,:,:,:), POINTER :: flux0, flux0po,          &
    flux0pi, q2grp0, t_xs, a_xs, psii, psij, psik, jb_in, jb_out,      &
    kb_in, kb_out, flkx, flky, flkz

  REAL(r_knd), DIMENSION(:,:,:,:,:), POINTER :: q2grpm, fluxm, s_xs

  REAL(r_knd), DIMENSION(:,:,:,:,:,:), POINTER :: qtot
#else
  REAL(r_knd), ALLOCATABLE, DIMENSION(:) :: fmin, fmax, pop

  REAL(r_knd), ALLOCATABLE, DIMENSION(:,:,:,:) :: flux0, flux0po,      &
    flux0pi, q2grp0, t_xs, a_xs, psii, psij, psik, jb_in, jb_out,      &
    kb_in, kb_out, flkx, flky, flkz

  REAL(r_knd), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: q2grpm, fluxm, s_xs

  REAL(r_knd), ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: qtot
#endif

  REAL(r_knd), DIMENSION(:,:,:,:,:,:), POINTER :: ptr_in, ptr_out

  CONTAINS


  SUBROUTINE solvar_allocate ( ierr )

!-----------------------------------------------------------------------
!
! Allocate solution arrays.
!
!-----------------------------------------------------------------------

    INTEGER(i_knd), INTENT(OUT) :: ierr
!_______________________________________________________________________
!
!   Allocate ptr_in/out if needed. If angcpy=1, only allocate for one
!   copy and point the other pointer at it. If angcpy=2, allocate for
!   two copies. Provide an initial condition of zero. This may be
!   changed in the future if necessary.
!_______________________________________________________________________

    ierr = 0

#ifdef SHM
    CALL plib_dbg_printf("solvar_allocate")
#endif

    NULLIFY( ptr_in, ptr_out )

    IF ( timedep == 1 ) THEN
      IF ( angcpy == 1 ) THEN
#ifdef SHM
        CALL shm_allocate(nang,nx,ny,nz,noct,ng, ptr_in, "ptr_in")
#else
        ALLOCATE( ptr_in(nang,nx,ny,nz,noct,ng), STAT=ierr )
#endif
      ELSE
#ifdef SHM
        CALL shm_allocate(nang,nx,ny,nz,noct,ng, ptr_in, "ptr_in")
        CALL shm_allocate(nang,nx,ny,nz,noct,ng, ptr_out, "ptr_out")
#else
        ALLOCATE( ptr_in(nang,nx,ny,nz,noct,ng),                       &
          ptr_out(nang,nx,ny,nz,noct,ng), STAT=ierr )
#endif
      END IF
    ELSE
      ALLOCATE( ptr_in(0,0,0,0,0,0), ptr_out(0,0,0,0,0,0), STAT=ierr )
    END IF
    IF ( ierr /= 0 ) RETURN

    IF ( timedep == 1 ) THEN
      ptr_in = zero
      IF ( angcpy == 1 ) THEN
        ptr_out => ptr_in
      ELSE
        ptr_out = zero
      END IF
    END IF
!_______________________________________________________________________
!
!   Allocate the flux moments arrays. Keep an old copy. If isotropic,
!   allocate fluxm as a dummy array to make passing contiguous pieces of
!   it in argument lists possible even in debug mode. There are better
!   ways to do this, but want to keep data structures simple for others
!   to change as they want easily.
!_______________________________________________________________________

#ifdef SHM
    CALL shm_allocate(nx,ny,nz,ng, flux0, "flux0")
    CALL shm_allocate(nx,ny,nz,ng, flux0po, "flux0po")
    CALL shm_allocate(nx,ny,nz,ng, flux0pi, "flux0pi")
#else
    ALLOCATE( flux0(nx,ny,nz,ng), flux0po(nx,ny,nz,ng),                &
      flux0pi(nx,ny,nz,ng), STAT=ierr )
    IF ( ierr /= 0 ) RETURN
#endif

    IF ( cmom > 1 ) THEN
#ifdef SHM
      CALL shm_allocate(cmom-1,nx,ny,nz,ng, fluxm, "fluxm")
#else
      ALLOCATE( fluxm(cmom-1,nx,ny,nz,ng), STAT=ierr )
      IF ( ierr /= 0 ) RETURN
#endif
    ELSE
#ifdef SHM
      CALL shm_allocate(1,nx,ny,nz,ng, fluxm, "fluxm")
#else
      ALLOCATE( fluxm(1,nx,ny,nz,ng), STAT=ierr )
      IF ( ierr /= 0 ) RETURN
#endif
    END IF

    flux0   = zero
    flux0po = zero
    flux0pi = zero
    fluxm   = zero
!_______________________________________________________________________
!
!   Allocate the source arrays. Do the same thing for q2grpm as was done
!   to fluxm above.
!_______________________________________________________________________
#ifdef SHM
    CALL shm_allocate(nx,ny,nz,ng, q2grp0, "q2grp0")
    CALL shm_allocate(cmom,ichunk,ny,nz,nc,ng, qtot, "qtot")
#else
    ALLOCATE( q2grp0(nx,ny,nz,ng), qtot(cmom,ichunk,ny,nz,nc,ng),      &
      STAT=ierr )
    IF ( ierr /= 0 ) RETURN
#endif

    IF ( cmom > 1 ) THEN
#ifdef SHM
      CALL shm_allocate(cmom-1,nx,ny,nz,ng, q2grpm, "q2grpm")
#else
      ALLOCATE( q2grpm(cmom-1,nx,ny,nz,ng), STAT=ierr )
      IF ( ierr /= 0 ) RETURN
#endif
    ELSE
#ifdef SHM
      CALL shm_allocate(1,nx,ny,nz,ng, q2grpm, "q2grpm")
#else
      ALLOCATE( q2grpm(0:0,nx,ny,nz,ng), STAT=ierr )
      IF ( ierr /= 0 ) RETURN
#endif
    END IF

    q2grp0 = zero
    q2grpm = zero
    qtot = zero
!_______________________________________________________________________
!
!   Allocate the cross section expanded to spatial mesh arrays
!_______________________________________________________________________

#ifdef SHM
    CALL shm_allocate(nx,ny,nz,ng, t_xs, "t_xs")
    CALL shm_allocate(nx,ny,nz,ng, a_xs, "a_xs")
    CALL shm_allocate(nx,ny,nz,nmom,ng, s_xs, "s_xs")
#else
    ALLOCATE( t_xs(nx,ny,nz,ng), a_xs(nx,ny,nz,ng),                    &
      s_xs(nx,ny,nz,nmom,ng), STAT=ierr )
    IF ( ierr /= 0 ) RETURN
#endif

    t_xs = zero
    a_xs = zero
    s_xs = zero
!_______________________________________________________________________
!
!   Working arrays
!_______________________________________________________________________

#ifdef SHM
    CALL shm_allocate(nang,ny,nz,ng, psii, "psii")
    CALL shm_allocate(nang,ichunk,nz,ng, psij, "psij")
    CALL shm_allocate(nang,ichunk,ny,ng, psik, "psik")
#else
    ALLOCATE( psii(nang,ny,nz,ng), psij(nang,ichunk,nz,ng),            &
      psik(nang,ichunk,ny,ng), STAT=ierr )
    IF ( ierr /= 0 ) RETURN
#endif

    psii = zero
    psij = zero
    psik = zero
!_______________________________________________________________________
!
!   PE boundary flux arrays
!_______________________________________________________________________

#ifdef SHM
    CALL shm_allocate(nang,ichunk,nz,ng, jb_in, "jb_in")
    CALL shm_allocate(nang,ichunk,nz,ng, jb_out, "jb_out")
    CALL shm_allocate(nang,ichunk,ny,ng, kb_in, "kb_in")
    CALL shm_allocate(nang,ichunk,ny,ng, kb_out, "kb_out")
#else
    ALLOCATE( jb_in(nang,ichunk,nz,ng), jb_out(nang,ichunk,nz,ng),     &
      kb_in(nang,ichunk,ny,ng), kb_out(nang,ichunk,ny,ng), STAT=ierr )
    IF ( ierr /= 0 ) RETURN
#endif

    jb_in  = zero
    jb_out = zero
    kb_in  = zero
    kb_out = zero
!_______________________________________________________________________
!
!   Leakage arrays
!_______________________________________________________________________

#ifdef SHM
    CALL shm_allocate(nx+1,ny,nz,ng, flkx, "flkx")
    CALL shm_allocate(nx,ny+1,nz,ng, flky, "flky")
    CALL shm_allocate(nx,ny,nz+1,ng, flkz, "flkz")
#else
    ALLOCATE( flkx(nx+1,ny,nz,ng), flky(nx,ny+1,nz,ng),                &
      flkz(nx,ny,nz+1,ng), STAT=ierr )
    IF ( ierr /= 0 ) RETURN
#endif

    flkx = zero
    flky = zero
    flkz = zero
!_______________________________________________________________________
!
!   Flux extremes and particle population spectrum
!_______________________________________________________________________

#ifdef SHM
    CALL shm_allocate(ng, fmin, "fmin")
    CALL shm_allocate(ng, fmax, "fmax")
    CALL shm_allocate(ng, pop, "pop")
#else
    ALLOCATE( fmin(ng), fmax(ng), pop(ng), STAT=ierr )
    IF ( ierr /= 0 ) RETURN
#endif

    fmin = zero
    fmax = zero
    pop = zero
!_______________________________________________________________________
!_______________________________________________________________________

#ifdef SHM
  CALL shm_barrier
#endif
  CALL plib_dbg_rootprintf("solvar_allocate done")

  END SUBROUTINE solvar_allocate


  SUBROUTINE solvar_deallocate

!-----------------------------------------------------------------------
!
! Deallocate solve_module arrays.
!
!-----------------------------------------------------------------------
!_______________________________________________________________________

#ifdef SHM
    CALL shm_deallocate(ptr_in)
    IF ( angcpy==2 .OR. timedep==0 ) CALL shm_deallocate(ptr_out)
    CALL shm_deallocate(flux0)
    CALL shm_deallocate(flux0po)
    CALL shm_deallocate(flux0pi)
    CALL shm_deallocate(fluxm)

    CALL shm_deallocate(q2grp0)
    CALL shm_deallocate(q2grpm)
    CALL shm_deallocate(qtot)

    CALL shm_deallocate(t_xs)
    CALL shm_deallocate(a_xs)
    CALL shm_deallocate(s_xs)

    CALL shm_deallocate(psii)
    CALL shm_deallocate(psij)
    CALL shm_deallocate(psik)

    CALL shm_deallocate(jb_in)
    CALL shm_deallocate(jb_out)
    CALL shm_deallocate(kb_in)
    CALL shm_deallocate(kb_out)

    CALL shm_deallocate(flkx)
    CALL shm_deallocate(flky)
    CALL shm_deallocate(flkz)

    CALL shm_deallocate(fmin)
    CALL shm_deallocate(fmax)
    CALL shm_deallocate(pop)
#else
    DEALLOCATE( ptr_in )
    IF ( angcpy==2 .OR. timedep==0 ) DEALLOCATE( ptr_out )
    DEALLOCATE( flux0, flux0po, flux0pi, fluxm )
    DEALLOCATE( q2grp0, q2grpm, qtot )
    DEALLOCATE( t_xs, a_xs, s_xs )
    DEALLOCATE( psii, psij, psik )
    DEALLOCATE( jb_in, jb_out, kb_in, kb_out )
    DEALLOCATE( flkx, flky, flkz )
    DEALLOCATE( fmin, fmax, pop )
#endif
!_______________________________________________________________________
!_______________________________________________________________________

  END SUBROUTINE solvar_deallocate


END MODULE solvar_module
