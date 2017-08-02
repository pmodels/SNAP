/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 * (C) 2016 by Argonne National Laboratory.
 *     See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#ifdef SHM_PIP

#include <pip.h>

static int wrank = 0, pip_rank = 0;
#define DEBUG
#ifdef DEBUG
#define dbg_print(str,...) do {                                             \
    fprintf(stdout, "wrank %d "str, wrank, ## __VA_ARGS__); fflush(stdout); \
    } while (0)
#else
#define dbg_print(str,...) do {} while (0)
#endif

#define fail_print(str,...) do {                                              \
        fprintf(stderr, "PLIBPIP error at %s: wrank %d "str, __FUNCTION__,    \
                wrank, ## __VA_ARGS__);                                       \
        fflush(stderr);                                                       \
        MPI_Abort(MPI_COMM_WORLD, -1);                                        \
    } while (0)

static int shm_pip_alloc_cnt = 0;
static MPI_Comm shm_pip_comm = MPI_COMM_NULL;
static pip_barrier_t pip_barrier, *pip_barrier_p = NULL;

static inline void shm_pip_init(MPI_Comm comm)
{
    int pip_size = 0;
    int err = MPI_SUCCESS;
    MPI_Aint pip_barrier_p_aint = 0;

    shm_pip_comm = comm;

    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    MPI_Comm_rank(shm_pip_comm, &pip_rank);
    MPI_Comm_rank(shm_pip_comm, &pip_size);

    /* Root initializes barrier and bcast address. */
    if (pip_rank == 0) {
        pip_barrier_init(&pip_barrier, pip_size);
        pip_barrier_p = &pip_barrier;
        pip_barrier_p_aint = (MPI_Aint) pip_barrier_p;
    }

    err = MPI_Bcast(&pip_barrier_p_aint, 1, MPI_AINT, 0, shm_pip_comm);
    if (err != MPI_SUCCESS) {
        fail_print("bcast fails\n");
    }
    pip_barrier_p = (pip_barrier_t *) pip_barrier_p_aint;

    /* Make sure no one access the pip barrier before initialization. */
    MPI_Barrier(shm_pip_comm);
}

static inline void shm_pip_destroy(void)
{
    if (shm_pip_alloc_cnt > 0) {
        fail_print("allocate_count %d > 0!\n", shm_pip_alloc_cnt);
    }
}

/* TODO: PIP way to bcast instead of MPI ? need pipid. */
static inline void shm_pip_bcast_fp(void **ptr)
{
    int err = MPI_SUCCESS;

    MPI_Aint ptr_addr = (MPI_Aint) (*ptr);

    err = MPI_Bcast(&ptr_addr, 1, MPI_AINT, 0, shm_pip_comm);
    if (err != MPI_SUCCESS) {
        fail_print("bcast fails\n");
    }

    *ptr = (void *) ptr_addr;
    dbg_print("PLIBPIP bcast fp %p on shm_pip_comm\n", *ptr);
}

static inline void shm_pip_barrier(void)
{
    pip_barrier_wait(pip_barrier_p);
}

static inline void shm_pip_allocate(int *size, void **ptr, const char *str)
{
    void *c_ptr = NULL;

    if (pip_rank == 0) {
        c_ptr = malloc(*size);
        if (c_ptr == NULL) {
            fail_print("malloc %ld bytes fails\n", *size);
        }

        dbg_print("PLIBPIP malloc %s ptr %p, size 0x%x, allocate_count=%d\n",
                  str, c_ptr, *size, shm_pip_alloc_cnt);
    }

    shm_pip_bcast_fp(&c_ptr);

    *ptr = c_ptr;
    shm_pip_alloc_cnt++;
}

static inline void shm_pip_deallocate(void **ptr)
{
    shm_pip_alloc_cnt--;

    if (pip_rank == 0) {
        dbg_print("PLIBPIP free ptr %p, allocate_count=%d\n", *ptr, shm_pip_alloc_cnt);
        free(*ptr);
    }
}
#endif

/**
 * Public interfaces exposed to FORTRAN routines
 */
void plib_shm_init_(MPI_Comm * comm_shm_ptr)
{
#ifdef SHM_PIP
    shm_pip_init(*comm_shm_ptr);
#endif
}

void plib_shm_destroy_(void)
{
#ifdef SHM_PIP
    shm_pip_destroy();
#endif
}

void plib_shm_allocate_(int *size, void **ptr, const char *str)
{
#ifdef SHM_PIP
    return shm_pip_allocate(size, ptr, str);
#endif
}

void plib_shm_deallocate_(void **ptr)
{
#ifdef SHM_PIP
    return shm_pip_deallocate(ptr);
#endif
}

void plib_shm_barrier_(void)
{
#ifdef SHM_PIP
    return shm_pip_barrier();
#endif
}
