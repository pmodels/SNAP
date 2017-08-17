/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 * (C) 2016 by Argonne National Laboratory.
 *     See COPYRIGHT in top-level directory.
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <unistd.h>

#ifdef SHM_MPIWIN

#define SHM_SIZE (2UL<<32)
static MPI_Aint shm_off = 0;
static char *shm_base_ptr = NULL;
static MPI_Win shm_win = MPI_WIN_NULL;
static int wrank = 0;
static size_t pagesize = 0;
static MPI_Comm shm_mpiwin_comm = MPI_COMM_NULL;

//#define DEBUG
#ifdef DEBUG
#define dbg_print(str,...) do {                                             \
    fprintf(stdout, "wrank %d "str, wrank, ## __VA_ARGS__); fflush(stdout);    \
    } while (0)
#else
#define dbg_print(str,...) do {} while (0)
#endif


static inline void shm_mpiwin_init(MPI_Comm comm_shm)
{
    int rank = 0, shm_rank = 0;

    shm_mpiwin_comm = comm_shm;

    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    MPI_Comm_rank(shm_mpiwin_comm, &shm_rank);

    pagesize = getpagesize();

    /* Only master PIP allocates, the children PIPs query the start address. */
    if (shm_rank == 0) {
        MPI_Win_allocate_shared(SHM_SIZE, 1, MPI_INFO_NULL, shm_mpiwin_comm,
                                &shm_base_ptr, &shm_win);
    }
    else {
        MPI_Aint r_size = 0;
        int r_disp_unit = 0;

        MPI_Win_allocate_shared(0, 1, MPI_INFO_NULL, shm_mpiwin_comm, &shm_base_ptr, &shm_win);
        MPI_Win_shared_query(shm_win, 0, &r_size, &r_disp_unit, &shm_base_ptr);
    }
    dbg_print("init shm size 0x%lx baseptr %p on comm 0x%x, win 0x%x\n",
              SHM_SIZE, shm_base_ptr, shm_mpiwin_comm, shm_win);
}

#define PLIB_PIP_ALIGN(val, align) (((size_t)(val) + (align) - 1) & ~((align) - 1))

static inline void shm_mpiwin_barrier(void)
{
    MPI_Barrier(shm_mpiwin_comm);
}

static inline void shm_mpiwin_allocate(int *size, void **ptr, const char *str)
{
    *ptr = shm_base_ptr + shm_off;
/*  shm_off += PLIB_PIP_ALIGN(*size, pagesize); */
    shm_off += *size;
    if (shm_off >= SHM_SIZE) {
        fprintf(stderr, "allocate %s shm ptr %p, size 0x%x, shm_off=0x%lx -- overflow (max 0x%lx)\n",
                str, *ptr, *size, shm_off, SHM_SIZE);
        fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    dbg_print("PLIBMPIWIN allocate %s shm ptr %p, size 0x%x, shm_off=0x%lx\n",
              str, *ptr, *size, shm_off);
}

static inline void shm_mpiwin_deallocate(void **ptr)
{
    /* do nothing */
}

static inline void shm_mpiwin_destroy(void)
{
    MPI_Win_free(&shm_win);
    dbg_print("free win 0x%x\n", shm_win);
    if (wrank == 0) {
        fprintf(stdout, "Total allocated memory on rank 0 %ld bytes %ld MB\n",
                shm_off, shm_off / 1024 / 1024);
        fflush(stdout);
    }
}
#endif /* SHM_MPIWIN */

/**
 * Public interfaces exposed to FORTRAN routines
 */
void plib_shm_init_(MPI_Comm * comm_shm_ptr)
{
#ifdef SHM_MPIWIN
    shm_mpiwin_init(*comm_shm_ptr);
#endif
}

void plib_shm_destroy_(void)
{
#ifdef SHM_MPIWIN
    shm_mpiwin_destroy();
#endif
}

void plib_shm_allocate_(int *size, void **ptr, const char *str)
{
#ifdef SHM_MPIWIN
    return shm_mpiwin_allocate(size, ptr, str);
#endif
}

void plib_shm_deallocate_(void **ptr)
{
#ifdef SHM_MPIWIN
    return shm_mpiwin_deallocate(ptr);
#endif
}

void plib_shm_barrier_(void)
{
#ifdef SHM_MPIWIN
    return shm_mpiwin_barrier();
#endif
}
