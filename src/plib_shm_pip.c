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

static MPI_Aint SHM_SIZE = (2UL << 29); /*512MB */
static MPI_Aint shm_off = 0;
static char *shm_base_ptr = NULL;

static int wrank = 0, pip_rank = 0;
#ifdef SHM_PIP_ALIGN
static size_t pagesize = 0;
#endif

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

#define info_print(str,...) do {                                              \
        fprintf(stdout, "PLIBPIP info at %s: wrank %d "str, __FUNCTION__,    \
                wrank, ## __VA_ARGS__);                                       \
        fflush(stdout);                                                       \
    } while (0)

static int shm_pip_alloc_cnt = 0;
static MPI_Comm shm_pip_comm = MPI_COMM_NULL;

#ifdef SHM_PIP_PBARRIER
static pip_barrier_t pip_barrier, *pip_barrier_p = NULL;
#endif

#ifdef SHM_PIP_MEMPOOL_WIN
MPI_Win shm_win = MPI_WIN_NULL;
#endif

static inline void shm_pip_readenv(void)
{
    char *val = NULL;
    val = getenv("PLIB_SHM_MEMPOOL_SIZE");
    if (val != NULL && strlen(val) > 0) {
        MPI_Aint sz = atol(val);
        if (sz > 0)
            SHM_SIZE = sz;
    }
}

static inline void shm_pip_init(MPI_Comm comm)
{
    int pip_size = 0;
    int err = MPI_SUCCESS;
    MPI_Aint pip_barrier_p_aint = 0;
    MPI_Aint shm_base_ptr_aint = 0;

    shm_pip_comm = comm;

    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    MPI_Comm_rank(shm_pip_comm, &pip_rank);
    MPI_Comm_size(shm_pip_comm, &pip_size);

#ifdef SHM_PIP_ALIGN
    pagesize = getpagesize();
#endif

#ifdef SHM_PIP_MEMPOOL
    shm_pip_readenv();

#ifdef SHM_PIP_MEMPOOL_WIN
    if (pip_rank == 0) {
        MPI_Win_allocate(SHM_SIZE, 1, MPI_INFO_NULL, shm_pip_comm, &shm_base_ptr, &shm_win);
        shm_base_ptr_aint = (MPI_Aint) shm_base_ptr;
    }
    else {
        MPI_Win_allocate(0, 1, MPI_INFO_NULL, shm_pip_comm, &shm_base_ptr, &shm_win);
    }
#else
    /* Only master PIP allocates, the children PIPs query the start address. */
    if (pip_rank == 0) {
        shm_base_ptr = malloc(SHM_SIZE);
        memset(shm_base_ptr, 0, SHM_SIZE);
        shm_base_ptr_aint = (MPI_Aint) shm_base_ptr;
    }
#endif
    err = MPI_Bcast(&shm_base_ptr_aint, 1, MPI_AINT, 0, shm_pip_comm);
    if (err != MPI_SUCCESS) {
        fail_print("bcast fails\n");
    }

    shm_base_ptr = (char *) shm_base_ptr_aint;
#endif

#ifdef SHM_PIP_PBARRIER
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
    dbg_print("PLIBPIP init pip_barrier_p=%p, pip_size=%d, shm_base_ptr=%p, size=0x%lx\n",
              pip_barrier_p, pip_size, shm_base_ptr, SHM_SIZE);
    /* Make sure no one access the pip barrier before initialization. */
    MPI_Barrier(shm_pip_comm);
#endif /* end of SHM_PIP_PBARRIER */

    if (wrank == 0) {
#ifdef SHM_PIP_ALIGN
        info_print("SHM_PIP_ALIGN enabled\n");
#endif
#ifdef SHM_PIP_PBARRIER
        info_print("SHM_PIP_PBARRIER enabled\n");
#endif
#ifdef SHM_PIP_MEMPOOL
        info_print("SHM_PIP_MEMPOOL enabled\n");
        info_print("SHM_SIZE 0x%lx (%ld MB)\n", SHM_SIZE, SHM_SIZE / 1024 / 1024);
#ifdef SHM_PIP_MEMPOOL_WIN
        info_print("SHM_PIP_MEMPOOL_WIN enabled\n");
#endif
#endif

    }
}

static inline void shm_pip_destroy(void)
{
    if (shm_pip_alloc_cnt > 0) {
        fail_print("allocate_count %d > 0!\n", shm_pip_alloc_cnt);
    }

#ifdef SHM_PIP_MEMPOOL
#ifdef SHM_PIP_MEMPOOL_WIN
    MPI_Win_free(&shm_win);
#else
    if (pip_rank == 0) {
        free(shm_base_ptr);
    }
#endif
#endif
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
//    dbg_print("PLIBPIP bcast fp %p on shm_pip_comm\n", *ptr);
}

static inline void shm_pip_barrier(void)
{
#ifdef SHM_PIP_PBARRIER
    pip_barrier_wait(pip_barrier_p);
#else
    MPI_Barrier(shm_pip_comm);
#endif
}

#define PLIB_PIP_ALIGN(val, align) (((size_t)(val) + (align) - 1) & ~((align) - 1))

static inline void shm_pip_allocate(size_t * size, void **ptr, const char *str)
{
    void *c_ptr = NULL;
#ifdef SHM_PIP_ALIGN
    size_t sz_align = PLIB_PIP_ALIGN(*size, pagesize);
#else
    size_t sz_align = *size;
#endif

#ifdef SHM_PIP_MEMPOOL
    if (shm_off + sz_align > SHM_SIZE) {
        fail_print("malloc %ld bytes fails -- overflow limit 0x%lx\n", sz_align, SHM_SIZE);
    }
    c_ptr = shm_base_ptr + shm_off;
    shm_off += sz_align;
#else
    if (pip_rank == 0) {
        c_ptr = malloc(sz_align);
        if (c_ptr == NULL) {
            fail_print("malloc %ld bytes fails\n", sz_align);
        }

        dbg_print("PLIBPIP malloc %s ptr %p, size 0x%lx(align 0x%lx), allocate_count=%d\n",
                  str, c_ptr, *size, sz_align, shm_pip_alloc_cnt);
    }

    shm_pip_bcast_fp(&c_ptr);
#endif
    *ptr = c_ptr;
    shm_pip_alloc_cnt++;
}

static inline void shm_pip_deallocate(void **ptr)
{
    shm_pip_alloc_cnt--;

#ifdef SHM_PIP_MEMPOOL
    /* Do nothing */
#else
    if (pip_rank == 0) {
        dbg_print("PLIBPIP free ptr %p, allocate_count=%d\n", *ptr, shm_pip_alloc_cnt);
        free(*ptr);
    }
#endif
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

void plib_shm_allocate_(size_t * size, void **ptr, const char *str)
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
