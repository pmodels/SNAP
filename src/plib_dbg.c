/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 * (C) 2016 by Argonne National Laboratory.
 *     See COPYRIGHT in top-level directory.
 */
#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

static int wrank = 0;

#define fail_print(str,...) do {                                              \
        fprintf(stderr, "PLIBDBG error at %s: wrank %d "str, __FUNCTION__,    \
                wrank, ## __VA_ARGS__);                                       \
        fflush(stderr);                                                       \
        MPI_Abort(MPI_COMM_WORLD, -1);                                        \
    } while (0)

#define info_print(str,...) do {                                              \
        fprintf(stdout, "PLIBDBG info at %s: wrank %d "str, __FUNCTION__,    \
                wrank, ## __VA_ARGS__);                                       \
        fflush(stdout);                                                       \
    } while (0)

void plib_dbg_printf_(const char *str)
{
#ifdef ENABLE_DBG
    int wrank = 0;

    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);

    fprintf(stdout, "wrank %d %s\n", wrank, str);
    fflush(stdout);
#endif
}

void plib_dbg_rootprintf_(const char *str)
{
#ifdef ENABLE_ROOT_DBG
    int wrank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    if (wrank == 0) {
        fprintf(stdout, "%s\n", str);
        fflush(stdout);
    }
#endif
}

#include <sched.h>
#include <unistd.h>
#ifdef OPENMP
#include <omp.h>
#endif

void plib_dbg_check_thread_cpubind_(void)
{
#if defined(PLIB_CHECK_CPUBIND) && defined(OPENMP)
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
#pragma omp parallel
    {
        int i;
        int nthreads, tid;
        cpu_set_t cpuset;
        int ncpu = sysconf(_SC_NPROCESSORS_ONLN);

        tid = omp_get_thread_num();
        nthreads = omp_get_num_threads();

        CPU_ZERO(&cpuset);
        /* FIXME: segfault at pthread_getaffinity_np ? */
        if (pthread_getaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset)) {
            perror("pthread_getaffinity_np");
            fail_print("pthread_getaffinity_np failed\n");
        }

        for (i = 0; i < sizeof(cpuset) * 8; i++) {
            if (CPU_ISSET(i, &cpuset)) {
                info_print("thread %d/%d is bound on core %d/%d\n", tid, nthreads, i, ncpu);
            }
        }
    }
#endif
}

void plib_dbg_check_cpubind_(void)
{
#ifdef PLIB_CHECK_CPUBIND
    cpu_set_t cpuset;
    int i;
    int ncpu = sysconf(_SC_NPROCESSORS_ONLN);

    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);

    CPU_ZERO(&cpuset);
    if (sched_getaffinity(getpid(), sizeof(cpuset), &cpuset)) {
        perror("sched_getaffinity");
        fail_print("sched_getaffinity failure\n");
    }

    for (i = 0; i < sizeof(cpuset) * 8; i++) {
        if (CPU_ISSET(i, &cpuset)) {
            info_print("I am bound on core %d/%d\n", i, ncpu);
        }
    }
#endif
}

void plib_dbg_breakpoint_(void)
{
    fprintf(stdout, "reach breakpoint\n");
    fflush(stdout);
}
