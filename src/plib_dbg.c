/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 * (C) 2016 by Argonne National Laboratory.
 *     See COPYRIGHT in top-level directory.
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void plib_dbg_printf_(const char *str) {
#ifdef ENABLE_DBG
    int wrank = 0;

    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);

    fprintf(stdout, "wrank %d %s\n", wrank, str);
    fflush(stdout);
#endif
}

void plib_dbg_rootprintf_(const char *str) {
#ifdef ENABLE_ROOT_DBG
    int wrank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    if (wrank == 0) {
        fprintf(stdout, "%s\n", str);
        fflush(stdout);
    }
#endif
}

void plib_dbg_breakpoint_(void)
{
    fprintf(stdout, "reach breakpoint\n");
    fflush(stdout);
}
