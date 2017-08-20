#
# Build SNAP
#
OPT = yes
FFLAGS =
FFLAG2 =
DEFS =
PP =
SHM_OBJ=

PIPDIR=
MPIDIR=

CFLAGS=
LDFLAGS=
PREFIX=

MPI = yes
OPENMP = yes
# If enable PIP, must set also the PIP-embeded mpif90
SHM_PIP = yes
SHM_MPIWIN = no
SHM = no

FORTRAN = mpif90
CC = icc
#FORTRAN = mpifort
#FORTRAN = ftn
#FORTRAN = mpiifort

TARGET = gsnap
#TARGET = isnap
#TARGET = ksnap

#HASWELL = yes
KNL=yes

ifeq ($(MPI),yes)
  DEFS := $(DEFS) MPI
else
  ifeq ($(TARGET),gsnap)
    FORTRAN = gfortran
  else
    FORTRAN = ifort
  endif
endif

ifeq ($(SHM_PIP),yes)
  SHM = yes
  SHM_MPIWIN = no
  CFLAGS := $(CFLAGS) -DSHM_PIP -I$(PIPDIR)/include
  LDFLAGS := $(LDFLAGS) -L$(PIPDIR)/lib -Wl,-rpath -Wl,$(PIPDIR)/lib -lpip
else
  ifeq ($(SHM_MPIWIN),yes)
    SHM = yes
    SHM_PIP = no
    CFLAGS := $(CFLAGS) -DSHM_MPIWIN
  endif
endif

# Disable openmp if SHM is enabled
ifeq ($(SHM),yes)
  OPENMP = no
  DEFS := $(DEFS) SHM
  
  ifeq ($(OPT),no)
    CFLAGS := $(CFLAGS) -O0 -g
  else
    CFLAGS := $(CFLAGS) -O3
  endif
endif

ifeq ($(OPENMP),yes)
  DEFS := $(DEFS) OPENMP
  ifeq ($(TARGET),gsnap)
    OMPFLAG = -fopenmp
  else
    # Assuming Intel 15, if 14, use -openmp
    OMPFLAG = -qopenmp
  endif
else
  OMPFLAG =
endif

PPFLAGS = $(foreach def,$(DEFS),-D$(def))

ifeq ($(TARGET),gsnap)
  FFLAGS = -O3 $(OMPFLAG)
  CFLAGS := $(CFLAGS) -O3
  FFLAG2 =
  PP = cpp
  ifeq ($(OPT),no)
    TARGET = gsnap_db
    FFLAGS = -O0 $(OMPFLAG) -g -fbounds-check -fbacktrace
    CFLAGS := $(CFLAGS) -O0 -g
    FFLAG2 =
  endif
  
  ifeq ($(KNL),yes)
	FFLAGS := $(FFLAGS) -xMIC-AVX512 -qopt-report=1
	CFLAGS := $(CFLAGS) -xMIC-AVX512 -qopt-report=1
  endif
endif

ifeq ($(TARGET),isnap)
  FFLAGS = -O3 $(OMPFLAG) -ip -align array32byte -qno-opt-dynamic-align -fno-fnalias -fp-model fast -fp-speculation fast
  CFLAGS := $(CFLAGS) -O3
  ifeq ($(HASWELL),yes)
    FFLAGS := $(FFLAGS) -xcore-avx2
  else
    ifeq ($(KNL),yes)
      FFLAGS := $(FFLAGS) -xMIC-AVX512 -qopt-report=1
      CFLAGS := $(CFLAGS) -xMIC-AVX512 -qopt-report=1
    else
      FFLAGS := $(FFLAGS) -xHost
    endif
  endif
  FFLAG2 =
  PP = fpp
  ifeq ($(OPT),no)
    TARGET = isnap_db
    FFLAGS = -O0 $(OMPFLAG) -g -check bounds -traceback -warn unused
    CFLAGS := $(CFLAGS) -O0 -g
    FFLAG2 =
  endif
endif

ifeq ($(TARGET),ksnap)
  FFLAGS = -O3 $(OMPFLAG) -xmic-avx512 -ip -align array64byte -qopt-streaming-cache-evict=0 -qno-opt-dynamic-align -fp-model fast -fp-speculation fast -fno-alias -fno-fnalias
  FFLAG2 =
  PP = fpp
  ifeq ($(OPT),no)
    TARGET = ksnap_db
    FFLAGS = -O0 $(OMPFLAG) -g -check bounds -traceback -warn unused
    FFLAG2 =
  endif
endif


# avoid ifort to create temp large array on stack, which might be bad in PiP.
FFLAGS := $(FFLAGS) -heap-arrays

.SUFFIXES:	.F90 .f90 .c .o

OBJS = global.o snap_main.o utils.o version.o plib.o geom.o sn.o \
       data.o control.o input.o setup.o dealloc.o translv.o solvar.o \
       outer.o expxs.o inner.o sweep.o octsweep.o dim1_sweep.o dim3_sweep.o \
       output.o time.o mms.o analyze.o thrd_comm.o mkba_sweep.o plib_dbg.o 

ifeq ($(SHM_PIP),yes)
  SHM_OBJ := plib_shm_pip.o
else
  ifeq ($(SHM_MPIWIN),yes)
    SHM_OBJ := plib_shm_mpiwin.o
  endif
endif
OBJS := $(OBJS) $(SHM_OBJ)

SRCS = global.f90 snap_main.f90 utils.f90 version.f90 plib.F90 geom.f90 \
       sn.f90 data.f90 control.f90 input.f90 setup.f90 dealloc.f90 \
       translv.f90 solvar.f90 outer.f90 expxs.f90 inner.f90 sweep.f90 \
       octsweep.f90 dim1_sweep.f90 dim3_sweep.f90 output.f90 time.F90 mms.f90 \
       analyze.f90 thrd_comm.f90 mkba_sweep.f90 plib_dbg.c

ifeq ($(SHM_PIP),yes)
  SRCS := $(OBJS) plib_shm_pip.c
else
  ifeq ($(SHM_MPIWIN),yes)
    SRCS := $(OBJS) plib_shm_mpiwin.c
  endif
endif

$(TARGET) :	$(OBJS)
		$(FORTRAN) $(FFLAGS) $(FFLAG2) $(LDFLAGS) -o $(PREFIX)/$@ $(OBJS)

version.o :	global.o
geom.o :	global.o plib.o $(SHM_OBJ)
sn.o :		global.o plib.o
data.o :	global.o plib.o
control.o :	global.o plib.o
time.o :	global.o
plib.o :	global.o time.o plib_dbg.o $(SHM_OBJ)
thrd_comm.o:	global.o plib.o
expxs.o :	global.o geom.o sn.o data.o
solvar.o :	global.o plib.o geom.o sn.o data.o control.o
mms.o:		global.o geom.o sn.o data.o control.o plib.o
dealloc.o :	global.o sn.o data.o mms.o geom.o solvar.o control.o
utils.o :	global.o dealloc.o plib.o control.o
analyze.o:	global.o control.o geom.o data.o solvar.o plib.o
input.o :	global.o plib.o geom.o sn.o data.o control.o utils.o time.o
setup.o :	global.o geom.o sn.o data.o control.o mms.o utils.o time.o \
		plib.o
output.o :	global.o geom.o sn.o data.o control.o utils.o solvar.o \
		mms.o time.o plib.o
snap_main.o :	global.o utils.o version.o input.o setup.o output.o time.o \
		dealloc.o plib.o control.o
dim1_sweep.o :	global.o geom.o sn.o data.o control.o plib.o
dim3_sweep.o :	global.o plib.o geom.o sn.o data.o control.o thrd_comm.o
mkba_sweep.o :	global.o plib.o geom.o sn.o data.o control.o thrd_comm.o
octsweep.o :	global.o geom.o sn.o data.o control.o solvar.o \
		thrd_comm.o dim1_sweep.o dim3_sweep.o mkba_sweep.o
sweep.o :	global.o data.o geom.o control.o octsweep.o solvar.o plib.o \
		thrd_comm.o
inner.o :	global.o geom.o sn.o data.o control.o solvar.o sweep.o \
		time.o plib.o thrd_comm.o
outer.o :	global.o geom.o sn.o data.o solvar.o control.o inner.o \
		time.o plib.o expxs.o thrd_comm.o
translv.o :	global.o plib.o geom.o sn.o data.o control.o utils.o \
		solvar.o expxs.o outer.o time.o thrd_comm.o analyze.o

plib_%.o:	plib_%.c
	$(CC) $(CFLAGS) -c $<

#
# Fortran rules
#
%.o:	%.F90
	$(PP) -P $(PPFLAGS) $< >$*.f90
	$(FORTRAN) $(FFLAGS) -c $*.f90

%.o:	%.f90
	$(FORTRAN) $(FFLAGS) -c $<

#
# Cleanup
#
clean:
	rm -f *.o *.mod *.bc

#
# Count lines of code
#
count:
	rm -f Lines
	for file in $(SRCS); do ./LineCount $$file Lines; done
	awk -f ./LineReport < Lines >> Lines
	cat Lines

#
# Link compiled files only. No recompile.
#
link:
	$(FORTRAN) $(FFLAGS) $(FFLAG2) $(LDFLAGS) -o $(PREFIX)/$(TARGET) $(OBJS)


