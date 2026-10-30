#
# MonteGrappa Makefile
#
#   make                       serial montegrappa              -> bin/montegrappa
#   make version=STEMPERING    serial montegrappa + stempering -> bin/montegrappa (needs GSL)
#   make version=MPI           parallel tempering              -> bin/montegrappa_mpi
#   make version=ASAN          serial, AddressSanitizer+UBSan  -> bin/montegrappa_asan
#   make grappino | mgp2pdb | dcd2pdb | mhistogram
#   make version=ASAN grappino                                 -> bin/grappino_asan
#   make all                   everything above except ASAN
#
# Each version is compiled in its own obj/<version> directory, so switching
# version never links objects built with different preprocessor flags.
# Paths can be overridden from the command line, e.g.
#   make version=MPI MPICC=/usr/lib64/openmpi/bin/mpicc
#   make version=STEMPERING GSL_PREFIX=/opt/local
#

CC       ?= cc
MPICC    ?= mpicc
OPTFLAG  ?= -O2
WARN     ?= -Wall
CFLAGS   ?= $(OPTFLAG) -g $(WARN)
LDLIBS   ?= -lm

GSL_PREFIX ?= $(shell gsl-config --prefix 2>/dev/null || echo /usr/local)
GSL_CFLAGS ?= -I$(GSL_PREFIX)/include
GSL_LIBS   ?= -L$(GSL_PREFIX)/lib -lgsl -lgslcblas

# OpenMP (experimental): make OPENMP=1 (needs a compiler supporting -fopenmp)
ifeq ($(OPENMP),1)
  CFLAGS += -fopenmp
  LDLIBS += -fopenmp
endif

SRCDIR = src
BINDIR = bin
CPPFLAGS += -Iinclude

version ?= SERIAL

MC_SRC = geometry io mc local_move memory misc potential optimizepot checkpoint cells observables trajectory mginfo montegrappa

ifeq ($(version),MPI)
  CCV       = $(MPICC)
  CPPFLAGS += -D ACTIVE_MPI
  MC_SRC   += MPIfunc
  MC_BIN    = $(BINDIR)/montegrappa_mpi
else ifeq ($(version),STEMPERING)
  CCV       = $(CC)
  CPPFLAGS += -D STEMPERING $(GSL_CFLAGS)
  LDLIBS   += $(GSL_LIBS)
  MC_SRC   += stempering memory1 adjust_st do_mhistogram
  MC_BIN    = $(BINDIR)/montegrappa
else ifeq ($(version),ASAN)
  CCV       = $(CC)
  CFLAGS    = -O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined $(WARN)
  LDLIBS   += -fsanitize=address,undefined
  MC_BIN    = $(BINDIR)/montegrappa_asan
else
  CCV       = $(CC)
  MC_BIN    = $(BINDIR)/montegrappa
endif

# no fused multiply-add: the same arithmetic on all platforms (also with CFLAGS from the environment)
override CFLAGS += -ffp-contract=off

OBJDIR = obj/$(version)
MC_OBJ = $(MC_SRC:%=$(OBJDIR)/%.o)

# the tools are built with the sanitizers too in the ASAN version (objects and names of their own)
ifeq ($(version),ASAN)
  TOOLDIR = $(OBJDIR)
  TOOLSUF = _asan
else
  TOOLDIR = obj/TOOLS
  TOOLSUF =
endif
GRAPPINO_BIN = $(BINDIR)/grappino$(TOOLSUF)
GRAPPINO_SRC = geometry io memory misc potential energy pdb rotamers grappino optimizepot cells mginfo
GRAPPINO_OBJ = $(GRAPPINO_SRC:%=$(TOOLDIR)/%.o)
MGP2PDB_BIN  = $(BINDIR)/mgp2pdb$(TOOLSUF)
MGP2PDB_SRC  = mgp2pdb io memory geometry misc cells
MGP2PDB_OBJ  = $(MGP2PDB_SRC:%=$(TOOLDIR)/%.o)
DCD2PDB_BIN  = $(BINDIR)/dcd2pdb$(TOOLSUF)

.PHONY: montegrappa grappino mgp2pdb dcd2pdb mhistogram all clean cleanobj check check-asan check-stat FORCE

montegrappa: $(MC_BIN)

# always relinked: SERIAL and STEMPERING both produce bin/montegrappa
$(MC_BIN): $(MC_OBJ) FORCE | $(BINDIR)
	$(CCV) $(CFLAGS) $(LDFLAGS) $(MC_OBJ) -o $@ $(LDLIBS)

FORCE:

$(OBJDIR)/%.o: $(SRCDIR)/%.c | $(OBJDIR)
	$(CCV) $(CPPFLAGS) $(CFLAGS) -MMD -MP -c $< -o $@

grappino: $(GRAPPINO_BIN)
mgp2pdb:  $(MGP2PDB_BIN)
dcd2pdb:  $(DCD2PDB_BIN)

$(GRAPPINO_BIN): $(GRAPPINO_OBJ) | $(BINDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@ $(LDLIBS)

$(MGP2PDB_BIN): $(MGP2PDB_OBJ) | $(BINDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@ $(LDLIBS)

$(DCD2PDB_BIN): $(TOOLDIR)/dcd2pdb.o | $(BINDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@ $(LDLIBS)

obj/TOOLS/%.o: $(SRCDIR)/%.c | obj/TOOLS
	$(CC) $(CPPFLAGS) $(CFLAGS) -MMD -MP -c $< -o $@

mhistogram: | $(BINDIR)
	$(MAKE) -C src/mhistogram CC="$(CC)" GSL_CFLAGS="$(GSL_CFLAGS)" GSL_LIBS="$(GSL_LIBS)"

$(BINDIR) $(OBJDIR) obj/TOOLS:
	mkdir -p $@

all:
	$(MAKE) version=STEMPERING
	$(MAKE) version=MPI
	$(MAKE) grappino mgp2pdb dcd2pdb mhistogram

# tests (see tests/run_tests.sh)
check:
	$(MAKE) montegrappa grappino mgp2pdb dcd2pdb
	if command -v $(MPICC) >/dev/null 2>&1; then $(MAKE) version=MPI; fi
	if command -v gsl-config >/dev/null 2>&1; then $(MAKE) version=STEMPERING MC_BIN=$(BINDIR)/montegrappa_st; fi
	tests/run_tests.sh

check-asan:
	$(MAKE) version=ASAN montegrappa grappino
	tests/run_tests.sh --asan errors regression echeck restart cells tools avbmc

check-stat:
	$(MAKE) montegrappa
	tests/run_tests.sh stat

cleanobj:
	rm -rf obj/SERIAL obj/STEMPERING obj/MPI obj/ASAN obj/TOOLS
	rm -f src/mhistogram/*.o

clean: cleanobj
	rm -f $(BINDIR)/montegrappa $(BINDIR)/montegrappa_mpi $(BINDIR)/montegrappa_asan $(BINDIR)/montegrappa_st \
	      $(BINDIR)/grappino $(BINDIR)/grappino_asan $(BINDIR)/mgp2pdb $(BINDIR)/mgp2pdb_asan \
	      $(BINDIR)/dcd2pdb $(BINDIR)/dcd2pdb_asan $(BINDIR)/mhistogram

-include $(wildcard obj/*/*.d)
