#
# MonteGrappa Makefile
#
#   make                       serial montegrappa              -> bin/montegrappa
#   make version=STEMPERING    serial montegrappa + stempering -> bin/montegrappa (needs GSL)
#   make version=MPI           parallel tempering              -> bin/montegrappa_mpi
#   make version=ASAN          serial, AddressSanitizer+UBSan  -> bin/montegrappa_asan
#   make grappino | mgp2pdb | mhistogram
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

MC_SRC = geometry io mc local_move memory misc potential optimizepot montegrappa

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

OBJDIR = obj/$(version)
MC_OBJ = $(MC_SRC:%=$(OBJDIR)/%.o)

GRAPPINO_SRC = geometry io memory misc potential energy pdb rotamers grappino optimizepot
GRAPPINO_OBJ = $(GRAPPINO_SRC:%=obj/TOOLS/%.o)
MGP2PDB_SRC  = mgp2pdb io memory geometry misc
MGP2PDB_OBJ  = $(MGP2PDB_SRC:%=obj/TOOLS/%.o)

.PHONY: montegrappa grappino mgp2pdb mhistogram all clean cleanobj FORCE

montegrappa: $(MC_BIN)

# always relinked: SERIAL and STEMPERING both produce bin/montegrappa
$(MC_BIN): $(MC_OBJ) FORCE | $(BINDIR)
	$(CCV) $(CFLAGS) $(MC_OBJ) -o $@ $(LDLIBS)

FORCE:

$(OBJDIR)/%.o: $(SRCDIR)/%.c | $(OBJDIR)
	$(CCV) $(CPPFLAGS) $(CFLAGS) -MMD -MP -c $< -o $@

grappino: $(BINDIR)/grappino
mgp2pdb:  $(BINDIR)/mgp2pdb

$(BINDIR)/grappino: $(GRAPPINO_OBJ) | $(BINDIR)
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

$(BINDIR)/mgp2pdb: $(MGP2PDB_OBJ) | $(BINDIR)
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

obj/TOOLS/%.o: $(SRCDIR)/%.c | obj/TOOLS
	$(CC) $(CPPFLAGS) $(CFLAGS) -MMD -MP -c $< -o $@

mhistogram: | $(BINDIR)
	$(MAKE) -C src/mhistogram CC="$(CC)" GSL_CFLAGS="$(GSL_CFLAGS)" GSL_LIBS="$(GSL_LIBS)"

$(BINDIR) $(OBJDIR) obj/TOOLS:
	mkdir -p $@

all:
	$(MAKE) version=STEMPERING
	$(MAKE) version=MPI
	$(MAKE) grappino mgp2pdb mhistogram

cleanobj:
	rm -rf obj/SERIAL obj/STEMPERING obj/MPI obj/ASAN obj/TOOLS
	rm -f src/mhistogram/*.o

clean: cleanobj
	rm -f $(BINDIR)/montegrappa $(BINDIR)/montegrappa_mpi $(BINDIR)/montegrappa_asan \
	      $(BINDIR)/grappino $(BINDIR)/mgp2pdb $(BINDIR)/mhistogram

-include $(wildcard obj/*/*.d)
