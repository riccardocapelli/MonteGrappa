CC         ?= cc
GSL_PREFIX ?= $(shell gsl-config --prefix 2>/dev/null || echo /usr/local)
GSL_CFLAGS ?= -I$(GSL_PREFIX)/include
GSL_LIBS   ?= -L$(GSL_PREFIX)/lib -lgsl -lgslcblas
CFLAGS     ?= -Wall -O3 -funroll-loops -finline-functions
CPPFLAGS   += $(GSL_CFLAGS)

OBJ = mhistogram.o memory.o misc.o do_mhistogram.o io.o

../../bin/mhistogram: $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $@ $(GSL_LIBS) -lm

%.o: %.c mhistogram.h do_mhistogram.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@
