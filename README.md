===============================================

                     ,/k.
                    /  ih,              
               ,-' ,  `:7b                      
             _.-/   '  /b.`.4p,
          --   ,    ,-' ^6x, `.'^=._


                 WELCOME TO
               MONTEGRAPPA 2.0 

===============================================


This file will help you to compile and install Montegrappa v2.0 .
The package contains five programs:

	- "montegrappa", the main Monte Carlo engine 
	- "grappino", a tool to create the input files for montegrappa.
	- "mhistogram", a data analysis tool
	- "mgp2pdb", a simple converter from .pol files to .pdb
	- "dcd2pdb", a converter of dcd trajectories (trjformat dcd) to .pdb

The code needs GSL libraries and the MPI environment to be correctly 
installed on your machine. This is a requirement to enable STEMPERING
and PARALLEL TEMPERING functions. A basic version of Montegrappa can be
compiled without GSL and MPI.

GSL is found with gsl-config; MPI with the mpicc wrapper in your PATH.
Both can be set on the command line, e.g.

	$ make version=STEMPERING GSL_PREFIX=/opt/local
	$ make version=MPI MPICC=/usr/lib64/openmpi/bin/mpicc

The whole distribution can be compiled typing

	$ make all

The binary files "montegrappa" (with STEMPERING), "montegrappa_mpi",
"grappino", "mhistogram", "mgp2pdb" and "dcd2pdb" will be created in
the ./bin directory.


OPENMP (experimental): compile with

	$ make OPENMP=1

(the compiler must support -fopenmp; on macOS the default clang does
not) and set the variable

	$ export OMP_NUM_THREADS=N

where N is the number of threads you want to use. Performances are
system and parameters dependent: always test your configuration before
running long simulations.

================ CUSTOM INSTALLATION ====================================

Each version is compiled in its own directory obj/<version>, so there
is no need to clean between versions; "make cleanobj" removes all of
them.

(1) MonteGrappa Single Core, no STEMPERING (no PARALLEL TEMPERING)

	$ make

An executable called "montegrappa" will be created in the ./bin directory.

(2) MonteGrappa Single Core, with STEMPERING (no PARALLEL TEMPERING)

    - SYSTEM REQUIREMENTS: GSL libraries

	$ make version=STEMPERING

As well as the (1) version, an executable called "montegrappa" will be
created in the ./bin directory.

(3) MonteGrappa Multi Core, with PARALLEL TEMPERING (no STEMPERING)

    - SYSTEM REQUIREMENTS: MPI environment

	$ make version=MPI

An executable called "montegrappa_mpi" will be created in the ./bin directory.

(4) MonteGrappa with AddressSanitizer and UndefinedBehaviorSanitizer,
    to look for memory errors (slow, for development)

	$ make version=ASAN

An executable called "montegrappa_asan" will be created in the ./bin directory.

(5) grappino, mhistogram and mgp2pdb

The tools can be compiled alone, with the commands

	$ make grappino
	$ make mhistogram
	$ make mgp2pdb
	$ make dcd2pdb


(6) Tests

	$ make check          (about 1 min; also builds montegrappa_mpi if mpicc is found,
	                       and montegrappa_st for simulated tempering if gsl-config is)
	$ make check-asan     (memory errors, about 1.5 min)
	$ make check-stat     (detailed balance, a few minutes)

See tests/README.md.


(7) Manual

The reference manual of all programs, input files and keywords is
docs/manual.tex; compile it with

	$ cd docs && latexmk -pdf manual.tex


================ NEW IN 2.0 =============================================

The changes from version 1.3 are listed in CHANGELOG.md.
