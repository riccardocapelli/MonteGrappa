===============================================

                     ,/k.
                    /  ih,              
               ,-' ,  `:7b                      
             _.-/   '  /b.`.4p,
          --   ,    ,-' ^6x, `.'^=._


                 WELCOME TO
               MONTEGRAPPA 1.3 

===============================================


This file will help you to compile and install Montegrappa v1.3 .
The package contains three programs:

	- "montegrappa", the main Monte Carlo engine 
	- "grappino", a tool to create the input files for montegrappa.
	- "mhistogram", a data analysis tool
	- "mgp2pdb", a simple converter from .pol files to .pdb

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
"grappino", "mhistogram" and "mgp2pdb" will be created in the ./bin
directory.


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


================ NEW IN THIS REVISION ===================================

Results:
 - Random numbers come from xoshiro256** instead of the C library rand():
   trajectories are the same on every platform, but differ from those of
   previous versions for the same seed. Each MPI replica uses its own
   independent stream (previously all replicas used the same numbers).
 - Several moves did not satisfy detailed balance or did not update the
   energy correctly (movebias, lpivot, flip/mflip with d_cloose, backrub
   with several chains, cluster moves, H fields with several chains).
   Simulations using them should be repeated.
 - A missing Temp keyword is now an error (it silently ran at T=0).

Restart:
 - A checkpoint is written by default every 10 minutes (keyword
   tcheckpoint, in seconds); "ncheckpoint N" writes it every N steps
   instead. With MPI the default is every 100000 steps.
 - The checkpoint (checkpoint.dat, or checkpoint_proc<rank>.dat with MPI)
   contains the full state of the simulation. Restart with

	$ montegrappa checkpoint potential.pot parameters.par

   in the same directory: the run continues exactly as if it had not
   been interrupted, appending to its trajectory and energy files. The
   original polymer file must still be there. Simulated tempering state
   and frames recorded for the optimization of the potential are not
   saved.

New keywords of the parameter file:
 - echeck          after every move, compare the energy with a full
                   recalculation and report differences (slow, to test moves)
 - cells N         cell list for the energy: 1 always, 0 never, -1 automatic
                   (default; used for large systems when it is faster)
 - tcheckpoint N   seconds between checkpoints (default 600)
 - op_tchi2 T      fictitious temperature of the Metropolis search of the
                   potential optimization, on the reduced chi2 (default 0:
                   only steps that lower chi2 are accepted, as before); the
                   best matrix found is returned
 - flag keywords accept an explicit 0 (e.g. "shell 0")

Output:
 - With MPI, each replica writes its exchanges to exchange_proc<rank>.dat
   and the acceptance of exchanges is printed at the end.
