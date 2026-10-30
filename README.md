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

	$ make check          (about 1 min; also rebuilds montegrappa_mpi if mpicc is found)
	$ make check-asan     (memory errors, about 1.5 min)
	$ make check-stat     (detailed balance, a few minutes)

See tests/README.md.


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
 - The multiple-histogram equations (simulated tempering and mhistogram)
   used the number of samples of one temperature for all of them: the
   free energies were wrong whenever the temperatures had different
   numbers of samples. With keepall, the number of temperatures was
   overwritten by the number of stored histograms.
 - disentangle did not work: the energy of an overlap was not kept in
   the contacts, so a move removing some overlaps but not all was always
   rejected and the energy kept the overlaps after they were removed. An
   overlap now costs 1E4 (EOVERLAP), kept like any pair energy.
 - grappino counted as native the contacts at exactly imin along the
   chain, which montegrappa excludes (they never interacted).
 - Pivot moves in a periodic box could skip pairs whose minimum image
   changed with the rotation (chains longer than about half the box).

Speed:
 - Optimization of the potential: a step changes one element of the
   matrix, so only the frames that contain that pair of types change their
   weight; the sums of the reweighting are updated with them (and computed
   again from scratch every 1000 steps). The cost of a step is proportional
   to the frames of the pair instead of all the frames: for 1PGB a few
   percent of them.
 - Pivot moves do not recompute the pairs of residues that they rotate
   rigidly together (about a third of the pairs for 1PGB, 15% faster).
   A residue is taken as rigid only if its sidechain atoms, and the atoms
   its hydrogen bonds refer to, are checked to have followed the rigid
   motion; the trajectories are the same as before (up to rounding).
 - The energy of a pair of residues reads the parameters of each pair of
   atoms once and shifts the atoms of a residue to the periodic image once,
   with the same arithmetic: in all 1PGB is 30% faster than before these
   changes, with the same trajectory.
 - Pairs of residues are skipped when the spheres containing their atoms
   (centred in their mean position) are farther than the interaction
   range: most of the pairs that passed the old test, on the backbone atoms
   plus the extents of the sidechains, had no atom in range. 1PGB: 20%
   faster again, same trajectory.
 - After a move only the contact lists that it changed are copied between
   the trial and the backup copy of the system (they are recorded as they
   change), and the parameters of the pairs of types are kept in one array:
   another 6%.
 - Each entry of a contact list knows the position of the reverse entry in
   the list of its partner, so a residue is removed from the lists of its
   neighbours without searching them, and the sidechains are rebuilt from
   offsets computed once for each rotamer instead of their sines and
   cosines: another 4-5%, same trajectory.
 - The cell list keeps a compact copy of all the atoms of each residue
   (positions and types), updated with its sphere; with it the energy of a
   pair of residues computes the distances of an atom from all the atoms
   of the other residue together, and skips the parameters of the pairs of
   atoms farther than any interaction. Same trajectory, 1PGB 7% faster
   (the potentials with hydrogen bonds keep the previous function). In
   all, the Monte Carlo of 1PGB is 2 times faster than at the start of
   these changes.
 - Cell list (systems of many chains): a residue is in the cell of the
   centre of the sphere containing its atoms instead of its backbone atom,
   so the cells only need to span the radius of those spheres (about half
   the extent of the sidechains) and contain fewer candidates; the
   candidates are put in order through a bitmap instead of an insertion
   sort. Same trajectory, 16 all-atom TTR(105-115) chains in a periodic
   box 2.2 times faster.

Restart:
 - A checkpoint is written by default every 10 minutes (keyword
   tcheckpoint, in seconds); "ncheckpoint N" writes it every N steps
   instead. With MPI the default is every 100000 steps.
 - The checkpoint (checkpoint.dat, or checkpoint_proc<rank>.dat with MPI)
   contains the full state of the simulation. Restart with

	$ montegrappa checkpoint potential.pot parameters.par

   in the same directory: the run continues exactly as if it had not
   been interrupted. The checkpoint holds the length of the trajectory,
   energy, observables and exchange files, which are cut back to it, so
   nothing written after the checkpoint is repeated; the log is
   continued. Each run writes a checkpoint when it starts, so with nrun>1
   a restart never goes back to the previous run. The original polymer
   file must still be there. Simulated tempering state and frames
   recorded for the optimization of the potential are not saved.
 - kill -USR1 (soft exit) writes a checkpoint at the end of the current
   step and stops the program.

Periodic box:
 - In the [ global ] chapter of the potential file,
	boxtype p
	boxsize L
   puts the system in a periodic cubic box of side L (A), with minimum
   image interactions (boxtype c and s are boxes with walls, as before).
   L must exceed twice the range of interaction of two residues; shells
   cannot be used with it. Coordinates are not wrapped into the box.
 - Any number of chains (the limit was 25).

New keywords of the parameter file:
 - echeck          after every move, compare the energy with a full
                   recalculation and report differences (slow, to test moves)
 - cells N         cell list for the energy: 1 always, 0 never, -1 automatic
                   (default; used for large systems when it is faster)
 - tcheckpoint N   seconds between checkpoints (default 600)
 - tmginfo N       seconds between updates of the file mginfo (default 60, 0 =
                   not written), the progress of the simulation as the mdinfo
                   file of Amber: run, step, energy, acceptance and speed, or
                   step and chi2 of the optimization of the potential; time
                   elapsed and left in the run and in the whole simulation,
                   and expected end (from the duration of the runs already
                   completed; before the first one ends, Monte Carlo only).
                   With MPI it is written by rank 0
 - nprintobs N     every N steps, write observables of aggregation to obs.dat
                   (obs_run<run>_proc<rank>.dat as for energies): number of
                   clusters of chains, largest and mean cluster size, pairs of
                   chains in contact, parallel and antiparallel pairs (cosine of
                   the end-to-end vectors above 0.8), pairs in register and
                   mean register, nematic order parameter
 - op_tchi2 T      fictitious temperature of the Metropolis search of the
                   potential optimization, on the reduced chi2 (default 0:
                   only steps that lower chi2 are accepted, as before); the
                   best matrix found is returned
 - op_neff F       smallest effective fraction of frames, (sum w)^2/sum w^2
                   divided by the number of frames, allowed in the
                   reweighting of the potential optimization, relative to
                   its value with the original matrix (0 <= F < 1, default
                   0 = no limit): steps that make the sample less
                   representative are rejected
 - op_allreplicas  with parallel tempering, every replica records frames
                   (not only rank 0) and rank 0 reweights all of them to
                   op_T with MBAR (Shirts and Chodera, J. Chem. Phys. 129,
                   124105 (2008)): the weight of a frame is
                   exp(-Enew/op_T) / sum_k N_k exp(f_k - Eold/T_k), with the
                   free energies f_k of the temperatures solved
                   self-consistently (printed in the log). The hot replicas
                   bring states that the replica at op_T does not visit,
                   e.g. the unfolded ones when it stays native.
                   The contacts of the frames are now stored sparse (about
                   1200 pairs of types per frame for 1PGB instead of all
                   those met): 12 bytes per pair and frame, held on rank 0
                   for all the replicas, twice while they are reordered by
                   pair of types for the reweighting. For 1PGB, 8 replicas
                   and 50000 frames each: about 9 GB on rank 0, 13 GB in
                   all, and 36 ms per step of op_itermax
 - restraints with a negative sigma in the op_file are not fitted: their
   reduced chi2 (with |sigma|) is a validation, written as third column of
   chi2.dat and reweighted at the end of each optimization
 - record_native   now works (it was read and ignored): the starting
                   structure is the first frame of the reweighting; with
                   always_restart it is the native one at every iteration
 - grappino use_nativedist: the native distance of a pair of types is the
   largest one below r_native (it was the last one found)
 - flag keywords accept an explicit 0 (e.g. "shell 0")
 - avbmc N         aggregation-volume-bias move (Chen and Siepmann, J. Phys.
                   Chem. B 105, 11275 (2001)), frequency as the other moves:
                   a chain is moved, with a random orientation, into the
                   bonded region of another chain (centre of mass within
                   avb_rin of its centre of mass) or out of it, with the
                   acceptance that keeps detailed balance. It helps chains in
                   a dilute box find each other. Periodic box only, of side
                   larger than 2*avb_rin.
   avb_rin R       radius (A) of the bonded region (required)
   avb_pin P       probability of a move into the bonded region (default 0.5)
 - cbmc N          configurational-bias regrowth of the ends of the chains
                   (Siepmann and Frenkel, Mol. Phys. 75, 59 (1992)), frequency
                   as the other moves: the last n backbone atoms of one end
                   (n random up to nmax_cbmc) are grown again one at a time,
                   with their bond lengths and angles; the dihedral of each one
                   is chosen among ntrial_cbmc random trials with the Boltzmann
                   weight of its backbone pairs and dihedral potential, and the
                   move is accepted with the ratio of the Rosenbluth weights of
                   the new and old conformations, corrected for the rest of the
                   energy. Dihedrals that cannot rotate (move 0) are kept. It
                   rearranges the ends also in compact states, where pivot
                   moves are rejected by overlaps.
   nmax_cbmc N     largest number of backbone atoms regrown (default 6)
   ntrial_cbmc N   trials for each dihedral (default 8)
 - trjformat F     format of the trajectory: pdb (default) or dcd. The dcd
                   file <traj>.dcd is binary, single precision, with the unit
                   cell in a periodic box (read by VMD, MDAnalysis, mdtraj);
                   the names of its atoms are in <traj>_top.pdb, the initial
                   structure in the same order. "dcd2pdb traj_top.pdb traj.dcd
                   [first [last [stride]]]" converts it to pdb. A restart
                   continues the dcd file as it does the pdb one.
 - movebias without nmul_local is an error (it did nothing)
 - Moves that change the backbone rebuild all the sidechains built from what
   they moved: a flip near the end of a chain did not rebuild its last oxygen,
   whose rotamer refers to the nitrogen two atoms before, and a later move
   rebuilt it as a hidden, irreversible change. echeck now also checks that
   every sidechain matches its rotamer. Runs with sidechains differ from
   before.

Checkpoints written by earlier versions of this revision cannot be read
(more move types).

Output:
 - The coordinates in the pdb trajectories are in columns 31-54, as in the
   pdb format (they were one column to the left).
 - With MPI, each replica writes its exchanges to exchange_proc<rank>.dat
   (exchange_run<run>_proc<rank>.dat with more than one run, as the energy
   files, and a restart appends to it) and the acceptance of exchanges is
   printed at the end.
