# TTR(105-115) aggregation

Structure-based models of the amyloid fibril of the transthyretin fragment
105-115 (YTIAALLSPYS), used to test the aggregation machinery of MonteGrappa
(periodic box, many chains, cluster moves, observables).

## Structure

`protofilament.pdb` is model 1 of PDB 2M5N (Fitzpatrick et al., PNAS 110,
5468, 2013): 16 chains in two sheets of 8 strands. Within a sheet strands are
parallel and in register, 4.7 A apart; the two sheets are antiparallel and face
each other at 8.3 A (steric zipper of IAALL). With `nprintobs` the native
structure gives 14 parallel pairs in register and 15 antiparallel pairs.

## Models

Both use `atomtypes go_homo`: atoms of all chains share the types of chain 0,
so native contacts (also between chains) attract any pair of chains.

    ../../bin/grappino ttr_ca.in     # CA model, r_native 7.5 A -> ttr_ca.pol/pot
    ../../bin/grappino ttr_aa.in     # all heavy atoms, rotamers, r_native 4.5 A -> ttr_aa.pol/pot

`ttr_aa.in` uses `pdb_rot`, so that the experimental structure (sidechains
included, with the CB of the PDB) is one of the rotamers and is reproduced
exactly. For a periodic box of side 120 A (16 chains, about 15 mM, the
experiments used 12 mM) add to the `[ global ]` chapter of the potential

    boxtype p
    boxsize 120

## Results so far

| model | native fibril | assembly from dispersed chains |
|---|---|---|
| CA | stable up to T ~ 2, melts between 2.5 and 3 | 16 chains aggregate into parallel in-register sheets stacked in parallel, with lower energy than the native fibril: the CA model cannot tell the antiparallel pairing of the sheets |
| all-atom | stable at T = 1.2 over 5e6 steps, melts at T = 1.4 (runs of 2e5 steps are too short to tell) | no aggregation for T >= 1.3 in 3e7 steps; at T = 1.0 the chains form clusters of 6-16 chains with parallel in-register pairs, antiparallel pairs are few and never in register. The native fibril has lower energy (pair energy -770 at T = 1.2 against -550 at T = 1.0): the all-atom model prefers the native architecture, the assembly is kinetically trapped |

`status.sh` prints the last observables of the all-atom assembly runs in
`runs/`.
