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

`ttr_aa.in` uses `pdb_rot` and `cb_pdb`, so that the experimental structure is
reproduced exactly (the default CB of grappino is an achiral pseudo-CB in the
backbone plane). For a periodic box of side 120 A (16 chains, about 15 mM, the
experiments used 12 mM) add to the `[ global ]` chapter of the potential

    boxtype p
    boxsize 120

## Results so far

| model | native fibril | assembly from dispersed chains |
|---|---|---|
| CA | stable up to T ~ 2, melts between 2.5 and 3 | 16 chains aggregate into parallel in-register sheets stacked in parallel, with lower energy than the native fibril: the CA model cannot tell the antiparallel pairing of the sheets |
| all-atom | stable up to T = 2.0, melts between 2.0 and 2.5 | running (T = 1.8, 2.0) |

`status.sh` prints the last observables of the all-atom assembly runs in
`runs/`.
