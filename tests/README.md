# MonteGrappa tests

    make check          # fast tests (about 20 s)
    make check-asan     # the same with AddressSanitizer/UBSan (about 1.5 min)
    make check-stat     # detailed balance test (a few minutes)

or directly `tests/run_tests.sh [--asan] [--keep] [group ...]`. The exit
status is the number of failed tests; failed runs are kept in a temporary
directory, whose path is printed.

## Groups

| group | what it checks |
|---|---|
| `errors` | wrong input (no `Temp`, wrong `nchains`, unknown `op_minim`, too long file name, missing checkpoint) stops with a message |
| `regression` | fixed-seed runs (CA Go model, full atom with rotamers, three chains with cluster moves, two chains with H fields and all moves, hydrogen bonds with splice, potential optimisation) against `tests/reference`, with a small numerical tolerance |
| `echeck` | every move type on five systems with the keyword `echeck`, with the automatic and the forced cell list: the incremental energy, the backup copy `oldp` and the cell list must agree with a full recalculation after every move |
| `restart` | a run restarted from its checkpoint gives the same energy file, trajectory and final structure, bit by bit |
| `cells` | cell list on and off give identical trajectories |
| `tools` | `grappino` reproduces the CA model of tutorial 1; `mgp2pdb` converts a two-chain system |
| `mpi` | `montegrappa_mpi` on one rank: run and exact restart (multi-rank runs are not tested) |
| `stat` | on a 30-residue Go fragment at T=0.6, the BGS local move gives the same mean energy as pivot moves within 3 sigma |

## References

`tests/reference` holds the outputs of the regression runs. When a change
alters trajectories on purpose (a new move, a fix of the physics), check the
change and regenerate them:

    tests/run_tests.sh --update-reference regression

The comparison tolerates differences in the last digits, since the C math
library differs between platforms; a different trajectory still fails.

## Input data (`tests/data`)

| file | system |
|---|---|
| `ca.pol`, `ca.pot` | 1YPC, CA Go model (tutorial 1) |
| `polyala.pol`, `polyala.pot`, `polyala.op` | 15-residue polyalanine, full atom with rotamers (tutorial 4) |
| `polyala_ang.pot` | the same with harmonic bond angles |
| `polyala_hb.pot` | the same with hydrogen bonds (O acceptors, CB donors) and splice |
| `three.pol`, `homo.pot` | three homopolymer chains |
| `two.pol`, `two_hf.pot` | two polyalanine chains with angles and H fields |
| `fragment.pol`, `fragment.pot` | first 30 residues of 1YPC, CA Go model |
| `chain1100.pol`, `chain1100.pot` | 1100-residue homopolymer (cell list) |
| `many.pol` (with `homo.pot`) | 40 homopolymer chains |
