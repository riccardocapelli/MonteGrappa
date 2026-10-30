# Statistical all-atom potential

Square-well pair potential between reduced atom types, backbone hydrogen bonds
and statistical dihedrals, derived from high-resolution protein structures.
It is meant for transferable all-atom models (e.g. aggregation of amyloid
peptides), as an alternative to the structure-based Go models.

## Ingredients

- **Atom types**: `lib/atomtypes.reduced.lib`, 13 types (N, CA, C, O, H;
  aliphatic, aromatic and polar carbons; polar O and N; charged O and N; S).
- **Amide hydrogens**: grappino `hb_backbone` builds the H on the bisector of
  C(i-1)-N-CA (1.01 A from N, keyword `hb_dist`; not for prolines and first
  residues) and writes the backbone hydrogen bonds (H donors, O acceptors);
  montegrappa uses them with the keyword `hb` (angular factor on the H-O well).
- **Pairs**: DFIRE-like potentials of mean force between the types, turned
  into square wells (hard core at the 1e-4 quantile of the distances, well
  over the first attractive shell). Repulsive pairs only keep their hard core;
  the H-O well is the hydrogen bond, and its geometric correlates (N-O, C-H,
  N-H, H-H) only keep their hard cores.
- **Hydrophobic term**: the atom-level statistics hardly see the hydrophobic
  effect in the apolar pairs, so Cal, Car and S get a well of depth EPS out
  to RW (`--hp EPS RW`).
- **Dihedrals**: -ln P of phi and psi (all residues but glycine, glycine),
  grappino `dihfile` (phi on CA except prolines, psi on C), scale `e_dihtab`.

## Files

| file | content |
|---|---|
| `fetch_set.py` | builds the training set from RCSB (X-ray <= 1.6 A, Rfree < 0.22, 60-400 residues, one chain per 30% identity cluster) and downloads it |
| `training_set.txt` | the 1500 chains used (1457 could be read) |
| `contact_stats.py` | histograms of the distances between types (pairs more than `imin` 6 backbone atoms apart) |
| `stats_reduced.npz` | its output on the training set |
| `derive_potential.py` | square wells from the histograms, with the scales `--hb`, `--polar` and the hydrophobic well `--hp` |
| `torsion_stats.py` | tables of phi and psi |
| `dihedrals.dat` | its output on the training set (bins of 5 degrees) |

To rebuild from scratch (about 1 GB of downloads):

    fetch_set.py fetch training_set.txt CACHE
    contact_stats.py ../../lib/atomtypes.reduced.lib training_set.txt CACHE stats_reduced.npz
    torsion_stats.py training_set.txt CACHE dihedrals.dat

## A model

    derive_potential.py stats_reduced.npz pairs.dat --hp 0.3 5.0

and in the input of grappino:

    atomtypes       lib/atomtypes.reduced.lib
    backbone_atoms  3 N CA C
    locked_atoms    N
    imin            6
    r_hardcore      1.5
    potential       pairs
    pairsfile       pairs.dat
    hb_backbone
    dihfile         dihedrals.dat

with `hb` in the parameter file of montegrappa. A pot file that contains
types the polymer does not have is fine (those pairs are ignored).

## Calibration (in progress)

On the TTR(105-115) protofilament (16 chains, `examples/ttr105-115`), runs of
1e6 steps from the native fibril:

| EPS (`--hp`) | `--hb` | fibril at T = 1.0 | melts at |
|---|---|---|---|
| 0 (only statistical wells) | 1 | melts already at T = 0.5 | |
| 0.15 | 1, 2 | melts | |
| 0.2 | 2 | frays (largest cluster 13) | |
| 0.25 | 2 | stable | |
| 0.3 | 1 | stable | ~1.25 |
| 0.3 | 2 | stable | ~1.35 |

With EPS = 0.3 a single chain at T = 1.0 is disordered (about 20% of residues
in the alpha region, 27% in the beta region). 16 dispersed chains in a box of
120 A form amorphous clusters of 4-8 chains in 1e7 steps at T = 1.0, with
an energy much higher than that of the fibril: assembly into the fibril is
not yet observed. The dihedral term, from mostly helical globular proteins,
disfavours the fibril by about 3 kT per chain.
