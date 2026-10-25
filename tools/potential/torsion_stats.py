#!/usr/bin/env python3
"""
Statistical potential of the backbone dihedrals phi and psi from the chains of
the training set.

    torsion_stats.py training_set.txt CACHEDIR dihedrals.dat [--plot FILE]

For each class of residues (a: all but glycine; b: glycine), the distributions
of phi and psi in bins of BIN degrees, smoothed with a (1,2,1)/4 kernel on the
circle, give E = -ln P/Pmax in kT (0 at the most populated bin, at most EMAX).
Only residues whose peptide bonds with the neighbours are intact (C-N < 1.5 A).
The output is the [ dihedral_table ] chapter of montegrappa, i.e. the lines
    i  phi_a  phi_b  psi_a  psi_b
where bin i covers [-180+i*BIN, -180+(i+1)*BIN); grappino (keyword dihfile)
writes it in the potential file with the weights of the atoms: phi on CA, psi
on C, table a for every residue but glycine, b for glycine, no phi for proline.
Phi and psi are taken as independent (no correlation in the Ramachandran plane).
"""
import os, sys
import numpy as np
from contact_stats import read_chain

BIN = 5.
EMAX = 6.


def dihedral(a, b, c, d):
    b1, b2, b3 = b - a, c - b, d - c
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    x = np.dot(n1, n2)
    y = np.dot(np.cross(n1, n2), b2) / np.linalg.norm(b2)
    return np.degrees(np.arctan2(y, x))


def bonded(r1, r2):
    return "C" in r1 and "N" in r2 and np.linalg.norm(r2["N"] - r1["C"]) < 1.5


def table(angles):
    nbin = int(round(360. / BIN))
    h, _ = np.histogram(angles, bins=nbin, range=(-180., 180.))
    h = h + 1.						# pseudocount
    h = 0.25 * np.roll(h, 1) + 0.5 * h + 0.25 * np.roll(h, -1)
    e = -np.log(h / h.max())
    return np.minimum(e, EMAX)


def main(setfile, cache, out, plot=None):
    ang = {("phi", "a"): [], ("phi", "b"): [], ("psi", "a"): [], ("psi", "b"): []}
    for line in open(setfile):
        if line.startswith("#") or not line.strip():
            continue
        pdb, chain = line.split()[:2]
        fname = os.path.join(cache, pdb + ".pdb.gz")
        if not os.path.exists(fname):
            continue
        res = read_chain(fname, chain)
        for k in range(1, len(res) - 1):
            resn, at = res[k]
            prev, nxt = res[k-1][1], res[k+1][1]
            if not all(a in at for a in ("N", "CA", "C")) or not bonded(prev, at) or not bonded(at, nxt):
                continue
            c = "b" if resn == "GLY" else "a"
            ang[("phi", c)].append(dihedral(prev["C"], at["N"], at["CA"], at["C"]))
            ang[("psi", c)].append(dihedral(at["N"], at["CA"], at["C"], nxt["N"]))
    e = {k: table(np.array(v)) for k, v in ang.items()}
    nbin = len(e[("phi", "a")])
    with open(out, "w") as f:
        f.write(f"# torsion_stats.py {setfile}: {len(ang[('phi','a')])} residues (a, not glycine), {len(ang[('phi','b')])} glycines\n")
        f.write(f"# bin {BIN} degrees; E = -ln P/Pmax (kT, at most {EMAX})\n")
        f.write("# i  phi_a  phi_b  psi_a  psi_b\n")
        # one more line for dih=180 exactly, which montegrappa puts in bin nbin
        for i in range(nbin + 1):
            j = i % nbin
            f.write(f"{i:3d}\t{e[('phi','a')][j]:.4f}\t{e[('phi','b')][j]:.4f}\t{e[('psi','a')][j]:.4f}\t{e[('psi','b')][j]:.4f}\n")
    print(f"{len(ang[('phi','a')])} residues, {len(ang[('phi','b')])} glycines", file=sys.stderr)
    if plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        x = -180. + (np.arange(nbin) + 0.5) * BIN
        fig, ax = plt.subplots(figsize=(5, 3.5))
        for k, v in e.items():
            ax.plot(x, v, label=f"{k[0]} {'glycine' if k[1]=='b' else 'other'}")
        ax.set_xlabel("dihedral (deg)"); ax.set_ylabel("E / kT"); ax.legend(fontsize=7)
        fig.tight_layout(); fig.savefig(plot, dpi=120)


if __name__ == "__main__":
    args = sys.argv[1:]
    plot = None
    if "--plot" in args:
        k = args.index("--plot"); plot = args[k + 1]; del args[k:k + 2]
    if len(args) != 3:
        sys.exit(__doc__)
    main(*args, plot=plot)
