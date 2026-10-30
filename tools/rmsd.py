#!/usr/bin/env python3
"""RMSD of the CA atoms of each frame of a MonteGrappa pdb trajectory from a
reference structure, after optimal superposition (Kabsch), written as an xvg
file in nm like gmx rms: rmsd.py traj.pdb ref.pdb out.xvg"""
import sys
import numpy as np


def frames(fname):
    cur = []
    for line in open(fname):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            cur.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        elif line.startswith("ENDMDL"):
            if cur:
                yield np.array(cur)
            cur = []
    if cur:
        yield np.array(cur)


def rmsd(a, b):
    a = a - a.mean(0)
    b = b - b.mean(0)
    u, s, vt = np.linalg.svd(a.T @ b)
    d = np.sign(np.linalg.det(u @ vt))
    s[-1] *= d
    e = (a * a).sum() + (b * b).sum() - 2 * s.sum()
    return np.sqrt(max(e, 0) / len(a))


ref = next(frames(sys.argv[2]))
with open(sys.argv[3], "w") as out:
    out.write('@    title "RMSD"\n@    xaxis  label "frame"\n@    yaxis  label "RMSD (nm)"\n@TYPE xy\n')
    out.write('@ subtitle "C-alpha after lsq fit to C-alpha"\n')
    for i, x in enumerate(frames(sys.argv[1])):
        out.write("%12.7f %12.7f\n" % (i - 1, rmsd(x, ref) / 10))
