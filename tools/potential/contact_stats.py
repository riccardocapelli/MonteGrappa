#!/usr/bin/env python3
"""
Distance statistics between reduced atom types in the chains of the training set.

    contact_stats.py TYPES.lib training_set.txt CACHEDIR stats.npz

For each chain: heavy atoms of the first model (alternate location blank or A,
selenomethionine read as methionine), plus the amide H built as grappino does
(hb_backbone: on the bisector of C(i-1)-N-CA, 1.01 A from N; not for prolines
and first residues). For each pair of types, histogram of the distances up to
RMAX between atoms whose backbone atoms (N, CA, C; H belongs to N, O to C,
the sidechain to CA) are more than IMIN apart along the chain, the rule of
montegrappa for non-bonded pairs; and the number of atoms of each
type. Only intra-chain pairs: crystal contacts are not used; disulphide bonds
(S-S closer than 2.5 A) are covalent and are not counted.
"""
import gzip, os, sys
import numpy as np

RMAX = 8.0
DR = 0.05
IMIN = 6		# as imin of montegrappa: pairs whose backbone atoms are <= IMIN apart are excluded
NH_DIST = 1.01
OWNER = {"N": 0, "H": 0, "CA": 1, "C": 2, "O": 2, "OXT": 2}	# backbone atom each atom belongs to (default CA)
STYPE = -1		# type of the cysteine S (set from the type library)


def read_types(fname):
    rules = []
    for line in open(fname):
        f = line.split()
        if len(f) >= 3 and f[0].lstrip("-").isdigit():
            rules.append((int(f[0]), f[1], f[2]))
    def typeof(atom, aa):
        t = -1
        for it, a, r in rules:			# the last matching line wins, as in grappino
            if (a == atom or a == "*") and (r == aa or r == "*"):
                t = it
        return t
    return typeof, max(r[0] for r in rules) + 1


def read_chain(fname, chain):
    """residues of the chain, in order: list of (resname, {atom: xyz})"""
    res, keys = {}, []
    with gzip.open(fname, "rt") as f:
        for l in f:
            if l.startswith("ENDMDL"):
                break
            rec = l[:6]
            if rec not in ("ATOM  ", "HETATM") or l[21] != chain:
                continue
            resn = l[17:20].strip()
            if rec == "HETATM" and resn != "MSE":
                continue
            if l[16] not in " A":
                continue
            atom = l[12:16].strip()
            elem = l[76:78].strip() if len(l) > 77 else ""
            if elem in ("H", "D") or (not elem and atom[0] in "HD"):
                continue
            if resn == "MSE":
                resn = "MET"
                if atom == "SE":
                    atom = "SD"
            key = (int(l[22:26]), l[26])
            if key not in res:
                res[key] = (resn, {})
                keys.append(key)
            res[key][1].setdefault(atom, np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]))
    return [res[k] for k in keys]


def add_amide_h(residues):
    for k in range(1, len(residues)):
        resn, at = residues[k]
        prev = residues[k-1][1]
        if resn == "PRO" or not all(a in at for a in ("N", "CA")) or "C" not in prev:
            continue
        if np.linalg.norm(at["N"] - prev["C"]) > 1.5:		# chain break
            continue
        u1 = at["N"] - prev["C"]; u1 /= np.linalg.norm(u1)
        u2 = at["N"] - at["CA"]; u2 /= np.linalg.norm(u2)
        x = u1 + u2
        at["H"] = at["N"] + NH_DIST * x / np.linalg.norm(x)


def main(typesfile, setfile, cache, out):
    typeof, ntypes = read_types(typesfile)
    global STYPE
    STYPE = typeof("SG", "CYS")
    nbin = int(round(RMAX / DR))
    hist = np.zeros((ntypes, ntypes, nbin))
    natoms = np.zeros(ntypes)
    nres = nchain = 0
    untyped = {}
    for line in open(setfile):
        if line.startswith("#") or not line.strip():
            continue
        pdb, chain = line.split()[:2]
        fname = os.path.join(cache, pdb + ".pdb.gz")
        if not os.path.exists(fname):
            continue
        residues = read_chain(fname, chain)
        if len(residues) < 30:
            continue
        add_amide_h(residues)
        xyz, typ, ires = [], [], []
        for k, (resn, at) in enumerate(residues):
            for a, x in at.items():
                t = typeof(a, resn)
                if t >= 0:
                    xyz.append(x); typ.append(t); ires.append(3*k + OWNER.get(a, 1))
                else:
                    untyped[(resn, a)] = untyped.get((resn, a), 0) + 1
        xyz = np.array(xyz); typ = np.array(typ); ires = np.array(ires)
        np.add.at(natoms, typ, 1)
        d = np.sqrt(((xyz[:, None, :] - xyz[None, :, :])**2).sum(-1))
        sel = (d < RMAX) & (np.abs(ires[:, None] - ires[None, :]) > IMIN)
        sel &= ~((typ[:, None] == STYPE) & (typ[None, :] == STYPE) & (d < 2.5))	# disulphide bonds are covalent
        i, j = np.nonzero(np.triu(sel, 1))
        b = (d[i, j] / DR).astype(int)
        ti, tj = typ[i], typ[j]
        np.add.at(hist, (ti, tj, b), 1)
        np.add.at(hist, (tj, ti, b), 1)		# symmetric (i!=j atoms: each pair once per order)
        nres += len(residues); nchain += 1
    # the diagonal of the type matrix got each pair twice
    for t in range(ntypes):
        hist[t, t] /= 2
    np.savez(out, hist=hist, natoms=natoms, dr=DR, rmax=RMAX, imin=IMIN, nchain=nchain, nres=nres)
    if untyped:
        print("atoms without a type (ignored):", " ".join(f"{r}:{a}({n})" for (r, a), n in sorted(untyped.items(), key=lambda x: -x[1])[:20]), file=sys.stderr)
    print(f"{nchain} chains, {nres} residues, {int(natoms.sum())} atoms; pairs within {RMAX} A: {int(np.triu(hist.sum(-1)).sum())}", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(*sys.argv[1:])
