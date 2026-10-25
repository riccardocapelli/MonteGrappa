#!/usr/bin/env python3
"""
Square-well pair potential between reduced atom types from the statistics of
contact_stats.py.

    derive_potential.py stats.npz pairs.dat [--plot DIR] [--hb S] [--polar S] [--hp EPS RW]

For each pair of types, the potential of mean force is
    w(r) = -ln[ N(r) / Nref(r) ],   Nref(r) = N(rc) (r/rc)^alpha
(DFIRE reference state, Zhou and Zhou, Protein Sci. 11, 2714 (2002): an ideal
gas in a finite volume, alpha = 1.61, normalised in the shell [rc-1, rc]).
The square well that replaces it:
  - hard core r0: the distance below which a fraction QHARD (1e-4) of the pairs is
    found (rounded down to 0.05 A);
  - depth e: -ln of the mean of N/Nref between r0 and rw.
  - outer radius rw: the end of the first attractive shell, i.e. the first
    distance after the minimum of w(r) where w(r) becomes positive (at most
    RWMAX);
Pairs whose first shell is not attractive (e >= 0) only keep their hard core:
their positive mean comes mostly from the rise of w(r) next to the hard core,
i.e. from the excluded volume that the hard core already accounts for, and
summed over the many close contacts of a compact structure it would dominate
the energy. Pairs listed in HB_CORRELATES (N-O, H-C, H-N, H-H) only keep their
hard core: their wells are geometric consequences of the backbone hydrogen
bond, which is the H-O well with its angular factor (keyword hb). Energies are in units of kT at the temperature of the crystal
structures (T=1 in montegrappa is then about room temperature).

The atom-level statistics hardly see the hydrophobic effect in the apolar pairs
(Cal-Cal is about -0.05 kT): it is mostly in the depletion of polar atoms next
to apolar ones, i.e. in the repulsive terms that are dropped. So the pairs of
apolar types (APOLAR) get instead a well of depth EPS out to RW (--hp), with
their statistical hard core. The other wells can be scaled: --hb for the H-O
hydrogen bond, --polar for all the remaining attractive pairs. The scales are
meant to be calibrated (stability of a fibril, behaviour of the monomer).

Output lines: "type1 type2 e rw r0", the format of the [ pairs ] chapter.
"""
import sys
import numpy as np

ALPHA = 1.61
RC = 8.0
QHARD = 0.0001
RWMAX = 5.5
HB_CORRELATES = {tuple(sorted(x)) for x in (("N", "O"), ("C", "H"), ("N", "H"), ("H", "H"))}
MINCOUNT = 200			# fewer pairs in the contact shell: no well (only the hard core)
NAMES = "N CA C O H Cal Car Cpo Opo Npo Oneg Npos S".split()
APOLAR = {"Cal", "Car", "S"}


def well(x, r, dr):
    tot = x.sum()
    if tot == 0:
        return None
    c = np.cumsum(x) / tot
    r0 = np.floor(r[np.searchsorted(c, QHARD)] / 0.05) * 0.05
    shell = (r > RC - 1) & (r <= RC)
    nref = x[shell].sum() / ((r[shell] / RC) ** ALPHA).sum()
    ratio = x / (nref * (r / RC) ** ALPHA)
    with np.errstate(divide="ignore"):
        w = -np.log(ratio)
    m = (r > r0) & (r < RWMAX) & np.isfinite(w)
    if not m.any():
        return None
    imin = np.nonzero(m)[0][np.argmin(w[m])]
    if w[imin] < 0:
        k = imin
        while k + 1 < len(r) and r[k + 1] < RWMAX and w[k + 1] < 0:
            k += 1
        rw = r[k] + dr / 2
    else:
        rw = r0 + 1.5
    sel = (r > r0) & (r < rw)
    if x[sel].sum() < MINCOUNT:
        return (0.0, 0.0, r0, x[sel].sum(), w)
    e = -np.log(x[sel].sum() / (nref * ((r[sel] / RC) ** ALPHA)).sum())
    return (e, rw, r0, x[sel].sum(), w)


def main(statsfile, outfile, plotdir=None, hb=1., polar=1., hp=None):
    d = np.load(statsfile)
    h, dr = d["hist"], float(d["dr"])
    nt = h.shape[0]
    r = (np.arange(h.shape[2]) + 0.5) * dr
    out = open(outfile, "w")
    out.write(f"# derive_potential.py {statsfile}: {int(d['nchain'])} chains, {int(d['nres'])} residues\n")
    out.write(f"# scales: hb {hb}, polar {polar}, hydrophobic {'none' if hp is None else f'{hp[0]} out to {hp[1]} A'}\n")
    out.write("# type1 type2   e(kT)   r_well   r_hardcore   (contacts in the well)\n")
    curves = {}
    for a in range(nt):
        for b in range(a, nt):
            res = well(h[a, b], r, dr)
            if res is None:
                continue
            e, rw, r0, n, w = res
            curves[(a, b)] = w
            names = tuple(sorted((NAMES[a], NAMES[b]))) if nt == len(NAMES) else None
            if e >= 0 or (names and names in HB_CORRELATES):
                e, rw = 0.0, 0.0
            if names == ("H", "O"):
                e *= hb
            elif names and hp is not None and set(names) <= APOLAR:
                e, rw = -hp[0], max(hp[1], r0 + 0.5)
            else:
                e *= polar
            name = f"{NAMES[a]}-{NAMES[b]}" if nt == len(NAMES) else ""
            out.write(f"{a:3d}\t{b:3d}\t\t{e:9.4f}\t{rw:6.2f}\t{r0:6.2f}\t# {name} {int(n)}\n")
    out.close()
    if plotdir:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        for a in range(nt):
            fig, ax = plt.subplots(figsize=(5, 3.5))
            for b in range(nt):
                key = (min(a, b), max(a, b))
                if key in curves:
                    ax.plot(r, curves[key], label=NAMES[b] if nt == len(NAMES) else str(b), lw=1)
            ax.set_xlim(1.5, 8); ax.set_ylim(-2, 3); ax.axhline(0, color="k", lw=0.5)
            ax.set_xlabel("r (A)"); ax.set_ylabel("w(r) / kT")
            ax.set_title(NAMES[a] if nt == len(NAMES) else str(a)); ax.legend(fontsize=6, ncol=2)
            fig.tight_layout(); fig.savefig(f"{plotdir}/pmf_{a}.png", dpi=120); plt.close(fig)


if __name__ == "__main__":
    args = sys.argv[1:]
    opt = {}
    for key, n in (("--plot", 1), ("--hb", 1), ("--polar", 1), ("--hp", 2)):
        if key in args:
            k = args.index(key); opt[key] = args[k + 1:k + 1 + n]; del args[k:k + 1 + n]
    if len(args) != 2:
        sys.exit(__doc__)
    main(args[0], args[1], opt["--plot"][0] if "--plot" in opt else None,
         float(opt.get("--hb", [1])[0]), float(opt.get("--polar", [1])[0]),
         tuple(map(float, opt["--hp"])) if "--hp" in opt else None)
