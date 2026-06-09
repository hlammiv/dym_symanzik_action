#!/usr/bin/env python3
"""
LLR plaquette heatmap: from each beta2 slice's reconstructed rho(E), reweight to
the *equilibrium* <E> = <plaquette>/3 + 1 as a function of beta1, and assemble a
heatmap over (beta1, beta2) with the freezing line overlaid.

Unlike a direct Monte Carlo scan, this is metastability-free: <E>(beta1) =
sum E rho(E) e^{-beta1 E} / sum rho(E) e^{-beta1 E}.

Usage: llr_heatmap.py <dir-with-b2_*.out> [--nplaq N] [--b1max B] [--plot out.png]
"""
import sys, glob, argparse, importlib.util, os
import numpy as np

_here = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("llr_reconstruct",
            os.path.join(_here, "llr_reconstruct.py"))
R = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(R)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dir")
    ap.add_argument("--nplaq", type=int, default=96)
    ap.add_argument("--b1max", type=float, default=5.5)
    ap.add_argument("--plot", default="results/llr_plaquette_heatmap.png")
    args = ap.parse_args()

    b1 = np.linspace(0.0, args.b1max, 160)
    b2s, rows, fo = [], [], {}
    for f in sorted(glob.glob(os.path.join(args.dir, "b2_*.out"))):
        b2 = float(os.path.basename(f)[3:-4])
        E, a, _ = R.parse(f)
        if len(E) < 3: continue
        a = R.smooth(a, 3); Ef, af, lnrho = R.reconstruct(E, a)
        Edens = Ef/(args.nplaq*3.0) + 1.0                 # <E> per plaquette
        row = np.empty_like(b1)
        for i, b in enumerate(b1):
            lw = lnrho - b*Ef; lw -= lw.max(); w = np.exp(lw)
            row[i] = np.sum(Edens*w)/np.sum(w)
        b2s.append(b2); rows.append(row)
        r = R.find_betac(Ef, lnrho, a)
        if r and r[0] > 0: fo[b2] = r[0]

    b2s = np.array(b2s); heat = np.array(rows)
    o = np.argsort(b2s); b2s = b2s[o]; heat = heat[o]
    # smooth the beta2 axis by linear interpolation onto a finer grid
    b2f = np.linspace(b2s.min(), b2s.max(), 120)
    heatf = np.vstack([np.interp(b2f, b2s, heat[:, i]) for i in range(len(b1))]).T

    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    plt.figure(figsize=(7.5, 5.5))
    pc = plt.pcolormesh(b1, b2f, heatf, shading="auto", cmap="viridis", vmin=0, vmax=1)
    plt.colorbar(pc, label=r"$\langle E\rangle$  (plaquette,  1=disordered, 0=frozen)")
    fb = sorted(fo)
    if fb:
        plt.plot([fo[b] for b in fb], fb, "w.-", lw=1.5, ms=7, label="freezing line (first-order)")
        plt.legend(loc="upper right", framealpha=0.6)
    plt.xlabel(r"$\beta_1$"); plt.ylabel(r"$\beta_2$")
    plt.title(r"S1080 improved-action plaquette $\langle E\rangle(\beta_1,\beta_2)$  (LLR, 2$^4$)")
    plt.tight_layout(); plt.savefig(args.plot, dpi=120)
    print(f"plot -> {args.plot}")
    print("freezing line (beta1_c, beta2):", ", ".join(f"({fo[b]:.2f},{b:+.2f})" for b in fb))

if __name__ == "__main__":
    main()
