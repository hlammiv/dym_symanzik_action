#!/usr/bin/env python3
"""
LLR order-parameter heatmap over (beta1, beta2), combining BOTH LLR axes.

Each LLR run reconstructs an effective density of states rho~(E) at a fixed
background coupling, with E = the *windowed* action:

  axis-1 file (fixed beta2):  rho~_{b2}(S1)  ->  phi1(b1) = <S1>/(NPLAQ*3)+1
                              reweight in beta1  (fills a ROW at that b2)
  axis-2 file (fixed beta1):  rho~_{b1}(S2)  ->  phi2(b2) = <S2>/(NRECT*3)+1
                              reweight in beta2  (fills a COLUMN at that b1)

Both phi are normalized order parameters (0 = frozen, 1 = disordered) for the
SAME freezing transition. axis-1 has full leverage in beta1 but is nearly blind
to beta2 (S2 barely fluctuates at fixed S1, esp. in the frozen branch), and
vice-versa. Averaging the two fields per cell uses whichever axis actually
resolves the transition in that region -- so the negative-beta2 line (where
axis-1 saturates) is corrected by axis-2.

Metastability-free: <E>(b) = sum E rho~(E) e^{-b E} / sum rho~(E) e^{-b E}.

Usage: llr_heatmap.py <dir-or-glob> [<dir-or-glob>...]
       [--b1max B] [--b2lim L] [--plot out.png]
Files are auto-routed by their AXIS: header (axis-1 vs axis-2); the background
coupling is read from the LLR(...) header line.
"""
import sys, glob, argparse, importlib.util, os
import numpy as np

_here = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("llr_reconstruct",
            os.path.join(_here, "llr_reconstruct.py"))
R = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(R)


def read_bg(fname):
    """Background coupling from the LLR(...) header (beta2 if axis-1, beta1 if
    axis-2). Layout: 'LLR(...): group D Nt Nx axis bg | ...'."""
    for line in open(fname):
        if line.startswith("LLR("):
            t = line.split()
            try: return float(t[6])
            except (IndexError, ValueError): return None
        if line.startswith("ANE:"): break
    return None


def opcurve(E, a, norm, grid):
    """Order parameter phi(coupling) = <E>/(norm*3)+1 reweighted along `grid`
    from a single slice's rho~(E). P(E) ~ rho~(E) e^{-c E}."""
    Ef, af, lnrho = R.reconstruct(E, a)
    dens = Ef/(norm*3.0) + 1.0
    out = np.empty_like(grid)
    for i, c in enumerate(grid):
        lw = lnrho - c*Ef; lw -= lw.max(); w = np.exp(lw)
        out[i] = np.sum(dens*w)/np.sum(w)
    return out


def collect(paths):
    """Expand dirs/globs to a flat list of .out files."""
    files = []
    for p in paths:
        if os.path.isdir(p): files += glob.glob(os.path.join(p, "*.out"))
        else:               files += glob.glob(p)
    return sorted(set(files))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("paths", nargs="+", help="dir(s) or glob(s) of llr .out files")
    ap.add_argument("--b1max", type=float, default=5.5)
    ap.add_argument("--b2lim", type=float, default=1.8, help="|beta2| extent of the y-axis")
    ap.add_argument("--n1", type=int, default=200)
    ap.add_argument("--n2", type=int, default=200)
    ap.add_argument("--plot", default="results/llr_plaquette_heatmap.png")
    ap.add_argument("--plot-agree", default="results/llr_heatmap_agreement.png",
                    help="mean(color)/disagreement(whiteness) consensus map")
    ap.add_argument("--spread-ref", type=float, default=0.5,
                    help="|phi1-phi2| at which a cell is fully whited-out")
    args = ap.parse_args()

    b1g = np.linspace(0.0, args.b1max, args.n1)
    b2g = np.linspace(-args.b2lim, args.b2lim, args.n2)

    ax1, ax2 = [], []          # (bg, curve), (bg, curve)
    fo1, fo2 = {}, {}          # first-order transition points per axis
    for f in collect(args.paths):
        E, a, norm, axis = R.parse(f)
        if E is None or len(E) < 3 or not norm: continue
        bg = read_bg(f)
        if bg is None: continue
        a = R.smooth(a, 3)
        Ef, af, lnrho = R.reconstruct(E, a)
        if axis == 1:                                   # fixed beta2 -> row over beta1
            ax1.append((bg, opcurve(E, a, norm, b1g)))
            r = R.find_betac(Ef, lnrho, a)
            if r and r[0] > 0: fo1[bg] = r[0]           # (beta1_c, beta2=bg)
        else:                                           # fixed beta1 -> column over beta2
            ax2.append((bg, opcurve(E, a, norm, b2g)))
            r = R.find_betac(Ef, lnrho, a)
            if r is not None: fo2[bg] = r[0]            # (beta1=bg, beta2_c)

    # --- field 1: axis-1, defined on rows (fixed beta2), interp across beta2 ---
    F1 = np.full((args.n2, args.n1), np.nan)
    if ax1:
        ax1.sort(key=lambda t: t[0]); b2s = np.array([b for b, _ in ax1]); M = np.array([c for _, c in ax1])
        for i in range(args.n1):
            F1[:, i] = np.interp(b2g, b2s, M[:, i], left=np.nan, right=np.nan)
    # --- field 2: axis-2, defined on columns (fixed beta1), interp across beta1 ---
    F2 = np.full((args.n2, args.n1), np.nan)
    if ax2:
        ax2.sort(key=lambda t: t[0]); b1s = np.array([b for b, _ in ax2]); C = np.array([c for _, c in ax2])
        for j in range(args.n2):
            F2[j, :] = np.interp(b1g, b1s, C[:, j], left=np.nan, right=np.nan)

    # --- per-cell average of whichever axes resolve the cell ---
    stack = np.stack([F1, F2])
    heat = np.nanmean(stack, axis=0)               # NaN only where BOTH missing
    ncov = np.sum(~np.isnan(stack), axis=0)        # 0,1,2 axes contributing

    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    plt.figure(figsize=(8, 6))
    pc = plt.pcolormesh(b1g, b2g, heat, shading="auto", cmap="viridis", vmin=0, vmax=1)
    plt.colorbar(pc, label=r"$\langle\phi\rangle$  (order param,  1=disordered, 0=frozen)")
    # hatch cells covered by only one axis (less reliable)
    plt.contourf(b1g, b2g, (ncov == 1).astype(float), levels=[0.5, 1.5],
                 hatches=["////"], colors="none", alpha=0)
    # transition points: axis-1 (beta1_c, beta2) and axis-2 (beta1, beta2_c)
    if fo1:
        b = sorted(fo1); plt.plot([fo1[x] for x in b], b, "w.-", lw=1.4, ms=8,
                                  label=r"axis-1 $\beta_{1c}(\beta_2)$")
    if fo2:
        b = sorted(fo2); plt.plot(b, [fo2[x] for x in b], "r.-", lw=1.4, ms=8,
                                  label=r"axis-2 $\beta_{2c}(\beta_1)$")
    if fo1 or fo2: plt.legend(loc="upper right", framealpha=0.7)
    plt.xlabel(r"$\beta_1$"); plt.ylabel(r"$\beta_2$")
    plt.title(r"S1080 improved-action order parameter $\langle\phi\rangle(\beta_1,\beta_2)$"
              "\n(LLR, both axes averaged; hatched = single-axis only)")
    plt.tight_layout(); plt.savefig(args.plot, dpi=120)

    # --- consensus map: hue = mean phi, whiteness = disagreement between axes ---
    # confidence = 1 - |phi1-phi2|/spread_ref  (clamped). Cells covered by both axes
    # blend toward white as the two scans diverge; single-axis cells are shown
    # faded (no cross-check possible); uncovered cells are pure white.
    spread = np.abs(F1 - F2)                              # nan unless both present
    conf = np.clip(1.0 - spread/args.spread_ref, 0.0, 1.0)
    conf = np.where(ncov == 2, conf, np.where(ncov == 1, 0.30, 0.0))
    rgb = plt.get_cmap("viridis")(np.clip(np.nan_to_num(heat), 0, 1))[..., :3]
    w = conf[..., None]
    img = w*rgb + (1.0 - w)*np.ones_like(rgb)             # blend toward white
    fig2 = plt.figure(figsize=(8, 6)); ax = fig2.add_subplot(111)
    ax.imshow(img, origin="lower", aspect="auto",
              extent=[b1g[0], b1g[-1], b2g[0], b2g[-1]])
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(0, 1))
    fig2.colorbar(sm, ax=ax, label=r"mean $\langle\phi\rangle$  (hue; 1=disordered, 0=frozen)")
    if fo1:
        b = sorted(fo1); ax.plot([fo1[x] for x in b], b, "k.-", lw=1.2, ms=7,
                                 label=r"axis-1 $\beta_{1c}(\beta_2)$")
    if fo2:
        b = sorted(fo2); ax.plot(b, [fo2[x] for x in b], "r.-", lw=1.2, ms=7,
                                 label=r"axis-2 $\beta_{2c}(\beta_1)$")
    if fo1 or fo2: ax.legend(loc="upper right", framealpha=0.7)
    ax.set_xlabel(r"$\beta_1$"); ax.set_ylabel(r"$\beta_2$")
    ax.set_title(r"S1080 consensus map: hue = mean of both LLR axes,"
                 "\nwhiteness = inter-axis disagreement (white = the two scans differ)")
    fig2.tight_layout(); fig2.savefig(args.plot_agree, dpi=120)

    print(f"plot -> {args.plot}")
    print(f"plot -> {args.plot_agree}")
    print(f"axis-1 slices (fixed beta2): {len(ax1)};  axis-2 slices (fixed beta1): {len(ax2)}")
    if fo1: print("axis-1 line (beta1_c, beta2):",
                  ", ".join(f"({fo1[b]:.2f},{b:+.2f})" for b in sorted(fo1)))
    if fo2: print("axis-2 line (beta1, beta2_c):",
                  ", ".join(f"({b:.2f},{fo2[b]:+.2f})" for b in sorted(fo2)))


if __name__ == "__main__":
    main()
