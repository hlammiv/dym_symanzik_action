#!/usr/bin/env python3
"""
LLR milestone 4: reconstruct the density of states from the a_n(E) curve and
locate the first-order transition.

Input: the stdout of `llr` (the ANE: / NPLAQ: lines), passed as a file or stdin.

Pipeline (paper Eq. 2, 11):
  1. a_n = d ln rho / dE in each window -> ln rho(E) by integrating a_n.
  2. P_beta(E) = rho(E) e^{-beta E}  (beta is the beta1-conjugate coupling).
  3. beta1_c = the beta where the two maxima of ln P_beta(E) have equal height
     (equal free energy of the two phases) -- equivalently the Maxwell equal-area
     construction on the a_n(E) hump.
  4. latent heat = separation of the two peaks (in S1, and in <E>=S1/NPLAQ/3+1).

Usage: llr_reconstruct.py ladder.out [--nplaq N] [--plot results/llr.png]
"""
import sys, argparse
import numpy as np
_trap = getattr(np, "trapezoid", getattr(np, "trapz", None))

def parse(fname):
    """Returns (E, a_n, norm, axis): norm = NPLAQ (axis 1) or NRECT (axis 2),
    axis = which action is windowed (1=S1/beta1-driven, 2=S2/beta2-driven)."""
    Es, ans = [], []
    nplaq = nrect = None; axis = 1
    f = open(fname) if fname and fname != "-" else sys.stdin
    for line in f:
        if   line.startswith("NPLAQ:"): nplaq = int(line.split()[1])
        elif line.startswith("NRECT:"): nrect = int(line.split()[1])
        elif line.startswith("AXIS:"):  axis  = int(line.split()[1])
        elif line.startswith("ANE:") and "unreached" not in line:
            t = line.split(); Es.append(float(t[1])); ans.append(float(t[2]))
    E = np.array(Es); a = np.array(ans)
    o = np.argsort(E)                 # ascending E (frozen -> disordered)
    norm = nrect if axis == 2 else nplaq
    return E[o], a[o], norm, axis

def smooth(a, w=3):
    """Light moving-average to suppress single-point noise at the frozen floor
    (the deepest windows are under-resolved and occasionally spike)."""
    if w <= 1 or len(a) < w: return a
    k = w//2; s = a.astype(float).copy()
    for i in range(len(a)):
        s[i] = np.mean(a[max(0,i-k):min(len(a),i+k+1)])
    return s

def reconstruct(E, a, npts=4000):
    Ef = np.linspace(E[0], E[-1], npts)
    af = np.interp(Ef, E, a)
    # ln rho(E) = integral of a dE  (rho fixed up to an arbitrary constant)
    lnrho = np.concatenate([[0.0], np.cumsum(0.5*(af[1:]+af[:-1])*np.diff(Ef))])
    return Ef, af, lnrho

def peaks(Ef, lnrho, beta):
    """At coupling beta: split E at the barrier (interior min of g=lnrho-beta*E)
    into frozen (low E) / disordered (high E) basins, and return the height
    difference g(frozen peak) - g(disordered peak) and the two peak energies.
    Returns None if there is no interior barrier (single phase at this beta).
    The frozen peak may sit at the low-E boundary (we sample down to the floor)."""
    g = lnrho - beta*Ef
    j = 1 + int(np.argmin(g[1:-1]))            # barrier (interior minimum)
    if not (g[j] < g[j-1] and g[j] < g[j+1]):  # not a genuine well -> single phase
        return None
    gf = g[:j+1].max();  Ef_f = Ef[:j+1][np.argmax(g[:j+1])]
    gd = g[j:].max();    Ef_d = Ef[j:][np.argmax(g[j:])]
    return gf - gd, Ef_f, Ef_d

def find_betac(Ef, lnrho, a):
    """beta_c where the two ln P peaks have equal height (equal free energy,
    paper Eq. 11). The interior barrier only exists over the resolved coexistence
    range; g(frozen)-g(disord) is ~linear in beta there, so fit it and take the
    zero crossing (extrapolating past the spinodal where the barrier closes).
    Returns (beta_c, E_frozen, E_disord, beta_lo, beta_hi). beta_hi is the most
    bimodal resolved beta (for plotting)."""
    betas, diffs, Efs, Eds = [], [], [], []
    for b in np.linspace(a.min()+1e-3, a.max()-1e-3, 800):
        r = peaks(Ef, lnrho, b)
        if r is None: continue
        betas.append(b); diffs.append(r[0]); Efs.append(r[1]); Eds.append(r[2])
    if len(betas) < 3:
        return None
    betas = np.array(betas); diffs = np.array(diffs)
    sl, ic = np.polyfit(betas, diffs, 1)       # g_f - g_d  vs  beta
    bc = -ic/sl
    k = int(np.argmin(np.abs(betas - bc)))     # nearest resolved beta for peak Es
    return bc, Efs[k], Eds[k], betas[0], betas[-1]

def _rwmean(Ef, lnrho, beta):
    g = lnrho - beta*Ef; g -= g.max(); w = np.exp(g)
    return np.sum(Ef*w)/np.sum(w)

def transition(Ef, lnrho, a, b1max=6.0, npts=600):
    """beta1_c from the steepest drop of <S1>(beta1) -- robust for every slice,
    including where a_n goes negative. Latent heat from find_betac (the a_n
    equal-height extrapolation), which resolves where the back-bend is sharp
    enough; None otherwise (the frozen peak is suppressed by floor under-
    resolution, so a P(S1) double-peak read-off is unreliable). Returns
    (beta1_c, E_frozen, E_disord, latent_or_None, resolved)."""
    b1 = np.linspace(0.0, b1max, npts)
    Em = np.array([_rwmean(Ef, lnrho, b) for b in b1])
    bc = b1[int(np.argmax(np.abs(np.gradient(Em, b1))))]
    r = find_betac(Ef, lnrho, a)
    if r is not None:
        return bc, r[1], r[2], (r[2]-r[1]), True
    return bc, None, None, None, False

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("file", nargs="?", default="-")
    ap.add_argument("--nplaq", type=int, default=None)
    ap.add_argument("--plot", default=None)
    ap.add_argument("--smooth", type=int, default=3, help="moving-average window on a_n (1=off)")
    args = ap.parse_args()

    E, a, norm, axis = parse(args.file)
    if args.nplaq: norm = args.nplaq
    if len(E) < 3:
        sys.exit("need >=3 reached windows")
    a = smooth(a, args.smooth)
    cc = "beta2" if axis == 2 else "beta1"     # the windowed coupling
    en = "S2" if axis == 2 else "S1"           # the windowed energy

    Ef, af, lnrho = reconstruct(E, a)
    bc, Elo, Ehi, latent, bimodal = transition(Ef, lnrho, a)
    print(f"{cc}_c (steepest drop of <{en}>) = {bc:.4f}")
    if bimodal:
        print(f"coexistence {en} : E- = {Elo:.2f} (frozen)  E+ = {Ehi:.2f} (disordered)")
        print(f"latent heat (Delta {en})  : {latent:.2f}")
        if norm:
            print(f"  in <E> units (norm={norm}) : Delta<E> = {latent/(norm*3.0):.4f}")
    else:
        print(f"P({en}) not bimodal at {cc}_c -> transition weak/unresolved (no latent heat)")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        ax[0].plot(Ef, af, '-', lw=1.5); ax[0].plot(E, a, 'o', ms=3)
        ax[0].axhline(bc, color='r', ls='--', lw=1, label=f"$\\beta_{{1c}}$={bc:.3f}")
        if bimodal:
            ax[0].axvline(Elo, color='gray', ls=':'); ax[0].axvline(Ehi, color='gray', ls=':')
        ax[0].set_xlabel(f"E = {en}"); ax[0].set_ylabel(r"$a_n(E)=\beta_{micro}$")
        ax[0].set_title("DoS slope"); ax[0].legend()
        g = lnrho - bc*Ef; P = np.exp(g - g.max())
        ax[1].fill_between(Ef, P, alpha=0.3); ax[1].plot(Ef, P, lw=1.5)
        ax[1].set_xlabel(f"E = {en}"); ax[1].set_ylabel(r"$P_{\beta_{1c}}(E)$")
        ax[1].set_title(f"P at beta1_c={bc:.3f}")
        plt.tight_layout(); plt.savefig(args.plot, dpi=110)
        print(f"plot -> {args.plot}")

if __name__ == "__main__":
    main()
