#!/usr/bin/env python3
"""
Assemble the LLR freezing line beta1_c(beta2) from a set of per-slice ladder
outputs (one `llr` run per beta2, files named b2_<value>.out).

For each slice: reconstruct rho(E) (via llr_reconstruct), locate beta1_c and the
latent heat. Slices whose a_n(E) shows no back-bend are flagged crossover. Plots
beta1_c(beta2) and Delta<E>(beta2) and marks where the first-order signal ends.

Usage: llr_transition_line.py <dir-with-b2_*.out> [--nplaq N] [--plot out.png]
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
    ap.add_argument("--smooth", type=int, default=3)
    ap.add_argument("--plot", default=None)
    args = ap.parse_args()

    fo, dE, cross = {}, {}, []
    for f in sorted(glob.glob(os.path.join(args.dir, "b2_*.out"))):
        b2 = float(os.path.basename(f)[3:-4])
        E, a, npl = R.parse(f)
        if len(E) < 3: continue
        a = R.smooth(a, args.smooth)
        Ef, af, lnrho = R.reconstruct(E, a)
        res = R.find_betac(Ef, lnrho, a)
        if res and res[0] > 0:
            fo[b2] = res[0]; dE[b2] = (res[2]-res[1])/(npl*3.0)
        else:
            cross.append(b2)

    bx = sorted(fo)
    print("# beta2   beta1_c   Delta<E>")
    for b in bx: print(f"  {b:+.3f}  {fo[b]:8.3f}  {dE[b]:8.3f}")
    if cross: print("# crossover (no resolved back-bend): " +
                    ", ".join(f"{b:+.2f}" for b in sorted(cross)))
    # endpoint = boundary of the contiguous first-order region on the high-beta2 side
    if cross and bx:
        cpos = [c for c in cross if c > min(bx)]
        endpoint = (max(bx) + min(cpos))/2 if cpos else None
        if endpoint is not None:
            print(f"# first-order signal ends near beta2 ~ {endpoint:.2f}")

    if args.plot:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        ax[0].plot([fo[b] for b in bx], bx, 'o-', ms=6)   # beta1_c on x, beta2 on y
        ax[0].set_xlabel(r'$\beta_{1c}$'); ax[0].set_ylabel(r'$\beta_2$')
        ax[0].set_title('freezing line')
        db = sorted(dE)
        ax[1].plot(db, [dE[b] for b in db], 's-', ms=6, color='green')
        ax[1].set_xlabel(r'$\beta_2$'); ax[1].set_ylabel(r'$\Delta\langle E\rangle$')
        ax[1].set_title('latent heat')
        plt.tight_layout(); plt.savefig(args.plot, dpi=120)
        print(f"# plot -> {args.plot}")

if __name__ == "__main__":
    main()
