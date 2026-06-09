#!/usr/bin/env python3
"""
Plaquette (and rectangle) order-parameter heatmap over (beta1,beta2) from the
2D-LLR reconstructed rho(S1,S2). ONE reconstruction -> the whole phase diagram:
metastability-free, and -- unlike the 1D-stitched map -- NO spectator problem
and NO inter-axis disagreement wedge, because both couplings reweight the same
joint density.

  phi_plaq(b1,b2) = <S1>/(NPLAQ*3) + 1   (0=frozen, 1=disordered)
  <S1> = sum_cells S1 e^{lnrho - b1 S1 - b2 S2} / sum_cells e^{...}

Usage: llr2d_heatmap.py run.out [run2.out ...] [--obs plaq|rect] [--plot out.png]
"""
import sys, argparse, importlib.util, os, numpy as np
_here=os.path.dirname(os.path.abspath(__file__))
_spec=importlib.util.spec_from_file_location("llr2d_reconstruct",
        os.path.join(_here,"llr2d_reconstruct.py"))
R=importlib.util.module_from_spec(_spec); _spec.loader.exec_module(R)

def norms_from_header(files):
    """(NPLAQ, NRECT, Vsites) from the 'D Nt Nx' in the LLR2D header."""
    for f in files:
        for l in open(f):
            if l.startswith("LLR2D"):
                t=l.split()
                try:
                    D=int(t[2]); Nt=int(t[3]); Nx=int(t[4])
                    V=Nt*Nx**(D-1); return V*D*(D-1)//2, V*D*(D-1), V
                except Exception: pass
            if l.startswith("ANE2:"): break
    return 1536, 3072, 256   # 4^4 fallback

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("file",nargs="+")
    ap.add_argument("--obs",choices=["plaq","rect","chi","chirect"],default="plaq",
                    help="plaq/rect order parameter, or chi/chirect = its connected susceptibility")
    ap.add_argument("--b1max",type=float,default=4.0)
    ap.add_argument("--b2lim",type=float,default=1.6)
    ap.add_argument("--n",type=int,default=180)
    ap.add_argument("--chi-sigma",type=float,default=1.5,
                    help="Gaussian smoothing (in pixels) of <S> before differentiating for chi")
    ap.add_argument("--smooth",type=int,default=140,
                    help="interpolate lnrho onto an NxN (S1,S2) mesh before reweighting "
                         "(0=off, use raw cell point-masses); essential for susceptibilities")
    ap.add_argument("--plot",default="results/llr2d_plaquette_heatmap.png")
    args=ap.parse_args()

    if len(args.file)==1: cells,geo=R.parse(args.file[0])
    else: cells,geo,_=R.parse_multi(args.file)
    lnrho=R.integrate(cells,R.grid_index(cells,geo))
    NPLAQ,NRECT,Vsites=norms_from_header(args.file)
    is_rect = args.obs in ("rect","chirect")
    is_chi  = args.obs in ("chi","chirect")
    Naction = (NRECT if is_rect else NPLAQ)
    h1=geo["step1"]/2.0; h2=geo["step2"]/2.0                  # non-overlapping tiles

    b1=np.linspace(0,args.b1max,args.n); b2=np.linspace(-args.b2lim,args.b2lim,args.n)
    # CONTINUOUS-cell <S>: each cell is an exp-linear density over its tile (uses the
    # measured gradient), integrated analytically -> the most accurate <S>(b1,b2).
    Smap=np.empty((args.n,args.n))
    for ii,bb1 in enumerate(b1):
        for jj,bb2 in enumerate(b2):
            S1,_,S2,_ = R.cell_moments(cells,lnrho,h1,h2,bb1,bb2)
            Smap[jj,ii]= S2 if is_rect else S1
    if is_chi:
        # NOTE: the *direct* Var(S) (point-mass OR continuous-cell) facets -- a finite
        # noisy cell set gives <S> tiny kinks that the variance amplifies into spurious
        # pseudo-coexistence ridges, and that does NOT wash out with density/statistics.
        # So chi = -d<S>/d(coupling) * V/N^2 from a lightly-smoothed <S>, which discards
        # exactly those kinks and leaves the clean physical ridge.
        from scipy.ndimage import gaussian_filter
        sm=gaussian_filter(Smap,sigma=args.chi_sigma)
        d=np.gradient(sm,b2,axis=0) if is_rect else np.gradient(sm,b1,axis=1)
        phi=-d*Vsites/(Naction*Naction)
    else:
        phi=Smap/(Naction*3.0)+1.0                               # 0 frozen, 1 disordered

    # freezing line (steepest drop of <S1>), overlaid
    line=[(x,R.betac_steep(cells,lnrho,'b1',x,-args.b2lim,0.55)) for x in np.linspace(1.4,args.b1max,16)]
    b1c0=R.betac_steep(cells,lnrho,'b2',0.0,0.3,3.0)

    import matplotlib;matplotlib.use("Agg");import matplotlib.pyplot as plt
    plt.figure(figsize=(8,6))
    which = "rect" if is_rect else "plaq"
    if is_chi:
        cc = "beta_2" if is_rect else "beta_1"
        vmax=np.percentile(phi[phi>0],99.0) if np.any(phi>0) else None
        pc=plt.pcolormesh(b1,b2,phi,shading="auto",cmap="inferno",vmin=0,vmax=vmax)
        plt.colorbar(pc,label=rf"$\chi_{{{which}}} = -\partial\langle S\rangle/\partial\{cc}\cdot V/N^2$")
        linecol="cyan"; title=rf"S1080 {which} susceptibility from 2D LLR $\rho(S_1,S_2)$ (4$^4$)"
    else:
        pc=plt.pcolormesh(b1,b2,phi,shading="auto",cmap="viridis",vmin=0,vmax=1)
        plt.colorbar(pc,label=rf"$\langle\phi_{{{which}}}\rangle$  (1=disordered, 0=frozen)")
        linecol="white"; title=rf"S1080 {which} order parameter from 2D LLR $\rho(S_1,S_2)$ (4$^4$)"
    xx,yy=zip(*line); plt.plot([b1c0]+list(xx),[0.0]+list(yy),".-",color=linecol,lw=1.6,ms=6,label="freezing line (2D LLR)")
    plt.legend(loc="upper right",framealpha=.6)
    plt.xlabel(r"$\beta_1$"); plt.ylabel(r"$\beta_2$")
    plt.title(title)
    plt.tight_layout(); plt.savefig(args.plot,dpi=120)
    print("plot ->",args.plot)
    print(f"NPLAQ={NPLAQ} NRECT={NRECT} V={Vsites}; freezing line beta2_c:",
          ", ".join(f"({x:.1f},{y:+.2f})" for x,y in line))

if __name__=="__main__":
    main()
