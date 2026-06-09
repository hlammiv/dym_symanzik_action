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

    # reweight over a continuous interpolated rho (smooth moments) or raw cells.
    cE1,cE2=cells[:,0],cells[:,1]
    if args.smooth>0:
        try:
            from scipy.interpolate import griddata
            g1=np.linspace(cE1.min(),cE1.max(),args.smooth)
            g2=np.linspace(cE2.min(),cE2.max(),args.smooth)
            G1,G2=np.meshgrid(g1,g2)
            LR=griddata((cE1,cE2),lnrho,(G1,G2),method="linear")
            m=~np.isnan(LR); E1=G1[m]; E2=G2[m]; LRv=LR[m]
        except Exception as e:
            print(f"(smooth off: {e})"); E1,E2,LRv=cE1,cE2,lnrho
    else:
        E1,E2,LRv=cE1,cE2,lnrho
    S = E2 if is_rect else E1

    b1=np.linspace(0,args.b1max,args.n); b2=np.linspace(-args.b2lim,args.b2lim,args.n)
    # first moment <S>(b1,b2) -- smooth in the couplings even from a coarse grid.
    m1map=np.empty((args.n,args.n))
    for ii,bb1 in enumerate(b1):
        lw0=LRv-bb1*E1
        for jj,bb2 in enumerate(b2):
            lw=lw0-bb2*E2; lw-=lw.max(); w=np.exp(lw)
            m1map[jj,ii]=np.sum(S*w)/np.sum(w)
    if is_chi:
        # chi = -d<S>/d(its coupling) * V/Naction^2  (fluctuation-dissipation).
        # We differentiate a LIGHTLY-SMOOTHED <S> rather than take the direct
        # variance: the direct Var(S) over discrete cells is spiky (every pair of
        # cells makes a pseudo-coexistence ridge), whereas <S> is smooth, so its
        # derivative gives the clean susceptibility ridge.
        from scipy.ndimage import gaussian_filter
        sm=gaussian_filter(m1map,sigma=args.chi_sigma)
        d=np.gradient(sm,b2,axis=0) if is_rect else np.gradient(sm,b1,axis=1)
        phi=-d*Vsites/(Naction*Naction)
    else:
        phi=m1map/(Naction*3.0)+1.0                          # 0=frozen, 1=disordered

    # freezing line (steepest drop of <S1>), overlaid
    line=[(x,R.betac_steep(cells,lnrho,'b1',x,-args.b2lim,0.55)) for x in np.linspace(1.4,args.b1max,16)]
    b1c0=R.betac_steep(cells,lnrho,'b2',0.0,0.3,3.0)

    import matplotlib;matplotlib.use("Agg");import matplotlib.pyplot as plt
    plt.figure(figsize=(8,6))
    which = "rect" if is_rect else "plaq"
    if is_chi:
        cc = "beta_1" if not is_rect else "beta_2"
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
