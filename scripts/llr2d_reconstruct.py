#!/usr/bin/env python3
"""
2D-LLR reconstruction (genuine 2D): least-squares-integrate the per-cell gradient
a=(a1,a2)=(dlnrho/dS1,dlnrho/dS2) into lnrho(S1,S2) on the tiled band, then reweight
to ANY (beta1,beta2) with NO spectator and locate the freezing line.

Steps:
  1. parse ANE2 cells; recover (i,j) grid indices from the header geometry.
  2. lnrho by least squares over grid edges: for adjacent reached cells a,b,
       lnrho_b - lnrho_a ~ 0.5 (grad_a+grad_b).(x_b-x_a)     (curl-noise tolerant)
     (one cell gauge-fixed to 0). This is the 2D analogue of integrating a_n.
  3. effective free energy along S1 at (b1,b2):
       F(E1;b1,b2) = -b1 E1 + logsumexp_j [ lnrho(E1,E2_j) - b2 E2_j ]
     i.e. the transverse (S2) fan is properly summed -- this is what 1D LLR could
     not do. The transition is where F has two equal-height maxima (frozen E1 vs
     disordered E1); scan the chosen coupling for the equal-height point.

Usage: llr2d_reconstruct.py run.out [--plot out.png] [--b1 2 3 4]
"""
import sys, argparse, numpy as np
from numpy.linalg import lstsq

def logsumexp(x):
    x=np.asarray(x,float); m=x.max(); return m+np.log(np.sum(np.exp(x-m)))

def parse(fname):
    cells=[]; geo={}
    for l in open(fname):
        if l.startswith("LLR2D"):
            for tok in l.replace("|"," ").split():
                if "=" in tok:
                    k,v=tok.split("=",1)
                    try: geo[k]=float(v)
                    except: pass
            # E1=[a,b] captured separately
            import re
            m=re.search(r"E1=\[([-\d.eE]+),([-\d.eE]+)\]", l)
            if m: geo["E1top"]=float(m.group(1)); geo["E1bot"]=float(m.group(2))
        if l.startswith("ANE2:"):
            t=l.split()
            cells.append([float(t[1]),float(t[2]),float(t[3]),float(t[4])])
    return np.array(cells), geo

def grid_index(cells, geo):
    """Assign (i,j): i from E1, j from E2 offset off the ridge."""
    E1top=geo["E1top"]; step1=geo["step1"]; nperp=int(geo["nperp"]); step2=geo["step2"]
    c0=geo["c0"]; c1=geo["c1"]; c2=geo["c2"]
    idx=[]
    for E1,E2,a1,a2 in cells:
        i=int(round((E1top-E1)/step1))
        ridge=c0+c1*E1+c2*E1*E1
        j=int(round((E2-ridge)/step2))+nperp
        idx.append((i,j))
    return idx

def integrate(cells, idx):
    """Least-squares lnrho from gradient over (i,j)-adjacent edges."""
    n=len(cells); pos={ij:k for k,ij in enumerate(idx)}
    rows=[]; rhs=[]
    for k,(i,j) in enumerate(idx):
        for (ni,nj) in [(i+1,j),(i,j+1)]:           # forward neighbours only (each edge once)
            if (ni,nj) in pos:
                k2=pos[(ni,nj)]
                x1=cells[k][:2]; x2=cells[k2][:2]
                g1=cells[k][2:]; g2=cells[k2][2:]
                d=0.5*(g1+g2)@(x2-x1)                # lnrho[k2]-lnrho[k] ~ d
                r=np.zeros(n); r[k2]=1; r[k]=-1; rows.append(r); rhs.append(d)
    # gauge fix: lnrho at the most-disordered reached cell = 0
    kref=int(np.argmax(cells[:,0]))                  # max E1 (least negative)
    r=np.zeros(n); r[kref]=1; rows.append(r); rhs.append(0.0)
    A=np.array(rows); b=np.array(rhs)
    lnrho,_,_,_=lstsq(A,b,rcond=None)
    return lnrho

def Fmarg(cells, lnrho, b1, b2):
    """F(E1)= -b1 E1 + logsumexp_j[lnrho - b2 E2] over transverse cells at each E1."""
    E1=cells[:,0]; E2=cells[:,1]
    levels=sorted(set(np.round(E1,2)))
    Es=[]; F=[]
    for e in levels:
        m=np.round(E1,2)==e
        Es.append(e); F.append(-b1*e + logsumexp(lnrho[m]-b2*E2[m]))
    return np.array(Es), np.array(F)

def peakdiff(Es,F):
    """height(frozen peak) - height(disordered peak); None if not bimodal.
    Es ascending (frozen=most negative=index 0)."""
    j=1+int(np.argmin(F[1:-1]))
    if not (F[j]<F[j-1] and F[j]<F[j+1]): return None
    return F[:j+1].max()-F[j:].max()

def betac(cells,lnrho, fix, val, lo, hi, npts=300):
    """Maxwell equal-height locator (clean, but needs a resolved barrier)."""
    xs=np.linspace(lo,hi,npts); xk=[]; dk=[]
    for x in xs:
        b1 = val if fix=='b1' else x
        b2 = val if fix=='b2' else x
        Es,F=Fmarg(cells,lnrho,b1,b2); d=peakdiff(Es,F)
        if d is not None: xk.append(x); dk.append(d)
    if len(xk)<3: return None
    sl,ic=np.polyfit(xk,dk,1); return -ic/sl

def S1mean(cells,lnrho,b1,b2):
    E1=cells[:,0];E2=cells[:,1]
    lw=lnrho-b1*E1-b2*E2; lw-=lw.max(); w=np.exp(lw)
    return np.sum(E1*w)/np.sum(w)

def betac_steep(cells,lnrho, fix, val, lo, hi, npts=500):
    """Steepest drop of <S1> over the full 2D reweight -- robust at small volume
    where the barrier is shallow (matches the 1D code's locator)."""
    xs=np.linspace(lo,hi,npts)
    O=np.array([S1mean(cells,lnrho, val if fix=='b1' else x, val if fix=='b2' else x) for x in xs])
    g=np.gradient(O,xs); return xs[int(np.argmax(np.abs(g)))]

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("file"); ap.add_argument("--plot",default=None)
    args=ap.parse_args()
    cells,geo=parse(args.file)
    idx=grid_index(cells,geo)
    lnrho=integrate(cells,idx)
    print(f"cells: {len(cells)}   E1 {cells[:,0].min():.0f}..{cells[:,0].max():.0f}   E2 {cells[:,1].min():.0f}..{cells[:,1].max():.0f}")
    b1c=betac_steep(cells,lnrho,'b2',0.0,0.3,3.0)
    print(f"beta2=0  -> beta1_c = {b1c:.3f}   (MC: ~1.25)")
    for b1,truth in [(2.0,-0.23),(3.0,-0.55),(4.0,None)]:
        b2c=betac_steep(cells,lnrho,'b1',b1,-1.6,0.5)
        print(f"beta1={b1} -> beta2_c = {b2c:.3f}   (MC: {truth})")
    if args.plot:
        import matplotlib;matplotlib.use("Agg");import matplotlib.pyplot as plt
        plt.figure(figsize=(7,5.5))
        sc=plt.scatter(cells[:,0],cells[:,1],c=lnrho,s=40,cmap="viridis")
        plt.colorbar(sc,label="ln rho (arb const)")
        plt.xlabel("S1");plt.ylabel("S2");plt.title("2D-LLR reconstructed ln rho(S1,S2)")
        plt.tight_layout();plt.savefig(args.plot,dpi=120);print("plot ->",args.plot)

if __name__=="__main__":
    main()
