/*
 * 2D LLR for discrete Yang-Mills -- joint density of states rho(S1,S2)
 * =================================================================
 * The 1D LLR (llr.cpp) windows ONE action and lets the other ride a fixed
 * background coupling. When that background is strong it PINS the spectator and
 * breaks the S1-S2 coupling that drives the joint first-order jump, so a single
 * axis mislocates the freezing line (axis-1 hooks at -beta2; axis-2 is too
 * shallow at large beta1 -- both confirmed against direct Metropolis).
 *
 * Here we window BOTH actions at once: constrain (S1,S2) into a 2D cell
 * [E1+/-hw1] x [E2+/-hw2] and Robbins-Monro solve the gradient
 *     a = (a1,a2) = (d lnrho/dS1, d lnrho/dS2)
 * via two coupled recursions that pin <S1>->E1 and <S2>->E2 simultaneously.
 * The tilt weight inside the cell is exp(-a1*S1 - a2*S2) (negated-trace
 * convention, as in 1D). No spectator, no pinning -> the joint jump is
 * representable.
 *
 * The physical (S1,S2) support is a thin ribbon (S1,S2 strongly correlated:
 * ordered -> both saturated, disordered -> both ~0), so we tile a BAND around
 * the curved ridge E2(E1) rather than the full rectangle; off-ribbon cells fail
 * to seed and are skipped.
 *
 * Output: "ANE2: E1 E2 a1 a2 sd1 sd2 nrep" lines -> 2D reconstruction.
 */

extern "C" {
#include "group.h"
#include "lattice.h"
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <random>
#include "timer.h"

#define NHIT 20
typedef unsigned int uint;

extern "C" void getpos(unsigned idx, unsigned int *pos);

std::default_random_engine *rnd;
std::vector<group_t> smallgroup;
static std::uniform_int_distribution<> *randgrp;
double S1;                 // running plaquette action (sum_plaq ReTr)
double S2;                 // running rectangle action (sum_rect ReTr)
double g_drift1 = 0, g_drift2 = 0;

/* ---- physics helpers: identical to llr.cpp (same staples/rectangles) ---- */
double recompute_S1()
{
  double s = 0;
#pragma omp parallel for reduction(+: s)
  for (unsigned i = 0; i < V; ++i)
    for (unsigned d1 = 0; d1 < D; ++d1)
      for (unsigned d2 = d1+1; d2 < D; ++d2)
        s += ReTr[plaquette(i, d1, d2)];
  return s;
}

void get_staples(unsigned int i, unsigned int d, group_t* st)
{
  unsigned int i1 = step(i, d, 1); int k=0;
  for(unsigned int d1=0; d1<D; ++d1) if(d1 != d)
  {
    group_t g = a[i1*D+d1];
    g = mult[g][inv[a[step(i,d1,1)*D+d]]];
    g = mult[g][inv[a[i*D+d1]]];
    st[k++] = g;
    g = inv[a[step(i1,d1,-1)*D+d1]];
    int i2 = step(i,d1,-1);
    g = mult[g][inv[a[i2*D+d]]];
    g = mult[g][a[i2*D+d1]];
    st[k++] = g;
  }
}

void get_rect_loops(unsigned int i, unsigned int d, group_t* st)
{
  unsigned int i1 = step(i, d, 1); int k=0;
  for(unsigned int d1 = 0; d1 < D; ++d1) if(d1 != d)
  {
    group_t g = id;
    g = mult[g][a[step(i,d,1)*D+d1]];
    int i2 = step(step(i,d,1),d1,1);
    g = mult[g][a[i2*D+d1]];
    g = mult[g][inv[a[step(step(i2,d1,1),d,-1)*D+d]]];
    g = mult[g][inv[a[step(i,d1,1)*D+d1]]];
    g = mult[g][inv[a[i*D+d1]]];
    st[k++] = g;

    g = id;
    g = mult[g][a[step(i,d,1)*D+d]];
    i2 = step(step(i,d,1),d,1);
    g = mult[g][a[i2*D+d1]];
    g = mult[g][inv[a[step(step(i2,d1,1),d,-1)*D+d]]];
    g = mult[g][inv[a[step(i,d1,1)*D+d]]];
    g = mult[g][inv[a[i*D+d1]]];
    st[k++] = g;

    g = id;
    g = mult[g][a[step(i,d,1)*D+d1]];
    i2 = step(step(i,d,1),d1,1);
    g = mult[g][inv[a[step(i2,d,-1)*D+d]]];
    g = mult[g][inv[a[step(step(i2,d,-1),d,-1)*D+d]]];
    g = mult[g][inv[a[step(i,d,-1)*D+d1]]];
    g = mult[g][a[step(i,d,-1)*D+d]];
    st[k++] = g;

    g = id;
    g = mult[g][inv[a[step(step(i,d,1),d1,-1)*D+d1]]];
    i2 = step(step(i,d,1),d1,-1);
    g = mult[g][inv[a[step(i2,d1,-1)*D+d1]]];
    g = mult[g][inv[a[step(step(i2,d1,-1),d,-1)*D+d]]];
    g = mult[g][a[step(step(i,d1,-1),d1,-1)*D+d1]];
    g = mult[g][a[step(i,d1,-1)*D+d1]];
    st[k++] = g;

    g = id;
    g = mult[g][a[step(i,d,1)*D+d]];
    i2 = step(step(i,d,1),d,1);
    g = mult[g][inv[a[step(i2,d1,-1)*D+d1]]];
    g = mult[g][inv[a[step(step(i2,d1,-1),d,-1)*D+d]]];
    g = mult[g][inv[a[step(i,d1,-1)*D+d]]];
    g = mult[g][a[step(i,d1,-1)*D+d1]];
    st[k++] = g;

    g = id;
    g = mult[g][inv[a[step(step(i,d,1),d1,-1)*D+d1]]];
    i2 = step(step(i,d,1),d1,-1);
    g = mult[g][inv[a[step(i2,d,-1)*D+d]]];
    g = mult[g][inv[a[step(step(i2,d,-1),d,-1)*D+d]]];
    g = mult[g][a[step(step(i,d,-1),d1,-1)*D+d1]];
    g = mult[g][a[step(i,d,-1)*D+d]];
    st[k++] = g;
  }
}

double RECT_MULT = 6.0;
double sum_m1()
{
  double s = 0; const int nrect = 6*(D-1);
  for (unsigned i = 0; i < V; ++i)
  for (unsigned d = 0; d < D; ++d) {
    group_t rect[6*(D-1)]; get_rect_loops(i, d, rect);
    group_t b = a[i*D+d];
    for (int k = 0; k < nrect; ++k) s += ReTr[mult[b][rect[k]]];
  }
  return s;
}
double recompute_S2() { return sum_m1()/RECT_MULT; }

void calibrate_rect_mult()
{
  const int nrect = 6*(D-1);
  unsigned i = V/2, d = 0;
  group_t rect[6*(D-1)]; get_rect_loops(i, d, rect);
  group_t b = a[i*D+d];
  double m1_old = 0; for (int k=0;k<nrect;++k) m1_old += ReTr[mult[b][rect[k]]];
  double T_old = sum_m1();
  group_t bnew = mult[b][ smallgroup[ smallgroup.size()/2 ] ];
  double m1_new = 0; for (int k=0;k<nrect;++k) m1_new += ReTr[mult[bnew][rect[k]]];
  a[i*D+d] = bnew;
  double T_new = sum_m1();
  a[i*D+d] = b;
  RECT_MULT = (T_new - T_old)/(m1_new - m1_old);
  printf("RECT_MULT calibrated: %.6f\n", RECT_MULT);
}

inline void resync_E() { S1 = recompute_S1(); S2 = recompute_S2(); }

/* ---- 2D constrained update: tilt (a1,a2) on (S1,S2), box [E1lo,E1hi]x[E2lo,E2hi] ----
 * driveIn: greedily reduce the normalized distance to the cell centre so the
 * config enters the 2D box; otherwise Metropolis with weight exp(-a1 dS1 - a2 dS2)
 * restricted to the box. Serial over sites (S1,S2 are global running sums). */
void update_constrained2d(double a1, double a2,
                          double E1lo, double E1hi, double E2lo, double E2hi,
                          bool driveIn, long *hitp, long *accp)
{
  std::uniform_real_distribution<> rand01(0., 1.);
  const int nstaple = 2*(D-1), nrect = 6*(D-1);
  const double E1mid = 0.5*(E1lo+E1hi), E2mid = 0.5*(E2lo+E2hi);
  const double hw1 = 0.5*(E1hi-E1lo), hw2 = 0.5*(E2hi-E2lo);
  long lhit = 0, lacc = 0;

  for (unsigned d = 0; d < D; d++)
  for (uint i = 0; i < V; i++)
  {
    group_t staples[2*(D-1)]; get_staples(i, d, staples);
    group_t rect[6*(D-1)];    get_rect_loops(i, d, rect);
    group_t b = a[i*D+d];
    double r1 = 0, m1 = 0;
    for (int j=0;j<nstaple;++j) r1 += ReTr[mult[b][staples[j]]];
    for (int k=0;k<nrect;  ++k) m1 += ReTr[mult[b][rect[k]]];

    for (int h = 0; h < NHIT; ++h)
    {
      group_t bnew = mult[b][ smallgroup[ smallgroup.size()*rand01(rnd[i]) ] ];
      double r1new = 0, m1new = 0;
      for (int j=0;j<nstaple;++j) r1new += ReTr[mult[bnew][staples[j]]];
      for (int k=0;k<nrect;  ++k) m1new += ReTr[mult[bnew][rect[k]]];
      double dS1 = r1new - r1, dS2 = m1new - m1;
      double S1n = S1 + dS1, S2n = S2 + dS2;
      bool inbox_now = (S1>=E1lo && S1<=E1hi && S2>=E2lo && S2<=E2hi);
      bool inbox_new = (S1n>=E1lo && S1n<=E1hi && S2n>=E2lo && S2n<=E2hi);
      bool accept = false;
      if (driveIn && !inbox_now) {
        double dc = hypot((S1 -E1mid)/hw1, (S2 -E2mid)/hw2);
        double dn = hypot((S1n-E1mid)/hw1, (S2n-E2mid)/hw2);
        accept = (dn < dc);                       // move toward the cell centre
      } else if (inbox_new) {
        double w = exp(-a1*dS1 - a2*dS2);
        accept = (w > rand01(rnd[i]));
      }
      ++lhit;
      if (accept) { a[i*D+d]=bnew; S1=S1n; S2=S2n; b=bnew; r1=r1new; m1=m1new; ++lacc; }
    }
  }
  *hitp += lhit; *accp += lacc;
}

/* Drive (S1,S2) into the 2D box (greedy). Returns false if it cannot. */
bool seed_cell(double E1lo,double E1hi,double E2lo,double E2hi,
               double a1,double a2, unsigned CAP)
{
  long h=0, ac=0; unsigned drive=0;
  auto inbox = [&]{ return S1>=E1lo&&S1<=E1hi&&S2>=E2lo&&S2<=E2hi; };
  while (!inbox() && drive < CAP) {
    update_constrained2d(a1,a2,E1lo,E1hi,E2lo,E2hi,true,&h,&ac); ++drive;
  }
  resync_E();
  return inbox();
}

/* 2D Robbins-Monro: two coupled recursions pin <S1>->E1 and <S2>->E2.
 * a1 += 12/(W1^2 (m+1)) (<S1>-E1) ;  a2 += 12/(W2^2 (m+1)) (<S2>-E2). */
void rm_solve2d(double E1,double E2,double hw1,double hw2,
                double a1_0,double a2_0,unsigned K,unsigned NRM,
                double *a1_out,double *a2_out)
{
  double E1lo=E1-hw1,E1hi=E1+hw1,E2lo=E2-hw2,E2hi=E2+hw2;
  double W1=2*hw1, W2=2*hw2, a1=a1_0, a2=a2_0;
  for (unsigned m = 0; m < NRM; ++m) {
    double s1=0,s2=0; long hit=0,acc=0;
    for (unsigned k = 0; k < K; ++k) {
      update_constrained2d(a1,a2,E1lo,E1hi,E2lo,E2hi,false,&hit,&acc);
      s1 += S1; s2 += S2;
    }
    a1 += (12.0/(W1*W1*(m+1))) * (s1/K - E1);
    a2 += (12.0/(W2*W2*(m+1))) * (s2/K - E2);
    double d1 = fabs(S1 - recompute_S1()), d2 = fabs(S2 - recompute_S2());
    if (d1 > g_drift1) g_drift1 = d1;
    if (d2 > g_drift2) g_drift2 = d2;
    resync_E();
  }
  *a1_out = a1; *a2_out = a2;
}

/* ---- config I/O (for process-based parallel-over-cells: "level 2") ---- */
void save_cfg(const char *p) {
  FILE *f = fopen(p, "wb");
  if (f) { fwrite(a, sizeof(*a), (size_t)V*D, f); fclose(f); }
}
bool load_cfg(const char *p) {
  FILE *f = fopen(p, "rb");
  if (!f) return false;
  size_t n = fread(a, sizeof(*a), (size_t)V*D, f); fclose(f);
  return n == (size_t)V*D;
}

/* Common lattice setup: group tables, proposal pool, hot config, RECT_MULT. */
void setup_lattice(const char *grp, int iseed)
{
  rnd = new std::default_random_engine[V];
  load_group(grp);
  static std::uniform_int_distribution<> rg(0, P-1); randgrp = &rg;
  double min_retr = 0, nn_retr = 10;
  for (uint i=0;i<P;++i) if (ReTr[i]<min_retr) min_retr=ReTr[i];
  for (uint i=0;i<P;++i) if (ReTr[i]>min_retr && ReTr[i]<nn_retr) nn_retr=ReTr[i];
  for (uint i=0;i<P;++i) if (ReTr[i]<nn_retr+1e-6 && ReTr[i]>min_retr) smallgroup.push_back(i);
  a = (group_t*) malloc(sizeof(*a)*V*D);
  for (unsigned i=0;i<V;++i) rnd[i].seed(iseed+i);
  for (unsigned i=0;i<V*D;++i) a[i] = (*randgrp)(rnd[0]);
  step(0,0,1);
  calibrate_rect_mult();
  if (fabs(RECT_MULT-6.0) > 1e-3)
    printf("WARNING: RECT_MULT=%.3f != 6 -- rectangles wrap; need lattice extent >= 3\n", RECT_MULT);
}

/* ---- cell grid: a band around the curved ridge E2 = c0+c1*E1+c2*E1^2 ---- */
struct Cell { double E1, E2; std::vector<double> a1s, a2s; bool reached=false; };

int mode_full(int argc, char *argv[])
{
  if (argc < 18) {
    fprintf(stderr,
      "usage: %s group D Nt Nx E1top E1bot step1 hw1 c0 c1 c2 step2 hw2 nperp a0 seed [K] [NRM] [R]\n"
      "  windows BOTH S1 and S2; tiles a band around the curved ridge\n"
      "  E2 = c0 + c1*E1 + c2*E1^2, offset by (j-nperp)*step2 (j=0..2*nperp).\n", argv[0]);
    return 1;
  }
  const char *grp = argv[1];
  D=atoi(argv[2]); Nt=atoi(argv[3]); Nx=atoi(argv[4]);
  double E1top=atof(argv[5]), E1bot=atof(argv[6]), step1=atof(argv[7]), hw1=atof(argv[8]);
  double c0=atof(argv[9]), c1=atof(argv[10]), c2=atof(argv[11]);
  double step2=atof(argv[12]), hw2=atof(argv[13]);
  int nperp=atoi(argv[14]);
  double a0=atof(argv[15]);
  int iseed=atoi(argv[16]);
  unsigned K   = (argc>17)?(unsigned)atoi(argv[17]):40;
  unsigned NRM = (argc>18)?(unsigned)atoi(argv[18]):70;
  unsigned R   = (argc>19)?(unsigned)atoi(argv[19]):1;
  #define RIDGE(E1) (c0 + c1*(E1) + c2*(E1)*(E1))

  int N1 = (int)((E1top - E1bot)/step1 + 0.5) + 1;
  int N2 = 2*nperp + 1;
  printf("LLR2D(grp,D,Nt,Nx): %s %d %d %d | E1=[%g,%g] step1=%g hw1=%g | ridge c0=%g c1=%g c2=%g step2=%g hw2=%g nperp=%d | a0=%g seed=%d K=%u NRM=%u R=%u | %dx%d cells\n",
         grp,D,Nt,Nx,E1top,E1bot,step1,hw1,c0,c1,c2,step2,hw2,nperp,a0,iseed,K,NRM,R,N1,N2);
  V = Nt; for (unsigned d=1; d<D; ++d) V *= Nx;
  printf("NPLAQ: %u\nNRECT: %u\n", V*D*(D-1)/2, V*D*(D-1));
  setup_lattice(grp, iseed);

  // grid[i][j]: E1 = E1top - i*step1 ; E2 = RIDGE(E1) + (j-nperp)*step2
  std::vector<std::vector<Cell>> grid(N1, std::vector<Cell>(N2));
  for (int i=0;i<N1;++i) for (int j=0;j<N2;++j) {
    double E1 = E1top - i*step1;
    grid[i][j].E1 = E1;
    grid[i][j].E2 = RIDGE(E1) + (j-nperp)*step2;
  }

  const unsigned CAP = 4000;
  for (unsigned r=0;r<R;++r) {
    for (unsigned i=0;i<V;++i) rnd[i].seed(iseed + r*1000003u + i);
    // ---- cold pass: frozen start, walk ridge from frozen (i=N1-1) up to disordered (i=0)
    for (unsigned i=0;i<V*D;++i) a[i]=id;
    resync_E();
    int creach=0;
    for (int s=0;s<N1;++s) {
      int i = N1-1-s;                         // frozen -> disordered
      // seed the ridge centre cell (j=nperp) from current config, solve, then fan out transverse
      Cell &c0 = grid[i][nperp];
      if (!seed_cell(c0.E1-hw1,c0.E1+hw1,c0.E2-hw2,c0.E2+hw2,a0,a0,CAP)) continue;
      double b1,b2; rm_solve2d(c0.E1,c0.E2,hw1,hw2,a0,a0,K,NRM,&b1,&b2);
      c0.a1s.push_back(b1); c0.a2s.push_back(b2); c0.reached=true; ++creach;
      // snapshot the ridge config to re-seed transverse cells
      std::vector<group_t> ridgecfg(a, a+V*D);
      for (int dir=-1; dir<=1; dir+=2)
        for (int t=1;t<=nperp;++t) {
          int j = nperp + dir*t;
          Cell &c = grid[i][j];
          // restore ridge config (perpendicular drive from the ridge each time)
          for (unsigned q=0;q<V*D;++q) a[q]=ridgecfg[q];
          resync_E();
          if (!seed_cell(c.E1-hw1,c.E1+hw1,c.E2-hw2,c.E2+hw2,b1,b2,CAP)) break;
          double q1,q2; rm_solve2d(c.E1,c.E2,hw1,hw2,b1,b2,K,NRM,&q1,&q2);
          c.a1s.push_back(q1); c.a2s.push_back(q2); c.reached=true;
        }
      // restore ridge config for the next ridge step
      for (unsigned q=0;q<V*D;++q) a[q]=ridgecfg[q];
      resync_E();
    }
    // ---- hot pass: disordered start, walk ridge from disordered (i=0) to frozen
    for (unsigned i=0;i<V*D;++i) a[i]=(*randgrp)(rnd[0]);
    resync_E();
    int hreach=0;
    for (int i=0;i<N1;++i) {
      Cell &c0 = grid[i][nperp];
      if (!seed_cell(c0.E1-hw1,c0.E1+hw1,c0.E2-hw2,c0.E2+hw2,a0,a0,CAP)) continue;
      double b1,b2; rm_solve2d(c0.E1,c0.E2,hw1,hw2,a0,a0,K,NRM,&b1,&b2);
      c0.a1s.push_back(b1); c0.a2s.push_back(b2); c0.reached=true; ++hreach;
      std::vector<group_t> ridgecfg(a, a+V*D);
      for (int dir=-1; dir<=1; dir+=2)
        for (int t=1;t<=nperp;++t) {
          int j = nperp + dir*t;
          Cell &c = grid[i][j];
          for (unsigned q=0;q<V*D;++q) a[q]=ridgecfg[q];
          resync_E();
          if (!seed_cell(c.E1-hw1,c.E1+hw1,c.E2-hw2,c.E2+hw2,b1,b2,CAP)) break;
          double q1,q2; rm_solve2d(c.E1,c.E2,hw1,hw2,b1,b2,K,NRM,&q1,&q2);
          c.a1s.push_back(q1); c.a2s.push_back(q2); c.reached=true;
        }
      for (unsigned q=0;q<V*D;++q) a[q]=ridgecfg[q];
      resync_E();
    }
    printf("repeat %u: cold ridge %d/%d, hot ridge %d/%d\n", r, creach, N1, hreach, N1);
  }

  printf("\n# E1 E2 a1 a2 sd1 sd2 nrep\n");
  int nreached=0;
  for (int i=0;i<N1;++i) for (int j=0;j<N2;++j) {
    Cell &c = grid[i][j];
    if (!c.reached) continue;
    nreached++;
    double m1=0,m2=0; int n=c.a1s.size();
    for (int k=0;k<n;++k){m1+=c.a1s[k];m2+=c.a2s[k];} m1/=n; m2/=n;
    double v1=0,v2=0; for (int k=0;k<n;++k){v1+=(c.a1s[k]-m1)*(c.a1s[k]-m1);v2+=(c.a2s[k]-m2)*(c.a2s[k]-m2);}
    double sd1=n>1?sqrt(v1/(n-1)):0, sd2=n>1?sqrt(v2/(n-1)):0;
    printf("ANE2: %9.2f %9.2f %10.6f %10.6f %9.5f %9.5f %d\n", c.E1,c.E2,m1,m2,sd1,sd2,n);
  }
  printf("reached %d/%d cells\n", nreached, N1*N2);
  printf("max |incr-recompute|: S1 %.3e  S2 %.3e\n", g_drift1, g_drift2);
  free(a); delete[] rnd;
  return 0;
}

/* ============================ LEVEL-2 PARALLELISM ============================
 * Split the single-run cost across processes: a cheap SEQUENTIAL seeding pass
 * drives a config into every cell and dumps it to disk (+ a manifest); then N
 * independent SOLVE workers each load a disjoint range of cells and run the
 * (expensive) RM in parallel. Equivalent to mode_full -- rm_solve2d starts from
 * the drive-in config either way -- but the RM now parallelizes over cells.
 * ========================================================================== */

/* seed: drive-in every reachable cell (bidirectional ridge walk), save its
 * config, and emit a manifest "CELL: idx i j E1 E2 hw1 hw2 path" on stdout. */
int mode_seed(int argc, char *argv[])
{
  if (argc < 18) {
    fprintf(stderr, "usage: %s seed group D Nt Nx E1top E1bot step1 hw1 c0 c1 c2 step2 hw2 nperp a0 seed cfgdir\n", argv[0]);
    return 1;
  }
  const char *grp=argv[1];
  D=atoi(argv[2]); Nt=atoi(argv[3]); Nx=atoi(argv[4]);
  double E1top=atof(argv[5]),E1bot=atof(argv[6]),step1=atof(argv[7]),hw1=atof(argv[8]);
  double c0=atof(argv[9]),c1=atof(argv[10]),c2=atof(argv[11]);
  double step2=atof(argv[12]),hw2=atof(argv[13]);
  int nperp=atoi(argv[14]);
  double a0=atof(argv[15]);
  int iseed=atoi(argv[16]);
  const char *cfgdir=argv[17];
  int N1=(int)((E1top-E1bot)/step1+0.5)+1, N2=2*nperp+1;
  V=Nt; for(unsigned d=1;d<D;++d) V*=Nx;
  printf("LLR2D(seed): %s %d %d %d | E1=[%g,%g] step1=%g hw1=%g | ridge c0=%g c1=%g c2=%g step2=%g hw2=%g nperp=%d\n",
         grp,D,Nt,Nx,E1top,E1bot,step1,hw1,c0,c1,c2,step2,hw2,nperp);
  setup_lattice(grp, iseed);

  std::vector<std::vector<Cell>> grid(N1, std::vector<Cell>(N2));
  std::vector<std::vector<bool>> saved(N1, std::vector<bool>(N2,false));
  for(int i=0;i<N1;++i) for(int j=0;j<N2;++j){
    double E1=E1top-i*step1; grid[i][j].E1=E1; grid[i][j].E2=c0+c1*E1+c2*E1*E1+(j-nperp)*step2;
  }
  const unsigned CAP=4000; int idx=0; char path[1024];
  auto try_save=[&](int i,int j){
    Cell&c=grid[i][j];
    if(saved[i][j]) return;
    if(!(S1>=c.E1-hw1&&S1<=c.E1+hw1&&S2>=c.E2-hw2&&S2<=c.E2+hw2)) return;
    snprintf(path,sizeof(path),"%s/cell_%d_%d.cfg",cfgdir,i,j);
    save_cfg(path); saved[i][j]=true;
    printf("CELL: %d %d %d %.2f %.2f %.2f %.2f %s\n", idx++,i,j,c.E1,c.E2,hw1,hw2,path);
  };
  auto fan=[&](int i){
    std::vector<group_t> ridgecfg(a,a+V*D);
    for(int dir=-1;dir<=1;dir+=2) for(int t=1;t<=nperp;++t){
      int j=nperp+dir*t; Cell&c=grid[i][j];
      for(unsigned q=0;q<V*D;++q) a[q]=ridgecfg[q]; resync_E();
      if(!seed_cell(c.E1-hw1,c.E1+hw1,c.E2-hw2,c.E2+hw2,a0,a0,CAP)) break;
      try_save(i,j);
    }
    for(unsigned q=0;q<V*D;++q) a[q]=ridgecfg[q]; resync_E();
  };
  // cold pass (frozen -> disordered)
  for(unsigned q=0;q<V*D;++q) a[q]=id; resync_E();
  for(int s=0;s<N1;++s){ int i=N1-1-s; Cell&cc=grid[i][nperp];
    if(!seed_cell(cc.E1-hw1,cc.E1+hw1,cc.E2-hw2,cc.E2+hw2,a0,a0,CAP)) continue;
    try_save(i,nperp); fan(i); }
  // hot pass (disordered -> frozen)
  for(unsigned q=0;q<V*D;++q) a[q]=(*randgrp)(rnd[0]); resync_E();
  for(int i=0;i<N1;++i){ Cell&cc=grid[i][nperp];
    if(!seed_cell(cc.E1-hw1,cc.E1+hw1,cc.E2-hw2,cc.E2+hw2,a0,a0,CAP)) continue;
    try_save(i,nperp); fan(i); }
  fprintf(stderr,"SEED: saved %d cells to %s\n", idx, cfgdir);
  free(a); delete[] rnd; return 0;
}

/* solve: load cells [lo,hi) from the manifest and RM-solve each (parallel-safe:
 * independent processes, disjoint cell ranges). Emits ANE2 lines. */
int mode_solve(int argc, char *argv[])
{
  if (argc < 12) {
    fprintf(stderr, "usage: %s solve group D Nt Nx manifest lo hi a0 K NRM seed\n", argv[0]);
    return 1;
  }
  const char *grp=argv[1];
  D=atoi(argv[2]); Nt=atoi(argv[3]); Nx=atoi(argv[4]);
  const char *manifest=argv[5];
  int lo=atoi(argv[6]), hi=atoi(argv[7]);
  double a0=atof(argv[8]);
  unsigned K=(unsigned)atoi(argv[9]), NRM=(unsigned)atoi(argv[10]);
  int iseed=atoi(argv[11]);
  V=Nt; for(unsigned d=1;d<D;++d) V*=Nx;
  setup_lattice(grp, iseed);
  FILE *mf=fopen(manifest,"r");
  if(!mf){ fprintf(stderr,"solve: cannot open manifest %s\n",manifest); return 1; }
  char line[2048]; int done=0;
  while(fgets(line,sizeof(line),mf)){
    if(strncmp(line,"CELL:",5)) continue;
    int idx,i,j; double E1,E2,hw1,hw2; char path[1024];
    if(sscanf(line,"CELL: %d %d %d %lf %lf %lf %lf %1023s",&idx,&i,&j,&E1,&E2,&hw1,&hw2,path)!=8) continue;
    if(idx<lo||idx>=hi) continue;
    if(!load_cfg(path)){ fprintf(stderr,"solve: missing cfg %s\n",path); continue; }
    resync_E();
    double b1,b2; rm_solve2d(E1,E2,hw1,hw2,a0,a0,K,NRM,&b1,&b2);
    printf("ANE2: %9.2f %9.2f %10.6f %10.6f %9.5f %9.5f %d\n", E1,E2,b1,b2,0.0,0.0,1);
    fflush(stdout); ++done;
  }
  fclose(mf);
  fprintf(stderr,"SOLVE[%d,%d): %d cells, drift S1 %.3e S2 %.3e\n", lo,hi,done,g_drift1,g_drift2);
  free(a); delete[] rnd; return 0;
}

int main(int argc, char *argv[])
{
  if (argc > 1 && strcmp(argv[1],"seed")==0)  return mode_seed (argc-1, argv+1);
  if (argc > 1 && strcmp(argv[1],"solve")==0) return mode_solve(argc-1, argv+1);
  return mode_full(argc, argv);   // default: single-process full grid
}
