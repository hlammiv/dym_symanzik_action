/*
 * LLR for discrete Yang-Mills -- MILESTONES 1-3
 * ============================================
 * m1: incremental energy (S1) tracking + single-window constrained Metropolis.
 * m2: Robbins-Monro solver for the LLR coefficient a_n in a window.
 * m3: tile the whole action range into overlapping windows and SEQUENTIALLY
 *     SEED each window from the previous one's config, so deep (frozen) windows
 *     that a hot-start drive-in cannot reach become reachable -- producing the
 *     full a_n(E) curve.
 *
 * Energy E = S1 = sum_plaq ReTr U_plaq. Constrained ensemble in window
 * [Ewin +/- hw]: weight ~ exp(-a*S1 - beta2*S2) restricted to S1 in [Elo,Ehi].
 * RM step (negated-trace convention, see m2): a += 12/(W^2 (m+1)) * (<S1>-Ewin),
 * W = window width = 2*hw. a_n = d ln rho / dS1.
 *
 * Output: "ANE: Ewin a_n sd nrep" lines (the input to milestone-4 reconstruction).
 * Updates are serial over sites (the window constraint is global).
 */

extern "C" {
#include "group.h"
#include "lattice.h"
}

#include <stdio.h>
#include <stdlib.h>
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
double S1;            // running plaquette action (sum_plaq ReTr)
double S2;            // running rectangle action (sum_rect ReTr)
int    AXIS = 1;      // 1: window S1 (beta1-driven); 2: window S2 (beta2-driven)
double BG  = 0;       // fixed background coupling (beta2 if AXIS=1, beta1 if AXIS=2)
double g_maxdrift = 0;

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

double RECT_MULT = 6.0;   // how many times each rectangle is counted in sum_m1
                          // (calibrated at startup; see calibrate_rect_mult)

/* sum over all links of m1(link) = the rectangle staple sum. */
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
/* Rectangle action S2 = sum over distinct rectangles ReTr = sum_m1 / RECT_MULT. */
double recompute_S2() { return sum_m1()/RECT_MULT; }

/* Determine RECT_MULT = d(sum_m1)/d(m1 of one link) under a single link flip,
 * so that the incremental dS2 = dm1 matches d(recompute_S2). */
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
  a[i*D+d] = b;                                  // restore
  double mu = (T_new - T_old)/(m1_new - m1_old);
  RECT_MULT = mu;
  printf("RECT_MULT calibrated: %.6f\n", mu);
}
inline double Ecur()           { return (AXIS == 1) ? S1 : S2; }   // windowed energy
inline double recompute_Ewin() { return (AXIS == 1) ? recompute_S1() : recompute_S2(); }
inline void   resync_E()       { S1 = recompute_S1(); S2 = recompute_S2(); }

/* Constrained update windowing the chosen action (AXIS): tilt a on the windowed
 * action, fixed coupling BG on the other; both S1 and S2 tracked incrementally. */
void update_constrained(double a_tilt, double Elo, double Ehi,
                        bool driveIn, long *hitp, long *accp)
{
  std::uniform_real_distribution<> rand01(0., 1.);
  const int nstaple = 2*(D-1), nrect = 6*(D-1);
  const double Emid = 0.5*(Elo + Ehi);
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
      double dEw  = (AXIS == 1) ? dS1 : dS2;     // change of the windowed action
      double dEbg = (AXIS == 1) ? dS2 : dS1;     // change of the background action
      double Ec   = (AXIS == 1) ? S1  : S2;
      double Enew = Ec + dEw;
      bool accept = false;
      if (driveIn && (Ec < Elo || Ec > Ehi)) {
        accept = fabs(Enew - Emid) < fabs(Ec - Emid);
      } else if (Enew >= Elo && Enew <= Ehi) {
        double w = exp(-a_tilt*dEw - BG*dEbg);
        accept = (w > rand01(rnd[i]));
      }
      ++lhit;
      if (accept) { a[i*D+d] = bnew; S1 += dS1; S2 += dS2; b = bnew; r1 = r1new; m1 = m1new; ++lacc; }
    }
  }
  *hitp += lhit; *accp += lacc;
}

/* Drive the windowed energy into [Elo,Ehi] (greedy). Works when the config
 * starts within ~one window-step of the target (sequential seeding). */
bool seed_window(double Elo, double Ehi, double a_guess)
{
  long h=0, ac=0; unsigned drive=0, CAP=20000;
  while ((Ecur() < Elo || Ecur() > Ehi) && drive < CAP) {
    update_constrained(a_guess, Elo, Ehi, true, &h, &ac); ++drive;
  }
  resync_E();
  return (Ecur() >= Elo && Ecur() <= Ehi);
}

/* RM-solve a_n in window centred at Ewin (config assumed already in-window). */
double rm_solve(double Ewin, double hw, double a0, unsigned K, unsigned NRM)
{
  double Elo = Ewin - hw, Ehi = Ewin + hw, W = 2*hw, a = a0;
  for (unsigned m = 0; m < NRM; ++m) {
    double sum = 0; long hit = 0, acc = 0;
    for (unsigned k = 0; k < K; ++k) {
      update_constrained(a, Elo, Ehi, false, &hit, &acc);
      sum += Ecur();
    }
    double g = sum/K - Ewin;
    a += (12.0/(W*W*(m+1))) * g;
    double drift = fabs(Ecur() - recompute_Ewin());
    if (drift > g_maxdrift) g_maxdrift = drift;
    resync_E();
  }
  return a;
}

/* One ladder pass. dir=+1: hot start, walk down (Etop->Ebot); dir=-1: cold
 * (frozen) start, walk up (Ebot->Etop). Disordering (the up pass) is entropically
 * easy, so it reaches the frozen/gap windows the down pass stalls before -- the
 * two passes meet in the gap. Pushes converged a_n into an[]; returns #reached. */
int run_pass(int dir, double a0, double Etop, double stepE,
             double hw, unsigned K, unsigned NRM, int M,
             std::vector<std::vector<double>> &an)
{
  if (dir > 0) for (unsigned i = 0; i < V*D; ++i) a[i] = (*randgrp)(rnd[0]);  // hot
  else         for (unsigned i = 0; i < V*D; ++i) a[i] = id;                  // cold (frozen)
  resync_E();
  double a_prev = a0;
  int reached = 0;
  for (int s = 0; s < M; ++s) {
    int n = (dir > 0) ? s : (M-1-s);
    double Ewin = Etop - n*stepE, Elo = Ewin - hw, Ehi = Ewin + hw;
    if (!seed_window(Elo, Ehi, a_prev)) break;           // stalled; rest unreached this pass
    an[n].push_back(rm_solve(Ewin, hw, a_prev, K, NRM));
    a_prev = an[n].back();                                // smooth guess for next window
    ++reached;
  }
  return reached;
}

int main(int argc, char *argv[])
{
  if (argc < 13) {
    fprintf(stderr, "usage: %s group D Nt Nx axis bg Etop Ebot step hw a0 seed [K] [NRM] [R]\n"
                    "  axis=1: window S1 (beta1-driven), bg=beta2 ; axis=2: window S2 (beta2-driven), bg=beta1\n",
            argv[0]);
    return 1;
  }
  const char *groupfilename = argv[1];
  D  = atoi(argv[2]);  Nt = atoi(argv[3]);  Nx = atoi(argv[4]);
  AXIS         = atoi(argv[5]);   // 1: window S1 ; 2: window S2
  BG           = atof(argv[6]);   // fixed coupling (beta2 if axis 1, beta1 if axis 2)
  double Etop  = atof(argv[7]);   // highest (disordered) window centre
  double Ebot  = atof(argv[8]);   // lowest (frozen) window centre
  double stepE = atof(argv[9]);   // spacing between window centres (>0)
  double hw    = atof(argv[10]);  // window half-width (overlap if hw > stepE/2)
  double a0    = atof(argv[11]);
  int iseed    = atoi(argv[12]);
  unsigned K   = (argc > 13) ? (unsigned)atoi(argv[13]) : 50;
  unsigned NRM = (argc > 14) ? (unsigned)atoi(argv[14]) : 80;
  unsigned R   = (argc > 15) ? (unsigned)atoi(argv[15]) : 1;

  int M = (int)((Etop - Ebot)/stepE + 0.5) + 1;
  printf("LLR(grp,D,Nt,Nx,axis,bg): %s %d %d %d %d %g | Etop=%g Ebot=%g step=%g hw=%g a0=%g seed=%d K=%u NRM=%u R=%u | %d windows\n",
         groupfilename, D, Nt, Nx, AXIS, BG, Etop, Ebot, stepE, hw, a0, iseed, K, NRM, R, M);

  V = Nt; for (unsigned d = 1; d < D; ++d) V *= Nx;
  printf("NPLAQ: %u\n", V*D*(D-1)/2);     // <E1> = S1/NPLAQ/3 + 1
  printf("NRECT: %u\n", V*D*(D-1));       // <E2> = S2/NRECT/3 + 1
  printf("AXIS: %d\n", AXIS);             // which action is windowed (energy = S<AXIS>)
  rnd = new std::default_random_engine[V];
  load_group(groupfilename);
  static std::uniform_int_distribution<> rg(0, P-1); randgrp = &rg;
  double min_retr = 0, nn_retr = 10;
  for (uint i = 0; i < P; ++i) if (ReTr[i] < min_retr) min_retr = ReTr[i];
  for (uint i = 0; i < P; ++i) if (ReTr[i] > min_retr && ReTr[i] < nn_retr) nn_retr = ReTr[i];
  for (uint i = 0; i < P; ++i) if (ReTr[i] < nn_retr + 1e-6 && ReTr[i] > min_retr) smallgroup.push_back(i);
  a = (group_t*) malloc(sizeof(*a) * V * D);
  for (unsigned i = 0; i < V; ++i) rnd[i].seed(iseed + i);
  for (unsigned i = 0; i < V*D; ++i) a[i] = (*randgrp)(rnd[0]);
  step(0, 0, 1);
  calibrate_rect_mult();                  // fixes RECT_MULT so dS2 = d(recompute_S2)
  if (AXIS == 2 && fabs(RECT_MULT - 6.0) > 1e-3)
    printf("WARNING: RECT_MULT=%.3f != 6 -- the 1x2 rectangles wrap on this lattice\n"
           "  (extent too small in some direction); the S2 action is degenerate and the\n"
           "  S2 window is invalid. Use a lattice with extent >= 3 for axis=2.\n", RECT_MULT);

  std::vector<std::vector<double>> an(M);   // an[n] = a_n samples over repeats

  for (unsigned r = 0; r < R; ++r) {
    for (unsigned i = 0; i < V; ++i) rnd[i].seed(iseed + r*1000003u + i);
    int rd = run_pass(+1, a0, Etop, stepE, hw, K, NRM, M, an);   // hot,  down
    int ru = run_pass(-1, a0, Etop, stepE, hw, K, NRM, M, an);   // cold, up
    printf("repeat %u: down reached %d/%d, up reached %d/%d\n", r, rd, M, ru, M);
  }

  printf("\n# Ewin   a_n        sd        nrep\n");
  for (int n = 0; n < M; ++n) {
    double Ewin = Etop - n*stepE;
    if (an[n].empty()) { printf("ANE: %8.2f   --(unreached)\n", Ewin); continue; }
    double mean = 0; for (double x : an[n]) mean += x; mean /= an[n].size();
    double var = 0; for (double x : an[n]) var += (x-mean)*(x-mean);
    double sd = an[n].size()>1 ? sqrt(var/(an[n].size()-1)) : 0.0;
    printf("ANE: %8.2f  %9.5f  %9.5f  %lu\n", Ewin, mean, sd, an[n].size());
  }
  printf("max |incr - recompute| E(axis%d) : %.3e\n", AXIS, g_maxdrift);

  free(a); delete[] rnd;
  return 0;
}
