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
double S1;            // running plaquette action (LLR energy)
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

void update_constrained(double a_tilt, double beta2, double Elo, double Ehi,
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
      double dS1 = r1new - r1, S1new = S1 + dS1;
      bool accept = false;
      if (driveIn && (S1 < Elo || S1 > Ehi)) {
        accept = fabs(S1new - Emid) < fabs(S1 - Emid);
      } else if (S1new >= Elo && S1new <= Ehi) {
        double w = exp(-a_tilt*dS1 - beta2*(m1new - m1));
        accept = (w > rand01(rnd[i]));
      }
      ++lhit;
      if (accept) { a[i*D+d] = bnew; S1 = S1new; b = bnew; r1 = r1new; m1 = m1new; ++lacc; }
    }
  }
  *hitp += lhit; *accp += lacc;
}

/* Drive the current config's S1 into [Elo,Ehi] (greedy). Works when the config
 * starts within ~one window-step of the target (sequential seeding). */
bool seed_window(double beta2, double Elo, double Ehi, double a_guess)
{
  long h=0, ac=0; unsigned drive=0, CAP=20000;
  while ((S1 < Elo || S1 > Ehi) && drive < CAP) {
    update_constrained(a_guess, beta2, Elo, Ehi, true, &h, &ac); ++drive;
  }
  S1 = recompute_S1();
  return (S1 >= Elo && S1 <= Ehi);
}

/* RM-solve a_n in window centred at Ewin (config assumed already in-window). */
double rm_solve(double beta2, double Ewin, double hw, double a0,
                unsigned K, unsigned NRM)
{
  double Elo = Ewin - hw, Ehi = Ewin + hw, W = 2*hw, a = a0;
  for (unsigned m = 0; m < NRM; ++m) {
    double sum = 0; long hit = 0, acc = 0;
    for (unsigned k = 0; k < K; ++k) {
      update_constrained(a, beta2, Elo, Ehi, false, &hit, &acc);
      sum += S1;
    }
    double g = sum/K - Ewin;
    a += (12.0/(W*W*(m+1))) * g;
    double drift = fabs(S1 - recompute_S1());
    if (drift > g_maxdrift) g_maxdrift = drift;
    S1 = recompute_S1();
  }
  return a;
}

int main(int argc, char *argv[])
{
  if (argc < 12) {
    fprintf(stderr, "usage: %s group D Nt Nx beta2 Etop Ebot step hw a0 seed [K] [NRM] [R]\n", argv[0]);
    return 1;
  }
  const char *groupfilename = argv[1];
  D  = atoi(argv[2]);  Nt = atoi(argv[3]);  Nx = atoi(argv[4]);
  beta2        = atof(argv[5]);
  double Etop  = atof(argv[6]);   // highest (disordered) window centre
  double Ebot  = atof(argv[7]);   // lowest (frozen) window centre
  double stepE = atof(argv[8]);   // spacing between window centres (>0)
  double hw    = atof(argv[9]);   // window half-width (overlap if hw > stepE/2)
  double a0    = atof(argv[10]);
  int iseed    = atoi(argv[11]);
  unsigned K   = (argc > 12) ? (unsigned)atoi(argv[12]) : 50;
  unsigned NRM = (argc > 13) ? (unsigned)atoi(argv[13]) : 80;
  unsigned R   = (argc > 14) ? (unsigned)atoi(argv[14]) : 1;

  int M = (int)((Etop - Ebot)/stepE + 0.5) + 1;
  printf("LLR-m3(grp,D,Nt,Nx,beta2): %s %d %d %d %g | Etop=%g Ebot=%g step=%g hw=%g a0=%g seed=%d K=%u NRM=%u R=%u | %d windows\n",
         groupfilename, D, Nt, Nx, beta2, Etop, Ebot, stepE, hw, a0, iseed, K, NRM, R, M);

  V = Nt; for (unsigned d = 1; d < D; ++d) V *= Nx;
  rnd = new std::default_random_engine[V];
  load_group(groupfilename);
  static std::uniform_int_distribution<> rg(0, P-1); randgrp = &rg;
  double min_retr = 0, nn_retr = 10;
  for (uint i = 0; i < P; ++i) if (ReTr[i] < min_retr) min_retr = ReTr[i];
  for (uint i = 0; i < P; ++i) if (ReTr[i] > min_retr && ReTr[i] < nn_retr) nn_retr = ReTr[i];
  for (uint i = 0; i < P; ++i) if (ReTr[i] < nn_retr + 1e-6 && ReTr[i] > min_retr) smallgroup.push_back(i);
  a = (group_t*) malloc(sizeof(*a) * V * D);
  for (unsigned i = 0; i < V; ++i) rnd[i].seed(iseed + i);
  step(0, 0, 1);

  std::vector<std::vector<double>> an(M);   // an[n] = a_n samples over repeats

  for (unsigned r = 0; r < R; ++r) {
    for (unsigned i = 0; i < V; ++i) rnd[i].seed(iseed + r*1000003u + i);
    for (unsigned i = 0; i < V*D; ++i) a[i] = (*randgrp)(rnd[0]);   // hot start
    S1 = recompute_S1();
    double a_prev = a0;
    int reached = 0;
    for (int n = 0; n < M; ++n) {
      double Ewin = Etop - n*stepE, Elo = Ewin - hw, Ehi = Ewin + hw;
      if (!seed_window(beta2, Elo, Ehi, a_prev)) {
        printf("repeat %u: stop at window %d (Ewin=%.1f): seeding failed (S1=%.1f)\n",
               r, n, Ewin, S1);
        break;
      }
      double a_n = rm_solve(beta2, Ewin, hw, a_prev, K, NRM);
      an[n].push_back(a_n);
      a_prev = a_n;     // smooth initial guess for the next (lower) window
      ++reached;
    }
    printf("repeat %u: reached %d / %d windows\n", r, reached, M);
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
  printf("max |incr - recompute| S1 : %.3e\n", g_maxdrift);

  free(a); delete[] rnd;
  return 0;
}
