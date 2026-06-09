/*
 * LLR for discrete Yang-Mills -- MILESTONES 1-2
 * ============================================
 * Incremental energy tracking + single-window constrained Metropolis (m1) and a
 * Robbins-Monro solver for the LLR coefficient a_n (m2).
 *
 * Energy:  E = S1 = sum over all plaquettes of ReTr U_plaq  (the beta1-conjugate
 * action). Constrained ensemble at window [Ewin +/- delta/2]:
 *     weight ~ exp(-a*S1 - beta2*S2),  restricted to S1 in [Elo, Ehi].
 *
 * RM iteration solves a_n so that <<dE>> = <S1> - Ewin = 0 (the reweighted
 * distribution becomes flat in the window, i.e. a_n = d ln rho / dS1):
 *     a^(m+1) = a^(m) + (12 / (delta^2 (m+1))) * (<S1>_a - Ewin).
 * NOTE the + sign: with the negated-trace convention S1 decreases as a grows,
 * so g = <S1> - Ewin has dg/da < 0 and the stable RM step adds g (verified by
 * g -> 0 in the trajectory output). [Paper Eq. 6 uses - with E >= 0.]
 *
 * Updates are SERIAL over sites (the window constraint is on the global S1).
 * Real LLR parallelism is over independent (interval x repeat) tasks, later.
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
double S1;   // running plaquette action (the LLR energy)

/* From-scratch S1 = sum_{plaquettes} ReTr U_plaq (each plaquette once). */
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

/* Plaquette staples touching link (i,d) (same as dym-mod-metro). */
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

/* Rectangle staples touching link (i,d) (same as dym-mod-metro). */
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

/* One serial sweep of the constrained update, maintaining S1 incrementally.
 * driveIn=true: greedily move S1 toward the window centre (init only). */
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

      double dS1 = r1new - r1;
      double S1new = S1 + dS1;
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

static std::uniform_int_distribution<> *randgrp;

/* (Re)initialize a hot config and drive S1 into the window. Returns false if
 * the greedy drive-in could not reach the window (deep windows need sequential
 * seeding -- milestone 3). */
bool init_into_window(double a_tilt, double beta2, double Elo, double Ehi)
{
  for (unsigned i = 0; i < V*D; i++) a[i] = (*randgrp)(rnd[0]);
  S1 = recompute_S1();
  long h = 0, ac = 0; unsigned drive = 0, CAP = 20000;
  while ((S1 < Elo || S1 > Ehi) && drive < CAP) {
    update_constrained(a_tilt, beta2, Elo, Ehi, true, &h, &ac); ++drive;
  }
  S1 = recompute_S1();
  return (S1 >= Elo && S1 <= Ehi);
}

int main(int argc, char *argv[])
{
  if (argc < 10) {
    fprintf(stderr, "usage: %s group D Nt Nx beta2 Ewin delta a0 seed [K] [NRM] [R]\n", argv[0]);
    return 1;
  }
  const char *groupfilename = argv[1];
  D  = atoi(argv[2]);  Nt = atoi(argv[3]);  Nx = atoi(argv[4]);
  beta2        = atof(argv[5]);
  double Ewin  = atof(argv[6]);
  double delta = atof(argv[7]);
  double a0    = atof(argv[8]);
  int iseed    = atoi(argv[9]);
  unsigned K   = (argc > 10) ? (unsigned)atoi(argv[10]) : 200;   // sweeps / RM step
  unsigned NRM = (argc > 11) ? (unsigned)atoi(argv[11]) : 400;   // RM iterations
  unsigned R   = (argc > 12) ? (unsigned)atoi(argv[12]) : 1;     // repeats
  double Elo = Ewin - delta/2, Ehi = Ewin + delta/2;
  printf("LLR-m2(grp,D,Nt,Nx,beta2,Ewin,delta,a0,seed,K,NRM,R): %s %d %d %d %g %g %g %g %d %u %u %u\n",
         groupfilename, D, Nt, Nx, beta2, Ewin, delta, a0, iseed, K, NRM, R);
  printf("window S1 in [%.3f, %.3f]\n", Elo, Ehi);

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
  printf("small group size: %lu\n\n", smallgroup.size());

  const double coef0 = 12.0/(delta*delta);
  std::vector<double> a_finals;
  double maxdrift = 0;

  for (unsigned r = 0; r < R; ++r)
  {
    for (unsigned i = 0; i < V; ++i) rnd[i].seed(iseed + r*1000003u + i);
    if (!init_into_window(a0, beta2, Elo, Ehi)) {
      printf("repeat %u: drive-in failed to reach window (S1=%.3f) -- skipping\n", r, S1);
      continue;
    }

    double a = a0;
    printf("repeat %u: RM trajectory (m, a, <S1>, g=<S1>-Ewin, acc)\n", r);
    for (unsigned m = 0; m < NRM; ++m)
    {
      double sum = 0; long hit = 0, acc = 0;
      for (unsigned k = 0; k < K; ++k) {
        update_constrained(a, beta2, Elo, Ehi, false, &hit, &acc);
        sum += S1;
      }
      double meanS1 = sum/K, g = meanS1 - Ewin;
      a += (coef0/(m+1)) * g;                    // Robbins-Monro step (see header)

      double drift = fabs(S1 - recompute_S1());  // ongoing m1 sanity check
      if (drift > maxdrift) maxdrift = drift;
      S1 = recompute_S1();

      if (m < 6 || (m+1)%50 == 0 || m == NRM-1)
        printf("  %4u  a=%.5f  <S1>=%9.3f  g=%+8.3f  acc=%.3f\n",
               m, a, meanS1, g, hit ? (double)acc/hit : 0.0);
    }
    a_finals.push_back(a);
    printf("repeat %u: converged a_n = %.5f\n\n", r, a);
  }

  printf("=== result ===\n");
  if (a_finals.empty()) { printf("no successful repeats\n"); }
  else {
    double mean = 0; for (double x : a_finals) mean += x; mean /= a_finals.size();
    double var = 0; for (double x : a_finals) var += (x-mean)*(x-mean);
    double sd = a_finals.size()>1 ? sqrt(var/(a_finals.size()-1)) : 0.0;
    printf("window Ewin=%.3f delta=%.3f : a_n = %.5f +/- %.5f  (%lu repeats)\n",
           Ewin, delta, mean, sd, a_finals.size());
  }
  printf("max |incr - recompute| S1 over RM steps : %.3e\n", maxdrift);

  free(a); delete[] rnd;
  return 0;
}
