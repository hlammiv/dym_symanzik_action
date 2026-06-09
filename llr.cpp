/*
 * LLR for discrete Yang-Mills -- MILESTONE 1
 * =========================================
 * Incremental energy tracking + a single-window constrained Metropolis update.
 *
 * Energy:  E = S1 = sum over all plaquettes of ReTr U_plaq  (the beta1-conjugate
 * action). The rectangle term S2 (beta2-conjugate) enters only as a fixed
 * background in the sampling weight.
 *
 * Constrained ensemble at window [Elo, Ehi] = [Ewin - delta/2, Ewin + delta/2]:
 *     weight ~ exp(-a*S1 - beta2*S2),  restricted to S1 in [Elo, Ehi]
 * where `a` is the LLR tilt that replaces beta1 (Robbins-Monro will solve for it
 * later; here it is a fixed input).
 *
 * A local link move changes S1 by exactly dS1 = r1new - r1 (the change of the
 * plaquettes touching that link), so S1 is maintained incrementally instead of
 * recomputed globally. This file's job is to verify two things:
 *   (1) the incremental S1 matches a from-scratch global recompute, and
 *   (2) the constrained update confines S1 to the window -- including a window
 *       placed inside the first-order coexistence gap that unconstrained runs
 *       never reach.
 *
 * Updates here are SERIAL over sites: the window constraint is on the global S1,
 * so a checkerboard parallel sweep is not conditionally independent. (Real LLR
 * parallelism is over independent (interval x repeat) tasks, added later.)
 */

extern "C" {
#include "group.h"
#include "lattice.h"
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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
 * driveIn=true: when S1 is outside the window, greedily move it toward the
 * window centre (initialization only); once inside, the hard window applies. */
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
        /* funnel toward the window (no detailed balance; init only) */
        accept = fabs(S1new - Emid) < fabs(S1 - Emid);
      } else if (S1new >= Elo && S1new <= Ehi) {
        /* hard window + LLR tilt a on S1, beta2 background on S2 */
        double w = exp(-a_tilt*dS1 - beta2*(m1new - m1));
        accept = (w > rand01(rnd[i]));
      }
      ++lhit;
      if (accept) { a[i*D+d] = bnew; S1 = S1new; b = bnew; r1 = r1new; m1 = m1new; ++lacc; }
    }
  }
  *hitp += lhit; *accp += lacc;
}

int main(int argc, char *argv[])
{
  if (argc < 10) {
    fprintf(stderr, "usage: %s group D Nt Nx beta2 Ewin delta a seed [K] [N]\n", argv[0]);
    return 1;
  }
  const char *groupfilename = argv[1];
  D  = atoi(argv[2]);
  Nt = atoi(argv[3]);
  Nx = atoi(argv[4]);
  beta2       = atof(argv[5]);
  double Ewin = atof(argv[6]);
  double delta= atof(argv[7]);
  double a_tilt = atof(argv[8]);
  int iseed   = atoi(argv[9]);
  unsigned K = (argc > 10) ? (unsigned)atoi(argv[10]) : 100;
  unsigned N = (argc > 11) ? (unsigned)atoi(argv[11]) : 100;
  double Elo = Ewin - delta/2, Ehi = Ewin + delta/2;
  printf("LLR-m1(grp,D,Nt,Nx,beta2,Ewin,delta,a,seed): %s %d %d %d %e %e %e %e %d\n",
         groupfilename, D, Nt, Nx, beta2, Ewin, delta, a_tilt, iseed);
  printf("window S1 in [%.3f, %.3f]\n", Elo, Ehi);

  V = Nt; for (unsigned d = 1; d < D; ++d) V *= Nx;
  rnd = new std::default_random_engine[V];
  for (int i = 0; i < V; ++i) rnd[i].seed(iseed + i);
  load_group(groupfilename);
  std::uniform_int_distribution<> randgrp(0, P-1);

  double min_retr = 0, nn_retr = 10;
  for (uint i = 0; i < P; ++i) if (ReTr[i] < min_retr) min_retr = ReTr[i];
  for (uint i = 0; i < P; ++i) if (ReTr[i] > min_retr && ReTr[i] < nn_retr) nn_retr = ReTr[i];
  for (uint i = 0; i < P; ++i) if (ReTr[i] < nn_retr + 1e-6 && ReTr[i] > min_retr) smallgroup.push_back(i);
  printf("small group size: %lu\n", smallgroup.size());

  a = (group_t*) malloc(sizeof(*a) * V * D);
  for (unsigned i = 0; i < V*D; i++) a[i] = randgrp(rnd[0]);
  step(0, 0, 1);
  S1 = recompute_S1();
  printf("hot-start S1 = %.3f  (plaquettes = %u)\n", S1, V*D*(D-1)/2);

  /* Drive S1 into the window. */
  long h0 = 0, a0 = 0;
  unsigned drive = 0, DRIVE_CAP = 20000;
  while ((S1 < Elo || S1 > Ehi) && drive < DRIVE_CAP) {
    update_constrained(a_tilt, beta2, Elo, Ehi, true, &h0, &a0);
    ++drive;
  }
  if (S1 < Elo || S1 > Ehi)
    printf("WARNING: failed to reach window in %u sweeps (S1=%.3f)\n", DRIVE_CAP, S1);
  else
    printf("reached window after %u drive-in sweeps (S1=%.3f)\n", drive, S1);
  S1 = recompute_S1();   /* resync after drive-in */

  /* Constrained sampling + verification. */
  long hit = 0, acc = 0;
  double sumS1 = 0, minS1 = 1e30, maxS1 = -1e30, maxdrift = 0;
  int outside = 0;
  for (unsigned n = 0; n < N; n++) {
    for (unsigned k = 0; k < K; k++) update_constrained(a_tilt, beta2, Elo, Ehi, false, &hit, &acc);
    double S1rec = recompute_S1();
    double drift = fabs(S1 - S1rec);
    if (drift > maxdrift) maxdrift = drift;
    S1 = S1rec;                                  /* resync (report drift first) */
    sumS1 += S1; if (S1 < minS1) minS1 = S1; if (S1 > maxS1) maxS1 = S1;
    if (S1 < Elo - 1e-6 || S1 > Ehi + 1e-6) ++outside;
  }

  printf("\n=== verification ===\n");
  printf("window           : [%.3f, %.3f]\n", Elo, Ehi);
  printf("<S1>             : %.3f\n", sumS1/N);
  printf("S1 range sampled : [%.3f, %.3f]  (must be within window)\n", minS1, maxS1);
  printf("measurements out : %d / %u\n", outside, N);
  printf("max |incr - recompute| over K-sweep blocks : %.3e\n", maxdrift);
  printf("acceptance       : %.3f\n", hit ? (double)acc/hit : 0.0);

  free(a); delete[] rnd;
  return 0;
}
