/*
 * Discrete Yang-Mills -- Metropolis Monte Carlo
 * ============================================
 *
 * Pure lattice gauge theory whose gauge group is a *finite* subgroup of SU(N)
 * (loaded from a multiplication table; see group.c / groups/). Each link
 * U_{i,mu} is stored as an index `a[i*D + mu]` into that table.
 *
 * Action (Symanzik-improved):
 *
 *     S = beta0 + beta1 * sum_plaquettes  ReTr U_plaq
 *              + beta2 * sum_rectangles  ReTr U_rect
 *
 * beta1 multiplies the standard 1x1 Wilson plaquette term; beta2 multiplies the
 * 1x2 rectangle term (tree-level Symanzik improvement). A local link update only
 * needs the *staples* (the open plaquette/rectangle loops that touch that link):
 * the change in S from U -> U' is
 *
 *     dS = beta1 * sum_j ReTr(U'·st_j - U·st_j) + beta2 * sum_k ReTr(U'·re_k - U·re_k)
 *
 * Update: checkerboard (even/odd) Metropolis. For each direction and parity we
 * sweep all sites of that parity in parallel (OpenMP) -- sites of one parity
 * only read links belonging to the other parity, so there are no data races.
 * Each link gets NHIT Metropolis hits; a proposal multiplies the link by a group
 * element close to the identity (`smallgroup`).
 *
 * Per-measurement output is one line tagged "GMES:" (parsed by scripts/):
 *     GMES: <placeholder> <Re Polyakov> <Im Polyakov> <mean plaquette ReTr> <Wilson loops...>
 */

extern "C" {
#include "group.h"
#include "lattice.h"
}

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <random>
#include "timer.h"

#define K       1000   // decorrelation sweeps between measurements
#define N       1000   // number of measurements
#define NTHERM  0      // thermalization sweeps before the first measurement
                       // !! currently 0 -- raise (e.g. to K) for production runs
#define NHIT    20     // Metropolis hits per link per visit

typedef unsigned int uint;

extern "C" void getpos(unsigned idx, unsigned int *pos);   // defined in lattice.c
void update();

std::default_random_engine *rnd;   // one independent PRNG per site (for OpenMP)
std::vector<group_t> smallgroup;   // group elements nearest the identity (proposals)
char *site_parity = NULL;          // site_parity[i] = (sum of coords) & 1 (checkerboard)

unsigned long hit = 0, acc = 0;    // running Metropolis hit/accept counters

/* The 2(D-1) staples of the 1x1 plaquette term touching link (i,d). */
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

/* The 6(D-1) staples of the 1x2 rectangle term touching link (i,d): for each
 * transverse direction d1 there are six rectangles (the link sits in each of the
 * three positions along the long edge, for both orientations). */
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

/* One full Metropolis sweep: every link gets NHIT hits, visited in checkerboard
 * order so each parity can be updated in parallel without races. */
void update()
{
  std::uniform_real_distribution<> rand01(0., 1.);
  const int nstaple = 2 * (D - 1);
  const int nrect   = 6 * (D - 1);

  for (unsigned d = 0; d < D; d++)
  for (int parity = 0; parity < 2; ++parity)
  {
    long lhit = 0, lacc = 0;
#pragma omp parallel for schedule(static) reduction(+: lacc, lhit)
    for (uint i = 0; i < V; i++)
    {
      if (site_parity[i] != parity) continue;        // only this parity's links

      group_t staples[2 * (D - 1)];   get_staples(i, d, staples);
      group_t rect[6 * (D - 1)];      get_rect_loops(i, d, rect);

      /* Action contribution of the current link, cached across hits and only
       * recomputed when a hit is accepted (b changes). */
      group_t b = a[i*D+d];
      double r1 = 0, m1 = 0;
      for (int j = 0; j < nstaple; ++j) r1 += ReTr[mult[b][staples[j]]];
      for (int k = 0; k < nrect;   ++k) m1 += ReTr[mult[b][rect[k]]];

      for (int h = 0; h < NHIT; ++h)
      {
        /* Propose b -> b * g, with g a group element near the identity. */
        group_t bnew = mult[b][ smallgroup[ smallgroup.size() * rand01(rnd[i]) ] ];

        double r1new = 0, m1new = 0;
        for (int j = 0; j < nstaple; ++j) r1new += ReTr[mult[bnew][staples[j]]];
        for (int k = 0; k < nrect;   ++k) m1new += ReTr[mult[bnew][rect[k]]];

        double oldact = beta0 + beta1*r1    + beta2*m1;
        double newact = beta0 + beta1*r1new + beta2*m1new;
        /* Accept with prob min(1, exp(-dS)); skip exp when dS<=0 (auto-accept),
         * but always draw rand01 so the PRNG stream is unchanged. */
        double probrat = (newact <= oldact) ? 1.0 : exp(-(newact - oldact));

        ++lhit;
        if (probrat > rand01(rnd[i])) {                // accept
          a[i*D+d] = bnew;
          b = bnew; r1 = r1new; m1 = m1new;
          ++lacc;
        }
      }
    }
    hit += lhit;
    acc += lacc;
  }
}

int main(int argc, char *argv[]) {
	/* Theory parameters come from the command line. */
	if (argc < 9) {
		fprintf(stderr, "usage: %s group D Nt Nx beta0 beta1 beta2 seed\n", argv[0]);
		return 1;
	}
	const char *groupfilename = argv[1];
	D     = atoi(argv[2]);
	Nt    = atoi(argv[3]);
	Nx    = atoi(argv[4]);
	beta0 = atof(argv[5]);
	beta1 = atof(argv[6]);
	beta2 = atof(argv[7]);
	int iseed = atoi(argv[8]);
	printf("PARAMS(grp,D,Nt,Nx,beta0,beta1,beta2,seed): %s %d %d %d %e %e %e %d\n",
	       groupfilename, D, Nt, Nx, beta0, beta1, beta2, iseed);

	/* Volume = Nt * Nx^(D-1) (integer product; avoids pow() rounding). */
	V = Nt;
	for (unsigned d = 1; d < D; ++d) V *= Nx;

	/* One PRNG per site so the parallel sweep is reproducible regardless of
	 * thread count. */
	rnd = new std::default_random_engine[V];
	for (int i = 0; i < V; ++i) rnd[i].seed(iseed + i);

	/* Load the multiplication table. */
	load_group(groupfilename);
	std::uniform_int_distribution<> randgrp(0, P-1);
	printf("read the group\n"); fflush(stdout);
	printf("id: %i\n", id);
	printf("check stuff: (a*b)*c: %d a*(b*c): %d\n", mult[mult[1][2]][3], mult[1][mult[2][3]]);

	/* Proposal pool: group elements at the smallest non-identity ReTr (i.e.
	 * "closest" to the identity), so single hits make small, ergodic moves. */
	double min_retr = 0;
	double nn_retr  = 10;
	for (uint i = 0; i < P; ++i) if (ReTr[i] < min_retr) min_retr = ReTr[i];
	for (uint i = 0; i < P; ++i) if (ReTr[i] > min_retr && ReTr[i] < nn_retr) nn_retr = ReTr[i];
	for (uint i = 0; i < P; ++i) if (ReTr[i] < nn_retr + 1e-6 && ReTr[i] > min_retr) smallgroup.push_back(i);
	printf("min_retr: %e nn_retr: %e \n", min_retr, nn_retr);
	printf("small group size: %lu\n", smallgroup.size());

	/* Initialize the gauge field to a random ("hot") start. */
	a = (group_t*) malloc(sizeof(*a) * V * D);
	for (unsigned i = 0; i < V*D; i++)
		a[i] = randgrp(rnd[0]);

	step(0, 0, 1);   // prime the nearest-neighbour table

	/* Precompute each site's checkerboard parity once (was recomputed every
	 * sweep, and the old inline version assumed D==4). */
	site_parity = (char*) malloc(V);
	for (unsigned i = 0; i < V; ++i) {
		unsigned pos[D];
		getpos(i, pos);
		int s = 0;
		for (unsigned d = 0; d < D; ++d) s += pos[d];
		site_parity[i] = s & 1;
	}
	printf("init done\n"); fflush(stdout);

	for (unsigned k = 0; k < NTHERM; k++) update();
	printf("thermo done\n"); fflush(stdout);

	for (unsigned n = 0; n < N; n++)
	{
		timer tm;
		tm.start("update");
		for (unsigned k = 0; k < K; k++) update();
		tm.stop();

		tm.start("meas");
		double rep, imp; getpoly(&rep, &imp);          // Polyakov loop (Re, Im)

		/* Mean plaquette ReTr over all V*D(D-1)/2 plaquettes. */
		double simpleplaq = 0;
#pragma omp parallel for reduction(+: simpleplaq)
		for (unsigned int i = 0; i < V; ++i)
		for (unsigned int d1 = 0; d1 < D; ++d1)
		for (unsigned int d2 = d1+1; d2 < D; ++d2)
			simpleplaq += ReTr[plaquette(i, d1, d2)];
		simpleplaq /= V*D*(D-1)/2;

		/* Nx x Nt grid of (spatial x temporal) Wilson loops. */
		double wloop[Nx][Nt];
		for (unsigned int k = 0; k < Nx; ++k)
		for (unsigned int l = 0; l < Nt; ++l)
		{
			double lwloop = 0;
#pragma omp parallel for reduction(+: lwloop)
			for (unsigned int i = 0; i < V; ++i)
			for (unsigned int d1 = 1; d1 < D; ++d1)
				lwloop += ReTr[wilson(i, d1, 0, k+1, l+1)];
			/* average over the V*(D-1) loops summed (V sites x (D-1) spatial
			 * directions), each a (k+1)x(l+1) spatial-temporal Wilson loop). */
			wloop[k][l] = lwloop / (V*(D-1));
		}

		/* GMES line: field [1] is a legacy placeholder; scripts read [4]=plaquette. */
		printf("GMES: %e %e %e %e", 999.0, rep, imp, simpleplaq);
		for (unsigned int k = 0; k < Nx; ++k)
		for (unsigned int l = 0; l < Nt; ++l)
			printf(" %e", wloop[k][l]);
		printf("\n");
		fflush(stdout);
		printf("ACC: %f\n", ((double)acc)/hit);
		tm.stop();
	}

	free(a);
	free(site_parity);
	delete[] rnd;
	return 0;
}
