/*
 * Discrete Yang-Mills -- Metropolis MC that also saves SU(3) configurations.
 *
 * Identical Monte Carlo to dym-mod-metro.cpp (Symanzik-improved action:
 * S = beta0 + beta1*sum_plaq ReTr U + beta2*sum_rect ReTr U), but after each
 * measurement it writes the configuration as SU(3) matrices in Kentucky (NERSC)
 * byte order for Wilson-flow scale setting.
 *
 * Requires a group file with the 18-real (3x3 complex) defining-rep matrices
 * appended after the multiplication table (e.g. groups/mys1080-v4). The output
 * directory/prefix is the last command-line argument.
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
#include <fstream>
#include "timer.h"

/* Defaults for the optional CLI sweep counts (override via argv 10/11/12). */
#define K_DEFAULT       200    // decorrelation sweeps between saved configurations
#define N_DEFAULT       1000   // number of configurations to save
#define NTHERM_DEFAULT  0      // thermalization sweeps (set >0 before saving physics)
#define NHIT            20     // Metropolis hits per link per visit

typedef unsigned int uint;

extern "C" void getpos(unsigned idx, unsigned int *pos);   // defined in lattice.c
void update();

std::default_random_engine *rnd;   // one PRNG per site (reproducible under OpenMP)
std::vector<group_t> smallgroup;   // proposal pool (elements near the identity)
char *site_parity = NULL;          // checkerboard parity per site
unsigned long hit = 0, acc = 0;

const char *outprefix = "./";      // output directory/prefix (CLI arg)

/* Defining-rep matrix of each group element: data[part + 2*(c + 3*r)]. */
struct su3 { double data[18]; double& operator()(int r, int c, int part) { return data[part+2*(c+3*r)]; } };
std::vector<su3> groupsu3;

/* Byte-swap an array of doubles (little<->big endian for the NERSC format). */
void switchend(unsigned char *buffer, int length)
{
  unsigned char *pos = buffer, save[8];
  for(int j=0; j<length; j++) {
    for(int i=0; i<8; i++) save[i] = pos[i];
    for(int i=7; i>-1; i--, pos++) *pos = save[i];
  }
}

/* Write the current configuration (as SU(3) matrices) to <outprefix>lat-...num%04d. */
void savekyconfig(int n)
{
  double *buf = new double[(size_t)V*4*18];
  uint pos[4];
  for(uint i=0; i<V; ++i)
  {
    getpos(i, pos);
    for(int d=0; d<D; ++d)
    {
      su3 &m = groupsu3[a[i*D+d]];
      for(int part=0; part<2; ++part) for(int c=0; c<3; ++c) for(int r=0; r<3; ++r)
        buf[pos[1]+Nx*(pos[2]+Nx*(pos[3]+Nx*(pos[0]+Nt*(part+2*(r+3*(c+3*((d+3)%D)))))))] = m(r, c, part);
    }
  }
  switchend((unsigned char*)buf, V*4*18);

  char name[512];
  snprintf(name, sizeof(name), "%slat-beta1%.3f-beta2%.3f-nt%02d-nx%02d-num%04d",
           outprefix, beta1, beta2, Nt, Nx, n);
  FILE* f = fopen(name, "w");
  if(!f) { fprintf(stderr, "cannot open %s for writing\n", name); abort(); }
  fwrite(buf, V*4, sizeof(su3), f);
  fclose(f);
  delete [] buf;
}

/* --- Monte Carlo (mirrors dym-mod-metro.cpp) --- */

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

/* The 6(D-1) staples of the 1x2 rectangle term touching link (i,d). */
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

/* One checkerboard Metropolis sweep with NHIT hits per link. */
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
      if (site_parity[i] != parity) continue;

      group_t staples[2 * (D - 1)];   get_staples(i, d, staples);
      group_t rect[6 * (D - 1)];      get_rect_loops(i, d, rect);

      group_t b = a[i*D+d];
      double r1 = 0, m1 = 0;
      for (int j = 0; j < nstaple; ++j) r1 += ReTr[mult[b][staples[j]]];
      for (int k = 0; k < nrect;   ++k) m1 += ReTr[mult[b][rect[k]]];

      for (int h = 0; h < NHIT; ++h)
      {
        group_t bnew = mult[b][ smallgroup[ smallgroup.size() * rand01(rnd[i]) ] ];

        double r1new = 0, m1new = 0;
        for (int j = 0; j < nstaple; ++j) r1new += ReTr[mult[bnew][staples[j]]];
        for (int k = 0; k < nrect;   ++k) m1new += ReTr[mult[bnew][rect[k]]];

        double oldact = beta0 + beta1*r1    + beta2*m1;
        double newact = beta0 + beta1*r1new + beta2*m1new;
        double probrat = (newact <= oldact) ? 1.0 : exp(-(newact - oldact));

        ++lhit;
        if (probrat > rand01(rnd[i])) {
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
	if (argc < 9) {
		fprintf(stderr, "usage: %s group D Nt Nx beta0 beta1 beta2 seed [outprefix] [K] [N] [Ntherm]\n", argv[0]);
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
	if (argc >= 10) outprefix = argv[9];
	/* Optional sweep counts; default to K_DEFAULT / N_DEFAULT / NTHERM_DEFAULT. */
	unsigned K      = (argc > 10) ? (unsigned)atoi(argv[10]) : K_DEFAULT;
	unsigned N      = (argc > 11) ? (unsigned)atoi(argv[11]) : N_DEFAULT;
	unsigned NTHERM = (argc > 12) ? (unsigned)atoi(argv[12]) : NTHERM_DEFAULT;
	printf("PARAMS(grp,D,Nt,Nx,beta0,beta1,beta2,seed): %s %d %d %d %e %e %e %d\n",
	       groupfilename, D, Nt, Nx, beta0, beta1, beta2, iseed);
	printf("RUN(K,N,Ntherm,NHIT): %u %u %u %d\n", K, N, NTHERM, NHIT);
	printf("saving configs to: %s\n", outprefix);

	V = Nt;
	for (unsigned d = 1; d < D; ++d) V *= Nx;
	rnd = new std::default_random_engine[V];
	for (int i = 0; i < V; ++i) rnd[i].seed(iseed + i);

	load_group(groupfilename);

	/* Read the appended defining-rep matrices (skip order, ReTr, ImTr, table). */
	{
		std::ifstream ingroup(groupfilename);
		int num; ingroup >> num;
		double tmp; int itmp;
		for (int i = 0; i < 2*num; ++i) ingroup >> tmp;
		for (int i = 0; i < num*num; ++i) ingroup >> itmp;
		groupsu3.resize(num);
		for (int i = 0; i < num; ++i) {
			double *pmat = groupsu3[i].data;
			for (int j = 0; j < 18; ++j, ++pmat) ingroup >> *pmat;
		}
		if (ingroup.fail()) {
			fprintf(stderr, "group file '%s' has no appended SU(3) matrices; "
			                "cannot save configurations\n", groupfilename);
			abort();
		}
	}

	std::uniform_int_distribution<> randgrp(0, P-1);
	printf("read the group\n"); fflush(stdout);
	printf("id: %i\n", id);
	printf("check stuff: (a*b)*c: %d a*(b*c): %d\n", mult[mult[1][2]][3], mult[1][mult[2][3]]);

	/* Proposal pool: elements at the smallest non-identity ReTr (as in dym-mod-metro). */
	double min_retr = 0, nn_retr = 10;
	for (uint i = 0; i < P; ++i) if (ReTr[i] < min_retr) min_retr = ReTr[i];
	for (uint i = 0; i < P; ++i) if (ReTr[i] > min_retr && ReTr[i] < nn_retr) nn_retr = ReTr[i];
	for (uint i = 0; i < P; ++i) if (ReTr[i] < nn_retr + 1e-6 && ReTr[i] > min_retr) smallgroup.push_back(i);
	printf("small group size: %lu\n", smallgroup.size());

	a = (group_t*) malloc(sizeof(*a) * V * D);
	for (unsigned i = 0; i < V*D; i++) a[i] = randgrp(rnd[0]);
	step(0, 0, 1);

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
		double rep, imp; getpoly(&rep, &imp);
		double simpleplaq = 0;
#pragma omp parallel for reduction(+: simpleplaq)
		for (unsigned int i = 0; i < V; ++i)
		for (unsigned int d1 = 0; d1 < D; ++d1)
		for (unsigned int d2 = d1+1; d2 < D; ++d2)
			simpleplaq += ReTr[plaquette(i, d1, d2)];
		simpleplaq /= V*D*(D-1)/2;
		printf("GMES: %e %e %e %e\n", 999.0, rep, imp, simpleplaq); fflush(stdout);
		printf("ACC: %f\n", ((double)acc)/hit);

		savekyconfig(n);
		tm.stop();
	}

	free(a);
	free(site_parity);
	delete[] rnd;
	return 0;
}
