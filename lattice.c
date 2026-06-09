/*
 * Lattice geometry and gauge-invariant loops.
 *
 * Sites are indexed 0..V-1 with coordinate 0 = Euclidean time (extent Nt) and
 * coordinates 1..D-1 = space (extent Nx each); both directions are periodic.
 * The gauge field `a[i*D + mu]` holds the group-element index of link (i, mu).
 */
#include "lattice.h"

#include <stdlib.h>
#include <stdio.h>

unsigned D, Nt, Nx, V;
double beta0, beta1, beta2, beta3, beta4;
group_t *a;
unsigned int *nn = NULL;   // precomputed neighbour table: nn[i*2D + 2*d + (s>0)]

void init_nn();

/* Index of the site one step from i along direction d (s = +1 or -1).
 * Uses a precomputed neighbour table (built lazily on first call). */
unsigned step(unsigned i, unsigned d, int s) {
  if(abs(s) != 1) {printf("Not implemented ... exiting\n"); abort();}
  if(nn == NULL) init_nn();
  return nn[i*2*D+d*2+(1+s)/2];
}

/* Walk s steps along direction d (s may exceed 1 in magnitude). */
unsigned bigstep(unsigned i, unsigned d, int s) {
	unsigned idx=0;
	do{
		if(s>0){i=step(i,d,1);idx++;}
		if(s<0){i=step(i,d,-1);idx--;}
	}while(idx != s);
	return i;
}

/* coords -> linear index (coord 0 is time, period Nt; the rest space, period Nx). */
unsigned int getidx(unsigned int *pos)
{
  int idx = 0;
  for(int d=D-1; d>0; --d) idx = (pos[d]+Nx)%Nx + Nx*idx;
  return (pos[0] + Nt)%Nt + Nt*idx;
}

/* linear index -> coords. */
void getpos(unsigned idx, unsigned int *pos)
{
  pos[0] = idx % Nt; idx /= Nt;
  for(int d=1; d<D; ++d) { pos[d] = idx % Nx; idx /= Nx; }
}

/* Build the neighbour table once: for every site and direction, the indices of
 * the two adjacent sites (minus and plus). */
void init_nn()
{
  printf("Initializing nearest neighbours [%d:%d]... \n", Nt, Nx);
  nn = malloc(V*2*D*sizeof(unsigned int));
  unsigned int pos[D];
  for(int i=0; i<V; ++i)
  {
    getpos(i, pos);
    for(int d=0; d<D; ++d)
    {
      pos[d]--;
      nn[2*D*i+2*d+0] = getidx(pos);   // neighbour at -1
      pos[d] += 2;
      nn[2*D*i+2*d+1] = getidx(pos);   // neighbour at +1
      pos[d]--;
    }
  }
}

/* 1x1 Wilson plaquette in the (d1,d2) plane at site i:
 * U_{d1}(i) U_{d2}(i+d1) U_{d1}^{-1}(i+d2) U_{d2}^{-1}(i). */
group_t plaquette(unsigned i, unsigned d1, unsigned d2) {
	group_t g = id;
	g = mult[g][a[i*D+d1]];
	g = mult[g][a[step(i,d1,1)*D+d2]];
	g = mult[g][inv[a[step(i,d2,1)*D+d1]]];
	g = mult[g][inv[a[i*D+d2]]];
	return g;
}

/* Product of temporal links around the (periodic) time direction through i. */
group_t polyakov(unsigned i) {
	group_t r = id;
	unsigned i0 = i;
	do {
		r = mult[r][a[i*D+0]];
		i = step(i,0,1);
	} while (i0 != i);
	return r;
}

/* Volume-averaged Polyakov loop (real and imaginary parts of its trace),
 * averaged over the V/Nt spatial sites of the t=0 timeslice. */
void getpoly(double *re, double *im)
{
  double lre=0, lim=0;
#pragma omp parallel for reduction(+: lre, lim)
  for(int i=0; i<V; ++i)
  {
    unsigned int pos[D];
    getpos(i, pos);
    if(pos[0] != 0) continue;
    group_t p = polyakov(i);
    lre += ReTr[p];
    lim += ImTr[p];
  }
  *re = lre/(V/Nt);
  *im = lim/(V/Nt);
}

/* nd1 x nd2 rectangular Wilson loop in the (d1,d2) plane starting at site i. */
group_t wilson(unsigned i, unsigned d1, unsigned d2, unsigned nd1, unsigned nd2)
{
	group_t g = id;
	for( int j=0; j<nd1; ++j) { g = mult[g][a[i*D+d1]];      i=step(i,d1,1); }
	for( int j=0; j<nd2; ++j) { g = mult[g][a[i*D+d2]];      i=step(i,d2,1); }
	for( int j=0; j<nd1; ++j) { i=step(i,d1,-1); g = mult[g][inv[a[i*D+d1]]]; }
	for( int j=0; j<nd2; ++j) { i=step(i,d2,-1); g = mult[g][inv[a[i*D+d2]]]; }
	return g;
}
