#include "lattice.h"

#include <stdlib.h>
#include <stdio.h>

unsigned D, Nt, Nx, V;
double beta0, beta1, beta2, beta3, beta4;
group_t *a;
unsigned int *nn = NULL;

void init_nn();
unsigned step(unsigned i, unsigned d, int s) {
  if(abs(s) != 1) {printf("Not implemented ... exiting\n"); abort();}
  if(nn == NULL) init_nn();
  return nn[i*2*D+d*2+(1+s)/2];
#if 0
	unsigned under = 1;
	for (unsigned i = 0; i < D-d-1; i++) under *= L;
	return (under*L)*(i/(under*L)) + (i+under*s+abs(s)*under*L)%(under*L);
#endif
}

unsigned bigstep(unsigned i, unsigned d, int s) {
	unsigned idx=0;

	do{
		if(s>0){i=step(i,d,1);idx++;}
		if(s<0){i=step(i,d,-1);idx--;}
//		printf("%d %d %d %d\n",i,d,s,idx);
	}while(idx != s);
//	printf("%d %d %d %d\n",i,d,s,idx);
	return i;
}

// 0 is the temp dir up to Nt and all else is spatial
unsigned int getidx(unsigned int *pos)
{
  int idx = 0;
  for(int d=D-1; d>0; --d) idx = (pos[d]+Nx)%Nx + Nx*idx; 
//printf("pos: %d %d %d %d -> %d (Nt:%d idx:%d Nx:%d)\n", pos[0], pos[1], pos[2], pos[3], (pos[0] + Nt)%Nt + Nt*idx, Nt, idx, Nx);
  return (pos[0] + Nt)%Nt + Nt*idx;
}

void getpos(unsigned idx, unsigned int *pos)
{
  pos[0] = idx % Nt; idx /= Nt;
  for(int d=1; d<D; ++d) { pos[d] = idx % Nx; idx /= Nx; }
}

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
      nn[2*D*i+2*d+0] = getidx(pos);
      pos[d] += 2;
      nn[2*D*i+2*d+1] = getidx(pos);
      pos[d]--;    
    }
    printf("% 4d: ", i);
    for(int j=0; j<2*D; ++j) printf("% 4d ", nn[2*D*i+j]);
    printf("\n");
  }
}

group_t plaquette(unsigned i, unsigned d1, unsigned d2) {
	group_t g = id;
//printf("(i,d1,d2): %d %d %d\n", i, d1, d2);
//printf("[0123]: %d %d %d %d\n", a[i*D+d1], a[step(i,d1,1)*D+d2], inv[a[step(i,d2,1)*D+d1]], inv[a[i*D+d2]]);
	g = mult[g][a[i*D+d1]];
	g = mult[g][a[step(i,d1,1)*D+d2]];
	g = mult[g][inv[a[step(i,d2,1)*D+d1]]];
	g = mult[g][inv[a[i*D+d2]]];
	return g;
}

group_t polyakov(unsigned i) {
	group_t r = id;
	unsigned i0 = i;
	do {
		r = mult[r][a[i*D+0]];
		i = step(i,0,1);
	} while (i0 != i);
	return r;
}

void getpoly(double *re, double *im)
{
  //*re = *im = 0.;
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

group_t wilson(unsigned i, unsigned d1, unsigned d2, unsigned nd1, unsigned nd2)
{
	group_t g = id;
	
	for( int j=0; j<nd1; ++j)
	{
		g = mult[g][a[i*D+d1]];
		i=step(i,d1,1);
	}
	for( int j=0; j<nd2; ++j)
	{	
		g = mult[g][a[i*D+d2]];
		i=step(i,d2,1);
	}
	for( int j=0; j<nd1; ++j)
	{
		i=step(i,d1,-1);
        	g = mult[g][inv[a[i*D+d1]]];
	}
	for( int j=0; j<nd2; ++j)
	{
        	i=step(i,d2,-1);
		g = mult[g][inv[a[i*D+d2]]];
	}
        return g;
}
