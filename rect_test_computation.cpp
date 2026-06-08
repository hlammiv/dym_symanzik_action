/*
 * Discrete Yang-Mills
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
#include "lattice.h"

#define K 1000      // Decorrelation time
#define N 10000   // Number of samples

unsigned D, Nt, Nx, V;
double beta0, beta1, beta2, beta3, beta4;
group_t *a;
unsigned int *nn = NULL; //creates a pointer with value 0

unsigned P;
double ReTr[PMAX];
double ImTr[PMAX];
group_t mult[PMAX][PMAX];
group_t inv[PMAX];
group_t id;

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

typedef unsigned int uint;

void update();
// void get_rect_loops(unsigned int i, unsigned int d, group_t* st);
// void get_staples(unsigned int i, unsigned int d, group_t* st);
double S_inv(unsigned i, unsigned d);

std::default_random_engine *rnd;

std::vector<group_t> smallgroup;

unsigned long hit = 0, acc = 0;


int main(int argc, char *argv[]) {
    int iseed = 0; //time(NULL) + getpid();

	/* Theory parameters come from the command line. */
	if (argc < 8) {
		fprintf(stderr, "usage: %s group D Nt Nx beta0 beta1 beta2\n", argv[0]);
		return 1;
	}
	const char *groupfilename = argv[1];
	D = atoi(argv[2]);
	Nt = atoi(argv[3]);
	Nx = atoi(argv[4]);
	beta0 = atof(argv[5]);
	beta1 = atof(argv[6]);
	beta2 = atof(argv[7]);
	if(argc == 9) iseed = atoi(argv[8]);
printf("PARAMS(grp,D,Nt,Nx,beta0,beta1,beta2,seed): %s %d %d %d %e %e %e %d\n", groupfilename, D, Nt, Nx, beta0, beta1, beta2, iseed);

	rnd = new std::default_random_engine[V];
	for(int i=0; i<V; ++i) rnd[i].seed(iseed+i);
	/* Load the multiplication table. */
	load_group(groupfilename);
	std::uniform_int_distribution<> randgrp(0,P-1);
printf("read the group\n"); fflush(stdout);
printf("id: %i\n",id);
printf("check stuff: (a*b)*c: %d a*(b*c): %d\n", mult[mult[1][2]][3], mult[1][mult[2][3]]);

	// select the group elems of S1080 that are close to identity
	double min_retr=0;
	double nn_retr=10;
	for(uint i=0; i<P; ++i) if(ReTr[i] < min_retr) min_retr=ReTr[i];
	for(uint i=0; i<P; ++i) if(ReTr[i] > min_retr && ReTr[i] < nn_retr) nn_retr=ReTr[i];
	for(uint i=0; i<P; ++i) if(ReTr[i] < nn_retr+1e-6 && ReTr[i] > min_retr) smallgroup.push_back(i);
	printf("min_retr: %e nn_retr: %e \n",min_retr,nn_retr);
	printf("small group size: %lu\n", smallgroup.size());

	/* Initialize the gauge field. */
	a = (group_t*) malloc(sizeof(unsigned) * V * D);
	for (unsigned i = 0; i < V*D; i++){
//			a[i] = randgrp(rnd[0]);
//			a[i] = rand()%P;
			a[i] = id;
//			if(i%2==0){a[i] = rand()%P;}else{a[i] = id;}
	}
step(0,0,1); //prime the nn table
printf("init done\n"); fflush(stdout);

//	for (unsigned i = 0; i < V*D; i++)
//                        printf(" %03d", a[i]);
//                printf("\n");


	for (unsigned k = 0; k < K*0; k++) update();
printf("thermo done\n"); fflush(stdout);

  for (unsigned n = 0; n < N; n++) 
  {
    timer tm;
    tm.start("update");
    for (unsigned k = 0; k < K; k++) update();
    tm.stop();
#if 0
    for (unsigned i = 0; i < V*D; i++)
      printf(" %03d", a[i]);
    printf("\n");
    fflush(stdout);
#endif
    // update done -- print plaq
    tm.start("meas");
//    double plaq = 0;
//#pragma omp parallel for reduction(+: plaq)
//    for(unsigned int i=0; i<V; ++i) for(unsigned int d=0; d<D; ++d) plaq += S_inv(i,d);
    double rep, imp; getpoly(&rep, &imp);
    // plaq computed is actually 4*action
    double simpleplaq = 0;
    double wloop[Nx][Nt];

    for(unsigned int i=0; i<Nx; ++i)
    for(unsigned int j=0; j<Nt; ++j)
    {
	    wloop[i][j]=0;
    }

#pragma omp parallel for reduction(+: simpleplaq)
    for(unsigned int i=0; i<V; ++i)
    for(unsigned int d1=0; d1<D; ++d1)
    for(unsigned int d2=d1+1; d2<D; ++d2)
    {
	    //printf("%e \n",ReTr[wilson(i,d1,d2,3,2)]);
	simpleplaq += ReTr[plaquette(i, d1, d2)];
    }
   
    for(unsigned int k=0; k<Nx; ++k)
    for(unsigned int l=0; l<Nt; ++l)
    {
	double lwloop=0;
	#pragma omp parallel for reduction(+: lwloop)
	for(unsigned int i=0; i<V; ++i)
    	for(unsigned int d1=1; d1<D; ++d1)
    	{
    		lwloop+=ReTr[wilson(i,d1,0,k+1,l+1)];
    	}
    	wloop[k][l]=lwloop/(V*D*(D-1)/2);
    }
   
    simpleplaq /= V*D*(D-1)/2;
    double cor;
    printf("GMES: %e %e %e %e", 999.0, rep, imp, simpleplaq);
//	printf("\n");
    for(unsigned int k=0; k<Nx; ++k)
    {
//	printf("WL %d: ",k);
    	for(unsigned int l=0; l<Nt; ++l) printf(" %e",wloop[k][l]);
//	printf("\n");
    }

    printf("\n");
    fflush(stdout);
    printf("ACC: %f\n", ((double)acc)/hit);

    
    printf("CONFIGS: ");
      for (unsigned i = 0; i < V*D; i++)
                        printf(" %03d", a[i]);
                printf("\n");

    tm.stop();
  }

	free(a);
	return 0;
}
//NEED TO ADD NEW ACTION TERMS
double S_inv(unsigned i, unsigned d) {
	double r1 = 0., r2 = 0.;
	for (unsigned dp = 0; dp < D; dp++) {
		if (dp == d) continue;
		group_t p1 = plaquette(i, d, dp);
		group_t p2 = plaquette(step(i, dp, -1), d, dp);
//printf("(i,d,dp): %d %d %d\n", i, d, dp);
//printf("p1: %d p2: %d trp1: %e trp2: %e\n", p1, p2, ReTr[p1], ReTr[p2]);
		r1 += ReTr[p1] + ReTr[p2];
		r2 += ReTr[mult[p1][p1]] + ReTr[mult[p2][p2]];
	}
	return beta0 + beta1 * r1 + beta2 * r2;
}

void get_staples(unsigned int i, unsigned int d, group_t* st)
{
  unsigned int i1 = step(i, d, 1); int k=0;
  for(unsigned int d1=0; d1<D; ++d1) if(d1 != d)
  {
//printf("(i,d,d1): %d %d %d\n", i, d, d1);
//printf("[123]: %d %d %d\n", a[i1*D+d1], inv[a[step(i,d1,1)*D+d]], inv[a[i*D+d1]]);
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
    //printf("(i,d1,d2): %d %d %d\n", i, d1, d2);
    //printf("[0123]: %d %d %d %d\n", a[i*D+d1], a[step(i,d1,1)*D+d2], inv[a[step(i,d2,1)*D+d1]], inv[a[i*D+d2]]);
	  // g = mult[g][a[i*D+d]];
	  g = mult[g][a[step(i,d,1)*D+d1]];
    int i2 = step(step(i,d,1),d1,1);
    g = mult[g][a[i2*D+d1]];
    g = mult[g][inv[a[step(step(i2,d1,1),d,-1)*D+d]]];
    g = mult[g][inv[a[step(i,d1,1)*D+d1]]];
	  g = mult[g][inv[a[i*D+d1]]];
    st[k++] = g;
    // printf("a product: %d\n", g);

    g = id;
    // g = mult[g][a[i*D+d]]
    g = mult[g][a[step(i,d,1)*D+d]];
    i2 = step(step(i,d,1),d,1);
    g = mult[g][a[i2*D+d1]];
    g = mult[g][inv[a[step(step(i2,d1,1),d,-1)*D+d]]];
    g = mult[g][inv[a[step(i,d1,1)*D+d]]];
    g = mult[g][inv[a[i*D+d1]]];
    st[k++] = g;
    // printf("a product: %d\n", g);

    g = id;
    // g = mult[g][a[i*D+d]]
    g = mult[g][a[step(i,d,1)*D+d1]];
    i2 = step(step(i,d,1),d1,1);
    g = mult[g][inv[a[step(i2,d,-1)*D+d]]];
    g = mult[g][inv[a[step(step(i2,d,-1),d,-1)*D+d]]];
    g = mult[g][inv[a[step(i,d,-1)*D+d1]]];
    g = mult[g][a[step(i,d,-1)*D+d]];
    st[k++] = g;
    // printf("a product: %d\n", g);

    g = id;
    // g = mult[g][a[i*D+d]]
    g = mult[g][inv[a[step(step(i,d,1),d1,-1)*D+d1]]];
    i2 = step(step(i,d,1),d1,-1);
    g = mult[g][inv[a[step(i2,d1,-1)*D+d1]]];
    g = mult[g][inv[a[step(step(i2,d1,-1),d,-1)*D+d]]];
    g = mult[g][a[step(step(i,d1,-1),d1,-1)*D+d1]];
    g = mult[g][a[step(i,d1,-1)*D+d1]];
    st[k++] = g;
    // printf("a product: %d\n", g);

    g = id;
    // g = mult[g][a[i*D+d]]
    g = mult[g][a[step(i,d,1)*D+d]];
    i2 = step(step(i,d,1),d,1);
    g = mult[g][inv[a[step(i2,d1,-1)*D+d1]]];
    g = mult[g][inv[a[step(step(i2,d1,-1),d,-1)*D+d]]];
    g = mult[g][inv[a[step(i,d1,-1)*D+d]]];
    g = mult[g][a[step(i,d1,-1)*D+d1]];
    st[k++] = g;
    // printf("a product: %d\n", g);

    g = id;
    // g = mult[g][a[i*D+d]]
    g = mult[g][inv[a[step(step(i,d,1),d1,-1)*D+d1]]];
    i2 = step(step(i,d,1),d1,-1);
    g = mult[g][inv[a[step(i2,d,-1)*D+d]]];
    g = mult[g][inv[a[step(step(i2,d,-1),d,-1)*D+d]]];
    g = mult[g][a[step(step(i,d,-1),d1,-1)*D+d1]];
    g = mult[g][a[step(i,d,-1)*D+d]];
    st[k++] = g;
  }
}

extern "C" void getpos(unsigned idx, unsigned int *pos);
void update() 
{
  timer tm(false); tm.start("update");
  std::uniform_int_distribution<> randsmallgrp(0,smallgroup.size()-1);
  std::uniform_real_distribution<> rand01(0.,1.);
  double kept = 0;
  double total = 0;
  for (unsigned d = 0; d < D; d++) 
  for (int parity = 0; parity < 2; ++parity)
  {
    int lhit=0, lacc=0;
#pragma omp parallel for schedule(static), reduction(+: lacc, lhit)
  for (uint i = 0; i < V; i++) 
  {
    uint pos[4+00];
    getpos(i, pos); if( (pos[0]+pos[1]+pos[2]+pos[3])%2 != parity) continue;
    group_t staples[2*(D-1)+00]; 
    get_staples(i, d, staples);
    group_t rectangles[6 * (D - 1) + 00];
    get_rect_loops(i, d, rectangles);
    int nhit = 20;
    for(int h=0; h<nhit; ++h)
    {
      double r1=0, r2=0, r1new=0, r2new=0, m1 = 0, m1new = 0;
      group_t b = a[i*D+d];
//      group_t bnew = mult[b][smallgroup[randsmallgrp(rnd[i])]];
      group_t bnew = mult[b][smallgroup[smallgroup.size()*rand01(rnd[i])]];
      for(int j=0; j<2*(D-1); ++j) 
      {
        group_t p1 = mult[b][staples[j]];
//      group_t p2 = mult[p1][p1];
        r1 += ReTr[p1];

        p1 = mult[bnew][staples[j]];
//      p2 = mult[p1][p1];
        r1new += ReTr[p1];
      }
      for(int k = 0; k < 6 * (D-1); ++k){
        group_t rect = mult[b][rectangles[k]];
        m1 += ReTr[rect];

        rect = mult[bnew][rectangles[k]];
        m1new += ReTr[rect];
      }
      double oldact = beta0 + beta1*r1 + beta2*m1; 
      double newact = beta0 + beta1*r1new + beta2*m1new;
      double changeact = newact - oldact;
      double probrat = exp(-(newact-oldact));
      // printf("Action change: %e \n", changeact);
//#pragma omp atomic
//      hit++;
      lhit++;
      total += 1;
      if(probrat > rand01(rnd[i])) {
        a[i*D+d] = bnew;
        // printf("true"); 
//#pragma omp atomic
//	acc++; 
	      lacc++; 
//	if(hit%(100*V*D)==0) printf("ACC: %f\n", ((double)acc)/hit);
      }
    }
  }

  hit += lhit;
  acc += lacc;
  }
  double acceptrate = 40/total;
  printf("acceptance rate: %e", acceptrate);
  tm.stop();
  
}

unsigned int getidx(unsigned int *pos)
{
  int idx = 0;
  for(int d=D-1; d>0; --d) idx = (pos[d]+Nx)%Nx + Nx*idx; 
//printf("pos: %d %d %d %d -> %d (Nt:%d idx:%d Nx:%d)\n", pos[0], pos[1], pos[2], pos[3], (pos[0] + Nt)%Nt + Nt*idx, Nt, idx, Nx);
  return (pos[0] + Nt)%Nt + Nt*idx;
}

void init_nn()
{
  printf("Initializing nearest neighbours [%d:%d]... \n", Nt, Nx);
  // nn = malloc(sizeof(unsigned int) * V * D * 2);
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

void load_group(const char *fn) {
	FILE *fin = fopen(fn, "r");
	/* File format for group specification is as follows:
	 *
	 *  <order>
	 *  ReTr0 ReTr1 ReTr2 ...
	 *  ImTr0 ImTr1 ImTr2 ...
	 *  0x0 0x1 0x2 ...
	 *  1x0 1x1 1x2 ...
	 *  ...
	 *
	 * For a group of order P, there are thus 1 + P + P^2 entries. Whitespace is ignored.
	 */
	fscanf(fin, "%d", &P);
	if (P > PMAX) {
		fprintf(stderr, "Order of group too large: %d > %d\n", P, PMAX);
		abort();
	}
	for (unsigned n = 0; n < P; n++)
		fscanf(fin, "%lf", &ReTr[n]);
        for (unsigned n = 0; n < P; n++)
                fscanf(fin, "%lf", &ImTr[n]);
	for (unsigned n = 0; n < P; n++)
		for (unsigned m = 0; m < P; m++)
			fscanf(fin, "%d", &mult[n][m]);
	fclose(fin);

	// Find the identity.
	char id_found = 0;
	for (unsigned n = 0; n < P; n++) {
		char is_id = 1;
		for (unsigned m = 0; m < P; m++)
			if (mult[n][m] != m || mult[m][n] != m)
				is_id = 0;
		if (is_id) {
			id_found = 1;
			id = n;
			break;
		}
	}
	if (!id_found) {
		fprintf(stderr, "Group does not have identity element\n");
		abort();
	}

	// Find inverses.
	for (unsigned n = 0; n < P; n++) {
		char inv_found = 0;
		for (unsigned m = 0; m < P; m++)
			if (mult[n][m] == id && mult[m][n] == id) {
				inv[n] = m;
				inv_found = 1;
				break;
			}
		if (!inv_found) {
			fprintf(stderr, "Group does not have inverses\n");
			abort();
		}
	}
}