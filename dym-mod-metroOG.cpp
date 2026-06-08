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

#define K 100      // Decorrelation time
#define N 100   // Number of samples

typedef unsigned int uint;

void update();
double S_inv(unsigned i, unsigned d);

std::default_random_engine *rnd;

std::vector<group_t> smallgroup;

unsigned long hit = 0, acc = 0;

int main(int argc, char *argv[]) {
	/* Seed the PRNG. */
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

	V = Nt*pow(Nx,D-1);
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
			a[i] = randgrp(rnd[0]);
//			a[i] = rand()%P;
//			a[i] = id;
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

extern "C" void getpos(unsigned idx, unsigned int *pos);
void update() 
{
  timer tm(false); tm.start("update");
  std::uniform_int_distribution<> randsmallgrp(0,smallgroup.size()-1);
  std::uniform_real_distribution<> rand01(0.,1.);
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
    int nhit = 20;
    for(int h=0; h<nhit; ++h)
    {
      double r1=0, r2=0, r1new=0, r2new=0;
      group_t b = a[i*D+d];
//      group_t bnew = mult[b][smallgroup[randsmallgrp(rnd[i])]];
      group_t bnew = mult[b][smallgroup[smallgroup.size()*rand01(rnd[i])]];
      for(int j=0; j<2*(D-1); ++j) 
      {
        group_t p1 = mult[b][staples[j]];
//      group_t p2 = mult[p1][p1];
        r1 += ReTr[p1];
//      r2 += ReTr[p2];
        r2 += ReTr[p1]*ReTr[p1];
        r2 += ImTr[p1]*ImTr[p1];

        p1 = mult[bnew][staples[j]];
//      p2 = mult[p1][p1];
        r1new += ReTr[p1];
//      r2new += ReTr[p2];
        r2new += ReTr[p1]*ReTr[p1];
        r2new += ImTr[p1]*ImTr[p1];

      }
      double oldact = beta0 + beta1*r1 + beta2*r2; 
      double newact = beta0 + beta1*r1new + beta2*r2new; 
      double probrat = exp(-(newact-oldact));
//#pragma omp atomic
//      hit++;
      lhit++;
      if(probrat > rand01(rnd[i])) { 
	a[i*D+d] = bnew; 
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

  tm.stop();
}

