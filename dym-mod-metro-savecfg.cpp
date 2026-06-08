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
#include <fstream>

#include "timer.h"

#define K 200      // Decorrelation time
#define N 1000   // Number of samples

typedef unsigned int uint;

void update();
double S_inv(unsigned i, unsigned d);

std::default_random_engine *rnd;

std::vector<group_t> smallgroup;

struct su3 { double data[18]; double& operator()(int r, int c, int part) { return data[part+2*(c+3*r)]; } };
std::vector<su3> groupsu3;

// switch endianess in the array
void switchend(unsigned char *buffer, int length)
{
  int i, j;
  unsigned char *pos;
  unsigned char save[8];

  pos = buffer;

  for(j=0; j<length; j++)
  {
    for(i=0; i<8; i++) save[i] = pos[i];
    for(i=7; i>-1; i--, pos++) *pos=save[i];
  }
}

extern "C" void getpos(unsigned idx, unsigned int *pos);
void savekyconfig(char* name)
{
  double *buf = new double[V*4*sizeof(su3)];
  uint pos[4];
  for(uint i=0; i<V; ++i)
  {
    getpos(i, pos);
    for(int d=0; d<D; ++d)
    {
      su3 &m = groupsu3[a[i*D+d]];
      for(int part = 0; part<2; ++part) for(int c=0; c<3; ++c) for(int r=0; r<3; ++r)
        buf[pos[1]+Nx*(pos[2]+Nx*(pos[3]+Nx*(pos[0]+Nt*(part+2*(r+3*(c+3*((d+3)%D)))))))] = m(r, c, part);
    }
  }

  switchend((unsigned char*)buf, V*4*18);

  FILE* f = fopen(name, "w");
  fwrite(buf, V*4, sizeof(su3), f);
  fclose(f);
  delete [] buf;
}

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

{
  std::ifstream ingroup(groupfilename);
  int num;
  ingroup >> num;
  double tmp; int itmp;
  for(int i=0; i<2*num; ++i) ingroup >> tmp;
  for(int i=0; i<num*num; ++i) ingroup >> itmp;
  su3 mat;
  groupsu3.resize(num);
  for(int i=0; i<num; ++i)
  {
    double * pmat = groupsu3[i].data;
    for(int j=0; j<18; ++j, ++pmat) ingroup >> *pmat;
  }
  ingroup.close();
}

	std::uniform_int_distribution<> randgrp(0,P-1);
printf("read the group\n"); fflush(stdout);
printf("id: %i\n",id);
printf("check stuff: (a*b)*c: %d a*(b*c): %d\n", mult[mult[1][2]][3], mult[1][mult[2][3]]);

	// select the group elems of S1080 that are close to identity
        for(uint i=0; i<P; ++i) if(fabs(ReTr[i]) > 1.6 && fabs(ReTr[i]) < 3) smallgroup.push_back(i);
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
#pragma omp parallel for reduction(+: simpleplaq)
    for(unsigned int i=0; i<V; ++i)
    for(unsigned int d1=0; d1<D; ++d1)
    for(unsigned int d2=d1+1; d2<D; ++d2)
      simpleplaq += ReTr[plaquette(i, d1, d2)];
    simpleplaq /= V*D*(D-1)/2;
    printf("GMES: %e %e %e %e\n", 999.0, rep, imp, simpleplaq); fflush(stdout);
    printf("ACC: %f\n", ((double)acc)/hit);
    char cfgname[100];
    sprintf(cfgname, "/lustre/hlamm/wilson_flow/lat-beta1%.3f-beta2%3f-nt%02d-nx%02d-num%04d", beta1, beta2, Nt, Nx, n);
    savekyconfig(cfgname);
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

#if 0
void update_hb() 
{
  for (unsigned d = 0; d < D; d++) 
  for (unsigned i = 0; i < V; i++) 
  {
    group_t staples[2*(D-1)]; 
    get_staples(i, d, staples);
    double pr[P];
    for (unsigned b = 0; b < P; b++) 
    {
      double r1=0, r2=0;
      for(int j=0; j<2*(D-1); ++j) 
      {
	group_t p1 = mult[b][staples[j]];
//printf("p1: %d trp1: %e\n", p1, ReTr[p1]);
	group_t p2 = mult[p1][p1];
	r1 += ReTr[p1];
	r2 += ReTr[p2];
      }
      double act = beta0 + beta1*r1 + beta2*r2; 
//a[i*D+d] = b;
//printf("act: %e old: %e\n", act, S_inv(i,d));
      pr[b] = exp(-act);
    }


    // Pick a random number (normalized to Pc).
    std::discrete_distribution<> heatbath(pr,pr+P);
    a[i*D+d] = heatbath(rnd);
  }

}
#endif
