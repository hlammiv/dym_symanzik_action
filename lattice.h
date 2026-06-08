#ifndef LATTICE_H
#define LATTICE_H

#include "group.h"

//extern unsigned D, L, V;
extern unsigned D, Nt, Nx, V;
extern double beta0, beta1, beta2, beta3, beta4;
extern group_t *a;

unsigned step(unsigned i, unsigned d, int s);
unsigned bigstep(unsigned i, unsigned d, int s);
group_t plaquette(unsigned i, unsigned d1, unsigned d2);
group_t polyakov(unsigned i);
void getpoly(double* re, double *im);
group_t wilson(unsigned i, unsigned d1, unsigned d2, unsigned nd1, unsigned nd2);

#endif /* ndef LATTICE_H */
