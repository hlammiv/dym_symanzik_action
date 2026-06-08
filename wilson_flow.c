/*
 * wilson_flow.c             by Szabolcs Borsanyi, 2012
 * ----------------------------------------------------
 * Wilson/Symanzik flow calculator for lattice gauge configuartions.
 *
 * The program calculates the wilson flow by implementing
 * Luscher's integrator in arXiv:1006.4518, appendix C.
 * In the present version the field strength and the topological charge
 * is written as a function of the  flow time 't'.
 *
 * Use the w0_scale.c program to extract the lattice spacing.
 *
 * The gauge configurations are given on the command line, one by one,
 * there will be an output file for each of them.
 *
 * This is a standalone c99 program, a reference implementation.
 * You are encouraged to read it, not only run it.
 *
 * With the gcc compiler you may need to include the -c99 flag
 *
 *  gcc -std=c99 wilson_flow.c  -lm
 *
 * or
 *
 *  gcc -std=c99 -fopenmp -DHOST_BIG_ENDIAN   -DENABLE_OMP  wilson_flow.c  -lm
 *  ./a.out
 * will print out the help pages.
 *
 * openMP parallelism is supported, check the ENABLE_OMP macro.
 * There is no MPI parallelism here, but you could still use this
 * implementation in 'production' on e.g. a cluster of multi-core 
 * workstations in farming mode, where several single-node
 * jobs analyze their share of a configuration set.
 * But be warned that this program has not been optimized for speed.
 *
 * For example:
 * gcc -std=c99 -fopenmp  -DENABLE_OMP  wilson_flow.c  -lm

 * Byte orders: If the program compiles and links without any setting,
 * than you seem to have the BSD extensions. If not, define the
 * HOST_LITTLE_ENDIAN or HOST_BIG_ENDIAN macro to 1 and then the program
 * will not need those extensions.
 *
 * For example: 
 * gcc -std=c99 -fopenmp -DHOST_BIG_ENDIAN   wilson_flow.c  -lm
 *
 * gcc -std=c99 -fopenmp -DHOST_LITTLE_ENDIAN   wilson_flow.c  -lm
 *
 * Lattice format:
 *  For now, only the NERSC format is understood.
 *
 *  look for help()
 *
 * Last modified: 7 March 2010
 */

/* Headers {{{ */
#define _BSD_SOURCE

/* the standard C99 headers */
#include <complex.h>
#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <float.h>
/* POSIX, should be available mostly everywhere */
#include <getopt.h>

/* HOST Byte order */

/* If you know, the byte order of your workstation, or compute node,
 * you can set it here. */

// #define HOST_BIG_ENDIAN 1
// #define HOST_LITTLE_ENDIAN 1

/* Or perhaps your precompiler knows */

#if defined(BYTE_ORDER) && defined(BIG_ENDIAN) && defined(LITTLE_ENDIAN) && !defined(HOST_LITTLE_ENDIAN) && !defined(HOST_BIG_ENDIAN)
#if BYTE_ORDER == BIG_ENDIAN
#define HOST_BIG_ENDIAN 1
#elif BYTE_ORDER == LITTLE_ENDIAN
#define HOST_LITTLE_ENDIAN 1
#endif
#endif


/* openMP */

/* switch openMP support on and off */
// #define ENABLE_OMP 1

#ifdef ENABLE_OMP
#include <omp.h>
#endif

/* Headers }}}*/

int help()/*{{{*/
{
  puts("wilson_flow  (by the Budapest-Marseille-Wuppertal collaboration)");
  puts("----------------------------------------------------------------");
  puts("Usage:");
  puts("       wilson_flow [options] conf1 conf2 ... ");
  puts("");
  puts("This program integrates the wilson flow into flow.conf1 flow.cont2 ... files.");
  puts("Options:");
  puts("        -e<eps>         integration flow time step");
  puts("        -t<t-max>       finish when flow time reaches this");
  puts("        -x<xi>          (renormalized) anisotropy");
  puts("        -s              use the tree level Symanzik discretization");
  puts("        -f<format>      format of the input configurations");
  puts("");
  puts("Input format:");
  puts(" -f nersc     see https://qcdlattices.bnl.gov/formatDescription/index.html");
  puts("");
  puts("Output format:");
  puts(" one line for each flow time:");
  puts("WFLOW <t> <Et> <Es> <topological_charge>");
  puts("");
  puts("Questions? Ask borsanyi@uni-wuppertal.de");
  return EXIT_FAILURE;
}/*}}}*/

/* program parameters {{{ */
double Xi=1.0;          /* Anisotropy in the flow equation: Xi=a_s/a_t */
/* program parameters }}} */

/* data structures {{{ */

/* Linear algebra structures */
/* The links are stored in C99's native double precision complex format.
 * We only store the upper two rows of the links using the 'su3'
 * type.
 * Other matrices have the 'mat' type. Antihermitian traceless matrices
 * come about very frequently, they are stored in 'ahm' typed objects.
 */

typedef struct mat {
  complex double c[3][3];  /* 3 x 3 format */
} mat;

typedef struct su3 {
  complex double c[2][3];  /* two upper rows */
} su3;

typedef struct ahm {
  /* Antihermitian traceless matrix: only some components are stored */
  /* See M_A() for the conversion. Such matrices are usually a result
   * of a projection, which we implement in A_M() */
  complex double c01;
  complex double c02;
  complex double c12;
  double i00;
  double i11;
} ahm;

/* The following pointers are used to access lattice fields.
 * The sites are layed out as a big one-dimensional array.
 * Only one field component (e.g. link direction) is stored in one array.
 * The whole configuration is stored in a "latsu3 U[4]".
 * Use the neighbour tables neigh_b, neigh_f to find the lattice neighbours.
 *
 * These pointers are initialized with the
 * kernel_U_new(), kernel_M_new(), kernel_A_new() calls, 
 * and deleted by
 * kernel_U_delete(), kernel_M_delete(), kernel_A_delete() calls, 
 *
 * The kernel_U_new4, kernel_U_delete4 and similar functions create and
 * delete four such fields at a time.
 */
typedef struct su3 *latsu3;  /* link   field, typically referred to as "U" */
typedef struct mat *latmat;  /* matrix field, typically referred to as "M" */
typedef struct ahm *latahm;  /* algebra field, ypically referred to as "A" */


/* data structures }}} */

/* matrix algebra {{{ */

static inline void U_reunitarize(su3 *U)/*{{{*/
{
  /* The two upper rows of an su3 matrix is corrected such that
   * the resulting two rows correspond to an su3 matrix.
   * This function is often used after restoring configurations from
   * single precision storage. */
#define rdot(A,B) (creal((A)[0])*creal((B)[0])+cimag((A)[0])*cimag((B)[0])+ creal((A)[1])*creal((B)[1])+cimag((A)[1])*cimag((B)[1])+ creal((A)[2])*creal((B)[2])+cimag((A)[2])*cimag((B)[2]))
#define cdot(A,B) (A)[0]*conj((B)[0])+(A)[1]*conj((B)[1])+(A)[2]*conj((B)[2])
  double n;
  complex double c;
  n=rdot(U->c[0],U->c[0]);
  n=1./sqrt(n);
  U->c[0][0]*=n; 
  U->c[0][1]*=n;
  U->c[0][2]*=n;
  c=cdot(U->c[1],U->c[0]);
  U->c[1][0]-=c*U->c[0][0];
  U->c[1][1]-=c*U->c[0][1];
  U->c[1][2]-=c*U->c[0][2];
  n=rdot(U->c[1],U->c[1]);
  n=1.0/sqrt(n);
  U->c[1][0]*=n;
  U->c[1][1]*=n;
  U->c[1][2]*=n;
#undef rdot
#undef cdot
}/*}}}*/

static inline void M_U(mat *M,const su3 *U)/*{{{*/
{
  /* reconstruction of an su3 matrix */
  M->c[0][0]=U->c[0][0];
  M->c[0][1]=U->c[0][1];
  M->c[0][2]=U->c[0][2];
  M->c[1][0]=U->c[1][0];
  M->c[1][1]=U->c[1][1];
  M->c[1][2]=U->c[1][2];
  M->c[2][0]=conj(M->c[0][1]*M->c[1][2]-M->c[0][2]*M->c[1][1]);
  M->c[2][1]=conj(M->c[0][2]*M->c[1][0]-M->c[0][0]*M->c[1][2]);
  M->c[2][2]=conj(M->c[0][0]*M->c[1][1]-M->c[1][0]*M->c[0][1]);
}/*}}}*/

static inline void U_M(su3 *U,const mat *M)/*{{{*/
{       
  U->c[0][0]=M->c[0][0];
  U->c[0][1]=M->c[0][1];
  U->c[0][2]=M->c[0][2];
  U->c[1][0]=M->c[1][0];
  U->c[1][1]=M->c[1][1];
  U->c[1][2]=M->c[1][2];
  U_reunitarize(U);
}       
/*}}}*/

static inline void M_M(mat *D,const mat *S)/*{{{*/
{
  /* plain matrix copy */
  memcpy(D,S,sizeof(mat));
}/*}}}*/

static inline void M_ct1(mat *M,complex double c)/*{{{*/
{
  /* initialization of a complex unit matrix */
  M->c[0][0]=c;
  M->c[0][1]=0;
  M->c[0][2]=0;
  M->c[1][0]=0;
  M->c[1][1]=c;
  M->c[1][2]=0;
  M->c[2][0]=0;
  M->c[2][1]=0;
  M->c[2][2]=c;
}/*}}}*/

static inline complex double c_ahm(complex double a,complex double b)/*{{{*/
{
  return 0.5*(a-conj(b));
}/*}}}*/

static inline double cnorm(complex double a)/*{{{*/
{
  register double r=creal(a);
  register double i=cimag(a);
  return r*r+i*i;
}/*}}}*/

static inline void A_M(ahm *A,const mat *M)/*{{{*/
{ /* projection: A = ( M - adj(M) ) / 2  - 1*Im(Tr(M))/3   */
  double tr;
  A->c01=c_ahm(M->c[0][1],M->c[1][0]);
  A->c02=c_ahm(M->c[0][2],M->c[2][0]);
  A->c12=c_ahm(M->c[1][2],M->c[2][1]);
  A->i00=cimag(M->c[0][0]);
  A->i11=cimag(M->c[1][1]);
  tr=(A->i00+A->i11+cimag(M->c[2][2]))/3.;
  A->i00-=tr;
  A->i11-=tr;
}/*}}}*/

static inline void M_A(mat *M,const ahm *A)/*{{{*/
{
  /* An antihermitian traceless matrix 'A' is converted
   * to the full matrix format, where the [2][2] element is reconstructed
   * so that the resulting matrix becomes traceless */
  double i22=-(A->i00+A->i11);
  M->c[0][0]=I*A->i00;
  M->c[1][1]=I*A->i11;
  M->c[2][2]=I*i22;
  M->c[0][1]=A->c01;
  M->c[1][0]=-conj(A->c01);
  M->c[0][2]=A->c02;
  M->c[2][0]=-conj(A->c02);
  M->c[1][2]=A->c12;
  M->c[2][1]=-conj(A->c12);
}/*}}}*/

static inline void A_StA(ahm *R,double s,const ahm *A)/*{{{*/
{
  const double *a=(const void*)A;
  double *r=(void*)R;
  r[0]=s*a[0];
  r[1]=s*a[1];
  r[2]=s*a[2];
  r[3]=s*a[3];
  r[4]=s*a[4];
  r[5]=s*a[5];
  r[6]=s*a[6];
  r[7]=s*a[7];
}/*}}}*/

static inline void Ap_StA(ahm *R,double s,const ahm *A)/*{{{*/
{
  const double *a=(const void*)A;
  double *r=(void*)R;
  r[0]+=s*a[0];
  r[1]+=s*a[1];
  r[2]+=s*a[2];
  r[3]+=s*a[3];
  r[4]+=s*a[4];
  r[5]+=s*a[5];
  r[6]+=s*a[6];
  r[7]+=s*a[7];
}/*}}}*/

/* Matrix exponentialization: void M_expA(mat *M,const ahm *A) {{{ */

/* The following set of somewhat sophisticated functions serve one
 * single purpose: to calculate the exponential function of
 * an antihermitian traceless matrix.
 * The user should access this through the  M_expA() call.
 *
 * The formulas correspond to Morningstar&Peardon hep-lat/0311018
 */

static inline void A_add_poly(mat *M,const ahm *A,complex double a,complex double b)/*{{{*/
  /* mat+=a*ahm + b* ahm^2 */
{   
  double i22=-(A->i00+A->i11);
  /* mat+=a*ahm */
  M->c[0][0]+=I*a*A->i00;
  M->c[0][1]+=a*A->c01;
  M->c[0][2]+=a*A->c02;
  M->c[1][0]-=a*conj(A->c01);
  M->c[1][1]+=I*a*A->i11;
  M->c[1][2]+=a*A->c12;
  M->c[2][0]-=a*conj(A->c02);
  M->c[2][1]-=a*conj(A->c12);
  M->c[2][2]+=I*a*i22;
  /* mat+=b*ahm */
  /* diagonal first */
  double s01,s02,s12;
  s01=cnorm(A->c01);
  s02=cnorm(A->c02);
  s12=cnorm(A->c12);
  M->c[0][0]-=(A->i00*A->i00+s01+s02)*b;
  M->c[1][1]-=(A->i11*A->i11+s01+s12)*b;
  M->c[2][2]-=(i22*i22+s02+s12)*b;
  complex double cc;
  cc=A->c02*conj(A->c12)+I*i22*A->c01;
  M->c[0][1]-=cc*b;
  M->c[1][0]-=conj(cc)*b;
  cc=A->c01*A->c12-I*A->i11*A->c02;
  M->c[0][2]+=cc*b;
  M->c[2][0]+=conj(cc)*b;
  cc=A->c02*conj(A->c01)+I*A->i00*A->c12;
  M->c[1][2]-=cc*b;
  M->c[2][1]-=conj(cc)*b;
}/*}}}*/

static inline void A_trace23(double *traa, double *traaa,const ahm *A)/*{{{*/
  /* traa= retr(ahm^2), traaa= imtr(ahm^3) */
  /* some piece of maple code:

     A:=Matrix([[I*im00, re01+I*im01, re02+I*im02], [-re01+I*im01, I*im11, re12+I*im12], [-re02+
     I*im02, -re12+I*im12, -I*(im00+im11)]]);
     with(LinearAlgebra):
     with(linalg):
     collect(simplify(evalc( Im(trace(A.A.A)))), {im11, im00});
   */
{
  /* element squares */
  double v00= A->i00*A->i00;
  double v11= A->i11*A->i11;
  double v01= cimag(A->c01)*cimag(A->c01);
  double v02= cimag(A->c02)*cimag(A->c02);
  double v12= cimag(A->c12)*cimag(A->c12);

  double u01= creal(A->c01)*creal(A->c01);
  double u02= creal(A->c02)*creal(A->c02);
  double u12= creal(A->c12)*creal(A->c12);

  /* traa= re tr(A^2) = -2*(a00^2 + a11^2 + a00*a11 + |a01|^2 + |a02|^2 + |a12|^2) */
  /* it is just the matrix normsq */
  *traa=-2*( u01+v01+u02+v02+u12+v12+v00+v11+A->i00*A->i11);

  /* traaa= im tr(A^3) = 3term0 + 3term1 + 6term */

  /* 3term0= 3 im00 (im11^2 + |a12|^2 - |a01|^2) */
  *traaa=3.*A->i00*(v11+u12+v12-u01-v01);

  /* 3term1= 3 im11 (im00^2 + |a02|^2 - |a01|^2) */
  *traaa+=3.*A->i11*(v00+u02+v02-u01-v01);

  /* 6term = 6 (- re02 im12 re01 + im02 re12 re01 - re02 re12 im01 - im02 im12 im01 ) */
  *traaa+=6.*(
      -creal(A->c02)*cimag(A->c12)*creal(A->c01)
      +cimag(A->c02)*creal(A->c12)*creal(A->c01)
      -creal(A->c02)*creal(A->c12)*cimag(A->c01)
      -cimag(A->c02)*cimag(A->c12)*cimag(A->c01)
      );
}/*}}}*/

static inline void emr_e2r_r ( complex double *emu, complex double *e2u, double u )/*{{{*/
{
  /* exp(-iu), exp(i*2*u) */
  complex double eu=cexp(I*u);
  *emu=conj(eu);
  *e2u=eu*eu;
}
/*}}}*/

/* a-> cos(a), ksi0(a)= sin(a)/a, ksi1(a)= cos(a)/a^2 - sin(a)/a^3 */
/* for small 'a'-s taylor expansion used in ksi1 (see doc/expa.ps) */
static inline void cosr_ksir_r( double *cosa, double *ksi0a, double *ksi1a, double a )/*{{{*/
{
#ifndef KSI1_TAYLOR_DEFINES
#define KSI1_TAYLOR_DEFINES
  /* taylor[6]={
     -1.0/3.0,
     +1.0/30.0,
     -1.0/840.0,
     +1.0/45360.0,
     -1.0/3991680.0,
     +1.0/518918400.0}
   */
#define C_EXPA_TAYLOR0 -0.3333333333333333333333333333333333333333
#define C_EXPA_TAYLOR1 0.03333333333333333333333333333333333333333
#define C_EXPA_TAYLOR2 -0.001190476190476190476190476190476190476190
#define C_EXPA_TAYLOR3 0.00002204585537918871252204585537918871252205
#define C_EXPA_TAYLOR4 -0.2505210838544171877505210838544171877505e-6
#define C_EXPA_TAYLOR5 0.1927085260418593751927085260418593751927e-8
  /* use taylor below these numbers */
#define C_EXPA_TAYLOR_BOUND  0.50
#endif
  double a2=a*a;
  complex double ce=cexp(I*a);
  *cosa=creal(ce);
  if(a2>0.0025) *ksi0a=cimag(ce)/a;
  else *ksi0a=1-(1-(1-a2/42.)*a2/20.)*a2/6;
  if(a2>(C_EXPA_TAYLOR_BOUND)*(C_EXPA_TAYLOR_BOUND)) {
    *ksi1a=(creal(ce)-*ksi0a)/a2;
  } else {
    /* ksi1-taylor */
    double tay;
    tay= C_EXPA_TAYLOR5;
    tay= C_EXPA_TAYLOR4+a2*tay;
    tay= C_EXPA_TAYLOR3+a2*tay;
    tay= C_EXPA_TAYLOR2+a2*tay;
    tay= C_EXPA_TAYLOR1+a2*tay;
    tay= C_EXPA_TAYLOR0+a2*tay;
    *ksi1a=tay;
  }
  /* make it safe : if a<DBL_MIN */
  if(a<DBL_MIN) {
    *cosa=1.;
    *ksi0a=1.;
    *ksi1a=C_EXPA_TAYLOR0;
  }
}/*}}}*/

/* for exp: exp(A)=f0+f1*A+f2*A^2 with A ahm_matrix */
/* traa= re tr (A^2); traaa= im tr (A^3) */
static inline void exp_trig0( complex double f[3], double traa, double traaa )/*{{{*/
{
#define C_SQRT3  1.7320508075688772935274463415058723669430
  complex double c, emu, e2u;
  double signc0, c0, c1o3, c0max, fround;
  double theta;
  double cosw, ksi_0w, ksi_1w;
  double sinth, costh;
  double sqrtc1o3;
  double u, w, u2, w2, oodenom, rr;
  /* c0, c1/3, c2 - with roundoff handling */
  c0= traaa*(-1./3.);
  signc0= copysign(1.,c0);
  c0*=signc0;
  c1o3= traa*(-1./6.);
  sqrtc1o3= sqrt(c1o3);
  /* make it safe: if (c1o3<0) sqrtc1o3=0 */
  if(c1o3<DBL_MIN) sqrtc1o3=0.;
  c0max=2.*c1o3*sqrtc1o3;
  fround= c0/c0max;
  /* make it safe: fround cannot be larger than 1.0 */
  fround= fmin(fround,1.);
  /* u, w, theta */
  theta = acos ( fround );
  rr= theta*(1./3.);
  sinth=sin(rr);
  costh=cos(rr);
  u = sqrtc1o3*costh;
  w = sqrtc1o3 *C_SQRT3 *sinth;
  emr_e2r_r (&emu, &e2u, u);

  /* speedup */
  cosr_ksir_r (&cosw, &ksi_0w, &ksi_1w, w);
  u2= u*u;
  w2= w*w;
  oodenom=1./( 9.*u2-w2 );
  /* h[] */
  rr= u2-w2;
  f[0]= rr*e2u;
  rr= w2+3.*u2;
  c=8.*u2*cosw + I*2.*u*rr*ksi_0w;
  f[0]+=emu*c;
  rr=2.*u;
  f[1]= rr*e2u;
  rr=3.*u2-w2;
  c=-2.*u*cosw+I*rr*ksi_0w;
  f[1]+=emu*c;
  f[2]= e2u;
  c=-cosw-3.*I*u*ksi_0w;
  f[2]+=emu*c;
  f[0]*=oodenom;
  f[1]*=oodenom;
  f[2]*=oodenom;
  /* f0=h0, f1=-i h1, f2=-h2 */
  f[1]=-I*f[1];
  f[2]=-f[2];
  /* signc0 */
  f[0]=creal(f[0])+I*cimag(f[0])*signc0;
  f[1]=creal(f[1])+I*cimag(f[1])*signc0;
  f[2]=creal(f[2])+I*cimag(f[2])*signc0;
  /* make it safe : if c0max<eps -> c0, c1 =0 -> A=0 -> f[0]=1, f[1]=f[2]=0 */
  if(c0max<DBL_MIN) {
    f[0]=1.;
    f[1]=0.;
    f[2]=0.;
  }
#undef C_SQRT3
}/*}}}*/

static inline void M_expA(mat *M,const ahm *A)/*{{{*/
{
  double traa,traaa;
  complex double f[3];
  A_trace23(&traa,&traaa,A);
  exp_trig0(f,traa,traaa);
  M_ct1(M,f[0]);
  A_add_poly(M,A,f[1],f[2]);
}/*}}}*/

/* Matrix exponentialization }}} */

static inline void M_MtM(mat *M,const mat *A,const mat * restrict B)/*{{{*/
  /* A and M may be the same */

{
  for(int a=0;a<3;a++) {
    complex double a0=A->c[a][0];
    complex double a1=A->c[a][1];
    complex double a2=A->c[a][2];
    M->c[a][0]=a0*B->c[0][0]+a1*B->c[1][0]+a2*B->c[2][0];
    M->c[a][1]=a0*B->c[0][1]+a1*B->c[1][1]+a2*B->c[2][1];
    M->c[a][2]=a0*B->c[0][2]+a1*B->c[1][2]+a2*B->c[2][2];
  }
}/*}}}*/

static double retrace_MtM(const mat *A,const mat *B)/*{{{*/
{
  /* the real part of the trace of the product of A and B matrices */
  /* A and M may be the same */
  double r=0;
  r+=creal(
      A->c[0][0]*B->c[0][0]+A->c[0][1]*B->c[1][0]+A->c[0][2]*B->c[2][0]
      );
  r+=creal(
      A->c[1][0]*B->c[0][1]+A->c[1][1]*B->c[1][1]+A->c[1][2]*B->c[2][1]
      );
  r+=creal(
      A->c[2][0]*B->c[0][2]+A->c[2][1]*B->c[1][2]+A->c[2][2]*B->c[2][2]
      );
  return r;
}/*}}}*/

static inline void U_UtU(su3 *R,const su3*A,const su3*B)/*{{{*/
{
  /* product of su3 matrices: R = A * B
   * the three matrices may point to the same object */
  complex double B00=B->c[0][0];
  complex double B01=B->c[0][1];
  complex double B02=B->c[0][2];
  complex double B10=B->c[1][0];
  complex double B11=B->c[1][1];
  complex double B12=B->c[1][2];
  complex double B20=conj(B01*B12-B02*B11); /* reconstruction */
  complex double B21=conj(B02*B10-B00*B12);
  complex double B22=conj(B00*B11-B10*B01);
  complex double a0,a1,a2;
  a0=A->c[0][0];                /* A and R may point to the same object */
  a1=A->c[0][1];
  a2=A->c[0][2];
  R->c[0][0]=a0*B00+a1*B10+a2*B20;
  R->c[0][1]=a0*B01+a1*B11+a2*B21;
  R->c[0][2]=a0*B02+a1*B12+a2*B22;
  a0=A->c[1][0];
  a1=A->c[1][1];
  a2=A->c[1][2];
  R->c[1][0]=a0*B00+a1*B10+a2*B20;
  R->c[1][1]=a0*B01+a1*B11+a2*B21;
  R->c[1][2]=a0*B02+a1*B12+a2*B22;
  /* A's bottom row is not reconstructed */
  /* The bottom row is not calculated. */
}/*}}}*/

static inline void U_UatU(su3 *R,const su3*A,const su3*B)/*{{{*/
{
  /* product of su3 matrices: R = adj(A) * B
   * the three matrices may point to the same object */
  complex double B00=B->c[0][0];
  complex double B01=B->c[0][1];
  complex double B02=B->c[0][2];
  complex double B10=B->c[1][0];
  complex double B11=B->c[1][1];
  complex double B12=B->c[1][2];
  complex double B20=conj(B01*B12-B02*B11);
  complex double B21=conj(B02*B10-B00*B12);
  complex double B22=conj(B00*B11-B10*B01);
  complex double A00=conj(A->c[0][0]);
  complex double A10=conj(A->c[0][1]);
  complex double A20=conj(A->c[0][2]);
  complex double A01=conj(A->c[1][0]);
  complex double A11=conj(A->c[1][1]);
  complex double A21=conj(A->c[1][2]);
  complex double A02=conj(A10*A21-A20*A11);
  complex double A12=conj(A20*A01-A00*A21);
  //    complex double A22=conj(A00*A11-A01*A10); /* this we don't need */
  R->c[0][0]=A00*B00+A01*B10+A02*B20;
  R->c[0][1]=A00*B01+A01*B11+A02*B21;
  R->c[0][2]=A00*B02+A01*B12+A02*B22;
  R->c[1][0]=A10*B00+A11*B10+A12*B20;
  R->c[1][1]=A10*B01+A11*B11+A12*B21;
  R->c[1][2]=A10*B02+A11*B12+A12*B22;
}/*}}}*/

static inline void U_UtUa(su3 *R,const su3*A,const su3*B)/*{{{*/
{
  /* product of su3 matrices: R = A * adj(B)
   * the three matrices may point to the same object */
  complex double B00=conj(B->c[0][0]);
  complex double B10=conj(B->c[0][1]);
  complex double B20=conj(B->c[0][2]);
  complex double B01=conj(B->c[1][0]);
  complex double B11=conj(B->c[1][1]);
  complex double B21=conj(B->c[1][2]);
  complex double B02=conj(B10*B21-B20*B11);
  complex double B12=conj(B20*B01-B00*B21);
  complex double B22=conj(B00*B11-B01*B10);
  complex double a0,a1,a2;
  a0=A->c[0][0];
  a1=A->c[0][1];
  a2=A->c[0][2];
  R->c[0][0]=a0*B00+a1*B10+a2*B20;
  R->c[0][1]=a0*B01+a1*B11+a2*B21;
  R->c[0][2]=a0*B02+a1*B12+a2*B22;
  a0=A->c[1][0];
  a1=A->c[1][1];
  a2=A->c[1][2];
  R->c[1][0]=a0*B00+a1*B10+a2*B20;
  R->c[1][1]=a0*B01+a1*B11+a2*B21;
  R->c[1][2]=a0*B02+a1*B12+a2*B22;
}/*}}}*/

static inline complex double U_trace(const su3 *U)/*{{{*/
{
  /* Trace of an su3 matrix.
   * Only the [2][2] component has to be reconstructed */
  return 
    U->c[0][0]+U->c[1][1]+
    conj(U->c[0][0]*U->c[1][1]-U->c[1][0]*U->c[0][1]);
}/*}}}*/

static inline void Mp_StU(mat *R,double f,const su3 * B)/*{{{*/
{
  /* R = R + f * B */
  /* A typical use the add up tha path products to build a staple. */
  complex double B00=B->c[0][0];
  complex double B01=B->c[0][1];
  complex double B02=B->c[0][2];
  complex double B10=B->c[1][0];
  complex double B11=B->c[1][1];
  complex double B12=B->c[1][2];
  complex double B20=conj(B01*B12-B02*B11);
  complex double B21=conj(B02*B10-B00*B12);
  complex double B22=conj(B00*B11-B10*B01);
  R->c[0][0]+=f*B00;
  R->c[0][1]+=f*B01;
  R->c[0][2]+=f*B02;
  R->c[1][0]+=f*B10;
  R->c[1][1]+=f*B11;
  R->c[1][2]+=f*B12;
  R->c[2][0]+=f*B20;
  R->c[2][1]+=f*B21;
  R->c[2][2]+=f*B22;
}/*}}}*/

static inline void M_0(mat *R)/*{{{*/
{
  /* reset a matrix */
  R->c[0][0]=0;
  R->c[0][1]=0;
  R->c[0][2]=0;
  R->c[1][0]=0;
  R->c[1][1]=0;
  R->c[1][2]=0;
  R->c[2][0]=0;
  R->c[2][1]=0;
  R->c[2][2]=0;
}/*}}}*/

static inline void M_projahm(mat *M)/*{{{*/
{ 
  /* projection to anti hermitian, 
   * but not traceless matrix: A = ( M - adj(M) ) / 2  
   * We cannot use the ahm type here (that would have no trace)
   * */
  ahm a;
  a.c01=c_ahm(M->c[0][1],M->c[1][0]);
  a.c02=c_ahm(M->c[0][2],M->c[2][0]);
  a.c12=c_ahm(M->c[1][2],M->c[2][1]);
  M->c[0][1]=a.c01;
  M->c[0][2]=a.c02;
  M->c[1][2]=a.c12;
  M->c[1][0]=-conj(a.c01);
  M->c[2][0]=-conj(a.c02);
  M->c[2][1]=-conj(a.c12);
  M->c[0][0]=I*cimag(M->c[0][0]);
  M->c[1][1]=I*cimag(M->c[1][1]);
  M->c[2][2]=I*cimag(M->c[2][2]);
}/*}}}*/

static inline void M_StM(mat *D,double f,const mat *S)/*{{{*/
{
  /* multiply a matrix with a scalar: D = f * S */
  complex double m00=S->c[0][0];
  complex double m01=S->c[0][1];
  complex double m02=S->c[0][2];
  complex double m10=S->c[1][0];
  complex double m11=S->c[1][1];
  complex double m12=S->c[1][2];
  complex double m20=S->c[2][0];
  complex double m21=S->c[2][1];
  complex double m22=S->c[2][2];
  D->c[0][0]=f*(m00);
  D->c[0][1]=f*(m01);
  D->c[0][2]=f*(m02);
  D->c[1][0]=f*(m10);
  D->c[1][1]=f*(m11);
  D->c[1][2]=f*(m12);
  D->c[2][0]=f*(m20);
  D->c[2][1]=f*(m21);
  D->c[2][2]=f*(m22);
}/*}}}*/

static inline complex double M_trace(mat *D)/*{{{*/
{
  return D->c[0][0]+D->c[1][1]+D->c[2][2];
}/*}}}*/

/* matrix algebra }}} */

/* lattice and layout {{{ */

enum { dir_x=3,dir_y=2,dir_z=1,dir_t=0 };

/* The 'coord' is simply an array ouf four coordinates.
 * In this program, lattice sites are simply indexed by a
 * single number, the site index. The 'coord' type is seldom used
 * to identifie a site, it is mainly for convenience, where four
 * integers are needed.
 */

typedef int coord[4]; 

/* lattice geometry */

/* The following variables and tables are all initialized by
 * the init_geometry(int lattice_dimension[4]) call. It only needs the
 * lattice size. 
 * You must call init_geometry() before any 
 * lattice field constructor or kernel.
 * So init_geometry() is typically called when you have
 * read in a configuration header.
 */
int Ls;
coord L;
coord Lmul;
/* Neighbour tables.
 * These tables implement the periodic boundary conditions as well as
 * the lattice layout */
coord *neigh_f;  /* neigh_f[i][dir_y]  tells the i-th site's forward
		    forward neighbour in the y direction */
coord *neigh_b; /* neigh_b[i][dir_t]  tells the i-th site's forward
		   backward neighbour in the t direction */

coord *latcoord; /* latcoord[i][3] tells the 3rd coordinate of the site */

char Last_lattice[512];  /* when we read in a configuration, we store its
			    name here for later reference */

/* alloc, release {{{ */

static size_t allocated; /* a record of the previous allocations */

static void *alloc(size_t s) /* s bytes */
{
  void *tmp=malloc(s);
  if(!tmp) {
    fprintf(stderr,
	"Memory allocation error. This attempt: %u bytes."
	"Allocated sofar %u Megabytes",
	(unsigned)s,(unsigned)(allocated/(1024*1024)));
    exit(EXIT_FAILURE);
  }
  allocated+=s;
  return tmp;
}

static void release(void *tmp,size_t s)
{
  free(tmp);
  allocated-=s;
}

/* alloc, release }}} */

/* Lattice field constructuors and destructors {{{ */
void kernel_U_new4(latsu3 U[4])/*{{{*/
{
  size_t one=Ls*sizeof(struct su3);
  char *start=alloc(4*one);
  U[0]=(void*)(start+0*one);
  U[1]=(void*)(start+1*one);
  U[2]=(void*)(start+2*one);
  U[3]=(void*)(start+3*one);
}/*}}}*/

void kernel_M_new4(latmat U[4])/*{{{*/
{
  size_t one=Ls*sizeof(struct mat);
  char *start=alloc(4*one);
  U[0]=(void*)(start+0*one);
  U[1]=(void*)(start+1*one);
  U[2]=(void*)(start+2*one);
  U[3]=(void*)(start+3*one);
}/*}}}*/

void kernel_U_delete4(latsu3 U[4])/*{{{*/
{
  size_t one=Ls*sizeof(struct su3);
  if(U[0]) release(U[0],one*4);
  U[0]=NULL;
  U[1]=NULL;
  U[2]=NULL;
  U[3]=NULL;
}/*}}}*/

void kernel_M_delete4(latmat U[4])/*{{{*/
{
  size_t one=Ls*sizeof(struct mat);
  if(U[0]) release(U[0],one*4);
  U[0]=NULL;
  U[1]=NULL;
  U[2]=NULL;
  U[3]=NULL;
}/*}}}*/

void kernel_A_new4(latahm U[4])/*{{{*/
{
  size_t one=Ls*sizeof(struct ahm);
  char *start=alloc(4*one);
  U[0]=(void*)(start+0*one);
  U[1]=(void*)(start+1*one);
  U[2]=(void*)(start+2*one);
  U[3]=(void*)(start+3*one);
}/*}}}*/

void kernel_A_delete4(latahm U[4])/*{{{*/
{
  size_t one=Ls*sizeof(struct ahm);
  if(U[0]) release(U[0],one*4);
  U[0]=NULL;
  U[1]=NULL;
  U[2]=NULL;
  U[3]=NULL;
}/*}}}*/

latsu3 kernel_U_new(void)/*{{{*/
{
  return alloc(Ls*sizeof(struct su3));
}/*}}}*/

latmat kernel_M_new(void)/*{{{*/
{
  return alloc(Ls*sizeof(struct mat));
}/*}}}*/

latahm kernel_A_new(void)/*{{{*/
{
  return alloc(Ls*sizeof(struct ahm));
}/*}}}*/

void kernel_U_delete(latsu3 ptr)/*{{{*/
{
  release(ptr,Ls*sizeof(struct su3));
}/*}}}*/

void kernel_M_delete(latmat ptr)/*{{{*/
{
  release(ptr,Ls*sizeof(struct mat));
}/*}}}*/

void kernel_A_delete(latahm ptr)/*{{{*/
{
  release(ptr,Ls*sizeof(struct ahm));
}/*}}}*/

/* double precision new and delete }}} */

/* lattice macros + openMP details {{{ */

#define SITE_NULL {0,0,0,0}
#define for_coords(s,i) for(int s[4]=SITE_NULL,i=0;s[0]<L[0];s[0]++) \
				     for(s[1]=0;s[1]<L[1];s[1]++) \
for(s[2]=0;s[2]<L[2];s[2]++) \
for(s[3]=0;s[3]<L[3];s[3]++,i++)

#define for_sites(i) for(int i=0;i<Ls;i++)


#ifdef ENABLE_OMP
#define for_rows(i) for(int nth=omp_get_num_threads(),nr=Ls/L[3],sep=nr/nth*L[3],tid=omp_get_thread_num(),i=tid*sep,endi=tid<nth-1?(tid+1)*sep:Ls;i<endi;i+=L[3])
#else
#define for_rows(i) for(int i=0;i<Ls;i+=L[3])
#endif
#define for_within_row(j) for(int j=0;j<L[3];j++)




#define for_dir(dir) for(int dir=0;dir<4;dir++)

#define site_step(s,step0,step1,step2,step3) { \
  (L[0]+s[0]+step0)%L[0], \
  (L[1]+s[1]+step1)%L[1], \
  (L[2]+s[2]+step2)%L[2], \
  (L[3]+s[3]+step3)%L[3]}

#define site_nind(s,step0,step1,step2,step3) ( \
    ((L[0]+s[0]+step0)%L[0])*Lmul[0]+ \
    ((L[1]+s[1]+step1)%L[1])*Lmul[1]+ \
    ((L[2]+s[2]+step2)%L[2])*Lmul[2]+ \
    ((L[3]+s[3]+step3)%L[3])*Lmul[3])


static inline int site2index(const coord s)
{
  return Lmul[0]*s[0]+Lmul[1]*s[1]+Lmul[2]*s[2]+Lmul[3]*s[3];
}

/* lattice macros }}} */

void init_geometry(int nlat[4])/*{{{*/
{
  if(
      L[dir_t]==nlat[dir_t] &&
      L[dir_z]==nlat[dir_z] &&
      L[dir_y]==nlat[dir_y] &&
      L[dir_x]==nlat[dir_x]
    ) return; /* same as before */
  /* We clean up a bit to forget previous layouts.
   * This is safe even if we are here the first time. */
  if(neigh_f) free(neigh_f);
  if(neigh_b) free(neigh_b);
  if(latcoord) free(latcoord);
  L[dir_t]=nlat[dir_t];
  L[dir_z]=nlat[dir_z];
  L[dir_y]=nlat[dir_y];
  L[dir_x]=nlat[dir_x];
  /* helper's for calculating the site index from coordinates */
  Lmul[3]=1;
  Lmul[2]=L[3]*Lmul[3];
  Lmul[1]=L[2]*Lmul[2];
  Lmul[0]=L[1]*Lmul[1];
  /* the lattice volume */
  Ls=L[0]*L[1]*L[2]*L[3];
  neigh_f=malloc(Ls*4*sizeof(int));
  neigh_b=malloc(Ls*4*sizeof(int));
  latcoord=malloc(Ls*4*sizeof(int));
  if(neigh_b==NULL || neigh_f==NULL || latcoord==NULL) {
    perror("Could not allocate the auxillary lattices");
    exit(EXIT_FAILURE);
  }
  /* store each site's four coordinates */
  for_coords(s,i) {
    latcoord[i][0]=s[0];
    latcoord[i][1]=s[1];
    latcoord[i][2]=s[2];
    latcoord[i][3]=s[3];
  }
  /* store each site's eight neigbours */
  for_coords(s,i)
    for(int dir=0;dir<4;dir++)
    {
      coord z={s[0],s[1],s[2],s[3]};
      z[dir]=(s[dir]+1)%L[dir];
      neigh_f[i][dir]=site2index(z);
      z[dir]=(L[dir]+s[dir]-1)%L[dir];
      neigh_b[i][dir]=site2index(z);
    }
}/*}}}*/

/* lattice }}} */

/* IO and lattice formats {{{ */

/* a pontier-to-function type that specifies the input file format */

typedef int (*configuration_loader)(latsu3 U[4],const char *filename);

/* Endianness {{{ */
/* these routines read n floating point numbers from a file,
 * and convert it to double precision host format
 *
 * To write this part of the code we used some glibc tools
 * to convert between host, big and little endian numbers in 32 and 64 bits.
 * Should these be absent on your system, you can at least
 * replace be32toh by ntohl and htobe32 by htonl. (man htonl)
 * If you know the endianness of your architecture,
 * you can define HOST_LITTLE_ENDIAN or HOST_BIG_ENDIAN 
 * and then these BSD extensions
 * are no longer necessary.
 */

uint32_t byte_swap32(uint32_t x)
{
  /* feel free to replace this with something more efficient */
  union U { unsigned char a[4]; uint32_t x; } u,v;
  u.x=x;
  v.a[0]=u.a[3];
  v.a[1]=u.a[2];
  v.a[2]=u.a[1];
  v.a[3]=u.a[0];
  return v.x;
}

uint64_t byte_swap64(uint64_t x)
{
  /* feel free to replace this with something more efficient */
  union U { unsigned char a[8]; uint64_t x; } u,v;
  u.x=x;
  v.a[0]=u.a[7];
  v.a[1]=u.a[6];
  v.a[2]=u.a[5];
  v.a[3]=u.a[4];
  v.a[4]=u.a[3];
  v.a[5]=u.a[2];
  v.a[6]=u.a[1];
  v.a[7]=u.a[0];
  return v.x;
}

#ifdef HOST_BIG_ENDIAN

#define local_be32toh(x) (x)
#define local_be64toh(x) (x)
#define local_le32toh(x) byte_swap32(x)
#define local_le64toh(x) byte_swap64(x)
#define local_hto32le(x) byte_swap32(x)

#elif defined(HOST_LITTLE_ENDIAN)

#define local_le32toh(x) (x)
#define local_le64toh(x) (x)
#define local_be32toh(x) byte_swap32(x)
#define local_be64toh(x) byte_swap64(x)
#define local_hto32le(x) (x)
#else
#include <endian.h>
#define local_le32toh(x) le32toh(x)
#define local_le64toh(x) le64toh(x)
#define local_be32toh(x) be32toh(x)
#define local_be64toh(x) be64toh(x)
#define local_hto32le(x) htole32(x)
#endif


int read_n_float_32big(double *array,int n,FILE *F)/*{{{*/
{
  union U { uint32_t u; float f; };
  union U *input=malloc(sizeof(union U)*n);
  union U output;
  int r;
  if(!input) { fprintf(stderr,"Allocation failure.\n"); exit(EXIT_FAILURE); }
  r=fread(input,sizeof(union U),n,F);
  for(int i=0;i<r;i++) {
    output.u=local_be32toh(input[i].u);
    array[i]=output.f;
  }
  free(input);
  return r;
}/*}}}*/

int read_n_float_64big(double *array,int n,FILE *F)/*{{{*/
{
  union U { uint64_t u; double f; };
  union U *input=malloc(sizeof(union U)*n);
  union U output;
  int r;
  if(!input) { fprintf(stderr,"Allocation failure.\n"); exit(EXIT_FAILURE); }
  r=fread(input,sizeof(union U),n,F);
  for(int i=0;i<r;i++) {
    output.u=local_be64toh(input[i].u);
    array[i]=output.f;
  }
  free(input);
  return r;
}/*}}}*/

int read_n_float_32little(double *array,int n,FILE *F)/*{{{*/
{
  union U { uint32_t u; float f; };
  union U *input=malloc(sizeof(union U)*n);
  union U output;
  int r;
  if(!input) { fprintf(stderr,"Allocation failure.\n"); exit(EXIT_FAILURE); }
  r=fread(input,sizeof(union U),n,F);
  for(int i=0;i<r;i++) {
    output.u=local_le32toh(input[i].u);
    array[i]=output.f;
  }
  free(input);
  return r;
}/*}}}*/

int read_n_float_64little(double *array,int n,FILE *F)/*{{{*/
{
  union U { uint64_t u; double f; };
  union U *input=malloc(sizeof(union U)*n);
  union U output;
  int r;
  if(!input) { fprintf(stderr,"Allocation failure.\n"); exit(EXIT_FAILURE); }
  r=fread(input,sizeof(union U),n,F);
  for(int i=0;i<r;i++) {
    output.u=local_le64toh(input[i].u);
    array[i]=output.f;
  }
  free(input);
  return r;
}/*}}}*/

/* Endianness }}} */

/* Link storage format {{{ */


/* Here we collected the basic link read and write routines
 * they accept the byte-order I/O routines as an argument
 * (e.g. read_n_float_32big)
 * Each of these functions do one matrix at a time. 
 * These functions return 0 on success, and -1 on error.
 *
 */

int read_matrix_2x3( mat *m,int (*read_n_float)(double *,int,FILE *),FILE *F)/*{{{*/
{
  /* our su3 format happens to be identical to this storage format */
  su3 u;
  if(read_n_float((double*)&u.c[0][0],12,F)!=12) return -1;
  M_U(m,&u);
  return 0;
}/*}}}*/

int read_matrix_3x3(mat *m,int (*read_n_float)(double *,int,FILE *),FILE *F)/*{{{*/
{
  /* our matrix format happens to be identical to this storage format */
  if(read_n_float((double*)&m->c[0][0],18,F)!=18) return -1;
  return 0;
}/*}}}*/

/* Link storage format }}}  */

int load_configuration_nersc(latsu3 U[4],const char *filename)/*{{{*/
{
  FILE *F;
  int e=0;                                      /* error code */
  int nlat[4]={0,0,0,0};                        /* geometry in the header */
  int store_order[4]={dir_x,dir_y,dir_z,dir_t}; /* link ordering */
  char line[128];                               /* buffer for one header line */
  /* for the storage of some header info: */
  char float_type[16];                        
  char data_type[24];
  unsigned checksum=0;
  double plaquette=0,link_trace=0;
  int (*read_n_float)(double *array,int n,FILE *F)=read_n_float_32big;
  int (*read_matrix)(mat *m,
      int (*read_n_float)(double *array,int n,FILE *F),
      FILE *F)=read_matrix_2x3;
  /* temporary variables for header values: */
  double p;     
  int pi,pj;
  /* these will be used later for tests */
  double true_link_trace=0.;
  double true_plaquette=0;
  /* function headers: */
  double normalized_plaquette(latsu3 U[4]);  /* to check the plaquette */
  /* we reset the pointers to Null, because if anything goes wrong, a
   * kernel_U_delete will be attepmted, and that is safe to call with NULL */
  for(int mu=0;mu<4;mu++) U[mu]=NULL;
  /* We open the file. Is it there at all ? */
  F=fopen(filename,"r");
  if(!F) {
    fprintf(stderr,"Could not open %s:%s\n",filename,strerror(errno));
    return -1;
  }
  /* NERSC header {{{ */
  /* the first line we check separately */
  if(fgets(line,sizeof(line),F)==NULL || strcmp(line,"BEGIN_HEADER\n")!=0) {
    fprintf(stderr,"Missing or invalid header in %s\n",filename);
    e=-1;
    goto finish;
  }
  /* then we read in lines in a loop until the header ends */
  while(fgets(line,sizeof(line),F)) {
    if(strcmp(line,"END_HEADER\n")==0) 
      break; /* binary part starts, we exit from this loop */
    if(sscanf(line,"CHECKSUM = %x",&checksum)==1) continue;
    if(sscanf(line,"LINK_TRACE = %lf",&link_trace)==1) continue;
    if(sscanf(line,"PLAQUETTE = %lf",&plaquette)==1) continue;
    if(sscanf(line,"HDR_VERSION = %lf",&p)==1) {
      if(p==1.0) continue;
      else {
	fprintf(stderr,"Unsupported header version in %s\n",filename);
	e=-1;
	goto finish;
      }
    }
    if(sscanf(line,"STORAGE_FORMAT = %lf",&p)==1) {
      if(p==1.0) continue;
      else {
	fprintf(stderr,"Unsupported format version in %s\n",filename);
	e=-1;
	goto finish;
      }
    }
    if(sscanf(line,"DIMENSION_%d = %d",&pi,&pj)==2) {
      if(pi>0 && pi<5 && pj>0) nlat[store_order[pi-1]]=pj;
      else {
	fprintf(stderr,"Invalid lattice dimension in %s\n",filename);
	e=-1;
	goto finish;
      }
    }
    if(sscanf(line,"FLOATING_POINT = %15s",float_type)==1) {
      if(!strcmp(float_type,"IEEE32BIG"))
	read_n_float=read_n_float_32big;
      else if(!strcmp(float_type,"IEEE32LITTLE"))
	read_n_float=read_n_float_32little;
      else if(!strcmp(float_type,"IEEE64BIG"))
	read_n_float=read_n_float_64big;
      else if(!strcmp(float_type,"IEEE64LITTLE"))
	read_n_float=read_n_float_64little;
      else {
	fprintf(stderr,"Invalid float format in %s\n",filename);
	e=-1;
	goto finish;
      }
    }
    if(sscanf(line,"DATATYPE = %23s",data_type)==1) {
      for(int i=0;data_type[i];i++) data_type[i]=toupper(data_type[i]);
      if(!strcmp(data_type,"4D_SU3_GAUGE"))
	read_matrix=read_matrix_2x3;
      else if(!strcmp(data_type,"4D_SU3_GAUGE_3X3"))
	read_matrix=read_matrix_3x3;
      else {
	fprintf(stderr,"Invalid matrix format in %s\n",filename);
	fprintf(stderr,"Format found %s\n",data_type);
	e=-1;
	goto finish;
      }
    }
  }
  if(nlat[0]*nlat[1]*nlat[2]*nlat[3]==0 || checksum==0) {
    fprintf(stderr,"Some header information are missing from %s\n",filename);
    e=-1;
    goto finish;
  }
  /* NERSC header }}} */
  init_geometry(nlat);
  kernel_U_new4(U);
  /* The main read loop */
  for_sites(s) {
    for(int dir=0;dir<4;dir++) {
      int mu=store_order[dir];
      mat m;
      if(read_matrix(&m,read_n_float,F)) {
	fprintf(stderr,"Short read in configuration %s\n", filename);
	e=-1;
	goto finish;
      }
      U_M(&U[mu][s],&m);
      true_link_trace+=creal(M_trace(&m));
    }
  }
  /* checks {{{ */
  true_link_trace/=12.*Ls;
  true_plaquette=normalized_plaquette(U);
  if(!isfinite(true_link_trace) || !isfinite(true_plaquette)) {
    fprintf(stderr,"The configuration is full of nans. Is the byte-order correct?\n");
#ifdef  HOST_BIG_ENDIAN
    fprintf(stderr,"HOST_BIG_ENDIAN is set\n" ) ;
#endif
#ifdef  HOST_LITTLE_ENDIAN
    fprintf(stderr,"HOST_LITTLE_ENDIAN is set\n" ) ;
#endif
    exit(EXIT_FAILURE);
  }
  if(fabs(true_link_trace-link_trace)/(true_link_trace+link_trace) > 1e-5) {
    fprintf(stderr,"The link_trace is wrong in %s\n", filename);
    e=-1;
    goto finish;
  } 
  //fprintf(stderr,"Measured plaquette = %g\n", true_plaquette);
fprintf(stdout,"Measured plaquette = %g\n", true_plaquette);
  if(fabs(true_plaquette-plaquette)/(true_plaquette+plaquette) > 1e-5) {
    fprintf(stderr,"The plaquette is wrong in %s\n", filename);
    e=-1;
    goto finish;
  } 
  /* checks }}} */
finish:
  if(e) {
    kernel_U_delete4(U);
    Last_lattice[0]=0;
  } else {
    strncpy(Last_lattice,filename,sizeof(Last_lattice)-1);
    Last_lattice[sizeof(Last_lattice)-1]=0;
  }
  return e;
}/*}}}*/

configuration_loader format_selector(const char *format_name)/*{{{*/
{
  if(!strcmp(format_name,"nersc")) return load_configuration_nersc;
  /* add your favourite format here 
   * else if(!strcmp(format_name,"my-format")) return load_configuration_my;
   * */
  else {
    fprintf(stderr,"Unsupported lattice format: %s\n",format_name);
    exit(EXIT_FAILURE);
  }
}/*}}}*/

/* IO }}} */

/* kernels {{{ */

/* Kernels are basic lattice operations that involve a
 * loop over lattice sites. */
/* As a general rule, the field in the first argument is modified,
 * while the others are inputs. 
 * (Exception: kernel_flow_update has more outputs )
 * The naming convention is similar to that of the individiual matrix routines.
 */

/* a compact list of what they do 
 *
 * kernel_staple_plaquette(latmat S[4],latsu3 U[4])
 *              staple in all direction.
 *              The path is ordered such that 
 *                  tr(U[mu]*S[mu]) is gauge invariant for any mu
 *
 * kernel_staple_symanzik(latmat S[4],latsu3 U[4])
 *              staple for the Symanzik action in all direction.
 *              The path is ordered such that 
 *                  tr(U[mu]*S[mu]) is gauge invariant for any mu
 *
 * kernel_staple(latmat S[4],latsu3 U[4])
 *              either kernel_staple or kernel_staple_symanzik,
 *              depending on the program parameters
 *
 * void kernel_tU_expStA(latsu3 U,double f,latahm A)
 *              U = exp ( f * A ) * U           (for all sites)
 *
 * void kernel_A_StProjM(latahm R,double f,latmat B)      
 *              A = f * Project_antihermitian_traceless(B)
 *
 * void kernel_Ap_StProjM(latahm R,double f,latmat B)      
 *              A = A + f * Project_antihermitian_traceless(B)
 *
 * void kernel_A_0(latahm D)
 *              reset D at all sites
 *
 * void kernel_M_0(latmat D)
 *              reset D at all sites
 *
 * void kernel_M_UtM(latmat R,latsu3 U,latmat B)
 *              M = reconstruct(U) * B  (for all sites)
 * 
 * void kernel_flow_update(latsu3 U,latahm A,latmat S,double f1,double f2)
 *              A = A + f1 * Project_antihermitian_traceless(U*S)
 *              U = exp( f2 * A ) * U
 *              S may be modified 
 *
 * void kernel_fmunu(latmat Fmn,latsu3 U[4],int mu,int nu)
 *              Clover representation of the Fmunu matrix,
 *              for one choiche of mu and nu (=0..3), for all sites.
 *              See Fig. 1 in arXiv:1006.4518
 *
 *
 *
 */

void kernel_staple_plaquette(latmat S[4],latsu3 U[4])/*{{{*/
{
  const double coeff[4]={Xi*Xi,1.0,1.0,1.0};
#ifdef ENABLE_OMP
#pragma omp parallel 
#endif
  {
    for_rows(r) for_dir(mu) {
      for_within_row(x) M_0(&S[mu][r+x]);
      for_dir(nu) if(nu!=mu) for_within_row(x) {
	int s=x+r;
	su3 u;
	su3 v;
	int s1,s2;
	/* upper staple */
	s1=neigh_f[s][mu];
	s2=neigh_f[s][nu];
	U_UtUa(&u,&U[nu][s1],&U[mu][s2]);
	U_UtUa(&v,&u,&U[nu][s]);
	Mp_StU(&S[mu][s],coeff[nu],&v);
	/* lower staple */
	s1=neigh_b[s1][nu];
	s2=neigh_b[s][nu];
	U_UtU(&u,&U[mu][s2],&U[nu][s1]);
	U_UatU(&v,&u,&U[nu][s2]);
	Mp_StU(&S[mu][s],coeff[nu],&v);
      }
    }
  }
}/*}}}*/

void kernel_staple_symanzik(latmat S[4],latsu3 U[4])/*{{{*/
{
  double fac_1x1[4]={Xi*Xi,1.0,1.0,1.0};
  double fac_1x2[4]={Xi*Xi,1.0,1.0,1.0};
  for(int dir=0;dir<4;dir++) fac_1x1[dir]*=5./3.;
  for(int dir=0;dir<4;dir++) fac_1x2[dir]*=-1./12.;
#if ENABLE_OMP
#pragma omp parallel
#endif
  {
    for_rows(r) for_dir(mu) {
      for_within_row(x) M_0(&S[mu][r+x]);
      for_dir(nu) if(nu!=mu) 
	for_within_row(x) {
	  int s=r+x;
	  su3 u;
	  su3 v;
	  int sf0=neigh_f[s][mu];
	  int s0f=neigh_f[s][nu];
	  int sfb=neigh_b[sf0][nu];
	  int s0b=neigh_b[s][nu];
	  int sF0=neigh_f[sf0][mu];
	  int s0F=neigh_f[s0f][nu];
	  int sff=neigh_f[sf0][nu];
	  int sbf=neigh_b[s0f][mu];
	  int sb0=neigh_b[s][mu];
	  int sFb=neigh_f[sfb][mu];
	  int s0B=neigh_b[s0b][nu];
	  int sfB=neigh_b[sfb][nu];
	  int sbb=neigh_b[s0b][mu];
	  /* upper staple */
	  U_UtUa(&u,&U[nu][sf0],&U[mu][s0f]);
	  U_UtUa(&v,&u,&U[nu][s]);
	  Mp_StU(&S[mu][s],fac_1x1[nu],&v);
	  /* 1x2 staples */
	  /* mu,nu,-mu,-mu,-nu */
	  U_UtU(&u,&U[mu][sf0],&U[nu][sF0]);
	  U_UtUa(&v,&u,&U[mu][sff]);
	  U_UtUa(&u,&v,&U[mu][s0f]);
	  U_UtUa(&v,&u,&U[nu][s]);
	  Mp_StU(&S[mu][s],fac_1x2[nu],&v);
	  /* nu,nu,-mu,-nu,-nu */
	  U_UtU(&u,&U[nu][sf0],&U[nu][sff]);
	  U_UtUa(&v,&u,&U[mu][s0F]);
	  U_UtUa(&u,&v,&U[nu][s0f]);
	  U_UtUa(&v,&u,&U[nu][s]);
	  Mp_StU(&S[mu][s],fac_1x2[nu],&v);
	  /* nu,-mu,-mu,-nu,mu */
	  U_UtUa(&u,&U[nu][sf0],&U[mu][s0f]);
	  U_UtUa(&v,&u,&U[mu][sbf]);
	  U_UtUa(&u,&v,&U[nu][sb0]);
	  U_UtU(&v,&u,&U[mu][sb0]);
	  Mp_StU(&S[mu][s],fac_1x2[nu],&v);
	  /* lower staple */
	  U_UtU(&u,&U[mu][s0b],&U[nu][sfb]);
	  U_UatU(&v,&u,&U[nu][s0b]);
	  Mp_StU(&S[mu][s],fac_1x1[nu],&v);
	  /* mu,-nu,-mu,-mu,nu */
	  U_UtUa(&u,&U[mu][sf0],&U[nu][sFb]);
	  U_UtUa(&v,&u,&U[mu][sfb]);
	  U_UtUa(&u,&v,&U[mu][s0b]);
	  U_UtU(&v,&u,&U[nu][s0b]);
	  Mp_StU(&S[mu][s],fac_1x2[nu],&v);
	  /* -nu,-nu,-mu,nu,nu */
	  U_UtU(&u,&U[mu][s0B],&U[nu][sfB]);
	  U_UtU(&v,&u,&U[nu][sfb]);
	  U_UatU(&u,&v,&U[nu][s0B]);
	  U_UtU(&v,&u,&U[nu][s0b]);
	  Mp_StU(&S[mu][s],fac_1x2[nu],&v);
	  /* -nu,-mu,-mu,nu,mu */
	  U_UtU(&u,&U[mu][sbb],&U[mu][s0b]);
	  U_UtU(&v,&u,&U[nu][sfb]);
	  U_UatU(&u,&v,&U[nu][sbb]);
	  U_UtU(&v,&u,&U[mu][sb0]);
	  Mp_StU(&S[mu][s],fac_1x2[nu],&v);
	}
    }
  }
}/*}}}*/

void (*kernel_staple)(latmat S[4],latsu3 U[4])=kernel_staple_plaquette;

void kernel_tU_expStA(latsu3 U,double f,latahm A)/*{{{*/
{
#ifdef ENABLE_OMP
#pragma omp parallel
#endif
  {
    for_rows(s)  {
      for_within_row(x) {
	mat m1,m2,m3;
	ahm sta;
	A_StA(&sta,f,&A[s+x]);
	M_expA(&m1,&sta);
	M_U(&m2,&U[s+x]);
	M_MtM(&m3,&m1,&m2);
	U_M(&U[s+x],&m3);
      }
    }
  }
}/*}}}*/

void kernel_A_StProjM(latahm R,double b,latmat B)/*{{{*/
{
#ifdef ENABLE_OMP
#pragma omp parallel
#endif
  {
    for_rows(r) for_within_row(x) { int s=x+r;
      ahm a;
      A_M(&a,&B[s]);
      A_StA(&R[s],b,&a);
    }
  }
}/*}}}*/

void kernel_Ap_StProjM(latahm R,double b,latmat B)/*{{{*/
{
#ifdef ENABLE_OMP
#pragma omp parallel
#endif
  {
    for_rows(r) for_within_row(x) { int s=x+r;
      ahm a;
      A_M(&a,&B[s]);
      Ap_StA(&R[s],b,&a);
    }
  }
}/*}}}*/

void kernel_A_0(latahm D)/*{{{*/
{
  memset(D,0,sizeof(ahm)*Ls);
}/*}}}*/

void kernel_M_0(latmat D)/*{{{*/
{
  memset(D,0,sizeof(mat)*Ls);
}/*}}}*/

void kernel_M_UtM(latmat R,latsu3 U,latmat B)/*{{{*/
{
#ifdef ENABLE_OMP
#pragma omp parallel
#endif
  {
    for_rows(r) for_within_row(x) { int s=x+r;
      mat a,b;
      M_U(&a,&U[s]);
      M_M(&b,&B[s]);
      M_MtM(&R[s],&a,&b);
    }
  }
}/*}}}*/

void kernel_flow_update(latsu3 U,latahm A,latmat S,double f1,double f2)/*{{{*/
{
#if 0
  /* less optimized version */
  kernel_M_UtM(S,U,S);
  kernel_Ap_StProjM(A,f1,S);
  kernel_tU_expStA(U,f2,A);
#else
  /* more optimized version */
  /* but S is not changed here */
#ifdef ENABLE_OMP
#pragma omp parallel
#endif
  {
    for_rows(r) for_within_row(x) { int s=x+r;
      mat a,b,p;
      mat m1,m2,m3;
      ahm ap;
      M_U(&a,&U[s]);
      M_M(&b,&S[s]);
      M_MtM(&p,&a,&b);
      A_M(&ap,&p);
      Ap_StA(&A[s],f1,&ap);
      A_StA(&ap,f2,&A[s]);
      M_expA(&m1,&ap);
      M_U(&m2,&U[s]);
      M_MtM(&m3,&m1,&m2);
      U_M(&U[s],&m3);
    }
  }
#endif
}/*}}}*/

void kernel_fmunu(latmat Fmn,latsu3 U[4],int mu,int nu)/*{{{*/
{
#ifdef ENABLE_OMP
#pragma omp parallel
#endif
  {
    for_rows(r) for_within_row(x) {
      int s=r+x;
      int s1,s2,s3;
      mat m;
      su3 u,v;
      /* 1st */
      s1=neigh_f[s][mu];
      s2=neigh_f[s][nu];
      U_UtU(&u,&U[mu][s],&U[nu][s1]);
      U_UtUa(&v,&u,&U[mu][s2]);
      U_UtUa(&u,&v,&U[nu][s]);
      M_U(&m,&u);
      /* 2nd */
      s1=neigh_b[s][mu];
      s2=neigh_b[s1][nu];
      s3=neigh_b[s][nu];
      U_UtU(&u,&U[nu][s2],&U[mu][s1]);
      U_UatU(&v,&u,&U[mu][s2]);
      U_UtU(&u,&v,&U[nu][s3]);
      Mp_StU(&m,1.0,&u);
      /* 3rd */
      s2=neigh_b[s][mu];
      s1=neigh_f[s2][nu];
      U_UtUa(&u,&U[nu][s],&U[mu][s1]);
      U_UtUa(&v,&u,&U[nu][s2]);
      U_UtU(&u,&v,&U[mu][s2]);
      Mp_StU(&m,1.0,&u);
      /* 4th */
      s1=neigh_b[s][nu];
      s2=neigh_f[s1][mu];
      U_UatU(&u,&U[nu][s1],&U[mu][s1]);
      U_UtU(&v,&u,&U[nu][s2]);
      U_UtUa(&u,&v,&U[mu][s]);
      Mp_StU(&m,1.0,&u);
      /* projecticon and factors */
      M_projahm(&m);
      M_StM(&Fmn[s],0.25,&m);
    }
  }
}/*}}}*/

void kernel_plaquette(double *plaq_t,double *plaq_s,latsu3 U[4])/*{{{*/
{
  double sum_s=0;
  double sum_t=0;
#if ENABLE_OMP
#pragma omp parallel reduction(+ :sum_s,sum_t)
#endif
  {
    for_rows(r) for_within_row(x) { int s=x+r;
      for(int mu=1;mu<4;mu++) for(int nu=mu+1;nu<4;nu++) {
	int s1,s2;
	su3 u,v;
	s1=neigh_f[s][mu];
	s2=neigh_f[s][nu];
	U_UtUa(&u,&U[nu][s1],&U[mu][s2]);
	U_UtUa(&v,&u,&U[nu][s]);
	U_UtU(&u,&v,&U[mu][s]);
	sum_s+=creal(U_trace(&u));
      }
    }
    for_rows(r) for_within_row(x) { int s=x+r;
      for(int mu=1;mu<4;mu++) {
	int nu=dir_t;
	int s1,s2;
	su3 u,v;
	s1=neigh_f[s][mu];
	s2=neigh_f[s][nu];
	U_UtUa(&u,&U[nu][s1],&U[mu][s2]);
	U_UtUa(&v,&u,&U[nu][s]);
	U_UtU(&u,&v,&U[mu][s]);
	sum_t+=creal(U_trace(&u));
      }
    }
  }
  *plaq_t=sum_t/(3.*Ls);
  *plaq_s=sum_s/(3.*Ls);
}/*}}}*/

complex double kernel_trace_MtM(latmat M1,latmat M2)/*{{{*/
{
  complex double n=0;
#ifdef ENABLE_OMP
#pragma omp parallel reduction(+ : n)
#endif
  {
    for_rows(r) for_within_row(x) { int s=x+r;
      mat c;
      M_MtM(&c,&M1[s],&M2[s]);
      n+=M_trace(&c);
    }
  }
  return n/Ls;
}/*}}}*/

double kernel_retrace_MtM(latmat M1,latmat M2)/*{{{*/
{
  double n=0;
#ifdef ENABLE_OMP
#pragma omp parallel reduction(+ : n)
#endif
  {
    for_rows(r) for_within_row(x) { int s=x+r;
      n+=retrace_MtM(&M1[s],&M2[s]);
    }
  }
  return n/Ls;
}/*}}}*/

/* kernels }}} */

/* high level tasks {{{ */

/* I think most of the readers of this program already have some
 * implementation of the kernels and the basic linear algebra on the
 * hardware that he or she controls. Then, these high level routines
 * should be easily portable to the user's favourite environment.
 */

void stout_smearing(latsu3 U[4],double rho)/*{{{*/
{
  /* The Wilson flow is 'in principle' a series of infinitesimal 
   * stout smeraings. Performing lot's of stout smearing is, however,
   * a very sub-optimal integrator (like the Euler method for differential
   * equations). Neventheless, we give the plain (in-place) 
   * stout_smearing routine for optional use and experimentation. */
  latmat S[4];
  latahm A;
  kernel_M_new4(S);
  A=kernel_A_new();
  kernel_staple_plaquette(S,U);
  for(int mu=0;mu<4;mu++) {
    kernel_M_UtM(S[mu],U[mu],S[mu]);
    kernel_A_StProjM(A,-rho,S[mu]);
    kernel_tU_expStA(U[mu],1.0,A);
  }
  kernel_A_delete(A);
  kernel_M_delete4(S);
}/*}}}*/

void stout_step_rk(latsu3 U[4],double eps)/*{{{*/
{
  /* This is the actual implementation of
   * Luscher's integrator in arXiv:1006.4518, appendix C */
  latmat S[4];
  latahm A[4];
  kernel_M_new4(S);
  kernel_A_new4(A);
  for(int mu=0;mu<4;mu++) kernel_A_0(A[mu]);
  kernel_staple(S,U);
  for(int mu=0;mu<4;mu++) 
    kernel_flow_update(U[mu],A[mu],S[mu],+17./36.*eps,-9./17.);
  kernel_staple(S,U);
  for(int mu=0;mu<4;mu++) 
    kernel_flow_update(U[mu],A[mu],S[mu],-8./9.*eps,1.0);
  kernel_staple(S,U);
  for(int mu=0;mu<4;mu++) 
    kernel_flow_update(U[mu],A[mu],S[mu],3./4.*eps,-1.0);
  kernel_A_delete4(A);
  kernel_M_delete4(S);
}/*}}}*/

double normalized_plaquette(latsu3 U[4])/*{{{*/
{
  /* nothing more than adding up the plaquette compontents */
  double pt,ps;
  kernel_plaquette(&pt,&ps,U);
  return (pt+ps)/6.;
}/*}}}*/

void generic_fmunufmunu(double *temporal,double *spatial,double *charge,latsu3 U[4])/*{{{*/
{
  /* this routine calculates the temporal and spatial component of the
   * field strength operator. We also calculate the topological charge. */
  latmat fmn1,fmn2;
  fmn1=kernel_M_new();
  fmn2=kernel_M_new();
  /* We calculate Fmunu and Fmunu_tilde, and get the traces of the
   * products that one can form. */
  /* epsilon_0123= 1 */
  kernel_fmunu(fmn1,U,0,1);
  kernel_fmunu(fmn2,U,2,3);
  *charge=kernel_retrace_MtM(fmn1,fmn2);
  *temporal=-kernel_retrace_MtM(fmn1,fmn1);
  *spatial=-kernel_retrace_MtM(fmn2,fmn2);
  /* epsilon_0231= 1 */
  kernel_fmunu(fmn1,U,0,2);
  kernel_fmunu(fmn2,U,3,1);
  *charge+=kernel_retrace_MtM(fmn1,fmn2);
  *temporal-=kernel_retrace_MtM(fmn1,fmn1);
  *spatial-=kernel_retrace_MtM(fmn2,fmn2);
  /* epsilon_0312= 1 */
  kernel_fmunu(fmn1,U,0,3);
  kernel_fmunu(fmn2,U,1,2);
  *charge+=kernel_retrace_MtM(fmn1,fmn2);
  *temporal-=kernel_retrace_MtM(fmn1,fmn1);
  *spatial-=kernel_retrace_MtM(fmn2,fmn2);
  /* norm */
  *charge*=-Ls*.02533029591058444286; /* -1/4pi^2 */
  kernel_M_delete(fmn1);
  kernel_M_delete(fmn2);
}/*}}}*/

void wflow_destructive(latsu3 U[4],double xi,double tmax,double dt,const char *latname)/*{{{*/
{
  /* This function implements the outer loop for the wilson flow
   * integrator. It controlls the output as the flow progresses.
   * The lattice links U[4] are overwritten.
   *
   * t goes from zero to tmax in dt steps.
   * latname is the file name of the configuration.
   * xi is the anisotropy that goes into the gauge action that defines
   * the flow.
   *
   */
  /* The output consists in five coulmns:
   * 1) a label "WFLOW"
   * 2) flow time "t"
   * 3) the temperal component of the field strength "Et"
   * 4) the spatial component of the field strength "Es"
   * 5) the topological charge
   *
   * Et and Es are given separately so that one can set the
   * spatial and temporal lattice spacing on an anisotropic lattice.
   *
   * The dimensionless wilson flow of the field strength is then
   * f(t) = t*t*(Es+Et)
   * The w0 scale is defined as
   * w0 = sqrt(t) where t * df/dt = 0.3
   * Our best estimate (Budapest-Marseille-Wuppertal collab, 2012): 
   *      w0 = 0.1755 fm
   */
  char flow_name[512];
  double t=0;                 /* flow time */
  double Et,Es,charge;        /* this is the output */
  double old_value=0,new_value=0; /* t^2*E, so thet we know where we are */
  double der_value=0;         /* t * d( t^2E) / dt */
  int stepi=0;
  FILE *O;
  snprintf(flow_name,sizeof(flow_name),"flow.%s",latname);
  O=fopen(flow_name,"w");
  if(!O) {
    fprintf(stderr,"Could not create %s: %s",flow_name,strerror(errno));
    exit(EXIT_FAILURE);
  }
  fprintf(O,"#WFLOW time gact4i gactij topcharge\n");
  /* tmax == 0 parameter tells that the user has no idea how long
   * the flow should be. In this case we try to figure out ourselves,
   * such that flow is long enough to determine the w0 scale. */
  for(stepi=0;tmax==0 || t<tmax-dt/2;stepi++) {
    stout_step_rk(U,dt); 
    t+=dt;
    /* fmunufmunu does not know about the anisotropy ! */
    generic_fmunufmunu(&Et,&Es,&charge,U);
    /* the variables necessary for scale setting: used if tmax==0 only */
    old_value=new_value;
    new_value=t*t*(Et*xi*xi+Es); /* the anisotropic formula */
    der_value=t*(new_value-old_value)/dt;
    /* The output to is printed and is directed to a file as well. */
    fprintf(O,"WFLOW %g %g %g %g\n",t,Et,Es,charge);
    printf("WFLOW %g %g %g %g\n",t,Et,Es,charge);
    /* Now, for the tmax == 0 case, we figure out if we can stop: */
    /* we never stop before t=1 (the lattice spacing scale) */
    /* The stepi%20==0 condition suppresses the fluctuation in the
     * flow length when whole ensemble is analyzed with tmax=0 */
    if(tmax==0 && t>1 && stepi%20==0 ) { 
      if(new_value > 0.45 && der_value > 0.35) break;
    }
  }
  fclose(O);
}/*}}}*/

/* high level tasks }}} */

void run_info()/*{{{*/
{
#ifdef ENABLE_OMP
  int nthreads=1;
#pragma omp parallel
  {
    nthreads = omp_get_num_threads();
  }
  printf("Running with OMP number of threads = %d\n", nthreads);
#endif
}/*}}}*/

int main(int argc,char **argv)/*{{{*/
{
  /* some defaults */
  double xi=1.;        /* Anisotropy of the flow equation as/at. */
  double epsilon=0.02; /* 0.01 is safe, 0.02 is mostly OK, >0.1 is not OK */
  double tmax=0;       /* 0: the program figures out how long the flow
			  should be for a w0-scale determination. */
  run_info();
  /* the format that we understand: */
  configuration_loader load_configuration=load_configuration_nersc;
  /* parameter parsing {{{ */
  if ( argc < 2 ) /* need command line arguments **/ { help() ; exit(0) ; }
  for(int o;(o=getopt(argc,argv,"he:t:x:sf:"))!=-1;) switch(o) {
    case 'e': epsilon=atof(optarg); break;
    case 't': tmax=atof(optarg); break;
    case 'x': xi=atof(optarg); break; /* smearing anisotropy */
    case 's': kernel_staple=kernel_staple_symanzik; break;
    case 'f': load_configuration=format_selector(optarg); break;
    case 'h': return help();
    default: printf("-h gives a help\n"); return EXIT_FAILURE;
  }
  /* parameter parsing }}} */
  /* parameter validation {{{ */
  if(Xi<1) 
    fprintf(stderr,"Warning: The anisotropy in -x <Xi> is Xi=as/at."
	" Is your number Xi=%g correct?",Xi);
  if(tmax<0) {
    fprintf(stderr,"Error: Negative flow time? I cannot unsmear, sorry!\n");
    return EXIT_FAILURE;
  }
  if(epsilon>0.1 || epsilon <= 0.0) {
    fprintf(stderr,"Error: Invalid flow time step.\n");
    return EXIT_FAILURE;
  }
  /* parameter validation }}} */
  /* loop over configurations, specified as arguments {{{ */
  for(int i=optind;i<argc;i++) {
    /* load, analyze and forget one configuration */
    latsu3 U[4];
    if(load_configuration(U,argv[i])) continue;
    printf("%s\t%dx%dx%dx%d\n",argv[i],L[dir_x],L[dir_y],L[dir_z],L[dir_t]);
    wflow_destructive(U,xi,tmax,epsilon,argv[i]);
    kernel_U_delete4(U);
  }
  /* loop over configurations, specified as arguments }}} */
  return 0;
}/*}}}*/
