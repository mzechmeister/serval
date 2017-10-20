#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <dlfcn.h>

#define MIN(a,b) ((a<b)?a:b)

/*
   gcc -c  -Wall -O2 -ansi -pedantic -fPIC cbspline.c; gcc -o cbspline.so -shared cbspline.o
*/


void cbsplcoeff(double *x, double *G, long *kk, int n) {
/*
      kk,pp = divmod(x,1)
      qq = 1-pp
      G[0,:]= qq**3
      G[1,:]= (3*pp**3 - 6*pp**2 + 4) # (-3*pp**2 (1-pp)  -3*pp**2 + 4)
      G[2,:]= (3*qq**3 - 6*qq**2 + 4) #= (-3*pp**3 + 3*pp**2 + 3*pp + 1)
      G[3,:]= pp**3
      G /= 6
*/
   long i, i4;
   double pp, qq, tmp;
   for (i=0; i<n; ++i) {
      i4 = 4 * i;
      kk[i] = (long)x[i];    /* integer part */
      pp = x[i] - kk[i];     /* fractional part */
      qq = 1 - pp;
      tmp = pp*pp; pp *= tmp;
      G[i4+3] = pp/6;
      G[i4+1] = (3*pp +4)/6 - tmp;
      tmp = qq*qq; qq *= tmp;
      G[i4]   = qq/6;
      G[i4+2] = (3*qq +4)/6 - tmp;
   }
}


int rhs_fill(double *rhs, double *G, double *w, long *kk, int n) {
/*      for i,kki in enumerate(kk):
         rhs[kk[i]:kk[i]+4] += G[:,i]*wy[i] # rhs[-1] = 0 !!!!
*/
   int  j;
   long i, in;

   for (i=0; i<n; ++i) {
      in = i * 4;
      for (j=0; j<4; ++j)
         rhs[kk[i]+j] += w[i] * G[in+j];
   }

   return 1;
}


int lhsbnd_fill(double *lhsbnd, double *G, double *w, long *kk, int n, int bw, int nk, double lam) {
  /* fill  lhs = w*G*G
   lhsbnd - band matrix ((bw x kkmax))
   n - number of data points
   G - matrix with bspline coefficients of data points ((bw/2 x n))
   kk - reference knots of data points ((n))
   w - weights ((n))

   python:
   bklind = 3+np.arange(4)+(6*np.arange(4))[:,np.newaxis]
   for i,kki in enumerate(kk): lhsbnd.flat[kk[i]*7+bklind] += w[i]*G[:,i]*G[:,i,np.newaxis]
   */
   int  j, k;
   long i, in, kkk, kkmax=0;  /* kkmax = dim lhsbnd*/
   double tmp;

   for (i=0; i<n; ++i) {
      in = i * 4;
      kkk = kk[i] + 3*nk;
      if (kk[i]>kkmax) kkmax = kk[i];
      /*printf("%ld %ld %ld\n",i,kk[i],in);*/

      for (j=0; j<4; ++j)
         for (k=0; k<=j; ++k)
            lhsbnd[kkk-nk*(j-k)+j] += w[i] * G[in+j] * G[in+k];
      /*printf("%ld %ld %f\n",nk*bw,k,lhsbnd[kkk+(nk-1)*k+j]);*/
   }
   kkmax += 3;

/*
array([[ 1., -2.,  1.,  0.,  0.],
       [-2.,  5., -4.,  1.,  0.],
       [ 1., -4.,  6., -4.,  1.],
       [ 0.,  1., -4.,  5., -2.],
       [ 0.,  0.,  1., -2.,  1.]]) */
   if (lam>0.0) {
      /* main diag */
      lhsbnd[0+3*nk] += lam;   lhsbnd[(kkmax-1)+3*nk] += lam;
      lhsbnd[1+3*nk] += 5*lam; lhsbnd[(kkmax-2)+3*nk] += 5*lam;
      tmp = 6*lam; for (i=2; i<kkmax-2; ++i) lhsbnd[i+3*nk] += tmp;
      /* 1st subdiag */
      lhsbnd[1+2*nk] += -2*lam; lhsbnd[(kkmax-1)+2*nk] += -2*lam;
      tmp = -4*lam; for (i=2; i<kkmax-1; ++i) lhsbnd[i+2*nk] += tmp;
      /* 2nd subdiag */
      for (i=2; i<kkmax; ++i)   lhsbnd[i+1*nk] += lam;
   }
   for (j=0; j<bw/2; ++j)  /* fill the symmetric part */
      for (i=bw/2-j; i<kkmax; ++i) lhsbnd[(i-bw/2+j)+(bw-1-j)*nk] = lhsbnd[i+j*nk];

   return 1;
}


int bandsol(double *a, double *r, int n, int nd) {
   int i, j, k;
   double aa;
/* Bandsolver with Gauss algorithm */
/* taken from REDUCE */
/*
   bandsol solve a sparse system of linear equations with band-diagonal matrix.
   Band is assumed to be symmetrix relative to the main diagonal. Usage:
   CALL_EXTERNAL('bandsol.so', 'bandsol', a, r, n, nd)
   where a is 2D array [n,m] where n - is the number of equations and nd
           is the width of the band (3 for tri-diagonal system),
           nd is always an odd number. The main diagonal should be in a(*,nd/2)
           The first lower subdiagonal should be in a(1:n-1,nd-2-1), the first
           upper subdiagonal is in a(0:n-2,nd/2+1) etc. For example:
                  / 0 0 X X X \
                  | 0 X X X X |
                  | X X X X X |
                  | X X X X X |
              A = | X X X X X |
                  | X X X X X |
                  | X X X X X |
                  | X X X X 0 |
                  \ X X X 0 0 /
         r is the array of RHS of size n.
*/

   /* Forward sweep */
   for (i=0; i<n-1; ++i) {
      aa = a[i+n*(nd/2)];
      r[i] /= aa;
      for (j=0; j<nd; ++j) a[i+j*n] /= aa;
      for (j=1; j<MIN(nd/2+1,n-i); ++j) {
         aa = a[i+j+n*(nd/2-j)];
         r[i+j] -= r[i] * aa;
         for(k=0; k<n*(nd-j); k+=n) a[i+j+k] -= a[i+k+n*j] * aa;
      }
   }

   /* Backward sweep */
   r[n-1] /= a[n-1+n*(nd/2)];
   for (i=n-1; i>0; --i) {
      for (j=1; j<=MIN(nd/2,i); ++j)
         r[i-j] -= r[i] * a[i-j+n*(nd/2+j)];
      r[i-1] /= a[i-1+n*(nd/2)];
   }
   r[0] /= a[n*(nd/2)];

   return 0;
}


/* wrapper to bandsol argv for FL */
int bandsol_obsolete(double *a, double *b, long *k, long *n){
   void *handle;
   char *error;
   int (*cosine)(int, void*);
   int  rr =1;
   void  *argv[4];
     /*printf("wrapi in %f %ld %ld\n",a[1],k[0],n[0]);*/
   argv[0]=a;
   argv[1]=b;
   argv[2]=k;
   argv[3]=n;

   /* handle = dlopen("/home/raid0/zechmeister/idl/test_load_shared_lib/bla.so", RTLD_LAZY);*/
   handle = dlopen("/home/san0/zechmeister/idl/test_load_shared_lib/bandsol.so.linux.x86_64.64", RTLD_LAZY);
    if (!handle) {
        fprintf(stderr, "%s\n", dlerror());
        exit(EXIT_FAILURE);
    }

   dlerror();    /* Clear any existing error */

   /* Writing: cosine = (double (*)(double)) dlsym(handle, "cos");
       would seem more natural, but the C99 standard leaves
       casting from "void *" to a function pointer undefined.
       The assignment used below is the POSIX.1-2003 (Technical
       Corrigendum 1) workaround; see the Rationale for the
       POSIX specification of dlsym(). */

    *(void **) (&cosine) = dlsym(handle, "bandsol");

/*int (*cosine)(int, void *) = (int (*)(int, void *)) dlysm(handle, "bandsol");
*/
   if ((error = dlerror()) != NULL)  {
        fprintf(stderr, "%s\n", error);
   /*    exit(EXIT_FAILURE);
*/    }

    /*printf("wrapo in : %f \n", a[0]);*/

    rr=  ((* cosine) (1, argv));
/*    printf("out : %d %d %f\n",rr,( int) cosine(1, argv), a[0]);
*/
    dlclose(handle);

   return rr;
}
