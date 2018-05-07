#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <dlfcn.h>

#define MIN(a,b) ((a<b)?a:b)
#define MAX(a,b) ((a>b)?a:b)

/*
   gcc -c  -Wall -O2 -ansi -pedantic -fPIC cbspline.c; gcc -o cbspline.so -shared cbspline.o
*/

int cholsol(double *A, double *b, int n) {
  /* solves Ax = b with cholesky decomposition
     x is returned in b and L in A = LL^T
     returns 1 if Matrix is not positive definite
   */
   int i, j, k;
   double sum;

   /* Find the Cholesky factorization of A* A = R* R */
   for (i=0; i<n; ++i) {
      for (j=0; j<i; ++j) {
         sum = A[i+n*j];
         for (k=0; k<j; ++k)
            sum -= A[i+n*k] * A[j+n*k];
         A[i+n*j] = sum / A[j+n*j];
      }
      sum = A[i+n*i];
      for (k=0; k<i; ++k)
       sum -= A[i+n*k] * A[i+n*k];

      if (sum >= 0.0) A[i+n*i] = sqrt(sum);
      else return 1;
   }

   /* Solve the lower-triangular system R* y = A* b */
   b[0] = b[0]/A[0];                /* forward substitution */
   for (i=1; i<n; ++i) {
      sum = 0.;
      for (k=0; k<i; ++k)
        sum += A[i+n*k] * b[k];
      b[i] = (b[i]-sum) / A[i+n*i];
   }

   /* Solve the upper-triangular system Rx = y for x. */
   b[n-1] = b[n-1] / A[n-1+n*(n-1)];    /* backward substitution */
   for (i=n-2; i>=0; --i) {
      sum = 0.;
      for (k=i+1; k<n; ++k)
        sum += A[n*i+k] * b[k];
      b[i] = (b[i]-sum) / A[n*i+i];
   }
   return 0;
}

int cholbnd(double *Abnd, double *b, int n, int bw) {
  /* Cholesky solver for band matrices in lower form
   *
   * solves Ax = b with cholesky decomposition
     x is returned in b and L in A = LL^T
     returns 1 if Matrix is not positive definite

   bandsolver index = i,j (i+n*j) =>    i,j (i+ n*(j-i))
python
import numpy as np
bw=3
n=12
A=np.zeros((bw+1,n))*0j +np.nan
for i in range(n):
   for j in range(n):
      if 0<=i-j<=bw: A[i-j,i] = i+j*1j

print A[:,:6]; print A[:,-6:]

A=np.zeros((bw+1,n))*0j +np.nan
for i in range(n):
   for j in range(n):
      if 0<=j-i<=bw: A[j-i,i] = i+j*1j

print A[:,:6]; print A[:,-6:]

   */
   int i, j, k;
   double sum;

   /* Find the Cholesky factorization of A* A = R* R */
   for (i=0; i<n; ++i) {
      for (j=MAX(0,i-bw); j<i; ++j) {
         sum = Abnd[j+n*(i-j)];
         for (k=MAX(0,i-bw); k<j; ++k)
            sum -= Abnd[k+n*(i-k)] * Abnd[k+n*(j-k)];
         Abnd[j+n*(i-j)] = sum / Abnd[j]; /* Li,j = 1/L_jj (A_ij - sum_k(L_ik*L_jk) ) */
      }
      sum = Abnd[i];
      for (k=MAX(0,i-bw); k<i; ++k)
         sum -= Abnd[k+n*(i-k)] * Abnd[k+n*(i-k)];
      if (sum >= 0.0) Abnd[i] = sqrt(sum); /* L_jj = sqrt(A_jj - sum_k L_jk^2) */
      else return 1;
   }

   /* Solve the lower-triangular system R* y = A* b */
   b[0] = b[0] / Abnd[0];                /* forward substitution */
   for (i=1; i<n; ++i) {
      sum = 0.;
      for (k=MAX(0,i-bw); k<i; ++k)
         sum += Abnd[k+n*(i-k)] * b[k];
      b[i] = (b[i]-sum) / Abnd[i];
   }

   /* Solve the upper-triangular system Rx = y for x. */
   /* x_i = (y_i - sum L_ik x_j) / L_ii  */
   b[n-1] = b[n-1] / Abnd[n-1];        /* backward substitution */
   for (i=n-2; i>=0; --i) {
      sum = 0.;
      for (k=i+1; k<MIN(n,i+bw+1); ++k) 
         sum += Abnd[n*(k-i)+i] * b[k];
      b[i] = (b[i]-sum) / Abnd[i];
   }
   return 0;
}

int cholbnd_upper(double *Abnd, double *b, int n, int bw) {
  /* Cholesky solver for band matrices upper form
   *
   * solves Ax = b with cholesky decomposition
     x is returned in b and L in A = LL^T
     returns 1 if Matrix is not positive definite

   bandsolver index = i,j (i+n*j) =>    i,j (i+ n*(bw-j+i))
python
import numpy as np
bw=3
n=12
A=np.zeros((bw+1,n))*0j +np.nan
for i in range(n):
   for j in range(n):
      if 0<=bw+j-i<=bw: A[bw+j-i,i] =   i+j*1j
      if 0<=j-i<=bw: A[j-i,i] = i+j*1j

print A[:,:6]; print A[:,-6:]

   */
   int i, j, k;
   double sum;

   /* Find the Cholesky factorization of A* A = R* R */
   for (i=0; i<n; ++i) {
      for (j=MAX(0,i-bw); j<i; ++j) {
         sum = Abnd[i+n*(bw+j-i)];
         for (k=MAX(0,i-bw); k<j; ++k)
            sum -= Abnd[i+n*(bw+k-i)] * Abnd[j+n*(bw+k-j)];
         Abnd[i+n*(bw+j-i)] = sum / Abnd[j+n*bw]; /* Li,j = 1/L_jj (A_ij - sum_k(L_ik*L_jk) ) */
      }

      sum = Abnd[i+n*bw];
      for (k=MAX(0,i-bw); k<i; ++k)
         sum -= Abnd[i+n*(bw+k-i)] * Abnd[i+n*(bw+k-i)];
      if (sum >= 0.0) Abnd[i+n*bw] = sqrt(sum); /* L_jj = sqrt(A_jj - sum_k L_jk^2) */
      else return 1;
   }

   /* Solve the lower-triangular system R* y = A* b */
   b[0] = b[0] / Abnd[0+n*bw];                /* forward substitution */
   for (i=1; i<n; ++i) {
      sum = 0.;
      for (k=MAX(0,i-bw); k<i; ++k)
         sum += Abnd[i+n*(bw+k-i)] * b[k];
      b[i] = (b[i]-sum) / Abnd[i+n*bw];
   }

   /* Solve the upper-triangular system Rx = y for x. */
   /* x_i = (y_i - sum L_ik x_j) / L_ii  */
   b[n-1] = b[n-1] / Abnd[n-1+n*bw];        /* backward substitution */
   for (i=n-2; i>=0; --i) {
      sum = 0.;
      for (k=i+1; k<MIN(n,i+bw+1); ++k) 
         sum += Abnd[n*(bw+i-k)+k] * b[k];
      b[i] = (b[i]-sum) / Abnd[n*bw+i];
   }
   return 0;
}

void cbspl_Bk(double *x, double *G, long *kk, int n, long fix) {
/*
      kk,pp = divmod(x,1)
      qq = 1-pp
      G[0,:] = qq**3
      G[1,:] = 3*pp**3 - 6*pp**2 + 4
      G[2,:] = 3*qq**3 - 6*qq**2 + 4
      G[3,:] = pp**3
      G /= 6
*/
   long i, i4;
   double pp, qq, tmp;
   for (i=0; i<n; ++i) {
      i4 = 4 * i;
      kk[i] = (long)x[i];    /* integer part */
      if (fix && kk[i] > fix) kk[i] = fix;
      pp = x[i] - kk[i];     /* fractional part */
      qq = 1 - pp;
      tmp = pp*pp; pp *= tmp;
      G[i4+3] = 1./6*pp;
      G[i4+1] = 3./6*pp - tmp + 4./6;
      tmp = qq*qq; qq *= tmp;
      G[i4]   = 1./6*qq;
      G[i4+2] = 3./6*qq - tmp + 4./6;
   }
}


int rhs_fill(double *rhs, double *G, double *w, long *kk, int n) {
   /* fill the right hand side of the equation system
      rhs_k = sum_i w_i B_ki
      for i,kki in enumerate(kk):
         rhs[kk[i]:kk[i]+4] += G[:,i]*wy[i] # rhs[-1] = 0 !!!!
   */
   int j;
   long i, i4;

   for (i=0; i<n; ++i) {
      i4 = i * 4;
      for (j=0; j<4; ++j)
         rhs[kk[i]+j] += w[i] * G[i4+j];
   }
   return 0;
}

int lhsbnd_fill(double *lhsbnd, double *G, double *w, long *kk, int n, int nk) {
  /* fill the left hand side of the equation system
   lhs = w * G * G
   lhsbnd - band matrix ((4 x nk)),  last dim is a dummy
   n - number of data points
   G - matrix with bspline coefficients of data points ((4 x n))
   w - weights ((n))
   kk - reference knots of data points ((n))
   nk - number of knots

   BTB_jk = sum_i B_j(x_i) * B_k(x_i)
   */
   long i, i4;
   int j, k;

   for (i=0; i<n; ++i) {
      i4 = i * 4;
      for (k=0; k<4; ++k)
         for (j=k; j<4; ++j)
            lhsbnd[nk*(j-k)+kk[i]+k] += w[i] * G[i4+k] * G[i4+j];
   }
   return 0;
}


int lhsbnd_fill_upper(double *lhsbnd, double *G, double *w, long *kk, int n, int nk) {
   /* fill the left hand side of the equation system
   lhs = w * G * G
   lhsbnd - band matrix ((7 x nk)),  last dim is a dummy
   n - number of data points
   G - matrix with bspline coefficients of data points ((4 x n))
   w - weights ((n))
   kk - reference knots of data points ((n))
   nk - number of knots

   */
   int j, k;
   long i, i4, kkk;

   for (i=0; i<n; ++i) {
      i4 = i * 4;
      kkk = kk[i] + 3*nk;
      /*printf("%ld %ld %ld\n",i,kk[i],in);*/

      for (j=0; j<4; ++j)
         for (k=0; k<=j; ++k)
            lhsbnd[kkk-nk*(j-k)+j] += w[i] * G[i4+j] * G[i4+k];
      /*printf("%ld %ld %f\n",nk*bw,k,lhsbnd[kkk+(nk-1)*k+j]);*/
   }
   return 0;
}

int add_DTD(double *lhsbnd, int nk, double lam, int pord) {
   /* add penalty and make band matrix symmetric

   lhsbnd - band matrix ((7 x nk)),  last dim is a dummy
   nk - number of knots
   lam - penalty
   pord - penalty order
   */
   int k;

   if (lam > 0.) {
      if (pord==0){
         /* main diag */
         for (k=0; k<nk; ++k) lhsbnd[k] += lam;
      }
      if (pord==1){
         /* main diag */
         lhsbnd[0] += lam;
         for (k=1; k<nk-1; ++k) lhsbnd[k] += 2 * lam;
         lhsbnd[k] += lam;
         /* 1st subdiag */
         for (k=0; k<nk-1; ++k) lhsbnd[k+1*nk] -= lam;
      }
      if (pord==2){
         /* [ 1 -2  1  0  0
             -2  5 -4  1  0
              1 -4  6 -4  1
              0  1 -4  5 -2
              0  0  1 -2  1] */
         /* main diag */
         lhsbnd[0] += lam;
         lhsbnd[1] += 5*lam;
         for (k=2; k<nk-2; ++k) lhsbnd[k] += 6*lam;
         lhsbnd[k] += 5*lam;
         lhsbnd[++k] += lam;
         /* 1st subdiag */
         lhsbnd[1*nk] += -2*lam;
         for (k=1; k<nk-2; ++k) lhsbnd[k+1*nk] += -4*lam;
         lhsbnd[k+1*nk] += -2*lam;
         /* 2nd subdiag */
         for (k=0; k<nk-2; ++k) lhsbnd[k+2*nk] += lam;
      }
   }
   return 0;
}

int add_DTD_upper(double *lhsbnd, int nk, double lam, int pord) {
   /* add penalty and make band matrix symmetric

   lhsbnd - band matrix ((7 x nk)),  last dim is a dummy
   nk - number of knots
   lam - penalty
   pord - penalty order
   */
   int j, k, bw=3;

   if (lam > 0.) {
      if (pord==0){
         /* main diag */
         for (k=0; k<nk; ++k) lhsbnd[k+3*nk] += lam;
      }
      if (pord==1){
         /* main diag */
         lhsbnd[0+3*nk] += lam;
         for (k=1; k<nk-1; ++k) lhsbnd[k+3*nk] += 2 * lam;
         lhsbnd[k+3*nk] += lam;
         /* 1st subdiag */
         for (k=1; k<nk; ++k) lhsbnd[k+2*nk] -= lam;
      }
      if (pord==2){
         /* [ 1 -2  1  0  0
             -2  5 -4  1  0
              1 -4  6 -4  1
              0  1 -4  5 -2
              0  0  1 -2  1] */
         /* main diag */
         lhsbnd[0+3*nk] += lam;
         lhsbnd[1+3*nk] += 5*lam;
         for (k=2; k<nk-2; ++k) lhsbnd[k+3*nk] += 6*lam;
         lhsbnd[k+3*nk] += 5*lam;
         lhsbnd[++k +3*nk] += lam;
         /* 1st subdiag */
         lhsbnd[1+2*nk] += -2*lam;
         for (k=2; k<nk-1; ++k) lhsbnd[k+2*nk] += -4*lam;
         lhsbnd[k+2*nk] += -2*lam;
         /* 2nd subdiag */
         for (k=2; k<nk; ++k) lhsbnd[k+1*nk] += lam;
      }
   }

   /* fill the symmetric part */
   for (j=0; j<bw; ++j)
      for (k=bw-j; k<nk; ++k) lhsbnd[(k-bw+j)+(2*bw-j)*nk] = lhsbnd[k+j*nk];

   return 0;
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
