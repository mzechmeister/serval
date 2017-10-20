#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
gcc -c  -Wall -O2 -ansi -pedantic -fPIC polyregression.c; gcc -o polyregression.so -shared polyregression.o

   wrapper for poly regression
   converted from original python call
   polynomial.polyfit(x2[ind]-globvar.wcen, y2[ind]/fmod[ind], deg-1, w=fmod[ind]/yerr2[ind], full=True)
*/

int cholsol(double *A, double *b, int n) {
  /* solves Ax = b with cholesky decomposition
     x is returned in b and L in A = LL^T
     returns 1 if Matrix is not positive definite
   */
   int i, j, k;
   double sum;
/*   double A[] = {25,  15,  -5, 15 , 18 ,  0, -5 ,  0 , 11};
   double b[]={25,          15      ,    -5}; 
   double A[] = {4,  2,  14, 2 , 17 ,  -5, 14 ,  -5 , 83};
   double b[]={14,          -101      ,   155};
   n=3;*/
   /* Find the Cholesky factorization of A* A = R* R */
   for (i=0; i<n; ++i) {
      for (j=0; j<i; ++j) {
         sum = A[i+n*j];
         for (k=0; k<j; ++k) sum -= A[i+n*k]*A[j+n*k];
         A[i+n*j] = sum/A[j+n*j];
      }
      sum = A[i+n*i];
      for (k=0; k<i; ++k) sum -= A[i+n*k]*A[i+n*k];

      if (sum >= 0.0) A[i+n*i] = sqrt(sum);
      else return 1;
 /* printf("%f %f %f\n",A[i],A[i+1*n],A[i+2*n]); */
   }
/*   for (i=1; i<n; ++i) A[0:i-1,i] = 0.
   if not keyword_set(b) then return, 1
*/
   /* Solve the lower-triangular system R* y = A* b */
   b[0] = b[0]/A[0];                /* forward substitution */
   /* printf("y%f\n", b[0]); */
   for (i=1; i<n; ++i) {
      sum = 0.;
      for (k=0; k<i; ++k) sum += A[i+n*k]*b[k]; /* b[i] = (b[i]-A[i,0:i-1]#b[0:i-1])/A[i,i] */
      b[i] = (b[i]-sum)/A[i+n*i];
   /* printf("y%f\n", b[i]);*/
   }

   /* Solve the upper-triangular system Rx = y for x.*/
   b[n-1] = b[n-1]/A[n-1+n*(n-1)]       ; /* backward substitution */
   /* printf("e%f\n", b[n-1]);  */
   for (i=n-2; i>=0; --i) {
     sum = 0.;
     for (k=i+1; k<n; ++k) sum += A[n*i+k]*b[k]; /* A[i,i+1:n-1]#b[i+1:n-1] */
     b[i] = (b[i]-sum)/A[n*i+i];
   /*printf("e %f\n", b[i]); */
   }
  return 0;
}

/* double chisq(double *wave, double *flux, double *ferr, double *fmod, _Bool *ind, int n, double wcen, int deg, double *p) {*/
double chisq(double *wave, double *flux, double *ferr, double *fmod, double ind, int n, double wcen, int deg, double *p) {
   double xi, yi, wi, res, poly, sum=0.0; /* solve A*p = b */
   int i, k;
   for (i=0; i<n; ++i) {
/*     if (ind[i]) {*/
     if (fmod[i]>ind) {
         xi = wave[i] - wcen;
         yi = flux[i] / fmod[i];
         wi = fmod[i] / ferr[i];
         poly = p[deg-1];
         for (k=deg-2; k>=0; --k) poly = p[k] + xi*poly ;
         res = wi * (yi-poly);
         sum += res * res;
     }
   }
   /*printf("%f\n",sum);*/
   return sum;
}

/* void polyfit(double *wave, double *flux, double *ferr, double *fmod, _Bool *ind, int n, double wcen, int deg, double *p, double *lhs, double *covar, double chisq) {*/
double polyfit(double *wave, double *flux, double *ferr, double *fmod, double ind, int n, double wcen, int deg, double *p, double *lhs, double *covar) {
   int i,k;
   /* fit polynomial:  flux = fmod * poly(wave)
      minimize:        total( ((flux - fmod*poly)/ferr )^2 )
      reformulate:     total( ((flux/fmod - poly)*fmod/ferr )^2 )
                       total( w*(y - x))^2 )
   */

   double xi, yi, wi, *b=p , *A=covar, mom;   /* solve A*p = b */

   for (k=0; k<2*deg-1; ++k)  A[k] = 0.;
   for (k=0; k<deg; ++k)      b[k] = 0.;

   for (i=0; i<n; ++i) {                               /* create design matrix */
/*      if (ind[i]) { */
     if (fmod[i]>ind) {
         xi = wave[i]-wcen;
         yi = flux[i]/fmod[i];
         wi = fmod[i]/ferr[i]; /* squared errors */
         wi *= wi;

         mom = wi * yi;                   /* fill up RHS */
         b[0] += mom;
         for (k=1; k<deg; ++k) b[k] += (mom *= xi);       /* wi*yi*xi^k */
         mom = wi;
         A[0] += mom;
         for (k=1; k<2*deg-1; ++k) A[k] += (mom *= xi);    /* higher matrix momments wi*xi^k */
/* printf("%f %f %d\n",b[0],flux[i],i); */
      }
   }

   for (k=0; k<deg; ++k)
      for (i=0; i<deg; ++i) lhs[k+deg*i] = A[k+i];

   i = cholsol(lhs, b, deg);
   if (i) return -1.0; /* */
   return chisq(wave, flux, ferr, fmod, ind, n, wcen, deg, p);
}


void interpol1D(double *xn, double *yn, double *x, double *y, int nn, int n) {
   /* xn and x must be sorted */
   double slope; /*, *y=malloc(n * sizeof(double));  */
   int i, k=0;
   for (i=0; i<n; ++i) {
       if (x[i]<=xn[0])
          y[i] = yn[0]; /* no extrapolation */
       else if (x[i]>=xn[nn-1])
          y[i] = yn[nn-1];
       else {
          while (k<nn-1 && xn[k+1]<x[i]) ++k;
          /* 4. Calculate the slope of regions that each x_new value falls in. */
          slope = (yn[k+1]-yn[k]) / (xn[k+1]-xn[k]);
          /* 5. Calculate the actual value for each entry in x_new. */
          y[i] = slope*(x[i]-xn[k]) + yn[k];
  /*printf("%d %f %f %f %f %f\n",i,x[i], y[i],xn[k],yn[k],yn[k+1]);*/
       }
   }
}

/* http://rosettacode.org/wiki/Cholesky_decomposition#C */
