#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "utility.h"



/*===========================================================================
* NUMERICAL RECIPES
*===========================================================================*/
#define NR_END 1
#define FREE_ARG char*
void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
   fprintf(stderr,"Numerical Recipes run-time error...\n");
   fprintf(stderr,"%s\n",error_text);
   fprintf(stderr,"...now exiting to system...\n");
   exit(1);
}
double *vector(long nl, long nh)
/* allocate a vector with subscript range v[nl..nh] */
{
   double *v;
   
   v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
   if (!v) nrerror("allocation failure in vector()");
   return v-nl+NR_END;
}
void free_vector(double *v, long nl, long nh)
/* free a vector allocated with vector() */
{
   free((FREE_ARG) (v+nl-NR_END));
}
void polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
   int i,m,ns=1;
   double den,dif,dift,ho,hp,w;
   double *c,*d;
   
   dif=fabs(x-xa[1]);
   c=vector(1,n);
   d=vector(1,n);
   for (i=1;i<=n;i++) {
      if ( (dift=fabs(x-xa[i])) < dif) {
         ns=i;
         dif=dift;
      }
      c[i]=ya[i];
      d[i]=ya[i];
   }
   *y=ya[ns--];
   for (m=1;m<n;m++) {
      for (i=1;i<=n-m;i++) {
         ho=xa[i]-x;
         hp=xa[i+m]-x;
         w=c[i+1]-d[i];
         if ( (den=ho-hp) == 0.0) nrerror("Error in routine polint");
         den=w/den;
         d[i]=hp*den;
         c[i]=ho*den;
      }
      *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
   }
   free_vector(d,1,n);
   free_vector(c,1,n);
}
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

/* ran3: random generator */
float ran3(int *idum)
{
   static int inext,inextp;
   static long ma[56];
   static int iff=0;
   long mj,mk;
   int i,ii,k;
   
   if (*idum < 0 || iff == 0) {
      iff=1;
      mj=MSEED-(*idum < 0 ? -*idum : *idum);
      mj %= MBIG;
      ma[55]=mj;
      mk=1;
      for (i=1;i<=54;i++) {
         ii=(21*i) % 55;
         ma[ii]=mk;
         mk=mj-mk;
         if (mk < MZ) mk += MBIG;
         mj=ma[ii];
      }
      for (k=1;k<=4;k++)
         for (i=1;i<=55;i++) {
            ma[i] -= ma[1+(i+30) % 55];
            if (ma[i] < MZ) ma[i] += MBIG;
         }
            inext=0;
      inextp=31;
      *idum=1;
   }
   if (++inext == 56) inext=1;
   if (++inextp == 56) inextp=1;
   mj=ma[inext]-ma[inextp];
   if (mj < MZ) mj += MBIG;
   ma[inext]=mj;
   return mj*FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

#define NRANSI
#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 50

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
   int *v;
   
   v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
   if (!v) nrerror("allocation failure in ivector()");
   return v-nl+NR_END;
}
void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
   free((FREE_ARG) (v+nl-NR_END));
}


void indexx(unsigned long n, double arr[], unsigned long indx[])
{
   unsigned long i,indxt,ir=n,itemp,j,k,l=1;
   int jstack=0,*istack;
   double a;
   
   istack=ivector(1,NSTACK);
   for (j=1;j<=n;j++) indx[j]=j;
   for (;;) {
      if (ir-l < M) {
         for (j=l+1;j<=ir;j++) {
            indxt=indx[j];
            a=arr[indxt];
            for (i=j-1;i>=1;i--) {
               if (arr[indx[i]] <= a) break;
               indx[i+1]=indx[i];
            }
            indx[i+1]=indxt;
         }
         if (jstack == 0) break;
         ir=istack[jstack--];
         l=istack[jstack--];
      } else {
         k=(l+ir) >> 1;
         SWAP(indx[k],indx[l+1]);
         if (arr[indx[l+1]] > arr[indx[ir]]) {
            SWAP(indx[l+1],indx[ir])
         }
         if (arr[indx[l]] > arr[indx[ir]]) {
            SWAP(indx[l],indx[ir])
         }
         if (arr[indx[l+1]] > arr[indx[l]]) {
            SWAP(indx[l+1],indx[l])
         }
         i=l+1;
         j=ir;
         indxt=indx[l];
         a=arr[indxt];
         for (;;) {
            do i++; while (arr[indx[i]] < a);
            do j--; while (arr[indx[j]] > a);
            if (j < i) break;
            SWAP(indx[i],indx[j])
         }
         indx[l]=indx[j];
         indx[j]=indxt;
         jstack += 2;
         if (jstack > NSTACK) nrerror("NSTACK too small in indexx.");
         if (ir-i+1 >= j-l) {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
         } else {
            istack[jstack]=j-1;
            istack[jstack-1]=l;
            l=i;
         }
      }
   }
   free_ivector(istack,1,NSTACK);
}
#undef M
#undef NSTACK
#undef SWAP
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software ,2kB. */
