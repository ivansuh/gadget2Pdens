#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

/* define some MIN and MAX funtions */
#define MIN(A,B) ((A)<(B)?(A):(B))
#define MAX(A,B) ((A)>(B)?(A):(B))

/* faster way to calculate power of 2, 3 and 4 */
#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))
#define pow4(x) ((x)*(x)*(x)*(x))
#define pow5(x) ((x)*(x)*(x)*(x)*(x))
#define pow6(x) ((x)*(x)*(x)*(x)*(x)*(x))

float   ran3        (int *idum);
void    indexx      (unsigned long n, double arr[], unsigned long indx[]);
int    *ivector     (long nl, long nh);
void    free_ivector(int *v, long nl, long nh);
