#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utility.h"

#define GADGET_LUNIT 1.0e-3
#define GADGET_MUNIT 1.0e10
#define GADGET           2
#define BYTESWAP         0
#define TRUE    1
#define FALSE   0

#define GADGET_SKIP    ReadUInt(fd,&blklen,SWAPBYTES);
#define NFILL   (256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8)
#define IDX(ix,iy,iz)  ((iz)*NGRID*NGRID + (iy)*NGRID + (ix))

int ReadInt       (FILE *fptr,int *n,int swap);
int ReadFloat     (FILE *fptr,float *n, int swap);
int ReadDouble    (FILE *fptr,double *n,int swap);
int ReadLong      (FILE *fptr,long *n,int swap);
int ReadUInt      (FILE *fptr,unsigned int *n,int swap);
int ReadLongLong  (FILE *fptr,long long *n,int swap);


int allocate_memory(void);

int load_snapshot(char *fname, int files);

int unit_conversion(void);

int write_ascii(char *fname, double *rdens);

int reordering();

void get_density(int NGRID);

void assign(int NGRID, float *dens);

void interp(int NGRID, float *dens);

struct io_header_1
{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[NFILL];  /* fills to 256 Bytes */
} header1;

unsigned int blklen;
int     NumPart, Ngas;
double	TotMass;

struct particle_data 
{
  float  	Pos[3];
  float		Vel[3];
  float  	Mass;
  float		Rho;
  int		ID;

  int    	Type;

  float  	U, Temp, Ne;
} *P;


double  	Time, Redshift;



/* Here we load a snapshot file. It can be distributed
 * onto several files (for files>1).
 * The particles are brought back into the order
 * implied by their ID's.
 * A unit conversion routine is called to do unit
 * conversion, and to evaluate the gas temperature.
 */
int main(argc,argv)
int argc;
char **argv;
{
  char		indata[200],outdata[200];
  int  		files;
  int  		NGRID, ilen, i, j, idummy;
  double        *rdens; 
  float		*dens;
  long unsigned *index;
  long unsigned Npart, ipart;
  double        x_fac;
  float         MaxDens, MinDens, CurDens;

  if(argc<4)
     {
      fprintf(stderr,"usage: %s Number_of_files Gadget_file NGRID\n",*argv);
      exit(1);
     }
     
  printf("=====================================================================\n");
  printf("  Read Gadget binary and calculate density at each particle position \n");
  printf("               (write particles to the ascii file)                   \n");
  printf("=====================================================================\n");   

  MaxDens = -1000.;
  MinDens =  1000.;
  NGRID = (int) atoi(argv[3]);

  files = (int) atoi(argv[1]);                               /* number of files per snapshot */
  
  ilen = strlen(argv[2]);
   for(i=ilen-1; i>=0; i--)
      /* chop filename at last '/' */
      if(argv[2][i] == '/') break;
   for(j=i+1; j<ilen; j++)
      indata[j-(i+1)] = argv[2][j];
   indata[j-(i+1)] = '\0';
   sprintf(outdata,"Pdens-%d_%s",NGRID,indata);
   fprintf(stderr,"\n-> writing result to file %s\n\n",outdata);
  
  
/* Reads snapshot (already converts (v, m) to (km/sec, Msun/h) */
  load_snapshot(indata, files);
  
  reordering();  /* call this routine only if your ID's are set properly */

  /* get density at particle positions */
  get_density(NGRID);
  
  /* scale coordinates to physical coordinates */
for(ipart = 1; ipart <= NumPart; ipart++)
    {
     P[ipart].Pos[0]  *= (header1.BoxSize*GADGET_LUNIT);
     P[ipart].Pos[1]  *= (header1.BoxSize*GADGET_LUNIT);
     P[ipart].Pos[2]  *= (header1.BoxSize*GADGET_LUNIT);

     if(P[ipart].Rho > MaxDens) MaxDens = P[ipart].Rho;
     if(P[ipart].Rho < MinDens) MinDens = P[ipart].Rho;
    }
  /* write DATA to outfile */
  fprintf(stderr,"  o writing output (minDens=%g maxDens=%g)...\n",MinDens,MaxDens);
  write_ascii(outdata,rdens);
}



/* here the particle data is at your disposal 
 */
int write_ascii(char *fname,double *rdens)
{
   FILE *fd;             /* File handle */
   char   buf[200];

   int i,j,k;            /* Dummy variables */


	sprintf(buf,"%s",fname);
        if(!(fd=fopen(buf,"w")))
	{
	  printf("can't open file `%s`\n",buf);
	  exit(0);
	}
        printf("Writing ascii data to file `%s`\n",buf);
#ifdef HEADER
        for(i=0;i<6;i++) 
	        fprintf(fd,"# np[%d]=%9d\t nTotal[%d]=%9d\t mass[%d]=%g\n",i,header1.npart[i],i,header1.npartTotal[i],i,header1.mass[i]*GADGET_MUNIT);

      	  fprintf(fd,"\n# Expansion factor: %lf\t Redshift: %lf\n",header1.time,header1.redshift);
      	  fprintf(fd,"# BoxSize (Mpc/h): %lf\n",header1.BoxSize*GADGET_LUNIT);
      	  fprintf(fd,"# Omega0: %lf\t OmegaLambda: %lf \t HubbleParam: %lf\n\n",header1.Omega0,header1.OmegaLambda,header1.HubbleParam);
#endif

#ifdef VELOCITY
	for(i=1;i<=NumPart;i++)
	  {
        fprintf(fd,"%g  %g  %g  %lf  %lf  %lf  %g %d\n",P[i].Pos[0],P[i].Pos[1],P[i].Pos[2],P[i].Vel[0],P[i].Vel[1],P[i].Vel[2],P[i].Rho,P[i].ID-1);
      }
        printf("done.\n");
#else
for(i=1;i<=NumPart;i++)
	  {
        fprintf(fd,"%lf  %lf  %lf  %lf %d\n",P[i].Pos[0],P[i].Pos[1],P[i].Pos[2],P[i].Rho,P[i].ID-1);
      }
        printf("done.\n");
#endif
return 0;
}


int unit_conversion(void)
{
  float  x_min=1E40, x_max=-1E40, y_min=1E40, y_max=-1E40, z_min=1E40, z_max=-1E40;

  int i;

/* Convert Gadget unit to physical units */
  printf("Converting Gadget Units to physical units.\n");
  header1.BoxSize *=GADGET_LUNIT;
  
  for(i=1;i<=NumPart;i++)
    {
     P[i].Pos[0] *=GADGET_LUNIT;  
     P[i].Pos[1] *=GADGET_LUNIT;  
     P[i].Pos[2] *=GADGET_LUNIT;  
     
     if(P[i].Pos[0] > x_max) x_max=P[i].Pos[0];
        if(P[i].Pos[1] > y_max) y_max=P[i].Pos[1];
        if(P[i].Pos[2] > z_max) z_max=P[i].Pos[2];
        if(P[i].Pos[0] < x_min) x_min=P[i].Pos[0];
        if(P[i].Pos[1] < y_min) y_min=P[i].Pos[1];
        if(P[i].Pos[2] < z_min) z_min=P[i].Pos[2];
  
  /* convert to peculiar velocities in km/sec */
     P[i].Vel[0] *= sqrt(header1.time);
     P[i].Vel[1] *= sqrt(header1.time);
     P[i].Vel[2] *= sqrt(header1.time);
  /* convert to Msun */
     P[i].Mass *= GADGET_MUNIT; 
     
    } 
    printf("ranges:\n");
    printf("x = %g  %g\n",x_min,x_max);
    printf("y = %g  %g\n",y_min,y_max);
    printf("z = %g  %g\n",z_min,z_max);
    
    printf("done.\n");
    return 0;
}

/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.
 */
int load_snapshot(char *fname, int files)
{
  FILE *fd;
  char   buf[200],DATA[512];
  int    i,j,k,l,m,dummy,ntot_withmasses;
  int    t,n,off,pc_new,pc_sph,pc,ic,c,SWAPBYTES;
  float  x_min=1E40, x_max=-1E40, y_min=1E40, y_max=-1E40, z_min=1E40, z_max=-1E40;

//#define SKIP fread(&dummy, sizeof(dummy), 1, fd);

  for(i=0, pc=1; i<files; i++, pc=pc_new)
    {
      if(files>1)
	sprintf(buf,"%s.%d",fname,i);
      else
	sprintf(buf,"%s",fname);

      if(!(fd=fopen(buf,"r")))
	{
	  printf("can't open file `%s`\n",buf);
	  exit(0);
	}
#if (BYTESWAP==1)
      SWAPBYTES = TRUE;
      printf("o reading BYTESWAPed `%s' ...\n",buf); fflush(stdout);
#else
      SWAPBYTES = FALSE;
      printf("o reading `%s' ...\n",buf); fflush(stdout);
#endif
      
      /*================= read in GADGET IO header =================*/
#if (GADGET==2)
      GADGET_SKIP;
      fread(DATA,sizeof(char),4,fd);
      DATA[4] = '\0';
      fprintf(stderr,"\n o reading %s\n",DATA);
      GADGET_SKIP;
      GADGET_SKIP;
#endif
      
      GADGET_SKIP;
      ReadInt(fd, &(header1.npart[0]), SWAPBYTES);
      ReadInt(fd, &(header1.npart[1]), SWAPBYTES);
      ReadInt(fd, &(header1.npart[2]), SWAPBYTES);
      ReadInt(fd, &(header1.npart[3]), SWAPBYTES);
      ReadInt(fd, &(header1.npart[4]), SWAPBYTES);
      ReadInt(fd, &(header1.npart[5]), SWAPBYTES);
      ReadDouble(fd, &(header1.mass[0]), SWAPBYTES);
      ReadDouble(fd, &(header1.mass[1]), SWAPBYTES);
      ReadDouble(fd, &(header1.mass[2]), SWAPBYTES);
      ReadDouble(fd, &(header1.mass[3]), SWAPBYTES);
      ReadDouble(fd, &(header1.mass[4]), SWAPBYTES);
      ReadDouble(fd, &(header1.mass[5]), SWAPBYTES);
      ReadDouble(fd, &(header1.time), SWAPBYTES);
      ReadDouble(fd, &(header1.redshift), SWAPBYTES);
      ReadInt(fd, &(header1.flag_sfr), SWAPBYTES);
      ReadInt(fd, &(header1.flag_feedback), SWAPBYTES);
      ReadInt(fd, &(header1.npartTotal[0]), SWAPBYTES);
      ReadInt(fd, &(header1.npartTotal[1]), SWAPBYTES);
      ReadInt(fd, &(header1.npartTotal[2]), SWAPBYTES);
      ReadInt(fd, &(header1.npartTotal[3]), SWAPBYTES);
      ReadInt(fd, &(header1.npartTotal[4]), SWAPBYTES);
      ReadInt(fd, &(header1.npartTotal[5]), SWAPBYTES);
      ReadInt(fd, &(header1.flag_cooling), SWAPBYTES);
      ReadInt(fd, &(header1.num_files), SWAPBYTES);
      ReadDouble(fd, &(header1.BoxSize), SWAPBYTES);
      ReadDouble(fd, &(header1.Omega0), SWAPBYTES);
      ReadDouble(fd, &(header1.OmegaLambda), SWAPBYTES);
      ReadDouble(fd, &(header1.HubbleParam), SWAPBYTES);
      for(ic=0; ic<NFILL; ic++)
        {
         c = fgetc(fd);
         if (c == EOF)
            return(FALSE);
         header1.fill[ic]   = c;
         header1.fill[ic+1] = '\0';
        }
      GADGET_SKIP;
     
      
      printf("   Header:\n");
      for(l=0;l<6;l++) 
	printf("   np[%d]=%9d\t nTotal[%d]=%9d\t mass[%d]=%g\n",i,header1.npart[l],i,header1.npartTotal[l],i,header1.mass[l]*GADGET_MUNIT);

      printf("   \nExpansion factor: %lf\t Redshift: %lf\n",header1.time,header1.redshift);
      printf("   flag_sfr		: %d\n", header1.flag_sfr);
      printf("   flag_feedback	: %d\n", header1.flag_feedback);
      printf("   flag_cooling	: %d\n", header1.flag_cooling);
      printf("   num_files		: %d\n", header1.num_files);
      printf("   BoxSize (Mpc/h)	: %lf\n",header1.BoxSize*GADGET_LUNIT);
      printf("   Omega0: %lf\t OmegaLambda: %lf \t HubbleParam: %lf\n",header1.Omega0,header1.OmegaLambda,header1.HubbleParam);

//      exit(0); /* For testing purposes */
      if(files==1)
	{
	  for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
	    NumPart+= header1.npart[k];
	  Ngas= header1.npart[0];
	}
      else
	{
	  for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
	    NumPart+= header1.npartTotal[k];
	  Ngas= header1.npartTotal[0];
	}

      for(k=0, ntot_withmasses=0; k<5; k++)
	{
	  if(header1.mass[k]==0)
	    ntot_withmasses+= header1.npart[k];
	}

      if(i==0)
	allocate_memory();

      /*================= read in GADGET particles =================*/
#if (GADGET==2)
      GADGET_SKIP;
      fread(DATA,sizeof(char),4,fd);
      DATA[4] = '\0';
      fprintf(stderr,"\n o reading %s\n",DATA);
      GADGET_SKIP;
      GADGET_SKIP;
#endif
      GADGET_SKIP; 
      for(k=0,pc_new=pc;k<6;k++)
	  {
	     for(n=0;n<header1.npart[k];n++)
	       {
	        ReadFloat(fd, &(P[pc_new].Pos[0]), SWAPBYTES);
            ReadFloat(fd, &(P[pc_new].Pos[1]), SWAPBYTES);
            ReadFloat(fd, &(P[pc_new].Pos[2]), SWAPBYTES);
            
	    /* convert to [0,1] */
            P[pc_new].Pos[0] /= header1.BoxSize;
            P[pc_new].Pos[1] /= header1.BoxSize;
            P[pc_new].Pos[2] /= header1.BoxSize;
            
            if(P[pc_new].Pos[0] > x_max) x_max=P[pc_new].Pos[0];
            if(P[pc_new].Pos[1] > y_max) y_max=P[pc_new].Pos[1];
            if(P[pc_new].Pos[2] > z_max) z_max=P[pc_new].Pos[2];
            if(P[pc_new].Pos[0] < x_min) x_min=P[pc_new].Pos[0];
            if(P[pc_new].Pos[1] < y_min) y_min=P[pc_new].Pos[1];
            if(P[pc_new].Pos[2] < z_min) z_min=P[pc_new].Pos[2];
            
            
            if(n==0)
               printf("   P[%d].Pos = %f  %f  %f\n", pc_new, P[pc_new].Pos[0], P[pc_new].Pos[1], P[pc_new].Pos[2]);
	    pc_new++;
	    }
	}
      GADGET_SKIP;
      printf("   ranges:\n");
      printf("   x = %g  %g\n",x_min,x_max);
      printf("   y = %g  %g\n",y_min,y_max);
      printf("   z = %g  %g\n",z_min,z_max);

      /*================= read in GADGET velocities =================*/
#if (GADGET==2)
      GADGET_SKIP;
      fread(DATA,sizeof(char),4,fd);
      DATA[4] = '\0';
      fprintf(stderr,"\n o reading %s\n",DATA);
      GADGET_SKIP;
      GADGET_SKIP;
#endif
      GADGET_SKIP;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	     ReadFloat(fd, &(P[pc_new].Vel[0]), SWAPBYTES);
         ReadFloat(fd, &(P[pc_new].Vel[1]), SWAPBYTES);
         ReadFloat(fd, &(P[pc_new].Vel[2]), SWAPBYTES);
             
	     /* convert to peculiar velocities in km/sec */
             P[pc_new].Vel[0] *= sqrt(header1.time);
             P[pc_new].Vel[1] *= sqrt(header1.time);
             P[pc_new].Vel[2] *= sqrt(header1.time);
            
             if(n==0)
               printf("   P[%d].Vel = %f  %f  %f\n", pc_new, P[pc_new].Vel[0], P[pc_new].Vel[1], P[pc_new].Vel[2]); 
	     pc_new++;
	    }
	}
      GADGET_SKIP;
    
      /*================= read in GADGET id's =================*/
#if (GADGET==2)
      GADGET_SKIP;
      fread(DATA,sizeof(char),4,fd);
      DATA[4] = '\0';
      fprintf(stderr,"\n o reading %s\n",DATA);
      GADGET_SKIP;
      GADGET_SKIP;
#endif
      GADGET_SKIP;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      ReadInt(fd, &(P[pc_new].ID), SWAPBYTES);
	      pc_new++;
	    }
	}
      GADGET_SKIP;

      /*================ read in GADGET masses ================*/
      if(ntot_withmasses>0)
#if (GADGET==2)
         GADGET_SKIP;
         fread(DATA,sizeof(char),4,fd);
         DATA[4] = '\0';
         fprintf(stderr,"\n o reading %s\n",DATA);
         GADGET_SKIP;
         GADGET_SKIP;
#endif
	GADGET_SKIP;
      for(k=0, pc_new=pc; k<6; k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      P[pc_new].Type=k;

	      if(header1.mass[k]==0)
		    ReadFloat(fd, &(P[pc_new].Mass), SWAPBYTES);
	      else
		P[pc_new].Mass= header1.mass[k];

	      TotMass += P[pc_new].Mass;

	      pc_new++;
	    }
	}
      if(ntot_withmasses>0)
	GADGET_SKIP;
      

      if(header1.npart[0]>0)
	{
#if (GADGET==2)
      GADGET_SKIP;
      fread(DATA,sizeof(char),4,fd);
      DATA[4] = '\0';
      fprintf(stderr,"\n o reading %s\n",DATA);
      GADGET_SKIP;
      GADGET_SKIP;
#endif
	  GADGET_SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      ReadFloat(fd, &(P[pc_sph].U), SWAPBYTES);
	      pc_sph++;
	    }
	  GADGET_SKIP;

#if (GADGET==2)
         GADGET_SKIP;
         fread(DATA,sizeof(char),4,fd);
         DATA[4] = '\0';
         fprintf(stderr,"\n o reading %s\n",DATA);
         GADGET_SKIP;
         GADGET_SKIP;
#endif
	  GADGET_SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      ReadFloat(fd, &(P[pc_sph].Rho), SWAPBYTES);
	      pc_sph++;
	    }
	  GADGET_SKIP;

	  if(header1.flag_cooling)
	    {
#if (GADGET==2)
         GADGET_SKIP;
         fread(DATA,sizeof(char),4,fd);
         DATA[4] = '\0';
         fprintf(stderr,"\n o reading %s\n",DATA);
         GADGET_SKIP;
         GADGET_SKIP;
#endif
	      GADGET_SKIP;
	      for(n=0, pc_sph=pc; n<header1.npart[0];n++)
		{
		  ReadFloat(fd, &(P[pc_sph].Ne), SWAPBYTES);
		  pc_sph++;
		}
	      GADGET_SKIP;
	    }
	  else
	    for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	      {
		P[pc_sph].Ne= 1.0;
		pc_sph++;
	      }
	}

      fclose(fd);
    }


  Time= header1.time;
  Redshift= header1.time;
  return 0;
}


/* this routine allocates the memory for the 
 * particle data.
 */
int allocate_memory(void)
{
  printf("allocating memory...\n");

  if(!(P=malloc(NumPart*sizeof(struct particle_data))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
  
  P--;   /* start with offset 1 */

  printf("allocating memory...done\n");
return 0;
}




/* This routine brings the particles back into
 * the order of their ID's.
 * NOTE: The routine only works if the ID's cover
 * the range from 1 to NumPart !
 * In other cases, one has to use more general
 * sorting routines.
 */
int reordering(void)
{
  int i,j;
  int idsource, idsave, dest;
  struct particle_data psave, psource;


/*  for(j=1; j<=NumPart; j++) 
      P[j].ID+=1;
*/
  printf("reordering....\n");

  for(i=1; i<=NumPart; i++)
    {
      if(P[i].ID != i)
	{
	  psource= P[i];
	  idsource=P[i].ID;
	  dest=P[i].ID;

	  do
	    {
	      idsave=P[dest].ID;
	      psave= P[dest];

	      P[dest]= psource;
	      P[dest].ID= idsource;
	      
	      if(dest == i) {
		break;
		}
	      psource= psave;
	      idsource=idsave;

	      dest=idsource;
	    }
	  while(1);
	}
    }

  printf("done.\n");
return 0;
}

/*=============================================================================
* get_density()    calculate the density at each particle position
*                  using a grid with NGRID^3 grid points
*
*                  particle positions are assumed to lie in
*                  the range [0,1] !
*
* 
*     INPUT:    x[npart], y[npart], z[npart], L1DIM
*
*     OUTPUT:   P[npart].Rho
*
*============================================================================*/
void get_density(int NGRID)
{
   double rho_mean;
   float *dens;
   int   ix, iy, iz;
   
   rho_mean = (double)NumPart / (double) (NGRID*NGRID*NGRID);
   
   dens   = (float *) calloc(NGRID*NGRID*NGRID, sizeof(float*));

   /* assign particles to the grid */
   fprintf(stderr,"  o assigning particles to grid (mass density)...");
   assign(NGRID, dens);
   fprintf(stderr,"done\n");
   
   /* get density in terms of mean density */
   for(ix=0; ix<NGRID; ix++)
      for(iy=0; iy<NGRID; iy++)
         for(iz=0; iz<NGRID; iz++)
            dens[IDX(ix,iy,iz)] /= rho_mean;
   
   /* interpolate density contrast back to particle positions */ 
   fprintf(stderr,"  o interpolating to particle positions...");
   interp(NGRID, dens);
   fprintf(stderr,"done\n");
   free(dens);
}

/*==================================================================
* assign:   assign particles to the grid points using a TSC scheme
*==================================================================*/
void assign(int NGRID, float *dens)
{
   int           ix, iy, iz, ixp1, iyp1, izp1, ixm1, iym1, izm1;
   long unsigned ipart;
   float         rrx, rry, rrz, www;
   float         hx, hy, hz;
   float         hx0, hy0, hz0, hxp1, hyp1, hzp1, hxm1, hym1, hzm1;
   
   for(ix=0; ix<NGRID; ix++)
      for(iy=0; iy<NGRID; iy++)
         for(iz=0; iz<NGRID; iz++)
            dens[IDX(ix,iy,iz)]=0.0;
   
   for(ipart=0; ipart<NumPart; ipart++)
     {
      /* coordinates in grid units */
      rrx = P[ipart+1].Pos[0] * (float)(NGRID);
      rry = P[ipart+1].Pos[1] * (float)(NGRID);
      rrz = P[ipart+1].Pos[2] * (float)(NGRID);
      
      /* particle weight */

      www = 1.0;
      
      /* index of nearest grid point */
      ix  = (int)(rrx+0.5);
      iy  = (int)(rry+0.5);
      iz  = (int)(rrz+0.5);
      
      /* distance to nearest grid point */
      hx  = rrx - (float)ix;
      hy  = rry - (float)iy;
      hz  = rrz - (float)iz;
      
      /* keep track of peridoc boundaries */
      ix=(int)fmod(ix,NGRID);
      iy=(int)fmod(iy,NGRID);
      iz=(int)fmod(iz,NGRID);
      
      /* calculate TSC weights */
      hx0=0.75 - hx*hx;
      hxp1=0.5* pow2(0.5 + hx);
      hxm1=0.5* pow2(0.5 - hx);
      hy0=0.75 - hy*hy;
      hyp1=0.5* pow2(0.5 + hy);
      hym1=0.5* pow2(0.5 - hy);
      hz0= 0.75 - hz*hz;
      hzp1=0.5* pow2(0.5 + hz);
      hzm1=0.5* pow2(0.5 - hz);
      
      /* keep track of peridoc boundaries */
      ixp1=(int)fmod(ix+1,NGRID);
      iyp1=(int)fmod(iy+1,NGRID);
      izp1=(int)fmod(iz+1,NGRID);
      ixm1=(int)fmod(ix-1+NGRID,NGRID);
      iym1=(int)fmod(iy-1+NGRID,NGRID);
      izm1=(int)fmod(iz-1+NGRID,NGRID);
      
      /* assign particle according to weights to 27 neighboring nodes */
      dens[IDX(ixm1,iym1,izm1)] += hxm1*hym1 *hzm1 * www;
      dens[IDX(ix,  iym1,izm1)] += hx0 *hym1 *hzm1 * www;
      dens[IDX(ixp1,iym1,izm1)] += hxp1*hym1 *hzm1 * www;
      dens[IDX(ixm1,  iy,izm1)] += hxm1*hy0  *hzm1 * www;
      dens[IDX(  ix,  iy,izm1)] += hx0 *hy0  *hzm1 * www;
      dens[IDX(ixp1,  iy,izm1)] += hxp1*hy0  *hzm1 * www;
      dens[IDX(ixm1,iyp1,izm1)] += hxm1*hyp1 *hzm1 * www;
      dens[IDX(  ix,iyp1,izm1)] += hx0 *hyp1 *hzm1 * www;
      dens[IDX(ixp1,iyp1,izm1)] += hxp1*hyp1 *hzm1 * www;
      dens[IDX(ixm1,iym1,  iz)] += hxm1*hym1 *hz0 * www;
      dens[IDX(  ix,iym1,  iz)] += hx0 *hym1 *hz0 * www;
      dens[IDX(ixp1,iym1,  iz)] += hxp1*hym1 *hz0 * www;
      dens[IDX(ixm1,  iy,  iz)] += hxm1*hy0  *hz0 * www;
      dens[IDX(  ix,  iy,  iz)] += hx0 *hy0  *hz0 * www;
      dens[IDX(ixp1,  iy,  iz)] += hxp1*hy0  *hz0 * www;
      dens[IDX(ixm1,iyp1,  iz)] += hxm1*hyp1 *hz0 * www;
      dens[IDX(  ix,iyp1,  iz)] += hx0 *hyp1 *hz0 * www;
      dens[IDX(ixp1,iyp1,  iz)] += hxp1*hyp1 *hz0 * www;
      dens[IDX(ixm1,iym1,izp1)] += hxm1*hym1 *hzp1 * www;
      dens[IDX(  ix,iym1,izp1)] += hx0 *hym1 *hzp1 * www;
      dens[IDX(ixp1,iym1,izp1)] += hxp1*hym1 *hzp1 * www;
      dens[IDX(ixm1,  iy,izp1)] += hxm1*hy0  *hzp1 * www;
      dens[IDX(  ix,  iy,izp1)] += hx0 *hy0  *hzp1 * www;
      dens[IDX(ixp1,  iy,izp1)] += hxp1*hy0  *hzp1 * www;
      dens[IDX(ixm1,iyp1,izp1)] += hxm1*hyp1 *hzp1 * www;
      dens[IDX(  ix,iyp1,izp1)] += hx0 *hyp1 *hzp1 * www;
      dens[IDX(ixp1,iyp1,izp1)] += hxp1*hyp1 *hzp1 * www;
     }
}

/*==============================================================
* interp:  interpolate the density the grid points back to the 
*          particle positions using again a TSC scheme
*==============================================================*/
void interp(int NGRID, float *dens)
{
   int           ix, iy, iz, ixp1, iyp1, izp1, ixm1, iym1, izm1;
   long unsigned ipart;
   float         rrx, rry, rrz;
   float         hx, hy, hz;
   float         hx0, hy0, hz0, hxp1, hyp1, hzp1, hxm1, hym1, hzm1;
   float         ac;
   
   for(ipart=0; ipart<NumPart; ipart++)
     {
      /* coordinates in grid units */
      rrx = P[ipart+1].Pos[0] * (float)(NGRID);
      rry = P[ipart+1].Pos[1] * (float)(NGRID);
      rrz = P[ipart+1].Pos[2] * (float)(NGRID);
      
      /* index of nearest grid point */
      ix  = (int)(rrx+0.5);
      iy  = (int)(rry+0.5);
      iz  = (int)(rrz+0.5);
      
      /* distance to nearest grid point */
      hx  = rrx - (float)ix;
      hy  = rry - (float)iy;
      hz  = rrz - (float)iz;
      
      /* keep track of peridoc boundaries */
      ix=(int)fmod(ix,NGRID);
      iy=(int)fmod(iy,NGRID);
      iz=(int)fmod(iz,NGRID);
      
      /* calculate TSC weights */
      hx0=0.75 - hx*hx;
      hxp1=0.5* pow2(0.5 + hx);
      hxm1=0.5* pow2(0.5 - hx);
      hy0=0.75 - hy*hy;
      hyp1=0.5* pow2(0.5 + hy);
      hym1=0.5* pow2(0.5 - hy);
      hz0= 0.75 - hz*hz;
      hzp1=0.5* pow2(0.5 + hz);
      hzm1=0.5* pow2(0.5 - hz);
      
      /* keep track of peridoc boundaries */
      ixp1=(int)fmod(ix+1,NGRID);
      iyp1=(int)fmod(iy+1,NGRID);
      izp1=(int)fmod(iz+1,NGRID);
      ixm1=(int)fmod(ix-1+NGRID,NGRID);
      iym1=(int)fmod(iy-1+NGRID,NGRID);
      izm1=(int)fmod(iz-1+NGRID,NGRID);
      
      
      ac =dens[IDX(ixm1,iym1,izm1)]* hxm1*hym1 *hzm1
         +dens[IDX(ix  ,iym1,izm1)]* hx0 *hym1 *hzm1
         +dens[IDX(ixp1,iym1,izm1)]* hxp1*hym1 *hzm1
         +dens[IDX(ixm1,iy  ,izm1)]* hxm1*hy0  *hzm1
         +dens[IDX(ix  ,iy  ,izm1)]* hx0 *hy0  *hzm1
         +dens[IDX(ixp1,iy  ,izm1)]* hxp1*hy0  *hzm1
         +dens[IDX(ixm1,iyp1,izm1)]* hxm1*hyp1 *hzm1
         +dens[IDX(ix  ,iyp1,izm1)]* hx0 *hyp1 *hzm1
         +dens[IDX(ixp1,iyp1,izm1)]* hxp1*hyp1 *hzm1
         +dens[IDX(ixm1,iym1,iz  )]* hxm1*hym1 *hz0
         +dens[IDX(ix  ,iym1,iz  )]* hx0 *hym1 *hz0
         +dens[IDX(ixp1,iym1,iz  )]* hxp1*hym1 *hz0
         +dens[IDX(ixm1,iy  ,iz  )]* hxm1*hy0  *hz0
         +dens[IDX(ix  ,iy  ,iz  )]* hx0 *hy0  *hz0;
      ac+=dens[IDX(ixp1,iy  ,iz  )]* hxp1*hy0  *hz0
         +dens[IDX(ixm1,iyp1,iz  )]* hxm1*hyp1 *hz0
         +dens[IDX(ix  ,iyp1,iz  )]* hx0 *hyp1 *hz0
         +dens[IDX(ixp1,iyp1,iz  )]* hxp1*hyp1 *hz0
         +dens[IDX(ixm1,iym1,izp1)]* hxm1*hym1 *hzp1
         +dens[IDX(ix  ,iym1,izp1)]* hx0 *hym1 *hzp1
         +dens[IDX(ixp1,iym1,izp1)]* hxp1*hym1 *hzp1
         +dens[IDX(ixm1,iy  ,izp1)]* hxm1*hy0  *hzp1
         +dens[IDX(ix  ,iy  ,izp1)]* hx0 *hy0  *hzp1
         +dens[IDX(ixp1,iy  ,izp1)]* hxp1*hy0  *hzp1
         +dens[IDX(ixm1,iyp1,izp1)]* hxm1*hyp1 *hzp1
         +dens[IDX(ix  ,iyp1,izp1)]* hx0 *hyp1 *hzp1
         +dens[IDX(ixp1,iyp1,izp1)]* hxp1*hyp1 *hzp1;
      P[ipart+1].Rho = (double) ac;
      
     }
}


/****************************************************************************/
/************    UTILITY ROUTINES TO READ BYTESWAPPED DATA    ***************/
/****************************************************************************/


/*
 Read a string of n characters
 */
int ReadString(FILE *fptr,char *s,int n)
{
   int i,c;
   
   if(sizeof(char) != 1)
     {
      fprintf(stderr,"ReadString: sizeof(char)=%ld and not 1\n",sizeof(char));
      exit(0);
     }
   
   s[0] = '\0';
   for (i=0;i<n;i++) {
      c = fgetc(fptr);
      if (c == EOF)
         return(FALSE);
      s[i] = c;
      s[i+1] = '\0';
   }
   return(TRUE);
}

/*
 Read a possibly byte swapped integer
 */
int ReadInt(FILE *fptr,int *n,int swap)
{
   unsigned char *cptr,tmp;
   
   if(sizeof(int) != 4)
     {
      fprintf(stderr,"ReadInt: sizeof(int)=%ld and not 4\n",sizeof(int));
      exit(0);
     }
   
   if (fread(n,4,1,fptr) != 1)
      return(FALSE);
   if (swap) {
      cptr = (unsigned char *)n;
      tmp     = cptr[0];
      cptr[0] = cptr[3];
      cptr[3] = tmp;
      tmp     = cptr[1];
      cptr[1] = cptr[2];
      cptr[2] = tmp;
   }
   return(TRUE);
}

/*
 Read a possibly byte swapped unsigned integer
 */
int ReadUInt(FILE *fptr,unsigned int *n,int swap)
{
   unsigned char *cptr,tmp;
   
   if(sizeof(int) != 4)
     {
      fprintf(stderr,"ReadInt: sizeof(int)=%ld and not 4\n",sizeof(int));
      exit(0);
     }
   
   if (fread(n,4,1,fptr) != 1)
      return(FALSE);
   if (swap) {
      cptr = (unsigned char *)n;
      tmp     = cptr[0];
      cptr[0] = cptr[3];
      cptr[3] = tmp;
      tmp     = cptr[1];
      cptr[1] = cptr[2];
      cptr[2] = tmp;
   }
   return(TRUE);
}

/*
 Read a possibly byte swapped long integer
 */
int ReadLong(FILE *fptr,long *n,int swap)
{
   unsigned char *cptr,tmp;
   
   if(sizeof(long) == 4)
     {
      if (fread(n,4,1,fptr) != 1)
         return(FALSE);
      if (swap) {
         cptr = (unsigned char *)n;
         tmp     = cptr[0];
         cptr[0] = cptr[3];
         cptr[3] = tmp;
         tmp     = cptr[1];
         cptr[1] = cptr[2];
         cptr[2] = tmp;
      }
     }
   else if(sizeof(long) == 8)
     {
      if (fread(n,8,1,fptr) != 1)
         return(FALSE);
      if (swap) {
         cptr = (unsigned char *)n;
         tmp     = cptr[0];
         cptr[0] = cptr[7];
         cptr[7] = tmp;
         tmp     = cptr[1];
         cptr[1] = cptr[6];
         cptr[6] = tmp;
         tmp     = cptr[2];
         cptr[2] = cptr[5];
         cptr[5] = tmp;
         tmp     = cptr[3];
         cptr[3] = cptr[4];
         cptr[4] = tmp;
      }
     }
   else
     {
      fprintf(stderr,"ReadLong: something wrong...cannot read long\n");
      exit(0);
     }
   
   
   
   return(TRUE);
}

/*
 Read a possibly byte swapped long long integer
 */
int ReadLongLong(FILE *fptr,long long *n,int swap)
{
   unsigned char *cptr,tmp;
   
   if (fread(n,8,1,fptr) != 1)
      return(FALSE);
   if (swap) {
      cptr = (unsigned char *)n;
      tmp     = cptr[0];
      cptr[0] = cptr[7];
      cptr[7] = tmp;
      tmp     = cptr[1];
      cptr[1] = cptr[6];
      cptr[6] = tmp;
      tmp     = cptr[2];
      cptr[2] = cptr[5];
      cptr[5] = tmp;
      tmp     = cptr[3];
      cptr[3] = cptr[4];
      cptr[4] = tmp;
   }
   return(TRUE);
}

/*
 Read a possibly byte swapped double precision number
 Assume IEEE
 */
int ReadDouble(FILE *fptr,double *n,int swap)
{
   unsigned char *cptr,tmp;
   
   if(sizeof(double) != 8)
     {
      fprintf(stderr,"ReadDouble: sizeof(double)=%ld and not 8\n",sizeof(double));
      exit(0);
     }
   
   if (fread(n,8,1,fptr) != 1)
      return(FALSE);
   if (swap) {
      cptr = (unsigned char *)n;
      tmp     = cptr[0];
      cptr[0] = cptr[7];
      cptr[7] = tmp;
      tmp     = cptr[1];
      cptr[1] = cptr[6];
      cptr[6] = tmp;
      tmp     = cptr[2];
      cptr[2] = cptr[5];
      cptr[5] = tmp;
      tmp     = cptr[3];
      cptr[3] = cptr[4];
      cptr[4] = tmp;
   }
   
   return(TRUE);
}

/*
 Read a possibly byte swapped floating point number
 Assume IEEE format
 */
int ReadFloat(FILE *fptr,float *n, int swap)
{
   unsigned char *cptr,tmp;
   
   if(sizeof(float) != 4)
     {
      fprintf(stderr,"ReadFloat: sizeof(float)=%ld and not 4\n",sizeof(float));
      exit(0);
     }
   
   if (fread(n,4,1,fptr) != 1)
      return(FALSE);
   if (swap) 
     {
      cptr = (unsigned char *)n;
      tmp     = cptr[0];
      cptr[0] = cptr[3];
      cptr[3] = tmp;
      tmp     = cptr[1];
      cptr[1] = cptr[2];
      cptr[2] = tmp;
     }
   return(TRUE);
}
