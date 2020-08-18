/* This program calculates the Poisson potential as a function of a given coordinate 
from the charges, dipoles and an xyz positions file where cartesian coordinates are given in Anstroms. */

/* The equation for the Poisson potential can be found for example in:
J. Phys. Chem. C, 115, 16613 (2011) - equation 6 */

/* To execute the program, you just have to launch it with an input file giving
the following information (program < inputfile)
first line: name of the positions' file (> name)
second line: are there dipoles ? 0 if no dipoles, 1 if dipoles (> DIP)
supplementary line (if there are dipoles): name of the dipoles' file /space/ number of columns (> name2)
(3 if just dipoles, 4 if indices in the first column)
third line: number of configurations and number of independent time zones you want (> nconfigs, ntimezones)
fifth line: direction along which the ionic densities will be calculated (x, y or z)
following lines: length of the box in the three cartesian directions (> Lbox)
following line: number of boxes (> nbox)
following line: number of different species (> diffspecies)
following lines: number of atoms/ions for each species (> nspecies)
following line: number of species with non fluctuating charge (> diffmove)
following lines: type and charge for each species with non fluctuating charge (> type, charge)
following line: initial condition for the charge potential psiq(z0) (> psiq0) in volts
following line: initial condition for the dipole potential psimu(z0) (> psimu0) in volts 
(0.0 if there are no dipoles)
following line: number of configurations you want to skip at the beginning of the files (> nskip). 
if there are atoms with gaussian charge densities:
add one line: file containing the fixed atoms charges for each configuration (> name3)
add one line: number of lines with comments to skip for each config (> ncomments)
add one line: width of the gaussian for the metallic atoms in Anstroms (> eta)
add one line: shift to integrate before and after the physical ends of the cell (> shift). */

/* Pay attention to the file dipoles.out: is it indices and dipoles or just dipoles ? */

/* Two solutions for the point charges: if all charges are equal, give the charge; if they differ, put eta = 0.*/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#define L 1000
#define epsilon_zero 8.85418782e-12
/* epsilon_zero is the permittivity of free space in C2.N-1.m-2 */

FILE *in,*in2,*in3;
char name[100],name2[100],name3[100],dirxyz[5];
/* Possibility to have up to 1000 atom types */
char type[1001][5];
int nions,*nspecies,diffspecies,nbox,nconfigs,dir,nskip,ncomments,diffmove,nmove,nshiftbins,ncol,DIP,ntimezones;
double psimu0,psiq0,eta,shift,binwidth,Lbox[4];
double *coord,*qdensity,*dipdensity,*intqdensity,*intintqdensity,*intdipdensity,*poisson,*charge,*wallq,*wallpos;
double **pos,**dip,**Charge,**Dip;

void read_input();
void calc_poisson();
void calc_density();
void zero_matrices();
double erfunc(double x);

char **cmatrix(int nl, int nc);
int **imatrix(int nl, int nc);
double *dvector(int n);
int *ivector(int n);
double **dmatrix(int nl, int nc);


int main ( void )
{ 
 char ligne[L],ligne2[L],ligne3[L];
 int i,j,k,t,nsteps;
 FILE *out;

 printf("Hello\n");
 read_input();

 nsteps=(nconfigs-nskip)/ntimezones;
 printf("Time zones of: %d steps.\n",nsteps);

 /* Opening positions, dipoles and wall charges if all files are used */
 in=fopen(name,"r");
 if(DIP==1) in2=fopen(name2,"r");
 if((diffspecies-diffmove)!=0)  {in3=fopen(name3,"r");}

 /* Skipping useless lines */
 if(nskip!=0)
 {
  for(i=1;i<=nskip;i++) 
     {
      /* Two lines not used at the beginning of each config in the xyz file */
      fgets(ligne,L,in); 
      fgets(ligne,L,in);
      /* ncomments lines not used at the beginning of each config in the wall charges file */ 
      for(k=1;k<=ncomments;k++) fgets(ligne,L,in3);
      for(j=1;j<=nions;j++) 
         {
          fgets(ligne,L,in); 
          if(DIP==1) fgets(ligne2,L,in2);
          if(j>nmove) fgets(ligne3,L,in3);
         }
     }
 }

 /* For each configuration, calculation of the density */	
 for(t=1;t<=ntimezones;t++)
    {
     for(i=((t-1)*nsteps+nskip+1);i<=(t*nsteps+nskip);i++)
        {
         if(i%100==0) printf("%d steps over %d\n",i,nconfigs);
         zero_matrices(); 
         calc_density(); 
         for(j=1;j<=nbox;j++) {Charge[j][t]+=qdensity[j]; Dip[j][t]+=dipdensity[j];}
        }
    }
 
 /* average of the quantities over the configurations */
 for(t=1;t<=ntimezones;t++)
    {
     for(i=1;i<=nbox;i++) 
        {
         Charge[i][t]/=nsteps; 
         Dip[i][t]/=nsteps;
         Charge[i][0]+=Charge[i][t];
         Dip[i][0]+=Dip[i][t];
        }
    }
 
 /* Average over the time zones */
 for(i=1;i<=nbox;i++) 
    {
     Charge[i][0]/=ntimezones;
     Dip[i][0]/=ntimezones;
    }

 /* Calculation of the Poisson Potential and derived quantities from the averages charge and dipole densities */ 
 calc_poisson();

 out=fopen("poisson.dat","w");
 for(i=1;i<=nbox;i++) fprintf(out,"%e %e\n",coord[i]*1e10,poisson[i]);
 fclose(out);
 
 out=fopen("charge.dat","w");
 for(i=1;i<=nbox;i++) fprintf(out,"%e %e\n",coord[i]*1e9,Charge[i][0]/(0.1*0.1*0.1));
 fclose(out); 

 out=fopen("int2charge.dat","w");
 for(i=1;i<=nbox;i++) fprintf(out,"%e %e\n",coord[i]*1e10,intintqdensity[i]);
 fclose(out); 

 if(DIP==1) 
    {
     out=fopen("dipole.dat","w");
     for(i=1;i<=nbox;i++) fprintf(out,"%e %e\n",coord[i]*1e9,Dip[i][0]);
     fclose(out);
    }

 out=fopen("intdipole.dat","w");
 for(i=1;i<=nbox;i++) fprintf(out,"%e %e\n",coord[i]*1e10,intdipdensity[i]);
 fclose(out);

 fclose(in);
 if(DIP==1) fclose(in2);
 if((diffspecies-diffmove)!=0) {fclose(in3);}
 
 return 0;
}


/*Reads the input file and allocates the vectors and matrices*/
void read_input()
{
 int i;
 
 scanf("%s",name);					/* name of the positions' file */
 scanf("%d",&DIP);
 if(DIP==1) scanf("%s %d",name2,&ncol);			/* name of the dipoles' file and number of columns in it*/
 scanf("%d %d",&nconfigs,&ntimezones);  		/* number of configurations and number of time zones */	
 scanf("%s",dirxyz);	    				/* direction of divisions */ 
 scanf("%lf",&Lbox[1]);	   				/* length of the simulation box x */
 scanf("%lf",&Lbox[2]);	    				/* length of the simulation box y */
 scanf("%lf",&Lbox[3]);	     				/* length of the simulation box z */
 scanf("%d",&nbox);	     				/* number of boxes in the chosen direction */
 scanf("%d",&diffspecies);  				/* number of different species */	
 nspecies=ivector(diffspecies+1);
 for(i=1;i<=diffspecies;i++) scanf("%d",&nspecies[i]); 	/* number of ions for each species */
 scanf("%d",&diffmove);  				/* number of species allowed to move */	
 charge=dvector(diffmove+1);
 for(i=1;i<=diffmove;i++) scanf("%s %lf",type[i],&charge[i]);
 scanf("%lf",&psiq0);	
 scanf("%lf",&psimu0);					/* initial conditions */
 scanf("%d",&nskip);					/* number of configurations to skip */

 /* If there are fixed atoms with gaussian and fluctuating charges */
 if((diffspecies-diffmove)!=0) {scanf("%s",name3);	/* name of the wall charges' file */
 				scanf("%d",&ncomments);	/* number of lines with comments to skip for each config */
 				scanf("%lf",&eta);	/* width of the gaussian charge for wall atoms */
			        scanf("%lf",&shift);}	/* shift to integrate before and after the ends of the cell */
 if((diffspecies-diffmove)==0) shift=0;			/* if no gaussian species, no need to integrate before and after the limits of the cell */

 if(eta==0.0) printf("Point charges for the walls !\n"); 
 
 for(i=1;i<=diffspecies;i++) nions+=nspecies[i];	/* Calculation of the total number of atoms/ions */
 for(i=1;i<=diffmove;i++) nmove+=nspecies[i];		/* Calculation of the total number of moving atoms/ions */
 printf("Reading OK !\n"); 

 /* Transformation of the direction into a given number */
 if(strcmp(dirxyz,"x")==0) dir=1;
 if(strcmp(dirxyz,"y")==0) dir=2;
 if(strcmp(dirxyz,"z")==0) dir=3;

 binwidth=Lbox[dir]/nbox;
 nshiftbins=0;

 /* Calculates the number of bins added */
 if((diffspecies-diffmove)!=0)
  {
   nshiftbins=2*floor(shift/binwidth);
   shift=0.5*nshiftbins*binwidth;
   printf("Adding %d bins.\n",nshiftbins);
   nbox+=nshiftbins;
   Lbox[dir]+=2.0*shift;
  }

 /* Allocation of the matrices and vectors */
 coord=dvector(nbox+1);
 qdensity=dvector(nbox+1); intqdensity=dvector(nbox+1); poisson=dvector(nbox+1);
 dipdensity=dvector(nbox+1); 
 pos=dmatrix(4,nions+1), dip=dmatrix(4,nions+1); 
 Charge=dmatrix(nbox+1,ntimezones+1); Dip=dmatrix(nbox+1,ntimezones+1);
 wallq=dvector(nions+1); wallpos=dvector(nions+1);
 intintqdensity=dvector(nbox+1); 
 intdipdensity=dvector(nbox+1); 
 printf("Allocation OK !\n");
}


/* Calculates the Poisson potential and writes it in poisson.dat */
void calc_poisson()
{
 int i,j;

 /* First step : calculation of the integral of the charge density (intqdensity) 
    using trapezian method */
 intqdensity[1]=0;
 for(i=2;i<=nbox;i++) intqdensity[i]=intqdensity[i-1]+(Charge[i][0]+Charge[i-1][0])*(coord[i]-coord[i-1])/2;

 /* Second step : calculation of the poisson potential using trapezian method */
 /* Conversion of coord, intqdensity and dipdensity in SI units. */
 for(i=1;i<=nbox;i++)
    {
     intqdensity[i]*=1.6e-19/(1e-10*1e-10);
     coord[i]*=(1e-10);
     Dip[i][0]*=1.6e-19/(1e-10*1e-10);
    }

 /* Integrals of intqdensity and dipole density */
 intintqdensity[1]=psiq0; intdipdensity[1]=psimu0;
 for(i=2;i<=nbox;i++) 
    {
     intintqdensity[i]=intintqdensity[i-1]+(-1)*(intqdensity[i]+intqdensity[i-1])*(coord[i]-coord[i-1])/(2*epsilon_zero);
     intdipdensity[i]=intdipdensity[i-1]+(Dip[i][0]+Dip[i-1][0])*(coord[i]-coord[i-1])/(2*epsilon_zero);
    }

 /* Poisson potential */
 for(i=1;i<=nbox;i++) poisson[i]=intintqdensity[i]+intdipdensity[i];
}


/* Calculates the charge and dipole densities */
void calc_density()
{
 char ligne[L],ligne2[L],ligne3[L];   
 char temp[5];
 int i,j,k;
 double *qtype;

 qtype=dvector(nions+1);

 /* Calculation of the coordinates */
 for(i=1;i<=nbox;i++) coord[i]=(i-0.5)*binwidth-shift;

 /* Jump the two lines at the beginning of the xyz file */
 fgets(ligne,L,in); 
 fgets(ligne,L,in); 

 /* Reading of positions, dipoles and charges */
 for(j=1;j<=nions;j++) 
  {
   fgets(ligne,L,in); 
   if(DIP==1) fgets(ligne2,L,in2);
   /* If the input files are written with D instead of e for the scientific notation, these loops change it.*/
   for(i=0;i<=strlen(ligne);i++)
      {
       if(ligne[i] == 'D') ligne[i] = 'e';
      }   
   if(DIP==1)
     {
      for(i=0;i<=strlen(ligne2);i++)
         {
          if(ligne2[i] == 'D') ligne2[i] = 'e';
         }
     }
   /* Read positions and find/assign charge */
   sscanf(ligne,"%s %lf %lf %lf",temp,&pos[1][j],&pos[2][j],&pos[3][j]);
   for(i=1;i<=diffmove;i++)
      {
       if(strcmp(type[i],temp)==0) qtype[j]=charge[i];
      }
   /* Read dipoles */
   if(DIP==1)
	{
   	 if(ncol==3) sscanf(ligne2,"%lf %lf %lf",&dip[1][j],&dip[2][j],&dip[3][j]); 
   	 if(ncol==4) sscanf(ligne2,"%*d %lf %lf %lf",&dip[1][j],&dip[2][j],&dip[3][j]); 
	}
   if(DIP==0) 
	{
	 dip[1][j]=dip[2][j]=dip[3][j]=0.0;
   	}
   /* Read charges */
   if(j>nmove)  
	{
	  fgets(ligne3,L,in3);      
	  for(i=0;i<=strlen(ligne3);i++)
	    {
	     if(ligne3[i] == 'D') ligne3[i] = 'e';
	    }
	  sscanf(ligne3,"%le",&wallq[j]); 
	}
  }

 /* Calculation of charge and dipole densities */ 
 for(j=1;j<=nions;j++) 
  {
   k=1;
   while((pos[dir][j]>(coord[k]+0.5*binwidth))&&(k<nbox)) k++;  /* Find the box in which the ion is. */
   dipdensity[k]+=dip[dir][j];
   if(j<=nmove) qdensity[k]+=qtype[j]; 
   if((j>nmove)&&(eta!=0.0)) {wallpos[j]=coord[k];} 	/* If it is a fixed atom, remember where it is. */
   if((j>nmove)&&(eta==0.0)) {qdensity[k]+=wallq[j];} 	/* Fixed atom with a point charge. */
  }

 /* Calculation of the charge density due the walls */
 if(((diffspecies-diffmove)!=0)&&(eta!=0))
   {
    for(k=1;k<=nbox;k++)
       {
        for(j=nmove+1;j<=nions;j++) qdensity[k]+=0.5*wallq[j]*(erfunc(eta*(wallpos[j]-coord[k]+binwidth*0.5))-erfunc(eta*(wallpos[j]-coord[k]-binwidth*0.5)));
       }
   }

 /* Division by the volume to obtain the density */
 for(k=1;k<=nbox;k++) 
    {
     qdensity[k]/=(Lbox[1]*Lbox[2]*Lbox[3]/nbox); dipdensity[k]/=(Lbox[1]*Lbox[2]*Lbox[3]/nbox);
    }
 
}


/* Put all the matrices to zero */
void zero_matrices()
{
 int i;

 for(i=0;i<=nbox;i++)
  {
   qdensity[i]=0; dipdensity[i]=0;
   poisson[i]=0;
   intqdensity[i]=0;
  }

}


/* Erfunc function, used in the calculation of the gaussian shapes of the charge density */
/* Taken from Abramowitz and Stegun, p.299, section 7.1.28, Dover books 1965 */
/* See also Hastings, Approximations for digital computers, 1955 */
double erfunc(double x)
{
 float xe,er,ye;

 xe=fabs(x);
 ye=1.0-(1.0/pow(1.0+xe*(0.0705230745+xe*(0.0422820123+xe*(0.0092705272+xe*(0.0001520143+xe*(0.0002765672+xe*0.0000430638))))),16));
 if(x>=0) er=fabs(ye);
 if(x<0) er=-fabs(ye);

 return er;
}


/********************************************************************************************/
/* Dynamic allocation */

char **cmatrix(int nl, int nc)
{
 int i;
 char **m;
 m=(char **) malloc(nl*sizeof(char*));
 if (m) {m[0]=(char *) malloc(nl*nc*sizeof(char));
         if (m[0]==NULL) return NULL;
         for (i=1;i<nl;i++) m[i]=m[i-1]+nc;
        }
 return m;
}


int **imatrix(int nl, int nc)
{
 int i;
 int **m;
 m=(int **) malloc(nl*sizeof(int*));
 if (m) {m[0]=(int *) malloc(nl*nc*sizeof(int));
         if (m[0]==NULL) return NULL;
         for (i=1;i<nl;i++) m[i]=m[i-1]+nc;
        }
 return m;
}

double *dvector(int n)
{
 return (double *) malloc(n*sizeof(double));
}

int *ivector(int n)
{
 return (int *) malloc(n*sizeof(int));
}

double **dmatrix(int nl, int nc)
{
 int i;
 double **m;
 m=(double **) malloc(nl*sizeof(double*));
 if (m) {m[0]=(double *) malloc(nl*nc*sizeof(double));
         if (m[0]==NULL) return NULL;
         for (i=1;i<nl;i++) m[i]=m[i-1]+nc;
        }
 return m;
}






