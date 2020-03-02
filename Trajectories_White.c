#include <stdio.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include "noises.c"


/* DA COMPILARE CON gcc hist_cai_norm2.c `pkg-config --cflags --libs gsl`

/* Programma per generare istogrammi della norma2 della differenza tra la soluzione deterministica e quelle stocastiche */


int main(){

  // COSTANTI GENERALI DEL PROBLEMA DA RISOLVERE

  const double k=0.03;

  // DEFINIZIONE DELLE PRINCIPALI VARIABILI USATE

  double T, dt, sigma;
  double *X;
  int i, seed, N;
  FILE *fp;

  seed=26;

  sigma=0.15;

  dt=0.001;
  T=180;
  N=T/dt;  

  printf("sigma=%g, seed=%d.\n", sigma, seed);

  X=(double*)malloc((N+1)*sizeof(double));

  fp=fopen("..//File txt//White//trajectories_white.txt", "w");
  if(fp!=NULL){
		  
	Lansky_white(X, k, sigma, T, N, seed);
	fprintf(fp, "%d\n", N);	
	for(i=0; i<N+1; i++)
		fprintf(fp, "%g %g\n", dt*i, X[i]);

  }

  else  printf("Errore nell'apertura di trajectories_white.txt");

  fclose(fp);

  free(X);

return 0;

}
