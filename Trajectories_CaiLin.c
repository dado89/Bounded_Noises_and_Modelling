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

  double T, dt, T_aux, dt_aux, tau_aux;
  double *X, *Noise, *time, *Noise_aux, *time_aux;
  double tau_c, tau, B, delta, x0;
  int i, j, seed, N, N_aux;
  FILE *fp;

  delta=1;
  
  B=k/2;
  
  seed=4; 

  tau_c=1/k;

  tau=0.1*tau_c;
  
  dt_aux=0.001;
  tau_aux=1;
  T_aux=50;
  N_aux=T_aux/dt_aux;

  dt=0.001;
  T=180;
  N=T/dt;  


  Noise=(double*)malloc(2*N*sizeof(double));
  time=(double*)malloc(2*N*sizeof(double));
  X=(double*)malloc(2*N*sizeof(double));

  Noise_aux=(double*)malloc(2*(N_aux)*sizeof(double));
  time_aux=(double*)malloc(2*(N_aux)*sizeof(double));

  printf("B=%g, coeff=%g, delta=%g.\n", B, tau/tau_c, delta);

  fp=fopen("..//File txt//Cai//trajectories_cai.txt", "w");
  if(fp!=NULL){
		  
	  	N_aux=Cai_Lin(Noise_aux, time_aux, tau_aux, B, delta, T_aux, dt_aux, 0, 15); 
/*54:molto>0,  1,11,15,52:poco>0, 3,5,7,51:poco<0, 4,6,8,12:molto<0*/
		x0=Noise_aux[N_aux];

		printf("x0=%g\n", x0);

		N=Cai_Lin(Noise, time, tau, B, delta, T, dt, x0, seed);
	
		Lansky(X, Noise, time, k, N);

		fprintf(fp, "%d\n", N);

	  	for(i=0; i<N+1; i++){
			fprintf(fp, "%g %g\n", time[i], X[i]);
	 	  }
	}

  else  printf("Errore nell'apertura di trajectories_cai.txt");

  fclose(fp);


				// LIBERAZIONE MEMORIA

  free(X);
  free(Noise);
  free(time);
  free(Noise_aux);
  free(time_aux);

return 0;

}	
