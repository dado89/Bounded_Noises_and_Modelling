#include <stdio.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include "noises.c"


/* DA COMPILARE CON gcc hist_cai0509.c `pkg-config --cflags --libs gsl`

/* Programma per generare istogrammi dei tempi d'ingresso in [0.5, 1] e [0.9, 1].*/


int main(){

  // COSTANTI GENERALI DEL PROBLEMA DA RISOLVERE

  const double k=0.03;

  // DEFINIZIONE DELLE PRINCIPALI VARIABILI USATE

  double T, dt, T_aux, dt_aux;
  double *X, *Noise, *time, *Noise_aux, *time_aux;
  double tau_c, tau, tau_aux, B, delta, x0;
  double x, y;
  int i, j, i1, i2, seed, seed1, M, N, N_aux;
  FILE *fp;

  M=20000;

  delta=10;
  B=k/2;
  tau_c=1/k;

  tau=10*tau_c;

  dt_aux=0.001;
  tau_aux=1;
  T_aux=10;
  N_aux=T_aux/dt_aux;

  dt=0.003;
  T=120;
  N=T/dt;  

  printf("delta=%g, B=%g, coeff=%g, T=%g.\n", delta, B, tau/tau_c, T);

  Noise=(double*)malloc(2*N*sizeof(double));
  time=(double*)malloc(2*N*sizeof(double));
  X=(double*)malloc(2*N*sizeof(double));
  Noise_aux=(double*)malloc(2*N_aux*sizeof(double));
  time_aux=(double*)malloc(2*N_aux*sizeof(double));


  fp=fopen("..//File txt//Cai//hist_cai0509.txt", "w");

  if(fp!=NULL){
	  for(j=0; j<M; j++){
	
		seed=3*j + 4;
		seed1=3*j + 2;
		  
	  	N_aux=Cai_Lin(Noise_aux, time_aux, tau_aux, B, delta, T_aux, dt_aux, 0, seed1);
		x0=Noise_aux[N_aux];

		N=Cai_Lin(Noise, time, tau, B, delta, T, dt, x0, seed);
	
		Lansky(X, Noise, time, k, N);

		i1=indice(X, 0.5, 0, N);
		i2=indice(X, 0.9, i1, N);
		x=time[i1];
		y=time[i2];

		fprintf(fp, "%g %g\n", x, y);

	  }				/* chiuso for in j */
	}				/* chiuso if */
  else  printf("Errore nell'apertura di hist_cai_0509.txt");

  fclose(fp);


				// LIBERAZIONE MEMORIA

  free(X);
  free(Noise);
  free(time);
  free(Noise_aux);
  free(time_aux);

return 0;

}
