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

  double sol, T, dt, T_aux, dt_aux, x;
  double *X, *Noise, *time, *Noise_aux, *time_aux;
  double tau_c, tau, tau_aux, B, delta, x0;
  int a, i, j, seed, seed1, M, N, N_aux;
  FILE *fp;

  M=2000;

  delta=10;
  B=k;
  tau_c=1/k;

  tau=tau_c;

  dt_aux=0.001;
  tau_aux=1;
  T_aux=10;
  N_aux=T_aux/dt_aux;

  dt=0.002;
  T=120;
  N=T/dt;  

  printf("M=%d, delta=%g, B=%g, coeff=%g\n", M, delta, B, tau/tau_c);

  Noise=(double*)malloc(2*N*sizeof(double));
  time=(double*)malloc(2*N*sizeof(double));
  X=(double*)malloc(2*N*sizeof(double));

  Noise_aux=(double*)malloc(2*(N_aux)*sizeof(double));
  time_aux=(double*)malloc(2*(N_aux)*sizeof(double));

  a=7;

  printf("a=%d\n", a);
  fp=fopen("..//File txt//Cai//hist_cai_norminf.txt", "w");

  if(fp!=NULL){
	  for(j=0; j<M; j++){
	
		seed=27*j + a;
		seed1= 27*j + 8 + a;
		  
	  	N_aux=Cai_Lin(Noise_aux, time_aux, tau_aux, B, delta, T_aux, dt_aux, 0, seed1);
		x0=Noise_aux[N_aux];
		N=Cai_Lin(Noise, time, tau, B, delta, T, dt, x0, seed);
	
		Lansky(X, Noise, time, k, N);
		
		for(i=0; i<N+1; i++){
		  sol = 1 - exp(- k*time[i]);
		  X[i] -= sol;
	 	}
	
		x=norminf(X, N);
		fprintf(fp, "%g\n", x);

	  }				/* chiuso for in j */
	}				/* chiuso if */
  else  printf("Errore nell'apertura di hist_cai_norm2.txt");

  fclose(fp);


				// LIBERAZIONE MEMORIA

  free(X);
  free(Noise);
  free(time);
  free(Noise_aux);
  free(time_aux);

return 0;

}
