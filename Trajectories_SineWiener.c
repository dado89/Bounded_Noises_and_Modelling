#include <stdio.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include "noises.c"


/* DA COMPILARE CON gcc trajectories_sine.c `pkg-config --cflags --libs gsl`

/* Programma per generare 3 soluzioni del modello Lansky, col rumore Sine-Wiener */


int main(){

  // COSTANTI GENERALI DEL PROBLEMA DA RISOLVERE

  const double k=0.03;

  // DEFINIZIONE DELLE PRINCIPALI VARIABILI USATE

  double T, dt;
  double *X, *Noise, *time;
  double tau_c, tau, B, w, c;
  int i, j, seed, N;
  FILE *fp;

  B=k/2;
  
  tau_c=1/k;

  tau=0.1*tau_c;
  c=10*sqrt(tau);

  dt=0.001;
  T=180;
  N=T/dt;  


  X=(double*)malloc((N+1)*sizeof(double));
  Noise=(double*)malloc((N+1)*sizeof(double));
  time=(double*)malloc((N+1)*sizeof(double));


  for(i=0; i<N+1; i++)
	time[i]=i*dt;

  seed=17;

  fp=fopen("..//File txt//Sine//trajectories_sine.txt", "w");
  if(fp!=NULL){
		  
		w=c*Brown1(2*seed);
		
		printf("B=%g, coeff=%g, x0=%g, seed=%d.\n", B, tau/tau_c, B*sin(sqrt(2/tau)*w), seed);
		Sine_Wien(Noise, tau, B, dt, N, w, seed);
	
		Lansky(X, Noise, time, k, N);
		
		fprintf(fp, "%d\n", N);

		for(i=0; i<N+1; i++){
			fprintf(fp, "%g %g\n", time[i], X[i]);
	 	  }
	}

  else  printf("Errore nell'apertura di trajectories_sine.txt");

  fclose(fp);


				// LIBERAZIONE MEMORIA

  free(X);
  free(Noise);
  free(time);

return 0;

}	
