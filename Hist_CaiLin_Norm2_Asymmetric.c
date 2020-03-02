#include <stdio.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include "noises_deg.c"


/* DA COMPILARE CON gcc hist_cai0509.c `pkg-config --cflags --libs gsl`

/* Programma per generare istogrammi dei tempi d'ingresso in [0.5, 1] e [0.9, 1].*/

/*
RAT1
  const double k1=0.839;
  const double k2=0.033;
  const double C0=165.5;

RAT2
  const double k1=0.046;
  const double k2=0.062;
  const double C0=86.3;

RAT5
  const double k1=0.044;
  const double k2=0.016;
  const double C0=86.4;
*/

int main(){

  // COSTANTI GENERALI DEL PROBLEMA DA RISOLVERE

  const double k1=0.839;
  const double k2=0.033;
  const double C0=165.5;

  // DEFINIZIONE DELLE PRINCIPALI VARIABILI USATE

  double T, dt, delta, B, tau, tau1, tau2, dt_aux, tau_aux, T_aux, x0, Discr, lambda1, lambda2, esp1, esp2, t, x, y;
  double *C, *Q, *time, *Noise, *Noise_aux, *time_aux, *SolC, *SolQ;
  int i, j, N, N_aux, M, seed1, seed, i1, i2;
  FILE *fp;

  M=20000;
  delta=1;
  B=1;

  Discr=sqrt( 4*k1*k1 + k2*k2 );
  lambda1=-(2*k1 + k2 - Discr)/2;
  lambda2=-(2*k1 + k2 + Discr)/2;
  tau1=-1/lambda2;
  tau2=-1/lambda1;

  dt=0.002;
  T=120;
  N=T/dt;
  tau=10*tau2;

  printf("M=%d\ntau=%g, B=%g*k2, T=%g\n", M, tau, B/k2, T);
  printf("dt=%g\n", dt);
  
  dt_aux=0.001;
  tau_aux=1;
  T_aux=10;
  N_aux=T_aux/dt_aux;


  Noise=(double*)malloc((2*N +1)*sizeof(double));
  time=(double*)malloc((N+1)*sizeof(double));
  C=(double*)malloc((N+1)*sizeof(double));
  Q=(double*)malloc((N+1)*sizeof(double));
  SolC=(double*)malloc((N+1)*sizeof(double));
  SolQ=(double*)malloc((N+1)*sizeof(double));

  Noise_aux=(double*)malloc((N_aux +1)*sizeof(double));
  time_aux=(double*)malloc((N_aux +1)*sizeof(double));

  for(i=0; i<N+1; i++)
	time[i]=i*dt;

  for(i=0; i<N+1; i++){
	t=time[i];
	esp1=exp(lambda1*t);
	esp2=exp(lambda2*t);
	SolC[i]=C0/(2*Discr) * (k2* (esp1 - esp2) + Discr*(esp1 + esp2) );
	SolQ[i]=C0*k1/Discr * (esp1 - esp2);
	}

  fp=fopen("..//..//File txt//DeGaetano//Cai//hist_cai_norm2_deg.txt", "w");
  if(fp!=NULL){
	  for(j=0; j<M; j++){

	  seed=3*j +7;
	  seed1=4*j+8;
	  
	  Cai_Lin(Noise_aux, tau_aux, B, delta, T_aux, N_aux, 0, seed1); 
	  x0=Noise_aux[N_aux];
	  Cai_Lin(Noise, tau, B, delta, T, 2*N, x0, seed);

	  DeGaetano(C, Q, Noise, time, k1, k2, C0, N);

	  for(i=0; i<N+1; i++){
		  C[i] -= SolC[i];
		  Q[i] -= SolQ[i];
	 	}
	  x=norm2(C, time, N);
	  y=norm2(Q, time, N);
	  fprintf(fp, "%g %g\n", x, y);
	}
  }

  else  printf("Errore nell'apertura di hist_cai_norm2_deg.txt");

  fclose(fp);


				// LIBERAZIONE MEMORIA

  free(C);
  free(Q);
  free(SolC);
  free(SolQ);
  free(Noise);
  free(time);
  free(Noise_aux);
  free(time_aux);

return 0;

}	
