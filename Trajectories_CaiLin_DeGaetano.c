#include <stdio.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include "noises_deg.c"


// DA COMPILARE CON gcc trajectories_cai_deg.c `pkg-config --cflags --libs gsl`

int main(){

  // COSTANTI GENERALI DEL PROBLEMA DA RISOLVERE

  const double k1=0.044;
  const double k2=0.016;
  const double C0=86.4;

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

  // DEFINIZIONE DELLE PRINCIPALI VARIABILI USATE

  double T, dt, delta, B, tau, tau1, tau2, dt_aux, tau_aux, T_aux, x0, Discr, lambda1, lambda2, t, esp1, esp2, meanx0;
  double *C, *Q, *Noise, *times, *Noise_aux, *SolC, *SolQ;
  int i, j, M, N, N_aux, seed1, seed;
  FILE *fp;

  Discr=sqrt( 4*k1*k1 + k2*k2 );
  lambda1=-(2*k1 + k2 - Discr)/2;
  lambda2=-(2*k1 + k2 + Discr)/2;
  tau1=-1/lambda2;
  tau2=-1/lambda1;

printf("tau1=%g, tau2=%g\n", tau1, tau2);

  M=500;
  delta=1;
  B=0.8*k2;
  tau=10*tau2;

  dt=0.008;
  T=180;
  N=T/dt; 

  printf("M di Matlab = %d\n", M);
  printf("B=%g*k2=%g, delta=%g\n", B/k2, B, delta);
  printf("tau=%g\n", tau);
 
  dt_aux=0.001;
  tau_aux=1;
  T_aux=10;
  N_aux=T_aux/dt_aux;

  Noise=(double*)malloc((2*N +1)*sizeof(double));
  times=(double*)malloc((N+1)*sizeof(double));
  Noise_aux=(double*)malloc((N_aux+1)*sizeof(double));

  C=(double*)malloc((N+1)*sizeof(double));
  Q=(double*)malloc((N+1)*sizeof(double));
  SolC=(double*)malloc((N+1)*sizeof(double));
  SolQ=(double*)malloc((N+1)*sizeof(double));


  for(i=0; i<N+1; i++)
	times[i]=i*dt;

  for(i=0; i<N+1; i++){
	t=times[i];
	esp1=exp(lambda1*t);
	esp2=exp(lambda2*t);
	SolC[i]=C0/(2*Discr) * (k2* (esp1 - esp2) + Discr*(esp1 + esp2) );
	SolQ[i]=C0*k1/Discr * (esp1 - esp2);
	}

  meanx0=0;
srand(time(NULL));

  fp=fopen("..//..//File txt//DeGaetano//Cai//trajectories_cai_deg.txt", "w");
  if(fp!=NULL){

	  for(i=0; i<N+1; i++)
		fprintf(fp, "%.4f\n", times[i]);
	  for(i=0; i<N+1; i++)
		fprintf(fp, "%g\n", SolC[i]);
	  for(i=0; i<N+1; i++)
		fprintf(fp, "%g\n", SolQ[i]);

	  for(j=0; j<M; j++){
//		seed=3*(j+M) +8;
//		seed1=3*j + 5;

		seed=rand();
		seed1=rand();

		Cai_Lin(Noise_aux, tau_aux, B, delta, T_aux, N_aux, 0, seed1);
		x0=Noise_aux[N_aux];
		meanx0 += x0;
		Cai_Lin(Noise, tau, B, delta, T, 2*N, x0, seed);
		DeGaetano(C, Q, Noise, dt, k1, k2, C0, N);

	    for(i=0; i<N+1; i++)
		  fprintf(fp, "%g\n", C[i]);
	    for(i=0; i<N+1; i++)
		  fprintf(fp, "%g\n", Q[i]);
	  }
  }

  else  printf("Errore nell'apertura di trajectories_cai_deg.txt");

  fclose(fp);

  meanx0 /= M;
printf("media di x0 rispetto a k2: %g\n", meanx0/k2);


				// LIBERAZIONE MEMORIA

  free(C);
  free(Q);
  free(SolC);
  free(SolQ);
  free(Noise);
  free(times);
  free(Noise_aux);

return 0;

}
