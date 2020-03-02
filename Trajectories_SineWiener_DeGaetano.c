#include <stdio.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include "noises_deg.c"


// DA COMPILARE CON gcc trajectories_cai_deg.c `pkg-config --cflags --libs gsl`

int main(){

  // COSTANTI GENERALI DEL PROBLEMA DA RISOLVERE

  const double k1=0.839;
  const double k2=0.033;
  const double C0=165.5;

  // DEFINIZIONE DELLE PRINCIPALI VARIABILI USATE

  double T, dt, B, dt_brown, tau, tau1, tau2, x0, Discr, lambda1, lambda2, t, w, esp1, esp2;
  double *C, *Q, *time, *Noise, *SolC, *SolQ, *Brown;
  int i, j, N, M, seed1, seed;
  FILE *fp;

  Discr=sqrt( 4*k1*k1 + k2*k2 );
  lambda1=-(2*k1 + k2 - Discr)/2;
  lambda2=-(2*k1 + k2 + Discr)/2;
  tau1=-1/lambda2;
  tau2=-1/lambda1;

  M=25;
  dt=0.001;
  T=180;
  N=T/dt; 

  B=k2;
  tau=tau2;
  dt_brown=0.1;
  
  Noise=(double*)malloc((2*N +1)*sizeof(double));
  time=(double*)malloc((N+1)*sizeof(double));
  C=(double*)malloc((N+1)*sizeof(double));
  Q=(double*)malloc((N+1)*sizeof(double));
  SolC=(double*)malloc((N+1)*sizeof(double));
  SolQ=(double*)malloc((N+1)*sizeof(double));
  Brown=(double*)malloc(11*sizeof(double));

  for(i=0; i<N+1; i++)
	time[i]=i*dt;

  for(i=0; i<N+1; i++){
	t=time[i];
	esp1=exp(lambda1*t);
	esp2=exp(lambda2*t);
	SolC[i]=C0/(2*Discr) * (k2* (esp1 - esp2) + Discr*(esp1 + esp2) );
	SolQ[i]=C0*k1/Discr * (esp1 - esp2);
	}

  printf("M di Matlab = %d\n", M);
  printf("B=%g*k2=%g, tau=%g\n", B/k2, B, tau);

  fp=fopen("..//..//File txt//DeGaetano//Sine//trajectories_sine_deg.txt", "w");
  if(fp!=NULL){

	  for(i=0; i<N+1; i++)
		fprintf(fp, "%.4f\n", time[i]);
	  for(i=0; i<N+1; i++)
		fprintf(fp, "%g\n", SolC[i]);
	  for(i=0; i<N+1; i++)
		fprintf(fp, "%g\n", SolQ[i]);  

	  for(j=0; j<M; j++){
		seed=3*(j+M) +8;
		seed1=3*j + 5;

		Brown1(dt_brown, Brown, seed1);
		w=10*Brown[10];
		Sine_Wien(Noise, tau, B, dt/2, 2*N, w, seed);
		DeGaetano(C, Q, Noise, time, k1, k2, C0, N);

	    for(i=0; i<N+1; i++)
		  fprintf(fp, "%g\n", C[i]);
	    for(i=0; i<N+1; i++)
		  fprintf(fp, "%g\n", Q[i]);
	  }
  }

  else  printf("Errore nell'apertura di trajectories_cai_deg.txt");

  fclose(fp);

/*
  fp=fopen("..//..//File txt//DeGaetano//Sine//trajectories_sine_deg.txt", "w");
  if(fp!=NULL){
  	for(i=0; i<N+1; i++){
		fprintf(fp, "%g %g %g %g %g\n", time[i], C[i], Q[i], SolC[i], SolQ[i]);
	 	  }
	}

  else  printf("Errore nell'apertura di trajectories_sine_deg.txt");

  fclose(fp);
*/

				// LIBERAZIONE MEMORIA

  free(C);
  free(Q);
  free(SolC);
  free(SolQ);
  free(Noise);
  free(time);
  free(Brown);

return 0;

}
