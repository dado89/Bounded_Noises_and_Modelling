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

  double T, dt, q, B, tau, tau1, tau2, dt_aux, tau_aux, T_aux, x0, Discr, lambda1, lambda2, t, esp1, esp2;
  double *C, *Q, *Noise, *time, *Noise_aux, *SolC, *SolQ;
  int i, j, M, N, N_aux, seed1, seed;
  FILE *fp;

  Discr=sqrt( 4*k1*k1 + k2*k2 );
  lambda1=-(2*k1 + k2 - Discr)/2;
  lambda2=-(2*k1 + k2 + Discr)/2;
  tau_c1=-1/lambda2;
  tau_c2=-1/lambda1;

  M=10;
  q=0.5;
  B=0.5*k2;
  tau=0.1*tau1;

  dt=0.0005;
  T=180;
  N=T/dt; 

  printf("M di Matlab = %d\n", M);
  printf("B=%g*k2=%g, q=%g\n", B/k2, B, q);
  printf("tau=%g\n", tau);
 
  dt_aux=0.001;
  tau_aux=1;
  T_aux=10;
  N_aux=T_aux/dt_aux;

  Noise=(double*)malloc((2*N +1)*sizeof(double));
  time=(double*)malloc((N+1)*sizeof(double));
  Noise_aux=(double*)malloc((N_aux+1)*sizeof(double));

  C=(double*)malloc((N+1)*sizeof(double));
  Q=(double*)malloc((N+1)*sizeof(double));
  SolC=(double*)malloc((N+1)*sizeof(double));
  SolQ=(double*)malloc((N+1)*sizeof(double));


  for(i=0; i<N+1; i++)
	time[i]=i*dt;

  for(i=0; i<N+1; i++){
	t=time[i];
	esp1=exp(lambda1*t);
	esp2=exp(lambda2*t);
	SolC[i]=C0/(2*Discr) * (k2* (esp1 - esp2) + Discr*(esp1 + esp2) );
	SolQ[i]=C0*k1/Discr * (esp1 - esp2);
	}

  fp=fopen("..//..//File txt//DeGaetano//Tsallis//trajectories_tsb_deg.txt", "w");
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

		Tsall_Borl(Noise_aux, tau_aux, B, q, T_aux, N_aux, 0, seed1);
		x0=Noise_aux[N_aux];
		Tsall_Borl(Noise, tau, B, q, T, 2*N, x0, seed);
		DeGaetano(C, Q, Noise, time, k1, k2, C0, N);

	    for(i=0; i<N+1; i++)
		  fprintf(fp, "%g\n", C[i]);
	    for(i=0; i<N+1; i++)
		  fprintf(fp, "%g\n", Q[i]);
	  }
  }

  else  printf("Errore nell'apertura di trajectories_tsb_deg.txt");

  fclose(fp);


				// LIBERAZIONE MEMORIA

  free(C);
  free(Q);
  free(SolC);
  free(SolQ);
  free(Noise);
  free(time);
  free(Noise_aux);

return 0;

}	
