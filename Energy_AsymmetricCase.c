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


  // DEFINIZIONE DELLE PRINCIPALI VARIABILI USATE

  double T, dt, delta, B, tau, tau1, tau2, x0, En, alpha;
  double Zr, sx, x, Mean;
  double *C, *Q, *Noise, *times, *P, *U;
  int i, N, N_u, K, ind, seed;
  FILE *fp;

  N_u=20;

  dt=0.001;
  T=100;
  N=T/dt;

  K=100;
 
  Noise=(double*)malloc((2*N +1)*sizeof(double));
  times=(double*)malloc((N+1)*sizeof(double));

  C=(double*)malloc((N+1)*sizeof(double));
  Q=(double*)malloc((N+1)*sizeof(double));
  P=(double*)malloc(4*sizeof(double));			// k1, k2, C0
  U=(double*)malloc(2*N_u*sizeof(double));


  B=1;

  Zr=4;
  delta=1;
  tau=5;
  x0=k2/5;

  alpha=alpha_Zr1(delta, Zr, &Mean);
  sx=1/(2+alpha);

  seed=91;

  Cai_Lin(Noise, tau, B, delta, T, 2*N, x0, seed);
  for(i=0; i<(2*N +1); i++){
	x=phi1(Noise[i], alpha);
	Noise[i]=k2 * (x - Mean)/(Mean - sx);	
   }

  DeGaetano(C, Q, Noise, dt, k1, k2, C0, N);

 /* for(i=0; i<N+1; i++)
	times[i]=i*dt;*/



  U[0]=0.5;
  for(i=1; i<N_u; i++)
	U[i]=5*i;

  for(i=0; i<N_u; i++){
	ind=U[i]/dt;
	U[N_u + i]=C[ind];
	}


  alpha=0.3086;
  delta=1.42601;
  tau=93.256;
  x0=0.95;

  dt=0.01;

  P[0]=alpha;
  P[1]=delta;
  P[2]=tau;
  P[3]=x0;

  srand(time(NULL));
//  srand(5);
for(i=0; i<4; i++)
seed=rand();

  En=Energy_min_asymm(P, U, N_u, K, dt, k1, k2, C0);
printf("Energia=%g\n", En);


  free(C);
  free(Q);
  free(Noise);
  free(times);
  free(P);
  free(U);

return 0;

}
