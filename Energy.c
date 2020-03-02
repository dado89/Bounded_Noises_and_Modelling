#include <stdio.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include "noises_deg.c"


// DA COMPILARE CON gcc trajectories_cai_deg.c `pkg-config --cflags --libs gsl`

int main(){

  // Pharmacokinetic parameters of the model (dissolution rates and initial concentration)

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

  // Declaration of main variables used

  double T, dt, delta, B, tau, tau1, tau2, x0, En;
  double *C, *Q, *Noise, *times, *P, *U;
  int i, N, N_u, K, ind, seed;
  FILE *fp;

  N_u=20;

  dt=0.001;
  T=100;
  N=T/dt;
 
  // Allocate memory
	
  Noise=(double*)malloc((2*N +1)*sizeof(double));
  times=(double*)malloc((N+1)*sizeof(double));

  C=(double*)malloc((N+1)*sizeof(double));
  Q=(double*)malloc((N+1)*sizeof(double));
  P=(double*)malloc(4*sizeof(double));			// k1, k2, C0
  U=(double*)malloc(2*N_u*sizeof(double));

  // Solve the system of equations and compute energy
  B=k2;
  delta=1;
  tau=5;
  x0=k2/5;

  K=100;

  for(i=0; i<N+1; i++)
	times[i]=i*dt;

  seed=90;

  Cai_Lin(Noise, tau, B, delta, T, 2*N, x0, seed);
  DeGaetano(C, Q, Noise, dt, k1, k2, C0, N);

  U[0]=0.5;
  for(i=1; i<N_u; i++)
	U[i]=5*i;

  for(i=0; i<N_u; i++){
	ind=U[i]/dt;
	U[N_u + i]=C[ind];
  }

  for(i=0; i<N_u; i++){
	ind=U[i]/dt;
	U[N_u + i]=Q[ind];
  
  B=0.133791*k2;
  delta=0.366798;
  tau=90.442;
  x0=0.102216*B;

  dt=0.005;

  P[0]=B;
  P[1]=delta;
  P[2]=tau;
  P[3]=x0;

  srand(time(NULL));
//  srand(5);
for(i=0; i<4; i++)
seed=rand();

En=Energy_min_dt(P, U, N_u, K, dt, k1, k2, C0);
printf("Energia=%g\n", En);


  free(C);
  free(Q);
  free(Noise);
  free(times);
  free(P);
  free(U);

return 0;

}	
