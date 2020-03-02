#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include "noises_deg.c"


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

  const double k1=0.044;
  const double k2=0.016;
  const double C0=86.4;

  FILE *fp;
  double *Noise, *C, *Q, *P, *U, *times, *SolC, *SolQ, *P_in;
  double *step_param, *domain_param, *K_param;
  double dt, tau, B, delta, T, T0, Tf, x0, en1, en2, x, En_fin;
  double t, esp1, esp2, Discr, lambda1, lambda2, tau1, tau2;
  int i, j, N, M, seed, K, N_u, ind;

  dt=0.001;
  T=100;
  N=T/dt;
  N_u=20;

  K=100;

  Noise=(double*)malloc((2*N +1)*sizeof(double));
  C=(double*)malloc((N+1)*sizeof(double));
  Q=(double*)malloc((N+1)*sizeof(double));
  SolC=(double*)malloc((N+1)*sizeof(double));
  SolQ=(double*)malloc((N+1)*sizeof(double));
  times=(double*)malloc((N+1)*sizeof(double));
  P=(double*)malloc(4*sizeof(double));					// B, delta, tau, x0
  P_in=(double*)malloc(4*sizeof(double));		
  step_param=(double*)malloc(3*sizeof(double));			// B, delta, tau
  domain_param=(double*)malloc(6*sizeof(double));		// B, delta, tau
  K_param=(double*)malloc(3*sizeof(double));			// k1, k2, C0
  U=(double*)malloc(2*N_u*sizeof(double));

  K_param[0]=k1;
  K_param[1]=k2;
  K_param[2]=C0;

// Ordine informazioni: B, delta, tau.
  domain_param[0]=0;
  domain_param[1]=0;
  domain_param[2]=0.1;
  domain_param[3]=k2;
  domain_param[4]=3;
  domain_param[5]=100;

  step_param[0]=0.1*k2;	//B
  step_param[1]=0.05;		//delta
  step_param[2]=1;			//tau

  for(i=0; i<(N+1); i++)
	times[i]=i*dt;

//  srand(time(NULL));

  B=k2;
  delta=1;
  tau=5;
  x0=k2/5;

  seed=91;

  Cai_Lin(Noise, tau, B, delta, T, 2*N, x0, seed);
  DeGaetano(C, Q, Noise, dt, k1, k2, C0, N);
/*
  Discr=sqrt( 4*k1*k1 + k2*k2 );
  lambda1=-(2*k1 + k2 - Discr)/2;
  lambda2=-(2*k1 + k2 + Discr)/2;
  tau1=-1/lambda2;
  tau2=-1/lambda1;

  for(i=0; i<N+1; i++){
	t=times[i];
	esp1=exp(lambda1*t);
	esp2=exp(lambda2*t);
	SolC[i]=C0/(2*Discr) * (k2*(esp1 - esp2) + Discr*(esp1 + esp2) );
	SolQ[i]=C0*k1/Discr * (esp1 - esp2);
	}
*/

  U[0]=0.5;
  for(i=1; i<N_u; i++)
	U[i]=5*i;

  for(i=0; i<N_u; i++){
	ind=U[i]/dt;
	U[N_u + i]=C[ind];
//	U[N_u + i]=Q[ind];
  }



  T0=3;
  Tf=0.5;
  SA(U, N_u, P_in, K, T0, Tf, step_param, domain_param, K_param, P, &En_fin);

  printf("\nValori inziali: B=%g*k2, delta=%g, tau=%g, x0=%g*B\n", P_in[0]/k2, P_in[1], P_in[2], P_in[3]/P_in[0]);
  printf("Valori finali:  B=%g*k2, delta=%g, tau=%g, x0=%g*B\n", P[0]/k2, P[1], P[2], P[3]/P[0]);
  //printf("\nEnergia finale=%g\n", En_fin);

  free(C);
  free(Q);
  free(SolC);
  free(SolQ);
  free(times);
  free(Noise);
  free(P);
  free(P_in);
  free(U);
  free(step_param);
  free(domain_param);
  free(K_param);

return 0;
}







/*
  srand(1);
  fp=fopen("..//..//File txt//random_unif.txt", "w");

  for(i=0; i<1000000; i++){
	en=random_unif();
	fprintf(fp, "%g\n", en);
  }

  fclose(fp);
*/
