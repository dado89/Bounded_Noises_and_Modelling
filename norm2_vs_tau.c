#include <stdio.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include "noises.c"


/* DA COMPILARE CON gcc norm2_vs_tau.c `pkg-config --cflags --libs gsl`

/* Programma per generare traiettorie della norma2 della differenza tra la soluzione deterministica e quelle stocastiche,
	in funzione del tempo di autocorrelazione */


int main(){

  // COSTANTI GENERALI DEL PROBLEMA DA RISOLVERE

  const double k=0.03;

  // DEFINIZIONE DELLE PRINCIPALI VARIABILI USATE

  double T, sol, dt, dt_aux, tau, tau_aux, T_aux, sigma;
  double *X_cai, *Cai_noise, *X_tsb, *Tsb_noise, *X_sine, *Sine_noise, *X_white, *Noise_aux; 	/* size N+1*/
  double *time_cai, *time_tsb, *time_sw, *time_aux;
  double *Normcai, *Normtsb, *Normsine, *Normwhite;							/* size M */
  double *Time_val, *Sigma_val, *Cai, *Tsb, *Sine, *White; 	  					/* size R */
  double delta, q, B, x0, w0;
  int i, j, h, N_cai, N_tsb, N_aux, N, seed, seed1, M, R;
  FILE *fp;

  M=1000;		
  R=16;		

  T=120; 	
  dt=0.003;
  N=T/dt;

  B=k/2;
  delta=-0.5;
  q=0;

  dt_aux=0.001;
  tau_aux=1;
  T_aux=10;
  N_aux=T_aux/dt_aux;


  // ALLOCAZIONE DELLA MEMORIA

  time_aux=(double*)malloc(2*N_aux*sizeof(double));
  time_cai=(double*)malloc(2*N*sizeof(double));
  time_tsb=(double*)malloc(2*N*sizeof(double));
  time_sw=(double*)malloc((N+1)*sizeof(double));
  
  Noise_aux=(double*)malloc(2*N_aux*sizeof(double));
  Cai_noise=(double*)malloc(2*N*sizeof(double));			// variabile soluzione del noise
  Tsb_noise=(double*)malloc(2*N*sizeof(double));
  Sine_noise=(double*)malloc((N+1)*sizeof(double));

  X_cai=(double*)malloc(2*N*sizeof(double));
  X_tsb=(double*)malloc(2*N*sizeof(double));
  X_sine=(double*)malloc((N+1)*sizeof(double));
  X_white=(double*)malloc((N+1)*sizeof(double));

  Normcai=(double*)malloc(M*sizeof(double));
  Normtsb=(double*)malloc(M*sizeof(double));
  Normsine=(double*)malloc(M*sizeof(double));
  Normwhite=(double*)malloc(M*sizeof(double));

  Cai=(double*)malloc(R*sizeof(double));
  Tsb=(double*)malloc(R*sizeof(double));
  Sine=(double*)malloc(R*sizeof(double));
  White=(double*)malloc(R*sizeof(double));
  Time_val=(double*)malloc(R*sizeof(double));
  Sigma_val=(double*)malloc(R*sizeof(double));


  for(i=0; i<N+1; i++)
	time_sw[i]=i*dt;
 
	// Inizializzazione di Time_val ai valori 
	// 0.01, 0.02, 0.05, 0.08, 0.1, 0.2, 0.5, 0.8, 1, 2, 5, 8, 10, 20, 50, 80
  for(i=0;i<4;i++){
	j=4*i;
	Time_val[j]=0.01*pow(10,i);
	Time_val[j+1]=0.02*pow(10,i);
	Time_val[j+2]=0.05*pow(10,i);
	Time_val[j+3]=0.08*pow(10,i);
  }

	// Inizializzazione di Sigma_val ai valori
	// 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22
  for(i=0; i<10; i++)
	Sigma_val[i]=0.01*(i+1);

  for(i=0; i<6; i++)
	Sigma_val[10+i]=0.1 + 0.02*(i+1);


	// Iniziamo il vero programma

  for(i=0; i<R; i++){
	
		tau=Time_val[i]/k;
		sigma=Sigma_val[i];

		for(j=0; j<M; j++){

			seed=5*j+2;
			seed1= 5*j + 1;

		  	N_aux=Cai_Lin(Noise_aux, time_aux, tau_aux, B, delta, T_aux, dt_aux, 0, seed1);
			x0=Noise_aux[N_aux];
		  	N_cai=Cai_Lin(Cai_noise, time_cai, tau, B, delta, T, dt, x0, seed);

		  	N_aux=Tsall_Borl(Noise_aux, time_aux, tau_aux, B, q, T_aux, dt_aux, 0, seed1);
			x0=Noise_aux[N_aux];
		  	N_tsb=Tsall_Borl(Tsb_noise, time_tsb, tau, B, q, T, dt, x0, seed);

			w0=Brown1(seed1)*sqrt(100*tau);
			Sine_Wien(Sine_noise, tau, B, dt, N, w0, seed);


		  	Lansky_white(X_white, k, sigma, T, N, seed);
			Lansky(X_cai, Cai_noise, time_cai, k, N_cai);
		  	Lansky(X_tsb, Tsb_noise, time_tsb, k, N_tsb);
		  	Lansky(X_sine, Sine_noise, time_sw, k, N);

			for(h=0; h<N_cai+1; h++){
			  sol = 1 - exp(-k*time_cai[h]);
			  X_cai[h] -= sol;}

			for(h=0; h<N_tsb+1; h++){
			  sol = 1 - exp(-k*time_tsb[h]);
			  X_tsb[h] -= sol;}

			for(h=0; h<N+1; h++){
			  sol = 1 - exp(-k*time_sw[h]);
			  X_sine[h] -= sol;
			  X_white[h] -= sol;}

			Normcai[j]=norm2(X_cai, time_cai, N_cai);
			Normtsb[j]=norm2(X_tsb, time_tsb, N_tsb);
			Normsine[j]=norm2(X_sine, time_sw, N);
			Normwhite[j]=norm2(X_white, time_sw, N);

		  } 	// fine for in j. A tau fissato (= i fissato), abbiamo calcolato M norme con quel tau

		Cai[i]=mean(Normcai, M);
		Tsb[i]=mean(Normtsb, M);
		Sine[i]=mean(Normsine, M);		
		White[i]=mean(Normwhite, M);

	}


  fp=fopen("..//File txt//norm2_vs_tau.txt", "w");

  if(fp!=NULL)
	for(i=0; i<R; i++)
		fprintf(fp, "%g %g %g %g %g %g\n", Time_val[i], Sigma_val[i], Cai[i], Tsb[i], Sine[i], White[i]);
  else printf("Errore nell'apertura di norm2_vs_tau.txt");

  fclose(fp);


				// LIBERAZIONE MEMORIA

  free(time_aux);
  free(time_cai);
  free(time_tsb);
  free(time_sw);

  free(Noise_aux);
  free(Cai_noise);
  free(Tsb_noise);
  free(Sine_noise);

  free(X_cai);
  free(X_tsb);
  free(X_sine);
  free(X_white);

  free(Normcai);
  free(Normtsb);
  free(Normsine);
  free(Normwhite);

  free(Cai);
  free(Tsb);
  free(Sine);
  free(White);
  free(Time_val);
  free(Sigma_val);

return 0;

}
