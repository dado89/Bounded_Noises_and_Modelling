#include <stdio.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include "noises.c"


/* DA COMPILARE CON gcc norm2_vs_tau.c `pkg-config --cflags --libs gsl`

/* Programma per generare traiettorie del FPT della differenza tra la soluzione deterministica e quelle stocastiche,
	in funzione del tempo di autocorrelazione */


int main(){

  // COSTANTI GENERALI DEL PROBLEMA DA RISOLVERE

  const double k=0.03;

  // DEFINIZIONE DELLE PRINCIPALI VARIABILI USATE

  double T, sol, dt, dt_aux, tau, tau_aux, T_aux, sigma;
  double *X_cai, *Cai_noise, *X_tsb, *Tsb_noise, *X_sine, *Sine_noise, *X_white, *Noise_aux; 			/* size N+1*/
  double *time_cai, *time_tsb, *time_sw, *time_aux;
  double *FTP05cai, *FTP05tsb, *FTP05sine, *FTP05white, *FTP09cai, *FTP09tsb, *FTP09sine, *FTP09white;	/* size M */
  double *Time_val, *Sigma_val, *Cai05, *Tsb05, *Sine05, *White05, *Cai09, *Tsb09, *Sine09, *White09;	/* size R */
  double delta, q, B, x0, w0;
  int i, j, h, i1, i2, N_cai, N_tsb, N_aux, N, seed, seed1, M, R;
  FILE *fp;

  M=1500;		
  R=16;		

  T=75; 	
  dt=0.003;
  N=T/dt;

  B=2*k/3;
  delta=5;
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

  FTP05cai=(double*)malloc(M*sizeof(double));
  FTP05tsb=(double*)malloc(M*sizeof(double));
  FTP05sine=(double*)malloc(M*sizeof(double));
  FTP05white=(double*)malloc(M*sizeof(double));

  FTP09cai=(double*)malloc(M*sizeof(double));
  FTP09tsb=(double*)malloc(M*sizeof(double));
  FTP09sine=(double*)malloc(M*sizeof(double));
  FTP09white=(double*)malloc(M*sizeof(double));

  Cai05=(double*)malloc(R*sizeof(double));
  Tsb05=(double*)malloc(R*sizeof(double));
  Sine05=(double*)malloc(R*sizeof(double));
  White05=(double*)malloc(R*sizeof(double));

  Cai09=(double*)malloc(R*sizeof(double));
  Tsb09=(double*)malloc(R*sizeof(double));
  Sine09=(double*)malloc(R*sizeof(double));
  White09=(double*)malloc(R*sizeof(double));

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
			seed1= 5*(j+1);

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


			i1=indice(X_cai, 0.5, 0, N_cai);
			i2=indice(X_cai, 0.9, i1, N_cai);
			FTP05cai[j]=time_cai[i1];
			FTP09cai[j]=time_cai[i2];

			i1=indice(X_tsb, 0.5, 0, N_tsb);
			i2=indice(X_tsb, 0.9, i1, N_tsb);
			FTP05tsb[j]=time_tsb[i1];
			FTP09tsb[j]=time_tsb[i2];

			i1=indice(X_sine, 0.5, 0, N);
			i2=indice(X_sine, 0.9, i1, N);
			FTP05sine[j]=time_sw[i1];
			FTP09sine[j]=time_sw[i2];

			i1=indice_white(X_white, 0.5, 0, N);
			i2=indice_white(X_white, 0.9, i1, N);
			FTP05white[j]=time_sw[i1];
			FTP09white[j]=time_sw[i2];

		  } 	// fine for in j. A tau fissato (= i fissato), abbiamo calcolato M norme con quel tau

		Cai05[i]=mean(FTP05cai, M);
		Tsb05[i]=mean(FTP05tsb, M);
		Sine05[i]=mean(FTP05sine, M);		
		White05[i]=mean(FTP05white, M);

		Cai09[i]=mean(FTP09cai, M);
		Tsb09[i]=mean(FTP09tsb, M);
		Sine09[i]=mean(FTP09sine, M);		
		White09[i]=mean(FTP09white, M);

	}


  fp=fopen("..//File txt//FPT_vs_tau.txt", "w");

  if(fp!=NULL)
	for(i=0; i<R; i++) {
		fprintf(fp, "%g %g ", Time_val[i], Sigma_val[i]);
		fprintf(fp, "%g %g %g %g ", Cai05[i], Tsb05[i], Sine05[i], White05[i]);
		fprintf(fp, "%g %g %g %g\n", Cai09[i], Tsb09[i], Sine09[i], White09[i]);
			}
  else printf("Errore nell'apertura di FPT_vs_tau.txt");

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

  free(FTP05cai);
  free(FTP05tsb);
  free(FTP05sine);
  free(FTP05white);

  free(FTP09cai);
  free(FTP09tsb);
  free(FTP09sine);
  free(FTP09white);

  free(Cai05);
  free(Tsb05);
  free(Sine05);
  free(White05);

  free(Cai09);
  free(Tsb09);
  free(Sine09);
  free(White09);

  free(Time_val);
  free(Sigma_val);

return 0;

}
