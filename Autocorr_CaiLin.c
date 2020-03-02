#include <stdio.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include "noises.c"

// DA COMPILARE CON gcc autocorr_cailin.c `pkg-config --cflags --libs gsl`

// Stores values of an histogram (to be later plotted in Matlab) for the statonary density of 
// the Cai-Lin noise with the specified parameters.

int main() {

	double *time, *Noise;
	double *Autocorr, *time_autocorr;
	double delta, B, tau, dt, T, s, t;
	int N, R, M, seed, i, j, is, it;
	FILE *fp;
	
	dt=0.001;
	T=2.5;
	tau=1;
	B=1;
	delta=-0.5;
	
	N=T/dt;
	M=200000;
	R=20;

	Noise=(double*)malloc(2*N*sizeof(double));
	time=(double*)malloc(2*N*sizeof(double));
	Autocorr=(double*)malloc(R*sizeof(double));
	time_autocorr=(double*)malloc(R*sizeof(double));

	
	s=0.25;
	for(j=0; j<R; j++){
		time_autocorr[j]=s+ (j*0.1);
		Autocorr[j]=0;
	  }

	for(i=0; i<M; i++){
		seed=2*i+5;
		N=Cai_Lin(Noise, time, tau, B, delta, T, dt, 0, seed);
		is=indice(time, s, 0, N);		

		for(j=0; j<R; j++){
			t=time_autocorr[j];
			it=indice(time, t, 0, N);
			Autocorr[j] += Noise[is]*Noise[it];
	      }
	  }

	for(j=0; j<R; j++)
		Autocorr[j]/=M;

	fp=fopen("..//File txt//autocorr_cailin.txt", "w");
	for(j=0; j<R; j++){
		fprintf(fp, "%g %g\n", time_autocorr[j], Autocorr[j] );
	  }
	fclose(fp);

	free(Noise);
	free(time);
	free(Autocorr);
	free(time_autocorr);
}
