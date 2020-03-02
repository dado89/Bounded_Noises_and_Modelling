#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
//#include "noises.h"



int sign(double x) {
	if(x>0) return 1;
	else if (x<0) return -1;
	else return 0;
  }

int max(int a, int b){
	if(a>b)	return a;
	else return b;	
  }




// Restituisce N=ultimo indice del vettore X  time. I vettori vanno quindi da 0 a N.

int Cai_Lin_fals(double *X, double *time, double tau, double B, double delta, double T, double dt, double x0, int seed){

  double dW, dt1;
  double x, gammaquad, gamma, drift, diff, mu, c;
  int i;

  const gsl_rng_type * R;
  gsl_rng * r;
  gsl_rng_env_setup ();
  R = gsl_rng_default;
  r = gsl_rng_alloc (R);
  gsl_rng_set(r,seed); 

  if(delta<=-1) {
	printf("Valore di delta non valido. Inserire un delta maggiore di -1\n");
	return -1;}

  // DICHIARAZIONE DELLE VARIABILI
  
  mu=1/tau;
  gammaquad=mu/(1+delta);
  gamma=sqrt( gammaquad );
    
  X[0]=x0;
  time[0]=0;
  i=0;
  
  c=sqrt(dt);
  while(time[i]<T){
	x = X[i];
	drift = -mu*x;
	dW=gsl_ran_gaussian_ziggurat(r,1);
	diff=sqrt(B*B - x*x);
	X[i+1] = X[i] + drift*dt + gamma*diff*c*dW - 0.5*gammaquad*x*dt*(dW*dW -1);
	dt1=dt;
	while(fabs(X[i+1]) > B){
		dt1 = 0.1*dt1;
		X[i+1] = X[i] + drift*dt1 + gamma*diff*sqrt(dt1)*dW - 0.5*gammaquad*x*dt1*(dW*dW -1);
		}
	time[i+1] = time[i] + dt1;
	i=i+1;
   }

  return i;
  gsl_rng_free(r);

}

void Brown1(double dt, double *B, int seed){
  
  int i, N;
  double c, dW;

  const gsl_rng_type * R;
  gsl_rng * r;
  gsl_rng_env_setup ();
  R = gsl_rng_default;
  r = gsl_rng_alloc (R);
  gsl_rng_set(r,seed); 

  N=1/dt;
  c=sqrt(dt);
  B[0]=0;
  for(i=0; i<N; i++){
	dW=gsl_ran_gaussian_ziggurat(r,c);
	B[i+1] = B[i] + dW;
	}
  gsl_rng_free(r);
}

/* Se un MB al tempo t0 vale a e al tempo t1 vale b, il seguente programma restituisce un possibile
valore del MB al tempo (t0+t1)/2 */

double Half_bridge(double b, double dt, int seed){
  double dt_temp, x, c;
  double *B;
  int N;

  c=sqrt(dt);
  dt_temp=0.5;
  N=1/dt_temp;
  B=(double*)malloc((N+1)*sizeof(double));
  Brown1(dt_temp, B, seed);
  x = 0.5*b + c*(B[1] - 0.5*B[2]);
  free(B);
  return x;
}


double F_cai(double x, double dt, double B, double mu, double gamma, double gammaquad, double w, int seed){
  double drift, diff, w1;
  double y, y1;

  drift = -mu*x;
  diff=sqrt(B*B - x*x);
  
  y = x + drift*dt + gamma*diff*w - 0.5*gammaquad*x*(w*w - dt);
  while(fabs(y) > B){
	w1= Half_bridge(w, dt, 2*seed +5);
	y1= F_cai(x,  dt/2, B, mu, gamma, gammaquad, w1, 2*seed + 7);
	y = F_cai(y1, dt/2, B, mu, gamma, gammaquad, w-w1, 2*seed + 19);
   }
  return y;
}
  


void Cai_Lin(double *X, double tau, double B, double delta, double T, int N, double x0, int seed){

  double dW;
  double dt, gammaquad, gamma, mu, c;
  int i;

  const gsl_rng_type * R;
  gsl_rng * r;
  gsl_rng_env_setup ();
  R = gsl_rng_default;
  r = gsl_rng_alloc (R);
  gsl_rng_set(r,seed); 

  if(delta<=-1) {
	printf("Valore di delta non valido. Inserire un delta maggiore di -1\n");
	exit(1); }

  dt=T/N;
    
  mu=1/tau;
  gammaquad=mu/(1+delta);
  gamma=sqrt( gammaquad );
  c=sqrt(dt);

  X[0]=x0;

  for(i=0; i<N; i++){
	dW=gsl_ran_gaussian_ziggurat(r,c);
	X[i+1] = F_cai(X[i], dt, B, mu, gamma, gammaquad, dW, seed+i);
  }

  gsl_rng_free(r);
}



int Tsall_Borl_fals(double *X, double *time, double tau, double B, double q, double T, double dt, double x0, int seed){
  
  double dW, theta, dt1;
  double x, gamma, drift, c, mu;
  int i;

  const gsl_rng_type * R;
  gsl_rng * r;
  gsl_rng_env_setup ();
  R = gsl_rng_default;
  r = gsl_rng_alloc (R);
  gsl_rng_set(r,seed); 

  if(q>=1) {
	printf("Valore di q non valido. Inserire un q minore di 1\n");
	return -1;}

  theta=(2.5 - 1.5*q) * tau;
  mu=1/theta;
  gamma=B*sqrt((1-q)*mu);
    
  // CALCOLO DELLA SOLUZIONE
  
  X[0]=x0;
  time[0]=0;
  i=0;

  c=sqrt(dt);
  while(time[i]<T){
	x = X[i];
	drift = mu * (B*B*x)/(B*B - x*x);
	dW=gsl_ran_gaussian_ziggurat(r,1);
	X[i+1] = X[i] - drift*dt + gamma*c*dW;
	dt1=dt;
	while(fabs(X[i+1]) >= B){
		dt1 = 0.1*dt1;
		X[i+1] = X[i] - drift*dt1 + gamma*sqrt(dt1)*dW;
		}
	time[i+1] = time[i] + dt1;
	i=i+1;
   }

  return i;
  gsl_rng_free(r);

}


double F_tsb(double x, double dt, double B, double mu, double gamma, double w, int seed){
  double drift, diff, w1, y, y1;

  drift = -mu * (B*B*x)/(B*B - x*x);
  diff=gamma*w;
  
  y = x + drift*dt + diff;
  while(fabs(y) > B){
	w1= Half_bridge(w, dt, 2*seed +5);
	y1= F_tsb(x,  dt/2, B, mu, gamma, w1, 2*seed + 7);
	y = F_tsb(y1, dt/2, B, mu, gamma, w-w1, 2*seed + 19);
   }
  return y;
}

void Tsall_Borl(double *X, double tau, double B, double q, double T, int N, double x0, int seed){
  
  double dW, theta;
  double dt, gamma, c, mu;
  int i;

  const gsl_rng_type * R;
  gsl_rng * r;
  gsl_rng_env_setup ();
  R = gsl_rng_default;
  r = gsl_rng_alloc (R);
  gsl_rng_set(r,seed); 

  if(q>=1) {
	printf("Valore di q non valido. Inserire un q minore di 1\n");
	exit(1); }

  dt=T/N;

  theta=(2.5 - 1.5*q) * tau;
  mu=1/theta;
  gamma=B*sqrt((1-q)*mu);
  c=sqrt(dt);
  
  X[0]=x0;

  for(i=0; i<N; i++){
	dW=gsl_ran_gaussian_ziggurat(r,c);
	X[i+1] = F_tsb(X[i], dt, B, mu, gamma, dW, seed+i);
   }

  gsl_rng_free(r);
}


void Sine_Wien(double *X, double tau, double B, double dt, int N, double W0, int seed){

  // DICHIARAZIONE DELLE VARIABILI
  
  double x, gamma, c, dW;
  int i;
  
  const gsl_rng_type * R;
  gsl_rng * r;
  gsl_rng_env_setup ();
  R = gsl_rng_default;
  r = gsl_rng_alloc (R);
  gsl_rng_set(r,seed); 

  // CALCOLO DEL VETTORE
  
  x=W0;
  gamma=sqrt(2/tau);
  c=sqrt(dt);
  X[0]=B*sin(gamma*x);
  for(i=0; i<N; i++){
	dW=gsl_ran_gaussian_ziggurat(r,c);
	X[i+1]=B*sin(gamma*(x+dW));
	x+=dW;
	}
  
  gsl_rng_free(r);

}


void F_Deg(double *X, double *Y, double c, double q, double k1, double k2){
	*X= -k1*c + k1*q;
	*Y= k1*c - (k1+k2)*q;
}	


void DeGaetano(double *C, double *Q, double *Noise, double *time, double k1, double k2, double C0, int N){

//C, Q e time vanno da 0 a N. Noise va da 0 a 2N.
	double h, No1, No2, No3;
	double H1c, H1q, H2c, H2q, H3c, H3q, H4c, H4q, Z2c, Z2q, Z3c, Z3q, Z4c, Z4q; 
	int i;

	C[0]=C0;
	Q[0]=0;

	for(i=0; i<N; i++){
		h=time[i+1] - time[i];
		No1=Noise[2*i];
		No2=Noise[2*i + 1];
		No3=Noise[2*(i+1)];
		F_Deg(&H1c, &H1q, C[i], Q[i], k1, k2+No1);
		Z2c=C[i] + 0.5*h*H1c;
		Z2q=Q[i] + 0.5*h*H1q;
		F_Deg(&H2c, &H2q, Z2c, Z2q, k1, k2+ No2);
		Z3c=C[i] + 0.5*h*H2c;
		Z3q=Q[i] + 0.5*h*H2q;		
		F_Deg(&H3c, &H3q, Z3c, Z3q, k1, k2+No2);
		Z4c=C[i] + h*H3c;
		Z4q=Q[i] + h*H3q;	
		F_Deg(&H4c, &H4q, Z4c, Z4q, k1, k2+No3);

		C[i+1]= C[i] + h/6*(H1c + 2*H2c + 2*H3c + H4c);
		Q[i+1]= Q[i] + h/6*(H1q + 2*H2q + 2*H3q + H4q);
	}

}

void DeGaetano_white(double *C, double *Q, double k1, double k2, double C0, double sigma, double h, int N, int seed){
	
	double H1c, H1q, H2c, H2q, H3c, H3q, H4c, H4q, Z2c, Z2q, Z3c, Z3q, Z4c, Z4q, Z1c, Z1q, c, dW; 
	int i;	

	const gsl_rng_type * R;
	gsl_rng * r;
	gsl_rng_env_setup ();
	R = gsl_rng_default;
	r = gsl_rng_alloc (R);
	gsl_rng_set(r,seed); 

	c=sqrt(h);

	C[0]=C0;
	Q[0]=0;

	for(i=0; i<N; i++){
		F_Deg(&H1c, &H1q, C[i], Q[i], k1, k2);
		Z2c=C[i] + 0.5*h*H1c;
		Z2q=Q[i] + 0.5*h*H1q;
		F_Deg(&H2c, &H2q, Z2c, Z2q, k1, k2);
		Z3c=C[i] + 0.5*h*H2c;
		Z3q=Q[i] + 0.5*h*H2q;		
		F_Deg(&H3c, &H3q, Z3c, Z3q, k1, k2);
		Z4c=C[i] + h*H3c;
		Z4q=Q[i] + h*H3q;	
		F_Deg(&H4c, &H4q, Z4c, Z4q, k1, k2);

		C[i+1]= C[i] + h/6*(H1c + 2*H2c + 2*H3c + H4c);
		F_Deg(&Z1c, &Z1q, C[i]+ h*H1c, Q[i] + h*H1q,  k1, k2);
		dW=gsl_ran_gaussian_ziggurat(r,1);
		Q[i+1]= Q[i] + h/2*(H1q + Z1q) + sigma*Q[i]*c*dW + 0.5*sigma*sigma*Q[i]*h*(dW*dW -1);
	}
}



double F_Nonlin1(double y, double J, double nu_t){
	double output;
	output = J - (1 + nu_t)*y/(1+y);
	return output;
	}	


void Nonlin1(double *Y, double *Noise, double J, double y0, int N, double dt){

// Y va da 0 a N. Noise va da 0 a 2N.

  double y, No1, No2, No3;
  double H1, H2, H3, H4, Z1, Z2, Z3; 
  int i;

  Y[0]=y0;

  for(i=0; i<N; i++){
	y=Y[i];
	No1=Noise[2*i];
	No2=Noise[2*i + 1];
	No3=Noise[2*(i+1)];
	H1=F_Nonlin1(y, J, No1);
	Z1=y+ (0.5*dt*H1);
	H2=F_Nonlin1(Z1, J, No2);
	Z2=y + (0.5*dt*H2);	
	H3=F_Nonlin1(Z2, J, No2);
	Z3=y + (dt*H3);
	H4=F_Nonlin1(Z3, J, No3);

	Y[i+1]= y + (dt/6)*(H1 + 2*H2 + 2*H3 + H4);
  }
}



void Nonlin1_white(double *Y, double J, double sigma, double y0, int N, double dt, int seed){

// Y va da 0 a N. Noise va da 0 a 2N.

  double y, c, drift, diff, z, c1, c2, c3, dW;
  int i;

  const gsl_rng_type * R;
  gsl_rng * r;
  gsl_rng_env_setup ();
  R = gsl_rng_default;
  r = gsl_rng_alloc (R);
  gsl_rng_set(r,seed); 

  Y[0]=y0;
  c=sqrt(dt);

  for(i=0; i<N; i++){
	y=Y[i];
	c1=1+y;
	c2=c1*c1;
	c3=c1*c2;
	dW=gsl_ran_gaussian_ziggurat(r, c);
	drift= J - y/c1;
	diff=sigma*y/c1;
	z=sigma*sigma*y/c3;
	
	Y[i+1]= y + drift*dt + diff*dW + 0.5*z*(dW*dW - dt);
  }

  gsl_rng_free(r);
}


/* Risolve il problema dX = (k+noise)(1-X)dt con schema di Heun */

void Lansky(double *X, double *Noise, double *time, double k, int N){

	int i;
	double x, x_temp, dt;

	X[0]=0;

	for(i=0; i<N; i++){
		x=X[i];
		dt=time[i+1]-time[i];
		x_temp = x + dt*(1-x)*(k + Noise[i]);		
		X[i+1] = x + 0.5*dt* ( (1-x)*(k + Noise[i]) + (1-x_temp)*(k + Noise[i+1]) );
	  }
  }



/* Calcola una approssimazione della traiettoria di
d(phi) = K*(1-phi)*dt + sigma*(1-phi)*dW.
*/

void Lansky_white(double *X, double k, double sigma, double T, int N, int seed){

  // DICHIARAZIONE DELLE VARIABILI
  double x, dW, dt, c;
  int i;

  const gsl_rng_type * R;
  gsl_rng * r;
  gsl_rng_env_setup ();
  R = gsl_rng_default;
  r = gsl_rng_alloc (R);
  gsl_rng_set(r,seed); 

  dt=T/N;
  c=sqrt(dt);

  x=0;
  X[0]=0;
  for(i=0; i<N; i++){
	dW=gsl_ran_gaussian_ziggurat(r,c);
	x = x + k*(1-x)*dt + sigma*(1-x)*dW - 0.5*sigma*sigma*(1-x)*(dW*dW - dt);
	X[i+1]=x;
	}

  gsl_rng_free(r);
}


/* Come prima, ma restituisce anche il vettore di N componenti con i valori di (k*dt + sigma*dW[i]) */

void Lansky_white_k(double *X, double *k_val, double k, double sigma, double T, int N, int seed){

  // DICHIARAZIONE DELLE VARIABILI
  double x, c, dW, dt;
  int i;

  const gsl_rng_type * R;
  gsl_rng * r;
  gsl_rng_env_setup ();
  R = gsl_rng_default;
  r = gsl_rng_alloc (R);
  gsl_rng_set(r,seed); 

  // INIZIALIZZAZIONE DEI DATI E ALLOCAZIONE DELLA MEMORIA

  dt=T/N;
  c=sqrt(dt);
  x=0;
  X[0]=0;
  for(i=0; i<N; i++){
  	dW=gsl_ran_gaussian_ziggurat(r,c);
	x = x + k*(1-x)*dt + sigma*(1-x)*dW - 0.5*sigma*sigma*(1-x)*(dW*dW - dt);
	X[i+1]=x;
	k_val[i]=k*dt + sigma*dW;
	}

  gsl_rng_free(r);
}




/* TROVA IL PRIMO INDICE DEL VETTORE X IL CUI VALORE E' MAGGIORE O UGUALE A t€[0,1] */

int indice(double *X, double t, int sx, int dx){
	
	int cx;

	if(sx==dx) return sx;
	else{
		cx=(sx+dx)/2;
		if(t<=X[cx]) return indice(X,t,sx,cx);
		else return indice(X,t,cx+1,dx);
	}

}




int indice_white(double *X, double t, int sx, int dx){
	int ind=sx;
	while(ind<dx && X[ind]<t)
		ind++;
	return ind;
  }

int indice_deg_seq(double *X, double t, int sx, int dx){
	int i;
	i=0;
	while(X[i]>t && i<dx)
		i++;
	return i;
  }


double norm2(double *X, double *time, int N){
	double dt, sum=0;
	int i;
	for(i=0; i<N; i++){
		dt=time[i+1]-time[i];
		sum += X[i]*X[i];
	  }
	return sqrt(sum*dt);
  }


double norminf(double *X, int N){
	double max=0;
	int i;
	for(i=0; i<N+1; i++){
		if(fabs(X[i])>max)
			max=fabs(X[i]);
	  }
	return max;
  }


double mean(double *X, int M){
	double mean=0;
	int i;
	for(i=0; i<M; i++)
		mean+=X[i];
	  
	mean=mean/M;
	return mean;
  }

