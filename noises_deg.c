#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include "noises_deg.h"

#define PI 3.141592654


double random_unif(double a, double b){
  double x;

  x = rand();
  x = x/RAND_MAX;
  x = a + (b-a)*x;
  return x;
}


int sign(double x) {
	if(x>0) return 1;
	else if (x<0) return -1;
	else return 0;
  }

int max(int a, int b){
	if(a>b)	return a;
	else return b;	
  }


double Mean_value(double A, double delta, double a){
  double G;
  G = pow(1+A,a) * gsl_sf_hyperg_2F1(-a, 1+delta, 2+ 2*delta, 2/(1+A));
  return G;
}

double phi(double x, double A, double a){
	double y;
	y=pow(A + x, a);
	return y;
	}

// Calcola la Zr corrispondente alla trasformazione
double A_Zr(double Zr, double delta, double a, double *M){
  double A, sx, dx;
  double Zr_temp, mean;

  Zr_temp=Zr +2;
  A=1.04;

  if(a<0){
	while(Zr_temp > Zr){
	  A+=0.00005;
	  sx=pow(A+1, a);
	  dx=pow(A-1,a);
	  mean=pow(1+A,a) * gsl_sf_hyperg_2F1(-a, 1+delta, 2+ 2*delta, 2/(1+A));
	  Zr_temp=(dx - mean)/(mean - sx);
	 }
   }
  else{
	while(Zr_temp > Zr){
	  A+=0.00005;
	  sx=pow(A-1, a);
	  dx=pow(A+1,a);
	  mean=pow(1+A,a) * gsl_sf_hyperg_2F1(-a, 1+delta, 2+ 2*delta, 2/(1+A));
	  Zr_temp=(dx - mean)/(mean - sx);
	 }
   }

  *M=mean;
  return A;
}

// Calcola la norma2 discreta della differenza tra la serie di misurazioni in U e la traiettoria DG di seme int e parametri P
// U ha 2*N_u elementi, da 0 a (2*N_u -1): posti pari tempi, posti dispari corrispondenti valori misurati.
// Ordine dei parametri in P: B, delta, tau, x0.

double Square_difference_dt(double *P, double *U, int N_u, double dt, double k1, double k2, double C0, int seed){

  double *Noise, *C, *Q;
  double T, x0, delta, tau, B, t, sum, y;
  int N, i, ind;

  T=U[N_u -1] +1;
  N=T/dt;
  Noise=(double*)malloc((2*N +1)*sizeof(double));
  C=(double*)malloc((N+1)*sizeof(double));
  Q=(double*)malloc((N+1)*sizeof(double));

  B=P[0];
  delta=P[1];
  tau=P[2];
  x0=P[3];

  Cai_Lin(Noise, tau, B, delta, T, 2*N, x0, seed);
  DeGaetano(C, Q, Noise, dt, k1, k2, C0, N);

  sum=0;
  for(i=0; i<N_u; i++){
	t=U[i];
	ind=t/dt;
	y= C[ind] - U[N_u + i];
//	y= Q[ind] - U[N_u + i];
	sum+= y*y;
  }
  sum=sqrt(sum);

  free(Noise);
  free(C);
  free(Q);

 return sum;
}

double Energy_min_dt(double *P, double *U, int N_u, int K, double dt, double k1, double k2, double C0){
  
  int i, seed_temp;
  double min, valore;

  seed_temp=rand();
  min=Square_difference_dt(P, U, N_u, dt, k1, k2, C0, seed_temp);

  for(i=0; i<K-1; i++){
	seed_temp=rand();
	valore=Square_difference_dt(P, U, N_u, dt, k1, k2, C0, seed_temp);
	if(valore<min)
		min=valore;
  }

  return min;
}


double Energy_mean_dt(double *P, double *U, int N_u, int K, double dt, double k1, double k2, double C0){
  
  int i, seed_temp;
  double mean, valore;

  mean=0;

  for(i=0; i<K; i++){
	seed_temp=rand();
	valore=Square_difference_dt(P, U, N_u, dt, k1, k2, C0, seed_temp);
	mean+= valore;
  }

  mean=mean/K;
  return mean;
}



Cooling(double *T, double R){
	*T = (*T)*R;	
}

// P: B, delta, tau, x0.   Step, domain: B, delta, tau.
void P_generation(double *P_old, double *P_new, double *step, double *domain){

  double min, max, B;
  int i;

  for(i=0; i<3; i++){
		P_new[i] = P_old[i] + step[i]*random_unif(-1,1);
		min=domain[i];
		max=domain[i+3];
		if(P_new[i]<min)
			P_new[i]=min;
		else if(P_new[i]>max)
			P_new[i]=max;
	 }	
  B=P_new[0];
  P_new[3]=P_old[3] + 0.01*B*random_unif(-1,1);

  if(P_new[3]<-B)
	P_new[3]= -0.999*B;
  else if(P_new[3]>B)
	P_new[3]=0.999*B;

}

// Ordine informazioni nei vettori P, P_fin, step_par, domain_par: B, delta, tau, (x0). 
void SA(double *U, int N_u, double *P_iniz, int K, double T0, double Tf, double *step_par, double *domain_par, double *Fixed_par, double *P_fin, double *E_fin){

  double Temp, min, max, k1, k2, C0;
  double E, E_new, Delta, s, dt;
  double *P, *P_new;
  int i, seed;

  P=(double*)malloc(4*sizeof(double));
  P_new=(double*)malloc(4*sizeof(double));

  k1=Fixed_par[0];
  k2=Fixed_par[1];
  C0=Fixed_par[2];

  dt=0.005;

// Generazione del P inziale
  srand(5);  						// DA CAMBIARE

  for(i=0; i<3; i++){
	min = domain_par[i];
	max = domain_par[i+3];
	P[i] = random_unif(min,max);
  }
  P[3]=random_unif(-P[0], P[0]);

  for(i=0;i<4; i++)
	P_iniz[i]=P[i];

  E=Energy_min_dt(P, U, N_u, K, dt, k1, k2, C0);
printf("E=%g\n", E);

  Temp=T0;

  while(Temp>Tf){

// Generazione di P_new vicino a P, e calcolo di E_new
	P_generation(P, P_new, step_par, domain_par);
	E_new=Energy_min_dt(P_new, U, N_u, K, dt, k1, k2, C0);
	Delta=E_new - E; 
//printf("Delta=%g, T=%g\n", Delta, Temp);

// Accettazione di P_new nei due casi

	if(Delta<0){
		for(i=0; i<4;i++){
			P[i]=P_new[i];
		 }
		E=E_new; 
		printf("E=%g, T=%g\n", E, Temp);
	 }
	else{
		s=random_unif(0,1);
		if(exp(-Delta/Temp)>s){
			for(i=0; i<4;i++){
				P[i]=P_new[i];
			 }
		E=E_new; 
		printf("E=%g, T=%g\n", E, Temp);
		 }
	 }
	if(Temp>2)
	Cooling(&Temp, 0.99);
	else if(Temp>1)
	Cooling(&Temp, 0.995);
	else
	Cooling(&Temp, 0.995);

  }//end while

  for(i=0; i<4; i++)
	P_fin[i]=P[i];

//  *(E_fin) =E;

  free(P);
  free(P_new);
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


void Sine_Wien(double *X, double tau, double B, double T, int N, double W0, int seed){

  // DICHIARAZIONE DELLE VARIABILI
  
  double dt, x, gamma, c, dW;
  int i;
  
  const gsl_rng_type * R;
  gsl_rng * r;
  gsl_rng_env_setup ();
  R = gsl_rng_default;
  r = gsl_rng_alloc (R);
  gsl_rng_set(r,seed); 

  // CALCOLO DEL VETTORE

  dt=T/N;
  
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


void DeGaetano(double *C, double *Q, double *Noise, double dt, double k1, double k2, double C0, int N){

//C, Q e time vanno da 0 a N. Noise va da 0 a 2N.
	double h, No1, No2, No3;
	double H1c, H1q, H2c, H2q, H3c, H3q, H4c, H4q, Z2c, Z2q, Z3c, Z3q, Z4c, Z4q; 
	int i;

	C[0]=C0;
	Q[0]=0;
	h=dt;

	for(i=0; i<N; i++){
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

void DeGaetano_white(double *C, double *Q, double *k, double k1, double k2, double C0, double sigma, double h, int N, int seed){
	
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
		C[i+1]= C[i] + h*H1c;
		dW=gsl_ran_gaussian_ziggurat(r,1);
		Q[i+1]= Q[i] + h*H1q + sigma*Q[i]*c*dW + 0.5*sigma*sigma*Q[i]*h*(dW*dW -1);
		k[i] = k2 + sigma*c*dW/h + 0.5*sigma*sigma*(dW*dW -1);

/*
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
*/
	}
k[N]=0;
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

int indice_bisez_cresc(double *X, double t, int sx, int dx){
	int cx;
	if(sx==dx) return sx;
	else{
		cx=(sx+dx)/2;
		if(t<=X[cx]) return indice_bisez_cresc(X,t,sx,cx);
		else return indice_bisez_cresc(X,t,cx+1,dx);
	}
}



int indice_bisez_decr(double *X, double t, int sx, int dx){
	int cx;
	if(sx==dx) return sx;
	else{
		cx=(sx+dx)/2;
		if(t<=X[cx]) return indice_bisez_decr(X,t,cx+1,dx);
		else return indice_bisez_decr(X,t,sx,cx);
	}
}


int indice_cresc(double *X, double t, int sx, int dx){
	int ind=sx;
	while(ind<dx && X[ind]<t)
		ind++;
	return ind;
  }

int indice_decr(double *X, double t, int sx, int dx){
	int i;
	i=sx;
	while(i<dx && X[i]>t)
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












double Square_difference_asymm(double *P, double A, double mean, double *U, int N_u, double dt, double k1, double k2, double C0, int seed){

  double *Noise, *C, *Q;
  double T, B, Zr, delta, a, tau, x0, t, sum, y;
  double sx, x;
  int N, i, ind;

  T=U[N_u -1] +1;
  N=T/dt;
  Noise=(double*)malloc((2*N +1)*sizeof(double));
  C=(double*)malloc((N+1)*sizeof(double));
  Q=(double*)malloc((N+1)*sizeof(double));

  B=1;

  Zr=P[0];
  delta=P[1];
  a=P[2];
  tau=P[3];
  x0=P[4];

  if(a<0)
	sx=pow(A+1, a);
  else
	sx=pow(A-1, a);

  Cai_Lin(Noise, tau, B, delta, T, 2*N, x0, seed);
  for(i=0; i<(2*N +1); i++){
	x=phi(Noise[i], A, a);
	Noise[i]=k2 * (x - mean)/(mean - sx);	
   }	
  DeGaetano(C, Q, Noise, dt, k1, k2, C0, N);

  sum=0;
  for(i=0; i<N_u; i++){
	t=U[i];
	ind=t/dt;
	y= C[ind] - U[N_u + i];
//	y= Q[ind] - U[N_u + i];
	sum+= y*y;
  }
  sum=sqrt(sum);

  free(Noise);
  free(C);
  free(Q);

 return sum;
}

// Elementi in P: Zr, delta, a, tau, x0
double Energy_min_asymm(double *P, double *U, int N_u, int K, double dt, double k1, double k2, double C0, FILE *fp){
  
  int i, seed;
  double A, delta, a, Zr, mean;
  double min, valore;

  Zr=P[0];
  delta=P[1];
  a=P[2];

  A=A_Zr(Zr, delta, a, &mean);

  seed=rand();
  min=Square_difference_asymm(P, A, mean, U, N_u, dt, k1, k2, C0, seed);

  for(i=0; i<K-1; i++){
	seed=rand();
	valore=Square_difference_asymm(P, A, mean, U, N_u, dt, k1, k2, C0, seed);
	if(valore<min)
		min=valore;
  }

  return min;
}

// Elementi in P: Zr, delta, a, tau, x0
double Energy_mean_asymm(double *P, double *U, int N_u, int K, double dt, double k1, double k2, double C0, FILE *fp){
  
  int i, seed;
  double A, delta, a, Zr, Tot;
  double mean, valore;

  Zr=P[0];
  delta=P[1];
  a=P[2];

  A=A_Zr(Zr, delta, a, &mean);

  Tot=0;

  for(i=0; i<K; i++){
	seed=rand();
	valore=Square_difference_asymm(P, A, mean, U, N_u, dt, k1, k2, C0, seed);
	Tot += valore;
  }

  Tot=Tot/K;
  return Tot;
}


// Elementi in P, Step, domain: Zr, delta, a, tau, x0.
void P_generation_asymm(double *P_old, double *P_new, double *step, double *domain, FILE *fp){

  double min, max, u;
  int i;

  for(i=0; i<5; i++){
		u=random_unif(-1,1);
		if(i==3){
			if(P_old[i]>30)
			  P_new[i] = P_old[i] + step[i]*u;
			else if(P_old[i]>15)
			  P_new[i] = P_old[i] + 0.7*step[i]*u;
			else if (P_old[i]>6)
			  P_new[i] = P_old[i] + 0.25*step[i]*u;
			else 
			  P_new[i] = P_old[i] + 0.1*step[i]*u;

			min=domain[i];
			max=domain[i+5];
			if(P_new[i]<min)
				P_new[i]=min;
			else if(P_new[i]>max)
				P_new[i]=max;
		  }

		else{
			P_new[i] = P_old[i] + step[i]*u;
			min=domain[i];
			max=domain[i+5];
			if(P_new[i]<min)
				P_new[i]=min;
			else if(P_new[i]>max)
				P_new[i]=max;
		  }
	 }	

}

// Ordine informazioni nei vettori P_iniz, P_fin, step_par, domain_par: Zr, delta, a, tau, x0. 
void SA_asymm(double *U, int N_u, double *P_iniz, int K, double T0, double Tf, double *step_par, double *domain_par, double *Fixed_par, double *P_fin, double *E_fin, double R, int N_step){

  double Temp, min, max, k1, k2, C0, B;
  double E, E_new, Delta, s, dt;
  double *P, *P_new;
  int i, seed, current_step;
  FILE *fp;

  P=(double*)malloc(5*sizeof(double));
  P_new=(double*)malloc(5*sizeof(double));

  k1=Fixed_par[0];
  k2=Fixed_par[1];
  C0=Fixed_par[2];

  dt=0.01;

//  fp=fopen("..//..//File txt/random_unif.txt", "w");

  srand(time(NULL)); 

// Generazione del P iniziale

  for(i=0; i<5; i++){
	min = domain_par[i];
	max = domain_par[i+5];
	P[i] = random_unif(min,max);
  }

  for(i=0;i<5; i++){
	P_iniz[i] = P[i];
	P_fin[i] = P[i];
	}

  E=Energy_min_asymm(P, U, N_u, K, dt, k1, k2, C0, fp);
  *(E_fin)=E;

  Temp=T0;
  current_step=0;

  while(Temp>Tf){

// Generazione di P_new vicino a P, e calcolo di E_new
	P_generation_asymm(P, P_new, step_par, domain_par, fp);
	E_new=Energy_min_asymm(P_new, U, N_u, K, dt, k1, k2, C0, fp);
	Delta=E_new - E; 

// Accettazione di P_new nei due casi

	if(Delta<0){
		for(i=0; i<5;i++){
			P[i]=P_new[i];
		 }
		E=E_new; 
//		printf("E=%g, T=%g, Zr=%g, d=%g, a=%g, tau=%g, x0=%g\n", E, Temp, P[0], P[1], P[2], P[3], P[4]);
		if(E< *(E_fin)){
			*(E_fin)=E;
			for(i=0; i<5; i++)
				P_fin[i]=P[i];
			  }
	 }

	else{
		s=random_unif(0,1);
		if(exp(-Delta/Temp)>s){
			for(i=0; i<5; i++)
				P[i]=P_new[i];
		E=E_new;
		}
	 }

	current_step++;
	if(current_step==N_step){
		Cooling(&Temp, R);
		current_step=0;
	  }

  }//end while

//  fclose(fp);


  free(P);
  free(P_new);
}


