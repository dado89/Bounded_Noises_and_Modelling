#include <stdio.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>

int sign(double x);
int max(int a, int b);
double norm2(double *X, double *time, int N);
double norminf(double *X, int N);
double mean(double *X, int M);



double phi1(double x, double alpha);
double phi2(double x, double a1, double a2);
double alpha_Zr1(double delta, double Zr, double *M);
double alpha2_Zr2(double delta, double a1, double Zr, double *M);

double Square_difference_dt(double *P, double *U, int N_u, double dt, double k1, double k2, double C0, int seed);
double Energy_min_dt(double *P, double *U, int N_u, int K, double dt, double k1, double k2, double C0);
double Energy_mean_dt(double *P, double *U, int N_u, int K, double dt, double k1, double k2, double C0);



void Brown1(double dt, double *B, int seed);
double Half_bridge(double b, double dt, int seed);
double F_cai(double x, double dt, double B, double mu, double gamma, double gammaquad, double w, int seed);
double F_tsb(double x, double dt, double B, double mu, double gamma, double w, int seed);

void Cai_Lin(double *X, double tau, double B, double delta, double T, int N, double x0, int seed);
void Tsall_Borl(double *X, double tau, double B, double q, double T, int N, double x0, int seed);
void Sine_Wien(double *X, double tau, double B, double dt, int N, double W0, int seed);



void F_Deg(double *X, double *Y, double c, double q, double k1, double k2);
void DeGaetano(double *C, double *Q, double *Noise, double dt, double k1, double k2, double C0, int N);
void DeGaetano_white(double *C, double *Q, double *k, double k1, double k2, double C0, double sigma, double h, int N, int seed);

void Lansky_white(double *X, double k, double sigma, double T, int N, int seed);
void Lansky_white_k(double *X, double *k_val, double k, double sigma, double T, int N, int seed);
void Lansky(double *X, double *Noise, double *time, double k, int N);



int indice(double *X, double t, int sx, int dx);
int indice_cresc(double *X, double t, int sx, int dx);
int indice_decr(double *X, double t, int sx, int dx);



double F_Nonlin1(double y, double J, double nu_t);
void Nonlin1(double *Y, double *Noise, double J, double y0, int N, double dt);
void Nonlin1_white(double *Y, double J, double sigma, double y0, int N, double dt, int seed);


int Cai_Lin_fals(double *X, double *time, double tau, double B, double delta, double T, double dt, double x0, int seed);
int Tsall_Borl_fals(double *X, double *time, double tau, double B, double q, double T, double dt, double x0, int seed);


