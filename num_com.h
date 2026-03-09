#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Phi_X(x) log10(x+2)
//the function of Picard's method
#define F(x) x*exp(x)-1
//the function of the Newton's method



void picard(double x0,double epsilon,double (*phi_x)(double x));
void picard_aitken(double x0,double epsilon,double (*phi_x)(double x));
double derivative(double x,double epsilon,double (*f)(double x));
void newton(double x0,double epsilon,double (*f)(double x));
void newton_down(double x0,double epsilon,double (*f)(double x));
void single_secant(double x0,double epsilon,double (*f)(double x));
void double_secant(double x0,double epsilon,double (*f)(double x));
double dichotomy(double x1,double x2,double epsilon,double (*f)(double x));
void up_triangle(double U[n][n],double B[n],int n);
void down_triangle(double L[n][n],double B[n],int n);