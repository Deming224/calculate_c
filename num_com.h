#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void picard(double x0,double epsilon,double (*phi_x)(double x));
void picard_aitken(double x0,double epsilon,double (*phi_x)(double x));
double derivative(double x,double epsilon,double (*f)(double x));
void newton(double x0,double epsilon,double (*f)(double x));
void newton_down(double x0,double epsilon,double (*f)(double x));
void single_secant(double x0,double epsilon,double (*f)(double x));
void double_secant(double x0,double x00,double epsilon,double (*f)(double x));
double dichotomy(double x1,double x2,double epsilon,double (*f)(double x));

void up_triangle(double **U, double *B, int n);
void down_triangle(double **L, double *B, int n);

void gauss_elimination(double **A,double *B,int n);
void gauss_elimination_pivot(double **A, double *B, int n);
void chasing_method(double *a, double *b, double *c, double *d, int n);