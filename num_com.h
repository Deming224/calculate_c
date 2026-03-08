#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Phi_X(x) pow(x+1,1.0/3.0)
//the function of Picard's method
#define F(x) x*exp(x)-1
//the function of the Newton's method

void picard(double x0,double epsilon);
void picard_aitken(double x0,double epsilon);
double derivative(double x,double epsilon);
void Newton(double x0,double epsilon);
double dichotomy(double x1,double x2,double epsilon);