#include "num_com.h"

void picard_aitken(double x0,double epsilon)
{
    double x,x1,x2;
    x0 = 0;//the initial value
    x = x0;
    for(int k=1;;k++)
    {
        x1 = x;
        x2 = ((x1*(Phi_X(Phi_X(x1))))-((Phi_X(x1))*(Phi_X(x1))))/((Phi_X(Phi_X(x1)))-(2*(Phi_X(x1)))+x1);//the formula of Aitken's method
        if(fabs(x2-x1)<=epsilon)
        {
            printf("%lf\n",x2);
            printf("k=%d\n",k);
            break;
        }
        else
        {
            x = x2;
            printf("%lf\n",x);
        }
    }
    printf("calulation is done");
    system("pause");
}

void picard(double x0,double epsilon)
{
     double x,x1,x2;
    x0 = 0;//the initial value
    x = x0;
    for(int k=1;;k++)
    {
        x1 = x;
        x2 = Phi_X(x);//the formula of Picard's method
        if(fabs(x2-x1)<=epsilon)
        {
            printf("%lf\n",x2);
            printf("k=%d\n",k);
            break;
        }
        else
        {
            x = x2;
            printf("%lf\n",x);
        }
    }
    printf("calulation is done");
    system("pause");
}

double derivative(double x,double epsilon)
{
    double dx = 0.0001;
    double dy_0,dy_1;
    do
    {
        dy_0 = ((F(x+dx)) - (F(x)))/dx;
        dx = 0.5*dx;
        dy_1 = ((F(x+dx))-(F(x)))/dx;
    }
    while(fabs(dy_0-dy_1)>epsilon);
    return dy_1;
}

void Newton(double x0,double epsilon)
{
    double k,x,x1,x2;
    x = x0;
    for(int i=1;;i++)
    {
        x1 = x;
        k = derivative(x1,epsilon);//epsilon needs to be changed as your wish
        x2 = x1-((F(x1))/k);//the formula of Newton's method
        if(fabs(x2-x1)<=epsilon)
        {
            printf("%lf\n",x2);
            printf("k=%d\n",i);
            break;
        }
        else
        {
            x = x2;
            printf("%lf\n",x);
        }
    }
    printf("calculation is done");
    system("pause");
}

double dichotomy(double x1,double x2,double epsilon)
{
    double x0;
    for(int k;;k++)
    {
        if((F(x1))*(F(x2))<0)
        {
            x0 = (x1+x2)/2;
            if((x2-x1)<=epsilon)
            {
                break;
            }
            else
            {
                if((F(x1)*(F(x0)))<0)
                {
                    x2 = x0;
                }
                else if((F(x0))*(F(x2))<0)
                {
                    x1 = x0;
                }
            }       
        }
        else
        {
            x0 = 0;
            break;
        }
    }
    return x0;
}