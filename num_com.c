#include "num_com.h"

//埃德金加速的Picard迭代法
void picard_aitken(double x0,double epsilon,double (*phi_x)(double x))
{
    double x,x1,x2;
    x = x0;
    for(int k=1;;k++)
    {
        x1 = x;
        x2 = ((x1*((*phi_x)(phi_x(x1))))-(((*phi_x)(x1))*((*phi_x)(x1))))/(((*phi_x)((*phi_x)(x1)))-(2*((*phi_x)(x1)))+x1);//埃德金加速的Picard迭代法迭代公式
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
    getchar();
}

//Picard迭代法（第二版）
void picard(double x0,double epsilon,double (*phi_x)(double x))
{
    double x,x1,x2;
    x = x0;
    for(int k=1;;k++)
    {
        x1 = x;
        x2 = (*phi_x)(x);//Picard迭代法迭代公式
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
    getchar();
}

//求导方法能否更加精简？
double derivative(double x,double epsilon,double (*f)(double x))
{
    double dx = 0.0001;
    double dy_0,dy_1;
    do
    {
        dy_0 = (((*f)(x+dx)) - ((*f)(x)))/dx;
        dx = 0.5*dx;
        dy_1 = (((*f)(x+dx))-((*f)(x)))/dx;
    }
    while(fabs(dy_0-dy_1)>epsilon);
    return dy_1;
}

//牛顿下山法
void newton_down(double x0,double epsilon,double (*f)(double x))
{
    double k,x,x1,x2,lamda;
    x = x0;
    for(int i=1;;i++)
    {
        x1 = x;
        k = derivative(x1,epsilon,f);
        lamda = 1;
        x2 = x1-(((*f)(x1))/k)*lamda;//牛顿下山法公式
        if(fabs(x2-x1)<=epsilon)
        {
            printf("%lf\n",x2);
            printf("k=%d\n",i);
            break;
        }
        else
        {
            if(((*f)(x2))*((*f)(x1))>0)
                lamda = 0.5*lamda;//0.5为缩放因子
                x2 = x1-(((*f)(x1))/k)*lamda;
            }
            x = x2;
            printf("%lf\n",x);
        }
    printf("calculation is done");
    getchar();
}

//牛顿迭代法（第二版）
void newton(double x0,double epsilon,double (*f)(double x))
{
    double k,x,x1,x2;
    x = x0;
    for(int i=1;;i++)
    {
        x1 = x;
        k = derivative(x1,epsilon,f);
        x2 = x1-(((*f)(x1))/k);//牛顿法公式
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
    getchar();
}

//二分法还有更加简洁的写法吗？
double dichotomy(double x1,double x2,double epsilon,double (*f)(double x))
{
    double x0;
    for(int k;;k++)
    {
        if(((*f)(x1))*((*f)(x2))<0)
        {
            x0 = (x1+x2)/2;
            if((x2-x1)<=epsilon)
            {
                break;
            }
            else
            {
                if(((*f)(x1)*((*f)(x0)))<0)
                {
                    x2 = x0;
                }
                else if(((*f)(x0))*((*f)(x2))<0)
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

//单点弦截法
void single_secant(double x0,double epsilon,double (*f)(double x))
{
    double k,x1,x2;
    k = derivative(x0,epsilon,f);
    x1 = x0-(((*f)(x0))/k);
    for(int i=1;;i++)
    {
        k = (((*f)(x1))-((*f)(x0)))/(x1-x0);
        x2 = x1-(((*f)(x1))/k);//单点弦截法公式
        if(fabs(x2-x1)<=epsilon)
        {
            printf("%lf\n",x2);
            printf("k=%d\n",i);
            break;
        }
        else
        {
            x1 = x2;
            printf("%lf\n",x1);
        }
    }
    printf("calculation is done");
    getchar();
}

//两点弦截法
void double_secant(double x0,double epsilon,double (*f)(double x))
{
    double k,x1,x2;
    k = derivative(x0,epsilon,f);
    x1 = x0-(((*f)(x0))/k);
    for(int i=1;;i++)
    {
        k = (((*f)(x1))-((*f)(x0)))/(x1-x0);
        x2 = x1-(((*f)(x1))/k);//两点弦截法迭代公式
        if(fabs(x2-x1)<=epsilon)
        {
            printf("%lf\n",x2);
            printf("k=%d\n",i);
            break;
        }
        else
        {
            x0 = x1;
            x1 = x2;
            printf("%lf\n",x2);
        }
    }
    printf("calculation is done");
    getchar();
}

//上三角矩阵的求解方法（线性方程组）
//矩阵U(n*n)为系数矩阵，B(1*n)为方程右侧常数矩阵
void up_triangle(double U[n][n],double B[n],int n)//n为矩阵阶数，自行输入
{
    double *x[n];
    double ux = 0;
    x[n] = B[n]/U[n][n];
    for(int i=n-1;i>=0;i--)
    {
        for(int j=i+1;j<=n;j++)
        {
            ux+ = U[i][j]*x[j];
        }
        x[i] = (B[i]-ux)/U[i][i];
    }
    for(int k=0;k<=n-1;k++)
    {
        printf("x[%d]=%lf\n",k+1,x[k]);
    }
    getchar();
}

//下三角矩阵的求解方法（线性方程组）
//矩阵U(n*n)为系数矩阵，B(1*n)为方程右侧常数矩阵
void down_triangle(double U[n][n],double B[n],int n)//n为矩阵阶数，自行输入
{
    double *x[n];
    double ux = 0;
    x[0] = B[0]/U[0][0];
    for(int i=2;i<=n;i++)
    {
        for(int j=1;j<=i-1;j++)
        {
            ux+ = U[i][j]*x[j];
        }
        x[i] = (B[i]-ux)/U[i][i];
    }
    for(int k=0;k<=n-1;k++)
    {
        printf("x[%d]=%lf\n",k+1,x[k]);
    }
    getchar();
}