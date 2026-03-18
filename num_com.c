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
        x2 = x1 - (((*f)(x1)) / k) * lamda; // 牛顿下山法公式

        if (fabs(x2 - x1) <= epsilon)
        {
            printf("%lf\n", x2);
            printf("k=%d\n", i);
            break;
        }
        else
        {
            if (((*f)(x2)) * ((*f)(x1)) > 0) {
                lamda *= 0.5; // 0.5 为缩放因子
                x2 = x1 - (((*f)(x1)) / k) * lamda;
            }
            x = x2;
            printf("%lf\n", x);
        }
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
void double_secant(double x0,double x00,double epsilon,double (*f)(double x))
{
    double k,x1,x2,x3;
    x1 = x0;
    k = derivative(x1,epsilon,f);
    x2 = x1-(((*f)(x1))/k);
    printf("x1=%lf,x2=%lf\n",x1,x2);
    for(int i=1;;i++)
    {
        k = (((*f)(x2))-((*f)(x1)))/(x2-x1);
        x3 = x2-(((*f)(x2))/k);//两点弦截法迭代公式
        if(fabs(x3-x2)<=epsilon)
        {
            printf("%lf\n",x3);
            printf("k=%d\n",i);
            break;
        }
        else
        {
            x1 = x2;
            x2 = x3;
            printf("%lf\n",x2);
        }
    }
    printf("calculation is done");
    getchar();
}

//秦九韶算法
;

//上三角矩阵的求解方法（线性方程组）
//U[n][n] 为系数矩阵（按行存放），B[n] 为右侧常数向量
void up_triangle(double **U, double *B, int n)
{
    double *x = (double*)malloc(n * sizeof(double));
    if (!x) 
    {
        fprintf(stderr, "out of memory\n");
        return;
    }
    for (int i = n - 1; i >= 0; --i) 
    {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) {
            sum += U[i][j] * x[j];
        }
        x[i] = (B[i] - sum) / U[i][i];
    }
    for (int k = 0; k < n; ++k) 
    {
        printf("x[%d] = %lf\n", k + 1, x[k]);
    }
    free(x);
    getchar();
}

//下三角矩阵的求解方法（线性方程组）
//L[n][n] 为系数矩阵（按行存放），B[n] 为右侧常数向量
void down_triangle(double **L, double *B, int n)
{
    double *x = (double*)malloc(n * sizeof(double));
    if (!x) 
    {
        fprintf(stderr, "out of memory\n");
        return;
    }
    for (int i = 0; i < n; ++i) 
    {
        double sum = 0.0;
        for (int j = 0; j < i; ++j) {
            sum += L[i][j] * x[j];
        }
        x[i] = (B[i] - sum) / L[i][i];
    }
    for (int k = 0; k < n; ++k) {
        printf("x[%d] = %lf\n", k + 1, x[k]);
    }
    free(x);
    getchar();
}

//高斯消元法（未换主元版）
void gauss_elimination(double **A,int n)
{
    double *x = (double*)malloc(n * sizeof(double));
    if(!x)
    {
        fprintf(stderr,"out of memory\n");
    }
    for(int k=1;k<=n-1;k++)
    {
        if(A[k][k]==0)
        {
            printf("pivot element equals error");
            break;
        }
        else
        {
            for(int i=k+1;i<=n;i++)
            {
                A[i][k] = A[i][k]/A[k][k];
                for(int j=k+1;k<=n+1;j++)
                {
                    A[i][j] = A[i][j]-(A[i][k]*A[k][j]);
                }
            }
        }
    }
    x[n-1] = A[n-1][n]/A[n-1][n-1];
    for(int i=n-2;i>=0;i--)
    {
        double sum=0;
        for(int j=i+1;j<=n;j++)
        {
            sum += A[i][j]*x[j];
        }
        x[i] = (A[i][n]-sum)/A[i][i];
    }
    for(int k=0;k<=n-1;k++)
    {
        printf("%lf\n",x[k]);
    }
    free(x);
    getchar();
}

int findmax(double **A,int n,int k)//用于高斯列主元消元法，查找返回行序数
{
    int t = 0;
    double max = A[0][k];
    for(int i=0;i<=n-1;i++)
    {
        if(A[i][k]>=max)
        {
            max = A[i][k];
            t = i;
        }
    }
    return t;
}

//高斯消元法（列主元版）
void gauss_elimination_pivot(double **A, int n)
{
    double *x = (double*)malloc(n * sizeof(double));
    if(!x)
    {
        fprintf(stderr,"out of memory\n");
    }
    for(int k=0;k<=n-2;k++)
    {
        int t = findmax(*A,n,k);
        if(t != k)
        {
            double swap;
            for(int i=0;i<=n-1;i++)
            {
                swap = A[t][i];
                A[t][i] = A[k][i];
                A[k][i] = swap;
            }
        }
        for(int i=k+1;i<=n+1;i++)
        {
            A[i][k] = A[i][k]/A[k][k];
            for(int j=k+1;j<=n;j++)
            {
                A[i][j] = A[i][j]-(A[i][k]*A[k][j]);
            }
        }
    }
    x[n-1] = A[n-1][n]/A[n-1][n-1];
    for(int i=n-2;i>=0;i--)
    {
        double sum=0;
        for(int j=i+1;j<=n;j++)
        {
            sum += A[i][j]*x[j];
        }
        x[i] = (A[i][n]-sum)/A[i][i];
    }
    for(int k=0;k<=n-1;k++)
    {
        printf("%lf\n",x[k]);
    }
    free(x);
    getchar();
}

//Guass_Jordan
void gauss_jordan()
{
    ;
}

//Doolite分解
void doolite(double **A,int n,int m)//n为行数，m为列数
{
    if(A[0][0] == 0)
    {
        printf("pivot element equals error");
    }
    else
    {
        for(int i=1;i<=n-1;i++)
        {
            A[i][1] = A[i][1]/A[1][1];
        }
    }
    for(int k=1;k<=n-2;k++)
    {
        if(A[k][k] == 0)
        {
            printf("pivot element equals error");
            break;
        }
        else
        {
            for(int j=k;j<=m-1;j++)
            {
                static double sum=0;
                for(int s=0;s<=k-2;s++)
                {
                    sum += A[k][s]*A[s][j];
                }
                A[k][j] = A[k][j]-sum;
            }
            for(int i=k+1;i<=n-1;i++)
            {
                static double sum=0;
                for(int s=0;s<=k-2;s++)
                {
                    sum += A[i][s]/A[s][k];
                }
                A[i][k] = (A[i][k]-sum)/A[k][k];
            }
            if(A[0][0] == 0)
            {
                printf("pivot element equals error");
            }
            else
            {
                for(int i=1;i<=n-1;i++)
                {
                    A[i][1] = A[i][1]/A[1][1];
                }
            }
            for(int k=1;k<=n-2;k++)
            {
                if(A[k][k] == 0)
                {
                    printf("pivot element equals error");
                    break;
                }
                else
                {
                    for(int j=k;j<=m-1;j++)
                    {
                        static double sum=0;
                        for(int s=0;s<=k-2;s++)
                        {
                            sum += A[k][s]*A[s][j];
                        }
                        A[k][j] = A[k][j]-sum;
                    }
                    for(int i=k+1;i<=n-1;i++)
                    {
                        static double sum=0;
                        for(int s=0;s<=k-2;s++)
                        {
                            sum += A[i][s]/A[s][k];
                        }
                        A[i][k] = (A[i][k]-sum)/A[k][k];
                    }
                }
            }
        }
    }
    return **A;
}

//LL分解
void ll(double **A,int n)
{
    double **l;
    l = (double**)malloc(n * sizeof(double*));
    for(int i = 0; i < n; i++)
    {
        l[i] = (double*)malloc(n * sizeof(double));
    }
    for(int i = 0; i < n; i++)
    {
        A[i][0] = A[i][0]/l[0][0];
    }
    for(int j=1;j<=n-2;j++)
    {
        double sum = 0;
        for(int k=0;k<=j-1;k++)
        {
            sum += (l[j][k]*l[j][k]);
        }
        l[j][j] = sqrt(A[j][j]-sum);
        A[j][j] = l[j][j];
        int i = j+1;
        double sum_0 = 0;
        for(int k=0;k<=j-1;k++)
        {
            sum_0 += (l[i][k]*l[j][k]);
        }
        l[i][j] = (A[i][j]-sum_0)/A[j][j];
        A[i][j] = l[i][j];
    }
    double sum_1 = 0;
    for(int k=0;k<=n-2;k++)
    {
        sum_1 += (l[n-1][k]*l[n-1][k]);
    }
    l[n-1][n-1] = sqrt(A[n][n]-sum_1);
    A[n-1][n-1] = l[n-1][n-1];
    free(l);
    return **A;
}

//LDL分解
void ldl(double **A, int n)
{
    for(int i=1;i<=n-1;i++)
    {
        for(int j=0;j<=i-1;j++)
        {
            double sum_0 =0;
            for(int k=0;k<=j-1;k++)
            {
                sum_0 += A[k][i]*A[j][k];
            }
            A[j][i] = A[i][j]-sum_0;
        }
        for(int j=0;j<=i-1;j++)
        {
            A[i][j] = A[j][i]/A[j][j];
        }
        double sum_1 = 0;
        for(int k=0;k<=i-2;k++)
        {
            sum_1 += A[k][i]*A[i][k];
        }
        A[i][i] = A[i][i]-sum_1;
    }
    return **A;
}

//追赶法（求解对角阵）
void chasing_method(double *a, double *b, double *c, double *d, int n)
{
    double *x = (double*)malloc(n * sizeof(double));
    if (!x) 
    {
        fprintf(stderr, "out of memory\n");
        return;
    }
    double *l = (double*)malloc((n-1) * sizeof(double));
    double *y = (double*)malloc((n-1) * sizeof(double));
    if (!l || !y) 
    {
        fprintf(stderr, "out of memory\n");
        free(l);
        free(y);
        return;
    }
    double *beta = (double*)malloc(n * sizeof(double));
    if (!beta) 
    {
        fprintf(stderr, "out of memory\n");
        free(beta);
        return;
    }
    beta[0] = b[0];
    y[0] = d[0];
    for(int i=1;i<n;i++)
    {
        l[i-1] = a[i]/beta[i-1];
        beta[i] = b[i]-(l[i-1]*c[i-1]);
        y[i] = d[i]-(l[i-1]*y[i-1]);
    }
    x[n-1] = y[n-1]/beta[n-1];
    for(int j=n-2;j>=0;j--)
    {
        x[j] = (y[j]-(c[j]*x[j+1]))/beta[j];
    }
    for (int k = 0; k < n; ++k)
    {
        printf("x[%d] = %lf\n",k+1,x[k]);
    }
    free(x);
    free(l);
    free(y);
    free(beta);
    getchar();
}