#include <omp.h>
#include <math.h>
#include <stdio.h>
#include "project.h"
int cg(int n, int level, int kmax, dtype eps, dtype **u, dtype **v, dtype **bf, dtype **bg)
{
    int k = 0;
    dtype rho = 0, bn = 0, rho1, beta, alpha, ptw;
    dtype h = 1. / n;
    dtype **pu, **pv, **wu, **wv;
    bn=0;
    rho=0;
    #pragma omp parallel for reduction(+:bn)
    for(int i = 0; i < n + 1; ++i)
        for(int j = 0; j < n; ++j)
        {
            bn += bf[i][j] * bf[i][j];
            // bn+=1;
        }
    
    #pragma omp parallel for reduction(+:bn)
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n + 1; ++j)
        {
            bn += bg[i][j] * bg[i][j];
            // bn+=1;
        }
    // printf("bn: %lf\n",bn);
    //update u
    // j = 0;
    #pragma omp parallel for
    for(int i = 1; i< n;i++)
    {
        bf[i][0] += (u[i][1] + u[i + 1][0] + u[i - 1][0] - u[i][0] * 3) / (h * h);
    }
    #pragma omp parallel for
    for(int j = 1; j < n - 1; j++)
    {
        for(int i = 1; i < n; i++)
        {
            bf[i][j] += (u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1] - 4 * u[i][j]) / (h * h);
        }
    }
    // j = n - 1;
    #pragma omp parallel for
    for(int i = 1; i< n;i++)
    {
        bf[i][n-1] += (u[i][n-2] + u[i + 1][n-1] + u[i - 1][n-1] - 3 * u[i][n - 1]) / (h * h);
    }

    //update v (j,i)
    // j = 0;
    #pragma omp parallel for
    for(int i = 1; i< n;i++)
    {
        bg[0][i] += (v[1][i] + v[0][i + 1] + v[0][i - 1] - 3 * v[0][i]) / (h * h);
    }
    #pragma omp parallel for
    for(int j = 1; j < n - 1; j++)
    {
        for(int i = 1; i < n; i++)
        {
            bg[j][i] += (v[j][i+1] + v[j][i-1] + v[j+1][i] + v[j-1][i] - 4 * v[j][i]) / (h * h);
        }
    }
    // j = n - 1;
    #pragma omp parallel for
    for(int i = 1; i< n;i++)
    {
        bg[n-1][i] += (v[n-2][i] + v[n-1][i + 1] + v[n-1][i - 1] - 3 * v[n-1][i]) / (h * h);
    }

    #pragma omp parallel for reduction(+:rho)
    for(int i = 0; i < n + 1; ++i)
        for(int j = 0; j < n; ++j)
        {
            rho += bf[i][j] * bf[i][j];
        }
    #pragma omp parallel for reduction(+:rho)
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n + 1; ++j)
        {
            rho += bg[i][j] * bg[i][j];
        }
    
    pu = new_2darray(n + 1, n);pv = new_2darray(n, n + 1);
    wu = new_2darray(n + 1, n);wv = new_2darray(n, n + 1);
    #pragma omp parallel for 
    for(int i = 0; i < n + 1; ++i)
        for(int j = 0; j < n; ++j)
        {
            pu[i][j] = bf[i][j];
        }
    #pragma omp parallel for 
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n + 1; ++j)
        {
            pv[i][j] = bg[i][j];
        }
    while ( rho > (eps * eps) * bn && k < kmax)
    {
        // printf("CG iterateion %d, rho/bn=%lf, rho=%lf, bn=%lf\n", k, rho/bn, rho, bn);
        print(u, n + 1, n, "u");
        print(v, n, n + 1, "v");
        if (k)
        {
            beta = rho / rho1;
    #pragma omp parallel for 
            for(int i = 0; i < n + 1; ++i)
                for(int j = 0; j < n; ++j)
                {
                    pu[i][j] = bf[i][j] + beta * pu[i][j];
                }
    #pragma omp parallel for 
            for(int i = 0; i < n; ++i)
                for(int j = 0; j < n + 1; ++j)
                {
                    pv[i][j] = bg[i][j] + beta * pv[i][j];
                }
        }
        //wu
        // j = 0;
    #pragma omp parallel for 
        for(int i = 1; i< n;i++)
        {
            wu[i][0] = -(pu[i][1] + pu[i + 1][0] + pu[i - 1][0] - pu[i][0] * 3) / (h * h);
        }
    #pragma omp parallel for 
        for(int j = 1; j < n - 1; j++)
        {
            for(int i = 1; i < n; i++)
            {
                wu[i][j] = -(pu[i+1][j] + pu[i-1][j] + pu[i][j+1] + pu[i][j-1] - 4 * pu[i][j]) / (h * h);
            }
        }
        // j = n - 1;
    #pragma omp parallel for 
        for(int i = 1; i< n;i++)
        {
            wu[i][n-1] = -(pu[i][n-2] + pu[i + 1][n-1] + pu[i - 1][n-1] - 3 * pu[i][n - 1]) / (h * h);
        }

        //wv
        // j = 0;
    #pragma omp parallel for 
        for(int i = 1; i< n;i++)
        {
            wv[0][i] = -(pv[1][i] + pv[0][i + 1] + pv[0][i - 1] - 3 * pv[0][i]) / (h * h);
        }
    #pragma omp parallel for 
        for(int j = 1; j < n - 1; j++)
        {
            for(int i = 1; i < n; i++)
            {
                wv[j][i] = -(pv[j][i+1] + pv[j][i-1] + pv[j+1][i] + pv[j-1][i] - 4 * pv[j][i]) / (h * h);
            }
        }
        // j = n - 1;
    #pragma omp parallel for 
        for(int i = 1; i< n;i++)
        {
            wv[n-1][i] = -(pv[n-2][i] + pv[n-1][i + 1] + pv[n-1][i - 1] - 3 * pv[n-1][i]) / (h * h);
        }

        ptw = 0;
    #pragma omp parallel for reduction(+:ptw)
        for(int i = 0; i < n + 1; ++i)
            for(int j = 0; j < n; ++j)
            {
                ptw += pu[i][j] * wu[i][j];
            }
    #pragma omp parallel for reduction(+:ptw)
        for(int i = 0; i < n; ++i)
            for(int j = 0; j < n + 1; ++j)
            {
                ptw += pv[i][j] * wv[i][j];
            }
        alpha = rho / ptw;
    #pragma omp parallel for 
        for(int i = 0; i < n + 1; ++i)
            for(int j = 0; j < n; ++j)
            {
                u[i][j] += alpha * pu[i][j];
            }
    #pragma omp parallel for 
        for(int i = 0; i < n; ++i)
            for(int j = 0; j < n + 1; ++j)
            {
                v[i][j] += alpha * pv[i][j];
            }

    #pragma omp parallel for 
        for(int i = 0; i < n + 1; ++i)
            for(int j = 0; j < n; ++j)
            {
                bf[i][j] -= alpha * wu[i][j];
            }
    #pragma omp parallel for 
        for(int i = 0; i < n; ++i)
            for(int j = 0; j < n + 1; ++j)
            {
                bg[i][j] -= alpha * wv[i][j];
            }
        
        rho1 = rho;
        rho = 0;
    #pragma omp parallel for reduction(+:rho)
        for(int i = 0; i < n + 1; ++i)
            for(int j = 0; j < n; ++j)
            {
                rho += bf[i][j] * bf[i][j];
            }
    #pragma omp parallel for reduction(+:rho)
        for(int i = 0; i < n; ++i)
            for(int j = 0; j < n + 1; ++j)
            {
                rho += bg[i][j] * bg[i][j];
            }
        k++;
    }
    delete_2darray(pu, n + 1);
    delete_2darray(pv, n);
    delete_2darray(wu, n + 1);
    delete_2darray(wv, n);
}