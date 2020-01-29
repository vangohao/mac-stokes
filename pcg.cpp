#include <omp.h>
#include <math.h>
#include <stdio.h>
#include "project.h"
int pcg(int n, int level, int mgiter, int mgv0, int mgv1, int kmax, dtype eps, dtype **u, dtype **v, dtype **bf, dtype **bg)
{
    int k = 0;
    dtype rho = 0, bn = 0, mu, mu1, mu_tmp, beta, alpha, ptw;
    dtype h = 1. / n;
    dtype **pu, **pv, **wu, **wv, **zu, ** zv;
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

    // #pragma omp parallel for reduction(+:rho)
    for(int i = 0; i < n + 1; ++i)
        for(int j = 0; j < n; ++j)
        {
            rho += bf[i][j] * bf[i][j];
        }
    // #pragma omp parallel for reduction(+:rho)
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n + 1; ++j)
        {
            rho += bg[i][j] * bg[i][j];
        }
    bn = rho; //FIXME
    
    pu = new_2darray(n + 1, n);pv = new_2darray(n, n + 1);
    wu = new_2darray(n + 1, n);wv = new_2darray(n, n + 1);
    zu = new_2darray(n + 1, n);zv = new_2darray(n, n + 1);

    while ( rho > (eps * eps) * bn && k < kmax)
    {
        printf("CG iterateion %d, rho/bn=%lf, rho=%lf, bn=%lf\n", k, rho/bn, rho, bn);
        print(u, n + 1, n, "u");
        print(v, n, n + 1, "v");

        // #pragma omp parallel for reduction(+:rho)
        // for(int i = 0; i < n + 1; ++i)
        //     for(int j = 0; j < n; ++j)
        //     {
        //         zu[i][j] = bf[i][j];
        //     }
        // #pragma omp parallel for reduction(+:rho)
        // for(int i = 0; i < n; ++i)
        //     for(int j = 0; j < n + 1; ++j)
        //     {
        //         zv[i][j] = bg[i][j];
        //     }
        vcycle_precondition(n, level, mgiter, mgv0, mgv1, zu, zv, bf, bg);
        

        mu_tmp = 0.;
        // #pragma omp parallel for reduction(+:rho)
        for(int i = 0; i < n + 1; ++i)
            for(int j = 0; j < n; ++j)
            {
                mu_tmp += zu[i][j] * bf[i][j];
            }
        // #pragma omp parallel for reduction(+:rho)
        for(int i = 0; i < n; ++i)
            for(int j = 0; j < n + 1; ++j)
            {
                mu_tmp += zv[i][j] * bg[i][j];
            }

        if (k == 0)
        {
            #pragma omp parallel for 
            for(int i = 0; i < n + 1; ++i)
                for(int j = 0; j < n; ++j)
                {
                    pu[i][j] = zu[i][j];
                }
            #pragma omp parallel for 
            for(int i = 0; i < n; ++i)
                for(int j = 0; j < n + 1; ++j)
                {
                    pv[i][j] = zv[i][j];
                }
            mu = mu_tmp;
        }
        else
        {
            mu1 = mu;
            mu = mu_tmp;
            beta = mu / mu1;
    #pragma omp parallel for 
            for(int i = 0; i < n + 1; ++i)
                for(int j = 0; j < n; ++j)
                {
                    pu[i][j] = zu[i][j] + beta * pu[i][j];
                }
    #pragma omp parallel for 
            for(int i = 0; i < n; ++i)
                for(int j = 0; j < n + 1; ++j)
                {
                    pv[i][j] = zv[i][j] + beta * pv[i][j];
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
        alpha = mu / ptw;
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
    delete_2darray(zu, n + 1);
    delete_2darray(zv, n);
}