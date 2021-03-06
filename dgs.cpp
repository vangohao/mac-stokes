#include <iostream>
#include "project.h"
using namespace std;
int GS_uv(int n, int level, dtype ** u, dtype ** v, dtype ** f, dtype **g)
{
    dtype h = 1. / n;
    //update u
    // j = 0;
    for(int i = 1; i< n;i++)
    {
        u[i][0] = (f[i][0] * h * h + u[i][1] + u[i + 1][0] + u[i - 1][0]) / 3.;
    }
    for(int j = 1; j < n - 1; j++)
    {
        for(int i = 1; i < n; i++)
        {
            u[i][j] = (f[i][j] * h * h + u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1]) / 4.;
        }
    }
    // j = n - 1;
    for(int i = 1; i< n;i++)
    {
        u[i][n-1] = (f[i][n-1] * h * h + u[i][n-2] + u[i + 1][n-1] + u[i - 1][n-1]) / 3.;
    }

    //update v (j,i)
    // j = 0;
    for(int i = 1; i< n;i++)
    {
        v[0][i] = (g[0][i] * h * h + v[1][i] + v[0][i + 1] + v[0][i - 1]) / 3.;
    }
    for(int j = 1; j < n - 1; j++)
    {
        for(int i = 1; i < n; i++)
        {
            v[j][i] = (g[j][i] * h * h + v[j][i+1] + v[j][i-1] + v[j+1][i] + v[j-1][i]) / 4.;
        }
    }
    // j = n - 1;
    for(int i = 1; i< n;i++)
    {
        v[n-1][i] = (g[n-1][i] * h * h + v[n-2][i] + v[n-1][i + 1] + v[n-1][i - 1]) / 3.;
    }
}
int dgs_iteration(int n, int level, int mgiter, int mgv0, int mgv1, dtype ** u, dtype ** v, dtype ** p, dtype ** f, dtype **g, dtype ** d, dtype * b, dtype * t, dtype * l, dtype * r)
{
    dtype h = 1. / n;
    dtype tmp, delta;
    //update u
    // j = 0;
    for(int i = 1; i< n;i++)
    {
        u[i][0] = (f[i][0] * h * h - (p[i][0] - p[i - 1][0]) * h /* + b[i] * h */ + u[i][1] + u[i + 1][0] + u[i - 1][0]) / 3.;
    }
    for(int j = 1; j < n - 1; j++)
    {
        for(int i = 1; i < n; i++)
        {
            u[i][j] = (f[i][j] * h * h - (p[i][j] - p[i -1][j]) * h + u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1]) / 4.;
        }
    }
    // j = n - 1;
    for(int i = 1; i< n;i++)
    {
        u[i][n-1] = (f[i][n-1] * h * h - (p[i][n-1] - p[i - 1][n-1]) * h/*  + t[i] * h */ + u[i][n-2] + u[i + 1][n-1] + u[i - 1][n-1]) / 3.;
    }

    //update v (j,i)
    // j = 0;
    for(int i = 1; i< n;i++)
    {
        v[0][i] = (g[0][i] * h * h - (p[0][i] - p[0][i - 1]) * h /* + l[i] * h */ + v[1][i] + v[0][i + 1] + v[0][i - 1]) / 3.;
    }
    for(int j = 1; j < n - 1; j++)
    {
        for(int i = 1; i < n; i++)
        {
            v[j][i] = (g[j][i] * h * h - (p[j][i] - p[j][i -1]) * h + v[j][i+1] + v[j][i-1] + v[j+1][i] + v[j-1][i]) / 4.;
        }
    }
    // j = n - 1;
    for(int i = 1; i< n;i++)
    {
        v[n-1][i] = (g[n-1][i] * h * h - (p[n-1][i] - p[n-1][i - 1]) * h /* + r[i] * h  */+ v[n-2][i] + v[n-1][i + 1] + v[n-1][i - 1]) / 3.;
    }

    dtype res;
    //calculate internal residual
    for(int i = 1; i < n - 1; i++)
    {
        for(int j = 1; j < n - 1; j++)
        {
            res = - (u[i+1][j] - u[i][j]) / h - (v[i][j+1] - v[i][j]) / h + d[i][j];

            delta = res * h / 4;
            u[i][j] -= delta;
            u[i+1][j] += delta;
            v[i][j] -= delta;
            v[i][j+1] += delta;

            tmp = res / 4.;
            p[i][j] += res;
            p[i+1][j] -= tmp;
            p[i-1][j] -= tmp;
            p[i][j+1] -= tmp;
            p[i][j-1] -= tmp;
        }
    }

    //calculate edge residual (i, n - 1)
    for(int i = 1; i < n - 1; i++)
    {
        res = - (u[i+1][n-1] - u[i][n-1]) / h - (v[i][n] - v[i][n-1]) / h;
        delta = res * h / 3;
        u[i][n-1] -= delta;
        u[i+1][n-1] += delta;
        v[i][n-1] -= delta;
        tmp = res / 3.;
        p[i][n-1] += tmp * 4.;
        // p[i][n-1] += res;
        p[i+1][n-1] -= tmp;
        p[i-1][n-1] -= tmp;
        p[i][n-2] -= tmp;
    }

    //calculate edge residual (i, 0)
    for(int i = 1; i < n - 1; i++)
    {
        res = - (u[i+1][0] - u[i][0]) / h - (v[i][1] - v[i][0]) / h;
        delta = res * h / 3;
        u[i][0] -= delta;
        u[i+1][0] += delta;
        v[i][1] += delta;
        tmp = res / 3.;
        p[i][0] += tmp * 4.;
        // p[i][0] += res;
        p[i+1][0] -= tmp;
        p[i-1][0] -= tmp;
        p[i][1] -= tmp;
    }

    //calculate edge residual
    for(int j = 1; j < n - 1; j++)
    {
        res = - (u[n][j] - u[n-1][j]) / h - (v[n-1][j+1] - v[n-1][j]) / h;
        delta = res * h / 4;
        u[n-1][j] -= delta;
        v[n-1][j] -= delta;
        v[n-1][j+1] += delta;
        tmp = res / 3.;
        p[n-1][j] += tmp * 4.;
        // p[n-1][j] += res;
        p[n-2][j] -= tmp;
        p[n-1][j+1] -= tmp;
        p[n-1][j-1] -= tmp;
    }

    //calculate edge residual
    for(int j = 1; j < n - 1; j++)
    {
        res = - (u[1][j] - u[0][j]) / h - (v[0][j+1] - v[0][j]) / h;
        delta = res * h / 4;
        u[1][j] += delta;
        v[0][j] -= delta;
        v[0][j+1] += delta;
        tmp = res / 3.;
        p[0][j] += tmp * 4.;
        // p[0][j] += res;
        p[1][j] -= tmp;
        p[0][j+1] -= tmp;
        p[0][j-1] -= tmp;
    }

    // vertux 
    res = - (u[1][0] - u[0][0]) / h - (v[0][1] - v[0][0]) / h;
    delta = res * h / 2;
    u[1][0] += delta;
    v[0][1] += delta;
    p[0][0] += 2 * res;
    // p[0][0] += res;
    p[1][0] -= res / 2;
    p[0][1] -= res / 2;
    // vertux 
    res = - (u[n][n-1] - u[n-1][n-1]) / h - (v[n-1][n] - v[n-1][n-1]) / h;
    delta = res * h / 2;
    u[n-1][n-1] -= delta;
    v[n-1][n-1] -= delta;
    p[n-1][n-1] += 2 * res;
    // p[n-1][n-1] += res;
    p[n-2][n-1] -= res / 2;
    p[n-1][n-2] -= res / 2;
    // vertux 
    res = - (u[n][0] - u[n-1][0]) / h - (v[n-1][1] - v[n-1][0]) / h;
    delta = res * h / 2;
    u[n-1][0] -= delta;
    v[n-1][1] += delta;
    p[n-1][0] += 2 * res;
    // p[n-1][0] += res;
    p[n-2][0] -= res / 2;
    p[n-1][1] -= res / 2;
    // vertux 
    res = - (u[1][n-1] - u[0][n-1]) / h - (v[0][n] - v[0][n-1]) / h;
    delta = res * h / 2;
    u[1][n-1] += delta;
    v[0][n-1] -= delta;
    p[0][n-1] += 2 * res;
    // p[0][n-1] += res;
    p[1][n-1] -= res / 2;
    p[0][n-2] -= res / 2;

    return 1;
}