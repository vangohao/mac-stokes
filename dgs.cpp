#include <iostream>
#include "project.h"
using namespace std;
int dgs_iteration(int n, dtype ** u, dtype ** v, dtype ** p, dtype ** f, dtype **g, dtype * b, dtype * t, dtype * l, dtype * r, dtype ** res)
{
    dtype h = 1. / n;
    dtype tmp, delta;
    //update u
    // j = 0;
    for(int i = 1; i< n;i++)
    {
        u[i][0] = (checkf(level, i,0) * h * h - (p[i][0] - p[i - 1][0]) * h + b[i] * h + u[i][1] + u[i + 1][0] + u[i - 1][0]) / 3.;
    }
    for(int j = 1; j < n - 1; j++)
    {
        for(int i = 1; i < n; i++)
        {
            u[i][j] = (checkf(level, i,j) * h * h - (p[i][j] - p[i -1][j]) * h + u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1]) / 4.;
        }
    }
    // j = n - 1;
    for(int i = 1; i< n;i++)
    {
        u[i][n-1] = (checkf(level, i,n-1) * h * h - (p[i][n-1] - p[i - 1][n-1]) * h + t[i] * h + u[i][n-2] + u[i + 1][n-1] + u[i - 1][n-1]) / 3.;
    }

    //update v (j,i)
    // j = 0;
    for(int i = 1; i< n;i++)
    {
        v[0][i] = (checkg(level, 0, i) * h * h - (p[0][i] - p[0][i - 1]) * h + l[i] * h + v[1][i] + v[0][i + 1] + v[0][i - 1]) / 3.;
    }
    for(int j = 1; j < n - 1; j++)
    {
        for(int i = 1; i < n; i++)
        {
            v[j][i] = (checkg(level, j, i) * h * h - (p[j][i] - p[j][i -1]) * h + v[j][i+1] + v[j][i-1] + v[j+1][i] + v[j-1][i]) / 4.;
        }
    }
    // j = n - 1;
    for(int i = 1; i< n;i++)
    {
        v[n-1][i] = (checkg(level, n-1, i) * h * h - (p[n-1][i] - p[n-1][i - 1]) * h + r[i] * h + v[n-2][i] + v[n-1][i + 1] + v[n-1][i - 1]) / 3.;
    }

    //calculate internal residual
    for(int i = 1; i < n - 1; i++)
    {
        for(int j = 1; j < n - 1; j++)
        {
            res[i][j] = - (u[i+1][j] - u[i][j]) / h - (v[i][j+1] - v[i][j]) / h;
        }
    }

    //update internal velocity
    for(int i = 1; i < n - 1;i++)
    {
        for(int j = 1; j < n - 1; j++)
        {
            delta = res[i][j] * h / 4;
            u[i][j] -= delta;
            u[i+1][j] += delta;
            v[i][j] -= delta;
            v[i][j+1] += delta;
        }
    }

    //update internal pressure
    for(int i = 1; i < n - 1;i++)
    {
        for(int j = 1; j < n - 1; j++)
        {
            tmp = res[i][j] / 4.;
            p[i][j] += res[i][j];
            p[i+1][j] -= tmp;
            p[i-1][j] -= tmp;
            p[i][j+1] -= tmp;
            p[i][j-1] -= tmp;
        }
    }

    //calculate edge residual (i, n - 1)
    for(int i = 1; i < n - 1; i++)
    {
        res[i][n-1] = - (u[i+1][n-1] - u[i][n-1]) / h - (v[i][n] - v[i][n-1]) / h;
    }

    //update edge velocity (i, n - 1)
    for(int i = 1; i < n - 1;i++)
    {
        delta = res[i][n-1] * h / 3;
        u[i][n-1] -= delta;
        u[i+1][n-1] += delta;
        v[i][n-1] -= delta;
    }

    //update internal pressure
    for(int i = 1; i < n - 1;i++)
    {
        tmp = res[i][n - 1] / 3.;
        p[i][n-1] += tmp * 4.;
        p[i+1][n-1] -= tmp;
        p[i-1][n-1] -= tmp;
        p[i][n-2] -= tmp;
    }

    //calculate edge residual (i, 0)
    for(int i = 1; i < n - 1; i++)
    {
        res[i][0] = - (u[i+1][0] - u[i][0]) / h - (v[i][1] - v[i][0]) / h;
    }

    //update edge velocity (i, 0)
    for(int i = 1; i < n - 1;i++)
    {
        delta = res[i][0] * h / 3;
        u[i][0] -= delta;
        u[i+1][0] += delta;
        v[i][1] += delta;
        // v[i][0] -= delta;
    }

    //update edge pressure
    for(int i = 1; i < n - 1;i++)
    {
        tmp = res[i][0] / 3.;
        p[i][0] += tmp * 4.;
        p[i+1][0] -= tmp;
        p[i-1][0] -= tmp;
        p[i][1] -= tmp;
    }

    //calculate edge residual
    for(int j = 1; j < n - 1; j++)
    {
        res[n-1][j] = - (u[n][j] - u[n-1][j]) / h - (v[n-1][j+1] - v[n-1][j]) / h;
    }

    //update edge velocity
    for(int j = 1; j < n - 1; j++)
    {
        delta = res[n-1][j] * h / 4;
        u[n-1][j] -= delta;
        v[n-1][j] -= delta;
        v[n-1][j+1] += delta;
    }

    //update edge pressure
    for(int j = 1; j < n - 1; j++)
    {
        tmp = res[n-1][j] / 3.;
        p[n-1][j] += tmp * 4.;
        p[n-2][j] -= tmp;
        p[n-1][j+1] -= tmp;
        p[n-1][j-1] -= tmp;
    }

    //calculate edge residual
    for(int j = 1; j < n - 1; j++)
    {
        res[0][j] = - (u[1][j] - u[0][j]) / h - (v[0][j+1] - v[0][j]) / h;
    }

    //update edge velocity
    for(int j = 1; j < n - 1; j++)
    {
        delta = res[0][j] * h / 4;
        u[1][j] += delta;
        v[0][j] -= delta;
        v[0][j+1] += delta;
    }

    //update edge pressure
    for(int j = 1; j < n - 1; j++)
    {
        tmp = res[0][j] / 3.;
        p[0][j] += tmp * 4.;
        p[1][j] -= tmp;
        p[0][j+1] -= tmp;
        p[0][j-1] -= tmp;
    }

    // vertux 
    res[0][0] = - (u[1][0] - u[0][0]) / h - (v[0][1] - v[0][0]) / h;
    delta = res[0][0] * h / 2;
    u[1][0] += delta;
    v[0][1] += delta;
    p[0][0] += 2 * res[0][0];
    p[1][0] -= res[0][0] / 2;
    p[0][1] -= res[0][0] / 2;
    // vertux 
    res[n-1][n-1] = - (u[n][n-1] - u[n-1][n-1]) / h - (v[n-1][n] - v[n-1][n-1]) / h;
    delta = res[n-1][n-1] * h / 2;
    u[n-1][n-1] -= delta;
    v[n-1][n-1] -= delta;
    p[n-1][n-1] += 2 * res[n-1][n-1];
    p[n-2][n-1] -= res[n-1][n-1] / 2;
    p[n-1][n-2] -= res[n-1][n-1] / 2;
    // vertux 
    res[n-1][0] = - (u[n][0] - u[n-1][0]) / h - (v[n-1][1] - v[n-1][0]) / h;
    delta = res[n-1][0] * h / 2;
    u[n-1][0] -= delta;
    v[n-1][1] += delta;
    p[n-1][0] += 2 * res[n-1][0];
    p[n-2][0] -= res[n-1][0] / 2;
    p[n-1][1] -= res[n-1][0] / 2;
    // vertux 
    res[0][n-1] = - (u[1][n-1] - u[0][n-1]) / h - (v[0][n] - v[0][n-1]) / h;
    delta = res[0][n-1] * h / 2;
    u[1][n-1] += delta;
    v[0][n-1] -= delta;
    p[0][n-1] += 2 * res[0][n-1];
    p[1][n-1] -= res[0][n-1] / 2;
    p[0][n-2] -= res[0][n-1] / 2;

    return 1;
}