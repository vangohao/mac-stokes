#include "project.h"
#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
using namespace std;
#define sqr(x) ((x)*(x))
int residual(int n, int level, dtype ** u, dtype ** v, dtype ** p, dtype ** f, dtype **g, dtype ** d, dtype * b, dtype * t, dtype * l, dtype * r, dtype ** rf, dtype ** rg, dtype ** rdiv, dtype * r0, dtype * r1)
{
    dtype h = 1. / n;
    dtype tmp, delta;
    //update u
    // j = 0;
    for(int i = 1; i< n;i++)
    {
        rf[i][0] = f[i][0] + (- (p[i][0] - p[i - 1][0]) * h /* + b[i] * h */ + u[i][1] + u[i + 1][0] + u[i - 1][0] - u[i][0] * 3) / (h * h);
    }
    for(int j = 1; j < n - 1; j++)
    {
        for(int i = 1; i < n; i++)
        {
            rf[i][j] = f[i][j] + ( - (p[i][j] - p[i -1][j]) * h + u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1] - 4 * u[i][j]) / (h * h);
        }
    }
    // j = n - 1;
    for(int i = 1; i< n;i++)
    {
        rf[i][n-1] = f[i][n-1] + (- (p[i][n-1] - p[i - 1][n-1]) * h /* + t[i] * h */ + u[i][n-2] + u[i + 1][n-1] + u[i - 1][n-1] - 3 * u[i][n - 1]) / (h * h);
    }

    //update v (j,i)
    // j = 0;
    for(int i = 1; i< n;i++)
    {
        rg[0][i] = g[0][i] + (- (p[0][i] - p[0][i - 1]) * h /* + l[i] * h */ + v[1][i] + v[0][i + 1] + v[0][i - 1] - 3 * v[0][i]) / (h * h);
    }
    for(int j = 1; j < n - 1; j++)
    {
        for(int i = 1; i < n; i++)
        {
            rg[j][i] = g[j][i] + (- (p[j][i] - p[j][i -1]) * h + v[j][i+1] + v[j][i-1] + v[j+1][i] + v[j-1][i] - 4 * v[j][i]) / (h * h);
        }
    }
    // j = n - 1;
    for(int i = 1; i< n;i++)
    {
        rg[n-1][i] = g[n-1][i] + (- (p[n-1][i] - p[n-1][i - 1]) * h /* + r[i] * h */ + v[n-2][i] + v[n-1][i + 1] + v[n-1][i - 1] - 3 * v[n-1][i]) / (h * h);
    }

    //calculate divergence residual
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            rdiv[i][j] = - (u[i+1][j] - u[i][j]) / h - (v[i][j+1] - v[i][j]) / h + d[i][j];
        }
    }

    double r0_sq = 0;
    for(int i = 0; i <= n; i++)
        for(int j = 0; j < n; j++)
        {
            r0_sq += sqr(rf[i][j]) + sqr(rg[j][i]);
        }
    *r0 = sqrt(r0_sq);

    double r1_sq = 0;
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
        {
            r1_sq += sqr(rdiv[i][j]);
        }
    *r1 = sqrt(r1_sq);

    return 1;
}