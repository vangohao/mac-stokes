#include <limits.h>
#include <stdio.h>
#include "project.h"
int pcg_uzawa_iteration(int n, int level, int mgiter, int mgv0, int mgv1, dtype ** u, dtype ** v, dtype ** p, dtype ** f, dtype **g, dtype ** d, dtype * b, dtype * t, dtype * l, dtype * r)
{
    // printf("begin UZAWA\n");
    dtype alpha = 1.;
    dtype h = 1. / n;
    dtype tmp, delta;
    dtype ** bf = new_2darray(n + 1, n);
    dtype ** bg = new_2darray(n, n + 1);
    dtype ** btu = new_2darray(n, n);
    //update u
    for(int i = 1; i< n;i++)
    {
        bf[i][0] = f[i][0] + (- (p[i][0] - p[i - 1][0]) * h) / (h * h);
    }
    for(int j = 1; j < n - 1; j++)
    {
        for(int i = 1; i < n; i++)
        {
            bf[i][j] = f[i][j] + ( - (p[i][j] - p[i -1][j]) * h) / (h * h);
        }
    }
    // j = n - 1;
    for(int i = 1; i< n;i++)
    {
        bf[i][n-1] = f[i][n-1] + (- (p[i][n-1] - p[i - 1][n-1]) * h) / (h * h);
    }

    //update v (j,i)
    for(int i = 1; i< n;i++)
    {
        bg[0][i] = g[0][i] + (- (p[0][i] - p[0][i - 1]) * h) / (h * h);
    }
    for(int j = 1; j < n - 1; j++)
    {
        for(int i = 1; i < n; i++)
        {
            bg[j][i] = g[j][i] + (- (p[j][i] - p[j][i -1]) * h) / (h * h);
        }
    }
    // j = n - 1;
    for(int i = 1; i< n;i++)
    {
        bg[n-1][i] = g[n-1][i] + (- (p[n-1][i] - p[n-1][i - 1]) * h) / (h * h);
    }
    // cg(n, level, 2 * n * (n - 1), 1e-9, u, v, bf, bg);
    int ll = 0, nn=n;
    while(nn % 2 == 0) {l ++; nn /= 2;}
    //uzawa
    printf("pcg level, %d %d %d %d\n", level, mgiter, mgv0, mgv1);
    pcg(n, level, mgiter, mgv0, mgv1, 1000, 1e-3, u, v, bf, bg);
    print(u, n + 1, n, "u");
    print(v, n, n + 1, "v");
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            btu[i][j] = (u[i+1][j] - u[i][j]) / h + (v[i][j+1] - v[i][j]) / h - d[i][j];
        }
    }
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            p[i][j] -= alpha * btu[i][j];
        }
    }
    print(p, n, n, "p");
    delete_2darray(bf, n + 1);
    delete_2darray(bg, n);
    delete_2darray(btu, n);
}