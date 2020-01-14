#include "project.h"
#include <iostream>
#include <stdlib.h>
#include <math.h>
using namespace std;
int initproblem(int n, dtype ** u, dtype ** v, dtype ** p, dtype ** f, dtype **g, dtype * b, dtype * t, dtype * l, dtype * r, dtype ** u_exact, dtype ** v_exact)
{
    dtype h = 1. / n;
    //f, g 总共要开辟 2n+1 * 2n+1 的空间
    for(int i = 0; i < 2 * n + 1; i ++)
    {
        for(int j = 0; j < 2 * n + 1; j++)
        {
            dtype x = i * h / 2;
            dtype y = j * h / 2;
            f[i][j] = -4 * M_PIl * M_PIl * (2 * cos(2 * M_PIl * x) - 1) * sin(2 * M_PIl * y) + x * x;
            g[i][j] = 4 * M_PIl * M_PIl * (2 * cos(2 * M_PIl * y) - 1) * sin(2 * M_PIl * x);
            u_exact[i][j] = (1 - cos(2 * M_PIl * x)) * sin(2 * M_PIl * y);
            v_exact[i][j] = -(1 - cos(2 * M_PIl * y)) * sin(2 * M_PIl * x);
        }
    }
    for(int i = 0; i < n + 1; i ++)
    {
        b[i] = 0;
        t[i] = 0;
        l[i] = 0;
        r[i] = 0;
    }
}