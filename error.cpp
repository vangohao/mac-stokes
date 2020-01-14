#include "project.h"
#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
using namespace std;
#define sqr(x) ((x)*(x))
int error(int n, dtype ** u, dtype ** v, dtype ** u_exact, dtype ** v_exact, dtype * en)
{
    dtype h = 1. / n;
    dtype en_sq = 0;
    for(int i = 1; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            en_sq += sqr(u[i][j] - u_exact[2 * i][2 * j + 1]);
            en_sq += sqr(v[j][i] - v_exact[2 * j + 1][2 * i]);
        }
    }
    *en = sqr(en_sq * h);
    return 0;
}