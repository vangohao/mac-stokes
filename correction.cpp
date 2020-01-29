#include <iostream>
#include "project.h"
using namespace std;
int correction_uv(int n, int level, dtype ** ur, dtype ** vr, dtype ** u, dtype ** v)
{
    for(int i = 0; i<=n ;i++)
    {
        for(int j = 0; j<n; j++)
        {
            ur[i][j] += u[i][j];
            vr[j][i] += v[j][i];
        }
    }
}

int correction(int n, int level, dtype ** ur, dtype ** vr, dtype ** pr, dtype ** u, dtype ** v, dtype ** p)
{
    correction_uv(n, level, ur, vr, u, v);
    for(int i = 0; i<n ;i++)
    {
        for(int j = 0; j<n; j++)
        {
            pr[i][j] += p[i][j];
        }
    }
}