#include <iostream>
#include "project.h"
using namespace std;
int correction(int n, int level, dtype ** ur, dtype ** vr, dtype ** pr, dtype ** u, dtype ** v, dtype ** p)
{
    for(int i = 0; i<=n ;i++)
    {
        for(int j = 0; j<n; j++)
        {
            ur[i][j] += u[i][j];
            vr[j][i] += v[j][i];
        }
    }
    for(int i = 0; i<n ;i++)
    {
        for(int j = 0; j<n; j++)
        {
            pr[i][j] += p[i][j];
        }
    }
}