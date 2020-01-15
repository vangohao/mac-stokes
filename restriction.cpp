#include <iostream>
#include "project.h"
using namespace std;
int restriction(int n, int level, dtype ** ur, dtype ** vr, dtype ** pr, dtype ** u, dtype ** v, dtype ** p)
{
    // for u
    for(int j = 0; j < n; j++)
    {
        // FIXME
        ur[0][j] = (2 * u[0][2 * j] + u[1][2 * j] 
                + 2 * u[0][2 * j + 1] + u[1][2 * j + 1]) / 8.;
        ur[n][j] = (2 * u[2 * n][2 * j] + u[2 * n - 1][2 * j] 
                + 2 * u[2 * n][2 * j + 1] + u[2 * n - 1][2 * j + 1]) / 8.;
    }
    for(int i = 1; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            ur[i][j] = (u[2 * i - 1][2 * j] + 2 * u[2 * i][2 * j] + u[2 * i + 1][2 * j] 
                + u[2 * i - 1][2 * j + 1] + 2 * u[2 * i][2 * j + 1] + u[2 * i + 1][2 * j + 1]) / 8.;
        }
    }

    //for v
    for(int j = 0; j < n; j++)
    {
        // FIXME
        vr[j][0] = (2 * v[2 * j][0] + v[2 * j][1] 
                + 2 * v[2 * j + 1][0] + v[2 * j + 1][1]) / 8.;
        vr[j][n] = (2 * v[2 * j][2 * n] + v[2 * j][2 * n - 1] 
                + 2 * v[2 * j + 1][2 * n] + v[2 * j + 1][2 * n - 1]) / 8.;
    }
    for(int i = 1; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            vr[j][i] = (v[2 * j][2 * i - 1] + 2 * v[2 * j][2 * i] + v[2 * j][2 * i + 1] 
                + v[2 * j + 1][2 * i - 1] + 2 * v[2 * j + 1][2 * i] + v[2 * j + 1][2 * i + 1]) / 8.;
        }
    }

    //for p
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            pr[i][j] = (p[2 * i][2 * j] + p[2 * i + 1][2 * j] + p[2 * i][2 * j + 1] + p[2 * i + 1][2 * j + 1]) / 4.;
        }
    }
    
    return 1;
}
