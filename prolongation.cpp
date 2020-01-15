#include <iostream>
#include "project.h"
using namespace std;
int prolongation(int n, int level, dtype ** ur, dtype ** vr, dtype ** pr, dtype ** u, dtype ** v, dtype ** p)
{
    //for u
    for(int i = 0; i < n; i++)
    {
        if (i == 0)
        {
            u[2 * i][0] = ur[i][0] / 2;
        }
        u[2 * i + 2][0] = ur[i + 1][0] / 2;
        u[2 * i + 1][0] = (u[2 * i][0] + u[2 * i + 2][0]) / 2;
    }
    for(int i = 0; i < n; i++)
    {
        if (i == 0)
        {
            u[2 * i][2 * n - 1] = ur[i][n - 1] / 2;
        }
        u[2 * i + 2][2 * n - 1] = ur[i + 1][n - 1] / 2;
        u[2 * i + 1][2 * n - 1] = (u[2 * i][2 * n - 1] + u[2 * i + 2][2 * n - 1]) / 2;
    }
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n - 1; j++)
        {
            if (i == 0)
            {
                u[2 * i][2 * j + 1] = ur[i][j] * 0.75 + ur[i][j + 1] * 0.25;
                u[2 * i][2 * j + 2] = ur[i][j] * 0.25 + ur[i][j + 1] * 0.75;
            }
            u[2 * i + 2][2 * j + 1] = ur[i + 1][j] * 0.75 + ur[i + 1][j + 1] * 0.25;
            u[2 * i + 2][2 * j + 2] = ur[i + 1][j] * 0.25 + ur[i + 1][j + 1] * 0.75;
            u[2 *i + 1][2 * j + 1] = (u[2 * i][2 *j + 1] + u[2 * i + 2][2 * j + 1]) / 2.;
            u[2 *i + 1][2 * j + 2] = (u[2 * i][2 *j + 2] + u[2 * i + 2][2 * j + 2]) / 2.;
        }
    }

    //for v
    for(int i = 0; i < n; i++)
    {
        if (i == 0)
        {
            v[0][2 * i] = vr[0][i] / 2;
        }
        v[0][2 * i + 2] = vr[0][i + 1] / 2;
        v[0][2 * i + 1] = (v[0][2 * i] + v[0][2 * i + 2]) / 2;
    }
    for(int i = 0; i < n; i++)
    {
        if (i == 0)
        {
            v[2 * n - 1][2 * i] = vr[n - 1][i] / 2;
        }
        v[2 * n - 1][2 * i + 2] = vr[n - 1][i + 1] / 2;
        v[2 * n - 1][2 * i + 1] = (v[n - 1][2 * i] + v[n - 1][2 * i + 2]) / 2;
    }
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n - 1; j++)
        {
            if (i == 0)
            {
                v[2 * j + 1][2 * i] = vr[j][i] * 0.75 + vr[j + 1][i] * 0.25;
                v[2 * j + 2][2 * i] = vr[j][i] * 0.25 + vr[j + 1][i] * 0.75;
            }
            v[2 * j + 1][2 * i + 2] = vr[j][i + 1] * 0.75 + vr[j + 1][i + 1] * 0.25;
            v[2 * j + 2][2 * i + 2] = vr[j][i + 1] * 0.25 + vr[j + 1][i + 1] * 0.75;
            v[2 * j + 1][2 *i + 1] = (v[2 *j + 1][2 * i] + v[2 * j + 1][2 * i + 2]) / 2.;
            v[2 * j + 2][2 *i + 1] = (v[2 *j + 2][2 * i] + v[2 * j + 2][2 * i + 2]) / 2.;
        }
    }

    //for p
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
        {
            p[2 * i][2 * j] = pr[i][j];
            p[2 * i + 1][2 * j] = pr[i][j];
            p[2 * i][2 * j + 1] = pr[i][j];
            p[2 * i + 1][2 * j + 1] = pr[i][j];
        }
}