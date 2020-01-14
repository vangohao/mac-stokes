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
            u[2 * i][0] = u[i][0] / 2;
        }
        u[2 * i + 2][0] = u[i + 1][0] / 2;
        u[2 * i + 1][0] = (u[2 * i][0] + u[2 * i + 2][0]) / 2;
        if (i == n)
        {
            u[2 * i][n] = u[i][n] / 2;
        }
        u[2 * i + 2][n] = u[i + 1][n] / 2;
        u[2 * i + 1][n] = (u[2 * i][n] + u[2 * i + 2][n]) / 2;
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
            v[2 * i][0] = v[i][0] / 2;
        }
        v[2 * i + 2][0] = v[i + 1][0] / 2;
        v[2 * i + 1][0] = (v[2 * i][0] + v[2 * i + 2][0]) / 2;
        if (i == n)
        {
            v[2 * i][n] = v[i][n] / 2;
        }
        v[2 * i + 2][n] = v[i + 1][n] / 2;
        v[2 * i + 1][n] = (v[2 * i][n] + v[2 * i + 2][n]) / 2;
    }
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n - 1; j++)
        {
            if (i == 0)
            {
                v[2 * i][2 * j + 1] = vr[i][j] * 0.75 + vr[i][j + 1] * 0.25;
                v[2 * i][2 * j + 2] = vr[i][j] * 0.25 + vr[i][j + 1] * 0.75;
            }
            v[2 * i + 2][2 * j + 1] = vr[i + 1][j] * 0.75 + vr[i + 1][j + 1] * 0.25;
            v[2 * i + 2][2 * j + 2] = vr[i + 1][j] * 0.25 + vr[i + 1][j + 1] * 0.75;
            v[2 *i + 1][2 * j + 1] = (v[2 * i][2 *j + 1] + v[2 * i + 2][2 * j + 1]) / 2.;
            v[2 *i + 1][2 * j + 2] = (v[2 * i][2 *j + 2] + v[2 * i + 2][2 * j + 2]) / 2.;
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