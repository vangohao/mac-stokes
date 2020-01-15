#include "project.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
using namespace std;
int main(int argc, char* argv[])
{
    int n = argc > 1 ? atoi(argv[1]) : 64;
    int level = argc > 2 ? atoi(argv[2]): 6;
    dtype** u = new_2darray(n + 1, n);
    dtype** v = new_2darray(n, n + 1);
    dtype** p = new_2darray(n, n);
    dtype** f = new_2darray(2*n+1, 2*n+1);
    dtype** g = new_2darray(2*n+1, 2*n+1);
    dtype** u_exact = new_2darray(2*n+1, 2*n+1);
    dtype** v_exact = new_2darray(2*n+1, 2*n+1);
    dtype** rf = new_2darray(n+1, n);
    dtype** rg = new_2darray(n, n+1);
    dtype** rdiv = new_2darray(n, n);

    dtype * b = new_vector(2 * n + 1);
    dtype * t = new_vector(2 * n + 1);
    dtype * l = new_vector(2 * n + 1);
    dtype * r = new_vector(2 * n + 1);

    dtype** f_main = new_2darray(n + 1, n);
    dtype** g_main = new_2darray(n, n + 1);
    dtype * b_main = new_vector(n + 1);
    dtype * t_main = new_vector(n + 1);
    dtype * l_main = new_vector(n + 1);
    dtype * r_main = new_vector(n + 1);

    initproblem(n, u, v, p, f, g, b, t, l, r, u_exact, v_exact);
    dtype r0, r0_div;

    for(int i = 0; i<=n ; i++)
    {
        for(int j = 0; j< n; j++)
        {
            f_main[i][j] = f[2 * i][2 * j + 1];
            g_main[j][i] = g[2 * j + 1][2 * i];
        }
        b_main[i] = b[2 * i];
        t_main[i] = t[2 * i];
        l_main[i] = l[2 * i];
        r_main[i] = r[2 * i];
    }
    for(int i = 1; i< n; i++)
    {
        f_main[i][0] += b[2 * i] * n;
        f_main[i][n - 1] += t[2 * i] * n;
        g_main[0][i] += l[2 * i] * n;
        g_main[n - 1][i] += r[2 * i] * n;
    }

    residual(n, 0, u, v, p, f_main, g_main, b_main, t_main, l_main, r_main, rf, rg, rdiv, &r0, &r0_div);
    dtype e0;
    error(n, u, v, u_exact, v_exact, &e0);
    cout<<"initial residual: "<<r0<<" "<<r0_div<<endl;

    dtype res, res_div;
    res = r0;
    res_div = r0_div;
    int cnt = 0;
    
        print(f_main, n + 1, n, "main f");
        print(g_main, n, n + 1, "main g");
    while (res / r0 > 1e-8 && cnt < 10000)
    {
        vcycle(n, level, 1, &dgs_iteration, u, v, p, f, g, b, t, l, r);
        residual(n, 0, u, v, p, f_main, g_main, b_main, t_main, l_main, r_main, rf, rg, rdiv, &res, &res_div);
        cout<<"iteration "<<cnt<<", residual "<<res<<", "<<res_div<<endl;
        print(rf, n + 1, n,"residual f");
        print(rg, n, n + 1,"residual g");
        cnt ++;
    }
    dtype en;
    error(n, u, v, u_exact, v_exact, &en);

    cout<<"result error:" << en<<" "<<"initial error:"<<e0<<endl;
    cout<<"iteration cnt:" <<cnt<<endl;
    return 0;
}