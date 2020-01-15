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

    dtype * b = new_vector(n + 1);
    dtype * t = new_vector(n + 1);
    dtype * l = new_vector(n + 1);
    dtype * r = new_vector(n + 1);

    initproblem(n, u, v, p, f, g, b, t, l, r, u_exact, v_exact);
    dtype r0, r0_div;

    residual(n, 0, u, v, p, f, g, b, t, l, r, rf, rg, rdiv, &r0, &r0_div);
    dtype e0;
    error(n, u, v, u_exact, v_exact, &e0);
    cout<<"initial residual: "<<r0<<" "<<r0_div<<endl;

    dtype res, res_div;
    res = r0;
    res_div = r0_div;
    int cnt = 0;
    while (res / r0 > 1e-8 && cnt < 100)
    {
        vcycle(n, level, 1, &dgs_iteration, u, v, p, f, g, b, t, l, r);
        residual(n, 0, u, v, p, f, g, b, t, l, r, rf, rg, rdiv, &res, &res_div);
        cout<<"iteration "<<cnt<<", residual "<<res<<", "<<res_div<<endl;
        cnt ++;
    }
    dtype en;
    error(n, u, v, u_exact, v_exact, &en);

    cout<<"result error:" << en<<" "<<"initial error:"<<e0<<endl;
    return 0;
}