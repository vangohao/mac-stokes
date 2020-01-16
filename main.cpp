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

    dtype * b = new_vector(2 * n + 1);
    dtype * t = new_vector(2 * n + 1);
    dtype * l = new_vector(2 * n + 1);
    dtype * r = new_vector(2 * n + 1);
    
    dtype ** d = new_2darray(n, n);

    initproblem(n, u, v, p, f, g, b, t, l, r, u_exact, v_exact);
    dtype r0, r0_div;

    int iter = argc > 3 ?atoi(argv[3]) : 3;
    int maxcnt = argc > 4 ?atoi(argv[4]) : 10000;
    int cnt;

    dtype e0;
    error(n, u, v, u_exact, v_exact, &e0);

    cnt = vcycle(n, level, iter, maxcnt, &dgs_iteration, u, v, p, f, g, b, t, l, r);
    
    dtype en;
    error(n, u, v, u_exact, v_exact, &en);

    cout<<"result error:" << en<<" "<<"initial error:"<<e0<<endl;
    cout<<"iteration cnt:" <<cnt<<endl;
    return 0;
}