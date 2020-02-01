#include "project.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

using namespace std;
int pcg_level;
int main(int argc, char* argv[])
{
    smoother_type* smt[4] = {&dgs_iteration, &uzawa_iteration, &inexact_uzawa_iteration, &pcg_uzawa_iteration};
    int smtnum = argc > 1? atoi(argv[1]) : 0;
    int n = argc > 2 ? atoi(argv[2]) : 64;
    int level = argc > 3 ? atoi(argv[3]): 1;
    if (smtnum == 3)
    {
        pcg_level = level;
        level = 1;
    }
    int v0 = argc > 4 ?atoi(argv[4]) : 3;
    int v1 = argc > 5 ?atoi(argv[5]) : 3;
    int mgiter = argc > 6 ?atoi(argv[6]) : 1; //only available for pcg_uzawa
    int maxcnt = argc > 7 ?atoi(argv[7]) : 10000;
    if (smtnum != 3 && mgiter != 1)
    {
        fprintf(stderr, "only method 3 can set mgiter!\n");
        mgiter = 1;
    }
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

    int cnt;

    dtype e0;
    error(n, u, v, u_exact, v_exact, &e0);

    clock_t start, end;
    time_t wstart,wend;
    double cpu_time_used, wall_time_used;

    start = clock();
    wstart = time(NULL);

    cnt = vcycle(n, level, mgiter, v0, v1, maxcnt, smt[smtnum], u, v, p, f, g, b, t, l, r);

    end = clock();
    wend = time(NULL);

    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    wall_time_used = difftime(wend, wstart);
    
    dtype en;
    error(n, u, v, u_exact, v_exact, &en);

    cout<<"result error:" << en<<" "<<"initial error:"<<e0<<endl;
    cout<<"iteration cnt:" <<cnt<<endl;
    cout<<"CPU time:"<<cpu_time_used<<endl;
    cout<<"Wall time:"<<wall_time_used<<endl;
    return 0;
}