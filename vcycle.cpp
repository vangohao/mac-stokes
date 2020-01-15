#include "project.h"
#include <iostream>
using namespace std;
int vcycle(int n, int level, int iter, smoother_type smoother, dtype ** u, dtype ** v, dtype ** p, dtype ** f, dtype **g, dtype * b, dtype * t, dtype * l, dtype * r)
{
    dtype** res = new_2darray(n, n);
    dtype *** u_cycle = (dtype ***) malloc(level * sizeof(dtype **));
    dtype *** v_cycle = (dtype ***) malloc(level * sizeof(dtype **));
    dtype *** p_cycle = (dtype ***) malloc(level * sizeof(dtype **));
    dtype *** u_cycle_pro = (dtype ***) malloc(level * sizeof(dtype **));
    dtype *** v_cycle_pro = (dtype ***) malloc(level * sizeof(dtype **));
    dtype *** p_cycle_pro = (dtype ***) malloc(level * sizeof(dtype **));
    dtype *** f_cycle = (dtype ***) malloc(level * sizeof(dtype **));
    dtype *** g_cycle = (dtype ***) malloc(level * sizeof(dtype **));
    dtype *** rf_cycle = (dtype ***) malloc(level * sizeof(dtype **));
    dtype *** rg_cycle = (dtype ***) malloc(level * sizeof(dtype **));
    dtype *** rdiv_cycle = (dtype ***) malloc(level * sizeof(dtype **));
    
    dtype r0, r0_div;
    u_cycle[0] = u;
    v_cycle[0] = v;
    p_cycle[0] = p;
    
    u_cycle_pro[0] = new_2darray(n + 1, n);
    v_cycle_pro[0] = new_2darray(n, n + 1);
    p_cycle_pro[0] = new_2darray(n, n);
    f_cycle[0] = new_2darray(n+1, n);
    g_cycle[0] = new_2darray(n, n+1);
    rf_cycle[0] = new_2darray(n+1, n);
    rg_cycle[0] = new_2darray(n, n+1);
    rdiv_cycle[0] = new_2darray(n, n);
    for(int lvl = 1; lvl< level; lvl++)
    {
        u_cycle[lvl] = new_2darray((n >> (lvl)) + 1, n >> (lvl));
        v_cycle[lvl] = new_2darray(n >> (lvl), (n >> (lvl)) + 1);
        p_cycle[lvl] = new_2darray(n >> (lvl), n >> (lvl));
        u_cycle_pro[lvl] = new_2darray((n >> (lvl)) + 1, n >> (lvl));
        v_cycle_pro[lvl] = new_2darray(n >> (lvl), (n >> (lvl)) + 1);
        p_cycle_pro[lvl] = new_2darray(n >> (lvl), n >> (lvl));
        f_cycle[lvl] = new_2darray((n >> (lvl)) + 1, n >> (lvl));
        g_cycle[lvl] = new_2darray(n >> (lvl), (n >> (lvl)) + 1);
        rf_cycle[lvl] = new_2darray((n >> (lvl)) + 1, n >> (lvl));
        rg_cycle[lvl] = new_2darray(n >> (lvl), (n >> (lvl)) + 1);
        rdiv_cycle[lvl] = new_2darray(n >> (lvl), n >> (lvl));
    }
    for(int i = 0; i<=n ; i++)
    {
        for(int j = 0; j< n; j++)
        {
            f_cycle[0][i][j] = f[2 * i][2 * j + 1];
            g_cycle[0][j][i] = g[2 * j + 1][2 * i];
        }
    }

    (*smoother)(n, 0, u, v, p, f_cycle[0], g_cycle[0], b, t, l, r);
    for(int lvl = 1; lvl < level; lvl ++)
    {
        print(u_cycle[lvl - 1], (n >> (lvl - 1)) + 1, n >> (lvl - 1));
        residual(n, lvl-1 , u_cycle[lvl - 1], v_cycle[lvl - 1], p_cycle[lvl - 1], 
        f_cycle[lvl - 1], g_cycle[lvl - 1], b, t, l, r, rf_cycle[lvl - 1], rg_cycle[lvl - 1], rdiv_cycle[lvl - 1], &r0, &r0_div);
        restriction(n >> (lvl), lvl, f_cycle[lvl], g_cycle[lvl], p_cycle[lvl], rf_cycle[lvl - 1], rg_cycle[lvl - 1], p_cycle[lvl - 1]);
        for(int i = 0; i<iter; i++)
        {
            (*smoother)(n >> (lvl), lvl, u_cycle[lvl], v_cycle[lvl], p_cycle[lvl], f_cycle[lvl], g_cycle[lvl], b, t, l, r);
        }
    }
    for(int lvl = level - 1; lvl > 0; lvl --) 
    {
        prolongation(n >> (lvl), lvl, u_cycle[lvl], v_cycle[lvl], p_cycle[lvl], u_cycle_pro[lvl - 1], v_cycle_pro[lvl - 1], p_cycle_pro[lvl - 1]);
        correction(n >> (lvl -1), lvl-1, u_cycle[lvl - 1], v_cycle[lvl - 1], p_cycle[lvl - 1], u_cycle_pro[lvl - 1], v_cycle_pro[lvl - 1], p_cycle_pro[lvl - 1]);
    }
}