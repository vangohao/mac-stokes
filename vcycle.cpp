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
    dtype ** b_cycle = (dtype **) malloc(level * sizeof(dtype *));
    dtype ** t_cycle = (dtype **) malloc(level * sizeof(dtype *));
    dtype ** l_cycle = (dtype **) malloc(level * sizeof(dtype *));
    dtype ** r_cycle = (dtype **) malloc(level * sizeof(dtype *));

    
    dtype r0, r0_div;
    u_cycle[0] = u;
    v_cycle[0] = v;
    p_cycle[0] = p;
    
    // u_cycle_pro[0] = new_2darray(n + 1, n);
    // v_cycle_pro[0] = new_2darray(n, n + 1);
    // p_cycle_pro[0] = new_2darray(n, n);
    // f_cycle[0] = new_2darray(n+1, n);
    // g_cycle[0] = new_2darray(n, n+1);
    // rf_cycle[0] = new_2darray(n+1, n);
    // rg_cycle[0] = new_2darray(n, n+1);
    // rdiv_cycle[0] = new_2darray(n, n);
    for(int lvl = 0; lvl< level; lvl++)
    {
        if (lvl>0)
        {
            u_cycle[lvl] = new_2darray((n >> (lvl)) + 1, n >> (lvl));
            v_cycle[lvl] = new_2darray(n >> (lvl), (n >> (lvl)) + 1);
            p_cycle[lvl] = new_2darray(n >> (lvl), n >> (lvl));
        }
        u_cycle_pro[lvl] = new_2darray((n >> (lvl)) + 1, n >> (lvl));
        v_cycle_pro[lvl] = new_2darray(n >> (lvl), (n >> (lvl)) + 1);
        p_cycle_pro[lvl] = new_2darray(n >> (lvl), n >> (lvl));
        f_cycle[lvl] = new_2darray((n >> (lvl)) + 1, n >> (lvl));
        g_cycle[lvl] = new_2darray(n >> (lvl), (n >> (lvl)) + 1);
        rf_cycle[lvl] = new_2darray((n >> (lvl)) + 1, n >> (lvl));
        rg_cycle[lvl] = new_2darray(n >> (lvl), (n >> (lvl)) + 1);
        rdiv_cycle[lvl] = new_2darray(n >> (lvl), n >> (lvl));
        b_cycle[lvl] = new_vector((n >> (lvl)) + 1);
        t_cycle[lvl] = new_vector((n >> (lvl)) + 1);
        l_cycle[lvl] = new_vector((n >> (lvl)) + 1);
        r_cycle[lvl] = new_vector((n >> (lvl)) + 1);
    }
    for(int i = 0; i<=n ; i++)
    {
        for(int j = 0; j< n; j++)
        {
            f_cycle[0][i][j] = f[2 * i][2 * j + 1];
            g_cycle[0][j][i] = g[2 * j + 1][2 * i];
        }
        b_cycle[0][i] = b[2 * i];
        t_cycle[0][i] = t[2 * i];
        l_cycle[0][i] = l[2 * i];
        r_cycle[0][i] = r[2 * i];
    }
    for(int i = 1; i< n; i++)
    {
        f_cycle[0][i][0] += b[2 * i] * n;
        f_cycle[0][i][n - 1] += t[2 * i] * n;
        g_cycle[0][0][i] += l[2 * i] * n;
        g_cycle[0][n - 1][i] += r[2 * i] * n;
    }

    (*smoother)(n, 0, u, v, p, f_cycle[0], g_cycle[0], b_cycle[0], t_cycle[0], l_cycle[0], r_cycle[0]);
    print(u_cycle[0], n + 1, n, "u_cycle");
    print(v_cycle[0], n, n+1, "v_cycle");
    for(int lvl = 1; lvl < level; lvl ++)//TODO
    {
        residual(n, lvl-1 , u_cycle[lvl - 1], v_cycle[lvl - 1], p_cycle[lvl - 1], 
        f_cycle[lvl - 1], g_cycle[lvl - 1], b, t, l, r, rf_cycle[lvl - 1], rg_cycle[lvl - 1], rdiv_cycle[lvl - 1], &r0, &r0_div);
        restriction(n >> (lvl), lvl, f_cycle[lvl], g_cycle[lvl], p_cycle[lvl], rf_cycle[lvl - 1], rg_cycle[lvl - 1], p_cycle[lvl - 1]);
        for(int i = 0; i<iter; i++)
        {
            (*smoother)(n >> (lvl), lvl, u_cycle[lvl], v_cycle[lvl], p_cycle[lvl], f_cycle[lvl], g_cycle[lvl], b, t, l, r);
        }
        print(u_cycle[lvl], (n >> (lvl)) + 1, n >> (lvl), "u_cycle");
    }
    for(int lvl = level - 1; lvl > 0; lvl --) 
    {
        prolongation(n >> (lvl), lvl, u_cycle[lvl], v_cycle[lvl], p_cycle[lvl], u_cycle_pro[lvl - 1], v_cycle_pro[lvl - 1], p_cycle_pro[lvl - 1]);
        correction(n >> (lvl -1), lvl-1, u_cycle[lvl - 1], v_cycle[lvl - 1], p_cycle[lvl - 1], u_cycle_pro[lvl - 1], v_cycle_pro[lvl - 1], p_cycle_pro[lvl - 1]);
    }

    for(int lvl = 0; lvl< level; lvl++)
    {
        if (lvl>0)
        {
            delete[] u_cycle[lvl];
            delete[] v_cycle[lvl];
            delete[] p_cycle[lvl];
        }
        delete[] u_cycle_pro[lvl];
        delete[] v_cycle_pro[lvl];
        delete[] p_cycle_pro[lvl];
        delete[] f_cycle[lvl] ;
        delete[] g_cycle[lvl] ;
        delete[] rf_cycle[lvl];
        delete[] rg_cycle[lvl];
        delete[] rdiv_cycle[lvl] ;
        delete[] b_cycle[lvl];
        delete[] t_cycle[lvl];
        delete[] l_cycle[lvl];
        delete[] r_cycle[lvl];
    }

    delete[] res ;
    delete[] u_cycle ;
    delete[] v_cycle ;
    delete[] p_cycle ;
    delete[] u_cycle_pro ;
    delete[] v_cycle_pro ;
    delete[] p_cycle_pro ;
    delete[] f_cycle ;
    delete[] g_cycle ;
    delete[] rf_cycle ;
    delete[] rg_cycle ;
    delete[] rdiv_cycle ;
    delete[] b_cycle ;
    delete[] t_cycle ;
    delete[] l_cycle ;
    delete[] r_cycle ;
}