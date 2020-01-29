#include "project.h"
#include <iostream>
using namespace std;
int vcycle(int n, int level, int mgiter, int mgv0, int mgv1, int maxcnt, smoother_type smoother, dtype ** u, dtype ** v, dtype ** p, dtype ** f, dtype **g, dtype * b, dtype * t, dtype * l, dtype * r)
{
    dtype *** u_cycle = (dtype ***) malloc(level * sizeof(dtype **));
    dtype *** v_cycle = (dtype ***) malloc(level * sizeof(dtype **));
    dtype *** p_cycle = (dtype ***) malloc(level * sizeof(dtype **));
    dtype *** d_cycle = (dtype ***) malloc(level * sizeof(dtype **));
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
    dtype res, res_div;

    u_cycle[0] = u;
    v_cycle[0] = v;
    p_cycle[0] = p;
    
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
        d_cycle[lvl] = new_2darray(n >> (lvl), n >> (lvl));
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

    residual(n, 0, u, v, p, f_cycle[0], g_cycle[0], d_cycle[0], b_cycle[0], t_cycle[0], l_cycle[0], r_cycle[0], rf_cycle[0], rg_cycle[0], rdiv_cycle[0], &r0, &r0_div);
    cout<<"initial residual: "<<r0<<" "<<r0_div<<endl;

    res = r0;
    res_div = r0_div;
    int cnt = 0;
    
    print(f_cycle[0], n + 1, n, "cycle[0] f");
    print(g_cycle[0], n, n + 1, "cycle[0] g");
    while (res / r0 > 1e-8 && cnt < maxcnt)
    {
        for(int lvl = 0; lvl< level; lvl++)
        {
            if (lvl>0)
            {
                clear(u_cycle[lvl], (n >> (lvl)) + 1, n >> (lvl));
                clear(v_cycle[lvl], n >> (lvl), (n >> (lvl)) + 1);
                clear(p_cycle[lvl], n >> (lvl), n >> (lvl));
            }
        }

        residual(n  , 0 , u_cycle[0], v_cycle[0], p_cycle[0], f_cycle[0], g_cycle[0], d_cycle[0], b, t, l, r, rf_cycle[0], rg_cycle[0], rdiv_cycle[0], &res, &res_div);
        print(rf_cycle[0], (n+1), n , "rf_before_cycle");
        print(rg_cycle[0], (n), n+1 , "rg_before_cycle");
        for(int i = 0; i<mgv0; i++)
        {
            (*smoother)(n, level, mgiter, mgv0, mgv1, u, v, p, f_cycle[0], g_cycle[0], d_cycle[0], b_cycle[0], t_cycle[0], l_cycle[0], r_cycle[0]);
        }
        print(u_cycle[0], n + 1, n, "u_cycle");
        print(v_cycle[0], n, n+1, "v_cycle");
        for(int lvl = 1; lvl < level; lvl ++)
        {
            residual(n >> (lvl - 1), lvl-1 , u_cycle[lvl - 1], v_cycle[lvl - 1], p_cycle[lvl - 1], f_cycle[lvl - 1], g_cycle[lvl - 1], d_cycle[lvl - 1], b, t, l, r, rf_cycle[lvl - 1], rg_cycle[lvl - 1], rdiv_cycle[lvl - 1], &res, &res_div);
            restriction(n >> (lvl), lvl, f_cycle[lvl], g_cycle[lvl], d_cycle[lvl], rf_cycle[lvl - 1], rg_cycle[lvl - 1], rdiv_cycle[lvl - 1]);
            print(rf_cycle[lvl-1], (n >>(lvl-1)) + 1, n >> (lvl-1), "rf_cycle");
            print(rg_cycle[lvl-1], (n >>(lvl-1)), (n >> (lvl-1)) + 1, "rg_cycle");
            print(f_cycle[lvl], (n >>(lvl)) + 1, n >> (lvl), "f_cycle");
            print(g_cycle[lvl], (n >>(lvl)), (n >> (lvl)) + 1, "g_cycle");
            for(int i = 0; i<mgv0; i++)
            {
                (*smoother)(n >> (lvl), level, mgiter, mgv0, mgv1, u_cycle[lvl], v_cycle[lvl], p_cycle[lvl], f_cycle[lvl], g_cycle[lvl], d_cycle[lvl], b, t, l, r);
            }
            residual(n >> (lvl), lvl , u_cycle[lvl], v_cycle[lvl], p_cycle[lvl], f_cycle[lvl], g_cycle[lvl], d_cycle[lvl], b, t, l, r, rf_cycle[lvl], rg_cycle[lvl], rdiv_cycle[lvl], &res, &res_div);
            print(rf_cycle[lvl], (n >>(lvl)) + 1, n >> (lvl), "rf_cycle_rough");
            print(rg_cycle[lvl], (n >>(lvl)), (n >> (lvl)) + 1, "rg_cycle_rough"); 
            print(u_cycle[lvl], (n >> (lvl)) + 1, n >> (lvl), "u_cycle");
        }
        for(int lvl = level - 1; lvl > 0; lvl --) 
        {
            prolongation(n >> (lvl), lvl, u_cycle[lvl], v_cycle[lvl], p_cycle[lvl], u_cycle_pro[lvl - 1], v_cycle_pro[lvl - 1], p_cycle_pro[lvl - 1]);
            print(u_cycle_pro[lvl-1], (n >> (lvl-1)) + 1, n >> (lvl-1), "u_cycle_pro");
            print(v_cycle_pro[lvl-1], (n >> (lvl-1)) , (n >> (lvl-1))+1, "v_cycle_pro");
            print(u_cycle[lvl-1], (n >> (lvl-1)) + 1, n >> (lvl-1), "u_cycle_bef");
            print(v_cycle[lvl-1], (n >> (lvl-1)) , (n >> (lvl-1))+1, "v_cycle_bef");
            correction(n >> (lvl -1), lvl-1, u_cycle[lvl - 1], v_cycle[lvl - 1], p_cycle[lvl - 1], u_cycle_pro[lvl - 1], v_cycle_pro[lvl - 1], p_cycle_pro[lvl - 1]);
            print(u_cycle[lvl-1], (n >> (lvl-1)) + 1, n >> (lvl-1), "u_cycle_correct");
            print(v_cycle[lvl-1], (n >> (lvl-1)) , (n >> (lvl-1))+1, "v_cycle_correct");
            for(int i = 0; i<mgv1; i++)
            {
                (*smoother)(n >> (lvl - 1), level, mgiter, mgv0, mgv1, u_cycle[lvl - 1], v_cycle[lvl - 1], p_cycle[lvl - 1], f_cycle[lvl - 1], g_cycle[lvl - 1], d_cycle[lvl - 1], b, t, l, r);
            }
            
        }

        residual(n  , 0 , u_cycle[0], v_cycle[0], p_cycle[0], f_cycle[0], g_cycle[0], d_cycle[0], b, t, l, r, rf_cycle[0], rg_cycle[0], rdiv_cycle[0], &res, &res_div);
        print(rf_cycle[0], (n+1), n , "rf_after_cycle");
        print(rg_cycle[0], (n), n+1 , "rg_after_cycle");


        cout<<"iteration "<<cnt<<", residual "<<res<<", "<<res_div<<endl;
        cnt ++;
    }
    
    for(int lvl = 0; lvl< level; lvl++)
    {
        if (lvl>0)
        {
            delete_2darray(u_cycle[lvl], (n >> (lvl)) + 1);
            delete_2darray(v_cycle[lvl], n >> (lvl));
            delete_2darray(p_cycle[lvl], n >> (lvl));
        }
        delete_2darray(u_cycle_pro[lvl], (n >> (lvl)) + 1);
        delete_2darray(v_cycle_pro[lvl], n >> (lvl));
        delete_2darray(p_cycle_pro[lvl], n >> (lvl));
        delete_2darray(d_cycle[lvl], n >> (lvl));
        delete_2darray(f_cycle[lvl], (n >> (lvl)) + 1);
        delete_2darray(g_cycle[lvl], n >> (lvl));
        delete_2darray(rf_cycle[lvl], (n >> (lvl)) + 1);
        delete_2darray(rg_cycle[lvl], n >> (lvl));
        delete_2darray(rdiv_cycle[lvl], n >> (lvl));
        free(b_cycle[lvl]);
        free(t_cycle[lvl]);
        free(l_cycle[lvl]);
        free(r_cycle[lvl]);
    }

    free(u_cycle);
    free(v_cycle);
    free(p_cycle);
    free(u_cycle_pro);
    free(v_cycle_pro);
    free(p_cycle_pro);
    free(f_cycle);
    free(g_cycle);
    free(d_cycle);
    free(rf_cycle);
    free(rg_cycle);
    free(rdiv_cycle);
    free(b_cycle);
    free(t_cycle);
    free(l_cycle);
    free(r_cycle);
    return cnt;
}