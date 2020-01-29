#include "project.h"
#include <iostream>
using namespace std;
int vcycle_precondition(int n, int level, int mgiter, int mgv0, int mgv1, dtype ** u, dtype ** v, dtype ** f, dtype **g)
{
    cout<<"VCYCLE_PRECONDITION"<<endl;
    dtype *** u_cycle = (dtype ***) malloc(level * sizeof(dtype **));
    dtype *** v_cycle = (dtype ***) malloc(level * sizeof(dtype **));
    dtype *** u_cycle_pro = (dtype ***) malloc(level * sizeof(dtype **));
    dtype *** v_cycle_pro = (dtype ***) malloc(level * sizeof(dtype **));
    dtype *** f_cycle = (dtype ***) malloc(level * sizeof(dtype **));
    dtype *** g_cycle = (dtype ***) malloc(level * sizeof(dtype **));
    dtype *** rf_cycle = (dtype ***) malloc(level * sizeof(dtype **));
    dtype *** rg_cycle = (dtype ***) malloc(level * sizeof(dtype **));

    dtype r0;
    dtype res;

    u_cycle[0] = u;
    v_cycle[0] = v;
    
    for(int lvl = 0; lvl< level; lvl++)
    {
        if (lvl>0)
        {
            u_cycle[lvl] = new_2darray((n >> (lvl)) + 1, n >> (lvl));
            v_cycle[lvl] = new_2darray(n >> (lvl), (n >> (lvl)) + 1);
        }
        u_cycle_pro[lvl] = new_2darray((n >> (lvl)) + 1, n >> (lvl));
        v_cycle_pro[lvl] = new_2darray(n >> (lvl), (n >> (lvl)) + 1);
        f_cycle[lvl] = new_2darray((n >> (lvl)) + 1, n >> (lvl));
        g_cycle[lvl] = new_2darray(n >> (lvl), (n >> (lvl)) + 1);
        rf_cycle[lvl] = new_2darray((n >> (lvl)) + 1, n >> (lvl));
        rg_cycle[lvl] = new_2darray(n >> (lvl), (n >> (lvl)) + 1);
    }
    for(int i = 0; i<=n ; i++)
    {
        for(int j = 0; j< n; j++)
        {
            f_cycle[0][i][j] = f[i][j];
            g_cycle[0][j][i] = g[j][i];
        }
    }

    residual_uv(n, 0, u, v, f_cycle[0], g_cycle[0], rf_cycle[0], rg_cycle[0], &r0);
    // cout<<"precondition initial residual: "<<r0<<endl;

    res = r0;
    int cnt = 0;
    
    print(f_cycle[0], n + 1, n, "cycle[0] f");
    print(g_cycle[0], n, n + 1, "cycle[0] g");
    while (res / r0 > 1e-8 && cnt < mgiter)
    {
        for(int lvl = 0; lvl< level; lvl++)
        {
            if (lvl>0)
            {
                clear(u_cycle[lvl], (n >> (lvl)) + 1, n >> (lvl));
                clear(v_cycle[lvl], n >> (lvl), (n >> (lvl)) + 1);
            }
        }

        residual_uv(n, 0, u_cycle[0], v_cycle[0], f_cycle[0], g_cycle[0], rf_cycle[0], rg_cycle[0], &res);
        print(rf_cycle[0], (n+1), n , "rf_before_cycle");
        print(rg_cycle[0], (n), n+1 , "rg_before_cycle");
        for(int i = 0; i<mgv0; i++)
        {
            GS_uv(n, 0, u, v, f_cycle[0], g_cycle[0]);
        }
        print(u_cycle[0], n + 1, n, "u_cycle");
        print(v_cycle[0], n, n+1, "v_cycle");
        for(int lvl = 1; lvl < level; lvl ++)
        {
            residual_uv(n >> (lvl - 1), lvl-1 , u_cycle[lvl - 1], v_cycle[lvl - 1], f_cycle[lvl - 1], g_cycle[lvl - 1], rf_cycle[lvl - 1], rg_cycle[lvl - 1], &res);
            restriction_uv(n >> (lvl), lvl, f_cycle[lvl], g_cycle[lvl], rf_cycle[lvl - 1], rg_cycle[lvl - 1]);
            print(rf_cycle[lvl-1], (n >>(lvl-1)) + 1, n >> (lvl-1), "rf_cycle");
            print(rg_cycle[lvl-1], (n >>(lvl-1)), (n >> (lvl-1)) + 1, "rg_cycle");
            print(f_cycle[lvl], (n >>(lvl)) + 1, n >> (lvl), "f_cycle");
            print(g_cycle[lvl], (n >>(lvl)), (n >> (lvl)) + 1, "g_cycle");
            for(int i = 0; i<mgv0; i++)
            {
                GS_uv(n >> (lvl), lvl, u_cycle[lvl], v_cycle[lvl], f_cycle[lvl], g_cycle[lvl]);
            }
            residual_uv(n >> (lvl), lvl , u_cycle[lvl], v_cycle[lvl], f_cycle[lvl], g_cycle[lvl], rf_cycle[lvl], rg_cycle[lvl], &res);
            print(rf_cycle[lvl], (n >>(lvl)) + 1, n >> (lvl), "rf_cycle_rough");
            print(rg_cycle[lvl], (n >>(lvl)), (n >> (lvl)) + 1, "rg_cycle_rough"); 
            print(u_cycle[lvl], (n >> (lvl)) + 1, n >> (lvl), "u_cycle");
        }
        for(int lvl = level - 1; lvl > 0; lvl --) 
        {
            prolongation_uv(n >> (lvl), lvl, u_cycle[lvl], v_cycle[lvl], u_cycle_pro[lvl - 1], v_cycle_pro[lvl - 1]);
            print(u_cycle_pro[lvl-1], (n >> (lvl-1)) + 1, n >> (lvl-1), "u_cycle_pro");
            print(v_cycle_pro[lvl-1], (n >> (lvl-1)) , (n >> (lvl-1))+1, "v_cycle_pro");
            print(u_cycle[lvl-1], (n >> (lvl-1)) + 1, n >> (lvl-1), "u_cycle_bef");
            print(v_cycle[lvl-1], (n >> (lvl-1)) , (n >> (lvl-1))+1, "v_cycle_bef");
            correction_uv(n >> (lvl -1), lvl-1, u_cycle[lvl - 1], v_cycle[lvl - 1], u_cycle_pro[lvl - 1], v_cycle_pro[lvl - 1]);
            print(u_cycle[lvl-1], (n >> (lvl-1)) + 1, n >> (lvl-1), "u_cycle_correct");
            print(v_cycle[lvl-1], (n >> (lvl-1)) , (n >> (lvl-1))+1, "v_cycle_correct");
            for(int i = 0; i<mgv1; i++)
            {
                GS_uv(n >> (lvl - 1), lvl - 1, u_cycle[lvl - 1], v_cycle[lvl - 1], f_cycle[lvl - 1], g_cycle[lvl - 1]);
            }
            
        }

        residual_uv(n, 0 , u_cycle[0], v_cycle[0], f_cycle[0], g_cycle[0], rf_cycle[0], rg_cycle[0], &res);
        print(rf_cycle[0], (n+1), n , "rf_after_cycle");
        print(rg_cycle[0], (n), n+1 , "rg_after_cycle");


        // cout<<"precondition iteration "<<cnt<<", residual "<<res<<endl;
        cnt ++;
    }
    
    for(int lvl = 0; lvl< level; lvl++)
    {
        if (lvl>0)
        {
            delete_2darray(u_cycle[lvl], (n >> (lvl)) + 1);
            delete_2darray(v_cycle[lvl], n >> (lvl));
        }
        delete_2darray(u_cycle_pro[lvl], (n >> (lvl)) + 1);
        delete_2darray(v_cycle_pro[lvl], n >> (lvl));
        delete_2darray(f_cycle[lvl], (n >> (lvl)) + 1);
        delete_2darray(g_cycle[lvl], n >> (lvl));
        delete_2darray(rf_cycle[lvl], (n >> (lvl)) + 1);
        delete_2darray(rg_cycle[lvl], n >> (lvl));
    }

    free(u_cycle);
    free(v_cycle);
    free(u_cycle_pro);
    free(v_cycle_pro);
    free(f_cycle);
    free(g_cycle);
    free(rf_cycle);
    free(rg_cycle);
    return cnt;
}