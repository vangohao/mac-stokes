#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "project.h"
using namespace std;
dtype **new_2darray(int n, int m)
{
    dtype ** array = (dtype **)malloc(sizeof(dtype*) * n);
    for(int i = 0; i < n; i++)
    {
        array[i] = (dtype*)calloc(m, sizeof(dtype));
    }
    return array;
}
dtype *new_vector(int n)
{
    return (dtype*)calloc(n, sizeof(dtype));
}
void print(dtype ** a, int n, int m, const char * name)
{
    // cout<<"start printing "<<name<<":"<<endl;
    // for(int i = 0; i < n; i++)
    // {
    //     for(int j = 0; j < m; j++)
    //     {
    //         cout<<a[i][j]<<" ";
    //     }
    //     cout<<endl;
    // }
}
void clear(dtype ** a, int n, int m)
{
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < m; j++)
        {
            a[i][j] = 0;
        }
    }
}