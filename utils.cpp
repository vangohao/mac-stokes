#include <stdio.h>
#include <stdlib.h>
#include "project.h"
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