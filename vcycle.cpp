#include "project.h"
#include <iostream>
using namespace std;
int vcycle(int n, int level, int iter, smoother_type smoother, dtype ** u, dtype ** v, dtype ** p, dtype ** f, dtype **g, dtype * b, dtype * t, dtype * l, dtype * r)
{
    dtype** res = new_2darray(n, n);
    (*smoother)(n, u, v, p, f, g, b, t, l, r, res);
}