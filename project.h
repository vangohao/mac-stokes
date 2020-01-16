typedef double dtype;

typedef int smoother_type(int n, int level, dtype ** u, dtype ** v, dtype ** p, dtype ** f, dtype **g, dtype ** d, dtype * b, dtype * t, dtype * l, dtype * r) ;
int residual(int n, int level, dtype ** u, dtype ** v, dtype ** p, dtype ** f, dtype **g, dtype ** d, dtype * b, dtype * t, dtype * l, dtype * r, dtype ** rf, dtype ** rg, dtype ** rdiv, dtype * r0, dtype * r1);
int vcycle(int n, int level, int iter, smoother_type smoother, dtype ** u, dtype ** v, dtype ** p, dtype ** f, dtype **g, dtype * b, dtype * t, dtype * l, dtype * r);
int dgs_iteration(int n, int level, dtype ** u, dtype ** v, dtype ** p, dtype ** f, dtype **g, dtype ** d, dtype * b, dtype * t, dtype * l, dtype * r);
int initproblem(int n, dtype ** u, dtype ** v, dtype ** p, dtype ** f, dtype **g, dtype * b, dtype * t, dtype * l, dtype * r, dtype ** u_exact, dtype ** v_exact);
int prolongation(int n, int level, dtype ** ur, dtype ** vr, dtype ** pr, dtype ** u, dtype ** v, dtype ** p);
int restriction(int n, int level, dtype ** ur, dtype ** vr, dtype ** pr, dtype ** u, dtype ** v, dtype ** p);
int error(int n, dtype ** u, dtype ** v, dtype ** u_exact, dtype ** v_exact, dtype * en);
int correction(int n, int level, dtype ** ur, dtype ** vr, dtype ** pr, dtype ** u, dtype ** v, dtype ** p);


dtype **new_2darray(int n, int m);
dtype *new_vector(int n);
void print(dtype ** a, int n, int m, const char * name);

// #define checkf(level, i, j) ((f)[i << (level + 1)][(j << (level + 1)) + (1 << (level))])
// #define checkg(level, i, j) ((g)[(i << (level+1)) + (1 << (level))][j << level])
#define checkf(level,i,j)  (f[i][j])
#define checkg(level,i,j)  (g[i][j])