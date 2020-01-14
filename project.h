typedef double dtype;
typedef int smoother_type(int n, dtype ** u, dtype ** v, dtype ** p, dtype ** f, dtype **g, dtype * b, dtype * t, dtype * l, dtype * r, dtype ** res) ;

int residual(int n, int level, dtype ** u, dtype ** v, dtype ** p, dtype ** f, dtype **g, dtype * b, dtype * t, dtype * l, dtype * r, dtype ** rf, dtype ** rg, dtype ** rdiv, dtype * r0, dtype * r1);
int vcycle(int n, int level, int iter, smoother_type smoother, dtype ** u, dtype ** v, dtype ** p, dtype ** f, dtype **g, dtype * b, dtype * t, dtype * l, dtype * r);
int dgs_iteration(int n, dtype ** u, dtype ** v, dtype ** p, dtype ** f, dtype **g, dtype * b, dtype * t, dtype * l, dtype * r, dtype ** res);
int initproblem(int n, dtype ** u, dtype ** v, dtype ** p, dtype ** f, dtype **g, dtype * b, dtype * t, dtype * l, dtype * r, dtype ** u_exact, dtype ** v_exact);
int prolongation(int n, int level, dtype ** ur, dtype ** vr, dtype ** pr, dtype ** u, dtype ** v, dtype ** p);
int restriction(int n, int level, dtype ** ur, dtype ** vr, dtype ** pr, dtype ** u, dtype ** v, dtype ** p);
int error(int n, dtype ** u, dtype ** v, dtype ** u_exact, dtype ** v_exact, dtype * en);


dtype **new_2darray(int n, int m);
dtype *new_vector(int n);

#define checkf(level, i, j) ((f)[i << level][(j << level) + (1 << (level -1))])
#define checkg(level, i, j) ((g)[(i << level) + (1 << (level -1))][j << level])