typedef double dtype;

typedef int smoother_type(int n, int level, int mgiter, int mgv0, int mgv1, dtype ** u, dtype ** v, dtype ** p, dtype ** f, dtype **g, dtype ** d, dtype * b, dtype * t, dtype * l, dtype * r) ;
int residual(int n, int level, dtype ** u, dtype ** v, dtype ** p, dtype ** f, dtype **g, dtype ** d, dtype * b, dtype * t, dtype * l, dtype * r, dtype ** rf, dtype ** rg, dtype ** rdiv, dtype * r0, dtype * r1);
int vcycle(int n, int level, int mgiter, int mgv0, int mgv1, int maxcnt, smoother_type smoother, dtype ** u, dtype ** v, dtype ** p, dtype ** f, dtype **g, dtype * b, dtype * t, dtype * l, dtype * r);
int dgs_iteration(int n, int level, int mgiter, int mgv0, int mgv1, dtype ** u, dtype ** v, dtype ** p, dtype ** f, dtype **g, dtype ** d, dtype * b, dtype * t, dtype * l, dtype * r);
int initproblem(int n, dtype ** u, dtype ** v, dtype ** p, dtype ** f, dtype **g, dtype * b, dtype * t, dtype * l, dtype * r, dtype ** u_exact, dtype ** v_exact);
int prolongation(int n, int level, dtype ** ur, dtype ** vr, dtype ** pr, dtype ** u, dtype ** v, dtype ** p);
int restriction(int n, int level, dtype ** ur, dtype ** vr, dtype ** pr, dtype ** u, dtype ** v, dtype ** p);
int error(int n, dtype ** u, dtype ** v, dtype ** u_exact, dtype ** v_exact, dtype * en);
int correction(int n, int level, dtype ** ur, dtype ** vr, dtype ** pr, dtype ** u, dtype ** v, dtype ** p);
int cg(int n, int level, int kmax, dtype epsb, dtype **u, dtype **v, dtype **rf, dtype **rg);
int uzawa_iteration(int n, int level, int mgiter, int mgv0, int mgv1, dtype ** u, dtype ** v, dtype ** p, dtype ** f, dtype **g, dtype ** d, dtype * b, dtype * t, dtype * l, dtype * r);
int inexact_uzawa_iteration(int n, int level, int mgiter, int mgv0, int mgv1, dtype ** u, dtype ** v, dtype ** p, dtype ** f, dtype **g, dtype ** d, dtype * b, dtype * t, dtype * l, dtype * r);
int pcg(int n, int level, int mgiter, int mgv0, int mgv1, int kmax, dtype eps, dtype **u, dtype **v, dtype **bf, dtype **bg);
int residual_uv(int n, int level, dtype ** u, dtype ** v, dtype ** f, dtype **g, dtype ** rf, dtype ** rg, dtype * r0);
int prolongation_uv(int n, int level, dtype ** ur, dtype ** vr, dtype ** u, dtype ** v);
int restriction_uv(int n, int level, dtype ** ur, dtype ** vr, dtype ** u, dtype ** v);
int correction_uv(int n, int level, dtype ** ur, dtype ** vr, dtype ** u, dtype ** v);
int GS_uv(int n, int level, dtype ** u, dtype ** v, dtype ** f, dtype **g);
int vcycle_precondition(int n, int level, int mgiter, int mgv0, int mgv1, dtype ** u, dtype ** v, dtype ** f, dtype **g);
int pcg_uzawa_iteration(int n, int level, int mgiter, int mgv0, int mgv1, dtype ** u, dtype ** v, dtype ** p, dtype ** f, dtype **g, dtype ** d, dtype * b, dtype * t, dtype * l, dtype * r);

dtype **new_2darray(int n, int m);
dtype *new_vector(int n);
void delete_2darray(dtype ** a, int n);
void print(dtype ** a, int n, int m, const char * name);
void clear(dtype ** a, int n, int m);

// #define checkf(level, i, j) ((f)[i << (level + 1)][(j << (level + 1)) + (1 << (level))])
// #define checkg(level, i, j) ((g)[(i << (level+1)) + (1 << (level))][j << level])
#define checkf(level,i,j)  (f[i][j])
#define checkg(level,i,j)  (g[i][j])