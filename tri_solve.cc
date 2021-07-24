#include "tri_solve.hh"

#include <cstdio>
#include <cstdlib>

// Tell the compiler about the existence of LAPACK's tridiagonal solver
extern "C" {
    void dgtsv_(int *n,int *nrhs,double *dl,double *d,double *du,double *b,int *ldb,int *info);
}

/** The tridiagonal solver constructor allocates memory for the tridiagonal matrix
 * coefficients and the source term.
 * \param[in] m_ the dimension of the source term. */
tri_solve::tri_solve(int m_) : m(m_), a(new double[4*m]),
    b(a+m), c(b+m), d(c+m) {}

/** The destructor frees the dynamically allocated memory. */
tri_solve::~tri_solve() {
    delete [] a;
}

/** Solves the tridiagonal system using LAPACK's routine, placing the solution
 * in the d array. This routine also overwrites the values of a, b, and c. */
void tri_solve::solve_lapack() {
    int info,nrhs=1;

    // Call the LAPACK routine and check for successful completion
    dgtsv_(&m,&nrhs,a+1,b,c,d,&m,&info);
    if(info!=0) {
        fputs("LAPACK routine failed\n",stderr);
        exit(1);
    }
}

/** Solves the tridiagonal system using the Thomas algorithm, placing the
 * answer in the d array. */
void tri_solve::solve() {

    // Forward sweep
    *c/=*b;
    for(int i=1;i<m-1;i++) c[i]/=b[i]-a[i]*c[i-1];

    // Backward sweep
    *d/=*b;
    for(int i=1;i<m;i++) d[i]=(d[i]-a[i]*d[i-1])/(b[i]-a[i]*c[i-1]);

    // Solution
    for(int i=m-2;i>=0;i--) d[i]-=c[i]*d[i+1];
}
