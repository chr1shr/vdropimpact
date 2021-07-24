#include <cstdio>
#include <cstdlib>

#include "tri_solve.hh"

int main() {

    // Number of grid points
    const int m=32;

    // Grid spacing constants
    const double dx=1./(m+1.),xsp=1/dx,xxsp=xsp*xsp;

    // Initialize the matrix coefficients for a Poisson problem on [0,1]
    // with zero Dirichlet BCs. (Note t.a[0] and t.c[n-1] will be ignored.)
    tri_solve t(m);
    for(int i=0;i<m;i++) {
        t.a[i]=t.c[i]=xxsp;
        t.b[i]=-2*xxsp;
        t.d[i]=1;
    }

    // Solve the tridiagonal matrix problem
    t.solve();

    // Print the solution
    puts("0 0");
    for(int i=0;i<m;i++) printf("%g %g\n",(i+1)*dx,t.d[i]);
    puts("1 0");
}
