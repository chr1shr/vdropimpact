#ifndef TRI_SOLVE_HH
#define TRI_SOLVE_HH

/** A class for representing and solving a tridiagonal linear system using
 * LAPACK. */
class tri_solve {
    public:
        /** The dimension of the system. */
        int m;
        /** An array of lower-diagonal terms in the matrix. */
        double* const a;
        /** An array of diagonal terms in the matrix. */
        double* const b;
        /** An array of upper-diagonal terms in the matrix. */
        double* const c;
        /** The source vector. */
        double* const d;
        tri_solve(int m_);
        ~tri_solve();
        void solve();
        void solve_lapack();
};

#endif
