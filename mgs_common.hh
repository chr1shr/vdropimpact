#ifndef VDI_MGS_COMMON_HH
#define VDI_MGS_COMMON_HH

#include "tgmg.hh"

class fluid_2d;

/** \brief A base class describing common constants in all linear systems to
 * be solved by the multigrid library. */
struct mgs_common {
    /** The number of gridpoints in the x direction. */
    const int m;
    /** The number of gridpoints in the y direction. */
    const int n;
    /** Total number of gridpoints. */
    const int mn;
    /** The total number of multigrid solves that have been performed. This is
     * used by the initial guess function, since the code only starts using
     * linear interpolation after two solves have been performed so that there
     * are enough solutions to interpolate from. */
    int scount;
    /** Periodicity in the x direction. */
    static const bool x_prd=false;
    /** Periodicity in the y direction. */
    static const bool y_prd=false;
    /** The mode to use for the Gauss-Seidel smoothing. (0=default) */
    static const char gs_mode=0;
    /** The simulation time corresponding to the solution stored in the z_prev
     * array. */
    double t_prev;
    /** The simulation time corresponding to the solution in the z array. */
    double t0;
    /** An array for storing a previous solution, which is used to initialize a
     * better starting guess for the multigrid solve, using linear
     * interpolation. */
    double* const z_prev;
    /** An array for computing the solution. */
    double* const z;
    mgs_common(int m_,int n_);
    /** The class destructor frees the dynamically allocated memory that was
     * allocated for the solutions. */
    ~mgs_common() {
        delete [] z;
        delete [] z_prev;
    }
    void initial_guess(double t);
    inline void reset_counters() {
        wc_time=0;
        tp.avg_iters(true);
    }
    /** The wall clock time counter. */
    double wc_time;
    /** The multigrid helper class that is used for predicting how many
     * V-cycles to do before testing for convergence. */
    tgmg_predict tp;
};

#endif
