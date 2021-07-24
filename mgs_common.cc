#include "mgs_common.hh"
#include "tgmg.hh"

#include <cstdio>
#include <cstring>
#include <limits>

/** Sets up the geometry constants and initializes the solution array.
 * \param[in] (m_,n_) the dimensions of the grid. */
mgs_common::mgs_common(int m_,int n_) :
    m(m_), n(n_), mn(m*n), scount(0), z_prev(new double[mn]),
    z(new double[mn]), wc_time(0.), tp() {}

/** Computes an initial guess for the multigrid solve, by interpolating in time
 * from previous solutions.
 * \param[in] t the current time. */
void mgs_common::initial_guess(double t) {

    // If fewer than two solutions have been performed, there is not enough
    // information to interpolate in time. Just store the results if available,
    // and move on.
    if(scount<2) {
        if(scount==1) {
            memcpy(z_prev,z,mn*sizeof(double));
            t0=t;
        } else t_prev=t;
        scount++;
    } else {

        // Check for a pathological case when this routine is called for
        // matching times, which make it impossible to perform
        if(fabs(t0-t_prev)<std::numeric_limits<double>::epsilon()*10) {
            fputs("Multigrid routine called for duplicate times\n",stderr);
            return;
        }

        // Interpolate in time
        double a=1./(t0-t_prev),b=(t-t_prev)*a,c=(t-t0)*a;
#pragma omp parallel for
        for(int ij=0;ij<mn;ij++) {
            double zz=z[ij];
            z[ij]=b*zz-c*z_prev[ij];
            z_prev[ij]=zz;
        }

        // Store the times of the previous frames
        t_prev=t0;
        t0=t;
    }
}

// Explicit instantiation
#include "tgmg.cc"
template class tgmg_base<tgmg_level<double,double>,double,double>;
