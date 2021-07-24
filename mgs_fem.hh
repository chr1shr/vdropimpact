#ifndef VDI_MGS_FEM_HH
#define VDI_MGS_FEM_HH

#include <cstdio>
#include <cstdlib>

#include "common.hh"
#include "tgmg.hh"
#include "mgs_common.hh"

class fluid_2d;

/** \brief A class describing the linear system arising in the finite-element
 * pressure projection. */
struct mgs_fem : public mgs_common {
    /** The ratio of vertical and horizontal grid spacings, which appears in
     * the finite-element stencils. */
    const double dydx;
    /** The ratio of horizontal and vertical grid spacings, which appears in
     * the finite-element stencils. */
    const double dxdy;
    /** The central term in the finite-element stencil. */
    const double fm;
    /** The reciprocal of the central term in the finite-element stencil, which
     * is stored separately because it is used frequently. */
    const double fm_inv;
    /** The vertical term in the finite-element stencil. */
    const double fey;
    /** Half of the vertical term in the finite-element stencil. */
    const double hey;
    /** The horizontal term in the finite-element stencil. */
    const double fex;
    /** Half of the horizontal term in the finite-element stencil. */
    const double hex;
    /** The corner term in the finite-element stencil. */
    const double fc;
    /** Threshold on L_2 norm of residual to terminate the multigrid solve. */
    const double acc;
    mgs_fem(fluid_2d &f);
    ~mgs_fem() {}
    mgs_fem(int m_,int n_,double dx,double dy,double *b);
    inline bool bottom(int ij) {return ij<m;}
    inline double a_dl(int i,int ij) {return !bottom(ij)&&ij>=m&&i>0?fc:0;}
    inline double a_dr(int i,int ij) {return !bottom(ij)&&ij>=m&&i<m-1?fc:0;}
    inline double a_ul(int i,int ij) {return !bottom(ij)&&ij<mn-m&&i>0?fc:0;}
    inline double a_ur(int i,int ij) {return !bottom(ij)&&ij<mn-m&&i<m-1?fc:0;}
    inline double a_dc(int i,int ij) {return !bottom(ij)&&ij>=m?(i>0&&i<m-1?fey:hey):0;}
    inline double a_uc(int i,int ij) {return !bottom(ij)&&ij<mn-m?(i>0&&i<m-1?fey:hey):0;}
    inline double a_cl(int i,int ij) {return !bottom(ij)&&i>0?(ij>=m&&ij<mn-m?fex:hex):0;}
    inline double a_cr(int i,int ij) {return !bottom(ij)&&i<m-1?(ij>=m&&ij<mn-m?fex:hex):0;}
    inline double a_cc(int i,int ij) {
        return (i>0&&i<m-1?1:0.5)*(ij>=m&&ij<mn-m?1:0.5)*fm;
    }
    inline double inv_cc(int i,int ij,double v) {
        return (i>0&&i<m-1?1:2)*(ij>=m&&ij<mn-m?1:2)*fm_inv*v;
    }
    double mul_a(int i,int ij);
        /** Solves the linear system using multigrid V-cycles. The initial guess is
     * based on linear interpolation in time from the previous solutions.
     * \param[in] t the current simulation time. */
    inline void solve_v_cycle(double t) {
        initial_guess(t);
        double t0=wtime();
        if(!mg.solve_v_cycle(tp)) failure();
        wc_time+=wtime()-t0;
    }
    /** Solves the linear system using multigrid V-cycles. */
    inline void solve_v_cycle() {
        double t0=wtime();
        if(!mg.solve_v_cycle(tp)) failure();
        wc_time+=wtime()-t0;
    }
    /** Prints an error message in the case when the multigrid method failed to
     * converge. */
    inline void failure() {
        fatal_error("Multigrid routine for finite-element problem failed",1);
    }
    /** The multigrid solver. */
    tgmg<mgs_fem,double,double> mg;
};

#endif
