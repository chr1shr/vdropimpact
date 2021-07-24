#ifndef VDI_MGS_MAC_HH
#define VDI_MGS_MAC_HH

#include <cstdio>
#include <cstdlib>

#include "common.hh"
#include "tgmg.hh"
#include "mgs_common.hh"

class fluid_2d;

struct mgs_mac : public mgs_common {
    /** The inverse grid spacing squared in the x direction. */
    const double xxsp;
    /** The inverse grid spacing squared in the y direction. */
    const double yysp;
    /** Pre-computed reciprocals of diagonal matrix entries. */
    const double cc0,cc1,cc2,cc3,cc4,cc5;
    /** Threshold on L_2 norm of residual to terminate the multigrid solve. */
    double acc;
    mgs_mac(fluid_2d &f);
    ~mgs_mac() {}
    mgs_mac(int m_,int n_,double dx,double dy,double *b);
    inline double a_dl(int i,int ij) {return 0;}
    inline double a_dr(int i,int ij) {return 0;}
    inline double a_ul(int i,int ij) {return 0;}
    inline double a_ur(int i,int ij) {return 0;}
    inline double a_dc(int i,int ij) {return ij>=m?yysp:0;}
    inline double a_uc(int i,int ij) {return ij<mn-m?yysp:0;}
    inline double a_cl(int i,int ij) {return i>0?xxsp:0;}
    inline double a_cr(int i,int ij) {return i<m-1?xxsp:0;}
    inline double a_cc(int i,int ij) {
        return (i>0&&i<m-1?-2*xxsp:-xxsp)
               +(ij<m?-3*yysp:ij<mn-m?-2*yysp:-yysp);
    }
    inline double inv_cc(int i,int ij,double v) {
        return v*(i>0&&i<m-1?(ij<m?cc0:ij<mn-m?cc1:cc2)
                            :(ij<m?cc3:ij<mn-m?cc4:cc5));
    }
    /** Performs multiplication by the linear system, excluding the diagonal
     * term.
     * \param[in] i the horizontal grid index.
     * \param[in] ij the overall grid index.
     * \return The result of the multiplication. */
    inline double mul_a(int i,int ij) {
        double *w=z+ij;
        return xxsp*(i==0?w[1]:i<m-1?w[-1]+w[1]:w[-1])
               +yysp*(ij<m?w[m]:ij<mn-m?w[-m]+w[m]:w[-m]);
    }
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
        fatal_error("Multigrid routine for MAC projection failed",1);
    }
    /** The multigrid solver. */
    tgmg<mgs_mac,double,double> mg;
};

#endif
