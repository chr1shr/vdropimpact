#ifndef VDI_VISCO_IMPL_HH
#define VDI_VISCO_IMPL_HH

#include <cstdlib>
#include <cmath>
#include <limits>

#include "tgmg.hh"
#include "mgs_common.hh"

class fluid_2d;

/** \brief Class describing the linear system arising from the the liquid
 * viscous term using an implicit Crank--Nicolson discretization. */
struct visco_impl : public mgs_common {
    /** Symmetric about x = 0. */
    const bool x_sym;
    /** Variable to differentiate between linear systems for horizontal and
     * vertical velocities. */
    const bool horizontal_vel;
    /** Threshold on L_2 norm of residual to terminate the multigrid solve. */
    const double acc;
    /** The viscosity factors. */
    double vfx,vfy;
    visco_impl(fluid_2d &f,bool horizontal_vel_);
    ~visco_impl() {}
    void setup_and_solve(double vfx_,double vfy_,double t);
    inline double a_dl(int i,int ij) {return 0.;}
    inline double a_dr(int i,int ij) {return 0.;}
    inline double a_ul(int i,int ij) {return 0.;}
    inline double a_ur(int i,int ij) {return 0.;}
    inline double a_dc(int i,int ij) {return ij>=m?-vfy:0.;}
    inline double a_uc(int i,int ij) {return ij<mn-m?-vfy:0.;}
    inline double a_cl(int i,int ij) {return i>0?-vfx:0.;}
    inline double a_cr(int i,int ij) {return i<m-1?-vfx:0.;}
    inline double a_cc(int i,int ij) {return x_sym?
        1+(i>0?2*vfx:(horizontal_vel?3*vfx:vfx))+(ij>=m?2*vfy:vfy)
        :1+2*vfx+(ij>=m?2*vfy:vfy);}
    inline double inv_cc(int i,int ij,double v) {return v/a_cc(i,ij);}
    /** Calculates the result of multiplying the matrix A by the solution
     * vector at one grid point. The diagonal term in the matrix A is omitted.
     * \param[in] i the horizontal index of the grid point.
     * \param[in] ij the overall index of the grid point.
     * \return The result of the multiplication. */
    inline double mul_a(int i,int ij) {
        double *w=z+ij;
        return -vfx*(i>0?(i<m-1?w[1]+w[-1]
                               :w[-1])
                        :w[1])
               -vfy*(ij>=m?(ij<mn-m?w[m]+w[-m]
                                   :w[-m])
                          :w[m]);
    }
    /** The multigrid solver. */
    tgmg<visco_impl,double,double> mg;
};

#endif
