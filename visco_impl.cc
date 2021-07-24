#include <cstdio>
#include <limits>

#include "visco_impl.hh"
#include "fluid_2d.hh"
#include "tgmg.hh"

/** The constructor sets many internal constants from the parent fluid_2d
 * class.
 * \param[in] f a reference to the parent fluid_2d class. */
visco_impl::visco_impl(fluid_2d &f,bool horizontal_vel_) : mgs_common(f.m,f.n),
    x_sym(f.x_sym), horizontal_vel(horizontal_vel_),
    acc(tgmg_accuracy(double(1.),1e3)), vfx(0.), mg(*this,f.src,z) {}

/** Sets up the linear system and solves it using the multigrid method.
 * \param[in] (vfx,vfy) prefactors in the linear system in the x and y
 *                      directions, which depend on the simulation timestep
 *                      used.
 * \param[in] t the current simulation time, used to construct an initial guess
 *              for the multigrid method via linear interpolation in time. */
void visco_impl::setup_and_solve(double vfx_,double vfy_,double t) {

    // If the timestep has changed from the previous call, then set up the
    // multigrid hierarchy again
    double t0=wtime();
    if(fabs(vfx_-vfx)>1e2*std::numeric_limits<double>::epsilon()*vfx) {
        vfx=vfx_;
        vfy=vfy_;
        mg.setup();
    }

    // Create initial guess using time interpolation, and solve using V-cycles
    initial_guess(t);
    if(!mg.solve_v_cycle(tp)) {
        fputs("Multigrid routine for implicit viscosity failed\n",stderr);
        exit(1);
    }
    wc_time+=wtime()-t0;
}

// Explicit instantiation
#include "tgmg.cc"
template class tgmg<visco_impl,double,double>;
template void tgmg_base<visco_impl,double,double>::output(char const*,double*,double,double,double,double);
template void tgmg_base<visco_impl,double,double>::clear_z();
template void tgmg_base<visco_impl,double,double>::output_res(char const*,double,double,double,double);
