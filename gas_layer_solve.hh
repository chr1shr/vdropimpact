#ifndef VDI_GAS_LAYER_SOLVE_HH
#define VDI_GAS_LAYER_SOLVE_HH

#include <limits>

#include "fields.hh"
#include "tri_solve.hh"
#include "gas_layer.hh"

/** A scale for rounding error for the nonlinear pressure problem. */
const double pressure_r_error=std::numeric_limits<double>::epsilon()*1e8;

/** Threshold on the mean squared residual per gridpoint for the Newton solver.
 */
const double newton_thresh=pressure_r_error*pressure_r_error;

/** The maximum number of Newton steps to try before bailing out. */
const int max_newton_steps=256;

/** \brief Structure describing the physical parameters for the gas layer. */
struct gl_params {
    /** The reciprocal exponent in the relationship linking density to pressure. */
    const double alpha;
    /** The viscosity of the gas. */
    const double mu;
    /** The ambient pressure. */
    const double Pamb;
    /** The initial height of the drop. */
    const double h0;
    /** The drop radius. */
    const double R;
    /** The initial velocity of the drop. */
    const double V;
    /** The surface tension at the liquid-gas interface. */
    const double sigma;
    gl_params(double alpha_,double mu_,double Pamb_,double h0_,double R_,double V_,double sigma_) :
        alpha(alpha_), mu(mu_), Pamb(Pamb_), h0(h0_), R(R_), V(V_),
        sigma(sigma_) {}
    gl_params(gl_params &g) : alpha(g.alpha), mu(g.mu), Pamb(g.Pamb), h0(g.h0),
        R(g.R), V(g.V), sigma(g.sigma) {}
};

/** \brief A class for solving for the evolving height of the gas layer via a
 * partial differential equation, and using the result to compute pressure and
 * tangential stress. */
class gas_layer_solve: public gas_layer, public gl_params {
    public:
        /** A constant to distinguish between the (0,bx) and (-bx,bx) cases. */
        const bool x_sym;
        /** The reciprocal of the grid spacing. */
        const double xsp;
        /** The squared reciprocal of the grid spacing. */
        const double xxsp;
        /** An array for tracking the flow through the bottom boundary. */
        double* const fbd;
        /** An array of the droplet height. */
        double* const h;
        /** An array of the x derivative of the droplet height. */
        double* const hx;
        /** An array of the t derivative of the droplet height. */
        double* const ht;
        /** The previous pressure. */
        double* const pg_old;
        /** The current pressure. */
        double* const pg;
        void step_forward(double dt);
        gas_layer_solve(int m_,double ax_,double bx_,gl_params &glp,bool x_sym_);
        virtual ~gas_layer_solve();
        void update_gas_pressure(double dt);
        void update_height(double dt);
        virtual void init_fields();
        virtual void set_field_pointer(field* fm_) {fm=fm_;fb=fm-ml;}
        virtual void output_profile(const char *filename,float *buf,int sn,double ay,bool height);
        virtual void output_profile_double(const char *filename,float *buf,int sn,bool height);
    protected:
        virtual void calc_p_and_tau(double t,double dt);
        virtual void load_restart_internal(double time,FILE *inf);
        virtual void save_restart_internal(FILE *outf);
    private:
        /** A pointer to the (0,0) gridpoint in the fluid_2d field data
         * structure from which to read fluid velocities. */
        field* fm;
        /** A pointer to the (0,-1) gridpoint in the fluid_2d field data
         * structure from which to read fluid velocities. */
        field* fb;
        /** The tridiagonal solver. */
        tri_solve t;
        void calculate_hx();
        void calculate_liquid_fields();
        void compute_frame_acceleration(double dt);
        void implicit_step(double dt);
        double newton_step(double dtinv);
        inline double h_ref(int i);
        inline double eno2(double p0,double p1,double p2,double p3);
};

#endif
