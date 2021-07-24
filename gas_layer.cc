#include "gas_layer.hh"

#include <cstring>

/** Initializes the base gas layer class, which contains grid constants and memory
 * for pressure and shear stress that are common to both gas layer models.
 * \param[in] m_ the number of grid cells.
 * \param[in] (ax_,bx_) the horizontal coordinate range to compute over.
 * \param[in] Vframe_ the initial downward frame velocity. */
gas_layer::gas_layer(int m_,double ax_,double bx_,double Vframe_):
    m(m_), me(m_+1), ml(m_+4), ax(ax_), bx(bx_), dx((bx-ax)/m),
    Vframe(Vframe_), p(new double[me]), tau(new double[ml]),
    p_old(new double[me]), tau_old(new double[ml]) {}

/** The class destructor frees the dynamically allocated memory. */
gas_layer::~gas_layer() {
    delete [] tau_old;
    delete [] p_old;
    delete [] tau;
    delete [] p;
}

/** Applies the boundary conditions on tau, by imposing tau=0 at the left and
 * right boundaries. This is the same for gas_layer_solve and gas_layer_data. */
void gas_layer::apply_tau_bc() {
    *tau=-tau[3];
    tau[1]=-tau[2];
    tau[ml-2]=-tau[ml-3];
    tau[ml-1]=-tau[ml-4];
}

/** Performs updates of the gas layer that are common between the "data" and
 * "solve" options.
 * \param[in] time the current simulation time.
 * \param[in] dt the current simulation timestep. */
void gas_layer::update(double time,double dt) {

    // Store the previous pressure and shear stress
    memcpy(p_old,p,me*sizeof(double));
    memcpy(tau_old,tau,ml*sizeof(double));

    // Compute the new pressure and shear stress, and apply the boundary
    // conditions to the shear stress that are common for both models
    calc_p_and_tau(time,dt);
    apply_tau_bc();
}
