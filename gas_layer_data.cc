#include "common.hh"
#include "gas_layer_data.hh"

#include <cstdio>
#include <cstring>

/** Initializes the class for representing the pressures and shear stresses due the gas layer,
 * for the case when it is set from a data file.
 * \param[in] m_ the number of grid cells.
 * \param[in] (ax_,bx_) the horizontal coordinate range to compute over.
 * \param[in] (m_dat,n_dat) the dimensions of the data files.
 * \param[in] t_end the end time of the information in the data files.
 * \param[in] path_pressure the path to the pressure data file.
 * \param[in] path_tau the path to the shear stress data file. */
gas_layer_data::gas_layer_data(int m_,double ax_,double bx_,int m_dat,int n_dat,
        double t_end,const char* path_pressure,const char* path_tau) : gas_layer(m_,ax_,bx_,0),
    p_array(m_dat,n_dat,ax_,bx_,0,t_end), tau_array(m_dat,n_dat,ax_,bx_,0,t_end) {

    // Read in the two data files
    read_array(m_dat,n_dat,p_array.u,path_pressure);
    read_array(m_dat,n_dat,tau_array.u,path_tau);
}

/** Updates the pressure and shear stress arrays.
 * \param[in] t the current simulation time.
 * \param[in] dt the size of the timestep. */
void gas_layer_data::calc_p_and_tau(double t,double dt) {

    // Update the pressure
    for(int i=0;i<me;i++) p[i]=p_array.f(ax+i*dx,t+dt);

    // Update the shear stress
    const double axh=ax+0.5*dx;
    for(int i=0;i<m;i++) tau[i+2]=tau_screened(axh+i*dx,t+dt);
}

/** Reads a text file into a bicubic interpolation array.
 * \param[in] u a pointer the array in the bicubic interpolation class.
 * \param[in] path the path of the file to read from. */
void gas_layer_data::read_array(int m_dat,int n_dat,double *u,const char* path) {
    FILE *fp=safe_fopen(path,"r");
    for(int i=0;i<m_dat;i++) for(int j=0;j<n_dat;j++)
        if(fscanf(fp,"%lg",u+i+m_dat*j)!=1)
            fatal_error("Data file import error",1);
    fclose(fp);
}

/** Calculates the shear stress using the bicubic interpolation of the provided
 * data, screening out some values to zero if they are outside the x range of
 * [25%,75%] of the simulation domain.
 * \param[in] (x,t) the position and time at which to interpolate. */
double gas_layer_data::tau_screened(double x,double t) {
    return x>=ax+0.75*(bx-ax)||x<=ax+0.25*(bx-ax)?0.:tau_array.f(x,t);
}
