#ifndef VDI_GAS_LAYER_DATA_HH
#define VDI_GAS_LAYER_DATA_HH

#include "bi_interp.hh"
#include "gas_layer.hh"

/** \brief A class for setting pressure and tangential stress in the gas layer
 * from a data file. */
class gas_layer_data : public gas_layer {
    public:
        gas_layer_data(int m_,double ax_,double bx_,int m_dat,int n_dat,double t_end,const char* path_pressure,const char* path_tau);
        virtual ~gas_layer_data() {}
        virtual void init_fields() {update(0,0);}
        virtual void set_field_pointer(field *fm_) {}
    protected:
        virtual void calc_p_and_tau(double t,double dt);
        virtual void load_restart_internal(double time,FILE *fp) {update(time,0);}
        virtual void save_restart_internal(FILE *fp) {}
    private:
        void read_array(int m_dat,int n_dat,double *u,const char* path);
        double tau_screened(double x,double t);
        /** A bicubic interpolation of the pressure at the bottom boundary. */
        bicubic_interp p_array;
        /** A bicubic interpolation of the tangential stress at the bottom
         * boundary. */
        bicubic_interp tau_array;
};

#endif
