#ifndef VDI_GAS_LAYER_HH
#define VDI_GAS_LAYER_HH

#include <cstdio>

#include "common.hh"
#include "fields.hh"

/** \brief A class for updating the pressure and shear stress in the gas layer.
 */
class gas_layer {
    public:
        /** The number of grid cells. */
        const int m;
        /** The number of grid cell corners. */
        const int me;
        /** The number of grid cells, including ghost cells. */
        const int ml;
        /** The lower coordinate bound. */
        const double ax;
        /** The upper coordinate bound. */
        const double bx;
        /** The grid spacing. */
        const double dx;
        /** The current downward frame velocity. */
        double Vframe;
        /** The pressure at the current step. */
        double* const p;
        /** The shear stress at the current step. */
        double* const tau;
        /** The pressure at the previous step. */
        double* const p_old;
        /** The shear stress at the previous step. */
        double* const tau_old;
        gas_layer(int m_,double ax_,double bx_,double Vframe_);
        virtual ~gas_layer();
        void update(double t,double dt);
        /** Loads the gas layer information from a restart file.
         * \param[in] time the current simulation time.
         * \param[in] inf the file handle to read from. */
        inline void load_restart(double time,FILE *inf) {
            safe_fread(&Vframe,sizeof(double),1,inf,"Vframe");
            load_restart_internal(time,inf);
        }
        /** Saves the gas layer information to a restart file.
         * \param[in] outf the file handle to save to. */
        inline void save_restart(FILE *outf) {
            fwrite(&Vframe,sizeof(double),1,outf);
            save_restart_internal(outf);
        }
        virtual void init_fields() = 0;
        virtual void set_field_pointer(field *fm_) = 0;
        virtual void output_profile(const char *filename,float *buf,int sn,double ay,bool height) {}
        virtual void output_profile_double(const char *filename,float *buf,int sn,bool height) {}
    protected:
        virtual void calc_p_and_tau(double t,double dt) = 0;
        virtual void load_restart_internal(double time,FILE *inf) = 0;
        virtual void save_restart_internal(FILE *outf) = 0;
        void apply_tau_bc();
};

#endif
