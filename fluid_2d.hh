#ifndef VDI_FLUID_2D_HH
#define VDI_FLUID_2D_HH

#include <cstdio>

#include "common.hh"
#include "fields.hh"
#include "fileinfo.hh"
#include "tgmg.hh"
#include "mgs_mac.hh"
#include "mgs_fem.hh"
#include "visco_impl.hh"

#ifdef _OPENMP
#include "omp.h"
#else
#include <ctime>
#endif

/** \brief A class to compute the flow field in the liquid. */
class fluid_2d : public fileinfo {
    public:
        /** The total number of grid cells. */
        const int mn;
        /** The number of corner grid points in the horizontal direction. */
        const int me;
        /** The number of corner grid points in the vertical direction. */
        const int ne;
        /** The memory step length, taking into account ghost point allocation. */
        const int ml;
        /** The grid spacing in the x direction. */
        const double dx;
        /** The grid spacing in the y direction. */
        const double dy;
        /** The inverse grid spacing in the x direction. */
        const double xsp;
        /** The inverse grid spacing in the y direction. */
        const double ysp;
        /** The square inverse grid spacing in the x direction. */
        const double xxsp;
        /** The square inverse grid spacing in the y direction. */
        const double yysp;
        /** The regular timestep choice. */
        const double dt_reg;
        /** An array containing the simulation fields. */
        field* const fbase;
        /** A pointer to the (0,0) grid cell in the field array. */
        field* const fm;
        /** An array containing the tracer positions. */
        double* const tm;
        /** An array for the source terms used during the algebraic
         * multigrid solve.*/
        double* const src;
        /** A pointer to within the source term array to store extra BC
         * information needed in the MAC solve. */
        double* const s_bc;
        /** The current simulation time. */
        double time;
        /** The current frame number. */
        int f_num;
        fluid_2d(const char* infile);
        ~fluid_2d();
        void solve();
        void step_forward(double dt);
        void init_fields();
        void load_restart(int sn);
        void load_last_restart();
        void save_restart(int sn);
        void write_files(int k);
        void init_tracers();
        void update_tracers(double dt);
        void output(const char *prefix,const int mode,int sn,bool ghost=false);
        void output_tracers(const int sn);
        void set_nif_mode(int mode,double w=0);
        /** Chooses a timestep size that is the largest value smaller than dt_reg,
         * such that a given interval length is a perfect multiple of this
         * timestep.
         * \param[in] interval the interval length to consider.
         * \param[out] adt the timestep size.
         * \return The number of timesteps the fit into the interval. */
        inline int timestep_select(double interval,double &adt) {
            int l=static_cast<int>(interval/dt_reg)+1;
            adt=interval/l;
            return l;
        }
        /** Opens a file in the output directory with a specific filename.
         * \param[in] prefix the field name to use as the filename prefix.
         * \param[in] sn the current frame number to append to the filename.
         * \param[in] mode the cstdio fopen mode to use. */
        inline FILE* odir_open(const char* prefix,int sn,const char* mode) {
            char *bufc=reinterpret_cast<char*>(buf);
            sprintf(bufc,"%s/%s.%d",filename,prefix,sn);
            return safe_fopen(bufc,mode);
        }
    private:
        /** The multigrid solver for the marker-and-cell (MAC) correction
         * problem. */
        mgs_mac ms_mac;
        /** The multigrid solver for the finite-element elliptic problem used
         * in the pressure projection. */
        mgs_fem ms_fem;
        /** The multigrid solver for the implicit solution of the horizontal
         * velocity component. */
        visco_impl* vix;
        /** The multigrid solver for the implicit solution of the vertical
         * velocity component. */
        visco_impl* viy;
        /** Temporary storage for used during the output routine. */
        float *buf;
        double t_gl,t_bc;
        void set_boundaries(double dt);
        void smooth_line(field *fp,int d);
        void adjust_frame();
        void set_inviscid_boundaries(double dt);
        void set_viscous_boundaries();
        void simpsons_velocity_bc(double x,double y,double prefac,double &u,double &v);
        void simpsons_terms(double xms,double xps,double ysq,double p,double &fsu,double &fsv);
        void zero_edge_velocities();
        double mono_diff(double &f0,double &f1,double &f2,double &f3,double &f4);
        inline double delta_lim(double &f0,double &f1,double &f2);
        inline double delta_f(double &f0,double &f1,double &f2);
        inline double min(double a,double b) {return a<b?a:b;}
        inline double max(double a,double b) {return a>b?a:b;}
        void centered_diff(field *fp,double &uxx,double &uyy,double &vxx,double &vyy);
        void godunov_choose();
        inline void godunov_set(double &u0,double &u1,double &v0,double &v1);
        inline void godunov_set_tang(double &u0,double &u1,double &v0,double &v1);
        inline bool restart_exists(int sn);
        /** Computes the vorticity at the lower left corner of a grid cell, based
         * on centered finite differences of the velocity.
         * \param[in] fp a pointer to the grid cell. */
        inline double vorticity(field *fp) {
            return 0.5*(xsp*(fp->v-fp[-1].v+fp[-ml].v-fp[-ml-1].v)
                       -ysp*(fp->u+fp[-1].u-fp[-ml].u-fp[-ml-1].u));
        }
#ifdef _OPENMP
        inline double wtime() {return omp_get_wtime();}
#else
        inline double wtime() {return double(clock())*(1./CLOCKS_PER_SEC);}
#endif
};

#endif
