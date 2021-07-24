#ifndef VDI_BI_INTERP_HH
#define VDI_BI_INTERP_HH

/** \brief A class to perform bicubic interpolation of a field. */
class bicubic_interp {
    public:
        /** The number of grid points in the x direction. */
        const int m;
        /** The number of grid points in the y direction. */
        const int n;
        /** The lower x coordinate. */
        const double ax;
        /** The lower y coordinate. */
        const double ay;
        /** The reciprocal of the x grid spacing. */
        const double xsp;
        /** The reciprocal of the y grid spacing. */
        const double ysp;
        /** The grid of values to interpolate. */
        double* const u;
        bicubic_interp(int m_,int n_,double ax_,double bx_,double ay_,double by_);
        ~bicubic_interp();
        double f(double x,double y);
        double f_grad_f(double x,double y,double &fx,double &fy);
    private:
        /** Temporary workspace for computing the interpolation. */
        double a[16];
        /** The index of grid point currently referenced in the temporary workspace. */
        int ijc;
        void table_setup(int i,int j,int ij);
        void compute_x(int i,double *up,double &c0,double &c1,double &c2,double &c3);
        void grid_index(double &x,double &y);
        void fill_ad(double *ap,double c1,double c2,double c3);
        void fill_au(double *ap,double c0,double c1,double c2);
        void fill_a(double *ap,double c0,double c1,double c2,double c3);
        /** Calculates a cubic interpolant value.
         * \param[in] ap a pointer to the cubic coefficients.
         * \param[in] y the function argument. */
        inline double yl(double *ap,double y) {
            return *ap+y*(ap[1]+y*(ap[2]+y*ap[3]));
        }
        /** Calculates the first derivative of the cubic interpolant.
         * \param[in] ap a pointer to the cubic coefficients.
         * \param[in] y the function argument.*/
        inline double dyl(double *ap,double y) {
            return ap[1]+y*(2*ap[2]+3*y*ap[3]);
        }
};

#endif
