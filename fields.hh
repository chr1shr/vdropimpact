#ifndef VDI_FIELDS_HH
#define VDI_FIELDS_HH

#include <cmath>

/** \brief Data structure for storing the fields at grid points. */
struct field {
    /** The horizontal velocity. */
    double u;
    /** The vertical velocity. */
    double v;
    /** The pressure. */
    double p;
    /** Gradient of velocity, stored at the cell-center. */
    double ux,uy,vx,vy;
    /** Velocities extrapolated to edges. */
    double ud,vd,ul,vl,ur,vr,uu,vu;
    /** Temporary storage for intermediate quantities. */
    double c0,c1,c2,c3;
    /** Sets the velocity components and pressure to zero. */
    inline void clear_main() {
        u=v=p=0;
    }
    /** Packs the primary fields into an array.
     * \param[in,out] pp a pointer to the array. */
    inline void pack(double *&pp) {
        *(pp++)=u;*(pp++)=v;*(pp++)=p;
    }
    /** Unpacks the primary fields from an array.
     * \param[in,out] pp a pointer to the array. */
    inline void unpack(double *&pp) {
        u=*(pp++);v=*(pp++);p=*(pp++);
    }
};

#endif
