#ifndef VDI_FILEINFO_HH
#define VDI_FILEINFO_HH

#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>

#include "gas_layer.hh"

/** The size of the temporary buffer for parsing the input file. */
const int fileinfo_buf_size=512;

/** \brief A class for parsing the drop impact parameters from a text
 * configuration file. */
struct fileinfo {
    public:
        /** The number of grid cells in the horizontal direction. */
        int m;
        /** The number of grid cells in the vertical direction. */
        int n;
        /** The types of file output to save. */
        unsigned int fflags;
        /** The number of tracers. */
        int ntrace;
        /** The number of frames. */
        int nframes;
        /** Restart file output frequency. */
        int restart_freq;
        /** The non-inertial frame mode. (0: inactive, 1: central velocity, 2:
         * averaging window.) */
        int nif_mode;
        /** A boolean variable to determine if the viscous term will be
         * calculated implicitly. */
        bool implicit_visc;
        /** A boolean variable to determine whether u, v is
         * antisymmetric, symmetric about the left boundary. */
        bool x_sym;
        /** A boolean variable to enable machine-readable timing information.
         */
        bool mr_time_output;
        /** The lower bound in the x direction. */
        double ax;
        /** The lower bound in the y direction. */
        const double ay;
        /** The upper bound in the x direction. */
        double bx;
        /** The upper bound in the y direction. */
        double by;
        /** The grid spacing in the x direction. */
        double dx;
        /** The grid spacing in the y direction. */
        double dy;
        /** The inverse grid spacing in the x direction. */
        double xsp;
        /** The inverse grid spacing in the y direction. */
        double ysp;
        /** The dynamic viscosity of the fluid. */
        double mu;
        /** The reciprocal of the dynamic viscosity of the fluid. */
        double muinv;
        /** The density of the fluid. */
        double rho;
        /** The reciprocal of the density of the fluid. */
        double rhoinv;
        /** The padding factor for the timestep for the physical terms, which
         * should be smaller than 1. */
        double tmult;
        /** The simulation duration. */
        double t_end;
        /** The output directory filename. */
        char *filename;
        /** A pointer to the gas layer model. */
        gas_layer* gl;
        fileinfo(const char* infile);
        /** The class destructor frees the dynamically allocated memory. */
        ~fileinfo() {
            delete gl;
            delete [] filename;
        }
    protected:
        /** The lower index for non-inertial frame mode 2. */
        int anif;
        /** The upper index for non-inertial frame mode 2. */
        int bnif;
    private:
        void set_nif_range(double w);
        /** Finds the next token in a string and interprets it as a double
         * precision floating point number. If none is availble, it gives an
         * error message.
         * \param[in] ln the current line number. */
        inline double next_double(int ln) {
            return atof(next_token(ln));
        }
        /** Finds the next token in a string, interprets it as a double
         * precision floating point number, and checks that there are no
         * subsequent values.
         * \param[in] ln the current line number. */
        inline double final_double(int ln) {
            double temp=next_double(ln);
            check_no_more(ln);
            return temp;
        }
        /** Finds the next token in a string, interprets it as an integer, and
         * checks that there are no subsequent values.
         * \param[in] ln the current line number. */
        inline int final_int(int ln) {
            int temp=atoi(next_token(ln));
            check_no_more(ln);
            return temp;
        }
        /** Tests to see if two strings are equal.
         * \param[in] p1 a pointer to the first string.
         * \param[in] p2 a pointer to the second string.
         * \return True if they are equal, false otherwise. */
        inline bool se(const char *p1,const char *p2) {
            return strcmp(p1,p2)==0;
        }
        char* next_token(int ln);
        void check_no_more(int ln);
        void check_invalid(double val,const char *p);
};

#endif
