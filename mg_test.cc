#include <cstdio>
#include <cstdlib>

#include "tgmg.hh"
#include "mgs_mac.hh"
#include "mgs_fem.hh"

// Grid dimensions
const double ax=-1,bx=1,ay=-1,by=1;

// Safe grid limit
const int grid_limit=32768;

// Minimum time in seconds to perform V-cycles for, in order to get an accurate
// timing result
const double v_cycle_time_min=10;

// Template for setting up the multigrid hierarchy and solving the problem
template<class T>
void setup_and_do(T& msu,bool solve,double dx,double dy) {
    msu.mg.setup();

    if(solve) {
        msu.solve_v_cycle();

        // Output the fields using the Gnuplot binary matrix format
        msu.mg.output_b("b.0",ax,dx,ay,dy);
        msu.mg.output_z("z.0",ax,dx,ay,dy);
        msu.mg.output_res("r.0",ax,dx,ay,dy);
    } else {
        double t0=wtime(),t1;
        int k=0;

        // To get accurate timing statistics, perform V-cycles until a certain
        // time has elasped
        do {
            ++k;
            msu.mg.v_cycle();
        } while((t1=wtime())<t0+v_cycle_time_min);

        // Calculate and print timing statistics
        t1-=t0;
        printf("%d V-cycles performed in %.2f s\n%g ms per V-cycle\n",k,t1,1e3*t1/k);
    }
}

int main(int argc,char **argv) {

    // Check for the right number of command-line arguments. If it's incorrect
    // then print a syntax message.
    if(argc<3||argc>4) {
        fputs("Syntax: ./mg_test <case> <x gridpoints> [<y gridpoints>]\n\n"
              "Case=0 : MAC correction (solve)\n"
              "Case=1 : MAC correction (time V-cycle)\n"
              "Case=2 : Finite-element projection (solve)\n"
              "Case=3 : Finite-element projection (time V-cycle)\n\n"
              "If y gridpoints aren't specified, then the code uses a square grid.\n",stderr);
        return 1;
    }

    // Check case number
    int ca=atoi(argv[1]);
    if(ca<0||ca>3) {
        fputs("Case number should be between 0 and 3\n",stderr);
        return 1;
    }

    // Check horizontal grid size
    int m=atoi(argv[2]);
    if(m<=8||m>grid_limit) {
        fputs("x gridpoints out of range\n",stderr);
        return 1;
    }

    // Check vertical grid size
    int n;
    if(argc==3) n=m;
    else {
        n=atoi(argv[3]);
        if(n<=8||n>grid_limit) {
            fputs("y gridpoints out of range\n",stderr);
            return 1;
        }
    }

    // Allocate memory and set up grids
    int ij,mn=m*n,i1=m/5,j1=n/3,i2=(6*m)/7,j2=(2*n)/3;
    double dx=(bx-ax)/(m-1),dy=(by-ay)/(n-1),*b=new double[mn],*z=new double[mn];
    for(ij=0;ij<mn;ij++) b[ij]=0;
    for(ij=0;ij<mn;ij++) z[ij]=0;
    b[i1+j1*m]=1;b[i2+j2*m]=-1;

    // Solve the multigrid problem
    if(ca<2) {
        mgs_mac msu(m,n,dx,dy,b);
        setup_and_do(msu,ca==0,dx,dy);
    } else {
        mgs_fem msu(m,n,dx,dy,b);
        setup_and_do(msu,ca==2,dx,dy);
    }

    // Delete dynamically allocated memory
    delete [] z;
    delete [] b;
}
