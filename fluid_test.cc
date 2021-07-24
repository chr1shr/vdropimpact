#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>

#include "fluid_2d.hh"

int main(int argc, char **argv) {

    // Check for the correct number of command-line arguments
    if(argc!=2&&(argc!=3||strcmp(argv[1],"-r")!=0)) {
        fputs("Syntax: ./fluid_test [-r] <input_file>\n\n"
              "The input file should have a '.cfg' suffix\n",stderr);
        return 1;
    }

    // Read in the parameters from the input file
    fluid_2d f2d(argv[argc-1]);

    // Create output directory
    mkdir(f2d.filename,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

    // Initialize the fluid tracers and the timestep
    argc==2?f2d.init_fields():f2d.load_last_restart();

    // Run the simulation
    f2d.solve();
}
