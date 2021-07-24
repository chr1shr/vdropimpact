#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/stat.h>

#include "common.hh"
#include "gp_matrix.hh"

// Number of grid points for doing the quadratic fit
const int m=8;

// Maximum number of files to check for
const int files_max=16384;

int main(int argc,char **argv) {

    // Check for the correct number of command-line arguments
    if(argc!=2) {
        fputs("Syntax: ./h_analysis <output_dir>\n\n"
              "The output directory should have a '.odr' suffix\n",stderr);
        return 1;
    }

    // Check that the filename ends in '.odr'
    int l=strlen(argv[1]),i;
    if(l<4) fatal_error("Filename is too short",1);
    if(l>224) fatal_error("Filename is too long",1);
    const char* ip=argv[1]+l-4;
    if(*ip!='.'||ip[1]!='o'||ip[2]!='d'||ip[3]!='r')
        fatal_error("Filename must end in '.odr'",1);

    // Assemble output filename by replacing '.odr' with '.han'
    char buf[256],*fp=buf+l-3;
    memcpy(buf,argv[1],l-3);
    *fp='h';fp[1]='a';fp[2]='n';fp[3]=0;

    // Loop over the height files
    struct stat st;
    FILE *of=safe_fopen(buf,"w");
    for(i=0;i<=files_max;i++) {

        // Quit if the height file doens't exist
        sprintf(buf,"%s/height.%d",argv[1],i);
        if(stat(buf,&st)!=0) break;

        // Loop over the first few field values to do quadratic interpolation
        gp_matrix gp(buf);
        double xx,sxx=0,sx4=0,sf=0,sfxx=0;
        for(int k=0;k<m;k++) {
            xx=gp.x[k]*gp.x[k];
            sf+=gp.f[k];
            sfxx+=gp.f[k]*xx;
            sxx+=xx;
            sx4+=xx*xx;
        }

        // Print coefficients of quadratic fit
        double idet=1./(sx4*m-sxx*sxx);
        fprintf(of,"%d %g %g\n",i,idet*(sx4*sf-sxx*sfxx),idet*(-sxx*sf+m*sfxx));
    }
    printf("%s: %d files processed\n",argv[1],i);
    fclose(of);
}
