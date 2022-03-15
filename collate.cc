#include <cstdio>
#include <cstdlib>
#include <cfloat>
#include <cstring>
#include <sys/stat.h>

#include "common.hh"
#include "fileinfo.hh"

int main(int argc,char **argv) {

    // Check for the correct number of command-line arguments
    if(argc!=3) {
        fputs("Syntax: ./collate [hei|pgl] <output_dir>\n\n"
              "The output directory should have a '.odr' suffix.\n",stderr);
        return 1;
    }

    // Check for the correct command-line argument
    int fld=0;
    if(strcmp(argv[1],"pgl")==0) fld=1;
    else if(strcmp(argv[1],"hei")!=0)
        fatal_error("First command-line argument should be 'hei' or 'pgl'",1);

    // Check that the filename ends in '.odr'
    const char* fn=argv[argc-1];
    int l=strlen(fn),i;
    if(l<4) fatal_error("Filename is too short",1);
    if(l>224) fatal_error("Filename is too long",1);
    const char* ip=fn+l-4;
    if(*ip!='.'||ip[1]!='o'||ip[2]!='d'||ip[3]!='r')
        fatal_error("Filename must end in '.odr'",1);

    // Read in the configuration file to get the frame information
    char buf[256],*fp=buf+l-3;
    memcpy(buf,fn,l-3);
    *fp='c';fp[1]='f';fp[2]='g';fp[3]=0;
    fileinfo fi(buf);
    double fdt=fi.t_end/fi.nframes;

    // Assemble the output filename by replacing '.odr' with the suffix
    *fp=argv[1][0];fp[1]=argv[1][1];fp[2]=argv[1][2];
    FILE *outf=safe_fopen(buf,"wb");

    // Set up output matrix
    int q=fld==0?fi.m:fi.m+1;
    float *o=new float[(l+1)*(fi.nframes+2)],*op=o;
    double *uu=new double[l+1],offset=fld==0?0.5:0;
    *(op++)=q;
    for(int i=0;i<fi.m;i++) *(op++)=fi.ax+(i+offset)*(fi.bx-fi.ax)/fi.m;
   
    // Loop over the height files
    struct stat st;
    FILE *inf;
    for(i=0;i<=fi.nframes;i++) {

        // Quit if the pressure file doens't exist
        sprintf(buf,"%s/%sfull.%d",fn,fld==0?"h":"pg",i);
        if(stat(buf,&st)!=0) break;

        // Read in the header
        inf=safe_fopen(buf,"rb");
        safe_fread(&l,sizeof(int),1,inf,"grid size");
        if(l!=q) fatal_error("Grids change size",1);
        safe_fread(uu,sizeof(double),2,inf,"grid dimensions");

        // Read in the data
        safe_fread(uu,sizeof(double),q,inf,"height data");

        // Convert data into single precision
        *(op++)=i*fdt;
        for(int j=0;j<q;j++) *(op++)=uu[j];
        fclose(inf);
    }

    fwrite(o,sizeof(float),(q+1)*(i+1),outf);
    fclose(outf);
    printf("# %s: %d files processed, fdt = %g µs\n",fn,i,fdt*1e6);

}
