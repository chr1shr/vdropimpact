#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/stat.h>

#include "common.hh"
#include "gp_matrix.hh"

// Maximum number of files to check for
const int files_max=16384;

int main(int argc,char **argv) {

    // Check for the correct number of command-line arguments
    if(argc!=2) {
        fputs("Syntax: ./p_strip <output_dir>\n\n"
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

    // Loop over the height files
    struct stat st;
    char buf[256];
    for(i=0;i<=files_max;i++) {

        // Quit if the height file doens't exist
        sprintf(buf,"%s/p.%d",argv[1],i);
        if(stat(buf,&st)!=0) break;

        // Load in the full pressure file
        gp_matrix gp(buf);

        // Strip out the first line of the file
        sprintf(buf,"%s/pbase.%d",argv[1],i);
        FILE *fp=safe_fopen(buf,"wb");
        float g=gp.m;
        fwrite(&g,sizeof(float),1,fp);
        fwrite(gp.x,sizeof(float),gp.m,fp);
        g=0;
        fwrite(&g,sizeof(float),1,fp);
        fwrite(gp.f,sizeof(float),gp.m,fp);
        fclose(fp);
    }
    printf("%s: %d files processed\n",argv[1],i);
}
