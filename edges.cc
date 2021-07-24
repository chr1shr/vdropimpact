#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>

#include "common.hh"
#include "gp_matrix.hh"

/** Calculates a Lagrange extrapolation from field values at positions
 * -k,..,0 to the point 0.5.
 * \param[in] f a pointer to the zero position.
 * \param[in] d a memory displacement factor between the zero and one
 *              positions.
 * \param[in] k the order of the extrapolation (between 1 and 3).
 * \return The extrapolated value. */
double quad(float *f,int d,int k) {
    return k==1?*f*2-f[d]:(k==2?*f*1.875-1.25*f[d]+0.375*f[d<<1]
                   :(*f-f[d])*2.1875+1.3125*f[d<<1]-0.3125*f[3*d]);
}

int main(int argc,char **argv) {

    // Check for the correct number of command-line arguments
    if(argc!=3) {
        fputs("Usage: ./edges <output_dir> <suffix>\n",stderr);
        return 1;
    }

    // Read the frame number and open the velocity components
    char buf[256];
    int fr=atoi(argv[2]);
    sprintf(buf,"%s/u.%d",argv[1],fr);
    gp_matrix gpu(buf);
    sprintf(buf,"%s/v.%d",argv[1],fr);
    gp_matrix gpv(buf);

    // Check the the dimensions of the two velocity components agree
    int m=gpu.m,n=gpu.n;
    if(gpu.m!=m||gpu.m!=m) {
        fputs("Size mismatch between velocity components\n",stderr);
        return 1;
    }

    // Open the file to save the edge velocity information
    sprintf(buf,"%s/ed.%d",argv[1],fr);
    FILE *fp=safe_fopen(buf,"w");

    // Extrapolate the velocity to the edges, using three different
    // extrapolation orders
    float *fu=gpu.f,*fv=gpv.f;
    for(int i=0,mn=gpu.mn,o=m>n?m:n;i<o;i++) {
        if(i<m) fprintf(fp,"%g",gpu.x[i]);else fputs("NaN",fp);
        if(i<n) fprintf(fp," %g",gpu.y[i]);else fputs(" NaN",fp);
        for(int k=1;k<=3;k++) {
            if(i<m) fprintf(fp," %g %g %g %g",quad(fu+i,m,k),quad(fv+i,m,k),
                    quad(fu+(mn-m+i),-m,k),quad(fv+(mn-m+i),-m,k));
            else fputs(" NaN NaN NaN NaN",fp);
            if(i<m) fprintf(fp," %g %g %g %g",quad(fu+m*i,1,k),quad(fv+m*i,1,k),
                    quad(fu+(m*i+(m-1)),-1,k),quad(fv+(m*i+(m-1)),-1,k));
            else fputs(" NaN NaN NaN NaN",fp);
        }
        fputs("\n",fp);
    }
    fclose(fp);
}
