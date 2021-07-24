#include <cstdio>
#include <cstdlib>
#include <cfloat>
#include <cstring>
#include <sys/stat.h>

#include "common.hh"
#include "fileinfo.hh"
#include "gas_layer_solve.hh"
#include "gp_matrix.hh"

// Maximum number of files to check for
const int files_max=16384;

// The height of the evanescent field (m)
float h_ev=500e-9;

// The length scale of the smoothing Gaussian (m)
float glen=1.2e-6;

/** The number of standard deviations to consider when doing Gaussian
 * smoothing. */
float gaussian_fac=7.;

/** Function to calculate the cubic Hermite interpolant over the interval
 * [0,1].
 * \param[in] x the position at which to evaluate the interpolant.
 * \param[in] (t0,t1,t2,t3) the coefficients of the Hermite basis functions. */
double hermite(double x,double t0,double t1,double t2,double t3) {
    return t0*x*x*(3-2*x)+x*(1-x)*(-t1*(1-x)+t2*x)+t3*(2*x*x*x-3*x*x+1);
}

/** Applies Gaussian smoothing to a Gnuplot field
 * \param[in] gp a reference to the Gnuplot field to smooth.
 * \param[in] q an array in which to write the smoothed values. */
void smooth(gp_matrix &gp,double *q) {
    float sls=glen/gp.dx,fac=0.5/(sls*sls),v=0.5;
    int &m=gp.m,cut=int(sls*gaussian_fac)+1;

    // Compute the normalizing factor for the weights
    double *sp=new double[2*cut+1],*sp2=sp+cut;
    for(int i=1;i<=cut;i++) v+=sp2[i]=exp(-fac*i*i);
    v=0.5/v;

    // Normalize and symmetrize
    *sp2=v;
    for(int i=1;i<=cut;i++) sp2[-i]=sp2[i]*=v;

    // Loop over the grid points and calculated the smoothed value
    for(int i=0;i<gp.m;i++) {
        v=0;
        for(int j=i-cut;j<=i+cut;j++)
            v+=sp2[j-i]*gp.f[j<0?-1-j:(j>=m?2*m-1-j:j)];
        q[i]=v;
    }
    delete [] sp;
}

int main(int argc,char **argv) {

    // Check for the correct number of command-line arguments
    if(argc<2||(argc>3&&strcmp(argv[1],"-r")!=0)) {
        fputs("Syntax: ./tt_analysis [-r] <output_dir>\n\n"
              "The output directory should have a '.odr' suffix. The '-r' flag uses the\n"
              "raw values without smoothing and cubic interpolation.\n",stderr);
        return 1;
    }

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
    gas_layer_solve *gls=reinterpret_cast<gas_layer_solve*>(fi.gl);
    double fdt=fi.t_end/fi.nframes,tmin0=gls->h0/gls->V,tmin0_ev=0;

    // Assemble the output filename by replacing '.odr' with '.tta'
    *fp='t';fp[1]='t';fp[2]=argc==3?'b':'a';fp[3]=0;

    // Open output file
    struct stat st;
    FILE *of=safe_fopen(buf,"w");

    // Loop over the height files
    int gimin=-1,m;
    double gfmin=DBL_MAX,*q=NULL;
    bool efield=true,hsearch=true;
    for(i=0;i<=files_max;i++) {

        // Quit if the height file doens't exist
        sprintf(buf,"%s/height.%d",fn,i);
        if(stat(buf,&st)!=0) break;

        // Loop over the first few field values to do quadratic interpolation
        gp_matrix gp(buf);

        // Create temporary array
        if(i==0) q=new double[m=gp.m];
        else if(m!=gp.m) fatal_error("Grids change size",1);

        // Apply smoothing if needed
        if(argc==2) smooth(gp,q);
        else for(int j=0;j<gp.m;j++) q[j]=gp.f[j];

        // Loop over the grid points to find the minimum
        double fmin=gp.f[0],xmin=gp.x[0];
        for(int k=0;k<gp.m-1;k++) {
            double *f=q+k,
                   t0=f[1],
                   t1=k==gp.m-1?(f[1]-*f):0.5*(f[2]-*f),
                   t2=k==0?(f[1]-*f):0.5*(f[1]-f[-1]),
                   t3=*f,
                   a=6*(t3-t0)+3*(t1+t2),
                   b=6*t0-2*t1-4*t2-6*t3,
                   c=t2,
                   det=b*b-4*a*c;

            // Check if the value at the grid point is a new minimum
            if(f[1]<fmin) fmin=f[1],xmin=gp.x[k+1];

            // Check if the cubic interpolant from this grid point to the next
            // has any internal minima
            if(argc==2&&det>=0) {
                det=sqrt(det);
                double xplus=0.5*(-b+det),xminus=0.5*(-b-det),fv;
                if(0<xplus&&xplus<a) {
                    xplus/=a;
                    fv=hermite(xplus,t0,t1,t2,t3);
                    if(fv<fmin) {fmin=fv;xmin=(1-xplus)*gp.x[k]+xplus*gp.x[k+1];}
                }
                if(0<xminus&&xminus<a) {
                    xminus/=a;
                    fv=hermite(xminus,t0,t1,t2,t3);
                    if(fv<fmin) {fmin=fv;xmin=(1-xminus)*gp.x[k]+xminus*gp.x[k+1];}
                }
            }
        }

        // Check for lift-off
        if(hsearch&&fmin<gfmin) {
            gfmin=fmin;gimin=i;

            // Calculate the time origin via the method of the Kolinski et al.,
            // which uses the moment when the drop is experimentally detectable
            // at a height of h_ev
            if(efield&&fmin<h_ev) {
                tmin0_ev=i*fdt-xmin*xmin/(2*gls->R*gls->V)+h_ev/gls->V;
                efield=false;
            }
        }

        // Print the minimum that was found to the output file
        fprintf(of,"%d %g %.10g %.10g\n",i,i*fdt,xmin,fmin);
        if(fmin>gfmin*1.005) hsearch=false;
    }
    delete [] q;
    fclose(of);

    // Print information about liftoff
    double gtmin=gimin*fdt;
    printf("# %s: %d files processed\n%g %g %g %g %g\n",fn,i,fi.mu*fi.rhoinv*1e6,
           gfmin,gtmin-tmin0,gtmin-tmin0_ev,100.*float(gimin)/float(i-1));

}
