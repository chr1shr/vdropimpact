#include "fileinfo.hh"
#include "gas_layer_data.hh"
#include "gas_layer_solve.hh"

#include <cmath>

/** The class constructor reads parameters from an input file, converts some
 * physical units into simulation units, and performs a number of consistency
 * checks.
 * \param[in] infile the name of the input file to read. */
fileinfo::fileinfo(const char* infile) : m(-1), fflags(0), ntrace(0),
    nframes(-1), restart_freq(0), nif_mode(0), implicit_visc(false),
    x_sym(false), mr_time_output(false), ay(0.), tmult(-1) {
    double L_nd=-1,h0_nd=-1,t_end_nd=-1,nul=-1,R=-1,V=-1,
           alpha=-1,sigma=-1,mug=-1,Pamb=-1,nif_r=-1;
    int m_dat=-1,n_dat=-1;
    char *path_pressure=NULL,*path_tau=NULL,glo=0;

    // Check that the filename ends in '.cfg'
    int l=strlen(infile),ln=1;
    if(l<4) fatal_error("Filename is too short",1);
    const char* ip=infile+l-4;
    if(*ip!='.'||ip[1]!='c'||ip[2]!='f'||ip[3]!='g')
        fatal_error("Filename must end in '.cfg'",1);

    // Assemble output filename by replacing '.cfg' with '.odr'
    filename=new char[l+1];
    memcpy(filename,infile,l-3);
    char *fp=filename+l-3;
    *fp='o';fp[1]='d';fp[2]='r';fp[3]=0;

    // Open the input file and read
    FILE *f=safe_fopen(infile,"r");
    char *buf=new char[fileinfo_buf_size],*bp;
    while(!feof(f)) {
        if(fgets(buf,fileinfo_buf_size,f)==NULL) break;

        // Locate comments and remove by replacing comment character
        // with a null character
        bp=buf;
        while((*bp)!=0) {
            if(*bp=='#') {*bp=0;break;}
            bp++;
        }

        // Separate entries in the file by tabs and spaces; if no entries are
        // available then skip this line
        bp=strtok(buf," \t\n");
        if(bp==NULL) {ln++;continue;}

        // Look for a known keyword, and read in any extra values
        if(se(bp,"grid_points")) {
            m=atoi(next_token(ln));
            n=final_int(ln);
            if(m<=0||n<=0) fatal_error("Grid dimensions must be positive",1);
        } else if(se(bp,"L_nd")) L_nd=final_double(ln);
        else if(se(bp,"h0_nd")) h0_nd=final_double(ln);
        else if(se(bp,"t_end_nd")) t_end_nd=final_double(ln);
        else if(se(bp,"tracers")) ntrace=final_int(ln);
        else if(se(bp,"frames")) nframes=final_int(ln);
        else if(se(bp,"tmult")) tmult=final_double(ln);
        else if(se(bp,"restart_freq")) restart_freq=final_int(ln);
        else if(se(bp,"nul")) nul=final_double(ln);
        else if(se(bp,"nul_cSt")) nul=1e-6*final_double(ln);
        else if(se(bp,"rhol")) rho=final_double(ln);
        else if(se(bp,"R")) R=final_double(ln);
        else if(se(bp,"V")) V=final_double(ln);
        else if(se(bp,"gamma")) alpha=1./final_double(ln);
        else if(se(bp,"alpha")) alpha=final_double(ln);
        else if(se(bp,"sigma")) sigma=final_double(ln);
        else if(se(bp,"mug")) mug=final_double(ln);
        else if(se(bp,"Pamb")) Pamb=final_double(ln);
        else if(se(bp,"x_sym")) {
            x_sym=true;
            check_no_more(ln);
        } else if(se(bp,"implicit_visc")) {
            implicit_visc=true;
            check_no_more(ln);
        } else if(se(bp,"mr_time_output")) {
            mr_time_output=true;
            check_no_more(ln);
        } else if(se(bp,"gas_layer_data")) {
            glo=1;
            m_dat=atoi(next_token(ln));
            n_dat=atoi(next_token(ln));
            char *p=next_token(ln);
            path_pressure=new char[strlen(p)+1];
            strcpy(path_pressure,p);
            p=next_token(ln);
            path_tau=new char[strlen(p)+1];
            strcpy(path_tau,p);
            check_no_more(ln);
        } else if(se(bp,"gas_layer_model")) {
            glo=2;
            check_no_more(ln);
        } else if(se(bp,"nif_center")) {
            nif_mode=1;
            check_no_more(ln);
        } else if(se(bp,"nif_range")) {
            nif_mode=2;
            nif_r=final_double(ln);
        } else if(se(bp,"output")) {
            fflags=0;
            bp=next_token(ln);
            while(bp!=NULL) {
                if(se(bp,"u")) fflags|=1;
                else if(se(bp,"v")) fflags|=2;
                else if(se(bp,"p")) fflags|=4;
                else if(se(bp,"pb")) fflags|=8;
                else if(se(bp,"pbfull")) fflags|=16;
                else if(se(bp,"w")) fflags|=32;
                else if(se(bp,"h")) fflags|=64;
                else if(se(bp,"fbd")) fflags|=128;
                else if(se(bp,"pg")) fflags|=256;
                else if(se(bp,"hfull")) fflags|=512;
                else if(se(bp,"fbdfull")) fflags|=1024;
                else if(se(bp,"pgfull")) fflags|=2048;
                else fatal_error("Output type not understood",1);
                bp=strtok(NULL," \t\n");
            }
        } else {
            fprintf(stderr,"Keyword '%s' not understood on line %d\n",bp,ln);
            exit(1);
        }
        ln++;
    }
    delete [] buf;
    fclose(f);

    // Check that all basic class parameters have been specified
    if(m==-1) fatal_error("Grid dimensions not set",1);
    check_invalid(L_nd,"L_nd");
    check_invalid(h0_nd,"h0_nd");
    check_invalid(t_end_nd,"t_end_nd");
    check_invalid(tmult,"tmult");
    check_invalid(nframes,"frames");
    check_invalid(restart_freq,"restart_freq");

    // Check the viscosity and density, and calculate class constants
    check_invalid(nul,"nul");
    check_invalid(rho,"rhol");
    mu=(nul*rho);muinv=1./mu;rhoinv=1./rho;

    // Calculate dimensional quantities from the supplied non-dimensional ones
    double st=mug/(rho*V*R),st3=pow(st,1./3.),st23=st3*st3,
           L=L_nd*R*st3,h0=h0_nd*R*st23;
    t_end=t_end_nd*R*st23/V;

    // Set grid dimensions
    bx=L;by=L*n/m;
    if(x_sym) ax=0;
    else ax=-L,by*=2;

    // Check that all intermediate quantities are valid
    if(glo==0) fatal_error("Either gas_layer_data or gas_layer_model must be set",1);
    else if(glo==1) {

        // If we are using data in the gas layer, then allocate the
        // gas_layer_data object and pass in the data filenames
        gl=new gas_layer_data(m,ax,bx,m_dat,n_dat,t_end,path_pressure,path_tau);
        delete [] path_tau;
        delete [] path_pressure;
    } else {

        // If we are solving in the gas layer, then check that all of the
        // parameters controlling the gas layer have been specified
        check_invalid(alpha,"alpha");
        check_invalid(mug,"mug");
        check_invalid(Pamb,"Pamb");
        check_invalid(R,"R");
        check_invalid(V,"V");
        check_invalid(sigma,"sigma");
        gl_params glp(alpha,mug,Pamb,h0,R,V,sigma);
        gl=new gas_layer_solve(m,ax,bx,glp,x_sym);
    }

    // If the non-inertial frame using a velocity range is in use, then calculate
    // the range of grid points that
    if(nif_mode==2) {

        // Compute the index ranges
        bnif=int((nif_r-ax)*xsp+0.5);
        if(bnif>=m) bnif=m-1;
        if(x_sym) anif=0;
        else {
            anif=int((-nif_r-ax)*xsp+0.5);
            if(anif<0) anif=0;
        }

        // Check that the velocity window is has at least one grid point in it
        if(bnif-anif<=0) fatal_error("Velocity window for non-inertial frame is too small",1);
    }
}

/** Checks that there are no subsequent values.
 * \param[in] ln the current line number. */
void fileinfo::check_no_more(int ln) {
    if(strtok(NULL," \t\n")!=NULL) {
        fprintf(stderr,"Too many arguments at input line %d\n",ln);
        exit(1);
    }
}

/** Finds the next token in a string, and if none is availble, gives an error
 * message.
 * \param[in] ln the current line number. */
char* fileinfo::next_token(int ln) {
    char *temp=strtok(NULL," \t\n");
    if(temp==NULL) {
        fprintf(stderr,"Not enough arguments at input line %d\n",ln);
        exit(1);
    }
    return temp;
}

/** Checks that a parameter is a valid positive value.
 * \param[in] val the value to check.
 * \param[in] p the name of the value. */
void fileinfo::check_invalid(double val,const char *p) {
    if(val<0) {
        fprintf(stderr,"Value of %s either invalid or not set\n",p);
        exit(1);
    }
}
