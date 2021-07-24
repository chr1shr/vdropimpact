#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <sys/stat.h>

#include "common.hh"
#include "fluid_2d.hh"

/** This source code file contains routines in the fluid_2d class that are
 * splash-specific or are related to input & output. */

/** Initializes the simulation fields. */
void fluid_2d::init_fields() {

    // Set the time and frame number
    f_num=0;time=0;

    // Set the velocity and pressure fields. Since the velocity field in the
    // inviscid boundaries is time-integrated it is initialized to be zero here
    // also.
#pragma omp parallel for
    for(field *fp=fbase;fp<fbase+(n+4)*ml;fp++) fp->clear_main();

    // Clear the pressure and shear stress arrays used in the gas layer
    // calculation
    gl->init_fields();

    // Initialize the tracers
    if(ntrace>0) init_tracers();

    // Now that the primary grid points are set up, initialize the ghost
    // points according to the boundary conditions
    set_boundaries(0);
}

/** Loads a restart file containing all of the simulation information.
 * \param[in] sn the frame number to read from. */
void fluid_2d::load_restart(int sn) {

    // Open the file and read the simulation time
    FILE *inf=odir_open("res",sn,"rb");
    safe_fread(&time,sizeof(double),1,inf,"simulation time");
    f_num=sn;

    // Read the primary fields
    double *lbuf=new double[3*ml],*lp,*le=lbuf+3*ml;
    field *fp=fbase;
    for(int j=0;j<n+4;j++) {
        safe_fread(lbuf,sizeof(double),3*ml,inf,"primary fields");
        lp=lbuf;
        while(lp<le) (fp++)->unpack(lp);
    }
    delete [] lbuf;

    // Load the tracers
    safe_fread(tm,sizeof(double),ntrace<<1,inf,"tracer information");

    // Load any required gas_layer fields
    gl->load_restart(time,inf);
    fclose(inf);
}

/** Searches for the last restart file available. */
void fluid_2d::load_last_restart() {
    if(restart_freq==0) puts("# Restart frequency not set, starting from frame 0");
    if(restart_exists(restart_freq)) {
        int lr=restart_freq;

        // See if any subsequent restart files exist
        while(restart_exists(lr+restart_freq)) {
            lr+=restart_freq;
            if(lr>=nframes) {
                fputs("# Final restart file found; nothing to do\n",stderr);
                exit(1);
            }
        }
        printf("# Starting from restart file at frame %d\n",lr);
        load_restart(lr);
        return;
    } else puts("# No restart file found, starting from frame 0");
    init_fields();
}

/** Checks to see a restart file exists.
 * \param[in] sn the frame number to check. */
inline bool fluid_2d::restart_exists(int sn) {
    char *bufc=reinterpret_cast<char*>(buf);
    sprintf(bufc,"%s/res.%d",filename,sn);
    struct stat st;
    return stat(bufc,&st)==0;
}

/** Saves a restart file. This contains all of the simulation information
 * required to restart from this point.
 * \param[in] sn the frame number to read from. */
void fluid_2d::save_restart(int sn) {

    // Open the file and output the simulation time
    FILE *outf=odir_open("res",sn,"wb");
    fwrite(&time,sizeof(double),1,outf);

    // Output the primary fields
    double *lbuf=new double[3*ml],*lp,*le=lbuf+3*ml;
    field *fp=fbase;
    for(int j=0;j<n+4;j++) {
        lp=lbuf;
        while(lp<le) (fp++)->pack(lp);
        fwrite(lbuf,sizeof(double),3*ml,outf);
    }
    delete [] lbuf;

    // Output the tracers
    fwrite(tm,sizeof(double),ntrace<<1,outf);

    // Output any required gas_layer fields
    gl->save_restart(outf);
    fclose(outf);
}

/** Outputs the tracer positions in a binary format that can be read by
 * Gnuplot. */
void fluid_2d::output_tracers(const int sn) {

    // Output the tracer positions in batches of 128 floats
    FILE *outf=odir_open("trace",sn,"wb");
    int j,tbatch=(ntrace<<1)/128,tres=(ntrace<<1)%128;
    float *fp,*fe=buf+128;
    double *tp=tm,*te=tm+(ntrace<<1);
    for(j=0;j<tbatch;j++) {
        fp=buf;
        while(fp<fe) *(fp++)=*(tp++);
        fwrite(buf,sizeof(float),128,outf);
    }

    // Output the remaining tracer positions, if any
    if(tres>0) {
        fp=buf;
        do {*(fp++)=*(tp++);} while(tp<te);
        fwrite(buf,sizeof(float),tres,outf);
    }

    // Close the file
    fclose(outf);
}

/** Writes a selection of simulation fields to the output directory.
 * \param[in] k the frame number to append to the output. */
void fluid_2d::write_files(int k) {
    if(fflags&1) output("u",0,k,true);
    if(fflags&2) output("v",1,k,true);
    if(fflags&4) output("p",2,k);
    if(fflags&8) output("w",3,k);
    if(fflags&16) gl->output_profile(filename,buf,k,ay,true);
    if(fflags&32) gl->output_profile(filename,buf,k,ay,false);
    if(fflags&64) gl->output_profile_double(filename,buf,k,true);
    if(fflags&128) gl->output_profile_double(filename,buf,k,false);
    if(ntrace>0) output_tracers(k);
    if(restart_freq>0&&k>0&&k%restart_freq==0) save_restart(k);
}

/** Outputs a 2D array to a file in a format that can be read by Gnuplot.
 * \param[in] prefix the field name to use as the filename prefix.
 * \param[in] mode the code of the field to print.
 * \param[in] sn the current frame number to append to the filename.
 * \param[in] ghost whether to output the ghost regions or not. */
void fluid_2d::output(const char *prefix,const int mode,const int sn,const bool ghost) {

    // Determine whether to output a cell-centered field or not
    bool cen=mode>=0&&mode<=1;
    int li=ghost?ml:(cen?m:m+1),
        lj=ghost?n+4:(cen?n:n+1);
    double disp=(cen?0.5:0)-(ghost?2:0);

    // Output the first line of the file
    FILE *outf=odir_open(prefix,sn,"wb");
    int i,j;
    float *bp=buf+1,*be=bp+li;
    *buf=li;
    for(i=0;i<li;i++) *(bp++)=ax+(i+disp)*dx;
    fwrite(buf,sizeof(float),li+1,outf);

    // Output the field values to the file
    field *fr=ghost?fbase:fm;
    for(j=0;j<lj;j++,fr+=ml) {
        field *fp=fr;
        *buf=ay+(j+disp)*dy;bp=buf+1;
        switch(mode) {
            case 0: while(bp<be) *(bp++)=(fp++)->u;break;
            case 1: while(bp<be) *(bp++)=(fp++)->v;break;
            case 2: while(bp<be) *(bp++)=(fp++)->p;break;
            case 3: while(bp<be) *(bp++)=vorticity(fp++);
        }
        fwrite(buf,sizeof(float),li+1,outf);
    }

    // Close the file
    fclose(outf);
}

/** Modifies the vertical velocity in the simulation to adjust the frame in
 * use. */
void fluid_2d::adjust_frame() {
    field *fb=fm-ml;
    double Vaccdt;

    // Mode 1: use the velocity under the droplet. These formulae are based
    // a Lagrange interpolant of the four cell-centered grid points
    if(nif_mode==1) {
        int m2=m/2;
        Vaccdt=x_sym?0.5625*(fb->v+fm->v)-0.0625*(fb[1].v+fm[1].v)
                    :0.28125*(fb[m2-1].v+fb[m2].v+fm[m2-1].v+fm[m2].v)
                    -0.03125*(fb[m2-2].v+fb[m2+1].v+fm[m2-2].v+fm[m2+1].v);

    // Mode 2: use an average velocity over a range
    } else {
        Vaccdt=0;
        for(int i=anif;i<bnif;i++) Vaccdt+=fm[i].v+fb[i].v;
        Vaccdt*=0.5/(bnif-anif);
    }

    // Adjust the vertical velocity field
#pragma omp parallel for
    for(field *fp=fbase;fp<fbase+(n+4)*ml;fp++) fp->v-=Vaccdt;

    // Adjust the frame velocity to match
    gl->Vframe-=Vaccdt;
}

/** Sets the fields in the ghost regions where the effect of viscosity is small. */
void fluid_2d::set_inviscid_boundaries(double dt) {
    double prefac=dx*rhoinv*dt/(3*M_PI);

    // Set right ghost values and left ghost values if x_sym==false.
#pragma omp parallel for
    for(int j=0;j<n;j++) {
        field *fp=fm+j*ml+m;

        // Set inner ghost layer on right
        simpsons_velocity_bc(bx+dx*0.5,ay+(j+0.5)*dy,prefac,fp->u,fp->v);

        // Set outer ghost layer on right
        simpsons_velocity_bc(bx+dx*1.5,ay+(j+0.5)*dy,prefac,fp[1].u,fp[1].v);

        // For symmetric case, copy the values to the left side
        if(!x_sym) {
            fp[-m-1].u=-fp->u;
            fp[-m-1].v=fp->v;
            fp[-m-2].u=-fp[1].u;
            fp[-m-2].v=fp[1].v;
        }
    }

    // Set top ghost values
#pragma omp parallel for
    for(int i=0;i<ml;i++) {
        field *fp=fbase+(n+2)*ml+i;

        // Set inner ghost layer on top
        simpsons_velocity_bc(ax+(i-1.5)*dx,by+dy*0.5,prefac,fp->u,fp->v);

        // Set outer ghost layer on top
        simpsons_velocity_bc(ax+(i-1.5)*dx,by+dy*1.5,prefac,fp[ml].u,fp[ml].v);
    }
}

/** Sets the fields in the ghost regions where the effect of viscosity is to be
 * understood. */
void fluid_2d::set_viscous_boundaries() {

    /** Set left boundary if x_sym==true */
    for(int j=0;j<n;j++) {
        field *fp=fm+j*ml;
        if(x_sym) {
            fp[-1].u=-fp->u;
            fp[-1].v=fp->v;
            fp[-2].u=-fp[1].u;
            fp[-2].v=fp[1].v;
        }
    }

    // Set bottom horizontal velocity values according to the tangential stress
    // condition for the liquid at the liquid-air interface
    int i=0;
    for(field *fp=fbase,*fe=fp+ml;fp<fe;fp++,i++) {

        // Compute du/dy at y=0
        double un_edge=muinv*gl->tau[i];

        // Set the inner ghost layer
        fp[ml].u=fp[2*ml].u-dy*un_edge;

        // Set the outer ghost layer
        fp->u=fp[3*ml].u-3*dy*un_edge;
    }

    // Set bottom vertical velocity using incompressibility
    for(field *fp=fbase+1,*fe=fp+m+1;fp<fe;fp++) {

        // Calculate the velocity derivatives, since the current values are for
        // the previous time step. NOTE: it's not clear the
        // monotonicity-limited derivative is useful here, and this could be
        // replaced with centered differences.
        fp[ml].ux=0.5*xsp*(fp[ml+1].u-fp[ml-1].u);
        fp->ux=0.5*xsp*(fp[1].u-fp[-1].u);

        // Calculate v from incompressibility
        fp[ml].v=fp[2*ml].v+dy*fp[ml].ux;
        fp->v=fp[ml].v+dy*fp->ux;
    }

//    // Optional smoothing of ghost values at the extreme edge of the domain
//    if(x_sym) {
//        int l=9*m/10;
//        smooth_line(fbase+2+l,m-l+1);
//        smooth_line(fbase+2+l+ml,m-l+1);
//    }
}

/** Smooths a horizontal line of the velocity field.
 * \param[in] fp a pointer to the start of the line.
 * \param[in] m the number of grid points in the line to smooth. */
void fluid_2d::smooth_line(field *fp,int d) {
    double *R=new double[2*d],*Rp=R;
    for(int i=0;i<d;i++) {
        *(Rp++)=fp[i].u;
        *(Rp++)=fp[i].v;
    }
    Rp=R;
    for(int i=1;i<d-1;i++,Rp+=2) {
        fp[i].u=0.25*(*Rp+Rp[4])+0.5*Rp[2];
        fp[i].v=0.25*(Rp[1]+Rp[5])+0.5*Rp[3];
    }
    delete [] R;
}

/** Carries out composite Simpson's rule to compute u at a specified location (x,y).
 * \param[in] (x,y) the point at which u is computed.
 * \param[in] prefac the prefactor to apply.
 * \param[in,out] (u,v) the velocities to adjust. */
void fluid_2d::simpsons_velocity_bc(double x,double y,double prefac,double &u,double &v) {
    double su,sv,fac,fsu,fsv,xms=x-ax,xps=x+ax,ysq=y*y,*pp=gl->p;

    // Calculate the approximate integral
    simpsons_terms(xms,xps,ysq,*pp,su,sv);
    xms-=dx;xps+=dx;pp++;
    for(int j=1;j<me-1;j++,xms-=dx,xps+=dx,pp++) {
        simpsons_terms(xms,xps,ysq,*pp,fsu,fsv);
        fac=(j&1)?4.:2.;
        su+=fac*fsu;sv+=fac*fsv;
    }
    simpsons_terms(xms,xps,ysq,*pp,fsu,fsv);
    su+=fsu;sv+=fsv;

    /** Forward Euler step to evolve the velocities. The Simpson's integration
     * result is multiplied by the prefactor containing the integral factor
     * dx/3. An additional factor of 2y is incorporated into the u calculation.*/
    u+=2*y*prefac*su;
    v+=prefac*sv;
}

/** Calculates the individual function values that are used in the composite
 * Simpson's rule for the velocity at the position (x,y).
 * \param[in] (xms,xps) the displacements between x, and s and -s where the
 *                      integrand is evaluated.
 * \param[in] ysq the square of the y coordinate.
 * \param[in] p the pressure.
 * \param[out] (fsu,fsv) the individual terms that are used in the Simpson's rule
 *                       integration. */
void fluid_2d::simpsons_terms(double xms,double xps,double ysq,double p,double &fsu,double &fsv) {
    double rm=xms*xms+ysq;rm=1/(rm*rm);
    if(x_sym) {
        double rp=xps*xps+ysq;rp=1/(rp*rp);
        fsu=p*(xms*rm+xps*rp);
        fsv=-p*((xms*xms-ysq)*rm+(xps*xps-ysq)*rp);
    } else {
        fsu=p*xms*rm;
        fsv=-p*(xms*xms-ysq)*rm;
    }
}
