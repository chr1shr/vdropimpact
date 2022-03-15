#include "common.hh"
#include "gas_layer_solve.hh"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

/** The gas layer equation solver initializes constants allocates memory for
 * the fluid height and pressure.
 * \param[in] m_ the number of gridpoints.
 * \param[in] (ax_,bx_) the horizontal coordinate range to compute over.
 * \param[in] glp the physical parameters for the gas layer.
 * \param[in] x_sym_ whether to use the symmetry boundary conditions. */
gas_layer_solve::gas_layer_solve(int m_,double ax_,double bx_,gl_params &glp,bool x_sym_) :
    gas_layer(m_,ax_,bx_,glp.V), gl_params(glp), x_sym(x_sym_), xsp(1./dx),
    xxsp(xsp*xsp), fbd(new double[m]), h(new double[m]), hx(new double[me]),
    ht(new double[m]), pg_old(new double[me]), pg(new double[me]), t(me) {}

/** The class destructor frees the dynamically allocated memory. */
gas_layer_solve::~gas_layer_solve() {
    delete [] pg;
    delete [] pg_old;
    delete [] ht;
    delete [] hx;
    delete [] h;
    delete [] fbd;
}

/** Initializes the interface height and velocity based on the drop height and
 * radius and the ambient pressure. */
void gas_layer_solve::init_fields() {

    // Initialize the array to track flow through the bottom boundary
    for(int i=0;i<m;i++) fbd[i]=0;

    // Initialize the drop height and x-derivative
    for(int i=0;i<m;i++){
        double x=ax+(i+0.5)*dx;
        h[i]=h0+0.5/R*x*x;
    }
    calculate_hx();

    // Initialize the drop velocity
    for(int i=0;i<m;i++) ht[i]=-V;

    // Initially specify the pressure to be Pamb everywhere. Since the code
    // solves for P-Pamb, initialize the pressure array to be zero.
    for(int i=0;i<me;i++) pg[i]=0;

    // Initialize the liquid pressure and shear stress
    calculate_liquid_fields();
}

/** Calculates the derivative of the height field. The derivative is stored at
 * cell edges, and so this function uses a second-order centered
 * finite-difference formula. */
void gas_layer_solve::calculate_hx() {
    for(int i=1;i<m;i++) hx[i]=xsp*(h[i]-h[i-1]);
    *hx=x_sym?0:hx[1];
    hx[m]=hx[m-1];
}

/** Performs a Newton step to solve the nonlinear pressure equation.
 * \param[in] dtinv the inverse timestep.
 * \return The average square residual before the Newton step was performed. */
double gas_layer_solve::newton_step(double dtinv) {
    double br,ressq=0.,*a=t.a,*b=t.b,*c=t.c,*d=t.d,px,pxx,hsq,hcu,h_ns,ht_ns;

    // Assemble the matrix, first setting the terms corresponding to the
    // boundaries. Note that an extra minus sign appears on the d terms due
    // to the minus in the standard Newton step.

    // Apply left boundary condition at x=0 -- from logbook on 14/1/18.
    double hcu_0=h[0]*h[0]*h[0];
    *b=(x_sym?
        -12*mu*(alpha*dtinv*h[0]+ht[0])+2*xxsp*hcu_0*(pg[1]-pg[0])
        -2*xxsp*hcu_0*pg[0]-2*xxsp*hcu_0*Pamb
        :1);
    *c=(x_sym?
        2*xxsp*hcu_0*pg[0]+2*xxsp*hcu_0*Pamb
        :0);
    *d=(x_sym?
        12*mu*(alpha*dtinv*h[0]*(pg[0]-pg_old[0])+ht[0]*pg[0]+ht[0]*Pamb)
        -2*xxsp*hcu_0*pg[0]*(pg[1]-pg[0]) -2*xxsp*hcu_0*Pamb*(pg[1]-pg[0])
        :-*pg);
    a[me-1]=0;b[me-1]=1;d[me-1]=-pg[me-1];

    // Loop over the other lines; update at the interior grid points
    for(int i=x_sym?0:1;i<me-1;i++) {

        // Compute h and ht in staggered grid
        h_ns=i>0?0.5*(h[i-1]+h[i]):h[i];
        ht_ns=i>0?0.5*(ht[i-1]+ht[i]):ht[i];

        // Compute useful constants: the centered difference of pressure, the
        // height squared, and the height cubed
        px=i>0?0.5*xsp*(pg[i+1]-pg[i-1]):0;
        pxx=i>0?xxsp*(pg[i+1]-2*pg[i]+pg[i-1])
               :xxsp*(pg[i+1]-2*pg[i]+pg[i+1]);
        hsq=h_ns*h_ns;
        hcu=hsq*h_ns;

        // Compute the large bracketed term in the equation
        br=3*hsq*hx[i]*px
          +hcu*pxx
          -12*mu*alpha*h_ns*dtinv
          -12*mu*ht_ns;

        // Compute the matrix entries. These remain unchanged in the staggered
        // grid, since hcu, hsq and hx are already on the same grid as the
        // pressure.
        a[i]=-xsp*alpha*hcu*px+pg[i]*(-1.5*xsp*hsq*hx[i]+hcu*xxsp)+Pamb*(-1.5*xsp*hsq*hx[i]+hcu*xxsp);
        b[i]=br-2*pg[i]*hcu*xxsp-2*Pamb*hcu*xxsp;
        c[i]= xsp*alpha*hcu*px+pg[i]*( 1.5*xsp*hsq*hx[i]+hcu*xxsp)+Pamb*( 1.5*xsp*hsq*hx[i]+hcu*xxsp);
        if(i==0) c[i]+=a[i];
        d[i]=-(
                12*mu*alpha*h_ns*dtinv*pg_old[i]
                +alpha*hcu*px*px+pg[i]*br
                +Pamb*(-12*mu*ht_ns+3*hsq*hx[i]*px+hcu*pxx));
        ressq+=d[i]*d[i];
    }

    // Solve using the tridiagonal solver, and add the Newton step to the gas
    // pressure. Return the square residual.
    t.solve_lapack();
    for(int i=0;i<me;i++) pg[i]+=d[i];
    return ressq/me;
}

/** Updates the droplet height.
 * \param[in] dt the timestep to use. */
void gas_layer_solve::update_height(double dt) {

    // Compute the drop height at the new time
    for(int i=0;i<m;i++) {

        // Compute the horizontal velocity. This also incorporates a factor of
        // 1/(2 dx), which is required to properly scale the output of the eno2
        // function.
        double uc=0.25*xsp*(fm[i].u+fb[i].u),

        // Compute the vertical velocity
               vc=0.5*(fm[i].v+fb[i].v),

        // Compute the x derivative of h using the ENO2 method
               dhdx=uc>0?eno2(h_ref(i+1),h[i],h_ref(i-1),h_ref(i-2))
                        :-eno2(h_ref(i-1),h[i],h_ref(i+1),h_ref(i+2));

        // Compute the partial t derivative of h
        ht[i]=-uc*dhdx+vc-Vframe;

        // Update the height
        h[i]+=dt*ht[i];

        // Update the flow across the boundary
        fbd[i]+=dt*vc;
    }

    // Update the x-derivative of the height field
    calculate_hx();
}

/** Calculates the ENO derivative using a sequence of values at
 * four gridpoints.
 * \param[in] (p0,p1,p2,p3) the sequence of values to use.
 * \return The computed derivative. */
inline double gas_layer_solve::eno2(double p0,double p1,double p2,double p3) {
    return fabs(p0-2*p1+p2)>fabs(p1-2*p2+p3)?3*p1-4*p2+p3:p0-p2;
}

/** Calculate the height, interpreting out-of-bounds values using linear
 * interpolation.
 * \param[in] i the index at which to evaluate the height, between -2 and m+1.
 * \return The height. */
inline double gas_layer_solve::h_ref(int i) {
    return i>=0?(i<m?h[i]:i==m?2*h[m-1]-h[m-2]:3*h[m-1]-2*h[m-2])
               :(x_sym?(i==-1?*h:h[1])
                       :(i==-1?2*(*h)-h[1]:3*(*h)-2*h[1]));
}

/** Updates the gas pressure.
 * \param[in] dt the timestep to use. */
void gas_layer_solve::update_gas_pressure(double dt) {

    // Copy the previous pressure into the pg_old array, since it appears in
    // the equation that must be solved (apart from in the special case when
    // alpha is zero)
    memcpy(pg_old,pg,me*sizeof(double));

    // Step forward the pressure using the implicit timestepper
    implicit_step(dt);
}

/** Performs an implicit step of the pressure in the gas layer.
 * \param[in] dt the timestep to use. */
void gas_layer_solve::implicit_step(double dt) {
    double dtinv=1./dt;
    int nstep=0,naccept=0;
    do {

        // Check for too many Newton steps, indicating non-convergence of the
        // nonlinear rootfinding procedure
        if(++nstep>max_newton_steps) fatal_error("Too many Newton steps\n",1);

        // Perform a Newton step, and check if the error is below the
        // threshold. Finish the routine once two successive steps were
        // acceptable.
        if(newton_step(dtinv)<newton_thresh) naccept++;
        else naccept=0;
    } while(naccept<2);
}

/** A function to compute the fields at the next time step. */
void gas_layer_solve::calc_p_and_tau(double t,double dt) {

    // Update the height of the liquid droplet
    update_height(dt);

    // Solve for the pressure in the gas using gas_layer_solve
    update_gas_pressure(dt);

    // Calculate the liquid pressure (by adding in surface tension) and the
    // shear stress
    calculate_liquid_fields();
}

/** Computes the liquid pressure from the gas pressure and surface tension. */
void gas_layer_solve::calculate_liquid_fields() {
    const double fac=0.5*sigma*xxsp,fac2=sigma/R;

    // Set the liquid pressure on the left edge. See logbook on 07/02/18 for
    // derivation of the symmetric case.
    if(x_sym) {
        *p=*pg+2*fac*(h[1]-*h)-fac2;
        p[1]=pg[1]+fac*(h[2]-h[1])-fac2;
    } else {
        *p=*pg;p[1]=pg[1];
    }

    // Set the liquid pressure in the interior of the grid
    for(int i=2;i<m-1;i++) p[i]=pg[i]+fac*(h[i+1]-h[i]-h[i-1]+h[i-2])-fac2;

    // Set the liquid pressure on the right edge
    p[m-1]=pg[m-1];p[m]=pg[m];

    // Compute the shear stress on the bottom boundary
    for(int i=0;i<m;i++) tau[i]=-0.5*xsp*mu*(p[i+1]-p[i])*h[i];
}

/** Outputs a linear profile in the Gnuplot matrix binary format.
 * \param[in] filename the filename of the output directory.
 * \param[in] buf a pointer to temporary space to assemble the output filename.
 * \param[in] sn the current frame number to append to the filename.
 * \param[in] ay the y position of the gas layer.
 * \param[in] fld which field to output (0: height, 1: flow across the
 *                boundary, 2: gas pressure relative to ambient pressure). */
void gas_layer_solve::output_profile(const char *filename,float *buf,int sn,double ay,int fld) {

    // Assemble the output filename and open the output file
    char *bufc=reinterpret_cast<char*>(buf);
    sprintf(bufc,fld==0?"%s/height.%d":
                (fld==1?"%s/fbd.%d":"%s/pg.%d"),filename,sn);
    FILE *outf=safe_fopen(bufc,"wb");

    // Output the first line of the file
    int l=fld==2?me:m;
    float *bp=buf+1,*be=bp+l;
    *buf=l;
    for(int i=0;i<l;i++) *(bp++)=ax+(i+(fld==2?0:0.5))*dx;
    fwrite(buf,sizeof(float),l+1,outf);

    // Output the field values to the file
    *buf=ay;bp=buf+1;
    double *hp=fld==0?h:(fld==1?fbd:pg);
    while(bp<be) *(bp++)=*(hp++);
    fwrite(buf,sizeof(float),l+1,outf);

    // Close the file
    fclose(outf);
}

/** Outputs a linear profile in a custom binary format at full double
 * precision.
 * \param[in] filename the filename of the output directory.
 * \param[in] buf a pointer to temporary space to assemble the output filename.
 * \param[in] sn the current frame number to append to the filename.
 * \param[in] fld which field to output (0: height, 1: flow across the
 *                boundary, 2: gas pressure relative to ambient pressure). */
void gas_layer_solve::output_profile_double(const char *filename,float *buf,int sn,int fld) {

    // Assemble the output filename and open the output file
    char *bufc=reinterpret_cast<char*>(buf);
    sprintf(bufc,fld==0?"%s/hfull.%d":
                (fld==1?"%s/fbdfull.%d":"%s/pgfull.%d"),filename,sn);
    FILE *outf=safe_fopen(bufc,"wb");

    // Output coordinate information
    int l=fld==2?me:m;
    fwrite(&l,sizeof(int),1,outf);
    fwrite(&ax,sizeof(double),2,outf);

    // Output the field values
    fwrite(fld==0?h:(fld==1?fbd:pg),sizeof(double),l,outf);

    // Close the file
    fclose(outf);
}

/** Loads the required fields from a restart file.
 * \param[in] time the current simulation time.
 * \param[in] inf a file handle to read from. */
void gas_layer_solve::load_restart_internal(double time,FILE *inf) {
    safe_fread(p,sizeof(double),me,inf,"p array");
    safe_fread(tau+2,sizeof(double),m,inf,"tau array");
    apply_tau_bc();
    safe_fread(h,sizeof(double),m,inf,"h array");
    calculate_hx();
    safe_fread(ht,sizeof(double),m,inf,"ht array");
    safe_fread(pg,sizeof(double),me,inf,"pg array");
}

/** Loads the required fields from a restart file.
 * \param[in] time the current simulation time.
 * \param[in] outf a file handle to read from. */
void gas_layer_solve::save_restart_internal(FILE *outf) {
    fwrite(p,sizeof(double),me,outf);
    fwrite(tau+2,sizeof(double),m,outf);
    fwrite(h,sizeof(double),m,outf);
    fwrite(ht,sizeof(double),m,outf);
    fwrite(pg,sizeof(double),me,outf);
}
