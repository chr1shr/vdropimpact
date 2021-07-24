#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

#include "common.hh"
#include "fluid_2d.hh"
#include "gas_layer.hh"

/** The class constructor sets up constants the control the geometry and the
 * simulation, dynamically allocates memory for the fields, and calls the
 * routine to initialize the fields.
 * \param[in] infile the name of the configuration file. */
fluid_2d::fluid_2d(const char* infile) : fileinfo(infile),
    mn(m*n), me(m+1), ne(n+1), ml(m+4), dx((bx-ax)/m), dy((by-ay)/n),
    xsp(1./dx), ysp(1./dy), xxsp(xsp*xsp), yysp(ysp*ysp),
    dt_reg(implicit_visc?dx*tmult:dx*dx*rho*muinv*tmult*0.25),
    fbase(new field[ml*(n+4)]), fm(fbase+2*ml+2), tm(new double[ntrace<<1]),
    src(new double[me*ne]), s_bc(src+mn), ms_mac(*this), ms_fem(*this),
    vix(implicit_visc?new visco_impl(*this,true):NULL),
    viy(implicit_visc?new visco_impl(*this,false):NULL),
    buf(new float[m>123?m+5:128]) {
    gl->set_field_pointer(fm);
}

/** The class destructor frees the dynamically allocated memory. */
fluid_2d::~fluid_2d() {
    delete [] buf;
    if(viy!=NULL) {
        delete viy;
        delete vix;
    }
    delete [] src;
    delete [] tm;
    delete [] fbase;
}

/** Carries out the simulation for a specified time interval using the direct
 * simulation method, periodically saving the output. */
void fluid_2d::solve() {
    double t0,t1,t2,adt;
    int l=timestep_select(t_end/nframes,adt);t_gl=t_bc=0;

    // Save header file
    if(f_num==0) {
        write_files(0);
        puts("# Output frame 0");
    }
    t0=wtime();

    // Loop over the output frames
    for(int k=f_num+1;k<=nframes;k++) {

        // Perform the simulation steps
        for(int j=0;j<l;j++) step_forward(adt);

        // Output the fields
        t1=wtime();
        write_files(k);

        // Print diagnostic information
        t2=wtime();
        if(mr_time_output) {

            // Print machine-readable timing information
            printf("%d %d %.8g %.8g  %d %.8g %d %.8g",k,l,t1-t0,t2-t1,
                   ms_mac.tp.vcount_extra(),ms_mac.wc_time,
                   ms_fem.tp.vcount_extra(),ms_fem.wc_time);
            ms_mac.reset_counters();
            ms_fem.reset_counters();
            if(vix==NULL) printf("  %.8g %.8g\n",t_gl,t_bc);
            else {
                printf(" %d %.8g %d %.8g  %.8g %.8g\n",vix->tp.vcount_extra(),vix->wc_time,
                        viy->tp.vcount_extra(),viy->wc_time,t_gl,t_bc);
                vix->reset_counters();
                viy->reset_counters();
            }
            t_gl=t_bc=0;
        } else {

            // Print human-readable timing information
            printf("# Output frame %d [%d, %.8g s, %.8g s] {MAC %.2f, FEM %.2f",
                   k,l,t1-t0,t2-t1,ms_mac.tp.avg_iters(),ms_fem.tp.avg_iters());
            if(vix==NULL) puts("}");
            else printf(", VX %.2f, VY %.2f}\n",vix->tp.avg_iters(),viy->tp.avg_iters());
        }
        t0=t2;
    }
}

/** Steps the simulation fields forward.
 * \param[in] dt the time step to use. */
void fluid_2d::step_forward(double dt) {
    int j;

    // Move the tracers. Step the gas pressure and shear stress forward. Update
    // the time to the half-timestep.
    update_tracers(dt);
    double t0=wtime();
    gl->update(time,dt);
    t_gl+=wtime()-t0;
    time+=0.5*dt;

    // Do the pre-computation step to evaluate the tangential derivatives
#pragma omp parallel for
    for(j=0;j<n;j++) for(int i=0;i<m;i++) {

        // Create references to the fields in the neighboring gridpoints
        field *fp=fm+(ml*j+i),&f=*fp,&fl=fp[-1],&fr=fp[1],&fd=fp[-ml],&fu=fp[ml];

        // Calculate monotonicity-limited derivatives
        f.uy=ysp*mono_diff(fp[-2*ml].u,fd.u,f.u,fu.u,fp[2*ml].u);
        f.vy=ysp*mono_diff(fp[-2*ml].v,fd.v,f.v,fu.v,fp[2*ml].v);
        f.ux=xsp*mono_diff(fp[-2].u,fl.u,f.u,fr.u,fp[2].u);
        f.vx=xsp*mono_diff(fp[-2].v,fl.v,f.v,fr.v,fp[2].v);

        // Extrapolate to bottom
        f.ud=f.u-0.5*(dy+dt*f.v)*f.uy;
        f.vd=f.v-0.5*(dy+dt*f.v)*f.vy;

        // Extrapolate to left
        f.ul=f.u-0.5*(dx+dt*f.u)*f.ux;
        f.vl=f.v-0.5*(dx+dt*f.u)*f.vx;

        // Extrapolate to right
        f.ur=f.u+0.5*(dx-dt*f.u)*f.ux;
        f.vr=f.v+0.5*(dx-dt*f.u)*f.vx;

        // Extrapolate to top
        f.uu=f.u+0.5*(dy-dt*f.v)*f.uy;
        f.vu=f.v+0.5*(dy-dt*f.v)*f.vy;
    }

    // Select which velocities to use at each boundary, using the Godunov
    // upwinding procedure
    zero_edge_velocities();
#pragma omp parallel for
    for(int j=0;j<n;j++) {
        for(int i=0;i<m;i++) {
            field *fp=fm+(i+ml*j),&f=*fp,&fl=fp[-1],&fr=fp[1],&fd=fp[-ml],&fu=fp[ml];
            if(i<m-1) {
                godunov_set_tang(f.ur,fr.ul,f.vr,fr.vl);
                if(i==0) fl.ur=f.ul;
            } else {fr.ul=f.ur;fr.vl=f.vr;}
            if(j<n-1) {
                godunov_set_tang(f.vu,fu.vd,f.uu,fu.ud);
                if(j==0) fd.vu=f.vd;
            } else {fu.vd=f.vu;fu.ud=f.uu;}
        }
    }

    // Store the tangential stability information in the temporary q arrays
#pragma omp parallel for
    for(j=0;j<n;j++) for(int i=0;i<m;i++) {
        field *fp=fm+(ml*j+i),&f=*fp;
        f.c0=0.5*xsp*(fp[1].ul+f.ul)*(f.ur-fp[-1].ur);
        f.c1=0.5*xsp*(fp[1].ul+f.ul)*(fp[1].vl-f.vl);
        f.c2=0.5*ysp*(fp[ml].vd+f.vd)*(fp[ml].ud-f.ud);
        f.c3=0.5*ysp*(fp[ml].vd+f.vd)*(f.vu-fp[-ml].vu);
    }

    // Now that the tangential-stability derivatives are calculated, compute
    // edge velocities at a half-timestep
    double viscfac=vix==NULL?rhoinv*mu:rhoinv*mu*0.5;
#pragma omp parallel for
    for(j=0;j<n;j++) for(int i=0;i<m;i++) {
        double uxx,uyy,vxx,vyy,px,py,f_x,f_y;

        // Create references to the fields in the neighboring gridpoints
        field *fp=fm+(ml*j+i),&f=*fp;

        // Calculate pressure gradient -- this is cell-centered
        px=0.5*xsp*(fp[ml+1].p+fp[1].p-fp[ml].p-f.p);
        py=0.5*ysp*(fp[ml+1].p-fp[1].p+fp[ml].p-f.p);

        // Calculate force as the sum of pressure gradient and viscous stress
        // term
        centered_diff(fp,uxx,uyy,vxx,vyy);
        f_x=-px*rhoinv+viscfac*(uxx+uyy);
        f_y=-py*rhoinv+viscfac*(vxx+vyy);

        // Extrapolate to bottom
        f.ud=f.u+0.5*(-dy*f.uy+dt*(-f.c0-f.v*f.uy+f_x));
        f.vd=f.v+0.5*(-dy*f.vy+dt*(-f.c1-f.v*f.vy+f_y));

        // Extrapolate to left
        f.ul=f.u+0.5*(-dx*f.ux+dt*(-f.u*f.ux-f.c2+f_x));
        f.vl=f.v+0.5*(-dx*f.vx+dt*(-f.u*f.vx-f.c3+f_y));

        // Extrapolate to right
        f.ur=f.u+0.5*(dx*f.ux+dt*(-f.u*f.ux-f.c2+f_x));
        f.vr=f.v+0.5*(dx*f.vx+dt*(-f.u*f.vx-f.c3+f_y));

        // Extrapolate to top
        f.uu=f.u+0.5*(dy*f.uy+dt*(-f.c0-f.v*f.uy+f_x));
        f.vu=f.v+0.5*(dy*f.vy+dt*(-f.c1-f.v*f.vy+f_y));

        // Store force for reference later
        f.c0=viscfac*(uxx+uyy);
        f.c1=viscfac*(vxx+vyy);
    }

    // Select which velocities to use at each boundary, using the Godunov
    // upwinding procedure
    zero_edge_velocities();
#pragma omp parallel for
    for(int j=0;j<n;j++) {
        for(int i=0;i<m;i++) {
            field *fp=fm+(i+ml*j),&f=*fp,&fr=fp[1],&fu=fp[ml];
            if(i<m-1) godunov_set(f.ur,fr.ul,f.vr,fr.vl);
            else {fr.ul=f.ur;fr.vl=f.vr;}
            if(j<n-1) godunov_set(f.vu,fu.vd,f.uu,fu.ud);
            else {fu.vd=f.vu;fu.ud=f.uu;}
        }
    }

    // Compute source term for the MAC solve
#pragma omp parallel for
    for(j=0;j<n;j++) {
        field *fp=fm+j*ml,*fe=fp+m;
        double *srp=src+j*m;
        while(fp<fe) {
            *(srp++)=xsp*(fp[1].ul-fp->ul)+ysp*(fp[ml].vd-fp->vd);
            fp++;
        }
    }

    // Add contribution from Dirichlet condition at the bottom, and solve for
    // the q field
    const double macfac=0.125*rhoinv*dt;
    for(int i=0;i<m;i++)
        src[i]-=2*yysp*(s_bc[i]=macfac*(gl->p_old[i]+gl->p_old[i+1]+gl->p[i]+gl->p[i+1]));
    ms_mac.solve_v_cycle(time);

    // Update the edge velocities based on the calculation of q
#pragma omp parallel for
    for(j=0;j<n;j++) {
        field *fp=fm+j*ml,*fe=fp+m;
        double *sop=ms_mac.z+j*m;
        if(j!=0) fp->vd-=(*sop-sop[-m])*ysp;
        fp++;sop++;
        while(fp<fe) {
            fp->ul-=(*sop-sop[-1])*xsp;
            if(j!=0) fp->vd-=(*sop-sop[-m])*ysp;
            fp++;sop++;
        }
    }
    for(int i=0;i<m;i++) fm[i].vd-=2*ysp*(ms_mac.z[i]-s_bc[i]);

    // Update the simulation time
    time+=0.5*dt;

#pragma omp parallel for
    for(j=0;j<n;j++) for(int i=0;i<m;i++) {

        // Create references to the fields in the neighboring gridpoints
        field *fp=fm+(ml*j+i),&f=*fp,&fr=fp[1],&fu=fp[ml];

        f.c0=f.u+dt*(-0.5*(xsp*(fr.ul+f.ul)*(fr.ul-f.ul)+ysp*(fu.vd+f.vd)*(fu.ud-f.ud))+f.c0);
        f.c1=f.v+dt*(-0.5*(xsp*(fr.ul+f.ul)*(fr.vl-f.vl)+ysp*(fu.vd+f.vd)*(fu.vd-f.vd))+f.c1);
    }

    // Reset ALL the ghost points according to the boundary conditions. The
    // field variables at the interior grid points remain unchanged. This is
    // needed because the implicit viscous solve can reference all boundary
    // values.
    t0=wtime();
    set_boundaries(dt);
    t_bc+=wtime()-t0;

    // Calculate the implicit viscous term if needed
    if(vix!=NULL) {

        // Compute constants
        const double vfx=viscfac*xxsp*dt,vfy=viscfac*yysp*dt;

        // Compute source term for the linear system solve for the horizontal
        // velocity component
#pragma omp parallel for
        for(j=0;j<n;j++) {
            field *fp=fm+ml*j;
            double *srp=src+j*m;
            for(int i=0;i<m;i++,srp++,fp++) {
                *srp=fp->c0;

                // Add a contribution from the Dirichlet boundary condition at
                // the left, right, and top edges
                if(!x_sym&&i==0) *srp+=vfx*fp[-1].u;
                if(i==m-1) *srp+=vfx*fp[1].u;
                if(j==n-1) *srp+=vfy*fp[ml].u;

                // Add a contribution from the Neumann boundary condition on
                // the bottom edge, using dudy at the current time
                if(j==0) {
                    double dudy=muinv*0.5*(gl->tau_old[i+2]+gl->tau[i+2]);
                    *srp-=vfy*dy*dudy;
                }
            }
        }

        // Solve for the horizontal velocity component
        vix->setup_and_solve(vfx,vfy,time);

#pragma omp parallel for
        for(j=0;j<n;j++) {
            field *fp=fm+ml*j;
            double *srp=src+j*m;
            for(int i=0;i<m;i++,srp++,fp++) {
                *srp=fp->c1;

                // Add a contribution from the Dirichlet boundary condition at
                // the left, right, and top edges
                if(!x_sym&&i==0) *srp+=vfx*fp[-1].v;
                if(i==m-1) *srp+=vfx*fp[1].v;
                if(j==n-1) *srp+=vfy*fp[ml].v;

                // Add a contribution from the Neumann boundary condition on
                // the bottom edge, using dudx+dvdy=0 at the current time
                if(j==0) {
                    double dvdy=-0.5*xsp*(fp[1].u-fp[-1].u);
                    *srp-=vfy*dy*dvdy;
                }
            }
        }

        // Solve for the vertical velocity component
        viy->setup_and_solve(vfx,vfy,time);

        // Copy the two solutions back into the main data structure
#pragma omp parallel for
        for(j=0;j<n;j++) for(int i=0;i<m;i++) {
            fm[ml*j+i].c0=vix->z[m*j+i];
            fm[ml*j+i].c1=viy->z[m*j+i];
        }
    }

    // Calculate the source term for the finite-element projection
#pragma omp parallel for
    for(j=1;j<ne;j++) {
        field *fp=fm+j*ml,*fs=fp,*fpen=fp+m,*fe=fpen+1;
        double *srp=src+j*me,spx,spy;
        while(fp<fe) {
            spx=spy=0;

            // Bottom left
            if(fp!=fs) {spx+=fp[-ml-1].c0;spy+=fp[-ml-1].c1;}

            // Top left
            if(fp!=fs&&j<n) {spx+=fp[-1].c0;spy-=fp[-1].c1;}

            // Bottom right
            if(fp!=fpen) {spx-=fp[-ml].c0;spy+=fp[-ml].c1;}

            // Top right
            if(fp!=fpen&&j<n) {spx-=fp->c0;spy-=fp->c1;}

            *(srp++)=0.5*(dy*spx+dx*spy)*rho/dt;
            fp++;
        }
    }

    // Add flux terms to FEM projection
    if(!x_sym) for(j=1;j<ne;j++) {
        if(j<n) src[j*me]+=0.25*dy*rho/dt*(fm[j*ml].c0+fm[j*ml-1].u);
        src[j*me]+=0.25*dy*rho/dt*(fm[(j-1)*ml].c0+fm[(j-1)*ml-1].u);
    }
    for(j=1;j<ne;j++) {
        if(j<n) src[m+j*me]-=0.25*dy*rho/dt*(fm[j*ml+m-1].c0+fm[j*ml+m].u);
        src[m+j*me]-=0.25*dy*rho/dt*(fm[(j-1)*ml+m-1].c0+fm[(j-1)*ml+m].u);
    }
    for(int i=0;i<me;i++) {
        if(i>0) src[n*me+i]-=0.25*dx*rho/dt*(fm[(n-1)*ml+i-1].c1+fm[n*ml+i-1].v);
        if(i<m) src[n*me+i]-=0.25*dx*rho/dt*(fm[(n-1)*ml+i].c1+fm[n*ml+i].v);
    }

    // Fill in pressure values on bottom boundary for the source term calculation
    for(int i=0;i<me;i++) src[i]=ms_fem.a_cc(i,i)*gl->p[i];

    // Solve the finite-element problem
    ms_fem.solve_v_cycle(time);

    // Update u and v based on c0, c1, and the computed pressure
#pragma omp parallel for
    for(j=0;j<n;j++) {
        field *fp=fm+j*ml,*fe=fp+m;
        double *sop=ms_fem.z+j*me;
        while(fp<fe) {
            fp->u=fp->c0-0.5*dt*rhoinv*xsp*(sop[me+1]+sop[1]-sop[me]-*sop);
            fp->v=fp->c1-0.5*dt*rhoinv*ysp*(sop[me+1]-sop[1]+sop[me]-*sop);
            (fp++)->p=*(sop++);
        }
        fp->p=*sop;
    }
    field *fp=fm+n*ml,*fe=fp+me;
    double *sop=ms_fem.z+n*me;
    while(fp<fe) (fp++)->p=*(sop++);

    // Since the velocity field changed, update the viscous boundaries
    t0=wtime();
    set_viscous_boundaries();
    t_bc+=wtime()-t0;

    // If the non-inertial frame is in use, then modify the velocity here
    if(nif_mode) adjust_frame();
}

/** Calculates the second-order derivatives of the velocity components using
 * centered differences.
 * \param[in] fp a pointer to the gridpoint to consider.
 * \param[out] (uxx,uyy) the derivatives of horizontal velocity.
 * \param[out] (vxx,vyy) the derivatives of vertical velocity. */
void fluid_2d::centered_diff(field *fp,double &uxx,double &uyy,double &vxx,double &vyy) {
    uyy=yysp*(fp[-ml].u-2*fp->u+fp[ml].u);
    vyy=yysp*(fp[-ml].v-2*fp->v+fp[ml].v);
    uxx=xxsp*(fp[-1].u-2*fp->u+fp[1].u);
    vxx=xxsp*(fp[-1].v-2*fp->v+fp[1].v);
}

/** Calculates the fourth-order monotonicity-limited derivative of a function.
 * \param[in] (f0,f1,f2,f3,f4) five values of the function at consecutive gridpoints.
 * \return The derivative, without normalization due to the grid spacing. */
double fluid_2d::mono_diff(double &f0,double &f1,double &f2,double &f3,double &f4) {
    double tdrc=f3-f1,s=(1/6.)*(4*tdrc-delta_f(f0,f1,f2)-delta_f(f2,f3,f4)),
           t=delta_lim(f1,f2,f3);
    return tdrc>0?min(s,t):max(s,-t);
}

/** Calculates the delta_f field at the field fp, which is an intermediate
 * quantity in the fourth-order monotonicity-limited derivative.
 * \param[in] f0 the value at left/down.
 * \param[in] f1 the value at center.
 * \param[in] f2 the value at right/upper.
 * \return The delta_f value. */
inline double fluid_2d::delta_f(double &f0,double &f1,double &f2) {
    double drc=0.5*(f2-f0);
    return drc>0?min(drc,delta_lim(f0,f1,f2)):max(drc,-delta_lim(f0,f1,f2));
}

/** Calculates the delta_lim field at the field fp, which is an intermediate
 * quantity in the fourth-order monotonicity-limited derivative.
 * \param[in] f0 the value at left/down.
 * \param[in] f1 the value at center.
 * \param[in] f2 the value at right/upper.
 * \return The delta_lim value. */
inline double fluid_2d::delta_lim(double &f0,double &f1,double &f2) {
    double drm=f1-f0,drp=f2-f1;
    return drm>0?(drp>0?2*min(drm,drp):0)
            :(drp>0?0:-2*max(drm,drp));
}

/** Selects which velocities to use at an edge based on the Godunov upwinding
 * criterion.
 * \param[in] (u0,u1) the velocity components transverse to the edge.
 * \param[in] (v0,v1) the velocity components parallel to the edge. */
inline void fluid_2d::godunov_set(double &u0,double &u1,double &v0,double &v1) {
    if(u0>-u1) {
        if(u0>0) {u1=u0;v1=v0;return;}
    } else if(u1<0) return;
    u1=0;
    v1=0.5*(v0+v1);
}

/** Selects which velocities to use at an edge based on the Godunov upwinding
 * criterion. This also stores the additional velocity term that is needed to
 * compute the tangential stability term, which has a slightly different
 * formula.
 * \param[in] (u0,u1) the velocity components transverse to the edge.
 * \param[in] (v0,v1) the velocity components parallel to the edge. */
inline void fluid_2d::godunov_set_tang(double &u0,double &u1,double &v0,double &v1) {
    if(u0>-u1) {
        if(u0>0) {u1=u0;v1=v0;return;}
    } else if(u1<0) {u0=u1;return;}
    v1=0.5*(v0+v1);
    u0=0.5*(u0+u1);
    u1=0;
}

/** Sets the extrapolated edge velocities to zero. For the splashing problem,
 * it is only needed at the x=0 boundary for the symmetrized case. */
void fluid_2d::zero_edge_velocities() {
    if(x_sym) for(field *fp=fm,*fe=fp+n*ml;fp<fe;fp+=ml) fp->ul=0;
}

/** Sets the fields in the ghost regions according to the boundary conditions. */
void fluid_2d::set_boundaries(double dt) {

    // Set right ghost values (and left ghost values if x_sym is false)
    set_inviscid_boundaries(dt);

    // Set the bottom ghost values (and left ghost values if x_sym is true)
    set_viscous_boundaries();
}

/** Sets up the fluid tracers by initializing them at random positions. */
void fluid_2d::init_tracers() {
    for(double *tp=tm;tp<tm+(ntrace<<1);) {

        // Create a random position vector within the simulation region
        *(tp++)=ax+(bx-ax)/RAND_MAX*static_cast<double>(rand());
        *(tp++)=ay+(by-ay)/RAND_MAX*static_cast<double>(rand());
    }
}

/** Moves the tracers according to the bilinear interpolation of the fluid
 * velocity.
 * \param[in] dt the timestep to use. */
void fluid_2d::update_tracers(double dt) {
    int i,j;
    double x,y,*tp=tm,*te=tm+(ntrace<<1);
    field *fp;
    while(tp<te) {

        // Find which grid cell the tracer is in
        x=(*tp-ax)*xsp-0.5;y=(tp[1]-ay)*ysp-0.5;
        i=int(x);if(i<-1) i=-1;if(i>m-1) i=m-1;
        j=int(y);if(j<-1) j=-1;if(j>n-1) j=n-1;

        // Compute the tracer's fractional position with the grid cell
        x-=i;y-=j;

        // Update tracer's position
        fp=fm+(i+ml*j);
        *(tp++)+=dt*((1-y)*(fp->u*(1-x)+fp[1].u*x)+y*(fp[ml].u*(1-x)+fp[ml+1].u*x));
        *(tp++)+=dt*((1-y)*(fp->v*(1-x)+fp[1].v*x)+y*(fp[ml].v*(1-x)+fp[ml+1].v*x));
    }
}
