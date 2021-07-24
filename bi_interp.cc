#include "bi_interp.hh"
#include <cstdio>

/** Constructs the bicubic interpolation class, initializing geometry
 * constants, the workspace, and links to the fields to interpolate
 * \param[in] (m_,n_) the number of grid points in the horizontal and
 *                    vertical directions of the fields to interpolate.
 * \param[in] (ax_,bx_) the lower and upper x-coordinate field bounds.
 * \param[in] (ay_,by_) the lower and upper y-coordinate field bounds. */
bicubic_interp::bicubic_interp(int m_,int n_,double ax_,double bx_,double ay_,double by_)
    : m(m_), n(n_), ax(ax_), ay(ay_), xsp((m-1)/(bx_-ax)), ysp((n-1)/(by_-ay)),
    u(new double[m*n]), ijc(-1) {
    *a=a[1]=a[2]=a[3]=0;
    a[4]=a[5]=a[6]=a[7]=0;
    a[8]=a[9]=a[10]=a[11]=0;
    a[12]=a[13]=a[14]=a[15]=0;
}

/** The class destructor frees the dynamically allocated memory. */
bicubic_interp::~bicubic_interp() {
    delete [] u;
}

/** Calculates the grid square for a given position at which to interpolate, and
 * sets up the workspace.
 * \param[in,out] (x,y) the position to use, which is mapped into a fractional
 *                      position in the grid square. */
void bicubic_interp::grid_index(double &x,double &y) {

    // Compute which grid square the given point is within
    x=(x-ax)*xsp;y=(y-ay)*ysp;
    int i=static_cast<int>(x),j=static_cast<int>(y),ij;

    // If the grid square is out of range, then map it to the nearest one
    // available
    if(i<0) i=0;else if(i>m-2) i=m-2;
    if(j<0) j=0;else if(j>n-2) j=n-2;

    // Compute index of the grid square to interpolate in. If this different
    // from the previous call then set up the temporary workspace with values
    // needed in the interpolation.
    ij=i+m*j;
    if(ijc!=ij) table_setup(i,j,ij);
    x-=i;y-=j;
}

/** Sets up the table of coefficients of the bicubic interpolation function in
 * a grid square.
 * \param[in] (i,j) the (x,y)-indices of the lower left of the grid square.
 * \param[in] ij the index of the grid square. */
void bicubic_interp::table_setup(int i,int j,int ij) {
    ijc=ij;
    double *up=u+ij;
    double c00,c01,c02,c03;
    double c10,c11,c12,c13;
    double c20,c21,c22,c23;
    double c30,c31,c32,c33;

    // Compute the contributions of the bottom and top rows of the grid square
    compute_x(i,up,c01,c11,c21,c31);
    compute_x(i,up+m,c02,c12,c22,c32);
    if(j==0) {

        // If this is the zeroth row, then compute the contribution from the
        // row above grid square (but not below). Use all the contributions to
        // make the table based on quadratic interpolation from the three rows.
        compute_x(i,up+2*m,c03,c13,c23,c33);
        fill_ad(a,c01,c02,c03);
        fill_ad(a+4,c11,c12,c13);
        fill_ad(a+8,c21,c22,c23);
        fill_ad(a+12,c31,c32,c33);
    } else if(j==n-2) {

        // If this is the topmost row, then compute the contribution from the
        // row below the grid square (but not above). Use all the contributions
        // to make the table based on quadratic interpolation from the three
        // rows.
        compute_x(i,up-m,c00,c10,c20,c30);
        fill_au(a,c00,c01,c02);
        fill_au(a+4,c10,c11,c12);
        fill_au(a+8,c20,c21,c22);
        fill_au(a+12,c30,c31,c32);
    } else {

        // If this row is in the bulk of the grid, then compute contributions
        // from the rows below and above the grid square. Use all the
        // contributions to make the table based on cubic interpolation from
        // the four rows.
        compute_x(i,up-m,c00,c10,c20,c30);
        compute_x(i,up+2*m,c03,c13,c23,c33);
        fill_a(a,c00,c01,c02,c03);
        fill_a(a+4,c10,c11,c12,c13);
        fill_a(a+8,c20,c21,c22,c23);
        fill_a(a+12,c30,c31,c32,c33);
    }
}

/** Compute the interpolation contributions from a row.
 * \param[in] i the horizontal index of the grid points to consider.
 * \param[in] ij the grid index to consider.
 * \param[out] (c0,c1,c2,c3) the coefficients of 1, x, x^2, x^3 in the
 *                           contribution. */
void bicubic_interp::compute_x(int i,double *up,double &c0,double &c1,double &c2,double &c3) {
    c0=*up;
    if(i==0) {

        // If this is the zeroth column then do quadratic interpolation
        // using three field values
        c1=-1.5*(*up)+2*up[1]-0.5*up[2];
        c2=0.5*(*up)-up[1]+0.5*up[2];
        c3=0;
    } else if(i==m-2) {

        // If this is the last column, then do quadratic interpolation
        // using three field values
        c1=-0.5*up[-1]+0.5*up[1];
        c2=0.5*up[-1]-*up+0.5*up[1];
        c3=0;
    } else {

        // If this is a central column then do cubic interpolation
        c1=-0.5*up[-1]+0.5*up[1];
        c2=up[-1]-2.5*(*up)+2*up[1]-0.5*up[2];
        c3=-0.5*up[-1]+1.5*(*up)-1.5*up[1]+0.5*up[2];
    }
}

/** Calculates a bicubic interpolation at a point.
 * \param[in] (x,y) the point to consider.
 * \return The interpolated value. */
double bicubic_interp::f(double x,double y) {
    grid_index(x,y);
    return yl(a,y)+x*(yl(a+4,y)+x*(yl(a+8,y)+x*yl(a+12,y)));
}

/** Calculates a bicubic interpolation at a point, and its gradient.
 * \param[in] (x,y) the point to consider.
 * \param[out] (fx,fy) the components of the gradient.
 * \return The interpolated value. */
double bicubic_interp::f_grad_f(double x,double y,double &fx,double &fy) {
    grid_index(x,y);
    fx=xsp*(yl(a+4,y)+x*(2*yl(a+8,y)+3*x*yl(a+12,y)));
    fy=ysp*(dyl(a,y)+x*(dyl(a+4,y)+x*(dyl(a+8,y)+x*dyl(a+12,y))));
    return yl(a,y)+x*(yl(a+4,y)+x*(yl(a+8,y)+x*yl(a+12,y)));
}

/** Fills four table entries based on quadratic interpolation of three values
 * (assuming the lower, fourth one is out of range).
 * \param[in] ap a pointer to the table entries to fill.
 * \param[in] (c1,c2,c3) the three values to use. */
void bicubic_interp::fill_ad(double *ap,double c1,double c2,double c3) {
    *ap=c1;
    ap[1]=-1.5*c1+2*c2-0.5*c3;
    ap[2]=0.5*c1-c2+0.5*c3;
    ap[3]=0;
}

/** Fills four table entries based on quadratic interpolation of three values
 * (assuming the upper, fourth one is out of range).
 * \param[in] ap a pointer to the table entries to fill.
 * \param[in] (c0,c1,c2) the three values to use. */
void bicubic_interp::fill_au(double *ap,double c0,double c1,double c2) {
    *ap=c1;
    ap[1]=-0.5*c0+0.5*c2;
    ap[2]=0.5*c0-c1+0.5*c2;
    ap[3]=0;
}

/** Fills four table entries based on cubic interpolation of four values.
 * \param[in] ap a pointer to the table entries to fill.
 * \param[in] (c0,c1,c2,c3) the four values to use. */
void bicubic_interp::fill_a(double *ap,double c0,double c1,double c2,double c3) {
    *ap=c1;
    ap[1]=-0.5*c0+0.5*c2;
    ap[2]=c0-2.5*c1+2*c2-0.5*c3;
    ap[3]=-0.5*c0+1.5*c1-1.5*c2+0.5*c3;
}
