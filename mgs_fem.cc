#include "mgs_fem.hh"
#include "fluid_2d.hh"
#include "tgmg.hh"

/** The constructor sets many internal constants from the parent fluid_2d
 * class.
 * \param[in] f a reference to the parent fluid_2d class. */
mgs_fem::mgs_fem(fluid_2d &f) : mgs_common(f.me,f.ne), dydx(f.dy/f.dx),
    dxdy(f.dx/f.dy), fm(4./3.*(dxdy+dydx)), fm_inv(1.0/fm),
    fey(1./3.*(-2*dxdy+dydx)), hey(0.5*fey), fex(1./3.*(-2*dydx+dxdy)),
    hex(0.5*fex), fc(-1./6.*(dxdy+dydx)), acc(tgmg_accuracy(fm,1e8)),
    mg(*this,f.src,z) {
    mg.setup();
    mg.clear_z();
}

/** An alternative constructor that independently sets up the multigrid
 * class for testing purposes.
 * \param[in] (m_,n_) the dimensions of the grid.
 * \param[in] (dx,dy) the grid spacings in the x and y directions,
 *            respectively.
 * \param[in] b a pointer to the source term array. */
mgs_fem::mgs_fem(int m_,int n_,double dx,double dy,double *b) :
    mgs_common(m_,n_), dydx(dy/dx), dxdy(dx/dy), fm(4./3.*(dxdy+dydx)),
    fm_inv(1.0/fm), fey(1./3.*(-2*dxdy+dydx)), hey(0.5*fey),
    fex(1./3.*(-2*dydx+dxdy)), hex(0.5*fex), fc(-1./6.*(dxdy+dydx)),
    acc(tgmg_accuracy(fm,1e8)), mg(*this,b,z) {}

/** Performs multiplication by the linear system, excluding the diagonal
 * term.
 * \param[in] i the horizontal grid index.
 * \param[in] ij the overall grid index.
 * \return The result of the multiplication. */
double mgs_fem::mul_a(int i,int ij) {

    // Return zero if on bottom row, to apply Dirichlet pressure boundary
    // condition
    if(ij<m) return 0;

    // Deal with the top boundary
    double *w=z+ij,*wm=w-m;
    return ij>=mn-m?(i==0?fc*wm[1]+hey*(*wm)+hex*w[1]:
                     i==m-1?fc*wm[-1]+hey*(*wm)+hex*w[-1]:
                     fc*(wm[-1]+wm[1])+fey*(*wm)+hex*(w[-1]+w[1])):

    // Deal with the middle rows
                    (i==0?fc*(wm[1]+w[m+1])+hey*(*wm+w[m])+fex*w[1]:
                    i==m-1?fc*(wm[-1]+w[m-1])+hey*(*wm+w[m])+fex*w[-1]:
                    fc*(wm[-1]+wm[1]+w[m-1]+w[m+1])+fey*(*wm+w[m])+fex*(w[-1]+w[1]));
}

// Explicit instantiation
#include "tgmg.cc"
template class tgmg<mgs_fem,double,double>;
template void tgmg_base<mgs_fem,double,double>::output(char const*,double*,double,double,double,double);
template void tgmg_base<mgs_fem,double,double>::output_res(char const*,double,double,double,double);
template void tgmg_base<mgs_fem,double,double>::clear_z();
