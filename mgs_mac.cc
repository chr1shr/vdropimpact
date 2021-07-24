#include "mgs_mac.hh"
#include "fluid_2d.hh"
#include "tgmg.hh"

/** The constructor sets many internal constants from the parent fluid_2d
 * class.
 * \param[in] f a reference to the parent fluid_2d class. */
mgs_mac::mgs_mac(fluid_2d &f) : mgs_common(f.m,f.n), xxsp(f.xxsp), yysp(f.yysp),
    cc0(-1./(2*xxsp+3*yysp)), cc1(-1./(2*xxsp+2*yysp)), cc2(-1./(2*xxsp+yysp)),
    cc3(-1./(xxsp+3*yysp)), cc4(-1./(xxsp+2*yysp)), cc5(-1./(xxsp+yysp)),
    acc(tgmg_accuracy(2*(xxsp+yysp)*f.dt_reg,1e6)),
    mg(*this,f.src,z) {
    mg.setup();
    mg.clear_z();
}

/** An alternative constructor that independently sets up the multigrid
 * class for testing purposes.
 * \param[in] (m_,n_) the dimensions of the grid.
 * \param[in] (dx,dy) the grid spacings in the x and y directions,
 *            respectively.
 * \param[in] b a pointer to the source array. */
mgs_mac::mgs_mac(int m_,int n_,double dx,double dy,double *b) :
    mgs_common(m_,n_), xxsp(1./(dx*dx)), yysp(1./(dy*dy)),
    cc0(-1./(2*xxsp+3*yysp)), cc1(-1./(2*xxsp+2*yysp)), cc2(-1./(2*xxsp+yysp)),
    cc3(-1./(xxsp+3*yysp)), cc4(-1./(xxsp+2*yysp)), cc5(-1./(xxsp+yysp)),
    acc(tgmg_accuracy(2*(xxsp+yysp),1e6)), mg(*this,b,z) {}

// Explicit instantiation
#include "tgmg.cc"
template class tgmg<mgs_mac,double,double>;
template void tgmg_base<mgs_mac,double,double>::output(char const*,double*,double,double,double,double);
template void tgmg_base<mgs_mac,double,double>::clear_z();
template void tgmg_base<mgs_mac,double,double>::output_res(char const*,double,double,double,double);
