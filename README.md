# VDropImpact: A simulation of viscous effects in early-time drop impact dynamics
This repository contains a C++ code for simulating a reduced model of the
initial dynamics of when a drop of liquid impacts a horizontal surface. In
particular it can be used to explore the effects of liquid viscosity on this
process. It supports the following scientific publication:

- Shruti Mishra, Shmuel M. Rubinstein, and Chris H. Rycroft, *Computing
  the viscous effect in early drop impact dynamics*, J. Fluid Mech **945**,
  A13 (2022). [doi:10.1017/jfm.2022.445](https://doi.org/10.1017/jfm.2022.445).
  [arXiv:2108.11497](https://arxiv.org/abs/2108.11497).

## Background
Understanding how liquid drops impact and splash on contact with solid surfaces
has been the subject of scientific research for over a century. Experiments
performed by Worthington in 1877 demonstrated the wide variety of patterns
formed by drops falling on a plate, even before the invention a flash
photography [1]. Since then, drop impact has been the subject of a wide range
of theoretical, computational, and experimental studies [2].

In this project we studied the early-time dynamics of a liquid drop as it first
approaches the solid surface. The drop decelerates and deforms, creating a
*dimple* of trapped gas. The drop then spreads outward, forming a thin gas
layer between the surface and the liquid, with thickness on order of 200 nm.
The early-time dynamics and formation of this layer have been analyzed by
coupling a one-dimensional model for the gas layer height to a two-dimensional
liquid model [3,4]. Since liquid viscosity was assumed not to be important,
on the scales of interest, [potential flow
theory](https://en.wikipedia.org/wiki/Potential_flow) was used to simulate the
liquid domain. This considerably simplifies the analysis since entire
two-dimensional flow can be determined via fields described only on the
one-dimensional interface.

Experiments by Kolinski *et al.* [5] showed that the inviscid assumption used in
this previous work was inaccurate. As the thin gas layer is formed, it tilts
upward, lifting off the surface. The time at which lift-off happens is
dependent on the liquid viscosity. For liquid viscosities from 1 cSt to
100 cSt, Kolinski *et al.* observed an approximate square root dependence of
lift-off time on viscosity.

In this study, we developed a reduced model to examine these viscous effects in
early-time drop impact dynamics. We coupled a one-dimensional gas layer model
to a two-dimensional simulation of the liquid using the incompressible
[Navier–Stokes equations](https://en.wikipedia.org/wiki/Navier–Stokes_equations).
Full implementation details of the Navier–Stokes solver are given by Yu *et
al.* [6] and Rycroft *et al.* [7]. This is more computationally expensive than
the previous work [3,4] since the full two-dimensional flow must be resolved.
However, with this simulation we can successfully recreate and analyze the
viscosity dependence seen by Kolinski *et al.*, and explore the role of surface
tension, impact velocity, and drop geometry.

## Compiling the code
The code is written in C++ and uses the OpenMP library for multithreading. It has
been tested on Linux, MacOS, and Windows (via [Cygwin](https://www.cygwin.com)).

- The code requires on TGMG, a C++ library for solving linear systems using the
  geometric multigrid method [8,9], which is available as a
  [separate repository on GitHub](http:/github.com/chr1shr/tgmg).

- The code outputs data in a binary format that can be read by the freeware
  plotting program [Gnuplot](http://www.gnuplot.info). The code uses
  utils-gp, a collection of tools for processing and analyzing Gnuplot output
  files. This is available as a [separate repository on
  GitHub](http://github.com/chr1shr/utils-gp).

- The utils-gp repository requires [libpng](http://www.libpng.org/pub/png/) for
  making for full functionality, but this dependency can be omitted.

- The code requires [LAPACK](http://www.netlib.org/lapack/) for solving a
  tridiagonal linear system arising in the gas layer description. LAPACK is
  often installed by default on many new operating systems, and is available
  via many software package management systems.

By default the code assumes that the **vdropimpact**, **tgmg**, **utils-gp**
repositories are placed in the same parent directory.

To compile the code it is necessary to create a common configuration file
called **config.mk** in the parent directory, which can be used by all three
repositories. Several templates are provided in the **config** directory. To
use, copy one of the templates into the parent directory. From the
**vdropimpact** directory, on a Linux computer, type
```Shell
cp config/config.mk.linux ../config.mk
```
On a Mac using GCC 11 installed via [MacPorts](http://www.macports.org), type
```Shell
cp config/config.mk.mac_mp ../config.mk
```
On a Mac using GCC installed via [Homebrew](http://brew.sh), type
```Shell
cp config/config.mk.mac_hb ../config.mk
```
On a Windows computer with Cygwin installed, type
```Shell
cp config/config.mk.win_cw ../config.mk
```
After this, the code can be compiled by typing
```Shell
make
```
This will build several executables such as **fluid_test** and **h_analysis**.

## Configuring and running a simulation
The main program for running a drop impact simulation is called **fluid_test**.
To run this program, it must be supplied with a text file containing the
simulation configuration. Several samples are provided in the **sims**
directory. A good place to start is **id_vis32.cfg**, which simulates the
initial deceleration of a drop as it approaches a surface—see Sec. 4.1 of the
paper. This simulation can be run by typing
```Shell
OMP_NUM_THREADS=<n> ./fluid_test sims/id_vis32.cfg
```
where `<n>` is replaced with the number of OpenMP threads to use.

The configuration file is divided into sections. A pound symbol in front of any
line marks it as a comment. The first section controls geometrical constants:
```
# Geometry
grid_points         2048 256
L_nd                30
h0_nd               15
t_end_nd            50
```
`grid_points` sets the number of grid points in the horizontal and vertical
directions in the liquid domain. The code uses the same grid spacing in the *x*
and *y* directions, and therefore the grid size specified here affects the
aspect ratio of the simulation.

The next three parameters are non-dimensionalized and set the size of
the simulation, the initial height of the drop above the surface, and the
simulation duration. Formulae for linking them to physical units are provided
in Sec. 2.4 of the paper.

The next section of the input file contains parameters relating to the
computation:
```
# Miscellaneous constants
tracers             512
frames              250
tmult               8e-3
#restart_freq        50
```
To visualize the flow in the liquid domain, the code randomly initializes
tracers throughout the domain that move with the liquid. The `tracers` option
sets how many tracers to use. Tracers have no effect on the simulation itself.

`frames` controls how many snapshots of the simulation fields to output over
the simulation duration. `tmult` controls the size of the timestep. When the
viscous term in the Navier–Stokes equation is handled implicitly, the timestep
is chosen as proportional to the grid spacing &Delta;*x*. Otherwise the timestep
is proportional to &Delta;*x*&sup2;.

The code can output restart files that save a snapshot of all simulation fields
at full precision. The code can be resumed from a restart file by passing the
`-r` flag option to *fluid_test*. The `restart_freq` option specifies how
often to save restart files.

The next section sets the physical parameters for the drop, and are specified
in physical units:
```
# Liquid/drop parameters:
# nul - liquid kin. viscosity (m^2/s)
# nul_cSt - liquid kin. viscosity (centistokes)
# rhol - liquid density (kg m^{-3})
# R - drop radius (m)
# V - initial drop velocity (m/s)
nul_cSt             10
rhol                997.96
R                   1.5e-3
V                   0.45
```
Two different options `nul` and `nul_cSt` can be used to specify the liquid
kinematic viscosity, using two different units. The next section sets the
physical parameters for the gas layer and the liquid--gas interface in physical 
units:
```
# Gas parameters:
# gamma - exponent in equation of state
# alpha - 1/exponent in equation of state
# sigma - surface tension (N/m)
# mug - gas dyn. viscosity (Pa s)
# Pamb - ambient gas pressure (Pa)
gamma               1.4
sigma               0
mug                 1.820775e-5
Pamb                1e5
```
The next section controls different boolean simulation options:
```
# Extra flags
x_sym
implicit_visc
gas_layer_model
```
As described in Sec. 2.2 of the paper the simulation domain is symmetric
about *x*=0. The `x_sym` option halves the simulation domain from [-*L*,*L*]
to [0,*L*] and uses symmetry boundary conditions at *x*=0.

The `implicit_visc` option handles the viscous term in the Navier–Stokes
equations using an implicit Crank–Nicolson discretization. This removes a
timestep restriction that scales quadratically with the grid spacing. The
`gas_layer_model` calculates the fields in the gas layer using the PDE model
given in Sec. 2.2 of the paper, following the numerical methods described in
Sec. 3.3. There is also another option called `gas_layer_data` used for testing
purposes, where the boundary is updated according to a pre-computed table of
data.

An additional option called `mr_time_output` can be specified that makes the
command-line utility output timing statistics in a machine-readable format.
This is useful for analyzing code performance (see the **t_perf.pl** script).

The options `x_sym`, `implicit_visc`, and `gas_layer_model` are used for all
results in the paper.

The next section enables the usage of a non-inertial frame. 
```
# Non-inertial frame:
# nif_center - use droplet velocity at x=0
# nif_range <w> - use average droplet velocity in |x|<w (w is given in m)
#nif_center
#nif_range           1e-4
```
By default, the fluid simulation is performed in an inertial frame moving
downward with the initial drop velocity. This means that at *t*=0 in the
simulation, the fluid simulation velocity is zero, and the fluid velocity
points upward once the drop begins to decelerate. Hence, fluid enters the
domain from the bottom boundary.

The code has the capability to use a non-inertial frame, where the frame is
slowed to take into account the deceleration of the drop, which limits
the fluid entering the bottom boundary. This creates a fictious vertical
acceleration in the fluid domain. Several procedures are available for
setting the velocity of the frame: `nif_center` tracks the velocity at *x*=0,
and `nif_range` uses the average velocity over a range |*x*|&le;w.

In the paper, the non-inertial frame is not used since it has a minimal effect
on the results.

The final section of the code controls the output of the simulation:
```
# Types of file output:
# (u,v) - velocity components
# p - pressure
# h - height
# w - vorticity
# fbd - flow across boundary
output              u v p h
```
The simulation creates an output directory, where the **cfg** suffix is
replaced with **odr**. Each of the requested fields is saved to the directory
with a numerical suffix to indicate the frame number.

The data is output in Gnuplot's binary matrix format. Full details on this
format can be found by typing "help binary matrix nonuniform" within Gnuplot;
[see here for additional explanation](https://github.com/chr1shr/utils-gp#the-binary-matrix-format).
For an *m*&times;*n* grid, this format consists of an (*m*+1)&times;(*n*+1)
array of single-precision floating point numbers, where the extra row and
column contain the coordinate information. Since the simulation internally
uses double precision, this output format truncates the precision by
half. But this halves the required disk storage, and single precision is
usually sufficient for plotting and analysis. The one-dimensional fields for
the `h` and `fbd` options are output as *m*&times;1 fields in the Gnuplot
format.

For the example configuration file, the height profiles can be plot using
the following Gnuplot commands. To begin, set the axis labels and switch off
the key:
```Gnuplot
unset key
set xlabel 'x (micrometer)'
set ylabel 'h (micrometer)'
```
The height profiles can be plot using the for-loop functionality in Gnuplot
with the following command:
```Gnuplot
a=1e6
plot for [k=0:250:20] sprintf('id_vis32.odr/height.%d',k) u ($1*a):($3*a) matrix binary
```
Here, the variable `a=1e6` is used to convert the values in the file from
meters into micrometers. This will produce a plot shown below:

![Plots of the height profile over time](https://people.seas.harvard.edu/~chr/software/vdi_graph1.png)

This shows the overall behavior of the simulation, where initially parabolic
height profile deforms as it approaches the surface. To see the thin gas
layer in more detail, the axes can be changed:
```Gnuplot
set xrange [150:550]
set yrange [0:1]
plot for [k=0:250:10] sprintf('id_vis32.odr/height.%d',k) u ($1*a):($3*a) matrix binary
```
This produces the plot shown below:

![Zoomed-in plots of the height profile over time](https://people.seas.harvard.edu/~chr/software/vdi_graph2.png)

The lift-off behavior can be observed at approximately *x*=330 &mu;m. Even
though lift-off can be clearly seen, the grid resolution specified in
**id_vis32.cfg** is relatively coarse and in the paper it is only used to study
the initial dynamics in Sec. 4.1. See the other examples **lo_vis10.cfg**,
**ld_vis32.cfg**, and **ld_vis100.cfg** for configuring a higher-resolution
simulation for better accuracy in computing the lift-off behavior.

## Analysis codes and utilities
In addition to the main **fluid_test** program, a variety of utility functions
and scripts are provided that were used during the project.

- **mg_test** – This program tests the multigrid method applied to the MAC and
  finite-element linear systems that are part of the incompressible
  Navier–Stokes solver in the liquid. The code can record the time to perform
  multigrid V-cycles, and can output solution fields in the Gnuplot binary
  matrix format.

- **edges** – The velocity fields that are outputted by the simulation are
  cell-centered. Hence at the edges of the simulation, the grid points are
  inset by half a grid spacing from the physical domain. This utility
  interpolates the velocities to the physical domain, which can be useful
  in determining how well boundary conditions are being imposed. This is
  subject of research in projection methods for incompressible Navier–Stokes
  equations [10,11].

- **tri_test** - This code tests the tridiagonal matrix solver on a simple
  linear system describing a finite-difference discretization of a
  one-dimensional Poisson problem.

- **h_analysis** – This utility can be run on completed simulation output
  directory (run with the `x_sym` option) to track droplet height and curvature
  near *x*=0. It creates a file with a **han** suffix. This can subsequently be
  passed to the script **min_hei.pl** to extract the height *H*<sup>&ast;</sup>
  [3,4] where the drop reaches a stagnation point.

- **tt_analysis** – This utility can be run on a completed simulation output
  directory to compute the trajectory of the leading tip via the global minimum
  of the height profile. It outputs a file with a **tta** suffix containing the
  trajectory of the leading tip. It also prints information about the lift-off
  time. By default the utility applies Gaussian smoothing to the height profile.
  If the program is run with the `-r` flag then Gaussian smoothing is switched
  off, and the resulting trajectory is saved to a file with a **ttb** suffix.

- **p_strip** – This utility processes all of the *m*&times;*n* pressure
  fields, **p**.&lt;n&gt; in the output directory, extracts the pressure field
  at the base at *y*=0, and saves the resulting information as *m*&times;1
  pressure field.

- **collate** - This utility takes the one-dimensional height fields or gas
  pressure fields, and collates them into a two-dimensional Gnuplot binary
  matrix format in the (*x*,*t*) space.

- **t_perf.pl** – This script can be run on the terminal output of a simulation
  that was run with the `mr_time_output` option that created machine-readable
  timing information. The script prints statistics about the time taken in
  different parts of the code, similar to Table 4 in Appendix B.

- **sims/create_cfg.pl** – In the paper, many sweeps of simulations were
  performed using ranges of different viscosities. This script can be used
  systematically generate **cfg** files for a sweep of simulations.
  The script makes use of the **template.cfg** file, performs
  search-and-replace on keywords, and outputs a set of **cfg** files.

## Code structure
The code is structured around several C++ classes:

- **fluid_2d** – This is the main class containing the core routines for the
  incompressible Navier–Stokes solver. The core numerical routines are defined
  in the file **fluid_2d.cc**. Additional routines for doing input/output and
  routines that are specific to the gas layer coupling, are defined in the file
  **fluid_2d_sio.cc**.

- **gas_layer** – This class handles time evolution of the gas layer, and computes
  the pressure and tangential stress that are required for coupling to the liquid.
  It is a pure virtual class that has two different options.
  **gas_layer_model** updates the gas layer according to the partial
  differential equation in Sec. 2.2 of the paper. **gas_layer_data** is an
  alternative that loads a precomputed table of pressures and tangential
  stresses to apply to the liquid; it was used during development of the code but
  is not referenced in the paper.

- **fileinfo** – This class parses the text configuration files to extract the
  simulation parameters. It also converts various non-dimensional constants
  supplied in the configuration file in dimensional constants to be used in the
  simulation.

- **bicubic_interp** – This class performs bicubic interpolation on a regular
  two-dimensional grid of data. It is only used by the **gas_layer_data**
  class.

- **tri_solve** – This class represents a tridiagonal linear system and contains
  a routine for solving it via the LAPACK function `dgtsv`.

The code is commented in the style of [Doxygen](https://www.doxygen.nl/), where
each class and function has a special comment block beginning with `/**` that
describes it.

The simulation method requires three linear systems to be solved during each
timestep. The code contains classes that describe these linear systems, which
are used by the TGMG library. There is a base class called **mgs_common**
that contains definitions and functions that are common to all classes.
The three linear systems are then defined via three classes:

- **mgs_mac** describes the linear system for the marker-and-cell (MAC)
  projection that is used an intermediate step in the Navier–Stokes update. 

- **mgs_fem** describes the linear system for the approximate projection
  implemented using a finite-element discretization [12] using a bilinear hat
  function basis [6,7].

- **visco_impl** describes the linear system arising from the liquid viscous
  term using an implicit Crank–Nicolson discretization. The code creates two
  separate versions of this class for solving for the *x* and *y* velocity
  components. The linear systems for the two components are near-identical
  apart from a small difference in the boundary condition implementation.

## Known issues
The code has been tested over a wide range of different parameter choices.
However, parameter choices that are outside those reported in the paper may
cause the code to fail. This usually manifests as the multigrid method
of Newton method (inside the **gas_layer_solve** class) failing to converge
to the required tolerance.

Two known regimes with problems are follows:

- If the initial velocity is low, and the liquid viscosity is low, then
  the height profile may develop capillary waves and progressively sharp
  features as shown in Fig. 15 of the paper.

- For very low viscosities (*i.e.* below 2 cSt), a numerical instability
  emerges at the far end of the grid.

## Bibliography
1. Arthur M. Worthington, *XXVIII. On the forms assumed by drops of liquids
   falling vertically on a horizontal plate*, P. R. Soc. London **25**, 261–272
   (1877).
   [doi:10.1098/rspl.1876.0048](https://doi.org/10.1098/rspl.1876.0048)

2. Christophe Josserand and Sigurdur T. Thoroddsen, *Drop impact on a solid
   surface*, Annu. Rev. Fluid Mech. **48**, 365–391 (2016).
   [doi:10.1146/annurev-fluid-122414-034401](https://doi.org/10.1146/annurev-fluid-122414-034401)

3. Shreyas Mandre, Madhav Mani, and Michael P. Brenner, *Precursors to splashing
   of liquid droplets on a solid surface*, Phys. Rev. Lett.,
   102(13):134502, 2009.
   [doi:10.1103/PhysRevLett.102.134502](https://doi.org/10.1103/PhysRevLett.102.134502)

4. Mahdav Mani, Shreyas Mandre, and Michael P. Brenner, *Events before droplet
   splashing on a solid surface*, J. Fluid Mech. **647**, 163–185
   (2010).
   [doi:10.1017/S0022112009993594](https://doi.org/10.1017/S0022112009993594)

5. John M. Kolinski, L. Mahadevan, and Shmuel M. Rubinstein, *Lift-off
   instability during the impact of a drop on a solid surface*, Phys. Rev.
   Lett. **112**, 134501 (2014).
   [doi:10.1103/PhysRevLett.112.134501](https://doi.org/10.1103/PhysRevLett.112.134501)

6. Jiun-Der Yu, Shinri Sakai, and James A. Sethian, *A coupled level set
   projection method applied to ink jet simulation*, Interface. Free Bound.
   **5**, 459–482 (2003).
   [doi:10.4171/IFB/87](https://doi.org/10.4171/IFB/87)

7. Chris H. Rycroft, Chen-Hung Wu, Yue Yu, and Ken Kamrin, *Reference map
   technique for incompressible fluid–structure interaction*, J. Fluid Mech.
   **898**, A9 (2020).
   [doi:10.1017/jfm.2020.353](https://doi.org/10.1017/jfm.2020.353)

8. James W. Demmel, *Applied Numerical Linear Algebra*, SIAM (1997).
   [doi:10.1137/1.9781611971446](https://doi.org/10.1137/1.9781611971446)

9. William L. Briggs, Van Emden Henson, and Steve F. McCormick, *A Multigrid
   Tutorial, Second Edition*, SIAM (2000).
   [doi:10.1137/1.9780898719505](https://doi.org/10.1137/1.9780898719505)

10. David L. Brown, Ricardo Cortez, and Michael L. Minion, *Accurate projection
   methods for the incompressible Navier–Stokes equations*, J. Comput. Phys.
   **168**, 464–499 (2001).
   [doi:10.1006/jcph.2001.6715](https://doi.org/10.1006/jcph.2001.6715)

11. Robert Saye, *Interfacial gauge method for incompressible fluid dynamics*,
   Sci. Advances **2**, e1501869 (2016).
   [doi:10.1126/sciadv.1501869](https://doi.org/10.1126/sciadv.1501869)

12. Ann S. Almgren, John B. Bell, and William G. Szymczak, *A numerical method
   for the incompressible Navier–Stokes equations based on an approximate
   projection*, SIAM J. Sci. Comput. **17**, 358–369 (1996).
   [doi:10.1137/S1064827593244213](https://doi.org/10.1137/S1064827593244213)
