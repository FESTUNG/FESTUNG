# Gallery

## Data visualization

![Elementwise constant, linear, and quadratic approximation.](doxygen/images/visualization.png "Elementwise constant, linear, and quadratic approximation.")

The routine [visualizeDataLagr.m](visualizeDataLagr_8m.html) writes simulation data such as solutions or coefficients to VTK- and Tecplot files, which can be visualized using, e.g., [Paraview](http://www.paraview.org/) or [Tecplot 360](http://www.tecplot.com/).  

## Shallow-water equations

Author: [Hennes Hajduk](http://www.mathematik.tu-dortmund.de/de/personen/person/Hennes+Hajduk.html)

References: [[HHAR2018]](@ref HHAR2018)

![Galveston Bay model domain with computational mesh, velocity magnitude, and surface elevation after a 10 days simulation.](doxygen/images/sweBVE.png "Galveston Bay model domain with computational mesh.")

### Model equations
\f[ \frac{\partial H}{\partial t} + \nabla\cdot(\vec{u}\,H) ~=~ 0 \f]
\f[ \frac{\partial (u\,H)}{\partial t} + \nabla\cdot(u\,\vec{u}\,H)
+ g \,H \,\frac{\partial \xi}{\partial x} + \tau_\mathrm{bf}\,u\,H 
 -  f_\mathrm{c} \,v\,H ~=~ H\,f_x + F_x \f]
\f[ \frac{\partial (v\,H)}{\partial t} + \nabla\cdot(v\,\vec{u}\,H)
+ g \,H \,\frac{\partial \xi}{\partial y} + \tau_\mathrm{bf}\,v\,H 
 +  f_\mathrm{c} \,u\,H ~=~ H\,f_y + F_y \f]

### Unknowns
\f$H=\xi - b\f$ (total height of water), with \f$\xi\f$ free surface elevation, \f$b\f$ bathymetry; \f$ \vec{u} = [u,v]^\mathrm{T}\f$ (depth-averaged horizontal velocity). 

### Parameters
 \f$g\f$ (acceleration due to gravity), \f$\tau_\mathrm{bf}\f$ (friction coefficient), \f$f_\mathrm{c}\f$ (Coriolis constant), and force terms \f$\vec{f} =
[f_x,f_y]^\mathrm{T}\f$ and \f$\vec{F} = [F_x,F_y]^\mathrm{T}\,\f$.

### Features
* Well-balanced DG scheme of orders \f$p = 0,\ldots,4\f$ implemented.
* Lax-Friedrichs Riemann-solver for flux approximation.
* SSP-Runge-Kutta methods used for time-stepping, cf. [rungeKuttaExplicit.m](rungeKuttaExplicit_8m.html).
* Slope-limiting framework incorporated, cf. [applySlopeLimiterTaylor.m](applySlopeLimiterTaylor_8m.html).
* In total five different boundary types are implemented.
* Multiple model and numerical parameters, such as the choice between linear and nonlinear bottom-friction or the possibility to turn on/off Riemann solver usage for boundary conditions.
* Flexible coding strategies are used to allow easy modifications.
* Additionally to the usual FESTUNG features an interface to ADCIRC for grid-generation and configuration is provided.
* Unstructured and spherical grids supported.
* Matrix-assembly is only required during pre-processing which makes the solver extremely fast.

## Coupled free-surface and subsurface flow

Author: Balthasar Reuter, Andreas Rupp

References: [[RRAK2019]](@ref RRAK2019)

![Coupled free-surface and subsurface flow in a vertical slice with flow from left to right. Arrows indicate velocity direction and magnitude in the free-flow domain, colors and isolines indicate hydraulic head in the subsurface domain.](doxygen/images/coupledSweDarcy.png "Coupled free-surface and subsurface flow in a vertical slice.")

![Computational grid for the coupled simulation.](doxygen/images/coupledSweDarcyGrid.png "Computational grid for the coupled simulation")

### Model equations
Flow in the free-surface domain is modeled by the primite hydrostatic equations (three-dimensional shallow water equations), restricted to a vertical slice:
\f[ \partial_t h(t, x^1) + \partial_{x^1} \left( \int_{z_\mathrm{b}(x^1)}^{\xi (t, x^1)} u^1(t, \vec{x}) \, \mathrm{d} x^2 \right) = 0\f]
\f[ \partial_t u^1(t, \vec{x}) + \nabla\cdot \left( u^1(t, \vec{x}) \, \vec{u}(t, \vec{x}) \right) + g\,\partial_{x^1} h(t, \vec{x})  - \nabla\cdot \left( \mathbf{D}(t,\vec{x}) \, \nabla u^1(t, \vec{x}) \right) = f(t, \vec{x}) - g \, \partial_{x^1} z_\mathrm{b} (t, x^1) \f]
\f[ \nabla\cdot \vec{u}(t, \vec{x}) = 0 \f]

Flow in the subsurface domain is modeled by the equation for groundwater flow, based on Darcy's law:
\f[S_0 \partial_t \tilde{h} - \nabla \cdot \left(\mathbf{K} \, \nabla \tilde{h}\right) = f\f]

Coupling at the interface is done using conservation of mass and continuity of pressure:
\f[\vec{u}(t,\vec{x}) \cdot \vec{\nu} = \mathbf{K}(t,\vec{x}) \, \nabla \tilde{h}(t,\vec{x}) \cdot \tilde{\vec{\nu}} \f]
\f[\tilde{h}(t,\vec{x}) = \xi(t,x^1) + \frac{1}{2\,g} \left( u^{1} \right)^2\f]

### Unknowns
\f$h=\xi - z_\mathrm{b}\f$ (total height of water), with \f$\xi\f$ free surface elevation, \f$z_\mathrm{b}\f$ bathymetry; \f$ \vec{u} = [u^1,u^2]^\mathrm{T}\f$ (horizontal and vertical velocity components), \f$\tilde{h}\f$ (hydraulic head). 

### Parameters
 \f$g\f$ (acceleration due to gravity), \f$\mathbf{D}\f$ (horizontal and vertical diffusivities), \f$S_0\f$ (specific storativity), \f$\mathbf{K}\f$ (hydraulic conductivity), and force terms \f$f\f$ and \f$\tilde{f}\f$.

### Features
* LDG scheme of orders \f$p = 0,\ldots,4\f$ implemented.
* Lax-Friedrichs Riemann-solver for approximation of horizontal fluxes.
* SSP-Runge-Kutta methods used for time-stepping in the free-flow domain, cf. [rungeKuttaExplicit.m](rungeKuttaExplicit_8m.html).
* Implicit Euler method used for time-stepping in the subsurface domain.
* Dirichlet- and Neumann-type boundary conditions in subsurface domain.
* Land, open-sea, inflow and outflow boundary conditions in free-flow domain.
* Time-dependent top boundary with grid adaptation.
* Unstructured and spherical grids supported.
