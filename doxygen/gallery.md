# Gallery

## Data visualization

![Elementwise constant, linear, and quadratic approximation.](doxygen/images/visualization.png "Elementwise constant, linear, and quadratic approximation.")

The routine [visualizeDataLagr.m](visualizeDataLagr_8m.html) writes simulation data such as solutions or coefficients to VTK- and Tecplot files, which can be visualized using, e.g., [Paraview](http://www.paraview.org/) or [Tecplot 360](http://www.tecplot.com/).  

## Shallow-water equations

![Galveston Bay model domain with computational mesh after a 10 days simulation.](doxygen/images/sweBathymetry.png "Galveston Bay model domain with computational mesh.")

[comment]:<> ([velocity magnitude](doxygen/images/sweVelocity.png), [surface elevation](doxygen/images/sweBathymetry.png))

\f[ \frac{\partial H}{\partial t} + \nabla\cdot(\vec{u}\,H) ~=~ 0 \f]
\f[ \frac{\partial (u\,H)}{\partial t} + \nabla\cdot(u\,\vec{u}\,H)
+ g \,H \,\frac{\partial \xi}{\partial x} + \tau_\mathrm{bf}\,u\,H 
 -  f_\mathrm{c} \,v\,H ~=~ H\,f_x + F_x \f]
\f[ \frac{\partial (v\,H)}{\partial t} + \nabla\cdot(v\,\vec{u}\,H)
+ g \,H \,\frac{\partial \xi}{\partial y} + \tau_\mathrm{bf}\,v\,H 
 +  f_\mathrm{c} \,u\,H ~=~ H\,f_y + F_y \f]

### Unknowns: 
\f$H=\xi - b\f$ (total height of water), with \f$\xi\f$ free surface elevation, \f$b\f$ bathymetry; \f$ \vec{u} = [u,v]^\mathrm{T}\f$ (depth-averaged horizontal velocity). 

### Parameters:
 \f$g\f$ (acceleration due to gravity), \f$\tau_\mathrm{bf}\f$ (friction coefficient), \f$f_\mathrm{c}\f$ (Coriolis constant), and force terms \f$\vec{f} =
[f_x,f_y]^\mathrm{T}\f$ and \f$\vec{F} = [F_x,F_y]^\mathrm{T}\,\f$.

### Features:
* Well-balanced DG scheme of orders \f$p = 0,\hdots,4\f$ implemented.
* Lax-Friedrichs Riemann-solver for flux approximation.
* SSP-Runge-Kutta methods used for time-stepping, cf. [rungeKuttaExplicit.m](rungeKuttaExplicit_8m.html).
* Slope-limiting framework incorporated, cf. [applySlopeLimiterTaylor.m](applySlopeLimiterTaylor_8m.html).
* Additionally to the usual FESTUNG features an interface to ADCIRC for grid-generation and configuration is provided.
* Unstructured and spherical grids supported.
* Matrix-assembly is only required during pre-processing which makes the solver extremely fast.
* Flexible coding strategies are used to allow easy modifications.
* Multiple model and numerical parameters, such as the choice between linear and nonlinear bottom-friction or the possibility to turn on/off Riemann solver usage for boundary conditions.
* In total five different boundary types are implemented.
