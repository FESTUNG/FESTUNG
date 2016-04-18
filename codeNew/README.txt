Implementation of a MATLAB / GNU Octave solver for 2D Shallow-Water Equations

by Hennes Hajduk

SYSTEM OF EQUATIONS

This code solves the following system of 2-Dimensional Shallow-Water Equations using the Local Discontinuous Galerkin Method.

	c_t + div B(c) = M(c)		(*)
	
with the primary unknowns c = (c1, c2, c3)^T = (H, uH, vH)^T, where uH = u*H and vH = v*H,

    B(c) = [ c2														 , c3														 ;
						 c2*c2/c1 + 0.5 * gConst * c1^2, c2*c3/c1											 ;
						 c3*c2/c1					   					 , c3*c3/c1 + 0.5 * gConst * c1^2 ]
			 
	M(c) = [ F0;
					 F1 - gConst * (\partial_{x^1} z_b) c1 - \tau_{bf} * c2 + fc * c3;
					 F2 - gConst * (\partial_{x^2} z_b) c1 - \tau_{bf} * c3 - fc * c2 ].

The physical meaning of H is the total height of the water column. 
u and v are velocity components in x- and y-direction while uH and vH are corresponding momenta.
The given parameters and functions are as follows:
gConst is the acceleration due to gravity, 
\tau_{bf} = tauConst * (c2^2 + c3^2)^0.5 / c1 is the bottom friction coefficient
and fc is a given function representing the Coriolis coefficient.
z_b is the function parametrizing the sea bed (z_b should be negative or zero).
F0, F1 and F2 are given source terms, F0 = 0 is usually used, but it is possible to use a source term here as well.
(*) is a system of three equations, the first line of (*) is called the continuity equation, the second and third are called momentum equations.

BOUNDARY CONDITIONS

The code supports so called land boundary conditions, where with \nu being the outer normal unit vector of the boundary one has

	[u, v] * \nu = 0    on land boundaries    => [c2, c3] * \nu = [uH, vH] * \nu = H * [u, v] * \nu = 0	on land boundaries,

and open sea boundary conditions, where the so called free surface elevation \xi = H + z_b is prescribed

	\xi = \xi_{OS}		on open sea boundaries, 

from which one gets a prescribed height H_{OS} = \xi_{OS} - z_b.

UNKNOWNS

In each time step of the algorithm the main unknowns H, uH and vH are computed.
Further unknowns are the velocities u and v which can be approximated by projecting the quotients uH/H and vH/H into the discrete space. This is later done in the routine getDGvelocities.
Also instead of the height of water, it is usually more interesting to know the free surface elevation \xi. Since z_b is given, it can be easily calculated from H.

TIME STEPPING

Time-stepping is done by either the forward Euler or a semi-implicit combination of forward and backward Euler scheme, in which all non-linearities are discretized explicitly and all linear terms are discretized implicitly with regard to time.
This scheme combines the advantages of forward and backward schemes, since we do not need to solve a non linear system of equations in every time step.
Furthermore, one hopes to achieve better stability properties than in the forward Euler scheme, while not being forced to use extremely small time steps to have a numerically stable solution.
The routine Euler performs one step of the time stepping loop.

RIEMANN SOLVER

It is in the nature of DG methods that one has to use certain stabilizing tools to approximate non-linear fluxes over element edges.
The so called Lax-Friedrichs solver is used throughout this implementation one every interior edge and optionally on open sea boundary edges. Here, instead of using for example the left or right limits of the flux B(c) * \nu on the boundary, one inserts

	0.5 * (B(c^-) + B(c^+)) * \nu + 0.5 * \lambda * (c^- - c^+)
	
where

	\lambda = abs( ( (H^0.5 * uH)^- + (H^0.5 * uH)^+ ) * nu^{1} + ( (H^0.5 * vH)^- + (H^0.5 * vH)^+ ) * nu^{2} ) / ((H^-)^0.5 + (H^+)^0.5) + (0.5 * gConst * (H^- + H^+))^0.5.

instead of B(c) * \nu into every edge integral.

NUMERICAL PARAMETERS

- The user needs to specify the refinement level, i.e. how many times the initial grid should be refined by quartering of each element. The variable that specifies this is refinement. It is set to arg1 in case main is handed at least one input argument. refinement = 0 uses the initial grid that domainHierarchy produces.
- Furthermore the polynomial degree of the discrete space needs to be specified by the variable p, which is set to arg2 in case main is handed at least two input arguments.
- If the user wants to use a Riemann solver on Open Sea Boundaries the variable OSRiem needs to be set to true. In case main is handed at least three input arguments OSRiem is set to arg3.
- To distinguish between the two variants of the Euler method the variable scheme has to be set either to 'explicit' or 'semi-implicit'. In case main is handed four input arguments scheme is set to arg4.
  If the routine main is handed between zero and three arguments, default values are set for the remaining parameters.
( Since the user can select one of two available schemes and decide whether or not to use a Riemann solver stabilisation on open sea boundaries, four different numerical variants of the algorithm are supported. All these variants are implemented in the routine Euler. )

DIVISION OF ALGORITHM

In order to keep the main routine simple and easy to adapt for changes the following auxiliary routines provide a further division of the basic algorithm: solutionProcessing, problemVariables, problemGrid and problemTimeStepControl. These routines help to specify the problem parameters such as initial and boundary conditions and the grid that is used. Their content could just as easily be inserted in main, however the implemented approach makes main a lot shorter. Each of those routines features a short description of its task. 

PROBLEM DATA

To make problem data input like problem parameters and grid data easier, the routines solutionProcessing, problemVariables, problemGrid and problemTimeStepControl have an input argument called problem. Right now the only accepted value for this is 1. If further problems are to be solved, it is possible to extend those routines and add different output arguments that correspond to the problem at hand. They can be called by adding another problem specification corresponding to different values of the variable problem. The only thing that needs to be changed in the routine main are the input arguments of those routines.

VARIABLES INPUT

The algorithm needs the following initial and boundary conditions as well as parameters

- zbAlg: a function handle with two arguments x1, x2 (not time-dependent) parametrizing the sea bed, should be of non-positive value,
- xi0, u0, v0: function handles with two arguments x1, x2 (not time-dependent) parametrizing the initial conditions of the free surface elevation and velocities,
- xiOSAlg function-handle with three arguments x1, x2, t (time-dependent) specifying the free surface elevation on Open Sea Boundaries,
- fcAlg: a function handle with two arguments x1, x2 (not time-dependent) parametrizing the Coriolis coefficient,
- gConst: a variable specifying the acceleration due to gravity,
- tauConst: a variable specifying the constant on which the bottom friction depends,
- F0Alg, F1Alg and F2Alg: function handles with three arguments x1, x2, t (time-dependent) which specify the right hand sides of continuity equation (if needed) and the momentum equations.
- Also function handles for xi, u and v are set. If the analytical solution is available, these are used for error estimation. If the solution is unknown they are set to an arbitrary value, since they are unused.

All these are set in the auxiliary routine problemVariables.

NON-ZERO RIGHT HAND SIDES OF THE CONTINUITY EQUATION

It may be handy to use non-zero right-hand sides of the continuity equation for some applications or test cases. If this is intended, one needs to set the variable F0 to true and specify a function_handle F0Alg to represent the right-hand-side function of the continuity equation. It is handy to do this in the auxiliary routine problemVariables.

GRID INPUT

Right now the code is adapted for convergence analysis. The routine domainHierarchy is used to build an initial unstructured grid. This grid consists of elements of diameter at most initialSize. The variable refinement specifies how many times the initial grid will be refined by quartering of every element.
For different scenarios the domain must be adapted to the application. For that purpose circular and all kinds of 2-dimensional polygonal grids can be created by domainCircle, domainPolygon and domainSquare. 
If the grid and grid information should be visualized one has to set the variables isVisGrid to true in the auxiliary routine solutionProcessing. The grid and land boundary conditions are specified in the routine problemGrid.
A grid needs to be created according to the problem that wants to be simulated and the user has to identify the edges of the grid, where there shall be a land boundary condition. This is done in the routine problemGrid which also accepts an input argument called structure. If structure is false an unstructured grid is created. If structure is true a Friedrichs-Keller triangulation is used.
Since the boundary edges that are not land boundaries are open sea boundaries, the code is able to automatically specify open sea boundaries. 

TIMESTEP CONTROL INPUT

For every simulation the end time tEnd and the number of time steps numSteps need to be specified. This is done in the auxiliary routine problemTimeStepControl. It also sets the size of the constant time increment and the variable output (see VISUALIZATION).

EXACT LINEAR DG REPRESENTATION OF z_b

For applications z_b is only known from linear interpolation allowing to use an exact representation of z_b. Therefore zbExact, the corresponding coefficients of the projection into the linear discrete space, is computed using exactly three degrees of freedom per element. Before that the basis functions must be evaluated in each quadrature point. This is done by calling computeBasesOnQuad(3), then finding the coefficients of zbExact and finally calling computeBasesOnQuad(N) for all further usage, where N is the number of local degrees of freedom. The only aim of this approach is to get an exact DG representation zbExact of zbAlg for the locally constant discrete space. For all higher orders of polynomials this exact representation is guaranteed anyway. For error calculation and visualization zbDG is used which holds the DG coefficients of an approximation of zbAlg which is only exact for p > 0.

ASSEMBLY

The different contributions of the discrete versions of the PDE are assembled in various routines. Some of those contributions are time-independent and need to be assembled only once. Others depend on the solution and therefore need to be build in every time step. 
The contributions of all linear terms of the PDE as well as of those in which an unknown appears in a polynomial form can be assembled by computing integrals of non-solution dependent functions using representative tensors only of the size of local degrees of freedom. This approach is not possible for all other non-linearities. However the idea of how to solve this problem is to compute similar local tensors and evaluate the functions that are to be integrated in each quadrature point multiplied by the corresponding quadrature weight. Therefore those tensors have one extra dimension for all the quadrature points.
In every integral there appears a test function \varphi_{ki} as well as another function \varphi_{kj}, which is obtained by factoring out one component of the solution and inserting the linear combination there for the solution. The remaining non-linearity can then be evaluated in every quadrature point and multiplied with the matrix that holds the entries of the local tensor only for the corresponding quadrature point. 
These local tensors are PhiPhi2D and gradPhiPhi2D for area integral contributions as well as PhiPhidiag and PhiPhioffdiag for boundary edge integrals respectively.
All routines which assemble global matrices and vectors feature a short description.
For various non-linearities the unknowns need to be evaluated in quadrature points of different kind. Since all non-linearities use the previous values of the unknowns, these can be computed easily only once for every time step. This is done in the routine computeSolutionOnQuad.
Since some of the matrices, which are used in the case of no Riemann solver stabilisation on open sea boundaries, reappear in case of such a stabilisation, they are used there, too. However some are scaled with a factor of 0.5. The additional matrices and vectors that are needed for the stabilisation are computed only if needed in the routine additionalOSRiem. These matrices and vectors get the suffix 'Riem' attached to their name. This only indicates that they are not needed in case of no additional stabilisation. There is no more logic to their name, as for instance every matrix or vector that uses Dirichlet boundary conditions features the suffix 'Riem' at its name.

CORRCECTION OF HEIGHT

It is vital for the code to work properly, that H is positive and bounded from below by a positive constant. This property is checked for every element in its vertices. If it does not hold the routine heightCorrection determines a local linear function which is added to H such that after the correction the values in every vertex of every element are at least as high as the user specified tolerance minTol. If such a correction is used, a warning is printed which indicates at what time step a correction is needed and what its norm is.

VISUALIZATION

The variables isVisu, isVisuH, isVisv, isVisvH, isVisxi specify which unknowns are to be visualized. These variables are set in the auxiliary routine solutionProcessing.
One can use the variable output to plot the solution only at time steps that are zero modulo output. output is set in the auxiliary routine problemTimestepControl in such a way that exactly numPlots (or numSteps if numSteps < numPlots) plots are created (not counting the plot of the initial conditions). numPlots is an input argument of problemTimeStepControl.
Visualization is only supported for discrete spaces of polynomial order 2 or lower.

CODE VERIFICATION

To verify the code by comparing analytical and numerical solutions set the variable isSolAvail to true in the auxiliary routine solutionProcessing. Then for the last time step of the calculation the L2-Errors for each unknown are estimated.
Right now an analytical test case is used in main. It is called by the problem identification 1. The order of convergence is approximately p+1 if p = 0,1 is the local polynomial order of the discrete space. Since very small time increments are needed for optimal convergence for higher polynomial orders, those have not been explored to a full scale. However a convergence order of at least 1 could be observed for all higher polynomial orders and moderate time increments. Smaller time increments should yield an optimal convergence order of p+1.

FURTHER ROUTINES

- One can call the routine errorCalculation with input arguments p, refinement, OSRiem and scheme. It estimates the L2-errors and convergence orders for \xi, u, v, uH, vH by simulating the given application starting with the coarse initial grid and refining it successively. The input arguments are the same as the user input variables in main and do not need to be specified since default values are set. It also saves the errors and convergence orders for each unknown as well as the other workspace variables into "'p = ' <p> ', OS-Stab = <OSRiem>, scheme = ' <scheme>".
- The routine errorCalculation uses another routine called orderOfConvergence, which accepts an array of numerical errors and corresponding maximal sizes of grid elements and returns the convergence orders estimated from to consecutive refinement levels.
- A script called simulate is available, which calculates the errors and convergence orders of the algorithm for all four different values of numerical parameters for a given polynomial . It also estimates the runtime of each routine.