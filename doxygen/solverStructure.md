Solver structure
================


The solver structure of FESTUNG is built around the perception that almost
every solver for a PDE problem can be subdivided into three major steps:

1. A setup phase: define problem parameters, read configuration files, 
   allocate grid data structures, initialize solution vectors, etc.;
2. A solver phase: assemble and solve a linear system (possibly repeatedly
   for different time levels), or apply iterative methods for non-linear
   problems, etc.;
3. A post-processing phase: evaluate errors, write visualization output, etc.

In FESTUNG, these phases are split into some more fine-grained intermediate
steps, allowing to implement solvers for numerical problems in a 
well-structured, read- and maintainable, yet efficient way, as depicted in 
the following flow chart:

![Generic solver formulation, including optional sub-stepping](doxygen/images/solver-structure.png)

Each of these steps is implemented as a MATLAB function, named 
`configureProblem`, `preprocessProblem`, etc.
All these step functions that make up the implementation of a problem solver
are put together in a subfolder.
For an (empty) example, see the folder `template`, which holds a fully
working problem definition with do-nothing steps.

The implemented solver is then executed from the FESTUNG main directory
using the command 

    $ main('template')

where `'template'` is to be replaced by the name of the folder.
