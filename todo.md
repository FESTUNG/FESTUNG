# Notes for the implementation of a hybridized discontinuous Galerkin scheme in FESTUNG

## TODO

* Find a proper naming scheme for functions and variables
** Should I call surfaces of elements __edges__ or __faces__?
** I should give basis functions on edges/faces a different name than phi.This emphasizes the difference between phi (2D) and basis functions on edges/faces (1D). I think I should stick with the definitions in our other papers and call them mu. I don't think that this letter is used already.
** I need to set up a proper data structure to access edges/faces in an easy and efficient way.
** I need to set up a data structure/function to get the edge/face transformations in an easy/efficient way.
** 