function [minTol, isVisGrid, isVisu, isVisuH, isVisv, isVisvH, isVisxi, isSolAvail] = solutionProcessing(problem)
%% defines miscellaneous parameters on how the solution should be processed, like the minimal height in each vertex,
%% what unknowns are to be visualized and if an analytical solution is available for error calculation
switch problem
  case 1 % analytical test case
    minTol      = 0.0001;            % values of c1 in each vertex must be at least as high as this tolerance 
    isVisGrid   = false;            % visualization of grid
    isVisu      = false;            % visualization of x-velocity component
    isVisuH     = true;            % visualization of x-momentum component
    isVisv      = false;            % visualization of y-velocity component
    isVisvH     = true;            % visualization of y-momentum component
    isVisxi     = true;            % visualization of free surface elevation
    isSolAvail  = true;             % if the solution is available one can compute the error of the computation
                                    % the user must then create function handles H, u and v parametrizing the height and both velocity components
  case 2 % bahamas
    minTol      = 0.001;            % values of c1 in each vertex must be at least as high as this tolerance 
    isVisGrid   = false;            % visualization of grid
    isVisu      = true;             % visualization of x-velocity component
    isVisuH     = false;            % visualization of x-momentum component
    isVisv      = true;             % visualization of y-velocity component
    isVisvH     = false;            % visualization of y-momentum component
    isVisxi     = true;             % visualization of free surface elevation
    isSolAvail  = false;            % if the solution is available one can compute the error of the computation
                                    % the user must then create function handles H, u and v parametrizing the height and both velocity components
  otherwise
    error('Unknown problem.');
end % switch
if minTol <= 0
  warning('minTol should be positive.');
end % if
end %function