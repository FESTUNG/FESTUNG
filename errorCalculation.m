function [error, order] = errorCalculation(arg1, arg2, arg3, arg4)
%% approximates residuals of approximation for different levels of refinement
%% and calculates numerical order of convergence
%% specification of refinement, order of discrete ansatz space, usage of stabilisation on open sea boundaries and scheme
%% either by dafault values that can be changed here
refinement  = 3;
p           = 1;
OSRiem      = true;
scheme      = 'semi-implicit';
%% or by user input
if nargin > 0
  refinement = arg1;
end % if
if nargin > 1
  p = arg2;
end % if
if nargin > 2
  OSRiem = arg3;
end % if
if nargin > 3
  scheme = arg4;
end % if
error = zeros(refinement+1,5);
h = 0.3*0.5.^[0 : refinement];
for k = 0 : refinement
  [error(k+1,1),error(k+1,2),error(k+1,3), error(k+1,4), error(k+1,5)] = main(k, p, OSRiem, scheme);
end % for
order = cell(5,1);
order{1} = orderOfConvergence(error(:,1),h);
order{2} = orderOfConvergence(error(:,2),h);
order{3} = orderOfConvergence(error(:,3),h);
order{4} = orderOfConvergence(error(:,4),h);
order{5} = orderOfConvergence(error(:,5),h);
%% clearing all unnecessary variables
clear h; clear k;
clear arg1; clear arg2; clear arg3; clear arg4;
%% saving error and convergence tables
if OSRiem
  save(['p = ' num2str(p) ', OS-Stab = true, scheme = ' scheme ]);
else
  save(['p = ' num2str(p) ', OS-Stab = false, scheme = ' scheme ]);
end % if
disp(error);
disp([order{1}, order{2}, order{3}, order{4}, order{5}]);
end % function
