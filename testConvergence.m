% Helper routine for testConvergence() that carries out convergence tests
% for problem 'advection_implicit'.

%===============================================================================
%> @file
%>
%> @brief Helper routine for testConvergence() that carries out convergence tests
%>        for problem 'advection_implicit'.
%===============================================================================
%>
%> @brief Helper routine for testConvergence() that carries out convergence tests
%>        for problem 'advection_implicit'.
%>
%> It solves a given testcase with continuously refined mesh size and/or
%> time step size and saves the obtained error estimates to compute
%> experimental orders of convergence.
%>
%> Mesh size is given by @f$h_j = \frac{1}{3\cdot 2^j}@f$, with @f$j@f$
%> any of the values in `hLevel`.
%> Number of time steps is given by @f$N_\mathrm{steps} = 10 \cdot 2^j@f$,
%> with @f$j@f$ any of the values in `tLevel`.
%>
%> @param  testcase     The name of the testcase to be solved. See
%>                      getTestcase() for available options.
%> @param  pLevel       Vector with polynomial approximation orders to be tested.
%> @param  hLevel       Vector with mesh refinement levels to be tested.
%> @param  tLevel       Vector with time step refinement levels to be tested.
%>
%> @retval err          A cell array of vectors, with errors from the
%>                      computations on all space/time levels for each
%>                      polynomial degree.
%> @retval eoc          A cell array with the corresponding experimental
%>                      orders of convergence.
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Balthasar Reuter, Florian Frank, Vadym Aizinger
%> @author Balthasar Reuter, 2017.
%>
%> @par License
%> @parblock
%> This program is free software: you can redistribute it and/or modify
%> it under the terms of the GNU General Public License as published by
%> the Free Software Foundation, either version 3 of the License, or
%> (at your option) any later version.
%>
%> This program is distributed in the hope that it will be useful,
%> but WITHOUT ANY WARRANTY; without even the implied warranty of
%> MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%> GNU General Public License for more details.
%>
%> You should have received a copy of the GNU General Public License
%> along with this program.  If not, see <http://www.gnu.org/licenses/>.
%> @endparblock
%
function [err, eoc] = testConvergence(testcase, pLevel, hLevel, tLevel)
if nargin < 2
  pLevel = 0:4;
end % if
if nargin < 3
  hLevel = 1:6;
end % if
if nargin < 4
  tLevel = 1;
end % if

isSpatConv = length(hLevel) > 1 && length(tLevel) == 1;
isTimeSpatConv = length(hLevel) == length(tLevel);
assert(xor(isSpatConv, isTimeSpatConv), 'Time and space convergence levels must fit to each other!');
nLevel = max(length(hLevel), length(tLevel));

err = cell(size(pLevel));
eoc = cell(size(pLevel));

for ip = 1 : length(pLevel)
  for level = 1 : nLevel
    pd = struct('isVisSol', false, 'isVisGrid', false, 'isConvergence', true, 'testcase', testcase);
    pd.p = pLevel(ip);
    if isSpatConv || isTimeSpatConv
      pd.hmax = 2^-hLevel(level) / 3;
    end % if
    if ~isSpatConv
      pd.numSteps = 10 * 2^tLevel(level);
    end % if
    try
      pd = main('advection_implicit', pd);
      err{ip}(level) = pd.error;
      N = length(err{ip});
      if isSpatConv
        eoc{ip} = [ 0, log(err{ip}(1:N-1) ./ err{ip}(2:N)) ./ ...
                        ((hLevel(2:N) - hLevel(1:N-1)) * log(2)) ];
      else
        eoc{ip} = [ 0, log(err{ip}(1:N-1) ./ err{ip}(2:N)) ./ ...
                        ((tLevel(2:N) - tLevel(1:N-1)) * log(2)) ];
      end % if
      printConvergence(err, eoc, pLevel);
    catch ME
      stack = arrayfun(@(a) sprintf('  In %s (line %d)\n', a.file, a.line), ME.stack, 'UniformOutput', false);
      warning('%s: %s\n%s---', ME.identifier, ME.message, [stack{:}]);
      disp(ME.stack)
    end % try
  end % for level
end % for ip


end % function

function printConvergence(err, eoc, pLevel)
fprintf('Err        EOC\n');
fprintf('================\n');

for ip = 1 : length(err)
  N = length(err{ip});
  fprintf('---- p = %d -----\n', pLevel(ip)); 
  for i = 1 : N
    fprintf('%6.2e  %6.3f   \n', err{ip}(i), eoc{ip}(i));
  end % for i
end % for p
end % function