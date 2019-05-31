% Helper routine for testConvergence.m that carries out convergence tests
% for problem  'swe_2dv'.

%===============================================================================
%> @file
%>
%> @brief Helper routine for testConvergence.m that carries out convergence tests
%>        for problem  @link swe_2dv @endlink.
%===============================================================================
%>
%> @brief Helper routine for testConvergence.m that carries out convergence tests
%>        for problem  @link swe_2dv @endlink.
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
%>                      @link swe_2dv/getTestcase.m @endlink for available options.
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
function [err, eoc] = doConvergenceTest(testcase, pLevel, hLevel, tLevel)
if nargin < 2
  pLevel = 0:4;
end % if
if nargin < 3
  hLevel = 1:6;
end % if
if nargin < 4
  tLevel = mat2cell((2.^pLevel)' * (250/4 * 2.^(2*hLevel)), ones(1, length(pLevel)), length(hLevel));
end % if

if ~iscell(tLevel)
  tLevel = cellfun(@(c) tLevel, cell(size(pLevel)), 'UniformOutput', false);
end % if
nTimeLevel = cell2mat(cellfun(@(c) length(c), tLevel, 'UniformOutput', false));

isSpatConv = length(hLevel) > 1 && all(1  == cell2mat(cellfun(@(c) length(c), tLevel, 'UniformOutput', false)));
isTimeSpatConv = all(length(hLevel) == nTimeLevel);
assert(xor(isSpatConv, isTimeSpatConv), 'Time and space convergence levels must fit to each other!');
nLevel = max(length(hLevel), max(nTimeLevel));

err = cell(size(pLevel));
eoc = cell(size(pLevel));

for ip = 1 : length(pLevel)
  for level = 1 : nLevel
    pd = struct('isVisSol', false, 'isVisGrid', false, 'testcase', testcase, 'tEnd', 5);
    pd.p = pLevel(ip);
    if isSpatConv || isTimeSpatConv
      pd.numElem = [2^hLevel(level), 2^(hLevel(level)-1)];
    end % if
    if ~isSpatConv
      pd.numSteps = tLevel{ip}(level);
    end % if
    try
      pd = main('swe_2dv', pd);
      err{ip}(level, :) = pd.error;
      N = size(err{ip}, 1);
      if isSpatConv || isTimeSpatConv
        eoc{ip} = [ zeros(1,3); log(err{ip}(1:N-1, :) ./ err{ip}(2:N, :)) ./ ...
                        repmat((hLevel(2:N) - hLevel(1:N-1))' * log(2), 1, 3) ];
      else
        eoc{ip} = [ zeros(1,3); log(err{ip}(1:N-1, :) ./ err{ip}(2:N, :)) ./ ...
                        repmat((tLevel{ip}(2:N) - tLevel{ip}(1:N-1))' * log(2), 1, 3) ];
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
fprintf('Err(h)    EOC(h)   Err(u)    EOC(u)   Err(w)    EOC(w)\n');
fprintf('=======================================================\n');

for ip = 1 : length(err)
  N = size(err{ip}, 1);
  fprintf('------------------------ p = %d ------------------------\n', pLevel(ip)); 
  for i = 1 : N
    fprintf('%6.2e %6.2f    %6.2e %6.2f    %6.2e %6.2f\n', ...
      err{ip}(i,1), eoc{ip}(i,1), err{ip}(i,2), eoc{ip}(i,2), err{ip}(i,3), eoc{ip}(i,3));
  end % for i
end % for p
end % function
