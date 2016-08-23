% Modifies a discrete function such that its vertex values in each element do 
% not fall below a certain threshold.
%
%===============================================================================
%> @file correctMinValueExceedanceDisc.m
%>
%> @brief Modifies a discrete function such that its vertex values in each 
%>        element do not fall below a certain threshold.
%===============================================================================
%>
%> @brief Modifies a discrete function such that its vertex values in each 
%>        element do not fall below a certain threshold.
%>
%> This routine computes the vertex values of a discrete function and changes 
%> them to a treshold value for every vertex where they are below that treshold.
%>
%> @param dataDisc     A representation of the discrete function, e.g., as 
%>                     computed by <code>projectFuncCont2DataDisc()</code>.
%> @param corrSys      The values of the first three basis functions in each of 
%>                     the vertices on the reference element (only used for 
%>                     polynomial order at least 1). @f$[3 \times 3]@f$
%> @param nStep        The current iteration number of the main loop. 
%> @param minTol       The treshold values for each local vertex.
%>                     @f$[K \times 3]@f$ or scalar.
%> @param corrTol      The maximal value for which computation continues without
%>                     throwing an error.
%> @retval dataDisc    The possibly modified representation of the discrete 
%>                     function. @f$[K \times N]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Hennes Hajduk, Balthasar Reuter, Florian Frank,
%>                      Vadym Aizinger
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
%>
function dataDisc = correctMinValueExceedanceDisc(dataDisc, corrSys, nStep, minTol, corrTol)
if nargin < 3
	nStep = 0;
end % if
if nargin < 4
	minTol = 0;
end % if
if nargin < 5
	corrTol = Inf;
end % if

dataV0T = projectDataDisc2DataLagr(dataDisc, 1);
corr = max(minTol - dataV0T, 0);
maxCorr = max(corr(:));

if maxCorr > corrTol
  [indx, indy] = find(corr == maxCorr, 1);
  error([ 'Unknown at local vertex ' num2str(indy) ' of element ' num2str(indx) ...
          ' in step ' num2str(nStep) ' is ' num2str(dataV0T(indx, indy)) ...
          ' (below the minimum tolerance of ' num2str(corrTol) ').' ]);
end % if

if any(corr(:))
  [indx, indy] = find(corr == maxCorr, 1);
  warning([ 'A maximum value of ' num2str(maxCorr) ...
            ' had to be added to unknown at local vertex ' num2str(indy) ...
            ' of element ' num2str(indx) ' in step ' num2str(nStep) '.' ]);
  if size(dataDisc,2) == 1
    dataDisc = dataDisc + max(corr, [], 2) / phi(1,1/3,1/3);
  else
    dataDisc(:,1:3) = dataDisc(:,1:3) + corr / corrSys;
  end % if
end % if
end % function
