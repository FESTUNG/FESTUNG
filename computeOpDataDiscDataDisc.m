% Computes the sum of two discrete data sets that may be of different order.

%===============================================================================
%> @file computeOpDataDiscDataDisc.m
%>
%> @brief Computes a given operation on two discrete data sets that may be of  
%>        different order.
%===============================================================================
%>
%> @brief Computes a given operation on two discrete data sets that may be of  
%>        different order.
%>
%> This routine prolongs both data sets to the be of the same order and
%> performs a given operation on them.
%>
%> @param op          A function handle that expects two arguments, e.g.,
%>                    <code>plus()</code>.
%> @param dataDisc1   A representation of the discrete function, e.g., as 
%>                    computed by <code>projectFuncCont2DataDisc()</code>
%> @param dataDisc2   A representation of the discrete function, e.g., as 
%>                    computed by <code>projectFuncCont2DataDisc()</code>
%> @retval ret        A representation of the discrete function representing the
%>                    exact sum of the two input functions in the higher order 
%>                    space.
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
function ret = computeOpDataDiscDataDisc(op, dataDisc1, dataDisc2)
[K, N1] = size(dataDisc1);
[~, N2] = size(dataDisc2);

validateattributes(dataDisc1, {'numeric'}, {'size', [K NaN]}, mfilename, 'dataDisc1')
validateattributes(dataDisc2, {'numeric'}, {'size', [K NaN]}, mfilename, 'dataDisc2')

ret = op([dataDisc1, zeros(K,N2-N1)], [dataDisc2, zeros(K,N1-N2)]);
end % function