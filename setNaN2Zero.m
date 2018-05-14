% Technical routine required to elimate possible NaN entries in any input by 
% replacing them with zeros.

%===============================================================================
%> @file
%>
%> @brief Technical routine required to elimate possible NaN entries in any 
%>        input by replacing them with zeros.
%===============================================================================
%>
%> @brief Technical routine required to elimate possible NaN entries in any 
%>        input by replacing them with zeros.
%>
%> This routine can be used for any application where the algorithm produces 
%> NaN entries in a tensor of any order alongside realistic values. Those NaN 
%> entries are typically required to be set to zero.
%> An example of applicability in our framework is the divison of zero by zero.
%> While these NaN entries are supposed to be eliminated by multiplying with a
%> logical field that has zero entries in a way that each NaN is multiplied by
%> zero, MATLAB returns NaN as the product of zero and NaN. Thus, this routine
%> is necessary to modify data that possibly includes NaN entries.

%> @param  input        Double-valued tensor of arbitary order.
%>
%> @retval ret          The same tensor with zeros instead of NaN's.
%>
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Hennes Hajduk, Balthasar Reuter, Florian Frank, 
%>                      Vadym Aizinger
%>                      
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
function ret = setNaN2Zero(input)
ret = max(input, 0) + min(input,0);
end % function
