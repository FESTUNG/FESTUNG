% Routine that uses the information of the centeroid values of neighbouring
% elements to determine on which elements the gradient of a quantity is
% below a certain tolerance.

%===============================================================================
%> @file
%>
%> @brief Routine that uses the information of the centeroid values of 
%>        neighbouring elements to determine on which elements the gradient of a
%>        quantity is below a certain tolerance.
%===============================================================================
%>
%> @brief Routine that uses the information of the centeroid values of 
%>        neighbouring elements to determine on which elements the gradient of a
%>        quantity is below a certain tolerance.
%>
%> This routine computes the differnece between the maximum and minimum 
%> centeroidal values of each of the elements that share either of its vertices
%> and deactivetes the elements for which this differnece is small.
%> 
%> There are two types of methods available: element-based computes the
%> difference between the maximum of any neighbouring centeroidal value and
%> the minimum of any such value. vertex-based one uses maxima and minima
%> corresponding to a certain vertex and uses the maximum of these
%> differences on each element. 
%>
%> @param  minMaxV0T  Two matrices with minimum or maximum centroid values,
%>                    respectively, of the patch of elements surrounding each
%>                    vertex of each element as computed by 
%>                    <code>computeMinMaxV0TElementPatch()</code>
%>                    @f$[2 \times 1 \mathrm{cell}]@f$
%> @param  maskTol    The tolerance below which the slope on an element has to
%>                    be such that that element is considered inactive. 
%> @param  maskType   String that determines which of the available methods
%>                    should be used.
%> @retval mask       <code>logical</code> array that marks each triangles
%>                    as either active (true) or inactive (false).
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Hennes Hajduk, Balthasar Reuter, Florian Frank, Vadym Aizinger
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
function mask = computeMask(minMaxV0T, maskTol, maskType)
if strcmp(maskType, 'element-based')
  mask = max(minMaxV0T{2},[],2) - min(minMaxV0T{1},[],2) >= maskTol;
elseif strcmp(maskType, 'vertex-based')
  mask = max(minMaxV0T{2} - minMaxV0T{1},[],2) >= maskTol;
else
  error('Invalid type of mask specified.')
end % if
end % function