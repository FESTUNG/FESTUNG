% Assembles a matrix containing integrals of products of two edge basis
% functions.

%===============================================================================
%> @file assembleMatEdgeMuMu.m
%>
%> @brief Assembles a matrix containing integrals of products of two edge basis 
%>        functions. This corresponds to a mass matrix.
%===============================================================================
%>
%> @brief Assembles a mass matrix @f$M_{\lambda}@f$
%>        containing integrals of products of two edge basis functions.
%>
%> TODO
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  markE0T    A marker indicating whether an edge should be 
%>                    recognized or not. @f$[K \times 3]@f$
%> @param refElemMuMu Local matrix @f$\hat{\mathsf{M}}_{\mu}@f$ as provided
%>                    by <code>integrateRefElemMuMu()</code>.
%>                    @f$[\bar{N} \times \bar{N}]@f$
%> @retval ret        The assembled matrix @f$[\bar{K}\bar{N} \times \bar{K}\bar{N}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Balthasar Reuter, Florian Frank, Vadym Aizinger
%> @author Alexander Jaust, 2017
%> @author Balthasar Reuter, 2017
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
% TODO Check this
% Compute global mas matrix of lambda
% If it is an internal edge, the entry has to be scaled by 2*alpha
% check the element transformation
function ret = assembleMatEdgeMuMu(g, markE0T, refEdgeMuMu)
Nmu = size(refEdgeMuMu, 1);
Kedge = g.numE;
validateattributes(refEdgeMuMu, {'numeric'}, {'size', [Nmu Nmu]}, mfilename, 'refEdgeMuMu');
ret = sparse(Kedge * Nmu, Kedge * Nmu);
for n = 1 : 3
  Kkn = g.areaE0T( :, n ) .*  markE0T(:, n) ;
  ret = ret + kron( sparse( g.E0T(:, n), g.E0T(:, n), Kkn, Kedge, Kedge ), refEdgeMuMu );
end
end % function
