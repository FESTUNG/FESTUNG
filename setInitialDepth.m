% Sets the initial values for the bathymetry.

%===============================================================================
%> @file sweInverse/setInitialDepth.m
%>
%> @brief Sets the initial values for the bathymetry.
%===============================================================================
%>
%> @brief Sets the initial values for the bathymetry.
%>
%> @param  pd           A struct with problem parameters and precomputed
%>                      fields, as provided by configureProblem() and 
%>                      preprocessProblem(). @f$[\text{struct}]@f$
%>
%> @retval ret					The discrete data representing the bathymetry.
%>											@f$[K \times N]@f$
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
function ret = setInitialDepth(pd)
K = pd.K;
N = pd.N;
ret = zeros(K, N);
switch pd.initType
	case 'exact'
		ret = pd.zbExact;
	case 'constant'
		ret(:,1) = -pd.constValue / sqrt(2);
  case 'linear'
    ret = projectFuncCont2DataDisc(pd.g, @(x1,x2) 0.001* x1 - 2, 2*pd.p, pd.refElemPhiPhi, pd.basesOnQuad);
%     ret = projectFuncCont2DataDisc(pd.g, @(x1,x2) -x1/1000-2, 2*pd.p, pd.refElemPhiPhi, pd.basesOnQuad);
%     ret = projectFuncCont2DataDisc(pd.g, @(x1,x2) -0.000175*sqrt(x1.^2+x2.^2)+25, 2*pd.p, pd.refElemPhiPhi, pd.basesOnQuad);
  case 'error'
    [~, depth] = domainADCIRC(['fort_' pd.name '10%.14'],['fort_' pd.name '.17'],5);
    zbCont = @(x1,x2) evaluateFuncFromVertexValues(pd.g, -depth, x1,x2);
    
    assert(max( max( abs(depth(pd.g.V0T) + zbCont(pd.g.coordV0T(:,:,1), pd.g.coordV0T(:,:,2))) ) ) < 1.e-5, ...
           'Bathymetry incorrectly constructed!');
    ret = projectFuncCont2DataDisc(pd.g, zbCont, 2*pd.p, pd.refElemPhiPhi, pd.basesOnQuad);
  otherwise
		error('Unknown type of data initialisation.');
end % switch
if ~strcmp(pd.initType, 'exact')
%   ret = ret + setInitialBoundaryDepth(pd, constValue);
end % if
end % function

function ret = setInitialBoundaryDepth(pd)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% copied from slope limiting %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pd.g.markV0TbdrOS = ismember(pd.g.V0T, pd.g.V0E(pd.g.E0T(pd.g.markE0TbdrOS), :));

xiV0T = -pd.depth(pd.g.V0T) + pd.constValue;

% all vertex indices for which values have to be set
nonUniqueVertices = pd.g.V0T(pd.g.markV0TbdrOS);
% the vertex values for triangles with boundary edges and (possibly)
% NaN if it is not a boundary vertex
xiV0T = xiV0T(pd.g.markV0TbdrOS);

% add the values of all contributing vertices via matrix-vector product
pd.vertInd2VertIndUniqueOS = bsxfun(@eq, (1:pd.g.numV).', nonUniqueVertices.');

xiV = pd.vertInd2VertIndUniqueOS * setNaN2Zero(xiV0T);
pd.xiVCountOS = pd.vertInd2VertIndUniqueOS * double(~isnan(xiV0T));

% xiV contains the sum of the vertex values of all triangles and is
% divided by the number of contributing elements in dataVCountOS
xiV = xiV ./ pd.xiVCountOS;
pd.xiV0Tos = xiV(pd.g.V0T);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sys = pd.sysMinValueCorrection.';
ret = zeros(pd.K, pd.N);
for k = 1:pd.g.numT
  ret(k,1:3) = sys \ setNaN2Zero(pd.xiV0Tos(k,:).');
end % for
end % function
