% Performs all post-processing steps, such as error estimates of the final
% solution, etc.

%===============================================================================
%> @file
%>
%> @brief Performs all post-processing tasks, such as error estimates of the 
%>        final solution, etc.
%===============================================================================
%>
%> @brief Performs all post-processing tasks, such as error estimates of the 
%>        final solution, etc.
%>
%> This routine is called after the main loop.
%>
%> It can include things such as error estimates, an output operation of
%> the final solution, etc.
%>
%> @param  problemData  A struct with problem parameters and solution
%>                      vectors. @f$[\text{struct}]@f$
%>
%> @retval problemData  A struct with all necessary parameters and definitions
%>                      for the problem description and precomputed fields.
%>                      @f$[\text{struct}]@f$
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
function problemData = postprocessProblem(problemData)
%% Close waitbar.
if problemData.isWaitbar
  close(problemData.waitbar)
end % if

%% Save stations.
if problemData.isVisStations
  if isfield(problemData, 'dataElev')
    data = problemData.dataElev;
    save(['output' filesep problemData.name '_elev.mat'], 'data')
    fprintf(['Data written to output' filesep problemData.name '_elev.mat\n'])
  end % if
  if isfield(problemData, 'dataVel')
    data = problemData.dataVel;
    save(['output' filesep problemData.name '_vel.mat'], 'data')
    fprintf(['Data written to output' filesep problemData.name '_vel.mat\n'])
  end % if
end % if

%% Compute error if analytical solution available.
p = problemData.p;
qOrd = max(2*p,1);
if problemData.isSolutionAvail
  
  problemData.err = zeros(6,1);
  
  % Continuous solution
  t = problemData.t;
  zbEnd = @(x1,x2) problemData.zbCont(x1,x2,t);
  hEnd = @(x1,x2) problemData.xiCont(x1,x2,t) - problemData.zbCont(x1,x2);
  uEnd = @(x1,x2) problemData.uCont(x1,x2,t);
  vEnd = @(x1,x2) problemData.vCont(x1,x2,t);
  
  % Error in bathymetry (zb)
  problemData.err(1) = computeL2Error(problemData.g, problemData.zbDisc, zbEnd, 2*p+1, problemData.basesOnQuad);
	
	% Error in water height (H)
  hDisc = problemData.cDisc(:,:,1) - problemData.zbDisc;
  problemData.err(2) = computeL2Error(problemData.g, problemData.cDisc(:,:,1), hEnd, 2*p+1, problemData.basesOnQuad);
  
  % Error in primary variable uH
  problemData.err(3) = computeL2Error(problemData.g, problemData.cDisc(:,:,2), @(x1,x2) uEnd(x1,x2) .* hEnd(x1,x2), 2*p+1, problemData.basesOnQuad);
  
  % Error in primary variable vH
  problemData.err(4) = computeL2Error(problemData.g, problemData.cDisc(:,:,3), @(x1,x2) vEnd(x1,x2) .* hEnd(x1,x2), 2*p+1, problemData.basesOnQuad);
  
  % Error in x-velocity (u)
  dataQ0T = (problemData.cDisc(:,:,2) * problemData.basesOnQuad.phi2D{qOrd}.') ./ (hDisc * problemData.basesOnQuad.phi2D{qOrd}.');
  dataDisc = projectDataQ0T2DataDisc(dataQ0T, 2*p, problemData.refElemPhiPhi, problemData.basesOnQuad);
  problemData.err(5) = computeL2Error(problemData.g, dataDisc, uEnd, 2*p+1, problemData.basesOnQuad);
  
  % Error in y-velocity (v)
  dataQ0T = (problemData.cDisc(:,:,3) * problemData.basesOnQuad.phi2D{qOrd}.') ./ (hDisc * problemData.basesOnQuad.phi2D{qOrd}.');
  dataDisc = projectDataQ0T2DataDisc(dataQ0T, 2*p, problemData.refElemPhiPhi, problemData.basesOnQuad);
  problemData.err(6) = computeL2Error(problemData.g, dataDisc, vEnd, 2*p+1, problemData.basesOnQuad);
end % if
end % function
