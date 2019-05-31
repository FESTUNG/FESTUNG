% Basic routine to apply a slope limiter in any way for the depth-averaged
% shallow water equations.

%===============================================================================
%> @file
%>
%> @brief Basic routine to apply a slope limiter in any way for the depth-
%>        averaged shallow water equations.
%===============================================================================
%>
%> @brief Basic routine to apply a slope limiter in any way for the depth-
%>        averaged shallow water equations.
%>
%> The appropriate slope limiter for the 2D shallow water equations is
%> called depending on the value of slopeLimName.
%>
%> @param  pd					  The main struct of the code with problem parameters,
%>                      precomputed fields, and solution data structures, as 
%>											provided by initializeProblem() or postprocessSubStep().
%>                      @f$[\text{struct}]@f$
%> @param  slopeLimName A string, specifying the appropriate slope limiter.
%>											So far, 'segregated' and 'sequential' are supported.
%> @param  nSubStep     The current iteration number of the iterateSubSteps 
%>                      loop.
%>
%> @retval pd					  The input struct enriched with the solution data after
%>                      the slope limiting step. @f$[\text{struct}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2018 Florian Frank, Balthasar Reuter, Vadym Aizinger
%> 
%> @author Hennes Hajduk, 2018.
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
function pd = applySlopeLimiter(pd, slopeLimName, nSubStep)

time = getdefault(pd.tLvls, nSubStep+1, pd.t + pd.dt);

if strcmp(slopeLimName(end-4:end), '+norm')
  slopeLimName = slopeLimName(1:end-5);
  isNormLimiter = true;
else
  isNormLimiter = false;
end % if

switch slopeLimName
	case 'componentwise'
		pd = applyComponentwiseSlopeLimiter(pd, time);
  case 'vector'
    pd = applyVectorSlopeLimiter(pd, time);
	case 'sequentialSegregated'
		pd = applySequentialSegregatedSlopeLimiter(pd, time);
  case 'sequentialVector'
		pd = applySequentialVectorSlopeLimiter(pd, time);
	otherwise
		error('Unsupported slope limiter.');
end % switch

if isNormLimiter
  pd = applyNormSlopeLimiter(pd, time);
end % if

end % function

function pd = applyComponentwiseSlopeLimiter(pd, time)

for i = 1 : length(pd.slopeLimList)
  switch pd.slopeLimList{i}
    case 'height'
      if pd.isRivCont
        pd.xiV0Triv = pd.g.markV0TbdrRI .* computeFuncContV0T(pd.g, @(x1, x2) pd.xiRivCont(x1, x2, time));
      end % if
      if pd.isOSCont
        pd.xiV0Tos = pd.g.markV0TbdrOS .* computeFuncContV0T(pd.g, @(x1, x2) pd.xiOSCont(x1, x2, time));
      end % if
      xiV0T = pd.ramp(time/86400) * (pd.xiV0Triv + pd.xiV0Tos);
      pd.cDisc(:,:,1) = applySlopeLimiterDisc(pd.g, pd.cDisc(:,:,1), pd.g.markV0TbdrD, xiV0T, pd.globM, pd.globMDiscTaylor, pd.basesOnQuad, pd.typeSlopeLim);
    case 'momentum'
      if pd.isRivCont
        hV0Triv = pd.xiV0Triv - pd.zbV0T;
        pd.uHV0Triv = computeFuncContV0T(pd.g, @(x1, x2) pd.uRivCont(x1, x2, time)) .* hV0Triv;
        pd.vHV0Triv = computeFuncContV0T(pd.g, @(x1, x2) pd.vRivCont(x1, x2, time)) .* hV0Triv;
      end % if
      pd.cDisc(:,:,2) = applySlopeLimiterDisc(pd.g, pd.cDisc(:,:,2), pd.g.markV0TbdrRI, pd.ramp(time/86400) * pd.uHV0Triv, pd.globM, pd.globMDiscTaylor, pd.basesOnQuad, pd.typeSlopeLim);
      pd.cDisc(:,:,3) = applySlopeLimiterDisc(pd.g, pd.cDisc(:,:,3), pd.g.markV0TbdrRI, pd.ramp(time/86400) * pd.vHV0Triv, pd.globM, pd.globMDiscTaylor, pd.basesOnQuad, pd.typeSlopeLim);
    otherwise
      error('Componentwise slope limiting only for primary variables. For limiting the velocity use sequential limiting.')
  end % switch
end % for
end % function

function pd = applyVectorSlopeLimiter(pd, time)

for i = 1 : length(pd.slopeLimList)
  switch pd.slopeLimList{i}
    case 'height'
      if pd.isRivCont
        pd.xiV0Triv = pd.g.markV0TbdrRI .* computeFuncContV0T(pd.g, @(x1, x2) pd.xiRivCont(x1, x2, time));
      end % if
      if pd.isOSCont
        pd.xiV0Tos = pd.g.markV0TbdrOS .* computeFuncContV0T(pd.g, @(x1, x2) pd.xiOSCont(x1, x2, time));
      end % if
      xiV0T = pd.ramp(time/86400) * (pd.xiV0Triv + pd.xiV0Tos);
      pd.cDisc(:,:,1) = applySlopeLimiterDisc(pd.g, pd.cDisc(:,:,1), pd.g.markV0TbdrD, xiV0T, pd.globM, pd.globMDiscTaylor, pd.basesOnQuad, pd.typeSlopeLim);
    case 'momentum'
      if pd.isRivCont
        hV0Triv = pd.xiV0Triv - pd.zbV0T;
        pd.uHV0Triv = computeFuncContV0T(pd.g, @(x1, x2) pd.uRivCont(x1, x2, time)) .* hV0Triv;
        pd.vHV0Triv = computeFuncContV0T(pd.g, @(x1, x2) pd.vRivCont(x1, x2, time)) .* hV0Triv;
      end % if
      
      pd.cDisc = applyVectorSlopeLimiterDisc( pd.g, pd.cDisc, pd.zbDisc, pd.g.markV0TbdrRI, {pd.ramp(time/86400) * pd.uHV0Triv; pd.ramp(time/86400) * pd.vHV0Triv}, ...
                                              pd.globM, pd.globMDiscTaylor, pd.basesOnQuad, 'direction', pd.variant, pd.angle );
    otherwise
      error('Componentwise slope limiting only for primary variables. For limiting the velocity use sequential limiting.')
  end % switch
end % for
end % function

function pd = applySequentialSegregatedSlopeLimiter(pd, time)

if pd.isRivCont
  pd.xiV0Triv = pd.g.markV0TbdrRI .* computeFuncContV0T(pd.g, @(x1, x2) pd.xiRivCont(x1, x2, time));
  pd.uV0Triv = computeFuncContV0T(pd.g, @(x1, x2) pd.uRivCont(x1, x2, time));
  pd.vV0Triv = computeFuncContV0T(pd.g, @(x1, x2) pd.vRivCont(x1, x2, time)) ;
end % if
if pd.isOSCont
  pd.xiV0Tos = pd.g.markV0TbdrOS .* computeFuncContV0T(pd.g, @(x1, x2) pd.xiOSCont(x1, x2, time));
end % if

dataV0T = { pd.ramp(time/86400) * (pd.xiV0Triv + pd.xiV0Tos); pd.ramp(time/86400) * pd.uV0Triv; pd.ramp(time/86400) * pd.vV0Triv };
markV0Tbdr = { pd.g.markV0TbdrD, pd.g.markV0TbdrRI, pd.g.markV0TbdrRI };
pd.cDisc = applySequentialSlopeLimiterDisc(pd.g, pd.cDisc, pd.zbDisc, markV0Tbdr, dataV0T, pd.globM, pd.globMDiscTaylor, pd.basesOnQuad, pd.typeSlopeLim);

end % function

function pd = applySequentialVectorSlopeLimiter(pd, time)

if pd.isRivCont
  pd.xiV0Triv = pd.g.markV0TbdrRI .* computeFuncContV0T(pd.g, @(x1, x2) pd.xiRivCont(x1, x2, time));
  pd.uV0Triv = computeFuncContV0T(pd.g, @(x1, x2) pd.uRivCont(x1, x2, time));
  pd.vV0Triv = computeFuncContV0T(pd.g, @(x1, x2) pd.vRivCont(x1, x2, time)) ;
end % if
if pd.isOSCont
  pd.xiV0Tos = pd.g.markV0TbdrOS .* computeFuncContV0T(pd.g, @(x1, x2) pd.xiOSCont(x1, x2, time));
end % if

dataV0T = { pd.ramp(time/86400) * (pd.xiV0Triv + pd.xiV0Tos); pd.ramp(time/86400) * pd.uV0Triv; pd.ramp(time/86400) * pd.vV0Triv };
markV0Tbdr = {pd.g.markV0TbdrD, pd.g.markV0TbdrRI};
pd.cDisc = applySequentialVectorSlopeLimiterDisc(pd.g, pd.cDisc, pd.zbDisc, markV0Tbdr, dataV0T, pd.globM, pd.globMDiscTaylor, pd.basesOnQuad, pd.typeSlopeLim, pd.variant, pd.angle);

end % function

function pd = applyNormSlopeLimiter(pd, time) % requires the computation of boundary data and storing in pd by one of the other routines

for i = 1 : length(pd.slopeLimList)
  switch pd.slopeLimList{i}
    case 'momentum'
      pd.cDisc = applyVectorSlopeLimiterDisc( pd.g, pd.cDisc, pd.zbDisc, pd.g.markV0TbdrRI, {pd.ramp(time/86400) * pd.uHV0Triv; pd.ramp(time/86400) * pd.vHV0Triv}, ...
                                              pd.globM, pd.globMDiscTaylor, pd.basesOnQuad, 'norm-momentum', pd.variant, []);
    case 'velocity'
      pd.cDisc = applyVectorSlopeLimiterDisc( pd.g, pd.cDisc, pd.zbDisc, pd.g.markV0TbdrRI, {pd.ramp(time/86400) * pd.uV0Triv; pd.ramp(time/86400) * pd.vV0Triv }, ...
                                              pd.globM, pd.globMDiscTaylor, pd.basesOnQuad, 'norm-velocity', pd.variant, []);
    otherwise
      continue;
  end % switch
end % for
end % function
