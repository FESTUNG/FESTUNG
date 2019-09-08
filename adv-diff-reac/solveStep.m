% Second step of the four-part algorithm in the main loop.

%===============================================================================
%> @file
%>
%> @brief Second step of the four-part algorithm in the main loop.
%===============================================================================
%>
%> @brief Second step of the four-part algorithm in the main loop.
%>
%> The main loop repeatedly executes four steps until the parameter
%> <code>problemData.isFinished</code> becomes <code>true</code>.
%> These four steps are:
%>
%>  1. preprocessStep()
%>  2. solveStep()
%>  3. postprocessStep()
%>  4. outputStep()
%> 
%> This routine is executed second in each loop iteration and is intended to
%> produce the solution at the next step, e.g., at a new time-level.
%>
%> @param  problemData  A struct with problem parameters, precomputed
%>                      fields, and solution data structures (either filled
%>                      with initial data or the solution from the previous
%>                      loop iteration), as provided by configureProblem()  
%>                      and preprocessProblem(). @f$[\text{struct}]@f$
%> @param  nStep        The current iteration number of the main loop. 
%>
%> @retval problemData  The input struct enriched with solution data at
%>                      the next step. @f$[\text{struct}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Balthasar Reuter, Florian Frank, Vadym Aizinger
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
function problemData = solveStep(problemData, nStep) %#ok<INUSD>
K = problemData.g.numT;
N = problemData.N;
sysU = reshape(problemData.uDisc', [], 1);

% System matrix and right hand side vector
sysA = sparse(K * N, K * N);
sysV = problemData.globF;
sysW = problemData.globM;

% System matrix and right-hand side from advection-reaction
sysAadv = -problemData.globAadv{1} - problemData.globAadv{2} + problemData.globAreac ...
            + problemData.globBadv;
sysVadv = -problemData.globJadv;
if problemData.deltaAdvReac == 1
  sysA = sysAadv;
  sysV = sysV + sysVadv;
else
  sysV = sysV - sysA * sysU + sysVadv;
end % if

if problemData.isIP
  % System matrix and right-hand side from diffusion using IP discretization
  sysAdiff = problemData.globAIP - problemData.globBIP + problemData.symparam * problemData.globBsym ...
              + problemData.penparam * problemData.globBjmp;
  sysVdiff = problemData.symparam * problemData.globJsym - problemData.globJN ...
              + problemData.penparam * problemData.globJjmp;

  % Build linear system
  if problemData.deltaDiff == 1
    sysA = sysA + sysAdiff;
    sysV = sysV + sysVdiff;
  else
    sysV = sysV - sysAdiff * sysU + sysVdiff;
  end % if
else
  % System matrix and right-hand side from diffusion using IP discretization
  sysAdiffUU = problemData.penparam * problemData.globBjmp;
  sysAdiffUQ = { -problemData.globAu{1} + problemData.globBu{1} + problemData.globBuD{1}, ...
                 -problemData.globAu{2} + problemData.globBu{2} + problemData.globBuD{2} };
  sysAdiffQ = { -problemData.globAq{1} + problemData.globBq{1} + problemData.globBqN{1}; ...
                -problemData.globAq{2} + problemData.globBq{2} + problemData.globBqN{2} };
  sysVdiff = -problemData.globJN + problemData.penparam * problemData.globJjmp;

  % Build and solve linear system
  if problemData.deltaDiff == 1
    sysW = [ sparse(2 * K * N, 3 * K * N); ...
             sparse(K * N, 2 * K * N), problemData.globM ];
    sysA = [ problemData.globM    , sparse(K * N, K * N) , sysAdiffQ{1} ; ...
             sparse(K * N, K * N) , problemData.globM    , sysAdiffQ{2} ; ...
             sysAdiffUQ{1}        , sysAdiffUQ{2}        , sysAdiffUU ];
    sysV = [ -problemData.globJD{1} ; -problemData.globJD{2} ; sysV + sysVdiff ];
    sysU = [ zeros(2 * K * N, 1) ; sysU ];
  else
    sysQ = { problemData.globM \ (sysAdiffQ{1} * sysU - problemData.globJD{1}); ...
             problemData.globM \ (sysAdiffQ{2} * sysU - problemData.globJD{2}) };
    sysV = sysV - sysAdiffUU * sysU - sysAdiffUQ{1} * sysQ{1} - sysAdiffUQ{2} * sysQ{2} + sysVdiff;
  end % if
end % if

% Solve linear system
sysU = (sysW + problemData.dt * sysA) \ (sysW * sysU + problemData.dt * sysV);

% Extract U for implicit LDG
if ~problemData.isIP && problemData.deltaDiff == 1
  sysU = sysU(2 * K * N + 1 : 3 * K * N);
end % if

% % System matrix from advection-reaction
% sysA = -problemData.globAadv{1} - problemData.globAadv{2} + problemData.globAreac + problemData.globBadv;

% if problemData.isIP
%   sysA = sysA + problemData.globAIP - problemData.globBIP ...
%           + problemData.symparam * problemData.globBsym + problemData.penparam * problemData.globBjmp;
%   sysU = sysU + problemData.dt * (problemData.sysF - problemData.globM \ ( ...
%           sysA * sysU + problemData.globJadv - problemData.symparam * problemData.globJsym + ...
%           problemData.globJN - problemData.penparam * problemData.globJjmp ));
% else
%   % Compute Q from previous time level
%   sysQ = cell(2,1);
%   for m = 1 : 2
%     sysAq = problemData.globAq{m} - problemData.globBq{m} - problemData.globBqN{m};
%     sysQ{m} =  problemData.globM \ ( sysAq * sysU - problemData.globJD{m} );
%   end % for
  
%   % Compute U at next time level
%   sysAu = sysA + problemData.penparam * problemData.globBjmp;
%   sysAq = { -problemData.globAu{1} + problemData.globBu{1} + problemData.globBuD{1}, ...
%             -problemData.globAu{2} + problemData.globBu{2} + problemData.globBuD{2} };
%   sysU = sysU + problemData.dt * (problemData.sysF - problemData.globM \ ( ...
%           sysAu * sysU + sysAq{1} * sysQ{1} + sysAq{2} * sysQ{2} + ...
%           problemData.globJadv + problemData.globJN - problemData.penparam * problemData.globJjmp ));
% end % if

problemData.uDisc = reshape(sysU, problemData.N, K)';
end % function
