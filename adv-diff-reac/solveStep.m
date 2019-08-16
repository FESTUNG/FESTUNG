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
sysU = reshape(problemData.uDisc', [], 1);

% System matrix from advection-reaction
sysA = -problemData.globAadv{1} - problemData.globAadv{2} + problemData.globAreac + problemData.globBadv;

if problemData.isIP
  sysA = sysA + problemData.globAIP - problemData.globBIP ...
          + problemData.symparam * problemData.globBsym + problemData.penparam * problemData.globBjmp;
  sysU = sysU + problemData.dt * (problemData.sysF - problemData.globM \ ( ...
          sysA * sysU + problemData.globJadv - problemData.symparam * problemData.globJsym + ...
          problemData.globJN - problemData.penparam * problemData.globJjmp ));
else
  % Compute Q from previous time level
  sysQ = cell(2,1);
  for m = 1 : 2
    sysAq = problemData.globAq{m} - problemData.globBq{m} - problemData.globBqN{m};
    sysQ{m} =  problemData.globM \ ( sysAq * sysU - problemData.globJD{m} );
  end % for
  
  % Compute U at next time level
  sysAu = sysA + problemData.penparam * problemData.globBjmp;
  sysAq = { -problemData.globAu{1} + problemData.globBu{1} + problemData.globBuD{1}, ...
            -problemData.globAu{2} + problemData.globBu{2} + problemData.globBuD{2} };
  sysU = sysU + problemData.dt * (problemData.sysF - problemData.globM \ ( ...
          sysAu * sysU + sysAq{1} * sysQ{1} + sysAq{2} * sysQ{2} + ...
          problemData.globJadv + problemData.globJN - problemData.penparam * problemData.globJjmp ));
end % if

problemData.uDisc = reshape(sysU, problemData.N, problemData.g.numT)';
end % function