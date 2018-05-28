% Returns true if the current environment is GNU Octave.

%===============================================================================
%> @file
%>
%> @brief Returns true if the current environment is GNU Octave.
%===============================================================================
%>
%> @brief Returns true if the current environment is GNU Octave.
%>
%> This function uses the recommended way of distinguishing between MATLAB and
%> GNU Octave: 
%> https://octave.org/doc/interpreter/How-to-Distinguish-Between-Octave-and-Matlab.html
%> 
%> This is usually used when working around some bugs of Octave.
%> 
%> 
%> @retval  retval  true, if Octave is used, false otherwise (i.e., MATLAB)
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2018 Balthasar Reuter, Florian Frank, Vadym Aizinger
%> @author Balthasar Reuter, 2018
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
function retval = isOctave
  persistent cacheval;  % speeds up repeated calls

  if isempty (cacheval)
    cacheval = (exist ('OCTAVE_VERSION', 'builtin') > 0);
  end

  retval = cacheval;
end % function