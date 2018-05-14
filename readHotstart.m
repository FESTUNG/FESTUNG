% Reads data from an existing hotstart file.

%===============================================================================
%> @file
%>
%> @brief Reads data from an existing hotstart file.
%===============================================================================
%>
%> @brief Reads data from an existing hotstart file.
%>
%> This function allows to load hotstart data written by outputHotstart().
%> The stored data is returned as a struct.
%>
%> @param  filename   The name of the hotstart file.
%> @retval  data      A struct containing all data from the hotstart file.
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
function data = readHotstart(filename)
validateattributes(filename, {'char'}, {'nonempty'});
assert(exist(filename, 'file') == 2, 'Specified hotstart file does not exist');
data = struct();
load(filename);
end % function

