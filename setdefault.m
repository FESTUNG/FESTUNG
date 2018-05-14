% Set a default value for a non-existing field in a struct.

%===============================================================================
%> @file
%>
%> @brief Set a default value for a non-existing field in a struct.
%===============================================================================
%>
%> @brief Set a default value for a non-existing field in a struct.
%>
%> This function expects a struct and a field name and checks if the
%> specified field exists. If not, it takes the given default value and
%> inserts the field into the struct with the given value.
%>
%> @param s        The struct to be checked and modified.
%> @param field    The name of the field to be checked
%> @param value    The default value to be used if the given field does not
%>                 exist.
%> @retval s       The struct, supplemented with the field.
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
function s = setdefault(s, field, value)
validateattributes(s, {'struct'}, {})
validateattributes(field, {'char'}, {'nonempty'})
if ~isfield(s, field)
  s.(field) = value;
end % if
end % function