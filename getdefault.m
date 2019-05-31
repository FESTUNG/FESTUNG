% Returns the entry of a given datastructure associated with the given
% index if it exists, or a default value otherwise.

%===============================================================================
%> @file
%>
%> @brief Returns the entry of a given datastructure associated with the 
%>        given index if it exists, or a default value otherwise.
%===============================================================================
%>
%> @brief Returns the entry of a given datastructure associated with the 
%>        given index if it exists, or a default value otherwise.
%>
%> This is tested to work with vectors, matrices, cells, and struct.
%>
%> The index can be given as a linear or multi-dimensional index (for
%> matrices, cell arrays) or as a string (for fields in a struct).
%>
%> @par Example
%> @code
%> v = [0, 1, 2];
%> for i = 1:4
%>   fprintf('%d ', getdefault(v, i, -1))
%> end
%> A = rand(5);
%> fprintf('%d\n', getdefault(A, [7 3], 100))
%> @endcode
%> produces the output '0 1 2 -1 100'.
%> 
%> @param  v       The datastructure in which to search for the associated entry.
%> @param  idx     The index in the datastructure v.
%> @param  default The default value to return if the associated entry does not exist.
%>
%> This file is part of FESTUNG
%>
%> @copyright 2016 Balthasar Reuter, Florian Frank, Vadym Aizinger
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
function ret = getdefault(v, idx, default)
switch class(v)
  
  case { 'single', 'double', 'char' }
    if isValidIdx(v, idx)
      ret = v(linearizeIdx(idx, size(v)));
    else
      ret = default;
    end % if
    
  case 'cell'
    if isValidIdx(v, idx)
      ret = v{linearizeIdx(idx, size(v))};
    else
      ret = default;
    end % if
    
  case 'struct'
    if isfield(v, idx)
      ret = v.(idx);
    else
      ret = default;
    end % if
    
  otherwise
    error('Unsupported class type for v')
    
end % switch
end % function
%
% Helper function to test if a valid index
function isValid = isValidIdx(v, idx)
isValid = (isscalar(idx) && idx <= numel(v)) || ...
          (numel(idx) == numel(size(v)) && all(idx <= size(v)));
end % function
%
% Helper function to get a linearized index
function idx = linearizeIdx(idx, n)
idx = cellfun(@(c) c(1), num2cell(idx), 'UniformOutput', false);
idx = sub2ind(n, idx{:});
end % function

