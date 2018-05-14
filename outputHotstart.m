% Stores an arbitrary number of variables to a hotstart file with given
% filename.

%===============================================================================
%> @file
%>
%> @brief Stores an arbitrary number of variables to a hotstart file with 
%>        given filename.
%===============================================================================
%>
%> @brief Stores an arbitrary number of variables to a hotstart file with 
%>        given filename.
%>
%> It can be used to save the state of the algorithm and later continue
%> the computation.
%>
%> The stored data can be loaded using readHotstart().
%>
%> All variables and values that are to be stored must be given as pairs of
%> name and value.
%>
%> @par Example
%> @parblock
%> a = 5;
%> b = rand(10);
%> c = 1:9;
%>
%> @code
%> outputHotstart('hotstartfile.mat', 'a', a, 'b', b, 'c', c);
%> @endcode
%> @endparblock
%>
%> @param  filename     The name of the hotstart file.
%> @param  varargin     An arbitrary number of variables that is to be
%>                      stored. Must be specified as pairs of name and
%>                      value.
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
function outputHotstart(filename, varargin)
if nargin < 2
  error('No variable to save specified!')
end % if

%% Ensure target directory exists and file has an extension
validateattributes(filename, {'char'}, {'nonempty'})
[dirname, ~, ext] = fileparts(filename);
if ~isempty(dirname) && ~isdir(dirname)
  mkdir(dirname);
end % if
if isempty(ext)
  filename = [filename '.mat'];
end % if

%% Convert given data into struct
labels = varargin(1:2:end);
assert(length(labels) == (nargin-1) / 2, 'Not all variables given as pairs');
data = cell2struct(varargin(2:2:end), labels, 2);
validateattributes(data, {'struct'}, {'size', [1 1]})

%% Save struct
save(filename, 'data');
disp(['Hotstart data written to ' filename]);
end % function

