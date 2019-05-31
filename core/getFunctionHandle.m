% Returns a function handle for an arbitrary function file outside the
% search path.

%===============================================================================
%> @file
%>
%> @brief Returns a function handle for an arbitrary function file outside 
%>        the search path.
%===============================================================================
%>
%> @brief Returns a function handle for an arbitrary function file outside 
%>        the search path.
%>
%> Assuming a function file 'func.m' exists in folder 'path/to' relative to
%> the current working directory. Then, getFunctionHandle('path/to/func')
%> creates a function handle for this function file without adding the 
%> respective path to the path variable.
%>
%> This is helpful, e.g., if a function with the same name exists in the
%> current path and should not be shadowed.
%> 
%> Please note that the called function should not rely on any other
%> function files, as they might not be in the search path.
%>
%> If this function is to be called more than once, then this is 
%> significantly faster than changing the working directory every time for 
%> the function call or adding and removing the directory from the search
%> path.
%>
%> For functions that are called only once, execin() might be an
%> alternative.
%>
%> @par Example
%> @parblock
%> For a function 'func.m' in folder 'path/to', which expects two arguments,
%> the following is an example on how to use getFunctionHandle():
%>
%> @code
%> h = getFunctionHandle('path/to/func');
%> val = h(1,2);
%> @endcode
%> 
%> This gives the same result as the following snippet of code:
%>
%> @code
%> p = pwd;
%> cd('path/to');
%> val = func(1,2);
%> cd(p);
%> @endcode
%> @endparblock
%>
%> @par Speed comparison
%> @parblock
%> Using MATLAB R2016a on an Intel Core i7-4790 gives the following runtime
%> measurements for the different strategies:
%>
%> @code
%> n = 1e4;
%> a = rand(1);
%> b = rand(1);
%> c = rand(1);
%> dirname = 'test';
%> funcname = 'test';
%> curDir = pwd;
%> 
%> % Create test function
%> mkdir(dirname);
%> f = fopen([dirname filesep funcname '.m'], 'w');
%> fprintf(f, 'function ret = test(a,b,c)\nret=a+b+c;\nend\n');
%> fclose(f);
%> 
%> % (1) Reference: 0.004480 seconds
%> cd(dirname)
%> tic
%> for i = 1 : n
%>   d = test(a, b, c);
%> end
%> fprintf('Reference: %f\n', toc)
%> cd(curDir);
%> 
%> % (2) getFunctionHandle: 0.029120 seconds
%> tic
%> f = getFunctionHandle([dirname filesep funcname]);
%> for i = 1 : n
%>   d = f(a,b,c);
%> end
%> fprintf('getFunctionHandle: %f\n', toc)
%> 
%> % (3) addpath: 423.274278 seconds
%> tic
%> for i = 1 : n
%>   p = addpath(dirname);
%>   d = test(a,b,c);
%>   path(p);
%> end
%> fprintf('addpath: %f\n', toc)
%> 
%> % (4) cd: 19.549983 seconds
%> tic
%> for i = 1 : n
%>   cd(dirname);
%>   d = test(a,b,c);
%>   cd(curDir);
%> end
%> fprintf('cd: %f\n', toc)
%> @endcode
%> @endparblock
%>
%> @param filename A string specifying the name including path of the
%>                 function to execute. It may contain the suffix '.m' or
%>                 leave it out.
%> @retval h       A function handle for the given function
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
function h = getFunctionHandle(filename)
% Extract path and name of function
validateattributes(filename, {'char'}, {'nonempty'})
[funcPath, funcName, funcExt] = fileparts(filename);
if isempty(funcExt)
  funcExt = '.m'; 
end % if
if exist(fullfile(funcPath, [funcName funcExt]), 'file') ~= 2
  error('Specified file ''%s'' does not exist', filename);
end % if

% Create function handle
% Originally, this was done by changing the directory, creating the handle
% and changing back to the original directory. However, for unknown reasons
% this leads to Octave not finding the correct functions. 
% Adding the directory of the desired function to the search path, creating
% the handle and restoring the original path is a lot slower (10x on
% Octave, 150x on MATLAB) but works on both.
oldpath = addpath(funcPath);
% curDir = pwd;
try
%   cd(funcPath);
  h = str2func(funcName);
  path(oldpath);
%   cd(curDir);
catch ME
  path(oldpath);
%   cd(curDir);
  fprintf('%s\n', ME.message);
  error('Could not create function handle for ''%s''.', filename);
end % try
end % function
