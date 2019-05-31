% Executes an arbitrary function outside the search path. This is slower
% than getFunctionHandle, hence the latter is preferred whenever possible.

%===============================================================================
%> @file
%>
%> @brief Executes an arbitrary function outside the search path. This is 
%>        slower than getFunctionHandle, hence the latter is preferred 
%>        whenever possible.
%===============================================================================
%>
%> @brief Executes an arbitrary function outside the search path. This is 
%>        slower than getFunctionHandle, hence the latter is preferred 
%>        whenever possible.
%>
%> Assuming a function file 'func.m' exists in folder 'path/to' relative to
%> the current working directory. Then, execin('path/to/func', arg1, arg2, ...)
%> calls this function with the specified arguments without having to
%> change the directory first or having to add the directory to the
%> path-variable.
%>
%> If this function is to be called more than once, then it is 
%> significantly faster (factor 1000!) to use getFunctionHandle() instead.
%>
%> Internally, execin() changes the working directory to the specified
%> path, calls the function and returns to the original working directory.
%> Thus, the called function can also call other functions located in its
%> directory, which is a difference to getFunctionHandle().
%>
%> This function is inspired by and a simplified version for our needs of
%> Jiro Doke's execin() function at Matlab File Exchange: 
%> https://www.mathworks.com/matlabcentral/fileexchange/8518-execin
%>
%> @par Example
%> @parblock
%> For a function 'func.m' in folder 'path/to', which expects two arguments,
%> the following is an example on how to use getFunctionHandle():
%>
%> @code
%> val = execin('path/to/func', 1, 2);
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
%> @param filename A string specifying the name including path of the
%>                 function to execute. It may contain the suffix '.m' or
%>                 leave it out.
%> @param varargin An arbitrary number of input parameters that are to be
%>                 passed to the called function.
%> @retval vararbout An arbitrary number of return values as given by the
%>                   called function.
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
function varargout = execin(filename, varargin)
% Extract path and name of function
validateattributes(filename, {'char'}, {'nonempty'})
[funcPath, funcName, funcExt] = fileparts(filename);
if isempty(funcExt)
  funcExt = '.m'; 
end % if
if exist(fullfile(funcPath, [funcName funcExt]), 'file') ~= 2
  error('Specified file ''%s'' does not exist', filename);
end % if

curDir = pwd;
try
  cd(funcPath);
  [varargout{1:nargout}] = feval(funcName, varargin{:});
  cd(curDir);
catch ME
  cd(curDir);
  fprintf('%s\n', ME.message);
  error('Could not execute function ''%s''.', filename);
end % try
end % function
