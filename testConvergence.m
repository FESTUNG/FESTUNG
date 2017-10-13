function [err, eoc] = testConvergence(problemName, varargin)
oldpath = addpath([pwd filesep 'core']);
fn_testConvergence = getFunctionHandle([problemName filesep 'testConvergence']);
[err, eoc] = fn_testConvergence(varargin{:});
path(oldpath);
end % function
