function [err, eoc] = testConvergence(testcase, pLevel, hLevel, tLevel)
if nargin < 2
  pLevel = 0:4;
end % if
if nargin < 3
  hLevel = 1:6;
end % if
if nargin < 4
  tLevel = 1;
end % if

isSpatConv = length(hLevel) > 1 && length(tLevel) == 1;
isTimeSpatConv = length(hLevel) == length(tLevel);
assert(xor(isSpatConv, isTimeSpatConv), 'Time and space convergence levels must fit to each other!');
nLevel = max(length(hLevel), length(tLevel));

err = cell(size(pLevel));
eoc = cell(size(pLevel));

for ip = 1 : length(pLevel)
  for level = 1 : nLevel
    pd = struct('isVisSol', false, 'isVisGrid', false, 'isConvergence', true, 'testcase', testcase);
    pd.p = pLevel(ip);
    if isSpatConv || isTimeSpatConv
      pd.hmax = 2^-hLevel(level) / 3;
    end % if
    if ~isSpatConv
      pd.numSteps = 10 * 2^tLevel(level);
    end % if
    try
      pd = main('hdg_advection', pd);
      err{ip}(level) = pd.error;
      N = length(err{ip});
      if isSpatConv
        eoc{ip} = [ 0, log(err{ip}(1:N-1) ./ err{ip}(2:N)) ./ ...
                        ((hLevel(2:N) - hLevel(1:N-1)) * log(2)) ];
      else
        eoc{ip} = [ 0, log(err{ip}(1:N-1) ./ err{ip}(2:N)) ./ ...
                        ((tLevel(2:N) - tLevel(1:N-1)) * log(2)) ];
      end % if
      printConvergence(err, eoc, pLevel);
    catch ME
      stack = arrayfun(@(a) sprintf('  In %s (line %d)\n', a.file, a.line), ME.stack, 'UniformOutput', false);
      warning('%s: %s\n%s---', ME.identifier, ME.message, [stack{:}]);
      disp(ME.stack)
    end % try
  end % for level
end % for ip


end % function

function printConvergence(err, eoc, pLevel)
fprintf('Err        EOC\n');
fprintf('================\n');

for ip = 1 : length(err)
  N = length(err{ip});
  fprintf('---- p = %d -----\n', pLevel(ip)); 
  for i = 1 : N
    fprintf('%6.2e  %6.3f   \n', err{ip}(i), eoc{ip}(i));
  end % for i
end % for p
end % function