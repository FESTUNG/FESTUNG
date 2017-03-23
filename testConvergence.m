% problem = 'darcyVert';
% p = [0; 1; 2];
% testcase = 'convergence';
% numElem = [16, 16; 24, 24; 36, 36; 54, 54; 81, 81];
% numSteps = [100; 400; 1600; 6400; 25600];
% testcase = 'convergence2';
% numElem = 3 * 2.^(0:6)' * [1, 1];
% numSteps = ones(size(numElem,1),1);

problem = 'sweVert';
p = [0; 1; 2];
% testcase = 'convergence';
% numElem = [2, 1; 4, 2; 8, 4; 16, 8; 32, 16];
% numSteps = [50; 100; 200; 400; 800];
% testcase = 'convergence2';
% numElem = [16, 16; 24, 24; 36, 36];
% numSteps = [100; 400; 1600];
testcase = 'convergence3';
numElem = 3 * 2.^(0:6)' * [1, 1];
numSteps = ones(size(numElem,1),1);
% testcase = 'sin h u';
% numElem = [2, 1; 4, 2; 8, 4; 16, 8; 32, 16];
% numSteps= [8; 16; 32; 64; 128];

if iscell(numSteps)
  for ip = p
    assert(isequal(size(numElem, 1), length(numSteps{ip+1})), 'numElem and numSteps must be same size')
  end
else
  assert(isequal(size(numElem, 1), length(numSteps)), 'numElem and numSteps must be same size')
end % if

err = {}; 
conv = {};
for ip = 1 : size(p, 1)
  for i = 1 : size(numElem, 1)
    pd = struct;
    pd.testcase = testcase;
    pd.p = p(ip);
    pd.numElem = numElem(i, :);
    if iscell(numSteps)
      pd.numSteps = numSteps{ip}(i);
    else
      pd.numSteps = numSteps(i);
    end % if
    try
      pd = main(problem, pd); 
      err{ip}(i,:) = pd.error;  %#ok<SAGROW>
      [N, n] = size(err{ip});
      conv{ip} = [zeros(1,n); ...
                  log(err{ip}(1:N-1,:) ./ err{ip}(2:N,:)) ./ ...
                    repmat(log(numElem(2:N,1) ./ numElem(1:N-1,1)), 1, n)];  %#ok<SAGROW>
      disp(conv{ip}); 
    catch e
      warning('%s: %s', e.identifier, e.message);
    end % try
  end % for i
end % for p

for i = 1 : size(err{1}, 2)
  fprintf('Err(%d)     C(%d)    ', i, i);
end
fprintf('\n');

for ip = 1 : length(err)
  [N, n] = size(err{ip});
  fprintf(repmat('-', 1, ceil(21 * n/2 - 8)));
  fprintf(' p = %d ', p(ip)); 
  fprintf([repmat('-', 1, ceil(21 * n/2 - 8)), '\n']);
  for i = 1 : N
    for j = 1 : n
      fprintf('%6.2e  %6.3f   ', err{ip}(i,j), conv{ip}(i,j));
    end % for j
    fprintf('\n');
  end % for i
end % for p
