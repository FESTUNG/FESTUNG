% problem = 'darcyVert';
% p = [0; 1; 2];
% testcase = 'convergence';
% numElem = [16, 16; 24, 24; 36, 36; 54, 54];%; 81, 81];
% numSteps = [100; 400; 1600; 6400];%; 25600];
% testcase = 'convergence2';
% numElem = 3 * 2.^(0:6)' * [1, 1];
% numSteps = ones(size(numElem,1),1);

% problem = 'sweVert';
% p = [0; 1; 2; 3];
% testcase = 'convergence';
% tEnd = 0.1;
% level = 0:4;
% isCoupling = false;
% numElem = 2.^level(:) * [2, 1];
% numSteps = { 10 * 2.^level(:); 40 * 2.^level(:); 160 * 2.^level(:); 640 * 2.^level(:) };

% problem = 'sweVert';
% p = [0; 1; 2];
% testcase = 'convergence';
% tEnd = 0.1;
% numElem = [16,16;24,24;36,36];
% numSteps = [100;400;1600];

problem = 'sweVert';
p = [0;1;2];
testcase = 'utbest_sinus';
tEnd = 86.4;
level = 0:3;
numElem = 2.^level(:) * [2, 1];
dt = { [1, 0.5, 0.5, 0.5] ; ...
        [0.5, 0.25, 0.25, 0.25]; ...
        [0.25, 0.125, 0.03125, 0.0025] ...
     };
numSteps = cellfun(@(c) ceil(tEnd ./ c), dt, 'UniformOutput', false);
isCoupling = false;



% problem = 'darcyVert_sweVert';
% % problem = 'darcyVert';
% p = [0; 1; 2];
% testcase = 'coupling';
% tEnd = 0.1;
% level = 0:5;
% numElem = 2.^level(:) * [2, 1];
% numSteps = { 1 * 2.^level(:); 4 * 4.^level(:); 16 * 8.^level(:) };

% problem = 'sweVert';
% p = [0; 1; 2];
% testcase = 'convergence';
% tEnd = 2;
% nLevel = 3;
% numElem = 2.^(0:nLevel-1)' * [2, 1];
% numSteps = { 2.^(0:nLevel-1)' * 5400; 2.^(0:nLevel-1)' * 10800; 2.^(0:nLevel-1)' * 43200 };
% numSteps = cellfun(@(c) max(tEnd/86400 * c, 1), numSteps, 'UniformOutput', false);
% testcase = 'convergence';
% tEnd = 0.1;
% level = 0:6;
% numElem = [16, 16; 24, 24; 36, 36];
% numSteps = [100; 400; 1600] / 10;
% numElem = 2.^level(:) * [2, 1];
% numSteps = { 10 * 2.^level(:); 40 * 2.^level(:); 160 * 2.^level(:) };
% testcase = 'convergence3';
% numElem = 3 * 2.^(0:6)' * [1, 1];
% numSteps = ones(size(numElem,1),1);
% numElem = [16, 16; 24, 24; 36, 36];
% numSteps = [100; 400; 1600];
% testcase = 'convergence3';
% numElem = [2, 1; 4, 2; 8, 4; 16, 8; 32, 16];
% numSteps= [8; 16; 32; 64; 128];
% testcase = 't-sin u';
% numElem = repmat([4, 4], 4, 1);
% numSteps = [10; 20; 40; 80];
% testcase = 't-sin h';
% numElem = repmat([4, 4], 4, 1);
% numSteps = [10; 20; 40; 80];
% testcase = 't-sin w';
% numElem = repmat([4, 4], 4, 1);
% numSteps = [10; 20; 40; 80];

if iscell(numSteps)
  for ip = 1 : size(p, 1)
    assert(isequal(size(numElem, 1), length(numSteps{ip})), 'numElem and numSteps must be same size')
  end
else
  assert(isequal(size(numElem, 1), length(numSteps)), 'numElem and numSteps must be same size')
end % if

err = {}; 
conv = {};
for ip = 1 : length(p)
  for i = 1 : size(numElem, 1)
    pd = struct;
    pd.isVisSol = false;
    pd.isVisGrid = false;
%     pd.isCouplingDarcy = isCoupling;
%     pd.isCouplingSWE = isCoupling;
    pd.testcase = testcase;
    pd.tEnd = tEnd;
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
      save([problem '_' testcase '_' num2str(isCoupling) '.mat'], 'err', 'conv');
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

% 
% for i = 1 : size(err_utbest{1}, 2)
%   fprintf('Err(%d)     C(%d)    ', i, i);
% end
% fprintf('\n');
% 
% for ip = 1 : length(err_utbest)
%   [N, n] = size(err_utbest{ip});
%   fprintf(repmat('-', 1, ceil(21 * n/2 - 8)));
%   fprintf(' p = %d ', p(ip)); 
%   fprintf([repmat('-', 1, ceil(21 * n/2 - 8)), '\n']);
%   for i = 1 : N
%     for j = 1 : n
%       fprintf('%6.2e  %6.3f   ', err_utbest{ip}(i,j), conv_utbest{ip}(i,j));
%     end % for j
%     fprintf('\n');
%   end % for i
% end % for p
% 
