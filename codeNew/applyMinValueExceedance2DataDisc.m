function dataDisc = applyMinValueExceedance2DataDisc(g, dataDisc, corrSys, nStep, minTol, correctionTolerance)
if nargin < 4
	nStep = [];
end % if
if nargin < 5
	minTol = 0;
end % if
if nargin < 6
	correctionTolerance = Inf;
end % if
[K, N] = size(dataDisc);
p = (sqrt(8*N+1)-3)/2;
dataV0T = zeros(K,3);
for i = 1:N
  dataV0T = dataV0T + dataDisc(:,i) * phi(i, [0, 1, 0], [0, 0, 1]);
end % for
corr = max(minTol - dataV0T, 0);
maxCorr = max(corr(:));
if maxCorr > correctionTolerance
  [indx, indy] = find(corr == maxCorr, 1);
  error([ 'Unknown at node ' num2str(g.V0T(indx, indy)) ' and step ' num2str(nStep) ' is ' num2str(dataV0T(indx, indy)) ...
          ', which is below the tolerated minimal value.' ]);
end % if
if p == 0
  if ~isequal(corr, zeros(K,3))
    [indx, indy] = find(corr == maxCorr, 1);
    dataDisc = dataDisc + max(corr, [], 2) / phi(1, 1/3, 1/3);
    warning([ 'A maximal value of ' num2str(maxCorr) ' had to be added to unknown at node ' num2str(g.V0T(indx, indy)) ...
              ' and step ' num2str(nStep) '.' ]);
  end %  if
else
  if ~isequal(corr, zeros(K,3))
    [indx, indy] = find(corr == maxCorr, 1);
    dataDisc(:,1:3) = dataDisc(:,1:3) + corr / corrSys; % TODO evtl speichern
    warning([ 'A maximal value of ' num2str(maxCorr) ' had to be added to unknown at node ' num2str(g.V0T(indx, indy)) ...
              ' and step ' num2str(nStep) '.' ]);
  end % if
end % if
end % function
