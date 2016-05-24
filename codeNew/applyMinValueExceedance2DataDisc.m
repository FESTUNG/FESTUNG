function dataDisc = applyMinValueExceedance2DataDisc(g, dataDisc, corrSys, nStep, minTol, correctionTolerance)
if nargin < 3
	nStep = [];
end % if
if nargin < 4
	minTol = 0;
end % if
if nargin < 5
	correctionTolerance = Inf;
end % if
[K, N] = size(dataDisc);
p = (sqrt(8*N+1)-3)/2;
if p == 0
  dataLagrange = dataDisc * phi(1,1/3,1/3);
else
  dataLagrange = zeros(K,3);
  for i = 1:N
    dataLagrange = dataLagrange + dataDisc(:,i) * phi(i, [0, 1, 0], [0, 0, 1]);
  end % for
end % if
corr = max(minTol - dataLagrange, 0);
maxCorr = max(abs(corr(:)));
[indx, indy] = find(abs(corr) == maxCorr);
if maxCorr > correctionTolerance
  error([ 'Unknown at node ' num2str(g.V0T(indx, indy)) ' and step ' num2str(nStep) ' is ' num2str(abs(dataLagrange(indx, indy))) ...
          ', which is below the tolerated minimal value.' ]);
end % if
if p == 0
  if ~isequal(corr, zeros(K,1))
    dataDisc = dataDisc + corr / phi(1,1/3,1/3);
    warning([ 'A maximal value of ' num2str(maxCorr) ' had to be added to unknown at node ' num2str(g.V0T(indx, indy)) ...
              ' and step ' num2str(nStep) '.' ]);
  end %  if
else
  if ~isequal(corr, zeros(K,3))
    dataDisc(:,1:3) = dataDisc(:,1:3) + corr / corrSys;
    warning([ 'A maximal value of ' num2str(maxCorr) ' had to be added to unknown at node ' num2str(g.V0T(indx, indy)) ...
              ' and step ' num2str(nStep) '.' ]);
  end % if
end % if
end % function
