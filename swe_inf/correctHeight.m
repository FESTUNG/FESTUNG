function cDG = correctHeight(p, cDG, minTol, corrSys, nStep)
[K, N] = size(cDG(:,:,1));
if p == 0
  HLagrange = cDG(:,:,1) * phi(1,1/3,1/3);
else
  HLagrange = zeros(K,3);
  for i = 1:N
    HLagrange = HLagrange + cDG(:,i,1) * phi(i, [0, 1, 0], [0, 0, 1]);
  end % for
end % if
corr = max(minTol - HLagrange, 0);
if p == 0
  if ~isequal(corr, zeros(K,1))
    warning([ 'Correction of varaible c1 in vertices is needed at timestep' ' ' num2str(nStep) '. Norm of correction values:' ' ' num2str(norm(corr)) ] );
    cDG(:,:,1) = cDG(:,:,1) + corr / phi(1,1/3,1/3);
  end %  if
else
  if ~isequal(corr, zeros(K,3))
    warning([ 'Correction of varaible c1 in vertices is needed at timestep' ' ' num2str(nStep) '. Norm of correction values:' ' ' num2str(norm(corr)) ] );
    cDG(:,1:3,1) = cDG(:,1:3,1) + corr / corrSys;
  end % if
end % if
end % function
