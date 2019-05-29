function ret = applyAlgebraicLimiter(g, dataV, eps, q)
numerator = zeros(g.numV,1);
denominator = zeros(g.numV,1);

% nominator(g.V0E(:,1)) = (dataV(g.V0E(:,1)) - dataV(g.V0E(:,2)));
% denominator(g.V0E(:,1)) = abs(dataV(g.V0E(:,1)) - dataV(g.V0E(:,2)));
% 
% nominator(g.V0E(:,2)) = nominator(g.V0E(:,2)) + (dataV(g.V0E(:,2)) - dataV(g.V0E(:,1)));
% denominator(g.V0E(:,2)) = denominator(g.V0E(:,2)) + abs(dataV(g.V0E(:,2)) - dataV(g.V0E(:,1)));

for i = 1 : g.numE
  v1 = g.V0E(i,1);  v2 = g.V0E(i,2);
  numerator(v1) = numerator(v1) + (dataV(v1) - dataV(v2));
  denominator(v1) = denominator(v1) + abs(dataV(v1) - dataV(v2));
  numerator(v2) = numerator(v2) + (dataV(v2) - dataV(v1));
  denominator(v2) = denominator(v2) + abs(dataV(v2) - dataV(v1));
end % for

numerator = abs(numerator) + eps;
denominator = denominator + eps;
alpha = 1 - (numerator ./ denominator).^q;
alphaMat = repmat(alpha, 1, g.numV);
ret = min(alphaMat, alphaMat');
ret = 1 - ret;
ret = ret - diag(ret);
end % function