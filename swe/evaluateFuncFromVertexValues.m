function cEval = evaluateFuncFromVertexValues(g, vertexValues, X1, X2)
%% evaluates a piecewise linear continuous function in a set of matching physical points X1, X2, where the number of rows in X1, X2 is K 
%% and the points (X1(k,r), X2(k,r)) are in element k
%% Note that due to the backtransformation to the reference element and repeated use of piecewise multiplications with bsxfun this method can be quite
%% time consuming and is therefore meant to be used only for preprocessing purposes contrary to using it in every time step.
lenX = size(X1,2);
vertexValues = vertexValues(:);
FinvX1 = 0.5 ./ (g.areaT * ones(1,lenX)) .* ( bsxfun(@times, g.B(:,2,2), X1) - bsxfun(@times, g.B(:,1,2), X2) - (g.B(:,2,2).*g.coordV0T(:,1,1) - g.B(:,1,2).*g.coordV0T(:,1,2) ) * ones(1,lenX) );
FinvX2 = 0.5 ./ (g.areaT * ones(1,lenX)) .* (-bsxfun(@times, g.B(:,2,1), X1) + bsxfun(@times, g.B(:,1,1), X2) + (g.B(:,2,1).*g.coordV0T(:,1,1) - g.B(:,1,1).*g.coordV0T(:,1,2) ) * ones(1,lenX) );
cEval  = bsxfun(@times, vertexValues(g.V0T(:,1)), 1 - FinvX1 - FinvX2) ...
       + bsxfun(@times, vertexValues(g.V0T(:,2)),     FinvX1         ) ...
       + bsxfun(@times, vertexValues(g.V0T(:,3)),             FinvX2 );
end % function
