
%
function dataDisc = projectFuncCont2FaceDataDisc(g, funcCont, ord, refFacePhiPhi, basesOnGamma)
validateattributes(funcCont, {'function_handle'}, {}, mfilename, 'funcCont');
validateattributes(basesOnGamma, {'struct'}, {}, mfilename, 'basesOnGamma');
ord = max(ord,1);  [Q1, W] = quadRule1D(ord);
N = size(refFacePhiPhi, 1);

rhs = zeros( g.numE, N );

% localIdx = g.E0E([1:g.numE],1);
% adjTri = g.T0E([1:g.numE]', 1);

% F1 = @(X1, X2) g.B(adjTri,1,1)*X1 + g.B(adjTri,1,2)*X2 + g.coordV0T(adjTri,1,1)*ones(size(X1));
% F2 = @(X1, X2) g.B(adjTri,2,1)*X1 + g.B(adjTri,2,2)*X2 + g.coordV0T(adjTri,1,2)*ones(size(X1));
% x1 = zeros(g.numE, N);
% x2 = zeros(g.numE, N);
% for i=1:g.numE
%     gammaMap( localIdx(i), Q1 );
%     [x1(i,:), x2(i,:)] = gammaMap( localIdx(i), Q1 );
% end
% 
% rhs = funcCont( F1(x1, x2), F2(x1, x2) ) * (repmat(W.', 1, N) .* basesOnGamma.phi1D{ord}(:,1:N));

%%TODO restructure the code such that we do not need the loop here
for i = 1:g.numE
    localIdx = g.E0E(i,1);
    adjTri = g.T0E(i, 1);

    F1 = @(X1, X2) g.B(adjTri,1,1)*X1 + g.B(adjTri,1,2)*X2 + g.coordV0T(adjTri,1,1)*ones(size(X1));
    F2 = @(X1, X2) g.B(adjTri,2,1)*X1 + g.B(adjTri,2,2)*X2 + g.coordV0T(adjTri,1,2)*ones(size(X1));

    [x1, x2] = gammaMap( localIdx, Q1 );

    rhs(i, :)  = funcCont( F1(x1, x2), F2(x1, x2) ) * (repmat(W.', 1, N) .* basesOnGamma.phi1D{ord}(:,1:N));
%     dataDisc(i, :) = rhs / refFacePhiPhi(:,:);
end % for
dataDisc = rhs / refFacePhiPhi(:,:);
end % function
