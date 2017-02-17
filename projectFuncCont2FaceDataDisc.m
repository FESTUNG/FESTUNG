
%
function dataDisc = projectFuncCont2FaceDataDisc(g, funcCont, ord, refFacePhiPhi, basesOnGamma)
validateattributes(funcCont, {'function_handle'}, {}, mfilename, 'funcCont');
validateattributes(basesOnGamma, {'struct'}, {}, mfilename, 'basesOnGamma');
ord = max(ord,1);  [Q1, W] = quadRule1D(ord);
N = size(refFacePhiPhi, 1);
K = g.numE;
rhs = zeros( g.numE, N );

rhsTest = zeros( K, N );
F1 = @(X1, X2) g.B(:,1,1)*X1 + g.B(:,1,2)*X2 + g.coordV0T(:,1,1)*ones(size(X1));
F2 = @(X1, X2) g.B(:,2,1)*X1 + g.B(:,2,2)*X2 + g.coordV0T(:,1,2)*ones(size(X1));
for n = 1:3
    [x1, x2] = gammaMap( n, Q1 );
    rhsTest(g.E0T(:, n),:) = W .* funcCont( F1(x1, x2), F2(x1, x2) ) * basesOnGamma.phi1D{ord}(:,:);    
end
dataDisc = rhsTest / refFacePhiPhi(:,:);
end % function
