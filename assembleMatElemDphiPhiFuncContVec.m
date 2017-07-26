% TODO
function ret = assembleMatElemDphiPhiFuncContVec(g, refElemDphiPhiPerQuad, funcCont1, funcCont2, qOrd)
K = g.numT;
[N, ~, R] = size(refElemDphiPhiPerQuad{1});
validateattributes(refElemDphiPhiPerQuad, {'cell'}, {'size', [2 1]}, mfilename, 'refElemDphiPhiPerQuad');
validateattributes(refElemDphiPhiPerQuad{1}, {'numeric'}, {'size', [N N R]}, mfilename, 'refElemDphiPhiPerQuad{1}');
validateattributes(refElemDphiPhiPerQuad{2}, {'numeric'}, {'size', [N N R]}, mfilename, 'refElemDphiPhiPerQuad{2}');
validateattributes(funcCont1, {'function_handle'}, {'size', [1 1]}, mfilename, 'funcCont1');
validateattributes(funcCont2, {'function_handle'}, {'size', [1 1]}, mfilename, 'funcCont2');

if nargin < 5, p = (sqrt(8*N+1)-3)/2;  qOrd = max(2*p, 1); end
[Q1, Q2, ~] = quadRule2D(qOrd);

% Assemble matrix
ret = { zeros(K*N, N), zeros(K*N,N) };
for r = 1 : R
  valOnQuad1 = funcCont1(g.mapRef2Phy(1, Q1(r), Q2(r)), g.mapRef2Phy(2, Q1(r), Q2(r)));
  ret{1} = ret{1} + kron(g.B(:,2,2) .* valOnQuad1, refElemDphiPhiPerQuad{1}(:, :, r)) ...
                  - kron(g.B(:,2,1) .* valOnQuad1, refElemDphiPhiPerQuad{2}(:, :, r));
        
  valOnQuad2 = funcCont2(g.mapRef2Phy(1, Q1(r), Q2(r)), g.mapRef2Phy(2, Q1(r), Q2(r)));
  ret{2} = ret{2} - kron(g.B(:,1,2) .* valOnQuad2, refElemDphiPhiPerQuad{1}(:, :, r)) ...
                  + kron(g.B(:,1,1) .* valOnQuad2, refElemDphiPhiPerQuad{2}(:, :, r));
end %for
ret{1} = kronVec(speye(K,K), ret{1});
ret{2} = kronVec(speye(K,K), ret{2});

% ret = { sparse(K*N, K*N); sparse(K*N, K*N) };
% for r = 1 : R
%   valOnQuad1 = funcCont1(g.mapRef2Phy(1, Q1(r), Q2(r)), g.mapRef2Phy(2, Q1(r), Q2(r)));
%   ret{1} = ret{1} ...
%           + kron(spdiags(g.B(:,2,2) .* valOnQuad1, 0, K, K), refElemDphiPhiPerQuad{1}(:, :, r) ) ...
%           - kron(spdiags(g.B(:,2,1) .* valOnQuad1, 0, K, K), refElemDphiPhiPerQuad{2}(:, :, r) );
%         
%   valOnQuad2 = funcCont2(g.mapRef2Phy(1, Q1(r), Q2(r)), g.mapRef2Phy(2, Q1(r), Q2(r)));
%   ret{2} = ret{2} ...
%           - kron(spdiags(g.B(:,1,2) .* valOnQuad2, 0, K, K), refElemDphiPhiPerQuad{1}(:, :, r) ) ...
%           + kron(spdiags(g.B(:,1,1) .* valOnQuad2, 0, K, K), refElemDphiPhiPerQuad{2}(:, :, r) );
% end %for
end
