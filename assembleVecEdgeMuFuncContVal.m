%
function ret = assembleVecEdgeMuFuncContVal(g, markE0Tbdr, funcCont, Nlambda, basesOnGamma )
K = g.numT;  
KEdge = g.numE;


% Determine quadrature rule
p = Nlambda - 1;
qOrd = 2*p+1;  
% [~, W] = quadRule1D(qOrd);
[Q, W] = quadRule1D(qOrd);

% [Q1, Q2] = gammaMap(n, Q);
% funcOnQuad = funcCont(Q2X1(Q1, Q2), Q2X2(Q1, Q2));
% Kkn = markE0Tbdr(:, n) .* g.areaE0T(:,n);

% Assemble vector
ret = zeros(KEdge, Nlambda);
for iEl = 1 : K
    % Determine mapping to physical element
    Q2X1 = @(X1,X2) g.B(iEl,1,1)*X1 + g.B(iEl,1,2)*X2 + g.coordV0T(iEl,1,1)*ones(size(X1));
    Q2X2 = @(X1,X2) g.B(iEl,2,1)*X1 + g.B(iEl,2,2)*X2 + g.coordV0T(iEl,1,2)*ones(size(X1));
    
    for n = 1 : 3
        [Q1, Q2] = gammaMap(n, Q);
        funcOnQuad = funcCont(Q2X1(Q1, Q2), Q2X2(Q1, Q2));
        Kkn = markE0Tbdr(iEl, n) .* g.areaE0T(iEl,n);
        for i = 1 : Nlambda
%             W'
%             basesOnGamma.phi1D{qOrd}(:,i)
%             W' .* basesOnGamma.phi1D{qOrd}(:,i)
            integral = funcOnQuad * ( W' .* basesOnGamma.phi1D{qOrd}(:,i) );
            ret(g.E0T(iEl, n),i) = ret(g.E0T(iEl, n),i) + Kkn .* integral;
        end % for
    end
end % for
ret = reshape(ret',KEdge*Nlambda,1);
end

% % Check function arguments that are directly used
% validateattributes(markE0Tbdr, {'logical'}, {'size', [g.numT 3]}, mfilename, 'markE0Tbdr');
% validateattributes(funcCont, {'function_handle'}, {}, mfilename, 'funcCont');
% validateattributes(valOnQuad, {'numeric'}, {'size', [g.numT 3 length(W)]}, mfilename, 'valOnQuad');
% validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
% 
% if nargin > 6
%   ret = assembleVecEdgePhiIntFuncContVal_withAreaE0Tbdr(g, funcCont, valOnQuad, N, basesOnQuad, areaE0Tbdr, qOrd);
% else
%   ret = assembleVecEdgePhiIntFuncContVal_noAreaE0Tbdr(g, markE0Tbdr, funcCont, valOnQuad, N, basesOnQuad, qOrd);
% end % if
% 
% end % function
% %
% %===============================================================================
% %> @brief Helper function for the case that assembleVecEdgePhiIntFuncContVal()
% %> was called with a precomputed field areaE0Tbdr.
% %
% function ret = assembleVecEdgePhiIntFuncContVal_withAreaE0Tbdr(g, funcCont, valOnQuad, N, basesOnQuad, areaE0Tbdr, qOrd)
% % Determine quadrature rule
% [Q, W] = quadRule1D(qOrd);
% 
% % Determine mapping to physical element
% Q2X1 = @(X1,X2) g.B(:,1,1)*X1 + g.B(:,1,2)*X2 + g.coordV0T(:,1,1)*ones(size(X1));
% Q2X2 = @(X1,X2) g.B(:,2,1)*X1 + g.B(:,2,2)*X2 + g.coordV0T(:,1,2)*ones(size(X1));
% 
% % Assemble vector
% ret = zeros(g.numT, N);
% for n = 1 : 3
%   [Q1, Q2] = gammaMap(n, Q);
%   funcOnQuad = funcCont(Q2X1(Q1, Q2), Q2X2(Q1, Q2));
%   for i = 1 : N
%     integral = (funcOnQuad .* squeeze((valOnQuad(:, n, :) < 0) .* valOnQuad(:, n, :))) * ( W' .* basesOnQuad.phi1D{qOrd}(:,i,n));
%     ret(:,i) = ret(:,i) + areaE0Tbdr{n} .* integral;
%   end % for
% end % for
% 
% ret = reshape(ret',g.numT*N,1);
% end % function
% %
% %===============================================================================
% %> @brief Helper function for the case that assembleVecEdgePhiIntFuncContVal()
% %> was called with no precomputed field areaE0Tbdr.
% %
% function ret = assembleVecEdgePhiIntFuncContVal_noAreaE0Tbdr(g, markE0Tbdr, funcCont, valOnQuad, N, basesOnQuad, qOrd)
% % Determine quadrature rule
% [Q, W] = quadRule1D(qOrd);
% 
% % Determine mapping to physical element
% Q2X1 = @(X1,X2) g.B(:,1,1)*X1 + g.B(:,1,2)*X2 + g.coordV0T(:,1,1)*ones(size(X1));
% Q2X2 = @(X1,X2) g.B(:,2,1)*X1 + g.B(:,2,2)*X2 + g.coordV0T(:,1,2)*ones(size(X1));
% 
% % Assemble vector
% ret = zeros(g.numT, N);
% for n = 1 : 3
%   [Q1, Q2] = gammaMap(n, Q);
%   funcOnQuad = funcCont(Q2X1(Q1, Q2), Q2X2(Q1, Q2));
%   Kkn = markE0Tbdr(:, n) .* g.areaE0T(:,n);
%   for i = 1 : N
%     integral = (funcOnQuad .* squeeze((valOnQuad(:, n, :) < 0) .* valOnQuad(:, n, :))) * ( W' .* basesOnQuad.phi1D{qOrd}(:,i,n));
%     ret(:,i) = ret(:,i) + Kkn .* integral;
%   end % for
% end % for
% 
% ret = reshape(ret',g.numT*N,1);
% end % function
