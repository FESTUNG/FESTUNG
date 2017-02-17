%
function ret = assembleVecEdgeMuFuncContValNeumann(g, markE0Tbdr, cDiscLamba, Nlambda, basesOnGamma )
K = g.numT;  
KEdge = g.numE;


% % Determine quadrature rule
% p = Nlambda - 1;
% qOrd = 2*p+1;  
% % [~, W] = quadRule1D(qOrd);
% [Q, W] = quadRule1D(qOrd);
% 
% % [Q1, Q2] = gammaMap(n, Q);
% % funcOnQuad = funcCont(Q2X1(Q1, Q2), Q2X2(Q1, Q2));
% % Kkn = markE0Tbdr(:, n) .* g.areaE0T(:,n);
% 
% % Assemble vector
% ret = zeros(KEdge, Nlambda);
% for iEl = 1 : K
%     % Determine mapping to physical element
%     Q2X1 = @(X1,X2) g.B(iEl,1,1)*X1 + g.B(iEl,1,2)*X2 + g.coordV0T(iEl,1,1)*ones(size(X1));
%     Q2X2 = @(X1,X2) g.B(iEl,2,1)*X1 + g.B(iEl,2,2)*X2 + g.coordV0T(iEl,1,2)*ones(size(X1));
%     
%     for n = 1 : 3
%         [Q1, Q2] = gammaMap(n, Q);
%         funcOnQuad = funcCont(Q2X1(Q1, Q2), Q2X2(Q1, Q2));
%         Kkn = markE0Tbdr(iEl, n) .* g.areaE0T(iEl,n);
%         for i = 1 : Nlambda
% %             W'
% %             basesOnGamma.phi1D{qOrd}(:,i)
% %             W' .* basesOnGamma.phi1D{qOrd}(:,i)
%             integral = funcOnQuad * ( W' .* basesOnGamma.phi1D{qOrd}(:,i) );
%             ret(g.E0T(iEl, n),i) = ret(g.E0T(iEl, n),i) + Kkn .* integral;
%         end % for
%     end
% end % for
% ret = reshape(ret',KEdge*Nlambda,1);
end

