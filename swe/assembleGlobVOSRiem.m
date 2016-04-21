function globVOSRiem = assembleGlobVOSRiem(g, markE0TOS, HOS, PhiPhidiag, lambdaOSRiem, areaE0Tbdr)
%% assembles global open sea boundary edge contributions of \int_E \varphi_i * \lambda * \varphi_j
%% by building a matrix and 
%% assembles global open sea boundary edge contributions of \int_E \varphi_i * \lambda * HOS
%% by building a vector
%% here \lambda = abs( ( (H^0.5 * uH)^- + (H^0.5 * uH)^+ ) * nu^{1} + ( (H^0.5 * vH)^- + (H^0.5 * vH)^+ ) * nu^{2} ) / ((H^-)^1.5 + (H^+)^1.5)
%%              + ( 0.5 * gConst * ((H^-)^1.5 + (H^+)^1.5) / ((H^-)^0.5 + (H^+)^0.5) )^0.5
%% with H^- being the unknown height and H^+ the perscribed open sea boundary height
%% and uH^- = uH^+, vH^- = vH^+ respectively being the values of the unknown momenta from the interior
%% multiplying globVOSRiem{1} with H corresponds to open sea boundary edge integrals resulting from using a Lax-Friedrichs solver
%% globVOSRiem{2} corresponds to open sea boundary edge integrals resulting from using a Lax-Friedrichs solver
global gPhi1D
K = g.numT; N = size(PhiPhidiag, 1);
globVOSRiem = cell(2,1);
globVOSRiem{1} = sparse(K*N, K*N); 
globVOSRiem{2} = zeros(K, N);
p = (sqrt(8*N+1)-3)/2;  qOrd = max(2*p, 1); [~, W] = quadRule1D(qOrd);
for n = 1 : 3
  for r = 1 : length(W)
		if nargin > 5
			lambdaOSRiem{n}(:,r) = lambdaOSRiem{n}(:,r) .* areaE0Tbdr{n};
		else
			lambdaOSRiem{n}(:,r) = lambdaOSRiem{n}(:,r) .* g.areaE0T(:,n) .* markE0TOS(:,n);
		end % if
		globVOSRiem{1} = globVOSRiem{1} + 0.5 * kron( spdiags( lambdaOSRiem{n}(:,r), 0, K, K ), PhiPhidiag(:,:,n,r) );
  end % for
  globVOSRiem{2} = globVOSRiem{2} + ( HOS{n} .* lambdaOSRiem{n} .* ( ones(K,1) * W ) ) * gPhi1D{qOrd}(:, :, n);
end % for
globVOSRiem{2} = reshape(0.5 * globVOSRiem{2}.',K*N,1);
end % function
