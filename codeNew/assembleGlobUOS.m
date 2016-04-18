function globUOS = assembleGlobUOS(g, markE0Tbdr, PhiPhidiag, HOS, gCdiag, areaE0Tbdr)
%% assembles global open sea boundary edge contributions of \int_E \varphi_i * ( uH * nu^{1} + vH * nu^{2} ) / HOS * \varphi_j
%% multiplying globUOS with uH (in the first momentum equation), vH (in the second momentum equation)
%% corresponds to open sea boundary edge integrals in weak form of div( [uH*uH, vH*uH; uH*vH, vH*vH] / HOS ) in PDE
%% here the values of the unknown momenta from the interior are used
K = g.numT; N = size(PhiPhidiag, 1); 
globUOS = sparse(K*N, K*N);
p = (sqrt(8*N+1)-3)/2; qOrd = 2*p+1;  [~, W] = quadRule1D(qOrd);
for n = 1 : 3
	for r = 1 : length(W)
		if nargin > 5
	    globUOS = globUOS + kron( spdiags( areaE0Tbdr{n} .* setNaN2Zero( ( gCdiag{2,n}(:,r) .* g.nuE0T(:,n,1) + gCdiag{3,n}(:,r) .* g.nuE0T(:,n,2) ) ...
																																			./ HOS{n}(:,r) ), 0, K, K ), PhiPhidiag(:,:,n,r) );
		else
	    globUOS = globUOS + kron( spdiags( g.areaE0T(:,n) .* markE0Tbdr(:,n) .* setNaN2Zero( ( g.nuE0T(:,n,1) .* gCdiag{2,n}(:,r) + ... 
																																														 g.nuE0T(:,n,2) .* gCdiag{3,n}(:,r) ) ...
																																														./ HOS{n}(:,r) ), 0, K, K ), PhiPhidiag(:,:,n,r) );
		end % if
  end % for
end % for
end % function
