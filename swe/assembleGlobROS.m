function globROS = assembleGlobROS(g, markE0Tbdr, HOS, N, gConst, areaNuE0Tbdr)
%% assembles global open sea boundary edge contributions of -0.5 * gConst * \int_E \varphi_i * HOS^2 * nu
%% globROS corresponds to open sea boundary edge integrals in weak form of -0.5 * gConst * HOS^2 in PDE
global gPhi1D
K = g.numT; p = (sqrt(8*N+1)-3)/2;
qOrd = 2*p+1;  [~, W] = quadRule1D(qOrd); 
globROS = cell(2, 1);
globROS{1} = zeros(K, N); globROS{2} = zeros(K, N);
for n = 1 : 3
  HOSnWphi = bsxfun(@times, HOS{n} .* HOS{n}, W) * gPhi1D{qOrd}(:, :, n);
	if nargin > 5
	  globROS{1} = globROS{1} - bsxfun(@times, HOSnWphi, areaNuE0Tbdr{n,1});
		globROS{2} = globROS{2} - bsxfun(@times, HOSnWphi, areaNuE0Tbdr{n,2});
	else
	  globROS{1} = globROS{1} - bsxfun(@times, HOSnWphi, g.areaE0T(:,n) .* g.nuE0T(:,n,1) .* markE0Tbdr(:,n));
		globROS{2} = globROS{2} - bsxfun(@times, HOSnWphi, g.areaE0T(:,n) .* g.nuE0T(:,n,2) .* markE0Tbdr(:,n));
	end % if
end % for
globROS{1} = reshape(0.5 * gConst * globROS{1}.',K*N,1);  
globROS{2} = reshape(0.5 * gConst * globROS{2}.',K*N,1);
end % function
