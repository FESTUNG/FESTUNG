function globM = assembleGlobM(g, hatM)
%% assembles global contributions of \int_T \varphi_i * \varphi_j
%% globM is the mass matrix
K = g.numT;
globM = 2 * kron(spdiags(g.areaT, 0, K, K), hatM);
end % function
