function ret = generateFlipArray(g)
K = g.numT;
ret = true( K, 3 );
% for iT = 1:K
%     for iE = 1:3
%         edgeNr = g.E0T(iT, iE);
%         if (g.T0E(edgeNr, 2) == iT)
%             ret(iT, iE) = 0;
%         end
%     end
% 
% end

for n=1:3
    ret(:, n) = g.T0E( g.E0T(:, n), 2) ~= (1:K)';
end

end % function