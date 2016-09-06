function [ fLagr ] = projectDG2LagrangeSub( fDG )

[K, N] = size(fDG);
fLagr = zeros(K, N);

switch N
    case 1, L1 = 0.5;           L2 = 0.5;
    case 4, L1 = [0, 1, 1, 0];  L2 = [0, 0, 1, 1];
    case 9, L1 = [0, 1, 1, 0, 0.5, 1, 0.5, 0, 0.5];
            L2 = [0, 0, 1, 1, 0, 0.5, 1, 0.5, 0.5];
end  % switch

for i = 1 : N
    fLagr = fLagr + fDG(:, i) * phi(i, L1, L2);
end  % for

end  % function