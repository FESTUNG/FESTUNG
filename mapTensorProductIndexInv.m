function [m, n] = mapTensorProductIndexInv(i)
p = ceil(sqrt(i));
local_ind = p-(i-(p-1).^2);
m = p + min(local_ind, 0);
n = p - max(local_ind, 0);
end

