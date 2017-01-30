function ret = mapTensorProductIndex(m, n)
p = max(m,n);
ret = (p - 1) .* (p - 1) + p - m + n;
end

