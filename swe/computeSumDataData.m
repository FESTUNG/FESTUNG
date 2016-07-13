function ret = computeSumDataData(dataDisc1, dataDisc2)
[K, N1] = size(dataDisc1);
[~, N2] = size(dataDisc2);

validateattributes(dataDisc1, {'numeric'}, {'size', [K NaN]}, mfilename, 'dataDisc1')
validateattributes(dataDisc2, {'numeric'}, {'size', [K NaN]}, mfilename, 'dataDisc2')

ret = [dataDisc1, zeros(K,N2-N1)] + [dataDisc2, zeros(K,N1-N2)];
end % function