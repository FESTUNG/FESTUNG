% TODO
function invA = blkinv(A, blockSize)
validateattributes(A, {'numeric'}, {'square'}, 'A')
validateattributes(blockSize, {'numeric'}, {'scalar', '>', 0}, 'blockSize')

n = size(A,1);
numBlocks = ceil(n / blockSize);
idxEnd = 0;
    
invA = cell(numBlocks, 1);
for idxBlock = 1 : numBlocks - 1
  idxStart = (idxBlock - 1) * blockSize + 1;
  idxEnd = idxBlock * blockSize;
  invA{idxBlock} = A(idxStart : idxEnd, idxStart : idxEnd) \ speye(blockSize);
end % for
invA{numBlocks} = A(idxEnd + 1 : end, idxEnd + 1 : end) \ speye(n - (numBlocks - 1) * blockSize);
invA = blkdiag(invA{:});
end % function