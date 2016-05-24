function ret = readHotStart(g, file)
cd('output');
K = g.numT;
fileID = fopen(file);
ret = cell2mat(textscan(fileID, '%f'));
N = length(ret) / K;
assert(ismember(N, [1 3 6 10 15]), 'Hot start file has wrong number of values.');
ret = reshape(ret, N, K).';
fclose(fileID);
cd ..
end % function
