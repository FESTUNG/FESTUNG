function ret = laguerrePolZeros(i, beta)
% % switch i
% %   case 1,  p = 1;%ret = ones(size(X));
% %   case 2,  p = [-beta 1];%ret = 1 - beta*X;
% %   case 3,  p = [1/2*beta^2 -2*beta 1];%ret = 1 - 2*beta*X + 1/2*(beta*X).^2;
% %   case 4,  p = [-1/6*beta^3 3/2*beta^2 -3*beta 1];%ret = 1 - 3*beta*X + 3/2*(beta*X).^2 - 1/6*(beta*X).^3;  
% %   case 5,  p = [1/24*beta^4 -2/3*beta^3 3*beta^2 -4*beta 1];%ret = 1 - 4*beta*X + 3*(beta*X).^2 - 2/3*(beta*X).^3 + 1/24*(beta*X).^4;
% %   case 6,  p = [-1/120*beta^5 5/24*beta^4 -5/3*beta^3 5*beta^2 -5*beta 1];%ret = 1 - 5*beta*X + 5*(beta*X).^2 - 5/3*(beta*X).^3 + 5/24*(beta*X).^4 - 1/120*(beta*X).^5;
% %   case 7,  p = [1/720*beta^6 -1/20*beta^5 5/8*beta^4 -10/3*beta^3 7.5*beta^2 -6*beta 1]; 
% %   otherwise
% %         ret = zeros(size(X));
% %         for k = (0:1:(i-1))
% %             ret = ret + ((-beta).^k)./factorial(k).*nchoosek((i-1),k).* X.^k;
% %         end
% % end % switch
% % ret = roots(p);

% switch i
%   case 1,  p = 1;
%     ret = roots(p);
%   case 2,  ret = [ 1];
%   case 3,  ret = [ 3.414213562373094; 0.585786437626904];
%   case 4,  ret = [ 6.289945082937476; 2.294280360279044; 0.415774556783478];  
%   case 5,  ret = [ 9.395070912301122; 4.536620296921130; 1.745761101158350; 0.322547689619392];
%   case 6,  ret = [12.640800844275793; 7.085810005858806; 3.596425771040726; 1.413403059106518; 0.263560319718140];
%   case 7,  ret = [15.982873980601580; 9.837467418382694; 5.775143569104473; 2.992736326059316; 1.188932101672622; 0.222846604179261];
%   otherwise
%     ret = zeros(i-1,1);
%     for k = (0:1:(i-1))
%       ret(i-k) = ((-beta).^k)./factorial(k).*nchoosek((i-1),k);
%     end
%     ret = roots(ret);
%     return;
% end % switch
% ret = 1/beta*ret;

alphaVec = 0:i-2;
alphaVec = (2*alphaVec+1);
betaVec  = 1:i-2;
A = diag(alphaVec) + diag(betaVec,1) + diag(betaVec,-1);
ret = eig(A)/beta;
end % function