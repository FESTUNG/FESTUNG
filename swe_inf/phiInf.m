%orthonormal in X2-Richtung
%orthogonal in X1-Richtung mit Faktor 1/beta ||phi||=1/beta
%phi(x) = p(x1)*p(x2)

function ret = phiInf(i, p, X1, X2, beta) % i \in {1,...,(p+1)*(pInf+1)}
leg = mod(i, p+1);
lag = ceil(i/(p+1));
if leg == 0
  leg = p+1;
end % if
ret = laguerrePol(lag, X1, beta) .* legendre01(leg, X2);
end % function

% function ret = phiInf(i, X1, X2, beta)
% switch i
%   case 1,  ret = laguerrePol(1, X1, beta).*ones(size(X2));
%   
%   case 2,  ret = laguerrePol(2, X1, beta).*ones(size(X2));
%   case 3,  ret = laguerrePol(1, X1, beta).*sqrt(12).*(X2-0.5);
%   case 4,  ret = laguerrePol(2, X1, beta).*sqrt(12).*(X2-0.5);
%     
%   case 5,  ret = laguerrePol(3, X1, beta).*ones(size(X2));
%   case 6,  ret = laguerrePol(3, X1, beta).*sqrt(12).*(X2-0.5);
%   case 7,  ret = laguerrePol(1, X1, beta).*sqrt(180).*(X2.^2-X2+1/6);
%   case 8,  ret = laguerrePol(2, X1, beta).*sqrt(180).*(X2.^2-X2+1/6);
%   case 9,  ret = laguerrePol(3, X1, beta).*sqrt(180).*(X2.^2-X2+1/6);  
%   
%   case 10,  ret = laguerrePol(4, X1, beta).*ones(size(X2));
%   case 11,  ret = laguerrePol(4, X1, beta).*sqrt(12).*(X2-0.5);
%   case 12,  ret = laguerrePol(4, X1, beta).*sqrt(180).*(X2.^2-X2+1/6);
%   case 13,  ret = laguerrePol(1, X1, beta).*sqrt(2800) ...
%       .*(X2.^3 - 3/2*X2.^2 + 3/5*X2 - 1/20);
%   case 14,  ret = laguerrePol(2, X1, beta).*sqrt(2800) ...
%       .*(X2.^3 - 3/2*X2.^2 + 3/5*X2 - 1/20);
%   case 15,  ret = laguerrePol(3, X1, beta).*sqrt(2800) ...
%       .*(X2.^3 - 3/2*X2.^2 + 3/5*X2 - 1/20);
%   case 16,  ret = laguerrePol(4, X1, beta).*sqrt(2800) ...
%       .*(X2.^3 - 3/2*X2.^2 + 3/5*X2 - 1/20);
%     
%   case 17,  ret = laguerrePol(5, X1, beta).*ones(size(X2));         
%   case 18,  ret = laguerrePol(5, X1, beta).*sqrt(12).*(X2-0.5);
%   case 19,  ret = laguerrePol(5, X1, beta).*sqrt(180).*(X2.^2-X2+1/6);
%   case 20,  ret = laguerrePol(5, X1, beta).*sqrt(2800) ...
%       .*(X2.^3 - 3/2*X2.^2 + 3/5*X2 - 1/20);  
%   case 21,  ret = laguerrePol(1, X1, beta).*sqrt(44100) ...
%       .*(X2.^4 - 2.*X2.^3 + 9/7*X2.^2 - 2/7*X2 + 1/70);
%   case 22,  ret = laguerrePol(2, X1, beta).*sqrt(44100) ...
%       .*(X2.^4 - 2.*X2.^3 + 9/7*X2.^2 - 2/7*X2 + 1/70);
%   case 23,  ret = laguerrePol(3, X1, beta).*sqrt(44100) ...
%       .*(X2.^4 - 2.*X2.^3 + 9/7*X2.^2 - 2/7*X2 + 1/70);
%   case 24,  ret = laguerrePol(4, X1, beta).*sqrt(44100) ...
%       .*(X2.^4 - 2.*X2.^3 + 9/7*X2.^2 - 2/7*X2 + 1/70);
%   case 25,  ret = laguerrePol(5, X1, beta).*sqrt(44100) ...
%       .*(X2.^4 - 2.*X2.^3 + 9/7*X2.^2 - 2/7*X2 + 1/70);
% end % switch
% end % function


