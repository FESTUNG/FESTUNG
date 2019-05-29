
function ret = gradPhiInf(i, p, m, X1, X2, beta)
leg = mod(i, p+1);
lag = ceil(i/(p+1));
if leg == 0
  leg = p+1;
end % if
ret = laguerrePol(lag, X1, beta) .* legendre01(leg, X2);

switch m
  case 1, ret = laguerrePolPrime(lag, X1, beta) .* legendre01(leg, X2);
  case 2, ret = laguerrePol(lag, X1, beta) .* legendre01Prime(leg, X2);
end %switch

% switch m
%   case 1
%     switch i
%       case 1,  ret = zeros(size(X2));
%         
%         
%       case 2,  ret = -beta.*ones(size(X2));
%       case 3,  ret = zeros(size(X2));
%       case 4,  ret = -beta.*sqrt(12).*(X2-0.5);      
%         
%       case 5,  ret = (-2*beta+beta.^2*X1).*ones(size(X2));
%       case 6,  ret = (-2*beta+beta.^2*X1).*sqrt(12).*(X2-0.5);        
%       case 7,  ret = zeros(size(X2));
%       case 8,  ret = -beta.*sqrt(180).*(X2.^2-X2+1/6);
%       case 9,  ret = (-2*beta+beta.^2*X1).*sqrt(180).*(X2.^2-X2+1/6);                
%         
%       case 10,  ret = (-3*beta+3*beta^2*X1-0.5*beta^3*X1.^2).*ones(size(X2));
%       case 11,  ret = (-3*beta+3*beta^2*X1-0.5*beta^3*X1.^2).*sqrt(12).*(X2-0.5);   
%       case 12,  ret = (-3*beta+3*beta^2*X1-0.5*beta^3*X1.^2).*sqrt(180).*(X2.^2-X2+1/6);    
%       case 13,  ret = zeros(size(X2));
%       case 14,  ret = -beta.*sqrt(2800).*(X2.^3 - 3/2*X2.^2 + 3/5*X2 - 1/20);
%       case 15,  ret = (-2*beta+beta.^2*X1).*sqrt(2800).*(X2.^3 - 3/2*X2.^2 + 3/5*X2 - 1/20);
%       case 16,  ret = (-3*beta+3*beta^2*X1-0.5*beta^3*X1.^2).*sqrt(2800).*(X2.^3 - 3/2*X2.^2 + 3/5*X2 - 1/20);               
%         
%       case 17,  ret = (-4*beta+6*beta^2*X1-2*beta^3*X1.^2+1/6*beta^4*X1.^3).*ones(size(X2));         
%       case 18,  ret = (-4*beta+6*beta^2*X1-2*beta^3*X1.^2+1/6*beta^4*X1.^3).*sqrt(12).*(X2-0.5);
%       case 19,  ret = (-4*beta+6*beta^2*X1-2*beta^3*X1.^2+1/6*beta^4*X1.^3).*sqrt(180).*(X2.^2-X2+1/6);
%       case 20,  ret = (-4*beta+6*beta^2*X1-2*beta^3*X1.^2+1/6*beta^4*X1.^3).*sqrt(2800).*(X2.^3 - 3/2*X2.^2 + 3/5*X2 - 1/20);
%       case 21,  ret = zeros(size(X2));
%       case 22,  ret = -beta.*sqrt(44100).*(X2.^4 - 2.*X2.^3 + 9/7*X2.^2 - 2/7*X2 + 1/70);
%       case 23,  ret = (-2*beta+beta.^2*X1).*sqrt(44100) ...
%           .*(X2.^4 - 2.*X2.^3 + 9/7*X2.^2 - 2/7*X2 + 1/70);
%       case 24,  ret = (-3*beta+3*beta^2*X1-0.5*beta^3*X1.^2).*sqrt(44100) ...
%           .*(X2.^4 - 2.*X2.^3 + 9/7*X2.^2 - 2/7*X2 + 1/70);
%       case 25,  ret = (-4*beta+6*beta^2*X1-2*beta^3*X1.^2+1/6*beta^4*X1.^3).*sqrt(44100) ...
%           .*(X2.^4 - 2.*X2.^3 + 9/7*X2.^2 - 2/7*X2 + 1/70);
%     end
%   case 2
%     switch i
%       case 1,  ret = zeros(size(X2));
%         
%       case 2,  ret = zeros(size(X2));
%       case 3,  ret = laguerrePol(1, X1, beta).*sqrt(12);
%       case 4,  ret = laguerrePol(2, X1, beta).*sqrt(12);        
%         
%       case 5,  ret = zeros(size(X2));
%       case 6,  ret = laguerrePol(3, X1, beta).*sqrt(12);
%       case 7,  ret = laguerrePol(1, X1, beta).*sqrt(180).*(2*X2-1);
%       case 8,  ret = laguerrePol(2, X1, beta).*sqrt(180).*(2*X2-1);
%       case 9,  ret = laguerrePol(3, X1, beta).*sqrt(180).*(2*X2-1);
%                
%       case 10,  ret = zeros(size(X2));
%       case 11,  ret = laguerrePol(4, X1, beta).*sqrt(12);  
%       case 12,  ret = laguerrePol(4, X1, beta).*sqrt(180).*(2*X2-1);
%       case 13,  ret = laguerrePol(1, X1, beta).*sqrt(2800) ...
%           .*(3*X2.^2 - 3*X2 + 3/5);
%       case 14,  ret = laguerrePol(2, X1, beta).*sqrt(2800) ...
%           .*(3*X2.^2 - 3*X2 + 3/5);
%       case 15,  ret = laguerrePol(3, X1, beta).*sqrt(2800) ...
%           .*(3*X2.^2 - 3*X2 + 3/5);
%       case 16,  ret = laguerrePol(4, X1, beta).*sqrt(2800) ...
%           .*(3*X2.^2 - 3*X2 + 3/5);
%         
%       case 17,  ret = zeros(size(X2));              
%       case 18,  ret = laguerrePol(5, X1, beta).*sqrt(12);   
%       case 19,  ret = laguerrePol(5, X1, beta).*sqrt(180).*(2*X2-1);
%       case 20,  ret = laguerrePol(5, X1, beta).*sqrt(2800) ...
%           .*(3*X2.^2 - 3*X2 + 3/5);
%       case 21,  ret = laguerrePol(1, X1, beta).*sqrt(44100) ...
%           .*(4*X2.^3 - 6.*X2.^2 + 18/7*X2 - 2/7);
%       case 22,  ret = laguerrePol(2, X1, beta).*sqrt(44100) ...
%           .*(4*X2.^3 - 6.*X2.^2 + 18/7*X2 - 2/7);
%       case 23,  ret = laguerrePol(3, X1, beta).*sqrt(44100) ...
%           .*(4*X2.^3 - 6.*X2.^2 + 18/7*X2 - 2/7);
%       case 24,  ret = laguerrePol(4, X1, beta).*sqrt(44100) ...
%           .*(4*X2.^3 - 6.*X2.^2 + 18/7*X2 - 2/7);
%       case 25,  ret = laguerrePol(5, X1, beta).*sqrt(44100) ...
%           .*(4*X2.^3 - 6.*X2.^2 + 18/7*X2 - 2/7);
%     end % switch
% end % switch
end % function
