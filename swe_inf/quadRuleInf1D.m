
function [Q, W] = quadRuleInf1D(qOrd, beta)

numExp = 2;

betaInt = beta/2*numExp; %Quadraturregel fuer p(x)*exp(-beta*x)

NQ = ceil((qOrd-1)/2);%0-indiziert => bei NQ=0: 1 Quadraturpunkt
% a=(2*(0:1:NQ)+1)/(betaInt);
% b=(1:1:NQ).^2/(betaInt);
% if length(a)== 1
%   Q = a;
% else
%   A = diag(a, 0) + diag(sqrt(b), -1) + diag(sqrt(b), 1);
%   A
%   Q = eig(A);
% end
Q = laguerrePolZeros(NQ+2,betaInt);
W = (Q./((NQ+1)^2*laguerrePol(NQ+1, Q, betaInt).^2).*exp(betaInt*Q))';
%exp kann man weglassen wenn man bei den Integralen exp(-...) weglaesst
%dann ist sie aber nicht mehr allgemein mit "betaInt" nutzbar!!!
Q=Q';

end % function
