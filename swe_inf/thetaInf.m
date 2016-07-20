%wenn nn==2 || np==2 => Grenze zum endlichen
function [XP1, XP2] = thetaInf(nn, np, X1, X2)
switch nn
  case 1
    switch np
      case 1,  XP1 = X1; XP2 = X2;%Irrelevant
      case 2,  XP1 = zeros(size(X1));  XP2 = X2;
      case 3,  XP1 = X1;  XP2 = zeros(size(X1));
    end % switch
  case 2
    switch np
      case 1,  XP1 = 1-X2;  XP2 = X2;
      case 2,  XP1 = zeros(size(X1));  XP2 = 1-X2;
      case 3,  XP1 = X2;  XP2 = zeros(size(X1));
    end % switch
  case 3
    switch np
      case 1,  XP1 = X1;  XP2 = ones(size(X1));
      case 2,  XP1 = zeros(size(X1));  XP2 = X1;
      case 3,  XP1 = zeros(size(X1)); XP2 = zeros(size(X1));%Irrelevant
    end % switch
end % switch
end % function
