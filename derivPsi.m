function [ ret ] = derivPsi( i , X )

switch i
    case 1, ret = zeros(size(X));
    case 2, ret = -2 * sqrt(3) * ones(size(X));
    case 3, ret = 6 * sqrt(5) * (2 * X - 1);
end  % switch

end  % function