function [ ret ] = psi1D( i , X )

switch i
    case 1, ret = ones(size(X));
    case 2, ret = sqrt(3) * (1 - 2 * X);
    case 3, ret = sqrt(5) * ( (6 * X - 6) .* X + 1 );
end  % switch            

end  % function