% cQ0E0T = { H, U, V }
function lambda = computeLaxFriedrichsCoefficient(cQ0E0T, nuQ0E0T, gConst)
lambda = setNaN2Zero(abs(cQ0E0T{2} .* nuQ0E0T{1} + cQ0E0T{3} .* nuQ0E0T{2})) + sqrt(gConst * cQ0E0T{1});
end % function
