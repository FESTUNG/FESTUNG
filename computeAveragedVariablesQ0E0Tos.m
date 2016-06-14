function cRoePikeQ0E0T = computeAveragedVariablesQ0E0Tos(cQ0E0Tint, cQ0E0Text, hQ0E0Tint, hQ0E0Text, markQ0E0Tbdr, averagingType)
cRoePikeQ0E0T = { zeros(size(hQ0E0Tint)); zeros(size(hQ0E0Tint)); zeros(size(hQ0E0Tint)) };

switch averagingType
  case 'full-harmonic'
    hL = sqrt(hQ0E0Tint(markQ0E0Tbdr));
    hR = sqrt(hQ0E0Text(markQ0E0Tbdr));
    hIntExtInv = 1 ./ (hL .* hR);
    
    cRoePikeQ0E0T{1}(markQ0E0Tbdr) = (hL .* hQ0E0Tint(markQ0E0Tbdr) + hR .* hQ0E0Text(markQ0E0Tbdr)) ./ (hL + hR);
    cRoePikeQ0E0T{2}(markQ0E0Tbdr) = cQ0E0Tint{2}(markQ0E0Tbdr) .* hIntExtInv;
    cRoePikeQ0E0T{3}(markQ0E0Tbdr) = cQ0E0Tint{3}(markQ0E0Tbdr) .* hIntExtInv;

  case 'semi-harmonic'
    hL = sqrt(hQ0E0Tint(markQ0E0Tbdr));
    hR = sqrt(hQ0E0Text(markQ0E0Tbdr));
    hIntExtInv = 1 ./ (hL .* hR);
    
    cRoePikeQ0E0T{1}(markQ0E0Tbdr) = 0.5 * (hQ0E0Tint(markQ0E0Tbdr) + hQ0E0Text(markQ0E0Tbdr));
    cRoePikeQ0E0T{2}(markQ0E0Tbdr) = cQ0E0Tint{2}(markQ0E0Tbdr) .* hIntExtInv;
    cRoePikeQ0E0T{3}(markQ0E0Tbdr) = cQ0E0Tint{3}(markQ0E0Tbdr) .* hIntExtInv;

  case 'mean'
    cRoePikeQ0E0T{1}(markQ0E0Tbdr) = 0.5 * (hQ0E0Tint(markQ0E0Tbdr) + hQ0E0Text(markQ0E0Tbdr));
    
    hIntExtInv = cRoePikeQ0E0T{1}(markQ0E0Tbdr) ./ (hQ0E0Tint(markQ0E0Tbdr) .* hQ0E0Text(markQ0E0Tbdr));
    cRoePikeQ0E0T{2}(markQ0E0Tbdr) = cQ0E0Tint{2}(markQ0E0Tbdr) .* hIntExtInv;
    cRoePikeQ0E0T{3}(markQ0E0Tbdr) = cQ0E0Tint{3}(markQ0E0Tbdr) .* hIntExtInv;

  otherwise
    error('Unknown averaging type for interior edges')
end % switch
end

