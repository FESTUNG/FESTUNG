function cRoePikeQ0E0T = computeAveragedVariablesQ0E0Tint(cQ0E0Tint, cQ0E0Text, hQ0E0Tint, hQ0E0Text, averagingType)
cRoePikeQ0E0T = cell(3,1);

switch averagingType
  case 'full-harmonic'
    hL = sqrt(hQ0E0Tint);
    hR = sqrt(hQ0E0Text);
    hIntExtInv = 1 ./ (hR .* hQ0E0Tint + hL .* hQ0E0Text);
    
    cRoePikeQ0E0T{1} = (hL .* hQ0E0Tint + hR .* hQ0E0Text) ./ (hL + hR);
    cRoePikeQ0E0T{2} = (hL .* cQ0E0Tint{2} + hR .* cQ0E0Text{2}) .* hIntExtInv;
    cRoePikeQ0E0T{3} = (hL .* cQ0E0Tint{3} + hR .* cQ0E0Text{3}) .* hIntExtInv;

  case 'semi-harmonic'
    hL = sqrt(hQ0E0Tint);
    hR = sqrt(hQ0E0Text);
    hIntExtInv = 1 ./ (hR .* hQ0E0Tint + hL .* hQ0E0Text);
    
    cRoePikeQ0E0T{1} = 0.5 * (hQ0E0Tint + hQ0E0Text);
    cRoePikeQ0E0T{2} = (hR .* cQ0E0Tint{2} + hL .* cQ0E0Text{2}) .* hIntExtInv;
    cRoePikeQ0E0T{3} = (hL .* cQ0E0Tint{3} + hR .* cQ0E0Text{3}) .* hIntExtInv;

  case 'mean'
    hIntExtInv = 0.5 ./ (hQ0E0Tint + hQ0E0Text);
    
    cRoePikeQ0E0T{1} = 0.5 * (hQ0E0Tint + hQ0E0Text);
    cRoePikeQ0E0T{2} = (hQ0E0Text .* cQ0E0Tint{2} + hQ0E0Tint .* cQ0E0Text{2}) .* hIntExtInv;
    cRoePikeQ0E0T{3} = (hQ0E0Text .* cQ0E0Tint{3} + hQ0E0Tint .* cQ0E0Text{3}) .* hIntExtInv;

  otherwise
    error('Unknown averaging type for interior edges')
end % switch
end

