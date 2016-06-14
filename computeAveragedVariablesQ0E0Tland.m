function cRoePikeQ0E0T = computeAveragedVariablesQ0E0Tland(cQ0E0Tint, cQ0E0Text, hQ0E0Tint, hQ0E0Text, markQ0E0Tbdr, averagingType)
hIntExtInv = 0.5 ./ hQ0E0Tint(markQ0E0Tbdr);
cRoePikeQ0E0T = { zeros(size(hQ0E0Tint)); zeros(size(hQ0E0Tint)); zeros(size(hQ0E0Tint)) };
cRoePikeQ0E0T{1}(markQ0E0Tbdr) = hQ0E0Tint(markQ0E0Tbdr);
cRoePikeQ0E0T{2}(markQ0E0Tbdr) = (cQ0E0Tint{2}(markQ0E0Tbdr) + cQ0E0Text{2}(markQ0E0Tbdr)) .* hIntExtInv;
cRoePikeQ0E0T{3}(markQ0E0Tbdr) = (cQ0E0Tint{3}(markQ0E0Tbdr) + cQ0E0Text{3}(markQ0E0Tbdr)) .* hIntExtInv;
end

