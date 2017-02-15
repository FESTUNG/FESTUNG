function rhs = evalExactUHNu(g, t, HeightTimesVelocityIndex, UHexact, basesOnQuad)

% AR: Achtung! Hier wird davon ausgegangen, dass das 1D Gitter wie folgt
% aufgebaut ist (lexikographische Nummerierung):
% | --1-- | --2-- | --3-- | ...
% Diese Annahme an die Nummerierung der Elemente wird in den Tests auch
% nicht verletzt. Ich weiß aber nicht, ob das immer so ist :(

numBases = size(basesOnQuad.phi1D, 2);
rhs = zeros(g.numT * numBases, 1);

for i = 1 : 2
    switch HeightTimesVelocityIndex(i)
        case 2 % right boundaries
            rhs(end - numBases + 1: end) = +1 * UHexact(t,g.x_max) * basesOnQuad.phi0D(:, 2);
        case 4 % left boundaries
            rhs(1: numBases) = -1 * UHexact(t,g.x_min) * basesOnQuad.phi0D(:, 1);
    end
end