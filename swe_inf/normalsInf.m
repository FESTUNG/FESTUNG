%zeigen immer nach RECHTS von infDir
function nuInf = normalsInf(infDir, nV)
nuInf = zeros(nV,2);
nuInf(:,1) = infDir(:,2);
nuInf(:,2) = -infDir(:,1);
end%function