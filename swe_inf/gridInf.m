%Kante 2 entspricht Grenze zum
%endlichen, Kante 3 an 1. Knoten der endlichen Kante
%Nummerierung der Elemente nach sortierter Kantennummer
function gInf = gridInf(pd)
g = pd.g;
p = pd.p;
pInf = pd.pInf;
beta = pd.beta;
markE0TOS = pd.g.markE0TOS;

%% numInfElem
gInf.numInfElem = sum(sum(markE0TOS)); %number of elements
if gInf.numInfElem == 0
  return;
end%if
%% E0Inf
gInf.E0Inf = []; %edges neighboring the infinite elements
for i = 1 : 3
  gInf.E0Inf = cat(1,gInf.E0Inf,g.E0T(find(markE0TOS(:,i)),i));
end%for
gInf.E0Inf = sort(gInf.E0Inf);
%% Inf0E
gInf.Inf0E = sparse(zeros(g.numE,0));
gInf.Inf0E(gInf.E0Inf) = 1:gInf.numInfElem;
%% V0Inf 
gInf.V0Inf = zeros(gInf.numInfElem,2); %index of vertex of Inf-element
gInf.numV = 2;
gInf.V0Inf(1,:) = [1 2];
elem = 1;
for i = 1:gInf.numInfElem
  found = 0;
  for k = 1:gInf.numInfElem
    if g.V0E(gInf.E0Inf(elem),2) == g.V0E(gInf.E0Inf(k),1)
      if gInf.V0Inf(k,1) ~= 0 %nur moeglich wenn k=1
        gInf.V0Inf(elem,2) = gInf.V0Inf(k,1);
        gInf.numV = gInf.numV - 1;
        break;
      end%if
      found = 1;
      gInf.V0Inf(k,1) = gInf.numV;
      gInf.numV = gInf.numV + 1;
      gInf.V0Inf(k,2) = gInf.numV;
      elem = k;
      break;
    end%if
  end%for
  %kein Element gefunden mit gleichem Knoten
  if i<gInf.numInfElem && ~found
    elem = find(~gInf.V0Inf(:,1),1);
    gInf.V0Inf(elem,1) = gInf.numV + 1;
    gInf.V0Inf(elem,2) = gInf.numV + 2;
    gInf.numV = gInf.numV + 2;
  end%if
end%for
%% VT0VInf
gInf.VT0VInf = zeros(gInf.numV,1); %Vertex triangle of vertex inf
gInf.VT0VInf(gInf.V0Inf(:,1)) = g.V0E(gInf.E0Inf(1:gInf.numInfElem),1);
gInf.VT0VInf(gInf.V0Inf(:,2)) = g.V0E(gInf.E0Inf(1:gInf.numInfElem),2);
%% markBdr
gInf.markBdr = zeros(gInf.numV,1);
gInf.markBdr(setdiff(gInf.V0Inf(:,1),gInf.V0Inf(:,2))) = 1;
gInf.markBdr(setdiff(gInf.V0Inf(:,2),gInf.V0Inf(:,1))) = 1;
%% markE0IL
%Rand der unendlichen Elemente entsprechen markE0IL
gInf.markE0IL = zeros(gInf.numInfElem, 3);
gInf.markE0IL(:,1) = ismember(gInf.V0Inf(:,2),find(gInf.markBdr));
gInf.markE0IL(:,3) = ismember(gInf.V0Inf(:,1),find(gInf.markBdr));
%% coordV
gInf.coordV = g.coordV(gInf.VT0VInf,:);
%% infDir
gInf.infDir = zeros(gInf.numV,2);
for i = 1 : gInf.numV
  z2 = zeros(g.numE,1);% Kanten von Knoten i
%   g.V2E(:,gInf.VT0VInf(i));
  z2(g.V2E(find(g.V2E(:,gInf.VT0VInf(i))),gInf.VT0VInf(i))) = g.V2E(find(g.V2E(:,gInf.VT0VInf(i))),gInf.VT0VInf(i));
  e1 = 0;
  e2 = 0;
  for j = 1:3
    z = zeros(g.numE,1);% Kanten auf Rand
%     z(g.E0T(find(markE0TL(:,j)+markE0TOS(:,j)),j)) = g.E0T(find(markE0TL(:,j)+markE0TOS(:,j)),j);
%     %eBV: 1 wenn Kante auf Rand und am Knoten i ist sonst 0
%     eBV = z2(g.E0T(find(markE0TL(:,j)+markE0TOS(:,j)),j)) == z(g.E0T(find(markE0TL(:,j)+markE0TOS(:,j)),j));
%     edges = z2(g.E0T(find(markE0TL(:,j)+markE0TOS(:,j)),j));
    z(g.E0T(find(markE0TOS(:,j)),j)) = g.E0T(find(markE0TOS(:,j)),j);
    %eBV: 1 wenn Kante auf OS-Rand und am Knoten i ist sonst 0
    eBV = z2(g.E0T(find(markE0TOS(:,j)),j)) == z(g.E0T(find(markE0TOS(:,j)),j));
    edges = z2(g.E0T(find(markE0TOS(:,j)),j));
    if size(find(eBV),1) == 1
      if e1 == 0
        e1 = edges(find(eBV));
      else
        e2 = edges(find(eBV));
      end%if
    elseif size(find(eBV),1) == 2
      f = find(eBV);
      e1 = edges(f(1));
      e2 = edges(f(2));
    end%if
  end%for
  if e2 == 0
    gInf.infDir(i,:) = g.nuE(e1,:);
  else
    gInf.infDir(i,:) = 0.5 * (g.nuE(e1,:) + g.nuE(e2,:));
  end % if
  gInf.infDir(i,:) = gInf.infDir(i,:)/norm(gInf.infDir(i,:));
end%for
%% nuE0Inf
gInf.nuE0Inf = zeros(gInf.numInfElem, 3, 2);
gInf.nuE0Inf(:,1,1) = -gInf.infDir(gInf.V0Inf(:,2),2);
gInf.nuE0Inf(:,1,2) = gInf.infDir(gInf.V0Inf(:,2),1);
gInf.nuE0Inf(:,2,1) = -g.nuE(gInf.E0Inf, 1);
gInf.nuE0Inf(:,2,2) = -g.nuE(gInf.E0Inf, 2);
gInf.nuE0Inf(:,3,1) = gInf.infDir(gInf.V0Inf(:,1),2);
gInf.nuE0Inf(:,3,2) = -gInf.infDir(gInf.V0Inf(:,1),1);
%% trafoParam
%for transformation from reference to physical element
[gInf.M, gInf.H, gInf.coordAxInf, gInf.Alpha, gInf.Beta, gInf.xHat] = compTrafoParamInf(gInf);
%% detQuadPoints
qOrd = max(2*p,1);
qOrdInf = max(2*pInf,1);
[Q1, Q2, ~] = quadRuleInf2D(qOrd, qOrdInf, beta);
gInf.detQuadPoints = zeros(gInf.numInfElem,size(Q1,2));
for nq = 1 : size(Q1,2)
  gInf.detQuadPoints(:,nq) = (cos(gInf.Alpha)-sin(gInf.Beta).*2.*gInf.M.*(Q2(nq)-0.5)).*cos(gInf.Beta).*(gInf.H + 2*Q1(nq)*gInf.M) + ...
    sin(gInf.Beta).*(gInf.H+2*Q1(nq)*gInf.M).*(sin(gInf.Alpha) + 2*cos(gInf.Beta).*gInf.M*(Q2(nq)-0.5));
end%for
%gInf.detQuadPoints = abs(gInf.detQuadPoints);
%% detE0I
gInf.detE0I = zeros(gInf.numInfElem, 3);
gInf.detE0I(:,2) = g.areaE(gInf.E0Inf);
X0 = zeros(gInf.numInfElem,2);
X1 = zeros(gInf.numInfElem,2);
X0(:,1) = gInf.coordV(gInf.V0Inf(:,2),1);
X0(:,2) = gInf.coordV(gInf.V0Inf(:,2),2);
X1(:,1) = gInf.coordV(gInf.V0Inf(:,1),1);
X1(:,2) = gInf.coordV(gInf.V0Inf(:,1),2);
Y = X0-X1;
normY = sqrt(Y(:,1).^2 + Y(:,2).^2);
Y(:,1) = Y(:,1)./normY;
Y(:,2) = Y(:,2)./normY;
K = gInf.xHat + gInf.coordAxInf + Y.*(0.5*cat(2,gInf.H,gInf.H) + cat(2,gInf.M,gInf.M)) - X0;
gInf.detE0I(:,1) = sqrt(K(:,1).^2 + K(:,2).^2);
K = gInf.xHat + gInf.coordAxInf - Y.*(0.5*cat(2,gInf.H,gInf.H) + cat(2,gInf.M,gInf.M)) - X1;
gInf.detE0I(:,3) = sqrt(K(:,1).^2 + K(:,2).^2);
%% markE0IE0I
gInf.markE0IE0I = cell(3,3);
for nn = 1 : 3
  for np = 1 : 3
    gInf.markE0IE0I{nn,np} = sparse(gInf.numInfElem, gInf.numInfElem);
    if nn == 2 || np == 2 || nn == np%es kann nur 1-3 eine gemeinsame Kante sein
      continue;
    end%if
    if nn == 1 %np == 3
      [~,kH] = ismember(gInf.V0Inf(:,2),gInf.V0Inf(:,1));
      [kn,~,kp] = find(kH);
      for i = 1:length(kn)
        gInf.markE0IE0I{nn,np}(kn(i),kp(i)) = 1;
      end % for
    end%if
    if nn == 3 %np == 1
      gInf.markE0IE0I{nn,np}=gInf.markE0IE0I{1,3}';
    end%if  
  end % for
end % for
%% markE0Iint
gInf.markE0Iint = zeros(gInf.numInfElem,3);
gInf.markE0Iint(:,1) = ~gInf.markE0IL(:,1);
gInf.markE0Iint(:,3) = ~gInf.markE0IL(:,3);
%% markE0Iaux
gInf.markE0Iaux = cell(3,3); % auxiliary cell of vectors of length K, needed in some routines
for nn = 1:3
  for np = 1:3
    gInf.markE0Iaux{nn,np} = gInf.markE0IE0I{nn,np} * ones(gInf.numInfElem,1);
  end
end

gInf.T0I = g.T0E(gInf.E0Inf,1);
%% markE0TE0I
gInf.markE0TE0I = cell(3,3);
for nn = 1:3
  for np = 1:3
    if nn==2
      gInf.markE0TE0I{nn,np} = sparse(gInf.numInfElem, g.numT);
      if sum(markE0TOS(gInf.T0I,np))>0
        gInf.markE0TE0I{nn,np}(sub2ind(size(gInf.markE0TE0I{nn,np}),find(markE0TOS(gInf.T0I,np)),gInf.T0I(find(markE0TOS(gInf.T0I,np))))) = 1;
      end%if
    elseif np == 2
      gInf.markE0TE0I{nn,np} = sparse(g.numT, gInf.numInfElem);
      if sum(markE0TOS(gInf.T0I,nn))>0
        gInf.markE0TE0I{nn,np}(sub2ind(size(gInf.markE0TE0I{nn,np}),gInf.T0I(find(markE0TOS(gInf.T0I,nn))),find(markE0TOS(gInf.T0I,nn)))) = 1;
      end%if
    else
      gInf.markE0TE0I{nn,np} = sparse(1,1);
    end%if
  end%for
end%for
if sum(sum(markE0TOS(gInf.T0I,:),2)>1)>0   %  at least one Triangle has more than 1 OS-Edge
  infElems = find(sum(markE0TOS(gInf.T0I,:),2)>1);
  OSEdges = g.E0T(gInf.T0I(infElems),:).*markE0TOS(gInf.T0I(infElems),:);
  edges = zeros(g.numE,1);
  for i = 1:length(infElems)
    s = sort(OSEdges(i,:));
    for j = 1:3
      if s(j)>0 && edges(s(j))==0
        edges(s(j))=1;
        markLocalEdge = (OSEdges(i,:)==[s(j),s(j),s(j)]);
        gInf.markE0TE0I{2,1}(sub2ind(size(gInf.markE0TE0I{2,1}),infElems(i),gInf.T0I(infElems(i)))) = markLocalEdge(1);
        gInf.markE0TE0I{2,2}(sub2ind(size(gInf.markE0TE0I{2,2}),infElems(i),gInf.T0I(infElems(i)))) = markLocalEdge(2);
        gInf.markE0TE0I{2,3}(sub2ind(size(gInf.markE0TE0I{2,3}),infElems(i),gInf.T0I(infElems(i)))) = markLocalEdge(3);
        gInf.markE0TE0I{1,2}(sub2ind(size(gInf.markE0TE0I{1,2}),gInf.T0I(infElems(i)),infElems(i))) = markLocalEdge(1);
        gInf.markE0TE0I{3,1}(sub2ind(size(gInf.markE0TE0I{3,1}),gInf.T0I(infElems(i)),infElems(i))) = markLocalEdge(3);
        break;
      end%if
    end%for
  end%for
end%if

%% markE0ITaux
gInf.markE0ITaux= cell(3,3); % entry {3,3} corresponds to: nn=np=2, fin = negative
for nn = 1:3
  for np = 1:3
    if nn==2
      gInf.markE0ITaux{nn,np} = gInf.markE0TE0I{nn,np} * ones(g.numT,1);
    elseif np == 2
      gInf.markE0ITaux{nn,np} = gInf.markE0TE0I{nn,np} * ones(gInf.numInfElem,1);
    end%if
  end%for
end%for
gInf.markE0ITaux{3,3} = gInf.markE0TE0I{2,2}' * ones(gInf.numInfElem,1);

end%function