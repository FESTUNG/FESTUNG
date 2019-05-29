function basesOnQuadInf = computeBasesOnQuadInf(N, NInf, pd)
beta = pd.beta;
%global gPhiInf2D gGradPhiInf2D gPhiInf1D gThetaPhiInf1D gThetaPhiLinkT2I21D%PhiLink for the linking edge
p    = (sqrt(8*N+1)-3)/2;
pInf = NInf/(p+1)-1;
if pInf > 0, requiredOrdersInf = [2*pInf, 2*pInf+1]; else requiredOrdersInf = 1; end
if p > 0,    requiredOrders    = [2*p   , 2*p+1   ]; else requiredOrders    = [1 1]; end
basesOnQuadInf.gPhiInf2D = cell(max(requiredOrdersInf),1);  basesOnQuadInf.gGradPhiInf2D  = cell(max(requiredOrdersInf),1);
basesOnQuadInf.gPhiInf1D = cell(max(requiredOrdersInf),1);  basesOnQuadInf.gThetaPhiInf1D = cell(max(requiredOrdersInf),1);
basesOnQuadInf.gThetaPhiLinkT2I21D = cell(max(requiredOrdersInf),1);
for it = 1 : length(requiredOrdersInf)
  ordInf = requiredOrdersInf(it);
  ord    = requiredOrders(it);
  [Q1, ~, ~] = quadRuleInf2D(ord, ordInf, beta);
  basesOnQuadInf.gPhiInf2D{ordInf}      = zeros(length(Q1), NInf);
  basesOnQuadInf.gGradPhiInf2D{ordInf}  = zeros(length(Q1), NInf, 2);
  [Q1, Q2, ~] = quadRuleInf2D(ord, ordInf, beta);
  for i = 1 : NInf
    basesOnQuadInf.gPhiInf2D{ordInf}(:, i) = phiInf(i, p, Q1, Q2, beta);
  end % for  
  for m = 1 : 2
    for i = 1 : NInf
      basesOnQuadInf.gGradPhiInf2D{ordInf}(:, i, m) = gradPhiInf(i, p, m, Q1, Q2, beta);
    end % for
  end % for
  [Q, ~] = quadRuleInf1D(ordInf, beta);
  basesOnQuadInf.gThetaPhiInf1D{ordInf} = zeros(length(Q), NInf, 3, 3);
  basesOnQuadInf.gThetaPhiLinkT2I21D{ordInf} = zeros(length(Q),NInf);
  basesOnQuadInf.gPhiInf1D{ordInf} = zeros(length(Q), NInf, 3);
  for nn = 1 : 3 %nn=2 entspricht Kante zum Dreieck
    if nn == 2
      [Q, ~] = quadRule1D(ord);
    else
      [Q, ~] = quadRuleInf1D(ordInf, beta);
    end % if 
    [Q1, Q2] = gammaMapInf(nn, Q);
    for i = 1 : NInf
      for r = 1 : length(Q1)
        basesOnQuadInf.gPhiInf1D{ordInf}(r, i, nn) = phiInf(i, p, Q1(r), Q2(r), beta);
      end % for
    end
    for np = 1 : 3
      if np == 2 || nn == 2
        [Q, ~] = quadRule1D(ord);%Inf);
        [Q1, Q2] = gammaMap(nn, Q);  
        [QP1,QP2] = theta(nn, np, Q1, Q2);
      else
        [Q, ~] = quadRuleInf1D(ordInf, beta);
        [Q1, Q2] = gammaMapInf(nn, Q);
        [QP1,QP2] = thetaInf(nn, np, Q1, Q2);
      end % if
      if nn == 2
        for i = 1 : N
          for r = 1 : length(QP1)
            basesOnQuadInf.gThetaPhiInf1D{ordInf}(r, i, nn, np) = phi(i, QP1(r), QP2(r));
          end % for
        end % for
        if np == 2
          for i = 1:NInf
            for r = 1 : length(QP1)
              basesOnQuadInf.gThetaPhiLinkT2I21D{ordInf}(r,i) = phiInf(i, p, QP1(r), QP2(r), beta);
            end % for
          end % for
        end % if
      else
        if np == 2
          for i = 1 : NInf
            for r = 1 : length(QP1)
              basesOnQuadInf.gThetaPhiInf1D{ordInf}(r, i, nn, np) = phiInf(i, p, QP1(r), QP2(r), beta);
            end % for
          end % for
        else 
          for i = 1 : NInf
            basesOnQuadInf.gThetaPhiInf1D{ordInf}(:, i, nn, np) = phiInf(i, p, QP1, QP2, beta);
          end % for
        end % if
      end % if
    end % for
  end % for
end % for
end % function
