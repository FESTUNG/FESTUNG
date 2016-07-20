function [Q1, Q2, W] = quadRuleInf2D(qOrd, qOrdInf, beta)

[Q1_1D,W1_1D]=quadRuleInf1D(qOrdInf, beta);
[Q2_1D,W2_1D]=quadRule1D(qOrd);

W = kron(W1_1D,W2_1D);
E1 = ones(1,length(Q2_1D'));
E2 = ones(1,length(Q1_1D'));

Q1=kron(Q1_1D,E1);
Q2=kron(E2,Q2_1D);
end % function
