function [ddq] = ivd_MCg_pcc_n4(mu,L,Iz,gp,q,dq,tau)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
[M,C,g] = MCg_pcc_n4(mu,L,Iz,gp,q,dq);

ddq = M \ (tau - C*dq - g);

end