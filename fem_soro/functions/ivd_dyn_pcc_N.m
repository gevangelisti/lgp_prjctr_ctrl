function [tau,g] = ivd_dyn_pcc_N(pb,q,dq,ddq)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

MCg_func = @(q,dq) MCg_pcc_n4(pb.mu,pb.L,pb.Iz,pb.g,q,dq);
[M,C,g] = MCg_func(q,dq);

tau = M*ddq + C*dq + g + diag(pb.d)*dq + diag(pb.k)*q;

end