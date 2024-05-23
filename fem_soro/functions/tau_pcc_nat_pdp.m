function [tau] = tau_pcc_nat_pdp(t,q,dq,q_des,dq_des,ddq_des,Kd,Ncc,mu,L,Iz,gp,k_pcc,d_pcc)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[M,C,g] = MCg_pcc_n4(mu,L,Iz,gp,q,dq);
[~,~,ge] = MCg_pcc_n4(mu,L,Iz,gp,q-q_des(t),dq);

g = g + diag(k_pcc)*q;
ge = ge + diag(k_pcc)*(q-q_des(t));
d = diag(d_pcc)*dq;
de = diag(d_pcc)*(dq-dq_des(t));
tau = M*ddq_des(t) + C*dq_des(t) ...
    + (eye(Ncc)-heaviside((q-q_des(t))'*g)*(q-q_des(t))*(q-q_des(t))'./(1e-3 + norm((q-q_des(t)))^2) )*g - ge ...
    + (eye(Ncc)-heaviside((dq-dq_des(t))'*d)*(dq-dq_des(t))*(dq-dq_des(t))'./(1e-3 + norm((dq-dq_des(t)))^2) )*d - de - Kd*(dq-dq_des(t));

end