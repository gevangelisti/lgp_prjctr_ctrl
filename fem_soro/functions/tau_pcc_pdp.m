function [tau] = tau_pcc_pdp(t,q,dq,q_des,dq_des,ddq_des,Kp,Kd,mu,L,Iz,gp,k_pcc,d_pcc)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[M,C,g] = MCg_pcc_n4(mu,L,Iz,gp,q,dq);

tau = M*ddq_des(t) + C*dq_des(t) + g + diag(k_pcc)*q + diag(d_pcc)*dq ...
    - Kp * (q - q_des(t)) - Kd * (dq - dq_des(t));

end