function [LLT_Sigma_sq,lambda_isq,...
    rho_g_sq,P_g_inv2,P_g_sq,...
    L_Sigma_sq_dmp,lambda_isq_dmp] = get_hyp_params(N,psi)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% kinetic energy length scale parameters
% Lambda parameter ordering = Lower Triangle:
% lambda_11(1), ... lambda_11(N), lambda_21(1), ... lambda_21(N), lambda_22(1) ... lambda_22(N), ...
L_Sigma_sq = get_cholesky_hyp_matrix(N,2,psi);
LLT_Sigma_sq = L_Sigma_sq*L_Sigma_sq';
offset = 3 + N*(N+1)/2;
lambda_isq = 1/psi(offset)^2;

% potential energy hyperparameters
offset = offset + 1;
rho_g_sq = psi(offset)^2;
P_g_inv2 = diag(1./(psi(offset+1:offset+N).^2));
P_g_sq = diag(psi(offset+1:offset+N).^2);

% dissipation hyperparameters
offset = offset + N;
L_Sigma_sq_dmp = zeros(N);
lambda_isq_dmp = zeros(N);
if length(psi) > offset
    L_Sigma_sq_dmp = get_cholesky_hyp_matrix(N,offset,psi);
    offset = offset + 1 + N*(N+1)/2;
    lambda_isq_dmp = diag(1./psi(offset:end)).^2;
end

end

function L_Sigma_sq = get_cholesky_hyp_matrix(N,offset_prev,psi)

Sigma = zeros(N);
for i = 1:N
    offset = offset_prev + i*(i-1)/2;
    for j = 1:i-1
        Sigma(i,j) = psi(offset+j);
    end
    Sigma(i,i) = psi(offset+i);
end
L_Sigma_sq = Sigma.^2;

end