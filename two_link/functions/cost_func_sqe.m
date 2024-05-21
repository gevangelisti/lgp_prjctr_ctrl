function [cost] = cost_func_sqe(...
    X, z, y_D, X_, z_, tau_, Z, Upsilon, pb, mu_D, D, N, eps_sq, x)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

sigma_n = sqrt(Z(1,1));
sigma_tau = sqrt(Upsilon(1,1));
psi = [sigma_n, sigma_tau, [x(1:8), 0, x(9:end)]];
[LLT_Sigma_sq,lambda_isq,rho_g_sq,P_g_inv2,P_g_sq,L_Sigma_sq_dmp,lambda_isq_dmp] = get_hyp_params(N,psi);

%%

% compute covariance matrices
K_y_plus_Pi = calc_K_yy_mex(X, z, LLT_Sigma_sq, lambda_isq, rho_g_sq, P_g_inv2, L_Sigma_sq_dmp, lambda_isq_dmp, Z, Upsilon, pb, D, N);
K_y0 = calc_K_y0_mex(X, rho_g_sq, P_g_inv2, D, N);
K_0_inv = (1/rho_g_sq) * [1, zeros(1,N); zeros(N,1), P_g_sq];
A = K_y_plus_Pi + eps_sq*eye(length(y_D)) - K_y0*K_0_inv*K_y0';

% cholesky factorization
try
    L_A = chol(A,'lower');
catch
    cost = Inf;
    return
end

% compute inversion using factorization
y_min_mu = y_D - mu_D;
A_inv_times_y_min_mu = L_A'\(L_A\y_min_mu);
K_times_a = K_0_inv * (K_y0' * A_inv_times_y_min_mu);

tau_tr = [y_D(1:2:end),y_D(2:2:end)];
[tau_est_tr,~,~,~,~,~,~,~] = ...
    calc_lgp_ests_mex(X, z, LLT_Sigma_sq, lambda_isq, rho_g_sq, P_g_inv2, L_Sigma_sq_dmp, lambda_isq_dmp, [X; X_], [z; z_], [tau_tr; tau_], pb, A_inv_times_y_min_mu, K_times_a, D, N);
cost = sum(sum((tau_est_tr-[tau_tr; tau_]).^2));

end