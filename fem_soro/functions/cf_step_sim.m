function [cost] = cf_step_sim(X, z, y_D, Z, Upsilon, pb, t_fem, q_fem, a_pcc, t_fem_2, q_fem_2, mu_D, D, N, x)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% D = size(X,1); % # of training pairs
% N = size(X,2)/2; % # of DOFs
% x: optimization variables (= GP hyperparameters)
% n = 2
% 1 to n(n+1)/2=3: sigma_ij (signal variances of mass matrix components)
% n(n+1)/2+1=4 to 3*n(n+1)/2=9: Lambda_ij (length scales)
% 3*n(n+1)/2+1=10: rho_g
% 3*n(n+1)/2+2=11 to 3*n(n+1)/2+1+n=12: P_g

% pb.d = pb.d + x(end-2*N+1:end-N)';
% pb.k = pb.k + x(end-N+1:end)';
% 
% mu_D = mu_D + repmat(x(end-2*N+1:end-N)', [D,1]) .* reshape(X(:,N+1:end),[N*D,1]) + repmat(x(end-N+1:end)', [D,1]) .* reshape(X(:,1:N),[N*D,1]);
% x = x(1:end-2*N);

%% psi
sigma_n = sqrt(Z(1,1));
sigma_tau = sqrt(Upsilon(1,1));
psi = [sigma_n, sigma_tau, x];
[LLT_Sigma_sq,lambda_isq,rho_g_sq,P_g_inv2,~] = get_hyp_params(N,psi);
eps_sq = 0;

%% Computations

% compute covariance matrices
K_y_plus_Pi = calc_Kyy_pcc_mex(X, z, LLT_Sigma_sq, lambda_isq, rho_g_sq, P_g_inv2, Z, Upsilon, pb, D, N);
K_y0 = calc_K_y0_pcc_mex(X, rho_g_sq, P_g_inv2, D, N);
k_0_inv = 1/(2*rho_g_sq);
A = K_y_plus_Pi + eps_sq*eye(length(y_D)) - k_0_inv*(K_y0*K_y0');

% cholesky factorization
try
    L_A = chol(A,'lower');
catch
    cost = Inf;
    return
end

% compute inversion using factorization
y_min_mu = y_D-mu_D;
A_inv_times_y_min_mu = L_A'\(L_A\y_min_mu); % the backslash operator recognizes triangular systems!
K_times_a = k_0_inv * (K_y0' * A_inv_times_y_min_mu);

% compute estimates
[q_pcc_lgp_sim,~] = sim_lgp_step_syf_timed(t_fem, a_pcc, X, z, LLT_Sigma_sq, lambda_isq, rho_g_sq, P_g_inv2, A_inv_times_y_min_mu, K_times_a, pb, D, N, 20);
[q_pcc_lgp_sim_2,~] = sim_lgp_step_syf_timed(t_fem_2, a_pcc/2, X, z, LLT_Sigma_sq, lambda_isq, rho_g_sq, P_g_inv2, A_inv_times_y_min_mu, K_times_a, pb, D, N, 20);

% cost
if sum(sum(isnan(q_pcc_lgp_sim))) == 0 && sum(sum(diff(q_pcc_lgp_sim') == 0)) == 0 && sum(sum(isnan(q_pcc_lgp_sim_2))) == 0 && sum(sum(diff(q_pcc_lgp_sim_2') == 0)) == 0
    if (size(q_pcc_lgp_sim) == size(q_fem)) .* (size(q_pcc_lgp_sim_2) == size(q_fem_2))
        cost = sum(sum((q_fem-q_pcc_lgp_sim).^2)) + sum(sum((q_fem_2-q_pcc_lgp_sim_2).^2));
    else
        cost = Inf;
        w = warning('query','last');
        id = w.identifier;
        warning('off',id)
    end
else
    cost = Inf;
end

end 