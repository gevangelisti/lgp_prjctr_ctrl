function [K_y0] = calc_K_y0(X_tr, rho_g_sq, P_g_inv2, D, N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Initializations
K_y0 = zeros(N*D,N+1);
k_yg_0 = zeros(N*D,1);
k_yg_nabla_chi__0_mat = zeros(N*D,N);

%%

for d = 1:D
    P_g_inv2_chi = P_g_inv2 * X_tr(d,1:N)';
    k_g_chi0 = rho_g_sq * exp( -1/2 * X_tr(d,1:N) * P_g_inv2_chi );
    k_yg_0(1+(d-1)*N:d*N) = -k_g_chi0 * P_g_inv2_chi;
    k_yg_nabla_chi__0_mat(1+(d-1)*N:d*N,:) = k_g_chi0 * (-P_g_inv2_chi*P_g_inv2_chi' + P_g_inv2);
end

K_y0(:,1) = k_yg_0;
K_y0(:,2:end) = k_yg_nabla_chi__0_mat;

end