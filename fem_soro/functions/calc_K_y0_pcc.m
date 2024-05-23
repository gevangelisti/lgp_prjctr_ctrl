function [k_yg_0] = calc_K_y0_pcc(X_tr, rho_g_sq, P_g_inv2, D, N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Initialization
k_yg_0 = zeros(N*D,1);

%%
for d = 1:D
    P_g_inv2_chi = P_g_inv2 * X_tr(d,1:N)';
    k_g_chi0 = rho_g_sq * exp( -1/2 * X_tr(d,1:N) * P_g_inv2_chi );
    k_yg_0(1+(d-1)*N:d*N) = -2 * k_g_chi0 * P_g_inv2_chi;
end

end