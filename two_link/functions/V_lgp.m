function [g_est, G_est] = V_lgp(chi_, X, rho_g_sq, P_g_inv2, a, s, pb, D, N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Parametric model estimates
hat_nabla_chi_g_GP = @(q) ...
    [(pb.m1h*pb.l1h/2 + pb.m2h*pb.l1h)*pb.gh*sin(q(1)) + pb.m2h*pb.l2h/2*pb.gh*sin(q(1)+q(2)); ...
    pb.m2h*pb.l2h/2*pb.gh*sin(q(1)+q(2))];
hatV = @(q) ...
    (pb.m1h*pb.l1h/2 + pb.m2h*pb.l1h)*pb.gh*(1-cos(q(1))) + pb.m2h*pb.l2h/2*pb.gh*(1-cos(q(1)+q(2)));

%% Lagrange GP computations

% Initializations
k_yg = zeros(N*D,1);
nabla_chi__k_ygT = zeros(N,N*D);

for d = 1:D
    % Set input/state values
    chi = X(d,1:N)';

    % 
    chi_min_chi_ = chi-chi_;
    chi__min_chi = -chi_min_chi_;

    % g
    chi_min_chi_T_Pginv2 = chi_min_chi_' * P_g_inv2;
    k_g_chichi_times_P_g_inv2 = rho_g_sq * exp(-1/2*chi_min_chi_T_Pginv2*chi_min_chi_) * P_g_inv2;

    % set vecs
    index_vec = 1+(d-1)*N:d*N;
    chi_min_chi_T_Pginv2 = chi_min_chi_' * P_g_inv2;
    nabla_chi__k_ygT(:,index_vec) = rho_g_sq * exp(-1/2*chi_min_chi_T_Pginv2*chi_min_chi_) * P_g_inv2 * (-chi_min_chi_*chi_min_chi_T_Pginv2 + eye(N));
    k_yg(index_vec) = k_g_chichi_times_P_g_inv2 * chi__min_chi;
end

%% Compute estimate

% g
chi_T_Pginv2 = chi_' * P_g_inv2;
nabla_chi__k_0gT = -rho_g_sq * exp( -1/2 * chi_T_Pginv2 * chi_ ) * P_g_inv2 * [chi_, chi_*chi_T_Pginv2 - eye(N)];
g_est = hat_nabla_chi_g_GP(chi_) + nabla_chi__k_ygT * a - nabla_chi__k_0gT * s;

% G
chi_T_Pginv2 = chi_' * P_g_inv2;
k_g_0chi_ = rho_g_sq * exp( -1/2 * chi_T_Pginv2 * chi_ );
k_0gT = k_g_0chi_ * [1, chi_T_Pginv2];
G_est = hatV(chi_) + k_yg' * a - k_0gT * s;

end