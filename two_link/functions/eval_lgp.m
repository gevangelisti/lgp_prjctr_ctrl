function [tau_est,ddq_est,M_est,d2_f_est,nabla_chi_f_est,nabla_chi_g_est,f_est,g_est,D_est,P_dsp_est] = ...
    eval_lgp(...
    chi_, xi_, tau_, ddq_, X, ddq_meas, ...
    LLT_Sigma_sq, lambda_isq, rho_g_sq, P_g_inv2, L_Sigma_sq_dmp, lambda_isq_dmp, ...
    a, s, pb, D, N ...
    )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Parametric model estimates
hatM = @(q) ...
    [pb.alphah + 2*pb.betah*cos(q(2)), ...
    pb.deltah + pb.betah*cos(q(2)); ...
    pb.deltah + pb.betah*cos(q(2)), ...
    pb.deltah];
hat_M_dot = @(q,dq) ...
    [0, ...
    -pb.betah*sin(q(2)).*(2*dq(1)+dq(2)); ...
    0, ...
    -pb.betah*sin(q(2)).*dq(1)];
hat_nabla_chi_f_GP = @(q,dq) ...
    [0; ...
    -pb.betah*sin(q(2)).*dq(1).*(dq(1)+dq(2))];
hat_nabla_chi_g_GP = @(q) ...
    [(pb.m1h*pb.l1h/2 + pb.m2h*pb.l1h)*pb.gh*sin(q(1)) + pb.m2h*pb.l2h/2*pb.gh*sin(q(1)+q(2)); ...
    pb.m2h*pb.l2h/2*pb.gh*sin(q(1)+q(2))];
hatV = @(q) ...
    (pb.m1h*pb.l1h/2 + pb.m2h*pb.l1h)*pb.gh*(1-cos(q(1))) + pb.m2h*pb.l2h/2*pb.gh*(1-cos(q(1)+q(2)));
hatD = @(dq) [pb.dh(1)+pb.dh(2)*abs(dq(1)), 0; 0, pb.dh(1)+pb.dh(2)*abs(dq(2))];

%% Lagrange GP computations

% Initializations
M_GP = zeros(N);
d2_f_GP_dxi_dchi_ = zeros(N);
nabla_chi__k_yfT_pos_part = zeros(N,N*D);
nabla_chi__k_yfT_neg_part = zeros(N,N*D);
nabla_chi__k_ygT = zeros(N,N*D);
k_yf_pos_part_vec = zeros(N*D,1);
k_yf_neg_part_vec = zeros(N*D,1);
k_yg = zeros(N*D,1);
D_GP = zeros(N);

% Initializations
mu = zeros(N,1,N,N);
Xi_Mu = zeros(N);
Alpha_Mu = zeros(N);
mat_xi__alpha_Mu_ = zeros(N);
mat_xi__z_Mu_ = zeros(N);
Omega = zeros(N,N,N,N);
xiT_Omega = zeros(1,N,N,N);
Omega_xi = zeros(N,1,N,N);
mat_xi_xi__xiTOmega = zeros(N);
mat_xi__alpha_xiTOmega = zeros(N);
alphaT_Omega = zeros(1,N,N,N);
mat_xi_xi__alphaTOmega = zeros(N);
mat_xi_xi__Mu = zeros(N);
mat_Mu_xi__z = zeros(N);
mat_Omegaxi_xi__xi = zeros(N);

%
xi__xi_T = xi_ * xi_';
diag_xi_ = diag(xi_);

for d = 1:D
    % Set input/state values
    chi = X(d,1:N)';
    xi = X(d,N+1:2*N)';
    alpha = a(1+(d-1)*N:d*N);
    z = ddq_meas(d,:)';

    % 
    chi_min_chi_ = chi-chi_;
    LLT = exp(-lambda_isq * (chi_min_chi_' * chi_min_chi_)) * LLT_Sigma_sq;
    chi__min_chi = -chi_min_chi_;
    diag_chi__min_chi = diag(chi__min_chi);
    dchi__dchi_kf = zeros(N);
    dkf_dchi = zeros(N,1);
    %
    xi_min_xi_ = xi-xi_;
    d_mat_dmp = -1/2*lambda_isq_dmp .* [xi_min_xi_(1)^2, 0; xi_min_xi_'*xi_min_xi_, xi_min_xi_(2)^2];
    L_dmp = L_Sigma_sq_dmp .* [exp(d_mat_dmp(1,1)), 0; exp(d_mat_dmp(2,1)), exp(d_mat_dmp(2,2))];
    LLT_dmp = L_dmp * L_dmp';

    % Pre-computations
    for i = 1:N
        for j = 1:i-1
            %
            Gamma_km_ = lambda_isq * [ones(N,j),zeros(N,N-j)];
            zi_sj = LLT(i,:)' .* LLT(:,j);
            Mu_ij = 2 * diag_chi__min_chi * Gamma_km_;
            mu_ij = Mu_ij * zi_sj;
            Omega_ij = Mu_ij * diag(zi_sj) * (-Mu_ij') + 2 * diag( Gamma_km_ * zi_sj );
            %
            mu(:,:,i,j) = mu_ij;
            mu(:,:,j,i) = mu_ij;
            %
            Xi_Mu(i,j) = mu_ij' * xi;
            Xi_Mu(j,i) = Xi_Mu(i,j);
            %
            Alpha_Mu(i,j) = mu_ij' * alpha;
            Alpha_Mu(j,i) = Alpha_Mu(i,j);
            %
            Omega(:,:,i,j) = Omega_ij;
            Omega(:,:,j,i) = Omega_ij;
            %
            xiT_Omega(:,:,i,j) = xi' * Omega_ij;
            xiT_Omega(:,:,j,i) = xiT_Omega(:,:,i,j);
            Omega_xi(:,:,i,j) = xiT_Omega(:,:,i,j)';
            Omega_xi(:,:,j,i) = Omega_xi(:,:,i,j);
            %
            alphaT_Omega(:,:,i,j) = alpha' * Omega_ij;
            alphaT_Omega(:,:,j,i) = alphaT_Omega(:,:,i,j);
            %
            double_xi_ij = 2 * xi(i) * xi(j) * xi_(i) * xi_(j);
            dkf_dchi = dkf_dchi + double_xi_ij * mu_ij;
            %
            dchi__dchi_kf = dchi__dchi_kf + double_xi_ij * Omega_ij;
        end
        %
        Gamma_km_ = lambda_isq * [ones(N,i),zeros(N,N-i)];
        zi_si = LLT(i,:)' .* LLT(:,i);
        Mu_ii = 2 * diag_chi__min_chi * Gamma_km_;
        mu_ii = Mu_ii * zi_si;
        Omega_ii = Mu_ii * diag(zi_si) * (-Mu_ii') + 2 * diag(Gamma_km_ * zi_si);
        %
        mu(:,:,i,i) = mu_ii;
        %
        Xi_Mu(i,i) = mu_ii' * xi;
        %
        Alpha_Mu(i,i) = mu_ii' * alpha;
        %
        Omega(:,:,i,i) = Omega_ii;
        xiT_Omega(:,:,i,i) = xi' * Omega_ii;
        Omega_xi(:,:,i,i) = xiT_Omega(:,:,i,i)';
        %
        alphaT_Omega(:,:,i,i) = alpha' * Omega_ii;
        %
        xi_ii = xi(i) * xi(i) * xi_(i) * xi_(i);
        dkf_dchi = dkf_dchi + xi_ii * mu_ii;
        %
        dchi__dchi_kf = dchi__dchi_kf + xi_ii * Omega_ii;
    end
    xi_alpha_T = (xi_ .* alpha)';
    xi__z = xi_ .* z;
    xi_xi__T = (xi .* xi_)';
    for i = 1:N
        %
        Mu_i_T = squeeze(mu(:,:,:,i))';
        mat_xi__alpha_Mu_(i,:) = -xi_alpha_T * Mu_i_T;
        mat_xi__z_Mu_(i,:) = -xi__z' * Mu_i_T;
        %
        xiT_Omega_i = squeeze(xiT_Omega(:,:,:,i));
        mat_xi_xi__xiTOmega(i,:) = xi_xi__T * xiT_Omega_i;
        mat_xi__alpha_xiTOmega(i,:) = xi_alpha_T * xiT_Omega_i;
        %
        mat_xi_xi__alphaTOmega(i,:) = xi_xi__T * squeeze(alphaT_Omega(:,:,:,i));
        %
        mat_xi_xi__Mu(i,:) = xi_xi__T * Mu_i_T;
        %
        mat_Mu_xi__z(:,i) = Mu_i_T' * xi__z;
        %
        mat_Omegaxi_xi__xi(:,i) = squeeze(Omega_xi(:,:,i,:)) * xi_xi__T';
    end
    
    % M
    z_alphaT = z * alpha';
    d2_zT_kHess_alpha_dxi_2 = LLT .* (z_alphaT + z_alphaT');
    xi_alphaT = xi * alpha';
    xi_alpha_sym_sum = xi_alphaT + xi_alphaT';
    d2_xiT_d2kdchidxi_alpha_dxi_2 = Xi_Mu .* xi_alpha_sym_sum;
    xi_xiT = xi * xi';
    d2_nablachikT_alpha_dxi_2 = xi_xiT .* Alpha_Mu;
    % nabla_xi nabla_chi^T f
    d2_zT_kHess_alpha_dxi_dchi_ = diag(z) * mat_xi__alpha_Mu_ + diag(alpha) * mat_xi__z_Mu_;
    d2_xiT_d2kdchidxi_alpha_dxi_dchi_ = diag(alpha) * mat_xi_xi__xiTOmega + diag(xi) * mat_xi__alpha_xiTOmega;
    d2_nablachikT_alpha_dxi_dchi_ = diag(xi) * mat_xi_xi__alphaTOmega;
    % nabla_chi f
    nabla_chi__d2_kf_dxi2_times_zT = mat_Mu_xi__z * diag_xi_;
    nabla_chi__d2_kf_dxidchi_times_xiT = mat_Omegaxi_xi__xi * diag_xi_;
    % g
    chi_min_chi_T_Pginv2 = chi_min_chi_' * P_g_inv2;
    k_g_chichi_times_P_g_inv2 = rho_g_sq * exp(-1/2*chi_min_chi_T_Pginv2*chi_min_chi_) * P_g_inv2;
    % k_yf
    d2_kf_dxi2 = LLT .* xi__xi_T;
    d2_kf_dxidchi = diag_xi_ * mat_xi_xi__Mu;

    % M = nabla_xi nabla_xi^T f
    d2_cfTalpha_dxi_2 = d2_zT_kHess_alpha_dxi_2 + d2_xiT_d2kdchidxi_alpha_dxi_2 - d2_nablachikT_alpha_dxi_2;
    % nabla_xi nabla_chi^T f
    d2_cfTalpha_dxi_dchi_ = d2_zT_kHess_alpha_dxi_dchi_ + d2_xiT_d2kdchidxi_alpha_dxi_dchi_ - d2_nablachikT_alpha_dxi_dchi_;

    % update sums
    M_GP = M_GP + d2_cfTalpha_dxi_2;
    d2_f_GP_dxi_dchi_ = d2_f_GP_dxi_dchi_ + d2_cfTalpha_dxi_dchi_;
    D_GP = D_GP + LLT_dmp .* xi_alpha_sym_sum;

    % set vecs
    index_vec = 1+(d-1)*N:d*N;
    nabla_chi__k_yfT_pos_part(:,index_vec) = nabla_chi__d2_kf_dxi2_times_zT + nabla_chi__d2_kf_dxidchi_times_xiT;
    nabla_chi__k_yfT_neg_part(:,index_vec) = dchi__dchi_kf;
    nabla_chi__k_ygT(:,index_vec) = k_g_chichi_times_P_g_inv2 * (chi__min_chi*chi_min_chi_T_Pginv2 + eye(N));
    k_yf_pos_part_vec(index_vec) = d2_kf_dxi2 * z + d2_kf_dxidchi * xi;
    k_yf_neg_part_vec(index_vec) = dkf_dchi;
    k_yg(index_vec) = k_g_chichi_times_P_g_inv2 * chi__min_chi;
end

%% Compute estimates

% M = nabla_xi nabla_xi^T f
M_prior = hatM(chi_);
M_est = M_prior + 1/2 * M_GP;

% nabla_xi nabla_chi^T f
d2_f_est = hat_M_dot(chi_,xi_) + 1/2 * d2_f_GP_dxi_dchi_;

% nabla_chi f
nabla_chi__k_yfT = 1/2 * nabla_chi__k_yfT_pos_part - 1/4 * nabla_chi__k_yfT_neg_part;
nabla_chi_f_est = hat_nabla_chi_f_GP(chi_,xi_) + nabla_chi__k_yfT * a;

% nabla_chi g
chi_T_Pginv2 = chi_' * P_g_inv2;
k_g_0chi_ = rho_g_sq * exp( -1/2 * chi_T_Pginv2 * chi_ );
nabla_chi__k_0gT = -k_g_0chi_ * P_g_inv2 * [chi_, chi_*chi_T_Pginv2 - eye(N)];
nabla_chi_g_est = hat_nabla_chi_g_GP(chi_) + nabla_chi__k_ygT * a - nabla_chi__k_0gT * s;

% f
k_yf = 1/2 * k_yf_pos_part_vec - 1/4 * k_yf_neg_part_vec;
f_est = 1/2 * xi_' * M_prior * xi_ + k_yf' * a;

% g
k_0gT = k_g_0chi_ * [1, chi_T_Pginv2];
g_est = hatV(chi_) + k_yg' * a - k_0gT * s;

% D
D_est = hatD(xi_) + 1/2 * D_GP;

% P_dsp
P_dsp_est = xi_' * D_est * xi_;

% tau
tau_est = M_est * ddq_ + d2_f_est * xi_ - nabla_chi_f_est + nabla_chi_g_est + D_est * xi_;

% ddq
ddq_est = M_est \ (tau_ - d2_f_est * xi_ + nabla_chi_f_est - nabla_chi_g_est - D_est * xi_);

end