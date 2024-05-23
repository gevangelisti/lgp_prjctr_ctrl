function [ddq_est] = ...
    ddq_lgp_pcc(...
    chi_, xi_, tau_, X, ddq_meas, LLT_Sigma_sq, lambda_isq, rho_g_sq, P_g_inv2, a, s, pb, D, N ...
    )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Parametric model estimates
[M_prior,C_prior,g_prior] = MCg_pcc_n4(pb.mu,pb.L,pb.Iz,pb.g,chi_,xi_);

%% Lagrange GP computations

% Initializations
M_GP = zeros(N);
d2_f_GP_dxi_dchi_ = zeros(N);
nabla_chi__k_yfT_pos_part = zeros(N,N*D);
nabla_chi__k_yfT_neg_part = zeros(N,N*D);
nabla_chi__k_ygT = zeros(N,N*D);

% Initializations
mu = zeros(N,1,N,N);
mu_n = zeros(N,1,N,N);
Xi_Mu = zeros(N);
Alpha_Mu = zeros(N);
mat_xi__alpha_Mu_ = zeros(N);
mat_xi__z_Mu_ = zeros(N);
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
diag_xi_ = diag(xi_);

for d = 1:D
    % Set input/state values
    chi = X(d,1:N)';
    xi = X(d,N+1:2*N)';
    alpha = a(1+(d-1)*N:d*N);
    z = ddq_meas(d,:)';

    % 
    chi_min_chi_ = chi-chi_;
    chi_p_chi_ = chi+chi_;
    LLT = exp(-lambda_isq * (chi_min_chi_' * chi_min_chi_)) * LLT_Sigma_sq;
    LLT_n = exp(-lambda_isq * (chi_p_chi_' * chi_p_chi_)) * LLT_Sigma_sq;
    chi__min_chi = -chi_min_chi_;
    diag_chi__min_chi = diag(chi__min_chi);
    diag_chi__p_chi = diag(chi_p_chi_);
    dchi__dchi_kf = zeros(N);
    dkf_dchi = zeros(N,1);

    % Pre-computations
    for i = 1:N
        for j = 1:i-1
            %
            Gamma_km_ = lambda_isq * [ones(N,j),zeros(N,N-j)];
            zi_sj = LLT(i,:)' .* LLT(:,j);
            zi_sj_n = LLT_n(i,:)' .* LLT_n(:,j);
            Mu_ij = 2 * diag_chi__min_chi * Gamma_km_;
            Mu_ij_n = -2 * diag_chi__p_chi * Gamma_km_;
            mu_ij = Mu_ij * zi_sj;
            mu_ij_n = Mu_ij_n * zi_sj_n;
            Omega_ij = Mu_ij * diag(zi_sj) * (-Mu_ij') + 2 * diag( Gamma_km_ * zi_sj );
            Omega_ij_n = Mu_ij_n * diag(zi_sj_n) * Mu_ij_n' - 2 * diag( Gamma_km_ * zi_sj_n );
            %
            mu(:,:,i,j) = mu_ij;
            mu_n(:,:,i,j) = mu_ij_n;
            mu(:,:,j,i) = mu_ij;
            mu(:,:,j,i) = mu_ij_n;
            %
            mu_ij_sum = mu_ij+mu_ij_n;
            Xi_Mu(i,j) = mu_ij_sum' * xi;
            Xi_Mu(j,i) = Xi_Mu(i,j);
            %
            Alpha_Mu(i,j) = mu_ij_sum' * alpha;
            Alpha_Mu(j,i) = Alpha_Mu(i,j);
            %
            Omega_ij_sum = Omega_ij + Omega_ij_n;
            xiT_Omega(:,:,i,j) = xi' * Omega_ij_sum;
            xiT_Omega(:,:,j,i) = xiT_Omega(:,:,i,j);
            Omega_xi(:,:,i,j) = xiT_Omega(:,:,i,j)';
            Omega_xi(:,:,j,i) = Omega_xi(:,:,i,j);
            %
            alphaT_Omega(:,:,i,j) = alpha' * Omega_ij_sum;
            alphaT_Omega(:,:,j,i) = alphaT_Omega(:,:,i,j);
            %
            double_xi_ij = 2 * xi(i) * xi(j) * xi_(i) * xi_(j);
            dkf_dchi = dkf_dchi + double_xi_ij * mu_ij_sum;
            %
            dchi__dchi_kf = dchi__dchi_kf + double_xi_ij * Omega_ij_sum;
        end
        %
        Gamma_km_ = lambda_isq * [ones(N,i),zeros(N,N-i)];
        zi_si = LLT(i,:)' .* LLT(:,i);
        zi_si_n = LLT_n(i,:)' .* LLT_n(:,i);
        Mu_ii = 2 * diag_chi__min_chi * Gamma_km_;
        Mu_ii_n = -2 * diag_chi__p_chi * Gamma_km_;
        mu_ii = Mu_ii * zi_si;
        mu_ii_n = Mu_ii_n * zi_si_n;
        Omega_ii = Mu_ii * diag(zi_si) * (-Mu_ii') + 2 * diag(Gamma_km_ * zi_si);
        Omega_ii_n = Mu_ii_n * diag(zi_si_n) * Mu_ii_n' - 2 * diag( Gamma_km_ * zi_si_n );
        %
        mu(:,:,i,i) = mu_ii;
        mu_n(:,:,i,i) = mu_ii_n;
        %
        mu_ii_sum = mu_ii + mu_ii_n;
        Xi_Mu(i,i) = mu_ii_sum' * xi;
        %
        Alpha_Mu(i,i) = mu_ii_sum' * alpha;
        %
        Omega_ii_sum = Omega_ii + Omega_ii_n;
        xiT_Omega(:,:,i,i) = xi' * Omega_ii_sum;
        Omega_xi(:,:,i,i) = xiT_Omega(:,:,i,i)';
        %
        alphaT_Omega(:,:,i,i) = alpha' * Omega_ii_sum;
        %
        xi_ii = xi(i) * xi(i) * xi_(i) * xi_(i);
        dkf_dchi = dkf_dchi + xi_ii * mu_ii_sum;
        %
        dchi__dchi_kf = dchi__dchi_kf + xi_ii * Omega_ii_sum;
    end
    xi_alpha_T = (xi_ .* alpha)';
    xi__z = xi_ .* z;
    xi_xi__T = (xi .* xi_)';
    for i = 1:N
        %
        Mu_i_T = squeeze(mu(:,:,:,i))';
        Mu_i_T_n = squeeze(mu_n(:,:,:,i))';
        mat_xi__alpha_Mu_(i,:) = xi_alpha_T * (-Mu_i_T + Mu_i_T_n);
        mat_xi__z_Mu_(i,:) = xi__z' * (-Mu_i_T + Mu_i_T_n);
        %
        xiT_Omega_i = squeeze(xiT_Omega(:,:,:,i));
        mat_xi_xi__xiTOmega(i,:) = xi_xi__T * xiT_Omega_i;
        mat_xi__alpha_xiTOmega(i,:) = xi_alpha_T * xiT_Omega_i;
        %
        mat_xi_xi__alphaTOmega(i,:) = xi_xi__T * squeeze(alphaT_Omega(:,:,:,i));
        %
        mat_xi_xi__Mu(i,:) = xi_xi__T * (Mu_i_T + Mu_i_T_n);
        %
        mat_Mu_xi__z(:,i) = (Mu_i_T + Mu_i_T_n)' * xi__z;
        %
        mat_Omegaxi_xi__xi(:,i) = squeeze(Omega_xi(:,:,i,:)) * xi_xi__T';
    end
    
    % M
    z_alphaT = z * alpha';
    LLT_sym = LLT + LLT_n;
    d2_zT_kHess_alpha_dxi_2 = LLT_sym .* (z_alphaT + z_alphaT');
    xi_alphaT = xi * alpha';
    d2_xiT_d2kdchidxi_alpha_dxi_2 = Xi_Mu .* (xi_alphaT + xi_alphaT');
    d2_nablachikT_alpha_dxi_2 = (xi * xi') .* Alpha_Mu;
    % nabla_xi nabla_chi^T f
    d2_zT_kHess_alpha_dxi_dchi_ = diag(z) * mat_xi__alpha_Mu_ + diag(alpha) * mat_xi__z_Mu_;
    d2_xiT_d2kdchidxi_alpha_dxi_dchi_ = diag(alpha) * mat_xi_xi__xiTOmega + diag(xi) * mat_xi__alpha_xiTOmega;
    d2_nablachikT_alpha_dxi_dchi_ = diag(xi) * mat_xi_xi__alphaTOmega;
    % nabla_chi f
    nabla_chi__d2_kf_dxi2_times_zT = mat_Mu_xi__z * diag_xi_;
    nabla_chi__d2_kf_dxidchi_times_xiT = mat_Omegaxi_xi__xi * diag_xi_;
    % g
    chisT_Pginv2 = chi_min_chi_' * P_g_inv2;
    chi_p_chi_ = chi+chi_;
    chis_pl_T_Pginv2 = chi_p_chi_' * P_g_inv2;
    nabla_chi__nabla_chi_k_g = rho_g_sq * ( ...
        exp( -1/2 * chisT_Pginv2 * chi_min_chi_ ) * P_g_inv2 * (chi__min_chi*chisT_Pginv2 + eye(N)) + ...
        exp( -1/2 * chis_pl_T_Pginv2 * chi_p_chi_ ) * P_g_inv2 * (chi_p_chi_*chis_pl_T_Pginv2 - eye(N)) ...
        );

    % M = nabla_xi nabla_xi^T f
    d2_cfTalpha_dxi_2 = d2_zT_kHess_alpha_dxi_2 + d2_xiT_d2kdchidxi_alpha_dxi_2 - d2_nablachikT_alpha_dxi_2;
    % nabla_xi nabla_chi^T f
    d2_cfTalpha_dxi_dchi_ = d2_zT_kHess_alpha_dxi_dchi_ + d2_xiT_d2kdchidxi_alpha_dxi_dchi_ - d2_nablachikT_alpha_dxi_dchi_;

    % update sums
    M_GP = M_GP + d2_cfTalpha_dxi_2;
    d2_f_GP_dxi_dchi_ = d2_f_GP_dxi_dchi_ + d2_cfTalpha_dxi_dchi_;

    % set vecs
    index_vec = 1+(d-1)*N:d*N;
    nabla_chi__k_yfT_pos_part(:,index_vec) = nabla_chi__d2_kf_dxi2_times_zT + nabla_chi__d2_kf_dxidchi_times_xiT;
    nabla_chi__k_yfT_neg_part(:,index_vec) = dchi__dchi_kf;
    nabla_chi__k_ygT(:,index_vec) = nabla_chi__nabla_chi_k_g;
end

%% Compute estimates

% M = nabla_xi nabla_xi^T f
M_est = M_prior + 1/2 * M_GP;

% nabla_xi nabla_chi^T f
d2_f_est = C_prior + 1/2 * d2_f_GP_dxi_dchi_;

% nabla_chi f
nabla_chi__k_yfT = 1/2 * nabla_chi__k_yfT_pos_part - 1/4 * nabla_chi__k_yfT_neg_part;
nabla_chi_f_est = nabla_chi__k_yfT * a;

% nabla_chi g
k_g_0chi_ = rho_g_sq * exp( -1/2 * chi_' * P_g_inv2 * chi_ );
nabla_chi__k_0gT = -2 * k_g_0chi_ * P_g_inv2 * chi_;
nabla_chi_g_est = g_prior + nabla_chi__k_ygT * a - nabla_chi__k_0gT * s;

% ddq
ddq_est = M_est \ (tau_ - d2_f_est * xi_ + nabla_chi_f_est - nabla_chi_g_est - diag(pb.d)*xi_ - diag(pb.k)*chi_);

end