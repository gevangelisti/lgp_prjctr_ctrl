function [tau,M_est,C_est,D_est,g_est,Sigma_lgp,g_est_e,D_est_de,e,de] = eval_pcc_lgp_var_nat_pdp(t,chi_,xi_,a_des,omega_d,...
    K1,K2,K3,Kd,K_y_plus_Pi,...
    X, ddq_meas, LLT_Sigma_sq, lambda_isq, rho_g_sq, P_g_inv2, a, s, pb, D, N)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

q_des = @(t) a_des.*sin(omega_d*t);
dq_des = @(t) omega_d*a_des.*cos(omega_d*t);
ddq_des = @(t) -omega_d^2*a_des.*sin(omega_d*t);

[M_prior,C_prior,g_prior] = MCg_pcc_n4(pb.mu,pb.L,pb.Iz,pb.g,chi_,xi_);
e = chi_-q_des(t);
[~,~,g_prior_e] = MCg_pcc_n4(pb.mu,pb.L,pb.Iz,pb.g,e,xi_);

%% Lagrange GP computations

% Initializations
M_GP = zeros(N);
d2_f_GP_dxi_dchi_ = zeros(N);
sumk_xik_dchik_d2_f_GP_dxi_2 = zeros(N);
nabla_chi__k_ygT = zeros(N,N*D);
nabla_e_k_ygT = zeros(N,N*D);

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
Xi__Mu_ = zeros(N);
xiT_Omega_xi_ = zeros(N);
alphaT_Omega_xi_ = zeros(N);

for d = 1:D
    % Set input/state values
    chi = X(d,1:N)';
    xi = X(d,N+1:2*N)';
    alpha = a(1+(d-1)*N:d*N);
    z = ddq_meas(d,:)';

    % 
    chi_min_chi_ = chi-chi_;
    chi_min_e = chi-e;
    chi_p_chi_ = chi+chi_;
    LLT = exp(-lambda_isq * (chi_min_chi_' * chi_min_chi_)) * LLT_Sigma_sq;
    LLT_n = exp(-lambda_isq * (chi_p_chi_' * chi_p_chi_)) * LLT_Sigma_sq;
    chi__min_chi = -chi_min_chi_;
    diag_chi__min_chi = diag(chi__min_chi);
    diag_chi__p_chi = diag(chi_p_chi_);

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
            Xi__Mu_(i,j) = -mu_ij_sum' * xi_;
            Xi__Mu_(j,i) = Xi__Mu_(i,j);
            %
            xiT_Omega_xi_(i,j) = xiT_Omega(:,:,i,j) * xi_;
            xiT_Omega_xi_(j,i) = xiT_Omega_xi_(i,j);
            %
            alphaT_Omega_xi_(i,j) = alphaT_Omega(:,:,i,j) * xi_;
            alphaT_Omega_xi_(j,i) = alphaT_Omega_xi_(i,j);
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
        Xi__Mu_(i,i) = -mu_ii_sum' * xi_;
        %
        xiT_Omega_xi_(i,i) = xiT_Omega(:,:,i,i) * xi_;
        %
        alphaT_Omega_xi_(i,i) = alphaT_Omega(:,:,i,i) * xi_;
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
    z_alphaT_sumT = z_alphaT + z_alphaT';
    LLT_sym = LLT + LLT_n;
    d2_zT_kHess_alpha_dxi_2 = LLT_sym .* z_alphaT_sumT;
    xi_alphaT = xi * alpha';
    xi_alphaT_sumT = xi_alphaT + xi_alphaT';
    d2_xiT_d2kdchidxi_alpha_dxi_2 = Xi_Mu .* xi_alphaT_sumT;
    xi_xiT = xi * xi';
    d2_nablachikT_alpha_dxi_2 = xi_xiT .* Alpha_Mu;
    % nabla_xi nabla_chi^T f
    d2_zT_kHess_alpha_dxi_dchi_ = diag(z) * mat_xi__alpha_Mu_ + diag(alpha) * mat_xi__z_Mu_;
    d2_xiT_d2kdchidxi_alpha_dxi_dchi_ = diag(alpha) * mat_xi_xi__xiTOmega + diag(xi) * mat_xi__alpha_xiTOmega;
    d2_nablachikT_alpha_dxi_dchi_ = diag(xi) * mat_xi_xi__alphaTOmega;
    sumk_xik_dchik_d2_zT_kHess_alpha_dxi_2 = Xi__Mu_ .* z_alphaT_sumT;
    sumk_xik_dchik_d2_xiT_d2kdchidxi_alpha_dxi_2 = xiT_Omega_xi_ .* xi_alphaT_sumT;
    sumk_xik_dchik_d2_nablachikT_alpha_dxi_2 = xi_xiT .* alphaT_Omega_xi_;
    % g
    chisT_Pginv2 = chi_min_chi_' * P_g_inv2;
    chi_p_chi_ = chi+chi_;
    chis_pl_T_Pginv2 = chi_p_chi_' * P_g_inv2;
    nabla_chi__nabla_chi_k_g = rho_g_sq * ( ...
        exp( -1/2 * chisT_Pginv2 * chi_min_chi_ ) * P_g_inv2 * (chi__min_chi*chisT_Pginv2 + eye(N)) + ...
        exp( -1/2 * chis_pl_T_Pginv2 * chi_p_chi_ ) * P_g_inv2 * (chi_p_chi_*chis_pl_T_Pginv2 - eye(N)) ...
        );
    % g(e)
    chisT_Pginv2 = chi_min_e' * P_g_inv2;
    chi_p_e = chi+e;
    chis_pl_T_Pginv2 = chi_p_e' * P_g_inv2;
    nabla_e_nabla_chi_k_g = rho_g_sq * ( ...
        exp( -1/2 * chisT_Pginv2 * chi_min_e ) * P_g_inv2 * (-chi_min_e*chisT_Pginv2 + eye(N)) + ...
        exp( -1/2 * chis_pl_T_Pginv2 * chi_p_e ) * P_g_inv2 * (chi_p_e*chis_pl_T_Pginv2 - eye(N)) ...
        );

    % M = nabla_xi nabla_xi^T f
    d2_cfTalpha_dxi_2 = d2_zT_kHess_alpha_dxi_2 + d2_xiT_d2kdchidxi_alpha_dxi_2 - d2_nablachikT_alpha_dxi_2;
    % nabla_xi nabla_chi^T f
    d2_cfTalpha_dxi_dchi_ = d2_zT_kHess_alpha_dxi_dchi_ + d2_xiT_d2kdchidxi_alpha_dxi_dchi_ - d2_nablachikT_alpha_dxi_dchi_;
    sumk_xik_dchik_d2_cfTalpha_dxi_2 = sumk_xik_dchik_d2_zT_kHess_alpha_dxi_2 + sumk_xik_dchik_d2_xiT_d2kdchidxi_alpha_dxi_2 - sumk_xik_dchik_d2_nablachikT_alpha_dxi_2;
    
    % update sums
    M_GP = M_GP + d2_cfTalpha_dxi_2;
    d2_f_GP_dxi_dchi_ = d2_f_GP_dxi_dchi_ + d2_cfTalpha_dxi_dchi_;
    sumk_xik_dchik_d2_f_GP_dxi_2 = sumk_xik_dchik_d2_f_GP_dxi_2 + sumk_xik_dchik_d2_cfTalpha_dxi_2;

    % set vecs
    index_vec = 1+(d-1)*N:d*N;
    nabla_chi__k_ygT(:,index_vec) = nabla_chi__nabla_chi_k_g;
    nabla_e_k_ygT(:,index_vec) = nabla_e_nabla_chi_k_g;
end

%% Compute estimates

% M = nabla_xi nabla_xi^T f
M_est = M_prior + 1/2 * M_GP;

% nabla_xi nabla_chi^T f
C_est = C_prior + 1/4 * (d2_f_GP_dxi_dchi_ + sumk_xik_dchik_d2_f_GP_dxi_2 - d2_f_GP_dxi_dchi_');

% nabla_chi g
k_g_0chi_ = rho_g_sq * exp( -1/2 * chi_' * P_g_inv2 * chi_ );
nabla_chi__k_0gT = -2 * k_g_0chi_ * P_g_inv2 * chi_;
g_est = g_prior + nabla_chi__k_ygT * a - nabla_chi__k_0gT * s;

% g(e)
k_g_0e = rho_g_sq * exp( -1/2 * e' * P_g_inv2 * e );
nabla_e_k_0gT = -2 * k_g_0e * P_g_inv2 * e;
g_est_e = g_prior_e + nabla_e_k_ygT * a - nabla_e_k_0gT * s;
g_est = g_est + diag(pb.k)*chi_;
g_est_e = g_est_e + diag(pb.k)*e;

d = diag(pb.d)*xi_;
D_est = diag(pb.d);

de_ = xi_-dq_des(t);
de = diag(pb.d)*de_;
D_est_de = diag(pb.d);

%%
Sigma_lgp_temp = calc_lgp_var(chi_, xi_, zeros(N,1), X, ddq_meas, LLT_Sigma_sq, lambda_isq, rho_g_sq, P_g_inv2, D, N, K_y_plus_Pi);
tau_ = M_est*ddq_des(t) + C_est*dq_des(t) ... 
    + (eye(N)-heaviside_man(e'*g_est)*(e*e')./(1e-3 + norm(e)^2) )*g_est ...
    - g_est_e - K_vfb(K1,K2,K3,Sigma_lgp_temp,e) ...
    + (eye(N)-heaviside_man(de_'*d)*(de_*de_')./(1e-3 + norm(de_)^2) )*d ...
    - de - Kd*de_ - K_vfb(K1,K2,K3,Sigma_lgp_temp,de);

ddq_ = M_est \ (tau_ - C_est * xi_ - g_est - diag(pb.d)*xi_);

%%
Sigma_lgp = calc_lgp_var(chi_, xi_, ddq_, X, ddq_meas, LLT_Sigma_sq, lambda_isq, rho_g_sq, P_g_inv2, D, N, K_y_plus_Pi);

tau = M_est*ddq_des(t) + C_est*dq_des(t) ... 
    + (eye(N)-heaviside_man(e'*g_est)*(e*e')./(1e-3 + norm(e)^2) )*g_est ...
    - g_est_e - K_vfb(K1,K2,K3,Sigma_lgp,e) ... 
    + (eye(N)-heaviside_man(de_'*d)*(de_*de_')./(1e-3 + norm(de_)^2) )*d ...
    - de - Kd*de_ - K_vfb(K1,K2,K3,Sigma_lgp,de);

end

function [x_amp] = K_vfb(K1,K2,K3,Sigma,x)

x_amp = K1 * (x - ( ( K3*(K2+Sigma)*K3 + K1 )\(K1*x) ) );

end

function [y] = heaviside_man(x)

if x > 0
    y = 1;
elseif x < 0
    y = 0;
else
    y = 1/2;
end

end

function [tau_gp_var_est] = calc_lgp_var(q, dq, ddq, X_tr, z_meas, LLT_Sigma_sq, lambda_isq, rho_g_sq, P_g_inv2, D, N, K)

% Var
var_mat = Cov_mat_func([q;dq], [q;dq], ddq, ddq, N, lambda_isq, rho_g_sq, P_g_inv2, LLT_Sigma_sq);

% Cov
K_sum = zeros(D,N,N);
parfor j = 1:D
    K_sum(j,:,:) = Cov_mat_func([q;dq], X_tr(j,:)', ddq, z_meas(j,:)', N, lambda_isq, rho_g_sq, P_g_inv2, LLT_Sigma_sq);
end
cov_mat = zeros(N,N*D);
for j = 1:D
    cov_mat(:,1+(j-1)*N:j*N) = K_sum(j,:,:);
end

% Full
tau_gp_var_est = var_mat - cov_mat*(K\cov_mat');

end

function [Cov_yy_] = Cov_mat_func(x, x_, z, z_, N, lambda_isq, rho_g_sq, P_g_inv2, LLT_Sigma_sq)

% Set input/state values
chi = x(1:N);
xi = x(N+1:2*N);
chi_ = x_(1:N);
xi_ = x_(N+1:2*N);

% Initializations
chi_min_chi_ = chi-chi_;
chi__min_chi = -chi_min_chi_;
diag_chi__min_chi = diag(chi__min_chi);
LLT = exp(-lambda_isq * (chi_min_chi_' * chi_min_chi_)) * LLT_Sigma_sq;
chi_p_chi_ = chi+chi_;
diag_chi__p_chi = diag(chi_p_chi_);
LLT_n = exp(-lambda_isq * (chi_p_chi_' * chi_p_chi_)) * LLT_Sigma_sq;
mu = zeros(N,1,N,N);
mu_n = zeros(N,1,N,N);
V__Mu_ = zeros(N);
diag_z_xi__Mu__v_ = zeros(N,1);
V_Mu = zeros(N);
diag_v_Mu_xi_z_ = zeros(N,1);
v_Omega_v_ = zeros(N);
Omega_v_ = zeros(N,1,N,N);
Omega_v = zeros(N,1,N,N);
diag_v_Omegav__xi_xi_ = zeros(N,1);
mat_z_xi__Mu_ = zeros(N);
mat_Omegav__xi_xi_ = zeros(N);
mat_xi_xi__Omegav = zeros(N);
mat_Mu_xi_z_ = zeros(N);
dchi_dchi__kf = zeros(N);

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
        % 2
        mu(:,:,i,j) = mu_ij;
        mu_n(:,:,i,j) = mu_ij_n;
        mu(:,:,j,i) = mu_ij;
        mu(:,:,j,i) = mu_ij_n;
        V__Mu_(i,j) = (-mu_ij+mu_ij_n)' * xi_;
        V__Mu_(j,i) = V__Mu_(i,j);
        % 3
        V_Mu(i,j) = (mu_ij+mu_ij_n)' * xi;
        V_Mu(j,i) = V_Mu(i,j);
        % 4
        Omega_ij_sum = Omega_ij + Omega_ij_n;
        Omega_v_(:,:,i,j) = Omega_ij_sum * xi_;
        Omega_v_(:,:,j,i) = Omega_v_(:,:,i,j);
        v_Omega_v_(i,j) = xi' * Omega_v_(:,:,i,j);
        v_Omega_v_(j,i) = v_Omega_v_(i,j);
        % 6
        Omega_v(:,:,i,j) = Omega_ij_sum * xi;
        Omega_v(:,:,j,i) = Omega_v(:,:,i,j);
        % 9
        dchi_dchi__kf = dchi_dchi__kf + 2 * xi(i) * xi(j) * xi_(i) * xi_(j) * Omega_ij_sum;
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
    % 2
    mu(:,:,i,i) = mu_ii + mu_ii_n;
    mu_n(:,:,i,i) = mu_ii_n;
    V__Mu_(i,i) = (-mu_ii+mu_ii_n)' * xi_;
    % 3
    V_Mu(i,i) = (mu_ii+mu_ii_n)' * xi;
    % 4
    Omega_ii_sum = Omega_ii + Omega_ii_n;
    Omega_v_(:,:,i,i) = Omega_ii_sum * xi_;
    v_Omega_v_(i,i) = xi' * Omega_v_(:,:,i,i);
    % 6
    Omega_v(:,:,i,i) = Omega_ii_sum * xi;
    % 9
    dchi_dchi__kf = dchi_dchi__kf + xi(i) * xi(i) * xi_(i) * xi_(i) * Omega_ii_sum;
end
z_xi__T = (z .* xi_)';
xi_z_ = xi .* z_;
xi_xi_ = xi .* xi_;
for i = 1:N
    %
    mu_i = squeeze(mu(:,:,:,i));
    mu_i_n = squeeze(mu_n(:,:,:,i));
    % 5
    mat_z_xi__Mu_(i,:) = z_xi__T * (-mu_i+mu_i_n)';
    % 2
    diag_z_xi__Mu__v_(i) = mat_z_xi__Mu_(i,:) * xi_;
    % 7
    mat_Mu_xi_z_(:,i) = (mu_i+mu_i_n) * xi_z_;
    % 3
    diag_v_Mu_xi_z_(i) = xi' * mat_Mu_xi_z_(:,i);
    % 8
    mat_Omegav__xi_xi_(:,i) = squeeze(Omega_v_(:,:,:,i)) * xi_xi_;
    % 4
    diag_v_Omegav__xi_xi_(i) = xi' * mat_Omegav__xi_xi_(:,i);
    % 6
    mat_xi_xi__Omegav(i,:) = xi_xi_' * squeeze(Omega_v(:,:,i,:))';
end
diag_z_xi__Mu__v_ = diag(diag_z_xi__Mu__v_);
diag_v_Mu_xi_z_ = diag(diag_v_Mu_xi_z_);
diag_v_Omegav__xi_xi_ = diag(diag_v_Omegav__xi_xi_);

% 1 - nablas: xi, xi', xi, xi'
LLT_sym = LLT + LLT_n;
C_double = LLT_sym .* (z_*z') + diag(LLT_sym * (z.*z_));

% 2 - nablas: xi, xi', xi, chi'
C_double = C_double + (xi_*z').*V__Mu_ + diag_z_xi__Mu__v_;

% 3 - nablas: xi, xi', chi, xi'
C_double = C_double + (z_*xi').*V_Mu + diag_v_Mu_xi_z_;

% 4 - nablas: xi, xi', chi, chi'
C_double = C_double + (xi_*xi').*v_Omega_v_ + diag_v_Omegav__xi_xi_;

% 5 - nablas: xi, xi, chi'
diag_xi_ = diag(xi_);
C_double = C_double - diag_xi_ * mat_z_xi__Mu_;

% 6 - nablas: xi, chi, chi'
C_double = C_double - diag_xi_ * mat_xi_xi__Omegav;

% 7 - nablas: chi, xi', xi'
diag_xi = diag(xi);
C_double = C_double - mat_Mu_xi_z_ * diag_xi;

% 8 - nablas: chi, chi', xi'
C_double = C_double - mat_Omegav__xi_xi_ * diag_xi;

% 9 - nablas: chi, chi'
C = 1/2 * C_double + 1/4 * dchi_dchi__kf;

% k_g
chisT_Pginv2 = chi_min_chi_' * P_g_inv2;
chi_p_chi_ = chi+chi_;
chis_pl_T_Pginv2 = chi_p_chi_' * P_g_inv2;
C = C + rho_g_sq * ( ...
    exp( -1/2 * chisT_Pginv2 * chi_min_chi_ ) * P_g_inv2 * (chi__min_chi*chisT_Pginv2 + eye(N)) + ...
    exp( -1/2 * chis_pl_T_Pginv2 * chi_p_chi_ ) * P_g_inv2 * (chi_p_chi_*chis_pl_T_Pginv2 - eye(N)) ...
    );

% Final sum
Cov_yy_ = C;

end