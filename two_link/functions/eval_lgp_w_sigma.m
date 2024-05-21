function [M_est,C_est,D_est,g_est,G_est,G_est_e,Sigma_lgp,g_est_e,D_est_de,e,de] = eval_lgp_w_sigma(t_, chi_, xi_, ...
    a_d, omega_d, K1, K2, K3, Kp, Kd, K_y_plus_Pi, ...
    X_tr, ddq_meas, ...
    LLT_Sigma_sq, lambda_isq, rho_g_sq, P_g_inv2, L_Sigma_sq_dmp, lambda_isq_dmp, ...
    a, s, pb, D, N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

q_d = @(t) a_d*sin(omega_d*t)*ones(N,1);
dq_d = @(t) a_d*omega_d*cos(omega_d*t)*ones(N,1);
ddq_d = @(t) -omega_d^2*a_d*sin(omega_d*t)*ones(N,1);

%% Parametric model estimates
hatM = @(q) ...
    [pb.alphah + 2*pb.betah*cos(q(2)), ...
    pb.deltah + pb.betah*cos(q(2)); ...
    pb.deltah + pb.betah*cos(q(2)), ...
    pb.deltah];
hatC = @(q,dq) ...
    [-pb.betah*sin(q(2))*dq(2), ...
    -pb.betah*sin(q(2))*(dq(1)+dq(2)); ...
    pb.betah*sin(q(2))*dq(1), ...
    0];
hat_nabla_chi_g_GP = @(q) ...
    [(pb.m1h*pb.l1h/2 + pb.m2h*pb.l1h)*pb.gh*sin(q(1)) + pb.m2h*pb.l2h/2*pb.gh*sin(q(1)+q(2)); ...
    pb.m2h*pb.l2h/2*pb.gh*sin(q(1)+q(2))];
hatD = @(dq) [pb.dh(1)+pb.dh(2)*abs(dq(1)), 0; 0, pb.dh(1)+pb.dh(2)*abs(dq(2))];
hatV = @(q) ...
    (pb.m1h*pb.l1h/2 + pb.m2h*pb.l1h)*pb.gh*(1-cos(q(1))) + pb.m2h*pb.l2h/2*pb.gh*(1-cos(q(1)+q(2)));

%% Lagrange GP computations
e = chi_-q_d(t_);
de = xi_-dq_d(t_);

% Initializations
M_GP = zeros(N);
d2_f_GP_dxi_dchi_ = zeros(N);
sumk_xik_dchik_d2_f_GP_dxi_2 = zeros(N);
nabla_chi__k_ygT = zeros(N,N*D);
D_GP = zeros(N);
nabla_chi__k_ygT_e = zeros(N,N*D);
D_GP_de = zeros(N);
k_yg = zeros(N*D,1);
k_yg_e = zeros(N*D,1);

% Initializations
mu = zeros(N,1,N,N);
Xi_Mu = zeros(N);
Alpha_Mu = zeros(N);
mat_xi__alpha_Mu_ = zeros(N);
mat_xi__z_Mu_ = zeros(N);
Omega = zeros(N,N,N,N);
xiT_Omega = zeros(1,N,N,N);
mat_xi_xi__xiTOmega = zeros(N);
mat_xi__alpha_xiTOmega = zeros(N);
alphaT_Omega = zeros(1,N,N,N);
mat_xi_xi__alphaTOmega = zeros(N);
Xi__Mu_ = zeros(N);
xiT_Omega_xi_ = zeros(N);
alphaT_Omega_xi_ = zeros(N);

for d = 1:D
    % Set input/state values
    chi = X_tr(d,1:N)';
    xi = X_tr(d,N+1:2*N)';
    alpha = a(1+(d-1)*N:d*N);
    z = ddq_meas(d,:)';

    % 
    chi_min_chi_ = chi-chi_;
    chi_min_e = chi-e;
    LLT = exp(-lambda_isq * (chi_min_chi_' * chi_min_chi_)) * LLT_Sigma_sq;
    chi__min_chi = -chi_min_chi_;
    diag_chi__min_chi = diag(chi__min_chi);
    %
    xi_min_xi_ = xi-xi_;
    d_mat_dmp = -1/2*lambda_isq_dmp .* [xi_min_xi_(1)^2, 0; xi_min_xi_'*xi_min_xi_, xi_min_xi_(2)^2];
    L_dmp = L_Sigma_sq_dmp .* [exp(d_mat_dmp(1,1)), 0; exp(d_mat_dmp(2,1)), exp(d_mat_dmp(2,2))];
    LLT_dmp = L_dmp * L_dmp';
    %
    xi_min_de = xi-de;
    d_mat_dmp_de = -1/2*lambda_isq_dmp .* [xi_min_de(1)^2, 0; xi_min_de'*xi_min_de, xi_min_de(2)^2];
    L_dmp_de = L_Sigma_sq_dmp .* [exp(d_mat_dmp_de(1,1)), 0; exp(d_mat_dmp_de(2,1)), exp(d_mat_dmp_de(2,2))];
    LLT_dmp_de = L_dmp_de * L_dmp_de';

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
            %
            alphaT_Omega(:,:,i,j) = alpha' * Omega_ij;
            alphaT_Omega(:,:,j,i) = alphaT_Omega(:,:,i,j);
            %
            Xi__Mu_(i,j) = -mu_ij' * xi_;
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
        %
        alphaT_Omega(:,:,i,i) = alpha' * Omega_ii;
        %
        Xi__Mu_(i,i) = -mu_ii' * xi_;
        %
        xiT_Omega_xi_(i,i) = xiT_Omega(:,:,i,i) * xi_;
        %
        alphaT_Omega_xi_(i,i) = alphaT_Omega(:,:,i,i) * xi_;
    end
    xi_alpha_T = (xi_ .* alpha)';
    xi__z_T = (xi_ .* z)';
    xi_xi__T = (xi .* xi_)';
    for i = 1:N
        %
        neg_mu_i_T = -squeeze(mu(:,:,:,i))';
        mat_xi__alpha_Mu_(i,:) = xi_alpha_T * neg_mu_i_T;
        mat_xi__z_Mu_(i,:) = xi__z_T * neg_mu_i_T;
        %
        xiT_Omega_i = squeeze(xiT_Omega(:,:,:,i));
        mat_xi_xi__xiTOmega(i,:) = xi_xi__T * xiT_Omega_i;
        mat_xi__alpha_xiTOmega(i,:) = xi_alpha_T * xiT_Omega_i;
        %
        mat_xi_xi__alphaTOmega(i,:) = xi_xi__T * squeeze(alphaT_Omega(:,:,:,i));
    end
    
    % M
    z_alphaT = z * alpha';
    z_alphaT_sumT = z_alphaT + z_alphaT';
    d2_zT_kHess_alpha_dxi_2 = LLT .* z_alphaT_sumT;
    xi_alphaT = xi * alpha';
    xi_alphaT_sumT = xi_alphaT + xi_alphaT';
    d2_xiT_d2kdchidxi_alpha_dxi_2 = Xi_Mu .* xi_alphaT_sumT;
    xi_xiT = xi * xi';
    d2_nablachikT_alpha_dxi_2 = xi_xiT .* Alpha_Mu;
    % C
    d2_zT_kHess_alpha_dxi_dchi_ = diag(z) * mat_xi__alpha_Mu_ + diag(alpha) * mat_xi__z_Mu_;
    d2_xiT_d2kdchidxi_alpha_dxi_dchi_ = diag(alpha) * mat_xi_xi__xiTOmega + diag(xi) * mat_xi__alpha_xiTOmega;
    d2_nablachikT_alpha_dxi_dchi_ = diag(xi) * mat_xi_xi__alphaTOmega;
    sumk_xik_dchik_d2_zT_kHess_alpha_dxi_2 = Xi__Mu_ .* z_alphaT_sumT;
    sumk_xik_dchik_d2_xiT_d2kdchidxi_alpha_dxi_2 = xiT_Omega_xi_ .* xi_alphaT_sumT;
    sumk_xik_dchik_d2_nablachikT_alpha_dxi_2 = xi_xiT .* alphaT_Omega_xi_;
    % g
    index_vec = 1+(d-1)*N:d*N;
    chi_min_chi_T_Pginv2 = chi_min_chi_' * P_g_inv2;
    nabla_chi__k_ygT(:,index_vec) = rho_g_sq * exp(-1/2*chi_min_chi_T_Pginv2*chi_min_chi_) * P_g_inv2 * (-chi_min_chi_*chi_min_chi_T_Pginv2 + eye(N));
    chi_min_eT_Pginv2 = chi_min_e' * P_g_inv2;
    nabla_chi__k_ygT_e(:,index_vec) = rho_g_sq * exp(-1/2*chi_min_eT_Pginv2*chi_min_e) * P_g_inv2 * (-chi_min_e*chi_min_eT_Pginv2 + eye(N));
    chi_min_chi_T_Pginv2 = chi_min_chi_' * P_g_inv2;
    k_g_chichi_times_P_g_inv2 = rho_g_sq * exp(-1/2*chi_min_chi_T_Pginv2*chi_min_chi_) * P_g_inv2;
    k_yg(index_vec) = k_g_chichi_times_P_g_inv2 * chi__min_chi;
    chi_min_eT_Pginv2 = chi_min_e' * P_g_inv2;
    k_g_chietimes_P_g_inv2 = rho_g_sq * exp(-1/2*chi_min_eT_Pginv2*chi_min_e) * P_g_inv2;
    k_yg_e(index_vec) = k_g_chietimes_P_g_inv2 * (-chi_min_e);

    % M
    d2_cfTalpha_dxi_2 = d2_zT_kHess_alpha_dxi_2 + d2_xiT_d2kdchidxi_alpha_dxi_2 - d2_nablachikT_alpha_dxi_2;
    % C
    d2_cfTalpha_dxi_dchi_ = d2_zT_kHess_alpha_dxi_dchi_ + d2_xiT_d2kdchidxi_alpha_dxi_dchi_ - d2_nablachikT_alpha_dxi_dchi_;
    sumk_xik_dchik_d2_cfTalpha_dxi_2 = sumk_xik_dchik_d2_zT_kHess_alpha_dxi_2 + sumk_xik_dchik_d2_xiT_d2kdchidxi_alpha_dxi_2 - sumk_xik_dchik_d2_nablachikT_alpha_dxi_2;
    
    % Calcs
    M_GP = M_GP + d2_cfTalpha_dxi_2;
    d2_f_GP_dxi_dchi_ = d2_f_GP_dxi_dchi_ + d2_cfTalpha_dxi_dchi_;
    sumk_xik_dchik_d2_f_GP_dxi_2 = sumk_xik_dchik_d2_f_GP_dxi_2 + sumk_xik_dchik_d2_cfTalpha_dxi_2;
    D_GP = D_GP + LLT_dmp .* xi_alphaT_sumT;
    D_GP_de = D_GP_de + LLT_dmp_de .* xi_alphaT_sumT;
end

% M
M_est = hatM(chi_) + 1/2 * M_GP;

% C
C_est = hatC(chi_,xi_) + 1/4 * (d2_f_GP_dxi_dchi_ + sumk_xik_dchik_d2_f_GP_dxi_2 - d2_f_GP_dxi_dchi_');

% g
chi_T_Pginv2 = chi_' * P_g_inv2;
nabla_chi__k_0gT = -rho_g_sq * exp( -1/2 * chi_T_Pginv2 * chi_ ) * P_g_inv2 * [chi_, chi_*chi_T_Pginv2 - eye(N)];
g_est = hat_nabla_chi_g_GP(chi_) + nabla_chi__k_ygT * a - nabla_chi__k_0gT * s;

% g(e)
chi_T_Pginv2 = e' * P_g_inv2;
nabla_chi__k_0gT_e = -rho_g_sq * exp( -1/2 * chi_T_Pginv2 * e ) * P_g_inv2 * [e, e*chi_T_Pginv2 - eye(N)];
g_est_e = hat_nabla_chi_g_GP(e) + nabla_chi__k_ygT_e * a - nabla_chi__k_0gT_e * s;

% G
chi_T_Pginv2 = chi_' * P_g_inv2;
k_g_0chi_ = rho_g_sq * exp( -1/2 * chi_T_Pginv2 * chi_ );
k_0gT = k_g_0chi_ * [1, chi_T_Pginv2];
G_est = hatV(chi_) + k_yg' * a - k_0gT * s;

% G(e)
e_T_Pginv2 = e' * P_g_inv2;
k_g_0e = rho_g_sq * exp( -1/2 * e_T_Pginv2 * e );
k_0gT = k_g_0e * [1, e_T_Pginv2];
G_est_e = hatV(e) + k_yg_e' * a - k_0gT * s;

% D
D_est = hatD(xi_) + 1/2 * D_GP;
%d_est = D_est * xi_;

% d(de)
D_est_de = hatD(de) + 1/2 * D_GP_de;
d_est_de = D_est_de * de;

%%
Sigma_lgp_temp = calc_lgp_var(chi_, xi_, zeros(N,1), X_tr, ddq_meas, LLT_Sigma_sq, L_Sigma_sq_dmp, lambda_isq, lambda_isq_dmp, rho_g_sq, P_g_inv2, D, N, K_y_plus_Pi);

tau_ = M_est * ddq_d(t_) + C_est * dq_d(t_) ...
    + (eye(N)-heaviside_man(e'*g_est)*(e*e')./(1e-3 + norm(e)^2)) * g_est ...
    - g_est_e - K_vfb(K1,K2,K3,Sigma_lgp_temp,e) ...
    - Kd*de - K_vfb(K1,K2,K3,Sigma_lgp_temp,de);
ddq_ = M_est \ (tau_ - C_est * xi_ - g_est);

%%
%tau_gp = M_est * ddq_d(t_) + C_est * dq_d(t_) + g_est + D_est * xi_ - Kp*(chi_-q_d(t_)) - Kd*(xi_-dq_d(t_));
Sigma_lgp = calc_lgp_var(chi_, xi_, ddq_, X_tr, ddq_meas, LLT_Sigma_sq, L_Sigma_sq_dmp, lambda_isq, lambda_isq_dmp, rho_g_sq, P_g_inv2, D, N, K_y_plus_Pi);

% tau_gp = M_est * ddq_d(t_) + C_est * dq_d(t_) ...
%     + (eye(N)-heaviside_man(e'*g_est)*(e*e')./(1e-3 + norm(e)^2)) * g_est ...
%     - g_est_e  - K_vfb(K1,K2,K3,Sigma_lgp,e) - Kp*e ...
%     + (eye(N)-heaviside_man(de'*d_est)*(de*de')./(1e-3 + norm(de)^2)) * d_est ...
%     - d_est_de - K_vfb(K1,K2,K3,Sigma_lgp,de) - Kd*de;

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

function [tau_gp_var_est] = calc_lgp_var(q, dq, ddq, X_tr, z_meas, LLT_Sigma_sq, L_Sigma_sq_dmp, lambda_isq, lambda_isq_dmp, rho_g_sq, P_g_inv2, D, N, K)

% Var
var_mat = Cov_mat_func([q;dq], [q;dq], ddq, ddq, N, lambda_isq, lambda_isq_dmp, rho_g_sq, P_g_inv2, LLT_Sigma_sq, L_Sigma_sq_dmp);

% Cov
K_sum = zeros(D,N,N);
parfor j = 1:D
    K_sum(j,:,:) = Cov_mat_func([q;dq], X_tr(j,:)', ddq, z_meas(j,:)', N, lambda_isq, lambda_isq_dmp, rho_g_sq, P_g_inv2, LLT_Sigma_sq, L_Sigma_sq_dmp);
end
cov_mat = zeros(N,N*D);
for j = 1:D
    cov_mat(:,1+(j-1)*N:j*N) = K_sum(j,:,:);
end

% Full
tau_gp_var_est = var_mat - cov_mat*(K\cov_mat');

end

function [Cov_yy_] = Cov_mat_func(x, x_, z, z_, N, lambda_isq, lambda_isq_dmp, rho_g_sq, P_g_inv2, LLT_Sigma_sq, L_Sigma_sq_dmp)

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
xi_min_xi_ = xi-xi_;
d_mat_dmp = -1/2*lambda_isq_dmp .* [xi_min_xi_(1)^2, 0; xi_min_xi_'*xi_min_xi_, xi_min_xi_(2)^2];
L_dmp = L_Sigma_sq_dmp .* [exp(d_mat_dmp(1,1)), 0; exp(d_mat_dmp(2,1)), exp(d_mat_dmp(2,2))];
LLT_dmp = L_dmp * L_dmp';
mu = zeros(N,1,N,N);
V__Mu_ = zeros(N);
diag_z_xi__Mu__v_ = zeros(N);
V_Mu = zeros(N);
diag_v_Mu_xi_z_ = zeros(N);
Omega = zeros(N,N,N,N);
v_Omega_v_ = zeros(N);
Omega_v_ = zeros(N,1,N,N);
Omega_v = zeros(N,1,N,N);
diag_v_Omegav__xi_xi_ = zeros(N);
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
        Mu_ij = 2 * diag_chi__min_chi * Gamma_km_;
        mu_ij = Mu_ij * zi_sj;
        Omega_ij = Mu_ij * diag(zi_sj) * (-Mu_ij') + 2 * diag( Gamma_km_ * zi_sj );
        % 2
        mu(:,:,i,j) = mu_ij;
        mu(:,:,j,i) = mu_ij;
        V__Mu_(i,j) = -mu_ij' * xi_;
        V__Mu_(j,i) = V__Mu_(i,j);
        % 3
        V_Mu(i,j) = mu_ij' * xi;
        V_Mu(j,i) = V_Mu(i,j);
        % 4
        Omega(:,:,i,j) = Omega_ij;
        Omega(:,:,j,i) = Omega_ij;
        Omega_v_(:,:,i,j) = Omega_ij * xi_;
        Omega_v_(:,:,j,i) = Omega_v_(:,:,i,j);
        v_Omega_v_(i,j) = xi' * Omega_v_(:,:,i,j);
        v_Omega_v_(j,i) = v_Omega_v_(i,j);
        % 6
        Omega_v(:,:,i,j) = Omega_ij * xi;
        Omega_v(:,:,j,i) = Omega_v(:,:,i,j);
        % 9
        dchi_dchi__kf = dchi_dchi__kf + 2 * xi(i) * xi(j) * xi_(i) * xi_(j) * Omega_ij;
    end
    %
    Gamma_km_ = lambda_isq * [ones(N,i),zeros(N,N-i)];
    zi_si = LLT(i,:)' .* LLT(:,i);
    Mu_ii = 2 * diag_chi__min_chi * Gamma_km_;
    mu_ii = Mu_ii * zi_si;
    Omega_ii = Mu_ii * diag(zi_si) * (-Mu_ii') + 2 * diag(Gamma_km_ * zi_si);
    % 2
    mu(:,:,i,i) = mu_ii;
    V__Mu_(i,i) = -mu_ii' * xi_;
    %
    %mu_d(:,:,i,i) = mu_ii_d;
    % 3
    V_Mu(i,i) = mu_ii' * xi;
    % 4
    Omega(:,:,i,i) = Omega_ii;
    Omega_v_(:,:,i,i) = Omega_ii * xi_;
    v_Omega_v_(i,i) = xi' * Omega_v_(:,:,i,i);
    % 6
    Omega_v(:,:,i,i) = Omega_ii * xi;
    % 9
    dchi_dchi__kf = dchi_dchi__kf + xi(i) * xi(i) * xi_(i) * xi_(i) * Omega_ii;
end
z_xi__T = (z .* xi_)';
xi_z_ = xi .* z_;
xi_xi_ = xi .* xi_;
for i = 1:N
    %
    mu_i = squeeze(mu(:,:,:,i));
    % 5
    mat_z_xi__Mu_(i,:) = z_xi__T * (-mu_i');
    % 2
    diag_z_xi__Mu__v_(i,i) = mat_z_xi__Mu_(i,:) * xi_;
    % 7
    mat_Mu_xi_z_(:,i) = mu_i * xi_z_;
    % 3
    diag_v_Mu_xi_z_(i,i) = xi' * mat_Mu_xi_z_(:,i);
    % 8
    mat_Omegav__xi_xi_(:,i) = squeeze(Omega_v_(:,:,:,i)) * xi_xi_;
    % 4
    diag_v_Omegav__xi_xi_(i,i) = xi' * mat_Omegav__xi_xi_(:,i);
    % 6
    mat_xi_xi__Omegav(i,:) = xi_xi_' * squeeze(Omega_v(:,:,i,:))';
end

%% kinetic kernel k_T
% 1 - nablas: xi, xi', xi, xi'
C_double = LLT .* (z_*z') + diag(LLT * (z.*z_));

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

%% gravitational kernel k_G
chisT_Pginv2 = chi_min_chi_' * P_g_inv2;
C = C + rho_g_sq * exp( -1/2 * chisT_Pginv2 * chi_min_chi_ ) * P_g_inv2 * (chi__min_chi*chisT_Pginv2 + eye(2));

%% dissipation kernel k_F
C = C + LLT_dmp .* (xi_*xi');

%% Final sum
Cov_yy_ = C;

end