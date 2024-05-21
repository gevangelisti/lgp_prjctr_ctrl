function [K_sum] ...
    = calc_K_yy(X_tr, z_meas, LLT_Sigma_sq, lambda_isq, rho_g_sq, P_g_inv2, L_Sigma_sq_dmp, lambda_isq_dmp, Z, Upsilon, pb, D, N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%%
% X_tr: Dx2N
% z_meas
% psi

%% 
hatM = @(q) [pb.alphah + 2*pb.betah*cos(q(2)), pb.deltah + pb.betah*cos(q(2)); pb.deltah + pb.betah*cos(q(2)), pb.deltah];

%%
K_sum = zeros(N*D);

for i = 1:D
    i_vec = 1+(i-1)*N:i*N;
    for j = 1:i-1
        Cov_yy_ = Cov_mat_func(X_tr(i,:)',X_tr(j,:)',z_meas(i,:)',z_meas(j,:)', N, lambda_isq, lambda_isq_dmp, rho_g_sq, P_g_inv2, LLT_Sigma_sq, L_Sigma_sq_dmp);
        j_vec = 1+(j-1)*N:j*N;
        K_sum(i_vec,j_vec) = Cov_yy_;
        K_sum(j_vec,i_vec) = Cov_yy_';
    end
    Var_y = Cov_mat_func(X_tr(i,:)',X_tr(i,:)',z_meas(i,:)',z_meas(i,:)', N, lambda_isq, lambda_isq_dmp, rho_g_sq, P_g_inv2, LLT_Sigma_sq, L_Sigma_sq_dmp);
    M_0 = hatM(X_tr(i,1:N)');
    Sigma_n = M_0*Z*M_0 + (ones(N)-eye(N)).*Z.*LLT_Sigma_sq + diag(LLT_Sigma_sq*diag(Z)); 
    K_sum(i_vec,i_vec) = Var_y + Sigma_n + Upsilon;
end

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