function [M,D_] = MD_lgp(chi_, xi_, X_tr, ddq_meas, LLT_Sigma_sq, lambda_isq, L_Sigma_sq_dmp, lambda_isq_dmp, a, pb, D, N)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Parametric model estimate
hatM = @(q) ...
    [pb.alphah + 2*pb.betah*cos(q(2)), ...
    pb.deltah + pb.betah*cos(q(2)); ...
    pb.deltah + pb.betah*cos(q(2)), ...
    pb.deltah];
hatD = @(dq) [pb.dh(1)+pb.dh(2)*abs(dq(1)), 0; 0, pb.dh(1)+pb.dh(2)*abs(dq(2))];

%% Lagrange GP computations

% Initializations
M_GP = zeros(N);
D_GP = zeros(N);

% Initializations
mu = zeros(N,1,N,N);
Xi_Mu = zeros(N);
Alpha_Mu = zeros(N);

for d = 1:D
    % Set input/state values
    chi = X_tr(d,1:N)';
    xi = X_tr(d,N+1:2*N)';
    alpha = a(1+(d-1)*N:d*N);
    z = ddq_meas(d,:)';

    % 
    chi_min_chi_ = chi-chi_;
    LLT = exp(-lambda_isq * (chi_min_chi_' * chi_min_chi_)) * LLT_Sigma_sq;
    chi__min_chi = -chi_min_chi_;
    diag_chi__min_chi = diag(chi__min_chi);
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
            %
            mu(:,:,i,j) = mu_ij;
            mu(:,:,j,i) = mu_ij;
            %
            Xi_Mu(i,j) = mu_ij' * xi;
            Xi_Mu(j,i) = Xi_Mu(i,j);
            %
            Alpha_Mu(i,j) = mu_ij' * alpha;
            Alpha_Mu(j,i) = Alpha_Mu(i,j);
        end
        %
        Gamma_km_ = lambda_isq * [ones(N,i),zeros(N,N-i)];
        zi_si = LLT(i,:)' .* LLT(:,i);
        Mu_ii = 2 * diag_chi__min_chi * Gamma_km_;
        mu_ii = Mu_ii * zi_si;
        %
        mu(:,:,i,i) = mu_ii;
        %
        Xi_Mu(i,i) = mu_ii' * xi;
        %
        Alpha_Mu(i,i) = mu_ii' * alpha;
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
    
    % M
    d2_cfTalpha_dxi_2 = d2_zT_kHess_alpha_dxi_2 + d2_xiT_d2kdchidxi_alpha_dxi_2 - d2_nablachikT_alpha_dxi_2;

    % Calc
    M_GP = M_GP + d2_cfTalpha_dxi_2;
    D_GP = D_GP + LLT_dmp .* xi_alphaT_sumT;
end

% M
M = hatM(chi_) + 1/2 * M_GP;

% D
D_ = hatD(xi_) + 1/2 * D_GP;

end