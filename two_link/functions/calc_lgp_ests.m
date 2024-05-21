function [tau_est,ddq_est,m_ij_est,d2_f_dxi_dchi__est,df_dchi__est,dg_dchi__est,f_est,g_est,d_ij_est,P_dsp_est] = ...
    calc_lgp_ests(X, z, LLT_Sigma_sq, lambda_isq, rho_g_sq, P_g_inv2, L_Sigma_sq_dmp, lambda_isq_dmp, X_, ddq_, tau_, pb, A_inv_times_y_min_mu, K_times_a, D, N)
%UNTITLED15 Summary of this function goes here
%   Detailed explanation goes here

%% Initializations
T = size(X_,1);
tau_est = zeros(T,N);
ddq_est = zeros(T,N);
m_ij_est = zeros(T,N*(N+1)/2);
d2_f_dxi_dchi__est = zeros(T,N^2);
df_dchi__est = zeros(T,N);
dg_dchi__est = zeros(T,N);
f_est = zeros(T,1);
g_est = zeros(T,1);
P_dsp_est = zeros(T,1);
d_ij_est = zeros(T,N*(N+1)/2);

%% compute estimates at query points
for t = 1:T

    % compute GP estimates
    [tau_est_t,ddq_est_t,M_est,d2_f_est,nabla_chi_f_est,nabla_chi_g_est,f_kin,g_pot,D_est,dsp_est] = ...
        eval_lgp(...
        X_(t,1:N)', X_(t,N+1:end)', tau_(t,:)', ddq_(t,:)', X, z, ...
        LLT_Sigma_sq, lambda_isq, rho_g_sq, P_g_inv2, L_Sigma_sq_dmp, lambda_isq_dmp, ...
        A_inv_times_y_min_mu, K_times_a, pb, D, N ...
        );

    % f, g & F
    f_est(t) = f_kin;
    g_est(t) = g_pot;
    P_dsp_est(t) = dsp_est;
    
    % nabla_chi f
    df_dchi__est(t,:) = nabla_chi_f_est;
    
    % nabla_chi g
    dg_dchi__est(t,:) = nabla_chi_g_est;
    
    % nabla_xi^2 f & D
    for j = 1:N
        m_ij_est(t,j) = M_est(1,j);
        d_ij_est(t,j) = D_est(1,j);
    end
    for i = 2:N
        for j = i:N
            m_ij_est(t,(2-i+2*N)*(i-1)/2+j-i+1) = M_est(i,j);
            d_ij_est(t,(2-i+2*N)*(i-1)/2+j-i+1) = D_est(i,j);
        end
    end
    
    % nabla_xi nabla_chi^T f
    for i = 1:N
        for j = 1:N
            d2_f_dxi_dchi__est(t,j+(i-1)*N) = d2_f_est(i,j);
        end
    end
    
    % tau
    tau_est(t,:) = tau_est_t;

    % ddq
    ddq_est(t,:) = ddq_est_t;
end

end