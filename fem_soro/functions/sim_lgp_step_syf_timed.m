function [q_pcc_lgp_sim,dq_pcc_lgp_sim] = sim_lgp_step_syf_timed(t_fem, a_pcc, X_tr, ddq_tr, LLT_Sigma_sq, lambda_isq, rho_g_sq, P_g_inv2, A_inv_times_y_min_mu, K_times_a, est_params, D, Ncc, t_abort)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
tstart = tic;
ode_fun_GP = @(t,x) ode_fun_GP_timed(t, x, tstart, t_abort, a_pcc, X_tr, ddq_tr, LLT_Sigma_sq, lambda_isq, rho_g_sq, P_g_inv2, A_inv_times_y_min_mu, K_times_a, est_params, D, Ncc);

[~,x_pcclgp_sim] = ode15s(ode_fun_GP,t_fem,zeros(2*Ncc,1));
q_pcc_lgp_sim = x_pcclgp_sim(:,1:Ncc)';
dq_pcc_lgp_sim = x_pcclgp_sim(:,Ncc+1:end)';

end

function [dx] = ode_fun_GP_timed(t,x,tstart,t_abort, a_pcc, X_tr, ddq_tr, LLT_Sigma_sq, lambda_isq, rho_g_sq, P_g_inv2, A_inv_times_y_min_mu, K_times_a, est_params, D, Ncc)

diff_t = t_abort - toc(tstart);
if diff_t < 0
    dx = zeros(2*Ncc,1);
else
    dx = [...
        x(Ncc+1:2*Ncc);
        ddq_lgp_pcc_mex(x(1:Ncc), x(Ncc+1:2*Ncc), a_pcc, X_tr, ddq_tr, LLT_Sigma_sq, lambda_isq, rho_g_sq, P_g_inv2, A_inv_times_y_min_mu, K_times_a, est_params, D, Ncc) ...
        ];
end

end