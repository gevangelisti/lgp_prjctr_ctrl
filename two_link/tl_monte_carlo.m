%% Open-source implementation of Lagrangian-Gaussian Processes presented in:
% [1] Evangelisti et al., (2024). Exponentially Stable Projector-based Control of Lagrangian Systems with Gaussian Processes, ArXiv
% [2] Evangelisti and Hirche, (2022). Physically Consistent Learning of Conservative Lagrangian Systems with Gaussian Processes, Conference on Decision and Control (CDC), https://doi.org/10.1109/CDC51059.2022.9993123 
% [3] Evangelisti and Hirche, (2024). Data-Driven Momentum Observers With Physically Consistent Gaussian Processes, Transactions on Robotics (T-RO), https://doi.org/10.1109/TRO.2024.3366818 
% 
% Note: Run tl_main first!
%
% Citation: If you use this code within your work/paper, please cite:
%{
@unpublished{evangeli2024control,
	author = {Evangelisti, Giulio and Della Santina, Cosimo and Hirche, Sandra},
	title = {Exponentially Stable Projector-based Control of Lagrangian Systems with Gaussian Processes},
	note = {arXiv},
	year = {2024}}
@inproceedings{evangeli2022learning,
	author = {Evangelisti, Giulio and Hirche, Sandra},
	booktitle = {2022 IEEE 61st Conference on Decision and Control (CDC)},
	pages = {4078-4085},
	title = {Physically Consistent Learning of Conservative Lagrangian Systems with Gaussian Processes},
	year = {2022}}
@article{evangeli2024momentum,
	author = {Evangelisti, Giulio and Hirche, Sandra},
	journal = {IEEE Transactions on Robotics},
	pages = {1-14},
	title = {Data-Driven Momentum Observers With Physically Consistent Gaussian Processes},
	year = {2024}}
%}

%%
rng(0);
close all
do_sim = 0;

%%
N_mc_runs = 1e2;
omega_freq_vec = 0.1:0.1:10;
N_freqs = length(omega_freq_vec);

mean_e_pdp_perf = zeros(N_mc_runs,N_freqs); mean_de_pdp_perf = zeros(N_mc_runs,N_freqs); e_pdp_L2_perf = zeros(N_mc_runs,N_freqs);
tau_pdp_max_perf = zeros(N_mc_runs,N_freqs); tau_pdp_mean_perf = zeros(N_mc_runs,N_freqs); tau_pdp_L2_perf = zeros(N_mc_runs,N_freqs);

mean_e_nat_pdp_perf = zeros(N_mc_runs,N_freqs); mean_de_nat_pdp_perf = zeros(N_mc_runs,N_freqs); e_nat_pdp_L2_perf = zeros(N_mc_runs,N_freqs);
tau_nat_pdp_max_perf = zeros(N_mc_runs,N_freqs); tau_nat_pdp_mean_perf = zeros(N_mc_runs,N_freqs); tau_nat_pdp_L2_perf = zeros(N_mc_runs,N_freqs);

mean_e_lgp_pdp_perf = zeros(N_mc_runs,N_freqs); mean_de_lgp_pdp_perf = zeros(N_mc_runs,N_freqs); e_lgp_pdp_L2_perf = zeros(N_mc_runs,N_freqs);
tau_lgp_pdp_max_perf = zeros(N_mc_runs,N_freqs); tau_lgp_pdp_mean_perf = zeros(N_mc_runs,N_freqs); tau_lgp_pdp_L2_perf = zeros(N_mc_runs,N_freqs);

mean_e_lgp_nat_pdp_perf = zeros(N_mc_runs,N_freqs); mean_de_lgp_nat_pdp_perf = zeros(N_mc_runs,N_freqs); e_lgp_nat_pdp_L2_perf = zeros(N_mc_runs,N_freqs);
tau_lgp_nat_pdp_max_perf = zeros(N_mc_runs,N_freqs); tau_lgp_nat_pdp_mean_perf = zeros(N_mc_runs,N_freqs); tau_lgp_nat_pdp_L2_perf = zeros(N_mc_runs,N_freqs);

mean_e_lgp_nat_pdp_v_perf = zeros(N_mc_runs,N_freqs); mean_de_lgp_nat_pdp_v_perf = zeros(N_mc_runs,N_freqs); e_lgp_nat_pdp_v_L2_perf = zeros(N_mc_runs,N_freqs);
tau_lgp_nat_pdp_v_max_perf = zeros(N_mc_runs,N_freqs); tau_lgp_nat_pdp_v_mean_perf = zeros(N_mc_runs,N_freqs); tau_lgp_nat_pdp_v_L2_perf = zeros(N_mc_runs,N_freqs);

if do_sim == 0
    load('monte_carlo_tl.mat')
else
    for id_omega = N_freqs:-1:1
    
        q_d = @(t) pi/2*sin(omega_freq_vec(id_omega)*t)*ones(N,1);
        dq_d = @(t) pi/2*omega_freq_vec(id_omega)*cos(omega_freq_vec(id_omega)*t)*ones(N,1);
        ddq_d = @(t) -omega_freq_vec(id_omega)^2*pi/2*sin(omega_freq_vec(id_omega)*t)*ones(N,1);
        if omega_freq_vec(id_omega) ~= 0
            tend = 4*2*pi/omega_freq_vec(id_omega);
        else
            tend = 4*2*pi/omega_freq_vec(id_omega+1);
        end
    
        parfor n_mc = 1:N_mc_runs
            
            x0 = pi/4*(2*rand(2*N,1)-1);
            
            % Simulation with classical PD+
            tau_pdp = @(t,q,dq) hatM(q)*ddq_d(t) + hatC(q,dq)*dq_d(t) + hatG(q) + hatD(dq)*dq - Kp*(q-q_d(t)) - Kd*(dq-dq_d(t));
            ode_fun = @(t,x) [x(3); x(4); M(x(1:2)) \ (tau_pdp(t,x(1:2),x(3:4)) - C(x(1:2),x(3:4))*x(3:4) - G(x(1:2)) - Dmp(x(3:4))*x(3:4))];
            [t_pdp,x_pdp] = ode45(ode_fun,linspace(0,tend,F_s*tend+1),x0);
    
            % Simulation with parametric nat-PD+
            tau_sc_param = @(t,q,dq) hatM(q)*ddq_d(t) + hatC(q,dq)*dq_d(t) ...
                + (eye(N)-heaviside((q-q_d(t))'*hatG(q))*(q-q_d(t))*(q-q_d(t))'./(1e-3 + norm((q-q_d(t)))^2) )*hatG(q) - hatG(q-q_d(t)) ...
                + (eye(N)-heaviside((dq-dq_d(t))'*hatD(dq)*dq)*(dq-dq_d(t))*(dq-dq_d(t))'./(1e-3 + norm((dq-dq_d(t)))^2) )*hatD(dq)*dq - hatD(dq-dq_d(t))*(dq-dq_d(t)) - Kd*(dq-dq_d(t));
            ode_fun_sc = @(t,x) [x(3); x(4); M(x(1:2)) \ (tau_sc_param(t,x(1:2),x(3:4)) - C(x(1:2),x(3:4))*x(3:4) - G(x(1:2)) - Dmp(x(3:4))*x(3:4))];
            [t_nat_pdp,x_scp] = ode45(ode_fun_sc,linspace(0,tend,F_s*tend+1),x0);
            
            % L-GP-PD+
            tau_lgp_pdp_const = @(t,q,dq) tau_lgp_pdp_mex(t, q, dq, pi/2, omega_freq_vec(id_omega), Kp_s, Kd_s, X_tr, ddq_meas, ...
                LLT_Sigma_sq, lambda_isq, rho_g_sq, P_g_inv2, L_Sigma_sq_dmp, lambda_isq_dmp, A_inv_times_y_min_mu, K_times_a, est_params, D, N);
            ode_robot = @(t,x) [x(3); x(4); M(x(1:2)) \ (tau_lgp_pdp_const(t,x(1:2),x(3:4)) - C(x(1:2),x(3:4))*x(3:4) - G(x(1:2)) - Dmp(x(3:4))*x(3:4))];
            [t_lgp_pdp,x_lgp_pdp] = ode45(ode_robot,linspace(0,tend,F_s*tend+1),x0);
    
            % Proposed structure-preserving tracking control
            tau_lgp_nat_pdp_const = @(t,q,dq) tau_lgp_nat_pdp_mex(t, q, dq, pi/2, omega_freq_vec(id_omega), 1*Kp_s, 1*Kd_s, X_tr, ddq_meas, ...
                LLT_Sigma_sq, lambda_isq, rho_g_sq, P_g_inv2, L_Sigma_sq_dmp, lambda_isq_dmp, A_inv_times_y_min_mu, K_times_a, est_params, D, N);
            ode_robot_nat = @(t,x) [x(3); x(4); M(x(1:2)) \ (tau_lgp_nat_pdp_const(t,x(1:2),x(3:4)) - C(x(1:2),x(3:4))*x(3:4) - G(x(1:2)) - Dmp(x(3:4))*x(3:4))];
            [t_lgp_nat_pdp,x_lgp_nat_pdp] = ode45(ode_robot_nat,linspace(0,tend,F_s*tend+1),x0);
            
            % Proposed uncertainty-adaptive variant
            tau_lgp_nat_pdp_var = @(t,q,dq) tau_lgp_var_nat_pdp_mex(t, q, dq, pi/2, omega_freq_vec(id_omega), ...
                1e2*eye(2), 2e-2*eye(2), 1/sqrt(2e-2*(1/1e0-1/1e2))*eye(2), 1*Kp_s, 1*Kd_s, A, ... %K_y_plus_Pi, ...
                X_tr, ddq_meas, LLT_Sigma_sq, lambda_isq, rho_g_sq, P_g_inv2, L_Sigma_sq_dmp, lambda_isq_dmp, A_inv_times_y_min_mu, K_times_a, est_params, D, N);
            ode_robot_nat_var = @(t,x) [x(3); x(4); M(x(1:2)) \ (tau_lgp_nat_pdp_var(t,x(1:2),x(3:4)) - C(x(1:2),x(3:4))*x(3:4) - G(x(1:2)) - Dmp(x(3:4))*x(3:4))];
            [t_lgp_var_nat_pdp,x_lgp_var_nat_pdp] = ode45(ode_robot_nat_var,linspace(0,tend,F_s*tend+1),x0);
            
            % Compute errors
            e1_pdp = x_pdp(:,1)-pi/2*sin(omega_freq_vec(id_omega)*t_pdp);
            e2_pdp = x_pdp(:,2)-pi/2*sin(omega_freq_vec(id_omega)*t_pdp);
            e1_nat_pdp = x_scp(:,1)-pi/2*sin(omega_freq_vec(id_omega)*t_nat_pdp);
            e2_nat_pdp = x_scp(:,2)-pi/2*sin(omega_freq_vec(id_omega)*t_nat_pdp);
            e1_lgp_pdp = x_lgp_pdp(:,1)-pi/2*sin(omega_freq_vec(id_omega)*t_lgp_pdp);
            e2_lgp_pdp = x_lgp_pdp(:,2)-pi/2*sin(omega_freq_vec(id_omega)*t_lgp_pdp);
            e1_lgp_nat_pdp = x_lgp_nat_pdp(:,1)-pi/2*sin(omega_freq_vec(id_omega)*t_lgp_nat_pdp);
            e2_lgp_nat_pdp = x_lgp_nat_pdp(:,2)-pi/2*sin(omega_freq_vec(id_omega)*t_lgp_nat_pdp);
            
            e1_lgp_nat_pdp_v = x_lgp_var_nat_pdp(:,1)-pi/2*sin(omega_freq_vec(id_omega)*t_lgp_var_nat_pdp);
            e2_lgp_nat_pdp_v = x_lgp_var_nat_pdp(:,2)-pi/2*sin(omega_freq_vec(id_omega)*t_lgp_var_nat_pdp);
            de1_pdp = x_pdp(:,3)-pi/2*omega_freq_vec(id_omega)*cos(omega_freq_vec(id_omega)*t_pdp);
            de2_pdp = x_pdp(:,4)-pi/2*omega_freq_vec(id_omega)*cos(omega_freq_vec(id_omega)*t_pdp);
            de1_nat_pdp = x_scp(:,3)-pi/2*omega_freq_vec(id_omega)*cos(omega_freq_vec(id_omega)*t_nat_pdp);
            de2_nat_pdp = x_scp(:,4)-pi/2*omega_freq_vec(id_omega)*cos(omega_freq_vec(id_omega)*t_nat_pdp);
            de1_lgp_pdp = x_lgp_pdp(:,3)-pi/2*omega_freq_vec(id_omega)*cos(omega_freq_vec(id_omega)*t_lgp_pdp);
            de2_lgp_pdp = x_lgp_pdp(:,4)-pi/2*omega_freq_vec(id_omega)*cos(omega_freq_vec(id_omega)*t_lgp_pdp);
            de1_lgp_nat_pdp = x_lgp_nat_pdp(:,3)-pi/2*omega_freq_vec(id_omega)*cos(omega_freq_vec(id_omega)*t_lgp_nat_pdp);
            de2_lgp_nat_pdp = x_lgp_nat_pdp(:,4)-pi/2*omega_freq_vec(id_omega)*cos(omega_freq_vec(id_omega)*t_lgp_nat_pdp);
            de1_lgp_nat_pdp_v = x_lgp_var_nat_pdp(:,3)-pi/2*omega_freq_vec(id_omega)*cos(omega_freq_vec(id_omega)*t_lgp_var_nat_pdp);
            de2_lgp_nat_pdp_v = x_lgp_var_nat_pdp(:,4)-pi/2*omega_freq_vec(id_omega)*cos(omega_freq_vec(id_omega)*t_lgp_var_nat_pdp);
            
            % Compute controller norms 
            tau_pdp_norm_vec = zeros(size(t_pdp));
            for i = 1:length(t_pdp)
                tau_pdp_norm_vec(i) = norm(tau_pdp(t_pdp(i),x_pdp(i,1:N)',x_pdp(i,N+1:end)'));
            end
            tau_pdp_max = max(tau_pdp_norm_vec);
            tau_pdp_mean = mean(tau_pdp_norm_vec);
                    
            tau_nat_pdp_norm_vec = zeros(size(t_nat_pdp));
            for i = 1:length(t_nat_pdp)
                tau_nat_pdp_norm_vec(i) = norm(tau_sc_param(t_nat_pdp(i),x_scp(i,1:N)',x_scp(i,N+1:end)'));
            end
            tau_nat_pdp_max = max(tau_nat_pdp_norm_vec);
            tau_nat_pdp_mean = mean(tau_nat_pdp_norm_vec);
            
            tau_lgp_pdp_norm_vec = zeros(size(t_lgp_pdp));
            for i = 1:length(t_lgp_pdp)
                tau_lgp_pdp_norm_vec(i) = norm(tau_lgp_pdp_const(t_lgp_pdp(i),x_lgp_pdp(i,1:2)',x_lgp_pdp(i,3:4)'));
            end
            tau_lgp_pdp_max = max(tau_lgp_pdp_norm_vec);
            tau_lgp_pdp_mean = mean(tau_lgp_pdp_norm_vec);
    
            tau_lgp_nat_pdp_norm_vec = zeros(size(t_lgp_nat_pdp));
            for i = 1:length(t_lgp_nat_pdp)
                tau_lgp_nat_pdp_norm_vec(i) = norm(tau_lgp_nat_pdp_const(t_lgp_nat_pdp(i),x_lgp_nat_pdp(i,1:2)',x_lgp_nat_pdp(i,3:4)'));
            end
            tau_lgp_nat_pdp_max = max(tau_lgp_nat_pdp_norm_vec);
            tau_lgp_nat_pdp_mean = mean(tau_lgp_nat_pdp_norm_vec);
            
            tau_lgp_nat_pdp_var_norm_vec = zeros(size(t_lgp_var_nat_pdp));
            for i = 1:length(t_lgp_var_nat_pdp)
                tau_lgp_nat_pdp_var_norm_vec(i) = norm(tau_lgp_nat_pdp_var(t_lgp_var_nat_pdp(i),x_lgp_var_nat_pdp(i,1:2)',x_lgp_var_nat_pdp(i,3:4)'));
            end
            tau_lgp_nat_pdp_var_max = max(tau_lgp_nat_pdp_var_norm_vec);
            tau_lgp_nat_pdp_var_mean = mean(tau_lgp_nat_pdp_var_norm_vec);
            
            % Evaluate
            %T_ss = 0;
            if omega_freq_vec(id_omega) ~= 0
                T_ss = 2*2*pi/omega_freq_vec(id_omega);
            else
                T_ss = 2*2*pi/omega_freq_vec(id_omega+1);
            end
            
            if length(t_pdp) == floor(F_s*tend+1)
                ind = (t_pdp >= T_ss);
                tau_pdp_L2_perf(n_mc,id_omega) = sqrt(trapz(t_pdp(ind), tau_pdp_norm_vec(ind).^2));
                mean_e_pdp_perf(n_mc,id_omega) = mean(sqrt(e1_pdp(ind).^2 + e2_pdp(ind).^2));
                mean_de_pdp_perf(n_mc,id_omega) = mean(sqrt(de1_pdp(ind).^2 + de2_pdp(ind).^2));
                e_pdp_L2_perf(n_mc,id_omega) = sqrt(trapz(t_pdp(ind), e1_pdp(ind).^2+e2_pdp(ind).^2+de1_pdp(ind).^2+de2_pdp(ind).^2));
                tau_pdp_max_perf(n_mc,id_omega) = max(tau_pdp_norm_vec(ind));
                tau_pdp_mean_perf(n_mc,id_omega) = mean(tau_pdp_norm_vec(ind));
            end
            
            if length(t_nat_pdp) == floor(F_s*tend+1)
                ind = (t_nat_pdp >= T_ss);
                tau_nat_pdp_L2_perf(n_mc,id_omega) = sqrt(trapz(t_nat_pdp(ind), tau_nat_pdp_norm_vec(ind).^2));
                mean_e_nat_pdp_perf(n_mc,id_omega) = mean(sqrt(e1_nat_pdp(ind).^2 + e2_nat_pdp(ind).^2));
                mean_de_nat_pdp_perf(n_mc,id_omega) = mean(sqrt(de1_nat_pdp(ind).^2 + de2_nat_pdp(ind).^2));
                e_nat_pdp_L2_perf(n_mc,id_omega) = sqrt(trapz(t_nat_pdp(ind), e1_nat_pdp(ind).^2+e2_nat_pdp(ind).^2+de1_nat_pdp(ind).^2+de2_nat_pdp(ind).^2));
                tau_nat_pdp_max_perf(n_mc,id_omega) = max(tau_nat_pdp_norm_vec(ind));
                tau_nat_pdp_mean_perf(n_mc,id_omega) = mean(tau_nat_pdp_norm_vec(ind));
            end
            
            if length(t_lgp_pdp) == floor(F_s*tend+1)
                ind = (t_lgp_pdp >= T_ss);
                tau_lgp_pdp_L2_perf(n_mc,id_omega) = sqrt(trapz(t_lgp_pdp(ind), tau_lgp_pdp_norm_vec(ind).^2));
                mean_e_lgp_pdp_perf(n_mc,id_omega) = mean(sqrt(e1_lgp_pdp(ind).^2 + e2_lgp_pdp(ind).^2));
                mean_de_lgp_pdp_perf(n_mc,id_omega) = mean(sqrt(de1_lgp_pdp(ind).^2 + de2_lgp_pdp(ind).^2));
                e_lgp_pdp_L2_perf(n_mc,id_omega) = sqrt(trapz(t_lgp_pdp(ind), e1_lgp_pdp(ind).^2+e2_lgp_pdp(ind).^2+de1_lgp_pdp(ind).^2+de2_lgp_pdp(ind).^2));
                tau_lgp_pdp_max_perf(n_mc,id_omega) = max(tau_lgp_pdp_norm_vec(ind));
                tau_lgp_pdp_mean_perf(n_mc,id_omega) = mean(tau_lgp_pdp_norm_vec(ind));
            end
            
            if length(t_lgp_nat_pdp) == floor(F_s*tend+1)
                ind = (t_lgp_nat_pdp >= T_ss);
                tau_lgp_nat_pdp_L2_perf(n_mc,id_omega) = sqrt(trapz(t_lgp_nat_pdp(ind), tau_lgp_nat_pdp_norm_vec(ind).^2));
                mean_e_lgp_nat_pdp_perf(n_mc,id_omega) = mean(sqrt(e1_lgp_nat_pdp(ind).^2 + e2_lgp_nat_pdp(ind).^2));
                mean_de_lgp_nat_pdp_perf(n_mc,id_omega) = mean(sqrt(de1_lgp_nat_pdp(ind).^2 + de2_lgp_nat_pdp(ind).^2));
                e_lgp_nat_pdp_L2_perf(n_mc,id_omega) = sqrt(trapz(t_lgp_nat_pdp(ind), e1_lgp_nat_pdp(ind).^2+e2_lgp_nat_pdp(ind).^2+de1_lgp_nat_pdp(ind).^2+de2_lgp_nat_pdp(ind).^2));
                tau_lgp_nat_pdp_max_perf(n_mc,id_omega) = max(tau_lgp_nat_pdp_norm_vec(ind));
                tau_lgp_nat_pdp_mean_perf(n_mc,id_omega) = mean(tau_lgp_nat_pdp_norm_vec(ind));
            end
            
            if length(t_lgp_var_nat_pdp) == floor(F_s*tend+1)
                ind = (t_lgp_var_nat_pdp >= T_ss);
                tau_lgp_nat_pdp_v_L2_perf(n_mc,id_omega) = sqrt(trapz(t_lgp_var_nat_pdp(ind), tau_lgp_nat_pdp_var_norm_vec(ind).^2));
                mean_e_lgp_nat_pdp_v_perf(n_mc,id_omega) = mean(sqrt(e1_lgp_nat_pdp_v(ind).^2 + e2_lgp_nat_pdp_v(ind).^2));
                mean_de_lgp_nat_pdp_v_perf(n_mc,id_omega) = mean(sqrt(de1_lgp_nat_pdp_v(ind).^2 + de2_lgp_nat_pdp_v(ind).^2));
                e_lgp_nat_pdp_v_L2_perf(n_mc,id_omega) = sqrt(trapz(t_lgp_var_nat_pdp(ind), e1_lgp_nat_pdp_v(ind).^2+e2_lgp_nat_pdp_v(ind).^2+de1_lgp_nat_pdp_v(ind).^2+de2_lgp_nat_pdp_v(ind).^2));
                tau_lgp_nat_pdp_v_max_perf(n_mc,id_omega) = max(tau_lgp_nat_pdp_var_norm_vec(ind));
                tau_lgp_nat_pdp_v_mean_perf(n_mc,id_omega) = mean(tau_lgp_nat_pdp_var_norm_vec(ind));
            end
        end
    end
end

%%

if omega_freq_vec(1) == 0
    omega_freq_vec = omega_freq_vec(2:end);
    mean_e_pdp_perf = mean_e_pdp_perf(:,2:end);
    mean_e_lgp_pdp_perf = mean_e_lgp_pdp_perf(:,2:end);
    mean_e_nat_pdp_perf = mean_e_nat_pdp_perf(:,2:end);
    mean_e_lgp_nat_pdp_perf = mean_e_lgp_nat_pdp_perf(:,2:end);
    mean_e_lgp_nat_pdp_v_perf = mean_e_lgp_nat_pdp_v_perf(:,2:end);
    mean_de_pdp_perf = mean_de_pdp_perf(:,2:end);
    mean_de_lgp_pdp_perf = mean_de_lgp_pdp_perf(:,2:end);
    mean_de_nat_pdp_perf = mean_de_nat_pdp_perf(:,2:end);
    mean_de_lgp_nat_pdp_perf = mean_de_lgp_nat_pdp_perf(:,2:end);
    mean_de_lgp_nat_pdp_v_perf = mean_de_lgp_nat_pdp_v_perf(:,2:end);
    tau_pdp_mean_perf = tau_pdp_mean_perf(:,2:end);
    tau_lgp_pdp_mean_perf = tau_lgp_pdp_mean_perf(:,2:end);
    tau_nat_pdp_mean_perf = tau_nat_pdp_mean_perf(:,2:end);
    tau_lgp_nat_pdp_mean_perf = tau_lgp_nat_pdp_mean_perf(:,2:end);
    tau_lgp_nat_pdp_v_mean_perf = tau_lgp_nat_pdp_v_mean_perf(:,2:end);
    e_pdp_L2_perf = e_pdp_L2_perf(:,2:end);
    e_lgp_pdp_L2_perf = e_lgp_pdp_L2_perf(:,2:end);
    e_nat_pdp_L2_perf = e_nat_pdp_L2_perf(:,2:end);
    e_lgp_nat_pdp_L2_perf = e_lgp_nat_pdp_L2_perf(:,2:end);
    e_lgp_nat_pdp_v_L2_perf = e_lgp_nat_pdp_v_L2_perf(:,2:end);
end

ax_fs = 13;

fig1 = figure(1); clf;
subplot(121)
sh_alpha = 0.15;
end_pdp = N_mc_runs - floor(sum(sum(mean_e_pdp_perf==0))/N_mc_runs);
end_lgp_pdp = N_mc_runs - floor(sum(sum(mean_e_lgp_pdp_perf==0))/N_mc_runs);
fill([omega_freq_vec(1:end_pdp), fliplr(omega_freq_vec(1:end_pdp))], [mean(e_pdp_L2_perf(:,1:end_pdp))-std(e_pdp_L2_perf(:,1:end_pdp)), fliplr(mean(e_pdp_L2_perf(:,1:end_pdp))+std(e_pdp_L2_perf(:,1:end_pdp)))],c(3,:),'EdgeColor','none','FaceAlpha',sh_alpha);
hold on
hnd(1) = semilogy(omega_freq_vec(1:end_pdp),mean(e_pdp_L2_perf(:,1:end_pdp)),'LineWidth',mlw,'Color',c(3,:));
mlgp_pdp = mean(e_lgp_pdp_L2_perf(:,1:end_lgp_pdp-2));
slgp_pdp = std(e_lgp_pdp_L2_perf(:,1:end_lgp_pdp-2));
for i = 1:2
    ind_i = (e_lgp_pdp_L2_perf(:,end_lgp_pdp-2+i)~=0);
    mlgp_pdp(end_lgp_pdp-2+i) = mean(e_lgp_pdp_L2_perf(ind_i,end_lgp_pdp-2+i));
    slgp_pdp(end_lgp_pdp-2+i) = std(e_lgp_pdp_L2_perf(ind_i,end_lgp_pdp-2+i));
end
fill([omega_freq_vec(1:end_lgp_pdp), fliplr(omega_freq_vec(1:end_lgp_pdp))], [mlgp_pdp-slgp_pdp, fliplr(mlgp_pdp+slgp_pdp)],c(4,:),'EdgeColor','none','FaceAlpha',sh_alpha);
hnd(2) = semilogy(omega_freq_vec(1:end_lgp_pdp),mlgp_pdp,'LineWidth',mlw,'Color',c(4,:));
min_bnd = mean(e_nat_pdp_L2_perf)-std(e_nat_pdp_L2_perf);
min_bnd(min_bnd<0) = 1e-10;
fill([omega_freq_vec, fliplr(omega_freq_vec)], [min_bnd, fliplr(mean(e_nat_pdp_L2_perf)+std(e_nat_pdp_L2_perf))],c(5,:),'EdgeColor','none','FaceAlpha',sh_alpha);
hnd(3) = semilogy(omega_freq_vec,mean(e_nat_pdp_L2_perf),'LineWidth',mlw,'Color',c(5,:));
fill([omega_freq_vec, fliplr(omega_freq_vec)], [mean(e_lgp_nat_pdp_L2_perf)-std(e_lgp_nat_pdp_L2_perf), fliplr(mean(e_lgp_nat_pdp_L2_perf)+std(e_lgp_nat_pdp_L2_perf))],c(7,:),'EdgeColor','none','FaceAlpha',sh_alpha);
hnd(4) = semilogy(omega_freq_vec,mean(e_lgp_nat_pdp_L2_perf),'LineWidth',mlw,'Color',c(7,:));
fill([omega_freq_vec, fliplr(omega_freq_vec)], [mean(e_lgp_nat_pdp_v_L2_perf)-std(e_lgp_nat_pdp_v_L2_perf), fliplr(mean(e_lgp_nat_pdp_v_L2_perf)+std(e_lgp_nat_pdp_v_L2_perf))],c(1,:),'EdgeColor','none','FaceAlpha',sh_alpha);
hnd(5) = semilogy(omega_freq_vec,mean(e_lgp_nat_pdp_v_L2_perf),'LineWidth',mlw,'Color',c(1,:));
xlabel('$\omega$ (1/s)', 'Interpreter', 'Latex','FontSize',fs)
ylabel('$||[\mbox{\boldmath $e$}^T \, s\dot{\mbox{\boldmath $e$}}^T]||_{\mathcal{L}_2}$ (rad)', 'Interpreter', 'Latex','FontSize',fs)
grid
set(gca,'FontSize',ax_fs,'YScale','log')
ylim([2e-2 2.1e1])
xlim([0.1 10])

subplot(122)
fill([omega_freq_vec(1:end_pdp), fliplr(omega_freq_vec(1:end_pdp))], [mean(tau_pdp_L2_perf(:,1:end_pdp))-std(tau_pdp_L2_perf(:,1:end_pdp)), fliplr(mean(tau_pdp_L2_perf(:,1:end_pdp))+std(tau_pdp_L2_perf(:,1:end_pdp)))],c(3,:),'EdgeColor','none','FaceAlpha',sh_alpha);
hold on
hnd(1) = semilogy(omega_freq_vec(1:end_pdp),mean(tau_pdp_L2_perf(:,1:end_pdp)),'LineWidth',mlw,'Color',c(3,:));
mlgp_pdp = mean(tau_lgp_pdp_L2_perf(:,1:end_lgp_pdp-2));
slgp_pdp = std(tau_lgp_pdp_L2_perf(:,1:end_lgp_pdp-2));
for i = 1:2
    ind_i = (tau_lgp_pdp_L2_perf(:,end_lgp_pdp-2+i)~=0);
    mlgp_pdp(end_lgp_pdp-2+i) = mean(tau_lgp_pdp_L2_perf(ind_i,end_lgp_pdp-2+i));
    slgp_pdp(end_lgp_pdp-2+i) = std(tau_lgp_pdp_L2_perf(ind_i,end_lgp_pdp-2+i));
end
fill([omega_freq_vec(1:end_lgp_pdp), fliplr(omega_freq_vec(1:end_lgp_pdp))], [mlgp_pdp-slgp_pdp, fliplr(mlgp_pdp+slgp_pdp)],c(4,:),'EdgeColor','none','FaceAlpha',sh_alpha);
hnd(2) = semilogy(omega_freq_vec(1:end_lgp_pdp),mlgp_pdp,'LineWidth',mlw,'Color',c(4,:));
min_bnd = mean(tau_nat_pdp_L2_perf)-std(tau_nat_pdp_L2_perf);
min_bnd(min_bnd<0) = 1e-10;
fill([omega_freq_vec, fliplr(omega_freq_vec)], [min_bnd, fliplr(mean(tau_nat_pdp_L2_perf)+std(tau_nat_pdp_L2_perf))],c(5,:),'EdgeColor','none','FaceAlpha',sh_alpha);
hnd(3) = semilogy(omega_freq_vec,mean(tau_nat_pdp_L2_perf),'LineWidth',mlw,'Color',c(5,:));
fill([omega_freq_vec, fliplr(omega_freq_vec)], [mean(tau_lgp_nat_pdp_L2_perf)-std(tau_lgp_nat_pdp_L2_perf), fliplr(mean(tau_lgp_nat_pdp_L2_perf)+std(tau_lgp_nat_pdp_L2_perf))],c(7,:),'EdgeColor','none','FaceAlpha',sh_alpha);
hnd(4) = semilogy(omega_freq_vec,mean(tau_lgp_nat_pdp_L2_perf),'LineWidth',mlw,'Color',c(7,:));
fill([omega_freq_vec, fliplr(omega_freq_vec)], [mean(tau_lgp_nat_pdp_v_L2_perf)-std(tau_lgp_nat_pdp_v_L2_perf), fliplr(mean(tau_lgp_nat_pdp_v_L2_perf)+std(tau_lgp_nat_pdp_v_L2_perf))],c(1,:),'EdgeColor','none','FaceAlpha',sh_alpha);
hnd(5) = semilogy(omega_freq_vec,mean(tau_lgp_nat_pdp_v_L2_perf),'LineWidth',mlw,'Color',c(1,:));
legend(hnd,'SoA (PD+)','new PD+ with $\mbox{\boldmath $\mu_{\tau}$}$ (L-GP-PD+)','ours no L-GP (nat-PD+)','ours $\mbox{\boldmath $\mu_{\tau}$}$ (L-GP nat-PD+)','ours $\mbox{\boldmath $\mu_{\tau}$},\mbox{\boldmath $\Sigma_{\tau}$}$ (L-GP var-nat-PD+)','Interpreter','Latex','Location','Best','FontSize',fs)%,'Orientation','Horizontal')%,'Orientation','Horizontal')
xlabel('$\omega$ (1/s)', 'Interpreter', 'Latex','FontSize',fs)
ylabel('$||\mbox{\boldmath $\tau$}||_{L_2}$ (Nm)', 'Interpreter', 'Latex','FontSize',fs)
grid
set(gca,'FontSize',ax_fs,'YScale','log')
xlim([0.1 4])

axes('position',[.79875 .25 .1 .25])
box on % put box around new pair of axes
iOI = (omega_freq_vec >= 4); % range of t near perturbation
hold on
semilogy(omega_freq_vec(iOI),mean(tau_nat_pdp_L2_perf(:,iOI)),'LineWidth',mlw,'Color',c(5,:))
semilogy(omega_freq_vec(iOI),mean(tau_lgp_nat_pdp_L2_perf(:,iOI)),'LineWidth',mlw,'Color',c(7,:))
semilogy(omega_freq_vec(iOI),mean(tau_lgp_nat_pdp_v_L2_perf(:,iOI)),'LineWidth',mlw,'Color',c(1,:))
set(gca,'FontSize',ax_fs,'YScale','log')
grid
axis tight
