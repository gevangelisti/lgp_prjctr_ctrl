%% Open-source implementation of Lagrangian-Gaussian Processes presented in:
% [1] Evangelisti et al., (2024). Exponentially Stable Projector-based Control of Lagrangian Systems with Gaussian Processes, ArXiv
% [2] Evangelisti and Hirche, (2022). Physically Consistent Learning of Conservative Lagrangian Systems with Gaussian Processes, Conference on Decision and Control (CDC), https://doi.org/10.1109/CDC51059.2022.9993123 
% [3] Evangelisti and Hirche, (2024). Data-Driven Momentum Observers With Physically Consistent Gaussian Processes, Transactions on Robotics (T-RO), https://doi.org/10.1109/TRO.2024.3366818 
% 
% Note: spatial_v2 must be in session path!
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

%% Parametrization
do_sim = 0;
t_end = 6*pi;
F_s = 1e3;
omega_d = 1;
Ncc = 4;
N = 100;

t_cl = linspace(0,t_end,F_s*t_end+1);
T_frame = t_end/omega_d;
a_des = pi/180*[2.5;15;60;135];
q_des = @(t) a_des.*sin(omega_d*t);
da_des_dt = @(t) zeros(Ncc,1);
dq_des = @(t) omega_d*a_des.*cos(omega_d*t) + da_des_dt(t).*sin(omega_d*t);
d2a_des_dt2 = @(t) zeros(Ncc,1);
ddq_des = @(t) -omega_d^2*a_des.*sin(omega_d*t);

x0 = [kron(-Ncc/N*a_des,ones(N/Ncc,1)); zeros(N,1)];
Kp = 1*eye(Ncc);
Kd = 1*eye(Ncc);

C_pcc = [];
for i = 1:Ncc
    C_pcc = blkdiag(C_pcc,ones(1,N/Ncc));
end

%% Simulations
if do_sim == 0
    load('soro_sim_data.mat')
else
    %% State-of-the-Art parametric PD+
    ode_fun_fem_soro_pdp = @(t,x) [...
        x(N+1:end); ...
        FDab_wDK( ...
        soro, x(1:N), x(N+1:end), -k*eye(N)*x(1:N) - d*eye(N)*x(N+1:end) + ...
        kron(tau_pcc_pdp(t,C_pcc*x(1:N),C_pcc*x(N+1:end),q_des,dq_des,ddq_des,Kp,Kd,mu,L,Iz,9.81,k_pcc,d_pcc),ones(N/Ncc,1))) ... 
        ];
    disp('Starting standard PD+ simulation ...')
    [~,x_pdp] = ode15s(ode_fun_fem_soro_pdp,t_cl,x0);
    disp('done')
    showmotion(plnr_so_ro(N),t_cl,x_pdp(:,1:N)')

    %% New L-GP-PD+
    ode_fun_fem_soro_lgp_pdp = @(t,x) [...
        x(N+1:end); ...
        FDab_wDK( ...
        soro, x(1:N), x(N+1:end), -k*eye(N)*x(1:N) - d*eye(N)*x(N+1:end) + ...
        kron(tau_pcc_lgp_pdp_mex(t,C_pcc*x(1:N),C_pcc*x(N+1:end),a_des,omega_d,Kp,Kd,...
        X_tr, ddq_tr', LLT_Sigma_sq, lambda_isq, rho_g_sq, P_g_inv2, A_inv_times_y_min_mu, K_times_a, est_params, D, Ncc),...
        ones(N/Ncc,1))) ... 
        ];
    disp('Starting static L-GP-PD+ simulation ...')
    [~,x_lgp_pdp] = ode15s(ode_fun_fem_soro_lgp_pdp,t_cl,x0);
    disp('done')
    showmotion(plnr_so_ro(N),t_cl,x_lgp_pdp(:,1:N)')

    %% Natural parametric nat-PD+
    ode_fun_fem_soro_nat_pdp = @(t,x) [...
        x(N+1:end); ...
        FDab_wDK( ...
        soro, x(1:N), x(N+1:end), -k*eye(N)*x(1:N) - d*eye(N)*x(N+1:end) + ...
        kron(tau_pcc_nat_pdp(t,C_pcc*x(1:N),C_pcc*x(N+1:end),q_des,dq_des,ddq_des,Kd,Ncc,mu,L,Iz,9.81,k_pcc,d_pcc),ones(N/Ncc,1))) ... 
        ];
    disp('Starting static parametric PD+ simulation ...')
    [~,x_nat_pdp] = ode15s(ode_fun_fem_soro_nat_pdp,t_cl,x0);
    disp('done')
    showmotion(plnr_so_ro(N),t_cl,x_nat_pdp(:,1:N)')

    %% Natural L-GP nat-PD+
    ode_fun_fem_soro_lgp_nat_pdp = @(t,x) [...
        x(N+1:end); ...
        FDab_wDK( ...
        soro, x(1:N), x(N+1:end), -k*eye(N)*x(1:N) - d*eye(N)*x(N+1:end) + ...
        kron(tau_pcc_lgp_nat_pdp_mex(t,C_pcc*x(1:N),C_pcc*x(N+1:end),a_des,omega_d,Kd,...
        X_tr, ddq_tr', LLT_Sigma_sq, lambda_isq, rho_g_sq, P_g_inv2, A_inv_times_y_min_mu, K_times_a, est_params, D, Ncc),...
        ones(N/Ncc,1))) ... 
        ];
    disp('Starting static nat L-GP-PD+ simulation ...')
    [~,x_lgp_nat_pdp] = ode15s(ode_fun_fem_soro_lgp_nat_pdp,t_cl,x0);
    disp('done')
    showmotion(plnr_so_ro(N),t_cl,x_lgp_nat_pdp(:,1:N)')
    
    %% Natural variance-adaptive var-nat L-GP-PD+
    k1 = 1e1;
    k2 = 1e-3;
    k_min = 1e-1;
    k3 = 1/sqrt(k2*(1/1e-1-1/k1));
    ode_fun_fem_soro_lgp_nat_pdp_var = @(t,x) [...
        x(N+1:end); ...
        FDab_wDK( ...
        soro, x(1:N), x(N+1:end), -k*eye(N)*x(1:N) - d*eye(N)*x(N+1:end) + ...
        kron(tau_pcc_lgp_var_nat_pdp_mex(t,C_pcc*x(1:N),C_pcc*x(N+1:end),a_des,omega_d,...
        k1*eye(Ncc),k2*eye(Ncc),k3*eye(Ncc),Kd,K_y_plus_Pi,...
        X_tr, ddq_tr', LLT_Sigma_sq, lambda_isq, rho_g_sq, P_g_inv2, A_inv_times_y_min_mu, K_times_a, est_params, D, Ncc),...
        ones(N/Ncc,1))) ... 
        ];
    disp('Starting variance-adaptive nat-L-GP-PD+ simulation ...')
    [~,x_lgp_nat_pdp_var] = ode15s(ode_fun_fem_soro_lgp_nat_pdp_var,t_cl,x0);
    disp('done')
    showmotion(plnr_so_ro(N),t_cl,x_lgp_nat_pdp_var(:,1:N)')

    %%
    [q_pdp,dq_pdp,~] = get_states_from_fem(N,Ncc,x_pdp(:,1:N)',x_pdp(:,N+1:end)',soro,k,d,0);
    [q_lgp_pdp,dq_lgp_pdp,~] = get_states_from_fem(N,Ncc,x_lgp_pdp(:,1:N)',x_lgp_pdp(:,N+1:end)',soro,k,d,0);
    [q_nat_pdp,dq_nat_pdp,~] = get_states_from_fem(N,Ncc,x_nat_pdp(:,1:N)',x_nat_pdp(:,N+1:end)',soro,k,d,0);
    [q_lgp_nat_pdp,dq_lgp_nat_pdp,~] = get_states_from_fem(N,Ncc,x_lgp_nat_pdp(:,1:N)',x_lgp_nat_pdp(:,N+1:end)',soro,k,d,0);
    [q_lgp_nat_pdp_var,dq_lgp_nat_pdp_var,~] = get_states_from_fem(N,Ncc,x_lgp_nat_pdp_var(:,1:N)',x_lgp_nat_pdp_var(:,N+1:end)',soro,k,d,0);
end

%% Data parsing
q_des_cl = q_des(t_cl);
dq_des_cl = dq_des(t_cl);
e_pdp = q_pdp - q_des(t_cl);
de_pdp = dq_pdp - dq_des(t_cl);
e_lgp_pdp = q_lgp_pdp - q_des(t_cl);
de_lgp_pdp = dq_lgp_pdp - dq_des(t_cl);
e_nat_pdp = q_nat_pdp - q_des(t_cl);
de_nat_pdp = dq_nat_pdp - dq_des(t_cl);
e_lgp_nat_pdp = q_lgp_nat_pdp - q_des(t_cl);
de_lgp_nat_pdp = dq_lgp_nat_pdp - dq_des(t_cl);
e_lgp_nat_pdp_var = q_lgp_nat_pdp_var - q_des(t_cl);
de_lgp_nat_pdp_var = dq_lgp_nat_pdp_var - dq_des(t_cl);

%% Compute errors
disp('Analyzing performance metrics ...')

% Compute L_2-error ||e,de||_L2 trapezoidal numerical integration
e_pdp_L2 = sqrt(trapz(t_cl, sum(e_pdp.^2)+sum(de_pdp.^2)));
e_lgp_pdp_c_L2 = sqrt(trapz(t_cl, sum(e_lgp_pdp.^2)+sum(de_lgp_pdp.^2)));

%% Compute controller norms
tau_pdp_norm_vec = zeros(size(t_cl)); k_pdp_norm_vec = zeros(size(t_cl));
for i = 1:length(t_cl)
    tau_i = tau_pcc_pdp(t_cl(i),q_pdp(:,i),dq_pdp(:,i),q_des,dq_des,ddq_des,Kp,Kd,mu,L,Iz,9.81,k_pcc,d_pcc);
    tau_pdp_norm_vec(i) = norm(tau_i);
    k_pdp_norm_vec(i) = norm(kron(tau_i,ones(N/Ncc,1)) - ID( soro, kron(q_des(t_cl(i))/(N/Ncc),ones(N/Ncc,1)), kron(dq_des(t_cl(i))/(N/Ncc),ones(N/Ncc,1)), kron(ddq_des(t_cl(i))/(N/Ncc),ones(N/Ncc,1)) ) );
end
tau_pdp_max = max(tau_pdp_norm_vec);
tau_pdp_mean = mean(tau_pdp_norm_vec);

%%
tau_lgp_pdp_norm_vec = zeros(size(t_cl)); k_lgp_pdp_norm_vec = zeros(size(t_cl));
M_pdp_mat_t = zeros(Ncc,Ncc,length(t_cl)); M_eig_vec_pdp_t = zeros(Ncc,length(t_cl)); C_pdp_mat_t = zeros(Ncc,Ncc,length(t_cl));
g_pdp_est_q_vec_t = zeros(Ncc,length(t_cl));
for i = 1:length(t_cl)
    [tau_i,M_pdp_mat_t(:,:,i),C_pdp_mat_t(:,:,i),g_pdp_est_q_vec_t(:,i)] = eval_pcc_lgp_pdp_mex(t_cl(i),q_lgp_pdp(:,i),dq_lgp_pdp(:,i),a_des,omega_d,Kp,Kd,...
    X_tr, ddq_tr', LLT_Sigma_sq, lambda_isq, rho_g_sq, P_g_inv2, A_inv_times_y_min_mu, K_times_a, est_params, D, Ncc);
    tau_lgp_pdp_norm_vec(i) = norm(tau_i);
    k_lgp_pdp_norm_vec(i) = norm(kron(tau_i,ones(N/Ncc,1)) - ID( soro, kron(q_des(t_cl(i))/(N/Ncc),ones(N/Ncc,1)), kron(dq_des(t_cl(i))/(N/Ncc),ones(N/Ncc,1)), kron(ddq_des(t_cl(i))/(N/Ncc),ones(N/Ncc,1)) ) );
    M_eig_vec_pdp_t(:,i) = sort(eig(M_pdp_mat_t(:,:,i)));
end
tau_lgp_pdp_max = max(tau_lgp_pdp_norm_vec);
tau_lgp_pdp_mean = mean(tau_lgp_pdp_norm_vec);

%%
tau_nat_pdp_norm_vec = zeros(size(t_cl)); k_nat_pdp_norm_vec = zeros(size(t_cl));
for i = 1:length(t_cl)
    tau_i = tau_pcc_nat_pdp(t_cl(i),q_nat_pdp(:,i),dq_nat_pdp(:,i),q_des,dq_des,ddq_des,Kd,Ncc,mu,L,Iz,9.81,k_pcc,d_pcc);
    tau_nat_pdp_norm_vec(i) = norm(tau_i);
    k_nat_pdp_norm_vec(i) = norm(kron(tau_i,ones(N/Ncc,1)) - ID( soro, kron(q_des(t_cl(i))/(N/Ncc),ones(N/Ncc,1)), kron(dq_des(t_cl(i))/(N/Ncc),ones(N/Ncc,1)), kron(ddq_des(t_cl(i))/(N/Ncc),ones(N/Ncc,1)) ) );
end
tau_nat_pdp_max = max(tau_nat_pdp_norm_vec);
tau_nat_pdp_mean = mean(tau_nat_pdp_norm_vec);

%%
tau_lgp_nat_pdp_norm_vec = zeros(size(t_cl)); k_lgp_nat_pdp_norm_vec = zeros(size(t_cl));
M_mat_t = zeros(Ncc,Ncc,length(t_cl)); M_eig_vec_t = zeros(Ncc,length(t_cl)); C_mat_t = zeros(Ncc,Ncc,length(t_cl));
g_est_e_vec_t = zeros(Ncc,length(t_cl)); g_est_q_vec_t = zeros(Ncc,length(t_cl));
for i = 1:length(t_cl)
    [tau_i,M_mat_t(:,:,i),C_mat_t(:,:,i),D_est_i,g_est_q_vec_t(:,i),g_est_e_vec_t(:,i),D_est_de_i] = eval_pcc_lgp_nat_pdp(t_cl(i),q_lgp_nat_pdp(:,i),dq_lgp_nat_pdp(:,i),a_des,omega_d,Kd,...
    X_tr, ddq_tr', LLT_Sigma_sq, lambda_isq, rho_g_sq, P_g_inv2, A_inv_times_y_min_mu, K_times_a, est_params, D, Ncc);
    tau_lgp_nat_pdp_norm_vec(i) = norm(tau_i);
    k_lgp_nat_pdp_norm_vec(i) = norm(kron(tau_i,ones(N/Ncc,1)) - ID( soro, kron(q_des(t_cl(i))/(N/Ncc),ones(N/Ncc,1)), kron(dq_des(t_cl(i))/(N/Ncc),ones(N/Ncc,1)), kron(ddq_des(t_cl(i))/(N/Ncc),ones(N/Ncc,1)) ) );
    M_eig_vec_t(:,i) = sort(eig(M_mat_t(:,:,i)));
end
tau_lgp_nat_pdp_max = max(tau_lgp_nat_pdp_norm_vec);
tau_lgp_nat_pdp_mean = mean(tau_lgp_nat_pdp_norm_vec);

%%
tau_lgp_nat_pdp_var_norm_vec = zeros(size(t_cl)); Sigma_t = zeros(Ncc,Ncc,length(t_cl)); k_lgp_nat_pdp_var_norm_vec = zeros(size(t_cl));
M_v_mat_t = zeros(Ncc,Ncc,length(t_cl)); M_v_eig_vec_t = zeros(Ncc,length(t_cl)); C_v_mat_t = zeros(Ncc,Ncc,length(t_cl));
de_v_vec_t = zeros(Ncc,length(t_cl)); g_v_est_e_vec_t = zeros(Ncc,length(t_cl)); g_v_est_q_vec_t = zeros(Ncc,length(t_cl));
for i = 1:length(t_cl)
    [tau_i,M_v_mat_t(:,:,i),C_v_mat_t(:,:,i),D_est_i,g_v_est_q_vec_t(:,i),Sigma_i,g_v_est_e_vec_t(:,i),D_est_de_i,e_i,de_i] = eval_pcc_lgp_var_nat_pdp_mex(t_cl(i),q_lgp_nat_pdp_var(:,i),dq_lgp_nat_pdp_var(:,i),a_des,omega_d,...
        k1*eye(Ncc),k2*eye(Ncc),k3*eye(Ncc),Kd,K_y_plus_Pi,...
        X_tr, ddq_tr', LLT_Sigma_sq, lambda_isq, rho_g_sq, P_g_inv2, A_inv_times_y_min_mu, K_times_a, est_params, D, Ncc);
    tau_lgp_nat_pdp_var_norm_vec(i) = norm(tau_i);
    Sigma_t(:,:,i) = Sigma_i;
    k_lgp_nat_pdp_var_norm_vec(i) = norm(kron(tau_i,ones(N/Ncc,1)) - ID( soro, kron(q_des(t_cl(i))/(N/Ncc),ones(N/Ncc,1)), kron(dq_des(t_cl(i))/(N/Ncc),ones(N/Ncc,1)), kron(ddq_des(t_cl(i))/(N/Ncc),ones(N/Ncc,1)) ) );
    M_v_eig_vec_t(:,i) = sort(eig(M_v_mat_t(:,:,i)));
end
tau_lgp_nat_pdp_var_max = max(tau_lgp_nat_pdp_var_norm_vec);
tau_lgp_nat_pdp_var_mean = mean(tau_lgp_nat_pdp_var_norm_vec);

%%
G_est_e_t = zeros(size(t_cl)); G_v_est_e_t = zeros(size(t_cl)); G_est_q_t = zeros(size(t_cl));
for i = 2:length(t_cl)
    G_v_est_e_t(i) = G_v_est_e_t(i-1) + trapz(t_cl(i-1:i), diag(de_lgp_nat_pdp_var(:,i-1:i)'*g_v_est_e_vec_t(:,i-1:i)));
    G_est_e_t(i) = G_est_e_t(i-1) + trapz(t_cl(i-1:i), diag(de_lgp_nat_pdp(:,i-1:i)'*g_est_e_vec_t(:,i-1:i)));
    G_est_q_t(i) = G_est_q_t(i-1) + trapz(t_cl(i-1:i), diag(dq_lgp_nat_pdp_var(:,i-1:i)'*g_v_est_q_vec_t(:,i-1:i)));
end
G_est_e_t = G_est_e_t-min(G_est_e_t) + 1e-6;
G_v_est_e_t = G_v_est_e_t-min(G_v_est_e_t) + 1e-6;
G_est_q_t = G_est_q_t-min(G_est_q_t);

%% Evaluate steady-state metrics
T_ss = 2*pi;
ind = (t_cl >= T_ss);

tau_pdp_ss_L2 = sqrt(trapz(t_cl(ind), tau_pdp_norm_vec(ind).^2));
e_pdp_ss = sqrt(sum(e_pdp(:,ind).^2));
de_pdp_ss = sqrt(sum(de_pdp(:,ind).^2));
e_pdp_ss_L2 = sqrt(trapz(t_cl(ind), sum(e_pdp(:,ind).^2)+sum(de_pdp(:,ind).^2)));
tau_pdp_max_ss = max(tau_pdp_norm_vec(ind));
tau_pdp_mean_ss = mean(tau_pdp_norm_vec(ind));

ind = (t_cl >= T_ss);
tau_lgp_pdp_ss_L2 = sqrt(trapz(t_cl(ind), tau_lgp_pdp_norm_vec(ind).^2));
e_lgp_pdp_ss = sqrt(sum(e_lgp_pdp(:,ind).^2));
de_lgp_pdp_ss = sqrt(sum(de_lgp_pdp(:,ind).^2));
e_lgp_pdp_ss_L2 = sqrt(trapz(t_cl(ind), sum(e_lgp_pdp(:,ind).^2)+sum(de_lgp_pdp(:,ind).^2)));
tau_lgp_pdp_max_ss = max(tau_lgp_pdp_norm_vec(ind));
tau_lgp_pdp_mean_ss = mean(tau_lgp_pdp_norm_vec(ind));

ind = (t_cl >= T_ss);
tau_nat_pdp_ss_L2 = sqrt(trapz(t_cl(ind), tau_nat_pdp_norm_vec(ind).^2));
e_nat_pdp_ss = sqrt(sum(e_nat_pdp(:,ind).^2));
de_nat_pdp_ss = sqrt(sum(de_nat_pdp(:,ind).^2));
e_nat_pdp_ss_L2 = sqrt(trapz(t_cl(ind), sum(e_nat_pdp(:,ind).^2)+sum(de_nat_pdp(:,ind).^2)));
tau_nat_pdp_max_ss = max(tau_nat_pdp_norm_vec(ind));
tau_nat_pdp_mean_ss = mean(tau_nat_pdp_norm_vec(ind));

ind = (t_cl >= T_ss);
tau_lgp_nat_pdp_ss_L2 = sqrt(trapz(t_cl(ind), tau_lgp_nat_pdp_norm_vec(ind).^2));
e_lgp_nat_pdp_ss = sqrt(sum(e_lgp_nat_pdp(:,ind).^2));
de_lgp_nat_pdp_ss = sqrt(sum(de_lgp_nat_pdp(:,ind).^2));
e_lgp_nat_pdp_ss_L2 = sqrt(trapz(t_cl(ind), sum(e_lgp_nat_pdp(:,ind).^2)+sum(de_lgp_nat_pdp(:,ind).^2)));
tau_lgp_nat_pdp_max_ss = max(tau_lgp_nat_pdp_norm_vec(ind));
tau_lgp_nat_pdp_mean_ss = mean(tau_lgp_nat_pdp_norm_vec(ind));

ind = (t_cl >= T_ss);
tau_lgp_nat_pdp_v_ss_L2 = sqrt(trapz(t_cl(ind), tau_lgp_nat_pdp_var_norm_vec(ind).^2));
e_lgp_nat_pdp_v_ss = sqrt(sum(e_lgp_nat_pdp_var(:,ind).^2));
de_lgp_nat_pdp_v_ss = sqrt(sum(de_lgp_nat_pdp_var(:,ind).^2));
e_lgp_nat_pdp_v_ss_L2 = sqrt(trapz(t_cl(ind), sum(e_lgp_nat_pdp_var(:,ind).^2)+sum(de_lgp_nat_pdp_var(:,ind).^2)));
tau_lgp_nat_pdp_var_max_ss = max(tau_lgp_nat_pdp_var_norm_vec(ind));
tau_lgp_nat_pdp_var_mean_ss = mean(tau_lgp_nat_pdp_var_norm_vec(ind));

disp('done')

%%
fprintf('-----------------------------------------------------------------------------------------------------------------------------\n')
fprintf('\t\t\t\t|\tPD+\t\t\tL-GP PD+\tnat-PD+\t\tL-GP nat-PD+\tL-GP var-nat-PD+\n')
fprintf('-----------------------------------------------------------------------------------------------------------------------------\n')
fprintf('||K_p||\t\t\t|\t%.0f\t\t\t%.0f\t\t\t%.0f\t\t\t%.0f\n', norm(Kp), norm(Kp), 0, 0)
fprintf('||K_d||\t\t\t|\t%.0f\t\t\t%.0f\t\t\t%.0f\t\t\t%.0f\n', norm(Kd), norm(Kd), 0, 0)
fprintf('-----------------------------------------------------------------------------------------------------------------------------\n')
fprintf('t >= 2pi:\t\t|\n')
fprintf('||[tau]||_L^2\t|\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t\t%.2f\n', tau_pdp_ss_L2, tau_lgp_pdp_ss_L2, tau_nat_pdp_ss_L2, tau_lgp_nat_pdp_ss_L2, tau_lgp_nat_pdp_v_ss_L2)
fprintf('max(||tau||)\t|\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t\t%.2f\n', tau_pdp_max_ss, tau_lgp_pdp_max_ss, tau_nat_pdp_max_ss, tau_lgp_nat_pdp_max_ss, tau_lgp_nat_pdp_var_max_ss)
fprintf('mean(||tau||)\t|\t%.3f\t\t%.3f\t\t%.3f\t\t%.3f\t\t\t%.3f\n', tau_pdp_mean_ss, tau_lgp_pdp_mean_ss, tau_nat_pdp_mean_ss, tau_lgp_nat_pdp_mean_ss, tau_lgp_nat_pdp_var_mean_ss)
fprintf('||[e, de]||_L^2\t|\t%.3f\t\t%.3f\t\t%.3f\t\t%.3f\t\t\t%.3f\n', e_pdp_ss_L2, e_lgp_pdp_ss_L2, e_nat_pdp_ss_L2, e_lgp_nat_pdp_ss_L2, e_lgp_nat_pdp_v_ss_L2)
fprintf('max(||e(t)||\t|\t%.3f\t\t%.3f\t\t%.3f\t\t%.3f\t\t\t%.3f\n', max(e_pdp_ss), max(e_lgp_pdp_ss), max(e_nat_pdp_ss), max(e_lgp_nat_pdp_ss), max(e_lgp_nat_pdp_v_ss))
fprintf('max(||de||\t\t|\t%.3f\t\t%.3f\t\t%.3f\t\t%.3f\t\t\t%.3f\n', max(de_pdp_ss), max(de_lgp_pdp_ss), max(de_nat_pdp_ss), max(de_lgp_nat_pdp_ss), max(de_lgp_nat_pdp_v_ss))
fprintf('mean(||e(t)||\t|\t%.3f\t\t%.3f\t\t%.3f\t\t%.3f\t\t\t%.3f\n', mean(e_pdp_ss), mean(e_lgp_pdp_ss), mean(e_nat_pdp_ss), mean(e_lgp_nat_pdp_ss), mean(e_lgp_nat_pdp_v_ss))
fprintf('mean(||de||\t\t|\t%.3f\t\t%.3f\t\t%.3f\t\t%.4f\t\t\t%.4f\n', mean(de_pdp_ss), mean(de_lgp_pdp_ss), mean(de_nat_pdp_ss), mean(de_lgp_nat_pdp_ss), mean(de_lgp_nat_pdp_v_ss))
fprintf('-----------------------------------------------------------------------------------------------------------------------------\n')

%%
c = get(gca,'colororder');
close all
fs = 14; % font size
ax_fs = 14; % axis font size
mlw = 1.5; % main line width
msz = 12; % marker size
fig1 = figure(1); clf
tp_end = 7;

% q
subplot(2,3,1)
hold on
plot(t_cl,180/pi*sqrt(sum(q_pdp.^2)),'LineWidth',mlw,'Color',c(3,:))
plot(t_cl,180/pi*sqrt(sum(q_lgp_pdp.^2)),'LineWidth',mlw,'Color',c(4,:))
plot(t_cl,180/pi*sqrt(sum(q_nat_pdp.^2)),'LineWidth',mlw,'Color',c(5,:))
plot(t_cl,180/pi*sqrt(sum(q_lgp_nat_pdp.^2)),'LineWidth',mlw,'Color',c(7,:))
plot(t_cl,180/pi*sqrt(sum(q_lgp_nat_pdp_var.^2)),'LineWidth',mlw,'Color',c(1,:))
plot(t_cl,180/pi*sqrt(sum(q_des_cl.^2)),'--','LineWidth',mlw,'Color',c(2,:))
xlim([0 tp_end])
grid
set(gca,'FontSize',ax_fs)
ylabel('$||\mbox{\boldmath $q$}||$ (deg)','Interpreter','Latex','FontSize',fs)

% dq
subplot(2,3,2)
semilogy(t_cl,180/pi*sqrt(sum(dq_pdp.^2)),'LineWidth',mlw,'Color',c(3,:))
hold on
semilogy(t_cl,180/pi*sqrt(sum(dq_lgp_pdp.^2)),'LineWidth',mlw,'Color',c(4,:))
semilogy(t_cl,180/pi*sqrt(sum(dq_nat_pdp.^2)),'LineWidth',mlw,'Color',c(5,:))
semilogy(t_cl,180/pi*sqrt(sum(dq_lgp_nat_pdp.^2)),'LineWidth',mlw,'Color',c(7,:))
semilogy(t_cl,180/pi*sqrt(sum(dq_lgp_nat_pdp_var.^2)),'LineWidth',mlw,'Color',c(1,:))
semilogy(t_cl,180/pi*sqrt(sum(dq_des_cl.^2)),'--','LineWidth',mlw,'Color',c(2,:))
xlim([0 tp_end]); ylim(180/pi*[min(sqrt(sum(dq_lgp_nat_pdp_var(:,2:end).^2))) max(sqrt(sum(dq_lgp_nat_pdp_var(:,2:end).^2)))])
grid; yticks([1e0 1e2]);
set(gca,'FontSize',ax_fs)
ylabel('$||\dot{\mbox{\boldmath $q$}}||$ (deg/s)','Interpreter','Latex','FontSize',fs)

% e
subplot(2,3,4)
semilogy(t_cl,180/pi*sqrt(sum(e_pdp.^2)),'LineWidth',mlw,'Color',c(3,:))
hold on
semilogy(t_cl,180/pi*sqrt(sum(e_lgp_pdp.^2)),'LineWidth',mlw,'Color',c(4,:))
semilogy(t_cl,180/pi*sqrt(sum(e_nat_pdp.^2)),'LineWidth',mlw,'Color',c(5,:))
semilogy(t_cl,180/pi*sqrt(sum(e_lgp_nat_pdp.^2)),'LineWidth',mlw,'Color',c(7,:))
semilogy(t_cl,180/pi*sqrt(sum(e_lgp_nat_pdp_var.^2)),'LineWidth',mlw,'Color',c(1,:))
xlim([0 tp_end])
grid
set(gca,'FontSize',ax_fs)
ylabel('$||\mbox{\boldmath $e$}||$ (deg)','Interpreter','Latex','FontSize',fs)
xlabel('$t$ (s)','Interpreter','Latex','FontSize',fs)
yticks([1e0 1e2])

% de
subplot(2,3,5)
semilogy(t_cl,180/pi*sqrt(sum(de_pdp.^2)),'LineWidth',mlw,'Color',c(3,:))
hold on
semilogy(t_cl,180/pi*sqrt(sum(de_lgp_pdp.^2)),'LineWidth',mlw,'Color',c(4,:))
semilogy(t_cl,180/pi*sqrt(sum(de_nat_pdp.^2)),'LineWidth',mlw,'Color',c(5,:))
semilogy(t_cl,180/pi*sqrt(sum(de_lgp_nat_pdp.^2)),'LineWidth',mlw,'Color',c(7,:))
semilogy(t_cl,180/pi*sqrt(sum(de_lgp_nat_pdp_var.^2)),'LineWidth',mlw,'Color',c(1,:))
xlim([0 tp_end])
grid
set(gca,'FontSize',ax_fs)
ylabel('$||\dot{\mbox{\boldmath $e$}}||$ (deg/s)','Interpreter','Latex','FontSize',fs)
xlabel('$t$ (s)','Interpreter','Latex','FontSize',fs)
yticks([1e0 1e2])

% tau
tau_ref_norm = zeros(size(t_cl));
for i = 1:length(t_cl)
    tau_ref_norm(i) = norm(ID( soro, kron(q_des(t_cl(i))/(N/Ncc),ones(N/Ncc,1)), kron(dq_des(t_cl(i))/(N/Ncc),ones(N/Ncc,1)), kron(ddq_des(t_cl(i))/(N/Ncc),ones(N/Ncc,1)) ) );
end
subplot(2,3,3)
semilogy(t_cl, tau_pdp_norm_vec, 'LineWidth', mlw,'Color',c(3,:))
hold on
semilogy(t_cl, tau_lgp_pdp_norm_vec, 'LineWidth',mlw,'Color',c(4,:))
semilogy(t_cl, tau_nat_pdp_norm_vec, 'LineWidth',mlw,'Color',c(5,:))
semilogy(t_cl, tau_lgp_nat_pdp_norm_vec, 'LineWidth', mlw,'Color',c(7,:))
semilogy(t_cl, tau_lgp_nat_pdp_var_norm_vec, 'LineWidth', mlw,'Color',c(1,:))
semilogy(t_cl,tau_ref_norm./sqrt(N/Ncc),'--','LineWidth',mlw,'Color',c(2,:))
xlim([0 tp_end])
ylim([1e-1 1e1])
grid
set(gca,'FontSize',ax_fs)
ylabel('$||\mbox{\boldmath $\tau$}||$ (Nm)','Interpreter','Latex','FontSize',fs)
legend('SoA (PD+)','new PD+ with $\mbox{\boldmath $\mu_{\tau}$}$ (L-GP-PD+)','ours no L-GP (nat-PD+)','ours $\mbox{\boldmath $\mu_{\tau}$}$ (L-GP nat-PD+)','ours $\mbox{\boldmath $\mu_{\tau}$},\mbox{\boldmath $\Sigma_{\tau}$}$ (L-GP var-nat-PD+)','Ref','Interpreter','Latex','FontSize',fs,'Location','Best','Orientation','Horizontal')%,'Orientation','Horizontal')

% desired force error
subplot(2,3,6)
semilogy(t_cl,k_pdp_norm_vec,'LineWidth',mlw,'Color',c(3,:))
hold on
semilogy(t_cl,k_lgp_pdp_norm_vec,'LineWidth',mlw,'Color',c(4,:))
semilogy(t_cl,k_nat_pdp_norm_vec,'LineWidth',mlw,'Color',c(5,:))
semilogy(t_cl,k_lgp_nat_pdp_norm_vec,'LineWidth',mlw,'Color',c(7,:))
semilogy(t_cl,k_lgp_nat_pdp_var_norm_vec,'LineWidth',mlw,'Color',c(1,:))
xlim([0 tp_end])
ylim([5e-1 2e1])
grid
set(gca,'FontSize',ax_fs)
ylabel('$||\mbox{\boldmath $\tau$}-\mbox{\boldmath $\tau$}^*||$ (Nm)','Interpreter','Latex','FontSize',fs)
xlabel('$t$ (s)','Interpreter','Latex','FontSize',fs)
yticks([1e0 1e1])

%%
d_lower_est = min(est_params.d);
m_0 = max(M_v_eig_vec_t(Ncc,:)); m_lower = 1e3*min(M_v_eig_vec_t(1,:));
kappa_G_lower_est = 2;

fig3 = figure(3); clf; hold on;
subplot(311)
set(gca,'yscale','log','FontSize',ax_fs)
hold on
subplot(312)
set(gca,'yscale','log','FontSize',ax_fs)
hold on
sh_alpha = 0.15;

ctr = 1; c = colororder; num_sim = 3; hndls = [];

if exist('lgp_error_fac') == 0
    lgp_error_fac = 0.913155517794163;
end

t_mat = {}; x_t_mat = {}; x_bnd_t_mat = {}; V_t_mat = {};

for w_var = [0,1e-1]

    if w_var == 0
        w_nat_vec = 0:1;
    else
        w_nat_vec = 1;
    end

    for w_nat = w_nat_vec

        a0f = @(alph) alph * (Kp(1,1)+w_var);
        a1f = @(alph,the) (Kp(1,1)+w_var) - the/2 - alph * (w_nat*d_lower_est + Kd(1,1)+w_var - alph*(m_0+m_lower));
        delta_m = m_0 - m_lower;
        phif = @(alph,eps,the) 2 * (w_nat*d_lower_est + Kd(1,1)+w_var - a1f(alph,the)*eps + a0f(alph)) - (eps+alph)*(m_0+m_lower);
        af = @(alph,eps,the) -a0f(alph) + a1f(alph,the)*eps;
        upsf = @(alph,eps,the) af(alph,eps,the) - (eps+alph)*delta_m/4 * (1 + sqrt(1 + 4*eps^2*alph^2/(eps+alph)^2));
        
        kappa = @(alph,eps) (Kp(1,1)+w_var) + eps*(Kd(1,1)+w_var+w_nat*d_lower_est) - eps*alph*(m_lower+m_0);
        mu_l = @(alph,eps) (w_nat*kappa_G_lower_est+kappa(alph,eps)+m_lower)/2 - 1/2*sqrt((w_nat*kappa_G_lower_est+kappa(alph,eps)-m_lower)^2+(2*eps*m_lower)^2);
        rho_func = @(alph,eps,the,phi) sqrt((eps/the + 1/phi)/(2*mu_l(alph,eps)));
        
        if w_nat
            u_min = 2/3; 
            if w_var == 0.1
                lb_the = 1;
            else
                lb_the = 0.25;
            end
        else
            u_min = 0.1;
            lb_the = 1;
        end

        rng(0)
        x_test = ga(@(x)rho_func(x(1),x(2),x(3),phif(x(1),x(2),x(3))) + 1/x(1),3,[],[],[],[],[1e-2,0.1,lb_the],[1,10,10],...
            @(x) nl_ineqs_beta_soro(x(1),x(2),x(3),phif(x(1),x(2),x(3)),Kd(1,1)+w_var,w_nat*d_lower_est,Kp(1,1)+w_var,w_nat*kappa_G_lower_est,m_lower,m_0,u_min));
        
        alpha_ct = x_test(1); eps_ct = x_test(2); the_ct  = x_test(3); phi_ct = phif(x_test(1),x_test(2),x_test(3));
        mu_lower = mu_l(alpha_ct,eps_ct);
        rho_ct = rho_func(alpha_ct,eps_ct,the_ct,phi_ct);
        kappa_ct = kappa(alpha_ct,eps_ct);
        
        %
        if w_var == 0
            if w_nat
                % Proposed structure-preserving tracking control
                t_n = t_cl; q_n = q_lgp_nat_pdp; dq_n = dq_lgp_nat_pdp;
                e_n = e_lgp_nat_pdp; de_n = de_lgp_nat_pdp;
                norm_e_sq = sum(e_lgp_nat_pdp.^2) + sum(de_lgp_nat_pdp.^2);
            else
                % PD+
                t_n = t_cl; q_n = q_lgp_pdp; dq_n = dq_nat_pdp;
                e_n = e_lgp_pdp; de_n = de_lgp_pdp;
                norm_e_sq = sum(e_lgp_pdp.^2) + sum(de_lgp_pdp.^2);
            end
        else
            % Proposed uncertainty-adaptive variant
            t_n = t_cl; q_n = q_lgp_nat_pdp_var; dq_n = dq_lgp_nat_pdp_var;
            e_n = e_lgp_nat_pdp_var; de_n = de_lgp_nat_pdp_var;
            norm_e_sq = sum(e_lgp_nat_pdp_var.^2) + sum(de_lgp_nat_pdp_var.^2);
        end
        
        V_lgp_nat_pdp_traj = zeros(size(t_n));
        mu_l_temp = zeros(1,length(t_n)); mu_u_temp = zeros(1,length(t_n)); alpha_t = zeros(size(t_n));
        
        for i = 1:length(t_n)
            q_i = q_n(:,i); dq_i = dq_n(:,i);
            e_i = e_n(:,i); de_i = de_n(:,i);
            
            if w_var == 0.1
                g_est = g_v_est_q_vec_t(:,i); g_est_e = g_v_est_e_vec_t(:,i); G_est_e = G_v_est_e_t(i); M_est = M_v_mat_t(:,:,i); C_est = C_v_mat_t(:,:,i);
                m_eig_vec_i = M_v_eig_vec_t(:,i);
                Sigma_lgp = Sigma_t(:,:,i);
            else
                if w_nat
                    g_est = g_est_q_vec_t(:,i); g_est_e = g_est_e_vec_t(:,i); G_est_e = G_est_e_t(i); M_est = M_mat_t(:,:,i); C_est = C_mat_t(:,:,i);
                    m_eig_vec_i = M_eig_vec_t(:,i);
                else
                    g_est = g_pdp_est_q_vec_t(:,i); M_est = M_pdp_mat_t(:,:,i); C_est = C_pdp_mat_t(:,:,i); m_eig_vec_i = M_eig_vec_pdp_t(:,i);
                end
            end
            D_est = diag(est_params.d); D_est_de = D_est;
            
            V_lgp_nat_pdp_traj(i) = 1/2*[e_i;de_i]' * [kappa_ct*eye(Ncc), eps_ct*M_est; eps_ct*M_est, M_est] * [e_i;de_i];
            if w_nat
                V_lgp_nat_pdp_traj(i) = V_lgp_nat_pdp_traj(i) + G_est_e;
            end
            
            if w_nat
                mu_l_temp(i) = 2*G_est_e/norm_e_sq(i) + (kappa_ct+m_eig_vec_i(1))/2 - 1/2*sqrt((kappa_ct-m_eig_vec_i(1))^2+(2*eps_ct*m_eig_vec_i(1))^2);
                mu_u_temp(i) = 2*G_est_e/norm_e_sq(i) + (kappa_ct+m_eig_vec_i(Ncc))/2 + 1/2*sqrt((kappa_ct-m_eig_vec_i(Ncc))^2+(2*eps_ct*m_eig_vec_i(Ncc))^2);
            else
                mu_l_temp(i) = (kappa_ct+m_eig_vec_i(1))/2 - 1/2*sqrt((kappa_ct-m_eig_vec_i(1))^2+(2*eps_ct*m_eig_vec_i(1))^2);
                mu_u_temp(i) = (kappa_ct+m_eig_vec_i(Ncc))/2 + 1/2*sqrt((kappa_ct-m_eig_vec_i(Ncc))^2+(2*eps_ct*m_eig_vec_i(Ncc))^2);
            end
            
            eta_e_i = heaviside(e_i'*g_est)*e_i'*g_est;
            d_est = D_est*dq_i;
            eta_de_i = heaviside(de_i'*d_est)*de_i'*d_est;
            omega_e_i = heaviside(e_i'*g_est)*de_i'*((e_i*e_i')./(1e-3+norm(e_i)^2))*g_est;
            omega_de_i = heaviside(de_i'*d_est)*e_i'*((de_i*de_i')./(1e-3+norm(de_i)^2))*d_est;
        
            % time-variant alpha calc
            D_tilde = w_nat*(D_est_de - min(diag(D_est_de))*eye(Ncc));
            R_mat = [zeros(Ncc), eps_ct/2*(D_tilde-C_est'); eps_ct/2*(D_tilde-C_est), D_tilde];
            kp = Kp(1,1);
            if w_var == 0.1
                K1 = k1*eye(Ncc); K2 = k2*eye(Ncc); K3 = k3*eye(Ncc);
                Kp_v = K1 * (eye(Ncc) - ( ( K3*(K2+Sigma_lgp)*K3 + K1 )\K1 ) );
                kp = min(eig(Kp + Kp_v));
                Kp_tilde = Kp + Kp_v - kp*eye(Ncc);
                R_mat = R_mat + [eps_ct*Kp_tilde, (1+eps_ct)/2*Kp_tilde; (1+eps_ct)/2*Kp_tilde, Kp_tilde];
            end
            r_min = min(eig(R_mat));
            gam_t = w_nat*min(diag(D_est_de)) + kp - phi_ct/2;
            b_t = kp + eps_ct*(w_nat*min(diag(D_est_de)) + kp) - kappa_ct;
            if w_nat == 1
                X = eps_ct * (e_i'*g_est_e+eta_e_i+omega_de_i) + eta_de_i + omega_e_i;
            else
                X = 0;
            end
        
            u0 = (eps_ct * (kp-the_ct/2) + gam_t - eps_ct*min(eig(M_est)))/2; u1 = (kappa_ct+min(eig(M_est)))/2;
            p0 = gam_t - eps_ct * (kp-the_ct/2 + min(eig(M_est))); p1 = kappa_ct - min(eig(M_est)); q = 2*eps_ct*min(eig(M_est));
            AA = u1;
            if w_nat
                AA = AA + G_est_e/(norm(e_i)^2+norm(de_i)^2);
            end
            BB = u0 + r_min + X/(norm(e_i)^2+norm(de_i)^2);
            A2 = 4*AA^2 - p1^2 - q^2; A1 = 2*(q*b_t - p0*p1 - 4*AA*BB); A0 = 4*BB^2 - p0^2 - b_t^2;
            alpha_min = 1/(2*A2) * (-A1 - sqrt(A1^2 - 4*A2*A0)); a_t = eps_ct*(kp - the_ct/2)-alpha_min*kappa_ct;
            u_min = (a_t+gam_t-(eps_ct+alpha_min)*min(eig(M_est)))/2; upsil_min = u_min - 1/2*sqrt((a_t - gam_t + (eps_ct+alpha_min)*min(eig(M_est)))^2 + (2*eps_ct*alpha_min*min(eig(M_est))-b_t)^2);
        
            u0 = (eps_ct * (kp-the_ct/2) + gam_t - eps_ct*max(eig(M_est)))/2; u1 = (kappa_ct+max(eig(M_est)))/2;
            p0 = gam_t - eps_ct * (kp-the_ct/2 + max(eig(M_est))); p1 = kappa_ct - max(eig(M_est)); q = 2*eps_ct*max(eig(M_est));
            AA = u1;
            if w_nat
                AA = AA + G_est_e/(norm(e_i)^2+norm(de_i)^2);
            end
            BB = u0 + r_min + X/(norm(e_i)^2+norm(de_i)^2);
            A2 = 4*AA^2 - p1^2 - q^2; A1 = 2*(q*b_t - p0*p1 - 4*AA*BB); A0 = 4*BB^2 - p0^2 - b_t^2;
            alpha_max = 1/(2*A2) * (-A1 - sqrt(A1^2 - 4*A2*A0)); a_t = eps_ct*(kp - the_ct/2)-alpha_max*kappa_ct;
            u_max = (a_t+gam_t-(eps_ct+alpha_max)*max(eig(M_est)))/2; upsil_max = u_max - 1/2*sqrt((a_t - gam_t + (eps_ct+alpha_max)*max(eig(M_est)))^2 + (2*eps_ct*alpha_max*max(eig(M_est))-b_t)^2);
        
            if upsil_min < upsil_max
                alpha_t(i) = alpha_min;
            else
                alpha_t(i) = alpha_max;
            end
            %
        end
        
        c_temp = c(num_sim+1-ctr,:);
        
        figure(3)
        subplot(311)
        fill([t_n, fliplr(t_n)], [1/2*(mu_l_temp.*norm_e_sq), fliplr(1/2*(mu_u_temp.*norm_e_sq))], c_temp,'EdgeColor','none','FaceAlpha',sh_alpha);
        V_t_mat{end+1} = V_lgp_nat_pdp_traj;
        t_mat{end+1} = t_n;
        
        % time-variant alpha
        alpha_t_int = zeros(1,length(t_n));
        for i = 2:length(t_n)
            alpha_t_int(i) = alpha_t_int(i-1) + trapz(t_n(i-1:i), alpha_t(i-1:i)');
        end
        if w_var == 0.1
            lgp_error_fac = 2 * max(mu_l_temp) / (1/phi_ct + eps_ct/the_ct);
        end
        V_border_ct = lgp_error_fac * (1/phi_ct + eps_ct/the_ct) / 4; %./ (4*alpha_t);
        x_V_ct_bound_time_var = sqrt(2*V_border_ct./mu_l_temp) + sqrt(2*(V_lgp_nat_pdp_traj(1)-V_border_ct)./mu_l_temp) .* exp(-alpha_t_int);
        
        subplot(312)
        fill([t_n, fliplr(t_n)], [2.5e-3*ones(size(t_n)), fliplr(x_V_ct_bound_time_var)],c_temp,'EdgeColor','none','FaceAlpha',sh_alpha);
        x_bnd_t_mat{end+1} = x_V_ct_bound_time_var;
        
        subplot(311)
        grid on;
        ylabel('$V(\mbox{\boldmath $e$}(t),\dot{\mbox{\boldmath $e$}}(t),t)$ (J)','Interpreter','Latex','FontSize',fs)
        
        subplot(312)
        x_t_mat{end+1} = sqrt(norm_e_sq);
        hndls(ctr) = semilogy(t_n,sqrt(norm_e_sq),'LineWidth',mlw, 'Color', c_temp);
        ctr = ctr + 1;
        grid on;
        ylabel('$||[\mbox{\boldmath $e$}^\top \, s\dot{\mbox{\boldmath $e$}}^\top]||$ (rad)','Interpreter','Latex','FontSize',fs)
        
        subplot(3,1,3)
        semilogy(t_n, alpha_t, 'LineWidth', mlw, 'Color', c_temp)
        hold on
        ylim([2e-2 1e2])
    
    end

end

%%
tp_end = t_end;
figure(3); subplot(311);
for id = 1:3
    semilogy(t_mat{id},V_t_mat{id},'LineWidth',mlw, 'Color',c(num_sim+1-id,:));
end
ylim([1e-8 120]); xlim([0 t_end]); yticks([1e-5 1e0])

subplot(312)
for id = 1:3
    semilogy(t_mat{id},x_bnd_t_mat{id},'--','LineWidth',mlw, 'Color',c(num_sim+1-id,:))
    semilogy(t_mat{id},x_t_mat{id},'LineWidth',mlw, 'Color',c(num_sim+1-id,:))
end

if ctr == 4
    legend(hndls,'L-GP-PD+','L-GP nat-PD+','L-GP var-nat-PD+','Interpreter','Latex','FontSize',fs,'Location','Best','Orientation','Horizontal')
end
ylim([7e-3 5e3]); xlim([0 tp_end]); yticks([1e-2 1e0 1e2])

subplot(3,1,3)
set(gca,'yscale','log','FontSize',ax_fs); grid on;
ylabel('$\alpha(t)$ (1/s)','Interpreter','Latex','FontSize',fs)
xlim([0 tp_end])
xlabel('$t$ (s)','Interpreter','Latex','FontSize',fs)

%%
fig11 = figure(11); clf; hold on; grid; set(gca,'FontSize',12);

t_des = linspace(-pi/2,pi/2,1e3); lxy_fem = zeros(2,N,length(t_des));
for i = 1:length(t_des)
    lxy_fem(:,:,i) = fem_li_xy_n100(1/N,kron(q_des(t_des(i))./(N/Ncc), ones(N/Ncc,1)));
end
plot(squeeze(lxy_fem(1,end,:)),squeeze(lxy_fem(2,end,:)),'--','Color',c(2,:),'LineWidth',mlw)

t_des = linspace(0,2*pi,2^5+1); lxy_fem = zeros(2,N,length(t_des));
for i = 1:length(t_des)
    lxy_fem(:,:,i) = fem_li_xy_n100(1/N,kron(q_des(t_des(i))./(N/Ncc), ones(N/Ncc,1)));
    plot(lxy_fem(1,:,i),lxy_fem(2,:,i),'Color',c(1,:))
    if i == 1
        text(lxy_fem(1,end,i)+0.015,lxy_fem(2,end,i),'$t=0,\pi$','Interpreter','Latex','FontSize',fs)
    elseif i == 3
        text(lxy_fem(1,end,i)+0.015,lxy_fem(2,end,i)+0.015,'$t=\frac{\pi}{8},\frac{7\pi}{8}$','Interpreter','Latex','FontSize',fs)
    elseif i == 9
        text(lxy_fem(1,end,i)-0.125,lxy_fem(2,end,i),'$t=\frac{\pi}{2}$','Interpreter','Latex','FontSize',fs)
    elseif i == 19
        text(lxy_fem(1,end,i)+0.015,lxy_fem(2,end,i)-0.015,'$t=\frac{9\pi}{8},\frac{15\pi}{8}$','Interpreter','Latex','FontSize',fs)
    elseif i == 27
        text(lxy_fem(1,end,i)-0.175,lxy_fem(2,end,i),'$t=\frac{3\pi}{2}$','Interpreter','Latex','FontSize',fs)
    end
end
xlabel('$x$ (m)','Interpreter','Latex','FontSize',fs)
ylabel('$y$ (m)','Interpreter','Latex','FontSize',fs)
