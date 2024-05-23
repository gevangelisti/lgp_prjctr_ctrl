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

%%
clear
close all
clc

%% Parametrization
train = 0;
rng(0)
N = 100;
Ncc = 4;
best = 1;
len_sim = 1e5+1;
tsim = 4;
d = 5;
k = 10;
a = 20*5/N;

%% Simulate
x0 = zeros(2*N,1);
soro = plnr_so_ro(N);

t_fem = linspace(0,tsim,len_sim);
[qout_full,dqout_full] = fem_sim_so_ro(soro,N,k,d,x0,t_fem,a);
[qout_full_2,dqout_full_2] = fem_sim_so_ro(soro,N,k,d,x0,t_fem,a/2);
showmotion(plnr_so_ro(N),t_fem,qout_full)

%%
c = get(gca,'colororder');
close all
fs = 14; % font size
ax_fs = 14; % axis font size
mlw = 1.5; % main line width
msz = 12; % marker size

figure
hold on
ind_vec = 1:(len_sim-1)/4:len_sim;
t_p_vec = t_fem(ind_vec);

lxy_fem = zeros(2,N,length(ind_vec));

for i = 1:length(t_p_vec)
    lxy_fem(:,:,i) = fem_li_xy_n100(1/N,qout_full(:,ind_vec(i)));
    plot(lxy_fem(1,:,i),lxy_fem(2,:,i),'LineWidth',mlw)
end
legend(['$t = ',num2str(t_p_vec(1)),'$ s'],['$t = ',num2str(t_p_vec(2)),'$ s'],['$t = ',num2str(t_p_vec(3)),'$ s'],['$t = ',num2str(t_p_vec(4)),'$ s'],['$t = ',num2str(t_p_vec(5)),'$ s'], 'Interpreter', 'Latex','FontSize',ax_fs,'Location','Best')
xlim([0 1])
ylim([-0.05 0.45])
xlabel('x (m)', 'Interpreter', 'Latex','FontSize',ax_fs)
ylabel('y (m)', 'Interpreter', 'Latex','FontSize',ax_fs)
set(gca,'FontSize',ax_fs)

%% Simulate PCC model
a_pcc = a*ones(Ncc,1);
L = 1/Ncc*ones(Ncc,1);
mu = 1/Ncc*[1.25;0.75;1.25;0.75];
Iz = 1/12*mu.*L.^2;
d_pcc = d/(N/Ncc)*ones(Ncc,1);
k_pcc = k/(N/Ncc)*[1.25;0.75;1.25;0.75];

test_sim_pcc;
t_init = t;
q_init = q_mat;
title('Initial parametrization PCC model')

%% Simulate measurement
num_intrvl = 100;
start_ind = 10;
ind_meas = 1:(len_sim-1)/num_intrvl:len_sim;
ind_meas = ind_meas(start_ind:end);
t_meas = t_fem(ind_meas);
[q_meas,dq_meas,ddq_meas] = get_states_from_fem(N,Ncc,qout_full,dqout_full,soro,k,d,a,ind_meas);
sigma_tau = 0.01;
tau_meas = a * (1 + sigma_tau*randn(numel(q_meas),1));

t_id = t;
q_id = q_mat;

%% Train & Test L-GP

% Set training data set
ns = 4; ss = 8; se = 0;
ind_tr = 1:(len_sim-1)/num_intrvl:len_sim;
ind_tr = ind_tr(1+ss:end-se);
t_tr = t_fem(ind_tr);
[q_tr,dq_tr,ddq_tr] = get_states_from_fem(N,Ncc,qout_full,dqout_full,soro,k,d,a,ind_tr);
q_tr = q_tr(:,1:ns:end);
dq_tr = dq_tr(:,1:ns:end);
X_tr = [q_tr',dq_tr'];
ddq_tr = ddq_tr(:,1:ns:end);
tau_tr = a * (1 + sigma_tau*randn(numel(q_tr),1));
y_D = tau_tr;
D = size(X_tr,1);

ind_v = ind_tr(1)+(ind_tr(2)-ind_tr(1))/2:(len_sim-1)/(1*num_intrvl):len_sim;
[q_v,dq_v,ddq_v] = get_states_from_fem(N,Ncc,qout_full,dqout_full,soro,k,d,a,ind_v);
X_v = [q_v',dq_v'];
tau_v = a * ones(size(ddq_v'));

% Set initial parameter estimates
est_params.mu = 1/Ncc*ones(Ncc,1);
est_params.L = 1/Ncc*ones(Ncc,1);
est_params.d = d_pcc;
est_params.k = k_pcc;
est_params.g = 9.81;
est_params.Iz = 1/12/Ncc^3*ones(Ncc,1);

% Compute mean a-priori estimates
mu_D = zeros(Ncc*D,1);
g_est_real = zeros(size(ddq_tr));
for i = 1:D
    [mu_D(1+(i-1)*Ncc:i*Ncc),g_est_real(:,i)] = ivd_dyn_pcc_N(est_params,q_tr(:,i),dq_tr(:,i),ddq_tr(:,i));
end

% Noise matrices
sigma_n = 0;
eps_sq = 0;
Z = diag(sigma_n^2*ones(1,Ncc));
Upsilon = diag(sigma_tau^2*ones(1,Ncc));
    
% Train
options = optimoptions('ga','UseParallel', true);
nonlcon = [];
[q_fem_full,~,~] = get_states_from_fem(N,Ncc,qout_full,dqout_full,soro,k,d,a);
[q_fem_full_2,~,~] = get_states_from_fem(N,Ncc,qout_full_2,dqout_full_2,soro,k,d,a/2);
ds = 100;

cf_test = @(x) cf_step_sim(X_tr, ddq_tr', y_D, Z, Upsilon, est_params, t_fem(1:ds:end), q_fem_full(:,1:ds:end), a_pcc, t_fem(1:ds:end), q_fem_full_2(:,1:ds:end), mu_D, D, Ncc, x);

lb = [zeros(1,10),     1e0*ones(1,6)];
ub = [0.25*ones(1,10), 2e1, 2e1, 4e1, 2e1, 4e1, 2e1];
if train
    xp_test = ga(cf_test,length(lb),[],[],[],[],lb,ub,nonlcon,[],options);
else
    load('soro_sim_data.mat','xp_test');
end

% Parameter parsing
psi = [sigma_n, sigma_tau, xp_test];
[LLT_Sigma_sq,lambda_isq,rho_g_sq,P_g_inv2,P_g_sq,~,~] = get_hyp_params(Ncc,psi);

% Pre-computations: compute covariance matrices
K_y_plus_Pi = calc_Kyy_pcc_mex(X_tr, ddq_tr', LLT_Sigma_sq, lambda_isq, rho_g_sq, P_g_inv2, Z, Upsilon, est_params, D, Ncc);
K_y0 = calc_K_y0_pcc_mex(X_tr, rho_g_sq, P_g_inv2, D, Ncc);
k_0_inv = 1/(2*rho_g_sq);
A = K_y_plus_Pi + eps_sq*eye(length(y_D)) - k_0_inv*(K_y0*K_y0');

% compute inversion using Cholesky factorization
L_A = chol(A,'lower');
A_inv_times_y_min_mu = L_A'\(L_A\(y_D-mu_D)); % the backslash operator recognizes triangular systems!
K_times_a = k_0_inv * (K_y0' * A_inv_times_y_min_mu);

%%
ode_fun_GP = @(t,x) [...
    x(Ncc+1:2*Ncc);
    ddq_lgp_pcc_mex(x(1:Ncc), x(Ncc+1:2*Ncc), a_pcc, X_tr, ddq_tr', LLT_Sigma_sq, lambda_isq, rho_g_sq, P_g_inv2, A_inv_times_y_min_mu, K_times_a, est_params, D, Ncc) ...
    ];

[~,x_pcclgp_sim] = ode15s(ode_fun_GP,t_fem,zeros(2*Ncc,1));

%%
figure
qfem_i_rs2id = zeros(size(q_id));

for i = 1:Ncc
    subplot(2,Ncc,i)
    hold on
    plot(t_id,180/pi*q_id(:,i),'LineWidth',mlw,'Color',c(3,:))
    plot(t_fem,180/pi*x_pcclgp_sim(:,i),'LineWidth',mlw,'Color',c(1,:))
    plot(t_fem,180/pi*q_fem_full(i,:),'--','LineWidth',mlw,'Color',c(2,:))
    grid
    title(['$n=',num2str(i),'$'],'Interpreter','Latex')
    xlabel('$t$ (s)','Interpreter','Latex')
    ylabel(['$q_',num2str(i),'$ (deg)'],'Interpreter','Latex')
    set(gca,'FontSize',ax_fs)

    subplot(2,Ncc,Ncc+i)
    hold on
    qfem_i_rs2id(:,i) = interp1(t_fem,q_fem_full(i,:),t_id);
    plot(t_id,180/pi*(q_id(:,i)-qfem_i_rs2id(:,i)),'Color',c(3,:),'LineWidth',mlw)
    plot(t_fem,180/pi*(x_pcclgp_sim(:,i)-q_fem_full(i,:)'),'Color',c(1,:),'LineWidth',mlw)
    grid
    xlabel('$t$ (s)', 'Interpreter', 'Latex','FontSize',ax_fs)
    ylabel(['$e_',num2str(i),'$ (deg)'],'Interpreter','Latex','FontSize',ax_fs)
    set(gca,'FontSize',ax_fs)
end
subplot(244)
legend('PCC','L-GP','FEM', 'Interpreter', 'Latex','Location','Best')