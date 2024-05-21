%% Open-source implementation of Lagrangian-Gaussian Processes presented in:
% [1] Evangelisti et al., (2024). Exponentially Stable Projector-based Control of Lagrangian Systems with Gaussian Processes, ArXiv
% [2] Evangelisti and Hirche, (2022). Physically Consistent Learning of Conservative Lagrangian Systems with Gaussian Processes, Conference on Decision and Control (CDC), https://doi.org/10.1109/CDC51059.2022.9993123 
% [3] Evangelisti and Hirche, (2024). Data-Driven Momentum Observers With Physically Consistent Gaussian Processes, Transactions on Robotics (T-RO), https://doi.org/10.1109/TRO.2024.3366818 
% 
% Note: Run start_lgp before first execution!
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
clc
close all
clear

%% Parametrization
rng(0);
N = 2;
tend = 20;
F_s = 1e3;
constrain_M_0 = 0;

% system parameters
m = [1, 1]; % masses [kg]
l = [1, 1]; % massless rod lengths [m]
g = 10; % force field (=gravitational) acceleration [m/s^2]
d = [1, 1]; % dampers [kg m^2/s]
Iz1 = 1/3*m(1)*l(1)^2;
Iz2 = 1/3*m(2)*l(2)^2;
alpha = Iz1 + Iz2 + m(1)*(l(1)/2)^2 + m(2)*(l(1)^2 + (l(2)/2)^2);
beta = m(2)*l(1)*l(2)/2;
delta = Iz2 + m(2)*(l(2)/2)^2;

%% noise
sigma_n = pi/180; % acceleration measurement std [1/s^2]
sigma_tau = 0.1; % torque measurement std [Nm]

% covariance matrices
Z = diag(sigma_n^2*ones(1,2));
Upsilon = diag(sigma_tau^2*ones(1,2));

%% Real system parametrization & analytical descriptions

true_params.m1 = m(1);
true_params.m2 = m(2);
true_params.l1 = l(1);
true_params.l2 = l(2);
true_params.d = d;
true_params.g = g;
true_params.alpha = alpha;
true_params.beta = beta;
true_params.delta = delta;

% forward kinematics
xp1 = @(q1) l(1)*cos(q1);
yp1 = @(q1) l(1)*sin(q1);
xp2 = @(q1,q2) xp1(q1) + l(2)*cos(q1+q2);
yp2 = @(q1,q2) yp1(q1) + l(2)*sin(q1+q2);

% Lagrange derivatives
L_dq1_dq1 = @(q1,q2) alpha + 2*beta*cos(q2);
L_dq1_dq2 = @(q1,q2) delta + beta*cos(q2);
L_dq2_dq2 = @(q1,q2) delta*ones(size(q1));
L_dq1_q1 = @(q1,q2,dq1,dq2) zeros(size(q1));
L_dq1_q2 = @(q1,q2,dq1,dq2) -beta*sin(q2).*(2*dq1+dq2);
L_dq2_q1 = @(q1,q2,dq1,dq2) zeros(size(q1));
L_dq2_q2 = @(q1,q2,dq1,dq2) -beta*sin(q2).*dq1;
T_q1 = @(q1,q2,dq1,dq2) zeros(size(q1));
T_q2 = @(q1,q2,dq1,dq2) -beta*sin(q2).*dq1.*(dq1+dq2);
V_q1 = @(q1,q2) (m(1)*l(1)/2 + m(2)*l(1))*g*sin(q1) + m(2)*l(2)/2*g*sin(q1+q2);
V_q2 = @(q1,q2) m(2)*l(2)/2*g*sin(q1+q2);
V = @(q) ...
    (m(1)*l(1)/2 + m(2)*l(1))*g*(1-cos(q(1))) + ...
    m(2)*l(2)/2*g*(1-cos(q(1)+q(2)));

% Dynamics
M = @(q) [L_dq1_dq1(q(1),q(2)), L_dq1_dq2(q(1),q(2)); L_dq1_dq2(q(1),q(2)), L_dq2_dq2(q(1),q(2))];
C = @(q,dq) [-beta*sin(q(2))*dq(2), -beta*sin(q(2))*(dq(1)+dq(2)); beta*sin(q(2))*dq(1), 0];
G = @(q) [V_q1(q(1),q(2)); V_q2(q(1),q(2))];
Dmp = @(dq) [d(1)+d(2)*abs(dq(1)), 0; 0, d(1)+d(2)*abs(dq(2))];

%% Estimated model parametrization & analytical descriptions
m1h = 1.5;
m2h = 0.5;
l1h = 1.5;
l2h = 0.5;
dh = [0.5, 1.5];
gh = g;
Iz1h = 1/3*m1h*l1h^2;
Iz2h = 1/3*m2h*l2h^2;
alphah = Iz1h + Iz2h + m1h*(l1h/2)^2 + m2h*(l1h^2 + (l2h/2)^2);
betah = m2h*l1h*l2h/2;
deltah = Iz2h + m2h*(l2h/2)^2;

est_params.m1h = m1h;
est_params.m2h = m2h;
est_params.l1h = l1h;
est_params.l2h = l2h;
est_params.dh = dh;
est_params.gh = gh;
est_params.alphah = alphah;
est_params.betah = betah;
est_params.deltah = deltah;

% Lagrange derivatives
hat_d2_f_GP_dxidchi = @(dq,q) [0, -betah*sin(q(2)).*(2*dq(1)+dq(2)); 0, -betah*sin(q(2)).*dq(1)];
hat_nabla_chi_f_GP = @(dq,q) [0; -betah*sin(q(2)).*dq(1).*(dq(1)+dq(2))];
hat_nabla_chi_g_GP = @(q) [(m1h*l1h/2 + m2h*l1h)*gh*sin(q(1)) + m2h*l2h/2*gh*sin(q(1)+q(2)); m2h*l2h/2*gh*sin(q(1)+q(2))];

% Dynamics
hatM = @(q) [alphah + 2*betah*cos(q(2)), deltah + betah*cos(q(2)); deltah + betah*cos(q(2)), deltah];
hatC = @(q,dq) [-betah*sin(q(2))*dq(2), -betah*sin(q(2))*(dq(1)+dq(2)); betah*sin(q(2))*dq(1), 0];
hatG = @(q) [(m1h*l1h/2 + m2h*l1h)*gh*sin(q(1)) + m2h*l2h/2*gh*sin(q(1)+q(2)); m2h*l2h/2*gh*sin(q(1)+q(2))];
hatD = @(dq) [dh(1)+dh(2)*abs(dq(1)), 0; 0, dh(1)+dh(2)*abs(dq(2))];

%% Desired trajectory
amp1 = pi/2;

qd1 = @(t) amp1*sin(t);
dqd1 = @(t) amp1*cos(t);
ddqd1 = @(t) -amp1*sin(t);

amp2 = amp1;
omega2 = 1;
qd2 = @(t) amp2*sin(omega2*t);
dqd2 = @(t) amp2*omega2*cos(omega2*t);
ddqd2 = @(t) -amp2*omega2^2*sin(omega2*t);

q_d = @(t) [qd1(t); qd2(t)];
dq_d = @(t) [dqd1(t); dqd2(t)];
ddq_d = @(t) [ddqd1(t); ddqd2(t)];

%% Classical PD+ control law
Kp = 10*eye(2);
Kd = 10*eye(2);
tau_pdp = @(t,q,dq) hatM(q)*ddq_d(t) + hatC(q,dq)*dq_d(t) + hatG(q) + hatD(dq)*dq - Kp*(q-q_d(t)) - Kd*(dq-dq_d(t));

%% Parametric model-based structure-preserving controller
ukp = 1; ukd = 1;
tau_sc_param = @(t,q,dq) hatM(q)*ddq_d(t) + hatC(q,dq)*dq_d(t) ...
    + (eye(N)-heaviside((q-q_d(t))'*hatG(q))*(q-q_d(t))*(q-q_d(t))'./(1e-3 + norm((q-q_d(t)))^2) )*hatG(q) - hatG(q-q_d(t)) - ukp*Kp*(q-q_d(t)) ... 
    + (eye(N)-heaviside((dq-dq_d(t))'*hatD(dq)*dq)*(dq-dq_d(t))*(dq-dq_d(t))'./(1e-3 + norm((dq-dq_d(t)))^2) )*hatD(dq)*dq - hatD(dq-dq_d(t))*(dq-dq_d(t)) - ukd*Kd*(dq-dq_d(t));

%% Simulation with classical PD+
ode_fun = @(t,x) [...
    x(3); ...
    x(4); ...
    M(x(1:2)) \ (tau_pdp(t,x(1:2),x(3:4)) - C(x(1:2),x(3:4))*x(3:4) - G(x(1:2)) - Dmp(x(3:4))*x(3:4)) ...
    ];
ode_fun_sc = @(t,x) [x(3); x(4); M(x(1:2)) \ (tau_sc_param(t,x(1:2),x(3:4)) - C(x(1:2),x(3:4))*x(3:4) - G(x(1:2)) - Dmp(x(3:4))*x(3:4))];

% Simulate
x0 = [pi/4;pi/4;zeros(2,1)];
[~,x] = ode45(ode_fun,linspace(0,tend,F_s*tend+1),x0);
[~,x_scp] = ode45(ode_fun_sc,linspace(0,tend,F_s*tend+1),x0);
t = linspace(0,tend,F_s*tend+1);

% trajectories
q1 = x(:,1);
q2 = x(:,2);
dq1 = x(:,3);
dq2 = x(:,4);
u1 = zeros(size(t));
u2 = zeros(size(t));
u1_sc = zeros(size(t));
u2_sc = zeros(size(t));
ddq1 = zeros(size(q1));
ddq2 = zeros(size(q1));
for i = 1:length(q1)
    u = tau_pdp(t(i),x(i,1:2)',x(i,3:4)');
    u1(i) = u(1);
    u2(i) = u(2);
    ddq = M([q1(i),q2(i)]) \ (u - C([q1(i),q2(i)],[dq1(i),dq2(i)])*[dq1(i);dq2(i)] - G([q1(i),q2(i)]) - Dmp([dq1(i);dq2(i)])*[dq1(i);dq2(i)]);
    ddq1(i) = ddq(1);
    ddq2(i) = ddq(2);
    u = tau_sc_param(t(i),x_scp(i,1:2)',x_scp(i,3:4)');
    u1_sc(i) = u(1);
    u2_sc(i) = u(2);
end

%% Simulate measurement

% training points
q_vec = linspace(-1,1,5);
dq_vec = linspace(-1,1,3);
q_vec_g = linspace(-1.25,1.25,6);
[q1_nd_1,q2_nd_1,dq1_nd_1,dq2_nd_1,ddq1_nd_1,ddq2_nd_1] = ndgrid(q_vec,q_vec,-1,1,4,4);
[q1_nd_2,q2_nd_2,dq1_nd_2,dq2_nd_2,ddq1_nd_2,ddq2_nd_2] = ndgrid(0,0,dq_vec,dq_vec,0,0);
[q1_nd_g,q2_nd_g,dq1_nd_g,dq2_nd_g,ddq1_nd_g,ddq2_nd_g] = ndgrid(q_vec_g,q_vec_g,1.5,0,0,0);
q1_tr = [q1_nd_1(:); q1_nd_2(:); q1_nd_g(:)];
q2_tr = [q2_nd_1(:); q2_nd_2(:); q2_nd_g(:)];
dq1_tr = [dq1_nd_1(:); dq1_nd_2(:); dq1_nd_g(:)];
dq2_tr = [dq2_nd_1(:); dq2_nd_2(:); dq2_nd_g(:)];
ddq1_tr = [ddq1_nd_1(:); ddq1_nd_2(:); ddq1_nd_g(:)];
ddq2_tr = [ddq2_nd_1(:); ddq2_nd_2(:); ddq2_nd_g(:)];
ddq1_meas = ddq1_tr;
ddq2_meas = ddq2_tr;

% compute torques
q_tr = [q1_tr,q2_tr];
dq_tr = [dq1_tr,dq2_tr];
u1_tr = zeros(1,length(q1_tr));
u2_tr = zeros(1,length(q2_tr));
for i = 1:length(q1_tr)
    tau = M([q1_tr(i),q2_tr(i)]) * [ddq1_tr(i);ddq2_tr(i)] + ...
        C([q1_tr(i),q2_tr(i)],[dq1_tr(i),dq2_tr(i)]) * [dq1_tr(i);dq2_tr(i)] + ...
        G([q1_tr(i),q2_tr(i)]) + ...
        Dmp([dq1_tr(i);dq2_tr(i)]) * [dq1_tr(i);dq2_tr(i)];
    u1_tr(i) = tau(1);
    u2_tr(i) = tau(2);
end

if eig(Z) == zeros(2,1)
    L_zeta = zeros(2);
else
    L_zeta = chol(Z)';
end
if eig(Upsilon) == zeros(2,1)
    L_upsilon = zeros(2);
else
    L_upsilon = chol(Upsilon)';
end
for i = 1:length(q1_tr)
    % ddq measurement
    zeta = L_zeta*randn(2,1);
    ddq1_meas(i) = ddq1_meas(i) + zeta(1);
    ddq2_meas(i) = ddq2_meas(i) + zeta(2);
    % tau measurement
    upsilon = L_upsilon*randn(2,1);
    u1_tr(i) = u1_tr(i) + upsilon(1);
    u2_tr(i) = u2_tr(i) + upsilon(2);
end

ddq_meas = [ddq1_meas, ddq2_meas];

D = length(ddq1_meas) - numel(q1_nd_g);
q1_v = q1_tr(D+1:end);
q1_tr = q1_tr(1:D);
q2_v = q2_tr(D+1:end);
q2_tr = q2_tr(1:D);
dq1_v = dq1_tr(D+1:end);
dq1_tr = dq1_tr(1:D);
dq2_v = dq2_tr(D+1:end);
dq2_tr = dq2_tr(1:D);
ddq1_v = ddq1_tr(D+1:end);
ddq1_tr = ddq1_tr(1:D);
ddq2_v = ddq2_tr(D+1:end);
ddq2_tr = ddq2_tr(1:D);
ddq1_meas_v = ddq1_meas(D+1:end);
ddq1_meas = ddq1_meas(1:D);
ddq2_meas_v = ddq2_meas(D+1:end);
ddq2_meas = ddq2_meas(1:D);
ddq_meas_v = ddq_meas(D+1:end,:);
ddq_meas = ddq_meas(1:D,:);
q_tr = q_tr(1:D,:);
dq_tr = dq_tr(1:D,:);
q_v = [q1_v,q2_v];
dq_v = [dq1_v,dq2_v];
u1_v = u1_tr(D+1:end);
u1_tr = u1_tr(1:D);
u2_v = u2_tr(D+1:end);
u2_tr = u2_tr(1:D);

mu_D = zeros(N*D,1);

X_tr = [q_tr,dq_tr];
X_v = [q_v,dq_v];
for i = 1:D
    mu_D(1+(i-1)*N:i*N) = hatM(X_tr(i,1:2)')*ddq_meas(i,:)' + ...
        hatC(X_tr(i,1:2)',X_tr(i,3:4)')*X_tr(i,3:4)' + ...
        hatG(X_tr(i,1:2)') + hatD(X_tr(i,3:4)')*X_tr(i,3:4)';
end

y_D = zeros(N*length(ddq1_meas),1);
y_D(1:2:end) = u1_tr;
y_D(2:2:end) = u2_tr;

%% Hyperparameter optimization
eps_sq = sqrt(0.1);
options = optimoptions('fmincon','Display','Iter');
nonlcon = [];

x0_test =   [1e0*ones(1,3), 1e1, 1e1,    1,    1,  1e0*ones(1,2), 1e0*ones(1,2)];
lb =        [   zeros(1,3), 1e0, 1e0, 1e-1, 1e-1,     zeros(1,2),   0*ones(1,2)];
ub =        [1e1*ones(1,3), 1e5, 1e2,  1e1,  1e1,  1e1*ones(1,2), 1e5*ones(1,2)];
cf_test = @(x) cost_func_sqe(X_tr, ddq_meas, y_D, X_v, ddq_meas_v, [u1_v', u2_v'], Z, Upsilon, est_params, mu_D, D, N, eps_sq, x);

disp('Starting local least-squares-based hyperparameter optimization ...')
xp_test = fmincon(cf_test,x0_test,[],[],[],[],lb,ub,nonlcon,options);
disp('Finished local least-squares-based hyperparameter optimization')

%% Output final cost function value
cf_test(xp_test)

%% Parameter parsing
psi = [sigma_n, sigma_tau, [xp_test(1:8), 0, xp_test(9:end)]];
[LLT_Sigma_sq,lambda_isq,rho_g_sq,P_g_inv2,P_g_sq,L_Sigma_sq_dmp,lambda_isq_dmp] = get_hyp_params(N,psi);

%% Pre-computations

% compute covariance matrices
K_y_plus_Pi = calc_K_yy_mex(X_tr, ddq_meas, LLT_Sigma_sq, lambda_isq, rho_g_sq, P_g_inv2, L_Sigma_sq_dmp, lambda_isq_dmp, Z, Upsilon, est_params, D, N);
K_y0 = calc_K_y0_mex(X_tr, rho_g_sq, P_g_inv2, D, N);
K_0_inv = (1/rho_g_sq) * [1, zeros(1,N); zeros(N,1), P_g_sq];
A = K_y_plus_Pi + eps_sq*eye(length(y_D)) - K_y0*K_0_inv*K_y0';

% compute inversion using Cholesky factorization
L_A = chol(A,'lower');
A_inv_times_y_min_mu = L_A'\(L_A\(y_D-mu_D)); % the backslash operator recognizes triangular systems!
K_times_a = K_0_inv * (K_y0' * A_inv_times_y_min_mu);

%%
tau_tr = [y_D(1:2:end),y_D(2:2:end)];
[tau_est_tr,~,~,~,~,~,~] = ...
    calc_lgp_ests_mex(X_tr, ddq_meas, LLT_Sigma_sq, lambda_isq, rho_g_sq, P_g_inv2, L_Sigma_sq_dmp, lambda_isq_dmp, X_tr, ddq_meas, tau_tr, est_params, A_inv_times_y_min_mu, K_times_a, D, N);
opt_mean_tau_error = mean(mean(abs(tau_est_tr-tau_tr)))

%% Plots

fs = 14; % font size
ax_fs = 10; % axis font size
mlw = 1.5; % main line width
msz = 12; % marker size

c = get(gca,'colororder');
close all

%% static GP-based PD+ control law
Kp_s = Kp;
Kd_s = Kd;
Kp_v = Kp;
Kd_v = Kp;
tau_lgp_pdp_const = @(t,q,dq) tau_lgp_pdp_mex(t, q, dq, pi/2, 1, Kp_s, Kd_s, X_tr, ddq_meas, ...
    LLT_Sigma_sq, lambda_isq, rho_g_sq, P_g_inv2, L_Sigma_sq_dmp, lambda_isq_dmp, A_inv_times_y_min_mu, K_times_a, est_params, D, N);
tau_lgp_nat_pdp_const = @(t,q,dq) tau_lgp_nat_pdp_mex(t, q, dq, pi/2, 1, ukp*Kp_s, ukd*Kd_s, X_tr, ddq_meas, ... 
    LLT_Sigma_sq, lambda_isq, rho_g_sq, P_g_inv2, L_Sigma_sq_dmp, lambda_isq_dmp, A_inv_times_y_min_mu, K_times_a, est_params, D, N);
tau_lgp_nat_pdp_var = @(t,q,dq) tau_lgp_var_nat_pdp_mex(t, q, dq, pi/2, 1, ...
    1e2*eye(2), 2e-2*eye(2), 1/sqrt(2e-2*(1/1e0-1/1e2))*eye(2), ukp*Kp_s, ukd*Kd_s, A, ... 
    X_tr, ddq_meas, LLT_Sigma_sq, lambda_isq, rho_g_sq, P_g_inv2, L_Sigma_sq_dmp, lambda_isq_dmp, A_inv_times_y_min_mu, K_times_a, est_params, D, N);

%% PD+: static-GP-gain simulation
ode_robot = @(t,x) [x(3); x(4); M(x(1:2)) \ (tau_lgp_pdp_const(t,x(1:2),x(3:4)) - C(x(1:2),x(3:4))*x(3:4) - G(x(1:2)) - Dmp(x(3:4))*x(3:4))];
disp('Starting static L-GP-PD+ simulation ...')
[t_lgp_pdp_c,x_lgp_pdp_c] = ode45(ode_robot,linspace(0,tend,F_s*tend+1),x0);%[0,tend],x0);
disp('done')

% trajectories
q1_lgp_pdp_c = x_lgp_pdp_c(:,1); q2_lgp_pdp_c = x_lgp_pdp_c(:,2); dq1_lgp_pdp_c = x_lgp_pdp_c(:,3); dq2_lgp_pdp_c = x_lgp_pdp_c(:,4);

%% Proposed structure-preserving tracking control
ode_robot_nat = @(t,x) [x(3); x(4); M(x(1:2)) \ (tau_lgp_nat_pdp_const(t,x(1:2),x(3:4)) - C(x(1:2),x(3:4))*x(3:4) - G(x(1:2)))];
disp('Starting static nat L-GP-PD+ simulation ...')
[t_lgp_nat_pdp_c,x_lgp_nat_pdp_c] = ode45(ode_robot_nat,linspace(0,tend,F_s*tend+1),x0);%[0,tend],x0);
disp('done')

% trajectories
q1_lgp_nat_pdp_c = x_lgp_nat_pdp_c(:,1); q2_lgp_nat_pdp_c = x_lgp_nat_pdp_c(:,2); dq1_lgp_nat_pdp_c = x_lgp_nat_pdp_c(:,3); dq2_lgp_nat_pdp_c = x_lgp_nat_pdp_c(:,4);

%% Test
ode_robot_nat_var = @(t,x) [x(3); x(4); M(x(1:2)) \ (tau_lgp_nat_pdp_var(t,x(1:2),x(3:4)) - C(x(1:2),x(3:4))*x(3:4) - G(x(1:2)))];
disp('Starting variance-adaptive nat-L-GP-PD+ simulation ...')
[t_lgp_var_nat_pdp,x_lgp_var_nat_pdp] = ode45(ode_robot_nat_var,linspace(0,tend,F_s*tend+1),x0);
disp('done')

% trajectories
q1_lgp_nat_pdp_v = x_lgp_var_nat_pdp(:,1); q2_lgp_nat_pdp_v = x_lgp_var_nat_pdp(:,2); dq1_lgp_nat_pdp_v = x_lgp_var_nat_pdp(:,3); dq2_lgp_nat_pdp_v = x_lgp_var_nat_pdp(:,4);

%% Compute errors
disp('Analyzing performance metrics ...')
t_pdp = t;
e1_pdp = q1-qd1(t_pdp)';
e2_pdp = q2-qd2(t_pdp)';
e1_nat_pdp = x_scp(:,1)-qd1(t_pdp)';
e2_nat_pdp = x_scp(:,2)-qd2(t_pdp)';
e1_lgp_pdp = q1_lgp_pdp_c-qd1(t_lgp_pdp_c);
e2_lgp_pdp = q2_lgp_pdp_c-qd2(t_lgp_pdp_c);
e1_lgp_nat_pdp = q1_lgp_nat_pdp_c-qd1(t_lgp_nat_pdp_c);
e2_lgp_nat_pdp = q2_lgp_nat_pdp_c-qd2(t_lgp_nat_pdp_c);
e1_lgp_nat_pdp_v = q1_lgp_nat_pdp_v-qd1(t_lgp_var_nat_pdp);
e2_lgp_nat_pdp_v = q2_lgp_nat_pdp_v-qd2(t_lgp_var_nat_pdp);
de1_pdp = dq1-dqd1(t_pdp)';
de2_pdp = dq2-dqd2(t_pdp)';
de1_nat_pdp = x_scp(:,3)-dqd1(t_pdp)';
de2_nat_pdp = x_scp(:,4)-dqd2(t_pdp)';
de1_lgp_pdp = dq1_lgp_pdp_c-dqd1(t_lgp_pdp_c);
de2_lgp_pdp = dq2_lgp_pdp_c-dqd2(t_lgp_pdp_c);
de1_lgp_nat_pdp = dq1_lgp_nat_pdp_c-dqd1(t_lgp_nat_pdp_c);
de2_lgp_nat_pdp = dq2_lgp_nat_pdp_c-dqd2(t_lgp_nat_pdp_c);
de1_lgp_nat_pdp_v = dq1_lgp_nat_pdp_v-dqd1(t_lgp_var_nat_pdp);
de2_lgp_nat_pdp_v = dq2_lgp_nat_pdp_v-dqd2(t_lgp_var_nat_pdp);

%% Compute controller norms 
tau_pdp_norm_vec = zeros(size(t_pdp));
for i = 1:length(t_pdp)
    tau_pdp_norm_vec(i) = norm(tau_pdp(t_pdp(i),[q1(i);q2(i)],[dq1(i);dq2(i)]));
end
tau_pdp_max = max(tau_pdp_norm_vec);
tau_pdp_mean = mean(tau_pdp_norm_vec);

tau_gp_pdp_c_norm_vec = zeros(size(t_lgp_pdp_c));
for i = 1:length(t_lgp_pdp_c)
    tau_gp_pdp_c_norm_vec(i) = norm(tau_lgp_pdp_const(t_lgp_pdp_c(i),x_lgp_pdp_c(i,1:2)',x_lgp_pdp_c(i,3:4)'));
end
tau_gp_pdp_c_max = max(tau_gp_pdp_c_norm_vec);
tau_gp_pdp_c_mean = mean(tau_gp_pdp_c_norm_vec);

tau_nat_pdp_norm_vec = sqrt(u1_sc.^2 + u2_sc.^2);
tau_nat_pdp_max = max(tau_nat_pdp_norm_vec);
tau_nat_pdp_mean = mean(tau_nat_pdp_norm_vec);

tau_lgp_nat_pdp_norm_vec = zeros(size(t_lgp_nat_pdp_c));
for i = 1:length(t_lgp_nat_pdp_c)
    tau_lgp_nat_pdp_norm_vec(i) = norm(tau_lgp_nat_pdp_const(t_lgp_nat_pdp_c(i),x_lgp_nat_pdp_c(i,1:2)',x_lgp_nat_pdp_c(i,3:4)'));
end
tau_lgp_nat_pdp_max = max(tau_lgp_nat_pdp_norm_vec);
tau_lgp_nat_pdp_mean = mean(tau_lgp_nat_pdp_norm_vec);

tau_lgp_nat_pdp_var_norm_vec = zeros(size(t_lgp_var_nat_pdp));
for i = 1:length(t_lgp_var_nat_pdp)
    tau_lgp_nat_pdp_var_norm_vec(i) = norm(tau_lgp_nat_pdp_var(t_lgp_var_nat_pdp(i),x_lgp_var_nat_pdp(i,1:2)',x_lgp_var_nat_pdp(i,3:4)'));
end
tau_lgp_nat_pdp_var_max = max(tau_lgp_nat_pdp_var_norm_vec);
tau_lgp_nat_pdp_var_mean = mean(tau_lgp_nat_pdp_var_norm_vec);

%% Evaluate steady-state metrics
T_ss = 10;
ind = (t_pdp >= T_ss);
tau_pdp_ss_L2 = sqrt(trapz(t_pdp(ind), tau_pdp_norm_vec(ind).^2));
e_pdp_ss = sqrt(e1_pdp(ind).^2 + e2_pdp(ind).^2);
de_pdp_ss = sqrt(de1_pdp(ind).^2 + de2_pdp(ind).^2);
e_pdp_ss_L2 = sqrt(trapz(t_pdp(ind), e1_pdp(ind).^2+e2_pdp(ind).^2+de1_pdp(ind).^2+de2_pdp(ind).^2));
tau_pdp_max_ss = max(tau_pdp_norm_vec(ind));
tau_pdp_mean_ss = mean(tau_pdp_norm_vec(ind));

ind = (t_lgp_pdp_c >= T_ss);
tau_lgp_pdp_ss_L2 = sqrt(trapz(t_lgp_pdp_c(ind), tau_gp_pdp_c_norm_vec(ind).^2));
e_lgp_pdp_c_ss = sqrt(e1_lgp_pdp(ind).^2 + e2_lgp_pdp(ind).^2);
de_lgp_pdp_c_ss = sqrt(de1_lgp_pdp(ind).^2 + de2_lgp_pdp(ind).^2);
e_lgp_pdp_c_ss_L2 = sqrt(trapz(t_lgp_pdp_c(ind), e1_lgp_pdp(ind).^2+e2_lgp_pdp(ind).^2+de1_lgp_pdp(ind).^2+de2_lgp_pdp(ind).^2));
tau_lgp_pdp_max_ss = max(tau_gp_pdp_c_norm_vec(ind));
tau_lgp_pdp_mean_ss = mean(tau_gp_pdp_c_norm_vec(ind));

ind = (t_pdp >= T_ss);
tau_nat_pdp_ss_L2 = sqrt(trapz(t_pdp(ind), tau_nat_pdp_norm_vec(ind).^2));
e_nat_pdp_ss = sqrt(e1_nat_pdp(ind).^2 + e2_nat_pdp(ind).^2);
de_nat_pdp_ss = sqrt(de1_nat_pdp(ind).^2 + de2_nat_pdp(ind).^2);
e_nat_pdp_ss_L2 = sqrt(trapz(t_pdp(ind), e1_nat_pdp(ind).^2+e2_nat_pdp(ind).^2+de1_nat_pdp(ind).^2+de2_nat_pdp(ind).^2));
tau_nat_pdp_max_ss = max(tau_nat_pdp_norm_vec(ind));
tau_nat_pdp_mean_ss = mean(tau_nat_pdp_norm_vec(ind));

ind = (t_lgp_nat_pdp_c >= T_ss);
tau_lgp_nat_pdp_ss_L2 = sqrt(trapz(t_lgp_nat_pdp_c(ind), tau_lgp_nat_pdp_norm_vec(ind).^2));
e_lgp_nat_pdp_ss = sqrt(e1_lgp_nat_pdp(ind).^2 + e2_lgp_nat_pdp(ind).^2);
de_lgp_nat_pdp_ss = sqrt(de1_lgp_nat_pdp(ind).^2 + de2_lgp_nat_pdp(ind).^2);
e_lgp_nat_pdp_ss_L2 = sqrt(trapz(t_lgp_nat_pdp_c(ind), e1_lgp_nat_pdp(ind).^2+e2_lgp_nat_pdp(ind).^2+de1_lgp_nat_pdp(ind).^2+de2_lgp_nat_pdp(ind).^2));
tau_lgp_nat_pdp_max_ss = max(tau_lgp_nat_pdp_norm_vec(ind));
tau_lgp_nat_pdp_mean_ss = mean(tau_lgp_nat_pdp_norm_vec(ind));

ind = (t_lgp_var_nat_pdp >= T_ss);
tau_lgp_nat_pdp_v_ss_L2 = sqrt(trapz(t_lgp_var_nat_pdp(ind), tau_lgp_nat_pdp_var_norm_vec(ind).^2));
e_lgp_nat_pdp_v_ss = sqrt(e1_lgp_nat_pdp_v(ind).^2 + e2_lgp_nat_pdp_v(ind).^2);
de_lgp_nat_pdp_v_ss = sqrt(de1_lgp_nat_pdp_v(ind).^2 + de2_lgp_nat_pdp_v(ind).^2);
e_lgp_nat_pdp_v_ss_L2 = sqrt(trapz(t_lgp_var_nat_pdp(ind), e1_lgp_nat_pdp_v(ind).^2+e2_lgp_nat_pdp_v(ind).^2+de1_lgp_nat_pdp_v(ind).^2+de2_lgp_nat_pdp_v(ind).^2));
tau_lgp_nat_pdp_var_max_ss = max(tau_lgp_nat_pdp_var_norm_vec(ind));
tau_lgp_nat_pdp_var_mean_ss = mean(tau_lgp_nat_pdp_var_norm_vec(ind));

disp('done')

%
fprintf('-----------------------------------------------------------------------------------------------------------------------------\n')
fprintf('\t\t\t\t|\tPD+\t\t\tL-GP PD+\tnat-PD+\t\tL-GP nat-PD+\tL-GP var-nat-PD+\n')
fprintf('-----------------------------------------------------------------------------------------------------------------------------\n')
fprintf('||K_p||\t\t\t|\t%.0f\t\t\t%.0f\t\t\t%.0f\t\t\t%.0f\n', norm(Kp), norm(Kp_s), 0, 0)
fprintf('||K_d||\t\t\t|\t%.0f\t\t\t%.0f\t\t\t%.0f\t\t\t%.0f\n', norm(Kd), norm(Kd_s), norm(Kd_s), norm(Kd_s))
fprintf('-----------------------------------------------------------------------------------------------------------------------------\n')
fprintf('t >= 10s:\t\t|\n')
fprintf('||[tau]||_L^2\t|\t%.3f\t\t%.3f\t\t%.3f\t\t%.3f\t\t\t%.3f\n', tau_pdp_ss_L2, tau_lgp_pdp_ss_L2, tau_nat_pdp_ss_L2, tau_lgp_nat_pdp_ss_L2, tau_lgp_nat_pdp_v_ss_L2)
fprintf('max(||tau||)\t|\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t\t%.2f\n', tau_pdp_max_ss, tau_lgp_pdp_max_ss, tau_nat_pdp_max_ss, tau_lgp_nat_pdp_max_ss, tau_lgp_nat_pdp_var_max_ss)
fprintf('mean(||tau||)\t|\t%.3f\t\t%.3f\t\t%.3f\t\t%.3f\t\t\t%.3f\n', tau_pdp_mean_ss, tau_lgp_pdp_mean_ss, tau_nat_pdp_mean_ss, tau_lgp_nat_pdp_mean_ss, tau_lgp_nat_pdp_var_mean_ss)
fprintf('||[e, de]||_L^2\t|\t%.3f\t\t%.3f\t\t%.3f\t\t%.3f\t\t\t%.3f\n', e_pdp_ss_L2, e_lgp_pdp_c_ss_L2, e_nat_pdp_ss_L2, e_lgp_nat_pdp_ss_L2, e_lgp_nat_pdp_v_ss_L2)
fprintf('max(||e(t)||\t|\t%.3f\t\t%.3f\t\t%.3f\t\t%.3f\t\t\t%.3f\n', max(e_pdp_ss), max(e_lgp_pdp_c_ss), max(e_nat_pdp_ss), max(e_lgp_nat_pdp_ss), max(e_lgp_nat_pdp_v_ss))
fprintf('max(||de||\t\t|\t%.3f\t\t%.3f\t\t%.3f\t\t%.3f\t\t\t%.3f\n', max(de_pdp_ss), max(de_lgp_pdp_c_ss), max(de_nat_pdp_ss), max(de_lgp_nat_pdp_ss), max(de_lgp_nat_pdp_v_ss))
fprintf('mean(||e(t)||\t|\t%.3f\t\t%.3f\t\t%.3f\t\t%.3f\t\t\t%.3f\n', mean(e_pdp_ss), mean(e_lgp_pdp_c_ss), mean(e_nat_pdp_ss), mean(e_lgp_nat_pdp_ss), mean(e_lgp_nat_pdp_v_ss))
fprintf('mean(||de||\t\t|\t%.3f\t\t%.3f\t\t%.3f\t\t%.3f\t\t\t%.3f\n', mean(de_pdp_ss), mean(de_lgp_pdp_c_ss), mean(de_nat_pdp_ss), mean(de_lgp_nat_pdp_ss), mean(de_lgp_nat_pdp_v_ss))
fprintf('-----------------------------------------------------------------------------------------------------------------------------\n')

%%
fig1 = figure(1); clf

subplot(121)
hold on
set(gca,'FontSize',ax_fs)

fc = 6;
rc = 0.2;
for i = [0:fc/2, fc+(1:fc/2)]
    xp1_i = xp1(qd1(i/fc*pi));
    xp2_i = xp2(qd1(i/fc*pi),qd2(i/fc*pi));
    yp1_i = yp1(qd1(i/fc*pi));
    yp2_i = yp2(qd1(i/fc*pi),qd2(i/fc*pi));
    plot([0 xp1_i], [0 yp1_i],'Color','k','LineWidth',mlw)
    plot([xp1_i xp2_i], [yp1_i yp2_i],'Color','k','LineWidth',mlw)
    plot(xp2_i, yp2_i,'x','MarkerSize',msz,'Color','k','LineWidth',mlw);
    pa = atan2(yp2_i,xp2_i-0.75);
    if i == 0
        text(xp2_i-3.75*rc*cos(pa),yp2_i+rc,['$t=0,\pi$'],'Interpreter','Latex','FontSize',fs)
    elseif i == 1
        text(xp2_i+0.05,yp2_i-4*rc*sin(pa),['$t=\frac{\pi}{6},\frac{5\pi}{6}$'],'Interpreter','Latex','FontSize',fs)
    elseif i == 2
        text(xp2_i-2.75*rc*cos(pa),yp2_i,['$t=\frac{\pi}{3},\frac{2\pi}{3}$'],'Interpreter','Latex','FontSize',fs)
    elseif i == 3
        text(xp2_i-0.25*rc*cos(pa),yp2_i-3*rc*sin(pa),['$t=\frac{\pi}{2}$'],'Interpreter','Latex','FontSize',fs)
    elseif i == 7
        text(xp2_i+0.05,yp2_i-4*rc*sin(pa),['$t=\frac{7\pi}{6},\frac{11\pi}{6}$'],'Interpreter','Latex','FontSize',fs)
    elseif i == 8
        text(xp2_i-2.75*rc*cos(pa),yp2_i,['$t=\frac{4\pi}{3},\frac{5\pi}{3}$'],'Interpreter','Latex','FontSize',fs)
    elseif i == 9
        text(xp2_i-0.25*rc*cos(pa),yp2_i-3*rc*sin(pa),['$t=\frac{3\pi}{2}$'],'Interpreter','Latex','FontSize',fs)    
    end
end

t_plot = linspace(0,2*pi,1e2);
plot(xp2(qd1(t_plot),qd2(t_plot)), yp2(qd1(t_plot),qd2(t_plot)),'--','LineWidth',mlw,'Color',c(2,:))
grid
xlabel('$x$ (m)','Interpreter','Latex','FontSize',fs)
ylabel('$y$ (m)','Interpreter','Latex','FontSize',fs)
title('Reference Trajectory','Interpreter','Latex','FontSize',fs)

tp = 6*pi;
tp_ctc = 3.35*pi;

subplot(122)
hold on
set(gca,'FontSize',ax_fs)
h2 = plot(q1(t_pdp<tp_ctc), dq1(t_pdp<tp_ctc),'LineWidth',mlw,'Color',c(3,:));
h3 = plot(q1_lgp_nat_pdp_v(t_lgp_var_nat_pdp<tp), dq1_lgp_nat_pdp_v(t_lgp_var_nat_pdp<tp),'LineWidth',mlw,'Color',c(1,:));
h1 = plot(qd1(t_plot), dqd1(t_plot), '--','LineWidth',mlw,'Color',c(2,:));
h4 = plot(q1_tr, dq1_tr, '+','MarkerSize',msz,'LineWidth',mlw,'Color',c(4,:));
h5 = plot(q1_v, dq1_v, '+','MarkerSize',msz,'LineWidth',mlw,'Color',c(5,:));

grid
xlabel('$q_1$ (rad)','Interpreter','Latex','FontSize',fs)
ylabel('$\dot{q}_1$ (rad/s)','Interpreter','Latex','FontSize',fs)
title('Tracking Control','Interpreter','Latex','FontSize',fs)
lh = legend([h1,h2,h3,h4,h5], 'Ref','PD+','var-nat','Train','Valid','Interpreter','Latex','FontSize',fs,'Location','Best');

e_norm_gp_pdp = sqrt(e1_lgp_pdp.^2+e2_lgp_pdp.^2);
de_norm_gp_pdp = sqrt(de1_lgp_pdp.^2+de2_lgp_pdp.^2);
e_norm_pdp = sqrt(e1_pdp.^2+e2_pdp.^2);
de_norm_pdp = sqrt(de1_pdp.^2+de2_pdp.^2);
e_norm_nat_pdp = sqrt(e1_nat_pdp.^2+e2_nat_pdp.^2);
de_norm_nat_pdp = sqrt(de1_nat_pdp.^2+de2_nat_pdp.^2);
e_norm_gp_nat_pdp = sqrt(e1_lgp_nat_pdp.^2+e2_lgp_nat_pdp.^2);
de_norm_gp_nat_pdp = sqrt(de1_lgp_nat_pdp.^2+de2_lgp_nat_pdp.^2);
e_norm_gp_nat_pdp_v = sqrt(e1_lgp_nat_pdp_v.^2+e2_lgp_nat_pdp_v.^2);
de_norm_gp_nat_pdp_v = sqrt(de1_lgp_nat_pdp_v.^2+de2_lgp_nat_pdp_v.^2);

fig2 = figure(2); clf;
subplot(131)
semilogy(t_pdp(t_pdp<tp), e_norm_pdp(t_pdp<tp),'LineWidth',mlw,'Color',c(3,:))
hold on
semilogy(t_lgp_pdp_c(t_lgp_pdp_c<tp), e_norm_gp_pdp(t_lgp_pdp_c<tp),'LineWidth',mlw,'Color',c(4,:))
semilogy(t(t<tp), e_norm_nat_pdp(t<tp),'LineWidth',mlw,'Color',c(5,:))
semilogy(t_lgp_nat_pdp_c(t_lgp_nat_pdp_c<tp), e_norm_gp_nat_pdp(t_lgp_nat_pdp_c<tp),'LineWidth',mlw,'Color',c(7,:))
semilogy(t_lgp_var_nat_pdp(t_lgp_var_nat_pdp<tp), e_norm_gp_nat_pdp_v(t_lgp_var_nat_pdp<tp),'LineWidth',mlw,'Color',c(1,:))
set(gca,'FontSize',ax_fs)
xlabel('$t$ (s)','Interpreter','Latex','FontSize',fs)
ylabel('$\left\|\mbox{\boldmath $e$}(t)\right\|_2$ (rad)','Interpreter','Latex','FontSize',fs)
title('Position Error $\mbox{\boldmath $e$}=\mbox{\boldmath $q$}-\mbox{\boldmath $q$}_{d}$','Interpreter','Latex','FontSize',fs)
grid
xlim([0 tp])
ylim([1e-3 1.11+0.09])

subplot(132)
semilogy(t_pdp(t_pdp<tp), de_norm_pdp(t_pdp<tp),'LineWidth',mlw,'Color',c(3,:))
hold on
semilogy(t_lgp_pdp_c(t_lgp_pdp_c<tp), de_norm_gp_pdp(t_lgp_pdp_c<tp),'LineWidth',mlw,'Color',c(4,:))
semilogy(t(t<tp), de_norm_nat_pdp(t<tp),'LineWidth',mlw,'Color',c(5,:))
semilogy(t_lgp_nat_pdp_c(t_lgp_nat_pdp_c<tp), de_norm_gp_nat_pdp(t_lgp_nat_pdp_c<tp),'LineWidth',mlw,'Color',c(7,:))
semilogy(t_lgp_var_nat_pdp(t_lgp_var_nat_pdp<tp), de_norm_gp_nat_pdp_v(t_lgp_var_nat_pdp<tp),'LineWidth',mlw,'Color',c(1,:))
set(gca,'FontSize',ax_fs)
xlabel('$t$ (s)','Interpreter','Latex','FontSize',fs)
ylabel('$\left\|\dot{\mbox{\boldmath $e$}}(t)\right\|_2$ (rad/s)','Interpreter','Latex','FontSize',fs)
title('Velocity Error $\dot{\mbox{\boldmath $e$}}=\dot{\mbox{\boldmath $q$}}-\dot{\mbox{\boldmath $q$}}_{d}$','Interpreter','Latex','FontSize',fs)
grid
xlim([0 tp])
ylim([1e-3 pi])

subplot(133)
plot(t, tau_pdp_norm_vec, 'LineWidth',mlw,'Color',c(3,:))
hold on
plot(t_lgp_pdp_c, tau_gp_pdp_c_norm_vec, 'LineWidth', mlw,'Color',c(4,:))
plot(t, tau_nat_pdp_norm_vec, 'LineWidth',mlw,'Color',c(5,:))
plot(t_lgp_nat_pdp_c, tau_lgp_nat_pdp_norm_vec, 'LineWidth', mlw,'Color',c(7,:))
plot(t_lgp_var_nat_pdp, tau_lgp_nat_pdp_var_norm_vec, 'LineWidth', mlw,'Color',c(1,:))
set(gca,'FontSize',ax_fs)
xlim([0 tp])
ylim([0 31])
xlabel('$t$ (s)','Interpreter','Latex','FontSize',fs)
ylabel('$\left\|\mbox{\boldmath $\tau$}(t)\right\|_2$ (Nm)','Interpreter','Latex','FontSize',fs)
title('Actuation','Interpreter','Latex','FontSize',fs)
grid
legend('SoA (PD+)','new PD+ with $\mbox{\boldmath $\mu_{\tau}$}$ (L-GP-PD+)','ours no L-GP (nat-PD+)','ours $\mbox{\boldmath $\mu_{\tau}$}$ (L-GP nat-PD+)','ours $\mbox{\boldmath $\mu_{\tau}$},\mbox{\boldmath $\Sigma_{\tau}$}$ (L-GP var-nat-PD+)','Interpreter','Latex','FontSize',fs,'Location','NorthEast')%,'Orientation','Horizontal')
