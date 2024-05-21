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
N_mc_runs = 4;
tend = 12;

%%
qd_bar = pi/2*sqrt(2);
e_bar = pi/4*sqrt(2);
q_bar = qd_bar+e_bar;
dq_bar = qd_bar+e_bar;
de_bar = 2*pi;

%%
fq2e = 2;
vec_q = linspace(-q_bar,q_bar,2e2);
[Q1,Q2] = meshgrid(vec_q,vec_q);
Lambda_max_est = zeros(size(Q1));
Lambda_min_est = zeros(size(Q1));
vec_dq = linspace(-dq_bar,dq_bar,2e2);
[DQ1,DQ2] = meshgrid(vec_dq,vec_dq);
Lambda_D_max_est = zeros(size(Q1));
Lambda_D_min_est = zeros(size(Q1));
GV_est = zeros(size(Q1));
Lambda_KG_max_est = zeros(size(Q1));
Lambda_KG_min_est = zeros(size(Q1));
for i = 1:size(Q1,1)
    for j = 1:size(Q1,2)
        
        q_g = [Q1(i,j);Q2(i,j)]./fq2e;
        [g_est,GV_est(i,j)] = V_lgp_mex(q_g, X_tr, rho_g_sq, P_g_inv2, A_inv_times_y_min_mu, K_times_a, est_params, D, N);

        % g
        Lambda_KG_min_est(i,j) = min(g_est./q_g);
        Lambda_KG_max_est(i,j) = max(g_est./q_g);
        
        [M_est, D_est] = MD_lgp_mex(fq2e*q_g, [DQ1(i,j);DQ2(i,j)]./fq2e, X_tr, ddq_meas, LLT_Sigma_sq, lambda_isq, L_Sigma_sq_dmp, lambda_isq_dmp, A_inv_times_y_min_mu, est_params, D, N);
        % M
        eig_est = eig(M_est);
        Lambda_min_est(i,j) = min(eig_est);
        Lambda_max_est(i,j) = max(eig_est);
        % D
        eig_est = eig(D_est);
        Lambda_D_min_est(i,j) = min(eig_est);
        Lambda_D_max_est(i,j) = max(eig_est);
    end
end

norm_q = sqrt(Q1.^2+Q2.^2);
ind = (norm_q <= q_bar);
m_lower = min(min(Lambda_min_est(ind)));
[M_0, ~] = MD_lgp_mex(zeros(N,1), zeros(N,1), X_tr, ddq_meas, LLT_Sigma_sq, lambda_isq, L_Sigma_sq_dmp, lambda_isq_dmp, A_inv_times_y_min_mu, est_params, D, N);
m_0 = max(eig(M_0));
m_1 = 0*1e-2/q_bar;
m_2 = 0*1e-2/q_bar^2;

d_lower_est = min(min(Lambda_D_min_est));
d_upper_est = max(max(Lambda_D_max_est));

kappa_G_upper_est = max(max(2*GV_est./(norm_q./fq2e).^2));
kappa_G_lower_est = min(min(2*GV_est./(norm_q./fq2e).^2));

k_G_lower_est = kappa_G_lower_est;
k_G_upper_est = kappa_G_upper_est;

mu_hat = @(q) m_0 + m_1*norm(q) + m_2*norm(q)^2;

fig4 = figure(4); clf; hold on;

for w_var = 0:1
    
    fig3 = figure(3+2*w_var);
    clf
    subplot(211)
    set(gca,'yscale','log','FontSize',ax_fs)
    hold on
    subplot(212)
    set(gca,'yscale','log','FontSize',ax_fs)
    hold on
    sh_alpha = 0.15;

    a0f = @(alph) alph * (Kp_s(1,1)+w_var);
    a1f = @(alph,the) (Kp_s(1,1)+w_var) - the/2 - alph * (d_lower_est + Kd_s(1,1)+w_var - alph*(m_0+m_lower));
    delta_m = m_0 - m_lower;
    phif = @(alph,eps,the) 2 * (d_lower_est + Kd_s(1,1)+w_var - a1f(alph,the)*eps + a0f(alph)) - (eps+alph)*(m_0+m_lower);
    af = @(alph,eps,the) -a0f(alph) + a1f(alph,the)*eps;
    upsf = @(alph,eps,the) af(alph,eps,the) - (eps+alph)*delta_m/4 * (1 + sqrt(1 + 4*eps^2*alph^2/(eps+alph)^2));
    
    kappa = @(alph,eps) (Kp_s(1,1)+w_var) + eps*(Kd_s(1,1)+w_var+d_lower_est) - eps*alph*(m_lower+m_0);
    mu_l = @(alph,eps) (kappa_G_lower_est+kappa(alph,eps)+m_lower)/2 - 1/2*sqrt((kappa_G_lower_est+kappa(alph,eps)-m_lower)^2+(2*eps*m_lower)^2);
    rho_func = @(alph,eps,the,phi) sqrt((eps/the + 1/phi)/(2*mu_l(alph,eps)));
    
    rng(0)
    x_test = ga(@(x)rho_func(x(1),x(2),x(3),phif(x(1),x(2),x(3))) + 1/x(1),3,[],[],[],[],[0.1,0.1,1],[1,10,10],...
        @(x) nl_ineqs_beta(x(1),x(2),x(3),phif(x(1),x(2),x(3)),Kd_s(1,1)+w_var,d_lower_est,Kp_s(1,1)+w_var,kappa_G_lower_est,m_lower,m_0));
    
    alpha_ct = x_test(1); eps_ct = x_test(2); the_ct  = x_test(3); phi_ct = phif(x_test(1),x_test(2),x_test(3));
    mu_lower = mu_l(alpha_ct,eps_ct);
    rho_ct = rho_func(alpha_ct,eps_ct,the_ct,phi_ct);
    
    kappa_ct = kappa(alpha_ct,eps_ct);
    
    %%
    rng(1)
    x0_mat = pi/3*randn(2*N,2*N_mc_runs) + repmat([zeros(2,1);pi/2*ones(2,1)], [1,2*N_mc_runs]);
    x_t_mat = cell(1,N_mc_runs); t_mat = cell(1,N_mc_runs); V_t_mat = cell(1,N_mc_runs);
    
    c_cntr = 1;
    
    N_mc_run_order = [4,1,2,3];
    for n_mc = N_mc_run_order
    
        x0 = x0_mat(:,n_mc);
        
        if w_var == 0
            % Proposed structure-preserving tracking control
            tau_lgp_nat_pdp_const = @(t,q,dq) tau_lgp_nat_pdp_mex(t, q, dq, pi/2, 1, 1*Kp_s, 1*Kd_s, X_tr, ddq_meas, ...
                LLT_Sigma_sq, lambda_isq, rho_g_sq, P_g_inv2, L_Sigma_sq_dmp, lambda_isq_dmp, A_inv_times_y_min_mu, K_times_a, est_params, D, N);
            ode_robot_nat = @(t,x) [x(3); x(4); M(x(1:2)) \ (tau_lgp_nat_pdp_const(t,x(1:2),x(3:4)) - C(x(1:2),x(3:4))*x(3:4) - G(x(1:2)) - Dmp(x(3:4))*x(3:4))];
            [t_lgp_nat_pdp,x_lgp_nat_pdp] = ode45(ode_robot_nat,linspace(0,tend,F_s/10*tend+1),x0);
            t_mat{n_mc} = t_lgp_nat_pdp;
            t_n = t_lgp_nat_pdp;
            x_n = x_lgp_nat_pdp;
            e1_lgp_nat_pdp = x_lgp_nat_pdp(:,1)-pi/2*sin(t_lgp_nat_pdp);
            e2_lgp_nat_pdp = x_lgp_nat_pdp(:,2)-pi/2*sin(t_lgp_nat_pdp);
            de1_lgp_nat_pdp = x_lgp_nat_pdp(:,3)-pi/2*cos(t_lgp_nat_pdp);
            de2_lgp_nat_pdp = x_lgp_nat_pdp(:,4)-pi/2*cos(t_lgp_nat_pdp);
            norm_e_sq = e1_lgp_nat_pdp.^2+e2_lgp_nat_pdp.^2+de1_lgp_nat_pdp.^2+de2_lgp_nat_pdp.^2;
        else
            % Proposed uncertainty-adaptive variant
            tau_lgp_nat_pdp_var = @(t,q,dq) tau_lgp_var_nat_pdp_mex(t, q, dq, pi/2, 1, ...
                1e2*eye(2), 2e-2*eye(2), 1/sqrt(2e-2*(1/1e0-1/1e2))*eye(2), 1*Kp_s, 1*Kd_s, A, ... 
                X_tr, ddq_meas, LLT_Sigma_sq, lambda_isq, rho_g_sq, P_g_inv2, L_Sigma_sq_dmp, lambda_isq_dmp, A_inv_times_y_min_mu, K_times_a, est_params, D, N);
            ode_robot_nat_var = @(t,x) [x(3); x(4); M(x(1:2)) \ (tau_lgp_nat_pdp_var(t,x(1:2),x(3:4)) - C(x(1:2),x(3:4))*x(3:4) - G(x(1:2)) - Dmp(x(3:4))*x(3:4))];
            [t_lgp_var_nat_pdp,x_lgp_var_nat_pdp] = ode45(ode_robot_nat_var,linspace(0,tend,F_s/10*tend+1),x0);
            t_mat{n_mc} = t_lgp_var_nat_pdp;
            t_n = t_lgp_var_nat_pdp;
            x_n = x_lgp_var_nat_pdp;
            e1_lgp_nat_pdp_v = x_lgp_var_nat_pdp(:,1)-pi/2*sin(t_lgp_var_nat_pdp);
            e2_lgp_nat_pdp_v = x_lgp_var_nat_pdp(:,2)-pi/2*sin(t_lgp_var_nat_pdp);
            de1_lgp_nat_pdp_v = x_lgp_var_nat_pdp(:,3)-pi/2*cos(t_lgp_var_nat_pdp);
            de2_lgp_nat_pdp_v = x_lgp_var_nat_pdp(:,4)-pi/2*cos(t_lgp_var_nat_pdp);
            norm_e_sq = e1_lgp_nat_pdp_v.^2+e2_lgp_nat_pdp_v.^2+de1_lgp_nat_pdp_v.^2+de2_lgp_nat_pdp_v.^2;
        end
    
        % Plots
        GV_gp_est_traj = zeros(size(t_n));
        V_lgp_nat_pdp_traj = zeros(size(t_n));
        dVdt_lgp_nat_pdp_traj = zeros(size(t_n));
        dVdt_bound1_deriv_th1 = zeros(size(t_n));
        dVdt_bound_w_cross_prod_deriv_th1 = zeros(size(t_n));
        dVdt_bound_w_taueb_deriv_th1 = zeros(size(t_n));
        ct = zeros(size(t_n));
        mu_l_temp = zeros(size(t_n)); mu_u_temp = zeros(size(t_n));
        alpha_t = zeros(size(t_n));
    
        gamma_1_t = zeros(size(t_n)); gamma_2_t = zeros(size(t_n)); gamma_3_t = zeros(size(t_n));
        xi_bar_t = zeros(size(t_n)); dxi_bar_t = zeros(size(t_n)); xi_min_t = zeros(size(t_n)); dxi_min_t = zeros(size(t_n));
        eps_t = 0.3; the_t = 0.5; dth_t = 2; vt_t = 2.5e-1; gam_t = 0.15; del_t = 0.15;%7.5e-1;
    
        lambda_mat_test_ct = zeros(size(t_n)); C_norm_vec = zeros(size(t_n));
    
        for i = 1:length(t_n)
            if w_var == 0
                e_i = [e1_lgp_nat_pdp(i);e2_lgp_nat_pdp(i)];
                de_i = [de1_lgp_nat_pdp(i);de2_lgp_nat_pdp(i)];
            else
                e_i = [e1_lgp_nat_pdp_v(i);e2_lgp_nat_pdp_v(i)];
                de_i = [de1_lgp_nat_pdp_v(i);de2_lgp_nat_pdp_v(i)];
            end
    
            [M_est,C_est,D_est,g_est,G_est,G_est_e,Sigma_lgp,g_est_e,D_est_de,e,de] = eval_lgp_w_sigma_mex(t_n(i), x_n(i,1:2)', x_n(i,3:4)', pi/2, 1, ...
                1e2*eye(2), 2e-2*eye(2), 1/sqrt(2e-2*(1/1e0-1/1e2))*eye(2), 1*Kp_s, 1*Kd_s, A, ... 
                X_tr, ddq_meas, LLT_Sigma_sq, lambda_isq, rho_g_sq, P_g_inv2, L_Sigma_sq_dmp, lambda_isq_dmp, A_inv_times_y_min_mu, K_times_a, est_params, D, N);
            
            d_est_de = D_est_de * de;
            GV_gp_est_traj(i) = G_est_e + 1/2*e_i'*Kp_s*e_i;
            V_lgp_nat_pdp_traj(i) = 1/2*[e_i;de_i]' * [kappa_ct*eye(2), eps_ct*M_est; eps_ct*M_est, M_est] * [e_i;de_i] + G_est_e;
            
            q_i = x_n(i,1:2)';
            dq_i = x_n(i,3:4)';
            
            mu_l_temp(i) = 2*G_est_e/(norm(e_i)^2+norm(de_i)^2) + (kappa_ct+min(eig(M_est)))/2 - 1/2*sqrt((kappa_ct-min(eig(M_est)))^2+(2*eps_ct*min(eig(M_est)))^2);
            mu_u_temp(i) = 2*G_est_e/(norm(e_i)^2+norm(de_i)^2) + (kappa_ct+max(eig(M_est)))/2 + 1/2*sqrt((kappa_ct-max(eig(M_est)))^2+(2*eps_ct*max(eig(M_est)))^2);
            
            eta_e_i = heaviside(e_i'*g_est)*e_i'*g_est;
            d_est = D_est*dq_i;
            eta_de_i = heaviside(de_i'*d_est)*de_i'*d_est;
            omega_e_i = heaviside(e_i'*g_est)*de_i'*((e_i*e_i')./(1e-3+norm(e_i)^2))*g_est;
            omega_de_i = heaviside(de_i'*d_est)*e_i'*((de_i*de_i')./(1e-3+norm(de_i)^2))*d_est;
    
            % time-variant alpha calc
            D_tilde = D_est_de - min(diag(D_est_de))*eye(2);
            R_mat = [zeros(2), eps_ct/2*(D_tilde-C_est'); eps_ct/2*(D_tilde-C_est), D_tilde];
            kp = Kp_s(1,1);
            if w_var
                K1 = 1e2*eye(2); K2 = 2e-2*eye(2); K3 = 1/sqrt(2e-2*(1/1e0-1/1e2))*eye(2);
                Kp_v = K1 * (eye(2) - ( ( K3*(K2+Sigma_lgp)*K3 + K1 )\K1 ) );
                kp = min(eig(Kp_s + Kp_v));
                Kp_tilde = Kp_s + Kp_v - kp*eye(N);
                R_mat = R_mat + [eps_ct*Kp_tilde, (1+eps_ct)/2*Kp_tilde; (1+eps_ct)/2*Kp_tilde, Kp_tilde];
            end
            r_min = min(eig(R_mat));
            gam_t = min(diag(D_est_de)) + kp - phi_ct/2;
            b_t = kp + eps_ct*(min(diag(D_est_de)) + kp) - kappa_ct;
            X = eps_ct * (e_i'*g_est_e+eta_e_i+omega_de_i) + eta_de_i + omega_e_i;
            
            u0 = (eps_ct * (kp-the_ct/2) + gam_t - eps_ct*min(eig(M_est)))/2; u1 = (kappa_ct+min(eig(M_est)))/2;
            p0 = gam_t - eps_ct * (kp-the_ct/2 + min(eig(M_est))); p1 = kappa_ct - min(eig(M_est)); q = 2*eps_ct*min(eig(M_est));
            AA = u1 + G_est_e/(norm(e_i)^2+norm(de_i)^2); BB = u0 + r_min + X/(norm(e_i)^2+norm(de_i)^2);
            A2 = 4*AA^2 - p1^2 - q^2; A1 = 2*(q*b_t - p0*p1 - 4*AA*BB); A0 = 4*BB^2 - p0^2 - b_t^2;
            alpha_min = 1/(2*A2) * (-A1 - sqrt(A1^2 - 4*A2*A0)); a_t = eps_ct*(kp - the_ct/2)-alpha_min*kappa_ct;
            u_min = (a_t+gam_t-(eps_ct+alpha_min)*min(eig(M_est)))/2; upsil_min = u_min - 1/2*sqrt((a_t - gam_t + (eps_ct+alpha_min)*min(eig(M_est)))^2 + (2*eps_ct*alpha_min*min(eig(M_est))-b_t)^2);
            
            u0 = (eps_ct * (kp-the_ct/2) + gam_t - eps_ct*max(eig(M_est)))/2; u1 = (kappa_ct+max(eig(M_est)))/2;
            p0 = gam_t - eps_ct * (kp-the_ct/2 + max(eig(M_est))); p1 = kappa_ct - max(eig(M_est)); q = 2*eps_ct*max(eig(M_est));
            AA = u1 + G_est_e/(norm(e_i)^2+norm(de_i)^2); BB = u0 + r_min + X/(norm(e_i)^2+norm(de_i)^2);
            A2 = 4*AA^2 - p1^2 - q^2; A1 = 2*(q*b_t - p0*p1 - 4*AA*BB); A0 = 4*BB^2 - p0^2 - b_t^2;
            alpha_max = 1/(2*A2) * (-A1 - sqrt(A1^2 - 4*A2*A0)); a_t = eps_ct*(kp - the_ct/2)-alpha_max*kappa_ct;
            u_max = (a_t+gam_t-(eps_ct+alpha_max)*max(eig(M_est)))/2; upsil_max = u_max - 1/2*sqrt((a_t - gam_t + (eps_ct+alpha_max)*max(eig(M_est)))^2 + (2*eps_ct*alpha_max*max(eig(M_est))-b_t)^2);
    
            if upsil_min < upsil_max
                alpha_t(i) = alpha_min;
            else
                alpha_t(i) = alpha_max;
            end
            
        end
        
        if n_mc <= 7
            c_temp = c(N_mc_runs-c_cntr+1,:);
        else 
            c_temp = rand(1,3);
        end
    
        figure(3+2*w_var)
        subplot(211)
        fill([t_n', fliplr(t_n')], [1/2*(mu_l_temp.*norm_e_sq)', fliplr(1/2*(mu_u_temp.*norm_e_sq)')],c_temp,'EdgeColor','none','FaceAlpha',sh_alpha);
        semilogy(t_n,1/2*(mu_l_temp.*norm_e_sq)','--','LineWidth',mlw, 'Color',c_temp)
        semilogy(t_n,1/2*(mu_u_temp.*norm_e_sq)','--','LineWidth',mlw, 'Color',c_temp)
        V_t_mat{n_mc} = V_lgp_nat_pdp_traj;
    
        alpha_t_int = zeros(size(t_n));
        for i = 2:length(t_n)
            alpha_t_int(i) = alpha_t_int(i-1) + trapz(t_n(i-1:i), alpha_t(i-1:i)');
        end
        if w_var == 0 && n_mc == N_mc_run_order(1)
            lgp_error_fac = 2 * min(mu_l_temp) / (1/phi_ct + eps_ct/the_ct);
        end
        V_border_ct = lgp_error_fac * (1/phi_ct + eps_ct/the_ct) / 4;
        x_V_ct_bound_time_var = sqrt(2*V_border_ct./mu_l_temp) + sqrt(2*(V_lgp_nat_pdp_traj(1)-V_border_ct)./mu_l_temp) .* exp(-alpha_t_int);
        
        subplot(212)
        if w_var == 0
            fill([t_n', fliplr(t_n')], [8e-3*ones(size(t_n')), fliplr(x_V_ct_bound_time_var')],c_temp,'EdgeColor','none','FaceAlpha',sh_alpha);
        else
            fill([t_n', fliplr(t_n')], [2.5e-3*ones(size(t_n')), fliplr(x_V_ct_bound_time_var')],c_temp,'EdgeColor','none','FaceAlpha',sh_alpha);
        end
        semilogy(t_n,x_V_ct_bound_time_var,'--','LineWidth',mlw, 'Color',c_temp)
        x_t_mat{n_mc} = sqrt(norm_e_sq);
    
        figure(4)
        subplot(2,1,1)
        if w_var
            semilogy(t_n, alpha_t, 'LineWidth', mlw, 'Color', c_temp)
        else
            semilogy(t_n, alpha_t, '--', 'LineWidth',mlw, 'Color',c_temp)
        end
        hold on

        subplot(2,1,2)
        if w_var
            semilogy(t_n, x_V_ct_bound_time_var, 'LineWidth', mlw, 'Color', c_temp)
        else
            semilogy(t_n, x_V_ct_bound_time_var, '--', 'LineWidth',mlw, 'Color',c_temp)
        end
        hold on

        c_cntr = c_cntr + 1;
        
    end
    %%
    figure(3+2*w_var)
    subplot(211)
    grid on;
    xlim([0 tend])
    ylabel('$V(\mbox{\boldmath $e$}(t),\dot{\mbox{\boldmath $e$}}(t),t)$ (J)','Interpreter','Latex','FontSize',fs)
    c_cntr = 1;
    for n_mc = N_mc_run_order
        semilogy(t_mat{n_mc},V_t_mat{n_mc},'LineWidth',mlw, 'Color',c(N_mc_runs-c_cntr+1,:)); c_cntr = c_cntr + 1;
    end
    
    subplot(212)
    c_cntr = 1;
    for n_mc = N_mc_run_order
        semilogy(t_mat{n_mc},x_t_mat{n_mc},'LineWidth',mlw, 'Color',c(N_mc_runs-c_cntr+1,:))
        c_cntr = c_cntr + 1;
    end
    grid on;
    xlim([0 tend])
    ylabel('$||[\mbox{\boldmath $e$}^\top \, s\dot{\mbox{\boldmath $e$}}^\top]||$ (rad)','Interpreter','Latex','FontSize',fs)
    xlabel('$t$ (s)','Interpreter','Latex','FontSize',fs)
    
end

figure(4)
subplot(2,1,1)
set(gca,'yscale','log','FontSize',ax_fs); grid on; xlim([0 tend]); yticks([1e0 1e2]);
ylabel('$\alpha(t)$ (1/s)','Interpreter','Latex','FontSize',fs)

subplot(2,1,2)
set(gca,'yscale','log','FontSize',ax_fs); grid on; xlim([0 tend]); ylim([1e-1 6e0]); yticks([1e-1 1e0]);
xlabel('$t$ (s)','Interpreter','Latex','FontSize',fs)
ylabel('$\varrho(t)$ (rad)','Interpreter','Latex','FontSize',fs)
