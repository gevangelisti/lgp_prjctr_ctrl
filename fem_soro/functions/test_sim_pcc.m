
ivd_func = @(x) ivd_MCg_pcc_n4(mu,L,Iz,9.81,x(1:Ncc),x(Ncc+1:2*Ncc),a_pcc-diag(k_pcc)*x(1:Ncc)-diag(d_pcc)*x(Ncc+1:2*Ncc));
xy_func = @(cosc_vec,q_vec,sinc_vec) pcc_xy_n4(L(1),L(2),L(3),L(4),cosc_vec(1),cosc_vec(2),cosc_vec(3),cosc_vec(4),q_vec(1),q_vec(2),q_vec(3),sinc_vec(1),sinc_vec(2),sinc_vec(3),sinc_vec(4));


% Simulate
ode_fun = @(t,x) [x(Ncc+1:2*Ncc); ivd_func(x)];
x0 = zeros(2*Ncc,1);
if Ncc <= 4
    [t,x] = ode45(ode_fun,[0 tsim],x0);
else
    [t,x] = ode15s(ode_fun,[0 tsim],x0);
end

% Plot
plot_soro_cs_ws;