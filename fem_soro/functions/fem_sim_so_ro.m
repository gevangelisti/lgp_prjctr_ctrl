function [qout,dqout] = fem_sim_so_ro(soro,N,k,d,x0,tsim,a)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ode_fun = @(t,x) [...
    x(N+1:end); ...
    FDab_wDK( ...
    soro, x(1:N), x(N+1:end), -k*eye(N)*x(1:N) - d*eye(N)*x(N+1:end) + ...
    a*ones(N,1)) ...
    ];

[~,x] = ode15s(ode_fun,tsim,x0);

qout = x(:,1:N)';
dqout = x(:,N+1:end)';

end