function [c,ceq] = nl_ineqs_beta_soro(alph,eps,the,phi,...
    kd,d,kp,kapl,m_l,m_u,u_min)%eps_t)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% c <= 0
% ceq = 0

%eps = (kd-phi/2 + beta*(kp-m_l)) / (kp+m_l*(1+2*beta^2)-beta*kd-the/2);
%kappa = kp + eps*kd - eps*beta*(m_l+m_u);

%%
%a_l = eps * (kp - the/2) - beta*kappa;
%d_l = kd - phi/2;

%ups = (a_l+d_l-(eps+beta)*m_u)/2;
%ups_min = ups - sqrt(ups^2 + (eps*beta*(m_u-m_l))^2/4 - a_l*(d_l-(eps+beta)*m_u));

kappa = kp + eps*(d+kd-alph*(m_l+m_u));
a = eps*(kp-the/2) - alph*kappa;
%a1 = kp - the/2 - alph * (d + kd - alph*(m_l+m_u));

delta_m = m_u - m_l;
ups_min = a - (eps+alph)*delta_m/4 * (1 + sqrt(1 + 4*eps^2*alph^2/(eps+alph)^2));

%%
ceq = [];

c = [-a, u_min-ups_min, 1.01*eps-sqrt((kapl+kappa)/m_u), -kappa, -eps, -phi];%, alph/eps-2];

%%

end